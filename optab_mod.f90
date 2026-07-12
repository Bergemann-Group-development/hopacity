MODULE optab_mod
  USE os_mod
  USE prec_mod
  USE const_mod
  USE spectrum_mod, only: spectrum_type
  USE op_math_mod, only: searchsorted, pca_grid_2d
  USE interpol_mod, only: hermite_coefs, pchip_mono
  implicit none

  !=============================================================================
  ! Externally-read binned opacity table (special case): 2D (log temp, log rho)
  ! interpolation, one independent lookup per wavelength bin. Node derivatives
  ! along both axes are precomputed once in build_derivatives (called from init)
  ! so interpolate is a pure Hermite evaluation with stored slopes.
  !=============================================================================
  type :: optab_type
    real(4), dimension(:), allocatable :: ltemp, lrho
    real(4), dimension(:,:,:), allocatable :: chi, src
    ! precomputed monotone node slopes: d/dlogT and d/dlogrho of chi and src
    real(4), dimension(:,:,:), allocatable :: dchi_t, dchi_r, dsrc_t, dsrc_r
    integer :: lunit, ntemp, nrho, nbin
    logical :: initialised=.false.
  contains
    procedure :: init
    procedure :: build_derivatives => optab_build_derivatives
    procedure :: interpolate
  end type

  !=============================================================================
  ! Orthogonal PCA coordinate transform of (lnT, lnNe):
  !   c = rot^T ((lnT,lnNe) - mu),   coord_k(i) = c0(k) + (i-1)*dc(k)
  ! Ties a sampled opacity table's (temp, ne) node grid to the physical
  ! coordinates it was built from, so a query point can be rotated into the
  ! same node-index space (see sample_interpolate). Built once, in
  ! pca_opac_grid, and carried from the composition patch that owns the grid
  ! down to every rt_patch that samples it (composition_type%pca).
  !=============================================================================
  type :: pca_transform_type
    real(8) :: mu(2)=0            ! mean of (lnT, lnNe)
    real(8) :: rot(2,2)=0          ! PC directions (unit eigenvectors as columns)
    real(8) :: c0(2)=0, dc(2)=1    ! min coordinate & spacing along each PC axis
  end type

  !=============================================================================
  ! Runtime background-opacity sampling table: 3D (wavelength, temp, log ne).
  ! The temp axis is regularly spaced (linear) and the ne axis is regularly
  ! spaced in log(ne); the wavelength axis is irregular. Monotone PCHIP node
  ! slopes along every dimension are precomputed and stored, so interpolate is
  ! an O(1) separable tensor-product cubic Hermite evaluation per query.
  !=============================================================================
  type :: sample_optab_type
    integer :: nnu=0, ntemp=0, nne=0
    type(pca_transform_type) :: pca       ! (lnT, lnNe) -> table-axis coordinates
    real(8), dimension(:), allocatable :: lam                 ! (nnu) irregular wavelength nodes
    real(4), dimension(:,:,:), allocatable :: val             ! (nnu,ntemp,nne) chi/chi500
    real(4), dimension(:,:,:), allocatable :: dl, dt, de      ! precomputed node slopes
    logical :: initialised=.false.
  contains
    procedure :: build_derivatives => sample_build_derivatives
    procedure :: interpolate => sample_interpolate
    procedure :: dealloc => sample_dealloc
  end type

  type(optab_type), public, save :: optab

contains

  !===============================================================================
  ! initialises opacity table.
  !===============================================================================
  SUBROUTINE init (self, opac_table, sp)
    class(optab_type) :: self
    class(spectrum_type) :: sp
    character(len=*), intent(in) :: opac_table
    integer :: i, lunit, iostat
    !---------------------------------------------------------------------------
    print '(A,A)', 'Reading opacity table: ', trim(opac_table)

    call os%check_file(trim(opac_table)//'.txt')
    open(file=trim(opac_table)//'.txt', newunit=lunit, form='formatted')
    read(lunit, *) self%ntemp, self%nrho, self%nbin
    close(lunit)

    allocate(self%lrho(self%nrho))
    allocate(self%ltemp(self%ntemp))
    allocate(self%chi(self%nbin, self%ntemp, self%nrho))
    allocate(self%src(self%nbin, self%ntemp, self%nrho))

    call os%check_file(trim(opac_table)//'.bin')
    open(file=trim(opac_table)//'.bin', newunit=lunit, form='unformatted', access='stream')
    read(lunit) self%ltemp, self%lrho, self%chi, self%src
    close(lunit)

    call self%build_derivatives()

    sp%nnu = self%nbin
    sp%name = 'bins'
    call sp%reset(lam=[(dble(i), i=1, self%nbin)])
    self%initialised = .true.
  END SUBROUTINE init

  !===============================================================================
  ! Precompute monotone node slopes of the binned table along the (irregular)
  ! log-temperature and log-density axes, for both chi and src.
  !===============================================================================
  SUBROUTINE optab_build_derivatives (self)
    class(optab_type) :: self
    !---------------------------------------------------------------------------
    allocate(self%dchi_t, self%dchi_r, self%dsrc_t, self%dsrc_r, mold=self%chi)
    call node_slopes_dim(self%ltemp, self%chi, 2, self%dchi_t)
    call node_slopes_dim(self%lrho , self%chi, 3, self%dchi_r)
    call node_slopes_dim(self%ltemp, self%src, 2, self%dsrc_t)
    call node_slopes_dim(self%lrho , self%src, 3, self%dsrc_r)
  END SUBROUTINE optab_build_derivatives

  !===============================================================================
  ! Interpolate the binned table at the requested (temp, rho). Separable 2D
  ! tensor-product cubic Hermite (collapse log-temp then log-rho) using the
  ! precomputed node slopes; mixed partials are ignored. chi and eta are
  ! shaped (nbin, ncells...) with the bin as the first, contiguous index.
  !===============================================================================
  SUBROUTINE interpolate (self, chi, eta, temp, rho)
    class(optab_type) :: self
    real(4), contiguous, dimension(:,:,:,:), target :: chi, eta
    real(4), contiguous, dimension(:,:,:), target :: temp, rho
    real(4), dimension(:), pointer :: chip, etap
    real(4), dimension(:), pointer :: temp1, rho1
    integer :: c, f, it, ir, ncells, ntot, nb
    real(4) :: lt, lr, ht, hr, tt, tr, pt, pr
    real(4) :: tc0, tc1, tcm0, tcm1, rc0, rc1, rcm0, rcm1
    real(4) :: v0, v1, g0, g1, cval, sval
    !---------------------------------------------------------------------------
    temp1(1:size(temp)) => temp(:,:,:)
    rho1 (1:size(rho))  => rho (:,:,:)
    nb = self%nbin
    ncells = size(temp1)
    ntot = size(chi)
    chip(1:ntot) => chi
    etap(1:ntot) => eta

    do c=1, ncells
      lt = log(temp1(c))
      lr = log(rho1(c))
      it = min(max(searchsorted(self%ltemp, lt, outbounds_ok=.true.), 1), self%ntemp-1)
      ir = min(max(searchsorted(self%lrho , lr, outbounds_ok=.true.), 1), self%nrho -1)
      ht = self%ltemp(it+1) - self%ltemp(it)
      hr = self%lrho (ir+1) - self%lrho (ir)
      tt = (lt - self%ltemp(it)) / ht
      tr = (lr - self%lrho (ir)) / hr
      call hermite_coefs(tt, ht, pt, tc0, tc1, tcm0, tcm1)
      call hermite_coefs(tr, hr, pr, rc0, rc1, rcm0, rcm1)
      do f=1, nb
        ! collapse log-temp for the two rho nodes ir, ir+1 (chi)
        v0 = tc0*self%chi(f,it,ir  ) + tc1*self%chi(f,it+1,ir  ) &
           + tcm0*self%dchi_t(f,it,ir  ) + tcm1*self%dchi_t(f,it+1,ir  )
        v1 = tc0*self%chi(f,it,ir+1) + tc1*self%chi(f,it+1,ir+1) &
           + tcm0*self%dchi_t(f,it,ir+1) + tcm1*self%dchi_t(f,it+1,ir+1)
        g0 = (1.-pt)*self%dchi_r(f,it,ir  ) + pt*self%dchi_r(f,it+1,ir  )
        g1 = (1.-pt)*self%dchi_r(f,it,ir+1) + pt*self%dchi_r(f,it+1,ir+1)
        cval = rc0*v0 + rc1*v1 + rcm0*g0 + rcm1*g1
        ! same for src
        v0 = tc0*self%src(f,it,ir  ) + tc1*self%src(f,it+1,ir  ) &
           + tcm0*self%dsrc_t(f,it,ir  ) + tcm1*self%dsrc_t(f,it+1,ir  )
        v1 = tc0*self%src(f,it,ir+1) + tc1*self%src(f,it+1,ir+1) &
           + tcm0*self%dsrc_t(f,it,ir+1) + tcm1*self%dsrc_t(f,it+1,ir+1)
        g0 = (1.-pt)*self%dsrc_r(f,it,ir  ) + pt*self%dsrc_r(f,it+1,ir  )
        g1 = (1.-pt)*self%dsrc_r(f,it,ir+1) + pt*self%dsrc_r(f,it+1,ir+1)
        sval = rc0*v0 + rc1*v1 + rcm0*g0 + rcm1*g1
        chip(f + (c-1)*nb) = cval
        etap(f + (c-1)*nb) = sval * cval
      end do
    end do
  END SUBROUTINE interpolate

  !===============================================================================
  ! Precompute monotone node slopes of the 3D sampling table along wavelength
  ! (dl, irregular), temperature (dt, regular) and log ne (de, regular).
  !===============================================================================
  SUBROUTINE sample_build_derivatives (self)
    class(sample_optab_type) :: self
    real(4), dimension(:), allocatable :: lam4, taxis, eaxis
    integer :: i
    !---------------------------------------------------------------------------
    allocate(self%dl, self%dt, self%de, mold=self%val)
    call os%mem_register(storage_size(self%dl), 3*product(shape(self%dl)), 'opac_slopes')

    allocate(lam4(self%nnu), taxis(self%ntemp), eaxis(self%nne))
    lam4  = real(self%lam)
    taxis = [(real(self%pca%c0(1) + (i-1)*self%pca%dc(1)), i=1, self%ntemp)]
    eaxis = [(real(self%pca%c0(2) + (i-1)*self%pca%dc(2)), i=1, self%nne)]

    call node_slopes_dim(lam4 , self%val, 1, self%dl)
    call node_slopes_dim(taxis, self%val, 2, self%dt)
    call node_slopes_dim(eaxis, self%val, 3, self%de)
  END SUBROUTINE sample_build_derivatives

  !===============================================================================
  ! Sample the 3D table into back(fs:fe,x,y,z), multiplying by the local chi500.
  ! Separable tensor-product cubic Hermite: collapse wavelength -> temp -> log ne
  ! using the precomputed node slopes (mixed partials ignored). The regular
  ! temp/ne axes use direct floor indexing; the irregular wavelength axis uses a
  ! per-frequency searchsorted precompute.
  !===============================================================================
  SUBROUTINE sample_interpolate (self, temp, ne, chi500, back, lam_q, fs, fe)
    class(sample_optab_type) :: self
    integer, intent(in) :: fs, fe
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(4), dimension(:,:,:), intent(in) :: chi500
    real(4), dimension(fs:,:,:,:), intent(inout) :: back
    real(8), dimension(:), intent(in) :: lam_q
    integer, dimension(:), allocatable :: kl
    real(4), dimension(:), allocatable :: pl, wc0, wc1, wcm0, wcm1
    real(8) :: rd1, rd2, ft, fe_, d1, d2, p1, p2
    real(4) :: hl, tl, tt, te, pt, pe
    real(4) :: tc0, tc1, tcm0, tcm1, nc0, nc1, ncm0, ncm1
    real(4) :: v00, v10, v01, v11, gt00, gt10, gt01, gt11
    real(4) :: ge00, ge10, ge01, ge11, f0, f1, ge0, ge1, kappa
    real(4) :: dc1_4, dc2_4
    integer :: nx, ny, nz, x, y, z, f, k, k1, it, je
    !---------------------------------------------------------------------------
    nx = size(chi500,1); ny = size(chi500,2); nz = size(chi500,3)
    rd1   = 1.d0 / self%pca%dc(1)
    rd2   = 1.d0 / self%pca%dc(2)
    dc1_4 = real(self%pca%dc(1))
    dc2_4 = real(self%pca%dc(2))

    ! per-frequency wavelength bracket and Hermite value coefficients (once)
    allocate(kl(fs:fe), pl(fs:fe), wc0(fs:fe), wc1(fs:fe), wcm0(fs:fe), wcm1(fs:fe))
    do f=fs, fe
      k = min(max(searchsorted(self%lam, lam_q(f), outbounds_ok=.true.), 1), self%nnu-1)
      kl(f) = k
      hl = real(self%lam(k+1) - self%lam(k))
      tl = real((lam_q(f) - self%lam(k)) / (self%lam(k+1) - self%lam(k)))
      call hermite_coefs(tl, hl, pl(f), wc0(f), wc1(f), wcm0(f), wcm1(f))
    end do

    do z=1, nz
      do y=1, ny
        do x=1, nx
          ! rotate (lnT, lnNe) into PC coordinates, then index each regular PC axis
          d1 = log(temp(x,y,z)) - self%pca%mu(1)
          d2 = log(ne(x,y,z))   - self%pca%mu(2)
          p1 = self%pca%rot(1,1)*d1 + self%pca%rot(2,1)*d2
          p2 = self%pca%rot(1,2)*d1 + self%pca%rot(2,2)*d2
          ! PC axis 1: direct index + Hermite coefficients (per cell)
          ft = (p1 - self%pca%c0(1)) * rd1 + 1.d0
          it = min(max(floor(ft), 1), self%ntemp-1)
          tt = real(ft - it)
          call hermite_coefs(tt, dc1_4, pt, tc0, tc1, tcm0, tcm1)
          ! PC axis 2
          fe_ = (p2 - self%pca%c0(2)) * rd2 + 1.d0
          je = min(max(floor(fe_), 1), self%nne-1)
          te = real(fe_ - je)
          call hermite_coefs(te, dc2_4, pe, nc0, nc1, ncm0, ncm1)

          do f=fs, fe
            k  = kl(f)
            k1 = k + 1
            ! collapse wavelength at the 4 (temp,ne) corners
            v00 = wc0(f)*self%val(k,it  ,je  ) + wc1(f)*self%val(k1,it  ,je  ) &
                + wcm0(f)*self%dl(k,it  ,je  ) + wcm1(f)*self%dl(k1,it  ,je  )
            v10 = wc0(f)*self%val(k,it+1,je  ) + wc1(f)*self%val(k1,it+1,je  ) &
                + wcm0(f)*self%dl(k,it+1,je  ) + wcm1(f)*self%dl(k1,it+1,je  )
            v01 = wc0(f)*self%val(k,it  ,je+1) + wc1(f)*self%val(k1,it  ,je+1) &
                + wcm0(f)*self%dl(k,it  ,je+1) + wcm1(f)*self%dl(k1,it  ,je+1)
            v11 = wc0(f)*self%val(k,it+1,je+1) + wc1(f)*self%val(k1,it+1,je+1) &
                + wcm0(f)*self%dl(k,it+1,je+1) + wcm1(f)*self%dl(k1,it+1,je+1)
            ! carry temp-slope and ne-slope across the collapsed wavelength axis
            gt00 = (1.-pl(f))*self%dt(k,it  ,je  ) + pl(f)*self%dt(k1,it  ,je  )
            gt10 = (1.-pl(f))*self%dt(k,it+1,je  ) + pl(f)*self%dt(k1,it+1,je  )
            gt01 = (1.-pl(f))*self%dt(k,it  ,je+1) + pl(f)*self%dt(k1,it  ,je+1)
            gt11 = (1.-pl(f))*self%dt(k,it+1,je+1) + pl(f)*self%dt(k1,it+1,je+1)
            ge00 = (1.-pl(f))*self%de(k,it  ,je  ) + pl(f)*self%de(k1,it  ,je  )
            ge10 = (1.-pl(f))*self%de(k,it+1,je  ) + pl(f)*self%de(k1,it+1,je  )
            ge01 = (1.-pl(f))*self%de(k,it  ,je+1) + pl(f)*self%de(k1,it  ,je+1)
            ge11 = (1.-pl(f))*self%de(k,it+1,je+1) + pl(f)*self%de(k1,it+1,je+1)
            ! collapse temp for the two ne rows
            f0  = tc0*v00 + tc1*v10 + tcm0*gt00 + tcm1*gt10
            f1  = tc0*v01 + tc1*v11 + tcm0*gt01 + tcm1*gt11
            ge0 = (1.-pt)*ge00 + pt*ge10
            ge1 = (1.-pt)*ge01 + pt*ge11
            ! collapse log ne
            kappa = nc0*f0 + nc1*f1 + ncm0*ge0 + ncm1*ge1
            back(f,x,y,z) = kappa * chi500(x,y,z)
          end do
        end do
      end do
    end do
  END SUBROUTINE sample_interpolate

  !===============================================================================
  SUBROUTINE sample_dealloc (self)
    class(sample_optab_type) :: self
    !---------------------------------------------------------------------------
    if (allocated(self%dl)) &
      call os%mem_register(storage_size(self%dl), -3*product(shape(self%dl)), 'opac_slopes')
    if (allocated(self%val)) call os%mem_register(storage_size(self%val), -product(shape(self%val)), 'opac_table')
    if (allocated(self%val)) deallocate(self%val)
    if (allocated(self%lam)) deallocate(self%lam)
    if (allocated(self%dl))  deallocate(self%dl, self%dt, self%de)
    self%initialised = .false.
  END SUBROUTINE sample_dealloc

  !===============================================================================
  ! Monotone PCHIP node slopes of a 3D array f along dimension dim, given the
  ! node coordinates axis. Reuses pchip_mono evaluated at the nodes: at a node
  ! the interpolant derivative equals the node slope. Works for regular and
  ! irregular axes (the Fritsch-Carlson weights collapse to the plain harmonic
  ! mean on a uniform grid automatically).
  !===============================================================================
  SUBROUTINE node_slopes_dim (axis, f, dim, slope)
    real(4), dimension(:), intent(in) :: axis
    real(4), dimension(:,:,:), intent(in) :: f
    integer, intent(in) :: dim
    real(4), dimension(:,:,:), intent(out) :: slope
    real(4), dimension(:), allocatable :: scratch
    integer :: n, a, b
    !---------------------------------------------------------------------------
    n = size(axis)
    allocate(scratch(n))
    select case (dim)
    case (1)
      do b=1, size(f,3)
        do a=1, size(f,2)
          call pchip_mono(axis, f(:,a,b), axis, scratch, sorted=.true., newdydx=slope(:,a,b))
        end do
      end do
    case (2)
      do b=1, size(f,3)
        do a=1, size(f,1)
          call pchip_mono(axis, f(a,:,b), axis, scratch, sorted=.true., newdydx=slope(a,:,b))
        end do
      end do
    case (3)
      do b=1, size(f,2)
        do a=1, size(f,1)
          call pchip_mono(axis, f(a,b,:), axis, scratch, sorted=.true., newdydx=slope(a,b,:))
        end do
      end do
    end select
  END SUBROUTINE node_slopes_dim

  !===============================================================================
  ! Build the node grid for the 3D sampling table (sample_optab_type) in the
  ! orthogonal PCA basis of (lnT, lnNe) over all cells, rather than an
  ! axis-aligned box in (T, lnNe). Since T and Ne are strongly correlated the
  ! atmosphere hugs a thin locus; a PC-rotated grid tracks it, dropping the
  ! unphysical off-diagonal corners that crash the EOS and cutting the node
  ! count. ldtemp/ldne are the requested lnT/lnNe resolutions, projected onto
  ! the PC axes for spacing. Callers pass temp/ne as read (real4, any shape)
  ! and get back ready-to-use (T, ne) node tables plus the transform needed
  ! later to rotate a query point into table coordinates (see sample_interpolate).
  !===============================================================================
  SUBROUTINE pca_opac_grid(temp, ne, ldtemp, ldne, npad, tab_temp, tab_ne, pca, n1, n2)
    real(4), dimension(:,:,:,:), intent(in) :: temp, ne
    real(8), intent(in) :: ldtemp, ldne
    integer, intent(in) :: npad
    real(8), dimension(:,:,:), allocatable, intent(out) :: tab_temp, tab_ne
    type(pca_transform_type), intent(out) :: pca
    integer, intent(out) :: n1, n2
    real(8), dimension(:), allocatable :: lnT, lnNe
    real(8) :: p1, p2
    integer :: i, j, ncell
    !---------------------------------------------------------------------------
    ! ---- lnT, lnNe computed once and reused for the whole PCA fit -------------
    ncell = size(temp)
    allocate(lnT(ncell), lnNe(ncell))
    lnT  = reshape(log(dble(temp)), [ncell])
    lnNe = reshape(log(dble(ne)),   [ncell])

    ! ---- PCA-rotated grid over (lnT, lnNe), npad nodes past the data hull ------
    call pca_grid_2d(lnT, lnNe, ldtemp, ldne, npad, pca%mu, pca%rot, pca%c0, pca%dc, n1, n2)
    deallocate(lnT, lnNe)

    ! ---- build the node grid in PC space, map back to (T, ne) -----------------
    allocate(tab_temp(n1, n2, 1), tab_ne(n1, n2, 1))
    do j=1, n2
    do i=1, n1
      p1 = pca%c0(1) + (i-1)*pca%dc(1)
      p2 = pca%c0(2) + (j-1)*pca%dc(2)
      tab_temp(i,j,1) = exp(pca%mu(1) + pca%rot(1,1)*p1 + pca%rot(1,2)*p2)   ! exp(lnT)
      tab_ne  (i,j,1) = exp(pca%mu(2) + pca%rot(2,1)*p1 + pca%rot(2,2)*p2)   ! exp(lnNe)
    end do
    end do
  END SUBROUTINE pca_opac_grid

END MODULE optab_mod
