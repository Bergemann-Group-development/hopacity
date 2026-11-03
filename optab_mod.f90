MODULE optab_mod
  USE hop_env_mod, only: os, trace
  USE console_mod
  USE prec_mod
  USE const_mod
  USE spectrum_mod, only: spectrum_type
  USE op_math_mod, only: searchsorted, pca_grid_2d
  USE interpol_mod, only: hermite_coefs, pchip_mono, pchip_slopes
  implicit none

  !=============================================================================
  !> Externally-read binned opacity table (special case): 2D (log temp, log rho)
  !> interpolation, one independent lookup per wavelength bin. Node derivatives
  !> along both axes are precomputed once in build_derivatives (called from init)
  !> so interpolate is a pure Hermite evaluation with stored slopes.
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
  !> Orthogonal PCA coordinate transform of (lnT, lnNe):
  !>   c = rot^T ((lnT,lnNe) - mu),   coord_k(i) = c0(k) + (i-1)*dc(k)
  !> Ties a sampled opacity table's (temp, ne) node grid to the physical
  !> coordinates it was built from, so a query point can be rotated into the
  !> same node-index space (see sample_interpolate). Built once, in
  !> pca_opac_grid, and carried from the composition patch that owns the grid
  !> down to every rt_patch that samples it (composition_type%pca).
  !=============================================================================
  type :: pca_transform_type
    real(8) :: mu(2)=0            !< mean of (lnT, lnNe)
    real(8) :: rot(2,2)=0          !< PC directions (unit eigenvectors as columns)
    real(8) :: c0(2)=0, dc(2)=1    !< min coordinate & spacing along each PC axis
  end type

  !=============================================================================
  !> Runtime background-opacity sampling table: 3D (wavelength, temp, log ne).
  !> The temp axis is regularly spaced (linear) and the ne axis is regularly
  !> spaced in log(ne); the wavelength axis is irregular. Monotone PCHIP node
  !> slopes along every dimension are precomputed and stored, so interpolate is
  !> an O(1) separable tensor-product cubic Hermite evaluation per query.
  !=============================================================================
  type :: sample_optab_type
    integer :: nnu=0, ntemp=0, nne=0
    type(pca_transform_type) :: pca       !< (lnT, lnNe) -> table-axis coordinates
    real(8), dimension(:), allocatable :: lam                 !< (nnu) irregular wavelength nodes
    real(4), dimension(:,:,:), allocatable :: val             !< (nnu,ntemp,nne) /chi500
    real(4), dimension(:,:,:), allocatable :: dl, dt, de      !< precomputed node slopes
    ! Second, scattering component with its own node slopes, /chi500.
    !   background table: val is the true absorption and scat the coherent
    !     scattering, kept DISJOINT and summed only after interpolation.  Both
    !     interpolants are non-negative (pchip_mono keeps each inside its local
    !     node range), so scat <= val + scat = chi_total holds by construction --
    !     that is what makes the interpolated single-scattering albedo <= 1
    !     without a clamp.
    !   line-list table: scat is the (1-eps)-weighted part of the line opacity
    !     already contained in val, because the foreground redirects that
    !     fraction of the line source from Bnu to Jnu without altering chi.
    ! Allocated on the background path always, on the line path only when line
    ! scattering is active.
    real(4), dimension(:,:,:), allocatable :: scat
    real(4), dimension(:,:,:), allocatable :: sl, st, se
    logical :: initialised=.false.
  contains
    procedure :: build_derivatives => sample_build_derivatives
    procedure :: interpolate => sample_interpolate
    procedure :: interpolate_ll => sample_interpolate_ll
    procedure :: dealloc => sample_dealloc
  end type

  type(optab_type), public, save :: optab

contains

  !===============================================================================
  !> initialises opacity table.
  !===============================================================================
  SUBROUTINE init (self, opac_table, sp)
    class(optab_type) :: self
    class(spectrum_type) :: sp
    character(len=*), intent(in) :: opac_table
    integer :: i, lunit
    !---------------------------------------------------------------------------
    if (console_on(ch_hop, 1)) write(console_unit,'(A,A)') 'Reading opacity table: ', trim(opac_table)

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
  !> Precompute monotone node slopes of the binned table along the (irregular)
  !> log-temperature and log-density axes, for both chi and src.
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
  !> Interpolate the binned table at the requested (temp, rho). Separable 2D
  !> tensor-product cubic Hermite (collapse log-temp then log-rho) using the
  !> precomputed node slopes; mixed partials are ignored. chi and eta are
  !> shaped (nbin, ncells...) with the bin as the first, contiguous index.
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
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('optab_mod%interpolate', itimer=itimer)
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
    call trace%end(itimer)
  END SUBROUTINE interpolate

  !===============================================================================
  !> Precompute monotone node slopes of the 3D sampling table along wavelength
  !> (dl, irregular), temperature (dt, regular) and log ne (de, regular).
  !> The wavelength slopes are taken on the real(8) nodes: the spectrum keeps its
  !> nodes lam_tolerance apart, but in real(4) two nodes 0.002 A apart at 3.8 um
  !> (the continua of two nearly degenerate levels) coincide. The PC axes are
  !> centred on the data and spaced dc apart, which real(4) resolves.
  !===============================================================================
  SUBROUTINE sample_build_derivatives (self)
    class(sample_optab_type) :: self
    real(4), dimension(:), allocatable :: taxis, eaxis
    integer :: i
    !---------------------------------------------------------------------------
    allocate(self%dl, self%dt, self%de, mold=self%val)
    call os%mem_register(storage_size(self%dl), 3*product(shape(self%dl)), 'opac_slopes')

    allocate(taxis(self%ntemp), eaxis(self%nne))
    taxis = [(real(self%pca%c0(1) + (i-1)*self%pca%dc(1)), i=1, self%ntemp)]
    eaxis = [(real(self%pca%c0(2) + (i-1)*self%pca%dc(2)), i=1, self%nne)]

    call lam_slopes(self%val, self%dl)
    call node_slopes_dim(taxis, self%val, 2, self%dt)
    call node_slopes_dim(eaxis, self%val, 3, self%de)

    if (allocated(self%scat)) then
      allocate(self%sl, self%st, self%se, mold=self%scat)
      call os%mem_register(storage_size(self%sl), 3*product(shape(self%sl)), 'opac_slopes')
      call lam_slopes(self%scat, self%sl)
      call node_slopes_dim(taxis, self%scat, 2, self%st)
      call node_slopes_dim(eaxis, self%scat, 3, self%se)
    endif
  contains
    subroutine lam_slopes (f, slope)
      real(4), dimension(:,:,:), intent(in) :: f
      real(4), dimension(:,:,:), intent(out) :: slope
      real(8), dimension(:), allocatable :: f8, dx, dy, dydx
      integer :: n, a, b
      n = self%nnu
      allocate(f8(n), dx(n-1), dy(n-1), dydx(n))
      do b=1, size(f,3)
        do a=1, size(f,2)
          f8 = f(:,a,b)
          call pchip_slopes(self%lam, f8, 1, n, dx, dy, dydx)
          slope(:,a,b) = real(dydx)
        end do
      end do
    end subroutine lam_slopes
  END SUBROUTINE sample_build_derivatives

  !===============================================================================
  !> Sample the 3D table into back(fs:fe,x,y,z), multiplying by the local chi500.
  !> Separable tensor-product cubic Hermite: collapse wavelength -> temp -> log ne
  !> using the precomputed node slopes (mixed partials ignored). The regular
  !> temp/ne axes use direct floor indexing; the irregular wavelength axis uses a
  !> per-frequency searchsorted precompute.
  !===============================================================================
  SUBROUTINE sample_interpolate (self, temp, ne, chi500, back, scat, lam_q, fs, fe)
    class(sample_optab_type) :: self
    integer, intent(in) :: fs, fe
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(4), dimension(:,:,:), intent(in) :: chi500
    ! back receives the TOTAL extinction and scat its scattering part.  The two
    ! table components are disjoint (true absorption / scattering), interpolated
    ! separately with identical weights and summed here -- so scat <= back holds
    ! pointwise, and the single-scattering albedo the transfer sees can never
    ! exceed 1 no matter how coarse the table is.
    real(4), dimension(fs:,:,:,:), intent(inout) :: back, scat
    real(8), dimension(:), intent(in) :: lam_q
    integer, dimension(:), allocatable :: kl
    real(4), dimension(:), allocatable :: pl, wc0, wc1, wcm0, wcm1
    real(8) :: rd1, rd2, ft, fe_, d1, d2, p1, p2
    real(4) :: hl, tl, tt, te, pt, pe
    real(4) :: tc0, tc1, tcm0, tcm1, nc0, nc1, ncm0, ncm1
    real(4) :: kappa, kappa_scat
    real(4) :: dc1_4, dc2_4
    integer :: nx, ny, nz, x, y, z, f, k, k1, it, je
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('optab_mod%sample_interpolate', itimer=itimer)
    ! the wavelength bracket below needs two nodes; one would be read out of bounds
    if (self%nnu < 2) error stop 'sample_interpolate: a sampled opacity table needs at least two wavelengths'
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
            kappa      = collapse(self%val , self%dl, self%dt, self%de)
            kappa_scat = collapse(self%scat, self%sl, self%st, self%se)
            back(f,x,y,z) = (kappa + kappa_scat) * chi500(x,y,z)
            scat(f,x,y,z) =          kappa_scat  * chi500(x,y,z)
          end do
        end do
      end do
    end do
    call trace%end(itimer)

  contains
    !> separable tensor-product cubic Hermite collapse (wavelength -> temp ->
    !> log ne) at the corner (k,it,je), reusing the host's precomputed Hermite
    !> coefficients; one call per table component
    real(4) function collapse(val, dl, dt, de) result(kap)
      real(4), dimension(:,:,:), intent(in) :: val, dl, dt, de
      real(4) :: v00, v10, v01, v11, gt00, gt10, gt01, gt11
      real(4) :: ge00, ge10, ge01, ge11, f0, f1, ge0, ge1
      ! collapse wavelength at the 4 (temp,ne) corners
      v00 = wc0(f)*val(k,it  ,je  ) + wc1(f)*val(k1,it  ,je  ) &
          + wcm0(f)*dl(k,it  ,je  ) + wcm1(f)*dl(k1,it  ,je  )
      v10 = wc0(f)*val(k,it+1,je  ) + wc1(f)*val(k1,it+1,je  ) &
          + wcm0(f)*dl(k,it+1,je  ) + wcm1(f)*dl(k1,it+1,je  )
      v01 = wc0(f)*val(k,it  ,je+1) + wc1(f)*val(k1,it  ,je+1) &
          + wcm0(f)*dl(k,it  ,je+1) + wcm1(f)*dl(k1,it  ,je+1)
      v11 = wc0(f)*val(k,it+1,je+1) + wc1(f)*val(k1,it+1,je+1) &
          + wcm0(f)*dl(k,it+1,je+1) + wcm1(f)*dl(k1,it+1,je+1)
      ! carry temp-slope and ne-slope across the collapsed wavelength axis
      gt00 = (1.-pl(f))*dt(k,it  ,je  ) + pl(f)*dt(k1,it  ,je  )
      gt10 = (1.-pl(f))*dt(k,it+1,je  ) + pl(f)*dt(k1,it+1,je  )
      gt01 = (1.-pl(f))*dt(k,it  ,je+1) + pl(f)*dt(k1,it  ,je+1)
      gt11 = (1.-pl(f))*dt(k,it+1,je+1) + pl(f)*dt(k1,it+1,je+1)
      ge00 = (1.-pl(f))*de(k,it  ,je  ) + pl(f)*de(k1,it  ,je  )
      ge10 = (1.-pl(f))*de(k,it+1,je  ) + pl(f)*de(k1,it+1,je  )
      ge01 = (1.-pl(f))*de(k,it  ,je+1) + pl(f)*de(k1,it  ,je+1)
      ge11 = (1.-pl(f))*de(k,it+1,je+1) + pl(f)*de(k1,it+1,je+1)
      ! collapse temp for the two ne rows, then log ne
      f0  = tc0*v00 + tc1*v10 + tcm0*gt00 + tcm1*gt10
      f1  = tc0*v01 + tc1*v11 + tcm0*gt01 + tcm1*gt11
      ge0 = (1.-pt)*ge00 + pt*ge10
      ge1 = (1.-pt)*ge01 + pt*ge11
      kap = nc0*f0 + nc1*f1 + ncm0*ge0 + ncm1*ge1
    end function collapse
  END SUBROUTINE sample_interpolate

  !===============================================================================
  !> Sample the LTE line-list table into the foreground work arrays, ADDING the
  !> interpolated value to both xhat (work(...,1)) and ehat (work(...,2)) -- in LTE
  !> the line's two working forms are identical (val stores chi_line/((1-beta)*chi500)).
  !> The query wavelength is Doppler-shifted per cell by the projected velocity:
  !> lam_query = lam_q(f)*(1 - vproj) (equivalently log-lam - vproj, the same sign
  !> convention as the line kernels' llc = llam0 - vproj_c). Because the shift is
  !> per cell, the wavelength bracket is computed inside the (x,y,z) loop rather
  !> than precomputed per frequency. Separable tensor-product cubic Hermite over
  !> (wavelength -> temp -> log ne), reusing the linear-lambda node slopes dl.
  !===============================================================================
  SUBROUTINE sample_interpolate_ll (self, temp, ne, chi500, work, lam_q, vproj, fs, fe, chi_scat)
    class(sample_optab_type) :: self
    integer, intent(in) :: fs, fe
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(4), dimension(:,:,:), intent(in) :: chi500, vproj
    real(4), dimension(fs:,:,:,:,:), intent(inout) :: work
    real(8), dimension(:), intent(in) :: lam_q
    ! line_scattering: accumulate the interpolated scattering component of the
    ! line opacity here (only when the table carries it)
    real(4), dimension(fs:,:,:,:), intent(inout), optional :: chi_scat
    real(8) :: rd1, rd2, ft, fe_, d1, d2, p1, p2, lamq, vf
    real(4) :: hl, tl, pl, tt, te, pt, pe
    real(4) :: wc0, wc1, wcm0, wcm1
    real(4) :: tc0, tc1, tcm0, tcm1, nc0, nc1, ncm0, ncm1
    real(4) :: kappa, add
    real(4) :: dc1_4, dc2_4
    integer :: nx, ny, nz, x, y, z, f, k, k1, it, je
    logical :: do_scat
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('optab_mod%sample_interpolate_ll', itimer=itimer)
    ! the wavelength bracket below needs two nodes; one would be read out of bounds
    if (self%nnu < 2) error stop 'sample_interpolate_ll: a sampled opacity table needs at least two wavelengths'
    nx = size(chi500,1); ny = size(chi500,2); nz = size(chi500,3)
    do_scat = present(chi_scat) .and. allocated(self%scat)
    rd1   = 1.d0 / self%pca%dc(1)
    rd2   = 1.d0 / self%pca%dc(2)
    dc1_4 = real(self%pca%dc(1))
    dc2_4 = real(self%pca%dc(2))

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

          ! Per-cell Doppler shift. The detailed line kernel places a line at
          ! observed log-lambda llam0 - vproj (llc = llam0 - vproj), so the rest
          ! table (built at vproj=0, peak at llam0) must be queried at
          ! llam(f) + vproj, i.e. lambda*(1 + vproj) to first order.
          vf = 1.d0 + real(vproj(x,y,z), 8)

          do f=fs, fe
            lamq = lam_q(f) * vf
            k  = min(max(searchsorted(self%lam, lamq, outbounds_ok=.true.), 1), self%nnu-1)
            k1 = k + 1
            hl = real(self%lam(k1) - self%lam(k))
            tl = real((lamq - self%lam(k)) / (self%lam(k1) - self%lam(k)))
            call hermite_coefs(tl, hl, pl, wc0, wc1, wcm0, wcm1)
            kappa = collapse(self%val, self%dl, self%dt, self%de)
            add = kappa * chi500(x,y,z)
            work(f,x,y,z,1) = work(f,x,y,z,1) + add
            work(f,x,y,z,2) = work(f,x,y,z,2) + add
            if (do_scat) then
              kappa = collapse(self%scat, self%sl, self%st, self%se)
              chi_scat(f,x,y,z) = chi_scat(f,x,y,z) + kappa * chi500(x,y,z)
            endif
          end do
        end do
      end do
    end do
    call trace%end(itimer)

  contains
    !> separable tensor-product cubic Hermite collapse (wavelength -> temp ->
    !> log ne) at the corner (k,it,je), reusing the host's precomputed Hermite
    !> coefficients; one call per table component
    real(4) function collapse(val, dl, dt, de) result(kap)
      real(4), dimension(:,:,:), intent(in) :: val, dl, dt, de
      real(4) :: cv00, cv10, cv01, cv11, cgt00, cgt10, cgt01, cgt11
      real(4) :: cge00, cge10, cge01, cge11, cf0, cf1, cge0, cge1
      cv00 = wc0*val(k,it  ,je  ) + wc1*val(k1,it  ,je  ) &
           + wcm0*dl(k,it  ,je  ) + wcm1*dl(k1,it  ,je  )
      cv10 = wc0*val(k,it+1,je  ) + wc1*val(k1,it+1,je  ) &
           + wcm0*dl(k,it+1,je  ) + wcm1*dl(k1,it+1,je  )
      cv01 = wc0*val(k,it  ,je+1) + wc1*val(k1,it  ,je+1) &
           + wcm0*dl(k,it  ,je+1) + wcm1*dl(k1,it  ,je+1)
      cv11 = wc0*val(k,it+1,je+1) + wc1*val(k1,it+1,je+1) &
           + wcm0*dl(k,it+1,je+1) + wcm1*dl(k1,it+1,je+1)
      cgt00 = (1.-pl)*dt(k,it  ,je  ) + pl*dt(k1,it  ,je  )
      cgt10 = (1.-pl)*dt(k,it+1,je  ) + pl*dt(k1,it+1,je  )
      cgt01 = (1.-pl)*dt(k,it  ,je+1) + pl*dt(k1,it  ,je+1)
      cgt11 = (1.-pl)*dt(k,it+1,je+1) + pl*dt(k1,it+1,je+1)
      cge00 = (1.-pl)*de(k,it  ,je  ) + pl*de(k1,it  ,je  )
      cge10 = (1.-pl)*de(k,it+1,je  ) + pl*de(k1,it+1,je  )
      cge01 = (1.-pl)*de(k,it  ,je+1) + pl*de(k1,it  ,je+1)
      cge11 = (1.-pl)*de(k,it+1,je+1) + pl*de(k1,it+1,je+1)
      cf0  = tc0*cv00 + tc1*cv10 + tcm0*cgt00 + tcm1*cgt10
      cf1  = tc0*cv01 + tc1*cv11 + tcm0*cgt01 + tcm1*cgt11
      cge0 = (1.-pt)*cge00 + pt*cge10
      cge1 = (1.-pt)*cge01 + pt*cge11
      kap = nc0*cf0 + nc1*cf1 + ncm0*cge0 + ncm1*cge1
    end function collapse
  END SUBROUTINE sample_interpolate_ll

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
    if (allocated(self%sl)) &
      call os%mem_register(storage_size(self%sl), -3*product(shape(self%sl)), 'opac_slopes')
    if (allocated(self%scat)) call os%mem_register(storage_size(self%scat), -product(shape(self%scat)), 'opac_table')
    if (allocated(self%scat)) deallocate(self%scat)
    if (allocated(self%sl))   deallocate(self%sl, self%st, self%se)
    self%initialised = .false.
  END SUBROUTINE sample_dealloc

  !===============================================================================
  !> Monotone PCHIP node slopes of a 3D array f along dimension dim, given the
  !> node coordinates axis. Reuses pchip_mono evaluated at the nodes: at a node
  !> the interpolant derivative equals the node slope. Works for regular and
  !> irregular axes (the Fritsch-Carlson weights collapse to the plain harmonic
  !> mean on a uniform grid automatically).
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
  !> Build the node grid for the 3D sampling table (sample_optab_type) in the
  !> orthogonal PCA basis of (lnT, lnNe) over all cells, rather than an
  !> axis-aligned box in (T, lnNe). Since T and Ne are strongly correlated the
  !> atmosphere hugs a thin locus; a PC-rotated grid tracks it, dropping the
  !> unphysical off-diagonal corners that crash the EOS and cutting the node
  !> count. ldtemp/ldne are the requested lnT/lnNe resolutions, projected onto
  !> the PC axes for spacing. Callers pass temp/ne as read (real4, any shape)
  !> and get back ready-to-use (T, ne) node tables plus the transform needed
  !> later to rotate a query point into table coordinates (see sample_interpolate).
  !===============================================================================
  SUBROUTINE pca_opac_grid(temp, ne, ldtemp, ldne, npad, tab_temp, tab_ne, pca, n1, n2, use_pca)
    real(4), dimension(:,:,:,:), intent(in) :: temp, ne
    real(8), intent(in) :: ldtemp, ldne
    integer, intent(in) :: npad
    logical, intent(in), optional :: use_pca   !< .false. -> axis-aligned (lnT,lnNe) grid
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
    ! (use_pca=.false. bypasses the rotation for an axis-aligned diagnostic grid)
    call pca_grid_2d(lnT, lnNe, ldtemp, ldne, npad, pca%mu, pca%rot, pca%c0, pca%dc, n1, n2, use_pca)
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
