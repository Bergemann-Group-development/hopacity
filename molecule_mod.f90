MODULE molecule_mod
  USE os_mod
  USE trace_mod
  USE prec_mod
  USE op_math_mod
  USE ion_mod
  USE handyfuncs_mod
  implicit none

  type, extends(particle_type) :: molec_type
    integer :: nions
    class(ion_pointer_t), allocatable :: ions(:)
    real(dp) :: dc(5) ! dissociation_constants
    real(dp), dimension(:), allocatable :: de ! dissociation_energies
    real(dp), dimension(:,:,:), allocatable :: dis
  contains
    procedure :: read_molec_data
    procedure :: get_dis
    procedure :: calc_f
    procedure :: calc_e
    procedure :: allocate => allocate_molec
    procedure :: deallocate => deallocate_molec
    procedure :: generate_isotopic_combinations
  end type

  logical, public :: irwin_dis=.true.

contains

  !===============================================================================
  SUBROUTINE read_molec_data(self, lunit, iostat)
    class(molec_type), intent(inout) :: self
    integer, intent(in) :: lunit
    integer, intent(out) :: iostat
    integer :: defix, n_irwin, n_de
    character(64) :: text
    !---------------------------------------------------------------------------
    read(lunit, *, iostat=iostat) text
    if (iostat < 0) return
    do while (text(1:1) == '!')
      read(lunit, *, iostat=iostat) text
    end do
    backspace lunit
    if (allocated(self%de)) deallocate(self%de)
    if (allocated(self%irwin)) deallocate(self%irwin)
    read(lunit, *, iostat=iostat) self%id, defix, n_irwin
    n_de = merge(6, 1, defix==0)
    allocate(self%de(n_de))
    allocate(self%irwin(n_irwin))
    read(lunit, *) self%de
    read(lunit, *) self%dc ! always 5 numbers
    read(lunit, *) self%irwin
    self%irwin(:) = self%irwin(n_irwin:1:-1)
  END SUBROUTINE read_molec_data

  !===============================================================================
  ! Computes the equilibrium constant log(P(AB)/(P(A)*P(B))) for molecules
  !===============================================================================
  SUBROUTINE get_dis(self, temp, teta, e_kT, dis)
    class(molec_type), intent(in) :: self
    real(dp), intent(in), dimension(:,:,:) :: temp, teta, e_kT
    real(dp), intent(out), dimension(:,:,:) :: dis
    real(dp) :: awgt15, scalar, N_tot, const_term
    real(dp), parameter :: const_base = log((2*PI*UU/HH**2)**(-1.5) * BK**(-2.5))
    integer :: j, nm1
    !---------------------------------------------------------------------------
    if (irwin_dis .or. all(self%dc == 0)) then
      nm1 = self%nions - 1
      const_term = const_base*nm1

      dis = 0
      do j=1, self%nions
        associate(ion => self%ions(j)%p)
          awgt15 = -1.5*log(ion%awgt)
          dis = dis + awgt15 - ion%part
        end associate
      end do

      scalar = const_term + 1.5*log(self%awgt)
      dis = dis + scalar
      dis = dis + self%part + self%de(1)*e_kT - (2.5*nm1)*log(temp)
    else
      call poly3d(teta, self%dc, dis)
      dis = -ln10*dis
    endif
  END SUBROUTINE get_dis

  !===============================================================================
  FUNCTION calc_f(self, fepe) result(f)
    class(molec_type), intent(in) :: self
    real(dp), dimension(:,:,:) :: fepe
    real(dp), dimension(:,:,:), allocatable :: f
    type(ion_type), pointer :: ion
    integer :: j, x, y, z
    !---------------------------------------------------------------------------
    allocate(f, source=self%dis)
    do j=1, self%nions
      ion => self%ions(j)%p
      where (ion%f > 0)
        f = f + log(ion%f)
      elsewhere
        f = - huge(1.d0)
      endwhere
      ! fepe needed to convert from pressure to f, first fepe cancels out
      if (j > 1) f = f - fepe
    end do
    f = min(f, 0.d0)
    f = exp(f)
  END FUNCTION calc_f

  !===============================================================================
  FUNCTION calc_e(self, teta, temp, dxi) result(e)
    class(molec_type), intent(in) :: self
    class(ion_type), pointer :: ion
    real(dp), dimension(:,:,:) :: teta, temp, dxi
    real(dp), dimension(:,:,:), allocatable :: de, e
    integer :: i
    !---------------------------------------------------------------------------
    allocate(de, e, mold=self%n)

    if (size(self%de) > 1) then
      ! rotational and vibrational energy has to be subtracted
      call poly3d(teta, self%de(2:6), de)
      de = self%de(1) - de*temp/EK
    else
      de = self%de(1)
    endif

    do i=1, self%nions
      ion => self%ions(i)%p

      ! why only consider ionisation energy from H?
      if (trim(ion%id) == 'H I') then
        de = de + ion%xion - dxi
      endif
    end do

    e = - de * self%n
  END FUNCTION calc_e

  !=========================================================================
  ! Molecular + ionization equilibrium, solved cell by cell.
  !
  ! The caller supplies the molecules and the solve elements: one neutral
  ! "root" ion per molecule-forming element (ions).  The routine runs at
  ! fixed electron density, so within one cell the Saha ratios between the
  ! ionization stages of every element, the electron pressure, and the
  ! dissociation constants are all constants.  With yi = ln f_i for the
  ! neutral fraction of each solve element and ye = ln fe, every molecular
  ! fraction is an explicit function of the unknowns,
  !
  !   ln f_mol = C_m + sum_i nu_im*yi - (n_m - 1)*ye ,
  !
  ! where the per-cell constant C_m absorbs the dissociation constant, the
  ! electron pressure, and the Saha factors of charged constituents (H2+).
  ! The nonlinear system is number conservation for each solve element plus
  ! charge conservation,
  !
  !   F_i = f_i*phi_i + sum_m nu_im*f_mol - A_i           = 0
  !   F_e = xnecno + sum_i q_i*f_i + sum_m q_m*f_mol - fe = 0 ,
  !
  ! with phi_i (total nuclei per neutral) and q_i (electrons per neutral)
  ! given by the fixed Saha ladder of element i.  It is solved by a damped
  ! Newton method in log space: analytic Jacobian, step capped in ln f,
  ! Armijo backtracking on |F|^2.  One iteration costs nions+nmol+1
  ! exponentials plus an (nions+1)^2 linear solve; there are no logarithms,
  ! pointer walks, 3D-array accesses, or array-syntax statements inside the
  ! iteration, and the working set is O((nions+1)^2) scalars.
  !=========================================================================
  SUBROUTINE molecular_equilibrium(molecules, ions, tmolim, temp, teta, e_kT, &
                                   pe, xnecno, fe)
    class(molec_type), target, intent(inout) :: molecules(:)
    type(ion_pointer_t), intent(in) :: ions(:)
    real(dp), intent(in) :: tmolim
    real(dp), dimension(:,:,:), intent(in) :: temp, teta, e_kT, pe, xnecno
    real(dp), dimension(:,:,:), intent(inout) :: fe
    class(molec_type), pointer :: mol
    type(ion_type), pointer :: pt
    integer :: i, j, k, m, iq, nx, ny, nz, np1, nions, nmol, mxc
    integer :: x, y, z, iter, nfail, iH, iH2, qlo, qhi
    ! fixed stoichiometry: per-constituent root index / charge offset, and
    ! per-molecule lists of distinct solve elements with their multiplicity
    integer, dimension(:),   allocatable :: ncons, nele, nm1, qmol
    integer, dimension(:,:), allocatable :: cidx, cq, eidx, enu
    ! per-cell constants
    real(dp), dimension(:,:), allocatable :: rat        ! stage/neutral ratios
    real(dp), dimension(:),   allocatable :: fil, ce    ! phi_i, q_i
    real(dp), dimension(:),   allocatable :: abund, Cm
    real(dp) :: pe_c, lnpe, fe_c, xne_c, b, s, fH, amax
    ! Newton state
    real(dp), dimension(:),   allocatable :: yv, yt, Fv, Ft, stp
    real(dp), dimension(:),   allocatable :: fi, fit, fm, fmt
    real(dp), dimension(:,:), allocatable :: Amat
    real(dp) :: norm2, norm2t, errf, lam, smax, worst
    logical  :: converged
    integer,  parameter :: niter_max = 100
    real(dp), parameter :: tol       = 1d-10      ! scaled residual target
    real(dp), parameter :: tol_stall = 1d-6       ! acceptable when stalled
    real(dp), parameter :: step_max  = 3d0        ! cap on |d ln f| per step
    real(dp), parameter :: lam_min   = 1d0/1024d0 ! smallest damping factor
    integer :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('molecule_mod::molecular_equilibrium', itimer=itimer)

    nx = size(temp, 1)
    ny = size(temp, 2)
    nz = size(temp, 3)
    nmol  = size(molecules)
    nions = size(ions)
    np1   = nions + 1

    ! molecule partition functions and equilibrium constants; the scalar
    ! irwin_partf call keeps the elemental expansion out of the hot files
    do m=1, nmol
      mol => molecules(m)
      call mol%allocate(shape(temp))
      do z=1, nz
        do y=1, ny
          do x=1, nx
            mol%part(x,y,z) = mol%irwin_partf(temp(x,y,z))
          end do
        end do
      end do
      call mol%get_dis(temp, teta, e_kT, mol%dis)
      ! molecules default to zero; only cells below tmolim get a solved value
      do z=1, nz
        do y=1, ny
          do x=1, nx
            mol%f(x,y,z) = 0
          end do
        end do
      end do
    end do

    ! charge extent of the ionization ladders (H- gives qlo = -1)
    qlo = 0
    qhi = 0
    do i=1, nions
      k = 0
      pt => ions(i)%p
      do while (associated(pt%above))
        k = k + 1
        pt => pt%above
      end do
      qhi = max(qhi, k)
      k = 0
      pt => ions(i)%p
      do while (associated(pt%below))
        k = k - 1
        pt => pt%below
      end do
      qlo = min(qlo, k)
    end do

    ! stoichiometry, resolved once: every constituent is mapped to its
    ! neutral root (solve index) and the charge it carries relative to it
    mxc = 0
    do m=1, nmol
      mxc = max(mxc, molecules(m)%nions)
    end do
    allocate(ncons(nmol), nele(nmol), nm1(nmol), qmol(nmol))
    allocate(cidx(mxc,nmol), cq(mxc,nmol), eidx(mxc,nmol), enu(mxc,nmol))
    do m=1, nmol
      mol => molecules(m)
      ncons(m) = mol%nions
      nm1(m)   = mol%nions - 1
      qmol(m)  = 0
      nele(m)  = 0
      do j=1, mol%nions
        call constituent_root(mol%ions(j)%p, i, iq)
        cidx(j,m) = i
        cq(j,m)   = iq
        qmol(m)   = qmol(m) + iq
        do k=1, nele(m)
          if (eidx(k,m) == i) exit
        end do
        if (k > nele(m)) then
          nele(m)   = k
          eidx(k,m) = i
          enu(k,m)  = 1
        else
          enu(k,m)  = enu(k,m) + 1
        endif
      end do
    end do

    ! H and H2 drive the initial guess (harmless to skip if absent)
    iH = 0
    iH2 = 0
    do i=1, nions
      if (trim(ions(i)%p%id) == 'H I') iH = i
    end do
    do m=1, nmol
      if (trim(molecules(m)%id) == 'H2') iH2 = m
    end do

    allocate(rat(qlo:qhi, nions), fil(nions), ce(nions), abund(nions), Cm(nmol))
    allocate(yv(np1), yt(np1), Fv(np1), Ft(np1), stp(np1))
    allocate(fi(nions), fit(nions), fm(nmol), fmt(nmol), Amat(np1,np1))

    do i=1, nions
      abund(i) = ions(i)%p%abund
    end do

    nfail = 0
    worst = 0

    do z=1, nz
      do y=1, ny
        do x=1, nx
          if (temp(x,y,z) >= tmolim) cycle ! molecules frozen out (=0)

          pe_c  = pe(x,y,z)
          xne_c = xnecno(x,y,z)
          ! fe(x,y,z) is a fresh Saha sum every call, but in very cold/thin
          ! cells an ionization-ladder ratio upstream can underflow/round to
          ! something non-positive; log() of that is NaN, which this routine
          ! writes straight back into fe(x,y,z), permanently poisoning every
          ! later outer iteration for the same cell. Floor it here so a
          ! single bad upstream sample can't cascade.
          fe_c  = max(fe(x,y,z), tiny(1d0))
          lnpe  = log(pe_c)

          ! --- fixed per-cell data --------------------------------------
          ! stage/neutral ratios from the Saha ladder, then
          ! phi_i = sum_q rat_q  and  q_i = sum_q q*rat_q
          do i=1, nions
            do k=qlo, qhi
              rat(k,i) = 0
            end do
            rat(0,i) = 1
            k = 0
            pt => ions(i)%p
            do while (associated(pt%above))
              rat(k+1,i) = rat(k,i) * pt%saha(x,y,z)
              k = k + 1
              pt => pt%above
            end do
            k = 0
            pt => ions(i)%p
            do while (associated(pt%below))
              rat(k-1,i) = rat(k,i) / pt%below%saha(x,y,z)
              k = k - 1
              pt => pt%below
            end do
            fil(i) = 0
            ce(i)  = 0
            do k=qlo, qhi
              fil(i) = fil(i) + rat(k,i)
              ce(i)  = ce(i)  + k*rat(k,i)
            end do
          end do

          ! ln K of each molecule in f-units, up to the -(n-1)*ye term
          do m=1, nmol
            b = molecules(m)%dis(x,y,z) + nm1(m)*lnpe
            do j=1, ncons(m)
              if (cq(j,m) /= 0) b = b + log(rat(cq(j,m), cidx(j,m)))
            end do
            Cm(m) = b
          end do

          ! --- initial guess ----------------------------------------------
          ! Saha neutral fractions, refined for H by solving the H/H2
          ! quadratic  phi_H*f_H + 2K*f_H^2 = A_H  (stable root), with the
          ! electron budget adjusted for the hydrogen thus removed
          do i=1, nions
            yv(i) = log(abund(i)/fil(i))
          end do
          yv(np1) = log(fe_c)
          if (iH > 0 .and. iH2 > 0) then
            b  = exp(min(Cm(iH2) - yv(np1), 600d0)) ! K = f_H2/f_H^2
            s  = sqrt(fil(iH)**2 + 8d0*b*abund(iH))
            fH = 2d0*abund(iH) / (fil(iH) + s)
            fe_c = fe_c + (fH - abund(iH)/fil(iH))*ce(iH)
            fe_c = max(fe_c, 1d-3*xne_c, tiny(1d0))
            yv(iH)  = log(fH)
            yv(np1) = log(fe_c)
          endif

          ! Beyond H2, the "assume everything stays atomic" guess can be
          ! wildly wrong in cold/thin cells where a molecule's equilibrium
          ! constant so overwhelmingly favors association that its implied
          ! ln(f_mol) is already deep positive (>>300) before Newton takes a
          ! single step: the resulting f_mol saturates at the exp() overflow
          ! guard, its Jacobian entries dwarf every physically-relevant one
          ! (which span only ~1e-12 to 1), and the linear solve's precision
          ! is destroyed for every other element from iteration one. A few
          ! Gauss-Seidel sweeps redistribute each over-saturated molecule's
          ! excess back onto its constituent elements' guesses (minimal-norm
          ! correction, weighted by stoichiometry), capping ln(f_mol) at what
          ! its scarcest constituent could actually supply -- not at 0 (i.e.
          ! f_mol<=1), which is no tighter than the old, already-rejected
          ! "f_mol<=1 by construction" cap and just as unphysical for trace
          ! elements -- so the starting point lands in a regime the damped
          ! Newton loop can actually refine.
          do iter=1, 5
            do m=1, nmol
              b = Cm(m) - nm1(m)*yv(np1)
              do k=1, nele(m)
                b = b + enu(k,m)*yv(eidx(k,m))
              end do
              amax = huge(1d0)
              do k=1, nele(m)
                amax = min(amax, log(abund(eidx(k,m))) - log(dble(enu(k,m))))
              end do
              if (b > amax) then
                s = 0
                do k=1, nele(m)
                  s = s + enu(k,m)**2
                end do
                do k=1, nele(m)
                  yv(eidx(k,m)) = yv(eidx(k,m)) - (b-amax)*enu(k,m)/s
                end do
              endif
            end do
          end do

          ! --- damped Newton to convergence -------------------------------
          call eval(yv, fi, fm, Fv, norm2)
          converged = .false.

          do iter=1, niter_max
            errf = abs(Fv(np1)) / max(fe_c, xne_c)
            do i=1, nions
              errf = max(errf, abs(Fv(i))/abund(i))
            end do
            if (errf < tol) then
              converged = .true.
              exit
            endif

            ! analytic Jacobian dF/d(ln f) at the current state
            do j=1, np1
              do i=1, np1
                Amat(i,j) = 0
              end do
            end do
            do i=1, nions
              Amat(i,i)   = fi(i)*fil(i)
              Amat(np1,i) = fi(i)*ce(i)
            end do
            Amat(np1,np1) = -fe_c
            do m=1, nmol
              do k=1, nele(m)
                i = eidx(k,m)
                Amat(i,np1) = Amat(i,np1) - (enu(k,m)*nm1(m))*fm(m)
                do j=1, nele(m)
                  Amat(i,eidx(j,m)) = Amat(i,eidx(j,m)) + (enu(k,m)*enu(j,m))*fm(m)
                end do
              end do
              if (qmol(m) /= 0) then
                Amat(np1,np1) = Amat(np1,np1) - (qmol(m)*nm1(m))*fm(m)
                do k=1, nele(m)
                  Amat(np1,eidx(k,m)) = Amat(np1,eidx(k,m)) + (qmol(m)*enu(k,m))*fm(m)
                end do
              endif
            end do

            ! Newton step, capped in ln f to stay in the basin. Also cap the
            ! implied step in each molecule's ln(f_mol) argument 'a': a
            ! molecule couples several elements (nu up to ~4), so its 'a' can
            ! swing far more per iteration than any single y component does,
            ! even though every component individually respects step_max.
            ! Left unchecked, that let a handful of very cold/thin cells
            ! (M-dwarf-like) run 'a' up past +300 in a couple of iterations,
            ! saturating f_mol at an astronomically unphysical value whose
            ! Jacobian entries then swamped the physically-relevant ones
            ! (spanning only ~1e-12 to 1) in the linear solve, wrecking its
            ! conditioning so badly that Newton stopped making any progress
            ! at all on the unrelated elements.
            do i=1, np1
              stp(i) = Fv(i)
            end do
            call lineq(Amat, stp, np1)
            smax = 0
            do i=1, np1
              smax = max(smax, abs(stp(i)))
            end do
            do m=1, nmol
              s = nm1(m)*stp(np1)
              do k=1, nele(m)
                s = s - enu(k,m)*stp(eidx(k,m))
              end do
              smax = max(smax, abs(s))
            end do
            if (smax > step_max) then
              b = step_max/smax
              do i=1, np1
                stp(i) = stp(i)*b
              end do
            endif

            ! Armijo backtracking on |F|^2
            lam = 1
            do
              do i=1, np1
                yt(i) = yv(i) - lam*stp(i)
              end do
              call eval(yt, fit, fmt, Ft, norm2t)
              if (norm2t <= (1d0 - 1d-4*lam)*norm2 .or. lam <= lam_min) exit
              lam = 0.5d0*lam
            end do

            ! accept the (possibly damped) step
            smax = 0
            do i=1, np1
              smax = max(smax, abs(yt(i) - yv(i)))
              yv(i) = yt(i)
              Fv(i) = Ft(i)
            end do
            do i=1, nions
              fi(i) = fit(i)
            end do
            do m=1, nmol
              fm(m) = fmt(m)
            end do
            norm2 = norm2t

            ! no movement in ln f: rounding level reached, stop iterating
            if (smax < 1d-13) exit
          end do

          if (.not. converged) then
            errf = abs(Fv(np1)) / max(fe_c, xne_c)
            do i=1, nions
              errf = max(errf, abs(Fv(i))/abund(i))
            end do
            converged = errf < tol_stall
          endif
          if (.not. converged) then
            nfail = nfail + 1
            worst = max(worst, errf)
          endif

          ! --- write back: fe, every ladder stage, molecular fractions ----
          fe(x,y,z) = fe_c
          do i=1, nions
            pt => ions(i)%p
            pt%f(x,y,z) = fi(i)
            k = 0
            do while (associated(pt%above))
              k = k + 1
              pt%above%f(x,y,z) = fi(i)*rat(k,i)
              pt => pt%above
            end do
            k = 0
            pt => ions(i)%p
            do while (associated(pt%below))
              k = k - 1
              pt%below%f(x,y,z) = fi(i)*rat(k,i)
              pt => pt%below
            end do
          end do
          do m=1, nmol
            molecules(m)%f(x,y,z) = fm(m)
          end do
        end do
      end do
    end do

    if (nfail > 0) then
      print *, 'molecular_equilibrium: ', nfail, '/', nx*ny*nz, ' points not converged after', &
        niter_max, ' iterations; worst scaled residual:', worst
    endif

    call trace%end(itimer)

  contains
    !---------------------------------------------------------------------------
    ! Residual F and its ingredients at state yy = (ln f_1..ln f_nions, ln fe).
    ! Also updates the host fe_c to exp(yy(np1)), so that after the accepted
    ! line-search evaluation fe_c always matches the current state.
    !---------------------------------------------------------------------------
    subroutine eval(yy, fatom, fmol, resid, nrm)
      real(dp), intent(in)  :: yy(:)
      real(dp), intent(out) :: fatom(:), fmol(:), resid(:), nrm
      real(dp) :: a
      integer  :: ii, mm, kk
      !-------------------------------------------------------------------------
      fe_c = exp(yy(np1))
      do ii=1, nions
        fatom(ii) = exp(yy(ii))
        resid(ii) = fatom(ii)*fil(ii) - abund(ii)
      end do
      resid(np1) = xne_c - fe_c
      do ii=1, nions
        resid(np1) = resid(np1) + fatom(ii)*ce(ii)
      end do
      do mm=1, nmol
        a = Cm(mm) - nm1(mm)*yy(np1)
        do kk=1, nele(mm)
          a = a + enu(kk,mm)*yy(eidx(kk,mm))
        end do
        ! No physical solution ever needs f_mol anywhere near 1 (even H2 stays
        ! below ~0.5); the old min(a,0) cap enforcing "f_mol <= 1" was a
        ! spurious plateau that Newton's search path legitimately crosses
        ! while still converging, and its zero slope there is inconsistent
        ! with the analytic Jacobian below, breaking descent and stalling
        ! the iteration. Clamp only far out, purely as an overflow guard, so
        ! fmol stays smooth (exp's own derivative matches the Jacobian)
        ! throughout the whole physically reachable range.
        fmol(mm) = exp(min(a, 300d0))
        do kk=1, nele(mm)
          resid(eidx(kk,mm)) = resid(eidx(kk,mm)) + enu(kk,mm)*fmol(mm)
        end do
        resid(np1) = resid(np1) + qmol(mm)*fmol(mm)
      end do
      nrm = 0
      do ii=1, np1
        nrm = nrm + resid(ii)**2
      end do
    end subroutine eval

    !---------------------------------------------------------------------------
    ! Find the solve-list index of the neutral root of constituent ion c, and
    ! the charge that c carries relative to it
    !---------------------------------------------------------------------------
    subroutine constituent_root(c, idx, charge)
      type(ion_type), pointer, intent(in) :: c
      integer, intent(out) :: idx, charge
      type(ion_type), pointer :: p2
      !-------------------------------------------------------------------------
      p2 => c
      charge = 0
      do ! search downwards, covering neutral and positive-ion constituents
        idx = solve_index(p2)
        if (idx > 0) return
        if (.not. associated(p2%below)) exit
        p2 => p2%below
        charge = charge + 1
      end do
      p2 => c
      charge = 0
      do while (associated(p2%above)) ! covers negative-ion constituents
        p2 => p2%above
        charge = charge - 1
        idx = solve_index(p2)
        if (idx > 0) return
      end do
      print *, 'molecular_equilibrium: no solve element for constituent ', c%id
      stop
    end subroutine constituent_root

    !---------------------------------------------------------------------------
    function solve_index(p2) result(idx)
      type(ion_type), pointer, intent(in) :: p2
      integer :: idx, ii
      !-------------------------------------------------------------------------
      do ii=1, nions
        if (associated(ions(ii)%p, p2)) then
          idx = ii
          return
        endif
      end do
      idx = 0
    end function solve_index
  END SUBROUTINE molecular_equilibrium

  !===============================================================================
  SUBROUTINE allocate_molec (self, m)
    class(molec_type) :: self
    integer :: m(3), mold(3), s, p
    !---------------------------------------------------------------------------
    if (all(m == self%m)) return
    self%m = m
    s = storage_size(self%f)
    if (self%allocated) then
      mold = shape(self%f)
      p = product(mold)
      deallocate(self%f, self%dis, self%part, self%n)
      call os%mem_register(s, -p, label='mol%f')
      call os%mem_register(s, -p, label='mol%n')
      call os%mem_register(s, -p, label='mol%dis')
      call os%mem_register(s, -p, label='mol%part')
    endif
    p = product(m)
    allocate(self%f(m(1),m(2),m(3)), source=0.d0)
    allocate(self%n, self%dis, self%part, source=self%f)
    call os%mem_register(s, p, label='mol%f')
    call os%mem_register(s, p, label='mol%n')
    call os%mem_register(s, p, label='mol%dis')
    call os%mem_register(s, p, label='mol%part')
    self%allocated = .true.
  END SUBROUTINE allocate_molec

  !===============================================================================
  SUBROUTINE deallocate_molec (self)
    class(molec_type) :: self
    integer :: m(3), s, p
    !---------------------------------------------------------------------------
    m = self%m
    p = product(m)
    s = storage_size(self%f)
    if (self%allocated) then
      deallocate(self%f, self%dis)
      call os%mem_register(s, -p, label='mol%f')
      call os%mem_register(s, -p, label='mol%dis')
      if (.not. (self%in_opac .or. self%in_linelist)) then
        deallocate(self%n, self%part)
        call os%mem_register(s, -p, label='mol%n')
        call os%mem_register(s, -p, label='mol%part')
      endif
    endif
    self%allocated = .false.
  END SUBROUTINE deallocate_molec

  !===============================================================================
  SUBROUTINE generate_isotopic_combinations(self)
    class(molec_type), intent(inout) :: self
    class(ion_type), pointer :: ion
    class(isotope_type), pointer :: iso, iso2
    integer :: i, j, k
    integer, allocatable :: current_iso(:), order(:)
    real(8) :: total_mass, total_frac
    character(16) :: mass_id
    logical :: done
    !---------------------------------------------------------------------------
    ! Estimate number of isotopic combinations (upper bound)
    self%niso = 1
    do i=1, self%nions
      ion => self%ions(i)%p
      if (associated(ion%isotope)) then
        self%niso = self%niso * (ion%niso+1)
      endif
    end do

    allocate(self%isotope(self%niso))

    ! Initialize index tracking for each ion
    allocate(current_iso(self%nions))
    allocate(order(self%nions))
    call argsort_int([(nint(self%ions(i)%p%awgt), i=1, self%nions)], order)
    current_iso = [(merge(0,1,associated(self%ions(order(i))%p%isotope)), i=1, self%nions)]
    k = 1
    done = .false.

    ! Generate all isotopic combinations
    do while (.not. done)
      ! Compute total mass and fraction
      total_mass = 0.0d0
      total_frac = 1.0d0
      mass_id = ''

      do i=1, self%nions
        ion => self%ions(order(i))%p
        if (.not. associated(ion%isotope)) then
          total_mass = total_mass + ion%awgt
          write(mass_id, '(A, I0)') trim(mass_id), nint(ion%awgt)
        else
          if (current_iso(i)==0) then
            total_mass = total_mass + ion%awgt
            write(mass_id, '(A, I0)') trim(mass_id), 0
          else
            iso => ion%isotope(current_iso(i))
            total_mass = total_mass + iso%awgt
            total_frac = total_frac * iso%frac
            write(mass_id, '(A, I0)') trim(mass_id), nint(iso%awgt)
          endif
        endif
      end do

      ! Store computed isotope in molecule
      iso => self%isotope(k)
      iso%awgt = total_mass
      iso%frac = total_frac
      read(mass_id, *) iso%iso
      k = k + 1

      ! Generate next combination
      do i=self%nions, 1, -1
        ion => self%ions(order(i))%p
        current_iso(i) = current_iso(i) + 1
        if (current_iso(i) <= ion%niso) exit
        current_iso(i) = merge(0,1,associated(ion%isotope))
        if (i == 1) done = .true.
      end do
    end do
    ! sum up duplicates
    ! this isn't quite perfect yet as it generates extra isotopes
    ! for H2, C2, etc. with frac > 1, as e.g. (0)12 and 120 are summed.
    ! It's fine as long as we don't use them though.
    ! Also haven't tested this for triatomic molecules yet which do not contain H.
    do i=1, self%niso
      iso => self%isotope(i)
      do k=i+1, self%niso
        iso2 => self%isotope(k)
        if (iso%awgt == iso2%awgt .and. iso%frac == iso2%frac) then
          iso%frac = iso%frac + iso2%frac
          iso2%frac = iso%frac
        endif
      end do
    end do
    deallocate(current_iso)
  END SUBROUTINE generate_isotopic_combinations

END MODULE molecule_mod
