MODULE molecule_mod
  USE hop_env_mod, only: os, trace
  USE console_mod
  USE prec_mod
  USE op_math_mod
  USE ion_mod
  USE handyfuncs_mod
  implicit none

  type, extends(particle_type) :: molec_type
    integer :: nions
    class(ion_pointer_t), allocatable :: ions(:)
    !> atomic number of each constituent, parallel to ions(:)
    integer, allocatable :: z(:)
    real(dp) :: dc(5) !< dissociation_constants
    !> Total nuclear spin statistical weight, prod over constituent nuclei of
    !> (2I+1). This is a pure convention marker for %part: g_nuclear > 1 means
    !> the tabulated partition function counts nuclear spin states, as modern
    !> molecular line lists do (the BT2/Barber H2O g column carries the 3:1
    !> ortho/para weight, so its Q must carry the matching factor 4 = 2*2*1).
    !> The line-opacity path therefore uses %part as-is, while get_dis divides
    !> it out to reach the nuclear-spin-free convention the atomic partition
    !> functions use. The factor cancels identically in the equilibrium
    !> constant anyway (it appears once in Q_mol and once in prod Q_atom), so
    !> this only keeps the two sides on the same footing; it is never physics.
    real(dp) :: g_nuclear = 1
    ! The molecule's ionisation energy in eV is its %xion, the limit the line-list
    ! Cowley Stark and Unsold van der Waals estimates are hydrogenic in (prepare_stark,
    ! prepare_vanderwaals). It is the optional 6th header field; a molecule without it
    ! keeps xion = 0, so its lines get neither estimate, and
    ! molecule_damping_estimates = .false. (bbtrans_mod) withholds both from all
    ! molecules. The dissociation energy is no such limit: it ends the nuclear-motion
    ! levels, not the Rydberg series of the electron.
    real(dp), dimension(:), allocatable :: de !< dissociation_energies
    real(dp), dimension(:,:,:), allocatable :: dis
  contains
    procedure :: read_molec_data
    procedure :: get_dis
    procedure :: calc_e
    procedure :: allocate => allocate_molec
    procedure :: deallocate => deallocate_molec
    procedure :: generate_isotopic_combinations
    procedure :: canonical_iso
  end type

  logical, public :: irwin_dis=.true.

contains

  !===============================================================================
  SUBROUTINE read_molec_data(self, lunit, iostat)
    class(molec_type), intent(inout) :: self
    integer, intent(in) :: lunit
    integer, intent(out) :: iostat
    integer :: defix, n_irwin, n_de, ios
    character(256) :: line
    !---------------------------------------------------------------------------
    ! Pull the header in as a whole record. The trailing nuclear-spin weight and
    ! polynomial low-temperature limit are optional, and parsing them out of an
    ! internal file lets a missing field fail cleanly -- a list-directed read
    ! straight off the unit would instead run on into the next record looking
    ! for the item.
    do
      read(lunit, '(a)', iostat=iostat) line
      if (iostat /= 0) return
      if (len_trim(line) > 0 .and. line(1:1) /= '!') exit
    end do
    if (allocated(self%de)) deallocate(self%de)
    if (allocated(self%irwin)) deallocate(self%irwin)
    self%g_nuclear = 1
    self%irwin_tmin = 1000
    self%xion = 0
    read(line, *, iostat=ios) self%id, defix, n_irwin, self%g_nuclear, self%irwin_tmin, &
                              self%xion
    if (ios /= 0) then
      self%xion = 0
      read(line, *, iostat=ios) self%id, defix, n_irwin, self%g_nuclear, self%irwin_tmin
    endif
    if (ios /= 0) then
      self%irwin_tmin = 1000
      read(line, *, iostat=ios) self%id, defix, n_irwin, self%g_nuclear
    endif
    if (ios /= 0) then
      ! 3-field header: partition function is in the nuclear-spin-free
      ! convention, so nothing to divide out
      self%g_nuclear = 1
      read(line, *, iostat=iostat) self%id, defix, n_irwin
      if (iostat /= 0) return
    endif
    ! ionisation stage, the Z of the damping estimates: 2 for a molecular ion (H2+)
    self%j = merge(2, 1, index(self%id, '+') > 0)
    n_de = merge(6, 1, defix==0)
    allocate(self%de(n_de))
    allocate(self%irwin(n_irwin))
    read(lunit, *) self%de
    read(lunit, *) self%dc ! always 5 numbers
    read(lunit, *) self%irwin
    self%irwin(:) = self%irwin(n_irwin:1:-1)
  END SUBROUTINE read_molec_data

  !===============================================================================
  !> Computes the equilibrium constant log(P(AB)/(P(A)*P(B))) for molecules
  !===============================================================================
  SUBROUTINE get_dis(self, temp, teta, e_kT, dis)
    class(molec_type), intent(in) :: self
    real(dp), intent(in), dimension(:,:,:) :: temp, teta, e_kT
    real(dp), intent(out), dimension(:,:,:) :: dis
    real(dp) :: awgt15, scalar, const_term
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

      ! %part may count nuclear spin states (see g_nuclear); the constituent
      ! ion%part above never do, so drop the factor here to put both sides of
      ! the equilibrium constant in the same convention
      scalar = const_term + 1.5*log(self%awgt) - log(self%g_nuclear)
      dis = dis + scalar
      dis = dis + self%part + self%de(1)*e_kT - (2.5*nm1)*log(temp)
    else
      call poly3d(teta, self%dc, dis)
      dis = -ln10*dis
    endif
  END SUBROUTINE get_dis

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
  !> Molecular + ionization equilibrium, solved cell by cell.
  !>
  !> The caller supplies the molecules and the solve elements: one neutral
  !> "root" ion per molecule-forming element (ions).  The routine runs at
  !> fixed electron density, so within one cell the Saha ratios between the
  !> ionization stages of every element, the electron pressure, and the
  !> dissociation constants are all constants.  With yi = ln f_i for the
  !> neutral fraction of each solve element and ye = ln fe, every molecular
  !> fraction is an explicit function of the unknowns,
  !>
  !>   ln f_mol = C_m + sum_i nu_im*yi - (n_m - 1)*ye ,
  !>
  !> where the per-cell constant C_m absorbs the dissociation constant, the
  !> electron pressure, and the Saha factors of charged constituents (H2+).
  !> The nonlinear system is number conservation for each solve element plus
  !> charge conservation,
  !>
  !>   F_i = f_i*phi_i + sum_m nu_im*f_mol - A_i           = 0
  !>   F_e = xnecno + sum_i q_i*f_i + sum_m q_m*f_mol - fe = 0 ,
  !>
  !> with phi_i (total nuclei per neutral) and q_i (electrons per neutral)
  !> given by the fixed Saha ladder of element i.  It is solved by a damped
  !> Newton method in log space: analytic Jacobian, step capped in ln f,
  !> Armijo backtracking on |F|^2.  One iteration costs nions+nmol+1
  !> exponentials plus an (nions+1)^2 linear solve; there are no logarithms,
  !> pointer walks, 3D-array accesses, or array-syntax statements inside the
  !> iteration, and the working set is O((nions+1)^2) scalars.
  !=========================================================================
  SUBROUTINE molecular_equilibrium(molecules, ions, tmolim, temp, teta, e_kT, &
                                   pe, xnecno, fe, ysave, same_state)
    class(molec_type), target, intent(inout) :: molecules(:)
    type(ion_pointer_t), intent(in) :: ions(:)
    real(dp), intent(in) :: tmolim
    real(dp), dimension(:,:,:), intent(in) :: temp, teta, e_kT, pe, xnecno
    real(dp), dimension(:,:,:), intent(inout) :: fe
    ! caller-owned per-cell seeds, carried across calls: the ne iteration
    ! re-solves the same cells at nearly the same state, so each cell's own
    ! previous solution is the best possible Newton start
    real(dp), dimension(:,:,:,:), allocatable, intent(inout), optional :: ysave
    ! temp (and hence teta, e_kT) is identical to the previous call's, so
    ! temp-only setup (mol%part, mol%dis) may be reused where valid
    logical, intent(in), optional :: same_state
    class(molec_type), pointer :: mol
    type(ion_type), pointer :: pt
    integer :: i, j, k, m, iq, nx, ny, nz, np1, nions, nmol, mxc
    integer :: x, y, z, iter, nfail, qlo, qhi
    ! fixed stoichiometry: per-constituent root index / charge offset, and
    ! per-molecule lists of distinct solve elements with their multiplicity
    integer, dimension(:),   allocatable :: ncons, nele, nm1, qmol
    integer, dimension(:,:), allocatable :: cidx, cq, eidx, enu
    ! transpose of the molecule->element map, for the scalar init solves
    integer, dimension(:),   allocatable :: nmele
    integer, dimension(:,:), allocatable :: meli, melnu
    ! per-cell constants
    real(dp), dimension(:,:), allocatable :: rat        ! stage/neutral ratios
    real(dp), dimension(:),   allocatable :: fil, ce    ! phi_i, q_i
    real(dp), dimension(:),   allocatable :: abund, Cm, avail, lockf, yprev
    integer,  dimension(:),   allocatable :: limel, limnu, plim
    real(dp) :: pe_c, lnpe, fe_c, xne_c, b, s, flock, sa, dsy, bb, dy
    ! scratch for the scalar solves: hoisted per-molecule constants and
    ! the last d H/d y, shared between solve_element and solve_pair
    real(dp), dimension(:), allocatable :: bcon, b1c
    integer,  dimension(:), allocatable :: n2p
    real(dp) :: dH_last
    ! Newton state
    real(dp), dimension(:),   allocatable :: yv, yt, Fv, Ft, stp
    real(dp), dimension(:),   allocatable :: fi, fit, fm, fmt
    real(dp), dimension(:,:), allocatable :: Amat
    ! last converged state, seeding the next cell's warm start, and the
    ! previous row's converged states (the row-below neighbor is closer
    ! than the row-wrapped previous cell at line starts)
    real(dp), dimension(:),   allocatable :: ywarm
    real(dp), dimension(:,:), allocatable :: ywrow
    logical,  dimension(:),   allocatable :: wrow_ok
    logical :: have_warm, use_save, have_seed, skip_setup, same_alloc
    real(dp) :: norm2, norm2t, errf, worst
    real(dp) :: y1, worsty
    integer  :: sweep
    logical  :: converged
    integer,  parameter :: niter_max = 100
    integer,  parameter :: nround    = 2          ! allocation auction rounds
    integer,  parameter :: nsweep    = 32         ! Gauss-Seidel sweeps/round
    real(dp), parameter :: tol       = 1d-10      ! scaled residual target
    real(dp), parameter :: tol_stall = 1d-6       ! acceptable when stalled
    real(dp), parameter :: step_max  = 3d0        ! cap on |d ln f| per step
    real(dp), parameter :: lam_min   = 1d0/1024d0 ! smallest damping factor
    integer :: itimer=0
    character(len=256) :: msg
    !---------------------------------------------------------------------------
    call trace%begin ('molecule_mod%molecular_equilibrium', itimer=itimer)

    nx = size(temp, 1)
    ny = size(temp, 2)
    nz = size(temp, 3)
    nmol  = size(molecules)
    nions = size(ions)
    np1   = nions + 1

    ! mol%part and mol%dis depend only on temp (teta and e_kT are functions
    ! of it), so when the caller certifies temp is unchanged they can be
    ! reused from the previous call -- provided every constituent ion's
    ! partition function is also temp-only (irwin or fixed); ne-dependent
    ! level dissolution in ion%part voids the cache, which such ions signal
    ! via state_dependent_part (H I under mhd_eos)
    skip_setup = .false.
    if (present(same_state)) then
      if (same_state) then
        skip_setup = .true.
        molcheck: do m=1, nmol
          mol => molecules(m)
          if (.not. all(mol%m == shape(temp))) then
            skip_setup = .false.
            exit molcheck
          endif
          do j=1, mol%nions
            pt => mol%ions(j)%p
            if (.not. (pt%fixed_part .or. allocated(pt%irwin)) &
                .or. pt%state_dependent_part) then
              skip_setup = .false.
              exit molcheck
            endif
          end do
        end do molcheck
      endif
    endif

    ! molecule partition functions and equilibrium constants; the scalar
    ! irwin_partf call keeps the elemental expansion out of the hot files
    do m=1, nmol
      mol => molecules(m)
      call mol%allocate(shape(temp))
      if (.not. skip_setup) then
        do z=1, nz
          do y=1, ny
            do x=1, nx
              mol%part(x,y,z) = mol%irwin_partf(temp(x,y,z))
            end do
          end do
        end do
        call mol%get_dis(temp, teta, e_kT, mol%dis)
      endif
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

    ! transpose of the molecule->element map: for every solve element the
    ! molecules that contain it, and with which multiplicity (used by the
    ! scalar Gauss-Seidel solves of the initial guess)
    allocate(nmele(nions), meli(nmol,nions), melnu(nmol,nions))
    do i=1, nions
      nmele(i) = 0
    end do
    do m=1, nmol
      do k=1, nele(m)
        i = eidx(k,m)
        nmele(i) = nmele(i) + 1
        meli(nmele(i), i)  = m
        melnu(nmele(i), i) = enu(k,m)
      end do
    end do

    allocate(rat(qlo:qhi, nions), fil(nions), ce(nions), abund(nions), Cm(nmol))
    allocate(avail(nions), lockf(nmol), limel(nmol), limnu(nmol), plim(nmol))
    allocate(yprev(np1))
    allocate(yv(np1), yt(np1), Fv(np1), Ft(np1), stp(np1))
    allocate(fi(nions), fit(nions), fm(nmol), fmt(nmol), Amat(np1,np1))
    allocate(ywarm(np1), bcon(nmol), b1c(nmol), n2p(nmol))
    allocate(ywrow(np1,nx), wrow_ok(nx))
    have_warm = .false.
    wrow_ok = .false.
    use_save = present(ysave)
    if (use_save) then
      if (allocated(ysave)) then
        if (size(ysave,1) /= np1 .or. size(ysave,2) /= nx .or. &
            size(ysave,3) /= ny  .or. size(ysave,4) /= nz) deallocate(ysave)
      endif
      if (.not. allocated(ysave)) then
        allocate(ysave(np1,nx,ny,nz))
        ysave = huge(1d0) ! sentinel: no seed stored yet
      endif
    endif

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

          ! ln K of each molecule in f-units, up to the -(n-1)*ye term.
          ! A charged constituent's ladder ratio can underflow to exactly 0
          ! in extreme cells; the floor keeps log() finite (the molecule
          ! then simply drops out) instead of injecting -Inf into Cm.
          do m=1, nmol
            b = molecules(m)%dis(x,y,z) + nm1(m)*lnpe
            do j=1, ncons(m)
              if (cq(j,m) /= 0) b = b + log(max(rat(cq(j,m), cidx(j,m)), 1d-300))
            end do
            Cm(m) = b
          end do

          ! --- warm start --------------------------------------------------
          ! In a smooth atmosphere the previous cell's converged state is
          ! usually already inside Newton's basin, so try the damped Newton
          ! directly from it, and only fall back to the full auction +
          ! relaxation seed where that fails (first cold cell, sharp
          ! fronts, branch flips).  A Newton that reaches tol lands on the
          ! same root regardless of its seed, so this changes nothing but
          ! the iteration count.
          ! one warm attempt from the best available seed: this cell's own
          ! previous solution, else the row-below neighbor, else the
          ! previously solved cell
          converged = .false.
          have_seed = .false.
          ! the ysave test must stay nested: Fortran does not short-circuit
          ! .and., so folding use_save in would touch an unallocated ysave
          if (use_save) then
            if (ysave(1,x,y,z) < 0.5d0*huge(1d0)) then
              yv(1:np1) = ysave(1:np1,x,y,z)
              have_seed = .true.
            endif
          endif
          if (.not. have_seed .and. wrow_ok(x)) then
            yv(1:np1) = ywrow(:,x)
            have_seed = .true.
          endif
          if (.not. have_seed .and. have_warm) then
            yv(1:np1) = ywarm
            have_seed = .true.
          endif
          if (have_seed) call newton_solve(40, converged)

          if (.not. converged) then
          ! eval() tracks fe_c with the state, so a failed warm start
          ! leaves it clobbered; re-seed it from the caller's Saha sum
          fe_c = max(fe(x,y,z), tiny(1d0))

          ! --- initial guess ----------------------------------------------
          ! Stoichiometric allocation of super-saturated molecules.  In
          ! cold cells a molecule's equilibrium so overwhelmingly favors
          ! association (dissociation energies of 100..500 kT) that
          ! essentially all of its scarcest constituent ends up bound.
          ! Which molecule binds which element is a branch choice that
          ! neither relaxation nor the damped Newton below can repair once
          ! it is made wrongly: with the atomic fractions underflowed,
          ! moving budget from one molecule to another passes through
          ! states of larger |F| (a genuine non-monotone barrier), so an
          ! Armijo-damped iteration can never cross it -- it stalls with
          ! some residual pinned at a finite fraction of an abundance.
          !
          ! The allocation is found as a small auction, run in nround
          ! rounds.  Within a round, molecules repeatedly bid for their
          ! constituents: a molecule's bid is the excess of its ln f_mol --
          ! evaluated at the current element depths -- over what its
          ! scarcest remaining budget could supply.  The highest bidder
          ! locks (1 - 1e-3 of, so budgets stay positive) the scarcest
          ! remaining budget; its abundant constituents are re-seeded from
          ! their leftover budgets, while the limiting element is seeded by
          ! inverting ln f_lock = c_m + sum_k nu_k y_k, which is essential:
          ! seeded from its 1e-3 leftover instead, the first sweep below
          ! would solve the molecule's abundant partner against a molecule
          ! term that is orders of magnitude too large and re-pin the wrong
          ! branch (the classic f_CO = A_C / f_CO = A_O seesaw).
          !
          ! One round is not enough, because early bids are evaluated
          ! against not-yet-deepened partners: SiO outbids H2O for O while
          ! H is still at full atomic strength, steals Si's budget, and
          ! starves SiS even though S has nowhere else to go.  The second
          ! round therefore re-ranks all bids using each element's depth
          ! from the end of the previous round (excluding, for each
          ! molecule, the element it locked itself, whose depth it caused).
          ! Crucially the carried depths are the ones refined by the
          ! Gauss-Seidel sweeps below, not the raw auction seeds: an
          ! auction-internal inversion can be badly stale (computed before
          ! a partner deepened) or reflect a starved crumb of a leftover
          ! budget, and one such fake price handed to the next round makes
          ! the wrong molecule win an element outright.  The sweeps produce
          ! self-consistent depths for whatever allocation the round chose,
          ! which is exactly what a price should be.
          do m=1, nmol
            plim(m) = 0
          end do
          do i=1, np1
            yprev(i) = huge(1d0)
          end do
          yv(np1) = log(fe_c)

          do iter=1, nround
            ! fresh budgets and atomic seeds each round (the -600 floor
            ! covers phi_i overflowing to Inf in many-stage ladders of
            ! hot, thin cells, where log() would seed y = -Inf and strand
            ! every refinement below at -Inf as well)
            do i=1, nions
              avail(i) = abund(i)
              yv(i) = max(log(abund(i)/fil(i)), -600d0)
            end do
            do m=1, nmol
              lockf(m) = -1d0
            end do

            do ! auction: highest bidder locks its scarcest budget
              iq = 0
              b = 0
              do m=1, nmol
                if (lockf(m) >= 0d0) cycle
                s = Cm(m) - nm1(m)*yv(np1)
                flock = huge(1d0)
                do k=1, nele(m)
                  i = eidx(k,m)
                  y1 = yv(i)
                  if (i /= plim(m)) y1 = min(y1, yprev(i))
                  s = s + enu(k,m)*y1
                  flock = min(flock, log(avail(i)) - log(dble(enu(k,m))))
                end do
                s = s - flock
                if (s > b) then
                  iq = m
                  b = s
                endif
              end do
              if (iq == 0) exit

              flock = huge(1d0)
              do k=1, nele(iq)
                if (avail(eidx(k,iq))/enu(k,iq) < flock) then
                  flock = avail(eidx(k,iq))/enu(k,iq)
                  limel(iq) = eidx(k,iq)
                  limnu(iq) = enu(k,iq)
                endif
              end do
              flock = (1d0 - 1d-3)*flock
              lockf(iq) = flock
              if (flock <= 0d0) cycle
              do k=1, nele(iq)
                i = eidx(k,iq)
                avail(i) = avail(i) - enu(k,iq)*flock
                yv(i) = min(yv(i), max(log(avail(i)/fil(i)), -600d0))
              end do
              b = Cm(iq) - nm1(iq)*yv(np1)
              do k=1, nele(iq)
                if (eidx(k,iq) /= limel(iq)) b = b + enu(k,iq)*yv(eidx(k,iq))
              end do
              i = limel(iq)
              s = (log(lockf(iq)) - b) / limnu(iq)
              yv(i) = min(yv(i), max(s, -600d0))
            end do

            ! a re-ranked auction that reproduces the previous round's
            ! allocation would only re-derive the depths the sweeps below
            ! already refined for it: restore those and stop instead
            if (iter > 1) then
              same_alloc = .true.
              do m=1, nmol
                if (alloc_code(m) /= plim(m)) same_alloc = .false.
              end do
              if (same_alloc) then
                yv(1:np1) = yprev
                exit
              endif
            endif

          ! Refine with Gauss-Seidel sweeps in which every element solves
          ! its own scalar conservation equation exactly, holding the
          ! other elements and ye fixed:
          !
          !   H(y) = phi_i e^y + sum_m nu_im f_m(y)  =  A_i .
          !
          ! Solved as ln H(y) = ln A_i: ln of a positive sum of increasing
          ! exponentials is increasing and convex with slope in [1, nu_max],
          ! so undamped Newton converges globally and fast (the step cap
          ! only guards exp overflow).  An earlier heuristic here only
          ! redistributed the excess of oversaturated molecules among their
          ! constituents; with several molecules sharing an element it
          ! could push that element's y so low that e^y and all its
          ! molecular terms underflowed to zero, handing the damped Newton
          ! below an exactly singular Jacobian row and a residual pinned at
          ! -A_i (the persistent "scaled residual 1.0" stalls).  The exact
          ! scalar solves cannot underflow that way, because the molecular
          ! term itself must balance A_i.  Charge conservation re-solves ye
          ! once per sweep (scalar Newton, floored like the seed fe).
          ! Two-element budget exchanges (Si moving between SiO and SiS)
          ! relax with a contraction factor that can approach 1, hence the
          ! generous sweep count; each sweep is a handful of scalar ops.
            do sweep=1, nsweep
              worsty = 0
              do i=1, nions
                y1 = yv(i)
                call solve_element(i)
                worsty = max(worsty, abs(yv(i) - y1))
              end do

              ! Two elements coupled through a molecule are then re-solved
              ! jointly: single-element relaxation exchanges budget between
              ! two molecules of a shared element (Si moving between SiO
              ! and SiS) at under an e-fold per sweep across gaps of
              ! hundreds, and it is exactly this one-at-a-time structure
              ! that creates the two locally stable branches of the CO
              ! seesaw.  The pair subsystem itself has an SPD Jacobian --
              ! a strictly monotone 2D field -- hence a unique root and no
              ! branch ambiguity at all; the nested Schur solve below finds
              ! it from any starting point.
              do m=1, nmol
                if (nele(m) /= 2) cycle
                b = yv(eidx(1,m))
                s = yv(eidx(2,m))
                call solve_pair(eidx(1,m), eidx(2,m))
                worsty = max(worsty, abs(yv(eidx(1,m))-b), abs(yv(eidx(2,m))-s))
              end do

            ! ye from charge conservation at the current element guesses.
            ! The atomic part is constant in ye and the charged-molecule
            ! part has a known slope, so the scalar equation ye = ln s(ye)
            ! is solved exactly with Newton (one step when no molecule is
            ! charged, since s is then constant).  The former 0.5-damped
            ! fixed point here only shed half the ye error per sweep, so
            ! the 1d-10 exit below was gated at ~30 sweeps per round even
            ! in cells whose element solves had long converged.
            sa = xne_c
            do i=1, nions
              sa = sa + exp(yv(i))*ce(i)
            end do
            ! the ye-independent part of each charged molecule's ln f
            ! argument, Cm(m) + sum_j enu(j,m)*yv(eidx(j,m)), is constant
            ! across the ye Newton (only yv(np1) moves below); hoist it out
            ! of the k-loop, reusing bcon whose contents are dead here
            do m=1, nmol
              if (qmol(m) == 0) cycle
              bb = Cm(m)
              do j=1, nele(m)
                bb = bb + enu(j,m)*yv(eidx(j,m))
              end do
              bcon(m) = bb
            end do
            b = yv(np1)
            do k=1, 20
              s   = sa
              dsy = 0
              do m=1, nmol
                if (qmol(m) == 0) cycle
                bb  = qmol(m)*exp(min(bcon(m) - nm1(m)*yv(np1), 600d0))
                s   = s   + bb
                dsy = dsy - nm1(m)*bb
              end do
              s = max(s, 1d-3*xne_c, tiny(1d0))
              ! Newton on g(ye) = ln s(ye) - ye: g' = s'/s - 1 <= -1 for
              ! positively charged molecules; the min() keeps the slope
              ! magnitude >= 1 should a negatively charged one ever appear
              dy = (log(s) - yv(np1)) / (1d0 - min(dsy/s, 0d0))
              dy = max(min(dy, 60d0), -60d0)
              yv(np1) = yv(np1) + dy
              if (abs(dy) < 1d-12) exit
            end do
            worsty = max(worsty, abs(yv(np1) - b))
            if (worsty < 1d-7) exit
            end do
            fe_c = exp(yv(np1))

            ! remember this round's outcome: sweep-refined depths become
            ! the prices of the next round's auction, and each molecule is
            ! exempted from the price of the element it locked itself
            do m=1, nmol
              plim(m) = alloc_code(m)
            end do
            yprev = yv(1:np1)
          end do

          ! --- damped Newton to convergence -------------------------------
          call newton_solve(niter_max, converged)
          endif ! warm start failed

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
            ! don't re-seed the next call/row from a stalled state
            if (use_save) ysave(1,x,y,z) = huge(1d0)
            wrow_ok(x) = .false.
          else
            ! seed for the next cell's warm start and this cell's next call
            ywarm      = yv
            ywrow(:,x) = yv
            have_warm = .true.
            wrow_ok(x) = .true.
            if (use_save) ysave(1:np1,x,y,z) = yv
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

    if (nfail > 0 .and. console_level(ch_hop) >= 1) then
      write(msg,'(a,i0,a,i0,a,i0,a,es10.3)') 'molecular_equilibrium: ', nfail, '/', nx*ny*nz, &
        ' points not converged after ', niter_max, ' iterations; worst scaled residual ', worst
      call console_warn (ch_hop, msg)
    endif

    call trace%end(itimer)

  contains
    !---------------------------------------------------------------------------
    !> The molecule's current auction allocation encoded as one integer: the
    !> element it locked (limel) when it holds a lock, else 0.  Single-sourced
    !> so the round-to-round comparison and the store into plim cannot drift.
    !---------------------------------------------------------------------------
    integer function alloc_code(mm)
      integer, intent(in) :: mm
      if (lockf(mm) >= 0d0) then
        alloc_code = limel(mm)
      else
        alloc_code = 0
      endif
    end function alloc_code

    !---------------------------------------------------------------------------
    !> Damped Newton on the full coupled system from the current yv, for at
    !> most nit iterations: analytic Jacobian, step capped in ln f, Armijo
    !> backtracking on |F|^2.  Leaves yv, fi, fm, Fv, norm2, fe_c consistent
    !> with the last accepted state and mutates the shared Newton scratch
    !> (Amat, stp, yt, fit, fmt, errf, norm2t); ok reports scaled-residual
    !> convergence.
    !---------------------------------------------------------------------------
    subroutine newton_solve(nit, ok)
      integer, intent(in)  :: nit
      logical, intent(out) :: ok
      real(dp) :: ss, smx, lm
      integer  :: ii, jj, kk, mm, itn
      logical  :: bad
      !-------------------------------------------------------------------------
      call eval(yv, fi, fm, Fv, norm2)
      ok = .false.

      do itn=1, nit
        errf = abs(Fv(np1)) / max(fe_c, xne_c)
        do ii=1, nions
          errf = max(errf, abs(Fv(ii))/abund(ii))
        end do
        if (errf < tol) then
          ok = .true.
          exit
        endif

        ! analytic Jacobian dF/d(ln f) at the current state
        do jj=1, np1
          do ii=1, np1
            Amat(ii,jj) = 0
          end do
        end do
        do ii=1, nions
          Amat(ii,ii)   = fi(ii)*fil(ii)
          Amat(np1,ii) = fi(ii)*ce(ii)
        end do
        Amat(np1,np1) = -fe_c
        do mm=1, nmol
          do kk=1, nele(mm)
            ii = eidx(kk,mm)
            Amat(ii,np1) = Amat(ii,np1) - (enu(kk,mm)*nm1(mm))*fm(mm)
            do jj=1, nele(mm)
              Amat(ii,eidx(jj,mm)) = Amat(ii,eidx(jj,mm)) + (enu(kk,mm)*enu(jj,mm))*fm(mm)
            end do
          end do
          if (qmol(mm) /= 0) then
            Amat(np1,np1) = Amat(np1,np1) - (qmol(mm)*nm1(mm))*fm(mm)
            do kk=1, nele(mm)
              Amat(np1,eidx(kk,mm)) = Amat(np1,eidx(kk,mm)) + (qmol(mm)*enu(kk,mm))*fm(mm)
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
        do ii=1, np1
          stp(ii) = Fv(ii)
        end do
        call lineq(Amat, stp, np1)
        ! lineq has no singularity check: a zero pivot yields NaN. Such
        ! a step must never reach the state (a NaN trial would pass the
        ! backtracking loop at lam_min and poison yv, and through the
        ! write-back fe(x,y,z), permanently); keep the last finite
        ! state and let the stall test in the caller judge it instead.
        bad = .false.
        do ii=1, np1
          if (stp(ii) /= stp(ii)) bad = .true.
        end do
        if (bad) exit
        smx = 0
        do ii=1, np1
          smx = max(smx, abs(stp(ii)))
        end do
        do mm=1, nmol
          ss = nm1(mm)*stp(np1)
          do kk=1, nele(mm)
            ss = ss - enu(kk,mm)*stp(eidx(kk,mm))
          end do
          smx = max(smx, abs(ss))
        end do
        if (smx > step_max) then
          ss = step_max/smx
          do ii=1, np1
            stp(ii) = stp(ii)*ss
          end do
        endif

        ! Armijo backtracking on |F|^2
        lm = 1
        do
          do ii=1, np1
            yt(ii) = yv(ii) - lm*stp(ii)
          end do
          call eval(yt, fit, fmt, Ft, norm2t)
          if (norm2t <= (1d0 - 1d-4*lm)*norm2 .or. lm <= lam_min) exit
          lm = 0.5d0*lm
        end do

        ! a non-finite trial slipped through at lam_min: restore fe_c
        ! (eval clobbered it) and keep the last finite state
        if (norm2t /= norm2t) then
          fe_c = exp(yv(np1))
          exit
        endif

        ! accept the (possibly damped) step
        smx = 0
        do ii=1, np1
          smx = max(smx, abs(yt(ii) - yv(ii)))
          yv(ii) = yt(ii)
          Fv(ii) = Ft(ii)
        end do
        do ii=1, nions
          fi(ii) = fit(ii)
        end do
        do mm=1, nmol
          fm(mm) = fmt(mm)
        end do
        norm2 = norm2t

        ! no movement in ln f: rounding level reached, stop iterating
        if (smx < 1d-13) exit
      end do
    end subroutine newton_solve

    !---------------------------------------------------------------------------
    !> Solve element ii's scalar conservation equation over yv(ii), all other
    !> unknowns held fixed:  H(y) = phi_ii e^y + sum_m nu_m f_m(y) = A_ii .
    !> Newton on ln H - ln A: ln of a positive sum of increasing exponentials
    !> is increasing and convex with slope in [1, nu_max], so this converges
    !> globally; the step cap only guards exp overflow.
    !---------------------------------------------------------------------------
    subroutine solve_element(ii)
      integer, intent(in) :: ii
      real(dp) :: yy1, HH, dHH, bb, ss, ddy
      integer  :: kk, mm, jj, itn
      !-------------------------------------------------------------------------
      ! each containing molecule's ln f_mol argument minus this element's
      ! own contribution is constant across the Newton iteration
      do kk=1, nmele(ii)
        mm = meli(kk,ii)
        bb = Cm(mm) - nm1(mm)*yv(np1) - melnu(kk,ii)*yv(ii)
        do jj=1, nele(mm)
          bb = bb + enu(jj,mm)*yv(eidx(jj,mm))
        end do
        bcon(kk) = bb
      end do
      yy1 = yv(ii)
      do itn=1, 50
        HH  = fil(ii)*exp(yy1)
        dHH = HH
        do kk=1, nmele(ii)
          ss  = exp(min(bcon(kk) + melnu(kk,ii)*yy1, 600d0))
          HH  = HH  + melnu(kk,ii)*ss
          dHH = dHH + melnu(kk,ii)**2*ss
        end do
        if (HH < 1d-290) then
          ddy = -60d0 ! everything underflowed: march up blind
        else
          ddy = (log(HH) - log(abund(ii))) * (HH/dHH)
          ddy = max(min(ddy, 60d0), -60d0)
        endif
        yy1 = yy1 - ddy
        if (abs(ddy) < 1d-9) exit
      end do
      yv(ii) = yy1
      dH_last = dHH ! d H/d y at the solution, reused as J22 by solve_pair
    end subroutine solve_element

    !---------------------------------------------------------------------------
    !> Jointly solve the two conservation equations of elements i1 and i2
    !> (everything else fixed): the inner loop keeps F_i2 = 0 exactly via
    !> solve_element, the outer Newton moves yv(i1) along the constraint with
    !> the Schur-complement slope J11 - J12^2/J22 of the pair's SPD Jacobian.
    !---------------------------------------------------------------------------
    subroutine solve_pair(i1, i2)
      integer, intent(in) :: i1, i2
      real(dp) :: H1, J11, J12, J22, bb, ss, ddy, slope
      integer  :: kk, mm, jj, itn
      !-------------------------------------------------------------------------
      ! each of i1's molecules' ln f_mol argument minus the contributions
      ! of the two solved elements is constant across the outer iteration;
      ! n2p records i2's multiplicity in it (0 when i2 is not a member)
      do kk=1, nmele(i1)
        mm = meli(kk,i1)
        bb = Cm(mm) - nm1(mm)*yv(np1)
        n2p(kk) = 0
        do jj=1, nele(mm)
          if (eidx(jj,mm) == i1) cycle
          if (eidx(jj,mm) == i2) then
            n2p(kk) = enu(jj,mm)
          else
            bb = bb + enu(jj,mm)*yv(eidx(jj,mm))
          endif
        end do
        b1c(kk) = bb
      end do

      do itn=1, 40
        call solve_element(i2)
        ! solve_element left d H2/d y2 at the solution in dH_last
        J22 = dH_last

        H1  = fil(i1)*exp(yv(i1))
        J11 = H1
        J12 = 0
        do kk=1, nmele(i1)
          ss  = exp(min(b1c(kk) + melnu(kk,i1)*yv(i1) + n2p(kk)*yv(i2), 600d0))
          H1  = H1  + melnu(kk,i1)*ss
          J11 = J11 + melnu(kk,i1)**2*ss
          J12 = J12 + melnu(kk,i1)*n2p(kk)*ss
        end do

        if (H1 < 1d-290) then
          ddy = -60d0
        else
          ! the Schur complement can cancel to roundoff when one shared
          ! molecule dominates both rows; the floor then just yields a
          ! large step in the right direction, harmlessly capped
          slope = J11 - J12*J12/max(J22, 1d-300)
          slope = max(slope, 1d-12*J11)
          ddy = (log(H1) - log(abund(i1))) * (H1/slope)
          ddy = max(min(ddy, 60d0), -60d0)
        endif
        yv(i1) = yv(i1) - ddy
        if (abs(ddy) < 1d-9) exit
      end do
      call solve_element(i2)
    end subroutine solve_pair

    !---------------------------------------------------------------------------
    !> Residual F and its ingredients at state yy = (ln f_1..ln f_nions, ln fe).
    !> Also updates the host fe_c to exp(yy(np1)), so that after the accepted
    !> line-search evaluation fe_c always matches the current state.
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
    !> Find the solve-list index of the neutral root of constituent ion c, and
    !> the charge that c carries relative to it
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
      call console_error ('molecular_equilibrium: no solve element for constituent '//trim(c%id))
      error stop
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
  !> Isotopologue labels follow Turbospectrum's mass part: three digits per constituent,
  !> constituents in ascending atomic number, mass numbers ascending among equal atomic
  !> numbers, 000 for the natural mixture.  12C14N is 012014, 1H12C is 000012 since
  !> hydrogen has no isotope table, and 12C13C is 012013.  The label is stored as an
  !> integer, which drops the leading zeros; the number of constituents restores them.
  !> A constituent without an isotope table is always 000.  Combinations that differ
  !> only in the order of equal atoms share a label, and their fractions are summed.
  !===============================================================================
  SUBROUTINE generate_isotopic_combinations(self)
    class(molec_type), intent(inout) :: self
    class(ion_type), pointer :: ion
    class(isotope_type), pointer :: iso
    integer :: i, k, n
    integer, allocatable :: current_iso(:), order(:), mass(:), label(:)
    real(8), allocatable :: awgt(:), frac(:)
    logical, allocatable :: first(:)
    logical :: done
    !---------------------------------------------------------------------------
    n = 1
    do i=1, self%nions
      ion => self%ions(i)%p
      if (associated(ion%isotope)) n = n * (ion%niso+1)
    end do
    allocate(label(n), awgt(n), frac(n), first(n))

    allocate(current_iso(self%nions), order(self%nions), mass(self%nions))
    call argsort_int(self%z, order)
    current_iso = 0
    k = 0
    done = .false.
    do while (.not. done)
      k = k + 1
      awgt(k) = 0
      frac(k) = 1
      do i=1, self%nions
        ion => self%ions(order(i))%p
        if (current_iso(i) == 0) then
          awgt(k) = awgt(k) + ion%awgt
          mass(i) = 0
        else
          iso => ion%isotope(current_iso(i))
          awgt(k) = awgt(k) + iso%awgt
          frac(k) = frac(k) * iso%frac
          mass(i) = nint(iso%awgt)
        endif
      end do
      label(k) = pack_label(self, mass)

      ! next combination; constituents without an isotope table stay natural
      do i=self%nions, 1, -1
        ion => self%ions(order(i))%p
        if (associated(ion%isotope)) then
          current_iso(i) = current_iso(i) + 1
          if (current_iso(i) <= ion%niso) exit
          current_iso(i) = 0
        endif
        if (i == 1) done = .true.
      end do
      if (self%nions == 0) done = .true.
    end do

    first = [(.not. any(label(:i-1) == label(i)), i=1, n)]
    self%niso = count(first)
    allocate(self%isotope(self%niso))
    k = 0
    do i=1, n
      if (.not. first(i)) cycle
      k = k + 1
      iso => self%isotope(k)
      iso%iso = label(i)
      iso%awgt = awgt(i)
      iso%frac = sum(frac, mask=label == label(i))
    end do
  END SUBROUTINE generate_isotopic_combinations

  !===============================================================================
  !> The isotopologue label of a line-list header in the convention above.  An
  !> explicit mass for a constituent without an isotope table that equals its
  !> natural mass number (hydrogen written 001) reads as the natural mixture, and the
  !> masses of equal atoms are put in ascending order.
  !===============================================================================
  INTEGER FUNCTION canonical_iso(self, iso) result(label)
    class(molec_type), intent(in) :: self
    integer, intent(in) :: iso
    integer, allocatable :: order(:), mass(:)
    integer :: i, rest
    !---------------------------------------------------------------------------
    allocate(order(self%nions), mass(self%nions))
    call argsort_int(self%z, order)
    rest = iso
    do i=self%nions, 1, -1
      mass(i) = mod(rest, 1000)
      rest = rest / 1000
    end do
    if (rest /= 0) then
      label = -1
      return
    end if
    do i=1, self%nions
      associate(ion => self%ions(order(i))%p)
        if (.not. associated(ion%isotope) .and. mass(i) == nint(ion%awgt)) mass(i) = 0
      end associate
    end do
    label = pack_label(self, mass)
  END FUNCTION canonical_iso

  !===============================================================================
  !> Masses given in ascending atomic number, packed three digits each.
  !===============================================================================
  INTEGER FUNCTION pack_label(self, mass) result(label)
    class(molec_type), intent(in) :: self
    integer, intent(inout) :: mass(:)
    integer, allocatable :: order(:), zs(:)
    integer(8) :: packed
    integer :: i, j, t
    !---------------------------------------------------------------------------
    allocate(order(self%nions))
    call argsort_int(self%z, order)
    allocate(zs(self%nions))
    zs = self%z(order)
    ! insertion sort of the masses within each run of equal atomic number
    do i=2, size(mass)
      j = i
      do while (j > 1)
        if (zs(j-1) /= zs(j) .or. mass(j-1) <= mass(j)) exit
        t = mass(j); mass(j) = mass(j-1); mass(j-1) = t
        j = j - 1
      end do
    end do
    packed = 0
    do i=1, size(mass)
      packed = packed * 1000 + mass(i)
      if (packed > huge(label)) then
        call console_error('isotopologue label of '//trim(self%id)//' exceeds the integer range')
        error stop
      end if
    end do
    label = int(packed)
  END FUNCTION pack_label


END MODULE molecule_mod
