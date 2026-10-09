MODULE composition_mod
  USE hop_env_mod, only: io, os, trace
  USE console_mod
  USE prec_mod
  USE op_math_mod
  USE handyfuncs_mod
  USE element_mod
  USE ion_mod, only: ion_type, particle_type, debye_shield
  USE level_mod, only: level_type, wcalc
  USE bbtrans_mod, only:tolerance
  USE bftrans_mod, only:bftrans_type, H_occupation_probability
  USE molecule_mod, only:molec_type, irwin_dis, molecular_equilibrium
  USE hprofl_mod, only:hline_t, hydrogen_line, add_h_opacities
  USE hvcs_mod, only:hydrogen_vcs_table, ensure_vcs_loaded
  USE anstee_mod, only:ensure_anstee_loaded
  USE line_types_mod, only: linelist_type
  USE linelist_reader_mod, only: ll_file_t, ll_block_t, ll_open, ll_close, &
                                 ll_next_block, ll_block_range, ll_whole_block, &
                                 ll_read_slice, ll_read_ranges, ll_index_records, ll_has_records
  USE linelist_mod, only: build_linelist, finalize_linelist, add_particle_opacity, &
                          append_h_line, bij_const, merge_ll_diag, &
                          voigt_peak, in_any_interval, print_ll_particle, &
                          ll_profile, ll_nvisit, ll_npass, ll_nvoigt, ll_nmarch
  USE spectrum_mod, only: spectrum_type, window_file
  USE optab_mod, only: pca_transform_type
  USE osmet_mod, only:osmet, compute_epsilon_stage
  USE omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none
  private

  ! One species' answer for one abundance set: either the populations a
  ! molecular re-solve produced, or -- for an element outside the network, whose
  ! stages simply follow its abundance at fixed (T, ne) -- the scalar that
  ! multiplies the base populations. Both are resolved when the sets are
  ! prepared, so the opacity loop only ever matches an id against this table.
  type :: set_population_t
    integer :: i_set=0
    character(len=6) :: id=''
    real(dp) :: scale=1d0
    real(dp), dimension(:,:,:), pointer :: n => null()
  end type set_population_t

  type, public :: composition_type
    logical :: initialised, ll_initialised=.false., loaded=.false.
    integer :: nel, nmolec, nline
    real(dp) :: xmy, suma, grph, tmolim, molh, FeH
    type(element_type), allocatable :: elements(:)
    type(molec_type),   allocatable :: molecules(:)
    type(bftrans_type), allocatable :: trans(:)
    type(level_type), dimension(:), allocatable :: input_H_levels
    type(particle_type) :: peachmet
    real(dp), dimension(:,:,:), allocatable :: temp, ne, rho, zero
    real(dp), dimension(:,:,:), allocatable :: e_kT, pe, dxi, toth, tothi
    real(dp), dimension(:,:,:), allocatable :: saha, lsaha, stark_pe
    ! per-cell warm-start seeds for molecular_equilibrium, carried across
    ! calls (the ne iteration re-solves the same cells with nearly the
    ! same state); allocated and maintained by molecular_equilibrium
    real(dp), dimension(:,:,:,:), allocatable :: moleq_seed
    ! electron fraction per hydrogen nucleus the last molecular equilibrium
    ! converged on. It is not ne/toth: toth is rebuilt from the model's own rho,
    ! while this is the fe that balanced the molecules, and re-entering the
    ! solver with anything else moves a converged state (see
    ! prepare_abundance_sets)
    real(dp), dimension(:,:,:), allocatable :: moleq_fe
    ! set by set_aux_data: temp identical to the previous call's, so
    ! temp-only derived quantities (molecular part/dis) may be reused
    logical :: same_temp = .false.
    ! converged occupation-weighted H I log-partition (mhd_eos), carried
    ! across thermal_equilibrium calls as a warm start for the inner loop
    real(dp), dimension(:,:,:), allocatable :: part_hocc
    real(dp), dimension(:,:,:), allocatable :: teta, pg, E, p6, shield, fo, vbar2h, t3nhe, eh
    real(4), dimension(:,:,:), allocatable :: chi500, vmic
    ! velocities in units of c, already masked by use_vel; the rotated copies each
    ! carry their own, which is what makes them the right home for the RT Doppler shift
    real(4), dimension(:,:,:), allocatable :: vx, vy, vz
    ! partial derivatives of ln(T) and of Pe (linear, not log) per code length
    ! unit, evaluated on the PCHIP interpolant of the input snapshot
    ! (set_gradients), plus a real(4) copy of Pe itself; consumed by the
    ! chain-rule 'hermite' RT solvers, which assume chi/Pe linear along the ray.
    ! Rotated copies carry only the z pair (their columns already run along the
    ! ray); only the unrotated composition holds the horizontal partials
    real(4), dimension(:,:,:), allocatable :: dlnt_dx, dlnt_dy, dlnt_dz
    real(4), dimension(:,:,:), allocatable :: dpe_dx, dpe_dy, dpe_dz
    real(4), dimension(:,:,:), allocatable :: pe4
    ! PCA transform of (lnT, lnNe) used by the sampled opacity table grid;
    ! ferried from init_tabulated_comp down to each rt_patch's opac table
    type(pca_transform_type) :: pca
    real(dp), allocatable :: set_da(:,:)       ! (atomic number, output set)
    logical, allocatable :: set_element(:,:)   ! explicit manifest entries
    type(set_population_t), allocatable :: set_pop(:)
  contains
    procedure :: init => injon
    procedure :: load_state_variables
    procedure :: set_gradients
    procedure :: iterate_thermal_eq
    procedure :: thermal_equilibrium
    procedure :: compute_energy_and_gas_pressure
    procedure :: set_aux_data
    procedure :: init_H_levels
    procedure :: compute_shield
    procedure :: Hoccupation_partition
    procedure :: molecule_by_id
    procedure :: element_by_id
    procedure :: particle_by_id
    procedure :: ion_by_id
    procedure :: boltzmann_factor
    procedure :: continuum_opacities
    procedure :: add_absmet_opacities
    procedure :: opacity
    procedure :: add_hrest
    procedure :: electron_scattering
    procedure :: rayleigh_scattering
    procedure :: reassign_pointer
    procedure :: get_totn
    procedure :: compute_grams_per_hydrogen
    procedure :: compute_damping_aux
    procedure :: overwrite_abundance
    procedure :: init_linelist
    procedure :: mask_linelist
    procedure :: mask_departure_species
    procedure :: add_linelist_opacities
    procedure :: merge_line_diagnostics
    procedure :: read_molec_file
    procedure :: parse_mol_id
    procedure :: print_opacity_contributions
    procedure :: get_table_sp
    procedure :: prepare_abundance_sets
    procedure :: set_population
    procedure :: population_ratio
    procedure :: write_linelists
  end type

  ! The levels and bound-bound transitions of the model atom a departure set belongs
  ! to, read from its table-format atom file (departure_atom).  One instance per file
  ! is shared by every patch and every ray-aligned copy.
  type, public :: departure_atom_t
    character(len=512) :: file=''
    real(8), allocatable :: ev(:)        ! 100*stage + energy above the stage's ground [eV]
    real(8), allocatable :: level_ev(:)  ! energy above the atom's ground [eV]
    integer, allocatable :: line_lo(:), line_up(:)  ! level indices of each transition
    real(8), allocatable :: line_lam(:)  ! vacuum wavelength of each transition [AA]
    type(departure_atom_t), pointer :: next => null()
  end type departure_atom_t

  type, public :: departure_pointer_t
    character(len=6) :: key=''  ! species key: element symbol ('Mg') or molecule formula ('H2O')
    logical :: is_set=.false.   ! an unset key must never match anything
    real(4), pointer :: mem(:,:,:,:) => null()
    type(departure_atom_t), pointer :: atom => null()
    integer :: nlevel=0
  contains
    procedure :: set_key => departure_set_key
    procedure :: matches => departure_matches
    procedure :: assign_levels => departure_assign_levels
    procedure :: assign_h_levels => departure_assign_h_levels
  end type departure_pointer_t

  type(departure_atom_t), pointer, save :: departure_atoms => null()
  public :: departure_atom

  character(128), public, save :: exclude_elements='', only_elements=''
  ! feed the Hummer-Mihalas occupation-weighted H I partition function into the
  ! ionization and molecular equilibrium (pressure ionization consistent with
  ! the H_occupation_probability opacity treatment)
  logical, public :: mhd_eos=.true.
  logical, public :: helium_scattering=.true.
  real(dp), public :: debug_opac=0
  real(dp), public :: ref_lam=5000.d0
  character(len=1024), public :: abundances
  logical, public :: relative_abundances=.false.
  logical, public :: read_departures=.false.
  ! match line-list lines without level indices to the model atom as Turbospectrum's
  ! line-list conversion does (departure_assign_levels); without it only lines with
  ! both indices are in NLTE
  logical, public :: match_nlte_levels=.false.
  logical, public :: sample_linelist=.false.
  ! split LTE line-list opacity into true absorption and scattering with the
  ! stage-correct van Regemorter epsilon (osmet_mod%compute_epsilon_stage);
  ! NLTE lines are untouched -- their scattering lives in the rate equations
  logical, public :: line_scattering=.false.
  real(4), public :: opac_table_vmic=-1

  character(len=3), parameter, public :: &
    element_ids(92) = [ &
    'H ','He','Li','Be','B ','C ','N ','O ','F ', &
    'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar', &
    'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co', &
    'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
    'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh', &
    'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
    'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu', &
    'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
    'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl', &
    'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
    'Pa','U ']

  character(len=5), parameter :: roman(0:10) = [ &
    ' -   ', ' I   ', ' II  ', ' III ', ' IV  ', ' V   ', &
    ' VI  ', ' VII ', ' VIII', ' IX  ', ' X   ' ]

contains

  !===============================================================================
  ! Set the species key of a departure set.  An empty id leaves the set unmatched.
  !-------------------------------------------------------------------------------
  SUBROUTINE departure_set_key(self, id)
    class(departure_pointer_t) :: self
    character(len=*), intent(in) :: id
    !---------------------------------------------------------------------------
    self%key = adjustl(id)
    self%is_set = len_trim(self%key) > 0
  END SUBROUTINE departure_set_key

  !===============================================================================
  ! True when this departure set belongs to the species of ptc.  Ion ids are
  ! '<Sym> <ROMAN>' (H- is 'H -', so it too carries the blank) while molecule ids
  ! never contain a blank, so the token before the first blank is the element
  ! symbol and a blank-free id is a whole molecule.  The match is at species level:
  ! key 'Mg' covers Mg I, Mg II, ... -- dep%atom%ev spans every stage, tagged 100*stage --
  ! but never MgH or Mn I.
  !-------------------------------------------------------------------------------
  LOGICAL FUNCTION departure_matches(self, ptc) RESULT(hit)
    class(departure_pointer_t), intent(in) :: self
    class(particle_type), intent(in) :: ptc
    character(len=6) :: pid
    integer :: ib
    !---------------------------------------------------------------------------
    hit = .false.
    if (.not. self%is_set) return
    pid = adjustl(ptc%id)
    ib = index(trim(pid), ' ')
    if (ib > 0) pid = pid(1:ib-1)
    hit = uppercase(trim(pid)) == uppercase(trim(self%key))
  END FUNCTION departure_matches

  !===============================================================================
  ! The levels and transitions of the table-format model atom in file, read on the
  ! first request and shared afterwards.  The format is the one atom_mod%read_pandas
  ! reads; level and transition indices are 0-based in the file.
  !-------------------------------------------------------------------------------
  FUNCTION departure_atom(file) RESULT(atom)
    character(len=*), intent(in) :: file
    type(departure_atom_t), pointer :: atom
    integer :: lunit, nlevel, nline, idx, l
    integer, allocatable :: stage(:)
    real(8) :: g, lam_air, ground_ev
    !---------------------------------------------------------------------------
    !$omp critical (departure_atom_cache)
    atom => departure_atoms
    do while (associated(atom))
      if (trim(atom%file) == trim(file)) exit
      atom => atom%next
    end do
    if (.not. associated(atom)) then
      call os%check_file(file)
      allocate(atom)
      atom%file = file
      open(newunit=lunit, file=file, form='formatted', status='old', action='read')
      call NextEntry(lunit)
      read(lunit, *)                     ! id
      call NextEntry(lunit)
      read(lunit, *)                     ! abnd awgt titles
      call NextEntry(lunit)
      read(lunit, *)                     ! abundance, atomic weight
      call NextEntry(lunit)
      read(lunit, *)                     ! nlevel nline nbf ncol titles
      call NextEntry(lunit)
      read(lunit, *) nlevel, nline
      call NextEntry(lunit)
      read(lunit, *)                     ! level header
      allocate(atom%ev(nlevel), atom%level_ev(nlevel), stage(nlevel))
      do l=1, nlevel
        call NextEntry(lunit)
        read(lunit, *) idx, atom%level_ev(l), g, stage(l)
      end do
      ground_ev = atom%level_ev(1)
      do l=1, nlevel
        if (l > 1) ground_ev = merge(atom%level_ev(l), ground_ev, stage(l) > stage(l-1))
        atom%ev(l) = 100*stage(l) + atom%level_ev(l) - ground_ev
      end do
      allocate(atom%line_lo(nline), atom%line_up(nline), atom%line_lam(nline))
      if (nline > 0) then
        call NextEntry(lunit)
        read(lunit, *)                   ! transition header
      end if
      do l=1, nline
        call NextEntry(lunit)
        read(lunit, *) idx, atom%line_lo(l), atom%line_up(l), lam_air
        atom%line_lo(l) = atom%line_lo(l) + 1
        atom%line_up(l) = atom%line_up(l) + 1
        atom%line_lam(l) = air2vac(lam_air)
      end do
      close(lunit)
      atom%next => departure_atoms
      departure_atoms => atom
    end if
    !$omp end critical (departure_atom_cache)
  END FUNCTION departure_atom

  !===============================================================================
  ! The departure-set levels of the H I lines.  As in Turbospectrum (hydropac.f), the
  ! principal quantum numbers n_lo/n_up, which lev_i/lev_j start as, are the indices of
  ! the atom's levels, so the atom must have one level per n, in order from n=1.  Each
  ! index is checked against the line's energy (lower_ev, and lower_ev + h nu0 for the
  ! upper end) to TS's relative tolerance; a mismatch stops the run.  A line with an n
  ! beyond the atom's H I levels is set to -1 and stays in LTE.
  !-------------------------------------------------------------------------------
  SUBROUTINE departure_assign_h_levels(self, hll)
    class(departure_pointer_t), intent(in) :: self
    type(linelist_type), intent(inout) :: hll
    real(8), parameter :: ts_ev_rel=1d-3
    integer :: l
    !---------------------------------------------------------------------------
    do l = 1, hll%nline
      if (hll%lev_i(l) < 0) cycle
      if (.not. (bound(hll%lev_i(l)) .and. bound(hll%lev_j(l)))) then
        hll%lev_i(l) = -1
        hll%lev_j(l) = -1
        cycle
      end if
      call check_energy(hll%lev_i(l), real(hll%lower_ev(l), 8))
      call check_energy(hll%lev_j(l), hll%lower_ev(l) + H_E*hll%nu0(l))
    end do
  CONTAINS
    ! level n is an H I level of the atom
    LOGICAL FUNCTION bound(n)
      integer, intent(in) :: n
      bound = n <= size(self%atom%ev)
      if (bound) bound = int(self%atom%ev(n)/100) == 1
    END FUNCTION bound

    SUBROUTINE check_energy(n, e)
      integer, intent(in) :: n
      real(8), intent(in) :: e
      real(8) :: e_atom
      character(len=160) :: msg
      e_atom = self%atom%ev(n) - 100
      if (abs(e_atom - e) > ts_ev_rel*max(e_atom, e)) then
        write(msg, '(a,i0,a,f9.4,a,i0,a,f9.4,a)') 'H I line level n=', n, ' at ', e, &
          ' eV does not match atom level ', n, ' at ', e_atom, ' eV'
        call console_error(trim(msg))
        error stop 'assign_h_levels: the H atom''s levels are not one per n'
      end if
    END SUBROUTINE check_energy
  END SUBROUTINE departure_assign_h_levels

  !===============================================================================
  ! Give each line of ptc's line lists the departure-set levels of its two ends.  A
  ! level index from the line list is kept and must name a level of the line's own
  ! stage.  Without match_nlte_levels a line is in NLTE only when the list gives both
  ! indices.  With it, the ends without an index (0) are matched as Turbospectrum's
  ! convert_lte_to_nlte.py assigns indices: a line without any index first takes both
  ! levels of the first atom transition of its stage within ts_dlam of its wavelength
  ! whose lower level agrees with lower_ev to the factor ts_ev_factor.  An end still
  ! unmatched takes the nearest level of the stage within that factor, the last of
  ! equally near ones as in TS.  Both tolerances are TS's.  The energies are compared
  ! above the atom's ground, the line's lower_ev shifted by the atom's energy of the
  ! stage's first level, so the window is wider for ions.  An end that finds no level
  ! is set to -1, which leaves the line in LTE, as does an end left at 0.
  !-------------------------------------------------------------------------------
  SUBROUTINE departure_assign_levels(self, ptc)
    class(departure_pointer_t), intent(in) :: self
    class(particle_type), intent(inout) :: ptc
    real(8), parameter :: ts_dlam=0.02d0, ts_ev_factor=1.0005d0
    real(8) :: ground_ev, e_lo, e_up
    integer :: kiso, l, t, lo
    !---------------------------------------------------------------------------
    associate(atom => self%atom)
    ground_ev = 0
    lo = findloc(int(atom%ev/100), ptc%j, 1)
    if (lo > 0) ground_ev = atom%level_ev(lo)
    do kiso = 0, ubound(ptc%ll, 1)
      associate(sll => ptc%ll(kiso))
      do l = 1, sll%nline
        if (sll%lev_i(l) > 0) call check_level(sll%lev_i(l))
        if (sll%lev_j(l) > 0) call check_level(sll%lev_j(l))
        if (.not. match_nlte_levels) cycle
        e_lo = ground_ev + sll%lower_ev(l)
        e_up = e_lo + H_E*sll%nu0(l)
        if (sll%lev_i(l) == 0 .and. sll%lev_j(l) == 0) then
          do t = 1, size(atom%line_lo)
            lo = atom%line_lo(t)
            if (int(atom%ev(lo)/100) /= ptc%j) cycle
            if (abs(atom%line_lam(t) - exp(sll%llam0(l))) < ts_dlam .and. &
                within(atom%level_ev(lo), e_lo)) then
              sll%lev_i(l) = lo
              sll%lev_j(l) = atom%line_up(t)
              exit
            end if
          end do
        end if
        if (sll%lev_i(l) == 0) sll%lev_i(l) = nearest_within(e_lo)
        if (sll%lev_j(l) == 0) sll%lev_j(l) = nearest_within(e_up)
      end do
      end associate
    end do
    end associate
  CONTAINS
    SUBROUTINE check_level(lev)
      integer, intent(in) :: lev
      if (lev > size(self%atom%ev)) &
        error stop 'assign_levels: a line-list level index exceeds the atom''s levels'
      if (int(self%atom%ev(lev)/100) /= ptc%j) &
        error stop 'assign_levels: a line-list level index names a level of another stage'
    END SUBROUTINE check_level

    ! TS's acceptance window: e lies within the factor ts_ev_factor of level_ev
    LOGICAL FUNCTION within(level_ev, e)
      real(8), intent(in) :: level_ev, e
      within = level_ev*ts_ev_factor >= e .and. e >= level_ev/ts_ev_factor
    END FUNCTION within

    INTEGER FUNCTION nearest_within(e) RESULT(lev)
      real(8), intent(in) :: e
      integer :: k
      lev = -1
      do k = 1, size(self%atom%ev)
        if (int(self%atom%ev(k)/100) /= ptc%j) cycle
        if (.not. within(self%atom%level_ev(k), e)) cycle
        if (lev < 0) then
          lev = k
        else if (abs(self%atom%level_ev(k) - e) <= abs(self%atom%level_ev(lev) - e)) then
          lev = k
        end if
      end do
    END FUNCTION nearest_within
  END SUBROUTINE departure_assign_levels

  !===============================================================================
  SUBROUTINE injon(self, absdat_file, abund_file, FeH, &
    isotope_file, molec_file, tmolim, external_abund, alpha_enhancement, &
    r_process_enhancement)
    class(composition_type), target :: self
    character(len=*), intent(in) :: absdat_file, abund_file
    character(len=*), intent(in), optional :: isotope_file, molec_file
    real(4), intent(in), optional :: tmolim
    real(8), intent(in) :: FeH
    real(4), intent(in), optional :: external_abund(:), alpha_enhancement, r_process_enhancement
    integer :: absdat, lunit, iostat
    integer :: i, j, k, l, m, nl, nk, niso
    integer :: lvl_count, skip_back
    real(dp) :: abund, xion, g2, dump, tecor, decor
    character(len=3) :: id
    character(len=6) :: key
    character(len=65) :: text
    character(len=128) :: text_line
    integer :: g0, nqfix, nqtemp, ifish, ncor
    real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c, scale
    real(dp), dimension(:,:), allocatable :: xkap
    integer :: nkomp,nextl,nutzl,nextt,nutzt,null
    integer :: ilogl, kvadl, minex, maxex, nnu, iscal
    integer :: ilogt, kvadt, minet, maxet, ntt, iteta, istim

    class(element_type),    pointer :: elements(:)
    class(element_type),    pointer :: element, H
    class(molec_type),      pointer :: molecule
    class(ion_type),        pointer :: ion, Hm, HI
    class(isotope_type),    pointer :: iso
    class(level_type),      pointer :: lvl
    class(bftrans_type), pointer :: trans
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%init', itimer=itimer)

106 FORMAT(4(F10.3,F10.4))
    call os%check_file(absdat_file)
    open(file=absdat_file, newunit=absdat)

    read(absdat, *) self%nel
    allocate(self%elements(self%nel))

    elements => self%elements

    read(absdat, *) elements%id
    read(absdat, *) elements%awgt
    read(absdat, *) elements%nion

    do i=1, self%nel
      element => elements(i)
      element%np = findloc(element_ids, element%id(:3), dim=1)

      if (i>1) then
        allocate(element%ions(element%nion))
      else
        allocate(element%ions(0:element%nion))
      endif
    end do

    call os%check_file(abund_file)
    open(file=abund_file, newunit=lunit, iostat=iostat)
    call nextentry(lunit)
    do while (iostat==0)
      read(lunit, *, iostat=iostat) id, abund
      call self%element_by_id(id, element, verbose=.true.)
      if(associated(element)) element%abund = abund
    end do
    close(lunit)

    if (elements(1)%abund == 1) then
      elements%abund = log10(elements%abund) + 12.0
    endif

    self%FeH = merge(FeH, 0.d0, FeH > -99)
    do i=1, self%nel
      element => elements(i)
      select case (element%id)
       case ('H', 'He', 'Li')
        ! do not scale with FeH
       case default
        element%abund = element%abund + self%FeH
      end select
      if (present(alpha_enhancement)) then
        select case (element%id)
         case ('O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti')
          element%abund = element%abund + alpha_enhancement
        end select
      endif
      if (present(r_process_enhancement)) then
        select case (element%id)
         case ('Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', &
           'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  &
           'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', &
           'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', &
           'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Th', 'U')
          element%abund = element%abund + r_process_enhancement
        end select
      endif
    end do

    if (present(external_abund)) then
      ! do not overwrite external abundances with FeH
      do i=1, size(external_abund)
        call self%element_by_id(element_ids(i), element, verbose=.true.)
        if(associated(element)) element%abund = external_abund(i)
      end do
    endif

    if (len_trim(abundances) > 0) then
      do i=1, self%nel
        element => elements(i)
        j = index(abundances, trim(element%id)//'=')
        if (j > 0) then
          l = len(abundances)
          l = j + index(abundances(j:l), ' ') - 1
          j = j + len_trim(element%id) + 1
          read(abundances(j:l), *) abund
          if (relative_abundances) then
            element%abund = element%abund + abund
          else
            element%abund = abund
          endif
          if (console_on(ch_hop, 1)) &
            write(console_unit,'(A, A, A, F8.3)') 'Changing abundance of ', trim(element%id), ' to ', element%abund
        endif
      end do
    endif

    elements%abund = 10**(elements%abund - elements(1)%abund)
    call self%compute_grams_per_hydrogen()

    do i=1, self%nel
      element => elements(i)

      do j=1, elements(i)%nion
        ion => element%ions(j)
        ion%j = j
        ion%dxi_fac = j

        !g0 = THE STATISTICAL WEIGHT OF THE GROUND LEVEL
        !nk = THE NUMBER OF ELECTRON CONFIGURATIONS CONSIDERED.

        read(absdat, *, iostat=iostat) g0
        backspace(absdat)

        if (g0 > 0) then
          read(absdat, *, iostat=iostat) g0, nk
          ion%nk = nk+1
          skip_back = 0
          do k=1, nk
            read(absdat, *) dump, dump, dump, nl
            skip_back = skip_back + 1 + ceiling(nl / 4.0)

            if (nl > 0) then
              ion%nk = ion%nk + nl
              read(absdat, 106) (xion, g2, m=1,nl)
            endif

          end do

          do k=1, skip_back
            backspace absdat
          end do

          allocate(ion%levels(ion%nk))
          ion%levels%xl=0
          ion%levels(1)%g = g0
          ion%levels(1)%xion = 0.d0

          lvl_count = 1
          do k=1, nk
            lvl_count = lvl_count + 1
            lvl => ion%levels(lvl_count)

            !xion = THE IONIZATION ENERGY IN ELECTRON VOLTS
            !g2   = THE STATISTICAL WEIGHT (2L+1)*(2J+1)
            !xl   = THE LOWEST QUANTUM NUMBER OF THE ASYMPTOTIC
            !       (HYDROGENIC) PARTOF THE PARTITION FUNCTION,
            !nl   = THE NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSION
            !       FOR THE 'MIDDLE PART' OF THE PARTITION FUNCTION ('QPRIME').

            read(absdat, *) lvl%xion, lvl%g, lvl%xl, nl

            if (nl > 0) then
              read(absdat, 106) (ion%levels(m)%xion, ion%levels(m)%g, m=lvl_count+1, lvl_count+nl)
              lvl_count = lvl_count + nl
            endif
          end do
          ion%xion = ion%levels(2)%xion

        else ! IRWIN FORMAT
          allocate(ion%irwin(6))
          read(absdat, *) g0, ion%xion, ion%irwin
          ion%irwin(:) = ion%irwin(6:1:-1)
        endif
      end do
    end do

    do i=1, self%nel
      element => elements(i)
      read(absdat, *) element%active, (element%ions(j)%active, j=1, element%nion)
    end do

    read(absdat, *) nqfix ! number of ions with constant partition FUNCTION

    do m=1, nqfix
      read(absdat, *) i, j, elements(i)%ions(j)%parco

      !IQFIX(I,J) = 0 MEANS THAT THE PARTITION FUNCTION IS CONSIDERED TO BE CONSTANT.
      elements(i)%ions(j)%fixed_part = .true.
    end do

    read(absdat, *) nqtemp ! number of ions with PF(temp) instead of PF(temp, pe)

    do m=1, nqtemp
      call console_error ('tabulated partition functions are not implemented')
      error stop
    end do

    ! IFISH IS A PARAMETER FOR THE CHOICE OF THE ASYMPTOTIC PARTITION
    ! FUNCTION. IFISH=0 MEANS THAT THE ASYMPTOTIC PART WILL BE EVALU-
    ! ATED FOLLOWING BASCHEK ET AL., ABH. HAMB. VIII,26 (1966). IFISH
    ! =1 MEANS THAT IT WILL BE EVALUATED FOLLOWING FISCHEL AND SPARKS
    ! ASTROPHYS. J. 164, 356 (1971).

    read(absdat, *) ifish

    ! WHEN MOLH IS GREATER THAN ZERO THE MOLECULAR FORMATION WILL BE COMPUTED
    ! IN SUBR. MOLEQ (ONLY H2 AND H2+), ELSE MORE COMPLETE MOLECULAR
    ! FORMATION WILL BE EVALUATED IN SUBR. MOL.
    read(absdat, *) self%tmolim, self%molh

    if (present(tmolim)) then
      if (tmolim > -98) then
        self%tmolim = tmolim
      endif
    endif

    read(absdat, '(A9)') text
    if (text(1:9)=='* CORONAL') then
      ! READS CORONAL ION BALANCE PARAMETERS
      read(absdat, *) ncor, tecor, decor

      do m=1, ncor
        read(absdat, *) i, j, acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
        element => elements(findloc(elements%np, i, dim=1))

        if (j <= element%nion) then
          ion => element%ions(j)

          ion%acolc = acolc
          ion%tcolc = tcolc
          ion%aradc = aradc
          ion%adic  = adic
          ion%bdic  = bdic
          ion%t0c   = t0c
          ion%t1c   = t1c
        endif
      end do
    else
      backspace absdat
    endif

    do i=1, self%nel
      element => elements(i)
      do j=lbound(element%ions, 1), ubound(element%ions, 1)
        ion => element%ions(j)

        if (j >= 0 .and. j <= 10) then
          ion%id = trim(element%id(:3)) // trim(roman(j))
        else
          error stop 'Ion id construction needs to be updated for j > 10'
        end if

        if(j < ubound(element%ions, 1)) ion%above => element%ions(j+1)
        if(j > 1) ion%below => element%ions(j-1)
        ion%awgt = element%awgt
        ion%abund => element%abund
      end do
    end do

    call self%element_by_id('H', H)
    Hm => H%ions(0) ! H- is ion number 0
    HI => H%ions(1)
    HI%below => Hm
    HI%nk = 15

    Hm%fixed_part = .true.
    Hm%parco = 1
    Hm%xion = 0.754
    Hm%dxi_fac = 2
    ! --- Different Hm ionisation in Oslo legacy code ---!
    !Hm%xion = 0.747
    ! --------------------------------------------------!
    RYD = HI%xion ! overwriting RYD constant
    ! --------------------------------------------------!
    self%peachmet%id = 'Me I' ! pseudo-ion for combination of minor metals

    call self%read_molec_file(molec_file)

    !101 FORMAT(7X,I3,5(9X,I1))
    !999 FORMAT(4(9X,I1),8X,I2,9X,I1)

    !997 FORMAT(6(F10.0)) ! lam format
    !998 FORMAT(6(F10.5)) ! xkap format
    !102 FORMAT(6F10.4)
    !105 FORMAT(F6.3,F4.0,F5.1,I5)

    !read(absdat,101) nkomp,nextl,nutzl,nextt,nutzt,null
    read(absdat, '(A)') text_line
    call march_line(text_line, 10, nkomp)
    call march_line(text_line, 10, nextl)
    call march_line(text_line, 10, nutzl)
    call march_line(text_line, 10, nextt)
    call march_line(text_line, 10, nutzt)
    call march_line(text_line, 10, null)

    allocate(self%trans(nkomp))

    self%elements(1)%ions%in_opac = .true.
    self%elements(2)%ions%in_opac = .true.
    call self%molecule_by_id('H2', molecule)
    if (associated(molecule)) molecule%in_opac = .true.

    do i=1, nkomp
      trans => self%trans(i)

      read(absdat,'(A65)') text
      text = uppercase(text)
      trans%label = text(:20)
      trans%key = text(1:5)
      trans%peach = index(trans%label, 'PEACH') > 0
      trans%lam3 = trans%peach

      select case(trim(trans%key))
       case('H- BF')
        trans%ion => Hm
       case('H BF')
        trans%ion => HI
        ! hydrogen energy level
        trans%k = extract_int(text(1:10))
        trans%lam3 = .true.
       case('H FF')
        trans%ion => HI
        trans%lam3 = .true.
       case('H- FF')
        trans%ion => HI
       case('H2+')
        trans%ion => HI
        trans%ion2 => HI%above
       case('H+H')
        trans%ion => HI
        trans%ion2 => HI
       case('H+HE')
        trans%ion => HI
        call self%particle_by_id('He I', trans%ion2)
       case('H2+H')
        trans%ion => HI
        call self%particle_by_id('H2', trans%ion2)
       case('H2+H2')
        call self%particle_by_id('H2', trans%ion)
        call self%particle_by_id('H2', trans%ion2)
       case('H2+HE')
        call self%particle_by_id('H2', trans%ion)
        call self%particle_by_id('He I', trans%ion2)
       case('H2-')
        call self%particle_by_id('H2', trans%ion)
       case('HE-')
        call self%particle_by_id('He I', trans%ion)
       case ('H2O-')
        call self%particle_by_id('H2O', trans%ion)
        trans%molecular = .true.
       case ('MEIFF')
        trans%ion => self%peachmet
        trans%regular = .true.
        trans%peach = .true.
       case default
        if (text(1:1) == '/') then
          trans%key = text(2:6)
          id = text(2:3)
          call self%element_by_id(id, element)
          read(absdat, *) j, xion, trans%level%g
          trans%level%xion = xion * hc_e
          trans%ion => element%ions(j)
          trans%regular = .true.
        else
          id = text(1:2)
          j = scan(id,'_')
          if (j > 1) id = id(1:j-1)
          trans%key = id // text(3:4)
          j = scan(trans%key, 'IV', back=.true.)
          if (j > 0) trans%key = trans%key(1:j)
          j = scan(trans%key, '-')
          trans%regular = j == 0
          if (j > 0) then
            key = trans%key(1:j-1)
          else
            key = trans%key
          endif
          call self%molecule_by_id(key, molecule)
          trans%molecular = associated(molecule)
          trans%ion => molecule
          if (trans%regular) then
            if (.not. trans%molecular) then
              call self%ion_by_id(trans%key, ion)
              if (associated(ion)) then
                trans%ion => ion
                if (trans%peach) trans%ion => ion%above
              endif
            endif
            trans%level%xion = 0
            trans%level%g = 1
          elseif (.not. trans%molecular) then
            call self%ion_by_id(key//' I', ion)
            trans%ion => ion
          endif
        endif
      end select

      if (.not. associated(trans%ion)) then
        call trans%deactivate()
        if (console_on(ch_hop, 2)) write(console_unit,'(2a)') trim(trans%key), ' deactivated'
      else
        trans%ion%in_opac = .true.
      endif
      if (associated(trans%ion2)) trans%ion2%in_opac = .true.

      !read(absdat,999) ilogl,kvadl,minex,maxex,nnu
      read(absdat, '(A)') text_line
      call march_line(text_line, 10, ilogl)
      call march_line(text_line, 10, kvadl)
      call march_line(text_line, 10, minex)
      call march_line(text_line, 10, maxex)
      call march_line(text_line, 12, nnu)

      trans%nnu = nnu
      trans%ilogl = ilogl > 0
      trans%kvadl = kvadl > 0
      trans%minex = minex > 0
      trans%maxex = maxex > 0

      if (trans%key == 'H BF') then
        trans%maxex  = H_occupation_probability
        !trans%loglam = H_occupation_probability
      endif

      allocate(trans%lam(nnu))
      allocate(trans%ldisk(nnu))

      trans%ldisk = .false.

      read(absdat, *, iostat=iostat) trans%lam

      if (iostat /= 0) then
        call console_error ('cannot read wavelength values for transition '//trim(trans%key))
        error stop
      endif

      if (trans%loglam) then
        allocate(trans%lnlam, source=log(trans%lam))
      endif

      trans%minl = trans%lam(1)
      trans%maxl = trans%lam(nnu)


      do j=1, nnu-1
        if (trans%lam(j+1) - trans%lam(j) <= 10.) then
          trans%ldisk(j) = .true.
        endif
      end do

      !read(absdat,999) ilogt,kvadt,minet,maxet,ntt,iteta
      read(absdat, '(A)') text_line
      call march_line(text_line, 10, ilogt)
      call march_line(text_line, 10, kvadt)
      call march_line(text_line, 10, minet)
      call march_line(text_line, 10, maxet)
      call march_line(text_line, 10, ntt)
      call march_line(text_line, 10, iteta)
      call march_line(text_line, 10, iscal)
      call march_line(text_line, 10, istim)

      if (iscal>0) read(absdat, *) scale

      trans%ntt = ntt
      trans%ilogt = ilogt > 0
      trans%kvadt = kvadt > 0
      trans%minet = minet > 0
      trans%maxet = maxet > 0
      trans%iteta = iteta > 0
      trans%pedep = iscal > 0
      trans%istim = istim > 0

      if (ntt > 1) then
        allocate(trans%temp(ntt))
        allocate(trans%tt  (ntt))
        read(absdat, *, iostat=iostat) trans%temp

        if (iostat /= 0) then
          call console_error ('cannot read temperature values for transition '//trim(trans%key))
          error stop
        endif

        if (iteta > 0) then
          allocate(trans%teta(ntt))
          trans%teta = trans%temp
          trans%temp = temp_to_theta / trans%teta(ntt:1:-1)
          trans%tt = trans%teta
        else
          trans%tt = trans%temp
        endif

        trans%mint = trans%temp(1)
        trans%maxt = trans%temp(ntt)
      endif

      allocate(xkap(nnu, ntt))

      do j=1, ntt
        read(absdat, '(6E10.3)') (xkap(k,j),k=1,nnu)
      end do

      if (iscal > 0) then
        xkap = xkap * scale
      endif

      call trans%set_table(xkap)
      deallocate(xkap)
    end do

    close(absdat)

    if (present(isotope_file)) then
      if (isotope_file /= '') then
        call os%check_file(isotope_file)
        open(file=isotope_file, newunit=lunit)
        read(lunit, *, iostat=iostat) text
        do while (iostat==0)
          do while (text(1:1) == '!')
            read(lunit, *, iostat=iostat) text
          end do
          backspace lunit
          read(lunit, *, iostat=iostat) id, niso
          call self%element_by_id(id, element)
          if (associated(element)) then
            element%niso = niso
            allocate(element%isotope(element%niso))
            do i=1, element%niso
              iso => element%isotope(i)
              read(lunit, *) iso%iso, iso%awgt, iso%frac
            end do
            do j=1, element%nion
              element%ions(j)%isotope => element%isotope
              element%ions(j)%niso = element%niso
            end do
          endif
          read(lunit, *, iostat=iostat) text
        end do
        close(lunit)
        do i=1, self%nmolec
          molecule => self%molecules(i)
          call molecule%generate_isotopic_combinations()
        end do
      endif
    endif

    if (allocated(HI%levels)) then
      allocate(self%input_H_levels, source=HI%levels)
    endif

    self%initialised = .true.
    call trace%end(itimer)
  contains

    !===========================================
    subroutine march_line(line, n, val)
      character(len=*) :: line
      integer :: i, n, val
      !---------------------------------------
      i = scan(line(1:n), '0123456789')
      read(line(i:n), *, iostat=iostat) val
      if (iostat /= 0) then
        call console_error ('cannot read fixed format line for transition '//trim(trans%key)//': '//trim(line))
        error stop
      endif
      line = line(n+1:)
    end subroutine march_line

  END SUBROUTINE injon

  !===============================================================================
  SUBROUTINE mask_linelist (self)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    integer :: i
    character(len=512) :: drop, only
    logical :: in_drop, no_only
    ! ----------------------------------------------------------------------------
    drop = ''
    only = ''
    drop(2:len(drop)-1) = trim_between_words(exclude_elements)
    only(2:len(only)-1) = trim_between_words(only_elements)
    do i=1, self%nel
      element => self%elements(i)
      in_drop = index(drop, ' '//trim(element%id)//' ') > 0
      no_only = index(only, ' '//trim(element%id)//' ') == 0 .and. only_elements /= ''
      if (in_drop .or. no_only) then
        element%ions%in_linelist = .false.
        if (element_has_lines(element)) then
          if (console_on(ch_hop, 1)) write(console_unit,'(3a)') 'removing ', trim(element%id), ' from linelist'
        endif
      endif
    end do
    do i=1, self%nmolec
      mol => self%molecules(i)
      in_drop = index(drop, ' '//trim(mol%id)//' ') > 0
      no_only = index(only, ' '//trim(mol%id)//' ') == 0 .and. only_elements /= ''
      if (in_drop .or. no_only) then
        mol%in_linelist = .false.
        if (has_lines(mol)) then
          if (console_on(ch_hop, 1)) write(console_unit,'(3a)') 'removing ', trim(mol%id), ' from linelist'
        endif
      endif
    end do
  END SUBROUTINE mask_linelist

  !===============================================================================
  ! Whether the line-list block of species id (an element symbol or a molecule
  ! formula, the first word of a block header) is to be read at all. This is the
  ! test mask_linelist applies to the composition, hoisted to read time so the
  ! blocks of unwanted species are never parsed: with only_elements naming the one
  ! species a run synthesises, the master list is skipped a block at a time.
  ! Both spellings are tried, since element ids are capitalised ('Fe') while
  ! molecule formulae keep the list's own case ('MgH').
  !===============================================================================
  FUNCTION linelist_species_wanted (id) RESULT (wanted)
    character(len=*), intent(in) :: id
    character(len=512) :: drop, only
    character(len=8) :: sym, cap
    logical :: wanted
    !-----------------------------------------------------------------------------
    wanted = .true.
    if (exclude_elements == '' .and. only_elements == '') return
    sym = adjustl(id)
    cap = capitalize(trim(sym))
    drop = ''
    only = ''
    drop(2:len(drop)-1) = trim_between_words(exclude_elements)
    only(2:len(only)-1) = trim_between_words(only_elements)
    if (index(drop, ' '//trim(sym)//' ') > 0 .or. &
        index(drop, ' '//trim(cap)//' ') > 0) then
      wanted = .false.
    else if (only_elements /= '') then
      wanted = index(only, ' '//trim(sym)//' ') > 0 .or. &
               index(only, ' '//trim(cap)//' ') > 0
    endif
  END FUNCTION linelist_species_wanted

  !===============================================================================
  ! Disable the line list for every species named in ids (element symbol or
  ! molecule formula), at the same per-element / per-molecule granularity as the
  ! departure_matches association (an element key covers all its ion stages).
  ! Used to strip the departure-carrying (NLTE) species out of tab_comp before
  ! building the sampled LTE line-list table, so those lines are instead handled
  ! in detail per cell with their departure coefficients.
  !===============================================================================
  SUBROUTINE mask_departure_species (self, ids)
    class(composition_type), target :: self
    character(len=*), intent(in), dimension(:) :: ids
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    integer :: i, di
    logical :: hit
    logical, save :: first_time = .true.
    ! ----------------------------------------------------------------------------
    do i=1, self%nel
      element => self%elements(i)
      hit = .false.
      do di=1, size(ids)
        if (uppercase(trim(element%id)) == uppercase(trim(ids(di)))) then
          hit = .true.; exit
        endif
      end do
      if (hit) then
        if (element_has_lines(element) .and. console_on(ch_hop, 1) .and. first_time) &
          write(console_unit,'(3a)') 'sample_linelist: ', trim(element%id), ' has departures -> excluded from LTE table'
        element%ions%in_linelist = .false.
      endif
    end do
    do i=1, self%nmolec
      mol => self%molecules(i)
      hit = .false.
      do di=1, size(ids)
        if (uppercase(trim(mol%id)) == uppercase(trim(ids(di)))) then
          hit = .true.; exit
        endif
      end do
      if (hit) then
        if (has_lines(mol) .and. console_on(ch_hop, 1) .and. first_time) &
          write(console_unit,'(3a)') 'sample_linelist: ', trim(mol%id), ' has departures -> excluded from LTE table'
        mol%in_linelist = .false.
      endif
    end do
    first_time = .false.
  END SUBROUTINE mask_departure_species

  !===============================================================================
  ! Whether a particle, or any ion stage of an element, holds line-list lines.
  !-------------------------------------------------------------------------------
  LOGICAL FUNCTION has_lines (p)
    class(particle_type), intent(in) :: p
    !-----------------------------------------------------------------------------
    has_lines = .false.
    if (allocated(p%ll)) has_lines = any(p%ll%nline > 0)
  END FUNCTION has_lines

  LOGICAL FUNCTION element_has_lines (element)
    class(element_type), intent(in) :: element
    integer :: j
    !-----------------------------------------------------------------------------
    element_has_lines = any([(has_lines(element%ions(j)), &
                              j=lbound(element%ions,1), ubound(element%ions,1))])
  END FUNCTION element_has_lines

  !===============================================================================
  ! init_linelist -- structure-of-arrays reader for the bulk linelist.
  ! Walks the file's per-species blocks through linelist_reader_mod, which addresses
  ! and bisects a block's records without reading what lies between them, and builds
  ! each particle's SoA linelist (particle%ll) via linelist_mod from the byte ranges
  ! the windows select. Hydrogen lines (Stark profiles) are read by add_h_block into
  ! the hydrogen extension of the H I store, which handles both the embedded H I
  ! block layout and the dedicated 'Hlinedata' column layout (selected by the
  ! hlinedata flag); its block is read whole and masked there.
  !-------------------------------------------------------------------------------
  SUBROUTINE init_linelist (self, linelist_file, dlam, lam, intervals)
    class(composition_type), target :: self
    class(particle_type), pointer :: particle
    character(len=*) :: linelist_file
    real(dp), intent(in) :: dlam
    real(dp), intent(in), dimension(:), optional :: lam
    type(interval_t), intent(in), dimension(:), optional :: intervals
    character(len=:), allocatable :: buf
    integer(8), allocatable :: lstart(:), lend(:)
    character(20) :: ion_id
    character(8) :: sym
    class(ion_type), pointer :: HI
    type(ll_file_t) :: llf
    type(ll_block_t) :: blk
    real(dp) :: lam_min, lam_max, vbar2e4
    integer :: iostat, nslice, nline, iso, iso_idx, nbefore, nkh, iv
    integer :: i, j, nactive, ntotal
    integer(8) :: b0, b1, b_after
    integer(8), allocatable :: r0(:), r1(:)
    logical :: isotope, ok, is_mol, lam_mask, hlinedata, use_intervals, have_block
    logical :: expect_lines
    ! atomic number of a Turbospectrum block (0 for a molecule), passed on to
    ! build_linelist; unallocated, it is an absent argument for a tumult block
    integer, allocatable :: ts_iel
    integer, save :: itimer=0, it_read=0, it_build=0, it_sort=0
    character(len=256) :: msg
    !-----------------------------------------------------------------------------
    if (linelist_file == '' .and. .not. self%ll_initialised) return

    ! every particle gets a per-isotope SoA linelist array ll(0:niso) (0 = main),
    ! before the hydrogen delegation so the diagnostic writer always sees an
    ! allocated array; idempotent so repeated calls (one per file) accumulate into it
    do i=1, self%nel
      do j=lbound(self%elements(i)%ions,1), ubound(self%elements(i)%ions,1)
        call ensure_ll_allocated(self%elements(i)%ions(j))
      end do
    end do
    do i=1, self%nmolec
      call ensure_ll_allocated(self%molecules(i))
    end do

    ! a dedicated 'Hlinedata' file uses the alternate H column layout (add_h_block)
    hlinedata = endswith(linelist_file, 'Hlinedata')
    call trace%begin ('composition_type%init_linelist', itimer=itimer)
    self%ll_initialised = .true.

    lam_mask = present(lam)
    use_intervals = present(intervals)
    if (lam_mask .and. .not. use_intervals) then
      lam_min = lam(1) - dlam
      lam_max = lam(size(lam)) + dlam
    else
      lam_min = 0
      lam_max = 0
    endif

    ! Walk the file block by block, reading only what each block contributes.
    ! linelist_reader_mod addresses a block's records without touching what lies
    ! between them, so a block of an unwanted species costs its two header lines
    ! and a narrow window costs its own records rather than the whole file. The
    ! per-interval slice is widened by one record on either side, and the mask
    ! below still decides membership in vacuum, so the selection is the one a
    ! reader of the whole file would make.
    call ll_open (llf, linelist_file)
    if (llf%unit == 0) then
      call trace%end(itimer)
      return
    end if

    if (use_intervals) allocate (r0(size(intervals)), r1(size(intervals)))

    call trace%begin ('init_linelist_build', itimer=it_build)
    ntotal = 0
    nactive = 0
    nkh = 0
    expect_lines = .false.
    do
      call ll_next_block (llf, blk, have_block)
      if (.not. have_block) exit
      nline = blk%nrec
      ntotal = ntotal + nline
      if (nline == 0) cycle
      if (allocated(ts_iel)) deallocate(ts_iel)
      if (blk%ts) then
        allocate(ts_iel)
        call ts_species (blk%siso, blk%ion, ion_id, ts_iel, iso, ok)
        if (.not. ok) then
          call console_warn (ch_hop, 'Turbospectrum species '//trim(blk%siso)//' not in the composition')
          cycle
        end if
        if (ts_iel > 0) call ensure_anstee_loaded()
      else
        read (blk%id, *, iostat=iostat) ion_id
        if (iostat /= 0) cycle
        ! the isotope token may carry Turbospectrum's species prefix (26.056)
        i = index(blk%siso, '.', back=.true.)
        read (blk%siso(i+1:), *, iostat=iostat) iso
        if (iostat /= 0) cycle
      end if
      ion_id = trim_between_words(ion_id)
      ! species switches are applied here rather than after the read, so a block of
      ! an unwanted species costs only its header (see linelist_species_wanted)
      read(ion_id, *) sym
      if (.not. linelist_species_wanted(sym)) cycle
      nullify(particle)
      call self%particle_by_id(ion_id(1:6), particle, verbose=.true., is_molecule=is_mol)
      ok = associated(particle)
      if (ok .and. is_mol .and. iso > 0) then
        select type (particle)
        class is (molec_type)
          iso = particle%canonical_iso(iso)
        end select
      end if
      isotope = iso /= 0

      iso_idx = 0
      if (ok .and. isotope) then
        if (associated(particle%isotope)) then
          iso_idx = findloc(particle%isotope%iso, iso, 1)
          ok = iso_idx > 0
        else
          ok = iso == nint(particle%awgt)
          isotope = .not. ok
        endif
      endif

      if (.not. ok) then
        if (isotope) then
          write(msg,'(2a,i0,a)') trim(ion_id(1:6)), ' isotope ', iso, ' not found'
          call console_warn (ch_hop, msg)
        endif
        cycle
      endif

      ! the bisection's own count of records in the requested windows, kept apart
      ! from the selection below so that an empty selection can be told from empty
      ! windows (see the check after the loop)
      if (.not. lam_mask) then
        expect_lines = .true.
      else if (use_intervals) then
        do iv = 1, size(intervals)
          if (ll_has_records (llf, blk, vac2air(intervals(iv)%blue_vac), &
                              vac2air(intervals(iv)%red_vac))) expect_lines = .true.
        end do
      else
        if (ll_has_records (llf, blk, vac2air(lam_min), vac2air(lam_max))) expect_lines = .true.
      end if

      if (particle%id == 'H I') then
        if (blk%ts) then
          call console_error ('init_linelist: '//trim(linelist_file)//' holds hydrogen lines in '// &
            'the Turbospectrum layout, which carries no Stark-profile data; read them from Hlinedata')
          error stop 'init_linelist: Turbospectrum hydrogen block'
        end if
        ! hydrogen keeps its Stark-profile kernel (add_h_opacities), and is read
        ! whole: its Stark wings reach far beyond any window margin
        call ll_whole_block (blk, b0, b1)
        call read_slice (b0, b1)
        nbefore = particle%ll(0)%nline
        call add_h_block(particle, buf, lstart(1:nslice), lend(1:nslice), nslice, &
                         lam_mask, lam_min, lam_max, intervals, use_intervals)
        nactive = nactive + (particle%ll(0)%nline - nbefore)
        cycle
      endif

      if (isotope) then
        vbar2e4 = particle%isotope(iso_idx)%compute_vbar2(1e4)
      else
        vbar2e4 = particle%compute_vbar2(1e4)
      endif

      particle%in_linelist = .true.
      nbefore = particle%ll(iso_idx)%nline
      if (.not. lam_mask) then
        call ll_whole_block (blk, b0, b1)
        call read_slice (b0, b1)
        call build_linelist(particle%ll(iso_idx), buf, lstart(1:nslice), &
          lend(1:nslice), nslice, particle%awgt, particle%j, particle%xion, is_mol, &
          vbar2e4, lam_mask, lam_min, lam_max, ts_iel=ts_iel)
      else if (use_intervals) then
        ! one byte range per interval, read into a single slice: windows far apart
        ! cost only their own records, while build_linelist still sees one ascending
        ! run of them and tests every record against every interval, as it does when
        ! it is handed the whole block.  b_after keeps two ranges from overlapping,
        ! so no record is added twice.
        b_after = blk%dstart
        do iv = 1, size(intervals)
          call ll_block_range (llf, blk, vac2air(intervals(iv)%blue_vac), &
                               vac2air(intervals(iv)%red_vac), b_after, r0(iv), r1(iv))
          b_after = max(b_after, r1(iv))
        end do
        call trace%begin ('init_linelist_read', itimer=it_read)
        call ll_read_ranges (llf, r0, r1, size(intervals), buf)
        call ll_index_records (buf, lstart, lend, nslice)
        call trace%end(it_read)
        if (nslice > 0) &
          call build_linelist(particle%ll(iso_idx), buf, lstart(1:nslice), &
            lend(1:nslice), nslice, particle%awgt, particle%j, particle%xion, is_mol, &
            vbar2e4, lam_mask, lam_min, lam_max, intervals, ts_iel=ts_iel)
      else
        call ll_block_range (llf, blk, vac2air(lam_min), vac2air(lam_max), &
                             blk%dstart, b0, b1)
        if (b1 > b0) then
          call read_slice (b0, b1)
          call build_linelist(particle%ll(iso_idx), buf, lstart(1:nslice), &
            lend(1:nslice), nslice, particle%awgt, particle%j, particle%xion, is_mol, &
            vbar2e4, lam_mask, lam_min, lam_max, ts_iel=ts_iel)
        end if
      end if
      nactive = nactive + (particle%ll(iso_idx)%nline - nbefore)
    end do
    call ll_close (llf)
    if (allocated(r0)) deallocate(r0, r1)
    if (allocated(buf)) deallocate(buf)
    if (allocated(lstart)) deallocate(lstart, lend)
    call trace%end(it_build)

    ! sort every touched particle's linelist by wavelength (march prerequisite)
    call trace%begin ('init_linelist_sort', itimer=it_sort)
    do i=1, self%nel
      do j=lbound(self%elements(i)%ions,1), ubound(self%elements(i)%ions,1)
        do iso_idx=0, ubound(self%elements(i)%ions(j)%ll, 1)
          call finalize_linelist(self%elements(i)%ions(j)%ll(iso_idx))
        end do
      end do
    end do
    do i=1, self%nmolec
      do iso_idx=0, ubound(self%molecules(i)%ll, 1)
        call finalize_linelist(self%molecules(i)%ll(iso_idx))
      end do
    end do
    call trace%end(it_sort)

    if (console_on(ch_hop, 1)) &
      write(console_unit,'(I0, a, I0, a, a)') nactive, ' of ', ntotal, &
        ' lines read from linelist ', trim(linelist_file)

    ! Records lie in the requested windows but none was selected: the selection is
    ! broken, and a run would go on to a line-free spectrum without saying so.
    ! Windows that hold no record of this list (a molecular list over a clean
    ! stretch, say) are legitimate and only noted.
    if (nactive == 0 .and. ntotal > 0) then
      if (expect_lines) then
        call console_error ('init_linelist: '//trim(linelist_file)//' has records in the '// &
          'requested wavelength windows, but none was selected')
        error stop 'init_linelist: no lines selected from a list that covers the windows'
      else
        call console_warn (ch_hop, 'no line of '//trim(linelist_file)//' lies in the requested windows')
      end if
    end if
    self%nline = nactive
    call self%ion_by_id('H I', HI)
    if (associated(HI)) HI%nk = max(HI%nk, nkh)
    call self%mask_linelist()
    call trace%end(itimer)

  contains
    ! A Turbospectrum species code: Z.AAA for an atom, or two-digit atomic numbers in
    ! ascending order followed by their three-digit masses for a molecule (0607.012014
    ! is 12C14N, 000 the natural mixture).  Gives the composition id, the atomic
    ! number (0 for a molecule) and the mass part as an isotope label.
    subroutine ts_species (code, ion, id, iel, label, found)
      character(len=*), intent(in) :: code
      integer, intent(in) :: ion
      character(len=*), intent(out) :: id
      integer, intent(out) :: iel, label
      logical, intent(out) :: found
      character(len=32) :: head, tail
      integer :: dot, k, m, nz, ios
      integer, allocatable :: z(:), order(:)
      !-------------------------------------------------------------------------
      found = .false.
      id = ''
      iel = 0
      dot = index(code, '.')
      if (dot == 0) dot = len_trim(code) + 1
      head = adjustl(code(:dot-1))
      tail = code(dot+1:)
      label = 0
      ios = 0
      if (tail /= '') read (tail, *, iostat=ios) label
      if (ios /= 0) return
      if (len_trim(head) <= 2) then
        read (head, *, iostat=ios) iel
        if (ios /= 0 .or. iel < 1 .or. iel > size(element_ids) .or. &
            ion < 1 .or. ion > ubound(roman, 1)) return
        id = trim(element_ids(iel))//roman(ion)
        found = .true.
        return
      end if
      if (ion /= 1) return
      if (mod(len_trim(head), 2) == 1) head = '0'//head
      nz = len_trim(head) / 2
      allocate (z(nz), order(nz))
      do k = 1, nz
        read (head(2*k-1:2*k), *, iostat=ios) z(k)
        if (ios /= 0) return
      end do
      call argsort_int (z, order)
      z = z(order)
      do m = 1, self%nmolec
        associate (mol => self%molecules(m))
          if (mol%nions /= nz .or. scan(mol%id, '+-') > 0) cycle
          call argsort_int (mol%z, order)
          if (all(mol%z(order) == z)) then
            id = mol%id
            found = .true.
            return
          end if
        end associate
      end do
    end subroutine ts_species

    ! the byte range [b0,b1) of the open list as buf, with its record boundaries;
    ! the offsets are into the slice, which is what build_linelist and add_h_block
    ! index with
    subroutine read_slice (b0, b1)
      integer(8), intent(in) :: b0, b1
      !-------------------------------------------------------------------------
      call trace%begin ('init_linelist_read', itimer=it_read)
      call ll_read_slice (llf, b0, b1, buf)
      call ll_index_records (buf, lstart, lend, nslice)
      call trace%end(it_read)
    end subroutine read_slice

    ! allocate the per-isotope SoA linelist array ll(0:niso) once (0 = main isotope)
    subroutine ensure_ll_allocated(p)
      class(particle_type), intent(inout) :: p
      integer :: niso
      !-------------------------------------------------------------------------
      if (allocated(p%ll)) return
      niso = 0
      if (associated(p%isotope)) niso = size(p%isotope)
      allocate(p%ll(0:niso))
    end subroutine ensure_ll_allocated

    ! hydrogen block -> particle%ll(0), with the principal quantum numbers that select
    ! the Stark profile (append_h_line). Two column layouts are accepted, selected by
    ! hlinedata:
    !   hlinedata=.true.  (dedicated 'Hlinedata' file):
    !       lam_air, i, j, lower_ev, upper_ev, loggf, label, grad_in
    !     -- principal quantum numbers i,j are read directly from the file.
    !   hlinedata=.false. (H I block embedded in a metal file):
    !       lam_air, lower_ev, loggf, gvdw, grad, gstk, i, j
    !     -- i,j are the principal quantum numbers when both are given, else derived
    !        from the level energies with the ionisation energy of hydrogen.
    ! The broadening columns are inert: hydrogen_line supplies the widths.
    subroutine add_h_block(particle, buf, lstart, lend, n, mask, lam_min_arg, lam_max_arg, &
                           intervals_arg, use_intervals_arg)
      class(particle_type), pointer :: particle
      character(len=*), intent(in) :: buf
      integer(8), intent(in) :: lstart(:), lend(:)
      integer, intent(in) :: n
      logical, intent(in) :: mask
      real(dp), intent(in) :: lam_min_arg, lam_max_arg
      type(interval_t), intent(in), optional :: intervals_arg(:)
      logical, intent(in), optional :: use_intervals_arg
      real(dp) :: lam_air_h, lower_ev_h, upper_ev_h, loggf_h, gv_h, ga_h, gq_h, r_i, r_j
      real(dp) :: lam0v, nu0, upper_ev
      character(len=150) :: label
      integer :: ii, ios, i_h, j_h
      logical :: use_intervals_local
      !-------------------------------------------------------------------------
      use_intervals_local = .false.
      if (present(use_intervals_arg)) use_intervals_local = use_intervals_arg
      if (hydrogen_vcs_table) call ensure_vcs_loaded()

      do ii=1, n
        gv_h = 0; ga_h = 0; gq_h = 0; r_i = 0; r_j = 0
        if (hlinedata) then
          read(buf(lstart(ii):lend(ii)), *, iostat=ios) lam_air_h, i_h, j_h, &
                                          lower_ev_h, upper_ev_h, loggf_h, label, ga_h
        else
          read(buf(lstart(ii):lend(ii)), *, iostat=ios) lam_air_h, lower_ev_h, loggf_h, &
                                          gv_h, ga_h, gq_h, r_i, r_j
        endif
        ! ios < 0 is a record that ends before the optional trailing columns
        if (ios > 0) then
          call console_error ('cannot read hydrogen line: '//buf(lstart(ii):lend(ii)))
          error stop
        endif
        lam0v = air2vac(lam_air_h)
        if (mask) then
          if (use_intervals_local .and. present(intervals_arg)) then
            if (.not. in_any_interval(lam0v, intervals_arg)) cycle
          else
            if (lam0v < lam_min_arg .or. lam0v > lam_max_arg) cycle
          end if
        endif
        nu0 = HZ_TO_AA / lam0v
        if (.not. hlinedata) then
          if (r_i > 0 .and. r_j > 0) then
            i_h = nint(r_i)
            j_h = nint(r_j)
          else
            upper_ev = lower_ev_h + H_E*nu0
            i_h = nint(sqrt(RYD / (RYD-lower_ev_h)))
            j_h = nint(sqrt(RYD / (RYD-upper_ev)))
          endif
        endif
        call append_h_line(particle%ll(0), log(lam0v), nu0, lower_ev_h, loggf_h, &
                           gv_h, ga_h, gq_h, i_h, j_h)
        nkh = max(j_h, nkh)
        particle%in_linelist = .true.
      end do
    end subroutine add_h_block

  END SUBROUTINE init_linelist

  !===============================================================================
  ! add_linelist_opacities -- line opacity orchestrator. The H I lines use the
  ! Stark-profile kernel add_h_opacities; metals and molecules the linelist_mod kernel.
  !-------------------------------------------------------------------------------
  SUBROUTINE add_linelist_opacities(self, sp, all_chi, all_eta, dlam, vproj, departures, &
                                    only_departures, one_m_eps, chi_scat, beta, i0, i1, llam_in, &
                                    i_set, window_lam, chi_lines, eta_lines, lte_of)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(molec_type), pointer :: mol
    class(departure_pointer_t), pointer :: dep
    type(hline_t) :: hline
    class(departure_pointer_t), intent(in), dimension(:), target, optional :: departures
    ! only_departures=.true. skips every particle that has no matching departure
    ! set, i.e. processes only NLTE species -- the complement of the masked
    ! tab_comp used to build the sampled LTE line table (sample_linelist).
    logical, intent(in), optional :: only_departures
    ! lte_of: with only_departures, also processes the species that these departure
    ! sets belong to, in LTE (no departures passed): the LTE companion set of a run
    ! whose sampled table leaves them out.
    class(departure_pointer_t), intent(in), dimension(:), optional :: lte_of
    ! line_scattering: precomputed scattering fractions 1-eps per stage group
    ! (last index 1=neutral incl. molecules, 2=ionized) and the accumulator the
    ! (1-eps)-weighted LTE line opacity is added to; hydrogen stays thermal
    ! chi_scat must enter zeroed.  The kernels accumulate it in working units
    ! (pre-stimulated-emission, like all_chi); with beta present it is converted
    ! to physical 1/cm on the way out, so callers never see the raw form.
    real(4), intent(in),    dimension(:,:,:,:,:), optional :: one_m_eps
    real(4), intent(inout), dimension(:,:,:,:),   optional :: chi_scat
    real(4), intent(in),    dimension(:,:,:,:),   optional :: beta
    integer, intent(in), optional :: i0, i1
    integer, intent(in), optional :: i_set
    ! window_lam: [blue, red] wavelength in Angstrom of the window_file window this
    ! output patch lies in (see below)
    real(dp), intent(in), optional :: window_lam(2)
    ! chi_lines, eta_lines: receive the line increments added to all_chi and all_eta
    ! over i0:i1, so a caller can add them again later with the same result
    real(4), intent(out), dimension(:,:,:,:), optional :: chi_lines, eta_lines
    class(spectrum_type) :: sp
    real(4), intent(inout), dimension(:,:,:,:) :: all_chi, all_eta
    real(8), intent(in) :: dlam
    real(8) :: llam_lines(2)
    real(4), intent(in), dimension(:,:,:), contiguous, target, optional :: vproj
    real(4), dimension(:,:,:), allocatable :: tothi, stark_pe, ne, temp, e_kT
    real(4), dimension(:,:,:), allocatable :: inv_df, vbar2, fo, t3nhe, ni, nj
    real(4), dimension(:,:,:), allocatable :: vmic_c
    real(4), dimension(:,:,:,:), allocatable, target :: chi_l, eta_l
    real(4), dimension(:,:,:,:), allocatable :: dummy_dep, logref
    real(4), dimension(:,:,:), pointer, contiguous :: vpr
    real(4), dimension(:,:,:), pointer :: bi, bj
    real(dp), dimension(:,:,:), pointer :: set_n
    real(dp) :: set_scale
    real(8), dimension(:), allocatable, target :: llam_mem
    real(8), dimension(:), pointer :: llam
    real(8), intent(in), dimension(:), optional, target :: llam_in
    real(4) :: log_tol_over_voigt_peak
    real(dp) :: lam0, nu0, Bij
    integer :: i, j, k, l, f, f0, m(3), nf, lo, hi, io, nn, set_index, n_lo, n_up
    logical :: departures_present, do_nlte, skip_lte
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('composition_type%add_linelist_opacities', itimer=itimer)
    m = shape(self%ne)
    nf = sp%nnu
    lo = 1;  hi = m(1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    nn = hi - lo + 1
    set_index = 0
    if (present(i_set)) set_index = i_set
    allocate(tothi(nn,m(2),m(3)))
    allocate(stark_pe, ne, temp, e_kT, inv_df, vbar2, fo, t3nhe, &
             ni, nj, mold=tothi)
    allocate(chi_l(nf, nn, m(2), m(3)), eta_l(nf, nn, m(2), m(3)))
    allocate(dummy_dep(1,1,1,1), source=0.)
    chi_l = 0.
    eta_l = 0.

    departures_present = present(departures)
    skip_lte = .false.
    if (present(only_departures)) skip_lte = only_departures
    if (present(vproj)) then
      allocate(vpr, source=vproj(lo:hi,:,:))
    else
      allocate(vpr, mold=tothi)
      vpr = 0
    endif

    ne   = real(self%ne(lo:hi,:,:))
    temp = real(self%temp(lo:hi,:,:))
    e_kT = real(self%e_kT(lo:hi,:,:))
    tothi = real(self%tothi(lo:hi,:,:))
    stark_pe = real(self%stark_pe(lo:hi,:,:))
    allocate(vmic_c, source=self%vmic(lo:hi,:,:))

    if (present(llam_in)) then
      llam => llam_in
    else
      allocate(llam_mem(nf))
      llam_mem = log(sp%lam)
      llam => llam_mem
    endif

    ! This patch's line neighbourhood in log(lambda), handed to the kernels. dlam is the
    ! selection half-width in Angstrom, so dividing by the bluest wavelength gives the
    ! widest log-space equivalent. Only a segmented (multi-window) grid needs this: there
    ! the shared line list is disjoint, and the forward march would otherwise clamp another
    ! window's lines onto this patch's edge. A patch covering part of a window takes the
    ! neighbourhood of the whole window, so splitting a window into several patches visits
    ! the same lines as one patch would. A contiguous grid keeps the whole-list visit.
    if (present(window_lam)) then
      llam_lines = log(window_lam) + [-1d0, 1d0] * dlam / window_lam(1)
    else if (window_file /= '') then
      llam_lines = [llam(1) - dlam / exp(llam(1)), llam(nf) + dlam / exp(llam(1))]
    else
      llam_lines = [-huge(1d0), huge(1d0)]
    end if

    ! Precompute once here (shared by all particles); guard non-positive continuum.
    allocate(logref(nf, nn, m(2), m(3)))
    log_tol_over_voigt_peak = log(tolerance / voigt_peak)
    where (all_chi(:, lo:hi, :, :) > 0.)
      logref = log_tol_over_voigt_peak + log(all_chi(:, lo:hi, :, :))
    elsewhere
      logref = -huge(1.)
    end where

    ! --- hydrogen: Stark profiles (add_h_opacities) on the H I line list ---
    ! n_lo/n_up select the profile; in NLTE lev_i/lev_j are the departure-set levels
    ! (departure_assign_h_levels), whose coefficients scale the LTE populations of
    ! levels n_lo/n_up.
    element => self%elements(1)
    ion => element%ions(1)
    call get_departure_coefficients(ion, dep)
    if (ion%in_linelist .and. allocated(ion%ll) .and. &
        .not. (skip_lte .and. .not. associated(dep) .and. .not. in_lte_of(ion))) then
      associate (hll => ion%ll(0))
      if (hll%nline > 0) then
        fo = self%fo(lo:hi,:,:)
        vbar2 = self%vbar2h(lo:hi,:,:)
        t3nhe = self%t3nhe(lo:hi,:,:)
        inv_df = element%inverse_doppler_factor(temp, self%vmic(lo:hi,:,:))
      endif

      if (associated(dep)) call dep%assign_h_levels(hll)
      do l=1, hll%nline
        if (.not. hll%active(l)) cycle
        lam0 = exp(hll%llam0(l))
        if (dlam == 0) then
          f0 = findloc(llam, hll%llam0(l), 1)
        else
          f0 = find_closest(sp%lam, lam0)
          if (abs(sp%lam(f0) - lam0) > dlam) f0 = 0
        endif
        if (f0 == 0) cycle
        nu0 = HZ_TO_AA / lam0
        n_lo = hll%n_lo(l)
        n_up = hll%n_up(l)
        ni = ion%levels(n_lo)%n(lo:hi,:,:)
        nj = ion%levels(n_up)%n(lo:hi,:,:)
        call self%set_population(ion%id, set_index, set_n, set_scale)
        if (associated(set_n)) then
          where (ion%n(lo:hi,:,:) > 0d0)
            ni = ni * real(set_n(lo:hi,:,:)/ion%n(lo:hi,:,:),4)
            nj = nj * real(set_n(lo:hi,:,:)/ion%n(lo:hi,:,:),4)
          endwhere
        else
          ni = ni * real(set_scale,4)
          nj = nj * real(set_scale,4)
        endif
        if (H_occupation_probability) then
          ni = ni * ion%levels(n_up)%w(lo:hi,:,:)
          nj = nj * ion%levels(n_up)%w(lo:hi,:,:)
        endif
        do_nlte = associated(dep) .and. hll%lev_i(l)>0 .and. hll%lev_j(l)>0
        if (do_nlte) then
          bi => dep%mem(:,:,:,hll%lev_i(l))
          bj => dep%mem(:,:,:,hll%lev_j(l))
        endif
        hline = hydrogen_line(n_lo, n_up, nu0)
        ! f = gf / g_lower, with g_lower = 2 n_lo**2
        Bij = bij_const * 10**dble(hll%loggf(l)) / (2*n_lo**2 * nu0)
        call add_h_opacities(hline, Bij, f0, sp%nu, vpr, inv_df, temp, &
          ni, nj, ne, tothi, vbar2, fo, t3nhe, bi, bj, chi_l, eta_l, all_chi, do_nlte, &
          hll%used(l), hll%strength(l), lo, hi)
      end do
      end associate
    endif

    ! --- metals and molecules: SoA kernel ---
    do i=2, self%nel
      element => self%elements(i)
      if (.not. any(element%ions%in_linelist)) cycle
      do j=lbound(element%ions,1), ubound(element%ions,1)
        call particle_add_opacities_soa(element%ions(j))
      end do
    end do
    do i=1, self%nmolec
      call particle_add_opacities_soa(self%molecules(i))
    end do

    do k=1, m(3)
      do j=1, m(2)
        do i=lo, hi
          io = i - lo + 1
          do f=1, nf
            all_eta(f,i,j,k) = all_eta(f,i,j,k) + eta_l(f,io,j,k)
            all_chi(f,i,j,k) = all_chi(f,i,j,k) + chi_l(f,io,j,k)
          end do
        end do
      end do
    end do
    if (present(chi_lines)) chi_lines = chi_l
    if (present(eta_lines)) eta_lines = eta_l
    deallocate(vpr)

    ! working units -> physical extinction for the scattering component
    if (present(chi_scat) .and. present(beta)) then
      do k=1, m(3)
        do j=1, m(2)
          do i=lo, hi
            do f=1, nf
              chi_scat(f,i,j,k) = chi_scat(f,i,j,k) * (1 - beta(f,i,j,k))
            end do
          end do
        end do
      end do
    endif

    if (ll_profile) then
      write(console_unit,'(a,3(a,i0),a,f7.2,a,f6.2)') 'll_profile (cumulative):', &
        ' visit=', ll_nvisit, ' pass=', ll_npass, ' voigt=', ll_nvoigt, &
        ' | voigt/pass=', real(ll_nvoigt)/max(ll_npass,1_8), &
        ' pass%=', 100.*real(ll_npass)/max(ll_nvisit,1_8)
    endif
    call trace%end(itimer)

  contains
    ! --------------------------------------------------------------------------
    subroutine get_departure_coefficients(ptc, dep)
      class(particle_type), target :: ptc
      class(departure_pointer_t), pointer :: dep
      integer :: di
      !-------------------------------------------------------------------------
      dep => null()
      if (.not. departures_present) return

      do di=1, size(departures)
        if (departures(di)%matches(ptc)) then
          dep => departures(di)
          exit
        endif
      end do
    end subroutine get_departure_coefficients

    ! --------------------------------------------------------------------------
    logical function in_lte_of(ptc)
      class(particle_type), target :: ptc
      integer :: di
      !-------------------------------------------------------------------------
      in_lte_of = .false.
      if (.not. present(lte_of)) return
      do di=1, size(lte_of)
        if (lte_of(di)%matches(ptc)) then
          in_lte_of = .true.
          return
        endif
      end do
    end function in_lte_of

    ! --------------------------------------------------------------------------
    subroutine particle_add_opacities_soa(ptc)
      class(particle_type), target :: ptc
      !-------------------------------------------------------------------------
      if (.not. ptc%in_linelist) return
      if (.not. allocated(ptc%n)) error stop 'particle density not allocated'

      call get_departure_coefficients(ptc, dep)
      ! LTE species live in the sampled table, except those named in lte_of
      if (skip_lte .and. .not. associated(dep) .and. .not. in_lte_of(ptc)) return

      ! the per-isotope population base / Doppler width construction and the kernel
      ! dispatch now live in linelist_mod (add_particle_opacity); we own only the patch
      ! work arrays, the departure lookup and the level assignment, and hand dep%mem in.
      call self%set_population(ptc%id, set_index, set_n, set_scale)
      if (associated(dep)) then
        ! NLTE species: departure-coefficient source, no scattering split
        call dep%assign_levels(ptc)
        call add_particle_opacity(ptc, temp, vmic_c, llam, chi_l, eta_l, &
          all_chi, logref, vpr, e_kT, ne, tothi, stark_pe, dep%mem, &
          .true., sp%equidistant, i0=lo, i1=hi, &
          population=set_n, population_scale=set_scale, llam_lines=llam_lines)
      elseif (present(one_m_eps) .and. present(chi_scat)) then
        ! LTE species with the stage-appropriate scattering fraction
        ! (molecules keep the default j=1 -> neutral p-function)
        call add_particle_opacity(ptc, temp, vmic_c, llam, chi_l, eta_l, &
          all_chi, logref, vpr, e_kT, ne, tothi, stark_pe, dummy_dep, &
          .false., sp%equidistant, &
          one_m_eps=one_m_eps(:,:,:,:, merge(1, 2, ptc%j <= 1)), chi_scat=chi_scat, &
          i0=lo, i1=hi, population=set_n, population_scale=set_scale, llam_lines=llam_lines)
      else
        call add_particle_opacity(ptc, temp, vmic_c, llam, chi_l, eta_l, &
          all_chi, logref, vpr, e_kT, ne, tothi, stark_pe, dummy_dep, &
          .false., sp%equidistant, i0=lo, i1=hi, &
          population=set_n, population_scale=set_scale, llam_lines=llam_lines)
      endif
    end subroutine particle_add_opacities_soa
  END SUBROUTINE add_linelist_opacities

  !===============================================================================
  ! Accumulate the line-list diagnostics (used, strength, reach) of src, a copy of
  ! self made after init_linelist, into self. The kernels write them per copy, which
  ! keeps the hot loops free of shared writes; merging every copy into the composition
  ! they were made from gives the maximum over all of them.
  !-------------------------------------------------------------------------------
  SUBROUTINE merge_line_diagnostics(self, src)
    class(composition_type), target :: self
    class(composition_type), intent(in) :: src
    integer :: i, j
    !-----------------------------------------------------------------------------
    do i=1, self%nel
      do j=lbound(self%elements(i)%ions,1), ubound(self%elements(i)%ions,1)
        call merge_particle(self%elements(i)%ions(j), src%elements(i)%ions(j))
      end do
    end do
    do i=1, self%nmolec
      call merge_particle(self%molecules(i), src%molecules(i))
    end do
  contains
    subroutine merge_particle (dst_p, src_p)
      class(particle_type), intent(inout) :: dst_p
      class(particle_type), intent(in) :: src_p
      integer :: k
      !---------------------------------------------------------------------------
      if (.not. allocated(dst_p%ll) .or. .not. allocated(src_p%ll)) return
      do k=lbound(dst_p%ll,1), ubound(dst_p%ll,1)
        call merge_ll_diag(dst_p%ll(k), src_p%ll(k))
      end do
    end subroutine merge_particle
  END SUBROUTINE merge_line_diagnostics

  !===============================================================================
  SUBROUTINE reassign_pointer(self)
    class(composition_type), target :: self
    class(element_type),     pointer :: element
    class(molec_type),       pointer :: molecule
    class(ion_type),         pointer :: ion
    class(level_type),       pointer :: lvl
    class(bftrans_type),  pointer :: trans
    integer :: i, j, m
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%reassign_pointer')

    do i=1, self%nel
      element => self%elements(i)

      do j=lbound(element%ions,1), ubound(element%ions,1)
        ion => element%ions(j)
        ion%abund => element%abund

        if (j > lbound(element%ions,1)) then
          ion%below => element%ions(j-1)
        endif
        if(j < ubound(element%ions,1)) then
          ion%above => element%ions(j+1)
        endif

        do m=1, ion%nmolec
          call self%molecule_by_id(ion%molecules(m)%p%id, molecule)
          ion%molecules(m)%p => molecule
        end do
      end do
    end do

    do i=1, size(self%molecules)
      molecule => self%molecules(i)

      do j=1, molecule%nions
        ion => molecule%ions(j)%p
        call self%ion_by_id(ion%id, ion)
        molecule%ions(j)%p => ion
      end do
    end do

    do i=1, size(self%trans)
      trans => self%trans(i)
      if (associated(trans%ion)) then
        if (trans%ion%id(1:2) == 'Me') then
          trans%ion => self%peachmet
        else
          call self%particle_by_id(trans%ion%id, trans%ion)
          if (associated(trans%ion2)) then
            call self%particle_by_id(trans%ion2%id, trans%ion2)
          endif
        endif
      endif
    end do
    call trace%end ()
  END SUBROUTINE reassign_pointer

  !===============================================================================
  SUBROUTINE load_state_variables(self, temp, ne, rho, use_rho, use_ne, vmic, vx, vy, vz, &
                                  set_da, set_element)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(molec_type), pointer :: mol
    type(spectrum_type) :: ref_sp
    real(dp), intent(in), dimension(:,:,:) :: temp, ne, rho
    logical, intent(in), optional :: use_rho, use_ne
    real(4), intent(in), dimension(:,:,:), optional :: vmic, vx, vy, vz
    real(8), intent(in), dimension(:,:), optional :: set_da
    logical, intent(in), dimension(:,:), optional :: set_element
    real(dp), dimension(:,:,:), allocatable :: ro
    real(4), dimension(:,:,:,:), pointer :: chi500
    real(dp) :: nudebug(1)
    integer :: i,j, m(3)
    integer :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%load_state_variables', itimer=itimer)

    if (.not. self%initialised) then
      call console_error ('load_state_variables: composition has not been initialised yet')
      error stop
    endif

    if (.not. (option(use_rho, .true.) .or. option(use_ne, .true.))) then
      call console_error ('load_state_variables needs at least one of rho or ne')
      error stop
    endif

    ! NOTE: this cleanup is partial (chi500/fo/vbar2h/t3nhe are not freed, and
    ! ion_type retains partial state after the deallocate loop below), so a
    ! second call on the same instance aborts.  Every caller -- dispatch and
    ! the standalone hop_c_api -- copies a never-loaded composition and loads
    ! exactly once per copy.
    if (self%loaded) then
      deallocate(self%temp)
      deallocate(self%ne, self%rho, self%pe, self%pg, self%E, self%eh)
      deallocate(self%toth, self%teta, self%dxi, self%lsaha)
      deallocate(self%e_kT, self%p6, self%stark_pe, self%tothi, self%zero)
      if (allocated(self%vmic)) deallocate(self%vmic)
      if (allocated(self%vx)) deallocate(self%vx, self%vy, self%vz)
    endif

    allocate(self%temp, self%ne, self%rho, self%pe, self%pg, self%E, self%eh, mold=temp)
    self%temp = -1 ! sentinel: never matches a physical temp in same_temp
    self%same_temp = .false.
    allocate(self%toth, self%teta, self%dxi, self%lsaha, self%fo, self%vbar2h, self%t3nhe, mold=temp)
    allocate(self%e_kT, self%p6, ro, self%stark_pe, self%tothi, self%zero, mold=temp)
    if (present(vmic)) allocate(self%vmic, source=vmic)
    if (present(vx)) allocate(self%vx, source=vx)
    if (present(vy)) allocate(self%vy, source=vy)
    if (present(vz)) allocate(self%vz, source=vz)

    self%zero = 0.d0
    self%E = 0

    self%ne = ne
    if (option(use_ne, .true.)) then
      call self%thermal_equilibrium(temp, ne, self%rho)
    else
      call self%iterate_thermal_eq(temp, rho, self%ne, tol=1d-6)
      call self%thermal_equilibrium(temp, self%ne, self%rho)
    endif

    if (option(use_rho, .true.)) then
      self%rho = rho
    endif

    self%toth = self%rho / self%grph

    do i=1, self%nel
      element => self%elements(i)
      do j=lbound(element%ions,1), ubound(element%ions,1)
        ion => element%ions(j)
        ion%n = ion%f * self%toth
      end do
    end do

    do i=1, self%nmolec
      mol => self%molecules(i)
      mol%n = mol%f * self%toth
    end do

    call self%set_aux_data(temp, self%ne)
    call self%compute_energy_and_gas_pressure(self%eh, self%E, self%pg)
    call self%compute_damping_aux()
    call self%init_H_levels()

    if (present(set_da) .and. present(set_element)) &
      call self%prepare_abundance_sets(set_da, set_element)

    do i=1, self%nel
      element => self%elements(i)
      do j=lbound(element%ions,1), ubound(element%ions,1)
        ion => element%ions(j)
        call ion%deallocate()
      end do
    end do

    do i=1, self%nmolec
      mol => self%molecules(i)
      call mol%deallocate()
    end do

    m = shape(self%temp)
    allocate(self%chi500(m(1),m(2),m(3)))
    chi500(1:1,1:m(1),1:m(2),1:m(3)) => self%chi500
    call ref_sp%reset(lam=[air2vac(ref_lam)])
    ! chi500 is the *total* extinction: it is the tau500 scale and the
    ! normalisation of every sampled opacity table
    call self%opacity(ref_sp, chi500, include_absmet=.false., include_llist=.false., &
      sorted=.true.)

    if (debug_opac > 0) then
      nudebug = HZ_TO_AA / air2vac(debug_opac)
      call self%print_opacity_contributions(nudebug)
    endif

    self%loaded = .true.
    call trace%end(itimer)
  END SUBROUTINE load_state_variables

  !=========================================================================
  !> Store the ln(T) and Pe gradients captured from the snapshot interpolant
  !> (per code length unit; the Pe gradients are of Pe itself, not its log)
  !> and cache a real(4) copy of Pe, so the chain-rule RT solvers never
  !> convert precision in hot code.  Must be called after
  !> load_state_variables, which finalizes self%pe.  The horizontal partials
  !> are optional -- rotated compositions carry only the z pair.
  !=========================================================================
  SUBROUTINE set_gradients(self, dlnt_dz, dpe_dz, dlnt_dx, dlnt_dy, dpe_dx, dpe_dy)
    class(composition_type) :: self
    real(4), intent(in), dimension(:,:,:) :: dlnt_dz, dpe_dz
    real(4), intent(in), dimension(:,:,:), optional :: dlnt_dx, dlnt_dy, dpe_dx, dpe_dy
    !---------------------------------------------------------------------------
    if (.not. self%loaded) then
      call console_error ('set_gradients must be called after load_state_variables')
      error stop
    endif
    if (allocated(self%dlnt_dz)) deallocate(self%dlnt_dz, self%dpe_dz)
    if (allocated(self%dlnt_dx)) deallocate(self%dlnt_dx, self%dlnt_dy)
    if (allocated(self%dpe_dx))  deallocate(self%dpe_dx, self%dpe_dy)
    if (allocated(self%pe4))     deallocate(self%pe4)
    allocate(self%dlnt_dz, source=dlnt_dz)
    allocate(self%dpe_dz,  source=dpe_dz)
    if (present(dlnt_dx)) allocate(self%dlnt_dx, source=dlnt_dx)
    if (present(dlnt_dy)) allocate(self%dlnt_dy, source=dlnt_dy)
    if (present(dpe_dx))  allocate(self%dpe_dx,  source=dpe_dx)
    if (present(dpe_dy))  allocate(self%dpe_dy,  source=dpe_dy)
    allocate(self%pe4(size(dlnt_dz,1), size(dlnt_dz,2), size(dlnt_dz,3)))
    self%pe4 = real(self%pe, 4)
  END SUBROUTINE set_gradients

  !=========================================================================
  SUBROUTINE thermal_equilibrium(self, temp, ne, rho)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion, HI, HeI
    real(dp), intent(in),  dimension(:,:,:) :: temp, ne
    real(dp), intent(out), dimension(:,:,:) :: rho
    real(dp), dimension(:,:,:), allocatable :: fe, xnen, xnecno, toth, hpart
    type(ion_pointer_t), dimension(:), allocatable :: mions
    integer :: i, j, nions, itimer=0, itimer2=0, itimer3=0, x, y, z
    logical :: converged
    logical, save :: hocc_warned = .false.
    character(len=256) :: msg
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%thermal_equilibrium', itimer=itimer)
    if ((mhd_eos .neqv. H_occupation_probability) .and. .not. hocc_warned) then
      hocc_warned = .true.
      call console_warn (ch_hop, 'mhd_eos and H_occupation_probability differ;'// &
                   ' EOS and opacity treat H level dissolution inconsistently')
    endif
    allocate(fe, xnen, xnecno, toth, mold=temp)
    call self%set_aux_data(temp, ne)
    call self%ion_by_id('H I', HI)
    call self%ion_by_id('He I', HeI)

    fe = 0
    xnecno = 0
    if (mhd_eos) then
      allocate(hpart, mold=self%temp)
      HI%state_dependent_part = .true.   ! veto temp-only moleq dis caching
      ! no Debye lowering on the H I -> H II step: the charged-perturber factor
      ! in wcalc is itself a continuum-lowering model, and applying both
      ! double-counts.  H- keeps its Debye term (dxi_fac=2): its binding is not
      ! described by the H-atom occupation weights.
      HI%dxi_fac = 0
    endif

    do i=1, self%nel
      element => self%elements(i)
      call trace%begin ('composition_type%partition_functions', itimer=itimer2)
      call element%compute_partition_functions(self%temp, self%e_kT, self%ne)
      call trace%end(itimer2)

      call trace%begin ('composition_type%ion_equilibrium', itimer=itimer3)
      if (mhd_eos .and. i == 1) then
        ! hydrogen: occupation-weighted partition from the previous call as a
        ! warm start, refined after molecular equilibrium below
        if (allocated(self%part_hocc)) then
          if (all(shape(self%part_hocc) == shape(self%temp))) HI%part = self%part_hocc
        endif
      endif
      call element%compute_ion_equilibrium(self%e_kT, self%lsaha, self%dxi)
      call trace%end(itimer3)

      do j=lbound(element%ions, 1), ubound(element%ions, 1)
        ion => element%ions(j)

        ! number of electrons from ion
        xnen = ion%f * (j-1)
        fe = fe + xnen

        if (all(element%ions%nmolec == 0)) then
          ! number of electrons from elements without molecules
          xnecno = xnecno + xnen
        endif
      end do
    end do

    if (self%nmolec > 0) then
      ! solve elements: the neutral root of every molecule-forming element
      nions = 0
      do i=1, self%nel
        if (self%elements(i)%ions(1)%nmolec > 0) nions = nions + 1
      end do
      allocate(mions(nions))
      nions = 0
      do i=1, self%nel
        if (self%elements(i)%ions(1)%nmolec > 0) then
          nions = nions + 1
          mions(nions)%p => self%elements(i)%ions(1)
        endif
      end do
      call molecular_equilibrium(self%molecules, mions, self%tmolim, self%temp, &
                                 self%teta, self%e_kT, self%pe, xnecno, fe, &
                                 ysave=self%moleq_seed, same_state=self%same_temp)
    endif

    ! total hydrogen number density, computed from scratch
    toth = ne / fe
    ! mass density computed from scratch
    rho = toth * self%grph

    if (mhd_eos) then
      ! Feed the Hummer-Mihalas occupation-weighted H I partition function into
      ! the ionization balance.  Overwriting HI%part (rather than a dxi shim)
      ! makes H-, the Saha ladder, and the molecular equilibrium (get_dis uses
      ! ion%part) all see the same dissolved partition.  The update is
      ! evaluated here, after molecular equilibrium, where toth is physical
      ! (H-, H2 and metal donors balanced; the pre-moleq fe can pass through
      ! zero at cool dense states).  If the partition moved, hydrogen and the
      ! molecules are re-solved; with the part_hocc warm start installed in the
      ! Saha loop above the tolerance is met on the first check in all but the
      ! first call, so the re-solve is cold-start only.
      do i=1, 6
        HI%n  = max(HI%f  * toth, 0d0)
        HeI%n = max(HeI%f * toth, 0d0)
        call self%Hoccupation_partition(HI%n, HeI%n, hpart)
        hpart = log(max(hpart, 1d-40))
        ! tolerance well below the ne-secant's rho tolerance, so the partition
        ! lag never limits the rho(T,ne) round-trip fidelity; NaN-safe: cells
        ! the outer ne iteration has driven to garbage must not keep the loop
        ! from exiting
        converged = .not. any(abs(hpart - HI%part) >= 1d-5)
        HI%part = hpart
        if (converged) exit
        ! re-solve hydrogen, tracking its electron contribution
        ! (H- has j=0 and contributes -f); xnecno excludes H
        call h_electrons(xnen)
        fe = fe - xnen
        call self%elements(1)%compute_ion_equilibrium(self%e_kT, self%lsaha, self%dxi)
        call h_electrons(xnen)
        fe = fe + xnen
        if (self%nmolec > 0) then
          call molecular_equilibrium(self%molecules, mions, self%tmolim, self%temp, &
                                     self%teta, self%e_kT, self%pe, xnecno, fe, &
                                     ysave=self%moleq_seed, same_state=self%same_temp)
        endif
        toth = ne / fe
        rho = toth * self%grph
      end do
      if (allocated(self%part_hocc)) then
        if (.not. all(shape(self%part_hocc) == shape(self%temp))) &
          deallocate(self%part_hocc)
      endif
      if (.not. allocated(self%part_hocc)) allocate(self%part_hocc, mold=self%temp)
      self%part_hocc = HI%part
    endif

    ! the fe the molecules balanced at, kept for the abundance-set re-solve
    if (allocated(self%moleq_fe)) then
      if (.not. all(shape(self%moleq_fe) == shape(fe))) deallocate(self%moleq_fe)
    endif
    if (.not. allocated(self%moleq_fe)) allocate(self%moleq_fe, mold=fe)
    self%moleq_fe = fe

    do z=1, size(rho,3)
      do y=1, size(rho,2)
        do x=1, size(rho,1)
          if (rho(x,y,z) /= rho(x,y,z)) then
            write(msg,'(a,3(1x,i0))') 'NaN density in thermal_equilibrium at x,y,z =', x, y, z
            call console_warn (ch_hop, msg)
          endif
        end do
      end do
    end do

    call trace%end(itimer)

  contains
    ! hydrogen's electron contribution per H nucleus (H- has j=0: -f)
    subroutine h_electrons(s)
      real(dp), dimension(:,:,:), intent(out) :: s
      integer :: jj
      s = 0
      associate(H_el => self%elements(1))
        do jj=lbound(H_el%ions,1), ubound(H_el%ions,1)
          s = s + H_el%ions(jj)%f * (jj-1)
        end do
      end associate
    end subroutine h_electrons
  END SUBROUTINE thermal_equilibrium

  !=========================================================================
  SUBROUTINE set_aux_data(self, temp, ne)
    class(composition_type), target :: self
    real(dp), intent(in),  dimension(:,:,:) :: temp, ne
    real(dp), parameter :: fo_const = 2*PI*(4./15.)**two_third*EC ! 1.25e-9
    real(dp), parameter :: dxi_const = (EC**3 * sqrt(8*PI)) / (temp_to_theta*BK*EE) ! 4.98e-4
    !---------------------------------------------------------------------------
    if (all(shape(self%temp) == shape(temp))) then
      self%same_temp = all(self%temp == temp)
    else
      self%same_temp = .false.
    endif
    self%ne = ne
    self%temp = temp
    self%teta = temp_to_theta / temp
    self%e_kT = EK / temp
    self%pe = temp * BK * ne
    self%lsaha = log(2 * op_saha_const * temp*sqrt(temp)) - log(ne)
    self%fo = fo_const * ne**two_third ! Holtsmark normal field strength

    if (debye_shield) then
      self%dxi = dxi_const * self%teta*sqrt(self%pe)
      ! --- this only shows up in the Oslo legacy code ---!
      !self%dxi=min(self%dxi,0.747*0.5)
      ! --------------------------------------------------!
    else
      self%dxi = 0
    endif
  END SUBROUTINE set_aux_data

  !=========================================================================
  ! find a state variables consistent with given temp & rho
  !=========================================================================
  SUBROUTINE iterate_thermal_eq(self, temp, rho, ne, tol)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    real(dp), dimension(:,:,:), intent(in) :: temp, rho
    real(dp), dimension(:,:,:), intent(out) :: ne
    real(dp), intent(in) :: tol ! stop criterion
    logical, dimension(:,:,:), allocatable :: ne_converged
    real(dp), dimension(:,:,:), allocatable :: a, b, c, x, xold, xprev2, xprev3, fx, fa, fb, ro, pg
    real(dp), dimension(:,:,:), allocatable :: lrho
    integer :: i, nconverged, arr_size, m(3), ix, iy, iz
    real(dp) :: nh, tt, sc, nel
    integer, parameter :: nmax = 40 ! max number of iteration
    character(len=256) :: msg
    !---------------------------------------------------------------------------
    m = shape(temp)
    arr_size = size(temp)

    allocate(ne_converged(m(1),m(2),m(3)), source=.false.) ! cant use mold= on different kinds
    allocate(x, xold, xprev2, xprev3, a, b, c, fx, fa, fb, ro, pg, lrho, mold=temp)

    if (ne(1,1,1) == 0.0d0) then
      ! Charge-neutrality guess from independent Saha quadratics
      ! ne^2/(n_donor - ne) = S(T) for the three donor channels that take
      ! over from each other as T drops: hydrogen, the abundant chi~8 eV
      ! metals (Mg+Si+Fe), and the easily ionized chi~5 eV species
      ! (Na+Ca+K).  The former hydrogen-only ionization-fraction formula
      ! underflowed to the tiny() floor below ~2500 K, which put the
      ! log-space secant start several hundred e-folds from the root and
      ! made the first thermal_equilibrium evaluations run at absurd
      ! states; its 0.1 fraction cap also started hot cells a factor ~10
      ! low.  Effective abundances and ionization energies are lumped;
      ! summing the channels independently overcounts ne slightly, which
      ! is harmless in an initial guess.
      do iz=1, m(3)
        do iy=1, m(2)
          do ix=1, m(1)
            nh = rho(ix,iy,iz) / self%grph
            tt = temp(ix,iy,iz)
            sc = op_saha_const * tt*sqrt(tt)
            nel =       donor_ne(       nh, sc*exp(max(-1.578d5/tt, -600d0)))
            nel = nel + donor_ne(8.0d-5*nh, sc*exp(max(-0.920d5/tt, -600d0)))
            nel = nel + donor_ne(7.0d-6*nh, sc*exp(max(-0.580d5/tt, -600d0)))
            ne(ix,iy,iz) = min(max(nel, 1d-80*nh), 1.5d0*nh)
          end do
        end do
      end do
    endif

    lrho = log(rho)

    call self%thermal_equilibrium(temp, ne, ro)

    where (ro==ro .and. ro>0) ! we skip NaN and non-positive values (log() would turn the
                               ! latter into NaN one step later, uncaught by ro==ro alone)
      fa = log(ro) - lrho
      ne_converged = abs(fa) < tol
      a = log(ne)
      b = a - 0.69*fa
      ne = exp(b)
    end where

    call self%thermal_equilibrium(temp, ne, ro)

    where(ro==ro .and. ro>0) ! we skip NaN and non-positive values
      fb = log(ro) - lrho
      ne_converged = abs(fb) < tol
    end where

    nconverged = count(ne_converged)

    i = 0
    x = b
    xprev2 = -huge(1d0) ! sentinels: no cycle possible against these on the first checks
    xprev3 = -huge(1d0)
    !write(*,fmt="(5x,'it [',I2.2,']',1x,'n-converged:',1x,I8,1X,'/ total size: ',I8)") i, nconverged, arr_size

    do while (nconverged < arr_size .and. i <= nmax)
      xold = x
      where (.not. ne_converged .and. ro==ro .and. ro>0)
        c = fb-fa ! log(ro2) - log(ro1)
        c = (b-a)/c
        x = a - c*fa
        x = max(x, min(a,b)-2.d0)
        x = min(x, max(a,b)+2.d0)
      end where

      ! The secant step above, once clamped to +-2 around the bracket, can
      ! fall into an exact repeating cycle when f is nearly flat locally
      ! (fb-fa tiny): the clamp keeps landing back on a value from one of
      ! the last few iterations instead of narrowing in (observed both
      ! period-2 and period-3 cycles in practice). Detect either per cell
      ! and fall back to bisection, which is guaranteed to shrink the
      ! bracket regardless of which endpoint ends up updated below.
      where (.not. ne_converged .and. ro==ro .and. ro>0 .and. &
             (abs(x-xprev2) < tol .or. abs(x-xprev3) < tol))
        x = 0.5d0*(a+b)
      end where
      xprev3 = xprev2
      xprev2 = xold

      ne_converged = abs(x-xold) < tol
      nconverged = count(ne_converged)
      if (arr_size == nconverged) exit ! all good !

      where (.not. ne_converged .and. ro==ro)
        ne = exp(x)
      endwhere

      call self%thermal_equilibrium(temp, ne, ro)
      where (ro==ro .and. ro>0) ! we skip NaN and non-positive values
        fx = log(ro) - lrho
        where (abs(a-x) < abs(b-x))
          b = x
          fb = fx
        else where
          a = x
          fa = fx
        end where
      end where
      i = i+1
    end do

    if (nconverged /= arr_size) then
      write(msg,'(a,i0,a,i0,a,i0,a)') 'electron density: ', nconverged, ' of ', arr_size, &
        ' points converged after ', i, ' iterations'
      call console_warn (ch_hop, msg)
    endif

    where (.not. ne_converged)
      self%ne = ieee_value(self%ne, ieee_quiet_nan)
      self%pg = ieee_value(self%pg, ieee_quiet_nan)
      self%E  = ieee_value(self%E,  ieee_quiet_nan)
    end where

    deallocate(ne_converged)
    deallocate(x,xold,a,b,fx,fa,fb)

  contains
    !---------------------------------------------------------------------------
    ! Stable root of the Saha quadratic ne^2 = s*(a - ne) for one donor
    ! species of number density a and Saha function s, branch-selected so
    ! that neither 4*a/s nor the subtraction can overflow or cancel
    !---------------------------------------------------------------------------
    function donor_ne(a, s) result(nel)
      real(dp), intent(in) :: a, s
      real(dp) :: nel
      !-------------------------------------------------------------------------
      if (s <= 0d0 .or. a <= 0d0) then
        nel = 0d0
      else if (a > 1d8*s) then
        nel = sqrt(a*s)                      ! strongly neutral limit
      else
        nel = 2d0*a / (1d0 + sqrt(1d0 + 4d0*a/s))
      endif
    end function donor_ne
  END SUBROUTINE iterate_thermal_eq

  !===============================================================================
  ! Hummer-Mihalas charged-perturber shielding parameter; depends on (temp, ne),
  ! so it must be refreshed whenever the state changes
  !-------------------------------------------------------------------------------
  SUBROUTINE compute_shield(self)
    class(composition_type), target :: self
    !---------------------------------------------------------------------------
    if (allocated(self%shield)) then
      if (.not. all(shape(self%shield) == shape(self%temp))) deallocate(self%shield)
    endif
    if (.not. allocated(self%shield)) allocate(self%shield, mold=self%temp)
    self%shield = 0.08989 * self%ne**0.16667d0 / sqrt(self%temp) ! the shielding parameter
    self%shield = (1.d0 + self%shield)**3.15d0 / self%ne
  END SUBROUTINE compute_shield

  !===============================================================================
  ! Occupation-weighted H I partition function U_w = sum_i g_i w_i exp(-E_i/kT)
  ! (linear, not log).  Refreshes the (temp, ne)-dependent inputs itself, so it
  ! is safe to call at any point of the ne iteration; the Boltzmann factors
  ! lvl%f are reused when temp is unchanged (same_temp).
  !-------------------------------------------------------------------------------
  SUBROUTINE Hoccupation_partition(self, nh, nhe, part)
    class(composition_type), target :: self
    real(dp), intent(in), dimension(:,:,:) :: nh, nhe
    real(dp), intent(out), dimension(:,:,:) :: part
    class(ion_type), pointer :: HI
    class(level_type), pointer :: lvl
    logical :: fresh
    integer :: i, itimer=0, itimer2=0
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%compute_occupation_weights', itimer=itimer)
    call self%ion_by_id('H I', HI)
    if (.not. allocated(HI%levels)) allocate(HI%levels(HI%nk))
    call self%compute_shield()

    part = 0.
    do i=1, HI%nk
      lvl => HI%levels(i)
      if (allocated(lvl%f)) then
        if (.not. all(shape(lvl%f) == shape(self%temp))) deallocate(lvl%f)
      endif
      fresh = .not. allocated(lvl%f)
      if (fresh) then
        allocate(lvl%f, mold=self%temp)
        call os%mem_register(storage_size(lvl%f), product(HI%m), 'lvl%f')
        lvl%i = i
        lvl%xion = RYD * (1. - 1./i**2)
        lvl%g = 2 * i**2
      endif
      if (fresh .or. .not. self%same_temp) then
        lvl%f = lvl%g * exp(-self%e_kT*lvl%xion)
      endif
      if (allocated(lvl%w)) then
        if (.not. all(shape(lvl%w) == shape(self%temp))) deallocate(lvl%w)
      endif
      if (.not. allocated(lvl%w)) then
        allocate(lvl%w, mold=self%temp)
        call os%mem_register(storage_size(lvl%w), product(HI%m), 'lvl%w')
      endif
      call trace%begin ('composition_type%compute_occupation_weights%wcalc', itimer=itimer2)
      call wcalc(nh, nhe, self%ne, self%shield, dble(i), lvl%w)
      call trace%end (itimer2)
      part = part + lvl%w * lvl%f
    end do
    call trace%end (itimer)
  END SUBROUTINE Hoccupation_partition

  !===============================================================================
  SUBROUTINE init_H_levels(self)
    class(composition_type), target :: self
    class(ion_type), pointer :: HI, HeI
    class(level_type), pointer :: lvl
    real(dp), dimension(:,:,:), allocatable :: h1, part
    integer, dimension(:,:,:), allocatable :: nniv
    integer :: i, itimer=0, itimer2=0
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%init_H_levels', itimer=itimer)
    call self%ion_by_id('H I', HI)
    call self%ion_by_id('HE I', HeI)

    call self%compute_shield()

    !HI%nk = 15
    ! the levels of Hoccupation_partition (f, w) are freed with HI%levels, so
    ! they leave the memory registry here as well
    if (allocated(HI%levels)) then
      do i=1, size(HI%levels)
        lvl => HI%levels(i)
        if (allocated(lvl%f)) call os%mem_register(storage_size(lvl%f), -size(lvl%f), 'lvl%f')
        if (allocated(lvl%w)) call os%mem_register(storage_size(lvl%w), -size(lvl%w), 'lvl%w')
        if (allocated(lvl%n)) call os%mem_register(storage_size(lvl%n), -size(lvl%n), 'lvl%n')
      end do
      deallocate(HI%levels)
    endif
    allocate(HI%levels(HI%nk))
    allocate(part, mold=self%temp)

    do i=1, HI%nk
      lvl => HI%levels(i)
      allocate(lvl%n, mold=self%temp)
      call os%mem_register(storage_size(lvl%n), product(HI%m), 'lvl%n')
      lvl%i = i
      lvl%xion = RYD * (1. - 1./i**2)
      lvl%g = 2 * i**2
      lvl%n = lvl%g * exp(-self%e_kT*lvl%xion)
    end do

    part = exp(-HI%part)
    if (H_occupation_probability) then
      part = 0.
      do i=1, HI%nk
        lvl => HI%levels(i)
        allocate(lvl%w, mold=self%temp)
        call os%mem_register(storage_size(lvl%w), product(HI%m), 'lvl%w')
        call trace%begin ('composition_type%init_H_levels%wcalc', itimer=itimer2)
        call wcalc(HI%n, HeI%n, self%ne, self%shield, dble(i), lvl%w)
        call trace%end (itimer2)
        part = part + lvl%w * lvl%n
        ! we chose not to apply the occupation probability to the hydrogen levels here
        ! as we avoid divisions in the opacity calculations later this way
      end do
      part = 1.d0 / part
    elseif (self%dxi(1,1,1) > 0) then
      allocate(h1,   source = sqrt(RYD/self%dxi) )
      allocate(nniv, source = int(h1)            )
      do i=1, HI%nk
        lvl => HI%levels(i)
        where (nniv < i)
          lvl%n = 0
        elsewhere (nniv == i)
          lvl%n = lvl%n * (h1 - nniv)
        end where
      end do
    endif
    do i=1, HI%nk
      lvl => HI%levels(i)
      lvl%n = HI%n * lvl%n * part
    end do
    call trace%end (itimer)
  END SUBROUTINE init_H_levels

  !===============================================================================
  SUBROUTINE compute_damping_aux(self)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion, HeI, HI
    class(molec_type), pointer :: H2, H2p
    integer :: i
    !---------------------------------------------------------------------------
    call self%molecule_by_id('H2',  H2)
    call self%molecule_by_id('H2+', H2p)
    call self%ion_by_id('H I', HI)
    call self%ion_by_id('HE I', HeI)

    self%vbar2h = HI%compute_vbar2(real(self%temp))
    self%t3nhe = (self%temp/10000.)**0.3 * HeI%n
    self%stark_pe = self%ne * self%temp**one_sixth
    self%tothi = HI%n ! total neutral hydrogen density

    if (associated(H2) .and. associated(H2p)) then
      self%p6 = LOG10((HI%n + 0.41336*HeI%n + 0.85*(H2%n + H2p%n)) * BK * self%temp)
    else
      self%p6 = LOG10((HI%n + 0.41336*HeI%n) * BK * self%temp)
    endif

    if (allocated(self%peachmet%n)) deallocate(self%peachmet%n)
    allocate(self%peachmet%n, source=self%zero)
    do i=1, self%nel
      element => self%elements(i)
      if (element%nion >= 2) then
        ion => element%ions(2)
        select case (trim(element%id))
         case ('H', 'He', 'C', 'Mg', 'Si')
          ! do not add
         case default
          self%peachmet%n = self%peachmet%n + ion%n
        end select
      end if
    end do
  END SUBROUTINE compute_damping_aux

  !=============================================================================
  ! Boltzmann factor beta = exp(-h nu / kT), clamped to the smallest positive
  ! real(4) so that ehat = eta/beta can never divide by zero.  beta is the one
  ! quantity that converts between physical extinction and the working
  ! (xhat, ehat) pair, so every routine below takes it as an argument rather
  ! than recomputing the exponential.
  !=============================================================================
  SUBROUTINE boltzmann_factor(self, sp, beta, temp)
    class(composition_type), target :: self
    class(spectrum_type), intent(in) :: sp
    real(4), dimension(:,:,:,:), contiguous, intent(out) :: beta
    real(dp), dimension(:,:,:), intent(in), optional, target :: temp
    real(dp), dimension(:,:,:), pointer :: t
    integer :: x, y, z, f, nf, m(3)
    real(8) :: hk_T
    !---------------------------------------------------------------------------
    if (present(temp)) then
      t => temp
    else
      t => self%temp
    endif
    m = shape(t)
    nf = size(beta, 1)
    if (size(beta, 2) /= m(1)) &
      error stop 'boltzmann_factor: beta and temp span different cells'
    ! h/(k T) is frequency-independent: hoist the division out of the frequency
    ! loop so the inner loop is one multiply plus the exponential.
    do z=1, m(3)
      do y=1, m(2)
        do x=1, m(1)
          hk_T = h_k / t(x,y,z)
          do f=1, nf
            beta(f,x,y,z) = max(real(exp(-sp%nu(f) * hk_T), 4), tiny(1.0))
          end do
        end do
      end do
    end do
  END SUBROUTINE boltzmann_factor

  !=============================================================================
  ! OPACITY PIPELINE
  !
  ! Three primitives -- continuum_opacities, add_absmet_opacities and
  ! add_linelist_opacities -- accumulate onto the same two buffers, and
  ! `opacity` below is a thin wrapper that runs all three for the callers that
  ! want one call.  The anisotropic rt_patch foreground drives the same
  ! primitives directly, interleaved with its bound-free/bound-bound kernels, so
  ! both paths share one implementation.
  !
  ! Two accumulators, in two different unit systems, deliberately:
  !
  !   xhat, ehat  working units.  xhat = chi + eta, ehat = eta/beta, so that a
  !               line adds R*b_i to xhat and R*b_j to ehat -- two positive
  !               population-weighted sums, with the subtractive
  !               chi = xhat - ehat*beta formed once at the end instead of per
  !               line.  In LTE the two coincide, xhat = ehat = chi/(1-beta),
  !               and that is also the pre-stimulated-emission accumulation.
  !
  !   sigma       coherent scattering, in physical 1/cm.  It NEVER enters the
  !               (xhat, ehat) transform: a pure scatterer has no thermal source
  !               of its own, and its emission depends on Jnu, which the
  !               transfer solution supplies later.  Keeping it a separate,
  !               purely additive accumulator is what makes
  !                   sigma <= chi_total = chi_true + sigma
  !               structural rather than something to clamp.
  !=============================================================================
  ! Continuum true absorption in the LTE working form, plus the continuum
  ! scattering opacity:
  !
  !   xhat = ehat = chi_true_cont / (1 - beta)   [working units]
  !   sigma       = electron + Rayleigh          [1/cm]
  !
  ! The working form is the *cheaper* output, not the more expensive one: the
  ! opacity tables do not carry stimulated emission, so their raw sum already is
  ! xhat.  Only transitions flagged istim arrive stim-corrected and have to be
  ! lifted by 1/(1-beta); they go through a separate buffer so the dominant term
  ! never round-trips through a small 1-beta.
  !=============================================================================
  SUBROUTINE continuum_opacities(self, sp, xhat, sigma, beta, sorted, i0, i1)
    class(composition_type), target :: self
    class(spectrum_type), intent(in) :: sp
    class(bftrans_type), pointer :: trans
    class(ion_type), pointer :: HeI
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: xhat
    real(4), dimension(:,:,:,:), contiguous, intent(inout), optional :: sigma
    real(4), dimension(:,:,:,:), contiguous, intent(in), target, optional :: beta
    logical, optional :: sorted
    integer, intent(in), optional :: i0, i1        ! optional cell range on dim 1 (threading)
    real(4), dimension(:,:,:,:), pointer, contiguous :: bet
    real(4), dimension(:,:,:,:), allocatable, target :: bet_mem
    real(4), dimension(:,:,:,:), allocatable :: tmp
    logical :: any_istim
    integer :: i, x, y, z, f, nf, nx, ny, nz, lo, hi, nn, boff
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%continuum_opacities', itimer=itimer)
    call self%ion_by_id('HE I', HeI)
    nx = size(self%temp, 1)
    ny = size(self%temp, 2)
    nz = size(self%temp, 3)
    nf = size(sp%nu, 1)
    lo = 1; hi = nx
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    nn = hi - lo + 1
    if (size(xhat, 2) /= nx) error stop 'continuum_opacities: xhat must be whole-tile'
    if (present(sigma)) then
      if (size(sigma, 2) /= nx) error stop 'continuum_opacities: sigma must be whole-tile'
    endif

    ! the tables below carry no stimulated emission, so their sum IS xhat = ehat.
    xhat(:,lo:hi,:,:) = 0
    call self%add_hrest(xhat, sp%lam, lo, hi)
    any_istim = .false.
    do i=1, size(self%trans)
      trans => self%trans(i)
      if (trans%istim) then
        any_istim = .true.
      else
        call trans%add_absorption(xhat, sp%lam, self%temp, self%e_kT, self%teta, &
          self%pe, self%dxi, ne=self%ne, nhe=HeI%n, shield=self%shield, sorted=sorted, i0=lo, i1=hi)
      endif
    end do

    if (any_istim) then
      boff = lo - 1
      if (present(beta)) then
        bet => beta
        if (size(bet, 2) == nx) then
          boff = 0
        else if (size(bet, 2) /= nn) then
          error stop 'continuum_opacities: beta is neither whole-tile nor chunk-local'
        endif
      else
        allocate(bet_mem(nf, nn, ny, nz))
        call self%boltzmann_factor(sp, bet_mem, self%temp(lo:hi,:,:))
        bet => bet_mem
      endif
      allocate(tmp(nf, nn, ny, nz))
      tmp = 0
      do i=1, size(self%trans)
        trans => self%trans(i)
        if (trans%istim) then
          call trans%add_absorption(tmp, sp%lam, self%temp, self%e_kT, self%teta, &
            self%pe, self%dxi, ne=self%ne, nhe=HeI%n, shield=self%shield, &
            sorted=sorted, i0=lo, i1=hi)
        endif
      end do
      do z=1, nz
        do y=1, ny
          do x=lo, hi
            do f=1, nf
              xhat(f,x,y,z) = xhat(f,x,y,z) &
                + tmp(f,x-lo+1,y,z) / (1 - bet(f,x-boff,y,z))
            end do
          end do
        end do
      end do
      deallocate(tmp)
    endif

    if (present(sigma)) then
      sigma(:,lo:hi,:,:) = 0
      call self%electron_scattering(sigma, lo, hi)
      call self%rayleigh_scattering(sp%lam, sigma, lo, hi)
    endif

    call trace%end(itimer)
  END SUBROUTINE continuum_opacities

  !=============================================================================
  ! Add the sampled metal-line (absmet) opacity to the same two accumulators.
  ! osmet%add_line_opacities returns *physical* extinction, so its thermal part
  ! is lifted into the working form by 1/(1-beta):
  !
  !   xhat  += eps*chi_met/(1 - beta)      sigma += (1-eps)*chi_met
  !
  ! split governs the eps split, and follows line_scattering just as it does for
  ! the line list -- absmet opacities *are* line opacities.  With split=.false.
  ! metal lines are pure absorbers (eps = 1), which is what line_scattering=F
  ! means everywhere else in the code.
  !=============================================================================
  SUBROUTINE add_absmet_opacities(self, sp, xhat, beta, sigma, split)
    class(composition_type), target :: self
    class(spectrum_type), intent(in) :: sp
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: xhat
    real(4), dimension(:,:,:,:), contiguous, intent(in) :: beta
    real(4), dimension(:,:,:,:), contiguous, intent(inout), optional :: sigma
    logical, intent(in), optional :: split
    real(4), dimension(:,:,:,:), allocatable :: chi_met, epsln
    logical :: do_split
    integer :: i, j, k, f, nf, m(3)
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    if (.not. osmet%initialised) return
    call trace%begin ('composition_type%add_absmet_opacities', itimer=itimer)
    m = shape(self%temp)
    nf = size(xhat, 1)
    do_split = .false.
    if (present(split)) do_split = split .and. present(sigma)

    allocate(chi_met, mold=xhat)
    chi_met = 0
    call osmet%add_line_opacities(sp%lam, self%p6, self%temp, self%rho, chi_met)

    if (do_split) then
      allocate(epsln, mold=xhat)
      call osmet%compute_epsilon(sp%nu, self%temp, self%ne, epsln)
      do k=1, m(3)
        do j=1, m(2)
          do i=1, m(1)
            do f=1, nf
              xhat (f,i,j,k) = xhat (f,i,j,k) &
                + epsln(f,i,j,k)*chi_met(f,i,j,k) / (1 - beta(f,i,j,k))
              sigma(f,i,j,k) = sigma(f,i,j,k) &
                + (1 - epsln(f,i,j,k))*chi_met(f,i,j,k)
            end do
          end do
        end do
      end do
      deallocate(epsln)
    else
      do k=1, m(3)
        do j=1, m(2)
          do i=1, m(1)
            do f=1, nf
              xhat(f,i,j,k) = xhat(f,i,j,k) + chi_met(f,i,j,k) / (1 - beta(f,i,j,k))
            end do
          end do
        end do
      end do
    endif
    deallocate(chi_met)
    call trace%end(itimer)
  END SUBROUTINE add_absmet_opacities

  !=============================================================================
  ! One call over the whole pipeline, for the callers that want the finished
  ! numbers: the standalone/Python API, chi500, and the opacity-table builders.
  !
  !   chi   total extinction chi_true + scat   [1/cm]   (disjoint=.false., default)
  !         true absorption alone              [1/cm]   (disjoint=.true.)
  !   scat  electron + Rayleigh, plus -- when line_scat -- the (1-eps) fraction
  !         of the absmet and of the LTE line-list opacity.  Lines carrying
  !         departure coefficients are left out: their scattering lives in the
  !         NLTE rate equations.
  !
  ! The eps split is a *partition*: scat is a share of an extinction that is
  ! already whole, never an addition to it.  So chi at disjoint=.false. is
  ! invariant under line_scat -- flipping the split moves opacity between the
  ! two outputs and must not change their sum.
  !
  ! The default keeps the published meaning of hop.opacity() and of the
  ! opa/scat tables: chi is the total, scat the scattering part *of* it.
  ! disjoint=.true. is for the sampled background table, where tabulating the
  ! two components separately and summing after interpolation is what makes the
  ! interpolated albedo scat/chi_total <= 1 by construction.
  !
  ! beta is optional; supply it (the rt path has it stored) to skip nf*ncell
  ! exponentials.
  !=============================================================================
  SUBROUTINE opacity (self, sp, chi, scat, include_absmet, include_llist, dlam, &
                      sorted, line_scat, vproj, departures, beta, disjoint)
    class(composition_type), target :: self
    class(spectrum_type), intent(in) :: sp
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
    real(4), dimension(:,:,:,:), contiguous, intent(inout), optional, target :: scat
    logical, optional :: include_absmet, include_llist, sorted, line_scat, disjoint
    real(8), intent(in), optional :: dlam
    real(4), dimension(:,:,:), intent(in), contiguous, target, optional :: vproj
    class(departure_pointer_t), intent(in), dimension(:), target, optional :: departures
    real(4), dimension(:,:,:,:), contiguous, intent(in), target, optional :: beta
    real(4), dimension(:,:,:,:), pointer, contiguous :: bet, sig
    real(4), dimension(:,:,:,:), allocatable, target :: bet_mem, sig_mem
    real(4), dimension(:,:,:,:), allocatable :: ehat, chi_scat
    real(4), dimension(:,:,:,:,:), allocatable :: one_m_eps
    logical :: do_absmet, do_llist, do_split, do_disjoint
    real(8) :: dl
    integer :: i, j, k, f, nf, m(3)
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%opacity', itimer=itimer)
    m = shape(self%temp)
    nf = size(chi, 1)
    dl = 100d0
    if (present(dlam)) dl = dlam
    do_absmet = .false.
    if (present(include_absmet)) do_absmet = include_absmet .and. osmet%initialised
    do_llist = self%ll_initialised
    if (present(include_llist)) do_llist = include_llist .and. self%ll_initialised
    do_split = line_scattering
    if (present(line_scat)) do_split = line_scat
    do_split = do_split .and. present(scat)
    do_disjoint = .false.
    if (present(disjoint)) do_disjoint = disjoint .and. present(scat)

    if (present(beta)) then
      bet => beta
    else
      allocate(bet_mem, mold=chi)
      call self%boltzmann_factor(sp, bet_mem)
      bet => bet_mem
    endif
    if (present(scat)) then
      sig => scat
    else
      allocate(sig_mem, mold=chi)
      sig => sig_mem
    endif

    ! chi is the working xhat from here until the single back-conversion below
    call self%continuum_opacities(sp, chi, sigma=sig, beta=bet, sorted=sorted)
    if (do_absmet) call self%add_absmet_opacities(sp, chi, bet, sigma=sig, split=do_split)

    if (do_llist) then
      allocate(ehat, mold=chi)
      ehat = chi                       ! in LTE the continuum's two forms coincide
      if (do_split) then
        allocate(chi_scat, mold=chi)
        allocate(one_m_eps(nf, m(1), m(2), m(3), 2))
        chi_scat = 0
        call compute_epsilon_stage(1, sp%nu, self%temp, self%ne, one_m_eps(:,:,:,:,1))
        call compute_epsilon_stage(2, sp%nu, self%temp, self%ne, one_m_eps(:,:,:,:,2))
        call self%add_linelist_opacities(sp, chi, ehat, dl, vproj, departures=departures, &
          one_m_eps=one_m_eps, chi_scat=chi_scat, beta=bet)
        do k=1, m(3)
          do j=1, m(2)
            do i=1, m(1)
              do f=1, nf
                sig(f,i,j,k) = sig(f,i,j,k) + chi_scat(f,i,j,k)
              end do
            end do
          end do
        end do
      else
        call self%add_linelist_opacities(sp, chi, ehat, dl, vproj, departures=departures)
      endif
      ! (xhat, ehat) -> physical true absorption, once.  The line kernels add the
      ! whole LTE line to xhat and only *record* its (1-eps) share in chi_scat, so
      ! with the split on that share is still sitting in chi while it is also being
      ! handed to sig above.  Take it out here: the eps split is a partition of the
      ! line opacity, not an addition to it, exactly as add_absmet_opacities does
      ! for the sampled metal lines.  chi_scat is already physical (1/cm), the
      ! kernels' working form having been undone with beta on the way out.
      if (do_split) then
        do k=1, m(3)
          do j=1, m(2)
            do i=1, m(1)
              do f=1, nf
                chi(f,i,j,k) = chi(f,i,j,k) - ehat(f,i,j,k)*bet(f,i,j,k) &
                                            - chi_scat(f,i,j,k)
              end do
            end do
          end do
        end do
      else
        do k=1, m(3)
          do j=1, m(2)
            do i=1, m(1)
              do f=1, nf
                chi(f,i,j,k) = chi(f,i,j,k) - ehat(f,i,j,k)*bet(f,i,j,k)
              end do
            end do
          end do
        end do
      endif
    else
      do k=1, m(3)
        do j=1, m(2)
          do i=1, m(1)
            do f=1, nf
              chi(f,i,j,k) = chi(f,i,j,k) * (1 - bet(f,i,j,k))
            end do
          end do
        end do
      end do
    endif

    if (.not. do_disjoint) then
      do k=1, m(3)
        do j=1, m(2)
          do i=1, m(1)
            do f=1, nf
              chi(f,i,j,k) = chi(f,i,j,k) + sig(f,i,j,k)
            end do
          end do
        end do
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE opacity

  !===============================================================================
  !> Precompute the population changes for molecular abundance sets.
  !>
  !> A set's offset elements are re-solved through molecular_equilibrium at the
  !> model's own (T, ne), warm-started from the base solution, and the resulting
  !> populations are stored for the solve elements and the molecules. Elements
  !> outside the molecular network are represented only by set_da and are scaled
  !> where the opacity is formed: at fixed (T, ne) an offset leaves their
  !> ionisation fractions alone, so the scaling is exact.
  !>
  !> What this does not carry is the continuum. The background is computed once
  !> for the first set and shared (see rt_patch_mod's i_set==1 guard), while a
  !> real run at that abundance re-forms it: measured on the CH G band of a cool
  !> giant, a 0.2 dex C offset moves the continuum by 2.1%, against a line signal
  !> of 1.8%, and what survives continuum normalisation is 11.5% of that signal.
  !> The error is a wrong derivative -- it does not shrink with the offset -- so
  !> a ladder on a molecule-forming element in a cool star carries it at any
  !> rung. Restoring it means a per-set background, which the shared-continuum
  !> design deliberately trades away.
  !===============================================================================
  SUBROUTINE prepare_abundance_sets(self, da, present_element)
    class(composition_type), target :: self
    type(ion_pointer_t), allocatable :: mions(:)
    real(dp), allocatable :: fe(:,:,:), fe0(:,:,:), xnecno(:,:,:), seed(:,:,:,:)
    real(8), intent(in) :: da(:,:)
    logical, intent(in) :: present_element(:,:)
    real(dp) :: base_abund(92), frac_change
    integer :: i, j, k, s, nions, nstore, ipop
    logical, allocatable :: molecular_set(:)
    !-----------------------------------------------------------------------------
    if (size(da,1) /= 92 .or. any(shape(da) /= shape(present_element))) &
      error stop 'prepare_abundance_sets: invalid set layout shape'
    allocate(self%set_da, source=da)
    allocate(self%set_element, source=present_element)
    allocate(molecular_set(size(da,2)), source=.false.)
    base_abund = 0d0
    do i=1,self%nel
      if (self%elements(i)%np > 0) base_abund(self%elements(i)%np) = self%elements(i)%abund
      if (self%elements(i)%ions(1)%nmolec > 0 .and. self%elements(i)%np > 0) then
        do s=1,size(da,2)
          molecular_set(s) = molecular_set(s) .or. present_element(self%elements(i)%np,s)
        enddo
      endif
    enddo
    ! Table size: a molecular set stores populations for every solve element's
    ! ions and every molecule, and any set stores one scale entry per ion of each
    ! offset element that forms no molecules.
    nions = count([(self%elements(i)%ions(1)%nmolec > 0, i=1,self%nel)])
    nstore = self%nmolec
    do i=1,self%nel
      if (self%elements(i)%ions(1)%nmolec <= 0) cycle
      nstore = nstore + size(self%elements(i)%ions)
    enddo
    nstore = nstore*count(molecular_set)
    do s=1,size(da,2)
      do i=1,self%nel
        if (self%elements(i)%np <= 0) cycle
        if (self%elements(i)%ions(1)%nmolec > 0) cycle
        if (present_element(self%elements(i)%np,s)) nstore = nstore + size(self%elements(i)%ions)
      enddo
    enddo
    if (nstore == 0) return
    allocate(self%set_pop(nstore))
    ipop = 0

    ! the scale entries: at fixed (T, ne) an element outside the molecular
    ! network keeps its ionisation fractions, so every stage follows its
    ! abundance by the same factor
    do s=1,size(da,2)
      do i=1,self%nel
        if (self%elements(i)%np <= 0) cycle
        if (self%elements(i)%ions(1)%nmolec > 0) cycle
        if (.not. present_element(self%elements(i)%np,s)) cycle
        do j=lbound(self%elements(i)%ions,1),ubound(self%elements(i)%ions,1)
          ipop = ipop + 1
          self%set_pop(ipop)%i_set = s
          self%set_pop(ipop)%id = self%elements(i)%ions(j)%id
          self%set_pop(ipop)%scale = 10d0**da(self%elements(i)%np,s)
        enddo
      enddo
    enddo
    if (.not. any(molecular_set)) return

    allocate(mions(nions))
    k = 0
    do i=1,self%nel
      if (self%elements(i)%ions(1)%nmolec <= 0) cycle
      k = k + 1
      mions(k)%p => self%elements(i)%ions(1)
    enddo
    allocate(fe, fe0, xnecno, mold=self%temp)
    ! the base solve's own converged fe: re-entering the Newton with it and an
    ! unchanged composition is a fixed point, which is what makes a da = 0 set
    ! reproduce the base state exactly.
    if (.not. allocated(self%moleq_fe)) &
      error stop 'prepare_abundance_sets: no converged molecular equilibrium to re-enter'
    fe0 = self%moleq_fe

    ! Electrons from the elements that form no molecules. At fixed (T, ne) an
    ! offset elsewhere leaves their ionisation fractions alone, so this is the
    ! same for every set and is built once.
    xnecno = 0d0
    do i=1,self%nel
      if (any(self%elements(i)%ions%nmolec /= 0)) cycle
      do j=lbound(self%elements(i)%ions,1),ubound(self%elements(i)%ions,1)
        xnecno = xnecno + self%elements(i)%ions(j)%f * (j-1)
      enddo
    enddo

    do s=1,size(da,2)
      if (.not. molecular_set(s)) cycle
      do i=1,self%nel
        if (self%elements(i)%np > 0 .and. present_element(self%elements(i)%np,s)) &
          self%elements(i)%abund = base_abund(self%elements(i)%np) * 10d0**da(self%elements(i)%np,s)
      enddo

      fe = fe0
      if (allocated(self%moleq_seed)) allocate(seed, source=self%moleq_seed)
      call molecular_equilibrium(self%molecules, mions, self%tmolim, self%temp, &
                                 self%teta, self%e_kT, self%pe, xnecno, fe, &
                                 ysave=seed, same_state=self%same_temp)

      do i=1,self%nel
        if (self%elements(i)%ions(1)%nmolec <= 0) cycle
        do j=lbound(self%elements(i)%ions,1),ubound(self%elements(i)%ions,1)
          ipop = ipop + 1
          self%set_pop(ipop)%i_set = s
          self%set_pop(ipop)%id = self%elements(i)%ions(j)%id
          allocate(self%set_pop(ipop)%n, source=self%elements(i)%ions(j)%f*self%toth)
        enddo
      enddo
      do i=1,self%nmolec
        ipop = ipop + 1
        self%set_pop(ipop)%i_set = s
        self%set_pop(ipop)%id = self%molecules(i)%id
        allocate(self%set_pop(ipop)%n, source=self%molecules(i)%f*self%toth)
      enddo

      frac_change = maxval(abs(fe-fe0) / max(abs(fe0), tiny(1d0)))
      if (frac_change > 1d-3 .and. console_on(ch_hop, 2)) &
        write(console_unit,'(a,i0,a,es10.3)') 'abundance set ', s, &
          ' changes the electron fraction at fixed ne by ', frac_change

      ! Back to the base state for the next set and for the RT. The number
      ! densities n were formed before this routine and the solver never touches
      ! them, so the fractions they came from are recovered exactly.
      do i=1,self%nel
        if (self%elements(i)%np > 0) self%elements(i)%abund = base_abund(self%elements(i)%np)
        if (self%elements(i)%ions(1)%nmolec <= 0) cycle
        do j=lbound(self%elements(i)%ions,1),ubound(self%elements(i)%ions,1)
          self%elements(i)%ions(j)%f = self%elements(i)%ions(j)%n / self%toth
        enddo
      enddo
      do i=1,self%nmolec
        self%molecules(i)%f = self%molecules(i)%n / self%toth
      enddo
      if (allocated(seed)) deallocate(seed)
    enddo
  END SUBROUTINE prepare_abundance_sets

  !===============================================================================
  SUBROUTINE set_population(self, id, i_set, population, scale)
    class(composition_type), target :: self
    character(len=*), intent(in) :: id
    integer, intent(in) :: i_set
    real(dp), dimension(:,:,:), pointer, intent(out) :: population
    real(dp), intent(out) :: scale
    integer :: i
    !-----------------------------------------------------------------------------
    ! A table match, and nothing else: this runs once per species per set inside
    ! the opacity loop, so the species -> element resolution (and the abundance
    ! it implies) is done when the sets are prepared, not here.
    population => null()
    scale = 1d0
    if (.not. allocated(self%set_pop)) return
    do i=1,size(self%set_pop)
      if (self%set_pop(i)%i_set == i_set .and. &
          trim(self%set_pop(i)%id) == trim(id)) then
        population => self%set_pop(i)%n
        scale = self%set_pop(i)%scale
        return
      endif
    enddo
  END SUBROUTINE set_population

  !===============================================================================
  SUBROUTINE population_ratio(self, id, i_set, ratio)
    class(composition_type), target :: self
    character(len=*), intent(in) :: id
    integer, intent(in) :: i_set
    real(dp), intent(out), dimension(:,:,:) :: ratio
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    real(dp), dimension(:,:,:), pointer :: pop
    real(dp), allocatable :: base(:,:,:), changed(:,:,:)
    real(dp) :: scale
    integer :: j
    !-----------------------------------------------------------------------------
    ratio = 1d0
    call self%element_by_id(id, element, verbose=.false.)
    if (associated(element)) then
      allocate(base, changed, mold=ratio)
      base = 0d0; changed = 0d0
      do j=lbound(element%ions,1),ubound(element%ions,1)
        base = base + element%ions(j)%n
        call self%set_population(element%ions(j)%id, i_set, pop, scale)
        if (associated(pop)) then
          changed = changed + pop
        else
          changed = changed + scale*element%ions(j)%n
        endif
      enddo
      where (base > 0d0) ratio = changed/base
      return
    endif
    call self%molecule_by_id(id, mol, verbose=.false.)
    if (associated(mol)) then
      call self%set_population(mol%id, i_set, pop, scale)
      if (associated(pop)) then
        where (mol%n > 0d0) ratio = pop/mol%n
      else
        ratio = scale
      endif
    endif
  END SUBROUTINE population_ratio

  !===============================================================================
  SUBROUTINE get_totn(self, id, abund, totn, el, mol, i_set)
    class(composition_type), target :: self
    class(element_type), intent(out), optional, pointer :: el
    class(molec_type), intent(out), optional, pointer :: mol
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    character(len=*), intent(in) :: id
    real(dp), intent(in) :: abund
    real(dp), dimension(:,:,:) :: totn
    integer, intent(in), optional :: i_set
    real(dp), allocatable :: ratio(:,:,:)
    integer :: j
    !---------------------------------------------------------------------------
    call trace%begin ('composition_type%get_totn')
    element => null()
    molec   => null()

    call self%element_by_id(id, element)
    if (associated(element)) then
      totn = 0
      do j=1, element%nion
        totn = totn + element%ions(j)%n
      end do
    else
      call self%molecule_by_id(id, molec)
      if (associated(molec)) then
        totn = molec%n
      else
        totn = 10.0**(abund-12.0) * self%toth
      endif
    endif

    if (present(i_set)) then
      allocate(ratio, mold=totn)
      call self%population_ratio(id, i_set, ratio)
      totn = totn * ratio
    endif

    if (present(el) .and. present(mol)) then
      el => element
      mol => molec
    endif
    call trace%end ()
  END SUBROUTINE get_totn

  !===============================================================================
  SUBROUTINE compute_energy_and_gas_pressure(self, eh, e, pg)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    class(ion_type), pointer :: ion, HI, Hm
    real(dp), intent(out), dimension(:,:,:) :: eh, pg
    real(dp), intent(inout), dimension(:,:,:) :: e
    real(dp), dimension(:,:,:), allocatable :: xerg, ejon, emol
    integer :: i, j
    !---------------------------------------------------------------------------
    allocate(xerg, ejon, emol, mold=pg)
    call self%ion_by_id('H I', HI)
    call self%ion_by_id('H -', Hm)

    ejon = 0
    pg = self%ne
    eh = - (HI%xion - self%dxi)*(HI%n + Hm%n) - (Hm%xion - 2*self%dxi)*Hm%n

    do i=1, self%nel
      element => self%elements(i)
      do j=1, element%nion
        ion => element%ions(j)
        pg = pg + ion%n
        if (i>1 .and. associated(ion%above)) then
          xerg = ion%xion - self%dxi*j
          ! energy of ionization (ejon). hydrogen is not included.
          ejon = ejon + ion%above%n * xerg
        endif
      end do
    end do

    emol = 0
    do i=1, self%nmolec
      mol => self%molecules(i)
      pg = pg + mol%n
      emol = emol + mol%calc_e(self%teta, self%temp, self%dxi)
    end do

    ! The gas pressure is Pg = NkT, where N is (Ne + Ni + Nmol)
    pg = pg * BK * self%temp

    !ehm = (Hm%xion + HI%xion) * Hm%f
    where (.not.ieee_is_nan(E))
      E = EE*(ejon + eh + emol)
      E = (1.5*pg + E) / self%rho
    end where
  END SUBROUTINE compute_energy_and_gas_pressure

  !===============================================================================
  SUBROUTINE molecule_by_id(self, id, molecule, verbose)
    class(composition_type), target :: self
    character(len=*), intent(in) :: id
    class(molec_type), pointer, intent(out) :: molecule
    integer :: idx
    logical, optional :: verbose
    !---------------------------------------------------------------------------
    idx = findloc(self%molecules%id, id, 1)
    if (idx == 0) then
      molecule => null()
      if (option(verbose)) then
        call console_warn (ch_hop, 'molecule '//trim(id)//' not found')
      endif
    else
      molecule => self%molecules(idx)
    endif
  END SUBROUTINE molecule_by_id

  !===============================================================================
  SUBROUTINE element_by_id(self, id, element, verbose)
    class(composition_type), target :: self
    character(len=*), intent(in) :: id
    class(element_type), pointer, intent(out) :: element
    integer :: idx
    logical, optional :: verbose
    !---------------------------------------------------------------------------
    idx = findloc(self%elements%id, capitalize(id), 1)

    element => null()
    if (idx == 0) then
      if (option(verbose)) then
        call console_warn (ch_hop, 'element '//trim(id)//' not found')
      endif
    else
      element => self%elements(idx)
    endif
  END SUBROUTINE element_by_id

  !===============================================================================
  SUBROUTINE ion_by_id(self, id, ion, verbose, element_out)
    class(composition_type), target :: self
    character(len=*), intent(in) :: id
    class(ion_type), pointer, intent(out) :: ion
    class(element_type), intent(out), pointer, optional :: element_out
    class(element_type), pointer :: element
    logical, optional :: verbose
    character(3) :: eid, iid
    character(6) :: trim_id
    integer :: idx, iostat
    !---------------------------------------------------------------------------
    read(id, *, iostat=iostat) eid, iid
    eid = capitalize(eid)
    iid = uppercase(iid)
    ion => null()
    element => null()
    if (iostat /= 0) then
      if (option(verbose)) then
        call console_warn (ch_hop, 'ion '//trim(id)//' not found')
      endif
    else
      trim_id = trim(eid)//' '//trim(iid)
      call self%element_by_id(eid, element)
      if (associated(element)) then
        idx = findloc(element%ions%id, trim_id, 1)
        if (lbound(element%ions,1) /= 1) then
          idx = idx + lbound(element%ions,1) - 1
        endif
        if (idx < lbound(element%ions,1) .or. idx>ubound(element%ions,1)) then
          call console_warn (ch_hop, 'ion '//trim(trim_id)//' not found')
        else
          ion => element%ions(idx)
        endif
      endif
      if (present(element_out)) then
        element_out => element
      endif
    endif
  END SUBROUTINE ion_by_id

  !===============================================================================
  SUBROUTINE particle_by_id(self, id, particle, verbose, np, ion, mol, is_molecule)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    character(len=*), intent(in) :: id
    class(particle_type), pointer, intent(out) :: particle
    logical, intent(in), optional :: verbose
    integer, intent(out), optional :: np, ion, mol
    logical, intent(out), optional :: is_molecule
    logical :: is_mol
    character(3) :: id1, id2
    character(6) :: trim_id
    integer :: idx, imol, iele
    !---------------------------------------------------------------------------
    element => null()
    particle => null()
    is_mol = .false.
    read(id, *) id1
    iele = findloc(self%elements%id, capitalize(id1), 1)
    imol = findloc(self%molecules%id, id1, 1)

    if (iele > 0 .and. .not. imol > 0) then
      element => self%elements(iele)
      is_mol = .false.
    elseif (imol > 0 .and. .not. iele > 0) then
      particle => self%molecules(imol)
      is_mol = .true.
    elseif (imol > 0 .and. iele > 0) then
      if (findloc(self%elements%id, id1, 1) > 0) then
        element => self%elements(iele)
        is_mol = .false.
      else
        particle => self%molecules(imol)
        is_mol = .true.
      endif
    endif

    if (associated(element)) then
      read(id, *) id1, id2
      id1 = capitalize(id1)
      id2 = uppercase(id2)
      trim_id = trim(id1)//' '//trim(id2)
      idx = findloc(element%ions%id, trim_id, 1)
      if (idx > 0) then
        idx = idx + (lbound(element%ions,1)-1)
        particle => element%ions(idx)
        if (present(np)) np = element%np
        if (present(ion)) ion = idx
      endif
    endif
    if (present(mol) .and. is_mol) mol = imol
    if (present(is_molecule)) is_molecule = is_mol
    if (option(verbose) .and. .not. associated(particle)) then
      call console_warn (ch_hop, 'particle '//trim(id)//' not recognised')
    endif
  END SUBROUTINE particle_by_id

  !=========================================================================
  SUBROUTINE compute_grams_per_hydrogen(self)
    class(composition_type), target :: self
    !---------------------------------------------------------------------
    ! Abundances should not be in log units when they enter here!
    ! XMY = GRAMS OF STELLAR MATTER / GRAMS OF HYDROGEN
    !---------------------------------------------------------------------
    self%suma = sum(self%elements%abund)
    self%xmy  = sum(self%elements%abund*self%elements%awgt) / self%elements(1)%awgt
    self%grph = XMH*self%xmy
  END SUBROUTINE compute_grams_per_hydrogen

  !=========================================================================
  SUBROUTINE overwrite_abundance(self, species, abnd, awgt, relative_abundance)
    class(composition_type), target :: self
    class(massive_type), pointer :: species
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    class(ion_type), pointer :: ion
    real(dp), intent(in) :: abnd, awgt
    logical, intent(in) :: relative_abundance
    real(dp) :: abund, a
    integer :: j
    !---------------------------------------------------------------------
    select type (species)
     class is (element_type)
      element => species
     class is (molec_type)
      mol => species
      do j=1, mol%nions
        ion => mol%ions(j)%p
        if (ion%id(1:2) /= 'H ') then
          call self%element_by_id(ion%id(1:2), element, verbose=.false.)
          exit
        endif
      end do
     class default
      call console_warn (ch_hop, 'set_abundance: particle not recognised')
      return
    end select

    abund = log10(element%abund) + 12.0
    if (relative_abundance) then
      a = abund + abnd
    else
      a = abnd
    endif
    if (abund /= a) then
      if (console_on(ch_hop, 1)) write(console_unit,'(A, A, A, F5.2, A, F5.2)') &
        'Changing abundance of ', trim(element%id), ' from ', abund, ' to ', a
      element%abund = 10.0**(a - 12.0)
      !if (.not. associated(mol)) element%awgt = awgt
      call self%compute_grams_per_hydrogen()
    endif
  END SUBROUTINE overwrite_abundance

  !=========================================================================
  SUBROUTINE add_hrest(self, chi, lam, i0, i1)
    class(composition_type), target :: self
    class(ion_type), pointer :: HI
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    real(dp), dimension(:), intent(in) :: lam
    integer, intent(in), optional :: i0, i1        ! optional cell range on dim 1 (threading)
    real(dp), dimension(:), allocatable :: lam3
    real(dp) :: hrest, fakt, e_kT, dxi
    real(dp), parameter :: nniv2 = 15.5**2
    integer :: x,y,z, xlo, xhi
    !---------------------------------------------------------------------
    call self%ion_by_id('H I', HI)
    allocate(lam3, mold=lam)
    lam3 = lam**3
    xlo = 1; xhi = size(self%temp,1)
    if (present(i0)) xlo = i0
    if (present(i1)) xhi = i1
    do z=1, size(self%temp,3)
      do y=1, size(self%temp,2)
        do x=xlo, xhi
          e_kT = self%e_kT(x,y,z)
          dxi = self%dxi(x,y,z)
          fakt = exp(-e_kT*(HI%xion - dxi) - HI%part(x,y,z))
          fakt = 2.0898E-26 / (2*RYD*e_kT) * fakt * HI%n(x,y,z)
          if (RYD >= nniv2 * dxi) then
            !NNIV=15
            !XNIV=15.
            !IF(H(1).LT.XNIV)NNIV=INT(H(1))
            !ADDF=EXP(TETA31/((FLOAT(NNIV)+0.5)**2)-UMC)-1.
            !IF(H(1).LT.XNIV+0.5)ADDF=0.
            !HREST(NTP)=EXPJ*ADDF
            hrest = exp(e_kT*(RYD/nniv2 - dxi))
            hrest = fakt * (hrest - 1)
            chi(:,x,y,z) = chi(:,x,y,z) + hrest*lam3
          endif
        end do
      end do
    end do
  END SUBROUTINE add_hrest

  !=========================================================================
  SUBROUTINE electron_scattering(self, chi, i0, i1)
    class(composition_type), target :: self
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    integer, intent(in), optional :: i0, i1        ! optional cell range on dim 1 (threading)
    real(dp), parameter :: thomson = 8./3.*PI*(EC**2/(EM*CC**2))**2
    real(dp) :: els
    integer :: x, y, z, lo, hi
    !---------------------------------------------------------------------
    lo = 1; hi = size(chi, 2)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    do z=1, size(chi, 4)
      do y=1, size(chi, 3)
        do x=lo, hi
          !els = 4.8206e-9 * self%pe / self%temp
          els = thomson * self%ne(x,y,z)
          chi(:,x,y,z) = chi(:,x,y,z) + els
        end do
      end do
    end do
  END SUBROUTINE electron_scattering

  !=========================================================================
  SUBROUTINE rayleigh_scattering(self, lam, chi, i0, i1)
    class(composition_type), target :: self
    class(molec_type), pointer :: H2
    class(ion_type), pointer :: HI, HE
    real(dp), dimension(:), intent(in) :: lam
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    integer, intent(in), optional :: i0, i1        ! optional cell range on dim 1 (threading)
    real(dp), dimension(:), allocatable :: xray2, HeIscat
    real(dp) :: xr, rayh, rayh2, nHI, nH2, nHE
    logical :: do_H2
    integer :: f, x, y, z, xlo, xhi
    !---------------------------------------------------------------------
    allocate(xray2, HeIscat,  mold=lam)
    xlo = 1; xhi = size(chi, 2)
    if (present(i0)) xlo = i0
    if (present(i1)) xhi = i1

    ! rayleigh scattering
    call self%ion_by_id('H I', HI)
    call self%ion_by_id('HE I', HE)
    call self%molecule_by_id('H2', H2)
    do_H2 = associated(H2)

    xray2 = 0.
    !where (lam > 911.27) xray2 = 1 / max(1026., lam)**2
    where (lam > 1400) xray2 = 1 / max(1026., lam)**2

    if (helium_scattering) then
      ! Scattering by HeI in cm**2/atom, from Bues & Wehrse, 1976, A&A 51, 461
      HeIscat = 0.66520e-24*4.*(500./lam)**4
    else
      HeIscat = 0
    endif

    do z=1, size(chi, 4)
      do y=1, size(chi, 3)
        do x=xlo, xhi
          nHI = HI%n(x,y,z)
          nHE = HE%n(x,y,z)
          if (do_H2) then
            nH2 = H2%n(x,y,z)
          else
            nH2 = 0
          endif

          do f=1, size(lam)
            xr = xray2(f)
            if (xr > 0.) then
              rayh  = (5.799e-13 + xr*(1.422e-6 + xr*2.784)) * nHI ! Dalgarno, quoted by Gingerich (1964)
              rayh2 = ( 8.14e-13 + xr*( 1.28e-6 + xr*1.61 )) * nH2 ! Dalgarno & Williams (1962)

              ! add sumsca to sumsca from lines
              chi(f,x,y,z) = chi(f,x,y,z) + (rayh + rayh2)*xr*xr
            endif

            chi(f,x,y,z) = chi(f,x,y,z) + HeIscat(f) * nHE
          end do
        end do
      end do
    end do
  END SUBROUTINE rayleigh_scattering

  !===============================================================================
  SUBROUTINE read_molec_file (self, molec_file)
    class(composition_type), target :: self
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    type(molec_type) :: molecule
    type(ion_type), pointer :: ion
    character(len=*) :: molec_file
    integer :: i, j, lunit, iostat
    !---------------------------------------------------------------------------
    self%nmolec = 0
    if (self%tmolim <= 0) return

    do i=1, self%nel
      element => self%elements(i)
      do j=lbound(element%ions,1), ubound(element%ions,1)
        ion => element%ions(j)
        if (allocated(ion%molecules)) then
          deallocate(ion%molecules)
          ion%nmolec = 0
        endif
      end do
    end do

    call os%check_file(molec_file)
    open(file=molec_file, newunit=lunit, iostat=iostat)
    do while (iostat==0)
      call molecule%read_molec_data(lunit, iostat)
      if (iostat==0) self%nmolec = self%nmolec + 1
    end do

    allocate(self%molecules(self%nmolec))
    rewind(lunit)

    do i=1, self%nmolec
      mol => self%molecules(i)
      call mol%read_molec_data(lunit, iostat)
      call self%parse_mol_id(mol)
      do j=1, mol%nions
        ion => mol%ions(j)%p
        ion%nmolec = ion%nmolec + 1
      end do
    end do

    close(lunit)

    do i=1, self%nel
      element => self%elements(i)
      do j=lbound(element%ions,1), ubound(element%ions,1)
        ion => element%ions(j)
        if (ion%nmolec > 0) then
          allocate(ion%molecules(ion%nmolec))
        endif
        ion%mcount = 0
      end do
    end do

    do i=1, self%nmolec
      mol => self%molecules(i)
      do j=1, mol%nions
        ion => mol%ions(j)%p
        ion%mcount = ion%mcount + 1
        ion%molecules(ion%mcount)%p => mol
      end do
    end do
  END SUBROUTINE read_molec_file

  !===============================================================================
  SUBROUTINE parse_mol_id (self, mol)
    class(composition_type) :: self
    class(molec_type) :: mol
    class(ion_type), pointer :: ion
    integer :: i, j, lenid, count, nions
    character(len=2) :: elem
    character :: c, next
    !-----------------------------------------------------------------------
    lenid = len_trim(mol%id)
    nions = 0
    mol%awgt = 0.0

    ! --- First pass: count how many atomic ions are present ---
    i = 1
    do while (i <= lenid)
      c = mol%id(i:i)
      if (i < lenid) then
        next = mol%id(i+1:i+1)
      else
        next = ' '
      end if
      if (c >= 'A' .and. c <= 'Z') then
        nions = nions + 1
        ! check next char for lowercase
        if (next >= 'a' .and. next <= 'z') then
          i = i + 1
        end if
      elseif (c >= '2' .and. c <= '9') then
        nions = nions + (ichar(c) - ichar('1'))
      elseif (c == '+' .or. c == '-') then
        ! handled later
      end if
      i = i + 1
    end do

    if (allocated(mol%ions)) deallocate(mol%ions)
    if (allocated(mol%z)) deallocate(mol%z)
    allocate(mol%ions(nions), mol%z(nions))
    mol%nions = nions

    ! --- Second pass: fill ions ---
    nions = 0
    i = 1
    do while (i <= lenid)
      c = mol%id(i:i)
      if (i < lenid) then
        next = mol%id(i+1:i+1)
      else
        next = ' '
      end if
      if (c >= 'A' .and. c <= 'Z') then
        elem = mol%id(i:i)
        if (next >= 'a' .and. next <= 'z') then
          elem = mol%id(i:i+1)
          i = i + 1
          if (i < lenid) then
            next = mol%id(i+1:i+1)
          else
            next = ' '
          end if
        end if

        count = 1
        if (next >= '2' .and. next <= '9') then
          count = ichar(next) - ichar('0')
          i = i + 1
        end if

        ! look up the ion
        call self%ion_by_id(trim(elem)//' I', ion)
        if (associated(ion)) then
          do j=1, count
            nions = nions + 1
            mol%ions(nions)%p => ion
            mol%z(nions) = findloc(element_ids, elem, dim=1)
            mol%awgt = mol%awgt + ion%awgt
          end do
        else
          call console_error ('element '//trim(elem)//' not found for molecule '//trim(mol%id))
          error stop
        endif
      elseif (c == '+') then
        mol%ions(nions)%p => ion%above
      elseif (c == '-') then
        mol%ions(nions)%p => ion%below
      end if
      i = i + 1
    end do
  END SUBROUTINE parse_mol_id

  !===============================================================================
  SUBROUTINE get_table_sp(self, sp_in, sp_out)
    class(composition_type) :: self
    type(spectrum_type), intent(in) :: sp_in
    type(spectrum_type), intent(out) :: sp_out
    real(dp), dimension(:), allocatable :: tab_lam
    logical, dimension(:), allocatable :: needed
    integer :: i, j, nnu, pos
    !---------------------------------------------------------------------------
    nnu = sum(self%trans%nnu)
    allocate(tab_lam(nnu))
    i=1
    do j=1, size(self%trans)
      associate(trans=>self%trans(j))
        tab_lam(i:i+trans%nnu-1) = trans%lam
        i = i + trans%nnu
      end associate
    end do
    tab_lam = max(1d0, tab_lam) ! remove zero values
    call sort_asc(tab_lam)
    call remove_close_entries(tab_lam, threshold=0.5d0)
    nnu = size(tab_lam)
    allocate(needed(nnu), source=.false.)
    do i=1, sp_in%nnu
      pos = searchsorted(tab_lam, sp_in%lam(i), outbounds_ok=.true.)
      needed(max(1, pos)) = .true.
      needed(min(pos+1, nnu)) = .true.
    end do
    call sp_out%reset(lam=pack(tab_lam, needed))
  END SUBROUTINE get_table_sp

  !===============================================================================
  SUBROUTINE print_opacity_contributions(self, nu)
    class(composition_type), target :: self
    class(bftrans_type), pointer :: trans
    class(ion_type), pointer :: HeI
    real(dp), dimension(:), intent(in) :: nu
    real(4), dimension(:,:,:,:), allocatable :: stim, tchi
    real(4), dimension(:,:,:,:,:), allocatable :: absk, scat
    real(dp), dimension(:), allocatable :: lam, hnu_e
    real(4) :: absk_tot, scat_tot
    character(len=1024) :: text
    integer, dimension(:), allocatable :: bin, nel, order
    character(6), dimension(:), allocatable :: key
    integer :: f, x, y, z, i, ib, nb, nx, ny, nz, nf, nt, lunit
    logical :: first_time = .true.
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('composition_type%print_opacity_contributions', itimer=itimer)
    nx = size(self%temp, 1)
    ny = size(self%temp, 2)
    nz = size(self%temp, 3)
    nf = size(nu, 1)
    nt = size(self%trans)
    allocate(stim(nf,nx,ny,nz))
    allocate(tchi(nf,nx,ny,nz))
    allocate(bin(nt))
    allocate(key(nt))
    allocate(nel(nt), order(nt))

    do i=1, nt
      trans=> self%trans(i)
      if (associated(trans%ion)) then
        nel(i) = findloc([element_ids, 'Me'], trans%ion%id(1:2), 1)
      else
        nel(i) = 0
      endif
      if (nel(i) == 0) nel(i) = nt+1
    end do
    call argsort_int(nel, order)

    nb = 1
    key(1) = 'Hrest'
    do i=1, nt
      trans => self%trans(order(i))
      ib = findloc(key(1:nb), trim(trans%key), 1)
      if (ib == 0) then
        nb = nb + 1
        ib = nb
        key(ib) = trim(trans%key)
      endif
      bin(i) = ib
    end do

    allocate(absk(nf, nx, ny, nz, nb), source=0.0)
    allocate(scat(nf, nx, ny, nz, 2), source=0.0)
    call self%ion_by_id('HE I', HeI)
    allocate(lam, source=CC / nu * 1e8)
    allocate(hnu_e, source=h_e * nu)

    do z=1, nz
      do y=1, ny
        do x=1, nx
          do f=1, nf
            stim(f,x,y,z) = 1. - exp(-hnu_e(f) * self%e_kT(x,y,z))
          end do
        end do
      end do
    end do

    call self%add_hrest(absk(:,:,:,:,1), lam)
    absk(:,:,:,:,1) = absk(:,:,:,:,1) * stim

    call self%electron_scattering(scat(:,:,:,:,1))
    call self%rayleigh_scattering(lam, scat(:,:,:,:,2))

    ! add absorption from all tables, which do not include stimulated emission
    do i=1, nt
      trans => self%trans(order(i))
      ib = bin(i)
      tchi = 0
      call trans%add_absorption(tchi, lam, self%temp, self%e_kT, self%teta, &
        self%pe, self%dxi, ne=self%ne, nhe=HeI%n, shield=self%shield, sorted=.true.)
      if (.not. trans%istim) tchi = tchi * stim
      absk(:,:,:,:,ib) = absk(:,:,:,:,ib) + tchi
    end do

    call io%lock%set('opacity_contributions')
    open(newunit=lunit, file=trim(io%run_dir)//'opacity_contributions.txt', access='append')
    if (first_time) then
      first_time = .false.
      write(lunit, '(6A11, 300(A7))') 'Temp', 'Pe', 'Pg', 'NH1', 'ABS/Rho', 'SCAT/Rho', &
        'Escat', 'Rscat', (trim(key(i)), i=1, nb)
    endif

    do z=1, nz
      associate(absk=>absk(1,1,1,z,:),   &
        scat=>scat(1,1,1,z,:), &
        rho=>self%rho(1,1,z), &
        t=>self%temp(1,1,z), &
        pe=>self%pe(1,1,z), &
        nh=>self%tothi(1,1,z), &
        pg=>self%pg(1,1,z))

        absk_tot = sum(absk)
        scat_tot = sum(scat)
        write(text, '(6ES11.3, 300(F7.2))') t, pe, pg, nh, &
          absk_tot/rho, scat_tot/rho, (scat(i)/scat_tot*100, i=1, 2), (absk(i)/absk_tot*100, i=1, nb)
      end associate
      call replace_substring(text, '-0.00', '     ')
      call replace_substring(text, ' 0.00', '     ')
      write(lunit, '(A)') trim(text)
    end do
    close(lunit)
    call io%lock%unset('opacity_contributions')
    call trace%end(itimer)
  END SUBROUTINE print_opacity_contributions
  !===============================================================================
  ! Write the used line-list lines of every ion and molecule to unit, in the
  ! linelist_out.txt format (linelist_mod%print_ll_particle).
  !-------------------------------------------------------------------------------
  SUBROUTINE write_linelists (self, unit)
    class(composition_type), target :: self
    integer, intent(in) :: unit
    integer :: i, j
    !-----------------------------------------------------------------------------
    if (allocated(self%elements)) then
      do i=1, self%nel
        if (.not. allocated(self%elements(i)%ions)) cycle
        if (.not. any(self%elements(i)%ions%in_linelist)) cycle
        do j=lbound(self%elements(i)%ions,1), ubound(self%elements(i)%ions,1)
          call write_particle(self%elements(i)%ions(j))
        end do
      end do
    end if
    if (allocated(self%molecules)) then
      do i=1, self%nmolec
        call write_particle(self%molecules(i))
      end do
    end if
  contains
    subroutine write_particle (ptc)
      class(particle_type), target :: ptc
      integer, allocatable :: iso_values(:)
      integer :: k, niso
      !---------------------------------------------------------------------------
      if (.not. allocated(ptc%ll)) return   ! no SoA linelists were built
      niso = 0
      if (associated(ptc%isotope)) niso = size(ptc%isotope)
      allocate(iso_values(niso))
      iso_values = [(ptc%isotope(k)%iso, k=1, niso)]
      select type (ptc)
      class is (molec_type)
        call print_ll_particle(unit, ptc%id, iso_values, ptc%ll, iso_digits=3*ptc%nions)
      class default
        call print_ll_particle(unit, ptc%id, iso_values, ptc%ll)
      end select
    end subroutine write_particle
  END SUBROUTINE write_linelists

END MODULE composition_mod
