!===============================================================================
!> C interface to the HOP EOS/opacity engine, for use from Python (ctypes).
!>
!> Design notes:
!> - A single composition_type instance per process.  HOP configuration lives
!>   partly in module-level globals (abundances, tolerance, ref_lam, ...), so
!>   two differently-configured instances cannot coexist anyway; a second
!>   hop_init without hop_free returns an error.
!> - All functions return an integer status: 0 = ok, > 0 = error (see body).
!>   Hard failures inside HOP (missing files, non-convergence) call stop/error
!>   stop and terminate the process -- the Python wrapper pre-validates what it
!>   can.
!> - Strings are NUL-terminated C strings.  Numeric arrays are allocated by the
!>   caller and filled here; nothing crosses the boundary by reference after a
!>   call returns.
!> - EOS quantities cross as float64 (HOP real(dp)); per-wavelength opacity
!>   buffers as float32 (HOP real(4)).
!===============================================================================
MODULE hop_c_api
  USE iso_c_binding
  USE composition_mod
  USE spectrum_mod, only: spectrum_type, lam_tolerance
  USE bbtrans_mod,  only: tolerance, disable_ga, disable_stark, disable_vdw, ga_fallback
  USE bftrans_mod,  only: H_occupation_probability
  USE ion_mod,      only: debye_shield, particle_type
  USE molecule_mod, only: irwin_dis
  USE hprofl_mod,   only: hydrogen_BPO
  USE hvcs_mod,     only: hydrogen_vcs_table, vcs_data_dir
  USE osmet_mod,    only: eosmet
  USE handyfuncs_mod, only: lowercase
  USE prec_mod,     only: dp
  implicit none
  private

  ! master is the initialized-but-never-loaded composition; each hop_load_state
  ! solves in a fresh copy of it (comp), mirroring the dispatch pattern
  ! (tab_comp = composition; reassign_pointer): load_state_variables is not
  ! re-entrant on one instance (ion_type keeps partially-deallocated state).
  type(composition_type), allocatable, save :: master, comp
  logical, save :: abundances_set = .false.
  ! verbose is the one option that must survive being set before hop_init: the
  ! messages it gates (abundance overrides, ...) are printed *during* injon
  integer, save :: init_verbose = 0

CONTAINS

  !=============================================================================
  ! Helpers
  !=============================================================================
  FUNCTION c2f_string(cs) RESULT(fs)
    character(kind=c_char), intent(in) :: cs(*)
    character(len=:), allocatable :: fs
    integer :: i, n
    !---------------------------------------------------------------------------
    n = 0
    do while (cs(n+1) /= c_null_char .and. n < 4096)
      n = n + 1
    end do
    allocate(character(len=n) :: fs)
    do i = 1, n
      fs(i:i) = cs(i)
    end do
  END FUNCTION c2f_string

  !-----------------------------------------------------------------------------
  ! Species lookup shared by hop_keep_species/hop_get_species.  Guards the bare
  ! element id case ('Fe' without an ionization stage): particle_by_id would
  ! crash on its unguarded second read, so only blank-free ids that name a
  ! molecule are allowed through.
  !-----------------------------------------------------------------------------
  SUBROUTINE species_lookup (c, nm, particle, status)
    class(composition_type), target :: c
    character(len=*), intent(in) :: nm
    class(particle_type), pointer, intent(out) :: particle
    integer(c_int), intent(out) :: status
    !---------------------------------------------------------------------------
    particle => null()
    status = 0
    if (index(trim(nm), ' ') == 0) then
      if (c%nmolec == 0) then
        status = 1
        return
      else if (findloc(c%molecules%id, trim(nm), 1) == 0) then
        status = 1   ! unknown molecule, or an element id without ' <stage>'
        return
      end if
    end if
    call c%particle_by_id(trim(nm), particle)
    if (.not. associated(particle)) status = 1
  END SUBROUTINE species_lookup

  !=============================================================================
  ! Lifecycle
  !=============================================================================
  FUNCTION hop_init(absdat_file, abund_file, isotope_file, molec_file, &
                    feh, tmolim, alpha_enh, r_process_enh) RESULT(status) &
                    BIND(c, name='hop_init')
    character(kind=c_char), intent(in) :: absdat_file(*), abund_file(*)
    character(kind=c_char), intent(in) :: isotope_file(*), molec_file(*)
    real(c_double), value :: feh
    real(c_float),  value :: tmolim, alpha_enh, r_process_enh
    integer(c_int) :: status
    character(len=:), allocatable :: absdat_f, abund_f, iso_f, molec_f
    !---------------------------------------------------------------------------
    if (allocated(master)) then
      status = 1   ! already initialized; call hop_free first
      return
    end if
    absdat_f = c2f_string(absdat_file)
    abund_f  = c2f_string(abund_file)
    iso_f    = c2f_string(isotope_file)
    molec_f  = c2f_string(molec_file)

    allocate(master)
    master%verbose = init_verbose           ! no default initializer in the type
    master%initialised = .false.
    if (.not. abundances_set) abundances = ''   ! module global, no initializer

    ! an unallocated allocatable actual makes the optional dummy absent (F2008)
    if (len(iso_f) == 0)   deallocate(iso_f)
    if (len(molec_f) == 0) deallocate(molec_f)
    call master%init(absdat_f, abund_f, feh, isotope_file=iso_f, molec_file=molec_f, &
      tmolim=tmolim, alpha_enhancement=alpha_enh, r_process_enhancement=r_process_enh)
    status = 0
  END FUNCTION hop_init

  FUNCTION hop_free() RESULT(status) BIND(c, name='hop_free')
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    if (allocated(comp)) deallocate(comp)
    if (allocated(master)) deallocate(master)
    status = 0
  END FUNCTION hop_free

  !=============================================================================
  ! Module-global options.  String/real options that influence injon
  ! (abundances, exclude_elements, ...) must be set BEFORE hop_init.
  !=============================================================================
  FUNCTION hop_set_option_str(name, value) RESULT(status) &
           BIND(c, name='hop_set_option_str')
    character(kind=c_char), intent(in) :: name(*), value(*)
    integer(c_int) :: status
    character(len=:), allocatable :: nm, val
    !---------------------------------------------------------------------------
    nm  = lowercase(c2f_string(name))
    val = c2f_string(value)
    status = 0
    select case (nm)
    case ('abundances')
      abundances = val
      abundances_set = .true.
    case ('exclude_elements')
      exclude_elements = val
    case ('only_elements')
      only_elements = val
    case ('vcs_data_dir')
      vcs_data_dir = val
    case default
      status = 1
    end select
  END FUNCTION hop_set_option_str

  FUNCTION hop_set_option_real(name, value) RESULT(status) &
           BIND(c, name='hop_set_option_real')
    character(kind=c_char), intent(in) :: name(*)
    real(c_double), value :: value
    integer(c_int) :: status
    character(len=:), allocatable :: nm
    !---------------------------------------------------------------------------
    nm = lowercase(c2f_string(name))
    status = 0
    select case (nm)
    case ('ref_lam')
      ref_lam = value
    case ('debug_opac')
      debug_opac = value
    case ('tolerance')
      tolerance = real(value, 4)
    case ('opac_table_vmic')
      opac_table_vmic = real(value, 4)
    case ('eosmet')
      eosmet = value
    case ('lam_tolerance')
      lam_tolerance = value
    case default
      status = 1
    end select
  END FUNCTION hop_set_option_real

  FUNCTION hop_set_option_int(name, value) RESULT(status) &
           BIND(c, name='hop_set_option_int')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: value
    integer(c_int) :: status
    character(len=:), allocatable :: nm
    !---------------------------------------------------------------------------
    nm = lowercase(c2f_string(name))
    status = 0
    select case (nm)
    case ('mhd_eos')
      mhd_eos = value /= 0
    case ('scattering_absorption')
      scattering_absorption = value /= 0
    case ('helium_scattering')
      helium_scattering = value /= 0
    case ('relative_abundances')
      relative_abundances = value /= 0
    case ('line_scattering')
      line_scattering = value /= 0
    case ('sample_linelist')
      sample_linelist = value /= 0
    case ('h_occupation_probability')
      H_occupation_probability = value /= 0
    case ('debye_shield')
      debye_shield = value /= 0
    case ('irwin_dis')
      irwin_dis = value /= 0
    case ('disable_ga')
      disable_ga = value /= 0
    case ('disable_stark')
      disable_stark = value /= 0
    case ('disable_vdw')
      disable_vdw = value /= 0
    case ('ga_fallback')
      ga_fallback = value /= 0
    case ('hydrogen_bpo')
      hydrogen_BPO = value /= 0
    case ('hydrogen_vcs_table')
      hydrogen_vcs_table = value /= 0
    case ('switch')
      switch = value
    case ('verbose')
      init_verbose = value
      if (allocated(master)) master%verbose = value
    case default
      status = 1
    end select
  END FUNCTION hop_set_option_int

  !=============================================================================
  ! EOS
  !=============================================================================
  FUNCTION hop_load_state(n, temp, ne, rho, vmic, use_rho, use_ne) &
           RESULT(status) BIND(c, name='hop_load_state')
    integer(c_int), value :: n, use_rho, use_ne
    real(c_double), intent(in)    :: temp(n), rho(n)
    real(c_double), intent(inout) :: ne(n)
    real(c_float),  value         :: vmic
    integer(c_int) :: status
    real(dp), allocatable :: t3(:,:,:), ne3(:,:,:), rho3(:,:,:)
    real(4),  allocatable :: vm3(:,:,:)
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(master)) then
      status = 2   ! hop_init not called
      return
    end if
    if (.not. master%initialised) then
      status = 2
      return
    end if
    if (n < 1) then
      status = 3
      return
    end if

    allocate(t3(n,1,1), ne3(n,1,1), rho3(n,1,1))
    t3(:,1,1)   = temp(1:n)
    ne3(:,1,1)  = ne(1:n)
    rho3(:,1,1) = rho(1:n)

    ! solve in a fresh copy of the never-loaded master (dispatch pattern)
    if (allocated(comp)) deallocate(comp)
    allocate(comp)
    comp = master
    call comp%reassign_pointer()

    if (vmic >= 0.0) then
      allocate(vm3(n,1,1))
      vm3 = vmic
      call comp%load_state_variables(t3, ne3, rho3, &
        use_rho = use_rho /= 0, use_ne = use_ne /= 0, vmic = vm3)
    else
      call comp%load_state_variables(t3, ne3, rho3, &
        use_rho = use_rho /= 0, use_ne = use_ne /= 0)
    end if

    ne(1:n) = comp%ne(:,1,1)
  END FUNCTION hop_load_state

  FUNCTION hop_get_n() RESULT(n) BIND(c, name='hop_get_n')
    integer(c_int) :: n
    !---------------------------------------------------------------------------
    n = 0
    if (allocated(comp)) then
      if (comp%loaded) n = size(comp%temp)
    end if
  END FUNCTION hop_get_n

  FUNCTION hop_get_field(name, n, out) RESULT(status) &
           BIND(c, name='hop_get_field')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: n
    real(c_double), intent(out) :: out(n)
    integer(c_int) :: status
    character(len=:), allocatable :: nm
    !---------------------------------------------------------------------------
    status = 0
    if (hop_get_n() /= n .or. n < 1) then
      status = 3   ! no state loaded or size mismatch
      return
    end if
    nm = lowercase(c2f_string(name))
    select case (nm)
    case ('temp')
      out(1:n) = reshape(comp%temp, [n])
    case ('ne')
      out(1:n) = reshape(comp%ne, [n])
    case ('rho')
      out(1:n) = reshape(comp%rho, [n])
    case ('pe')
      out(1:n) = reshape(comp%pe, [n])
    case ('pg')
      out(1:n) = reshape(comp%pg, [n])
    case ('e')
      out(1:n) = reshape(comp%E, [n])
    case ('eh')
      out(1:n) = reshape(comp%eh, [n])
    case ('chi500')
      out(1:n) = reshape(real(comp%chi500, dp), [n])
    case default
      status = 1
    end select
  END FUNCTION hop_get_field

  !-----------------------------------------------------------------------------
  ! Mark a species ('Fe I', 'CH', ...) so its number density survives the solve.
  ! Call after hop_init and before hop_load_state; the mark is set on the master
  ! instance, so every subsequent solve retains the species.  Mirrors dispatch's
  ! save_species pre-init marking: in_opac is read only by the deallocate keep-
  ! checks, so this does not change any opacity.
  !-----------------------------------------------------------------------------
  FUNCTION hop_keep_species(name) RESULT(status) BIND(c, name='hop_keep_species')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int) :: status
    class(particle_type), pointer :: particle
    !---------------------------------------------------------------------------
    if (.not. allocated(master)) then
      status = 2   ! hop_init not called
      return
    end if
    call species_lookup(master, c2f_string(name), particle, status)
    if (status == 0) particle%in_opac = .true.
  END FUNCTION hop_keep_species

  !-----------------------------------------------------------------------------
  ! Number density [1/cm3] of a species in the loaded state.  The species must
  ! have been kept (hop_keep_species) before the hop_load_state call -- or be
  ! one HOP retains anyway (continuum-opacity and line-list species).
  !-----------------------------------------------------------------------------
  FUNCTION hop_get_species(name, n, out) RESULT(status) &
           BIND(c, name='hop_get_species')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: n
    real(c_double), intent(out) :: out(n)
    integer(c_int) :: status
    class(particle_type), pointer :: particle
    !---------------------------------------------------------------------------
    if (hop_get_n() /= n .or. n < 1) then
      status = 3   ! no state loaded or size mismatch
      return
    end if
    call species_lookup(comp, c2f_string(name), particle, status)
    if (status /= 0) return
    if (.not. allocated(particle%n)) then
      status = 4   ! density not retained; call hop_keep_species before loading
      return
    end if
    out(1:n) = reshape(particle%n, [n])
  END FUNCTION hop_get_species

  FUNCTION hop_get_grph() RESULT(grph) BIND(c, name='hop_get_grph')
    real(c_double) :: grph
    !---------------------------------------------------------------------------
    grph = 0
    if (allocated(master)) grph = master%grph
  END FUNCTION hop_get_grph

  !=============================================================================
  ! Continuum opacity.  lam is vacuum Angstrom, ascending.  chi/sig buffers are
  ! float32 of size nlam*n (Fortran order: wavelength fastest), filled with
  ! chi in 1/cm (per-volume extinction; divide by rho for per-gram).
  !=============================================================================
  FUNCTION hop_absorption(nlam, lam, chi, include_absmet) RESULT(status) &
           BIND(c, name='hop_absorption')
    integer(c_int), value :: nlam, include_absmet
    real(c_double), intent(in) :: lam(nlam)
    real(c_float), intent(inout), target :: chi(*)
    integer(c_int) :: status
    type(spectrum_type) :: sp
    real(4), pointer, contiguous :: chi4(:,:,:,:)
    real(8), allocatable :: lam8(:)
    integer :: n
    !---------------------------------------------------------------------------
    status = 0
    n = hop_get_n()
    if (n < 1 .or. nlam < 1) then
      status = 3
      return
    end if
    allocate(lam8(nlam))
    lam8 = lam(1:nlam)
    call sp%reset(lam=lam8)
    sp%sorted = .true.
    chi4(1:nlam, 1:n, 1:1, 1:1) => chi(1:nlam*n)
    call comp%absorption(sp, chi4, include_absmet = include_absmet /= 0, &
                         sorted = .true.)
  END FUNCTION hop_absorption

  FUNCTION hop_scattering(nlam, lam, sig) RESULT(status) &
           BIND(c, name='hop_scattering')
    integer(c_int), value :: nlam
    real(c_double), intent(in) :: lam(nlam)
    real(c_float), intent(inout), target :: sig(*)
    integer(c_int) :: status
    real(4), pointer, contiguous :: sig4(:,:,:,:)
    real(8), allocatable :: lam8(:)
    integer :: n
    !---------------------------------------------------------------------------
    status = 0
    n = hop_get_n()
    if (n < 1 .or. nlam < 1) then
      status = 3
      return
    end if
    allocate(lam8(nlam))
    lam8 = lam(1:nlam)
    sig4(1:nlam, 1:n, 1:1, 1:1) => sig(1:nlam*n)
    sig4 = 0.0
    call comp%electron_scattering(sig4)
    call comp%rayleigh_scattering(lam8, sig4)
  END FUNCTION hop_scattering

END MODULE hop_c_api
