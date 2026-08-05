!===============================================================================
!> C interface to the HOP EOS/opacity engine, for use from Python (cffi).
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
  USE hprofl_mod,   only: hydrogen_BPO, holtsmark
  USE hvcs_mod,     only: hydrogen_vcs_table, vcs_data_dir
  USE linelist_mod, only: ll_profile, voigt_kernel
  USE osmet_mod,    only: eosmet, compute_epsilon_stage
  USE handyfuncs_mod, only: lowercase
  USE prec_mod,     only: dp
  USE const_mod,    only: h_k
#ifdef _OPENMP
  USE omp_lib,      only: omp_set_num_threads, omp_get_max_threads
#endif
  implicit none
  private

  ! master is the initialized-but-never-loaded composition; each hop_load_state
  ! solves in a fresh copy of it (comp), mirroring the dispatch pattern
  ! (tab_comp = composition; reassign_pointer): load_state_variables is not
  ! re-entrant on one instance (ion_type keeps partially-deallocated state).
  type(composition_type), allocatable, save :: master, comp
  logical, save :: abundances_set = .false.
  ! verbose is the one option that must survive being set before hop_init: the
  ! messages it gates (abundance overrides, ...) are printed *during* injon.
  ! The default is -1 rather than dispatch's 0, so that a library driven from
  ! Python says nothing at all unless asked; 0 restores what a dispatch run
  ! prints at its own default (the linelist line count), 1 and up add the rest.
  integer, save :: init_verbose = -1

  ! Message for the most recent nonzero status, readable with hop_last_error.  The
  ! status codes alone cannot say *which* option was rejected or *why* a value would
  ! not parse, and the wrapper cannot reconstruct that.  This does not cover the
  ! stop/error stop sites inside HOP itself, which still kill the process.
  character(len=512), save :: last_error = ''

  !> Every option hop_set_option accepts, in the order hop_option_names reports them.
  !> Kept adjacent to the select case below -- the two must agree, and
  !> tests/test_hop_abi.py checks that they do.
  character(len=*), parameter :: option_names = &
    'abundances exclude_elements only_elements vcs_data_dir ' // &
    'ref_lam debug_opac tolerance opac_table_vmic eosmet lam_tolerance ' // &
    'mhd_eos helium_scattering relative_abundances ' // &
    'line_scattering sample_linelist h_occupation_probability debye_shield ' // &
    'irwin_dis disable_ga disable_stark disable_vdw ga_fallback hydrogen_bpo ' // &
    'hydrogen_vcs_table ll_profile voigt_kernel switch verbose'

CONTAINS

  !=============================================================================
  ! Error reporting
  !=============================================================================
  !> Record msg and return code, so every failure path is one call.
  SUBROUTINE fail(msg, code, status)
    character(len=*), intent(in) :: msg
    integer(c_int), intent(in)   :: code
    integer(c_int), intent(out)  :: status
    !---------------------------------------------------------------------------
    last_error = msg
    status = code
  END SUBROUTINE fail

  !> Copy a Fortran string into a caller-provided C buffer, NUL-terminated.
  SUBROUTINE f2c_string(fs, buf, buflen, status)
    character(len=*), intent(in) :: fs
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), intent(in)  :: buflen
    integer(c_int), intent(out) :: status
    integer :: i, n
    !---------------------------------------------------------------------------
    status = 0
    n = len_trim(fs)
    if (buflen < 1) then
      status = 5
      return
    end if
    if (n + 1 > buflen) then
      status = 5           ! truncated; the caller still gets what fits
      n = buflen - 1
    end if
    do i = 1, n
      buf(i) = fs(i:i)
    end do
    buf(n+1) = c_null_char
  END SUBROUTINE f2c_string

  !> Message for the most recent nonzero status ('' if none).  Always succeeds
  !> unless the buffer is too small, in which case the text is truncated and the
  !> status is 5.
  FUNCTION hop_last_error(buf, buflen) RESULT(status) &
           BIND(c, name='hop_last_error')
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), value :: buflen
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    call f2c_string(last_error, buf, buflen, status)
  END FUNCTION hop_last_error

  !> Space-separated list of the names hop_set_option accepts, so a wrapper can
  !> report or validate them instead of hardcoding its own copy.
  FUNCTION hop_option_names(buf, buflen) RESULT(status) &
           BIND(c, name='hop_option_names')
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), value :: buflen
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    call f2c_string(option_names, buf, buflen, status)
  END FUNCTION hop_option_names

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
        call fail("unknown species '" // trim(nm) // "': no molecules loaded, and an " // &
                  "element id needs an ionization stage (e.g. 'Fe I')", 1, status)
        return
      else if (findloc(c%molecules%id, trim(nm), 1) == 0) then
        ! unknown molecule, or an element id without ' <stage>'
        call fail("unknown species '" // trim(nm) // "': not a molecule; an element " // &
                  "id needs an ionization stage (e.g. 'Fe I')", 1, status)
        return
      end if
    end if
    call c%particle_by_id(trim(nm), particle)
    if (.not. associated(particle)) &
      call fail("unknown species '" // trim(nm) // "'", 1, status)
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
      call fail('already initialized; call hop_free first', 1, status)
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
  !> Parse helpers.  Each leaves status alone on success and sets 4 on a value that
  !> will not read as the option's type, so the caller can report which it was.
  SUBROUTINE parse_real(s, out, status)
    character(len=*), intent(in) :: s
    real(dp), intent(out) :: out
    integer(c_int), intent(inout) :: status
    integer :: ios
    !---------------------------------------------------------------------------
    read(s, *, iostat=ios) out
    if (ios /= 0) call fail("cannot read '" // trim(s) // "' as a number", 4, status)
  END SUBROUTINE parse_real

  SUBROUTINE parse_int(s, out, status)
    character(len=*), intent(in) :: s
    integer, intent(out) :: out
    integer(c_int), intent(inout) :: status
    integer :: ios
    !---------------------------------------------------------------------------
    read(s, *, iostat=ios) out
    if (ios /= 0) call fail("cannot read '" // trim(s) // "' as an integer", 4, status)
  END SUBROUTINE parse_int

  !> Accepts what a namelist and a Python bool both stringify to.
  SUBROUTINE parse_logical(s, out, status)
    character(len=*), intent(in) :: s
    logical, intent(out) :: out
    integer(c_int), intent(inout) :: status
    character(len=:), allocatable :: t
    !---------------------------------------------------------------------------
    t = lowercase(trim(adjustl(s)))
    select case (t)
    case ('1', 't', 'true', '.true.', 'yes', 'on')
      out = .true.
    case ('0', 'f', 'false', '.false.', 'no', 'off')
      out = .false.
    case default
      call fail("cannot read '" // trim(s) // "' as a logical (use true/false or 1/0)", &
                4, status)
    end select
  END SUBROUTINE parse_logical

  !=============================================================================
  ! The one option setter.  Name and value both cross as strings and the parse
  ! happens here, against the option's own type, so no caller has to know which
  ! type bucket an option belongs to, and an unknown name stays distinguishable
  ! from a value that will not parse (see the status codes below).
  !
  ! Options that influence injon (abundances, exclude_elements, ...) must still
  ! be set BEFORE hop_init.
  !
  ! status: 0 ok, 1 unknown name, 4 value will not parse as that option's type
  !=============================================================================
  FUNCTION hop_set_option(name, value) RESULT(status) &
           BIND(c, name='hop_set_option')
    character(kind=c_char), intent(in) :: name(*), value(*)
    integer(c_int) :: status
    character(len=:), allocatable :: nm, val
    real(dp) :: rv
    integer  :: iv
    logical  :: lv
    !---------------------------------------------------------------------------
    nm  = lowercase(c2f_string(name))
    val = c2f_string(value)
    status = 0
    select case (nm)
    ! --- strings ------------------------------------------------------------
    case ('abundances')
      abundances = val
      abundances_set = .true.
    case ('exclude_elements')
      exclude_elements = val
    case ('only_elements')
      only_elements = val
    case ('vcs_data_dir')
      vcs_data_dir = val
    ! --- reals --------------------------------------------------------------
    case ('ref_lam')
      call parse_real(val, rv, status); if (status == 0) ref_lam = rv
    case ('debug_opac')
      call parse_real(val, rv, status); if (status == 0) debug_opac = rv
    case ('tolerance')
      call parse_real(val, rv, status); if (status == 0) tolerance = real(rv, 4)
    case ('opac_table_vmic')
      call parse_real(val, rv, status); if (status == 0) opac_table_vmic = real(rv, 4)
    case ('eosmet')
      call parse_real(val, rv, status); if (status == 0) eosmet = rv
    case ('lam_tolerance')
      call parse_real(val, rv, status); if (status == 0) lam_tolerance = rv
    ! --- logicals -----------------------------------------------------------
    case ('mhd_eos')
      call parse_logical(val, lv, status); if (status == 0) mhd_eos = lv
    case ('helium_scattering')
      call parse_logical(val, lv, status); if (status == 0) helium_scattering = lv
    case ('relative_abundances')
      call parse_logical(val, lv, status); if (status == 0) relative_abundances = lv
    case ('line_scattering')
      call parse_logical(val, lv, status); if (status == 0) line_scattering = lv
    case ('sample_linelist')
      call parse_logical(val, lv, status); if (status == 0) sample_linelist = lv
    case ('h_occupation_probability')
      call parse_logical(val, lv, status); if (status == 0) H_occupation_probability = lv
    case ('debye_shield')
      call parse_logical(val, lv, status); if (status == 0) debye_shield = lv
    case ('irwin_dis')
      call parse_logical(val, lv, status); if (status == 0) irwin_dis = lv
    case ('disable_ga')
      call parse_logical(val, lv, status); if (status == 0) disable_ga = lv
    case ('disable_stark')
      call parse_logical(val, lv, status); if (status == 0) disable_stark = lv
    case ('disable_vdw')
      call parse_logical(val, lv, status); if (status == 0) disable_vdw = lv
    case ('ga_fallback')
      call parse_logical(val, lv, status); if (status == 0) ga_fallback = lv
    case ('hydrogen_bpo')
      call parse_logical(val, lv, status); if (status == 0) hydrogen_BPO = lv
    case ('hydrogen_vcs_table')
      call parse_logical(val, lv, status); if (status == 0) hydrogen_vcs_table = lv
    case ('ll_profile')
      call parse_logical(val, lv, status); if (status == 0) ll_profile = lv
    ! --- integers -----------------------------------------------------------
    case ('voigt_kernel')
      call parse_int(val, iv, status); if (status == 0) voigt_kernel = iv
    case ('switch')
      call parse_int(val, iv, status); if (status == 0) switch = iv
    case ('verbose')
      call parse_int(val, iv, status)
      if (status == 0) then
        init_verbose = iv
        if (allocated(master)) master%verbose = iv
      end if
    case default
      call fail("unknown option '" // trim(nm) // "' (see hop_option_names)", 1, status)
    end select
  END FUNCTION hop_set_option

  !=============================================================================
  ! EOS
  !=============================================================================
  FUNCTION hop_load_state(n, temp, ne, rho, vmic, use_rho, use_ne) &
           RESULT(status) BIND(c, name='hop_load_state')
    integer(c_int), value :: n, use_rho, use_ne
    ! assumed-size, not temp(n): -fc-prototypes renders an explicit-shape dummy as a
    ! by-value scalar (`const double temp`) and only gets `(*)` right as a pointer, so
    ! the generated header is only trustworthy if every array dummy here is `(*)`
    real(c_double), intent(in)    :: temp(*), rho(*)
    real(c_double), intent(inout) :: ne(*)
    real(c_float),  value         :: vmic
    integer(c_int) :: status
    real(dp), allocatable :: t3(:,:,:), ne3(:,:,:), rho3(:,:,:)
    real(4),  allocatable :: vm3(:,:,:)
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(master)) then
      call fail('hop_load_state: no instance; call hop_init first', 2, status)
      return
    end if
    if (.not. master%initialised) then
      call fail('hop_load_state: instance not initialised', 2, status)
      return
    end if
    if (n < 1) then
      call fail('hop_load_state: n must be >= 1', 3, status)
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

    ! vmic is always passed, negative ("unset") included: dispatch's table paths
    ! do the same (tab_vmic = opac_table_vmic, which defaults to -1), and the line
    ! opacity reads comp%vmic unguarded -- leaving it unallocated segfaults there.
    ! A negative value is thermal-only broadening, see inverse_doppler_factor.
    allocate(vm3(n,1,1))
    vm3 = vmic
    call comp%load_state_variables(t3, ne3, rho3, &
      use_rho = use_rho /= 0, use_ne = use_ne /= 0, vmic = vm3)

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
    real(c_double), intent(out) :: out(*)
    integer(c_int) :: status
    character(len=:), allocatable :: nm
    !---------------------------------------------------------------------------
    status = 0
    if (hop_get_n() /= n .or. n < 1) then
      call fail('hop_get_field: no state loaded, or n differs from the loaded size', &
                3, status)
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
      call fail("hop_get_field: unknown field '" // trim(nm) // "'", 1, status)
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
      call fail('hop_keep_species: no instance; call hop_init first', 2, status)
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
    real(c_double), intent(out) :: out(*)
    integer(c_int) :: status
    class(particle_type), pointer :: particle
    !---------------------------------------------------------------------------
    if (hop_get_n() /= n .or. n < 1) then
      call fail('hop_get_species: no state loaded, or n differs from the loaded size', &
                3, status)
      return
    end if
    call species_lookup(comp, c2f_string(name), particle, status)
    if (status /= 0) return
    if (.not. allocated(particle%n)) then
      call fail('species density not retained; call hop_keep_species before ' // &
                'hop_load_state', 4, status)
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

  !-----------------------------------------------------------------------------
  ! Microturbulence [km/s] of the loaded state.  vmic enters only the line
  ! broadening (composition%add_linelist_opacities2), never the EOS, so it can
  ! be changed on a solved state instead of re-solving it.  A negative value is
  ! thermal-only broadening.
  !-----------------------------------------------------------------------------
  FUNCTION hop_set_vmic(n, vmic) RESULT(status) BIND(c, name='hop_set_vmic')
    integer(c_int), value :: n
    real(c_float), value :: vmic
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    status = 0
    if (hop_get_n() /= n .or. n < 1) then
      call fail('hop_set_vmic: no state loaded, or n differs from the loaded size', &
                3, status)
      return
    end if
    comp%vmic = vmic
  END FUNCTION hop_set_vmic

  !=============================================================================
  ! Opacity.  lam is vacuum Angstrom, ascending.  The chi and scat buffers are
  ! float32 of size nlam*n (Fortran order: wavelength fastest), filled in 1/cm
  ! (per-volume extinction; divide by rho for per-gram).  See
  ! composition_type%opacity for what each of them contains: chi is the total
  ! extinction (continuum + absmet + line list), scat its scattering part.
  ! line_scat < 0 leaves the split to the module-global line_scattering, 0/1
  ! override it; include_llist is honoured only once a line list has been read.
  !=============================================================================
  FUNCTION hop_opacity(nlam, lam, chi, scat, include_absmet, include_llist, dlam, &
                       line_scat) RESULT(status) BIND(c, name='hop_opacity')
    integer(c_int), value :: nlam, include_absmet, include_llist, line_scat
    real(c_double), intent(in) :: lam(*)
    real(c_float), intent(inout), target :: chi(*), scat(*)
    real(c_double), value :: dlam
    integer(c_int) :: status
    type(spectrum_type) :: sp
    real(4), pointer, contiguous :: chi4(:,:,:,:), scat4(:,:,:,:)
    real(8), allocatable :: lam8(:)
    integer :: n
    !---------------------------------------------------------------------------
    status = 0
    n = hop_get_n()
    if (n < 1 .or. nlam < 1) then
      call fail('hop_opacity: no state loaded, or nlam < 1', 3, status)
      return
    end if
    allocate(lam8(nlam))
    lam8 = lam(1:nlam)
    call sp%reset(lam=lam8)
    sp%sorted = .true.
    chi4 (1:nlam, 1:n, 1:1, 1:1) => chi (1:nlam*n)
    scat4(1:nlam, 1:n, 1:1, 1:1) => scat(1:nlam*n)
    if (line_scat < 0) then
      call comp%opacity(sp, chi4, scat4, include_absmet = include_absmet /= 0, &
        include_llist = include_llist /= 0, dlam = dlam, sorted = .true.)
    else
      call comp%opacity(sp, chi4, scat4, include_absmet = include_absmet /= 0, &
        include_llist = include_llist /= 0, dlam = dlam, sorted = .true., &
        line_scat = line_scat /= 0)
    end if
  END FUNCTION hop_opacity

  !=============================================================================
  ! Threaded Opacity
  !=============================================================================
  FUNCTION hop_opacity_omp(nlam, lam, chi, scat, include_absmet, include_llist, dlam, &
                           line_scat, nthreads) RESULT(status) BIND(c, name='hop_opacity_omp')
    integer(c_int), value :: nlam, include_absmet, include_llist, line_scat, nthreads
    real(c_double), intent(in) :: lam(*)
    real(c_float), intent(inout), target :: chi(*), scat(*)
    real(c_double), value :: dlam
    integer(c_int) :: status

    type(spectrum_type) :: sp
    real(4), pointer, contiguous :: chi4(:,:,:,:), scat4(:,:,:,:)
    real(8), allocatable :: lam8(:)
    integer :: n

    real(4), dimension(:,:,:,:), allocatable, target :: bet, sig
    real(4), dimension(:,:,:,:), allocatable :: ehat, chi_scat
    real(4), dimension(:,:,:,:,:), allocatable :: one_m_eps
    logical :: do_absmet, do_llist, do_split, sorted
    real(8) :: dl
    integer :: i, j, k, f, nf, m(3)
    real(4) :: s
    integer :: nch_l, nch_c, nthr, w, i0, i1

    status = 0
    n = hop_get_n()
    if (n < 1 .or. nlam < 1) then
      call fail('hop_opacity_omp: no state loaded, or nlam < 1', 3, status)
      return
    end if

#ifdef _OPENMP
    if (nthreads > 0) call omp_set_num_threads(nthreads)
    nthr = omp_get_max_threads()
#else
    nthr = 1
#endif
    ! increase the chunk multiplier to allow dynamic scheduling
    nch_l = min(n, max(1, nthr * 64))
    nch_c = min(n, max(1, nthr * 4))
#ifndef _OPENMP
    nch_l = 1                   
    nch_c = 1
#endif

    allocate(lam8(nlam))
    lam8 = lam(1:nlam)
    call sp%reset(lam=lam8)
    sp%sorted = .true.
    chi4 (1:nlam, 1:n, 1:1, 1:1) => chi (1:nlam*n)
    scat4(1:nlam, 1:n, 1:1, 1:1) => scat(1:nlam*n)

    m = shape(comp%temp)
    if (m(1) /= n) then
      call fail('hop_opacity_omp: state is not shaped (n,1,1)', 3, status)
      return
    end if
    nf = size(chi4, 1)
    dl = dlam
    do_absmet = .false.
    do_llist = (include_llist /= 0) .and. comp%ll_initialised
    do_split = (line_scat /= 0)
    sorted = .true.

    allocate(bet, mold=chi4)
    allocate(sig, mold=chi4)
    if (do_llist) then
      allocate(ehat, mold=chi4)
      if (do_split) then
        allocate(chi_scat, mold=chi4)
        allocate(one_m_eps(nf, m(1), m(2), m(3), 2))
      endif
    endif

    ! --- everything that precedes the line opacity, over cell chunks ----------
    !$omp parallel do schedule(dynamic) default(shared) private(w, i0, i1)
    do w = 1, nch_c
      call chunk_range(w, nch_c, n, i0, i1)
      if (i1 < i0) cycle
      call comp%boltzmann_factor(sp, bet(:,i0:i1,:,:), comp%temp(i0:i1,:,:))
      call comp%continuum_opacities(sp, chi4, sigma=sig, beta=bet, sorted=sorted, i0=i0, i1=i1)
      if (do_llist) then
        ehat(:,i0:i1,:,:) = chi4(:,i0:i1,:,:)
        if (do_split) then
          chi_scat(:,i0:i1,:,:) = 0
          call compute_epsilon_stage(1, sp%nu, comp%temp(i0:i1,:,:), comp%ne(i0:i1,:,:), one_m_eps(:,i0:i1,:,:,1))
          call compute_epsilon_stage(2, sp%nu, comp%temp(i0:i1,:,:), comp%ne(i0:i1,:,:), one_m_eps(:,i0:i1,:,:,2))
        endif
      endif
    end do
    !$omp end parallel do

    if (do_absmet) call comp%add_absmet_opacities(sp, chi4, bet, sigma=sig, split=do_split)

    ! --- line opacity --------------------------------------------------------
    if (do_llist) then
      if (do_split) then
        !$omp parallel do schedule(dynamic) default(shared) private(w, i0, i1)
        do w = 1, nch_l
          call chunk_range(w, nch_l, n, i0, i1)
          if (i1 < i0) cycle
          call comp%add_linelist_opacities2(sp, chi4, ehat, dl, one_m_eps=one_m_eps, chi_scat=chi_scat, beta=bet, i0=i0, i1=i1)
        end do
        !$omp end parallel do
      else
        !$omp parallel do schedule(dynamic) default(shared) private(w, i0, i1)
        do w = 1, nch_l
          call chunk_range(w, nch_l, n, i0, i1)
          if (i1 < i0) cycle
          call comp%add_linelist_opacities2(sp, chi4, ehat, dl, i0=i0, i1=i1)
        end do
        !$omp end parallel do
      endif
    endif

    !$omp parallel do collapse(3) schedule(static) default(shared) private(f, s)
    do k=1, m(3)
      do j=1, m(2)
        do i=1, m(1)
          if (do_llist) then
            if (do_split) then
              do f=1, nf
                s = sig(f,i,j,k) + chi_scat(f,i,j,k)
                chi4 (f,i,j,k) = chi4(f,i,j,k) - ehat(f,i,j,k)*bet(f,i,j,k) + s
                scat4(f,i,j,k) = s
              end do
            else
              do f=1, nf
                s = sig(f,i,j,k)
                chi4 (f,i,j,k) = chi4(f,i,j,k) - ehat(f,i,j,k)*bet(f,i,j,k) + s
                scat4(f,i,j,k) = s
              end do
            endif
          else
            do f=1, nf
              s = sig(f,i,j,k)
              chi4 (f,i,j,k) = chi4(f,i,j,k) * (1 - bet(f,i,j,k)) + s
              scat4(f,i,j,k) = s
            end do
          endif
        end do
      end do
    end do
    !$omp end parallel do

  CONTAINS
    !> Cell range [i0, i1] of chunk w of nch over n cells, balanced to +-1 cell.
    !> Empty (i1 < i0) when there are fewer cells than chunks.
    PURE SUBROUTINE chunk_range(w, nch, n, i0, i1)
      integer, intent(in)  :: w, nch, n
      integer, intent(out) :: i0, i1
      integer :: base, rem
      !-------------------------------------------------------------------------
      base = n / nch
      rem  = mod(n, nch)
      i0 = (w-1)*base + min(w-1, rem) + 1
      i1 =  w   *base + min(w,   rem)
    END SUBROUTINE chunk_range
  END FUNCTION hop_opacity_omp

  !-----------------------------------------------------------------------------
  ! Read one line-list file into the master instance, keeping only lines inside
  ! [lam(1)-dlam, lam(nlam)+dlam] (pass nlam=0 for no wavelength mask).  Mirrors
  ! the composition%init_linelist2 calls dispatch makes at setup, one per file in
  ! linelist_folder; repeated calls accumulate.  Must be called after hop_init and
  ! BEFORE hop_load_state: the solve works on a fresh copy of master, and
  ! particle%ll is allocatable, so the copy carries the lines along.
  !-----------------------------------------------------------------------------
  FUNCTION hop_init_linelist(linelist_file, dlam, nlam, lam) RESULT(status) &
           BIND(c, name='hop_init_linelist')
    character(kind=c_char), intent(in) :: linelist_file(*)
    real(c_double), value :: dlam
    integer(c_int), value :: nlam
    real(c_double), intent(in) :: lam(*)
    integer(c_int) :: status
    real(8), allocatable :: lam8(:)
    character(len=:), allocatable :: file_f
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(master)) then
      call fail('hop_init_linelist: no instance; call hop_init first', 2, status)
      return
    end if
    if (.not. master%initialised) then
      call fail('hop_init_linelist: instance not initialised', 2, status)
      return
    end if
    file_f = c2f_string(linelist_file)
    if (len(file_f) == 0) then
      call fail('hop_init_linelist: empty file name', 3, status)
      return
    end if
    if (nlam > 0) then
      allocate(lam8(nlam))
      lam8 = lam(1:nlam)
      call master%init_linelist2(file_f, dlam, lam8)
    else
      call master%init_linelist2(file_f, dlam)
    end if
  END FUNCTION hop_init_linelist

  !-----------------------------------------------------------------------------
  ! Griem's tabulated quasistatic Stark profile S(beta, p) for the hydrogen line
  ! (n, m), evaluated over a vector of beta.  HOP itself calls holtsmark() one
  ! beta at a time from stark1 (hprofl_mod); this batches it for the hydrogen-
  ! Stark study code, which fits a component pattern to S.  Needs no instance --
  ! holtsmark reads only its own tabulated data.
  !-----------------------------------------------------------------------------
  FUNCTION hop_holtsmark(nbeta, beta, p, n, m, out) RESULT(status) &
           BIND(c, name='hop_holtsmark')
    integer(c_int), value :: nbeta, n, m
    real(c_double), intent(in)  :: beta(*)
    real(c_double), value       :: p
    real(c_double), intent(out) :: out(*)
    integer(c_int) :: status
    integer :: i
    !---------------------------------------------------------------------------
    status = 0
    if (nbeta < 1) then
      call fail('hop_holtsmark: nbeta must be >= 1', 3, status)
      return
    end if
    do i = 1, nbeta
      out(i) = real(holtsmark(beta(i), p, n, m), c_double)
    end do
  END FUNCTION hop_holtsmark

  FUNCTION hop_scattering(nlam, lam, sig) RESULT(status) &
           BIND(c, name='hop_scattering')
    integer(c_int), value :: nlam
    real(c_double), intent(in) :: lam(*)
    real(c_float), intent(inout), target :: sig(*)
    integer(c_int) :: status
    real(4), pointer, contiguous :: sig4(:,:,:,:)
    real(8), allocatable :: lam8(:)
    integer :: n
    !---------------------------------------------------------------------------
    status = 0
    n = hop_get_n()
    if (n < 1 .or. nlam < 1) then
      call fail('hop_scattering: no state loaded, or nlam < 1', 3, status)
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
