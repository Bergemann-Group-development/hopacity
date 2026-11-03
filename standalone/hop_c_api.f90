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
  USE console_mod
  USE composition_mod
  USE spectrum_mod, only: spectrum_type, lam_tolerance
  USE bbtrans_mod,  only: tolerance, line_reject_frac, disable_rad, disable_stk, &
                          disable_vdw, molecule_damping_estimates
  USE bftrans_mod,  only: H_occupation_probability
  USE ion_mod,      only: debye_shield, particle_type
  USE molecule_mod, only: irwin_dis
  USE hprofl_mod,   only: hydrogen_BPO
  USE hvcs_mod,     only: hydrogen_profiles, vcs_data_dir
  USE anstee_mod,   only: anstee_file
  USE linelist_mod, only: ll_profile, voigt_kernel
  USE op_math_mod,  only: interval_t, merge_intervals
  USE interpol_mod, only: optical_depth, pchip_integral
  USE qlg,          only: qlg_init
  USE osmet_mod,    only: osmet, eosmet, compute_epsilon_stage
  USE handyfuncs_mod, only: lowercase, vac2air
  USE prec_mod,     only: dp
  USE const_mod,    only: h_k
  USE atmos_reader_mod, only: atmos_data_t, read_atmos_file, atmos_format_from_name
#ifdef _OPENMP
  USE omp_lib,      only: omp_set_num_threads, omp_get_max_threads, &
                          omp_pause_resource_all, omp_pause_soft
#endif
  implicit none
  private
  ! The C symbols are what Python binds; the same entry points are also USEd by the
  ! hop_prune.x program in this directory, so it goes through exactly the calls the
  ! Python wrapper makes.  Every BIND(C) procedure is public: gfortran 16 gives a
  ! private one local linkage, which drops it from the library's exports.
  public :: hop_init, hop_free, hop_set_option, hop_set_external_abundances, &
            hop_option_names, hop_last_error, hop_omp_release
  public :: hop_init_linelist, hop_init_linelist_windows, hop_keep_species, &
            hop_linelist_count, hop_get_linelist, hop_linelist_species, hop_write_linelist
  public :: hop_load_state, hop_set_vmic, hop_get_n, hop_get_field, hop_get_species, &
            hop_get_grph, hop_opacity, hop_opacity_omp, hop_scattering
  public :: hop_atmos_read, hop_atmos_get, hop_atmos_free, hop_optical_depth, hop_pchip_integral

  !> master is the initialized-but-never-loaded composition; each hop_load_state
  !> solves in a fresh copy of it (comp), mirroring the dispatch pattern
  !> (tab_comp = composition; reassign_pointer): load_state_variables is not
  !> re-entrant on one instance (ion_type keeps partially-deallocated state).
  type(composition_type), allocatable, save :: master, comp
  logical, save :: abundances_set = .false.
  !> verbose is HOP's console level (see console_mod) and the one option that must
  !> survive being set before hop_init: the messages it gates (abundance
  !> overrides, ...) are printed *during* init. The default is -1 rather than
  !> dispatch's 1, so that a library driven from Python says nothing at all
  !> unless asked; 0 adds warnings, 1 the setup summaries a dispatch run prints
  !> by default (the linelist line count, ...), 2 and up the diagnostics.
  integer, save :: init_verbose = -1
  !> Abundances an atmosphere file supplies (MARCS carries log eps for Z=1..92),
  !> applied by hop_init between the [Fe/H]-scaled standard set and the explicit
  !> 'abundances' overrides -- the order composition_type%init applies them in a
  !> dispatch run with use_abnd=T.  Set before hop_init; cleared by hop_free.
  real(4), allocatable, save :: external_abund(:)

  !> Message for the most recent nonzero status, readable with hop_last_error.  The
  !> status codes alone cannot say *which* option was rejected or *why* a value would
  !> not parse, and the wrapper cannot reconstruct that.  This does not cover the
  !> stop/error stop sites inside HOP itself, which still kill the process.
  character(len=512), save :: last_error = ''
  integer, parameter :: MIN_CHUNK_CELLS = 1

  !> The model atmosphere most recently parsed by hop_atmos_read, served field by
  !> field to hop_atmos_get.  Independent of the composition instance.
  type(atmos_data_t), allocatable, save :: atmos_file_data

  !> The per-line strength collection, filled by hop_linelist_count and served to
  !> hop_get_linelist / hop_linelist_species.  Collecting on the count call is what
  !> makes the pair consistent: the caller has to ask how many lines there are before
  !> it can size its buffers, and everything it then reads describes that same walk.
  real(8), allocatable, save :: ll_lam(:), ll_evi(:), ll_loggf(:), ll_strength(:)
  integer, allocatable, save :: ll_ispec(:)
  character(len=:), allocatable, save :: absmet_loaded   !< the file osmet holds
  character(len=8), allocatable, save :: ll_species(:)
  integer, save :: ll_nspecies = 0

  !> Every option hop_set_option accepts, in the order hop_option_names reports them.
  !> Kept adjacent to the select case below -- the two must agree, and
  !> tests/test_hop_abi.py checks that they do.
  character(len=*), parameter :: option_names = &
    'abundances exclude_elements only_elements vcs_data_dir anstee_file absmet_file ' // &
    'ref_lam debug_opac tolerance line_reject_frac opac_table_vmic eosmet lam_tolerance ' // &
    'mhd_eos helium_scattering relative_abundances ' // &
    'line_scattering sample_linelist h_occupation_probability debye_shield ' // &
    'irwin_dis disable_rad disable_stk disable_vdw molecule_damping_estimates hydrogen_bpo ' // &
    'hydrogen_profiles ll_profile voigt_kernel verbose'

CONTAINS

  !=============================================================================
  ! Error reporting
  !=============================================================================
  !> Record msg and return code, so every failure path is one call.
  !=============================================================================
  SUBROUTINE fail(msg, code, status)
    character(len=*), intent(in) :: msg
    integer(c_int), intent(in)   :: code
    integer(c_int), intent(out)  :: status
    !---------------------------------------------------------------------------
    last_error = msg
    status = code
  END SUBROUTINE fail

  !=============================================================================
  !> Copy a Fortran string into a caller-provided C buffer, NUL-terminated.
  !=============================================================================
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

  !=============================================================================
  !> Message for the most recent nonzero status ('' if none).  Always succeeds
  !> unless the buffer is too small, in which case the text is truncated and the
  !> status is 5.
  !=============================================================================
  FUNCTION hop_last_error(buf, buflen) RESULT(status) &
           BIND(c, name='hop_last_error')
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), value :: buflen
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    call f2c_string(last_error, buf, buflen, status)
  END FUNCTION hop_last_error

  !=============================================================================
  !> Space-separated list of the names hop_set_option accepts, so a wrapper can
  !> report or validate them instead of hardcoding its own copy.
  !=============================================================================
  FUNCTION hop_option_names(buf, buflen) RESULT(status) &
           BIND(c, name='hop_option_names')
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), value :: buflen
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    call f2c_string(option_names, buf, buflen, status)
  END FUNCTION hop_option_names

  !=============================================================================
  ! Per-line strengths -- what the last hop_opacity call worked out about every line
  ! it touched.
  !
  ! composition_type%opacity already computes, for each line, the maximum over the
  ! solved cells of its line-centre extinction over the continuum extinction
  ! (linelist_mod: ll%strength(l) = max(..., voigt_peak*Rij / cont(f0,i,j,k))).  That
  ! number is the only line ranking that compares across species -- it carries the
  ! abundance, the ionization fraction, the dissociation equilibrium and the partition
  ! function, none of which a gf-and-Boltzmann estimate has -- and until now the only
  ! way to see it was the linelist_out.txt diagnostic a full run writes.  These three
  ! functions hand back the same numbers from a solved state in memory, so a caller
  ! can rank the lines of a window at an arbitrary (T, rho) without a run and without
  ! keeping a file.
  !
  ! Usage: hop_load_state, then hop_opacity over the window of interest (the strengths
  ! are a by-product of that call and are zero before it), then hop_linelist_count to
  ! collect and size, then hop_get_linelist / hop_linelist_species to read.  A line is
  ! reported when it is active and was actually used, which is the same set
  ! linelist_out.txt writes.
  !=============================================================================
  !> Write the reported lines as a line list in the *input* format, so a caller can
  !> use them as the linelist_file of a later run. The records are the ones
  !> linelist_out.txt carries, written by the same routines a full run writes its dump
  !> with. The trailing strength column is a diagnostic both readers ignore; the reach
  !> column after it lets a run select lines per output patch (linelist_mod).
  !>
  !> What makes this useful: a curve of growth repeats one atmosphere at many
  !> abundances, and every rung then re-reads and re-selects from the whole master
  !> list. Selecting once per atmosphere -- at the strongest rung, so the result is a
  !> superset of what any weaker one needs -- and handing the cropped list to every
  !> rung makes the selection a property of the atmosphere rather than of each run.
  !>
  !> Which lines appear follows how the list was *loaded*, not the wavelengths of the
  !> opacity call: a line is reported when it was active and used, and the load is what
  !> decides the pool.  Crop with hop_init_linelist_windows, therefore --
  !> hop_init_linelist masks on the hull of its wavelengths and would put the whole
  !> span in the file.
  !>
  !> Model-atom transitions are not included: a standalone HOP state has none.
  !=============================================================================
  FUNCTION hop_write_linelist(path) RESULT(status) BIND(c, name='hop_write_linelist')
    character(kind=c_char), intent(in) :: path(*)
    integer(c_int) :: status
    integer :: lunit, iostat
    character(len=:), allocatable :: file
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(comp)) then
      call fail('hop_write_linelist: no composition; call hop_init first', 3, status)
      return
    end if
    call c2f_string(path, file)
    open(newunit=lunit, file=file, form='formatted', action='write', &
         status='replace', iostat=iostat)
    if (iostat /= 0) then
      call fail('hop_write_linelist: cannot open the output file', 3, status)
      return
    end if
    call comp%write_linelists(lunit)
    close(lunit)
  END FUNCTION hop_write_linelist

  !-----------------------------------------------------------------------------
  FUNCTION hop_linelist_count() RESULT(nout) BIND(c, name='hop_linelist_count')
    integer(c_int) :: nout
    !---------------------------------------------------------------------------
    call collect_lines(nout)
  END FUNCTION hop_linelist_count

  !-----------------------------------------------------------------------------
  !> The collection itself.  lam_air is Angstrom (as in the line-list files), evi the
  !> lower-level energy [eV], strength the ratio above, and ispec a 0-based index into
  !> the list hop_linelist_species reports.
  !-----------------------------------------------------------------------------
  FUNCTION hop_get_linelist(n, lam_air, evi, loggf, strength, ispec) RESULT(status) &
           BIND(c, name='hop_get_linelist')
    integer(c_int), value :: n
    real(c_double), intent(out) :: lam_air(*), evi(*), loggf(*), strength(*)
    integer(c_int), intent(out) :: ispec(*)
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(ll_lam)) then
      call fail('hop_get_linelist: nothing collected; call hop_linelist_count first', &
                3, status)
      return
    end if
    if (n /= size(ll_lam)) then
      call fail('hop_get_linelist: n differs from the collected line count', 3, status)
      return
    end if
    if (n < 1) return
    lam_air (1:n) = ll_lam(1:n)
    evi     (1:n) = ll_evi(1:n)
    loggf   (1:n) = ll_loggf(1:n)
    strength(1:n) = ll_strength(1:n)
    ispec   (1:n) = ll_ispec(1:n) - 1        ! 0-based for the C caller
  END FUNCTION hop_get_linelist

  !-----------------------------------------------------------------------------
  !> The species ids the ispec column indexes, newline separated in index order.  Ids
  !> carry the ionization stage ('Fe I', 'Fe II', 'TiO'), so they cannot be separated
  !> by blanks the way hop_option_names does it.
  !-----------------------------------------------------------------------------
  FUNCTION hop_linelist_species(buf, buflen) RESULT(status) &
           BIND(c, name='hop_linelist_species')
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), value :: buflen
    integer(c_int) :: status
    character(len=:), allocatable :: text, longer
    integer :: i
    !---------------------------------------------------------------------------
    if (.not. allocated(ll_species) .or. ll_nspecies < 1) then
      allocate(character(len=0) :: text)
    else
      allocate(text, source=trim(ll_species(1)))
      do i = 2, ll_nspecies
        allocate(longer, source=text // new_line('a') // trim(ll_species(i)))
        call move_alloc(longer, text)
      end do
    end if
    call f2c_string(text, buf, buflen, status)
  END FUNCTION hop_linelist_species

  !-----------------------------------------------------------------------------
  !> Walk the solved composition and gather every reported line.  The order is the one
  !> m3d_patch_t%print_linelist writes linelist_out.txt in -- each element's ions, then
  !> the molecules -- so the two agree line for line, which is what makes the dump a
  !> usable check on this.
  !-----------------------------------------------------------------------------
  SUBROUTINE collect_lines(nout)
    integer(c_int), intent(out) :: nout
    integer :: i, j, pass, ntot
    !---------------------------------------------------------------------------
    nout = 0
    if (allocated(ll_lam)) deallocate(ll_lam, ll_evi, ll_loggf, ll_strength, ll_ispec)
    if (allocated(ll_species)) deallocate(ll_species)
    ll_nspecies = 0
    if (.not. allocated(comp)) return
    if (.not. allocated(comp%elements)) return

    ! pass 1 counts so the buffers can be sized exactly, pass 2 fills them
    do pass = 1, 2
      ntot = 0
      ll_nspecies = 0
      do i = 1, comp%nel
        if (.not. allocated(comp%elements(i)%ions)) cycle
        do j = lbound(comp%elements(i)%ions, 1), ubound(comp%elements(i)%ions, 1)
          call take_particle(comp%elements(i)%ions(j), pass, ntot)
        end do
      end do
      if (allocated(comp%molecules)) then
        do i = 1, comp%nmolec
          call take_particle(comp%molecules(i), pass, ntot)
        end do
      end if
      if (pass == 1) then
        if (ntot < 1) return
        allocate(ll_lam(ntot), ll_evi(ntot), ll_loggf(ntot), ll_strength(ntot), &
                 ll_ispec(ntot))
        allocate(ll_species(max(ll_nspecies, 1)))
      end if
    end do
    nout = ntot
  END SUBROUTINE collect_lines

  !-----------------------------------------------------------------------------
  !> One particle's reported lines, from its line lists (one per isotope).
  !-----------------------------------------------------------------------------
  SUBROUTINE take_particle(ptc, pass, ntot)
    class(particle_type), intent(in) :: ptc
    integer, intent(in) :: pass
    integer, intent(inout) :: ntot
    integer :: k, l, ispec
    !---------------------------------------------------------------------------
    ispec = 0
    if (.not. allocated(ptc%ll)) return
    do k = lbound(ptc%ll, 1), ubound(ptc%ll, 1)
      associate(sll => ptc%ll(k))
        if (.not. allocated(sll%strength)) cycle
        do l = 1, sll%nline
          if (.not. (sll%active(l) .and. sll%used(l))) cycle
          call emit(ptc%id, ispec, pass, ntot, vac2air(exp(sll%llam0(l))), &
                    dble(sll%lower_ev(l)), dble(sll%loggf(l)), dble(sll%strength(l)))
        end do
      end associate
    end do
  END SUBROUTINE take_particle

  !-----------------------------------------------------------------------------
  !> Record one line.  ispec is the particle's index in ll_species, claimed lazily so a
  !> species that reports nothing takes no slot.
  !-----------------------------------------------------------------------------
  SUBROUTINE emit(id, ispec, pass, ntot, lam_air, evi, loggf, strength)
    character(len=*), intent(in) :: id
    integer, intent(inout) :: ispec, ntot
    integer, intent(in) :: pass
    real(8), intent(in) :: lam_air, evi, loggf, strength
    !---------------------------------------------------------------------------
    ntot = ntot + 1
    if (ispec == 0) then
      ll_nspecies = ll_nspecies + 1
      ispec = ll_nspecies
      if (pass == 2) ll_species(ispec) = id
    end if
    if (pass /= 2) return
    ll_lam(ntot) = lam_air
    ll_evi(ntot) = evi
    ll_loggf(ntot) = loggf
    ll_strength(ntot) = strength
    ll_ispec(ntot) = ispec
  END SUBROUTINE emit

  !=============================================================================
  ! Helpers
  !=============================================================================
  !> A subroutine rather than a function with a deferred-length result: built with
  !> -fno-realloc-lhs, gfortran hands such a result to callers in this module with
  !> length 0.
  !=============================================================================
  SUBROUTINE c2f_string(cs, fs)
    character(kind=c_char), intent(in) :: cs(*)
    character(len=:), allocatable, intent(out) :: fs
    integer :: i, n
    !---------------------------------------------------------------------------
    n = 0
    do while (n < 4096)
      if (cs(n+1) == c_null_char) exit
      n = n + 1
    end do
    allocate(character(len=n) :: fs)
    do i = 1, n
      fs(i:i) = cs(i)
    end do
  END SUBROUTINE c2f_string

  !-----------------------------------------------------------------------------
  !> Species lookup shared by hop_keep_species/hop_get_species.  Guards the bare
  !> element id case ('Fe' without an ionization stage): particle_by_id would
  !> crash on its unguarded second read, so only blank-free ids that name a
  !> molecule are allowed through.
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
  !> Lifecycle
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
    call c2f_string(absdat_file, absdat_f)
    call c2f_string(abund_file, abund_f)
    call c2f_string(isotope_file, iso_f)
    call c2f_string(molec_file, molec_f)

    ! QLG-bot's lookup tables. Must happen before any opacity or rate is accumulated:
    ! the kernel is PURE and cannot check, so a missed call is a silently zero profile
    ! rather than an error. Idempotent, a few ms, and independent of everything below.
    call qlg_init

    allocate(master)
    call console_set_level (init_verbose, ch_hop)
    master%initialised = .false.
    if (.not. abundances_set) abundances = ''   ! module global, no initializer

    ! an unallocated allocatable actual makes the optional dummy absent (F2008)
    if (len(iso_f) == 0)   deallocate(iso_f)
    if (len(molec_f) == 0) deallocate(molec_f)
    ! external_abund is passed as an unallocated actual -> absent optional dummy
    call master%init(absdat_f, abund_f, feh, isotope_file=iso_f, molec_file=molec_f, &
      tmolim=tmolim, external_abund=external_abund, alpha_enhancement=alpha_enh, &
      r_process_enhancement=r_process_enh)
    status = 0
  END FUNCTION hop_init

  !-----------------------------------------------------------------------------
  !> Atmosphere-provided abundances, log eps indexed by atomic number, for the next
  !> hop_init: what a dispatch run applies with atmos_params use_abnd=T.  Elements
  !> beyond n keep the standard value.  n=0 clears a previous vector.
  !-----------------------------------------------------------------------------
  FUNCTION hop_set_external_abundances(n, abund) RESULT(status) &
           BIND(c, name='hop_set_external_abundances')
    integer(c_int), value :: n
    real(c_double), intent(in) :: abund(*)
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    status = 0
    if (allocated(master)) then
      call fail('hop_set_external_abundances: set before hop_init', 1, status)
      return
    end if
    if (allocated(external_abund)) deallocate(external_abund)
    if (n < 0 .or. n > 92) then
      call fail('hop_set_external_abundances: n must be 0..92', 3, status)
      return
    end if
    if (n == 0) return
    allocate(external_abund(n))
    external_abund = real(abund(1:n), 4)
  END FUNCTION hop_set_external_abundances

  FUNCTION hop_free() RESULT(status) BIND(c, name='hop_free')
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    if (allocated(comp)) deallocate(comp)
    if (allocated(master)) deallocate(master)
    if (allocated(external_abund)) deallocate(external_abund)
    status = hop_omp_release()
    status = 0
  END FUNCTION hop_free

  FUNCTION hop_omp_release() RESULT(status) BIND(c, name='hop_omp_release')
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    status = 0
#ifdef _OPENMP
    status = int(omp_pause_resource_all(omp_pause_soft), c_int)
#endif
  END FUNCTION hop_omp_release

  !=============================================================================
  ! Module-global options.  String/real options that influence injon
  ! (abundances, exclude_elements, ...) must be set BEFORE hop_init.
  !=============================================================================
  !> Parse helpers.  Each leaves status alone on success and sets 4 on a value that
  !> will not read as the option's type, so the caller can report which it was.
  !=============================================================================
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

  !=============================================================================
  !> Accepts what a namelist and a Python bool both stringify to.
  !=============================================================================
  SUBROUTINE parse_logical(s, out, status)
    character(len=*), intent(in) :: s
    logical, intent(out) :: out
    integer(c_int), intent(inout) :: status
    character(len=:), allocatable :: t
    !---------------------------------------------------------------------------
    allocate(t, source=trim(adjustl(s)))
    t = lowercase(t)
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
  !> The one option setter.  Name and value both cross as strings and the parse
  !> happens here, against the option's own type, so no caller has to know which
  !> type bucket an option belongs to, and an unknown name stays distinguishable
  !> from a value that will not parse (see the status codes below).
  !>
  !> Options that influence injon (abundances, exclude_elements, ...) must still
  !> be set BEFORE hop_init.
  !>
  !> status: 0 ok, 1 unknown name, 4 value will not parse as that option's type
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
    call c2f_string(name, nm)
    nm = lowercase(nm)
    call c2f_string(value, val)
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
    case ('hydrogen_profiles')
      hydrogen_profiles = lowercase(val)
    case ('anstee_file')
      anstee_file = val
    case ('absmet_file')
      ! The sampled metal-line table (dispatch's &opac_table_params absmet_file), read
      ! at once.  osmet is process-wide and outlives hop_free, so later calls may
      ! repeat the loaded file but cannot replace it.
      if (osmet%initialised) then
        if (val /= absmet_loaded) call fail("absmet_file: '" // absmet_loaded // &
          "' is already loaded and cannot be replaced", 1, status)
      else
        inquire(file=val, exist=lv)
        if (lv) then
          call osmet%init(val)
          allocate(absmet_loaded, source=val)
        else
          call fail("absmet_file: '" // val // "' not found", 1, status)
        end if
      end if
    ! --- reals --------------------------------------------------------------
    case ('ref_lam')
      call parse_real(val, rv, status); if (status == 0) ref_lam = rv
    case ('debug_opac')
      call parse_real(val, rv, status); if (status == 0) debug_opac = rv
    case ('tolerance')
      call parse_real(val, rv, status); if (status == 0) tolerance = real(rv, 4)
    case ('line_reject_frac')
      call parse_real(val, rv, status); if (status == 0) line_reject_frac = real(rv, 4)
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
    case ('disable_rad')
      call parse_logical(val, lv, status); if (status == 0) disable_rad = lv
    case ('disable_stk')
      call parse_logical(val, lv, status); if (status == 0) disable_stk = lv
    case ('disable_vdw')
      call parse_logical(val, lv, status); if (status == 0) disable_vdw = lv
    case ('molecule_damping_estimates')
      call parse_logical(val, lv, status); if (status == 0) molecule_damping_estimates = lv
    case ('hydrogen_bpo')
      call parse_logical(val, lv, status); if (status == 0) hydrogen_BPO = lv
    case ('ll_profile')
      call parse_logical(val, lv, status); if (status == 0) ll_profile = lv
    ! --- integers -----------------------------------------------------------
    case ('voigt_kernel')
      ! parsed and stored, but nothing reads it: QLG-bot is the only kernel compiled in.
      ! Kept so that callers still setting it keep working. See linelist_mod's kernel note.
      call parse_int(val, iv, status); if (status == 0) voigt_kernel = iv
    case ('verbose')
      call parse_int(val, iv, status)
      if (status == 0) then
        init_verbose = iv
        call console_set_level (iv, ch_hop)
      end if
    case default
      call fail("unknown option '" // trim(nm) // "' (see hop_option_names)", 1, status)
    end select
  END FUNCTION hop_set_option

  !=============================================================================
  !> EOS
  !=============================================================================
  FUNCTION hop_load_state(n, temp, ne, rho, vmic, use_rho, use_ne) &
           RESULT(status) BIND(c, name='hop_load_state')
    integer(c_int), value :: n, use_rho, use_ne
    ! assumed-size: gen_bindings.py maps a dummy with any dimension to a pointer and a
    ! `value` scalar to a C scalar; everything in this interface stays one of the two
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
    call c2f_string(name, nm)
    nm = lowercase(nm)
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
  !> Mark a species ('Fe I', 'CH', ...) so its number density survives the solve.
  !> Call after hop_init and before hop_load_state; the mark is set on the master
  !> instance, so every subsequent solve retains the species.  Mirrors dispatch's
  !> save_species pre-init marking: in_opac is read only by the deallocate keep-
  !> checks, so this does not change any opacity.
  !-----------------------------------------------------------------------------
  FUNCTION hop_keep_species(name) RESULT(status) BIND(c, name='hop_keep_species')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int) :: status
    class(particle_type), pointer :: particle
    character(len=:), allocatable :: nm
    !---------------------------------------------------------------------------
    if (.not. allocated(master)) then
      call fail('hop_keep_species: no instance; call hop_init first', 2, status)
      return
    end if
    call c2f_string(name, nm)
    call species_lookup(master, nm, particle, status)
    if (status == 0) particle%in_opac = .true.
  END FUNCTION hop_keep_species

  !-----------------------------------------------------------------------------
  !> Number density [1/cm3] of a species in the loaded state.  The species must
  !> have been kept (hop_keep_species) before the hop_load_state call -- or be
  !> one HOP retains anyway (continuum-opacity and line-list species).
  !-----------------------------------------------------------------------------
  FUNCTION hop_get_species(name, n, out) RESULT(status) &
           BIND(c, name='hop_get_species')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: n
    real(c_double), intent(out) :: out(*)
    integer(c_int) :: status
    class(particle_type), pointer :: particle
    character(len=:), allocatable :: nm
    !---------------------------------------------------------------------------
    if (hop_get_n() /= n .or. n < 1) then
      call fail('hop_get_species: no state loaded, or n differs from the loaded size', &
                3, status)
      return
    end if
    call c2f_string(name, nm)
    call species_lookup(comp, nm, particle, status)
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
  !> Microturbulence [km/s] of the loaded state.  vmic enters only the line
  !> broadening (composition%add_linelist_opacities), never the EOS, so it can
  !> be changed on a solved state instead of re-solving it.  A negative value is
  !> thermal-only broadening.
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
  !> Opacity.  lam is vacuum Angstrom, ascending.  The chi and scat buffers are
  !> float32 of size nlam*n (Fortran order: wavelength fastest), filled in 1/cm
  !> (per-volume extinction; divide by rho for per-gram).  See
  !> composition_type%opacity for what each of them contains: chi is the total
  !> extinction (continuum + absmet + line list), scat its scattering part.
  !> line_scat < 0 leaves the split to the module-global line_scattering, 0/1
  !> override it; include_llist is honoured only once a line list has been read.
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
  !> Threaded Opacity
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
    real(8), allocatable :: lam8(:), llam(:)
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
    nch_l = min(max(1, n / MIN_CHUNK_CELLS), max(1, nthr * 64))
    nch_c = min(max(1, n / MIN_CHUNK_CELLS), max(1, nthr * 4))
#ifndef _OPENMP
    nch_l = 1                   
    nch_c = 1
#endif

    allocate(lam8(nlam))
    lam8 = lam(1:nlam)
    call sp%reset(lam=lam8)
    sp%sorted = .true.
    allocate(llam(nlam))
    llam = log(sp%lam)
    chi4 (1:nlam, 1:n, 1:1, 1:1) => chi (1:nlam*n)
    scat4(1:nlam, 1:n, 1:1, 1:1) => scat(1:nlam*n)

    m = shape(comp%temp)
    if (m(1) /= n) then
      call fail('hop_opacity_omp: state is not shaped (n,1,1)', 3, status)
      return
    end if
    nf = size(chi4, 1)
    dl = dlam
    do_absmet = (include_absmet /= 0) .and. osmet%initialised
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
          call comp%add_linelist_opacities(sp, chi4, ehat, dl, one_m_eps=one_m_eps, chi_scat=chi_scat, beta=bet, i0=i0, i1=i1, llam_in=llam)
        end do
        !$omp end parallel do
      else
        !$omp parallel do schedule(dynamic) default(shared) private(w, i0, i1)
        do w = 1, nch_l
          call chunk_range(w, nch_l, n, i0, i1)
          if (i1 < i0) cycle
          call comp%add_linelist_opacities(sp, chi4, ehat, dl, i0=i0, i1=i1, llam_in=llam)
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
                ! the line's (1-eps) share is already inside chi4 (the kernels add
                ! the whole line and only record the share in chi_scat), so chi4
                ! takes the continuum scattering alone while scat4 reports both --
                ! the split partitions the extinction, it does not add to it
                s = sig(f,i,j,k) + chi_scat(f,i,j,k)
                chi4 (f,i,j,k) = chi4(f,i,j,k) - ehat(f,i,j,k)*bet(f,i,j,k) &
                                               + sig(f,i,j,k)
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
    !=============================================================================
    !> Cell range [i0, i1] of chunk w of nch over n cells, balanced to +-1 cell.
    !> Empty (i1 < i0) when there are fewer cells than chunks.
    !=============================================================================
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

  !=============================================================================
  !> Read a line list keeping only the lines inside a set of windows, rather than
  !> everything between the first and last wavelength.  hop_init_linelist masks on
  !> [lam(1)-dlam, lam(nlam)+dlam], so a request for fifteen 2.5 A windows spread
  !> over 1300 A loads all 1300 A of it: measured on the master list, 18943 lines
  !> where the windows hold 484.  The windows are merged and expanded by dlam here
  !> exactly as a full run does it (build_linelist_intervals), so the kept set is the
  !> one that run would keep, wing contributors included.
  !>
  !> blue_vac/red_vac are vacuum Angstrom, like interval_t itself; the reader masks
  !> on the vacuum line centre.
  !=============================================================================
  FUNCTION hop_init_linelist_windows(linelist_file, dlam, nwin, blue_vac, red_vac) &
           RESULT(status) BIND(c, name='hop_init_linelist_windows')
    character(kind=c_char), intent(in) :: linelist_file(*)
    real(c_double), value :: dlam
    integer(c_int), value :: nwin
    real(c_double), intent(in) :: blue_vac(*), red_vac(*)
    integer(c_int) :: status
    type(interval_t), allocatable :: merged(:)
    real(8), allocatable :: windows(:,:), lam8(:)
    character(len=:), allocatable :: file_f
    integer :: i, nmerged
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(master)) then
      call fail('hop_init_linelist_windows: no instance; call hop_init first', 2, status)
      return
    end if
    if (.not. master%initialised) then
      call fail('hop_init_linelist_windows: instance not initialised', 2, status)
      return
    end if
    call c2f_string(linelist_file, file_f)
    if (len(file_f) == 0) then
      call fail('hop_init_linelist_windows: empty file name', 3, status)
      return
    end if
    if (nwin < 1) then
      call fail('hop_init_linelist_windows: no windows given', 3, status)
      return
    end if
    allocate(windows(nwin, 2))
    windows(:,1) = blue_vac(1:nwin)
    windows(:,2) = red_vac(1:nwin)
    call merge_intervals(windows, nwin, dlam, merged, nmerged)
    if (nmerged < 1) then
      call fail('hop_init_linelist_windows: the windows merged to nothing', 3, status)
      return
    end if
    ! the lam range is the hull of the merged intervals; the intervals do the masking
    allocate(lam8(2))
    lam8 = [minval(merged(1:nmerged)%blue_vac), maxval(merged(1:nmerged)%red_vac)]
    call master%init_linelist(file_f, dlam, lam8, merged(1:nmerged))
  END FUNCTION hop_init_linelist_windows

  !-----------------------------------------------------------------------------
  !> Read one line-list file into the master instance, keeping only lines inside
  !> [lam(1)-dlam, lam(nlam)+dlam] (pass nlam=0 for no wavelength mask).  Mirrors
  !> the composition%init_linelist calls dispatch makes at setup, one per file in
  !> linelist_folder; repeated calls accumulate.  Must be called after hop_init and
  !> BEFORE hop_load_state: the solve works on a fresh copy of master, and
  !> particle%ll is allocatable, so the copy carries the lines along.
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
    call c2f_string(linelist_file, file_f)
    if (len(file_f) == 0) then
      call fail('hop_init_linelist: empty file name', 3, status)
      return
    end if
    if (nlam > 0) then
      allocate(lam8(nlam))
      lam8 = lam(1:nlam)
      call master%init_linelist(file_f, dlam, lam8)
    else
      call master%init_linelist(file_f, dlam)
    end if
  END FUNCTION hop_init_linelist

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

  !=============================================================================
  ! Model-atmosphere files, parsed by atmos_reader_mod exactly as dispatch.x
  ! reads them, at the file's full resolution and in its own orientation.
  !=============================================================================
  !> Parse a file; format 'unknown' infers it from the name as dispatch.x does.
  !> snap is the Stagger snapshot number and use_aux /= 0 lets a Stagger snapshot
  !> take T and n_e from its .aux file (both ignored by the other formats).
  !> dims returns nx, ny, nz, the ghost zones, has_pe, has_rho, has_feh and the
  !> number of abundances the file carries (92 for MARCS, else 0).
  !=============================================================================
  FUNCTION hop_atmos_read(atmos_format, atmos_file, snap, use_aux, dims) RESULT(status) &
           BIND(c, name='hop_atmos_read')
    character(kind=c_char), intent(in) :: atmos_format(*), atmos_file(*)
    integer(c_int), value :: snap, use_aux
    integer(c_int), intent(out) :: dims(*)
    integer(c_int) :: status
    character(len=:), allocatable :: fmt, file
    character(len=512) :: msg
    integer :: ierr
    !---------------------------------------------------------------------------
    status = 0
    call c2f_string(atmos_format, fmt)
    fmt = lowercase(fmt)
    call c2f_string(atmos_file, file)
    if (fmt == 'unknown') then
      deallocate(fmt)
      allocate(fmt, source=trim(atmos_format_from_name(file)))
    end if
    if (allocated(atmos_file_data)) deallocate(atmos_file_data)
    allocate(atmos_file_data)
    call read_atmos_file(fmt, file, atmos_file_data, ierr, msg, snap=int(snap), &
                         use_aux=use_aux /= 0)
    if (ierr /= 0) then
      deallocate(atmos_file_data)
      call fail(trim(msg), int(ierr), status)
      return
    end if
    associate (a => atmos_file_data)
      dims(1:8) = [a%nx, a%ny, a%nz, a%gz, merge(1, 0, a%has_pe), merge(1, 0, a%has_rho), &
                   merge(1, 0, a%has_feh), 0]
      if (allocated(a%abund)) dims(8) = size(a%abund)
    end associate
  END FUNCTION hop_atmos_read

  !=============================================================================
  !> One field of the parsed atmosphere as float32, x fastest: ltemp, lpe, lrho
  !> (natural logs, cgs), vx, vy, vz [cm/s], vmic [km/s] (nx*ny*nz values), xx, yy,
  !> zz (the axes), feh (1) or abund (log eps for Z=1..92, MARCS only).
  !=============================================================================
  FUNCTION hop_atmos_get(name, n, out) RESULT(status) BIND(c, name='hop_atmos_get')
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: n
    real(c_float), intent(out) :: out(*)
    integer(c_int) :: status
    character(len=:), allocatable :: nm
    integer :: ncell
    !---------------------------------------------------------------------------
    status = 0
    if (.not. allocated(atmos_file_data)) then
      call fail('hop_atmos_get: no atmosphere read; call hop_atmos_read first', 2, status)
      return
    end if
    call c2f_string(name, nm)
    nm = lowercase(nm)
    associate (a => atmos_file_data)
      ncell = a%nx * a%ny * a%nz
      select case (nm)
      case ('ltemp', 'lpe', 'lrho', 'vx', 'vy', 'vz', 'vmic')
        if (n /= ncell) then
          call fail('hop_atmos_get: n differs from nx*ny*nz', 3, status)
          return
        end if
      case ('xx')
        if (n /= a%nx) call fail('hop_atmos_get: n differs from nx', 3, status)
      case ('yy')
        if (n /= a%ny) call fail('hop_atmos_get: n differs from ny', 3, status)
      case ('zz')
        if (n /= a%nz) call fail('hop_atmos_get: n differs from nz', 3, status)
      case ('feh')
        if (n /= 1) call fail('hop_atmos_get: feh is one value', 3, status)
      case ('abund')
        if (.not. allocated(a%abund)) then
          call fail('hop_atmos_get: the file carries no abundances', 3, status)
        else if (n /= size(a%abund)) then
          call fail('hop_atmos_get: n differs from the abundance count', 3, status)
        end if
      case default
        call fail("hop_atmos_get: unknown field '" // nm // "'", 1, status)
      end select
      if (status /= 0) return
      select case (nm)
      case ('ltemp'); out(1:n) = reshape(a%ltemp, [n])
      case ('lpe');   out(1:n) = reshape(a%lpe, [n])
      case ('lrho');  out(1:n) = reshape(a%lrho, [n])
      case ('vx');    out(1:n) = reshape(a%vx, [n])
      case ('vy');    out(1:n) = reshape(a%vy, [n])
      case ('vz');    out(1:n) = reshape(a%vz, [n])
      case ('vmic');  out(1:n) = reshape(a%vmic, [n])
      case ('xx');    out(1:n) = a%xx
      case ('yy');    out(1:n) = a%yy
      case ('zz');    out(1:n) = a%zz
      case ('feh');   out(1) = real(a%feh, c_float)
      case ('abund'); out(1:n) = a%abund
      end select
    end associate
  END FUNCTION hop_atmos_get

  FUNCTION hop_atmos_free() RESULT(status) BIND(c, name='hop_atmos_free')
    integer(c_int) :: status
    !---------------------------------------------------------------------------
    status = 0
    if (allocated(atmos_file_data)) deallocate(atmos_file_data)
  END FUNCTION hop_atmos_free

  !=============================================================================
  !> Optical depth of ncol columns of nz cells, each top first and depth fastest
  !> (C order (ncol, nz)), with the integrator and top rule of dispatch.x
  !> (interpol_mod optical_depth): chi [1/cm], ln rho and the height z [cm].
  !> Needs no hop_init.
  !=============================================================================
  FUNCTION hop_optical_depth(ncol, nz, chi, lnrho, z, tau) RESULT(status) &
           BIND(c, name='hop_optical_depth')
    integer(c_int), value :: ncol, nz
    real(c_double), intent(in) :: chi(*), lnrho(*), z(*)
    real(c_double), intent(out) :: tau(*)
    integer(c_int) :: status
    integer :: j, i0
    !---------------------------------------------------------------------------
    status = 0
    if (nz < 2) then
      call fail('hop_optical_depth: a column needs two cells or more', 3, status)
      return
    end if
    do j=1, ncol
      i0 = (j-1) * nz
      call optical_depth(chi(i0+1:i0+nz), lnrho(i0+1:i0+nz), z(i0+1:i0+nz), &
                         tau(i0+1:i0+nz))
    end do
  END FUNCTION hop_optical_depth

  !=============================================================================
  !> Cumulative PCHIP integrals of y over x from the first point, for ncol series
  !> of n points, each contiguous (C order (ncol, n)) -- the quadrature of
  !> optical_depth (interpol_mod pchip_integral).  Needs no hop_init.
  !=============================================================================
  FUNCTION hop_pchip_integral(ncol, n, x, y, cum) RESULT(status) &
           BIND(c, name='hop_pchip_integral')
    integer(c_int), value :: ncol, n
    real(c_double), intent(in) :: x(*), y(*)
    real(c_double), intent(out) :: cum(*)
    integer(c_int) :: status
    integer :: j, i0
    !---------------------------------------------------------------------------
    status = 0
    if (n < 2) then
      call fail('hop_pchip_integral: a series needs two points or more', 3, status)
      return
    end if
    do j=1, ncol
      i0 = (j-1) * n
      call pchip_integral(x(i0+1:i0+n), y(i0+1:i0+n), cum(i0+1:i0+n))
    end do
  END FUNCTION hop_pchip_integral

END MODULE hop_c_api
