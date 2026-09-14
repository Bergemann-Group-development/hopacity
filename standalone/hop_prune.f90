!===============================================================================
!> hop_prune.x input.nml -- select the transitions a synthesis on one atmosphere
!> needs, and write them as a line list in the ordinary input format.
!>
!> The program reads the same namelist groups its consuming dispatch.x run reads
!> (&atmos_params, &composition_params, &opac_table_params, &linelist_params,
!> &spectrum_params) plus a small &prune_params group, so one set of settings
!> describes both.  It parses the atmosphere with atmos_reader_mod, samples its
!> depths, solves the EOS on those states through the standalone HOP C API, loads
!> the configured line lists over the requested windows expanded by dlam, computes
!> the opacity over the run's own wavelength grid so HOP marks the transitions that
!> pass its negligibility test, and writes those with hop_write_linelist.
!>
!> A run that reads the written file as its linelist_file selects the same
!> transitions it would have kept from the master list, so a ladder of runs on one
!> atmosphere (a curve of growth, an abundance fit) reads a few hundred lines each
!> instead of re-reading and re-selecting the whole master list.  Prune at the
!> strongest abundance the ladder will use: strength grows with abundance, so the
!> top rung's selection covers every weaker one.
!>
!> The output is written to a process-specific temporary file and renamed into
!> place, so an interrupted or failed run never leaves a truncated line list behind.
!> A <output_file>.provenance.nml sidecar echoes every input group the selection
!> depended on, plus a &prune_result group with the derived sampling and windows.  Invalid input and failures inside HOP exit with status 1 and a
!> one-line diagnostic on stderr; HOP's own error stops also terminate the process
!> with a nonzero status.
!===============================================================================
PROGRAM hop_prune
  USE iso_c_binding, only: c_char, c_null_char, c_int
  USE iso_fortran_env, only: error_unit, output_unit
  USE hop_c_api, only: hop_init, hop_free, hop_set_option, hop_set_external_abundances, &
                       hop_init_linelist_windows, hop_load_state, hop_set_vmic, &
                       hop_opacity, hop_write_linelist, hop_linelist_count, hop_last_error
  USE hop_env_mod, only: io
  USE atmos_reader_mod, only: atmos_data_t, read_atmos_file, atmos_format_from_name
  USE spectrum_mod, only: spectrum_type, daa, aa_red, aa_blue, nnu, give_lam, in_log, R, &
                          lam_file, window_file, in_air
  USE composition_mod, only: abundances, relative_abundances, exclude_elements, &
                             only_elements, opac_table_vmic, mhd_eos, switch, &
                             line_scattering, helium_scattering, ref_lam, debug_opac
  USE bbtrans_mod, only: tolerance, disable_ga, disable_stark, disable_vdw, ga_fallback
  USE bftrans_mod, only: H_occupation_probability
  USE ion_mod, only: debye_shield
  USE molecule_mod, only: irwin_dis
  USE hprofl_mod, only: hydrogen_BPO
  USE hvcs_mod, only: hydrogen_vcs_table, vcs_data_dir
  USE linelist_mod, only: ll_profile, voigt_kernel
  USE osmet_mod, only: eosmet
  USE const_mod, only: bk
  USE handyfuncs_mod, only: lowercase
  implicit none

  ! --- &atmos_params, the read_external_mod list; only the named few are used ----
  character(len=512) :: atmos_file=''
  character(len=64) :: atmos_format='unknown', amr=''
  integer :: dims=1, verbose=0, patches_per_mpi=1, nx=-1, ny=-1, nz=96, snap=0
  real(8) :: FeH=-99.9
  real(4) :: vmic=-1, timestamp=-1.0, prune_tol=0.0
  logical :: save_atmos=.true., use_rho=.true., use_ne=.true., pchip_interpolation=.true., &
             log_interpolation=.true., save_eos=.false., use_abnd=.false., use_vel=.true., &
             read_departures=.false., prune_top=.false.
  namelist /atmos_params/ atmos_file, dims, verbose, patches_per_mpi, FeH, &
    save_atmos, atmos_format, vmic, use_rho, use_ne, nx, ny, nz, snap, pchip_interpolation, &
    log_interpolation, save_eos, amr, use_abnd, timestamp, use_vel, read_departures, &
    prune_top, prune_tol

  ! --- &composition_params, the composition_patch_mod list (group in read_composition,
  !     whose 'verbose' would otherwise clash with the atmos_params one) -------------
  character(len=256) :: absdat_file='./input/TS_absdat.dat'
  character(len=256) :: abund_file='./input/abund/abund_magg22'
  character(len=256) :: isotope_file='./input/isotopes_blp'
  character(len=256) :: molec_file='./input/molec_data_BC.txt'
  character(len=256) :: save_species=''
  real(4) :: tmolim=-99, alpha_enhancement=0.0, r_process_enhancement=0.0
  integer :: comp_verbose=0

  ! --- &opac_table_params: only opac_table_vmic influences the selection ------------
  character(len=256) :: opac_table='', absmet_file=''
  logical :: pca_opac=.true.
  real(8) :: dtemp=250., ldne=0.25, ldrho=0.25, ldtemp=0.05
  namelist /opac_table_params/ opac_table, opac_table_vmic, pca_opac, &
    dtemp, ldne, ldrho, ldtemp, absmet_file

  ! --- &linelist_params ---------------------------------------------------------------
  real(8) :: dlam=100.
  character(len=512) :: linelist_folder=''
  character(len=256) :: linelist_file(128)
  namelist /linelist_params/ dlam, exclude_elements, linelist_file, linelist_folder, &
    tolerance, only_elements, line_scattering, ll_profile, voigt_kernel

  ! --- &prune_params ------------------------------------------------------------------
  ! sample_count depths spread evenly over the model's index range (0 = every depth).
  ! Spanning the whole range is what matters: a strength is a maximum over the
  ! sampled states and the maximum of a high-excitation line sits deep, that of a
  ! resonance or molecular line at the top.
  character(len=512) :: output_file=''
  integer :: sample_count=8
  namelist /prune_params/ output_file, sample_count

  type(atmos_data_t) :: atmos
  type(spectrum_type) :: sp
  character(len=512) :: nml_file, msg, tmp_file
  character(len=1024) :: abundances_in
  character(len=16) :: fmt
  integer :: nml_unit, echo_unit, iostat, ierr, i, nfiles, nsample, nwin, nkept
  integer, allocatable :: sample(:)
  real(8), allocatable :: temp(:), ne(:), rho(:), blue(:), red(:)
  real(4), allocatable :: chi(:), scat(:)
  real(4) :: vmic_eff
  logical :: with_rho, with_ne
  !-----------------------------------------------------------------------------
  linelist_file = ''
  abundances = ''                  ! module global without an initializer
  if (command_argument_count() /= 1) call die('usage: hop_prune.x input.nml')
  call get_command_argument(1, nml_file)
  open(newunit=nml_unit, file=trim(nml_file), status='old', action='read', iostat=iostat)
  if (iostat /= 0) call die('cannot open ' // trim(nml_file))
  ! HOP's own namelist readers (spectrum_type%init) read io%input and echo the
  ! group to io%namelist; the echo goes to a scratch unit here.
  io%input = nml_unit
  open(newunit=echo_unit, status='scratch', form='formatted')
  io%namelist = echo_unit

  ! --- the namelist groups; the atmosphere and the output are mandatory --------------
  call read_group('atmos_params', required=.true.)
  call read_group('composition_params')
  call read_group('opac_table_params')
  call read_group('linelist_params')
  call read_group('prune_params', required=.true.)
  abundances_in = abundances       ! hop_init resets the global unless set as an option
  if (output_file == '') call die('prune_params: output_file is required')
  if (atmos_file == '') call die('atmos_params: atmos_file is required')
  if (linelist_folder /= '') &
    call die('linelist_params: linelist_folder is not supported; list the files in linelist_file')
  nfiles = count(linelist_file /= '')
  if (nfiles == 0) call die('linelist_params: linelist_file is required')
  do i = 1, nfiles
    call check_exists(linelist_file(i))
  end do

  ! --- the wavelength grid and its windows, as the run builds them --------------------
  call sp%init()
  if (sp%nnu < 1) call die('spectrum_params: no wavelength grid; give R with ' // &
                           'window_file or aa_blue/aa_red')
  if (sp%n_segments > 0) then
    nwin = sp%n_segments
    allocate(blue(nwin), red(nwin))
    blue = sp%segment_blue_vac
    red  = sp%segment_red_vac
  else
    nwin = 1
    allocate(blue(1), red(1))
    blue = minval(sp%lam)
    red  = maxval(sp%lam)
  end if

  ! --- the atmosphere ---------------------------------------------------------------
  fmt = lowercase(trim(atmos_format))
  if (fmt == 'unknown') then
    fmt = atmos_format_from_name(atmos_file)
    if (fmt == '') call die('cannot infer atmos_format from ' // trim(atmos_file) // &
                            '; set it in atmos_params')
  end if
  call read_atmos_file(fmt, atmos_file, atmos, ierr, msg)
  if (ierr /= 0) call die(trim(msg))
  if (atmos%nx * atmos%ny > 1) &
    call die('hop_prune.x samples 1D atmospheres only; no sampling policy is ' // &
             'defined for 3D cubes yet')
  ! the namelist overrides read_atmos_mod applies before the EOS sees the model
  if (FeH < -99 .and. atmos%has_feh) FeH = atmos%feh
  if (vmic >= 0) atmos%vmic = vmic
  with_rho = use_rho .and. atmos%has_rho
  with_ne  = use_ne .and. atmos%has_pe
  if (.not. (with_rho .or. with_ne)) &
    call die('the atmosphere provides neither rho nor Pe under the use_rho/use_ne settings')
  ! read_external_mod: the table vmic is the atmosphere's mean when any of it is set
  vmic_eff = opac_table_vmic
  if (any(atmos%vmic > 0)) vmic_eff = sum(atmos%vmic) / size(atmos%vmic)

  ! --- sampled depths -----------------------------------------------------------------
  call sample_depths(atmos%nz, sample_count, sample)
  nsample = size(sample)
  allocate(temp(nsample), ne(nsample), rho(nsample))
  ne = 0; rho = 0
  ! the same conversions read_atmos_mod applies to the 1D state
  temp = exp(dble(atmos%ltemp(1,1,sample)))
  if (with_rho) rho = exp(dble(atmos%lrho(1,1,sample)))
  if (with_ne)  ne  = exp(dble(atmos%lpe(1,1,sample))) / (bk * temp)

  ! --- HOP: composition, line lists, state, opacity -----------------------------------
  call check_status(hop_set_option(cstr('verbose'), cstr(itoa(comp_verbose))), 'verbose')
  call check_status(hop_set_option(cstr('abundances'), cstr(abundances_in)), 'abundances')
  if (use_abnd .and. allocated(atmos%abund)) &
    call check_status(hop_set_external_abundances(int(size(atmos%abund), c_int), &
                      dble(atmos%abund)), 'hop_set_external_abundances')
  call check_exists(absdat_file)
  call check_exists(abund_file)
  call check_status(hop_init(cstr(absdat_file), cstr(abund_file), cstr(isotope_file), &
                    cstr(molec_file), FeH, tmolim, alpha_enhancement, r_process_enhancement), &
                    'hop_init')
  do i = 1, nfiles
    call check_status(hop_init_linelist_windows(cstr(linelist_file(i)), dlam, &
                      int(nwin, c_int), blue, red), 'hop_init_linelist_windows')
  end do
  call check_status(hop_load_state(int(nsample, c_int), temp, ne, rho, vmic_eff, &
                    merge(1_c_int, 0_c_int, with_rho), merge(1_c_int, 0_c_int, with_ne)), &
                    'hop_load_state')
  call check_status(hop_set_vmic(int(nsample, c_int), vmic_eff), 'hop_set_vmic')
  allocate(chi(sp%nnu * nsample), scat(sp%nnu * nsample))
  call check_status(hop_opacity(int(sp%nnu, c_int), sp%lam, chi, scat, 0_c_int, 1_c_int, &
                    dlam, -1_c_int), 'hop_opacity')
  nkept = hop_linelist_count()

  ! --- output: temporary file, then an atomic rename ----------------------------------
  write(tmp_file, '(A,A,I0,A)') trim(output_file), '.', getpid(), '.partial'
  call check_status(hop_write_linelist(cstr(tmp_file)), 'hop_write_linelist')
  call rename(trim(tmp_file), trim(output_file), iostat)
  if (iostat /= 0) call die('cannot rename ' // trim(tmp_file) // ' to ' // trim(output_file))
  call write_provenance(trim(output_file) // '.provenance.nml')
  if (verbose >= 0) write(output_unit, '(A,I0,A,A)') 'hop_prune: ', nkept, &
    ' transitions written to ', trim(output_file)
  iostat = hop_free()

CONTAINS

  !-----------------------------------------------------------------------------
  !> Read one namelist group from the input; a missing optional group keeps the
  !> defaults, a malformed one is an error.
  !-----------------------------------------------------------------------------
  SUBROUTINE read_group(name, required)
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: required
    integer :: ios
    !---------------------------------------------------------------------------
    rewind(nml_unit)
    select case (name)
    case ('atmos_params');        read(nml_unit, nml=atmos_params, iostat=ios)
    case ('composition_params');  call read_composition(ios)
    case ('opac_table_params');   read(nml_unit, nml=opac_table_params, iostat=ios)
    case ('linelist_params');     read(nml_unit, nml=linelist_params, iostat=ios)
    case ('prune_params');        read(nml_unit, nml=prune_params, iostat=ios)
    end select
    if (ios > 0) call die('cannot read &' // name // ' from ' // trim(nml_file))
    if (ios < 0 .and. present(required)) then
      if (required) call die('&' // name // ' is missing from ' // trim(nml_file))
    end if
  END SUBROUTINE read_group

  SUBROUTINE read_composition(ios, echo_unit)
    integer, intent(out), optional :: ios
    integer, intent(in), optional :: echo_unit
    integer :: verbose
    namelist /composition_params/ verbose, absdat_file, abund_file, &
      isotope_file, molec_file, H_occupation_probability, switch, mhd_eos, debye_shield, &
      tmolim, irwin_dis, save_species, ref_lam, debug_opac, abundances, relative_abundances, &
      alpha_enhancement, r_process_enhancement, disable_ga, disable_stark, disable_vdw, &
      helium_scattering, hydrogen_BPO, hydrogen_vcs_table, vcs_data_dir, ga_fallback, eosmet
    !---------------------------------------------------------------------------
    verbose = comp_verbose
    if (present(echo_unit)) then
      write(echo_unit, nml=composition_params)
      return
    end if
    read(nml_unit, nml=composition_params, iostat=ios)
    comp_verbose = verbose
  END SUBROUTINE read_composition

  !-----------------------------------------------------------------------------
  !> n depth indices spread evenly over 1..nz (all of them for n <= 0 or n >= nz),
  !> the sampling tumult.linestrength.read_marcs established: rounded linspace
  !> over the full index range, duplicates removed.
  !-----------------------------------------------------------------------------
  SUBROUTINE sample_depths(nz, n, idx)
    integer, intent(in) :: nz, n
    integer, allocatable, intent(out) :: idx(:)
    integer :: m, i, k, last
    !---------------------------------------------------------------------------
    m = n
    if (m <= 0 .or. m > nz) m = nz
    allocate(idx(m))
    k = 0
    last = 0
    do i = 0, m-1
      if (m == 1) then
        idx(1) = 1
        k = 1
        exit
      end if
      if (nint(dble(nz-1) * dble(i) / dble(m-1)) + 1 /= last) then
        k = k + 1
        idx(k) = nint(dble(nz-1) * dble(i) / dble(m-1)) + 1
        last = idx(k)
      end if
    end do
    idx = idx(1:k)
  END SUBROUTINE sample_depths

  !-----------------------------------------------------------------------------
  !> The inputs as read, echoed group by group the way a run echoes its own, plus a
  !> &prune_result group with what the program derived from them: the format, the
  !> depths it sampled, the vmic and [Fe/H] the EOS saw, the vacuum windows and the
  !> number of transitions written.
  !-----------------------------------------------------------------------------
  SUBROUTINE write_provenance(path)
    character(len=*), intent(in) :: path
    character(len=16) :: atmos_format_used
    integer :: u, nz_model, nsampled, nnu_grid, ntransitions
    integer, allocatable :: sampled_depths(:)
    real(8), allocatable :: blue_vac(:), red_vac(:)
    real(4) :: vmic_used
    real(8) :: feh_used
    logical :: use_rho_used, use_ne_used, use_abnd_used
    ! spectrum_type%init reads this group from io%input and echoes it to io%namelist,
    ! which is a scratch unit here, so it is echoed again from the same module variables
    namelist /spectrum_params/ daa, aa_red, aa_blue, nnu, give_lam, in_log, R, &
      lam_file, window_file, in_air
    namelist /prune_result/ atmos_format_used, nz_model, nsampled, sampled_depths, &
      use_rho_used, use_ne_used, use_abnd_used, vmic_used, feh_used, blue_vac, red_vac, &
      nnu_grid, ntransitions
    !---------------------------------------------------------------------------
    atmos_format_used = fmt
    nz_model = atmos%nz
    nsampled = nsample
    sampled_depths = sample
    use_rho_used = with_rho
    use_ne_used = with_ne
    use_abnd_used = use_abnd .and. allocated(atmos%abund)
    vmic_used = vmic_eff
    feh_used = FeH
    blue_vac = blue
    red_vac = red
    nnu_grid = sp%nnu
    ntransitions = nkept
    open(newunit=u, file=path, status='replace', action='write')
    write(u, nml=atmos_params)
    call read_composition(echo_unit=u)
    write(u, nml=opac_table_params)
    write(u, nml=linelist_params)
    write(u, nml=spectrum_params)
    write(u, nml=prune_params)
    write(u, nml=prune_result)
    close(u)
  END SUBROUTINE write_provenance

  !-----------------------------------------------------------------------------
  FUNCTION cstr(s) RESULT(c)
    character(len=*), intent(in) :: s
    character(kind=c_char, len=len_trim(s)+1) :: c
    !---------------------------------------------------------------------------
    c = trim(s) // c_null_char
  END FUNCTION cstr

  FUNCTION itoa(i) RESULT(s)
    integer, intent(in) :: i
    character(len=16) :: s
    !---------------------------------------------------------------------------
    write(s, '(I0)') i
  END FUNCTION itoa

  SUBROUTINE check_exists(path)
    character(len=*), intent(in) :: path
    logical :: exists
    !---------------------------------------------------------------------------
    inquire(file=trim(path), exist=exists)
    if (.not. exists) call die('file not found: ' // trim(path))
  END SUBROUTINE check_exists

  SUBROUTINE check_status(status, what)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: what
    character(kind=c_char) :: buf(512)
    character(len=512) :: text
    integer :: i, ios
    !---------------------------------------------------------------------------
    if (status == 0) return
    ios = hop_last_error(buf, 512_c_int)
    text = ''
    do i = 1, 512
      if (buf(i) == c_null_char) exit
      text(i:i) = buf(i)
    end do
    call die(what // ' failed: ' // trim(text))
  END SUBROUTINE check_status

  SUBROUTINE die(text)
    character(len=*), intent(in) :: text
    !---------------------------------------------------------------------------
    write(error_unit, '(A)') 'hop_prune: ' // text
    call exit(1)
  END SUBROUTINE die

END PROGRAM hop_prune
