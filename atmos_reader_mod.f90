!===============================================================================
!> Model-atmosphere file readers with no DISPATCH dependency.
!>
!> Both consumers of an atmosphere file share these routines: read_atmos_mod copies
!> the parsed data into the runtime input_atmos_t (mesh, patches, AMR, pruning of the
!> top planes), and hop_prune.x reduces it to the thermodynamic states it hands to
!> standalone HOP.  Everything a format-specific parser knows lands in atmos_data_t;
!> everything a consumer decides (namelist overrides of vmic or [Fe/H], which of
!> rho/pe the EOS uses, messages) stays with the consumer.
!>
!> The thermodynamic columns are stored as natural logarithms in the same real(4)
!> precision read_atmos_mod has always used, so the runtime copy is exact and the
!> binary formats that tabulate logarithms convert with one multiplication.
!>
!> Failures are reported through ierr/msg rather than stop, so a standalone caller
!> can turn them into an exit status; read_atmos_mod stops on them itself.
!===============================================================================
MODULE atmos_reader_mod
  USE const_mod, only: bk
  USE handyfuncs_mod, only: nextentry, count_words, lowercase
  implicit none
  private

  type, public :: atmos_data_t
    character(len=16) :: format = ''
    integer :: nx=0, ny=0, nz=0
    ! ln T [K], ln P_e [dyn/cm2], ln rho [g/cm3]; a column a format does not carry
    ! is flagged absent and left unset.
    real(4), allocatable :: ltemp(:,:,:), lpe(:,:,:), lrho(:,:,:)
    logical :: has_pe=.false., has_rho=.false.
    real(4), allocatable :: vmic(:,:,:)      ! [km/s], the file's own values
    real(4), allocatable :: vz(:,:,:)        ! [cm/s] where the format has it, else 0
    real(4), allocatable :: zz(:)            ! depth coordinate as the file gives it
    logical :: has_feh=.false.
    real(8) :: feh=0                         ! model [Fe/H] when the file states one
    real(4) :: radius=0                      ! [cm], MARCS (1 for plane-parallel)
    real(4), allocatable :: abund(:)         ! MARCS: log eps for Z=1..92
  contains
    procedure :: allocate => allocate_atmos_data
  end type atmos_data_t

  public :: read_atmos_file, read_marcs_file, read_text_1d_file, read_multi_1d_file
  public :: atmos_format_from_name

CONTAINS

  !=============================================================================
  SUBROUTINE allocate_atmos_data(self, nx, ny, nz)
    class(atmos_data_t), intent(inout) :: self
    integer, intent(in) :: nx, ny, nz
    !---------------------------------------------------------------------------
    self%nx = nx; self%ny = ny; self%nz = nz
    allocate(self%zz(nz), source=0.)
    allocate(self%vz(nx,ny,nz), source=0.)
    allocate(self%vmic(nx,ny,nz), source=0.)
    allocate(self%ltemp(nx,ny,nz), self%lpe(nx,ny,nz), self%lrho(nx,ny,nz))
  END SUBROUTINE allocate_atmos_data

  !=============================================================================
  !> The format read_external_mod infers from a file name when atmos_format is
  !> 'unknown'; '' when nothing matches.
  !=============================================================================
  FUNCTION atmos_format_from_name(atmos_file) RESULT(fmt)
    character(len=*), intent(in) :: atmos_file
    character(len=16) :: fmt
    !---------------------------------------------------------------------------
    if (index(atmos_file, '.mod') > 0) then
      fmt = 'marcs'
    elseif (index(atmos_file, 'atmos.') > 0) then
      fmt = 'multi1d'
    elseif (index(atmos_file, '.dat') > 0) then
      fmt = 'siu_bin'
    elseif (index(atmos_file, '.txt') > 0) then
      fmt = 'text'
    elseif (index(lowercase(atmos_file), 'stagger_grid2') > 0) then
      fmt = 'stagger2'
    elseif (index(lowercase(atmos_file), 'stagger') > 0) then
      fmt = 'stagger'
    elseif (index(atmos_file, '.atmos') > 0) then
      fmt = 'multi'
    elseif (index(atmos_file, 'data/') > 0) then
      fmt = 'save'
    else
      fmt = ''
    endif
  END FUNCTION atmos_format_from_name

  !=============================================================================
  !> Dispatch on the format names read_external_mod accepts.  Formats this module
  !> does not parse yet (the 3D snapshots, SIU) return ierr=2.
  !=============================================================================
  SUBROUTINE read_atmos_file(atmos_format, atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_format, atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    !---------------------------------------------------------------------------
    select case (lowercase(trim(atmos_format)))
    case ('marcs')
      call read_marcs_file(atmos_file, atmos, ierr, msg)
    case ('text')
      call read_text_1d_file(atmos_file, atmos, ierr, msg)
    case ('multi1d')
      call read_multi_1d_file(atmos_file, atmos, ierr, msg)
    case default
      ierr = 2
      msg = "atmos_format '" // trim(atmos_format) // "' is not supported by the " // &
            "standalone reader (marcs, text, multi1d)"
    end select
  END SUBROUTINE read_atmos_file

  !=============================================================================
  SUBROUTINE open_existing(atmos_file, lunit, ierr, msg, form, access)
    character(len=*), intent(in) :: atmos_file
    integer, intent(out) :: lunit, ierr
    character(len=*), intent(out) :: msg
    character(len=*), intent(in), optional :: form, access
    logical :: exists
    integer :: iostat
    !---------------------------------------------------------------------------
    ierr = 0
    msg = ''
    inquire(file=atmos_file, exist=exists)
    if (.not. exists) then
      ierr = 1
      msg = 'File not found - ' // trim(atmos_file)
      return
    endif
    if (present(form)) then
      open(newunit=lunit, file=atmos_file, form=form, access=access, iostat=iostat)
    else
      open(newunit=lunit, file=atmos_file, iostat=iostat)
    endif
    if (iostat /= 0) then
      ierr = 1
      msg = 'Cannot open - ' // trim(atmos_file)
    endif
  END SUBROUTINE open_existing

  !=============================================================================
  !> MARCS .mod: header (Teff, log g, vmic, [Fe/H], radius, 92 abundances), then
  !> the T/Pe table and the density table, both keyed by the depth index.
  !=============================================================================
  SUBROUTINE read_marcs_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    integer :: lunit, nz, kk, k, iostat
    real(8) :: FeH_marcs
    real(4) :: teff, flux, vmic_marcs, g, rmass, alpha, rluminosity
    real(4) :: xx, yy, zz, palfa, pny, py, pbeta
    real(4) :: lgTauR, lgTau5, pg, prad, pturb, xkapr, emu, vconv, fconv, rhox
    real(4), dimension(:), allocatable :: temp, pe, rho
    character(512) :: fileID
    !---------------------------------------------------------------------------
    call open_existing(atmos_file, lunit, ierr, msg)
    if (ierr /= 0) return
    atmos%format = 'marcs'
    allocate(atmos%abund(92))
    read(lunit,*, iostat=iostat) fileID
    read(lunit,*, iostat=iostat) teff                  ! Teff [K].
    read(lunit,*, iostat=iostat) flux                  ! Flux [erg/cm2/s]
    read(lunit,*, iostat=iostat) g                     ! Surface gravity [cm/s2]
    read(lunit,*, iostat=iostat) vmic_marcs            ! Microturbulence parameter [km/s]
    read(lunit,*, iostat=iostat) rmass                 ! No mass for plane-parallel models
    read(lunit,*, iostat=iostat) FeH_marcs, alpha      ! Metallicity [Fe/H] and [alpha/Fe]
    read(lunit,*, iostat=iostat) atmos%radius          ! radius [cm] (plane-parallel: 1)
    read(lunit,*, iostat=iostat) rluminosity           ! Luminosity [Lsun] for a radius of 1 cm
    read(lunit,*, iostat=iostat) palfa, pny, py, pbeta ! Convection parameters
    read(lunit,*, iostat=iostat) xx, yy, zz            ! X, Y and Z
    read(lunit,*, iostat=iostat)                       ! skip comment
    read(lunit,*, iostat=iostat) atmos%abund
    read(lunit,*, iostat=iostat) nz
    read(lunit,*, iostat=iostat)                       ! skip comment
    read(lunit,*, iostat=iostat)                       ! skip comment
    if (iostat /= 0 .or. nz < 1) then
      close(lunit)
      ierr = 3
      msg = 'Cannot parse the MARCS header of ' // trim(atmos_file)
      return
    endif

    call atmos%allocate(1,1,nz)
    allocate(temp(nz), pe(nz), rho(nz))
    atmos%has_feh = .true.
    atmos%feh = FeH_marcs
    atmos%vmic(1,1,:) = vmic_marcs

    do kk=1, nz
      read(lunit,*, iostat=iostat) k, lgTauR, lgTau5, atmos%zz(k), temp(k), pe(k), pg, prad, pturb
      if (iostat /= 0 .or. kk /= k) then
        close(lunit)
        ierr = 3
        msg = 'model and loop indices differ'
        return
      endif
    end do

    read(lunit,*, iostat=iostat)                       ! skip comment

    do kk=1, nz
      read(lunit,*, iostat=iostat) k, lgTauR, xkapr, rho(k), emu, vconv, fconv, rhox
      if (iostat /= 0 .or. kk /= k) then
        close(lunit)
        ierr = 3
        msg = 'model and loop indices differ'
        return
      endif
    end do
    close(lunit)

    atmos%ltemp(1,1,:) = log(temp)
    atmos%lpe  (1,1,:) = log(pe)
    atmos%lrho (1,1,:) = log(rho)
    atmos%has_pe = .true.
    atmos%has_rho = .true.
  END SUBROUTINE read_marcs_file

  !=============================================================================
  !> Text 1D: name, "nz [FeH]", then one "z T Pe rho vmic" row per depth.  A
  !> column of zeros marks Pe or rho as absent, and the logarithms are taken only
  !> where the value is positive, as read_atmos_mod always has.
  !=============================================================================
  SUBROUTINE read_text_1d_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    integer :: nwords, nz, lunit, i, iostat
    character(128) :: text
    real(4), dimension(:), allocatable :: temp, pe, rho
    !---------------------------------------------------------------------------
    call open_existing(atmos_file, lunit, ierr, msg)
    if (ierr /= 0) return
    atmos%format = 'text'
    call nextentry(lunit)
    read(lunit, *) ! skip name
    call nextentry(lunit)

    read(lunit, '(A)') text
    call count_words(text, nwords)
    if (nwords > 1) then
      read(text, *, iostat=iostat) nz, atmos%feh
      atmos%has_feh = .true.
    else
      read(text, *, iostat=iostat) nz
    endif
    if (iostat /= 0 .or. nz < 1) then
      close(lunit)
      ierr = 3
      msg = 'Cannot parse the depth count of ' // trim(atmos_file)
      return
    endif

    call atmos%allocate(1,1,nz)
    allocate(temp(nz), pe(nz), rho(nz))

    do i=1, nz
      call nextentry(lunit)
      read(lunit, *, iostat=iostat) atmos%zz(i), temp(i), pe(i), rho(i), atmos%vmic(1,1,i)
      if (iostat /= 0) then
        close(lunit)
        ierr = 3
        write(msg, '(A,I0,A)') 'Cannot parse depth ', i, ' of ' // trim(atmos_file)
        return
      endif
    end do
    close(lunit)

    atmos%has_pe = pe(1) /= 0
    atmos%has_rho = rho(1) /= 0
    where (temp > 0) atmos%ltemp(1,1,:) = log(temp)
    where (rho  > 0) atmos%lrho (1,1,:) = log(rho)
    where (pe   > 0) atmos%lpe  (1,1,:) = log(pe)
  END SUBROUTINE read_text_1d_file

  !=============================================================================
  !> MULTI 1D atmos.* : name, scale name, log g, nz, then "z T ne vz vmic" rows.
  !> The file carries electron density, stored here as ln P_e = ln(n_e k T); there
  !> is no mass density.
  !=============================================================================
  SUBROUTINE read_multi_1d_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    integer :: nz, lunit, i, iostat
    character(128) :: text
    real(4), dimension(:), allocatable :: temp, ne
    !---------------------------------------------------------------------------
    call open_existing(atmos_file, lunit, ierr, msg)
    if (ierr /= 0) return
    atmos%format = 'multi1d'
    call nextentry(lunit)
    read(lunit, *) ! skip name
    call nextentry(lunit)

    read(lunit, '(A)') text ! scale name
    text = trim(text)
    if (lowercase(text(1:1))=='m') then
      close(lunit)
      ierr = 3
      msg = 'Mass scale not supported in multi 1D atm format'
      return
    endif
    call nextentry(lunit)
    read(lunit, '(A)') text ! logg
    call nextentry(lunit)
    read(lunit, *, iostat=iostat) nz
    if (iostat /= 0 .or. nz < 1) then
      close(lunit)
      ierr = 3
      msg = 'Cannot parse the depth count of ' // trim(atmos_file)
      return
    endif

    call atmos%allocate(1,1,nz)
    allocate(temp(nz), ne(nz))

    do i=1, nz
      call nextentry(lunit)
      read(lunit, *, iostat=iostat) atmos%zz(i), temp(i), ne(i), atmos%vz(1,1,i), atmos%vmic(1,1,i)
      if (iostat /= 0) then
        close(lunit)
        ierr = 3
        write(msg, '(A,I0,A)') 'Cannot parse depth ', i, ' of ' // trim(atmos_file)
        return
      endif
    end do
    close(lunit)

    atmos%ltemp(1,1,:) = log(temp)
    atmos%lpe  (1,1,:) = log(ne * BK * temp)
    atmos%has_pe = .true.
    atmos%has_rho = .false.
  END SUBROUTINE read_multi_1d_file

END MODULE atmos_reader_mod
