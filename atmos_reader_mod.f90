!===============================================================================
!> Model-atmosphere file readers with no DISPATCH dependency.
!>
!> Every consumer of an atmosphere file shares these routines: read_atmos_mod copies
!> the parsed data into the runtime input_atmos_t (mesh, patches, AMR, pruning of the
!> top planes), hop_prune.x reduces it to the thermodynamic states it hands to
!> standalone HOP, and libhop's hop_atmos_* entry points hand it to Python
!> (tumult.atmos_reader) at the file's full resolution.  Everything a format-specific parser knows lands in atmos_data_t;
!> everything a consumer decides (namelist overrides of vmic or [Fe/H], which of
!> rho/pe the EOS uses, messages) stays with the consumer.
!>
!> The thermodynamic columns are stored as natural logarithms in the same real(4)
!> precision read_atmos_mod has always used, so the runtime copy is exact and the
!> binary formats that tabulate logarithms convert with one multiplication.
!> Arrays keep the file's own orientation and units after the format's conversion
!> to cgs (Stagger code units, km/s, momenta): read_atmos_mod%thick2thin decides
!> which end is the top.
!>
!> Failures are reported through ierr/msg rather than stop, so a standalone caller
!> can turn them into an exit status; read_atmos_mod stops on them itself.
!===============================================================================
MODULE atmos_reader_mod
  USE const_mod, only: bk
  USE handyfuncs_mod, only: nextentry, count_words, lowercase, endswith, replace_char
  implicit none
  private

  type, public :: atmos_data_t
    character(len=16) :: format = ''
    integer :: nx=0, ny=0, nz=0
    integer :: gz=0         !< ghost zones at each end in z (Stagger)
    real(4) :: timestamp=-1 !< CO5BOLD snapshot time
    !> ln T [K], ln P_e [dyn/cm2], ln rho [g/cm3]; a column a format does not carry
    !> is flagged absent and left unset.
    real(4), allocatable :: ltemp(:,:,:), lpe(:,:,:), lrho(:,:,:)
    logical :: has_pe=.false., has_rho=.false.
    real(4), allocatable :: vmic(:,:,:)                     !< [km/s], the file's own values
    real(4), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:) !< [cm/s] where the format has them, else 0
    real(4), allocatable :: xx(:), yy(:)                    !< horizontal mesh [cm] (3D formats), else 0
    real(4), allocatable :: zz(:)                           !< depth coordinate as the file gives it
    logical :: has_feh=.false.
    real(8) :: feh=0                 !< model [Fe/H] when the file states one
    real(4) :: radius=0              !< [cm], MARCS (1 for plane-parallel)
    real(4), allocatable :: abund(:) !< MARCS: log eps for Z=1..92
  contains
    procedure :: allocate => allocate_atmos_data
  end type atmos_data_t

  public :: read_atmos_file, read_marcs_file, read_text_1d_file, read_multi_1d_file
  public :: read_stagger_file, read_stagger2_file, read_multi_snapshot_file
  public :: read_co5bold_file, read_siu_bin_file
  public :: atmos_format_from_name

  ! Stagger code units
  real(4), parameter :: uul = 1.0e+8      ! length
  real(4), parameter :: uur = 1.0e-7      ! density
  real(4), parameter :: luul = log(uul)
  real(4), parameter :: luur = log(uur)
  real(4), parameter :: uut = 1.0e2       ! time
  real(4), parameter :: uuv = uul / uut   ! velocity

CONTAINS

  !=============================================================================
  SUBROUTINE allocate_atmos_data(self, nx, ny, nz)
    class(atmos_data_t), intent(inout) :: self
    integer, intent(in) :: nx, ny, nz
    !---------------------------------------------------------------------------
    self%nx = nx; self%ny = ny; self%nz = nz
    allocate(self%xx(nx), source=0.)
    allocate(self%yy(ny), source=0.)
    allocate(self%zz(nz), source=0.)
    allocate(self%vx(nx,ny,nz), source=0.)
    allocate(self%vy(nx,ny,nz), source=0.)
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
  !> Dispatch on the format names read_external_mod accepts ('save' is a run
  !> folder, not an atmosphere file, and returns ierr=2).  snap is the Stagger
  !> snapshot number, use_aux whether a Stagger snapshot takes T and n_e from its
  !> .aux file when there is one; both default as in &atmos_params.
  !=============================================================================
  SUBROUTINE read_atmos_file(atmos_format, atmos_file, atmos, ierr, msg, snap, use_aux)
    character(len=*), intent(in) :: atmos_format, atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    integer, intent(in), optional :: snap
    logical, intent(in), optional :: use_aux
    integer :: snap_
    logical :: aux_
    !---------------------------------------------------------------------------
    snap_ = 0
    if (present(snap)) snap_ = snap
    aux_ = .true.
    if (present(use_aux)) aux_ = use_aux
    select case (lowercase(trim(atmos_format)))
    case ('marcs')
      call read_marcs_file(atmos_file, atmos, ierr, msg)
    case ('text')
      call read_text_1d_file(atmos_file, atmos, ierr, msg)
    case ('multi1d')
      call read_multi_1d_file(atmos_file, atmos, ierr, msg)
    case ('stagger')
      call read_stagger_file(atmos_file, snap_, aux_, atmos, ierr, msg)
    case ('stagger2')
      call read_stagger2_file(atmos_file, atmos, ierr, msg)
    case ('multi', 'm3dis')
      call read_multi_snapshot_file(atmos_file, atmos, ierr, msg)
    case ('co5bold')
      call read_co5bold_file(atmos_file, atmos, ierr, msg)
    case ('siu_bin')
      call read_siu_bin_file(atmos_file, atmos, ierr, msg)
    case default
      ierr = 2
      msg = "atmos_format '" // trim(atmos_format) // "' is not supported by the " // &
            "standalone reader (marcs, text, multi1d, stagger, stagger2, multi, " // &
            "m3dis, co5bold, siu_bin)"
    end select
  END SUBROUTINE read_atmos_file

  !=============================================================================
  SUBROUTINE open_existing(atmos_file, lunit, ierr, msg, form, access, recl, convert)
    character(len=*), intent(in) :: atmos_file
    integer, intent(out) :: lunit, ierr
    character(len=*), intent(out) :: msg
    character(len=*), intent(in), optional :: form, access, convert
    integer, intent(in), optional :: recl
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
    if (present(recl)) then
      open(newunit=lunit, file=atmos_file, form=form, access=access, recl=recl, &
           convert=convert, status='old', iostat=iostat)
    elseif (present(convert)) then
      open(newunit=lunit, file=atmos_file, form=form, convert=convert, iostat=iostat)
    elseif (present(access)) then
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
    character(4096) :: fileID
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
  !> Text 1D: name, "nz [FeH]", then one "z T Pe rho vmic [vz]" row per depth.  A
  !> column of zeros marks Pe or rho as absent, and the logarithms are taken only
  !> where the value is positive, as read_atmos_mod always has.  The first row
  !> decides whether the rows carry the vertical velocity vz [km/s], positive
  !> outward whatever the orientation of the depth column; it is stored in cm/s.
  !=============================================================================
  SUBROUTINE read_text_1d_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    integer :: nwords, nz, lunit, i, iostat
    character(128) :: text
    character(256) :: row
    logical :: has_vz
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
      read(lunit, '(A)', iostat=iostat) row
      if (iostat == 0) then
        if (i == 1) then
          call replace_char(row, achar(9), ' ')
          call count_words(row, nwords)
          has_vz = nwords >= 6
        endif
        if (has_vz) then
          read(row, *, iostat=iostat) atmos%zz(i), temp(i), pe(i), rho(i), atmos%vmic(1,1,i), &
                                      atmos%vz(1,1,i)
        else
          read(row, *, iostat=iostat) atmos%zz(i), temp(i), pe(i), rho(i), atmos%vmic(1,1,i)
        endif
      endif
      if (iostat /= 0) then
        close(lunit)
        ierr = 3
        write(msg, '(A,I0,A)') 'Cannot parse depth ', i, ' of ' // trim(atmos_file)
        return
      endif
    end do
    close(lunit)

    atmos%vz = atmos%vz * 1e5
    atmos%has_pe = pe(1) /= 0
    atmos%has_rho = rho(1) /= 0
    where (temp > 0) atmos%ltemp(1,1,:) = log(temp)
    where (rho  > 0) atmos%lrho (1,1,:) = log(rho)
    where (pe   > 0) atmos%lpe  (1,1,:) = log(pe)
  END SUBROUTINE read_text_1d_file

  !=============================================================================
  !> MULTI 1D atmos.* : name, scale name, log g, nz, then "z T ne vz vmic" rows.
  !> The file carries electron density, stored here as ln P_e = ln(n_e k T); there
  !> is no mass density.  vz [km/s] is positive outward, as in MULTI, and is
  !> stored in cm/s.
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

    atmos%vz = atmos%vz * 1e5
    atmos%ltemp(1,1,:) = log(temp)
    atmos%lpe  (1,1,:) = log(ne * BK * temp)
    atmos%has_pe = .true.
    atmos%has_rho = .false.
  END SUBROUTINE read_multi_1d_file

  !=============================================================================
  !> CO5BOLD stream file: time, dimensions, mesh, ln T, ln rho, ln P_e and the
  !> momenta, which are divided by the density here.
  !=============================================================================
  SUBROUTINE read_co5bold_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    real(4), dimension(:,:,:), allocatable :: over_rho
    integer :: lunit, nx, ny, nz
    !---------------------------------------------------------------------------
    call open_existing(atmos_file, lunit, ierr, msg, form='unformatted', access='stream')
    if (ierr /= 0) return
    atmos%format = 'co5bold'
    read(lunit) atmos%timestamp
    read(lunit) nx, ny, nz
    call atmos%allocate(nx, ny, nz)
    read(lunit) atmos%xx, atmos%yy, atmos%zz
    read(lunit) atmos%ltemp
    read(lunit) atmos%lrho
    read(lunit) atmos%lpe
    read(lunit) atmos%vx ! this is px
    read(lunit) atmos%vy ! this is py
    read(lunit) atmos%vz ! this is pz
    close(lunit)

    allocate(over_rho, source=exp(-atmos%lrho))
    atmos%vx = atmos%vx * over_rho
    atmos%vy = atmos%vy * over_rho
    atmos%vz = atmos%vz * over_rho
    atmos%has_pe = .true.
    atmos%has_rho = .true.
  END SUBROUTINE read_co5bold_file

  !=============================================================================
  !> The original Stagger grid: <file>.msh (mesh) and <file>_<snap>.dat (density,
  !> momenta, energy, T) in code units, big-endian, with y the vertical axis and
  !> five ghost zones at each end.  With use_aux and an .aux file, T and n_e come
  !> from the .aux file (has_pe); otherwise T from the .dat file and no P_e.
  !=============================================================================
  SUBROUTINE read_stagger_file(stagger_file, snap, use_aux, atmos, ierr, msg)
    character(len=*), intent(in) :: stagger_file
    integer, intent(in) :: snap
    logical, intent(in) :: use_aux
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    character(len=4096) :: mesh_file, snap_file, data_file, aux_file
    real(4), allocatable, dimension(:) :: xdn, dxidxup, dxidxdn, dx, dxdn
    real(4), allocatable, dimension(:) :: ydn, dyidyup, dyidydn, dy, dydn
    real(4), allocatable, dimension(:) :: zdn, dzidzup, dzidzdn, dz, dzdn
    real(4), dimension(:,:,:), allocatable :: ppx, ppy, ppz, lne, rho
    integer :: iostat, lunit, recl, nx, ny, nz, i, j, k
    logical :: with_aux
    real(4) :: const, over_rho
    !---------------------------------------------------------------------------
    mesh_file = trim(stagger_file) // '.msh'
    write(snap_file, '(A,A1,I5.5)') trim(stagger_file), '_', snap
    write(data_file, '(A,A)') trim(snap_file), '.dat'
    write(aux_file,  '(A,A)') trim(snap_file), '.aux'

    ! mesh
    call open_existing(mesh_file, lunit, ierr, msg, form='unformatted', &
                       convert='big_endian')
    if (ierr /= 0) return
    atmos%format = 'stagger'
    read(unit=lunit, iostat=iostat) nx, nz, ny ! y is z in stagger
    call atmos%allocate(nx, ny, nz)
    atmos%gz = 5 ! stagger has 5 ghost zones on each side
    allocate(xdn, dxidxup, dxidxdn, dx, dxdn, mold=atmos%xx)
    allocate(ydn, dyidyup, dyidydn, dy, dydn, mold=atmos%yy)
    allocate(zdn, dzidzup, dzidzdn, dz, dzdn, mold=atmos%zz)
    read (unit=lunit, iostat=iostat) dx, dxdn, atmos%xx, xdn, dxidxup, dxidxdn
    read (unit=lunit, iostat=iostat) dz, dzdn, atmos%zz, zdn, dzidzup, dzidzdn ! y is z in stagger
    read (unit=lunit, iostat=iostat) dy, dydn, atmos%yy, ydn, dyidyup, dyidydn
    close(unit=lunit)
    atmos%xx = uul * atmos%xx
    atmos%yy = uul * atmos%yy
    atmos%zz = uul * atmos%zz

    ! data and aux
    allocate(ppx(nx,ny,nz), ppy(nx,ny,nz), ppz(nx,ny,nz), rho(nx,ny,nz), lne(nx,ny,nz))
    recl = kind(1.0)*nx*ny*nz
    inquire(file=aux_file, exist=with_aux)
    with_aux = use_aux .and. with_aux

    call open_existing(data_file, lunit, ierr, msg, form='unformatted', access='direct', &
                       recl=recl, convert='big_endian')
    if (ierr /= 0) return
    ! 1: density, 2-4: momentum per unit volume in x, y, z, 5: internal energy per
    ! unit volume, 6: temperature
    read (unit=lunit, iostat=iostat, rec=1) rho
    read (unit=lunit, iostat=iostat, rec=2) ppx
    read (unit=lunit, iostat=iostat, rec=3) ppz ! y is z in stagger
    read (unit=lunit, iostat=iostat, rec=4) ppy ! z is y in stagger
    if (.not. with_aux) then
      read (unit=lunit, iostat=iostat, rec=6) atmos%ltemp ! this is temp
      atmos%ltemp = log(atmos%ltemp)
    endif
    close(unit=lunit)

    if (with_aux) then
      call open_existing(aux_file, lunit, ierr, msg, form='unformatted', access='direct', &
                         recl=recl, convert='big_endian')
      if (ierr /= 0) return
      ! 1: ln total pressure (gas + rad), 2: ln Rosseland extinction coeff.,
      ! 3: ln T, 4: ln n_e, 5: ln Planck-averaged extinction coeff., 6: ln
      ! extinction coeff. at 500 nm
      read (unit=lunit, iostat=iostat, rec=3) atmos%ltemp
      read (unit=lunit, iostat=iostat, rec=4) lne
      close(unit=lunit)
    else
      atmos%lpe = 0
    endif

    ! "un"-staggering vectors relative to scalars
    call xup_set(ppx, atmos%vx)

    ! transposing data (y to z)
    call transpose_stagger(rho)
    call transpose_stagger(atmos%ltemp)
    if (with_aux) call transpose_stagger(lne)

    call transpose_stagger(atmos%vx)
    call transpose_stagger(ppy)
    call transpose_stagger(ppz)

    call yup_set(ppy, atmos%vy)
    call zup_set(ppz, atmos%vz)

    const = -3.0 * luul + log(BK)
    do k=1, nz
      do j=1, ny
        do i=1, nx
          over_rho = uuv / rho(i,j,k)
          atmos%lrho(i,j,k) = log(rho(i,j,k)) + luur
          atmos%vx  (i,j,k) = atmos%vx(i,j,k) * over_rho
          atmos%vy  (i,j,k) = atmos%vy(i,j,k) * over_rho
          atmos%vz  (i,j,k) = atmos%vz(i,j,k) * over_rho
          if (with_aux) then
            atmos%lpe(i,j,k) = lne(i,j,k) + const + atmos%ltemp(i,j,k)
          endif
        end do
      end do
    end do
    atmos%has_pe = with_aux
    atmos%has_rho = .true.
  END SUBROUTINE read_stagger_file

  !=============================================================================
  !> MULTI3D/M3DIS snapshot: <base>.mesh (formatted nx ny nz and the axes [cm])
  !> and <base>.atmos (stream: n_e, T, vx, vy, vz [cm/s], rho).  atmos_file may
  !> name either file or the common base.
  !=============================================================================
  SUBROUTINE read_multi_snapshot_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    character(len=:), allocatable :: base
    integer :: lunit, nx, ny, nz, nb
    real(4), dimension(:,:,:), allocatable :: ne, temp, rho
    !---------------------------------------------------------------------------
    nb = len_trim(atmos_file)
    if (endswith(atmos_file(1:nb), '.atmos')) then
      nb = nb - 6
    elseif (endswith(atmos_file(1:nb), '.mesh')) then
      nb = nb - 5
    endif
    allocate(base, source=atmos_file(1:nb))
    call open_existing(base // '.mesh', lunit, ierr, msg)
    if (ierr /= 0) return
    atmos%format = 'multi'
    read(lunit, *) nx, ny, nz
    call atmos%allocate(nx, ny, nz)
    allocate(temp, ne, rho, mold=atmos%lrho)
    read(lunit, *) atmos%xx
    read(lunit, *) atmos%yy
    read(lunit, *) atmos%zz
    close(lunit)

    call open_existing(base // '.atmos', lunit, ierr, msg, form='unformatted', access='stream')
    if (ierr /= 0) return
    read(lunit) ne
    read(lunit) temp
    read(lunit) atmos%vx
    read(lunit) atmos%vy
    read(lunit) atmos%vz
    read(lunit) rho
    close(lunit)

    atmos%lpe = log(ne * BK * temp)
    atmos%ltemp = log(temp)
    atmos%lrho  = log(rho)
    atmos%has_pe = .true.
    atmos%has_rho = .true.
  END SUBROUTINE read_multi_snapshot_file

  !=============================================================================
  !> Stagger grid 2 (Lagae): mesh.<id> (formatted, each axis as count then
  !> values) beside atm3d.<id> (stream: n_e, T, vx, vy, vz [km/s], rho).
  !=============================================================================
  SUBROUTINE read_stagger2_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    character(len=:), allocatable :: mesh_file
    integer :: lunit, nx, ny, nz, pos
    real(4), dimension(:), allocatable :: xx, yy, zz
    real(4), dimension(:,:,:), allocatable :: ne, temp, rho
    !---------------------------------------------------------------------------
    pos = index(atmos_file, 'atm3d.')
    if (pos < 1) then
      ierr = 3
      msg = "a stagger2 file name contains 'atm3d.' - " // trim(atmos_file)
      return
    endif
    allocate(mesh_file, source=atmos_file(1:pos-1) // 'mesh.' // trim(atmos_file(pos+6:)))
    call open_existing(mesh_file, lunit, ierr, msg)
    if (ierr /= 0) return
    atmos%format = 'stagger2'
    read(lunit, *) nx
    allocate(xx(nx))
    read(lunit, *) xx
    read(lunit, *) ny
    allocate(yy(ny))
    read(lunit, *) yy
    read(lunit, *) nz
    allocate(zz(nz))
    read(lunit, *) zz
    close(lunit)

    call atmos%allocate(nx, ny, nz)
    atmos%xx = xx
    atmos%yy = yy
    atmos%zz = zz

    allocate(temp, ne, rho, mold=atmos%lrho)
    call open_existing(atmos_file, lunit, ierr, msg, form='unformatted', access='stream')
    if (ierr /= 0) return
    read(lunit) ne
    read(lunit) temp
    read(lunit) atmos%vx
    read(lunit) atmos%vy
    read(lunit) atmos%vz
    read(lunit) rho
    close(lunit)

    atmos%vx = atmos%vx * 1e5
    atmos%vy = atmos%vy * 1e5
    atmos%vz = atmos%vz * 1e5
    atmos%lpe = log(ne * BK * temp)
    atmos%ltemp = log(temp)
    atmos%lrho  = log(rho)
    atmos%has_pe = .true.
    atmos%has_rho = .true.
  END SUBROUTINE read_stagger2_file

  !=============================================================================
  !> SIU binary 1D model: log10 tau_500 as depth, log10 of T, P_e and rho.
  !=============================================================================
  SUBROUTINE read_siu_bin_file(atmos_file, atmos, ierr, msg)
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: msg
    integer :: lunit, nz, nel, nmol, nbreak, recl
    real(4), parameter :: ln10 = log(10.0)
    integer, dimension(:), allocatable :: at_idx, at_nion, break
    real(4), dimension(:), allocatable :: kappa500, pgas
    !---------------------------------------------------------------------------
    call open_existing(atmos_file, lunit, ierr, msg, form='unformatted', access='stream')
    if (ierr /= 0) return
    atmos%format = 'siu_bin'
    read(lunit) recl ! skip the record length in stream access
    read(lunit) nz, nel, nmol, nbreak

    call atmos%allocate(1,1,nz)
    allocate(kappa500(nz), pgas(nz))
    allocate(at_idx(nel), at_nion(nel))
    allocate(break(nbreak))

    read(lunit) at_idx           ! Atomic number of the elements
    read(lunit) at_nion          ! Number of ionisation stages considered for each element
    read(lunit) break            ! break points in tau scale for MARCS
    read(lunit) atmos%zz         ! log10(tau500)
    read(lunit) atmos%ltemp      ! log10(temperature)
    read(lunit) kappa500         ! log10(kappa500)
    read(lunit) atmos%lpe        ! log10(electron pressure)
    read(lunit) pgas             ! p gas MARCS
    read(lunit) pgas             ! p gas MARCS (not sure why this appears twice)
    read(lunit) atmos%lrho       ! log10(density)
    ! the data goes on, but we have everything we need
    close(lunit)

    ! convert from log10() to ln
    atmos%ltemp = atmos%ltemp * ln10
    atmos%lpe   = atmos%lpe   * ln10
    atmos%lrho  = atmos%lrho  * ln10
    atmos%has_pe = .true.
    atmos%has_rho = .true.
  END SUBROUTINE read_siu_bin_file

  !=============================================================================
  !> f is centered on (i-.5,j,k), fifth order stagger
  !=============================================================================
  SUBROUTINE xup_set (f, fp)
    real(4), dimension(:,:,:), intent(in) :: f
    real(4), dimension(:,:,:), intent(out) :: fp
    integer :: j, k
    !---------------------------------------------------------------------------
    do k=1, size(f, 3)
      do j=1, size(f, 2)
        call destagger_vector(f(:,j,k), fp(:,j,k))
      end do
    end do
  END SUBROUTINE xup_set

  !=============================================================================
  SUBROUTINE yup_set (f, fp)
    real(4), dimension(:,:,:), intent(in) :: f
    real(4), dimension(:,:,:), intent(out) :: fp
    integer :: i, k
    !---------------------------------------------------------------------------
    do k=1, size(f, 3)
      do i=1, size(f, 1)
        call destagger_vector(f(i,:,k), fp(i,:,k))
      end do
    end do
  END SUBROUTINE yup_set

  !=============================================================================
  SUBROUTINE zup_set (f, fp)
    real(4), dimension(:,:,:), intent(in) :: f
    real(4), dimension(:,:,:), intent(out) :: fp
    integer :: i, j
    !---------------------------------------------------------------------------
    do j=1, size(f, 2)
      do i=1, size(f, 1)
        call destagger_vector(f(i,j,:), fp(i,j,:))
      end do
    end do
  END SUBROUTINE zup_set

  !=============================================================================
  !> Periodic fifth-order interpolation from the staggered (i-1/2) to the centred
  !> position.
  !=============================================================================
  SUBROUTINE destagger_vector(f, fout)
    real, dimension(:), intent(in)  :: f
    real, dimension(:), intent(out) :: fout
    integer :: i, n, j, idx
    real(4), parameter :: c = 3./256.
    real(4), parameter :: b = -25./256.
    real(4), parameter :: a = 0.5-b-c
    real(4) :: wts(-2:3)=[c,b,a,a,b,c]
    !---------------------------------------------------------------------------
    n = size(f,1)
    fout = 0
    do i=1, n
      do j=-2,3
        idx = i + j
        if (idx<1) then
          idx = n + idx
        elseif (idx>n) then
          idx = idx-n
        endif
        fout(i) = fout(i) + wts(j) * f(idx)
      end do
    end do
  END SUBROUTINE destagger_vector

  !=============================================================================
  !> Swap the y and z axes of every x-plane (Stagger's vertical axis is y).
  !=============================================================================
  SUBROUTINE transpose_stagger(f)
    real, dimension(:,:,:), intent(inout)  :: f
    integer :: i
    !---------------------------------------------------------------------------
    do i=1, size(f, 1)
      f(i, :, :) = transpose(f(i, :, :))
    end do
  END SUBROUTINE transpose_stagger

END MODULE atmos_reader_mod
