!===============================================================================
!> Conversion of MULTI-format model atoms (atom.*) to the tabular ("pandas")
!> format that atom_mod%read_pandas reads.
!>
!> dispatch.x reads every model atom through read_pandas: a MULTI atom is parsed
!> here, written as <run_dir>/atoms/<id>.txt and read back, so the copy in the run
!> directory is the atom the run used.  atom_convert.x (HOP/standalone) is the same
!> conversion as a program for the Python layer.
!>
!> The parsing follows the MULTI conventions dispatch has always applied:
!>  - level energies are given in cm^-1; a line's vacuum wavelength follows from
!>    its level energies and is written as an air wavelength, which read_pandas
!>    converts back;
!>  - a line with io >= 2 components becomes one line per component with positive
!>    weight: the component's wavelength is shifted by dq [cm^-1], its gf is the
!>    line's gf times the weight and its radiative damping is the component's.
!>    The first component keeps the line's position in the table, the others
!>    follow the MULTI lines in order, and parent() names the MULTI line of every
!>    written line;
!>  - hydrogenic continua (qmax >= 0) are tabulated on the nodes MULTI defines:
!>    nq points linear in wavelength from the edge to qmax, alpha0 (lam/lam0)^3;
!>  - a continuum whose upper level is not the ground level of the next ion is
!>    written with the optional j column of the bound-free table;
!>  - collision records on several TEMP grids are interpolated onto the union of
!>    the grids with gencol's scheme (cubic Hermite with monotone slopes, in log
!>    for OHM, OHMEGA, CE, CI and CH0, constant outside the record's grid), which
!>    is exact on each record's own nodes.
!>
!> The continuum header's alpha0 has no column in the table format; read_pandas
!> takes a continuum's threshold cross section from the first tabulated node.
!===============================================================================
MODULE atom_convert_mod
  USE const_mod, only: hc_e, h_e, HZ_TO_AA, AA_TO_CM
  USE handyfuncs_mod, only: NextEntry, print_failed_line, vac2air, uppercase, lowercase
  USE interpol_mod, only: hermite_derivative, hermite_interpolation
  USE console_mod
  implicit none
  private

  type, public :: classic_atom_t
    character(len=20) :: id = ''
    real(8) :: abnd = 0, awgt = 0
    integer :: nlevel = 0, nline = 0, nbf = 0, ncol = 0
    real(8), allocatable :: ev(:), g(:)               ! [eV]
    integer, allocatable :: ion(:)
    character(len=64), allocatable :: label(:)
    ! bound-bound records as written, components already split
    integer, allocatable :: li(:), lj(:), nq(:), parent(:)
    real(8), allocatable :: lam_vac(:), loggf(:), grad(:), gvdw(:), gstk(:), qmax(:), q0(:)
    ! bound-free nodes, each continuum from its threshold to shorter wavelengths
    integer, allocatable :: bf_i(:), bf_j(:)
    real(8), allocatable :: bf_lam(:), bf_alpha(:)  ! [AA], [cm^2]
    logical :: excited_bf = .false.
    ! collision records on one temperature grid
    real(8), allocatable :: temp(:), col_val(:,:)
    character(len=8), allocatable :: col_key(:)
    integer, allocatable :: col_l(:), col_u(:)
  end type classic_atom_t

  ! a MULTI line with its components (none for a single line)
  type :: line_rec_t
    integer :: i, j, nq, io
    real(8) :: f, qmax, q0, grad, gvdw, gstk
    real(8), allocatable :: dq(:), wt(:), da(:)
  end type line_rec_t

  character(len=8), dimension(7), parameter :: processes = &
    [character(len=8) :: 'OHMEGA', 'CI', 'CH', 'CE', 'OHM', 'CH0', 'CHI']
  character(len=8), dimension(5), parameter :: log_processes = &
    [character(len=8) :: 'OHMEGA', 'OHM', 'CE', 'CI', 'CH0']

  public :: read_classic_atom, write_pandas_atom, convert_classic_atom
CONTAINS

!===============================================================================
!> Convert the MULTI atom in_file to the table format in out_file.
!===============================================================================
SUBROUTINE convert_classic_atom (in_file, out_file)
  character(len=*), intent(in) :: in_file, out_file
  type(classic_atom_t) :: atom
  !-----------------------------------------------------------------------------
  call read_classic_atom(in_file, atom)
  call write_pandas_atom(atom, out_file, in_file)
END SUBROUTINE convert_classic_atom

!===============================================================================
!> Parse the MULTI atom in atom_file.
!===============================================================================
SUBROUTINE read_classic_atom (atom_file, atom)
  character(len=*), intent(in) :: atom_file
  type(classic_atom_t), intent(out) :: atom
  type(line_rec_t), allocatable :: lines(:)
  integer :: latom, ios, nline, ncont, nrfix, i
  logical :: exists, abort
  character(len=256) :: msg
  !-----------------------------------------------------------------------------
  inquire(file=atom_file, exist=exists)
  if (.not. exists) then
    call console_error('model atom '//trim(atom_file)//' not found')
    error stop 'atom_convert_mod%read_classic_atom'
  endif
  open(newunit=latom, file=atom_file, form='formatted', status='old', action='read')
  abort = .false.

  call NextEntry(latom)
  read(latom, fmt='(a)') atom%id
  atom%id = uppercase(adjustl(atom%id))

  call NextEntry(latom)
  read(latom, *, iostat=ios) atom%abnd, atom%awgt
  call print_failed_line(latom, ios)

  call NextEntry(latom)
  read(latom, *, iostat=ios) atom%nlevel, nline, ncont, nrfix
  call print_failed_line(latom, ios)

  allocate(atom%ev(atom%nlevel), atom%g(atom%nlevel), atom%ion(atom%nlevel), &
           atom%label(atom%nlevel))
  do i=1, atom%nlevel
    call NextEntry(latom)
    read(latom, *, iostat=ios) atom%ev(i), atom%g(i), atom%label(i), atom%ion(i)
    call print_failed_line(latom, ios)
    atom%ev(i) = atom%ev(i) * hc_e
  end do
  do i=2, atom%nlevel
    if (atom%ev(i) < atom%ev(i-1)) then
      call console_error('level energies not increasing at '//trim(atom%label(i))// &
                         ' in '//trim(atom_file))
      abort = .true.
    endif
  end do

  allocate(lines(nline))
  call read_lines(latom, atom, lines, abort)
  call split_components(atom, lines, abort)
  call read_continua(latom, atom, ncont, abort)
  call read_collisions(latom, atom, abort)
  close(latom)

  if (abort) then
    write(msg, '(a)') 'cannot convert model atom '//trim(atom_file)
    call console_error(msg)
    error stop 'atom_convert_mod%read_classic_atom'
  endif
END SUBROUTINE read_classic_atom

!===============================================================================
!> The bound-bound records and their components.
!===============================================================================
SUBROUTINE read_lines (latom, atom, lines, abort)
  integer, intent(in) :: latom
  type(classic_atom_t), intent(in) :: atom
  type(line_rec_t), intent(inout) :: lines(:)
  logical, intent(inout) :: abort
  integer :: i, k, ii, jj, ios
  real(8) :: fterm, gw, gq
  character(len=512) :: text
  character(len=256) :: msg
  !-----------------------------------------------------------------------------
  do i=1, size(lines)
    associate (lin => lines(i))
      call NextEntry(latom)
      ! the closing PROFILE word is not read: every line has a Voigt profile, and a
      ! hydrogen atom takes the Stark profiles
      read(latom, *, iostat=ios) jj, ii, lin%f, lin%nq, lin%qmax, lin%q0, lin%io, &
        lin%grad, lin%gvdw, lin%gstk
      call print_failed_line(latom, ios)
      lin%i = min(ii, jj)
      lin%j = max(ii, jj)
      if (lin%j > atom%nlevel .or. lin%i < 1) then
        write(msg, '(a,3(1x,i0))') 'transition to a non-existing level:', i, lin%i, lin%j
        call console_error(msg)
        abort = .true.
      endif
      if (lin%io >= 2) then
        allocate(lin%dq(lin%io), lin%wt(lin%io), lin%da(lin%io))
        do k=1, lin%io
          read(latom, '(a)') text
          read(text, *, iostat=ios) lin%dq(k), lin%wt(k), fterm, lin%da(k), gw, gq
          call print_failed_line(latom, ios, text)
        end do
      endif
    end associate
  end do
END SUBROUTINE read_lines

!===============================================================================
!> One written line per MULTI line or per component of positive weight.
!===============================================================================
SUBROUTINE split_components (atom, lines, abort)
  type(classic_atom_t), intent(inout) :: atom
  type(line_rec_t), intent(in) :: lines(:)
  logical, intent(inout) :: abort
  integer :: i, k, n, first, nextra
  real(8) :: lam0
  character(len=256) :: msg
  !-----------------------------------------------------------------------------
  n = size(lines)
  do i=1, size(lines)
    if (lines(i)%io >= 2) n = n + max(count(lines(i)%wt > 0) - 1, 0)
  end do
  atom%nline = n
  allocate(atom%li(n), atom%lj(n), atom%nq(n), atom%parent(n), atom%lam_vac(n), &
           atom%loggf(n), atom%grad(n), atom%gvdw(n), atom%gstk(n), atom%qmax(n), &
           atom%q0(n))

  nextra = size(lines)
  do i=1, size(lines)
    associate (lin => lines(i))
      lam0 = HZ_TO_AA * h_e / (atom%ev(lin%j) - atom%ev(lin%i))
      call set_line(i, i, lam0, lin%f, lin%grad)
      if (lin%io < 2) cycle
      first = findloc(lin%wt > 0, .true., 1)
      if (first == 0) then
        write(msg, '(a,i0,a)') 'line ', i, ' has no component of positive weight'
        call console_error(msg)
        abort = .true.
        cycle
      endif
      do k=first, lin%io
        if (lin%wt(k) <= 0) cycle
        ! dq is the component's wavenumber offset [cm^-1]; build_ll places it at
        ! lam0*exp(-dq*lam0) as well
        if (k == first) then
          call set_line(i, i, lam0*exp(-lin%dq(k)*lam0*AA_TO_CM), lin%f*lin%wt(k), lin%da(k))
        else
          nextra = nextra + 1
          call set_line(nextra, i, lam0*exp(-lin%dq(k)*lam0*AA_TO_CM), lin%f*lin%wt(k), &
                        lin%da(k))
        endif
      end do
    end associate
  end do
CONTAINS
  SUBROUTINE set_line (l, m, lam, f, grad)
    integer, intent(in) :: l, m
    real(8), intent(in) :: lam, f, grad
    atom%parent(l)  = m
    atom%li(l)      = lines(m)%i
    atom%lj(l)      = lines(m)%j
    atom%lam_vac(l) = lam
    atom%loggf(l)   = log10(atom%g(lines(m)%i) * f)
    atom%grad(l)    = grad
    atom%gvdw(l)    = lines(m)%gvdw
    atom%gstk(l)    = lines(m)%gstk
    atom%nq(l)      = lines(m)%nq
    atom%qmax(l)    = lines(m)%qmax
    atom%q0(l)      = lines(m)%q0
  END SUBROUTINE set_line
END SUBROUTINE split_components

!===============================================================================
!> The bound-free records: tabulated (qmax < 0) or hydrogenic.  A second record
!> for the same lower level is ignored, as dispatch always did.
!===============================================================================
SUBROUTINE read_continua (latom, atom, ncont, abort)
  integer, intent(in) :: latom, ncont
  type(classic_atom_t), intent(inout) :: atom
  logical, intent(inout) :: abort
  integer :: c, k, ii, jj, lo, up, ground, nnu, n, ios
  integer, allocatable :: ilo(:), iup(:), nodes(:)
  real(8), allocatable :: lam(:,:), alpha(:,:)
  real(8) :: alpha0, lam_min, lam0, dlam
  logical, allocatable :: done(:)
  character(len=256) :: msg
  !-----------------------------------------------------------------------------
  allocate(ilo(ncont), iup(ncont), nodes(ncont), source=0)
  allocate(done(atom%nlevel), source=.false.)
  allocate(lam(1,ncont), alpha(1,ncont))
  do c=1, ncont
    call NextEntry(latom)
    read(latom, *, iostat=ios) jj, ii, alpha0, nnu, lam_min
    call print_failed_line(latom, ios)
    lo = min(ii, jj)
    up = max(ii, jj)
    if (atom%ion(lo) == atom%ion(up)) then
      write(msg, '(a,i0,a,2(1x,i0))') 'bound-free transition ', c, &
        ' connects two levels of the same ion:', lo, up
      call console_error(msg)
      abort = .true.
    endif
    if (done(lo)) then
      write(msg, '(a,i0,a)') 'multiple bound-free transitions from lower state ', lo, &
        '; ignoring the later one'
      call console_warn(ch_m3d, msg)
      if (lam_min < 0) then
        do k=1, nnu
          call NextEntry(latom)
          read(latom, '(a)')
        end do
      endif
      cycle
    endif
    done(lo) = .true.
    if (nnu > size(lam, 1)) call grow(nnu)
    ilo(c) = lo
    iup(c) = up
    nodes(c) = nnu
    if (lam_min < 0) then
      do k=1, nnu
        call NextEntry(latom)
        read(latom, *, iostat=ios) lam(k,c), alpha(k,c)
        call print_failed_line(latom, ios)
      end do
    else
      lam0 = HZ_TO_AA * h_e / (atom%ev(up) - atom%ev(lo))
      dlam = (lam0 - lam_min) / (nnu-1)
      do k=1, nnu
        lam(k,c) = lam0 - (k-1)*dlam
        alpha(k,c) = (lam(k,c) / lam0)**3 * alpha0
      end do
    endif
  end do

  n = sum(nodes)
  atom%nbf = n
  allocate(atom%bf_i(n), atom%bf_j(n), atom%bf_lam(n), atom%bf_alpha(n))
  n = 0
  do c=1, ncont
    if (nodes(c) == 0) cycle
    ground = findloc(atom%ion(ilo(c):), atom%ion(ilo(c))+1, 1) + ilo(c) - 1
    if (iup(c) /= ground) atom%excited_bf = .true.
    atom%bf_i(n+1:n+nodes(c)) = ilo(c)
    atom%bf_j(n+1:n+nodes(c)) = iup(c)
    atom%bf_lam(n+1:n+nodes(c)) = lam(1:nodes(c),c)
    atom%bf_alpha(n+1:n+nodes(c)) = alpha(1:nodes(c),c)
    n = n + nodes(c)
  end do
CONTAINS
  SUBROUTINE grow (m)
    integer, intent(in) :: m
    real(8), allocatable :: tmp(:,:)
    allocate(tmp(m,ncont), source=0d0)
    tmp(1:size(lam,1),:) = lam
    call move_alloc(tmp, lam)
    allocate(tmp(m,ncont), source=0d0)
    tmp(1:size(alpha,1),:) = alpha
    call move_alloc(tmp, alpha)
  END SUBROUTINE grow
END SUBROUTINE read_continua

!===============================================================================
!> The GENCOL section, onto one temperature grid.
!===============================================================================
SUBROUTINE read_collisions (latom, atom, abort)
  integer, intent(in) :: latom
  type(classic_atom_t), intent(inout) :: atom
  logical, intent(inout) :: abort
  character(len=20) :: crout, txt
  character(len=8) :: key
  character(len=256) :: msg
  integer :: ios, ncoll, ngrids, max_ntemp, ntemp, igrid, i, l, u, n, k
  integer, allocatable :: grid_of(:), ntemps(:)
  real(8), allocatable :: tgrids(:,:), values(:,:), union(:), slope(:), y(:)
  logical :: same
  !-----------------------------------------------------------------------------
  call NextEntry(latom)
  read(latom, fmt='(a)') crout
  crout = lowercase(adjustl(crout))
  select case (crout)
  case ('end')
    return
  case ('gencol')
  case default
    call console_error('collision routine '//trim(crout)//' not supported')
    abort = .true.
    return
  end select

  ! first pass: count the records and the temperature grids
  ncoll = 0
  ngrids = 0
  max_ntemp = 1
  do
    call NextEntry(latom)
    read(latom, '(a)', iostat=ios) key
    if (ios < 0) then
      call console_error('read_collisions: keyword "END" not found')
      abort = .true.
      return
    endif
    key = uppercase(adjustl(key))
    if (scan(key(1:1), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ') == 0) cycle
    if (key == 'END') exit
    if (key == 'TEMP') then
      call NextEntry(latom)
      read(latom, *) ntemp
      max_ntemp = max(max_ntemp, ntemp)
      ngrids = ngrids + 1
    elseif (any(processes == key)) then
      ncoll = ncoll + 1
    else
      call console_error('read_collisions: unrecognized keyword '//trim(key))
      abort = .true.
      return
    endif
  end do

  allocate(tgrids(max_ntemp, ngrids), ntemps(ngrids), values(max_ntemp, ncoll))
  allocate(atom%col_key(ncoll), atom%col_l(ncoll), atom%col_u(ncoll), grid_of(ncoll))

  ! second pass from the routine's keyword
  rewind(latom)
  do
    call NextEntry(latom)
    read(latom, *) txt
    if (lowercase(adjustl(txt)) == crout) exit
  end do
  igrid = 0
  n = 0
  do
    call NextEntry(latom)
    read(latom, fmt='(a8)') key
    key = uppercase(adjustl(key))
    if (key == 'END') exit
    call NextEntry(latom)
    if (key == 'TEMP') then
      igrid = igrid + 1
      read(latom, *, iostat=ios) ntemps(igrid), tgrids(1:ntemps(igrid), igrid)
      call print_failed_line(latom, ios)
      if (ntemps(igrid) < 2) then
        call console_error('read_collisions: ntemp < 2')
        abort = .true.
        return
      endif
    elseif (any(processes == key)) then
      read(latom, *, iostat=ios) l, u, values(1:ntemps(igrid), n+1)
      call print_failed_line(latom, ios)
      if (l == u) then
        write(msg, '(3a,2(1x,i0),a)') 'll==ul for collision process ', trim(key), ':', &
          l, u, '; ignored'
        call console_warn(ch_m3d, msg)
        cycle
      endif
      n = n + 1
      if (any(values(1:ntemps(igrid), n) <= 0)) then
        write(msg, '(2a,2(1x,i0))') 'non-positive values for collision process ', &
          trim(key), l, u
        call console_error(msg)
        abort = .true.
      endif
      if (min(l,u) < 1 .or. max(l,u) > atom%nlevel) then
        write(msg, '(3a,3(a,i0))') 'level index outside range (1,nk) for ', trim(key), &
          ':', ' ll=', l, ' ul=', u, ' nk=', atom%nlevel
        call console_error(msg)
        abort = .true.
      endif
      atom%col_key(n) = key
      atom%col_l(n) = min(l, u)
      atom%col_u(n) = max(l, u)
      grid_of(n) = igrid
    endif
  end do
  atom%ncol = n

  ! one grid if every record's grid equals the first one used
  same = .true.
  do i=2, n
    k = grid_of(i)
    if (ntemps(k) /= ntemps(grid_of(1))) then
      same = .false.
    elseif (any(tgrids(1:ntemps(k), k) /= tgrids(1:ntemps(k), grid_of(1)))) then
      same = .false.
    endif
  end do
  if (n == 0) then
    allocate(atom%temp(0), atom%col_val(0,0))
    return
  endif
  if (same) then
    k = ntemps(grid_of(1))
    allocate(atom%temp(k), source=tgrids(1:k, grid_of(1)))
    allocate(atom%col_val(k, n))
    atom%col_val = values(1:k, 1:n)
    return
  endif

  call union_grid(tgrids, ntemps, grid_of(1:n), union)
  allocate(atom%temp(size(union)), source=union)
  allocate(atom%col_val(size(union), n))
  allocate(slope(max_ntemp), y(max_ntemp))
  do i=1, n
    k = ntemps(grid_of(i))
    associate (t => tgrids(1:k, grid_of(i)))
      if (any(log_processes == atom%col_key(i))) then
        y(1:k) = log(values(1:k, i))
      else
        y(1:k) = values(1:k, i)
      endif
      call hermite_derivative(slope(1:k), t, y(1:k))
      do l=1, size(union)
        call hermite_interpolation(atom%col_val(l,i), union(l), t, y(1:k), slope(1:k))
      end do
      if (any(log_processes == atom%col_key(i))) atom%col_val(:,i) = exp(atom%col_val(:,i))
    end associate
  end do
END SUBROUTINE read_collisions

!===============================================================================
!> Sorted union of the temperature grids that records use.
!===============================================================================
SUBROUTINE union_grid (tgrids, ntemps, used, union)
  real(8), intent(in) :: tgrids(:,:)
  integer, intent(in) :: ntemps(:), used(:)
  real(8), allocatable, intent(out) :: union(:)
  real(8), allocatable :: t(:)
  integer :: k, n, i, m
  logical :: keep(size(ntemps))
  !-----------------------------------------------------------------------------
  keep = .false.
  keep(used) = .true.
  allocate(t(sum(ntemps, mask=keep)))
  n = 0
  do k=1, size(ntemps)
    if (.not. keep(k)) cycle
    t(n+1:n+ntemps(k)) = tgrids(1:ntemps(k), k)
    n = n + ntemps(k)
  end do
  ! insertion sort; the grids have a few tens of nodes
  do i=2, n
    do m=i, 2, -1
      if (t(m-1) <= t(m)) exit
      t(m-1:m) = [t(m), t(m-1)]
    end do
  end do
  m = 1
  do i=2, n
    if (t(i) /= t(m)) then
      m = m + 1
      t(m) = t(i)
    endif
  end do
  allocate(union(m), source=t(1:m))
END SUBROUTINE union_grid

!===============================================================================
!> Write atom in the table format.  Indices are 0-based, as pandas writes them;
!> source names the converted file in the header comment.  Each column has the
!> precision tumult.atom.Atom.write gives it, so a converted atom and one written
!> from Python agree; the run that converts a MULTI atom uses the rounded values.
!===============================================================================
SUBROUTINE write_pandas_atom (atom, out_file, source)
  type(classic_atom_t), intent(in) :: atom
  character(len=*), intent(in) :: out_file, source
  integer :: u, i, k, wl, wk
  ! field widths without the separating blank
  integer, parameter :: w_awgt=8, w_n=6, w_ev=11, w_g=4, w_ion=3, w_lam=14, w_loggf=8, &
    w_grad=10, w_gvdw=16, w_gstk=11, w_nq=4, w_q=6, w_alpha=10, w_col=11
  character(len=*), parameter :: awgt = '1x,f8.3', ev = '1x,f11.6', g = '1x,i4', ion = '1x,i3', &
    lam = '1x,f14.3', loggf = '1x,f8.3', grad = '1x,es10.3', gvdw = '1x,g16.8', &
    gstk = '1x,g11.3', nq = '1x,i4', q = '1x,f6.1', alpha = '1x,es10.3', col = '1x,es11.4'
  character(len=8) :: idx, lvl, lvl0, key
  !-----------------------------------------------------------------------------
  open(newunit=u, file=out_file, form='formatted', status='replace', action='write')
  ! the MULTI comments are not copied: they describe the MULTI layout and units
  write(u, '(a)') '* converted from the MULTI-format atom '//trim(source)
  write(u, '(a)') trim(atom%id)
  write(u, '(a)') title('abnd', 5)//title('awgt', w_awgt)
  write(u, '(f6.2,'//awgt//')') atom%abnd, atom%awgt
  write(u, '(a)') title('nlevel', 5)//title('nline', w_n)//title('nbf', w_n)//title('ncol', w_n)
  write(u, '(i6,3(1x,i6))') atom%nlevel, atom%nline, atom%nbf, atom%ncol

  ! the level indices i and j share one width; idx is the table's own index
  wl = ndigits(atom%nlevel-1)
  write(lvl, '(a,i0)') '1x,i', wl
  write(lvl0, '(a,i0)') 'i', wl

  call set_index(atom%nlevel)
  write(u, '(a)') repeat(' ', ndigits(atom%nlevel-1))//title('ev', w_ev)//title('g', w_g)// &
    title('ion', w_ion)//' label'
  do i=1, atom%nlevel
    write(u, '('//trim(idx)//','//ev//','//g//','//ion//',3a)') i-1, atom%ev(i), &
      nint(atom%g(i)), atom%ion(i), " '", trim(atom%label(i)), "'"
  end do

  if (atom%nline > 0) then
    call set_index(atom%nline)
    write(u, '(a)') repeat(' ', ndigits(atom%nline-1))//title('i', wl)//title('j', wl)// &
      title('lam_air', w_lam)//title('evi', w_ev)//title('loggf', w_loggf)// &
      title('grad', w_grad)//title('gvdw', w_gvdw, 4)//title('gstk', w_gstk, 4)// &
      title('nq', w_nq)//title('qmax', w_q)//title('q0', w_q)//' iso'
    do i=1, atom%nline
      write(u, '('//trim(idx)//',2('//trim(lvl)//'),'//lam//','//ev//','//loggf//','// &
        grad//','//gvdw//','//gstk//','//nq//',2('//q//'),a)') i-1, &
        atom%li(i)-1, atom%lj(i)-1, vac2air(atom%lam_vac(i)), atom%ev(atom%li(i)), &
        atom%loggf(i), atom%grad(i), atom%gvdw(i), atom%gstk(i), atom%nq(i), &
        atom%qmax(i), atom%q0(i), '   0'
    end do
  endif

  ! the bound-free wavelengths are vacuum wavelengths, under the column name
  ! read_pandas has always used; this table has no index column
  if (atom%nbf > 0) then
    if (atom%excited_bf) then
      write(u, '(a)') title('i', wl-1)//title('lam_air', w_lam)//title('alpha', w_alpha)// &
        title('j', wl)
      do i=1, atom%nbf
        write(u, '('//trim(lvl0)//','//lam//','//alpha//','//trim(lvl)//')') atom%bf_i(i)-1, &
          atom%bf_lam(i), atom%bf_alpha(i), atom%bf_j(i)-1
      end do
    else
      write(u, '(a)') title('i', wl-1)//title('lam_air', w_lam)//title('alpha', w_alpha)
      do i=1, atom%nbf
        write(u, '('//trim(lvl0)//','//lam//','//alpha//')') atom%bf_i(i)-1, atom%bf_lam(i), &
          atom%bf_alpha(i)
      end do
    endif
  endif

  ! the temperatures are the titles of the rate columns
  if (atom%ncol > 0) then
    wk = max(len('key'), maxval(len_trim(atom%col_key)))
    write(key, '(a,i0)') '1x,a', wk
    call set_index(atom%ncol)
    write(u, '(a,*('//col//'))') repeat(' ', ndigits(atom%ncol-1))//title('key', wk)// &
      title('i', wl)//title('j', wl), atom%temp
    do i=1, atom%ncol
      write(u, '('//trim(idx)//','//trim(key)//',2('//trim(lvl)//'),*('//col//'))') i-1, &
        adjustr(atom%col_key(i)(1:wk)), atom%col_l(i)-1, atom%col_u(i)-1, &
        (atom%col_val(k,i), k=1, size(atom%temp))
    end do
  endif
  close(u)
CONTAINS
  ! edit descriptor of the index column of a table with n rows
  SUBROUTINE set_index (n)
    integer, intent(in) :: n
    write(idx, '(a,i0)') 'i', ndigits(n-1)
  END SUBROUTINE set_index
  ! name right-aligned over a field of width w after its separating blank; a G edit
  ! leaves trail blanks after a fixed-point value
  FUNCTION title (name, w, trail) result (s)
    character(len=*), intent(in) :: name
    integer, intent(in) :: w
    integer, intent(in), optional :: trail
    character(len=w+1) :: s
    integer :: t
    t = 0
    if (present(trail)) t = trail
    s = repeat(' ', w+1-t-len(name))//name
  END FUNCTION title
  INTEGER FUNCTION ndigits (n)
    integer, intent(in) :: n
    integer :: m
    ndigits = 1
    m = n
    do while (m >= 10)
      m = m/10
      ndigits = ndigits+1
    end do
  END FUNCTION ndigits
END SUBROUTINE write_pandas_atom

END MODULE atom_convert_mod
