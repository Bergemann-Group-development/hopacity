!===============================================================================
!> linelist_reader_mod -- random access to the line-list text files.
!>
!> A master line list is a sequence of per-species blocks, each a two-line header
!> ("'Fe I'" / "'26' 1234") followed by that many records, ascending in air
!> wavelength.  A run keeps the records inside a handful of windows and ignores
!> the rest, so reading the whole file to find them costs time proportional to
!> the list rather than to the selection: the masters run to a gigabyte, and the
!> newline index alone over 23 million records costs seconds and hundreds of
!> megabytes before the first record is parsed.
!>
!> The records of a block are almost always written with one fixed width, which
!> makes the k-th record addressable as dstart + k*width and the block bisectable
!> without reading anything between the probes.  The width is probed from the
!> first record and then verified (the block must end on a record boundary, the
!> last record must still begin with a number, and what follows must be the next
!> header), so a list whose records differ in width -- optional trailing columns,
!> say -- is detected rather than mis-addressed, and falls back to a bisection
!> that hunts record boundaries by newline.  Both paths return the same records.
!>
!> The file is never held in memory: a block walk reads two header lines per
!> block, a bisection reads one record per probe, and only the selected byte
!> range is read in full.  This module does no parsing beyond the leading
!> wavelength of a record, which is what the bisection compares; the callers
!> (composition_mod's init_linelist, the hop_ll_* C API) parse the slice.
!>
!> Wavelengths here are as written in the file, air Angstrom.  Selection in
!> vacuum stays with the caller, which is why a range is widened by one record on
!> either side: the air-space bisection bounds the work and never the selection.
!===============================================================================
MODULE linelist_reader_mod
  USE hop_env_mod, only: os, trace
  USE console_mod
  USE prec_mod
  implicit none
  private

  !> One species block: its header, its records' byte range, and how they are addressed.
  type, public :: ll_block_t
    character(len=80) :: id = ''      !< first header line, verbatim
    character(len=32) :: siso = ''    !< isotope token of the second header line
    integer           :: nrec = 0     !< record count the header declares
    integer(8)        :: dstart = 0   !< first byte of the first record
    integer(8)        :: dend = 0     !< one past the last byte of the block
    integer           :: width = 0    !< bytes per record incl. newline; 0 when ragged
    real(dp)          :: lam_first = 0d0, lam_last = 0d0   !< air Angstrom
  end type ll_block_t

  !> An open line-list file, positioned at the next block header.
  type, public :: ll_file_t
    integer    :: unit = 0            !< 0 when closed; newunit= hands out negatives
    integer(8) :: fsize = 0
    integer(8) :: pos = 1             !< byte offset of the next block header
    character(len=:), allocatable :: path
  end type ll_file_t

  public :: ll_open, ll_close, ll_next_block, ll_block_range, ll_whole_block
  public :: ll_read_slice, ll_read_ranges, ll_index_records

  integer, parameter :: HDR_CHUNK = 4096   !< enough for two headers and one record
  integer, parameter :: PROBE_LEN = 32     !< bytes read to parse a record's wavelength
  character, parameter :: NL = char(10)

CONTAINS

!===============================================================================
!> Open a line list for random access.  The file stays open until ll_close.
!===============================================================================
SUBROUTINE ll_open (self, path)
  type(ll_file_t), intent(out) :: self
  character(len=*), intent(in) :: path
  integer :: iostat
  !-----------------------------------------------------------------------------
  call os%check_file (path)
  self%path = trim(path)
  open (newunit=self%unit, file=self%path, access='stream', form='unformatted', &
        status='old', action='read', iostat=iostat)
  if (iostat /= 0) then
    call console_error ('linelist_reader: cannot open '//trim(self%path))
    self%unit = 0
    return
  end if
  inquire (unit=self%unit, size=self%fsize)
  self%pos = 1
END SUBROUTINE ll_open

!===============================================================================
!> Close the file.  Safe on a handle that never opened.
!===============================================================================
SUBROUTINE ll_close (self)
  type(ll_file_t), intent(inout) :: self
  !-----------------------------------------------------------------------------
  if (self%unit /= 0) close (self%unit)
  self%unit = 0
  self%pos = 1
END SUBROUTINE ll_close

!===============================================================================
!> The next block, or ok=.false. at end of file.  Reads two header lines and one
!> record; the block's extent follows from the declared record count once the
!> record width is known, and is verified before it is trusted (see the module
!> header).  A block that fails the check keeps width=0 and its end is found by
!> searching for the line that starts the next header.
!===============================================================================
SUBROUTINE ll_next_block (self, blk, ok)
  type(ll_file_t), intent(inout) :: self
  type(ll_block_t), intent(out) :: blk
  logical, intent(out) :: ok
  character(len=:), allocatable :: chunk
  integer :: n, i1, i2, i3, iostat
  integer(8) :: dend
  logical :: fixed
  !-----------------------------------------------------------------------------
  ok = .false.
  if (self%unit == 0 .or. self%pos > self%fsize) return

  n = int(min(int(HDR_CHUNK,8), self%fsize - self%pos + 1))
  allocate (character(len=n) :: chunk)
  read (self%unit, pos=self%pos) chunk

  i1 = index(chunk, NL)                              ! end of the id line
  if (i1 == 0) return
  i2 = index(chunk(i1+1:), NL)                       ! end of the iso/count line
  if (i2 == 0) return
  i2 = i1 + i2
  blk%id = chunk(1:i1-1)
  read (chunk(i1+1:i2-1), *, iostat=iostat) blk%siso, blk%nrec
  if (iostat /= 0 .or. blk%nrec < 0) return

  blk%dstart = self%pos + i2                         ! first byte after the newline
  if (blk%nrec == 0) then
    blk%dend = blk%dstart
    blk%width = 0
    self%pos = blk%dend
    ok = .true.
    return
  end if

  i3 = index(chunk(i2+1:), NL)                       ! end of the first record
  if (i3 == 0) then
    blk%width = 0
  else
    blk%width = i3                                   ! record length incl. its newline
  end if

  blk%lam_first = lam_at (self, blk%dstart)

  dend = 0
  if (blk%width > 0) then
    dend = blk%dstart + int(blk%nrec,8)*int(blk%width,8)
    call check_fixed_width (self, blk, dend, fixed)
    if (.not. fixed) blk%width = 0
  end if

  if (blk%width > 0) then
    blk%dend = dend
    blk%lam_last = lam_at (self, blk%dstart + int(blk%nrec-1,8)*int(blk%width,8))
  else
    blk%dend = ragged_block_end (self, blk%dstart)
    blk%lam_last = lam_at (self, last_record_start (self, blk%dstart, blk%dend))
  end if

  ! a block whose records descend breaks every bisection here; say so once and
  ! leave the caller to read it whole (ll_block_range then returns the full range)
  if (blk%lam_last < blk%lam_first) then
    call console_warn (ch_hop, 'linelist_reader: '//trim(blk%id)//' in '//trim(self%path)// &
             ' is not ascending in wavelength; reading it whole')
    blk%width = -1
  end if

  self%pos = blk%dend
  ok = .true.
END SUBROUTINE ll_next_block

!===============================================================================
!> Whether dstart + nrec*width really is the end of the block, and where that end
!> is.  Three endings are legitimate: another block follows, in which case the
!> computed end must land on a record boundary and on the quote that opens the
!> next header; the file ends with the block's own newline; or the file ends
!> without one, which leaves the last record a byte short of the stride and so
!> puts the computed end one past the file.  Anything else means the records are
!> not one width and the caller falls back to hunting record boundaries.  The
!> last record must also still begin with a number, which is what catches a
!> stride that happens to land on a boundary by coincidence.
!===============================================================================
SUBROUTINE check_fixed_width (self, blk, dend, ok)
  type(ll_file_t), intent(in) :: self
  type(ll_block_t), intent(in) :: blk
  integer(8), intent(inout) :: dend
  logical, intent(out) :: ok
  character(len=1) :: c
  real(dp) :: lam
  integer :: iostat
  !-----------------------------------------------------------------------------
  ok = .false.
  if (dend <= self%fsize) then
    read (self%unit, pos=dend-1) c                   ! the block's last newline
    if (c /= NL) return
    read (self%unit, pos=dend) c                     ! the next header line
    if (c /= "'") return
  else if (dend == self%fsize + 1) then
    read (self%unit, pos=self%fsize) c               ! the file ends on that newline
    if (c /= NL) return
  else if (dend == self%fsize + 2) then
    read (self%unit, pos=self%fsize) c               ! no newline after the last record
    if (c == NL) return
    dend = self%fsize + 1
  else
    return
  end if
  lam = lam_at (self, blk%dstart + int(blk%nrec-1,8)*int(blk%width,8), iostat)
  ok = iostat == 0
END SUBROUTINE check_fixed_width

!===============================================================================
!> End of a block whose records are not one width: the start of the next header
!> line, found by searching for a newline followed by a quote, else end of file.
!===============================================================================
FUNCTION ragged_block_end (self, dstart) RESULT (dend)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: dstart
  integer(8) :: dend, pos
  character(len=:), allocatable :: chunk
  integer, parameter :: SCAN_CHUNK = 1024*1024
  integer :: n, i
  !-----------------------------------------------------------------------------
  pos = dstart
  do while (pos <= self%fsize)
    n = int(min(int(SCAN_CHUNK,8), self%fsize - pos + 1))
    allocate (character(len=n) :: chunk)
    read (self%unit, pos=pos) chunk
    i = index(chunk, NL//"'")
    if (i > 0) then
      dend = pos + i                                 ! one past that newline
      deallocate (chunk)
      return
    end if
    deallocate (chunk)
    if (n < SCAN_CHUNK) exit
    pos = pos + n - 1                                ! overlap by one for a split pair
  end do
  dend = self%fsize + 1
END FUNCTION ragged_block_end

!===============================================================================
!> Start of the last record of [dstart, dend): the byte after the last newline
!> that is not the block's own trailing one.
!===============================================================================
FUNCTION last_record_start (self, dstart, dend) RESULT (start)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: dstart, dend
  integer(8) :: start, from
  character(len=:), allocatable :: chunk
  integer :: n, i, k
  !-----------------------------------------------------------------------------
  start = dstart
  from = max(dstart, dend - 512_8)
  n = int(dend - from)
  if (n <= 1) return
  allocate (character(len=n) :: chunk)
  read (self%unit, pos=from) chunk
  k = 0
  do i = n-1, 1, -1                                  ! skip the trailing newline itself
    if (chunk(i:i) == NL) then
      k = i
      exit
    end if
  end do
  if (k > 0) start = from + k
END FUNCTION last_record_start

!===============================================================================
!> The air wavelength a record starts with.
!===============================================================================
FUNCTION lam_at (self, pos, iostat) RESULT (lam)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: pos
  integer, intent(out), optional :: iostat
  real(dp) :: lam
  character(len=:), allocatable :: probe
  integer :: n, ios
  !-----------------------------------------------------------------------------
  lam = 0d0
  ios = -1
  n = int(min(int(PROBE_LEN,8), self%fsize - pos + 1))
  if (n > 0) then
    allocate (character(len=n) :: probe)
    read (self%unit, pos=pos) probe
    read (probe, *, iostat=ios) lam
  end if
  if (present(iostat)) iostat = ios
END FUNCTION lam_at

!===============================================================================
!> The whole block as a byte range, for a caller that wants every record.
!===============================================================================
SUBROUTINE ll_whole_block (blk, b0, b1)
  type(ll_block_t), intent(in) :: blk
  integer(8), intent(out) :: b0, b1
  !-----------------------------------------------------------------------------
  b0 = blk%dstart
  b1 = blk%dend
END SUBROUTINE ll_whole_block

!===============================================================================
!> The byte range of the records of blk between lam_lo and lam_hi (air Angstrom),
!> widened by one record on either side so a caller masking in vacuum still sees
!> the neighbours of its own edges.  b1 is exclusive, and b0 is clamped to
!> b_after, the end of a range the caller has already taken from this block, so
!> that windows a record apart do not hand it the same record twice.  An empty
!> selection returns b0 == b1.
!===============================================================================
SUBROUTINE ll_block_range (self, blk, lam_lo, lam_hi, b_after, b0, b1)
  type(ll_file_t), intent(in) :: self
  type(ll_block_t), intent(in) :: blk
  real(dp), intent(in) :: lam_lo, lam_hi
  integer(8), intent(in) :: b_after
  integer(8), intent(out) :: b0, b1
  integer :: k0, k1
  !-----------------------------------------------------------------------------
  if (blk%nrec == 0) then
    b0 = blk%dstart
    b1 = blk%dstart
    return
  end if

  if (blk%width < 0) then                            ! unsorted: the caller masks it all
    b0 = blk%dstart
    b1 = blk%dend
  else if (blk%width > 0) then
    k0 = bisect_fixed (self, blk, lam_lo, .false.)   ! first record with lam >= lam_lo
    k1 = bisect_fixed (self, blk, lam_hi, .true.)    ! first record with lam >  lam_hi
    k0 = max(0, k0 - 1)
    k1 = min(blk%nrec, k1 + 1)
    if (k1 <= k0) then
      b0 = blk%dstart
      b1 = blk%dstart
      return
    end if
    b0 = blk%dstart + int(k0,8)*int(blk%width,8)
    b1 = blk%dstart + int(k1,8)*int(blk%width,8)
  else
    b0 = bisect_ragged (self, blk%dstart, blk%dend, lam_lo, .false.)
    b1 = bisect_ragged (self, b0,         blk%dend, lam_hi, .true.)
    b0 = widen_back (self, blk%dstart, b0)
    b1 = widen_forward (self, b1, blk%dend)
  end if

  if (b0 < b_after) b0 = min(b_after, b1)
END SUBROUTINE ll_block_range

!===============================================================================
!> Index of the first record of a fixed-width block whose wavelength is >= lam
!> (strictly > lam with above=.true.); nrec when there is none.
!===============================================================================
FUNCTION bisect_fixed (self, blk, lam, above) RESULT (k)
  type(ll_file_t), intent(in) :: self
  type(ll_block_t), intent(in) :: blk
  real(dp), intent(in) :: lam
  logical, intent(in) :: above
  integer :: k, lo, hi, mid
  real(dp) :: probe
  !-----------------------------------------------------------------------------
  lo = 0
  hi = blk%nrec
  do while (lo < hi)
    mid = (lo + hi) / 2
    probe = lam_at (self, blk%dstart + int(mid,8)*int(blk%width,8))
    if (probe < lam .or. (above .and. probe == lam)) then
      lo = mid + 1
    else
      hi = mid
    end if
  end do
  k = lo
END FUNCTION bisect_fixed

!===============================================================================
!> Byte offset of the first record in [lo, hi) whose wavelength is >= lam
!> (strictly > lam with above=.true.), for a block that is not one width; hi when
!> there is none.  Each probe backs up from the midpoint to the record boundary.
!===============================================================================
FUNCTION bisect_ragged (self, lo_in, hi_in, lam, above) RESULT (pos)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: lo_in, hi_in
  real(dp), intent(in) :: lam
  logical, intent(in) :: above
  integer(8) :: pos, lo, hi, mid, start, nxt
  real(dp) :: probe
  !-----------------------------------------------------------------------------
  lo = lo_in
  hi = hi_in
  do while (lo < hi)
    mid = lo + (hi - lo) / 2
    start = record_start (self, lo, mid)
    probe = lam_at (self, start)
    if (probe < lam .or. (above .and. probe == lam)) then
      nxt = next_record (self, start, hi)
      if (nxt >= hi) then
        lo = hi
      else
        lo = nxt
      end if
    else
      hi = start
    end if
  end do
  pos = lo
END FUNCTION bisect_ragged

!===============================================================================
!> Start of the record containing byte pos, not before lo.
!===============================================================================
FUNCTION record_start (self, lo, pos) RESULT (start)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: lo, pos
  integer(8) :: start, from
  character(len=:), allocatable :: chunk
  integer :: n, i, k
  !-----------------------------------------------------------------------------
  start = lo
  from = max(lo, pos - 512_8)
  n = int(pos - from + 1)
  if (n <= 0) return
  allocate (character(len=n) :: chunk)
  read (self%unit, pos=from) chunk
  k = 0
  do i = n, 1, -1
    if (chunk(i:i) == NL) then
      k = i
      exit
    end if
  end do
  if (k > 0) start = from + k
END FUNCTION record_start

!===============================================================================
!> Start of the record after the one at pos, or hi when it is the last.
!===============================================================================
FUNCTION next_record (self, pos, hi) RESULT (nxt)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: pos, hi
  integer(8) :: nxt, at
  character(len=:), allocatable :: chunk
  integer :: n, i
  !-----------------------------------------------------------------------------
  at = pos
  do while (at < hi)
    n = int(min(512_8, hi - at))
    allocate (character(len=n) :: chunk)
    read (self%unit, pos=at) chunk
    i = index(chunk, NL)
    if (i > 0) then
      nxt = at + i
      deallocate (chunk)
      return
    end if
    deallocate (chunk)
    at = at + n
  end do
  nxt = hi
END FUNCTION next_record

!===============================================================================
!> One record earlier, not before lo.
!===============================================================================
FUNCTION widen_back (self, lo, pos) RESULT (out)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: lo, pos
  integer(8) :: out
  !-----------------------------------------------------------------------------
  out = lo
  if (pos > lo) out = record_start (self, lo, pos - 2_8)
END FUNCTION widen_back

!===============================================================================
!> One record later, not past hi.
!===============================================================================
FUNCTION widen_forward (self, pos, hi) RESULT (out)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: pos, hi
  integer(8) :: out
  !-----------------------------------------------------------------------------
  out = hi
  if (pos < hi) out = next_record (self, pos, hi)
END FUNCTION widen_forward

!===============================================================================
!> Read the bytes [b0, b1) into buf.
!===============================================================================
SUBROUTINE ll_read_slice (self, b0, b1, buf)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: b0, b1
  character(len=:), allocatable, intent(out) :: buf
  integer(8) :: n
  !-----------------------------------------------------------------------------
  n = range_len (self, b0, b1)
  allocate (character(len=n) :: buf)
  if (n > 0) read (self%unit, pos=b0) buf
END SUBROUTINE ll_read_slice

!===============================================================================
!> Read several byte ranges into one buffer, in the order given.  The ranges must
!> ascend and each must start and end on a record boundary, which is what
!> ll_block_range returns: the result is then a valid, still ascending run of
!> records with the bytes between the ranges left out, and a caller can treat a
!> set of windows as one slice of the block.
!===============================================================================
SUBROUTINE ll_read_ranges (self, r0, r1, nr, buf)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: r0(:), r1(:)
  integer, intent(in) :: nr
  character(len=:), allocatable, intent(out) :: buf
  integer(8) :: total, at, n
  integer :: k
  !-----------------------------------------------------------------------------
  total = 0
  do k = 1, nr
    total = total + range_len (self, r0(k), r1(k))
  end do
  allocate (character(len=total) :: buf)
  at = 1
  do k = 1, nr
    n = range_len (self, r0(k), r1(k))
    if (n > 0) then
      read (self%unit, pos=r0(k)) buf(at:at+n-1)
      at = at + n
    end if
  end do
END SUBROUTINE ll_read_ranges

!===============================================================================
!> Bytes of [b0, b1) that lie in the file.
!===============================================================================
FUNCTION range_len (self, b0, b1) RESULT (n)
  type(ll_file_t), intent(in) :: self
  integer(8), intent(in) :: b0, b1
  integer(8) :: n
  !-----------------------------------------------------------------------------
  n = max(0_8, min(b1, self%fsize + 1) - b0)
END FUNCTION range_len

!===============================================================================
!> Record boundaries within buf, as byte offsets into it: lstart(k):lend(k) is
!> record k without its newline (and without a carriage return before it).
!===============================================================================
SUBROUTINE ll_index_records (buf, lstart, lend, n)
  character(len=*), intent(in) :: buf
  integer(8), allocatable, intent(inout) :: lstart(:), lend(:)
  integer, intent(out) :: n
  integer(8) :: b0, bnl, e, fsize
  integer :: rel
  !-----------------------------------------------------------------------------
  fsize = len(buf)
  if (.not. allocated(lstart)) allocate (lstart(1024), lend(1024))
  n = 0
  b0 = 1
  do
    if (b0 > fsize) exit
    rel = index(buf(b0:fsize), NL)
    if (rel == 0) then
      bnl = fsize + 1                                ! last record, no newline
    else
      bnl = b0 + rel - 1
    end if
    n = n + 1
    if (n > size(lstart)) then
      call grow_i8 (lstart)
      call grow_i8 (lend)
    end if
    lstart(n) = b0
    e = bnl - 1
    if (e >= b0) then
      if (buf(e:e) == char(13)) e = e - 1            ! drop CR (CRLF files)
    end if
    lend(n) = e
    b0 = bnl + 1
  end do

CONTAINS
  subroutine grow_i8 (a)
    integer(8), allocatable, intent(inout) :: a(:)
    integer(8), allocatable :: t(:)
    allocate (t(2_8*size(a)))
    t(1:size(a)) = a
    call move_alloc (t, a)
  end subroutine grow_i8
END SUBROUTINE ll_index_records

END MODULE linelist_reader_mod
