!===============================================================================
!> Reports a namelist group that failed to read and stops the run.
!>
!> Each key given in the input that the group does not declare is listed, with
!> the closest declared key when the edit distance makes a typo likely, and the
!> declared keys are listed when some unknown key has no such match.  An error
!> without unknown keys (a malformed value, say) is reported with the compiler
!> message.
!>
!> A namelist group cannot be passed to a procedure, so the caller writes the
!> group to a scratch unit and the declared keys are parsed from that text; the
!> READ_NAMELIST macro in read_namelist.inc does this.
!===============================================================================
MODULE namelist_check_mod
  USE iso_fortran_env, only: error_unit, iostat_eor
  USE console_mod, only: console_error
  implicit none
  private
  integer, parameter:: name_len=64
  public:: namelist_error
CONTAINS

!===============================================================================
!> Report the failed read of a group from unit input and stop.  scratch holds the
!> group as written by write(scratch, nml=group) and is closed here.
!===============================================================================
SUBROUTINE namelist_error (input, scratch, iomsg)
  integer, intent(in):: input, scratch
  character(len=*), intent(in):: iomsg
  character(len=:), allocatable:: declared_text, given_text, group
  character(len=name_len), allocatable:: declared(:), given(:)
  character(len=256):: hint
  logical:: list_declared
  integer:: i, n_unknown
  !-----------------------------------------------------------------------------
  call read_unit (scratch, declared_text)
  close (scratch)
  allocate (group, source=trim(group_name(declared_text)))
  call collect_keys (declared_text, declared)
  call group_text (input, group, given_text)
  call collect_keys (given_text, given)
  list_declared = .false.
  n_unknown = 0
  do i=1, size(given)
    if (any(declared == given(i))) cycle
    n_unknown = n_unknown + 1
    hint = suggestion(given(i), declared)
    if (hint == '') then
      call console_error ('unknown key '//trim(given(i))//' in &'//group)
      list_declared = .true.
    else
      call console_error ('unknown key '//trim(given(i))//' in &'//group// &
        ' (did you mean '//trim(hint)//'?)')
    end if
  end do
  if (n_unknown == 0) call console_error ('cannot read &'//group//': '//trim(iomsg))
  if (list_declared) call write_keys (group, declared)
  error stop 'namelist error'
END SUBROUTINE namelist_error

!===============================================================================
!> All records of a unit from the start, joined by blanks
!===============================================================================
SUBROUTINE read_unit (unit, text)
  integer, intent(in):: unit
  character(len=:), allocatable, intent(out):: text
  character(len=:), allocatable:: record
  integer:: iostat
  !-----------------------------------------------------------------------------
  allocate (character(len=0):: text)
  rewind (unit)
  do
    call read_record (unit, record, iostat)
    if (iostat /= 0) exit
    call append (text, ' '//record)
  end do
END SUBROUTINE read_unit

!===============================================================================
!> One record of any length; iostat is nonzero at the end of the file
!===============================================================================
SUBROUTINE read_record (unit, record, iostat)
  integer, intent(in):: unit
  character(len=:), allocatable, intent(out):: record
  integer, intent(out):: iostat
  character(len=256):: chunk
  integer:: n
  !-----------------------------------------------------------------------------
  allocate (character(len=0):: record)
  do
    read (unit, '(a)', advance='no', size=n, iostat=iostat) chunk
    call append (record, chunk(1:n))
    if (iostat == iostat_eor) then
      iostat = 0
      exit
    end if
    if (iostat /= 0) exit
  end do
END SUBROUTINE read_record

!===============================================================================
!> The input text of a group, from after '&group' to its terminating '/', with
!> comments removed; empty when the group is not found
!===============================================================================
SUBROUTINE group_text (unit, group, text)
  integer, intent(in):: unit
  character(len=*), intent(in):: group
  character(len=:), allocatable, intent(out):: text
  character(len=:), allocatable:: record
  character:: quote
  integer:: iostat, i, i0, n
  logical:: inside
  !-----------------------------------------------------------------------------
  allocate (character(len=0):: text)
  inside = .false.
  quote = ' '
  n = len(group) + 1
  rewind (unit)
  do
    call read_record (unit, record, iostat)
    if (iostat /= 0) exit
    i0 = 1
    if (.not. inside) then
      record = adjustl(record)
      if (lower(record(1:min(n,len(record)))) /= '&'//group) cycle
      if (len(record) > n) then
        if (is_name_char(record(n+1:n+1))) cycle
      end if
      inside = .true.
      i0 = n + 1
    end if
    do i=i0, len(record)
      if (quote /= ' ') then
        if (record(i:i) == quote) quote = ' '
      else if (record(i:i) == '"' .or. record(i:i) == "'") then
        quote = record(i:i)
      else if (record(i:i) == '!') then
        exit
      else if (record(i:i) == '/') then
        call append (text, ' '//record(i0:i-1))
        return
      end if
    end do
    call append (text, ' '//record(i0:i-1))
  end do
END SUBROUTINE group_text

!===============================================================================
!> The lower-case group name after the leading '&'
!===============================================================================
FUNCTION group_name (text) result(name)
  character(len=*), intent(in):: text
  character(len=name_len):: name
  integer:: i0, i1
  !-----------------------------------------------------------------------------
  i0 = index(text, '&') + 1
  i1 = i0
  do while (i1 <= len(text))
    if (.not. is_name_char(text(i1:i1))) exit
    i1 = i1 + 1
  end do
  name = lower(text(i0:i1-1))
END FUNCTION group_name

!===============================================================================
!> The distinct lower-case names that precede '=' outside quoted strings, without
!> subscripts or derived-type components
!===============================================================================
SUBROUTINE collect_keys (text, keys)
  character(len=*), intent(in):: text
  character(len=name_len), allocatable, intent(out):: keys(:)
  character(len=name_len), allocatable:: grown(:)
  character(len=name_len):: key
  character:: quote
  integer:: i, j, last, depth, n
  !-----------------------------------------------------------------------------
  allocate (keys(8))
  n = 0
  quote = ' '
  do i=1, len(text)
    if (quote /= ' ') then
      if (text(i:i) == quote) quote = ' '
      cycle
    end if
    if (text(i:i) == '"' .or. text(i:i) == "'") quote = text(i:i)
    if (text(i:i) /= '=') cycle
    j = skip_blanks_back(text, i-1)
    if (j > 0) then
      if (text(j:j) == ')') then
        depth = 0
        do while (j > 0)
          if (text(j:j) == ')') depth = depth + 1
          if (text(j:j) == '(') depth = depth - 1
          if (depth == 0) exit
          j = j - 1
        end do
        j = skip_blanks_back(text, j-1)
      end if
    end if
    last = j
    do while (j > 0)
      if (.not. (is_name_char(text(j:j)) .or. text(j:j) == '%')) exit
      j = j - 1
    end do
    if (j == last) cycle
    key = lower(text(j+1:last))
    if (index(key, '%') > 0) key = key(1:index(key, '%')-1)
    if (any(keys(1:n) == key)) cycle
    if (n == size(keys)) then
      allocate (grown(2*n))
      grown(1:n) = keys
      call move_alloc (grown, keys)
    end if
    n = n + 1
    keys(n) = key
  end do
  allocate (grown(n))
  grown = keys(1:n)
  call move_alloc (grown, keys)
END SUBROUTINE collect_keys

!===============================================================================
!> The declared keys closest to an unknown key, or '' when none is close enough.
!> The threshold grows with the key length: 1 edit up to 5 characters, 2 up to
!> 8, 3 beyond.
!===============================================================================
FUNCTION suggestion (key, declared) result(hint)
  character(len=*), intent(in):: key
  character(len=name_len), intent(in):: declared(:)
  character(len=256):: hint
  integer:: i, best, d(size(declared))
  !-----------------------------------------------------------------------------
  hint = ''
  if (size(declared) == 0) return
  do i=1, size(declared)
    d(i) = edit_distance(trim(key), trim(declared(i)))
  end do
  best = minval(d)
  if (best > max(1, min(3, len_trim(key)/3))) return
  do i=1, size(declared)
    if (d(i) /= best) cycle
    if (hint /= '') hint = trim(hint)//' or '
    hint = trim(hint)//trim(declared(i))
  end do
END FUNCTION suggestion

!===============================================================================
!> Damerau-Levenshtein distance in its optimal-string-alignment form: insertions,
!> deletions, substitutions and swaps of adjacent characters each count as one
!===============================================================================
INTEGER FUNCTION edit_distance (a, b) result(dist)
  character(len=*), intent(in):: a, b
  integer:: d(0:len(a),0:len(b)), i, j, cost
  !-----------------------------------------------------------------------------
  do i=0, len(a)
    d(i,0) = i
  end do
  do j=0, len(b)
    d(0,j) = j
  end do
  do j=1, len(b)
    do i=1, len(a)
      cost = merge(0, 1, a(i:i) == b(j:j))
      d(i,j) = min(d(i-1,j) + 1, d(i,j-1) + 1, d(i-1,j-1) + cost)
      if (i > 1 .and. j > 1) then
        if (a(i:i) == b(j-1:j-1) .and. a(i-1:i-1) == b(j:j)) &
          d(i,j) = min(d(i,j), d(i-2,j-2) + 1)
      end if
    end do
  end do
  dist = d(len(a),len(b))
END FUNCTION edit_distance

!===============================================================================
!> The declared keys, wrapped to lines of about 80 characters
!===============================================================================
SUBROUTINE write_keys (group, keys)
  character(len=*), intent(in):: group
  character(len=name_len), intent(in):: keys(:)
  character(len=160):: line
  integer:: i, n
  !-----------------------------------------------------------------------------
  call console_error ('valid keys in &'//group//':')
  line = ''
  n = 1
  do i=1, size(keys)
    if (n + len_trim(keys(i)) > 78) then
      write (error_unit, '(a)') trim(line)
      line = ''
      n = 1
    end if
    line(n+1:) = ' '//trim(keys(i))
    n = n + 1 + len_trim(keys(i))
  end do
  write (error_unit, '(a)') trim(line)
  flush (error_unit)
END SUBROUTINE write_keys

!===============================================================================
!> Append piece to an allocated string; growing by assignment would need
!> reallocation on assignment, which the builds disable
!===============================================================================
SUBROUTINE append (text, piece)
  character(len=:), allocatable, intent(inout):: text
  character(len=*), intent(in):: piece
  character(len=:), allocatable:: grown
  !-----------------------------------------------------------------------------
  allocate (character(len=len(text)+len(piece)):: grown)
  grown(1:len(text)) = text
  grown(len(text)+1:) = piece
  call move_alloc (grown, text)
END SUBROUTINE append

!===============================================================================
INTEGER FUNCTION skip_blanks_back (text, i) result(j)
  character(len=*), intent(in):: text
  integer, intent(in):: i
  !-----------------------------------------------------------------------------
  j = i
  do while (j > 0)
    if (text(j:j) /= ' ' .and. text(j:j) /= achar(9)) exit
    j = j - 1
  end do
END FUNCTION skip_blanks_back

!===============================================================================
LOGICAL FUNCTION is_name_char (c)
  character, intent(in):: c
  !-----------------------------------------------------------------------------
  is_name_char = (c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z') .or. &
                 (c >= '0' .and. c <= '9') .or. c == '_'
END FUNCTION is_name_char

!===============================================================================
FUNCTION lower (s) result(t)
  character(len=*), intent(in):: s
  character(len=len(s)):: t
  integer:: i
  !-----------------------------------------------------------------------------
  t = s
  do i=1, len(s)
    if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') t(i:i) = achar(iachar(s(i:i)) + 32)
  end do
END FUNCTION lower

END MODULE namelist_check_mod
