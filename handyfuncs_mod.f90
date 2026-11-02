MODULE handyfuncs_mod
  USE op_math_mod
  USE hop_env_mod, only: os
  USE console_mod
  implicit none

  interface truncate_array
    module procedure truncate_string_array, truncate_float_array, &
      truncate_dble_array, truncate_int_array, truncate_bool_array
  end interface truncate_array

  interface check_array
    module procedure check_1d_array_dp, check_2d_array_dp, check_3d_array_dp, check_4d_array_dp
    module procedure check_1d_array_sp, check_2d_array_sp, check_3d_array_sp, check_4d_array_sp
  end interface check_array

contains

  !====================================================================================
  SUBROUTINE print_failed_line(lunit, iostat, text)
    integer, intent(in) :: lunit, iostat
    character(len=512) :: line_text
    character(len=4096) :: filename
    character(len=*), intent(in), optional :: text
    ! ----------------------------------------------------------------------
    if (iostat > 0) then
      if (present(text)) then
        read(text,'(A)') line_text
      else
        inquire(unit=lunit, name=filename)
        backspace lunit
        read(lunit,'(A)') line_text
      endif
      write(filename,'(a,i0,2a)') 'failed to read a line (iostat=', iostat, ') in ', trim(filename)
      call console_error (trim(filename)//': '//trim(line_text))
      error stop
    endif
  END SUBROUTINE print_failed_line

  !====================================================================================
  SUBROUTINE NextEntry(fileunit, comment)
    ! ***********************************************************************
    ! --- Reads an input line, returns when it is not a
    ! --- comment or an empty line.
    ! --- Feb 21 2005:
    ! ---   Created by Jorrit Leenaarts (j.leenaarts@astro.uu.nl)
    ! --- May 28 2008 :
    !       jorritl: updated, will now also check for empty lines,
    !       and allow data in the first column of the file.
    ! ***********************************************************************
    integer,intent(in) :: fileunit
    character(len=132) :: instring
    character(len=1)   :: token1, cmt
    character(len=*),intent(in), optional :: comment
    !---------------------------------------------------------------------------------
    if (present(comment)) then
      cmt = comment
    else
      cmt = '*'
    endif
    do
      ! read the line
      read(unit=fileunit,fmt='(A72)') instring
      instring=adjustl(instring)
      token1=instring(1:1)

      ! exit if line is not a comment and not empty
      if (token1 .ne. cmt .and. token1 .ne. ' ') then
        backspace(fileunit)
        exit
      endif
    end do
  END SUBROUTINE NextEntry

  !====================================================================================
  SUBROUTINE unequal_intervals(n, i, ni, start, end)
    integer, intent(in)  :: n, i, ni
    integer, intent(out) :: start, end
    integer :: nv, md, j
    !---------------------------------------------------------------------------
    nv    = n / ni
    md    = mod(n, ni)
    j     = mod(abs(i)-1, ni) + 1
    start = (j-1)*nv + min(md, j-1) + 1
    end   =  j   *nv + min(md, j)
  END SUBROUTINE unequal_intervals

  !===============================================================================
  SUBROUTINE replace_char (string, out, in)
    character(len=*), intent(inout) :: string
    character, intent(in) :: out, in
    integer :: i, n
    !-----------------------------------------------------------------------------
    n = len(string)
    do i=1, n
      if (string(i:i) == out) then
        string(i:i) = in
      endif
    end do
  END SUBROUTINE replace_char

  !===============================================================================
  SUBROUTINE replace_substring (string, out, in)
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: out, in
    character(:), allocatable :: tmp
    integer :: i, j, n, no, ni
    !-----------------------------------------------------------------------------
    n = len(string)
    no = len(out)
    ni = len(in)
    allocate(tmp, mold=string)
    i = 1
    j = 1
    do while (i<=n-no+1)
      if (string(i:i+no-1) == out) then
        tmp(j:j+ni-1) = in
        j = j + ni
        i = i + no
      else
        tmp(j:j) = string(i:i)
        i = i + 1
        j = j + 1
      endif
    end do
    string = tmp
  END SUBROUTINE replace_substring

  !===============================================================================
  SUBROUTINE count_words (string, nwords, sep)
    character(len=*), intent(in) :: string
    integer, intent(out) :: nwords
    character, intent(in), optional :: sep
    character :: seperator
    integer :: i, n
    logical :: new_word
    !-----------------------------------------------------------------------------
    if (present(sep)) then
      seperator = sep
    else
      seperator = ' '
    endif

    n = len(string)
    nwords = 0
    new_word = .true.
    do i=1, n
      if (string(i:i) == seperator) then
        new_word = .true.
      elseif(new_word) then
        nwords = nwords + 1
        new_word = .false.
      endif
    end do
  END SUBROUTINE count_words

  !===============================================================================
  FUNCTION trim_between_words(string) result(res)
    character(len=*) :: string
    character(len=len(string)) :: res
    character :: char
    integer :: i, j
    !---------------------------------------------------------------------------
    res = ''
    char = ' '
    j = 1
    call replace_char(string, "'", ' ')
    call replace_char(string, ",", ' ')
    do i=1, len(string)
      if (string(i:i) == ' ') then
        if (char == ' ') then
          cycle
        else
          res(j:j) = string(i:i)
          j = j + 1
        endif
      else
        res(j:j) = string(i:i)
        j = j + 1
      endif
      char = string(i:i)
    end do
  END FUNCTION

  !===============================================================================
  FUNCTION string_in_text (text, str) result(in_text)
    character(len=*) :: str
    character(len=*) :: text
    character(len=256) :: dummy
    integer :: i, j
    logical :: in_text
    ! ----------------------------------------------------------------------------
    in_text = .false.
    dummy = ''
    j = 0
    do i=1, len(text)
      if (text(i:i) == ' ') then
        if (dummy(1:j) == str) exit
        j = 0
      else
        j = j + 1
        dummy(j:j) = text(i:i)
      endif
    end do
    if (dummy(1:j) == str) in_text = .true.
  END FUNCTION

  !===============================================================================
  FUNCTION extract_int (str) result(val)
    character(len=*) :: str
    integer :: i, o, val
    ! ----------------------------------------------------------------------------
    i = scan(str, '0123456789')
    o = scan(str, '0123456789', back=.true.)
    read(str(i:o), *) val
  END FUNCTION

    !===============================================================================
  FUNCTION hash_string(text) result(hashed)
    character(*), intent(in) :: text
    integer :: hashed ! 32-bit integers assumed
    integer, parameter :: magic_numb = int(Z'5D7A9F43', kind=kind(magic_numb))  ! hexadecimal
    integer :: i, j
    !-----------------------------------------------------------------------------
    hashed = 0
    do i = 1, len_trim(text)
      j = mod(i-1, 4) * 8
      hashed = IEOR( hashed, ISHFT( ICHAR( text(i:i) ), j ) )
    end do
    hashed = abs( IEOR( hashed, magic_numb ) )
  END FUNCTION hash_string

  !===============================================================================
  SUBROUTINE argsort_int(a, order)
    integer, intent(in) :: a(:)
    integer, intent(out) :: order(:)
    integer, allocatable :: scr(:)
    integer :: i, j, n
    !-----------------------------------------------------------------------------
    allocate(scr, source=a)
    n = size(a, 1)
    order = 1

    do i = 2, n
      j = i - 1
      do while (j >= 1)
        if (scr(j) <= a(i)) exit
        scr(j + 1) = scr(j)
        order(j+1) = order(j)
        j = j - 1
      end do
      scr(j + 1) = a(i)
      order(j+1) = i
    end do
  END SUBROUTINE argsort_int

  !===============================================================================
  ! In-place ascending sort (recursive quicksort, middle-element pivot).  Both call
  ! sites (spectrum_mod%sort, composition_mod) are one-time setup sorts, not hot
  ! loops.  The middle-element pivot keeps the already-sorted and
  ! concatenated-sorted-runs inputs they pass off the O(n^2) path.
  RECURSIVE SUBROUTINE sort_asc(a)
    real(8), intent(inout) :: a(:)
    !-----------------------------------------------------------------------------
    call qs(1, size(a))
  contains
    recursive subroutine qs(lo, hi)
      integer, intent(in) :: lo, hi
      integer  :: i, j
      real(8)  :: pivot, t
      if (lo >= hi) return
      pivot = a((lo+hi)/2)
      i = lo; j = hi
      do
        do while (a(i) < pivot); i = i + 1; end do
        do while (a(j) > pivot); j = j - 1; end do
        if (i <= j) then
          t = a(i); a(i) = a(j); a(j) = t
          i = i + 1; j = j - 1
        end if
        if (i > j) exit
      end do
      if (lo < j) call qs(lo, j)
      if (i < hi) call qs(i, hi)
    end subroutine qs
  END SUBROUTINE sort_asc

  !====================================================================================
  PURE FUNCTION lowercase(string)
    character(len=*),intent(in) :: string
    character(len=len(string)) :: lowercase
    character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
    integer :: i,k
    !---------------------------------------------------------------------------------
    lowercase=' '
    lowercase=string
    do i=1, len(string)
      k = index(ucase,string(i:i))
      if (k > 0) then
        lowercase(i:i)=lcase(k:k)
      else
        lowercase(i:i)=string(i:i)
      end if
    end do
  END FUNCTION lowercase

  !====================================================================================
  PURE FUNCTION uppercase(string)
    ! converts the first ss characters of string to uppercase and returns them as a new string.
    implicit none
    character(len=*),intent(in) :: string
    character(len=len(string)) :: uppercase
    character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
    integer :: i,k
    !---------------------------------------------------------------------------------
    uppercase=' '
    uppercase=string
    do i=1, len(string)
      k=index(lcase,string(i:i))
      if (k > 0) then
        uppercase(i:i)=ucase(k:k)
      else
        uppercase(i:i)=string(i:i)
      end if
    end do
  END FUNCTION uppercase

  !====================================================================================
  PURE FUNCTION capitalize(string)
    character(len=*),intent(in) :: string
    character(len=len(string)) :: capitalize
    !---------------------------------------------------------------------------------
    capitalize = lowercase(trim(string))
    capitalize(1:1) = uppercase(capitalize(1:1))
  END FUNCTION capitalize

  !===============================================================================
  SUBROUTINE remove_close_entries(array, threshold)
    real(8), intent(inout), allocatable :: array(:)
    real(8), intent(in) :: threshold
    real(8), allocatable :: scratch(:)
    integer :: i, new_size
    !-----------------------------------------------------------------------------
    new_size = size(array)
    i = 1
    do while (i < new_size)
      if (abs(array(i+1) - array(i)) < threshold) then
        array(i+1:new_size-1) = array(i+2:new_size)
        new_size = new_size - 1
      else
        i = i + 1
      end if
    end do

    allocate(scratch, source=array(1:new_size))
    deallocate(array)
    allocate(array, source=scratch)
    deallocate(scratch)
  END SUBROUTINE remove_close_entries

  !===============================================================================
  elemental FUNCTION vac2air (vac) result(air)
    real(8), intent(in) :: vac
    real(8) :: air
    !-----------------------------------------------------------------------------
    air = vac/airvac_corr(vac)
  END FUNCTION vac2air

  !===============================================================================
  elemental FUNCTION air2vac (air) result(vac)
    real(8), intent(in) :: air
    real(8) :: vac, w
    !-----------------------------------------------------------------------------
    w = air*airvac_corr(air)
    w = air*airvac_corr(w)
    vac = air*airvac_corr(w)
  END FUNCTION air2vac

  !===============================================================================
  elemental FUNCTION airvac_corr (air) result(n)
    real(8), intent(in) :: air
    real(8) :: n, s, s2
    !-----------------------------------------------------------------------------
    s = 1d4 / air
    s2 = s*s
    n = 1.00008336624212083d0 + 0.02408926869968d0 / (130.1065924522d0 - s2) &
      + 0.0001599740894897d0 / (38.92568793293d0 - s2)
  END FUNCTION airvac_corr

  !===============================================================================
  FUNCTION rindex(text, str) result(index)
    character(len=*), intent(in) :: text, str
    integer :: index
    integer :: i, ntext, nstr
    !-----------------------------------------------------------------------------
    ntext = len_trim(text)
    nstr = len_trim(str)
    index = 0

    if (nstr > ntext) return

    do i = ntext - nstr + 1, 1, -1
      if (text(i:i+nstr-1) == str) then
        index = i
        exit
      end if
    end do
  END FUNCTION rindex

  !===============================================================================
  FUNCTION startswith(text, prefix)
    character(len=*), intent(in) :: text, prefix
    logical :: startswith
    integer :: text_len, prefix_len
    !-----------------------------------------------------------------------------
    text_len = len_trim(text)
    prefix_len = len_trim(prefix)
    if (prefix_len > text_len) then
      startswith = .false.
    else
      startswith = (text(1:prefix_len) == prefix)
    end if
  END FUNCTION startswith

  !===============================================================================
  FUNCTION endswith(text, suffix)
    character(len=*), intent(in) :: text, suffix
    logical :: endswith
    integer :: text_len, suffix_len
    !-----------------------------------------------------------------------------
    text_len = len_trim(text)
    suffix_len = len_trim(suffix)
    if (suffix_len > text_len) then
      endswith = .false.
    else
      endswith = (text(text_len - suffix_len + 1:text_len) == suffix)
    end if
  END FUNCTION endswith

  !====================================================================================
  SUBROUTINE truncate_string_array(array, nkeep)
    character(len=*), intent(inout), dimension(:), allocatable :: array
    character(len=len(array(1))), dimension(:), allocatable :: tmp
    integer, intent(in) :: nkeep
    ! ----------------------------------------------------------------------
    allocate(tmp, source=array(1:nkeep))
    deallocate(array)
    call move_alloc(tmp, array)
  END SUBROUTINE truncate_string_array

  !====================================================================================
  SUBROUTINE truncate_float_array(array, nkeep)
    real(4), intent(inout), dimension(:), allocatable :: array
    real(4), dimension(:), allocatable :: tmp
    integer, intent(in) :: nkeep
    ! ----------------------------------------------------------------------
    allocate(tmp(nkeep), source=array(1:nkeep))
    deallocate(array)
    call move_alloc(tmp, array)
  END SUBROUTINE truncate_float_array

  !====================================================================================
  SUBROUTINE truncate_dble_array(array, nkeep)
    real(8), intent(inout), dimension(:), allocatable :: array
    real(8), dimension(:), allocatable :: tmp
    integer, intent(in) :: nkeep
    ! ----------------------------------------------------------------------
    allocate(tmp(nkeep), source=array(1:nkeep))
    deallocate(array)
    call move_alloc(tmp, array)
  END SUBROUTINE truncate_dble_array

  !====================================================================================
  SUBROUTINE truncate_int_array(array, nkeep)
    integer, intent(inout), dimension(:), allocatable :: array
    integer, dimension(:), allocatable :: tmp
    integer, intent(in) :: nkeep
    ! ----------------------------------------------------------------------
    allocate(tmp(nkeep), source=array(1:nkeep))
    deallocate(array)
    call move_alloc(tmp, array)
  END SUBROUTINE truncate_int_array

  !====================================================================================
  SUBROUTINE truncate_bool_array(array, nkeep)
    logical, intent(inout), dimension(:), allocatable :: array
    logical, dimension(:), allocatable :: tmp
    integer, intent(in) :: nkeep
    ! ----------------------------------------------------------------------
    allocate(tmp(nkeep), source=array(1:nkeep))
    deallocate(array)
    call move_alloc(tmp, array)
  END SUBROUTINE truncate_bool_array

  ! ==============================================================================
  SUBROUTINE report_invalid (what, label, idx)
    character(len=*), intent(in) :: what, label
    integer, intent(in) :: idx(:)
    character(len=256) :: msg
    ! ----------------------------------------------------------------------
    write(msg,'(4a,*(1x,i0))') what, ' in ', trim(label), ' at position', idx
    call console_error (msg)
  END SUBROUTINE report_invalid

  SUBROUTINE check_array_core_dp(arr, arr_shape, label, abort, check_zero, check_negative, ok)
    real(8), intent(in), dimension(:) :: arr
    integer, intent(in), dimension(:) :: arr_shape
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    
    logical :: has_issue, do_check_zero, do_check_negative
    integer :: i, d, remainder
    integer, dimension(size(arr_shape)) :: idx
    ! ----------------------------------------------------------------------
    has_issue = .false.
    do_check_zero = .false.
    do_check_negative = .false.
    if (present(check_zero)) do_check_zero = check_zero
    if (present(check_negative)) do_check_negative = check_negative
    do i = 1, size(arr)
      if (arr(i) /= arr(i) .or. arr(i) > inf8 &
      .or. (do_check_zero .and. arr(i) == 0.d0) &
      .or. (do_check_negative .and. arr(i) < 0.d0)) &
      then
        has_issue = .true.
        remainder = i - 1
        do d = 1, size(arr_shape)
          idx(d) = mod(remainder, arr_shape(d)) + 1
          remainder = remainder / arr_shape(d)
        end do

        if (arr(i) /= arr(i)) then
          call report_invalid ('NaN', label, idx)
        elseif (arr(i) > inf8) then
          call report_invalid ('Infinity', label, idx)
        elseif (do_check_zero .and. arr(i) == 0.d0) then
          call report_invalid ('Zero', label, idx)
        elseif (do_check_negative .and. arr(i) < 0.d0) then
          call report_invalid ('Negative value', label, idx)
        endif
      end if
    end do

    if (present(ok)) ok = .not. has_issue
    if (has_issue .and. abort) then
      call console_error ('aborting due to invalid values in '//trim(label))
      error stop
    end if
  END SUBROUTINE check_array_core_dp

  SUBROUTINE check_array_core_sp(arr, arr_shape, label, abort, check_zero, check_negative, ok)
    real(4), intent(in), dimension(:) :: arr
    integer, intent(in), dimension(:) :: arr_shape
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    
    logical :: has_issue, do_check_zero, do_check_negative
    integer :: i, d, remainder
    integer, dimension(size(arr_shape)) :: idx
    ! ----------------------------------------------------------------------
    has_issue = .false.
    do_check_zero = .false.
    do_check_negative = .false.
    if (present(check_zero)) do_check_zero = check_zero
    if (present(check_negative)) do_check_negative = check_negative
    do i = 1, size(arr)
      if (arr(i) /= arr(i) .or. arr(i) > inf4 &
        .or. (do_check_zero .and. arr(i) == 0.0) &
        .or. (do_check_negative .and. arr(i) < 0.0)) then
        has_issue = .true.
        remainder = i - 1
        do d = 1, size(arr_shape)
          idx(d) = mod(remainder, arr_shape(d)) + 1
          remainder = remainder / arr_shape(d)
        end do

        if (arr(i) /= arr(i)) then
          call report_invalid ('NaN', label, idx)
        elseif (arr(i) > inf4) then
          call report_invalid ('Infinity', label, idx)
        elseif (do_check_zero .and. arr(i) == 0.0) then
          call report_invalid ('Zero', label, idx)
        elseif (do_check_negative .and. arr(i) < 0.0) then
          call report_invalid ('Negative value', label, idx)
        endif
      end if
    end do

    if (present(ok)) ok = .not. has_issue
    if (has_issue .and. abort) then
      call console_error ('aborting due to invalid values in '//trim(label))
      error stop
    end if
  END SUBROUTINE check_array_core_sp

  ! ==============================================================================
  SUBROUTINE check_1d_array_dp(arr, label, abort, check_zero, check_negative, ok)
    real(8), intent(in), dimension(:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_dp(arr, shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_1d_array_dp

  ! ==============================================================================
  SUBROUTINE check_2d_array_dp(arr, label, abort, check_zero, check_negative, ok)
    real(8), intent(in), dimension(:,:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_dp(reshape(arr, [size(arr)]), shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_2d_array_dp

  ! ==============================================================================
  SUBROUTINE check_3d_array_dp(arr, label, abort, check_zero, check_negative, ok)
    real(8), intent(in), dimension(:,:,:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_dp(reshape(arr, [size(arr)]), shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_3d_array_dp

  ! ==============================================================================
  SUBROUTINE check_4d_array_dp(arr, label, abort, check_zero, check_negative, ok)
    real(8), intent(in), dimension(:,:,:,:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_dp(reshape(arr, [size(arr)]), shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_4d_array_dp

  ! ==============================================================================
  SUBROUTINE check_1d_array_sp(arr, label, abort, check_zero, check_negative, ok)
    real(4), intent(in), dimension(:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_sp(arr, shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_1d_array_sp

  ! ==============================================================================
  SUBROUTINE check_2d_array_sp(arr, label, abort, check_zero, check_negative, ok)
    real(4), intent(in), dimension(:,:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_sp(reshape(arr, [size(arr)]), shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_2d_array_sp

  ! ==============================================================================
  SUBROUTINE check_3d_array_sp(arr, label, abort, check_zero, check_negative, ok)
    real(4), intent(in), dimension(:,:,:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_sp(reshape(arr, [size(arr)]), shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_3d_array_sp

  ! ==============================================================================
  SUBROUTINE check_4d_array_sp(arr, label, abort, check_zero, check_negative, ok)
    real(4), intent(in), dimension(:,:,:,:) :: arr
    character(len=*), intent(in) :: label
    logical, intent(in) :: abort
    logical, intent(in), optional :: check_zero, check_negative
    logical, optional, intent(out) :: ok
    ! ----------------------------------------------------------------------
    call check_array_core_sp(reshape(arr, [size(arr)]), shape(arr), label, abort, check_zero, check_negative, ok)
  END SUBROUTINE check_4d_array_sp


  ! ==============================================================================
  PURE SUBROUTINE symmetric_log_grid(q0, qmax, q)
    ! Build a symmetric log-spaced grid with preferred-number mantissas
    ! (1, 1.5, 2, 3, 4, 6, 8) per decade, mirrored about zero.
    ! The endpoints are snapped to the nearest preferred mantissa.
    ! Output q has the shape  [-q_n, ..., -q_1, 0, q_1, ..., q_n].
    ! ----------------------------------------------------------------------
    real(8),              intent(in)  :: q0, qmax
    real(8), allocatable, intent(out) :: q(:)

    real(8), parameter :: mantissas(*) = &
        [1.0_8, 1.5_8, 2.0_8, 3.0_8, 4.0_8, 6.0_8, 8.0_8]
    integer,  parameter :: nref = size(mantissas)

    real(8) :: ref(nref)
    real(8) :: lq0, lqm, r0, rm, d0, dm
    integer  :: istart, iend, d, n_pos, n_total, k
    real(8), allocatable :: pos(:)

    ! Reference mantissas in log space
    ref = log10(mantissas)

    ! Decade / remainder split.
    lq0 = log10(q0)
    d0  = floor(lq0)
    r0  = lq0 - d0

    lqm = log10(qmax)
    dm  = floor(lqm)
    rm  = lqm - dm

    ! Closest mantissa indices.
    istart = minloc(abs(ref - r0), dim=1)
    iend   = minloc(abs(ref - rm), dim=1)

    ! Build the positive side.
    if (nint(d0) == nint(dm)) then
        n_pos = iend - istart + 1
        allocate(pos(n_pos))
        pos = 10.0_8**(d0 + ref(istart:iend))
    else
        n_pos = (nref - istart + 1) &
              + (nint(dm) - nint(d0) - 1)*nref &
              + iend
        allocate(pos(n_pos))

        k = 0
        ! Tail of the first decade.
        pos(k+1 : k+nref-istart+1) = 10.0_8**(d0 + ref(istart:nref))
        k = k + nref - istart + 1

        ! Full intermediate decades.
        do d = nint(d0) + 1, nint(dm) - 1
            pos(k+1 : k+nref) = 10.0_8**(real(d, 8) + ref)
            k = k + nref
        end do

        ! Head of the last decade.
        pos(k+1 : k+iend) = 10.0_8**(dm + ref(1:iend))
    end if

    ! Mirror about zero.
    n_total = 2*n_pos + 1
    allocate(q(n_total))
    q = [ -pos(n_pos:1:-1), 0.0_8, pos ]
  END SUBROUTINE symmetric_log_grid

END MODULE handyfuncs_mod
