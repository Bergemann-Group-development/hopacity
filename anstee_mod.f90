!===============================================================================
!> anstee_mod -- hydrogen-collision cross-sections for Turbospectrum line lists.
!>
!> A Turbospectrum line list leaves the van der Waals broadening of an atomic line
!> to Turbospectrum's anstee.f whenever its fdamp column is below 20, i.e. not an
!> already packed ABO cross-section.  anstee.f looks the line up in two places, in
!> this order, and the later one wins:
!>
!>   * the ABO tables of Anstee & O'Mara (1995, s-p), Barklem & O'Mara (1997, p-d)
!>     and Barklem, O'Mara & Ross (1998, d-f), interpolated in the effective
!>     principal quantum numbers of the two levels.  They apply to neutral lines
!>     whose orbital symbols (columns 8-9 of the record) name one of those pairs.
!>   * individually computed values for about 42,000 lines (Barklem and
!>     collaborators, among them all Fe II and Cr II lines of Barklem &
!>     Aspelund-Johansson 2005), matched by wavelength within 0.05 AA and lower
!>     level energy within 0.02 eV.  These apply whatever the orbital symbols.
!>
!> A neutral d-p line that is in neither gets a fit in n* of its upper p level
!> alone.  Only when all of this fails does Turbospectrum fall back on the fdamp
!> column (log gamma_6 when negative, Unsold fudge factor when positive), which is
!> then what HOP does too.  The module reproduces that lookup, so that a
!> Turbospectrum list means the same broadening here as there; the caller packs
!> the result as an ABO value (linelist_mod%build_linelist), which is also what
!> linelist_out.txt then states.
!>
!> The data file is Turbospectrum's DATA/Anstee-1802-2014.dat (GPL-3,
!> github.com/bertrandplez/Turbospectrum_NLTE), read once on first use.
!===============================================================================
MODULE anstee_mod
  USE prec_mod
  USE console_mod
  implicit none
  private

  character(len=4096), public :: anstee_file = './input/Anstee-1802-2014.txt'
  public :: ensure_anstee_loaded, anstee_cross_section

  !> Rydberg energy and the wavelength-energy product, as anstee.f uses them
  real(dp), parameter :: RYD_EV = 13.595_dp, HC_EV_AA = 12399.8_dp
  !> anstee.f converts the cm^-1 energies of the special lines with 8067 cm^-1/eV
  real(dp), parameter :: CM_PER_EV = 8067._dp
  !> record counts of the special-line sections, in file order: general, Fe II, Cr II
  integer, parameter :: NSPEC(3) = [5024, 24188, 13165]
  !> lines skipped before each section
  integer, parameter :: NSKIP(3) = [14, 5, 6]

  !> table pair k (sigma, alpha) indexed (10*n*_upper-orbital, 10*n*_lower-orbital)
  real(dp), save :: tab(10:50, 10:50, 6) = 0

  !> one section of special lines, sorted by wavelength; pos is the file position
  !> within the section, which breaks ties the way anstee.f's first match does
  type special_t
    integer,  allocatable :: iel(:), ion(:), pos(:)
    real(dp), allocatable :: lam(:), elo(:), sigma(:), alpha(:)
  end type
  type(special_t), save :: spec(3)
  logical, save :: loaded = .false.

CONTAINS

!===============================================================================
!> Read the data file.  Idempotent; call from serial setup before the first lookup.
!===============================================================================
SUBROUTINE ensure_anstee_loaded ()
  integer :: lunit, ios, k
  !-----------------------------------------------------------------------------
  if (loaded) return
  open (newunit=lunit, file=trim(anstee_file), status='old', form='formatted', iostat=ios)
  if (ios /= 0) then
    call console_error ('anstee_mod: cannot open '//trim(anstee_file))
    error stop 'anstee_mod: missing ABO data file (anstee_file)'
  end if
  ! the s-p, p-d and d-f pairs; rows are the lower orbital in anstee.f's order
  call read_table (1, 10, 30, 13, 30)
  call read_table (2, 10, 30, 13, 30)
  call read_table (3, 13, 30, 23, 40)
  call read_table (4, 13, 30, 23, 40)
  call read_table (5, 23, 40, 33, 50)
  call read_table (6, 23, 40, 33, 50)
  do k = 1, 3
    call read_special (spec(k), NSKIP(k), NSPEC(k))
  end do
  close (lunit)
  loaded = .true.
CONTAINS
  !-----------------------------------------------------------------------------
  !> One table: comment lines starting with '*', then rows i0..i1 of columns j0..j1.
  SUBROUTINE read_table (k, i0, i1, j0, j1)
    integer, intent(in) :: k, i0, i1, j0, j1
    character(len=256) :: line
    integer :: i
    do
      read (lunit, '(a)', iostat=ios) line
      call check ()
      if (line(1:1) /= '*') exit
    end do
    read (line, *, iostat=ios) tab(j0:j1, i0, k)
    call check ()
    do i = i0+1, i1
      read (lunit, *, iostat=ios) tab(j0:j1, i, k)
      call check ()
    end do
  END SUBROUTINE read_table
  !-----------------------------------------------------------------------------
  SUBROUTINE read_special (s, nskip, n)
    type(special_t), intent(inout) :: s
    integer, intent(in) :: nskip, n
    integer, allocatable :: idx(:)
    integer :: i
    real(dp) :: ecm
    allocate (s%iel(n), s%ion(n), s%pos(n), s%lam(n), s%elo(n), s%sigma(n), s%alpha(n))
    do i = 1, nskip
      read (lunit, *, iostat=ios)
      call check ()
    end do
    do i = 1, n
      read (lunit, '(i2,3x,i2,f10.3,f11.3,59x,f8.2,f8.3)', iostat=ios) &
        s%iel(i), s%ion(i), s%lam(i), ecm, s%sigma(i), s%alpha(i)
      call check ()
      s%elo(i) = ecm / CM_PER_EV
      s%pos(i) = i
    end do
    allocate (idx(n))
    call sort_index (s%lam, idx)
    s%iel = s%iel(idx); s%ion = s%ion(idx); s%pos = s%pos(idx)
    s%lam = s%lam(idx); s%elo = s%elo(idx)
    s%sigma = s%sigma(idx); s%alpha = s%alpha(idx)
  END SUBROUTINE read_special
  !-----------------------------------------------------------------------------
  SUBROUTINE check ()
    if (ios /= 0) then
      call console_error ('anstee_mod: '//trim(anstee_file)//' ends early or is malformed')
      error stop 'anstee_mod: unreadable ABO data file'
    end if
  END SUBROUTINE check
END SUBROUTINE ensure_anstee_loaded

!===============================================================================
!> anstee.f's cross-section sigma (atomic units, at 10^4 m/s) and velocity
!> exponent alpha for one line, and whether it found one (else Turbospectrum
!> would use the line list's fdamp).  iel is the atomic number, ion the stage
!> (1 = neutral), lam_air the wavelength [AA] and elo the lower level energy [eV]
!> as the record states them, eion the stage's ionisation energy [eV], and
!> lev_lo/lev_up the orbital symbols of columns 8-9.
!===============================================================================
SUBROUTINE anstee_cross_section (iel, ion, lam_air, elo, eion, lev_lo, lev_up, &
                                 sigma, alpha, found)
  integer,          intent(in)  :: iel, ion
  real(dp),         intent(in)  :: lam_air, elo, eion
  character(len=1), intent(in)  :: lev_lo, lev_up
  real(dp),         intent(out) :: sigma, alpha
  logical,          intent(out) :: found
  character(len=2) :: pair
  real(dp) :: eup, nlo, nup
  integer :: k
  !-----------------------------------------------------------------------------
  if (.not. loaded) error stop 'anstee_cross_section: call ensure_anstee_loaded first'
  sigma = 0
  alpha = 0
  found = .false.
  eup = HC_EV_AA / lam_air + elo
  if (eup >= eion) return
  if (ion < 1 .or. ion > 2) return

  ! the ABO tables, for neutral lines; the table rows are the lower orbital of
  ! the pair (s of s-p, p of p-d, d of d-f) whichever level it is
  if (ion == 1) then
    pair = lower(lev_lo)//lower(lev_up)
    nlo = nstar (elo)
    nup = nstar (eup)
    select case (pair)
    case ('sp');  call from_table (1, nlo, nup, 1.0_dp, 3.0_dp, 1.3_dp, 3.0_dp)
    case ('ps');  call from_table (1, nup, nlo, 1.0_dp, 3.0_dp, 1.3_dp, 3.0_dp)
    case ('pd');  call from_table (3, nlo, nup, 1.3_dp, 3.0_dp, 2.3_dp, 4.0_dp)
    case ('dp');  call from_table (3, nup, nlo, 1.3_dp, 3.0_dp, 2.3_dp, 4.0_dp)
    case ('df');  call from_table (5, nlo, nup, 2.3_dp, 4.0_dp, 3.3_dp, 5.0_dp)
    case ('fd');  call from_table (5, nup, nlo, 2.3_dp, 4.0_dp, 3.3_dp, 5.0_dp)
    end select
  end if

  ! individually computed lines take precedence; Fe II and Cr II have their own lists
  if (iel == 24 .and. ion == 2) then
    k = 3
  else if (iel == 26 .and. ion == 2) then
    k = 2
  else
    k = 1
  end if
  call from_special (spec(k))

  ! the last resort for neutral d-p lines: the upper p level alone (Barklem 1997)
  if (ion == 1 .and. .not. found .and. lower(lev_lo) == 'd' .and. lower(lev_up) == 'p') then
    if (nstar(elo) < 2.0_dp) then
      nup = nstar (eup)
      found = .true.
      if (nup >= 2.6_dp .and. nup <= 2.7_dp) then
        sigma = 762._dp + 960._dp*(nup-2.6_dp);  alpha = 0.270_dp - 0.04_dp*(nup-2.6_dp)
      else if (nup > 2.7_dp .and. nup <= 2.8_dp) then
        sigma = 858._dp + 1110._dp*(nup-2.7_dp); alpha = 0.266_dp + 0.02_dp*(nup-2.7_dp)
      else if (nup > 2.8_dp .and. nup <= 2.9_dp) then
        sigma = 969._dp + 610._dp*(nup-2.8_dp);  alpha = 0.268_dp - 0.10_dp*(nup-2.8_dp)
      else if (nup > 2.9_dp .and. nup <= 3.0_dp) then
        sigma = 1030._dp + 440._dp*(nup-2.9_dp); alpha = 0.258_dp - 0.04_dp*(nup-2.9_dp)
      else
        found = .false.
      end if
    end if
  end if
CONTAINS
  !-----------------------------------------------------------------------------
  !> Effective principal quantum number of a level at energy e [eV].
  PURE FUNCTION nstar (e) RESULT (n)
    real(dp), intent(in) :: e
    real(dp) :: n
    n = 1._dp / sqrt((eion - e) / RYD_EV)
  END FUNCTION nstar
  !-----------------------------------------------------------------------------
  !> Bilinear interpolation in table pair k at (n1 = row orbital, n2 = column
  !> orbital), when both lie within the table's ranges.
  SUBROUTINE from_table (k, n1, n2, lo1, hi1, lo2, hi2)
    integer,  intent(in) :: k
    real(dp), intent(in) :: n1, n2, lo1, hi1, lo2, hi2
    real(dp) :: f1, f2
    integer :: i1, i2
    if (n1 < lo1 .or. n1 > hi1 .or. n2 < lo2 .or. n2 > hi2) return
    i1 = int(10*n1);  f1 = 10*n1 - i1
    i2 = int(10*n2);  f2 = 10*n2 - i2
    sigma = (tab(i2,i1,k)*(1-f2) + tab(i2+1,i1,k)*f2) * (1-f1) &
          + (tab(i2,i1+1,k)*(1-f2) + tab(i2+1,i1+1,k)*f2) * f1
    alpha = (tab(i2,i1,k+1)*(1-f2) + tab(i2+1,i1,k+1)*f2) * (1-f1) &
          + (tab(i2,i1+1,k+1)*(1-f2) + tab(i2+1,i1+1,k+1)*f2) * f1
    found = .true.
  END SUBROUTINE from_table
  !-----------------------------------------------------------------------------
  !> The special line of section s within 0.05 AA and 0.02 eV, the earliest in
  !> the file when several qualify.
  SUBROUTINE from_special (s)
    type(special_t), intent(in) :: s
    integer :: lo, hi, mid, i, best
    lo = 1
    hi = size(s%lam) + 1
    do while (lo < hi)                       ! first entry with lam >= lam_air - 0.05
      mid = (lo + hi) / 2
      if (s%lam(mid) < lam_air - 0.05_dp) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    best = 0
    do i = lo, size(s%lam)
      if (s%lam(i) > lam_air + 0.05_dp) exit
      if (s%iel(i) /= iel .or. s%ion(i) /= ion) cycle
      if (abs(elo - s%elo(i)) >= 0.02_dp) cycle
      if (best == 0) then
        best = i
      else if (s%pos(i) < s%pos(best)) then
        best = i
      end if
    end do
    if (best > 0) then
      sigma = s%sigma(best)
      alpha = s%alpha(best)
      found = .true.
    end if
  END SUBROUTINE from_special
  !-----------------------------------------------------------------------------
  PURE FUNCTION lower (c) RESULT (l)
    character(len=1), intent(in) :: c
    character(len=1) :: l
    l = c
    if (c >= 'A' .and. c <= 'Z') l = achar(iachar(c) + 32)
  END FUNCTION lower
END SUBROUTINE anstee_cross_section

!===============================================================================
!> Indices that sort a ascending (a stable bottom-up merge sort).
!===============================================================================
SUBROUTINE sort_index (a, idx)
  real(dp), intent(in)  :: a(:)
  integer,  intent(out) :: idx(:)
  integer, allocatable :: tmp(:)
  integer :: n, width, lo, mid, hi, i, j, k
  !-----------------------------------------------------------------------------
  n = size(a)
  idx = [(i, i=1, n)]
  allocate (tmp(n))
  width = 1
  do while (width < n)
    lo = 1
    do while (lo <= n)
      mid = min(lo + width - 1, n)
      hi  = min(lo + 2*width - 1, n)
      i = lo;  j = mid + 1;  k = lo
      do while (i <= mid .and. j <= hi)
        if (a(idx(j)) < a(idx(i))) then
          tmp(k) = idx(j);  j = j + 1
        else
          tmp(k) = idx(i);  i = i + 1
        end if
        k = k + 1
      end do
      do while (i <= mid);  tmp(k) = idx(i);  i = i + 1;  k = k + 1;  end do
      do while (j <= hi);   tmp(k) = idx(j);  j = j + 1;  k = k + 1;  end do
      lo = lo + 2*width
    end do
    idx = tmp
    width = 2*width
  end do
END SUBROUTINE sort_index

END MODULE anstee_mod
