MODULE hvcs_mod
  !  Tabulated VCS (Vidal-Cooper-Smith) hydrogen line profiles of M. Lemke
  !  (1997, A&AS 122, 285), pre-convolved with thermal Doppler broadening,
  !  Ali & Griem resonance broadening and radiative damping. Ported from the
  !  SIU/lineform routines READ_HYDPROF and HYDL_AG_CON
  !  (T. Gehren / F. Grupp / M. Bergemann) and the QPOL quadratic interpolator.
  !
  !  Provides an OPTIONAL alternative to the analytic Griem+VCS Stark profiles
  !  in hprofl_mod, selected at runtime with hydrogen_vcs_table. It covers the
  !  Balmer (n=2), Paschen (n=3) and Brackett (n=4) series for the line members
  !  present in the grids; Lyman and any unlisted member fall back to the
  !  analytic path in hprofl_mod.
  !
  !  The tabulated profile is the WHOLE (already Doppler+resonance-convolved)
  !  line shape, so it replaces the analytic Voigt-core + Stark-wing
  !  construction. phi_vcs returns a pure lineshape normalised to unit area in
  !  frequency, matching the convention of hvoigt/stark1; the caller applies the
  !  populations (nlo*Bij).
  USE const_mod, only: cc
  implicit none
  private

  !  VCS table grid dimensions (from the SIU lf_param.inc); the product
  !  NNHW*NNHT*NNHE == 2288 is asserted against every profile block on read.
  integer, parameter :: NNHW   = 44   ! wavelength (delta-lambda) points
  integer, parameter :: NNHT   = 4    ! temperature nodes
  integer, parameter :: NNHE   = 13   ! electron-density nodes
  integer, parameter :: NNHH   = 13   ! hydrogen-density (resonance) nodes
  integer, parameter :: NNHMAX = 18   ! max lines per series
  integer, parameter :: NSER   = 3    ! Balmer, Paschen, Brackett

  !  Per-series grid + axes (log10 quantities), loaded once at setup.
  real(4), save :: grid_vcs(NNHW,NNHH,NNHT,NNHE,NNHMAX,NSER) = 0.0
  real(4), save :: nw_vcs(NNHW,NSER) = 0.0   ! log10 delta-lambda [AA]
  real(4), save :: nt_vcs(NNHT,NSER) = 0.0   ! log10 T
  real(4), save :: ne_vcs(NNHE,NSER) = 0.0   ! log10 Ne
  real(4), save :: nh_vcs(NNHH,NSER) = 0.0   ! log10 N(HI)
  integer, save :: nlow(NNHMAX,NSER)  = 0
  integer, save :: nup (NNHMAX,NSER)  = 0
  integer, save :: nhhnr(NNHMAX,NSER) = 0
  integer, save :: nlines(NSER)       = 0
  logical, save :: vcs_loaded = .false.

  !  Runtime controls (exposed in the linelist_params namelist).
  logical, public           :: hydrogen_vcs_table = .true.
  character(len=256), public :: vcs_data_dir = './input/vcs'

  real(8), parameter :: AA_PER_S = cc*1.d8   ! speed of light [AA/s]

  !  Thread-private single-cell memo: the (NH,T,Ne) interpolation of the grid is
  !  reused across all frequencies of one cell/line. Safe under OpenMP; degrades
  !  to a plain cache when OpenMP is off.
  real(4),    save :: cache_lp(NNHW) = 0.0, cache_lognorm = 0.0
  integer,    save :: cache_n = -1, cache_m = -1
  real(4),    save :: cache_t = -1.0, cache_e = -1.0, cache_h = -1.0
  logical,    save :: cache_ok = .false.
  !$omp threadprivate(cache_lp, cache_lognorm, cache_n, cache_m, &
  !$omp               cache_t, cache_e, cache_h, cache_ok)

  public :: ensure_vcs_loaded, vcs_supported, phi_vcs

contains

  !====================================================================
  !  Load all three VCS grid files once. Idempotent; call during (serial)
  !  setup, before any profile evaluation.
  SUBROUTINE ensure_vcs_loaded(dir)
    character(len=*), intent(in), optional :: dir
    character(len=256) :: d
    !------------------------------------------------------------------
    if (vcs_loaded) return
    d = vcs_data_dir
    if (present(dir)) d = dir
    call read_vcs_grid(1, trim(d)//'/vcs_balmer.dat')
    call read_vcs_grid(2, trim(d)//'/vcs_paschen.dat')
    call read_vcs_grid(3, trim(d)//'/vcs_bracket.dat')
    vcs_loaded = .true.
  END SUBROUTINE ensure_vcs_loaded

  !====================================================================
  SUBROUTINE read_vcs_grid(s, fname)
    integer, intent(in) :: s
    character(len=*), intent(in) :: fname
    integer :: lunit, ios, k, i, j, kk, l
    real(4) :: dw0,ddw,ne0,dne,nh0,dnh,lgt0,dlgt
    !------------------------------------------------------------------
    open(newunit=lunit, file=fname, status='old', form='formatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'hvcs_mod: cannot open VCS grid file ', trim(fname)
      error stop 'hvcs_mod: missing VCS data file'
    endif

    read(lunit,*) dw0,ddw,ne0,dne,nh0,dnh,lgt0,dlgt

    !  One block per line: nlow nup nhhnr, then grid_vcs(i,j,kk,l). nhhnr is read
    !  before the implied-do that uses it as a bound (as in SIU read_hydprof).
    k = 1
    do
      read(lunit,*,iostat=ios) nlow(k,s), nup(k,s), nhhnr(k,s), &
        ((((grid_vcs(i,j,kk,l,k,s), i=1,NNHW), j=1,nhhnr(k,s)), kk=1,NNHT), l=1,NNHE)
      if (ios /= 0) exit
      if (nhhnr(k,s) > NNHH) then
        write(*,*) 'hvcs_mod: nhhnr exceeds NNHH in ', trim(fname), nhhnr(k,s)
        error stop 'hvcs_mod: VCS grid dimension mismatch'
      endif
      k = k + 1
      if (k > NNHMAX) exit
    enddo
    nlines(s) = k - 1
    close(lunit)

    if (nlines(s) < 1) then
      write(*,*) 'hvcs_mod: no line blocks parsed in ', trim(fname)
      error stop 'hvcs_mod: empty/unreadable VCS grid'
    endif

    !  Axis vectors (all log10).
    do i=1,NNHW ; nw_vcs(i,s)  = dw0  + (i-1)*ddw  ; enddo
    do kk=1,NNHT; nt_vcs(kk,s) = lgt0 + (kk-1)*dlgt; enddo
    do l=1,NNHE ; ne_vcs(l,s)  = ne0  + (l-1)*dne  ; enddo
    do j=1,NNHH ; nh_vcs(j,s)  = nh0  + (j-1)*dnh  ; enddo
  END SUBROUTINE read_vcs_grid

  !====================================================================
  !  True when the tabulated path is enabled, loaded, and the line (n,m) has a
  !  grid entry.
  LOGICAL FUNCTION vcs_supported(n, m)
    integer, intent(in) :: n, m
    !------------------------------------------------------------------
    vcs_supported = hydrogen_vcs_table .and. vcs_loaded .and. (vcs_row(n,m) > 0)
  END FUNCTION vcs_supported

  !====================================================================
  !  Row index of line (n,m) in its series grid, or 0 if absent. Paschen/Brackett
  !  grids do not start at the alpha member, so we match (nlow,nup) rather than
  !  assume row = m-n.
  INTEGER FUNCTION vcs_row(n, m) result(r)
    integer, intent(in) :: n, m
    integer :: s, k
    !------------------------------------------------------------------
    r = 0
    s = n - 1                       ! 2->1 (Balmer), 3->2 (Paschen), 4->3 (Brackett)
    if (s < 1 .or. s > NSER) return
    do k=1, nlines(s)
      if (nlow(k,s) == n .and. nup(k,s) == m) then
        r = k
        return
      endif
    enddo
  END FUNCTION vcs_row

  !====================================================================
  !  Tabulated VCS profile at fractional frequency detuning q for line (n,m).
  !  Returned in the same convention as hvoigt/stark1: a lineshape normalised to
  !  unit area in the FRACTIONAL detuning q = (nu-nu0)/nu0 (i.e. phi_nu*nu0), so
  !  it plugs directly into add_h_opacities/hprofile where the caller applies
  !  nlo*Bij. (A per-Hz normalisation would be a factor nu0 ~1e15 too small.)
  !
  !  The heavy (NH,T,Ne) interpolation onto the delta-lambda grid is memoised per
  !  cell/line, so a whole frequency scan costs one grid interpolation plus a
  !  cheap quadratic wavelength interpolation per point.
  REAL(8) FUNCTION phi_vcs(n, m, q, nu0, temp, ne, nh) result(phi)
    integer, intent(in) :: n, m
    real(8), intent(in) :: q, nu0
    real(4), intent(in) :: temp, ne, nh
    real(8) :: lam0_AA, dlam_AA, deltaw, p, dwmin
    !------------------------------------------------------------------
    if (cache_n /= n .or. cache_m /= m .or. temp /= cache_t .or. &
        ne /= cache_e .or. nh /= cache_h) then
      call vcs_cell_profile(n, m, temp, ne, nh, cache_lp, cache_lognorm, cache_ok)
      cache_n = n ; cache_m = m
      cache_t = temp ; cache_e = ne ; cache_h = nh
    endif

    if (.not. cache_ok) then
      phi = 0.d0
      return
    endif

    !  Line centre wavelength and the requested detuning in Angstrom.
    lam0_AA = AA_PER_S / nu0
    dlam_AA = abs(q) * lam0_AA
    !  delta-lambda grid is identical across series; use Balmer (s=1) axis.
    dwmin  = 10.d0**nw_vcs(1,1)
    deltaw = log10(max(dlam_AA, dwmin))
    !  Quadratic (qpol) interpolation along the wavelength axis, matching SIU
    !  hydl_ag_con. Beyond the outer table edge (~200 AA) qpol extrapolates
    !  linearly in log-log (its goto-50 branch), continuing the asymptotic Stark
    !  power law without the divergence a cubic extrapolation would introduce.
    p = qpol(real(deltaw,4), nw_vcs(1,1), cache_lp, NNHW)

    !  10**(p-lognorm) is phi_lambda [1/AA, unit area in lambda]. The Jacobian
    !  |dlambda/dq| = lam0 converts to phi_q [unit area in q], matching hvoigt.
    phi = 10.d0**(dble(p) - dble(cache_lognorm)) * lam0_AA
  END FUNCTION phi_vcs

  !====================================================================
  !  Interpolate the stored grid over (log NH, log T, log Ne) for one cell,
  !  yielding the log10 profile lp(:) on the delta-lambda grid, plus lognorm =
  !  log10 of its area in delta-lambda (used to renormalise to unit area).
  SUBROUTINE vcs_cell_profile(n, m, temp, ne, nh, lp, lognorm, ok)
    integer, intent(in)  :: n, m
    real(4), intent(in)  :: temp, ne, nh
    real(4), intent(out) :: lp(NNHW), lognorm
    logical, intent(out) :: ok
    integer :: s, r, i, kk, l, nn
    real(4) :: tlog, nelog, nhlog, area, dl
    real(4) :: vt(NNHT), ve(NNHE), vh(NNHH)
    !------------------------------------------------------------------
    ok = .false.
    lp = 0.0 ; lognorm = 0.0
    r = vcs_row(n, m)
    if (r == 0) return
    s  = n - 1
    nn = nhhnr(r,s)

    tlog  = log10(max(temp, 1.0))
    nelog = log10(max(ne, 1.0))
    nhlog = log10(max(nh, 1.0))
    !  Electron densities below the table floor are Doppler dominated (SIU: use
    !  the 10^10 cm^-3 entry). Clamp the other axes to avoid extrapolation.
    if (nelog < ne_vcs(1,s)) nelog = ne_vcs(1,s)
    nelog = min(nelog, ne_vcs(NNHE,s))
    tlog  = min(max(tlog,  nt_vcs(1,s)), nt_vcs(NNHT,s))
    if (nn > 1) nhlog = min(max(nhlog, nh_vcs(1,s)), nh_vcs(nn,s))

    do i=1, NNHW
      do l=1, NNHE
        do kk=1, NNHT
          if (nn > 1) then
            !  Ali & Griem resonance broadening (NH axis, dim 2) present. Copy the
            !  slice to a contiguous buffer: dim 2 is not contiguous from (i,1,...).
            vh(1:nn) = grid_vcs(i,1:nn,kk,l,r,s)
            vt(kk) = qpol4(nhlog, nh_vcs(1,s), vh, nn)
          else
            vt(kk) = grid_vcs(i,1,kk,l,r,s)
          endif
        enddo
        ve(l) = qpol4(tlog, nt_vcs(1,s), vt, NNHT)
      enddo
      lp(i) = qpol4(nelog, ne_vcs(1,s), ve, NNHE)
    enddo

    !  Area under the (symmetric) profile in delta-lambda [AA].
    area = 10.0**lp(1) * (10.0**nw_vcs(1,s))          ! core [0, first node]
    do i=1, NNHW-1
      dl = 10.0**nw_vcs(i+1,s) - 10.0**nw_vcs(i,s)
      area = area + 0.5*(10.0**lp(i) + 10.0**lp(i+1))*dl
    enddo
    area = 2.0*area                                   ! red + blue wings
    lognorm = log10(area)
    ok = .true.
  END SUBROUTINE vcs_cell_profile

  !====================================================================
  !  Quadratic (3-point) interpolation, ported verbatim from SIU qpol.f90 so the
  !  tabulated profiles are reproduced with the same scheme. a(:) must be sorted
  !  ascending; array-element actual arguments (sequence association) supply the
  !  per-series column.
  REAL(4) FUNCTION qpol(x, a, f, no) result(res)
    integer, intent(in) :: no
    real(4), intent(in) :: x, a(no), f(no)
    integer :: i, k, nn
    real(4) :: u, v, w, y, z
    !------------------------------------------------------------------
    k = no - 1
    if (x > a(no)) goto 50
    k = 1
    if (x < a(1) .or. no == 2) goto 50
    nn = no
10  if (k == nn-1) goto 30
    i = (k + nn) / 2
    if (a(i) >= x) goto 20
    k = i
    goto 10
20  nn = i
    goto 10
30  if (k == no-1) k = k - 1
    u = x - a(k)
    v = x - a(k+1)
    w = x - a(k+2)
    y = (f(k)   - f(k+1)) / (v-u)
    z = (f(k+1) - f(k+2)) / (w-v)
    if (y*z > 0.0) goto 40
    res = f(k) + u * (y + v*(y-z)/(w-u))
    return
40  res = f(k) + u*y / (1.0 - v*(y-z)/(f(k)-f(k+2)))
    return
50  res = f(k) + (f(k+1)-f(k)) * (x-a(k)) / (a(k+1)-a(k))
  END FUNCTION qpol

  !====================================================================
  !  Piecewise cubic HERMITE interpolation, C1 continuous. Used for the grid
  !  axes (NH,T,Ne). Each node k carries a slope d(k) = the derivative at a(k) of
  !  the quadratic through the three nearest nodes (one-sided at the ends); that
  !  slope is the SAME whichever of the two bracketing intervals is evaluated, so
  !  value AND first derivative are continuous at every node. Unlike the 3-point
  !  qpol or a 4-point Lagrange (whose stencil shifts at each node -> only C0,
  !  leaving slope kinks that print as wiggly wings), this scheme has no node
  !  kinks. a(:) ascending; inputs on the grid axes are pre-clamped to
  !  [a(1),a(no)] so no extrapolation occurs. Falls back to qpol for <3 nodes.
  REAL(4) FUNCTION qpol4(x, a, f, no) result(res)
    integer, intent(in) :: no
    real(4), intent(in) :: x, a(no), f(no)
    integer :: j, kk
    real(4) :: h, t, t2, t3, dl, dr, h00, h10, h01, h11
    !------------------------------------------------------------------
    if (no < 3) then
      res = qpol(x, a, f, no)
      return
    endif
    !  Bracketing interval [a(j),a(j+1)] (largest j with a(j) <= x, in [1,no-1]).
    j = 1
    do kk = 1, no-1
      if (x >= a(kk)) j = kk
    enddo
    !  Node slopes at the two interval ends (shared with the neighbouring
    !  intervals -> C1), then the cubic Hermite basis on t in [0,1].
    dl = node_slope(a, f, no, j)
    dr = node_slope(a, f, no, j+1)
    h  = a(j+1) - a(j)
    t  = (x - a(j)) / h
    t2 = t*t
    t3 = t2*t
    h00 =  2.0*t3 - 3.0*t2 + 1.0
    h10 =      t3 - 2.0*t2 + t
    h01 = -2.0*t3 + 3.0*t2
    h11 =      t3 -     t2
    res = h00*f(j) + h10*h*dl + h01*f(j+1) + h11*h*dr
  END FUNCTION qpol4

  !====================================================================
  !  Slope at node k: derivative at a(k) of the quadratic through the three
  !  nearest nodes (centred for interior nodes, one-sided at the ends). Single-
  !  valued per node, which is what gives qpol4 its C1 continuity.
  REAL(4) FUNCTION node_slope(a, f, no, k) result(d)
    integer, intent(in) :: no, k
    real(4), intent(in) :: a(no), f(no)
    integer :: i0
    real(4) :: d12, d23, d123, xe
    !------------------------------------------------------------------
    i0   = min(max(k-1, 1), no-2)              ! left index of the 3-node triplet
    d12  = (f(i0+1) - f(i0))   / (a(i0+1) - a(i0))
    d23  = (f(i0+2) - f(i0+1)) / (a(i0+2) - a(i0+1))
    d123 = (d23 - d12) / (a(i0+2) - a(i0))
    xe   = a(k)
    d    = d12 + d123*((xe - a(i0)) + (xe - a(i0+1)))
  END FUNCTION node_slope

END MODULE hvcs_mod
