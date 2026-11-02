MODULE hvcs_mod
  !  Tabulated hydrogen line profiles, selected at runtime with hydrogen_profiles:
  !
  !  'stehle'  Stark profiles of Stehle & Hutcheon (1999, A&AS 140, 93; Model
  !            Microfield Method for electrons and ions), convolved with thermal
  !            Doppler broadening, for Balmer and Paschen lines up to n=20 and 21,
  !            with the Lemke tables for the Brackett series. Converted from
  !            Korg.jl's copy of CDS VI/98A onto the Lemke grid; beyond a line's
  !            last tabulated density the profile is kept fixed in dnu/F0
  !            (Holtsmark scaling). They contain no resonance broadening and no
  !            radiative damping.
  !  'lemke'   VCS (Vidal-Cooper-Smith) profiles of M. Lemke (1997, A&AS 122,
  !            285), convolved with thermal Doppler broadening, radiative damping
  !            and, along an N(HI) axis for H alpha to H delta, Ali & Griem
  !            resonance broadening.
  !  'analytic' no table; the Griem+VCS Stark profiles of hprofl_mod.
  !
  !  The reader and interpolation are ported from the SIU/lineform routines
  !  READ_HYDPROF and HYDL_AG_CON (T. Gehren / F. Grupp / M. Bergemann) and the
  !  QPOL quadratic interpolator. Lyman lines and members missing from a table
  !  use the analytic path in hprofl_mod; hprofl_mod adds the broadening a
  !  table lacks (vcs_resonance, vcs_radiative).
  !
  !  The tabulated profile is the whole Doppler-convolved line shape, so it
  !  replaces the analytic Voigt-core + Stark-wing
  !  construction. vcs_prepare builds the profile for one (line, cell) and
  !  vcs_eval samples it at a fractional frequency detuning q = (nu-nu0)/nu0,
  !  returning a lineshape normalised to unit area in q (phi_nu*nu0), matching
  !  the convention of hvoigt/stark1. The tables are symmetric about a single
  !  line centre; hprofl_mod applies the fine-structure component pattern on
  !  top (vcs_fine), as it does in the analytic path, and the populations
  !  (nlo*Bij).
  USE const_mod, only: cc
  USE interpol_mod, only: pchip_slopes, hermite_value
  USE console_mod
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
  !  Longest axis qpol4 is ever handed (NH, T or Ne) -- sizes its scratch.
  integer, parameter :: NNAXMAX = max(NNHH, max(NNHT, NNHE))

  !  Per-series grid + axes (log10 quantities), loaded once at setup.
  real(4), save :: grid_vcs(NNHW,NNHH,NNHT,NNHE,NNHMAX,NSER) = 0.0
  real(4), save :: nw_vcs(NNHW,NSER) = 0.0   ! log10 delta-lambda [AA]
  real(4), save :: nt_vcs(NNHT,NSER) = 0.0   ! log10 T
  real(4), save :: ne_vcs(NNHE,NSER) = 0.0   ! log10 Ne
  real(4), save :: nh_vcs(NNHH,NSER) = 0.0   ! log10 N(HI)
  !  Node slopes along the T axis (qpol4's) of the lines without an N(HI) axis,
  !  where the T interpolation reads the grid directly.
  real(4), save :: tslope_vcs(NNHW,NNHT,NNHE,NNHMAX,NSER) = 0.0
  !  Derived from nw_vcs at load time, so the frequency scan never pays for them:
  !  the delta-lambda nodes in linear AA, and the reciprocal of the (uniform)
  !  log10 node spacing used to bracket without a search.
  real(4), save :: dw_lin(NNHW,NSER)  = 0.0  ! 10**nw_vcs
  real(4), save :: nw_invdx(NSER)     = 0.0  ! 1/ddw
  real(8), save :: dw_first(NSER)     = 0.d0 ! 10**nw_vcs(1,s), the table's inner edge
  real(8), save :: dw_last(NSER)      = 0.d0 ! 10**nw_vcs(NNHW,s), its outer edge
  integer, save :: nlow(NNHMAX,NSER)  = 0
  integer, save :: nup (NNHMAX,NSER)  = 0
  integer, save :: nhhnr(NNHMAX,NSER) = 0
  integer, save :: nlines(NSER)       = 0
  logical, save :: radiative(NSER)    = .false.  ! the series' table includes radiative damping
  logical, save :: vcs_loaded = .false.

  !  Runtime controls (exposed in the composition_params namelist).
  character(len=16), public :: hydrogen_profiles = 'stehle'
  character(len=4096), public :: vcs_data_dir = './input/vcs'

  real(8), parameter :: AA_PER_S = cc*1.d8   ! speed of light [AA/s]
  real(8), parameter :: LN10 = 2.30258509299404568402d0

  !  One cell's worth of prepared profile: the (NH,T,Ne) interpolation of the
  !  grid, which is constant over a whole frequency scan. Callers that sweep
  !  frequency hold one of these as a LOCAL (vcs_prepare once, vcs_eval per
  !  point) -- keeping it off module scope is what keeps the inner loop clear of
  !  threadprivate/TLS lookups, which on macOS gfortran cost a
  !  pthread_getspecific per reference. The (line, cell) stamp makes prepare
  !  idempotent, so a caller that keeps one per cell and line (rt_patch's hvcs)
  !  builds each only once.
  type, public :: vcs_profile_t
    real(4) :: lp(NNHW) = 0.0     ! log10 profile on the delta-lambda grid
    real(8) :: lognorm  = 0.d0    ! log10 of its area in delta-lambda
    real(8) :: lam0     = 0.d0    ! line centre wavelength [AA]
    real(8) :: q_tail   = 0.d0    ! |q| of the outer table edge; beyond it a power law
    real(8) :: ta = 0.d0, ts = 0.d0   ! ... phi = exp(ta + ts*ln|q|), see vcs_tail
    integer :: n = -1, m = -1     ! stamp: line ...
    real(8) :: nu0 = -1.d0        ! ...
    real(4) :: t = -1.0, e = -1.0, h = -1.0   ! ... and cell it was built for
    logical :: ok = .false.
  end type

  public :: ensure_vcs_loaded, vcs_supported, vcs_prepare, vcs_eval, vcs_tail, &
            vcs_resonance, vcs_radiative

contains

  !====================================================================
  !  Load the grid files hydrogen_profiles selects, once. Idempotent; call during
  !  (serial) setup, before any profile evaluation.
  SUBROUTINE ensure_vcs_loaded(dir)
    character(len=*), intent(in), optional :: dir
    character(len=len(vcs_data_dir)) :: d
    !------------------------------------------------------------------
    if (vcs_loaded) return
    d = vcs_data_dir
    if (present(dir)) d = dir
    select case (trim(hydrogen_profiles))
    case ('analytic')
      return
    case ('stehle')
      call read_vcs_grid(1, trim(d)//'/stehle_balmer.dat')
      call read_vcs_grid(2, trim(d)//'/stehle_paschen.dat')
      call read_vcs_grid(3, trim(d)//'/vcs_bracket.dat')
      radiative = [.false., .false., .true.]
    case ('lemke')
      call read_vcs_grid(1, trim(d)//'/vcs_balmer.dat')
      call read_vcs_grid(2, trim(d)//'/vcs_paschen.dat')
      call read_vcs_grid(3, trim(d)//'/vcs_bracket.dat')
      radiative = .true.
    case default
      call console_error ("hvcs_mod: hydrogen_profiles must be 'stehle', 'lemke' or 'analytic', not '" &
                          //trim(hydrogen_profiles)//"'")
      error stop 'hvcs_mod: unknown hydrogen_profiles'
    end select
    vcs_loaded = .true.
  END SUBROUTINE ensure_vcs_loaded

  !====================================================================
  SUBROUTINE read_vcs_grid(s, fname)
    integer, intent(in) :: s
    character(len=*), intent(in) :: fname
    integer :: lunit, ios, k, i, j, kk, l
    real(4) :: dw0,ddw,ne0,dne,nh0,dnh,lgt0,dlgt
    real(4) :: dx(NNHT), dy(NNHT)
    !------------------------------------------------------------------
    open(newunit=lunit, file=fname, status='old', form='formatted', iostat=ios)
    if (ios /= 0) then
      call console_error ('hvcs_mod: cannot open VCS grid file '//trim(fname))
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
        call console_error ('hvcs_mod: nhhnr exceeds NNHH in '//trim(fname))
        error stop 'hvcs_mod: VCS grid dimension mismatch'
      endif
      k = k + 1
      if (k > NNHMAX) exit
    enddo
    nlines(s) = k - 1
    close(lunit)

    if (nlines(s) < 1) then
      call console_error ('hvcs_mod: no line blocks parsed in '//trim(fname))
      error stop 'hvcs_mod: empty/unreadable VCS grid'
    endif

    !  Axis vectors (all log10).
    do i=1,NNHW ; nw_vcs(i,s)  = dw0  + (i-1)*ddw  ; enddo
    do kk=1,NNHT; nt_vcs(kk,s) = lgt0 + (kk-1)*dlgt; enddo
    do l=1,NNHE ; ne_vcs(l,s)  = ne0  + (l-1)*dne  ; enddo
    do j=1,NNHH ; nh_vcs(j,s)  = nh0  + (j-1)*dnh  ; enddo

    !  Delta-lambda quantities the hot paths would otherwise recompute: the nodes
    !  in linear AA (area integral), the table's inner edge (profile clamp) and
    !  the reciprocal node spacing (bracketing, see qpol_w).
    do i=1,NNHW ; dw_lin(i,s) = 10.0**nw_vcs(i,s) ; enddo
    dw_first(s) = 10.d0**nw_vcs(1,s)
    dw_last(s)  = 10.d0**nw_vcs(NNHW,s)
    nw_invdx(s) = 1.0/ddw

    do k=1, nlines(s)
      if (nhhnr(k,s) > 1) cycle
      do l=1,NNHE
        do i=1,NNHW
          call pchip_slopes(nt_vcs(:,s), grid_vcs(i,1,:,l,k,s), 1, NNHT, &
                            dx(1:NNHT-1), dy(1:NNHT-1), tslope_vcs(i,:,l,k,s))
        enddo
      enddo
    enddo
  END SUBROUTINE read_vcs_grid

  !====================================================================
  !  True when the tabulated path is enabled, loaded, and the line (n,m) has a
  !  grid entry.
  LOGICAL FUNCTION vcs_supported(n, m)
    integer, intent(in) :: n, m
    !------------------------------------------------------------------
    vcs_supported = vcs_loaded .and. (vcs_row(n,m) > 0)
  END FUNCTION vcs_supported

  !====================================================================
  !  True when the table of a supported line (n,m) carries resonance broadening
  !  along its N(HI) axis.
  LOGICAL FUNCTION vcs_resonance(n, m)
    integer, intent(in) :: n, m
    !------------------------------------------------------------------
    vcs_resonance = nhhnr(vcs_row(n,m), n-1) > 1
  END FUNCTION vcs_resonance

  !====================================================================
  !  True when the table of series n includes radiative damping.
  LOGICAL FUNCTION vcs_radiative(n)
    integer, intent(in) :: n
    !------------------------------------------------------------------
    vcs_radiative = radiative(n-1)
  END FUNCTION vcs_radiative

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
  !  Build (or reuse) the prepared profile for line (n,m) in one cell. Idempotent
  !  for a repeated (line, cell) stamp, so callers may call it unconditionally.
  SUBROUTINE vcs_prepare(n, m, nu0, temp, ne, nh, c)
    integer, intent(in) :: n, m
    real(8), intent(in) :: nu0
    real(4), intent(in) :: temp, ne, nh
    type(vcs_profile_t), intent(inout) :: c
    real(4) :: lognorm
    !------------------------------------------------------------------
    if (c%n == n .and. c%m == m .and. c%nu0 == nu0 .and. c%t == temp .and. &
        c%e == ne .and. c%h == nh) return
    call vcs_cell_profile(n, m, temp, ne, nh, c%lp, lognorm, c%ok)
    c%lognorm = dble(lognorm)
    c%lam0 = AA_PER_S / nu0          ! line centre wavelength [AA]
    c%q_tail = dw_last(1) / c%lam0   ! vcs_eval uses the Balmer axis for all series
    if (c%ok) call tail_coef(c)
    c%n = n ; c%m = m ; c%nu0 = nu0
    c%t = temp ; c%e = ne ; c%h = nh
  END SUBROUTINE vcs_prepare

  !====================================================================
  !  Evaluate a prepared profile at fractional frequency detuning q. This is the
  !  frequency-scan inner loop: everything that does not depend on q lives in
  !  vcs_prepare, leaving one log10 and one exp per point.
  PURE REAL(8) FUNCTION vcs_eval(c, q) result(phi)
    type(vcs_profile_t), intent(in) :: c
    real(8), intent(in) :: q
    real(8) :: dlam_AA, deltaw, p
    !------------------------------------------------------------------
    if (.not. c%ok) then
      phi = 0.d0
      return
    endif

    !  Requested detuning in Angstrom, clamped to the table's inner edge. The
    !  delta-lambda grid is identical across series; use the Balmer (s=1) axis.
    dlam_AA = abs(q) * c%lam0
    deltaw = log10(max(dlam_AA, dw_first(1)))
    !  Quadratic (qpol) interpolation along the wavelength axis, matching SIU
    !  hydl_ag_con. Beyond the outer table edge (~200 AA) qpol extrapolates
    !  linearly in log-log (its goto-50 branch), continuing the asymptotic Stark
    !  power law without the divergence a cubic extrapolation would introduce.
    p = dble(qpol_w(real(deltaw,4), 1, c%lp))

    !  10**(p-lognorm) is phi_lambda [1/AA, unit area in lambda], written as an
    !  exp because pow(10,x) is ~45% dearer for the same result. The Jacobian
    !  |dlambda/dq| = lam0 converts to phi_q [unit area in q], matching hvoigt.
    phi = exp(LN10*(p - c%lognorm)) * c%lam0
  END FUNCTION vcs_eval

  !====================================================================
  !  The power law that vcs_eval continues past the outer table edge (qpol's
  !  linear log-log extrapolation through the last two nodes), with the area
  !  normalisation and the Jacobian lam0 folded in: ln(phi_q) = ta + ts*ln|q|.
  SUBROUTINE tail_coef(c)
    type(vcs_profile_t), intent(inout) :: c
    real(8) :: slope, a
    !------------------------------------------------------------------
    a = nw_vcs(NNHW-1,1)
    slope = dble(c%lp(NNHW) - c%lp(NNHW-1)) / dble(nw_vcs(NNHW,1) - nw_vcs(NNHW-1,1))
    c%ts = slope
    c%ta = LN10*(dble(c%lp(NNHW-1)) - slope*a - c%lognorm) + (1.d0 + slope)*log(c%lam0)
  END SUBROUTINE tail_coef

  !====================================================================
  !  vcs_eval beyond the outer table edge (|q| > q_tail), from ln|q|. Callers
  !  that sweep many cells at the same frequencies can share ln|q| between them.
  PURE REAL(8) FUNCTION vcs_tail(c, lnq) result(phi)
    type(vcs_profile_t), intent(in) :: c
    real(8), intent(in) :: lnq
    !------------------------------------------------------------------
    if (c%ok) then
      phi = exp(c%ta + c%ts*lnq)
    else
      phi = 0.d0
    endif
  END FUNCTION vcs_tail

  !====================================================================
  !  Interpolate the stored grid over (log NH, log T, log Ne) for one cell,
  !  yielding the log10 profile lp(:) on the delta-lambda grid, plus lognorm =
  !  log10 of its area in delta-lambda (used to renormalise to unit area).
  SUBROUTINE vcs_cell_profile(n, m, temp, ne, nh, lp, lognorm, ok)
    integer, intent(in)  :: n, m
    real(4), intent(in)  :: temp, ne, nh
    real(4), intent(out) :: lp(NNHW), lognorm
    logical, intent(out) :: ok
    integer :: s, r, i, kk, l, nn, jt, je, le0, le1, ns
    real(4) :: tlog, nelog, nhlog, area, dl, plin, plin_next, ht, tt, he, te
    real(4) :: vt(NNHT), ve(NNHE), vh(NNHH), dx(4), dy(4), dydx(4)
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

    !  qpol4's bracketing and Hermite parameters on the T and Ne axes, the same
    !  for every wavelength. The Ne interpolation reads the node slopes at je and
    !  je+1, which depend on the nodes le0..le1 only, so only those columns are
    !  interpolated in T. Without an N(HI) axis that interpolation uses the T
    !  slopes stored at load time.
    jt = bracket(tlog, nt_vcs(:,s))
    ht = nt_vcs(jt+1,s) - nt_vcs(jt,s)
    tt = (tlog - nt_vcs(jt,s)) / ht
    je = bracket(nelog, ne_vcs(:,s))
    he = ne_vcs(je+1,s) - ne_vcs(je,s)
    te = (nelog - ne_vcs(je,s)) / he
    le0 = max(1, je-1)
    le1 = min(NNHE, je+2)
    ns = le1 - le0 + 1
    do i=1, NNHW
      do l=le0, le1
        if (nn > 1) then
          do kk=1, NNHT
            !  Ali & Griem resonance broadening (NH axis, dim 2) present. Copy the
            !  slice to a contiguous buffer: dim 2 is not contiguous from (i,1,...).
            vh(1:nn) = grid_vcs(i,1:nn,kk,l,r,s)
            vt(kk) = qpol4(nhlog, nh_vcs(1,s), vh, nn)
          enddo
          ve(l) = qpol4(tlog, nt_vcs(1,s), vt, NNHT)
        else
          ve(l) = hermite_value(tt, ht, grid_vcs(i,1,jt,l,r,s), grid_vcs(i,1,jt+1,l,r,s), &
                                tslope_vcs(i,jt,l,r,s), tslope_vcs(i,jt+1,l,r,s))
        endif
      enddo
      call pchip_slopes(ne_vcs(le0:le1,s), ve(le0:le1), je-le0+1, je-le0+2, &
                        dx(1:ns-1), dy(1:ns-1), dydx(1:ns))
      lp(i) = hermite_value(te, he, ve(je), ve(je+1), dydx(je-le0+1), dydx(je-le0+2))
    enddo

    !  Area under the (symmetric) profile in delta-lambda [AA]. The nodes in
    !  linear AA come from dw_lin (built at load), and each 10**lp is carried
    !  over to the next trapezium rather than raised twice.
    plin = 10.0**lp(1)
    area = plin * dw_lin(1,s)                         ! core [0, first node]
    do i=1, NNHW-1
      dl = dw_lin(i+1,s) - dw_lin(i,s)
      plin_next = 10.0**lp(i+1)
      area = area + 0.5*(plin + plin_next)*dl
      plin = plin_next
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
  !  qpol specialised to the delta-lambda axis of series s, which read_vcs_grid
  !  lays down as a UNIFORM log10 ladder (nw_vcs(i,s) = dw0 + (i-1)*ddw). The
  !  bracketing index is therefore one multiply instead of qpol's binary search
  !  -- worth having, as this runs once per frequency point. The arithmetic that
  !  follows is qpol verbatim, and the two bracketings agree: qpol's search
  !  returns the largest k with a(k) < x, which differs from the truncation only
  !  when x lands exactly on a node, and there both stencils collapse to f(node).
  PURE REAL(4) FUNCTION qpol_w(x, s, f) result(res)
    real(4), intent(in) :: x
    integer, intent(in) :: s
    real(4), intent(in) :: f(NNHW)
    integer :: k
    real(4) :: u, v, w, y, z
    !------------------------------------------------------------------
    if (x > nw_vcs(NNHW,s)) then
      !  Linear (log-log) extrapolation past the outer edge, qpol's goto-50.
      k = NNHW-1
      res = f(k) + (f(k+1)-f(k)) * (x-nw_vcs(k,s)) / (nw_vcs(k+1,s)-nw_vcs(k,s))
      return
    else if (x < nw_vcs(1,s)) then
      k = 1
      res = f(k) + (f(k+1)-f(k)) * (x-nw_vcs(k,s)) / (nw_vcs(k+1,s)-nw_vcs(k,s))
      return
    endif

    k = int((x - nw_vcs(1,s)) * nw_invdx(s)) + 1
    k = min(max(k, 1), NNHW-2)                        ! qpol's "if (k==no-1) k=k-1"
    u = x - nw_vcs(k,s)
    v = x - nw_vcs(k+1,s)
    w = x - nw_vcs(k+2,s)
    y = (f(k)   - f(k+1)) / (v-u)
    z = (f(k+1) - f(k+2)) / (w-v)
    if (y*z > 0.0) then
      res = f(k) + u*y / (1.0 - v*(y-z)/(f(k)-f(k+2)))
    else
      res = f(k) + u * (y + v*(y-z)/(w-u))
    endif
  END FUNCTION qpol_w

  !====================================================================
  !  qpol4's bracketing interval [a(j),a(j+1)]: the largest j in [1,size(a)-1]
  !  with a(j) <= x, or 1.
  PURE INTEGER FUNCTION bracket(x, a) result(j)
    real(4), intent(in) :: x, a(:)
    integer :: kk
    !------------------------------------------------------------------
    j = 1
    do kk = 1, size(a)-1
      if (x >= a(kk)) j = kk
    enddo
  END FUNCTION bracket

  !====================================================================
  !  Piecewise cubic HERMITE interpolation, C1 continuous. Used for the grid
  !  axes (NH,T,Ne). Node slopes are the canonical monotone (Fritsch-Carlson
  !  harmonic-mean) rule shared with the rest of the code via interpol_mod
  !  (pchip_slopes), of which only the two bracketing ones are actually computed;
  !  the slope at each node is single-valued, so value AND first derivative are
  !  continuous at every node (no wiggly wings from stencil
  !  shifts). a(:) ascending; inputs on the grid axes are pre-clamped to
  !  [a(1),a(no)] so no extrapolation occurs. Falls back to qpol for <3 nodes.
  REAL(4) FUNCTION qpol4(x, a, f, no) result(res)
    integer, intent(in) :: no
    real(4), intent(in) :: x, a(no), f(no)
    integer :: j, kk
    !  Fixed-size scratch: automatic arrays here would be a malloc/free pair per
    !  call, and vcs_cell_profile makes thousands of calls per cell.
    real(4) :: h, t, dx(NNAXMAX), dy(NNAXMAX), dydx(NNAXMAX)
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
    !  Canonical monotone node slopes, then the cubic Hermite basis on t in [0,1]
    !  using the two bracketing node slopes. Only those two are asked for, via
    !  pchip_slopes' imin/imax window: a node slope is a local quantity, so the
    !  window changes what is computed but not the values, and this is the
    !  innermost call of the grid interpolation below.
    call pchip_slopes(a, f, j, j+1, dx(1:no-1), dy(1:no-1), dydx(1:no))
    h  = a(j+1) - a(j)
    t  = (x - a(j)) / h
    res = hermite_value(t, h, f(j), f(j+1), dydx(j), dydx(j+1))
  END FUNCTION qpol4

END MODULE hvcs_mod
