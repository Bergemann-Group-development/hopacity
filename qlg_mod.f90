!=======================================================================================
! qlg_mod.f90 -- the four deployed members of the QLG family of Voigt approximations.
!
! Companion code to "A branchless single-precision Voigt approximation for spectral
! synthesis" (Hoppe).  Every kernel here is the compiled form of a member of Eq. (5) of
! that paper: the two analytic limits, plus a correction whose coefficients are
! polynomials in x = 2r - 1 with r = a/Gamma the damping ratio.  The constants below are
! the ones in the paper's appendix table of coefficients, in that table's convention.
!
! The name is the family's three ingredients: the Lorentzian generator L = 1/(1+u), the
! Gaussian generator G = 2^-u, and the envelope q = r(1-r) that multiplies every fitted
! correction.  q vanishes at r = 0 and r = 1, so the Gaussian and Lorentzian limits are
! exact by construction and no fitted constant can perturb them.
!
!   qlg_toy_sweep   3 constants, 2.2e-2   -- the smallest useful member, Eq. (6)
!   qlg_bot_sweep   6 constants, 5.7e-3   -- Synspec's cost, a much smaller flux error
!   qlg_mid_sweep  14 constants, 9.2e-4   -- flux error below 1e-3 on all three stars
!   qlg_top_sweep  56 constants, 3.2e-5   -- a replacement for Humlicek w4
!
! A fifth public routine, qlg_ref_f64, is the fp64 Humlicek-Weideman reference the four
! are validated against.  It is reference grade and correspondingly slow: it is here to
! check a member, and to build QLG-bot's table, not to be called in a synthesis.
!
! CALLING.  All four share one signature and evaluate a whole line window into an array:
!
!     call qlg_mid_sweep (a, v0, dv, nw, K)      ! K(f) = F(a, v0 + (f-1)*dv), f = 1..nw
!
! with a, v0, dv and K all real(4).  K is the profile in the normalisation K(0,0) = 1 for
! a = 0, i.e. K = Re w(v + i a) with w the Faddeeva function, so a caller multiplies by
! its own 1/(sqrt(pi) Delta_nu_D) as it would for any Voigt routine.  The window form is
! not a convenience: the per-line setup -- the half-width, r, and every coefficient
! polynomial -- is evaluated once per call and amortised over the window, which is where
! most of these members' advantage lies.  A one-point-per-call wrapper throws that away.
!
! INITIALISATION.  qlg_bot_sweep reads two tables and requires
!
!     call qlg_init
!
! once before its first use.  The call is idempotent and costs a few milliseconds, so a
! caller that cannot easily tell which member will be selected may simply always make it.
! The other three members are arithmetic: they need no tables and no initialisation.
!
! fp32.  All four are single precision throughout, and are fitted to be so: the
! coefficient polynomials are evaluated in x = 2r - 1 rather than in r, which keeps their
! constants O(1) and the single-precision loss below 1e-5 (Appendix C of the paper).
!
! SELF-CONTAINED, DELIBERATELY.  This file compiles alone, with no dependency on anything
! else.  Two consequences worth knowing before tidying it:
!
!   - fexp2, the 2^-u kernel, lives in this module rather than in a utility module of its
!     own because a cross-module call inside the point loop does not inline without -flto,
!     and one measured 4.07 ns/point against 1.64 for the same arithmetic.
!   - qlg_init builds its table from the Weideman-32 reference qlg_ref_f64 at the bottom
!     of this file, so the module carries its own reference rather than depending on one.
!
! COEFFICIENT PROVENANCE, for the record: QLG-toy and QLG-mid and QLG-top are the fits
! stored as `Q5r_H0.2_d111_w1` (fit/botshapes.json), `G0_Gm1_L2_G2_G3_L4_G4|d1|whiting0`
! (fit/midarith_scored.json) and `N13|alloc43` (fit/toparith_final.json) in the paper's
! working repository, transcribed here without change.
!
!   gfortran -O3 -mcpu=native -c qlg_mod.f90
!=======================================================================================
MODULE qlg
  implicit none
  private

  public :: qlg_init
  public :: qlg_toy_sweep, qlg_bot_sweep, qlg_mid_sweep, qlg_top_sweep
  public :: qlg_ref_f64
  public :: qlg_bot_phi, qlg_bot_profile          ! local addition, see the bottom of this file
  public :: qlg_bot_coef, qlg_bot_eval             ! local addition, as is qlg_bot_coef_t

  ! QLG-bot's per-line setup for one damping a (local addition, see qlg_bot_coef)
  type, public :: qlg_bot_coef_t
    real(4) :: ih, bL, gG, cZH, cL2
  end type qlg_bot_coef_t

  !-------------------------------------------------------------------------------------
  ! Scalar constants of Eq. (5): ln2, sqrt(pi), 1/pi, sqrt(ln2/pi), sqrt(ln2).
  !-------------------------------------------------------------------------------------
  real(4), parameter :: LN2 = 0.6931471806, RPI = 1.7724538509
  real(4), parameter :: OVERPI = 0.3183098862, GNORM = 0.4697186
  real(4), parameter :: SQL2 = 0.8325546112
  real(8), parameter :: DLN2 = 0.69314718055994531d0
  real(8), parameter :: ISQPI8 = 0.564189583547756287d0   ! 1/sqrt(pi)
  real(4), parameter :: ISQPI  = 0.5641895835              ! 1/sqrt(pi), fp32 (local addition)

  !-------------------------------------------------------------------------------------
  ! minimax p2(f) ~ 2^-f on f in [0,1); G = 2^-u as p2(u-n) * 2^-n. See fexp2.
  !-------------------------------------------------------------------------------------
  real(4), parameter :: EXP2_P(6) = [ &
    +9.999999404e-01, -6.931422353e-01, +2.401721776e-01, &
    -5.527971312e-02, +9.188604541e-03, -9.387863684e-04 ]
  real(4), parameter :: U_CAP = 127.0
  integer, parameter :: F32_BIAS = 127, F32_MANT = 23

  !-------------------------------------------------------------------------------------
  ! QLG-bot's two tables (Appendix C.3).  Both are indexed on xi = |v|/Gamma = sqrt(u) at
  ! XSCL = 800 intervals per unit xi, over 0 <= xi <= NX/XSCL = 11.27, and both are read
  ! at ONE snapped nearest-node index, so the second table costs one load and one FMA.
  !
  !   TGX(i)  = 2^-u_i                      the Gaussian generator
  !   TH2X(i) = Z_0.2(u_i) (1 + u_i)^2      the anchored slice, divided by its L^2 wing
  !
  ! with u_i = (i/XSCL)^2.  The upper limit is xi^2 = 127, past which 2^-u has underflowed
  ! in fp32, so clamping the index gives the correct G = 0 in the far wing and no min() on
  ! u is needed in the point loop.
  !
  ! TH2X is stored DIVIDED by L^2 and the factor restored arithmetically after the read.
  ! Stored bare it would be a far-wing failure: past the index clamp the table freezes at
  ! its last entry, which for the bare slice is an O(1) constant added to a decaying wing,
  ! while the divided form's frozen value is multiplied by an L^2 that keeps decaying.
  !
  ! TWO SEPARATE RANK-1 ARRAYS.  Do not interleave them into one rank-2 table to save an
  ! index computation: that was tried and costs 17 %, because it halves the vector width.
  !-------------------------------------------------------------------------------------
  integer, parameter :: NX   = 9019
  real(4), parameter :: XSCL = 800.0
  real(4), parameter :: RNX  = real(NX)
  real(4), save      :: TGX(0:NX)
  real(4), save      :: TH2X(0:NX)
  logical, save      :: qlg_done = .false.

  !-------------------------------------------------------------------------------------
  ! Weideman (1994) w = 32 rational coefficients, fp64.  Used ONLY by qlg_init, to build
  ! TH2X; no member evaluates them.
  !-------------------------------------------------------------------------------------
  real(8), parameter :: W32_LD = 4.75682846001088411d+00
  real(8), parameter :: W32_AD(32) = [ &
    -1.30317978630500875d-12, +3.74107755472269789d-12, &
    +8.03033538228436012d-12, -2.15435373143760700d-11, &
    -5.54423561143486540d-11, +1.16582546820233235d-10, &
    +4.15374271614691897d-10, -5.23102143724924582d-10, &
    -3.20801514078875819d-09, +8.12488896568210460d-10, &
    +2.37975567593842971d-08, +2.29304391194392872d-08, &
    -1.48130789124396946d-07, -4.18407637053084035d-07, &
    +4.25583313792292626d-07, +4.40153173159092505d-06, &
    +6.82103194402497841d-06, -2.14096192018698432d-05, &
    -1.30754492546239731d-04, -2.45329802700202065d-04, &
    +3.92591360700705167d-04, +4.51954110534923473d-03, &
    +1.90061557848454077d-02, +5.73044035298371918d-02, &
    +1.40607162268937713d-01, +2.95444510715087261d-01, &
    +5.46013972063934205d-01, +9.01925489364799993d-01, &
    +1.34554416923454490d+00, +1.82566962963248147d+00, &
    +2.26353729990026764d+00, +2.57225340812456960d+00 ]

contains

  !=====================================================================================
  ! G = 2^-u, branchless: a degree-5 minimax polynomial on the fractional part times the
  ! integer part spliced straight into the exponent field.  Used by QLG-toy, QLG-mid and
  ! QLG-top.  The clamp at U_CAP = 127 is what makes the splice safe; past it 2^-u has
  ! underflowed in fp32 anyway.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE fexp2(u, E)
    real(4), intent(in) :: u
    real(4), intent(out) :: E
    real(4) :: uc, fr, p
    integer :: j, ne
    !-----------------------------------------------------------------
    uc = min(u, U_CAP)
    ne = int(uc)
    fr = uc - real(ne)
    p = EXP2_P(6)
    do j = 5, 1, -1
      p = p*fr + EXP2_P(j)
    end do
    E = p * transfer(ishft(F32_BIAS - ne, F32_MANT), 1.0)
  END SUBROUTINE fexp2

  !=====================================================================================
  ! QLG-toy -- 3 basis functions at x-degree 0, so the three coefficients are plain
  ! numbers.  Eq. (6) of the paper, on Whiting's half-width.  The smallest member for
  ! which both counts are minimal at its accuracy; stated in the paper to be read, and
  ! kept here so that it can be run.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE qlg_toy_sweep (a, v0, dv, n, K)
    real(4), intent(in)  :: a, v0, dv
    integer, intent(in)  :: n
    real(4), intent(out) :: K(n)
    integer :: f
    real(4) :: hw, ih, ih2, r, q, sc, qs, bL, gG
    real(4) :: c0_2, c2_2
    real(4) :: v, u, L, G, L2, G2
    !-----------------------------------------------------------------
    hw = 0.5346*a + sqrt(0.2166*a*a + LN2)
    ih = 1.0 / hw
    r  = a * ih
    q  = r * (1.0 - r)
    sc = RPI * ih
    qs = q * sc
    ih2 = ih * ih
    bL  = (OVERPI*r) * sc
    gG  = (GNORM*(1.0 - r)) * sc
    c0_2    = qs*(1.413611936)
    gG  = gG + qs*(-0.6357795021)
    c2_2    = qs*(-0.8508010787)

    do f = 1, n
      v = v0 + real(f-1, 4)*dv
      u = v*v * ih2
      L = 1.0 / (1.0 + u)
      call fexp2(u, G)
      L2 = L*L
      G2 = G*G
      K(f) = bL*L + gG*G + L2*(c0_2 + c2_2*G2)
    end do
  END SUBROUTINE qlg_toy_sweep

  !=====================================================================================
  ! QLG-bot -- 3 basis functions {q G, q Z_0.2, q L^2} at x-degree 1, 6 constants, on the
  ! sqrt-free rational half-width of Eq. (4).  The tabulated member: G and Z_0.2/L^2 are
  ! read from TGX and TH2X at one shared snapped index, so call qlg_init first.
  ! Maximum relative error over the validation domain, compiled in fp32: 5.6923e-3.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE qlg_bot_sweep (a, v0, dv, n, K)
    real(4), intent(in)  :: a, v0, dv
    integer, intent(in)  :: n
    real(4), intent(out) :: K(n)
    integer :: f, ig
    real(4) :: d, ih, r, q, sc, qs, x, bL, gG
    real(4) :: cZH, cL2
    real(4) :: v, xs, u, L, G, W, p
    !-----------------------------------------------------------------
    d  = 1.0 + a*0.60602209212341318
    ih = d / (a*d + SQL2)
    r  = a * ih
    q  = r * (1.0 - r)
    sc = RPI * ih
    qs = q * sc
    x  = 2.0*r - 1.0
    bL = (OVERPI*r) * sc
    gG = (GNORM*(1.0 - r)) * sc + qs*(-0.064924393191013452 + x*(-0.11019326278794868))
    cZH = qs*(0.38494195930078284 + x*(-0.10567135279610407))
    cL2 = qs*(0.04054415746208978 + x*(0.16241029253712097))

    do f = 1, n
      v  = v0 + real(f-1, 4)*dv
      xs = abs(v) * ih
      u  = xs * xs
      L  = 1.0 / (1.0 + u)

      ig = nint(min(xs*XSCL, RNX))       ! one snapped index, both tables
      G  = TGX(ig)
      W  = TH2X(ig)

      p = (cL2 + cZH*W)
      p = bL                           + L*p
      K(f) = G*(gG) + L*p
    end do
  END SUBROUTINE qlg_bot_sweep

  !=====================================================================================
  ! QLG-mid -- 7 basis functions at x-degree 1, 14 constants, on Whiting's half-width.
  ! Arithmetic: G comes from fexp2, so there is no table and no gather.  Exactly one
  ! multiply by G stands between fexp2 and the result, which is what keeps a G-coupled
  ! correction as cheap as a pure-L one of the same L-power count.
  ! Maximum relative error, compiled in fp32: 9.2002e-4.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE qlg_mid_sweep (a, v0, dv, n, K)
    real(4), intent(in)  :: a, v0, dv
    integer, intent(in)  :: n
    real(4), intent(out) :: K(n)
    integer :: f
    real(4) :: ih, r, q, sc, qs, x, bL, gG, cGm1, cL2, cG2, cG3, cL4, cG4
    real(4) :: v, u, iL, L, G, t
    !-----------------------------------------------------------------
    ih = 1.0 / (0.5346*a + sqrt(0.2166*a*a + LN2))
    r  = a * ih
    q  = r * (1.0 - r)
    sc = RPI * ih
    qs = q * sc
    x  = 2.0*r - 1.0
    bL = (OVERPI*r) * sc
    gG   = (GNORM*(1.0 - r)) * sc + qs*(-0.7895388603 + x*(0.3148033917))
    cGm1 = qs*(-0.0426536724 + x*(0.09756947309))
    cL2  = qs*(1.187075973 + x*(0.114530772))
    cG2  = qs*(-5.72019577 + x*(11.42426872))
    cG3  = qs*(-10.74044037 + x*(15.57648754))
    cL4  = qs*(32.14587021 + x*(-58.42788315))
    cG4  = qs*(-16.07919693 + x*(30.8985672))

    do f = 1, n
      v  = v0 + real(f-1, 4)*dv
      u  = v*v*ih*ih
      iL = 1.0 + u
      L  = 1.0 / iL
      call fexp2(u, G)   ! kern.f90::fexp2, no table, no gather

      t = cL4 + cG4*G
      t = cG3*G + L*t
      t = (cL2 + cG2*G) + L*t
      t = bL + L*t
      K(f) = G*(gG + cGm1*iL) + L*t
    end do
  END SUBROUTINE qlg_mid_sweep

  !=====================================================================================
  ! QLG-top -- 13 basis functions at x-degrees 2 to 4, 56 constants, on Whiting's
  ! half-width.  Arithmetic, as QLG-mid is.  The shapes are grouped by their power of G
  ! into three Horner chains in L -- pure L, G L, G^2 -- so the whole correction costs one
  ! evaluation of G, one of G^2, and thirteen FMAs.
  ! Maximum relative error, compiled in fp32: 3.1790e-5.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE qlg_top_sweep (a, v0, dv, n, K)
    real(4), intent(in)  :: a, v0, dv
    integer, intent(in)  :: n
    real(4), intent(out) :: K(n)
    integer :: f
    real(4) :: hw, ih, ih2, r, m, q, sq, x, bL, gG
    real(4) :: c2, c3, c4, c5
    real(4) :: d1, d2, d3, d4, d5
    real(4) :: e0, e1, e2
    real(4) :: v, u, L, G, P, G2
    !-----------------------------------------------------------------
    hw = 0.5346*a + sqrt(0.2166*a*a + LN2)
    ih = 1.0 / hw
    ih2 = ih * ih
    r  = a * ih
    m  = 1.0 - r
    q  = r * m
    sq = RPI * ih
    x  = 2.0*r - 1.0
    bL = sq * r * OVERPI
    gG = sq*(GNORM*m + q*(-0.8646222353 + x*(1.196227074 + x*(-1.946778655 + x*(0.8700489998)))))
    c2 = sq*q*(1.141897202 + x*(0.1391619146 + x*(0.001768314862)))
    c3 = sq*q*(2.27953124 + x*(-2.901926041 + x*(-0.3081531227 + x*(0.2136156261 + x*(-0.007055087946)))))
    c4 = sq*q*(6.46969986 + x*(-0.1679614633 + x*(-14.35427475 + x*(-5.158498764))))
    c5 = sq*q*(-0.985999763 + x*(-223.0961151 + x*(380.4154358 + x*(0.8698176742 + x*(-0.1589891016)))))
    d1 = sq*q*(-0.8064422011 + x*(-3.389981747 + x*(17.56308556 + x*(-11.88365173 + x*(0.2252746075)))))
    d2 = sq*q*(3.915934801 + x*(11.83234501 + x*(-70.72115326 + x*(48.13945007 + x*(-1.714033365)))))
    d3 = sq*q*(0.4049704671 + x*(28.79962349 + x*(-55.43619156 + x*(7.028235912))))
    d4 = sq*q*(-4.503441811 + x*(77.21445465 + x*(-117.9185104)))
    d5 = sq*q*(0.7967261076 + x*(111.5775375 + x*(-191.4626923)))
    e0 = sq*q*(-0.5365843177 + x*(0.85891819 + x*(2.652786493 + x*(-3.128877878 + x*(0.3580500185)))))
    e1 = sq*q*(-1.964906335 + x*(0.6071469188 + x*(11.63689137 + x*(-8.979340553 + x*(-0.1078078002)))))
    e2 = sq*q*(-5.386292934 + x*(-2.671384335 + x*(39.87611389 + x*(-27.97250366 + x*(1.403038502)))))

    do f = 1, n
      v = v0 + real(f-1, 4)*dv
      u = v*v * ih2
      L = 1.0 / (1.0 + u)
      call fexp2(u, G)
      P  = G*L
      G2 = G*G
      K(f) = bL*L + gG*G + L*L*(c2 + L*(c3 + L*(c4 + L*(c5)))) + P*(d1 + L*(d2 + L*(d3 + L*(d4 + L*(d5))))) + G2*(e0 + L*(e1 + L*(e2)))
    end do
  END SUBROUTINE qlg_top_sweep

  !=====================================================================================
  ! Build QLG-bot's two tables.  MUST be called once before any opacity or rate is
  ! accumulated: qlg_bot_sweep and qlg_bot_phi are PURE and cannot check, so an omitted
  ! call is not an error but a silently zero profile.  Idempotent, and a few
  ! milliseconds.  Entry points that must make it: SRC/composition_patch_mod.f90 after
  ! reading &linelist_params, and the standalone library's hop_init.
  ! Both tables are built in fp64 and rounded once, on the same grid u_i = (i/XSCL)^2.
  !-------------------------------------------------------------------------------------
  SUBROUTINE qlg_init
    integer :: i
    real(8) :: x
    !-----------------------------------------------------------------
    if (qlg_done) return
    do i = 0, NX
      x = real(i,8) / real(XSCL,8)
      TGX(i)  = real(2.0d0**(-x*x), 4)
      TH2X(i) = real(qlg_wr02(x*x), 4)
    end do
    qlg_done = .true.
  END SUBROUTINE qlg_init

  !=====================================================================================
  ! Z_0.2(u)/L^2 in fp64 -- the integrand of TH2X, and the definition of the basis
  ! function Z_r* at r* = 0.2 (Eq. 7 of the paper): the exact residual left by the two
  ! analytic limits at one fixed damping ratio, divided by the correction envelope q.
  ! The four anchor constants belong to that r* and are quoted to full fp64.  Evaluated
  ! against qlg_ref_f64 below, which is this file's own reference-grade Voigt in double.
  !-------------------------------------------------------------------------------------
  FUNCTION qlg_wr02 (u) result (w)
    real(8), intent(in) :: u
    real(8), parameter  :: A2  = 0.18741813853432507d0   ! a* : a/Gamma_exact(a) = 0.2
    real(8), parameter  :: GAM = 0.93709069267162530d0   ! Gamma_exact(a*)
    real(8), parameter  :: KP  = 0.81924853602069690d0   ! K(a*, 0) = erfcx(a*)
    real(8), parameter  :: ETA = 0.14697983683153040d0   ! a*/(sqrt(pi) Gamma^2 K_peak)
    real(8), parameter  :: QRS = 0.16d0                  ! r*(1 - r*)
    real(8) :: v, L, G, P, Kv(1), w
    !-----------------------------------------------------------------
    v = GAM * sqrt(u)
    L = 1.0d0 / (1.0d0 + u)
    G = exp(-DLN2 * u)
    call qlg_ref_f64 (A2, v, 0.0d0, 1, Kv)
    P = KP * (ETA*L + (1.0d0 - ETA)*G)
    w = (Kv(1) - P) / (QRS * KP) * (1.0d0 + u)**2
  END FUNCTION qlg_wr02

  !=====================================================================================
  ! Humlicek-Weideman, fp64: the two-region hybrid Schreier (2011) recommends -- Weideman's
  ! (1994) N = 32 rational approximation inside |v| + a < 15, and Humlicek's (1982) w4
  ! region-1 asymptote outside it -- and the reference every accuracy figure in the paper
  ! is measured against.  It is NOT w4: that is the kernel SMElib and Turbospectrum deploy
  ! and the one QLG-top replaces, and this hybrid is far more accurate and far dearer than
  ! either.  Public, so that a caller can check a member against it, and used
  ! internally by qlg_init to build TH2X.  It is not one of the four members -- it costs
  ! roughly an order of magnitude more than any of them and no point loop should call it.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE qlg_ref_f64 (a, v0, dv, n, K)
    real(8), intent(in)  :: a, v0, dv
    integer, intent(in)  :: n
    real(8), intent(out) :: K(n)
    integer :: f, j
    real(8) :: dr, nr0, v, av, di, q, zr, zi, pr, pi_, t
    real(8) :: ir, ii, sr, si, tr, ti, ur
    !-----------------------------------------------------------------------------------
    dr  = W32_LD + a
    nr0 = W32_LD - a
    do f = 1, n
      v  = v0 + real(f-1,8)*dv
      av = abs(v)
      if (av + a >= 15.0d0) then
        tr = a*a - av*av + 0.5d0
        ti = -2.0d0*a*av
        q  = 1.0d0 / (tr*tr + ti*ti)
        ur = ( a*tr + (-av)*ti) * q
        K(f) = ISQPI8 * ur
      else
        di = -v
        q  = 1.0d0 / (dr*dr + di*di)
        zr = (nr0*dr + v*di) * q
        zi = (v*dr - nr0*di) * q
        pr = W32_AD(1); pi_ = 0.0d0
        do j = 2, 32
          t   = pr*zr - pi_*zi + W32_AD(j)
          pi_ = pr*zi + pi_*zr
          pr  = t
        end do
        ir = dr*q;  ii = -di*q
        sr = ir*ir - ii*ii
        si = 2.0d0*ir*ii
        K(f) = 2.0d0*(pr*sr - pi_*si) + ISQPI8*ir
      end if
    end do
  END SUBROUTINE qlg_ref_f64

  !#####################################################################################
  !#####################################################################################
  !##
  !##  LOCAL ADDITION -- NOT PART OF THE PAPER'S MODULE.
  !##
  !##  Everything above this banner is qlg_mod.f90 as it ships with the paper, verbatim,
  !##  apart from the header-block lines marked "local addition" (the ISQPI parameter,
  !##  the public statements and the qlg_bot_coef_t type the routines below need).  Everything below is OpenMULTI3D's, and a
  !##  future re-copy from the paper's source should re-apply it rather than lose it.
  !##
  !##  WHY.  The module's four members are window sweeps, and OpenMULTI3D's hot path --
  !##  linelist_mod's add_opacity_equi -- calls qlg_bot_sweep that way, which is where
  !##  QLG-bot's advantage lies.  But the code has three colder Voigt paths that evaluate
  !##  one point at a time on abscissae that are not uniform ramps:
  !##
  !##    hprofl_mod        hydrogen fine-structure components, per component per point
  !##                      (through the split qlg_bot_coef / qlg_bot_eval form below)
  !##    rt_patch_mod      the NLTE rate integrals, over line%q
  !##    linelist_mod      add_opacity, the non-uniform-grid twin of add_opacity_equi,
  !##                      and the NLTE rate accumulation
  !##
  !##  qlg_bot_phi takes the Doppler-unit offset v; qlg_bot_profile takes the relative
  !##  detuning q with 1/df and applies the nu/nu0 Jacobian.
  !##
  !##  COST.  A one-point call redoes QLG-bot's per-line setup -- the half-width, r, and
  !##  the three coefficient polynomials, ~15 flops -- which the sweep form amortises over
  !##  a window.  It will not inline across the module boundary without -flto (see this
  !##  file's header, on fexp2).  A caller that evaluates many points at one damping on a
  !##  non-uniform abscissa, as the hydrogen profiles do per cell, splits the call instead:
  !##  qlg_bot_coef(a) once, then qlg_bot_eval(c, v) per point.  The setup there repeats
  !##  qlg_bot_sweep's, which stays verbatim above.  Do not "optimise" the hot path into
  !##  these; call qlg_bot_sweep directly, as it does.
  !##
  !##  qlg_init MUST have been called first, as for qlg_bot_sweep itself: these read the
  !##  same TGX/TH2X.  SRC/composition_patch_mod.f90 and the standalone library's
  !##  hop_init call it.
  !##
  !#####################################################################################
  !#####################################################################################

  !=====================================================================================
  ! K(a,v)/sqrt(pi) at one point -- the convention of the retired voigt_function2, whose
  ! call sites pass a = gamma/(4 pi dnud) and v = (nu - nu0)/dnud, both already scaled by
  ! the inverse Doppler width.  Peaks at 1/sqrt(pi) for a = v = 0, which is linelist_mod's
  ! voigt_peak and the bound its tolerance pre-filter relies on.
  !
  ! ELEMENTAL is legal over the PURE qlg_bot_sweep, and is used: rt_patch_mod and
  ! hprofl_mod both apply these to array sections.
  !-------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION qlg_bot_phi (a, v) result (h)
    real(4), intent(in) :: a, v
    real(4) :: h, K1(1)
    !-----------------------------------------------------------------
    call qlg_bot_sweep (a, v, 0.0, 1, K1)
    h = K1(1) * ISQPI
  END FUNCTION qlg_bot_phi

  !=====================================================================================
  ! qlg_bot_phi split at the damping: qlg_bot_coef holds qlg_bot_sweep's per-line setup
  ! for damping a, and qlg_bot_eval(c, v) is the profile at one point. It equals
  ! qlg_bot_phi(a, v) up to rounding: the compiler may contract the two forms into
  ! fused multiply-adds differently, which is why qlg_bot_phi still calls the sweep.
  !-------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION qlg_bot_coef (a) result (c)
    real(4), intent(in) :: a
    type(qlg_bot_coef_t) :: c
    real(4) :: d, r, q, sc, qs, x
    !-----------------------------------------------------------------
    d  = 1.0 + a*0.60602209212341318
    c%ih = d / (a*d + SQL2)
    r  = a * c%ih
    q  = r * (1.0 - r)
    sc = RPI * c%ih
    qs = q * sc
    x  = 2.0*r - 1.0
    c%bL = (OVERPI*r) * sc
    c%gG = (GNORM*(1.0 - r)) * sc + qs*(-0.064924393191013452 + x*(-0.11019326278794868))
    c%cZH = qs*(0.38494195930078284 + x*(-0.10567135279610407))
    c%cL2 = qs*(0.04054415746208978 + x*(0.16241029253712097))
  END FUNCTION qlg_bot_coef

  ELEMENTAL FUNCTION qlg_bot_eval (c, v) result (h)
    type(qlg_bot_coef_t), intent(in) :: c
    real(4), intent(in) :: v
    real(4) :: h, xs, u, L, G, W, p
    integer :: ig
    !-----------------------------------------------------------------
    xs = abs(v) * c%ih
    u  = xs * xs
    L  = 1.0 / (1.0 + u)
    ig = nint(min(xs*XSCL, RNX))
    G  = TGX(ig)
    W  = TH2X(ig)
    p = (c%cL2 + c%cZH*W)
    p = c%bL + L*p
    h = (G*(c%gG) + L*p) * ISQPI
  END FUNCTION qlg_bot_eval

  !=====================================================================================
  ! The convention of the retired voigt_profile: the same profile, but on the relative
  ! frequency q = (nu - nu0)/nu0 with the projected Doppler shift vpr subtracted and the
  ! inverse Doppler width idf applied here rather than by the caller, and normalised per
  ! unit frequency rather than per Doppler width -- hence the trailing idf*(q + 1), the
  ! nu/nu0 Jacobian.  Kept as a separate routine rather than folded into the call sites
  ! because that tail is exactly what distinguishes the NLTE rate integrals' convention
  ! from the opacity path's, and it is easy to drop.
  !-------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION qlg_bot_profile (a, q, vpr, idf) result (h)
    real(4), intent(in) :: a, q, vpr, idf
    real(4) :: h, v, K1(1)
    !-----------------------------------------------------------------
    v = abs(q - vpr) * idf
    call qlg_bot_sweep (a, v, 0.0, 1, K1)
    h = K1(1) * ISQPI * idf * (q + 1.0)
  END FUNCTION qlg_bot_profile

END MODULE qlg
