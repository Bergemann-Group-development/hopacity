!=======================================================================================
! voigt_av_mod -- branchless, vectorisable fp32 Voigt function K(a,v) = Re[w(v + i a)].
!
! This is a Fortran port of voigt_av.c (fitted in voigt.ipynb against scipy's wofz):
!
!   K = (sqrt(ln2)/hw) * E * g * P_A(r,x)  +  a/(sqrt(pi)(hw^2+v^2)) * (1 + (hw/(v+hw)) P_B(r,x))
!
! with hw = hwhm(a) from a polynomial front-end in rho = a/(sqrt(ln2)+a), r = a/hw,
! x = v/(v+hw), E = p(min(v^2,18)/32)^32 ~ exp(-v^2) (minimax poly then five squarings),
! and g = min(1, 100(1-x)^2). P_A and P_B are monomial tables in xh = 2x-1 (rows) by
! mapped r (columns); the far wing is the exact Lorentz asymptote.
!
! Why it exists: the incumbent voigt_function2 (bbtrans_mod) is a four-way data-dependent
! branch tree containing a complex(8) division and a real exp(). Since `a` is fixed per
! line while `v` sweeps across the profile window, the branch taken varies point-to-point
! *within* the frequency loop, so that loop cannot vectorise -- it costs ~170 cycles per
! evaluation and is ~32% of a full SED run's wall time.
!
! The split below mirrors the two-phase structure the C header calls out:
!   voigt_av_setup  ~130 flops, once per (cell, line)   -- collapses the r-direction of
!                   the coefficient tables into cA/cB by Horner
!   voigt_av_eval   ~48 flops + 2 divisions per point   -- branchless except two min()
!
! Accuracy (fp32 vs double wofz, v in [0,200], relative where K > 1e-3):
!   a=1e-6: 2.6e-4   a=1e-2: 9.4e-5   a=0.21: 1.4e-4   a=1.0: 1.9e-4   a=86: 1.8e-4
! i.e. uniformly ~2e-4. For comparison voigt_function2 is better below a ~ 1e-3 (1.7e-7
! at a=1e-6) but degrades to 3.3e-2 at (a=0.21, v=3.0) and 6.0e-2 at (a=1.0, v=2.5) --
! both on the `u = v+a > 3.2` seam where it switches to the far-wing asymptote early --
! and 1-2e-2 across the core for a > 1.4. So for damped lines this is the more accurate
! of the two; its only regression is the deep Lorentz wing at a <~ 1e-5, where the
! absolute error ~1e-8 sits far below any caller's opacity tolerance.
!
! NOTE ON NORMALISATION: this returns K = H(a,v), whereas voigt_function2 returns
! H(a,v)/sqrt(pi). Callers substituting one for the other must fold in voigt_peak.
!=======================================================================================
MODULE voigt_av_mod
  implicit none
  private

  public :: voigt_av_t, voigt_av_setup, voigt_av_eval, voigt_av_block, voigt_av

  ! per-`a` state; built once per line, consumed by every point of that line's window
  type :: voigt_av_t
    real(4) :: cA(10)      ! P_A monomial coefficients in xh, r-direction already collapsed
    real(4) :: cB(12)      ! P_B ditto
    real(4) :: hw          ! half width at half maximum of K(a,.)
    real(4) :: hw2         ! hw*hw
    real(4) :: q           ! sqrt(ln2)/hw
    real(4) :: arpi        ! a/sqrt(pi)
  end type

  real(4), parameter :: SQL2 = 0.8325546112   ! sqrt(ln 2)
  real(4), parameter :: RPI  = 1.7724538509   ! sqrt(pi)

  real(4), parameter :: P_HWHM(6) = [ &
    +6.082167625e-01, -5.944246650e-01, -2.180101573e-01, &
    -2.297060937e-01, -3.986304104e-01, +8.325546384e-01 ]

  real(4), parameter :: EXP_P(7) = [ &
    +1.000000000e+00, -9.999999404e-01, +4.999980927e-01, -1.666474938e-01, &
    +4.157148302e-02, -8.083330467e-03, +1.051157247e-03 ]

  ! NA/NB keep the C row/column sense (NA[10][5] -> NA(10,5)), i.e. the *coefficient*
  ! index runs contiguously and the Horner index is the slow one. That is deliberate:
  ! the r-direction collapse in voigt_av_setup is independent across coefficients, so
  ! laying them out contiguously lets that collapse vectorise as whole-array statements.
  ! It matters -- profile windows average only ~7.6 points wide, so a scalar setup would
  ! be a ~35% tax on the per-point work rather than a negligible prologue.
  real(4), parameter :: NA(10,5) = reshape([ &
    -2.841740549e-01, +3.051739335e-01, +4.145147204e-01, -7.996425629e-01, +2.906724811e-01, &
    -7.885839939e-01, -9.768277407e-01, +1.374700189e+00, +3.868988454e-01, -1.132824540e+00, &
    +9.390527725e+00, -1.686799812e+01, -5.940451622e+00, +1.453922081e+01, -4.597479820e+00, &
    -2.331789207e+01, +5.113978577e+01, -9.649711609e+00, -5.084135056e+01, +2.910062027e+01, &
    +3.456605530e+01, -6.646157837e+01, +1.826349449e+01, +6.226839066e+01, -5.100187683e+01, &
    -2.664812565e+00, +3.581811142e+01, +1.443977833e+01, -3.797651672e+01, -1.475073051e+01, &
    -6.289942741e+00, +1.028932953e+02, -8.747564697e+01, -1.005940170e+02, +8.826245117e+01, &
    +6.660982513e+01, -1.024929123e+02, -6.492916107e+01, +1.030761490e+02, +5.133430481e+00, &
    -6.415218830e+00, -1.060191574e+02, +1.169886093e+02, +1.072529602e+02, -1.016141739e+02, &
    -5.501876068e+01, +1.130189590e+02, +7.080995083e+00, -1.069713516e+02, +4.540658569e+01 ], [10,5], order=[2,1])

  real(4), parameter :: NB(12,7) = reshape([ &
    +3.548822165e+00, -8.368619919e+00, +8.148386002e+00, -3.699010134e+00, +8.234729618e-02, +2.067441344e-01, +8.166331798e-02, &
    +3.515838385e-01, -7.093361616e-01, +5.674955994e-02, +1.574641019e-01, +2.679947913e-01, +2.967274785e-01, -4.223370850e-01, &
    -4.408294296e+01, +1.394360352e+02, -1.559816437e+02, +6.414968109e+01, -1.235230088e+00, -1.768387318e+00, -5.161782503e-01, &
    +1.128925323e+02, -3.631920776e+02, +3.891399231e+02, -1.390268860e+02, -5.439591408e-01, -1.906112671e+00, +2.655450821e+00, &
    -1.026814194e+02, +2.599161072e+02, -1.828707428e+02, +1.412961769e+01, +4.386076927e+00, +5.116020203e+00, +1.982548714e+00, &
    -3.981498718e+01, +2.593503723e+02, -4.074282227e+02, +1.900587769e+02, -9.019944072e-01, +6.237632275e+00, -7.591190338e+00, &
    +1.645718689e+02, -4.675358887e+02, +3.900527344e+02, -6.928217316e+01, -7.081229687e+00, -6.905075550e+00, -3.742149830e+00, &
    -8.793880463e+01, +3.058077812e+00, +2.774294739e+02, -1.967609100e+02, +4.222378731e+00, -1.096558189e+01, +1.113395405e+01, &
    -8.026850128e+01, +4.016867981e+02, -5.183965454e+02, +1.835954895e+02, +5.498904705e+00, +4.544492245e+00, +3.234449148e+00, &
    +9.395753479e+01, -2.565823669e+02, +2.116531372e+02, -4.567527390e+01, -4.965072155e+00, +9.485515594e+00, -8.033469200e+00, &
    -7.839607716e+00, +6.462976933e+00, -3.202336788e+00, +8.511311531e+00, -1.645035982e+00, -1.190568209e+00, -1.049151540e+00, &
    -1.272360420e+01, +2.649419022e+01, -8.601893425e+00, -6.155086517e+00, +1.941450000e+00, -3.149172783e+00, +2.247422218e+00 ], [12,7], order=[2,1])

contains

  !=====================================================================================
  ! Per-`a` setup: collapse the r-direction of the coefficient tables. Call once per
  ! line, before sweeping its frequency window. ~130 flops.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE voigt_av_setup (a, s)
    real(4),          intent(in)  :: a
    type(voigt_av_t), intent(out) :: s
    integer :: k
    real(4) :: rho, hw, invhw, r, rh
    !-----------------------------------------------------------------------------------
    rho = a / (SQL2 + a)
    hw  = P_HWHM(1)
    do k = 2, 6
      hw = hw*rho + P_HWHM(k)
    end do
    hw = hw + a

    invhw = 1.0 / hw
    r     = a * invhw
    rh    = min((r - 0.01)*2.0204061 - 1.0, 1.0)

    ! Horner in rh, all coefficients at once. Written as one expression per table with
    ! literal column indices on purpose: with the Horner step in a DO loop, NA(:,k) is a
    ! variable-offset slice and gfortran emits a __builtin_memcpy from the constant table
    ! that clobbers memory and blocks vectorisation entirely. Spelled out like this the
    ! columns are compile-time addresses and the collapse SLP-vectorises to ~30 vector
    ! FMAs instead of ~112 scalar ones -- which matters, because this runs once per
    ! (cell, line) for an average of only ~7.6 evaluated points.
    s%cA = (((NA(:,5)*rh + NA(:,4))*rh + NA(:,3))*rh + NA(:,2))*rh + NA(:,1)
    s%cB = (((((NB(:,7)*rh + NB(:,6))*rh + NB(:,5))*rh + NB(:,4))*rh + NB(:,3))*rh &
             + NB(:,2))*rh + NB(:,1)

    s%hw   = hw
    s%hw2  = hw*hw
    s%q    = SQL2 * invhw
    s%arpi = a / RPI
  END SUBROUTINE voigt_av_setup

  !=====================================================================================
  ! Per-point evaluation, given the setup state. Branchless apart from two min(), which
  ! lower to NEON fmin -- this is what lets the caller's frequency loop vectorise.
  ! ~48 flops + 2 divisions. v must be non-negative (K is even in v).
  !-------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION voigt_av_eval (s, v) result (K)
    type(voigt_av_t), intent(in) :: s
    real(4),          intent(in) :: v
    real(4) :: K
    integer :: j
    real(4) :: inv, x, xm, v2, t, p, E, g, b2, xh, PA, PB
    !-----------------------------------------------------------------------------------
    inv = 1.0 / (v + s%hw)
    x   = v * inv
    xm  = s%hw * inv
    v2  = v * v

    ! E ~ exp(-v^2) as p(min(v^2,18)/32)^32: minimax poly then five squarings
    t = min(v2, 18.0) * 0.03125
    p = EXP_P(7)
    do j = 6, 1, -1
      p = p*t + EXP_P(j)
    end do
    E = p*p
    E = E*E
    E = E*E
    E = E*E
    E = E*E

    g  = min(1.0, 100.0*xm*xm) * s%q
    b2 = s%arpi / (s%hw2 + v2)
    xh = 2.0*x - 1.0

    PA = s%cA(10)
    do j = 9, 1, -1
      PA = PA*xh + s%cA(j)
    end do
    PB = s%cB(12)
    do j = 11, 1, -1
      PB = PB*xh + s%cB(j)
    end do

    K = E*g*PA + b2*(1.0 + xm*PB)
  END FUNCTION voigt_av_eval

  !=====================================================================================
  ! Whole-window sweep: K(f) = voigt_av_eval(s, |v0 + (f-1)*dv|), f = 1..n.
  !
  ! This is the form the opacity kernels should use. Writing the per-point body out
  ! inline (rather than calling the elemental above) is deliberate: it leaves gfortran a
  ! straight-line branchless loop body with no call, which is what it needs to emit NEON.
  ! The abscissa is computed as v0 + (f-1)*dv rather than accumulated as v = v + dv --
  ! the latter is a loop-carried dependency and blocks vectorisation entirely.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE voigt_av_block (s, v0, dv, n, K)
    type(voigt_av_t), intent(in)  :: s
    real(4),          intent(in)  :: v0, dv
    integer,          intent(in)  :: n
    real(4),          intent(out) :: K(n)
    integer :: f, j
    real(4) :: v, inv, x, xm, v2, t, p, E, g, b2, xh, PA, PB
    !-----------------------------------------------------------------------------------
    do f = 1, n
      v   = abs(v0 + real(f-1)*dv)
      inv = 1.0 / (v + s%hw)
      x   = v * inv
      xm  = s%hw * inv
      v2  = v * v

      t = min(v2, 18.0) * 0.03125
      p = EXP_P(7)
      do j = 6, 1, -1
        p = p*t + EXP_P(j)
      end do
      E = p*p
      E = E*E
      E = E*E
      E = E*E
      E = E*E

      g  = min(1.0, 100.0*xm*xm) * s%q
      b2 = s%arpi / (s%hw2 + v2)
      xh = 2.0*x - 1.0

      PA = s%cA(10)
      do j = 9, 1, -1
        PA = PA*xh + s%cA(j)
      end do
      PB = s%cB(12)
      do j = 11, 1, -1
        PB = PB*xh + s%cB(j)
      end do

      K(f) = E*g*PA + b2*(1.0 + xm*PB)
    end do
  END SUBROUTINE voigt_av_block

  !=====================================================================================
  ! Convenience one-shot form (setup + eval). Only for callers that evaluate a single
  ! point per `a` -- in a profile sweep, hoist voigt_av_setup instead.
  !-------------------------------------------------------------------------------------
  PURE FUNCTION voigt_av (a, v) result (K)
    real(4), intent(in) :: a, v
    real(4) :: K
    type(voigt_av_t) :: s
    !-----------------------------------------------------------------------------------
    call voigt_av_setup (a, s)
    K = voigt_av_eval (s, abs(v))
  END FUNCTION voigt_av

END MODULE voigt_av_mod
