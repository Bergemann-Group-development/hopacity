  ! exp(-t) and 1 - exp(-t) in real(4), included into the contains section of each
  ! module that uses them, so that gfortran can inline them into the calling loop
  ! (it does not inline across modules).  Included by SRC/rt_solvers_mod.f90,
  ! HOP/bftrans_mod.f90 and HOP/level_mod.f90; the Makefiles list these includers.

  !===============================================================================
  !> exp(-t) and 1 - exp(-t), neither formed from the other.
  !>
  !> The reduction t = n*ln2 + r (ln2 split in two, as in cephes expf, |r| <= ln2/2)
  !> gives E = s*(1-m) and u0 = (1-s) + s*m, with s = 2^-n spliced into the exponent
  !> field and m = 1-exp(-r) = r*q(r), q a relative-error minimax fit of degree 5.
  !> E is formed as the product s*(1-m), so no intermediate is subnormal unless E
  !> is, and a flush-to-zero mode (ifort's default) cannot drop s*m from it.  The
  !> compensated forms need parentheses honoured; ifort builds pass
  !> -assume protect_parens for that.
  !> t is clamped at 87, so E does not fall below exp(-87) ~ 1.6e-38; this suits
  !> callers in which exp(-t) multiplies quantities of order one.  The exponent
  !> field limits t to t >= -87.
  !===============================================================================
  pure subroutine exp_neg (t, E, u0)
    real(4), intent(in) :: t
    real(4), intent(out) :: E, u0
    ! ln2 split so that n*EXP_C1 is exact; EXP_Q is q, lowest order first
    real(4), parameter :: EXP_LOG2E = 1.44269504, EXP_C1 = 0.693359375, EXP_C2 = -2.12194440e-4
    real(4), parameter :: EXP_Q(0:5) = [ &
      +1.000000000e+00, -5.000000000e-01, +1.666666716e-01, -4.167019576e-02, &
      +8.353385143e-03, -1.351550105e-03 ]
    real(4) :: tc, n, r, m, s
    !-----------------------------------------------------------------------------
    tc = min(t, 87.0)
    n = anint(tc*EXP_LOG2E)
    r = (tc - n*EXP_C1) - n*EXP_C2
    m = r*(EXP_Q(0) + r*(EXP_Q(1) + r*(EXP_Q(2) + r*(EXP_Q(3) + r*(EXP_Q(4) + r*EXP_Q(5))))))
    s = transfer(ishft(127 - int(n), 23), 1.0)
    E = s*(1.0 - m)
    u0 = (1.0 - s) + s*m
  end subroutine exp_neg

  !===============================================================================
  !> exp_neg for any t, where exp(-t) may multiply large quantities: E is 0 above
  !> t = 87, by bit mask so that the calling loop stays branch-free, and t is
  !> clamped at -87 below, where E saturates at exp(87) ~ 6e37.
  !===============================================================================
  pure subroutine exp_neg0 (t, E, u0)
    real(4), intent(in) :: t
    real(4), intent(out) :: E, u0
    !-----------------------------------------------------------------------------
    call exp_neg(max(t, -87.0), E, u0)
    E = transfer(iand(transfer(E, 1), merge(-1, 0, t <= 87.0)), 1.0)
  end subroutine exp_neg0

  !===============================================================================
  !> exp(x) by exp_neg0: 0 below x = -87, saturating at exp(87) above x = 87.
  pure real(4) function exp4 (x)
    real(4), intent(in) :: x
    real(4) :: u0
    !-----------------------------------------------------------------------------
    call exp_neg0(-x, exp4, u0)
  end function exp4

  !===============================================================================
  !> y = exp4(y) elementwise.  Kept out of line, so that exp4 is inlined into
  !> this loop also where the caller is too large for gfortran's early inliner.
  pure subroutine exp4_array (y)
!GCC$ ATTRIBUTES noinline :: exp4_array
    real(4), intent(inout) :: y(:)
    integer :: i
    !-----------------------------------------------------------------------------
    do i=1, size(y)
      y(i) = exp4(y(i))
    end do
  end subroutine exp4_array
