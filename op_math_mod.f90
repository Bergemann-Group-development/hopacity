MODULE op_math_mod
  USE prec_mod
  implicit none

  interface poly1d
    module procedure poly3, poly4, poly5, poly6
  end interface poly1d

  interface find_closest
    module procedure find_closest4, find_closest8
  end interface find_closest

  interface searchsorted
    module procedure searchsorted4, searchsorted8
  end interface searchsorted

  interface BinarySearch
    module procedure BinarySearch4, BinarySearch8
  end interface BinarySearch
contains

  !====================================================================================
  PURE FUNCTION option(bool, default) result(out)
    logical, intent(in), optional :: bool, default
    logical :: out
    !---------------------------------------------------------------------------------
    if (present(bool)) then
      out = bool
    else
      if (present(default)) then
        out = default
      else
        out = .false.
      endif
    endif
  END FUNCTION option

  !====================================================================================
  pure FUNCTION find_closest4(arr, val) result(i)
    real(4), intent(in) :: arr(:), val
    integer :: i
    !----------------------------------------------------------------------------
    i = searchsorted4(arr, val, outbounds_ok=.true.)
    if (i == 0) then
      i = 1
    elseif (i < size(arr)) then
      if (val - arr(i) > arr(i+1) - val) then
        i = i + 1
      endif
    endif
  end FUNCTION find_closest4

  !====================================================================================
  pure FUNCTION find_closest8(arr, val) result(i)
    real(8), intent(in) :: arr(:), val
    integer :: i
    !----------------------------------------------------------------------------
    i = searchsorted8(arr, val, outbounds_ok=.true.)
    if (i == 0) then
      i = 1
    elseif (i < size(arr)) then
      if (val - arr(i) > arr(i+1) - val) then
        i = i + 1
      endif
    endif
  end FUNCTION find_closest8

  !====================================================================================
  ! Finds index of element after which value should be inserted to maintain order.
  !====================================================================================
  pure FUNCTION searchsorted4(arr, val, outbounds_ok) result(i)
    real(4), intent(in) :: arr(:), val
    logical, optional, intent(in) :: outbounds_ok
    logical :: ob_ok, above, below
    integer :: i, n
    !----------------------------------------------------------------------------
    if (present(outbounds_ok)) then
      ob_ok = outbounds_ok
    else
      ob_ok = .false.
    endif

    n = size(arr)
    below = val < arr(1)
    above = val > arr(n)

    if (below .or. above) then
      if (ob_ok) then
        if (below) i = 0
        if (above) i = n
      else
        i = -1
      endif
    else
      i = BinarySearch4(arr, 1, n, val)
    endif
  end FUNCTION searchsorted4

  !====================================================================================
  ! Finds index of element after which value should be inserted to maintain order.
  !====================================================================================
  pure FUNCTION searchsorted8(arr, val, outbounds_ok) result(i)
    real(8), intent(in) :: arr(:), val
    logical, optional, intent(in) :: outbounds_ok
    logical :: ob_ok, above, below
    integer :: i, n
    !----------------------------------------------------------------------------
    if (present(outbounds_ok)) then
      ob_ok = outbounds_ok
    else
      ob_ok = .false.
    endif

    n = size(arr)
    below = val < arr(1)
    above = val > arr(n)

    if (below .or. above) then
      if (ob_ok) then
        if (below) i = 0
        if (above) i = n
      else
        i = -1
      endif
    else
      i = BinarySearch8(arr, 1, n, val)
    endif
  end FUNCTION searchsorted8

  !====================================================================================
  pure recursive FUNCTION BinarySearch4(arr, l, r, x) result(idx)
    real(4), intent(in) :: arr(:)
    real(4), intent(in) :: x
    integer, intent(in) :: l, r
    integer :: idx, mid
    !----------------------------------------------------------------------------
    if (r - l > 1) then
      mid = l + (r - l) / 2

      if(arr(mid) >= x) then
        idx = BinarySearch4(arr, l, mid, x)

      elseif(arr(mid) < x) then
        idx = BinarySearch4(arr, mid, r, x)

      endif
    else
      idx = l
    endif
  end FUNCTION BinarySearch4

  !====================================================================================
  pure recursive FUNCTION BinarySearch8(arr, l, r, x) result(idx)
    real(8), intent(in) :: arr(:)
    real(8), intent(in) :: x
    integer, intent(in) :: l, r
    integer :: idx, mid
    !----------------------------------------------------------------------------
    if (r - l > 1) then
      mid = l + (r - l) / 2

      if(arr(mid) >= x) then
        idx = BinarySearch8(arr, l, mid, x)

      elseif(arr(mid) < x) then
        idx = BinarySearch8(arr, mid, r, x)

      endif
    else
      idx = l
    endif
  end FUNCTION BinarySearch8

  !===============================================================================
  PURE FUNCTION linspace (x1, x2, n) result(xx)
    real(4), intent(in) :: x1, x2
    integer, intent(in) :: n
    real(4), dimension(n) :: xx
    real(4) :: dx
    integer :: i
    !-----------------------------------------------------------------------------
    dx = (x2 - x1) / real(n-1)
    xx = [(x1 + i*dx, i=0, n-1)]
  END FUNCTION linspace

  !===============================================================================
  PURE FUNCTION normspace (n) result(xx)
    integer, intent(in) :: n
    real(4), dimension(n) :: xx
    real(4) :: dx
    integer :: i
    !-----------------------------------------------------------------------------
    dx = 1. / real(n)
    xx = [(i*dx, i=0, n-1)]
  END FUNCTION normspace

  !====================================================================================
  pure SUBROUTINE get_interval(arr, lower, upper, lo, up, limit)
    real(8), intent(in) :: arr(:), lower, upper
    integer, intent(out) ::lo, up
    logical, intent(in) :: limit
    integer :: n
    !----------------------------------------------------------------------------
    lo = searchsorted(arr, lower, outbounds_ok=.true.) + 1
    up = searchsorted(arr, upper, outbounds_ok=.true.)
    n = size(arr)
    if (lo > 1) then
      do while(arr(lo-1) == lower)
        lo = lo - 1
        if (lo == 1) exit
      end do
    endif
    if (up < n) then
      do while(arr(up+1) == upper)
        up = up + 1
        if (up >= n) exit
      end do
    endif
    if (limit) then
      lo = max(1, lo)
      up = min(n, up)
    endif
  end SUBROUTINE get_interval

  !====================================================================================
  elemental real(dp) FUNCTION poly3(x, a, b, c) result(poly)
    real(dp), intent(in) :: x, a, b, c
    !----------------------------------------------------------------------------
    poly = a*x*x + b*x + c
  end FUNCTION poly3

  !====================================================================================
  elemental real(dp) FUNCTION poly4(x, a, b, c, d) result(poly)
    real(dp), intent(in) :: x, a, b, c, d
    !----------------------------------------------------------------------------
    poly = x*poly3(x, a, b, c) + d
  end FUNCTION poly4

  !====================================================================================
  elemental real(dp) FUNCTION poly5(x, a, b, c, d, e) result(poly)
    real(dp), intent(in) :: x, a, b, c, d, e
    !----------------------------------------------------------------------------
    poly = x*poly4(x, a, b, c, d) + e
  end FUNCTION poly5

  !====================================================================================
  elemental real(dp) FUNCTION poly6(x, a, b, c, d, e, f) result(poly)
    real(dp), intent(in) :: x, a, b, c, d, e, f
    !----------------------------------------------------------------------------
    poly = x*poly5(x, a, b, c, d, e) + f
  end FUNCTION poly6

  !====================================================================================
  SUBROUTINE poly3d(x, p, res)
    real(dp), intent(in), dimension(:,:,:) :: x
    real(dp), intent(in), dimension(:) :: p
    real(dp), intent(out), dimension(:,:,:) :: res
    integer :: i
    !----------------------------------------------------------------------------
    res = p(1)
    do i = 2, size(p)
      res = res * x + p(i)
    end do
  END SUBROUTINE poly3d

  !===============================================================================
  elemental real(dp) FUNCTION solve_quadratic(p, q) result(x)
    real(dp), intent(in) :: p, q
    real(dp) :: cam, root
    !----------------------------------------------------------------------------
    if (p > 0) then
      root = sqrt(abs(p*p - 4.d0*q))
      x = - 2.d0*q / (p + root)
    else
      cam = - 0.5d0*p
      root = sqrt(abs(cam*cam - q))
      x = cam + root
    endif
  end FUNCTION solve_quadratic

  !===============================================================================
  !  finds solution of system of linear equations
  !  with gaussian elimination with pivoting
  !: lineq  90-06-05  new routine: (martin j stift)
  !===============================================================================
  SUBROUTINE lineq(a,b,n)
    integer :: n
    real(dp), dimension(n,n) :: a
    real(dp), dimension(n)   :: b

    real(dp), dimension(n)  :: c
    integer, dimension(n) :: icol
    integer :: i, ibeg, jbeg, ima, jma, j, it, n1, jmin, i1, imin
    real(dp)  :: ama, temp, quot
    !-----------------------------------------------------------------------------
    ! initialize column count and starting points

    do i = 1,n
      icol(i) = i
    end do

    ibeg = 1
    jbeg = 1

    ! determine pivot
    do while (ibeg.lt.n)
      ama = abs(a(ibeg,jbeg))
      ima = ibeg
      jma = jbeg
      do i = ibeg,n
        do j = jbeg,n
          if(abs(a(i,j)).le.ama) cycle
          ama = abs(a(i,j))
          ima = i
          jma = j
        end do
      end do

      ! order matrix depending on pivot

      do i = 1,n
        temp = a(i,jma)
        a(i,jma) = a(i,jbeg)
        a(i,jbeg) = temp
      end do

      do j = jbeg,n
        temp = a(ima,j)
        a(ima,j) = a(ibeg,j)
        a(ibeg,j) = temp
      end do

      temp = b(ima)
      b(ima) = b(ibeg)
      b(ibeg) = temp
      it = icol(jbeg)
      icol(jbeg) = icol(jma)
      icol(jma)= it
      ibeg = ibeg + 1
      jbeg = jbeg + 1

      ! eliminate

      imin = ibeg - 1
      jmin = jbeg - 1
      do i = ibeg,n
        quot = a(i,jmin) / a(imin,jmin)
        do j = jbeg,n
          a(i,j) = a(i,j) - quot * a(imin,j)
        end do
        b(i) = b(i) - quot * b(imin)
      end do
    end do
    !if (ibeg.lt.n) goto 20

    ! determine coefficients

    b(n) = b(n) / a(n,n)
    n1 = n - 1
    do i = n1,1,-1
      i1 = i + 1
      do j = n,i1,-1
        b(i) = b(i) - a(i,j) * b(j)
      end do
      b(i) = b(i) / a(i,i)
    end do

    ! reorder coefficients

    do i = 1,n
      c(icol(i))=b(i)
    end do

    do i = 1,n
      b(i) = c(i)
    end do

  end SUBROUTINE lineq
  !---------------------------------------------------------------------

  !===============================================================================
  !  solves the equation  system  a*x=b.
  !  when newmat=true, the system is rearranged into u*x=l*b, where u
  !  is upper and l is lower triangular. these are then reused in later
  !  calls with newmat=false and new right hand sides b. the solution
  !  vector is returned in b. no pivoting, i.e. the matrix a is assumed
  !  to have nonzero diagonal elements.
  !
  !  coded by: a. nordlund (feb-1979)
  !
  !  this is a modified version of eqsyst which tests for zero elements
  !  below the diagonal and also stops at the last non-zero element above
  !  the diagonal. considerable savings are obtained for loose matrices.
  !
  !  this is a column oriented version (m. carlsson jan-1986)
  !  temporary scalars bl, all, alm and bk are used to show the compiler
  !  that there is no vector dependency in the innermost do-loop
  !===============================================================================
  SUBROUTINE eqsyst(a,b,n)
    integer :: n
    real(8) :: a(n,n), b(n)
    real(8) :: all, bl, alm, bk
    integer :: lastn(n), l, k, m
    !-----------------------------------------------------------------------------
    ! find the last non-zero element in each column
    do l=1,n
      do k=n,l+1,-1
        if(a(k,l)/=0.0) exit
      end do
      lastn(l)=k
    end do

    !   column loop: eliminate elements below the diagonal in column l.
    do l=1,n-1
      !  store -a(k,l)/a(l,l) in element a(k,l)
      !  multiply right hand side with -a(k,l)/a(l,l)
      all=a(l,l)
      bl=b(l)
      do k=l+1,lastn(l)
        a(k,l)=-a(k,l)/all
        b(k)=b(k)+a(k,l)*bl
      end do
      !  add fraction -a(k,l)/a(l,l) of row l to row k.
      !  in each column go through all rows
      do m=l+1,n
        if(a(l,m)/=0.0) then
          alm=a(l,m)
          lastn(m)=max(lastn(l),lastn(m))
          do k=l+1,lastn(l)
            a(k,m)=a(k,m)+a(k,l)*alm
          end do
        end if
      end do
    end do

    ! backsubstitute
    do k=n,1,-1
      bk=b(k)
      do l=k+1,n
        bk=bk-a(k,l)*b(l)
      end do
      b(k)=bk/a(k,k)
    end do
  END SUBROUTINE eqsyst

  !===============================================================================
  SUBROUTINE lineq_no_pivot (a, b, n)
    integer :: n
    real(8):: a(n,n), b(n)
    integer :: i, j, k
    !-----------------------------------------------------------------------------
    do k=1,n-1
      do i=k+1,n
        a(i,k) = a(i,k)/a(k,k)
      end do
      do j=k+1,n
        if (a(k,j) /= 0) then
          do i=k+1,n
            a(i,j) = a(i,j) - a(i,k)*a(k,j)
          end do
        endif
      end do
    end do
    do i=1,n
      do j=1,i-1
        b(i) = b(i) - a(i,j)*b(j)
      end do
    end do
    do i=n,1,-1
      do j=i+1,n
        b(i) = b(i) - a(i,j)*b(j)
      end do
      b(i) = b(i)/a(i,i)
    end do
  END SUBROUTINE lineq_no_pivot

  !===============================================================================
  SUBROUTINE AINV2(A,M)
    implicit none
    integer :: M
    real(8), dimension(M,M) :: A
    integer, dimension(M) :: IND
    real(8), dimension(M) :: C
    real(8) :: AMAX,STO,W1,W
    integer :: I, I4, MM, J, ISTO, K, J1, I1, J2, I2
    !***** THIS SUBROUTINE EVALUATES THE INVERSE OF A
    !***** SQUARE M*M MATRIX A
    !***** THE DIMENSIONS OF THE ARRAYS MAY BE INTEGER VARIABLES WHEN USED*********
    !***** ON THE IBM 7094,BUT THEY MUST BE INTEGER CONSTANTS ON THE IBM 1130******
    !-----------------------------------------------------------------------------
    AMAX=0.0

    do I=1,M
      IND(I)=I
      if(abs(A(I,1))-AMAX > 0) then
        AMAX=abs(A(I,1))
        I4=I
      endif
    end do

    MM=M-1

    do J=1,MM
      if(I4-J > 0) then
        ISTO=IND(J)
        IND(J)=IND(I4)
        IND(I4)=ISTO
        do K=1,M
          STO=A(I4,K)
          A(I4,K)=A(J,K)
          A(J,K)=STO
        end do
      endif
      AMAX=0.0
      J1=J+1
      do I=J1,M
        A(I,J)=A(I,J)/A(J,J)
        do K=J1,M
          A(I,K)=A(I,K)-A(I,J)*A(J,K)
          IF (K-J1 <= 0) then
            IF(abs(A(I,K))-AMAX > 0) then
              AMAX=abs(A(I,K))
              I4=I
            endif
          endif
        end do
      end do
    end do

    do I1=1,MM
      I=M+1-I1
      I2=I-1
      do J1=1,I2
        J=I2+1-J1
        J2=J+1
        W1=-A(I,J)
        IF(I2-J2 >= 0) then
          do K=J2,I2
            W1=W1-A(K,J)*C(K)
          end do
        endif
        C(J)=W1
      end do
      do K=1,I2
        A(I,K)=C(K)
      end do
    end do

    do I1=1,M
      I=M+1-I1
      I2=I+1
      W=A(I,I)
      do J=1,M
        if (I-J < 0) then
          W1=0.0
        elseif(I-J == 0) then
          W1=1.0
        else
          W1=A(I,J)
        endif
        if (I1-1 > 0) then
          do K=I2,M
            W1=W1-A(I,K)*A(K,J)
          end do
        endif
        C(J)=W1
      end do
      do J=1,M
        A(I,J)=C(J)/W
      end do
    end do

    do I=1,M
      do while (IND(I)-I /=0)
        J=IND(I)
        do K=1,M
          STO=A(K,I)
          A(K,I)=A(K,J)
          A(K,J)=STO
        end do
        ISTO=IND(J)
        IND(J)=J
        IND(I)=ISTO
      end do
    end do
    return
  END SUBROUTINE

  !===============================================================================
  ! Orthogonal PCA basis of 2D scatter data (x1,x2), e.g. (lnT, lnNe) samples.
  ! Returns the mean, rotation into principal-axis coordinates, and the origin/
  ! spacing/node count of a rectangular grid laid out along those axes with
  ! projected resolution d1,d2 (in x1,x2 units) and npad extra nodes past the
  ! data hull on each end. Caller passes x1,x2 pre-transformed (e.g. logged)
  ! so the transform is computed once rather than repeated per moment/extent.
  !===============================================================================
  SUBROUTINE pca_grid_2d(x1, x2, d1, d2, npad, mu, rot, c0, dc, n1, n2, use_pca)
    real(8), intent(in)  :: x1(:), x2(:)     ! sample coordinates, e.g. (lnT, lnNe)
    real(8), intent(in)  :: d1, d2           ! requested resolution along x1, x2
    integer, intent(in)  :: npad             ! extra nodes past the data hull, each end
    real(8), intent(out) :: mu(2), rot(2,2), c0(2), dc(2)
    integer, intent(out) :: n1, n2
    logical, intent(in), optional :: use_pca ! .false. -> axis-aligned (no rotation)
    real(8) :: a, b, d, tr, disc, l1, l2, tiny_b, vn
    real(8) :: pmin1, pmax1, pmin2, pmax2
    real(8), dimension(size(x1)) :: p1, p2
    integer :: n, np
    !-------------------------------------------------------------------------
    n = size(x1)

    ! ---- mean and 2x2 covariance of (x1, x2) ----------------------------------
    mu(1) = sum(x1) / n
    mu(2) = sum(x2) / n
    a = sum((x1 - mu(1))**2) / n                      ! var(x1)
    d = sum((x2 - mu(2))**2) / n                      ! var(x2)
    b = sum((x1 - mu(1)) * (x2 - mu(2))) / n          ! cov

    ! ---- closed-form symmetric 2x2 eigendecomposition -------------------------
    tr   = a + d
    disc = sqrt(max(((a - d)*0.5d0)**2 + b*b, 0.d0))
    l1   = tr*0.5d0 + disc                         ! larger eigenvalue -> PC1
    l2   = tr*0.5d0 - disc
    tiny_b = 1.d-30 * (abs(a) + abs(d) + 1.d-300)
    if (abs(b) > tiny_b) then
      rot(:,1) = [b, l1 - a];  vn = sqrt(sum(rot(:,1)**2));  rot(:,1) = rot(:,1)/vn
      rot(:,2) = [b, l2 - a];  vn = sqrt(sum(rot(:,2)**2));  rot(:,2) = rot(:,2)/vn
    else                                           ! already diagonal: axis-aligned
      if (a >= d) then
        rot(:,1) = [1.d0, 0.d0];  rot(:,2) = [0.d0, 1.d0]
      else
        rot(:,1) = [0.d0, 1.d0];  rot(:,2) = [1.d0, 0.d0]
      endif
    endif

    ! ---- bypass: force the identity rotation so the grid is an axis-aligned box
    ! in (x1,x2). PC1->x1 with spacing d1, PC2->x2 with spacing d2. Diagnostic
    ! switch to confirm the PCA rotation is not the source of interpolation error.
    if (present(use_pca)) then
      if (.not. use_pca) then
        rot(:,1) = [1.d0, 0.d0];  rot(:,2) = [0.d0, 1.d0]
      endif
    endif

    ! ---- spacing (projected resolution) and extent along each PC axis ---------
    dc(1) = sqrt((rot(1,1)*d1)**2 + (rot(2,1)*d2)**2)
    dc(2) = sqrt((rot(1,2)*d1)**2 + (rot(2,2)*d2)**2)

    p1 = rot(1,1)*(x1 - mu(1)) + rot(2,1)*(x2 - mu(2))
    p2 = rot(1,2)*(x1 - mu(1)) + rot(2,2)*(x2 - mu(2))
    pmin1 = minval(p1);  pmax1 = maxval(p1)
    pmin2 = minval(p2);  pmax2 = maxval(p2)

    np    = max(npad, 0)                            ! one node past the data hull each end
    c0(1) = pmin1 - np*dc(1)
    c0(2) = pmin2 - np*dc(2)
    n1    = ceiling((pmax1 - pmin1)/dc(1)) + 1 + 2*np
    n2    = ceiling((pmax2 - pmin2)/dc(2)) + 1 + 2*np
  END SUBROUTINE pca_grid_2d
END MODULE op_math_mod
