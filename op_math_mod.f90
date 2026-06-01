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

  interface pchip_mono
    module procedure pchip_mono4, pchip_mono8
  end interface pchip_mono
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
  FUNCTION interp1d(xx,yy,newx, order, sorted) result(newy)
    real(dp), dimension(:), intent(in) :: xx, yy, newx
    real(dp), dimension(:), allocatable :: newy
    integer, value :: order
    logical, value :: sorted
    integer :: i, j, nx, ny
    procedure(interp1d_linear), pointer :: interp
    !----------------------------------------------------------------------------
    nx = size(xx)
    ny = size(yy)

    if (nx /= ny) then
      print *, 'X array is of length', nx
      print *, 'but Y array is of lenght', ny
      stop
    endif

    interp => interp1d_linear

    !if (.not. present(sorted)) sorted = .false.
    !if (.not. present(order))  order  = 1

    if (order == 2) then
      if (nx < 3) then
        print *, 'too few points for quadratic interpolation'
        stop
      else
        interp => interp1d_quadratic
      endif
    elseif (order > 2) then
      print *, 'only linear and quadratic are supported'
      stop
    endif

    allocate(newy, mold=newx)

    if (sorted) then
      j = 1
      i = 1

      do while (i <= size(newy))
        if (j >= nx) then
          newy(i) = yy(ny)

        elseif (xx(j+1) < newx(i)) then
          j = j+1
          cycle

        elseif (xx(j) < newx(i)) then
          newy(i) = interp(j,xx,yy,newx(i))
        else
          newy(i) = yy(j)
        endif

        i = i+1
      end do
    else
      do i=1, size(newy)
        ! xx and yy are still sorted
        j = searchsorted(xx, newx(i))
        newy(i) = interp(j,xx,yy,newx(i))
      end do
    endif
  end FUNCTION interp1d

  !====================================================================================
  pure FUNCTION interp1d_linear(j,xx,yy,newx) result(newy)
    real(dp), dimension(:), intent(in) :: xx, yy
    real(dp), intent(in) :: newx
    integer, intent(in) :: j
    real(dp) :: m, newy
    !----------------------------------------------------------------------------
    m = (yy(j+1) - yy(j)) / (xx(j+1) - xx(j))
    newy = m * (newx - xx(j)) + yy(j)
  end FUNCTION interp1d_linear

  !====================================================================================
  pure FUNCTION interp1d_quadratic(j,xx,yy,newx) result(newy)
    real(dp), dimension(:), intent(in) :: xx, yy
    real(dp), intent(in) :: newx
    integer, intent(in) :: j
    real(dp) :: newy, w(3)
    integer :: j1, j2, j3
    !-----------------------------------------------------------------------------
    if (size(xx) - j > 1) then
      j1 = j
      j2 = j+1
      j3 = j+2
    else
      j1 = j-1
      j2 = j
      j3 = j+1
    endif

    w(1) = (newx-xx(j2)) * (newx-xx(j3)) / ((xx(j1)-xx(j2)) * (xx(j1)-xx(j3)))
    w(2) = (newx-xx(j1)) * (newx-xx(j3)) / ((xx(j2)-xx(j1)) * (xx(j2)-xx(j3)))
    w(3) = (newx-xx(j1)) * (newx-xx(j2)) / ((xx(j3)-xx(j1)) * (xx(j3)-xx(j2)))

    newy = sum(yy(j1:j3)*w(1:3))
  end FUNCTION interp1d_quadratic

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

  !===============================================================================
  SUBROUTINE interp1d_sorted (xx, yy, newx, newy, fill)
    real(4), dimension(:), intent(in) :: newx
    real(4), dimension(:), intent(in) :: xx, yy
    real(4), dimension(:), intent(out) :: newy
    real(4) :: m, left, right
    real(4), optional :: fill
    integer :: ii, i, nnew, nold
    !-----------------------------------------------------------------------------
    nnew = size(newx, 1)
    nold = size(xx, 1)

    if (present(fill)) then
      left = fill
      right = fill
    else
      left = yy(1)
      right = yy(nold)
    endif

    i = 1
    do ii=1, nnew
      do while (xx(i+1) < newx(ii))
        i = i + 1
        if (i >= nold) then
          i = nold
          exit
        endif
      end do

      if (i == nold) then
        newy(ii:nnew) = right
        exit
      elseif (xx(1) > newx(ii)) then
        newy(ii) = left
      else
        m = (yy(i+1) - yy(i)) /  (xx(i+1) - xx(i))
        newy(ii) = m * (newx(ii) - xx(i)) + yy(i)
      endif
    end do
  END SUBROUTINE interp1d_sorted

  !===============================================================================
  SUBROUTINE linterpolate_3d(f, fout, dx, dy, oldz, lower, ds, newz, pchip)
    real(4), intent(in),  dimension(:,:,:) :: f
    real(4), intent(out), dimension(:,:,:) :: fout
    real(8), intent(in) :: dx, dy, lower(2), ds(2)
    real(4), dimension(:), intent(in) :: oldz, newz
    logical, optional, intent(in) :: pchip
    real(4), dimension(:,:,:), allocatable :: fx
    real(4), dimension(:,:,:), allocatable :: fxy
    real(8) :: df(2), offset(2)
    integer :: i,j,k, nx, ny, nz, ox,oy,oz, zs, ze
    logical :: do_pchip
    !-----------------------------------------------------------------------------
    ox = size(f,1); nx = size(fout,1)
    oy = size(f,2); ny = size(fout,2)
    oz = size(f,3); nz = size(fout,3)
    do_pchip = option(pchip, .true.)
    allocate(fx(nx,oy,oz))
    allocate(fxy(nx,ny,oz))
    df = ds / [dx, dy]
    offset = lower / [dx, dy]

    zs = searchsorted(oldz, newz(1), outbounds_ok=.true.)
    ze = searchsorted(oldz, newz(nz), outbounds_ok=.true.)
    zs = max(1, zs-1)
    ze = min(size(oldz), ze+2)

    ! x-dimension interpolation
    do k = zs, ze
      do j = 1, oy
        if (do_pchip) then
          call pchip_mono4_reg_period(f(:,j,k), offset(1), df(1), fx(:,j,k))
        else
          call interp1d_reg_period(f(:,j,k), offset(1), df(1), fx(:,j,k))
        endif
      end do
      ! y-dimension interpolation
      do i = 1, nx
        if (do_pchip) then
          call pchip_mono4_reg_period(fx(i,:,k), offset(2), df(2), fxy(i,:,k))
        else
          call interp1d_reg_period(fx(i,:,k), offset(2), df(2), fxy(i,:,k))
        endif
      end do
    end do
    ! z-dimension interpolation
    do j = 1, ny
      do i = 1, nx
        if (do_pchip) then
          call pchip_mono4(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:), sorted=.true.)
        else
          call interp1d_sorted(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:))
        endif
      end do
    end do
  END SUBROUTINE linterpolate_3d

  !===============================================================================
  SUBROUTINE linterpolate_3d_rotated(f, fout, dx, dy, oldz, lower, ds, newz, vec, bot2pivot, pchip)
    real(4), intent(in),  dimension(:,:,:) :: f
    real(4), intent(out), dimension(:,:,:) :: fout
    real(4), dimension(:), intent(in) :: oldz, newz
    real(4), dimension(:,:,:), allocatable :: fx, fxy
    real(8), intent(in) :: dx, dy, lower(2), ds(2)
    real(8), intent(in) :: vec(3), bot2pivot
    logical, optional, intent(in) :: pchip
    real(8) :: delta(2), deltak(2), offset(2), df(2)
    integer :: i,j,k, nx, ny, nz, ox,oy,oz, zs, ze
    logical :: do_pchip
    !-----------------------------------------------------------------------------
    ox = size(f,1); nx = size(fout,1)
    oy = size(f,2); ny = size(fout,2)
    oz = size(f,3); nz = size(fout,3)
    do_pchip = option(pchip, .true.)
    allocate(fx(nx,oy,oz))
    allocate(fxy(nx,ny,oz))
    delta = vec(1:2) / (vec(3) * [dx, dy])
    offset = lower / [dx, dy]
    df = ds / [dx, dy]

    zs = searchsorted(oldz, newz(1), outbounds_ok=.true.)
    ze = searchsorted(oldz, newz(nz), outbounds_ok=.true.)
    zs = max(1, zs-1)
    ze = min(size(oldz), ze+2)

    do k = zs, ze
      deltak = (bot2pivot - oldz(k))*delta + offset
      ! x-dimension linear interpolation
      do j = 1, oy
        if (do_pchip) then
          call pchip_mono4_reg_period(f(:,j,k), deltak(1), df(1), fx(:,j,k))
        else
          call interp1d_reg_period(f(:,j,k), deltak(1), df(1), fx(:,j,k))
        endif
      end do
      ! y-dimension linear interpolation
      do i = 1, nx
        if (do_pchip) then
          call pchip_mono4_reg_period(fx(i,:,k), deltak(2), df(2), fxy(i,:,k))
        else
          call interp1d_reg_period(fx(i,:,k), deltak(2), df(2), fxy(i,:,k))
        endif
      end do
    end do
    ! z-dimension interpolation
    do j = 1, ny
      do i = 1, nx
        if (do_pchip) then
          call pchip_mono4(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:), sorted=.true.)
        else
          call interp1d_sorted(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:))
        endif
      end do
    end do
  END SUBROUTINE linterpolate_3d_rotated

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
  ! Monotonic piecewise cubic hermite interpolation
  ! Coded by Richard Hoppe, July 2023
  !===============================================================================
  SUBROUTINE pchip_mono4(xx, yy, newx, newy, sorted, newdydx)
    real(4), intent(in) :: xx(:)            ! X values of the data points
    real(4), intent(in) :: yy(:)            ! Y values of the data points
    real(4), intent(in) :: newx(:)          ! X values to interpolate at
    real(4), intent(out) :: newy(:)         ! Interpolated Y values
    real(4), intent(out), optional :: newdydx(:) ! Derivative at interpolation points
    real(4), save, allocatable :: dydx(:)   ! derivative
    real(4), save, allocatable :: dx(:)
    real(4), save, allocatable :: dy(:)
    logical, intent(in), optional :: sorted
    logical :: is_sorted, get_derivative
    integer :: n                 ! Number of data points
    integer :: n_interp          ! Number of interpolation points
    integer :: j                 ! Loop index for interpolation points
    integer :: i                 ! Loop index for data points
    integer :: imin, imax
    real(4) :: alpha, t1, t2, t3, p
    real(8) :: dydx8, d1, d2
    !$omp threadprivate(dydx, dx, dy)
    !-----------------------------------------------------------------------------
    n = size(xx)
    n_interp = size(newx)
    is_sorted = option(sorted, .false.)
    get_derivative = present(newdydx)

    if (allocated(dydx)) then
      if (size(dydx) /= n) then
        deallocate(dydx)
        deallocate(dx)
        deallocate(dy)
        allocate(dydx(n))
        allocate(dx(n-1))
        allocate(dy(n-1))
      endif
    else
      allocate(dydx(n))
      allocate(dx(n-1))
      allocate(dy(n-1))
    endif

    dydx = 0.
    dx = xx(2:n) - xx(1:n-1)
    dy = yy(2:n) - yy(1:n-1)

    imin = searchsorted(xx, minval(newx), outbounds_ok=.true.)
    imax = searchsorted(xx, maxval(newx), outbounds_ok=.true.)+1
    imin = max(imin, 1)
    imax = min(imax, n)

    if (n == 2) then
      dydx(1) = dy(1)/dx(1)
      dydx(2) = dydx(1)
    endif

    do i=max(imin, 2), min(imax, n-1)
      if (dy(i-1) == 0 .or. dy(i) == 0) cycle

      d1 = dy(i-1)/dx(i-1)
      d2 = dy(i)/dx(i)

      if (sign(1.d0, d1) == sign(1.d0, d2)) then
        ! harmonic mean derivative
        alpha = (1. + dx(i-1)/(dx(i) + dx(i-1))) / 3.0
        dydx8 = d1*d2 / ((1. - alpha)*d1 + alpha*d2)
        dydx(i) = real(dydx8)
      endif
    end do

    ! boundary slopes from quadratic polynomial using 2 points and 1 dydx
    if (imin == 1) then
      dydx(1) = 2.*dy(1)/dx(1) - dydx(2)
    endif

    if (imax == n) then
      dydx(n) = 2.*dy(n-1)/dx(n-1) - dydx(n-1)
    endif

    i = 1
    do j=1, n_interp
      if (is_sorted) then
        ! Find the interval containing newx(j)
        do i=i, n - 1
          if (newx(j) <= xx(i + 1)) then
            exit ! Found the interval
          end if
        end do
        i = min(i, n-1)
      else
        i = searchsorted(xx, newx(j), outbounds_ok=.true.)
      endif

      ! Calculate the interpolation parameter
      t1 = (newx(j) - xx(i)) / dx(i)

      t2 = t1*t1
      t3 = t2*t1
      p = 3.0*t2 - 2.0*t3

      newy(j) = (1.0 - p)*yy(i) + p*yy(i+1) + dx(i)*((t3 - 2.0*t2 + t1)*dydx(i) + (t3-t2)*dydx(i+1))

      if (get_derivative) then
        newdydx(j) =  6.0*(t1 - t2)          * dy(i)/dx(i) &
                   + (3.0*t2 - 4.0*t1 + 1.0) * dydx(i)     &
                   + (3.0*t2 - 2.0*t1)       * dydx(i+1)
      end if
    end do
  END SUBROUTINE pchip_mono4

  !===============================================================================
  ! Monotonic piecewise cubic hermite interpolation
  ! Coded by Richard Hoppe, July 2023
  !===============================================================================
  SUBROUTINE pchip_mono8(xx, yy, newx, newy, sorted)
    real(8), intent(in) :: xx(:)            ! X values of the data points
    real(8), intent(in) :: yy(:)            ! Y values of the data points
    real(8), intent(in) :: newx(:)          ! X values to interpolate at
    real(8), intent(out) :: newy(:)         ! Interpolated Y values
    real(8), save, allocatable :: dydx(:)   ! derivative
    real(8), save, allocatable :: dx(:)
    real(8), save, allocatable :: dy(:)
    logical, intent(in), optional :: sorted
    logical :: is_sorted
    integer :: n                 ! Number of data points
    integer :: n_interp          ! Number of interpolation points
    integer :: j                 ! Loop index for interpolation points
    integer :: i                 ! Loop index for data points
    integer :: imin, imax
    real(8) :: d1, d2, alpha, t1, t2, t3, p
    !$omp threadprivate(dydx, dx, dy)
    !-----------------------------------------------------------------------------
    n = size(xx)
    n_interp = size(newx)
    is_sorted = option(sorted, .false.)

    if (allocated(dydx)) then
      if (size(dydx) /= n) then
        deallocate(dydx)
        deallocate(dx)
        deallocate(dy)
        allocate(dydx(n))
        allocate(dx(n-1))
        allocate(dy(n-1))
      endif
    else
      allocate(dydx(n))
      allocate(dx(n-1))
      allocate(dy(n-1))
    endif

    dx = xx(2:n) - xx(1:n-1)
    dy = yy(2:n) - yy(1:n-1)

    imin = searchsorted(xx, minval(newx), outbounds_ok=.true.)
    imax = searchsorted(xx, maxval(newx), outbounds_ok=.true.)+1
    imin = max(imin, 1)
    imax = min(imax, n)

    if (n == 2) then
      dydx(1) = dy(1)/dx(1)
      dydx(2) = dydx(1)
    endif

    do i=max(imin, 2), min(imax, n-1)
      if (dy(i-1) == 0 .or. dy(i) == 0) cycle

      d1 = dy(i-1)/dx(i-1)
      d2 = dy(i)/dx(i)

      if (sign(1.d0, d1) == sign(1.d0, d2)) then
        ! harmonic mean derivative
        alpha = (1. + dx(i-1)/(dx(i) + dx(i-1))) / 3.0
        dydx(i) = 1. / ((1. - alpha)/d1 + alpha/d2)
      endif
    end do

    ! boundary slopes from quadratic polynomial using 2 points and 1 dydx
    if (imin == 1) then
      dydx(1) = 2.*dy(1)/dx(1) - dydx(2)
    endif

    if (imax == n) then
      dydx(n) = 2.*dy(n-1)/dx(n-1) - dydx(n-1)
    endif

    i = 1
    do j=1, n_interp
      if (is_sorted) then
        ! Find the interval containing newx(j)
        do i=i, n - 1
          if (newx(j) <= xx(i + 1)) then
            exit ! Found the interval
          end if
        end do
        i = min(i, n-1)
      else
        i = searchsorted(xx, newx(j), outbounds_ok=.true.)
      endif

      ! Calculate the interpolation parameter
      t1 = (newx(j) - xx(i)) / dx(i)

      t2 = t1*t1
      t3 = t2*t1
      p = 3.0*t2 - 2.0*t3

      newy(j) = (1.0 - p)*yy(i) + p*yy(i+1) + dx(i)*((t3 - 2.0*t2 + t1)*dydx(i) + (t3-t2)*dydx(i+1))
    end do
  END SUBROUTINE pchip_mono8
  
  !===============================================================================
  SUBROUTINE interp1d_reg_period (yy, new_x0, new_dx, newy)
    real(8), intent(in) :: new_x0    ! Starting x value in grid units
    real(8), intent(in) :: new_dx    ! Spacing of new x values in grid units
    real(4), dimension(:), intent(in) :: yy
    real(4), dimension(:), intent(out) :: newy
    real(4) :: ratio, delta
    integer :: i, i1, i2, nnew, nold
    real(8) :: x0
    !-----------------------------------------------------------------------------
    nnew = size(newy, 1)
    nold = size(yy, 1)
    x0 = new_x0 - new_dx

    do i=1, nnew
      ratio = x0 + i*new_dx
      i1 = floor(ratio)
      delta = ratio - i1
      i1 = modulo(i1, nold) + 1
      i2 = modulo(i1, nold) + 1
      newy(i) = (1.0 - delta) * yy(i1) + delta * yy(i2)
    end do
  END SUBROUTINE interp1d_reg_period

  !===============================================================================
  ! Monotonic piecewise cubic hermite interpolation
  ! Coded by Richard Hoppe, July 2023
  !===============================================================================
  SUBROUTINE pchip_mono4_reg_period(yy, new_x0, new_dx, newy)
    real(4), intent(in) :: yy(:)            ! Y values of the data points
    real(8), intent(in) :: new_x0           ! Starting x value in grid units
    real(8), intent(in) :: new_dx           ! Spacing of new x values in grid units
    real(4), intent(out) :: newy(:)         ! Interpolated Y values
    real(4), save, allocatable :: dydx(:)   ! derivative
    real(4), save, allocatable :: dy(:)
    integer :: n                 ! Number of data points
    integer :: n_interp          ! Number of interpolation points
    integer :: j                 ! Loop index for interpolation points
    integer :: i, i1, i2         ! Loop index for data points
    real(4) :: t1, t2, t3, p, ratio, newx, x0
    real(8) :: dydx8, d1, d2
    !$omp threadprivate(dydx, dy)
    !-----------------------------------------------------------------------------
    n = size(yy)
    n_interp = size(newy)

    if (new_x0 == 0 .and. new_dx - nint(new_dx) < 2.0 * epsilon(new_dx)) then
      ! newx is on the same grid as xx, just copy the values
      x0 = new_x0 - new_dx
      do j=1, n_interp
        newx = x0 + j*new_dx
        i = modulo(nint(newx), n) + 1
        newy(j) = yy(i)
      end do
      return
    endif

    if (allocated(dydx)) then
      if (size(dydx) /= n) then
        deallocate(dydx)
        deallocate(dy)
        allocate(dydx(n))
        allocate(dy(n))
      end if
    else
      allocate(dydx(n))
      allocate(dy(n))
    end if

    dydx = 0.
    dy(1:n-1) = yy(2:n) - yy(1:n-1)
    dy(n) = yy(1) - yy(n)

    do i=1, n
      i1 = modulo(i-2, n) + 1
      i2 = modulo(i-1, n) + 1
      if (dy(i1) == 0 .or. dy(i2) == 0) cycle

      d1 = dy(i1)
      d2 = dy(i2)

      if (sign(1.d0, d1) == sign(1.d0, d2)) then
        ! harmonic mean derivative
        dydx8 = 2*d1*d2 / (d1 + d2)
        dydx(i) = real(dydx8)
      endif
    end do

    x0 = new_x0 - new_dx
    do j=1, n_interp
      ratio = x0 + j*new_dx
      i = floor(ratio)
      i1 = modulo(i,   n) + 1
      i2 = modulo(i+1, n) + 1

      ! Calculate the interpolation parameter
      t1 = ratio - i

      t2 = t1*t1
      t3 = t2*t1
      p = 3.0*t2 - 2.0*t3

      newy(j) = (1.0 - p)*yy(i1) + p*yy(i2) + (t3 - 2.0*t2 + t1)*dydx(i1) + (t3-t2)*dydx(i2)
    end do
  END SUBROUTINE pchip_mono4_reg_period
END MODULE op_math_mod
