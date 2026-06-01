MODULE legacy_math

CONTAINS

  FUNCTION psemi(z,b)
    !
    !  thermal p-FUNCTION for semic empirical collision rates
    !  reference: sobel'man - 'atomic physics ii.  '
    !
    !:
    !: psemi  94-08-29  modifications: (mats carlsson)
    !:        multiple return changed to single return to avoid
    !:        compiler bug in digital fortran
    !:
    implicit real(8) (A-H,O-Z)
    parameter (n1=10)
    dimension  betref(n1),pnref(n1),pcref(n1)
    integer z
    data betref/-2.0,-1.699,-1.398,-1.0,-0.699,-0.398,0.0,0.301,0.602,1.0/
    data pnref/1.160,0.956,0.758,0.493,0.331,0.209,0.100,0.063,0.040,0.023/
    data pcref/1.160,0.977,0.788,0.554,0.403,0.290,0.214,0.201,0.200,0.200/
    !
    !  limit of high temperature: kt >> e
    !
    if (b .lt. 0.01) then
      psemi = 0.27566 * (0.577 + log(b))
      !
      !  intermediate temps (most important for eqm plasmas)
      !  spline interpolation onto logb grid:
      !
    else if(b .gt. 0.01 .and. b .lt. 10.0)then
      bb=log10(b)
      if (z .eq. 1) psemi = splin(bb,betref,pnref,n1,n1)
      if (z .gt. 1) psemi = splin(bb,betref,pcref,n1,n1)
      !
      !  limit of low temperature: kt << e
      !
    else  if(b .gt. 10.0)then
      if(z .eq. 1) psemi = 0.066 / sqrt(b)
      if(z .gt. 1) psemi = 0.200
    endif
    return
  end FUNCTION psemi
  !
  ! ****************************************************************
  !
  FUNCTION splin(t, x, y, n, nnn)
    !
    !  spline interpolation of t in a given table of points(x(i),y(i)).
    !
    ! y(i) is a FUNCTION of x(i) - x(n)>x(n-1)>....x(1) -or- x(1)> x(2)
    ! > ....> x(n) - when x is normalised x(min) = 0, xmax=1,
    ! norm(abs(x(i+1)-x(i))) < 1.0e-8,(i=1,n-1). the array x is divided
    ! into overlapping intervals. each interval consists of nnn points,
    ! the overlap is int(log(nnn)+3) points. ( 3<=nnn<=100,n).
    ! each time after a call when t changes of interval or after a call
    ! with another set of points, a new cubic spline of nnn points has
    ! to be computed, therefore the computing time depends strongly of
    ! the succession of calls.
    ! for example: successive calls, where t is randomly chosen will cost
    ! a lot of computing time - (t often changes of interval) --->
    ! take nnn then as large as possible. successive calls where t changes
    ! slightly cost less computing time, so it will be clear that your
    ! choice of nnn depends heavily upon its way of useage.
    !
    ! computing time:   for the computation of spline coefficients
    !        this routine needs about 21*n multiplications or divisions
    !        for the interpolation by means of the spline coefficients
    !        the routine needs 9 multiplications or divisions.
    !
    ! extrapolation: if t is not in [x(1),x(n)] t is linear extrapolated
    !
    ! a part of the algoritm for the computation of a natural spline
    ! is taken from t.n.e. greville, the linear system is solved by
    ! the method of successive overrelaxation. (error = 1.0e-6)
    !
    ! ref: t.n.e. greville, mathematics research center, u.s. army,
    !        university of wisconsin. mathematical methods.
    !
    ! this routine was implemented by e.b.j. van der zalm, sterrenwacht
    ! utrecht (may 7th 1981). adapted for vax\11 by paul kuin at oxford.
    ! written in real*4 by p.g. judge, oxford
    ! modified so that if the argument is outside the range of x-values
    ! then the first (or last) y-value is taken
    !
    ! input:         t        argument   (real)
    !                x        array of arguments        (real)
    !                y        array of FUNCTION values (real)
    !                n        length of the arrays x,y
    !                        if n = 2 ---> linear interpolation
    !                nnn        number of points of the cubic spline
    !                        if nnn > n   ---> nnn = n
    !                        if nnn > 100 ---> nnn = 100
    !                        if nnn = 2 linear interpolation adopted.
    !
    !                        if nnn is zero nnn=7  (altered by p judge)
    !                       (used to be if nnn omitted )
    ! output:        splin interpolated value at t
    !
    implicit real(8) (A-H,O-Z)
    dimension x(n), y(n), s2(100), s3(100), dely(100), h(100), &
      h2(100), delsqy(100), c(100), b(100), xx(100), yy(100)
    logical search, found, skip
    character*4 init
    data init/'noto'/
    save
    !  nn = number of points of the spline interval
    nn = min(7,n)
    if (nnn .ne. 0) then
      ! pgj        if (%loc(nnn).eq.0) go to 10
      nn = nnn
      nn = min(n,nnn,100)
    endif
    !  if t in [x(k)-error,x(k)+error] t = x(k)
    error = 1.0e-12
    if (t.ne.0) error = abs(1.0e-12*t)
    found = search(x,t,1,n,k,error)

    if (found) then
      splin = y(k)
      return
    endif

    if (n.eq.2.or.nn.eq.2) then
      ! if n=2 or nnn=2 linear interpolation is adopted
      splin = (y(k+1)-y(k))*(t-x(k))/(x(k+1)-x(k)) + y(k)
      return
    endif
    if (n.eq.1) then
      splin = y(1)
      return
    endif
    !       interval of the spline is computed
    iover = log10(float(nn))*3
    int = (k-iover)/(nn-2*iover) + 1
    i1 = (int-1)*(nn-2*iover)
    i2 = i1 + nn -1
    if (i1.lt.1) then
      i2 = nn
      i1 = 1
    endif
    if (i2.gt.n) then
      i2 = n
      i1 = n - nn + 1
    endif
    ishift = i1 - 1
    n1 = nn - 1
    ! initialization check
    if (init.eq.'okay') then
      ! interval check
      skip = .true.
      do i=1,nn
        if (xx(i).ne.x(i+ishift)) skip = .false.
        if (yy(i).ne.y(i+ishift)) skip = .false.
      end do
    else
      init = 'okay'
      skip = .false.
    endif

    if (.not. skip) then
      ! start of the computation of spline coefficients
      rmin = 1.0e37
      rmax = -1.0e37
      do i=1,nn
        if (y(i+ishift).lt.rmin)  rmin = y(i+ishift)
        if (y(i+ishift).gt.rmax)  rmax = y(i+ishift)
      end do
      ! computation of the norm factors xnorm and ynorm
      xnorm = 1/(x(nn+ishift)-x(1+ishift))
      ynorm = rmax - rmin
      if (ynorm.eq.0) ynorm = rmax
      if (ynorm.eq.0) ynorm = 1
      ynorm = 1 / ynorm
      do i=1,n1
        h(i) = (x(i+1+ishift)-x(i+ishift))*xnorm
        dely(i) = ((y(i+1+ishift)-y(i+ishift))/h(i))*ynorm
      end do
      do i=2,n1
        h2(i) = h(i-1) + h(i)
        b(i) = 0.5*h(i-1)/h2(i)
        delsqy(i) = (dely(i)-dely(i-1))/h2(i)
        s2(i) = 2.*delsqy(i)
        c(i) = 3.*delsqy(i)
      end do
      s2(1) = 0
      s2(n1+1) = 0
      ! solution of the linear system of the spline coefficients
      ! by successive overrelaxation.
      ! constants: eps = error criterion in the iterative solution.
      !            omega = relaxation coefficient.
      epsln = 1.0e-8
      omega = 1.0717968
      eta = epsln
      do while (eta-epsln >= 0)
        eta = 0.
        do i=2,n1
          w = (c(i)-b(i)*s2(i-1)-(0.5-b(i))*s2(i+1)-s2(i))*omega
          if (abs(w)-eta > 0) then
            eta = abs(w)
          endif
          s2(i) = s2(i) + w
        end do
      end do
      do i=1,n1
        s3(i) = (s2(i+1)-s2(i))/h(i)
      end do
      ! x and y stored in xx and yy for interval check
      do i=1,nn
        xx(i) = x(i+ishift)
        yy(i) = y(i+ishift)
      end do
    endif
    i = k - i1 + 1
    ht1 = (t-x(i+ishift))*xnorm
    ht2 = (t-x(i+ishift+1) )*xnorm
    if ((t-x(1))*xnorm.le.0) then
      ! extrapolation t < x(1)
      splin = y(1)
      return
    endif
    if ((t-x(n))*xnorm.le.0) then
      ! interpolation by means of the spline coefficients
      ss2 = s2(i) + ht1*s3(i)
      prod = ht1*ht2
      delsqs = (s2(i)+s2(i+1)+ss2)/6.
      splin = y(i+ishift) + (ht1*dely(i)+prod*delsqs)/ynorm
      return
    endif
    ! extrapolation t > x(n)
    splin =  y(n)
    return
  end

!
! **********************************************************************
!
  logical FUNCTION search(arr,x,ia,ib,k,err)
!
! * searches the point x within an error bound -err- in array arr      *
! * if that point is not found, it gives the                           *
! * interval in the array, where that point fits.                      *
! * the interval k if x in [arr(k),arr(k+1)), interval 1 if x in       *
! * (-inf,arr(2)), interval n-1 if x in [arr(n-1),inf)-(arr increasing)*
! * interval 1 if x in (inf,arr(2)),interval n-1 if x in [arr(n-1),    *
! * -inf)-(arr decreasing)                                             *
! * in the worst case a search costs 2*log2(n) cycles.                 *
! * this routine is very economic, when in a sequence of searches,     *
! * the differences between the points which are to be searched are    *
! * small.                                                             *
! * for example: a sequence of searches which goes along the array,    *
! * or a sequence of searches which remains in a small area of         *
! * the array.                                                         *
! *                                                                    *
! * input:     - arr         array with a non decreasing or a non      *
! *                          increasing fucntion.                      *
! *            - ia,ib       lower and upper limit of the array within *
! *                          is searched.                              *
! *            - x           value which you want to search            *
! *            - k           last index which was found                *
! *                          search                                    *
! *            - err         error bound                               *
! *                                                                    *
! * output:    - k           index of the point or interval of the     *
! *                          array which is found.                     *
! *            - search      true - x is found within the errorbound   *
! *                                  -err-                             *
! *                          false - x is not found.                   *
! *                                                                    *
! * send in by e.v.d.zalm, utrecht, sterrewacht, march 24th 1981       *
! **********************************************************************
    implicit real(8) (A-H,O-Z)
    dimension arr(ib)
!
    one=1.0
    search = .true.
    if (arr(ib)-arr(ia).lt.0.0) goto 100
5   if (x.lt.arr(ia).or.x.gt.arr(ib)) goto 25
    if (k.lt.ia.or.k.gt.ib-1) goto 25
    l = k
    ip = sign(one,x-arr(k))
    if (ip.lt.0.) goto 20
10  l = min(k+ip,ib)
    if (arr(l).ge.x) goto 30
    ip = ip * 2
    k = l
    goto 10
20  k = max(l+ip,ia)
    if (arr(k).le.x) goto 30
    ip = ip * 2
    l = k
    goto 20
25  k = ia
    l = ib
30  if ((l-k).le.1) goto 50
    i = int((k+l)/2.)
    if (arr(i).ge.x) goto 40
    k = i
    goto 30
40  l = i
    goto 30
100 if (x.lt.arr(ib).or.x.gt.arr(ia)) goto 250
    if (k.lt.ia.or.k.gt.ib-1) goto 250
    l = k
    ip = sign(one,arr(k)-x)
    if (ip.lt.0.) goto 200
110 l = min(k+ip,ib)
    if (arr(l).le.x) goto 300
    ip = ip * 2
    k = l
    goto 110
200 k = max(l+ip,ia)
    if (arr(k).ge.x) goto 300
    ip = ip * 2
    l = k
    goto 200
250 k = ia
    l = ib
300 if ((l-k).le.1) goto 50
    i = int((k+l)/2.)
    if (arr(i).le.x) goto 400
    k = i
    goto 300
400 l = i
    goto 300
50  if (abs(x-arr(k)).le.err) return
    if (abs(x-arr(l)).gt.err) goto 60
    k = l
    return
60  search = .false.
    return
  end
!
! **********************************************************************
!

END MODULE
