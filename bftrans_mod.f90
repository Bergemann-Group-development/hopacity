MODULE bftrans_mod
  USE prec_mod
  USE op_math_mod
  USE level_mod, only: level_type, wcalc, wcalc2
  USE ion_mod
  USE element_mod
  USE molecule_mod
  USE os_mod
  USE handyfuncs_mod
  USE trace_mod
  implicit none
  private

  type, public :: bftrans_type
    integer :: k, ntt, nnu, verbose
    type(level_type) :: level
    class(particle_type), pointer :: ion => null()
    class(particle_type), pointer :: ion2 => null()
    real(dp) :: mint, maxt, minl, maxl
    logical :: ilogl, kvadl, minex, maxex, loglam=.false.
    logical :: ilogt, kvadt, minet, maxet, iteta, pedep, istim
    real(dp), dimension(:), allocatable :: lam, temp, teta, tt, lnlam
    real(dp), dimension(:,:), allocatable :: xkap
    real(dp), dimension(:,:), allocatable :: log_xkap
    logical, dimension(:), allocatable :: ldisk!, tdisk
    logical :: peach=.false., regular=.false., molecular=.false.
    logical :: lam3=.false.
    logical :: active=.true.
    character(len=5) :: key
    character(len=20) :: label

  contains
    procedure :: interp_xkap
    procedure :: absorption_factor
    procedure :: add_absorption
    procedure :: load_lam_interp_weights
    procedure :: deactivate
  end type

  integer :: fs, fe
  logical, public :: H_occupation_probability=.true.
  integer, dimension(:), allocatable :: iw
  logical, dimension(:), allocatable :: quad
  real(dp), dimension(:,:), allocatable :: lagrange
  real(dp), dimension(:,:), allocatable, target :: pvt_xkap
  !$omp threadprivate(fs, fe, iw, quad, lagrange, pvt_xkap)

contains
  
  !====================================================================================
  PURE SUBROUTINE deactivate(self)
    class(bftrans_type), intent(inout) :: self
    !----------------------------------------------------------------------------------
    self%active = .false.
  END SUBROUTINE deactivate

  !===============================================================================
  SUBROUTINE load_lam_interp_weights(self, lam, active)
    class(bftrans_type), target :: self
    real(dp), dimension(:), intent(in) :: lam
    logical, intent(out) :: active
    real(dp), dimension(:), pointer, contiguous :: xx
    real(dp) :: x
    integer :: nf, f
    integer :: i,j,k
    logical :: quadratic, above, below, skip
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    active = .true.
    nf = size(lam)

    above = minval(lam, 1) > self%maxl
    below = maxval(lam, 1) < self%minl

    skip = (.not. self%maxex .and. above) .or. (.not. self%minex .and. below)
    if (skip) then
      active = .false.
      return
    endif

    call trace%begin('load_lam_interp_weights', itimer=itimer)

    fs = nf
    fe = 1

    do i=1, nf
      x = lam(i)
      above = x > self%maxl
      below = x < self%minl
      skip = (.not. self%maxex .and. above) .or. (.not. self%minex .and. below)
      if (.not. skip) then
        fs = min(i, fs)
        fe = max(i, fe)
        j = -1
        if (below) then
          j = 0
        elseif (above) then
          j = self%nnu
        elseif (i > 1) then
          if (fits_in(iw(i-1))) then
            j = iw(i-1)
          elseif (fits_in(iw(i-1)-1)) then
            j = iw(i-1)-1
          endif
        endif
        if (j<0) then
          do j=1, self%nnu-1
            if (fits_in(j)) exit
          end do
        endif

        quadratic = self%kvadl
        !if (trim(self%key) /= 'H BF') then ! H BF is quadratically extrapolated
        if (j<1 .or. j >= (self%nnu - 1)) then
          quadratic = .false.
        else
          quadratic = quadratic .and. .not. self%ldisk(j)
        endif
        !endif

        ! if there are not enough points, step back
        k = merge(2, 1, quadratic)
        j = min(j, self%nnu-k)
        j = max(j, 1)
        if (self%loglam) then
          x = log(x)
          xx => self%lnlam(j:j+k)
        else
          xx => self%lam(j:j+k)
        endif

        if (quadratic) then
          lagrange(1:1+k,i) = quadratic_weight(xx, x)
        else
          lagrange(1:1+k,i) = linear_weight(xx, x)
        endif

        iw(i) = j
        quad(i) = quadratic
      endif
    end do
    call trace%end(itimer)
  contains
    pure logical function fits_in(idx)
      integer, intent(in) :: idx
      fits_in = idx > 0 .and. idx < self%nnu
      if (fits_in) then
        fits_in = x >= self%lam(idx) .and. x < self%lam(idx+1)
      endif
    end function fits_in
  END SUBROUTINE load_lam_interp_weights

  !===============================================================================
  PURE FUNCTION quadratic_weight(xxx,newx) result(w)
    real(dp), dimension(3), intent(in) :: xxx
    real(dp), intent(in) :: newx
    real(dp) :: diff(3), w(3)
    !---------------------------------------------------------------------------
    diff = newx-xxx
    w(1) = diff(2) * diff(3) / ((xxx(1)-xxx(2)) * (xxx(1)-xxx(3)))
    w(2) = diff(1) * diff(3) / ((xxx(2)-xxx(1)) * (xxx(2)-xxx(3)))
    w(3) = diff(1) * diff(2) / ((xxx(3)-xxx(1)) * (xxx(3)-xxx(2)))
  end FUNCTION quadratic_weight

  !===============================================================================
  PURE FUNCTION linear_weight(xx,newx) result(w)
    real(dp), dimension(2), intent(in) :: xx
    real(dp), intent(in) :: newx
    real(dp) :: w(2), dx
    !---------------------------------------------------------------------------
    dx = 1 / (xx(2) - xx(1))
    w(1) = (xx(2) - newx) * dx
    w(2) = (newx - xx(1)) * dx
  end FUNCTION linear_weight

  !===============================================================================
  SUBROUTINE interp_xkap(self, chi, temp, teta, fakt, lam, ne, nhe, shield)
    class(bftrans_type) :: self
    class(level_type), pointer :: lvl
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
    real(dp), intent(in), dimension(:,:,:), target :: temp, teta
    real(dp), intent(inout), dimension(:,:,:) :: fakt
    real(dp), intent(in), dimension(:) :: lam
    real(dp), intent(in), dimension(:,:,:), optional :: ne, nhe, shield
    real(dp), dimension(:,:,:), pointer :: newt
    real(dp), dimension(:), allocatable :: ab, lam3
    real(dp), dimension(:,:), pointer :: xkap
    real(dp), dimension(3,self%ntt-1) :: dtdt
    real(dp), dimension(self%ntt-1) :: dtdt1
    real(dp) :: ab1, fakt1, t, diff(3), w(3), t3(3), dt, a, b, y1,y2,y3
    real(dp), dimension(self%nnu) :: yy
    real(dp), dimension(self%ntt) :: tt
    integer :: i, f, nx, ny, nz, x, y, z, it, it2, j1, j2, j3
    logical :: skip, do_log, do_exp, in_bounds, quadratic
    integer, save :: itimer1=0, itimer2=0, itimer3=0, itimer4=0, itimer5=0
    !---------------------------------------------------------------------------
    nx = size(temp, 1)
    ny = size(temp, 2)
    nz = size(temp, 3)
    allocate(ab(fs:fe))
    allocate(lam3(fs:fe))
    lam3 = lam**3

    if (self%ntt > 1) then
      call trace%begin('interp_xkap_1', itimer=itimer1)
      do_exp = .false.
      do_log = .false.
      if (self%ilogt) then
        pvt_xkap = self%log_xkap
        do_exp = .not. self%ilogl
      else
        pvt_xkap = self%xkap
        do_log = self%ilogl
      endif

      if (self%iteta) then
        newt => teta
      else
        newt => temp
      endif

      if (self%kvadt) then
        do it=1, self%ntt-1
          if (it+2 <= self%ntt) then
            t3 = self%tt(it:it+2)
          else
            t3 = self%tt(it-1:it+1)
          endif
          dtdt(1,it) = 1. / ((t3(1)-t3(2)) * (t3(1)-t3(3)))
          dtdt(2,it) = 1. / ((t3(2)-t3(1)) * (t3(2)-t3(3)))
          dtdt(3,it) = 1. / ((t3(3)-t3(1)) * (t3(3)-t3(2)))
        end do
      endif

      do it=1, self%ntt-1
        dtdt1(it) = 1. / (self%tt(it+1) - self%tt(it))
      end do

      do z=1, nz
        do y=1, ny
          do x=1, nx
            t = temp(x,y,z)
            in_bounds = .true.
            quadratic = .false.

            if (.not. self%minet) in_bounds = in_bounds .and. (t > self%mint)
            if (.not. self%maxet) in_bounds = in_bounds .and. (t < self%maxt)

            if (in_bounds) then
              t = newt(x,y,z)
              if (t < self%tt(1)) then
                it = 1
              elseif (t >= self%tt(self%ntt)) then
                it = self%ntt - 1
              else
                do it=1, self%ntt-1
                  if (self%tt(it+1) > t) exit
                end do
                quadratic = self%kvadt
              endif

              if (quadratic) then
                ! quadratic interpolation
                if (it+2 <= self%ntt) then
                  i = it
                  xkap => pvt_xkap(:,it:it+2)
                else
                  i = it-1
                  xkap => pvt_xkap(:,it-1:it+1)
                endif
                diff = t - self%tt(i:i+2)
                w(1) = diff(2) * diff(3) * dtdt(1,it)
                w(2) = diff(1) * diff(3) * dtdt(2,it)
                w(3) = diff(1) * diff(2) * dtdt(3,it)
                yy(:) = 0
                yy(:) = yy(:) + xkap(:,1) * w(1)
                yy(:) = yy(:) + xkap(:,2) * w(2)
                yy(:) = yy(:) + xkap(:,3) * w(3)
              else
                ! linear interpolation
                dt = (t - self%tt(it)) * dtdt1(it)
                if (t < self%tt(1)) then
                  do f=1, self%nnu
                    a = pvt_xkap(f,1)
                    b = pvt_xkap(f,2)
                    yy(f) = (b - a)*dt + a
                    yy(f) = min(yy(f), pvt_xkap(f,1))
                  end do
                elseif (t >= self%tt(self%ntt)) then
                  do f=1, self%nnu
                    a = pvt_xkap(f,self%ntt-1)
                    b = pvt_xkap(f,self%ntt)
                    yy(f) = (b - a)*dt + a
                    yy(f) = min(yy(f), pvt_xkap(f,self%ntt))
                  end do
                else
                  do f=1, self%nnu
                    a = pvt_xkap(f,it)
                    b = pvt_xkap(f,it+1)
                    yy(f) = (b - a)*dt + a
                  end do
                endif
              endif

              if (do_exp) then
                yy = exp(yy)
              elseif(do_log) then
                yy = max(yy, tiny(1.d0))
                yy = log(yy)
              endif

              i = 0
              ab = 0
              do f=fs,fe
                if (iw(f) > 0) then
                  if (i /= iw(f)) then
                    i = iw(f)
                    y1 = yy(i)
                    y2 = yy(i+1)
                    y3 = yy(min(i+2, self%nnu))
                  endif
                  if (i /= 0) then
                    ab(f) = sum([y1, y2, y3]*lagrange(:,f))
                  endif
                endif
              end do
              if (self%ilogl) ab = exp(ab)

              if (self%lam3) then
                ab(fs:fe) = ab(fs:fe) * lam3(fs:fe)
              endif

              fakt1 = fakt(x,y,z)
              chi(fs:fe,x,y,z) = chi(fs:fe,x,y,z) + fakt1*ab(fs:fe)
            endif
          end do
        end do
      end do
      call trace%end(itimer1)
    else
      call trace%begin('interp_xkap_2', itimer=itimer2)

      if(self%ilogl) then
        yy = self%log_xkap(:,1)
      else
        yy = self%xkap(:,1)
      endif

      do f=fs,fe
        i = iw(f)
        if (i<=0) then
          ab(f) = 0
        else
          if (quad(f)) then
            ab(f) = sum(yy(i:i+2) * lagrange(1:3,f))
          else
            ab(f) = sum(yy(i:i+1) * lagrange(1:2,f))
          endif
          if (self%ilogl) ab(f) = exp(ab(f))
        endif
      end do

      if (self%lam3) then
        ab(fs:fe) = ab(fs:fe) * lam3(fs:fe)
      endif

      if (H_occupation_probability .and. trim(self%key) == 'H BF') then
        lvl => self%ion%levels(self%k)
        lvl%i = self%k
        call trace%begin('wcalc2', itimer=itimer3)
        call wcalc2(lvl, chi, lam(fs:fe), ab, self%ion%n, nhe, ne, shield)
        call trace%end(itimer3)
      else
        do z=1, nz
          do y=1, ny
            do x=1, nx
              do f=fs,fe
                chi(f,x,y,z) = chi(f,x,y,z) + ab(f)*fakt(x,y,z)
              end do
            end do
          end do
        end do
      endif
      call trace%end(itimer2)
    endif
  END SUBROUTINE interp_xkap

  !===============================================================================
  SUBROUTINE absorption_factor(self, e_kT, teta, pe, dxi, fakt)
    class(bftrans_type), intent(in) :: self
    class(particle_type), pointer :: ion
    real(dp), dimension(:,:,:), intent(in) :: e_kT, teta, pe, dxi
    real(dp), dimension(:,:,:), intent(out) :: fakt
    real(dp), dimension(:,:,:), allocatable :: h1
    integer, dimension(:,:,:), allocatable :: nniv
    real(dp) :: scalar
    integer :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('absorption_factor', itimer=itimer)
    ion => self%ion
    select case(trim(self%key))
     case('H- BF')
      fakt = ion%n
     case('H- FF', 'H2-', 'HE-')
      fakt = ion%n * pe
     case('H BF')
      if (.not. H_occupation_probability) then
        scalar = 2.0898E-26 / (2*self%k**5)
        fakt = scalar * ion%levels(self%k)%n
      endif
     case('H FF')
      !fakt = 2.0898E-26 / 2 * ion%part * ion%n * exp(-31.30364*teta + 2.3026*dxi*teta) / (31.30364*teta)
      fakt = exp(-e_kT*(ion%xion - dxi) - ion%part)
      fakt = 2.0898E-26 / (2*RYD*e_kT) * fakt * ion%n
     case('H+H', 'H+HE', 'H2+H', 'H2+H2', 'H2+HE')
      fakt = ion%n * self%ion2%n
     case('H2+')
      fakt = 1d-40 * ion%n * self%ion2%n
     case default
      if (self%pedep) then
        fakt = pe * ion%n
        if (self%peach) then
          fakt = (1d-40 / EE) * e_kT * fakt
        endif
      elseif (self%molecular) then
        fakt = ion%n * ion%awgt / 6.023d23
      elseif (self%regular) then
        fakt = self%level%g * exp(-e_kT*self%level%xion - ion%part)
        fakt = fakt * ion%n
      endif
    end select
    call trace%end(itimer)
  END SUBROUTINE absorption_factor

  !===============================================================================
  SUBROUTINE add_absorption(self, chi, lam, temp, e_kT, teta, pe, dxi, ne, nhe, shield, sorted)
    class(bftrans_type), intent(in) :: self
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
    real(dp), dimension(:,:,:), intent(in) :: temp, e_kT, teta, pe, dxi
    real(dp), dimension(:,:,:), intent(in), optional :: ne, nhe, shield
    real(dp), dimension(:), intent(in) :: lam
    real(dp), dimension(:,:,:), allocatable :: fakt
    logical, optional :: sorted
    logical :: active
    real(dp) :: prev, post
    integer :: i,j,k,f,nf, m(2), p(4)=pp
    integer :: itimer=0
    !---------------------------------------------------------------------------
    if (.not. self%active) return
    call trace%begin('add_absorption', itimer=itimer)
    nf = size(lam)
    allocate(quad(nf))
    allocate(iw(nf), source=-1)
    allocate(lagrange(3,nf), source=0.d0)
    allocate(pvt_xkap, mold=self%xkap)
    allocate(fakt, mold=temp)
    fakt = huge(1d0)
    if (.not. allocated(self%ion%n)) then
      stop 'add_absorption: ion%n not allocated'
    endif
    call self%load_lam_interp_weights(lam, active)

    if (active) then
      !if (nf > 1) prev = chi(p(1), p(2), p(3), p(4))

      call self%absorption_factor(e_kT, teta, pe, dxi, fakt)
      call self%interp_xkap(chi, temp, teta, fakt, lam, ne, nhe, shield)

      !if (nf > 1) then
      !  post = chi(p(1),p(2),p(3),p(4))
      !  print '(A8, 2x, A12, 2x, A, 8(ES14.6))', self%key, self%label, self%ion%id, &
      !  fakt(p(1),p(2),p(3)), post - prev, post
      !endif
    endif

    deallocate(quad)
    deallocate(iw)
    deallocate(lagrange)
    deallocate(pvt_xkap)
    call trace%end(itimer)
  END SUBROUTINE add_absorption


END MODULE
