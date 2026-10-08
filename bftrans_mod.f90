MODULE bftrans_mod
  USE prec_mod
  USE op_math_mod
  USE level_mod, only: level_type, wcalc, wcalc2
  USE ion_mod
  USE element_mod
  USE molecule_mod
  USE hop_env_mod, only: os, trace
  USE handyfuncs_mod
  implicit none
  private

  type, public :: bftrans_type
    integer :: k, ntt, nnu
    type(level_type) :: level
    class(particle_type), pointer :: ion => null()
    class(particle_type), pointer :: ion2 => null()
    real(dp) :: mint, maxt, minl, maxl
    logical :: ilogl, kvadl, minex, maxex, loglam=.false.
    logical :: ilogt, kvadt, minet, maxet, iteta, pedep, istim
    real(dp), dimension(:), allocatable :: lam, temp, teta, tt, lnlam
    real(4), dimension(:,:), allocatable :: xkap   ! scaled table, see set_table
    real(dp) :: xkap_scale=1
    logical :: log_table=.false.
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
    procedure :: set_table
    procedure :: deactivate
  end type

  logical, public :: H_occupation_probability=.true.

contains
  
  !====================================================================================
  PURE SUBROUTINE deactivate(self)
    class(bftrans_type), intent(inout) :: self
    !----------------------------------------------------------------------------------
    self%active = .false.
  END SUBROUTINE deactivate

  !===============================================================================
  ! Store the absorption table xkap(nnu, ntt) in real(4), in the form in which
  ! interp_xkap interpolates it first: logarithmic for a logarithmic temperature
  ! interpolation, or for a table without temperature dependence that is
  ! interpolated logarithmically in wavelength.  The table is divided by its
  ! maximum, xkap_scale, which keeps cross sections far below the real(4) range
  ! (collision-induced absorption, ~1e-45) representable; interp_xkap multiplies
  ! xkap_scale into the double-precision absorption factor instead.  Logarithms
  ! of non-positive entries are taken of 1e-37.
  !===============================================================================
  SUBROUTINE set_table(self, xkap)
    class(bftrans_type), intent(inout) :: self
    real(dp), dimension(:,:), intent(in) :: xkap
    !---------------------------------------------------------------------------
    self%log_table = merge(self%ilogt, self%ilogl, self%ntt > 1)
    self%xkap_scale = maxval(xkap)
    if (self%xkap_scale <= 0) self%xkap_scale = 1
    if (allocated(self%xkap)) deallocate(self%xkap)
    allocate(self%xkap(size(xkap, 1), size(xkap, 2)))
    if (self%log_table) then
      self%xkap = real(log(merge(xkap, 1d-37, xkap > 0)) - log(self%xkap_scale), 4)
    else
      self%xkap = real(xkap / self%xkap_scale, 4)
    endif
  END SUBROUTINE set_table

  !===============================================================================
  SUBROUTINE load_lam_interp_weights(self, lam, active, fs, fe, iw, quad, lagrange)
    class(bftrans_type), target :: self
    real(dp), dimension(:), intent(in) :: lam
    logical, intent(out) :: active
    integer, intent(out) :: fs, fe
    integer, dimension(:), intent(inout) :: iw
    logical, dimension(:), intent(inout) :: quad
    real(4), dimension(:,:), intent(inout) :: lagrange
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
  SUBROUTINE interp_xkap(self, chi, temp, teta, fakt, lam, fs, fe, iw, quad, lagrange, &
                         ne, nhe, shield, i0, i1)
    class(bftrans_type), target :: self
    class(level_type), pointer :: lvl
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    real(dp), intent(in), dimension(:,:,:), target :: temp, teta
    real(dp), intent(inout), dimension(:,:,:) :: fakt
    real(dp), intent(in), dimension(:) :: lam
    integer, intent(in) :: fs, fe
    integer, dimension(:), intent(in) :: iw
    logical, dimension(:), intent(in) :: quad
    real(4), dimension(:,:), intent(in) :: lagrange
    real(dp), intent(in), dimension(:,:,:), optional :: ne, nhe, shield
    integer, intent(in), optional :: i0, i1
    real(dp), dimension(:,:,:), pointer :: newt
    real(4), dimension(:,:), pointer, contiguous :: tab
    real(4), dimension(:), allocatable :: ab, lam3
    real(4), dimension(self%nnu) :: yy
    real(4) :: w(3), dt, y1, y2, y3, fakt1
    real(dp), dimension(3,self%ntt-1) :: dtdt
    real(dp), dimension(self%ntt-1) :: dtdt1
    real(dp) :: t, diff(3), t3(3)
    integer :: i, f, nx, ny, nz, x, y, z, it, j0, j1, xlo, xhi, mem_offset
    logical :: do_log, do_exp, in_bounds, quadratic
    integer, save :: itimer1=0, itimer2=0, itimer3=0
    !---------------------------------------------------------------------------
    nx = size(temp, 1)
    xlo = 1; xhi = nx
    if (present(i0)) xlo = i0
    if (present(i1)) xhi = i1
    mem_offset = 0
    if (size(chi, 2) == xhi - xlo + 1) then
      mem_offset = xlo - 1                 ! chunk-local buffer: its cell 1 is xlo
    else if (size(chi, 2) /= nx) then
      error stop 'interp_xkap: chi is neither whole-tile nor the given cell range'
    end if
    ny = size(temp, 2)
    nz = size(temp, 3)
    allocate(ab(fs:fe))
    allocate(lam3(fs:fe))
    lam3(fs:fe) = real(lam(fs:fe)**3, 4)
    tab => self%xkap

    if (self%ntt > 1) then
      call trace%begin('interp_xkap_1', itimer=itimer1)
      ! the table is interpolated in temperature in its stored form, then
      ! converted to the form of the wavelength interpolation
      do_exp = self%log_table .and. .not. self%ilogl
      do_log = self%ilogl .and. .not. self%log_table

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

      ! The wavelength step below reads only the table rows j0:j1 of yy, so the
      ! per-cell temperature interpolation and log/exp are confined to them
      j0 = minval(iw(fs:fe), mask=iw(fs:fe) > 0)
      j1 = min(maxval(iw(fs:fe)) + 2, self%nnu)

      do z=1, nz
        do y=1, ny
          do x=xlo, xhi
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
                i = merge(it, it-1, it+2 <= self%ntt)
                diff = t - self%tt(i:i+2)
                w(1) = real(diff(2) * diff(3) * dtdt(1,it), 4)
                w(2) = real(diff(1) * diff(3) * dtdt(2,it), 4)
                w(3) = real(diff(1) * diff(2) * dtdt(3,it), 4)
                do f=j0, j1
                  yy(f) = tab(f,i)*w(1) + tab(f,i+1)*w(2) + tab(f,i+2)*w(3)
                end do
              else
                ! linear interpolation
                dt = real((t - self%tt(it)) * dtdt1(it), 4)
                do f=j0, j1
                  yy(f) = (tab(f,it+1) - tab(f,it))*dt + tab(f,it)
                end do
                ! extrapolation never exceeds the edge of the table
                if (t < self%tt(1)) then
                  yy(j0:j1) = min(yy(j0:j1), tab(j0:j1,1))
                elseif (t >= self%tt(self%ntt)) then
                  yy(j0:j1) = min(yy(j0:j1), tab(j0:j1,self%ntt))
                endif
              endif

              if (do_exp) then
                call exp4_array(yy(j0:j1))
              elseif(do_log) then
                yy(j0:j1) = log(max(yy(j0:j1), tiny(1.0)))
              endif

              i = 0
              do f=fs,fe
                if (iw(f) > 0) then
                  if (i /= iw(f)) then
                    i = iw(f)
                    y1 = yy(i)
                    y2 = yy(i+1)
                    y3 = yy(min(i+2, self%nnu))
                  endif
                  ab(f) = y1*lagrange(1,f) + y2*lagrange(2,f) + y3*lagrange(3,f)
                else
                  ab(f) = merge(-huge(1.0), 0.0, self%ilogl)
                endif
              end do
              if (self%ilogl) call exp4_array(ab)

              if (self%lam3) then
                ab = ab * lam3
              endif

              fakt1 = scaled_factor(fakt(x,y,z))
              chi(fs:fe,x-mem_offset,y,z) = chi(fs:fe,x-mem_offset,y,z) + fakt1*ab
            endif
          end do
        end do
      end do
      call trace%end(itimer1)
    else
      call trace%begin('interp_xkap_2', itimer=itimer2)

      yy = tab(:,1)
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
        ab = ab * lam3
      endif

      if (H_occupation_probability .and. trim(self%key) == 'H BF') then
        lvl => self%ion%levels(self%k)
        call trace%begin('wcalc2', itimer=itimer3)
        call wcalc2(lvl, chi(fs:fe, xlo-mem_offset:xhi-mem_offset, :, :), &
                    lam(fs:fe), ab*real(self%xkap_scale, 4), self%ion%n(xlo:xhi,:,:), &
                    nhe(xlo:xhi,:,:), ne(xlo:xhi,:,:), shield(xlo:xhi,:,:), &
                    wocc=lvl%w(xlo:xhi,:,:), nlvl=lvl%n(xlo:xhi,:,:))
        call trace%end(itimer3)
      else
        do z=1, nz
          do y=1, ny
            do x=xlo, xhi
              fakt1 = scaled_factor(fakt(x,y,z))
              do f=fs,fe
                chi(f,x-mem_offset,y,z) = chi(f,x-mem_offset,y,z) + ab(f)*fakt1
              end do
            end do
          end do
        end do
      endif
      call trace%end(itimer2)
    endif
  contains
    ! The absorption factor times xkap_scale, which the scaled table leaves out,
    ! capped at the real(4) maximum so that the conversion cannot overflow
    pure real(4) function scaled_factor(fk)
      real(dp), intent(in) :: fk
      scaled_factor = real(min(fk*self%xkap_scale, real(huge(1.0), dp)), 4)
    end function scaled_factor
  END SUBROUTINE interp_xkap

  include 'exp_neg.inc'

  !===============================================================================
  SUBROUTINE absorption_factor(self, e_kT, teta, pe, dxi, fakt, i0, i1)
    class(bftrans_type), intent(in) :: self
    class(particle_type), pointer :: ion
    real(dp), dimension(:,:,:), intent(in) :: e_kT, teta, pe, dxi
    real(dp), dimension(:,:,:), intent(inout) :: fakt
    integer, intent(in), optional :: i0, i1        ! optional cell range on dim 1 (threading)
    real(dp), dimension(:,:,:), allocatable :: h1
    integer, dimension(:,:,:), allocatable :: nniv
    real(dp) :: scalar
    integer :: lo, hi
    integer :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('absorption_factor', itimer=itimer)
    ion => self%ion
    lo = 1; hi = size(fakt, 1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    associate (fk   => fakt(lo:hi,:,:), ekt => e_kT(lo:hi,:,:), &
               pev  => pe  (lo:hi,:,:), dxv => dxi (lo:hi,:,:), &
               nion => ion%n(lo:hi,:,:))
    select case(trim(self%key))
     case('H- BF')
      fk = nion
     case('H- FF', 'H2-', 'HE-')
      fk = nion * pev
     case('H BF')
      if (.not. H_occupation_probability) then
        scalar = 2.0898E-26 / (2*self%k**5)
        fk = scalar * ion%levels(self%k)%n(lo:hi,:,:)
      endif
     case('H FF')
      !fakt = 2.0898E-26 / 2 * ion%part * ion%n * exp(-31.30364*teta + 2.3026*dxi*teta) / (31.30364*teta)
      fk = exp(-ekt*(ion%xion - dxv) - ion%part(lo:hi,:,:))
      fk = 2.0898E-26 / (2*RYD*ekt) * fk * nion
     case('H+H', 'H+HE', 'H2+H', 'H2+H2', 'H2+HE')
      fk = nion * self%ion2%n(lo:hi,:,:)
     case('H2+')
      fk = 1d-40 * nion * self%ion2%n(lo:hi,:,:)
     case default
      if (self%pedep) then
        fk = pev * nion
        if (self%peach) then
          fk = (1d-40 / EE) * ekt * fk
        endif
      elseif (self%molecular) then
        fk = nion * ion%awgt / 6.023d23
      elseif (self%regular) then
        fk = self%level%g * exp(-ekt*self%level%xion - ion%part(lo:hi,:,:))
        fk = fk * nion
      endif
    end select
    end associate
    call trace%end(itimer)
  END SUBROUTINE absorption_factor

  !===============================================================================
  SUBROUTINE add_absorption(self, chi, lam, temp, e_kT, teta, pe, dxi, ne, nhe, shield, sorted, i0, i1)
    class(bftrans_type), intent(in), target :: self
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
    real(dp), dimension(:,:,:), intent(in) :: temp, e_kT, teta, pe, dxi
    real(dp), dimension(:,:,:), intent(in), optional :: ne, nhe, shield
    integer, intent(in), optional :: i0, i1        ! optional cell range on dim 1 (threading)
    real(dp), dimension(:), intent(in) :: lam
    real(dp), dimension(:,:,:), allocatable :: fakt
    logical, optional :: sorted
    logical :: active
    real(dp) :: prev, post
    integer :: i,j,k,f,nf, m(2), p(4)=pp
    integer :: fs, fe, lo, hi
    integer, dimension(:), allocatable :: iw
    logical, dimension(:), allocatable :: quad
    real(4), dimension(:,:), allocatable :: lagrange
    integer :: itimer=0
    !---------------------------------------------------------------------------
    if (.not. self%active) return
    call trace%begin('add_absorption', itimer=itimer)
    nf = size(lam)
    allocate(quad(nf))
    allocate(iw(nf), source=-1)
    allocate(lagrange(3,nf), source=0.0)
    lo = 1; hi = size(temp, 1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    allocate(fakt, mold=temp)
    fakt(lo:hi,:,:) = huge(1d0)
    if (.not. allocated(self%ion%n)) then
      error stop 'add_absorption: ion%n not allocated'
    endif
    call self%load_lam_interp_weights(lam, active, fs, fe, iw, quad, lagrange)

    if (active) then
      !if (nf > 1) prev = chi(p(1), p(2), p(3), p(4))

      call self%absorption_factor(e_kT, teta, pe, dxi, fakt, lo, hi)
      call self%interp_xkap(chi, temp, teta, fakt, lam, fs, fe, iw, quad, lagrange, ne, nhe, shield, lo, hi)

      !if (nf > 1) then
      !  post = chi(p(1),p(2),p(3),p(4))
      !  print '(A8, 2x, A12, 2x, A, 8(ES14.6))', self%key, self%label, self%ion%id, &
      !  fakt(p(1),p(2),p(3)), post - prev, post
      !endif
    endif

    deallocate(quad)
    deallocate(iw)
    deallocate(lagrange)
    call trace%end(itimer)
  END SUBROUTINE add_absorption


END MODULE
