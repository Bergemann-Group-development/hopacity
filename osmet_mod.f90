MODULE osmet_mod
  USE hop_env_mod, only: os, trace
  USE console_mod
  USE prec_mod
  USE const_mod
  USE op_math_mod
  USE interpol_mod, only: pchip_mono
  implicit none

  real(dp) :: eosmet=0.5
  integer :: mnu=154000, mtem=17, mp6=18

  type :: osmet_type
    real(dp), dimension(:), pointer, contiguous :: lam
    real(dp), dimension(:), allocatable :: p6, temp
    real(4), dimension(:,:,:), pointer, contiguous :: ostabl
    logical, dimension(:), pointer, contiguous :: active
    integer :: lunit, np6, ntemp, nnu, maxrec
    logical :: initialised=.false.
  contains
    procedure :: init
    procedure :: add_line_opacities
    procedure :: compute_epsilon
  end type

  type(osmet_type), public :: osmet

contains

  !===============================================================================
  ! initialises osmet type.
  !===============================================================================
  SUBROUTINE init (self, absmet_file)
    class(osmet_type) :: self
    character(len=*), intent(in) :: absmet_file
    integer :: lunit, iostat, maxp6, maxtem
    real(4) :: orig_ostabl(mp6, mtem)
    real(4), allocatable :: p6(:), temp(:)
    real(4), allocatable :: all_ostabl(:,:,:)
    real(8), allocatable :: all_wnum(:)
    real(8) :: wnum
    integer :: nnu
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('osmet_type%init', itimer=itimer)
    if (console_on(ch_hop, 1)) write(console_unit,'(A, A)') 'Reading absmet file: ', trim(absmet_file)
    call os%check_file(absmet_file)
    open(file=absmet_file, newunit=lunit, form='unformatted')
    read(lunit, iostat=iostat) !textatom
    read(lunit, iostat=iostat) !textatom
    read(lunit, iostat=iostat) !ispec
    read(lunit, iostat=iostat) !xiturb
    if (iostat == -1) then
      close(lunit)
      open(file=absmet_file, newunit=lunit, form='unformatted', convert='big_endian')
      read(lunit, iostat=iostat) !textatom
      read(lunit, iostat=iostat) !textatom
      read(lunit, iostat=iostat) !ispec
      read(lunit, iostat=iostat) !xiturb
    endif
    read(lunit) maxp6
    allocate(self%p6(maxp6), p6(maxp6))
    read(lunit) p6
    self%p6 = log10(dble(p6))
    read(lunit) maxtem
    allocate(self%temp(maxtem), temp(maxtem))
    read(lunit) temp
    self%temp = dble(temp)
    read(lunit)
    read(lunit)
    read(lunit)
    read(lunit)
    ! single-pass read into pre-allocated max-size buffers
    allocate(all_ostabl(mp6, mtem, mnu))
    allocate(all_wnum(mnu))
    nnu = 0
    do
      read(lunit, iostat=iostat) wnum, orig_ostabl
      if (iostat < 0 .or. wnum == 0.) exit
      nnu = nnu + 1
      all_wnum(nnu) = wnum
      all_ostabl(:, :, nnu) = orig_ostabl
    end do
    mnu = nnu
    allocate(self%lam(mnu))
    allocate(self%active(mnu), source=.false.)
    call os%mem_register(storage_size(self%ostabl), maxp6*maxtem*mnu, label='osmet')
    allocate(self%ostabl(maxp6, maxtem, mnu))
    ! file is wavenumber-decreasing, i.e. wavelength-increasing: copy in file order
    do nnu = 1, mnu
      self%ostabl(:, :, nnu) = all_ostabl(:maxp6, :maxtem, nnu)
      self%lam(nnu) = CM_TO_AA / dble(all_wnum(nnu))
    end do
    deallocate(all_ostabl, all_wnum)
    close(lunit)
    self%initialised = .true.
    call trace%end(itimer)
  END SUBROUTINE init

  !===============================================================================
  SUBROUTINE add_line_opacities(self, lam, p6, temp, rho, chi)
    class(osmet_type), intent(inout) :: self
    real(8), dimension(:), intent(in) :: lam
    real(8), dimension(:,:,:), intent(in) :: p6, temp, rho
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    real(8) :: p, t
    real(4) :: a1, b1, c1, d1, a2, b2, c2, d2, r
    real(4) :: ab1, cd1, ab2, cd2, abcd1, abcd2
    real(4) :: dx, dy, dz, os
    integer :: ix, iy, iz, i, j, k, f, nx, ny, nz, nf
    !---------------------------------------------------------------------------
    nf = size(chi,1)
    nx = size(chi,2)
    ny = size(chi,3)
    nz = size(chi,4)

    do k=1, nz
      do j=1, ny
        do i=1, nx
          t = temp(i,j,k)
          p = p6(i,j,k)
          r = rho(i,j,k)
          ix = searchsorted(self%p6, p)
          iy = searchsorted(self%temp, t)

          if (ix<0 .or. iy<0) then
            cycle
          endif

          dx = (p - self%p6  (ix)) / (self%p6  (ix+1) - self%p6  (ix))
          dy = (t - self%temp(iy)) / (self%temp(iy+1) - self%temp(iy))

          do f=1, nf
            iz = searchsorted(self%lam, lam(f))
            if (iz < 0) cycle

            dz = (lam(f) - self%lam(iz)) / (self%lam(iz+1) - self%lam(iz))

            a1 = self%ostabl(ix,   iy,   iz)
            b1 = self%ostabl(ix+1, iy,   iz)
            c1 = self%ostabl(ix+1, iy+1, iz)
            d1 = self%ostabl(ix,   iy+1, iz)
            a2 = self%ostabl(ix,   iy,   iz+1)
            b2 = self%ostabl(ix+1, iy,   iz+1)
            c2 = self%ostabl(ix+1, iy+1, iz+1)
            d2 = self%ostabl(ix,   iy+1, iz+1)

            ab1 = (b1 - a1) * dx + a1
            cd1 = (c1 - d1) * dx + d1
            ab2 = (b2 - a2) * dx + a2
            cd2 = (c2 - d2) * dx + d2

            abcd1 = exp(ln10*((cd1 - ab1) * dy + ab1))
            abcd2 = exp(ln10*((cd2 - ab2) * dy + ab2))

            ! last interpolation along frequencies is done on linear opacities
            os = (abcd2 - abcd1) * dz + abcd1
            chi(f,i,j,k) = chi(f,i,j,k) + os*r
          end do
        end do
      end do
    end do
  END SUBROUTINE add_line_opacities

  !===============================================================================
  SUBROUTINE compute_epsilon(self, nu, temp, ne, epsln)
    class(osmet_type), intent(inout) :: self
    real(8), dimension(:), intent(in) :: nu
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(4), dimension(:,:,:,:), intent(inout) :: epsln
    real(8), dimension(:,:,:,:), allocatable, target :: beta, ecnst
    real(8), dimension(:), pointer :: vbeta, vecnst
    real(8), dimension(:), allocatable :: psemi1, psemi2, wnum3
    real(8) :: hkt, eosmet_m1, ecnst1
    integer :: i, j, k, f, nx, ny, nz, nf, n
    !---------------------------------------------------------------------------
    nf = size(epsln,1)
    nx = size(epsln,2)
    ny = size(epsln,3)
    nz = size(epsln,4)

    allocate(wnum3, source=(nu/CC)**3 / 20.6D0)
    allocate(beta(nf, nx, ny, nz))
    allocate(ecnst(nf, nx, ny, nz))
    do k=1, nz
      do j=1, ny
        do i=1, nx
          hkt = h_k / temp(i,j,k)
          beta(:,i,j,k) = hkt * nu
        end do
      end do
    end do
    n = nf*nx*ny*nz
    vbeta(1:n) => beta
    vecnst(1:n) => ecnst
    allocate(psemi1, mold=vbeta)
    allocate(psemi2, mold=vbeta)
    call psemi(1, vbeta, psemi1)
    call psemi(2, vbeta, psemi2)
    eosmet_m1 = 1.0D0 - eosmet
    do i=1, n
      vecnst(i) = eosmet/psemi2(i) + eosmet_m1/psemi1(i)
    end do

    do k=1, nz
      do j=1, ny
        do i=1, nx
          ecnst1 = sqrt(temp(i,j,k)) / ne(i,j,k)
          do f=1, nf
            epsln(f,i,j,k) = 1 / (ecnst1*ecnst(f,i,j,k)*wnum3(f) + 1)
          end do
        end do
      end do
    end do
  END SUBROUTINE compute_epsilon

  !===============================================================================
  ! Scattering fraction 1-eps of line opacity carried by one ionization stage
  ! (z=1 neutral p-function, z>=2 charged): the van Regemorter recipe of
  ! compute_epsilon without the eosmet blend -- the stage of a line-list line
  ! is known exactly, so no neutral/ion mixing fraction is needed.  Needs no
  ! absmet table, hence a plain module procedure.
  !===============================================================================
  SUBROUTINE compute_epsilon_stage(z, nu, temp, ne, one_m_eps)
    integer, intent(in) :: z
    real(8), dimension(:), intent(in) :: nu
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(4), dimension(:,:,:,:), intent(inout) :: one_m_eps
    real(8), dimension(:,:,:,:), allocatable, target :: beta
    real(8), dimension(:), pointer :: vbeta, vpsem
    real(8), dimension(:,:,:,:), allocatable, target :: psem
    real(8), dimension(:), allocatable :: wnum3
    real(8) :: hkt, ecnst1, x
    integer :: i, j, k, f, nx, ny, nz, nf, n
    !---------------------------------------------------------------------------
    nf = size(one_m_eps,1)
    nx = size(one_m_eps,2)
    ny = size(one_m_eps,3)
    nz = size(one_m_eps,4)

    allocate(wnum3, source=(nu/CC)**3 / 20.6D0)
    allocate(beta(nf, nx, ny, nz))
    allocate(psem(nf, nx, ny, nz))
    do k=1, nz
      do j=1, ny
        do i=1, nx
          hkt = h_k / temp(i,j,k)
          beta(:,i,j,k) = hkt * nu
        end do
      end do
    end do
    n = nf*nx*ny*nz
    vbeta(1:n) => beta
    vpsem(1:n) => psem
    call psemi(z, vbeta, vpsem)

    do k=1, nz
      do j=1, ny
        do i=1, nx
          ecnst1 = sqrt(temp(i,j,k)) / ne(i,j,k)
          do f=1, nf
            x = ecnst1 * wnum3(f) / psem(f,i,j,k)
            one_m_eps(f,i,j,k) = x / (1 + x)
          end do
        end do
      end do
    end do
  END SUBROUTINE compute_epsilon_stage

  !===============================================================================
  !  thermal p-FUNCTION for semic empirical collision rates
  !  reference: sobel'man - 'atomic physics ii.  '
  !===============================================================================
  SUBROUTINE psemi(z,b,out)
    integer, intent(in) :: z
    real(8), intent(in), dimension(:) :: b
    real(8), intent(out), dimension(:) :: out
    real(8), dimension(10), parameter :: betref = [-2.0,-1.699,-1.398,-1.0,-0.699,-0.398,0.0,0.301,0.602,1.0]
    real(8), dimension(10), parameter :: pnref  = [1.160,0.956,0.758,0.493,0.331,0.209,0.100,0.063,0.040,0.023]
    real(8), dimension(10), parameter :: pcref  = [1.160,0.977,0.788,0.554,0.403,0.290,0.214,0.201,0.200,0.200]
    real(8), allocatable :: b_packed(:), out_packed(:)
    logical, allocatable :: mask(:)
    !---------------------------------------------------------------------------
    if (z <= 0) then
      out = 0.0
      return
    endif

    allocate(mask, source=(b >= 0.01 .and. b <= 10.0))
    allocate(b_packed, source=pack(b, mask))
    allocate(out_packed, mold=b_packed)

    where (b < 0.01) out = 0.27566 * (0.577 + log(b))

    if (z == 1) then
      where (b > 10) out = 0.066 / sqrt(b)
      if (count(mask) > 0) then
        call pchip_mono(betref, pnref, log10(b_packed), out_packed)
      end if
    else
      where (b > 10) out = 0.200
      if (count(mask) > 0) then
        call pchip_mono(betref, pcref, log10(b_packed), out_packed)
      end if
    endif
    out = unpack(out_packed, mask, out)

    ! if (b .lt. 0.01) then
    !     !  limit of high temperature: kt >> e
    !     psemi = 0.27566 * (0.577 + log(b))
    ! elseif(b .gt. 0.01 .and. b .lt. 10.0)then
    !     !  intermediate temps (most important for eqm plasmas)
    !     !  spline interpolation onto logb grid:
    !     bb=log10(b)
    !     if (z .eq. 1) call pchip_mono(betref,pnref,bb,psemi)
    !     if (z .gt. 1) call pchip_mono(betref,pcref,bb,psemi)
    ! elseif(b .gt. 10.0)then
    !     !  limit of low temperature: kt << e
    !     if(z .eq. 1) psemi = 0.066 / sqrt(b)
    !     if(z .gt. 1) psemi = 0.200
    ! endif
  end SUBROUTINE psemi
END MODULE osmet_mod
