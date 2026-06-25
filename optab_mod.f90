MODULE optab_mod
  USE os_mod
  USE prec_mod
  USE const_mod
  USE spectrum_mod, only: spectrum_type
  USE op_math_mod, only: searchsorted, pchip_mono
  implicit none

  type :: optab_type
    real(4), dimension(:), allocatable :: ltemp, lrho
    real(4), dimension(:,:,:), allocatable :: chi, src
    integer :: lunit, ntemp, nrho, nbin
    logical :: initialised=.false.
  contains
    procedure :: init
    procedure :: interpolate
  end type

  type(optab_type), public, save :: optab

contains

  !===============================================================================
  ! initialises opacity table.
  !===============================================================================
  SUBROUTINE init (self, opac_table, sp)
    class(optab_type) :: self
    class(spectrum_type) :: sp
    character(len=*), intent(in) :: opac_table
    integer :: i, lunit, iostat
    !---------------------------------------------------------------------------
    print '(A,A)', 'Reading opacity table: ', trim(opac_table)

    call os%check_file(trim(opac_table)//'.txt')
    open(file=trim(opac_table)//'.txt', newunit=lunit, form='formatted')
    read(lunit, *) self%ntemp, self%nrho, self%nbin
    close(lunit)

    allocate(self%lrho(self%nrho))
    allocate(self%ltemp(self%ntemp))
    allocate(self%chi(self%nbin, self%ntemp, self%nrho))
    allocate(self%src(self%nbin, self%ntemp, self%nrho))

    call os%check_file(trim(opac_table)//'.bin')
    open(file=trim(opac_table)//'.bin', newunit=lunit, form='unformatted', access='stream')
    read(lunit) self%ltemp, self%lrho, self%chi, self%src
    close(lunit)

    sp%nnu = self%nbin
    sp%name = 'bins'
    call sp%reset(lam=[(dble(i), i=1, self%nbin)])
    self%initialised = .true.
  END SUBROUTINE init

  !===============================================================================
  SUBROUTINE interpolate (self, chi, eta, temp, rho)
    class(optab_type) :: self
    real(4), contiguous, dimension(:,:,:,:), target :: chi, eta
    real(4), contiguous, dimension(:,:,:), target :: temp, rho
    real(4), dimension(:), pointer :: chi1, eta1, temp1, rho1, chip, etap
    real(4), dimension(:), allocatable :: ltemp, lrho
    real(4), dimension(:,:), allocatable :: scratch
    integer :: i,j,k,t,r,f,m(4), n
    !---------------------------------------------------------------------------
    temp1(1:size(temp)) => temp(:,:,:)
    rho1 (1:size(rho))  => rho (:,:,:)
    m = shape(chi)
    n = product(m)

    allocate(ltemp, source=log(temp1))
    allocate(lrho , source=log(rho1))
    allocate(scratch(size(ltemp), self%nrho))

    chip(1:n) => chi
    etap(1:n) => eta

    do f=1, size(chi, 4)
      chi1(1:size(temp)) => chip(f:n:m(1))
      eta1(1:size(temp)) => etap(f:n:m(1))
      do i=1, self%nrho
        call pchip_mono(self%ltemp, self%chi(f, :, i), ltemp, scratch(:,i))
      end do
      do i=1, size(temp1)
        call pchip_mono(self%lrho, scratch(i,:), lrho(i:i), chi1(i:i))
      end do
      do i=1, self%nrho
        call pchip_mono(self%ltemp, self%src(f, :, i), ltemp, scratch(:,i))
      end do
      do i=1, size(temp1)
        call pchip_mono(self%lrho, scratch(i,:), lrho(i:i), eta1(i:i))
      end do
      eta1 = eta1 * chi1
    end do
  END SUBROUTINE interpolate
END MODULE optab_mod
