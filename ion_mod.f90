MODULE ion_mod
  USE os_mod
  USE prec_mod
  USE level_mod
  USE bbtrans_mod, only: bbtrans_type
  implicit none

  type, extends(state_type) :: massive_type
    real(dp) :: awgt
  contains
    procedure :: inverse_doppler_factor
    procedure :: compute_vbar2
  end type

  type, extends(massive_type) :: isotope_type
    integer  :: iso=0
    real(dp) :: frac=1.0
  end type isotope_type

  type, extends(isotope_type) :: particle_type
    real(dp), dimension(:,:,:), allocatable :: part
    integer :: j=1 ! ionisation stage (1 is neutral)
    integer :: m(3)=[0,0,0], nk, nline=0
    character(len=6) :: id
    integer :: niso=1
    class(isotope_type), dimension(:), pointer :: isotope => null()
    type(level_type), allocatable :: levels(:)
    type(level_type), allocatable :: nlte_levels(:)
    type(bbtrans_type), dimension(:), allocatable :: line
    real(dp), allocatable, dimension(:) :: irwin
    logical :: in_linelist=.false., in_opac=.false.
    logical :: allocated=.false.
  contains
    procedure :: irwin_partf
    procedure :: irwin_poly
  end type

  type particle_pointer_t
    class(particle_type), pointer :: p => null()
  end type

  type ion_pointer_t
    type(ion_type), pointer :: p => null()
  end type ion_pointer_t

  type, extends(particle_type) :: ion_type
    integer :: nmolec=0, mcount
    type(ion_type), pointer :: above => null()
    type(ion_type), pointer :: below => null()
    class(particle_pointer_t), allocatable :: molecules(:)
    integer :: active, dxi_fac
    real(dp) :: parco
    real(dp), pointer :: abund => null()
    real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
    real(dp), dimension(:,:,:), allocatable :: saha
    real(dp), dimension(:,:,:), allocatable :: fil, ne
    logical :: fixed_part = .false.
  contains
    procedure :: molecules_with
    procedure :: f_of_molecules_with
    procedure :: allocate => allocate_particle
    procedure :: deallocate => deallocate_particle
  end type

  logical, public :: debye_shield=.True.

contains

  !===============================================================================
  FUNCTION f_of_molecules_with(self, other, count_double) result(fsum)
    class(ion_type) :: self
    type(ion_type), optional :: other
    type(particle_type), pointer :: mol
    type(particle_pointer_t), allocatable :: mol_list(:)
    real(dp), dimension(:,:,:), allocatable :: fsum
    logical, optional :: count_double
    integer :: i
    !---------------------------------------------------------------------------
    allocate(fsum(self%m(1),self%m(2),self%m(3)), source=0.d0)

    if (present(other)) then
      call self%molecules_with(other, mol_list, count_double)
    else
      allocate(mol_list, source=self%molecules)
    end if

    do i=1, size(mol_list)
      mol => mol_list(i)%p
      fsum = fsum + mol%f
    end do
    deallocate(mol_list)
  END FUNCTION

  !===============================================================================
  SUBROUTINE molecules_with (self, other, mol_list, count_double)
    class(ion_type), target, intent(in) :: self
    type(ion_type), intent(in) :: other
    logical, intent(in), optional :: count_double
    type(particle_pointer_t), intent(out), allocatable :: mol_list(:)
    class(particle_type), pointer :: mol
    logical :: cd
    integer :: mi, mj, count
    !---------------------------------------------------------------------------
    if (present(count_double)) then
      cd = count_double
    else
      cd = .false.
    endif

    count = 0

    do mi=1, self%nmolec
      mol => self%molecules(mi)%p

      if (mi>1 .and. .not. cd) then
        if(mol%id == self%molecules(mi-1)%p%id) cycle
      endif

      do mj=1, other%nmolec
        if (mol%id == other%molecules(mj)%p%id) count = count + 1
      end do
    end do

    allocate(mol_list(count))

    count = 0
    do mi=1, self%nmolec
      mol => self%molecules(mi)%p

      if (mi>1 .and. .not. cd) then
        if(mol%id == self%molecules(mi-1)%p%id) cycle
      endif

      do mj=1, other%nmolec
        if (mol%id == other%molecules(mj)%p%id) then
          count = count + 1
          mol_list(count)%p => mol
        endif
      end do
    end do
  END SUBROUTINE molecules_with

  !===============================================================================
  SUBROUTINE allocate_particle (self, m)
    class(ion_type) :: self
    integer :: m(3), mold(3), s, p
    !---------------------------------------------------------------------------
    if (all(m == self%m)) return
    self%m = m
    s = storage_size(self%f)
    if (self%allocated) then
      mold = shape(self%f)
      p = product(mold)
      deallocate(self%f, self%saha, self%part, self%n)
      call os%mem_register(s, -p, label='ion%f')
      call os%mem_register(s, -p, label='ion%n')
      call os%mem_register(s, -p, label='ion%saha')
      call os%mem_register(s, -p, label='ion%part')
    endif
    p = product(m)
    allocate(self%f(m(1),m(2),m(3)))
    allocate(self%saha, self%part, self%n, mold=self%f)
    call os%mem_register(s, p, label='ion%f')
    call os%mem_register(s, p, label='ion%n')
    call os%mem_register(s, p, label='ion%saha')
    call os%mem_register(s, p, label='ion%part')
    self%allocated = .true.
  END SUBROUTINE allocate_particle

  !===============================================================================
  SUBROUTINE deallocate_particle (self)
    class(ion_type) :: self
    integer :: m(3), s, p
    !---------------------------------------------------------------------------
    m = self%m
    p = product(m)
    s = storage_size(self%f)
    if (self%allocated) then
      deallocate(self%f, self%saha)
      call os%mem_register(s, -p, label='ion%f')
      call os%mem_register(s, -p, label='ion%saha')
      if (.not. (self%in_opac .or. self%in_linelist)) then
        deallocate(self%n, self%part)
        call os%mem_register(s, -p, label='ion%n')
        call os%mem_register(s, -p, label='ion%part')
      endif
    endif
    self%allocated = .false.
  END SUBROUTINE deallocate_particle

  !===============================================================================
  elemental FUNCTION irwin_partf (self, temp) result(ulog)
    class(particle_type), intent(in) :: self
    real(dp), intent(in) :: temp
    real(dp) :: ulog1, ulog2, ulog, tt1, tt2
    !---------------------------------------------------------------------------
    if (temp < 1000.) then
      tt1=1000; tt2=1500
      ulog1 = self%irwin_poly(log(tt1))
      ulog2 = self%irwin_poly(log(tt2))
      ulog = (ulog2 - ulog1) * log(temp/tt1) / log(tt2/tt1) + ulog1
    elseif (temp <= 16000.) then
      ulog = self%irwin_poly(log(temp))
    else
      tt1=10000; tt2=16000
      ulog1 = self%irwin_poly(log(tt1))
      ulog2 = self%irwin_poly(log(tt2))
      ulog = (ulog2 - ulog1) * (temp-tt2)/(tt2-tt1) + ulog2
    endif
  END FUNCTION irwin_partf

  !===============================================================================
  pure FUNCTION irwin_poly (self, temp) result(u)
    class(particle_type), intent(in) :: self
    real(dp), intent(in) :: temp
    real(dp) :: u
    integer :: i
    !---------------------------------------------------------------------------
    u = self%irwin(1)
    do i=1, size(self%irwin)-1
      u=u*temp + self%irwin(i+1)
    end do
  END FUNCTION irwin_poly

  !====================================================================================
  pure elemental function inverse_doppler_factor(self, temp, vmic) result(inv_df)
    class(massive_type), intent(in) :: self
    real(4), intent(in)  :: temp, vmic
    real(4) :: dop, inv_df
    real(4), parameter :: const=2.0d0*BK/UU
    !---------------------------------------------------------------------------------
    dop = const*temp/self%awgt
    if (vmic > 0) dop = dop + (vmic*1e5)**2
    inv_df = CC / sqrt(dop)
  end function inverse_doppler_factor

  !====================================================================================
  pure elemental function compute_vbar2(self, temp) result(vbar2)
    class(massive_type), intent(in) :: self
    real(4),  intent(in)  :: temp
    real(4) :: vbar2
    real(4), parameter :: const=8.0d0*BK/(PI*UU), inv_MH=UU/MH
    !---------------------------------------------------------------------------------
    vbar2 = const*temp * (inv_MH + 1.0/self%awgt)
  end function compute_vbar2

END MODULE ion_mod
