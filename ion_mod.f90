!===============================================================================
!> Atoms, ions and molecules of HOP's equation of state, with their partition functions.
!===============================================================================
MODULE ion_mod
  USE hop_env_mod, only: os
  USE prec_mod
  USE level_mod
  USE line_types_mod, only: linelist_type
  implicit none

  !===============================================================================
  !> A state with an atomic weight awgt [u], for thermal velocities.
  !===============================================================================
  type, extends(state_type) :: massive_type
    real(dp) :: awgt
  contains
    procedure :: inverse_doppler_factor
    procedure :: compute_vbar2
  end type

  !===============================================================================
  !> An isotope: mass number iso and abundance fraction frac.
  !===============================================================================
  type, extends(massive_type) :: isotope_type
    integer  :: iso=0
    real(dp) :: frac=1.0
  end type isotope_type

  !===============================================================================
  !> A species (atom, ion or molecule) with its partition function, levels, line lists
  !> and isotopes.
  !===============================================================================
  type, extends(isotope_type) :: particle_type
    real(dp), dimension(:,:,:), allocatable :: part
    integer :: j=1 !< ionisation stage (1 is neutral)
    integer :: m(3)=[0,0,0], nk
    character(len=6) :: id
    integer :: niso=1
    class(isotope_type), dimension(:), pointer :: isotope => null()
    type(level_type), allocatable :: levels(:)
    !> line-list lines, one structure-of-arrays store per isotope (index 0 = main)
    type(linelist_type), allocatable :: ll(:)
    real(dp), allocatable, dimension(:) :: irwin
    !> lowest temperature the Irwin polynomial is fitted over. Below it the
    !> partition function is continued as a power law anchored on the polynomial
    !> at irwin_tmin and 1.5*irwin_tmin. Irwin's own atomic and diatomic fits
    !> start at 1000 K, but a polyatomic Q curves far too strongly there for a
    !> power law to follow: extrapolating one from 1000/1500 K undershoots the
    !> true Q by a factor of several at 300 K.  A fit that reaches lower
    !> declares it here and keeps the polynomial down to its own limit.
    real(dp) :: irwin_tmin = 1000
    logical :: in_linelist=.false., in_opac=.false.
    logical :: allocated=.false.
  contains
    procedure :: irwin_partf
    procedure :: irwin_poly
  end type

  !===============================================================================
  !> A pointer to a particle, for arrays of particles.
  !===============================================================================
  type particle_pointer_t
    class(particle_type), pointer :: p => null()
  end type

  !===============================================================================
  !> A pointer to an ion, for arrays of ions.
  !===============================================================================
  type ion_pointer_t
    type(ion_type), pointer :: p => null()
  end type ion_pointer_t

  !===============================================================================
  !> An ionisation stage of an element in HOP's equation of state, linked to the
  !> stages above and below it and to the molecules it forms.
  !===============================================================================
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
    logical :: fixed_part = .false.
    !> ion%part depends on more state than temp (H I under mhd_eos carries the
    !> ne-dependent occupation-weighted partition); vetoes temp-only caching
    logical :: state_dependent_part = .false.
  contains
    procedure :: allocate => allocate_particle
    procedure :: deallocate => deallocate_particle
  end type

  logical, public :: debye_shield=.True. !< lower the ionisation energies by Debye shielding

contains

  !===============================================================================
  !> Allocates the arrays f, n, saha and part with shape m, unless they have it
  !> already, and registers their memory.
  !===============================================================================
  SUBROUTINE allocate_particle (self, m)
    class(ion_type) :: self
    integer :: m(3) !< array shape
    integer :: mold(3), s, p
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
  !> Deallocates f and saha, and n and part unless the opacity or the line list still
  !> needs them.
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
  !> ln of the partition function from the Irwin polynomial in ln T.  Below irwin_tmin
  !> it continues as a power law through the values at irwin_tmin and 1.5 irwin_tmin,
  !> above 16000 K linearly in T through the values at 10000 and 16000 K.
  !===============================================================================
  elemental FUNCTION irwin_partf (self, temp) result(ulog)
    class(particle_type), intent(in) :: self
    real(dp), intent(in) :: temp !< temperature [K]
    real(dp) :: ulog1, ulog2, ulog, tt1, tt2
    !---------------------------------------------------------------------------
    if (temp < self%irwin_tmin) then
      tt1=self%irwin_tmin; tt2=1.5*self%irwin_tmin
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
  !> The Irwin polynomial, coefficients in self%irwin from the highest power down.
  !===============================================================================
  pure FUNCTION irwin_poly (self, temp) result(u)
    class(particle_type), intent(in) :: self
    real(dp), intent(in) :: temp !< ln of the temperature [K]
    real(dp) :: u
    integer :: i
    !---------------------------------------------------------------------------
    u = self%irwin(1)
    do i=1, size(self%irwin)-1
      u=u*temp + self%irwin(i+1)
    end do
  END FUNCTION irwin_poly

  !====================================================================================
  !> Inverse relative Doppler width c / sqrt(2kT/m + vmic^2).
  !====================================================================================
  pure elemental function inverse_doppler_factor(self, temp, vmic) result(inv_df)
    class(massive_type), intent(in) :: self
    real(4), intent(in)  :: temp !< temperature [K]
    real(4), intent(in)  :: vmic !< microturbulence [km/s]
    real(4) :: dop, inv_df
    real(4), parameter :: const=2.0d0*BK/UU
    !---------------------------------------------------------------------------------
    dop = const*temp/self%awgt
    if (vmic > 0) dop = dop + (vmic*1e5)**2
    inv_df = CC / sqrt(dop)
  end function inverse_doppler_factor

  !====================================================================================
  !> Mean relative speed squared of the particle and hydrogen, 8kT/(pi mu) [cm^2 s^-2],
  !> with mu the reduced mass, for van der Waals broadening.
  !====================================================================================
  pure elemental function compute_vbar2(self, temp) result(vbar2)
    class(massive_type), intent(in) :: self
    real(4),  intent(in)  :: temp !< temperature [K]
    real(4) :: vbar2
    real(4), parameter :: const=8.0d0*BK/(PI*UU), inv_MH=UU/MH
    !---------------------------------------------------------------------------------
    vbar2 = const*temp * (inv_MH + 1.0/self%awgt)
  end function compute_vbar2

END MODULE ion_mod
