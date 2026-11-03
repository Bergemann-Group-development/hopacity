!===============================================================================
!> Chemical elements of HOP's equation of state and their ionisation equilibrium.
!===============================================================================
MODULE element_mod
  USE prec_mod
  USE op_math_mod
  USE ion_mod
  implicit none

  !===============================================================================
  !> An element with its ionisation stages, abundance and isotopes.
  !===============================================================================
  type, extends(isotope_type) :: element_type
    integer :: np, nion, nmolec, niso=1, mcount
    real(dp) :: abund=-100
    type(ion_type), allocatable :: ions(:)
    class(isotope_type), allocatable :: isotope(:)
    integer :: active
    character(len=3) :: id
  contains
    procedure :: compute_partition_functions
    procedure :: compute_ion_equilibrium
    procedure :: print_info
  end type element_type

contains

  !===============================================================================
  !> Sets ion%part, the ln of the partition function, of every ionisation stage: a
  !> fixed value, the Irwin polynomial, or the sum g exp(-xion/kT) over the levels,
  !> with the Fischel & Sparks sum for levels that start a hydrogenic series.
  !===============================================================================
  SUBROUTINE compute_partition_functions(self, temp, e_kT, ne)
    class(element_type), target :: self
    real(dp), dimension(:,:,:), intent(in) :: temp !< temperature [K]
    real(dp), dimension(:,:,:), intent(in) :: e_kT !< 1 eV / kT
    real(dp), dimension(:,:,:), intent(in) :: ne   !< electron number density [cm^-3]
    type(ion_type), pointer :: ion
    type(level_type), pointer :: lvl
    real(dp), dimension(:,:,:), allocatable :: beta
    integer :: i, j, l, u
    !---------------------------------------------------------------------------
    l = lbound(self%ions,1)
    u = ubound(self%ions,1)

    allocate(beta, mold=temp)

    do j=l, u
      ion => self%ions(j)
      call ion%allocate(shape(ne))

      if (ion%fixed_part) then
        ion%part = log(ion%parco)
      elseif (allocated(ion%irwin)) then ! IRWIN FORMAT
        ion%part = ion%irwin_partf(temp)
      else
        ion%part = 0
        do i=1, ubound(ion%levels, 1)
          lvl => ion%levels(i)
          beta = lvl%g * exp(-e_kT*lvl%xion)
          if (lvl%xl > 0) then
            beta = beta * lvl%fischel_sparks(ion%j, e_kT, ne)
          endif
          ion%part = ion%part + beta
        end do
        ion%part = log(ion%part)
      endif
    end do
    deallocate(beta)
  END SUBROUTINE compute_partition_functions

  !===============================================================================
  !> Sets the Saha ratios ion%saha = n(j+1)/n(j) and the fractions ion%f = n(j)/n_H
  !> of the ionisation stages from the partition functions, with the ionisation
  !> energies lowered by dxi times ion%dxi_fac.
  !===============================================================================
  SUBROUTINE compute_ion_equilibrium(self, e_kT, lsaha, dxi)
    class(element_type), target :: self
    real(dp), dimension(:,:,:), intent(in) :: e_kT  !< 1 eV / kT
    real(dp), dimension(:,:,:), intent(in) :: lsaha !< ln(2 (2 pi m_e kT/h^2)^(3/2) / ne)
    real(dp), dimension(:,:,:), intent(in) :: dxi   !< ionisation energy lowering per unit charge [eV]
    real(dp), dimension(:,:,:), allocatable :: fil, partition_ratio, partp
    integer :: j, l, u
    !---------------------------------------------------------------------------
    allocate(fil, partition_ratio, partp, mold=dxi)
    l = lbound(self%ions,1)
    u = ubound(self%ions,1)

    fil = 1
    do j=u-1, l, -1
      associate(ion => self%ions(j))
        !partp = ions(j+1)%part
        !! --- this only shows up in the Oslo legacy code ---!
        !!partp = max(1d-3, partp)
        !! --------------------------------------------------!
        ion%saha = exp(ion%above%part - ion%part + lsaha - e_kT*(ion%xion - dxi*ion%dxi_fac))
        fil = 1 + fil * ion%saha
      end associate
    end do

    self%ions(l)%f = self%abund / fil 

    do j=l+1, u
      self%ions(j)%f = self%ions(j-1)%f * self%ions(j-1)%saha
    end do
    deallocate(fil, partition_ratio, partp)
  END SUBROUTINE compute_ion_equilibrium

  !===============================================================================
  !> Writes the abundance and the fractions, Saha ratios and partition functions of
  !> the ionisation stages in cell p to lunit.
  !===============================================================================
  SUBROUTINE print_info(self, lunit, p)
    class(element_type), target :: self
    integer, intent(in) :: lunit !< unit to write to
    integer, intent(in) :: p(:)  !< cell index
    type(ion_type), pointer :: ion
    character(32) :: info
    character(32) :: f_string
    character(48) :: saha_string
    character(64) :: part_string
    integer :: j
    !---------------------------------------------------------------------------
    write(info,'(A4, e12.4)') self%id, self%abund
    f_string=''
    part_string=''
    saha_string=''
    do j=1, self%nion
      ion => self%ions(j)
      write(f_string, '(A, f8.4)') trim(f_string), ion%f(p(1),p(2),p(3))
    end do

    do j=1, self%nion
      ion => self%ions(j)
      write(saha_string, '(A, es10.2)') trim(saha_string), ion%saha(p(1),p(2),p(3))
    end do

    do j=1, self%nion
      ion => self%ions(j)
      write(part_string, '(A, es10.2)') trim(part_string), 1/ion%part(p(1),p(2),p(3))
    end do
    write(lunit, '(A,A,A,A)'), trim(info), f_string, saha_string, part_string
  END SUBROUTINE print_info

END MODULE element_mod
