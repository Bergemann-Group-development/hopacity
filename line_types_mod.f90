!=======================================================================================
! line_types_mod -- the SoA bulk-linelist data type, factored into its own leaf module.
!
! linelist_type is a plain structure-of-arrays data container (parallel arrays over the
! lines of one particle/isotope). It lives here, depending only on prec_mod, so that
! particle_type (ion_mod) can embed a linelist_type component WITHOUT depending on the
! heavy linelist_mod (which holds the build/opacity/rate kernels). That inversion lets
! linelist_mod, in turn, USE ion_mod and operate directly on particle_type.
!
! All procedures acting on a linelist_type stay in linelist_mod.
!=======================================================================================
MODULE line_types_mod
  USE prec_mod
  implicit none
  private

  public :: linelist_type

  !-------------------------------------------------------------------------------------
  ! Single precision is sufficient for all line parameters except the log-wavelength
  ! sort key llam0, which must resolve sub-Doppler-width offsets (~1e-5 relative) used
  ! by the march and the profile abscissa; lam0 is dropped and reconstructed via exp(llam0).
  type :: linelist_type
    integer :: nline = 0
    real(dp), allocatable :: llam0(:)    ! log(lam0) (natural log) -- sort key for the march
    real(4),  allocatable :: nu0(:)      ! rest frequency [Hz]
    real(4),  allocatable :: lower_ev(:) ! lower level energy [eV]
    ! oscillator strength + Einstein prefactor folded into one log term, so
    ! Rij = inv_part*idf*exp(gpref - ex), with gpref = loggf*ln10 + log(bij_const/nu0).
    ! The isotope fraction is NOT folded in here; the caller applies it to the
    ! per-isotope population base (bulk path) or via append_line's gpref (model atoms).
    real(4),  allocatable :: gpref(:)
    ! broadening (already scaled by over_4pi/nu0, exactly as the AoS prepare_* routines)
    real(4),  allocatable :: ga(:)       ! radiative damping
    real(4),  allocatable :: gw1(:)      ! van der Waals coefficient
    real(4),  allocatable :: gw2(:)      ! van der Waals exponent
    real(4),  allocatable :: gq(:)       ! Stark constant (sign encodes recipe)
    logical,  allocatable :: unsold(:)   ! .true. -> Unsold vdW, .false. -> Barklem
    ! bookkeeping
    integer,  allocatable :: lev_i(:)    ! NLTE lower level index (filled at opacity time)
    integer,  allocatable :: lev_j(:)    ! NLTE upper level index
    logical,  allocatable :: active(:)   ! line is active
    logical,  allocatable :: used(:)     ! diagnostic: line contributed somewhere
    real(4),  allocatable :: strength(:) ! diagnostic: max peak-opacity/continuum estimate
    ! raw inputs kept only for the linelist_out.txt diagnostic dump
    real(4),  allocatable :: loggf(:), gv_in(:), ga_in(:), gq_in(:)
    ! --- optional extensions, used by model-atom linelists (atom_trans_t%ll) ---------
    ! ipart selects the 4th dim of the caller's lnbase array (population base per ion
    ! stage); unallocated means all lines use lnbase(:,:,:,1) (the bulk-linelist case).
    integer,  allocatable :: ipart(:)
    ! iws/iwe: per-line frequency window (local, 1-based); allocated -> "windowed" mode
    ! (the union grid of an NLTE run is per-line ascending but not globally sorted, so
    ! the marching cursor is replaced by a bisection inside each window).
    integer,  allocatable :: iws(:), iwe(:)
    integer,  allocatable :: icen(:)     ! window centre index (f0 seed; grid is symmetric)
    integer,  allocatable :: itrans(:)   ! target column in gamma_ij/gamma_ji (patch subset)
    real(4),  allocatable :: wA(:)       ! Einstein Aij (x component weight) for rates
    real(4),  allocatable :: cwt(:)      ! component weight (wphi measure; 1 for plain lines)
    real(4),  allocatable :: dej(:)      ! upper_ev - lower_ev = H_E*nu0 [eV] (preconditioning)
    integer,  allocatable :: iso_out(:)  ! isotope number for the linelist_out.txt grouping
  end type linelist_type

END MODULE line_types_mod
