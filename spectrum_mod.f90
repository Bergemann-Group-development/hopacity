!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE hop_env_mod, only: io, os
    USE const_mod
    USE handyfuncs_mod, only: air2vac, remove_close_entries, sort_asc
    USE op_math_mod, only: get_ivl=>get_interval
    USE trans_mod, only: transition_type
    implicit none

    !-------------------------------------------------------------------------------
    type :: spectrum_type
        ! total number of wavelengths
        integer :: nnu=0
        real(8), dimension(:), allocatable :: nu, wnu, lam, llam
        logical :: sorted=.false.
        ! .true. when the grid is uniform in log(lambda) (e.g. a constant-R grid);
        ! enables the arithmetic (march/bisection-free) opacity kernel add_opacity_equi
        logical :: equidistant=.false.
        character(len=32) :: name=''
    contains
      procedure :: init
      procedure :: reset
      procedure :: get_interval
      procedure :: sort
      procedure :: crop
      procedure :: save
    end type spectrum_type
    !---------------------------------------------------------------------

    public :: set_weights
    real(8), public :: lam_tolerance = 0.001
    logical :: give_lam=.false., in_air=.true.
    logical :: in_log=.false.
    real(8) :: daa=0, aa_red=0, aa_blue=0, R=0
    integer :: nnu=0
    character(len=256) :: lam_file=''
    ! A bare aa_red/aa_blue (no daa/nnu/R/lam_file) requests the atom-native output grid
    ! cropped to this vacuum window instead of a fresh uniform grid -- see init / crop.
    logical, public :: crop_output=.false.
    real(8), public :: crop_blue=0, crop_red=0

  contains
    !===============================================================================
    SUBROUTINE init(self)
      class(spectrum_type) :: self
      real(8), dimension(:), allocatable :: lam
      real(8) :: fac, blue, red
      integer :: iostat, i, lunit
      logical :: grid_in_vac
      namelist /spectrum_params/ daa, aa_red, aa_blue, nnu, give_lam, in_log, R, &
      lam_file, in_air
      namelist /sp_wavelengths/ lam
      !-----------------------------------------------------------------------------
      grid_in_vac = .false.
      rewind (io%input)
      read (io%input, spectrum_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, spectrum_params)
      end if
      call io%namelist_warning ('spectrum_params', iostat)
      write (io%namelist, spectrum_params)

      ! No explicit sampling grid (daa/nnu/R/lam_file). If a bare window is given, request
      ! cropping the atom-native output grid (built later from the model atom's transitions)
      ! to [aa_blue, aa_red]; otherwise the output grid spans the atom's full range.
      if (daa==0 .and. nnu==0 .and. R==0 .and. lam_file=='') then
        crop_output = aa_red /= 0 .and. aa_blue /= 0
        if (crop_output) then
          crop_blue = aa_blue
          crop_red  = aa_red
          if (in_air) then
            crop_blue = air2vac(aa_blue)
            crop_red  = air2vac(aa_red)
          endif
        endif
        return
      endif
      if (give_lam .and. nnu>0) then
        allocate(lam(nnu))
        rewind (io%input)
        read (io%input, sp_wavelengths, iostat=iostat)
        if (iostat > 0) then
          rewind (io%input)
          read (io%input, sp_wavelengths)
        end if
        call io%namelist_warning ('spectrum_params', iostat)
        write (io%namelist, sp_wavelengths)
      elseif (lam_file /= '') then
        nnu = 0
        call os%check_file(lam_file)
        open(newunit=lunit, file=trim(lam_file), form='FORMATTED')
        read(lunit,'(A)', iostat=iostat)
        do while (iostat == 0)
          nnu = nnu + 1
          read(lunit,'(A)', iostat=iostat)
        end do
        rewind(lunit)
        allocate(lam(nnu))
        read(lunit,*) lam
        close(lunit)
      elseif (aa_red /= 0 .and. aa_blue /= 0) then
        if (R /= 0) then
          ! constant-R grid: build it geometric (equidistant in log lambda) directly in
          ! vacuum, so the spacing stays uniform and the arithmetic opacity kernel
          ! (add_opacity_equi) applies. If aa_blue/aa_red are given in air, convert the
          ! endpoints to vacuum up front and skip the global air2vac below.
          blue = aa_blue
          red  = aa_red
          if (in_air) then
            blue = air2vac(aa_blue)
            red  = air2vac(aa_red)
            grid_in_vac = .true.
          endif
          fac = 1d0 + 1d0/R
          nnu = ceiling(log(red/blue) / log(fac)) + 1
          allocate(lam(nnu))
          lam(1) = blue
          do i=2, nnu
            lam(i) = lam(i-1) * fac
          end do
        else
          if (daa /= 0) then
            nnu = nint((aa_red - aa_blue) / daa) + 1
          else
            daa = (aa_red - aa_blue) / real(nnu-1)
          endif
          allocate(lam(nnu))
          lam = [(aa_blue + i*daa, i=0,nnu-1)]
        endif
      else
        stop 'Give either wavelength set or aa_red and aa_blue in spectrum_params'
      endif
      if (in_log) lam = exp(lam)
      if (in_air .and. .not. grid_in_vac) lam = air2vac(lam)
      call self%reset(lam=lam)
      allocate(self%wnu(nnu))
      call set_weights(self%wnu, self%nu)
      self%sorted=.true.
      deallocate(lam)
    END SUBROUTINE init

    !===============================================================================
    SUBROUTINE reset(self, lam, nu)
      class(spectrum_type) :: self
      real(8), intent(in), dimension(:), optional :: lam, nu
      real(8), allocatable, dimension(:) :: tmp_lam, tmp_nu
      !-----------------------------------------------------------------------------
      ! 1. Capture inputs into local temporaries immediately
      if (present(lam)) allocate(tmp_lam, source=lam)
      if (present(nu)) allocate(tmp_nu, source=nu)

      if (allocated(self%lam)) deallocate(self%lam)
      if (allocated(self%nu)) deallocate(self%nu)
      
      if (present(lam)) then
        allocate(self%lam, source=tmp_lam)
        allocate(self%nu, source=HZ_TO_AA / self%lam)
      endif
      if (present(nu)) then
        allocate(self%nu, source=tmp_nu)
        allocate(self%lam, source=HZ_TO_AA / self%nu)
      endif
      self%nnu = size(self%nu)
      self%equidistant = is_log_equidistant(self%lam)
    END SUBROUTINE reset

    !===============================================================================
    ! Test whether lam is uniform in log(lambda), i.e. a constant-R (geometric) grid.
    ! Each consecutive log spacing is compared to the GLOBAL average d0 = (llam(n)-
    ! llam(1))/(n-1) -- the very step the equi kernel uses -- which both bounds the
    ! cumulative index drift and separates geometric grids from arithmetic (daa) ones
    ! (whose local step varies as 1/lambda across the band, far above rtol). The
    ! tolerance is loose enough to admit the tiny (~1e-5) non-uniformity that the
    ! air->vacuum conversion imprints on an R grid. A contiguous subrange of a
    ! geometric grid stays geometric, so sub-spectra inherit the property.
    !===============================================================================
    PURE LOGICAL FUNCTION is_log_equidistant(lam) result(equi)
      real(8), intent(in) :: lam(:)
      integer :: n, i
      real(8) :: d0, d
      real(8), parameter :: rtol = 1d-4
      !-----------------------------------------------------------------------------
      n = size(lam)
      equi = .false.
      if (n < 3) return
      if (lam(1) <= 0d0 .or. lam(n) <= 0d0) return
      d0 = (log(lam(n)) - log(lam(1))) / real(n - 1, 8)
      if (d0 == 0d0) return
      do i = 2, n
        d = log(lam(i)) - log(lam(i-1))
        if (abs(d - d0) > rtol * abs(d0)) return
      end do
      equi = .true.
    END FUNCTION is_log_equidistant

    !===============================================================================
    SUBROUTINE set_weights (weight, nu)
      real(8), dimension(:), intent(out) :: weight
      real(8), dimension(:), intent(in) :: nu
      integer :: j, n
      !-----------------------------------------------------------------------------
      n = size(nu, 1)
      if (n==0) then
        return
      elseif (n==1) then
        print *, 'weight for single frequency wanted!'
        weight(1) = 1.0
        return
      elseif (size(weight) /= n) then
        stop 'weight and nu do not have the same number of elements!'
      endif
  
      weight(1) = 0.5 * ( nu(2) - nu(1) )
      do j=2,n-1
          weight(j) = 0.5 * ( nu(j+1) - nu(j-1) )
      end do
      weight(n) = 0.5 * ( nu(n) - nu(n-1) )
      weight = abs(weight)
    END SUBROUTINE set_weights

    !===============================================================================
    SUBROUTINE get_interval (self, trans, ired, iblue)
      class(spectrum_type), intent(in):: self
      class(transition_type), intent(in) :: trans
      integer, intent(out) :: ired, iblue
      !-----------------------------------------------------------------------------
      if (self%sorted) then
        call get_ivl(self%lam, trans%lam_min, trans%lam_max, iblue, ired, limit=.false.)
        if (iblue > 1) then
          if((trans%lam_min - self%lam(iblue-1)) < lam_tolerance) then
            iblue = iblue - 1
          endif
        endif
        if (ired < self%nnu) then
          if((self%lam(ired+1) - trans%lam_max) < lam_tolerance) then
            ired = ired + 1
          endif
        endif
      else
        ired = trans%ired
        iblue = trans%iblue
      endif
    END SUBROUTINE get_interval

    !===============================================================================
    SUBROUTINE sort (self)
      class(spectrum_type), intent(inout):: self
      !-----------------------------------------------------------------------------
      self%sorted = .true.
      call sort_asc(self%lam)
      ! discard duplicates, nnu contains number of unique elements on exit
      call remove_close_entries(self%lam, lam_tolerance)
      call self%reset(lam=self%lam)
    END SUBROUTINE

    !===============================================================================
    ! Restrict an atom-native output grid to the [crop_blue, crop_red] vacuum window
    ! requested via a bare aa_red/aa_blue (no daa/nnu/R). A no-op unless crop_output is
    ! set, so the internal (bb/bf/all) solve grids -- which need full coverage -- are
    ! left untouched; only the emitted output spectrum is trimmed. Call after the grid
    ! is assembled/sorted but before the per-transition [iblue,ired] windows are set,
    ! since those are recomputed against the (now smaller) grid.
    !===============================================================================
    SUBROUTINE crop (self)
      class(spectrum_type) :: self
      real(8), allocatable :: lam(:)
      logical, allocatable :: keep(:)
      !-----------------------------------------------------------------------------
      if (.not. crop_output) return
      if (.not. allocated(self%lam)) return
      allocate(keep(self%nnu))
      keep = self%lam >= crop_blue .and. self%lam <= crop_red
      if (count(keep) == self%nnu) return          ! nothing outside the window
      lam = pack(self%lam, keep)
      if (size(lam) == 0) then
        call io%warning('spectrum_type%crop: [aa_blue, aa_red] excludes the whole model '// &
                        'atom grid -- output spectrum will be empty')
        return
      endif
      call self%reset(lam=lam)
      ! keep the frequency weights in step with the trimmed grid only if they were built
      ! (the uniform-grid path); the atom-native output leaves self%wnu unallocated.
      if (allocated(self%wnu)) then
        deallocate(self%wnu)
        allocate(self%wnu(self%nnu))
        call set_weights(self%wnu, self%nu)
      endif
      self%sorted = .true.
    END SUBROUTINE crop

    !===============================================================================
    SUBROUTINE save (self, label)
      class(spectrum_type), intent(inout):: self
      character(len=*), intent(in) :: label
      character(len=256) :: filename
      integer :: lunit
      logical :: exists
      !-----------------------------------------------------------------------------
      filename = trim(io%run_dir) // 'save/' // trim(label) // '_lam.bin'
      inquire(file=filename, exist=exists)
      if (.not. exists) then
        open(newunit=lunit, file=filename, access='stream')
        write(lunit) self%lam
        close(lunit)
      endif
    END SUBROUTINE save

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  