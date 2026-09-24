!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE hop_env_mod, only: io, os
    USE console_mod
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
        ! Segment metadata for multi-window grids (A2)
        integer :: n_segments=0
        integer, dimension(:), allocatable :: segment_first, segment_last
        real(8), dimension(:), allocatable :: segment_blue_vac, segment_red_vac
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
    character(len=256) :: lam_file='', window_file=''
    ! window_file (A1/A2): Path to file with multi-window specification. One window per
    ! non-comment line: "blue_AA red_AA" in air or vacuum (apply in_air convention).
    ! Requires R > 0 for constant-R grid construction. Windows are validated, sorted,
    ! merged (overlapping/touching), then each generates a constant-R segment. Segment
    ! metadata is stored in spectrum_type%segment_* arrays (A2).
    ! Mutually exclusive with lam_file and give_lam; incompatible with R=0.
    !
    ! A bare aa_red/aa_blue (no daa/nnu/R/lam_file) requests the atom-native output grid
    ! cropped to this vacuum window instead of a fresh uniform grid -- see init / crop.
    logical, public :: crop_output=.false.
    real(8), public :: crop_blue=0, crop_red=0

  contains
    !===============================================================================
    ! Read and validate window_file, build segmented constant-R grid.
    ! Input: window_file path, in_air flag, R (sampling resolution)
    ! Output: lam (concatenated wavelengths in vacuum), segment metadata arrays
    !===============================================================================
    SUBROUTINE build_segmented_grid(wfile, in_air_flag, R_val, lam, nnu_out, self)
      character(len=*), intent(in) :: wfile
      logical, intent(in) :: in_air_flag
      real(8), intent(in) :: R_val
      real(8), dimension(:), allocatable, intent(out) :: lam
      integer, intent(out) :: nnu_out
      class(spectrum_type), intent(inout) :: self

      real(8), dimension(:,:), allocatable :: windows
      integer :: i, j, n_windows, nnu_seg, lunit, iostat, n_kept
      real(8) :: blue, red, fac, blue_vac, red_vac
      real(8), dimension(:), allocatable :: lam_seg, all_lam
      integer, dimension(:), allocatable :: seg_first, seg_last
      real(8), dimension(:), allocatable :: seg_blue, seg_red
      character(len=256) :: line, msg
      logical :: is_comment

      ! Read window file
      call os%check_file(wfile)
      open(newunit=lunit, file=trim(wfile), form='FORMATTED')

      ! First pass: count non-comment lines
      n_windows = 0
      do
        read(lunit, '(A)', iostat=iostat) line
        if (iostat /= 0) exit
        ! Skip empty lines and comments
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        n_windows = n_windows + 1
      end do

      if (n_windows == 0) then
        stop 'window_file is empty (all comments or blank lines)'
      endif

      ! Allocate and read windows
      allocate(windows(n_windows, 2))
      rewind(lunit)
      n_windows = 0
      do
        read(lunit, '(A)', iostat=iostat) line
        if (iostat /= 0) exit
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        n_windows = n_windows + 1
        read(line, *) windows(n_windows, 1), windows(n_windows, 2)
      end do
      close(lunit)

      ! Validate windows
      do i = 1, n_windows
        blue = windows(i, 1)
        red = windows(i, 2)

        if (.not. (blue > 0 .and. red > 0 .and. isfinite(blue) .and. isfinite(red))) then
          stop 'window_file: endpoints must be finite and positive'
        endif

        if (red <= blue) then
          stop 'window_file: red edge must be > blue edge'
        endif
      end do

      ! Sort windows by blue edge
      call sort_windows(windows, n_windows)

      ! Merge overlapping or touching windows
      call merge_overlapping_windows(windows, n_windows)

      if (n_windows == 0) then
        stop 'window_file: all windows were rejected during merging'
      endif

      ! Convert endpoints to vacuum if needed
      if (in_air_flag) then
        windows(:, 1) = air2vac(windows(:, 1))
        windows(:, 2) = air2vac(windows(:, 2))
      endif

      ! Allocate storage for segment metadata
      allocate(seg_first(n_windows))
      allocate(seg_last(n_windows))
      allocate(seg_blue(n_windows))
      allocate(seg_red(n_windows))
      allocate(all_lam(0))

      nnu_out = 0

      ! Build constant-R grid for each window and concatenate
      do i = 1, n_windows
        blue_vac = windows(i, 1)
        red_vac = windows(i, 2)

        ! Build constant-R grid for this window.
        ! Red-edge rule: start at blue edge and multiply by (1 + 1/R) until first sample
        ! at or beyond red edge. This ensures consistent spacing and matches the existing
        ! single-window constant-R behavior (spectrum_mod.f90:115-121).
        fac = 1d0 + 1d0/R_val
        nnu_seg = ceiling(log(red_vac/blue_vac) / log(fac)) + 1

        if (nnu_seg < 2) then
          write(msg,'(a,i0,2(a,f12.4),a,i0,a)') 'window ', i, ': blue=', blue_vac, ' red=', red_vac, &
            ' -> ', nnu_seg, ' points'
          call console_error (msg)
          stop 'window_file: each window must have at least 2 points'
        endif

        allocate(lam_seg(nnu_seg))
        lam_seg(1) = blue_vac
        do j = 2, nnu_seg
          lam_seg(j) = lam_seg(j-1) * fac
        end do

        ! Record segment boundaries
        seg_first(i) = nnu_out + 1
        seg_last(i) = nnu_out + nnu_seg
        seg_blue(i) = blue_vac
        seg_red(i) = red_vac

        ! Append to global wavelength array
        call append_to_array(all_lam, lam_seg)
        nnu_out = nnu_out + nnu_seg

        deallocate(lam_seg)
      end do

      ! Verify grid is strictly ascending (no gaps allowed, but concatenated segments are)
      do i = 2, size(all_lam)
        if (all_lam(i) <= all_lam(i-1)) then
          stop 'build_segmented_grid: grid not strictly ascending after concatenation'
        endif
      end do

      ! Set up spectrum segment metadata
      self%n_segments = n_windows
      allocate(self%segment_first, source=seg_first)
      allocate(self%segment_last, source=seg_last)
      allocate(self%segment_blue_vac, source=seg_blue)
      allocate(self%segment_red_vac, source=seg_red)

      ! Return the concatenated wavelength array
      allocate(lam, source=all_lam)

      deallocate(windows, seg_first, seg_last, seg_blue, seg_red, all_lam)
    END SUBROUTINE build_segmented_grid

    !===============================================================================
    ! Sort windows by blue edge in ascending order (simple bubble sort for small n)
    !===============================================================================
    SUBROUTINE sort_windows(windows, n)
      real(8), intent(inout) :: windows(:, :)
      integer, intent(in) :: n
      integer :: i, j
      real(8) :: tmp_blue, tmp_red

      do i = 1, n-1
        do j = i+1, n
          if (windows(j, 1) < windows(i, 1)) then
            tmp_blue = windows(i, 1)
            tmp_red = windows(i, 2)
            windows(i, 1) = windows(j, 1)
            windows(i, 2) = windows(j, 2)
            windows(j, 1) = tmp_blue
            windows(j, 2) = tmp_red
          endif
        end do
      end do
    END SUBROUTINE sort_windows

    !===============================================================================
    ! Merge overlapping or touching windows. Windows must be sorted by blue edge.
    ! Updates the windows array in place and reduces n.
    !===============================================================================
    SUBROUTINE merge_overlapping_windows(windows, n)
      real(8), intent(inout) :: windows(:, :)
      integer, intent(inout) :: n
      integer :: i, j
      real(8), parameter :: merge_tol = 1d-10

      i = 1
      do while (i < n)
        ! Check if window i and i+1 overlap or touch
        if (windows(i+1, 1) <= windows(i, 2) + merge_tol) then
          ! Merge: extend window i's red edge to the maximum of the two
          windows(i, 2) = max(windows(i, 2), windows(i+1, 2))
          ! Remove window i+1 by shifting
          do j = i+1, n-1
            windows(j, 1) = windows(j+1, 1)
            windows(j, 2) = windows(j+1, 2)
          end do
          n = n - 1
        else
          i = i + 1
        endif
      end do
    END SUBROUTINE merge_overlapping_windows

    !===============================================================================
    ! Append array b to array a, growing a in place
    !===============================================================================
    SUBROUTINE append_to_array(a, b)
      real(8), dimension(:), allocatable, intent(inout) :: a
      real(8), dimension(:), intent(in) :: b
      real(8), dimension(:), allocatable :: tmp
      integer :: na, nb

      na = size(a)
      nb = size(b)

      if (nb == 0) return

      if (na == 0) then
        if (allocated(a)) deallocate(a)
        allocate(a(nb))
        a = b
      else
        allocate(tmp(na + nb))
        tmp(1:na) = a
        tmp(na+1:na+nb) = b
        deallocate(a)
        allocate(a(na+nb))
        a = tmp
        deallocate(tmp)
      endif
    END SUBROUTINE append_to_array

    !===============================================================================
    ! Check if a real number is finite (not inf or nan)
    !===============================================================================
    PURE LOGICAL FUNCTION isfinite(x) result(res)
      real(8), intent(in) :: x
      res = x == x .and. abs(x) < huge(x)
    END FUNCTION isfinite

    !===============================================================================
    SUBROUTINE init(self)
      class(spectrum_type) :: self
      real(8), dimension(:), allocatable :: lam
      real(8) :: fac, blue, red
      integer :: iostat, i, lunit
      logical :: grid_in_vac
      namelist /spectrum_params/ daa, aa_red, aa_blue, nnu, give_lam, in_log, R, &
      lam_file, window_file, in_air
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

      ! No explicit sampling grid (daa/nnu/R/lam_file/window_file). If a bare window is given, request
      ! cropping the atom-native output grid (built later from the model atom's transitions)
      ! to [aa_blue, aa_red]; otherwise the output grid spans the atom's full range.
      if (daa==0 .and. nnu==0 .and. R==0 .and. lam_file=='' .and. window_file=='') then
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

      ! Validate that window_file is not used with incompatible options
      if (window_file /= '') then
        if (R <= 0) then
          stop 'window_file requires R > 0 for constant-R grid construction'
        endif
        if (lam_file /= '') then
          stop 'window_file and lam_file are mutually exclusive'
        endif
        if (give_lam) then
          stop 'window_file and give_lam are mutually exclusive'
        endif
        if (in_log) then
          stop 'window_file and in_log are mutually exclusive; grid is built in linear vacuum wavelength'
        endif
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
      elseif (window_file /= '') then
        call build_segmented_grid(window_file, in_air, R, lam, nnu, self)
        grid_in_vac = .true.
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
      if (self%n_segments > 0) then
        ! Trapezoidal weights over the concatenated grid would span the gaps between
        ! windows and are not quadrature weights for any of them. Each segment carries its
        ! own edge weights instead.
        do i = 1, self%n_segments
          call set_weights(self%wnu(self%segment_first(i):self%segment_last(i)), &
                           self%nu(self%segment_first(i):self%segment_last(i)))
        end do
      else
        call set_weights(self%wnu, self%nu)
      endif
      self%sorted=.true.
      deallocate(lam)
    END SUBROUTINE init

    !===============================================================================
    SUBROUTINE reset(self, lam, nu)
      class(spectrum_type) :: self
      real(8), intent(in), dimension(:), optional :: lam, nu
      real(8), allocatable, dimension(:) :: tmp_lam, tmp_nu
      integer :: i
      logical :: segments_valid
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

      ! Validate segment metadata invariant: the stored metadata is valid exactly when
      ! segment_last(n_segments) equals the installed grid size and the segment ranges
      ! tile 1..nnu without gaps or overlaps.
      if (self%n_segments > 0) then
        segments_valid = .false.
        if (self%nnu > 0 .and. allocated(self%segment_first) .and. allocated(self%segment_last)) then
          ! Check the main invariant: last segment ends at grid boundary
          if (self%segment_last(self%n_segments) == self%nnu) then
            ! First segment must start at 1
            if (self%segment_first(1) == 1) then
              ! Check for gaps or overlaps between segments
              segments_valid = .true.
              do i = 1, self%n_segments - 1
                if (self%segment_last(i) + 1 /= self%segment_first(i+1)) then
                  segments_valid = .false.
                  exit
                endif
              end do
            endif
          endif
        endif

        ! Invalidate if the invariant is violated
        if (.not. segments_valid) then
          if (allocated(self%segment_first)) deallocate(self%segment_first)
          if (allocated(self%segment_last)) deallocate(self%segment_last)
          if (allocated(self%segment_blue_vac)) deallocate(self%segment_blue_vac)
          if (allocated(self%segment_red_vac)) deallocate(self%segment_red_vac)
          call io%warning('spectrum_type%reset: grid size changed; segment metadata invalidated')
          self%n_segments = 0
        endif
      endif
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
        call console_warn (ch_hop, 'integration weight requested for a single frequency')
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
      integer :: lunit, i
      logical :: exists
      !-----------------------------------------------------------------------------
      filename = trim(io%run_dir) // 'save/' // trim(label) // '_lam.bin'
      inquire(file=filename, exist=exists)
      if (.not. exists) then
        open(newunit=lunit, file=filename, access='stream')
        write(lunit) self%lam
        close(lunit)
      endif
      ! A segmented grid is a concatenation of disjoint windows, so a reader cannot
      ! recover the window boundaries from the wavelengths alone without guessing at
      ! gaps. Write them out beside the grid.
      if (self%n_segments > 0) then
        filename = trim(io%run_dir) // 'save/' // trim(label) // '_segments.txt'
        inquire(file=filename, exist=exists)
        if (.not. exists) then
          open(newunit=lunit, file=filename, form='FORMATTED')
          write(lunit,'(a)') '# segment  first  last  blue_vac_AA  red_vac_AA'
          do i = 1, self%n_segments
            write(lunit,'(i9,2i7,2f14.5)') i, self%segment_first(i), self%segment_last(i), &
                                           self%segment_blue_vac(i), self%segment_red_vac(i)
          end do
          close(lunit)
        endif
      endif
    END SUBROUTINE save

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  