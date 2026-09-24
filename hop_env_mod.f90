!===============================================================================
!> Single point of contact between HOP and the DISPATCH framework.
!>
!> In a normal DISPATCH build this module simply re-exports the framework
!> singletons (io, os, trace), so behavior is identical to using the framework
!> modules directly.  When compiled with -DHOP_STANDALONE (see
!> HOP/standalone/Makefile) it instead provides minimal replacements covering
!> exactly the surface HOP uses -- os%check_file, os%mem_register,
!> trace%begin/end, io%run_dir, io%lock, io%namelist_warning, io%warning --
!> which makes HOP compilable as a standalone library with no framework
!> dependencies.
!===============================================================================
MODULE hop_env_mod
#ifndef HOP_STANDALONE
  USE io_mod,    only: io
  USE os_mod,    only: os
  USE trace_mod, only: trace
  implicit none
  public:: io, os, trace
#else
  USE omp_lib
  USE console_mod
  implicit none
  private

  type:: trace_t
  contains
    procedure, nopass:: begin => trace_begin
    procedure, nopass:: end   => trace_end
  end type

  type:: os_t
  contains
    procedure:: check_file
    procedure, nopass:: mem_register
  end type

  type:: lock_t
    integer(kind=omp_lock_kind):: lock
    logical:: initialized=.false.
  contains
    procedure:: set   => lock_set
    procedure:: unset => lock_unset
  end type

  type:: io_t
    integer:: input=-1                      ! namelist input unit (unused standalone)
    integer:: namelist=-1                   ! namelist echo unit (unused standalone)
    character(len=512):: run_dir='./'
    type(lock_t):: lock
  contains
    procedure:: namelist_warning
    procedure, nopass:: warning
  end type

  type(trace_t), public:: trace
  type(os_t),    public:: os
  type(io_t),    public:: io

CONTAINS

  SUBROUTINE trace_begin (id, set_verbose, itimer, detailed_timer)
    character(len=*):: id
    integer, optional:: set_verbose, itimer
    logical, optional:: detailed_timer
  END SUBROUTINE trace_begin

  SUBROUTINE trace_end (itimer, detailed_timer)
    integer, optional:: itimer
    logical, optional:: detailed_timer
  END SUBROUTINE trace_end

  SUBROUTINE check_file (self, name)
    class(os_t):: self
    character(*):: name
    logical:: exists
    !---------------------------------------------------------------------------
    inquire (file=name, exist=exists)
    if (.not. exists) then
      call console_error ('file not found: '//trim(name))
      error stop 'File not found'
    endif
  END SUBROUTINE check_file

  SUBROUTINE mem_register (bits, count, label)
    integer:: bits, count
    character(*), optional:: label
  END SUBROUTINE mem_register

  SUBROUTINE lock_set (self, label, itimer)
    class(lock_t):: self
    character(len=*), optional:: label
    integer, optional:: itimer
    !---------------------------------------------------------------------------
    if (.not. self%initialized) then
      !$omp critical (hop_env_lock_cr)
      if (.not. self%initialized) then
        call omp_init_lock (self%lock)
        self%initialized = .true.
      end if
      !$omp end critical (hop_env_lock_cr)
    end if
    call omp_set_lock (self%lock)
  END SUBROUTINE lock_set

  SUBROUTINE lock_unset (self, label)
    class(lock_t):: self
    character(len=*), optional:: label
    !---------------------------------------------------------------------------
    call omp_unset_lock (self%lock)
  END SUBROUTINE lock_unset

  !> Same contract as io_mod: silent on iostat 0, a missing group (iostat < 0) is
  !> reported only when error is present, a malformed one always.
  SUBROUTINE namelist_warning (self, namelist, iostat, error)
    class(io_t):: self
    character(len=*):: namelist
    integer, optional:: iostat
    logical, optional:: error
    integer:: iostat_l
    !---------------------------------------------------------------------------
    iostat_l = 0
    if (present(iostat)) iostat_l = iostat
    if (iostat_l == 0) return
    if (iostat_l < 0 .and. .not. present(error)) return
    if (iostat_l < 0) then
      call console_warn (ch_hop, 'namelist '//trim(namelist)//' missing')
    else
      call console_warn (ch_hop, 'problem reading namelist '//trim(namelist))
    end if
  END SUBROUTINE namelist_warning

  SUBROUTINE warning (text)
    character(len=*), optional:: text
    !---------------------------------------------------------------------------
    if (present(text)) call console_warn (ch_hop, text)
  END SUBROUTINE warning
#endif
END MODULE hop_env_mod
