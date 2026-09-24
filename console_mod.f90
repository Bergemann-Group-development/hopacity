!===============================================================================
!> Console output for OpenMULTI3D and HOP, with one verbosity level per channel.
!>
!> Levels:  -1 silent (errors only), 0 warnings and run progress (NLTE iteration
!> table, final summaries), 1 setup summaries, 2 per-patch/per-task diagnostics,
!> 3 trace.  Each channel follows the global level unless it has its own level;
!> console_inherit marks a channel without one.
!>
!> Errors go to stderr from every rank.  Warnings and informational output go to
!> console_unit (stdout); informational output comes from the master rank only.
!> Informational lines are written by the caller behind console_on(), which keeps
!> the Fortran format at the call site:
!>
!>   if (console_on(ch_m3d, 2)) write (console_unit,'(a,i0)') 'patch ', id
!>
!> The module has no DISPATCH dependency, so dispatch.x, libhop and hop_prune.x
!> share it; the master flag is set by the DISPATCH side once MPI is known.
!===============================================================================
MODULE console_mod
  USE iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  integer, parameter, public:: ch_atmos=1, ch_m3d=2, ch_hop=3
  integer, parameter, public:: console_inherit=-100
  integer, parameter:: n_channel=3
  integer, public, save:: console_unit=output_unit
  integer, save:: global_level=1, channel_level(n_channel)=console_inherit
  logical, save:: master=.true.
  public console_set_level, console_set_master, console_level, console_on, &
    console_warn, console_error
CONTAINS

!===============================================================================
!> Set the global level (channel absent) or one channel's level; console_inherit
!> returns a channel to the global level.
!===============================================================================
SUBROUTINE console_set_level (level, channel)
  integer, intent(in):: level
  integer, intent(in), optional:: channel
  !-----------------------------------------------------------------------------
  if (present(channel)) then
    channel_level(channel) = level
  else
    global_level = level
  end if
END SUBROUTINE console_set_level

!===============================================================================
SUBROUTINE console_set_master (is_master)
  logical, intent(in):: is_master
  !-----------------------------------------------------------------------------
  master = is_master
END SUBROUTINE console_set_master

!===============================================================================
!> Effective level of a channel
!===============================================================================
INTEGER FUNCTION console_level (channel)
  integer, intent(in):: channel
  !-----------------------------------------------------------------------------
  console_level = channel_level(channel)
  if (console_level == console_inherit) console_level = global_level
END FUNCTION console_level

!===============================================================================
!> True when the channel is at least as verbose as level, on the master rank
!> unless all_ranks is set (output that one task on any rank owns).
!===============================================================================
LOGICAL FUNCTION console_on (channel, level, all_ranks)
  integer, intent(in):: channel, level
  logical, intent(in), optional:: all_ranks
  logical:: any_rank
  !-----------------------------------------------------------------------------
  any_rank = .false.
  if (present(all_ranks)) any_rank = all_ranks
  console_on = (master .or. any_rank) .and. console_level(channel) >= level
END FUNCTION console_on

!===============================================================================
!> A warning from any rank, suppressed only when the channel is silent (-1)
!===============================================================================
SUBROUTINE console_warn (channel, message)
  integer, intent(in):: channel
  character(len=*), intent(in):: message
  !-----------------------------------------------------------------------------
  if (console_level(channel) < 0) return
  !$omp critical (console_cr)
  write (console_unit,'(2a)') 'WARNING: ', trim(message)
  flush (console_unit)
  !$omp end critical (console_cr)
END SUBROUTINE console_warn

!===============================================================================
!> An error from any rank, always written to stderr; the caller decides whether
!> to stop.
!===============================================================================
SUBROUTINE console_error (message)
  character(len=*), intent(in):: message
  !-----------------------------------------------------------------------------
  !$omp critical (console_cr)
  write (error_unit,'(2a)') 'ERROR: ', trim(message)
  flush (error_unit)
  !$omp end critical (console_cr)
END SUBROUTINE console_error

END MODULE console_mod
