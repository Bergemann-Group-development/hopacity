!===============================================================================
!> atom_convert.x atom.xx out.txt -- convert a MULTI-format model atom to the
!> table format dispatch.x reads, with the conversion a run applies to its own
!> MULTI atoms (HOP/atom_convert_mod.f90).  An unreadable atom ends the program
!> with a nonzero status and a diagnostic on stderr.
!===============================================================================
PROGRAM atom_convert
  USE iso_fortran_env, only: error_unit
  USE atom_convert_mod, only: convert_classic_atom
  implicit none
  character(len=4096) :: in_file, out_file
  !-----------------------------------------------------------------------------
  if (command_argument_count() /= 2) then
    write(error_unit, '(a)') 'usage: atom_convert.x atom.xx out.txt'
    error stop 2
  endif
  call get_command_argument(1, in_file)
  call get_command_argument(2, out_file)
  call convert_classic_atom(trim(in_file), trim(out_file))
END PROGRAM atom_convert
