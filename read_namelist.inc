! READ_NAMELIST(group) reads a namelist group from io%input.  io and
! namelist_error (namelist_check_mod) must be in scope.  A missing group keeps the defaults.  A group that fails to read is
! written to a scratch unit, so that namelist_error can list its declared keys,
! and the run stops with a report of the unknown keys in the input.
#define READ_NAMELIST(group) block; integer:: nml_iostat, nml_scratch; character(len=512):: nml_iomsg; rewind (io%input); read (io%input, group, iostat=nml_iostat, iomsg=nml_iomsg); if (nml_iostat > 0) then; open (newunit=nml_scratch, status='scratch', form='formatted'); write (nml_scratch, nml=group, delim='quote'); call namelist_error (io%input, nml_scratch, nml_iomsg); end if; end block
