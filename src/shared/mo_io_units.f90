MODULE mo_io_units

  IMPLICIT NONE
  
  PUBLIC
  
  ! This paramter is taken from /usr/include/stdio.h (ANSI C standard). If 
  ! problems with filename length appear, check the before mentioned file.
  
  INTEGER, PARAMETER :: filename_max = 1024
  
  ! Standard I/O-units
  
#ifdef hpux
  INTEGER, PARAMETER :: nerr  = 7     ! error output
#else
  INTEGER, PARAMETER :: nerr  = 0     ! error output
#endif
  INTEGER, PARAMETER :: nlog  = 1     ! standard log file unit
  INTEGER, PARAMETER :: nin   = 5     ! standard input
  INTEGER, PARAMETER :: nout  = 6     ! standard output  
  
  INTEGER, PARAMETER, PRIVATE :: none = -1  ! unit given back, when nothing 
                                            ! in the allowed range is available
  
END MODULE mo_io_units











