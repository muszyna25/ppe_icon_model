!> @brief profiling package
!!
!! @author Luis Kornblueh, MPIM
!!
!! @date June 2009
!!
!! @details Profileing package for IBM hpct package or NEC ftrace.
!!
!! @b License:
!!
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!!
!! 1. You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! 2. The code may not be re-distributed without the consent of the authors.
!! 3. The copyright notice and statement of authorship must appear in all
!!    copies.
!! 4. You accept the warranty conditions (see WARRANTY).
!! 5. In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!!
!! @b WARRANTY:
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!<
MODULE mo_profile

  IMPLICIT NONE

  PRIVATE

#if defined(__xlC__) && defined(__hpm__)
#include "f_hpm.h"
#endif

  PUBLIC :: trace_init
  PUBLIC :: trace_finalize
  PUBLIC :: trace_start
  PUBLIC :: trace_stop

CONTAINS

  SUBROUTINE trace_init(name, pe)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: pe

#if defined(__xlC__) && defined(__hpm__)
    call f_hpminit(pe, TRIM(name))
#endif

  END SUBROUTINE trace_init

  SUBROUTINE trace_finalize (pe)
    INTEGER, INTENT(in) :: pe

#if defined(__xlC__) && defined(__hpm__)
    call f_hpmterminate(pe)
#endif

  END SUBROUTINE trace_finalize

  SUBROUTINE trace_start(name, number)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: number

#if defined(__SX__) && defined(_FTRACE)
    call ftrace_region_begin(TRIM(name))
#endif
#if defined(__xlC__) && defined(__hpm__)
    CALL f_hpmstart(number, TRIM(name))
#endif

  END SUBROUTINE trace_start

  SUBROUTINE trace_stop(name, number)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: number

#if defined(__SX__) && defined(_FTRACE)
    call ftrace_region_end(TRIM(name))
#endif
#if defined(__xlC__) && defined(__hpm__)
    CALL f_hpmstop(number)
#endif

  END SUBROUTINE trace_stop

END MODULE mo_profile
