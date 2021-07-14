!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_backtrace

#ifdef __INTEL_COMPILER
  USE ifcore, ONLY: tracebackqq
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: util_backtrace

#ifdef __INTEL_COMPILER

CONTAINS
  SUBROUTINE util_backtrace()
    CALL tracebackqq
  END SUBROUTINE util_backtrace

#elif defined __xlC__

CONTAINS
  SUBROUTINE util_backtrace()
    CALL xl__trbk
  END SUBROUTINE util_backtrace

#else

  INTERFACE
    SUBROUTINE util_backtrace() BIND(C)
    END SUBROUTINE util_backtrace
  END INTERFACE

#endif

END MODULE mo_util_backtrace

