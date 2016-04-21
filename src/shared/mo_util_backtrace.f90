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

#if ! (defined (__INTEL_COMPILER) || defined (__xlC__))
  INTERFACE
    SUBROUTINE util_backtrace() BIND(C,NAME='util_backtrace')
    END SUBROUTINE util_backtrace
  END INTERFACE
#endif

  PUBLIC :: util_backtrace

#if (defined (__INTEL_COMPILER) || defined (__xlC__))
CONTAINS

  SUBROUTINE util_backtrace()
#if defined __INTEL_COMPILER
    CALL tracebackqq
#elif defined __xlC__
    CALL xl__trbk
#endif
  END SUBROUTINE util_backtrace
#endif

END MODULE mo_util_backtrace

