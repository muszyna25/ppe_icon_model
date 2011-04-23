MODULE mo_util_backtrace

#ifdef __INTEL_COMPILER
  USE ifcore, ONLY: tracebackqq
#endif

  IMPLICIT NONE

  PRIVATE

#if ! (defined (__SX__) || defined (__INTEL_COMPILER) || defined (__xlC__))
  INTERFACE
    SUBROUTINE util_backtrace() BIND(C,NAME='util_backtrace')
    END SUBROUTINE util_backtrace
  END INTERFACE
#endif

  PUBLIC :: util_backtrace

#if (defined (__SX__) || defined (__INTEL_COMPILER) || defined (__xlC__))
CONTAINS

  SUBROUTINE util_backtrace()
#if defined __INTEL_COMPILER
    CALL tracebackqq
#elif defined __xlC__
    CALL xl__trbk
#elif defined __SX__
    CALL mesput('Traceback: ', 11, 1)
#endif
  END SUBROUTINE util_backtrace
#endif

END MODULE mo_util_backtrace

