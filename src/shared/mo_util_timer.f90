MODULE mo_util_timer
  USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT, C_PTR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: util_cputime
  PUBLIC :: util_walltime
  PUBLIC :: util_gettimeofday
  PUBLIC :: util_init_real_time
  PUBLIC :: util_get_real_time_size
  PUBLIC :: util_read_real_time
  PUBLIC :: util_diff_real_time

  INTERFACE
    FUNCTION util_cputime(user_time, system_time) BIND(C) RESULT(f_result)
      IMPORT C_DOUBLE, C_INT
      REAL(C_DOUBLE), INTENT(OUT) :: user_time, system_time
      INTEGER(C_INT) :: f_result
    END FUNCTION util_cputime

    FUNCTION util_walltime() BIND(C) RESULT(f_result)
      IMPORT C_DOUBLE
      REAL(C_DOUBLE) :: f_result
    END FUNCTION util_walltime

    FUNCTION util_gettimeofday() BIND(C) RESULT(f_result)
      IMPORT C_DOUBLE
      REAL(C_DOUBLE) :: f_result
    END FUNCTION util_gettimeofday

    SUBROUTINE util_init_real_time() BIND(C)
    END SUBROUTINE util_init_real_time

    SUBROUTINE util_get_real_time_size(rt_size) BIND(C)
      IMPORT C_INT
      INTEGER(C_INT), INTENT(OUT) :: rt_size
    END SUBROUTINE util_get_real_time_size

    SUBROUTINE util_read_real_time(it) BIND(C)
      IMPORT C_PTR
      TYPE(C_PTR), VALUE :: it
    END SUBROUTINE util_read_real_time

    SUBROUTINE util_diff_real_time(it1, it2, t) BIND(C)
      IMPORT C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: it1, it2
      REAL(C_DOUBLE), INTENT(OUT) :: t
    END SUBROUTINE util_diff_real_time
  END INTERFACE

END MODULE mo_util_timer
