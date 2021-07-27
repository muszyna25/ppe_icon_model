MODULE mo_util_stride

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: util_stride_1d, util_stride_2d, util_get_ptrdiff

  INTERFACE
    SUBROUTINE util_stride_1d(f_out, elemsize, p1, p2) BIND(C)
      USE ISO_C_BINDING, ONLY: C_INT, C_PTR
      INTEGER(C_INT), INTENT(OUT) :: f_out
      INTEGER(C_INT), VALUE, INTENT(IN) :: elemsize
      TYPE(C_PTR), VALUE, INTENT(IN) :: p1, p2
    END SUBROUTINE util_stride_1d

    SUBROUTINE util_stride_2d(f_out, elemsize, p1, p2, p3) BIND(C)
      USE ISO_C_BINDING, ONLY: C_INT, C_PTR
      INTEGER(C_INT), INTENT(OUT) :: f_out(2)
      INTEGER(C_INT), VALUE, INTENT(IN) :: elemsize
      TYPE(C_PTR), VALUE, INTENT(IN) :: p1, p2, p3
    END SUBROUTINE util_stride_2d

    FUNCTION util_get_ptrdiff(a, b) RESULT(s) BIND(C,NAME='util_get_ptrdiff')
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      INTEGER(C_SIZE_T) :: s, a, b
    END FUNCTION util_get_ptrdiff
  END INTERFACE

END MODULE mo_util_stride

