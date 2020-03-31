# ACX_FC_ATTRIBUTE_CONTIGUOUS([ACTION-IF-SUCCESS],
#                             [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the Fortran compiler supports the CONTIGUOUS attribute. The
# result is either "yes" or "no".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_attribute_contiguous variable.
#
AC_DEFUN([ACX_FC_ATTRIBUTE_CONTIGUOUS],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_CACHE_CHECK([whether Fortran compiler supports the CONTIGUOUS attribute],
     [acx_cv_fc_attribute_contiguous],
     [AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[[      subroutine conftest_routine(x)
      implicit none
      real, contiguous, target, intent(in) :: x(:, :)
      real, contiguous, pointer :: p(:, :)
      integer :: i, j
      p => x
      do j = 1, size(p, 2)
      do i = 1, size(p, 1)
      p(i, j) = p(i, j) * 0.1
      end do
      end do
      end subroutine]])],
        [acx_cv_fc_attribute_contiguous=yes],
        [acx_cv_fc_attribute_contiguous=no])])
   AS_VAR_IF([acx_cv_fc_attribute_contiguous], [yes], [$1], [m4_default([$2],
     [AC_MSG_FAILURE([Fortran compiler does not support the CONTIGUOUS dnl
attribute])])])])
