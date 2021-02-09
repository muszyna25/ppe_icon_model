dnl
dnl Copyright  (C)  2019  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl ACX_FORTRAN_CHECK_INTEGER_PARAMETER([PARAMETER_VARIABLE],
dnl         [LOW-BOUND],[HIGH-BOUND],[INCREMENT],[EXTRA-DECL],[SHELL_VARIABLE],
dnl         [subst?],[define?])
dnl
dnl sets shell variable acx_cv_fortran_mangled_parameter
dnl (and SHELL_VARIABLE if specified) to value of PARAMETER_VARIABLE
dnl Fortran Parameter
dnl
AC_DEFUN([ACX_FORTRAN_CHECK_INTEGER_PARAMETER],
  [AS_VAR_PUSHDEF([acx_fortran_Param], [acx_cv_fortran_$1])dnl
   AC_CACHE_CHECK([value of $1],[acx_fortran_Param],
     [AC_LANG_PUSH([Fortran])
      ACX_FORTRAN_CHECK_INTEGER_EXPRESSION(,[$1],[$5],[$2],[$3],[$4],[acx_fortran_Param])
      AC_LANG_POP([Fortran])])
   m4_ifval([$6],
     [AS_VAR_COPY([$6],[acx_fortran_Param])
      m4_ifval([$7],[AC_SUBST([$6])])
      m4_ifval([$8],
        [AC_DEFINE_UNQUOTED([$6],AS_VAR_GET([acx_fortran_Param]),[$8])])])
   AS_VAR_POPDEF([acx_fortran_Param])])
dnl Local Variables:
dnl mode: autoconf
dnl End:
