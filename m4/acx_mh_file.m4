dnl acx_mh_file.m4 --- load machine and site specific configurations
dnl
dnl Copyright  (C)  2017  Luis Kornblueh <luis.kornblueh@mpimet.mpg.de>
dnl
dnl Version: 1
dnl Keywords:
dnl Author: Luis Kornblueh <luis.kornblueh@mpimet.mpg.de>
dnl Maintainer: Luis Kornblueh <luis.kornblueh@mpimet.mpg.de>
dnl URL: https://code.zmaw.de/projects/icon
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
dnl Neither the name of the Max PLanck Institute for Meteorology nor the
dnl names of its contributors may be used to endorse or promote products
dnl derived from this software without specific prior written permission.
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
dnl
AC_DEFUN([ACX_GET_MH],
         [AC_MSG_CHECKING([for machine dependent configuration])
          AC_REQUIRE([AC_PROG_SED]) 
          changequote(,)
          cat > confsed <<EOF
s/ *= */="/
/="/s/$/"/
/^#/{
    d
}
/[^LIBS].*=/{
    p
    s/=.*$//
    s/^/export /
    p
    d
}
/LIBS.*\\\/ {
    s/\\\/ /g
    s/"$//
    h
    D
    bnext
}
/^ *\-L.*\\\/{
    s/\\\/ /g
    s/^ *//g
    s/"$//
    H
    D
    bnext
}
/^ *\-L/{
    s/^ *//g
    s/"$//
    H
    x
    s/\n//g
    s/$/"/
    p
    a\\
export LIBS
    d
}

:next
EOF
          changequote([,])
          $SED -f confsed $1 > conftest 
          . ./conftest 
          /bin/rm -f confsed conftest
          AS_IF([test "$1" != 0],
                [AC_MSG_RESULT([$1])],
                [AC_MSG_RESULT(unavailable)])
])
dnl
dnl
dnl
AC_DEFUN([ACX_GET_MODULES_TO_LOAD],
         [AC_REQUIRE([AC_PROG_SED]) 
          changequote(,)
          cat > confsed <<EOF
s/ *= */="/
/="/s/$/"/
/^#/{
    d
}
/[^LIBS].*=/{
    p
    s/=.*$//
    s/^/export /
    p
    d
}
/LIBS.*\\\/ {
    s/\\\/ /g
    s/"$//
    h
    D
    bnext
}
/^ *\-L.*\\\/{
    s/\\\/ /g
    s/^ *//g
    s/"$//
    H
    D
    bnext
}
/^ *\-L/{
    s/^ *//g
    s/"$//
    H
    x
    s/\n//g
    s/$/"/
    p
    a\\
export LIBS
    d
}

:next
EOF
          changequote([,])
          $SED -f confsed $1 > conftest
          export fortran_compiler
          export ac_sitename
          export hostname
          export with_mpi
          export host
          load_modules=$($CONFIG_SHELL -c '. ./conftest; echo "$load_modules"')
          /bin/rm -f confsed conftest
])
dnl
dnl
dnl
AC_DEFUN([ACX_GET_PROFILE_TO_LOAD],
         [AC_REQUIRE([AC_PROG_SED]) 
          changequote(,)
          cat > confsed <<EOF
s/ *= */="/
/="/s/$/"/
/^#/{
    d
}
/[^LIBS].*=/{
    p
    s/=.*$//
    s/^/export /
    p
    d
}
/LIBS.*\\\/ {
    s/\\\/ /g
    s/"$//
    h
    D
    bnext
}
/^ *\-L.*\\\/{
    s/\\\/ /g
    s/^ *//g
    s/"$//
    H
    D
    bnext
}
/^ *\-L/{
    s/^ *//g
    s/"$//
    H
    x
    s/\n//g
    s/$/"/
    p
    a\\
export LIBS
    d
}

:next
EOF
          changequote([,])
          $SED -f confsed $1 > conftest
          export fortran_compiler
          export ac_sitename
          export hostname
          export with_mpi
          export host
          load_profile=$($CONFIG_SHELL -c '. ./conftest; echo "$load_profile"')
          /bin/rm -f confsed conftest
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://code.zmaw.de/projects/icon/wiki/ICON_licensees"
dnl license-default: "bsd"
dnl End:


