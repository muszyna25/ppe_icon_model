dnl AC_TEST_FORTTYPES tests to see if the following fortran datatypes are
dnl supported: INTEGER1, INTEGER2, INTEGER4, INTEGER8, REAL4, REAL8,
dnl            DOUBLE_COMPLEX
dnl
define(AC_TEST_FORTTYPES,dnl
   [
FIX_FILE=0
FORT_INT1=1
FORT_INT2=1
FORT_INT4=1
FORT_INT8=1
FORT_REAL4=1
FORT_REAL8=1
FORT_DOUBLE_COMPLEX=1
COUNT=13
AC_MSG_CHECKING(for integer*1)
cat > testfort.f <<EOF
        subroutine forttype( a )
        integer*1 a
        return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_INT1=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for integer*2)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        integer*2 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_INT2=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for integer*4)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        integer*4 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_INT4=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for integer*8)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        integer*8 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_INT8=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for integer*16)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        integer*16 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_INT16=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for real*4)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        real*4 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_REAL4=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for real*8)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        real*8 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_REAL8=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for real*16)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        real*16 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_REAL16=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for double complex)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        double complex a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_DOUBLE_COMPLEX=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for complex*8)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        complex*8 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_COMPLEX8=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for complex*16)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        complex*16 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_COMPLEX16=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
AC_MSG_CHECKING(for complex*32)
    cat > testfort.f <<EOF
        subroutine forttype( a )
        complex*32 a
	return
        end
EOF
   $FC -c testfort.f > /dev/null 2>&1
   if test ! -s testfort.o ; then
       AC_MSG_RESULT(no)
       FORT_COMPLEX32=0
       FIX_FILE=1
       COUNT=`expr ${COUNT} - 1`
   else
       AC_MSG_RESULT(yes)
   fi
   /bin/rm -f testfort.f testfort.o
dnl
   ])dnl
dnl
dnl AC_FORTRAN_GET_INTEGER_SIZE(var_for_size)
dnl
dnl sets var_for_size to the size.  Ignores if the size cannot be determined
dnl
define(AC_FORTRAN_GET_INTEGER_SIZE,
[AC_MSG_CHECKING([for size of Fortran INTEGER])
/bin/rm -f conftestval
/bin/rm -f conftestf.f conftestf.o conftestf
cat <<EOF > conftestf.f
      program conftest
      integer :: itest
      integer :: integer_byte
      integer_byte=bit_size(itest)/8	
      open(10,file="conftestval")
      write(10,*)integer_byte
      close(10)		
      end
EOF
Pac_CV_NAME=""
if $FC $FFLAGS -o conftestf conftestf.f >/dev/null 2>&1 ; then
    ./conftestf >/dev/null 2>&1
    Pac_CV_NAME=`cat conftestval`
fi
/bin/rm -f conftestf.f conftestf.o conftestf conftestval

if test -n "$Pac_CV_NAME" -a "$Pac_CV_NAME" != 0 ; then
    AC_MSG_RESULT($Pac_CV_NAME)
else
    AC_MSG_RESULT(unavailable)
fi
$1=$Pac_CV_NAME
])dnl
dnl
dnl AC_FORTRAN_GET_REAL_SIZE(var_for_size)
dnl
dnl sets var_for_size to the size.  Ignores if the size cannot be determined
dnl
define(AC_FORTRAN_GET_REAL_SIZE,
[AC_MSG_CHECKING([for size of Fortran REAL])
/bin/rm -f conftestval
/bin/rm -f conftestf.f conftestf.o conftestf
cat <<EOF > conftestf.f
      program conftest
      integer :: itest
      real    :: rtest
      integer :: integer_io_length, real_io_length
      integer :: integer_byte, real_byte
      inquire(iolength=integer_io_length) itest
      inquire(iolength=real_io_length) rtest	
      integer_byte=bit_size(itest)/8	
      real_byte = real_io_length/integer_io_length*integer_byte
      open(10,file="conftestval")
      write(10,*)real_byte
      close(10)		
      end
EOF
Pac_CV_NAME=""
if $FC $FFLAGS -o conftestf conftestf.f >/dev/null 2>&1 ; then
    ./conftestf >/dev/null 2>&1
    Pac_CV_NAME=`cat conftestval`
fi
/bin/rm -f conftestf.f conftestf.o conftestf conftestval

if test -n "$Pac_CV_NAME" -a "$Pac_CV_NAME" != 0 ; then
    AC_MSG_RESULT($Pac_CV_NAME)
else
    AC_MSG_RESULT(unavailable)
fi
$1=$Pac_CV_NAME
])dnl
dnl
dnl Get the format of Fortran names.  Uses F77, FFLAGS, and sets WDEF.
dnl If the test fails, sets NOF77 to 1, HAS_FORTRAN to 0
dnl
define(AC_GET_FORTNAMES,[
dnl Check for strange behavior of Fortran.  For example, some FreeBSD
dnl systems use f2c to implement f77, and the version of f2c that they
dnl use generates TWO (!!!) trailing underscores
dnl Currently, WDEF is not used but could be...
dnl
dnl Eventually, we want to be able to override the choices here and
dnl force a particular form.  This is particularly useful in systems
dnl where a Fortran compiler option is used to force a particular
dnl external name format (rs6000 xlf, for example).
cat > confftest.f <<EOF
      subroutine mpir_init_fop( a )
      integer a
      a = 1
      return
      end
EOF
$FC $FFLAGS -c confftest.f > /dev/null 2>&1
if test ! -s confftest.o ; then
    echo "Unable to test Fortran compiler"
    echo "(compiling a test program failed to produce an "
    echo "object file)."
    HAS_FORTRAN=0
elif test -z "$FORTRANNAMES" ; then
dnl We have to be careful here, since the name may occur in several
dnl forms.  We try to handle this by testing for several forms
dnl directly.
    if test $arch_CRAY ; then
dnl Cray doesn't accept -a ...
        nameform1=`strings confftest.o | grep mpir_init_fop_  | sed -n -e '1p'`
        nameform2=`strings confftest.o | grep MPIR_INIT_FOP   | sed -n -e '1p'`
        nameform3=`strings confftest.o | grep mpir_init_fop   | sed -n -e '1p'`
        nameform4=`strings confftest.o | grep mpir_init_fop__ | sed -n -e '1p'`
    else
        nameform1=`strings -a confftest.o | grep mpir_init_fop_  | sed -n -e '1p'`
        nameform2=`strings -a confftest.o | grep MPIR_INIT_FOP   | sed -n -e '1p'`
        nameform3=`strings -a confftest.o | grep mpir_init_fop   | sed -n -e '1p'`
        nameform4=`strings -a confftest.o | grep mpir_init_fop__ | sed -n -e '1p'`
    fi
    /bin/rm -f confftest.f confftest.o
    if test -n "$nameform4" ; then
	echo "Fortran externals are lower case and have 1 or 2 trailing underscores"
	FORTRANNAMES="FORTRANDOUBLEUNDERSCORE"
    elif test -n "$nameform1" ; then
        echo "Fortran externals have a trailing underscore and are lowercase"
	FORTRANNAMES="FORTRANUNDERSCORE"
    elif test -n "$nameform2" ; then
	echo "Fortran externals are uppercase"
	FORTRANNAMES="FORTRANCAPS"
    elif test -n "$nameform3" ; then
	echo "Fortran externals are lower case"
	FORTRANNAMES="FORTRANNOUNDERSCORE"
    else
	echo "Unable to determine the form of Fortran external names"
	echo "Make sure that the compiler $FC can be run on this system"
        HAS_FORTRAN=0
    fi
fi
])dnl
dnl
dnl AC_GET_MH
dnl
dnl Load machine specific compiler options
dnl
define(AC_GET_MH,
[AC_MSG_CHECKING([for machine dependent configuration])
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
sed -f confsed $1 > conftest 
. ./conftest 
/bin/rm -f confsed conftest
if test "$1" != 0 ; then
    AC_MSG_RESULT($1)
else
    AC_MSG_RESULT(unavailable)
fi
])dnl


