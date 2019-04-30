AC_INIT([icon], [2.6.0-rc], [icon@mpimet.mpg.de], [icon-2.6.0-rc])
AC_PREREQ([2.69])

AC_CONFIG_SRCDIR([src/drivers/icon.f90])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AC_PREFIX_DEFAULT([$(pwd)/install])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile depgen.f90.config])

AC_PROG_FC
AC_LANG([Fortran])
AC_FC_PP_SRCEXT([f90])
AC_FC_FREEFORM

ACX_FC_VENDOR
ACX_FC_PP_DEFINE
ACX_FC_MODULE_INC_FLAG
ACX_FC_MODULE_OUT_FLAG

ACX_FC_MODULE_FILE_NAMING(
  [FC_MODULE_UPPER=$acx_cv_fc_module_file_naming_upper
   FC_MODULE_EXT=$acx_cv_fc_module_file_naming_ext],
  [FC_MODULE_UPPER=no
   FC_MODULE_EXT=mod])
AC_SUBST([FC_MODULE_UPPER])
AC_SUBST([FC_MODULE_EXT])

ACX_FC_INCLUDE_ORDER(
  [FC_INC_ORDER=$acx_cv_fc_include_order],
  [FC_INC_ORDER=src,flg])
AC_SUBST([FC_INC_ORDER])

ACX_FC_HASH_INCLUDE_ORDER(
  [PP_INC_ORDER_f90=$acx_cv_fc_hash_include_order],
  [PP_INC_ORDER_f90=inc,flg])
AC_SUBST([PP_INC_ORDER_f90])

AC_MSG_CHECKING([whether the Fortran interface of mtime library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mtime
      character(len=max_date_str_len) :: test_date_string
      type(date), pointer :: test_date
      test_date => newdate(0,0,0)
      call dateToString(test_date, test_date_string)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([unable to compile test program using mtime: set FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether Fortran self 0.2 library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use self_object
      use self_vector_ref
      use self_vector
      use self_map
      use self_assert
      implicit none
      type(map) :: a
      call a%init()
]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([unable to compile test program using self 0.2 library: set FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether Fortran BLAS library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      real :: a(2,2), b(2,2), c(2,2)
      call dgemm('n','n',2,2,2,1.,a,2,b,2,0.,c,2)
]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([unable to compile test program using BLAS: set FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether Fortran LAPACK library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      real :: a(2,2), b(2,2), c(2,2), d(2), e(2), f(2)
      call dgeev('n','n',2,a,2,d,e,b,2,c,2,f,2,0)
]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([unable to compile test program using LAPACK: set FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether the Fortran 77 interface of NetCDF library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      include 'netcdf.inc'
      integer :: status
      status = nf_open('conftest.nc',0,0)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([unable to compile test program using NetCDF: set FCFLAGS, LDGLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether the Fortran 2003 interface of CDI 1.8.x library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: TSTEP_CONSTANT, cdiInqNatts, &
                        vlistInqVarTypeOfGeneratingProcess 
      implicit none
      integer :: a, b, c, d
! Appeared in 1.8.x
      a = cdiInqNatts(b, c, d)
! Disappeared in 1.9.x
      a = vlistInqVarTypeOfGeneratingProcess(b, c)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([unable to compile test program using CDI: set FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether the Fortran 2003 interface of CDI 1.8.x is patched])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: t_CdiIterator, cdiIterator_inqFiletype
      implicit none
      integer :: a
      type(t_CdiIterator) :: b
      a = cdiIterator_inqFiletype(b)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([the provided version of CDI 1.8.x does not have required patches])])

FCFLAGS="${FC_DEFINE}__ICON__ ${FC_DEFINE}__NO_JSBACH__ ${FC_DEFINE}__NO_ICON_TESTBED__ $FCFLAGS"

AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--enable-mpi], [enable MPI (parallelization) support @<:@default=yes@:>@])],
  [], [enable_mpi=yes])
AS_IF([test "x$enable_mpi" = "xyes"],
  [AC_MSG_CHECKING([whether Fortran MPI library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
       use mpi, only: mpi_init
       integer :: error
       call mpi_init(error)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_ERROR([unable to compile test program using MPI: \
either disable MPI support (--disable-mpi) or set FC to the MPI Fortran compiler wrapper, \
or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
  [FCFLAGS="${FC_DEFINE}NOMPI $FCFLAGS"])

AC_ARG_ENABLE([atmo],
  [AS_HELP_STRING([--enable-atmo], [enable the atmosphere component of ICON @<:@default=yes@:>@])],
  [], [enable_atmo=yes])
AS_IF([test "x$enable_atmo" = "xyes"], [],
  [AC_MSG_NOTICE([the atmosphere component is disabled])
  FCFLAGS="${FC_DEFINE}__NO_ICON_ATMO__ $FCFLAGS"])

AC_ARG_ENABLE([ocean],
  [AS_HELP_STRING([--enable-ocean], [enable the ocean component of ICON @<:@default=yes@:>@])],
  [], [enable_ocean=yes])
AS_IF([test "x$enable_ocean" = "xyes"], [],
  [AC_MSG_NOTICE([the ocean component is disabled])
  FCFLAGS="${FC_DEFINE}__NO_ICON_OCEAN__ $FCFLAGS"])

AC_PROG_CC
AC_LANG([C])
dnl Supporting ISO C89 means that the following header files are available:
dnl stdlib.h, stdio.h, math.h, time.h, ctype.h, string.h, signal.h, limits.h
dnl errno.h, ctype.h assert.h.
AS_IF([test "x$ac_cv_prog_cc_c89" = "xno"],
  [AC_MSG_ERROR([C compiler does not support ISO C89])])

dnl Some of the C source files require support for ISO C99. Supporting ISO C99
dnl also means that the following header files are available:
dnl stdbool.h.
AC_PROG_CC_C99
AS_IF([test "x$ac_cv_prog_cc_c99" = "xno"],
  [AC_MSG_ERROR([C compiler does not support ISO C99])])

dnl Check for headers that are protected with HAVE_* macros in the code
dnl (override AC_INCLUDES_DEFAULT of autoconf 2.69, which runs additional
dnl annoying checks; already fixed in the autoconf master branch).
AC_CHECK_HEADERS([execinfo.h ucontext.h link.h],,,[[// override default includes]])

dnl Check for headers that are not protected at all.
AC_CHECK_HEADERS([unistd.h fcntl.h pwd.h stdint.h sys/types.h sys/stat.h]dnl
  [sys/param.h sys/resource.h sys/utsname.h sys/time.h sys/times.h dirent.h],,
  [AC_MSG_ERROR([unable to find header file $ac_header])],[[// override default includes]])

dnl This should depend on the $acx_cv_fc_vendor
AS_VAR_APPEND([CPPFLAGS], [" -Df2cFortran"])

dnl Other header files that are protected in the code but not with the HAVE_* macros:
dnl cfortran.h: fortran.h malloc.h descrip.h
dnl util_hash.c: endian.h
dnl util_sysinfo.c: mcheck.h sys/systemcfg.h
dnl util_timer.c: sys/systemcfg.h asm/msr.h

dnl Header files that are included from non-mandatory files
dnl util_pbio.c: ffio.h pthread.h
dnl util_postgresql.c: libpq-fe.h
dnl util_signal.c: features.h strings.h fenv.h fptrap.h

AC_ARG_VAR([PYTHON], [Python interpreter command (used for Fortran dependency generation)])
ACX_CHECK_PROGS([PYTHON], [python python3], [for Python interpreter command],
  [acx_depgen_try='$acx_prog_candidate $srcdir/utils/mkhelper/depgen.py -h > /dev/null'
  AC_TRY_EVAL([acx_depgen_try])],
  [unknown])
AS_IF([test "x$PYTHON" = "xunknown"],
  [AC_MSG_ERROR([unable to find Python interpreter])])

AC_ARG_VAR([PERL], [Perl interpreter command (used to generate version information)])
ACX_CHECK_PROGS([PERL], [perl], [for Perl interpreter command],
  [acx_perl_try='$acx_prog_candidate $srcdir/utils/pvcs.pl --help > /dev/null'
  AC_TRY_EVAL([acx_perl_try])
  test $ac_status -eq 0 || test $ac_status -eq 1],
  [unknown])
AS_IF([test "x$PERL" = "xunknown"],
  [AC_MSG_ERROR([unable to find Perl interpreter])])

AC_PROG_INSTALL

AC_OUTPUT
