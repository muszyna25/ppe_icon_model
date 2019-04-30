AC_PREREQ([2.69])
AC_INIT([icon], [2.6.0-rc], [icon@mpimet.mpg.de], [icon-2.6.0-rc])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/drivers/icon.f90])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile depgen.f90.config depgen.c.config])

AC_PREFIX_DEFAULT(["`pwd`/install"])

AC_PROG_FC
AC_LANG([Fortran])
AC_FC_PP_SRCEXT([f90])

AC_SUBST([DEPGEN_DEFS])
AC_SUBST([DEPGEN_INTRINSIC_MODS],
  [[iso_c_binding,iso_fortran_env,ieee_exceptions,ieee_arithmetic,]]dnl
[[ieee_features,omp_lib,omp_lib_kinds,openacc]])
AC_SUBST([DEPGEN_EXTERNAL_MODS],
  [[mo_cdi,mpi,mtime,mtime_datetime,self_object,self_vector_ref,]]dnl
[[self_vector,self_map,self_assert,sct,mo_yac_finterface,yaxt]])

ACX_LANG_VENDOR

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FC_MOD_IN], [$acx_cv_fc_module_in_flag])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FC_MOD_OUT], [$acx_cv_fc_module_out_flag])])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FC_MOD_FILE_UPPER], [$acx_cv_fc_module_naming_upper])
   AC_SUBST([FC_MOD_FILE_EXT], [$acx_cv_fc_module_naming_ext])])

ACX_FC_INCLUDE_FLAG(
  [AC_SUBST([FC_INC_FLAG], [$acx_cv_fc_ftn_include_flag])])
ACX_FC_INCLUDE_ORDER(
  [FC_INC_ORDER=$acx_cv_fc_ftn_include_order],
  [FC_INC_ORDER=src,flg])
AC_SUBST([FC_INC_ORDER])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FC_INC_FLAG_PP_f90], [$acx_cv_fc_pp_include_flag])])
ACX_FC_INCLUDE_ORDER_PP(
  [FC_INC_ORDER_PP_f90=$acx_cv_fc_pp_include_order],
  [FC_INC_ORDER_PP_f90=inc,flg])
AC_SUBST([FC_INC_ORDER_PP_f90])

ACX_LANG_MACRO_FLAG([AC_SUBST([FC_PP_DEF], [$acx_cv_fc_macro_flag])])

AC_ARG_ENABLE([silent-rules],
  [AC_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=yes@:>@])],
dnl We pass this variable to the configure scripts of the bundled libraries,
dnl therefore, we need to make sure that it's either 'yes' or 'no'.
     [AS_IF([test x"$enablevar" != xno], [enable_silent_rules=yes])],
  [enable_silent_rules=yes])
AS_VAR_IF([enable_silent_rules], [yes],
  [DEFAULT_VERBOSITY=0], [DEFAULT_VERBOSITY=1])
AC_SUBST([DEFAULT_VERBOSITY])

AC_MSG_CHECKING([whether Fortran BLAS library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      real :: a(2,2), b(2,2), c(2,2)
      call dgemm('n','n',2,2,2,1.,a,2,b,2,0.,c,2)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to compile test program using BLAS: set FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether Fortran LAPACK library is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      real :: a(2,2), b(2,2), c(2,2), d(2), e(2), f(2)
      call dgeev('n','n',2,a,2,d,e,b,2,c,2,f,2,0)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to compile test program using LAPACK: set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether the Fortran 77 interface of NetCDF library dnl
is available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      include 'netcdf.inc'
      integer :: status
      status = nf_open('conftest.nc',0,0)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to compile test program using NetCDF: set dnl
FCFLAGS, LDGLAGS and LIBS accordingly])])

AC_SUBST([OPTION_FCFLAGS],
  ["${FC_PP_DEF}__ICON__ ${FC_PP_DEF}__NO_ICON_TESTBED__"])

AC_ARG_ENABLE([mpi],
  [AC_HELP_STRING([--enable-mpi],
     [enable MPI (parallelization) support @<:@default=yes@:>@])], [],
     [enable_mpi=yes])
AS_VAR_IF([enable_mpi], [yes],
  [AC_MSG_CHECKING([whether Fortran MPI library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
       use mpi, only: mpi_init
       integer :: error
       call mpi_init(error)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to compile test program using MPI: either dnl
disable MPI support (--disable-mpi) or set FC to the MPI Fortran compiler dnl
wrapper, or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
  [AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}NOMPI"])])

AC_ARG_ENABLE([openmp],
  [AC_HELP_STRING([--enable-openmp],
     [enable OpenMP (parallelization) support @<:@default=yes@:>@])], [],
     [enable_openmp=yes])
AS_VAR_IF([enable_openmp], [yes],
  [ACX_LANG_OPENMP_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openmp_flag"])
      ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
      AS_VAR_IF([acx_macro_defined], [yes],
        [ACX_LANG_MACRO_CHECK_VALUE([_OPENMP],
           [201511 201307 201107 200805 200505])
         AS_VAR_IF([acx_macro_value], [unknown],
           [AC_MSG_WARN([unable to detect supported OpenMP version])
            AS_VAR_APPEND([DEPGEN_DEFS], [" ${FC_PP_DEF}_OPENMP"])],
           [AS_IF([test "$acx_macro_value" -lt 201307],
              [AC_MSG_WARN([Fortran compiler does not support OpenMP 4.0 dnl
(201307): it is strongly recommended to disable OpenMP support dnl
(--disable-openmp)])])
            AS_VAR_APPEND([DEPGEN_DEFS],
              [" ${FC_PP_DEF}_OPENMP=$acx_macro_value"])])],
         [AC_MSG_FAILURE([Fortran compiler does not define the _OPENMP dnl
preprocessor macro])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenMP support by default: set dnl
FCFLAGS accordingly to disable it])])])

AC_ARG_ENABLE([atmo],
  [AC_HELP_STRING([--enable-atmo],
     [enable the atmosphere component of ICON @<:@default=yes@:>@])], [],
  [enable_atmo=yes])
AS_VAR_IF([enable_atmo], [yes], [],
  [AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_ATMO__"])
   AC_MSG_NOTICE([the atmosphere component is disabled])])

AC_ARG_ENABLE([ocean],
  [AC_HELP_STRING([--enable-ocean],
     [enable the ocean component of ICON @<:@default=yes@:>@])], [],
  [enable_ocean=yes])
AS_VAR_IF([enable_ocean], [yes], [],
  [AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_OCEAN__"])
   AC_MSG_NOTICE([the ocean component is disabled])])

AC_ARG_ENABLE([jsbach],
  [AC_HELP_STRING([--enable-jsbach],
     [enable the land component JSBACH @<:@default=yes@:>@])], [],
  [enable_jsbach=yes])
AC_SUBST([JSBACH_ENABLED])
AS_VAR_IF([enable_jsbach], [yes],
  [JSBACH_ENABLED=
   AC_CONFIG_FILES([externals/jsbach.mk])],
  [JSBACH_ENABLED='#'
   AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}__NO_JSBACH__"])
   AC_MSG_NOTICE([the land component JSBACH is disabled])])

AC_ARG_ENABLE([coupling],
  [AC_HELP_STRING([--enable-coupling],
     [enable the coupling with YAC @<:@default=yes@:>@])], [],
  [enable_coupling=yes])
AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}YAC_coupling"])],
  [AC_MSG_NOTICE([the coupling is disabled])])

AC_ARG_ENABLE([psrad],
  [AC_HELP_STRING([--enable-psrad],
     [enable the PSRAD component of ICON @<:@default=yes@:>@])], [],
  [enable_psrad=yes])
AS_VAR_IF([enable_psrad], [yes], [],
  [AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_PS_RAD__"])
   AC_MSG_NOTICE([the PSRAD component is disabled])])

AC_ARG_ENABLE([yaxt],
  [AC_HELP_STRING([--enable-yaxt],
     [enable the YAXT data exchange @<:@default=no@:>@])], [],
  [enable_yaxt=no])
AS_VAR_IF([enable_yaxt], [yes],
  [AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}HAVE_YAXT"])])

AC_ARG_ENABLE([sct],
  [AC_HELP_STRING([--enable-sct],
     [enable the SCT timer @<:@default=no@:>@])], [],
  [enable_sct=no])
AS_VAR_IF([enable_sct], [yes],
  [AC_MSG_CHECKING([whether the Fortran interface of SCT library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use sct, only: sct_init
      implicit none
      call sct_init(timer_max=512)]])],
     [AC_MSG_RESULT([yes])
      AS_VAR_APPEND([OPTION_FCFLAGS], [" ${FC_PP_DEF}__SCT__"])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to compile test program using SCT: dnl
set FCFLAGS, LDFLAGS and LIBS accordingly])])])

AC_ARG_WITH([external-mtime],
  [AC_HELP_STRING([--with-external-mtime],
     [use an external mtime library instead of the bundled one ]dnl
[@<:@default=no@:>@])], [],
  [with_external_mtime=no])

AC_ARG_WITH([external-self],
  [AC_HELP_STRING([--with-external-self],
     [use an external self library instead of the bundled one ]dnl
[@<:@default=no@:>@])], [],
  [with_external_self=no])

AC_ARG_WITH([external-cdi],
  [AC_HELP_STRING([--with-external-cdi],
     [use an external CDI library instead of the bundled one ]dnl
[@<:@default=no@:>@])], [],
  [with_external_cdi=no])

AC_ARG_WITH([external-yac],
  [AC_HELP_STRING([--with-external-yac],
     [use an external YAC library (when the coupling is enabled) instead ]dnl
[of the bundled one @<:@default=no@:>@])], [],
  [with_external_yac=no])

AC_ARG_WITH([external-yaxt],
  [AC_HELP_STRING([--with-external-yaxt],
     [use an external YAXT library (when the YAXT data exchange is ]dnl
[enabled) instead of the bundled one @<:@default=no@:>@])], [],
  [with_external_yaxt=no])

AC_SUBST([BUNDLED_FCFLAGS])
AC_SUBST([BUNDLED_LIBS])

dnl Set $ac_abs_top_srcdir and $ac_abs_top_builddir
_AC_SRCDIRS([.])

AS_VAR_IF([with_external_self], [no],
  [AS_IF([test -f "${srcdir}/externals/self/configure"],
     [AX_SUBDIRS_CONFIGURE([externals/self], [], [],
        [[--enable-silent-rules=$enable_silent_rules]],
        [[--enable-silent-rules], [--disable-silent-rules]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/self/src"])
      AS_VAR_APPEND([BUNDLED_LIBS], [" externals/self/src/libself.a"])],
     [AC_MSG_ERROR([unable to use the bundled version of self library: dnl
file ${srcdir}/externals/self/configure not found])])],
  [AC_MSG_CHECKING([whether Fortran self 0.2 library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use self_object
      use self_vector_ref
      use self_vector
      use self_map
      use self_assert
      implicit none
      type(map) :: a
      call a%init()]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to compile test program using self 0.2 dnl
library: enable usage of the bundled version of the library or set FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_IF([with_external_yac], [no],
     [AS_IF([test -f "${srcdir}/externals/yac/configure"],
        [AS_VAR_IF([with_external_mtime], [no],
           [yac_MTIME_CFLAGS="-I${ac_abs_top_srcdir}/externals/mtime/include"
            yac_MTIME_CLIBS="${ac_abs_top_builddir}/externals/mtime/src/dnl
.libs/libmtime.a"],
           [yac_MTIME_CFLAGS=
            yac_MTIME_CLIBS=])
         AX_SUBDIRS_CONFIGURE([externals/yac], [], [],
           [[--enable-lib-only=yes], [--enable-xml=yes],
            [--enable-netcdf=yes], [--enable-mpi=$enable_mpi],
            [--with-external-lapack=fortran], [--with-external-mtime=yes],
            [--enable-silent-rules=$enable_silent_rules],
            [MTIME_CFLAGS=$yac_MTIME_CFLAGS], [MTIME_CLIBS=$yac_MTIME_CLIBS],
            [NETCDF_CFLAGS=], [NETCDF_CLIBS=],
            [XML2_CFLAGS=], [XML2_CLIBS=],
            [MPI_CFLAGS=], [MPI_CLIBS=],
            [MPI_FCFLAGS=], [MPI_FCLIBS=],
            [FORTRAN_LAPACK_CLIBS=]],
           [[--enable-lib-only], [--disable-lib-only],
            [--enable-xml], [--disable-xml],
            [--enable-netcdf], [--disable-netcdf],
            [--enable-mpi], [--disable-mpi],
            [--with-external-lapack], [--without-external-lapack],
            [--with-external-mtime], [--without-external-mtime],
            [--enable-silent-rules], [--disable-silent-rules]])
         AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/yac/src"])
         AS_VAR_APPEND([BUNDLED_LIBS], [" externals/yac/src/libyac.a"])],
        [AC_MSG_ERROR([unable to use the bundled version of YAC library: dnl
file ${srcdir}/externals/yac/configure not found])])],
     [AC_MSG_CHECKING([whether the Fortran interface of YAC library is dnl
available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_yac_finterface, only: yac_finit
      implicit none
      character(len=1024) :: a, b
      call yac_finit(a, b)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to compile test program using YAC: dnl
set FCFLAGS, LDFLAGS and LIBS accordingly])])])])

AS_VAR_IF([with_external_mtime], [no],
  [AS_IF([test -f "${srcdir}/externals/mtime/configure"],
     [AX_SUBDIRS_CONFIGURE([externals/mtime], [], [],
        [[--enable-static=yes], [--enable-shared=no],
         [--enable-examples=no], [--enable-silent-rules=$enable_silent_rules]],
        [[--enable-static], [--disable-static],
         [--enable-shared], [--disable-shared],
         [--enable-examples], [--disable-examples],
         [--enable-silent-rules], [--disable-silent-rules]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/mtime/src"])
      AS_VAR_APPEND([BUNDLED_LIBS],
        [" externals/mtime/src/.libs/libmtime.a"])],
     [AC_MSG_ERROR([unable to use the bundled version of mtime library: dnl
file ${srcdir}/externals/mtime/configure not found])])],
  [AC_MSG_CHECKING(
     [whether the Fortran interface of mtime library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mtime
      character(len=max_date_str_len) :: test_date_string
      type(date), pointer :: test_date
      test_date => newdate(0,0,0)
      call dateToString(test_date, test_date_string)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to compile test program using mtime library: dnl
enable usage of the bundled version of the library or set FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])])

AC_ARG_ENABLE([grib2],
  [AC_HELP_STRING([--enable-grib2],
     [enable GRIB2 I/O @<:@default=no@:>@])],
dnl We pass this variable to the configure scripts of libcdi, therefore, we
dnl need to make sure that it's either 'yes' or 'no'.
  [AS_IF([test x"$enablevar" != xno], [enable_grib2=yes])],
  [enable_grib2=no])

AS_VAR_IF([with_external_cdi], [no],
  [AS_IF([test -f "${srcdir}/externals/cdi/configure"],
     [AX_SUBDIRS_CONFIGURE([externals/cdi], [], [],
dnl Currently, we don't have full control on choosing the GRIB2 backend (e.g.
dnl grib_api, eccodes) due to the limitations of the configure script of
dnl libcdi. Therefore, we set --with-grib_api=$enable_grib2 --with-eccodes=no
dnl unconditionally. This way, the configure script of libcdi enables support
dnl for grib2 regardless of whether grib_api or eccodes librariy is provided
dnl with the LDFLAGS and LIBS variables.
        [[--enable-static=yes], [--enable-shared=no],
         [--enable-silent-rules=$enable_silent_rules],
         [--enable-iso-c-interface=yes], [--enable-cf-interface=no],
         [--enable-cdi-app=no],  [--enable-grib=yes],
         [--enable-cgribex=yes], [--with-grib_api=$enable_grib2],
         [--with-eccodes=no], [--with-threads=no],
         [--with-netcdf=yes]],
        [[--enable-static], [--disable-static],
         [--enable-shared], [--disable-shared],
         [--enable-silent-rules], [--disable-silent-rules],
         [--enable-iso-c-interface], [--disable-iso-c-interface],
         [--enable-cf-interface], [--disable-cf-interface],
         [--enable-cdi-app], [--disable-cdi-app],
         [--enable-grib], [--disable-grib],
         [--enable-cgribex], [--disable-cgribex],
         [--with-grib_api], [--without-grib_api],
         [--with-eccodes], [--without-eccodes],
         [--with-threads], [--without-threads],
         [--with-netcdf], [--without-netcdf],
         [--enable-mpi], [--disable-mpi]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/cdi/src"])
      AS_VAR_APPEND([BUNDLED_LIBS],
        [" externals/cdi/src/.libs/libcdi_f2003.a externals/cdi/src/dnl
.libs/libcdi.a"])],
     [AC_MSG_ERROR([unable to use the bundled version of CDI library: dnl
file ${srcdir}/externals/cdi/configure not found])])],
  [AC_MSG_CHECKING([whether the Fortran 2003 interface of CDI 1.8.x library dnl
is available])
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
      AC_MSG_FAILURE([unable to compile test program using CDI library: dnl
enable usage of the bundled version of the library or set FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])
   AC_MSG_CHECKING(
     [whether the provided version of CDI 1.8.x is patched])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: t_CdiIterator, cdiIterator_inqFiletype
      implicit none
      integer :: a
      type(t_CdiIterator) :: b
      a = cdiIterator_inqFiletype(b)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([the provided version of CDI 1.8.x does not have dnl
required patches])])
   AS_VAR_IF([enable_grib2], [yes],
     [AC_MSG_CHECKING(
        [whether the provided version of CDI is built with GRIB2 support])
      icon_cdi_supports_grib2=unknown
      AS_VAR_IF([cross_compiling], [no],
        [AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: cdiHaveFiletype, CDI_FILETYPE_GRB2
      write(*, "(i0)") cdiHaveFiletype(CDI_FILETYPE_GRB2)]])],
           [icon_exec_result=`./conftest$ac_exeext 2>/dev/null`
            AS_IF([test $? -eq 0],
              [AS_IF([test $icon_exec_result -eq 1],
                 [icon_cdi_supports_grib2=yes],
                 [icon_cdi_supports_grib2=no])])])])
      AC_MSG_RESULT([$icon_cdi_supports_grib2])
      AS_CASE([$icon_cdi_supports_grib2],
        [no],
        [AC_MSG_FAILURE([the provided version of CDI does not support GRIB2])],
        [unknown],
        [AC_MSG_WARN([unable to check whether the provided version of dnl
CDI supports GRIB2])])])])

AS_VAR_IF([enable_yaxt], [yes],
  [AS_VAR_IF([with_external_yaxt], [no],
     [AS_IF([test -f "${srcdir}/externals/yaxt/configure"],
        [AX_SUBDIRS_CONFIGURE([externals/yaxt], [], [],
           [[--enable-static=yes], [--enable-shared=no],
            [--enable-dependency-tracking=yes],
            [--enable-silent-rules=$enable_silent_rules]],
           [[--enable-static], [--disable-static],
            [--enable-shared], [--disable-shared],
            [--enable-dependency-tracking], [--disable-dependency-tracking],
            [--enable-silent-rules], [--disable-silent-rules]])
         AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/yaxt/src"])
         AS_VAR_APPEND([BUNDLED_LIBS],
           [" externals/yaxt/src/.libs/libyaxt.a externals/yaxt/src/dnl
.libs/libyaxt_c.a"])],
        [AC_MSG_ERROR([unable to use the bundled version of YAXT library: dnl
file ${srcdir}/externals/yaxt/configure not found])])],
     [AC_MSG_CHECKING([whether the Fortran interface of YAXT library dnl
is available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use yaxt, only: xt_initialized, xt_initialize
      implicit none
      integer :: a
      if (.not. xt_initialized()) call xt_initialize(a)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to compile test program using YAXT library: dnl
enable usage of the bundled version of the library or set FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])])])

AC_PROG_CC
dnl Supporting ISO C89 means that the following header files are available:
dnl stdlib.h, stdio.h, math.h, time.h, ctype.h, string.h, signal.h, limits.h
dnl errno.h, ctype.h assert.h.
AS_VAR_IF([ac_cv_prog_cc_c89], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C89])])

dnl Some of the C source files require support for ISO C99. Supporting ISO C99
dnl also means that the following header files are available:
dnl stdbool.h.
AC_LANG([C])
AC_PROG_CC_C99
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])

dnl Check for headers that are protected with HAVE_* macros in the code
dnl (override AC_INCLUDES_DEFAULT of autoconf 2.69, which runs additional
dnl annoying checks; already fixed in the autoconf master branch).
AC_CHECK_HEADERS([execinfo.h ucontext.h link.h], [], [],
  [[// override default includes]])

dnl Check for headers that are not protected at all.
AC_CHECK_HEADERS([unistd.h fcntl.h pwd.h stdint.h sys/types.h sys/stat.h]dnl
  [sys/param.h sys/resource.h sys/utsname.h sys/time.h sys/times.h dirent.h],
  [],
  [AC_MSG_FAILURE([unable to find header file $ac_header])],
  [[// override default includes]])

dnl Other header files that are protected in the code but
dnl not with the HAVE_* macros:
dnl util_hash.c: endian.h
dnl util_sysinfo.c: mcheck.h sys/systemcfg.h
dnl util_timer.c: sys/systemcfg.h asm/msr.h

dnl Header files that are included from non-mandatory files
dnl util_pbio.c: ffio.h pthread.h
dnl util_postgresql.c: libpq-fe.h
dnl util_signal.c: features.h strings.h fenv.h fptrap.h

AC_ARG_VAR([PYTHON],
  [Python interpreter command (used for Fortran dependency generation)])
ACX_PROG_SEARCH([PYTHON], [python python3],
  [acx_depgen_try=dnl
'$acx_candidate $srcdir/utils/mkhelper/depgen.py -h > /dev/null'
   AC_TRY_EVAL([acx_depgen_try])],
  [PYTHON=$acx_candidate])

AC_ARG_VAR([PERL],
  [Perl interpreter command (used to generate version information)])
ACX_PROG_SEARCH([PERL], [perl],
  [acx_depgen_try='$acx_candidate $srcdir/utils/pvcs.pl --help > /dev/null'
   AC_TRY_EVAL([acx_depgen_try])
   test $ac_status -eq 0 || test $ac_status -eq 1],
  [PERL=$acx_candidate])

AC_PROG_INSTALL
AC_PROG_MKDIR_P

AC_OUTPUT
