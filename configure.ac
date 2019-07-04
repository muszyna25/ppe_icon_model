AC_PREREQ([2.69])
AC_INIT([icon], [2.6.0-rc], [icon@mpimet.mpg.de], [icon-2.6.0-rc])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/drivers/icon.f90])

AC_PREFIX_DEFAULT(["`pwd`/install"])

AS_IF([test "`cd $srcdir && pwd`" != "`pwd`" && dnl
test -f "$srcdir/config.status"],
  [AC_MSG_ERROR([source directory already configured; dnl
run "make distclean" there first])])

dnl First, we declare all the configure options and make sure that they have
dnl valid values.

AC_ARG_ENABLE([rpaths],
  [AC_HELP_STRING([--enable-rpaths],
     [add directories specified with -L flags in LDFLAGS and LIBS to the
runtime library search paths (RPATH) @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rpaths=yes])],
  [enable_rpaths=yes])

AC_ARG_ENABLE([silent-rules],
  [AC_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_silent_rules=yes])],
  [enable_silent_rules=yes])

m4_divert_once([HELP_ENABLE], [[
Model Features:]])

AC_ARG_ENABLE([atmo],
  [AC_HELP_STRING([--enable-atmo],
     [enable the atmosphere component of ICON @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_atmo=yes])],
  [enable_atmo=yes])

AC_ARG_ENABLE([ocean],
  [AC_HELP_STRING([--enable-ocean],
     [enable the ocean component of ICON @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_ocean=yes])],
  [enable_ocean=yes])

AC_ARG_ENABLE([jsbach],
  [AC_HELP_STRING([--enable-jsbach],
     [enable the land component JSBACH @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_jsbach=yes])],
  [enable_jsbach=yes])

AC_ARG_ENABLE([coupling],
  [AC_HELP_STRING([--enable-coupling],
     [enable the coupling @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_coupling=yes])],
  [enable_coupling=yes])

AC_ARG_ENABLE([psrad],
  [AC_HELP_STRING([--enable-psrad],
     [enable the radiation component PSRAD @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_psrad=yes])],
  [enable_psrad=yes])

AC_ARG_ENABLE([ecrad],
  [AC_HELP_STRING([--enable-ecrad],
     [enable usage of the ECMWF radiation scheme (ECRAD) @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_ecrad=yes])],
  [enable_ecrad=no])

AC_ARG_ENABLE([rte-rrtmgp],
  [AC_HELP_STRING([--enable-rte-rrtmgp],
     [enable usage of the RTE+RRTMGP toolbox for radiation calculations
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rte_rrtmgp=yes])],
  [enable_rte_rrtmgp=no])

AC_ARG_ENABLE([rttov],
  [AC_HELP_STRING([--enable-rttov],
     [enable usage of the radiative transfer model for TOVS
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rttov=yes])],
  [enable_rttov=no])

AC_ARG_ENABLE([art],
  [AC_HELP_STRING([--enable-art],
     [enable the aerosols and reactive trace component ART
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_art=yes])],
  [enable_art=no])

m4_divert_once([HELP_ENABLE], [[
Infrastructural Features:]])

AC_ARG_ENABLE([mpi],
  [AC_HELP_STRING([--enable-mpi],
     [enable MPI (parallelization) support @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mpi=yes])],
  [enable_mpi=yes])

AC_ARG_ENABLE([mpi-checks],
  [AC_HELP_STRING([--enable-mpi-checks],
     [enable configure-time checks of MPI library for known defects ]dnl
[@<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mpi_checks=yes])],
  [enable_mpi_checks=yes])

AC_ARG_ENABLE([active-target-sync],
  [AC_HELP_STRING([--enable-active-target-sync],
     [enable MPI active target mode (otherwise, passive target mode is used)
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_active_target_sync=yes])],
  [enable_active_target_sync=no])

AC_ARG_ENABLE([openmp],
  [AC_HELP_STRING([--enable-openmp],
     [enable OpenMP support @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_openmp=yes])],
  [enable_openmp=yes])

AC_ARG_ENABLE([openacc],
  [AC_HELP_STRING([--enable-openacc],
     [enable OpenACC support @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_openacc=yes])],
  [enable_openacc=no])

AC_ARG_ENABLE([grib2],
  [AC_HELP_STRING([--enable-grib2],
     [enable GRIB2 I/O @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_grib2=yes])],
  [enable_grib2=no])

AC_ARG_ENABLE([parallel-netcdf],
  [AC_HELP_STRING([--enable-parallel-netcdf],
     [enable usage of the parallel features of NetCDF
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_parallel_netcdf=yes])],
  [enable_parallel_netcdf=no])

AC_ARG_ENABLE([sct],
  [AC_HELP_STRING([--enable-sct],
     [enable the SCT timer @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_sct=yes])],
  [enable_sct=no])

AC_ARG_ENABLE([yaxt],
  [AC_HELP_STRING([--enable-yaxt],
     [enable the YAXT data exchange @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_yaxt=yes])],
  [enable_yaxt=no])

AC_ARG_ENABLE([serialization],
[  --enable-serialization  enable Serialbox2 serialization. The value of the
                          argument must be one of the following @<:@default=no@:>@:
                            read|yes  enable READ mode
                            perturb   enable READ & PERTURB mode
                            create    enable CREATE mode
                            no        disable serialization],
[AS_CASE(["$enableval"],
   [yes], [enable_serialization=read],
   [read|perturb|create|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-serialization='$enableval'; valid values are 'read', 'perturb', dnl
'create', 'yes' (same as 'read'), 'no'])])],
[enable_serialization=no])

m4_divert_once([HELP_ENABLE], [[
Optimization Features:]])

ICON_OPTIMIZATION_ARG_ENABLE

AC_ARG_ENABLE([loop-exchange],
  [AC_HELP_STRING([--enable-loop-exchange],
     [enable loop exchange @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_loop_exchange=yes])],
  [enable_loop_exchange=yes])

AC_ARG_ENABLE([intel-consistency],
  [AC_HELP_STRING([--enable-intel-consistency],
     [enable Intel compiler directives enforcing consistency
@<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_intel_consistency=yes])],
  [enable_intel_consistency=auto])

AC_ARG_ENABLE([vectorized-lrtm],
  [AC_HELP_STRING([--enable-vectorized-lrtm],
     [enable the parallelization-invariant version of LRTM
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_vectorized_lrtm=yes])],
  [enable_vectorized_lrtm=no])

AC_ARG_ENABLE([mixed-precision],
  [AC_HELP_STRING([--enable-mixed-precision],
     [enable mixed precision dycore @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mixed_precision=yes])],
  [enable_mixed_precision=no])

dnl Optional Packages:

AC_ARG_WITH([external-ecrad],
  [AC_HELP_STRING([--with-external-ecrad],
     [use external ECRAD library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_ecrad=yes])],
  [with_external_ecrad=no])

AC_ARG_WITH([external-rte-rrtmgp],
     [AC_HELP_STRING([--with-external-rte-rrtmgp],
        [use external RTE+RRTMGP library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_rte_rrtmgp=yes])],
  [with_external_rte_rrtmgp=no])

AC_ARG_WITH([external-sct],
  [AC_HELP_STRING([--with-external-sct],
     [use external SCT library (when the SCT timer is enabled)
@<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_sct=yes])],
  [with_external_sct=no])

AC_ARG_WITH([external-yaxt],
  [AC_HELP_STRING([--with-external-yaxt],
     [use external YAXT library (when the YAXT data exchange is enabled)
@<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_yaxt=yes])],
  [with_external_yaxt=no])

AC_ARG_WITH([external-cdi],
  [AC_HELP_STRING([--with-external-cdi],
     [use external CDI library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_cdi=yes])],
  [with_external_cdi=no])

AC_ARG_WITH([external-mtime],
  [AC_HELP_STRING([--with-external-mtime],
     [use external MTIME library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_mtime=yes])],
  [with_external_mtime=no])

AC_ARG_WITH([external-yac],
  [AC_HELP_STRING([--with-external-yac],
     [use external YAC library (when the coupling is enabled)
@<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_yac=yes])],
  [with_external_yac=no])

AC_ARG_WITH([external-tixi],
  [AC_HELP_STRING([--with-external-tixi],
     [use external TIXI library (when ART is enabled)
@<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_tixi=yes])],
  [with_external_tixi=no])

AC_ARG_WITH([external-self],
  [AC_HELP_STRING([--with-external-self],
     [use external SELF library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_self=yes])],
  [with_external_self=no])

ACX_BUILD_ENVIRONMENT

dnl Delay the cross-compilation check before we add RPATH flags for ICON.
ACX_COMPILER_CROSS_CHECK_DELAY

AC_LANG([C])
dnl Override the default flags set by AC_PROG_CC:
CFLAGS=${CFLAGS-}
dnl Flags to be passed to the configure scripts of the bundled libraries:
subdir_CFLAGS=$CFLAGS
AC_PROG_CC
dnl Supporting ISO C89 means that the following header files are available:
dnl stdlib.h, stdio.h, math.h, time.h, ctype.h, string.h, signal.h, limits.h
dnl errno.h, ctype.h assert.h.
AS_VAR_IF([ac_cv_prog_cc_c89], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C89])])

dnl Some of the C source files require support for ISO C99. Supporting ISO C99
dnl also means that the following header files are available:
dnl stdbool.h.
AC_PROG_CC_C99
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])

ACX_COMPILER_CC_VERSION

AS_VAR_IF([enable_rpaths], [yes], [ACX_SHLIB_CC_RPATH_FLAG])

dnl Check for headers that are protected with HAVE_* macros in the code
dnl (override AC_INCLUDES_DEFAULT of autoconf 2.69, which runs additional
dnl annoying checks; already fixed in the autoconf master branch).
AC_CHECK_HEADERS([execinfo.h ucontext.h link.h], [], [],
  [[// override default includes]])

dnl Check for headers that are not protected at all.
AC_CHECK_HEADERS([unistd.h fcntl.h pwd.h stdint.h sys/types.h sys/stat.h]dnl
  [sys/param.h sys/resource.h sys/utsname.h sys/time.h sys/times.h dirent.h],
  [], [AC_MSG_FAILURE([unable to find header file $ac_header])],
  [[// override default includes]])

dnl Other header files that are protected in the code but not with the HAVE_*
dnl macros:
dnl support/util_hash.c: endian.h
dnl support/util_sysinfo.c: mcheck.h sys/systemcfg.h
dnl support/util_timer.c: sys/systemcfg.h asm/msr.h

dnl Header files that are included from non-mandatory files
dnl support/util_signal.c: features.h strings.h fenv.h fptrap.h

ICON_OPTIMIZATION_SET_CFLAGS
ICON_OPTIMIZATION_CHECK_CFLAGS
AC_SUBST([ICON_OPTIM_CFLAGS], [$icon_optim_CFLAGS])
subdir_CFLAGS="$icon_optim_subdir_CFLAGS $subdir_CFLAGS"

AC_LANG([Fortran])
dnl Override the default flags set by AC_PROG_FC:
FCFLAGS=${FCFLAGS-}
dnl Flags to be passed to the configure scripts of the bundled libraries:
subdir_FCFLAGS=$FCFLAGS
AC_PROG_FC
ACX_FC_PP_SRCEXT([f90])
ACX_COMPILER_FC_VERSION

dnl Notes on RPATHs for the bundled libraries:
dnl
dnl 1. RPATH flags must enable the dynamic linker to locate all the libraries
dnl    listed in the LIBS at least when running an executable linked with the
dnl    primary compiler, i.e. the compiler that is checked first by the
dnl    configure script of a bundled library: if AC_PROG_CC is expanded first
dnl    then the primary compiler is the C compiler, if AC_PROG_FC is expanded
dnl    first, then the primary compiler is the Fortran compiler. The primary
dnl    compiler for a package that uses libtool is always the C compiler. This
dnl    requirement is induced by the standard cross-compilation check, which,
dnl    if unsuccessful, fails the configuration. Obviously, this requirement
dnl    must be satisfied for the rest of the compilers if they are used at
dnl    configure-time in non-optional checks for linking with a following run
dnl    of the resulting executable.
dnl
dnl 2. RPATH flags must not break at least the compilation with the secondary
dnl    or any other compiler required for building. Obviously, posing such a
dnl    loose constraint is possible only if we can disable the usage of the
dnl    secondary compiler for linking. For example, if the configure script of
dnl    a bundled library does not fail when the linking with the secondary
dnl    compiler is impossible and allows for generation of a makefile that does
dnl    not create shared libraries or executables (e.g. the configure script
dnl    has the option '--disable-shared'). In simple cases, this requirement
dnl    can be met by moving RPATH flags from LDFLAGS to the list used by the
dnl    primary compiler only (e.g. CFLAGS).
dnl
dnl 3. In some cases, RPATH flags for the secondary compiler do not necessarily
dnl    have to enable the dynamic linker to locate the libraries when running
dnl    an executable but they must not break the configure-time checks that
dnl    require linking with the secondary compiler. This means that sometimes
dnl    it is better not to set any RPATH flags for the secondary compiler at
dnl    all.
dnl
dnl 4. RPATH flags must be correctly interpreted by libtool. Especially if
dnl    libtool is used for linking by the configure script of a bundled
dnl    library. First, this means that such script is most likely unable to
dnl    process RPATH flags correctly if their prefixes differ from
dnl    '-Wl,-rpath -Wl,' or '-Wl,-rpath,'. For example, we cannot run such
dnl    configure script with FC='nagfor', because
dnl    LDFLAGS='-Wl,-Wl,,-rpath -Wl,-Wl,,<lib dir>', which are valid for the
dnl    compiler itself, are incorrectly transformed by libtool into
dnl    '-Wl,-Wl -Wl,"" -Wl,-rpath -Wl,-Wl -Wl,"" -Wl,<lib dir>'. Moreover,
dnl    NAG compiler doesn't get special treatment by libtool if it is used as
dnl    C compiler (i.e. CC='nagfor =C'), because ac_cv_c_compiler_gnu equals
dnl    'yes' in this case (see the implementation of libtool). The second
dnl    problem is that libtool relies on the name of the compiler executable.
dnl    This means that depending on whether we are using the MPI compiler
dnl    wrappers or not, we need to specify RPATH flags differently (unless the
dnl    bundled library has a patched version of libtool that can recognize
dnl    compiler even if it is called with a wrapper).
dnl
dnl 5. The configure script of a bundled library might have its non-standard
dnl    solutions for the problems above.
dnl
dnl Currently, the configure scripts of the bundled libraries need the
dnl following treatment regarding RPATH flags:
dnl
dnl 1. RTE+RRTMGP: does not use libtool, needs only Fortran compiler. RPATH
dnl    flags can be appended to LDFLAGS with the prefixes accepted by the
dnl    Fortran compiler.
dnl
dnl 2. ECRAD: does not use libtool, needs only Fortran compiler. RPATH flags
dnl    can be appended to LDFLAGS with the prefixes accepted by the Fortran
dnl    compiler.
dnl
dnl 3. SCT: libtool-based, no configure-time checks using libtool directly
dnl    (however, the configure script tries to find .la files and parse them),
dnl    runs link-and-run checks for C compiler and link-only checks for Fortran
dnl    compilers. RPATH flags can be appended to LDFLAGS if the linker prefixes
dnl    accepted by Fortran and C compilers are the same. It is possible to
dnl    build the static version of SCT library (--disable-shared) if the linker
dnl    prefix of the C compiler is different from the prefix of libtool (e.g.
dnl    CC='nagfor =C') but the linking of the post-build tests (i.e. tests from
dnl    the test suite that are run with 'make check' command) fail in this
dnl    case. If the linker prefixes accepted by C and Fortran compilers differ,
dnl    RPATH flags should be appended to CFLAGS only to give libtool a chance
dnl    to link the post-build tests. Obviously, in this case, post-build tests
dnl    that run executables linked with Fortran compiler might not get all
dnl    required RPATHS (some of them will be added by libtool though) and fail.
dnl
dnl 4. YAXT: libtool-based, runs link-and-run checks for C compiler directly
dnl    and link-and-run checks for both C and Fortran compilers using libtool.
dnl    Implements custom workarounds that transform FCFLAGS, which we can't
dnl    adopt since they don't seem to work as expected (at least with Libtool
dnl    2.4.6). Therefore, we should not try to modify FCFLAGS. RPATH flags can
dnl    be appended to LDFLAGS with the prefixes accepted by libtool (YAXT has a
dnl    patched version of libtool, which recognizes NAG compiler when it is
dnl    called with an MPI wrapper). However, this solution works only if the C
dnl    compiler accepts RPATH flags in the same form (i.e. CC='nagfor =C' is
dnl    not supported).
dnl
dnl 5. CDI: libtool-based, runs link-and-run checks for C compiler directly,
dnl    link-only checks for C compiler using libtool, and compile-only checks
dnl    for Fortran compiler using libtool (if MPI support is disabled in CDI).
dnl    RPATH flags can be appended to LDFLAGS with the linker prefixes accepted
dnl    by libtool (CDI has a patched version of libtool, which recognizes NAG
dnl    compiler when it is called with an MPI wrapper). However, this solution
dnl    works only if the C compiler accepts RPATH flags in the same form (i.e.
dnl    CC='nagfor =C' is not supported).
dnl
dnl 6. MTIME: libtool-based, runs link-and-run checks for C compiler directly,
dnl    no checks using libtool, no link checks for Fortran compiler. It is
dnl    possible to disable any linking at build-time (--disable-examples,
dnl    --disable-shared), therefore RPATH flags can be appended to LDFLAGS with
dnl    the prefixes of the C compiler even if they are different from the
dnl    prefixes accepted by libtool (e.g. CC='nagfor =C'). However, the linking
dnl    of post-build tests will fail in this case.
dnl
dnl 7. YAC: does not use libtool, runs link-and-run checks for both C and
dnl    Fortran compilers. If the original forms of RPATH flags for C and
dnl    Fortran compilers are the same, we set LDFLAGS, otherwise we append the
dnl    flags to CFLAGS and FCFLAGS.
dnl
dnl 8. TIXI: does not use libtool, does not run linking checks for Fortran
dnl    compiler. RPATH flags can be appended to LDFLAGS wit the prefixes
dnl    accepted by the C compiler.
dnl
dnl 9. SELF: does not use libtool, needs only Fortran compiler. RPATH flags
dnl    can be appended to LDFLAGS with the prefixes accepted by the Fortran
dnl    compiler.
dnl
dnl Another possible solution for enabling the building and testing of the
dnl bundled libraries is to rely on the fact that we do not install the bundled
dnl libraries and link them statically, therefore all we probably need is to
dnl export and extend LD_LIBRARY_PATH (or its analogue) here before running
dnl the configure scripts and in the makefile before running the tests (the
dnl latter can be done either by extending the BUILD_ENV variable or by
dnl introducing a separate TESTS_ENVIRONMENT variable like the one supported
dnl by Automake, see
dnl https://www.gnu.org/software/automake/manual/html_node/Scripts_002dbased-Testsuites.html).
dnl
dnl The detection of the environment variable that the dynamic linker refers to
dnl when searching for the shared libraries (e.g. LD_LIBRARY_PATH) is already
dnl implemented in ./m4/acx_shlib.m4 as an M4 macro ACX_SHLIB_PATH_VAR, which
dnl can be used in the following way:
dnl
dnl   ACX_SHLIB_PATH_VAR
dnl   for icon_libdir in $icon_libdirs; do
dnl     ASX_VAR_APPEND_UNIQ([$acx_cv_shlib_path_var],
dnl       ["$icon_libdir"], ["$PATH_SEPARATOR"])
dnl   done
dnl   eval export $acx_cv_shlib_path_var
dnl
dnl The drawback of this approach is that we modify the environment, which a
dnl user doesn't expect and might not want because the modification might
dnl affect not only the dynamic linker at runtime but the linker at build-time
dnl too.

icon_c_rpaths= ; icon_fc_rpaths= ; icon_lt_rpaths=
AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([icon_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   icon_libdirs=
   for icon_L_arg in $icon_L_args; do
     if test -d "$icon_L_arg"; then
       ASX_VAR_APPEND_UNIQ([icon_libdirs], ["$icon_L_arg"], [' '])
     fi
   done
   for icon_libdir in $icon_libdirs; do
     AS_VAR_APPEND([icon_c_rpaths], ["$acx_cv_c_rpath_flag$icon_libdir "])
     AS_VAR_APPEND([icon_fc_rpaths], ["$acx_cv_fc_rpath_flag$icon_libdir "])
     AS_VAR_APPEND([icon_lt_rpaths], ["-Wl,-rpath -Wl,$icon_libdir "])
   done
   AC_LANG_PUSH([C])
   ACX_SHLIB_RPATH_FLAGS_CHECK([$icon_c_rpaths])
   AC_LANG_POP([C])
   ACX_SHLIB_RPATH_FLAGS_CHECK([$icon_fc_rpaths])])

dnl Original LDFLAGS to be passed to the configure scripts in subdirectories:
subdir_LDFLAGS=$LDFLAGS

dnl Add RPATH flags:
LDFLAGS="$icon_fc_rpaths$LDFLAGS"

dnl Run the delayed cross-compilation check:
ACX_COMPILER_CROSS_CHECK_NOW

ACX_FC_LINE_LENGTH([unlimited])

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

ACX_FC_ATTRIBUTE_CONTIGUOUS(
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_FC_ATTRIBUTE_CONTIGUOUS"])],
  [:])

ICON_OPTIMIZATION_SET_FCFLAGS
ICON_OPTIMIZATION_CHECK_FCFLAGS
AC_SUBST([ICON_OPTIM_FCFLAGS], [$icon_optim_FCFLAGS])
AC_SUBST([ICON_OPTIM_OCEAN_FCFLAGS], [$icon_optim_ocean_FCFLAGS])
subdir_FCFLAGS="$icon_optim_subdir_FCFLAGS $subdir_FCFLAGS"

AC_SUBST([DEFAULT_VERBOSITY], [0])
AS_VAR_IF([enable_silent_rules], [no], [DEFAULT_VERBOSITY=1])

AS_VAR_APPEND([FCFLAGS],
  [" ${FC_PP_DEF}__ICON__ ${FC_PP_DEF}__NO_ICON_TESTBED__"])

AS_VAR_IF([enable_atmo], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_ATMO__"])])

AS_VAR_IF([enable_ocean], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_OCEAN__"])])

AS_VAR_IF([enable_jsbach], [yes],
  [icon_subdir_srcdir="$srcdir/externals/jsbach"
   icon_subdir_unique_file="src/base/mo_jsb_base.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the land component JSBACH dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable the component dnl
(--disable-jsbach)])])
   AC_SUBST([JSBACH_ENABLED], [''])],
  [AC_SUBST([JSBACH_ENABLED], ['#'])
   AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_JSBACH__"])])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}YAC_coupling"])])

AS_VAR_IF([enable_psrad], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_PS_RAD__"])])

AS_VAR_IF([enable_ecrad], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__ECRAD"])])

dnl Currently, the macro __NO_RTE_RRTMGP__ is not used anywhere in the code but
dnl is set in the configure script of icon-2.6.0-rc:
AS_VAR_IF([enable_rte_rrtmgp], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_RTE_RRTMGP__"])])

AS_VAR_IF([enable_rttov], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__USE_RTTOV"])])

AS_VAR_IF([enable_art], [yes],
  [icon_subdir_srcdir="$srcdir/src/art"
   icon_subdir_unique_file="chem_init/mo_art_chem_init_utils.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the ART component dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable the component dnl
(--disable-art)])])
   AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}__ICON_ART ${FC_PP_DEF}__ART_TMP_IFDEF"])])

AS_VAR_IF([enable_mpi], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NOMPI"])])

AS_VAR_IF([enable_active_target_sync], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}HAVE_SLOW_PASSIVE_TARGET_ONESIDED"])])

AS_VAR_IF([enable_parallel_netcdf], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_PARALLEL_NETCDF"])])

AS_VAR_IF([enable_sct], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__SCT__"])])

AS_VAR_IF([enable_yaxt], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_YAXT"])])

AS_IF([test "x$enable_serialization" != xno],
  [AC_SUBST([SB2_ENABLED], [''])
   AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}SERIALIZE"])
   AS_CASE(["$enable_serialization"],
     [read], [AS_VAR_APPEND([FCFLAGS],
                [" ${FC_PP_DEF}SERIALIZE_READ_REFERENCE"])],
     [perturb], [AS_VAR_APPEND([FCFLAGS],
                   [" ${FC_PP_DEF}SERIALIZE_PERTURB_REFERENCE"])],
     [create], [AS_VAR_APPEND([FCFLAGS],
                  [" ${FC_PP_DEF}SERIALIZE_CREATE_REFERENCE"])])],
  [AC_SUBST([SB2_ENABLED], ['#'])])

AS_VAR_IF([enable_loop_exchange], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__LOOP_EXCHANGE"])])

AS_IF([test "x$enable_intel_consistency" = xauto && dnl
test "x$acx_cv_fc_compiler_vendor" = xintel],
  [enable_intel_consistency=yes])
AS_VAR_IF([enable_intel_consistency], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}IFORT_CONSISTENCY_ENFORCE"])])

AS_VAR_IF([enable_vectorized_lrtm], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}LRTM_FULL_VECTORIZATION"])])

AS_VAR_IF([enable_mixed_precision], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}__MIXED_PRECISION ${FC_PP_DEF}__MIXED_PRECISION_2"])])

dnl One of the topological orderings of the ICON dependency graph:
dnl
dnl '-lself', '-licon_tixi', '-lyac', '-llapack -lblas (or other)',
dnl '-lmtime', '-lSerialboxFortranShared', '-lcdi_f2003 -lcdi',
dnl '-lyaxt -lyaxt_c', '-lsct', '-lradiance -lrttov10.2 (etc.)',
dnl '-lradiation -lifsrrtm -lutilities -lifsaux', '-lrrtmgp -lrte',
dnl '-lnetcdff', '-lgomp (or other)', '-lmpifort -lmpi (or other)'
dnl
dnl We check these dependencies in the reverse order.

dnl A list of library files (i.e. *.a files) of the enabled bundled libraries:
AC_SUBST([BUNDLED_LIBS], [''])
dnl A list of Fortran compiler flags specifying search paths of the module
dnl files of the bundled libraries:
AC_SUBST([BUNDLED_FCFLAGS], [''])

dnl Additional macro definitions to aid the dependency generator,
dnl e.g. compiler intrinsic macros:
AC_SUBST([DEPGEN_FCFLAGS], [''])

dnl A list of Fortran modules that are used in ICON but are not part of it:
dnl MPI:
AC_SUBST([DEPGEN_EXTERNAL_MODS], ['mpi'])
dnl RTE+RRTMGP:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS],
  [',mo_gas_optics_rrtmgp,mo_gas_concentrations,mo_optical_props'])
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS],
  [',mo_source_functions,mo_rte_kind,mo_fluxes,mo_rte_lw,mo_rte_sw'])
dnl ECRAD:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS],
  [',radiation_cloud,radiation_interface,radiation_flux,radiation_aerosol'])
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS],
  [',radiation_gas,radiation_thermodynamics,radiation_config'])
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',radiation_single_level'])
dnl RTTOV:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',mo_rttov_ifc'])
dnl SCT:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',sct'])
dnl YAXT:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',yaxt'])
dnl CDI:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',mo_cdi'])
dnl SERIALBOX2
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',m_serialize,utils_ppser'])
dnl MTIME:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',mtime,mtime_datetime'])
dnl YAC:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',mo_yac_finterface'])
dnl TIXI:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',tixi'])
dnl SELF:
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS],
  [',self_assert,self_map,self_object,self_vector,self_vector_ref'])

dnl A list of Fortran intrinsic modules:
AC_SUBST([DEPGEN_INTRINSIC_MODS],
  [[iso_c_binding,iso_fortran_env,ieee_exceptions,ieee_arithmetic,]]dnl
[[ieee_features,omp_lib,omp_lib_kinds,openacc]])

AS_VAR_IF([enable_mpi], [yes],
  [AC_MSG_CHECKING([whether Fortran 90 interface of MPI library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
       use mpi, only: mpi_init
       implicit none
       integer :: error
       call mpi_init(error)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using Fortran 90 dnl
interface of MPI: either disable MPI support (--disable-mpi) or set FC to dnl
the MPI Fortran compiler wrapper, or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])])

AS_VAR_IF([enable_openmp], [yes],
  [ACX_LANG_OPENMP_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openmp_flag"])
      ACX_LANG_MACRO_CHECK_VALUE([_OPENMP],
         [201511 201307 201107 200805 200505 200011 199911 199710])
      AS_VAR_IF([acx_macro_value], [unknown],
        [AC_MSG_WARN([unable to detect supported OpenMP version])
         AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENMP"])],
        [AS_IF([test "$acx_macro_value" -lt 201307 >&AS_MESSAGE_LOG_FD 2>&1],
           [AC_MSG_WARN([Fortran compiler does not support OpenMP 4.0 dnl
(201307): it is strongly recommended to disable OpenMP support dnl
(--disable-openmp)])])
         AS_VAR_APPEND([DEPGEN_FCFLAGS],
           [" ${FC_PP_DEF}_OPENMP=$acx_macro_value"])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenMP support by default: set dnl
FCFLAGS accordingly to disable it])])])

AS_VAR_IF([enable_openacc], [yes],
  [ACX_LANG_OPENACC_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openacc_flag"])
      ACX_LANG_MACRO_CHECK_VALUE([_OPENACC],
         [201711 201510 201308 201306 201111])
      AS_VAR_IF([acx_macro_value], [unknown],
        [AC_MSG_WARN([unable to detect supported OpenACC version])
         AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENACC"])],
        [AS_IF([test "$acx_macro_value" -lt 201711 >&AS_MESSAGE_LOG_FD 2>&1],
           [AC_MSG_WARN([Fortran compiler does not support OpenACC 2.6 dnl
(201711): it is strongly recommended to disable OpenACC support dnl
(--disable-openacc)])])
      AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FC_PP_DEF}_OPENACC=$acx_macro_value"])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenACC support by default: dnl
set FCFLAGS accordingly to disable it])])])

AC_MSG_CHECKING([whether the Fortran 77 interface of NetCDF library is dnl
available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      include 'netcdf.inc'
      integer :: status, ncid
      status = nf_open('conftest.nc', nf_nowrite, ncid)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using Fortran 77 interface dnl
of NetCDF library: set FCFLAGS, LDFLAGS and LIBS accordingly])])
AS_VAR_IF([enable_parallel_netcdf], [yes],
  [AC_MSG_CHECKING(
     [whether the provided version of NetCDF is built with parallel support])
   icon_netcdf_supports_parallel=unknown
   AS_VAR_IF([cross_compiling], [no],
     [AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mpi, only: mpi_init, MPI_COMM_WORLD, &
                   & MPI_INFO_NULL
      implicit none
      include 'netcdf.inc'
      integer :: result, tmp
      call mpi_init(tmp)
      result = nf_create_par('conftest.nc', &
                           & NF_MPIIO, &
                           & MPI_COMM_WORLD, &
                           & MPI_INFO_NULL, &
                           & tmp)
      if (result == NF_ENOPAR) then
        write(*,"(a)") "conftest_parallel_no"
      else
        write(*,"(a)") "conftest_parallel_yes"
      endif]])],
        [icon_exec_result=`./conftest$ac_exeext 2>/dev/null`
         rm -f conftest.nc
         AS_IF([test $? -eq 0],
           [AS_IF(
              [echo "$icon_exec_result" | dnl
grep 'conftest_parallel_no' >/dev/null 2>&1],
              [icon_netcdf_supports_parallel=no],
              [echo "$icon_exec_result" | dnl
grep 'conftest_parallel_yes' >/dev/null 2>&1],
              [icon_netcdf_supports_parallel=yes])])])])
   AC_MSG_RESULT([$icon_netcdf_supports_parallel])
   AS_IF([test x"$icon_netcdf_supports_parallel" = xno],
        [AC_MSG_FAILURE([the provided version of NetCDF does not support dnl
parallel features: either disable usage of the parallel features of NetCDF dnl
(--disable-parallel-netcdf) or set FCFLAGS, LDFLAGS and LIBS accordingly])],
        [test x"$icon_netcdf_supports_parallel" = xunknown],
        [AC_MSG_WARN([unable to check whether the provided version of dnl
NetCDF supports parallel features])])])

AS_VAR_IF([enable_rte_rrtmgp], [yes],
  [AS_VAR_IF([with_external_rte_rrtmgp], [yes],
     [AC_MSG_CHECKING([whether Fortran RTE+RRTMGP library is available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_optical_props, only: ty_optical_props_1scl, &
                                & ty_optical_props_2str
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mo_source_functions, only: ty_source_func_lw
      use mo_rte_kind, only: wp
      use mo_fluxes, only: ty_fluxes_broadband
      use mo_rte_lw, only: rte_lw
      use mo_rte_sw, only: rte_sw
      use mo_gas_concentrations, only: ty_gas_concs
      implicit none
      type(ty_optical_props_1scl) :: optical_props_lw
      type(ty_optical_props_2str) :: optical_props_sw
      logical :: top_at_1
      type(ty_source_func_lw) :: source
      real(wp), dimension(:,:), allocatable :: sfc_emis_spec, &
                                             & toa_flux, &
                                             & sfc_alb_spec
      type(ty_fluxes_broadband) :: fluxes
      real(wp), dimension(:), allocatable :: mu0
      type(ty_gas_concs), dimension(:), allocatable :: gas_conc_array
      real(wp), dimension(:,:,:), allocatable :: gas_conc_temp_3d
      character(len=128) :: error_msg

      error_msg = rte_lw(optical_props_lw, top_at_1, source, &
                       & sfc_emis_spec, fluxes, n_gauss_angles = 1)
      error_msg = rte_sw(optical_props_sw, top_at_1, mu0, &
                       & toa_flux, sfc_alb_spec, sfc_alb_spec, fluxes)
      error_msg = gas_conc_array(1)%set_vmr('h2o',gas_conc_temp_3d(:,:,1))]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using RTE+RRTMGP: dnl
either disable RTE+RRTMGP (--disable-rte-rrtmgp) or use the bundled version dnl
of RTE+RRTMGP (--without-external-rte-rrtmgp), or set FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])],
     [icon_subdir_srcdir="$srcdir/externals/rte-rrtmgp"
      icon_subdir_unique_file="rte/mo_rte_sw.F90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the the RTE+RRTMGP toolbox dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable dnl
RTE+RRTMGP (--disable-rte-rrtmgp) or provide an external version of dnl
RTE+RRTMGP (--with-external-rte-rrtmgp) by setting FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])
      ACX_CONFIG_SUBDIR([externals/rte-rrtmgp],
        [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
         [NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([openacc])]],
        [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
         [LDFLAGS="$icon_fc_rpaths$subdir_LDFLAGS"],
         [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
         [--enable-openacc="$enable_openacc"],
         [--enable-silent-rules="$enable_silent_rules"]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS],
        [" ${FC_MOD_IN}externals/rte-rrtmgp/mod"])
      BUNDLED_LIBS="externals/rte-rrtmgp/librte.a $BUNDLED_LIBS"
      BUNDLED_LIBS="externals/rte-rrtmgp/librrtmgp.a $BUNDLED_LIBS"])])

AS_VAR_IF([enable_ecrad], [yes],
  [AS_VAR_IF([with_external_ecrad], [yes],
     [AC_MSG_CHECKING([whether Fortran ECRAD library is available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use radiation_config, only: config_type
      use radiation_interface, only: setup_radiation
      implicit none
      type(config_type) :: conf
      call setup_radiation(conf)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using ECRAD: either dnl
disable ECRAD (--disable-ecrad) or use the bundled version of ECRAD dnl
(--without-external-ecrad), or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/ecrad"
      icon_subdir_unique_file="radiation/radiation_flux.F90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the the ECMWF radiation dnl
scheme ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable dnl
ECRAD (--disable-ecrad) or provide an external version of ECRAD dnl
(--with-external-ecrad) by setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      ACX_CONFIG_SUBDIR([externals/ecrad],
        [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
         [NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([single-precision])]],
        [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
         [LDFLAGS="$icon_fc_rpaths$subdir_LDFLAGS"],
         [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
         [--enable-silent-rules="$enable_silent_rules"]])
   AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/ecrad/mod"])
   BUNDLED_LIBS="externals/ecrad/libifsaux.a $BUNDLED_LIBS"
   BUNDLED_LIBS="externals/ecrad/libutilities.a $BUNDLED_LIBS"
   BUNDLED_LIBS="externals/ecrad/libifsrrtm.a $BUNDLED_LIBS"
   BUNDLED_LIBS="externals/ecrad/libradiation.a $BUNDLED_LIBS"])])

AS_VAR_IF([enable_rttov], [yes],
  [AC_MSG_CHECKING([whether Fortran RTTOV library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_rttov_ifc, only: rttov_init
      implicit none
      integer :: istatus, channels(8, 2), n_chans(2), p_pe, &
               & num_work_procs, p_io, p_comm_work
      integer, allocatable :: instruments(:,:)
      logical, allocatable :: addclouds(:)
      istatus = rttov_init(  &
        instruments     , &
        channels        , &
        n_chans         , &
        p_pe            , &
        num_work_procs  , &
        p_io            , &
        p_comm_work     , &
        appRegLim=.TRUE., &
        readCloud=addclouds)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using RTTOV: either dnl
disable usage of RTTOV (--disable-rttov) or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])])

AS_VAR_IF([enable_sct], [yes],
  [AS_VAR_IF([with_external_sct], [yes],
     [AC_MSG_CHECKING([whether the Fortran interface of SCT library is dnl
available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use sct, only: sct_init
      implicit none
      call sct_init(timer_max=512)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using SCT: either dnl
disable SCT timer (--disable-sct) or use the bundled version of SCT dnl
(--without-external-sct), or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/sct"
      icon_subdir_unique_file="src/sct_reporter.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
SCT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the SCT timer (--disable-sct) or provide an external dnl
version of SCT (--with-external-sct) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      sct_CFLAGS=$subdir_CFLAGS
      sct_LDFLAGS=$subdir_LDFLAGS
      AS_IF([test -n "$icon_c_rpaths"],
        [sct_test_failure_stage=
         AS_IF([test x"$icon_c_rpaths" = x"$icon_fc_rpaths"],
           [sct_LDFLAGS="$icon_c_rpaths$sct_LDFLAGS"],
           [AS_VAR_APPEND([sct_CFLAGS], [" $icon_c_rpaths"])
            sct_test_failure_stage=run])
         AS_IF([test x"$icon_c_rpaths" != x"$icon_lt_rpaths"],
           [sct_test_failure_stage=link])
         AS_IF([test -n "$sct_test_failure_stage"],
           [AC_MSG_WARN([unable to set valid RPATH flags for the bundled dnl
version of SCT: the test suite of SCT might fail at $sct_test_failure_stage dnl
time])])])
      ACX_CONFIG_SUBDIR([externals/sct],
        [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*|MPICC=*|MPIFC=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpi])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpirun-check])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([openmp])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([rpath])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([papi])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([hdf5])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([libpapi-prefix])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([libhdf5-prefix])]],
        [[FCFLAGS="$subdir_FCFLAGS"], [CFLAGS="$sct_CFLAGS"],
         [LDFLAGS="$sct_LDFLAGS"], [MPICC=], [MPIFC=], [--enable-static],
         [--disable-shared], [--disable-maintainer-mode],
         [--enable-mpi="$enable_mpi"], [--enable-openmp="$enable_openmp"],
         [--disable-rpath], [--disable-papi], [--enable-hdf5],
         [--enable-silent-rules="$enable_silent_rules"],
         [--without-libpapi-prefix], [--without-libhdf5-prefix]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/sct/src"])
      BUNDLED_LIBS="externals/sct/src/.libs/libsct.a $BUNDLED_LIBS"])])

AS_VAR_IF([enable_yaxt], [yes],
  [AS_VAR_IF([with_external_yaxt], [yes],
     [AC_MSG_CHECKING([whether the Fortran interface of YAXT library dnl
is available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use yaxt, only: xt_initialized, xt_initialize
      implicit none
      integer :: a
      if (.not. xt_initialized()) call xt_initialize(a)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using YAXT: either dnl
disable the YAXT data exchange (--disable-yaxt) or use the bundled version dnl
of YAXT (--without-external-yaxt), or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/yaxt"
      icon_subdir_unique_file="src/xt_idxlist.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
YAXT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the YAXT data exchange (--disable-yaxt) or provide an external dnl
version of YAXT (--with-external-yaxt) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      yaxt_LDFLAGS=$subdir_LDFLAGS
      AS_IF([test -n "$icon_lt_rpaths"],
        [AS_IF([test x"$icon_c_rpaths" = x"$icon_lt_rpaths"],
           [yaxt_LDFLAGS="$icon_lt_rpaths$yaxt_LDFLAGS"],
           [AC_MSG_WARN([unable to set valid RPATH flags for the bundled dnl
version of YAXT: the configure script of YAXT might fail])])])
      ACX_CONFIG_SUBDIR([externals/yaxt],
        [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([dependency-tracking])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([regard-for-quality])]],
        [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
         [LDFLAGS="$yaxt_LDFLAGS"], [--enable-static],
         [--disable-shared], [--enable-dependency-tracking],
         [--enable-silent-rules="$enable_silent_rules"],
         [--with-regard-for-quality="$enable_mpi_checks"]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/yaxt/src"])
      BUNDLED_LIBS="externals/yaxt/src/.libs/libyaxt.a dnl
externals/yaxt/src/.libs/libyaxt_c.a $BUNDLED_LIBS"])])

AS_VAR_IF([with_external_cdi], [yes],
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
      AC_MSG_FAILURE([unable to link a test program using CDI: either use dnl
the bundled version of CDI (--without-external-cdi) or set FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])
dnl
   AC_MSG_CHECKING([whether the provided version of CDI 1.8.x is patched])
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
dnl
   icon_cdi_formats="NetCDF-2:CDI_FILETYPE_NC2 NetCDF-4:CDI_FILETYPE_NC4 dnl
GRIB1:CDI_FILETYPE_GRB"
   AS_VAR_IF([enable_grib2], [yes],
     [AS_VAR_APPEND([icon_cdi_formats], [" GRIB2:CDI_FILETYPE_GRB2"])])
   for icon_cdi_format in $icon_cdi_formats; do
     icon_cdi_format_name=`echo "$icon_cdi_format" | cut -d: -f1`
     icon_cdi_format_ID=`echo "$icon_cdi_format" | cut -d: -f2`
     icon_cdi_format_result=unknown
     AC_MSG_CHECKING([whether the provided version of CDI is built with dnl
$icon_cdi_format_name support])
     AS_VAR_IF([cross_compiling], [no],
       [AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: cdiHaveFiletype, $icon_cdi_format_ID
      write(*, "(i0)") cdiHaveFiletype($icon_cdi_format_ID)]])],
           [icon_exec_result=`./conftest$ac_exeext 2>/dev/null`
            AS_IF([test $? -eq 0],
              [AS_IF([test $icon_exec_result -eq 1],
                 [icon_cdi_format_result=yes],
                 [icon_cdi_format_result=no])])])])
     AC_MSG_RESULT([$icon_cdi_format_result])
     AS_IF([test x"$icon_cdi_format_result" = xno],
        [AC_MSG_FAILURE([the provided version of CDI does not support dnl
$icon_cdi_format_name])],
        [test x"$icon_cdi_format_result" = xunknown],
        [AC_MSG_WARN([unable to check whether the provided version of CDI dnl
supports $icon_cdi_format_name])])
   done],
dnl
  [icon_subdir_srcdir="$srcdir/externals/cdi"
   icon_subdir_unique_file="src/cdi.h"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of CDI dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of CDI (--with-external-cdi) by setting FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])
   cdi_LDFLAGS=$subdir_LDFLAGS
   AS_IF([test -n "$icon_lt_rpaths"],
     [AS_IF([test x"$icon_c_rpaths" = x"$icon_lt_rpaths"],
        [cdi_LDFLAGS="$icon_lt_rpaths$cdi_LDFLAGS"],
        [AC_MSG_WARN([unable to set valid RPATH flags for the bundled dnl
version of CDI: the configure script of CDI might fail])])])
   ACX_CONFIG_SUBDIR([externals/cdi],
     [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([iso-c-interface])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([cf-interface])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([cdi-app])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([grib])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([cgribex])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpi])],
      [ACX_CONFIG_SUBDIR_PATTERN_WITH([grib_api])],
      [ACX_CONFIG_SUBDIR_PATTERN_WITH([eccodes])],
      [ACX_CONFIG_SUBDIR_PATTERN_WITH([threads])],
      [ACX_CONFIG_SUBDIR_PATTERN_WITH([netcdf])]],
     [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
      [LDFLAGS="$cdi_LDFLAGS"], [--enable-static],
      [--disable-shared], [--enable-silent-rules="$enable_silent_rules"],
      [--enable-iso-c-interface], [--disable-cf-interface],
      [--disable-cdi-app], [--enable-grib], [--enable-cgribex],
      [--disable-mpi],
dnl Currently, we don't have full control on choosing the GRIB2 backend (e.g.
dnl grib_api, eccodes) due to the limitations of the configure script of
dnl libcdi. Therefore, we set --with-grib_api=$enable_grib2 --without-eccodes
dnl unconditionally. This way, the configure script of libcdi enables support
dnl for grib2 regardless of whether grib_api or eccodes library is provided
dnl with the LDFLAGS and LIBS variables.
      [--with-grib_api="$enable_grib2"], [--without-eccodes],
      [--without-threads], [--with-netcdf],
dnl The configure script of CDI runs nc-config to check for the features of the
dnl provided NetCDF library. The script might not be available or might belong
dnl to another NetCDF installation. Therefore, we enforce enabling of NetCDF2
dnl and NetCDF4 features by setting the corresponding cache variables on the
dnl command line.
      [acx_cv_have_netcdf2=yes], [acx_cv_have_netcdf4=yes],
dnl Disable unwanted extra dependencies.
      [acx_cv_have_libnc_dap=no], [acx_cv_have_nc4hdf5=no]])
   AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/cdi/src"])
   BUNDLED_LIBS="externals/cdi/src/.libs/libcdi_f2003.a dnl
externals/cdi/src/.libs/libcdi.a $BUNDLED_LIBS"])

AS_IF([test "x$enable_serialization" != xno],
  [AC_MSG_CHECKING(
     [whether the Fortran interface of SERIALBOX2 library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use utils_ppser, only: ppser_set_mode
      implicit none
      call ppser_set_mode(1)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using SERIALBOX2: dnl
either disable the serialization (--disable-serialization) or set FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([with_external_mtime], [yes],
  [AC_MSG_CHECKING(
     [whether the Fortran interface of MTIME library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mtime, only: max_date_str_len, date, &
                     & newdate, dateToString
      implicit none
      character(len=max_date_str_len) :: test_date_string
      type(date), pointer :: test_date
      test_date => newdate(0,0,0)
      call dateToString(test_date, test_date_string)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using MTIME: either use dnl
the bundled version of MTIME (--without-external-mtime) or set FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])],
dnl
  [icon_subdir_srcdir="$srcdir/externals/mtime"
   icon_subdir_unique_file="include/mtime_calendar.h"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of MTIME dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of MTIME (--with-external-mtime) by setting CFLAGS, dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])
   mtime_LDFLAGS=$subdir_LDFLAGS
   AS_IF([test -n "$icon_c_rpaths"],
     [mtime_LDFLAGS="$icon_c_rpaths$mtime_LDFLAGS"
      AS_IF([test x"$icon_c_rpaths" != x"$icon_lt_rpaths"],
        [AC_MSG_WARN([unable to set valid RPATH flags for the bundled dnl
version of MTIME: the test suite of MTIME might fail at link time])])])
   ACX_CONFIG_SUBDIR([externals/mtime],
     [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([examples])]],
     [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
      [LDFLAGS="$mtime_LDFLAGS"], [--enable-static],
      [--disable-shared], [--enable-silent-rules="$enable_silent_rules"],
      [--disable-examples]])
   AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/mtime/src"])
   BUNDLED_LIBS="externals/mtime/src/.libs/libmtime.a $BUNDLED_LIBS"])

AC_MSG_CHECKING([whether Fortran BLAS library is available])
AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([dgemm])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using BLAS: set FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether Fortran LAPACK library is available])
AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([dgeev])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using LAPACK: set FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_IF([with_external_yac], [yes],
     [AC_MSG_CHECKING([whether the Fortran interface of YAC library is dnl
available])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_yac_finterface, only: yac_finit
      implicit none
      character(len=1024) :: a, b
      call yac_finit(a, b)]])],
        [AC_LINK_IFELSE([],
           [AC_MSG_RESULT([yes])],
           [AC_MSG_RESULT([no])
            AS_VAR_IF([with_external_mtime], [yes],
              [AC_MSG_FAILURE([unable to link a test program using YAC: dnl
either disable the coupling (--disable-coupling) or use the bundled version dnl
of YAC (--without-external-yac), or set LDFLAGS and LIBS accordingly])],
              [AC_MSG_FAILURE([unable to link a test program using YAC: dnl
you are trying to build with an external YAC and the bundled MTIME, so you dnl
need to make sure that the version of MTIME that YAC has been built with is dnl
available by setting LDFLAGS and LIBS accordingly, otherwise disable the dnl
coupling (--disable-coupling) or use the bundled version of YAC dnl
(--without-external-yac)])])])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to compile a test program using YAC: either dnl
disable the coupling (--disable-coupling) or use the bundled version of YAC dnl
(--without-external-yac), or set FCFLAGS accordingly])])],
     [icon_subdir_srcdir="$srcdir/externals/yac"
      icon_subdir_unique_file="src/grid.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of YAC dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the coupling (--disable-coupling) or provide an external version of dnl
YAC (--with-external-yac) by setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      yac_MTIME_CFLAGS= ; yac_MTIME_CLIBS=
      AS_VAR_IF([with_external_mtime], [no],
        [dnl Set $ac_abs_top_srcdir and $ac_abs_top_builddir:
         _AC_SRCDIRS([.])
         yac_MTIME_CFLAGS="-I${ac_abs_top_srcdir}/externals/mtime/include"
         yac_MTIME_CLIBS="${ac_abs_top_builddir}/externals/mtime/src/.libs/libmtime.a"])
      yac_CFLAGS=$subdir_CFLAGS
      yac_FCFLAGS=$subdir_FCFLAGS
      yac_LDFLAGS=$subdir_LDFLAGS
      AS_IF([test x"$icon_c_rpaths" = x"$icon_fc_rpaths"],
        [yac_LDFLAGS="$icon_c_rpaths$yac_LDFLAGS"],
        [AS_VAR_APPEND([yac_CFLAGS], [" $icon_c_rpaths"])
         AS_VAR_APPEND([yac_FCFLAGS], [" $icon_fc_rpaths"])])
      ACX_CONFIG_SUBDIR([externals/yac],
        [[FCFLAGS=*|CFLAGS=*|LDFLAGS=*],
         [MTIME_CFLAGS=*|MTIME_CLIBS=*|NETCDF_CFLAGS=*|NETCDF_CLIBS=*],
         [XML2_CFLAGS=*|XML2_CLIBS=*|MPI_CFLAGS=*|MPI_CLIBS=*],
         [MPI_FCFLAGS=*|MPI_FCLIBS=*|FORTRAN_LAPACK_CLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([lib-only])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([xml])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([netcdf])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpi])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpi-checks])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([external-lapack])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([external-mtime])]],
        [[FCFLAGS="$yac_FCFLAGS"], [CFLAGS="$yac_CFLAGS"],
         [LDFLAGS="$yac_LDFLAGS"], [MTIME_CFLAGS="$yac_MTIME_CFLAGS"],
         [MTIME_CLIBS=], [NETCDF_CFLAGS=], [NETCDF_CLIBS=],
         [XML2_CFLAGS=], [XML2_CLIBS=], [MPI_CFLAGS=], [MPI_CLIBS=],
         [MPI_FCFLAGS=], [MPI_FCLIBS=], [FORTRAN_LAPACK_CLIBS=],
         [--enable-silent-rules="$enable_silent_rules"], [--enable-lib-only],
         [--enable-xml], [--enable-netcdf], [--enable-mpi="$enable_mpi"],
         [--with-external-lapack=fortran], [--with-external-mtime],
         [--enable-mpi-checks="$enable_mpi_checks"],
         [acx_cv_c_lib_func_datetimeToString="$yac_MTIME_CLIBS"]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/yac/src"])
      BUNDLED_LIBS="externals/yac/src/libyac.a $BUNDLED_LIBS"])])

AS_VAR_IF([enable_art], [yes],
  [AS_VAR_IF([with_external_tixi], [yes],
     [AC_MSG_CHECKING([whether the Fortran interface of TIXI library is dnl
available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use tixi, only: tixiOpenDocument
      implicit none
      character(len=200) :: a
      integer :: b
      b = tixiOpenDocument(a, b)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using TIXI: either dnl
disable the ART component (--disable-art) or use the bundled version of dnl
TIXI (--without-external-tixi), or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])],
     [icon_subdir_srcdir="$srcdir/externals/tixi"
      icon_subdir_unique_file="src/tixi.f90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
TIXI library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the ART component (--disable-art) or provide an external version of dnl
TIXI (--with-external-tixi) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      ACX_CONFIG_SUBDIR([externals/tixi],
        [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
         [XML2_CFLAGS=*|XML2_CLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])]],
        [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
         [LDFLAGS="$icon_c_rpaths$subdir_LDFLAGS"], [XML2_CFLAGS=],
         [XML2_CLIBS=], [--enable-silent-rules="$enable_silent_rules"]])
      AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/tixi/src"])
      BUNDLED_LIBS="externals/tixi/src/libicon_tixi.a $BUNDLED_LIBS"])])

AS_VAR_IF([with_external_self], [yes],
  [AC_MSG_CHECKING([whether Fortran SELF 0.2 library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use self_object
      use self_vector_ref
      use self_vector
      use self_map, only: map
      use self_assert
      implicit none
      type(map) :: a
      call a%init()]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using SELF 0.2 dnl
library: either use the bundled version of SELF (--without-external-self) dnl
or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
  [icon_subdir_srcdir="$srcdir/externals/self"
   icon_subdir_unique_file="src/self_object.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of SELF dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of SELF (--with-external-self) by setting FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])
   ACX_CONFIG_SUBDIR([externals/self],
     [[CFLAGS=*|FCFLAGS=*|LDFLAGS=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])]],
     [[CFLAGS="$subdir_CFLAGS"], [FCFLAGS="$subdir_FCFLAGS"],
      [LDFLAGS="$icon_fc_rpaths$subdir_LDFLAGS"],
      [--enable-silent-rules="$enable_silent_rules"]])
   AS_VAR_APPEND([BUNDLED_FCFLAGS], [" ${FC_MOD_IN}externals/self/src"])
   BUNDLED_LIBS="externals/self/src/libself.a $BUNDLED_LIBS"])

AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python python3],
  [icon_prog_try=dnl
'$acx_candidate $srcdir/utils/mkhelper/depgen.py -h > /dev/null'
   AC_TRY_EVAL([icon_prog_try])],
  [PYTHON=$acx_candidate])

AC_ARG_VAR([PERL], [Perl interpreter command])dnl
ACX_PROG_SEARCH([PERL], [perl],
  [icon_prog_try='$acx_candidate $srcdir/utils/pvcs.pl --help > /dev/null'
   AC_TRY_EVAL([icon_prog_try])
   test $ac_status -eq 0 || test $ac_status -eq 1],
  [PERL=$acx_candidate])

AC_ARG_VAR([SB2PP], [Serialbox2 preprocessor command])
AC_MSG_CHECKING([for Serialbox2 preprocessor command])
AS_VAR_IF([enable_serialization], [no],
  [SB2PP=
   AC_MSG_RESULT([none needed])],
  [AC_LANG_ASSERT([Fortran])
   icon_sb2pp_result=unknown
   SB2PP=${SB2PP-'${PYTHON} -m pp_ser'}
dnl pp_ser.py can not process a code like this:
dnl   program main
dnl   end
dnl Therefore, we need to generate the source code manually:
   AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      program main
      implicit none
      !\$ser mode write
      call ppser_set_mode(1)
      end program]])])
   icon_prog_try=dnl
"$SB2PP -o conftest.$ac_ext conftest.$ac_ext >&AS_MESSAGE_LOG_FD"
   AC_TRY_EVAL([icon_prog_try])
   AS_IF([test $? -eq 0], [AC_LINK_IFELSE([], [icon_sb2pp_result=$SB2PP])])
   rm -f conftest.$ac_ext
   AC_MSG_RESULT([$icon_sb2pp_result])
   AS_VAR_IF([icon_sb2pp_result], [unknown],
     [AC_MSG_FAILURE([unable to find a valid Serialbox2 preprocessor: dnl
either disable the serialization (--disable-serialization) or make sure dnl
that path to the 'pp_ser' module is in the PYTHONPATH environment variable, dnl
or set SB2PP accordingly])])])

AC_PROG_INSTALL
AC_PROG_MKDIR_P
dnl The message generated by AC_PROG_LN_S is misleading in our case, therefore
dnl we do what the standard macro does but without the message:
AC_SUBST([LN_S], [$as_ln_s])
AC_PROG_SED

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  depgen.c.config
  depgen.f90.config
  icon.mk
])

dnl An additional step to generate dummy dependencies:
AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} -f icon.mk dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch icon.mk])

AC_OUTPUT
