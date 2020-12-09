AC_PREREQ([2.69])
AC_INIT([icon], [2.6.0], [], [icon-2.6.0])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([utils])
AC_CONFIG_SRCDIR([src/drivers/icon.f90])

dnl Temporary message for the users that is emitted on the first failure of the
dnl script. Please, remove it once all the current users are used to using the
dnl wrappers:
icon_usr_note="
********************************************************************************
Note that ICON building system has been updated. Now users are recommended to
run site-specific wrappers instead of this generic configure script.

You can find the wrappers in '${srcdir}/config/<organization>/'.

If you are running a wrapper and still see this message than either you are
using a wrong one or the wrapper is broken, or there is a problem with the new
configure script.
********************************************************************************
"
m4_define([AC_MSG_FAILURE],
[{ AS_MESSAGE([error: in `$ac_pwd':], 2)
AC_MSG_ERROR([$1
See `config.log' for more details
${icon_usr_note}], [$2]); }])
dnl End of the temporary message injection. Delete up until this line once the
dnl message is not needed anymore.

dnl We need to make sure that the source directory is clean when building
dnl out-of-source for two reasons:
dnl   1. Building systems of most of the bundled libraries are Automake-based
dnl      and, therefore, they rely on the standard non-selective search path
dnl      (i.e., VPATH) feature of make for out-of-source builds. This makes
dnl      having a clean source directory a hard requirement for them. Otherwise,
dnl      instead of generating object files in the build directory, make would
dnl      simply take them from the dirty source directory. For that reason, the
dnl      respective configure scripts will run a check similar to the one below
dnl      and fail if their source directories are dirty. It is better that we
dnl      run the check much earlier and fail the configuration with a more
dnl      informative error message.
dnl   2. The makefiles of ICON rely on the selective search path (i.e. vpath)
dnl      feature of GNU Make. Therefore, the problem with the object files that
dnl      is described above is irrelevant. However, icon.mk supports Fortran
dnl      source file preprocessing. The input and output files of the most of
dnl      the preprocessing steps are valid Fortran files. Thus, both the
dnl      original files in the source directory and the preprocessed files in
dnl      the build directory have the same suffix '.f90', which makes it
dnl      impossible to declare selective search paths that would work for the
dnl      original files but wouldn't work for the preprocessed files. This might
dnl      lead to the preprocessed files being taken from the dirty source
dnl      directory instead of being generated in the build directory. A
dnl      workaround for this, which had been implemented before, is to refer to
dnl      the preprocessed files by their absolute paths. However, it makes the
dnl      makefiles overcomplicated and error-prone. Instead, we simply forbid
dnl      the source directory to be dirty in the case of out-of-source build.
AS_IF([test "`cd $srcdir && pwd`" != "`pwd`" && dnl
test -f "$srcdir/config.status"],
  [AC_MSG_ERROR([source directory already configured; dnl
run "make distclean" there first])])

############################# OPTION DECLARATIONS ##############################

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

AC_ARG_ENABLE([concurrent-psrad],
  [AC_HELP_STRING([--enable-concurrent-psrad],
     [enable the concurrent mode of the radiation component PSRAD
@<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_concurrent_psrad=yes])],
  [enable_concurrent_psrad=yes])

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

AC_ARG_ENABLE([dace],
  [AC_HELP_STRING([--enable-dace],
     [enable the DACE modules for data assimilation @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_dace=yes])],
  [enable_dace=no])

AC_ARG_ENABLE([emvorado],
  [AC_HELP_STRING([--enable-emvorado],
     [enable the radar forward operator EMVORADO @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_emvorado=yes])],
  [enable_emvorado=no])

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
     [enable OpenMP support @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_openmp=yes])],
  [enable_openmp=no])

AC_ARG_ENABLE([gpu],
  [AC_HELP_STRING([--enable-gpu],
     [enable GPU support @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_gpu=yes])],
  [enable_gpu=no])

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

AC_ARG_ENABLE([cdi-pio],
  [AC_HELP_STRING([--enable-cdi-pio],
     [enable usage of the parallel features of CDI @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_cdi_pio=yes])],
  [enable_cdi_pio=no])

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

AC_ARG_ENABLE([claw],
[  --enable-claw           enable CLAW preprocessing. The value of the argument
                          must be one of the following @<:@default=no@:>@:
                            std|yes   enable STANDARD mode
                            validate  enable VALIDATE mode
                            no        disable CLAW preprocessing],
[AS_CASE(["$enableval"],
   [yes], [enable_claw=std],
   [std|validate|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-claw='$enableval'; valid values are 'std', 'validate', 'yes' dnl
(same as 'std'), 'no'])])],
[enable_claw=no])

AC_ARG_ENABLE([explicit-fpp],
  [AC_HELP_STRING([--enable-explicit-fpp],
     [enable explicit Fortran preprocessing @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_explicit_fpp=yes])],
  [enable_explicit_fpp=auto])

AC_ARG_ENABLE([serialization],
[  --enable-serialization  enable the Serialbox2 serialization. The value of
                          the argument must be one of the following
                          @<:@default=no@:>@:
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

AC_ARG_ENABLE([testbed],
  [AC_HELP_STRING([--enable-testbed],
     [enable ICON Testbed infrastructure @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_testbed=yes])],
  [enable_testbed=no])

m4_divert_once([HELP_ENABLE], [[
Optimization Features:]])

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

AC_ARG_WITH([external-cub],
  [AC_HELP_STRING([--with-external-cub],
     [use external CUB library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_cub=yes])],
  [with_external_cub=no])

############################ MAKEFILE CONDITIONALS #############################

AS_VAR_IF([enable_silent_rules], [yes],
  [AC_SUBST([DEFAULT_VERBOSITY], [0])],
  [AC_SUBST([DEFAULT_VERBOSITY], [1])])

AS_VAR_IF([enable_jsbach], [yes],
  [icon_subdir_srcdir="$srcdir/externals/jsbach"
   icon_subdir_unique_file="src/base/mo_jsb_base.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the land component JSBACH dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable the component dnl
(--disable-jsbach)])])
   AC_SUBST([JSBACH_ENABLED], [''])],
  [AC_SUBST([JSBACH_ENABLED], ['#'])])

icon_dace_enabled_or_needed=no
AS_IF([test x"$enable_dace" = xyes || test x"$enable_emvorado" = xyes],
  [icon_dace_enabled_or_needed=yes])

AS_VAR_IF([icon_dace_enabled_or_needed], [yes],
  [icon_subdir_srcdir="$srcdir/externals/dace_icon"
   icon_subdir_unique_file="src_for_icon/mo_fdbk_tables.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the DACE modules for data dnl
assimilation ($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable dnl
the modules (--disable-dace) and the radar forward operator EMVORADO dnl
(--disable-emvorado)])])
   AC_SUBST([DACE_ENABLED_OR_NEEDED], [''])],
  [AC_SUBST([DACE_ENABLED_OR_NEEDED], ['#'])])

AS_VAR_IF([enable_emvorado], [yes],
  [icon_subdir_srcdir="$srcdir/externals/emvorado"
   icon_subdir_unique_file="src_emvorado/radar_organize.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the radar forward operator dnl
EMVORADO ($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable the dnl
operator (--disable-emvorado)])])
   AC_SUBST([EMVORADO_ENABLED], [''])],
  [AC_SUBST([EMVORADO_ENABLED], ['#'])])

AS_VAR_IF([enable_art], [yes],
  [icon_subdir_srcdir="$srcdir/externals/art"
   icon_subdir_unique_file="chem_init/mo_art_chem_init_utils.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the aerosols and reactive dnl
trace component ART ($icon_subdir_unique_file) in '$icon_subdir_srcdir': dnl
disable the component (--disable-art)])])
   AC_SUBST([ART_ENABLED], [''])],
  [AC_SUBST([ART_ENABLED], ['#'])])

AS_VAR_IF([enable_mpi], [yes],
  [AC_SUBST([MPI_ENABLED], [''])],
  [AC_SUBST([MPI_ENABLED], ['#'])])

AS_VAR_IF([enable_openmp], [yes],
  [AC_SUBST([OPENMP_ENABLED], [''])],
  [AC_SUBST([OPENMP_ENABLED], ['#'])])

AS_VAR_IF([enable_gpu], [yes],
  [AC_SUBST([GPU_ENABLED], [''])],
  [AC_SUBST([GPU_ENABLED], ['#'])])

AS_IF([test x"$enable_claw" != xno],
  [AC_SUBST([CLAW_ENABLED], [''])],
  [AC_SUBST([CLAW_ENABLED], ['#'])])

AS_IF([test x"$enable_explicit_fpp" = xauto && dnl
test x"$enable_serialization" != xno],
  [AC_MSG_NOTICE([explicit Fortran preprocessing is enabled because the dnl
Serialbox2 serialization is requested: disable the explicit preprocessing dnl
if required (--disable-explicit-fpp)])
   enable_explicit_fpp=yes])
AS_VAR_IF([enable_explicit_fpp], [yes],
  [AC_SUBST([FPP_ENABLED], [''])],
  [AC_SUBST([FPP_ENABLED], ['#'])])

AS_IF([test x"$enable_serialization" != xno],
  [AC_SUBST([SB2_ENABLED], [''])],
  [AC_SUBST([SB2_ENABLED], ['#'])])

############################## ENVIRONMENT CHECKS ##############################

ACX_BUILD_ENVIRONMENT

AC_CANONICAL_HOST
ACX_HOST_FQDN([], [:])
AC_SUBST([host_fqdn], [$acx_cv_host_fqdn])

###################### FORTRAN COMPILER COMMAND AND FLAGS ######################
dnl Delay the cross-compilation check:
ACX_COMPILER_CROSS_CHECK_DELAY

dnl Override the default flags set by AC_PROG_FC:
FCFLAGS=${FCFLAGS-}
AC_PROG_FC
AC_LANG([Fortran])

ACX_FC_PP_SRCEXT([f90])
ACX_COMPILER_FC_VENDOR_SIMPLE
AC_SUBST([FC_VENDOR], [$acx_cv_fc_compiler_vendor])
ACX_COMPILER_FC_VERSION_SIMPLE
AC_SUBST([FC_VERSION], [$acx_cv_fc_compiler_version])

AC_ARG_VAR([ICON_FCFLAGS],
  [Fortran compiler flags to be appended to FCFLAGS when compiling and linking
ICON, e.g. optimization flags])dnl
ICON_FCFLAGS=${ICON_FCFLAGS-}

AC_ARG_VAR([ICON_LDFLAGS],
  [Fortran compiler flags to be appended to LDFLAGS when linking ICON, e.g.
to prevent configuration-time failures of the Libtool-based bundled
libraries])dnl
ICON_LDFLAGS=${ICON_LDFLAGS-}

AC_ARG_VAR([ICON_OCEAN_FCFLAGS],
  [Fortran compiler flags to be appended to FCFLAGS when compiling the ocean
component of ICON @<:@default=$ICON_FCFLAGS@:>@])dnl
ICON_OCEAN_FCFLAGS=${ICON_OCEAN_FCFLAGS-$ICON_FCFLAGS}

AC_ARG_VAR([ICON_BUNDLED_FCFLAGS],
  [Fortran compiler flags to be appended to FCFLAGS when compiling and linking
the bundled libraries @<:@default=$ICON_FCFLAGS@:>@])dnl
ICON_BUNDLED_FCFLAGS=${ICON_BUNDLED_FCFLAGS-$ICON_FCFLAGS}
subdir_FCFLAGS=$FCFLAGS

dnl Run the delayed cross-compilation check:
AS_VAR_IF([acx_save_cross_compiling], [no],
  [acx_compiler_cross_check_result=ok
   AC_LINK_IFELSE([AC_LANG_PROGRAM],
     [./conftest$ac_cv_exeext >/dev/null 2>&1 || dnl
acx_compiler_cross_check_result=fail],
     [acx_compiler_cross_check_result=fail])
   AS_VAR_IF([acx_compiler_cross_check_result], [fail],
     [AC_MSG_WARN([unable to link and run executables linked with the dnl
provided linker flags, which might result into a failure of the following dnl
standard cross compilation check with a misleading message: make sure that dnl
the dynamic linker is able to find all the shared libraries specified in dnl
LIBS, e.g. by extending the LD_LIBRARY_PATH environment variable (or its dnl
analogue)])])])
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

dnl Generate RPATH flags:
AC_SUBST([ICON_RPATH_LDFLAGS], [''])
AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([icon_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   for icon_L_arg in $icon_L_args; do
     ASX_VAR_APPEND_UNIQ([ICON_RPATH_LDFLAGS],
       ["$acx_cv_fc_rpath_flag$icon_L_arg"], [' '])
   done
   ACX_SHLIB_RPATH_FLAGS_CHECK([$ICON_RPATH_LDFLAGS])])

ACX_FC_ATTRIBUTE_CONTIGUOUS(
  [AS_VAR_APPEND([FCFLAGS],
  [" ${FC_PP_DEF}HAVE_FC_ATTRIBUTE_CONTIGUOUS"])], [:])

AS_VAR_APPEND([FCFLAGS],
  [" ${FC_PP_DEF}__ICON__"])

AS_VAR_IF([enable_atmo], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_ATMO__"])])

AS_VAR_IF([enable_ocean], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_OCEAN__"])])

AS_VAR_IF([enable_jsbach], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_JSBACH__"])])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}YAC_coupling"])])

AS_VAR_IF([enable_concurrent_psrad], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_PS_RAD__"])])

AS_VAR_IF([enable_ecrad], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__ECRAD"])])

AS_VAR_IF([enable_rte_rrtmgp], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_RTE_RRTMGP__"])])

AS_VAR_IF([enable_rttov], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__USE_RTTOV"])])

AS_VAR_IF([enable_dace], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__DACE__"])])

AS_VAR_IF([enable_emvorado], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}HAVE_RADARFWO ${FC_PP_DEF}NUDGING dnl
${FC_PP_DEF}TWOMOM_SB_NEW"])])

AS_VAR_IF([enable_art], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}__ICON_ART ${FC_PP_DEF}__ART_TMP_IFDEF"])])

AS_VAR_IF([enable_mpi], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NOMPI"])])

AS_VAR_IF([enable_active_target_sync], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}HAVE_SLOW_PASSIVE_TARGET_ONESIDED"])])

AS_VAR_IF([enable_parallel_netcdf], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_PARALLEL_NETCDF"])])

AS_VAR_IF([enable_cdi_pio], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_CDI_PIO"])])

AS_VAR_IF([enable_sct], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__SCT__"])])

AS_VAR_IF([enable_yaxt], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_YAXT"])])

AS_IF([test x"$enable_claw" != xno],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}_CLAW"])])

AS_IF([test x"$enable_serialization" != xno],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}SERIALIZE"])
   AS_CASE(["$enable_serialization"],
     [read], [AS_VAR_APPEND([FCFLAGS],
                [" ${FC_PP_DEF}SERIALIZE_READ_REFERENCE"])],
     [perturb], [AS_VAR_APPEND([FCFLAGS],
                   [" ${FC_PP_DEF}SERIALIZE_PERTURB_REFERENCE"])],
     [create], [AS_VAR_APPEND([FCFLAGS],
                  [" ${FC_PP_DEF}SERIALIZE_CREATE_REFERENCE"])])])

AS_VAR_IF([enable_testbed], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_TESTBED__"])])

AS_VAR_IF([enable_loop_exchange], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__LOOP_EXCHANGE"])])

AS_IF([test x"$enable_intel_consistency" = xauto && dnl
test x"$acx_cv_fc_compiler_vendor" = xintel],
  [AC_MSG_NOTICE([Intel compiler directives enforcing consistency are dnl
enabled because $FC is recognized as the Intel Fortran compiler: disable dnl
the directives if required (--disable-intel-consistency)])
   enable_intel_consistency=yes])
AS_VAR_IF([enable_intel_consistency], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}IFORT_CONSISTENCY_ENFORCE"])])

AS_VAR_IF([enable_vectorized_lrtm], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}LRTM_FULL_VECTORIZATION"])])

AS_VAR_IF([enable_mixed_precision], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}__MIXED_PRECISION"])])

######################### C COMPILER COMMAND AND FLAGS #########################
dnl Override the default flags set by AC_PROG_CC:
CFLAGS=${CFLAGS-}
AC_PROG_CC
AC_LANG([C])

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

ACX_COMPILER_CC_VENDOR_SIMPLE
ACX_COMPILER_CC_VERSION_SIMPLE

AC_ARG_VAR([ICON_CFLAGS],
  [C compiler flags to be appended to CFLAGS when compiling ICON,
e.g. optimization flags])dnl
ICON_CFLAGS=${ICON_CFLAGS-}

AC_ARG_VAR([ICON_BUNDLED_CFLAGS],
  [C compiler flags to be appended to CFLAGS when compiling and linking
the bundled libraries @<:@default=$ICON_CFLAGS@:>@])dnl
ICON_BUNDLED_CFLAGS=${ICON_BUNDLED_CFLAGS-$ICON_CFLAGS}
subdir_CFLAGS=$CFLAGS

dnl Check for headers that are protected with HAVE_* macros in the code
dnl (override AC_INCLUDES_DEFAULT of autoconf 2.69, which runs additional
dnl annoying checks; already fixed in the autoconf master branch):
AC_CHECK_HEADERS([execinfo.h ucontext.h link.h], [], [],
  [[// override default includes]])

dnl Check for headers that are not protected at all:
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

########################### CUDA COMPILER AND FLAGS ############################

dnl Don't use neither AS_IF nor AS_VAR_IF to prevent checking for NVCC when
dnl GPU support is disabled:
if test x"$enable_gpu" = xyes; then :
  ACX_PROG_NVCC
   AS_VAR_IF([acx_cv_prog_nvcc_cuda], [no],
     [AC_MSG_FAILURE([unable to find a working NVIDIA CUDA Compiler: either dnl
disable GPU support (--disable-gpu) or set NVCC accordingly])])
   ACX_PROG_NVCC_CXX11
   ACX_LANG_PUSH_CUDA
   AC_CHECK_HEADER([cuda.h], [],
     [AC_MSG_FAILURE([unable to find header file cuda.h: either disable GPU dnl
support (--disable-gpu) or set NVCFLAGS accordingly])],
     [[// override default includes]])
   ACX_LANG_POP_CUDA
   AC_LANG_PUSH([Fortran])
   ACX_FC_CUDA_COMPATIBLE([],
     [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
CUDA compiler: a possible reason is that it cannot find CUDA runtime dnl
libraries (e.g. '-lcudart'): set LDFLAGS and LIBS accordingly])])
   ACX_FC_CUDA_STDCXX_COMPATIBLE([],
     [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
CUDA compiler that require C++ standard library: a possible reason is that dnl
it cannot find C++ standard library (i.e. '-lstdc++'): set LDFLAGS and LIBS dnl
accordingly])])
   AC_LANG_POP([Fortran])
else
  NVCC=
  NVCFLAGS=
fi

################################ LIBRARY CHECKS ################################

dnl One of the topological orderings of the ICON dependency graph
dnl (run ./utils/icon_sorted_deps.sh):
dnl
dnl '-lself', '-licon_tixi', '-lyac', '-lmtime', '-lSerialboxFortranShared',
dnl '-lcdi_f2003 -lcdi', '-leccodes_f90', '-lyaxt -lyaxt_c', '-lsct',
dnl '-lradiance -lrttov10.2 (etc.)', '-llapack -lblas (or other)',
dnl '-lradiation -lifsrrtm -lutilities -lifsaux', '-lrrtmgp -lrte',
dnl '-lnetcdff', '-lhdf5hl_fortran -lhdf5_fortran', '-lz', '-lgomp (or other)',
dnl '-lmpifort -lmpi (or other)'
dnl
dnl We check these dependencies in the reverse order.

dnl Currently, we assume that the Fortran compiler implicitly links to libm and
dnl libc, which are required by ICON itself and some of the bundled libraries
dnl (including the mandatory ones).

dnl Most of the following checks are for Fortran compiler:
AC_LANG([Fortran])

dnl Additional macro definitions to aid the dependency generator,
dnl e.g. compiler intrinsic macros:
AC_SUBST([DEPGEN_FCFLAGS], [''])

dnl A list of library files (i.e. *.a files) of the enabled bundled libraries:
AC_SUBST([BUNDLED_LIBFILES], [''])

dnl A list of include search flags targeting directories of the bundled
dnl libraries:
bundled_includes=

AS_VAR_IF([enable_mpi], [yes],
  [AC_MSG_CHECKING([whether the Fortran 90 interface of MPI library is dnl
available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
       use mpi, only: mpi_init
       implicit none
       integer :: error
       call mpi_init(error)]])],
     [AC_MSG_RESULT([yes])
      AS_IF([test x"$icon_dace_enabled_or_needed" = xyes || dnl
test x"$enable_rttov" = xyes],
        [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_MPI_MOD"])])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interface of MPI library: either disable MPI support (--disable-mpi) or set dnl
FC to the MPI Fortran compiler wrapper, or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])])

dnl The configure scripts of the most of the bundled libraries expect
dnl MPI_LAUNCH to be set in the form of an absolute path, which we try to
dnl accommodate when required and if possible:
AS_IF([test -n "$MPI_LAUNCH"],
  [AS_CASE(["$MPI_LAUNCH"],
     [[[\\/]]* | ?:[[\\/]]*], [],
     [ACX_PROG_SEARCH_ABSPATH([$MPI_LAUNCH],
        [MPI_LAUNCH="$acx_prog_search_abspath"],
        [AC_MSG_WARN([unable to find the absolute path to $MPI_LAUNCH: the dnl
bundled libraries might pick other MPI launcher commands for their tests])])])])

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

AS_VAR_IF([enable_gpu], [yes],
  [ACX_LANG_OPENACC_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openacc_flag"])
      ACX_LANG_MACRO_CHECK_VALUE([_OPENACC],
         [201711 201510 201308 201306 201111])
      AS_VAR_IF([acx_macro_value], [unknown],
        [AC_MSG_WARN([unable to detect supported OpenACC version])
         AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENACC"])],
        [AS_IF([test "$acx_macro_value" -lt 201711 >&AS_MESSAGE_LOG_FD 2>&1],
           [AC_MSG_WARN([Fortran compiler does not support OpenACC 2.6 dnl
(201711): it is strongly recommended to disable GPU support (--disable-gpu)])])
      AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FC_PP_DEF}_OPENACC=$acx_macro_value"])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenACC support by default: dnl
either enable GPU support (--enable-gpu) or set FCFLAGS accordingly to dnl
disable OpenACC])])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_MSG_CHECKING([whether C ZLIB library is available via the Fortran dnl
ISO_C_BINDING interface])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      implicit none
      interface
      subroutine gzOpen() bind(c, name="gzopen")
      end subroutine gzOpen
      end interface
      call gzOpen]])],
     [AC_MSG_RESULT([yes])
      AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}WITH_ZLIB"])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using ZLIB library via dnl
the Fortran ISO_C_BINDING interface: either disable the radar forward dnl
operator EMVORADO (--disable-emvorado) or set LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_MSG_CHECKING([whether the Fortran 90 low- and high-level interfaces of dnl
HDF5 library are available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use hdf5, only: h5open_f, &
                    & h5fopen_f, &
                    & hid_t, &
                    & H5F_ACC_RDONLY_F
      use h5lt, only: h5ltget_attribute_string_f
      implicit none
      integer :: ierr
      integer(kind=hid_t) :: file_id
      character(len=500) :: catt_value
      call h5open_f(ierr)
      call h5fopen_f('conftest.h5', H5F_ACC_RDONLY_F, &
                   & file_id, ierr)
      call h5ltget_attribute_string_f(file_id, &
                                    & 'group', 'attr', &
                                    & catt_value, ierr)]])],
     [AC_MSG_RESULT([yes])
      AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HDF5_RADAR_INPUT"])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
low- and high-level interfaces of HDF5 library: set FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])])

AC_MSG_CHECKING([whether the Fortran 77 interface of NetCDF library is dnl
available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      implicit none
      include 'netcdf.inc'
      integer :: status, ncid
      status = nf_open('conftest.nc', nf_nowrite, ncid)]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using the Fortran 77 dnl
interface of NetCDF library: set FCFLAGS, LDFLAGS and LIBS accordingly])])
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
        [AC_MSG_FAILURE([the provided version of NetCDF does not support the dnl
parallel features: either disable usage of the parallel features of NetCDF dnl
(--disable-parallel-netcdf) or set FCFLAGS, LDFLAGS and LIBS accordingly])],
        [test x"$icon_netcdf_supports_parallel" = xunknown],
        [AC_MSG_WARN([unable to check whether the provided version of dnl
NetCDF supports the parallel features])])])
AC_MSG_CHECKING([whether the Fortran 90 interface of NetCDF library is dnl
available])
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use netcdf, only: nf90_open, nf90_nowrite
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', nf90_nowrite, ncid)]])],
  [AC_MSG_RESULT([yes])
   AS_VAR_IF([enable_emvorado], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NETCDF"])])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interface of NetCDF library: set FCFLAGS, LDFLAGS and LIBS accordingly])])

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
         AC_MSG_FAILURE([unable to link a test program using RTE+RRTMGP dnl
library: either disable usage of the RTE+RRTMGP toolbox for radiation dnl
calculations (--disable-rte-rrtmgp) or use the bundled version of RTE+RRTMGP dnl
library (--without-external-rte-rrtmgp), or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])],
     [icon_subdir_srcdir="$srcdir/externals/rte-rrtmgp"
      icon_subdir_unique_file="rte/mo_rte_sw.F90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the RTE+RRTMGP library dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable usage of dnl
the RTE+RRTMGP toolbox for radiation calculations (--disable-rte-rrtmgp) or dnl
provide an external version of RTE+RRTMGP library dnl
(--with-external-rte-rrtmgp) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      AC_ARG_VAR([ICON_RTE_RRTMGP_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of RTE+RRTMGP library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_RTE_RRTMGP_FCFLAGS=${ICON_RTE_RRTMGP_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_RTE_RRTMGP_CFLAGS],
        [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of RTE+RRTMGP library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_RTE_RRTMGP_CFLAGS=${ICON_RTE_RRTMGP_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_CONFIG_SUBDIR([externals/rte-rrtmgp],
        [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
         [NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([openacc])]],
        [[CFLAGS="$subdir_CFLAGS $ICON_RTE_RRTMGP_CFLAGS"],
         [FCFLAGS="$subdir_FCFLAGS $ICON_RTE_RRTMGP_FCFLAGS"],
         [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
         [--enable-openacc="$enable_gpu"],
         [--enable-silent-rules="$enable_silent_rules"]])
      AS_VAR_APPEND([bundled_includes],
        ["${FC_MOD_IN}externals/rte-rrtmgp/mod "])
      BUNDLED_LIBFILES="externals/rte-rrtmgp/librte.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/rte-rrtmgp/librrtmgp.a $BUNDLED_LIBFILES"])])

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
         AC_MSG_FAILURE([unable to link a test program using ECRAD library: dnl
either disable usage of the ECMWF radiation scheme (--disable-ecrad) or use dnl
the bundled version of ECRAD library (--without-external-ecrad), or set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/ecrad"
      icon_subdir_unique_file="radiation/radiation_flux.F90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the the ECRAD library dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable usage of dnl
the ECMWF radiation scheme (--disable-ecrad) or provide an external version dnl
of ECRAD library (--with-external-ecrad) by setting FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])
      AC_ARG_VAR([ICON_ECRAD_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of ECRAD library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_ECRAD_FCFLAGS=${ICON_ECRAD_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_ECRAD_CFLAGS],
        [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of ECRAD library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_ECRAD_CFLAGS=${ICON_ECRAD_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_CONFIG_SUBDIR([externals/ecrad],
        [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
         [NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([single-precision])]],
        [[CFLAGS="$subdir_CFLAGS $ICON_ECRAD_CFLAGS"],
         [FCFLAGS="$subdir_FCFLAGS $ICON_ECRAD_FCFLAGS"],
         [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
         [--enable-silent-rules="$enable_silent_rules"]])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/ecrad/mod "])
      BUNDLED_LIBFILES="externals/ecrad/libifsaux.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/ecrad/libutilities.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/ecrad/libifsrrtm.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/ecrad/libradiation.a $BUNDLED_LIBFILES"])])

AC_MSG_CHECKING([whether Fortran BLAS library is available])
AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([dgemm])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using BLAS library: set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_MSG_CHECKING([whether Fortran LAPACK library is available])
AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([dgeev])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unable to link a test program using LAPACK library: set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])

AS_VAR_IF([enable_rttov], [yes],
  [AC_MSG_CHECKING([whether Fortran RTTOV library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use rttov_types, only: rttov_coefs
      implicit none
#include "rttov_nullify_coef.interface"
      type(rttov_coefs),target:: coefs(1)
      call rttov_nullify_coef(coefs(1)%coef)]])],
     [AC_MSG_RESULT([yes])
dnl The only file in ICON (src/data_assimilation/interfaces/mo_rttov_ifc.f90)
dnl that uses the macro 'RTTOV12' defines it itself. Therefore, we define it
dnl only when the DACE modules are enabled because their source files need it.
      AS_VAR_IF([enable_dace], [yes],
        [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}RTTOV12"])])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using RTTOV library: dnl
either disable usage of the radiative transfer model for TOVS dnl
(--disable-rttov) or set FCFLAGS, LDFLAGS and LIBS accordingly])])dnl
dnl We need an extra check for the availability of the "god stuff (generalized
dnl optical depth, RTTOV smoothing)" when the DACE modules are enabled.
],
  [AS_VAR_IF([enable_dace], [yes],
     [dnl
dnl It is technically possible to build the DACE modules without RTTOV library
dnl but it is currently unclear whether that makes sense. Therefore, we fail
dnl the configuration when the modules are enabled but usage of the library is
dnl disabled:
      AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with the DACE modules for data assimilation but dnl
usage of the radiative transfer model for TOVS is disabled: either disable dnl
the DACE modules for data assimilation (--disable-dace) or enable usage of dnl
the radiative transfer model for TOVS (--enable-rttov)])
      AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NO_RTTOV"])])])

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
         AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of SCT library: either disable the SCT timer (--disable-sct) or dnl
use the bundled version of SCT library (--without-external-sct), or set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/sct"
      icon_subdir_unique_file="src/sct_reporter.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
SCT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the SCT timer (--disable-sct) or provide an external dnl
version of SCT library (--with-external-sct) by setting FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])
      AC_ARG_VAR([ICON_SCT_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of SCT library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_SCT_FCFLAGS=${ICON_SCT_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_SCT_CFLAGS],
        [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of SCT library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_SCT_CFLAGS=${ICON_SCT_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_CONFIG_SUBDIR([externals/sct],
        [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
         [MPICC=*|MPIFC=*|MPI_LAUNCH=*],
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
        [[CFLAGS="$subdir_CFLAGS $ICON_SCT_CFLAGS"],
         [FCFLAGS="$subdir_FCFLAGS $ICON_SCT_FCFLAGS"],
         [MPICC=], [MPIFC=], [--enable-static], [--disable-shared],
         [--disable-maintainer-mode], [--enable-mpi="$enable_mpi"],
         [--enable-openmp="$enable_openmp"], [--disable-rpath],
         [--disable-papi], [--enable-hdf5],
         [--enable-silent-rules="$enable_silent_rules"],
         [--without-libpapi-prefix], [--without-libhdf5-prefix],
         [MPI_LAUNCH="$MPI_LAUNCH"]])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/sct/src "])
      BUNDLED_LIBFILES="externals/sct/src/.libs/libsct.a $BUNDLED_LIBFILES"])])
dnl The bundled version of SCT, which we link statically, might require
dnl additional linker flags:
dnl   1) a flag enabling linking to librt;
dnl   2) a flag enabling the Fortran compiler to link C code that uses OpenMP.
AC_CONFIG_COMMANDS_PRE(
  [AS_CASE([" $subdirs_extra "],
     [*' externals/sct '*],
dnl First, we check whether we already have '-lrt' in the LIBS:
     [AS_CASE([" $LIBS "],
        [*' -lrt '*], [],
dnl If that is not the case, check whether the configure script of SCT "decided"
dnl to include '-lrt' in the LIBS of SCT:
        [ACX_CONFIG_SUBDIR_VAR([sct_LIBS], [externals/sct], [@LIBS@])
         AS_CASE([" $sct_LIBS "],
           [*' -lrt '*],
dnl If yes, prepend LIBS of ICON with '-lrt':
           [LIBS="-lrt $LIBS"])])
dnl Second, get the flags enabling the Fortran compiler to link C code that uses
dnl OpenMP.
      AS_VAR_IF([enable_openmp], [yes],
        [ACX_CONFIG_SUBDIR_VAR([sct_OPENMP_FC_C_FLAGS], [externals/sct],
           [@OPENMP_FC_C_FLAGS@])
dnl Check, whether we already have the flag in the linking command:
         AS_IF([test -n "$sct_OPENMP_FC_C_FLAGS"],
           [AS_CASE([" $FC $FCFLAGS $ICON_FCFLAGS $LDFLAGS $LIBS "],
dnl (note that the following does not properly account for multiple flags, which
dnl is, however, very unlikely):
              [*" $sct_OPENMP_FC_C_FLAGS "*], [],
dnl If that is not the case, prepend libs with the flag:
              [LIBS="$sct_OPENMP_FC_C_FLAGS $LIBS"])])])])])

AS_IF([test x"$enable_yaxt" = xyes || test x"$enable_cdi_pio" = xyes],
  [AS_VAR_IF([with_external_yaxt], [yes],
     [AC_MSG_CHECKING([whether the Fortran interface of YAXT library is dnl
available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use yaxt, only: xt_initialized, xt_initialize
      implicit none
      integer :: a
      if (.not. xt_initialized()) call xt_initialize(a)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of YAXT library: either disable the YAXT data exchange dnl
(--disable-yaxt) or use the bundled version of YAXT library dnl
(--without-external-yaxt), or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/yaxt"
      icon_subdir_unique_file="src/xt_idxlist.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
YAXT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the YAXT data exchange (--disable-yaxt) or provide an external dnl
version of YAXT library (--with-external-yaxt) by setting FCFLAGS, LDFLAGS dnl
and LIBS accordingly])])
      AC_ARG_VAR([ICON_YAXT_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of YAXT library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_YAXT_FCFLAGS=${ICON_YAXT_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_YAXT_CFLAGS],
        [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of YAXT library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_YAXT_CFLAGS=${ICON_YAXT_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_CONFIG_SUBDIR([externals/yaxt],
        [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|MPI_LAUNCH=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([dependency-tracking])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([regard-for-quality])]],
        [[CFLAGS="$subdir_CFLAGS $ICON_YAXT_CFLAGS"],
         [FCFLAGS="$subdir_FCFLAGS $ICON_YAXT_FCFLAGS"],
         [--enable-static], [--disable-shared], [--disable-maintainer-mode],
         [--enable-dependency-tracking],
         [--enable-silent-rules="$enable_silent_rules"],
         [--with-regard-for-quality="$enable_mpi_checks"],
         [MPI_LAUNCH="$MPI_LAUNCH"]])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/yaxt/src "])
      BUNDLED_LIBFILES="externals/yaxt/src/.libs/libyaxt.a dnl
externals/yaxt/src/.libs/libyaxt_c.a $BUNDLED_LIBFILES"])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_MSG_CHECKING([whether the Fortran 90 interfaces of GRIB-API library is dnl
available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use grib_api, only: grib_open_file, GRIB_SUCCESS
      implicit none
      integer :: igrib, ierr
      call grib_open_file(igrib, 'conftest.grb', 'w', ierr)
      if (ierr  /= GRIB_SUCCESS) igrib = 0]])],
     [AC_MSG_RESULT([yes])
      AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}GRIBAPI"])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interfaces of GRIB-API library: set FCFLAGS, LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([with_external_cdi], [yes],
dnl Check whether an external version of CDI is provided:
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
      AC_MSG_FAILURE([unable to link a test program using the Fortran 2003 dnl
interface of CDI 1.8.x library: either use the bundled version of CDI dnl
library (--without-external-cdi) or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
dnl Check whether the provided external version is compatible with ICON:
   AC_MSG_CHECKING([whether the provided version of CDI 1.8.x is patched])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: t_CdiIterator, cdiIterator_inqFiletype
      implicit none
      integer :: a
      type(t_CdiIterator) :: b
      a = cdiIterator_inqFiletype(b)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([the provided version of CDI 1.8.x does not have the dnl
required patches: use the bundled version of CDI library dnl
(--without-external-cdi)])])
dnl Check whether the provided external version supports the required file
dnl formats:
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
$icon_cdi_format_name: use the bundled version of CDI library dnl
(--without-external-cdi)])],
        [test x"$icon_cdi_format_result" = xunknown],
        [AC_MSG_WARN([unable to check whether the provided version of CDI dnl
supports $icon_cdi_format_name])])
   done
dnl Check whether the provided external version supports the parallel features:
   AS_VAR_IF([enable_cdi_pio], [yes],
     [AS_VAR_IF([with_external_yaxt], [no],
        [AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with the parallel features of an external CDI dnl
library and the bundled YAXT library: either disable the parallel features dnl
of CDI (--disable-cdi-pio) or use the bundled version of CDI library dnl
(--without-external-cdi), or provide the version of YAXT library dnl
(--with-external-yaxt) that the provided CDI library has been built with by dnl
setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      AC_MSG_CHECKING([whether the Fortran interface of CDI-PIO library is dnl
available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use yaxt, only: xt_idxlist
      implicit none
      include 'cdipio.inc'
      integer :: cdiCellGridID, gridtype, n_glb, &
               & max_cell_connectivity, grid_deco_part(2)
      type(xt_idxlist), allocatable :: reorder_idxlst_xt(:)
      type(xt_idxlist) :: null_idxlist
      cdiCellGridID = &
      &  cdiPioDistGridCreate(gridtype, n_glb, -1, -1, &
                            & max_cell_connectivity, &
                            & grid_deco_part, &
                            & reorder_idxlst_xt(1), &
                            & null_idxlist, null_idxlist)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of CDI-PIO library: either disable usage of the parallel features dnl
of CDI (--disable-cdi-pio) or set FCFLAGS, LDFLAGS and LIBS accordingly])])])],
dnl Configure the bundled version of CDI:
  [AS_VAR_IF([enable_cdi_pio], [yes],
     [AC_MSG_ERROR([the bundled version of CDI library does not currently dnl
support the parallel features: either disable usage of the parallel features dnl
of CDI (--disable-cdi-pio) or provide an external version of CDI library dnl
(--with-external-cdi) by setting FCFLAGS, LDFLAGS and LIBS accordingly])])
   icon_subdir_srcdir="$srcdir/externals/cdi"
   icon_subdir_unique_file="src/cdi.h"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of CDI dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of CDI library (--with-external-cdi) by setting FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])
   AC_ARG_VAR([ICON_CDI_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of CDI library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_CDI_FCFLAGS=${ICON_CDI_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_CDI_CFLAGS],
     [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of CDI library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_CDI_CFLAGS=${ICON_CDI_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_CONFIG_SUBDIR([externals/cdi],
     [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|MPI_LAUNCH=*|F77=*|CXX=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
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
     [[CFLAGS="$subdir_CFLAGS $ICON_CDI_CFLAGS"],
      [FCFLAGS="$subdir_FCFLAGS $ICON_CDI_FCFLAGS"],
      [--enable-static], [--disable-shared], [--disable-maintainer-mode],
      [--enable-silent-rules="$enable_silent_rules"],
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
      [acx_cv_have_libnc_dap=no], [acx_cv_have_nc4hdf5=no],
dnl Disable the compilers we do not need.
      [F77=no],
dnl Note that we need to set CXX to 'no' not only to disable the redundant
dnl checks but also to avoid the inconsistency in the generation of the
dnl 'postdeps' variable for Fortran (!!!) compiler inside the *.la files. We
dnl care about those files because they are used when linking the test suite of
dnl CDI (i.e. when running 'make check'). If you want to know more, take a look
dnl at the macro _LT_LANG_FC_CONFIG in libtool.m4 and note that it does not
dnl initialize the variable 'output_verbose_link_cmd' as it should but uses the
dnl value that is set for the C++ compiler. Therefore, by disabling CXX we get
dnl the same behaviour (still questionable but satisfactory for us) as if CDI
dnl didn't use the C++ compiler at all (e.g. like in YAXT).
      [CXX=no],
      [MPI_LAUNCH="$MPI_LAUNCH"]])
   AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/cdi/src "])
   BUNDLED_LIBFILES="externals/cdi/src/.libs/libcdi_f2003.a dnl
externals/cdi/src/.libs/libcdi.a $BUNDLED_LIBFILES"])

AS_IF([test x"$enable_serialization" != xno],
  [AC_MSG_CHECKING(
     [whether the Fortran interface of SERIALBOX2 library is available])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use utils_ppser, only: ppser_set_mode
      implicit none
      call ppser_set_mode(1)]])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_FAILURE([unable to link a test program using the Fortrain dnl
interface of SERIALBOX2 library: either disable the serialization dnl
(--disable-serialization) or set FCFLAGS, LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([with_external_mtime], [yes],
  [AC_MSG_CHECKING([whether the Fortran interface of MTIME library is dnl
available])
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
      AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of MTIME library: either use the bundled version of MTIME dnl
(--without-external-mtime) or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
  [icon_subdir_srcdir="$srcdir/externals/mtime"
   icon_subdir_unique_file="include/mtime_calendar.h"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of MTIME dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of MTIME library (--with-external-mtime) by setting CFLAGS, dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])
   AC_ARG_VAR([ICON_MTIME_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of MTIME library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_MTIME_FCFLAGS=${ICON_MTIME_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_MTIME_CFLAGS],
     [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of MTIME library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_MTIME_CFLAGS=${ICON_MTIME_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_CONFIG_SUBDIR([externals/mtime],
     [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([examples])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([openmp])]],
     [[CFLAGS="$subdir_CFLAGS $ICON_MTIME_CFLAGS"],
      [FCFLAGS="$subdir_FCFLAGS $ICON_MTIME_FCFLAGS"],
      [--enable-static], [--disable-shared], [--disable-maintainer-mode],
      [--enable-silent-rules="$enable_silent_rules"], [--disable-examples],
      [--enable-openmp="$enable_openmp"]])
   AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/mtime/src "])
   BUNDLED_LIBFILES="externals/mtime/src/.libs/libmtime.a $BUNDLED_LIBFILES"])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_IF([with_external_yac], [yes],
     [AS_VAR_IF([with_external_mtime], [no],
        [AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with an external YAC library and the bundled dnl
MTIME library: either disable the coupling (--disable-coupling) or use the dnl
bundled version of YAC library (--without-external-yac), or provide the dnl
version of MTIME library (--with-external-mtime) that the provided YAC dnl
library has been built with by setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      AC_MSG_CHECKING([whether the Fortran interface of YAC library is dnl
available])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_yac_finterface, only: yac_finit
      implicit none
      character(len=1024) :: a, b
      call yac_finit(a, b)]])],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of YAC library: either disable the coupling (--disable-coupling) dnl
or use the bundled version of YAC library (--without-external-yac), or set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])],
     [icon_subdir_srcdir="$srcdir/externals/yac"
      icon_subdir_unique_file="src/grid.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of YAC dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the coupling (--disable-coupling) or provide an external version of dnl
YAC library (--with-external-yac) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      yac_MTIME_CFLAGS= ; yac_MTIME_CLIBS=
      AS_VAR_IF([with_external_mtime], [no],
        [dnl Set $ac_abs_top_srcdir and $ac_abs_top_builddir:
         _AC_SRCDIRS([.])
         yac_MTIME_CFLAGS="-I${ac_abs_top_srcdir}/externals/mtime/include"
         yac_MTIME_CLIBS="${ac_abs_top_builddir}/externals/mtime/src/.libs/libmtime.a"])
      AC_ARG_VAR([ICON_YAC_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of YAC library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_YAC_FCFLAGS=${ICON_YAC_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_YAC_CFLAGS],
        [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of YAC library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_YAC_CFLAGS=${ICON_YAC_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_CONFIG_SUBDIR([externals/yac],
        [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|MPI_LAUNCH=*],
         [MTIME_CFLAGS=*|MTIME_CLIBS=*|NETCDF_CFLAGS=*|NETCDF_CLIBS=*],
         [XML2_CFLAGS=*|XML2_CLIBS=*|MPI_CFLAGS=*|MPI_CLIBS=*],
         [MPI_FCFLAGS=*|MPI_FCLIBS=*|FORTRAN_LAPACK_CLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([lib-only])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([xml])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([netcdf])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpi])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([mpi-checks])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([external-lapack])],
         [ACX_CONFIG_SUBDIR_PATTERN_WITH([external-mtime])]],
        [[CFLAGS="$subdir_CFLAGS $ICON_YAC_CFLAGS"],
         [FCFLAGS="$subdir_FCFLAGS $ICON_YAC_FCFLAGS"],
         [MTIME_CFLAGS="$yac_MTIME_CFLAGS"], [MTIME_CLIBS=],
         [NETCDF_CFLAGS=], [NETCDF_CLIBS=], [XML2_CFLAGS=], [XML2_CLIBS=],
         [MPI_CFLAGS=], [MPI_CLIBS=], [MPI_FCFLAGS=], [MPI_FCLIBS=],
         [FORTRAN_LAPACK_CLIBS=], [MPI_LAUNCH="$MPI_LAUNCH"],
         [--disable-maintainer-mode],
         [--enable-silent-rules="$enable_silent_rules"], [--enable-lib-only],
         [--enable-xml], [--enable-netcdf], [--enable-mpi="$enable_mpi"],
         [--with-external-lapack=fortran], [--with-external-mtime],
         [--enable-mpi-checks="$enable_mpi_checks"],
         [acx_cv_c_lib_func_datetimeToString="$yac_MTIME_CLIBS"]])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/yac/src "])
      BUNDLED_LIBFILES="externals/yac/src/libyac.a $BUNDLED_LIBFILES"])])

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
         AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of TIXI library: either disable the aerosols and reactive trace dnl
component (--disable-art) or use the bundled version of TIXI library dnl
(--without-external-tixi), or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
     [icon_subdir_srcdir="$srcdir/externals/tixi"
      icon_subdir_unique_file="src/tixi.f90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
TIXI library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the aerosols and reactive trace component (--disable-art) or provide dnl
an external version of TIXI library (--with-external-tixi) by setting dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])
      AC_ARG_VAR([ICON_TIXI_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of TIXI library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_TIXI_FCFLAGS=${ICON_TIXI_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_TIXI_CFLAGS],
        [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of TIXI library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_TIXI_CFLAGS=${ICON_TIXI_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_CONFIG_SUBDIR([externals/tixi],
        [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
         [XML2_CFLAGS=*|XML2_CLIBS=*],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
         [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])]],
        [[CFLAGS="$subdir_CFLAGS $ICON_TIXI_CFLAGS"],
         [FCFLAGS="$subdir_FCFLAGS $ICON_TIXI_FCFLAGS"],
         [XML2_CFLAGS=], [XML2_CLIBS=], [--disable-maintainer-mode],
         [--enable-silent-rules="$enable_silent_rules"]])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/tixi/src "])
      BUNDLED_LIBFILES="externals/tixi/src/libicon_tixi.a dnl
$BUNDLED_LIBFILES"])])

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
library: either use the bundled version of SELF library dnl
(--without-external-self) or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
  [icon_subdir_srcdir="$srcdir/externals/self"
   icon_subdir_unique_file="src/self_object.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of SELF dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of SELF library (--with-external-self) by setting FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])
   AC_ARG_VAR([ICON_SELF_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when compiling and
linking the bundled version of SELF library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_SELF_FCFLAGS=${ICON_SELF_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_SELF_CFLAGS],
     [Fortran compiler flags to be appended to CFLAGS when compiling and
linking the bundled version of SELF library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_SELF_CFLAGS=${ICON_SELF_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_CONFIG_SUBDIR([externals/self],
     [[BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])]],
     [[CFLAGS="$subdir_CFLAGS $ICON_SELF_CFLAGS"],
      [FCFLAGS="$subdir_FCFLAGS $ICON_SELF_FCFLAGS"],
      [--disable-maintainer-mode],
      [--enable-silent-rules="$enable_silent_rules"]])
   AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/self/src "])
   BUNDLED_LIBFILES="externals/self/src/libself.a $BUNDLED_LIBFILES"])

dnl Prepend the include search flags to FCFLAGS:
AS_IF([test -n "$bundled_includes"], [FCFLAGS="${bundled_includes}${FCFLAGS}"])

AS_VAR_IF([enable_gpu], [yes],
  [AS_VAR_IF([with_external_cub], [yes],
     [ACX_LANG_PUSH_CUDA
      AC_CHECK_HEADER([cub/device/device_select.cuh], [],
        [AC_MSG_FAILURE([unable to find header files of CUB library: either dnl
disable GPU support (--disable-gpu) or use the bundled version of CUB dnl
library (--without-external-cub) or set NVCFLAGS accordingly])],
        [[// override default includes]])
      ACX_LANG_POP_CUDA],
     [icon_subdir_srcdir="$srcdir/externals/cub"
      icon_subdir_unique_file="cub/device/device_select.cuh"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of CUB dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable dnl
GPU support (--disable-gpu) or provide an external version of CUB library dnl
(--with-external-cub) by setting NVCFLAGS accordingly])])
      AS_VAR_APPEND([NVCFLAGS], [" -I${srcdir}/externals/cub"])])])

##################### PREPROCESSORS AND AUXILIARY COMMANDS #####################

AC_PROG_INSTALL
AC_PROG_MKDIR_P
dnl The message generated by AC_PROG_LN_S is misleading in our case, therefore
dnl we do what the standard macro does but without the message:
AC_SUBST([LN_S], [$as_ln_s])

AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python python3],
  [icon_prog_try=dnl
'$acx_candidate $srcdir/utils/mkhelper/depgen.py -h >/dev/null'
   AC_TRY_EVAL([icon_prog_try])],
  [PYTHON=$acx_candidate])

AC_ARG_VAR([PERL], [Perl interpreter command])dnl
ACX_PROG_SEARCH([PERL], [perl],
  [icon_prog_try='$acx_candidate $srcdir/utils/pvcs.pl --help >/dev/null'
   AC_TRY_EVAL([icon_prog_try])
   test $ac_status -eq 0 || test $ac_status -eq 1],
  [PERL=$acx_candidate])

AC_ARG_VAR([CLAW], [CLAW preprocessor command (optional)])dnl
AC_ARG_VAR([CLAWFLAGS], [CLAW preprocessor flags (optional)])dnl
dnl Additional macro definitions to aid the dependency generator,
dnl e.g. CLAW intrinsic macros:
AC_SUBST([DEPGEN_CLAWFLAGS], ['-D_CLAW'])
AS_IF([test x"$enable_claw" != xno],
  [icon_subdir_srcdir="$srcdir/externals/omni-xmod-pool"
   icon_subdir_unique_file="cdi"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find precompiled module files for OMNI dnl
compiler required for CLAW preprocessing ($icon_subdir_unique_file) in dnl
'$icon_subdir_srcdir': disable CLAW preprocessing (--disable-claw)])])
   AS_VAR_IF([enable_openmp], [yes],
     [AC_MSG_WARN([CLAW preprocessing is not compatible with OpenMP: it is dnl
strongly recommended to disable OpenMP support (--disable-openmp)])])
dnl
   AS_VAR_IF([enable_gpu], [yes],
     [AS_VAR_APPEND([CLAWFLAGS], [' --target=gpu --directive=acc'])
      AS_VAR_APPEND([DEPGEN_CLAWFLAGS], [' -D_OPENACC'])],
     [AS_VAR_APPEND([CLAWFLAGS], [' --target=cpu --directive=none'])])
   AS_VAR_IF([enable_jsbach], [yes],
     [AS_VAR_APPEND([CLAWFLAGS],
        [" --model-config=$srcdir/externals/jsbach/config/icon_jsbach.toml"])])
   AS_VAR_IF([enable_claw], [validate],
     [AS_VAR_IF([enable_serialization], [no],
        [AC_MSG_ERROR([unsupported combination of the configure options: you dnl
are trying to configure with CLAW preprocessing in VALIDATE mode but the dnl
Serialbox2 serialization is disabled: either enable CLAW preprocessing in a dnl
different mode (--enable-claw=std) or enable the Serialbox2 serialization dnl
(--enable-serialization)])])
dnl Enable serialization support and update after kernels:
      AS_VAR_APPEND([CLAWFLAGS], [' -x=sca_serialization_enabled:true'])
      AS_VAR_APPEND([CLAWFLAGS],
        [' -x=sca_serialization_enabled_direction:write'])
      AS_VAR_APPEND([CLAWFLAGS],
        [' -x=sca_forward_update_enabled_direction:out'])],
     [dnl
dnl Disable serialization support and updates by kernel:
      AS_VAR_APPEND([CLAWFLAGS], [' -x=sca_serialization_enabled:false'])
      AS_VAR_APPEND([CLAWFLAGS], [' -x=sca_forward_update_enabled:false'])])
dnl Check for CLAW:
   AC_LANG_ASSERT([Fortran])
   AC_MSG_CHECKING([for CLAW preprocessor command])
   icon_claw_result=unknown
   CLAW=${CLAW-'clawfc'}
   AS_MKDIR_P([conftest.dir])
   cd conftest.dir
   AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      module conftest_module
      end module]])])
   icon_prog_try='$CLAW -f -o conftest.claw.$ac_ext conftest.$ac_ext'
   AC_TRY_EVAL([icon_prog_try])
   AS_IF([test $ac_status -eq 0 && test -f conftest.claw.$ac_ext && dnl
test -f conftest_module.xmod],
     [mv conftest.claw.$ac_ext conftest.$ac_ext
      AC_COMPILE_IFELSE([], [icon_claw_result=$CLAW])])
   cd ..
   rm -rf conftest.dir
   AC_MSG_RESULT([$icon_claw_result])
   AS_VAR_IF([icon_claw_result], [unknown],
     [AC_MSG_FAILURE([unable to find a valid CLAW preprocessor command: dnl
either disable CLAW preprocessing (--disable-claw) or set CLAW accordingly])])
dnl
   AC_CONFIG_FILES([claw.mk])],
  [CLAW=])

AC_ARG_VAR([FPP], [Fortran preprocessor command (optional)])dnl
AS_VAR_IF([enable_explicit_fpp], [yes],
  [ACX_PROG_FPP([FPP],
     [FPP=$acx_cv_prog_fpp],
     [AC_MSG_FAILURE([unable to find a valid Fortran preprocessor command: dnl
either disable explicit Fortran preprocessing (--disable-explicit-fpp) or dnl
set FPP accordingly])],
     [$FC_INC_FLAG_PP_f90], [$FC_PP_DEF], [utils/fpp-wrappers])],
  [FPP=])

AC_ARG_VAR([SB2PP], [Serialbox2 preprocessor command (optional)])dnl
AS_IF([test x"$enable_serialization" != xno],
  [AC_LANG_ASSERT([Fortran])
   AC_MSG_CHECKING([for Serialbox2 preprocessor command])
   icon_sb2pp_result=unknown
   SB2PP=${SB2PP-'${PYTHON} -m pp_ser'}
dnl pp_ser.py can not process a piece of code like this:
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
     [AC_MSG_FAILURE([unable to find a valid Serialbox2 preprocessor dnl
command: either disable the Serialbox2 serialization dnl
(--disable-serialization) or make sure that path to the 'pp_ser' module is dnl
in the PYTHONPATH environment variable, or set SB2PP accordingly])])],
   [SB2PP=])

dnl We have checked this command earlier but put the declaration here so
dnl it would be listed in the help message of the configure script at the very
dnl end, together with other optional commands.
AC_ARG_VAR([MPI_LAUNCH], [interactive MPI launcher command (optional)])dnl

#################################### OUTPUT ####################################

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  depgen.c.config
  depgen.f90.config
  deplist.config
  icon.mk
])

AC_CONFIG_FILES([run/collect.set-up.info], [chmod +x "$ac_file"])

dnl The runscript generator should call run/collect.set-up.info and generate
dnl run/set-up.info. This is currently not the case, therefore we delegate this
dnl task to the configure script (do not forget to update '*clean' rules in
dnl 'icon.mk.in' when the following becomes redundant):
AC_CONFIG_COMMANDS([run/set-up.info],
  [run/collect.set-up.info "$ac_file" 2>/dev/null || dnl
   AC_MSG_WARN([failed to create "$ac_file"])])

dnl An additional step to generate dummy dependencies:
AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} -f icon.mk dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch icon.mk])

AC_OUTPUT
