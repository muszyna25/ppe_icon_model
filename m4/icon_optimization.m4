# ICON_OPTIMIZATION_ARGS()
# -----------------------------------------------------------------------------
# Sets optimization related command-line arguments of the configure script:
#
#   '--enable-optimization' for choosing of one of the optimization levels:
#       'release'     default optimization for production runs;
#       'aggressive'  the most aggressive optimizations (that are known to
#                     work);
#       'precise'     optimization level ensuring binary reproducibility of the
#                     results;
#       'debug'       optimizations that are unlikely to interfere with
#                     debugging or profiling;
#       'test'        optimization level for the test suite;
#       'none'        all optimizations disabled.
#
#   '--enable-loop-exchange' for enabling loop exchange
#
#   '--enable-intel-consistency' for enabling Intel compiler directive
#                                enforcing consistency
#
#   '--enable-vectorized-lrtm' for enabling full vectorization of LRTM
#
#   '--enable-mixed-precision' for enabling mixed precision dycore
#
AC_DEFUN([ICON_OPTIMIZATION_ARGS],
  [AC_ARG_ENABLE([optimization],
[  --enable-optimization   prepend CFLAGS and FCFLAGS with additional
                          predefined compiler-specific set of flags. The value
                          of the argument must be one of the following
                          @<:@default=release@:>@:
                            release|yes  optimizations for production runs
                            aggressive   most aggressive optimizations
                            precise      optimizations ensuring reproducibility
                            debug        optimizations for debugging
                            test         optimizations for regression testing
                            none|no      no optimizations at all],
[AS_CASE(["$enableval"],
   [yes], [enable_optimization=release],
   [no], [enable_optimization=none],
   [release|aggressive|precise|debug|test|none], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-optimization='$enableval'; valid values are 'release', dnl
'aggressive', 'precise', 'debug', 'test', 'none', dnl
'yes' (same as 'release'), 'no' (same as 'none')])])],
[enable_optimization=release])
dnl
   AC_ARG_ENABLE([loop-exchange],
     [AC_HELP_STRING([--enable-loop-exchange],
        [enable loop exchange @<:@default=auto@:>@])], [],
     [enable_loop_exchange=auto])
dnl
   AC_ARG_ENABLE([intel-consistency],
     [AC_HELP_STRING([--enable-intel-consistency],
        [enable Intel compiler directives enforcing consistency
@<:@default=auto@:>@])], [],
     [enable_intel_consistency=auto])
dnl
   AC_ARG_ENABLE([vectorized-lrtm],
     [AC_HELP_STRING([--enable-vectorized-lrtm],
        [enable the parallelization-invariant version of LRTM
@<:@default=auto@:>@])], [],
     [enable_vectorized_lrtm=auto])
dnl
   AC_ARG_ENABLE([mixed-precision],
     [AC_HELP_STRING([--enable-mixed-precision],
        [enable mixed precision dycore @<:@default=auto@:>@])], [],
     [enable_mixed_precision=auto])])

# ICON_OPTIMIZATION_SET_FCFLAGS()
# -----------------------------------------------------------------------------
# Sets shell variables (see the list below) to be prepended to FCFLAGS when
# compiling different components of ICON. The variables are set according to
# the optimization level (see ICON_OPTIMIZATION_ARG_ENABLE) stored in the
# enable_optimization shell variable and to the version of the Fortran compiler
# in use.
#
# Currently, sets the following variables:
#
#     icon_optim_FCFLAGS
#         flags to be prepended to FCFLAGS when compiling ICON and its
#         components (including JSBACH and ART, excluding OCEAN);
#
#     icon_optim_ocean_FCFLAGS
#         flags to be prepended to FCFLAGS when compiling files ICON-OCEAN
#         (src/hamocc/%.f90, src/ocean/%.f90, src/sea_ice/%.f90)
#
#     icon_optim_subdir_FCFLAGS
#         flags to be prepended to FCFLAGS when running the configure scripts
#         of the bundled libraries.
#
AC_DEFUN([ICON_OPTIMIZATION_SET_FCFLAGS],
  [AC_REQUIRE([ACX_COMPILER_FC_VERSION])dnl
   AC_PROVIDE_IFELSE([ICON_OPTIMIZATION_ARGS], [],
     [m4_warn([syntax],
        [$0 should be called after ICON_OPTIMIZATION_ARGS])])dnl
dnl The following code is M4-quoted and is implemented using plain shell to
dnl be more maintainable by those who are less familiar with M4 syntax:
[
   icon_optim_error=
   icon_optim_FCFLAGS=
   icon_optim_ocean_FCFLAGS=
   icon_optim_subdir_FCFLAGS=
   case $enable_optimization in #(
     release)
       # release, default (example)
       icon_optim_FCFLAGS='-O2' ;; #(
     aggressive)
       case $acx_cv_fc_compiler_vendor in #(
         unknown)
           # aggressive, unknown compiler
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for unknown Fortran compiler" ;; #(
         *)
           # aggressive, default (example)
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for $acx_cv_fc_compiler_vendor Fortran compiler: the relevant \
set of extra FCFLAGS is not defined" ;;
       esac ;; #(
     precise)
       case $acx_cv_fc_compiler_vendor in #(
         unknown)
           # precise, unknown compiler
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for unknown Fortran compiler" ;; #(
         intel)
           # precise, Intel compiler (example)
           icon_optim_FCFLAGS='-O2 -fp-model precise'
           test x"$enable_intel_consistency" = xauto && \
             enable_intel_consistency=yes
           test x"$enable_vectorized_lrtm" = xauto && \
             enable_vectorized_lrtm=yes ;; #(
         *)
           # precise, default (example)
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for $acx_cv_fc_compiler_vendor Fortran compiler: the relevant \
set of extra FCFLAGS is not defined" ;;
       esac ;; #(
     debug)
       # debug, default (example)
       test x"$ac_cv_prog_fc_g" = xyes && icon_optim_FCFLAGS='-g' ;; #(
     test)
       case $acx_cv_fc_compiler_vendor in #(
         intel)
           # test, intel compiler (example)
           case $acx_cv_fc_compiler_version in #(
             1[0-6].*)
               # test, intel compiler, 10>=version<=16 (example)
               icon_optim_FCFLAGS='-O2' ;; #(
             *)
               # test, intel compiler, default (example)
               icon_optim_FCFLAGS='-O1' ;;
           esac ;; #(
         cray)
           # test, cray compiler, default (example)
           icon_optim_FCFLAGS='-O2' ;; #(
         *)
           # test, default (example)
           icon_optim_FCFLAGS='-O1' ;;
       esac ;; #(
     none)
       # none, keep empty
       ;; #(
     *)
       icon_optim_error="unexpected optimization level for Fortran compiler: \
'$enable_optimization'; valid values are 'release', 'aggressive', 'precise', \
'debug', 'test', or 'none'." ;;
   esac

   # Currently, we compile ICON-Ocean with the same flags:
   icon_optim_ocean_FCFLAGS=$icon_optim_FCFLAGS

   # Currently, we compile the bundled libraries with the same flags:
   icon_optim_subdir_FCFLAGS=$icon_optim_FCFLAGS

   # Currently, enable loop exchange by default,
   # regardless of the optimization level:
   test x"$enable_loop_exchange" = xauto && enable_loop_exchange=yes

   # Currently, disable mixed precision by default,
   # regardless of the optimization level:
   test x"$enable_mixed_precision" = xauto && enable_mixed_precision=no]
dnl The plain shell section ends here
dnl
   AS_IF([test -n "$icon_optim_error"],
     [AC_MSG_ERROR([$icon_optim_error])])
dnl
   AS_VAR_IF([enable_loop_exchange], [yes],
     [AS_VAR_APPEND([icon_optim_FCFLAGS],
        [" ${FC_PP_DEF}__LOOP_EXCHANGE"])])
dnl
   AS_VAR_IF([enable_intel_consistency], [yes],
     [AS_VAR_APPEND([icon_optim_FCFLAGS],
        [" ${FC_PP_DEF}IFORT_CONSISTENCY_ENFORCE"])])
dnl
   AS_VAR_IF([enable_vectorized_lrtm], [yes],
     [AS_VAR_APPEND([icon_optim_FCFLAGS],
        [" ${FC_PP_DEF}LRTM_FULL_VECTORIZATION"])])
dnl
   AS_VAR_IF([enable_mixed_precision], [yes],
     [AS_VAR_APPEND([icon_optim_FCFLAGS],
        [" ${FC_PP_DEF}__MIXED_PRECISION ${FC_PP_DEF}__MIXED_PRECISION_2"])])])

# ICON_OPTIMIZATION_SET_CFLAGS()
# -----------------------------------------------------------------------------
# Sets shell variables (see the list below) to be prepended to CFLAGS when
# compiling different components of ICON. The variables are set according to
# the optimization level (see ICON_OPTIMIZATION_ARG_ENABLE) stored in the
# enable_optimization shell variable and to the version of the C compiler in
# use.
#
# Currently, sets the following variables:
#
#     icon_optim_CFLAGS
#         flags to be prepended to CFLAGS when compiling ICON and its
#         components (including JSBACH and ART);
#
#     icon_optim_subdir_CFLAGS
#         flags to be prepended to CFLAGS when running the configure scripts of
#         the bundled libraries.
#
AC_DEFUN([ICON_OPTIMIZATION_SET_CFLAGS],
  [AC_REQUIRE([ACX_COMPILER_CC_VERSION])dnl
   AC_PROVIDE_IFELSE([ICON_OPTIMIZATION_ARGS], [],
     [m4_warn([syntax],
        [$0 should be called after ICON_OPTIMIZATION_ARGS])])dnl
dnl The following code is M4-quoted and is implemented using plain shell to
dnl be more maintainable by those who are less familiar with M4 syntax:
[
   icon_optim_error=
   icon_optim_CFLAGS=
   icon_optim_subdir_CFLAGS=
   case $enable_optimization in #(
     release)
       # release, default (example)
       icon_optim_CFLAGS='-O2' ;; #(
     aggressive)
       case $acx_cv_c_compiler_vendor in #(
         unknown)
           # aggressive, unknown compiler
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for unknown C compiler" ;; #(
         *)
           # aggressive, default (example)
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for $acx_cv_c_compiler_vendor C compiler: the relevant \
set of extra CFLAGS is not defined" ;;
       esac ;; #(
     precise)
       case $acx_cv_c_compiler_vendor in #(
         unknown)
           # precise, unknown compiler
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for unknown C compiler" ;; #(
         intel)
           # precise, Intel compiler (example)
           ;; #(
         *)
           # precise, default (example)
           icon_optim_error="unable to set '$enable_optimization' \
optimizations for $acx_cv_c_compiler_vendor C compiler: the relevant \
set of extra CFLAGS is not defined" ;;
       esac ;; #(
     debug)
       # debug, default (example)
       test x"$ac_cv_prog_cc_g" = xyes && icon_optim_CFLAGS='-g' ;; #(
     test)
       case $acx_cv_c_compiler_vendor in #(
         intel)
           # test, intel compiler (example)
           case $acx_cv_c_compiler_version in #(
             1[0-6].*)
               # test, intel compiler, 10>=version<=16 (example)
               icon_optim_CFLAGS='-O2' ;; #(
             *)
               # test, intel compiler, default (example)
               icon_optim_CFLAGS='-O1' ;;
           esac ;; #(
         cray)
           # test, cray compiler, default (example)
           icon_optim_CFLAGS='-O2' ;; #(
         *)
           # test, default (example)
           icon_optim_CFLAGS='-O1' ;;
       esac ;; #(
     none)
       # none, keep empty
       ;; #(
     *)
       icon_optim_error="unexpected optimization level for C compiler: \
'$enable_optimization'; valid values are 'release', 'aggressive', 'precise', \
'debug', 'test', or 'none'." ;;
   esac

   # Currently, we compile the bundled libraries with the same flags:
   icon_optim_subdir_CFLAGS=$icon_optim_CFLAGS]
dnl The plain shell section ends here
dnl
   AS_IF([test -n "$icon_optim_error"],
     [AC_MSG_ERROR([$icon_optim_error])])])

# ICON_OPTIMIZATION_CHECK_FCFLAGS()
# -----------------------------------------------------------------------------
# Checks whether the optimization flags set by ICON_OPTIMIZATION_SET_FCFLAGS
# and stored in the icon_optim_FCFLAGS shell variable are accepted by the
# Fortran compiler.
#
# If successful, does nothing, otherwise fails with an error message.
#
AC_DEFUN([ICON_OPTIMIZATION_CHECK_FCFLAGS],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_REQUIRE([ICON_OPTIMIZATION_SET_FCFLAGS])dnl
   _ICON_OPTIMIZATION_CHECK])

# ICON_OPTIMIZATION_CHECK_CFLAGS()
# -----------------------------------------------------------------------------
# Checks whether the optimization flags set by ICON_OPTIMIZATION_SET_CFLAGS
# and stored in the icon_optim_CFLAGS shell variable are accepted by the C
# compiler.
#
# If successful, does nothing, otherwise fails with an error message.
#
AC_DEFUN([ICON_OPTIMIZATION_CHECK_CFLAGS],
  [AC_LANG_ASSERT([C])dnl
   AC_REQUIRE([ICON_OPTIMIZATION_SET_CFLAGS])dnl
   _ICON_OPTIMIZATION_CHECK])

# _ICON_OPTIMIZATION_CHECK()
# -----------------------------------------------------------------------------
# Checks whether the optimization flags set by
# ICON_OPTIMIZATION_SET_[]_AC_LANG_PREFIX[]FLAGS and stored in the
# icon_optim_[]_AC_LANG_PREFIX[]FLAGS shell variable are accepted by the
# current compiler.
#
# If successful, does nothing, otherwise fails with an error message.
#
m4_define([_ICON_OPTIMIZATION_CHECK],
  [AS_IF([test x"$enable_optimization" != xnone],
     [AC_MSG_CHECKING([whether _AC_LANG compiler accepts flags for dnl
optimization level '$enable_optimization'])
      icon_optimization_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      _AC_LANG_PREFIX[]FLAGS=$icon_optim_[]_AC_LANG_PREFIX[]FLAGS
      icon_optimization_check_result=no
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
        [AC_LINK_IFELSE([], [icon_optimization_check_result=yes])])
      AC_MSG_RESULT([$icon_optimization_check_result])
      AS_VAR_IF([icon_optimization_check_result], [no],
        [AC_MSG_FAILURE([_AC_LANG compiler does not accept flags for dnl
optimization level '$enable_optimization': disable optimizations dnl
(--disable-optimization)])])
      _AC_LANG_PREFIX[]FLAGS=dnl
$icon_optimization_save_[]_AC_LANG_PREFIX[]FLAGS])])
