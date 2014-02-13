dnl Process this file with autoconf to produce a configure script.
dnl
dnl Luis Kornblueh
dnl Uwe Schulzweida
dnl
dnl  - add cross compiling for NEC (use --host=sx6)
dnl  - use autoconf version 2.59 
dnl  - added cross compiling for CRAY XT3
dnl  - added support for IBM AIX

dnl

AC_INIT(icon,1.2.1,luis.kornblueh@zmaw.de,icon-1.2.1)
AC_PREREQ(2.59)

dnl ######################################################################
dnl basic package definitions (by default ICON)
package="ICON"
if test -r "./config/set-up" ; then
. ./config/set-up
fi

dnl Checks for package options.
dnl ######################################################################

AC_PREFIX_DEFAULT(.)

AC_CONFIG_AUX_DIR(config)

AC_CANONICAL_SYSTEM

dnl ######################################################################

AC_ARG_WITH([fortran],
  [AS_HELP_STRING([--with-fortran=COMPILER],
    [select_available compiler (see README.compiler)])],
  [],
  [])


if test x$with_fortran != x; then
    fortran_compiler=$with_fortran
else
    fortran_compiler="default"
fi
make_command="make"

AC_ARG_WITH([openmp],
  [AS_HELP_STRING([--with-openmp],
    [enable OpenMP, if available for selected compiler])],
  [],
  [])

AC_ARG_WITH([setup],
  [AS_HELP_STRING([--with-setup=<name>],
    [use user specific setup (mh_setup)])],
    [],
    [])

if test x$with_setup != x ; then
  mh_setup=$with_setup
else
  mh_setup=""            
fi

AC_ARG_WITH([mtime],
  [AS_HELP_STRING([--without-mtime],
    [disable use of mtime])],
    [],
    [])

AS_IF([test "x$with_mtime" != "xno"],
      [use_mtime="yes"],
      [use_mtime="no"])

AC_ARG_WITH([messy],
  [AS_HELP_STRING([--with-messy],
    [enable Modular Earth Submodel System (MESSy) interface])],
  [],
  [])

AS_IF([test "x$with_messy" = "xyes"],
      [use_messy="yes"],
      [use_messy="no"])

AC_ARG_WITH([flags],
  [AS_HELP_STRING([--with-flags=<name>],
    [use specified group od flags (std, hiopt, noopt, etc])],
    [],
    [])

case $host in
  i*86-*-linux-*)  
    if test x$fortran_compiler != xintel ; then
      gcc_lib="$(gcc -print-libgcc-file-name)"
    fi
    ;;
  *)               
    gcc_lib=""                               
    ;;
esac



dnl ######################################################################

ac_hostname=$(hostname)
ac_hostname=${ac_hostname%%\.*}

AS_CASE([x"$ac_hostname"],
  [xplogintest1],
  [ac_sitename=dkrz.de],
  [xwizard*],
  [ac_sitename=dkrz.de; ac_hostname=wizard],
  [if [[ $(hostname) = $ac_hostname ]] ; then
   fqdn=$(host -n $ac_hostname | awk '/has address/ {print $1}')
else
   fqdn=$(hostname)
fi
  echo $fqdn
  if [[ -n "$fqdn" ]]; then
    ac_sitename=$(echo $fqdn | awk -F. '{i=NF-1; printf("%s.%s",$i,$NF);}')
  else
    if [[ -n "$PE_ENV" ]] ; then
      fqdn=$(echo ${PE_ENV} | awk '{print tolower($0)}')
      echo "detecting Cray environment... $fqdn"
    else
      fqdn=${HOSTNAME}
    fi
    ac_sitename=${fqdn}
  fi])

dnl ######################################################################

host_frag=""
case $host in
    sx*-*-*)             host_frag=config/mh-nec ;;
    *-ibm-aix*)          host_frag=config/mh-aix ;;
    *-*-linux-*)         host_frag=config/mh-linux ;;
    *-*-darwin*)         host_frag=config/mh-darwin ;;
    *-cray-*)            host_frag=config/mh-cray ;;
    *) echo This configuration is not supported.
       echo Please set a valid config/mh-* file up.
       ;;	 	
esac

if test ! -f $host_frag ; then
   echo "machine dependent configuration file $host_frag does not exist!"
   exit 1
fi

if test "$build" != "$host" ; then
   CROSS=yes
else
   CROSS=""
fi

if test -z "$CROSS" ; then
AC_CONFIG_HEADER(config/config.h)
fi

# Makefile is common and does not allow multiple builds
# Disable the the build_host folder
if test x$mh_setup != x ; then
  build_host=${host}-${mh_setup}
else
  build_host=$host
fi
#build_host=""

#------------------------------------
cat $host_frag > config/mh-config-use
mh_local_name=config/mh-${ac_sitename}-${ac_hostname}
echo "mh_local_name=$mh_local_name"
if test -f $mh_local_name ; then
  cat $mh_local_name >> config/mh-config-use
# else
#   echo "not found!"
fi
if test -f config/mh-override ; then
  cat config/mh-override >> config/mh-config-use
fi
host_frag=config/mh-config-use
#------------------------------------

AC_GET_MH($host_frag)

ac_libs=$LIBS
LIBS=""

if test "$CC" = mpcc_r ;  then
  CC=xlc_r
  FC=xlf2003_r
  F77=xlf_r
  use_mpcc_r=yes
fi

if test "x$load_modules" != x ;  then
  if test "$0" != "./conf__.tmp" ; then
    echo "$load_profile" > conf_load_module.tmp
    echo module purge >> conf_load_module.tmp
    echo module load "$load_modules" >> conf_load_module.tmp
    echo module list >> conf_load_module.tmp
    if test "x$config_target" = xblizzard ;  then
      # change configure to bash
      echo "sed '1s|/bin/sh|/bin/bash|' configure > conf__.tmp" >>  conf_load_module.tmp
    else
      cp configure conf__.tmp 
    fi
    echo "chmod +x conf__.tmp" >>  conf_load_module.tmp
    echo "./conf__.tmp $ac_configure_args" >> conf_load_module.tmp
    echo "rm conf__.tmp" >> conf_load_module.tmp
    chmod +x conf_load_module.tmp
    ./conf_load_module.tmp
    rm conf_load_module.tmp
    exit
  fi  
fi
  
dnl ######################################################################
dnl define the packages
use_components_file=config/icon_components
if test -f $use_components_file ; then
  AC_GET_MH($use_components_file)
fi
FFLAGS="${FFLAGS} ${DEFOPT}__${package}__"
CFLAGS="${CFLAGS} ${DEFCOPT}__${package}__"

for other_package in $icon_components ; do
  FFLAGS="${FFLAGS} ${DEFOPT}__${other_package}__"
  CFLAGS="${CFLAGS} ${DEFCOPT}__${other_package}__"

done
dnl ######################################################################


dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################

AC_PROG_RANLIB

dnl ######################################################################

if test -z "$CROSS" ; then
  AC_PROG_CC
fi

dnl ######################################################################
dnl Checks for LAPACK/BLAS

AC_ARG_WITH([lapack],
            [AC_HELP_STRING([--with-lapack=DIR],
                            [Alternative vendor optimized LAPACK and BLAS library])],,
                            withval=yes)

if test x$withval = xyes; then
  if test x$with_lapack != x; then
    LAPACKROOT=$with_lapack
  fi
fi

if test x$LAPACK_LIB_PATH != x; then
    LAPACK_LIB='-L$(LAPACK_LIB_PATH) '$LAPACK_LIB
elif test -n "$LAPACKROOT" ; then
    LAPACK_LIB='-L$(LAPACKROOT)/lib '$LAPACK_LIB
elif test -d lapack ; then
    LAPACKROOT="local"
    LAPACK_LIB="-L../lib -llapack -lblas"
else
    LAPACKROOT=""
    LAPACK_LIB=""
fi

dnl ######################################################################
dnl Checks for zlib.

ZLIB_LIB=""
ZLIB_INCLUDE=""

AC_ARG_WITH([zlib],
            [AC_HELP_STRING([--with-zlib=DIR],
                            [Required I/O package])],,
			    withval=yes)

if test x$withval = xyes; then
  if test x$with_zlib != x; then
    ZLIBROOT=$with_zlib
  fi
fi

if test -n "$ZLIBROOT" ; then
  if test -r "$ZLIBROOT/lib/libz.a" ; then
     ZLIB_LIB='-L$(ZLIBROOT)/lib -lz'
  else
     echo "cannot read $ZLIBROOT/lib/libz.a"  
     exit 1    
  fi
  if test -d "$ZLIBROOT/include" ; then
     ZLIB_INCLUDE='-I$(ZLIBROOT)/include'
  else
     echo "cannot find $ZLIBROOT/include"  
     exit 1    
  fi
else
  echo "Warning: the variable ZLIBROOT is not set in $host_frag"  
  echo "         set ZLIBROOT to the zlib root directory"
fi

dnl ######################################################################
dnl Checks for SZIP.

SZIP_LIB=""
SZIP_INCLUDE=""

AC_ARG_WITH([szip],
            [AC_HELP_STRING([--with-szip=DIR],
                            [Required I/O package])],,
			    withval=yes)

if test x$withval = xyes; then
  if test x$with_szip != x; then
    SZIPROOT=$with_szip
  fi
fi

if test -n "$SZIPROOT" ; then
  if test -r "$SZIPROOT/lib/libsz.a" ; then
     SZIP_LIB='-L$(SZIPROOT)/lib -lsz'
  else
     echo "cannot read $SZIPROOT/lib/libsz.a"  
     exit 1    
  fi
  if test -d "$SZIPROOT/include" ; then
     SZIP_INCLUDE='-I$(SZIPROOT)/include'
  else
     echo "cannot find $SZIPROOT/include"  
     exit 1    
  fi
else
  echo "Warning: the variable SZIPROOT is not set in $host_frag"  
  echo "         set SZIPROOT to the szip root directory"
fi

dnl ######################################################################
dnl Checks for HDF5.

HDF5_LIB=""
HDF5_INCLUDE=""

AC_ARG_WITH([hdf5],
            [AC_HELP_STRING([--with-hdf5=DIR],
                            [Required I/O package])],,
			    withval=yes)

if test x$withval = xyes; then
  if test x$with_hdf5 != x; then
    HDF5ROOT=$with_hdf5
  fi
fi

if test -n "$HDF5ROOT" ; then
  if test -r "$HDF5ROOT/lib/libhdf5.a" ; then
     HDF5_LIB='-L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5'
  else
     echo "cannot read $HDF5ROOT/lib/libhdf5.a"  
     exit 1    
  fi
  if test -d "$HDF5ROOT/include" ; then
     HDF5_INCLUDE='-I$(HDF5ROOT)/include'
  else
     echo "cannot find $HDF5ROOT/include"  
     exit 1    
  fi
else
  echo "Warning: the variable HDF5ROOT is not set in $host_frag"  
  echo "         set HDF5ROOT to the hdf5 root directory"
fi

dnl ######################################################################
dnl Checks for NETCDF.

NETCDF_LIB=""
NETCDF_INCLUDE=""

AC_ARG_WITH([netcdf],
            [AC_HELP_STRING([--with-netcdf=DIR],
                            [Required I/O package])],
                            [],
                            with_netcdf=yes)
AS_CASE([$with_netcdf],
        [no],[AC_MSG_CHECKING([for netcdf library])
              AC_MSG_RESULT([suppressed])],
        [yes],[AC_CHECK_FILE($NETCDFROOT/include/netcdf.inc,[],AC_MSG_ERROR([NETCDFROOT (currently: $NETCDFROOT) not properly defined]))],
        [*],[NETCDFROOT=$with_netcdf
             AC_CHECK_FILE($NETCDFROOT/include/netcdf.inc,[],AC_MSG_ERROR([NETCDFROOT (currently: $NETCDFROOT) not properly defined]))]
        )

if test -n "$NETCDFROOT" ; then
  if test  x$NETCDFLIBPATH = x; then
    NETCDFLIBPATH=$NETCDFROOT/lib
  fi
  if test -r "$NETCDFLIBPATH/libnetcdff.a" ; then
     NETCDF_LIB='-L$(NETCDFLIBPATH) -lnetcdff -lnetcdf'
  elif test -r "$NETCDFLIBPATH/libnetcdf.a" ; then
     NETCDF_LIB='-L$(NETCDFLIBPATH) -lnetcdf'
  elif test -r "$NETCDFLIBPATH/libnetcdff.so" ; then
     NETCDF_LIB='-L$(NETCDFLIBPATH)  -lnetcdff -lnetcdf'
  elif test -r "$NETCDFLIBPATH/libnetcdf.so" ; then
     NETCDF_LIB='-L$(NETCDFLIBPATH) -lnetcdf'
  else
     echo "cannot read $NETCDFROOT/lib/libnetcdf.a"
     exit 1
  fi
  if test x$gcc_lib != x ; then
     NETCDF_LIB="${NETCDF_LIB} ${gcc_lib}"
  fi

  if test -r "$NETCDFROOT/include/netcdf.inc" ; then
     NETCDF_INCLUDE='-I$(NETCDFROOT)/include'
  else
     echo "cannot read $NETCDFROOT/include/netcdf.inc"
     exit 1
  fi
else
  echo "error: the variable NETCDFROOT is not set in $host_frag"
  echo "       set NETCDFROOT to the netCDF root directory"
  exit 1
fi

dnl ######################################################################
dnl Checks for SCT.

SCT_LIB=""
SCT_INCLUDE=""

AC_ARG_WITH([sct],
            [AC_HELP_STRING([--with-sct=DIR],
                            [Optional SCT timer package])],
                            [],
                            with_sct=no)

AS_CASE([$with_sct],
      [no],[AC_MSG_CHECKING([for sct library])
            AC_MSG_RESULT([suppressed])],
      [yes],[AC_CHECK_FILE($SCTROOT/include/sct.mod,[],AC_MSG_ERROR([SCTROOT (currently: $SCTROOT) not properly defined to find sct.mod]))
             AC_CHECK_FILE($SCTROOT/lib/libsct.a,[],AC_MSG_ERROR([SCTROOT (currently: $SCTROOT) not properly defined to find libsct.a]))],
      [*],[SCTROOT=$with_sct
           AC_CHECK_FILE($SCTROOT/include/sct.mod,[],AC_MSG_ERROR([SCTROOT (currently: $SCTROOT) not properly defined to find sct.mod]))
           AC_CHECK_FILE($SCTROOT/lib/libsct.a,[],AC_MSG_ERROR([SCTROOT (currently: $SCTROOT) not properly defined to find libsct.a]))]
      )

if test -n "$SCTROOT" ; then
  if test -r "$SCTROOT/lib/libsct.a" ; then
    SCT_LIB='-L$(SCTROOT)/lib -lsct'
  else
    echo "cannot read $SCTROOT/lib/libsct.a"  
    exit 1
  fi
  if test -d "$SCTROOT/include" ; then
    SCT_INCLUDE='-I$(SCTROOT)/include'
  else
    echo "cannot find $SCTROOT/include"  
    exit 1
  fi

  FFLAGS="${FFLAGS} ${DEFOPT}__SCT__"
  CFLAGS="${CFLAGS} ${DEFCOPT}__SCT__"

dnl checks for additional libs needed for sct, eg. papi

  AC_ARG_WITH([sctlibs],
              [AC_HELP_STRING([--with-sctlibs=PATHS],
                              [Optional libs needed for SCT])],
                              [],
                              [])
  AS_IF([test "x$with_sctlibs" != ""],
        [SCT_LIB="${SCT_LIB} $with_sctlibs"],
        [])

else
  echo "Info: the variable SCTROOT is not set in $host_frag"  
  echo "      fallback to mo_real_timer is used"
  echo "      otherwise, set SCTROOT to the sct root directory"

  FFLAGS=$(echo ${FFLAGS} | sed "s/${DEFOPT}__SCT__//")
  CFLAGS=$(echo ${CFLAGS} | sed "s/${DEFCOPT}__SCT__//")
fi

dnl ######################################################################
dnl Checks for MPI.

case $host in
  powerpc-ibm-aix*)
    q64_result="no"
    mt_result="no"
    q64_result=$(echo $FFLAGS | awk -v b="-q64" '$0 ~ b { print "yes" }')
    mt_result=$(echo $FC | awk -v b="_r" '$0 ~ b { print "yes" } ')
    if [[ x$mt_result != xyes ]]; then
      echo "error: IBM MPI requires thread safe compiler versions."
      echo "       add _r to compiler selected by FC."
      exit 1
    fi

    if test x$q64_result = xyes; then
      MPI_INCLUDE_THREAD="thread64/"
    else
      MPI_INCLUDE_THREAD="thread/"
    fi
    ;;
  *-*-darwin*)
    if test -n "$MPIROOT" ; then
      mpifh=$(find $MPIROOT/include -name mpif.h)
      mpifh=${mpifh%\/mpif.h}
      mpifh=${mpifh##*/}
      if test -n $mpifh; then
        MPI_INCLUDE_THREAD="$mpifh/"
      fi
     fi
     ;;
  *)
    MPI_INCLUDE_THREAD=""
    ;;
esac


MPI_INCLUDE=""

withval=no
AC_ARG_WITH([mpi],
            [AC_HELP_STRING([--with-mpi=DIR],
                            [MPI library])],,)


#default is with mpi
if test x$with_mpi != xno ; then
  withval=yes
fi

if test x$withval = xno ; then
  FFLAGS="${FFLAGS} ${DEFOPT}NOMPI"
  CFLAGS="${CFLAGS} ${DEFCOPT}NOMPI"
  MPI_LIB=""
  MPI_VERSION=none
  MPIROOT=""
  use_mpcc_r=no
else
  if test x$with_mpi != x ; then
    if test $with_mpi != yes && $with_mpi != no  ; then
      MPIROOT=$with_mpi
    fi
  fi
  if test -n "$MPIROOT" ; then    
    AC_CHECK_FILE($MPIROOT/include/${MPI_INCLUDE_THREAD}mpif.h,[],AC_MSG_ERROR([MPIROOT not properly defined]))
    case $host in
    *-*-darwin*)
      MPI_LIB='-L$(MPIROOT)/lib/mpich-mp '$MPI_LIB
      ;;    
    *)
      MPI_LIB='-L$(MPIROOT)/lib '$MPI_LIB
      ;;   
    esac 
    MPI_I_INCLUDE="-I\$(MPIROOT)/include/${MPI_INCLUDE_THREAD}"
    case $host in
    powerpc-ibm-aix*)
      mod_file=$(find $MPIROOT -name mpi.mod -print 2>/dev/null | grep mpi.mod | grep ${MPI_INCLUDE_THREAD})
      ;;
    sx?-nec-superux)
      real_path=$(cd $MPIROOT; pwd -P)
      mod_file=$(find $real_path -name mpi.mod -print | grep -v modulew)
      mod_file="${MPIROOT}/lib${mod_file##*lib0}"
      ;;
    *-*-darwin*)
      real_path=$(cd $MPIROOT; pwd -P)
      mod_file=$(find $real_path/include $real_path/lib -name mpi.mod -print | awk 'NR==1')
      ;;
    *)
      real_path=$(cd $MPIROOT; pwd -P)
#      mod_file=$(find $real_path -name mpi.mod -print -quit)
      mod_file=$(find $real_path -name mpi.mod -print | awk 'NR==1')
      ;;
    esac
    if test -r $mod_file ; then
       mod_file=${mod_file#$real_path}
       mod_file=${mod_file#$MPIROOT}
       mod_file="\$(MPIROOT)"${mod_file}
       MPI_M_INCLUDE=${MODOPT}${mod_file%/*}
    else
       echo "cannot read $mod_file/mpi.mod"
       exit 1
    fi
    if test "x$MPI_M_INCLUDE" = "x$MPI_I_INCLUDE" ; then
      MPI_INCLUDE="$MPI_M_INCLUDE"
    else
      MPI_INCLUDE="$MPI_M_INCLUDE $MPI_I_INCLUDE"
    fi  
  else
    MPI="${DEFOPT}NOMPI"
    MPI_LIB=""
  fi
fi


dnl ######################################################################
dnl Checks for other flags 
PROFILE_LIB=""
PROFILE_INCLUDE=""

withval=no
AC_ARG_WITH([hpm],
            [],
            [],
             with_hpm=no)

if test x$withval = xyes; then
  FFLAGS="${FFLAGS} ${DEFOPT}__hpm__"
  PROFILE_LIB="$HPM_LIB"
  PROFILE_INCLUDE="$HPM_INCL"
fi


dnl ######################################################################
dnl Checks for header files.

if test -z "$CROSS" ; then

   AC_HEADER_STDC
   AC_CHECK_HEADERS_ONCE(fcntl.h fortran.h limits.h malloc.h netdb.h pwd.h unistd.h execinfo.h ucontext.h sys/param.h sys/time.h sys/unistd.h sys/utsname.h )

dnl Checks for typedefs, structures, and compiler characteristics.
   AC_C_CONST
   AC_HEADER_TIME

dnl Checks for library functions.
   AC_CHECK_FUNCS(getrusage gettimeofday sysconf uname valloc)

dnl Check for byte ordering

   AC_C_BIGENDIAN

dnl Check for Fortran naming convention

   AC_FC_WRAPPERS

fi # ! CROSS

test "x$exec_prefix" = xNONE && exec_prefix=build/${build_host}

my_srcdirs="support src"

if [[ "x$LAPACKROOT" = xlocal ]] ; then
    my_srcdirs="blas lapack ${my_srcdirs}"
    LAPACKROOT=""
fi


if [[ "x$use_mtime" = xyes ]] ; then
    my_srcdirs="externals/mtime/src ${my_srcdirs}"
fi

if [[ "x$use_messy" = xyes ]] ; then
    rm -f src/messy_interface
    ln -s -t src/ ../externals/messy_interface
    FFLAGS="${FFLAGS} ${DEFOPT}MESSY ${DEFOPT}MESSYTIMER"
    CFLAGS="${CFLAGS} ${DEFCOPT}MESSY ${DEFOPT}MESSYTIMER"
fi


ac_srcdirs=""
for dir in $my_srcdirs ; do
    ac_srcdirs="${ac_srcdirs} $dir"
done


if test "$use_mpcc_r" = yes ; then
  CC=mpcc_r
  FC=mpxlf2003_r
  F77=mpxlf_r
  MPIROOT=""
  MPI_INCLUDE=""
  MPI_LIB=""
fi

ac_objdirs=""
for dir in $my_srcdirs ; do
    ac_objdirs="${ac_objdirs} ${exec_prefix}/$dir"
done

#------------------------
if test x$with_openmp != x; then
    use_openmp=$with_openmp
else
    use_openmp="no"
fi

if [[ "x$target_os" = xsuperux ]] ; then  
    MODDIR=""	    
else
    MODDIR="${MODDIR}../module"
    FFLAGS="${MODOPT}../module ${FFLAGS}"
fi


if [[ "x$LDFLAGS" = "x" ]] ; then
  if [[ $(echo $FC | awk '{print index($0, "gfortran")}') -ne 0 ]] ; then  
    LDFLAGS=${FFLAGS//-xf95-cpp-input/}
  else
    LDFLAGS=$FFLAGS
  fi
fi

if [[ "x$use_openmp" = xyes ]] ; then
   FFLAGS="${FFLAGS} ${OMPFLAG}"   
   LDFLAGS="${LDFLAGS} ${OMPFLAG}"
fi

#------------------------
if test x$with_flags = x; then
    use_flags_group="std"
else
    use_flags_group=$with_flags
fi

if [[ "x$use_mtime" = xyes ]] ; then
   ac_libs="$ac_libs -lmtime"
fi

dnl the following prepends VARIANT_FLAGS to FFLAGS
dnl when given --with-flags=variant argument to configure
AS_IF([eval test x\$\{`echo ${use_flags_group} | tr a-z A-Z`_FLAGS+set\} = xset],
  [eval FFLAGS=\"\${FFLAGS} \${`echo ${use_flags_group} | tr a-z A-Z`_FLAGS}\"],
  [AS_IF([test x"$use_flags_group" != xstd],
     [AC_MSG_FAILURE([unknown with_flags parameter: $with_flags])])])


#------------------------

if test -z "$AR" ; then
  AR="ar"
fi  
if test -z "$AS" ; then
  AS="as"
fi  
if test -z "$ARFLAGS" ; then
  ARFLAGS="crv"
fi


if [[ "x$target_os" = xsuperux ]] ; then  
    arch=$(echo ${target_cpu:0:2} | tr "[:lower:]" "[:upper:]")
else
    arch=$target_cpu
fi

AC_SUBST(ARCH,$arch)
AC_SUBST(OS,$target_os)

AC_SUBST(SITE,$ac_sitename)
AC_SUBST(HOST,$ac_hostname)

AC_SUBST(NETCDFROOT,$NETCDFROOT)
AC_SUBST(NETCDFLIBPATH,$NETCDFLIBPATH)
AC_SUBST(NETCDF_LIB,$NETCDF_LIB)
AC_SUBST(NETCDF_INCLUDE,$NETCDF_INCLUDE)

AC_SUBST(HDF5ROOT,$HDF5ROOT)
AC_SUBST(HDF5_LIB,$HDF5_LIB)
AC_SUBST(HDF5_INCLUDE,$HDF5_INCLUDE)

AC_SUBST(SZIPROOT,$SZIPROOT)
AC_SUBST(SZIP_LIB,$SZIP_LIB)
AC_SUBST(SZIP_INCLUDE,$SZIP_INCLUDE)

AC_SUBST(ZLIBROOT,$ZLIBROOT)
AC_SUBST(ZLIB_LIB,$ZLIB_LIB)
AC_SUBST(ZLIB_INCLUDE,$ZLIB_INCLUDE)

AC_SUBST(SCTROOT,$SCTROOT)
AC_SUBST(SCT_LIB,$SCT_LIB)
AC_SUBST(SCT_INCLUDE,$SCT_INCLUDE)

AC_SUBST(MPIROOT,$MPIROOT)
AC_SUBST(MPI_LIB,$MPI_LIB)
AC_SUBST(MPI_INCLUDE,$MPI_INCLUDE)

AC_SUBST(LAPACKROOT,$LAPACKROOT)
AC_SUBST(LAPACK_LIB_PATH,$LAPACK_LIB_PATH)
AC_SUBST(LAPACK_LIB,$LAPACK_LIB)

AC_SUBST(METIS_LIB,$METIS_LIB)

AC_SUBST(PROFILE_LIB,$PROFILE_LIB)
AC_SUBST(PROFILE_INCLUDE,$PROFILE_INCLUDE)

AC_SUBST(AR,$AR)
AC_SUBST(AS,$AS)

AC_SUBST(CC,$CC)
AC_SUBST(FC,$FC)
AC_SUBST(F77,$F77)

AC_SUBST(ARFLAGS,$ARFLAGS)

AC_SUBST(CFLAGS,$CFLAGS)
AC_SUBST(FFLAGS,$FFLAGS)
AC_SUBST(FlibFLAGS,$FlibFLAGS)
AC_SUBST(F77FLAGS,$F77FLAGS)

AC_SUBST(LIB,$ac_libs)

AC_SUBST(MODDIR,$MODDIR)

AC_SUBST(LDFLAGS,$LDFLAGS)

AC_SUBST(SRCDIRS,$ac_srcdirs)

AC_SUBST(OBJDIRS,$ac_objdirs)


AC_OUTPUT(Makefile)

if test -n "$CROSS" ; then
   case "${host}" in	
       sx*-*-*)             ARCH=SX ;;
   esac
   cross_config=config/config.h.${ARCH}
   echo "config.status: creating config/config.h from $cross_config"
   if test ! -f $cross_config ; then
      echo "cross configuration file $cross_config does not exist!"
      exit 1
   fi   
   cp $cross_config config/config.h
fi

config/createMakefiles.pl --target=$build_host --srcdirs="${my_srcdirs}"

#-----------------------------------------------------------------------------
case $host in
  xt?-cray-*)
    f90version=$($(echo $FC -V) 2>&1 | awk 'NR==2')
    ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
    ;;
  *-*-solaris*)
    f90version=$($(echo $FC -V) 2>&1 | awk 'NR==1' | cut -f2- -d" ")
    ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
    ;;
  sx*-*-*)
    f90version=$($(echo $FC -V) 2>&1 | awk 'NR==1;NR==5' | tr -d "\n")
    ccversion=$($(echo sxc++ -V) 2>&1 | sed -e 's/.c.*$//' | sed -e 's/C++\/SX Compiler //' | sed -e 's/C\/C++ Compiler //' | awk 'NR==1;NR==3' | tr -d "\n")
    ;;
  *-ibm-aix*)
    f90version=$($(echo $FC -qversion) 2>&1 | tr -d "\n" | sed 's/Version: / \(/' | sed 's/$/\)/')
    ccversion=$($(echo $CC -qversion) 2>&1 | tr -d "\n" | sed 's/Version: / \(/' | sed 's/$/\)/') 
    ;;
  *-*-linux-*)
    case $fortran_compiler in
    default|gcc)
      f90version=$($(echo $FC --version) 2>&1 | awk 'NR==1')
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
    nag)
      f90version=$($(echo $FC -V) 2>&1 | awk 'NR==1')
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
    intel)
      f90version=$($(echo $FC -V) 2>&1 | awk 'NR==1')
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
    path)
      f90version=$($(echo $FC -v) 2>&1 | awk 'NR==1')
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
    pgi)
      f90version=$($(echo $FC -V) 2>&1 | awk 'NR==2')
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
    sun)
      f90version=$($(echo $FC -V) 2>&1 | awk 'NR==1' | cut -f2- -d" ")
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
    esac
    ;;
  *-*-darwin*)
      f90version=$($(echo $FC --version) 2>&1 | awk 'NR==1')
      ccversion=$($(echo $CC --version) 2>&1 | awk 'NR==1')
      ;;
  *)
    f90version="unknown"
    ccversion="unknown"
    ;;
esac
#-----------------------------------------------------------------------------
echo ""
echo "Selected compiler:"
echo ""
echo "    Fortran: ${f90version}"
echo "    C      : ${ccversion}"
echo ""
echo "    OpenMP : ${use_openmp}"
echo ""
#-----------------------------------------------------------------------------
if test x$MPIROOT = x || test x$CROSS = xyes; then
case $host in
powerpc-ibm-aix*)
  if test "$use_mpcc_r" = yes ; then
    mpi_version_number=$(lslpp -l ppe.poe | grep ppe.poe | uniq | awk '{print $2}')
    MPI_VERSION="IBM MPI ${mpi_version_number}"
  fi
  ;;
sx?-nec-superux)
  mpi_version_number=$(strings $MPIROOT/lib/libmpi.a | grep "MPI/SX: library Version" | awk '{print $4}')
  MPI_VERSION="SX MPI ${mpi_version_number}"
  ;;
*)
  MPI_VERSION=none
  ;;
esac
else
case $host in
*-*-darwin*)
  mpi_version_number=$(mpiexec --version | awk 'NR==2{print $2}')
  MPI_VERSION="mpich ${mpi_version_number} (MacPorts)"
  ;;		       
*)
  MPI_VERSION=${MPIROOT##*/}
esac
fi  
#-----------------------------------------------------------------------------
echo "Selected libraries:"
echo ""
echo "    MPI    : ${MPI_VERSION}"
echo "             installation  ${MPIROOT}"
if test x$NETCDFROOT = x; then
echo "    netCDF : none (error)"
else
netcdf_version=$(${NETCDFROOT}/bin/ncdump 2>&1 | grep "netcdf library version" | cut -f4 -d" ")
echo "    NETCDF : netcdf-${netcdf_version}"
echo "             installation  ${NETCDFROOT}"
fi  
if test x$with_netcdf4 != xno ; then
if test x$HDF5ROOT = x; then
echo "    HDF5   : none (error)"
else
hdf5_version=$(grep "define H5_VERSION" ${HDF5ROOT}/include/H5pubconf.h | cut -f2 -d"\"")
echo "    HDF5   : hdf5-${hdf5_version}"
echo "             installation  ${HDF5ROOT}"
fi  
if test x$SZIPROOT = x; then
echo "    SZIP   : none (error)"
else
szlib_version=$(grep "define SZLIB_VERSION" ${SZIPROOT}/include/szlib.h | cut -f2 -d"\"")
echo "    SZIP   : szip-${szlib_version}"
echo "             installation  ${SZIPROOT}"
fi  
if test x$ZLIBROOT = x; then
echo "    ZLIB   : none (error)"
else
zlib_version=$(grep "define ZLIB_VERSION" ${ZLIBROOT}/include/zlib.h | cut -f2 -d"\"")
echo "    ZLIB   : zlib-${zlib_version}"
echo "             installation  ${ZLIBROOT}"
fi
if test x$SCTROOT = x; then
echo "    SCT    : none"
else
sct_version=$(grep "define SCT_VERSION" ${SCTROOT}/include/sct.h)
echo "    SCT    : installation  ${SCTROOT}"
fi
fi
#############################
echo ""

if test "x$config_target" = "x"; then
  config_target="default"
fi

if test "x$with_mpi" = "xno" ; then
  use_mpi=no
else
  use_mpi=yes
fi

if test "x$load_module" != "x" ; then
  if test "$use_mpi" = "yes" ; then
    load_modules="$load_modules $load_mpi_modules"
  else
    load_modules="$load_modules $load_nompi_modules"
  fi
fi

# general command for running a progmam
if test "x$mpi_startrun" = "x" ; then
  if  test "x$MPIROOT" = "x" ; then
    mpi_startrun="mpiexec -n \\\$mpi_total_procs"
  else
    mpi_startrun="$MPIROOT/bin/mpiexec -n \\\$mpi_total_procs"
  fi
fi

sync_submit=${sync_submit:=$submit}

echo Creating set-up.info...
cat > config/set-up.info << EOF
#############################
# set-up.info created by configure
# do not modify this file
#############################
use_site="${ac_sitename}"
use_target="${config_target}"
use_compiler="${config_compiler}"
use_mpi="${use_mpi}"
use_openmp="${use_openmp}"
use_builddir="build/${build_host}"
use_mpi_root="${MPIROOT}"
use_mpi_procs_pernode="${use_mpi_procs_pernode}"
use_mpi_startrun="${mpi_startrun}"
use_load_modules="$load_modules"
use_load_profile="$load_profile"
use_netcdflibpath="$NETCDFLIBPATH"
use_make_command="$make_command"
use_submit="$submit"
use_sync_submit="$sync_submit"
use_shell="$use_shell"
use_icon_components="$icon_components"
#############################
EOF

echo Creating the build_command...
echo "$load_profile" > build_command
if test "x$load_modules" != x; then
   echo module purge >> build_command
   echo module load "$load_modules" >> build_command
   echo module list >> build_command
fi
if test "x${config_target}" = "xblizzard" ; then
   echo "export LDR_CNTRL=MAXDATA=0xD0000000@DSA"  >> build_command
fi
echo "$make_command" >> build_command
chmod +x build_command

#echo Creating the run scripts...
#if test -f ./make_runscripts ; then
#  ./make_runscripts
#fi

#############################
echo ""
cat config/set-up.info
echo ""
#############################

echo ""
echo "Configuration done. Run GNU make to build ICON."
echo ""
echo "build directory: build/${build_host}"
echo ""
