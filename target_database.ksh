#==============================================================================
#
# Database for target specific configuration parameters
# -----------------------------------------------------
#
# The database contains parameters which are used for the configuration
# of the Makefile. These parameters are defined for "targets", which are
# descriptive names referring typically to the computer system, the compiler
# or the model parallelization to be used.
#
#==============================================================================
#
# Parameters
# ---------- 
#
#   configureoption = option(s) for configure command
#
#------------------------------------------------------------------------------

db_status=0 # exit status of database, 0=ok, 1=target not known

# for specific targets
# --------------------
case ${target} in

    # Blizzard
    # --------
    blizz_nMnO)
        configureoption="--with-mpi=no  --with-openmp=no --with-flags=hiopt"
        ;;
    blizz_nMyO)
        configureoption="--with-mpi=no  --with-openmp=yes --with-flags=hiopt"
        ;;
    blizz_yMnO)
        configureoption="--with-openmp=no --with-flags=hiopt"
        ;;
    blizz_yMyO)
        configureoption="--with-openmp=yes --with-flags=hiopt"
        ;;

    # MPIPC
    # -----
    mpipc | mpipc_gcc)
        configureoption="--with-fortran=gcc"
        ;;
    mpipc_nag)
        configureoption="--with-fortran=nag"
        ;;
    mpipc_intel)
        configureoption="--with-fortran=intel"
        ;;
    mpipc_pgi)
        configureoption="--with-fortran=pgi"
        ;;
    mpipc_sun)
        configureoption="--with-fortran=sun"
        ;;

    # Squall
    # ------
    squall | squall_gcc)
        configureoption="--with-fortran=gcc"
        ;;
    squall_nag)
        configureoption="--with-fortran=nag"
        ;;
    squall_intel)
        configureoption="--with-fortran=intel"
        ;;
    squall_pgi)
        configureoption="--with-fortran=pgi"
        ;;
    squall_sun)
       configureoption="--with-fortran=sun"
        ;;

    # Tornado
    # -------
    tornado | tornado_gcc)
        configureoption="--with-fortran=gcc --with-loopexchange=yes"
        ;;
    tornado_nag)
        configureoption="--with-fortran=nag"
        ;;
    tornado_intel)
        configureoption="--with-fortran=intel"
        ;;
    tornado_pgi)
        configureoption="--with-fortran=pgi"
        ;;
    tornado_sun)
        configureoption="--with-fortran=sun"
        ;;

    # MPIMAC
    # ------
    mpimac | mpimac_gcc)
        configureoption="--with-fortran=gcc"
        ;;

    # DWD hpc 
    # -------
    hpc | hpc_noomp)
       configureoption="--with-fortran=sun-noomp"
        ;;
    hpc_sun_debug)
       configureoption="--with-fortran=sun-debug --without-mpi"
        ;;
    hpc_serial )
        configureoption="--with-fortran=sun --without-mpi"
        ;;

    # DWD SX9
    # -------
    sx9)
        configureoption="--host=sx9 --with-setup=sx9"
        ;;

    sx9omp)
        configureoption="--host=sx9 --with-setup=sx9omp"
        ;;

    sx9mpiomp)
        configureoption="--host=sx9 --with-setup=sx9mpiomp"
        ;;

    sx9ftromp)
        configureoption="--host=sx9 --with-setup=sx9ftromp"
        ;;

    # Default (for unspecified target)
    # -------
    default)
        configureoption="--with-mpi=no"
        ;;

    # Error trap, in case $target has unknown value
    *)
        db_status=1
        ;;

esac

if   [ "$db_status" == "0" ]
then
    echo "Parameters set for target=$target"
    echo "  configureoption = $configureoption"
else
    echo "No parameters set for target=$target --> exit"
    exit $db_status
fi
