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
        configureoption="--with-mpi=no --with-openmp=no"
        ;;
    blizz_nMyO)
        configureoption="--with-mpi=no --with-openmp=yes"
        ;;
    blizz_yMnO)
        configureoption="--with-openmp=no"
        ;;
    blizz_yMyO)
        configureoption="--with-openmp=yes"
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

    # Thunder
    # ------
    thunder | thunder_gcc)
        configureoption="--with-fortran=gcc"
        ;;
    thunder_nag)
        configureoption="--with-fortran=nag"
        ;;
    thunder_intel)
        # configureoption="--with-fortran=intel --with-mpi --with-openmp --with-flags=hiopt"
        configureoption="--with-fortran=intel --with-mpi --with-openmp"
        ;;
    thunder_pgi)
        configureoption="--with-fortran=pgi"
        ;;
    thunder_sun)
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
        configureoption="--host=sx9 --with-setup=sx9 --without-ocean"
        ;;

    sx9omp)
        configureoption="--host=sx9 --with-setup=sx9omp --without-ocean"
        ;;

    sx9mpiomp)
        configureoption="--host=sx9 --with-setup=sx9mpiomp --without-ocean"
        ;;

    sx9ftromp)
        configureoption="--host=sx9 --with-setup=sx9ftromp --without-ocean"
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
