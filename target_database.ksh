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

    # Generic
    # -------
    intel_hybrid)
        configureoption="--with-fortran=intel --with-openmp --with-yac"
        ;;
    intel)
        configureoption="--with-fortran=intel"
        ;;
    intel_openmp)
        configureoption="--with-fortran=intel --without-mpi --with-openmp"
        ;;
    cce)
        configureoption="--with-fortran=cray"
        ;;
    gcc)
        configureoption="--with-fortran=gcc"
        ;;
    nag)
        configureoption="--with-fortran=nag --with-yac"
        ;;
    nag_serial)
        configureoption="--with-fortran=nag --without-mpi"
        ;;
    nag_mtime)
        configureoption="--with-fortran=nag --enable-mtime-loop"
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
