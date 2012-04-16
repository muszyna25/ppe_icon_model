#!/bin/ksh
#==============================================================================
#
# This script configures and makes the ICON executables for a specific target. 
#
# The executables are made for the build target given by 
# the variable $target. Its value is determined as follows:
#
#   if the target name is given as argument, then:
#     target=$1
#   else:
#     target="default" 
#
# If the value of $target exists in the target database, then variables 
# used for the configuration are set, and the script proceeds with the 
# configure&make process.
#
# If the value of $target does not exist in the case table, 
# then the script stops.
#
# Marco Giorgetta, MPI-M, 2010-04-12
#
#==============================================================================
#==============================================================================
#
# Define useful functions to execute commands and check the exit status.
#
# Walter Sauf,     MPI-M, 2009-10-22
# Marco Giorgetta, MPI-M, 2010-04-12
#
#==============================================================================

check_error()
{
    # Check if the first parameter (return status) is not OK
    # Arguments:
    #   $1 = error status: 0 = OK, not 0 = ERROR
    #   $2 = error message

    if [ "${STATUS_FILE}" != "" ] 
    then
        echo "$1" > ${STATUS_FILE}
    fi

    if [ $1 != 0 ] 
    then
        echo "check_error()"
        echo "   ERROR : $2"

        exit $1
    fi
}

#------------------------------------------------------------------------------

return_ok()
{
    # Exit with status = 0 = OK
    # Arguments:
    #   $1 =  message

    echo "return_ok()"
    echo "$1"

    exit 0
}

#------------------------------------------------------------------------------

exec_and_check()
{
    # Execute command and check for error
    # Arguments:
    #   $1 = command to be executed

    $1

    check_error $? "exec_and_check(): $1"
}

#==============================================================================
# -------------------------
# Executable part of script
# -------------------------

scriptname="target_confmake.ksh"

echo "$scriptname: start"

#------------------------------------------------------------------------------
# Set target variable, which is used to select case specific parameters from a
# built in data base.

if   [ "$1" != "" ]          # --> from argument
then
    target="$1"
else                         # --> as default
    target="default"
fi

echo "$scriptname: target = $target"


# Get target specific parameters from database
#
#   loadmodule      = name of Fortran compiler module to be loaded
#   configureoption = option(s) for configure command

. ./target_database.ksh

#------------------------------------------------------------------------------
# Configure Makefile ...
exec_and_check "./configure $configureoption"
# clean...
exec_and_check "make clean"
# ... and make executables
#exec_and_check "$make_command"
#compilation on blizzard is too slow when called through exec_and_check
# include the build_command
. ./build_command
#------------------------------------------------------------------------------
error_status=$?

#------------------------------------------------------------------------------
# this is done in buildbot
#./make_runscripts
#------------------------------------------------------------------------------

if [ "${STATUS_FILE}" != "" ]
then
  echo "$error_status" > ${STATUS_FILE}
fi

if [ $error_status != 0 ]
then
  echo "check_error()"
  echo "   ERROR : build_command"
  exit $error_status
fi

#==============================================================================

# Finish script with status OK

return_ok "$scriptname: end"
