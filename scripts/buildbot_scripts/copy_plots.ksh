#!/bin/ksh
# This script is used to copy the eps plot files to mpipc91.
# On mpipc91 the eps-files are converted to jpg-file.
#-------------------------------------------------------------------
#
# First version by Walter Sauf (MPI-M, 2010-02-26)
#
# $Rev$:     Revision of last commit
# $Author$:  Author of last commit
# $Date$:    Date of last commit
#-------------------------------------------------------------------

REV=$1
SLAVE=$2
BUILDER=$3
BUILDER_NR=$4
DATE=`date '+%Y-%m-%d'`


#===============================================================
# functions
#===============================================================

copy_files ()
{

  for file in $FILES
  do
    dir_name=`echo $file | cut -d '/' -f2-10`
    name=`basename $dir_name`
    dir=`dirname $dir_name`
#    echo "DIR: $dir NAME: $name"
    mkdir -p ${BASE_DIR}/${dir}
#    echo "Copy:"
#    echo "Source: experiments/${dir}/${name}"
#    echo "Target: ${BASE_DIR}/${dir}/${name}"
    cp experiments/${dir}/${name} ${BASE_DIR}/${dir}/${name}
  done
}

#===============================================================
# Script begin
#===============================================================

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Save Plots for Buildbot web-Page
#---------------------------------------------------------------------
#---------------------------------------------------------------------

if [ "x${BUILDER}" = "x" ]
then
  echo "!!!! BUILDER not set "
  exit 1
fi

if [ -d /tmp/${BUILDER} ]
then
  rm -rf /tmp/${BUILDER} 
fi

echo "BB_SYSTEM=${SLAVE}"
echo "BB_SLAVE=${BUILDER}"

#REV=`svn info | grep Revision | cut -d ':' -f2`
echo "REV=${REV}"

mkdir /tmp/${BUILDER} 
echo "_COMPUTER_ ${SLAVE}" > /tmp/${BUILDER}/job_info.txt
echo "_BUILDER_ ${BUILDER}" >> /tmp/${BUILDER}/job_info.txt
echo "_REVISION_ ${REV}" >> /tmp/${BUILDER}/job_info.txt

find experiments -name '*.eps' -exec cp {} /tmp/${BUILDER}/. \;

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Save Plots for archive
#---------------------------------------------------------------------
#---------------------------------------------------------------------

BASE_DIR=/tmp/BuildBot/${BUILDER}/archive/${DATE}/buildbot/${BUILDER}/${BUILDER_NR}

if [ -d /tmp/BuildBot/${BUILDER}/archive ]
then
  rm -rf /tmp/BuildBot/${BUILDER}/archive
fi

mkdir -p ${BASE_DIR}

echo "Copy eps-files"
FILES=`find experiments -name '*.eps'`
copy_files

echo "Copy ps-files"
FILES=`find experiments -name '*.ps'`
copy_files
