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
BB_SYSTEM=$2
BB_SLAVE=$3

if [ "x${BB_SLAVE}" = "x" ]
then
  echo "!!!! BB_SLAVE not set "
  exit 1
fi

if [ -d /tmp/${BB_SLAVE} ]
then
  rm -rf /tmp/${BB_SLAVE} 
fi

echo "BB_SYSTEM=${BB_SYSTEM}"
echo "BB_SLAVE=${BB_SLAVE}"

#REV=`svn info | grep Revision | cut -d ':' -f2`
echo "REV=${REV}"

mkdir /tmp/${BB_SLAVE} 
echo "_COMPUTER_ ${BB_SYSTEM}" > /tmp/${BB_SLAVE}/job_info.txt
echo "_BUILDER_ ${BB_SLAVE}" >> /tmp/${BB_SLAVE}/job_info.txt
echo "_REVISION_ ${REV}" >> /tmp/${BB_SLAVE}/job_info.txt

find experiments -name '*.eps' -exec cp {} /tmp/${BB_SLAVE}/. \;
