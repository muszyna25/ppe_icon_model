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

. ../setting

echo "WORKING_PATH ${WORKING_PATH}"
SQURCE_PATH=${WORKING_PATH}/experiments/$1/plots

TARGET_PATH=`dirname ${WORKING_PATH}`
TARGET_PATH=`basename ${TARGET_PATH}`

echo "SQURCE_DIR $SQURCE_PATH"
echo "TARGET_DIR $TARGET_PATH"

cd ${SQURCE_PATH}
scp *.eps mpipc91.mpi.zmaw.de:/scratch/local1/m211098/master_ICON/public_html/model_compare/eps/${TARGET_PATH}/.
ssh mpipc91.mpi.zmaw.de /scratch/local1/m211098/master_ICON/public_html/model_compare/convert2jpg.ksh ${TARGET_PATH}