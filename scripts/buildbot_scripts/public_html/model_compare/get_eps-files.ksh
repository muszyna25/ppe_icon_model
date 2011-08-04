#!/bin/ksh
. /client/etc/profile.zmaw
module load svn
UPLOAD_DIR=/scratch/local1/m211098/master_ICON/upload

#=======================================================
#=======================================================


check_new_dir()
{
  DIRS=`ls /scratch/local1/m211098/master_ICON/upload/eps`
  for DIR in $DIRS
  do
    rm -rf ${UPLOAD_DIR}/eps_use/$DIR
    mv ${UPLOAD_DIR}/eps/$DIR ${UPLOAD_DIR}/eps_use/$DIR
  done
}

#=======================================================
#=======================================================

#
TARGET_PATH=/scratch/local1/m211098/master_ICON/public_html/model_compare
echo "=============  Begin  ==========================="
echo " "
date
echo " "
echo "================================================="

cd ${TARGET_PATH}

# This call is done at the begin of the script to save the jpg-file for the yesterday
# run. At the end when the new eps plots have been loaded and the eps-files has been 
# created out of the eps-files we are caolling the script ones again

${TARGET_PATH}/run_html.ksh

echo "STATUS: $?"

# Check if a new upload directory is there

check_new_dir

cp -r ${UPLOAD_DIR}/eps_use/* eps/.

DIRS=`ls ${UPLOAD_DIR}/eps_use`
for DIR in $DIRS
do
#    echo "DIR: $DIR"
  TARGET_DIR=`basename $DIR`
  echo "TARGET_DIR $TARGET_DIR"
  if [ -d eps/${TARGET_DIR} ]
  then
#     delete all old jpg files
    find eps/${TARGET_DIR} -name '*.jpg' -exec rm -f {} \; 2> /dev/null
  else
    echo "Create directory: eps/$TARGET_DIR"
    mkdir eps/$TARGET_DIR
  fi 



done

${TARGET_PATH}/convert2jpg.ksh
${TARGET_PATH}/run_html.ksh

echo "=============   End   ==========================="
echo " "
date
echo " "
echo "================================================="
