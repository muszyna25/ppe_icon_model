#!/bin/ksh

TARGET_DIR=/scratch/local1/m211098/master_ICON/public_html/model_compare
cd ${TARGET_DIR}
pwd

DIRS=`ls eps`
echo "DIRS $DIRS"

for DIR in $DIRS
do
  ./convert2jpg.ksh eps/$DIR
done

