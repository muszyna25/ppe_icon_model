#!/bin/ksh
LOCAL_PATH=$PWD

TARGET_DIR=/scratch/local1/m211098/master_ICON/public_html/model_compare/eps
DIRS=`ls ${TARGET_DIR}`

# echo "TARGET_DIR: $TARGET_DIR"

for DIR in $DIRS
do
  cd ${TARGET_DIR}/$DIR
  pwd
  FILES=`ls *.eps 2> /dev/null`
  basedir=${dir%/*}

  if [ "$FILES" != "" ]
  then
    rm *.jpg 2> /dev/null
    for FILE in $FILES
    do
      baseName=`echo $FILE | cut -d '.' -f1` 
#      echo "baseName $baseName  FILE: $FILE"
      convert ${baseName}.eps ${baseName}.jpg
    done
  fi

  rm *.eps 2> /dev/null
done

cd ${LOCAL_PATH}

