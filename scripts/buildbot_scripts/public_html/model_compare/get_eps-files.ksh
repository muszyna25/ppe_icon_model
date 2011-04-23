#!/bin/ksh
. /client/etc/profile.zmaw
module load svn

MACHINES="mpipc22 mpipc91 mpipc56 blizzard tornado squall"
#
TARGET_PATH=/scratch/local1/m211098/master_ICON/public_html/model_compare
echo "=============  Begin  ==========================="
echo " "
date
echo " "
echo "================================================="

cd ${TARGET_PATH}

${TARGET_PATH}/run_html.ksh
echo "STATUS: $?"

for MACHINE in $MACHINES
do
  DIRS=`cat dir.$MACHINE`
  for DIR in $DIRS
  do
#    echo "DIR: $DIR"
    TARGET_DIR=`basename $DIR`
    if [ -d eps/${TARGET_DIR} ]
    then
#     delete all old jpg files
      find eps/${TARGET_DIR} -name '*.jpg' -exec rm -f {} \; 2> /dev/null
    else
      echo "Create directory: eps/$TARGET_DIR"
      mkdir eps/$TARGET_DIR
    fi 

#    echo "TARGET_DIR: ${TARGET_DIR}"
    REMOTE_LISTS=`ssh $MACHINE ls $DIR/build/experiments 2> /dev/null`
    if [ "$REMOTE_LISTS" != "" ]
    then
#      echo "REMOTE_LISTS: $REMOTE_LISTS"
      for REMOTE_LIST in $REMOTE_LISTS
      do
#        echo "scp -r $MACHINE:$DIR/build/experiments/${REMOTE_LIST}/plots /tmp/experiments_$$"

        scp -r $MACHINE:$DIR/build/experiments/${REMOTE_LIST}/plots /tmp/experiments_$$ 2> /dev/null
        if [ $? -eq 0 ]
        then 
          cp /tmp/experiments_$$/*.eps ${TARGET_PATH}/eps/${TARGET_DIR}/. 2> /dev/null
        fi 
        rm -rf /tmp/experiments_$$
        echo "Transfer of $REMOTE_LIST from $MACHINE is done. TARGET_DIR: ${TARGET_DIR}"
      done
    else
      echo "!!!! REMOTE_LISTS of machine $MACHINE (TARGET_DIR: ${TARGET_DIR}) is emty"

    fi

    ssh $MACHINE svn info $DIR/build > .tmp
    REV_NR=`grep Revision .tmp | cut -d ' ' -f2`

    sed 's/_COMPUTER_NAME_/'$MACHINE'/g' job_info_template.txt > ${TARGET_PATH}/eps/${TARGET_DIR}/job_info.txt
    sed -i 's/_BUILDER_NAME_/'${TARGET_DIR}'/g' ${TARGET_PATH}/eps/${TARGET_DIR}/job_info.txt
    sed -i 's/_REV_NR_/'${REV_NR}'/g' ${TARGET_PATH}/eps/${TARGET_DIR}/job_info.txt
  done
done
${TARGET_PATH}/convert2jpg.ksh
${TARGET_PATH}/run_html.ksh

echo "=============   End   ==========================="
echo " "
date
echo " "
echo "================================================="
