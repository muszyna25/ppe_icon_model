#!/bin/ksh
TARGET_PATH=/scratch/local1/m211098/master_ICON/public_html/model_compare

# Change to actual working directory
cd ${TARGET_PATH}

# Get the actuel date
DATE=`date '+%d.%m.%Y'`

UPDATE=`date '+%d.%m.%Y %k:%M'`

# Get the date from last call of this script
SAVE_DATE=`cat .date`

# define the plotfile to show

EXPS="test_nh35_tri_JWs_bb"

NAMES="R2B04L31_VOR850_day09"

# Check if we are on a new day
if  [ "$SAVE_DATE" != "$DATE" ]
then
  echo "!!!!! Build new Infos"
# Delete plots from day before yesterday
  rm -rf eps_y/*

# rename the directory

  cp -r eps/* eps_y/.

# build the new html-page for yesterday data
  for NAME in $NAMES
  do 
    grep -v 'Yesterday Run' compare_${NAME}.html > compare_${NAME}_y.html
    sed -i 's/eps/eps_y/g' compare_${NAME}_y.html
  done

# Cleanup the eps-directory 
  DIRS=`ls eps`
  for DIR in $DIRS
  do
    if [ "$DIR" != "reference" ]
    then
      find eps/$DIR -name '*.jpg' -exec rm -f {} \;
    fi
  done
fi

# Save the actuale date for next script run
echo "$DATE" > .date

# delete actual html-page

DIRS=`ls eps`

# Make a run vor each directory

for NAME in $NAMES
do 
  rm -f compare.html

# create new html-page from compare_template.html and include the actuale date
  sed 's/_DATUM_/'$DATE'/g'  compare_template.html > compare.html
  sed -i 's/_UPDATE_/'"$UPDATE"'/g'  compare.html

  for EXP in $EXPS
  do 
    PLOT_FILE=${EXP}_${NAME}.jpg
  
# build the list of directors included in eps
    for DIR in $DIRS
    do 
# Chech if the specified plotfile exist
      if [ -a eps/$DIR/${PLOT_FILE} ]
      then
        echo "Found file eps/$DIR/${PLOT_FILE}"
#   replace the Missing info by including the imgage link
        sed -i 's/ Missing <\/td> <!-- '${DIR}' -->/ <img src="eps\/'${DIR}'\/'${PLOT_FILE}'"> <\/td>/g' compare.html
      fi
      COMPUTER_NAME=`grep _COMPUTER_ eps/${DIR}/job_info.txt | cut -d ' ' -f2`
      BUILDER_NAME=`grep _BUILDER_ eps/${DIR}/job_info.txt | cut -d ' ' -f2`
      REVISION_NR=`grep _REVISION_ eps/${DIR}/job_info.txt | cut -d ' ' -f2`

      sed -i 's/_COMPUTER_'${DIR}'/'${COMPUTER_NAME}'/g' compare.html
      sed -i 's/_BUILDER_'${DIR}'/'${BUILDER_NAME}'/g' compare.html
      sed -i 's/_REVISION_'${DIR}'/'${REVISION_NR}'/g' compare.html

    done
  done
  sed -i 's/compare_yesterday/compare_'${NAME}'_y/g' compare.html
  cp compare.html compare_${NAME}.html
done
