#!/bin/ksh

TARGET_PATH=/scratch/local1/m211098/master_ICON/public_html/model_compare
EXPs="EXPs_1 EXPs_2 EXPs_3 EXPs_4 EXPs_5 EXPs_6 EXPs_7"
#========================================================================
clean_eps_dir(){
  # Cleanup the eps-directory 
  DIRS=`ls eps`
  for DIR in $DIRS
  do
    if [ "$DIR" != "reference" ]
    then
      find eps/$DIR -name '*.jpg' -exec rm -f {} \;
    fi
  done
}
#========================================================================
switch_today2yesterday(){
   echo "!!!!! Switch Today to Yesterday"

# # Delete plots from day before yesterday
#     rm -rf eps_y/*
# 
# # rename the directory
#     cp -r eps/* eps_y/.
 
  for name in $NAMES
  do
# build the new html-page for yesterday data
    echo "name $name"
    grep -v 'Yesterday Run' compare_${name}.html > compare_${name}_y.html
    sed -i 's/eps/eps_y/g' compare_${name}_y.html
  done
# Cleanup the eps-directory 
  DIRS=`ls eps | grep -v "reference_"`
  echo "DIRS $DIRS"
  for Dir in $DIRS
  do
    echo "Dir $Dir"
    find eps/$Dir -name '*.jpg' -exec rm -f {} \;
  done

  echo "$DATE" > .date

}
#========================================================================
build_html_file(){

# Get a list of all directories in eps
  DIRS=`ls eps`
#  echo "DIRS: $DIRS"
  for name in $NAMES
  do
   echo "!!!!! Build a HTML-file: $name"

# delete old file
    rm -f compare.html

# create new html-page from compare_template.html and include the actuale date
    sed 's/_DATUM_/'$DATE'/g'  compare_template.html > compare.html
    sed -i 's/_EXP_/'$EXP'/g'  compare.html
    sed -i 's/_UPDATE_/'"$UPDATE"'/g'  compare.html

    sed -i 's/ Missing <\/td> <!-- reference -->/ Missing <\/td> <!-- reference_'"${name}"' -->/g' compare.html
    sed -i 's/_COMPUTER_reference/_COMPUTER_reference_'"${name}"'/g'  compare.html
    sed -i 's/_BUILDER_reference/_BUILDER_reference_'"${name}"'/g'  compare.html
    sed -i 's/_REVISION_reference/_REVISION_reference_'"${name}"'/g'  compare.html

    sed -i 's/_PUPDATE_reference/_PUPDATE_reference_'"${name}"'/g'  compare.html
    sed -i 's/_BUILDID_reference/_BUILDID_reference_'"${name}"'/g'  compare.html

    PLOT_FILE=${EXP}_${name}.jpg

# build the list of directors included in eps
    for Dir in $DIRS
    do 
# Chech if the specified plotfile exist
      if [ -a eps/$Dir/${PLOT_FILE} ]
      then
        echo "Found file eps/$Dir/${PLOT_FILE}"
#   replace the Missing info by including the imgage link
        sed -i 's/ Missing <\/td> <!-- '${Dir}' -->/ <img src="eps\/'${Dir}'\/'${PLOT_FILE}'"> <\/td>/g' compare.html
      else
        pwd
        echo "Could not find file eps/$Dir/${PLOT_FILE}"
      fi
      COMPUTER_NAME=`grep _COMPUTER_ eps/${Dir}/job_info.txt | cut -d ' ' -f2`
      BUILDER_NAME=`grep _BUILDER_ eps/${Dir}/job_info.txt | cut -d ' ' -f2`
      REVISION_NR=`grep _REVISION_ eps/${Dir}/job_info.txt | cut -d ' ' -f2`
      PLOT_UPDATE_D=`grep _UPDATE_ eps/${Dir}/job_info.txt | cut -d ' ' -f2`
      PLOT_UPDATE_T=`grep _UPDATE_ eps/${Dir}/job_info.txt | cut -d ' ' -f3`
      BUILD_ID=`grep _BUILDID_ eps/${Dir}/job_info.txt | cut -d ' ' -f2`
# 
      sed -i 's/_COMPUTER_'${Dir}'#/'${COMPUTER_NAME}'/g' compare.html
      sed -i 's/_BUILDER_'${Dir}'#/'${BUILDER_NAME}'/g' compare.html
      sed -i 's/_REVISION_'${Dir}'#/'${REVISION_NR}'/g' compare.html
      sed -i 's/_PUPDATE_'${Dir}'#/'${PLOT_UPDATE_D}'\ '${PLOT_UPDATE_T}'/g' compare.html
      sed -i 's/_BUILDID_'${Dir}'#/'${BUILD_ID}'/g' compare.html
# 
    done
   sed -i 's/compare_yesterday/compare_'${name}'_y/g' compare.html
     cp compare.html compare_${name}.html
  done
}

#========================================================================

create_html(){

  if  [ "$SAVE_DATE" != "$DATE" ]
  then
    switch_day=true
# Delete plots from day before yesterday
    day_of_month=`date --date=yesterday '+%d'`
    rm -f eps_${day_of_month}.tgz
    tar -czf eps_${day_of_month}.tgz eps
    rm -rf eps_y/*

# rename the directory
    cp -r eps/* eps_y/.
  fi

  if  [ "${switch_day}" == "true" ]
  then
    switch_today2yesterday
  fi

  build_html_file
}

#========================================================================



# Change to actual working directory
cd ${TARGET_PATH}

# Get the actuel date
DATE=`date '+%d.%m.%Y'`

UPDATE=`date '+%d.%m.%Y %k:%M'`

# Get the date from last call of this script
SAVE_DATE=`cat .date`

switch_day=false

for exp in $EXPs 
do 
   EXP=`cat ${exp} | cut -d ' ' -f1`
   NAMES=`cat ${exp} | cut -d ' ' -f2-`

   create_html
   SAVE_DATE=`cat .date`

done

exit 0
