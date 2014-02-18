#!/bin/ksh
MODEL_DIR1=$1
EXP1=$2
MODEL_DIR2=$3
EXP2=$4
TEST=$5
#TYPES="DOM01_ML"
DATES='040 050 100'
EXIT_STATUS=0
#create more beautiful output for variables MODEL_DIR[12]
cd $MODEL_DIR1; MODEL_DIR1=`pwd`
cd $MODEL_DIR2; MODEL_DIR2=`pwd`
#for TYPE in $TYPES; do 
cd ${MODEL_DIR1}/experiments/${EXP1}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP1}'_20080801T0'${DATE}'00Z.nc'
done
cdo -O mergetime $FILES ${EXP1}.nc 1> /dev/null 2> /dev/null
cd ${MODEL_DIR2}/experiments/${EXP2}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP2}'_20080801T0'${DATE}'00Z.nc'
done
cdo -O mergetime $FILES ${EXP2}.nc 1> /dev/null 2> /dev/null
cdo -O sub ${EXP2}.nc ../${EXP1}/${EXP1}.nc diff_${EXP2}-${EXP1}.nc 1> /dev/null 2> /dev/null
cdo remapdis,r360x20 diff_${EXP2}-${EXP1}.nc diff_${EXP2}-${EXP1}_r360x20.nc 1> /dev/null 2> /dev/null
if [ -f test$$.dat ]; then
  rm -f test$$.dat
fi
FILE1=${MODEL_DIR1}/experiments/${EXP1}/${EXP1}.nc
FILE2=${MODEL_DIR2}/experiments/${EXP2}/${EXP2}.nc
cdo diffv ${FILE1} ${FILE2} > test$$.dat
STATUS=`echo $?`
if [ $STATUS -ne 0 -o -s test$$.dat ]; then
echo -e "\033[31m${TEST} test FAILED:\033[00m"
echo -e "The variables in files ${FILE1} and ${FILE2} differ or the files cannot be compared"
  EXIT_STATUS=$((EXIT_STATUS + 1))
fi
#done # file types (atm_phy atm_dyn lnd_phy)
if [ -s test$$.dat ]; then
mv test$$.dat diff_${EXP2}-${EXP1}.dat
echo -e "\033[31mThe differences are in file diff_${EXP2}-${EXP1}.dat\033[00m"
else
rm -f test$$.dat
fi
if [ $EXIT_STATUS == 0 ]; then
echo -e "\033[32m${TEST} o.k.:\033[00m"
echo -e "The variables in files ${MODEL_DIR1}/experiments/${EXP1}/${EXP1}.nc and ${MODEL_DIR2}/experiments/${EXP2}/${EXP2}.nc do not differ"
fi
exit $EXIT_STATUS
