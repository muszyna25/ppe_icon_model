#!/bin/ksh
MODEL_DIR1=$1
EXP1=$2
MODEL_DIR2=$3
EXP2=$4
TEST=$5
TYPES="atm_phy atm_dyn lnd_phy"
DATES='040 045 050 055 100'
EXIT_STATUS=0
for TYPE in $TYPES; do 
cd ${MODEL_DIR1}/experiments/${EXP1}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP1}'_'${TYPE}'_19780101T0'${DATE}'00Z.nc'
done
cdo -O mergetime $FILES ${EXP1}_${TYPE}.nc 1> /dev/null 2> /dev/null
cd ${MODEL_DIR2}/experiments/${EXP2}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP2}'_'${TYPE}'_19780101T0'${DATE}'00Z.nc'
done
cdo -O mergetime $FILES ${EXP2}_${TYPE}.nc 1> /dev/null 2> /dev/null
cdo -O sub ${EXP2}_${TYPE}.nc ../${EXP1}/${EXP1}_${TYPE}.nc diff_${EXP2}-${EXP1}_${TYPE}.nc 1> /dev/null 2> /dev/null
REMAP_FILE=/pool/data/ICON/grids/private/r2b4_amip/remapdis_r2b4_amip_cell_t63grid.nc
if [ -f ${REMAP_FILE} ]; then
cdo remap,t63grid,${REMAP_FILE} diff_${EXP2}-${EXP1}_${TYPE}.nc diff_${EXP2}-${EXP1}_${TYPE}_t63.nc 1> /dev/null 2> /dev/null
fi
if [ -f test$$.dat ]; then
  rm -f test$$.dat
fi
FILE1=${MODEL_DIR1}/experiments/${EXP1}/${EXP1}_${TYPE}.nc
FILE2=${MODEL_DIR2}/experiments/${EXP2}/${EXP2}_${TYPE}.nc
cdo diff ${FILE1} ${FILE2} > test$$.dat
STATUS=`echo $?`
if [ $STATUS == 0 -a ! -s test$$.dat ]; then
echo -e "\033[32m${TEST} o.k.:\033[00m"
echo -e "The variables in files ${FILE1} and ${FILE2} do not differ"
else
echo -e "\033[31m${TEST} test FAILED:\033[00m"
echo -e "The variables in files ${FILE1} and ${FILE2} differ"
  EXIT_STATUS=$((EXIT_STATUS + 1))
fi
done # file types (atm_phy atm_dyn lnd_phy)
rm -f test$$.dat
exit $EXIT_STATUS
