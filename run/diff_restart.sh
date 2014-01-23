#!/bin/ksh
MODEL_DIR=$1
EXP1=$2
EXP2=$3
TYPES="atm_phy atm_dyn lnd_phy"
DATES='040 045 050 055 100'
for TYPE in $TYPES; do 
cd ${MODEL_DIR}/experiments/${EXP1}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP1}'_'${TYPE}'_000'${DATE}'00.nc'
done
cdo -O mergetime $FILES ${EXP1}_${TYPE}.nc
cd ${MODEL_DIR}/experiments/${EXP2}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP2}'_'${TYPE}'_000'${DATE}'00.nc'
done
cdo -O mergetime $FILES ${EXP2}_${TYPE}.nc
cdo -O sub ${EXP2}_${TYPE}.nc ../${EXP1}/${EXP1}_${TYPE}.nc diff_${EXP2}-${EXP1}_${TYPE}.nc
/scratch/local2/m218036/icon/t63_remap.sh  diff_${EXP2}-${EXP1}_${TYPE}.nc  diff_${EXP}_${TYPE}_t63.nc
cdo diff ${EXP2}_${TYPE}.nc ../${EXP1}/${EXP1}_${TYPE}.nc
done
