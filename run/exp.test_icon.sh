#!/bin/ksh
#################### functions ##########################################
function copy_experiment
{
MODEL_DIR=$1
EXP_ORI=$2
EXP_NEW=$3
cd ${MODEL_DIR}/experiments
if [ ! -d ${EXP_ORI} ]; then
  echo 'no experiment '${EXP_ORI}' found'
  exit 1
fi
if [ ! -f ${EXP_ORI}/restart.info -o ! -f ${EXP_ORI}/${EXP_ORI}_restart_atm_${RESTART_DATE}.nc ]; then
  echo -e "\033[31mThe run ${EXP_ORI} did not write any restart at ${RESTART_DATE}\033[00m"
  exit 1
fi
rm -rf ${EXP_NEW}
mkdir ${EXP_NEW}
cd ${EXP_NEW}
cp ../${EXP_ORI}/restart.info .
cp ../${EXP_ORI}/${EXP_ORI}_restart_atm_${RESTART_DATE}.nc .
ln -sf ${EXP_ORI}_restart_atm_${RESTART_DATE}.nc restart_atm_DOM01.nc
}
#------------------------------------------------------------------------
function print_usage
{
COMMAND=$1
echo "usage: ${COMMAND} [-e] [-h] [-m base|update|restart|nproma|ur|un|rn|urn] [-o yes|no] [-r <reference model path>]" 
echo '-e : experiment:'
echo '     r2b4_amip (hydrostatic amip experiment)'
echo '     mtest_nat_rce_cbl_120km_nwp (non-hydrostatic, radiative-convective equilibrium, nwp physics, torus grid'
echo '-h : display help'
echo '-m : test mode, either single tests like base, update, restart or nproma test or'
echo '     combined tests like ur (update and restart) possible'
echo '-o : ''yes'': overwrite existing experiments, that is the default'
echo '     ''no'': use existing experiments'
echo '-r : reference model path'
}
#------------------------------------------------------------------------
function diff_results
{
MODEL_DIR1=$1
EXP1=$2
MODEL_DIR2=$3
EXP2=$4
TEST=$5
#TYPES: defined in main script and used from main script 
#DATES: defined in main script and used from main script 
DIFF_STATUS=0
#create more beautiful output for variables MODEL_DIR[12]
cd $MODEL_DIR1; MODEL_DIR1=`pwd`
cd $MODEL_DIR2; MODEL_DIR2=`pwd`
for TYPE in $TYPES; do 
cd ${MODEL_DIR1}/experiments/${EXP1}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP1}'_'${TYPE}'_'${DATE}'.nc'
done
cdo -O mergetime $FILES ${EXP1}_${TYPE}.nc 1> /dev/null 2> /dev/null
cd ${MODEL_DIR2}/experiments/${EXP2}
FILES=''
for DATE in $DATES; do
FILES=$FILES' '${EXP2}'_'${TYPE}'_'${DATE}'.nc'
done
cdo -O mergetime $FILES ${EXP2}_${TYPE}.nc 1> /dev/null 2> /dev/null
cdo -O sub ${EXP2}_${TYPE}.nc ../${EXP1}/${EXP1}_${TYPE}.nc diff_${EXP2}-${EXP1}_${TYPE}.nc 1> /dev/null 2> /dev/null
if [ -z "${REMAP_FILE}" ]; then
cdo remapdis,${RESOL} diff_${EXP2}-${EXP1}_${TYPE}.nc diff_${EXP2}-${EXP1}_${TYPE}_${RESOL}.nc 1> /dev/null 2> /dev/null
else
if [ -f ${REMAP_FILE} ]; then
cdo remap,${RESOL},${REMAP_FILE} diff_${EXP2}-${EXP1}_${TYPE}.nc diff_${EXP2}-${EXP1}_${TYPE}_${RESOL}.nc 1> /dev/null 2> /dev/null
fi
fi
if [ -f test$$.dat ]; then
  rm -f test$$.dat
fi
FILE1=${MODEL_DIR1}/experiments/${EXP1}/${EXP1}_${TYPE}.nc
FILE2=${MODEL_DIR2}/experiments/${EXP2}/${EXP2}_${TYPE}.nc
cdo diff ${FILE1} ${FILE2} > test$$.dat
STATUS=`echo $?`
if [ $STATUS -ne 0 -o -s test$$.dat ]; then
echo -e "\033[31m${TEST} test FAILED:\033[00m"
echo -e "The variables in files ${FILE1} and ${FILE2} differ or the files cannot be compared"
  DIFF_STATUS=$((DIFF_STATUS + 1))
fi
done # file types 
if [ -s test$$.dat ]; then
mv test$$.dat diff_${EXP2}-${EXP1}.dat
echo -e "\033[31mThe differences are in file diff_${EXP2}-${EXP1}.dat\033[00m"
else
rm -f test$$.dat
fi
if [ $DIFF_STATUS == 0 ]; then
echo -e "\033[32m${TEST} o.k.:\033[00m"
echo -e "The variables in files ${MODEL_DIR1}/experiments/${EXP1}/${EXP1}_TYPE.nc and ${MODEL_DIR2}/experiments/${EXP2}/${EXP2}_TYPE.nc for TYPE=${TYPES} do not differ"
fi
echo $DIFF_STATUS
#exit $DIFF_STATUS
}
#################### main script ########################################
SCRIPT_DIR=`pwd`
MODEL_DIR=${SCRIPT_DIR}/../
EXPERIMENT='r2b4_amip'
MODE='rn'
REFERENCE=''
OVERWRITE='yes'
while getopts ":e:hm:o:r:" OPT; do
case $OPT in
 h  ) print_usage $0
      exit 1 ;;
 m  ) MODE=$OPTARG 
      if [ $MODE != 'base' -a $MODE != 'update' -a $MODE != 'restart' -a $MODE != 'nproma' -a $MODE != 'ur' -a $MODE != 'un' -a $MODE != 'ur' -a $MODE != 'rn' -a $MODE != 'urn' ]; then
      print_usage $0
      exit 1
      fi 
      ;;
 e  ) EXPERIMENT=$OPTARG
      if [ $EXPERIMENT != 'r2b4_amip' -a $EXPERIMENT != 'mtest_nat_rce_cbl_120km_nwp' ]; then
      print_usage $0
      exit 1
      fi
      ;;
 r  ) REFERENCE=$OPTARG
      if [ $MODE == 'update' -o $MODE == 'ur' -o $MODE == 'un' -o $MODE == 'urn' ]; then
        if [ ! -d ${REFERENCE}/run ]; then
	   echo 'you asked for an update test, but no reference model found'
	   echo "reference model ''${REFERENCE}'' does not exist"
	   exit 1
	 fi
      fi
      ;;
 o  ) OVERWRITE=$OPTARG
      if [ $OVERWRITE != 'yes' -a $OVERWRITE != 'no' ]; then
      print_usage $0
      exit 1
      fi
      ;;
esac
done
if [ ${EXPERIMENT} == 'r2b4_amip' ]; then
TYPES='atm_phy atm_dyn lnd_phy'
DATES='19780101T004000Z 19780101T004500Z 19780101T005000Z 19780101T005500Z 19780101T010000Z'
RESOL='t63grid'
REMAP_FILE='/pool/data/ICON/grids/private/r2b4_amip/remapdis_r2b4_amip_cell_t63grid.nc'
RESTART_DATE='19780101T003000Z'
fi
if [ ${EXPERIMENT} == 'mtest_nat_rce_cbl_120km_nwp' ]; then
TYPES='atm'
DATES='20080801T004000Z 20080801T005000Z 20080801T010000Z'
RESOL='r320x20'
REMAP_FILE=''
RESTART_DATE='20080801T003000Z'
fi
echo 'test mode is '$MODE
if [ $MODE == 'update' -o $MODE == 'ur' -o $MODE == 'un' -o $MODE == 'urn' ]; then
  if [ -z "${REFERENCE}" ]; then
     echo -e '\033[31mno reference model was specified, use -r option to specify one\033[00m'
     print_usage
     exit 1
  else
     echo "reference model is ${REFERENCE}"
  fi
# get standardized paths
  cd ${REFERENCE}; REFERENCE=`pwd`
  cd ${MODEL_DIR}; MODEL_DIR=`pwd`
  if [ ${REFERENCE} == ${MODEL_DIR} ]; then
     echo 'reference model and test model are the same'
     exit 1
  fi
fi
SCRIPT=${EXPERIMENT}_test
EXP1=${EXPERIMENT}_base
EXP2=${EXPERIMENT}_restart
EXP3=${EXPERIMENT}_nproma
EXIT_STATUS=0
cd ${MODEL_DIR}/
STATUS=$?
if [ ${STATUS} == 0 ]; then
echo "found model ${MODEL_DIR}"
else
echo "no model ${MODEL_DIR} found"
exit
fi
rm -f ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
echo "Output will be in exp.${EXPERIMENT}_$$.output of"
echo "$SCRIPT_DIR"
${MODEL_DIR}/make_runscripts 1> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output 2> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
if [ ! -d experiments ]; then
mkdir experiments
fi
#################### perform base run without restart ####################
if [ ${OVERWRITE} == 'yes' -o ! -d ${MODEL_DIR}/experiments/${EXP1} ]; then
cd run
if [ ! -f exp.${SCRIPT}.run ]; then
echo 'did not find base runscript 'exp.${SCRIPT}.run
exit 1
fi
RUN_SCRIPT=exp.${EXP1}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP1}/g ${RUN_SCRIPT}
echo 'Performing base run'
${SCRIPT_DIR}/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output 2>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
STATUS=$?
if [ $STATUS -ne 0 ]; then
EXIT_STATUS=$((EXIT_STATUS + 1))
fi
else
echo 'Found base run'
fi # OVERWRITE
#################### perfrom update test              ####################
if [ ${MODE} == 'update' -o ${MODE} == 'ur' -o ${MODE} == 'un' -o ${MODE} == 'urn' ]; then
cd ${REFERENCE}
if [ ${OVERWRITE} == 'yes' -o ! -d ${REFERENCE}/experiments/${EXP1} ] ;then 
${REFERENCE}/make_runscripts 1> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output 2> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
if [ ! -d experiments ]; then
mkdir experiments
fi
cd run
if [ ! -f exp.${SCRIPT}.run ]; then
echo 'reference model: '$REFERENCE
echo 'did not find base runscript for reference model 'exp.${SCRIPT}.run 
exit 1
fi
RUN_SCRIPT=exp.${EXP1}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP1}/g ${RUN_SCRIPT}
echo 'Performing update test (run reference model)'
${REFERENCE}/run/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output 2>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
else
echo 'found run of reference model (update test)'
fi # OVERWRITE
diff_results $MODEL_DIR $EXP1 $REFERENCE $EXP1 "update"
if [ $DIFF_STATUS == 0 ]; then
EXIT_STATUS=$(($EXIT_STATUS + 0))
else
EXIT_STATUS=$(($EXIT_STATUS + 1))
fi
fi # MODE
#################### perform restart test             ####################
if [ ${MODE} == 'restart' -o ${MODE} == 'ur' -o ${MODE} == 'rn' -o ${MODE} == 'urn' ]; then
if [ ${OVERWRITE} == 'yes' -o ! -d ${MODEL_DIR}/experiments/${EXP2} ]; then
copy_experiment ${MODEL_DIR} ${EXP1} ${EXP2}
if [ $? -ne 0 ]; then
  echo "could not get base experiment ${EXP1} or create new experiment ${EXP2}"
  exit 1
fi
echo 'model_dir= ' ${MODEL_DIR}
cd ${MODEL_DIR}/run/
RUN_SCRIPT=exp.${EXP2}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/restart:=\".false.\"/restart:=\".true.\"/g ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP2}/g ${RUN_SCRIPT}
echo 'Performing restart test (running restart)'
${SCRIPT_DIR}/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output 2>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
else
echo 'found restart run (restart test)'
fi # OVERWRITE
# compare restart with base run
diff_results $MODEL_DIR $EXP1 $MODEL_DIR $EXP2 "restart"
if [ $DIFF_STATUS == 0 ]; then
EXIT_STATUS=$(($EXIT_STATUS + 0))
else
EXIT_STATUS=$(($EXIT_STATUS + 1))
fi
fi # MODE
#################### perform nproma test             ####################
if [ ${MODE} == 'nproma' -o ${MODE} == 'un' -o ${MODE} == 'rn' -o ${MODE} == 'urn' ]; then
if [ ${OVERWRITE} == 'yes' -o ! -d ${MODEL_DIR}/experiments/${EXP3} ]; then
copy_experiment ${MODEL_DIR} ${EXP1} ${EXP3}
cd ${MODEL_DIR}/run/
RUN_SCRIPT=exp.${EXP3}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/restart:=\".false.\"/restart:=\".true.\"/g ${RUN_SCRIPT}
sed -i s/nproma=64/nproma=17/g ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP3}/g ${RUN_SCRIPT}
echo 'Performing nproma test (running with nproma=17)'
${SCRIPT_DIR}/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output 2>> ${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
else
echo 'found run with nproma=17 (nproma test)'
fi # OVERWRITE
# compare nproma run with restarted run
diff_results $MODEL_DIR $EXP2 $MODEL_DIR $EXP3 "nproma"
if [ $DIFF_STATUS == 0 ]; then
EXIT_STATUS=$(($EXIT_STATUS + 0))
else
EXIT_STATUS=$(($EXIT_STATUS + 1))
fi
fi # MODE
########################################################################
if [ $EXIT_STATUS -eq 0 ]; then
  echo -e "\033[32mtest mode ${MODE}: ${EXPERIMENT} passed all corresponding tests\033[00m"
  echo "OK" > ${SCRIPT_DIR}/finish.status
else
  echo -e "\033[31mtest mode ${MODE}: ${EXPERIMENT} did NOT pass the corresponding tests\033[00m"
  echo -e "\033[31mnumber of tests including base run that FAILED: $EXIT_STATUS\033[00m"
fi
OFILE=${SCRIPT_DIR}/exp.${EXPERIMENT}_$$.output
if [ -s $OFILE ]; then
echo "Detailed output in exp.${EXPERIMENT}_$$.output"
fi
exit $EXIT_STATUS
