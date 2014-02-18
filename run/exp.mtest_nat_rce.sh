#!/bin/ksh
#################### functions ######################
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
if [ ! -f ${EXP_ORI}/restart.info -o ! -f ${EXP_ORI}/torus_grid_r120000s64_restart_atm_20080801T003000Z.nc ]; then
  echo -e "\033[31mThe experiment ${EXP_ORI} did not write any restart at 20080801T003000Z\033[00m"
  exit 1
fi
rm -rf ${EXP_NEW}
mkdir ${EXP_NEW}
cd ${EXP_NEW}
cp ../${EXP_ORI}/restart.info .
cp ../${EXP_ORI}/torus_grid_r120000s64_restart_atm_20080801T003000Z.nc .
ln -sf torus_grid_r120000s64_restart_atm_20080801T003000Z.nc restart_atm_DOM01.nc
}
#################### main script ####################
SCRIPT_DIR=`pwd`
MODEL_DIR=${SCRIPT_DIR}/../
MODE='rn'
REFERENCE=''
OVERWRITE='yes'
while getopts ":hm:r:o:" OPT; do
case $OPT in
 h  ) echo "usage: $0 [-h] [-m base|update|restart|nproma|ur|un|rn|urn ] [-r <reference model path>] [-o yes|no ]"
      echo '-h : display help'
      echo '-m : test mode, either single tests like update, restart or nproma test or'
      echo '     combined tests like ur (update and restart) possible'
      echo '-r : reference model path'
      echo '-o : ''yes'': overwrite existing experiments, that is the default'
      echo '     ''no'': use existing experiments'
      exit 1 ;;
 m  ) MODE=$OPTARG 
      if [ $MODE != base -a $MODE != update -a $MODE != restart -a $MODE != nproma -a $MODE != ur -a $MODE != un -a $MODE != ur -a $MODE != rn -a $MODE != urn ]; then
      echo "usage: $0 [-h] [-m update|restart|nproma|ur|un|rn|urn ] [-r <reference model path>]" 
      exit 1
      fi 
      ;;
 r  ) REFERENCE=$OPTARG
      echo "REFERENCE=$REFERENCE"
      if [ $MODE == 'update' -o $MODE == 'ur' -o $MODE == 'un' -o $MODE == 'urn' ]; then
        if [ ! -d ${OPTARG}/run ]; then
	   echo 'you asked for an update test, but no reference model found'
	   echo "reference model ''${REFERENCE}'' does not exist"
	   exit 1
	 fi
      fi
      ;;
 o  ) OVERWRITE=$OPTARG ;;
esac
done
echo 'test mode is '$MODE
if [ $MODE == 'update' -o $MODE == 'ur' -o $MODE == 'un' -o $MODE == 'urn' ]; then
  if [ ! -n ${REFERENCE} ]; then
     echo 'no reference model was specified, use -r option to specify one'
     exit 1
  else
     echo "reference model is ${REFERENCE}"
  fi
  cd ${MODEL_DIR}
  MODEL_DIR=`pwd`
  if [ ${REFERENCE} == ${MODEL_DIR} ]; then
     echo 'reference model and test model are the same'
     exit 1
  fi
fi
SCRIPT=mtest_nat_rce_cbl_120km_nwp_test
EXP1=mtest_nat_rce_cbl_120km_nwp_test_base
EXP2=mtest_nat_rce_cbl_120km_nwp_test_restart
EXP3=mtest_nat_rce_cbl_120km_nwp_test_nproma
EXIT_STATUS=0
cd ${MODEL_DIR}/
STATUS=$?
if [ ${STATUS} == 0 ]; then
echo "found model ${MODEL_DIR}"
else
echo "no model ${MODEL_DIR} found"
exit
fi
rm -f ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
echo "Output will be in exp.mtest_nat_rce$$.output of"
echo "$SCRIPT_DIR"
${MODEL_DIR}/make_runscripts 1> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output 2> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
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
${SCRIPT_DIR}/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output 2>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
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
${REFERENCE}/make_runscripts 1> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output 2> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
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
${REFERENCE}/run/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output 2>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
else
echo 'found run of reference model (update test)'
fi # OVERWRITE
${SCRIPT_DIR}/diff_nat_rce_test.sh $MODEL_DIR $EXP1 $REFERENCE $EXP1 "update"
STATUS=$?
if [ $STATUS == 0 ]; then
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
${SCRIPT_DIR}/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output 2>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
else
echo 'found restart run (restart test)'
fi # OVERWRITE
# compare restart with base run
$SCRIPT_DIR/diff_nat_rce_test.sh $MODEL_DIR $EXP1 $MODEL_DIR $EXP2 "restart"
STATUS=$?
if [ $STATUS == 0 ]; then
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
${SCRIPT_DIR}/${RUN_SCRIPT} 1>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output 2>> ${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
else
echo 'found run with nproma=17 (nproma test)'
fi # OVERWRITE
# compare nproma run with restarted run
$SCRIPT_DIR/diff_nat_rce_test.sh $MODEL_DIR $EXP2 $MODEL_DIR $EXP3 "nproma"
STATUS=$?
if [ $STATUS == 0 ]; then
EXIT_STATUS=$(($EXIT_STATUS + 0))
else
EXIT_STATUS=$(($EXIT_STATUS + 1))
fi
fi # MODE
########################################################################
if [ $EXIT_STATUS -eq 0 ]; then
echo -e "\033[32mtest mode ${MODE}: amip model passed all corresponding tests\033[00m"
else
echo -e "\033[31mtest mode ${MODE}: amip model did NOT pass the corresponding tests\033[00m"
echo -e "\033[31mnumber of tests including base run that FAILED: $EXIT_STATUS\033[00m"
fi
OFILE=${SCRIPT_DIR}/exp.mtest_nat_rce$$.output
if [ -s $OFILE ]; then
echo "Detailed output in exp.mtest_nat_rce$$.output"
fi
exit $EXIT_STATUS
