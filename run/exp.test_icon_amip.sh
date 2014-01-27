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
if [ ! -f ${EXP_ORI}/restart.info -o ! -f ${EXP_ORI}/restart.r2b4_amip_19780101T003000Z_atm.nc ]; then
  echo "The experiment ${EXP_ORI} did not write any restart at 19780101T003000Z"
  exit 1
fi
rm -rf ${EXP_NEW}
mkdir ${EXP_NEW}
cd ${EXP_NEW}
cp ../${EXP_ORI}/restart.info .
cp ../${EXP_ORI}/restart.r2b4_amip_19780101T003000Z_atm.nc .
ln -sf restart.r2b4_amip_19780101T003000Z_atm.nc restart_atm.nc
}
#################### main script ####################
SCRIPT_DIR=`pwd`
MODEL_DIR=${SCRIPT_DIR}/../
SCRIPT=r2b4_amip_test
EXP1=r2b4_amip_base
EXP2=r2b4_amip_restart
EXP3=r2b4_amip_nproma
EXIT_STATUS=0
cd ${MODEL_DIR}/
STATUS=$?
if [ ${STATUS} == 0 ]; then
echo "found model ${MODEL_DIR}"
else
echo "no model ${MODEL_DIR} found"
exit
fi
${MODEL_DIR}/make_runscripts
if [ ! -d experiments ]; then
mkdir experiments
fi
cd experiments
# rm -rf ${EXP1} ${EXP2} ${EXP3}
cd ../run
if [ ! -f exp.${SCRIPT}.run ]; then
echo 'did not find base runscript 'exp.${SCRIPT}.run
exit
fi
# perform base run without restart 
RUN_SCRIPT=exp.${EXP1}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP1}/g ${RUN_SCRIPT}
${SCRIPT_DIR}/${RUN_SCRIPT}
# perform restart test
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
${SCRIPT_DIR}/${RUN_SCRIPT}
# compare restart with base run
$SCRIPT_DIR/diff_amip_test.sh $MODEL_DIR $EXP1 $EXP2
STATUS=$?
if [ $STATUS == 0 ]; then
EXIT_STATUS=$(($EXIT_STATUS + 0))
else
EXIT_STATUS=$(($EXIT_STATUS + 1))
fi
# perform nproma test
copy_experiment ${MODEL_DIR} ${EXP1} ${EXP3}
cd ${MODEL_DIR}/run/
RUN_SCRIPT=exp.${EXP3}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/restart:=\".false.\"/restart:=\".true.\"/g ${RUN_SCRIPT}
sed -i s/nproma=64/nproma=17/g ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP3}/g ${RUN_SCRIPT}
${SCRIPT_DIR}/${RUN_SCRIPT}
# compare nproma run with restarted run
$SCRIPT_DIR/diff_amip_test.sh $MODEL_DIR $EXP2 $EXP3
STATUS=$?
if [ $STATUS == 0 ]; then
EXIT_STATUS=$(($EXIT_STATUS + 0))
else
EXIT_STATUS=$(($EXIT_STATUS + 1))
fi
if [ $EXIT_STATUS -eq 0 ]; then
echo 'amip model passed rerun and nproma test'
else
echo 'amip model did not pass the rerun and/or nproma test'
fi
exit $EXIT_STATUS
