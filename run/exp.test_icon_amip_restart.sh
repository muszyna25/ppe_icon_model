#!/bin/ksh
SCRIPT_DIR=`pwd`
MODEL_DIR=${SCRIPT_DIR}/../
SCRIPT=r2b4_amip_restart_test
EXP1=r2b4_amip_without_restart
EXP2=r2b4_amip_restart
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
rm -rf ${EXP1} ${EXP2}
cd ../run
# run without restart first
RUN_SCRIPT=exp.${EXP1}.run
if [ ! -f exp.${SCRIPT}.run ]; then
echo 'did not find runscript 'exp.${SCRIPT}.run
exit
fi
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP1}/g ${RUN_SCRIPT}
${SCRIPT_DIR}/${RUN_SCRIPT}
cd ../experiments
mkdir ${EXP2}
cd ${EXP2}
cp ../${EXP1}/restart.info .
cp ../${EXP1}/restart.r2b4_amip_19780101T003000Z_atm.nc .
ln -s restart.r2b4_amip_19780101T003000Z_atm.nc restart_atm.nc
cd ${MODEL_DIR}/run
RUN_SCRIPT=exp.${EXP2}.run
cp -f exp.${SCRIPT}.run ${RUN_SCRIPT}
sed -i s/restart:=\".false.\"/restart:=\".true.\"/g ${RUN_SCRIPT}
sed -i s/${SCRIPT}/${EXP2}/g ${RUN_SCRIPT}
${SCRIPT_DIR}/${RUN_SCRIPT}
$SCRIPT_DIR/diff_restart.sh $MODEL_DIR $EXP1 $EXP2
STATUS=$?
if [ $STATUS == 0 ]; then
exit 0
else
exit 1
fi
