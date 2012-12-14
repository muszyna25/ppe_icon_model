#!/bin/ksh
#
# run as: many.error.s 20110101 00 20110102 00 1 R2B04 25
# -------------------------------------------------------
set -x
	
inidate=${1}
initime=${2}
verdate=${3}
vertime=${4}
ndays=${5}
res=${6}
expnum=${7}
echo "Arguments: inidate="${inidate}" initime="${initime}" verdate="${verdate}" vertime="${vertime}" ndays="${ndays}" res="${res}" expnum="${expnum}

#expnum=34

#              00h      24h      48h      120h
#set -A inidate 20110101 20110101 20110101 20110101
#set -A initime 00       00       00       00      
#set -A verdate 20110101 20110102 20110103 20110106
#set -A vertime 00       00       00       00

#set -A inidate 20110101 
#set -A initime 00       
#set -A verdate 20110102 
#set -A vertime 00       

#mkdir -p "/fe1-daten/"$USER"/plots/icon/nwp.exp"${expnum}
#scriptdir="/fe1-daten/"$USER"/metview/ICON/"
scriptdir="./"
cd ${scriptdir}

#metview=metview4_new
metview=metview4_dev
#metview=/usr/local/apps/Metview/metview4_expt

met_job=met.job.all.$nstart
\rm -rf $met_job

integer nt
while [[ $nt < ${#inidate[*]} ]]; do

  # -------------------------------------------------------
  
  set -A vars TQV  TQC  TQI   TCC  TQ1  TQ2  TQ3  PS                             \
              ACCSOB_S        ACCTHB_S       ACCSOB_T         ACCTHB_T           \
              ACCLHFL_S       ACCSHFL_S                                          \
              TOT_PREC        RAIN_GSP       SNOW_GSP         RAIN_CON  SNOW_CON \
              Z0  T_G  QV_S   T_2M  QV_2M    U_10M  V_10M                        \
              T_GT_tile_1     T_S_tile_1     W_I_tile_1                          \
              T_SNOW_tile_1   DZH_SNOW_tile_1                 H_SNOW_tile_1      \
              W_SNOW_tile_1   WTOT_SNOW_tile_1                WLIQ_SNOW_tile_1   \
              RHO_SNOW_tile_1 RHO_SNOW_MULT_tile_1
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $var sfc snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
  set -A vars TQV             TQC            TQI              TCC                \
              ACCSOB_S        ACCTHB_S       ACCSOB_T         ACCTHB_T           \
              ACCLHFL_S       ACCSHFL_S      TOT_PREC         PS                 \
              T_G             T_2M           U_10M            V_10M              \
              H_SNOW_tile_1   RHO_SNOW_tile_1 
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $var sfc diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $var sfc ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done

  # -------------------------------------------------------
 
#  set -A vars T_SO_tile_1  W_SO_tile_1  W_SO_ICE_tile_1  T_SNOW_MULT_tile_1
#  for var in ${vars[*]}
#  do
#    echo${metview} -b ${scriptdir}map.error $expnum $var lnd snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo${metview} -b ${scriptdir}map.error $expnum $var lnd ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo${metview} -b ${scriptdir}map.error $expnum $var lnd diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#  done
  
  # -------------------------------------------------------
  	
  set -A vars T  U  V  Q1  Q2  Q3  QV  QC  QI  CC  P  QR  QS #QR QS  QTVAR  O3  P               
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var ml diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var ml snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var ml ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $var ml diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $var ml snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $var ml ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
    
  # -------------------------------------------------------
  	
  set -A vars ttendcds qtendcds utendcds vtendcds \
              ttendts  qtendt   utendts  vtendts  \
              ttends   utends   vtends            \
              ewgd     nsgd                       \
              ttendsw  ttendlw  O3
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var ml snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $var ml snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
    
  # -------------------------------------------------------

  set -A vars T  QV  U  V  FI
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var pl diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var pl ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $var pl diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $var pl ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done

  set -A vars  Q1  Q2  Q3  QC  QI  CC # QR QS
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done

  # -------------------------------------------------------

  set -A vars T  U  V  P  Q1  Q2  Q3  QV  QC  QI  CC  # CRWC CSWC
#  for var in ${vars[*]}
#  do
#    echo ${metview} -b ${scriptdir}map.error   $expnum $var zl diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo ${metview} -b ${scriptdir}map.error   $expnum $var zl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo ${metview} -b ${scriptdir}map.error   $expnum $var zl ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#  done


  nt=nt+1
done  # end loop date

exit
