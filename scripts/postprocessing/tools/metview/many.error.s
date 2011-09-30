#!/bin/ksh
#
# run as: many.error.s 20110101 00 20110102 00 25
# -------------------------------------------------------
set -x
	
inidate=${1}
initime=${2}
verdate=${3}
vertime=${4}
expnum=${5}
echo "Arguments: inidate="${inidate}" initime="${initime}" verdate="${verdate}" vertime="${vertime}" expnum="${expnum}

#expnum=25

#              00h      24h      48h      120h
#set -A inidate 20110101 20110101 20110101 20110101
#set -A initime 00       00       00       00      
#set -A verdate 20110101 20110102 20110103 20110106
#set -A vertime 00       00       00       00

#set -A inidate 20110101 
#set -A initime 00       
#set -A verdate 20110102 
#set -A vertime 00       

mkdir -p "/fe1-daten/"$USER"/plots/icon/nwp.exp"${expnum}
scriptdir="/fe1-daten/"$USER"/metview/ICON/"

integer nt
while [[ $nt < ${#inidate[*]} ]]; do

  # -------------------------------------------------------
  
  set -A vars TQV  TQC  TQI   TCC  TQ1  TQ2  TQ3                                 \
              swflxsfc_acc    lwflxsfc_acc   swflxtoa_acc     lwflxtoa_acc       \
              LHFL_S_acc      SHFL_S_acc                                         \
              TOT_PREC        RAIN_GSP       SNOW_GSP         RAIN_CON  SNOW_CON \
              Z0  T_G  QV_S   T_2m  QV_2m    U_10m  V_10m                        \
              T_GT_tile_1     T_S_tile_1     W_I_tile_1                          \
              T_SNOW_tile_1   DZH_SNOW_tile_1                 H_SNOW_tile_1      \
              W_SNOW_tile_1   WTOT_SNOW_tile_1                WLIQ_SNOW_tile_1   \
              RHO_SNOW_tile_1 RHO_SNOW_MULT_tile_1
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}map.error $expnum $var sfc snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done

  set -A vars TQV             TQC            TQI              TCC                \
              swflxsfc_acc    lwflxsfc_acc   swflxtoa_acc     lwflxtoa_acc       \
              LHFL_S_acc      SHFL_S_acc     TOT_PREC                            \
              T_G             T_2m           U_10m            V_10m              \
              H_SNOW_tile_1   RHO_SNOW_tile_1
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}map.error $expnum $var sfc ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error $expnum $var sfc diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done
 
  # -------------------------------------------------------
 
  set -A vars T_SO_tile_1  W_SO_tile_1  W_SO_ICE_tile_1  T_SNOW_MULT_tile_1
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}map.error $expnum $var lnd snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error $expnum $var lnd ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error $expnum $var lnd diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done
  
  # -------------------------------------------------------
  	
  set -A vars T U V QV QC QI CC Q1 Q2 Q3
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}zonal.error $expnum $var ml snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}zonal.error $expnum $var ml ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}zonal.error $expnum $var ml diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var ml snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var ml ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var ml diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done
    
  # -------------------------------------------------------

  set -A vars T_P  QV_P  U_P  V_P  Z
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}zonal.error $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}zonal.error $expnum $var pl ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}zonal.error $expnum $var pl diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var pl ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var pl diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done

  set -A vars QC_P  QI_P  CC_P  Q1_P  Q2_P  Q3_P  Q4_P  Q5_P
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}zonal.error $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var pl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done

  # -------------------------------------------------------

  set -A vars T_Z  U_Z  V_Z  P_Z  QV_Z  QC_Z  QI_Z  CC_Z  Q1_Z  Q2_Z  Q3_Z  Q4_Z  Q5_Z
  for var in ${vars[*]}
  do
    metview -b ${scriptdir}map.error   $expnum $var zl snap  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var zl ctr   ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
    metview -b ${scriptdir}map.error   $expnum $var zl diff  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]}
  done


  nt=nt+1
done  # end loop date

exit
