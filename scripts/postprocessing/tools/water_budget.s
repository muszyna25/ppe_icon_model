#!/bin/ksh
# -------------------------------------------------------
# Calculate water budget from output
#
# run as: ./water_budget.s /scratch/ms/de/dei2/temp/dei2_146_2012010100_DOM01_ML_0001.grb
#
# tip: split GRIB file into typeOfLevel's
#      grib_copy dei2_146_2012010100_DOM01_ML_0001.grb dei2_146_2012010100_DOM01_ML_0001.[typeOfLevel].grb
#      ./water_budget.s /scratch/ms/de/dei2/temp/dei2_146_2012010100_DOM01_ML_0001.surface.grb
#
# Martin Koehler, May 2014
# -------------------------------------------------------

ml_file=${1}

module load cdo
LANG=''          # for math with 1.0 not 1,0

set -A date     `cdo infov ${ml_file} | grep slhf | awk '{print $3}'`
set -A time     `cdo infov ${ml_file} | grep slhf | awk '{print $4}'`
set -A lhf      `cdo output -fldmean -selname,slhf        ${ml_file}`
set -A tot_prec `cdo output -fldmean -selname,tp          ${ml_file}`   
set -A tqv      `cdo output -fldmean -selname,tciwv       ${ml_file}`   
set -A tql      `cdo output -fldmean -selname,param69.1.0 ${ml_file}`
set -A tqi      `cdo output -fldmean -selname,param70.1.0 ${ml_file}`
set -A tqr      `cdo output -fldmean -selname,tcolr       ${ml_file}`
set -A tqs      `cdo output -fldmean -selname,tcols       ${ml_file}`  

float lhf1
float pme
float dqt
float budget

echo '________________________________________________________________________________________________________________________________________________________________'
echo 'date         time        lhf            tot_prec     tqv           tql           tqi           tqr           tqs           P-E           del qt         budget'
echo '                         [mm]           [mm]         [mm]          [mm]          [mm]          [mm]          [mm]          [mm/12h]      [mm/12h]       [mm/12h]'
echo '________________________________________________________________________________________________________________________________________________________________'
for ((nn=0; nn<=20; nn+=1)) ; do

  lhf1=${lhf[$nn]}/2500800.0
  if [[ ${nn} -eq 0 ]] ; then
    pme=0.0
    dqt=0.0
    budget=0.0
  else
    pme=$(( ${tot_prec[$nn]}+${lhf[$nn]}/2500800.0 - ${tot_prec[$nn-1]} - ${lhf[$nn-1]}/2500800.0 ))
    dqt=$(( ${tqv[$nn]}+${tql[$nn]}+${tqi[$nn]}+${tqr[$nn]}+${tqs[$nn]} - ${tqv[$nn-1]} - ${tql[$nn-1]} - ${tqi[$nn-1]} - ${tqr[$nn-1]} - ${tqs[$nn-1]} ))
    budget=$(( ${pme} + ${dqt} ))
  fi

  printf '%-12s %-11s %-14s %-12s %-13s %-13s %-13s %-13s %-13s %-17s %-15s %-15s \n' \
    ${date[$nn]} ${time[$nn]} ${lhf1}     ${tot_prec[$nn]}                                  \
    ${tqv[$nn]}  ${tql[$nn]}  ${tqi[$nn]} ${tqs[$nn]}      ${tqr[$nn]} ${pme} ${dqt} ${budget}
done
echo '_________________________________________________________________________________________________________________________________________________________________'


exit

