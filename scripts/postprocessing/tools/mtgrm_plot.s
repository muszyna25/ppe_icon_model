#!/bin/ksh
# -------------------------------------------------------
# Meteogram plots
# -------------------------------------------------------

set ex

dir="/e/uwork/mkoehler/icon/experiments/exp38/"
iFile=${dir}"NWP_iconR2B06_DOM01_2011010100_0001_meteogram.nc"
mkdir -p ${dir}"meteo"

set -A varName P  T  PEXNER  QV  QC  QI  QR  QS    \
   REL_HUM  RHO  THETAV  U  V  CLC  TKVM  TKVH  W  \
   Phalf  t_so_1  w_so_1  w_so_ice_1
set -A iStation 1 2 3

for station in ${iStation[*]}
do
for var in ${varName[*]}
do
   oFile=${dir}"meteo/NWP_iconR2B06_DOM01_2011010100_0001_meteogram.st"${station}"."${var}
   ncl -n mtgrm_plot.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" \
      varName=\"${var}\" iStation=${station}
done	
done	

exit
