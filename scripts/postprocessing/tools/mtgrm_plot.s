#!/bin/ksh
# -------------------------------------------------------
# Meteogram plots
#
# run as: mtgrm_plot.s 2011010100 exp24
#
# attention:
#   scripts need to be local (mtgrm_plot.s, mtgrm_plot.ncl)
#
# Martin Koehler, Dec 2011
# -------------------------------------------------------

set ex

dates=${1}
expnum=${2}
echo "mtgrm_plot.s --- Arguments: dates="${dates}" expnum="${expnum}


dir="/e/uwork/mkoehler/icon/experiments/"${expnum}"/"
iFile=${dir}"NWP_iconR2B06_DOM01_"${dates}"_0001_meteogram.nc"
mkdir -p ${dir}"meteo"
oType="png"

set -A iStation 1 2 3 4 5 6 7 8 9

set -A varNameSfc \
  P_SFC    PL_Cov   LA_Ind  RO_Dept    Z0         qv_s       w_i_1    w_snow_1 \
  TCM      TCH      SHFL    LHFL       RUNOFF_S_1 RUNOFF_G_1 VIO3     HMO3     \
  t_snow_1 t_s_1    t_g     FRESHSNW_1 RHO_SNOW_1 H_SNOW_1   T2M      TD2M     \
  U10M     V10M     SOBT    SOBS       THBS       ALB        RAIN_GSP SNOW_GSP \
  RAIN_CON SNOW_CON 

set -A varName3D \
  P        T        PEXNER  QV         QC  QI   QR    QS       \
  REL_HUM  RHO      THETAV  U          V   CLC  TKVM  TKVH  W  \
  Phalf    t_so_1   w_so_1  w_so_ice_1

# -------------------------------------------------------

for station in ${iStation[*]}
do
for var in ${varNameSfc[*]}
do
   oFile=${dir}"meteo/NWP_iconR2B06_DOM01_"${dates}"_0001_meteogram.loc"${station}"."${var}
   ncl -n mtgrm_plot_sfc.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
      varName=\"${var}\" iStation=${station}
done	
done	

# -------------------------------------------------------

for station in ${iStation[*]}
do
for var in ${varName3D[*]}
do
   oFile=${dir}"meteo/NWP_iconR2B06_DOM01_"${dates}"_0001_meteogram.loc"${station}"."${var}
   ncl -n mtgrm_plot.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
      varName=\"${var}\" iStation=${station}
done	
done	

exit
