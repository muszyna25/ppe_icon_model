#!/bin/ksh

# --- setup -----------------------------------------------------------

turb='edmf'   # 'edmf' or 'turbdiff'
code='14217'
case='rico'   # 'cbl' or 'rico'

cdo=/fe1-daten/mkoehler/bin/cdo-1.6.0/src/cdo

expdir="/fe1-daten/mkoehler/icon-dev/experiments/nh_"${case}"_"${turb}
plotdir=${expdir}"/plots_rev"${code}
plotname="torus"
mkdir -p ${plotdir}

var="temp"
oType="eps"     # "png" or "eps"
title="RICO"

# --- horizontal mean -------------------------------------------------

file_tri=${expdir}"/nh_"${case}"_"${turb}"_tri_DOM01_ML_0001.nc"
file_avg=${expdir}"/nh_"${case}"_"${turb}"_tri_DOM01_ML_0001_avg.nc"
${cdo} fldmean ${file_tri} ${file_avg}

iFile=${file_avg}


# --- surface evolution plots -----------------------------------------

# for var in 'u_10m'
for var in 'clct' 'tqv' 'tqc' 'tqi' 'tot_prec'            \
  'rain_gsp' 'rain_con' 'snow_gsp' 'snow_con'             \
  'ashfl_s' 'alhfl_s' 'athb_s' 'athb_t' 'asob_s' 'asob_t' \
  'u_10m' 'v_10m' 't_2m' 't_g'
do
  echo "---------------------------------- plotting:"  $var " -----------------"
  oFile=${plotdir}/${var}_sfc_${case}"_"${turb}
  ncl -n plot_torus_sfc.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
     varName=\"${var}\" expnum=\"${title}\"
done


# --- vertical profile plots ------------------------------------------

# w has different levels can not plot it here
#for var in 'u'
for var in 'u' 'v' 'temp' 'theta_v' 'qv' 'qc' 'qi' 'qr' 'qs' 'clc'    \
  'tot_qv' 'tot_qc' 'tot_qi' 'qtvar'            \
  'ddt_temp_radsw' 'ddt_temp_radlw'  'ddt_temp_turb' 'ddt_temp_pconv' \
  'ddt_u_turb' 'ddt_u_pconv' 'ddt_v_turb' 'ddt_v_pconv' 'ddt_qv_turb' 'ddt_qv_conv' 
 # 'z_ifc' 'z_mc' 'pres' 
do
  echo "---------------------------------- plotting:"  $var " -----------------"
  oFile=${plotdir}/${var}_prof_${case}"_"${turb}
  ncl -n plot_torus_prof.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
    varName=\"${var}\" expnum=\"${title}\"
done


# --- vertical profile plots ------------------------------------------

#for var in 'qc'
for var in 'u' 'v' 'temp' 'theta_v' 'qv' 'qc' 'qi' 'qr' 'qs' 'clc'    \
  'tot_qv' 'tot_qc' 'tot_qi' 'qtvar'            \
  'ddt_temp_radsw' 'ddt_temp_radlw'  'ddt_temp_turb' 'ddt_temp_pconv' \
  'ddt_u_turb' 'ddt_u_pconv' 'ddt_v_turb' 'ddt_v_pconv' 'ddt_qv_turb' 'ddt_qv_conv' 
 # 'z_ifc' 'z_mc' 'pres' 
do
  echo "---------------------------------- plotting:"  $var " -----------------"
  oFile=${plotdir}/${var}_tz_${case}"_"${turb}
  ncl -n plot_torus_tz.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
    varName=\"${var}\" expnum=\"${title}\"
done


