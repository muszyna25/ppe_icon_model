#!/usr/bin/env sh

set -x

#==============================================================================
# CONFIGURATION
#==============================================================================
# input
PHC_POOL_DIR='/pool/data/SEP/'
     PHC_DIR=${PHC_DIR:-${PHC_POOL_DIR}}
   PHC_FILES='PHC__3.0__SO__1x1__annual.nc  PHC__3.0__TempO__1x1__annual.nc'
#==============================================================================
# output
              MODEL=${MODEL:-icon}
         RESOLUTION=${RESOLUTION:-R2B04}
             LEVELS=${LEVELS:-L40}
             TARGET='ts_phc3.0_annual'
TARGET_MODEL_NOMISS=${TARGET}_${MODEL}_${RESOLUTION}-nomiss.nc
TARGET_MODEL_OUTPUT=${TARGET}_${MODEL}_${RESOLUTION}_${LEVELS}.nc
  TARGET_MODEL_SURF=${TARGET}_${MODEL}_${RESOLUTION}_surf.nc
#==============================================================================
# internals
                   CDO=${CDO:-cdo-dev}
               THREADS=${THREADS:-8}

 ICON_GRID_DIR_DEFAULT='/pool/data/ICON/ocean_data/ocean_grid'
         ICON_GRID_DIR=${ICON_GRID_DIR:-${ICON_GRID_DIR_DEFAULT}}

MPIOM_GRID_DIR_DEFAULT='/pool/data/MPIOM/input'
        MPIOM_GRID_DIR=${MPIOM_GRID_DIR:-${MPIOM_GRID_DIR_DEFAULT}}
# grids/weights for horizontal interpolation
declare -A remapFields
    remapFields[icon]=$ICON_GRID_DIR
   remapFields[mpiom]=$MPIOM_GRID_DIR
           targetGrid=${remapFields[${MODEL}]}/cell_grid-${RESOLUTION}.nc
         targetWeight=${remapFields[${MODEL}]}/cell_weight-${RESOLUTION}.nc

declare -A remapOperator
  remapOperator[icon]=gencon
 remapOperator[mpiom]=genbic
declare -A gridSelect
     gridSelect[icon]=ifs2icon_cell_grid
    gridSelect[mpiom]=area
#==============================================================================
# vertical levels
declare -A remapLevels
remapLevels[L20]='10,30,50,75,110,155,215,295,400,535,700,895,1125,1400,1750,2200,2750,3400,4100,4800'
remapLevels[L40]='6,17,27,37,47,57,68.5,82.5,100,122.5,150,182.5,220,262.5,310,362.5,420,485,560,645,740,845,960,1085,1220,1365,1525,1700,1885,2080,2290,2525,2785,3070,3395,3770,4195,4670,5170,5720'
#==============================================================================
# tmpfiles
   PHC_TMP='_phc3.0-annual.nc'
PHC_MERGED='phc3.0-annual.nc'
PHC_NOMISS='phc3.0-annual-nomiss.nc'
 TEMPFILES="${PHC_TMP} ${PHC_MERGED}"
#==============================================================================
 FORCE=${FORCE:-0} # re-create interpolation weights/grids
#==============================================================================


#==============================================================================
# MAIN SCRIPT
#==============================================================================
# get annual 1x1 phc3.0 data form pool/SEP
# convert to potential temperature (operator adipot, cdo1.6.2 required)
cd ${PHC_DIR}
$CDO merge ${PHC_FILES} $OLDPWD/${PHC_TMP}
cd -
$CDO -adipot -setcode,-1 -chname,SO,s,TempO,t ${PHC_TMP} ${PHC_MERGED}

#==============================================================================
# filling of land points
# TODO: fillmis produces unreal values in the baltic sea
$CDO -r -settaxis,2000-01-01,0,years  -fillmiss ${PHC_MERGED} ${PHC_NOMISS}

#==============================================================================
# horiz. interpolation for ICON or MPIOM target grid
# some preprocessing:
# * if the targeWeight are not present, they have to be created by with the targeGrid
# * if the targeGrid is not present, it has to be selection from the file
#   provides with the GRID variable
if [ ! -f ${targetWeight} -o "$FORCE" -eq 1 ];then 
  if [ ! -f ${targeGrid} -o "$FORCE" -eq 1 ]; then
    if [[ ! -f ${GRID} ]]; then
      echo "GRID variablen has to be set correctly!"
      exit 1
    else
      $CDO -selname,${gridSelect[${MODEL}]} ${GRID} ${targetGrid}
    fi
  fi
  $CDO -P ${THREADS} ${remapOperator[${MODEL}]},${targetGrid} ${PHC_NOMISS} ${targetWeight}
fi
$CDO -P ${THREADS} -remap,${targetGrid},${targetWeight} ${PHC_NOMISS} ${TARGET_MODEL_NOMISS}

#==============================================================================
# select relaxation top layer field
$CDO -sellevel,10 ${TARGET_MODEL_NOMISS} ${TARGET_MODEL_SURF}

#==============================================================================
# vertical interpolation 
$CDO -intlevelx,${remapLevels[${LEVELS}]} ${TARGET_MODEL_NOMISS} ${TARGET_MODEL_OUTPUT}

#==============================================================================
# clean up
for file in ${TEMPFILES}; do
  [[ -f $file ]] && rm $file
done
#==============================================================================
