#!/usr/bin/env sh

set -x

#==============================================================================
# Usage:
#   MODEL=icon GRID=R2B04 PHC_DIR=. FORCE=1 LEV=L20 ./preproc4phc.sh
#   MODEL=mpiom GRID=GR15 PHC_DIR=. FORCE=1 LEV=L40 ./preproc4phc.sh
#==============================================================================
# CONFIGURATION
#==============================================================================
# input
PHC_POOL_DIR='/pool/data/SEP/data_sources/PHC3.0/DATA'
     PHC_DIR=${PHC_DIR:-${PHC_POOL_DIR}}
   PHC_FILES='PHC__3.0__SO__1x1__annual.nc  PHC__3.0__TempO__1x1__annual.nc'
#==============================================================================
# output
              MODEL=${MODEL:-icon}
               GRID=${GRID:-R2B04}
                LEV=${LEV:-L40}
             TARGET='ts_phc3.0_annual'
TARGET_MODEL_NOMISS=${TARGET}_${MODEL}_${GRID}-nomiss.nc
TARGET_MODEL_OUTPUT=${TARGET}_${MODEL}_${GRID}_${LEV}.nc
  TARGET_MODEL_SURF=${TARGET}_${MODEL}_${GRID}_surf.nc
#==============================================================================
# internals
                   CDO=${CDO:-cdo-dev}
               THREADS=${THREADS:-8}
#==============================================================================
# basic input directories
 ICON_GRID_DIR_DEFAULT='/pool/data/ICON/ocean_data/ocean_grid'
         ICON_GRID_DIR=${ICON_GRID_DIR:-${ICON_GRID_DIR_DEFAULT}}

MPIOM_GRID_DIR_DEFAULT="/pool/data/MPIOM/${GRID}"
        MPIOM_GRID_DIR=${MPIOM_GRID_DIR:-${MPIOM_GRID_DIR_DEFAULT}}
case "${MODEL}" in
 icon)
   remapOperator=genbil
   gridSelect=ifs2icon_cell_grid
   GRID_FILE=$(ls ${ICON_GRID_DIR}/icon${GRID}*etop*planet.nc | head -1)
  targetGrid=./cell_grid-${GRID}-${MODEL}.nc
  ${CDO} -f nc -selname,${gridSelect} ${GRID_FILE} ${targetGrid}
   ;;
 mpiom)
   remapOperator=genbil
   gridSelect=''
   GRID_FILE="${MPIOM_GRID_DIR}/${GRID}s.nc"
   targetGrid=./cell_grid-${GRID}-${MODEL}.nc
   ${CDO} -f nc -sethalo,1,1 -random,${GRID_FILE} ${targetGrid}
   ;;
 *) echo "Unsupported model! Use 'icon' or 'mpiom'."; exit 1;;
esac
# grids/weights for horizontal interpolation
targetWeight=./cell_weight-${remapOperator}-phc3-to-${MODEL}-${GRID}.nc
#==============================================================================
# vertical levels
case "${LEV}" in
  L10) remapLevels='20,65,135,260,475,805,1275,2000,3100,4900';;
  L20) remapLevels='10,30,50,75,110,155,215,295,400,535,700,895,1125,1400,1750,2200,2750,3400,4100,4800';;
  L40) remapLevels='6,17,27,37,47,57,68.5,82.5,100,122.5,150,182.5,220,262.5,310,362.5,420,485,560,645,740,845,960,1085,1220,1365,1525,1700,1885,2080,2290,2525,2785,3070,3395,3770,4195,4670,5170,5720';;
  L80) remapLevels='6,17,27,37,47,57,67,77,88,100,112,125,139,153,167,183,199,215,233,251,271,291,312,334,357,381,407,433,461,489,519,551,584,618,654,692,732,773,816,861,908,957,1008,1062,1119,1178,1239,1304,1371,1442,1516,1593,1674,1759,1848,1941,2038,2140,2246,2356,2472,2594,2721,2854,2993,3138,3290,3449,3616,3790,3972,4162,4362,4570,4788,5016,5255,5504,5765,6038';;
  *) echo "Unsupported Levels - use L10,L20,L40 or L80";exit 1;;
esac
#==============================================================================
# tmpfiles
   PHC_TMP='_phc3.0-annual.nc'
PHC_MERGED='phc3.0-annual.nc'
PHC_NOMISS='phc3.0-annual-nomiss.nc'
 TEMPFILES="${PHC_TMP}" # ${PHC_MERGED}"
#==============================================================================
 FORCE=${FORCE:-0} # re-create interpolation weights/grids
#==============================================================================


#==============================================================================
# MAIN SCRIPT
#==============================================================================
# get annual 1x1 phc3.0 data form pool/SEP
# convert to potential temperature (operator adipot, cdo1.6.2 required)
cd ${PHC_DIR}
$CDO -O merge ${PHC_FILES} $OLDPWD/${PHC_TMP}
cd -
$CDO -adipot -setcode,-1 -chname,SO,s,TempO,t ${PHC_TMP} ${PHC_MERGED}

#==============================================================================
# filling of land points
# TODO: fillmis produces unreal values in the baltic sea
$CDO -P ${THREADS} -O -r -settaxis,2000-01-01,0,years -fillmiss -ifthenelse -setmisstoc,0 ${PHC_MERGED} ${PHC_MERGED}  -remapbil,${PHC_MERGED} -fillmiss -fillmiss -sellonlatbox,20,28,64,66 ${PHC_MERGED} ${PHC_NOMISS}
#$CDO -O -r -settaxis,2000-01-01,0,years  -fillmiss ${PHC_MERGED} ${PHC_NOMISS}

#==============================================================================
# horiz. interpolation for ICON or MPIOM target grid
# some preprocessing:
# * if the targeWeight are not present, they have to be created by with the targeGrid
# * if the targeGrid is not present, it has to be selection from the file
#   provides with the GRID variable
if test ! -f "${targetWeight}" -o "$FORCE" -eq 1 ;then
  if test ! -f "${targetGrid}" ; then
    echo "GRID variable has to been set correctly!"
    exit 1
  fi
  $CDO -P ${THREADS} ${remapOperator},${targetGrid} ${PHC_NOMISS} ${targetWeight}
fi
$CDO -remap,${targetGrid},${targetWeight} ${PHC_NOMISS} ${TARGET_MODEL_NOMISS}

#==============================================================================
# select relaxation top layer field
$CDO -sellevel,10 ${TARGET_MODEL_NOMISS} ${TARGET_MODEL_SURF}

#==============================================================================
# vertical interpolation
$CDO -intlevelx,${remapLevels} ${TARGET_MODEL_NOMISS} ${TARGET_MODEL_OUTPUT}

#==============================================================================
# clean up
for file in ${TEMPFILES}; do
  [[ -f $file ]] && rm $file
done
#==============================================================================
