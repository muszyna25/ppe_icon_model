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
OMIP_POOL_DIR='/pool/data/MPIOM/setup/omip_365_era15'
     OMIP_DIR=${OMIP_DIR:-${OMIP_POOL_DIR}}
   OMIP_FILES='*.nc'
#==============================================================================
# output
              MODEL=${MODEL:-icon}
               GRID=${GRID:-R2B04}
             TARGET=${TARGET:-omip}
        TARGET_FILE='omip'
                LEV=${LEV:-forcing}
         _test_file="${OMIP_DIR}/runoff.nc"
TARGET_MODEL_NOMISS=${TARGET}_${MODEL}_${GRID}-nomiss.nc
TARGET_MODEL_OUTPUT=${TARGET}_${MODEL}_${GRID}_${LEV}.nc
  TARGET_MODEL_SURF=${TARGET}_${MODEL}_${GRID}_surf.nc
#==============================================================================
# internals
                   CDO=${CDO:-cdo}
                CDODEV=${CDODEV:-cdo-dev}
               THREADS=${THREADS:-8}
                 FORCE=${FORCE:-0}
#==============================================================================
# basic input directories
 ICON_GRID_DIR_DEFAULT='/pool/data/ICON/ocean_data/ocean_grid'
         ICON_GRID_DIR=${ICON_GRID_DIR:-${ICON_GRID_DIR_DEFAULT}}

MPIOM_GRID_DIR_DEFAULT="/pool/data/MPIOM/${GRID}"
        MPIOM_GRID_DIR=${MPIOM_GRID_DIR:-${MPIOM_GRID_DIR_DEFAULT}}
# remapping setup
remapOperator_general=genbil
remapOperator_special=genbic
case "${MODEL}" in
 icon)
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
targetWeight_general=./cell_weight-${remapOperator_general}-${TARGET}-to-${MODEL}-${GRID}.nc
targetWeight_special=./cell_weight-${remapOperator_special}-${TARGET}-to-${MODEL}-${GRID}.nc
#==============================================================================
# MAIN
#==============================================================================
# horiz. interpolation for ICON or MPIOM target grid
# some preprocessing:
# * if the targeWeight are not present, they have to be created by with the targeGrid
# * if the targeGrid is not present, it has to be selection from the file
#   provides with the GRID variable
#
# create the target land-sea-mask
lsmFile='lsm4omip.nc'
$CDO gtc,1.0e-4 ${OMIP_DIR}/land_sea_mask.ECMWF.nc ${lsmFile}

if test ! -f "${targetWeight}" -o "$FORCE" -eq 1 ;then
  if test ! -f "${targetGrid}" ; then
    echo "GRID variable has to been set correctly!"
    exit 1
  fi
#  $CDO -P ${THREADS} ${remapOperator_general},${targetGrid} -fillmiss -ifthen ${lsmFile} ${_test_file} ${targetWeight_general}
#  $CDO -P ${THREADS} ${remapOperator_special},${targetGrid} -fillmiss -ifthen ${lsmFile} ${_test_file} ${targetWeight_special}
fi
#==============================================================================
TMP='_omip.nc'
MERGED='omip.nc'
NOMISS='phc3.0-annual-nomiss.nc'
TEMPFILES="${TMP}" # ${MERGED}"
#==============================================================================
FORCE=${FORCE:-0} # re-create interpolation weights/grids

[[ -f jobs ]] && rm jobs
[[ -f $(ls -1 remapped_*.nc) ]] && rm remapped_*nc
[[ -f $(ls -1 fillmiss_*.nc) ]] && rm fillmiss_*nc
for file in $(ls ${OMIP_DIR}/${OMIP_FILES}); do 
  oFile=remapped_$(basename $file)
  _oFile=fillmiss_$(basename $file)
  echo -n "$CDO fillmiss -ifthen ${lsmFile}  $file ${_oFile};" >> jobs
  echo -n "$CDO remap,${targetGrid},${targetWeight_general} ${_oFile} ${oFile};" >> jobs
  echo "" >> jobs
done
cat jobs | parallel -j ${THREADS}
#==============================================================================
# clean up
for file in ${TEMPFILES}; do
  [[ -f $file ]] && rm $file
done
#==============================================================================
