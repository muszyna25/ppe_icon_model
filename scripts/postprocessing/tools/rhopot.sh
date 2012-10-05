# plot T,S and potential density variation to initial values from a list of ICON input files
#
# ==============================================================================
   DEBUG='TRUE'
  CDOOPT='-O'
     CDO="cdo-dev $CDOOPT"  # my developement version of cdo with a workong rhopot operator
ICONPLOT=$HOME/src/git/icon/scripts/postprocessing/tools/icon_plot.ncl

# ==============================================================================
# littel helper function fo debuggin
function call {
  if [ ! -z "$DEBUG" ];then
    echo "CALLING:'$@'"
  fi
  $@
}
# ==============================================================================
# ==============================================================================
# we supppose, that all files belong to the same experiment
  fileListPattern='RAM_test_oce_withIce2_*.nc'
     Temp_VarName='T'
 Salinity_VarName='S'
PotDensityVarName='rhopot'
      MaskVarName='wet_c'
   outputDataFile='output.nc'

# ==============================================================================
declare -a fileListArray
fileList=$(ls $fileListPattern)
i=0
for file in $fileList; do
  fileListArray[$i]=$file
  i=$((i+1))
done
numberOfFiles=$i
echo "$numberOfFiles number of files"
echo $fileList
# ==============================================================================
# handing of the initial values
# 1) get the initial values before any averading is done
initFile='initValues.nc'
call "$CDO -selname,$Temp_VarName,$Salinity_VarName -seltimestep,1 ${fileListArray[0]} $initFile"

# 2) compute rhopot and add it to the initial values file
    initRhopotFile='initRhopot.nc'
call "$CDO rhopot,0 $initFile $initRhopotFile"
initWithPotDensity='initValuesWithPotDensity.nc'
call "$CDO merge $initFile $initRhopotFile $initWithPotDensity"

# ==============================================================================
# creating a mask file
maskFile='mask.nc'
call "$CDO -selname,$MaskVarName ${fileListArray[0]} $maskFile"

# ==============================================================================
# Loop over all files in serial
for file in $fileList; do
  baseFilename=$(basename $file)
  # compute yearmean of temperature and salinity and then maskout the land points
  call "$CDO -div -yearmean -selname,$Temp_VarName,$Salinity_VarName $file $maskFile ${Temp_VarName}-${Salinity_VarName}_${baseFilename}"
  # compute the corresponding potential density
  call "$CDO rhopot,0 ${Temp_VarName}-${Salinity_VarName}_${baseFilename} ${PotDensityVarName}_${baseFilename}"
  # merge boith together
  call "$CDO merge ${Temp_VarName}-${Salinity_VarName}_${baseFilename} ${PotDensityVarName}_${baseFilename} merged_${baseFilename}"
  # substract the initial values from it
  call "$CDO sub merged_${baseFilename} $initWithPotDensity diff2init_${baseFilename}"
  # compute the fldmean
  call "$CDO fldmean diff2init_${baseFilename} fldmean_${baseFilename}"
done

# ==============================================================================
# Cat the files together
call "$CDO cat fldmean_${fileListPattern}  $outputDataFile"

# ==============================================================================
# Plot a  hovmoeller type graph
for varname in $Temp_VarName $Salinity_VarName $PotDensityVarName; do 
  call "nclsh $ICONPLOT  -varName=$varname -iFile=$outputDataFile -oFile=${varname}_$(basename $outputDataFile .nc) -oType=png -isIcon -DEBUG -hov=true"
done
