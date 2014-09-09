set -x
                ifileLowRes=$1
                 targetGrid=$2
intermediatLonLatResolution=${3:-r720x360}

#===============================================================================
# USAGE:
#   ./low2highResInit.sh <icon-ocean-output> <targetCellGrid> <intermediatenLonLatResolution>
#
# EXAMPLE:
#   ./low2highResInit.sh input/oce_omip_79km_cell_40levels_conservative_21510131T052000Z.nc_part_1+ r10x10 ../towards20km/cell_grid-OceanOnly_Icos_0079km_etopo05-icon.nc r720x360
#===============================================================================
[[ ! -f $1  ]] || [[ ! -f $2 ]] && echo "input data of target grid cannot be found" && exit 1

# select field from input
timemeanFileName=timmean_$(basename ${ifileLowRes})
cdo chname,t_acc,T,s_acc,S,h_acc,h -timmean -selname,t_acc,s_acc,h_acc ${ifileLowRes} ${timemeanFileName}
# put data in an intermediate lonlat grid
intermediateLonLatFileName=${intermediatLonLatResolution}_${timemeanFileName}
cdo -P 24 -remapycon,r720x360 ${timemeanFileName} ${intermediateLonLatFileName}
# fill the missing values with some interpolated data
nomissFileName=nomiss_${intermediateLonLatFileName}
cdo fillmiss2 ${intermediateLonLatFileName} ${nomissFileName}
# put things onto the target grid
cdo -P 12 remapycon,${targetGrid} ${nomissFileName} highRes_from_lowRes_init.nc

# rename the vertical dimension
ncrename -d depth,level   -v depth,level   -O highRes_from_lowRes_init.nc
ncrename -d depth_2,level -v depth_2,level -O highRes_from_lowRes_init.nc
ncrename -d Time,time     -v Time,time     -O highRes_from_lowRes_init.nc
