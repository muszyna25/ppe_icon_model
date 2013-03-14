#!/bin/ksh
#------------------------------------------------------------------------------
# generate jsbach initial file for icon grid
#------------------------------------------------------------------------------
set -e

grid=iconR2B04-ocean_etopo40_planet
#grid=iconR2B05-ocean_etopo05_planet

ntiles=11       # currently only works for 11 tiles
year=1850

gridfile=/pool/data/ICON/ocean_data/ocean_grid/${grid}.nc
jsbachfile=/pool/data/JSBACH/T127/jsbach_T127TP04_${ntiles}tiles_${year}.nc
soilpool=/pool/data/JSBACH

filename=${grid}_jsbach_${ntiles}tiles_${year}

# preparations
#--------------

# generate remapping matrices if not yest available
[[ -f rmp_t127_to_${grid} ]] || cdo gendis,${gridfile} ${jsbachfile} rmp_t127_to_${grid}
[[ -f rmp_0.5_to_${grid} ]]  || cdo gendis,${gridfile} ${soilpool}/bclapp.nc rmp_0.5_to_${grid}

# land sea mask on Gaussian grid
cdo selvar,slm ${jsbachfile} slm_t127.nc

# glacier mask on Gaussian grido
cdo selvar,glac ${jsbachfile} glac_t127.nc

# land sea mask on icon grid
cdo gtc,0 -selvar,cell_sea_land_mask ${gridfile} slm_${grid}.nc


rm -f ${filename}_*.tmp
varlist="$(cdo showvar ${jsbachfile}) bclapp fieldcap heatcapacity heatcond hydcond \
         moisture_pot poresize porosity soildepth wiltpoint"
#vg varlist="glac cover_type cover_fract"

for var in ${varlist}; do

    # extrapolate land values to the ocean on the source grid
    case ${var} in

        init_moist | roughness_length | orography_std_dev | fao )
            cdo -fillmiss -setmissval,0 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_t127_to_${grid} ;;
        albedo )
            cdo -fillmiss -setmissval,0.07 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_t127_to_${grid} ;;
        maxmoist )
            cdo -fillmiss -setmissval,1.e-13 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_t127_to_${grid} ;;
        bclapp | fieldcap | heatcapacity | heatcond | hydcond | moisture_pot | poresize | \
            porosity | soildepth | wiltpoint )
            cdo -fillmiss -setmissval,-9999 ${soilpool}/${var}.nc gauss_${var}.tmp
            rmp=rmp_0.5_to_${grid} ;;
        cover_fract | natural_veg )
            cdo splitlevel -selvar,${var} ${jsbachfile} ${var}
            cdo sub ${var}000001.nc glac_t127.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc
            rm -f gauss_${var}.tmp1
            cdo merge ${var}0000??.nc gauss_${var}.tmp1
            rm ${var}0000??.nc
            cdo setvar,${var} -fillmiss -ifnotthen glac_t127.nc -ifthen slm_t127.nc \
                gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_t127_to_${grid} ;;
        cover_type )
            cdo splitlevel -selvar,${var} ${jsbachfile} ${var}
            cdo add ${var}000001.nc glac_t127.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc
            rm -f gauss_${var}.tmp1
            cdo merge ${var}0000??.nc gauss_${var}.tmp1
            rm ${var}0000??.nc
            cdo setvar,${var} -fillmiss -ifthen slm_t127.nc gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_t127_to_${grid} ;;
        forest_fract | glac | snow )
            cdo setvar,${var} -fillmiss -ifthen slm_t127.nc -selvar,${var} ${jsbachfile} \
		gauss_${var}.tmp
            rmp=rmp_t127_to_${grid} ;;
        lai_clim | veg_fract | veg_ratio_max | roughness_length_oro )
            cdo setvar,${var} -fillmiss -sub -setmissval,-9e33 -setmisstoc,-9e33 -setmissval,0 \
                         -add slm_t127.nc -selvar,${var} ${jsbachfile} \
                         slm_t127.nc gauss_${var}.tmp
            rmp=rmp_t127_to_${grid} ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir | surf_temp )
            cdo selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_t127_to_${grid} ;;
    esac

    # do the remapping
    case ${var} in

 	elevation )
            cdo setvar,elevation -ifthen slm_${grid}.nc -selvar,cell_elevation ${gridfile} \
                       ${filename}_elevation.tmp ;;
	slf ) 
            rm -f ${filename}_${var}.tmp1 ${filename}_${var}.tmp ;;
	slm ) 
            cdo setvar,slm slm_${grid}.nc ${filename}_slm.tmp ;;
        * )
            cdo remap,${gridfile},${rmp} gauss_${var}.tmp ${filename}_${var}.tmp1
 	    cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc ${filename}_${var}.tmp1 \
		${filename}_${var}.tmp ;;
    esac
done

# treatment of glaciers
cdo gec,0.5 ${filename}_glac.tmp glac_${grid}.nc
cdo ltc,0.5 ${filename}_glac.tmp glac_inv_${grid}.nc
 
for var in ${varlist}; do
    case ${var} in

       cover_fract | natural_veg )
           cdo splitlevel ${filename}_${var}.tmp ${filename}_${var}_
	   typeset -Z2 nt=01
	   while [[ ${nt} -le ${ntiles} ]]; do
	       cdo mul ${filename}_${var}_0000${nt}.nc2 glac_inv_${grid}.nc \
		   ${filename}_${var}_${nt}.tmp
	       rm ${filename}_${var}_0000${nt}.nc2
	       (( nt = nt + 1 ))
	   done
	   cdo add ${filename}_${var}_01.tmp glac_${grid}.nc ${filename}_${var}_01.tmp1
           mv ${filename}_${var}_01.tmp1 ${filename}_${var}_01.tmp
           rm -f ${filename}_${var}.tmp
           cdo merge ${filename}_${var}_??.tmp ${filename}_${var}.tmp
           cdo -f nc setvar,${var} ${filename}_${var}.tmp ${filename}_${var}.tmp2
           rm ${filename}_${var}_??.tmp ;; 
       cover_type )
           cdo splitlevel ${filename}_${var}.tmp ${filename}_${var}_
	   typeset -Z2 nt=01
           cdo sub ${filename}_${var}_000001.nc2 glac_${grid}.nc ${filename}_${var}_01.tmp1
           mv ${filename}_${var}_01.tmp1 ${filename}_${var}_000001.nc2
           cdo nint ${filename}_${var}_000011.nc2 ${filename}_${var}_11.tmp1
           mv ${filename}_${var}_11.tmp1 ${filename}_${var}_000011.nc2
           rm -f ${filename}_${var}.tmp
           cdo merge ${filename}_${var}_0000??.nc2 ${filename}_${var}.tmp2
           rm ${filename}_${var}_0000??.nc2 ;; 
       glac )
	    cp glac_${grid}.nc ${filename}_${var}.tmp2 ;;
       lai_clim | bclapp | fieldcap | forest_fract | heatcapacity | heatcond | \
	   hydcond | init_moist | maxmoist | moisture_pot | poresize | porosity | \
	   snow | soildepth | veg_fract | veg_ratio_max | wiltpoint )
            cdo mul ${filename}_${var}.tmp glac_inv_${grid}.nc ${filename}_${var}.tmp2 ;;
       albedo )
            cdo setvar,${var} -ifthenelse glac_inv_${grid}.nc ${filename}_${var}.tmp \
		-mulc,0.7 glac_${grid}.nc ${filename}_${var}.tmp2 ;;
    esac
    [[ -f ${filename}_${var}.tmp2 ]] && mv ${filename}_${var}.tmp2  ${filename}_${var}.tmp
done

# generate three separate jsbach initial files, depending on the third dimension
for file in $(ls ${filename}_*.tmp); do

  # find out variables with time dimension
  if [[ $(ncdump -h ${file} | grep "time =") != "" ]]; then
    mv ${file} ${file}.time
  
  # find out variables with tile dimension
  elif [[ $(ncdump -h ${file} | grep "ntiles =") != "" ]]; then
    mv ${file} ${file}.ntiles

  # variables with two dimensions, only
  else
    mv ${file} ${file}.rest
  fi
done

# define uuid attribute
uuid=$(ncdump -h ${gridfile} | grep ':uuid = ' | cut -f2 -d'"')

rm -f ${grid}_jsbach_clim.nc ${grid}_jsbach_${ntiles}tiles_${year}.nc ${grid}_jsbach.nc
if [[ $(ls ${filename}_*.time | wc -l) != 0 ]]; then
  cdo merge ${filename}_*.time   ${grid}_jsbach_clim.nc
  rm  ${filename}_*.time
  ncatted -a uuid,global,o,c,${uuid} ${grid}_jsbach_clim.nc
fi
if [[ $(ls ${filename}_*.ntiles | wc -l) != 0 ]]; then
  cdo merge ${filename}_*.ntiles ${grid}_jsbach_${ntiles}tiles_${year}.nc
  rm ${filename}_*.ntiles
  ncatted -a uuid,global,o,c,${uuid} ${grid}_jsbach_${ntiles}tiles_${year}.nc
fi
if [[ $(ls ${filename}_*.rest | wc -l) != 0 ]]; then
  cdo merge ${filename}_*.rest ${grid}_jsbach.nc
  rm ${filename}_*.rest
  ncatted -a uuid,global,o,c,${uuid} ${grid}_jsbach.nc
fi

# clean up
rm  gauss_*.tmp
rm  ${filename}_*.tmp1
rm  slm_${grid}.nc  slm_t127.nc
rm  glac_t127.nc glac_${grid}.nc glac_inv_${grid}.nc

# test: remap to regular grid to check data with ncview
[[ -f rmp_${grid}_to_t255 ]] || cdo gencon,t255grid ${filename}.nc rmp_${grid}_to_t255
if [[ -f ${grid}_jsbach_clim.nc ]]; then
  cdo remap,t255grid,rmp_${grid}_to_t255 ${grid}_jsbach_clim.nc test_${grid}_jsbach_clim.nc
fi
if [[ -f ${grid}_jsbach_${ntiles}tiles_${year}.nc ]]; then
  cdo remap,t255grid,rmp_${grid}_to_t255 ${grid}_jsbach_${ntiles}tiles_${year}.nc \
                                    test_${grid}_jsbach_${ntiles}tiles_${year}.nc
fi
if [[ -f ${grid}_jsbach.nc ]]; then
  cdo remap,t255grid,rmp_${grid}_to_t255 ${grid}_jsbach.nc test_${grid}_jsbach.nc
fi
