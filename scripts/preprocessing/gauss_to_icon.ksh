#!/bin/ksh
#------------------------------------------------------------------------------
# generate jsbach initial boundary conditions and initial conditions for icon grid
# Veronika Gayler, MPI-M, Hamburg, 2013
#------------------------------------------------------------------------------
set -e

year=${year_cf}
[[ ${echam_fractional} = true ]] && fractional=fractional_ || fractional=""

# input files
jsbachfile=${srcdir}/input/jsbach/jsbach_${res_atm}${res_oce}_${fractional}${ntiles}tiles_${year}.nc
soilpool=/pool/data/JSBACH
ifile_topo=${srcdir}/input/echam6/${res_atm}${res_oce}_jan_surf.nc
filename=${grid}_jsbach_${ntiles}tiles_${year}

# preparations
#--------------

# generate remapping matrices if not yet available
[[ -f rmp_${res_atm}_to_${grid} ]] || cdo gendis,${gridfile} ${jsbachfile} rmp_${res_atm}_to_${grid}
[[ -f rmp_0.5_to_${grid} ]]  || cdo gendis,${gridfile} ${soilpool}/bclapp.nc rmp_0.5_to_${grid}

# land sea mask on Gaussian grid
cdo selvar,slm ${jsbachfile} slm_${res_atm}.nc

# glacier mask on Gaussian grido
cdo selvar,glac ${jsbachfile} glac_${res_atm}.nc

# land sea mask on icon grid
cdo gtc,0 -selvar,cell_sea_land_mask ${gridfile} slm_${grid}.nc

rm -f ${filename}_*.tmp
varlist="$(cdo showvar ${jsbachfile}) bclapp fieldcap heatcapacity heatcond hydcond \
         moisture_pot poresize porosity soildepth wiltpoint"

for var in ${varlist}; do

    # extrapolate land values to the ocean on the source grid
    case ${var} in

        init_moist | roughness_length | orography_std_dev | fao )
            cdo -fillmiss -setmissval,0 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        albedo )
            cdo -fillmiss -setmissval,0.07 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        maxmoist )
            cdo -fillmiss -setmissval,1.e-13 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        bclapp | fieldcap | heatcapacity | heatcond | hydcond | moisture_pot | poresize | \
            porosity | soildepth | wiltpoint )
            cdo -fillmiss -setmissval,-9999 ${soilpool}/${var}.nc gauss_${var}.tmp
            rmp=rmp_0.5_to_${grid} ;;
        cover_fract | natural_veg )
            cdo splitlevel -selvar,${var} ${jsbachfile} ${var}
            cdo sub ${var}000001.nc glac_${res_atm}.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc
            rm -f gauss_${var}.tmp1
            cdo merge ${var}0000??.nc gauss_${var}.tmp1
            rm ${var}0000??.nc
            cdo setvar,${var} -fillmiss -ifnotthen glac_${res_atm}.nc -ifthen slm_${res_atm}.nc \
                gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_${res_atm}_to_${grid} ;;
        cover_type )
            cdo splitlevel -selvar,${var} ${jsbachfile} ${var}
            cdo add ${var}000001.nc glac_${res_atm}.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc
            rm -f gauss_${var}.tmp1
            cdo merge ${var}0000??.nc gauss_${var}.tmp1
            rm ${var}0000??.nc
            cdo setvar,${var} -fillmiss -ifthen slm_${res_atm}.nc gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_${res_atm}_to_${grid} ;;
        forest_fract | glac | snow )
            cdo setvar,${var} -fillmiss -ifthen slm_${res_atm}.nc -selvar,${var} ${jsbachfile} \
		gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        lai_clim | veg_fract | veg_ratio_max | roughness_length_oro )
            cdo setvar,${var} -fillmiss -sub -setmissval,-9e33 -setmisstoc,-9e33 -setmissval,0 \
                         -add slm_${res_atm}.nc -selvar,${var} ${jsbachfile} \
                         slm_${res_atm}.nc gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir | surf_temp )
            cdo selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
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
        fao )
            cdo remapnn,${gridfile} gauss_${var}.tmp ${filename}_${var}.tmp1
 	    cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc ${filename}_${var}.tmp1 \
		${filename}_${var}.tmp ;;
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

#------------------------------------------------------------------------------
# generate boundary conditions depending on topography for icon grid
# Marco Giorgetta, MPI-M, Hamburg 2013-03-08
#------------------------------------------------------------------------------

# 1. remapping ECHAM data from the gaussian grid to the ICON grid at cell centers
# ------------------------------------------------------------------------------
# method = remapnn = nearest neighbor
#

cdo remapnn,${gridfile} ${ifile_topo} jan_surf_tmp.nc

# 2. land sea mask on the ICON grid
# ---------------------------------

# integer sea land mask from the icon grid file
# - ocean cell without land edge : -2
# - ocean cell with land edge    : -1
# - land cell with ocean edge    :  1
# - land cell without ocean edge :  2
#
cdo -f nc selname,cell_sea_land_mask ${gridfile} integer_slm.nc
#
# floating point sea mask
# - at land : 0.
# - at sea  : 1.
#
cdo -b F64 chname,cell_sea_land_mask,sea -ltc,0 integer_slm.nc ${filename}_sea.tmp
#
# floating point land mask
# - at land : 1.
# - at sea  : 0.
#
cdo -b F64 chname,cell_sea_land_mask,notsea -gtc,0 integer_slm.nc ${filename}_notsea.tmp

# 3. set all land parameters to 0 where at sea
# --------------------------------------------
# set surface land variables to zero at sea
# - at land : keep value
# - at sea  : set to 0
#
cdo mul jan_surf_tmp.nc ${filename}_notsea.tmp jan_surf.nc

# 4. separate variables transfered to new files
# ---------------------------------------------
#
# lake surface fractions
#
cdo chname,ALAKE,lake -selname,ALAKE jan_surf.nc ${filename}_lake.tmp
#
# land fraction (not sea and not lake)
#
cdo chname,notsea,land -sub ${filename}_notsea.tmp ${filename}_lake.tmp ${filename}_land.tmp
#
# surface physical properties
#
cdo chname,AZ0,z0  -selname,AZ0 jan_surf.nc ${filename}_z0.tmp
#
# orography
#
cdo chname,OROSTD,orostd -selname,OROSTD jan_surf.nc ${filename}_orostd.tmp
cdo chname,OROSIG,orosig -selname,OROSIG jan_surf.nc ${filename}_orosig.tmp
cdo chname,OROGAM,orogam -selname,OROGAM jan_surf.nc ${filename}_orogam.tmp
cdo chname,OROTHE,orothe -selname,OROTHE jan_surf.nc ${filename}_orothe.tmp
cdo chname,OROPIC,oropic -selname,OROPIC jan_surf.nc ${filename}_oropic.tmp
cdo chname,OROVAL,oroval -selname,OROVAL jan_surf.nc ${filename}_oroval.tmp
cdo chname,OROMEA,oromea -selname,OROMEA jan_surf.nc ${filename}_oromea.tmp

# generate four separate files for boundary conditions (bc) and one file for initial conditions (ic)

rm -f ${grid}_jsbach_clim.nc ${grid}_jsbach_${ntiles}tiles_${year}.nc ${grid}_jsbach.nc

# define surfaces
cdo merge ${filename}_slm.tmp ${filename}_glac.tmp ${filename}_sea.tmp ${filename}_notsea.tmp ${filename}_lake.tmp ${filename}_land.tmp ${filename}_cover_fract.tmp ${filename}_cover_type.tmp ${filename}_natural_veg.tmp ${filename}_veg_ratio_max.tmp bc_land_frac.nc
rm ${filename}_slm.tmp ${filename}_glac.tmp ${filename}_sea.tmp ${filename}_notsea.tmp ${filename}_lake.tmp ${filename}_land.tmp ${filename}_cover_fract.tmp ${filename}_cover_type.tmp ${filename}_natural_veg.tmp ${filename}_veg_ratio_max.tmp 
# physical properties
cdo merge ${filename}_lai_clim.tmp ${filename}_veg_fract.tmp ${filename}_roughness_length.tmp ${filename}_roughness_length_oro.tmp ${filename}_z0.tmp ${filename}_albedo.tmp ${filename}_albedo_veg_vis.tmp ${filename}_albedo_veg_nir.tmp ${filename}_albedo_soil_vis.tmp ${filename}_albedo_soil_nir.tmp ${filename}_forest_fract.tmp bc_land_phys.nc
rm ${filename}_lai_clim.tmp ${filename}_veg_fract.tmp ${filename}_orography_std_dev.tmp ${filename}_roughness_length.tmp ${filename}_roughness_length_oro.tmp ${filename}_z0.tmp ${filename}_albedo.tmp ${filename}_albedo_veg_vis.tmp ${filename}_albedo_veg_nir.tmp ${filename}_albedo_soil_vis.tmp ${filename}_albedo_soil_nir.tmp ${filename}_forest_fract.tmp
# soil properties
cdo merge ${filename}_fao.tmp ${filename}_maxmoist.tmp ${filename}_bclapp.tmp ${filename}_fieldcap.tmp ${filename}_heatcapacity.tmp ${filename}_heatcond.tmp ${filename}_hydcond.tmp ${filename}_moisture_pot.tmp ${filename}_poresize.tmp ${filename}_porosity.tmp ${filename}_soildepth.tmp ${filename}_wiltpoint.tmp bc_land_soil.nc
rm ${filename}_fao.tmp ${filename}_maxmoist.tmp ${filename}_bclapp.tmp ${filename}_fieldcap.tmp ${filename}_heatcapacity.tmp ${filename}_heatcond.tmp ${filename}_hydcond.tmp ${filename}_moisture_pot.tmp ${filename}_poresize.tmp ${filename}_porosity.tmp ${filename}_soildepth.tmp ${filename}_wiltpoint.tmp
# topography
cdo merge ${filename}_elevation.tmp ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp bc_land_sso.nc
rm ${filename}_elevation.tmp ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp
# initial soil conditions
cdo merge ${filename}_surf_temp.tmp ${filename}_init_moist.tmp ${filename}_snow.tmp ic_land_soil.nc
rm ${filename}_surf_temp.tmp ${filename}_init_moist.tmp ${filename}_snow.tmp
# uuid as an attribute
uuidOfHGrid=$(ncdump -h ${gridfile} | grep ':uuidOfHGrid = ' | cut -f2 -d'"')
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_frac.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_phys.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_soil.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_sso.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} ic_land_soil.nc

# clean up
rm gauss_*.tmp
rm ${filename}_*.tmp1
rm slm_${grid}.nc  slm_${res_atm}.nc
rm glac_${res_atm}.nc glac_${grid}.nc glac_inv_${grid}.nc
rm rmp_0.5_to_${grid} rmp_${res_atm}_to_${grid}
rm jan_surf_tmp.nc jan_surf.nc
rm integer_slm.nc

# test: remap to regular grid to check data with ncview
[[ -f rmp_${grid}_to_t255 ]] || cdo gencon,t255grid bc_land_frac.nc rmp_${grid}_to_t255
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_frac.nc test_bc_land_frac.nc
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_phys.nc test_bc_land_phys.nc
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_soil.nc test_bc_land_soil.nc
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_sso.nc test_bc_land_sso.nc
cdo remap,t255grid,rmp_${grid}_to_t255 ic_land_soil.nc test_ic_land_soil.nc
rm rmp_${grid}_to_t255

exit
