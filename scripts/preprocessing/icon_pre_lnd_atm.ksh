#!/bin/ksh
#------------------------------------------------------------------------------
# Script to construct from a given ICON grid file
# initial fields and boundary conditions consistent
# for the land model JSBACH and the ICON atmosphere.
# 1. compile and run jsbach_init_file.f90 via the script jsbach_init_file.ksh to construct
# a JSBACH input file on a gaussian grid, which can also be used for ECHAM6
# 2. run gauss_to_icon.ksh to interpolate the fields to the ICON grid
#
# jsbach_init_file.f90 was coded by Thomas Raddatz and Veronika Gayler at MPI-Met
# gauss_to_icon.ksh combines a script coded by Veronika Gayler (for land fields) and one
# coded by Marco Giorgetta (for atmosphere fields).
#
# Thomas Raddatz 4. April 2013

set -e

# control fortran compiler (use version 5.3.907)
#module load nag/5.3.907
# control cdo version (use version 1.5.9)
#module load cdo/1.5.9

# get the script jsbach_init_file.ksh and the program jsbach_init_file.f90
svn checkout -r 6285 https://svn.zmaw.de/svn/cosmos/branches/cosmos-landveg/contrib/initial_tarfiles
mv ./initial_tarfiles/jsbach_init_file.ksh .
mv ./initial_tarfiles/jsbach_init_file.f90 .

# switch off interactive mode in the script jsbach_init_file.ksh
export interactive_ctrl=1

# component model resolutions
export res_atm=T63         # atmosphere horizontal resolution
export res_oce=GR15        # ocean horizontal resolution

# ICON grid file
pool_grid=/pool/data/ICON/grids/private/r2b4_amip # pool with ICON grid data
export grid=r2b4_amip
export gridfile=${pool_grid}/${grid}.nc

# echam6 input data
export pool_atm=/pool/data/ECHAM6/input/r0002  # pool with atm input data
export pool_srf=/pool/data/JSBACH         # pool with srf input data

# script directory
startdir=$(pwd)

# parameters for jsbach_init_file.f90
export dynveg=true
export year_ct=1992; export year_cf=1992
export ntiles=11
export c3c4crop=true
export cmip5_pasture=true
export l5_soil_layers=true
export lroot_depth=${l5_soil_layers}
export echam_fractional=true
export srcdir=${startdir}
export wrkdir=${startdir}

# atmosphere input
mkdir -p input/echam6
cd input/echam6
cp -p ${pool_atm}/${res_atm}/${res_atm}${res_oce}_VGRATCLIM.nc  .
cp -p ${pool_atm}/${res_atm}/${res_atm}${res_oce}_VLTCLIM.nc    .
cp -p ${pool_atm}/${res_atm}/${res_atm}${res_oce}_jan_surf.nc   .
cp -p ${pool_atm}/${res_atm}/${res_atm}_TSLCLIM2.nc             .
cd ../..

# land (JSBACH) input
mkdir -p ${srcdir}/input/jsbach
cd ${srcdir}/input/jsbach

${srcdir}/jsbach_init_file.ksh

# remove the executable and modules
rm jsbach_init_file mo_kinds.mod mo_vegparams.mod

# ICON input
mkdir -p ${srcdir}/input/icon
cd ${srcdir}/input/icon

${srcdir}/gauss_to_icon.ksh

exit
