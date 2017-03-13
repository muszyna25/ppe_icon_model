#!/bin/bash
#--------------------------------------------------------------------------------------------------
#
# Atmosphere-Land-Ocean YAC Coupled Experiment
#
author_list="Marco Giorgetta, Rene Redler, Reiner Schnur, Stephan Lorenz"
#
#--------------------------------------------------------------------------------------------------
#
# This file describes a coupled Atmosphere-Land-Ocean experiment.
#
# based on atm_amip, ocean_omip_long, and yac_nh_r2b3
#
# The atmosphere is based on the non-hydrostatic dynamics and the ECHAM physics, and is intialized
# from analysis files and using transient boundary conditions for:
# - spectral solar irradiation
# - well mixed greenhouse gases CO2, CH4, N2O, CFCs
# - O3 concentration
# - tropospheric aerosol optical properties
# - stratospheric volcanic aerosol optical properties
#
# SST and sea ice are provided via coupling by the ocean
#
# The land ... JSBACH_lite ...
#
# The ocean is initialized with temperature and salinity from PHC3.0 
#
#--------------------------------------------------------------------------------------------------
# I. General setup
#--------------------------------------------------------------------------------------------------
#
# I.1 Variables provided or required by the scripting mechanism
# -------------------------------------------------------------
#
# EXPNAME                       = name of exp. in 'exp.<name>'
# basedir                       = base directory, where src/, run/ etc exist
# nproma                        = blocking length for array dimensioning and inner loop lengths
# second, minute, hour, day     = length of these intervals in [s]
#
# modelname_list[0]="atmo"  : needs to correspond with yac_fdef_comp and xml cmomponent name
# modelname_list[1]="ocean" : needs to correspond with yac_fdef_comp and xml cmomponent name
#
#--------------------------------------------------------------------------------------------------

# I.2 Set variables needed by the scripting mechanism
# ---------------------------------------------------

# Number of ocean processes
# -------------------------
mpi_oce_nodes=${mpi_oce_nodes:=((no_of_nodes/2))}   # default: half of requested nodes
((mpi_oce_procs=mpi_oce_nodes * mpi_procs_pernode))
if [ $mpi_oce_procs -lt 1 ] ; then
   mpi_oce_procs=1                                  # minimum: 1 process
fi

#-----------------------------------------------------------------------------

if [ $mpi_total_procs -lt `expr $mpi_oce_procs + 1` ] ; then
   echo "The coupled runs require at least 2 mpi procs. Exiting."
   check_error 0
   exit
fi

# horizontal grid(s)
# ------------------

# Variables atmo_dyn_grids and ocean_grids are used in namelist and linked to model directory
grids_folder=${icon_data_poolFolder}/coupled_input_temp/aloy_pre04
atmo_dyn_grids='ATMOOceGridWithCoast_158km.nc'
ocean_grids='OCEANGridNoLake_pre04_slo_0039km_bathy.edit2016.nc'

# start and end date+time
# -----------------------

start_date="1979-01-01T00:00:00Z"        # ISO-format date+time
  end_date="1979-01-03T00:00:00Z"        # ISO-format date+time

# time step length
# ----------------

# in ISO-format for atmosphere and ocean submodels:
mts_atm="PT10M"
mts_oce="PT10M"

# in seconds for YAC - must match those of the submodels:
dtime_atm=${dtime_atm:=600}
dtime_oce=${dtime_oce:=600}

# change timestep of ocean/coupled to multiple of atm. timestep
# dtime_oce must be multiple of dtime_atm
#(( dtime_oce=dtime_atm * 1 ))

# restart interval
# ----------------

checkpoint_interval="P200Y" #  mtime ISO-format

# high frequency output and restart intervals
restart_interval="P1D"
atm_output_interval="PT6H"
lnd_output_interval="PT6H"
oce_output_interval="PT6H"

# mid-frequency output and restart intervals
##  P6M ?
#restart_interval="P1Y"   
#atm_output_interval="P1D"
#lnd_output_interval="P1D"
#oce_output_interval="P1D"

## long-term low frequency output and restart intervals
#restart_interval="P2Y"   
#checkpoint_interval="P1Y"
#atm_output_interval="P1M"
#lnd_output_interval="P1M"
#oce_output_interval="P1M"

# ocean timestep:
#  - mts_oce can be increased during a run (up to 30 minutes tested after 1 year spinup)
#  - change timestep for coupler accordingly: dtime_oce
#  - dtime_oce must be multiple of dtime_atm
#mts_oce="PT30M"
#(( dtime_oce=dtime_atm * 3 ))

# standard long-term output switches - _min=yes only:
output_atm_min=yes
output_oce_min=yes
output_lnd_min=yes

atm_file_interval="P20Y"
lnd_file_interval="P20Y"
oce_file_interval="P20Y"

# full output switches - monthly atmospheric output:
atm_output_interval_3d="P1M"
atm_output_interval_2d="P1M"

# full output switches - set to "yes" if needed
output_atm_3d=no
output_atm_2d=no
output_phy_3d=no
output_aer_3d=no # "yes" needs lrad_aero_diag=.TRUE. in radiation_nml
output_oce_full=no
output_lnd_full=no

# coupling timestep for YAC (equals ocean timestep)
dtime_cpl=$dtime_oce

# namelist files
# --------------

atm_namelist=NAMELIST_${EXPNAME}_atm
lnd_namelist=NAMELIST_${EXPNAME}_lnd
oce_namelist=NAMELIST_${EXPNAME}_oce
cpl_xmlctrfn=XMLCTRFN_${EXPNAME}_cpl

# JSBACH settings
run_jsbach=yes
jsbach_with_hd=yes
jsbach_usecase=jsbach_lite
#
[[ $jsbach_with_hd == yes ]] && jsbach_usecase=${jsbach_usecase}_with_hd
ljsbach=$([ "${run_jsbach:=no}" == yes ] && echo .TRUE. || echo .FALSE. )

#--------------------------------------------------------------------------------------------------
# II. Atmosphere & land section
#--------------------------------------------------------------------------------------------------

# II.1 Define the atmosphere and land model configuration
# -------------------------------------------------------

# II.1.1 Atmospheric dynamics and physics
# ---------------------------------------

cat > ${atm_namelist} << EOF
&coupling_mode_nml
  coupled_mode = .TRUE.
/
&parallel_nml
 nproma           = ${nproma}
/
&grid_nml
 dynamics_grid_filename = "${atmo_dyn_grids}",
/
&run_nml
 num_lev          = 47           ! number of full levels
 modelTimeStep    = "$mts_atm"   ! model time step in ISO-format
 ltestcase        = .FALSE.      ! run testcase
 ldynamics        = .TRUE.       ! dynamics
 ltransport       = .TRUE.       ! transport
 ntracer          = 3            ! number of tracers
 iforcing         = 2            ! 0: none, 1: HS, 2: ECHAM, 3: NWP
 output           = 'nml'
 msg_level        = 15           ! level of details report during integration 
 restart_filename = "${EXPNAME}_restart_atm_<rsttime>.nc"
 activate_sync_timers = .TRUE.
 profiling_output     = 2

/
&extpar_nml
 itopo            = 1           ! 1: read topography from the grid file
 l_emiss          = .FALSE.
/
&initicon_nml
 init_mode        = 2           ! 2: initialize from IFS analysis
/
&dynamics_nml
 iequations       = 3           ! 3: ICONAM dynamics
/
&nonhydrostatic_nml
 ndyn_substeps    = 5           ! dtime/dt_dyn
 damp_height      = 50000.      ! [m]
 rayleigh_coeff   = 0.10
 vwind_offctr     = 0.2
 divdamp_fac      = 0.004
/
&interpol_nml
 rbf_scale_mode_ll = 1
/
&sleve_nml
 min_lay_thckn    = 40.         ! [m]
 top_height       = 83000.      ! [m]
 stretch_fac      = 0.9
 decay_scale_1    = 4000.       ! [m]
 decay_scale_2    = 2500.       ! [m]
 decay_exp        = 1.2
 flat_height      = 16000.      ! [m]
/
&diffusion_nml
/
&transport_nml
 ctracer_list     = 'vwi'       ! water vapour, cloud water, cloud ice
 ivadv_tracer     = 3,3,3
 itype_hlimit     = 3,4,4
 ihadv_tracer     = 52,2,2
/
&echam_phy_nml
 dt_rad           =  7200.      ! [s] radiation time step
 lrad             = .TRUE.
 lvdiff           = .TRUE.
 lconv            = .TRUE.
 lcond            = .TRUE.
 lgw_hines        = .TRUE.
 lssodrag         = .TRUE.
 lice             = .TRUE.
 ljsbach          = .TRUE.
 lebudget         = .TRUE.
 lmlo             = .FALSE.
 lamip            = .FALSE.
/
&radiation_nml
 irad_h2o         = 1           ! 1: prognostic vapor, liquid and ice
 irad_co2         = 2           ! 2: global constant; 4: from greenhouse gas scenario
 irad_ch4         = 2           ! 2: global constant; 4: from greenhouse gas scenario
 irad_n2o         = 2           ! 2: global constant; 4: from greenhouse gas scenario
 irad_o3          = 8           ! 8: horizontally and vertically variable
 irad_o2          = 2           ! 2: horizontally and vertically constant
 irad_cfc11       = 0           ! 0: no cfc
 irad_cfc12       = 0           ! 0: no cfc
 irad_aero        = 0           ! 0: no aerosol
 lrad_aero_diag   = .FALSE.     ! switch for diagnostics of the aerosol optical properties
 ighg             = 0           ! 0: select default values; 1: transient well mixed greenhouse gas concentrations
 izenith          = 4           ! 4: seasonal and diurnal cycle
 isolrad          = 6           ! 6: pre-industrial CMIP6 solar constant
 vmr_co2          = 284.317e-6  ! constant volume mixing ration
 vmr_ch4          = 808.249e-9  ! constant volume mixing ration
 vmr_n2o          = 273.021e-9  ! constant volume mixing ration
/
&psrad_nml
rad_perm          = 1           ! Integer for perturbing random number seeds
/
&echam_conv_nml
/
&echam_cloud_nml
/
&gw_hines_nml
/
&sea_ice_nml
/
EOF

# II.1.2 Land surface and soil
# ----------------------------

cat > ${lnd_namelist} <<EOF
&jsb_model_nml
  usecase         = "${jsbach_usecase}"
/
&jsb_srf_nml
  bc_filename     = 'bc_land_phys.nc'
  bc_sso_filename = 'bc_land_sso.nc'
  ic_filename     = 'ic_land_soil.nc'
/
&jsb_soil_nml
  active          = .TRUE.
  nsoil_energy    = 5
  nsoil_water     = 5
  bc_filename     = 'bc_land_soil.nc'
  ic_filename     = 'ic_land_soil.nc'
/
&jsb_veg_nml
  active          = .TRUE.
  bc_filename     = 'bc_land_phys.nc'
  ic_filename     = 'ic_land_soil.nc'
/
EOF
if [[ ${jsbach_with_hd} = yes ]]; then
cat >> ${lnd_namelist} <<EOF
&jsb_hd_nml
  active          = .TRUE.
  routing_scheme  = 'full'
  bc_filename     = 'bc_land_hd.nc'
  diag_water_budget = .TRUE.
  debug_hd        = .FALSE.
/
EOF
fi

#--------------------------------------------------------------------------------------------------

# II.2 Define the input
# ---------------------

# II.2.0 Grid files
# ------------------
# slo.2016/09 atmospere and ocean base masks, no longer copied from HGRIDDIR - not yet tested
#add_link_file $grids_folder/$atmo_dyn_grids                             ./$atmo_dyn_grids                
#add_link_file $grids_folder/$ocean_grids                                ./$ocean_grids                

# II.2.1 Model files
# ------------------

add_link_file ${basedir}/data/rrtmg_lw.nc                               ./
add_link_file ${basedir}/data/rrtmg_sw.nc                               ./
add_link_file ${basedir}/data/ECHAM6_CldOptProps.nc                     ./

# II.2.2 Namelist files
# ---------------------

add_required_file ${basedir}/run/${atm_namelist}                        ./
add_required_file ${basedir}/run/${lnd_namelist}                        ./

# II.2.3 Dictionary file for output variable names
# ------------------------------------------------

dict_file="dict.${EXPNAME}"
cat dict.iconam.mpim > ${dict_file}
add_required_file ${basedir}/run/${dict_file}                           ./

# II.2.4 Initial conditions
# -------------------------

INDATA=${icon_data_poolFolder}/input/r0003
#
# - atmosphere: ECMWF analysis, 1979-01-01T00:00:00Z
datadir=/pool/data/ICON/setup/ifs_iconremap_amip
add_link_file $datadir/ifs_remap_R2B4_00XX_AMIP_2012010100_setyear1979.nc ./ifs2icon_R2B04_DOM01.nc
#
# - land: source?, date+time?
datadir=$grids_folder/Land
add_link_file ${datadir}/ic_land_soil_1976.nc                           ./ic_land_soil.nc

# II.2.5 Boundary conditions
# --------------------------

# range of years for yearly files
# assume start_date and end_date have the format yyyy-...
start_year=$(( ${start_date%%-*} - 1 ))
end_year=$(( ${end_date%%-*} + 1 ))
#
# - well mixed greenhouse gases
# mbe/slo - attention: do not use, set to constant 1850 via namelist
datadir=$INDATA/global/atm
add_link_file $datadir/bc_greenhouse_rcp45_1765-2500.nc                 ./bc_greenhouse_gases.nc

# - ozone
# mbe/slo - attention: do not use, set to constant 1850 via namelist
datadir=$grids_folder/Atmos
add_link_file $datadir/bc_ozone_pre04_slo_0158km_CMIP6_picontrol.nc     ./bc_ozone.nc

# - sst and sic
datadir=$grids_folder/Atmos
add_link_file $datadir/bc_sic_pcmdi_cmip6_pre04_slo_0158km_1870-2015.nc ./bc_sic.nc
add_link_file $datadir/bc_sst_pcmdi_cmip6_pre04_slo_0158km_1870-2015.nc ./bc_sst.nc

# - land parameters
datadir=$grids_folder/Land
add_link_file $datadir/bc_land_frac_1976.nc                             ./bc_land_frac.nc
add_link_file $datadir/bc_land_phys_1976.nc                             ./bc_land_phys.nc
add_link_file $datadir/bc_land_soil_1976.nc                             ./bc_land_soil.nc
add_link_file $datadir/bc_land_sso_1976.nc                              ./bc_land_sso.nc
add_link_file $datadir/hd_para_icon_r2b4_pre04_vs1_5.nc                 ./bc_land_hd.nc

# slo.2016/09 read -2to2 cell_sea_land_mask, name hardcoded as r2b4_amip.nc - use same as for HD-model (lsf>0:=land)
add_link_file $datadir/slm-2to2_forHDmodel_158km.nc                     ./r2b4_amip.nc
add_link_file $datadir/slm-2to2_forHDmodel_158km.nc                     ./hd_mask.nc
add_link_file $datadir/runoff_a2o.r110.3.HD.corGibr.nc                  ./runoff_a2o.weights.nc
#add_link_file $datadir/runoff_a2o.r110.2.HD.nc                          ./runoff_a2o.weights.nc
#

#--------------------------------------------------------------------------------------------------

# II.3 Define the output
# ----------------------

# II.3.1 Parameters for all output files
# --------------------------------------

cat >> ${atm_namelist} << EOF
&io_nml
 output_nml_dict  = "${dict_file}"
 netcdf_dict      = "${dict_file}"
 itype_pres_msl   = 4
 lkeep_in_sync = .TRUE.          ! sync after each timestep
/
EOF

# II.3.3 Define output files
# --------------------------
#
# output_<xyz>=yes : yes --> output files for <xyz>, any other value --> no files for <xyz>

if [[ "$output_atm_min" == "yes" ]]; then
  cat >> ${atm_namelist} <<EOF
&output_nml
 output_filename  = "${EXPNAME}_atm_min"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 operation        = 'mean'
 output_grid      = .TRUE.
 output_start     = "${start_date}"
 output_end       = "${end_date}"
 output_interval  = "${atm_output_interval}"
 file_interval    = "${atm_file_interval}"
 include_last     = .FALSE.
 ml_varlist       = 'orog'    , 'ps'      ,
                    'frac_wtr', 'frac_ice', 'frac_lnd',
                    'sic'     , 'sit'     , 'qtop_icecl', 'qbot_icecl',
                    'clt'     , 'psl'     , 'tas'     , 'ts'      ,
                    'pr'      , 'prw'     , 'cllvi'   , 'clivi'   ,
                    'prlr'    , 'prls'    , 'prcr'    , 'prcs'    ,
                    'tauu'    , 'tauv'    , 'albedo'  ,
                    'hfls'    , 'hfss'    , 'evspsbl' ,
                    'hfls_wtr', 'hfls_ice', 'hfls_lnd',
                    'hfss_wtr', 'hfss_ice', 'hfss_lnd',
                    'rsns_wtr', 'rsns_ice', 'rsns_lnd',
                    'rlns_wtr', 'rlns_ice', 'rlns_lnd', 
                    'evspsbl_wtr', 'evspsbl_ice', 'evspsbl_lnd',
                    'ts_wtr'  , 'ts_ice'  , 'ts_lnd'  ,
                    'tauu_wtr', 'tauu_ice', 'tauu_lnd',
                    'tauv_wtr', 'tauv_ice', 'tauv_lnd',
                    'sfcwind' , 'uas'     , 'vas'
/
&dbg_index_nml
  idbg_mxmn        = 2                        ! initialize MIN/MAX  debug output
  idbg_val         = 0                        ! initialize one cell debug output
  idbg_slev        = 1                        ! initialize start level for debug output
  idbg_elev        = 2                        ! initialize end level for debug output
  dbg_lat_in       =  30.0                    ! latitude location of one cell debug output
  dbg_lon_in       = -30.0                    ! longitude location of one cell debug output
  str_mod_tst      ='InterFaceOce'            ! define modules to print out in debug mode
/
EOF
fi

#  Variables not allowed?
#  'rsns'    , 'rlns'    , 'rsnt'    , 'rlnt'    , 'rsdt'    ,

#  more variables for test output
#  'cosmu0'  , 'daylght_frc',
#  'visdffsfc','nirdffsfc','vissfc'  , 'nirsfc'  ,
#  'albedo_wtr', 'albedo_ice', 'albedo_lnd',
#  'albvisdir','albvisdif' , 'albnirdir' , 'albnirdif',
#  'albvisdir_ice', 'albvisdir_wtr', 'albvisdir_lnd',
#  'tauu_sso', 'tauv_sso', 'diss_sso', 
#  'sh_vdiff', 'qv_vdiff',
#  'ch_concloud',
#  'con_dtrl', 'con_dtri', 'con_iteqv',
#  'cld_dtrl', 'cld_dtri', 'cld_iteq' ,
#  'dew2'


if [[ "$output_atm_3d" == "yes" ]]; then
  cat >> ${atm_namelist} << EOF
&output_nml
 output_filename  = "${EXPNAME}_atm_3d"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 operation        = 'mean'
 output_grid      = .TRUE.
 output_start     = "${start_date}"
 output_end       = "${end_date}"
 output_interval  = "${atm_output_interval_3d}"
 file_interval    = "${atm_file_interval}"
 include_last     = .FALSE.
 ml_varlist       = 'ta','ua','va','wap','hus','hur','cl','clw','cli','rho',
                    'zg','pfull','ps'
/
EOF
fi


if [[ "$output_atm_2d" == "yes" ]]; then
# minimum list for table and quickplots
  cat >> ${atm_namelist} << EOF
&output_nml
 output_filename  = "${EXPNAME}_atm_2d"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 operation        = 'mean'
 output_grid      = .TRUE.
 output_start     = "${start_date}"
 output_end       = "${end_date}"
 output_interval  = "${atm_output_interval_2d}"
 file_interval    = "${atm_file_interval}"
 include_last     = .FALSE.
 ml_varlist       = 'orog'    , 'ps'      , 'psl'     ,
                    'cosmu0'  ,
                    'rsdt'    ,
                    'rsut'    , 'rsutcs'  , 'rlut'    , 'rlutcs'
                    'rsds'    , 'rsdscs'  , 'rlds'    , 'rldscs'
                    'rsus'    , 'rsuscs'  , 'rlus'    ,
                    'ts'      ,
                    'sic'     , 'sit'     ,
                    'albedo'  ,
                    'clt'     ,
                    'prlr'    , 'prls'    , 'prcr'    , 'prcs'    ,
                    'pr'      , 'prw'     , 'cllvi'   , 'clivi'   ,
                    'hfls'    , 'hfss'    , 'evspsbl' ,
                    'tauu'    , 'tauv'    ,
                    'tauu_sso', 'tauv_sso', 'diss_sso', 
                    'sfcwind' , 'uas'     , 'vas'     ,
                    'tas'     , 'dew2'
/
EOF
fi


if [[ "$output_phy_3d" == "yes" ]]; then
  cat >> ${atm_namelist} << EOF
&output_nml
 output_filename  = "${EXPNAME}_phy_3d"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 operation        = 'mean'
 output_grid      = .TRUE.
 output_start     = "${start_date}"
 output_end       = "${end_date}"
 output_interval  = "${atm_output_interval_3d}"
 file_interval    = "${atm_file_interval}"
 include_last     = .FALSE.
 ml_varlist       = 'tend_ta'      , 'tend_ta_dyn'  , 'tend_ta_phy'  ,
                    'tend_ta_rlw'  , 'tend_ta_rsw'  ,
                    'tend_ta_vdf'  , 'tend_ta_gwh'  , 'tend_ta_sso'  ,
                    'tend_ta_cnv'  , 'tend_ta_cld'  , 
                    'tend_ua'      , 'tend_ua_dyn'  , 'tend_ua_phy'  ,
                    'tend_ua_vdf'  , 'tend_ua_gwh'  , 'tend_ua_sso'  ,
                    'tend_ua_cnv'  , 
                    'tend_va'      , 'tend_va_dyn'  , 'tend_va_phy'  ,
                    'tend_va_vdf'  , 'tend_va_gwh'  , 'tend_va_sso'  ,
                    'tend_va_cnv'  ,
                    'tend_hus'     , 'tend_hus_dyn' , 'tend_hus_phy' ,
                    'tend_hus_cld' , 'tend_hus_cnv' , 'tend_hus_vdf' ,
                    'pfull','ps'
/
EOF
fi


if [[ "$output_aer_3d" == "yes" ]]; then
  cat >> ${atm_namelist} << EOF
&output_nml
 output_filename  = "${EXPNAME}_aer_3d"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 operation        = 'mean'
 output_grid      = .TRUE.
 output_start     = "${start_date}"
 output_end       = "${end_date}"
 output_interval  = "${atm_output_interval_3d}"
 file_interval    = "${atm_file_interval}"
 include_last     = .FALSE.
 ml_varlist       = 'aer_aod_533',  'aer_ssa_533',  'aer_asy_533' , 
                    'aer_aod_2325', 'aer_ssa_2325', 'aer_asy_2325', 
                    'aer_aod_9731',
                    'pfull','ps'
/
EOF
fi


if [[ "$output_lnd_min" == "yes" ]]; then
  cat >> ${atm_namelist} << EOF
&output_nml
 output_filename  = "${EXPNAME}_lnd_min"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 output_grid      = .TRUE.
 output_start     = "${start_date}"                  ! ISO-format date+time
 output_end       = "${end_date}"                    ! ISO-format date+time
 output_interval  = "${lnd_output_interval}"         ! ISO-format interval
 file_interval    = "${lnd_file_interval}"           ! ISO-format interval
 include_last     = .FALSE.
 ml_varlist       = 'fract',
                    't_srf', 'qsat_srf',
                    'lw_srf_down', 'swnir_srf_down', 'swpar_srf_down', 'swvis_srf_down',
                    'root_depth', 'soil_depth',
                    'runoff', 'drainage', 'pme_glacier', 'discharge', 'discharge_ocean'
/
EOF
fi

if [[ "$output_lnd_full" == "yes" ]]; then
  cat >> ${atm_namelist} << EOF
&output_nml
 output_filename  = "${EXPNAME}_lnd_full"
 filename_format  = "<output_filename>_<levtype_l>_<datetime2>"
 remap            = 0
 output_grid      = .TRUE.
 output_start     = "${start_date}"                  ! ISO-format date+time
 output_end       = "${end_date}"                    ! ISO-format date+time
 output_interval  = "${lnd_output_interval}"         ! ISO-format interval
 file_interval    = "${lnd_file_interval}"           ! ISO-format interval
 include_last     = .FALSE.
 ml_varlist       = 'fract',
                    't_srf', 'qsat_srf',
                    'alb_vis_srf', 'alb_nir_srf',
                    'canopy_cond', 'ws_root',
                    'albedo_srf', 't_air', 'q_air', 'ws',
                    'lw_srf_down', 'swnir_srf_down', 'swpar_srf_down', 'swvis_srf_down',
                    'root_depth', 'soil_depth',
                    'evapotrans', 'sensible_hflx', 'latent_hflx' 
                    't_soil', 'wsn_srf', 'wsr_srf', 
                    'sfract_srf', 'wfract_srf', 'sfract_soil', 'sfract_can', 'wfract_can', 'wfract_soil'
                    'rough_m_srf', 'rough_h_srf'
                    'runoff', 'drainage', 'pme_glacier', 'discharge', 'discharge_ocean'
/
EOF
fi

#--------------------------------------------------------------------------------------------------
# III. ocean and sea-ice section
#--------------------------------------------------------------------------------------------------

# III.1 ocean namelist
#---------------------

cat > ${oce_namelist} << EOF
&coupling_mode_nml
  coupled_mode = .TRUE.
/
&parallel_nml
 nproma        = ${nproma}
 p_test_run    = .FALSE.
 l_fast_sum    = .TRUE.
/
&grid_nml
 dynamics_grid_filename      = "${ocean_grids}",
 use_dummy_cell_closure      = .true.
 use_duplicated_connectivity = .false.
/
&dynamics_nml
 iequations  = -1                             ! -1: hydrost. ocean model
/
&run_nml
 modelTimeStep    = "$mts_oce"                ! model time step in ISO-format
 output      = 'nml'                          ! output mechanism via namelist
 activate_sync_timers = .TRUE.
 profiling_output     = 2
 timers_level         = 10
 debug_check_level    = 1
/
&output_nml
  output_start     = "${start_date}"                  ! start date in ISO-format
  output_end       = "${end_date}"                    ! end date in ISO-format
  output_interval  = "${oce_output_interval}"         ! interval in ISO-format
  file_interval    = "${oce_file_interval}"           ! interval in ISO-format
  output_grid      = .TRUE.
  output_filename  = "${EXPNAME}_oceanMonitor"
  filename_format  = "<output_filename>_<datetime2>"
  ml_varlist       =  'group:ocean_monitor'
/
&dbg_index_nml
  idbg_mxmn        = 2                        ! initialize MIN/MAX  debug output
  idbg_val         = 0                        ! initialize one cell debug output
  idbg_slev        = 1                        ! initialize start level for debug output
  idbg_elev        = 2                        ! initialize end level for debug output
  dbg_lat_in       =  30.0                    ! latitude location of one cell debug output
  dbg_lon_in       = -30.0                    ! longitude location of one cell debug output
  str_mod_tst      ='oceanCouplng'            ! define modules to print out in debug mode
/
&ocean_dynamics_nml
! 64 unevenly spaced levels used by MPIOM
 n_zlev             =   64
 dzlev_m(1:64)      =  12.0,   10.0,   10.0,   10.0,   10.0,   10.0,   10.0,   11.0,   11.0,   12.0,
                       13.0,   14.0,   15.0,   16.0,   17.0,   18.0,   20.0,   22.0,   24.0,   26.0,
                       29.0,   33.0,   36.0,   39.0,   42.0,   45.0,   48.0,   52.0,   56.0,   60.0,
                       64.0,   68.0,   72.0,   76.0,   80.0,   85.0,   90.0,   95.0,  100.0,  106.0,
                      112.0,  118.0,  124.0,  130.0,  138.0,  146.0,  154.0,  162.0,  170.0,  178.0,
                      186.0,  196.0,  206.0,  216.0,  226.0,  236.0,  246.0,  256.0,  266.0,  276.0,
                      286.0,  296.0,  306.0,  318.0
  l_edge_based                    = .FALSE.   ! edge- or cell-based mimetic discretization
  select_solver                   =   2       ! 1=gmres_oce_old
                                              ! 2=ocean_restart_gmres
                                              ! 3=mixed precisison restart
  solver_FirstGuess               =   1
  use_absolute_solver_tolerance   = .true.
  solver_tolerance                =   7.5E-14
  solver_max_iter_per_restart     =  17
  solver_tolerance_sp             =   1.0E-12 !
  solver_max_iter_per_restart_sp  =  26       !
  solver_max_restart_iterations   = 100       ! outer (restart solver)
  fast_performance_level          = 200       ! performance level 12: for cell-based; 5: default
  use_continuity_correction       = .TRUE.    ! height adjustment according to vertical velocity in dynamics
  cfl_check                       = .FALSE.
  cfl_write                       = .FALSE.
  i_bc_veloc_top                  =   4
  i_bc_veloc_bot                  =   1
/
&ocean_tracer_transport_nml
  FLUX_CALCULATION_HORZ = 5      ! 1=upwind, 2=central, 3=Lax-Friedrichs, 4=Miura, 5=FCT with Zalesak limiter (default)
  FLUX_CALCULATION_VERT = 7      ! 1=upwind 6=adpo; 7=upwind biased ppm (default); 8=FCT with zalesak limiter
  ! define low and high order methods to be used in horizontal flux corrected transport methods (flux_calculation_horz=4,5)
  fct_low_order_flux              =   1       ! horizontal low  order method: 1=upwind (def), no other implemented
  fct_high_order_flux             =   5       ! horizontal high order method: 1=upwind
                                              !                               2=central (def)
                                              !                               3=lax_friedrichs
                                              !                               4=miura_order1
                                              !                               5=??
  fct_limiter_horz                = 100 
  threshold_min_T                 =  -4.0     ! to avoid abort
/
&ocean_horizontal_diffusion_nml
  laplacian_form = 1                   ! 1=curlcurl-graddiv
  VelocityDiffusion_order = 2          ! 21=biharmonic+laplacian (for the laplacian leith)

  BiharmonicViscosity_scaling     =  1
  BiharmonicViscosity_reference   =  1.5E12  !  [m2/s] constant horizontal viscosity coefficient for velocity
  BiharmonicViscosity_background  =  0.0  ! [m2/s] constant horizontal viscosity coefficient for velocity

  TracerHorizontalDiffusion_scaling          = 1
  Temperature_HorizontalDiffusion_Background = 0.0
  Temperature_HorizontalDiffusion_Reference  = 120 !  40
  Salinity_HorizontalDiffusion_Background    = 0.0
  Salinity_HorizontalDiffusion_Reference     = 120 !  40
/
&ocean_vertical_diffusion_nml
  ! PPscheme_type: type of pp-scheme, 4 is recommended:
  ! 1=ICON-optimized; 2=MPIOM-type; 4=ICON-edge-based - factor e-3 for wind mixing needed
  PPscheme_type                            = 2
  VerticalViscosity_TimeWeight             = 0.005
  velocity_VerticalDiffusion_background    =   5.0E-5  ! [m2/s]  vertical background viscosity coefficient for velocity
  Temperature_VerticalDiffusion_background =   1.0E-5  ! [m2/s]  vertical background diffusion coefficient for temperature
  Salinity_VerticalDiffusion_background    =   1.0E-5  ! [m2/s]  vertical background diffusion coefficient for salinity
  tracer_convection_MixingCoefficient      =   0.1 ! max vertical tracer diffusion for convection used in case of instability
!  convection_InstabilityThreshold         =  -1.0E-6  ! used in update_ho_params - default=-5e-8
!  RichardsonDiffusion_threshold           =   0.0     ! used in update_ho_params - default=+5e-8
  tracer_RichardsonCoeff                   =   2.0E-3  ! factor for vertical diffusion coefficient in PP scheme
  velocity_RichardsonCoeff                 =   2.0E-3  ! factor for vertical viscosity coefficient in PP scheme
  bottom_drag_coeff                        =   3.0E-3  ! default=2.5E-3; active for i_bc_veloc_bot=1

  use_wind_mixing                 = .true.    ! true: use wind mixing scheme in MPIOM-type pp-scheme
  ! lambda_wind                     = 0.03
  tracer_TopWindMixing            = 5.0E-4 ! 4.0E-6 The used values maybe are too high 
  velocity_TopWindMixing          = 5.0E-4 ! 4.0E-6
/
&ocean_GentMcWilliamsRedi_nml
  GMRedi_configuration           =   0       ! 0=cartesian diffusion 1=GM-Redi: bolus advection + isopycnal diffusion
  GMREDI_COMBINED_DIAGNOSTIC    = .false.
/
&ocean_physics_nml
/                                           
&sea_ice_nml
  i_ice_therm                     =   1       ! 1=zero-layer (default), 2=Winton, 0/2: not allowed
  i_ice_dyn                       =   1       ! 1/0=switch on/off AWI ice dynamics; TODO: ice dynamics not yet ready for coupling
/
&ocean_forcing_nml
  iforc_oce                       =  14       ! ocean forcing: FORCING_FROM_COUPLED_FLUX  = 14
  type_surfRelax_Temp             =  -1       ! -1: use net surface heat flux from atmosphere as boundary condition - no relaxation
  forcing_enable_freshwater       = .TRUE.    ! enable freshwater flux
  forcing_windstress_u_type       =  22       ! 0: wind-stress set to zero (default) 1: read from (OMIP-) file; 2-100: untouched wind stress
  forcing_windstress_v_type       =  22
  limit_seaice                    = .TRUE.    ! true: mass and energy conserving maximum of sea ice
  seaice_limit                    =   0.7     ! 12m surface layer: 8.4m sea-ice limit
/
&ocean_initialConditions_nml
  initial_salinity_type           =   1       ! read S from initial_state.nc
  initial_temperature_type        =   1       ! read T from initial_state.nc
/                                    
&ocean_diagnostics_nml
   diagnostics_level              =   1
/
&io_nml
  lkeep_in_sync = .TRUE.          ! sync after each timestep
/
EOF

if [[ "$output_oce_min" == "yes" ]]; then
  cat >> ${oce_namelist} <<EOF

&output_nml
  filetype         =  4                               ! output format: 2=GRIB2, 4=NETCDFv2
  output_filename  = "${EXPNAME}_oce_min"
  filename_format  = "<output_filename>_<datetime2>"
  output_start     = "${start_date}"                  ! start date in ISO-format
  output_end       = "${end_date}"                    ! end date in ISO-format
  output_interval  = "${oce_output_interval}"         ! interval in ISO-format
  file_interval    = "${oce_file_interval}"           ! interval in ISO-format
  output_grid      = .TRUE.
  mode             =  2                               ! 1: forecast mode (relative t-axis),
                                                      ! 2: climate mode (absolute t-axis)
  include_last     = .false.
  m_levels         = "1...3,10"                       ! reduced output out of 64 levels
  filename_format  = "<output_filename>_<datetime2>"
  ml_varlist       =  'wet_c','basin_c','regio_c','lsm_ctr_c', 'Qtop', 'Qbot',
                      'mld','condep','h_acc','u_vint_acc','hi_acc','hs_acc','conc_acc',
                      't_acc','s_acc','u_acc','topBoundCond_windStress_u_acc',
                      'HeatFlux_Total_acc','atmos_fluxes_HeatFlux_ShortWave','atmos_fluxes_HeatFlux_LongWave',
                      'HeatFlux_ShortWave_acc',
                      'atmos_fluxes_HeatFlux_Latent','atmos_fluxes_HeatFlux_Sensible'
                      'FrshFlux_Runoff_acc','FrshFlux_Precipitation_acc','FrshFlux_Evaporation_acc','FrshFlux_SnowFall_acc',
                      'FrshFlux_TotalOcean_acc', 'oceWind_Speed_10m'
/
EOF
fi

# m_levels         = "1...3,10,21,32,44,51,56"        ! reduced output out of 64 levels
#                                                     ! (0-27m, 100m, 305.5, 829, 1986, 3092, 4147m) 

if [[ "$output_oce_full" == "yes" ]]; then
  cat >> ${oce_namelist} <<EOF

&output_nml
  filetype         =  4                               ! output format: 2=GRIB2, 4=NETCDFv2
  output_filename  = "${EXPNAME}_oce_full"
  filename_format  = "<output_filename>_<datetime2>"
  output_start     = "${start_date}"                  ! start date in ISO-format
  output_end       = "${end_date}"                    ! end date in ISO-format
  output_interval  = "${oce_output_interval}"         ! interval in ISO-format
  file_interval    = "${oce_file_interval}"           ! interval in ISO-format
  output_grid      = .TRUE.
  mode             =  2                               ! 1: forecast mode (relative t-axis),
                                                      ! 2: climate mode (absolute t-axis)
  include_last     = .false.
  filename_format  = "<output_filename>_<datetime2>"
  ml_varlist       =  'wet_c','basin_c','regio_c','mld','condep','h_acc','u_vint_acc','hi_acc','hs_acc','conc_acc',
                      'Qbot','Qtop','topBoundCond_windStress_u_acc','t_acc','s_acc','u_acc', 'v_acc', 'w_acc',
                      'HeatFlux_Total_acc','FrshFlux_Runoff_acc','FrshFlux_Precipitation_acc',
                      'FrshFlux_Evaporation_acc', 'FrshFlux_TotalIce_acc', 'FrshFlux_TotalSalt_acc',
                      'FrshFlux_TotalOcean_acc', 'FrshFlux_VolumeTotal_acc', 'FrshFlux_SnowFall_acc',
                      'atmos_fluxes_FrshFlux_VolumeIce','oceWind_Speed_10m'
/
EOF
fi

# III.2 Add standard ocean input files 
# ------------------------------------

# PHC3.0, 40km resolution:
datadir=$grids_folder/Ocean
add_required_file ${datadir}/ts_phc3.0_annual_pre04_slo_0039km_L64.nc initial_state.nc

#--------------------------------------------------------------------------------------------------
# IV. coupling section
#--------------------------------------------------------------------------------------------------

if [ $mpi_total_procs -lt 2 ] ; then
  check_error 0 "This setup requires at least 2 mpi processes. Exit"
fi

# IV.1 Split the number of procs in two for each of the dummy component
# ---------------------------------------------------------------------
oce_min_rank=`expr ${mpi_total_procs} - ${mpi_oce_procs}`
oce_max_rank=`expr ${oce_min_rank} + ${mpi_oce_procs} - 1`
oce_inc_rank=1
atm_min_rank=0
atm_max_rank=`expr ${oce_min_rank} - 1`
atm_inc_rank=1
#
# IV.2 Fill model list
# --------------------
#
namelist_list[0]="$atm_namelist"
modelname_list[0]="atmo"
modeltype_list[0]=1
minrank_list[0]=$atm_min_rank
maxrank_list[0]=$atm_max_rank
incrank_list[0]=$atm_inc_rank
#
namelist_list[1]="$oce_namelist"
modelname_list[1]="ocean"
modeltype_list[1]=2
minrank_list[1]=$oce_min_rank
maxrank_list[1]=$oce_max_rank
incrank_list[1]=$oce_inc_rank

# IV.3 Write YAC coupling xml file
# --------------------------------
#
atm_lag=1  #  time lag of counting of timesteps in atmosphere
#
cat > ${cpl_xmlctrfn} << EOF
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<coupling xmlns="http://www.w3schools.com"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xsi:schemaLocation="http://www.w3schools.com coupling.xsd">
   <redirect redirect_of_root="false" redirect_stdout="false"/>
   <components>
      <component id="1">
         <name>atmo</name>
         <model>ICON</model>
         <simulated>atmosphere</simulated>
         <transient_grid_refs>
            <transient_grid_ref collection_size="2" grid_ref="1" id="1" transient_ref="1"/>
            <transient_grid_ref collection_size="2" grid_ref="1" id="2" transient_ref="2"/>
            <transient_grid_ref collection_size="3" grid_ref="1" id="3" transient_ref="3"/>
            <transient_grid_ref collection_size="4" grid_ref="1" id="4" transient_ref="4"/>
            <transient_grid_ref collection_size="2" grid_ref="1" id="5" transient_ref="5"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="6" transient_ref="6"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="7" transient_ref="7"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="8" transient_ref="8"/>
            <transient_grid_ref collection_size="3" grid_ref="1" id="9" transient_ref="9"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="10" transient_ref="10"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="11" transient_ref="11"/>
         </transient_grid_refs>
      </component>
      <component id="2">
         <name>ocean</name>
         <model>ICON</model>
         <simulated>ocean</simulated>
         <transient_grid_refs>
            <transient_grid_ref collection_size="2" grid_ref="1" id="1" transient_ref="1"/>
            <transient_grid_ref collection_size="2" grid_ref="1" id="2" transient_ref="2"/>
            <transient_grid_ref collection_size="3" grid_ref="1" id="3" transient_ref="3"/>
            <transient_grid_ref collection_size="4" grid_ref="1" id="4" transient_ref="4"/>
            <transient_grid_ref collection_size="2" grid_ref="1" id="5" transient_ref="5"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="6" transient_ref="6"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="7" transient_ref="7"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="8" transient_ref="8"/>
            <transient_grid_ref collection_size="3" grid_ref="1" id="9" transient_ref="9"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="10" transient_ref="10"/>
            <transient_grid_ref collection_size="1" grid_ref="1" id="11" transient_ref="11"/>
         </transient_grid_refs>
      </component>
   </components>
   <transients>
      <transient id="1" transient_standard_name="surface_downward_eastward_stress"/>
      <transient id="2" transient_standard_name="surface_downward_northward_stress"/>
      <transient id="3" transient_standard_name="surface_fresh_water_flux"/>
      <transient id="4" transient_standard_name="total_heat_flux"/>
      <transient id="5" transient_standard_name="atmosphere_sea_ice_bundle"/>
      <transient id="6" transient_standard_name="sea_surface_temperature"/>
      <transient id="7" transient_standard_name="eastward_sea_water_velocity"/>
      <transient id="8" transient_standard_name="northward_sea_water_velocity"/>
      <transient id="9" transient_standard_name="ocean_sea_ice_bundle"/>
      <transient id="10" transient_standard_name="10m_wind_speed"/>
      <transient id="11" transient_standard_name="river_runoff"/>
   </transients>
   <grids>
      <grid alias_name="grid1" id="1"/>
   </grids>
   <dates>
      <start_date>+1800-01-01T00:00:00.000</start_date>
      <end_date>+2100-01-01T00:00:00.000</end_date>
      <calendar>proleptic-gregorian</calendar>
   </dates>
   <timestep_unit>second</timestep_unit>
   <couples>
      <couple>
         <component1 component_id="1"/>
         <component2 component_id="2"/>
         <transient_couple transient_id="1">
            <source component_ref="1" transient_grid_ref="1"/>
            <target transient_grid_ref="1"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation allow_extrapolation="true" extend_source_patch="1" gauss_order="7" method="patch_recovery" polynomial_order="3"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="2">
            <source component_ref="1" transient_grid_ref="2"/>
            <target transient_grid_ref="2"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation allow_extrapolation="true" extend_source_patch="1" gauss_order="7" method="patch_recovery" polynomial_order="3"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="3">
            <source component_ref="1" transient_grid_ref="3"/>
            <target transient_grid_ref="3"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="4">
            <source component_ref="1" transient_grid_ref="4"/>
            <target transient_grid_ref="4"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="5">
            <source component_ref="1" transient_grid_ref="5"/>
            <target transient_grid_ref="5"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="6">
            <source component_ref="2" transient_grid_ref="6"/>
            <target transient_grid_ref="6"/>
            <timestep>
               <source>${dtime_oce}</source>
               <target>${dtime_atm}</target>
               <coupling_period operation="none">${dtime_cpl}</coupling_period>
               <source_timelag>1</source_timelag>
               <target_timelag>${atm_lag}</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="7">
            <source component_ref="2" transient_grid_ref="7"/>
            <target transient_grid_ref="7"/>
            <timestep>
               <source>${dtime_oce}</source>
               <target>${dtime_atm}</target>
               <coupling_period operation="none">${dtime_cpl}</coupling_period>
               <source_timelag>1</source_timelag>
               <target_timelag>${atm_lag}</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="8">
            <source component_ref="2" transient_grid_ref="8"/>
            <target transient_grid_ref="8"/>
            <timestep>
               <source>${dtime_oce}</source>
               <target>${dtime_atm}</target>
               <coupling_period operation="none">${dtime_cpl}</coupling_period>
               <source_timelag>1</source_timelag>
               <target_timelag>${atm_lag}</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="9">
            <source component_ref="2" transient_grid_ref="9"/>
            <target transient_grid_ref="9"/>
            <timestep>
               <source>${dtime_oce}</source>
               <target>${dtime_atm}</target>
               <coupling_period operation="none">${dtime_cpl}</coupling_period>
               <source_timelag>1</source_timelag>
               <target_timelag>${atm_lag}</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation enforced_conservation="false" method="conservative" normalisation="FRACAREA" partial_coverage="true"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="11">
            <source component_ref="1" transient_grid_ref="11"/>
            <target transient_grid_ref="11"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="true" use_target_mask="true">
               <interpolation method="user_file" filename="runoff_a2o.weights.nc"/>
               <interpolation method="fixed_value" user_value="0.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
         <transient_couple transient_id="10">
            <source component_ref="1" transient_grid_ref="10"/>
            <target transient_grid_ref="10"/>
            <timestep>
               <source>${dtime_atm}</source>
               <target>${dtime_oce}</target>
               <coupling_period operation="average">${dtime_cpl}</coupling_period>
               <source_timelag>${atm_lag}</source_timelag>
               <target_timelag>1</target_timelag>
            </timestep>
            <interpolation_requirements use_source_mask="false" use_target_mask="true">
               <interpolation allow_extrapolation="true" extend_source_patch="1" gauss_order="7" method="patch_recovery" polynomial_order="1"/>
               <interpolation method="fixed_value" user_value="-999.0"/>
            </interpolation_requirements>
            <debug_mode at_source_after_interpolation="false" at_source_before_interpolation="false" at_target="false"/>
            <enforce_write_restart>false</enforce_write_restart>
            <enforce_write_weight_file filename="">false</enforce_write_weight_file>
         </transient_couple>
      </couple>
   </couples>
   <created date="19-04-2016 15:39" tool="YAC-CouplingGUI v.1.2.0"/>
</coupling>
EOF

# Mapping parameter:
#  <interpolation method="n-nearest_neighbor" n="1" weighted="ARITHMETIC_AVERAGE"/>
#  <interpolation allow_extrapolation="true" extend_source_patch="1" gauss_order="7" method="patch_recovery" polynomial_order="1"/>
#  <interpolation method="source_to_target_map"/>
#  <interpolation method="fixed_value" user_value="99.0"/>
#  <interpolation method="user_file" filename="runoff_a2o.weights.nc"/>
#  <enforce_write_weight_file filename="runoff_a2o_written.weights.nc">true</enforce_write_weight_file>

# IV.4 xsd and xml files for yac
# ------------------------------

add_required_file ${basedir}/run/${cpl_xmlctrfn}  ./coupling.xml
add_required_file ${basedir}/run/coupling.xsd  ./
add_required_file ${basedir}/run/component.xsd ./

# set grid filenames to undefined to avoid copy mechanism by exec.iconrun (not tested, change later):
#atmo_dyn_grids=""
#ocean_grids=""

