#! /bin/ksh
# ----------------------------------------------------------------------------
#PBS -q xc_norm_h
#PBS -l select=10:ompthreads=4
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -W umask=022
#PBS -o LOG.exp.run_ICON_06_R02B06N07_UPATMO_ifsinit_restarttest.run.run.o
# ----------------------------------------------------------------------------
# ============================================================================
#
# Cray batch script for the ICON model
#
# Real-data test case starting with interpolated IFS input.
# Setup with nest.
#
# Platform: Cray XC30
#           xce.dwd.de
# ---------------------------------------------------------------
# MERGING OF EXPERIMENTS run_ICON_05_R02B06N07_ifsinit_restarttest
# AND run_ICON_14_R2B6N7_oper_IAU_and_restarttest
# WITH MODIFICATIONS FOR AN UPPER-ATMOSPHERE SIMULATION.
#
# Please note: 
# * three runs are performed:
#   1)    simple run: all upper-atmosphere physics are switched on
#   2, 3) restart test: upper-atmosphere radiation has to be switched off, 
#         since it is inherently prone to the accumulation 
#         of round-off errors and will fail restart tests 
#         (at least under the typical compiler settings 
#         for operational use on the CRAY)
# 
# * in case that users would like to use parts of this script
#   as a template for upper-atmosphere simulations, 
#   all upper-atmosphere-related namelist switches are stated explicitly,
#   although the settings are often identical with the default values.
# ---------------------------------------------------------------
#
# 05/2014 : G. Zaengl, F. Prill, DWD
#
# ----------------------------------------------------------------------------
# ============================================================================

set -x

# OpenMP settings
export OMP_SCHEDULE="static"
export OMP_DYNAMIC="false"
# MPI: use DMAPP for off-node data movement
export MPICH_RMA_OVER_DMAPP=1
# disable core dumps
export RLIMIT_CORE=0
export ATP_MAX_CORES=0

# ----------------------------------------------------------------------------
# path definitions
# ----------------------------------------------------------------------------

# root directory for input data
DATAROOT=/lustre2/rwork0/routfor/test/icon/Checksuite_data/

# directory with input grids:
GRIDDIR=${DATAROOT}/GRF_R2B6N8_IFStest

# directory with input data from IFS:
IFSDATADIR=${DATAROOT}/Inidata_R2B6N7_IFStest

# directory with external parameter data:
EXTPARDIR=${DATAROOT}/Extpar_R2B6N8_IFStest

# base directory for ICON sources and binary:
ICONDIR=${PBS_O_WORKDIR}/../

# absolute path to directory with plenty of space:
EDIR="exp15"
EXPDIR=${TMPDIR}/${EDIR}/

# path to model binary, including the executable:
MODEL=${ICONDIR}/bin/icon

# ----------------------------------------------------------------------------
# copy input data: grids, external parameters, model
# ----------------------------------------------------------------------------

# the directory for the experiment will be created, if not already there
if [ ! -d $EXPDIR ]; then
    mkdir -p $EXPDIR
fi
cd ${EXPDIR}

# delete existing restart files
rm -rf *restart_*

# grid files
ln -sf ${GRIDDIR}/iconR2B05_DOM00.nc iconR2B05_DOM00.nc
ln -sf ${GRIDDIR}/iconR2B06_DOM01.nc iconR2B06_DOM01.nc
ln -sf ${GRIDDIR}/iconR2B07_DOM02.nc iconR2B07_DOM02.nc

# external parameter (from ExtPar tool)
ln -sf ${EXTPARDIR}/extpar_R2B06_DOM01.nc extpar_iconR2B06_DOM01.nc
ln -sf ${EXTPARDIR}/extpar_R2B07_DOM02.nc extpar_iconR2B07_DOM02.nc

# interpolated IFS data
ln -sf ${IFSDATADIR}/prepiconR2B06_DOM01_2012062100.nc ifs2icon_R2B06_DOM01.nc
ln -sf ${IFSDATADIR}/prepiconR2B07_DOM02_2012062100.nc ifs2icon_R2B07_DOM02.nc

# files needed for radiation
ln -sf ${ICONDIR}/data/ECHAM6_CldOptProps.nc .
ln -sf ${ICONDIR}/data/rrtmg_lw.nc .

# file(s) needed for upper atmosphere
ln -sf ${ICONDIR}/data/upatmo_gases_chemheat.nc .

# copy binary
cp -p $MODEL icon

# ----------------------------------------------------------------------------
# grid namelist settings
# ----------------------------------------------------------------------------

# the grid parameters
atmo_dyn_grids="iconR2B06_DOM01.nc iconR2B07_DOM02.nc"
atmo_rad_grids="iconR2B05_DOM00.nc"

# reconstruct the grid parameters in namelist form
dynamics_grid_filename=""
for gridfile in ${atmo_dyn_grids}; do
  dynamics_grid_filename="${dynamics_grid_filename} '${gridfile}',"
done
radiation_grid_filename=""
for gridfile in ${atmo_rad_grids}; do
  radiation_grid_filename="${radiation_grid_filename} '${gridfile}',"
done

# global timing
start_date="2012-06-21T00:00:00Z"
stop_date="2012-06-21T12:00:00Z"

# output timing
output_start_date=$stop_date
output_end_date="2012-06-22T00:00:00Z"
output_interval="P01D"

# time step
dtime=180

# tendency update periods for upper-atmosphere physics
dt_imf_dom1=$((${dtime}*2))
dt_imf_dom2=$((${dtime}*1))
dt_rad_dom1=$((${dtime}*3))
dt_rad_dom2=$((${dtime}*2))

# timing of computation of upper-atmosphere physics tendencies
start_date_imf=" "
stop_date_imf="2012-06-21T11:30:00Z"
start_date_rad="2012-06-21T00:30:00Z"
stop_date_rad=" "

# ----------------------------------------------------------------------------
# ============================================================================
#                     LOOP: MODEL IS EXECUTED THREE TIMES
# ============================================================================

# for error handling
SCRIPT_STATUS=0

for iter in 1 2 3 ; do

if [ $iter -eq 1 ] ; then

echo "##############################################################"
echo "###     Simple run with full upper-atmosphere physics      ###"
echo "##############################################################"

lrestart=".FALSE."
num_io_procs=1
num_restart_procs=0
restart_interval="P02D"
checkpoint_interval="P02D"
ofile_prfx="UPATMO"
imode_rad_dom1=1
imode_rad_dom2=1

elif [ $iter -eq 2 ] ; then

echo "##############################################################"
echo "###           Restart test: set up restart point           ###"
echo "##############################################################"

lrestart=".FALSE."
num_io_procs=1
num_restart_procs=2
restart_interval="P02D"
checkpoint_interval="PT11H"
ofile_prfx="UPATMO_REFF"
imode_rad_dom1=0
imode_rad_dom2=2

else

echo "##############################################################"
echo "###                Restart test: resumed run               ###"
echo "##############################################################"

lrestart=".TRUE."
num_io_procs=1
num_restart_procs=2
restart_interval="P02D"
checkpoint_interval="P02D"
ofile_prfx="UPATMO_REST"
imode_rad_dom1=0
imode_rad_dom2=2

fi

# ----------------------------------------------------------------------------
# create ICON master namelist
# ----------------------------------------------------------------------------

cat > icon_master.namelist << EOF

! master_nml: ----------------------------------------------------------------
&master_nml
 lrestart                   =                  ${lrestart}        ! .TRUE.=current experiment is resumed
/

! master_model_nml: repeated for each model ----------------------------------
&master_model_nml
 model_type                  =                          1         ! identifies which component to run (atmosphere,ocean,...)
 model_name                  =                      "ATMO"        ! character string for naming this component.
 model_namelist_filename     =              "NAMELIST_NWP"        ! file name containing the model namelists
 model_min_rank              =                          1         ! start MPI rank for this model
 model_max_rank              =                      65536         ! end MPI rank for this model
 model_inc_rank              =                          1         ! stride of MPI ranks
/

! master_time_control_nml: specification of date and time --------------------
&master_time_control_nml
 experimentStartDate         =               "$start_date"
 experimentStopDate          =                "$stop_date"
 restartTimeIntval           =         "$restart_interval"
 checkpointTimeIntval        =      "$checkpoint_interval"
/

EOF


# ----------------------------------------------------------------------------
# model namelists
# ----------------------------------------------------------------------------
# For a complete list see doc/Namelist_overview.pdf

cat > NAMELIST_NWP << EOF

! parallel_nml: MPI parallelization -------------------------------------------
&parallel_nml
 nproma                      =                          8         ! loop chunk length
 p_test_run                  =                     .FALSE.        ! .TRUE. means verification run for MPI parallelization
 num_io_procs                =            ${num_io_procs}         ! number of I/O processors
 num_restart_procs           =       ${num_restart_procs}         ! number of restart processors
 iorder_sendrecv             =                          3         ! sequence of MPI send/receive calls
/

! run_nml: general switches ---------------------------------------------------
&run_nml
 ltestcase                   =                     .FALSE.        ! idealized testcase runs
 num_lev                     =                   120, 110         ! number of full levels (atm.) for each domain
 lvert_nest                  =                      .TRUE.        ! vertical nesting
 dtime                       =                   ${dtime}         ! timestep in seconds
 ldynamics                   =                      .TRUE.        ! compute adiabatic dynamic tendencies
 ltransport                  =                      .TRUE.        ! compute large-scale tracer transport
 ntracer                     =                          5         ! number of advected tracers
 iforcing                    =                          3         ! forcing of dynamics and transport by parameterized processes
 msg_level                   =                         12         ! controls how much printout is written during runtime
 ltimer                      =                      .TRUE.        ! timer for monitoring the runtime of specific routines
 timers_level                =                         10         ! performance timer granularity
 output                      =                       "nml"        ! main switch for enabling/disabling components of the model output
/

! diffusion_nml: horizontal (numerical) diffusion ----------------------------
&diffusion_nml
 hdiff_order                 =                          5         ! order of nabla operator for diffusion
 itype_vn_diffu              =                          1         ! reconstruction method used for Smagorinsky diffusion
 itype_t_diffu               =                          2         ! discretization of temperature diffusion
 hdiff_efdt_ratio            =                         36.0       ! ratio of e-folding time to time step 
 hdiff_smag_fac              =                          0.025     ! scaling factor for Smagorinsky diffusion
 lhdiff_vn                   =                      .TRUE.        ! diffusion on the horizontal wind field
 lhdiff_temp                 =                      .TRUE.        ! diffusion on the temperature field
/

! dynamics_nml: dynamical core -----------------------------------------------
&dynamics_nml
 iequations                  =                          3         ! type of equations and prognostic variables
 idiv_method                 =                          1         ! method for divergence computation
 divavg_cntrwgt              =                          0.50      ! weight of central cell for divergence averaging
 lcoriolis                   =                      .TRUE.        ! Coriolis force
 ldeepatmo                   =                      .TRUE.        ! .TRUE. -> deep-atmosphere modification is applied 
                                                                  ! to the non-hydrostatic model equations 
                                                                  ! (default is .FALSE. -> standard shallow-atmosphere dynamics)
/

! extpar_nml: external data --------------------------------------------------
&extpar_nml
 extpar_filename             =  "<path>extpar_<gridfile>"         ! filename of external parameter input file
 itopo                       =                          1         ! topography (0:analytical)
 n_iter_smooth_topo          =                          1         ! iterations of topography smoother
 heightdiff_threshold        =                       3000.0       ! height difference between neighb. grid points
/

! initicon_nml: specify read-in of initial state ------------------------------
&initicon_nml
  zpbl1                      =                        500.0       ! bottom height of layer used for gradient computation
  zpbl2                      =                       1000.0       ! top height of layer used for gradient computation
  itype_vert_expol           =                          2         ! type of vertical extrapolation of initial data: 
                                                                  ! 1: linear extrapolation (default) 
                                                                  ! 2: blending with climatology 
                                                                  ! (intended for simulations with the upper-atmosphere configuration)
/

! grid_nml: horizontal grid --------------------------------------------------
&grid_nml
 dynamics_grid_filename      =  ${dynamics_grid_filename}         ! array of the grid filenames for the dycore
 radiation_grid_filename     = ${radiation_grid_filename}         ! array of the grid filenames for the radiation model
 dynamics_parent_grid_id     =                       0, 1         ! array of the indexes of the parent grid filenames
 lredgrid_phys               =               .TRUE.,.TRUE.        ! .true.=radiation is calculated on a reduced grid
 lfeedback                   =                      .FALSE.       ! specifies if feedback to parent grid is performed
/

! gridref_nml: grid refinement and nesting -----------------------------------
&gridref_nml
 grf_intmethod_e             =                          6         ! interpolation method for grid refinement
 grf_scalfbk                 =                          2         ! feedback method for dynamical scalar variables
 grf_tracfbk                 =                          2         ! feedback method for tracer variables
 denom_diffu_v               =                        150.0       ! Denominator for lateral boundary diffusion of velocity
/

! io_nml: general switches for model I/O -------------------------------------
&io_nml
 dt_diag                     =                      21600.0       ! diagnostic integral output interval
 itype_pres_msl              =                          2         ! method for computation of mean sea level pressure
/

! nonhydrostatic_nml: nonhydrostatic model -----------------------------------
&nonhydrostatic_nml
 iadv_rhotheta               =                          2         ! advection method for rho and rhotheta
 ivctype                     =                          2         ! type of vertical coordinate
 itime_scheme                =                          4         ! time integration scheme
 exner_expol                 =                          0.333     ! temporal extrapolation of Exner function
 vwind_offctr                =                          0.2       ! off-centering in vertical wind solver
 damp_height                 =            120000.0, 80000.0       ! height at which Rayleigh damping of vertical wind starts
 rayleigh_coeff              =                   10.0, 10.0       ! Rayleigh damping coefficient
 ndyn_substeps               =                          5         ! number of dynamical core substeps
 lhdiff_rcf                  =                      .TRUE.        ! .TRUE.=compute diffusion only at advection time steps
 divdamp_order               =                         24         ! order of divergence damping
 divdamp_fac                 =                          0.004     ! scaling factor for divergence damping 
 l_open_ubc                  =                     .FALSE.        ! .TRUE.=use open upper boundary condition
 igradp_method               =                          3         ! discretization of horizontal pressure gradient
 l_zdiffu_t                  =                      .TRUE.        ! specifies computation of Smagorinsky temperature diffusion
 thslp_zdiffu                =                          0.02      ! slope threshold (temperature diffusion)
 thhgtd_zdiffu               =                        125.0       ! threshold of height difference (temperature diffusion)
 htop_moist_proc             =                      22500.0       ! max. height for moist physics
 hbot_qvsubstep              =                      19500.0       ! height above which QV is advected with substepping scheme
/

! nwp_phy_nml: switches for the physics schemes ------------------------------
&nwp_phy_nml
 inwp_gscp                   =                          1         ! cloud microphysics and precipitation
 inwp_convection             =                          1         ! convection
 inwp_radiation              =                          1         ! radiation
 inwp_cldcover               =                          1         ! cloud cover scheme for radiation
 inwp_turb                   =                          1         ! vertical diffusion and transfer
 inwp_satad                  =                          1         ! saturation adjustment
 inwp_sso                    =                          1         ! subgrid scale orographic drag
 inwp_gwd                    =                          1         ! non-orographic gravity wave drag
 inwp_surface                =                          1         ! surface scheme
 latm_above_top              =            .FALSE., .FALSE.        ! take into account atmosphere above model top for radiation computation
 efdt_min_raylfric           =                       7200.0       ! minimum e-folding time of Rayleigh friction
 itype_z0                    =                          2         ! type of roughness length data
 icapdcycl                   =                          3         ! type of CAPE correction
 icpl_aero_conv              =                          1         ! simple coupling between autoconversion and Tegen aerosol climatology
 icpl_aero_gscp              =                          1         ! simple coupling between autoconversion and Tegen aerosol climatology
 dt_rad                      =                       1440.        ! time interval of radiation call
 dt_conv                     =                        420.        ! time interval of convection and cloud cover call
 dt_sso                      =                        840.        ! time interval of call of subgrid-scale orographic drag
 dt_gwd                      =                        840.        ! time interval of call of non-orographic gravity-wave drag
 lupatmo_phy                 =                      .TRUE.        ! .TRUE. -> switch on upper-atmosphere physics for NWP forcing
                                                                  ! (default is .FALSE.)
/

! sleve_nml: vertical level specification -------------------------------------
&sleve_nml
 min_lay_thckn               =                         20.0       ! layer thickness of lowermost layer
 max_lay_thckn               =                        700.0       ! max. layer thickness below htop_thcknlimit
 htop_thcknlimit             =                      40000.0       ! level below whick layer thickness is <= max_lay_thckn
 top_height                  =                     150000.0       ! height of model top
 stretch_fac                 =                          1.0       ! stretching factor to vary distribution of model levels
 decay_scale_1               =                       4000.0       ! decay scale of large-scale topography component
 decay_scale_2               =                       2500.0       ! decay scale of small-scale topography component
 decay_exp                   =                          1.2       ! exponent of decay function
 flat_height                 =                      16000.0       ! height above which the coordinate surfaces are flat
/

! radiation_nml: radiation scheme ---------------------------------------------
&radiation_nml
 irad_o3                     =                          7         ! ozone climatology
 irad_aero                   =                          6         ! aerosols
 albedo_type                 =                          2         ! type of surface albedo
 vmr_co2                     =                        390.e-06    ! volume mixing ratio of CO2 
 vmr_ch4                     =                       1800.e-09
 vmr_n2o                     =                        322.0e-09
 vmr_o2                      =                          0.20946
 vmr_cfc11                   =                        240.e-12
 vmr_cfc12                   =                        532.e-12 
/

! transport_nml: tracer transport ---------------------------------------------
&transport_nml
 ivadv_tracer                =              3, 3, 3, 3, 3         ! tracer specific method to compute vertical advection
 itype_hlimit                =           3, 4, 4, 4, 4, 0         ! type of limiter for horizontal transport
 ihadv_tracer                =          52, 2, 2, 2, 2, 0         ! tracer specific method to compute horizontal advection
 llsq_svd                    =                      .TRUE.        ! use SV decomposition for least squares design matrix
 beta_fct                    =                          1.005     ! factor of allowed over-/undershooting in monotonous limiter
/

! turbdiff_nml: turbulent diffusion -------------------------------------------
&turbdiff_nml
 tkhmin                      =                          0.75      ! scaling factor for minimum vertical diffusion coefficient
 tkmmin                      =                          0.75      ! scaling factor for minimum vertical diffusion coefficient
 pat_len                     =                        300.0       ! effective length scale of thermal surface patterns
 c_diff                      =                          0.2       ! length scale factor for vertical diffusion of TKE
 rat_sea                     =                          8.0       ! controls laminar resistance for sea surface
 frcsmot                     =                          0.2       ! vertical smoothing factor for TKE forcing terms
 imode_frcsmot               =                          2         ! restrict vertical smoothing to the tropics
 itype_sher                  =                          3         ! vertical and horizontal shear forcing
 ltkeshs                     =                      .TRUE.        ! coarse grid correction for horizontal shear production term
 a_hshr                      =                          2.0       ! length scale factor for separated horizontal shear mode
 icldm_turb                  =                          1         ! only grid scale condensation possible
/

! lnd_nml: land scheme switches -----------------------------------------------
&lnd_nml
 ntiles                      =                          3         ! number of tiles
 nlev_snow                   =                          3         ! number of snow layers
 lmulti_snow                 =                     .FALSE.        ! .TRUE. for use of multi-layer snow model
 idiag_snowfrac              =                         20         ! type of snow-fraction diagnosis
 lsnowtile                   =                      .TRUE.        ! .TRUE.=consider snow-covered and snow-free separately
 itype_root                  =                          2         ! root density distribution
 itype_heatcond              =                          2         ! type of soil heat conductivity
 itype_lndtbl                =                          3         ! table for associating surface parameters
 lseaice                     =                      .TRUE.        ! .TRUE. for use of sea-ice model
 llake                       =                      .TRUE.        ! .TRUE. for use of lake model
/

! upatmo_nml: upper-atmosphere switches ---------------------------------------
&upatmo_nml
 !----------------------------------------------
 ! deep-atmosphere dynamics: ldeepatmo = .TRUE.
 !----------------------------------------------
 lnontrad                    =                      .TRUE.        ! .FALSE. -> non-traditional deep-atmosphere terms 
                                                                  ! in components of momentum equation are switched off 
                                                                  ! (default is .TRUE.)
 lconstgrav                  =                     .FALSE.        ! .TRUE. -> gravitational acceleration is const. 
                                                                  ! like in case of the shallow atmosphere 
                                                                  ! (default is .FALSE.)
 lcentrifugal                =                     .FALSE.        ! .TRUE. -> explicit centrifugal acceleration is switched on 
                                                                  ! (default is .FALSE.)
 ldeepatmo2phys              =                     .FALSE.        ! .TRUE. -> the input fields to the ECHAM physics parameterizations 
                                                                  ! are modified for the deep atmosphere, if required 
                                                                  ! .FALSE. -> the input fields are computed in accordance 
                                                                  ! with the shallow-atmosphere approximation in any case (default)
 !------------------------------------------------------
 ! upper-atmosphere extrapolation: itype_vert_expol = 2
 !------------------------------------------------------
 expol_start_height          =                      70000.0       ! (m) height above which extrapolation (blending) 
                                                                  ! of initial data starts (default value is 70000 m)
 expol_blending_scale        =                      10000.0       ! (m) blending scale height (default value is 10000 m)
 expol_vn_decay_scale        =                      10000.0       ! (m) scale height for (exponential) decay of extrapolated 
                                                                  ! horizontal wind component (for stability reasons)
                                                                  ! (default value is 10000m)
 expol_temp_infty            =                        400.0       ! (K) climatological temperature of exosphere (z -> infinity)  
 lexpol_sanitycheck          =                      .TRUE.        ! .TRUE. -> apply sanity check to extrapolated fields 
                                                                  ! (default is .FALSE.)
 !----------------------------------------------------------------
 ! upper-atmosphere physics for NWP forcing: lupatmo_phy = .TRUE. 
 !----------------------------------------------------------------
 orbit_type                  =                           1        ! orbit model: 
                                                                  ! * 1: vsop87 -> standard and accurate model (default)
                                                                  ! * 2: kepler -> simple model, appropriate for idealized work      
 solvar_type                 =                           1        ! solar activity: 
                                                                  ! * 1: normal activity (default)
                                                                  ! * 2: low --,,--
                                                                  ! * 3: high --,,--
 solvar_data                 =                           2        ! solar activity data type:
                                                                  ! * 1: G. Rottman data
                                                                  ! * 2: J. Lean data (default)
 solcyc_type                 =                           2        ! solar cycle:
                                                                  ! * 1: standard cycle
                                                                  ! * 2: 27day cycle (default)
 nwp_grp_imf%imode           =                        1, 1        ! mode for physics group:
                                                                  ! (i)on drag, (m)olecular diffusion and (f)rictional heating
                                                                  ! * 0: switched off
                                                                  ! * 1: switched on (default)
                                                                  ! * 2: offline-mode: tendencies are computed,
                                                                  !      but not coupled to the dynamics
 nwp_grp_imf%dt              = ${dt_imf_dom1}, ${dt_imf_dom2}     ! (s) tendency update period (new tendencies are computed every dt)
                                                                  ! (default value is 300 s)
 nwp_grp_imf%t_start         =          "${start_date_imf}"       ! start time for physics group:
                                                                  ! * empty string: simulation start date is start time (default)
                                                                  ! * e.g., "2012-06-21T06:00:00Z": tendency computation would start 
                                                                  !   6 hours after the initial date and time 
                                                                  !   set with ini_datetime_string above. 
                                                                  !   Before the start time, tendencies are set to zero
 nwp_grp_imf%t_end           =           "${stop_date_imf}"       ! end time for physics group: 
                                                                  ! * empty string: simulation end date is end time (default)
                                                                  ! * e.g., "2012-06-21T18:00:00Z": tendency computation 
                                                                  !   would be stopped 18 hours after ini_datetime_string.
                                                                  !   After the end time, tendencies are set to zero
 nwp_grp_imf%start_height    =                        -999.0      ! (m) start height of physics group:
                                                                  ! * negative value: take default start heights
                                                                  !   defined in src/upper_atmosphere/mo_upatmo_impl_const (default)
                                                                  ! * e.g., 50000.0: all processes of the group start to compute
                                                                  !   tendencies above 50 km. Below, all processes are inactive
                                                                  !   and the tendencies are set to zero
 nwp_grp_rad%imode           = ${imode_rad_dom1}, ${imode_rad_dom2}  ! mode for physics group: (RAD)iation and chemical heating.
                                                                  ! The following variables, nwp_grp_rad%..., 
                                                                  ! have the same meaning as with nwp_grp_imf%...
 nwp_grp_rad%dt              = ${dt_rad_dom1}, ${dt_rad_dom2}     ! (default value is 600 s)
 nwp_grp_rad%t_start         =          "${start_date_rad}"
 nwp_grp_rad%t_end           =           "${stop_date_rad}"
 nwp_grp_rad%start_height    =                        -999.0

 nwp_gas_o3%imode            =                          2         ! mode for radiatively active gas: ozone (O3):
                                                                  ! * 0: zero gas concentration
                                                                  ! * 1: constant gas concentration (independent of space and time),
                                                                  !      specified by nwp_gas_o3%vmr
                                                                  ! * 2: external data; meridionally, vertically and monthly varying
                                                                  !      gas concentrations are read from a file 
                                                                  !      with name: nwp_extdat_gases%filename (default)
 nwp_gas_o3%vmr              =                           0.0      ! (m3 m-3) constant volume mixing ratio for ozone
                                                                  ! (for nwp_gas_o3%imode = 1, default value is 0 m3 m-3)
 nwp_gas_o3%fscale           =                           1.0      ! scaling factor the gas concentration is multiplied with 
                                                                  ! in each grid cell (default value is 1)
 nwp_gas_o2%imode            =                           2        ! mode for radiatively active gas: dioxygen (O2).
                                                                  ! The following variables, nwp_gas_o2%...,
                                                                  ! have the same meaning as with nwp_gas_o3%...
 nwp_gas_o2%vmr              =                           0.0
 nwp_gas_o2%fscale           =                           1.0 
 nwp_gas_o%imode             =                           2        ! mode for radiatively active gas: atomic oxygen (O)
 nwp_gas_o%vmr               =                           0.0 
 nwp_gas_o%fscale            =                           1.0 
 nwp_gas_co2%imode           =                           1        ! mode for radiatively active gas: carbon dioxide (CO2)
 nwp_gas_co2%vmr             =                         348.0e-6
 nwp_gas_co2%fscale          =                           1.0
 nwp_gas_no%imode            =                           2        ! mode for radiatively active gas: nitric oxide (NO)
 nwp_gas_no%vmr              =                           0.0
 nwp_gas_no%fscale           =                           1.0 
 nwp_extdat_gases%dt         =                       86400.0      ! (s) update period for time interpolation 
                                                                  ! of gas concentrations from external data (default value is 86400 s)
 nwp_extdat_gases%filename   =   "upatmo_gases_chemheat.nc"       ! name of file containing external data
 nwp_extdat_chemheat%dt      =                       86400.0      ! (s) update period for time interpolation 
                                                                  ! of temperature tendencies of chemical heating
                                                                  ! from external data (default value is 86400 s)
 nwp_extdat_chemheat%filename=   "upatmo_gases_chemheat.nc"       ! name of file containing external data
/

EOF

if [ $iter -eq 1 ] ; then

cat >> NAMELIST_NWP << EOF

! output_nml: specifies an output stream --------------------------------------
&output_nml
 filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
 dom                         =                          2         ! write domain 2
 output_start                =      "${output_start_date}"        ! start date and time of output
 output_end                  =        "${output_end_date}"        ! end date and time of output 
 output_interval             =        "${output_interval}"        ! time interval between two outputs
 steps_per_file              =                          5         ! number of output steps in one output file
 mode                        =                          2         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
 include_last                =                      .TRUE.        ! flag whether to include the last time step
 output_filename             =             "${ofile_prfx}"        ! file name base
 output_grid                 =                      .TRUE.        ! flag whether grid information is added to output.
 !
 ml_varlist                  =   'group:upatmo_tendencies', 'group:upatmo_rad_gases', 'upatmo_mdry', 'upatmo_amd', 'upatmo_cpair', 
'upatmo_grav', 'upatmo_sclrlw', 'upatmo_effrsw'
/

EOF

else

cat >> NAMELIST_NWP << EOF

! output_nml: specifies an output stream --------------------------------------
&output_nml
 filetype                    =                          4         ! output format: 2=GRIB2, 4=NETCDFv2
 dom                         =                          1         ! write domain 1
 output_start                =      "${output_start_date}"        ! start date and time of output
 output_end                  =        "${output_end_date}"        ! end date and time of output 
 output_interval             =        "${output_interval}"        ! time interval between two outputs
 steps_per_file              =                          5         ! number of output steps in one output file
 mode                        =                          2         ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
 include_last                =                      .TRUE.        ! flag whether to include the last time step
 output_filename             =             "${ofile_prfx}"        ! file name base
 output_grid                 =                      .TRUE.        ! flag whether grid information is added to output.
 !
 ml_varlist                  =                   'vn', 'w'
/

EOF

fi

# ----------------------------------------------------------------------------
# run the model!
# ----------------------------------------------------------------------------

# "aprun" command:
# -n xx   : number of MPI tasks
# -N xx   : number of MPI tasks/node
# -d  x   : number of threads/MPI task
# -j 2    : Hyperthreading enabled: 20 physical cores -> 40 "virtual" cores
# -m 6g   : 6G memory/task

aprun  -n 120 -N 12 -j 2 -d 4 -m 3g icon

SCRIPT_STATUS=$?

# in case of an error there is no reason to execute the next iteration
# or what follows after the for-loop
if [ $SCRIPT_STATUS -ne 0 ]; then
  echo " Iteration ${iter} resulted in error"
  echo " SCRIPT_STATUS: "$SCRIPT_STATUS
  exit $SCRIPT_STATUS
fi

# ----------------------------------------------------------------------------
done # iterations for restart testing
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# check model output
# ----------------------------------------------------------------------------

# absolute path to cdo
CDODIR="/e/uhome/dreinert/bin"

SCRIPT_STATUS=0
TEST_STATUS=0

# output filenames
FILE_SIMP="UPATMO_DOM02_ML_0001.nc"
FILE_REFF="UPATMO_REFF_DOM01_ML_0001.nc"
FILE_REST="UPATMO_REST_DOM01_ML_0001.nc"

# check, if files exist
if [ ! -f ${FILE_SIMP} ]; then
  echo " The file ${FILE_SIMP} does not exist"
  SCRIPT_STATUS=$(($SCRIPT_STATUS + 1))
  exit $SCRIPT_STATUS
fi

if [ ! -f ${FILE_REFF} ] || [ ! -f ${FILE_REST} ]; then
  echo " Both or either of the files ${FILE_REFF} and ${FILE_REST} do not exist"
  SCRIPT_STATUS=$(($SCRIPT_STATUS + 1))
  exit $SCRIPT_STATUS
fi

# difference file
DIFF_FILE="DIFF_DOM01_ML_0001.out"

echo "##############################################################"
echo "###       Restart test: check match of output files        ###"
echo "##############################################################"

${CDODIR}/cdo diffn ${FILE_REFF} ${FILE_REST} > ${DIFF_FILE}

SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
  echo " Command: 'cdo diffn' failed"
  echo " SCRIPT_STATUS: "$SCRIPT_STATUS
  exit $SCRIPT_STATUS
fi

if [ -s ${DIFF_FILE} ]; then
  echo " Restart test: file 1 and file 2 differ"
  echo " - file 1: ${FILE_REFF}"
  echo " - file 2: ${FILE_REST}"
  cat ${DIFF_FILE}
  rm  ${DIFF_FILE}
  #
  TEST_STATUS=$(($TEST_STATUS + 1))
fi

if [ $TEST_STATUS == 0 ]; then
  echo " "
  echo "Upper atmosphere test SUCCESSFUL"
  echo " "
else
  echo " "
  echo "Upper atmosphere test FAILED:"
  echo " "
fi

# remove data files
rm -rf ${FILE_SIMP} ${FILE_REFF} ${FILE_REST}

SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
  echo " Deletion of files ${FILE_SIMP}, ${FILE_REFF} and ${FILE_REST} failed"
  echo " SCRIPT_STATUS: "$SCRIPT_STATUS
  exit $SCRIPT_STATUS
fi

# remove restart files
SYMLINK_DOM1="restart_atm_DOM01.nc"
SYMLINK_DOM2="restart_atm_DOM02.nc"
REST_FILE_DOM1="iconR2B06_DOM01_restart_atm_20120621T110000Z.nc"
REST_FILE_DOM2="iconR2B07_DOM02_restart_atm_20120621T110000Z.nc"

rm -rf ${SYMLINK_DOM1} ${SYMLINK_DOM2}

SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
  echo " Deletion of symbolic links ${SYMLINK_DOM1} and ${SYMLINK_DOM2} failed"
  echo " SCRIPT_STATUS: "$SCRIPT_STATUS
  exit $SCRIPT_STATUS
fi

rm -rf ${REST_FILE_DOM1} ${REST_FILE_DOM2}

SCRIPT_STATUS=$?

if [ $SCRIPT_STATUS -ne 0 ]; then
  echo " Deletion of restart files ${REST_FILE_DOM1} and ${REST_FILE_DOM2} failed"
  echo " SCRIPT_STATUS: "$SCRIPT_STATUS
  exit $SCRIPT_STATUS
fi

echo " TEST_STATUS: "$TEST_STATUS
exit $TEST_STATUS
#
#-----------------------------------------------------------------------------
