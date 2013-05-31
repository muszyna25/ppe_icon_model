!>
!! The main program unit of the lake parameterization scheme FLake.
!! It contains all FLake procedures except for the procedures used 
!! to compute fluxes of momentum and of sensible and latent heat over lakes. 
!! (Note that the FLake flux-calculation routines can optionally be used 
!! to determine fluxes over lakes. By default, fluxes over lakes 
!! are computed in the same way as over the ocean).
!! Communication between FLake and the host atmospheric model (ICON) 
!! occurs through two subroutines, "flake_init" and "flake_interface".
!! In "flake_init", FLake prognostic variables are initialized  
!! and some consistency checks are performed.
!! In "flake_interface", a call to the FLake subroutine "flake_driver" is organized, 
!! where FLake variables are advanced one time step forward 
!! (see "flake_driver" for further comments).
!!
!! FLake (Fresh-water Lake) is a bulk lake parameterization scheme 
!! capable of predicting the water temperature 
!! and the characteristics of ice and snow in lakes of various depth 
!! on the time scales from a few hours to many years. 
!! The scheme is based on a parametric representation of the evolving
!! temperature profile, where the structure of the stratified layer between the
!! upper mixed layer and the basin bottom, the lake thermocline, is described
!! using the concept of self-similarity (assumed shape) of the temperature-depth curve.
!! The concept was put forward by Kitaigorodskii and Miropolsky (1970) to
!! describe the vertical temperature structure of the oceanic seasonal
!! thermocline. It has been successfully used in geophysical applications.
!! A bulk approach based on the assumed shape of the temperature-depth curve 
!! is also used to describe the vertical structure of the thermally active 
!! upper layer of bottom sediments and of the ice and snow cover.
!!
!! FLake incorporates the heat budget equations for the four layers
!! in question, viz., snow, ice, water and bottom sediments, developed with
!! due regard for the vertically distributed character of solar radiation
!! heating. The entrainment equation that incorporates the spin-up term
!! is used to compute the depth of a convectively-mixed layer.
!! A relaxation-type equation is used to compute the wind-mixed layer depth in
!! stable and neutral stratification, where a multi-limit formulation for the
!! equilibrium mixed-layer depth accounts for the effects of the earth's rotation,
!! of the surface buoyancy flux and of the static stability in the thermocline.
!! The equations for the mixed-layer depth are developed with due regard for
!! the volumetric character of the solar radiation heating.
!! Simple thermodynamic arguments are invoked to develop
!! the evolution equations for the ice and snow thickness.
!! The heat flux through the water-bottom sediment interface is computed,
!! using a parameterization proposed by Golosov et al. (1998).
!! The heat flux trough the air-water interface
!! (or through the air-ice or air-snow interface)
!! is provided by the driving atmospheric model.
!!
!! Disposable constants and parameters of FLake are estimated, using
!! independent empirical and numerical data. They should not be re-evaluated
!! when the scheme is applied to a particular lake. The only lake-specific
!! parameters are the lake depth, the optical characteristics of lake water,
!! the temperature at the lower boundary of the thermally active layer of bottom
!! sediments, and the depth of that layer.
!!
!! A detailed description of FLake is given in
!!
!! Mironov, D. V., 2008: 
!! Parameterization of lakes in numerical weather prediction. Description of a lake model. 
!! COSMO Technical Report, No. 11, Deutscher Wetterdienst, Offenbach am Main, Germany, 41 pp.
!
!! Mironov, D., E. Heise, E. Kourzeneva, B. Ritter, N. Schneider, and A. Terzhevik, 2010:
!! Implementation of the lake parameterisation scheme FLake
!! into the numerical weather prediction model COSMO.
!! Boreal Env. Res., 15, 218-230.
!!
!! A snow-ice module of FLake is discussed in
!!
!! Mironov, D., B. Ritter, J.-P. Schulz, M. Buchhold, M. Lange, and E. Machulskaya, 2012:
!! Parameterization of sea and lake ice in numerical weather prediction models
!! of the German Weather Service.
!! Tellus A, 64, 17330. doi:10.3402/tellusa.v64i0.17330
!!
!! Further information is available at the FLake web page http://lakemodel.net.
!!
!! In the present configuration, the bottom sediment module of FLake is switched off
!! and the heat flux at the water-bottom sediment interface is set to zero. 
!! Although the snowfall rate is provided by the driving atmospheric model, 
!! snow over lake ice is not considered explicitly.
!! The effect of snow is accounted for implicitly (parametrically) 
!! through changes in the surface albedo with respect to solar radiation.
!!
!!
!! @author Dmitrii Mironov, DWD.
!!
!! @par Revision History
!!
!! COSMO history:
!! Version    Date       Name
!! ---------- ---------- ----
!! 3.18       2006/03/03 Dmitrii Mironov
!!  Initial release
!! 3.21       2006/12/04 Dmitrii Mironov, Ulrich Schaettler
!!  General update of the scheme
!! V4_11        2009/11/30 Ekaterina Machulskaya, Jan-Peter Schulz
!!  Adaptations for multi-layer snow model
!!  Eliminate option for cold start (is put to INT2LM) (US)
!! V4_13        2010/05/11 Michael Gertz
!!  Adaptions to SVN
!! V4_15        2010/11/19 Ulrich Schaettler
!!  The field t_s_lake is removed again after adaptations in the SST Analysis
!! V4_16        2010/12/07 Dmitrii Mironov
!!  Adaptation to the two time level scheme, update of I/O (SUBROUTINE flake_init)
!! V4_21        2011/12/06 Dmitrii Mironov
!!  Changes for better vectorization
!! V4_23        2012/05/10 Burkhardt Rockel (CLM)
!!  Restrict some writing information to "debug output"
!!  Set  t_snow_mult = t_snow only for grid points depth_lk > 0.0,
!!   otherwise program may crash later on in organize_dynamics
!! V4_24        2012/06/22 Dmitrii Mironov
!!  Added consistency checks at the beginning, because values might be inconsistent
!!    due to grib packing
!! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!!  Replaced qx-variables by using them from the tracer module
!! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!!  Introduced MESSy interface
!!
!! Initial ICON release by Dmitrii Mironov, DWD (2013-MM-DD)
!! (adaptation of the COSMO code for ICON)
!! 
!! Modification by <name>, <organization> (YYYY-MM-DD)
!! - <brief description of modification>
!! 
!! 
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

! Lines embraced with "!_tmp>" and "!_tmp<" contain temporary parts of the code.
! Lines embraced/marked with "!_dev>" and "!_dev<" may be replaced
! as improved formulations are developed and tested.
! Lines embraced/marked with "!_cdm>" and "!_cdm<" are DM's comments that may be helpful to a user.
! Lines embraced/marked with "!_dbg>" and "!_dbg<" are used for debugging purposes only.
! Lines starting with "!_nu" are not used.

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

MODULE mo_flake

!===================================================================================================
  
! The following '#ifdef' statements make it possible to use "mo_flake" within both ICON and COSMO.

!_cdm>
! Notice that due to different data structure (organization of calls, etc.) 
! in ICON and COSMO the subroutines "flake_init" and "flake_interface" 
! in the two models are quite different.
! Since these two subroutines are just communication routines between FLake and the host
! atmospheric models, the actual "physical parameterizations" of FLake remain the same
! in ICON and COSMO.
!_cdm<
  
#ifdef __COSMO__       
  USE data_parameters, ONLY :             &
                            & ireals    , &  !< KIND-type parameter for real variables 
                            & iintegers      !< KIND-type parameter for standard integer variables 
#endif                 

#ifdef __ICON__
  USE mo_kind, ONLY:                   &
                   & ireals    => wp , &  !< KIND-type parameter for real variables
                   & iintegers => i4      !< KIND-type parameter for standard integer variables
#endif

!===================================================================================================

!_tmp>
! These physical constants seem to be needed for surface flux calculations only.
! Then, they are not required as in the ICON surface fluxes are computed outside FLake.
!_nu USE data_constants  , ONLY :   &
!_nu 
!_nu ! physical constants and related variables
!_nu ! -------------------------------------------
!_nu     lh_v,         & ! latent heat of vapourization
!_nu     lh_f,         & ! latent heat of fusion
!_nu     cp_d,         & ! specific heat of dry air at constant pressure
!_nu     rdocp           ! r_d / cp_d
!_nu 
!_tmp<

!_tmp>
! Variables/parameters from "data_modelconfig" are COSMO specific.
! They are not required within ICON.
!_nu USE data_modelconfig, ONLY :   &
!_nu 
!_nu ! horizontal and vertical sizes of the fields and related variables
!_nu ! --------------------------------------------------------------------
!_nu     ie         ,  & ! number of grid points in zonal direction
!_nu     je         ,  & ! number of grid points in meridional direction
!_nu     ke         ,  & ! number of grid points in vertical direction
!_nu     ke1        ,  & ! ke + 1
!_nu     ke_soil    ,  & ! number of layers in multi-layer soil model
!_nu     ke_snow    ,  & ! number of layers in multi-layer snow model
!_nu 
!_nu ! start- and end-indices for the computations in the horizontal layers
!_nu ! -----------------------------------------------------------------------
!_nu !    These variables give the start- and the end-indices of the
!_nu !    forecast for the prognostic variables in a horizontal layer.
!_nu !    Note, that the indices for the wind-speeds u and v differ from
!_nu !    the other ones because of the use of the staggered Arakawa-C-grid.
!_nu !
!_nu     istartpar  ,  & ! start index for computations in the parallel program
!_nu     iendpar    ,  & ! end index for computations in the parallel program
!_nu     jstartpar  ,  & ! start index for computations in the parallel program
!_nu     jendpar    ,  & ! end index for computations in the parallel program
!_nu 
!_nu ! variables for the time discretization and related variables
!_nu ! --------------------------------------------------------------
!_nu     dt         ,  & ! long time-step
!_nu     dt2             ! dt*2.
!_tmp<

!_tmp>
! Variables/arrays from "data_fields" are COSMO specific.
! Due to a different organization of calls in ICON,
! some of these variables should be passed to FLake subroutines
! (or defined elsewhere).
!_nu USE data_fields     , ONLY :   &
!_nu 
!_nu ! external parameter fields                                           (unit)
!_nu ! ----------------------------
!_nu     p0         ,    & ! base state pressure                           ( pa  )
!_nu     fc         ,    & ! coriolis-parameter                            ( 1/s )
!_nu     llandmask  ,    & ! landpoint mask                                (  -  )
!_nu 
!_nu ! external parameter fields of the lake model                         (unit)
!_nu ! ----------------------------
!_nu     fr_lake    ,    & ! lake fraction in a grid element [0,1]         (  -  )
!_nu     depth_lk   ,    & ! lake depth                                    (  m  )
!_nu     fetch_lk   ,    & ! wind fetch over lake                          (  m  )
!_nu     dp_bs_lk   ,    & ! depth of the thermally active layer
!_nu                       ! of bottom sediments                           (  m  )
!_nu     t_bs_lk    ,    & ! climatological temperature at the bottom of
!_nu                       ! the thermally active layer of sediments       (  K  )
!_nu     gamso_lk   ,    & ! attenuation coefficient for
!_nu                       ! solar radiation in lake water                 ( 1/m )
!_nu ! prognostic variables                                                (unit)
!_nu ! -----------------------
!_nu     u          ,    & ! zonal wind speed                              ( m/s )
!_nu     v          ,    & ! meridional wind speed                         ( m/s )
!_nu     t          ,    & ! temperature                                   (  k  )
!_nu     pp         ,    & ! deviation from the reference pressure         ( pa  )
!_nu 
!_nu ! four-dimensional fields (x,y,z,t)
!_nu ! ---------------------------------
!_nu     t_so       ,    & ! multi-layer soil temperature                  (  K  )
!_nu 
!_nu ! fields for surface values and soil model variables                  (unit )
!_nu ! -----------------------------------------------------
!_nu     ps         ,    & ! surface pressure                              ( pa  )
!_nu     t_s        ,    & ! temperature of the ground surface             (  K  )
!_nu     t_snow     ,    & ! temperature of the snow-surface               (  k  )
!_nu     t_snow_mult,    & ! temperature of the snow-surface               (  k  )
!_nu     t_ice      ,    & ! temperature at the snow-ice or
!_nu                       ! air-ice interface                             (  K  )
!_nu     h_snow     ,    & ! height of snow                                (  m  )
!_nu     h_ice      ,    & ! ice thickness                                 (  m  )
!_nu     t_g        ,    & ! weighted surface temperature                  (  K  )
!_nu     qv_s       ,    & ! specific water vapor content at the surface   (kg/kg)
!_nu 
!_nu ! fields computed by the turbulence scheme                            (unit )
!_nu ! -----------------------------------------------------
!_nu     tch        ,    & ! turbulent transfer coefficient for heat       ( -- )
!_nu     tcm        ,    & ! turbulent transfer coefficient for momentum   ( -- )
!_nu 
!_nu ! fields computed by the radiation scheme                              (unit )
!_nu ! -----------------------------------------------------
!_nu     sobs       ,    & ! solar radiation at the ground                 ( W/m2)
!_nu     thbs       ,    & ! thermal radiation at the ground               ( W/m2)
!_nu 
!_nu ! prognostic variables of the lake model                              (unit )
!_nu ! -----------------------------------------------------
!_nu     t_mnw_lk   ,    & ! mean temperature of the water column          (  K  )
!_nu     t_wml_lk   ,    & ! mixed-layer temperature                       (  K  )
!_nu     t_bot_lk   ,    & ! temperature at the water-bottom sediment
!_nu                       ! interface                                     (  K  )
!_nu     t_b1_lk    ,    & ! temperature at the bottom of the upper layer
!_nu                       ! of the sediments                              (  K  )
!_nu     c_t_lk     ,    & ! shape factor with respect to the
!_nu                       ! temperature profile in lake thermocline       (  -  )
!_nu     h_ml_lk    ,    & ! thickness of the mixed-layer                  (  m  )
!_nu     h_b1_lk           ! thickness of the upper layer
!_nu                       ! of bottom sediments                           (  m  )
!_nu 
!_tmp<

!===================================================================================================

  USE mo_data_flake, ONLY: &
    ! flake configure
    &  lflk_botsed_use   , & !< .TRUE. indicates that bottom-sediment scheme is used
    &  rflk_depth_bs_ref     !< Reference value of depth of thermally active layer 
                             !< of bottom sediments [m]

  USE mo_data_flake, ONLY:        &
    ! derived types
    &  opticpar_medium          , & !< derived type 
    ! optical characteristics of water, ice and snow 
    &  opticpar_water_ref       , & !< Water (reference)
!_nu    &  opticpar_water_trans     , & !< Transparent water (two-band)
!_nu    &  opticpar_whiteice_ref    , & !< White ice
!_nu    &  opticpar_blueice_ref     , & !< Blue ice
!_nu    &  opticpar_drysnow_ref     , & !< Dry snow
!_nu    &  opticpar_meltingsnow_ref , & !< Melting snow
    &  opticpar_ice_opaque      , & !< Opaque ice
    &  opticpar_snow_opaque         !< Opaque snow

  USE mo_data_flake, ONLY: &
    ! flake parameters
    &  c_cbl_1           , & !< Constant in the CBL entrainment equation [-]
    &  c_cbl_2           , & !< Constant in the CBL entrainment equation [-]
    &  c_sbl_ZM_n        , & !< Constant in the ZM1996 equation for the equilibrium SBL depth [-]
    &  c_sbl_ZM_s        , & !< Constant in the ZM1996 equation for the equilibrium SBL depth [-]
    &  c_sbl_ZM_i        , & !< Constant in the ZM1996 equation for the equilibrium SBL depth [-]
    &  c_relax_h         , & !< Constant in the relaxation equation for the SBL depth [-]
    &  c_relax_C         , & !< Constant in the relaxation equation for the shape factor 
                             !< with respect to the temperature profile in the thermocline [-]
    &  C_T_min           , & !< Minimum value of the shape factor C_T (thermocline) [-]
    &  C_T_max           , & !< Maximum value of the shape factor C_T (thermocline) [-]
    &  Phi_T_pr0_1       , & !< Constant in the expression for the T shape-function derivative [-]
    &  Phi_T_pr0_2       , & !< Constant in the expression for the T shape-function derivative [-]
    &  C_TT_1            , & !< Constant in the expression for C_TT (thermocline) [-]
    &  C_TT_2            , & !< Constant in the expression for C_TT (thermocline) [-]
    &  C_B1              , & !< Shape factor (upper layer of bottom sediments) [-]
    &  C_B2              , & !< Shape factor (lower layer of bottom sediments) [-]
    &  Phi_B1_pr0        , & !< B1 shape-function derivative [-]
    &  C_S_lin           , & !< Shape factor (linear temperature profile in the snow layer) [-]
    &  Phi_S_pr0_lin     , & !< S shape-function derivative (linear profile) [-]
    &  C_I_lin           , & !< Shape factor (linear temperature profile in the ice layer) [-]
    &  Phi_I_pr0_lin     , & !< I shape-function derivative (linear profile) [-]
    &  Phi_I_pr1_lin     , & !< I shape-function derivative (linear profile) [-]
    &  Phi_I_ast_MR      , & !< Constant in the MR2004 expression for I shape factor [-]
    &  C_I_MR            , & !< Constant in the MR2004 expression for I shape factor [-]
    &  H_Ice_max         , & !< Maximum ice thickness in the Mironov/Ritter ice model [m]
    &  h_Snow_min_flk    , & !< Minimum snow thickness [m]
    &  h_Ice_min_flk     , & !< Minimum ice thickness [m]
    &  h_ML_min_flk      , & !< Minimum mixed-layer depth [m]
    &  h_ML_max_flk      , & !< Maximum mixed-layer depth [m]
    &  H_B1_min_flk      , & !< Minimum thickness of the upper layer of bottom sediments [m]
    &  u_star_min_flk    , & !< Minimum value of the surface friction velocity [m s^{-1}]
    &  c_small_flk       , & !< A small number
    &  c_maxearg_flk     , & !< Maximum value of the EXP function argument [-]
    &  fr_lake_min           !< Minimum lake fraction within a host atmospheric model grid box [-]

!_cdm>
! Note that most physical constants are taken from the ICON module 
! "mo_physical_constants" rather than from "mo_data_flake" 
! (the respective lines are marked with "!_nu").
! FLake specific are in fact parameters in empirical approximation formulae
! for some physical properties of different media (e.g. snow).
!_cdm<
  USE mo_data_flake, ONLY: &
    ! flake parameters and physical constants
!_nu    &  tpl_grav          , & !< Acceleration due to gravity [m s^{-2}]
    &  tpl_T_r           , & !< Temperature of maximum density of fresh water [K]
!_nu    &  tpl_T_f           , & !< Fresh water freezing point [K]
    &  tpl_a_T           , & !< Constant in the fresh-water equation of state [K^{-2}]
!_nu    &  tpl_rho_w_r       , & !< Maximum density of fresh water [kg m^{-3}]
!_nu    &  tpl_rho_I         , & !< Density of ice [kg m^{-3}]
    &  tpl_rho_S_min     , & !< Minimum snow density [kg m^{-3}]
    &  tpl_rho_S_max     , & !< Maximum snow density [kg m^{-3}]
    &  tpl_Gamma_rho_S   , & !< Empirical parameter [kg m^{-4}] 
                             !< in the expression for the snow density
!_nu    &  tpl_L_f           , & !< Latent heat of fusion  [J kg^{-1}]
!_nu    &  tpl_c_w           , & !< Specific heat of water [J kg^{-1} K^{-1}]
!_nu    &  tpl_c_I           , & !< Specific heat of ice   [J kg^{-1} K^{-1}]
!_nu    &  tpl_c_S           , & !< Specific heat of snow  [J kg^{-1} K^{-1}]
    &  tpl_kappa_w       , & !< Molecular heat conductivity of water        [J m^{-1} s^{-1} K^{-1}]
!_nu    &  tpl_kappa_I       , & !< Molecular heat conductivity of ice          [J m^{-1} s^{-1} K^{-1}]
    &  tpl_kappa_S_min   , & !< Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
    &  tpl_kappa_S_max   , & !< Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
    &  tpl_Gamma_kappa_S     !< Empirical parameter in expression for the 
                             !< snow heat conductivity [J m^{-2} s^{-1} K^{-1}]

!_cdm>
! Most physical constants taken from the ICON module "mo_physical_constants" 
! rather than from the FLake module "mo_data_flake". 
! Note that values of some constants are slightly different 
! in "mo_physical_constants" and "mo_data_flake".
! FLake: tpl_grav    = 9.81,    ICON: grav  = 9.80665.
! FLake: tpl_rho_I   = 9.1E+02, ICON: rhoi = 917.0.
! FLake: tpl_L_f     = 3.3E+05, ICON: alf  = 3.337E+05.
! FLake: tpl_c_w     = 4.2E+03, ICON: clw  = 4192.664112.
! FLake: tpl_c_I     = 2.1E+03, ICON: ci   = 2106.0.
! FLake: tpl_c_S     = 2.1E+03, ICON: ci   = 2106.0.
! FLake: tpl_kappa_I = 2.29,    ICON: ki   = 2.1656.
!_cdm<
  USE mo_physical_constants, ONLY: &
    &  tpl_grav    => grav       , & !< Acceleration due to gravity [m s^{-2}]
    &  tpl_T_f     => tmelt      , & !< Fresh-water freezing point [K]
    &  tpl_rho_w_r => rhoh2o     , & !< Maximum density of fresh water [kg m^{-3}]
    &  tpl_rho_I   => rhoi       , & !< Density of ice [kg/m^3]
    &  tpl_L_f     => alf        , & !< Latent heat of fusion  [J kg^{-1}]
    &  tpl_c_w     => clw        , & !< Specific heat of water [J kg^{-1} K^{-1}]
    &  tpl_c_I     => ci         , & !< Specific heat of ice   [J kg^{-1} K^{-1}]
    &  tpl_c_S     => ci         , & !< Specific heat of snow  [J kg^{-1} K^{-1}]
    &  tpl_kappa_I => ki             !< Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]

  USE mo_exception, ONLY:       &
                    &  finish , &  !< external procedure, finishes model run and reports the reason
                    &  message     !< external procedure, sends a message (error, warning, etc.)

!===================================================================================================

!_cdm>
! COSMO stuff.
! Currently unused within ICON?
!_cdm<

!_nu USE data_parallel   , ONLY :   &
!_nu   my_cart_id       ! rank of this sub-domain in the Cartesian communicator

!===================================================================================================

!_tmp>
! COSMO specific stuff.
! Not needed within ICON.
!_nu USE data_runcontrol , ONLY :   &
!_nu 
!_nu ! controlling the physics
!_nu ! --------------------------
!_nu     lmulti_layer   , & ! run multi-layer soil model
!_nu     lmulti_snow    , & ! run multi-layer snow model
!_nu     lseaice        , & ! forecast with sea ice model
!_nu 
!_nu ! start and end of the forecast
!_nu ! --------------------------------
!_nu     ntstep       ,  & ! actual time step
!_nu     nold         ,  & ! corresponds to ntstep - 1
!_nu     nnow         ,  & ! corresponds to ntstep
!_nu     nnew         ,  & ! corresponds to ntstep + 1
!_nu 
!_nu ! controlling the dynamics
!_nu ! ---------------------------
!_nu     l2tls        ,  & ! time integration by two time level RK-scheme (.TRUE.)
!_nu                       ! or by three time level KW-scheme (.FALSE.)
!_nu 
!_nu ! 12. controlling verbosity of debug output
!_nu ! -----------------------------------------
!_nu     idbg_level,   & ! to control the verbosity of debug output
!_nu     ldebug_soi,   & ! if .TRUE., debug output for soil and surface
!_nu     lprintdeb_all   ! .TRUE.:  all tasks print debug output
!_nu                     ! .FALSE.: only task 0 prints debug output
!_nu
!_tmp<

!===================================================================================================

!_tmp>
! Within ICON, surface fluxes are computed elsewhere and passed to FLake routines.
! Then, SfcFlx_rhoair is not required.
!_nu USE src_flake_sfcflx, ONLY:    &
!_nu   SfcFlx_rhoair                   ! Function, returns the air density
!_tmp<

!===================================================================================================
  
!_cdm>
! Stuff from COSMO.
! Currently unused within ICON.
!_cdm<

!_nu USE environment,      ONLY: model_abort
  
!===================================================================================================

!_cdm>
! MESSY stuff from COSMO.
! Currently unused within ICON.
!_cdm<

!_nu #ifndef MESSY
!_nu   USE src_tracer,       ONLY: trcr_get, trcr_errorstr
!_nu #else
!_nu   USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!_nu   USE messy_main_tracer_bi,     ONLY: tracer_halt
!_nu   USE messy_main_tracer,        ONLY: get_tracer
!_nu #endif

!===================================================================================================

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! The variables declared below are accessible to all program units of the MODULE "mo_flake".
  ! These are basically the quantities computed by FLake.
  ! All variables declared below have a suffix "flk".

  ! FLake variables of type REAL

  ! Temperatures at the previous time step ("p"), and the updated temperatures ("n")
  REAL (KIND = ireals) ::            &
    &  T_mnw_p_flk,  T_mnw_n_flk   , & !< Mean temperature of the water column [K]
    &  T_snow_p_flk, T_snow_n_flk  , & !< Temperature at the air-snow interface [K]
    &  T_ice_p_flk,  T_ice_n_flk   , & !< Temperature at the snow-ice or air-ice interface [K]
    &  T_wML_p_flk,  T_wML_n_flk   , & !< Mixed-layer temperature [K]
    &  T_bot_p_flk,  T_bot_n_flk   , & !< Temperature at the water-bottom sediment interface [K]
    &  T_B1_p_flk,   T_B1_n_flk        !< Temperature at the bottom of the upper layer     
                                       !< of the sediments [K]

  ! Thickness of various layers at the previous time step ("p") and the updated values ("n")
  REAL (KIND = ireals) ::            &
    &  h_snow_p_flk, h_snow_n_flk  , & !< Snow thickness [m]
    &  h_ice_p_flk,  h_ice_n_flk   , & !< Ice thickness [m]
    &  h_ML_p_flk,   h_ML_n_flk    , & !< Thickness of the mixed-layer [m]
    &  H_B1_p_flk,   H_B1_n_flk        !< Thickness of the upper layer of bottom sediments [m]

  ! The shape factor(s) at the previous time step ("p") and the updated value(s) ("n")
  REAL (KIND = ireals) ::            &
    &  C_T_p_flk, C_T_n_flk        , & !< Shape factor (thermocline) [-]
    &  C_TT_flk                    , & !< Dimensionless parameter (thermocline) [-]
    &  C_Q_flk                     , & !< Shape factor with respect to heat flux (thermocline) [-]
    &  C_I_flk                     , & !< Shape factor (ice) [-]
    &  C_S_flk                         !< Shape factor (snow) [-]

  ! Derivatives of the shape functions
  REAL (KIND = ireals) ::            &
    &  Phi_T_pr0_flk               , & !< d\Phi_T(0)/d\zeta   (thermocline) [-]
    &  Phi_I_pr0_flk               , & !< d\Phi_I(0)/d\zeta_I (ice) [-]
    &  Phi_I_pr1_flk               , & !< d\Phi_I(1)/d\zeta_I (ice) [-]
    &  Phi_S_pr0_flk                   !< d\Phi_S(0)/d\zeta_S (snow) [-]

  ! Heat and radiation fluxes
  REAL (KIND = ireals) ::            &
    &  Q_snow_flk                  , & !< Heat flux through the air-snow interface [W m^{-2}]
    &  Q_ice_flk                   , & !< Heat flux through the snow-ice or air-ice 
                                       !< interface [W m^{-2}]
    &  Q_w_flk                     , & !< Heat flux through the ice-water or air-water 
                                       !< interface [W m^{-2}]
    &  Q_bot_flk                   , & !< Heat flux through the water-bottom sediment 
                                       !< interface [W m^{-2}]
    &  I_atm_flk                   , & !< Radiation flux at the lower boundary of 
                                       !< the atmosphere [W m^{-2}], i.e. the incident radiation flux 
                                       !< with no regard for the surface albedo.
    &  I_snow_flk                  , & !< Radiation flux through the air-snow interface [W m^{-2}]
    &  I_ice_flk                   , & !< Radiation flux through the snow-ice or air-ice 
                                       !< interface [W m^{-2}]
    &  I_w_flk                     , & !< Radiation flux through the ice-water or air-water 
                                       !< interface [W m^{-2}]
    &  I_h_flk                     , & !< Radiation flux through the mixed-layer-thermocline 
                                       !< interface [W m^{-2}]
    &  I_bot_flk                   , & !< Radiation flux through the water-bottom sediment 
                                       !< interface [W m^{-2}]
    &  I_intm_0_h_flk              , & !< Mean radiation flux over the mixed layer [W m^{-2}]
    &  I_intm_h_D_flk              , & !< Mean radiation flux over the thermocline [W m^{-2}]
    &  Q_star_flk                      !< A generalized heat flux scale [W m^{-2}]

  ! Velocity scales
  REAL (KIND = ireals) ::            &
    &  u_star_w_flk                , & !< Friction velocity in the surface layer 
                                       !< of lake water [m s^{-1}]
    &  w_star_sfc_flk                  !< Convective velocity scale based on
                                       !< a generalized heat flux scale [m s^{-1}]

  ! The rate of snow accumulation
  REAL (KIND = ireals) ::            &
    &  dMsnowdt_flk                    !< The rate of snow accumulation [kg m^{-2} s^{-1}]

!===================================================================================================

  PUBLIC ::                       &
         &  flake_init          , & ! procedure (FLake initialization)
         &  flake_interface         ! procedure (interface to FLake time stepping)

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

CONTAINS

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

!===================================================================================================
!  Initialize the lake parameterization scheme FLake
!---------------------------------------------------------------------------------------------------

  !> 
  !! This routine initializes prognostic variables of the lake model FLake.
  !! At present, the bottom-sediment module of FLake is switched off so that 
  !! the external parameters of the bottom-sediment module,
  !! namely, the depth of the thermally active layer of the sediment
  !! and the climatological temperature at the bottom of that layer,
  !! are not required.
  !! These parameters are set to their reference values and are not part of the input 
  !! (not read from external-parameter files).
  !! The same applies to the attenuation coefficient with respect to solar 
  !! radiation and to the wind fetch. Default values of these variables are used.
  !! Three FLake prognostic variables, namely, 
  !! the snow thickness (always zero as snow over lake ice is not considered 
  !! explicitly), thickness of the upper layer of bottom sediment and 
  !! the temperature at the outer edge of that layer 
  !! (these are equal to their reference values as the bottom-sediment module is 
  !! switched off), are not included into the ICON IO list.
  !! These variable are handled internally.
  !!
  !!
  !! @par Revision History
  !! Initial ICON release by Dmitrii Mironov, DWD (2013-MM-DD)
  !! (adaptation of the COSMO code for ICON)
  !! 
  !! Modification by <name>, <organization> (YYYY-MM-DD)
  !! - <brief description of modification>
  !!

  SUBROUTINE flake_init (                                       &
                     &  nflkgb,                                 &
                     &  fr_lake, depth_lk,                      &
                     &  fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk,  &
                     &  t_snow_p, h_snow_p,                     & 
                     &  t_ice_p, h_ice_p,                       & 
                     &  t_mnw_lk_p, t_wml_lk_p, t_bot_lk_p,     &
                     &  c_t_lk_p, h_ml_lk_p,                    & 
                     &  t_b1_lk_p, h_b1_lk_p,                   &           
                     &  t_scf_lk_p,                             &
                     &  t_snow_n, h_snow_n,                     & 
                     &  t_ice_n, h_ice_n,                       & 
                     &  t_mnw_lk_n, t_wml_lk_n, t_bot_lk_n,     &
                     &  c_t_lk_n, h_ml_lk_n,                    & 
                     &  t_b1_lk_n, h_b1_lk_n,                   &           
                     &  t_scf_lk_n                              &
                     &  )

    IMPLICIT NONE

    ! Procedure arguments

    ! Array dimension(s)

    INTEGER, INTENT(IN) ::         &
                        &  nflkgb    !< Array (vector) dimension
                                     !< (equal to the number of grid boxes within a block 
                                     !< where lakes are present)

    ! FLake external parameters

    REAL(KIND = ireals), DIMENSION(:), INTENT(IN) ::  &
                        &  fr_lake          , & !< lake fraction in a grid box [-]
                        &  depth_lk             !< lake depth [m]

    REAL(KIND = ireals), DIMENSION(:), INTENT(INOUT) ::  &
                        &  fetch_lk         , & !< wind fetch over lake [m]
                        &  dp_bs_lk         , & !< depth of the thermally active layer
                                                !< of bottom sediments [m]
                        &  t_bs_lk          , & !< climatological temperature at the bottom of
                                                !< the thermally active layer of sediments [K]
                        &  gamso_lk             !< attenuation coefficient of the lake water 
                                                !< with respect to solar radiation [m^{-1}]

    ! FLake prognostic variables
    ! (at the previous time step - "p", and the updated values - "n")

    REAL(KIND = ireals), DIMENSION(:), INTENT(INOUT) ::  &
                        &  t_snow_p         , & !< temperature of snow upper surface at previous time step [K]
                        &  h_snow_p         , & !< snow thickness at previous time step [m]
                        &  t_ice_p          , & !< temperature of ice upper surface at previous time step [K]
                        &  h_ice_p          , & !< ice thickness at previous time step [m]
                        &  t_mnw_lk_p       , & !< mean temperature of the water column at previous time step [K]
                        &  t_wml_lk_p       , & !< mixed-layer temperature at previous time step [K] 
                        &  t_bot_lk_p       , & !< temperature at the water-bottom sediment interface 
                                                !< at previous time step [K] 
                        &  c_t_lk_p         , & !< shape factor with respect to the temperature profile 
                                                !< in lake thermocline at previous time step [-] 
                        &  h_ml_lk_p        , & !< thickness of the mixed-layer at previous time step [m] 
                        &  t_b1_lk_p        , & !< temperature at the bottom of the upper layer
                                                !< of the sediments at previous time step [K]  
                        &  h_b1_lk_p        , & !< thickness of the upper layer of bottom sediments 
                                                !< at previous time step [m] 
                        &  t_scf_lk_p           !< lake surface temperature at previous time step [K]
                                                !< (i.e. the temperature at the air-water, air-ice 
                                                !< or air-snow interface) 
                       
    REAL(KIND = ireals), DIMENSION(:), INTENT(OUT) ::    &
                        &  t_snow_n         , & !< temperature of snow upper surface at new time step [K]
                        &  h_snow_n         , & !< snow thickness at new time step [m]
                        &  t_ice_n          , & !< temperature of ice upper surface at new time step [K]
                        &  h_ice_n          , & !< ice thickness at new time step [m]
                        &  t_mnw_lk_n       , & !< mean temperature of the water column at new time step [K]
                        &  t_wml_lk_n       , & !< mixed-layer temperature at new time step [K] 
                        &  t_bot_lk_n       , & !< temperature at the water-bottom sediment interface 
                                                !< at new time step [K] 
                        &  c_t_lk_n         , & !< shape factor with respect to the temperature profile 
                                                !< in lake thermocline at new time step [-] 
                        &  h_ml_lk_n        , & !< thickness of the mixed-layer at new time step [m] 
                        &  t_b1_lk_n        , & !< temperature at the bottom of the upper layer
                                                !< of the sediments at new time step [K]  
                        &  h_b1_lk_n        , & !< thickness of the upper layer of bottom sediments 
                                                !< at new time step [m] 
                        &  t_scf_lk_n           !< lake surface temperature at new time step [K]
                                                !< (i.e. the temperature at the air-water, air-ice 
                                                !< or air-snow interface) 

    ! Local variables

    INTEGER ::       &
            &  iflk  !< DO loop index

    CHARACTER(len=256) ::            &
                       &  nameerr  , &  !< name of procedure where an error occurs
                       &  texterr       !< error/warning message text

    LOGICAL ::             &
            &  lcallabort  !< logical switch, set .TRUE. if errors are encountered
                           !< (used to call model abort outside a DO loop)




!===================================================================================================
!_tmp>
! COSMO stuff, unused.
!_nu !  Local variables of type INTEGER
!_nu INTEGER (KIND = iintegers) :: &
!_nu   i, j, ks, nt              , & ! Loop indices
!_nu   istarts                   , & ! Start index for x-direction
!_nu   iends                     , & ! End   index for x-direction
!_nu   jstarts                   , & ! Start index for y-direction
!_nu   jends                     , & ! End   index for y-direction
!_nu   izdebug                   , & ! for debug output
!_nu   nztlev                    , & ! number of time-levels for prognostic variables
!_nu   nx                            ! Time-level for initialization (nnow or nnew)
!_tmp<
!===================================================================================================

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Logical switch, default value is .FALSE.
    lcallabort = .FALSE.

    ! Loop over grid boxes where lakes are (should be) present
    CheckFLakeExtPar: DO iflk=1, nflkgb
      ! Check lake-fraction and lake-depth fields
      IF( (fr_lake(iflk) < fr_lake_min) .OR. (depth_lk(iflk) < 0._ireals) ) THEN
        ! Lake fraction less than a minimum threshold value or negative lake depth is found
        ! Set logical switch
        lcallabort = .TRUE.
        ! Exit DO loop to call model abort
        EXIT CheckFLakeExtPar 
      END IF
    END DO CheckFLakeExtPar 

    ! Call model abort if errors are encountered
    IF( lcallabort ) THEN
      ! Send an error message
      WRITE(nameerr,*) "MODULE mo_flake, SUBROUTINE flake_init"
      WRITE(texterr,*) "Lake fraction ", fr_lake(iflk),                          &
                    &  " is less than a minimum threshold value ", fr_lake_min,  &
                    &  " or negative lake depth ", depth_lk(iflk),                &
                    &  " Call model abort."
      CALL message(TRIM(nameerr), TRIM(texterr))

      ! Call model abort
      WRITE(nameerr,*) "mo_flake:flake_init"
      WRITE(texterr,*) "error in lake fraction or lake depth"
      CALL finish(TRIM(nameerr), TRIM(texterr))
    END IF

    ! As different from lake fraction and lake depth 
    ! that are part from ICON (COSMO) IO (read from GRIB, NetCDF, etc., files), 
    ! other FLake external-parameter fields are handled internally by the host atmospheric model.
    ! Set those FLake external parameters to reference values.
    SetFLakeExtPar: DO iflk=1, nflkgb
      fetch_lk(iflk) = 1.0E+04_ireals        ! Use a constant long fetch
      dp_bs_lk(iflk) = rflk_depth_bs_ref     ! Reference value
      t_bs_lk (iflk) = tpl_T_r               ! Reference value
      gamso_lk(iflk) = opticpar_water_ref%extincoef_optic(1)
                                             ! Reference value
    END DO SetFLakeExtPar 

    ! Since a loss of accuracy may occur during IO
    ! (e.g. due to GRIB encoding/decoding),
    ! FLake prognostic variables must be checked for consistency.
    ! That is, the relations between
    ! "t_wml_lk", "t_mnw_lk", "t_bot_lk", "h_ml_lk", "c_t_lk", "h_ice" and "t_ice",
    ! suggested by the assumed shape (self-similarity) of the temperature-depth curve,
    ! must hold or should be enforced.

    GridBoxesWithLakes: DO iflk=1, nflkgb

      ! FLake variables at previous time step

      ! Limit the mixed-layer depth and the temperature-profile shape factor
      h_ml_lk_p(iflk) = MIN(depth_lk(iflk), MAX(h_ml_lk_p(iflk), 0._ireals))
      c_t_lk_p(iflk)  = MIN(C_T_max, MAX(c_t_lk_p(iflk), C_T_min))
      IF(h_ice_p(iflk) >= h_Ice_min_flk) THEN
        ! Ice-covered lakes
        h_ice_p(iflk)    = MIN(h_ice_p(iflk), H_Ice_max)
        t_ice_p(iflk)    = MIN(t_ice_p(iflk), tpl_T_f)
        t_wml_lk_p(iflk) = tpl_T_f
        t_mnw_lk_p(iflk) = MAX(tpl_T_f, MIN(t_mnw_lk_p(iflk), tpl_T_r))
        t_bot_lk_p(iflk) = MAX(tpl_T_f, MIN(t_bot_lk_p(iflk), tpl_T_r))
        ! Mixing down to the lake bottom or "inverse" stratification,
        ! reset mixed-layer depth, shape factor, mean temperature, and bottom temperature.
        IF( (h_ml_lk_p(iflk) >= depth_lk(iflk)-h_ML_min_flk) .OR.  &
            (t_bot_lk_p(iflk) <= t_wml_lk_p(iflk)) ) THEN
          h_ml_lk_p(iflk)  = 0._ireals
          c_t_lk_p(iflk)   = C_T_min
          t_mnw_lk_p(iflk) = t_wml_lk_p(iflk)
          t_bot_lk_p(iflk) = t_wml_lk_p(iflk)
        END IF
      ELSE
        ! Ice-free lakes 
        h_ice_p(iflk)    = 0.0_ireals
        t_ice_p(iflk)    = tpl_T_f
        t_wml_lk_p(iflk) = MAX(tpl_T_f, t_wml_lk_p(iflk))
        t_mnw_lk_p(iflk) = MAX(tpl_T_f, t_mnw_lk_p(iflk))
        t_bot_lk_p(iflk) = MAX(tpl_T_f, t_bot_lk_p(iflk))
        ! Avoid temperature of maximum density crossover
        IF( ((t_wml_lk_p(iflk) <= tpl_T_r) .AND. (t_bot_lk_p(iflk) > tpl_T_r)) .OR.  &
             ((t_wml_lk_p(iflk) >= tpl_T_r) .AND. (t_bot_lk_p(iflk) < tpl_T_r)) )    &
          t_bot_lk_p(iflk) = tpl_T_r
        ! Mixing down to the lake bottom or "inverse" stratification,
        ! reset mixed-layer depth, shape factor, mixed-layer temperature, and bottom temperature.
        IF( (h_ml_lk_p(iflk) >= depth_lk(iflk)-h_ML_min_flk) .OR.                               &
             ((t_wml_lk_p(iflk) <= tpl_T_r) .AND. (t_bot_lk_p(iflk) < t_wml_lk_p(iflk))) .OR.  &
             ((t_wml_lk_p(iflk) >= tpl_T_r) .AND. (t_bot_lk_p(iflk) > t_wml_lk_p(iflk))) ) THEN
          h_ml_lk_p(iflk)  = depth_lk(iflk)
          c_t_lk_p(iflk)   = C_T_min
          t_wml_lk_p(iflk) = t_mnw_lk_p(iflk)
          t_bot_lk_p(iflk) = t_mnw_lk_p(iflk)
        END IF
      END IF
      ! Snow over lake ice is not considered explicitly,
      ! set "h_snow" to zero and "t_snow" to the ice surface temperature.
      h_snow_p(iflk)   = 0.0_ireals
      t_snow_p(iflk)   = t_ice_p(iflk)
      ! Bottom-sediment module is switched off 
      ! and the respective FLake variables are handled internally, 
      ! i.e. they are not part of ICON (COSMO) IO. 
      ! Set the bottom-sediment FLake variables at reference values 
      t_b1_lk_p(iflk) = tpl_T_r          
      h_b1_lk_p(iflk) = rflk_depth_bs_ref
      ! Set the lake surface temperature
      IF (h_snow_p(iflk) >= h_Snow_min_flk) THEN
        ! Snow exists, use the snow surface temperature
        t_scf_lk_p(iflk) = t_snow_p(iflk)
      ELSE IF (h_ice_p(iflk) >= h_Ice_min_flk) THEN
        ! Ice exists but there is no snow, use the ice surface temperature
        t_scf_lk_p(iflk) = t_ice_p(iflk)
      ELSE
        ! No ice-snow cover, use the water surface temperature ( = mixed-layer temperature)
        t_scf_lk_p(iflk) = t_wml_lk_p(iflk)
      END IF

      ! FLake variables at new time step
      
      h_snow_n  (iflk) = h_snow_p  (iflk)
      t_snow_n  (iflk) = t_snow_p  (iflk)
      h_ice_n   (iflk) = h_ice_p   (iflk)
      t_ice_n   (iflk) = t_ice_p   (iflk)
      t_wml_lk_n(iflk) = t_wml_lk_p(iflk)
      t_mnw_lk_n(iflk) = t_mnw_lk_p(iflk)
      t_bot_lk_n(iflk) = t_bot_lk_p(iflk)
      c_t_lk_n  (iflk) = c_t_lk_p  (iflk)
      h_ml_lk_n (iflk) = h_ml_lk_p (iflk)
      t_b1_lk_n (iflk) = t_b1_lk_p (iflk) 
      h_b1_lk_n (iflk) = h_b1_lk_p (iflk) 
      t_scf_lk_n(iflk) = t_scf_lk_p(iflk)   

    END DO GridBoxesWithLakes 


!_tmp>
! COSMO stuff, unused.
! HOWEVER! Make sure that arrays for multi-layer snow and soil models
! and the like should/should not be initialized with the respective 
! FLake variables.
!_nu
!_nu ! The outermost loop is over the time levels
!_nu DO nt = 1, nztlev
!_nu DO j = jstarts, jends
!_nu !CDIR NOMOVE
!_nu DO i = istarts, iends
!_nu   Lake_points: IF(depth_lk(i,j) > 0.0_ireals) THEN
!_nu     ! At lake points, set "t_s", "t_g" and "t_so" to FLake values
!_nu     t_s(i,j,nt) = t_s(i,j,nx)
!_nu     t_g(i,j,nt) = t_s(i,j,nx)
!_nu     IF (lmulti_layer) THEN
!_nu       ! Save lake surface temperature in "t_so(:,:,0,:)"
!_nu       ! if the multi-layer soil model is used
!_nu       t_so(i,j,0,nt) = t_s(i,j,nx)
!_nu     ENDIF
!_nu     ! Set FLake variables at all time levels
!_nu     ! ("nnew, "nnow", and for three time level scheme "nold")
!_nu     ! to their values at "nx" (nx=nnew for two time level scheme, see above).
!_nu     t_wml_lk(i,j,nt) = t_wml_lk(i,j,nx)
!_nu     t_mnw_lk(i,j,nt) = t_mnw_lk(i,j,nx)
!_nu     t_bot_lk(i,j,nt) = t_bot_lk(i,j,nx)
!_nu     c_t_lk  (i,j,nt) = c_t_lk  (i,j,nx)
!_nu     h_ml_lk (i,j,nt) = h_ml_lk (i,j,nx)
!_nu     h_ice   (i,j,nt) = h_ice   (i,j,nx)
!_nu     t_ice   (i,j,nt) = t_ice   (i,j,nx)
!_nu     ! As snow over lake ice is not considered explicitly,
!_nu     ! "h_snow" is zero
!_nu     h_snow  (i,j,nt) = h_snow  (i,j,nx)
!_nu     ! and "t_snow" is equal to the ice surface temperature
!_nu     t_snow  (i,j,nt) = t_snow  (i,j,nx)
!_nu   ELSE Lake_points
!_nu     ! Set FLake variables at non-lake (sea and land) points at their reference
!_nu     ! values at all time levels. Then, all FLake variables at land and sea points
!_nu     ! remain equal to their initial values during the entire COSMO run
!_nu     ! independent of permutation of time indices.
!_nu     t_wml_lk(i,j,nt) = tpl_T_r
!_nu     t_mnw_lk(i,j,nt) = tpl_T_r
!_nu     t_bot_lk(i,j,nt) = tpl_T_r
!_nu     c_t_lk  (i,j,nt) = C_T_min
!_nu     h_ml_lk (i,j,nt) = 0.0_ireals
!_nu     ! If sea ice scheme is not used,
!_nu     ! "t_ice" and "h_ice" are only used if "llake" is set true.
!_nu     ! Initial values of "t_ice" and "h_ice" are set through the surface analysis.
!_nu     IF(.NOT.lseaice) THEN
!_nu       h_ice (i,j,nt) = h_ice (i,j,nx)
!_nu       t_ice (i,j,nt) = t_ice (i,j,nx)
!_nu     ENDIF
!_nu   ENDIF Lake_points
!_nu ! Variables are not part of the input (not read from GRIB files),
!_nu ! set them to their reference values at all points.
!_nu   t_b1_lk(i,j,nt) = tpl_T_r          ! Reference value, bottom-sediment module
!_nu                                      !     is switched off
!_nu   h_b1_lk(i,j,nt) = rflk_depth_bs_ref! Reference value, bottom-sediment module
!_nu                                      !     is switched off
!_nu END DO
!_nu END DO
!_nu END DO
!_nu 
!_nu ! If multi-layer snow model is used, set "t_snow_mult" at lake points
!_nu ! to the ice surface temperature at all vertical levels.
!_nu IF (lmulti_snow) THEN
!_nu   DO nt = 1, nztlev   ! DO loop over time levels
!_nu   DO ks = 0, ke_snow  ! DO loop over snow layers
!_nu     WHERE (depth_lk(:,:) > 0.0_ireals) t_snow_mult(:,:,ks,nt) = t_snow(:,:,nx)
!_nu   END DO
!_nu   END DO
!_nu ENDIF
!_nu 
!_nu ! If multi-layer soil model is used, set "t_so(:,:,:,:)" at lake points
!_nu ! to the fresh-water temperature of maximum density at all levels except 0.
!_nu ! Notice that the climatological soil temperature,
!_nu ! i.e. t_so(:,:,ke_soil+1,:), is left intact.
!_nu IF (lmulti_layer) THEN
!_nu   DO nt = 1, nztlev   ! DO loop over time levels
!_nu   DO ks = 1, ke_soil  ! DO loop over soil layers
!_nu     WHERE (depth_lk(:,:) > 0.0_ireals) t_so(:,:,ks,nt) = tpl_T_r
!_nu   END DO
!_nu   END DO
!_nu ENDIF
!_nu  
!_tmp<

!_dbg>
!_nu    WRITE(*,*) 'FLake: WARM START initialization completed. '
!_dbg<

    !===============================================================================================
    !  END calculations
    !----------------------------------------------------------------------------------------------

  END SUBROUTINE flake_init

!===================================================================================================
!  End of the lake parameterization scheme FLake initialization
!---------------------------------------------------------------------------------------------------

!_XXX

!===================================================================================================
!  FLake interface. Advances the lake model variables one time step ahead.
!---------------------------------------------------------------------------------------------------

  SUBROUTINE flake_interface 

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
! Description:
!
!  The FLake interface is a communication routine between "flake_driver" and 
!  COSMO (or any other model that uses FLake).
!  It assigns the FLake variables at the previous time step 
!  to their input values given by COSMO,
!  calls a number of routines to compute the heat and radiation fluxes,
!  calls "flake_driver",
!  and returns the updated FLake variables to COSMO.
!  The "flake_interface" does not contain any FLake physics. 
!  It only serves as a convenient means to organize calls of "flake_driver"
!  and of external routines that compute heat and radiation fluxes.
!  These calls are executed in a DO loop over all horizontal grid-points.
!  For sea and land grid-points, no action is taken so that
!  FLake variables remain equal to their reference values. 
!
!==============================================================================

!_dbg>
    WRITE(*,*) ' flake_interface is now dummy '
!_dbg<

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE flake_interface

!---------------------------------------------------------------------------------------------------
!  End of FLake interface
!===================================================================================================

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

END MODULE mo_flake

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890
