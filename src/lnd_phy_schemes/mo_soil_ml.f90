!>
!! soil vegetation atmosphere transfer: source module  "mo_soil_ml_v413.f90"
!!------------------------------------------------------------------------------
!!----------------------------------------------------------------------------- 
!! 
!! @par Description:
!!   The module "mo_soil_ml_v413.f90" performs calculations related to the
!!   parameterization of soil processes. It contains the soubroutine
!!   terra_multlay which is the combination of the former parts
!!   terra1_multlay.incf and terra2_multlay.incf of the LM.
!!   All parametric scalar and array data for this soil model routine are
!!   defined in the data module data_soil.f90.
!!
!!   All global variables of the model that are used by the soil model routine
!!   terra_multlay.incf is imported by USE statements below.
!! 
!!   The parameterization package has been provided by B. Ritter in a
!!   Plug-compatible Fortran77-Version, which is based on the EM/DM soil model
!!   by E. Heise. Some technical modifications have been done for the
!!   F90 and the parallel Version:
!!   Internal communication by common-blocks is replaced by module parameters,
!!   scalars and arrays defined in module data_soil. The plug compatible
!!   I/O lists of the subroutines have been replaced by the Module interface
!!   defined by Use lists below.
!!
!! @author E. Heise, R. Schrodin, B. Ritter  
!! @author E. Machulskaya, F. Ament, J. Helmert
!!
!! @par reference   This is an adaption of subroutine hydci_pp in file src_gscp.f90
!!  of the lm_f90 library (COSMO code). Equation numbers refer to
!!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!!
!! $Id: n/a$
!!
!! @par Revision History
!! implemented into ICON by K. Froehlich, E. Machulskaya, and J. Helmert (2010-11-XX)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
! <li> The code may not be re-distributed without the consent of the authors.
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
!!==============================================================================
!
!
!
! Current Code Owner: DWD, Reinhold Schrodin
!  phone:  +49  69  8062 2709
!  fax:    +49  69  8062 3721
!  email:  reinhold.schrodin@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 2.17       2002/05/08 Reinhold Schrodin
!  Initial release
! 2.18       2002/07/16 Reinhold Schrodin
!  Corrections and further developments
! 3.6        2003/12/11 Reinhold Schrodin
!  Reformulation of bare soil evaporation, transpiration, snow density, 
!  snow height, snow melting, subsoil runoff calculation
! 3.7        2004/02/18 Reinhold Schrodin / Ulrich Schaettler
!  Adapted use of nold to 2 time level scheme and some corrections
! 3.13       2004/12/03 Reinhold Schrodin, Bodo Ritter, Erdmann Heise
!  Snow albedo:   Snow age indicator calculation (freshsnow)
!  Transpiration: 2m temperature and 10m wind considered instead temperature
!                 and wind at lowest atmospheric main level
!  Runoff:        Reformulation of soil water runoff (gravitational part)
!  Soil water:    Soil water calculations restricted to ke_soil_hy,
!                 determined by soil level nearest 2.5 m, combined with
!                 a reformulation of the tridiagonal equation system
!  Snow melting:  Correction in case of pore volume overshooting of
!                 uppermost soil layer by infiltration of snow water
!  Precautions to avoid
!   - soil water content below air dryness point
!   - soil surface temperature problems at lateral LM boundaries.
!   - excessive turbulent fluxes at the soil surface due to mismatch
!     of atmospheric and soil variables in the initial state
!   - evaporation if soil water content is below air dryness point
! 3.16       2005/07/22 Reinhold Schrodin
!   - Modification of evaporation treatment at plant wilting point
!   - Soil water transport: Soil water diffusion term driven by the gradient
!                        of the total soil water concentration (water + ice)
!   - Combining terra1_multlay and terra2_multlay in terra_multlay
! 3.17       2005/12/12 Reinhold Schrodin
!   LME adaptations to GME soil model modifications (B. Ritter):
!   - Constant snow density (250 kg/m**3) replaced by prognostic snow density
!     formulation (min = 50 kg/m**3, max = 400 kg/m**3)
!   - Adjustment of soil temperatures in presence of snow
!   - Extended formulation of transfer coefficient limitation
! 3.18       2006/03/03 Ulrich Schaettler / Erdmann Heise
!  Adaptations to do some initializations also for restart runs
!  Field h_snow is now time dependent, because of use in the FLake model
!  Increase infiltration and reduce surface runoff (bugfix) (Erdmann Heise)
! 3.19       2006/04/25 Erdmann Heise
!  Remove spurious snow and set consistent t_snow and w_snow at the beginning
! 3.21       2006/12/04 Ulrich Schaettler
!  crsmin, rat_lam put to data_soil; eliminated variables that are not used
!  Use dt for ntstep=0 instead of dt/2
!  Additional use of graupel, if present (Thorsten Reinhardt)
!  Modifications to avoid soil water content below air dryness point
!  Gravitational runoff calculation modified in case of frozen soil layer
!                                                      (Reinhold Schrodin)
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introduction of tfv to consider different laminar resistance for heat
!  and water vapour; thus 'rat_lam' is not used here any longer.
! V4_4         2008/07/16 Ulrich Schaettler
!  Splitting of a loop in Section I.4.3b (m_styp is not defined for sea points 
!  and must not occur together with llandmask in the same IF-clause)
! V4_7         2008/12/12 Ulrich Schaettler
!  There were still some loops left with llandmask and m_styp in one line
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Inserted Call to collapse loops
! V4_10        2009/09/11 Christian Bollmann
!  Added compiler directive to use option _on_adb for NEC
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Lucio Torrisi
!  Implementation of multi-layer snow model (EM)
!  Use of an external parameter field for stomata resistance (JH)
!  Implementation of ground water as lower boundary of soil column and
!    soil moisture dependent heat conductivity of the soil (JH)
!  Save additional fluxes and stomata resistance to global memory for output (LT)
! V4_12        2010/05/11 Ulrich Schaettler, Ekaterina Machulskaya
!  Renamed t0 to t0_melt because of conflicting names
!  Renamed prs_min to rsmin2d because of conflicting names
!  Update of the new snow model (EM)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! 
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

MODULE mo_soil_ml
!
! Declarations:
!
! Modules used:

#ifdef __COSMO__ 

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &
! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    ke,           & ! number of grid points in vertical direction
    ke_soil,      & ! number of layers in multi-layer soil model
    ke_snow,      & ! number of layers in multi-layer soil model
    czmls,        & ! depth of the main level soil layers in m
                    ! (see organize_physics.f90)
! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt              ! long time-step
!   dt2             ! dt*2.

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! absolute zero for temperature
    r_d,          & ! gas constant for dry air    
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d -1  
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp,        & ! r_d / cp_d
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion         
    lh_s,         & ! latent heat of sublimation    
    g,            & ! acceleration due to gravity
    sigma,        & ! Boltzmann-constant

! 2. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    rho_w           ! density of liquid water  

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! base state pressure                           (Pa) 

! 2. external parameter fields                                        (unit)
! ----------------------------
    soiltyp_subs    ,    & ! type of the soil (keys 0-9)                     --
    plcov      ,    & ! fraction of plant cover                         --
    rootdp     ,    & ! depth of the roots                            ( m  )
    sai        ,    & ! surface area index                              --
    tai        ,    & ! transpiration area index                        --
    eai        ,    & ! earth area (evaporative surface area) index     --
    llandmask  ,    & ! landpoint mask                                  --
    rsmin2d    ,    & ! minimum stomata resistance                    ( s/m )


! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    qv         ,    & ! specific water vapor content                  (kg/kg)
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_snow    ,     & ! temperature of the snow-surface               (  K  )
    t_snow_mult,    & ! temperature of the snow-surface               (  K  )
    t_s       ,     & ! temperature of the ground surface             (  K  )
    t_g       ,     & ! weighted surface temperature                  (  K  )
    qv_s      ,     & ! specific humidity at the surface              (kg/kg)
    w_snow    ,     & ! water content of snow                         (m H2O)
    rho_snow  ,     & ! snow density                                  (kg/m**3)
    rho_snow_mult,  & ! snow density                                  (kg/m**3)
    h_snow    ,     & ! snow height                                   (  m  
    w_i       ,     & ! water content of interception water           (m H2O)
    t_so      ,     & ! soil temperature (main level)                 (  K  )
    w_so      ,     & ! total water conent (ice + liquid water)       (m H20)
    swi       ,     & ! soil wetness index                            (  1  )
    w_so_ice  ,     & ! ice content                                   (m H20)
    t_2m      ,     & ! temperature in 2m                             (  K  )
!    u_10m     ,     & ! zonal wind in 10m                             ( m/s )
!    v_10m     ,     & ! meridional wind in 10m                        ( m/s )
    uv_low     ,     & ! wind speed at lowest model level              ( m/s ) 
    freshsnow ,     & ! indicator for age of snow in top of snow layer(  -  )
    wliq_snow ,     & ! liquid water content in the snow              (m H2O)
    wtot_snow ,     & ! total (liquid + solid) water content of snow  (m H2O)
    dzh_snow          ! layer thickness between half levels in snow   (  m  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields of convective and grid-scale precipitation
    prr_con     ,   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con     ,   & ! precipitation rate of snow, convective        (kg/m2*s)
    prr_gsp     ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp     ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp     ,   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)

!   fields of the turbulence scheme
    tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
    tcm         ,   & ! turbulent transfer coefficient for momentum   ( -- )
    tfv         ,   & ! laminar reduction factor for evaporation      ( -- )

!   fields of the radiation
    sobs        ,   & ! solar radiation at the ground                 ( W/m2)
    thbs        ,   & ! thermal radiation at the ground               ( W/m2)
    pabs        ,   & ! photosynthetic active radiation               ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ---------------------------------------------------------------
    runoff_s    ,   & ! surface water runoff; sum over forecast      (kg/m2)
    runoff_g    ,   & ! soil water runoff; sum over forecast         (kg/m2)
    rstom       ,   & ! stomata resistance                           ( s/m )
    lhfl_bs     ,   & ! average latent heat flux from bare soil evap.( W/m2)
    lhfl_pl           ! average latent heat flux from plants         ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold  ,       & ! corresponds to ntstep - 1
    nnow  ,       & ! corresponds to ntstep 
    nnew  ,       & ! corresponds to ntstep + 1
    lmelt ,       & ! soil model with melting process
    lmelt_var,    & ! freezing temperature dependent on water content
    lmulti_snow,  & ! run the multi-layer snow model

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitation physics
    itype_trvg,   & ! type of vegetation transpiration parameterization
    itype_evsl,   & ! type of parameterization of bare soil evaporation
    itype_tran,   & ! type of surface to atmospher transfer
    itype_root,   & ! type of root density distribution
    itype_heatcond,&! type of soil heat conductivity
    itype_hydbound,&! type of hydraulic lower boundary condition
    lstomata,     & ! map of minimum stomata resistance

! 5. additional control variables
! --------------------------
     l2tls           ! forecast with 2-TL integration scheme 

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_io,            ONLY:  &
    lana_rho_snow   ! if .TRUE., take rho_snow-values from analysis file
                    ! else, it is set in the model

!------------------------------------------------------------------------------

USE data_parallel,      ONLY:  &
    my_cart_id        ! rank of this subdomain in the global communicator

!------------------------------------------------------------------------------

USE data_soil       ! All variables from data module "data_soil are used by
                    ! this module. These variables start with letter "c"

!------------------------------------------------------------------------------

! External routines used by this module

USE environment,     ONLY : collapse
USE meteo_utilities, ONLY : tgcom

#endif

#ifdef __ICON__ 
!
USE mo_kind,               ONLY: ireals=>wp,    &
                                 iintegers=>i4
USE mo_math_constants    , ONLY: pi
!
USE mo_physical_constants, ONLY: t0_melt => tmelt,& ! absolute zero for temperature
                                 r_v   => rv    , & ! gas constant for water vapour
                                 r_d   => rd    , & ! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & ! r_v/r_d - 1
                                 o_m_rdv        , & ! 1 - r_d/r_v
                                 rdv            , & ! r_d / r_v
                                 lh_v  => alv   , & ! latent heat of vapourization
                                 lh_s  => als   , & ! latent heat of sublimation
                                 lh_f  => alf   , & ! latent heat of fusion
                                 cp_d  => cpd   , & ! specific heat of dry air at constant press
                                 cpdr  => rcpd  , & ! (specific heat of dry air at constant press)^-1
                                 g     => grav  , & ! acceleration due to gravity
                                 sigma => stbo  , & ! Boltzmann-constant
                                 rdocp => rd_o_cpd  ! r_d / cp_d
!
USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --
!
USE mo_cuparameters,       ONLY: rho_w => rhoh2o ! density of liquid water (kg/m^3)
!
USE mo_phyparam_soil   
!
USE mo_lnd_nwp_config,     ONLY: nztlev, lmelt, lmelt_var, lmulti_snow,  &
  &                              itype_gscp, itype_trvg, itype_evsl,     &
  &                              itype_tran, itype_root, itype_heatcond, &
  &                              itype_hydbound, lstomata, l2tls,        &
  &                              lana_rho_snow, itype_subs
!
USE mo_lnd_nwp_config,     ONLY: t_tiles
!                           
USE mo_exception,          ONLY: message, finish, message_text
USE mo_run_config,         ONLY: msg_level
#endif



!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

!IMPLICIT NONE
!PRIVATE:: normalize, tgcom
PRIVATE:: tgcom

!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: terra_multlay,terra_multlay_init

!------------------------------------------------------------------------------
! Public variables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef COSMO
CHARACTER(132) :: message_text = ''
#endif


CONTAINS

#ifdef COSMO
SUBROUTINE message (name, text, all_print)
IMPLICIT NONE

CHARACTER (*) :: name, text
LOGICAL, INTENT(in), OPTIONAL :: all_print

LOGICAL :: lprint

IF (PRESENT(all_print)) THEN
  lprint = all_print
ELSE
  lprint = .FALSE.
ENDIF

IF (lprint .OR. my_cart_id==0) &
  WRITE (*,*) TRIM(name), TRIM(text)

END SUBROUTINE message

#endif



!==============================================================================
!+ Computation of the first part of the soil parameterization scheme
!------------------------------------------------------------------------------

!option! -pvctl _on_adb


  SUBROUTINE terra_multlay (         &   
                  ie               , & ! array dimensions
                  istartpar        , & ! start index for computations in the parallel program
                  iendpar          , & ! end index for computations in the parallel program
                  nsubs0           , & ! nsubs0=1 for single tile, nsubs0=2 for multi-tile
                  nsubs1           , & ! nsubs1=1 for single tile, nsubs1=#tiles+1 for multi-tile
                  ldiag_tg         , & ! if true: diagnose t_g and snow-cover fraction with tgcom
                  ke_soil, ke_snow , &
                  czmls            , & ! processing soil level structure 
                  dt               , &
!
                  soiltyp_subs     , & ! type of the soil (keys 0-9)                     --
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  eai              , & ! earth area (evaporative surface area) index     --
!                 llandmask        , & ! landpoint mask                                  --
                  rsmin2d          , & ! minimum stomata resistance                    ( s/m )
!                                  
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  p0               , & !!!! base state pressure                        ( Pa  ) 
!                 pp               , & ! deviation from the reference pressure         ( Pa  )
                  ps               , & ! surface pressure                              ( Pa  )
!                                  
                  t_snow_now       , & ! temperature of the snow-surface               (  K  )
                  t_snow_new       , & ! temperature of the snow-surface               (  K  )
!
                  t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_new  , & ! temperature of the snow-surface               (  K  )
!
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_s_new          , & ! temperature of the ground surface             (  K  )
!
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s             , & ! specific humidity at the surface              (kg/kg)
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  w_snow_new       , & ! water content of snow                         (m H2O)
!
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  rho_snow_new     , & ! snow density                                  (kg/m**3)
!
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  rho_snow_mult_new, & ! snow density                                  (kg/m**3)
!
                  h_snow           , & ! snow height                                   (  m  )
!
                  w_i_now          , & ! water content of interception water           (m H2O)
                  w_i_new          , & ! water content of interception water           (m H2O)
!
                  t_so_now         , & ! soil temperature (main level)                 (  K  )
                  t_so_new         , & ! soil temperature (main level)                 (  K  )
!
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
!
                  w_so_ice_now     , & ! ice content                                   (m H20)
                  w_so_ice_new     , & ! ice content                                   (m H20)
!
                  t_2m             , & ! temperature in 2m                             (  K  )
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow          , & ! snow-cover fraction                           (  -  )
!
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
!
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  wtot_snow_new    , & ! total (liquid + solid) water content of snow  (m H2O)
!
                  dzh_snow_now     , & ! layer thickness between half levels in snow   (  m  )
                  dzh_snow_new     , & ! layer thickness between half levels in snow   (  m  )
!
                  subsfrac         , &
!                                  
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  prg_gsp          , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
!                                  
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv              , & ! laminar reduction factor for evaporation      ( -- )
!                                  
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs             , & !!!! photosynthetic active radiation            ( W/m2)
!                                  
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g         , & ! soil water runoff; sum over forecast          (kg/m2)
!                 ntstep             & ! actual time step
                  pt_tiles         , & ! tiles structure
!
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2) 
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2) 
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2) 
                  zlhfl_snow         & ! latent   heat flux snow/air interface         (W/m2) 
                                     )


!-------------------------------------------------------------------------------
! Declarations for ICON (USE statements for COSMO)
!-------------------------------------------------------------------------------

IMPLICIT NONE



  INTEGER (KIND=iintegers), INTENT(IN)  ::  &
                  ie,                & ! array dimensions
                  istartpar,         & ! start index for computations in the parallel program
                  iendpar,           & ! end index for computations in the parallel program
                  nsubs0,nsubs1    , &
                  ke_soil, ke_snow      
  REAL    (KIND = ireals), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  czmls                ! processing soil level structure 
  LOGICAL, INTENT(IN) :: ldiag_tg      ! if .TRUE., use tgcom to diagnose t_g and snow-cover fraction
  REAL    (KIND = ireals), INTENT(IN)  ::  &
                  dt                    

  INTEGER (KIND=iintegers),DIMENSION(ie), INTENT(IN) :: & 
                  soiltyp_subs         ! type of the soil (keys 0-9)                     --
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  eai                  ! earth area (evaporative surface area) index     --
!  LOGICAL                , DIMENSION(ie), INTENT(IN) :: & 
!                  llandmask           ! landpoint mask                                  --
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                 rsmin2d               ! minimum stomata resistance                    ( s/m )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv                   ! specific water vapor content                  (kg/kg)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                  p0                   !!!! base state pressure                        ( Pa ) 
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) ::    &
                  ps                   ! surface pressure                              ( pa  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_snow_now           ! temperature of the snow-surface (K)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  t_snow_new
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_snow), INTENT(OUT) :: &
                  t_snow_mult_new      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_s_now              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  t_s_new              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie)  :: &
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s                 ! specific humidity at the surface              (kg/kg)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  w_snow_new       , & ! water content of snow                         (m H2O)
                  rho_snow_new         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(OUT) :: &
                  rho_snow_mult_new    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  h_snow               ! snow height  
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_i_now              ! water content of interception water           (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  w_i_new              ! water content of interception water           (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_now         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_new         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_2m             , & ! temperature in 2m                             (  K  )
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow              ! snow-cover fraction
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(OUT) :: &
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_new        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(OUT) :: &
                  dzh_snow_new         ! layer thickness between half levels in snow   (  m  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  subsfrac
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) ::    &
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp              ! precipitation rate of snow, grid-scale        (kg/m2*s)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) ::    &
                  prg_gsp              ! precipitation rate of graupel, grid-scale     (kg/m2*s)

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv                  ! laminar reduction factor for evaporation      ( -- )

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs                 !!!! photosynthetic active radiation            ( W/m2)

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g             ! soil water runoff; sum over forecast          (kg/m2)
  TYPE(t_tiles), TARGET,                  INTENT(IN) :: &
                  pt_tiles(nsubs1)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2) 
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2) 
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2) 
                  zlhfl_snow           ! latent   heat flux snow/air interface         (W/m2) 

!!$  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
!!$                  rstom            ! stomata resistance                           ( s/m )
!!$  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
!!$                  lhfl_bs          ! average latent heat flux from bare soil evap.( W/m2)
!!$  REAL    (KIND = ireals), DIMENSION(ie,ke_soil), INTENT(INOUT) :: &
!!$                  lhfl_pl          ! average latent heat flux from plants         ( W/m2)


!--------------------------------------------------------------------------------
! TERRA Declarations

! New declaration for ICON

!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------

  INTEGER (KIND=iintegers)              ::  &
    ierror                        ! error status variable

  CHARACTER (LEN=80)                    ::  &
    yerror
!
! Local parameters:
! ----------------

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals , & ! security constant
    zalfa  = 1.0_ireals    , & ! degree of impliciteness (1: full implicit,
                               !    (0.5: Cranck-Nicholson)
    rho_i  = 910._ireals       ! density of solid ice (soil model)  (kg/m**3)

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
!
!   Indices
!
!   nx             , & ! time-level for integration
    kso            , & ! loop index for soil moisture layers           
    ksn            , & ! loop index for snow layers
    k              , & ! loop index for snow layers
    ke_soil_hy     , & ! number of active soil moisture layers
    i              , & ! loop index in x-direction              
#ifndef __ICON__
!    im1, jm1       , & ! i-1, j-1   ! must be removed completely
#endif
    jb             , & ! loop index for soil-type               
    mstyp          , & ! soil type index
    msr_off        , & ! number of layers contributing to surface run off
    istarts        , & ! start index for x-direction      
    iends          , & ! end   index for x-direction     
    i250           , & ! number of layers above a soil depth of 2.50 m
    k10cm          , & ! index of half level closest to 0.1 m
    k100cm         , & ! index of half level closest to 1.0 m
    ns

  REAL    (KIND=ireals   ) ::  &
!
!   Timestep parameters
!
    zdt            , & ! integration time-step [s]
    zdelt_bound    , & ! for surface temperature problem at lateral boundaries
!   zdel_t_so      , & ! auxiliary variable
    zdtdrhw        , & ! timestep / density of water
    zrhwddt        , & ! density of water / timestep
    zroffdt        , & ! timestep for runoff calculation
    z1d2dt         , & ! 1./2*timestep
!
!   Connection to the atmosphere
!
    zuv            , & ! wind velocity in lowest atmospheric layer
    ztvs           , & ! virtual temperature at the surface
    zplow          , & ! pressure of lowest atmospheric layer 
    zqvlow         , & ! specific humidity of lowest atmospheric layer
    zrss           , & ! ice covered part of grid element
    zrww           , & ! water covered part of grid element
!
!   Snow parameters
!
    zsnow_rate     , & ! rate of snow fall            [kg/m**2 s]
    zdsn_new       , & ! snow age refresh increment   [-]
    zdsn_old       , & ! snow age decay increment     [-]
    zdz_snow(ie)    ! snow depth (not for snow heat content but snow
                       ! density calculation)

  REAL    (KIND=ireals   ) ::  &
!
!   Multi-snow layer parameters
    zsn_porosity   , & ! snow porosity                    (   -   )
    zp1            , &
    zfukt          , &
    zq0            , &
    zqbase(ie)  , &
    zdzh_old       , &
    zrefr(ie)   , & ! rate of liquid water refreezing
    zmelt(ie)   , & ! rate of snow melting
    ze_in          , &
    ze_out(ie)  , &
    zadd_dz        , &
    zrho_dry_old(ie)   , &
    zeta           , &
    zdens_old      , &
    zp(ie,ke_snow), &
    zcounter(ie), &
    ze_rad(ie)  , &
    zswitch(ie) , &
    fact1          , &
    fact2          , &
    tmp1           , &
    tmp2           , &
    tmp3           , &
    zf_snow_old(ie)    , &
    tmp_num        , &
    sum_weight(ie)      , &
    t_new  (ie,ke_snow) , &
    rho_new(ie,ke_snow) , &
    wl_new (ie,ke_snow) , &
    z_old  (ie,ke_snow) , & 
    dz_old (ie,ke_snow) , &
    weight         , &
!
!   Plant parameters
    zbeta          , & ! reduction factor for evaporation
    zevap          , & ! auxiliary variable
    zalpha              ! NP89 bare soil evaporation

  REAL    (KIND=ireals   ) ::  &
!
!   Local scalars for BATS-scheme
!
    zpsi0=.01_ireals,& ! air entry potential at water saturation (m)
    zbf1           , & ! auxiliary variable
    zbf2           , & ! auxiliary variable
    zdmax          , & ! auxiliary variable
    zd             , & ! auxiliary variable
    zck            , & ! auxiliary variable
    z1             , & ! depth of half level closest to 0.1 m
    znull          , & ! depth of half level closest to 1.0 m
    zfqmax         , & ! maximum sustainable flux in uppermost soil layer
    zevapor        , & ! evaporation rate of bare soil in BATS-scheme
    zrla           , & ! atmospheric resistance
    zcatm          , & ! transfer function CA
    zpar           , & ! PAR (interim value)
    zf_wat         , & ! soil water function for stomatal resistance
    zf_tem         , & ! temperature function for stomatal resistance
    zf_sat         , & ! saturation deficit function for stomatal resistance
    zepsat         , & ! saturation vapour pressure at near surface temperature
    zepke          , & ! near surface water vapour pressure
    zrstom         , & ! stomata resistance
    zrveg          , & ! sum of atmospheric and stomata resistance
    zedrstom       , & ! 1./zrstom
    zustar         , & ! friction velocity
    ztrabpf        , & ! area average transpiration rate
    ! variables for non-uniform root-distribution
    ztrfr          , & ! fraction of total transpiration associated with each layer
    zrootfr        , & ! normalized root fraction in each layer (=1 for z_soil=0.0)
    zrootdz        , & ! normalized root fraction * min(layer thickness,rootlength-top_of_layer)
!
!   Soil parameters
!
    zice           , & ! indicator of soil type ice
    zwqg           , & ! mean of fcap and pwp
    z4wdpv         , & ! 4*zwqg/porv         
!   zroot          , & ! fraction of layer filled by roots (Bucket scheme)
    zropart        , & ! fraction of layer filled by roots (Dickinson scheme)
    zrootfc        , & ! distributes transpiration to soil layers
!   zr_root        , & ! zroot/zbwt (Bucket scheme)
    ze_sum             ! sum of all contributions to evapotranspiration

  REAL    (KIND=ireals   ) ::  &
!
!   Thermal parameters
!
!    ztgt0          , & ! Indicator T_g > T_0
    zgstr          , & ! downward longwave radiation
    zalas          , & ! heat conductivity of snow
    zrnet_snow     , & ! net radiation at snow surface
    zfak           , & ! utility variable for implicit snow temperature forecast
    ztsnow_im      , & ! utility variable for implicit snow temperature forecast
    ztsnew         , & ! preliminary value of snow surface temperature
    ztsnownew      , & ! preliminary value of snow surface temperature
    zfor_snow      , & ! total forcing at snow surface
    zfr_melt       , & ! melting snow fraction
!    zdwsnm         , & ! utility variable for snow melt determination
    zdwgme         , & ! utility variable for snow melt determination
    zdelt_s        , & ! utility variable for snow melt determination
    ze_avail       , & ! utility variable for snow melt determination
    ze_total           ! utility variable for snow melt determination

  ! Some of these values have to be fields for the multi-layer snow model
  REAL    (KIND=ireals   ) ::  &
    zalas_mult    (ie,ke_snow),    & ! heat conductivity of snow
    ztsnownew_mult(ie,0:ke_snow),  & ! preliminary value of snow surface temperature
    zextinct      (ie,ke_snow),    & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(ie)               ! total forcing at snow surface

! REAL    (KIND=ireals   ) ::  &
!   zfor_total    (ie)               ! total forcing at surface

  REAL    (KIND=ireals   ) ::  &
!
!   Hydraulic parameters
!
    zwinstr        , & ! preliminary value of interception store
    zinfmx         , & ! maximum infiltration rate
    zwimax         , & ! maximum interception store
    zalf           , & ! utility variable
    zvers          , & ! water supply for infiltration
    zro_inf        , & ! surface runoff
    zdwsndtt       , & ! time rate of change of snow store
    zwsnstr        , & ! preliminary value of snow store
    zdwseps        , & ! artificial change of small amount of snow store
    zdwidtt        , & ! time rate of change of interception store
    zdwieps        , & ! artificial change of small amount of interception store
    zro_wi         , & ! surface runoff due to limited infiltration capacity
    zro_sfak       , & ! utility variable 
    zro_gfak       , & ! utility variable
    zfmb_fak       , & ! utility variable
    zdwg           , & ! preliminary change of soil water content
    zwgn           , & ! preliminary change of soil water content
    zredfu         , & ! utility variable for runoff determination
    zro            , & ! utility variable for runoff determination
    zro2           , & ! utility variable for runoff determination
    zkorr          , & ! utility variable for runoff determination
    zfr_ice        , & ! reduction factor for water transport
    zfr_ice_free   , & ! reduction factor for water transport
    zwso_new       , & ! preliminary value of soil water content
!    zwsnew         , & ! preliminary value of snow water equivalent
    zw_ovpv        , & ! utility variable
!
!   Implicit solution of thermal and hydraulic equations
!
    zakb           , & ! utility variable
    zzz            , & ! utility variable
    z1dgam1        , & ! utility variable
    zredm          , & ! utility variable
    zredm05        , & ! utility variable
    zredp05        , & ! utility variable
    zgam2m05       , & ! utility variable
    zgam2p05           ! utility variable

  REAL    (KIND=ireals   ) ::  &
!
!   Statement functions
!
    zsf_heav       , & ! Statement function: Heaviside function
    zsf_psat_iw    , & ! Saturation water vapour pressure over ice or water
                       ! depending on temperature "zstx"
    zsf_qsat       , & ! Specific humidity at saturation pressure
                       ! (depending on the saturation water vapour pressure
                       !  "zspsatx" and the air pressure "zspx")
    zsf_dqvdt_iw   , & ! Statement function: First derivative of specific
                       ! saturation humidity
                       ! with respect to temperature (depending on temperature
                       ! "zstx" and saturation specific humidity pressure
                       ! "zsqsatx")
    zstx           , & ! dummy argument for Stmt. function
    zspx           , & ! dummy argument for Stmt. function
    zspsatx        , & ! dummy argument for Stmt. function
    zsqsatx        , & ! dummy argument for Stmt. function
    z2iw           , & ! dummy argument for Stmt. function
    z4iw           , & ! dummy argument for Stmt. function
    z234iw         , & ! dummy argument for Stmt. function
    zqs            , & ! saturation humidty at T_s and ps
    zdqs           , & ! derivative of zqs with respect to T_s
    zqsnow         , & ! saturation humidty at T_snow and ps
    zdqsnow        , & ! derivative of zqs with respect to T_snow
!   Local scalars for transfer coefficient limitation calculations and for
!   treatment of possible problems at lateral boundaries
    zch_soil = 1.400E06_ireals ,& ! approximate average heat capacity of soil (J/m**3 K)
    zlim_dtdt = 2.5_ireals   ,&! maximum allowed temperature increment (K) in
                               ! one time step in  uppermost soil layer
    zdT_snow                            ,& ! maximum possible increment in snow layer before
                                           ! melting sets in
    zg1       (ie                  ) ,& ! heat flux between two uppermost soil layers (W/m**2)
    zlhfl     (ie                  ) ,& ! estimated       latent    heat flux surface (W/m**2)
    zshfl     (ie                  ) ,& ! estimated       sensible  heat flux surface (W/m**2)
    zthfl     (ie                  ) ,& ! estimated total turbulent heat flux surface (W/m**2)
    zradfl    (ie                  ) ,& ! total radiative flux at surface             (W/m**2)
    ze_melt   (ie                  ) ,& ! energy required for melting of snow         (J/m**2)
    zch_snow  (ie                  ) ,& ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
    zeb1      (ie                  ) ,& ! estimated energy budget of first soil layer (W/m**2)
    ztchv     (ie                  ) ,& ! transfer coefficient multiplied by wind velocity
    ztchv_max (ie                  ) ,& ! dito, but upper limit
    zrho_atm  (ie                  ) ,& ! air density of lowest atmospheric layer     (kg/m**3)
    zdt_atm   (ie                  ) ,& ! temperature difference between lowest
    zdq_atm   (ie                  ) ,& ! dito for moisture (kg/kg)
    !additional variables for root distribution
    zroota    (ie                  ) ,& ! root density profile parameter (1/m)
    zwrootdz  (ie         , ke_soil) ,& ! mean water content over root depth weighted by root density
    zrootdz_int (ie                ) ,& ! parameter needed to initialize the root density profile integral
    zwrootdz_int(ie                )    ! parameter needed to initialize the root water content integral
    
    
  INTEGER m_limit                          ! counter for application of limitation
  CHARACTER *7 yhc                         ! heating or cooling indicator


  REAL    (KIND=ireals   ) ::  &
!
!   Local scalars for water content dependent freezing/melting
!
    zliquid        , & ! utility variable
    zaa            , & ! utility variable
    znen           , & ! utility variable
    zargu          , & ! utility variable
!
!   Water transport
!
    zice_fr_ksom1  , & ! fractional ice content of actual layer - 1
    zice_fr_kso    , & ! fractional ice content of actual layer
    zice_fr_ksop1  , & ! fractional ice content of actual layer + 1
    zlw_fr_ksom1   , & ! fractional liquid water content of actual layer - 1
    zlw_fr_ksom1_new,& ! fractional liquid water content of actual layer -1
    zlw_fr_kso     , & ! fractional liquid water content of actual layer
    zlw_fr_kso_new , & ! fractional liquid water content of actual layer
    zlw_fr_ksop1   , & ! fractional liquid water content of actual layer + 1
    zlw_fr_ksom05  , & ! fractional liquid water content of actual layer - 1/2
    zdlw_fr_ksom05 , & ! hydraulic diffusivity coefficient at half level above
    zklw_fr_ksom05 , & ! hydraulic conductivity coefficient at half level above
    zklw_fr_kso_new, & ! hydraulic conductivity coefficient at main level
                       !    for actual runoff_g
    zlw_fr_ksop05  , & ! fractional liquid water content of actual layer + 1/2
    zdlw_fr_ksop05 , & ! hydraulic diffusivity coefficient at half level below
    zklw_fr_ksop05 , & ! hydraulic conductivity coefficient at half level below
    zinf           , & ! infiltration 
!
!   Snow density
    ztau_snow      , & ! 'ageing constant' for snow density          (-)
    zrho_snowe     , & ! updated density of existing snow ('ageing') kg/m**3
    zrho_snowf     , & ! density of fresh snow                       kg/m**3
    znorm          , & ! normalisation factor for weighted snow density mH2O
!
!
!   Freezing/melting of soil water/ice:
!
    zenergy        , & ! available melting/freezing energy per time step
    zdelwice       , & ! amount of melted soil ice/frozen soil water
    zdwi_max       , & ! maximum amount of freezing/melting soil water
    ztx                ! water content dependent freezing/melting temperature

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND = ireals) :: &
!
! Two-time level variables exchange with interface
!
    h_snow_now (ie)    , & ! snow height  (m)  
    h_snow_new (ie)    , & ! snow height  (m)  
!
!
! Model geometry
!
    zmls     (ke_soil+1)  , & ! depth of soil main level
    zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
    zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
    zdzms    (ke_soil+1)  , & ! distance between main levels
    zdz_snow_fl(ie)    , & ! snow depth for snow temperature gradient
!
! Multi-layer snow model
    zhh_snow (ie,ke_snow)    , & ! depth of the half level snow layers
    zhm_snow (ie,ke_snow)    , & ! depth of snow main levels
    zdzh_snow(ie,ke_snow)    , & ! layer thickness between half levels
    zdzm_snow(ie,ke_snow)    , & ! distance between main levels
!
! External parameters
!
    zbwt     (ie)      , & ! root depth (with artificial minimum value)
    zrock    (ie)      , & ! ice/rock-indicator: 0 for ice and rock
    zsandf   (ie)      , & ! mean fraction of sand (weight percent)
    zclayf   (ie)      , & ! mean fraction of clay (weight percent)
    zb_por   (ie)      , & ! pore size distribution index
    zpsis    (ie)      , & ! air entry potential (m)
    zw_m     (ie)          ! maximum of liquid water content  (m)

  INTEGER  (KIND=iintegers ) ::  &
    m_styp   (ie)          ! soil type

  REAL    (KIND=ireals   ) ::  &
!
! Connection to the atmosphere
!
    zrr      (ie)      , & ! total rain rate including formation of dew
    zrs      (ie)      , & ! total snow rate including formation of rime
    zesoil   (ie)      , & ! evaporation from bare soil
    zrhoch   (ie)      , & ! transfer coefficient*rho*g
    zth_low  (ie)      , & ! potential temperature of lowest layer
    zf_wi    (ie)      , & ! surface fraction covered by interception water
    ztmch    (ie)      , & ! heat transfer coefficient*density*velocity
    zep_s    (ie)      , & ! potential evaporation for t_s
    zep_snow (ie)      , & ! potential evaporation for t_snow
    zverbo   (ie)      , & ! total evapotranspiration
    zversn   (ie)      , & ! total evaporation from snow surface
    zthsoi   (ie)      , & ! thermal flux at soil surface
    zthsnw   (ie)      , & ! thermal flux at snow surface
    zfor_s   (ie)      , & ! total forcing at soil surface
    zgsb     (ie)      , & ! heat-flux through snow
    zrnet_s  (ie)      , & ! net radiation
    zsprs    (ie)      , & ! utility variable
!
! Tendencies
!
    zdwidt   (ie)      , & ! interception store tendency
    zdwsndt  (ie)      , & ! snow water content tendency
    zdtsdt   (ie)      , & ! tendency of zts
    zdtsnowdt(ie)      , & ! tendency of ztsnow
    zdtsnowdt_mult(ie,0:ke_snow)      , & ! tendency of ztsnow
    zdwgdt   (ie,ke_soil)  ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=ireals   ) ::  &
!
! Soil and plant parameters
!
    zroc(ie,ke_soil+1) , & ! heat capacity
    zfcap    (ie)      , & ! field capacity of soil
    zadp     (ie)      , & ! air dryness point
    zporv    (ie)      , & ! pore volume (fraction of volume)
    zdlam    (ie)      , & ! heat conductivity parameter
    zdw      (ie)      , & ! hydrological diff.parameter
    zdw1     (ie)      , & ! hydrological diff.parameter
    zkw      (ie)      , & ! hydrological cond.parameter
    zkw1     (ie)      , & ! hydrological cond.parameter
    zik2     (ie)      , & ! minimum infiltration rate
    zpwp     (ie)      , & ! plant wilting point  (fraction of volume)
    ztlpmwp  (ie)      , & ! turgor-loss-point minus plant wilting point
    zedb     (ie)      , & ! utility variable
!
! Hydraulic variables
!
    ztrang(ie,ke_soil) , & ! transpiration contribution by the different layers
    ztrangs(ie)        , & ! total transpiration (transpiration from all
                              !    soil layers)
    zwin     (ie)      , & ! water cont. of interception store   (m H20)
    zwsnow   (ie)      , & ! snow water equivalent               (m H20)
    zwsnew   (ie)      , & ! snow water equivalent               (m H20)
    zdwsnm   (ie)      , & ! utility variable for snow melt determination
    zw_fr    (ie,ke_soil+1),&!fractional total water content of soil layers
    zinfil   (ie)      , & ! infiltration rate
    zlw_fr   (ie,ke_soil+1), & ! fractional liqu. water content of soil layer
    ziw_fr   (ie,ke_soil+1), & ! fractional ice content of soil layer
    zwsnn    (ie)          , & ! new value of zwsnow
    zflmg    (ie,ke_soil+1), & ! flux of water at soil layer interfaces
    zrunoff_grav(ie,ke_soil+1) ! main level water gravitation
                                  !    flux (for runoff calc.)

  REAL    (KIND=ireals   ) ::  &
!
! Local arrays for BATS-scheme
!
    zk0di    (ie)      , & ! surface type dependent parameter
    zbedi    (ie)      , & ! surface type dependent parameter
    zsnull   (ie)      , & ! mean relative(zporv) water fraction of the
                              !    so-called total active soil layer (above 1m)
    zs1      (ie)      , & ! mean relative (zporv) water fraction  of
                              !    layers above 0.1m
    zf_rad   (ie)      , & ! radiation function for stomata resistance
    ztraleav (ie)      , & ! transpiration rate of dry leaves
!
! Root functions
!
    zwroot   (ie)      , & ! mean water content over root depth
    zropartw(ie,ke_soil)   ! fraction of layer filled by roots * w_fr

  REAL    (KIND=ireals   ) ::  &
!
! Thermal variables
!
    zts      (ie)      , & ! soil surface temperature
    ztsnow   (ie)      , & ! snow surface temperaure
    ztsnow_mult   (ie,0:ke_snow)      , & ! snow surface temperaure
! HEATCOND (soil moisture dependent heat conductivity of the soil)
    zalamtmp (ie)      , & ! heat conductivity
    zalam    (ie,ke_soil),&! heat conductivity
    zrocg    (ie)      , & ! volumetric heat capacity of bare soil
    zrocs    (ie)      , & ! heat capacity of snow
    ztsn     (ie)      , & ! new value of zts
    ztsnown  (ie)      , & ! new value of ztsnow
    ztsnown_mult  (ie,0:ke_snow)      , & ! new value of ztsnow
    znlgw1f  (ke_soil)    , & ! utility variable
    zaga(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zagb(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zagc(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zagd(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zage(ie,0:ke_soil+ke_snow+1),& ! utility variable
!
! Auxiliary variables
!
!    zw_snow_old(ie)    , & !
    zdqvtsnow(ie)      , & ! first derivative of saturation specific humidity
                              !    with respect to t_snow
    zts_pm   (ie)      , & ! indicator zts > < T_melt
    ztsnow_pm(ie)          ! indicator ztsnow > < T_melt

  LOGICAL     :: &
!
    ldebug                , & !
    limit_tch (ie)         ! indicator for flux limitation problem

  ! ground water as lower boundary of soil column
  REAL    (KIND=ireals) ::  &
    zdelta_sm, zdhydcond_dlwfr, zklw_fr_kso, zklw_fr_ksom1, zdlw_fr_kso 
  
  ! HEATCOND
  REAL    (KIND=ireals) ::          &
    hzalam   (ie,ke_soil+1),     & ! heat conductivity
    zthetas, zlamli, zlamsat, zlams, rsandf, zlamq, zlam0, zrhod, zlamdry,  &
    zsri, zKe, zthliq, zlamic

!  INTEGER (KIND=iintegers) :: i_loc, isub

!- End of header
!==============================================================================


!------------------------------------------------------------------------------
! Begin Subroutine terra_multlay              
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------


! Declaration of STATEMENT-FUNCTIONS

  zsf_heav     (zstx                    ) = 0.5_ireals+SIGN( 0.5_ireals, zstx )
  zsf_psat_iw  (zstx,z2iw   ,z4iw       )                                     &
                   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))
  zsf_qsat     (zspsatx, zspx           )                                     &
                   = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)
  zsf_dqvdt_iw (zstx,zsqsatx,z4iw,z234iw)                                     &
                   = z234iw*(1._ireals+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2


!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar

!>JH
  prg_gsp=0._ireals ! graupel not implemented yet 
!<JH

  ierror = 0
  yerror = '        '

! select timelevel and timestep for calculations
! In the new formulation a 2 timelevel scheme is used for the soil model.
! Therefore it is not necessary any more to distinguish between Leapfrog
! and Runge-Kutta scheme.
! But the first timestep must not be done with dt/2, as in the Leapfrog
! scheme
!  nx  = nnow

!!$  IF ( (ntstep == 0) .AND. (.NOT. l2tls) ) THEN
!!$    ! use the original dt and not dt/2
!!$    zdt = 2._ireals * dt
!!$  ELSE
    zdt = dt
!!$  ENDIF

  ! time step for run-off computation
  zroffdt = zdt


! Computation of derived constants

  zrhwddt = rho_w/zdt     ! density of liquid water/timestep
  zdtdrhw = zdt/rho_w     ! timestep/density of liquid water

! time constant for infiltration of water from interception store
  ctau_i        = MAX(ctau_i,zdt)

! grids for temperature and water content

  zzhls(1) = 2._ireals*czmls(1)   !depth of first half level
  zdzhs(1) = zzhls(1)      !layer thickness betw. half levels of uppermost layer
  zmls(1)  = czmls(1)      !depth of 1st main level
  zdzms(1) = czmls(1)      !layer thickness between soil surface and main level
                           ! of uppermost layer
 
  DO kso = 2,ke_soil+1
    zzhls(kso)  = zzhls(kso-1) + 2._ireals*(czmls(kso) -zzhls(kso-1))
    zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
    zmls(kso)  = czmls(kso)                ! depth of main levels
    zdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
  ENDDO


!!$  IF(itype_subs .EQ. 2) THEN    ! tiles
!!$    IF (ntstep == 0) THEN
!!$      DO ns = nsubs0+1, nsubs1, 2      ! 1 - mean, 2 - ocean, 3 - lake, 4 - no snow first
!!$          DO i = istarts, iends
!!$            IF (llandmask(i,ns)) THEN  ! for landpoints only  !check is not necessary
!!$  
!!$              zf_snow_old(i) = subsfrac(i,ns)/(subsfrac(i,ns-1)+subsfrac(i,ns))
!!$  
!!$              zwsnow(i) = w_snow_now(i,ns)*zf_snow_old(i) + &
!!$                w_snow_now(i,ns-1)*(1._ireals - zf_snow_old(i))
!!$              zf_snow(i) = MAX( 0.01_ireals, MIN(1.0_ireals,zwsnow(i)/cf_snow) ) &
!!$                &            *zsf_heav(zwsnow(i) - zepsi)
!!$  
!!$              w_snow_now(i,ns  ) = zwsnow(i)/MAX(zf_snow(i),zepsi)
!!$              w_snow_now(i,ns-1) = 0._ireals
!!$  
!!$              fact1 = subsfrac(i,ns-1)+subsfrac(i,ns)
!!$              subsfrac(i,ns-1) = (1._ireals - zf_snow(i))*fact1
!!$              subsfrac(i,ns  ) = zf_snow(i)              *fact1
!!$  
!!$  
!!$            END IF  ! land-points only
!!$          END DO
!!$     END DO
!!$  !    DO ns = nsubs0+1, nsubs1, 2      ! 1 - mean, 2 - ocean, 3 - lake, 4 - no snow first
!!$  !      DO kso = 1,ke_soil
!!$  !          DO i = istarts, iends
!!$  !            IF (llandmask(i,ns)) THEN  ! for landpoints only  !check is not necessary
!!$  !              fact1 = (MIN(1._ireals - zf_snow_old(i), 1._ireals - zf_snow(i)))/MAX((1._ireals - zf_snow(i)),zepsi)
!!$  !              fact2 = (MAX(zf_snow(i) - zf_snow_old(i), 0._ireals))/MAX(zf_snow(i),zepsi)
!!$  
!!$  !IF(MOD(i,340).eq.0 .and. w_snow_now(i,ns).gt.0) &
!!$  !print *,i,kso,'snow',w_snow_now(i,ns-1),w_snow_now(i,ns  ),zf_snow    (i),fact1,fact2
!!$  
!!$  !              tmp1 = t_so    (i,kso,nx,ns-1)
!!$  !              tmp2 = w_so    (i,kso,nx,ns-1)
!!$  !              tmp3 = w_so_ice(i,kso,nx,ns-1)
!!$  
!!$  !              t_so    (i,kso,nx,ns-1) = t_so_now(i,kso,ns-1)*fact1 + t_so_now(i,kso,ns)*(1._ireals - fact1)
!!$  !              w_so    (i,kso,nx,ns-1) = w_so(i,kso,nx,ns-1)*fact1 + w_so(i,kso,nx,ns)*(1._ireals - fact1)
!!$  !              w_so_ice(i,kso,nx,ns-1) = w_so_ice(i,kso,nx,ns-1)*fact1 + w_so_ice(i,kso,nx,ns)*(1._ireals - fact1)
!!$    
!!$  !              t_so    (i,kso,nx,ns) = tmp1 * fact2 + t_so_now(i,kso,ns)*(1._ireals - fact2)
!!$  !              w_so    (i,kso,nx,ns) = tmp2 * fact2 + w_so(i,kso,nx,ns)*(1._ireals - fact2)
!!$  !              w_so_ice(i,kso,nx,ns) = tmp3 * fact2 + w_so_ice(i,kso,nx,ns)*(1._ireals - fact2)
!!$  
!!$  !IF(MOD(i,340).eq.0 .and. w_snow_now(i,ns).gt.0) &
!!$  !print *,i,kso,'temp',t_so    (i,kso,nx,ns-1),t_so    (i,kso,nx,ns)
!!$  !            END IF  ! land-points only
!!$  !          END DO
!!$  !      END DO        ! soil layers
!!$  !    END DO
!!$    END IF
!!$  END IF
  !em>


!---loop over tiles---
!DO ns=nsubs0,nsubs1
!---------------------


! Prepare basic surface properties (for land-points only)
 
    DO i = istarts, iends

!>JH
!      IF(llandmask(i)) THEN        ! for land-points only

        mstyp       = soiltyp_subs(i)        ! soil type
        m_styp(i) = mstyp                     ! array for soil type
        zdw   (i)  = cdw0  (mstyp)
        zdw1  (i)  = cdw1  (mstyp)
        zkw   (i)  = ckw0  (mstyp)
        zkw1  (i)  = ckw1  (mstyp)
        zik2  (i)  = cik2  (mstyp)
        zporv(i)  = cporv(mstyp)              ! pore volume
        zpwp (i)  = cpwp (mstyp)              ! plant wilting point
        zadp (i)  = cadp (mstyp)              ! air dryness point
        zfcap(i)  = cfcap(mstyp)              ! field capacity
        zrock(i)  = crock(mstyp)              ! rock or ice indicator
        zrocg(i)  = crhoc(mstyp)              ! heat capacity
        zdlam(i)  = cala1(mstyp)-cala0(mstyp) ! heat conductivity parameter
        zbwt(i)   = MAX(0.001_ireals,rootdp(i))! Artificial minimum value
                                                ! for root depth
        zroota(i) = 3._ireals/zbwt(i)       ! root density profile parameter (1/m)
                                                ! zroota=0. creates the original TERRA_LM
                                                ! version with constant root density
                                                ! for root depth
        ! New arrays for BATS-scheme
        zk0di(i)  = ck0di(mstyp)              !
        zbedi(i)  = cbedi(mstyp)              !
        ! Arrays for soil water freezing/melting
        zsandf(i)   = csandf(mstyp)
        zclayf(i)   = cclayf(mstyp)
        zpsis(i)    = -zpsi0 * 10._ireals**(1.88_ireals-0.013_ireals*zsandf(i))
        zb_por(i)   = 2.91_ireals + .159_ireals*zclayf(i)
        zedb(i)     = 1._ireals/zb_por(i)

! print*,i,zporv(i),zclayf(i),zpsis(i),zb_por(i)

!<JH
!DR        IF(m_styp(i) < 9 ) llandmask(i) = .true.
!>JH      ENDIF

    ENDDO


  ! Set three-dimensional variables
  DO kso = 1, ke_soil
      DO i = istarts, iends
!        IF(llandmask(i)) THEN        ! for land-points only
          mstyp           = soiltyp_subs(i)        ! soil type
          zalam(i,kso)  = cala0(mstyp)              ! heat conductivity parameter
!        ENDIF
      ENDDO
  ENDDO


! For ntstep=nstart : Some preparations
! =====================================

!>JH  IF (ntstep == nstart) THEN
!   Determine constants clgk0 for BATS-scheme
  DO jb       = 1, 10
    clgk0(jb) = LOG10(MAX(zepsi,ck0di(jb)/ckrdi))
  END DO
!<JH  ENDIF

! To ensure t_s = t_so(1) also at gme2lm-modified lateral boundaries by
! modifying the soil temperature profile:
  DO kso   = 1,ke_soil+1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN             ! for land-points only
        ! Next lines have to be modified/deleted if the soil surface temperature
        ! t_so(i,0,nx) predicted by the heat conduction equation is used
          zdelt_bound = t_s_now(i) - t_so_now(i,1)
          t_so_now(i,kso) = t_so_now(i,kso) + zdelt_bound * (1._ireals - &
                           (zmls(kso) - zmls(1))/(zmls(ke_soil+1) - zmls(1)))
          IF(kso.EQ.1) t_so_now(i,0) = t_so_now(i,1)
!        ENDIF
      END DO
  END DO

    DO i = istarts, iends
        ztrangs(i)      = 0.0_ireals
        zsnull(i)       = 0.0_ireals
        zs1(i)          = 0.0_ireals
        zwroot(i)       = 0.0_ireals
        zdwidt (i)      = 0.0_ireals         ! Initialisation of all
        zdwsndt(i)      = 0.0_ireals         ! evaporation quantities
        zesoil (i)      = 0.0_ireals         ! to be equal to zero
        zrr    (i)      = 0.0_ireals         ! in first part formation of dew
        zrs    (i)      = 0.0_ireals         ! in first part formation of rime
        zw_fr(i,ke_soil+1)  = w_so_now(i,ke_soil+1)/zdzhs(ke_soil+1)
    END DO


  DO kso   = 1, ke_soil
      DO i = istarts, iends
        zw_fr(i,kso)    = w_so_now(i,kso)/zdzhs(kso)
        ztrang (i,kso)  = 0._ireals
        zropartw(i,kso) = 0._ireals
      END DO
  END DO

  !>JH WRITE(*,*) 'zw_fr: ',zw_fr(1,1,:)

! Determine the layer indices and total depth for water content averaging
! (in soil evaporation after Dickinson)
  k10cm  = 1
  k100cm = 1
  z1     = zzhls(1)
  znull  = zzhls(1)
  DO kso = 1, ke_soil
    IF (zmls(kso).le.0.1_ireals) THEN
      z1      = zzhls(kso)
      k10cm   = kso
    END IF
    IF (zmls(kso).le.1.0_ireals) THEN
      znull   = zzhls(kso)
      k100cm  = kso
    END IF
  END DO

! Determine average soil water content over 10 and 100 cm, respectively.
  DO kso   = 1, k100cm
      DO i = istarts, iends
!        IF (llandmask(i)) THEN   ! for land-points only
          IF (kso.le.k10cm) zs1(i) = zs1(i) + w_so_now(i,kso)
          zsnull(i)   = zsnull(i) + w_so_now(i,kso)
!        END IF
      END DO
  END DO

! JH Copy soil moisture and ice content for all soiltypes and levels on new state, untouched for rock and ice
 
  w_so_new(istarts:iends,ke_soil+1) = w_so_now(istarts:iends,ke_soil+1)
  w_so_ice_new(istarts:iends,:)     = w_so_ice_now(istarts:iends,:)
  

!<em subs  
!w_snow per fraction -> w_snow per unit
!  DO i = istarts, iends
!    IF (llandmask(i)) THEN   ! for land-points only
!      zf_snow  (i, 1) = MAX( 0.01_ireals, MIN(1.0_ireals,w_snow_now(i, 1)/cf_snow) )* &
!                        zsf_heav(w_snow_now(i, 1) - zepsi)                               !true snow fraction
!     zf_snow(i) = MAX(0.01_ireals, snowfrac(i))*zsf_heav(w_snow_now(i) - zepsi)
!      IF(pt_tiles(ns)%snow_tile .OR. pt_tiles(ns)%snowfree_tile) THEN
!        w_snow_now(i) = w_snow_now(i)*zf_snow(i, 1)
!      END IF
!    END IF
!  END DO
!em>

  IF (itype_heatcond == 2) THEN

! heat conductivity dependent on actual soil water content
! based on Peters-Lidard et al. (1998) and Johansen (1975)
! see also Block, Alexander (2007), Dissertation BTU Cottbus
    
    DO kso = 1, ke_soil+1
      DO i = istarts, iends
!       IF (llandmask(i)) THEN          ! land-points only
        zthetas = zporv(i)
  
        zlamli  = 0.57_ireals
        zlamic  = 2.2_ireals
        zthliq  = zthetas - (w_so_ice_now(i,kso)/zdzhs(kso))
        zlamq   = 7.7_ireals
  
        rsandf=zsandf(i)/100._ireals
  
        if (rsandf >= 0.2_ireals)  zlam0    =  2.0_ireals
        if (rsandf < 0.2_ireals)   zlam0    =  3.0_ireals
  
  
        zlams   = zlamq**rsandf* zlam0**(1._ireals-rsandf)
  
  !      zlamsat = (clamso(m_styp(i))**(1.0_ireals-zthetas)) * (zlamli**zthliq) &
  !                * (zlamic**(zthetas-zthliq))
        zlamsat = (zlams**(1.0_ireals-zthetas)) * (zlamli**zthliq) * (zlamic**(zthetas-zthliq))
  
        zrhod   = 2700.0_ireals*(1.0_ireals-zthetas)
        zlamdry = ( 0.135_ireals*zrhod + 64.7_ireals )                     &
                / ( 2700.0_ireals - 0.947_ireals*zrhod )
  
        zsri    = MIN(1.0_ireals, (w_so_now(i,kso)/zdzhs(kso)) / zthetas)
  
        IF ( zsri >= 0.1_ireals ) THEN
         IF ( t_so_now(i,kso) < t0_melt) THEN
          zKe = zsri
         ELSE
          zKe = LOG10(zsri) + 1.0_ireals
         ENDIF
         zKe = MAX(0.0_ireals, (MIN(1.0_ireals, zKe)) )
         hzalam(i,kso) = zKe*(zlamsat - zlamdry) + zlamdry
        ELSE
         hzalam(i,kso) = zlamdry
        ENDIF
!       END IF  ! land-points
      ENDDO
    ENDDO
  
    DO kso = 1, ke_soil
      DO i = istarts, iends
!       IF (llandmask(i)) THEN          ! land-points only
        ! mean heat conductivity
        zalam(i,kso) = hzalam(i,kso)*hzalam(i,kso+1)*(zmls(kso+1)-zmls(kso))    &
                       / ( hzalam(i,kso)*(zmls(kso+1)-zzhls(kso))                   &
                       +   hzalam(i,kso+1)*(zzhls(kso)-zmls(kso)) )
!       END IF  ! land-points
      ENDDO
    ENDDO
  
  ELSE

! heat conductivity based on assumption of a soil water content which is equal to the
! average between wilting point and field capacity
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          zwqg         = 0.5_ireals*(zfcap(i) + zpwp(i))
          z4wdpv       = 4._ireals*zwqg/zporv(i)
          ! heat conductivity
          zalamtmp(i) =              zdlam(i)                         &
                        * (0.25_ireals + 0.30_ireals*zdlam(i)           &
                        / (1._ireals+0.75_ireals*zdlam(i)))             &
                        * MIN (z4wdpv, 1.0_ireals + (z4wdpv-1.0_ireals)   &
                        *(1.0_ireals+0.35_ireals*zdlam(i))              &
                        /(1.0_ireals+1.95_ireals*zdlam(i)))
!        END IF  ! land-points
      ENDDO
    DO kso = 1, ke_soil
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            zalam(i,kso) = zalam(i,kso) + zalamtmp(i)
!          END IF  ! land-points
        ENDDO
    ENDDO
  
  ENDIF


! Initialisations and conversion of tch to tmch

  limit_tch(:) = .false.  !  preset problem indicator
  ldebug         = .false.

    DO i = istarts, iends
!      IF (llandmask(i)) THEN     ! for land-points only
#ifdef __ICON__
        zuv        = SQRT ( u(i)**2 + v(i)**2 )
#else
        im1        = MAX( 1, i-1)
        zuv        = 0.5_ireals*SQRT ( (u(i) + u(im1))**2 &
                               +(v(i) + v(im1))**2 )
#endif
        ztvs       = t_g (i)*(1.0_ireals + rvd_m_o*qv_s(i))

       
        !  'potential' temperature of lowest atmospheric layer
        zplow          = p0(i) ! + pp(i)
        zth_low (i)  =  t(i) *( (ps(i)/zplow)**rdocp )

        zdt_atm (i)  =  zth_low(i)-t_g(i)
        zdq_atm (i)  =  qv(i)-qv_s(i)

!   introduce an artificical upper boundary on transfer coefficients in cases
!   where an extreme cooling/heating of topmost soil layer may be caused due
!   to excessive sensible/latent heat fluxes (e.g. after creation of unbalanced
!   structures in the data assimilation or following massive changes in radiative
!   forcing due to infrequent radiation computations)

!   estimate current energy budget of topmost soil layer

!   heat flux between layers 1&2 based on current temperature profile
        zg1(i)= zalam(i,1)*(t_so_now(i,1)-t_so_now(i,2))/zdzms(2)

!   estimates of sensible and latent heat flux
        zrho_atm(i)=ps(i)/(r_d*ztvs)
        zshfl(i) = tch(i)*zuv*zrho_atm(i)*cp_d*zdt_atm(i)
        zlhfl(i) = tch(i)*zuv*zrho_atm(i)*lh_v*zdq_atm(i)
    
        IF (ABS(zshfl(i)) <= zepsi) zshfl(i)=SIGN(zepsi,zshfl(i))
    
        zthfl(i) = zshfl(i)+zlhfl(i)
    
!   net radiative fluxes at surface
        zradfl(i) = sobs(i)+thbs(i)

!   unconstrained estimated energy budget of topmost soil layer
        zeb1(i) = zshfl(i)+zlhfl(i)+zradfl(i)-zg1(i)

!   energy required to melt existing snow
        ze_melt(i)=w_snow_now(i)*rho_w*lh_f     ! (J/m**2)
!   heat capacity of snow layer
        zch_snow(i)=w_snow_now(i)*rho_w*chc_i   ! (J/(m**2 K))

!   constrain transfer coefficient, if energy budget  of topmost soil layer is:
!   a) negative & surface layer is unstable (i.e   upward directed turbulent heat flux)
!   b) positive & surface layer is stable   (i.e downward directed turbulent heat flux)

        IF (zeb1(i)<0.0_ireals .AND. zthfl(i)<0.0_ireals) THEN 
    ! cooling of 1st soil layer&upward SHF+LHF

          ztchv_max(i) =  &
               (-zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i))/zdt            &
                +zg1(i)-zradfl(i) )                                     &
               / zthfl(i) * tch(i) * zuv
        ELSEIF (zeb1(i)>0.0_ireals .AND. zthfl(i)>0.0_ireals) THEN 
    ! heating of 1st soil layer & downward SHF+LHF
    !   Note: The heat capacity of snow is only relevant for the difference 
    !         between the actual temperature and the melting point. The mean 
    !         snow temperature is set to the average of t_snow & t_so(1)
          IF(lmulti_snow) THEN
            zdT_snow=MIN(0._ireals, &
              &          t0_melt-0.5_ireals*(t_snow_mult_now(i,1)+t_so_now(i,1)))
          ELSE
            zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow_now(i)+t_so_now(i,1)))
          ENDIF
          ztchv_max(i) =  &
          ( (zlim_dtdt*zch_soil*zdzhs(1)+zdT_snow*zch_snow(i)+ze_melt(i))/zdt  &
            + zg1(i)-zradfl(i) )                                               &
                      / zthfl(i) * tch(i) * zuv
        ELSE                                                    
          ! unlimited transfer coefficient
          ztchv_max(i) = HUGE(1._ireals)
        ENDIF
                                                        ! required constraint as non-turbulent
                                                        ! energy budget components alone may
        ztchv_max(i) =MAX( ztchv_max(i) ,zepsi)     ! exceed the upper limit in the energy
                                                        ! budget of the 1st soil layer

        ztchv(i)    = tch(i)*zuv  ! transfer coefficient * velocity

        LIM: IF (ztchv(i) > ztchv_max(i)) THEN
          tch(i)=ztchv_max(i)/MAX(zuv,1.E-06_ireals)
!         IF (ntstep > 10) THEN          ! control print only after initial adaptation
!                                        ! to analysis state
            limit_tch(i) = .true.      ! set switch for later use
        END IF LIM

        ztmch(i) = tch(i)*zuv*g*zrho_atm(i)
!     ENDIF
    ENDDO

! counter for limitation of transfer coefficients
  m_limit = COUNT( limit_tch(:) ) 


! In debugging mode and if transfer coefficient occured for at least one grid point
  IF (m_limit > 0 .AND. ldebug) THEN
    WRITE(*,'(1X,A,/,2(1X,A,F10.2,A,/),1X,A,F10.2,/,1X,A,F10.3,/)')                  &
           'terra1: transfer coefficient had to be constrained',                     &
           'model time step                                 :', zdt     ,' seconds', &
           'max. temperature increment allowed per time step:',zlim_dtdt,' K',       &
           'upper soil model layer thickness                :', zdzhs(1)

      DO i = istarts, iends
#ifdef __ICON__
      IF (limit_tch(i)) THEN
        zuv        = SQRT (u(i)**2 + v(i)**2 )
#else
      im1 = MAX(1,i-1)
      IF (limit_tch(i)) THEN
        zuv        = 0.5_ireals*SQRT ( (u(i) + u(im1))**2 &
                               +(v(i) + v(im1))**2 )
#endif
        yhc       ='COOLING'
        IF (zeb1(i) > 0._ireals) Yhc='HEATING'
        END IF
      END DO
  ENDIF


! Update indicator for age of snow in top of snow layer

  DO i = istarts, iends
!    IF(llandmask(i)) THEN        ! for land-points only
      IF (w_snow_now(i) <=0.0_ireals) THEN
        ! if no snow exists, reinitialize age indicator
        freshsnow(i) = 1.0_ireals
      ELSE
        IF ( itype_gscp == 4 ) THEN
          zsnow_rate = prs_gsp(i)+prs_con(i)+prg_gsp(i) ! [kg/m**2 s]
        ELSE
          zsnow_rate = prs_gsp(i)+prs_con(i)              ! [kg/m**2 s]
        ENDIF

        ! linear decay equals 1.0 in 28 days
        zdsn_old   = zdt/86400._ireals/28._ireals

        ! linear growth rate equals 1.0 in 1 day for a constant
        ! snow rate of 5 mmH2O (kg/m**2) per day
        zdsn_new   = zdt*zsnow_rate*0.2_ireals   !

        ! reduce decay rate, if new snow is falling and as function of snow
        ! age itself
        zdsn_old   = (zdsn_old - zdsn_new)*freshsnow(i)
        zdsn_old   = MAX(zdsn_old,0._ireals)

        freshsnow(i) = freshsnow(i) + zdsn_new-zdsn_old

        freshsnow(i) = MIN(1._ireals,MAX(0._ireals,freshsnow(i)))
     END IF
!   END IF
  ENDDO


!------------------------------------------------------------------------------
! Section I.2: temperatures, water contents (in mH2O), surface pressure, 
!------------------------------------------------------------------------------

  IF(lmulti_snow) THEN

    DO ksn = 0,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN     ! for land-points only
            IF (w_snow_now(i) > 0.0_ireals) THEN
              ! existence of snow
              t_snow_mult_now(i,ksn) = MIN (t0_melt - zepsi, t_snow_mult_now(i,ksn) )
            ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
              ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
              t_snow_mult_now(i,ksn) = MAX (t0_melt + zepsi, t_s_now(i) )
            ELSE
              ! no snow and  t_snow < t0_melt
              ! --> t_snow = t_s
              t_snow_mult_now(i,ksn) = MIN (t0_melt - zepsi, t_s_now(i) )
            END IF
!          END IF
        ENDDO
    ENDDO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN     ! for land-points only
        IF (w_snow_now(i) > 0.0_ireals) THEN
          ! existence of snow
          ! --> no water in interception store and t_snow < t0_melt
          w_snow_now(i) = w_snow_now(i) + w_i_now(i)
          wtot_snow_now(i,1) = wtot_snow_now(i,1) + w_i_now(i)
          dzh_snow_now(i,1)  = dzh_snow_now(i,1)  + w_i_now(i) &
            &                      /rho_snow_mult_now(i,1)*rho_w
          w_i_now(i) = 0.0_ireals
          h_snow_now(i) = h_snow(i)
        ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
          ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
          t_s_now   (i) = t_snow_mult_now(i,ke_snow) 
          h_snow_now(i) = 0.0_ireals
        ELSE
          ! no snow and  t_snow < t0_melt
          ! --> t_snow = t_s and no water w_i in interception store
          t_s_now   (i) = t_snow_mult_now(i,ke_snow)
          w_i_now(i) = 0.0_ireals
          h_snow_now(i) = 0.0_ireals
        END IF
!      END IF
    ENDDO
  ELSE ! no multi-layer snow

    DO i = istarts, iends
!      IF (llandmask(i)) THEN     ! for land-points only
        IF (w_snow_now(i) > 0.0_ireals) THEN
          ! existence of snow
          ! --> no water in interception store and t_snow < t0_melt
          w_snow_now(i) = w_snow_now(i) + w_i_now(i)
          w_i_now(i) = 0.0_ireals
          t_snow_now(i) = MIN (t0_melt - zepsi, t_snow_now(i) )
        ELSE IF (t_snow_now(i) >= t0_melt) THEN
          ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
          t_s_now   (i) = MAX (t0_melt + zepsi, t_s_now(i) )
          t_snow_now(i) = t_s_now(i)
        ELSE
          ! no snow and  t_snow < t0_melt
          ! --> t_snow = t_s and no water w_i in interception store
          t_s_now   (i) = MIN (t0_melt - zepsi, t_s_now(i) )
          t_snow_now(i) = t_s_now(i)
          w_i_now(i) = 0.0_ireals
        END IF
!      END IF
    ENDDO

  ENDIF

  !Inizializations for the next sections
  IF (lmulti_snow) THEN
    ! some preparations for ksn==0 and ksn==1
    DO i = istarts, iends
!      IF (llandmask(i)) THEN   ! for land-points only
        ztsnow_mult   (i,0) = t_snow_mult_now(i,0)
        ztsnow_mult   (i,1) = t_snow_mult_now(i,1)

        zhh_snow(i,1) =  -h_snow_now(i) + dzh_snow_now(i,1)
        zhm_snow(i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._ireals

        zdzh_snow(i,1) = dzh_snow_now(i,1)
        zextinct (i,1) = 0.13_ireals*rho_snow_mult_now(i,1)+3.4_ireals

        ! set ztsnow to ztsnow_mult(ksn=1)
        ztsnow   (i) = t_snow_mult_now(i,1)

        ! initialize zdz_snow (from Section I.3) to zdzh_snow(i,1)
        zdz_snow (i) = dzh_snow_now(i,1) ! zdzh_snow
        zwsnow   (i) = w_snow_now(i)
!      END IF
    ENDDO

    DO  ksn = 2,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN   ! for land-points only
            ztsnow_mult(i,ksn) = t_snow_mult_now(i,ksn)

            zhh_snow   (i,ksn) = zhh_snow(i,ksn-1) + dzh_snow_now(i,ksn)
            zhm_snow   (i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals

            zdzh_snow  (i,ksn) = dzh_snow_now(i,ksn)
            zextinct   (i,ksn) = 0.13_ireals*rho_snow_mult_now(i,ksn)+3.4_ireals

            ! build sum over all layers for zdzh_snow in zdz_snow
            zdz_snow   (i)     = zdz_snow(i) + zdzh_snow(i,ksn)
 !         END IF
        ENDDO
    ENDDO
  ELSE

    ! set ztsnow to t_snow
    DO i = istarts, iends
!      IF (llandmask(i)) THEN   ! for land-points only
        ztsnow   (i) = t_snow_now(i)
        zwsnow   (i) = w_snow_now(i)
        zdz_snow (i) = zwsnow(i)*rho_w/rho_snow_now(i)
!      END IF
    ENDDO
  ENDIF

  DO i = istarts, iends
!    IF (llandmask(i)) THEN   ! for land-points only

    ! ztsnow is now set according to lmulti_snow (see above);
    ! therefore we need not take care of it in this loop
    ! ztsnow   (i) = t_snow(i,nx)
    zts      (i) = t_s_now   (i)
    zts_pm   (i) = zsf_heav(zts   (i) - t0_melt)
    ztsnow_pm(i) = zsf_heav(ztsnow(i) - t0_melt)
    zwin     (i) = w_i_now(i)

    ! moisture and potential temperature of lowest atmospheric layer
    zplow          = p0(i) ! + pp(i)
    zqvlow         = qv(i)
    zth_low (i)  =  t(i) *( (ps(i)/zplow)**rdocp )

    ! density*transfer coefficient*wind velocity
    zrhoch(i)    = ztmch(i)*(1._ireals/g) + zepsi

    ! saturation specific humidity for t_s and t_snow and first derivative
    z2iw        = zts_pm(i)*b2w + (1._ireals - zts_pm(i))*b2i
    z4iw        = zts_pm(i)*b4w + (1._ireals - zts_pm(i))*b4i
    z234iw      = z2iw*(b3 - z4iw)
    zqs         = zsf_qsat( zsf_psat_iw(zts(i), z2iw,z4iw), ps(i) )
    zdqs        = zqvlow - zqs
    IF (ABS(zdqs).LT.zepsi) zdqs = 0._ireals
    z2iw        = ztsnow_pm(i)*b2w + (1._ireals - ztsnow_pm(i))*b2i
    z4iw        = ztsnow_pm(i)*b4w + (1._ireals - ztsnow_pm(i))*b4i
    z234iw      = z2iw*(b3 - z4iw)
    zqsnow      = zsf_qsat(zsf_psat_iw(ztsnow(i),z2iw,z4iw), ps(i))
    zdqvtsnow(i)= zsf_dqvdt_iw(ztsnow(i), zqsnow, z4iw,z234iw)
    zdqsnow     = zqvlow - zqsnow
    IF (ABS(zdqsnow).LT.zepsi) zdqsnow = 0._ireals

    ! potential evaporation at T_snow and Ts
    zep_snow(i) = (1._ireals-ztsnow_pm(i))*                       &
                                      tfv(i)*zrhoch(i)*zdqsnow
    zep_s   (i) =                   tfv(i)*zrhoch(i)*zdqs
!    END IF
  ENDDO



!------------------------------------------------------------------------------
! Section I.3: heat conductivity, frozen fraction, snow and
!            water covered fraction, snow density and  height,
!            volumetric heat content
!------------------------------------------------------------------------------

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! snow and water covered fraction
!em        zrss = MAX( 0.01_ireals, MIN(1.0_ireals,zwsnow(i)/cf_snow) ) 
        zrww = MAX( 0.01_ireals, 1.0_ireals -                           &
                                 EXP(MAX( -5.0_ireals, - zwin(i)/cf_w) ) )
!em        zf_snow(i) = zrss*zsf_heav(zwsnow(i) - zepsi)
        zf_wi  (i) = zrww*zsf_heav(zwin  (i) - zepsi)

! BR 7/2005 prognostic  snow density
! US DMironov: for the FLake Model, h_snow has to be a prognostic variable but
!              here it is used only diagnostically. So the values are put
!              to every time level
        ! zdz_snow has been computed in the initializations above depending
        ! on lmulti_snow
        ! zdz_snow(i)=zwsnow(i)*rho_w/rho_snow(i,nx)
        h_snow_new(i) = zdz_snow(i)
        h_snow_now(i) = zdz_snow(i)
!        IF (.NOT. l2tls) h_snow(i,nold) = zdz_snow(i)

!       constrain snow depth and consider this constraint for the computation
!       of average snow density of snow layer
        zdz_snow (i) =  zdz_snow (i)/MAX(0.01_ireals,zf_snow(i))
        zdz_snow (i) =  MAX(cdsmin,zdz_snow(i))

!       limitation of snow depth to 1.5m for snow cover heat transfer
        zdz_snow_fl(i) = MIN(1.5_iREALs, zdz_snow(i))
        IF(.NOT. lmulti_snow) &
          zrocs(i) = chc_i*zdz_snow_fl(i)*rho_snow_now(i)
!BR 7/2005 End
!      END IF
    ENDDO



!------------------------------------------------------------------------------
! Section I.4: Hydrology, 1.Section
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section I.4.1: Evaporation from interception store and from snow cover,
  
  !----------------------------------------------------------------------------
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!

    DO i = istarts, iends
!      IF (llandmask(i)) THEN             ! land points only
        ! Evaporation from interception store if it contains water (wi>0) and
        ! if zep_s<0 indicates potential evaporation for temperature Ts
        ! amount of water evaporated is limited to total content of store
        zdwidt(i) = zsf_heav(-zep_s(i))      &
                      * MAX(-zrhwddt*zwin(i), zf_wi(i)*zep_s(i))
        ! Evaporation of snow, if snow exists (wsnow>0) and if zep_snow<0
        ! indicates potential evaporation for temperature t_snow
        zdwsndt(i) = zsf_heav(-zep_snow(i))  &
                       * MAX(-zrhwddt*zwsnow(i), zf_snow(i)*zep_snow(i))
        ! Formation of dew or rime, if zep_s > 0 . distinction between
        ! dew or rime is only controlled by sign of surface temperature
        ! and not effected by presence of snow !
        zrr(i)=zsf_heav(zep_s   (i))*zep_s   (i)*            zts_pm(i)
        zrs(i)=zsf_heav(zep_snow(i))*zep_snow(i)*(1.0_ireals-zts_pm(i))
!      END IF
    ENDDO
  
  
  !----------------------------------------------------------------------------
  ! Section I.4.2b: Bare soil evaporation, BATS version
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.2) THEN
    ! Calculation of bare soil evaporation after Dickinson (1984)
    ! Determination of mean water content relative to volume of voids
      DO i = istarts, iends
!        IF (llandmask(i)) THEN       ! land points only
          IF (zep_s(i) < 0.0_ireals) THEN   ! upwards directed potential
                                              ! evaporation
            zsnull(i) = zsnull(i)/(znull*zporv(i))
            ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
            zice   = zsf_heav(1.5_ireals - REAL(m_styp(i),ireals)) ! 1 only for ice
            zevap  = zrock(i) + zice                               ! 1 for all soil types
                                                                     ! but rock and ice (=0)
            zbeta  = 0.0_ireals
            IF (m_styp(i).ge.3) THEN ! Computations not for ice and rocks
              ! auxiliary quantities
              zbf1   = 5.5_ireals - 0.8_ireals* zbedi(i)*                &
                      (1.0_ireals + 0.1_ireals*(zbedi(i) - 4.0_ireals)*  &
                       clgk0(m_styp(i)) )
              zbf2   = (zbedi(i) - 3.7_ireals + 5.0_ireals/zbedi(i))/  &
                      (5.0_ireals + zbedi(i))
              zdmax  = zbedi(i)*cfinull*zk0di(i)/crhowm 
              zs1(i)  = zs1(i)/(z1*zporv(i))
              zd     = 1.02_ireals*zdmax*zs1(i)**(zbedi(i) + 2._ireals) * &
                                                (zsnull(i)/zs1(i))**zbf1
              zck    = (1.0_ireals + 1550.0_ireals*cdmin/zdmax)*zbf2
              ! maximum sustainable moisture flux in the uppermost surface
              ! layer in kg/(s*m**2)
              zfqmax = - rho_w*zck*zd*zsnull(i)/SQRT(znull*z1)
              zevapor= MAX(zep_s(i),zfqmax)
              IF(zw_fr(i,1)+zevapor*(1.0_ireals - zf_wi(i))   &
                                     *(1.0_ireals - zf_snow(i)) &
                                     *eai(i)/sai(i)* zdtdrhw/zdzhs(1) &
                                     .LE.zadp(i)) zevapor = 0._ireals
              zbeta  = zevapor/MIN(zep_s(i),-zepsi)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_ireals - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover
            zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
                          *(1.0_ireals - zf_wi  (i)) & ! not water covered
                          *(1.0_ireals - zf_snow(i)) & ! not snow covered
                          * eai(i)/sai(i) ! relative source surface
                                              ! of the bare soil
!            lhfl_bs(i) = lh_v * zesoil(i)
          END IF  ! upwards directed potential evaporation
!        END IF    ! land points
      END DO
  END IF ! BATS version


  !----------------------------------------------------------------------------
  ! Section I.4.2b: Bare soil evaporation, Noilhan and Platon, 1989
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.3) THEN
      DO i = istarts, iends
!        IF (llandmask(i)) THEN       ! land points only
          IF (zep_s(i) < 0.0_ireals) THEN   ! upwards directed potential
                                              ! evaporation
            zsnull(i) = zsnull(i)/(znull*zporv(i))
            ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
            zice   = zsf_heav(1.5_ireals - REAL(m_styp(i),ireals)) ! 1 only for ice
            zevap  = zrock(i) + zice                  ! 1 for all soil types
                                                        ! but rock and ice (=0)
            zbeta  = 0.0_ireals
            IF (m_styp(i).ge.3) THEN ! Computations not for ice and rocks

               IF (zw_fr(i,1)> zfcap(i)) THEN
                  zalpha = 1.0_ireals
               ELSE
                  zalpha = 0.5_ireals * (1 - COS ( 0.5_ireals * pi * &
                           (zw_fr(i,1) - zadp(i)) / ( zfcap(i) - zadp(i)) ) )
               ENDIF
               z2iw   = ztsnow_pm(i)*b2w + (1._ireals - ztsnow_pm(i))*b2i
               z4iw   = ztsnow_pm(i)*b4w + (1._ireals - ztsnow_pm(i))*b4i
               zqs    = zsf_qsat( zsf_psat_iw(zts(i), z2iw,z4iw), ps(i) )
               zevapor= MIN(0.0_ireals,zrhoch(i)*(qv(i)-zalpha*zqs))
               
               zbeta  = zevapor/MIN(zep_s(i),-zepsi)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_ireals - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover
            zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
                          *(1.0_ireals - zf_wi  (i)) & ! not water covered
                          *(1.0_ireals - zf_snow(i)) & ! not snow covered
                          * eai(i)/sai(i) ! relative source surface of the bare soil
            
          END IF  ! upwards directed potential evaporation
!        END IF    ! land points
      END DO
  END IF ! NP89


  !----------------------------------------------------------------------------
  ! Section I.4.3b: transpiration by plants, BATS version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.2) THEN   ! BATS version
    ! This version is based on Dickinson's (1984) BATS scheme, simplified by
    ! neglecting the water and energy transports between the soil and the plant
    ! canopy. This leads to a Monteith combination formula for the computation
    ! of plant transpiration.


  ! Root distribution

    IF (itype_root == 2) THEN
        DO i = istarts, iends
          zrootdz_int (i)= 0.0_ireals   !initialize the root density profile integral
          zwrootdz_int(i)= 0.0_ireals   !initialize the root water content integral
        END DO
      DO kso   = 1,ke_soil
          DO i = istarts, iends
!            IF (llandmask(i)) THEN ! land points only,
              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i) < 0.0_ireals) THEN  ! upwards directed potential evaporation
                  ! consider the effect of root depth & root density
                  zrootfr = EXP (-zroota(i)*zmls(kso)) ! root density
                  zrootdz = zrootfr*MIN(zdzhs(kso),MAX(0.0_ireals, zbwt(i)-(zmls(kso) &
                    &       -0.5_ireals*zdzhs(kso)) ) )
                  zrootdz_int(i)=zrootdz_int(i) + zrootdz
                  zwrootdz(i,kso)=zrootdz*(zw_fr(i,kso)-w_so_ice_now(i,kso)/zdzhs(kso))
                  zwrootdz_int(i)=zwrootdz_int(i) + zwrootdz(i,kso)
                END IF  ! negative potential evaporation only
              END IF  ! neither ice or rocks
!            END IF    ! land-points only
          END DO
      END DO


       ! Compute root zone integrated average of liquid water content

          DO i = istarts, iends
             zwrootdz_int(i)=zwrootdz_int(i)/MAX(zrootdz_int(i),zepsi)
          END DO

    ELSE

      DO kso   = 1,ke_soil
          DO i = istarts, iends
!            IF (llandmask(i)) THEN ! land points only,
              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i) < 0.0_ireals) THEN  ! upwards directed potential
                                                   ! evaporation
                  zropart  = MIN ( zdzhs(kso), MAX(0.0_ireals,                     &
                             zbwt(i) - (zmls(kso) - 0.5_ireals*zdzhs(kso))))
                  zropartw(i,kso) = zropart*(zw_fr(i,kso)-w_so_ice_now(i,kso)/zdzhs(kso))
                  zwroot(i) = zwroot(i) + zropartw(i,kso)/MAX(zepsi,zbwt(i))
                END IF  ! negative potential evaporation only
              END IF  ! neither ice or rocks
!            END IF    ! land-points only
          END DO
      END DO

    ENDIF


    ! Determination of the transfer functions CA, CF, and CV

      DO i = istarts, iends
!        IF (llandmask(i)) THEN ! land points only,
          IF (m_styp(i).ge.3) THEN ! neither ice or rocks
            IF (zep_s(i) < 0.0_ireals) THEN  ! upwards directed potential evaporation
              zuv        = SQRT (u_10m(i) **2 + v_10m(i)**2 )
              zcatm      = tch(i)*zuv           ! Function CA

              IF(itype_tran == 1) THEN
                zustar     = zuv*SQRT(tcm(i))
                zrla       = 1.0_ireals/MAX(cdash*SQRT(zustar),zepsi)
              ELSE
                zrla       = 0._ireals
              ENDIF

           
              ! to compute CV, first the stomatal resistance has to be determined
              ! this requires the determination of the F-functions:
              ! Radiation function
              zpar       = pabs(i)  !  PAR
              zf_rad(i)= MAX(0.0_ireals,MIN(1.0_ireals,zpar/cparcrit))
              ztlpmwp(i) = (zfcap(i) - zpwp(i))*(0.81_ireals +       &
                      0.121_ireals*ATAN(-86400._ireals*zep_s(i) - 4.75_ireals))

              ! Soil water function
              IF (itype_root == 2) THEN
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwrootdz_int(i) -  &
                                                zpwp(i))/ztlpmwp(i)))
              ELSE
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwroot(i) -  &
                                                zpwp(i))/ztlpmwp(i)))
              ENDIF

              ! Temperature function
!              IF (ntstep .EQ. 0 .AND. itype_tran .NE. 2) THEN
!                t_2m(i)=t(i)
!              ENDIF
              zf_tem     = MAX(0.0_ireals,MIN(1.0_ireals,4.0_ireals*     &
                           (t_2m(i)-t0_melt)*(ctend-t_2m(i))/(ctend-t0_melt)**2))

              ! Saturation deficit function (function not used, but computations
              ! necessary for determination of  slope of the saturation curve)
              z2iw       = zts_pm(i)*b2w + (1._ireals - zts_pm(i))*b2i
              z4iw       = zts_pm(i)*b4w + (1._ireals - zts_pm(i))*b4i
              zepsat     = zsf_psat_iw(t(i),z2iw,z4iw)
              zepke      = qv(i)*ps(i)/                   &
                                        (rdv + o_m_rdv*qv(i))
              zf_sat     = MAX(0.0_ireals,MIN(1.0_ireals,1.0_ireals -  &
                                              (zepsat - zepke)/csatdef))
              ! zf_sat paralysed:
              zf_sat     = 1.0_ireals

              IF (lstomata) THEN
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/MAX(40._ireals,rsmin2d(i)) -     &
                             1.0_ireals/crsmax)*zf_rad(i)*zf_wat*zf_tem*zf_sat
              ELSE
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/crsmin -     &
                             1.0_ireals/crsmax)*zf_rad(i)*zf_wat*zf_tem*zf_sat
              END IF

              zrstom     = 1.0_ireals/zedrstom              ! stomatal resistance
!              rstom(i) = zrstom
              zrveg      = zrla + zrstom
              ! Transpiration rate of dry leaves:
              ztraleav(i)=zep_s(i)*tai(i)/(sai(i)+zrveg*zcatm)
            END IF  ! upwards directed potential evaporation only
          END IF    ! m_styp > 2
!        END IF      ! land points
      END DO

    ! Consideration of water and snow coverage, distribution to the different
    ! soil layers

    DO     kso       = 1,ke_soil
        DO i         = istarts, iends
!          IF (llandmask(i)) THEN ! land points only,
            IF (m_styp(i).ge.3) THEN ! neither ice or rocks
              IF (zep_s(i) < 0.0_ireals) THEN    ! upwards potential evaporation
                ztrabpf  = ztraleav(i)*                   & ! plant covered part
                           (1.0_ireals - zf_wi(i))*       & ! not water covered
                           (1.0_ireals - zf_snow(i))        ! not snow covered

                ! for root distribution
                IF (itype_root == 2) THEN
                  ztrfr    = zwrootdz(i,kso)/(zrootdz_int(i)*zwrootdz_int(i))
                  ztrang(i,kso) = ztrabpf*ztrfr
                ELSE
                  zrootfc = zropartw(i,kso)/(zwroot(i) + zepsi)
                  ztrang(i,kso) = ztrabpf*zrootfc/MAX(zepsi,zbwt(i))
                  IF(zw_fr(i,kso)+ztrang(i,kso)*zdtdrhw/zdzhs(kso) &
                                    .LT.zpwp(i)) ztrang(i,kso) = 0._ireals
                ENDIF
!                lhfl_pl(i,kso)= lh_v * ztrang(i,kso)
                ztrangs(i)    = ztrangs(i) + ztrang(i,kso)
              END IF  ! upwards directed potential evaporation only
            END IF    ! m_styp > 2
!          END IF      ! land points
        END DO
    END DO          ! loop over soil layers

  END IF ! BATS version


  !----------------------------------------------------------------------------
  ! Section I.4.4: total evapotranspiration and 
  !              associated ficticious soil humidity qv_s
  !----------------------------------------------------------------------------

    DO i = istarts, iends
!      IF (llandmask(i)) THEN   ! land-points only
        ze_sum = zdwsndt(i  )  & ! evaporation of snow
               + zdwidt (i  )  & ! evaporation from interception store
               + zesoil (i  )  & ! evaporation from bare soil
               + ztrangs(i  )  & ! transpiration from all soil layers
               + zrr    (i  )  & ! formation of dew
               + zrs    (i  )    ! formation of rime
        qv_s(i) = qv (i) - ze_sum /(zrhoch(i) + zepsi)
!JH     qv_s(i,nnew) = qv_s(i,nx)

!     END IF     ! land points
    END DO

!------------------------------------------------------------------------------
! End of former module procedure terra1_multlay
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
! Former SUBROUTINE terra2_multlay (yerror, ierror)
!------------------------------------------------------------------------------

!   In the prognostic part II the equation of heat conduction and water
!   transport is solved for a multi-layer soil using the same vertical grid
!   Freezing/melting of soil water/ice is accounted for (optionally). A
!   simple one-layer snow model provides the snow surface temperature and
!   the snow water equivalent.


!------------------------------------------------------------------------------
! Section II.1: Initializations
!------------------------------------------------------------------------------

  ! Computation of derived constants
  z1d2dt = 1._ireals/zdt      ! 1./2*timestep

  ! Number of soil layers contributing to surface run-off
  msr_off  = 0

! Note:
!
! - Introduction of i250, i.e. index of last soil layer completely above
!   2.5m soil depth
! - Usage of i250 to compute subsoil runoff (runoff_g) as drainage flux
!   through bottom of this layer (as suggested by the Rhone-Aggregation
!   Experiment)and to control switching off of soil moisture gradient
!   related flux below (i.e. only sedimentation flux
!   allowed between i250 and i250+1)
  ! no. of layers above 2.5 m soil depth
  i250 = 0
  DO kso=1,ke_soil+1
    IF (zzhls(kso)<=2.5_ireals) i250=kso
    ke_soil_hy = i250
  END DO

  zdtdrhw = zdt/rho_w   ! timestep/density of liquid water

  ! time constant for infiltration of water from interception store
  ! must not be less than 2*time step
  ctau_i  = MAX( ctau_i, zdt )

  ! Utility variable to avoid IF-constructs
  DO kso     = 1,ke_soil
    znlgw1f(kso)  = 0.0_ireals
  END DO
  znlgw1f(1)         = 1.0_ireals

! Initialisations
  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN                 ! land-points
            zdtsnowdt_mult(i,ksn)  = 0.0_ireals
            zdtsdt(i)     = 0.0_ireals
!          END IF
        END DO
    END DO
  ELSE
      DO i = istarts, iends
!        IF (llandmask(i)) THEN                 ! land-points
          zdtsnowdt(i)  = 0.0_ireals
          zdtsdt(i)     = 0.0_ireals
!        END IF
      END DO
  END IF



!------------------------------------------------------------------------------
! Section II.2: Prepare basic surface properties and create some local
!               arrays of surface related quantities (for land-points only)
!               Initialise some fields
!------------------------------------------------------------------------------

    DO i = istarts, iends
!      IF (llandmask(i)) THEN                 ! land-points
        mstyp        = soiltyp_subs(i)
        m_styp(i)  = mstyp
        zsandf(i)  = csandf(mstyp)
        zclayf(i)  = cclayf(mstyp)
        zgsb   (i) = 0.0_ireals
        ztrangs(i) = 0.0_ireals
!      END IF
    END DO
  
  DO   kso = 1,ke_soil+1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN     ! land-points only
          ziw_fr(i,kso) = w_so_ice_now(i,kso)/zdzhs(kso)   ! ice frac.
          zlw_fr(i,kso) = zw_fr(i,kso) - ziw_fr(i,kso)  ! liquid water frac.
          zroc(i,kso)   = zrocg(i) + rho_w*zlw_fr(i,kso)*chc_w +          &
                                         rho_w*ziw_fr(i,kso)*chc_i
          IF (kso<=ke_soil) THEN
            ztrangs(i) = ztrangs(i) + ztrang(i,kso)
            zdwgdt(i,kso) = 0.0_ireals
          END IF
          zflmg(i,kso)  = 0.0_ireals
          zrunoff_grav(i,kso)  = 0.0_ireals
!        END IF
      END DO
  END DO      !soil layers



!------------------------------------------------------------------------------
! Section II.3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------

    DO i = istarts, iends
!      IF(llandmask(i))THEN     ! land-points only

        ! store forcing terms due to evapotranspiration, formation of dew
        ! and rime for later use
        zverbo(i) = zdwidt(i) + zesoil(i) + ztrangs(i) +              &
                        (1._ireals-zf_snow(i))*(zrr(i) + zrs(i))

        zversn(i) = zdwsndt(i) + zrs(i)                                 &
                                  + zsf_heav (zwsnow(i) - zepsi) * zrr(i)

        ! add grid scale and convective precipitation (and graupel, if present)
        ! to dew and rime
        zrr(i) = zrr(i) + prr_con(i) + prr_gsp(i)
        IF ( itype_gscp == 4 ) THEN
          zrs(i) = zrs(i) + prs_con(i) + prs_gsp(i) + prg_gsp(i)
        ELSE
          zrs(i) = zrs(i) + prs_con(i) + prs_gsp(i)
        ENDIF

!em subs
!!$        IF(pt_tiles(ns)%snowfree_tile) THEN
!!$          IF(subsfrac(i,pt_tiles(ns)%conjunct).gt.0._ireals) THEN
!!$            IF((1._ireals-ztsnow_pm(i))*zrs(i) > 0.0_ireals) zrs(i) = 0.0_ireals
!!$            IF((1._ireals-ztsnow_pm(i))*zrr(i) > 0.0_ireals) zrr(i) = 0.0_ireals
!!$          END IF
!!$        END IF 
!em>

        ! infiltration and surface run-off

        ! ice free fraction of first soil layer scaled by pore volume
        ! is used as reduction factor for infiltration rate
        zfr_ice_free     = 1._ireals-ziw_fr(i,1)/zporv(i)

        ! subtract evaporation from interception store to avoid negative
        ! values due to sum of evaporation+infiltration
        zwinstr = zwin(i) + zdwidt(i)*zdtdrhw
        zwinstr = MAX(0.0_ireals,zwinstr)

        ! maximum infiltration rate of the soil (rock/ice/water-exclusion
        zinfmx = zrock(i)*zfr_ice_free*csvoro &
                 *( cik1*MAX(0.5_ireals,plcov(i))*MAX(0.0_ireals,           &
                 zporv(i)-zw_fr(i,1))/zporv(i) + zik2(i) )

        ! to avoid pore volume water excess of the uppermost layer by 
        ! infiltration
        zinfmx = MIN(zinfmx, (zporv(i) - zw_fr(i,1))*zdzhs(1)*zrhwddt)

        ! to avoid infiltration at snow covered parts of soil surface
        zinfmx = zinfmx*(1._ireals - zf_snow(i))

        zwimax = cwimax_ml*(1._ireals + plcov(i)*5._ireals)
        zalf   = SQRT(MAX(0.0_ireals,1.0_ireals - zwinstr/zwimax))

        ! water supply from interception store (if Ts above freezing)
        zinf   = zfr_ice_free*zwinstr*rho_w/ctau_i

        ! possible contribution of rain to infiltration
        IF (zrr(i)-zepsi > 0.0_ireals) THEN
          zalf = MAX( zalf,                                                   &
                 (zrhwddt*MAX(0.0_ireals, zwimax-zwinstr) + zinf)/zrr(i) )
          zalf = MAX( 0.01_ireals, MIN(1.0_ireals, zalf) )

          ! if rain falls onto snow, all rain is considered for infiltration
          ! as no liquid water store is considered in the snowpack
          IF (zwsnow(i) > 0.0_ireals) zalf = 0.0_ireals
        END IF
        ! Increase infiltration and reduce surface runoff (bugfix)
        zalf = 0.0_ireals
        ! rain freezes on the snow surface
        IF (lmulti_snow .AND. zwsnow(i) > 0.0_ireals) zalf = 1.0_ireals
        ! add rain contribution to water supply for infiltration
        zvers = zinf + (1._ireals - zalf)*zrr(i)
        ! final infiltration rate limited by maximum value
        zinfil(i) = MIN(zinfmx,zvers)

        ! surface run-off (residual of potential minus actual infiltration)
        zro_inf       = zvers - zinfil(i)
        runoff_s(i) = runoff_s(i) + zro_inf*zroffdt

        ! change of snow water and interception water store
        ! (negligible residuals are added to the run-off)

        ! snow store
        zdwsndtt = zrs(i) + zdwsndt(i)
        zwsnstr  = zwsnow(i) + zdwsndtt*zdtdrhw
        zwsnstr  = MAX(0.0_ireals, zwsnstr) ! avoid negative values (security)
        zdwseps  = 0.0_ireals
        IF (zwsnstr > 0.0_ireals .AND. zwsnstr < zepsi) THEN
!         IF (ztsnow_pm(i) > 0.0_ireals) THEN
            zdwseps    = zwsnstr*zrhwddt
            runoff_s(i) = runoff_s(i) + zdwseps*zroffdt
            zdwsndtt   = zdwsndtt - zdwseps
!         END IF
        END IF
        zdwsndt(i) = zdwsndtt

        ! interception store
        zdwidtt  = zalf*zrr(i) + zdwidt(i)-zinf 
        zwinstr  = zwin(i) + zdwidtt*zdtdrhw
        zwinstr  = MAX(0.0_ireals, zwinstr) !avoid negative values (security)
        zdwieps  = 0.0_ireals
        IF (zwinstr > 0.0_ireals .AND. zwinstr < zepsi) THEN
          zdwieps    = zwinstr*zrhwddt
          runoff_s(i)= runoff_s(i) + zdwieps*zroffdt
          zdwidtt    = zdwidtt - zdwieps
          zwinstr    = 0.0_ireals
        END IF
        ! add excess over zwimax to runoff
        zro_wi       = zrhwddt*MAX( 0.0_ireals, zwinstr-zwimax )
        zdwidtt      = zdwidtt - zro_wi
        zdwidt(i)  = zdwidtt
        runoff_s(i)= runoff_s(i) + zro_wi*zroffdt
!      END IF            ! land-points only
    END DO



!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

! uppermost layer, kso = 1
    DO i = istarts, iends
      ! sedimentation and capillary transport in soil
      ! Note: The fractional liquid water content (concentration)  of each layer
      !       is normalized by the ice free fraction of each layer in order to
      !       obtain a representative concentration of liquid water in the
      !       'active' part of each soil layer
      !       Hydraulic diffusivity and conductivity coefficients are multiplied
      !       by a reduction factor depending on the maximum ice fraction of the
      !       adjacent layers in order to avoid the transport of liquid water
      !       in to the frozen part of the adjacent layer
!      IF (llandmask(i)) THEN      ! land-points only
        IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
          zice_fr_kso   = ziw_fr(i,1)
          zice_fr_ksop1 = ziw_fr(i,2)
          zlw_fr_kso    = zlw_fr(i,1)
          zlw_fr_ksop1  = zlw_fr(i,2)
  
          ! compute reduction factor for transport coefficients
          zfr_ice  = max (zice_fr_kso,zice_fr_ksop1)
          zredp05  = 1._ireals-zfr_ice/MAX(zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)
  
          ! interpolated scaled liquid water fraction at layer interface
          zlw_fr_ksop05  = 0.5_ireals*(zdzhs(2)*zlw_fr_kso+zdzhs(1)*zlw_fr_ksop1) &
                                               /zdzms(2)
          zdlw_fr_ksop05 = zredp05*zdw(i)*EXP(zdw1(i)*                        &
                           (zporv(i)-zlw_fr_ksop05)/(zporv(i)-zadp(i)) )
          zklw_fr_ksop05 = zredp05*zkw(i)*EXP(zkw1(i)*                        &
                           (zporv(i)-zlw_fr_ksop05)/(zporv(i)-zadp(i)) )
  
  
          ! coefficients for implicit flux computation
          z1dgam1     = zdt/zdzhs(1) 
          zgam2p05    = zdlw_fr_ksop05/zdzms(2)
          zaga(i,1) = 0._ireals
          zagb(i,1) = 1._ireals+zalfa*zgam2p05*z1dgam1
          zagc(i,1) = -zalfa * zgam2p05*z1dgam1
          zagd(i,1) = zlw_fr(i,1) + zinfil(i)*z1dgam1/rho_w  &
                       -zklw_fr_ksop05*z1dgam1                     &
                       +(1._ireals - zalfa)* zgam2p05*z1dgam1*(zlw_fr_ksop1 - zlw_fr_kso)  &
                       +                     zgam2p05*z1dgam1*(zice_fr_ksop1-zice_fr_kso)
  
          ! explicit part of soil surface water flux:
          zflmg (i,1) = - zinfil(i)! boundary value for soil water transport
        ENDIF
!      ENDIF
    END DO

! inner layers 2 <=kso<=ke_soil_hy-1
  DO   kso =2,ke_soil_hy-1
      DO i = istarts, iends
        ! sedimentation and capillary transport in soil
!        IF (llandmask(i)) THEN      ! land-points only
          IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
            zice_fr_ksom1 = ziw_fr(i,kso-1)
            zice_fr_kso   = ziw_fr(i,kso  )
            zice_fr_ksop1 = ziw_fr(i,kso+1)
            zlw_fr_ksom1  = zlw_fr(i,kso-1)
            zlw_fr_kso    = zlw_fr(i,kso  )
            zlw_fr_ksop1  = zlw_fr(i,kso+1)
            ! interpolated scaled liquid water content at interface to layer
            ! above and below
            zlw_fr_ksom05 = 0.5_ireals*(zdzhs(kso-1)*zlw_fr_kso+   &
                                   zdzhs(kso)*zlw_fr_ksom1)/zdzms(kso)
            zlw_fr_ksop05 = 0.5_ireals*(zdzhs(kso+1)*zlw_fr_kso+   &
                                   zdzhs(kso)*zlw_fr_ksop1)/zdzms(kso+1)
  
            ! compute reduction factor for coefficients
            zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
            zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
            zfr_ice          = max (zice_fr_kso,zice_fr_ksop1)
            zredp05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)
            zdlw_fr_ksom05= zredm05*zdw(i)*EXP( zdw1(i)*   &
                               (zporv(i)-zlw_fr_ksom05)/(zporv(i)-zadp(i)) )
            zdlw_fr_ksop05= zredp05*zdw(i)*EXP( zdw1(i)*   &
                               (zporv(i)-zlw_fr_ksop05)/(zporv(i)-zadp(i)) )
            zklw_fr_ksom05= zredm05*zkw(i)*EXP( zkw1(i)*   &
                               (zporv(i)-zlw_fr_ksom05)/(zporv(i)-zadp(i)) )
            zklw_fr_ksop05= zredp05*zkw(i)*EXP( zkw1(i)*   &
                               (zporv(i)-zlw_fr_ksop05)/(zporv(i)-zadp(i)) )
  
            ! coefficients for implicit flux computation
            z1dgam1 = zdt/zdzhs(kso)
            zgam2m05  = zdlw_fr_ksom05/zdzms(kso)
            zgam2p05  = zdlw_fr_ksop05/zdzms(kso+1)
            zaga (i,kso) = -zalfa*zgam2m05*z1dgam1
            zagc (i,kso) = -zalfa*zgam2p05*z1dgam1
            zagb (i,kso) = 1._ireals +zalfa*(zgam2m05+zgam2p05)*z1dgam1
            zagd (i,kso) = zlw_fr(i,kso)+                               &
                                  z1dgam1*(-zklw_fr_ksop05+zklw_fr_ksom05)+ &
                                  (1._ireals-zalfa)*z1dgam1*                &
                                  (zgam2p05*(zlw_fr_ksop1-zlw_fr_kso  )     &
                                  -zgam2m05*(zlw_fr_kso  -zlw_fr_ksom1)   ) &
                                 +z1dgam1*                                  &
                                  (zgam2p05*(zice_fr_ksop1-zice_fr_kso  )   &
                                  -zgam2m05*(zice_fr_kso-zice_fr_ksom1)   )
  
            !soil water flux, explicit part, for soil water flux investigations
            ! only)
            zflmg(i,kso) = rho_w &
              &              *(zdlw_fr_ksom05*(zlw_fr_kso+zice_fr_kso-zlw_fr_ksom1-zice_fr_ksom1) &
              &              / zdzms(kso) - zklw_fr_ksom05)
  
            IF(kso==ke_soil_hy-1) THEN
              zflmg(i,kso+1)=rho_w &
                &       *(zdlw_fr_ksop05*(zlw_fr_ksop1+zice_fr_ksop1-zlw_fr_kso-zice_fr_kso) &
                &       / zdzms(kso+1) - zklw_fr_ksop05)
            ENDIF
          ENDIF
!        ENDIF
      END DO
  END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN      ! land-points only
        IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
          ! lowest active hydrological layer ke_soil_hy-1
          zice_fr_ksom1 = ziw_fr(i,ke_soil_hy-1)
          zice_fr_kso   = ziw_fr(i,ke_soil_hy  )
          zlw_fr_ksom1  = zlw_fr(i,ke_soil_hy-1)
          zlw_fr_kso    = zlw_fr(i,ke_soil_hy  )
          zlw_fr_ksom05 = 0.5_ireals*(zdzhs(ke_soil_hy-1)*zlw_fr_kso+ &
                              zdzhs(ke_soil_hy)*zlw_fr_ksom1)/zdzms(ke_soil_hy)

          zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
          zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
          zdlw_fr_ksom05= zredm05*zdw(i)*EXP( zdw1(i)* &
                            (zporv(i)-zlw_fr_ksom05)/(zporv(i)-zadp(i)) )

          z1dgam1 = zdt/zdzhs(ke_soil_hy)
          zgam2m05  = zdlw_fr_ksom05/zdzms(ke_soil_hy)
          zklw_fr_ksom05= zredm05*zkw(i)*EXP( zkw1(i)* &
                            (zporv(i)-zlw_fr_ksom05)/(zporv(i)-zadp(i)) )
          zaga(i,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
          zagb(i,ke_soil_hy) = 1._ireals+ zalfa*zgam2m05*z1dgam1
          zagc(i,ke_soil_hy) = 0.0_ireals
          zagd(i,ke_soil_hy) = zlw_fr(i,ke_soil_hy)+z1dgam1*zklw_fr_ksom05 &
                            +(1._ireals-zalfa)*z1dgam1*                              &
                             zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)            &
                            +z1dgam1*                                         &
                             zgam2m05*(zice_fr_ksom1-zice_fr_kso )   
        ENDIF
!      ENDIF
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! generalized upper boundary condition
        IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
          zagc(i,1) = zagc(i,1)/zagb(i,1)
          zagd(i,1) = zagd(i,1)/zagb(i,1)
        ENDIF
!      END IF          ! land-points only
    END DO

  DO kso=2,ke_soil_hy-1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
            zzz = 1._ireals/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
            zagc(i,kso) = zagc(i,kso) * zzz
            zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
          ENDIF
!        END IF          ! land-points only
      END DO
  END DO                ! soil layers

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
           zage(i,ke_soil_hy) = (zagd(i,ke_soil_hy)-zaga(i,ke_soil_hy)*  &
                             zagd(i,ke_soil_hy-1))/                          &
                            (zagb(i,ke_soil_hy) - zaga(i,ke_soil_hy)*      &
                             zagc(i,ke_soil_hy-1))
        ENDIF
!      END IF          ! land-points only
    END DO

  DO kso = ke_soil_hy-1,1,-1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
            zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
            ! compute implicit part of new liquid water content and add existing
            ! ice content
            w_so_new(i,kso) = zage(i,kso)*zdzhs(kso) + w_so_ice_now(i,kso)
          END IF
!        END IF          ! land-points only
      END DO
  END DO                ! soil layers

!lowest active hydrological level
    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
          ! boundary values ensure that the calculation below leaves the climate
          ! layer water contents unchanged compute implicit part of new liquid
          ! water content and add existing ice content
          w_so_new(i,ke_soil_hy) = zage(i,ke_soil_hy)*zdzhs(ke_soil_hy) + &
                          w_so_ice_now(i,ke_soil_hy)
        END IF 
!      END IF          ! land-points only
    END DO

! to ensure vertical constant water concentration profile beginning at 
! layer ke_soil_hy for energetic treatment only
! soil water climate layer(s)

  IF (itype_hydbound == 3) THEN
    ! ground water as lower boundary of soil column
    DO kso = ke_soil_hy+1,ke_soil+1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
              w_so_new(i,kso) = zporv(i)*zdzhs(kso)
            END IF
!          END IF          ! land-points only
        END DO
    END DO
  ELSE
    DO kso = ke_soil_hy+1,ke_soil+1
       DO i = istarts, iends
!         IF (llandmask(i)) THEN          ! land-points only
           IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
             w_so_new(i,kso) = w_so_new(i,kso-1)*zdzhs(kso)/zdzhs(kso-1)
           END IF
!         END IF          ! land-points only
       END DO
    END DO
  ENDIF

! combine implicit part of sedimentation and capillary flux with explicit part
! (for soil water flux investigations only)
  DO kso = 2,ke_soil+1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type
            zice_fr_ksom1 = ziw_fr(i,kso-1)
            zice_fr_kso   = ziw_fr(i,kso)
            zlw_fr_ksom1_new= w_so_new(i,kso-1)/zdzhs(kso-1) - zice_fr_ksom1
            zlw_fr_kso_new  = w_so_new(i,kso  )/zdzhs(kso  ) - zice_fr_kso
            zlw_fr_ksom1  = w_so_now(i,kso-1)/zdzhs(kso-1) - zice_fr_ksom1
            zlw_fr_kso    = w_so_now(i,kso  )/zdzhs(kso  ) - zice_fr_kso
            !... additionally for runoff_g at lower level of lowest active water
            ! layer calculated with (upstream) main level soil water content
            ! compute reduction factor for transport coefficients
            zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
            zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
  
          ! interpolated liquid water content at interface to layer above
            zlw_fr_ksom05 =0.5_ireals*(zdzhs(kso)*zlw_fr_ksom1+zdzhs(kso-1)*zlw_fr_kso) &
                              /zdzms(kso)
            zdlw_fr_ksom05= zredm05*zdw(i)*EXP(zdw1(i) *  &
                            (zporv(i)-zlw_fr_ksom05)/(zporv(i)-zadp(i)) )
            zklw_fr_ksom05= zredm05*zkw(i) * EXP(zkw1(i)* &
                            (zporv(i)-zlw_fr_ksom05)/(zporv(i)-zadp(i)) )
  
            IF (kso> ke_soil_hy) zdlw_fr_ksom05=0.0_ireals   ! no flux gradient 
                                                      ! contribution below 2.5m
            IF (kso> ke_soil_hy) zklw_fr_ksom05=0.0_ireals   ! no gravitation flux below 2.5m
            zflmg(i,kso) =                              &
                     (1._ireals-zalfa) * zflmg(i,kso) + & ! explicit flux component
                                zalfa  * rho_w *          & ! implicit flux component
                   (zdlw_fr_ksom05 * (zlw_fr_kso_new+zice_fr_kso-zlw_fr_ksom1_new-zice_fr_ksom1) &
                     /zdzms(kso) - zklw_fr_ksom05)
            zredm = 1._ireals-zice_fr_kso/(zlw_fr_kso+zice_fr_kso)
            zklw_fr_kso_new = zredm*zkw(i) * EXP(zkw1(i)* &
                              (zporv(i) - zlw_fr_kso_new)/(zporv(i) - zadp(i)) )
  
            ! actual gravitation water flux
            IF (w_so_new(i,kso).LT.1.01_ireals*zadp(i)*zdzhs(kso)) THEN
              zklw_fr_kso_new=0._ireals
            ENDIF
            zrunoff_grav(i,kso) =  - rho_w * zklw_fr_kso_new
  
            ! ground water as lower boundary of soil column
            IF ((kso == ke_soil_hy+1).and.(itype_hydbound == 3)) THEN
               zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )
  
               zdlw_fr_kso = zredm05*zdw(i)*EXP(zdw1(i) *  &
                    (zporv(i)-zlw_fr_kso_new)/(zporv(i)-zadp(i)) )
               zklw_fr_kso = zredm05*zkw(i) * EXP(zkw1(i)* &
                    (zporv(i)-zlw_fr_kso_new)/(zporv(i)-zadp(i)) )
               zklw_fr_ksom1 = zredm05*zkw(i) * EXP(zkw1(i)* &
                    (zporv(i)-zlw_fr_ksom1_new)/(zporv(i)-zadp(i)) )
  
               zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
               zrunoff_grav(i,ke_soil_hy)=zrunoff_grav(i,ke_soil_hy)+ &
                  zdhydcond_dlwfr / &
                  (1.0_ireals-exp(-zdhydcond_dlwfr/zdlw_fr_kso*0.5_ireals*zdzms(ke_soil_hy+1)))* &
                  zdelta_sm
            ENDIF
          END IF 
!        END IF          ! land-points only
      END DO
  END DO


  DO  kso = 1,ke_soil
    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = zsf_heav(0.5_ireals + REAL(msr_off - kso,ireals))  ! 1.0 for 'surface runoff'
    zro_gfak = 1._ireals - zro_sfak                               ! 1.0 for 'ground runoff'

    IF (kso==i250) THEN  
      zfmb_fak = 1.0_ireals
    ELSE
      zfmb_fak = 0.0_ireals
    END IF

      DO i = istarts, iends
!        IF (llandmask(i)) THEN      ! land-points only
          ! sedimentation and capillary transport in soil
          IF (m_styp(i) >= 3) THEN   ! neither ice nor rock as soil type

            ! first runoff calculation without consideration of
            ! evapotranspiration
            !zdwg   =  zflmg(i,kso+1) - zflmg(i,kso)
            !zdwg calculated above by flux divergence has to be aequivalent with 
            zdwg =  (w_so_new(i,kso)/zdzhs(kso)-zw_fr(i,kso))*zdzhs(kso) &
                                                                   /zdtdrhw
            zdwg =  zdwg + zrunoff_grav(i,kso)*zfmb_fak
            zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,kso) -     &
                       zfcap(i))/MAX(zporv(i) - zfcap(i),zepsi)) )
            zredfu = zsf_heav(zdwg)*zredfu
            zro    = zdwg*zredfu
            zdwg   = zdwg*(1._ireals - zredfu)

            ! add evaporation (znlgw1f: first layer only)
            ! and transpiration (for each layer)
            zdwg   = zdwg + znlgw1f(kso) * zesoil(i) + ztrang (i,kso)
            zwgn   = zw_fr(i,kso) + zdtdrhw*zdwg/zdzhs(kso)
            zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zwgn - zporv(i))
            zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zadp(i) - zwgn )
            zdwgdt(i,kso)= zdwg + zkorr - zro2
            zro    = zro      + zro2
            runoff_s(i) = runoff_s(i) + zro*zro_sfak*zroffdt
            runoff_g(i) = runoff_g(i) + zro*zro_gfak*zroffdt
            ! runoff_g reformulation:
            runoff_g(i) = runoff_g(i) - (zrunoff_grav(i,kso) * zfmb_fak &
                                          + zkorr) * zroffdt
          END IF          ! ice/rock-exclusion
!        END IF   ! land-points only
      END DO
  END DO         ! end loop over soil layers



!------------------------------------------------------------------------------
! Section II.5: Soil surface heat flux (thermal forcing)
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)
        zgstr =   sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i))* &
                  zts(i) + zf_snow(i)*ztsnow_mult(i,1) )**4 + thbs(i)
        zthsnw(i) = - sigma*(1._ireals - Ctalb)*ztsnow_mult(i,1)**4 + zgstr
        zthsoi(i) = - sigma*(1._ireals - Ctalb)*zts(i)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux

        zrnet_s(i) = sobs(i) + zthsoi(i)
        zshfl_s(i) = cp_d*zrhoch(i) * (zth_low(i) - zts(i))
        zlhfl_s(i) = (zts_pm(i)*lh_v + (1._ireals-zts_pm(i))*lh_s)*zverbo(i) &
                     / MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
        zsprs  (i) = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i)*zrs(i) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs  (i) = - lh_f*zrs(i)
          zdwidt (i) = zdwidt (i) + zrs(i)
          zdwsndt(i) = zdwsndt(i) - zrs(i)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
          zwimax       = cwimax_ml*(1._ireals + plcov(i)*5._ireals)
          zwinstr      = zwin(i) + zdwidt(i)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zro        = (zwinstr - zwimax)*zrhwddt
            zdwidt(i)= zdwidt(i) - zro
            runoff_s(i)  = runoff_s(i) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF ((1._ireals-ztsnow_pm(i))*zrr(i) > 0.0_ireals) THEN
          zsprs  (i) = lh_f*zrr(i)
          zdwidt (i) = zdwidt (i) - zrr(i)
          zdwsndt(i) = zdwsndt(i) + zrr(i)
        END IF

!       Influence of heatflux through snow on total forcing:
        zdwsndt(i) = zdwsndt(i)*zsf_heav(zdwsndt(i) - zepsi/zdtdrhw)
        zwsnew(i)  = zwsnow(i) + zdwsndt(i)*zdtdrhw
        IF (zwsnew(i).GT.zepsi) THEN

          zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (zth_low(i)-csnow_tmin) &
                                                  /(t0_melt          -csnow_tmin)
          zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))

          IF(zdwsndt(i)-zrs(i)-zrr(i).GT.0.0_ireals) THEN

            wtot_snow_now(i,1) = max(wtot_snow_now(i,1) + zdwsndt(i)*zdtdrhw, &
              &                           0.0_ireals)

            zhm_snow(i,1) = zhm_snow(i,1) - (zdwsndt(i)-zrs(i)- &
                              zrr(i))*zdt/rho_i/2._ireals-  &
              zrs(i)*zdt/zrho_snowf/2._ireals- zrr(i)*zdt/rho_i/2._ireals
            zdzh_snow(i,1) = zdzh_snow(i,1) + (zdwsndt(i)-zrs(i)-zrr(i))*zdt/rho_i +  &
              zrs(i)*zdt/zrho_snowf + zrr(i)*zdt/rho_i

            rho_snow_mult_now(i,1) = max(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), &
              &                              0.0_ireals)
          ELSE

            wtot_snow_now(i,1) = max(wtot_snow_now(i,1) + (zrs(i)+zrr(i))*zdtdrhw, &
              &                          0.0_ireals)

            zhm_snow(i,1)  = zhm_snow(i,1) - zrs(i)*zdt/zrho_snowf/2._ireals- &
                               zrr(i)*zdt/rho_i/2._ireals
            zdzh_snow(i,1) = zdzh_snow(i,1) + zrs(i)*zdt/zrho_snowf + zrr(i)*zdt/rho_i

            IF(wtot_snow_now(i,1) .GT. 0._ireals) THEN
              rho_snow_mult_now(i,1) = max(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), &
                &                              0.0_ireals)

              wtot_snow_now(i,1) = max(wtot_snow_now(i,1) &
                &                      + (zdwsndt(i)-zrs(i)-zrr(i))*zdtdrhw,0.0_ireals)

              zhm_snow(i,1)  = zhm_snow(i,1) - (zdwsndt(i)-zrs(i)-zrr(i)) &
                &                *zdt/rho_snow_mult_now(i,1)/2._ireals
              zdzh_snow(i,1) = zdzh_snow(i,1) + (zdwsndt(i)-zrs(i)-zrr(i)) &
                &                *zdt/rho_snow_mult_now(i,1)
            ELSE
              rho_snow_mult_now(i,1) = 0.0_ireals
              zdzh_snow(i,1) = 0.0_ireals
            END IF

          END IF
        END IF
        h_snow_now(i) = 0.
        sum_weight(i) = 0.0_ireals
!      END IF          ! land-points only
    END DO
  DO ksn = 1,ke_snow  
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (zwsnew(i).GT.zepsi) THEN          
            h_snow_now(i) = h_snow_now(i) + zdzh_snow(i,ksn)
          END IF
!        END IF          ! land-points only
      END DO
  END DO    
              
  DO ksn = ke_snow,1,-1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(zwsnew(i) .GT. zepsi) THEN
            IF(zwsnow(i) .GT. zepsi) THEN
              dz_old(i,ksn) = zdzh_snow(i,ksn)
              z_old(i,ksn) = -sum_weight(i) - zdzh_snow(i,ksn)/2._ireals
              sum_weight(i) = sum_weight(i) + zdzh_snow(i,ksn)
              zhh_snow(i,ksn) = -h_snow_now(i)/ke_snow*(ke_snow-ksn)
            ELSE
              zhh_snow(i,ksn) = -h_snow_now(i)/ke_snow*(ke_snow-ksn)
            END IF
          END IF
!        END IF          ! land-points only
      END DO  
  END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        IF(zwsnew(i) .GT. zepsi) THEN
          IF(zwsnow(i) .GT. zepsi) THEN
            zhm_snow (i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._ireals
            zdzh_snow(i,1) = zhh_snow(i,1) + h_snow_now(i)            !layer thickness betw. half levels of uppermost snow layer
            zdzm_snow(i,1) = zhm_snow(i,1) + h_snow_now(i)            !layer thickness between snow surface and main level of uppermost layer
            IF(dz_old(i,1).ne.0..and.rho_snow_mult_now(i,1).ne.0.) THEN
              wliq_snow_now(i,1) = wliq_snow_now(i,1)/dz_old(i,1)
            END IF
          ELSE
            zhm_snow (i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._ireals
            zdzh_snow(i,1) = zhh_snow(i,1) + h_snow_now(i)
            zdzm_snow(i,1) = zhm_snow(i,1) + h_snow_now(i)
          END IF
        END IF
!      END IF          ! land-points only
    END DO
  DO ksn = 2,ke_snow
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(zwsnew(i) .GT. zepsi) THEN
            IF(zwsnow(i) .GT. zepsi) THEN
              zhm_snow  (i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
              zdzh_snow (i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,ksn ) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
              IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_now(i,ksn).ne.0.) THEN
                wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)/dz_old(i,ksn)
              END IF
            ELSE
              zhm_snow (i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
              zdzh_snow(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,ksn) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
            END IF
          END IF
!        END IF          ! land-points only
      END DO
  END DO

  DO ksn = ke_snow,1,-1
      DO i = istarts, iends
        t_new  (i,ksn) = 0.0_ireals
        rho_new(i,ksn) = 0.0_ireals
        wl_new (i,ksn) = 0.0_ireals
      END DO

    DO k = ke_snow,1,-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(zwsnew(i) .GT. zepsi .AND. zwsnow(i) .GT. zepsi) THEN

              weight = MAX(MIN(z_old(i,k)+dz_old(i,k)/2._ireals,zhm_snow(i,ksn) &
              &+zdzh_snow(i,ksn)/2._ireals)-   &
                       MAX(z_old(i,k)-dz_old(i,k)/2._ireals, &
                       zhm_snow(i,ksn)-zdzh_snow(i,ksn)/2._ireals),0._ireals) &
                       &/zdzh_snow(i,ksn)

              t_new  (i,ksn) = t_new  (i,ksn) + ztsnow_mult  (i,k      )*weight
              rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_now(i,k)*weight
              wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_now(i,k)*weight
            END IF
!          END IF          ! land-points only
        END DO
    END DO
  END DO

  DO ksn = ke_snow,1,-1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(zwsnew(i) .GT. zepsi) THEN
            IF(zwsnow(i) .GT. zepsi) THEN
              ztsnow_mult  (i,ksn      ) = t_new  (i,ksn)
              rho_snow_mult_now(i,ksn) = rho_new(i,ksn)
              wtot_snow_now    (i,ksn) = rho_new(i,ksn)*zdzh_snow(i,ksn)/rho_w
              wliq_snow_now    (i,ksn) = wl_new (i,ksn)*zdzh_snow(i,ksn)
            ELSE
              ztsnow_mult  (i,ksn      ) = t_s_now(i)
              rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,1)
              wtot_snow_now    (i,ksn) = zwsnew(i)/ke_snow
            END IF
          END IF
!        END IF          ! land-points only
      END DO
  END DO

  ! heat conductivity of snow as funtion of water content
  DO ksn = 1, ke_snow
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF (zwsnew(i).GT.zepsi) THEN
            zalas_mult(i,ksn) = 2.22_ireals*(rho_snow_mult_now(i,ksn)/rho_i)**1.88_ireals
          END IF
!        END IF          ! land-points only
      END DO
  END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        IF (zwsnew(i).GT.zepsi) THEN
          zgsb(i) = (zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+zalam(i,1)*zdzms(1))/ &
                      (-zhm_snow(i,ke_snow)+zdzms(1)) * &
                      (ztsnow_mult(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                      +zdzms(1))

        END IF

        ! total forcing for uppermost soil layer
        zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + zsprs(i) ) &
                         * (1._ireals - zf_snow(i)) &
                    + zf_snow(i) * (1._ireals-ztsnow_pm(i)) * zgsb(i)

        IF(zwsnew(i) .GT. zepsi) THEN
          zrnet_snow = sobs(i) * (1.0_ireals - EXP(-zextinct(i,1)*zdzm_snow(i,1))) &
                       + zthsnw(i)
        ELSE
          zrnet_snow = sobs(i) + zthsnw(i)
        END IF
        zshfl_snow(i) = zrhoch(i)*cp_d*(zth_low(i) - ztsnow_mult(i,1))
        zlhfl_snow(i) = lh_s*zversn(i)   
        zfor_snow_mult(i)  = (zrnet_snow + zshfl_snow(i) + zlhfl_snow(i) + zsprs(i))*zf_snow(i)

!      END IF          ! land-points only
    END DO

  ELSE

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)
        zgstr =   sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i))* &
                  zts(i) + zf_snow(i)*ztsnow(i) )**4 + thbs(i)
        zthsnw(i) = - sigma*(1._ireals - Ctalb)*ztsnow(i)**4 + zgstr
        zthsoi(i) = - sigma*(1._ireals - Ctalb)*zts(i)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux

        zrnet_s(i) = sobs(i) + zthsoi(i)
        zshfl_s(i) = cp_d*zrhoch(i) * (zth_low(i) - zts(i))
        zlhfl_s(i) = (zts_pm(i)*lh_v + (1._ireals-zts_pm(i))*lh_s)*zverbo(i) &
                     / MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
        zsprs  (i) = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i)*zrs(i) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs  (i) = - lh_f*zrs(i)
          zdwidt (i) = zdwidt (i) + zrs(i)
          zdwsndt(i) = zdwsndt(i) - zrs(i)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
          zwimax       = cwimax_ml*(1._ireals + plcov(i)*5._ireals)
          zwinstr      = zwin(i) + zdwidt(i)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zro        = (zwinstr - zwimax)*zrhwddt
            zdwidt(i)= zdwidt(i) - zro
            runoff_s(i)  = runoff_s(i) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF (zwsnow(i) == 0.0_ireals .AND.                            &
               (1._ireals-ztsnow_pm(i))*zrr(i) > 0.0_ireals) THEN
          zsprs  (i) = lh_f*zrr(i)
          zdwidt (i) = zdwidt (i) - zrr(i)
          zdwsndt(i) = zdwsndt(i) + zrr(i)
        END IF

!       Influence of heatflux through snow on total forcing:
        zwsnew(i)   = zwsnow(i) + zdwsndt(i)*zdtdrhw
        IF (zwsnew(i).GT.zepsi) THEN
!         heat conductivity of snow as funtion of water content
! BR      zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i)))
!
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
          zalas  = 2.22_ireals*(rho_snow_now(i)/rho_i)**1.88_ireals

! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
!            The water vapour transport associated conductivity is not included.

!        zalas   = 0.023_ireals+(2.290_ireals-0.023_ireals)* &
!                               (7.750E-05_ireals*rho_snow(i,nx) + &
!                                1.105E-06_ireals*prho_snow(i,nx)**2)

          zgsb(i) = zalas*(ztsnow(i) - zts(i))/zdz_snow_fl(i)
        END IF

        ! total forcing for uppermost soil layer
        zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) ) &
                         * (1._ireals - zf_snow(i)) + zsprs(i) &
                    + zf_snow(i) * (1._ireals-ztsnow_pm(i)) * zgsb(i)

!      END IF          ! land-points only
    END DO

 ENDIF ! lmulti_snow



!------------------------------------------------------------------------------
! Section II.6: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

  DO kso = 2,ke_soil
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          ! for heat conductivity: zalam is now 3D
          zakb = zalam(i,kso)/zroc(i,kso)
          zaga(i,kso) = -zalfa*zdt*zakb/(zdzhs(kso)*zdzms(kso))
          zagc(i,kso) = -zalfa*zdt*zakb/(zdzhs(kso)*zdzms(kso+1))
          zagb(i,kso) = 1._ireals - zaga(i,kso) - zagc(i,kso)
          zagd(i,kso) = t_so_now(i,kso) +                                     &
                 (1._ireals - zalfa)*( - zaga(i,kso)/zalfa*t_so_now(i,kso-1)+ &
                 (zaga(i,kso)/zalfa + zagc(i,kso)/zalfa)*t_so_now(i,kso) -  &
                  zagc(i,kso)/zalfa*t_so_now(i,kso+1)  )
!        END IF  ! land-points only
      END DO
  END DO        ! soil layers
  
  
    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! for heat conductivity: zalam is now 3D: here we need layer 1
        zakb = zalam(i,1)/zroc(i,1)
        zaga(i,  1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(1))
        zagc(i,  1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(2))
        zagb(i,  1) = 1._ireals - zaga(i,1) - zagc(i,1)
        zagd(i,  1) = t_so_now(i,1) + (1._ireals - zalfa)* (                 &
                        - zaga(i,1)/zalfa * t_s_now(i) +                     &
                        (zaga(i,1) + zagc(i,1))/zalfa * t_so_now(i,1) -    &
                         zagc(i,1)/zalfa * t_so_now(i,2)   )
        zaga(i,0)    = 0.0_ireals
        zagb(i,0)    = zalfa
        zagc(i,0)    = -zalfa
        zagd(i,0)    = zdzms(1) * zfor_s(i)/zalam(i,1)+(1._ireals-zalfa)* &
                        (t_so_now(i,1) - t_s_now(i))
        zaga(i,ke_soil+1) = 0.0_ireals
        zagb(i,ke_soil+1) = 1.0_ireals
        zagc(i,ke_soil+1) = 0.0_ireals
        zagd(i,ke_soil+1) = t_so_now(i,ke_soil+1)
  
!      END IF          ! land-points only
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        zagc(i,0) = zagc(i,0)/zagb(i,0)
        zagd(i,0) = zagd(i,0)/zagb(i,0)
!      END IF          ! land-points only
    END DO
  
  DO kso=1,ke_soil
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          zzz = 1._ireals/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
          zagc(i,kso) = zagc(i,kso) * zzz
          zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
!        END IF          ! land-points only
      END DO
  END DO                ! soil layers
  
    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        zage(i,ke_soil+1) = (zagd(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                               zagd(i,ke_soil))/                              &
                              (zagb(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                               zagc(i,ke_soil))
!      END IF          ! land-points only
    END DO
  
  DO kso = ke_soil,0,-1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
        ! The surface temperature computed by t_so(i,0,nnew)=zage(i,0) is
        ! presently unused
          t_so_new(i,kso) = zage(i,kso)
!        END IF          ! land-points only
      END DO
  END DO                ! soil layers

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        t_so_new(i,ke_soil+1) = zage(i,ke_soil+1) ! climate value, unchanged
!      END IF          ! land-points only
    END DO


  IF (lmulti_snow) THEN
  
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF(zwsnew(i) .GE. zepsi) THEN
            zrocs(i) = wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w*rho_w + &
              (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i*rho_i
            zakb      = zalas_mult(i,1)/zrocs(i)  !(chc_i*rho_snow(1,nx))
            zaga(i,1) = 0.0_ireals
            zagc(i,1) = -zdt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,2))
            zagb(i,1) = 1._ireals - zagc(i,1)
            zagd(i,1) = ztsnow_mult(i,1) + zdt/zdzh_snow(i,1)*zfor_snow_mult(i)/zrocs(i)
        
            zrocs(i) = wliq_snow_now(i,ke_snow)/wtot_snow_now(i,ke_snow)*chc_w*rho_w +&
              & (wtot_snow_now(i,ke_snow)-wliq_snow_now(i,ke_snow))/                    &
              & wtot_snow_now(i,ke_snow)*chc_i*rho_i
            zakb = zalas_mult(i,ke_snow)/zrocs(i)  !(chc_i*rho_snow(1,nx))
            zaga(i,ke_snow) = -zdt*zakb/(zdzh_snow(i,ke_snow)*zdzm_snow(i,ke_snow))
            zagb(i,ke_snow) = 1.0_ireals - zaga(i,ke_snow)
            zagc(i,ke_snow) = 0.0_ireals
            zagd(i,ke_snow) = ztsnow_mult(i,ke_snow) &
                                - zdt/zdzh_snow(i,ke_snow)*zgsb(i)/zrocs(i)
          END IF
!        END IF          ! land-points only
      END DO 
  
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (zwsnew(i) .GT. zepsi) THEN
            zagc(i,1) = zagc(i,1)/zagb(i,1)
            zagd(i,1) = zagd(i,1)/zagb(i,1)
          END IF
!        END IF          ! land-points only
      END DO
  
    DO ksn = 2, ke_snow-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF(zwsnew(i) .GT. zepsi) THEN
              zrocs(i) = wliq_snow_now(i,ksn)/wtot_snow_now(i,ksn)*chc_w*rho_w + &
                (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))/ &
                &wtot_snow_now(i,ksn)*chc_i*rho_i
              zakb = zalas_mult(i,ksn)/zrocs(i) !(chc_i*rho_snow(ksn,nx))
              zaga(i,ksn) = -zdt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn))
              zagc(i,ksn) = -zdt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn+1))
              zagb(i,ksn) = 1._ireals - zaga(i,ksn) - zagc(i,ksn)
              zagd(i,ksn) = ztsnow_mult(i,ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO                ! snow layers
  
    DO ksn=2,ke_snow-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF(zwsnew(i) .GT. zepsi) THEN
              zzz = 1._ireals/(zagb(i,ksn) - zaga(i,ksn)*zagc(i,ksn-1))
              zagc(i,ksn) = zagc(i,ksn) * zzz
              zagd(i,ksn) = (zagd(i,ksn) - zaga(i,ksn)*zagd(i,ksn-1)) * zzz
            END IF
!          END IF          ! land-points only
        END DO
    END DO                ! snow layers
  
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF(zwsnew(i) .GT. zepsi) THEN
            zage(i,ke_snow) = (zagd(i,ke_snow) - zaga(i,ke_snow)*       &
                                 zagd(i,ke_snow-1))/                        &
                                (zagb(i,ke_snow) - zaga(i,ke_snow)*       &
                                 zagc(i,ke_snow-1))
            ztsnown_mult(i,ke_snow) = zage(i,ke_snow)
          END IF
!        END IF          ! land-points only
      END DO
  
    DO ksn = ke_snow-1,1,-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF(zwsnew(i) .GT. zepsi) THEN
              zage(i,ksn)     = zagd(i,ksn) - zagc(i,ksn)*zage(i,ksn+1)
              ! The surface temperature computed by t_so(0,nnew)=zage(0) is
              ! presently unused
              ztsnown_mult(i,ksn) = zage(i,ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO                ! snow layers
  
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF(zwsnew(i) .GT. zepsi) THEN
            zrocs(i) = wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w*rho_w + &
              (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i*rho_i
            zswitch(i) = MAX(-zfor_snow_mult(i)/50./zrocs(i)*zdt*ke_snow, &
                                zgsb(i)/50./zrocs(i)*zdt*ke_snow)
            zswitch(i) = MAX(zswitch(i),1.E-03_ireals)
  
            IF(zwsnew(i) .LT. zswitch(i)) THEN
  
              ztsnow(i) = (ztsnow_mult(i,1)*zdzh_snow(i,1) &
              &+ ztsnow_mult(i,2)*zdzh_snow(i,2)) / &
                            (zdzh_snow(i,1) + zdzh_snow(i,2))
              ztsn  (i) = t_so_new(i,1)
              tmp_num = ztsnow(i) + zdt*2._ireals*(zfor_snow_mult(i) - zgsb(i))  &
                             /zrocs(i)/(zswitch(i)/rho_snow_mult_now(i,1)*rho_w) &
                             &- ( ztsn(i) - zts(i) )
              zalas  = 2.22_ireals*(rho_snow_mult_now(i,1)/rho_i)**1.88_ireals
  
              ztsnow_im    = - zrhoch(i) * (cp_d + zdqvtsnow(i) * lh_s)       &
                                           - zalas/(zdzh_snow(i,1) + zdzh_snow(i,2))
              zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i)/(zdzh_snow(i,1) &
              &+ zdzh_snow(i,2)))
              tmp_num = ztsnow(i) + (tmp_num-ztsnow(i))/zfak
  
              ztsnown_mult(i,1) = tmp_num
              ztsnown_mult(i,2) = tmp_num
            END IF
  
            ztsnown_mult(i,0) = ztsnown_mult(i,1)
          ELSE
            zswitch(i) = 0.0_ireals
          END IF
!        END IF          ! land-points only
      END DO
  
    DO ksn = 1, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF(zwsnew(i) .GT. zepsi .and. zwsnew(i) .LT. zswitch(i)) THEN
  
              IF(zfor_snow_mult(i)*zdt > zwsnew(i)*rho_w*lh_f) THEN
                ztsnown_mult(i,ksn) = t_so_now(i,0)
              ELSE IF(zfor_snow_mult(i) .GT. 0._ireals) THEN
                ztsnown_mult(i,ksn) = ztsnow_mult(i,ksn) + &
                  zfor_snow_mult(i)*zdt/(chc_i*wtot_snow_now(i,ksn))/rho_w/ke_snow
              END IF
            END IF
!          END IF  ! land-points only
        END DO
    END DO

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF(zwsnew(i) .GT. zepsi .and. zwsnew(i) .LT. zswitch(i)) THEN
  
            IF(zfor_snow_mult(i)*zdt > zwsnew(i)*rho_w*lh_f) THEN
              zfor_snow_mult(i) = zfor_snow_mult(i) - zwsnew(i)*rho_w*lh_f/zdt
              zdwsndt(i) = zdwsndt(i) - zwsnew(i)*rho_w/zdt
              zwsnew(i)  = 0._ireals
              ztsnown_mult(i,0) = t_so_now(i,0)
              zfor_s(i) = zfor_s(i) + zfor_snow_mult(i)
            ELSE IF(zfor_snow_mult(i) .GT. 0._ireals) THEN
              zfor_snow_mult(i) = 0._ireals
              ztsnown_mult(i,0) = ztsnown_mult(i,1)
            END IF
  
          END IF
!        END IF  ! land-points only
      END DO

  END IF


  IF(lmelt) THEN
    IF(.NOT.lmelt_var) THEN
      DO kso = 1,ke_soil
          DO i = istarts, iends
!            IF (llandmask(i)) THEN ! land points only,
              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                ! melting or freezing of soil water
                zenergy   = zroc(i,kso)*zdzhs(kso)*(t_so_new(i,kso) - t0_melt)
                zdwi_max  = - zenergy/(lh_f*rho_w)
                zdelwice  = zdwi_max
                zwso_new  = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
                IF (zdelwice.LT.0.0_ireals)                                    &
                          zdelwice = - MIN( - zdelwice,w_so_ice_now(i,kso))
                IF (zdelwice.GT.0.0_ireals)                                    &
                          zdelwice =   MIN(   zdelwice,zwso_new - w_so_ice_now(i,kso))
                w_so_ice_new(i,kso) = w_so_ice_now(i,kso) + zdelwice

  !             At this point we have 0.0 LE w_so_ice(i,kso,nnew) LE w_so(i,kso,nnew)
  !             If we have 0.0 LT w_so_ice(i,kso,nnew) LT w_so(i,kso,nnew)
  !             the resulting new temperature has to be equal to t0_melt.
  !             If not all energy available can be used to melt/freeze soil water, the
  !             following line corrects the temperature. It also applies to cases
  !             where all soil water is completely ice or water, respectively.

                t_so_new(i,kso) = t0_melt + (zdelwice - zdwi_max)*(lh_f*rho_w)/     &
                                              (zroc(i,kso)*zdzhs(kso))
              END IF  ! m_styp > 2
!            END IF  ! land-points only
          END DO
      END DO        ! soil layers
    ELSE IF(lmelt_var) THEN
      DO kso = 1,ke_soil
          DO i = istarts, iends
!            IF (llandmask(i)) THEN ! land points only,
              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                ztx      = t0_melt
                zw_m(i)     = zporv(i)*zdzhs(kso)
                IF(t_so_new(i,kso).LT.(t0_melt-zepsi)) THEN
                  zaa    = g*zpsis(i)/lh_f
                  zw_m(i) = zw_m(i)*((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)&
                    &         *zaa))**(-zedb(i))
                  zliquid= MAX(zepsi,w_so_now(i,kso) -  w_so_ice_now(i,kso))
                  znen   = 1._ireals-zaa*(zporv(i)*zdzhs(kso)/zliquid)**zb_por(i)
                  ztx    = t0_melt/znen
                ENDIF
                ztx      = MIN(t0_melt,ztx)
                zenergy  = zroc(i,kso)*zdzhs(kso)*(t_so_new(i,kso)-ztx)
                zdwi_max = - zenergy/(lh_f*rho_w)
                zdelwice = zdwi_max
                zwso_new  = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
                zargu = zwso_new - zw_m(i) - w_so_ice_now(i,kso)
                IF (zdelwice.LT.0.0_ireals) zdelwice =                           &
                        - MIN( - zdelwice,MIN(-zargu,w_so_ice_now(i,kso)))
                IF (zdelwice.GT.0.0_ireals) zdelwice =                           &
                          MIN(   zdelwice,MAX( zargu,0.0_ireals))
                w_so_ice_new(i,kso) = w_so_ice_now(i,kso) + zdelwice
  !             At this point we have 0.0 LE w_so_ice(i,kso,nnew) LE zwso_new
  !             If we have 0.0 LT w_so_ice(i,kso,nnew) LT zwso_new
  !             the resulting new temperature has to be equal to ztx.
  !             If not all energy available can be used to melt/freeze soil water,
  !             the following line corrects the temperature. It also applies
  !             to cases without any freezing/melting, in these cases the
  !             original temperature is reproduced.

                t_so_new(i,kso) = ztx + (zdelwice - zdwi_max)*       &
                                     (lh_f*rho_w)/(zroc(i,kso)*zdzhs(kso))
  
             END IF                   ! m_stpy > 2
!            END IF                    ! land-points only
          END DO
      ENDDO
    ENDIF   ! lmelt_var
  END IF ! lmelt


!------------------------------------------------------------------------------
! Section II.7: Energy budget and temperature prediction at snow-surface
!------------------------------------------------------------------------------

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! next line has to be changed if a soil surface temperature is
        ! predicted by the heat conduction equation
        zdtsdt (i) = (t_so_new(i,1) - zts(i))*z1d2dt
        ztsn   (i) =  t_so_new(i,1)
        IF(.NOT. lmulti_snow)  &
          ztsnown(i) = ztsn(i)       ! default setting
        zwsnew(i)       = zwsnow(i) + zdwsndt(i)*zdtdrhw

        ! forcing contributions for snow formation of dew and rime are
        ! contained in ze_ges, heat fluxes must not be multiplied by
        ! snow covered fraction

        IF(.NOT. lmulti_snow) THEN

          zrnet_snow    = sobs(i) + zthsnw(i)
          zshfl_snow(i) = zrhoch(i)*cp_d*(zth_low(i) - ztsnow(i))
          zlhfl_snow(i) = lh_s*zversn(i)
          zfor_snow     = zrnet_snow + zshfl_snow(i) + zlhfl_snow(i)

          ! forecast of snow temperature Tsnow
          IF (ztsnow(i) < t0_melt .AND. zwsnew(i) > zepsi) THEN
            ztsnown(i) = ztsnow(i) + zdt*2._ireals*(zfor_snow - zgsb(i))  &
                           /zrocs(i) - ( ztsn(i) - zts(i) )

            ! implicit formulation
! BR        zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnew(i)))
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
            zalas  = 2.22_ireals*(rho_snow_now(i)/rho_i)**1.88_ireals

            ztsnow_im    = - zrhoch(i) * (cp_d + zdqvtsnow(i) * lh_s)       &
                                         - zalas/zdz_snow_fl(i)
            zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i))
            ztsnown(i) = ztsnow(i) + (ztsnown(i)-ztsnow(i))/zfak
          END IF

          zdtsnowdt(i) = (ztsnown(i) - ztsnow(i))*z1d2dt
        ENDIF
!      END IF          ! land-points only
    END DO

  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (zwsnew(i) > zepsi) THEN
              zdtsnowdt_mult(i,ksn) = (ztsnown_mult(i,ksn) - ztsnow_mult(i,ksn))*z1d2dt
            END IF
!          ENDIF
        END DO
    ENDDO
  ENDIF


!------------------------------------------------------------------------------
! Section II.8: Melting of snow ,infiltration and surface runoff of snow water
!------------------------------------------------------------------------------

! If the soil surface temperature predicted by the equation of heat conduction
! is used instead of using T_s = T_so(1), the following section has to be
! adjusted accordingly.
!
! Basically this snow model uses heat fluxes to either heat the uppermost soil
! layer, if the snow surface temperature exceeds t0_melt, or to heat the snow, if
! the temperature of the uppermost soil layer exceeds t0_melt. Melting is considered
! after this process, if the snow temperature equals t0_melt AND the temperature of
! the uppermost soil layer exceeds t0_melt. The excess heat (t_so(1)-t0_melt) is used for
! melting snow. In cases this melting may be postponed to the next time step.

  IF (.NOT. lmulti_snow) THEN

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
    
          zwsnn  (i)  = zwsnow(i) + zdtdrhw*zdwsndt(i)
          zwsnew (i)  = zwsnn(i)
          ztsnownew     = ztsnown(i)
          ze_avail      = 0.0_ireals
          ze_total      = 0.0_ireals
          zfr_melt      = 0.0_ireals
    
          IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
            ! first case: T_snow > t0_melt: melting from above
            ! ----------
            IF (ztsnown(i) > t0_melt .AND. t_so_new(i,1) < t0_melt ) THEN
              zdwsnm(i)    = zwsnew(i)*.5_ireals*(ztsnown(i) - (t0_melt - zepsi))/ &
                               (.5_ireals* (zts(i) - (t0_melt - zepsi)) - lh_f/chc_i)
!em              zdwsnm(i)    = zdwsnm(i)*z1d2dt*rho_w 
              zdwsnm(i)    = zdwsnm(i)*z1d2dt*rho_w*zf_snow(i) 
              zdwsndt(i)   = zdwsndt (i) + zdwsnm(i)
              ztsnownew    = t0_melt - zepsi
              zdtsnowdt(i) = zdtsnowdt(i) + (ztsnownew - ztsnown(i))*z1d2dt
              runoff_s (i) = runoff_s(i) - zdwsnm(i)*zroffdt
            ENDIF ! melting from above
    
            IF (t_so_new(i,1) .gt. t0_melt) THEN
              !second case:  temperature of uppermost soil layer > t0_melt. First a
              !-----------   heat redistribution is performed. As a second step,
              !              melting of snow is considered.
              ! a) Heat redistribution
              ztsnew = t0_melt + zepsi
              ztsnownew      = ztsn(i) + ztsnown(i) - ztsnew +  &
                   2._ireals*(ztsn(i) - ztsnew)*zroc(i,1)*zdzhs(1)/zrocs(i)
              zdtsdt(i)    = zdtsdt(i) + (ztsnew - t_so_new(i,1))*z1d2dt
              zdtsnowdt(i) = zdtsnowdt(i) + (ztsnownew - ztsnown(i))*z1d2dt
              ! b) Melting of snow (if possible)
              IF (ztsnownew > t0_melt) THEN
                ze_avail     = 0.5_ireals*(ztsnownew - t0_melt)*zrocs(i)
                ze_total     = lh_f*zwsnew(i)*rho_w
                zfr_melt     = MIN(1.0_ireals,ze_avail/ze_total)
                zdtsnowdt(i)= zdtsnowdt(i) + (t0_melt - ztsnownew)*z1d2dt
                zdelt_s      = MAX(0.0_ireals,(ze_avail - ze_total)/(zroc(i,1)* &
                                                                        zdzhs(1)))
                zdtsdt(i)  = zdtsdt(i) + zdelt_s*z1d2dt
    
                ! melted snow is allowed to penetrate the soil (up to field
                ! capacity), if the soil type is neither ice nor rock (zrock = 0);
                ! else it contributes to surface run-off;
                ! fractional water content of the first soil layer determines
                ! a reduction factor which controls additional run-off
!em                zdwsnm(i)   = zfr_melt*zwsnew(i)*z1d2dt*rho_w  ! available water
                zdwsnm(i)   = zfr_melt*zwsnew(i)*z1d2dt*rho_w*zf_snow(i)  ! available water
                zdwsndt(i)  = zdwsndt (i) - zdwsnm(i)
                zdwgme        = zdwsnm(i)*zrock(i)             ! contribution to w_so
                zro           = (1._ireals - zrock(i))*zdwsnm(i)      ! surface runoff
                zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,1) -  &
                                zfcap(i))/MAX(zporv(i)-zfcap(i), zepsi)))
                zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._ireals - zredfu)
                zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                       ! for this fraction
    
                ! zro-, zdw_so_dt-correction in case of pore volume overshooting
                zw_ovpv = MAX(0._ireals, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                           zdwgdt(i,1) - zporv(i) * zdzhs(1) * zrhwddt)
                zro = zro + zw_ovpv
                zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv
    
    
                IF (zfr_melt > 0.9999_ireals) zdwsndt(i)= -zwsnow(i)*zrhwddt
                runoff_s (i)= runoff_s(i) + zro*zroffdt
    
              END IF   ! snow melting
            END IF     ! snow and/or soil temperatures
          END IF       ! points with snow cover only
!        END IF         ! land-points only
      END DO

  ELSE       ! new snow scheme

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          
          zwsnew(i) = zwsnow(i) + zdtdrhw*zdwsndt(i)
          zdwsnm(i) = 0.0_ireals

          ze_out  (i) = 0.0_ireals
          zqbase  (i) = 0.0_ireals
          zcounter(i) = 0.0_ireals

          ze_rad(i) = 0.0_ireals
          IF(zextinct(i,1).gt.0.0_ireals) ze_rad(i) = zf_snow(i) * sobs(i)

          ztsnownew_mult(i,0) = ztsnown_mult(i,0)
!        END IF         ! land-points only
      END DO

    DO ksn = 1,ke_snow

        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
    
              zrefr(i) = 0.0_ireals
              zmelt(i) = 0.0_ireals
              ztsnownew_mult(i,ksn) = ztsnown_mult(i,ksn)
    
              IF(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn).GT.zepsi .OR. &
                wtot_snow_now(i,ksn) - wliq_snow_now(i,ksn).GT.zepsi) THEN
                zrho_dry_old(i) = MAX(wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn), &
                  &                     zepsi)                                               &
                  &                 *rho_w/(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn))
              ELSE
                zrho_dry_old(i) = rho_w
              END IF
    
              ztsnownew_mult(i,ksn) = (ztsnown_mult(i,ksn)*wtot_snow_now(i,ksn) &
                &                       + t0_melt*zqbase(i)*zdt)/(zqbase(i)*zdt     &
                &                       + wtot_snow_now(i,ksn))
    
              IF(zextinct(i,ksn).eq.0.0_ireals) THEN
                ze_in = ze_out(i)
              ELSE
                IF(ksn.eq.ke_snow) THEN
                  ze_in = ze_out(i) + (ze_rad(i) - zcounter(i))
                ELSEIF(ksn.eq.1) then
                  ze_in = ze_rad(i) * (EXP (-zextinct(i,1)*zdzm_snow(i,1)) &
                    &     - EXP (-zextinct(i,1)*zdzh_snow(i,1)))
                  zcounter(i) = ze_rad(i) * (1._ireals-EXP(-zextinct(i,1)*zdzh_snow(i,1)))
                ELSE
                  ze_in = ze_out(i) + (ze_rad(i) - zcounter(i))  &
                    &     -(ze_rad(i) - zcounter(i))*EXP(-zextinct(i,ksn)*zdzh_snow(i,ksn))
                  zcounter(i) = ze_rad(i)-(ze_rad(i)-zcounter(i)) &
                    &             * EXP(-zextinct(i,ksn)*zdzh_snow(i,ksn))
                END IF
              END IF
    
              ztsnownew_mult(i,ksn) = ztsnownew_mult(i,ksn) &
                &                       + ze_in*zdt/(chc_i*wtot_snow_now(i,ksn))/rho_w
              wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) + zqbase(i)*zdt
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zqbase(i)*zdt
    
              zdzh_old = zdzh_snow(i,ksn)
              zdzh_snow(i,ksn) = zdzh_snow(i,ksn) + zqbase(i)*zdt
    
              rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*&
                &                                rho_w/zdzh_snow(i,ksn), &
                &                                0.0_ireals)
    
              IF(ztsnownew_mult(i,ksn) .GT. t0_melt) THEN
    
                IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
                  ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                    &      *z1d2dt*rho_w
                  zmelt(i) = 0.0_ireals
                ELSEIF(chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt)/lh_f <= &
                  wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn)) THEN
                  zmelt(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                    &          *z1d2dt/lh_f
                  ze_out(i) = 0.0_ireals
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt(i)*zdt
                ELSE
                  zmelt(i) = (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))*z1d2dt
                  ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                    &      *z1d2dt*rho_w - zmelt(i)*lh_f*rho_w
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt(i)*zdt
                END IF
                ztsnownew_mult(i,ksn) = t0_melt
    
              ELSE
!          T<0
                IF(wliq_snow_now(i,ksn) .GT. -chc_i*wtot_snow_now(i,ksn) &
                  & *(ztsnownew_mult(i,ksn) - t0_melt)/lh_f) THEN
                  zrefr(i) = -chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn) &
                    &          - t0_melt)*z1d2dt/lh_f
                  ztsnownew_mult(i,ksn)   = t0_melt
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zrefr(i)*zdt
                ELSE
                  zrefr(i) = wliq_snow_now(i,ksn)*z1d2dt
                  wliq_snow_now(i,ksn) = 0.0_ireals
                  ztsnownew_mult(i,ksn)   = ztsnownew_mult(i,ksn) + zrefr(i)*zdt*lh_f &
                    &                         /(chc_i*wtot_snow_now(i,ksn))
                END IF
                ze_out(i) = 0.0_ireals
    
              END IF
    
              zdtsnowdt_mult(i,ksn) = zdtsnowdt_mult(i,ksn) + &
                                        (ztsnownew_mult(i,ksn) - ztsnown_mult(i,ksn))*z1d2dt
              IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
                zqbase(i)           = wliq_snow_now(i,ksn)*z1d2dt
                wliq_snow_now(i,ksn) = 0.0_ireals
                wtot_snow_now(i,ksn) = 0.0_ireals
                zdzh_snow(i,ksn)    = 0.0_ireals
                rho_snow_mult_now(i,ksn)  = 0.0_ireals
              ELSE
                IF(zrefr(i) .GT. 0.0_ireals .OR. zmelt(i) .GT. 0.0_ireals) THEN
                  zadd_dz = 0.0_ireals
                  zadd_dz = MAX(zrefr(i),0._ireals)*(-1.0_ireals + 1.0_ireals/rho_i*rho_w)*zdt
                  zadd_dz = MAX(zmelt(i),0._ireals)*(-1.0_ireals/zrho_dry_old(i)*rho_w &
                    &       + 1.0_ireals)*zdt
                  zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) + zadd_dz
                  rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                    &                            /zdzh_snow(i,ksn),0.0_ireals)
                  IF(wtot_snow_now(i,ksn) .LE. 0.0_ireals) zdzh_snow(i,ksn) = 0.0_ireals
                  IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                    zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                    rho_snow_mult_now(i,ksn) = rho_w
                  END IF
                END IF
    
                zsn_porosity = 1._ireals - (rho_snow_mult_now(i,ksn)/rho_w -  &
                               wliq_snow_now(i,ksn)/zdzh_snow(i,ksn))/rho_i*rho_w - &
                               wliq_snow_now(i,ksn)/zdzh_snow(i,ksn)
                zsn_porosity = MAX(zsn_porosity,cwhc + 0.1_ireals)
                zp1 = zsn_porosity - cwhc
    
                IF (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) .GT. cwhc) THEN
                  zfukt             = (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) - cwhc)/zp1
                  zq0               = chcond * zfukt**3.0_ireals
                  zqbase(i)       = MIN(zq0*zdt,wliq_snow_now(i,ksn))
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zqbase(i)
                  wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) - zqbase(i)
    
                  zdzh_old = zdzh_snow(i,ksn)
                  zdzh_snow(i,ksn) = zdzh_snow(i,ksn) - zqbase(i)
                  zqbase(i)        = zqbase(i)*z1d2dt
    
                  IF(zdzh_snow(i,ksn) .LT. zepsi*0.01_ireals) THEN
                    wliq_snow_now(i,ksn) = 0.0_ireals
                    wtot_snow_now(i,ksn) = 0.0_ireals
                    zdzh_snow(i,ksn)    = 0.0_ireals
                    rho_snow_mult_now(i,ksn)  = 0.0_ireals
                  ELSE
                    rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                      &                            /zdzh_snow(i,ksn),0.0_ireals)
                    IF(wtot_snow_now(i,ksn) .LE. 0.0_ireals) zdzh_snow(i,ksn) = 0.0_ireals
                    IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                      zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                      rho_snow_mult_now(i,ksn) = rho_w
                    END IF
                  END IF
                ELSE
                  zqbase(i) = 0.0_ireals
                END IF
              END IF
    
            END IF       ! points with snow cover only
!          END IF         ! land-points only
        END DO
    END DO        ! snow layers

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
            zdwsnm(i) = zqbase(i)*rho_w       ! ksn == ke_snow
          END IF       ! points with snow cover only
!        END IF         ! land-points only
      END DO

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
            zdwsndt(i)  = zdwsndt(i) - zdwsnm(i)
    
              ! melted snow is allowed to penetrate the soil (up to field
              ! capacity), if the soil type is neither ice nor rock (zrock = 0);
              ! else it contributes to surface run-off;
              ! fractional water content of the first soil layer determines
              ! a reduction factor which controls additional run-off
    
            zdwgme        = zdwsnm(i)*zrock(i)             ! contribution to w_so
            zro           = (1._ireals - zrock(i))*zdwsnm(i)      ! surface runoff
            zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,1) -  &
                            zfcap(i))/MAX(zporv(i)-zfcap(i), zepsi)))
            zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._ireals - zredfu)
            zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                       ! for this fraction
    
            ! zro-, zdw_so_dt-correction in case of pore volume overshooting
            zw_ovpv = MAX(0._ireals, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                      zdwgdt(i,1) - zporv(i) * zdzhs(1) * zrhwddt)
            zro = zro + zw_ovpv
            zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv
    
            runoff_s(i) = runoff_s(i) + zro*zroffdt
          END IF       ! points with snow cover only
!        END IF         ! land-points only
      END DO
  
! snow densification due to gravity and metamorphism

    DO ksn = 2, ke_snow
      zp(:,ksn) = 0.0_ireals                         ! gravity, Pa
      DO k = ksn,1,-1
          DO i = istarts, iends
!            IF (llandmask(i)) THEN          ! land-points only
              zp(i,ksn) = zp(i,ksn) + rho_snow_mult_now(i,k)*g*zdzh_snow(i,ksn)
!            END IF         ! land-points only
          END DO
      END DO
    END DO
  
    DO ksn = 2, ke_snow 
        DO i = istarts, iends 
!          IF (llandmask(i)) THEN          ! land-points only 
            IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
              IF(rho_snow_mult_now(i,ksn) .LT. 600._ireals .AND. &
                rho_snow_mult_now(i,ksn) .NE. 0.0_ireals) THEN
                zdens_old = rho_snow_mult_now(i,ksn)
                zeta =         &! compactive viscosity of snow
                  ca2*EXP(19.3_ireals*rho_snow_mult_now(i,ksn)/rho_i)* &
                  EXP(67300._ireals/8.31_ireals/ztsnownew_mult(i,ksn))
                rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,ksn) + &
                  zdt*rho_snow_mult_now(i,ksn)*(csigma+zp(i,ksn))/zeta
                rho_snow_mult_now(i,ksn) = MIN(rho_snow_mult_now(i,ksn),rho_i)
                zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) * zdens_old/rho_snow_mult_now(i,ksn)
              END IF     
            END IF       ! points with snow cover only
!          END IF         ! land-points only
        END DO
    END DO

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only

            IF(ztsnownew_mult(i,0) .GT. t0_melt) THEN
              ztsnownew_mult(i,0) = t0_melt
              zdtsnowdt_mult(i,0) = zdtsnowdt_mult(i,0) +     &
                                      (ztsnownew_mult(i,0) - ztsnown_mult(i,0))*z1d2dt
            END IF
          END IF       ! points with snow cover only
!        END IF         ! land-points only
      END DO

  END IF

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN
    ! First for ksn == 0
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          t_snow_mult_new  (i,0) = t_snow_mult_now(i,0) + zdt*zdtsnowdt_mult(i,0)
          t_snow_new(i) = t_snow_mult_new (i,0)
!        ENDIF
      ENDDO
  
    DO ksn = 1, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            t_snow_mult_new  (i,ksn) = t_snow_mult_now(i,ksn) + &
              &                              zdt*zdtsnowdt_mult(i,ksn)
            dzh_snow_new     (i,ksn) = zdzh_snow(i,ksn)
            wtot_snow_new    (i,ksn) = wtot_snow_now(i,ksn)
            rho_snow_mult_new(i,ksn) = rho_snow_mult_now(i,ksn)
            wliq_snow_new    (i,ksn) = wliq_snow_now(i,ksn)
!          ENDIF
        ENDDO
    ENDDO
  ELSE
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          t_snow_new(i)  = t_snow_now(i) + zdt*zdtsnowdt(i)
!        ENDIF
      ENDDO
  ENDIF
  
    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        ! t_snow is computed above
        ! t_snow(i,nnew)  = t_snow(i,nx) + zdt*zdtsnowdt(i)
        t_so_new(i,1)  = t_so_now(i,1) + zdt*zdtsdt   (i)
  
        ! Next line has to be changed, if the soil surface temperature
        ! t_so(i,0,nnew) predicted by the heat conduction equation is used
        t_s_new   (i)    = t_so_new(i,1)
        t_so_new  (i,0)  = t_so_new(i,1)
        w_snow_new(i)  = w_snow_now(i) + zdt*zdwsndt  (i)/rho_w
        w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
  
        ! melting-point adjustment of snow: if snow temp is above freezing and a non-negligible
        ! amount of snow is available, then melt as much snow as needed to get snow temp
        ! back to t0_melt while conserving energy
        IF (t_snow_new(i) > t0_melt .AND. w_snow_new(i) > zepsi) THEN
          w_snow_new(i) = MAX(w_snow_new(i)*(1._ireals-(t_snow_new(i)-t0_melt) &
            *chc_i/lh_f), 0.0_ireals)
          t_snow_new(i) = t0_melt
        ELSE IF (w_snow_new(i) <= zepsi) THEN
          ! if the amount of snow is negligible, then just remove it
          w_snow_new(i) = 0.0_ireals
          t_snow_new(i) = t_so_new(i,0)
        ENDIF
        IF (w_i_new(i) <= zepsi) w_i_new(i) = 0.0_ireals
!      END IF          ! land-points only
    END DO


  ! Eliminate snow for multi-layer snow model, if w_snow = 0
  IF (lmulti_snow) THEN
    DO ksn = 1, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF (w_snow_new(i) <= zepsi) THEN
              t_snow_mult_new(i,ksn) = t_so_new(i,0)
              wliq_snow_new(i,ksn) = 0.0_ireals
              wtot_snow_new(i,ksn) = 0.0_ireals
              rho_snow_mult_new(i,ksn) = 0.0_ireals
              dzh_snow_new(i,ksn) = 0.0_ireals
            ENDIF
!          END IF          ! land-points only
        END DO
    END DO
  ENDIF


  IF(.NOT. lmulti_snow) THEN

      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only

!BR 7/2005 Update snow density
!
!     a) aging of existing snow
!
!     temperature dependence of relaxation/ageing constant

         ztau_snow = crhosmint+(crhosmaxt-crhosmint)*(t_snow_new(i)-csnow_tmin) &
                                                    /(t0_melt         -csnow_tmin)
         ztau_snow = MAX(0.005_ireals,MIN(crhosmaxt,ztau_snow)) ! JH change time constant to 200d!
         zrho_snowe= crhosmax_ml+(rho_snow_now(i)-crhosmax_ml)* &
                               EXP(-ztau_snow*zdt/86400._ireals)
!
!     b) density of fresh snow
!
         zrho_snowf= crhosminf+(crhosmaxf-crhosminf)* (zth_low(i)-csnow_tmin) &
                                                    /(t0_melt      -csnow_tmin)
         zrho_snowf= MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))
!
!     c) new snow density is weighted average of existing and new snow
!
         IF ( itype_gscp == 4 ) THEN
           znorm=MAX(w_snow_now(i)+(prs_gsp(i)+prs_con(i)+prg_gsp(i))      &
                     *zdtdrhw,zepsi)
           rho_snow_new(i)  = ( zrho_snowe*w_snow_now(i) + &
                              zrho_snowf*(prs_gsp(i)+prs_con(i)+prg_gsp(i)) &
                                 *zdtdrhw )    /znorm
         ELSE
           znorm=MAX(w_snow_now(i)+(prs_gsp(i)+prs_con(i) )      &
                     *zdtdrhw,zepsi)
           rho_snow_new(i)  = ( zrho_snowe*w_snow_now(i) + &
                              zrho_snowf*(prs_gsp(i)+prs_con(i) ) &
                                 *zdtdrhw )    /znorm
         ENDIF
         rho_snow_new(i) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow_new(i)))
!        END IF          ! land-points only
      END DO

  ELSE   ! new snow scheme

      DO i = istarts, iends 
!        IF (llandmask(i)) THEN  ! for landpoints only
          h_snow_new(i) = 0.0_ireals
          sum_weight(i) = 0.0_ireals
!        END IF          ! land-points only
      END DO  
    DO ksn = 1,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) .GT. zepsi) THEN
              h_snow_new(i) = h_snow_new(i) + zdzh_snow(i,ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO    
              
    DO ksn = ke_snow,1,-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) .GT. zepsi) THEN
              dz_old(i,ksn) = dzh_snow_new(i,ksn)
              z_old(i,ksn) = -sum_weight(i) - dzh_snow_new(i,ksn)/2._ireals
              sum_weight(i) = sum_weight(i) + dzh_snow_new(i,ksn)
              zhh_snow(i,ksn) = -h_snow_new(i)/ke_snow*(ke_snow-ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO   
    
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(w_snow_new(i) .GT. zepsi) THEN
            zhm_snow(i,1) = (-h_snow_new(i) + zhh_snow(i,1))/2._ireals
            dzh_snow_new(i,1) = zhh_snow(i,1) + h_snow_new(i)            !layer thickness betw. half levels of uppermost snow layer
            zdzm_snow(i,1     ) = zhm_snow(i,1) + h_snow_new(i)            !layer thickness between snow surface and main level of uppermost layer
            IF(dz_old(i,1).ne.0..and.rho_snow_mult_new(i,1).ne.0.) THEN
              wliq_snow_new(i,1) = wliq_snow_new(i,1)/dz_old(i,1)
            END IF
          END IF
!        END IF          ! land-points only
      END DO
    DO ksn = 2,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) .GT. zepsi) THEN
              zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
              dzh_snow_new(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,ksn     ) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
              IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_new(i,ksn).ne.0.) THEN
                wliq_snow_new(i,ksn) = wliq_snow_new(i,ksn)/dz_old(i,ksn)
              END IF
            END IF
!          END IF          ! land-points only
        END DO
    END DO

    DO ksn = ke_snow,1,-1
        DO i = istarts, iends
          t_new  (i,ksn) = 0.0_ireals
          rho_new(i,ksn) = 0.0_ireals
          wl_new (i,ksn) = 0.0_ireals
        END DO
    
      DO k = ke_snow,1,-1
          DO i = istarts, iends
!            IF (llandmask(i)) THEN  ! for landpoints only
              IF(w_snow_new(i) .GT. zepsi) THEN
    
                weight = MAX(MIN(z_old(i,k)+dz_old(i,k)/2._ireals,zhm_snow(i,ksn) &
                &+dzh_snow_new(i,ksn)/2._ireals)-   &
                         MAX(z_old(i,k)-dz_old(i,k)/2._ireals, &
                         zhm_snow(i,ksn)-dzh_snow_new(i,ksn)/2._ireals),0._ireals) &
                         &/dzh_snow_new(i,ksn)
    
                t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_new  (i,k)*weight
                rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_new(i,k)*weight
                wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_new(i,k)*weight
              END IF
!            END IF          ! land-points only
          END DO
      END DO
    END DO

    DO ksn = ke_snow,1,-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) > zepsi) THEN
              t_snow_mult_new  (i,ksn) = t_new  (i,ksn)
              rho_snow_mult_new(i,ksn) = rho_new(i,ksn)
              wtot_snow_new    (i,ksn) = rho_new(i,ksn)*dzh_snow_new(i,ksn)/rho_w
              wliq_snow_new    (i,ksn) = wl_new (i,ksn)*dzh_snow_new(i,ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO
    
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(w_snow_new(i) > zepsi) THEN
            rho_snow_new(i) = w_snow_new(i)/h_snow_new(i)*rho_w
          ELSE !JH
            rho_snow_new(i) = 250._ireals ! workaround need to be inspected!!
          END IF
!        END IF          ! land-points only
      END DO

  ENDIF ! lmulti_snow


  DO kso = 1,ke_soil
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          w_so_new(i,kso) = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
!        END IF  ! land-points only
      END DO
  END DO        ! soil layers

  ! Update of two-time level interface variables 
    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        h_snow(i) = h_snow_new(i)
!      END IF
    END DO
  
!---loop over tiles---
!END DO
!---------------------


!!$  IF(itype_subs .EQ. 2) THEN           ! tiles
!!$    DO ns = nsubs0+1, nsubs1, 2        ! 1 - mean, 2 - ocean, 3 - lake, 4 - no snow first
!!$        DO i = istarts, iends
!!$          IF (llandmask(i)) THEN  ! for landpoints only  !check is not necessary
!!$    
!!$            fact1 = subsfrac(i-1)+subsfrac(i)
!!$            w_snow_new(i  ) = (w_snow_new(i)  *subsfrac(i) + &
!!$                                     w_snow_new(i-1)*subsfrac(i-1))/MAX(fact1,zepsi)
!!$            w_snow_new(i-1) = 0._ireals
!!$
!!$            zf_snow_old(i) = subsfrac(i) / (subsfrac(i-1) + subsfrac(i))
!!$            zf_snow    (i) = MAX( 0.01_ireals, MIN(1.0_ireals,w_snow_new(i)/cf_snow) )* &
!!$                               zsf_heav(w_snow_new(i) - zepsi)
!!$
!!$            w_snow_new(i  ) = w_snow_new(i)/MAX(zf_snow(i),zepsi)
!!$
!!$            subsfrac(i-1) = (1._ireals - zf_snow(i))*fact1
!!$            subsfrac(i  ) = zf_snow(i)              *fact1
!!$          END IF  ! land-points only
!!$        END DO
!!$      DO kso = 1,ke_soil
!!$          DO i = istarts, iends
!!$            IF (llandmask(i)) THEN  ! for landpoints only  !check is not necessary
!!$
!!$              fact1 = (MIN(1._ireals - zf_snow_old(i), 1._ireals - zf_snow(i))) &
!!$                &     /MAX((1._ireals - zf_snow(i)),zepsi) 
!!$              fact2 = (MAX(zf_snow(i) - zf_snow_old(i), 0._ireals))/MAX(zf_snow(i), zepsi) 
!!$    
!!$              tmp1 = t_so_new    (i,kso-1)
!!$              tmp2 = w_so_new    (i,kso-1)
!!$              tmp3 = w_so_ice_new(i,kso-1)
!!$  
!!$              t_so_new    (i,kso-1) = t_so_new(i,kso-1)*fact1 &
!!$                &                         + t_so_new(i,kso)*(1._ireals - fact1)
!!$              w_so_new    (i,kso-1) = w_so_new(i,kso-1)*fact1 &
!!$                &                         + w_so_new(i,kso)*(1._ireals - fact1)
!!$              w_so_ice_new(i,kso-1) = w_so_ice_new(i,kso-1)*fact1 &
!!$                &                         + w_so_ice_new(i,kso)*(1._ireals - fact1)
!!$  
!!$              t_so_new    (i,kso) = tmp1*fact2 + t_so_new(i,kso)*(1._ireals-fact2)
!!$              w_so_new    (i,kso) = tmp2*fact2 + w_so_new(i,kso)*(1._ireals-fact2)
!!$              w_so_ice_new(i,kso) = tmp3*fact2 + w_so_ice_new(i,kso)*(1._ireals-fact2)
!!$
!!$            END IF  ! land-points only
!!$          END DO
!!$      END DO        ! soil layers
!!$    END DO
!!$ END IF


!  DO ns = nsubs0, nsubs1
  !em>
  
  ! computation of the temperature at the boundary soil/snow-atmosphere
    IF (ldiag_tg) THEN
      IF(lmulti_snow) THEN
        CALL tgcom ( t_g(:), t_snow_mult_new(:,1), t_s_new(:), &
                     w_snow_new(:), zf_snow(:), ie, cf_snow,  &
                     istarts, iends )
      ELSE
        CALL tgcom ( t_g(:), t_snow_new(:), t_s_new(:),       &
                     w_snow_new(:), zf_snow(:), ie, cf_snow, &
                     istarts, iends )
      ENDIF
    ENDIF
!  END DO


! Debug messages


#ifdef __ICON__
   IF (msg_level >= 14) THEN
!    DO ns = nsubs0, nsubs1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (w_snow_new(i) > zepsi .AND. (t_snow_new(i)<180. &
                & .OR. t_snow_new(i)>280.)) THEN 
!                & .OR. w_i_new(i)*1000. > 0.1_ireals ) THEN

              write(0,*) "SFC-DIAGNOSIS TERRA ",i,dt,nsubs1!,ntstep
              write(0,*)" nztlev ",               nztlev   
        !!$   write(0,*)" lmelt  ",               lmelt    
        !!$   write(0,*)" lmelt_var ",            lmelt_var
        !!$   write(0,*)" lmulti_snow ",          lmulti_snow 
        !!$   write(0,*)" itype_gscp ",           itype_gscp
        !!$   write(0,*)" itype_trvg ",           itype_trvg
        !!$   write(0,*)" itype_evsl ",           itype_evsl
        !!$   write(0,*)" itype_tran ",           itype_tran
        !!$   write(0,*)" itype_root ",           itype_root
        !!$   write(0,*)" itype_heatcond ",       itype_heatcond
        !!$   write(0,*)" itype_hydbound ",       itype_hydbound
        !!$   write(0,*)" lstomata  ",            lstomata
        !!$   write(0,*)" l2tls ",                l2tls  
        !!$   write(0,*)" lana_rho_snow ",        lana_rho_snow
        !!$   write(0,*)" itype_subs ",           itype_subs
        !!$   write(0,*) "zml_soil: ",czmls
              write(0,*) "t", t(i)
              write(0,*) "p0",p0(i)
              write(0,*) "qv",qv(i)
              write(0,*) "ps",ps(i)
              write(0,*) "t_g",t_g(i)
              write(0,*) "t_s",t_s_now(i),t_s_new(i)
              write(0,*) "t_snow",t_snow_now(i),t_snow_new(i)
              write(0,*) "w_snow",w_snow_now(i),w_snow_now(i)
              write(0,*) "h_snow",h_snow_now(i),h_snow_new(i)
              write(0,*) "qv_s",qv_s(i)
              write(0,*) "t_so",t_so_now(i,:),t_so_new(i,:)
              write(0,*) "w_so",w_so_now(i,:),w_so_new(i,:)
              write(0,*) "tch_t",tch(i)
              write(0,*) "tcm_t",tcm(i)
              write(0,*) " tfv_t",tfv(i)
              write(0,*) "zshfl,zlhfl,zradfl,zg1",zshfl(i),zlhfl(i),zradfl(i),zg1(i)
              write(0,*) "soiltyp_t",soiltyp_subs(i)
              write(0,*) "plcov_t",  plcov(i)
              write(0,*) "rootdp_t", rootdp(i)
              write(0,*) "sai_t",   sai(i) 
              write(0,*) "tai_t",   tai(i) 
              write(0,*) "eai_t",   eai(i) 
              write(0,*) "t_2m_t",  t_2m(i) 
              write(0,*) "u_10m_t", u_10m(i)   
              write(0,*) "v_10m_t", v_10m(i)    
              write(0,*) "sobs_t",  sobs(i)    
              write(0,*) "thbs_t",  thbs(i)     
              write(0,*) "pabs_t",  pabs(i)     
!              write(0,*) "llandmask_t",llandmask(i) 
            END IF
!          END IF
         END DO
!    END DO
  ENDIF
#endif



!------------------------------------------------------------------------------
! End of module procedure terra_multlay
!------------------------------------------------------------------------------


END SUBROUTINE terra_multlay


SUBROUTINE terra_multlay_init (                &   
                  ie,           & ! array dimensions
                  istartpar,    & ! start index for computations in the parallel program
                  iendpar,      & ! end index for computations in the parallel program
! 4. variables for the time discretization and related variables
!                  ke,           & ! nsubs0=1 for single tile, nsubs0=2 for multi-tile
!                  nsubs0ubs1,& ! nsubs1=1 for single tile, nsubs1=#tiles+1 for multi-tile
                  ke_soil, ke_snow      , &
                  czmls                 , & ! processing soil level structure 
                  soiltyp_subs      , & ! type of the soil (keys 0-9)                     --
                  rootdp       , & ! depth of the roots                            ( m  )
                  t_snow_now   , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_s_new          , & ! temperature of the ground surface             (  K  )
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  t_so_now         , & ! soil temperature (main level)                 (  K  )
                  t_so_new         , & ! soil temperature (main level)                 (  K  )
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_now     , & ! ice content                                   (m H20)
                  w_so_ice_new     , & ! ice content                                   (m H20)
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_now       & ! layer thickness between half levels in snow   (  m  )
                                  )

                                      

!-------------------------------------------------------------------------------
! Declarations for ICON (USE statements for COSMO)
!-------------------------------------------------------------------------------

IMPLICIT NONE


  INTEGER (KIND=iintegers), INTENT(IN)  ::  &
                  ie,           & ! array dimensions
                  istartpar,    & ! start index for computations in the parallel program
                  iendpar,      & ! end index for computations in the parallel program
!                  ke,           &
!                 nsubs0ubs1 , &
                  ke_soil, ke_snow      
  REAL    (KIND = ireals), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  czmls           ! processing soil level structure 

  INTEGER (KIND=iintegers),DIMENSION(ie), INTENT(IN) :: & 
                  soiltyp_subs      ! type of the soil (keys 0-9)                     --
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                  rootdp           ! depth of the roots                            ( m  )
!  LOGICAL                , DIMENSION(ie), INTENT(IN) :: & 
!                  llandmask        ! landpoint mask                                  --

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_snow_now              ! temperature of the snow-surface (K)
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_s_now              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  t_s_new              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_now         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_new         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )

!--------------------------------------------------------------------------------
! TERRA Declarations

! New declaration for ICON

!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------

  INTEGER (KIND=iintegers)              ::  &
    ierror                        ! error status variable

  CHARACTER (LEN=80)                    ::  &
    yerror
!
! Local parameters:
! ----------------

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals , & ! security constant
    zalfa  = 1.0_ireals    , & ! degree of impliciteness (1: full implicit,
                               !    (0.5: Cranck-Nicholson)
    rho_i  = 910._ireals       ! density of solid ice (soil model)  (kg/m**3)

! Local scalars:
! -------------

     INTEGER (KIND=iintegers) ::  &
       kso            , & ! loop index for soil moisture layers           
       ksn            , & ! loop index for snow layers
       i              , & ! loop index in x-direction              
       mstyp          , & ! soil type index
       istarts        , & ! start index for x-direction      
       iends          , & ! end   index for x-direction     
       ns

      REAL    (KIND=ireals   ) ::  &
       zdel_t_so          ! auxiliary variable

      REAL    (KIND=ireals   ) ::  &
      zpsi0=.01_ireals    ! air entry potential at water saturation (m)


      REAL    (KIND=ireals   ) ::  &
!
!   Statement functions
!
    zsf_heav       , & ! Statement function: Heaviside function
    zsf_psat_iw    , & ! Saturation water vapour pressure over ice or water
                       ! depending on temperature "zstx"
    zsf_qsat       , & ! Specific humidity at saturation pressure
                       ! (depending on the saturation water vapour pressure
                       !  "zspsatx" and the air pressure "zspx")
    zsf_dqvdt_iw   , & ! Statement function: First derivative of specific
                       ! saturation humidity
                       ! with respect to temperature (depending on temperature
                       ! "zstx" and saturation specific humidity pressure
                       ! "zsqsatx")
    zstx           , & ! dummy argument for Stmt. function
    zspx           , & ! dummy argument for Stmt. function
    zspsatx        , & ! dummy argument for Stmt. function
    zsqsatx        , & ! dummy argument for Stmt. function
    z2iw           , & ! dummy argument for Stmt. function
    z4iw           , & ! dummy argument for Stmt. function
    zroota    (ie                  ) ,& ! root density profile parameter (1/m)
    z234iw            ! dummy argument for Stmt. function


      REAL    (KIND=ireals   ) ::  &
       zaa                ! utility variable



     REAL    (KIND = ireals) :: &
!
    h_snow_now (ie)    , & ! snow height  (m)  
    h_snow_new (ie)    , & ! snow height  (m)  
    zmls     (ke_soil+1)  , & ! depth of soil main level
    zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
    zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
    zdzms    (ke_soil+1)  , & ! distance between main levels
!
! Multi-layer snow model
    zhh_snow (ie,ke_snow)    , & ! depth of the half level snow layers
    zhm_snow (ie,ke_snow)    , & ! depth of snow main levels
    zdzh_snow(ie,ke_snow)    , & ! layer thickness between half levels
!
! External parameters
!
    zbwt     (ie)      , & ! root depth (with artificial minimum value)
    zrock    (ie)      , & ! ice/rock-indicator: 0 for ice and rock
    zsandf   (ie)      , & ! mean fraction of sand (weight percent)
    zclayf   (ie)      , & ! mean fraction of clay (weight percent)
    zb_por   (ie)      , & ! pore size distribution index
    zpsis    (ie)      , & ! air entry potential (m)
    zw_m     (ie)          ! maximum of liquid water content  (m)

   INTEGER  (KIND=iintegers ) ::  &
    m_styp   (ie)      , & ! soil type
    zicount1 (ie)      , &
    zicount2 (ie)       

   REAL    (KIND=ireals   ) ::  & ! Soil and plant parameters
!
    zfcap    (ie)      , & ! field capacity of soil
    zadp     (ie)      , & ! air dryness point
    zporv    (ie)      , & ! pore volume (fraction of volume)
    zdlam    (ie)      , & ! heat conductivity parameter
    zdw      (ie)      , & ! hydrological diff.parameter
    zdw1     (ie)      , & ! hydrological diff.parameter
    zkw      (ie)      , & ! hydrological cond.parameter
    zkw1     (ie)      , & ! hydrological cond.parameter
    zik2     (ie)      , & ! minimum infiltration rate
    zpwp     (ie)      , & ! plant wilting point  (fraction of volume)
    zedb     (ie)          ! utility variable

     REAL    (KIND=ireals   ) ::  &
      zk0di    (ie)      , & ! surface type dependent parameter
      zbedi    (ie)          ! surface type dependent parameter

     REAL    (KIND=ireals   ) ::  &
       zrocg    (ie)        ,& ! volumetric heat capacity of bare soil
       zalam    (ie,ke_soil),& ! heat conductivity
       zw_snow_old(ie)      ,& !
       zrho_snow_old(ie)

!- End of header
!==============================================================================


!------------------------------------------------------------------------------
! Begin Subroutine terra_multlay              
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------


! Declaration of STATEMENT-FUNCTIONS

  zsf_heav     (zstx                    ) = 0.5_ireals+SIGN( 0.5_ireals, zstx )
  zsf_psat_iw  (zstx,z2iw   ,z4iw       )                                     &
                   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))
  zsf_qsat     (zspsatx, zspx           )                                     &
                   = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)
  zsf_dqvdt_iw (zstx,zsqsatx,z4iw,z234iw)                                     &
                   = z234iw*(1._ireals+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2


!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar

!>JH
!  prg_gsp=0._ireals ! graupel not implemented yet 
!<JH

  ierror = 0
  yerror = '        '


! grids for temperature and water content

  zzhls(1) = 2._ireals*czmls(1)   !depth of first half level
  zdzhs(1) = zzhls(1)      !layer thickness betw. half levels of uppermost layer
  zmls(1)  = czmls(1)      !depth of 1st main level
  zdzms(1) = czmls(1)      !layer thickness between soil surface and main level
                           ! of uppermost layer
 
  DO kso = 2,ke_soil+1
    zzhls(kso)  = zzhls(kso-1) + 2._ireals*(czmls(kso) -zzhls(kso-1))
    zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
    zmls(kso)  = czmls(kso)                ! depth of main levels
    zdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
  ENDDO



!---loop over tiles---
!DO ns=nsubs0ubs1
!---------------------


! Prepare basic surface properties (for land-points only)
 
    DO i = istarts, iends

!>JH
!      IF(llandmask(i)) THEN        ! for land-points only

        mstyp       = soiltyp_subs(i)        ! soil type
        m_styp(i) = mstyp                     ! array for soil type
        zdw   (i)  = cdw0  (mstyp)
        zdw1  (i)  = cdw1  (mstyp)
        zkw   (i)  = ckw0  (mstyp)
        zkw1  (i)  = ckw1  (mstyp)
        zik2  (i)  = cik2  (mstyp)
        zporv(i)  = cporv(mstyp)              ! pore volume
        zpwp (i)  = cpwp (mstyp)              ! plant wilting point
        zadp (i)  = cadp (mstyp)              ! air dryness point
        zfcap(i)  = cfcap(mstyp)              ! field capacity
        zrock(i)  = crock(mstyp)              ! rock or ice indicator
        zrocg(i)  = crhoc(mstyp)              ! heat capacity
        zdlam(i)  = cala1(mstyp)-cala0(mstyp) ! heat conductivity parameter
        zbwt(i)   = MAX(0.001_ireals,rootdp(i))! Artificial minimum value
                                                ! for root depth
        zroota(i) = 3._ireals/zbwt(i)       ! root density profile parameter (1/m)
                                                ! zroota=0. creates the original TERRA_LM
                                                ! version with constant root density
                                                ! for root depth
        ! New arrays for BATS-scheme
        zk0di(i)  = ck0di(mstyp)              !
        zbedi(i)  = cbedi(mstyp)              !
        ! Arrays for soil water freezing/melting
        zsandf(i)   = csandf(mstyp)
        zclayf(i)   = cclayf(mstyp)
        zpsis(i)    = -zpsi0 * 10._ireals**(1.88_ireals-0.013_ireals*zsandf(i))
        zb_por(i)   = 2.91_ireals + .159_ireals*zclayf(i)
        zedb(i)     = 1._ireals/zb_por(i)

! print*,i,zporv(i),zclayf(i),zpsis(i),zb_por(i)

!<JH
!DR        IF(m_styp(i) < 9 ) llandmask(i) = .true.
!>JH      ENDIF

    ENDDO
 

  ! Set three-dimensional variables
  DO kso = 1, ke_soil
      DO i = istarts, iends
!        IF(llandmask(i)) THEN        ! for land-points only
          mstyp           = soiltyp_subs(i)        ! soil type
          zalam(i,kso)  = cala0(mstyp)              ! heat conductivity parameter
!        ENDIF
      ENDDO
  ENDDO




! For ntstep=0 : Some preparations
! ================================

!!$  IF (ntstep == 0) THEN
    w_so_new(:,:) = w_so_now(:,:)

!   Provide for a soil moisture 1 % above air dryness point, reset soil
!   moisture to zero in case of ice and rock
    DO kso   = 1,ke_soil+1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN             ! for land-points only
            IF (m_styp(i).ge.3) THEN
              w_so_now (i,kso) = MAX(w_so_now(i,kso),                     &
                                           1.01_ireals*zadp(i)*zdzhs(kso) )
              w_so_new (i,kso) = MAX(w_so_new(i,kso),                     &
                                           1.01_ireals*zadp(i)*zdzhs(kso) )
            ELSE
              w_so_now(i,kso) = 0.0_ireals
              w_so_new(i,kso) = 0.0_ireals
            ENDIF
!          END IF   ! land-points
        END DO
    END DO


!   adjust temperature profile in lower soil layers, if temperature of first soil
!   layer was reduced due to the detection of snow (e.g. in the analysis)
!   loop over grid points
      DO i = istarts, iends
!        IF(llandmask(i)) THEN   ! for land-points only
          IF (w_snow_now(i) <=1.0E-6_ireals) THEN
            ! spurious snow is removed
            t_snow_now(i)=t_so_now(i,0)
!jh         t_snow_new(i)=t_so_now(i,0)
            w_snow_now(i)= 0.0_ireals
            IF(lmulti_snow) THEN
              t_snow_mult_now(i,0) = t_so_now(i,0)
              DO ksn = 1, ke_snow
                t_snow_mult_now(i,ksn) = t_so_now(i,0)
                wliq_snow_now(i,ksn) = 0.0_ireals
                wtot_snow_now(i,ksn) = 0.0_ireals
                rho_snow_mult_now(i,ksn) = 0.0_ireals
                dzh_snow_now(i,ksn) = 0.0_ireals
              END DO
            END IF
          ELSE
!           adjust soil temperatures in presence of snow
            zdel_t_so = MAX (0._ireals,t_so_now(i,0)-t0_melt)
            t_so_now(i,0)=MIN( t0_melt,t_so_now(i,0) )
            t_so_new(i,0)=MIN( t0_melt,t_so_now(i,0) )
            DO kso=1,ke_soil
              IF ( t_so_now(i,kso) > t0_melt) THEN
                 t_so_now(i,kso) = MAX( t0_melt,                  &
                                   t_so_now(i,kso) -zdel_t_so    &
                                   *(zmls(ke_soil+1)-zmls(kso))   &
                                   /(zmls(ke_soil+1)-zmls( 1 )) )
                 t_so_new(i,kso) = MAX( t0_melt,                  &
                                   t_so_now(i,kso) -zdel_t_so    &
                                   *(zmls(ke_soil+1)-zmls(kso))   &
                                   /(zmls(ke_soil+1)-zmls( 1 )) )
              ENDIF
              IF ( t_so_now(i,kso) <= t0_melt) EXIT
            ENDDO
          ENDIF

          t_s_now(i) = t_so_now(i,0)
          t_s_new(i) = t_so_new(i,0)

!         Set level 1 to level 0 for t_so for every landpoint
          t_so_now(i,1) = t_so_now(i,0)
          t_so_new(i,1) = t_so_now(i,0)

!        ENDIF    ! llandmask
      END DO


!   Initialization of soil ice content
!   ----------------------------------
!   Generally the combination lmelt=.true., lemlt_var=.true. has to be used.
!   lmelt_var=.false. (soil water freezing at 273.15K) should be used only
!   for special investigations because at model start a partial frozen layer
!   is only considered if the soil ice water content is read in from
!   a pre run.
!   ----------------------------------
    IF (lmelt .AND. .NOT. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
          DO i = istarts, iends
!            IF (llandmask(i)) THEN             ! for land-points only
              w_so_ice_now(i,kso) = 0.0_ireals
              w_so_ice_new(i,kso) = 0.0_ireals

              IF (t_so_now(i,kso) < t0_melt) THEN
                w_so_ice_now(i,kso) = w_so_now(i,kso)
                w_so_ice_new(i,kso) = w_so_now(i,kso)
              END IF ! t_so(kso) < t0_melt
!            END IF   ! land-points
          END DO
      END DO
    END IF           ! lmelt .AND. .NOT. lmelt_var
    IF(lmelt .AND. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
          DO i = istarts, iends
!            IF (llandmask(i)) THEN             ! for land-points only
              IF (t_so_now(i,kso) < (t0_melt-zepsi)) THEN
                zaa    = g*zpsis(i)/lh_f
                zw_m(i)     = zporv(i)*zdzhs(kso)
                zw_m(i)   = zw_m(i)*                                          &
                  ((t_so_now(i,kso) - t0_melt)/(t_so_now(i,kso)*zaa))**(-zedb(i))
                w_so_ice_now(i,kso) = MAX (0.0_ireals,w_so_now(i,kso) - zw_m(i))
                w_so_ice_new(i,kso) = MAX (0.0_ireals,w_so_now(i,kso) - zw_m(i))
              END IF
!            END IF   ! land-points
          END DO
      END DO
    END IF           ! lmelt .AND. lmelt_var


!   Initialization of snow density, if necessary
!   --------------------------------------------
    IF (.NOT. lana_rho_snow) THEN
        DO i = istarts, iends
!          IF (llandmask(i)) THEN             ! for land-points only
            IF(lmulti_snow) THEN
              DO ksn = 1, ke_snow
                rho_snow_mult_now(i,ksn) = 250.0_ireals    ! average initial density
                t_snow_mult_now  (i,ksn) = t_snow_now(i)
                wtot_snow_now    (i,ksn) = w_snow_now(i)/REAL(ke_snow,ireals)
                dzh_snow_now     (i,ksn) = w_snow_now(i)/REAL(ke_snow,ireals) &
                  &                            *rho_w/rho_snow_mult_now(i,ksn)
                wliq_snow_now    (i,ksn) = 0.0_ireals
              END DO
            ELSE
              rho_snow_now(i) = 250.0_ireals    ! average initial density
            ENDIF
!          ENDIF   ! land-points
        ENDDO
    ELSE
      IF(lmulti_snow) THEN
!The sum of wtot_snow, i.e. the "old" total water equivalent depth
        zw_snow_old  (:) = 0.0_ireals
        zrho_snow_old(:) = 0.0_ireals
        DO ksn = 1, ke_snow
          DO i = istarts, iends
            zw_snow_old(i) = zw_snow_old(i) + wtot_snow_now(i,ksn)
            zrho_snow_old(i) = zrho_snow_old(i) + dzh_snow_now(i,ksn)
          END DO
        END DO
!The "old" snow density
        DO i = istarts, iends
          zrho_snow_old(i) = zw_snow_old(i) / MAX(zrho_snow_old(i),1.E-09) *rho_w
        END DO 
        DO i = istarts, iends
          zicount1(i) = COUNT(wtot_snow_now(i,:).gt.zepsi)
          zicount2(i) = COUNT(dzh_snow_now(i,:).gt.zepsi)
        END DO
        DO ksn = 1, ke_snow
          DO i = istarts, iends
            IF(zicount1(i).EQ.ke_snow .AND. zicount2(i).EQ.ke_snow) THEN
              rho_snow_mult_now(i,ksn) = wtot_snow_now(i,ksn)/dzh_snow_now(i,ksn)*rho_w
              wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn)*w_snow_now(i)/zw_snow_old(i)
              dzh_snow_now (i,ksn) = dzh_snow_now (i,ksn)*w_snow_now(i)/zw_snow_old(i)
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)*w_snow_now(i)/zw_snow_old(i)
            ELSE
              IF(rho_snow_now(i) .EQ. 0._ireals) rho_snow_now(i) = 250._ireals
              rho_snow_mult_now(i,ksn) = rho_snow_now(i)
              wtot_snow_now    (i,ksn) = w_snow_now(i)/ke_snow
              dzh_snow_now     (i,ksn) = w_snow_now(i)*rho_w/rho_snow_mult_now(i,ksn)/ke_snow
              wliq_snow_now    (i,ksn) = 0.0_ireals
            END IF
          END DO
        END DO
      ELSE
          DO i = istarts, iends
!            IF (llandmask(i)) THEN             ! for land-points only
              IF(rho_snow_now(i) .EQ. 0._ireals) rho_snow_now(i) = 250._ireals
!            END IF
          END DO
      END IF
    ENDIF


!   Initialization of the local array of the grid in snow
!   -----------------------------------------------------
    IF(lmulti_snow) THEN
        DO i = istarts, iends
!          IF(llandmask(i)) THEN   ! for land-points only
            h_snow_new(i) = 0.0_ireals
            DO ksn = 1,ke_snow
              zdzh_snow(i,ksn) = dzh_snow_now(i,ksn)
              h_snow_new(i) = h_snow_new(i) + zdzh_snow(i,ksn)
            END DO
            h_snow_now(i) = h_snow_new(i)

            zhh_snow(i,1) = - h_snow_now(i) + dzh_snow_now(i,1)
            DO ksn = 2,ke_snow
              zhh_snow(i,ksn) = zhh_snow(i,ksn-1) + dzh_snow_now(i,ksn)
            END DO

            zhm_snow(i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._ireals
            DO ksn = 2,ke_snow
              zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
            END DO
!          ENDIF    ! llandmask
        END DO
    END IF

!!$ ENDIF             ! ntstep = 0

!---loop over tiles---
!END DO !ns=nsubs0,nsubs1
!---------------------

! End of timestep 0 preparations
! ==============================

END SUBROUTINE terra_multlay_init
 

!==============================================================================


SUBROUTINE tgcom (tg, ts, tb, ws, snowfrac, ie, cf_snow, istart, iend)

!-------------------------------------------------------------------------------
!
! Description:
!   Computation of the temperature tg at the boundary layer between the ground
!   and the atmosphere. Only 2-dimensional arrays can be passed to tgcom. It
!   must be called using the desired time level.
!
! Method:
!   For grid points above water and for grid points on land that are not
!   covered with snow:   tg = ground surface temperature tb
!   For snow covered land points, tg is a function of the temperature of the
!   the snow surface ts and the ground surface temperature tb:
!       tg = ts + exp( -rhde*ws ) * (tb-ts)
!   from Version 2.18 on replaced by
!       tg = ts + ( 1. - MIN(1.,ws/cf_snow)) * (tb -ts)
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie,              & ! dimensions of the fields
  istart, iend       ! start and end-indices of the computation

REAL (KIND=ireals), INTENT (INOUT)       ::    &
  tg (ie)    ! temperature at the boundary between ground and atmosphere

REAL (KIND=ireals), INTENT (IN)          ::    &
!DR  ts (ie), & ! temperature of the snow surface
  tb (ie), & ! temperature of the ground surface
  ws (ie)    ! water content of snow

REAL (KIND=ireals), INTENT (INOUT)          ::    &
  ts (ie),   & ! temperature of the snow surface
  snowfrac(ie) ! snow-cover fraction

!LOGICAL,  INTENT (IN)                    ::    &
!  llp (ie)   ! pattern of land- and sea-points

REAL (KIND=ireals), INTENT (IN)          ::    &
  cf_snow       ! factor for the computation

INTEGER :: i

!-------------------------------------------------------------------------------

! Begin subroutine tgcom

  DO i = istart, iend
    snowfrac(i) = MIN(1.0_ireals,ws(i)/cf_snow)
    tg(i) = ts(i) + (1.0_ireals - snowfrac(i))*(tb(i) - ts(i))
  ENDDO

END SUBROUTINE tgcom


!==============================================================================

!------------------------------------------------------------------------------
! End of module src_soil_multlay
!------------------------------------------------------------------------------




END MODULE mo_soil_ml

