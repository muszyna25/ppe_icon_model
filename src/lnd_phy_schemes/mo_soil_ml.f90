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
    je,           & ! number of grid points in meridional direction
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
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program
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
USE mo_lnd_nwp_config
!                           
USE mo_exception,       ONLY: message, finish, message_text
#endif



!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

!IMPLICIT NONE
PRIVATE:: normalize

!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: terra_multlay

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

SUBROUTINE terra_multlay (                &   
                  ie, je                , & ! array dimensions
                  istartpar,    & ! start index for computations in the parallel program
                  iendpar,      & ! end index for computations in the parallel program
                  jstartpar,    & ! start index for computations in the parallel program
                  jendpar,      & ! end index for computations in the parallel program
! 4. variables for the time discretization and related variables
                  ke,nsubs0,nsubs1      , &
                  ke_soil, ke_snow      , &
                  czmls                 , & ! processing soil level structure 
                  dt                    , &
!
                  soiltyp_subs      , & ! type of the soil (keys 0-9)                     --
                  plcov        , & ! fraction of plant cover                         --
                  rootdp       , & ! depth of the roots                            ( m  )
                  sai          , & ! surface area index                              --
                  tai          , & ! transpiration area index                        --
                  eai          , & ! earth area (evaporative surface area) index     --
                  landmask    , & ! landpoint mask                                  --
                  rsmin2d      ,  & ! minimum stomata resistance                    ( s/m )
!
                  u            , & ! zonal wind speed                              ( m/s )
                  v            , & ! meridional wind speed                         ( m/s )
                  t            , & ! temperature                                   (  k  )
                  qv           , & ! specific water vapor content                  (kg/kg)
                  p0           , & !!!! base state pressure                           (Pa) 
                  pp           , & ! deviation from the reference pressure         ( pa  )
                  ps           , & ! surface pressure                              ( pa  )
!
                  t_snow       , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult  , & ! temperature of the snow-surface               (  K  )
                  t_s          , & ! temperature of the ground surface             (  K  )
                 t_g          , & ! weighted surface temperature                  (  K  )
                 qv_s         , & ! specific humidity at the surface              (kg/kg)
                  w_snow       , & ! water content of snow                         (m H2O)
                  rho_snow     , & ! snow density                                  (kg/m**3)
                  rho_snow_mult, & ! snow density                                  (kg/m**3)
                  h_snow       , & ! snow height                                   (  m  
                  w_i          , & ! water content of interception water           (m H2O)
                  t_so         , & ! soil temperature (main level)                 (  K  )
                  w_so         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice     , & ! ice content                                   (m H20)
                  t_2m         , & ! temperature in 2m                             (  K  )
                  u_10m        , & ! zonal wind in 10m                             ( m/s )
                  v_10m        , & ! meridional wind in 10m                        ( m/s )
                  freshsnow    , & ! indicator for age of snow in top of snow layer(  -  )
                  wliq_snow    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow     , & ! layer thickness between half levels in snow   (  m  )
                  subsfrac     , &
!
                  prr_con      , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con      , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp      , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp      , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  prg_gsp      , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
!
                  tch          , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm          , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv          , & ! laminar reduction factor for evaporation      ( -- )
! 
                  sobs         , & ! solar radiation at the ground                 ( W/m2)
                  thbs         , & ! thermal radiation at the ground               ( W/m2)
                  pabs         , & !!!! photosynthetic active radiation               ( W/m2)
! 
                  runoff_s     , & ! surface water runoff; sum over forecast      (kg/m2)
                  runoff_g     , & ! soil water runoff; sum over forecast         (kg/m2)
!                 nstart       , & ! first time step of the forecast
                  ntstep       , & ! actual time step
                                   ! indices for permutation of three time levels
!                 nold         , & ! corresponds to ntstep - 1
                  nnow         , & ! corresponds to ntstep 
                  nnew           & ! corresponds to ntstep + 1
                                         )


!-------------------------------------------------------------------------------
! Declarations for ICON (USE statements for COSMO)
!-------------------------------------------------------------------------------

IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN)  ::  &
!                  nstart       , & ! first time step of the forecast
                  ntstep       , & ! actual time step
                                   ! indices for permutation of three time levels
!                 nold         , & ! corresponds to ntstep - 1
                  nnow         , & ! corresponds to ntstep 
                  nnew             ! corresponds to ntstep + 1
!                  nztlev              ! time step scheme 2,3 

  INTEGER (KIND=iintegers), INTENT(IN)  ::  &
                  ie, je           , & ! array dimensions
                  istartpar,    & ! start index for computations in the parallel program
                  iendpar,      & ! end index for computations in the parallel program
                  jstartpar,    & ! start index for computations in the parallel program
                  jendpar,      & ! end index for computations in the parallel program
                  ke,nsubs0,nsubs1 , &
                  ke_soil, ke_snow      
  REAL    (KIND = ireals), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  czmls           ! processing soil level structure 
  REAL    (KIND = ireals), INTENT(IN)  ::  &
                  dt                    

  INTEGER,                 DIMENSION(ie,je,nsubs1), INTENT(IN) :: & 
                  soiltyp_subs      ! type of the soil (keys 0-9)                     --
  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(IN) :: & 
                  plcov        , & ! fraction of plant cover                         --
                  rootdp       , & ! depth of the roots                            ( m  )
                  sai          , & ! surface area index                              --
                  tai          , & ! transpiration area index                        --
                  eai              ! earth area (evaporative surface area) index     --
  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(INOUT) :: & 
                  landmask        ! landpoint mask fractions                                 --
  REAL    (KIND = ireals), DIMENSION(ie,je), INTENT(IN) :: & 
                 rsmin2d          ! minimum stomata resistance                    ( s/m )

  REAL    (KIND = ireals), DIMENSION(ie,je,ke,nztlev), INTENT(IN) :: & 
                  u            , & ! zonal wind speed                              ( m/s )
                  v            , & ! meridional wind speed                         ( m/s )
                  t            , & ! temperature                                   (  k  )
                  qv               ! specific water vapor content                  (kg/kg)
  REAL    (KIND = ireals), DIMENSION(ie,je,ke), INTENT(IN) :: & 
                  p0               !!!! base state pressure                           (Pa) 
  REAL    (KIND = ireals), DIMENSION(ie,je,ke,nztlev), INTENT(IN) :: & 
                  pp               ! deviation from the reference pressure         ( pa  )

  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev), INTENT(IN) ::    &
                  ps               ! surface pressure                               ( pa  )

  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev,nsubs1), INTENT(INOUT) :: &
                  t_snow           ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,je,0:ke_snow,nztlev,nsubs1), INTENT(INOUT) :: &
                  t_snow_mult      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev,nsubs1), INTENT(INOUT) :: &
                  t_s              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev,nsubs1)  :: &
                  t_g          , & ! weighted surface temperature                  (  K  )
                  qv_s             ! specific humidity at the surface              (kg/kg)
  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev,nsubs1), INTENT(INOUT) :: &
                  w_snow       , & ! water content of snow                         (m H2O)
                  rho_snow         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,je,ke_snow,nztlev,nsubs1), INTENT(INOUT) :: &
                  rho_snow_mult    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev,nsubs1), INTENT(INOUT) :: &
                  h_snow           ! snow height                                   (  m   )  
  REAL    (KIND = ireals), DIMENSION(ie,je,nztlev,nsubs1), INTENT(INOUT) :: &
                  w_i              ! water content of interception water           (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,je,0:ke_soil+1,nztlev,nsubs1), INTENT(INOUT) :: &
                  t_so             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,je,ke_soil+1,nztlev,nsubs1), INTENT(INOUT) :: &
                  w_so         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(INOUT) :: &
                  t_2m         , & ! temperature in 2m                             (  K  )
                  u_10m        , & ! zonal wind in 10m                             ( m/s )
                  v_10m        , & ! meridional wind in 10m                        ( m/s )
                  freshsnow        ! indicator for age of snow in top of snow layer(  -  )
   REAL    (KIND = ireals), DIMENSION(ie,je,ke_snow,nztlev,nsubs1), INTENT(INOUT) :: &
                  wliq_snow    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,je,ke_snow,nztlev,nsubs1), INTENT(INOUT) :: &
                  dzh_snow         ! layer thickness between half levels in snow   (  m  )
  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(INOUT) :: &
                  subsfrac

  REAL    (KIND = ireals), DIMENSION(ie,je), INTENT(IN) ::    &
                  prr_con      , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con      , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp      , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          ! precipitation rate of snow, grid-scale        (kg/m2*s)
  REAL    (KIND = ireals), DIMENSION(ie,je), INTENT(INOUT) ::    &
                  prg_gsp          ! precipitation rate of graupel, grid-scale     (kg/m2*s)

  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(INOUT) :: &
                  tch          , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm          , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv              ! laminar reduction factor for evaporation      ( -- )
!
  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(IN) :: &
                  sobs         , & ! solar radiation at the ground                 ( W/m2)
                  thbs         , & ! thermal radiation at the ground               ( W/m2)
                  pabs             !!!! photosynthetic active radiation               ( W/m2)
!
  REAL    (KIND = ireals), DIMENSION(ie,je,nsubs1), INTENT(INOUT) :: &
                  runoff_s     , & ! surface water runoff; sum over forecast      (kg/m2)
                  runoff_g         ! soil water runoff; sum over forecast         (kg/m2)
!!$  REAL    (KIND = ireals), DIMENSION(ie,je), INTENT(INOUT) :: &
!!$                  rstom            ! stomata resistance                           ( s/m )
!!$  REAL    (KIND = ireals), DIMENSION(ie,je), INTENT(INOUT) :: &
!!$                  lhfl_bs          ! average latent heat flux from bare soil evap.( W/m2)
!!$  REAL    (KIND = ireals), DIMENSION(ie,je,ke_soil), INTENT(INOUT) :: &
!!$                  lhfl_pl          ! average latent heat flux from plants         ( W/m2)


!--------------------------------------------------------------------------------
! TERRA Declarations

! New declaration for ICON

  LOGICAL     llandmask(ie,je,nsubs1)        ! landpoint mask                                  --

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
    nx             , & ! time-level for integration
    kso            , & ! loop index for soil moisture layers           
    ksn            , & ! loop index for snow layers
    k              , & ! loop index for snow layers
    kso_zag        , & ! loop index for snow layers
    ke_soil_hy     , & ! number of active soil moisture layers
    i              , & ! loop index in x-direction              
    j              , & ! loop index in y-direction              
    im1, jm1       , & ! i-1, j-1
    jb             , & ! loop index for soil-type               
    mstyp          , & ! soil type index
    msr_off        , & ! number of layers contributing to surface run off
    istarts        , & ! start index for x-direction      
    iends          , & ! end   index for x-direction     
    jstarts        , & ! start index for y-direction     
    jends          , & ! end   index for y-directin      
    i250           , & ! number of layers above a soil depth of 2.50 m
    k10cm          , & ! index of half level closest to 0.1 m
    k100cm         , & ! index of half level closest to 1.0 m
!subs amf
    ns
!amf ---

  REAL    (KIND=ireals   ) ::  &
!
!   Timestep parameters
!
    zdt            , & ! integration time-step [s]
    zdelt_bound    , & ! for surface temperature problem at lateral boundaries
    zdel_t_so      , & ! auxiliary variable
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
    zdz_snow(ie,je)    ! snow depth (not for snow heat content but snow
                       ! density calculation)

  REAL    (KIND=ireals   ) ::  &
!
!   Multi-snow layer parameters
    zsn_porosity   , & ! snow porosity                    (   -   )
    zp1            , &
    zfukt          , &
    zq0            , &
    zqbase(ie,je)  , &
    zdzh_old       , &
    zrefr(ie,je)   , & ! rate of liquid water refreezing
    zmelt(ie,je)   , & ! rate of snow melting
    ze_in          , &
    ze_out         , &
    zflag          , &
    zadd_dz        , &
    zrho_dry_old(ie,je)   , &
    zeta           , &
    zdens_old      , &
    zp             , &
    zcounter       , &
    ze_rad         , &
    zswitch        , &
    fact1          , &
    fact2          , &
    tmp1           , &
    tmp2           , &
    tmp3           , &
    zf_snow_old(ie,je)    , &
!
!   Plant parameters
    zbeta          , & ! reduction factor for evaporation
    zevap              ! auxiliary variable

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
    zroot          , & ! fraction of layer filled by roots (Bucket scheme)
    zropart        , & ! fraction of layer filled by roots (Dickinson scheme)
    zrootfc        , & ! distributes transpiration to soil layers
    zr_root        , & ! zroot/zbwt (Bucket scheme)
    ze_sum             ! sum of all contributions to evapotranspiration

  REAL    (KIND=ireals   ) ::  &
!
!   Thermal parameters
!
    ztgt0          , & ! Indicator T_g > T_0
    zgstr          , & ! downward longwave radiation
    zrnet_s        , & ! net radiation
    zshfl_s        , & ! sensible heatflux at soil surface
    zlhfl_s        , & ! latent heatflux at soil surface
    zsprs          , & ! utility variable
    zalas          , & ! heat conductivity of snow
    zrnet_snow     , & ! net radiation at snow surface
    zfak           , & ! utility variable for implicit snow temperature forecast
    ztsnow_im      , & ! utility variable for implicit snow temperature forecast
    ztsnew         , & ! preliminary value of snow surface temperature
    ztsnownew      , & ! preliminary value of snow surface temperature
    zshfl_snow     , & ! sensible heatflux at snow surface
    zlhfl_snow     , & ! latent heatflux at snow surface
    zfor_snow      , & ! total forcing at snow surface
    zfr_melt       , & ! melting snow fraction
!    zdwsnm         , & ! utility variable for snow melt determination
    zdwgme         , & ! utility variable for snow melt determination
    zdelt_s        , & ! utility variable for snow melt determination
    ze_avail       , & ! utility variable for snow melt determination
    ze_total           ! utility variable for snow melt determination

  ! Some of these values have to be fields for the multi-layer snow model
  REAL    (KIND=ireals   ) ::  &
    zalas_mult    (ie,je,ke_snow),    & ! heat conductivity of snow
    ztsnownew_mult(ie,je,0:ke_snow),  & ! preliminary value of snow surface temperature
    zextinct      (ie,je,ke_snow),    & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(ie,je),            & ! total forcing at snow surface
    zfor_total    (ie,je)               ! total forcing at surface

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
    w_p            , & ! preliminary value of soil water content
!    zwsnew         , & ! preliminary value of snow water equivalent
    zw_ovpv        , & ! utility variable
!
!   Implicit solution of thermal and hydraulic equations
!
    zakb           , & ! utility variable
    zakb_m1        , & ! utility variable
    zakb_p1        , & ! utility variable
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
   zg1       (ie         ,je         ) ,& ! heat flux between two uppermost soil layers (W/m**2)
   zlhfl     (ie         ,je         ) ,& ! estimated       latent    heat flux surface (W/m**2)
   zshfl     (ie         ,je         ) ,& ! estimated       sensible  heat flux surface (W/m**2)
   zthfl     (ie         ,je         ) ,& ! estimated total turbulent heat flux surface (W/m**2)
   zradfl    (ie         ,je         ) ,& ! total radiative flux at surface             (W/m**2)
   ze_melt   (ie         ,je         ) ,& ! energy required for melting of snow         (J/m**2)
   zch_snow  (ie         ,je         ) ,& ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
   zeb1      (ie         ,je         ) ,& ! estimated energy budget of first soil layer (W/m**2)
   ztchv     (ie         ,je         ) ,& ! transfer coefficient multiplied by wind velocity
   ztchv_max (ie         ,je         ) ,& ! dito, but upper limit
   zrho_atm  (ie         ,je         ) ,& ! air density of lowest atmospheric layer     (kg/m**3)
   zdt_atm   (ie         ,je         ) ,& ! temperature difference between lowest
   zdq_atm   (ie         ,je         ) ,& ! dito for moisture (kg/kg)
   !additional variables for root distribution
   zroota    (ie         ,je         ) ,& ! root density profile parameter (1/m)
   zwrootdz  (ie         ,je, ke_soil) ,& ! mean water content over root depth weighted by root density
   zrootdz_int (ie       ,je         ) ,& ! parameter needed to initialize the root density profile integral
   zwrootdz_int(ie       ,je         )    ! parameter needed to initialize the root water content integral


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
! Model geometry
!
  zmls     (ke_soil+1)  , & ! depth of soil main level
  zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
  zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
  zdzms    (ke_soil+1)  , & ! distance between main levels
  zdz_snow_fl(ie,je)    , & ! snow depth for snow temperature gradient
!
! Multi-layer snow model
  zhh_snow (ie,je,ke_snow)    , & ! depth of the half level snow layers
  zhm_snow (ie,je,ke_snow)    , & ! depth of snow main levels
  zdzh_snow(ie,je,ke_snow)    , & ! layer thickness between half levels
  zdzm_snow(ie,je,ke_snow)    , & ! distance between main levels
!
! External parameters
!
  zbwt     (ie,je)      , & ! root depth (with artificial minimum value)
  zrock    (ie,je)      , & ! ice/rock-indicator: 0 for ice and rock
  zsandf   (ie,je)      , & ! mean fraction of sand (weight percent)
  zclayf   (ie,je)      , & ! mean fraction of clay (weight percent)
  zb_por   (ie,je)      , & ! pore size distribution index
  zpsis    (ie,je)      , & ! air entry potential (m)
  zw_m     (ie,je)          ! maximum of liquid water content  (m)

  INTEGER  (KIND=iintegers ) ::  &
  m_styp   (ie,je)          ! soil type

  REAL    (KIND=ireals   ) ::  &
!
! Connection to the atmosphere
!
  zrr      (ie,je)      , & ! total rain rate including formation of dew
  zrs      (ie,je)      , & ! total snow rate including formation of rime
  zesoil   (ie,je)      , & ! evaporation from bare soil
  zrhoch   (ie,je)      , & ! transfer coefficient*rho*g
  zf_snow  (ie,je)      , & ! surface fraction covered by snow
  zth_low  (ie,je)      , & ! potential temperature of lowest layer
  zf_wi    (ie,je)      , & ! surface fraction covered by interception water
  ztmch    (ie,je)      , & ! heat transfer coefficient*density*velocity
  zep_s    (ie,je)      , & ! potential evaporation for t_s
  zep_snow (ie,je)      , & ! potential evaporation for t_snow
  zverbo   (ie,je)      , & ! total evapotranspiration
  zversn   (ie,je)      , & ! total evaporation from snow surface
  zthsoi   (ie,je)      , & ! thermal flux at soil surface
  zthsnw   (ie,je)      , & ! thermal flux at snow surface
  zfor_s   (ie,je)      , & ! total forcing at soil surface
  zgsb     (ie,je)      , & ! heat-flux through snow
!
! Tendencies
!
  zdwidt   (ie,je)      , & ! interception store tendency
  zdwsndt  (ie,je)      , & ! snow water content tendency
  zdtsdt   (ie,je)      , & ! tendency of zts
  zdtsnowdt(ie,je)      , & ! tendency of ztsnow
  zdtsnowdt_mult(ie,je,0:ke_snow)      , & ! tendency of ztsnow
  zdwgdt   (ie,je,ke_soil)  ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=ireals   ) ::  &
!
! Soil and plant parameters
!
  zroc(ie,je,ke_soil+1) , & ! heat capacity
  zfcap    (ie,je)      , & ! field capacity of soil
  zadp     (ie,je)      , & ! air dryness point
  zporv    (ie,je)      , & ! pore volume (fraction of volume)
  zdlam    (ie,je)      , & ! heat conductivity parameter
  zdw      (ie,je)      , & ! hydrological diff.parameter
  zdw1     (ie,je)      , & ! hydrological diff.parameter
  zkw      (ie,je)      , & ! hydrological cond.parameter
  zkw1     (ie,je)      , & ! hydrological cond.parameter
  zik2     (ie,je)      , & ! minimum infiltration rate
  zpwp     (ie,je)      , & ! plant wilting point  (fraction of volume)
  ztlpmwp  (ie,je)      , & ! turgor-loss-point minus plant wilting point
  zedb     (ie,je)      , & ! utility variable
!
! Hydraulic variables
!
  ztrang(ie,je,ke_soil) , & ! transpiration contribution by the different layers
  ztrangs(ie,je)        , & ! total transpiration (transpiration from all
                            !    soil layers)
  zwin     (ie,je)      , & ! water cont. of interception store   (m H20)
  zwsnow   (ie,je)      , & ! snow water equivalent               (m H20)
  zwsnew   (ie,je)      , & ! snow water equivalent               (m H20)
  zdwsnm   (ie,je)      , & ! utility variable for snow melt determination
  zw_fr    (ie,je,ke_soil+1),&!fractional total water content of soil layers
  zinfil   (ie,je)      , & ! infiltration rate
  zlw_fr   (ie,je,ke_soil+1), & ! fractional liqu. water content of soil layer
  ziw_fr   (ie,je,ke_soil+1), & ! fractional ice content of soil layer
  zwsnn    (ie,je)          , & ! new value of zwsnow
  zflmg    (ie,je,ke_soil+1), & ! flux of water at soil layer interfaces
  zrunoff_grav(ie,je,ke_soil+1) ! main level water gravitation
                                !    flux (for runoff calc.)

  REAL    (KIND=ireals   ) ::  &
!
! Local arrays for BATS-scheme
!
  zk0di    (ie,je)      , & ! surface type dependent parameter
  zbedi    (ie,je)      , & ! surface type dependent parameter
  zsnull   (ie,je)      , & ! mean relative(zporv) water fraction of the
                            !    so-called total active soil layer (above 1m)
  zs1      (ie,je)      , & ! mean relative (zporv) water fraction  of
                            !    layers above 0.1m
  zf_rad   (ie,je)      , & ! radiation function for stomata resistance
  ztraleav (ie,je)      , & ! transpiration rate of dry leaves
!
! Root functions
!
  zwroot   (ie,je)      , & ! mean water content over root depth
  zropartw(ie,je,ke_soil)   ! fraction of layer filled by roots * w_fr

  REAL    (KIND=ireals   ) ::  &
!
! Thermal variables
!
  zts      (ie,je)      , & ! soil surface temperature
  ztsnow   (ie,je)      , & ! snow surface temperaure
  ztsnow_mult   (ie,je,0:ke_snow)      , & ! snow surface temperaure
! HEATCOND (soil moisture dependent heat conductivity of the soil)
  zalamtmp (ie,je)      , & ! heat conductivity
  zalam    (ie,je,ke_soil),&! heat conductivity
  zrocg    (ie,je)      , & ! volumetric heat capacity of bare soil
  zrocs    (ie,je)      , & ! heat capacity of snow
  ztsn     (ie,je)      , & ! new value of zts
  ztsnown  (ie,je)      , & ! new value of ztsnow
  ztsnown_mult  (ie,je,0:ke_snow)      , & ! new value of ztsnow
  znlgw1f  (ke_soil)    , & ! utility variable
  zaga(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zagb(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zagc(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zagd(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zage(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
!
! Auxiliary variables
!
  zw_snow_old(ie,je)    , & !
  zdqvtsnow(ie,je)      , & ! first derivative of saturation specific humidity
                            !    with respect to t_snow
  zts_pm   (ie,je)      , & ! indicator zts > < T_melt
  ztsnow_pm(ie,je)          ! indicator ztsnow > < T_melt

  LOGICAL     :: &
!
  ldebug                , & !
  limit_tch (ie,je)         ! indicator for flux limitation problem

  ! ground water as lower boundary of soil column
  REAL    (KIND=ireals) ::  &
    zdelta_sm, zdhydcond_dlwfr, zklw_fr_kso, zklw_fr_ksom1, zdlw_fr_kso 
  
  ! HEATCOND
  REAL    (KIND=ireals) ::          &
    hzalam   (ie,je,ke_soil+1),     & ! heat conductivity
    zthetas, zlamli, zlamsat, zlams, rsandf, zlamq, zlam0, zrhod, zlamdry,  &
    zsri, zKe, zthliq, zlamic

  INTEGER (KIND=iintegers) :: i_loc, j_loc, isub

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



!>JH
  prg_gsp=0._ireals ! graupel not implemented yet 
!<JH

! >JH NEW SECTION FOR DEFINITION OF LAND POINTS !
  llandmask=.FALSE.
  do ns=nsubs0,nsubs1
     !amf ---
     ! Prepare basic surface properties (for land-points only)
 
     DO   j = jstarts, jends
        DO i = istarts, iends
           IF(landmask(i,j,ns) > 0.5_ireals) THEN        ! for land-points only
              llandmask(i,j,ns) = .true.
           END IF
        ENDDO
     ENDDO
  END DO
! <JH




  ierror = 0
  yerror = '        '

! select timelevel and timestep for calculations
! In the new formulation a 2 timelevel scheme is used for the soil model.
! Therefore it is not necessary any more to distinguish between Leapfrog
! and Runge-Kutta scheme.
! But the first timestep must not be done with dt/2, as in the Leapfrog
! scheme
  nx  = nnow
  IF ( (ntstep == 0) .AND. (.NOT. l2tls) ) THEN
    ! use the original dt and not dt/2
    zdt = 2 * dt
  ELSE
    zdt = dt
  ENDIF

  ! time step for run-off computation
  zroffdt = zdt

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar
  jstarts = jstartpar
  jends   = jendpar

#ifdef NECSX
  CALL collapse(.TRUE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif

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


  !>JH WRITE(*,*) subsfrac(1,1,1),subsfrac(1,1,2),t_so(1,1,:,:,:),'cf_snow ',cf_snow


!!$  write(0,*) "SFC-DIAGNOSIS TERRA ",ke,dt,nsubs1
!!$  write(0,*) "zml_soil: ",czmls
!!$  write(0,*) "t", t(25,36,:,:)
!!$  write(0,*) "p0",p0(25,36,:)
!!$  write(0,*) "qv",qv(25,36,:,:)
!!$  write(0,*) "u",u(25,36,:,:)
!!$  write(0,*) "ps",ps(25,36,:)
!!$  write(0,*) "t_g",t_g(25,36,:,:)
!!$  write(0,*) "qv_s",qv_s(25,36,:,:)
!!$  write(0,*) "t_so",t_so(25,36,:,:,:)
!!$  write(0,*) "w_so",w_so(25,36,:,:,:)
!!$  write(0,*) "tch_t",tch(25,36,:)
!!$  write(0,*) "tcm_t",tcm(25,36,:)
!!$  write(0,*) " tfv_t",tfv(25,36,:)
!!$  write(0,*) "soiltyp_t",soiltyp_subs(25,36,:)
!!$  write(0,*) "plcov_t",  plcov(25,36,:)
!!$  write(0,*) "rootdp_t", rootdp(25,36,:)
!!$  write(0,*) "sai_t",   sai(25,36,:) 
!!$  write(0,*) "tai_t",   tai(25,36,:) 
!!$  write(0,*) "eai_t",   eai(25,36,:) 
!!$  write(0,*) "t_2m_t",  t_2m(25,36,:) 
!!$  write(0,*) "u_10m_t", u_10m(25,36,:)   
!!$  write(0,*) "v_10m_t", v_10m(25,36,:)    
!!$  write(0,*) "sobs_t",  sobs(25,36,:)    
!!$  write(0,*) "thbs_t",  thbs(25,36,:)     
!!$  write(0,*) "pabs_t",  pabs(25,36,:)     
!!$  write(0,*) "llandmask_t",llandmask(25,36,:)    



! >JH NEW SECTION FOR DEFINITION OF LAND POINTS !
do ns=nsubs0,nsubs1
!amf ---

! Prepare basic surface properties (for land-points only)
 
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF(llandmask(i,j)) THEN        ! for land-points only


!>JH
!      IF(llandmask(i,j,ns)) THEN        ! for land-points only


!subs        mstyp       = NINT(soiltyp_subs(i,j))        ! soil type
        mstyp       = soiltyp_subs(i,j,ns)        ! soil type
        IF(m_styp(i,j) < 9 ) llandmask(i,j,ns) = .true.
!>JH      ENDIF
    ENDDO
  ENDDO
END DO
!!!! <JH





!subs em
IF(itype_subs .EQ. 2) THEN    ! tiles
  IF (ntstep == 0) THEN
    DO ns = nsubs0+1, nsubs1, 2      ! 1 - mean, 2 - ocean, 3 - lake, 4 - no snow first
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j,ns)) THEN  ! for landpoints only  !check is not necessary

!    PRINT *,nsubs0,nsubs1,i,j,ns,'subsfrac(i,j,:)',subsfrac(i,j,ns-1),subsfrac(i,j,ns)

!            IF (w_snow(i,j,nx,ns).NE.w_snow(i,j,nx,ns-1)) &
!              & PRINT *,i,j,'w_snow(i,j,nx,ns).ne.w_snow(i,j,nx,ns-1)'
!            IF (subsfrac(i,j,ns-1)+subsfrac(i,j,ns).NE.1._ireals) &
!              & PRINT *,i,j,'subsfrac(i,j,ns-1)+subsfrac(i,j,ns).ne.1.'

!IF(MOD(i,340).eq.0 .and. MOD(j,186).eq.0 .and. w_snow(i,j,nx,ns).gt.0) &
!print *,i,j,'sub',subsfrac(i,j,ns-1),subsfrac(i,j,ns),w_snow(i,j,nx,ns),w_snow(i,j,nx,ns-1)
            zf_snow_old(i,j) = subsfrac(i,j,ns)/(subsfrac(i,j,ns-1)+subsfrac(i,j,ns))

            zwsnow(i,j) = w_snow(i,j,nx,ns)*zf_snow_old(i,j) + &
              w_snow(i,j,nx,ns-1)*(1._ireals - zf_snow_old(i,j))
            zf_snow(i,j) = MAX( 0.01_ireals, MIN(1.0_ireals,zwsnow(i,j)/cf_snow) ) &
              &            *zsf_heav(zwsnow(i,j) - zepsi)

            w_snow(i,j,nx,ns  ) = zwsnow(i,j)/MAX(zf_snow(i,j),zepsi)
            w_snow(i,j,nx,ns-1) = 0._ireals

            fact1 = subsfrac(i,j,ns-1)+subsfrac(i,j,ns)
            subsfrac(i,j,ns-1) = (1._ireals - zf_snow(i,j))*fact1
            subsfrac(i,j,ns  ) = zf_snow(i,j)              *fact1
!IF(subsfrac(i,j,ns-1).ne.0. .OR. subsfrac(i,j,ns-1).ne.1.) print *,i,j,'subsfrac(i,j,ns-1)',subsfrac(i,j,ns-1)
!IF(subsfrac(i,j,ns  ).ne.0. .OR. subsfrac(i,j,ns  ).ne.1.) print *,i,j,'subsfrac(i,j,ns  )',subsfrac(i,j,ns  )

!IF(MOD(i,340).eq.0 .and. MOD(j,186).eq.0 .and. w_snow(i,j,nx,ns).gt.0) &
!!IF(zf_snow(i,j).gt.0.1 .and. zf_snow(i,j).lt.0.9) &
!print *,i,j,'soil',subsfrac(i,j,ns-1),subsfrac(i,j,ns  ),zf_snow_old(i,j)

          END IF  ! land-points only
        END DO
      END DO
    END DO
!    DO ns = nsubs0+1, nsubs1, 2      ! 1 - mean, 2 - ocean, 3 - lake, 4 - no snow first
!      DO kso = 1,ke_soil
!        DO   j = jstarts, jends
!          DO i = istarts, iends
!            IF (llandmask(i,j,ns)) THEN  ! for landpoints only  !check is not necessary
!              fact1 = (MIN(1._ireals - zf_snow_old(i,j), 1._ireals - zf_snow(i,j)))/MAX((1._ireals - zf_snow(i,j)),zepsi)
!              fact2 = (MAX(zf_snow(i,j) - zf_snow_old(i,j), 0._ireals))/MAX(zf_snow(i,j),zepsi)

!IF(MOD(i,340).eq.0 .and. MOD(j,186).eq.0 .and. w_snow(i,j,nx,ns).gt.0) &
!print *,i,j,kso,'snow',w_snow(i,j,nx,ns-1),w_snow(i,j,nx,ns  ),zf_snow    (i,j),fact1,fact2

!              tmp1 = t_so    (i,j,kso,nx,ns-1)
!              tmp2 = w_so    (i,j,kso,nx,ns-1)
!              tmp3 = w_so_ice(i,j,kso,nx,ns-1)

!              t_so    (i,j,kso,nx,ns-1) = t_so(i,j,kso,nx,ns-1)*fact1 + t_so(i,j,kso,nx,ns)*(1._ireals - fact1)
!              w_so    (i,j,kso,nx,ns-1) = w_so(i,j,kso,nx,ns-1)*fact1 + w_so(i,j,kso,nx,ns)*(1._ireals - fact1)
!              w_so_ice(i,j,kso,nx,ns-1) = w_so_ice(i,j,kso,nx,ns-1)*fact1 + w_so_ice(i,j,kso,nx,ns)*(1._ireals - fact1)
  
!              t_so    (i,j,kso,nx,ns) = tmp1 * fact2 + t_so(i,j,kso,nx,ns)*(1._ireals - fact2)
!              w_so    (i,j,kso,nx,ns) = tmp2 * fact2 + w_so(i,j,kso,nx,ns)*(1._ireals - fact2)
!              w_so_ice(i,j,kso,nx,ns) = tmp3 * fact2 + w_so_ice(i,j,kso,nx,ns)*(1._ireals - fact2)

!IF(MOD(i,340).eq.0 .and. MOD(j,186).eq.0 .and. w_snow(i,j,nx,ns).gt.0) &
!print *,i,j,kso,'temp',t_so    (i,j,kso,nx,ns-1),t_so    (i,j,kso,nx,ns)
!            END IF  ! land-points only
!          END DO
!        END DO
!      END DO        ! soil layers
!    END DO
  END IF
END IF
!em>



!subs amf
!print *,'t_g    in TERRA, begin',t_g (340,186,nnow,:)
!print *,'t_snow in TERRA, begin',t_snow(340,186,nnow,:)
!print *,'w_snow in TERRA, begin',w_snow(340,186,nnow,:)
!print *,'t_s    in TERRA, begin',t_s (340,186,nnow,:)
!>JH CALL ij_local (388, 539, i_loc, j_loc, isub, ierror)
do ns=nsubs0,nsubs1
!amf ---

! Prepare basic surface properties (for land-points only)
 
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF(llandmask(i,j)) THEN        ! for land-points only


!>JH
!      IF(llandmask(i,j,ns)) THEN        ! for land-points only


!subs        mstyp       = NINT(soiltyp_subs(i,j))        ! soil type
        mstyp       = soiltyp_subs(i,j,ns)        ! soil type
        m_styp(i,j) = mstyp                     ! array for soil type
        zdw   (i,j)  = cdw0  (mstyp)
        zdw1  (i,j)  = cdw1  (mstyp)
        zkw   (i,j)  = ckw0  (mstyp)
        zkw1  (i,j)  = ckw1  (mstyp)
        zik2  (i,j)  = cik2  (mstyp)
        zporv(i,j)  = cporv(mstyp)              ! pore volume
        zpwp (i,j)  = cpwp (mstyp)              ! plant wilting point
        zadp (i,j)  = cadp (mstyp)              ! air dryness point
        zfcap(i,j)  = cfcap(mstyp)              ! field capacity
        zrock(i,j)  = crock(mstyp)              ! rock or ice indicator
        zrocg(i,j)  = crhoc(mstyp)              ! heat capacity
        zdlam(i,j)  = cala1(mstyp)-cala0(mstyp) ! heat conductivity parameter
!subs        zbwt(i,j)   = MAX(0.001_ireals,rootdp(i,j))! Artificial minimum value
        zbwt(i,j)   = MAX(0.001_ireals,rootdp(i,j,ns))! Artificial minimum value
                                                ! for root depth
        zroota(i,j) = 3._ireals/zbwt(i,j)       ! root density profile parameter (1/m)
                                                ! zroota=0. creates the original TERRA_LM
                                                ! version with constant root density
                                                ! for root depth
        ! New arrays for BATS-scheme
        zk0di(i,j)  = ck0di(mstyp)              !
        zbedi(i,j)  = cbedi(mstyp)              !
        ! Arrays for soil water freezing/melting
        zsandf(i,j)   = csandf(mstyp)
        zclayf(i,j)   = cclayf(mstyp)
        zpsis(i,j)    = -zpsi0 * 10**(1.88_ireals-0.013_ireals*zsandf(i,j))
        zb_por(i,j)   = 2.91_ireals + .159_ireals*zclayf(i,j)
        zedb(i,j)     = 1._ireals/zb_por(i,j)

! print*,i,j,zporv(i,j),zclayf(i,j),zpsis(i,j),zb_por(i,j)

!<JH
        IF(m_styp(i,j) < 9 ) llandmask(i,j,ns) = .true.
!>JH      ENDIF
    ENDDO
  ENDDO




  ! Set three-dimensional variables
  DO kso = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF(llandmask(i,j)) THEN        ! for land-points only
!subs          mstyp           = NINT(soiltyp_subs(i,j))        ! soil type
        IF(llandmask(i,j,ns)) THEN        ! for land-points only
          mstyp           = soiltyp_subs(i,j,ns)        ! soil type
          zalam(i,j,kso)  = cala0(mstyp)              ! heat conductivity parameter
        ENDIF
      ENDDO
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

! For ntstep=0 : Some preparations
! ================================

  IF (ntstep == 0) THEN
!subs    w_so(:,:,:,nnew) = w_so(:,:,:,nx)
    w_so(:,:,:,nnew,ns) = w_so(:,:,:,nx,ns)

!   Provide for a soil moisture 1 % above air dryness point, reset soil
!   moisture to zero in case of ice and rock
    DO kso   = 1,ke_soil+1
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN             ! for land-points only
          IF (llandmask(i,j,ns)) THEN             ! for land-points only
            IF (m_styp(i,j).ge.3) THEN
!subs              w_so (i,j,kso,:) = MAX(w_so(i,j,kso,:),                     &
              w_so (i,j,kso,:,ns) = MAX(w_so(i,j,kso,:,ns),                     &
                                           1.01_ireals*zadp(i,j)*zdzhs(kso) )
            ELSE
!subs              w_so(i,j,kso,:) = 0.0_ireals
              w_so(i,j,kso,:,ns) = 0.0_ireals
            ENDIF

          END IF   ! land-points
        END DO
      END DO
    END DO

  !>JH WRITE(*,*) 'W_SO: ',w_so (1,1,:,:,:),'zadp: ',zadp(1,1),'zdzhs: ',zdzhs(:)

!   adjust temperature profile in lower soil layers, if temperature of first soil
!   layer was reduced due to the detection of snow (e.g. in the analysis)
!   loop over grid points
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF(llandmask(i,j)) THEN   ! for land-points only
!subs          IF (w_snow(i,j,nx) <=1.0E-6_ireals) THEN
        IF(llandmask(i,j,ns)) THEN   ! for land-points only
          IF (w_snow(i,j,nx,ns) <=1.0E-6_ireals) THEN
            ! spurious snow is removed
!subs            t_snow(i,j,:)=t_so(i,j,0,nx)
!subs            w_snow(i,j,nx)= 0.0_ireals
            t_snow(i,j,:,ns)=t_so(i,j,0,nx,ns)
            w_snow(i,j,nx,ns)= 0.0_ireals
            IF(lmulti_snow) THEN
!subs              t_snow_mult(i,j,0,nx) = t_so(i,j,0,nx)
              t_snow_mult(i,j,0,nx,ns) = t_so(i,j,0,nx,ns)
              DO ksn = 1, ke_snow
!subs                t_snow_mult(i,j,ksn,nx) = t_so(i,j,0,nx)
!subs                wliq_snow(i,j,ksn,nx) = 0.0_ireals
!subs                wtot_snow(i,j,ksn,nx) = 0.0_ireals
!subs                rho_snow_mult(i,j,ksn,nx) = 0.0_ireals
!subs                dzh_snow (i,j,ksn,nx) = 0.0_ireals
                t_snow_mult(i,j,ksn,nx,ns) = t_so(i,j,0,nx,ns)
                wliq_snow(i,j,ksn,nx,ns) = 0.0_ireals
                wtot_snow(i,j,ksn,nx,ns) = 0.0_ireals
                rho_snow_mult(i,j,ksn,nx,ns) = 0.0_ireals
                dzh_snow (i,j,ksn,nx,ns) = 0.0_ireals
              END DO
            END IF
          ELSE
!           adjust soil temperatures in presence of snow
!subs            zdel_t_so = MAX (0._ireals,t_so(i,j,0,nx)-t0_melt)
!subs            t_so(i,j,0,:)=MIN( t0_melt,t_so(i,j,0,nx) )
            zdel_t_so = MAX (0._ireals,t_so(i,j,0,nx,ns)-t0_melt)
            t_so(i,j,0,:,ns)=MIN( t0_melt,t_so(i,j,0,nx,ns) )
            DO kso=1,ke_soil
!subs              IF ( t_so(i,j,kso,nx) > t0_melt) THEN
!subs                 t_so(i,j,kso,:) = MAX( t0_melt,                  &
!subs                                   t_so(i,j,kso,nx) -zdel_t_so    &
              IF ( t_so(i,j,kso,nx,ns) > t0_melt) THEN
                 t_so(i,j,kso,:,ns) = MAX( t0_melt,                  &
                                   t_so(i,j,kso,nx,ns) -zdel_t_so    &
                                   *(zmls(ke_soil+1)-zmls(kso))   &
                                   /(zmls(ke_soil+1)-zmls( 1 )) )
              ENDIF
!subs              IF ( t_so(i,j,kso,nx) <= t0_melt) EXIT
              IF ( t_so(i,j,kso,nx,ns) <= t0_melt) EXIT
            ENDDO
          ENDIF

!subs          t_s(i,j,:) = t_so(i,j,0,:)
          t_s(i,j,:,ns) = t_so(i,j,0,:,ns)

!         Set level 1 to level 0 for t_so for every landpoint
!subs          t_so(i,j,1,:) = t_so(i,j,0,nx)
          t_so(i,j,1,:,ns) = t_so(i,j,0,nx,ns)
        ENDIF    ! llandmask
      END DO
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
        DO   j = jstarts, jends
          DO i = istarts, iends
!subs            IF (llandmask(i,j)) THEN             ! for land-points only
!subs              w_so_ice(i,j,kso,:) = 0.0_ireals
!subs              IF (t_so(i,j,kso,nx) < t0_melt) THEN
!subs                w_so_ice(i,j,kso,:) = w_so(i,j,kso,nx)
            IF (llandmask(i,j,ns)) THEN             ! for land-points only
              w_so_ice(i,j,kso,:,ns) = 0.0_ireals
              IF (t_so(i,j,kso,nx,ns) < t0_melt) THEN
                w_so_ice(i,j,kso,:,ns) = w_so(i,j,kso,nx,ns)
              END IF ! t_so(kso) < t0_melt
            END IF   ! land-points
          END DO
        END DO
      END DO
    END IF           ! lmelt .AND. .NOT. lmelt_var
    IF(lmelt .AND. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
        DO   j = jstarts, jends
          DO i = istarts, iends
!subs            IF (llandmask(i,j)) THEN             ! for land-points only
!subs              IF (t_so(i,j,kso,nx) < (t0_melt-zepsi)) THEN
            IF (llandmask(i,j,ns)) THEN             ! for land-points only
              IF (t_so(i,j,kso,nx,ns) < (t0_melt-zepsi)) THEN
                zaa    = g*zpsis(i,j)/lh_f

! print*,i,j,zaa,zedb(i,j),zporv(i,j),t_so(i,j,0:4,1,:)

                zw_m(i,j)     = zporv(i,j)*zdzhs(kso)
                zw_m(i,j)   = zw_m(i,j)*                                          &
!subs                  ((t_so(i,j,kso,nx) - t0_melt)/(t_so(i,j,kso,nx)*zaa))**(-zedb(i,j))
!subs                w_so_ice(i,j,kso,:) = MAX (0.0_ireals,w_so(i,j,kso,nx) - zw_m(i,j))
                  ((t_so(i,j,kso,nx,ns) - t0_melt)/(t_so(i,j,kso,nx,ns)*zaa))**(-zedb(i,j))
                w_so_ice(i,j,kso,:,ns) = MAX (0.0_ireals,w_so(i,j,kso,nx,ns) - zw_m(i,j))
              END IF
            END IF   ! land-points
          END DO
        END DO
      END DO
    END IF           ! lmelt .AND. lmelt_var

!   Initialization of snow density, if necessary
!   --------------------------------------------
    IF (.NOT. lana_rho_snow) THEN
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN             ! for land-points only
          IF (llandmask(i,j,ns)) THEN             ! for land-points only
            IF(lmulti_snow) THEN
              DO ksn = 1, ke_snow
!subs                rho_snow_mult(i,j,ksn,nx) = 250.0_ireals    ! average initial density
!subs                wtot_snow    (i,j,ksn,nx) = w_snow(i,j,nx)/ke_snow
!subs                dzh_snow     (i,j,ksn,nx) = w_snow(i,j,nx)*rho_w/rho_snow_mult(i,j,ksn,nx)/ke_snow
!subs                wliq_snow    (i,j,ksn,nx) = 0.0_ireals
                rho_snow_mult(i,j,ksn,nx,ns) = 250.0_ireals    ! average initial density
                wtot_snow    (i,j,ksn,nx,ns) = w_snow(i,j,nx,ns)/REAL(ke_snow,ireals)
                dzh_snow     (i,j,ksn,nx,ns) = w_snow(i,j,nx,ns)/REAL(ke_snow,ireals) &
                  &                            *rho_w/rho_snow_mult(i,j,ksn,nx,ns)
                wliq_snow    (i,j,ksn,nx,ns) = 0.0_ireals
              END DO
            ELSE
!subs              rho_snow(i,j,nx) = 250.0_ireals    ! average initial density
              rho_snow(i,j,nx,ns) = 250.0_ireals    ! average initial density
            ENDIF
          ENDIF   ! land-points
        ENDDO
      ENDDO
    ELSE
      IF(lmulti_snow) THEN
        zw_snow_old(:,:) = 0.0_ireals
        DO ksn = 1, ke_snow
          DO   j = jstarts, jends
            DO i = istarts, iends
!subs              IF (llandmask(i,j)) THEN             ! for land-points only
              IF (llandmask(i,j,ns)) THEN             ! for land-points only
!subs                zw_snow_old(i,j) = zw_snow_old(i,j) + wtot_snow(i,j,ksn,nx)
                zw_snow_old(i,j) = zw_snow_old(i,j) + wtot_snow(i,j,ksn,nx,ns)
              END IF
            END DO
          END DO
        END DO
        DO ksn = 1, ke_snow
          DO   j = jstarts, jends
            DO i = istarts, iends
!subs              IF (llandmask(i,j)) THEN             ! for land-points only
!subs                IF(dzh_snow(i,j,ksn,nx) .LT. zepsi  .OR. ABS(zw_snow_old(i,j)-w_snow(i,j,nx)).GE.zepsi) THEN
!subs                  IF(rho_snow(i,j,nx) .EQ. 0._ireals) rho_snow(i,j,nx) = 250._ireals
!subs                  rho_snow_mult(i,j,ksn,nx) = rho_snow(i,j,nx)
!subs                  wtot_snow    (i,j,ksn,nx) = w_snow(i,j,nx)/ke_snow
!subs                  dzh_snow     (i,j,ksn,nx) = w_snow(i,j,nx)*rho_w/rho_snow_mult(i,j,ksn,nx)/ke_snow
!subs                  wliq_snow    (i,j,ksn,nx) = 0.0_ireals
              IF (llandmask(i,j,ns)) THEN             ! for land-points only
                IF(dzh_snow(i,j,ksn,nx,ns) .LT. zepsi  &
                  & .OR. ABS(zw_snow_old(i,j)-w_snow(i,j,nx,ns)).GE.zepsi) THEN
                  IF(rho_snow(i,j,nx,ns) .EQ. 0._ireals) rho_snow(i,j,nx,ns) = 250._ireals
                  rho_snow_mult(i,j,ksn,nx,ns) = rho_snow(i,j,nx,ns)
                  wtot_snow    (i,j,ksn,nx,ns) = w_snow(i,j,nx,ns)/REAL(ke_snow,ireals)
                  dzh_snow     (i,j,ksn,nx,ns) = w_snow(i,j,nx,ns)/REAL(ke_snow,ireals) &
                    &                            *rho_w/rho_snow_mult(i,j,ksn,nx,ns)
                  wliq_snow    (i,j,ksn,nx,ns) = 0.0_ireals
                ELSE
!subs             rho_snow_mult(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx)/MAX(dzh_snow(i,j,ksn,nx),0.0_ireals)*rho_w    ! initial density
                  rho_snow_mult(i,j,ksn,nx,ns) = wtot_snow(i,j,ksn,nx,ns)                 &
                                               & /MAX(dzh_snow(i,j,ksn,nx,ns),0.0_ireals) &
                                               & *rho_w  ! initial density
                END IF
              END IF
            END DO
          END DO
        END DO
      ELSE
        DO   j = jstarts, jends
          DO i = istarts, iends
!subs            IF (llandmask(i,j)) THEN             ! for land-points only
!subs              IF(rho_snow(i,j,nx) .EQ. 0._ireals) rho_snow(i,j,nx) = 250._ireals
            IF (llandmask(i,j,ns)) THEN             ! for land-points only
              IF(rho_snow(i,j,nx,ns) .EQ. 0._ireals) rho_snow(i,j,nx,ns) = 250._ireals
            END IF
          END DO
        END DO
      END IF
    ENDIF

!   Initialization of the local array of the grid in snow
!   --------------------------------------------
    IF(lmulti_snow) THEN
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF(llandmask(i,j)) THEN   ! for land-points only
          IF(llandmask(i,j,ns)) THEN   ! for land-points only
!subs            h_snow(i,j,nnew) = 0.0_ireals
            h_snow(i,j,nnew,ns) = 0.0_ireals
            DO ksn = 1,ke_snow
!subs              zdzh_snow(i,j,ksn) = dzh_snow(i,j,ksn,nx)
!subs              h_snow(i,j,nnew) = h_snow(i,j,nnew) + zdzh_snow(i,j,ksn)
              zdzh_snow(i,j,ksn) = dzh_snow(i,j,ksn,nx,ns)
              h_snow(i,j,nnew,ns) = h_snow(i,j,nnew,ns) + zdzh_snow(i,j,ksn)
            END DO
!subs            h_snow(i,j,nnow) = h_snow(i,j,nnew)
            h_snow(i,j,nnow,ns) = h_snow(i,j,nnew,ns)

!subs            zhh_snow(i,j,1) = - h_snow(i,j,nnow) + dzh_snow(i,j,1,nx)
            zhh_snow(i,j,1) = - h_snow(i,j,nnow,ns) + dzh_snow(i,j,1,nx,ns)
            DO ksn = 2,ke_snow
!subs              zhh_snow(i,j,ksn) = zhh_snow(i,j,ksn-1) + dzh_snow(i,j,ksn,nx)
              zhh_snow(i,j,ksn) = zhh_snow(i,j,ksn-1) + dzh_snow(i,j,ksn,nx,ns)
            END DO

!subs            zhm_snow(i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals
            zhm_snow(i,j,1) = (-h_snow(i,j,nnow,ns) + zhh_snow(i,j,1))/2._ireals
            DO ksn = 2,ke_snow
              zhm_snow(i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals
            END DO
          ENDIF    ! llandmask
        END DO
      END DO
    END IF

  ENDIF             ! ntstep = 0

! End of timestep 0 preparations
! ==============================




! To ensure t_s = t_so(1) also at gme2lm-modified lateral boundaries by
! modifying the soil temperature profile:
  DO kso   = 1,ke_soil+1
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN             ! for land-points only
        IF (llandmask(i,j,ns)) THEN             ! for land-points only
        ! Next lines have to be modified/deleted if the soil surface temperature
        ! t_so(i,j,0,nx) predicted by the heat conduction equation is used
!subs          zdelt_bound = t_s(i,j,nx) - t_so(i,j,1,nx)
!subs          t_so(i,j,kso,nx) = t_so(i,j,kso,nx) + zdelt_bound * (1._ireals - &
!subs                           (zmls(kso) - zmls(1))/(zmls(ke_soil+1) - zmls(1)))
!subs          IF(kso.EQ.1) t_so(i,j,0,nx) = t_so(i,j,1,nx)
          zdelt_bound = t_s(i,j,nx,ns) - t_so(i,j,1,nx,ns)
          t_so(i,j,kso,nx,ns) = t_so(i,j,kso,nx,ns) + zdelt_bound * (1._ireals - &
                           (zmls(kso) - zmls(1))/(zmls(ke_soil+1) - zmls(1)))
          IF(kso.EQ.1) t_so(i,j,0,nx,ns) = t_so(i,j,1,nx,ns)
        ENDIF
      END DO
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
        ztrangs(i,j)      = 0.0_ireals
        zsnull(i,j)       = 0.0_ireals
        zs1(i,j)          = 0.0_ireals
        zwroot(i,j)       = 0.0_ireals
        zdwidt (i,j)      = 0.0_ireals         ! Initialisation of all
        zdwsndt(i,j)      = 0.0_ireals         ! evaporation quantities
        zesoil (i,j)      = 0.0_ireals         ! to be equal to zero
        zrr    (i,j)      = 0.0_ireals         ! in first part formation of dew
        zrs    (i,j)      = 0.0_ireals         ! in first part formation of rime
!subs        zw_fr(i,j,ke_soil+1)  = w_so(i,j,ke_soil+1,nx)/zdzhs(ke_soil+1)
        zw_fr(i,j,ke_soil+1)  = w_so(i,j,ke_soil+1,nx,ns)/zdzhs(ke_soil+1)
    END DO
  END DO


  DO kso   = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        zw_fr(i,j,kso)    = w_so(i,j,kso,nx)/zdzhs(kso)
        zw_fr(i,j,kso)    = w_so(i,j,kso,nx,ns)/zdzhs(kso)
        ztrang (i,j,kso)  = 0._ireals
        zropartw(i,j,kso) = 0._ireals
      END DO
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
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN   ! for land-points only
!subs          IF (kso.le.k10cm) zs1(i,j) = zs1(i,j) + w_so(i,j,kso,nx)
!subs          zsnull(i,j)   = zsnull(i,j) + w_so(i,j,kso,nx)
        IF (llandmask(i,j,ns)) THEN   ! for land-points only
          IF (kso.le.k10cm) zs1(i,j) = zs1(i,j) + w_so(i,j,kso,nx,ns)
          zsnull(i,j)   = zsnull(i,j) + w_so(i,j,kso,nx,ns)
        END IF
      END DO
    END DO
  END DO

IF (itype_heatcond == 2) THEN

! heat conductivity dependent on actual soil water content
! based on Peters-Lidard et al. (1998) and Johansen (1975)
! see also Block, Alexander (2007), Dissertation BTU Cottbus
    
  DO kso = 1, ke_soil+1
   DO j = jstarts, jends
    DO i = istarts, iends
!subs     IF (llandmask(i,j)) THEN          ! land-points only
     IF (llandmask(i,j,ns)) THEN          ! land-points only
      zthetas = zporv(i,j)

      zlamli  = 0.57_ireals
      zlamic  = 2.2_ireals
!subs      zthliq  = zthetas - (w_so_ice(i,j,kso,nx)/zdzhs(kso))
      zthliq  = zthetas - (w_so_ice(i,j,kso,nx,ns)/zdzhs(kso))
      zlamq   = 7.7_ireals

      rsandf=zsandf(i,j)/100._ireals

      if (rsandf >= 0.2_ireals)  zlam0    =  2.0_ireals
      if (rsandf < 0.2_ireals)   zlam0    =  3.0_ireals


      zlams   = zlamq**rsandf* zlam0**(1._ireals-rsandf)

!      zlamsat = (clamso(m_styp(i,j))**(1.0_ireals-zthetas)) * (zlamli**zthliq) * (zlamic**(zthetas-zthliq))
      zlamsat = (zlams**(1.0_ireals-zthetas)) * (zlamli**zthliq) * (zlamic**(zthetas-zthliq))

      zrhod   = 2700.0_ireals*(1.0_ireals-zthetas)
      zlamdry = ( 0.135_ireals*zrhod + 64.7_ireals )                     &
              / ( 2700.0_ireals - 0.947*zrhod )

!subs      zsri    = MIN(1.0_ireals, (w_so(i,j,kso,nx)/zdzhs(kso)) / zthetas)
      zsri    = MIN(1.0_ireals, (w_so(i,j,kso,nx,ns)/zdzhs(kso)) / zthetas)

      IF ( zsri >= 0.1_ireals ) THEN
!subs       IF ( t_so(i,j,kso,nx) < t0_melt) THEN
       IF ( t_so(i,j,kso,nx,ns) < t0_melt) THEN
        zKe = zsri
       ELSE
        zKe = LOG10(zsri) + 1.0_ireals
       ENDIF
       zKe = MAX(0.0_ireals, (MIN(1.0_ireals, zKe)) )
       hzalam(i,j,kso) = zKe*(zlamsat - zlamdry) + zlamdry
      ELSE
       hzalam(i,j,kso) = zlamdry
      ENDIF
     END IF  ! land-points
    ENDDO
   ENDDO
  ENDDO

  DO kso = 1, ke_soil
   DO j = jstarts, jends
    DO i = istarts, iends
!subs     IF (llandmask(i,j)) THEN          ! land-points only
     IF (llandmask(i,j,ns)) THEN          ! land-points only
      ! mean heat conductivity
      zalam(i,j,kso) = hzalam(i,j,kso)*hzalam(i,j,kso+1)*(zmls(kso+1)-zmls(kso))    &
                     / ( hzalam(i,j,kso)*(zmls(kso+1)-zzhls(kso))                   &
                     +   hzalam(i,j,kso+1)*(zzhls(kso)-zmls(kso)) )
     END IF  ! land-points
    ENDDO
   ENDDO
  ENDDO

ELSE

! heat conductivity based on assumption of a soil water content which is equal to the
! average between wilting point and field capacity
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        zwqg         = 0.5_ireals*(zfcap(i,j) + zpwp(i,j))
        z4wdpv       = 4._ireals*zwqg/zporv(i,j)
        ! heat conductivity
        zalamtmp(i,j) =              zdlam(i,j)                         &
                      * (0.25_ireals + 0.30_ireals*zdlam(i,j)           &
                      / (1._ireals+0.75_ireals*zdlam(i,j)))             &
                      * MIN (z4wdpv, 1.0_ireals + (z4wdpv-1.0_ireals)   &
                      *(1.0_ireals+0.35_ireals*zdlam(i,j))              &
                      /(1.0_ireals+1.95_ireals*zdlam(i,j)))
      END IF  ! land-points
    ENDDO
  ENDDO
  DO kso = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN          ! land-points only
        IF (llandmask(i,j,ns)) THEN          ! land-points only
          zalam(i,j,kso) = zalam(i,j,kso) + zalamtmp(i,j)
        END IF  ! land-points
      ENDDO
    ENDDO
  ENDDO

ENDIF

#ifdef NECSX
  CALL collapse(.FALSE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif
! Initialisations and conversion of tch to tmch

  limit_tch(:,:) = .false.  !  preset problem indicator
  ldebug         = .false.

  DO   j = jstarts, jends
    jm1  = MAX( 1, j-1 )
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN     ! for land-points only
      IF (llandmask(i,j,ns)) THEN     ! for land-points only
        im1        = MAX( 1, i-1)
        zuv        = 0.5_ireals*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                               +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
!subs        ztvs       = t_g (i,j,nx)*(1.0_ireals + rvd_m_o*qv_s(i,j,nx))
        ztvs       = t_g (i,j,nx,ns)*(1.0_ireals + rvd_m_o*qv_s(i,j,nx,ns))
!JH IF(i.eq.1 .and. j.eq.1)  print *,ns,'ztvs',ztvs,t_g (i,j,nx,ns),rvd_m_o,qv_s(i,j,nx,ns),zuv,ps(i,j,nx),t(i,j,ke,nx)
       
        !  'potential' temperature of lowest atmospheric layer
        zplow          = p0(i,j,ke) + pp(i,j,ke,nx)
        zth_low (i,j)  =  t(i,j,ke,nx) *( (ps(i,j,nx)/zplow)**rdocp )
!>JH IF(zth_low (i,j).ge.350.) print *,ns,'zth_low (i,j).ge.350.',zth_low (i,j),i_global(i),j_global(j)
!subs        zdt_atm (i,j)  =  zth_low(i,j)-t_g(i,j,nx)
!subs        zdq_atm (i,j)  =  qv(i,j,ke,nx)-qv_s(i,j,nx)
        zdt_atm (i,j)  =  zth_low(i,j)-t_g(i,j,nx,ns)
        zdq_atm (i,j)  =  qv(i,j,ke,nx)-qv_s(i,j,nx,ns)

!   introduce an artificical upper boundary on transfer coefficients in cases
!   where an extreme cooling/heating of topmost soil layer may be caused due
!   to excessive sensible/latent heat fluxes (e.g. after creation of unbalanced
!   structures in the data assimilation or following massive changes in radiative
!   forcing due to infrequent radiation computations)

!   estimate current energy budget of topmost soil layer

!   heat flux between layers 1&2 based on current temperature profile
!subs    zg1(i,j)= zalam(i,j,1)*(t_so(i,j,1,nx)-t_so(i,j,2,nx))/zdzms(2)
    zg1(i,j)= zalam(i,j,1)*(t_so(i,j,1,nx,ns)-t_so(i,j,2,nx,ns))/zdzms(2)





  if (i.eq.40 .and. j.eq.1) then
print*, "SFC-DIAGNOSIS TERRA ",ke,dt,nsubs1,ntstep
print*," nztlev ",               nztlev   
print*," lmelt  ",               lmelt    
print*," lmelt_var ",            lmelt_var
print*," lmulti_snow ",          lmulti_snow 
print*," itype_gscp ",           itype_gscp
print*," itype_trvg ",           itype_trvg
print*," itype_evsl ",           itype_evsl
print*," itype_tran ",           itype_tran
print*," itype_root ",           itype_root
print*," itype_heatcond ",       itype_heatcond
print*," itype_hydbound ",       itype_hydbound
print*," lstomata  ",            lstomata
print*," l2tls ",                l2tls  
print*," lana_rho_snow ",        lana_rho_snow
print*," itype_subs ",           itype_subs
  print*, "zml_soil: ",czmls
  print*, "t", t(i,j,ke,:)
  print*, "p0",p0(i,j,ke)
  print*, "qv",qv(i,j,ke,:)
  print*, "u",u(i,j,ke,:)
  print*, "ps",ps(i,j,:)
  print*, "t_g",t_g(i,j,:,:)
  print*, "qv_s",qv_s(i,j,:,:)
  print*, "t_so",t_so(i,j,:,:,:)
  print*, "w_so",w_so(i,j,:,:,:)
  print*, "tch_t",tch(i,j,:)
  print*, "tcm_t",tcm(i,j,:)
  print*, " tfv_t",tfv(i,j,:)
  print*, "soiltyp_t",soiltyp_subs(i,j,:)
  print*, "plcov_t",  plcov(i,j,:)
  print*, "rootdp_t", rootdp(i,j,:)
  print*, "sai_t",   sai(i,j,:) 
  print*, "tai_t",   tai(i,j,:) 
  print*, "eai_t",   eai(i,j,:) 
  print*, "t_2m_t",  t_2m(i,j,:) 
  print*, "u_10m_t", u_10m(i,j,:)   
  print*, "v_10m_t", v_10m(i,j,:)    
  print*, "sobs_t",  sobs(i,j,:)    
  print*, "thbs_t",  thbs(i,j,:)     
  print*, "pabs_t",  pabs(i,j,:)     
  print*, "llandmask_t",llandmask(i,j,:) 

     end if

!   estimates of sensible and latent heat flux
    zrho_atm(i,j)=ps(i,j,nx)/(r_d*ztvs)
!subs    zshfl(i,j) = tch(i,j)*zuv*zrho_atm(i,j)*cp_d*zdt_atm(i,j)
!subs    zlhfl(i,j) = tch(i,j)*zuv*zrho_atm(i,j)*lh_v*zdq_atm(i,j)
    zshfl(i,j) = tch(i,j,ns)*zuv*zrho_atm(i,j)*cp_d*zdt_atm(i,j)
    zlhfl(i,j) = tch(i,j,ns)*zuv*zrho_atm(i,j)*lh_v*zdq_atm(i,j)

    IF (ABS(zshfl(i,j)) <= zepsi) zshfl(i,j)=SIGN(zepsi,zshfl(i,j))

    zthfl(i,j) = zshfl(i,j)+zlhfl(i,j)

!   net radiative fluxes at surface
!subs    zradfl(i,j) = sobs(i,j)+thbs(i,j)
    zradfl(i,j) = sobs(i,j,ns)+thbs(i,j,ns)

!   unconstrained estimated energy budget of topmost soil layer
    zeb1(i,j) = zshfl(i,j)+zlhfl(i,j)+zradfl(i,j)-zg1(i,j)

!   energy required to melt existing snow
!subs    ze_melt(i,j)=w_snow(i,j,nx)*rho_w*lh_f     ! (J/m**2)
    ze_melt(i,j)=w_snow(i,j,nx,ns)*rho_w*lh_f     ! (J/m**2)
!   heat capacity of snow layer
!subs    zch_snow(i,j)=w_snow(i,j,nx)*rho_w*chc_i   ! (J/(m**2 K))
    zch_snow(i,j)=w_snow(i,j,nx,ns)*rho_w*chc_i   ! (J/(m**2 K))

!   constrain transfer coefficient, if energy budget  of topmost soil layer is:
!   a) negative & surface layer is unstable (i.e   upward directed turbulent heat flux)
!   b) positive & surface layer is stable   (i.e downward directed turbulent heat flux)

  IF (zeb1(i,j)<0.0_ireals .AND. zthfl(i,j)<0.0_ireals) THEN 
    ! cooling of 1st soil layer&upward SHF+LHF

    ztchv_max(i,j) =  &
               (-zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i,j))/zdt            &
                +zg1(i,j)-zradfl(i,j) )                                     &
!subs               / zthfl(i,j) * tch(i,j) * zuv
               / zthfl(i,j) * tch(i,j,ns) * zuv
  ELSEIF (zeb1(i,j)>0.0_ireals .AND. zthfl(i,j)>0.0_ireals) THEN 
    ! heating of 1st soil layer & downward SHF+LHF
    !   Note: The heat capacity of snow is only relevant for the difference 
    !         between the actual temperature and the melting point. The mean 
    !         snow temperature is set to the average of t_snow & t_so(1)
    IF(lmulti_snow) THEN
!subs      zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow_mult(i,j,1,nx)+t_so(i,j,1,nx)))
      zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow_mult(i,j,1,nx,ns)+t_so(i,j,1,nx,ns)))
    ELSE
!subs      zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow(i,j,nx)+t_so(i,j,1,nx)))
      zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow(i,j,nx,ns)+t_so(i,j,1,nx,ns)))
    ENDIF
    ztchv_max(i,j) =  &
    ( (zlim_dtdt*zch_soil*zdzhs(1)+zdT_snow*zch_snow(i,j)+ze_melt(i,j))/zdt  &
      + zg1(i,j)-zradfl(i,j) )                                               &
!subs                / zthfl(i,j) * tch(i,j) * zuv
                / zthfl(i,j) * tch(i,j,ns) * zuv
  ELSE                                                    
    ! unlimited transfer coefficient
    ztchv_max(i,j) = HUGE(1._ireals)
  ENDIF

                                                  ! required constraint as non-turbulent
                                                  ! energy budget components alone may
  ztchv_max(i,j) =MAX( ztchv_max(i,j) ,zepsi)     ! exceed the upper limit in the energy
                                                  ! budget of the 1st soil layer

!subs  ztchv(i,j)    = tch(i,j)*zuv  ! transfer coefficient * velocity
  ztchv(i,j)    = tch(i,j,ns)*zuv  ! transfer coefficient * velocity

        LIM: IF (ztchv(i,j) > ztchv_max(i,j)) THEN
!subs          tch(i,j)=ztchv_max(i,j)/MAX(zuv,1.E-06_ireals)
          tch(i,j,ns)=ztchv_max(i,j)/MAX(zuv,1.E-06_ireals)
!         IF (ntstep > 10) THEN          ! control print only after initial adaptation
!                                        ! to analysis state
            limit_tch(i,j) = .true.      ! set switch for later use
        END IF LIM

!subs        ztmch(i,j) = tch(i,j)*zuv*g*zrho_atm(i,j)
        ztmch(i,j) = tch(i,j,ns)*zuv*g*zrho_atm(i,j)
      ENDIF
    ENDDO
  ENDDO

! counter for limitation of transfer coefficients
  m_limit = COUNT( limit_tch(:,:) ) 
!*DL: Reduce dayfile printout
! IF (m_limit > 0) THEN
!   WRITE(*,'(1X,A,I8,A,I12)')                        &
!   ' time step:',ntstep,                             &
!   ' no. of transfer coefficient limitations:',m_limit
! ENDIF


! In debugging mode and if transfer coefficient occured for at least one grid point
  IF (m_limit > 0 .AND. ldebug) THEN
    WRITE(*,'(1X,A,/,2(1X,A,F10.2,A,/),1X,A,F10.2,/,1X,A,F10.3,/)')                  &
           'terra1: transfer coefficient had to be constrained',                  &
           'model time step                                 :', zdt     ,' seconds', &
           'max. temperature increment allowed per time step:',zlim_dtdt,' K',       &
           'upper soil model layer thickness                :', zdzhs(1)
    DO j = jstarts, jends
    jm1  = MAX( 1, j-1 )
      DO i = istarts, iends
      im1 = MAX(1,i-1)
      IF (limit_tch(i,j)) THEN
        zuv        = 0.5_ireals*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                               +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
        yhc       ='COOLING'
        IF (zeb1(i,j) > 0._ireals) Yhc='HEATING'
!!$>JH        PRINT 7001,ntstep,i,j,my_cart_id,Yhc,ps(i,j,nx),  &
!!$              t(i,j,ke,nx), zrho_atm(i,j),     &
!!$!subs              w_snow(i,j,nx),zch_snow(i,j),ze_melt(i,j), &
!!$              w_snow(i,j,nx,ns),zch_snow(i,j),ze_melt(i,j), &
!!$              zth_low(i,j),zdt_atm(i,j),       &
!!$!subs              t_g   (i,j,nx),t_so(i,j,2,nx),t_so(i,j,3,nx), &
!!$!subs              qv(i,j,ke,nx)*1000._ireals,qv_s(i,j,nx)*1000._ireals, &
!!$!subs              w_so(i,j,1,nx)/zdzhs(1),w_so(i,j,2,nx)/zdzhs(2), &
!!$!subs              w_so(i,j,3,nx)/zdzhs(3),                               &
!!$              t_g   (i,j,nx,ns),t_so(i,j,2,nx,ns),t_so(i,j,3,nx,ns), &
!!$              qv(i,j,ke,nx)*1000._ireals,qv_s(i,j,nx,ns)*1000._ireals, &
!!$              w_so(i,j,1,nx,ns)/zdzhs(1),w_so(i,j,2,nx,ns)/zdzhs(2), &
!!$              w_so(i,j,3,nx,ns)/zdzhs(3),                               &
!!$              u(i,j,ke,nx),v(i,j,ke,nx),zuv, &
!!$!subs              sobs(i,j),thbs(i,j),   &
!!$              sobs(i,j,ns),thbs(i,j,ns),   &
!!$              zshfl(i,j),zlhfl(i,j),           &
!!$              zthfl(i,j),zg1(i,j),zeb1(i,j), &
!!$              ztchv(i,j),ztchv_max(i,j)
!!$     7001 FORMAT(1x,'time step:',I5,'   grid point:',3I4,2X,A8,/,&
!!$                 1X,'surface pressure   :',F12.2,' Pa'     ,/, &
!!$                 1X,'t(ke)              :',F12.2,' K '     ,/, &
!!$                 1X,'air density        :',F12.2,' kg/m**3',/, &
!!$                 1X,'snow water content :',F12.4,' mH2O   ',/, &
!!$                 1X,'snow heat capacity :',F12.2,' J/(m**2 K)',/, &
!!$                 1X,'snow melting energy:',F12.2,' J/m**2   ',/, &
!!$                 1X,'th_low             :',F12.2,' K '     ,/, &
!!$                 1X,'th_low-t_g         :',F12.2,' K '     ,/, &
!!$                 1X,'t_g                :',F12.2,' K '     ,/, &
!!$                 1X,'t_so(2)            :',F12.2,' K '     ,/, &
!!$                 1X,'t_so(3)            :',F12.2,' K '     ,/, &
!!$                 1X,'qv(ke)             :',F12.4,' g/kg'   ,/, &
!!$                 1X,'qv_s               :',F12.4,' g/kg'   ,/, &
!!$                 1X,'zwg_fr(1)          :',F12.2,' - '     ,/, &
!!$                 1X,'zwg_fr(2)          :',F12.2,' - '     ,/, &
!!$                 1X,'zwg_fr(3)          :',F12.2,' - '     ,/, &
!!$                 1X,'u(ke)              :',F12.2,' m/s'    ,/, &
!!$                 1X,'v(ke)              :',F12.2,' m/s'    ,/, &
!!$                 1X,'|u|                :',F12.2,' m/s'    ,/, &
!!$                 1X,'sobs               :',F12.2,' W/m**2' ,/, &
!!$                 1X,'thbs               :',F12.2,' W/m**2' ,/, &
!!$                 1X,'zshfl              :',F12.2,' W/m**2' ,/, &
!!$                 1X,'zlhfl              :',F12.2,' W/m**2' ,/, &
!!$                 1X,'zthfl              :',F12.2,' W/m**2' ,/, &
!!$                 1X,'zg1                :',F12.2,' W/m**2' ,/, &
!!$                 1X,'zeb1               :',F12.2,' W/m**2' ,/, &
!!$                 1X,'tch*|u| original   :',F12.6,' m/s'    ,/, &
!!$                 1X,'tch*|u| after limit:',F12.6,' m/s'    )
        END IF
      END DO
    END DO
  ENDIF

#ifdef NECSX
  CALL collapse(.TRUE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif
! Update indicator for age of snow in top of snow layer

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF(llandmask(i,j)) THEN        ! for land-points only
!subs        IF (w_snow(i,j,nx) <=0.0_ireals) THEN
      IF(llandmask(i,j,ns)) THEN        ! for land-points only
        IF (w_snow(i,j,nx,ns) <=0.0_ireals) THEN
          ! if no snow exists, reinitialize age indicator
!subs          freshsnow(i,j) = 1.0_ireals
          freshsnow(i,j,ns) = 1.0_ireals
        ELSE
          IF ( itype_gscp == 4 ) THEN
            zsnow_rate = prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j) ! [kg/m**2 s]
          ELSE
            zsnow_rate = prs_gsp(i,j)+prs_con(i,j)              ! [kg/m**2 s]
          ENDIF

          ! linear decay equals 1.0 in 28 days
          zdsn_old   = zdt/86400._ireals/28._ireals

          ! linear growth rate equals 1.0 in 1 day for a constant
          ! snow rate of 5 mmH2O (kg/m**2) per day
          zdsn_new   = zdt*zsnow_rate*0.2_ireals   !

          ! reduce decay rate, if new snow is falling and as function of snow
          ! age itself
!subs          zdsn_old   = (zdsn_old - zdsn_new)*freshsnow(i,j)
          zdsn_old   = (zdsn_old - zdsn_new)*freshsnow(i,j,ns)
          zdsn_old   = MAX(zdsn_old,0._ireals)

!subs          freshsnow(i,j) = freshsnow(i,j) + zdsn_new-zdsn_old
          freshsnow(i,j,ns) = freshsnow(i,j,ns) + zdsn_new-zdsn_old

!subs          freshsnow(i,j) = MIN(1._ireals,MAX(0._ireals,freshsnow(i,j)))
          freshsnow(i,j,ns) = MIN(1._ireals,MAX(0._ireals,freshsnow(i,j,ns)))
        END IF
      END IF
    ENDDO
  ENDDO

  !>JH WRITE(*,*) 'TERRA-DIAG END I.1 +++++++++> ',t_so(1,1,:,nnew,:)

!------------------------------------------------------------------------------
! Section I.2: temperatures, water contents (in mH2O), surface pressure, 
!------------------------------------------------------------------------------

  IF(lmulti_snow) THEN
    ! US: this loop does not vectorize until the ksn-loops are put outside
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN     ! for land-points only
        IF (llandmask(i,j,ns)) THEN     ! for land-points only

!subs          IF (w_snow(i,j,nx) > 0.0_ireals) THEN
          IF (w_snow(i,j,nx,ns) > 0.0_ireals) THEN
            ! existence of snow
            ! --> no water in interception store and t_snow < t0_melt
!subs            w_snow(i,j,nx) = w_snow(i,j,nx) + w_i(i,j,nx)
!subs            wtot_snow(i,j,1,nx) = wtot_snow(i,j,1,nx) + w_i(i,j,nx)
!subs            dzh_snow(i,j,1,nx)  = dzh_snow(i,j,1,nx)  + w_i(i,j,nx)/rho_snow_mult(i,j,1,nx)*rho_w
!subs            w_i   (i,j,nx) = 0.0_ireals
            w_snow(i,j,nx,ns) = w_snow(i,j,nx,ns) + w_i(i,j,nx,ns)
            wtot_snow(i,j,1,nx,ns) = wtot_snow(i,j,1,nx,ns) + w_i(i,j,nx,ns)
            dzh_snow(i,j,1,nx,ns)  = dzh_snow(i,j,1,nx,ns)  + w_i(i,j,nx,ns) &
              &                      /rho_snow_mult(i,j,1,nx,ns)*rho_w
            w_i   (i,j,nx,ns) = 0.0_ireals
            DO ksn = 0,ke_snow
!subs              t_snow_mult(i,j,ksn,nx) = MIN (t0_melt - zepsi, t_snow_mult(i,j,ksn,nx) )
              t_snow_mult(i,j,ksn,nx,ns) = MIN (t0_melt - zepsi, t_snow_mult(i,j,ksn,nx,ns) )
            END DO
!subs          ELSE IF (t_snow_mult(i,j,ke_snow,nx) >= t0_melt) THEN
          ELSE IF (t_snow_mult(i,j,ke_snow,nx,ns) >= t0_melt) THEN
            ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
!subs            t_s   (i,j,nx) = MAX (t0_melt + zepsi, t_s(i,j,nx) )
            t_s   (i,j,nx,ns) = MAX (t0_melt + zepsi, t_s(i,j,nx,ns) )
            DO ksn = 0,ke_snow
!subs              t_snow_mult(i,j,ksn,nx) = t_s(i,j,nx)
              t_snow_mult(i,j,ksn,nx,ns) = t_s(i,j,nx,ns)
            END DO
          ELSE
            ! no snow and  t_snow < t0_melt
            ! --> t_snow = t_s and no water w_i in interception store
!subs            t_s   (i,j,nx) = MIN (t0_melt - zepsi, t_s(i,j,nx) )
            t_s   (i,j,nx,ns) = MIN (t0_melt - zepsi, t_s(i,j,nx,ns) )
            DO ksn = 0,ke_snow
!subs              t_snow_mult(i,j,ksn,nx) = t_s(i,j,nx)
              t_snow_mult(i,j,ksn,nx,ns) = t_s(i,j,nx,ns)
            END DO
!subs            w_i   (i,j,nx) = 0.0_ireals
            w_i   (i,j,nx,ns) = 0.0_ireals
          END IF
        END IF
      ENDDO
    ENDDO

  ELSE ! no multi-layer snow

    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN     ! for land-points only
!subs          IF (w_snow(i,j,nx) > 0.0_ireals) THEN
        IF (llandmask(i,j,ns)) THEN     ! for land-points only
          IF (w_snow(i,j,nx,ns) > 0.0_ireals) THEN
            ! existence of snow
            ! --> no water in interception store and t_snow < t0_melt
!subs            w_snow(i,j,nx) = w_snow(i,j,nx) + w_i(i,j,nx)
!subs            w_i   (i,j,nx) = 0.0_ireals
!subs            t_snow(i,j,nx) = MIN (t0_melt - zepsi, t_snow(i,j,nx) )
            w_snow(i,j,nx,ns) = w_snow(i,j,nx,ns) + w_i(i,j,nx,ns)
            w_i   (i,j,nx,ns) = 0.0_ireals
            t_snow(i,j,nx,ns) = MIN (t0_melt - zepsi, t_snow(i,j,nx,ns) )
!subs          ELSE IF (t_snow(i,j,nx) >= t0_melt) THEN
          ELSE IF (t_snow(i,j,nx,ns) >= t0_melt) THEN
            ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
!subs            t_s   (i,j,nx) = MAX (t0_melt + zepsi, t_s(i,j,nx) )
!subs            t_snow(i,j,nx) = t_s(i,j,nx)
            t_s   (i,j,nx,ns) = MAX (t0_melt + zepsi, t_s(i,j,nx,ns) )
            t_snow(i,j,nx,ns) = t_s(i,j,nx,ns)
          ELSE
            ! no snow and  t_snow < t0_melt
            ! --> t_snow = t_s and no water w_i in interception store
!subs            t_s   (i,j,nx) = MIN (t0_melt - zepsi, t_s(i,j,nx) )
!subs            t_snow(i,j,nx) = t_s(i,j,nx)
!subs            w_i   (i,j,nx) = 0.0_ireals
            t_s   (i,j,nx,ns) = MIN (t0_melt - zepsi, t_s(i,j,nx,ns) )
            t_snow(i,j,nx,ns) = t_s(i,j,nx,ns)
            w_i   (i,j,nx,ns) = 0.0_ireals
          END IF
        END IF
      ENDDO
    ENDDO
  ENDIF

  !Inizializations for the next sections
  IF (lmulti_snow) THEN
    ! some preparations for ksn==0 and ksn==1
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN   ! for land-points only
!subs          ztsnow_mult   (i,j,0) = t_snow_mult(i,j,0,nx)
!subs          ztsnow_mult   (i,j,1) = t_snow_mult(i,j,1,nx)
        IF (llandmask(i,j,ns)) THEN   ! for land-points only
          ztsnow_mult   (i,j,0) = t_snow_mult(i,j,0,nx,ns)
          ztsnow_mult   (i,j,1) = t_snow_mult(i,j,1,nx,ns)

!subs          zhh_snow(i,j,1) =  -h_snow(i,j,nnow) + dzh_snow(i,j,1,nx)
!subs          zhm_snow(i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals
          zhh_snow(i,j,1) =  -h_snow(i,j,nnow,ns) + dzh_snow(i,j,1,nx,ns)
          zhm_snow(i,j,1) = (-h_snow(i,j,nnow,ns) + zhh_snow(i,j,1))/2._ireals

!subs          zdzh_snow(i,j,1) = dzh_snow(i,j,1,nx)
!subs          zextinct (i,j,1) = 0.13_ireals*rho_snow_mult(i,j,1,nx)+3.4_ireals
          zdzh_snow(i,j,1) = dzh_snow(i,j,1,nx,ns)
          zextinct (i,j,1) = 0.13_ireals*rho_snow_mult(i,j,1,nx,ns)+3.4_ireals

          ! set ztsnow to ztsnow_mult(ksn=1)
!subs          ztsnow   (i,j) = t_snow_mult(i,j,1,nx)
          ztsnow   (i,j) = t_snow_mult(i,j,1,nx,ns)

          ! initialize zdz_snow (from Section I.3) to zdzh_snow(i,j,1)
!subs          zdz_snow (i,j) = dzh_snow(i,j,1,nx) ! zdzh_snow
!subs          zwsnow   (i,j) = w_snow(i,j,nx)
          zdz_snow (i,j) = dzh_snow(i,j,1,nx,ns) ! zdzh_snow
          zwsnow   (i,j) = w_snow(i,j,nx,ns)
        END IF
      ENDDO
    ENDDO

    DO  ksn = 2,ke_snow
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN   ! for land-points only
!subs            ztsnow_mult(i,j,ksn) = t_snow_mult(i,j,ksn,nx)
          IF (llandmask(i,j,ns)) THEN   ! for land-points only
            ztsnow_mult(i,j,ksn) = t_snow_mult(i,j,ksn,nx,ns)

!subs            zhh_snow   (i,j,ksn) = zhh_snow(i,j,ksn-1) + dzh_snow(i,j,ksn,nx)
            zhh_snow   (i,j,ksn) = zhh_snow(i,j,ksn-1) + dzh_snow(i,j,ksn,nx,ns)
            zhm_snow   (i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals

!subs            zdzh_snow  (i,j,ksn) = dzh_snow(i,j,ksn,nx)
!subs            zextinct   (i,j,ksn) = 0.13_ireals*rho_snow_mult(i,j,ksn,nx)+3.4_ireals
            zdzh_snow  (i,j,ksn) = dzh_snow(i,j,ksn,nx,ns)
            zextinct   (i,j,ksn) = 0.13_ireals*rho_snow_mult(i,j,ksn,nx,ns)+3.4_ireals

            ! build sum over all layers for zdzh_snow in zdz_snow
            zdz_snow   (i,j)     = zdz_snow(i,j) + zdzh_snow(i,j,ksn)
          END IF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    ! set ztsnow to t_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN   ! for land-points only
!subs          ztsnow   (i,j) = t_snow(i,j,nx)
!subs          zwsnow   (i,j) = w_snow(i,j,nx)
!subs          zdz_snow (i,j) = zwsnow(i,j)*rho_w/rho_snow(i,j,nx)
        IF (llandmask(i,j,ns)) THEN   ! for land-points only
          ztsnow   (i,j) = t_snow(i,j,nx,ns)
          zwsnow   (i,j) = w_snow(i,j,nx,ns)
          zdz_snow (i,j) = zwsnow(i,j)*rho_w/rho_snow(i,j,nx,ns)
        END IF
      ENDDO
    ENDDO
  ENDIF

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN   ! for land-points only
      IF (llandmask(i,j,ns)) THEN   ! for land-points only
        ! ztsnow is now set according to lmulti_snow (see above);
        ! therefore we need not take care of it in this loop
        ! ztsnow   (i,j) = t_snow(i,j,nx)
!subs        zts      (i,j) = t_s   (i,j,nx)
        zts      (i,j) = t_s   (i,j,nx,ns)
        zts_pm   (i,j) = zsf_heav(zts   (i,j) - t0_melt)
        ztsnow_pm(i,j) = zsf_heav(ztsnow(i,j) - t0_melt)
!subs        zwin     (i,j) = w_i   (i,j,nx)
        zwin     (i,j) = w_i   (i,j,nx,ns)

        ! moisture and potential temperature of lowest atmospheric layer
        zplow          = p0(i,j,ke) + pp(i,j,ke,nx)
        zqvlow         = qv(i,j,ke,nx)
        zth_low (i,j)  =  t(i,j,ke,nx) *( (ps(i,j,nx)/zplow)**rdocp )

        ! density*transfer coefficient*wind velocity
        zrhoch(i,j)    = ztmch(i,j)*(1._ireals/g) + zepsi

IF(i.eq.1 .and. j.eq.1) &
!JH  print *,ns,'zplow',p0(i,j,ke) , pp(i,j,ke,nx),qv(i,j,ke)
        ! saturation specific humidity for t_s and t_snow and first derivative
        z2iw        = zts_pm(i,j)*b2w + (1._ireals - zts_pm(i,j))*b2i
        z4iw        = zts_pm(i,j)*b4w + (1._ireals - zts_pm(i,j))*b4i
        z234iw      = z2iw*(b3 - z4iw)
        zqs         = zsf_qsat( zsf_psat_iw(zts(i,j), z2iw,z4iw), ps(i,j,nx) )
        zdqs        = zqvlow - zqs
        IF (ABS(zdqs).LT.zepsi) zdqs = 0._ireals
        z2iw        = ztsnow_pm(i,j)*b2w + (1._ireals - ztsnow_pm(i,j))*b2i
        z4iw        = ztsnow_pm(i,j)*b4w + (1._ireals - ztsnow_pm(i,j))*b4i
        z234iw      = z2iw*(b3 - z4iw)
        zqsnow      = zsf_qsat(zsf_psat_iw(ztsnow(i,j),z2iw,z4iw), ps(i,j,nx))
        zdqvtsnow(i,j)= zsf_dqvdt_iw(ztsnow(i,j), zqsnow, z4iw,z234iw)
        zdqsnow     = zqvlow - zqsnow
        IF (ABS(zdqsnow).LT.zepsi) zdqsnow = 0._ireals

        ! potential evaporation at T_snow and Ts
        zep_snow(i,j) = (1._ireals-ztsnow_pm(i,j))*                       &
!subs                                          tfv(i,j)*zrhoch(i,j)*zdqsnow
!subs        zep_s   (i,j) =                   tfv(i,j)*zrhoch(i,j)*zdqs
                                          tfv(i,j,ns)*zrhoch(i,j)*zdqsnow
        zep_s   (i,j) =                   tfv(i,j,ns)*zrhoch(i,j)*zdqs
      END IF
    ENDDO
  ENDDO

  !>JH WRITE(*,*) 'TERRA-DIAG END I.2 +++++++++> ',t_so(1,1,:,nnew,:)


!------------------------------------------------------------------------------
! Section I.3: heat conductivity, frozen fraction, snow and
!            water covered fraction, snow density and  height,
!            volumetric heat content
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        ! snow and water covered fraction
        zrss = MAX( 0.01_ireals, MIN(1.0_ireals,zwsnow(i,j)/cf_snow) ) 
        zrww = MAX( 0.01_ireals, 1.0_ireals -                           &
                                 EXP(MAX( -5.0_ireals, - zwin(i,j)/cf_w) ) )
        zf_snow(i,j) = zrss*zsf_heav(zwsnow(i,j) - zepsi)
        zf_wi  (i,j) = zrww*zsf_heav(zwin  (i,j) - zepsi)

! BR 7/2005 prognostic  snow density
! US DMironov: for the FLake Model, h_snow has to be a prognostic variable but
!              here it is used only diagnostically. So the values are put
!              to every time level
        ! zdz_snow has been computed in the initializations above depending
        ! on lmulti_snow
        ! zdz_snow(i,j)=zwsnow(i,j)*rho_w/rho_snow(i,j,nx)
!subs        h_snow(i,j,nnew) = zdz_snow(i,j)
!subs        h_snow(i,j,nnow) = zdz_snow(i,j)
!subs        IF (.NOT. l2tls) h_snow(i,j,nold) = zdz_snow(i,j)
        h_snow(i,j,nnew,ns) = zdz_snow(i,j)
        h_snow(i,j,nnow,ns) = zdz_snow(i,j)
!        IF (.NOT. l2tls) h_snow(i,j,nold,ns) = zdz_snow(i,j)

!       constrain snow depth and consider this constraint for the computation
!       of average snow density of snow layer
        zdz_snow (i,j) =  zdz_snow (i,j)/MAX(0.01_ireals,zf_snow(i,j))
        zdz_snow (i,j) =  MAX(cdsmin,zdz_snow(i,j))

!       limitation of snow depth to 1.5m for snow cover heat transfer
        zdz_snow_fl(i,j) = MIN(1.5_iREALs, zdz_snow(i,j))
        IF(.NOT. lmulti_snow) &
!subs          zrocs(i,j) = chc_i*zdz_snow_fl(i,j)*rho_snow(i,j,nx)
          zrocs(i,j) = chc_i*zdz_snow_fl(i,j)*rho_snow(i,j,nx,ns)
!BR 7/2005 End
      END IF
    ENDDO
  ENDDO


  !>JH WRITE(*,*) 'TERRA-DIAG END I.3 +++++++++> ',t_so(1,1,:,nnew,:)


!------------------------------------------------------------------------------
! Section I.4: Hydrology, 1.Section
!------------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Section I.4.1: Evaporation from interception store and from snow cover,
  
  !----------------------------------------------------------------------------
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN             ! land points only
      IF (llandmask(i,j,ns)) THEN             ! land points only
        ! Evaporation from interception store if it contains water (wi>0) and
        ! if zep_s<0 indicates potential evaporation for temperature Ts
        ! amount of water evaporated is limited to total content of store
        zdwidt(i,j) = zsf_heav(-zep_s(i,j))      &
                      * MAX(-zrhwddt*zwin(i,j), zf_wi(i,j)*zep_s(i,j))
        ! Evaporation of snow, if snow exists (wsnow>0) and if zep_snow<0
        ! indicates potential evaporation for temperature t_snow
        zdwsndt(i,j) = zsf_heav(-zep_snow(i,j))  &
                       * MAX(-zrhwddt*zwsnow(i,j), zf_snow(i,j)*zep_snow(i,j))
        ! Formation of dew or rime, if zep_s > 0 . distinction between
        ! dew or rime is only controlled by sign of surface temperature
        ! and not effected by presence of snow !
        zrr(i,j)=zsf_heav(zep_s   (i,j))*zep_s   (i,j)*            zts_pm(i,j)
        zrs(i,j)=zsf_heav(zep_snow(i,j))*zep_snow(i,j)*(1.0_ireals-zts_pm(i,j))
      END IF
    ENDDO
  ENDDO
  

  !----------------------------------------------------------------------------
  ! Section I.4.2a: Bare soil evaporation, bucket version
  !----------------------------------------------------------------------------
  
  IF (itype_evsl.EQ.1) THEN   ! Bucket version
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN             ! land-points only
        IF (llandmask(i,j,ns)) THEN             ! land-points only
          IF (zep_s(i,j) < 0.0_ireals) THEN  ! upwards directed potential
                                             ! evaporation
            ! reduction factor for evaporation based on water content of
            ! first soil layer, air dryness point and field capacity
            zbeta = MAX( 0.0_ireals, MIN( 1.0_ireals,                        &
                    (zw_fr(i,j,1)-zadp(i,j)) / (zfcap(i,j)-zadp(i,j))))
            zbeta = zbeta**2

            ! if first soil layer is frozen, allow evaporation at potential
            ! rate; if soil type is rock, evaporation is not allowed at all
            zice        = zsf_heav(1.5_ireals - m_styp(i,j)) ! 1 only for ice
            zevap       = zrock(i,j) + zice          ! 1 for all, but rock
            zbeta  = zbeta + (1._ireals - zbeta)*zice
            zesoil(i,j) = zevap*zbeta*zep_s(i,j)      & ! evaporation
                          *(1._ireals - zf_wi  (i,j)) & ! not water covered
                          *(1._ireals - zf_snow(i,j)) & ! not snow covered
!subs                          * eai(i,j)/sai(i,j)        ! relative source surface
                          * eai(i,j,ns)/sai(i,j,ns)        ! relative source surface
                                                     !  of the bare soil
!            lhfl_bs(i,j) = lh_v * zesoil(i,j)
          END IF ! upwards directed potential evaporation
        END IF   ! land points
      END DO
    END DO
  END IF         ! Bucket version
  
  !----------------------------------------------------------------------------
  ! Section I.4.2b: Bare soil evaporation, BATS version
  !----------------------------------------------------------------------------
  IF (itype_evsl.EQ.2) THEN
    ! Calculation of bare soil evaporation after Dickinson (1984)
    ! Determination of mean water content relative to volume of voids
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN       ! land points only
        IF (llandmask(i,j,ns)) THEN       ! land points only
          IF (zep_s(i,j) < 0.0_ireals) THEN   ! upwards directed potential
                                              ! evaporation
            zsnull(i,j) = zsnull(i,j)/(znull*zporv(i,j))
            ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
            zice   = zsf_heav(1.5_ireals - m_styp(i,j)) ! 1 only for ice
            zevap  = zrock(i,j) + zice                  ! 1 for all soil types
                                                        ! but rock and ice (=0)
            zbeta  = 0.0_ireals
            IF (m_styp(i,j).ge.3) THEN ! Computations not for ice and rocks
              ! auxiliary quantities
              zbf1   = 5.5_ireals - 0.8_ireals* zbedi(i,j)*                &
                      (1.0_ireals + 0.1_ireals*(zbedi(i,j) - 4.0_ireals)*  &
                       clgk0(m_styp(i,j)) )
              zbf2   = (zbedi(i,j) - 3.7_ireals + 5.0_ireals/zbedi(i,j))/  &
                      (5.0_ireals + zbedi(i,j))
              zdmax  = zbedi(i,j)*cfinull*zk0di(i,j)/crhowm 
              zs1(i,j)  = zs1(i,j)/(z1*zporv(i,j))
              zd     = 1.02_ireals*zdmax*zs1(i,j)**(zbedi(i,j) + 2) *      &
                                                (zsnull(i,j)/zs1(i,j))**zbf1
              zck    = (1.0_ireals + 1550.0_ireals*cdmin/zdmax)*zbf2
              ! maximum sustainable moisture flux in the uppermost surface
              ! layer in kg/(s*m**2)
              zfqmax = - rho_w*zck*zd*zsnull(i,j)/SQRT(znull*z1)
              zevapor= MAX(zep_s(i,j),zfqmax)
              IF(zw_fr(i,j,1)+zevapor*(1.0_ireals - zf_wi(i,j))   &
                                     *(1.0_ireals - zf_snow(i,j)) &
!subs                                     *eai(i,j)/sai(i,j)* zdtdrhw/zdzhs(1) &
                                     *eai(i,j,ns)/sai(i,j,ns)* zdtdrhw/zdzhs(1) &
                                     .LE.zadp(i,j)) zevapor = 0._ireals
              zbeta  = zevapor/MIN(zep_s(i,j),-zepsi)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_ireals - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover
            zesoil(i,j) = zevap*zbeta*zep_s(i,j)       & ! evaporation
                          *(1.0_ireals - zf_wi  (i,j)) & ! not water covered
                          *(1.0_ireals - zf_snow(i,j)) & ! not snow covered
!subs                          * eai(i,j)/sai(i,j) ! relative source surface
                          * eai(i,j,ns)/sai(i,j,ns) ! relative source surface
                                              ! of the bare soil
!            lhfl_bs(i,j) = lh_v * zesoil(i,j)
          END IF  ! upwards directed potential evaporation
        END IF    ! land points
      END DO
    END DO
  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section I.4.3a: transpiration by plants, bucket version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.1) THEN   ! Bucket version
  ! for lower layers, even for water saturated soil and and maximum
  ! root extension, transpiration is limited to the potential evaporation rate
  ! the consideration of a root depth (zroot) allows the maximum
  ! transpiration, if the active soil layer is saturated and completely
  ! penetrated by roots

    DO kso = 1,ke_soil  ! loop over soil layers
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN     ! land-points only
          IF (llandmask(i,j,ns)) THEN     ! land-points only
            IF (zep_s(i,j) < 0.0_ireals  &    ! potential evaporation and
               .AND. m_styp(i,j) >= 3) THEN   ! neither ice nor rock
              ! turgor-loss-point minus plant wilting point
              ztlpmwp(i,j) = (zfcap(i,j) - zpwp(i,j))*(0.81_ireals +       &
                    0.121_ireals*ATAN(-86400._ireals*zep_s(i,j) - 4.75_ireals))
              ! determine whether neither surface nor lower boundary of
              ! second layer are below freezing point
!subs              ztgt0 = zsf_heav(t_so(i,j,kso,nx) - t0_melt)
              ztgt0 = zsf_heav(t_so(i,j,kso,nx,ns) - t0_melt)
              ! determine reduction factor beta for this layer
              zbeta   = ztgt0 * MAX(0.0_ireals, MIN(1.0_ireals,            &
!subs                        (zw_fr(i,j,kso) - w_so_ice(i,j,kso,nx)/zdzhs(kso)  &
                        (zw_fr(i,j,kso) - w_so_ice(i,j,kso,nx,ns)/zdzhs(kso)  &
                        -  zpwp(i,j)) / ztlpmwp(i,j)))
              zbeta   = zbeta**2

              ! consider the effect of root depth
              zroot = MIN ( zdzhs(kso), MAX(0.0_ireals,                    &
                         zbwt(i,j) - (zmls(kso) - 0.5_ireals*zdzhs(kso))))
              zr_root = zroot/MAX(zepsi,zbwt(i,j))
              ztrang(i,j,kso) =  zr_root*zbeta             & ! reduction
                                   * zep_s(i,j)            & ! transpiration
                                   * (1._ireals - zf_wi(i,j))     & ! non-water
                                   * (1._ireals - zf_snow(i,j))   & ! non-snow
!subs                                   * tai(i,j)/sai(i,j)       ! transp. surface
                                   * tai(i,j,ns)/sai(i,j,ns)       ! transp. surface
                                     ! relative source surface of the plants
!              lhfl_pl(i,j,kso)= lh_v * ztrang(i,j,kso)
              ztrangs(i,j)    = ztrangs(i,j) + ztrang(i,j,kso)
            END IF  ! upwards directed potential evaporation .AND. m_styp > 2
          END IF    ! land-points only
        END DO
      END DO
    END DO          ! loop over soil layers
  END IF            ! Bucket version


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
      DO   j = jstarts, jends
        DO i = istarts, iends
          zrootdz_int (i,j)= 0.0_ireals   !initialize the root density profile integral
          zwrootdz_int(i,j)= 0.0_ireals   !initialize the root water content integral
        END DO
      END DO
      DO kso   = 1,ke_soil
        DO   j = jstarts, jends
          DO i = istarts, iends
!subs            IF (llandmask(i,j)) THEN ! land points only,
            IF (llandmask(i,j,ns)) THEN ! land points only,
              IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i,j) < 0.0_ireals) THEN  ! upwards directed potential evaporation
                  ! consider the effect of root depth & root density
                  zrootfr = EXP (-zroota(i,j)*zmls(kso)) ! root density
                  zrootdz = zrootfr*MIN(zdzhs(kso),MAX(0.0_ireals, zbwt(i,j)-(zmls(kso) &
                    &       -0.5_ireals*zdzhs(kso)) ) )
                  zrootdz_int(i,j)=zrootdz_int(i,j) + zrootdz
!subs                  zwrootdz(i,j,kso)=zrootdz*(zw_fr(i,j,kso)-w_so_ice(i,j,kso,nx)/zdzhs(kso))
                  zwrootdz(i,j,kso)=zrootdz*(zw_fr(i,j,kso)-w_so_ice(i,j,kso,nx,ns)/zdzhs(kso))
                  zwrootdz_int(i,j)=zwrootdz_int(i,j) + zwrootdz(i,j,kso)
                END IF  ! negative potential evaporation only
              END IF  ! neither ice or rocks
            END IF    ! land-points only
          END DO
        END DO
      END DO

       ! Compute root zone integrated average of liquid water content
       DO   j = jstarts, jends
          DO i = istarts, iends
             zwrootdz_int(i,j)=zwrootdz_int(i,j)/MAX(zrootdz_int(i,j),zepsi)
          END DO
       END DO

    ELSE

      DO kso   = 1,ke_soil
        DO   j = jstarts, jends
          DO i = istarts, iends
!subs            IF (llandmask(i,j)) THEN ! land points only,
            IF (llandmask(i,j,ns)) THEN ! land points only,
              IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i,j) < 0.0_ireals) THEN  ! upwards directed potential
                                                   ! evaporation
                  zropart  = MIN ( zdzhs(kso), MAX(0.0_ireals,                     &
                             zbwt(i,j) - (zmls(kso) - 0.5_ireals*zdzhs(kso))))
!subs                  zropartw(i,j,kso) = zropart*(zw_fr(i,j,kso)-w_so_ice(i,j,kso,nx)/zdzhs(kso))
                  zropartw(i,j,kso) = zropart*(zw_fr(i,j,kso)-w_so_ice(i,j,kso,nx,ns)/zdzhs(kso))
                  zwroot(i,j) = zwroot(i,j) + zropartw(i,j,kso)/MAX(zepsi,zbwt(i,j))
                END IF  ! negative potential evaporation only
              END IF  ! neither ice or rocks
            END IF    ! land-points only
          END DO
        END DO
      END DO

    ENDIF

    ! Determination of the transfer functions CA, CF, and CV

    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN ! land points only,
        IF (llandmask(i,j,ns)) THEN ! land points only,
          IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
            IF (zep_s(i,j) < 0.0_ireals) THEN  ! upwards directed potential evaporation
              im1        = MAX( 1, i-1)
!subs              zuv        = SQRT (u_10m(i,j) **2 + v_10m(i,j)**2 )
!subs              zcatm      = tch(i,j)*zuv           ! Function CA
              zuv        = SQRT (u_10m(i,j,ns) **2 + v_10m(i,j,ns)**2 )
              zcatm      = tch(i,j,ns)*zuv           ! Function CA

              IF(itype_tran == 1) THEN
!subs                zustar     = zuv*SQRT(tcm(i,j))
                zustar     = zuv*SQRT(tcm(i,j,ns))
                zrla       = 1.0_ireals/MAX(cdash*SQRT(zustar),zepsi)
              ELSE
                zrla       = 0._ireals
              ENDIF

           
              ! to compute CV, first the stomatal resistance has to be determined
              ! this requires the determination of the F-functions:
              ! Radiation function
!subs              zpar       = pabs(i,j)  !  PAR
              zpar       = pabs(i,j,ns)  !  PAR
              zf_rad(i,j)= MAX(0.0_ireals,MIN(1.0_ireals,zpar/cparcrit))
              ztlpmwp(i,j) = (zfcap(i,j) - zpwp(i,j))*(0.81_ireals +       &
                      0.121_ireals*ATAN(-86400._ireals*zep_s(i,j) - 4.75_ireals))

              ! Soil water function
              IF (itype_root == 2) THEN
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwrootdz_int(i,j) -  &
                                                zpwp(i,j))/ztlpmwp(i,j)))
              ELSE
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwroot(i,j) -  &
                                                zpwp(i,j))/ztlpmwp(i,j)))
              ENDIF

              ! Temperature function
              if(ntstep.eq.0.and.itype_tran.ne.2) then
!subs                t_2m(i,j)=t(i,j,ke,nx)
                t_2m(i,j,ns)=t(i,j,ke,nx)
              endif
              zf_tem     = MAX(0.0_ireals,MIN(1.0_ireals,4.0_ireals*     &
!subs                           (t_2m(i,j)-t0_melt)*(ctend-t_2m(i,j))/(ctend-t0_melt)**2))
                           (t_2m(i,j,ns)-t0_melt)*(ctend-t_2m(i,j,ns))/(ctend-t0_melt)**2))

              ! Saturation deficit function (function not used, but computations
              ! necessary for determination of  slope of the saturation curve)
              z2iw       = zts_pm(i,j)*b2w + (1._ireals - zts_pm(i,j))*b2i
              z4iw       = zts_pm(i,j)*b4w + (1._ireals - zts_pm(i,j))*b4i
              zepsat     = zsf_psat_iw(t(i,j,ke,nx),z2iw,z4iw)
              zepke      = qv(i,j,ke,nx)*ps(i,j,nx)/                   &
                                        (rdv + o_m_rdv*qv(i,j,ke,nx))
              zf_sat     = MAX(0.0_ireals,MIN(1.0_ireals,1.0_ireals -  &
                                              (zepsat - zepke)/csatdef))
              ! zf_sat paralysed:
              zf_sat     = 1.0_ireals

              IF (lstomata) THEN
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/rsmin2d(i,j) -     &
                             1.0_ireals/crsmax)*zf_rad(i,j)*zf_wat*zf_tem*zf_sat
              ELSE
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/crsmin -     &
                             1.0_ireals/crsmax)*zf_rad(i,j)*zf_wat*zf_tem*zf_sat
              END IF

              zrstom     = 1.0_ireals/zedrstom              ! stomatal resistance
!              rstom(i,j) = zrstom
              zrveg      = zrla + zrstom
              ! Transpiration rate of dry leaves:
!subs              ztraleav(i,j)=zep_s(i,j)*tai(i,j)/(sai(i,j)+zrveg*zcatm)
              ztraleav(i,j)=zep_s(i,j)*tai(i,j,ns)/(sai(i,j,ns)+zrveg*zcatm)
            END IF  ! upwards directed potential evaporation only
          END IF    ! m_styp > 2
        END IF      ! land points
      END DO
    END DO

    ! Consideration of water and snow coverage, distribution to the different
    ! soil layers
    DO     kso       = 1,ke_soil
      DO   j         = jstarts, jends
        DO i         = istarts, iends
!subs          IF (llandmask(i,j)) THEN ! land points only,
          IF (llandmask(i,j,ns)) THEN ! land points only,
            IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
              IF (zep_s(i,j) < 0.0_ireals) THEN    ! upwards potential evaporation
                ztrabpf  = ztraleav(i,j)*                   & ! plant covered part
                           (1.0_ireals - zf_wi(i,j))*       & ! not water covered
                           (1.0_ireals - zf_snow(i,j))        ! not snow covered

                ! for root distribution
                IF (itype_root == 2) THEN
                  ztrfr    = zwrootdz(i,j,kso)/(zrootdz_int(i,j)*zwrootdz_int(i,j))
                  ztrang(i,j,kso) = ztrabpf*ztrfr
                ELSE
                  zrootfc = zropartw(i,j,kso)/(zwroot(i,j) + zepsi)
                  ztrang(i,j,kso) = ztrabpf*zrootfc/MAX(zepsi,zbwt(i,j))
                  IF(zw_fr(i,j,kso)+ztrang(i,j,kso)*zdtdrhw/zdzhs(kso) &
                                    .LT.zpwp(i,j)) ztrang(i,j,kso) = 0._ireals
                ENDIF
!                lhfl_pl(i,j,kso)= lh_v * ztrang(i,j,kso)
                ztrangs(i,j)    = ztrangs(i,j) + ztrang(i,j,kso)
              END IF  ! upwards directed potential evaporation only
            END IF    ! m_styp > 2
          END IF      ! land points
        END DO
      END DO
    END DO          ! loop over soil layers

  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section I.4.4: total evapotranspiration and 
  !              associated ficticious soil humidity qv_s
  !----------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN   ! land-points only
      IF (llandmask(i,j,ns)) THEN   ! land-points only
        ze_sum = zdwsndt(i,j  )  & ! evaporation of snow
               + zdwidt (i,j  )  & ! evaporation from interception store
               + zesoil (i,j  )  & ! evaporation from bare soil
               + ztrangs(i,j  )  & ! transpiration from all soil layers
               + zrr    (i,j  )  & ! formation of dew
               + zrs    (i,j  )    ! formation of rime
!subs        qv_s(i,j,nx  ) = qv (i,j,ke,nx) - ze_sum /(zrhoch(i,j) + zepsi)
!subs        qv_s(i,j,nnew) = qv_s(i,j,nx)
        qv_s(i,j,nx  ,ns) = qv (i,j,ke,nx) - ze_sum /(zrhoch(i,j) + zepsi)
        qv_s(i,j,nnew,ns) = qv_s(i,j,nx,ns)

 
     END IF     ! land points
    END DO
  END DO

!------------------------------------------------------------------------------
! End of former module procedure terra1_multlay
!------------------------------------------------------------------------------

  !>JH WRITE(*,*) 'TERRA-DIAG END I.4 +++++++++> ',t_so(1,1,:,nnew,:)

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
  do kso=1,ke_soil+1
  if (zzhls(kso)<=2.5_ireals) i250=kso
  ke_soil_hy = i250
  end do

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
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN                 ! land-points
          IF (llandmask(i,j,ns)) THEN                 ! land-points
            zdtsnowdt_mult(i,j,ksn)  = 0.0_ireals
            zdtsdt(i,j)     = 0.0_ireals
          END IF
        END DO
      END DO
    END DO
  ELSE
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN                 ! land-points
        IF (llandmask(i,j,ns)) THEN                 ! land-points
          zdtsnowdt(i,j)  = 0.0_ireals
          zdtsdt(i,j)     = 0.0_ireals
        END IF
      END DO
    END DO
  END IF

  !>JH WRITE(*,*) 'TERRA-DIAG END II.1 +++++++++> ',t_so(1,1,:,nnew,:)

!------------------------------------------------------------------------------
! Section II.2: Prepare basic surface properties and create some local
!               arrays of surface related quantities (for land-points only)
!               Initialise some fields
!------------------------------------------------------------------------------

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN                 ! land-points
!subs      mstyp        = NINT(soiltyp_subs(i,j))
    IF (llandmask(i,j,ns)) THEN                 ! land-points
      mstyp        = soiltyp_subs(i,j,ns)
      m_styp(i,j)  = mstyp
      zsandf(i,j)  = csandf(mstyp)
      zclayf(i,j)  = cclayf(mstyp)
      zgsb   (i,j) = 0.0_ireals
      ztrangs(i,j) = 0.0_ireals
    END IF
  END DO
END DO


DO   kso = 1,ke_soil+1
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN     ! land-points only
!subs        ziw_fr(i,j,kso) = w_so_ice(i,j,kso,nx)/zdzhs(kso)   ! ice frac.
      IF (llandmask(i,j,ns)) THEN     ! land-points only
        ziw_fr(i,j,kso) = w_so_ice(i,j,kso,nx,ns)/zdzhs(kso)   ! ice frac.
        zlw_fr(i,j,kso) = zw_fr(i,j,kso) - ziw_fr(i,j,kso)  ! liquid water frac.
        zroc(i,j,kso)   = zrocg(i,j) + rho_w*zlw_fr(i,j,kso)*chc_w +          &
                                       rho_w*ziw_fr(i,j,kso)*chc_i
        IF (kso<=ke_soil) THEN
          ztrangs(i,j) = ztrangs(i,j) + ztrang(i,j,kso)
          zdwgdt(i,j,kso) = 0.0_ireals
        END IF
        zflmg(i,j,kso)  = 0.0_ireals
        zrunoff_grav(i,j,kso)  = 0.0_ireals
      END IF
    END DO
  END DO
END DO      !soil layers

  !>JH WRITE(*,*) 'TERRA-DIAG END II.2 +++++++++> ',t_so(1,1,:,nnew,:)

!------------------------------------------------------------------------------
! Section II.3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF(llandmask(i,j))THEN     ! land-points only
      IF(llandmask(i,j,ns))THEN     ! land-points only

        ! store forcing terms due to evapotranspiration, formation of dew
        ! and rime for later use
        zverbo(i,j) = zdwidt(i,j) + zesoil(i,j) + ztrangs(i,j) +              &
                        (1._ireals-zf_snow(i,j))*(zrr(i,j) + zrs(i,j))

        zversn(i,j) = zdwsndt(i,j) + zrs(i,j)                                 &
                                  + zsf_heav (zwsnow(i,j) - zepsi) * zrr(i,j)

        ! add grid scale and convective precipitation (and graupel, if present)
        ! to dew and rime
        zrr(i,j) = zrr(i,j) + prr_con(i,j) + prr_gsp(i,j)
        IF ( itype_gscp == 4 ) THEN
          zrs(i,j) = zrs(i,j) + prs_con(i,j) + prs_gsp(i,j) + prg_gsp(i,j)
        ELSE
          zrs(i,j) = zrs(i,j) + prs_con(i,j) + prs_gsp(i,j)
        ENDIF

!em subs
        IF(ns.eq.2) THEN
          IF(subsfrac(i,j,ns+1).gt.0._ireals) THEN
            IF((1._ireals-ztsnow_pm(i,j))*zrs(i,j) > 0.0_ireals) zrs(i,j) = 0.0_ireals
            IF((1._ireals-ztsnow_pm(i,j))*zrr(i,j) > 0.0_ireals) zrr(i,j) = 0.0_ireals
          END IF
        ELSE
          IF(subsfrac(i,j,ns).gt.0._ireals) THEN
            IF((1._ireals-ztsnow_pm(i,j))*zrs(i,j) > 0.0_ireals) zrs(i,j) = &
              zrs(i,j)*(subsfrac(i,j,ns-1)+subsfrac(i,j,ns))/MAX(subsfrac(i,j,ns),zepsi)
            IF((1._ireals-ztsnow_pm(i,j))*zrr(i,j) > 0.0_ireals) zrr(i,j) = &
              zrr(i,j)*(subsfrac(i,j,ns-1)+subsfrac(i,j,ns))/MAX(subsfrac(i,j,ns),zepsi) 
          END IF
        END IF 

        ! infiltration and surface run-off

        ! ice free fraction of first soil layer scaled by pore volume
        ! is used as reduction factor for infiltration rate
        zfr_ice_free     = 1._ireals-ziw_fr(i,j,1)/zporv(i,j)

        ! subtract evaporation from interception store to avoid negative
        ! values due to sum of evaporation+infiltration
        zwinstr = zwin(i,j) + zdwidt(i,j)*zdtdrhw
        zwinstr = MAX(0.0_ireals,zwinstr)

        ! maximum infiltration rate of the soil (rock/ice/water-exclusion
        zinfmx = zrock(i,j)*zfr_ice_free*csvoro &
!subs                 *( cik1*MAX(0.5_ireals,plcov(i,j))*MAX(0.0_ireals,           &
                 *( cik1*MAX(0.5_ireals,plcov(i,j,ns))*MAX(0.0_ireals,           &
                 zporv(i,j)-zw_fr(i,j,1))/zporv(i,j) + zik2(i,j) )

        ! to avoid pore volume water excess of the uppermost layer by 
        ! infiltration
        zinfmx = MIN(zinfmx, (zporv(i,j) - zw_fr(i,j,1))*zdzhs(1)*zrhwddt)

        ! to avoid infiltration at snow covered parts of soil surface
        zinfmx = zinfmx*(1._ireals - zf_snow(i,j))

!subs        zwimax = cwimax_ml*(1._ireals + plcov(i,j)*5._ireals)
        zwimax = cwimax_ml*(1._ireals + plcov(i,j,ns)*5._ireals)
        zalf   = SQRT(MAX(0.0_ireals,1.0_ireals - zwinstr/zwimax))

        ! water supply from interception store (if Ts above freezing)
        zinf   = zfr_ice_free*zwinstr*rho_w/ctau_i

        ! possible contribution of rain to infiltration
        IF (zrr(i,j)-zepsi > 0.0_ireals) THEN
          zalf = MAX( zalf,                                                   &
                 (zrhwddt*MAX(0.0_ireals, zwimax-zwinstr) + zinf)/zrr(i,j) )
          zalf = MAX( 0.01_ireals, MIN(1.0_ireals, zalf) )

          ! if rain falls onto snow, all rain is considered for infiltration
          ! as no liquid water store is considered in the snowpack
          IF (zwsnow(i,j) > 0.0_ireals) zalf = 0.0_ireals
        END IF
        ! Increase infiltration and reduce surface runoff (bugfix)
        zalf = 0.0_ireals
        ! rain freezes on the snow surface
        IF (lmulti_snow .AND. zwsnow(i,j) > 0.0_ireals) zalf = 1.0_ireals
        ! add rain contribution to water supply for infiltration
        zvers = zinf + (1._ireals - zalf)*zrr(i,j)
        ! final infiltration rate limited by maximum value
        zinfil(i,j) = MIN(zinfmx,zvers)

        ! surface run-off (residual of potential minus actual infiltration)
        zro_inf       = zvers - zinfil(i,j)
!subs        runoff_s(i,j) = runoff_s(i,j) + zro_inf*zroffdt
        runoff_s(i,j,ns) = runoff_s(i,j,ns) + zro_inf*zroffdt

        ! change of snow water and interception water store
        ! (negligible residuals are added to the run-off)

        ! snow store
        zdwsndtt = zrs(i,j) + zdwsndt(i,j)
        zwsnstr  = zwsnow(i,j) + zdwsndtt*zdtdrhw
        zwsnstr  = MAX(0.0_ireals, zwsnstr) ! avoid negative values (security)
        zdwseps  = 0.0_ireals
        IF (zwsnstr > 0.0_ireals .AND. zwsnstr < zepsi) THEN
!         IF (ztsnow_pm(i,j) > 0.0_ireals) THEN
            zdwseps    = zwsnstr*zrhwddt
!subs            runoff_s(i,j) = runoff_s(i,j) + zdwseps*zroffdt
            runoff_s(i,j,ns) = runoff_s(i,j,ns) + zdwseps*zroffdt
            zdwsndtt   = zdwsndtt - zdwseps
!         END IF
        END IF
        zdwsndt(i,j) = zdwsndtt

        ! interception store
        zdwidtt  = zalf*zrr(i,j) + zdwidt(i,j)-zinf 
        zwinstr  = zwin(i,j) + zdwidtt*zdtdrhw
        zwinstr  = MAX(0.0_ireals, zwinstr) !avoid negative values (security)
        zdwieps  = 0.0_ireals
        IF (zwinstr > 0.0_ireals .AND. zwinstr < zepsi) THEN
          zdwieps    = zwinstr*zrhwddt
!subs          runoff_s(i,j)= runoff_s(i,j) + zdwieps*zroffdt
          runoff_s(i,j,ns)= runoff_s(i,j,ns) + zdwieps*zroffdt
          zdwidtt    = zdwidtt - zdwieps
          zwinstr    = 0.0_ireals
        END IF
        ! add excess over zwimax to runoff
        zro_wi       = zrhwddt*MAX( 0.0_ireals, zwinstr-zwimax )
        zdwidtt      = zdwidtt - zro_wi
        zdwidt(i,j)  = zdwidtt
!subs        runoff_s(i,j)= runoff_s(i,j) + zro_wi*zroffdt
        runoff_s(i,j,ns)= runoff_s(i,j,ns) + zro_wi*zroffdt
      END IF            ! land-points only
    END DO
  END DO

  !>JH WRITE(*,*) 'TERRA-DIAG END II.3 +++++++++> ',t_so(1,1,:,nnew,:)

!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

! uppermost layer, kso = 1
DO   j = jstarts, jends
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
!subs    IF (llandmask(i,j)) THEN      ! land-points only
    IF (llandmask(i,j,ns)) THEN      ! land-points only
      IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
        zice_fr_kso   = ziw_fr(i,j,1)
        zice_fr_ksop1 = ziw_fr(i,j,2)
        zlw_fr_kso  = zlw_fr(i,j,1)
        zlw_fr_ksop1= zlw_fr(i,j,2)

        ! compute reduction factor for transport coefficients
        zfr_ice          = max (zice_fr_kso,zice_fr_ksop1)
        zredp05          = 1._ireals-zfr_ice/MAX(zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)

        ! interpolated scaled liquid water fraction at layer interface
        zlw_fr_ksop05 = 0.5_ireals*(zdzhs(2)*zlw_fr_kso+zdzhs(1)*zlw_fr_ksop1) &
                                             /zdzms(2)
        zdlw_fr_ksop05 = zredp05*zdw(i,j)*EXP(zdw1(i,j)*                       &
                         (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
        zklw_fr_ksop05 = zredp05*zkw(i,j)*EXP(zkw1(i,j)*                       &
                         (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )


        ! coefficients for implicit flux computation
        z1dgam1 = zdt/zdzhs(1) 
        zgam2p05 = zdlw_fr_ksop05/zdzms(2)
        zaga(i,j,1) = 0._ireals
        zagb(i,j,1) = 1._ireals+zalfa*zgam2p05*z1dgam1
        zagc(i,j,1) = -zalfa * zgam2p05*z1dgam1
        zagd(i,j,1) = zlw_fr(i,j,1) + zinfil(i,j)*z1dgam1/rho_w  &
                     -zklw_fr_ksop05*z1dgam1                     &
                     +(1._ireals - zalfa)* zgam2p05*z1dgam1*(zlw_fr_ksop1 - zlw_fr_kso)  &
                          +         zgam2p05*z1dgam1*(zice_fr_ksop1-zice_fr_kso)

        ! explicit part of soil surface water flux:
        zflmg (i,j,1) = - zinfil(i,j)! boundary value for soil water transport
      ENDIF
    ENDIF
  END DO
END DO

! inner layers 2 <=kso<=ke_soil_hy-1
DO   kso =2,ke_soil_hy-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      ! sedimentation and capillary transport in soil
!subs      IF (llandmask(i,j)) THEN      ! land-points only
      IF (llandmask(i,j,ns)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zice_fr_ksom1 = ziw_fr(i,j,kso-1)
          zice_fr_kso   = ziw_fr(i,j,kso  )
          zice_fr_ksop1 = ziw_fr(i,j,kso+1)
          zlw_fr_ksom1  = zlw_fr(i,j,kso-1)
          zlw_fr_kso    = zlw_fr(i,j,kso  )
          zlw_fr_ksop1  = zlw_fr(i,j,kso+1)
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
          zdlw_fr_ksom05= zredm05*zdw(i,j)*EXP( zdw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zdlw_fr_ksop05= zredp05*zdw(i,j)*EXP( zdw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
          zklw_fr_ksom05= zredm05*zkw(i,j)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zklw_fr_ksop05= zredp05*zkw(i,j)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )

          ! coefficients for implicit flux computation
          z1dgam1 = zdt/zdzhs(kso)
          zgam2m05  = zdlw_fr_ksom05/zdzms(kso)
          zgam2p05  = zdlw_fr_ksop05/zdzms(kso+1)
          zaga (i,j,kso) = -zalfa*zgam2m05*z1dgam1
          zagc (i,j,kso) = -zalfa*zgam2p05*z1dgam1
          zagb (i,j,kso) = 1._ireals +zalfa*(zgam2m05+zgam2p05)*z1dgam1
          zagd (i,j,kso) = zlw_fr(i,j,kso)+                               &
                                z1dgam1*(-zklw_fr_ksop05+zklw_fr_ksom05)+ &
                                (1._ireals-zalfa)*z1dgam1*                &
                                (zgam2p05*(zlw_fr_ksop1-zlw_fr_kso  )     &
                                -zgam2m05*(zlw_fr_kso  -zlw_fr_ksom1)   ) &
                               +z1dgam1*                                  &
                                (zgam2p05*(zice_fr_ksop1-zice_fr_kso  )   &
                                -zgam2m05*(zice_fr_kso-zice_fr_ksom1)   )

          !soil water flux, explicit part, for soil water flux investigations
          ! only)
          zflmg(i,j,kso) = rho_w &
            &              *(zdlw_fr_ksom05*(zlw_fr_kso+zice_fr_kso-zlw_fr_ksom1-zice_fr_ksom1) &
            &              / zdzms(kso) - zklw_fr_ksom05)

          IF(kso==ke_soil_hy-1) THEN
            zflmg(i,j,kso+1)=rho_w &
              &              *(zdlw_fr_ksop05*(zlw_fr_ksop1+zice_fr_ksop1-zlw_fr_kso-zice_fr_kso) &
              &              / zdzms(kso+1) - zklw_fr_ksop05)
          ENDIF
        ENDIF
      ENDIF
    END DO
  END DO
END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN      ! land-points only
      IF (llandmask(i,j,ns)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          ! lowest active hydrological layer ke_soil_hy-1
          zice_fr_ksom1 = ziw_fr(i,j,ke_soil_hy-1)
          zice_fr_kso   = ziw_fr(i,j,ke_soil_hy  )
          zlw_fr_ksom1  = zlw_fr(i,j,ke_soil_hy-1)
          zlw_fr_kso    = zlw_fr(i,j,ke_soil_hy  )
          zlw_fr_ksom05 = 0.5_ireals*(zdzhs(ke_soil_hy-1)*zlw_fr_kso+ &
                              zdzhs(ke_soil_hy)*zlw_fr_ksom1)/zdzms(ke_soil_hy)

          zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
          zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
          zdlw_fr_ksom05= zredm05*zdw(i,j)*EXP( zdw1(i,j)* &
                            (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

          z1dgam1 = zdt/zdzhs(ke_soil_hy)
          zgam2m05  = zdlw_fr_ksom05/zdzms(ke_soil_hy)
          zklw_fr_ksom05= zredm05*zkw(i,j)*EXP( zkw1(i,j)* &
                            (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zaga(i,j,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
          zagb(i,j,ke_soil_hy) = 1._ireals+ zalfa*zgam2m05*z1dgam1
          zagc(i,j,ke_soil_hy) = 0.0_ireals
          zagd(i,j,ke_soil_hy) = zlw_fr(i,j,ke_soil_hy)+z1dgam1*zklw_fr_ksom05 &
                            +(1._ireals-zalfa)*z1dgam1*                              &
                             zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)            &
                            +z1dgam1*                                         &
                             zgam2m05*(zice_fr_ksom1-zice_fr_kso )   
        ENDIF
      ENDIF
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        ! generalized upper boundary condition
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zagc(i,j,1) = zagc(i,j,1)/zagb(i,j,1)
          zagd(i,j,1) = zagd(i,j,1)/zagb(i,j,1)
        ENDIF
      END IF          ! land-points only
    END DO
  END DO

DO kso=2,ke_soil_hy-1
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zzz = 1._ireals/(zagb(i,j,kso) - zaga(i,j,kso)*zagc(i,j,kso-1))
          zagc(i,j,kso) = zagc(i,j,kso) * zzz
          zagd(i,j,kso) = (zagd(i,j,kso) - zaga(i,j,kso)*zagd(i,j,kso-1)) * zzz
        ENDIF
      END IF          ! land-points only
    END DO
  END DO
END DO                ! soil layers

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
           zage(i,j,ke_soil_hy) = (zagd(i,j,ke_soil_hy)-zaga(i,j,ke_soil_hy)*  &
                             zagd(i,j,ke_soil_hy-1))/                          &
                            (zagb(i,j,ke_soil_hy) - zaga(i,j,ke_soil_hy)*      &
                             zagc(i,j,ke_soil_hy-1))
        ENDIF
      END IF          ! land-points only
   END DO
 END DO

 DO kso = ke_soil_hy-1,1,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zage(i,j,kso)     = zagd(i,j,kso) - zagc(i,j,kso)*zage(i,j,kso+1)
          ! compute implicit part of new liquid water content and add existing
          ! ice content
!subs          w_so(i,j,kso,nnew) = zage(i,j,kso)*zdzhs(kso) + w_so_ice(i,j,kso,nx)
          w_so(i,j,kso,nnew,ns) = zage(i,j,kso)*zdzhs(kso) + w_so_ice(i,j,kso,nx,ns)
        END IF
      END IF          ! land-points only
    END DO
  END DO
END DO                ! soil layers

!lowest active hydrological level
 DO   j = jstarts, jends
   DO i = istarts, iends
!subs     IF (llandmask(i,j)) THEN          ! land-points only
     IF (llandmask(i,j,ns)) THEN          ! land-points only
       IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
         ! boundary values ensure that the calculation below leaves the climate
         ! layer water contents unchanged compute implicit part of new liquid
         ! water content and add existing ice content
!subs         w_so(i,j,ke_soil_hy,nnew) = zage(i,j,ke_soil_hy)*zdzhs(ke_soil_hy) + &
!subs                         w_so_ice(i,j,ke_soil_hy,nx)
         w_so(i,j,ke_soil_hy,nnew,ns) = zage(i,j,ke_soil_hy)*zdzhs(ke_soil_hy) + &
                         w_so_ice(i,j,ke_soil_hy,nx,ns)
       END IF 
     END IF          ! land-points only
   END DO
 END DO

! to ensure vertical constant water concentration profile beginning at 
! layer ke_soil_hy for energetic treatment only
! soil water climate layer(s)

IF (itype_hydbound == 3) THEN
  ! ground water as lower boundary of soil column
  DO kso = ke_soil_hy+1,ke_soil+1
   DO   j = jstarts, jends
     DO i = istarts, iends
!subs       IF (llandmask(i,j)) THEN          ! land-points only
       IF (llandmask(i,j,ns)) THEN          ! land-points only
         IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
!subs           w_so(i,j,kso,nnew) = zporv(i,j)*zdzhs(kso)
           w_so(i,j,kso,nnew,ns) = zporv(i,j)*zdzhs(kso)
         END IF
       END IF          ! land-points only
     END DO
   END DO
  END DO
ELSE
  DO kso = ke_soil_hy+1,ke_soil+1
   DO   j = jstarts, jends
     DO i = istarts, iends
!subs       IF (llandmask(i,j)) THEN          ! land-points only
       IF (llandmask(i,j,ns)) THEN          ! land-points only
         IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
!subs           w_so(i,j,kso,nnew) = w_so(i,j,kso-1,nnew)*zdzhs(kso)/zdzhs(kso-1)
           w_so(i,j,kso,nnew,ns) = w_so(i,j,kso-1,nnew,ns)*zdzhs(kso)/zdzhs(kso-1)
         END IF
       END IF          ! land-points only
     END DO
   END DO
  END DO
ENDIF

! combine implicit part of sedimentation and capillary flux with explicit part
! (for soil water flux investigations only)
DO kso = 2,ke_soil+1
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zice_fr_ksom1 = ziw_fr(i,j,kso-1)
          zice_fr_kso   = ziw_fr(i,j,kso)
!subs          zlw_fr_ksom1_new= w_so(i,j,kso-1,nnew)/zdzhs(kso-1) - zice_fr_ksom1
!subs          zlw_fr_kso_new  = w_so(i,j,kso  ,nnew)/zdzhs(kso  ) - zice_fr_kso
!subs          zlw_fr_ksom1  = w_so(i,j,kso-1,nx)/zdzhs(kso-1) - zice_fr_ksom1
!subs          zlw_fr_kso    = w_so(i,j,kso  ,nx)/zdzhs(kso  ) - zice_fr_kso
          zlw_fr_ksom1_new= w_so(i,j,kso-1,nnew,ns)/zdzhs(kso-1) - zice_fr_ksom1
          zlw_fr_kso_new  = w_so(i,j,kso  ,nnew,ns)/zdzhs(kso  ) - zice_fr_kso
          zlw_fr_ksom1  = w_so(i,j,kso-1,nx,ns)/zdzhs(kso-1) - zice_fr_ksom1
          zlw_fr_kso    = w_so(i,j,kso  ,nx,ns)/zdzhs(kso  ) - zice_fr_kso
          !... additionally for runoff_g at lower level of lowest active water
          ! layer calculated with (upstream) main level soil water content
          ! compute reduction factor for transport coefficients
          zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
          zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)

        ! interpolated liquid water content at interface to layer above
          zlw_fr_ksom05 =0.5_ireals*(zdzhs(kso)*zlw_fr_ksom1+zdzhs(kso-1)*zlw_fr_kso) &
                            /zdzms(kso)
          zdlw_fr_ksom05= zredm05*zdw(i,j)*EXP(zdw1(i,j) *  &
                          (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zklw_fr_ksom05= zredm05*zkw(i,j) * EXP(zkw1(i,j)* &
                          (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

          IF (kso> ke_soil_hy) zdlw_fr_ksom05=0.0_ireals   ! no flux gradient 
                                                    ! contribution below 2.5m
          IF (kso> ke_soil_hy) zklw_fr_ksom05=0.0_ireals   ! no gravitation flux below 2.5m
          zflmg(i,j,kso) =                              &
                   (1._ireals-zalfa) * zflmg(i,j,kso) + & ! explicit flux component
                              zalfa  * rho_w *          & ! implicit flux component
                   (zdlw_fr_ksom05 * (zlw_fr_kso_new+zice_fr_kso-zlw_fr_ksom1_new-zice_fr_ksom1) &
                   /zdzms(kso) - zklw_fr_ksom05)
          zredm = 1._ireals-zice_fr_kso/(zlw_fr_kso+zice_fr_kso)
          zklw_fr_kso_new = zredm*zkw(i,j) * EXP(zkw1(i,j)* &
                            (zporv(i,j) - zlw_fr_kso_new)/(zporv(i,j) - zadp(i,j)) )

          ! actual gravitation water flux
!subs          IF (w_so(i,j,kso,nnew).LT.1.01_ireals*zadp(i,j)*zdzhs(kso)) zklw_fr_kso_new=0._ireals
          IF (w_so(i,j,kso,nnew,ns).LT.1.01_ireals*zadp(i,j)*zdzhs(kso)) zklw_fr_kso_new=0._ireals
          zrunoff_grav(i,j,kso) =  - rho_w * zklw_fr_kso_new

          ! ground water as lower boundary of soil column
          IF ((kso == ke_soil_hy+1).and.(itype_hydbound == 3)) THEN
             zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )

             zdlw_fr_kso = zredm05*zdw(i,j)*EXP(zdw1(i,j) *  &
                  (zporv(i,j)-zlw_fr_kso_new)/(zporv(i,j)-zadp(i,j)) )
             zklw_fr_kso = zredm05*zkw(i,j) * EXP(zkw1(i,j)* &
                  (zporv(i,j)-zlw_fr_kso_new)/(zporv(i,j)-zadp(i,j)) )
             zklw_fr_ksom1 = zredm05*zkw(i,j) * EXP(zkw1(i,j)* &
                  (zporv(i,j)-zlw_fr_ksom1_new)/(zporv(i,j)-zadp(i,j)) )

             zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
             zrunoff_grav(i,j,ke_soil_hy)=zrunoff_grav(i,j,ke_soil_hy)+ &
                  zdhydcond_dlwfr / &
                  (1.0_ireals-exp(-zdhydcond_dlwfr/zdlw_fr_kso*0.5_ireals*zdzms(ke_soil_hy+1)))* &
                  zdelta_sm
          ENDIF
        END IF 
      END IF          ! land-points only
    END DO
  END DO
END DO


  DO  kso = 1,ke_soil
    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = zsf_heav(0.5_ireals + msr_off - kso)  ! 1.0 for 'surface runoff'
    zro_gfak = 1._ireals - zro_sfak                  ! 1.0 for 'ground runoff'

    IF (kso==i250) THEN  
      zfmb_fak = 1.0_ireals
    ELSE
      zfmb_fak = 0.0_ireals
    END IF

    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN      ! land-points only
        IF (llandmask(i,j,ns)) THEN      ! land-points only
          ! sedimentation and capillary transport in soil
          IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type

            ! first runoff calculation without consideration of
            ! evapotranspiration
            !zdwg   =  zflmg(i,j,kso+1) - zflmg(i,j,kso)
            !zdwg calculated above by flux divergence has to be aequivalent with 
!subs            zdwg =  (w_so(i,j,kso,nnew)/zdzhs(kso)-zw_fr(i,j,kso))*zdzhs(kso) &
            zdwg =  (w_so(i,j,kso,nnew,ns)/zdzhs(kso)-zw_fr(i,j,kso))*zdzhs(kso) &
                                                                   /zdtdrhw
            zdwg =  zdwg + zrunoff_grav(i,j,kso)*zfmb_fak
            zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,j,kso) -     &
                       zfcap(i,j))/MAX(zporv(i,j) - zfcap(i,j),zepsi)) )
            zredfu = zsf_heav(zdwg)*zredfu
            zro    = zdwg*zredfu
            zdwg   = zdwg*(1._ireals - zredfu)

            ! add evaporation (znlgw1f: first layer only)
            ! and transpiration (for each layer)
            zdwg   = zdwg + znlgw1f(kso) * zesoil(i,j) + ztrang (i,j,kso)
            zwgn   = zw_fr(i,j,kso) + zdtdrhw*zdwg/zdzhs(kso)
            zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zwgn - zporv(i,j))
            zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zadp(i,j) - zwgn )
            zdwgdt(i,j,kso)= zdwg + zkorr - zro2
            zro    = zro      + zro2
!subs            runoff_s(i,j) = runoff_s(i,j) + zro*zro_sfak*zroffdt
!subs            runoff_g(i,j) = runoff_g(i,j) + zro*zro_gfak*zroffdt
            runoff_s(i,j,ns) = runoff_s(i,j,ns) + zro*zro_sfak*zroffdt
            runoff_g(i,j,ns) = runoff_g(i,j,ns) + zro*zro_gfak*zroffdt
            ! runoff_g reformulation:
!subs            runoff_g(i,j) = runoff_g(i,j) - (zrunoff_grav(i,j,kso) * zfmb_fak &
            runoff_g(i,j,ns) = runoff_g(i,j,ns) - (zrunoff_grav(i,j,kso) * zfmb_fak &
                                          + zkorr) * zroffdt
          END IF          ! ice/rock-exclusion
        END IF   ! land-points only
      END DO
    END DO
  END DO         ! end loop over soil layers

  !>JH WRITE(*,*) 'TERRA-DIAG END II.4 +++++++++> ',t_so(1,1,:,nnew,:)

!------------------------------------------------------------------------------
! Section II.5: Soil surface heat flux (thermal forcing)
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)
        zgstr =   sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i,j))* &
!subs                  zts(i,j) + zf_snow(i,j)*ztsnow_mult(i,j,1) )**4 + thbs(i,j)
                  zts(i,j) + zf_snow(i,j)*ztsnow_mult(i,j,1) )**4 + thbs(i,j,ns)
        zthsnw(i,j) = - sigma*(1._ireals - Ctalb)*ztsnow_mult(i,j,1)**4 + zgstr
        zthsoi(i,j) = - sigma*(1._ireals - Ctalb)*zts(i,j)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux
!subs        zrnet_s = (1._ireals - zf_snow(i,j))*(sobs(i,j)+zthsoi(i,j))
        zrnet_s = (1._ireals - zf_snow(i,j))*(sobs(i,j,ns)+zthsoi(i,j))
        zshfl_s = (1._ireals - zf_snow(i,j))*cp_d*zrhoch(i,j)*  &
                                                      (zth_low(i,j) - zts(i,j))
        zlhfl_s = (zts_pm(i,j)*lh_v + (1._ireals-zts_pm(i,j))*lh_s)*zverbo(i,j)
        zsprs      = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i,j)*zrs(i,j) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs        = - lh_f*zrs(i,j)
          zdwidt (i,j) = zdwidt (i,j) + zrs(i,j)
          zdwsndt(i,j) = zdwsndt(i,j) - zrs(i,j)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
!subs          zwimax       = cwimax_ml*(1._ireals + plcov(i,j)*5._ireals)
          zwimax       = cwimax_ml*(1._ireals + plcov(i,j,ns)*5._ireals)
          zwinstr      = zwin(i,j) + zdwidt(i,j)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zro        = (zwinstr - zwimax)*zrhwddt
            zdwidt(i,j)= zdwidt(i,j) - zro
!subs            runoff_s(i,j)  = runoff_s(i,j) + zro*zroffdt
            runoff_s(i,j,ns)  = runoff_s(i,j,ns) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF ((1._ireals-ztsnow_pm(i,j))*zrr(i,j) > 0.0_ireals) THEN
          zsprs   = lh_f*zrr(i,j)
          zdwidt (i,j) = zdwidt (i,j) - zrr(i,j)
          zdwsndt(i,j) = zdwsndt(i,j) + zrr(i,j)
        END IF

!       Influence of heatflux through snow on total forcing:
        zwsnew(i,j)   = zwsnow(i,j) + zdwsndt(i,j)*zdtdrhw
        IF (zwsnew(i,j).GT.zepsi) THEN

          zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (zth_low(i,j)-csnow_tmin) &
                                                  /(t0_melt          -csnow_tmin)
          zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))

          IF(zdwsndt(i,j)-zrs(i,j)-zrr(i,j).GT.0.0_ireals) THEN

!subs            wtot_snow(i,j,1,nx) = max(wtot_snow(i,j,1,nx) + zdwsndt(i,j)*zdtdrhw,0.0_ireals)
            wtot_snow(i,j,1,nx,ns) = max(wtot_snow(i,j,1,nx,ns) + zdwsndt(i,j)*zdtdrhw,0.0_ireals)

            zhm_snow(i,j,1) = zhm_snow(i,j,1) - (zdwsndt(i,j)-zrs(i,j)- &
                              zrr(i,j))*zdt/rho_i/2._ireals-  &
              zrs(i,j)*zdt/zrho_snowf/2._ireals- zrr(i,j)*zdt/rho_i/2._ireals
            zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_i +  &
              zrs(i,j)*zdt/zrho_snowf + zrr(i,j)*zdt/rho_i

!subs            rho_snow_mult(i,j,1,nx) = max(wtot_snow(i,j,1,nx)*rho_w/zdzh_snow(i,j,1),0.0_ireals)
            rho_snow_mult(i,j,1,nx,ns) = max(wtot_snow(i,j,1,nx,ns)*rho_w/zdzh_snow(i,j,1), &
              &                              0.0_ireals)
          ELSE

!subs            wtot_snow(i,j,1,nx) = max(wtot_snow(i,j,1,nx) + (zrs(i,j)+zrr(i,j))*zdtdrhw,0.0_ireals)
            wtot_snow(i,j,1,nx,ns) = max(wtot_snow(i,j,1,nx,ns) + (zrs(i,j)+zrr(i,j))*zdtdrhw, &
              &                          0.0_ireals)

            zhm_snow(i,j,1)  = zhm_snow(i,j,1) - zrs(i,j)*zdt/zrho_snowf/2._ireals- &
                               zrr(i,j)*zdt/rho_i/2._ireals
            zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + zrs(i,j)*zdt/zrho_snowf + zrr(i,j)*zdt/rho_i

!subs            IF(wtot_snow(i,j,1,nx) .GT. 0._ireals) THEN
!subs              rho_snow_mult(i,j,1,nx) = max(wtot_snow(i,j,1,nx)*rho_w/zdzh_snow(i,j,1),0.0_ireals)
            IF(wtot_snow(i,j,1,nx,ns) .GT. 0._ireals) THEN
              rho_snow_mult(i,j,1,nx,ns) = max(wtot_snow(i,j,1,nx,ns)*rho_w/zdzh_snow(i,j,1), &
                &                              0.0_ireals)

!subs              wtot_snow(i,j,1,nx) = max(wtot_snow(i,j,1,nx) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdtdrhw,0.0_ireals)
              wtot_snow(i,j,1,nx,ns) = max(wtot_snow(i,j,1,nx,ns) &
                &                      + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdtdrhw,0.0_ireals)

!subs              zhm_snow(i,j,1)  = zhm_snow(i,j,1) - (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_snow_mult(i,j,1,nx)/2.
!subs              zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_snow_mult(i,j,1,nx)
              zhm_snow(i,j,1)  = zhm_snow(i,j,1) - (zdwsndt(i,j)-zrs(i,j)-zrr(i,j)) &
                &                *zdt/rho_snow_mult(i,j,1,nx,ns)/2._ireals
              zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j)) &
                &                *zdt/rho_snow_mult(i,j,1,nx,ns)
            ELSE
!subs              rho_snow_mult(i,j,1,nx) = 0.0_ireals
              rho_snow_mult(i,j,1,nx,ns) = 0.0_ireals
              zdzh_snow(i,j,1) = 0.0_ireals
            END IF

          END IF
          zflag = 0._ireals
!subs          h_snow(i,j,nnow) = 0.
          h_snow(i,j,nnow,ns) = 0._ireals
          DO ksn = 1,ke_snow
!subs            h_snow(i,j,nnow) = h_snow(i,j,nnow) + zdzh_snow(i,j,ksn)
!subs            IF(wtot_snow(i,j,ksn,nx).LT.zepsi) zflag = 1._ireals
            h_snow(i,j,nnow,ns) = h_snow(i,j,nnow,ns) + zdzh_snow(i,j,ksn)
            IF(wtot_snow(i,j,ksn,nx,ns).LT.zepsi) zflag = 1._ireals
          END DO
          IF(zwsnow(i,j) .GT. zepsi) THEN
!subs            CALL normalize(h_snow(i,j,nnow),zdzm_snow(i,j,:),zdzh_snow(i,j,:), &
!subs              zhh_snow(i,j,:),wtot_snow(i,j,:,nx),ztsnow_mult(i,j,:),   &
!subs              wliq_snow(i,j,:,nx),rho_snow_mult(i,j,:,nx),i,j)
            CALL normalize(h_snow(i,j,nnow,ns),zdzm_snow(i,j,:),zdzh_snow(i,j,:), &
              zhh_snow(i,j,:),wtot_snow(i,j,:,nx,ns),ztsnow_mult(i,j,:),   &
              wliq_snow(i,j,:,nx,ns),rho_snow_mult(i,j,:,nx,ns),ke_snow)
          ELSE
            DO ksn = 1,ke_snow
!subs              rho_snow_mult(i,j,ksn,nx) = rho_snow_mult(i,j,1,nx)
!subs              ztsnow_mult(i,j,ksn) = t_s(i,j,nx)
!subs              wtot_snow(i,j,ksn,nx) = zwsnew(i,j)/ke_snow
!subs              zhh_snow(i,j,ksn) = -h_snow(i,j,nnow)/ke_snow*(ke_snow-ksn)
              rho_snow_mult(i,j,ksn,nx,ns) = rho_snow_mult(i,j,1,nx,ns)
              ztsnow_mult(i,j,ksn) = t_s(i,j,nx,ns)
              wtot_snow(i,j,ksn,nx,ns) = zwsnew(i,j)/REAL(ke_snow,ireals)
              zhh_snow(i,j,ksn) = -h_snow(i,j,nnow,ns)/REAL(ke_snow,ireals)* &
                         (REAL(ke_snow,ireals)-REAL(ksn,ireals))
            END DO

!subs            zhm_snow(i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals
            zhm_snow(i,j,1) = (-h_snow(i,j,nnow,ns) + zhh_snow(i,j,1))/2._ireals
            DO ksn = 2,ke_snow
              zhm_snow(i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals
            END DO

            !layer thickness betw. half levels of uppermost snow layer
!subs            zdzh_snow(i,j,1) = zhh_snow(i,j,1) + h_snow(i,j,nnow)
            zdzh_snow(i,j,1) = zhh_snow(i,j,1) + h_snow(i,j,nnow,ns)
            !layer thickness between soil surface and main level of uppermost layer
!subs            zdzm_snow(i,j,1) = zhm_snow(i,j,1) + h_snow(i,j,nnow)
            zdzm_snow(i,j,1) = zhm_snow(i,j,1) + h_snow(i,j,nnow,ns)

            DO ksn = 2,ke_snow
              zdzh_snow(i,j,ksn) = zhh_snow(i,j,ksn) - zhh_snow(i,j,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,j,ksn) = zhm_snow(i,j,ksn) - zhm_snow(i,j,ksn-1) ! layer thickness betw. main levels
            ENDDO
          END IF

!         heat conductivity of snow as funtion of water content
! BR      zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i,j)))
!
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
          DO ksn = 1, ke_snow
!subs            zalas_mult(i,j,ksn) = 2.22_ireals*(rho_snow_mult(i,j,ksn,nx)/rho_i)**1.88_ireals
            zalas_mult(i,j,ksn) = 2.22_ireals*(rho_snow_mult(i,j,ksn,nx,ns)/rho_i)**1.88_ireals
          END DO

! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
!            The water vapour transport associated conductivity is not included.

!        zalas   = 0.023_ireals+(2.290_ireals-0.023_ireals)* &
!                               (7.750E-05_ireals*rho_snow(i,j,nx) + &
!                                1.105E-06_ireals*prho_snow(i,j,nx)**2)

!em          zgsb(i,j) = zalas_mult(i,j,1)*(ztsnow_mult(i,j,1) - zts(i,j))/zdz_snow_fl(i,j)
          zgsb(i,j) = (zalas_mult(i,j,ke_snow)+zalam(i,j,1))/2._ireals * &
!subs                      (ztsnow_mult(i,j,ke_snow) - t_so(i,j,1,nx))/(zdzh_snow(i,j,ke_snow)+zdzhs(1))*2._ireals
                      (ztsnow_mult(i,j,ke_snow) - t_so(i,j,1,nx,ns))/(zdzh_snow(i,j,ke_snow) &
                      +zdzhs(1))*2._ireals
        END IF

        ! total forcing for uppermost soil layer
!em        zfor_s(i,j) = zrnet_s + zshfl_s + zlhfl_s + zsprs
        zfor_s(i,j) = zrnet_s + zshfl_s + zlhfl_s + zsprs*(1._ireals - zf_snow(i,j))       &
                      + zf_snow(i,j) * (1._ireals-ztsnow_pm(i,j)) * zgsb(i,j)

!em        IF (zwsnew(i,j) .GT. 5.0E-03_ireals) THEN
!em          IF(zextinct(i,j,1).ne.0.0) THEN
          IF(zwsnew(i,j) .GT. zepsi) THEN
!subs            zrnet_snow = sobs(i,j) * (1.0_ireals - EXP(-zextinct(i,j,1)*zdzm_snow(i,j,1))) + zthsnw(i,j)
            zrnet_snow = sobs(i,j,ns) * (1.0_ireals - EXP(-zextinct(i,j,1)*zdzm_snow(i,j,1))) &
              &          + zthsnw(i,j)
          ELSE
!subs            zrnet_snow = sobs(i,j) + zthsnw(i,j)
            zrnet_snow = sobs(i,j,ns) + zthsnw(i,j)
          END IF
          zshfl_snow = zrhoch(i,j)*cp_d*(zth_low(i,j) - ztsnow_mult(i,j,1))
          zlhfl_snow = lh_s*zversn(i,j)/MAX(zf_snow(i,j), zepsi)
          zfor_snow_mult(i,j)  = zf_snow(i,j)*(zrnet_snow + zshfl_snow + zlhfl_snow + zsprs)
!em          zfor_total(i,j) = zfor_s(i,j) + zfor_snow_mult(i,j)
!        END IF        ! points with snow cover only

      END IF          ! land-points only
    END DO
   END DO

  ELSE

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)
        zgstr =   sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i,j))* &
!subs                  zts(i,j) + zf_snow(i,j)*ztsnow(i,j) )**4 + thbs(i,j)
                  zts(i,j) + zf_snow(i,j)*ztsnow(i,j) )**4 + thbs(i,j,ns)
        zthsnw(i,j) = - sigma*(1._ireals - Ctalb)*ztsnow(i,j)**4 + zgstr
        zthsoi(i,j) = - sigma*(1._ireals - Ctalb)*zts(i,j)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux
!subs        zrnet_s = (1._ireals - zf_snow(i,j))*(sobs(i,j)+zthsoi(i,j))
        zrnet_s = (1._ireals - zf_snow(i,j))*(sobs(i,j,ns)+zthsoi(i,j))
        zshfl_s = (1._ireals - zf_snow(i,j))*cp_d*zrhoch(i,j)*  &
                                                      (zth_low(i,j) - zts(i,j))
        zlhfl_s = (zts_pm(i,j)*lh_v + (1._ireals-zts_pm(i,j))*lh_s)*zverbo(i,j)
        zsprs      = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i,j)*zrs(i,j) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs        = - lh_f*zrs(i,j)


!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!!$print *,ns,'zsprs, snow fall on soil with T>T0',zsprs,zrs(i,j)
!!$          zdwidt (i,j) = zdwidt (i,j) + zrs(i,j)
!!$          zdwsndt(i,j) = zdwsndt(i,j) - zrs(i,j)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
!subs          zwimax       = cwimax_ml*(1._ireals + plcov(i,j)*5._ireals)
          zwimax       = cwimax_ml*(1._ireals + plcov(i,j,ns)*5._ireals)
          zwinstr      = zwin(i,j) + zdwidt(i,j)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zro        = (zwinstr - zwimax)*zrhwddt
            zdwidt(i,j)= zdwidt(i,j) - zro
!subs            runoff_s(i,j)  = runoff_s(i,j) + zro*zroffdt
            runoff_s(i,j,ns)  = runoff_s(i,j,ns) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF (zwsnow(i,j) == 0.0_ireals .AND.                            &
               (1._ireals-ztsnow_pm(i,j))*zrr(i,j) > 0.0_ireals) THEN
          zsprs   = lh_f*zrr(i,j)
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!!$print *,ns,'zsprs, freezing of rain falling on soil with Ts < T0',zsprs,zrr(i,j)
!!$          zdwidt (i,j) = zdwidt (i,j) - zrr(i,j)
!!$          zdwsndt(i,j) = zdwsndt(i,j) + zrr(i,j)
        END IF

!       Influence of heatflux through snow on total forcing:
        zwsnew(i,j)   = zwsnow(i,j) + zdwsndt(i,j)*zdtdrhw
        IF (zwsnew(i,j).GT.zepsi) THEN
!         heat conductivity of snow as funtion of water content
! BR      zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i,j)))
!
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
!subs          zalas  = 2.22_ireals*(rho_snow(i,j,nx)/rho_i)**1.88_ireals
          zalas  = 2.22_ireals*(rho_snow(i,j,nx,ns)/rho_i)**1.88_ireals

! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
!            The water vapour transport associated conductivity is not included.

!        zalas   = 0.023_ireals+(2.290_ireals-0.023_ireals)* &
!                               (7.750E-05_ireals*rho_snow(i,j,nx) + &
!                                1.105E-06_ireals*prho_snow(i,j,nx)**2)

          zgsb(i,j) = zalas*(ztsnow(i,j) - zts(i,j))/zdz_snow_fl(i,j)
        END IF

        ! total forcing for uppermost soil layer
        zfor_s(i,j) = zrnet_s + zshfl_s + zlhfl_s + zsprs       &
                         + zf_snow(i,j) * (1._ireals-ztsnow_pm(i,j)) * zgsb(i,j)
      END IF          ! land-points only
    END DO
   END DO

 ENDIF ! lmulti_snow

  !>JH WRITE(*,*) 'TERRA-DIAG END II.5 +++++++++> ',t_so(1,1,:,nnew,:)

!------------------------------------------------------------------------------
! Section II.6: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

IF(lmulti_snow) THEN
  zswitch = 5.0E-03_ireals
ELSE
  zswitch = 1.0E+06_ireals
END IF

IF(lmulti_snow) THEN
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF(zwsnew(i,j) .LT. zswitch) THEN
          IF(zfor_snow_mult(i,j)*zdt > zwsnew(i,j)*rho_w*lh_f) THEN
            zfor_snow_mult(i,j) = zfor_snow_mult(i,j) - zwsnew(i,j)*rho_w*lh_f/zdt
            zwsnew(i,j) = 0._ireals
            DO ksn = 0, ke_snow
!subs              ztsnown_mult(i,j,ksn) = t_so(i,j,0,nx)
              ztsnown_mult(i,j,ksn) = t_so(i,j,0,nx,ns)
            END DO
            zfor_s(i,j) = zfor_s(i,j) + zfor_snow_mult(i,j)
          ELSEIF(zwsnew(i,j) .GE. zepsi .AND. zfor_snow_mult(i,j) .GT. 0._ireals) THEN
            DO ksn = 1, ke_snow
              ztsnown_mult(i,j,ksn) = ztsnow_mult(i,j,ksn) + &
!subs                zfor_snow_mult(i,j)*zdt/(chc_i*wtot_snow(i,j,ksn,nx))/rho_w/ke_snow
                zfor_snow_mult(i,j)*zdt/(chc_i*wtot_snow(i,j,ksn,nx,ns))/rho_w/REAL(ke_snow,ireals)
            END DO
            zfor_snow_mult(i,j) = 0._ireals
            ztsnown_mult(i,j,0) = ztsnown_mult(i,j,1)
          END IF
        END IF
      END IF  ! land-points only
    END DO
  END DO
END IF

DO kso = 2,ke_soil
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        ! for heat conductivity: zalam is now 3D
        zakb = zalam(i,j,kso)/zroc(i,j,kso)
        zaga(i,j,kso) = -zalfa*zdt*zakb/(zdzhs(kso)*zdzms(kso))
        zagc(i,j,kso) = -zalfa*zdt*zakb/(zdzhs(kso)*zdzms(kso+1))
        zagb(i,j,kso) = 1._ireals - zaga(i,j,kso) - zagc(i,j,kso)
!subs        zagd(i,j,kso) = t_so(i,j,kso,nx) +                                     &
!subs               (1._ireals - zalfa)*( - zaga(i,j,kso)/zalfa*t_so(i,j,kso-1,nx)+ &
!subs               (zaga(i,j,kso)/zalfa + zagc(i,j,kso)/zalfa)*t_so(i,j,kso,nx) -  &
!subs                zagc(i,j,kso)/zalfa*t_so(i,j,kso+1,nx)  )
        zagd(i,j,kso) = t_so(i,j,kso,nx,ns) +                                     &
               (1._ireals - zalfa)*( - zaga(i,j,kso)/zalfa*t_so(i,j,kso-1,nx,ns)+ &
               (zaga(i,j,kso)/zalfa + zagc(i,j,kso)/zalfa)*t_so(i,j,kso,nx,ns) -  &
                zagc(i,j,kso)/zalfa*t_so(i,j,kso+1,nx,ns)  )
      END IF  ! land-points only
    END DO
  END DO
END DO        ! soil layers

  !>JH WRITE(*,*) 'TERRA-DIAG END II.6 - 1 +++++++++> ',t_so(1,1,:,nnew,:)

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
    IF (llandmask(i,j,ns)) THEN          ! land-points only
      ! for heat conductivity: zalam is now 3D: here we need layer 1
      zakb = zalam(i,j,1)/zroc(i,j,1)
      zaga(i,j,  1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(1))
      zagc(i,j,  1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(2))
      zagb(i,j,  1) = 1._ireals - zaga(i,j,1) - zagc(i,j,1)
!subs      zagd(i,j,  1) = t_so(i,j,1,nx) + (1._ireals - zalfa)* (                 &
!subs                      - zaga(i,j,1)/zalfa * t_s(i,j,nx) +                     &
!subs                      (zaga(i,j,1) + zagc(i,j,1))/zalfa * t_so(i,j,1,nx) -    &
!subs                       zagc(i,j,1)/zalfa * t_so(i,j,2,nx)   )
      zagd(i,j,  1) = t_so(i,j,1,nx,ns) + (1._ireals - zalfa)* (                 &
                      - zaga(i,j,1)/zalfa * t_s(i,j,nx,ns) +                     &
                      (zaga(i,j,1) + zagc(i,j,1))/zalfa * t_so(i,j,1,nx,ns) -    &
                       zagc(i,j,1)/zalfa * t_so(i,j,2,nx,ns)   )
      zaga(i,j,0)    = 0.0_ireals
      zagb(i,j,0)    = zalfa
      zagc(i,j,0)    = -zalfa
      zagd(i,j,0)    = zdzms(1) * zfor_s(i,j)/zalam(i,j,1)+(1._ireals-zalfa)* &
!subs                      (t_so(i,j,1,nx) - t_s(i,j,nx))
                      (t_so(i,j,1,nx,ns) - t_s(i,j,nx,ns))
      zaga(i,j,ke_soil+1) = 0.0_ireals
      zagb(i,j,ke_soil+1) = 1.0_ireals
      zagc(i,j,ke_soil+1) = 0.0_ireals
!subs      zagd(i,j,ke_soil+1) = t_so(i,j,ke_soil+1,nx)
      zagd(i,j,ke_soil+1) = t_so(i,j,ke_soil+1,nx,ns)
!!$IF(i.eq.1.and.j.eq.1) WRITE(*,*) 'END II.6 - 2 +++++++++> zaga ',zaga(1,1,0),zaga(1,1,1),zaga(1,1,ke_soil+1)
!!$IF(i.eq.1.and.j.eq.1) WRITE(*,*) 'END II.6 - 2 +++++++++> zagb ',zagb(1,1,0),zagb(1,1,1),zagb(1,1,ke_soil+1)
!!$IF(i.eq.1.and.j.eq.1) WRITE(*,*) 'END II.6 - 2 +++++++++> zagc ',zagc(1,1,0),zagc(1,1,1),zagc(1,1,ke_soil+1)
!!$IF(i.eq.1.and.j.eq.1) WRITE(*,*) 'END II.6 - 2 +++++++++> zagd ',zagd(1,1,0),zagd(1,1,1),zagd(1,1,ke_soil+1)
!!$IF(i.eq.1.and.j.eq.1) WRITE(*,*) 'END II.6 - 2 +++++++++> zfor_s(i,j) ',zfor_s(i,j),zalam(i,j,1),zalfa,zdzms(1)
    END IF          ! land-points only
  END DO
END DO


DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
    IF (llandmask(i,j,ns)) THEN          ! land-points only
      zagc(i,j,0) = zagc(i,j,0)/zagb(i,j,0)
      zagd(i,j,0) = zagd(i,j,0)/zagb(i,j,0)
    END IF          ! land-points only
  END DO
END DO


DO kso=1,ke_soil
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        zzz = 1._ireals/(zagb(i,j,kso) - zaga(i,j,kso)*zagc(i,j,kso-1))
        zagc(i,j,kso) = zagc(i,j,kso) * zzz
        zagd(i,j,kso) = (zagd(i,j,kso) - zaga(i,j,kso)*zagd(i,j,kso-1)) * zzz
      END IF          ! land-points only
    END DO
  END DO
END DO                ! soil layers

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
    IF (llandmask(i,j,ns)) THEN          ! land-points only
      zage(i,j,ke_soil+1) = (zagd(i,j,ke_soil+1) - zaga(i,j,ke_soil+1)*       &
                             zagd(i,j,ke_soil))/                              &
                            (zagb(i,j,ke_soil+1) - zaga(i,j,ke_soil+1)*       &
                             zagc(i,j,ke_soil))
    END IF          ! land-points only
  END DO
END DO

  !>JH WRITE(*,*) 'TERRA-DIAG END II.6 - 2.1 +++++++++> ',t_so(1,1,:,nnew,:),'zaga ',zaga(1,1,:),'zagb ',zagb(1,1,:), &
  !  'zagc ',zagc(1,1,:),'zagd ',zagd(1,1,:), 'zage ',zage(1,1,:)

DO kso = ke_soil,0,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        zage(i,j,kso)     = zagd(i,j,kso) - zagc(i,j,kso)*zage(i,j,kso+1)
      ! The surface temperature computed by t_so(i,j,0,nnew)=zage(i,j,0) is
      ! presently unused
!subs        t_so(i,j,kso,nnew) = zage(i,j,kso)
        t_so(i,j,kso,nnew,ns) = zage(i,j,kso)
!IF(i.eq.340 .and. j.eq.186 .and. kso.eq.1) &

!>JH IF(i.eq.i_loc.and. j.eq.j_loc.and. kso.eq.1 .and. my_cart_id.eq.isub) &
!print *,ns,kso,'after sweep',t_so(i,j,kso,nnew,ns),zfor_s(i,j)

  if (i.eq.40 .and. j.eq.1) then
print*, "SFC-DIAGNOSIS TERRA II ",ke,dt,nsubs1,ntstep

  print*, "t_so",t_so(i,j,:,:,:)
  print*, "w_so",w_so(i,j,:,:,:)

end if


      END IF          ! land-points only
    END DO
  END DO
END DO                ! soil layers

 
DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
!subs      t_so(i,j,ke_soil+1,nnew) = zage(i,j,ke_soil+1) ! climate value, unchanged
    IF (llandmask(i,j,ns)) THEN          ! land-points only
      t_so(i,j,ke_soil+1,nnew,ns) = zage(i,j,ke_soil+1) ! climate value, unchanged
    END IF          ! land-points only
  END DO
END DO

  !>JH WRITE(*,*) 'TERRA-DIAG END II.6 - 3 +++++++++> ',t_so(1,1,:,nnew,:),'zage ',zage(1,1,:)

IF (lmulti_snow) THEN
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (zwsnew(i,j) .GE. zswitch &
          & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
!subs          zrocs(i,j) = wliq_snow(i,j,1,nx)/zdzh_snow(i,j,1)*rho_w*chc_w + &
!subs            (rho_snow_mult(i,j,1,nx)/rho_w -wliq_snow(i,j,1,nx)/zdzh_snow(i,j,1))*rho_w*chc_i
          zrocs(i,j) = wliq_snow(i,j,1,nx,ns)/zdzh_snow(i,j,1)*rho_w*chc_w + &
            (rho_snow_mult(i,j,1,nx,ns)/rho_w -wliq_snow(i,j,1,nx,ns)/zdzh_snow(i,j,1))*rho_w*chc_i
          zakb      = zalas_mult(i,j,1)/zrocs(i,j)  !(chc_i*rho_snow(1,nx))
          zaga(i,j,1) = 0.0_ireals
          zagc(i,j,1) = -zdt*zakb/(zdzh_snow(i,j,1)*zdzm_snow(i,j,2))
          zagb(i,j,1) = 1._ireals - zagc(i,j,1)
          zagd(i,j,1) = ztsnow_mult(i,j,1) + zdt/zdzh_snow(i,j,1)*zfor_snow_mult(i,j)/zrocs(i,j)

!subs          zrocs(i,j) = wliq_snow(i,j,ke_snow,nx)/zdzh_snow(i,j,ke_snow)*rho_w*chc_w + &
!subs            (rho_snow_mult(i,j,ke_snow,nx)/rho_w - &
!subs            wliq_snow(i,j,ke_snow,nx)/zdzh_snow(i,j,ke_snow))*rho_w*chc_i
          zrocs(i,j) = wliq_snow(i,j,ke_snow,nx,ns)/zdzh_snow(i,j,ke_snow)*rho_w*chc_w + &
            (rho_snow_mult(i,j,ke_snow,nx,ns)/rho_w - &
            wliq_snow(i,j,ke_snow,nx,ns)/zdzh_snow(i,j,ke_snow))*rho_w*chc_i
          zakb = zalas_mult(i,j,ke_snow)/zrocs(i,j)  !(chc_i*rho_snow(1,nx))
          zaga(i,j,ke_snow) = -zdt*zakb/(zdzh_snow(i,j,ke_snow)*zdzm_snow(i,j,ke_snow))
          zagb(i,j,ke_snow) = 1.0_ireals - zaga(i,j,ke_snow)
          zagc(i,j,ke_snow) = 0.0_ireals
          zagd(i,j,ke_snow) = ztsnow_mult(i,j,ke_snow) &
                              - zdt/zdzh_snow(i,j,ke_snow)*zgsb(i,j)*zf_snow(i,j)/zrocs(i,j)
        END IF
      END IF          ! land-points only
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (zwsnew(i,j) .GE. zswitch &
          & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
          zagc(i,j,1) = zagc(i,j,1)/zagb(i,j,1)
          zagd(i,j,1) = zagd(i,j,1)/zagb(i,j,1)
        END IF
      END IF          ! land-points only
    END DO
  END DO

  DO ksn = 2, ke_snow-1
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN          ! land-points only
        IF (llandmask(i,j,ns)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zswitch &
            & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
!subs            zrocs(i,j) = wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn)*rho_w*chc_w + &
!subs              (rho_snow_mult(i,j,ksn,nx)/rho_w - wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn))*rho_w*chc_i
            zrocs(i,j) = wliq_snow(i,j,ksn,nx,ns)/zdzh_snow(i,j,ksn)*rho_w*chc_w + &
              (rho_snow_mult(i,j,ksn,nx,ns)/rho_w - wliq_snow(i,j,ksn,nx,ns)/zdzh_snow(i,j,ksn)) &
              *rho_w*chc_i
            zakb = zalas_mult(i,j,ksn)/zrocs(i,j) !(chc_i*rho_snow(ksn,nx))
            zaga(i,j,ksn) = -zdt*zakb/(zdzh_snow(i,j,ksn)*zdzm_snow(i,j,ksn))
            zagc(i,j,ksn) = -zdt*zakb/(zdzh_snow(i,j,ksn)*zdzm_snow(i,j,ksn+1))
            zagb(i,j,ksn) = 1._ireals - zaga(i,j,ksn) - zagc(i,j,ksn)
            zagd(i,j,ksn) = ztsnow_mult(i,j,ksn)
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! snow layers

  DO ksn=2,ke_snow-1
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN          ! land-points only
        IF (llandmask(i,j,ns)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zswitch &
            & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
            zzz = 1._ireals/(zagb(i,j,ksn) - zaga(i,j,ksn)*zagc(i,j,ksn-1))
            zagc(i,j,ksn) = zagc(i,j,ksn) * zzz
            zagd(i,j,ksn) = (zagd(i,j,ksn) - zaga(i,j,ksn)*zagd(i,j,ksn-1)) * zzz
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! snow layers

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF(zwsnew(i,j) .GE. zswitch &
          & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
          zage(i,j,ke_snow) = (zagd(i,j,ke_snow) - zaga(i,j,ke_snow)*       &
                                          zagd(i,j,ke_snow-1))/                           &
                                          (zagb(i,j,ke_snow) - zaga(i,j,ke_snow)*       &
                                          zagc(i,j,ke_snow-1))
          ztsnown_mult(i,j,ke_snow) = MAX(zage(i,j,ke_snow),173.15_ireals)
        END IF
      END IF          ! land-points only
    END DO
  END DO

  DO ksn = ke_snow-1,1,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN          ! land-points only
        IF (llandmask(i,j,ns)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zswitch &
            & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
            zage(i,j,ksn)     = zagd(i,j,ksn) - zagc(i,j,ksn)*zage(i,j,ksn+1)
            ! The surface temperature computed by t_so(0,nnew)=zage(0) is
            ! presently unused
            ztsnown_mult(i,j,ksn) = MAX(zage(i,j,ksn),173.15_ireals)
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! snow layers

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF(zwsnew(i,j) .GE. zswitch &
          & .OR. (zfor_snow_mult(i,j) .LE. 0._ireals .AND. zwsnew(i,j) .GT. zepsi)) THEN
          ztsnown_mult(i,j,0) = ztsnown_mult(i,j,1)
        END IF
      END IF          ! land-points only
    END DO
  END DO
END IF

 !>JH WRITE(*,*) 'TERRA-DIAG END II.6 - 4 +++++++++> ',t_so(1,1,:,nnew,:)

IF(lmelt) THEN
  IF(.NOT.lmelt_var) THEN
    DO kso = 1,ke_soil
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN ! land points only,
          IF (llandmask(i,j,ns)) THEN ! land points only,
            IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
              ! melting or freezing of soil water
!subs              zenergy   = zroc(i,j,kso)*zdzhs(kso)*(t_so(i,j,kso,nnew) - t0_melt)
              zenergy   = zroc(i,j,kso)*zdzhs(kso)*(t_so(i,j,kso,nnew,ns) - t0_melt)
              zdwi_max  = - zenergy/(lh_f*rho_w)
              zdelwice  = zdwi_max
!subs              zwso_new  = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
              zwso_new  = w_so(i,j,kso,nx,ns) + zdt*zdwgdt(i,j,kso)/rho_w
              IF (zdelwice.LT.0.0_ireals)                                    &
!subs                        zdelwice = - MIN( - zdelwice,w_so_ice(i,j,kso,nx))
                        zdelwice = - MIN( - zdelwice,w_so_ice(i,j,kso,nx,ns))
              IF (zdelwice.GT.0.0_ireals)                                    &
!subs                        zdelwice =   MIN(   zdelwice,zwso_new - w_so_ice(i,j,kso,nx))
                        zdelwice =   MIN(   zdelwice,zwso_new - w_so_ice(i,j,kso,nx,ns))
!subs              w_so_ice(i,j,kso,nnew) = w_so_ice(i,j,kso,nx) + zdelwice
              w_so_ice(i,j,kso,nnew,ns) = w_so_ice(i,j,kso,nx,ns) + zdelwice
!             At this point we have 0.0 LE w_so_ice(i,j,kso,nnew) LE w_so(i,j,kso,nnew)
!             If we have 0.0 LT w_so_ice(i,j,kso,nnew) LT w_so(i,j,kso,nnew)
!             the resulting new temperature has to be equal to t0_melt.
!             If not all energy available can be used to melt/freeze soil water, the
!             following line corrects the temperature. It also applies to cases
!             where all soil water is completely ice or water, respectively.
!subs              t_so(i,j,kso,nnew) = t0_melt + (zdelwice - zdwi_max)*(lh_f*rho_w)/     &
              t_so(i,j,kso,nnew,ns) = t0_melt + (zdelwice - zdwi_max)*(lh_f*rho_w)/     &
                                            (zroc(i,j,kso)*zdzhs(kso))
            END IF  ! m_styp > 2
          END IF  ! land-points only
        END DO
      END DO
    END DO        ! soil layers
  ELSE IF(lmelt_var) THEN
    DO kso = 1,ke_soil
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN ! land points only,
          IF (llandmask(i,j,ns)) THEN ! land points only,
            IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
            !>JH WRITE(*,*) 'TERRA-DIAG END II.6 4.0 +++++++++> ',t_so(1,1,:,nnew,:)
              ztx      = t0_melt
              zw_m(i,j)     = zporv(i,j)*zdzhs(kso)
!subs              IF(t_so(i,j,kso,nnew).LT.(t0_melt-zepsi)) THEN
              IF(t_so(i,j,kso,nnew,ns).LT.(t0_melt-zepsi)) THEN
                zaa    = g*zpsis(i,j)/lh_f
!subs                zw_m(i,j) = zw_m(i,j)*((t_so(i,j,kso,nnew) - t0_melt)/(t_so(i,j,kso,nnew)*  &
                zw_m(i,j) = zw_m(i,j)*((t_so(i,j,kso,nnew,ns) - t0_melt)/(t_so(i,j,kso,nnew,ns)*  &
                         zaa))**(-zedb(i,j))
!subs                zliquid= MAX(zepsi,w_so(i,j,kso,nx) -  w_so_ice(i,j,kso,nx))
                zliquid= MAX(zepsi,w_so(i,j,kso,nx,ns) -  w_so_ice(i,j,kso,nx,ns))
                znen   = 1._ireals-zaa*(zporv(i,j)*zdzhs(kso)/zliquid)**zb_por(i,j)
                ztx    = t0_melt/znen
              ENDIF
              ztx      = MIN(t0_melt,ztx)
!subs              zenergy  = zroc(i,j,kso)*zdzhs(kso)*(t_so(i,j,kso,nnew)-ztx)
              zenergy  = zroc(i,j,kso)*zdzhs(kso)*(t_so(i,j,kso,nnew,ns)-ztx)
              zdwi_max = - zenergy/(lh_f*rho_w)
              zdelwice = zdwi_max
!subs              zwso_new  = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
!subs              zargu = zwso_new - zw_m(i,j) - w_so_ice(i,j,kso,nx)
              zwso_new  = w_so(i,j,kso,nx,ns) + zdt*zdwgdt(i,j,kso)/rho_w
              zargu = zwso_new - zw_m(i,j) - w_so_ice(i,j,kso,nx,ns)
              IF (zdelwice.LT.0.0_ireals) zdelwice =                           &
!subs                      - MIN( - zdelwice,MIN(-zargu,w_so_ice(i,j,kso,nx)))
                      - MIN( - zdelwice,MIN(-zargu,w_so_ice(i,j,kso,nx,ns)))
              IF (zdelwice.GT.0.0_ireals) zdelwice =                           &
                        MIN(   zdelwice,MAX( zargu,0.0_ireals))
!subs              w_so_ice(i,j,kso,nnew) = w_so_ice(i,j,kso,nx) + zdelwice
              w_so_ice(i,j,kso,nnew,ns) = w_so_ice(i,j,kso,nx,ns) + zdelwice
!             At this point we have 0.0 LE w_so_ice(i,j,kso,nnew) LE zwso_new
!             If we have 0.0 LT w_so_ice(i,j,kso,nnew) LT zwso_new
!             the resulting new temperature has to be equal to ztx.
!             If not all energy available can be used to melt/freeze soil water,
!             the following line corrects the temperature. It also applies
!             to cases without any freezing/melting, in these cases the
!             original temperature is reproduced.
!subs              t_so(i,j,kso,nnew) = ztx + (zdelwice - zdwi_max)*       &
            !>JH WRITE(*,*) 'TERRA-DIAG END II.6 4.1 +++++++++> ',t_so(1,1,:,nnew,ns),'lhf ',lh_f, 'zroc ',zroc(i,j,:)
              t_so(i,j,kso,nnew,ns) = ztx + (zdelwice - zdwi_max)*       &
                                   (lh_f*rho_w)/(zroc(i,j,kso)*zdzhs(kso))
!IF(i.eq.340 .and. j.eq.186 .and. kso.eq.1) &
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. kso.eq.1 .and. my_cart_id.eq.isub) &
! print *,ns,kso,'after lmelt',t_so(i,j,kso,nnew,ns)
            END IF                  ! m_stpy > 2
          END IF                    ! land-points only
        END DO
      END DO
    ENDDO
  ENDIF   ! lmelt_var
END IF ! lmelt

  !>JH WRITE(*,*) 'TERRA-DIAG END II.6 +++++++++> '

!------------------------------------------------------------------------------
! Section II.7: Energy budget and temperature prediction at snow-surface
!------------------------------------------------------------------------------
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        ! next line has to be changed if a soil surface temperature is
        ! predicted by the heat conduction equation
!subs        zdtsdt (i,j) = (t_so(i,j,1,nnew) - zts(i,j))*z1d2dt
!subs        ztsn   (i,j) =  t_so(i,j,1,nnew)
        zdtsdt (i,j) = (t_so(i,j,1,nnew,ns) - zts(i,j))*z1d2dt
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
! print *,ns,'before snow',t_so(i,j,1,nnew,ns),zdt*zdtsdt   (i,j)
        ztsn   (i,j) =  t_so(i,j,1,nnew,ns)
        IF(.NOT. lmulti_snow)  &
          ztsnown(i,j) = ztsn(i,j)       ! default setting
        zwsnew(i,j)       = zwsnow(i,j) + zdwsndt(i,j)*zdtdrhw

        ! forcing contributions for snow formation of dew and rime are
        ! contained in ze_ges, heat fluxes must not be multiplied by
        ! snow covered fraction

        IF(.NOT. lmulti_snow) THEN

!subs          zrnet_snow = sobs(i,j) + zthsnw(i,j)
          zrnet_snow = sobs(i,j,ns) + zthsnw(i,j)
          zshfl_snow = zrhoch(i,j)*cp_d*(zth_low(i,j) - ztsnow(i,j))
          zlhfl_snow = lh_s*zversn(i,j)
          zfor_snow  = zrnet_snow + zshfl_snow + zlhfl_snow

          ! forecast of snow temperature Tsnow
          IF (ztsnow(i,j) < t0_melt .AND. zwsnew(i,j) > zepsi) THEN
            ztsnown(i,j) = ztsnow(i,j) + zdt*2._ireals*(zfor_snow - zgsb(i,j))  &
                           /zrocs(i,j) - ( ztsn(i,j) - zts(i,j) )
            ! implicit formulation
! BR        zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnew(i,j)))
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
!subs            zalas  = 2.22_ireals*(rho_snow(i,j,nx)/rho_i)**1.88_ireals
            zalas  = 2.22_ireals*(rho_snow(i,j,nx,ns)/rho_i)**1.88_ireals

            ztsnow_im    = - zrhoch(i,j) * (cp_d + zdqvtsnow(i,j) * lh_s)       &
                                         - zalas/zdz_snow_fl(i,j)
            zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i,j))
            ztsnown(i,j) = ztsnow(i,j) + (ztsnown(i,j)-ztsnow(i,j))/zfak
          END IF

!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!!$print *,ns,'zfor_snow',zfor_snow,zrnet_snow , zshfl_snow , zlhfl_snow
!!$IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!!$print *,ns,'zshfl_snow', zshfl_snow , zrhoch(i,j),zth_low(i,j) , ztsnow(i,j)
!!$IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!!$print *,ns,'zlhfl_snow', zlhfl_snow , zversn(i,j)
!!$IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!!$print *,ns,'new snow temp',zfor_snow,ztsnown(i,j)
!!$          zdtsnowdt(i,j) = (ztsnown(i,j) - ztsnow(i,j))*z1d2dt
        ENDIF

      END IF          ! land-points only
    END DO
  END DO

  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
      DO   j = jstarts, jends
        DO i = istarts, iends
!subs          IF (llandmask(i,j)) THEN          ! land-points only
          IF (llandmask(i,j,ns)) THEN          ! land-points only
            IF (zwsnew(i,j) > zepsi) THEN
              zdtsnowdt_mult(i,j,ksn) = (ztsnown_mult(i,j,ksn) - ztsnow_mult(i,j,ksn))*z1d2dt
            END IF
          ENDIF
        END DO
      END DO
    ENDDO
  ENDIF

  !>JH WRITE(*,*) 'TERRA-DIAG END II.7 +++++++++> '

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

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
    IF (llandmask(i,j,ns)) THEN          ! land-points only

      zwsnn  (i,j)  = zwsnow(i,j) + zdtdrhw*zdwsndt(i,j)
      zwsnew (i,j)  = zwsnn(i,j)
      ztsnownew     = ztsnown(i,j)
      ze_avail      = 0.0_ireals
      ze_total      = 0.0_ireals
      zfr_melt      = 0.0_ireals

      IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
        ! first case: T_snow > t0_melt: melting from above
        ! ----------
!subs        IF (ztsnown(i,j) > t0_melt .AND. t_so(i,j,1,nnew) < t0_melt ) THEN
        IF (ztsnown(i,j) > t0_melt .AND. t_so(i,j,1,nnew,ns) < t0_melt ) THEN
          zdwsnm(i,j)   = zwsnew(i,j)*.5_ireals*(ztsnown(i,j) - (t0_melt - zepsi))/ &
                          (.5_ireals* (zts(i,j) - (t0_melt - zepsi)) - lh_f/chc_i)
          zdwsnm(i,j)   = zdwsnm(i,j)*z1d2dt*rho_w 
          zdwsndt(i,j)  = zdwsndt (i,j) + zdwsnm(i,j)
          ztsnownew       = t0_melt - zepsi
          zdtsnowdt(i,j)  = zdtsnowdt(i,j) + (ztsnownew - ztsnown(i,j))*z1d2dt
!subs          runoff_s (i,j)= runoff_s(i,j) - zdwsnm(i,j)*zroffdt
          runoff_s (i,j,ns)= runoff_s(i,j,ns) - zdwsnm(i,j)*zroffdt
        ENDIF ! melting from above

!subs        IF (t_so(i,j,1,nnew) .gt. t0_melt) THEN
        IF (t_so(i,j,1,nnew,ns) .gt. t0_melt) THEN
          !second case:  temperature of uppermost soil layer > t0_melt. First a
          !-----------   heat redistribution is performed. As a second step,
          !              melting of snow is considered.
          ! a) Heat redistribution
          ztsnew = t0_melt + zepsi
          ztsnownew      = ztsn(i,j) + ztsnown(i,j) - ztsnew +  &
               2._ireals*(ztsn(i,j) - ztsnew)*zroc(i,j,1)*zdzhs(1)/zrocs(i,j)
!subs          zdtsdt(i,j)    = zdtsdt(i,j) + (ztsnew - t_so(i,j,1,nnew))*z1d2dt
          zdtsdt(i,j)    = zdtsdt(i,j) + (ztsnew - t_so(i,j,1,nnew,ns))*z1d2dt
          zdtsnowdt(i,j) = zdtsnowdt(i,j) + (ztsnownew - ztsnown(i,j))*z1d2dt
          ! b) Melting of snow (if possible)
          IF (ztsnownew > t0_melt) THEN
            ze_avail     = 0.5_ireals*(ztsnownew - t0_melt)*zrocs(i,j)
            ze_total     = lh_f*zwsnew(i,j)*rho_w
            zfr_melt     = MIN(1.0_ireals,ze_avail/ze_total)
            zdtsnowdt(i,j)= zdtsnowdt(i,j) + (t0_melt - ztsnownew)*z1d2dt
            zdelt_s      = MAX(0.0_ireals,(ze_avail - ze_total)/(zroc(i,j,1)* &
                                                                    zdzhs(1)))
            zdtsdt(i,j)  = zdtsdt(i,j) + zdelt_s*z1d2dt
!IF(i.eq.340 .and. j.eq.186) &
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
! print *,ns,'in snow',t_so(i,j,1,nnew,ns),zdt*zdtsdt   (i,j)
            ! melted snow is allowed to penetrate the soil (up to field
            ! capacity), if the soil type is neither ice nor rock (zrock = 0);
            ! else it contributes to surface run-off;
            ! fractional water content of the first soil layer determines
            ! a reduction factor which controls additional run-off
            zdwsnm(i,j)   = zfr_melt*zwsnew(i,j)*z1d2dt*rho_w  ! available water
            zdwsndt(i,j)  = zdwsndt (i,j) - zdwsnm(i,j)
            zdwgme        = zdwsnm(i,j)*zrock(i,j)             ! contribution to w_so
            zro           = (1._ireals - zrock(i,j))*zdwsnm(i,j)      ! surface runoff
            zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,j,1) -  &
                            zfcap(i,j))/MAX(zporv(i,j)-zfcap(i,j), zepsi)))
            zdwgdt(i,j,1) = zdwgdt(i,j,1) + zdwgme*(1._ireals - zredfu)
            zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                   ! for this fraction

            ! zro-, zdw_so_dt-correction in case of pore volume overshooting
            zw_ovpv = MAX(0._ireals, zw_fr(i,j,1)* zdzhs(1) * zrhwddt +  &
                       zdwgdt(i,j,1) - zporv(i,j) * zdzhs(1) * zrhwddt)
            zro = zro + zw_ovpv
            zdwgdt(i,j,1)= zdwgdt(i,j,1) - zw_ovpv


            IF (zfr_melt > 0.9999_ireals) zdwsndt(i,j)= -zwsnow(i,j)*zrhwddt
!subs            runoff_s (i,j)= runoff_s(i,j) + zro*zroffdt
            runoff_s (i,j,ns)= runoff_s(i,j,ns) + zro*zroffdt
!           IF(runoff_s(i,j).LT.0.) THEN
!             WRITE (yerror,'(A,I5,A,I5,A,F10.5)')                        &
!              'Negative runoff_s for i= ',i,' and j = ',j,': ',runoff_s(i,j)
!             ierror = 1
!           ENDIF
!           IF(runoff_g(i,j).LT.0.) THEN
!             WRITE (yerror,'(A,I5,A,I5,A,F10.5)')                        &
!              'Negative runoff_g for i= ',i,' and j = ',j,': ',runoff_g(i,j)
!             ierror = 1
!           ENDIF
          END IF   ! snow melting
        END IF     ! snow and/or soil temperatures
      END IF       ! points with snow cover only
    END IF         ! land-points only
  END DO
END DO

ELSE       ! new snow scheme
DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
    IF (llandmask(i,j,ns)) THEN          ! land-points only

      zwsnn  (i,j)  = zwsnow(i,j) + zdtdrhw*zdwsndt(i,j)
      zwsnew (i,j)  = zwsnn(i,j)
      ze_avail      = 0.0_ireals
      ze_total      = 0.0_ireals
      zfr_melt      = 0.0_ireals
      zdwsnm(i,j)   = 0.0_ireals

      IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only

        ze_out      = 0.0_ireals
        zqbase(i,j) = 0.0_ireals
        zcounter    = 0.0_ireals

        IF(zextinct(i,j,1).eq.0.0_ireals) then
          ze_rad = 0.0_ireals
        ELSE
!subs          ze_rad = zf_snow(i,j) * sobs(i,j) !* exp(-zextinct(1)*zdzm_snow(1))
          ze_rad = zf_snow(i,j) * sobs(i,j,ns) !* exp(-zextinct(1)*zdzm_snow(1))
        END IF

        ztsnownew_mult(i,j,0) = ztsnown_mult(i,j,0)
        DO ksn = 1,ke_snow

          zrefr(i,j) = 0.0_ireals
          zmelt(i,j) = 0.0_ireals
          ztsnownew_mult(i,j,ksn) = ztsnown_mult(i,j,ksn)

!em          zrho_dry_old(i,j) = MAX(wtot_snow(i,j,ksn,nx) - wliq_snow(i,j,ksn,nx), 0._ireals)*   &
!em                         rho_w/(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx))
!subs          IF(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx).GT.zepsi .OR. &
!subs            wtot_snow(i,j,ksn,nx) - wliq_snow(i,j,ksn,nx).GT.zepsi) THEN
!subs            zrho_dry_old(i,j) = MAX(wtot_snow(i,j,ksn,nx) - wliq_snow(i,j,ksn,nx), zepsi)*   &
!subs                                rho_w/(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx))
          IF(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx,ns).GT.zepsi .OR. &
            wtot_snow(i,j,ksn,nx,ns) - wliq_snow(i,j,ksn,nx,ns).GT.zepsi) THEN
            zrho_dry_old(i,j) = MAX(wtot_snow(i,j,ksn,nx,ns) - wliq_snow(i,j,ksn,nx,ns), zepsi)*  &
                                rho_w/(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx,ns))
          ELSE
            zrho_dry_old(i,j) = rho_w
          END IF


!subs          ztsnownew_mult(i,j,ksn) = (ztsnown_mult(i,j,ksn)*wtot_snow(i,j,ksn,nx) + t0_melt*zqbase(i,j)*zdt)/  &
!subs            (zqbase(i,j)*zdt+wtot_snow(i,j,ksn,nx)) !+ zqbase(i,j)*zdt*lh_f/(chc_i*wtot_snow(ksn,nx))
          ztsnownew_mult(i,j,ksn) = (ztsnown_mult(i,j,ksn)*wtot_snow(i,j,ksn,nx,ns) &
            &                       + t0_melt*zqbase(i,j)*zdt)/(zqbase(i,j)*zdt     &
            &                       + wtot_snow(i,j,ksn,nx,ns))

          IF(zextinct(i,j,ksn).eq.0.0_ireals) THEN
            ze_in = ze_out
          ELSE
            IF(ksn.eq.ke_snow) THEN
              ze_in = ze_out + (ze_rad - zcounter)
            ELSEIF(ksn.eq.1) then
              ze_in = ze_rad * (EXP (-zextinct(i,j,1)*zdzm_snow(i,j,1)) &
                &     - EXP (-zextinct(i,j,1)*zdzh_snow(i,j,1)))
              zcounter = ze_rad * (1._ireals - EXP (-zextinct(i,j,1)*zdzh_snow(i,j,1)))
            ELSE
              ze_in = ze_out + (ze_rad - zcounter) -  &
                (ze_rad - zcounter) * EXP (-zextinct(i,j,ksn)*zdzh_snow(i,j,ksn))
              zcounter = ze_rad - (ze_rad - zcounter) * EXP (-zextinct(i,j,ksn)*zdzh_snow(i,j,ksn))
            END IF
          END IF

!subs          ztsnownew_mult(i,j,ksn) = ztsnownew_mult(i,j,ksn) + ze_in*zdt/(chc_i*wtot_snow(i,j,ksn,nx))/rho_w
!subs          wtot_snow(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx) + zqbase(i,j)*zdt
!subs          wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) + zqbase(i,j)*zdt
          ztsnownew_mult(i,j,ksn) = ztsnownew_mult(i,j,ksn) &
            &                       + ze_in*zdt/(chc_i*wtot_snow(i,j,ksn,nx,ns))/rho_w
          wtot_snow(i,j,ksn,nx,ns) = wtot_snow(i,j,ksn,nx,ns) + zqbase(i,j)*zdt
          wliq_snow(i,j,ksn,nx,ns) = wliq_snow(i,j,ksn,nx,ns) + zqbase(i,j)*zdt

          zdzh_old = zdzh_snow(i,j,ksn)
          zdzh_snow(i,j,ksn) = zdzh_snow(i,j,ksn) + zqbase(i,j)*zdt

!subs          rho_snow_mult(i,j,ksn,nx) = MAX(wtot_snow(i,j,ksn,nx)*rho_w/zdzh_snow(i,j,ksn),0.0_ireals)
          rho_snow_mult(i,j,ksn,nx,ns) = MAX(wtot_snow(i,j,ksn,nx,ns)*rho_w/zdzh_snow(i,j,ksn), &
            &                                0.0_ireals)

          IF(ztsnownew_mult(i,j,ksn) .GT. t0_melt) THEN

!subs            IF(wtot_snow(i,j,ksn,nx) .LE. wliq_snow(i,j,ksn,nx)) THEN
!subs              ze_out = chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt*rho_w
            IF(wtot_snow(i,j,ksn,nx,ns) .LE. wliq_snow(i,j,ksn,nx,ns)) THEN
              ze_out = chc_i*wtot_snow(i,j,ksn,nx,ns)*(ztsnownew_mult(i,j,ksn) - t0_melt) &
                &      *z1d2dt*rho_w
              zmelt(i,j) = 0.0_ireals
!subs            ELSEIF(chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)/lh_f .LE.  &
!subs              wtot_snow(i,j,ksn,nx)-wliq_snow(i,j,ksn,nx)) THEN
!subs              zmelt(i,j) = chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt/lh_f
            ELSEIF(chc_i*wtot_snow(i,j,ksn,nx,ns)*(ztsnownew_mult(i,j,ksn) - t0_melt)/lh_f .LE.  &
              wtot_snow(i,j,ksn,nx,ns)-wliq_snow(i,j,ksn,nx,ns)) THEN
              zmelt(i,j) = chc_i*wtot_snow(i,j,ksn,nx,ns)*(ztsnownew_mult(i,j,ksn) - t0_melt) &
                &          *z1d2dt/lh_f
              ze_out = 0.0_ireals
!subs              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) + zmelt(i,j)*zdt
              wliq_snow(i,j,ksn,nx,ns) = wliq_snow(i,j,ksn,nx,ns) + zmelt(i,j)*zdt
            ELSE
!subs              zmelt(i,j) = (wtot_snow(i,j,ksn,nx)-wliq_snow(i,j,ksn,nx))*z1d2dt
!subs              ze_out = chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt*rho_w - zmelt(i,j)*lh_f*rho_w
!subs              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) + zmelt(i,j)*zdt
              zmelt(i,j) = (wtot_snow(i,j,ksn,nx,ns)-wliq_snow(i,j,ksn,nx,ns))*z1d2dt
              ze_out = chc_i*wtot_snow(i,j,ksn,nx,ns)*(ztsnownew_mult(i,j,ksn) - t0_melt) &
                &      *z1d2dt*rho_w - zmelt(i,j)*lh_f*rho_w
              wliq_snow(i,j,ksn,nx,ns) = wliq_snow(i,j,ksn,nx,ns) + zmelt(i,j)*zdt
            END IF
            ztsnownew_mult(i,j,ksn) = t0_melt

          ELSE
!       T<0
!subs            IF(wliq_snow(i,j,ksn,nx) .GT. -chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)/lh_f) THEN
!subs              zrefr(i,j)            = -chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt/lh_f
            IF(wliq_snow(i,j,ksn,nx,ns) .GT. -chc_i*wtot_snow(i,j,ksn,nx,ns) &
              & *(ztsnownew_mult(i,j,ksn) - t0_melt)/lh_f) THEN
              zrefr(i,j)            = -chc_i*wtot_snow(i,j,ksn,nx,ns)*(ztsnownew_mult(i,j,ksn) &
                &                     - t0_melt)*z1d2dt/lh_f
              ztsnownew_mult(i,j,ksn)   = t0_melt
!subs              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) - zrefr(i,j)*zdt
              wliq_snow(i,j,ksn,nx,ns) = wliq_snow(i,j,ksn,nx,ns) - zrefr(i,j)*zdt
            ELSE
!subs              zrefr(i,j)            = wliq_snow(i,j,ksn,nx)*z1d2dt
!subs              wliq_snow(i,j,ksn,nx) = 0.0_ireals
!subs              ztsnownew_mult(i,j,ksn)   = ztsnownew_mult(i,j,ksn) + zrefr(i,j)*zdt*lh_f/(chc_i*wtot_snow(i,j,ksn,nx))
              zrefr(i,j)            = wliq_snow(i,j,ksn,nx,ns)*z1d2dt
              wliq_snow(i,j,ksn,nx,ns) = 0.0_ireals
              ztsnownew_mult(i,j,ksn)   = ztsnownew_mult(i,j,ksn) + zrefr(i,j)*zdt*lh_f &
                &                         /(chc_i*wtot_snow(i,j,ksn,nx,ns))
            END IF
            ze_out = 0.0_ireals

          END IF

          zdtsnowdt_mult(i,j,ksn) = zdtsnowdt_mult(i,j,ksn) + &
                                    (ztsnownew_mult(i,j,ksn) - ztsnown_mult(i,j,ksn))*z1d2dt
!subs          IF(wtot_snow(i,j,ksn,nx) .LE. wliq_snow(i,j,ksn,nx)) THEN
!subs            zqbase(i,j)           = wliq_snow(i,j,ksn,nx)*z1d2dt
!subs            wliq_snow(i,j,ksn,nx) = 0.0_ireals
!subs            wtot_snow(i,j,ksn,nx) = 0.0_ireals
          IF(wtot_snow(i,j,ksn,nx,ns) .LE. wliq_snow(i,j,ksn,nx,ns)) THEN
            zqbase(i,j)           = wliq_snow(i,j,ksn,nx,ns)*z1d2dt
            wliq_snow(i,j,ksn,nx,ns) = 0.0_ireals
            wtot_snow(i,j,ksn,nx,ns) = 0.0_ireals
            zdzh_snow(i,j,ksn)    = 0.0_ireals
!subs            rho_snow_mult(i,j,ksn,nx)  = 0.0_ireals
            rho_snow_mult(i,j,ksn,nx,ns)  = 0.0_ireals
          ELSE
            IF(zrefr(i,j) .GT. 0.0_ireals .OR. zmelt(i,j) .GT. 0.0_ireals) THEN
              zadd_dz = 0.0_ireals
              zadd_dz = MAX(zrefr(i,j),0._ireals)*(-1.0_ireals + 1.0_ireals/rho_i*rho_w)*zdt
              zadd_dz = MAX(zmelt(i,j),0._ireals)*(-1.0_ireals/zrho_dry_old(i,j)*rho_w &
                &       + 1.0_ireals)*zdt
              zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn) + zadd_dz
!subs              rho_snow_mult(i,j,ksn,nx) = MAX(wtot_snow(i,j,ksn,nx)*rho_w/zdzh_snow(i,j,ksn),0.0_ireals)
!subs              IF(wtot_snow(i,j,ksn,nx) .LE. 0.0_ireals) zdzh_snow(i,j,ksn) = 0.0_ireals
!subs              IF(rho_snow_mult(i,j,ksn,nx) .GT. rho_w) THEN
!subs                zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn)*rho_snow_mult(i,j,ksn,nx)/rho_w
!subs                rho_snow_mult(i,j,ksn,nx) = rho_w
              rho_snow_mult(i,j,ksn,nx,ns) = MAX(wtot_snow(i,j,ksn,nx,ns)*rho_w &
                &                            /zdzh_snow(i,j,ksn),0.0_ireals)
              IF(wtot_snow(i,j,ksn,nx,ns) .LE. 0.0_ireals) zdzh_snow(i,j,ksn) = 0.0_ireals
              IF(rho_snow_mult(i,j,ksn,nx,ns) .GT. rho_w) THEN
                zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn)*rho_snow_mult(i,j,ksn,nx,ns)/rho_w
                rho_snow_mult(i,j,ksn,nx,ns) = rho_w
              END IF
            END IF

!subs            zsn_porosity = 1._ireals - (rho_snow_mult(i,j,ksn,nx)/rho_w -  &
!subs                           wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn))/rho_i*rho_w - &
!subs                           wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn)
            zsn_porosity = 1._ireals - (rho_snow_mult(i,j,ksn,nx,ns)/rho_w -  &
                           wliq_snow(i,j,ksn,nx,ns)/zdzh_snow(i,j,ksn))/rho_i*rho_w - &
                           wliq_snow(i,j,ksn,nx,ns)/zdzh_snow(i,j,ksn)
            zsn_porosity = MAX(zsn_porosity,cwhc + 0.1_ireals)
            zp1 = zsn_porosity - cwhc

!subs            IF (wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn) .GT. cwhc) THEN
!subs              zfukt             = (wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn) - cwhc)/zp1
            IF (wliq_snow(i,j,ksn,nx,ns)/zdzh_snow(i,j,ksn) .GT. cwhc) THEN
              zfukt             = (wliq_snow(i,j,ksn,nx,ns)/zdzh_snow(i,j,ksn) - cwhc)/zp1
              zq0               = chcond * zfukt**3.0_ireals
!subs              zqbase(i,j)       = MIN(zq0*zdt,wliq_snow(i,j,ksn,nx))
!subs              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) - zqbase(i,j)
!subs              wtot_snow(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx) - zqbase(i,j)
              zqbase(i,j)       = MIN(zq0*zdt,wliq_snow(i,j,ksn,nx,ns))
              wliq_snow(i,j,ksn,nx,ns) = wliq_snow(i,j,ksn,nx,ns) - zqbase(i,j)
              wtot_snow(i,j,ksn,nx,ns) = wtot_snow(i,j,ksn,nx,ns) - zqbase(i,j)

              zdzh_old = zdzh_snow(i,j,ksn)
              zdzh_snow(i,j,ksn) = zdzh_snow(i,j,ksn) - zqbase(i,j)
              zqbase(i,j)        = zqbase(i,j)*z1d2dt

              IF(zdzh_snow(i,j,ksn) .LT. zepsi*0.01_ireals) THEN
!subs                wliq_snow(i,j,ksn,nx) = 0.0_ireals
!subs                wtot_snow(i,j,ksn,nx) = 0.0_ireals
                wliq_snow(i,j,ksn,nx,ns) = 0.0_ireals
                wtot_snow(i,j,ksn,nx,ns) = 0.0_ireals
                zdzh_snow(i,j,ksn)    = 0.0_ireals
!subs                rho_snow_mult(i,j,ksn,nx)  = 0.0_ireals
                rho_snow_mult(i,j,ksn,nx,ns)  = 0.0_ireals
              ELSE
!subs                rho_snow_mult(i,j,ksn,nx) = MAX(wtot_snow(i,j,ksn,nx)*rho_w/zdzh_snow(i,j,ksn),0.0_ireals)
!subs                IF(wtot_snow(i,j,ksn,nx) .LE. 0.0_ireals) zdzh_snow(i,j,ksn) = 0.0_ireals
!subs                IF(rho_snow_mult(i,j,ksn,nx) .GT. rho_w) THEN
!subs                  zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn)*rho_snow_mult(i,j,ksn,nx)/rho_w
!subs                  rho_snow_mult(i,j,ksn,nx) = rho_w
                rho_snow_mult(i,j,ksn,nx,ns) = MAX(wtot_snow(i,j,ksn,nx,ns)*rho_w &
                  &                            /zdzh_snow(i,j,ksn),0.0_ireals)
                IF(wtot_snow(i,j,ksn,nx,ns) .LE. 0.0_ireals) zdzh_snow(i,j,ksn) = 0.0_ireals
                IF(rho_snow_mult(i,j,ksn,nx,ns) .GT. rho_w) THEN
                  zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn)*rho_snow_mult(i,j,ksn,nx,ns)/rho_w
                  rho_snow_mult(i,j,ksn,nx,ns) = rho_w
                END IF
              END IF
            ELSE
              zqbase(i,j) = 0.0_ireals
            END IF
          END IF

        END DO        ! snow layers
        zdwsnm(i,j) = zqbase(i,j)*rho_w       ! ksn == ke_snow
      END IF       ! points with snow cover only
    END IF         ! land-points only
  END DO
END DO

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN          ! land-points only
    IF (llandmask(i,j,ns)) THEN          ! land-points only
      IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
        zdwsndt(i,j)  = zdwsndt(i,j) - zdwsnm(i,j)

          ! melted snow is allowed to penetrate the soil (up to field
          ! capacity), if the soil type is neither ice nor rock (zrock = 0);
          ! else it contributes to surface run-off;
          ! fractional water content of the first soil layer determines
          ! a reduction factor which controls additional run-off

        zdwgme        = zdwsnm(i,j)*zrock(i,j)             ! contribution to w_so
        zro           = (1._ireals - zrock(i,j))*zdwsnm(i,j)      ! surface runoff
        zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,j,1) -  &
                        zfcap(i,j))/MAX(zporv(i,j)-zfcap(i,j), zepsi)))
        zdwgdt(i,j,1) = zdwgdt(i,j,1) + zdwgme*(1._ireals - zredfu)
        zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                   ! for this fraction

        ! zro-, zdw_so_dt-correction in case of pore volume overshooting
        zw_ovpv = MAX(0._ireals, zw_fr(i,j,1)* zdzhs(1) * zrhwddt +  &
                  zdwgdt(i,j,1) - zporv(i,j) * zdzhs(1) * zrhwddt)
        zro = zro + zw_ovpv
        zdwgdt(i,j,1)= zdwgdt(i,j,1) - zw_ovpv

!subs        runoff_s(i,j) = runoff_s(i,j) + zro*zroffdt
        runoff_s(i,j,ns) = runoff_s(i,j,ns) + zro*zroffdt
      END IF       ! points with snow cover only
    END IF         ! land-points only
  END DO
END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN          ! land-points only
      IF (llandmask(i,j,ns)) THEN          ! land-points only
        IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
! snow densification due to gravity and metamorphism 
          zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (zth_low(i,j)-csnow_tmin) &
                       /(t0_melt           -csnow_tmin)
          zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))
          DO ksn = 2, ke_snow
!subs            IF (rho_snow_mult(i,j,ksn,nx) .LT. 6*zrho_snowf &
!subs              & .AND. rho_snow_mult(i,j,ksn,nx) .NE. 0.0_ireals) THEN
!subs              zdens_old = rho_snow_mult(i,j,ksn,nx)
            IF (rho_snow_mult(i,j,ksn,nx,ns) .LT. 6*zrho_snowf &
              & .AND. rho_snow_mult(i,j,ksn,nx,ns) .NE. 0.0_ireals) THEN
              zdens_old = rho_snow_mult(i,j,ksn,nx,ns)
              zeta =         &! compactive viscosity of snow
!subs                ca2*EXP(19.3_ireals*rho_snow_mult(i,j,ksn,nx)/rho_i)* &
                ca2*EXP(19.3_ireals*rho_snow_mult(i,j,ksn,nx,ns)/rho_i)* &
                EXP(67300._ireals/8.31_ireals/ztsnownew_mult(i,j,ksn))
              zp = 0.0_ireals                         ! gravity, Pa
              DO k = ksn,1,-1
!subs                zp = zp + rho_snow_mult(i,j,k,nx)*g*zdzh_snow(i,j,ksn)
                zp = zp + rho_snow_mult(i,j,k,nx,ns)*g*zdzh_snow(i,j,ksn)
              END DO
!subs              rho_snow_mult(i,j,ksn,nx) = rho_snow_mult(i,j,ksn,nx) + zdt*rho_snow_mult(i,j,ksn,nx)*(csigma+zp)/zeta
!subs              rho_snow_mult(i,j,ksn,nx) = MIN(rho_snow_mult(i,j,ksn,nx),rho_i)
!subs              zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn) * zdens_old/rho_snow_mult(i,j,ksn,nx)
              rho_snow_mult(i,j,ksn,nx,ns) = rho_snow_mult(i,j,ksn,nx,ns) &
                &                          + zdt*rho_snow_mult(i,j,ksn,nx,ns)*(csigma+zp)/zeta
              rho_snow_mult(i,j,ksn,nx,ns) = MIN(rho_snow_mult(i,j,ksn,nx,ns),rho_i)
              zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn) * zdens_old/rho_snow_mult(i,j,ksn,nx,ns)
            END IF
          END DO

          IF(ztsnownew_mult(i,j,0) .GT. t0_melt) THEN
            ztsnownew_mult(i,j,0) = t0_melt
            zdtsnowdt_mult(i,j,0) = zdtsnowdt_mult(i,j,0) +     &
                                    (ztsnownew_mult(i,j,0) - ztsnown_mult(i,j,0))*z1d2dt
          END IF
        END IF       ! points with snow cover only
      END IF         ! land-points only
    END DO
  END DO
END IF

  !>JH WRITE(*,*) 'TERRA-DIAG END II.8 +++++++++> '

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

IF (lmulti_snow) THEN
  ! First for ksn == 0
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN  ! for landpoints only
!subs        t_snow_mult  (i,j,0,nnew) = t_snow_mult(i,j,0,nx) + zdt*zdtsnowdt_mult(i,j,0)
      IF (llandmask(i,j,ns)) THEN  ! for landpoints only
        t_snow_mult  (i,j,0,nnew,ns) = t_snow_mult(i,j,0,nx,ns) + zdt*zdtsnowdt_mult(i,j,0)
      ENDIF
    ENDDO
  ENDDO

  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN  ! for landpoints only
!subs          t_snow_mult  (i,j,ksn,nnew) = t_snow_mult(i,j,ksn,nx) + zdt*zdtsnowdt_mult(i,j,ksn)
!subs          dzh_snow     (i,j,ksn,nnew) = zdzh_snow(i,j,ksn)
!subs          wtot_snow    (i,j,ksn,nnew) = wtot_snow(i,j,ksn,nx)
!subs          rho_snow_mult(i,j,ksn,nnew) = rho_snow_mult(i,j,ksn,nx)
!subs          wliq_snow    (i,j,ksn,nnew) = wliq_snow(i,j,ksn,nx)
        IF (llandmask(i,j,ns)) THEN  ! for landpoints only
          t_snow_mult  (i,j,ksn,nnew,ns) = t_snow_mult(i,j,ksn,nx,ns) + zdt*zdtsnowdt_mult(i,j,ksn)
          dzh_snow     (i,j,ksn,nnew,ns) = zdzh_snow(i,j,ksn)
          wtot_snow    (i,j,ksn,nnew,ns) = wtot_snow(i,j,ksn,nx,ns)
          rho_snow_mult(i,j,ksn,nnew,ns) = rho_snow_mult(i,j,ksn,nx,ns)
          wliq_snow    (i,j,ksn,nnew,ns) = wliq_snow(i,j,ksn,nx,ns)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN  ! for landpoints only
!subs        t_snow(i,j,nnew)  = t_snow(i,j,nx) + zdt*zdtsnowdt(i,j)
      IF (llandmask(i,j,ns)) THEN  ! for landpoints only
        t_snow(i,j,nnew,ns)  = t_snow(i,j,nx,ns) + zdt*zdtsnowdt(i,j)
      ENDIF
    ENDDO
  ENDDO
ENDIF

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN  ! for landpoints only
    IF (llandmask(i,j,ns)) THEN  ! for landpoints only
      ! t_snow is computed above
      ! t_snow(i,j,nnew)  = t_snow(i,j,nx) + zdt*zdtsnowdt(i,j)
!subs      t_so(i,j,1,nnew)  = t_so(i,j,1,nx) + zdt*zdtsdt   (i,j)
      t_so(i,j,1,nnew,ns)  = t_so(i,j,1,nx,ns) + zdt*zdtsdt   (i,j)
!IF(i.eq.340 .and. j.eq.186) &
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
! print *,ns,'after snow',t_so(i,j,1,nnew,ns),zdt*zdtsdt   (i,j)
      ! Next line has to be changed, if the soil surface temperature
      ! t_so(i,j,0,nnew) predicted by the heat conduction equation is used
!subs      t_s   (i,j,nnew)    = t_so(i,j,1,nnew)
!subs      t_so  (i,j,0,nnew)  = t_so(i,j,1,nnew)
!subs      w_snow(i,j,nnew)  = w_snow(i,j,nx) + zdt*zdwsndt  (i,j)/rho_w
!subs      w_i   (i,j,nnew)  = w_i   (i,j,nx) + zdt*zdwidt   (i,j)/rho_w
      t_s   (i,j,nnew,ns)    = t_so(i,j,1,nnew,ns)
      t_so  (i,j,0,nnew,ns)  = t_so(i,j,1,nnew,ns)
      w_snow(i,j,nnew,ns)  = w_snow(i,j,nx,ns) + zdt*zdwsndt  (i,j)/rho_w
      w_i   (i,j,nnew,ns)  = w_i   (i,j,nx,ns) + zdt*zdwidt   (i,j)/rho_w

!subs      IF (w_snow(i,j,nnew) <= zepsi) THEN
!subs        w_snow(i,j,nnew) = 0.0_ireals
!subs        t_snow(i,j,nnew) = t_so(i,j,0,nnew)
      IF (w_snow(i,j,nnew,ns) <= zepsi) THEN
        w_snow(i,j,nnew,ns) = 0.0_ireals
        t_snow(i,j,nnew,ns) = t_so(i,j,0,nnew,ns)
      ENDIF
!subs      IF (w_i(i,j,nnew) <= zepsi) w_i(i,j,nnew) = 0.0_ireals
      IF (w_i(i,j,nnew,ns) <= zepsi) w_i(i,j,nnew,ns) = 0.0_ireals
    END IF          ! land-points only
  END DO
END DO

! Eliminate snow for multi-layer snow model, if w_snow = 0
IF (lmulti_snow) THEN
  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
!subs        IF (llandmask(i,j)) THEN  ! for landpoints only
!subs          IF (w_snow(i,j,nnew) == 0.0_ireals) THEN
!subs            t_snow_mult(i,j,ksn,nnew) = t_so(i,j,0,nnew)
!subs            wliq_snow(i,j,ksn,nnew) = 0.0_ireals
!subs            wtot_snow(i,j,ksn,nnew) = 0.0_ireals
!subs            rho_snow_mult (i,j,ksn,nnew) = 0.0_ireals
!subs            dzh_snow (i,j,ksn,nnew) = 0.0_ireals
        IF (llandmask(i,j,ns)) THEN  ! for landpoints only
          IF (w_snow(i,j,nnew,ns) == 0.0_ireals) THEN
            t_snow_mult(i,j,ksn,nnew,ns) = t_so(i,j,0,nnew,ns)
            wliq_snow(i,j,ksn,nnew,ns) = 0.0_ireals
            wtot_snow(i,j,ksn,nnew,ns) = 0.0_ireals
            rho_snow_mult (i,j,ksn,nnew,ns) = 0.0_ireals
            dzh_snow (i,j,ksn,nnew,ns) = 0.0_ireals
          ENDIF
        END IF          ! land-points only
      END DO
    END DO
  END DO
ENDIF

IF(.NOT. lmulti_snow) THEN

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN  ! for landpoints only
    IF (llandmask(i,j,ns)) THEN  ! for landpoints only
!BR 7/2005 Update snow density
!
!     a) aging of existing snow
!
!    temperature dependence of relaxation/ageing constant
!subs     ztau_snow = crhosmint+(crhosmaxt-crhosmint)*(t_snow(i,j,nnew)-csnow_tmin) &
     ztau_snow = crhosmint+(crhosmaxt-crhosmint)*(t_snow(i,j,nnew,ns)-csnow_tmin) &
                                                /(t0_melt         -csnow_tmin)
     ztau_snow = MAX(crhosmint,MIN(crhosmaxt,ztau_snow))
!subs     zrho_snowe= crhosmax_ml+(rho_snow(i,j,nx)-crhosmax_ml)* &
     zrho_snowe= crhosmax_ml+(rho_snow(i,j,nx,ns)-crhosmax_ml)* &
                           EXP(-ztau_snow*zdt/86400._ireals)
!     b) density of fresh snow
!
     zrho_snowf= crhosminf+(crhosmaxf-crhosminf)* (zth_low(i,j)-csnow_tmin) &
                                                /(t0_melt      -csnow_tmin)
     zrho_snowf= MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))
!
!     c) new snow density is weighted average of existing and new snow
!
     IF ( itype_gscp == 4 ) THEN
!subs       znorm=MAX(w_snow(i,j,nx)+(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j))      &
       znorm=MAX(w_snow(i,j,nx,ns)+(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j))      &
                 *zdtdrhw,zepsi)
!subs       rho_snow(i,j,nnew)  = ( zrho_snowe*w_snow(i,j,nx) + &
       rho_snow(i,j,nnew,ns)  = ( zrho_snowe*w_snow(i,j,nx,ns) + &
                          zrho_snowf*(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j)) &
                             *zdtdrhw )    /znorm
     ELSE
!subs       znorm=MAX(w_snow(i,j,nx)+(prs_gsp(i,j)+prs_con(i,j) )      &
       znorm=MAX(w_snow(i,j,nx,ns)+(prs_gsp(i,j)+prs_con(i,j) )      &
                 *zdtdrhw,zepsi)
!subs       rho_snow(i,j,nnew)  = ( zrho_snowe*w_snow(i,j,nx) + &
       rho_snow(i,j,nnew,ns)  = ( zrho_snowe*w_snow(i,j,nx,ns) + &
                          zrho_snowf*(prs_gsp(i,j)+prs_con(i,j) ) &
                             *zdtdrhw )    /znorm
     ENDIF
!subs     rho_snow(i,j,nnew) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow(i,j,nnew)))
     rho_snow(i,j,nnew,ns) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow(i,j,nnew,ns)))
    END IF          ! land-points only
  END DO
END DO

ELSE   ! new snow scheme

DO   j = jstarts, jends
  DO i = istarts, iends
!subs    IF (llandmask(i,j)) THEN  ! for landpoints only
!subs      IF(w_snow(i,j,nnew) .GT. zepsi) THEN
!subs        h_snow(i,j,nnew) = 0.0_ireals
    IF (llandmask(i,j,ns)) THEN  ! for landpoints only
      IF(w_snow(i,j,nnew,ns) .GT. zepsi) THEN
        h_snow(i,j,nnew,ns) = 0.0_ireals
        DO ksn = 1,ke_snow
!subs          h_snow(i,j,nnew) = h_snow(i,j,nnew) + zdzh_snow(i,j,ksn)
          h_snow(i,j,nnew,ns) = h_snow(i,j,nnew,ns) + zdzh_snow(i,j,ksn)
        END DO
      ELSE
!subs        h_snow(i,j,nnew) = 0._ireals
        h_snow(i,j,nnew,ns) = 0._ireals
      END IF

      zflag = 0._ireals
      DO ksn = 1,ke_snow
!subs        IF(wtot_snow(i,j,ksn,nnew) .LT. zepsi) zflag = 1._ireals
        IF(wtot_snow(i,j,ksn,nnew,ns) .LT. zepsi) zflag = 1._ireals
      END DO
!subs      IF(w_snow(i,j,nx) .GT. zepsi .AND. w_snow(i,j,nnew) .GT. zepsi) THEN
      IF(w_snow(i,j,nx,ns) .GT. zepsi .AND. w_snow(i,j,nnew,ns) .GT. zepsi) THEN
!subs        CALL normalize(h_snow(i,j,nnew),zdzm_snow(i,j,:),dzh_snow(i,j,:,nnew),                  &
!subs                       zhh_snow(i,j,:),wtot_snow(i,j,:,nnew),               &
!subs                       t_snow_mult(i,j,:,nnew),wliq_snow(i,j,:,nnew),rho_snow_mult(i,j,:,nnew),i,j)
        CALL normalize(h_snow(i,j,nnew,ns),zdzm_snow(i,j,:),dzh_snow(i,j,:,nnew,ns),         &
          &            zhh_snow(i,j,:),wtot_snow(i,j,:,nnew,ns), t_snow_mult(i,j,:,nnew,ns), &
          &            wliq_snow(i,j,:,nnew,ns),rho_snow_mult(i,j,:,nnew,ns),ke_snow)
      END IF
    END IF          ! land-points only
  END DO
END DO

ENDIF ! lmulti_snow

DO kso = 1,ke_soil
  DO   j = jstarts, jends
    DO i = istarts, iends
!subs      IF (llandmask(i,j)) THEN  ! for landpoints only
!subs        w_so(i,j,kso,nnew) = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
      IF (llandmask(i,j,ns)) THEN  ! for landpoints only
        w_so(i,j,kso,nnew,ns) = w_so(i,j,kso,nx,ns) + zdt*zdwgdt(i,j,kso)/rho_w
!        print*,'w_so(i,j,kso,nnew,ns) ',kso,nnew,ns,w_so(25,36,kso,nnew,ns)
      END IF  ! land-points only
    END DO
  END DO
END DO        ! soil layers

!subs amf
enddo
!amf ---

!subs em 
IF(itype_subs .EQ. 2) THEN    ! tiles
  DO ns = nsubs0+1, nsubs1, 2      ! 1 - mean, 2 - ocean, 3 - lake, 4 - no snow first
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j,ns)) THEN  ! for landpoints only  !check is not necessary

!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
! print *,'w_snow in the end:',w_snow(i,j,nnew,ns-1),w_snow(i,j,nnew,ns  )
!dev          w_snow(i,j,nnew,ns  ) = w_snow(i,j,nnew,ns) + &
!dev                                  w_snow(i,j,nnew,ns-1)*subsfrac(i,j,ns-1)/MAX(subsfrac(i,j,ns),zepsi)
          fact1 = subsfrac(i,j,ns-1)+subsfrac(i,j,ns)
          w_snow(i,j,nnew,ns  ) = (w_snow(i,j,nnew,ns)*subsfrac(i,j,ns) + &
                                  w_snow(i,j,nnew,ns-1)*subsfrac(i,j,ns-1))/MAX(fact1,zepsi)
          w_snow(i,j,nnew,ns-1) = 0._ireals
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
! print *,'w_snow in the end-2:',w_snow(i,j,nnew,ns-1),w_snow(i,j,nnew,ns  )

          zf_snow_old(i,j) = subsfrac(i,j,ns) / (subsfrac(i,j,ns-1) + subsfrac(i,j,ns))
!dev          zf_snow    (i,j) = MAX( 0.01_ireals, MIN(1.0_ireals,w_snow(i,j,nnew,ns)*zf_snow_old(i,j)/cf_snow) )* &
          zf_snow    (i,j) = MAX( 0.01_ireals, MIN(1.0_ireals,w_snow(i,j,nnew,ns)/cf_snow) )* &
                             zsf_heav(w_snow(i,j,nnew,ns) - zepsi)

!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
! print *,'zf_snow_old, zf_snow',zf_snow_old(i,j),zf_snow    (i,j)

!dev          w_snow(i,j,nnew,ns  ) = w_snow(i,j,nnew,ns) * zf_snow_old(i,j)/MAX(zf_snow(i,j),zepsi)
          w_snow(i,j,nnew,ns  ) = w_snow(i,j,nnew,ns)/MAX(zf_snow(i,j),zepsi)

!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!print *,'w_snow in the end-3:',w_snow(i,j,nnew,ns-1),w_snow(i,j,nnew,ns  )

          subsfrac(i,j,ns-1) = (1._ireals - zf_snow(i,j))*fact1
          subsfrac(i,j,ns  ) = zf_snow(i,j)              *fact1
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. my_cart_id.eq.isub) &
!print *,'subsfrac:',subsfrac(i,j,ns-1),subsfrac(i,j,ns  )
        END IF  ! land-points only
      END DO
    END DO
    DO kso = 1,ke_soil
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j,ns)) THEN  ! for landpoints only  !check is not necessary

!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. kso.eq.1 .and. my_cart_id.eq.isub) &
!print *,i,j,kso,'temp before',t_so    (i,j,kso,nnew,ns-1),t_so    (i,j,kso,nnew,ns)
            fact1 = (MIN(1._ireals - zf_snow_old(i,j), 1._ireals - zf_snow(i,j))) &
              &     /MAX((1._ireals - zf_snow(i,j)),zepsi) 
            fact2 = (MAX(zf_snow(i,j) - zf_snow_old(i,j), 0._ireals))/MAX(zf_snow(i,j), zepsi) 

!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. kso.eq.1 .and. my_cart_id.eq.isub) &
! print *,'fact:',fact1,fact2
            tmp1 = t_so    (i,j,kso,nnew,ns-1)
            tmp2 = w_so    (i,j,kso,nnew,ns-1)
            tmp3 = w_so_ice(i,j,kso,nnew,ns-1)

            t_so    (i,j,kso,nnew,ns-1) = t_so(i,j,kso,nnew,ns-1)*fact1 &
              &                         + t_so(i,j,kso,nnew,ns)*(1._ireals - fact1)
            w_so    (i,j,kso,nnew,ns-1) = w_so(i,j,kso,nnew,ns-1)*fact1 &
              &                         + w_so(i,j,kso,nnew,ns)*(1._ireals - fact1)
            w_so_ice(i,j,kso,nnew,ns-1) = w_so_ice(i,j,kso,nnew,ns-1)*fact1 &
              &                         + w_so_ice(i,j,kso,nnew,ns)*(1._ireals - fact1)

            t_so    (i,j,kso,nnew,ns) = tmp1*fact2 + t_so(i,j,kso,nnew,ns)*(1._ireals - fact2)
            w_so    (i,j,kso,nnew,ns) = tmp2*fact2 + w_so(i,j,kso,nnew,ns)*(1._ireals - fact2)
            w_so_ice(i,j,kso,nnew,ns) = tmp3*fact2 + w_so_ice(i,j,kso,nnew,ns)*(1._ireals - fact2)

!print *,i,j,'t_snow in TERRA, end',t_snow(i,j,nnew,ns)
!print *,i,j,'w_snow in TERRA, end',w_snow(i,j,nnew,ns)
!print *,i,j,'t_s    in TERRA, end',t_s (i,j,nnew,ns)
!>JH IF(i.eq.i_loc .and. j.eq.j_loc .and. kso.eq.1 .and. my_cart_id.eq.isub) &
! print *,i,j,kso,'temp after',t_so    (i,j,kso,nnew,ns-1),t_so    (i,j,kso,nnew,ns)
          END IF  ! land-points only
        END DO
      END DO
    END DO        ! soil layers
  END DO
END IF

DO ns = nsubs0, nsubs1
!em>

! computation of the temperature at the boundary soil/snow-atmosphere
IF(lmulti_snow) THEN
!subs  CALL tgcom ( t_g(:,:,nnew), t_snow_mult(:,:,1,nnew), t_s(:,:,nnew), &
!subs               w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow,     &
!>JH  CALL tgcom ( t_g(:,:,nnew,ns), t_snow_mult(:,:,1,nnew,ns), t_s(:,:,nnew,ns), &
!               w_snow(:,:,nnew,ns), llandmask(:,:,ns), ie, je, cf_snow,     &
!<JH            istarts, iends, jstarts, jends )
ELSE
!subs  CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), t_s(:,:,nnew),        &
!subs               w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow,     &
!>JH  CALL tgcom ( t_g(:,:,nnew,ns), t_snow(:,:,nnew,ns), t_s(:,:,nnew,ns),        &
!               w_snow(:,:,nnew,ns), llandmask(:,:,ns), ie, je, cf_snow,     &
!<JH            istarts, iends, jstarts, jends )

ENDIF

!subs amf
!print *,ns,'in TERRA, t_g',w_snow(453,356,nnew,ns),t_snow(453,356,nnew,ns),t_s(453,356,nnew,ns),t_g(453,356,nnew,ns)
!print *,ns,'in TERRA, w_snow',MAXVAL(w_snow(:,:,nnew,ns)),MINVAL(w_snow(:,:,nnew,ns))
!print *,ns,'in TERRA, t_snow',MAXVAL(t_snow(:,:,nnew,ns)),MINVAL(t_snow(:,:,nnew,ns))
!print *,ns,'in TERRA, t_g',MAXVAL(t_g(:,:,nnew,ns)),MINVAL(t_g(:,:,nnew,ns))
enddo
!print *,'t_g    in TERRA, end',t_g (340,186,nnew,:)
!print *,'t_snow in TERRA, end',t_snow(340,186,nnew,:)
!print *,'w_snow in TERRA, end',w_snow(340,186,nnew,:)
!print *,'t_s    in TERRA, end',t_s (340,186,nnew,:)
!IF(ntstep.eq.2) STOP
!amf ---

#ifdef NECSX
  CALL collapse(.FALSE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif

  !>JH WRITE(*,*) 'TERRA-DIAG END II.9 END MODULE +++++++++> '

!------------------------------------------------------------------------------
! End of module procedure terra_multlay
!------------------------------------------------------------------------------



END SUBROUTINE terra_multlay

!==============================================================================

SUBROUTINE normalize(hsnow_new,dzm,dzh,zh,wtr,t,wl,rho,ke_snow)

!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------

INTEGER ::  ke_snow

REAL    (KIND=ireals   )                ::  &
  hsnow_new    ,  &
  dzm(ke_snow) ,  &
  dzh(ke_snow) ,  &
  zm(ke_snow)  ,  &
  zh(ke_snow)  ,  &
  wtr(ke_snow) ,  &
  wl(ke_snow)  ,  &
  rho(ke_snow) ,  &
  t(0:ke_snow)

! 
! Local scalars:
! -------------

INTEGER (KIND=iintegers) ::  &
!
!   Indices
!  
  ksn,    &       ! loop index for snow layers
  k               ! loop index for snow layers

REAL    (KIND=ireals   ) ::  &
!  
  z_old(ke_snow)     , &
  dz_old(ke_snow)    , &
  weight(ke_snow)    , & 
  t_new(0:ke_snow)     , &
  wtr_new(ke_snow)   , &
  wl_new(ke_snow)    , &
  rho_new(ke_snow)   , &
  sum_weight
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine normalize
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the
!------------------------------------------------------------------------------

  sum_weight = 0.0_ireals
  DO ksn = ke_snow,1,-1
    dz_old(ksn) = dzh(ksn)
    z_old(ksn) = -sum_weight - dzh(ksn)/2._ireals
    sum_weight = sum_weight + dzh(ksn)
  END DO

  DO ksn = 1,ke_snow
 !   zh(ksn) = -hsnow_new/REAL(ke_snow,ireals)*REAL((ke_snow-ksn),ireals))
     zh(ksn) =-hsnow_new/(REAL(ke_snow,ireals))*REAL(ke_snow-ksn,ireals)
  END DO

  zm(1) = (-hsnow_new + zh(1))/2._ireals
  DO ksn = 2,ke_snow
    zm(ksn) = (zh(ksn) + zh(ksn-1))/2._ireals
  END DO

  dzh(1) = zh(1) + hsnow_new      !layer thickness betw. half levels of uppermost snow layer
  dzm(1) = zm(1) + hsnow_new      !layer thickness between snow surface and main level of uppermost layer

  DO ksn = 2,ke_snow
    dzh(ksn) = zh(ksn) - zh(ksn-1) ! layer thickness betw. half levels
    dzm(ksn) = zm(ksn) - zm(ksn-1) ! layer thickness betw. main levels
  END DO

  DO ksn = 1, ke_snow
    IF(dz_old(ksn).ne.0._ireals.and.rho(ksn).ne.0._ireals) THEN
      wl(ksn) = wl(ksn)/dz_old(ksn)
      wtr(ksn) = wtr(ksn)/dz_old(ksn)
    END IF
  END DO

  DO ksn = ke_snow,1,-1
    t_new(ksn)   = 0.0_ireals
    rho_new(ksn) = 0.0_ireals
    wtr_new(ksn) = 0.0_ireals
    wl_new(ksn)  = 0.0_ireals

    DO k = ke_snow,1,-1
      weight(k) = 0.0_ireals

      weight(k) = MAX(MIN(z_old(k)+dz_old(k)/2._ireals,zm(ksn)+dzh(ksn)/2._ireals)-   &
                  MAX(z_old(k)-dz_old(k)/2._ireals,zm(ksn)-dzh(ksn)/2._ireals),0._ireals)/dzh(ksn)

      t_new(ksn)   = t_new(ksn)   + t(k)*weight(k)
      rho_new(ksn) = rho_new(ksn) + rho(k)*weight(k)
      wtr_new(ksn) = wtr_new(ksn) + wtr(k)*weight(k)
      wl_new(ksn)  = wl_new(ksn)  + wl(k)*weight(k)
    end do

  end do

  do ksn = 1,ke_snow
    t(ksn)   = t_new(ksn)
    rho(ksn) = rho_new(ksn)
    wtr(ksn) = wtr_new(ksn)
    wl(ksn)  = wl_new(ksn)
    wl(ksn) = wl(ksn)*dzh(ksn)
    wtr(ksn) = wtr(ksn)*dzh(ksn)
  end do

END SUBROUTINE normalize

!==============================================================================

!------------------------------------------------------------------------------
! End of module src_soil_multlay
!------------------------------------------------------------------------------




END MODULE mo_soil_ml
