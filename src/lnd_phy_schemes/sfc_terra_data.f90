!+ Data module for all parametric data in the soil model "terra"  
!------------------------------------------------------------------------------

MODULE sfc_terra_data

!------------------------------------------------------------------------------
!
! Description:
!  This module declares and initializes all parametric scalar and array      
!  data which are used in the soil model
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  Juergen.Helmert@dwd.de
!
! Old COSMO History: (from data_soil.f90)
! 1.1        1998/03/11 Guenther Doms
!  Initial release
! 1.30       1999/06/24 Erdmann Heise
! Implementation of variables for simplified BATS-scheme
! 1.33       1999/10/14 Matthias Raschendorfer
!  crsmin is now a namelist-parameter.
! 2.17       2002/05/08 Ulrich Schaettler
!  Additional parameters for soil water content dependent freezing/melting
! 2.18       2002/07/16 Reinhold Schrodin
!  Redefined variable cf_snow (for calculation of fractional snow coverage)
! 3.6        2003/12/11 Reinhold Schrodin
!  Adapted several variables for new multi-layer soil model
! 3.13       2004/12/03 Reinhold Schrodin
!  New variable cwimax_ml (maximum interception water content for multi-layer
!  soil model). Changed values for minimal and maximal density of snow
!  (crhosmin_ml: 100 => 250; crhosmax_ml: 400 => 250) (now consistent with GME)
! 3.17       2005/12/12 Reinhold Schrodin
!  New variables (crhosmin, crhosmaxf, crhosmin, crhosmaxt, csnow_tmin)
!  and changed variables (crhosmin_ml,crhosmax_ml) for ageing of snow density
!  calculation
! 3.18       2006/03/03 Ulrich Schaettler
!  Editorial changes
! 3.21       2006/12/04 Ulrich Schaettler
!  Put declaration of NL parameters crsmin and rat_lam to data_soil
! V3_23        2007/03/30 Matthias Raschendorfer
!  Moving 'rat_lam' into MODULE 'data_turbulence'.
!  Initialisation of 'crsmin' with a new default value of 150.0
! V4_11        2009/11/30 Ekaterina Machulskaya
!  Introduced parameters for the snow model
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced new logical lsoilinit_dfi to initialize soil variables after
!    a DFI forward launching
!  Changed the code owner
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V5_4e        2017-03-23 Ulrich Schaettler
!  Initial release for COSMO, including variables for TERRA_URB
! V5_4f        2017-09-01 Ulrich Schaettler, Valentin Clement
!  Introduced idiag_snowfrac as namelist variable (US)
!  Introduced OpenACC statements for data create / delete (VC)
! V5_4g        2017-11-13 Ulrich Schaettler
!  Adaptations to latest ICON updates (2017-11-08)
! V5_4h        2017-12-15 Ulrich Schaettler
!  Define zzhls, zdzhs, zdzms as global variables, which are computed in sfc_init
! V5_5         2018-02-23 Ulrich Schaettler
!  Modifications to run the full block of parameterizations on GPU
!   Removed soil and rock/ice lists
!   Added additional fields to ALLOC_WKARR and acc create/exit
! V5_6         2019-02-27 Ulrich Schaettler
!  Single precision version: replaced eps_soil where necessary
! V5_6a        2019-05-21 Jan-Peter Schulz
!  Introduce skin temperature approach
! @VERSION@    @DATE@     Ulrich Schaettler
!  Re-unification with ICON
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
#ifdef __COSMO__
USE kind_parameters, ONLY :   &
#elif __ICON__
USE mo_kind, ONLY:     &
#endif
    wp           ! KIND-type parameter for real variables

!==============================================================================

IMPLICIT NONE

PUBLIC           ! All constants and variables in this module are public

!==============================================================================

! Global (i.e. public) Declarations:

! 1. Data arrays for properties of different soil types (array index)     
! -------------------------------------------------------------------
 
  REAL  (KIND=wp)     ::  &
!   a) parameters describing the soil water budget
    cporv (10), &  !  pore volume (fraction of volume)
    cfcap (10), &  !  field capacity (fraction of volume)
    cpwp  (10), &  !  plant wilting point (fraction of volume)
    cadp  (10), &  !  air dryness point (fraction of volume)
    cik2  (10), &  !  minimum infiltration rate (kg/s*m**2)
    ckw0  (10), &  !  parameter for determination of hydr. conductivity (m/s)
    ckw1  (10), &  !  parameter for determination of hydr. conductivity (1)
    cdw0  (10), &  !  parameter for determination of hydr. diffusivity (m**2/s)
    cdw1  (10), &  !  parameter for determination of hydr. diffusivity (1)
    crock (10), &  !  rock/ice/water indicator (hydrological calculations 
                   !  only for crock=1)

!   b) parameters describing the soil heat budget
    cdz1  (10), &  !  top layer thickness (EFR-method)
    crhoc (10), &  !  soil heat capacity  (J/K*m**3)
    cala0 (10), &  !  parameters for the determination of
    cala1 (10), &  !      the soil heat conductivity (W/(K*m))
    csalb (10), &  !  solar albedo for dry soil                            
    csalbw(10), &  !  slope of solar albedo with respect to soil water content     

!   c) additional parameters for the BATS scheme (Dickinson)
    ck0di (10), &  !  (m/s)
    cbedi (10), &  !  (1)
    clgk0 (10), &  !  auxiliary variable

!   d) additional parameters for soil water content dependent freezing/melting
    csandf(10), &  !  mean fraction of sand (weight percent)
    cclayf(10)     !  mean fraction of clay (weight percent)
 
  !$acc declare copyin(clgk0)

  ! Initialization of soil type parameters except cdz1 
  ! (being calculated during execution)

  ! soil type:   ice    rock    sand    sandy   loam   clay      clay    peat    sea     sea  
  ! (by index)                          loam           loam                     water    ice

  DATA  cporv / 1.E-10_wp, 1.E-10_wp, 0.364_wp  , 0.445_wp  , 0.455_wp  , 0.475_wp  , 0.507_wp  , 0.863_wp  , 1.E-10_wp, 1.E-10_wp /
  DATA  cfcap / 1.E-10_wp, 1.E-10_wp, 0.196_wp  , 0.260_wp  , 0.340_wp  , 0.370_wp  , 0.463_wp  , 0.763_wp  , 1.E-10_wp, 1.E-10_wp /
  DATA  cpwp  / 0.0_wp   , 0.0_wp   , 0.042_wp  , 0.100_wp  , 0.110_wp  , 0.185_wp  , 0.257_wp  , 0.265_wp  , 0.0_wp   ,  0.0_wp   /
  DATA  cadp  / 0.0_wp   , 0.0_wp   , 0.012_wp  , 0.030_wp  , 0.035_wp  , 0.060_wp  , 0.065_wp  , 0.098_wp  , 0.0_wp   ,  0.0_wp   /
  DATA  crhoc / 1.92E6_wp, 2.10E6_wp, 1.28E6_wp , 1.35E6_wp , 1.42E6_wp , 1.50E6_wp , 1.63E6_wp , 0.58E6_wp , 4.18E6_wp, 1.92E6_wp /
  DATA  cik2  / 0.0_wp   , 0.0_wp   , 0.0035_wp , 0.0023_wp , 0.0010_wp , 0.0006_wp , 0.0001_wp , 0.0002_wp , 0.0_wp   ,  0.0_wp   /
#ifdef __COSMO__
  DATA  ckw0  / 0.0_wp   , 0.0_wp   , 479.E-7_wp, 943.E-8_wp, 531.E-8_wp, 764.E-9_wp, 17.E-9_wp , 58.E-9_wp , 0.0_wp   ,  0.0_wp   /
#endif
#ifdef __ICON__
  DATA  ckw0  / 0.0_wp   , 0.0_wp   , 479.E-7_wp, 943.E-8_wp, 531.E-8_wp, 764.E-9_wp, 85.E-9_wp , 58.E-9_wp , 0.0_wp   ,  0.0_wp   /
#endif
  DATA  ckw1  / 0.0_wp   , 0.0_wp   , -19.27_wp , -20.86_wp , -19.66_wp , -18.52_wp , -16.32_wp , -16.48_wp , 0.0_wp   ,  0.0_wp   /
  DATA  cdw0  / 0.0_wp   , 0.0_wp   , 184.E-7_wp, 346.E-8_wp, 357.E-8_wp, 118.E-8_wp, 442.E-9_wp, 106.E-9_wp, 0.0_wp   ,  0.0_wp   /
  DATA  cdw1  / 0.0_wp   , 0.0_wp   , -8.45_wp  , -9.47_wp  , -7.44_wp  , -7.76_wp  , -6.74_wp  , -5.97_wp  , 0.0_wp   ,  0.0_wp   /
  DATA  crock / 0.0_wp   , 0.0_wp   , 1.0_wp    , 1.0_wp    , 1.0_wp    , 1.0_wp    , 1.0_wp    , 1.0_wp    , 0.0_wp   ,  0.0_wp   /
  DATA  cala0 / 2.26_wp  , 2.41_wp  , 0.30_wp   , 0.28_wp   , 0.25_wp   , 0.21_wp   , 0.18_wp   , 0.06_wp   , 1.0_wp   ,  2.26_wp  /
  DATA  cala1 / 2.26_wp  , 2.41_wp  , 2.40_wp   , 2.40_wp   , 1.58_wp   , 1.55_wp   , 1.50_wp   , 0.50_wp   , 1.0_wp   ,  2.26_wp  /
  DATA  csalb / 0.70_wp  , 0.30_wp  , 0.30_wp   , 0.25_wp   , 0.25_wp   , 0.25_wp   , 0.25_wp   , 0.20_wp   , 0.07_wp  ,  0.70_wp  /
  DATA  csalbw/ 0.00_wp  , 0.00_wp  , 0.44_wp   , 0.27_wp   , 0.24_wp   , 0.23_wp   , 0.22_wp   , 0.10_wp   , 0.00_wp  ,  0.00_wp  /
  DATA  ck0di / 1.E-4_wp , 1.E-4_wp , 2.E-4_wp  , 2.E-5_wp  , 6.E-6_wp  , 2.E-6_wp  , 1.E-6_wp  , 1.5E-6_wp , 0.00_wp  ,  0.00_wp  /
  DATA  cbedi / 1.00_wp  , 1.00_wp  , 3.5_wp    , 4.8_wp    , 6.1_wp    , 8.6_wp    , 10.0_wp   , 9.0_wp    , 0.00_wp  ,  0.00_wp  /
  DATA  csandf/ 0.0_wp   , 0.0_wp   , 90._wp    , 65._wp    , 40._wp    , 35._wp    , 15._wp    , 90._wp    , 0.00_wp  ,  0.00_wp /
  DATA  cclayf/ 0.0_wp   , 0.0_wp   , 5.0_wp    , 10._wp    , 20._wp    , 35._wp    , 70._wp    , 5.0_wp    , 0.00_wp  ,  0.00_wp /
 

!==============================================================================
! Soiltype IDs
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: ist_seawtr = 9     ! ID of soiltype 'sea water'
  INTEGER, PARAMETER :: ist_seaice = 10    ! ID of soiltype 'sea ice'


! 2. Additional parameters for the soil model                             
! -------------------------------------------------------------------

  REAL  (KIND=wp)           ::  &
!==============================================================================

    csalb_p        = 0.15_wp  , & !  solar albedo of ground covered by plants
    csalb_snow     = 0.70_wp  , & !  solar albedo of ground covered by snow

#ifdef __COSMO__
    csalb_snow_min = 0.400_wp , & ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max = 0.700_wp , & ! max. solar albedo of snow for forest free surfaces
  ! for possible later use:
    csalb_snow_fe  = 0.200_wp , &  ! solar albedo of snow for surfaces with evergreen forest
    csalb_snow_fd  = 0.200_wp , &  ! solar albedo of snow for surfaces with deciduous forest
#elif __ICON__
    ! T.R. 2011-09-21 csalb_snow_min/max set to values used in GME
    csalb_snow_min = 0.500_wp , &
                           ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max = 0.850_wp , &
                           ! max. solar albedo of snow for forest free surfaces
    ! T.R. 2011-09-21 snow albedos for forests set to values used in GME
    csalb_snow_fe  = 0.270_wp , &  ! solar albedo of snow for surfaces with evergreen forest
    csalb_snow_fd  = 0.320_wp , &  ! solar albedo of snow for surfaces with deciduous forest
#endif

    ctalb          = 0.004_wp , & !  thermal albedo ( of all soil types )   
    cf_snow        = 0.0150_wp, & !  parameter for the calculation of the 
                                  !  fractional snow coverage
  ! for the multi-layer soil model
    cwhc       = 0.04_wp      , & !  water holding capacity of snow ()
    chcond     = 0.01_wp      , & !  saturation hydraulic conductivity of snow ()
    ca2        = 6.6E-07_wp   , & !  activation energy (for snow metamorphosis) (J)
    csigma     = 75._wp       , & !  snow metamorphosis, Pa

  ! cf_w changed from 0.0004 to 0.0010 (in agreement with GME)
    cf_w       = 0.0010_wp    , & !  parameter for the calculation of the
                                  !  fractional water coverage

    csvoro     = 1.0000_wp    , & !  parameter to estimate the subgrid-scale 
                                  !  variation of orography
    cik1       = 0.0020_wp    , & !  parameter for the determination of the 
                                  !  maximum infiltaration
#ifdef __COSMO__
    cwimax_ml                 , & !  maximum interception water content
              ! this is now a namelist variable in PHYCTL
    ctau_i     = 1000.0_wp    , & !  time constant for the drainage from the interception storage
#elif __ICON__
    ctau_i     = 7200.0_wp    , & !  time constant for the drainage from the interception storage
#endif
    cakw       = 0.8000_wp    , & !  parameter for averaging the water contents
                                  !  of the top and middle soil water layers to 
                                  !  calculate the hydraulic diffusivity and 
                                  !  conductiviy

    ctau1      = 1.0000_wp    , & !  first adjustment time period in EFR-method
    ctau2      = 5.0000_wp    , & !  second adjustment time period in EFR-method
    chc_i      = 2100.0_wp    , & !  heat capacity of ice     
    chc_w      = 4180.0_wp    , & !  heat capacity of water     

    cdzw12     = 0.1000_wp    , & !  thickness of upper soil water layer in 
                                  !  two-layer model         
    cdzw22     = 0.9000_wp    , & !  thickness of lower soil water layer in 
                                  !  two-layer model      
    cdzw13     = 0.0200_wp    , & !  thickness of upper soil water layer in 
                                  !  three-layer model
    cdzw23     = 0.0800_wp    , & !  thickness of middle soil water layer in 
                                  !  three-layer model 
    cdzw33     = 0.9000_wp        !  thickness of lower soil water layer in 
                                  !  three-layer model

  REAL  (KIND=wp)           ::  &
    cdsmin     = 0.0100_wp    , & !  minimum snow depth
    crhosmin   = 500.00_wp    , & !  minimum density of snow
    crhosmax   = 800.00_wp    , & !  maximum density of snow
    crhosmin_ml=  50.00_wp    , & !  minimum density of snow
    crhosmax_ml= 400.00_wp    , & !  maximum density of snow
    crhosminf  =  50.00_wp    , & !  minimum density of fresh snow
    crhosmaxf  = 150.00_wp    , & !  maximum density of fresh snow
    crhogminf  = 100.00_wp    , & !  minimum density of fresh graupel / convective snow
    crhogmaxf  = 200.00_wp    , & !  maximum density of fresh graupel / convective snow
#ifdef __COSMO__
    crhosmint  =   0.20_wp    , & !  minimum value of time constant for ageing of snow
#elif __ICON__
    crhosmint  =   0.125_wp    ,& !  value of time constant for ageing 
                                  !  of snow at csnow_tmin (8 days)
#endif
    crhosmaxt  =   0.40_wp    , & !  maximum value of time constant for ageing 
                                  !  of snow
    crhosmax_tmin = 200.00_wp , & ! maximum density of snow at csnow_tmin
    csnow_tmin = 258.15_wp    , & !  lower threshold temperature of snow for 
                                  !  ageing and fresh snow density computation 
                                  !  ( = 273.15-15.0)
    crhos_dw   = 300.00_wp    , & !  change of snow density with water content
    calasmin   = 0.2000_wp    , & !  minimum heat conductivity of snow (W/m K)
    calasmax   = 1.5000_wp    , & !  maximum heat conductivity of snow (W/m K)
    calas_dw   = 1.3000_wp    , & !  change of snow heat conductivity with
                                  !  water content                (W/(m**2) K)
   
    crhowm     =    0.8_wp    , & !  BATS (1)
    cdmin      =    0.25E-9_wp, & !  BATS (m**2/s)
    cfinull    =    0.2_wp    , & !  BATS (m)
    ckrdi      =    1.0E-5_wp , & !  BATS (m/s)
    cdash      =    0.05_wp   , & !  BATS ((m/s)**1/2)
    clai       =    3.0_wp    , & !  BATS
    cparcrit   =  100.0_wp    , & !  BATS (W/m**2)
    ctend      =  313.15_wp   , & !  BATS (K)
    csatdef    = 4000.0_wp    , & !  BATS (Pa)

    !Minimum and maximum value of stomatal resistance (s/m)
    !used by the Pen.-Mont. method for vegetation transpiration
    !(itype_trvg=2):
    crsmin     = 150.0_wp     , & !  BATS (s/m)
    crsmax     = 4000.0_wp        !  BATS (s/m)
    ! crsmax increased from 1000 to 4000 s/m (to reduce latent heat flux).

  ! Parameters for the skin temperature formulation

  REAL  (KIND=wp)           ::  &
    cskinc                    , & ! skin conductivity (W/m**2/K)
    cimpl                         ! stability parameter for the computation of the skin temperature

! 3. Additional variables for the soil geometry
! ---------------------------------------------

  ! these are allocated and computed in sfc_init
  REAL  (KIND=wp), ALLOCATABLE ::    &
    zzhls          (:)             , & ! depth of the half level soil layers in m
    zdzhs          (:)             , & ! layer thickness between half levels
    zdzms          (:)                 ! distance between main levels


! 4. Variables for TERRA_URB
! --------------------------

  ! Default urban fabric parameters are derived according to literature. 
  ! They are obtained/tested for 
  !      Toulouse, Basel (Wouters et al., 2015), 
  !      Paris (De Ridder et al., 2013; 
  !      Sarkar and De Ridder, 2010; 
  !      Demuzere et al., 2008), 

  REAL  (KIND=wp) ::       &
    ctalb_bm = 0.08_wp   , & !  default effective thermal albedo of building/road environment

    ! csalb_eff_uf = 0.80_wp,   & ! correction factor for effective albedo induced by the 
    !                             ! urban fabric (street canyons). 
    ! ! This value is based on observations and monte-carlo simulations, taking H/W ratio of 1.0 
    ! ! and roof fraction of 0.5, see Pawlak et al. 
    ! ! ????: http://nargeo.geo.uni.lodz.pl/~icuc5/text/P_4_6.pdf. csalb_eff_uf = (20% + 15%)/2

    csalb_bm = 0.213_wp  , & ! default short-wave albedo of building/road materials 
                             ! in the urban fabric
      ! csalb_bm is chosen in such a way that the effective albedo for a dense urban 
      ! environment csalb_bm * csalb_eff_uf is equal to 0.17

  ! Default surface area index of buildings/street environment. This is esimated from the 
  ! squared thermal inertia = 3.8E6 (for Paris, see De Ridder et al.,2013) estimated from 
  ! model simulations divided by the building material parameters below estimates 
  ! c_rhoc_bm * c_ala_bm 

    ! c_ai_uf = 2.0_wp,    &
    c_uf_h  = 15._wp     , & ! default height of building elements in urban fabric (metre)
    c_lnd_h = 0.01_wp    , & ! default height of natural soil elements in natural land (metre)
    c_htw   = 1.5_wp     , &
    c_roof  = 0.667_wp   , &

  ! default building/road 'material' properties...

    ! Value of 'concrete' is taken (see engineering-toolbox.com)
    c_rhoc_bm = 1.74E6_wp, & ! default specific heat times density of buildings, 

    !  Value of 'medium concrete' is taken (higher boundary, as average value for 
    !  concrete in,  see engineerintoolbox.com)
    c_ala_bm  = 0.87_wp  , & ! default building material heat conductivity: 
    c_isa_runoff = 1.0_wp, & ! default fraction of water exiting from the impervious 
                               ! surface leading to runoff. The remainder fraction is (potentially) 
                               ! infiltrated in the neighbouring natural soil (switched off by default). 
                               ! Some addaptation strategies (such as infiltration of roof water) 
                               ! will lead to values less than 1
    c_isa_delt = 0.12_wp , & ! The maximum wet-surface fraction delta_max, see Wouters et al., 2015
    c_isa_wmax = 1.31_wp , & ! Maximum amount of water that can be stored by impervious surfaces 
                               ! (estimate for urban areas, Toulouse centre), see Wouters et al., 2015

    ! constants for calculation of anthropogenic heat, see Flanner, 2009
    cb1        = 0.451_wp, & !
    cb2        = 0.8_wp  , & !
    csig       = 0.3_wp  , & ! #sigma =0.18
    cmu        = 0.5_wp  , & !
    cA1        = -0.5_wp , & ! #A1 = -0.3
    cff        = 2.0_wp  , & ! #ff = 1.9
    calph      = 10.0_wp , & !
    ceps       = 0.25_wp     !



! 5. Additional control variables
! -------------------------------

  LOGICAL                   ::  &
    lsoilinit_dfi = .FALSE.         ! initialize soil after dfi forward launching

! 6. Epsilons (security constants)
! --------------------------------

  REAL  (KIND=wp), PARAMETER :: &

    ! Avoid division by zero, e.g. x = y / MAX(z,eps_div).
    eps_div  = 1.0E-6_wp      , &
!!    eps_div  = repsilon       , &
!! RUS
!! eps_div is used in divisions to avoid division by zero, repsilon would
!! be appropriate therefore. However, the testsuite fails with repsilon (about 1E-30)
!! as it is 'used to' the (in this context) huge epsilon of 1E-6.

    ! Multi-purpose epsilon in soil model (former zepsi).
    eps_soil = 1.0E-6_wp      , &

    ! Small value to check if temperatures have exceeded a fixed threshold
    ! such as the freezing point.  In double precision (15 decimal digits)
    ! a value as small value such as 1.0E-6 can be used. In single
    ! precision (6-7 decimal digits), however, the value has to be larger
    ! in order not to vanish. The current formulation is save for
    ! temperatures up to 500K.
    eps_temp = MAX(1.0E-6_wp,500.0_wp*EPSILON(1.0_wp))


#ifdef __COSMO__
! 7. Taken from ICON (at least for the moment being: ICON declares these variables somewhere different)
! --------------------------------

  REAL  (KIND=wp)            ::  &
    max_toplaydepth  =  0.25_wp, & ! maximum depth of uppermost snow layer for multi-layer snow scheme (25 cm)
    tune_minsnowfrac = 0.125_wp    ! Minimum value to which the snow cover fraction is artificially reduced
                                   ! in case of melting show (in case of idiag_snowfrac = 20/30/40)

  INTEGER                    ::  &
    itype_interception = 1,      & ! type of plant interception
    idiag_snowfrac                 ! method for diagnosis of snow-cover fraction

  LOGICAL                    ::  &
    l2lay_rho_snow = .FALSE.       ! use two-layer snow density for single-layer snow model
#endif

! 8. Soil ice parameterization
! ----------------------------

  REAL  (KIND=wp), PARAMETER ::  &
    T_ref_ice  =  0.1_wp,        & !degC Soil ice parameterization
    T_star_ice =  0.01_wp,       & !degC according to K. Schaefer and Jafarov, E.,2016, doi:10.5194/bg-13-1991-2016
    b_clay     = -0.3_wp,        &
    b_silt     = -0.5_wp,        &
    b_sand     = -0.9_wp,        &
    b_org      = -1.0_wp


#ifdef ALLOC_WKARR
! Local arrays defined as allocatables here:
! ------------------------------------------

  INTEGER, ALLOCATABLE           ::  &
    m_styp         (:)             , & ! soil type
    ke_soil_hy_b   (:)             , &
    zicount1       (:)             , &
    zicount2       (:)

  REAL(KIND = wp), ALLOCATABLE    :: &

    ! Two-time level variables exchange with interface
    h_snow_now     (:)             , & ! snow height  (m)
    h_snow_new     (:)             , & ! snow height  (m)

    ! Model geometry
    zdz_snow_fl    (:)             , & ! snow depth for snow temperature gradient

    ! Multi-layer snow model
    zhh_snow       (:,:)           , & ! depth of the half level snow layers
    zhm_snow       (:,:)           , & ! depth of snow main levels
    zdzh_snow      (:,:)           , & ! layer thickness between half levels
    zdzm_snow      (:,:)           , & ! distance between main levels
    zdtsnowdt_mult (:,:)           , & ! tendency of ztsnow

    ! External parameters
    zbwt           (:)             , & ! root depth (with artificial minimum value)
    zrock          (:)             , & ! ice/rock-indicator: 0 for ice and rock
    zsandf         (:)             , & ! mean fraction of sand (weight percent)
    zclayf         (:)             , & ! mean fraction of clay (weight percent)
    zsiltf         (:)             , & ! mean fraction of silt (weight percent)
    zb_por         (:)             , & ! pore size distribution index
    zpsis          (:)             , & ! air entry potential (m)
    zw_m_org       (:)             , & ! maximum of liquid water content organic
    zw_m_soil      (:)             , & ! maximum of liquid water content mineral soil
    zw_m_up        (:)             , & ! maximum of liquid water content at temp -3 degC
    zw_m_low       (:,:)           , & ! maximum of liquid water content at temp -40 degC
    zw_m           (:)                 ! maximum of liquid water content  (m)

  REAL(KIND=wp), ALLOCATABLE     ::  &

    ! Connection to the atmosphere
    zrr            (:)             , & ! total rain rate including formation of dew
    zrs            (:)             , & ! total snow rate including formation of rime
    zesoil         (:)             , & ! evaporation from bare soil
    zsfc_frac_bs   (:)             , & ! evaporation from bare soil
    zrhoch         (:)             , & ! transfer coefficient*rho*g
    zth_low        (:)             , & ! potential temperature of lowest layer
    zf_wi          (:)             , & ! surface fraction covered by interception water
    ztmch          (:)             , & ! heat transfer coefficient*density*velocity
    zep_s          (:)             , & ! potential evaporation for t_s
    zep_snow       (:)             , & ! potential evaporation for t_snow
    zverbo         (:)             , & ! total evapotranspiration
    zversn         (:)             , & ! total evaporation from snow surface
    zthsoi         (:)             , & ! thermal flux at soil surface
    zthsnw         (:)             , & ! thermal flux at snow surface
    zfor_s         (:)             , & ! total forcing at soil surface
    zgsb           (:)             , & ! heat-flux through snow
    zrnet_s        (:)             , & ! net radiation
    zsprs          (:)             , & ! utility variable

    ! Tendencies
    zdwidt         (:)             , & ! interception store tendency
    zdwsndt        (:)             , & ! snow water content tendency
    zdtsdt         (:)             , & ! tendency of zts
    zdtsnowdt      (:)             , & ! tendency of ztsnow
    zdwgdt         (:,:)               ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=wp), ALLOCATABLE ::  &

    !   Interception variables
    zwinstr        (:)             , & ! preliminary value of interception store
    zinfmx         (:)             , & ! maximum infiltration rate
    zwimax         (:)             , & ! maximum interception store
    zvers          (:)             , & ! water supply for infiltration
    zwisnstr       (:)             , & ! water content of snow interception store (t+1) (mH2O)
    zwpnstr        (:)             , & ! water content of pond store (t+1) (mH2O)
    zfd_wi         (:)             , & ! surface fraction of intercepted water
    zf_pd          (:)             , & ! surface fraction covered by pond
    zwisn          (:)             , & ! water content of snow interception store (mH2O)
    zwpn           (:)             , & ! water content of pond
    zewi           (:)             , & ! evaporation from interception store
    zepd           (:)             , & ! evaporation from pond
    zept           (:)             , & ! potential evaporation
    zesn           (:)             , & ! evaporation from snow
    zdrr           (:)             , & ! formation of dew
    zrrs           (:)                 ! formation of rime

  REAL    (KIND=wp), ALLOCATABLE ::  &
    zg1            (:)             , & ! heat flux between two uppermost soil layers (W/m**2)
    zlhfl          (:)             , & ! estimated       latent    heat flux surface (W/m**2)
    zshfl          (:)             , & ! estimated       sensible  heat flux surface (W/m**2)
    zthfl          (:)             , & ! estimated total turbulent heat flux surface (W/m**2)
    zradfl         (:)             , & ! total radiative flux at surface             (W/m**2)
    ze_melt        (:)             , & ! energy required for melting of snow         (J/m**2)
    zch_snow       (:)             , & ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
    zeb1           (:)             , & ! estimated energy budget of first soil layer (W/m**2)
    ztchv          (:)             , & ! transfer coefficient multiplied by wind velocity
    ztchv_max      (:)             , & ! dito          , but upper limit
    zrho_atm       (:)             , & ! air density of lowest atmospheric layer     (kg/m**3)
    zdt_atm        (:)             , & ! temperature difference between lowest
    zdq_atm        (:)             , & ! dito for moisture (kg/kg)

    !additional variables for root distribution
    zroota         (:)             , & ! root density profile parameter (1/m)
    zwrootdz       (:,:)           , & ! mean water content over root depth weighted by root density
    zrootdz_int    (:)             , & ! parameter needed to initialize the root density profile integral
    zwrootdz_int   (:)             , & ! parameter needed to initialize the root water content integral
    zqhfl_s        (:)             , & ! moisture flux at soil/air interface
    zqhfl_snow     (:)                 ! moisture flux at snow/air interface

  REAL    (KIND=wp), ALLOCATABLE ::  &

    ! Soil and plant parameters
    zroc     (:,:)                 , & ! heat capacity
    zfcap    (:,:)                 , & ! field capacity of soil
    zadp     (:,:)                 , & ! air dryness point
    zporv    (:,:)                 , & ! pore volume (fraction of volume)
    zdlam    (:)                   , & ! heat conductivity parameter
    zdw      (:,:)                 , & ! hydrological diff.parameter
    zdw1     (:,:)                 , & ! hydrological diff.parameter
    zkw0     (:)                   , & ! hydrological cond.parameter at surface
    zkw      (:,:)                 , & ! hydrological cond.parameter
    zkwm     (:,:)                 , & ! hydrological cond.parameter (main levels)
    zkw1     (:,:)                 , & ! hydrological cond.parameter
    zik2     (:)                   , & ! minimum infiltration rate
    zpwp     (:,:)                 , & ! plant wilting point  (fraction of volume)
    ztlpmwp  (:)                   , & ! turgor-loss-point minus plant wilting point
    zedb     (:)                   , & ! utility variable
    zaa      (:)                   , & ! utility variable

    ! Hydraulic variables
    ztrang      (:,:)              , & ! transpiration contribution by the different layers
    ztrangs     (:)                , & ! total transpiration (transpiration from all
                                       !    soil layers)
    zwin        (:)                , & ! water cont. of interception store   (m H20)
    zwsnow      (:)                , & ! snow water equivalent               (m H20)
    zwsnew      (:)                , & ! snow water equivalent               (m H20)
    zdwsnm      (:)                , & ! utility variable for snow melt determination
    zw_fr       (:,:)              , & !fractional total water content of soil layers
    zinfil      (:)                , & ! infiltration rate
    zlw_fr      (:,:)              , & ! fractional liqu. water content of soil layer
    ziw_fr      (:,:)              , & ! fractional ice content of soil layer
    zwsnn       (:)                , & ! new value of zwsnow
    zflmg       (:,:)              , & ! flux of water at soil layer interfaces
    zrunoff_grav(:,:)                  ! main level water gravitation
                                       !    flux (for runoff calc.)

  REAL    (KIND=wp), ALLOCATABLE ::  &

    ! Local arrays for BATS-scheme
    zk0di       (:)                , & ! surface type dependent parameter
    zbedi       (:)                , & ! surface type dependent parameter
    zsnull      (:)                , & ! mean relative(zporv) water fraction of the
                                       !    so-called total active soil layer (above 1m)
    zs1         (:)                , & ! mean relative (zporv) water fraction  of
                                       !    layers above 0.1m
    zf_rad      (:)                , & ! radiation function for stomata resistance
    ztraleav    (:)                , & ! transpiration rate of dry leaves

    ! Root functions
    zwroot      (:)                , & ! mean water content over root depth
    zropartw    (:,:)                  ! fraction of layer filled by roots * w_fr

  REAL    (KIND=wp), ALLOCATABLE ::  &

    ! Thermal variables
    zts         (:)                , & ! soil surface temperature
    ztsk        (:)                , & ! skin temperature
    ztsnow      (:)                , & ! snow surface temperaure
    ztsnow_mult (:,:)              , & ! snow surface temperaure

    ! HEATCOND (soil moisture dependent heat conductivity of the soil)
    zalamtmp    (:,:)              , & ! heat conductivity
    zalam       (:,:)              , & ! heat conductivity
    zrocg       (:,:)              , & ! total volumetric heat capacity of soil
    zrocg_soil  (:,:)              , & ! volumetric heat capacity of bare soil
    zrocs       (:)                , & ! heat capacity of snow
    ztsn        (:)                , & ! new value of zts
    ztskn       (:)                , & ! new value of ztsk
    ztsnown     (:)                , & ! new value of ztsnow
    ztsnown_mult(:,:)              , & ! new value of ztsnow
    znlgw1f     (:)                    ! utility variable

  REAL(KIND=wp), ALLOCATABLE     ::  &

    !   Multi-snow layer parameters
    zqbase       (:)               , &
    zrefr        (:)               , & ! rate of liquid water refreezing
    zmelt        (:)               , & ! rate of snow melting
    ze_out       (:)               , &
    zrho_dry_old (:)               , &
    zp           (:,:)             , &
    zcounter     (:)               , &
    ze_rad       (:)               , &
    zswitch      (:)               , &
    tmp_num      (:)               , &
    sum_weight   (:)               , &
    t_new        (:,:)             , &
    rho_new      (:,:)             , &
    wl_new       (:,:)             , &
    z_old        (:,:)             , &
    dz_old       (:,:)             , &
    t_so_free_new(:,:)             , &
    t_so_snow_new(:,:)             , &
    sn_frac      (:)               , &
    zf_snow_lim  (:)               , &
    zdz_snow     (:)               , & ! snow depth (not for snow heat content but snow
                                       !             density calculation)
    zalas_mult    (:,:)            , & ! heat conductivity of snow
    ztsnownew_mult(:,:)            , & ! preliminary value of snow surface temperature
    zextinct      (:,:)            , & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(:)                  ! total forcing at snow surface

  REAL    (KIND=wp), ALLOCATABLE ::  &
    ! Auxiliary variables
    hzalam      (:,:)              , & ! heat conductivity
    zdqvtsnow   (:)                , & ! first derivative of saturation specific humidity
                                       !    with respect to t_snow
    zrho_snow   (:)                , & ! snow density used for computing heat capacity and conductivity
    zts_pm      (:)                , & ! indicator zts > < T_melt
    ztsk_pm     (:)                , & ! indicator ztsk > < T_melt
    ztfunc      (:)                , & ! smoothed transition function between T_melt and T_melt+2K
    ztsnow_pm   (:)                , & ! indicator ztsnow > < T_melt
    zeisa       (:)                    ! evaporation from impervious surfaces (puddles)

  REAL    (KIND=wp), ALLOCATABLE ::  &
    ! for terra_init
    zw_snow_old   (:)              , & !
    zrho_snow_old (:)              , & !
    h_snow_fg     (:)              , & !
    h_snow_incr   (:)                  !

  REAL    (KIND=wp), ALLOCATABLE ::  &
    ! Utility variables
    zaga        (:,:)              , &
    zagb        (:,:)              , &
    zagc        (:,:)              , &
    zagd        (:,:)              , &
    zage        (:,:)


  LOGICAL,           ALLOCATABLE  :: &
    limit_tch (:)                  , & ! indicator for flux limitation problem
    lzurban   (:)                  , & ! indicator for flux limitation problem
    l_redist  (:)

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE terra_wkarr_alloc (ke_soil, ke_snow, nproma, istat)

  INTEGER, INTENT(IN)  :: ke_soil, ke_snow, nproma
  INTEGER, INTENT(OUT) :: istat

  istat = 0

  ALLOCATE (m_styp         (nproma)          , & ! soil type
            ke_soil_hy_b   (nproma)          , & ! soil type
            zicount1       (nproma)          , & !
            zicount2       (nproma)          , & !
       STAT=istat)


  ALLOCATE (                                   &
            ! Two-time level variables exchange with interface
            h_snow_now     (nproma)          , & ! snow height  (m)
            h_snow_new     (nproma)          , & ! snow height  (m)

            ! Model geometry
            zdz_snow_fl    (nproma)          , & ! snow depth for snow temperature gradient

            ! Multi-layer snow model
            zhh_snow       (nproma,  ke_snow), & ! depth of the half level snow layers
            zhm_snow       (nproma,  ke_snow), & ! depth of snow main levels
            zdzh_snow      (nproma,  ke_snow), & ! layer thickness between half levels
            zdzm_snow      (nproma,  ke_snow), & ! distance between main levels
            zdtsnowdt_mult (nproma,0:ke_snow), & ! tendency of ztsnow

            ! External parameters
            zbwt           (nproma)          , & ! root depth (with artificial minimum value)
            zrock          (nproma)          , & ! ice/rock-indicator: 0 for ice and rock
            zsandf         (nproma)          , & ! mean fraction of sand (weight percent)
            zclayf         (nproma)          , & ! mean fraction of clay (weight percent)
            zsiltf         (nproma)          , & ! mean fraction of silt (weight percent)
            zb_por         (nproma)          , & ! pore size distribution index
            zpsis          (nproma)          , & ! air entry potential (m)
            zw_m_org       (nproma)          , & ! maximum of liquid water content organic
            zw_m_soil      (nproma)          , & ! maximum of liquid water content mineral soil
            zw_m_up        (nproma)          , & ! maximum of liquid water content at temp -3 degC
            zw_m_low       (nproma,ke_soil+1), & ! maximum of liquid water content at temp -40 degC
            zw_m           (nproma)          , & ! maximum of liquid water content  (m)
       STAT=istat)

  ALLOCATE (                                   &
            ! Connection to the atmosphere
            zrr            (nproma)          , & ! total rain rate including formation of dew
            zrs            (nproma)          , & ! total snow rate including formation of rime
            zesoil         (nproma)          , & ! evaporation from bare soil
            zsfc_frac_bs   (nproma)          , & ! evaporation from bare soil
            zrhoch         (nproma)          , & ! transfer coefficient*rho*g
            zth_low        (nproma)          , & ! potential temperature of lowest layer
            zf_wi          (nproma)          , & ! surface fraction covered by interception water
            ztmch          (nproma)          , & ! heat transfer coefficient*density*velocity
            zep_s          (nproma)          , & ! potential evaporation for t_s
            zep_snow       (nproma)          , & ! potential evaporation for t_snow
            zverbo         (nproma)          , & ! total evapotranspiration
            zversn         (nproma)          , & ! total evaporation from snow surface
            zthsoi         (nproma)          , & ! thermal flux at soil surface
            zthsnw         (nproma)          , & ! thermal flux at snow surface
            zfor_s         (nproma)          , & ! total forcing at soil surface
            zgsb           (nproma)          , & ! heat-flux through snow
            zrnet_s        (nproma)          , & ! net radiation
            zsprs          (nproma)          , & ! utility variable

            ! Tendencies
            zdwidt         (nproma)          , & ! interception store tendency
            zdwsndt        (nproma)          , & ! snow water content tendency
            zdtsdt         (nproma)          , & ! tendency of zts
            zdtsnowdt      (nproma)          , & ! tendency of ztsnow
            zdwgdt         (nproma,  ke_soil), & ! tendency of water content [kg/(m**3 s)]

            !   Interception variables
            zwinstr        (nproma)          , & ! preliminary value of interception store
            zinfmx         (nproma)          , & ! maximum infiltration rate
            zwimax         (nproma)          , & ! maximum interception store
            zvers          (nproma)          , & ! water supply for infiltration
            zwisnstr       (nproma)          , & ! water content of snow interception store (t+1) (mH2O)
            zwpnstr        (nproma)          , & ! water content of pond store (t+1) (mH2O)
            zfd_wi         (nproma)          , & ! surface fraction of intercepted water
            zf_pd          (nproma)          , & ! surface fraction covered by pond
            zwisn          (nproma)          , & ! water content of snow interception store (mH2O)
            zwpn           (nproma)          , & ! water content of pond
            zewi           (nproma)          , & ! evaporation from interception store
            zepd           (nproma)          , & ! evaporation from pond
            zept           (nproma)          , & ! potential evaporation
            zesn           (nproma)          , & ! evaporation from snow
            zdrr           (nproma)          , & ! formation of dew
            zrrs           (nproma)          , & ! formation of rime

            zg1            (nproma)          , & ! heat flux between two uppermost soil layers (W/m**2)
            zlhfl          (nproma)          , & ! estimated       latent    heat flux surface (W/m**2)
            zshfl          (nproma)          , & ! estimated       sensible  heat flux surface (W/m**2)
            zthfl          (nproma)          , & ! estimated total turbulent heat flux surface (W/m**2)
            zradfl         (nproma)          , & ! total radiative flux at surface             (W/m**2)
            ze_melt        (nproma)          , & ! energy required for melting of snow         (J/m**2)
            zch_snow       (nproma)          , & ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
            zeb1           (nproma)          , & ! estimated energy budget of first soil layer (W/m**2)
            ztchv          (nproma)          , & ! transfer coefficient multiplied by wind velocity
            ztchv_max      (nproma)          , & ! dito          , but upper limit
            zrho_atm       (nproma)          , & ! air density of lowest atmospheric layer     (kg/m**3)
            zdt_atm        (nproma)          , & ! temperature difference between lowest
            zdq_atm        (nproma)          , & ! dito for moisture (kg/kg)

            !additional variables for root distribution
            zroota         (nproma)          , & ! root density profile parameter (1/m)
            zwrootdz       (nproma, ke_soil) , & ! mean water content over root depth weighted by root density
            zrootdz_int    (nproma)          , & ! parameter needed to initialize the root density profile integral
            zwrootdz_int   (nproma)          , & ! parameter needed to initialize the root water content integral
            zqhfl_s        (nproma)          , & ! moisture flux at soil/air interface
            zqhfl_snow     (nproma)          , & ! moisture flux at snow/air interface
       STAT=istat)

  ALLOCATE (                                   &
            ! Soil and plant parameters
            zroc     (nproma,ke_soil+1)      , & ! heat capacity
            zfcap    (nproma,ke_soil+1)      , & ! field capacity of soil
            zadp     (nproma,ke_soil+1)      , & ! air dryness point
            zporv    (nproma,ke_soil+1)      , & ! pore volume (fraction of volume)
            zdlam    (nproma)                , & ! heat conductivity parameter
            zdw      (nproma,ke_soil+1)      , & ! hydrological diff.parameter
            zdw1     (nproma,ke_soil+1)      , & ! hydrological diff.parameter
            zkw0     (nproma)                , & ! hydrological cond.parameter at surface
            zkw      (nproma,ke_soil+1)      , & ! hydrological cond.parameter
            zkwm     (nproma,ke_soil+1)      , & ! hydrological cond.parameter (main levels)
            zkw1     (nproma,ke_soil+1)      , & ! hydrological cond.parameter
            zik2     (nproma)                , & ! minimum infiltration rate
            zpwp     (nproma,ke_soil+1)      , & ! plant wilting point  (fraction of volume)
            ztlpmwp  (nproma)                , & ! turgor-loss-point minus plant wilting point
            zedb     (nproma)                , & ! utility variable
            zaa      (nproma)                , & ! utility variable

            ! Hydraulic variables
            ztrang      (nproma,ke_soil)     , & ! transpiration contribution by the different layers
            ztrangs     (nproma)             , & ! total transpiration (transpiration from all
                                                 !    soil layers)
            zwin        (nproma)             , & ! water cont. of interception store   (m H20)
            zwsnow      (nproma)             , & ! snow water equivalent               (m H20)
            zwsnew      (nproma)             , & ! snow water equivalent               (m H20)
            zdwsnm      (nproma)             , & ! utility variable for snow melt determination
            zw_fr       (nproma,ke_soil+1)   , & !fractional total water content of soil layers
            zinfil      (nproma)             , & ! infiltration rate
            zlw_fr      (nproma,ke_soil+1)   , & ! fractional liqu. water content of soil layer
            ziw_fr      (nproma,ke_soil+1)   , & ! fractional ice content of soil layer
            zwsnn       (nproma)             , & ! new value of zwsnow
            zflmg       (nproma,ke_soil+1)   , & ! flux of water at soil layer interfaces
            zrunoff_grav(nproma,ke_soil+1)   , & ! main level water gravitation
                                                 !    flux (for runoff calc.)

            ! Local arrays for BATS-scheme
            zk0di       (nproma)             , & ! surface type dependent parameter
            zbedi       (nproma)             , & ! surface type dependent parameter
            zsnull      (nproma)             , & ! mean relative(zporv) water fraction of the
                                                 !    so-called total active soil layer (above 1m)
            zs1         (nproma)             , & ! mean relative (zporv) water fraction  of
                                                 !    layers above 0.1m
            zf_rad      (nproma)             , & ! radiation function for stomata resistance
            ztraleav    (nproma)             , & ! transpiration rate of dry leaves

            ! Root functions
            zwroot      (nproma)             , & ! mean water content over root depth
            zropartw    (nproma,ke_soil)     , & ! fraction of layer filled by roots * w_fr

            ! Thermal variables
            zts         (nproma)             , & ! soil surface temperature
            ztsk        (nproma)             , & ! skin temperature
            ztsnow      (nproma)             , & ! snow surface temperaure
            ztsnow_mult (nproma,0:ke_snow)   , & ! snow surface temperaure

            ! HEATCOND (soil moisture dependent heat conductivity of the soil)
            zalamtmp    (nproma,ke_soil)     , & ! heat conductivity
            zalam       (nproma,ke_soil)     , & ! heat conductivity
            zrocg       (nproma,ke_soil+1)   , & ! total volumetric heat capacity of soil
            zrocg_soil  (nproma,ke_soil+1)   , & ! volumetric heat capacity of bare soil
            zrocs       (nproma)             , & ! heat capacity of snow
            ztsn        (nproma)             , & ! new value of zts
            ztskn       (nproma)             , & ! new value of ztsk
            ztsnown     (nproma)             , & ! new value of ztsnow
            ztsnown_mult(nproma,0:ke_snow)   , & ! new value of ztsnow
            znlgw1f     (ke_soil)            , & ! utility variable
       STAT=istat)


  ALLOCATE (                                   &
            !   Multi-snow layer parameters
            zqbase       (nproma)            , &
            zrefr        (nproma)            , & ! rate of liquid water refreezing
            zmelt        (nproma)            , & ! rate of snow melting
            ze_out       (nproma)            , &
            zrho_dry_old (nproma)            , &
            zp           (nproma,ke_snow)    , &
            zcounter     (nproma)            , &
            ze_rad       (nproma)            , &
            zswitch      (nproma)            , &
            tmp_num      (nproma)            , &
            sum_weight   (nproma)            , &
            t_new        (nproma,ke_snow)    , &
            rho_new      (nproma,ke_snow)    , &
            wl_new       (nproma,ke_snow)    , &
            z_old        (nproma,ke_snow)    , &
            dz_old       (nproma,ke_snow)    , &
            t_so_free_new(nproma,0:ke_soil+1), &
            t_so_snow_new(nproma,0:ke_soil+1), &
            sn_frac      (nproma)            , &
            zf_snow_lim  (nproma)            , &
            zdz_snow     (nproma)            , & ! snow depth (not for snow heat content but snow
                                                 !             density calculation)
            zalas_mult    (nproma,  ke_snow) , & ! heat conductivity of snow
            ztsnownew_mult(nproma,0:ke_snow) , & ! preliminary value of snow surface temperature
            zextinct      (nproma,  ke_snow) , & ! solar radiation extinction coefficient in snow (1/m)
            zfor_snow_mult(nproma)           , & ! total forcing at snow surface

            ! Auxiliary variables
            hzalam      (nproma,ke_soil+1)   , & ! heat conductivity
            zdqvtsnow   (nproma)             , & ! first derivative of saturation specific humidity
                                                 !    with respect to t_snow
            zrho_snow   (nproma)             , & ! snow density used for computing heat capacity and conductivity
            zts_pm      (nproma)             , & ! indicator zts > < T_melt
            ztsk_pm     (nproma)             , & ! indicator ztsk > < T_melt
            ztfunc      (nproma)             , & ! smoothed transition function between T_melt and T_melt+2K
            ztsnow_pm   (nproma)             , & ! indicator ztsnow > < T_melt
            zeisa       (nproma)             , & ! evaporation from impervious surfaces (puddles)
       STAT=istat)

  ALLOCATE (                                 &
            ! for terra_init
            zw_snow_old  (nproma)             , &
            zrho_snow_old(nproma)             , &
            h_snow_fg    (nproma)             , &
            h_snow_incr  (nproma)             , &
       STAT=istat)

  ALLOCATE (                                 &
            ! Utility variables
            zaga        (nproma,0:ke_soil+ke_snow+1), &
            zagb        (nproma,0:ke_soil+ke_snow+1), &
            zagc        (nproma,0:ke_soil+ke_snow+1), &
            zagd        (nproma,0:ke_soil+ke_snow+1), &
            zage        (nproma,0:ke_soil+ke_snow+1), &
       STAT=istat)

  ALLOCATE (                                 &  ! logical
            limit_tch (nproma)             , &  ! indicator for flux limitation problem
            lzurban   (nproma)             , &  ! indicator for flux limitation problem
            l_redist  (nproma)             , &  !
       STAT=istat)

  !$acc enter data create(m_styp,ke_soil_hy_b)                      &
  !$acc create(zicount1, zicount2)                                  &
  !$acc create(h_snow_now, h_snow_new)                              &
  !$acc create(zdz_snow_fl, zhh_snow, zhm_snow, zdzh_snow)          &
  !$acc create(zdzm_snow, zdtsnowdt_mult, zbwt, zrock, zsandf)      &
  !$acc create(zclayf, zsiltf, zb_por, zpsis, zw_m)                 &
  !$acc create(zrr, zrs,zesoil, zsfc_frac_bs)                       &
  !$acc create(zw_m_org, zw_m_soil, zw_m_up, zw_m_low)              &
  !$acc create(zrhoch, zth_low, zf_wi, ztmch, zep_s, zep_snow)      &
  !$acc create(zverbo, zversn, zthsoi, zthsnw, zfor_s, zgsb)        &
  !$acc create(zrnet_s, zsprs, zdwidt, zdwsndt, zdtsdt)             &
  !$acc create(zdtsnowdt, zdwgdt, zwinstr, zinfmx, zwimax, zvers)   &
  !$acc create(zwisnstr, zwpnstr, zfd_wi, zf_pd, zwisn, zwpn)       &
  !$acc create(zewi, zepd, zept, zesn, zdrr, zrrs, zg1, zlhfl)      &
  !$acc create(zshfl, zthfl, zradfl, ze_melt, zch_snow, zeb1)       &
  !$acc create(ztchv, ztchv_max, zrho_atm, zdt_atm, zdq_atm)        &
  !$acc create(zroota, zwrootdz, zwrootdz_int, zrootdz_int)         &
  !$acc create(zqhfl_s, zqhfl_snow, zroc, zfcap, zadp, zporv)       &
  !$acc create(zdlam, zdw, zdw1, zkw, zkw1, zik2, zpwp, ztlpmwp)    &
  !$acc create(zkw0, zkwm)                                          &
  !$acc create(zedb, ztrang, ztrangs, zwin, zwsnow, zwsnew, zdwsnm) &
  !$acc create(zaa, zw_fr, zinfil, zlw_fr, ziw_fr, zwsnn, zflmg)    &
  !$acc create(zrunoff_grav, zk0di, zbedi, zsnull, zs1, zf_rad)     &
  !$acc create(ztraleav, zwroot, zropartw, zts, ztsk, ztsnow)       &
  !$acc create(ztsnow_mult, zalamtmp, zalam, zrocg, zrocg_soil)     &
  !$acc create(zrocs, ztsn, ztskn, ztsnown, ztsnown_mult, znlgw1f)  &
  !$acc create(zqbase, zrefr, zmelt, ze_out, zrho_dry_old)          &
  !$acc create(zp, zcounter, ze_rad, zswitch, tmp_num, sum_weight)  &
  !$acc create(t_new, rho_new, wl_new, dz_old, z_old)               &
  !$acc create(t_so_free_new, t_so_snow_new, sn_frac)               &
  !$acc create(zf_snow_lim, zdz_snow, zalas_mult, ztsnownew_mult)   &
  !$acc create(zextinct, zfor_snow_mult, hzalam, zdqvtsnow)         &
  !$acc create(zrho_snow, zts_pm, ztsk_pm, ztfunc, ztsnow_pm, zeisa)&
  !$acc create(zw_snow_old, zrho_snow_old, h_snow_fg, h_snow_incr)  &
  !$acc create(zaga, zagb, zagc, zagd, zage)                        &
  !$acc create(limit_tch, lzurban, l_redist)

END SUBROUTINE terra_wkarr_alloc

!==============================================================================
!==============================================================================

SUBROUTINE terra_wkarr_dealloc (istat)

  INTEGER, INTENT(OUT) :: istat

  istat = 0

  !$acc exit data delete(m_styp,ke_soil_hy_b)                       &
  !$acc delete(zicount1, zicount2)                                  &
  !$acc delete(h_snow_now, h_snow_new)                              &
  !$acc delete(zdz_snow_fl, zhh_snow, zhm_snow, zdzh_snow)          &
  !$acc delete(zdzm_snow, zdtsnowdt_mult, zbwt, zrock, zsandf)      &
  !$acc delete(zclayf, zsiltf, zb_por, zpsis, zw_m)                 &
  !$acc delete(zrr, zrs,zesoil, zsfc_frac_bs)                       &
  !$acc delete(zw_m_org, zw_m_soil, zw_m_up, zw_m_low)              &
  !$acc delete(zrhoch, zth_low, zf_wi, ztmch, zep_s, zep_snow)      &
  !$acc delete(zverbo, zversn, zthsoi, zthsnw, zfor_s, zgsb)        &
  !$acc delete(zrnet_s, zsprs, zdwidt, zdwsndt, zdtsdt)             &
  !$acc delete(zdtsnowdt, zdwgdt, zwinstr, zinfmx, zwimax, zvers)   &
  !$acc delete(zwisnstr, zwpnstr, zfd_wi, zf_pd, zwisn, zwpn)       &
  !$acc delete(zewi, zepd, zept, zesn, zdrr, zrrs, zg1, zlhfl)      &
  !$acc delete(zshfl, zthfl, zradfl, ze_melt, zch_snow, zeb1)       &
  !$acc delete(ztchv, ztchv_max, zrho_atm, zdt_atm, zdq_atm)        &
  !$acc delete(zroota, zwrootdz, zwrootdz_int, zrootdz_int)         &
  !$acc delete(zqhfl_s, zqhfl_snow, zroc, zfcap, zadp, zporv)       &
  !$acc delete(zdlam, zdw, zdw1, zkw, zkw1, zik2, zpwp, ztlpmwp)    &
  !$acc delete(zkw0, zkwm)                                          &
  !$acc delete(zedb, ztrang, ztrangs, zwin, zwsnow, zwsnew, zdwsnm) &
  !$acc delete(zaa, zw_fr, zinfil, zlw_fr, ziw_fr, zwsnn, zflmg)    &
  !$acc delete(zrunoff_grav, zk0di, zbedi, zsnull, zs1, zf_rad)     &
  !$acc delete(ztraleav, zwroot, zropartw, zts, ztsk, ztsnow)       &
  !$acc delete(ztsnow_mult, zalamtmp, zalam, zrocg, zrocg_soil)     &
  !$acc delete(zrocs, ztsn, ztskn, ztsnown, ztsnown_mult, znlgw1f)  &
  !$acc delete(zqbase, zrefr, zmelt, ze_out, zrho_dry_old)          &
  !$acc delete(zp, zcounter, ze_rad, zswitch, tmp_num, sum_weight)  &
  !$acc delete(t_new, rho_new, wl_new, dz_old, z_old)               &
  !$acc delete(t_so_free_new, t_so_snow_new, sn_frac)               &
  !$acc delete(zf_snow_lim, zdz_snow, zalas_mult, ztsnownew_mult)   &
  !$acc delete(zextinct, zfor_snow_mult, hzalam, zdqvtsnow)         &
  !$acc delete(zrho_snow, zts_pm, ztsk_pm, ztfunc, ztsnow_pm, zeisa)&
  !$acc delete(zw_snow_old, zrho_snow_old, h_snow_fg, h_snow_incr)  &
  !$acc delete(zaga, zagb, zagc, zagd, zage)                        &
  !$acc delete(limit_tch, lzurban, l_redist)

  DEALLOCATE (                                                                  &
            m_styp         , zicount1       , zicount2       ,                  &
            ke_soil_hy_b   ,                                                    &
       STAT=istat)


  DEALLOCATE (                                                                  &
            h_snow_now     , h_snow_new     , zdz_snow_fl    , zhh_snow       , &
            zhm_snow       , zdzh_snow      , zdzm_snow      , zdtsnowdt_mult , &
            zbwt           , zrock          , zsandf         , zclayf         , &
            zsiltf         , zb_por         , zpsis          , zw_m_org       , &
            zw_m_soil      , zw_m_up        , zw_m_low       , zw_m           , &
       STAT=istat)

  DEALLOCATE (                                                                  &
            zrr            , zrs            , zesoil         , zrhoch         , &
            zth_low        , zf_wi          , ztmch          , zep_s          , &
            zep_snow       , zverbo         , zversn         , zthsoi         , &
            zthsnw         , zfor_s         , zgsb           , zrnet_s        , &
            zsprs          , zdwidt         , zdwsndt        , zdtsdt         , &
            zdtsnowdt      , zdwgdt         , zwinstr        , zinfmx         , &
            zwimax         , zvers          , zwisnstr       , zwpnstr        , &
            zfd_wi         , zf_pd          , zwisn          , zwpn           , &
            zewi           , zepd           , zept           , zesn           , &
            zdrr           , zrrs           , zg1            , zlhfl          , &
            zshfl          , zthfl          , zradfl         , ze_melt        , &
            zch_snow       , zeb1           , ztchv          , ztchv_max      , &
            zrho_atm       , zdt_atm        , zdq_atm        , zroota         , &
            zwrootdz       , zrootdz_int    , zwrootdz_int   , zqhfl_s        , &
            zqhfl_snow     , zsfc_frac_bs   ,                                   &
       STAT=istat)

  DEALLOCATE (                                                                  &
            zroc           , zfcap          , zadp           , zporv          , &
            zdlam          , zdw            , zdw1           , zkw            , &
            zkw0           , zkwm           ,                                   &
            zkw1           , zik2           , zpwp           , ztlpmwp        , &
            zedb           , ztrang         , ztrangs        , zwin           , &
            zwsnow         , zwsnew         , zdwsnm         , zw_fr          , &
            zinfil         , zlw_fr         , ziw_fr         , zwsnn          , &
            zflmg          , zrunoff_grav   , zk0di          , zbedi          , &
            zsnull         , zs1            , zf_rad         , ztraleav       , &
            zwroot         , zropartw       , zts            , ztsnow         , &
            ztsnow_mult    , zalamtmp       , zalam          , zrocg          , &
            zrocg_soil     , zrocs          , ztsn           , ztsnown        , &
            ztsnown_mult   , znlgw1f        , zaa            ,                  &
            ztsk           , ztskn          ,                                   &
       STAT=istat)


  DEALLOCATE (                                                                  &
            zqbase         , zrefr          , zmelt          , ze_out         , &
            zrho_dry_old   , zp             , zcounter       , ze_rad         , &
            zswitch        , tmp_num        , sum_weight     , t_new          , &
            rho_new        , wl_new         , z_old          , dz_old         , &
            t_so_free_new  , t_so_snow_new  , sn_frac        , zf_snow_lim    , &
            zdz_snow       , zalas_mult     , ztsnownew_mult , zextinct       , &
            zfor_snow_mult , hzalam         , zdqvtsnow      , zrho_snow      , &
            zts_pm         , ztfunc         , ztsnow_pm      , zeisa          , &
            ztsk_pm        ,                                                    &
       STAT=istat)

  DEALLOCATE (                                                                  &
            ! for terra_init
            zw_snow_old    , zrho_snow_old  , h_snow_fg      , h_snow_incr    , &
       STAT=istat)

  DEALLOCATE (                                                                  &
            zaga           , zagb           , zagc           , zagd           , &
            zage           ,                                                    &
       STAT=istat)

  DEALLOCATE (                                                                  &
            limit_tch      , lzurban        , l_redist       ,                  &
       STAT=istat)

END SUBROUTINE terra_wkarr_dealloc

!==============================================================================
#endif
         !ALLOC_WKARR
!==============================================================================

END MODULE sfc_terra_data
