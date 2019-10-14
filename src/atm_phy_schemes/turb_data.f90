!+ Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------

MODULE turb_data

!------------------------------------------------------------------------------
!
! Description:
!  This module contains parameters that are used in the turbulence
!  parameterizations. With some of these parameters a tuning of the schemes
!  is possible.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  Matthias.Raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V5_4a        2016-05-10 Matthias Raschendorfer, Ulrich Schaettler
!  Initial release, based on data_turbulence with modifications for ICON
!  Note that this module is now also used for the old (not blocked) scheme
! V5_4b        2016-07-12 Ulrich Schaettler
!  Activate old defaults for some namelist variables for COSMO 
!  as long as non-blocked turbulence scheme is running per default
! V5_4c        2016-10-06 Ulrich Schaettler
!  Again use local memory if not running on GPUs 
!     (because on vectorization problem on CRAY)
! V5_4d        2016-12-12 Matthias Raschendorfer
!  Adapting some comments and reshuffling variables.
! V5_4e        2017-03-23 Ulrich Schaettler
!  Modified several Namelist defaults to the ICON defaults
! V5_4f        2017-09-01 Matthias Raschendorfer, Ulrich Schaettler
!  Unified settings for ICON and COSMO:
!   - imode_rat_sea:    1  (ICON value,  changes in COSMO)
!   - loutshs        TRUE  (COSMO value, changes in ICON)
!  To get an old-COSMO-like behaviour of the blocked turbulence scheme some 
!   hardcoded variables need different values. These are set by compiling 
!   with the pragma -DCOSMO_OLD
!  Removed wichfakt, securi from the old version (US)
! V5_4g        2017-11-13 Ulrich Schaettler
!  Implemented namelist variable loldtur (Default .FALSE.) and removed 
!  pragma -DCOSMO_OLD to be able to switch the behaviour during run time.
! V5_4h        2017-12-15 Xavier Lapillonne
!  Modifications to port turbulence scheme to GPU
! V5_5         2018-02-23 Ulrich Schaettler
!  Updated with ICON Version 7bcba73: 
!   - new (still internal) switch imode_tkesso (should replace ltkesso)
!   - modified some default values of namelist variables with ifdef ICON/COSMO
! V5_6         2019-02-27 Ulrich Schaettler
!  Updated with ICON Version d7e0252
!    (Set alpha1 = 0.75; was 1.0 before)
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
USE kind_parameters, ONLY : &
    wp           ! KIND-type parameter for real variables
#endif

#ifdef __ICON__
USE mo_kind,                ONLY: wp           ! KIND-type parameter for real variables

USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

USE mo_turbdiff_config,     ONLY: turbdiff_config
#endif

!==============================================================================

IMPLICIT NONE

PUBLIC

!==============================================================================
! Configuration parameters:
! ----------------------------------------------------------------------------

INTEGER, PARAMETER :: &
!
    mom=1,        & ! index for a momentum variable
    sca=2           ! index for a scalar   variable

INTEGER, PARAMETER :: &
!
!  Indexgrenzen:
!
   nscal=3,     & !aktive skalare Groessen 1-ter Ordnung ('tem', 'vap', 'liq')
   ninv=2,      & !daraus abgeleitete gegenueber vertikalen (feuchtadiabatischen)
                  !Verrueckungen invarianten Groessen ('tet_l', 'h2o_g')
   nvel=2,      & !aktive Geschwindigkeitskomponenten ('u_m', 'v_m')
   nmvar=nscal+nvel, &
   nred=ninv+nvel,   &
   naux=5,      & !number of auxilary variables
   ntyp=2,      & !Anzahl von Variablentypen (mom) und (sca)
   ntmax=3,     & !max. Anzahl der Zeitebenen fuer die TKE
!
!    Zeiger fuer die Variablen :
!
   u_m=1,       & !zonale Geschw.komp. im Massenzentrum
   v_m=2,       & !meridionale  ,,      ,,     ,,
   tet_l=3,     & !feucht-potentielle Temperatur
   tem_l=tet_l, & !Fluessigwasser-Temperatur
   h2o_g=4,     & !Gesamtwasseergehalt
   liq=5,       & !Fluessigwasser  ,,
   w_m=6,       & !vertikale Geschw.komp. im Massenzentrum
!
   tet=tet_l,   & !pot.Temperatur
   tem=tet,     & !Temperatur
   vap=h2o_g,   & !Wasserdampfmischungsverh.

   ndim=MAX(nmvar,naux)

!  Beachte:     u_m,v_m muessen in [1,nvel] liegen;
!         aber: tet_l,h2o_g in [nvel+1,nred]
!         und   tem (tet),vap,liq in [nvel+1,nmvar]


!==============================================================================
! Parameters that may be used for tuning and special configurations:
! ----------------------------------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs!

! 1. Numerical parameters:
!-----------------------------

REAL (KIND=wp)     ::        &
  impl_s       =  1.20_wp,   & ! implicit weight near the surface (maximal value)
  impl_t       =  0.75_wp,   & ! implicit weight near top of the atmosphere (maximal value)

  ! Minimal diffusion coefficients in [m^2/s] for vertical
  tkhmin       =  0.75_wp,   & ! scalar (heat) transport
  tkmmin       =  0.75_wp,   & ! momentum transport
#ifdef __COSMO__
  tkhmin_strat =  5.00_wp,   & ! scalar (heat) transport, enhanced value for stratosphere
  tkmmin_strat =  5.00_wp,   & ! momentum transport,      enhanced value for stratosphere
#endif
#ifdef __ICON__
  tkhmin_strat =  0.75_wp,   & ! scalar (heat) transport, enhanced value for stratosphere
  tkmmin_strat =  4.00_wp,   & ! momentum transport,      enhanced value for stratosphere
#endif

  ditsmot      =  0.00_wp,   & ! smoothing factor for direct time-step iteration

  tndsmot      =  0.00_wp,   & ! vertical smoothing factor for diffusion tendencies
  frcsmot      =  0.00_wp,   & ! vertical smoothing factor for TKE forcing (in ICON only in the tropics)
  tkesmot      =  0.15_wp,   & ! time smoothing factor for TKE and diffusion coefficients
  stbsmot      =  0.00_wp,   & ! time smoothing factor for stability function
  frcsecu      =  1.00_wp,   & ! security factor for TKE-forcing       (<=1)
  tkesecu      =  1.00_wp,   & ! security factor in  TKE equation      (out of [0; 1])
  stbsecu      =  0.00_wp,   & ! security factor in stability function (out of [0; 1])

  epsi         =  1.0E-6_wp    ! relative limit of accuracy for comparison of numbers

INTEGER            ::        &

  it_end       =  1            ! number of initialization iterations (>=0)

! 2. Parameters describing physical properties of the lower boundary 
!    of the atmosphere:
!------------------------------------------

REAL (KIND=wp)     ::        &
  rlam_mom     =  0.0_wp,    & ! scaling factor of the laminar boundary layer for momentum
  rlam_heat    =  1.0_wp,    & ! scaling factor of the laminar boundary layer for heat

  rat_lam      =  1.0_wp,    & ! ratio of laminar scaling factors for vapour and heat
  rat_can      =  1.0_wp,    & ! ratio of canopy height over z0m
  rat_sea      = 10.0_wp,    & ! ratio of laminar scaling factors for heat over sea and land

  z0m_dia      =  0.2_wp,    & ! roughness length of a typical synoptic station [m]

  alpha0       =  0.0123_wp, & ! Charnock-parameter
  alpha0_max   =  0.0335_wp, & ! upper limit of velocity-dependent Charnock-parameter
  alpha0_pert  =  0.0_wp,    & ! additive ensemble perturbation of Charnock-parameter

#ifdef __ICON__
  alpha1       =  0.7500_wp    ! parameter scaling the molecular roughness of water waves
#endif

#ifdef __COSMO__
  alpha1       =  1.0000_wp    ! parameter scaling the molecular roughness of water waves
#endif


!$acc declare copyin(alpha0,alpha0_max,alpha0_pert)

! 3. Parameters that should be external parameter fields being not yet 
!    available:
!------------------------------------------

REAL (KIND=wp)     ::        &
  c_lnd        = 2.0_wp,     & ! surface area density of the roughness elements over land
  c_sea        = 1.5_wp,     & ! surface area density of the waves over sea
  c_soil       = 1.0_wp,     & ! surface area density of the (evaporative) soil surface
  e_surf       = 1.0_wp        ! exponent to get the effective surface area


! 4. Parameters that should be dynamical fields being not yet available:
!------------------------------------------

REAL (KIND=wp)     ::        &
  z0_ice       =  0.001_wp     !roughness length of sea ice


! 5. Parameters for modelling turbulent diffusion:
!------------------------------------------

REAL (KIND=wp)     ::        &
  tur_len      = 500.0_wp,   & ! asymptotic maximal turbulent distance [m]
  pat_len      = 100.0_wp,   & ! effective length scale of subscale surface patterns over land [m]
                               ! (should be dependent on location)
  len_min      =  1.0E-6_wp, & ! minimal turbulent length scale [m]

  vel_min      =  0.01_wp,   & ! minimal velocity scale [m/s]

  akt          =  0.4_wp,    & ! von Karman-constant

  ! Length scale factors for pressure destruction of turbulent
  a_heat       =  0.74_wp,   & ! scalar (heat) transport
  a_mom        =  0.92_wp,   & ! momentum transport

  ! Length scale factors for dissipation of
  d_heat       =  10.1_wp,   & ! scalar (temperature) variance
  d_mom        =  16.6_wp,   & ! momentum variance

  ! Length scale factors for turbulent transport (vertical diffusion)
  c_diff       =  0.20_wp,   & ! of TKE

  ! Length scale factor for separate horizontal shear production
  a_hshr       =  1.00_wp,   & ! of TKE

  ! Length scale factor for the stability correction
  a_stab       =  0.00_wp,   & ! no stability correction so far

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):
  clc_diag     =  0.5_wp,    & !cloud cover at saturation
#ifdef __COSMO__
  q_crit       =  4.0_wp,    & !critical value for normalized over-saturation
#endif
#ifdef __ICON__
  q_crit       =  1.6_wp,    & !critical value for normalized over-saturation
#endif
  c_scld       =  1.0_wp       !factor for liquid water flux density in sub grid scale clouds

!==============================================================================
! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

LOGICAL :: &

  loldtur       =.FALSE., & ! use settings to simulate old ijk turbulence version
                            ! if .TRUE.: new ICON-like settings are used
  ltkesso       =.TRUE.,  & ! calculation SSO-wake turbulence production for TKE
  ltkecon       =.FALSE., & ! consider convective buoyancy production for TKE
  ltkeshs       =.TRUE. , & ! consider separ. horiz. shear production for TKE
  loutshs       =.TRUE. , & ! consider separ. horiz. shear production of TKE for output
  lnonloc       =.FALSE., & ! nonlocal calculation of vertical gradients used for turbul. diff.
  lprfcor       =.FALSE., & ! using the profile values of the lowest main level instead of
                            ! the mean value of the lowest layer for surface flux calulations

  ltmpcor       =.FALSE., & ! consideration of thermal TKE-sources in the enthalpy budget
  lcpfluc       =.FALSE., & ! consideration of fluctuations of the heat capacity of air

  lexpcor       =.FALSE., & ! explicit corrections of the implicit calculated turbul. diff.

! for semi-implicit vertical diffusion: 
  lsflcnd       =.TRUE. , & ! lower flux condition for vertical diffusion calculation
  ldynimp       =.FALSE., & ! dynamical calculation of implicit weights
  lprecnd       =.FALSE., & ! preconditioning of tridiagonal matrix
  lfreeslip     =.FALSE.    ! free-slip lower boundary condition (use for idealized runs only!)

! Notice that the following switches are provided by the parameter-list of 
! SUB 'turb_diffusion' or 'turb_transfer':

! lstfnct                   :calculation of stability function required
! lnsfdia                   :calculation of (synoptical) near-surface variables required
! lmomdif                   :calculation of complete gradient diffusion of horizontal momenum
! lscadif                   :calculation of complete gradient diffusion of scalar properties
! lturatm                   :running turbulence model between atmosph. layers (updating diffusion coefficients)
! ltursrf                   :running turbulence model at the surface layer (updating transfer coefficients
! lsfluse                   :use explicit heat flux densities at the suface
! ltkeinp                   :TKE present as input (at level k=ke1 for current time level 'ntur')
! lgz0inp                   :gz0 present as input

!==============================================================================
! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

INTEGER :: &

  imode_tran    =0, & ! mode of TKE-equation in transfer scheme             (compare 'imode_turb')
  imode_turb    =1, & ! mode of TKE-equation in turbulence scheme
                      !  0: diagnostic equation
                      !  1: prognostic equation (default)
                      !  2: prognostic equation (implicitly positive definit)
  icldm_tran    =2, & ! mode of cloud representation in transfer parametr.  (compare 'icldm_turb')
  icldm_turb    =2, & ! mode of cloud representation in turbulence parametr.
                      ! -1: ignoring cloud water completely (pure dry scheme)
                      !  0: no clouds considered (all cloud water is evaporated)
                      !  1: only grid scale condensation possible
                      !  2: also sub grid (turbulent) condensation considered
  itype_wcld    =2, & ! type of water cloud diagnosis within the turbulence scheme:
                      ! 1: employing a scheme based on relative humitidy
                      ! 2: employing a statistical saturation adjustment
  itype_sher    =0    ! type of shear production for TKE
                      ! 0: only vertical shear of horizontal wind
                      ! 1: previous plus horizontal shear correction
                      ! 2: previous plus shear from vertical velocity
#ifdef __ICON__
INTEGER :: &
  itype_diag_t2m=1    ! type of diagnostics of 2m-temperature and -dewpoint
                      ! 1: Considering a fictive surface roughness of a SYNOP lawn
                      ! 2: Considering the mean surface roughness of a grid box
                      !    and using an exponential roughness layer profile
#endif

! To reproduce the old ijk turbulence settings as good as possible, all these switches 
! have to be set to 1. If loldtur=.TRUE., the re-setting is done in organize_physics.

! These are the settings for the ICON-like setup of the physics
INTEGER :: &

  imode_stbcorr =1, & ! mode of correcting the stability function (related to 'stbsecu')
                      ! 1: always for strict.-non-stb. strat. using a restr. 'gama' in terms of prev. forc.
                      ! 2: only to avoid non-physic. solution or if current 'gama' is too large
  imode_tkemini =1, & ! mode of fixing a lower limit of q=2TKE**2
                      ! 1: by using 'vel_min' only
                      ! 2: by adapting to minimal diffusion coefficients
  imode_lamdiff =1, & ! mode of considering laminar diffusion at surface layer
                      ! 1: only when calculating the profile functions
                      ! 2: surface-layer diffusion coeff. always at least at laminar value
  ilow_def_cond =2, & ! type of the default condition at the lower boundary
                      ! 1: zero surface flux density
                      ! 2: zero surface value
  imode_calcirc =2, & ! mode of treating the circulation term (related to 'pat_len', imode_pat_len')
                      ! 1: explicit calculation of the flux convergence
                      ! 2: quasi implicit treatment by calculation of effective TKE-gradients
  imode_pat_len =2, & ! mode of determining the length scale of surface patterns (related to 'pat_len')
                      ! 1: by the constant value 'pat_len' only
                      ! 2: and the std. deviat. of SGS orography as a lower limit (only if 'd_pat' is pres.)
  imode_frcsmot =2, & ! if "frcsmot>0", apply smoothing of TKE source terms 
                      ! 1: globally or 
                      ! 2: in the tropics only (if 'trop_mask' is present) 
  imode_shshear =2, & ! mode of calculat. the separated horizontal shear mode (related to 'ltkeshs', 'a_hshr')
                      ! 0: with a constant lenght scale and based on 3D-shear and incompressibility
                      ! 1: with a constant lenght scale and considering the trace constraint for the 2D-strain tensor
                      ! 2: with a Ri-depend. length sclale correct. and the trace constraint for the 2D-strain tensor
  imode_tkesso =1, &  ! mode of calculat. the SSO source term for TKE production
                      ! 1: original implementation
                      ! 2: with a Ri-dependent reduction factor for Ri>1
  imode_tkvmini =2, & ! mode of calculating the minimal turbulent diff. coeffecients
                      ! 1: with a constant value
                      ! 2: with a stability dependent correction
  imode_rat_sea =1, & ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                      ! 1: constant ratio 'rat_sea' to land 
                      ! 2: with a correction for a strongly overheated SST
  imode_vel_min =2, & ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                      ! 1: with a constant value
                      ! 2: with a stability dependent correction
  imode_charpar =2    ! mode of estimating the Charnock-Parameter
                      ! 1: use a constant value 
                      ! 2: use a wind-dependent value with a constant lower bound

INTEGER :: &

  imode_syndiag =2, & ! mode of diagnostics at the synoptic near surface levels (related to 'itype_diag_t2m')
                      ! 1: direct interpolation of temperature and specific humidity
                      ! 2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                      !    allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
  imode_qvsatur =2, & ! mode of calculating the saturat. humidity
                      ! 1: old version using total pressure
                      ! 2: new version using partial pressure of dry air
  imode_stadlim =2, & ! mode of limitting statist. saturation adjustment (SUB 'turb_cloud')
                      ! 1: only absolut upper limit of stand. dev. of local oversatur. (sdsd)
                      ! 2: relative limit of sdsd and upper limit of cloud-water
  imode_trancnf =2, & ! mode of configuring the transfer-scheme (SUB 'turbtran')
                      ! 1: old version: start. with lamin. diffus.; with a lamin. correct. for profile-funct.;
                      !    interpol. T_s rather then Tet_l onto zero-level; calcul. only approx. Tet_l-grads.;
                      !    using an upper bound for TKE-forcing; without transmit. skin-layer depth to turbul.
                      ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                      !    interpol. Tet_l onto zero-level; calcul. Tet_l-gradients directly; 
                      !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                      ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                      !    for stable stratification
                      ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
  imode_tkediff =2, & ! mode of implicit TKE-Diffusion (related to 'c_diff')
                      ! 1: in terms of q=SQRT(2*TKE)) 
                      ! 2; in terms of TKE=0.5*TKE**2
  imode_adshear =2    ! mode of considering addit. shear by scale interaction (realt. to 'ltkesso', 'ltkeshs',
                      ! 'ltkecon')
                      ! 1: not consid. for stability functions
                      ! 2:  considered for stability functions

! Notice that the following selectors are provided by the parameter-list of 
! SUB 'turb_diffusion' or 'turb_transfer':

! iini                :type of initialization (0: no, 1: separate before the time loop
!                                                   , 2: within the first time step)
! itnd                :type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
!                                                     2: by adding to current profile before vertical diffusion
!                                                     3: by using corrected virtual vertical profiles

!==============================================================================
! Declarations of utility variables:

! Turbulence parameters which are computed during model run
!-------------------------------------------------------------------------------

REAL (KIND=wp)     ::        &
  ! do we need it as TARGET?
  ! these variables are set in SR turb_param
  c_tke,tet_g,rim, &
  c_m,c_h, b_m,b_h,  sm_0, sh_0, &
  d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
  a_3,a_5,a_6,                 &

  ! these parameters are physical constants which are either taken as
  ! they are or set to 0.0 in the turbulence for special applications
  tur_rcpv,         & ! cp_v/cp_d - 1
  tur_rcpl            ! cp_l/cp_d - 1 (where cp_l=cv_l)

! Definition of used data types
!-------------------------------

TYPE modvar !model variable
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
#else
     REAL (KIND=wp), POINTER                 ::         &
#endif
             av(:,:) => NULL(), & !atmospheric values
             sv(:)   => NULL(), & !surface     values (concentration of flux density)
             at(:,:) => NULL()    !atmospheric time tendencies
     LOGICAL                                 ::         &
             fc                !surface values are flux densities
END TYPE modvar

TYPE turvar !turbulence variables
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
#else
     REAL (KIND=wp), POINTER                 ::         &
#endif
             tkv(:,:) => NULL(), & !turbulent coefficient for vert. diff.
             tsv(:)   => NULL()    !turbulent velocity at the surface
END TYPE turvar

TYPE varprf !variable profile
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
#else
     REAL (KIND=wp), POINTER                 ::         &
#endif
             bl(:,:), & !variable at boundary model levels
             ml(:,:)    !variable at main     model levels
END TYPE varprf

#ifdef ALLOC_WKARR
! Turbulence parameters which are computed during model run
! ---------------------------------------------------------

! For any application of the blocked default turbulence code:
REAL (KIND=wp), ALLOCATABLE, TARGET :: &
     eprs(:,:),    &  ! surface Exner-factor

     len_scale(:,:)  ,& ! turbulent length-scale (m)

     rhon(:,:)       ,& ! bondary level air density (including surface level) (Kg/m3)
     frh(:,:)        ,& ! thermal forcing (1/s2) or thermal acceleration (m/s2)
     frm(:,:)        ,& ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)
     zaux(:,:,:)     ,& ! auxilary array containing thermodynamical properties
                        ! (dQs/dT,ex_fakt,cp_fakt,g_tet,g_vap) or various 
                        ! auxilary variables for calculation of implicit vertical diffusion
     zvari(:,:,:)       ! set of variables used in the turbulent 2-nd order equations
                        ! and later their effective vertical gradients

! For atmospheric turbulence or surface transfer only:

REAL (KIND=wp), ALLOCATABLE, TARGET :: &
     l_scal(:),    &  ! reduced maximal turbulent length scale due to horizontal grid spacing (m)
     fc_min(:),    &  ! minimal value for TKE-forcing (1/s2)

     grad  (:,:),  &  ! any vertical gradient
     hig   (:,:),  &  ! obere und untere Referenzhoehe bei der Bildung nicht-lokaler Gradienten

     prss(:,:),    &  ! surface pressurei (Pa)
     tmps(:,:),    &  ! surface temperature (K)
     vaps(:,:),    &  ! surface specific humidity
     liqs(:,:),    &  ! liquid water content at the surface

     diss_tar (:,:)   ! target for eddy dissipation rate (m2/s3)

! For atmpospheric turbulence or vertical diffusion (turbdiff) only:

REAL (KIND=wp), ALLOCATABLE, TARGET :: &
     hlp  (:,:)       ,& ! any 'help' variable
     dicke(:,:)       ,& ! any (effective) depth of model layers (m) or other auxilary variables

! For atmospheric turbulence only:

     can(:,:)         ,& ! auxilary array valid for the vertically resolved canopy

     ftm(:,:)         ,& ! mechan. forcing (1/s2) by pure turbulent shear 
     shv(:,:)         ,& ! velocity scale of the separated horiz. shear mode (m/s)

     hor_scale(:,:)   ,& ! effective hoprizontal length-scale used for sep. horiz. shear calc. (m)
     xri(:,:)         ,& ! a function of Ri-number used for tuning corrections (hyper-parameterizations) 


     lay(:)           ,& ! any variable at a specific layer
     lays(:,:)        ,& ! any (2-D) vector of variables at a specific layer

     src(:)           ,& ! arbitrary source term

     dzsm(:)          ,& ! effective depth of Prandtl-layer applied to momentum (m)
     dzsh(:)             ! effective depth of Prandtl-layer applied to scalars  (m)

! For surface transfer (turbtran) only:

REAL (KIND=wp), ALLOCATABLE, TARGET :: &
     rclc   (:,:),    & ! cloud cover

 tketens_tar(:,:)       ! target for turbulent transport of SQRT(TKE)

LOGICAL, ALLOCATABLE ::  &
     lo_ice   (:)       ! logical sea ice indicator

INTEGER, ALLOCATABLE ::  &
     lev      (:,:),  & ! eingrenzende Hoehenvieaus
     k_2d     (:)       ! index field of the upper level index to be used for near surface diagn.

REAL (KIND=wp), ALLOCATABLE, TARGET ::  &
     hk_2d    (:),    & ! mid level height above ground belonging to 'k_2d' (m)
     hk1_2d   (:),    & ! mid level height above ground of the previous layer (below) (m)

     h_top_2d (:),    & ! boundary level height of transfer layer (top  level) (m)
     h_atm_2d (:),    & ! mid      level heigth of transfer layer (atm. level) (m)
     h_can_2d (:),    & ! effective canopy height (m) used for calculation of roughness layer
                        ! resistance for momentum

     edh      (:),    & ! reciprocal of a layer depth (1/m)

     z0m_2d   (:),    & ! mean  roughness length (m)
     z0d_2d   (:),    & ! diag. roughness length (m)
     z2m_2d   (:),    & ! effective height of 2m  level (above the surface) (m)
     z10m_2d  (:)       ! effective height of 10m level (above the surface) (m)

REAL (KIND=wp), ALLOCATABLE, TARGET ::  &
     rat_m_2d (:),    & ! any surface layer ratio for momentum (like Re-number or Stability factor)
     rat_h_2d (:),    & ! any surface layer ratio for scalars  (like Re-number or Stability factor)
     fac_m_2d (:),    & ! surface layer profile factor for momentum
     fac_h_2d (:),    & ! surface layer profile factor for scalars

     frc_2d   (:),    & ! length scale fraction related to bottom and top of the transfer layer

     vel_2d   (:),    & ! wind speed (m/s) at the top of the transfer layer (lowest mid level of atm. model)
     tl_s_2d  (:),    & ! surface level value of conserved temperature (liquid water temperature) (K)
     qt_s_2d  (:),    & ! surface level value of conserved humidity    (total water)
     velmin   (:),    & ! modified 'vel_min' used for tuning corrections (hyper-parameterizations) (m/s)
     ratsea   (:),    & ! modified 'rat_sea' used for tuning corrections (hyper-parameterizations)

     dz_sg_m  (:),    & ! laminar resistance lenght for momentum (m)
     dz_sg_h  (:),    & ! laminar resistance length for scalars (m)
     dz_g0_m  (:),    & ! turbulent roughness layer resistance length for momentum (m)
     dz_g0_h  (:),    & ! turbulent roughness layer resistance length for scalars (m)
     dz_0a_m  (:),    & ! turbulent Prandtl-layer resistance length for momentum (m)
     dz_0a_h  (:),    & ! turbulent Prandtl-layer resistance length for scalars (m)
     dz_s0_h  (:),    & ! total roughness layer resistance lenght for scalars (m)
     dz_sa_h  (:)       ! total transfer resistance length for scalars (m)
#endif


#ifdef  __ICON__
! 7. Switches from the COSMO-Model: must only be defined for ICON
! ---------------------------------------------------------------

REAL(KIND=wp), POINTER :: &
  impl_weight(:)  ! implicit weights for tridiagonal solver

! Switches controlling turbulent diffusion:
! ------------------------------------------

INTEGER  :: &

    itype_tran   =2,       & ! type of surface-atmosphere transfer
    imode_circ   =2,       & ! mode of treating the circulation term
    ilow_dcond   =1          ! type of the default condition at the lower boundary
#endif

!==============================================================================

CONTAINS

!==============================================================================

! MR: Uli, the following SUB 'get_turbdiff' should be ported to the ICON-interface 
!     for TURBDIFF!!
#ifdef  __ICON__

! this subroutines sets the switches defined above during an ICON run

SUBROUTINE get_turbdiff_param (jg)

   INTEGER, INTENT(IN) :: jg !patch index

   impl_weight => turbdiff_config(jg)%impl_weight

   imode_tran   = turbdiff_config(jg)%imode_tran
   icldm_tran   = turbdiff_config(jg)%icldm_tran
   imode_turb   = turbdiff_config(jg)%imode_turb
   icldm_turb   = turbdiff_config(jg)%icldm_turb
   itype_sher   = turbdiff_config(jg)%itype_sher
   imode_frcsmot= turbdiff_config(jg)%imode_frcsmot
   imode_tkesso = turbdiff_config(jg)%imode_tkesso

   ltkesso      = turbdiff_config(jg)%ltkesso
   ltkecon      = turbdiff_config(jg)%ltkecon
   ltkeshs      = turbdiff_config(jg)%ltkeshs
   lexpcor      = turbdiff_config(jg)%lexpcor
   ltmpcor      = turbdiff_config(jg)%ltmpcor
   lprfcor      = turbdiff_config(jg)%lprfcor
   lnonloc      = turbdiff_config(jg)%lnonloc
   lfreeslip    = turbdiff_config(jg)%lfreeslip
   lcpfluc      = turbdiff_config(jg)%lcpfluc
   lsflcnd      = turbdiff_config(jg)%lsflcnd

   itype_wcld   = turbdiff_config(jg)%itype_wcld

   tur_len      = turbdiff_config(jg)%tur_len
   pat_len      = turbdiff_config(jg)%pat_len
   a_stab       = turbdiff_config(jg)%a_stab
   tkhmin       = turbdiff_config(jg)%tkhmin
   tkmmin       = turbdiff_config(jg)%tkmmin
   tkhmin_strat = turbdiff_config(jg)%tkhmin_strat
   tkmmin_strat = turbdiff_config(jg)%tkmmin_strat

   alpha0       = turbdiff_config(jg)%alpha0
   alpha0_max   = turbdiff_config(jg)%alpha0_max
   alpha0_pert  = turbdiff_config(jg)%alpha0_pert

   c_diff       = turbdiff_config(jg)%c_diff
   rlam_heat    = turbdiff_config(jg)%rlam_heat
   rlam_mom     = turbdiff_config(jg)%rlam_mom
   rat_sea      = turbdiff_config(jg)%rat_sea
   tkesmot      = turbdiff_config(jg)%tkesmot
   frcsmot      = turbdiff_config(jg)%frcsmot
   impl_s       = turbdiff_config(jg)%impl_s
   impl_t       = turbdiff_config(jg)%impl_t
   q_crit       = turbdiff_config(jg)%q_crit

   loutshs      = ltkeshs .OR. itype_sher > 0

END SUBROUTINE get_turbdiff_param
#endif

!==============================================================================

#ifdef ALLOC_WKARR
SUBROUTINE turb_wkarr_alloc (ke, kcm, nproma, istat)

  INTEGER, INTENT(IN)  :: ke, kcm, nproma
  INTEGER, INTENT(OUT) :: istat

  INTEGER :: kez1

  kez1 = ke+1

! For any application of the blocked default turbulence code:

  ALLOCATE ( eprs         (nproma,kez1:kez1)         , STAT=istat)

  ALLOCATE ( len_scale    (nproma,kez1)              , STAT=istat)

  ALLOCATE ( rhon         (nproma,kez1)              , STAT=istat)
  ALLOCATE ( frh          (nproma,kez1)              , STAT=istat)
  ALLOCATE ( frm          (nproma,kez1)              , STAT=istat)
  ALLOCATE ( zaux         (nproma,kez1,ndim)         , STAT=istat)
  ALLOCATE ( zvari        (nproma,kez1,ndim)         , STAT=istat)

! For turbdiff or turbtran only:

  ALLOCATE ( l_scal       (nproma)                   , STAT=istat)
  ALLOCATE ( fc_min       (nproma)                   , STAT=istat)

  ALLOCATE ( grad         (nproma,nmvar)             , STAT=istat)
  ALLOCATE ( hig          (nproma,2)                 , STAT=istat)

  ALLOCATE ( tmps         (nproma,kez1:kez1)         , STAT=istat)
  ALLOCATE ( vaps         (nproma,kez1:kez1)         , STAT=istat)
  ALLOCATE ( prss         (nproma,kez1:kez1)         , STAT=istat)
  ALLOCATE ( liqs         (nproma,kez1:kez1)         , STAT=istat)

  ALLOCATE ( diss_tar     (nproma,kez1)              , STAT=istat)

! For turbdiff or vertical diffusion:

  ALLOCATE ( hlp          (nproma,kez1)              , STAT=istat)
  ALLOCATE ( dicke        (nproma,kez1)              , STAT=istat)


! For turbdiff only:

  ALLOCATE ( can          (nproma,kcm:kez1)          , STAT=istat)

  ALLOCATE ( ftm          (nproma,kez1)              , STAT=istat)
  ALLOCATE ( shv          (nproma,kez1)              , STAT=istat)

  ALLOCATE ( hor_scale    (nproma,ke)                , STAT=istat)
  ALLOCATE ( xri          (nproma,ke)                , STAT=istat)

  ALLOCATE ( lay          (nproma)                   , STAT=istat)
  ALLOCATE ( lays         (nproma,2)                 , STAT=istat)

  ALLOCATE ( src          (nproma)                   , STAT=istat)

  ALLOCATE ( dzsm         (nproma)                   , STAT=istat)
  ALLOCATE ( dzsh         (nproma)                   , STAT=istat)

  ALLOCATE ( lev          (nproma,2)                 , STAT=istat)

! For turbtran only:

  ALLOCATE ( rclc         (nproma, kez1:kez1)        , STAT=istat)
  ALLOCATE ( tketens_tar  (nproma, kez1:kez1)        , STAT=istat)

  ALLOCATE ( lo_ice       (nproma)                   , STAT=istat)

  ALLOCATE ( k_2d         (nproma)                   , STAT=istat)

  ALLOCATE ( hk_2d        (nproma)                   , STAT=istat)
  ALLOCATE ( hk1_2d       (nproma)                   , STAT=istat)

  ALLOCATE ( h_top_2d     (nproma)                   , STAT=istat)
  ALLOCATE ( h_atm_2d     (nproma)                   , STAT=istat)
  ALLOCATE ( h_can_2d     (nproma)                   , STAT=istat)

  ALLOCATE ( edh          (nproma)                   , STAT=istat)

  ALLOCATE ( z0m_2d       (nproma)                   , STAT=istat)
  ALLOCATE ( z0d_2d       (nproma)                   , STAT=istat)
  ALLOCATE ( z2m_2d       (nproma)                   , STAT=istat)
  ALLOCATE ( z10m_2d      (nproma)                   , STAT=istat)

  ALLOCATE ( rat_m_2d     (nproma)                   , STAT=istat)
  ALLOCATE ( rat_h_2d     (nproma)                   , STAT=istat)
  ALLOCATE ( fac_h_2d     (nproma)                   , STAT=istat)
  ALLOCATE ( fac_m_2d     (nproma)                   , STAT=istat)

  ALLOCATE ( frc_2d       (nproma)                   , STAT=istat)

  ALLOCATE ( vel_2d       (nproma)                   , STAT=istat)
  ALLOCATE ( tl_s_2d      (nproma)                   , STAT=istat)
  ALLOCATE ( qt_s_2d      (nproma)                   , STAT=istat)
  ALLOCATE ( velmin       (nproma)                   , STAT=istat)
  ALLOCATE ( ratsea       (nproma)                   , STAT=istat)

  ALLOCATE ( dz_sg_m      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_sg_h      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_g0_m      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_g0_h      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_0a_m      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_0a_h      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_sa_h      (nproma)                   , STAT=istat)
  ALLOCATE ( dz_s0_h      (nproma)                   , STAT=istat)

  !$acc enter data create(eprs,len_scale,rhon,frh,frm,zaux,zvari,l_scal)
  !$acc enter data create(fc_min,grad,hig,tmps,vaps,prss,liqs,diss_tar,hlp)
  !$acc enter data create(dicke,can,ftm,shv)
  !$acc enter data create(hor_scale,xri,lay,lays,src,dzsm,dzsh,lev,rclc)
  !$acc enter data create(tketens_tar,lo_ice,k_2d,hk_2d,hk1_2d,h_top_2d)
  !$acc enter data create(h_atm_2d,h_can_2d,edh,z0m_2d,z0d_2d,z2m_2d)
  !$acc enter data create(z10m_2d,rat_m_2d,rat_h_2d,fac_h_2d,fac_m_2d,frc_2d)
  !$acc enter data create(vel_2d,tl_s_2d,qt_s_2d,velmin,ratsea,dz_sg_m,dz_sg_h)
  !$acc enter data create(dz_g0_m,dz_g0_h,dz_0a_m,dz_0a_h,dz_sa_h,dz_s0_h)

END SUBROUTINE turb_wkarr_alloc

!==============================================================================
!==============================================================================

SUBROUTINE turb_wkarr_dealloc (istat)

  INTEGER, INTENT(OUT) :: istat

! For any application of the blocked default turbulence code:
  !$acc exit data delete(eprs,len_scale,rhon,frh,frm,zaux,zvari,l_scal)
  !$acc exit data delete(fc_min,grad,hig,tmps,vaps,prss,liqs,diss_tar,hlp)
  !$acc exit data delete(dicke,can,ftm,shv)
  !$acc exit data delete(hor_scale,xri,lay,lays,src,dzsm,dzsh,lev,rclc)
  !$acc exit data delete(tketens_tar,lo_ice,k_2d,hk_2d,hk1_2d,h_top_2d)
  !$acc exit data delete(h_atm_2d,h_can_2d,edh,z0m_2d,z0d_2d,z2m_2d)
  !$acc exit data delete(z10m_2d,rat_m_2d,rat_h_2d,fac_h_2d,fac_m_2d,frc_2d)
  !$acc exit data delete(vel_2d,tl_s_2d,qt_s_2d,velmin,ratsea,dz_sg_m,dz_sg_h)
  !$acc exit data delete(dz_g0_m,dz_g0_h,dz_0a_m,dz_0a_h,dz_sa_h,dz_s0_h)

  DEALLOCATE ( eprs           , STAT=istat)

  DEALLOCATE ( len_scale      , STAT=istat)

  DEALLOCATE ( rhon           , STAT=istat)
  DEALLOCATE ( frh            , STAT=istat)
  DEALLOCATE ( frm            , STAT=istat)
  DEALLOCATE ( zaux           , STAT=istat)
  DEALLOCATE ( zvari          , STAT=istat)

! For turbdiff or turbtran only:

  DEALLOCATE ( l_scal         , STAT=istat)
  DEALLOCATE ( fc_min         , STAT=istat)

  DEALLOCATE ( grad           , STAT=istat)
  DEALLOCATE ( hig            , STAT=istat)

  DEALLOCATE ( tmps           , STAT=istat)
  DEALLOCATE ( vaps           , STAT=istat)
  DEALLOCATE ( prss           , STAT=istat)
  DEALLOCATE ( liqs           , STAT=istat)

  DEALLOCATE ( diss_tar       , STAT=istat)

! For turbdiff or vertical diffusion only:

  DEALLOCATE ( hlp            , STAT=istat)
  DEALLOCATE ( dicke          , STAT=istat)


! For turbdiff only:

  DEALLOCATE ( can            , STAT=istat)

  DEALLOCATE ( ftm            , STAT=istat)
  DEALLOCATE ( shv            , STAT=istat)


  DEALLOCATE ( hor_scale      , STAT=istat)
  DEALLOCATE ( xri            , STAT=istat)

  DEALLOCATE ( lay            , STAT=istat)
  DEALLOCATE ( lays           , STAT=istat)

  DEALLOCATE ( src            , STAT=istat)

  DEALLOCATE ( dzsm           , STAT=istat)
  DEALLOCATE ( dzsh           , STAT=istat)

  DEALLOCATE ( lev            , STAT=istat)

! For turbtran only:

  DEALLOCATE ( rclc           , STAT=istat)

  DEALLOCATE ( tketens_tar    , STAT=istat)

  DEALLOCATE ( lo_ice         , STAT=istat)

  DEALLOCATE ( k_2d           , STAT=istat)

  DEALLOCATE ( hk_2d          , STAT=istat)
  DEALLOCATE ( hk1_2d         , STAT=istat)

  DEALLOCATE ( h_top_2d       , STAT=istat)
  DEALLOCATE ( h_atm_2d       , STAT=istat)
  DEALLOCATE ( h_can_2d       , STAT=istat)

  DEALLOCATE ( edh            , STAT=istat)

  DEALLOCATE ( z0m_2d         , STAT=istat)
  DEALLOCATE ( z0d_2d         , STAT=istat)
  DEALLOCATE ( z2m_2d         , STAT=istat)
  DEALLOCATE ( z10m_2d        , STAT=istat)

  DEALLOCATE ( rat_m_2d       , STAT=istat)
  DEALLOCATE ( rat_h_2d       , STAT=istat)
  DEALLOCATE ( fac_h_2d       , STAT=istat)
  DEALLOCATE ( fac_m_2d       , STAT=istat)

  DEALLOCATE ( frc_2d         , STAT=istat)

  DEALLOCATE ( vel_2d         , STAT=istat)
  DEALLOCATE ( tl_s_2d        , STAT=istat)
  DEALLOCATE ( qt_s_2d        , STAT=istat)
  DEALLOCATE ( velmin         , STAT=istat)
  DEALLOCATE ( ratsea         , STAT=istat)

  DEALLOCATE ( dz_sg_m        , STAT=istat)
  DEALLOCATE ( dz_sg_h        , STAT=istat)
  DEALLOCATE ( dz_g0_m        , STAT=istat)
  DEALLOCATE ( dz_g0_h        , STAT=istat)
  DEALLOCATE ( dz_0a_m        , STAT=istat)
  DEALLOCATE ( dz_0a_h        , STAT=istat)
  DEALLOCATE ( dz_sa_h        , STAT=istat)
  DEALLOCATE ( dz_s0_h        , STAT=istat)

END SUBROUTINE turb_wkarr_dealloc
#endif

!==============================================================================

END MODULE turb_data
