!+ Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------

MODULE mo_data_turbdiff

!------------------------------------------------------------------------------
!
! Description:
!  This module contains constants, parameters switches and selectors
!  that are used in the turbulence MODULE 'turbulence_turbdiff'. 
!  With these variables a special configuration of the scheme is possible.
!  Each of the models 'COSMO' or 'ICON' may have its own specific MODULE
!  'turbulence_data'.
!  Notice that those control-variables, whose values may depend on the circumstances
!  of a single call of SUB 'organize_turbdiff' within one and the same model
!  are provided by the parameter-list of the SUB CALL.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  Matthias.Raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! @VERSION@    @DATE@     Matthias Raschendorfer
!  Initial release, based on the ICON-version of 'turbulence_turbdiff'
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

USE mo_kind, ONLY : &
    iintegers=> i4,    & ! KIND-type parameter for standard integer variables
    ireals   => wp       ! KIND-type parameter for real variables

USE mo_physical_constants, ONLY : &
 
! Physical constants and related variables:
! -------------------------------------------
 
    r_d     => rd,        & ! gas constant for dry air
    rdv,                  & ! r_d / r_v
    o_m_rdv,              & ! 1 - r_d/r_v
    rvd_m_o => vtmpc1,    & ! r_v/r_d - 1
    cp_d    => cpd,       & ! specific heat for dry air
    lh_v    => alv,       & ! evaporation heat
    rdocp   => rd_o_cpd,  & ! r_d / cp_d
    lhocp   => alvdcp,    & ! lh_v / cp_d
    rcpv,                 & ! cp_d/cp_v - 1
    rcpl,                 & ! cp_d/cp_l - 1
    con_m,                & ! kinematic vsicosity of dry air (m2/s)
    con_h,                & ! scalar conductivity of dry air (m2/s)
    t0_melt => tmelt,     & ! absolute zero for temperature  (K)
    grav,                 & ! acceleration due to gravity    (m/s2)
    p0ref,                & ! reference pressure for Exner-function

! Parameters for computing the saturation steam pressure (over water and ice):
! ------------------------------------------

    b3      => tmelt        ! absolute zero for temperature  (K)

USE mo_convect_tables, ONLY : &
 
    b1      => c1es,      & !
    b2w     => c3les,     & !                               
    b2i     => c3ies,     & !
    b4w     => c4les,     & !
    b4i     => c4ies,     & !
    b234w   => c5les        ! b2w * (b3 - b4w)
 
USE mo_physical_constants, ONLY : &

! Parameters for computing the rate of cloud cover based on relative humidity:
! ------------------------------------------

    ucl,                  & !
    uc1                     !

USE mo_math_constants, ONLY : &

    uc2     => sqrt3        !

USE mo_data_flake, ONLY : &

! Flake parameters:
! -------------------
 
    h_Ice_min_flk   ! Minimum ice thickness [m]

! Special parameter structures:
! -----------------------------------------------------

USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config

USE mo_turbdiff_config,      ONLY: turbdiff_config

! Switches controlling other parameterizations or the model configuration:
! ------------------------------------------

USE mo_lnd_nwp_config , ONLY :   &
 
    lseaice,      & ! forecast with sea ice model
    llake           ! forecast with lake model FLake

!==============================================================================

IMPLICIT NONE

PUBLIC

!==============================================================================

LOGICAL :: &

!Achtung: to be evaluated by ICON-switches <
  lconv    = .TRUE. , & ! convection scheme is active
  lsso     = .TRUE. , & ! SSO-Scheme is active
!>
 
  l3dturb  = .FALSE., & ! a run with 3D-(turbulent)-diffusion

  lscm     = .FALSE.    ! a SC-run

REAL (KIND=ireals), POINTER :: &

! Numerical constants:
!-----------------------------

  impl_weight(:) ! implicit weights for tridiagonal solver

!==============================================================================

! Parameters that may be used for tuning and special configurations:
!--------------------------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs!

!==============================================================================

REAL (KIND=ireals)     ::     &

! 1. Numerical parameters:
!-----------------------------

  impl_s     =  1.20_ireals,  & ! implicit weight near the surface (maximal value)
  impl_t     =  0.75_ireals,  & ! implicit weight near top of the atmosphere (maximal value)

  ! Minimal diffusion coefficients in [m^2/s] for vertical
  tkhmin     =  0.75_ireals,  & ! scalar (heat) transport
  tkmmin     =  0.75_ireals,  & ! momentum transport
  tkhmin_strat = 0.75_ireals, & ! scalar (heat) transport, enhanced value for stratosphere
  tkmmin_strat = 4.0_ireals,  & ! momentum transport,      enhanced value for stratosphere

  ditsmot    =  0.00_ireals,  & ! smoothing factor for direct time-step iteration

  tndsmot    =  0.00_ireals,  & ! vertical smoothing factor for diffusion tendencies
  frcsmot    =  0.00_ireals,  & ! vertical smoothing factor for TKE forcing
  tkesmot    =  0.15_ireals,  & ! time smoothing factor for TKE and diffusion coefficients
  stbsmot    =  0.00_ireals,  & ! time smoothing factor for stability function
  frcsecu    =  1.00_ireals,  & ! security factor for TKE-forcing       (<=1)
  tkesecu    =  1.00_ireals,  & ! security factor in  TKE equation      (out of [0; 1])
  stbsecu    =  0.00_ireals,  & ! security factor in stability function (out of [0; 1])

  epsi       =  1.0E-6_ireals   ! relative limit of accuracy for comparison of numbers

INTEGER            ::     &

  it_end     =  1           ! number of initialization iterations (>=0)

! 2. Parameters describing physical properties of the lower boundary 
!    of the atmosphere:
!------------------------------------------

REAL (KIND=ireals)     ::      &

  rlam_mom   =  0.0_ireals,    & ! scaling factor of the laminar boundary layer for momentum
  rlam_heat  =  1.0_ireals,    & ! scaling factor of the laminar boundary layer for heat

  rat_lam    =  1.0_ireals,    & ! ratio of laminar scaling factors for vapour and heat
  rat_can    =  1.0_ireals,    & ! ratio of canopy height over z0m
  rat_sea    = 10.0_ireals,    & ! ratio of laminar scaling factors for heat over sea and land

  z0m_dia    =  0.2_ireals,    & ! roughness length of a typical synoptic station [m]

  alpha0     =  0.0123_ireals, & ! Charnock-parameter
  alpha0_max =  0.0335_ireals, & ! upper limit of velocity-dependent Charnock-parameter
  alpha0_pert=  0.0_ireals,    & ! additive ensemble perturbation of Charnock-parameter
!<Achtung:
  alpha1     =  1.0000_ireals    ! parameter scaling the molek. roughness of water waves
!>

REAL (KIND=ireals)     ::      &

! 3. Parameters that should be external parameter fields being not yet 
!    available:
!------------------------------------------

  c_lnd      = 2.0_ireals,     & ! surface area density of the roughness elements over land
  c_sea      = 1.5_ireals,     & ! surface area density of the waves over sea
  c_soil     = 1.0_ireals,     & ! surface area density of the (evaporative) soil surface
  e_surf     = 1.0_ireals,     & ! exponent to get the effective surface area

! 4. Parameters that should be dynamical fields being not yet available:
!------------------------------------------

!!!DR in the long term, we should make use of tf_salt (see mo_physical_constants in ICON)
  zt_ice     = -1.7_ireals,    & !freezing temperature of sea ice
  z0_ice     =  0.001_ireals     !roughness length of sea ice

REAL (KIND=ireals)     ::     &

! 5. Parameters for modelling turbulent diffusion:
!------------------------------------------

  tur_len    = 500.0_ireals,  & ! asymptotic maximal turbulent distance [m]
  pat_len    = 100.0_ireals,  & ! effective length scale of subscale surface patterns over land [m]
                                ! (should be dependent on location)
  len_min    =  1.0E-6_ireals,& ! minimal turbulent length scale [m]

  vel_min    =  0.01_ireals,  & ! minimal velocity scale [m/s]

  akt        =  0.4_ireals,   & ! von Karman-constant

  ! Length scale factors for pressure destruction of turbulent
  a_heat     =  0.74_ireals,  & ! scalar (heat) transport
  a_mom      =  0.92_ireals,  & ! momentum transport

  ! Length scale factors for dissipation of
  d_heat     =  10.1_ireals,  & ! scalar (temperature) variance
  d_mom      =  16.6_ireals,  & ! momentum variance

  ! Length scale factors for turbulent transport (vertical diffusion)
  c_diff     =  0.20_ireals,  & ! of TKE

  ! Length scale factor for separate horizontal shear production
  a_hshr     =  1.00_ireals,  & ! of TKE

  ! Length scale factor for the stability correction
  a_stab     =  0.00_ireals,  & ! no stability correction so far

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):
  clc_diag   =  0.5_ireals,   & !cloud cover at saturation
  q_crit     =  1.6_ireals,   & !critical value for normalized over-saturation
  c_scld     =  1.0_ireals      !factor for liquid water flux density in sub grid scale clouds

!==============================================================================

! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

LOGICAL :: &

  ltkesso       =.TRUE. , & ! calculation SSO-wake turbulence production for TKE
  ltkecon       =.FALSE., & ! consider convective buoyancy production for TKE
  ltkeshs       =.TRUE. , & ! consider separ. horiz. shear production for TKE
  loutshs       =.TRUE., & ! consider separ. horiz. shear production of TKE for output
  
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

! Notice that the following switches are provided by the parameter-list of SUB 'organize_turbdiff':

! lstfnct                   :calculation of stability function required
! lnsfdia                   :calculation of (synoptical) near-surface variables required
! lmomdif                   :calculation of complete gradient diffusion of horizontal momenum
! lscadif                   :calculation of complete gradient diffusion of scalar properties
! lturatm                   :running turbulence model between atmosph. layers (updating diffusion coefficients)
! ltursrf                   :running turbulence model at the surface layer (updating transfer coefficients
! lsfluse                   :use explicit heat flux densities at the suface
! ltkeinp                   :TKE present as input (at level k=ke1 for current time level 'ntur')
! lgz0inp                   :gz0 present as input

INTEGER :: &

! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

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
  itype_sher    =0, & ! type of shear production for TKE
                      ! 0: only vertical shear of horizontal wind
                      ! 1: previous plus horizontal shear correction
                      ! 2: previous plus shear from vertical velocity
  itype_diag_t2m=1    ! type of diagnostics of 2m-temperature and -dewpoint
                      ! 1: Considering a fictive surface roughness of a SYNOP lawn
                      ! 2: Considering the mean surface roughness of a grid box
                      !    and using an exponential roughness layer profile

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
!< Achtung:
  ilow_def_cond =2, & ! type of the default condition at the lower boundary
                      ! 1: zero surface flux density
                      ! 2: zero surface value
  imode_calcirc =2, & ! mode of treating the circulation term (related to 'pat_len', imode_pat_len')
                      ! 1: explicit calculation of the flux convergence
                      ! 2: quasi implicit treatment by calculation of effective TKE-gradients
  imode_pat_len =2, & ! mode of determining the length scale of surface patterns (related to 'pat_len')
                      ! 1: by the constant value 'pat_len' only
                      ! 2: and the std. deviat. of SGS orography as a lower limit (only if 'd_pat' is pres.)
!Achtung: bislang "1"
  imode_frcsmot =2, & ! if "frcsmot>0", apply smoothing of TKE source terms 
                      ! 1: globally or 
                      ! 2: in the tropics only (if 'trop_mask' is present) 
  imode_shshear =2, & ! mode of calculat. the separated horizontal shear mode (related to 'ltkeshs', 'a_hshr')
                      ! 1: with a constant lenght scale 
                      ! 2: with a Ri-dependent length sclale correction
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
!>
!Notice that the following selectors are provided by the parameter-list of SUB 'organize_turbdiff:

! iini                :type of initialization (0: no, 1: separate before the time loop
!                                                   , 2: within the first time step)
! itnd                :type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
!                                                     2: by adding to current profile before vertical diffusion
!                                                     3: by using corrected virtual vertical profiles

!==============================================================================

!-----------------------------------------------------------------------------
 CONTAINS
!-----------------------------------------------------------------------------

 SUBROUTINE get_turbdiff_param (jg)

     INTEGER (KIND=iintegers), INTENT(IN) :: jg !patch index

     impl_weight => turbdiff_config(jg)%impl_weight

     lsso         =(atm_phy_nwp_config(jg)%inwp_sso.GT.0)
     lconv        =(atm_phy_nwp_config(jg)%inwp_convection.GT.0)

     imode_tran   = turbdiff_config(jg)%imode_tran
     icldm_tran   = turbdiff_config(jg)%icldm_tran
     imode_turb   = turbdiff_config(jg)%imode_turb
     icldm_turb   = turbdiff_config(jg)%icldm_turb
     itype_sher   = turbdiff_config(jg)%itype_sher
     imode_frcsmot= turbdiff_config(jg)%imode_frcsmot

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

     loutshs      = ltkeshs .OR. itype_sher > 0

 END SUBROUTINE get_turbdiff_param

!==============================================================================

END MODULE mo_data_turbdiff
