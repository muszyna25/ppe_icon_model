!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!-------------------------------------------------------------------------------

#ifdef __COSMO__
 MODULE turbulence_turbdiff
#endif
#ifdef __ICON__
 MODULE src_turbdiff
#endif

!-------------------------------------------------------------------------------
!
! Description:
!
! The module "turbdiff" calculates the tendencies for turbulent
! vertial transport of momentum and heat and the coefficients
! for turbulent diffusion and turbulent transfer as well.

! The clousure is made on lever 2.5 (Mellor/Yamada) using a prognostik
! TKE-equation and includes the formulation of a flow through a porous
! medium (roughnes layer)

! The turbulence model (with some Prandtl-layer approximations is used
! for the calculation of turbulent transfer between atmosphere and the
! lower boundary too.

! The module contains only the subroutines :
!
!   init_canopy, organize_turbdiff

! called from the basic driving routine of the model.

! 'organize_turbdiff' contains the further subroutines:
!
!   canpoy_source, turb_param, turbtran, turbdiff and stab_funct.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  matthias.raschendorfer@dwd.de
!
!--------------------------------------------------------------------------
! History of module 'src_turbdiff' in COSMO:
! (during this time the subroutines 'turbtran' and turbdiff' had been located
!  in separate indluce-files with their own histories):
!--------------------------------------------------------------------------

! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 Matthias Raschendorfer
!  Initial release
! 1.33       1999/10/14 Matthias Raschendorfer
!  USE 2 additional LOGICAL namelist-parameter controlling the physics
!  (ltmpcor,lprfcor).
!  USE 6 additional REAL-namelist-parameter controlling the turbulence:
!  (tur_len, a_heat, d_heat, a_mom, d_mom, c_diff, rat_can, c_lnd, c_see).
!  USE additional field: sai (for the calculation of the roughness length
!  for scalars in SUBROTINE turbtran).
!  Introduction of the LOGECAL parameter lstfnct in the parameter list of
!  SUBROUTINE turbdiff.
! 1.34       1999/12/10 Ulrich Schaettler
!  Consideration of the LAI in the calculation of the SAI
!  Put allocation of canopy fields to new module src_allocation
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names.
! 2.2        2000/08/18 Matthias Raschendorfer
!  USE of an additional REAL-namelist-parameter 'c_soil', defintion of
!  the additional 2-D fields 'eai' and 'tai'.
!  Definition of 'sai', 'eai' and 'tai' also in the case of the
!  previous surface scheme (itype_tran=1).
! 2.15       2002/03/20 Matthias Raschendorfer
!  Introduction of the roughness length for a typical SYNOP-station (z0m_dia)
! 2.19       2002/10/24 Ulrich Schaettler
!  Use l2tls for specifying the number of time levels
! 3.5        2003/09/02 Ulrich Schaettler
!  Adapted interface for routine exchg_boundaries
! 3.7        2004/02/18 Matthias Raschendorfer
!  Introduction of the parameter rat_sea
! 3.14       2005/01/25 Jochen Foerstner
!  Introduced prognostic treatment of TKE also for 3D turbulence scheme
!  (use lprog_tke)
! 3.16       2005/07/22 Matthias Raschendorfer
!  Introduction of the parameters tkesmot, wichfakt, securi, tkhmin, tkmmin
! 3.18       2006/03/03 Matthias Raschendorfer
!  Use rh_2m from data_fields
!  Eliminated use of exchg_boundaries (not necessary here)
! 3.21       2006/12/04 Dmitrii Mironov, Ulrich Schaettler
!  Changes for the FLake Model
!  Use tuning variables from data_turbulence
!  Use new NL variables from data_constants: clc_diag, q_crit
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introducing 'tfv' containing the laminar reduction factor for evaporation.
!  Introduction of some output variables for SCLM.
!  Moving some parameters from 'turb_param.incf' to MODULE data_turbulence.
!  Moving the deduced parameters from 'turb_param.incf' to SUBROUTINE turb_param.
!  'turb_param.incf' contains only statement functions and is renamed into 'stab_funct.incf'
!  Moving content of 'stab_funct.incf' into the SUBROUTINE stab_funct.
! V3_24        2007/04/26 Ulrich Schaettler
!  Introduced call to exchg_boundaries for wind-tendencies again;
!  it is necessary for imode_turb=2/3!
! V4_3         2008/02/25 Matthias Raschendorfer
!  Introduction of a 3D diagnostic field 'edr' for the eddy dissipotion rate.
!  Changing the treatment of parameter field dd(:,:) to achieve better vectorization
!  in SUBROUTINE 'stab_funct'.
! V4_4         2008/07/16 Ulrich Schaettler
!  Removed POINTERs from SR stab_funct, because these still prevented vectorization
! V4_8         2009/02/16 Ulrich Schaettler
!  Included itype_diag_t2m to add options for computation of 2m temperature
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Introduction of the INTEGERs im, jm for SCLM treatment.
!  Introduction of LOGICALs limpltkediff, ltkesso and removing lturhor.
!  Introduction of INTEGER itype_sher and of REALs a_hshr, edadlat, acrlat and a_stab.
!  Reference of ut_sso, vt_sso and lsso.
!  Modifications for new sea-ice model: eliminated l_ls_ice, introduced lseaice
!    (Jan-Peter Schulz)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN

!--------------------------------------------------------------------------
! History of module 'src_turbdiff' in ICON (for all including subroutines)
!  of some major changes by Matthias Raschendorfer:
!--------------------------------------------------------------------------

!              2010/09/30 Matthias Raschendorfer
!  Substitution of 'itype_diag_t2m' by 'itype_synd' being already used for that purpose
!   "1": alternative SYNOP-digansostics related to previous Lewis scheme (not included here)
!   "2": SYNOP-diagnostics according to current transfer scheme using SYNOP-z0 'z0d'
!   "3": like "2" but using 'z0d' only for 10m-wind and a specific roughness layer profile for
!        2m-temperature and -humidity.
!  Including the adiabatic lapse rate correction for 't_2m' for "itype_synd.EQ.3" as well.
!              2010/12/17 Matthias Raschendorfer
!  Introduction of a TKE-source term due to scale interaction with sub grid scale convection
!   using 'ltkecon' and the convective buoyant heat flux density in 'tket_conv'.
!              2010/12/30 Matthias Raschendorfer
!  Reorganization of the code for use in both models COSMO and ICON with various formal modifications:
!   MODULE 'src_turbdiff' CONTAINS now SUBs 'init_canopy', 'organize_turbdiff' and 'turb_cloud'.
!   The latter was before in 'meteo_utilities' with name 'cloud_diag' and contains a saturation
!   adjustment with due regard on turbulent fluctuations of thermodynamic model variables.
!   SUB 'organize_turbdiff' CONTAINS SUBs
!    'turbtran', 'turbdiff', 'turb_param', 'stab_funct', 'diag_level' and ('canopy_source').
!    Exept 'diag_level' they had been present before as well, partly in INCLUDE-files.
!   In accordance with ICON rules, SUB 'organize_turbdiff' is CALLed with the most needed parameters
!    in the SUB header. Some parameters are OPTIONAL, allowing to run the scheme in differen modes.
!              2011/02/18 Matthias Raschendorfer
!  Introduction of some minor formal modif. of some parts of the code mainly to achiev better vectorization
!   (results will be modyfied only because of numerical effects).
!  Introduction of OPTIONAL extra output fields 'tket_sso' and 'tket_hshr'.
!              2011/03/23 Matthias Raschendorfer
!  Correction of two bugs introduced together with the last modifications in SUB 'turbtran'
!   (related to SUB 'diag_level' and the 'tet_l'-gradient used for the flux output).
!  Substitution of run time allocations because of bad performance on some computers.
!              2011/06/17 Matthias Raschendorfer
!  Removing a bug related to the saving of the saturation fraction that had been introduced by
!  the ICON modifications.
!              2011/08/26 Matthias Raschendorfer
!  Discriminating the variable for standard deviation of oversaturation (input) and cloud fraction (output)
!   in the CALL of 'turb_cloud'.
!  Changing the the definition of LOGICAL 'lcircterm'.
!  Introduction of a preconditioning of the tridiagonal system for "imode_turb=3".
!              2011/09/23 Matthias Raschendorfer
!  Changing 'k' to 'ke' in the initialization section of 'turbdiff' in case of "itype_tran.EQ.3".
!  Formulating temp-grad for 'ltmpcor' for "k.EQ.ke1" similar to the other levels.
!  Removing a wrong multiplication with Exner-factor in case of "imode_turb.EQ.4".
!  Adding a missing multiplication by 'dicke()' in case of ".NOT.limpltkediff".
!  Dirscriminating between two effective values of effictive Prandtl-layer depth,
!   by introducing the arrays 'dz0(:,mom|sca)' and 'vh0(:,:)
!   removing wrong surface conditions in case "imode_turb.GE.3 .AND. itype_tran.EQ.2".
!              2011/09/28 Matthias Raschendorfer
!  Correction in formula for specific humidity at saturation by using partial pressure of dry air.
!              2011/12/08 Matthias Raschendorfer
!  Introduction of SUBs 'prep_impl_vert_diff' and 'calc_impl_vert_diff' substituting the code segment for
!   implicit vertical diffusion calcul. for TKE or all other variables, if 'imode_turb' is one of 3 or 4.
!   These SUBs allow the calculation of semi implicit diffusion for arbitrary variables defined on
!   full- or  half-levels.
!   Diffusion tendency (rather than flux-density) of pot. temp. is converted in that of ordinary temp.
!              2011/12/27 Matthias Raschendorfer
!  Correction of some bugs:
!   Including missing definition of 'km1' for the surface level in SUB calc_impl_vert_diff
!   Moving a wrong ')' in interpolation of 't' onto the lower boundary of model atmosphere.
!              2012/01/10 Matthias Raschendorfer
!  Correction of further bugs:
!   In case of "limpltkediff" 'ko' needs to be added by "1".
!   Coding mistake in n-loop bounds for surface layer gradients in 'turbdiff'.
!  Introduction of array 'rhoh' (density on main levels) to avoid a back interpolation from half levels.
!  Using a mass weighted density interpolation onto half levels as for all other variables.
!              2012/01/11 Matthias Raschendorfer
!  Some rearrangement and simplification:
!   Removing explicit vertical diffusion and expressing explicit (eg. moist) corrections
!    by corrected vertical profiles still using implicit formulation of vertical diffusion.
!   Now "imode_turb=2" is for semi-implicit diffusion within 'turbdiff' using a concentration condition
!    at the surface and "imode_turb=3" is for the same with a flux condition at the surface.
!    "imode_turb_turb=4" is obsolet now.
!   One loop for diffusion of all first order model variables.
!              2012/01/20 Matthias Raschendorfer
!  Reformulation of roughness layer drag by an implicit equation (keeping sign of wind components)
!  Expressing explicit circulation term (as far as it has divergence form) by a respective
!   q**2-profile correction within the implicit formulation of diffusion similar to the treatment
!   of the former explicit (moist) corrections
!  Introduction of effective flux profile calculation into SUB 'calc_impl_vert_diff'.
!  Setup of a positiv definit solution avoiding some of the former limitations.
!              2012/03/20 Matthias Raschendorfer
!  Rearrangement, modularization and revision of numerical treatment, such as:
!  Intoduction of the SUBs 'solve_turb_budgets' and 'adjust_satur_equil' in order to use the same code
!  for the same purpose in 'turbtran' and 'turbdiff.
!  Removing the explicit diffusion option.
!  Including the new driving SUB 'vert_grad_diff' oranizing vertical diffusion when called in a loop of
!  several full level variables and managing options how to treat non-gradient-fluxes and additional
!  explicit time tendencies numerically, including smoothing options and preconditioning the linear system.
!  Parameters 'imode_tran' and 'imode_turb' define whether the TKE-equation is solved in a dianostic (1)
!  or a prognostic (2) form.
!  New parameter 'lsflcnd' is a flag for using a lower flux condition for vertical diffusion.
!  Introducing the flags 'lturatm', 'ltursrf', 'lmomdif', 'lscadif' arranging what tasks of
!  'organize_turbdiff' are required.
!  Introduction of 'ltkeshs' and removing the case "itype_sher.EQ.3".
!              2014/07/28 Matthias Raschendorfer
!  Introduction of precalculated 'hdef2', 'hdiv', 'dwdx', 'dwdy' used for calculation of horizontal shear
!   including the scale separated non-turbulent part controlled by 'itype_sher'.
!  Some formal cleaning mainly in order to avoid code duplication
!   -> no influence of results, except due to non-associative multiplication using 'rprs'
!      and non-reversible calculation of the reciprocal: "z1/(temp/exner).NE.exner/temp"
!      both in SUB 'adjust_satur_equil'
!  Removing a bug in formular for 'rcld(:,ke1)' in SUB 'turbtran'
!   -> influence on near-surface temperature and - humidity
!  Simpler (physically identical) formulation of moist flux conversion
!   -> only numerical differences in the case of 'lexplcor'
!  Numerically more efficient formulation of Blackadar 'len_scale'
!   -> only numerical differences
!  Eliminating array 'wind', as 'u' and 'v' are already defined at mass positions.

!--------------------------------------------------------------------------
! History of the common module 'turbulence_turbdiff' in COSMO and ICON:
!--------------------------------------------------------------------------

! @VERSION@    @DATE@     Matthias Raschendorfer
! Blocked and further developed version of the turbulence model (including the calculation of transfer
!  resistances) and a new more general formulation of turbulent vertical diffusion using all turbulence data
!  from module 'turbulence_data' including various switches and selectors in order to configure the turbulence model
!  either as currently used in ICON or COSMO:
! Adopting other development within the ICON-version (namely by Guenther Zaengl) as switchable options
!  related to the following new selectors and switches:
!   imode_pat_len, imode_frcsmot, imode_shshear, imode_tkvmini, imode_rat_sea, imode_vel_min, imode_charpar,
!   lfreeslip.
! Rearranging the development by Matthias Raschendorfer that had not yet been transferred to COSMO as switchable
!  options related to the following switches:
!   lsflcnd, ldynimp, lprecnd, ltkeshs, loutshs
!  and selectors:
!   imode_stbcorr, itype_diag_t2m, imode_syndiag, imode_qvsatur, imode_stadlim, imode_trancnf,
!   imode_tkediff, imode_adshear, imode_tkemini, imode_lamdiff
!  and a partly new (more consistent) interpretation of:
!   imode_tran, imode_turb, icldm_tran, icldm_turb, itype_sher, and itype_diag_t2m
! Controlling numerical restrictions gradually namely by the parameters:
!  ditsmot, tndsmot, frcsmot, tkesmot, stbsmot, frcsecu, tkesecu, stbsecu
! Adopting the 3D-options from the COSMO-Version and correcting the horizontal limit of the turbulent length scale.
! Correction of some bugs (of the latest ICON-version) in particular related to the optional lower concentration condition.
! Using the arrays 'tvm', 'tvh' and 'tkm', allowing an easier formulation of transfer-resistances.
!
! Code Description:
! Language: Fortran 90
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".


! Modules used:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
    ireals,       & ! KIND-type parameter for real variables
    iintegers,    & ! KIND-type parameter for standard integer variables
!
! Physical constants and related variables:
! -------------------------------------------
!
    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat for dry air
    lh_v,         & ! evaporation heat
    rdocp,        & ! r_d / cp_d
    lhocp,        & ! lh_v / cp_d
    rcpv,         & ! cp_v/cp_d - 1
    rcpl,         & ! cp_l/cp_d - 1 (where cp_l=cv_l)
    con_m,        & ! kinematic vsicosity of dry air (m2/s)
    con_h,        & ! scalar conductivity of dry air (m2/s)
    t0_melt,      & ! absolute zero for temperature (K)
    grav,         & ! acceleration due to gravity
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    p0ref,        & ! reference pressure for Exner-function
!
    b1,           & ! variables for computing the saturation steam pressure
    b2w,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b234w,        & ! b2w * (b3 - b4w)
    uc1,          & ! variable for computing the rate of cloud cover in
    uc2,          & ! the unsaturated case
    ucl             !               -- " --

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Switches/selectors and parameters controlling other physical parameterizations
! or the overall model configuration:
!
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecast with lake model FLake
    lsso,         & ! SSO-Scheme is active
    lconv,        & ! confection scheme is active
!
    lscm,         & ! a SC-run
!
    h_Ice_min_flk   ! Minimum ice thickness [m]

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Numerical constants and parameters:
! -----------------------------------
!
    impl_weight,  & ! vertical column with implicit weights for tridiagonal solver
!
    impl_s,       & ! implicit weight near the surface (maximal value)
    impl_t,       & ! implicit weight near top of the atmosphere (minimal value)
!
    tkhmin,       & ! minimal diffusion coefficients for heat
    tkmmin,       & ! minimal diffusion coefficients for momentum
    tkhmin_strat, & ! additional minimal diffusion coefficients for heat for stratosphere
    tkmmin_strat, & ! additional minimal diffusion coefficients for momentum for stratosphere
!
    ditsmot,      & ! smoothing factor for direct time-step iteration
    tndsmot,      & ! vertical smoothing factor for diffusion tendencies
    frcsmot,      & ! vertical smoothing factor for TKE forcing
    tkesmot,      & ! time smoothing factor for TKE
    stbsmot,      & ! time smoothing factor for stability function
    frcsecu,      & ! security factor for TKE-forcing       (<=1)
    tkesecu,      & ! security factor in  TKE equation      (out of [0; 1])
    stbsecu,      & ! security factor in stability function (out of [0; 1])
!
    epsi,         & ! relative limit of accuracy for comparison of numbers
!
    it_end          ! number of initialization iterations (>=0)

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Parameters describing physical properties of the lower boundary
! of the atmosphere:
!---------------------------------------------------------------
!
    rlam_mom,     & ! scaling factor of the laminar boudary layer for momentum
    rlam_heat,    & ! scaling factor of the laminar boudary layer for heat
!
    rat_can,      & ! factor for the canopy height
    rat_sea,      & ! ratio of laminar scaling factors for heat over sea and land
    rat_lam,      & ! ratio of laminar scaling factors for vapour and heat
!
    z0m_dia,      & ! roughness length of a typical synoptic station
!
    alpha0,       & ! lower bound for Charnock-parameter
    alpha0_max,   & ! upper bound for Charnock-parameter
    alpha0_pert,  & ! ensemble perturbation for Charnock-parameter
    alpha1,       & ! parameter scaling the molek. roughness of water waves
!
    c_lnd,        & ! surface area index of the land exept the leaves
    c_soil,       & ! surface area index of the (evaporative) soil
    c_sea,        & ! surface area index of the waves over sea
    e_surf,       & ! exponent to get the effective surface area
!
    zt_ice,       & ! freezing temperature of sea ice
    z0_ice,       & ! roughness length of sea ice
!
    tur_len,      & ! maximal turbulent length scale [m]
    pat_len,      & ! lenth scale of subscale patterns over land [m]
    len_min,      & ! minimal turbulent length scale [m]
    vel_min,      & ! minimal velocity scale [m/s]
!
    akt,          & ! von Karman-constant
!
    a_h=>a_heat,  & ! factor for turbulent heat transport
    a_m=>a_mom,   & ! factor for turbulent momentum transport
    d_h=>d_heat,  & ! factor for turbulent heat dissipation
    d_m=>d_mom,   & ! factor for turbulent momentum dissipation
!
    c_diff,       & ! factor for turbulent diffusion of TKE
    a_hshr,       & ! factor for horizontal shear production of TKE
    a_stab,       & ! factor for stability correction of turbulent length scale
!
    clc_diag,     & ! cloud cover at saturation in statistical cloud diagnostic
    q_crit,       & ! critical value for normalized over-saturation
    c_scld          ! factor for liquid water flux density in sub grid scale clouds

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------
!
    ltkesso,      & ! consider SSO-wake turbulence production of TKE
    ltkecon,      & ! consider convective buoyancy production of TKE
    ltkeshs,      & ! consider separ. horiz. shear production of TKE
    loutshs,      & ! consider separ. horiz. shear production of TKE for output
!
    lnonloc,      & ! nonlocal calculation of vertical gradients used for turb. diff.
    lprfcor,      & ! using the profile values of the lowest main level instead of
!                   ! the mean value of the lowest layer for surface flux calulations
    ltmpcor,      & ! consideration of thermal TKE-sources in the enthalpy budget
    lcpfluc,      & ! consideration of fluctuations of the heat capacity of air
!
    lexpcor,      & ! explicit corrections of the implicit calculated turbul. diff.
!
    l3dturb,      & ! a model run with 3D-(turbulent)-diffusion
!
!   for semi-implicit vertical diffusion:
    lsflcnd,      & ! lower flux condition
    ldynimp,      & ! dynamical calculation of implicit weights
    lprecnd,      & ! preconditioning of tridiagonal matrix
    lfreeslip       ! free-slip lower boundary condition (enforeced zero-flux condition for
                    ! for all diffused variables, only for idealized test cases)

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------
!
    imode_tran,   & ! mode of TKE-equation in transfer scheme                    (compare 'imode_turb')
    imode_turb,   & ! mode of TKE-equation in turbulence scheme
                    !  0: diagnostic equation
                    !  1: prognostic equation (default)
                    !  2: prognostic equation (implicitly positive definit)
    icldm_tran,   & ! mode of water cloud representation in transfer parametr.   (compare 'icldm_tran)
    icldm_turb,   & ! mode of water cloud representation in turbulence parametr.
                    ! -1: ignoring cloud water completely (pure dry scheme)
                    !  0: no clouds considered (all cloud water is evaporated)
                    !  1: only grid scale condensation possible
                    !  2: also sub grid (turbulent) condensation considered
    itype_wcld,   & ! type of water cloud diagnosis within the turbulence scheme:
                    ! 1: employing a scheme based on relative humitidy
                    ! 2: employing a statistical saturation adjustment
    itype_sher,   & ! type of shear production for TKE
                    ! 0: only vertical shear of horizontal wind
                    ! 1: previous plus horizontal shear correction
                    ! 2: previous plus shear from vertical velocity
    imode_stbcorr,& ! mode of correcting the stability function (related to 'stbsecu')
                    ! 1: always for strict.-non-stb. strat. using a restr. gama in terms of prev. forc.
                    ! 2: only to avoid non-physic. solution or if current gama is too large
    itype_diag_t2m  ! type of diagnostics of 2m-temperature and -dewpoint
                    ! 1: Considering a fictive surface roughness of a SYNOP lawn
                    ! 2: Considering the mean surface roughness of a grid box
                    !    and using an exponential roughness layer profile

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
    ilow_def_cond,& !type of the default condition at the lower boundary
                    ! 1: zero surface flux density
                    ! 2: zero surface value
    imode_calcirc,& ! mode of treating the circulation term (related to 'pat_len', imode_pat_len')
                    ! 1: explicit calculation of the flux convergence
                    ! 2: quasi implicit treatment by calculation of effective TKE-gradients
    imode_pat_len,& ! mode of determining the length scale of surface patterns (related to 'pat_len')
                    ! 1: by the constant value 'pat_len' only
                    ! 2: and the std. deviat. of SGS orography as a lower limit (only if 'd_pat' is pres.)
    imode_frcsmot,& ! if frcsmot>0, apply smoothing of TKE source terms
                    ! 1: globally or
                    ! 2: in the tropics only (if 'trop_mask' is present)
    imode_shshear,& ! mode of calculat. the separated horizontal shear mode related to 'ltkeshs', 'a_hshr')
                    ! 1: with a constant lenght scale
                    ! 2: with a Ri-dependent length sclale correction
    imode_tkesso,&  ! mode of calculat. the SSO source term for TKE production
                    ! 1: original implementation
                    ! 2: with a Ri-dependent reduction factor for Ri>1
    imode_tkvmini,& ! mode of calculating the minimal turbulent diff. coeffecients
                    ! 1: with a constant value
                    ! 2: with a stability dependent correction
    imode_rat_sea,& ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                    ! 1: constant ratio compared to land surface
                    ! 2: with a correction for a strongly overheated SST
    imode_vel_min,& ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                    ! 1: with a constant value
                    ! 2: with a stability dependent correction
    imode_charpar   ! mode of estimating the Charnock-Parameter
                    ! 1: use a constant value
                    ! 2: use a wind-dependent value with a constant lower bound

#ifdef __COSMO__
USE turbulence_data, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
    imode_syndiag,& ! mode of diagnostics at the synoptic near surface levels (related to 'itype_diag_t2m')
                    ! 1: direct interpolation of temperature and specific humidity
                    ! 2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                    !    allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
    imode_qvsatur,& ! mode of calculating the saturat. humidity
                    ! 1: old version, using total pressure
                    ! 2: new version, using partial pressure of dry air
    imode_stadlim,& ! mode of mode of limitting statist. saturation adjustment
                    ! 1: only absolut upper limit of stand. dev. of local oversatur. (sdsd)
                    ! 2: relative limit of sdsd and upper limit of cloud-water
    imode_trancnf,& ! mode of configuring the transfer-scheme
                    ! 1: old version: start. with lamin. diffus.; with a lamin. correct. for profile-funct.;
                    !    interpol. T_s rather then Tet_l onto zero-level; calcul. only approx. Tet_l-grads.;
                    !    using an upper bound for TKE-forcing; without transmit. skin-layer depth to turbul.
                    ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                    !    interpol. Tet_l onto zero-level; calcul. Tet_l-gradients directly;
                    !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                    ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                    !    for stable stratification
                    ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
    imode_tkediff,& ! mode of implicit TKE-Diffusion (related to 'c_diff')
                    ! 1: in terms of q=SQRT(2*TKE))
                    ! 2; in terms of TKE=0.5*TKE**2
    imode_adshear,& ! mode of considering additional shear by scale interaction
                    ! 1: not considered for stability functions
                    ! 2:     considered for stability functions
    imode_tkemini,& ! mode of fixing a lower limit of q=2TKE**2
                    ! 1: by using 'vel_min' only
                    ! 2: by adapting to minimal diffusion coefficients
    imode_lamdiff   ! mode of considering laminar diffusion at surface layer
                    ! 1: only when calculating the profile functions
                    ! 2: surface-layer diffusion coeff. always at least at laminar value

!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
!
    lsclm, lsurflu, latmflu, i_cal, i_upd, i_mod, im, &
!
    UUA, VVA, UWA, VWA, WWA, UST, TWA, QWA, TTA, TQA, QQA, SHF, LHF, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
#endif
!SCLM---------------------------------------------------------------------------

!-------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC  :: init_canopy, organize_turbdiff, turb_cloud, vert_grad_diff, &
           modvar


INTEGER (KIND=iintegers), PARAMETER :: &
!
    mom=1,        & ! index for a momentum variable
    sca=2           ! index for a scalar   variable

REAL (KIND=ireals), PARAMETER :: &
!
    z0 = 0.0_ireals,&
    z1 = 1.0_ireals,&
    z2 = 2.0_ireals,&
    z3 = 3.0_ireals,&
    z4 = 4.0_ireals,&
    z5 = 5.0_ireals,&
    z6 = 6.0_ireals,&
    z7 = 7.0_ireals,&
    z8 = 8.0_ireals,&
    z9 = 9.0_ireals,&
    z10=10.0_ireals,&
!
    z1d2=z1/z2     ,&
    z1d3=z1/z3     ,&
    z2d3=z2/z3     ,&
    z3d2=z3/z2

REAL (KIND=ireals) :: xx

INTEGER (KIND=iintegers) :: &
!
    istat=0, ilocstat=0

LOGICAL :: &
!
    lerror=.FALSE.

TYPE modvar !model variable
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: av(:,:) => NULL() !atmospheric values
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: sv(:)   => NULL() !surface     values (concentration of flux density)
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: at(:,:) => NULL() !atmospheric time tendencies
     LOGICAL                                 :: fc                !surface values are flux densities
END TYPE modvar

TYPE turvar !turbulence variables
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: tkv(:,:) => NULL() !turbulent coefficient for vert. diff.
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: tsv(:)   => NULL() !effective surface layer depth
END TYPE turvar

TYPE varprf !variable profile
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: bl(:,:) !variable at boundary model levels
     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
          :: ml(:,:) !variable at main     model levels
END TYPE varprf
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------

!********************************************************************************
!********************************************************************************

!+ Module procedure init_canopy in "src_turbdiff" for initialization and allocation
!+ of special external parameters describing the surface canopy needed for for the
!+ description of surface-to-atmosphere transfer and within canopy diffusion:

SUBROUTINE init_canopy ( ie, ke, ke1, kcm, &
!
           i_stp, i_enp, icant, &
!
           l_hori, hhl, fr_land, plcov, &
!
           lai, sai, tai, eai, &
!
           d_pat, l_pat, h_can,  &
!
           c_big, c_sml, r_air )

!_________________________________________________________________________________
!
! Description:
!
!   In the module 'init_canopy' additional external parametr fields, which are
!   used in the new turbulence scheme 'turbdiff' (especially parameters for the
!   physical description of the roughness canopy) are covered by the appropriate values
!   by reading the refering parameter files and/or by making some diagnostic calculations
!   using the known parameters.
!   The 3-d fields of the canopy parameters are dynamically allocated using the maximum
!   canopy hight within the model domain.

! Method:
!
!   For the present there exists no concept of generating those additional data. Thus the
!   model will run either with an artificial canopy arcitecture or (for simplicity)
!   without any vertically resolved canopy. But even in the latter case at least the
!   allocation of the 3-d canopy data fields must be done, because the canopy concept
!   is incorporated in the tubulent diffusion scheme 'turbdiff'.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8236 1493
!  email:  mraschendorfer@dwd.d400.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================

!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!Formal Parameters:
!-------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------
!
    ie,    & ! number of grid points in zonal      direction
    ke,    & ! number of main model levels (start index is 1)
    ke1,   & ! number of half model levels (start index is 1)
    i_stp, & ! horizontal start-index (including the model bondary lines)
    i_enp    ! horizontal   end-index (including the model bondary lines)

INTEGER (KIND=iintegers), OPTIONAL, INTENT(IN) :: &
!
    icant    ! index for the used canopy-type
             ! 1: evapotransp.-fractions only based on plant-cover
             ! 2: based on a surface-area-index for all evapotransp.-types

INTEGER (KIND=iintegers), TARGET, INTENT(INOUT) :: &
!
    kcm             ! index of the lowest model layer higher than the canopy

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
!
! External parameter fields:
! ----------------------------
    fr_land         ! land portion of a grid point area             ( 1 )

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
!
    l_hori,       & ! horizontal grid spacing (m)
!
    plcov,        & ! fraction of plant cover                       ( 1 )
    lai,          & ! leaf area index                               ( 1 )
!
    h_can,        & ! hight of the vertically resolved canopy
    d_pat           ! external geometric dimension of circulation patterns

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: &
!
    sai,          & ! surface area index                            ( 1 )
    tai,          & ! transpiration area index                      ( 1 )
    eai,          & ! (evaporative) earth area index                ( 1 )
!
    l_pat           ! effective length scale of circulation patterns

REAL (KIND=ireals), DIMENSION(:,kcm:), OPTIONAL, INTENT(INOUT) :: &
!
    c_big,        & ! effective drag coefficient of canopy elements
                    ! larger than or equal to the turbulent length scale (1/m)
    c_sml           ! effective drag coefficient of canopy elements
                    ! smaller than the turbulent length scale            (1/m)

REAL (KIND=ireals), DIMENSION(:,kcm-1:), OPTIONAL, INTENT(INOUT) :: &
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)

! ----------------
! Local variables:
! ----------------

  INTEGER (KIND=iintegers) ::  &
    i,          & !  loop index
    kcp           !  buffer for the vertical index of the upper boudary of the canopy

  REAL    (KIND=ireals   ) ::  fakt

!-------------------------------------------------------------------------------
! Begin Subroutine init_canopy
!-------------------------------------------------------------------------------

  kcp=kcm !save current value of 'kcm', that might have been used for allocation before

  IF (PRESENT(h_can) .AND. PRESENT(hhl)) THEN
   ! h_can is a primary external parameter. The initial values of h_can are 0.
   ! If we don't change this, no canopy will be resolved in the vertical direction.

     kcm=ke
     DO WHILE (MAXVAL( h_can - hhl(:,kcm) + hhl(:,ke+1) ) .GT. z0)
        kcm=kcm-1
     END DO

   ! Up to now kcm points to the lowest layer being not a canopy layer.
   ! From now on kcm points the highest layer being     a canopy layer:

     kcm=kcm+1
  END IF

! Input of the external canopy-parameters:

! At this stage there is no concept of generating those Parameters (c_big, c_sml, r_air).
! They may be derived as functions of rbig, dbig, rsml, dsml, which either come from
! the primary external parameter files or may be derived from other primary external
! parameters like canopy-hight and -type:

! Provisional values for the canopy parameters:
  IF (kcp.LE.kcm) THEN
     !Uppermost canopy level has been determined before, so canopy fields are allocated yet
     IF (PRESENT(c_big)) c_big=z0 !cbig !isotr. drag-coeff. of big canopy-elem.
     IF (PRESENT(c_sml)) c_sml=z0 !csml ! ,,       ,,       ,, small     ,,
     IF (PRESENT(r_air)) r_air=z0 !log(1-rdrg) !log of the volume-fraction being not covered
  END IF

! Provisional values for pattern length array:
  IF (PRESENT(l_pat)) THEN
     DO i=i_stp, i_enp
        IF (fr_land(i) < z1d2) THEN
           l_pat(i)=z0
        ELSE
           IF (PRESENT(d_pat) .AND. imode_pat_len.EQ.2) THEN
              !Restriction of 'pat_len' by 'd_pat':
              l_pat(i)=MIN( pat_len, d_pat(i) )
           ELSE
              l_pat(i)=pat_len !should be a 2D external parameter field
           END IF
           l_pat(i)=l_hori(i)*l_pat(i)/(l_hori(i)+l_pat(i))
        END IF
     END DO
  END IF

! Effective values of the surface area indices:
  IF (PRESENT(sai) .AND. PRESENT(eai)   .AND. PRESENT(tai) .AND. &
      PRESENT(lai) .AND. PRESENT(plcov) .AND. PRESENT(icant) ) THEN

     DO i=i_stp, i_enp
        IF (fr_land(i) < z1d2) THEN
           sai(i)=c_sea
        ELSE
           tai(i)=MAX( 1.0E-6_ireals, lai(i) )
        END IF
     END DO

     IF (icant.EQ.1) THEN
        DO i=i_stp, i_enp
           IF (fr_land(i) >= z1d2) THEN
              sai(i)=tai(i)
              eai(i)=(z1-plcov(i))*sai(i)
              tai(i)=plcov(i)*tai(i)
            END IF
        END DO
     ELSE
        DO i=i_stp, i_enp
           IF (fr_land(i) >= z1d2) THEN
              tai(i)=plcov(i)*tai(i)  ! transpiration area index
              eai(i)=c_soil           ! evaporation area index
              sai(i)=c_lnd+tai(i)     ! surface area index
           END IF
        END DO
        IF (e_surf.NE.z1) THEN
           DO i=i_stp, i_enp
              fakt=EXP( e_surf*LOG( sai(i)) )/sai(i)
            ! Effective area indices by multiplication with the reduction factor fakt:
              sai(i)=fakt*sai(i)
              eai(i)=fakt*eai(i)
              tai(i)=fakt*tai(i)
           END DO
        END IF
     END IF

  END IF

END SUBROUTINE init_canopy

!********************************************************************************
!********************************************************************************

!+ Module procedure organize_turbdiff in "src_turbdiff" for organising the calls
!+ of surface-to-atmosphere trasnsfer and turbulent diffusion:

SUBROUTINE organize_turbdiff ( &
!
          iini, lturatm, ltursrf, lstfnct, lnsfdia, ltkeinp, lgz0inp, &
          itnd, lum_dif, lvm_dif, lscadif, lsrflux, lsfluse, lqvcrst, &
!
          dt_var,dt_tke, nprv,ntur,ntim, &
!
          ie, ke, ke1, kcm, &
!
          i_st, i_en, i_stp, i_enp, &
!
          l_hori, hhl, dp0, trop_mask, innertrop_mask, &
!
          fr_land, depth_lk, h_ice, gz0, sai, &
!
          d_pat, c_big, c_sml, r_air, &
!
          t_g, qv_s, ps, &
          u, v, w, t, qv, qc, prs, rho, epr, &
!
          ptr, ndtr, &
!
          tcm, tch, tvm, tvh, tfm, tfh, tfv, tkr, tkred_sfc, &
          tke, tkvm, tkvh, rcld, tkhm, tkhh, &
          hdef2, hdiv, dwdx, dwdy,           &
!
          edr, tket_sso, tket_conv, tket_hshr, &
          u_tens, v_tens, t_tens, &
          qv_tens, qc_tens, &
          tketens, tketadv, &
          qv_conv, ut_sso, vt_sso, &
!
          t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m, &
          shfl_s, lhfl_s, qvfl_s, umfl_s, vmfl_s, &
!
          ierrstat, errormsg, eroutine)


!-------------------------------------------------------------------------------
! Description:
!
! Organizes the CALL of 'turbtran' and 'turbdiff'
! and 'iinit'.

! Method:
!
! All tendency parameters are OPTIONAL (except 'tketens' in case of "lturatm=T". If they are missing,
!  calculated tendencies of SUB 'turbdiff' are automatically added to the related prognostic variables.
! It is also possible to use only one time level for TKE using "ntim=1" and thus "nprv=1=ntur".

! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8236 1493
!  email:  matthias.raschendorfer@dwd.de

!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!Formal Parameters:
!-------------------------------------------------------------------------------

! 0. Parameters controlling the call of 'organize_turbdiff':

LOGICAL, INTENT(IN) :: &

   lnsfdia,      & !calculation of (synoptical) near-surface variables required
   lsrflux,      & !calculation of surface flux densities in 'trubtran'

   lstfnct,      & !calculation of stability function required

   lturatm,      & !running turbulence model between atmosph. layers (updating diffusion coefficients)
   ltursrf,      & !running turbulence model at the surface layer (updating transfer coefficients
   lum_dif,      & !running vertical gradient diffusion of horizontal u-momenum
   lvm_dif,      & !running vertical gradient diffusion of horizontal v-momenum
   lscadif,      & !running vertical gradient diffusion of scalar properties

   lsfluse,      & !use explicit heat flux densities at the suface
   lqvcrst,      & !qv-flux-divergence reset requested (only if 'qv_conv' is present)

   ltkeinp,      & !TKE present as input (at level k=ke1 for current time level 'ntur')
   lgz0inp         !gz0 present as input

REAL (KIND=ireals), INTENT(IN) :: &

   dt_var,       & !time step for ordinary prognostic variables
   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER (KIND=iintegers), INTENT(IN) :: &

   iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
   itnd,         & !type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
                   !                               2: by adding to current profile before vertical diffusion
                   !                               3: by using corrected virtual vertical profiles
   nprv,         & !previous    time step index of tke
   ntur,         & !current new time step index of tke
   ntim            !number of tke time levels

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------
!
    ie,           & ! number of grid points in zonal      direction
    ke,           & ! index of the lowest main model level
    ke1             ! index of the lowest model half level (=ke+1)

INTEGER (KIND=iintegers), INTENT(INOUT) :: &
!
    kcm             ! level index of the upper canopy bound

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------
!   These variables give the start- and the end-indices of the
!   forecast for the prognostic variables in a horizontal layer.
!
! Horizontal indices:
! --------------------------------------------------------------------------

    i_st,    i_en,    & ! start and end index of mass points
    i_stp,   i_enp      ! start and end index of mass points including model boundary lines

! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------

REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
    dp0             ! pressure thickness of layer                   (pa )

REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
!
    l_hori,       & ! horizontal grid spacing (m)
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area             ( 1 )
    depth_lk,     & ! lake depth                                    ( m )
    sai             ! surface area index                            ( 1 )

REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
!
    d_pat           ! external geometric dimension of circulation patterns

REAL (KIND=ireals), DIMENSION(:,kcm:), TARGET, OPTIONAL, INTENT(IN) :: &
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml           ! effective drag coefficient of canopy elements
                    ! smaller than the turbulent length scale            (1/m)

REAL (KIND=ireals), DIMENSION(:,kcm-1:), TARGET, OPTIONAL, INTENT(IN) :: &
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)

! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------

REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
!
    h_ice           ! ice thickness                                 (  m  )

REAL (KIND=ireals), DIMENSION(:), TARGET, INTENT(IN) :: &
!
    ps,           & ! surface pressure                              ( pa  )
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
!
! Atmospheric model variables:
! ---------------------------------
!
     u,           & ! zonal wind speed       (at mass positions)    ( m/s )
     v,           & ! meridional wind speed  (at mass positions)    ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
     prs            ! atmospheric pressure                          ( pa  )

REAL (KIND=ireals), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(IN) :: &
!
     rho,         & ! total density of air                          (kg/m3)
     epr            ! exner pressure                                 (1)

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
     w              ! vertical wind speed (defined on half levels)  ( m/s )

TYPE (modvar),            OPTIONAL :: ptr(:) ! passive tracers
INTEGER (KIND=iintegers), OPTIONAL :: ndtr   ! number of tracers to be diffused

REAL (KIND=ireals), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------
!
     gz0,          & ! roughness length * g of the vertically not
                     ! resolved canopy                               (m2/s2)
!Achtung: Der g-Faktor ist ueberfluessig!

!    turbulent (transfer) velocity scales at the surface
     tvm,          & ! for momentum                                  ( m/s)
     tvh,          & ! for heat and moisture                         ( m/s)

     !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
     !vecolities  'tvm' and 'tvh' makes live much easier!!

!    turbulent transfer factors for laminar- and roughness-layer transfer
     tfm,          & ! of momentum                                     --
     tfh,          & ! of scalars                                      --
     tfv             ! of water vapor compared to heat                 --

REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL, INTENT(INOUT) :: &

!    turbulent transfer coefficients at the surface
!    (only as optional output for 'ltursrf'):
     tcm,          & ! for momentum                                  ( -- )
     tch,          & ! for scalars (heat and moisture)               ( -- )

!    reference surface diffusion coefficient
!    (only if 'ltursrf' and "imode_trancnf.GE.2"):
     tkr             ! l*Ustar                                       (m2/s)

REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
     tkred_sfc       ! reduction factor for minimum diffusion coefficients near the surface

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=ireals), DIMENSION(ie,ke1,ntim), TARGET, INTENT(INOUT) :: &
!
     tke             ! q:=SQRT(2*TKE); TKE='turbul. kin. energy'     ( m/s )
                     ! (defined on half levels)

REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
!
     tkvm,         & ! turbulent diffusion coefficient for momentum  (m2/s )
     tkvh            ! turbulent diffusion coefficient for heat      (m2/s )
                     ! (and other scalars)

REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
!
     rcld            ! standard deviation of the local oversaturation
                     ! (as input and output)
                     ! fractional cloud cover (in turbdiff)            --

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
!
     hdef2,        & ! horizontal deformation square at half levels  ( 1/s2 )
     hdiv,         & ! horizontal divergence                   ,,    ( 1/s )
!
     dwdx,         & ! zonal      derivative of vertical wind  ,,    ( 1/s )
     dwdy            ! meridional derivative of vertical wind  ,,    ( 1/s )

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
! Tendency fields for the prognostic variables:
! -----------------------------------------------
!
     u_tens,       & ! u-tendency                                    ( m/s2)
     v_tens,       & ! v-tendency                                    ( m/s2)
     t_tens,       & ! t-tendency                                    ( K/s )
     qv_tens,      & ! qv-tendency                                   ( 1/s )
     qc_tens         ! qc-tendency                                   ( 1/s )

REAL (KIND=ireals), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
!
     tketens,      & ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)
     tketadv         ! advection tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
!
     trop_mask,  &   ! mask-factor (1: within tropics; 0: within extra-tropics)
     innertrop_mask  ! used for vertical smoothing of TKE forcing terms
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: &
!
     qv_conv         ! qv-flux-convergence                            ( 1/s )

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
     ut_sso,       & ! u-tendency due to the SSO-Scheme              ( 1/s )
     vt_sso          ! v-tendency due to the SSO-Scheme              ( 1/s )

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(OUT) :: &
!
     edr,          & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
     tket_sso,     & ! TKE-tendency due to SSO wake production       (m2/s3)
     tket_hshr,    & ! TKE-tendency due to separ. horiz. shear       (m2/s3)
!
     tkhm,         & ! horizontal diffusion coefficient for momentum ( m2/s )
     tkhh            ! horizontal diffusion coefficient for scalars  ( m2/s )

REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
     tket_conv       ! TKE-tendency due to convective buoyancy       (m2/s3)

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: &
!
! Diagnostic near surface variables:
! -----------------------------------------------
!
     t_2m,         & ! temperature in 2m                             (  K  )
     qv_2m,        & ! specific water vapor content in 2m            (kg/kg)
     td_2m,        & ! dew-point in 2m                               (  K  )
     rh_2m,        & ! relative humidity in 2m                       (  %  )
     u_10m,        & ! zonal wind in 10m                             ( m/s )
     v_10m           ! meridional wind in 10m                        ( m/s )

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2)    (positive downward)
     lhfl_s,       & ! latent   heat flux at the surface             (W/m2)    (positive downward)
     qvfl_s,       & ! water vapor   flux at the surface             (kg/m2/s) (positive downward)
     umfl_s,       & ! u-momentum flux at the surface                (N/m2)    (positive downward)
     vmfl_s          ! v-momentum flux at the surface                (N/m2)    (positive downward)

INTEGER (KIND=iintegers), INTENT(INOUT) :: ierrstat

CHARACTER (LEN=*), INTENT(INOUT) :: eroutine
CHARACTER (LEN=*), INTENT(INOUT) :: errormsg

!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

INTEGER (KIND=iintegers), PARAMETER :: &
!
!    Indexgrenzen:
!
     nscal=3,     & !Skalare Groessen 1-ter Ordnung
     ninv=2,      & !daraus abgeleitete gegenueber vertikalen
                    !(feuchtadiabatischen) Verrueckungen
                    !invarianten Groessen
     nvel=2,      & !Geschwindigkeitskomponenten
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

!    Beachte:     u_m,v_m muessen in [1,nvel] liegen;
!           aber: tet_l,h2o_g in [nvel+1,nred]
!           und   tem (tet),vap,liq in [nvel+1,nmvar]

INTEGER (KIND=iintegers) :: &
     i, k,       & !horizontaler und vertikaler Laufindex
!
     nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes
     it_start,   & !Startindex der Iterationen
     kem,        & !ke oder ke1
     ntrac,      & !number of included passive tracers
     ndiff,      & !number of 1-st order model variables to be diffused
!
! Transformed horizontal indices:
!
! Zonal direction:
     istart,    iend,    & ! start and end index for the inner model domain
     istartpar, iendpar    ! start and end index including model boundary lines

LOGICAL ::  &
!
     lini,      & !initialization required
     ldovardif, & !berechne (teil-)implizite Vert.diff von Mod.var 1-ter Ordnung
     ldogrdcor, & !mache Gradientkorrektur bei Berechnung der vertikalen Diffusion
     lssintact    !trenne Skalen-Interaktionsterme vom mech. Forcing ab
!
!     Zwischenspeicher fuer rcpv, rcpl und lhocp
!     und Faktor zur Beruecksichtigung von Fluessigwasser:
!

REAL (KIND=ireals) :: &
!
     fr_tke, &      !z1/dt_tke
     zrcpv,zrcpl    !current value of 'rcpv' and 'rcpl' dependent on 'lcpfluc'

REAL (KIND=ireals), TARGET :: &
!
     c_tke,tet_g,rim, &
     c_m,c_h, b_m,b_h,  sm_0, sh_0, &
     d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
     a_3,a_5,a_6

REAL (KIND=ireals), TARGET :: &
!
     l_scal(ie),fc_min(ie), &
!
     tmps(ie,ke1:ke1), &  !surface temperature
     vaps(ie,ke1:ke1), &  !surface specific humidity
     prss(ie,ke1:ke1), &  !surface pressure
     eprs(ie,ke1:ke1), &  !surface Exner-factor
     liqs(ie,ke1:ke1)     !liquid water content at the surface

REAL (KIND=ireals) :: &
!
     grad(nmvar)  !vertical gradient

REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
!    pointer for density, exner factor and eddy dissipation rate:
     rhoh(:,:), exner(:,:), ediss(:,:), &
!
!    pointer for densities of sh-flux and qv-flux:
     shflu_s(:), qvflu_s(:)

REAL (KIND=ireals), TARGET :: &
!
!    memory-target for above pointer:
     rho_tar(ie,ke), exner_tar(ie,ke), diss_tar(ie,ke1), &
     shflu_s_tar(ie), qvflu_s_tar(ie)

!-------------------------------------------------------------------------------

 istat=0; ilocstat=0
 errormsg=''; eroutine='organize_turbdiff'; lerror=.FALSE.

!test: always 5 initial iterations
!it_end=5
!test

!Note:
!Since TKE may not be present at the boundary of the model domain, the turbulence model
! can be applied to the total model domain including the boundary, if 'i_stp', 'i_enp'
! do contain these over all boundary as well.

!All variables and their tendencies are defined at horizontal mass positions.

 IF (itnd.LT.0 .OR. itnd.GT.3) THEN
    ierrstat = 1002
    errormsg= &
    'ERROR *** Parameter ''itnd'' is out of range ***'
     lerror=.TRUE.; RETURN
 END IF
 IF (ltursrf .AND. imode_trancnf.GE.2 .AND. .NOT.PRESENT(tkr)) THEN
    ierrstat = 1002
    errormsg= &
    'ERROR *** Variable ''tkr'' needs to be present ***'
     lerror=.TRUE.; RETURN
 END IF

 ldogrdcor=(lexpcor .AND. lturatm)             !gradient correction has to be done
 ldovardif=(lum_dif .OR. lvm_dif .OR. lscadif) !any variable has to be diffused

 lssintact=((ltkesso.OR.ltkeshs.OR.ltkecon) .AND. imode_adshear.EQ.1)

 IF (PRESENT(ptr)) THEN !passive tracers are present
    !number of tracers
    IF (PRESENT(ndtr)) THEN
       ntrac = ndtr
    ELSE
       ntrac = UBOUND(ptr,1)
    END IF
 ELSE
    ntrac=0
 END IF
 ndiff=nmvar+ntrac !number of 1-st order variables used in the turbulence model
                   !note that cloud ice is treated like a passive trace here
 IF (lcpfluc) THEN
    zrcpv=rcpv
    zrcpl=rcpl
 ELSE
    zrcpv=z0
    zrcpl=z0
 END IF

 fr_tke=z1/dt_tke

 tet_g=grav/cp_d !adiabatic T-gradient

 kem=ke

 istart=i_st    ; iend=i_en
 istartpar=i_stp; iendpar=i_enp

 DO i=istartpar,iendpar
!Achtung: Korrektur durch Faktor 1/2 (wirkt bei sehr kleinen horiz. Gitterzellen
    l_scal(i)=MIN( z1d2*l_hori(i), tur_len )
!__________________________________________________________________________
!test: frm ohne fc_min-Beschraenkung: Bewirkt Unterschiede!
     fc_min(i)=(vel_min/MAX( l_hori(i), tur_len ))**2
!      fc_min(i)=z0
!__________________________________________________________________________
 END DO

 IF (iini.GT.0) THEN !an initialization run
    lini=.TRUE.
    IF (iini.EQ.1) THEN !separate initialization before the time loop
       it_start=1 !only 'it_end' iterations for initialization
                  !and an additional one at the first time loop
    ELSE !initialization within the first time step
       it_start=0 !"it_end+1" iterations for initializatzion
    END IF

!   Bestimmung der initialen Werte fuer die laminaren Reduktions-
!   koeffizienten und die Standardabw. des Saettigungsdef.:
!   Initializing some special variables:

    DO k=1, ke1
       DO i=istartpar,iendpar
          rcld(i,k)=z0 !no standard-deviat. of local over-saturation
       END DO
    END DO
    DO i=istartpar,iendpar
       tfh(i)=z1 !no roughness- and laminar-layer-resistance for scalars
       tfm(i)=z1 !no roughness- and laminar-layer-resistance for momentum
    END DO
    !Notice that the above initalizations will stay, if another turbulence-
    !or transfer-scheme is used!
 ELSE !not an initialization run
    lini=.FALSE.
    it_start=it_end !only one iteration
!test:
!it_start=1
!test
 END IF

 !Note:
 !A call with "iini=2" (at the first time step) provides the same result as a
 !  call with "iini=1" (before the   time loop) followed by a second
 !  call with "iini=0" (at the first time step).

 nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

 !Note:
 !It is also possible to use only one time level for TKE ("ntim=1" and thus "nprv=1=ntur").

 IF (PRESENT(rho)) THEN
    rhoh => rho
 ELSE
    rhoh => rho_tar
  ! ALLOCATE ( rhoh(ie,ke),   STAT=ilocstat ); istat = istat + ilocstat
 END IF

 IF (PRESENT(epr)) THEN
    exner => epr
 ELSE
    exner => exner_tar
  ! ALLOCATE ( exner(ie,ke),  STAT=ilocstat ); istat = istat + ilocstat
 END IF

 IF (PRESENT(edr)) THEN
    ediss => edr
 ELSE
    ediss => diss_tar
  ! ALLOCATE ( ediss(ie,ke1),  STAT=ilocstat ); istat = istat + ilocstat
 END IF
 IF (PRESENT(shfl_s)) THEN
    shflu_s => shfl_s
 ELSE
    shflu_s => shflu_s_tar
  ! ALLOCATE ( shflu_s(ie),    STAT=ilocstat ); istat = istat + ilocstat
 END IF
 IF (PRESENT(qvfl_s)) THEN
    qvflu_s => qvfl_s
 ELSE
    qvflu_s => qvflu_s_tar
  ! ALLOCATE ( qvflu_s(ie),    STAT=ilocstat ); istat = istat + ilocstat
 END IF

!Achtung: Korrektur: Konsistente Behandlung der Null-Fluss-Randbedingung
!als moeglicher 'default' (etwa fuer qc)
 IF (ilow_def_cond.EQ.2) THEN !zero surface value of liquid water
!DIR$ IVDEP
    DO i=istartpar, iendpar
       liqs(i,ke1)=z0
    END DO
 ELSE !constant liquid water within the transfer-layer
!DIR$ IVDEP
    DO i=istartpar, iendpar
       liqs(i,ke1)=qc(i,ke)
    END DO
 END IF

!-------------------------------------------------------------------------------

 IF (ltursrf) THEN
    IF (.NOT.lerror) CALL turbtran
 END IF
 IF (lturatm .OR. ldovardif) THEN
    IF (.NOT.lerror) CALL turbdiff
 END IF

!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------

!********************************************************************************

!+ Module procedure canopy_source in "src_turbdiff" for calculation
!+ of scalar source terms inside the model canopy

SUBROUTINE canopy_source

!_________________________________________________________________________________
!
! Description:
!

! Method:
!
!
! Current Code Owner: DWD,
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8236 1493
!  email:  matthias.raschendorfer@dwd.de
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================

IMPLICIT NONE

END SUBROUTINE canopy_source

!********************************************************************************

!+ Module procedure turb_param in "src_turbdiff" for computing some deduced parameters
!+ for turbulent transfer and - diffusion


SUBROUTINE turb_param

   IMPLICIT NONE

   INTEGER (KIND=iintegers) :: i !loop index

!     Belegung abgeleiteter Parameter:

      c_tke=d_m**z1d3

      c_m=z1-z1/(a_m*c_tke)-z6*a_m/d_m !=3*0.08
      c_h=z0 !kann auch als unabhaengiger Parameter aufgefasst werden

      b_m=z1-c_m
      b_h=z1-c_h

      d_0=d_m

      d_1=z1/a_h
      d_2=z1/a_m
      d_3=z9*a_h
      d_4=z6*a_m
      d_5=z3*(d_h+d_4)
      d_6=d_3+z3*d_4

      rim=z1/(z1+(d_m-d_4)/d_5) !1-Rf_c

      a_3=d_3/(d_2*d_m)        !
      a_5=d_5/(d_1*d_m)
      a_6=d_6/(d_2*d_m)

      sh_0=(b_h-d_4/d_m)/d_1 !stability-function for scalars  at neutr. strat.
      sm_0=(b_m-d_4/d_m)/d_2 !stability-function for momentum at neutr. strat.

END SUBROUTINE turb_param

!********************************************************************************

!+ Module procedure turbtran in "src_turbdiff" for computing the coefficients
!+ for turbulent transfer

SUBROUTINE turbtran

!     INCLUDE 'turbtran.incf'

!     TURBTRAN.incf

! Description:
!
!     Es werden die Transferkoeffizienten fuer den Austausch von Impuls,
!     sowie fuehlbarer und latenter Waerme bestimmt und die Modellwerte
!     fuer die bodennahen Messwerte (in 2m und 10m) berechnet.
!
! Method:
!
!     Hierzu wird der gesamte Bereich von den festen Oberflachen am
!     Unterrand des Modells bis hin zur untersten Hauptflaeche in
!     die drei Teilbereiche:
!
!     - laminare Grenzschicht (L-Schicht)
!     - turbulente Bestandesschicht (B-Schicht)
!     - turbulnte Prandtl-Schicht (P-Schicht)
!
!     aufgeteilt. Fuer jeden dieser Teilbereiche wird (getrennt nach
!     skalaren Eigenschaften und Impuls) ein zugehoeriger Transport-
!     widerstand berechnet, der gleich einer effektiven Widerstands-
!     laenge ( dz_(sg, g0, 0a)_(h,m) ) dividiert durch den Diffusions-
!     koeffizienten am Unterrand der P-Schicht (Niveau '0') ist.
!     Die Konzentrationen am Unterrand der B-Schicht, also im
!     Abstand der L-Schicht-Dicke entlang der festen Oberflaechen,
!     haben den Index 'g' (ground) und die Oberflaechenkonzentrationen
!     den Index 's'. Groessen fuer den Imoulst haben den Index 'm' (momentum)
!     und solche fuer skalare Eigenschaften 'h' (heat).
!     Der Widerstand der P-Schicht vom Nieveau '0' bis zum Niveau 'a'
!     (atmospheric) der untersten Hauptflaeche wird durch vertikale
!     Integration der Modellgleichungen in P-Schicht-Approximation
!     (vertikal konstante Flussdichten, turbulente Laengenskala lin. Funkt.
!      von der Hoehe) gewonnen.
!     Dabei wird das atmosphaerische Turbulenzschema aus der Subroutine
!     'turbdiff' benutzt, so dass alo keine empirischen Profilfunktionen
!     benutzt werden. Zur Vereinfachung der Integration  wird das Produkt
!     aus turbulenter Geschwindigkeitsskala 'q' und der Stabilitaetsfunktion
!     s(h,m), alo die stabilitaetsabhaengige turb. Geschwindigkeitsskala
!     'v' innerhalb der P-Schicht als linear angesehen.
!     Die turb. Laengenskala im Niveau '0' wird mit der Rauhigkeitslaenge 'z0'
!     (multipliziert mit der v.Kaman-Konstanten) gleichgesetzt. Formal werden
!     dann fuer das Nieveau '0' Vertikalgradienten und auch Diffusions-
!     koeffizienten abgeleitet.
!     Unter der Annahme, dass 'v' innerhalb der B-Schicht konstant bleibt,
!     ergibt sich die laminare Widerstandslaenge dz_sg als prop. zu 'z0'
!     und die Widerstandsstrecke durch die B-Schicht als prop. zu
!     'z0*ln(delta/z0)', wobei 'delta' die Dicke der L-Schicht ist, die der
!     Abstand von einer ebenen Wand sein soll in dem der turbulente
!     Diffusionskoeffizient fuer Impuls gleich dem molekularen ist.
!     Ferner wird angenommen, dass die Widerstaende durch die L- und
!     B-Schicht prop. zur effektiven Quellflaech der Bestandeselemente
!     zuzueglich der Grundflaeche des Erdbodens sind. Die Bestandesoberflaechen
!     werden durch den Wert 'sai' (surface area index) ausgedrueckt und setzt
!     sich aus dem Flaechenindex der transpirierenden Oberflaechen 'lai'
!     (leaf area index) und dem fuer die nicht transpirierenden Flaechen
!     zusammen. Im Falle nicht benetzter Oberlfaechen hat die latente Waerme
!     i.a. eine kleinere Quellflaeche als die fuehlbare Waerme, so dass die
!     Wiederstaende fuer beide Groessen unterschieden werden muessten.
!     Um dies zu vermeiden, wird nur der Widerstand fuer die fuehlbare Waerme
!     berechnet. Dafuer wird aber bei der Berechnung der effektiven
!     Oberflaechenkonzentration 'qv_s' der spez. Feuchtigkeit in Subroutine
!     'terra1' dieser Effekt beruecksichtigt.
!     Beim vertikalen Impulstransport ist aber noch die zusaetzliche
!     Impulssenke innerhalb der B-Schicht durch die Wirkung der Formreibungs-
!     kraft zu beruecksichtigen, was durch einen zusaetzlichen Flaechenindex
!     'dai' (drag area index) bewerkstelligt wird.
!
!     Die Vertikalprofile aller Eigenschaften innerhalb der P-Schicht ergeben
!     sich aus dem vertikal integrierten Turbulenzmodell in P-Schicht-
!     Approximation zu logarithmischen Funktionen, welche durch die
!     thermische Schichtung modifiziert sind. Wie bereits erwaehnt, ist die
!     Stabilitaetsfunktion nur noch von Konstanten des atmosphaerischen
!     Turbulenzmodells abhaengig. Das Transferschema ist somit auch automatisch
!     konsistent zum oben anschliessenden Turbulenzmodell formuliert.
!
!     Die Profilfunktionen innerhalb der B-Schicht ergeben sich aus der
!     Annahme eines Gleichgewichtes zwischen vertikalen Flussdichtedivergenzen
!     und Quellstaerken durch die laminaren Grenzschichten der Rauhigkeits-
!     elemente bei vertikal konstanten Bestandeseigenschaften zu exponentiellen
!     Funktionen. Durch die Bedingung eines glatten Ueberganges zwischen beiden
!     Profiltypen im Niveau '0' und der Bedingung, dass im Abstand einer
!     effektiven Bestandesdicke 'Hb' unterhalb des Nieveaus '0' die Bestandes-
!     profile in die Konzentration am Unterrand der B-Schicht (Niveau mit
!     Index 'g') uebergehen, ist das gesamte Transferschema geschlossen und
!     es kann auch der "drag area index" 'dai', sowie die Bestandeshoehe
!     'Hb' selbst eliminiert werden.
!
!     Zur Charakterisierung des Oberflaechentransfers werden dann nur die
!     externen Parameter 'z0', 'sai', 'lai' und je ein globaler Parameter
!     fuer den laminaren Grenzschichtwiderstand des skalaren - und des
!     Impulstransportes benoetigt '(lam_(h,m)'. Hieraus koennte auch eine
!     aequivalente Rauhigkeitslaenge fuer Skalare 'z0h' berechnet werden.
!     Die Oberfalaechenkonzentrationen (Niveau mit Index 's') fuer die skalaren
!     Groessen werden im Modul 'terra' berechnet. Fuer den Impuls gilt die
!     Haftbedingung. Im Grundniveau des atmosphaerischen Modells
!     (Niveau '0') verschwindet also der Wind i.a. nicht; dies ist erst
!     entlang der festen Oberfalechen der Fall. Die bodennahen synoptischen
!     Niveaus werden nun vom Niveau 'z=-Hb', also von der effektiven Bestandes-
!     grundflaeche (Umsatzniveau) aus gezaehlt. Ist z.B. 'Hb>2m', werden
!     die 2m-Werte entlang der exponentiellen Bestandesprofile ausgeweret.
!     Ist 'Hb<2m', wird das logarithmische Profil in der Hoehe '2m-Hb' entlang
!     dinnerhalb der P-Schicht ausgewertet.
!     Die resultierenden Transferkoeffizienten 'tc(h,m)' sind die Kehrwerte
!     des Gesamtwiderstandes von den festen Oberflaechen (Neviau 's') bis
!     zur untersten Modellhauptflaeche (Niveau 'a').
!     Die turbulenten Diffusionskoeffizienten 'tkv(h,m)' fuer den vertikalen
!     Index 'ke1', beziehen sich aber auf den Unterrand des atmosphaerischen
!     Modells (Niveau '0').
!     Mit Hilfe der Felder 'tf(mh)' werden noch Reduktionsfaktoren der
!     Transferkoeffizienten durch die Wirkung der L-Schicht uebertragen.
!     Diese koennen im Modul 'terra' benutzt werden, um ev. das effektive
!     'qv_s' so zu bestimmen, als gaebe es fuer fuehlbare und latente Waerme
!     unterschiedliche Parameter fuer den laminaren Transportwiderstand.
!     Zu beachten ist, dass im Falle eines vertikal vom atmosphaerischen Modell
!     aufgeloesten 'Makrobestandes' (z.B. Bebauung, Wald oder subskalige
!     Orographie) das Transferschema genauso wie im Falle eines nicht
!     aufgeloesten Bestandes angewendet wird. Allerdings beziehen sich die
!     den Bestand des Transferschemas charakterisierenden externen Parameter
!     dann auf den nicht vertikal aufgeloesten verbleibenden 'Mikrobestand',
!     der ev. allein durch niedrigen Bewuchs gebildet wird.
!     Im Transferschema eingearbeitet ist auch dei iterative Bestimmmung der
!     Rauhigkeitslaenge der Meeresoberflaeche gemaess einer modifizierten
!     Charnock-Formel, bei der die Wellenerzeugung bei verschwindenden
!     mittleren Wind mit hilfe der zur TKE ausgedrueckt wird.


! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  matthias.raschendorfer@dwd.de
!
! History of include-file 'turbtran.incf':

! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 Matthias Raschendorfer
!  Initial release
! 1.33       1999/10/14 Matthias Raschendorfer
!  Improved calculation of the vertical gradients at the z0-level, using a
!  simlified integration of the inverse of the turbulent diff. coeff. from
!  this level to the lowest model level.
!  Introduction of an optional correction of the vertical gradients controled
!  by the LOOGECAL parameter lprfcor.
!  Calculation of an seperate roughness parameter z0h for scalars as a
!  function of z0 and the effective surface area index sai.
!  Reformulation of the laminar resistance.
!  Reformulation of the Charnock-Formula using the additional parameter alpha1.
!  Using the redefined parameter c_m (including the factor 3).
!  Inclusion of the diagnosis of 2m- and 10m- values (former routine
!  'synop_diag').  Attention: this quite fundamental reformulation of the
!  surface scheme has not been implementet in turbdiff.incf yet.
! 1.34       1999/12/10 Matthias Raschendorfer
!  Reformulation of the consideration of a roughness canopy in the part 'synop_diag'.
! 1.37       2000/03/24 Matthias Raschendorfer
! Introduction of an canopy resistance for momentum
! 1.39       2000/05/03 Ulrich Schaettler
!  Global variable names lam_m and lam_h changed to rlam_m and rlam_h.
! 2.2        2000/08/18 Matthias Raschendorfer
!  Introduction of the molecular diffusion coefficients as minimal values
!  for the turbulent ones. Calculation of the SQRT of the wind energy in
!  the modified Charnock formula at the lowest model level
!  using the TKE at that level and not (as before) at the lower boundary.
!  Now sai is defined by (sai_before + 1).
! 2.3        2000/11/15 Guenther Doms
!  Some local variable have been redifined to allow for reproducible
!  results on MPP platforms.
! 2.12       2001/11/07 Matthias Raschendorfer
!  Limitation of 'gama' in the calulation of the Prandtl-layer resistance
!  in order to avoid unrealistic deformation of the Prandtl-layer profiles
!  which were the reason of some "jumps" in the 2m-temperature during periods
!  of stabilisation at the surface.
!  Introduction of a minimal roughness length for the initialization of 'gz0'
!  over water surfaces.
! 2.15       2002/03/20 Matthias Raschendorfer
!  Modified interpolation of the 10m wind vector with the help of the roughness length
!  of a typical SYNOP-station (z0m_d).
!  Multiplication of z0m with the laminar correction factor.
! 2.16       2002/03/28 Matthias Raschendorfer
!  Modification of the laminar limit without touching the value of z0m, as the actual
!  formulation became numerical unstable during a test assimilation run.
! 2.17       2002/05/08 Ulrich Schaettler
!  Optimizations for vectorization: split the big loop in Section 4
! 2.18       2002/07/16 Ulrich Schaettler
!  Added a test for variables dz_s0_h, fac_h_2d before a division
!  (because of problems on the NEC)
! 2.19       2002/10/24 Ulrich Schaettler
!  Deleted 2 lines of code that wrongly remained after update 2.17
!  Adaptations to 2-timelevel scheme (use of ntlev)
! 3.6        2003/12/11 Ulrich Schaettler
!  Optimizations for vectorization: modification of 2 DO WHILE loops
!  Modification of an IF...ELSEIF... Statement for vectorization in Section 4f
! 3.7        2004/02/18 Matthias Raschendorfer
!  Introduction of the parameter rat_sea
! 3.14       2005/01/25 Jochen Foerstner
!  Introduced new types of turbulent diffusion parameterizations
!  Replaced SIGN function
!  Adjusted upper boundary of DO-WHILE loops (by NEC)
! 3.16       2005/07/22 Matthias Raschendorfer
!  Some adaptations to the modifications in 'turbdiff.incf'
! 3.18       2006/03/03 Matthias Raschendorfer
!  Introduction of rh_2m and limitation of td_2m
! 3.19       2006/04/25 Jochen Foerstner / Matthias Raschendorfer
!  Application of the lower limit for the roughness length over sea
!  not only during the initialization of 'gz0'
! 3.21       2006/12/04 Dmitrii Mironov
!  Changes to use the FLake model, Ulrich Schaettler
!  Changed interface to subroutine cloud_diag from meteo_utilities
! V3_23        2007/03/30 Matthias Raschendorfer
!  Renaming some variables for better understanding.
!  Introducing the reduction factor for evaporation tfv.
!  Introduction of some output variables for SCLM.
!  Change of some parameter names.
!  Eliminated stab_funct.incf and turb_param.incf
!  (Substitution of file 'stab_funct.incf' by a SUBROUTINE)
!  New file 'statement_functs.incf' containing statement functions
!  Removing the final tfh- and tfm- calculation with respect to the lam. layer,
!  which was the wrong definition for its use in 'turbdiff.icnf'
! V3_25        2007/05/21 Ulrich Schaettler
!  Moved an IF-clause outside DO-Loops in line 1193ff
! V4_3         2008/02/25 Matthias Raschendorfer
!  Changing interpolation onto diagnostic levels,
!  in particular without an exponential canopy profile,
!  but with a dianostic Prandtl layer interpolation even for scalars,
!  using an adopted canopy layer resistance.
!  Calculation of the 3D diagnostic field 'edr' for the eddy dissipotion rate.
!  Changing the treatment of parameter field dd(:,:) to achiev better vectorisation
!  in SUBROUTINE 'stab_funct'.
! V4_8         2009/02/16 Ulrich Schaettler
!  Introduced itype_diag_t2m and "old" 2m temperature as an option
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Removing the horizontal loops for SCLM treatment.
!  Modifications for seaice model: eliminated l_ls_ice, introduced lseaice
!   (Jan-Peter Schulz)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN

!Note: The History of specific ICON-development
!
! Code Description:
! Language: Fortran 90
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================

! Declarations:

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &
!
              i,k,        & !Diskretis.index fuer lamda, phi, sigma-Level
              k1,k2,ks,   & !spezifische Level-Indices
              n,          & !Index fuer diverse Schleifen
              it_durch      !Durchgangsindex der Iterationen

!     Lokale real Variablen:

      REAL (KIND=ireals) :: &
!
!          Hilfsvariablen:
!
           wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,             & !  ,,         ,,     ,,      Faktoren
!
!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 (fh2,fm2):
!
           fh2,fm2, &
!
!     Platzh. fuer horiz. Geschw.-Komponenten, bel. Geschw. und Druck:
!
           vel1,vel2,velo, patm, &
!
!     Platzh. fuer Hoehendifferenzen und Laengaenskalen:
!
           dh,l_turb,lh,lm,z0d,z_surf,len1,len2, &
           dz_s0_m, dz_sa_m, &
           h_2m, h_10m, & !level heights (equal 2m and 10m)
           a_2m, a_10m, & !turbulent distance of 2m- and 10m-level (with respect to diag. roughness)
           a_atm, &       !turbulent distance of the atmosp. level (with respect to mean  roughness)
!
!     sonstiges:
!
           rin_m,rin_h, fr_sd_h, &
           xf,wf

! Local arrays:

      LOGICAL            ::  &
!
        lgz0ini, &   ! initialization of roughness lenght over water and ice
        lo_ice(ie)   ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=ireals), TARGET ::  &
!
!US for a better vectorization (hopefully):
!
        h_top_2d (ie),    & !boundary level height of transfer layer (top  level)
        h_atm_2d (ie),    & !main     level heigth of transfer layer (atm. level)
!
        edh      (ie),    & !reciprocal of a layer depth
!
        z0m_2d   (ie),    & !mean  roughness length
        z0d_2d   (ie),    & !diag. roughness length
        z2m_2d   (ie),    & !height of 2m  level (above the surface)
        z10m_2d  (ie),    & !height of 10m level (above the surface)
!
        hk_2d    (ie),    &
        hk1_2d   (ie),    &
        h_can_2d (ie),    &
!
        rat_m_2d (ie),    &
        rat_h_2d (ie),    &
        fac_h_2d (ie),    &
        fac_m_2d (ie),    &
        frc_2d   (ie),    & !length scale fraction
!
        vel_2d   (ie),    &
        ts_2d    (ie),    &
        qds_2d   (ie),    &
!
        dz_sg_m  (ie), dz_sg_h  (ie),    &
        dz_g0_m  (ie), dz_g0_h  (ie),    &
        dz_0a_m  (ie), dz_0a_h  (ie),    &
        dz_sa_h  (ie),    &
        dz_s0_h  (ie),    &
!
        velmin   (ie),    &
        ratsea   (ie)

     REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
          , CONTIGUOUS &
#endif
!
          :: &
        vel1_2d  (:),      &
        vel2_2d  (:),      &
        ta_2d    (:),      &
        qda_2d   (:),      &
!
        g_tet    (:),      &
        g_vap    (:),      &
        qsat_dT  (:),      &
!
        epr_2d   (:),      &
        rcl_2d   (:),      &
!
        l_tur_z0 (:),      &
!
        tvt      (:,:)

     REAL (KIND=ireals), TARGET :: &
!
        vari(ie,ke-1:ke1,ndim), & !conserved Prandtl-layer variables
!
           a(ie,ke1:ke1, ndim), & !special thermodynamic variables
!
        rclc(ie,ke1:ke1),       & !cloud cover
   len_scale(ie,ke1:ke1),       & !turbulent master length scale
 tketens_tar(ie,ke1:ke1)          !target for turbulent transport of SQRT(TKE)

     INTEGER (KIND=iintegers) ::  &
        k_2d     (ie)

     REAL (KIND=ireals) :: &
        frm(ie,ke1:ke1),    &
        frh(ie,ke1:ke1),    &
        rhon(ie,ke1:ke1)

!---- End of header ---------------------------------------------------------

! 1)  Vorbereitungen:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbtran'; lerror=.FALSE.

!     Unterste Hauptflaeche halbiert den Abstand zur
!     untersten Nebenflaeche:

      xf=z2
      wf=xf-z1
      xf=z1/xf

      IF (istat /= 0) THEN
         ierrstat = 1004
         errormsg= &
         'ERROR *** Allocation of space for meteofields failed ***'
         lerror=.TRUE.; RETURN
      ENDIF

      IF (PRESENT(tketens)) THEN
         tvt => tketens
      ELSE
         tvt => tketens_tar
         tvt=z0
      END IF

      vel1_2d => vari(:,ke ,u_m)
      vel2_2d => vari(:,ke ,v_m)
      ta_2d   => vari(:,ke ,tet_l)
      qda_2d  => vari(:,ke ,h2o_g)

      qsat_dT => a(:,ke1,3)
      g_tet   => a(:,ke1,4)
      g_vap   => a(:,ke1,5)

      epr_2d  => eprs(:,ke1)
      rcl_2d  => rclc(:,ke1)

      l_tur_z0 => len_scale(:,ke1)

!     Berechnung abgeleiteter Parameter:

      CALL turb_param

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren:

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

!DIR$ IVDEP
      DO i=istartpar, iendpar

         IF (fr_land(i) < z1d2) THEN
            ! Water point.
            IF (.NOT. lseaice) THEN
              ! Sea ice model is not used.
              ! Ice surface if SST is less than the salt water freezing temperature.
              lo_ice(i) = t_g(i) < t0_melt + zt_ice
            ELSE
              ! Sea ice model is used.
              ! Ice surface if ice is present.
              lo_ice(i) = h_ice(i) > z0
            END IF
            IF (llake) THEN
              ! Lake model is used.
              ! Ice surface if this is a lake point AND ice is present.
              IF ((depth_lk(i) > z0) .AND. (h_ice(i) >= h_Ice_min_flk)) &
              lo_ice(i) = .TRUE.
            END IF
         END IF

      END DO

! 3)  Berechnung einiger Hilfsgroessen und Initialisierung der Diffusionskoeff.:

      ! Hoehe des 2m- und 10m-Niveaus:
      h_2m  = z2 ; h_10m = z10

!DIR$ IVDEP
      DO i=istartpar, iendpar
         ! Dicke der Modell-Prandtl-Schicht
         h_top_2d(i) = hhl(i,ke)-hhl(i,ke1)
         h_atm_2d(i) = h_top_2d(i)*xf

         ! Surface-Exner-pressure:
         epr_2d(i) = zexner(ps(i))

         ! Scaling velocity = Wind at lowest full level:
         vel_2d(i) = MAX( vel_min, SQRT(u(i,ke)**2+v(i,ke)**2) )
      END DO

      IF (lprfcor) THEN
         ks=ke-1
      ELSE
         ks=ke
      END IF

      IF (.NOT.PRESENT(epr)) THEN
         DO k=ks, ke
!DIR$ IVDEP
            DO i=istartpar, iendpar
               exner(i,k)=zexner(prs(i,k))
            END DO
         END DO
      END IF

!DIR$ IVDEP
      DO i=istartpar, iendpar
         ratsea(i) = rat_sea
      END DO

      IF (imode_rat_sea.EQ.2) THEN
!<Tuning
!DIR$ IVDEP
         DO i=istartpar, iendpar
            IF (t_g(i) - t(i,ke) > 8._ireals) THEN
               ! Increase rat_sea for very large temperature differences between water and adjacent air
               ! in order to reduce excessive peaks in latent heat flux

               ratsea(i) = rat_sea*(1._ireals + 0.05_ireals*(t_g(i) - t(i,ke) - 8._ireals))
            END IF
         END DO
!>Tuning: This kind of correction can be substituded by a less ad-hoc approach.
      END IF

      IF (imode_vel_min.EQ.2) THEN
!<Tuning
!DIR$ IVDEP
         DO i=istartpar, iendpar
            ! stability-dependent minimum velocity serving as lower limit on surface TKE
            ! (parameterizes small-scale circulations developing over a strongly heated surface;
            ! tuned to get 1 m/s when the land surface is about 10 K warmer than the air in the
            ! lowest model level; nothing is set over water because this turned out to induce
            ! detrimental effects in NH winter)

            velmin(i) = MAX( vel_min, fr_land(i)*(t_g(i)/epr_2d(i) - t(i,ke)/exner(i,ke))/ &
                        LOG(2.e3_ireals*h_atm_2d(i)) )
         END DO
!>Tuning: his kind of correction can be substituded by a less ad-hoc approach.
      END IF

      IF (lini) THEN !only for initialization

         DO i=istartpar, iendpar

            lgz0ini=(.NOT.lgz0inp .AND. fr_land(i) < z1d2)

            IF (imode_trancnf.GE.2 .OR. (lgz0ini .AND. .NOT.lo_ice(i))) THEN
!              Einfachste Schaetzung der Schubspannung als Impusls-
!              flussdichte durch die Nebenflaeche ke mit Hilfe
!              einer diagnostischen TKE ohne Beruecksichtigung von
!              Feuchte-Effekten und mit neuchtralen Stabilitaets-
!              funktion:

               l_turb=h_top_2d(i) !approx. turb. length scale at level ke

!test: different turb. length scale
               l_turb=akt*MAX( len_min, l_turb/( z1+l_turb/l_scal(i) ) )
!  l_turb=akt*MIN( l_scal(i), hhl(i,ke)-hhl(i,ke1) )
!test

               dh=z1d2*(hhl(i,ke-1)-hhl(i,ke1))

               vel1=u(i,ke-1)
               vel2=u(i,ke  )
               grad(u_m)=(vel1-vel2)/dh

               vel1=v(i,ke-1)
               vel2=v(i,ke  )
               grad(v_m)=(vel1-vel2)/dh

               grad(tet_l)=(t(i,ke-1)-t(i,ke))/dh + tet_g

               fm2=MAX( grad(u_m)**2+grad(v_m)**2, fc_min(i) )
               fh2=grav*grad(tet_l)/t(i,ke)

               ! Vereinfachte Loesung mit Rf=Ri:
               IF (fh2.GE.(z1-rim)*fm2) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=l_turb*(sm_0-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2/(fm2-fh2)
                  lm=l_turb*(sm_0-(a_6+a_3)*fakt)
                  lh=l_turb*(sh_0-a_5*fakt)
               END IF

               val1=lm*fm2; val2=lh*fh2
               wert=MAX( val1-val2, rim*val1 )

               IF (ltkeinp) THEN
                  tke(i,ke,nvor)=tke(i,ke,ntur)
               ELSE
                  tke(i,ke,nvor)=SQRT(d_m*l_turb*wert)
               END IF

               val1=MAX ( con_m, tkmmin ); tkvm(i,ke)=lm*tke(i,ke,nvor)
               val2=MAX ( con_h, tkhmin ); tkvh(i,ke)=lh*tke(i,ke,nvor)

               IF (imode_tkemini.EQ.2) THEN
                  tke(i,ke,1)=tke(i,ke,1)*MAX( z1, val1/tkvm(i,ke), & !adapted tke
                                                   val2/tkvh(i,ke) )
               END IF
               tkvm(i,ke)=MAX(val1, tkvm(i,ke)) !corrected tkv
               tkvh(i,ke)=MAX(val2, tkvh(i,ke))

!Achtung: Die 'epsi'-Beschraenkung ist recht willkuerlich und fehlt in COSMO-Version!
               val2=MAX( epsi, tkvm(i,ke)*SQRT(fm2) ) !estimate of Ustar**2
               val1=SQRT(val2) !Ustar
            END IF

            IF (lgz0ini) THEN
              ! Iinitialization of roughness length for water or ice covered surface:

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice:

               IF ( lo_ice(i) ) THEN ! ice covered surface
                  gz0(i)=grav*z0_ice
               ELSE !water covered surface
!                 Schubspannung abhaengiger Wellenhoehe:
                  IF (imode_charpar.EQ.1) THEN !constant Charnock-Parameter
                     fakt=alpha0
                  ELSE IF (depth_lk(i) > z0) THEN
                     ! enhanced Charnock parameter over lakes, parameterizing a non-equlibrium wave spectrum
                     fakt = 0.1_ireals
                  ELSE
                     fakt=alpha0_char(vel_2d(i))
                     !Note: The argument of alpha0_char should be 'vel_10m', which might not yet be
                     !      present during initialization!
                  END IF
                  gz0(i)=MAX( grav*len_min, fakt*val2+alpha1*grav*con_m/val1 )
               END IF
            END IF

            IF (imode_trancnf.GE.2) THEN !new version of init. using estimated Ustar
               tkr(i)=l_turb*val1                          !l_0*Ustar
               rat_m_2d(i)= tkr(i)/tkvm(i,ke)              !Ustar/(q*Sm)_p
               rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke)*sm_0) !Ustar/(q*Sh)_p*Sh(0)/Sm(0)
            END IF

         END DO

      END IF !only for initialization

!DIR$ IVDEP
      DO i=istartpar, iendpar
         z0m_2d(i)  = gz0(i)/grav   !mean roughness length
         l_tur_z0(i)= akt*z0m_2d(i) !turbulent length scale
         frc_2d(i)  = z0m_2d(i)/(h_top_2d(i)+z0m_2d(i)) !length scale fraction
      END DO

      IF (lini) THEN

         IF (imode_trancnf.GE.2) THEN !new version of init. using estimated Ustar
            DO i=istartpar, iendpar
               tkvm(i,ke1)=tkvm(i,ke)*frc_2d(i)*(frc_2d(i)+(z1-frc_2d(i))*rat_m_2d(i))
               tkvh(i,ke1)=tkvh(i,ke)*frc_2d(i)*(frc_2d(i)+(z1-frc_2d(i))*rat_h_2d(i))
            END DO
         ELSE
            DO i=istartpar, iendpar
               tkvm(i,ke) =con_m; tkvh(i,ke) =con_h
               tkvm(i,ke1)=con_m; tkvh(i,ke1)=con_h
            END DO
         END IF

      ELSEIF (imode_trancnf.GE.4) THEN
         !Not for initialization, but for calculation the profile-factors
         !without an upper node based on previous values of transfer-velocity
         !and diffusion-coefficients (at the top of the roughness layer):

         DO i=istartpar, iendpar
            rat_m_2d(i)= tkr(i)/tkvm(i,ke1)              !Ustar/(q*Sm)_0
            rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(Sh(0)/Sm(0))
            tkr(i)=rat_m_2d(i) !saved profile factor for momentum
                               !to be optionally smoothed during direct iteration
         END DO

      END IF

      IF (imode_trancnf.EQ.2 .OR. imode_trancnf.EQ.3) THEN
         !Profile-factors by using the previous diffusion coefficients
         !without a laminar correction, but still based on the upper node:
!DIR$ IVDEP
         DO i=istartpar, iendpar
            rat_m_2d(i)=frc_2d(i)*tkvm(i,ke)/tkvm(i,ke1) !(q*Sm)_p/(q*Sm)_0
            rat_h_2d(i)=frc_2d(i)*tkvh(i,ke)/tkvh(i,ke1) !(q*Sh)_p/(q*Sh)_0
         END DO
      END IF

! 4)  Berechnung der Transferkoeffizienten:

      DO it_durch=it_start, it_end !Iterationen
!print *,"it_durch=",it_durch

!DIR$ IVDEP
         DO i=istartpar, iendpar

            z_surf= z0m_2d(i)/sai(i) !effektive Rauhigkeitslaenge

!           Laminare Korrektur der Diffusionskoeffizienten:
            tkvm(i,ke1)=MAX( con_m, tkvm(i,ke1) )
            tkvh(i,ke1)=MAX( con_h, tkvh(i,ke1) )

            fakt=z1+(z1-REAL(NINT(fr_land(i)),ireals))*(ratsea(i)-z1)

            rin_m=tkvm(i,ke1)/con_m
            rin_h=tkvh(i,ke1)/con_h
!rin_m=tkvm(i,ke1)/con_m+z1
!rin_m=tkvh(i,ke1)/con_h+z1

!           Effektiven Widerstandslaengen der Rauhigkeits-Schicht:

            dz_sg_m(i)=rlam_mom*z_surf
            dz_sg_h(i)=fakt*rlam_heat*z_surf*(rin_h/rin_m)

          ! ohne lam. Grenzschicht fuer Skalare:
            dz_g0_h(i)=z_surf*LOG(rin_m)

          ! inclusive lam. Grenzschicht fuer Skalare:
            dz_s0_h(i)=dz_sg_h(i)+dz_g0_h(i)

         END DO

!DIR$ IVDEP
         DO i=istartpar, iendpar

!           Effektiven Bestandeshoehe:
            IF (dz_sg_h(i).EQ.z0) THEN
               h_can_2d(i)=rat_can*z0m_2d(i)
            ELSE
               h_can_2d(i)=rat_can*dz_s0_h(i)*LOG(dz_s0_h(i)/dz_sg_h(i))
            END IF
         END DO

!DIR$ IVDEP
         DO i=istartpar, iendpar

          ! inclusive lam. Grenzschicht fuer Impuls:
            wert=z1d2*dz_sg_m(i)
            dz_s0_m=wert+SQRT(wert**2+h_can_2d(i)*dz_sg_m(i))

          ! ohne lam. Grenzschicht fuer Impuls:
            dz_g0_m(i)=dz_s0_m-dz_sg_m(i)

!------------------------------------------------------------------------------

!           Profilfakoren der turbulenten Prandtl-Schicht:

            fakt=z0m_2d(i)/h_top_2d(i)

            IF (imode_trancnf.LT.4) THEN
               !Profile-factors by employing previous values of the diffusion-coefficients
               !at the top fo the roughness-layer (0) and also at the upper bound of the
               !lowest atm. model layer (p) as an upper node:


               IF (imode_trancnf.EQ.1) THEN !first version
                  !Profile factors by using the previous diffusion coefficients
                  !including a laminar correction:
                  rat_m_2d(i)=frc_2d(i)*tkvm(i,ke)/tkvm(i,ke1)
                  rat_h_2d(i)=frc_2d(i)*tkvh(i,ke)/tkvh(i,ke1)
               END IF

!Achtung: Die Beschraenkung ist recht willkuerlich
               rat_m_2d(i)=MIN( z2, MAX( z1d2, rat_m_2d(i) ) ) !limitted (q*Sm)_p/(q*Sm)_0
               rat_h_2d(i)=MIN( z2, MAX( z1d2, rat_h_2d(i) ) ) !limitted (q*Sh)_p/(q*Sh)_0

               fac_m_2d(i)=(rat_m_2d(i)-z1)*fakt !non-stab. profile-factor for momentum
               fac_h_2d(i)=(rat_h_2d(i)-z1)*fakt !non-stab. profile-factor for scalars

            ELSE !Profile-factors without using the upper node
               fac_m_2d(i)=z1-rat_m_2d(i) !profile-factor for momentum
               fac_h_2d(i)=z1-rat_h_2d(i) !profile-factor for scalars
!Achtung: Analoge Beschraenkung wie bei "imode_trancnf.LT.4":
!fac_m_2d(i)=Min( fakt, MAX( -z1d2*fakt, fac_m_2d(i) ) )
!fac_h_2d(i)=Min( fakt, MAX( -z1d2*fakt, fac_h_2d(i) ) )
!Achtung: Neutrale Profile:
!fac_m_2d(i)=z0
!fac_h_2d(i)=z0
            END IF

         END DO

!DIR$ IVDEP
         DO i=istartpar, iendpar

            a_atm = h_atm_2d(i)+z0m_2d(i) !turbulente Distanz auf der untersten Hauptflaeche

!           Effektive Widerstandslaengen der turb. Prandtl-Schicht:

            IF (fac_m_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
               !non-stable strat. or based on linear interpolation of profile-function
               !for the velocity scale (q*Sm):
               IF (fac_m_2d(i).EQ.z1) THEN
                  dz_0a_m(i)=z0m_2d(i)*h_atm_2d(i)/a_atm
               ELSE
                  dz_0a_m(i)=z0m_2d(i)*LOG(a_atm/(z0m_2d(i)+fac_m_2d(i)*h_atm_2d(i))) &
                                      /(z1-fac_m_2d(i))
               END IF
            ELSE !based on hyperbolic interpolation of (q*Sm) for stable stratification
               dz_0a_m(i)=z0m_2d(i)*(LOG(a_atm/z0m_2d(i))-fac_m_2d(i)*h_atm_2d(i)/z0m_2d(i)) &
                                   /(z1-fac_m_2d(i))
            END IF
            IF (fac_h_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
               !non-stable strat. or using only linear interpolation of profile-function
               !for the velocity scale (q*Sh):
               IF (fac_h_2d(i).EQ.z1) THEN
                  dz_0a_h(i)=z0m_2d(i)*h_atm_2d(i)/a_atm
               ELSE
                  dz_0a_h(i)=z0m_2d(i)*LOG(a_atm/(z0m_2d(i)+fac_h_2d(i)*h_atm_2d(i))) &
                                      /(z1-fac_h_2d(i))
               END IF
            ELSE !hyperbolic interpolation of (q*Sh) for stable stratification
               dz_0a_h(i)=z0m_2d(i)*(LOG(a_atm/z0m_2d(i))-fac_h_2d(i)*h_atm_2d(i)/z0m_2d(i)) &
                                   /(z1-fac_h_2d(i))
            END IF

         END DO

!DIR$ IVDEP
         DO i=istartpar, iendpar

!           Effektive Widerstandslaengen von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):

            dz_sa_m    = dz_s0_m    + dz_0a_m(i)
            dz_sa_h(i) = dz_s0_h(i) + dz_0a_h(i)

!           Reduktionsfaktoren fuer die Bestandesschicht incl. lam. Grenzschicht:

            tfm(i)=dz_0a_m(i)/dz_sa_m
            tfh(i)=dz_0a_h(i)/dz_sa_h(i)

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den Faktor 'rat_lam'
!           gegenueber fuehlbarer Waerme vergroesserten laminaren Transpostwiderstandes:

            tfv(i)=z1/(z1+(rat_lam-z1)*dz_sg_h(i)/dz_sa_h(i))
         END DO

!        Berechnung der Erhaltungsgroessen in der Prandtl-Schicht:

!Achtung: <Korrketur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
         IF (icldm_tran.EQ.-1 .OR. ilow_def_cond.EQ.2) THEN
            !conserved values at the rigid surface are temperature and humidity
!DIR$ IVDEP
            DO i=istartpar, iendpar
               ts_2d(i)=t_g(i); qds_2d(i)=qv_s(i)
            END DO
         ELSE !conserved variables at the rigid surface depend on liquid water
!DIR$ IVDEP
            DO i=istartpar, iendpar
                ts_2d(i)=(t_g(i) - lhocp*liqs(i,ke1))/eprs(i,ke1)
               qds_2d(i)=qv_s(i) +       liqs(i,ke1)
            END DO
         END IF

         DO k=ks, ke
            DO i=istartpar, iendpar
               vari(i,k,u_m)=u(i,k)
               vari(i,k,v_m)=v(i,k)
            END DO
         END DO
!> Korrektur

         IF (icldm_tran.EQ.-1) THEN !no water phase change possible
            DO k=ks, ke
!DIR$ IVDEP
               DO i=istartpar, iendpar
                  vari(i,k,tet_l)=t(i,k)/exner(i,k)
                  vari(i,k,h2o_g)=qv(i,k)
               END DO
            END DO
         ELSE !water phase changes are possible
            DO k=ks, ke
!DIR$ IVDEP
               DO i=istartpar, iendpar
                  vari(i,k,tet_l)=(t(i,k) - lhocp*qc(i,k))/exner(i,k)
                  vari(i,k,h2o_g)=qv(i,k) +       qc(i,k)
               END DO
            END DO
         END IF

         IF (lprfcor) THEN
!DIR$ IVDEP
            DO i=istartpar, iendpar
               len1=z2*h_top_2d(i)
               len2=(h_top_2d(i)-h_atm_2d(i))**2 &
                   /((hhl(i,ke-1)+hhl(i,ke))*z1d2-hhl(i,ke1)-h_atm_2d(i))
               lm=len1-tfm(i)*h_atm_2d(i)-len2
               lh=len1-tfh(i)*h_atm_2d(i)-len2

               vari(i,ke,u_m  )=(len1*vari(i,ke  ,u_m  ) &
                                  -len2*vari(i,ke-1,u_m  ))/lm
               vari(i,ke,v_m  )=(len1*vari(i,ke  ,v_m  ) &
                                  -len2*vari(i,ke-1,v_m  ))/lm
               vari(i,ke,tet_l)=(len1*vari(i,ke  ,tet_l)-h_atm_2d(i)*tfh(i)* t_g(i)/epr_2d(i) &
                                  -len2*vari(i,ke-1,tet_l))/lh
               vari(i,ke,h2o_g)=(len1*vari(i,ke  ,h2o_g)-h_atm_2d(i)*tfh(i)*qv_s(i) &
                                 -len2*vari(i,ke-1,h2o_g))/lh
            END DO
         END IF

!        Thermodynamische Hilfsvariablen auf dem Unterrand der Prandtl-Schicht:

!DIR$ IVDEP
         DO i=istartpar, iendpar
            prss(i,ke1)=ps(i)
            tmps(i,ke1)=t_g(i)
            vaps(i,ke1)=qv_s(i)
         END DO

         CALL adjust_satur_equil ( khi=ke1, ktp=ke, &
!
              i_st=istartpar, i_en=iendpar, k_st=ke1, k_en=ke1,        &
!
              lcalrho=.TRUE.,  lcalepr=.FALSE., lcaltdv=.TRUE.,        &
              lpotinp=.FALSE., ladjout=.FALSE.,                        &
!
              icldmod=icldm_tran,                                      &
!
              zrcpv=zrcpv, zrcpl=zrcpl,                                &
!
!Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
!         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
              prs=prss, t=tmps, qv=vaps, qc=liqs,                      &
!
              fip=tfh,                                                 &
!
              exner=eprs, rcld=rcld(:,ke1:ke1), dens=rhon,             &
!
              qst_t=a(:,:,3), g_tet=a(:,:,4), g_h2o=a(:,:,5),          &
!
              tet_l=vari(:,ke:ke1,tet_l), q_h2o=vari(:,ke:ke1,h2o_g),  &
                                          q_liq=vari(:,ke:ke1,liq) )

!        Beachte:
!        'vari(:,ke1,tet_l)' und 'vari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen
!        am Unterrand der Prandtl-Schicht, waehrend  ta_2d' => 'vari(:,ke,tet_l)
!        und 'qda_2d  => 'vari(:,ke,h2o_g)' auf diese Groessen bzgl. der untersten
!        Hauptflaeche zeigen, welche zur Interpolation der Groessen an der Oberflaeche
!        auf jenen Unterrand der Prandtl-Schicht (Oberrand der Rauhigkeitsschicht)
!        benutzt werden.

!        Berechnung der benoetigten Vertikalgradienten und der TKE-Antriebe:

         !Vertikalgradienten des Horizontalwindes:
!DIR$ IVDEP
         DO i=istartpar, iendpar
            edh(i)=tfm(i)/dz_0a_m(i)
         END DO
         DO n=1, nvel
!DIR$ IVDEP
            DO i=istartpar, iendpar
               vari(i,ke1,n)=vari(i,ke,n)*edh(i)
            END DO
            !Beachte: Dies ist die Darstellung ohne Nutzung der unteren Randwerte der Prandtl-Schicht
         END DO

         !Scherungs-Antrieb der TKE:
         IF (itype_sher.EQ.2 .AND. PRESENT(dwdx) .AND. PRESENT(dwdy)) THEN
            !Einschliesslich der 3D-Korrektur durch den Vertikalwind bzgl. der mittleren Hangneigung:
!DIR$ IVDEP
            DO i=istartpar, iendpar
               frm(i,ke1)=MAX( (vari(i,ke1,u_m)+dwdx(i,ke1)*edh(i))**2 &
                              +(vari(i,ke1,v_m)+dwdy(i,ke1)*edh(i))**2 &
                              +hdef2(i,ke1)*edh(i)**2, fc_min(i) )
            END DO
            !Beachte: dwdx(ke1), dwdy(ke1) und hdef2(ke1) beziehen sich auf die vorlaeufige Schichtdicke 1m.
         ELSE
!DIR$ IVDEP
            DO i=istartpar, iendpar
               frm(i,ke1)=MAX( vari(i,ke1,u_m)**2+vari(i,ke1,v_m)**2, fc_min(i) )
            END DO
         END IF

         !Vertikalgradienten der dynamisch wirksamen Skalare:
!DIR$ IVDEP
         DO i=istartpar, iendpar
            edh(i)=z1/dz_0a_h(i)
         END DO
         DO n=nvel+1,nred
!DIR$ IVDEP
            DO i=istartpar, iendpar
               vari(i,ke1,n)=(vari(i,ke,n)-vari(i,ke1,n))*edh(i)
            END DO
         END DO
         !'vari(:,ke1,n)' enthaelt jetzt die Vertikalgradienten der Erhaltungsvariablen

         IF (imode_trancnf.EQ.1) THEN !old version of zero-level-gradients requested
            !Transformation of Tet_l-gradient into the old form following from interpolation
            !onto the zero-level in terms of T_l (rather than Tet_l) and correcting the
            !calculated T_l-Gradient by the adiabatic laps rate:
            k=ke1
            DO i=istartpar, iendpar
               wert=vari(i,k,liq)                                 !liquid water at zero-level
               val1=vari(i,k,h2o_g)-vari(i,k,liq)                 !spec. humid. at zero-level
               val2=eprs(i,k)*vari(i,k,tet_l)+lhocp*vari(i,k,liq) !temperature  at zero-level

               vari(i,k,tet_l) = vari(i,k,tet_l)  &
                               + ( (exner(i,k-1)-eprs(i,k))*vari(i,k-1,tet_l)*tfh(i)*edh(i) &
                               + tet_g*(z1-lhocp*wert/val2)/(z1+rvd_m_o*val1-wert) )/eprs(i,k)
            END DO
         END IF

         !Auftriebs-Antrieb der TKE:
!DIR$ IVDEP
         DO i=istartpar, iendpar
            frh(i,ke1)=g_tet(i)*vari(i,ke1,tet_l)+g_vap(i)*vari(i,ke1,h2o_g)
         END DO

!        Berechnung der Stabilitaetslaengen:

         IF (it_durch.EQ.it_start .AND. lini) THEN !Startinitialisierung

            DO i=istartpar, iendpar
               IF (frh(i,ke1).GE.(z1-rim)*frm(i,ke1)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  tkvm(i,ke1)=l_tur_z0(i)*(sm_0-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=tkvm(i,ke1)
               ELSE
                  fakt=frh(i,ke1)/(frm(i,ke1)-frh(i,ke1))
                  tkvm(i,ke1)=l_tur_z0(i)*(sm_0-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=l_tur_z0(i)*(sh_0-a_5*fakt)
               END IF

               val1=tkvm(i,ke1)*frm(i,ke1)
               val2=tkvh(i,ke1)*frh(i,ke1)
               wert=MAX( val1-val2, rim*val1 )

               IF (.NOT.ltkeinp) THEN !TKE not present as input
                  tke(i,ke1,nvor)=MAX( SQRT(d_m*l_tur_z0(i)*wert), vel_min )
               END IF

               IF (imode_tkemini.EQ.2) THEN
                  val1=con_m; val2=con_h
                  tke(i,ke1,nvor)=MAX( tke(i,ke1,nvor), val1/tkvm(i,ke1), &
                                                        val2/tkvh(i,ke1) )
               END IF

            END DO

         ELSE ! mit Hilfe der vorhergehenden TKE-Werte

!DIR$ IVDEP
            DO i=istartpar, iendpar
               wert=z1/tke(i,ke1,nvor)
               tkvm(i,ke1)=tkvm(i,ke1)*wert
               tkvh(i,ke1)=tkvh(i,ke1)*wert
            END DO

         END IF


! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes:

         CALL solve_turb_budgets( khi=ke1, it_s=it_durch, &
                                  i_st=istartpar, i_en=iendpar, k_st=ke1, k_en=ke1,   &
!
                                  lssintact=lssintact, lupfrclim=(imode_trancnf.EQ.1),&
!
                                  imode_stke=imode_tran, imode_vel_min=imode_vel_min, &
!
                                  fm2=frm, fh2=frh, ft2=frm,                          &
                                  lsm=tkvm(:,ke1:ke1), lsh=tkvh(:,ke1:ke1),           &
#ifdef SCLM
                                  grd=vari(:,ke1:ke1,:),                              &
#endif
                                  tls=len_scale(:,ke1:ke1), tvt=tvt(:,ke1:ke1),       &
                                  velmin=velmin(:)                                      )

!DIR$ IVDEP
         DO i=istartpar, iendpar
! 4h)       Bestimmung der durch Wirkung der L-Schicht
!           korrigierten Diffusionskoeffizienten
!           und der zugehoerigen Transferkoeffizienten:

!           Unkorrigierte Diffusionskoeffizienten:

            val1=con_m; tkvm(i,ke1)=tke(i,ke1,ntur)*tkvm(i,ke1)
            val2=con_h; tkvh(i,ke1)=tke(i,ke1,ntur)*tkvh(i,ke1)

            IF (imode_lamdiff.EQ.2) THEN !surface-layer diffusion coeff. always at least at laminar value
               IF (imode_tkemini.EQ.2) THEN
                   tke(i,ke1,ntur)=tke(i,ke1,ntur)*MAX( z1, val1/tkvm(i,ke1), &
                                                            val2/tkvh(i,ke1) )
               END IF
               tkvm(i,ke1)=MAX( val1, tkvm(i,ke1) )
               tkvh(i,ke1)=MAX( val2, tkvh(i,ke1) )
            END IF

!           Belegung der Felder fuer die Transferkoeffizienten:

!Achtung: <zum Vergleich mit alter Variante:
!tvm(i)=tkvm(i,ke1)*tfm(i)/(dz_0a_m(i)*vel_2d(i))
!tvh(i)=tkvh(i,ke1)*tfh(i)/(dz_0a_h(i)*vel_2d(i))
!Achtung: Modifikation tcm -> tvm; tch -> tvh: macht Unterschiede
            tvm(i)=tkvm(i,ke1)*tfm(i)/dz_0a_m(i)
            tvh(i)=tkvh(i,ke1)*tfh(i)/dz_0a_h(i)
         END DO

         IF (PRESENT(tcm) .AND. PRESENT(tch)) THEN
            DO i=istartpar, iendpar
!tcm(i)=tvm(i); tvm(i)=tcm(i)*vel_2d(i)
!tch(i)=tvh(i); tvh(i)=tch(i)*vel_2d(i)
               tcm(i)=tvm(i)/vel_2d(i)
               tch(i)=tvh(i)/vel_2d(i)
!>zum Vergleich
            END DO
         END IF

         IF (imode_trancnf.GE.4 .OR. (imode_trancnf.GE.2 .AND. it_durch.LT.it_end)) THEN
            DO i=istartpar, iendpar
               wert=l_tur_z0(i)*SQRT(tkvm(i,ke1)*SQRT(frm(i,ke1))) !updated l_0*Ustar
               IF (ditsmot.GT.z0) THEN
                  tkr(i)=ditsmot*tkr(i)*tkvm(i,ke1) + (z1-ditsmot)*wert
               ELSE
                  tkr(i)=wert
               END IF
            END DO
         END IF


!DIR$ IVDEP
         DO i=istartpar, iendpar
! 4i)       Einschraenkung von z0m_dia ueber Land:

            IF (fr_land(i) <= z1d2) THEN
               !Ueber See gibt es keinen synoptischen Garten
               z0d_2d(i)=z0m_2d(i)
            ELSE
               !Die Rauhigkeitslaenge einer SYNOP Station soll immer
               !kleiner als 10m bleiben:
               z0d_2d(i)=MIN( h_10m, z0m_dia )
            END IF
         END DO

         IF (it_durch.LT.it_end) THEN !at least one additional iteration will take place

            IF (imode_trancnf.EQ.2 .OR. imode_trancnf.EQ.3) THEN
               !new version of initializing the profile-factors using Ustar,
               !but still epressing this factor in terms of "(q*Sx)_p/(q*Sx)_0":
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  fakt=h_top_2d(i)/z0m_2d(i) !(l_p-l_0)/l_0; l_0=akt*z0m
                  rat_m_2d(i)=z1+fakt*(z1-(tkr(i)      / tkvm(i,ke1)     )) !(q*Sm)_p/(q*Sm)_0
                  rat_h_2d(i)=z1+fakt*(z1-(tkr(i)*sh_0)/(tkvh(i,ke1)*sm_0)) !(q*Sh)_p/(q*Sh)_0
               END DO

            ELSEIF (imode_trancnf.GE.4) THEN
               !new version of initializing the profile-factors and already expressing
               !them in terms of "Ustar/(q*Sh)_0*(Sh(0)/Sm(0))":
               DO i=istartpar,iendpar
                  rat_m_2d(i)= tkr(i)/tkvm(i,ke1)              !Ustar/(q*Sm)_0
                  rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(sh(0)/sm(0))
               END DO
            END IF

            IF (.NOT.ltkeinp) THEN
               nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
            END IF

         END IF

      END DO !Iteration

! 4j) Berechnung der Standardabweichnung des Saettigungsdefizites:

!k->ke1
!DIR$ IVDEP
      DO i=istartpar, iendpar
         rcld(i,ke1)=SQRT(l_tur_z0(i)*tkvh(i,ke1)*d_h)* &
!modif: previously vari' contained NOT gradient values
                       ABS(epr_2d(i)*qsat_dT(i)*vari(i,ke1,tet_l)-vari(i,ke1,h2o_g))
!modif
      ENDDO

! 4h) Berechnung der Enthalpie- und Impulsflussdichten sowie der EDR am Unterrand:

      IF ((lsrflux .AND. PRESENT(shfl_s)) .OR. lscm) THEN
!DIR$ IVDEP
         DO i=istartpar,iendpar
            shflu_s(i)=cp_d*rhon(i,ke1)*tkvh(i,ke1)*vari(i,ke1,tet_l)*epr_2d(i)
            !Note: shflu_s is positive downward and belogns to the T-equation!
         END DO
!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
         IF (lsclm) THEN
            IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
               !measured SHF has to be used for forcing:
               shflu_s(im)=SHF%mod(0)%val
            ELSEIF (lsurflu) THEN !SHF defined by explicit surface flux density
               SHF%mod(0)%val=shflu_s(im)
               SHF%mod(0)%vst=MAX(i_upd, SHF%mod(0)%vst) !SHF is at least updated
            END IF
         END IF
#endif
!SCLM --------------------------------------------------------------------------------
      END IF

      IF ((lsrflux .AND. PRESENT(qvfl_s)) .OR. lscm) THEN
!DIR$ IVDEP
         DO i=istartpar,iendpar
            qvflu_s(i)=rhon(i,ke1)*tkvh(i,ke1)*vari(i,ke1,h2o_g)
            !Note: qvflu_s is positive downward!
         END DO
      END IF
!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN
         IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN
            !measured LHF has to be used for forcing:
            qvflu_s(im)=LHF%mod(0)%val / lh_v
         ELSEIF (lsurflu) THEN !LHF defined by explicit surface flux density
            LHF%mod(0)%val=qvflu_s(im) * lh_v
            LHF%mod(0)%vst=MAX(i_upd, LHF%mod(0)%vst) !LHF is at least updated
         END IF
         !Note: LHF always is the latent heat flux connected with evaporation by definition,
         !      independent whether the surface is frozen or not!
      END IF
#endif
!SCLM --------------------------------------------------------------------------------

      !Note:
      !Both shflu_s and qvfl_s contain the flux densities of the related conserved variables now,
      ! which are equal to the usual definition only, if the cloud-water flux at the surface vanishes.

      IF (lsrflux .AND. PRESENT(umfl_s)) THEN
!DIR$ IVDEP
         DO i=istartpar,iendpar
            umfl_s(i)=rhon(i,ke1)*tkvm(i,ke1)*vari(i,ke1,u_m)
         END DO
      END IF
      IF (lsrflux .AND. PRESENT(vmfl_s)) THEN
!DIR$ IVDEP
         DO i=istartpar,iendpar
            vmfl_s(i)=rhon(i,ke1)*tkvm(i,ke1)*vari(i,ke1,v_m)
         END DO
      END IF

      IF (PRESENT(edr)) THEN
!DIR$ IVDEP
         DO i=istartpar,iendpar
            edr(i,ke1)=tke(i,ke1,ntur)**3/(d_m*l_tur_z0(i))
         END DO
      END IF

      !All these surface variables are calculated for the longer horizontal list (istartpar->iendpar),
      ! so horizontal boundary points may be included!

!----------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      IF (lnsfdia) THEN !diagnostics at near surface level required at this place

!DIR$ IVDEP
      DO i=istartpar, iendpar

         IF (itype_diag_t2m.EQ.2) THEN !using an exponetial rougness layer profile
           z2m_2d (i) = h_2m -h_can_2d(i) !2m ueber dem Bodenniveau des Bestandes
           z10m_2d(i) = h_10m-z0m_2d(i)   !Hoehe, in der turbulente Distanz 10m betraegt
         ELSE !using only a logarithmic profile above a SYNOP lawn
           z2m_2d (i) = h_2m
           z10m_2d(i) = h_10m
         END IF

         !Erste Belegung zweier benachbarter Modellniveaus:

         hk_2d(i)=h_atm_2d(i)
         hk1_2d(i)=z0
         k_2d(i)=ke

      ENDDO

!     Diagnose der 2m-Groessen:

      CALL diag_level(istartpar,iendpar, z2m_2d, k_2d, hk_2d, hk1_2d)

      IF (itype_diag_t2m.EQ.2) THEN !using an exponential rougness layer profile

         val2=z1/epsi

!DIR$ IVDEP
         DO i=istartpar, iendpar
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               IF (z2m_2d(i).LT.z0) THEN
!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

!modif: Vereinfachung der Sicherheitsabfrage:
                  val1=z2m_2d(i)/dz_s0_h(i)
                  IF (-val1.LE.val2) THEN
                    fakt=dz_s0_h(i)/dz_sa_h(i)
                    fakt=MIN( z1, MAX( z0, fakt*EXP(val1) ) )
                  ELSE
                    fakt=z0
                  ENDIF
!modif
               ELSE
!                 2m-Niveau liegt innerhalb der Modell_Prandtl-Schicht
!                 mit logarithmischen Vertikalprofilen:

                  IF (ABS(z1-fac_h_2d(i)) < epsi ) THEN
                     wert=z0m_2d(i)*z2m_2d(i)/(z2m_2d(i)+z0m_2d(i))
                  ELSEIF (fac_h_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
                     !non-stable strat. or using only linear interpolation of profile-function
                     !for the velocity scale (q*Sh):
                     wert=z0m_2d(i)*LOG((z2m_2d(i)+z0m_2d(i))/            &
                         (z0m_2d(i)+fac_h_2d(i)*z2m_2d(i)))/(z1-fac_h_2d(i))
                  ELSE !hyperbolic interpolation of (q*Sh) for stable stratification
                     wert=z0m_2d(i)*(LOG(z2m_2d(i)/z0m_2d(i)+z1)-fac_h_2d(i)*z2m_2d(i)/z0m_2d(i)) &
                                   /(z1-fac_h_2d(i))
                  END IF
                  fakt=(dz_s0_h(i)+wert)/dz_sa_h(i)
               END IF

               IF (imode_syndiag.EQ.1) THEN !direkte interpol. von temperatur und spezifischer Feuchte
                  tmps(i,ke1) = t_g(i) + (t(i,ke)-t_g(i))*fakt &
                              + tet_g*( (h_atm_2d(i)+h_can_2d(i) )*fakt-h_2m ) !Achtung: mit 'tet_g'-Korrektur
                  vaps(i,ke1)= qv_s(i) + (qv(i,ke)-qv_s(i))*fakt
               ELSE !interpolation von Erhaltungsvariablen
                  tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*ts_2d(i)/epr_2d(i)
                  vaps(i,ke1) = qds_2d(i) + fakt*(qda_2d(i)-qds_2d(i))
               END IF
               prss(i,ke1) = ps(i)

            END IF
         END DO

      ELSE !using only a logarithmic profile above a SYNOP lawn

!DIR$ IVDEP
         DO i=istartpar, iendpar
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z0d=z0d_2d(i)
               a_atm=h_atm_2d(i)+z0d
               a_2m=h_2m+z0d

!              Dimensionsloser Widerstand des Rauhigkeitsbestandes der SYNOP-Wiese,
!              wobei die turbulente Geschwindigkeitsskala und der Oberflaechenindex
!              gegenueber dem mittleren Rauhigkeitsbestand gleich bleiben:

               fr_sd_h=MAX( z0, dz_s0_h(i)/z0m_2d(i)+LOG(z0d/z0m_2d(i)) )

               IF (imode_trancnf.LT.4) THEN !profile-factor using an upper node
                  fac_h_2d(i)=(rat_h_2d(i)-z1)*z0d/h_top_2d(i)

!Achtung: Einfacher: fac_h_2d(i)=fac_h_2d(i)*z0d/z0m_2d(i)
!'fac_h_2d' sollte auch in diesem Falle wohl auch eine Profil-Konstante sein
! und nicht mit "z0d/z0m_2d" skaliert werden (also hier nicht mehr veraendert werden)
               END IF

               IF (fac_h_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sh):
                  IF (fac_h_2d(i).EQ.z1) THEN
                     val1=fr_sd_h+h_2m/a_2m
                     val2=fr_sd_h+h_atm_2d(i)/a_atm
                  ELSE
                     fakt=z1/(z1-fac_h_2d(i))
                     val1=fr_sd_h+LOG(a_2m /(z0d+fac_h_2d(i)*h_2m       ))*fakt
                     val2=fr_sd_h+LOG(a_atm/(z0d+fac_h_2d(i)*h_atm_2d(i)))*fakt
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sh) for stable stratification
                  fakt=z1/(z1-fac_h_2d(i))
                  val1=fr_sd_h+(LOG(a_2m /z0d)-fac_h_2d(i)*h_2m       /z0d)*fakt
                  val2=fr_sd_h+(LOG(a_atm/z0d)-fac_h_2d(i)*h_atm_2d(i)/z0d)*fakt
               END IF

               fakt=val1/val2

               IF (imode_syndiag.EQ.1) THEN
                  tmps(i,ke1) = t_g(i) + (t(i,ke)-t_g(i))*fakt &
                              + tet_g*(h_atm_2d(i)*fakt-h_2m)
                  vaps(i,ke1)= qv_s(i) + (qv(i,ke)-qv_s(i))*fakt
               ELSE
                  tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*ts_2d(i)/epr_2d(i)
                  vaps(i,ke1) = qds_2d(i) + fakt*(qda_2d(i)-qds_2d(i))
               END IF
               prss(i,ke1) = ps(i)

            END IF
         END DO

      END IF

!DIR$ IVDEP
      DO i=istartpar, iendpar
         IF (k_2d(i).LT.ke) THEN
!           2m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!           trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommenen
!           linearen Profile der progn. Modellvariablen eine logarith. Interpolation:

            k2=k_2d(i); k1=k2+1

            fakt=z1/(hk1_2d(i)+z0d_2d(i))
            wert=(h_2m    +z0d_2d(i))*fakt
            fakt=(hk_2d(i)+z0d_2d(i))*fakt
            fakt=LOG(wert)/LOG(fakt)
!test: different interpol. weight
!  fakt=(z2m_2d(i)-hk1_2d(i))/(hk_2d(i)-hk1_2d(i))
!test

            IF (imode_syndiag.EQ.1) THEN
               tmps(i,ke1)= t(i,k1)+fakt*( t(i,k2)- t(i,k1))
               vaps(i,ke1)=qv(i,k1)+fakt*(qv(i,k2)-qv(i,k1))
            ELSE
               val2=qv(i,k2)+qc(i,k2)
               val1=qv(i,k1)+qc(i,k1)

               vaps(i,ke1)=val1+fakt*(val2-val1)

               IF (.NOT.PRESENT(epr)) THEN
                  IF (k2.LT.ks) exner(i,k2)=zexner(prs(i,k2))
                  IF (k1.LT.ks) exner(i,k1)=zexner(prs(i,k1))
               END IF
               IF (icldm_tran.EQ.-1) THEN !no water phase change possible
                  val2=t(i,k2)/exner(i,k2)
                  val1=t(i,k1)/exner(i,k1)
               ELSE !water phase changes are possible
                  val2=(t(i,k2)-lhocp*qc(i,k2))/exner(i,k2)
                  val1=(t(i,k1)-lhocp*qc(i,k1))/exner(i,k1)
               END IF
               tmps(i,ke1)=val1+fakt*(val2-val1)

               rcl_2d(i)=rcld(i,k1)+fakt*(rcld(i,k2)-rcld(i,k1))
            END IF
            prss(i,ke1)=prs(i,k1)
         ELSEIF (z2m_2d(i).LE.z0) THEN
            rcl_2d(i)=rcld(i,ke1)
         ELSE
            fakt=z2m_2d(i)/h_atm_2d(i)
            rcl_2d(i)=rcld(i,ke1)+fakt*(rcld(i,ke)-rcld(i,ke1))
         END IF
         !Note: In the case "k_2d(i).EQ.ke" 'tmps' and 'prss' have already been
         !      caclculated above.
      END DO

      !Druck im 2m-Niveau:
!DIR$ IVDEP
      DO i=istartpar, iendpar
         wert=tmps(i,ke1)*(z1+rvd_m_o*vaps(i,ke1)) !angenaeherte virt. Temp.
         prss(i,ke1)=prss(i,ke1)                &  !Druck
                    *EXP(-(z2m_2d(i)-hk1_2d(i)) &
                    *grav/(r_d*wert))
      END DO

      IF (imode_syndiag.EQ.1) THEN
!DIR$ IVDEP
         DO i=istartpar, iendpar
             t_2m(i)=tmps(i,ke1)
            qv_2m(i)=vaps(i,ke1)
         END DO

      ELSE

!        Berechnung der zugehoerigen Modell- und Feuchtevariablen im 2m-Niveau
!        aus den Erhalturngsvariablen.

         CALL adjust_satur_equil ( khi=ke1, ktp=ke, &
!
              i_st=istartpar, i_en=iendpar, k_st=ke1, k_en=ke1,                &
!
              lcalrho=.FALSE., lcalepr=.TRUE., lcaltdv=.FALSE.,                &
              lpotinp=.TRUE. , ladjout=.TRUE.,                                 &
!
              icldmod=icldm_tran,                                              &
!
              zrcpv=zrcpv, zrcpl=zrcpl,                                        &
!
!Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
!         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
              prs=prss, t=tmps, qv=vaps, qc=liqs,                              &
!
              psf=ps,                                                          &
!
              exner=eprs, rcld=rclc,                                           &
!
              qst_t=a(:,:,3), g_tet=a(:,:,4), g_h2o=a(:,:,5),                  &
!
              tet_l=vari(:,ke:ke1,tet_l), q_h2o=vari(:,ke:ke1,h2o_g),          &
                                          q_liq=vari(:,ke:ke1,liq) )
!DIR$ IVDEP
         DO i=istartpar, iendpar
             t_2m(i)=vari(i,ke1,tet_l)
            qv_2m(i)=vari(i,ke1,h2o_g)
         END DO

      END IF

      IF (lfreeslip) THEN ! only for idealized dry runs with free-slip condition
!DIR$ IVDEP
         DO i=istartpar, iendpar
            qv_2m(i)=z0
            rh_2m(i)=z0
            td_2m(i)=z0
            u_10m(i)=z0
            v_10m(i)=z0
         END DO

      ELSE

!        Finale 2m-Diagnose:

!DIR$ IVDEP
         DO i=istartpar, iendpar
!Achtung: Macht minimale Unterschiede
            patm=prss(i,ke1)*qv_2m(i) &
                /(rdv+(z1-rdv)*qv_2m(i))          !Wasserdampfdruck
!patm=prss(i,ke1) &
!    /(rdv/MAX(qv_2m(i),epsi)+(z1-rdv))           !Wasserdampfdruck (alt)

            fakt=patm/zpsat_w( t_2m(i) )
            rh_2m(i)=100.0_ireals*MIN( fakt, z1 ) !relative Feuchte

            wert=LOG(patm/b1)
            td_2m(i)=MIN( (b2w*b3-b4w*wert) &
                      /(b2w-wert), t_2m(i) )      !Taupunktstemperatur
         END DO

!        Diagnose der 10m-Groessen:

         CALL diag_level(istartpar,iendpar, z10m_2d, k_2d, hk_2d, hk1_2d)

!DIR$ IVDEP
         DO i=istartpar, iendpar

            IF (k_2d(i).EQ.ke) THEN

!              10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z0d=z0d_2d(i)
               a_atm=h_atm_2d(i)+z0d
               a_10m=h_10m+z0d

               IF (imode_trancnf.LT.4) THEN !profile-factor using an upper node
                  fac_m_2d(i)=(rat_m_2d(i)-z1)*z0d/h_top_2d(i)
!Achtung: Einfacher: fac_m_2d(i)=fac_m_2d(i)*z0d/z0m_2d(i)
!'fac_m_2d' sollte auch in diesem Falle wohl auch eine Profil-Konstante sein
! und nicht mit "z0d/z0m_2d" skaliert werden (also hier nicht mehr veraendert werden)
               END IF

               IF (fac_m_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sm):
                  IF (fac_m_2d(i).EQ.z1) THEN
                     val1=h_10m/a_10m
                     val2=h_atm_2d(i)/a_atm
                  ELSE
                     val1=LOG(a_10m/(z0d+fac_m_2d(i)*h_10m))
                     val2=LOG(a_atm/(z0d+fac_m_2d(i)*h_atm_2d(i)))
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sm) for stable stratification
                  val1=(LOG(a_10m/z0d)-fac_m_2d(i)*h_10m      /z0d)
                  val2=(LOG(a_atm/z0d)-fac_m_2d(i)*h_atm_2d(i)/z0d)
               END IF

               fakt=val1/val2

               u_10m(i)=vel1_2d(i)*fakt
               v_10m(i)=vel2_2d(i)*fakt

            ELSE
!              10m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!              trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommen
!              lenearen Profile der progn. Modellvariablen eine logarithm. Interpolation:

               k2=k_2d(i); k1=k2+1

               fakt=z1/(hk1_2d(i)+z0d_2d(i))
               wert=(h_10m   +z0d_2d(i))*fakt
               fakt=(hk_2d(i)+z0d_2d(i))*fakt
               fakt=LOG(wert)/LOG(fakt)
!test: different interpol. weight
!  fakt=(z10m_2d(i)-hk1_2d(i))/(hk_2d(i)-hk1_2d(i))
!test
               val2=u(i,k2)
               val1=u(i,k1)
               u_10m(i)=val1+fakt*(val2-val1)

               val2=v(i,k2)
               val1=v(i,k1)
               v_10m(i)=val1+fakt*(val2-val1)
            END IF

         END DO

      END IF

      !Note: '(u, v)_10m' always belong to mass points!

      END IF !in case of ".NOT.lnsfdia" this kind of diagnostics is done at another place

!DIR$ IVDEP
      DO i=istartpar, iendpar

!        Diagnose von gz0 (fuer den naechsten Zeitschritt)
!        ueber Wasserflaechen mit der (angepassten) Charnock-Formel

         IF (fr_land(i) < z1d2) THEN

           ! Use ice surface roughness or open-water surface roughness
           ! according to lo_ice
            IF ( lo_ice(i) ) THEN
               ! Ice-covered grid box
               gz0(i)=grav*z0_ice
            ELSE !water covered surface
               velo=(tke(i,ke1,ntur)+tke(i,ke,nvor))*z1d2
!Achtung: Die 'epsi'-Beschraenkung ist recht willkuerlich und fehlt in COSMO-Version!
!Achtung: Modifikation tcm -> tvm: macht Unterschiede
               wert=MAX( epsi, tvm(i)*SQRT(vel_2d(i)**2+velo**2) ) !effective Ustar**2
               IF (imode_charpar.EQ.1) THEN !constant Charnock-Parameter
                  fakt=alpha0
               ELSE
                  IF (lini .AND. .NOT.lnsfdia) THEN
                     velo=vel_2d(i)
                  ELSE
                     velo=SQRT(u_10m(i)**2+v_10m(i)**2)
                  END IF
                  IF (depth_lk(i) > z0) THEN
                    ! enhanced Charnock parameter over lakes, parameterizing a non-equlibrium wave spectrum
                    fakt = 0.1_ireals
                  ELSE
                    fakt=alpha0_char(velo)
                  ENDIF
               END IF
               wert=MAX( grav*len_min, fakt*wert+grav*alpha1*con_m/SQRT(wert) )
               IF (ditsmot.GT.z0) THEN
                  gz0(i)=ditsmot*gz0(i)+(z1-ditsmot)*wert
               ELSE
                  gz0(i)=wert
               END IF
            END IF
         END IF
      END DO

END SUBROUTINE turbtran

!********************************************************************************

!+ Module procedure diag_level in "src_turbdiff" for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, zdia_2d, k_2d, hk_2d, hk1_2d)

   IMPLICIT NONE

   INTEGER (KIND=iintegers), INTENT(IN) :: &
!
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=ireals), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER (KIND=iintegers), INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                    !to be used of near surface diagnostics

   REAL (KIND=ireals), INTENT(INOUT) :: &
!
      hk_2d(:), & !mid level height above ground belonging to 'k_2d'
     hk1_2d(:)    !mid level height above ground of the previous layer (below)

   INTEGER (KIND=iintegers) :: i

   LOGICAL :: lcheck

   lcheck=.TRUE. !check whether a diagnostic level is above the current layer

   DO WHILE (lcheck) !loop while previous layer had to be checked
      lcheck=.FALSE. !check next layer ony, if diagnostic level is at least once
                     !above the current layer
      DO i=i_st,i_en
         IF (hk_2d(i)<zdia_2d(i) .AND. k_2d(i)>1) THEN !diagnostic level is above current layer
            lcheck=lcheck .OR. .TRUE. !for this point or any previous one the
                                      !diagnostic level is above the current layer
            k_2d(i)=k_2d(i)-1
            hk1_2d(i)=hk_2d(i)
            hk_2d(i)=(hhl(i,k_2d(i))+hhl(i,k_2d(i)+1))*z1d2-hhl(i,ke1)
          END IF
       END DO

   END DO

END SUBROUTINE diag_level

!********************************************************************************

!+ Module procedure turbdiff in "src_turbdiff" for computing the tendencies
!+ for vertical diffusion


SUBROUTINE turbdiff

!     INCLUDE 'turbdiff.incf'

!     TURBDIFF.incf

! Description:
!
!     Es werden die Diffusionskoeffizienten berechnet und ggf. Anteile
!     der zeitlichen Tendenzen der turbulenten Diffusion bestimmt
!     und zu den Tendenzfeldern hinzuaddiert.
!     Optional wird eine explizite oder (teil-)implizite Berechnung der
!     Diffusionstendenzen oder aber nur eine Berechnung der Diffusions-
!     koeffizienten durchgefuehrt. Im letzten Fall wird dann ein
!     implizit zu berechnender Anteil der Diffusionstendenzen an
!     anderer Stelle (slow_tendencies) bestimmt.
!     Allerdings koennen dann zusaetzliche explizite Korrekturtendenzen
!     hier in tubdiff bestimmt werden.

! Method:
!
!     Die Berechnung basiert auf einer Schliessung 2-ter Ordnung auf
!     dem level 2.5 (nach Mellor/Yamada). Demnach wird also eine
!     prognostische Gleichung fuer die TKE geloest.
!     Ausser der TKE-Advektion, die zusammen mit den Advektionstendenzen
!     der anderen prognostischen Variablen an anderer Stelle berechnet
!     wird, geschieht die gesamte TKE-Prognose in diesem Unterprogramm.

!     Die Formulierung des Schemas erfolgt mit thermodynamischen
!     Variablen, die bei feuchtadiabatischen vertikalen Verrueckungen
!     erhalten bleiben (pot. Fluessigw.temp. und Gesamtwassergehalt),
!     so dass der Kondesationseffekt auf subskalige Vertikalbewegungen
!     beruecksichtigt wird.
!     Die turbulenten Flussdichten der Erhaltungsgroessen werden in
!     solche der Modellvariablen konvertiert, so dass die thermodyn.
!     Kopplung der Flussdichten richtig erhalten bleibt.

!     Angeschlossen ist auch ein optionales statistisches Wolkenschema
!     (nach Sommeria und Deardorff), sub turb_cloud, welches auch
!     subskalige Bewoelkung mit Hilfe der ueber das Feld rcld ausge-
!     gebenen Standardabweichung des Saettigungsdefizites berechnet.

!     Das Turbulenzschema wurde so verallgemeinert, dass es auch bei
!     einer vertikal aufgeloesten Bestandesschicht gueltig ist, indem
!     idealisierend von der Durchstroemung eines poroesen Mediums
!     ausgegangen wird. Die Bilanzgleichungen 1-ter und 2-ter Ordnung
!     enthalten dann zusaetzliche Terme, welche die Wechselwirkungen
!     mit der Bestandes-Matrix beschreiben. Dies wirkt sich zum einen
!     auf die Stabilitaetsfunktionen und zum anderen vor allem auf die
!     TKE-Gleichung aus, welche einen auf den Formwiderstand der
!     Bestandeselemente zurueckzufuehrenden zusaetzlichen Quellterm
!     (Nachlaufturbulenz) enthaelt. Ausserdem werden die turbulenten
!     Flussdichtedivergenzen noch um einen Zusatzterm, welcher der
!     Reduktion des lufterfuellten Volumens im Gitterelement Rechnung
!     traegt, erweitert. Der Effekt des Formwiderstandes in der
!     Impulsgleichung ist ebenfalls beruecksichtigt. Die zusaetzlichen
!     Tendenzterme, die auf die Flussdichten zwischen Bestandes-Matrix
!     und umgebender Luft zurueckzufuehren sind (Bestandesquellen),
!     muessen noch in einem separaten Bestandesmodell parametrisiert
!     werden und sind nicht Gegenstand des Turbulenzmodells.
!     Schliesslich wird auch der Effekt der Transformation von Turbulenz
!     auf der dominierenden Skala in kleinskalige dissipative Turbulenz
!     durch Wirbelbrechen an Koerpern mit
!          Laengenskalen der Abmessungen << turbulente Laengenskala
!     beruecksichtigt; was sich durch eine (von der Laengenskala und
!     Volumendichte jener sehr kleinen Bestandeselemente aubhaengige)
!     Modifikation der Modellkonstanten ausdruecken laesst.

!     Es wird auch versucht den Effekt thermisch induzierter
!     Zirkulationen auf die TKE-Produktion zu beruecksichtigen
!     Hierdurch wird (vor allem) der Austauch in der naechtlichen
!     Grenzschicht erhoeht, was der Tendenz des alten Schemas,
!     in Bodennaehe zu kalte und nicht schnell genug anwachsende
!     Inversionen zu produzieren, entgegenwirkt.

!     Optional kann die Berechnung der vertikalen Gradienten durch eine
!     nicht-lokale Variante erfolgen. Hierbei werden die Gradienten
!     mit Profilen gebildet, die mit einem ueber die stabilitaets-
!     abhaengige Laengenskala gebildeten gleitenden Mittel behandelt
!     wurden.

!     Die Bildung der Anteile der Diffusionstendenzen, die numerisch
!     durch die Multiplikation mit einer Tridiagonalmatrix ausdrueckbar
!     sind, kann (neben der expliziten Variante) auch implizit erfolgen
!     (im Falle der Berechnung von nich-lokalen Gradienten und fuer die
!     TKE allerdings nur explizit).
!     Bei expliziter Rechnung ist, um auch bei Zeitschritten von
!     mehreren Minuten numerisch stabil zu bleiben, eine Limitierung
!     der Groesse der Diffusionskoeffezienten und eine im vertikalen Integral
!     quellenfreie numerische Glaettung der Vertikalprofile der Diffusions-
!     tendenzen erforderlich, sowie eine teilimplizite Behandlung der
!     Diffusionstendenzen in der untersten Modellschicht, erforderlich.

!     Die unteren Randwerte der turbulenten Flussdichten werden ueber
!     die Transferkoeffizienten zwischen Erdboden und unterster
!     Modellschicht (tcm und tch) bestimmt.
!     Optional koennen die Transferkoeffizienten auch mit diesem
!     Unterprogramm bestimmt werden, indem das Turbulenzmodell auch auf
!     das Niveau d+z0 angewandt wird, wobei vertikale Gradienten in
!     diesem Niveau mit Hilfe der Prandtl-Schicht-Hypothese berechnet
!     werden.
!     In diesem Zusammenhang wird auch die Wirkung der laminaren
!     Grenzschicht behandelt.

!     Turbulente Horizontaldiffusion (um ein 3-d-Schema zu erhalten)
!     ist noch nicht enthalten, kann aber integriert werden.
!     Uebergabevariablen:

! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  matthias.raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 Matthias Raschendorfer
!  Initial release
! 1.33       1999/10/14 Matthias Raschendorfer
!  Introduction of the LOGICAL parameters ltmpcor, lstfnct.
!  Rearranging the code in order to make it run faster.
!  Reformulation of the Charnock-Formula using the additional parameter alpha1.
!  Using the redefined parameter c_m (including the factor 3).
! 1.34       1999/12/10 Matthias Raschendofer
!  Introduction of minimal Diffusion Coefficients.
!  Modification in the Formulation of subgrid scale condensation.
!  Consideration of partial cloudiness on the TKE production due to thermal circulations.
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names.
! 2.2        2000/08/18 Matthias Raschendorfer
!  tkv(h,m) limited by the molecular diff.coef. con_(h,m) for small values.
! 2.3        2000/11/15 Guenther Doms
!  Some local variable have been redefined to allow for reproducible
!  results on MPP platforms.
! 2.15       2002/03/20 Matthias Raschendorfer
!  Some formal modifications to make the code more efficient on vector machines
! 2.17       2002/05/08 Ulrich Schaettler
!  Some more optimizations for vector machines
! 2.19       2002/10/24 Ulrich Schaettler
!  Adaptations to 2-timelevel scheme (use of ntlev)
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptation of the interface for exchg_boundaries
! 3.6        2003/12/11 Ulrich Schaettler
!  Only editorial changes
! 3.7        2004/02/18 Ulrich Schaettler
!  Increased dimension of kzdims
! 3.13       2004/12/03 Thorsten Reinhardt
!  Replaced SIGN-Function by IF-statements
! 3.16       2005/07/22 Matthias Raschendorfer
!  Introduction of the smoothing parameter 'tsmot'; removal of the parameter 'q_max'.
!  Changing the restriction in order to avoid a singularity in the stability function.
!  Modification of vertical diffusion in the TKE-equation
! 3.18       2006/03/03 Ulrich Schaettler
!  Eliminated call to exchg_boundaries for wind-tendencies, which is
!  not necessary here
! 3.19       2006/04/25 Matthias Raschendorfer / Jochen Foerstner
!  Correction of a bug related to the Exner-factor in the explicit correction
!  and some formal modifications
! 3.21       2006/12/04 Matthias Raschendorfer, Dmitrii Mironov, Ulrich Schaettler
!  Numerical security check for the second term of explicit vertical TKE diffusion
!  Introduction of the total dry option "icldm_turb.EQ.-1".
!  Changes to use the FLake model.
!  Changed interface to cloud_diag from meteo_utilities.
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introduction of some output variables for SCLM.
!  Change of some parameter names.
!  Eliminated stab_funct.incf and turb_param.incf
!  (Substitution of file 'stab_funct.incf' by a SUBROUTINE)
!  New file 'statement_functs.incf' containing statement functions
! V3_24        2007/04/26 Ulrich Schaettler
!  Introduced call to exchg_boundaries for wind-tendencies again;
!  it is necessary for imode_turb=2/3!
! V4_3         2008/02/25 Matthias Raschendorfer
!  Calculation of the 3D diagnostic field 'edr' for the eddy dissipotion rate.
!  Changing the treatment of parameter field dd(:,:) to achiev better vectorisation
!  in SUBROUTINE 'stab_funct'.
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!   Eliminated small loops to get the compiler vectorize over the correct loop
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Introduction of implicit vertical diffusion also for TKE and
!   correction of a bug related with the explicite TKE diffusion (by Oliver Fuhrer).
!  Removing the horizontal loops for SCLM treatment.
!  Adoptions with respect to the implicit vertical TKE diffusion.
!  Introduction of 3D and horizontal corrections for windshear production incl. metric terms.
!  Introduction of a separate horizontal shere mode
!   and wake turbulence terms due to the SSO-scheme.
!  Introduction of a stablity correction for turbulent length scale.
!  Modifications for seaice model: eliminated l_ls_ice, introduced lseaice
!   (Jan-Peter Schulz)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
!              2010/11/10 Matthias Raschendorfer
!  Calculation of surface heat fluxes for "itype_tran=3" in case of a SCLM run.
!
! Code Description:
! Language: Fortran 90
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================

! Declarations:

      IMPLICIT NONE

!     Lokale logical Variablen:

      LOGICAL ldotkedif, & !berechne (teil-)implizite Vert.diff von TKE
              lcircterm, & !Zirkulationsterm muss berechnet werden
              lcircdiff, & !Zirkulationsterm wird zusammen mit TKE-Diffusion bestimmt
              linisetup, & !initiales setup bei impliziter Vertikaldiffusion
              lnewvtype, & !neuer Variablentyp muss vorbereitet werden
              lsflucond    !untere Flussrandbedingung

!     Lokale Integer-Hilfsvariablen:

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &
!
              i,k,          & !Diskretis.index fuer lamda, phi, sigma
              ii,n,m,kk,    & !Indices fuer diverse Schleifen
              ku,k1,k2,     & !Schicht-Indices
              kgc,          & !oberer Schichtindex des Bereiches mit Gradientkorrektur
              it_durch,     & !Durchgangsindex der Iterationen
              nprim,        & !Erster  Index der diffundierenden Variablen
              nlast,        & !Letzter Index der diffundierenden Variablen
              ncorr,        & !Startindex der Variablen mit Gradientkorrektur
              igrdcon,      & !Index fuer Modus der Gradientberuecksichtigung
              itndcon,      & !Index fuer Modus der  Tendenzberuecksichtigung
              ivtype          !Index fuer Variablentyp

      INTEGER (KIND=iintegers) :: &
!
!             Eingrenzende Hoehenvieaus bei der Berechnung der
!             gemittelten Profile bei nicht-lokaler Gradient-Berechnung:
!
              lev(2)

!     Lokale real Variablen:

      REAL (KIND=ireals) :: &
!
!          Hilfsvariablen:
!
           wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,             & !  ,,         ,,     ,,     Faktoren
!
!     Platzh. fuer thermodynamische Hilfsgreossen
!
           virt,      & !z1+(Rv/Rd-1)*qv-qc
           flw_h2o_g, & !                 rc/(1+(lh_v/cp_d)*d_qsat/d_T)
           flw_tet_l    !exner*d_qsat/d_T*rc/(1+(lh_v/cp_d)*d_qsat/d_T)

      REAL (KIND=ireals) :: &
!
           fr_var       !1/dt_var

      REAL (KIND=ireals) :: &
!
!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2
!     (fh2,fm2):
!
           fh2,fm2, &
!
!     Platzh. fuer horiz. Geschw.-Komponenten und bel. Geschw.:
!
           vel1,vel2,velo, &
!
!     Platzh. fuer die Hoehe ueber Grund, Hoehendifferenzen, obere und
!     untere Hoehe, turbulente Laengaenskalen, Kohaerenzlaenge,
!     Dicke der laminaren Grenzschicht,sowie eine Laenge allgemein:
!
           h,hu,l_turb,lh,lm,kohae_len,len, &
!
           edh, & ! Kehrwert von Schichtdicken
!
!     Zwischenspeicher fuer
!
           thermik, & !(negative) Auftriebsproduktion an TKE
           phasdif, & !Temperaturtendenz durch Phasendiffusion
!
!
!<For_Tuning
!Achtung:
! x1,x2,x3, &
           x4, x4i, xri(ie,ke)
!>For_Tuning

! Local arrays:

!     Lokale Hilfsfelder:

      INTEGER (KIND=iintegers) :: &
!
           ivtp(ndiff) ! index of variable type

      REAL (KIND=ireals), TARGET :: &
!
           vari(ie,ke1,ndim) ,& !reduced set of variables in the first part of 'turbdiff'
                                !and later their effective vertical gradients
!
           dicke(ie,ke1),& !(effektive) Dicke einer Modellschicht
                              !bzw. sonstiges Hilfsfeld, bei k=ke1 steht
                              !eine effekt. Dicke der Modell-Prandtlsch.
!
           len_scale(ie,ke1),& !turbulent length-scale
           hor_scale(ie,ke) ,& !effective hoprizontal length-scale used for sep. horiz. shear calc.
!
           rhon(ie,ke1),& !Luftdichte auf Nebenflaechen (einschl. surface)

           shv(ie,ke1), & !velocity scale of the separated horiz. shear mode
!
           frh(ie,ke1),&  !thermal forcing (1/s2) or thermal acceleration (m/s2)
           frm(ie,ke1),&  !mechan. forcing (1/s2) or mechan. accelaration (m/s2)
           ftm(ie,ke1),&  !mechan. forcing (1/s2) by pure turbulent shear
!
           a(ie,ke1,ndim), &
!
!          a() enthaelt zu beginn die Werte der 5 thermodynamischen
!          Koeffizienten dQs/dT,ex_fakt,cp_fakt,g_tet,g_vap
!          auf den Positionen 1 bis 5 des letzten Index.
!
!          Im Falle der impliziten Berechnung von Flussdivergenzen
!          enthaelt a() hierfuer benoetigte Hilfsgroessen
!
!     3-d Hilfsfelder
!
           hlp(ie,ke1), can(ie,kcm:ke1), &
!
!     Hilfsfelder fuer eine und zwei Variablenschichten:
!
           lay(ie), lays(ie,2), vars(ie,2), &

           dzsm(ie), dzsh(ie), & !effektive Dicke der Prandtl-Schicht
!
           l_pat(ie),&    !effektive Laengenskala der therm. Inhomogenitaeten
                          !der Erdbodenoberflaeche
!
           src(ie),  &    ! arbitrary source term
!
!     Time increment and inverse time increment of ordinary prognostic variables:
!
           tinc(ndiff), tinv(ndiff), &
!
           hig(2) ! obere und untere Referenzhoehe bei der Bildung
                  ! nicht-lokaler Gradienten

      REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
           , CONTIGUOUS &
#endif
           :: &
!
!          Pointer fuer Tendenzfelder:
!
           utens(:,:), vtens(:,:), &
           ttens(:,:), qvtens(:,:), qctens(:,:), &
!
!          Pointer fuer Felder der Bestandes-Parameter:
!
           cbig(:,:), csml(:,:), rair(:,:)

      REAL (KIND=ireals), TARGET :: &
!
!          Seicher-Target fuer obige Pointer:
!
           cbig_tar(ie,kcm:ke1), csml_tar(ie,kcm:ke1), rair_tar(ie,kcm-1:ke1)

!          Note:
!          The following buffers wouldn't be necessary, if the related pointers above
!          were allowed to be allocated at run time:

      REAL (KIND=ireals), DIMENSION(:,:), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
           , CONTIGUOUS &
#endif
           :: &
!
            cur_prof, upd_prof, sav_prof, &
!++++
            expl_mom, impl_mom, invs_mom, &
            diff_dep, eff_flux
!++++

      LOGICAL ::  &
!
           ltend(ndiff), &  !calculation of tendencies required
           lsfli(ndiff), &  !surface value input is a flux density instead of a concentration
           leff_flux,    &  !calculation of effective flux density required
           lcalc_frcsmot    !local control variable if smoothing of TKE forcing terms needs to be calculated

      TYPE (modvar) :: dvar(ndiff)  !model variables to be diffused

      TYPE (turvar) :: vtyp(ntyp)   !variable types (momentum and scalars)

      TYPE (varprf) :: pvar(naux+1) !vertical variable profiles

!---- End of header ------------------------------------------------------------

!     nur in LM-Umgebung:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbdiff'; lerror=.FALSE.

!     Fuer die Turb.par. benutzter Variablensatz auf Hauptflaechen:
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)

!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     und hat einen der Werte u_m,v_m,tet_l,h2o_g,liq;
!                        bzw. tem,vap,liq.
!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)
!     Am Ende des Programms wird vari mit den (Co-Varianzen) der
!     Geschwindigkeitskomponenten ueberschrieben, die fuer die
!     Berechung der Horizontaldiffusion benoetigt werden.
!     vari() enthaelt spaeter auch die nmvar (nichtlokalen) vertikalen
!     Gradienten und auch die durch Wirkung der subskaligen Kondensation
!     veraenderten (effektiven) vertikalen Gradienten.
!     Zum Schluss enthaelt vari() fuer die turbulente Horizontaldiff.
!     benoetigte Komponenten des turbulenten Spannungstensors.

!print *,"in turbdiff c_diff=",c_diff," tkhmin=",tkhmin

!########################################################################

      fr_var=z1/dt_var

      IF (PRESENT(c_big)) THEN
         cbig => c_big
      ELSE
         cbig => cbig_tar
       ! ALLOCATE ( cbig(ie,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
      END IF
      IF (PRESENT(c_sml)) THEN
         csml => c_sml
      ELSE
         csml => csml_tar
       ! ALLOCATE ( csml(ie,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
      END IF
      IF (PRESENT(r_air)) THEN
         rair => r_air
      ELSE
         rair => rair_tar
       ! ALLOCATE ( rair(ie,kcm-1:ke1), STAT=ilocstat ); istat = istat + ilocstat
      END IF

      IF (istat /= 0) THEN
         ierrstat = 1004
         errormsg= &
         'ERROR *** Allocation of space for meteofields failed ***'
         lerror=.TRUE.; RETURN
      ENDIF

!print *,"vor init_canopy kcm=",kcm
      CALL init_canopy(ie=ie, ke=ke, ke1=ke1, kcm=kcm, &
           i_stp=istartpar, i_enp=iendpar, l_hori=l_hori, &
           fr_land=fr_land, &
           l_pat=l_pat, d_pat=d_pat, &
           c_big=cbig, c_sml=csml, r_air=rair)
!print *,"nach init_canopy kcm=",kcm

      ltend(u_m)=PRESENT(u_tens)
      IF (ltend(u_m)) THEN !calculation of tendencies required
         utens => u_tens !'utens' points to the tendency
      ELSE                 !update of ordinary prognostic variables required
         utens => u      !'utens' points to the prognostic variables
      END IF
      ltend(v_m)=PRESENT(v_tens)
      IF (ltend(v_m)) THEN
         vtens => v_tens
      ELSE
         vtens => v
      END IF
      ltend(tem)=PRESENT(t_tens)
      IF (ltend(tem)) THEN
         ttens => t_tens
      ELSE
         ttens => t
      END IF
      ltend(vap)=PRESENT(qv_tens)
      IF (ltend(vap)) THEN
         qvtens => qv_tens
      ELSE
         qvtens => qv
      END IF
      ltend(liq)=PRESENT(qc_tens)
      IF (ltend(liq)) THEN
         qctens => qc_tens
      ELSE
         qctens => qc
      END IF

      ! check if vertical smoothing of TKE forcing terms is needed
      IF (frcsmot > z0) THEN
        IF (.NOT. PRESENT(trop_mask)) THEN
          lcalc_frcsmot = .TRUE.
        ELSE IF (ANY(trop_mask(i_st:i_en) > z0)) THEN
          lcalc_frcsmot = .TRUE.
        ELSE
          lcalc_frcsmot = .FALSE.
        ENDIF
      ELSE
        lcalc_frcsmot = .FALSE.
      ENDIF

      lsfli(:)=.FALSE. !surface values are concentrations by default

      dvar(u_m)%av  => u  ; dvar(u_m)%at => utens  ; dvar(u_m)%sv => NULL()
      dvar(v_m)%av  => v  ; dvar(v_m)%at => vtens  ; dvar(v_m)%sv => NULL()

!Note: Use                                           var(u_m)%sv => u(:,ke)
!      and                                           var(v_m)%sv => v(:,ke)
!      in order to force a "free-slip condition"!

      dvar(tem)%av  => t  ; dvar(tem)%at => ttens  ; dvar(tem)%sv => t_g
      dvar(vap)%av  => qv ; dvar(vap)%at => qvtens ; dvar(vap)%sv => qv_s
      dvar(liq)%av  => qc ; dvar(liq)%at => qctens ; dvar(liq)%sv => NULL()

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN
         IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
            !measured SHF has to be used for forcing:
            lsfli(tem)=.TRUE.
         END IF
         IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN
            !measured LHF has to be used for forcing:
            lsfli(vap)=.TRUE.
         END IF
      END IF
      !Note: the measured SHF and LHF have to be present by shfl_s and qvfl_s!
#endif
!SCLM --------------------------------------------------------------------------------

      IF (lsfluse) THEN !use explicit heat flux densities at the surface
         lsfli(tem)=.TRUE.; lsfli(vap)=.TRUE.
      END IF

      IF ((lsfli(tem) .AND. .NOT.PRESENT(shfl_s)) .OR. &
          (lsfli(vap) .AND. .NOT.PRESENT(qvfl_s))) THEN
         ierrstat = 1004; lerror=.TRUE.
         errormsg='ERROR *** forcing with not present surface heat flux densities  ***'
         RETURN
      ENDIF

      IF (lsfli(tem)) dvar(tem)%sv => shfl_s
      IF (lsfli(vap)) dvar(vap)%sv => qvfl_s

      IF (PRESENT(ptr) .AND. ntrac .GE. 1) THEN !passive tracers are present
         DO m=1, ntrac
            n=liq+m
            dvar(n)%av => ptr(m)%av
            ltend(n)=ASSOCIATED(ptr(m)%at)
            IF (ltend(n)) THEN
               dvar(n)%at => ptr(m)%at
            ELSE
               dvar(n)%at => ptr(m)%av
            END IF
            IF (ASSOCIATED(ptr(m)%sv)) THEN
               dvar(n)%sv => ptr(m)%sv; lsfli(n)=ptr(m)%fc
            ELSE
               dvar(n)%sv => NULL()   ; lsfli(n)=.FALSE.
            END IF
         END DO
      END IF

      vtyp(mom)%tkv => tkvm ; vtyp(mom)%tsv => tvm
      vtyp(sca)%tkv => tkvh ; vtyp(sca)%tsv => tvh

      !Note:
      !If a tendency field of an ordinary prognostic variable is not present,
      !the related time step increment due to turbulent diffusion will be
      !added to the prognostic variable directly.

      DO n=1,ndiff
         IF (ltend(n)) THEN !calculation of tendencies required
            tinc(n)=z1        !no time increment multiplication for tendencies
            tinv(n)=fr_var    !division by time increment for variable increments
         ELSE               !update of prognostic variables required
            tinc(n)=dt_var    !time increment multiplication for tendencies
            tinv(n)=z1        !no division by time increment for variable increments
         END IF
         IF (n.LE.nvel) THEN
            ivtp(n)=mom
         ELSE
            ivtp(n)=sca
         END IF
!print *,"n=",n," tinc=",tinc(n)," tinv=",tinv(n)
      END DO

      IF (l3dturb .AND..NOT. (PRESENT(tkhm) .AND. PRESENT(tkhh))) THEN
         ierrstat = 1004; lerror=.TRUE.
         errormsg='ERROR *** 3D-diffusion with not present horiz. diff.coeffs. ***'
      END IF

!########################################################################

!-----------------------
   IF (lturatm) THEN
!-----------------------

! 0)  Berechnung der Erhaltungsvariablen (auf 'vari') samt des Bedeckungsgrades
!     und thermodynamischer Hilfgroessen, sowie der turbulenten Laengenskalen:

!-----------------------
!Achtung: Bei T-Gleichung in cv-Form
!         gesamte Thermodynamik ueberpruefen auch gegenueber satad

      lcircterm=(pat_len.GT.z0)
      ldotkedif=(c_diff .GT.z0)
      lcircdiff=(lcircterm .AND. imode_calcirc.EQ.2)

!     Berechnung abgeleiteter Parameter:

      CALL turb_param

!     Thermodynamische Hilfsvariablen auf Hauptflaechen:

      CALL adjust_satur_equil ( khi=1, ktp=1, &
!
           i_st=istartpar, i_en=iendpar, k_st=1, k_en=ke,           &
!
           lcalrho=.NOT.PRESENT(rho), lcalepr=.NOT.PRESENT(epr),    &
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        &
!
           icldmod=icldm_turb,                                      &
!
           zrcpv=zrcpv, zrcpl=zrcpl,                                &
!
           prs=prs, t=t,     qv=qv,    qc=qc,                       &
!
           psf=ps,                                                  &
!
           exner=exner, rcld=rcld, dens=rhoh,                       &
!
           r_cpd=a(:,:,2), qst_t=a(:,:,3), g_tet=a(:,:,4),          &
                                           g_h2o=a(:,:,5),          &
!
           tet_l=vari(:,:,tet_l), q_h2o=vari(:,:,h2o_g),            &
                                  q_liq=vari(:,:,liq) )

!     Thermodynamische Hilfsvariablen auf Unterrand der Prandtl-Schicht:

!DIR$ IVDEP
      DO i=istartpar, iendpar
         prss(i,ke1)=ps(i)
         tmps(i,ke1)=t_g(i)
         vaps(i,ke1)=qv_s(i)
      END DO

      CALL adjust_satur_equil ( khi=ke1, ktp=ke, &
!
           i_st=istartpar, i_en=iendpar, k_st=ke1, k_en=ke1,        &
!
           lcalrho=.TRUE., lcalepr=.TRUE.,                          &
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        &
!
           icldmod=icldm_turb,                                      &
!
           zrcpv=zrcpv, zrcpl=zrcpl,                                &
!
!Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
!         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
           prs=prss, t=tmps, qv=vaps, qc=liqs,                      &
!
           fip=tfh,                                                 &
!
           exner=a(:,ke1:ke1,1), rcld=rcld(:,ke1:ke1),              &
           dens=rhon(:,ke1:ke1),                                    &
!
           r_cpd=a(:,ke1:ke1,2), qst_t=a(:,ke1:ke1,3),              &
           g_tet=a(:,ke1:ke1,4), g_h2o=a(:,ke1:ke1,5),              &
!
           tet_l=vari(:,ke:ke1,tet_l), q_h2o=vari(:,ke:ke1,h2o_g),  &
                                       q_liq=vari(:,ke:ke1,liq) )

!     Beachte:
!     'vari(:,ke1,tet_l)' und 'vari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen
!      am Unterrand der Prandtl-Schicht (zero-level). Die Werte im Niveau 'ke' wurden dabei
!      zur Interpolation benutzt.
!     'a(:,ke1,1)' enthaelt den Exner-Faktor im zero-level. Das Feld 'a(:,:,1) wird im Folgenden
!      mit dem Exner-Faktor auf NF belegt.

!     Kommentar:
!     Der 2-te Aufruf von 'adjust_satur_equil' stellt die unteren Randwerte der thermodyn. Variablen
!      zur Verfuegung. Dies koennte in den 1-ten Aufruf integriert werden, wenn alle thermodyn.
!      Modell-Variablen bis "k=ke1" allociert waeren. Dies wuere Vieles vereinfachen!

      IF (imode_trancnf.EQ.1) THEN !old version of zero-level-values requested
         !Transformation of Tet_l at zero-level into the value following from the old
         !treatment of interpolation in terms of T_l (rather than Tet_l):
         DO i=istartpar, iendpar
            vari(i,ke1,tet_l) = vari(i,ke1,tet_l)  &
                            + ( (exner(i,ke)-a(i,ke1,1))*vari(i,ke,tet_l)*(z1-tfh(i)) ) &
                              / a(i,ke1,1)
         END DO
      END IF

!     Berechnung der horizontalen Windgeschwindigkeiten
!     im Massenzentrum der Gitterbox:

      DO k=1,ke
!DIR$ IVDEP
         DO i=istartpar,iendpar
            vari(i,k,u_m)=u(i,k)
            vari(i,k,v_m)=v(i,k)
         END DO
      END DO
!DIR$ IVDEP
      DO i=istartpar,iendpar
         vari(i,ke1,u_m)=vari(i,ke,u_m)*(z1-tfm(i))
         vari(i,ke1,v_m)=vari(i,ke,v_m)*(z1-tfm(i))
      END DO

!     Berechnung der Schichtdicken und der Dichte auf Nebenflaechen:

      DO k=1,ke
!DIR$ IVDEP
         DO i=istartpar,iendpar
            dicke(i,k)=hhl(i,k)-hhl(i,k+1)
         END DO
      END DO

!     Interpolation der thermodyn. Hilfsgroessen im Feld a(),
!     der Wolkendichte und der Luftdichte auf Nebenflaechen:

      CALL bound_level_interp( istartpar,iendpar, 2,ke, &
!___________________________________________________________________________
!test: mass weighted interpolation
!                              nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dicke)
!Achtung: Macht minimale Unterschiede
                               nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dp0)
!___________________________________________________________________________

      pvar(1)%bl => rcld     ; pvar(1)%ml => rcld  !NF-Werte wieder nach 'rcld'
      pvar(2)%bl => a(:,:,1) ; pvar(2)%ml => exner !NF-Werte nach 'a(:,:,1)', weil
                                                   !'exner' auf HF noch benoetigt wird.
      DO n=2,naux
         pvar(1+n)%bl => a(:,:,n) ; pvar(1+n)%ml => a(:,:,n)
      END DO
      !Note:
      !Internal order of level looping in 'bound_level_interp' allows to store the
      !'bl'-values (output) at the same place as the 'ml'-values (input).

      CALL bound_level_interp( istartpar,iendpar, 2,ke, &
                               nvars=naux+1, pvar=pvar, &
                               depth=dp0, rpdep=hlp)

!print *,"nach interpol"

      !Spezifische effektive Dicke der Prandtlschicht:

!DIR$ IVDEP
      DO i=istartpar,iendpar

!Achtung: < zum Vergleich mit alter Variante:
! velo=MAX( vel_min, SQRT(vari(i,ke,u_m)**2+vari(i,ke,v_m)**2) )
! tvm(i)=tcm(i)*velo
! tvh(i)=tch(i)*velo
!> zum Vergleich: in ICON werden zwischen turbtran und turbdiff 'u' und 'v' incrementiert!
!Achtung: Modifikation tcm -> tvm; tch -> tvh: macht Unterschiede
         dzsm(i)=tkvm(i,ke1)*tfm(i)/tvm(i)
         dzsh(i)=tkvh(i,ke1)*tfh(i)/tvh(i)
      END DO

      !Beachte: Auch wenn 'turbtran' nicht als Transfer-Schema benutzt wird, muessen 'tkvm', tkvh'
      !         und 'tke' und (falls es PRESENT ist) auch 'edr' fuer "k=ke1" belegt sein!

!     Berechnung der turbulenten Laengenscalen:

!DIR$ IVDEP
      DO i=istartpar,iendpar
         len_scale(i,ke1)=gz0(i)/grav
      END DO
      DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
!DIR$ IVDEP
         DO i=istartpar,iendpar
            IF (cbig(i,k).GT.z0) THEN
!              Die turbulente Laengenskala wird durch die Laengenskala
!              der lufterfuellten Zwischenraeume limitiert:

               l_turb=z1/(cbig(i,k)*sqrt(z1/exp(rair(i,k))-z1))
               len_scale(i,k)=MIN( dicke(i,k)+len_scale(i,k+1), l_turb )
            ELSE
               len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
            END IF
         END DO
      END DO
      DO k=kcm-1,1,-1
!DIR$ IVDEP
         DO i=istartpar,iendpar
            len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
         END DO
      END DO

!     Uebergang von der maximalen turbulenten Laengenskala zur
!     effektiven turbulenten Laengenskala:

      DO k=ke1,1,-1
!DIR$ IVDEP
         DO i=istartpar,iendpar
            lay(i)=l_scal(i)
         END DO
!DIR$ IVDEP
         DO i=istartpar,iendpar
            len_scale(i,k)=akt*MAX( len_min, &
                                  ! len_scale(i,k)/(z1+len_scale(i,k)/lay(i)) )
                                    lay(i)*len_scale(i,k)/(lay(i)+len_scale(i,k)) )
         END DO
      END DO

!     Initialisierung der Felder fuer tke,tkvh,tkvm:

      IF (lini) THEN  !nur beim allerersten Durchgang

!        Erste Schaetzwerte aus vereinfachtem TKE-Gleichgewicht:

         DO k=2,kem
            DO i=istartpar,iendpar

!              der Einfachheit halber nur lokale Berechnung der
!              vertikalen Gradienten:

               len=len_scale(i,k)
               edh=z2/(hhl(i,k+1)-hhl(i,k-1))

               grad(u_m  )=(vari(i,k,u_m  )-vari(i,k-1,u_m  ))*edh
               grad(v_m  )=(vari(i,k,v_m  )-vari(i,k-1,v_m  ))*edh
               grad(tet_l)=(vari(i,k,tet_l)-vari(i,k-1,tet_l))*edh
               grad(h2o_g)=(vari(i,k,h2o_g)-vari(i,k-1,h2o_g))*edh

               fh2=a(i,k,4)*grad(tet_l)+a(i,k,5)*grad(h2o_g)
               fm2=MAX( grad(u_m)**2+grad(v_m)**2, fc_min(i) )

               ! Vereinfachte Loesung mit Rf=Ri:
               IF (fh2.GE.(z1-rim)*fm2) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=len*(sm_0-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2/(fm2-fh2)
                  lm=len*(sm_0-(a_6+a_3)*fakt)
                  lh=len*(sh_0-a_5*fakt)
               END IF

               val1=lm*fm2
               val2=lh*fh2
               hlp(i,k)=MAX( val1-val2, rim*val1 )

               IF (ltkeinp) THEN
!Achtung: Bug: Auch in ICON falsch !
                  tke(i,k,1)=tke(i,k,ntur)
                ! tke(i,ke,1)=tke(i,ke,ntur)
               ELSE
                  tke(i,k,1)=MAX( vel_min, SQRT(d_m*len*hlp(i,k)) ) !Initialwert fuer SQRT(2TKE)
               END IF

               val1=MAX ( con_m, tkmmin ); tkvm(i,k)=lm*tke(i,k,1)
               val2=MAX ( con_h, tkhmin ); tkvh(i,k)=lh*tke(i,k,1)
               IF (imode_tkemini.EQ.2) THEN
                  tke(i,k,1)=tke(i,k,1)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                 val2/tkvh(i,k) )
               END IF
!Achtung: Bislang fehlte die Beschraenkung: Macht Unterschiede geg. ICON
               tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv
               tkvh(i,k)=MAX( val2, tkvh(i,k) )

!              Am Anfang konnte noch keine Diffusion von q=SQRT(2*TKE) berechnet werden:

               tketens(i,k)=z0

            END DO

         END DO

!DIR$ IVDEP
         DO i=istartpar,iendpar
            tke(i,1,1)=tke(i,2,1)
         END DO

         DO n=2,ntim
            DO k=1,kem
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  tke(i,k,n)=tke(i,k,1)
               END DO
            END DO
         END DO

      END IF

! 1)  Berechnung der benoetigten vertikalen Gradienten und
!     Abspeichern auf 'vari':

!     Am unteren Modellrand:

!DIR$ IVDEP
      DO i=istartpar,iendpar
         lays(i,mom)=z1/dzsm(i)
         lays(i,sca)=z1/dzsh(i)
      END DO
      DO n=1,nmvar
!DIR$ IVDEP
         DO i=istartpar,iendpar
            vari(i,ke1,n)=(vari(i,ke,n)-vari(i,ke1,n))*lays(i,ivtp(n))
         END DO
      END DO

!     An den darueberliegenden Nebenflaechen:

      IF (lnonloc) THEN

!        Berechnung nicht-lokaler Gradienten:

         DO n=1,nmvar

!           Berechnung vertikalen Integralfunktionen in hlp():

            DO i=istartpar,iendpar
               hlp(i,ke1)=z0
            END DO

            DO k=ke,2,-1
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  hlp(i,k)=hlp(i,k+1)+vari(i,k,n)*(hhl(i,k)-hhl(i,k+1))
               END DO
            END DO

            k1=1
            k2=2
!DIR$ IVDEP
            DO i=istartpar,iendpar
               lays(i,k1)=hhl(i,1)-hhl(i,2)
            END DO

            DO k=2,ke

!              Berechnung der nicht-lokalen Gradienten als mittlere
!              Differenzenquotienten ueber die stabilitaetsabhaengige
!              Laengenskala in tkvh() bzw tkvm():

!              Speichern der stab.abh. Laengenskala unter lay():

               IF (n.LE.nvel) THEN
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     lay(i)=tkvm(i,k)/tke(i,k,nvor)
                  END DO
               ELSE
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     lay(i)=tkvh(i,k)/tke(i,k,nvor)
                  END DO
               END IF

!              Bestimmung der nicht-lokalen Gradienten und
!              Zwischenspeichern derselben auf dem Feld dicke():

!DIR$ IVDEP
               DO i=istartpar,iendpar

                  lays(i,k2)=hhl(i,k)-hhl(i,k+1)

                  IF (lay(i).LE. &
                     z1d2*MIN( lays(i,k1), lays(i,k2)) ) THEN

!                    Die vertikalen Diffusionswege schneiden weder
!                    eine untere noch eine obere Hauptflaeche. Es
!                    koennen somit die lokalen Gradienten genommen
!                    werden. Bei sehr kleinen Diffusionslaengen, muss
!                    aus num. Gruenden sogar die lokale Berechnung
!                    gewaehlt werden:

                     dicke(i,k)=z2*(vari(i,k-1,n)-vari(i,k,n)) &
                                    /(hhl(i,k-1)-hhl(i,k+1))
                  ELSE

!                    Berechn. der benoetigten Referenzhoehen und -level:

                     h=hhl(i,k)
                     hu=MAX( h-lay(i), hhl(i,ke1) )
                     hig(1)=hu+lay(i)
                     hig(2)=h+lay(i)

                     kk=k
111                  IF (hhl(i,kk).GT.hu) THEN
                        kk=kk+1
                        GOTO 111
                     END IF
                     ku=kk
                     DO ii=1,2
112                     IF (kk.GT.1) THEN
                           IF (hhl(i,kk).LE.hig(ii)) THEN
                              kk=kk-1
                              GOTO 112
                           END IF
                        END IF
                        lev(ii)=kk+1
                     END DO

!                    Berechnung der gemittelten Differenzenquotienten
!                    als Ausdruck fuer die nicht-lokalen Gradienten:

                     wert=hlp(i,ku)-hlp(i,k) &
                         +hlp(i,lev(2))-hlp(i,lev(1)) &
                         +vari(i,ku-1,n)*(hu-hhl(i,ku)) &
                         -vari(i,lev(1)-1,n)*(hig(1)-hhl(i,lev(1)))&
                         +vari(i,lev(2)-1,n)*(hig(2)-hhl(i,lev(2)))

                     dicke(i,k)=wert/(lay(i)*(h-hu))
                  END IF
               END DO

               kk=k1
               k1=k2
               k2=kk

            END DO

!           Sichern der nicht-lokalen Gradienten im Feld vari():

            DO k=2,ke
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  vari(i,k,n)=dicke(i,k)
               END DO
            END DO

         END DO

!        Belegung von dicke() mit den Schichtdicken*rho/dt_tke
!        bzgl. Nebenflaechen:

         DO k=2,ke
!DIR$ IVDEP
            DO i=istartpar,iendpar
               dicke(i,k)=rhon(i,k)*z1d2*(hhl(i,k-1)-hhl(i,k+1))*fr_tke
            END DO
         END DO

      ELSE

!        Berechnung lokaler Gradienten:

         DO k=ke,2,-1
!DIR$ IVDEP
            DO i=istartpar,iendpar
               len=(hhl(i,k-1)-hhl(i,k+1))*z1d2
               hlp(i,k)=z1/len
               dicke(i,k)=rhon(i,k)*len*fr_tke
            END DO
         END DO

         DO n=1,nmvar
            DO k=ke,2,-1
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  vari(i,k,n)=(vari(i,k-1,n)-vari(i,k,n))*hlp(i,k)
               END DO
            END DO
         END DO

      END IF
!print *,"nach gradient"

!
! 2)  Berechnung der verallgemeinerten Antriebsfunktionen einschliesslich der
!     Korrekturen innerhalb der Rauhigkeitsschicht (samt der Windtendenz durch Formreibung)
!     und der Scherung durch nicht-turbulente subskalige Stroemungen:

!-------------------------
!     Thermal forcing:
!-------------------------

!Achtung:
      DO k=2,ke1 !'frh'(ke1) wird fuer Zirkulationsterm und Temperaturkorrektur benoetigt
!DIR$ IVDEP
         DO i=istartpar,iendpar
            frh(i,k)=a(i,k,4)*vari(i,k,tet_l) &
                    +a(i,k,5)*vari(i,k,h2o_g)
         END DO
      END DO
      !Note: a(:,:,4) and a(:,:,5) are free now.

!-------------------------
!     Total mechanical forcing:
!-------------------------

     !hdef2 = (d1v2+d2v1)^2 + (d1v1-d2v2)^2 !horizontal deformation square     (at half levels)
     !hdiv  = (d1v1+d2v2)                   !horizontal wind-divergence            ,,

     !dwdx !zonal      derivation of vertical wind                                 ,,
     !dwdy !meridional derivation of vertical wind                                 ,,

     !vel_div=hdiv+dzdw=0 !Incomressibility

     !itype_sher = 0 : only single column vertical shear
     !             1 : previous and additional 3D horiz. shear correction
     !             2 : previous and additional 3D vertc. shear correction

     !ltkeshs: consider separated non-turbulent horizontal shear mode for TKE forcing
     !loutshs: consider separated non-turbulent horizontal shear mode for output

!     Mechanical forcing by vertical shear:

      IF (itype_sher.EQ.2 .AND. (PRESENT(dwdx) .AND. PRESENT(dwdy) .AND. PRESENT(hdiv))) THEN
         !Include 3D-shear correction by the vertical wind (employing incomressibility):
         DO k=2,kem
!DIR$ IVDEP
            DO i=istartpar,iendpar
               frm(i,k)=MAX( (vari(i,k,u_m)+dwdx(i,k))**2+(vari(i,k,v_m)+dwdy(i,k))**2 &
                                                         +z3*hdiv(i,k)**2, fc_min(i) )
            END DO
         END DO
      ELSE
         !Load pure single column shear:
         DO k=2,kem
!DIR$ IVDEP
            DO i=istartpar,iendpar
               frm(i,k)=MAX( vari(i,k,u_m)**2+vari(i,k,v_m)**2, fc_min(i))
            END DO
         END DO
      END IF

!     Mechanical forcing by horizontal shear:

      IF (PRESENT(hdef2)) THEN

         IF (itype_sher.GE.1) THEN
            !Apply horizontal 3D-shear correction:
            DO k=2,kem
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  frm(i,k)=frm(i,k)+hdef2(i,k) !extended shear
               END DO
            END DO
         END IF

      END IF

      IF (lssintact) THEN !save pure turbulent shear
         DO k=2,kem
!DIR$ IVDEP
            DO i=istartpar,iendpar
               ftm(i,k)=frm(i,k)
            END DO
         END DO
      END IF

!For_Tuning>
!     Preparation for Richardson-number-dependent factor used for correcting
!     the minimum diffusion coefficient and the horizontal shear production term:

      DO k=2,ke
!DIR$ IVDEP
        DO i=istartpar,iendpar
!Achtung: Mit Hilfe von 'frm' und 'frh' auszudruecken: <
! x1 = z1d2*(hhl(i,k-1)-hhl(i,k+1))
! x2 = MAX( 1.e-6_ireals, ((u(i,k-1)-u(i,k))**2+(v(i,k-1)-v(i,k))**2)/x1**2 )          ! |du/dz|**2
! x3 = MAX( 1.e-5_ireals, grav/(z1d2*(t(i,k-1)+t(i,k)))*((t(i,k-1)-t(i,k))/x1+tet_g) ) !       N**2

! xri(i,k) = EXP(z2d3*LOG(x2/x3))  ! 1/Ri**(2/3)

          xri(i,k)=EXP( z2d3*LOG( MAX( 1.e-6_ireals, frm(i,k) ) / & !1/Ri**(2/3)
                                  MAX( 1.e-5_ireals, frh(i,k) ) ) )
!>
        END DO
      END DO
!>For_Tuning

      IF (PRESENT(hdef2)) THEN

         !Additional impact by separated horizontal shear:
         IF ((ltkeshs .OR. (loutshs .AND. PRESENT(tket_hshr))) .AND. PRESENT(hdiv)) THEN
            !Include separated horizontal shear mode:

            fakt=z1/(z2*sm_0)**2; wert=a_hshr*akt*z1d2

!DIR$ IVDEP
            DO i=istartpar,iendpar
               lay(i)=wert*l_hori(i) !uncorrected effective horizontal length scale
            END DO

            IF (imode_shshear.EQ.2) THEN
!>Tuning
              DO k=2,kem
!DIR$ IVDEP
                DO i=istartpar,iendpar
                  ! Factor for variable 3D horizontal-vertical length scale proportional to 1/SQRT(Ri),
                  ! decreasing to zero in the lowest two kilometers above ground

                  x4i = MIN( 1._ireals, 0.5e-3_ireals*(hhl(i,k)-hhl(i,ke1)) )
                  x4 = 3._ireals*x4i**2 - 2._ireals*x4i**3                   ! low-level reduction factor
                  hor_scale(i,k) = lay(i)*MIN( 5.0_ireals, MAX( 0.01_ireals, x4*xri(i,k) ) )
                END DO
              END DO
!>Tuning: This kind of correction can be substituded by a less ad-hoc approach.
            ELSE
              DO k=2,kem
!DIR$ IVDEP
                DO i=istartpar,iendpar
                  hor_scale(i,k) = lay(i)
                END DO
              END DO
            ENDIF

            DO k=2,kem
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  wert=fakt*hdiv(i,k)
                  shv(i,k)=hor_scale(i,k)*(SQRT(wert**2+hdef2(i,k))-wert) !strain velocity of the sep. hor. shear mode
                  hlp(i,k)=(shv(i,k))**3/hor_scale(i,k)                   !additional TKE-source by related shear
               END DO
            END DO
            IF (loutshs .AND. PRESENT(tket_hshr)) THEN
               !Load output variable for the TKE-source by separated horiz. shear:
               DO k=2,kem
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     tket_hshr(i,k)=hlp(i,k)
                  END DO
               END DO
            END IF
            IF (ltkeshs) THEN
               !Consider separated horizontal shear mode in mechanical forcing:
               DO k=2,kem
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     frm(i,k)=frm(i,k)+hlp(i,k)/tkvm(i,k) !extended shear
                  END DO
               END DO
               IF (l3dturb) THEN
                  ! Load related horizontal diffusion coefficients:
                  fakt=sh_0/sm_0
                  DO k=2,kem
!DIR$ IVDEP
                     DO i=istartpar,iendpar
                        tkhm(i,k)=hor_scale(i,k)*shv(i,k) !for momentum related to the sep. shear mode
                        tkhh(i,k)=fakt*tkhm(i,k)          !for scalars    ,,       ,,            ,,
                     END DO
                  END DO
               END IF
            END IF

         END IF

      END IF

      !Addition verallgemeinerter Scherterme durch nicht-turbulente subskalige Stroemung:

      DO k=2,kem

         IF (.NOT.lini) THEN !nicht bei der Initialisierung

            IF (lsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
               !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:

!              Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO-Tendenzen:

!DIR$ IVDEP
               DO i=istartpar,iendpar
!Achtung: horizontale Pos. beachten!
                 vel1=-(ut_sso(i,k)  *u(i,k  )+vt_sso(i,k)  *v(i,k  ))*dp0(i,k-1)
                 vel2=-(ut_sso(i,k-1)*u(i,k-1)+vt_sso(i,k-1)*v(i,k-1))*dp0(i,k)

                 src(i)=MAX( z0, (vel1+vel2)/(dp0(i,k-1)+dp0(i,k)) )

!                Beachte:
!                Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
!                erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
!                Obwohl vel1 und vel2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!
               END DO

               IF (PRESENT(tket_sso)) THEN
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     tket_sso(i,k)=src(i)
                  END DO
               END IF

!              Addition des Scherterms durch Nachlaufwirbel:

               IF (imode_tkesso == 1) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)
                  END DO
               ELSE IF (imode_tkesso == 2) THEN ! Reduce TKE production in the presence of large Richardson numbers
!DIR$ IVDEP
                  DO i=istartpar,iendpar
                     frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)*MIN(1.0_ireals,MAX(0.01_ireals,xri(i,k)))
                  END DO
               END IF

            END IF

            IF (lconv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
               !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

!              Addition des Scherterms durch konvektive Zirkulation:

!DIR$ IVDEP
               DO i=istartpar,iendpar
                  frm(i,k)=frm(i,k) + MAX( z0, tket_conv(i,k)/tkvm(i,k) )
               END DO

!              Beachte:  Obwohl tket_conv immer positiv sein muesste, wird zur Sicherheit MAX benutzt!
            END IF
         END IF

      END DO

!     Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht
!     (ausser Volumenterme, die zur Diffusion gehoeren):

      DO k=kcm,kem !von oben nach unten durch Rauhiggkeitsschicht

!Achtung: Neue Behandlung der Rauhigkeitsschicht einfuehren

!        Vertikalwind und Formreibungskoeff. auf Hauptflaechen:

!DIR$ IVDEP
         DO i=istartpar,iendpar
            lay(i)=z1d2*(w(i,k)+w(i,k+1))
            src(i)=z1d2*(cbig(i,k)  +csml(i,k) &
                        +cbig(i,k+1)+csml(i,k+1))
         END DO

!        Windbetrag auf der angrenzenden Hauptflaeche:

         IF (k.EQ.kcm) THEN
!DIR$ IVDEP
            DO i=istartpar,iendpar
               velo=z1d2*(w(i,k-1)+w(i,k))
               lays(i,1)=SQRT(u(i,k-1)**2+v(i,k-1)**2+velo**2)
             END DO
         END IF

!        Reduzierte Konstanten und implizite Horizontalwind-Tendenzen durch Formreibung:

!DIR$ IVDEP
         DO i=istartpar,iendpar

!           Windbetrag auf der aktuellen Hauptflaeche:
            lays(i,2)=SQRT(u(i,k)**2+v(i,k)**2+lay(i)**2)

!           Aufaddieren der Windtendenzen durch Fromreibung:
            wert=src(i)*lays(i,2)

            utens(i,k)=utens(i,k)-tinc(u_m)*wert*u(i,k)/(z1+dt_var*wert)
            vtens(i,k)=vtens(i,k)-tinc(v_m)*wert*v(i,k)/(z1+dt_var*wert)

!           Windbetrag auf Nebenflaechen:
            can(i,k)=(lays(i,2)*dp0(i,k-1)+lays(i,1)*dp0(i,k))/(dp0(i,k-1)+dp0(i,k))

!           Windbetrag unter Wirkung der Formreibung:
            velo=can(i,k)/(z1+can(i,k)*(cbig(i,k)+csml(i,k))*dt_var)

!           Addition des Scherterms durch Nachlaufwirbel an grossen Rauhigkeitselementen:
            frm(i,k)=frm(i,k)+cbig(i,k)*velo**3/tkvm(i,k)

!           Frequenz der kleinskaligen Rauhigkeitselemente:
            can(i,k)=csml(i,k)*can(i,k)

         END DO

      END DO  !k=kcm,kem !von oben nach unten druch Rauhigkeitsschicht

      !Optionale vertikale Glaettung des mechanischen Antriebs:

      IF (lcalc_frcsmot) THEN
         CALL vert_smooth ( &
              i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1, &
              disc_mom=dicke, cur_tend=frm, vertsmot=frcsmot, smotfac=trop_mask )
      END IF

!     Optionale vertikale Glaettung des thermischen Antriebs:

      IF (lcalc_frcsmot) THEN
         CALL vert_smooth ( &
              i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1, &
              disc_mom=dicke, cur_tend=frh, vertsmot=frcsmot, smotfac=trop_mask )
      END IF

!     Belegung von tkvh und tkvm mit den stabilitaetsabhaengigen Laengenmassen:

      DO k=2,kem
!DIR$ IVDEP
         DO i=istartpar,iendpar
            tkvh(i,k)=tkvh(i,k)/tke(i,k,nvor)
            tkvm(i,k)=tkvm(i,k)/tke(i,k,nvor)
         END DO
      END DO

!return

! 3)  Loesung der turbulenten Bilanzgleichungen (Bestimmung von TKE und der Stabilitaetsfunktionen)
!     und Berechnung der turbulenten Diffusionskoeffizienten:

      DO it_durch=it_start, it_end

        !Die Schleife wird nur bei der Initialisierung (d.h. beim ersten Aufruf) wiederholt,
        !um TKE-Gleichgewicht anzunaehern. Die resultierenden TKE-Werte der Zeitstufe 'ntur'
        !gehoeren in diesem Fall dann zur Zeitstufe der uebrigen prognostischen Variablen
        !('nold' bei "leap-frog" oder 'nnow' bei 2-Zeitebenen).
        !Fuer die folgenden Aufrufe wird die Schleife nur einmal durchlaufen und liefert TKE-Werte
        !die gegenueber den Vorgaengerwerten um einen Zeitschritt weiter in der Zukunft liegen,
        !also wieder zur Zeitstufe der uebrigen prognostischen Variablen gehoeren.

         CALL solve_turb_budgets (khi=1, it_s=it_durch, &
                                  i_st=istartpar, i_en=iendpar, k_st=2, k_en=kem, &
!
                                  lssintact=lssintact, lupfrclim=.FALSE.,         &
!
                                  imode_stke=imode_turb, imode_vel_min=1,         &
!
                                  fm2=frm, fh2=frh, ft2=ftm, lsm=tkvm, lsh=tkvh,  &
#ifdef SCLM
                                  grd=vari,                                       &
#endif
                                  fcd=can, tls=len_scale, tvt=tketens, avt=tketadv)

         IF (it_durch.LT.it_end .AND. .NOT.ltkeinp) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF

      END DO

!     Kein TKE-Gradient am Oberrand:

!DIR$ IVDEP
      DO i=istartpar,iendpar
         tke(i,1,ntur)=tke(i,2,ntur)
      END DO

      IF (iini.EQ.1) THEN !only for separate initialization before the time loop

         DO k=2, kem
!DIR$ IVDEP
            DO i=istartpar,iendpar
!Achtung: Bislang fehtlte die laminare Beschraenkung
               val1=MAX ( con_m, tkmmin ); tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
               val2=MAX ( con_h, tkhmin ); tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
               IF (imode_tkemini.EQ.2) THEN
                  tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                       val2/tkvh(i,k) )
               END IF
               tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv
               tkvh(i,k)=MAX( val2, tkvh(i,k) )

            END DO
         END DO

       ! IF (.NOT.PRESENT(rhoh))  DEALLOCATE ( rhoh , STAT=ilocstat )
       ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )

         RETURN !finish this subroutine

      END IF

!  4) Berechnung der effektiven turbulenten Vertikalgradienten,
!     Standardabweichnung des Saettigungsdefizites und der Diffusionskoeffizienten:

      IF (ltmpcor) THEN
      !  Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:
         DO k=2, ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
               frm(i,k)=a(i,k,1)*vari(i,k,tet_l)-tet_g  !vertical temperature gradient
            END DO
         END DO
         IF (icldm_turb.NE.-1) THEN !water phase changes are possible
            DO k=2, ke1
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  frm(i,k)=frm(i,k)+lhocp*vari(i,k,liq) !liquid water correction
               END DO
            END DO
         END IF

      !  Dies geschieht schon hier, weil im naechsten Schritt das Feld vari()
      !  durch die effiktiven Gradienten ueberschrieben wird.
      END IF

      DO k=2, ke1
!DIR$ IVDEP
         DO i=istartpar,iendpar
            a(i,k,5)=a(i,k,1)*a(i,k,3)                           !exner*d_qsat/d_T
            a(i,k,4)=c_scld*rcld(i,k)/(z1+rcld(i,k)*(c_scld-z1)) !effective cloud cover
            rcld(i,k)=SQRT(len_scale(i,k)*tkvh(i,k)*d_h)* & !standard deviation
                           ABS(a(i,k,5)*vari(i,k,tet_l)   & !of local
                              -vari(i,k,h2o_g))             !oversaturation
         END DO
      END DO

      IF (icldm_turb.NE.-1) THEN !consideration of water phase changes
         DO k=2, ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
             ! Effective vertical gradient of liquid water content:

               flw_h2o_g=a(i,k,4)/(z1+lhocp*a(i,k,3))  !weight of h2o_g-flux
               flw_tet_l=-flw_h2o_g*a(i,k,5)           !weight of tet_l-flux

               vari(i,k,liq)=  flw_tet_l*vari(i,k,tet_l) &                    ! eff_grad(liq)
                            +  flw_h2o_g*vari(i,k,h2o_g)

             ! Effective vertical gradient of water vapor content and pot. temper.:

               vari(i,k,vap)= vari(i,k,h2o_g)-vari(i,k,liq)                   ! eff_grad(vap)
               vari(i,k,tet)=(vari(i,k,tet_l)+vari(i,k,liq)*lhocp/a(i,k,1)) & ! eff_grad(tet)
                             *a(i,k,2)                                        !*(Cp/Cpd)

             ! Note:
             ! -flux(h2o_g)/(rho*K)=grad(h2o_g)=eff_grad(vap)+eff_grad(liq)
             ! -flux(tet_l)/(rho*K)=grad(tet_l)=eff_grad(tet)-eff_grad(liq)*lh_v/(cp_d*exner)
             ! Treating "lh_v/(cp_d*exner)" like a constnat, besides numerical effects,
             !  vertical gradient diffusion of non conserved variables without a moist gradient correction
             !  does not change the resulting diffusion of conserved variables. Thus the redistribution of
             !  water constituents can (or should) be left to a final (sub grid scale) saturation adjustment.
            END DO
         END DO
      ELSE !no water phase change possible
        !'vari(:,:,liq)' bleibt unveraendert, genauso wie auch
        !'vari(:,:,tet)'='vari(:,:,tet_l)' und 'vari(:,:,vap)'='vari(:,:,h2o_g)'.
         DO k=2, ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
               vari(i,k,tet)=vari(i,k,tet_l)*a(i,k,2) ! grad(tet)*(Cp/Cpd)
            END DO
         END DO
      END IF
      !Beachte:
      !Die eff_grad ergeben multipliziert mit 'rho*tkvh' die vertikalen Flussdichten
      ! von 'tet', 'vap' und 'liq' unter Beruecksichtigung der turbulenten Phasenuebergaenge.
      !'vari(:,:,tet)' ist der in der Teta-Gleichung benoetigte effective Gradient
      ! und entspricht der 'tet'-Flussdichte*(Cp/Cpd).
      !Die Indices 'tet', 'tem' und 'tet_l' haben den gleichen Wert
      ! so wie auch 'vap' und 'h2o_g' den gleichen Wert haben.
      !Die Unterscheidungen sollen nur den jeweiligen Inhalt verdeutlichen.
      !Im Falle "icldm_turb.EQ.EQ.-1" verschwindet der eff. Bedeckungsgrad "a(:,:,4)=0",
      ! troltzdem steht auch dann die Standardabweichung der Uebersaettigung in 'rcld'.

!     Beschraenkung der Diffusionskoeffizienten:

      val1=tkmmin; val2=tkhmin !default minimum diffusion values

      DO k=2, ke
!DIR$ IVDEP
         DO i=istartpar,iendpar
!           Berechn. der Diffusionskoeffizienten:

            IF (imode_tkvmini.EQ.2) THEN
!>Tuning
               ! Factor for variable minimum diffusion coefficient proportional to 1/SQRT(Ri);
               ! the namelist parameters tkhmin/tkmmin specify the value for Ri=1:
               fakt=MIN( z1, tkred_sfc(i)*(0.25_ireals+7.5e-3_ireals*(hhl(i,k)-hhl(i,ke1))) ) !low-level red.-fact.
               fakt=MIN( 2.5_ireals, MAX( 0.01_ireals, fakt*xri(i,k) ) )

               val1=tkmmin*fakt; val2=tkhmin*fakt

               IF (tkhmin_strat.GT.z0 .OR. tkmmin_strat.GT.z0) THEN
                  ! Enhanced diffusion in the stratosphere - needed primarily for momentum because 
                  ! there is otherwise too little dynamic coupling between adjacent model levels
                  fakt = MIN( z1, 2.e-4_ireals*MAX( z0, hhl(i,k) - 12500._ireals ) ) ! transition zone between 12.5 and 17.5 km
                  ! Wider transition zone in the tropics in order to avoid too strong diffusion in the tropopause region
                  x4  = z1-z1d3*trop_mask(i)*MIN(z1, 2.e-4_ireals*MAX(z0, 22500._ireals-hhl(i,k)) )
                  x4i = z1-z2d3*innertrop_mask(i)*MIN(z1, 2.e-4_ireals*MAX(z0, 27500._ireals-hhl(i,k)) )
                  fakt = fakt*MIN( x4*1.5_ireals, MAX( 0.25_ireals, SQRT(xri(i,k)) ) )
                  val1=MAX( val1, tkmmin_strat*MIN(x4,x4i)*fakt ) ; val2=MAX( val2, tkhmin_strat*x4*fakt )
               END IF
!>Tuning: This kind of correction can be substituted by a less ad-hoc approach.
! Remark (GZ): The enhanced stratospheric diffusion seems to parameterize a missing process outside the turbulence scheme,
! maybe momentum transports due to non-stationary gravity waves. This may also explain why we need a much larger
! minimum diffusion coefficient for momentum than for heat.
            END IF

!Achtung: Beschraenkung mit lam. diff.coef. fehlte bislang auch in ICON; macht ev. Unterschiede
            val1=MAX ( con_m, val1 ); tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
            val2=MAX ( con_h, val2 ); tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
            IF (imode_tkemini.EQ.2) THEN
               tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                    val2/tkvh(i,k) )
            END IF
            tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv's
            tkvh(i,k)=MAX( val2, tkvh(i,k) )

!test: ohne tk?min (wird aber bei Diffusion benutzt!)
!test: ohne tk?min und ohen lam. Diffkof. (wird aber bei Diffusion benutzt!)

!           tkvh und tkvm enthalten jetzt nicht mehr Diffusions-
!           laengen, sondern wieder Diffusionskoeffizienten in m^2/s!

         END DO
      END DO

      IF (l3dturb) THEN
         !Consider horizontal diffusion coefficients:
         IF (PRESENT(hdef2) .AND. PRESENT(hdiv) .AND. ltkeshs) THEN
            !Add isotropic turbulent part to that part due to the sep. horiz. shear mode:
            DO k=2,kem
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  tkhm(i,k)=tkhm(i,k)+tkvm(i,k)
                  tkhh(i,k)=tkhh(i,k)+tkvh(i,k)
               END DO
            END DO
         ELSE !no treatment of sep. horiz. shear mode has taken place
            !Load only the isotropic turbulent part:
            DO k=2,kem
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  tkhm(i,k)=tkvm(i,k); tkhh(i,k)=tkvh(i,k)
               END DO
            END DO
         END IF
      END IF

!  5) Untere Randwerte der Turbulenzgroessen:

      !Beachte: Auch wenn 'turbtran' nicht als Transfer-Schema benutzt wird,
      !         muessen 'tkvm', tkvh' und 'tke' und (falls es PRESENT ist)
      !         auch 'edr' fuer "k=ke1" belegt sein!

      IF (lsfli(tem)) THEN !use explicit shfl_s
!DIR$ IVDEP
         DO i=istartpar,iendpar
            vari(i,ke1,tet)=shfl_s(i)/(cp_d*rhon(i,ke1)*tkvh(i,ke1)*a(i,ke1,1))
         END DO
         !Note:'vari(tet)' belongs to potential temperature,
         !     and 'a(1)' contains Exner factor on half levels.
      END IF
      IF (lsfli(vap)) THEN !use explicit qvfl_s
!DIR$ IVDEP
         DO i=istartpar,iendpar
            vari(i,ke1,vap)=qvfl_s(i)/(rhon(i,ke1)*tkvh(i,ke1))
         END DO
      END IF
      !Note: "tkvh(ke1) > 0.0" is required!

! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperaturtendenzen
!     (ausser der Divergenz des Zirkulationsflusses):

      IF (ltmpcor) THEN
         IF (.NOT.PRESENT(edr)) THEN
!DIR$ IVDEP
            DO i=istartpar,iendpar
               ediss(i,ke1)=tke(i,ke1,ntur)**3/(d_m*len_scale(i,ke1))
            END DO
         END IF

         DO k=2, ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
               thermik=tkvh(i,k)*frh(i,k)
               phasdif=tkvh(i,k)*frm(i,k) &
                      *(zrcpv*vari(i,k,vap)+zrcpl*vari(i,k,liq))

               tketens(i,k)=len_scale(i,k)/a(i,k,2) &
                       *((ediss(i,k)+thermik)/cp_d+phasdif)
            END DO

!           Beachte:
!           In 'frm' steht an dieser Stelle der Temp.-Gradient!
!           Wegen der nachfolgenden Interpolation auf Hauptflaechen,
!           wird 'tketens' mit der turbulenten Laengenskala multipliziert.
         END DO

!DIR$ IVDEP
         DO i=i_st,i_en
            ttens(i,1)=ttens(i,1)+tinc(tem)*tketens(i,2) &
                         /(len_scale(i,1)+len_scale(i,2))

         END DO
         DO k=2,ke
!DIR$ IVDEP
            DO i=i_st,i_en
               ttens(i,k)=ttens(i,k)+tinc(tem)*(tketens(i,k)+tketens(i,k+1)) &
                            /(len_scale(i,k)+len_scale(i,k+1))
            END DO
         END DO
      END IF

! 7)  Bestimmung des Zirkulationstermes als zusaetzliche TKE-Flussdichte:
!Achtung: Zirkulationsterm revidieren:

      IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden

         DO k=2,ke1

            IF (k.LT.ke1) THEN
!              Interpolation des Druckes auf die naechst hoehere
!              Nebenflaeche:

!DIR$ IVDEP
               DO i=istartpar,iendpar
                  lay(i)=(prs(i,k)*dp0(i,k-1)+prs(i,k-1)*dp0(i,k)) &
                          /(dp0(i,k)+dp0(i,k-1))
               END DO
            ELSE
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  lay(i)=ps(i)
               END DO
            END IF

!DIR$ IVDEP
            DO i=istartpar,iendpar

!Achtung: Variation durch Wolkenbedeckung
               fakt=z1-z2*ABS(a(i,k,4)-z1d2)

               len=MAX( l_pat(i), SQRT(fakt*len_scale(i,k)*l_hori(i)) )

!              Berechnung der lokalen Kohaerenzlaenge fuer die Parameterisierung
!              des Zirkulationstermes (kohae_len=l_pat*exnr*grad(tet_v)*r/g):

               fakt=frh(i,k)*lay(i)/(rhon(i,k)*grav**2)
               kohae_len=len*SIGN(z1,fakt)*MIN( ABS(fakt), z1 )

!              Belegung von 'frh' mit der TKE-Flussdichte des Zirkulationstermes
!              (Zirkulationsfluss oder Drucktransport):

               frh(i,k)=tkvh(i,k)*kohae_len*frh(i,k)

!              Die Divergenz dieser Flussdichte ist gleichzeitig
!              eine weitere Quelle fuer thermische Energie.

            END DO

!           Addition des Teta-Gradienten, welcher zur Teta-Flussdichte
!           durch den Zirkulationsterm gehoert:

            IF (ltmpcor) THEN
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  vari(i,k,tet)=vari(i,k,tet) &
                                -frh(i,k)/(tkvh(i,k)*a(i,k,1)*a(i,k,2)*cp_d)
               END DO
            END IF

            !Die kinetische Energie der Zirkulation wird der thermischen Energie
            !entzogen!

         END DO

      END IF

! 8) Berechnung der Diffusionstendenz von q=SQRT(2*TKE) einschliesslich der
!    q-Tendenz durch den Zirkulationsterm:

!     Vorbereitung zur Bestimmung der zugehoerigen Incremente von TKE=(q**2)/2:

      cur_prof => hlp
      upd_prof => a(:,:,1)
      sav_prof => a(:,:,5)

      IF (ldotkedif .OR. lcircdiff) THEN

         expl_mom => a(:,:,2)
         impl_mom => a(:,:,3)
         invs_mom => a(:,:,4)

         ! Diffusions-Koeffizienten auf NF:
         DO k=2, ke1
            IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
!DIR$ IVDEP
               DO i=istartpar,iendpar
!___________________________________________________________________________
!test: TKE-Diffusion mit Stab.fnkt. fuer Skalare:
! sav_prof(i,k)=c_diff*tkvh(i,k)
                  sav_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)
!___________________________________________________________________________
               END DO
            ELSE !Diffusion in terms of q=SQRT(2*TKE)
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  sav_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)**2
               END DO
            END IF
         END DO

         DO k=3, ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
               expl_mom(i,k)=rhoh(i,k-1)*z1d2*(sav_prof(i,k-1)+sav_prof(i,k)) &
                                             /(hhl(i,k-1)-hhl(i,k))
            END DO
            !Beachte:
            !'expl_mom' bezieht sich auf HF, also die Fluss-Niveaus fuer die TKE (bzw. q-)-Diffusion.
            !Wegen der spaeteren Nutzung der SUBs 'prep_impl_vert_diff' und 'calc_impl_vert_diff'
            ! muss ein Fluss-Niveau (hier HF) ueber dem Variabl.-Niveau (hier NF) mit gleichem Index liegen.
!Achtung: In der COSMO-Version werden 'rhon', 'len_scale' und 'tke' einzeln auf HF interpoliert,
!         was geringe Unterschiede verursacht!
         END DO

      END IF

      IF (ldotkedif .OR. lcircterm) THEN
         IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
            DO k=2, ke1
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  sav_prof(i,k)=z1d2*tke(i,k,ntur)**2 !TKE
               END DO
            END DO
         ELSE !Diffusion in terms of q=SQRT(2*TKE)
            DO k=2, ke
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  sav_prof(i,k)=tke(i,k,ntur)         !q=SQRT(2*TKE)
                  dicke(i,k)=dicke(i,k)*tke(i,k,ntur) !related effective discretization momentum
               END DO
            END DO
            !Das Feld 'dicke' wird bei "k=ke1" nicht benoetigt und war zuvor dort auch nicht belegt!
            DO i=istartpar,iendpar
               sav_prof(i,ke1)=tke(i,ke1,ntur) !q at the surface layer
            END DO
         END IF
         !Beachte:
         !Das Feld 'tke' enthaelt nicht TKE sondern q!
      END IF

!     Aufnahme des Zirkulationstermes:

      IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden

!        Interpolation der Zirulationsflussdichte auf Hauptflaechen:

         !Beachte:
         !'frm' ist frei.
         !'frh' enthaelt bislang eine TKE-Flussdichte (bis auf den Faktor 'rhon'),
         ! deren Vertikalprofil im wesentlichen durch d_z(tet_v)**2 bestimmt ist,
         ! was zumindest in der Prandtl-Schicht prop. zu 1/len_scale ist.
         !Die Interpolation von "rhon*frh" auf Hauptflaechen erfolgt daher mit
         ! 'rhon*frh*len_scale'. Anderenfalls ist mit grossen Interpolationsfehlern zu rechnen.

         DO k=2,ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
!Achtung: In COSMO-Version ist "imode_tkediff=1".
!         Es wird in dieser Version aber 'frh' an dieser Stelle durch 'tke' dividiert,
!         was geringe Unterschiede verursacht:"
               frh(i,k)=rhon(i,k)*frh(i,k)*len_scale(i,k) !skalierte Flussdichte auf NF
!frh(i,k)=rhon(i,k)*frh(i,k)/tke(i,k,ntur)*len_scale(i,k) !skalierte Flussdichte auf NF

!if (k.lt.ke1) then
!   frh(i,k)=rhon(i,k)*grav*tkvh(i,k)*t(i,k)/t_g(i)/tke(i,k,ntur)*len_scale(i,k)
!else
!   frh(i,k)=rhon(i,k)*grav*tkvh(i,k)/tke(i,k,ntur)*len_scale(i,k)
!endif
            END DO
         END DO

!        Korrektur der TKE-Profile durch die Zirkulations-Tendenz:

         IF (imode_calcirc.EQ.1) THEN !expliziten Berechnung der Zirkulationstendenz

            cur_prof => sav_prof

            DO k=3,ke1
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  frm(i,k)=(frh(i,k)+frh(i,k-1))/(len_scale(i,k)+len_scale(i,k-1)) !interpolierte Flussdichte auf HF
               END DO
            END DO
!DIR$ IVDEP

!Achtung: In COSMO-Version ist "imode_tkediff=1".
!Da 'frh' bereits durch 'tke' dividiert wurde, muss fuer die exakte COSMO-Version der 'tke'-Faktor in 'dicke'
!wieder beseitigt werden:
            k=2
            DO i=istartpar,iendpar
               upd_prof(i,k)=sav_prof(i,k)+frm(i,k+1)/dicke(i,k)
!upd_prof(i,k)=sav_prof(i,k)+frm(i,k+1)*tke(i,k,ntur)/dicke(i,k)

!upd_prof(i,k)=sav_prof(i,k)
!tketens(i,k)=frm(i,k+1)*tke(i,k,ntur)/dicke(i,k)

            END DO
            DO k=3,ke
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  upd_prof(i,k)=sav_prof(i,k)-(frm(i,k)-frm(i,k+1))/dicke(i,k)
!upd_prof(i,k)=sav_prof(i,k)-(frm(i,k)-frm(i,k+1))*tke(i,k,ntur)/dicke(i,k)

!Achtung: In der COSMO-Version wird die explizite q-Tendenz durch den Zirkulationsterm
!erst nach der TKE-Diffusion hinzuaddiert (verursacht geringe Unterschiede):
!upd_prof(i,k)=sav_prof(i,k)
!tketens(i,k)=-(frm(i,k)-frm(i,k+1))*tke(i,k,ntur)/dicke(i,k)

               END DO
            END DO
            !Beachte:
            !'dicke' enthaelt bereits den Faktor "1/dt_tke" (sowie den Faktor 'tke' bei "imode_tkediff=1").
            !'upd_prof' enthaelt das um die Zirkulations-Tendenz aufdatierte TKE-Profil
            ! (oder q-Profile bei "imode_tkediff=1")

            !Zuschlag durch Volumenterm aus der Divergenzbildung:
            DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
               DO i=istartpar,iendpar
!Achtung: Korrektur
                ! upd_prof(i,k)=upd_prof(i,k)+dt_tke*frh(i,k)*z1d2*(rair(i,k-1)-rair(i,k+1)) &
                  upd_prof(i,k)=upd_prof(i,k)+frh(i,k)*z1d2*(rair(i,k-1)-rair(i,k+1)) &
                                                           /(len_scale(i,k)*dicke(i,k))
                  !'frh' enthaelt die Zirkultions-Flussdichte der TKE auf NF (skaliert mit 'len_scale').
               END DO
            END DO

            !Bereucksichtige Zirkulations-Tendenz:
             itndcon=1 !indem 'upd_prof' auf rechter Seite der impliz. Diff.-Gl. benutzt wird.
!Achtung: Um die COSMO-Version exakt nachzubilden, darf die Zirkulationstendenz nicht bei der
!impliziten Diffusions-Gleichung benutzt werden:
!itndcon=0
            !Fuer die expliziten Diff.-Anteile wird aber 'cur_prof' benutzt.

         ELSE !quasi implizite Berechnung der Zirkulatinstendenz (entspricht "lcircdiff=T")
!DIR$ IVDEP
            DO i=istartpar,iendpar
               cur_prof(i,2)=sav_prof(i,2)
            END DO
            DO k=3,ke1
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  wert=(frh(i,k)+frh(i,k-1))/((len_scale(i,k)+len_scale(i,k-1))*expl_mom(i,k))
                  cur_prof(i,k)=(cur_prof(i,k-1)-sav_prof(i,k-1)+wert)+sav_prof(i,k)
               END DO
            END DO
            !Beachte:
            !'cur_prof' enthaelt ein virtuelles TKE-Profil (oder q-Profile bei "imode_tkediff=1"),
            ! dessen Diffusions-Tendenz die Zirkulations-Tendenz einschliesst.

            !Bereucksichtige Zirkulations-Tendenz:
            itndcon=0 !indem 'cur_prof' auf rechter Seite der impliz. Diff.-Gl. benutzt wird.
            !Fuer die expliziten Diff.-Anteile wird ebenfalls 'cur_prof' benutzt.
         END IF

      ELSEIF (ldotkedif) THEN

         cur_prof => sav_prof

         itndcon=0 !'cur_prof' wird auf rechter Seite der impliz. Diff.-Gl.
                   ! und fuer explizite Diff.-Anteile benutzt.
      END IF

!     Aufdatieren des TKE-Profils durch die (erweiterte) Diffusions-Tendenz

      IF (ldotkedif .OR. lcircdiff) THEN

         !'frm', 'frh' und 'len_scale' sind frei.
         !In den Diffusionroutinen wird vorausgesetzt, dass ein Flussniveau mit gleichem
         !Vertikalindex wie ein Konzentrationsniveau gerade ueber letzterem liegt.
         !Die bisherige Hauptflaechenindizierung musste daher fuer die Uebergabefelder
         !der Routinen 'prep_impl_vert_diff' und 'calc_impl_vert_diff' um "1" verschoben werden.

         CALL prep_impl_vert_diff( lsflucond=.FALSE., ldynimpwt=ldynimp, lprecondi=lprecnd, &
                                   i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1,  &
!Achtung: q_Diff
              disc_mom=dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
!disc_mom=sav_prof*dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
              invs_fac=frh, scal_fac=frm, impl_weight=impl_weight )

!        Berechnung der vertikalen Diffusionstendenzen von TKE=z1d2*q**2:

         eff_flux => len_scale

         CALL calc_impl_vert_diff ( lsflucond=.FALSE.,lprecondi=lprecnd, &
                                    leff_flux=(kcm.LE.ke), itndcon=-itndcon, &
                                    i_st=istartpar, i_en=iendpar ,k_tp=1, k_sf=ke1, &
!Achtung: q_Diff
              disc_mom=dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
!disc_mom=sav_prof*dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
              invs_fac=frh, scal_fac=frm, cur_prof=cur_prof, upd_prof=upd_prof, eff_flux=eff_flux )

         !Beachte:
         !Bei "imode_diff=1" erfolgt q-Diffusion, so dass 'disc_mom' und 'expl_mom' den zusaetzlichen
         ! Faktor 'sav_prof'='tke'=q enthalten!
         !'upd_prof' enthaelt jetzt die mit der Diffusionstendenz aufdatierten (modifizierte) TKE-Werte.
         !Weil "itndcon<=0", bleiben auf 'cur_prof' die Eingangsprofile erhalten.
         !'eff_flux' enthaelt die effektiven Flussdichten (positiv abwaerts) der (semi-)impliziten
         ! Vertikaldiffusion.

         IF (lcircdiff) THEN !es wurden virtuelle Effektiv-Profile benutzt
            DO k=2,ke
!DIR$ IVDEP
               DO i=istartpar,iendpar
                  upd_prof(i,k)=sav_prof(i,k)+upd_prof(i,k)-cur_prof(i,k) !aufdatierte echte Profile
               END DO
            END DO
         END IF
!++++

!        Zuschlag durch Volumenterm innerhalb der Rauhigkeitsschicht:

         DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
            DO i=istartpar,iendpar
               wert=(eff_flux(i,k)*dp0(i,k)+eff_flux(i,k+1)*dp0(i,k-1))/(dp0(i,k)+dp0(i,k-1))
                   !effektive TKE-Flussdichte interpoliert auf die k-te Nebenflaeche,
                   !wobei sich 'eff_flux(:,k)' auf die (k-1)-te Hauptflaeche bezieht!
!Achtung: Korrektur
             ! upd_prof(i,k)=upd_prof(i,k)+dt_tke*wert*z1d2*(rair(i,k-1)-rair(i,k+1))/dicke(i,k)
               upd_prof(i,k)=upd_prof(i,k)+wert*z1d2*(rair(i,k-1)-rair(i,k+1))/dicke(i,k)
            END DO
         END DO
         !'upd_prof' enthaelt das mit dem Volunmenterm beaufschlagte aufdatierte Profil.

      END IF

!     Speichern der zugehoerigen q-Tendenzen:

      IF (ldotkedif .OR. lcircterm) THEN

         IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
            !'upd_prof' ist ein TKE-Profil:
            DO k=2,ke
!DIR$ IVDEP
               DO i=istartpar,iendpar
!___________________________________________________________________________
!test:
                ! tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - sav_prof(i,k) )*fr_tke/tke(i,k,ntur)
                  tketens(i,k)=( SQRT( 2*MAX( upd_prof(i,k), z0 ) ) - tke(i,k,ntur) )*fr_tke
!___________________________________________________________________________
               END DO
            END DO
         ELSE !Diffusion in terms of q=SQRT(2*TKE)
            !'upd_prof' ist ein q-Profil:
            DO k=2,ke
!DIR$ IVDEP
               DO i=istartpar,iendpar
!Achtung:
                 tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - tke(i,k,ntur) )*fr_tke

!Achtung: Bei der COSMO-Version wird erst hier die explizite Zirkulations-Tendenz hinzuaddiert:
!tketens(i,k)=MAX( -tke(i,k,ntur), tketens(i,k)+( upd_prof(i,k) - tke(i,k,ntur) ) )*fr_tke
               END DO
            END DO
         END IF
         !'tketens' enthaelt jetzt immer eine q-Tendenz!

         !Am Unterrand gibt es keine q-Tendenz durch Diffusionsterme:
!DIR$ IVDEP
         DO i=istartpar,iendpar
            tketens(i,ke1)=z0
         END DO

!        Optionale vertikale Glaettung der erweiterten Diffusionstendenz von q=SQRT(2*TKE):

         IF (tndsmot.GT.z0) THEN
            CALL vert_smooth ( &
                 i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1, &
                 disc_mom=dicke, cur_tend=tketens, vertsmot=tndsmot )
         END IF

      ELSE !keine q-Tendenzen, weder durch TKE-Diffusion noch durch den Zirkulationsterm

!        Zuruecksetzen der q-Tendenzen:

         DO k=2,ke1
!DIR$ IVDEP
            DO i=istartpar,iendpar
               tketens(i,k)=z0
            END DO
         END DO

      END IF

! 10) Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!     des Saettigungsdefizites:

!DIR$ IVDEP
      DO i=istartpar,iendpar
         rcld(i,1)=rcld(i,2)
      END DO
      DO k=2,kem-1
!DIR$ IVDEP
         DO i=istartpar,iendpar
            rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
         END DO
      END DO
!     Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
!     der Wert auf der entspr. Nebenflaeche beibehalten.

!-----------------------
   END IF !lturatm
!-----------------------

!########################################################################

!Achtung: ".NOT.lini" ergaenzt
!--------------------------------------------------
      IF ((ldovardif .OR. ldogrdcor) .AND. .NOT.lini) THEN !Vertikaldiffusion wird hier berechnet
!--------------------------------------------------
         !Note:
         !If ".NOT.ldovardif .AND. ldogrdcor", only a correction of pure vertical gradient diffusion
         ! due to sub grid scale condensation or non-local gradients is performed.

!        Berechnung der Luftdichte und des Exner-Faktors am Unterrand:

!DIR$ IVDEP
         DO i=istart,iend
            virt=z1+rvd_m_o*qv_s(i) !virtueller Faktor
            rhon(i,ke1)=ps(i)/(r_d*virt*t_g(i))
            eprs(i,ke1)=zexner(ps(i))
         END DO
         !Note:
         !In the turbulence model 'rhon(:,ke1)' belongs to the lower boundary of the
         !Prandtl-layer, rather than to the surface level.

!        Berechnung von Hilfsgroessen:

         IF (.NOT.lturatm) THEN !turbulence model was not running

            IF (.NOT.PRESENT(rho)) THEN
               DO k=1,ke
!DIR$ IVDEP
                  DO i=istart,iend
                     virt=z1+rvd_m_o*qv(i,k)-qc(i,k) !virtueller Faktor
                     rhoh(i,k)=prs(i,k)/(r_d*virt*t(i,k))
                  END DO
               END DO
            END IF

            CALL bound_level_interp( istart,iend, 2,ke, &
!___________________________________________________________________________
!test: mass weighted interpolation
!                              nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=hhl, auxil=hlp)
                               nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dp0)
!___________________________________________________________________________

            IF (lscadif .AND. .NOT.PRESENT(epr)) THEN
               DO k=1, ke
!DIR$ IVDEP
                  DO i=istart,iend
                     exner(i,k)=zexner(prs(i,k))
                  END DO
               END DO
            END IF

         END IF

!        Setzen von Steuerparametern:

         IF (ldogrdcor) THEN
            IF (lnonloc) THEN
               ncorr=1
            ELSE
               ncorr=nvel+1
            END IF
         ELSE
            ncorr=ndiff+1
         END IF

         ivtype=0

!-----------------------------------------------------------------
!        Berechnung der Vertikaldiffuion von Modellvariablen auf Hauptflaechen:
!-----------------------------------------------------------------

!        DO n=nprim, nlast !loop over all variables to be diffused
         DO n=1, ndiff !loop over all variables to be diffused potentially
         IF ( (lum_dif .AND. n.EQ.u_m)   .OR. &                   !u_m-diffusion or
              (lvm_dif .AND. n.EQ.v_m)   .OR. &                   !v_m-diffusion or
              (lscadif .AND. n.GT.nvel)  .OR. &                   !sca-diffusion or
            (ldogrdcor .AND. n.GE.ncorr .AND. n.LE.nmvar) ) THEN !gradien correction

            m=MIN(n,nmvar)

            IF (ivtype.EQ.0) THEN
               linisetup=.TRUE.; lnewvtype=.TRUE.
            ELSE
               linisetup=.FALSE.;lnewvtype=.FALSE.
            END IF

            IF (n.LE.nvel) THEN !a wind component
               IF (ivtype.EQ.sca) lnewvtype=.TRUE.

!Achtung:
               lsflucond=.FALSE. !no lower flux condition for momentum!
!lsflucond=lsflcnd !no lower flux condition for momentum!

               ivtype=mom
            ELSE
               IF (ivtype.EQ.mom) lnewvtype=.TRUE.

               lsflucond=lsflcnd !use chosen type of lower boundary condition

               ivtype=sca
            END IF

            IF (n.LT.ncorr .OR. n.GT.nmvar) THEN
               igrdcon=0 !keine Gradientkorrektur der Profile
            ELSEIF ( (.NOT.lscadif .AND. ivtype.EQ.sca) .OR. &
                     (.NOT.lum_dif .AND.      n.EQ.u_m) .OR. &
                     (.NOT.lvm_dif .AND.      n.EQ.v_m) ) THEN
               igrdcon=1 !nur Profile aus Gradientkorrektur
            ELSE
               igrdcon=2 !korrigierte Profile aus effektiven Gradienten
            END IF

            kgc=2 !uppermost level of gradient correction

            IF (.NOT.ltend(n)) THEN !tendency array not present
               itndcon=0 !no explicit tendency consideration
            ELSE
               itndcon=itnd !use chosen mode of tendency consideration
            END IF
!test: never expl tendency consideration
!itndcon=0
!test

            IF (lsfli(n)) THEN
               !Load effective surface layer gradients due to given flux values:

!DIR$ IVDEP
               DO i=istart,iend
                  vari(i,ke1,m)=dvar(n)%sv(i)/(rhon(i,ke1)*vtyp(ivtype)%tkv(i,ke1))
               END DO
               IF (n.EQ.tem) THEN !flux density is that of sensible heat
!DIR$ IVDEP
                  DO i=istart,iend
                     vari(i,ke1,m)=vari(i,ke1,m)/(cp_d*eprs(i,ke1))
                  END DO
               END IF
               !Note:
               !In this case not the current surface concentration but the current flux-density
               ! at the surface is used in 'vert_grad_diff'!
               !Hoewever, 'vari' contains vertical gradients at this place, and for ".NOT.lsflucond"
               ! a related surface concentration is recalculated in 'vert_grad_diff'.
               !'tkv(ke1)' needs to be >0, which is always the case, if calculated by 'turbtran'!
               !Thus "tkvh(ke1)=0.0" should not be forced, if "lsfli=.TRUE"!
               !For tracers it is "m=nmvar"!
               !In case of "lturatm=T" 'vari(ke1)' has already been loaded using shlf_s or qvfl_s!
            END IF

!           Belegung der Eingangsprofile und -Tendenzen:

            cur_prof => hlp

            DO k=1,ke
!DIR$ IVDEP
               DO i=istart,iend
                  cur_prof(i,k)=dvar(n)%av(i,k)
               END DO
            END DO

            IF (ASSOCIATED(dvar(n)%sv)) THEN !surface variable is present
!DIR$ IVDEP
               DO i=istart,iend
                  cur_prof(i,ke1)=dvar(n)%sv(i)
               END DO
            ELSEIF (n.LE.nvel .OR. ilow_def_cond.EQ.2) THEN
               !No-slip-condition for momentum or zero-concentr.-condition as a default:
!DIR$ IVDEP
               DO i=istart,iend
                  cur_prof(i,ke1)=z0
               END DO
            ELSE !enforce a zero flux condition as a default
!DIR$ IVDEP
               DO i=istart,iend
                  cur_prof(i,ke1)=cur_prof(i,ke)
               END DO
            END IF
            IF (itndcon.GT.0) THEN !explicit tendencies have to be considered
               DO k=1,ke
!DIR$ IVDEP
                  DO i=istart,iend
                     dicke(i,k)=dvar(n)%at(i,k)
                  END DO
               END DO
            END IF

            IF (n.EQ.tem) THEN !temperature needs to be transformed
               DO k=1,ke
!DIR$ IVDEP
                  DO i=istart,iend
                     cur_prof(i,k)=cur_prof(i,k)/exner(i,k) !potential temperature
                  END DO
               END DO
!DIR$ IVDEP
               DO i=istart,iend
                  cur_prof(i,ke1)=cur_prof(i,ke1)/eprs(i,ke1)
               END DO
               IF (itndcon.GT.0) THEN !explicit tendencies to be considered
                  DO k=1,ke
!DIR$ IVDEP
                     DO i=istart,iend
                        dicke(i,k)=dicke(i,k)/exner(i,k)
                     END DO
                  END DO
               END IF
            END IF

            IF (.NOT.(lsfluse .AND. lsflcnd)) THEN ! calculation of effective flux density required
               IF ( ( n.EQ.tem .AND. PRESENT(shfl_s) ) .OR. ( n.EQ.vap .AND. PRESENT(qvfl_s) ) ) THEN
                  leff_flux = .TRUE.
               ELSE
                  leff_flux = .FALSE.
               END IF
            ELSE
               leff_flux = .FALSE.
            END IF

!           Berechnung der vertikalen Diffusionstendenzen:

            CALL vert_grad_diff( kcm, kgc,                            &
!
                 i_st=istart, i_en=iend, k_tp=0, k_sf=ke1,            &
!
                 dt_var=dt_var, ivtype=ivtype, igrdcon=igrdcon, itndcon=itndcon, &
!
                 linisetup=linisetup, lnewvtype=lnewvtype,            &
                 lsflucond=lsflucond, lsfgrduse=lsfli(n),             &
                 ldynimpwt=ldynimp  , lprecondi=lprecnd,              &
                 leff_flux=leff_flux,                                 &
!
                 rho=rhoh, rho_n=rhon, hhl=hhl, r_air=rair,           &
!
                 tkv=vtyp(ivtype)%tkv, tsv=vtyp(ivtype)%tsv,          &
!
                 impl_weight=impl_weight,                             &
!
                 disc_mom=a(:,:,1), expl_mom=a(:,:,2),                &
                 impl_mom=a(:,:,3), invs_mom=a(:,:,4),                &
                 diff_dep=a(:,:,5), diff_mom=len_scale,               &
                 invs_fac=frh, scal_fac=frm,                          &
!
                 dif_tend=dicke, cur_prof=cur_prof, eff_flux=vari(:,:,m) )

            !Beachte:
            !'frh', 'frm' und 'len_scale' sind genauso wie 'a(:,:,1:5)' Hilfsspeicher in 'vert_grad_diff'.
            !Weil Fluesse ab "n>=liq=nmvar" nicht mehr benoetigt werden, bleibt 'vari' nur bis
            ! 'nmvar' dimensioniert und vari(nmvar) wird auch fuer "n>nmvar" benutzt.

!           Sichern der Tendenzen:

            IF (n.EQ.tem) THEN
               DO k=1,ke
!DIR$ IVDEP
                  DO i=istart,iend
                     dvar(n)%at(i,k)=dvar(n)%at(i,k)+exner(i,k)*dicke(i,k)*tinc(n)
                  END DO
               END DO
            ELSE
               DO k=1,ke
!DIR$ IVDEP
                  DO i=istart,iend
                     dvar(n)%at(i,k)=dvar(n)%at(i,k)+dicke(i,k)*tinc(n)
                  END DO
               END DO
            END IF

            IF (n.EQ.vap .AND. PRESENT(qv_conv)) THEN
               !qv-flux-convergence (always a tendency) needs to be adapted:
               DO k=1,ke
                  IF (lqvcrst) THEN
                     !by initializing 'qv_conv' with vertical qv-diffusion:
!DIR$ IVDEP
                     DO i=istart,iend
                        qv_conv(i,k)=dicke(i,k)
                     END DO
                  ELSE !by adding vertical qv-diffusion to 'qv_conv':
!DIR$ IVDEP
                     DO i=istart,iend
                        qv_conv(i,k)=qv_conv(i,k)+dicke(i,k)
                     END DO
                  END IF
               END DO
            END IF

         END IF !diffusion calculation requested
         END DO !1, ndiff

!-----------------------------------------------------------------

!Achtung:
!Ist cp-Fluss tatsaechlich der thermische Erdbodenantrieb?
!Was gilt im Falle der T-Gleichung in cv-Form?

!        Berechnung der effektiven Oberflaechenflussdichten:

!Achtung: "lscadif" ergaenzt
         IF (.NOT.(lsfluse .AND. lsflcnd) .AND. lscadif) THEN
            !effektive Oberfl.flussdichten wurden neu bestimmt

            IF (PRESENT(shfl_s) .OR. lscm) THEN
!DIR$ IVDEP
               DO i=istart, iend
                  shflu_s(i)=eprs(i,ke1)*cp_d*vari(i,ke1,tet)
               END DO
            END IF
            IF (PRESENT(qvfl_s) .OR. lscm) THEN
!DIR$ IVDEP
               DO i=istart, iend
                  qvflu_s(i)=vari(i,ke1,vap)
               END DO
            END IF

!---------------------------------------------------------------------------------------
#ifdef SCLM
            IF (lsclm .AND. latmflu) THEN
               !Berechnung der Enthalpieflussdichten:

               SHF%mod(0)%val=shflu_s(im)     ; SHF%mod(0)%vst=i_cal
               LHF%mod(0)%val=qvflu_s(im)*lh_v; LHF%mod(0)%vst=i_cal

               !Note:
               !IF ".NOT.latmflu", SHF and LHF either are loaded by the fluxes used for
               ! the soil budget (lertflu) or they have been loaded above by the explicit
               ! SHF and LHF at the surface (lsurflu).
               !SHF and LHF are positive downward and they may have been corrected with
               ! vertical integrated correction tendencies.
               !Thus they always refer to the used flux densities, which are only then equal
               ! to the explicit surface flux density, if a lower flux condition is used "lsflcnd=.TRUE.".
            END IF
#endif
!SCLM-----------------------------------------------------------------------------------

            !Bem: shflu_s und qvflu_s, sowie SHF und LHF sind positiv abwaerts!

         END IF

         IF (lum_dif .AND. PRESENT(umfl_s)) THEN
!DIR$ IVDEP
            DO i=istart, iend
               umfl_s(i)=vari(i,ke1,u_m)
            END DO
         END IF
         IF (lvm_dif .AND. PRESENT(vmfl_s)) THEN
!DIR$ IVDEP
            DO i=istart, iend
               vmfl_s(i)=vari(i,ke1,v_m)
            END DO
         END IF

!--------------------------------------------------
      END IF !Vertikaldiffusion wird hier berechnet
!--------------------------------------------------


! 11) Berechnung der Tendenzen infloge horizontaler turb. Diffusion
!     (und aufsummieren auf die Tendenzfelder):

! 16) Deallocierung lokaler dynamischer Felder:

    ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )

END SUBROUTINE turbdiff

!********************************************************************************

SUBROUTINE adjust_satur_equil ( khi, ktp, &
!
   i_st, i_en, k_st, k_en,            &
!
   lcalrho, lcalepr,                  &
   lcaltdv,                           &
   lpotinp, ladjout,                  &
!
   icldmod,                           &
!
   zrcpv, zrcpl,                      &
!
   prs, t, qv, qc,                    &
   psf, fip,                          &
!
   exner, rcld, dens, r_cpd,          &
!
   qst_t, g_tet, g_h2o, tet_l,        &
   q_h2o, q_liq )

INTEGER (KIND=iintegers), INTENT(IN) :: &
  khi,        & !usual start index of vertical dimension (highest level)
  ktp,        & !extra start index of vertical dimension (top level)
                !including the auxilary level
!
  i_st, i_en, & !horizontal start- and end-indices
  k_st, k_en    !vertical   start- and end-indices

LOGICAL, INTENT(IN) :: &
  lcalrho, &   !density calculation required
  lcalepr, &   !exner pressure calculation required
  lcaltdv, &   !calculation of special thermodynamic variables required
  lpotinp, &   !input temperature is a potential one
  ladjout      !output of adjusted variables insted of conserved ones

INTEGER (KIND=IINTEGERS), INTENT(IN) :: &
  icldmod      !mode of water cloud representation in transfer parametr.
               !-1: ignoring cloud water completely (pure dry scheme)
               ! 0: no clouds considered (all cloud water is evaporated)
               ! 1: only grid scale condensation possible
               ! 2: also sub grid (turbulent) condensation considered

REAL (KIND=ireals), DIMENSION(:,khi:), INTENT(IN) :: &
  prs,     &   !current pressure
  t, qv        !current temperature and water vapor content (spec. humidity)

REAL (KIND=ireals), DIMENSION(:,khi:), OPTIONAL, INTENT(IN) :: &
  qc           !current cloud water content

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
  psf, &       !surface pressure
  fip          !interpolation factor with respect to an auxilary level

REAL (KIND=ireals), DIMENSION(:,khi:), INTENT(INOUT) :: &
  exner, &     !current Exner-factor
  rcld         !inp: standard deviation of oversaturation
               !out: saturation fraction

REAL (KIND=ireals), TARGET, INTENT(INOUT) &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
  tet_l(:,ktp:), &     !inp: liquid water potent. temp. (only if 'fip' is present)
               !out: liquid water potent. temp. (or adjust. 't' , if "ladjout")
  q_h2o(:,ktp:), &     !inp: total  water content (only if 'fip' is present)
               !out: total  water content       (or adjust. 'qv', if "ladjout")
  q_liq(:,ktp:)        !out: liquid water content after adjustment

!REAL (KIND=ireals), DIMENSION(:,khi:), OPTIONAL, INTENT(INOUT) :: & !doesn't work on NAG-compiler
 REAL (KIND=ireals), DIMENSION(ie,k_st:k_en), OPTIONAL, INTENT(INOUT) :: &
  dens,  &     !current air density
  r_cpd        !cp/cpd

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, INTENT(OUT) :: &
  qst_t, &     !out: d_qsat/d_T
               !aux:                                    TARGET for 'qvap', if ".NOT.ladjout"
  g_tet, &     !out: g/tet-factor of tet_l-gradient in buoyancy term
               !aux: total  water content of ref. lev.; TARGET for 'temp', if ".NOT.ladjout"
  g_h2o        !out: g    -factpr of q_h2o-gradient in buoyancy term
               !aux: liq. wat. pot. temp. of ref. lev.; TARGET for 'virt'

REAL (KIND=ireals), INTENT(IN) :: &
  zrcpv,  &    !0-rcpv  switch for cp_v/cp_d - 1
  zrcpl        !0-rcpl  switch for cp_l/cp_d - 1

REAL (KIND=ireals) :: &
  pdry,  &     !corrected pot. temp. and partial pressure of dry air
  ccov,  &     !effective cloud cover
  mcor         !moist correction

REAL (KIND=ireals), DIMENSION(i_st:i_en,k_st:k_en) :: &
  rprs          !reduced pressure

REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
  temp(:,:), &  !corrected temperature
  qvap(:,:), &  !corrected water vapour content
  virt(:,:)     !reciprocal virtual factor

INTEGER (KIND=iintegers) :: &
  i, k, &
  icldtyp  !index for type of cloud diagnostics in SUB 'turb_cloud'

!-------------------------------------------------------------------------------

   !Note:
   !If 'qc' is not present, it is assumed that 't' and 'qv' already contain liquid water temperature
   ! and total water content.
   !In case of "k_st=khi=k_en", the vertical loop is only for one level.
   !Since a surface variable for 'qc' is not used in general, it is not forseen here as well, assuming
   ! a zero value at the surface. This implicats the surface variables to be identical with
   ! liquid water temperature and total water content (conserved variables for moist water conversions).
   !If 'fip' is resent, the conseved variables 'tet_l' and q_h2o' at level 'k' are an interpolation between these
   ! values and those of the level "k-1" using the weight 'fip'. This is used  in order to interpolate
   ! rigid surface values at level "k=ke1" and those valid for the lowermost atmospheric full level "k=ke"
   ! to an atmopspheric surface level at the top of the roughness layer formed by land use. In this case it is
   ! lways "k_st=ke1=k_en".
   !In this version, atmospheric ice is not included to the adjustment process.

   !Calculation of Exner-pressure:
   IF (lcalepr) THEN
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            exner(i,k)=zexner(prs(i,k))
         END DO
      END DO
   END IF

   !Conserved variables (with respect to phase change):
   IF (icldmod.EQ.-1 .OR. .NOT.PRESENT(qc)) THEN
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            q_h2o(i,k)=qv(i,k)
            tet_l(i,k)= t(i,k)
         END DO
      END DO
   ELSE !water phase changes are possible and 'qc' is present
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            q_h2o(i,k)=qv(i,k) +       qc(i,k) !tot. wat. cont.
            tet_l(i,k)= t(i,k) - lhocp*qc(i,k) !liq. wat. temp.
         END DO
      END DO
   END IF

   !Transformation in real liquid water temperature:
   IF (lpotinp) THEN
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            tet_l(i,k)=exner(i,k)*tet_l(i,k)
         END DO
      END DO
   END IF

   !Interpolation of conserved variables (with respect to phase change)
   !onto the zero-level of interest at the top of the roughness layer:
   IF (PRESENT(fip)) THEN
      k=k_en !only for the lowest level
!DIR$ IVDEP
      DO i=i_st,i_en
         q_h2o(i,k)=           q_h2o(i,k-1)*(z1-fip(i))+q_h2o(i,k)*fip(i)
         tet_l(i,k)=exner(i,k)*tet_l(i,k-1)*(z1-fip(i))+tet_l(i,k)*fip(i)
      END DO
      !Note:
      !'tet_l' at level "k-1" needs to be present and it is assumed to be already
      ! a potential (liquid water) temperature there, where it is still a pure
      ! liquid water temperature at level 'k'!
      !The roughness layer between the current level (at the rigid surface)
      ! and the desired level (at the top of the roughness layer) of index "k=ke1"
      ! is assumend to be a pure CFL layer without any mass. Consequentley
      ! the Exner pressure values at both levels are treated like being equal!
   END IF

   !Berechnung der effektiven Bedeckung mit Wasserwolken:

!Achtung: Korrektur: An COSMO-Version angepasste Bedeutung von "icldmod=-1":
!  IF (icldmod.LE.0 ) THEN
   IF (icldmod.EQ.0 .OR. (icldmod.EQ.-1 .AND. .NOT.PRESENT(qc))) THEN
      !Alles Wolkenwasser verdunstet oder wird ignoriert:

      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            rcld(i,k)=z0
           q_liq(i,k)=z0
         END DO
      END DO

   ELSEIF (icldmod.EQ.-1) THEN
    !Wolken sind vorhanden, sind aber an turbulenter Phasenumwandlungen unbeteiligt:

       DO k=k_st, k_en
!DIR$ IVDEP
          DO i=i_st,i_en
             rcld(i,k)=z0
            q_liq(i,k)=qc(i,k)
          END DO
       END DO

   ELSEIF (icldmod.EQ.1 .AND. PRESENT(qc)) THEN
      !Verwendung des vorhandenen skaligen Wolkenwassers:

      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            IF ( qc(i,k) .GT. z0 ) THEN
               rcld(i,k) = z1
            ELSE
               rcld(i,k) = z0
            END IF
            q_liq(i,k)=qc(i,k)
         END DO
      END DO

   ELSE !a special cloud diagnostics needs to be employed

      IF (icldmod.EQ.1) THEN !only grid scale clouds possible
         icldtyp=0
      ELSE !sub grid scale clouds possible
         icldtyp=itype_wcld !use specified type of cloud diagnostics
      END IF

      CALL turb_cloud( khi=khi,                            &
           istart=i_st, iend=i_en, kstart=k_st, kend=k_en, &
           icldtyp=icldtyp,                                &
           prs=prs, t=tet_l(:,khi:), qv=q_h2o(:,khi:),     &
           psf=psf,                                        &
           clcv=rcld, clwc=q_liq(:,khi:) )
   END IF

   !Berechnung der thermodynamischen Hilfsfelder:

   IF (ladjout) THEN !output of adjusted non conservative variables
      qvap => q_h2o
      temp => tet_l
   ELSE !output of conserved variables
      qvap => qst_t
      temp => g_tet
   END IF

   virt => g_h2o

   IF (ladjout .OR. lcaltdv .OR. lcalrho .OR. PRESENT(r_cpd)) THEN
      IF (.NOT.ladjout .AND. icldmod.LE.0) THEN !'temp' and 'vap' equal conserv. vars.
         DO k=k_st, k_en
!DIR$ IVDEP
            DO i=i_st,i_en
               temp(i,k)=tet_l(i,k)
               qvap(i,k)=q_h2o(i,k)
            END DO
         END DO
      ELSEIF (icldmod.GT.0) THEN !'temp' and 'qvap' my be different form conserv. vars.
         DO k=k_st, k_en
!DIR$ IVDEP
            DO i=i_st,i_en
               temp(i,k)=tet_l(i,k)+lhocp*q_liq(i,k) !corrected temperature
               qvap(i,k)=q_h2o(i,k)-      q_liq(i,k) !corrected water vapor
            END DO
         END DO
      END IF
      !Note: In the remaining case "ladjout .AND. icldmod.LE.0" 'temp' and 'qvap'
      !      already point to the conserved variables.
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            virt(i,k)=z1/(z1+rvd_m_o*qvap(i,k)-q_liq(i,k)) !rezipr. virtual factor
            rprs(i,k)=virt(i,k)*prs(i,k)                   !reduced pressure profile
         END DO
      END DO
   END IF

   IF (lcalrho .AND. PRESENT(dens)) THEN
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            dens(i,k)=rprs(i,k)/(r_d*temp(i,k))
         END DO
      END DO
   END IF

   IF (PRESENT(r_cpd)) THEN
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            r_cpd(i,k)=z1+zrcpv*qvap(i,k)+zrcpl*q_liq(i,k) !Cp/Cpd
         END DO
      END DO
   END IF

   IF (.NOT.ladjout) THEN !the potential (liquid water) temperature values are requested
      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            tet_l(i,k)=tet_l(i,k)/exner(i,k) !liquid water pot. temp.
         END DO
      END DO
   END IF

   IF (lcaltdv) THEN

      DO k=k_st, k_en
!DIR$ IVDEP
         DO i=i_st,i_en
            IF (imode_qvsatur.EQ.1) THEN
               qst_t(i,k)=zdqsdt_old( temp(i,k), zqvap_old( zpsat_w( temp(i,k) ), prs(i,k) ) )
                                                                !d_qsat/d_T (old version)
            ELSE
               pdry=(z1-qvap(i,k))*rprs(i,k)                    !partial pressure of dry air
               qst_t(i,k)=zdqsdt( temp(i,k), zqvap( zpsat_w( temp(i,k) ), pdry ) )
                                                                !d_qsat/d_T (new version)
            END IF
         END DO
      END DO
      IF (icldmod.EQ.-1) THEN !no consideration of water phase changes
         DO k=k_st, k_en
!DIR$ IVDEP
            DO i=i_st,i_en
               g_h2o(i,k)=grav*(rvd_m_o*virt(i,k))              !g    -factor of q_h2o-gradient
               g_tet(i,k)=grav*(exner(i,k)/temp(i,k))           !g/tet-factor of tet_l-gradient
            END DO
         END DO
      ELSE !water phase changes are possible
         DO k=k_st, k_en
!DIR$ IVDEP
            DO i=i_st,i_en
               ccov=c_scld*rcld(i,k)/(z1+rcld(i,k)*(c_scld-z1)) !resulting cloud cover
               mcor=ccov*(lhocp/temp(i,k)-(z1+rvd_m_o)*virt(i,k)) &
                        /(z1+qst_t(i,k)*lhocp)                  !moist correction

               g_h2o(i,k)=grav*(rvd_m_o*virt(i,k)+mcor)         !g    -factor of q_h2o-gradient
!modif: previously "z1/teta" was used instead of "exner/temp", which is not equivalent
               g_tet(i,k)=grav*(exner(i,k)/temp(i,k) &          !g/tet-factor of tet_l-gradient
                               -mcor*exner(i,k)*qst_t(i,k))
!modif
            END DO
         END DO
      END IF
   END IF

   !Die thermodynamischen Hilfsgroessen wurden hier unter Beruecksichtigung der diagnostizierten
   !Kondensationskorrektur gebildet, indem die entspr. korrigierten Werte fuer t, qv und ql benutzt wurden.
   !Beachte, dass es zu einer gewissen Inkosistenz kommt, wenn die Dichte nicht angepasst wird!

END SUBROUTINE adjust_satur_equil

!********************************************************************************

SUBROUTINE solve_turb_budgets ( khi, it_s, &
!
   i_st, i_en,                       &
   k_st, k_en,                       &
!
   lssintact, lupfrclim,             &
!
   imode_stke, imode_vel_min,        &
!
   fm2, fh2, ft2, lsm, lsh,          &
#ifdef SCLM
   grd,                              &
#endif
   fcd, tls, tvt, avt, velmin )

INTEGER (KIND=iintegers), INTENT(IN) :: &  !
!
  it_s,                   & !iteration step
!
  khi,                    & !start index of vertical dimension
  i_st, i_en,             & !horizontal start- and end-indices
  k_st, k_en,             & !vertical   start- and end-indices
!
  imode_stke,             & !mode of solving the TKE equation
  imode_vel_min             !mode of determination of minimum turbulent velocity scale

LOGICAL, INTENT(IN)  :: &
!
  lssintact, & !seperate treatment of non-turbulent shear (by scale interaction) requested
  lupfrclim    !enabling an upper limit for TKE-forcing

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, INTENT(IN) :: &
!
  fm2,  &  !squared frequency of mechanical forcing
  fh2,  &  !squared frequency of thermal    forcing
  ft2      !squared frequency of pure turbulent shear

#ifdef SCLM
REAL (KIND=ireals), DIMENSION(:,khi:,:), INTENT(IN) :: &
!
  grd      !vertical gradients (needed only for SUB 'turb_stat' during a SCM-run)
#endif

REAL (KIND=ireals), DIMENSION(:,kcm:), OPTIONAL, INTENT(IN) :: &
!
  fcd      !frequency of small scale canopy drag

REAL (KIND=ireals), DIMENSION(:,khi:), INTENT(INOUT) :: &
!
  tls, &   !turbulent master scale

  lsm, &   !turbulent length scale times value of stability function for momentum       [m]
  lsh      !turbulent length scale times value of stability function for scalars (heat) [m]

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, INTENT(INOUT) :: &
!
  tvt      !inp: turbulent transport of turbulent velocity scale
           !out: total     transport ...

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, OPTIONAL, INTENT(INOUT) :: &
!
  avt      !advective transport of turbulent velocity scale

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
!
  velmin ! location-dependent minimum velocity

INTEGER (KIND=iintegers) :: &
!
  i,i1,j,k, & !running indices for horizontal location loops or vertical levels
  n           !running-index for for selected horizontal points

REAL (KIND=ireals) :: &
!
  val1, val2, wert, fakt, & !auxilary values
  w1, w2, &                 !weighting factors
  q1, q2, q3, &             !turbulent velocity scales
!
  a11, a12, a22, a21, &
  a3, a5, a6, be1, be2, r_b_m, &
!
  d0, d1, d2, d3, d4, d5, d6, &
  d1_rec, d2_rec, &
!
  sm,sh, & !dimensionless stability functions
  det,   & !determinant of 2X2 matrix
!
  gm,gh, & !dimensionless forcing due to shear and buoyancy
  gama,  & !parameter of TKE-equilibrium
  gam0,  & !parameter of TKE-equilibrium (maximal value)
  tim2     !square of turbulent time scale

REAL (KIND=ireals), DIMENSION(i_st:i_en) :: &
!
  l_dis, & !dissipation length scale
  l_frc, & !forcing     length scale
  frc,   & !effective TKE-forcing (acceleration of turbulent motion)
  tvsm     !minimal turbulent velocity scale

REAL (KIND=ireals), DIMENSION(i_st:i_en), TARGET :: &
!
  tvs      !turbulent velocity scale

REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
  tvs0(:),   & !pointer for intermediate turbulent velocity scale
  fm2_e(:,:)   !pointer for the effictive mechanical forcing

REAL (KIND=ireals), DIMENSION(i_st:i_en,0:7), TARGET :: &
!
  dd       !local derived turbulence parameter

LOGICAL :: add_adv_inc, lvar_fcd, rogh_lay, alt_gama, corr

!------------------------------------------------------------------------------------------------------

  add_adv_inc=(PRESENT(avt) .AND. it_s.EQ.it_end)  !add advection increments
                                                   !(only at the last iteration step)
  lvar_fcd=PRESENT(fcd) !array for small-scale canpy drag is present

  alt_gama=(imode_stbcorr.EQ.1 .AND. .NOT.ltkeinp) !alternative gama-Berechnung

  IF (lssintact) THEN !seperate treatment of shear by scale interaction
     fm2_e => ft2 !effective shear is pure turbulent shear
  ELSE
     fm2_e => fm2 !effective shear is total mechanical shear
  END IF

! Stabilitaetskorrektur der turbulenten Laengenskala bei stabilier Schichtung:

  IF (a_stab.GT.z0) THEN
     DO k=k_st,k_en !von oben nach unten
!DIR$ IVDEP
        DO i=i_st, i_en
           wert=a_stab*SQRT(MAX( z0, fh2(i,k)) )
           tls(i,k)=tke(i,k,nvor)*tls(i,k)/(tke(i,k,nvor)+wert*tls(i,k))
        END DO
     END DO
  END IF

  r_b_m=z1/b_m

! Vorbelegung der Felder fuer Turbulenzparameter mit den Werten fuer die freie Atm.:
!DIR$ IVDEP
  DO i=i_st, i_en
     dd(i,0)=d_m
     dd(i,1)=d_1; dd(i,2)=d_2; dd(i,3)=d_3
     dd(i,4)=d_4; dd(i,5)=d_5; dd(i,6)=d_6
     dd(i,7)=rim
  END DO

  IF (PRESENT(velmin) .AND. imode_vel_min.EQ.2) THEN !nutze variirendes 'tvsm'
     DO i=i_st, i_en
        tvsm(i)=tkesecu*velmin(i) !effektiver variiernder Minimalwert fuer 'tvs'
     END DO
  ELSE
     DO i=i_st, i_en
        tvsm(i)=tkesecu*vel_min !effektiver konstanter Minimalwert fuer 'tvs'
     END DO
  END IF

!----------------------------------------------------------------------
  DO k=k_st, k_en !ueber alle Schichten beginnend mit der freien Atm.
!----------------------------------------------------------------------

     rogh_lay=(k.GE.kcm .AND. lvar_fcd) !innerhalb der Rauhigkeitsschicht

     IF (rogh_lay) THEN !innerhalb der Rauhigkeitsschicht
!Achtung: neue Behandlung der Rauhigkeitsschicht einfuehren

!       Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht
!       (ausser Volumenterme, die zur Diffusion gehoeren):

!DIR$ IVDEP
        DO i=i_st, i_en

!          Berechnung der modifizierten Modellparameter:

           wert=z3*tls(i,k)*fcd(i,k)/tke(i,k,nvor)

           dd(i,0)=d_m/(z1+d_m*wert)
           dd(i,1)=a_h/(z1+a_h*wert)
           dd(i,2)=a_m/(z1+z2*a_m*wert)
           dd(i,3)=z9*dd(i,1)
           dd(i,4)=z6*dd(i,2)
           dd(i,5)=z3*(d_h+dd(i,4))
           dd(i,6)=dd(i,3)+z3*dd(i,4)
           dd(i,1)=z1/dd(i,1)
           dd(i,2)=z1/dd(i,2)
           dd(i,7)=z1/(z1+(dd(i,0)-dd(i,4))/dd(i,5))
        END DO

     END IF

     !Gesamter TKE-Antrieb in [m/s2]:

!DIR$ IVDEP
     DO i=i_st, i_en
        l_dis(i)=tls(i,k)*dd(i,0)          !length scale of dissipation
        l_frc(i)=tls(i,k)*dd(i,4)*r_b_m    !length scale of forcing

        val1=lsm(i,k)*fm2_e(i,k); val2=lsh(i,k)*fh2(i,k)
!_______________________________________________________________
!test: Keine Beschraenkung von frc nach unten
        frc(i)=MAX( val1-val2, frcsecu*dd(i,7)*val1 )
! frc(i)=val1-val2
!_______________________________________________________________
     END DO
     !Beachte:
     !Bei "frcsecu=1" wird 'frc' so nach unten beschraenkt, dass die krit. Rf-Zahl (Rf=1-'dd(i,7)')
     ! im TKE-Quellgleichgewicht (ohne Transport und Zeittendenz) nicht ueberschritten wird.
     !Bei "frcsecu<1" wird 'frc' entsrechend weniger eingeschraenkt.
     !Bei "frcsecu=0" wird Rf nicht > 1 und die Summe der TKE-Quellterme wird nicht negativ.
     !Bei "frcsecu<0" kann "Rf>1" und die Summe der TKE-Quellterme auch negativ werden.

!________________________________________________________________
!test: keine Beschraenkung von frc nach oben (Auskommentierung)
!Achtung: Korrektur: Ermoeglicht obere 'frc'-Schranke im Transferschema (wie in COSMO-Version)
!    IF (frcsecu.GT.0) THEN
     IF (frcsecu.GT.0 .AND. lupfrclim) THEN
        DO i=i_st, i_en
           frc(i)=MIN( frc(i), frcsecu*tke(i,k,nvor)**2/l_frc(i)+(z1-frcsecu)*frc(i) )
        END DO
     END IF
!________________________________________________________________
     !Beachte:
     !'tke(i,k,nvor)**2/l_frc(i)' ist der Maximal-Wert fuer 'frc', mit dem die Abweichung
     ! vom TKE-Gleichgewicht den Wert besitzt, der mit der gegebenen 'tke' bei neutraler Schichtung
     ! nicht ueberschritten werden kann.
     !Bei "frcsecu<=0" wird 'frc' nicht nach oben beschraenkt.

     IF (lssintact) THEN !shear forcing by scale interaction needs to be added
!DIR$ IVDEP
        DO i=i_st, i_en
           frc(i)=frc(i)+lsm(i,k)*(fm2(i,k)-ft2(i,k)) !complete forcing
        END DO
     END IF

!    Berechnung der neuen SQRT(2TKE)-Werte:

     IF (.NOT.ltkeinp) THEN


        !Bestimmung der Zwischenwerte:

        IF (add_adv_inc) THEN !explizite Addition der Advektions-Inkremente
           tvs0 => tvs !die effektiven TKE-Vorgaengerwerte werden auf 'tvs' geschrieben
                       !und anschliessend durch die neuen Werte ueberschrieben
!DIR$ IVDEP
           DO i=i_st, i_en
              tvs0(i)=MAX( tvsm(i), tke(i,k,nvor)+avt(i,k)*dt_tke )
           END DO
           !Beachte:
           !Die Advektions-Inkremente werden explizit in 'tvs0' aufgenommen,
           ! damit die spaetere zeitliche Glaettung (tkesmot) nicht die Transportgeschindigkeit
           ! der Advektion beeinflusst (verlangsamt).
           !Da die Advektions-Inkremente auch von benachbarten 'tke'-Werten in horiz. Richtung abhaengen,
           ! werden sie auch nicht in die optionale Iteration gegen einen Gleichgewichtswert einbezogen
           ! und werden somit nur beim letzten Iterations-Schritt addiert.
        ELSE !benutze die Vorgaengerwerte
           tvs0 => tke(:,k,nvor)
        END IF

        !Integration von SQRT(2TKE):

        IF (imode_stke.EQ.1) THEN !1-st (former) type of prognostic solution
!DIR$ IVDEP
           DO i=i_st, i_en
              q1=l_dis(i)*fr_tke
              q2=MAX( z0, tvs0(i)+tvt(i,k)*dt_tke )+frc(i)*dt_tke
              tvs(i)=q1*(SQRT(z1+z4*q2/q1)-z1)*z1d2
           END DO
        ELSEIF (imode_stke.GE.2) THEN !2-nd (new) type of prognostic solution
!DIR$ IVDEP
           DO i=i_st, i_en
              fakt=z1/(z1+z2*dt_tke*tvs0(i)/l_dis(i))
              q1=fakt*(tvt(i,k)+frc(i))*dt_tke
              tvs(i)=q1+SQRT( q1**2+fakt*( tvs0(i)**2+z2*dt_tke*con_m*fm2(i,k) ) )
           END DO
        ELSEIF (imode_stke.EQ.-1) THEN !diagn. solution of station. TKE-equation
!DIR$ IVDEP
           DO i=i_st, i_en
              q2=l_dis(i)*(tvt(i,k)+frc(i))
              q3=l_dis(i)*con_m*fm2(i,k)
              IF (q2.GE.z0) THEN
               ! tvs(i)=SQRT(q2+q3/tvs0(i))
                 tvs(i)=EXP( z1d3*LOG( q2*tvs0(i)+q3 ) )
              ELSEIF (q2.LT.z0) THEN
                 tvs(i)=SQRT( q3/( tvs0(i)-q2/tvs0(i) ))
              ELSE
                 tvs(i)=EXP( z1d3*LOG(q3) )
              END IF
           END DO
        ELSE !standard diagnostic solution
!DIR$ IVDEP
           DO i=i_st, i_en
              tvs(i)=SQRT( l_dis(i)*MAX( frc(i), z0 ) )
           END DO
        END IF

        w1=tkesmot; w2=z1-tkesmot

!DIR$ IVDEP
       DO i=i_st, i_en
          q2=SQRT( l_frc(i)*MAX( frc(i), z0 ) )
!________________________________________________________________
!test: ohne tkesecu*vel_min als untere Schranke
          tke(i,k,ntur)=MAX( tvsm(i), tkesecu*q2, w1*tke(i,k,nvor)+w2*tvs(i) )
!         tke(i,k,ntur)=MAX(          tkesecu*q2, w1*tke(i,k,nvor)+w2*tvs(i) )
!________________________________________________________________

       END DO
       !'q2' ist ein Minimalwert fuer 'tke', mit dem die Abweichung vom TKE-Gleichgewicht
       !den Wert besitzt, der mit dem gegebenen 'frc' bei neutraler Schichtung nicht
       !ueberschritten werden kann.

     END IF

!    Berechnung der neuen stabilitaetsabhangigen Laengenskalen:


     IF (lstfnct) THEN

        w1=stbsmot; w2=z1-stbsmot

!DIR$ IVDEP
        DO i=i_st, i_en

           IF (rogh_lay .OR. (k.EQ.k_st .AND. i.EQ.i_st)) THEN
              d0=dd(i,0)
              d1=dd(i,1); d2=dd(i,2); d3=dd(i,3)
              d4=dd(i,4); d5=dd(i,5); d6=dd(i,6)

              d1_rec=z1/d1; d2_rec=z1/d2

              IF (ltkeinp) THEN
                 gam0=z1/d0 !TKE-equilibrium
              ELSE
                 !Obere Schranke fuer die Abweichung 'gama' vom TKE-Gleichgewicht:
                 gam0=stbsecu/d0+(z1-stbsecu)*b_m/d4
              END IF
           END IF

           tim2=(tls(i,k)/tke(i,k,ntur))**2 !Quadrat der turbulenten Zeitskala

           IF (imode_stbcorr.EQ.2 .OR. (imode_stbcorr.EQ.1 .AND. fh2(i,k).GE.z0)) THEN

              !Loesung bei vorgegebener 'tke':


              !Dimensionslose Antriebe der Turbulenz
              gh=fh2  (i,k)*tim2 !durch thermischen Auftrieb
              gm=fm2_e(i,k)*tim2 !durch mechanische Scherung

              IF (stbsecu.EQ.z0) THEN !Koeffizientenbelegung ohne Praeconditionierung

                 be1=b_h
                 be2=b_m

                 a11=d1+(d5-d4)*gh
                 a12=d4*gm
                 a21=(d6-d4)*gh
                 a22=d2+d3*gh+d4*gm

              ELSE !Koeffizientenbelegung mit Praeconditionierung:

                 wert=z1/tim2

                 a11=d1*wert+(d5-d4)*fh2(i,k)
                 a12=d4*fm2_e(i,k)
                 a21=(d6-d4)*fh2(i,k)
                 a22=d2*wert+d3*fh2(i,k)+d4*fm2_e(i,k)

                 be1=z1/MAX( ABS(a11), ABS(a12) )
                 be2=z1/MAX( ABS(a21), ABS(a22) )

                 a11=be1*a11; a12=be1*a12
                 a21=be2*a21; a22=be2*a22

                 be1=be1*wert*b_h; be2=be2*wert*b_m

              END IF

              det=a11*a22-a12*a21
              sh=be1*a22-be2*a12
              sm=be2*a11-be1*a21

!Achtung: Erweiterung: Ermoeglicht Korrektur, wie in COSMO-Variante
              IF (det.GT.z0 .AND. sh.GT.z0 .AND. sm.GT.z0) THEN !Loesung moeglich
!Achtung: test
            ! IF (imode_stbcorr.EQ.1 .OR. (det.GT.z0 .AND. sh.GT.z0 .AND. sm.GT.z0)) THEN !Loesung moeglich

                 det=z1/det
                 sh=sh*det
                 sm=sm*det

                 IF (imode_stbcorr.EQ.2) THEN
                    corr=(sm*gm-sh*gh.GT.gam0)
                 ELSE
                    corr=.FALSE.
                 END IF
              ELSE
                 corr=.TRUE.
              END IF
           ELSE
              corr=.TRUE.
           END IF

!          Korrektur mit beschraenktem 'gama':

           IF (corr) THEN !Korrektur noetig
!             Loesung bei der die TKE-Gleichung in Gleichgewichtsform eingesetzt ist,
!             wobei 'gama' die Abweichung vom Gleichgewicht darstellt,
!             die so beschraenkt wird, dass immer eine Loesung moeglich ist:

              IF (alt_gama) THEN
!Achtung: Korrektur bei alternativer Behandlung (wie in COMSO-Version)
                 gama=MIN( gam0, frc(i)*tim2/tls(i,k) )
               ! gama=frc(i)*tim2/tls(i,k)
              ELSE
                 gama=gam0
              END IF

              wert=d4*gama

              be1=(b_h-wert)*d1_rec
              be2=(b_m-wert)*d2_rec

              a3=d3*gama*d2_rec
              a5=d5*gama*d1_rec
              a6=d6*gama*d2_rec

              val1=(fm2_e(i,k)*be2+(a5-a3+be1)*fh2(i,k))/(z2*be1)
              val2=val1+SQRT(val1**2-(a6+be2)*fh2(i,k)*fm2_e(i,k)/be1)
              fakt=fh2(i,k)/(val2-fh2(i,k))

              sh=be1-a5*fakt
              sm=sh*(be2-a6*fakt)/(be1-(a5-a3)*fakt)
           END IF

!          Zeitliche Glaettung der Stabilitaetslaenge:
           lsh(i,k)=tls(i,k)*sh*w2+lsh(i,k)*w1
           lsm(i,k)=tls(i,k)*sm*w2+lsm(i,k)*w1

        END DO

     END IF

!------------------------------------------------------------------------------------
#ifdef SCLM
     IF (lsclm .AND. it_s.EQ.it_end) THEN
        CALL turb_stat(k=k, &
                       lsm=  lsm(im,k), lsh=lsh(im,k),      &
                       fm2=fm2_e(im,k), fh2=fh2(im,k),      &
                       tls=  tls(im,k), tvs=tke(im,k,ntur), &
                       tvt=  tvt(im,k), grd=grd(im,k,:),    &
                       d_m=dd(im,0), d_h=d_h, a_m=z1/dd(im,2))
     END IF
#endif
!SCLM--------------------------------------------------------------------------------

!    Sichern der TKE-Dissipation "edr=q**3/l_dis" u.a. als thermische Quelle:

     IF (PRESENT(edr).OR.ltmpcor) THEN
!DIR$ IVDEP
        DO i=i_st, i_en
           ediss(i,k)=tke(i,k,ntur)**3/(dd(i,0)*tls(i,k))
        END DO
        !Achtung: Dies ist der Wert, der im naechsten Prognoseschritt benutzt wird!
     END IF

!----------------------------------------------------------------------
  END DO !k
!----------------------------------------------------------------------

END SUBROUTINE solve_turb_budgets

!********************************************************************************

!------------------------------------------------------------------------------------
#ifdef SCLM
SUBROUTINE turb_stat (lsm, lsh, fm2, fh2, d_m, d_h, a_m, tls, tvs, tvt, grd, k)

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
    k         !vertical level index

REAL (KIND=ireals), INTENT(IN) :: &
!
    lsm,    & !turbulent length scale times value of stability function for momentum       [m]
    lsh,    & !turbulent length scale times value of stability function for scalars (heat) [m]
    fm2,    & !squared forcing frequency for momentum       [1/s2]
    fh2,    & !squared forcing frequency for scalars (heat) [1/s2]
    d_m,    & !dissipation parameter for momentum       [1]
    d_h,    & !dissipation parameter for scalars (heat) [1]
    a_m,    & !return-to-isotropy parameter for momentum [1]
    tls,    & !turbulent length scale   [m]
    tvs,    & !turbulent velocity scale [m/s]
    tvt,    & !turbulent transport of turbulent velocity scale [m/s2]
    grd(:)    !vector of vertical gradients of diffused (conserved) variables [{unit}/m]

REAL (KIND=ireals) ::  &
!
   tts,                  & !turbulent master time scale
   ts_d,                 & !dissipation time scale [s]
   ts_m,                 & !return-to-isotropy time scale for momentum [s]
   tvs2,                 & !tvs**2 [m2/s2]
   tkm, tkh,             & !turbulent diffusion coefficient for momentum and scalars (heat) [m2/s]
   x1, x2, x3,           & !auxilary TKE source terme values [m2/s3]
   cvar(nmvar+1,nmvar+1)   !covariance matrix  [{unit1}*{unit2}]

   cvar(tet_l,tet_l)=d_h*tls*lsh*grd(tet_l)**2
   cvar(h2o_g,h2o_g)=d_h*tls*lsh*grd(h2o_g)**2
   cvar(tet_l,h2o_g)=d_h*tls*lsh*grd(tet_l)*grd(h2o_g)

   tkm=lsm*tvs; tkh=lsh*tvs

   cvar(u_m  ,w_m)=-tkm*grd(u_m)
   cvar(v_m  ,w_m)=-tkm*grd(v_m)
   cvar(tet_l,w_m)=-tkh*grd(tet_l)
   cvar(h2o_g,w_m)=-tkh*grd(h2o_g)

   TKE_SCLM%mod(k)%val=tvs        ; TKE_SCLM%mod(k)%vst=i_cal

   !Achtung: TKE%mod(k)%val zeigt z.Z. noch auf die alte TKE-Zeitstufe.
   !Somit wird also der alte TKE-Wert mit dem neuen ueberschrieben,
   !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
   !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

   tvs2=tvs**2
   tts=tls/tvs
   ts_d=d_m*tts
   ts_m=a_m*tts

   x1=cvar(u_m,w_m)*grd(u_m)
   x2=cvar(v_m,w_m)*grd(v_m)
   x3=-tkh*fh2

   BOYPR%mod(k)%val=x3             ; BOYPR%mod(k)%vst=i_cal
   SHRPR%mod(k)%val=-(x1+x2)       ; SHRPR%mod(k)%vst=i_cal
   DISSI%mod(k)%val=tvs2/ts_d      ; DISSI%mod(k)%vst=i_cal
   TRANP%mod(k)%val=tvs*tvt        ; TRANP%mod(k)%vst=i_cal

   cvar(u_m,u_m)=z1d3*tvs2+ts_m*(-z4*x1+z2*x2-z2*x3)
   cvar(v_m,v_m)=z1d3*tvs2+ts_m*(+z2*x1-z4*x2-z2*x3)
   cvar(w_m,w_m)=tvs2-cvar(u_m,u_m)-cvar(v_m,v_m)

   UWA%mod(k)%val=cvar(u_m  ,w_m)  ; UWA%mod(k)%vst=i_cal
   VWA%mod(k)%val=cvar(v_m  ,w_m)  ; VWA%mod(k)%vst=i_cal
   TWA%mod(k)%val=cvar(tet_l,w_m)  ; TWA%mod(k)%vst=i_cal
   UST%mod(k)%val=SQRT(tkm*SQRT(fm2))
                                     UST%mod(k)%vst=i_cal
 ! LMO%mod(0)%val=-UST%mod(k)%val**3/x3
 !                                   LMO%mod(0)%vst=i_cal

   IF (cvar(u_m,u_m).GE.z0 .AND. cvar(v_m,v_m).GE.z0 .AND. &
       cvar(u_m,u_m)+cvar(v_m,v_m).LE.tvs2) THEN

      UUA%mod(k)%val=cvar(u_m,u_m)  ; UUA%mod(k)%vst=i_cal
      VVA%mod(k)%val=cvar(v_m,v_m)  ; VVA%mod(k)%vst=i_cal
      WWA%mod(k)%val=tvs2-UUA%mod(k)%val-VVA%mod(k)%val
                                      WWA%mod(k)%vst=i_cal
   END IF

   TTA%mod(k)%val=cvar(tet_l,tet_l); TTA%mod(k)%vst=i_cal
   QQA%mod(k)%val=cvar(h2o_g,h2o_g); QQA%mod(k)%vst=i_cal
   TQA%mod(k)%val=cvar(h2o_g,h2o_g); TQA%mod(k)%vst=i_cal

END SUBROUTINE turb_stat
#endif
!SCLM--------------------------------------------------------------------------------

!********************************************************************************

END SUBROUTINE organize_turbdiff

!********************************************************************************
!********************************************************************************

SUBROUTINE turb_cloud ( khi,            &
!
   istart, iend, kstart, kend,          &
!
   icldtyp,                             &
!
   prs, t, qv, qc,                      &
   psf,                                 &
!
   rcld, clcv, clwc )

!------------------------------------------------------------------------------
!
! Description:
!
!     This routine calculates the area fraction of a grid box covered
!     by stratiform (non-convective) clouds.
!     If subgrid-scale condensation is required, an additional
!     saturation adjustment is done.
!
! Method:
!
!     icldtyp = 0 : grid scale diagnosis of local oversaturation
!     Cloud water is estimated by grid scale saturation adjustment,
!     which equals the case "icldtyp = 2" in calse of no variance.
!     So cloud cover is either "0" or "1".
!
!     icldtyp = 1 : empirical diagnosis of local oversaturation
!     The fractional cloud cover clcv is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.

!     icldtyp = 2: statistical diagnosis of local oversaturation
!     A Gaussion distribution is assumed for the local oversaturation
!     dq = qt - qs where qt = qv + ql is the total water content and
!     qs is the saturation specific humidity. Using the standard deviation
!     rcld of this distribution (on input) and the conservative grid-scale
!     quantities qt and tl (liquid water temperature), a corrected liquid
!     water content is determined which contains also the contributions from
!     subgrid-scale clouds. A corresponding cloudiness is also calculated.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT(IN) :: &  ! indices used for allocation of arrays
  khi,               & ! start index of vertical dimension
  istart, iend,      & ! zonal      start and end index
  kstart, kend,      & ! vertical   start and end index
!
  icldtyp              ! type of cloud diagnostics

! Array arguments with intent(in):

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, INTENT(IN) :: &
  prs,   & !atmospheric pressure
  t, qv    !temperature and water vapor content (spec. humidity)

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, OPTIONAL, INTENT(IN) :: &
  qc,    & !cloud water content (if not present, 't' and 'qv' are conserved variables)
  rcld     !standard deviation of local oversaturation

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
  psf      !surface pressure

! Array arguments with intent(out):

REAL (KIND=ireals), DIMENSION(:,khi:), TARGET, INTENT(INOUT) :: &
  clcv     ! stratiform subgrid-scale cloud cover

REAL (KIND=ireals), DIMENSION(:,khi:), INTENT(OUT) :: &
  clwc     ! liquid water content of ""

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,k ! loop indices

REAL (KIND=ireals), PARAMETER :: &
!Achtung: Erweiterung
  asig_max = 0.001_ireals,   & ! max. absoute  standard deviation of local oversaturation
  rsig_max = 0.05_ireals,    & ! max. relative standard deviation of local oversaturation
  zclwfak  = 0.005_ireals,   & ! fraction of saturation specific humidity
  zuc      = 0.95_ireals       ! constant for critical relative humidity

REAL (KIND=ireals), DIMENSION(istart:iend) :: &
  qs, dq, gam, sig


REAL (KIND=ireals), DIMENSION(:,:), POINTER :: &
  qt, tl !total water and liquid water temperature

REAL (KIND=ireals), DIMENSION(istart:iend,kstart:kend), TARGET :: &
  qt_tar, tl_tar

REAL (KIND=ireals), DIMENSION(:,:), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
  sdsd !pointer for standard deviation of local oversaturation

REAL (KIND=ireals) :: &
  pdry, q, & !
  zsigma, zclc0, zq_max !

LOGICAL ::  &
  lsurpres    !surface pressure is present

!------------------------------------------------------------------
!Festlegung der Formelfunktionen fuer die Turbulenzparametrisierung:

!------------ End of header ----------------------------------------

! Note:
! If 'qc' is not present, 't' and 'qv' are assumed to be already liquid water temperature (tl)
!  and total water content (qt), which are conserved for moist conversions
! In case of "kstart=khi=kend", the vertical loop is only for one level.
! Since a surface variable for 'qc' is not used in general, it is not forseen here as well, assuming
!  't' and 'qv' to already contain the conserved variables.
! In this version, atmospheric ice is not included to the adjustment process.

! Begin Subroutine turb_cloud
! ---------------------------

  lsurpres = (PRESENT(psf))

  IF (PRESENT(rcld)) THEN !rcld contains standard deviation
     sdsd => rcld
  ELSE !clcv contains standard deviation and will be overwritten by cloud cover
     sdsd => clcv
  END IF

!Achtung: Korrektur: zsig_max -> epsi
!Achtung: test
! zclc0=MIN( clc_diag, z1-rsig_max )
  zclc0=MIN( clc_diag, z1-epsi )

  zq_max = q_crit*(z1/zclc0 - z1)

  IF (PRESENT(qc)) THEN
     qt => qt_tar
     tl => tl_tar

     DO k = kstart, kend
!DIR$ IVDEP
        DO i = istart, iend
           qt(i,k) = qv(i,k) +       qc(i,k) ! total water content
           tl(i,k) =  t(i,k) - lhocp*qc(i,k) ! liquid water temperature
        END DO
     END DO
  ELSE !'qv' and 't' already contain conserved variables
     qt => qv
     tl => t
  END IF

  !Note: 'qt' and 'tl' are not being changed in the following!

  DO k = kstart, kend

     !Calculation of saturation properties with respect to "t=tl" and "qv=qt":
!DIR$ IVDEP
     DO i = istart, iend
!mod_2011/09/28: zpres=patm -> zpres=pdry {
        IF (imode_qvsatur.EQ.1) THEN
           qs(i) = zqvap_old( zpsat_w( tl(i,k) ), prs(i,k) )       ! saturation mixing ratio (old version)
           gam(i) = z1/( z1 + lhocp*zdqsdt_old( tl(i,k), qs(i) ) ) ! slope factor (from old vers. of d_qsat/d_T)
        ELSE
           pdry=( z1-qt(i,k) )/( z1+rvd_m_o*qt(i,k) )*prs(i,k)     ! part. pressure of dry air
           qs(i) = zqvap( zpsat_w( tl(i,k) ), pdry )               ! saturation mixing ratio (new version)
           gam(i) = z1/( z1 + lhocp*zdqsdt( tl(i,k), qs(i) ) )     ! slope factor (from new vers. of d_qsat/d_T)
        END IF
!mod_2011/09/28: zpres=patm -> zpres=pdry }
         dq(i) = qt(i,k) - qs(i)                                   ! local oversaturation
     END DO

     IF (icldtyp.EQ.1) THEN
        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion
        IF (lsurpres) THEN  !surface pressure is present
!DIR$ IVDEP
           DO i = istart, iend
              zsigma = prs(i,k)/psf(i)

              ! critical relative humidity
              sig(i) = zuc - uc1 * zsigma * ( z1 - zsigma )  &
                                 * ( z1 + uc2*(zsigma-z1d2) )
           END DO
        ELSE !no pressure dependency of critical humidity (only near surface levels)
           DO i = istart, iend
              sig(i) = zuc
           END DO
        END IF
!DIR$ IVDEP
        DO i = istart, iend
           ! cloud cover
           clcv(i,k) = MAX( z0,  &
                            MIN( z1, zclc0 * ((ucl*qt(i,k)/qs(i)-sig(i))/(ucl-sig(i)))**2 ) )

           ! grid-volume water content
           clwc(i,k) = clcv(i,k)*qs(i)*zclwfak
        END DO
!DIR$ IVDEP
        DO i = istart, iend
           IF (dq(i).GT.z0) THEN
              ! adapted grid-volume water content
              clwc(i,k) = clwc(i,k) + (gam(i)*dq(i)-clwc(i,k))*(clcv(i,k)-zclc0)/(z1-zclc0)
           END IF
        END DO
     ELSE !cloud water diagnosis based on saturation adjustment
        IF ( icldtyp.EQ.2 ) THEN
           ! Statistical calculation of cloud cover and cloud water content
           ! using the standard deviation of the local oversaturation

!Achtung: Erweiterung, um COSMO-Version abzubilden
           IF (imode_stadlim.EQ.1) THEN !absolute upper 'sdsd'-limit
!DIR$ IVDEP
              DO i = istart, iend
                 sig(i) = MIN ( asig_max, sdsd(i,k) )
              END DO
           ELSE !relative upper 'sdsd'-limit
!DIR$ IVDEP
              DO i = istart, iend
                 sig(i) = MIN ( rsig_max*qs(i), sdsd(i,k) )
              END DO
           END IF
        ELSE !grid scale adjustment wihtout any variance
!DIR$ IVDEP
           DO i = istart, iend
              sig(i)=z0
           END DO
        END IF
        ! in case of sig=0, the method is similar to grid-scale
        ! saturation adjustment. Otherwise, a fractional cloud cover
        ! is diagnosed.
!DIR$ IVDEP
        DO i = istart, iend
           IF (sig(i).LE.z0) THEN
              IF (dq(i).LE.z0) THEN
                 clcv(i,k) = z0; clwc(i,k) = z0
              ELSE
                 clcv(i,k) = z1; clwc(i,k) = gam(i) * dq(i)
              ENDIF
           ELSE
              q = dq(i)/sig(i)
              clcv(i,k) = MIN ( z1, MAX ( z0, zclc0 * ( z1 + q/q_crit) ) )
              IF (q.LE.-q_crit) THEN !no clouds
                 clwc(i,k) = z0
              ELSEIF (q.GE.zq_max) THEN !grid-scale adjustment
                 clwc(i,k) = gam(i) * dq(i)
              ELSE !statistical adjustment
!Achtung: Erweiterung um COSMO-Variante
                 IF (imode_stadlim.EQ.1) THEN !no limit for cloud water
                    clwc(i,k) = gam(i) * sig(i) * zq_max &
                                       * ( (q + q_crit)/(zq_max + q_crit) )**2
                 ELSE !limiting "sig*zq_max" by 'qt'
                    clwc(i,k) = gam(i) * MIN( qt(i,k), sig(i)*zq_max ) &
                                       * ( (q + q_crit)/(zq_max + q_crit) )**2
                 END IF
              ENDIF
           END IF
        END DO
     ENDIF

  END DO

END SUBROUTINE turb_cloud

!********************************************************************************

PURE FUNCTION zexner (zpres) !Exner-factor

  REAL (KIND=ireals), INTENT(IN) :: zpres
  REAL (KIND=ireals) :: zexner

  zexner = EXP(rdocp*LOG(zpres/p0ref))

END FUNCTION zexner

PURE FUNCTION zpsat_w (ztemp) !satur. vapor pressure over water

  REAL (KIND=ireals), INTENT(IN) :: ztemp
  REAL (KIND=ireals) :: zpsat_w

  zpsat_w = b1*EXP(b2w*(ztemp-b3)/(ztemp-b4w))

END FUNCTION zpsat_w

PURE FUNCTION zqvap_old (zpvap, zpres) !satur. specif. humid. (old version)

  REAL (KIND=ireals), INTENT(IN) :: zpvap, &
                                    zpres !pressure
  REAL (KIND=ireals) :: zqvap_old

  zqvap_old = rdv*zpvap/(zpres-o_m_rdv*zpvap) !old form

END FUNCTION zqvap_old

PURE FUNCTION zqvap (zpvap, zpdry) !satur. specif. humid. (new version)

  REAL (KIND=ireals), INTENT(IN) :: zpvap, &
                                    zpdry !part. pressure of dry air
  REAL (KIND=ireals) :: zqvap

!mod_2011/09/28: zpres=patm -> zpres=pdry {
  zqvap = rdv*zpvap
  zqvap = zqvap/(zpdry+zqvap)
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zqvap

PURE FUNCTION zdqsdt_old (ztemp, zqsat) !d_qsat/d_temp (old version)

  REAL (KIND=ireals), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=ireals) :: zdqsdt_old

  zdqsdt_old=b234w*(1.0_ireals+rvd_m_o*zqsat) & !old form
                  *zqsat/(ztemp-b4w)**2         !old form

END FUNCTION zdqsdt_old

PURE FUNCTION zdqsdt (ztemp, zqsat) !d_qsat/d_tem (new version)

  REAL (KIND=ireals), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=ireals) :: zdqsdt

!mod_2011/09/28: zpres=patm -> zpres=pdry {
  zdqsdt=b234w*(z1-zqsat)*zqsat/(ztemp-b4w)**2
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zdqsdt

ELEMENTAL FUNCTION alpha0_char(u10)

  ! Wind-speed dependent specification of the Charnock parameter based on suggestions by
  ! Jean Bidlot and Peter Janssen, ECMWF
  REAL (KIND=ireals), INTENT(IN) :: u10 ! 10 m wind speed
  REAL (KIND=ireals), PARAMETER  :: a=6.e-3_ireals, b=5.5e-4_ireals, &
                                    c=4.e-5_ireals, d=6.e-5_ireals,  &
                                    u2=17.5_ireals, umax=40.0_ireals
  REAL (KIND=ireals) :: ulim, ured, alpha0_char

  ulim = MIN(u10,umax)
  ured = MAX(0._ireals, ulim-u2)
  alpha0_char = MIN(alpha0_max, MAX (alpha0, a + alpha0_pert + ulim*(b + c*ulim - d*ured)))

END FUNCTION alpha0_char

!********************************************************************************

SUBROUTINE vert_grad_diff ( kcm, kgc,                  &
!
          i_st, i_en, k_tp, k_sf,                      &
!
          dt_var, ivtype, igrdcon, itndcon,            &
          linisetup, lnewvtype,                        &
!++++
          lsflucond, lsfgrduse,                        &
          ldynimpwt, lprecondi, leff_flux,             &
!
          impl_weight,                                 &
!
          rho, rho_s, rho_n, hhl, r_air, tkv, tsv,     &
!
          disc_mom, expl_mom, impl_mom, invs_mom,      &
          diff_dep, diff_mom, invs_fac, scal_fac,      &
!
          dif_tend, cur_prof, eff_flux )
!++++


INTEGER (KIND=iintegers), INTENT(IN) :: &
!
! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------
!
    kcm,          &  ! level index of the upper roughness layer bound
    kgc,          &  ! index of the uppermost level with gradient correction
!
! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------
    i_st, i_en, &    ! start and end index for meridional direction
    k_tp, k_sf, &    ! vertical level indices for top and surface
!
    ivtype,     &    ! index for variable type
!
    igrdcon,    &    ! mode index for effective gradient consideration
                     ! 0: no consider., 1: use only the related profile correction
                     !                  2: use complete effective profile
    itndcon          ! mode index of current tendency  consideration
                     ! 0: no consider., 1: consider curr. tend. in implicit equation only
                     !                  2: add curr. tend. increment to current profile
                     !                  3: add related diffusion correction to current profile

REAL (KIND=ireals), INTENT(IN) :: &
!
    dt_var           ! time increment for vertical diffusion         ( s )

LOGICAL, INTENT(IN) :: &
!
    linisetup,  &    ! calculate initial setup
    lnewvtype,  &    ! calculate setup for a new variable type
    lsflucond,  &    ! surface flux condition
!++++
    lsfgrduse,  &    ! use explicit surface gradient
    ldynimpwt,  &    ! dynamical calculatin of implicit weights
    lprecondi        ! preconditioning of tridiagonal matrix
!++++


LOGICAL, INTENT(INOUT) :: &
!
    leff_flux        ! calculation of effective flux density required

                  ! DIMENSION(ie,ke)
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
    rho              ! air density at full levels                   (Kg/m3)

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl              ! half level height                              (m)

!Attention: Notice the start index of vertical dimension!
                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
!
    tkv              ! turbulent diffusion coefficient               (m/s2)

                  ! DIMENSION(ie)
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
!
    tsv, &           ! turbulent velocity at the surface              (m/s)
!++++
    impl_weight      ! profile of precalculated implicit weights
!++++

                  ! DIMENSION(ie,kcm-1:ke1)
REAL (KIND=ireals), DIMENSION(:,kcm-1:), OPTIONAL, TARGET, INTENT(IN) :: &
!
    r_air            ! log of air containing volume fraction

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
!
    rho_n            ! air density at half levels                   (Kg/m3)

                  ! DIMENSION(ie)
REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
!
    rho_s            ! air density at the surface                   (Kg/m3)

REAL (KIND=ireals), TARGET, INTENT(INOUT) &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
!   Auxilary arrays:
!
!++++
    disc_mom(:,:), &      ! prep_inp calc_inp: discretis. momentum (rho*dz/dt) on var. levels
    expl_mom(:,:), &      ! prep_inp         : diffusion  momentum (rho*K/dz))
                     ! prep_out calc_inp: explicit part of diffusion momentum
    impl_mom(:,:), &      ! prep_out calc_inp: implicit  part of diffusion momentum
    invs_mom(:,:), &      ! prep_out calc_inp: inverted momentum
    invs_fac(:,:), &      ! prep_out calc_inp: inversion factor
    scal_fac(:,:), &      ! prep_out calc_inp: scaling factor due to preconditioning
    diff_dep(:,:)         ! diffusion depth

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: &
!
    diff_mom         ! aux: saved comlete diffusion momentum (only in case of "itndcon.EQ.3")
!++++
                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), INTENT(INOUT) &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
!   Inp-out-variable:
!
    cur_prof(:,:), &      ! inp     : current   variable profile (will be overwritten!)
                     ! calc_inp: corrected variable profile
                     ! calc_out: current   variable profile including tendency increment
    eff_flux(:,:), &      ! inp     : effective gradient
                     ! out     : effective flux density (if "leff_flux=T")
                     ! aux     : downward flux density and related vertical profile increment
    dif_tend(:,:)         ! inp     : current time tendency of variable
                     ! calc_inp: updated   variable profile including tendency increment
                     ! calc_out: updated   variable profile including increments by tendency and diffusion
                     ! out     : pure (vertically smoothed) diffusion tendency

INTEGER (KIND=iintegers) :: &
!
    i,k, &
    k_lw, k_hi          ! vertical index of the lowest and highest level

REAL (KIND=ireals) :: &
!
    fr_var, &           ! 1/dt_var
    tkmin               ! effective minimal diffusion coefficient

REAL (KIND=ireals), DIMENSION(:,:), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
    rhon, rhoh

!-------------------------------------------------------------------------

  k_lw=k_sf-1; k_hi=k_tp+1

  fr_var=z1/dt_var

!++++
!  IF (ivtype.EQ.mom) THEN
!     tkmin=MAX( con_m, tkmmin )
!  ELSE
!     tkmin=MAX( con_h, tkhmin )
!  END IF
!++++

! Initial setup and adoptions for new variable type:

  IF (linisetup .OR. lnewvtype) THEN

     IF (linisetup .OR. .NOT.PRESENT(rho_n)) THEN
        DO k=k_hi,k_lw
!DIR$ IVDEP
           DO i=i_st,i_en
              expl_mom(i,k)=hhl(i,k)-hhl(i,k+1)
           END DO
        END DO
     END IF

     IF (PRESENT(rho_n)) THEN
        rhon => rho_n
     ELSE
!++++
        rhon => invs_mom ; rhoh => rho
!++++
        CALL bound_level_interp( i_st,i_en, k_hi+1,k_lw, &
                                 nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=expl_mom)

!DIR$ IVDEP
        DO i=i_st,i_en
           rhon(i,k_sf)=rho_s(i)
        END DO
     END IF

     IF (linisetup) THEN
!DIR$ IVDEP
        DO i=i_st,i_en
           disc_mom(i,k_hi)=rho(i,k_hi)*expl_mom(i,k_hi)*fr_var
        END DO
        DO k=k_hi+1,k_lw
!DIR$ IVDEP
           DO i=i_st,i_en
              disc_mom(i,k)=rho(i,k)*expl_mom(i,k)*fr_var
              diff_dep(i,k)=z1d2*(expl_mom(i,k-1)+expl_mom(i,k))
           END DO
        END DO
     END IF

     DO k=k_hi+1,k_lw
!DIR$ IVDEP
        DO i=i_st,i_en
!Achtung: Eventuell tkmin-Beschraenkung nur bei VDiff
         ! expl_mom(i,k)=MAX( tkmin, tkv(i,k) ) &
           expl_mom(i,k)=tkv(i,k) &
                        *rhon(i,k)/diff_dep(i,k)
        END DO
     END DO
!DIR$ IVDEP
     DO i=i_st,i_en
!Achtung: Einfuehrung von 'tsv': macht Unterschiede
        expl_mom(i,k_sf)=rhon(i,k_sf)*tsv(i)
        diff_dep(i,k_sf)=tkv(i,k_sf)/tsv(i)
     !  expl_mom(i,k_sf)=rhon(i,k_sf)*tkv(i,k_sf)/diff_dep(i,k_sf)
     !  diff_dep(i,k_sf)=dzs(i)
     END DO
     !Attention: 'tkmin' should be excluded for the surface level 'k_sf'!

     IF (itndcon.EQ.3) THEN
        DO k=k_hi+1,k_sf
!DIR$ IVDEP
           DO i=i_st,i_en
              diff_mom(i,k)=expl_mom(i,k)
           END DO
        END DO
     END IF

!    This manipulation enforces always a complete decoupling from the surface:
     IF (lfreeslip) THEN
        DO i=i_st,i_en
           expl_mom(i,k_sf)=z0
        END DO
     END IF

     CALL prep_impl_vert_diff( lsflucond, ldynimpwt, lprecondi, &
          i_st, i_en, k_tp=k_tp, k_sf=k_sf, &
          disc_mom=disc_mom, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
          invs_fac=invs_fac, scal_fac=scal_fac, impl_weight=impl_weight )

  END IF

! Optional correction of vertical profiles:

  IF (lsfgrduse .AND. igrdcon.NE.2) THEN !effective surface value from effective surface gradient
!DIR$ IVDEP
     DO i=i_st,i_en
        cur_prof(i,k_sf)=cur_prof(i,k_sf-1)-diff_dep(i,k_sf)*eff_flux(i,k_sf)
     END DO
     !Note that this would be done twice in case of "igrdcon.EQ.2"!
  END IF

  IF (igrdcon.EQ.1) THEN !only correction profiles
     DO k=k_sf,kgc,-1
!DIR$ IVDEP
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-cur_prof(i,k)
        END DO
     END DO
!DIR$ IVDEP
     DO i=i_st,i_en
        cur_prof(i,kgc-1)=z0
     END DO
     DO k=kgc,k_sf
!DIR$ IVDEP
        DO i=i_st,i_en
            cur_prof(i,k)=cur_prof(i,k-1)+(cur_prof(i,k)-diff_dep(i,k)*eff_flux(i,k))
        END DO
     END DO
  ELSEIF (igrdcon.EQ.2) THEN !effektive total profile
     DO k=kgc,k_sf
!DIR$ IVDEP
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-diff_dep(i,k)*eff_flux(i,k)
        END DO
     END DO
  END IF

  IF (itndcon.EQ.3) THEN !add diffusion correction from tendency to current profile
     !Related downward flux densities:
!DIR$ IVDEP
     DO i=i_st,i_en
        eff_flux(i,2)=-dif_tend(i,1)*disc_mom(i,1)*dt_var
     END DO
     DO k=k_hi+1,k_lw
!DIR$ IVDEP
        DO i=i_st,i_en
           eff_flux(i,k+1)=eff_flux(i,k)-dif_tend(i,k)*disc_mom(i,k)*dt_var
        END DO
     END DO
     !Virtual total vertical increment:
     DO k=k_hi+1,k_sf
!DIR$ IVDEP
        DO i=i_st,i_en
           eff_flux(i,k)=eff_flux(i,k)/diff_mom(i,k)+(cur_prof(i,k-1)-cur_prof(i,k))
        END DO
     END DO
     !Related corrected profile:
     DO k=k_hi+1,k_sf
!DIR$ IVDEP
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-eff_flux(i,k)
        END DO
     END DO
  END IF

  IF (itndcon.GE.1) THEN !calculate updated profile by adding tendency increment to current profile
     DO k=k_hi,k_lw
!DIR$ IVDEP
        DO i=i_st,i_en
            dif_tend(i,k)=cur_prof(i,k)+dif_tend(i,k)*dt_var
        END DO
     END DO
!DIR$ IVDEP
     DO i=i_st,i_en
        dif_tend(i,k_sf)=cur_prof(i,k_sf)
     END DO
  END IF

  IF (kcm.LE.k_lw) THEN
     leff_flux = .TRUE.
  END IF

  !Final solution of the semi-implicit diffusion equation:

  !Note:
  !'cur_prof' is the current profile including the gradient correction (if "igrdcon>0")
  !           including the virtual gradient correction of an explicit time tendency (if "itndcon=3").
  !'dif_tend' is only used, if "itndcon>0" and contains 'cur_prof' updated by the curr. tend. incr..

  CALL calc_impl_vert_diff ( lsflucond, lprecondi, leff_flux, itndcon, &
       i_st, i_en ,k_tp=k_tp, k_sf=k_sf, &
       disc_mom=disc_mom, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
       invs_fac=invs_fac, scal_fac=scal_fac, cur_prof=cur_prof, upd_prof=dif_tend, eff_flux=eff_flux )

  !Note:
  !'cur_prof' now contains the current profile updated by the current tendency increment (if "itndcon>0").
  !'dif_tend' now contains the final updated profile including vertical diffusion.

  !Calculation of time tendencies for pure vertical diffusion:

  DO k=k_hi,k_lw
!DIR$ IVDEP
     DO i=i_st,i_en
        dif_tend(i,k)=(dif_tend(i,k)-cur_prof(i,k))*fr_var
     END DO
  END DO

! Volume correction within the roughness layer:

  IF (PRESENT(r_air)) THEN
     DO k=kcm,k_lw  !r_air-gradient within the roughness layer
!DIR$ IVDEP
        DO i=i_st,i_en
           cur_prof(i,k)=(r_air(i,k)-r_air(i,k+1))/(dt_var*disc_mom(i,k))
        END DO
     END DO
     DO k=kcm,k_lw  !within the roughness layer
!DIR$ IVDEP
        DO i=i_st,i_en
           dif_tend(i,k)=dif_tend(i,k)                        &
                        +z1d2*(eff_flux(i,k+1)+eff_flux(i,k)) & !flux on full level
                             *cur_prof(i,k)                     !r_air-gradient
        END DO
     END DO
  END IF

! Optional conservative vertical smoothing of tendencies:

  IF (tndsmot.GT.z0) THEN
     CALL vert_smooth ( &
          i_st, i_en, k_tp=k_tp, k_sf=k_sf, &
          disc_mom=disc_mom, cur_tend=dif_tend, vertsmot=tndsmot )
  END IF

END SUBROUTINE vert_grad_diff

!********************************************************************************

!++++
SUBROUTINE prep_impl_vert_diff ( lsflucond, ldynimpwt, lprecondi, &
!
   i_st,i_en, k_tp, k_sf, &
!
   disc_mom, expl_mom, impl_mom, invs_mom, invs_fac, scal_fac, impl_weight )

!Achtung: l-Schleifen -> k-Schleifen
!Achtung: Vorzeichenwechsel fur impl. momentum ist uebersichtlicher

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for first  horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

LOGICAL, INTENT(IN) :: &
!
   lsflucond, & !flux condition at the surface
   ldynimpwt, & !dynamical calculatin of implicit weights
   lprecondi    !preconditioning of tridiagonal matrix

REAL (KIND=ireals), INTENT(IN) :: &
!
   impl_weight(:),& !profile of precalculated implicit weights
!
   disc_mom(:,:)    !discretis momentum (rho*dz/dt) on variable levels

REAL (KIND=ireals), INTENT(INOUT) :: &
!
   expl_mom(:,:)    !inp: diffusion momentum (rho*K/dz)
                    !out: explicit part of diffusion momentum

REAL (KIND=ireals), INTENT(OUT) :: &
!
   impl_mom(:,:), & !scaled implicit part of diffusion momentum
   invs_mom(:,:), & !inverted momentum
   invs_fac(:,:), & !inversion vactor
   scal_fac(:,:)    !scaling factor due to preconditioning

INTEGER (KIND=iintegers) :: &
!
   i,k, &           !horizontal and vertical coordiante indices
   m                !level increment dependent on type of surface condition

!--------------------------------------------------------------------------------------

   IF (lsflucond) THEN !a surface-flux condition
      m=2
   ELSE
      m=1
   END IF

!  Implicit and explicit weights:

   IF (ldynimpwt) THEN !dynamical determination of implicit weights
      DO k=k_tp+2, k_sf+1-m
!DIR$ IVDEP
         DO i=i_st, i_en
            impl_mom(i,k)=expl_mom(i,k) &
                            ! *MAX(MIN(expl_mom(i,k)/impl_mom(i,k), impl_s), impl_t)
                              *MAX(impl_s-z1d2*impl_mom(i,k)/expl_mom(i,k), impl_t)
!*z1
         END DO
      END DO
   ELSE !use precalculated implicit weights
      DO k=k_tp+2, k_sf+1-m
!DIR$ IVDEP
         DO i=i_st, i_en
!Achtung:
            impl_mom(i,k)=expl_mom(i,k)*impl_weight(k)
!impl_mom(i,k)=expl_mom(i,k)*1.00_ireals
         END DO
      END DO
   END IF

!Achtung: Korrektur: Um Konzentrations-Randbedingung richtig abzubilden, muss
!'expl_mom' bei "k=k_sf" das gesamte Diffusions-Moment enthalten:
!  DO k=k_tp+2, k_sf
   DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
      DO i=i_st, i_en
         expl_mom(i,k)=expl_mom(i,k)-impl_mom(i,k)
      END DO
   END DO
   !Notice that 'expl_mom' still contains the whole diffusion momentum at level 'k_sf'!

!  Inverse momentum vector:

   IF (lprecondi) THEN !apply symmetric preconditioning of tridiagonal matrix
      k=k_tp+1
!DIR$ IVDEP
      DO i=i_st, i_en
         scal_fac(i,k)=z1/SQRT(disc_mom(i,k)+impl_mom(i,k+1))
         invs_mom(i,k)=z1
      END DO
      DO k=k_tp+2, k_sf-m
!DIR$ IVDEP
         DO i=i_st, i_en
            scal_fac(i,k)=z1/SQRT(disc_mom(i,k)+impl_mom(i,k)+impl_mom(i,k+1))
         END DO
      END DO
      DO k=k_sf+1-m, k_sf-1 !only for a surface-flux condition and at level k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            scal_fac(i,k)=z1/SQRT(disc_mom(i,k)+impl_mom(i,k))
         END DO
      END DO
      DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            impl_mom(i,k)=scal_fac(i,k-1)*scal_fac(i,k)*impl_mom(i,k)
            invs_fac(i,k)=invs_mom(i,k-1)*impl_mom(i,k)
            invs_mom(i,k)=z1/( z1-invs_fac(i,k)*impl_mom(i,k) )
         END DO
      END DO
   ELSE !without preconditioning
!DIR$ IVDEP
      k=k_tp+1
      DO i=i_st, i_en
         invs_mom(i,k)=z1/(disc_mom(i,k)+impl_mom(i,k+1))
      END DO
      DO k=k_tp+2, k_sf-m
!DIR$ IVDEP
         DO i=i_st, i_en
            invs_fac(i,k)=invs_mom(i,k-1)*impl_mom(i,k)
            invs_mom(i,k)=z1/( disc_mom(i,k)+impl_mom(i,k+1) &
                              +impl_mom(i,k)*(z1-invs_fac(i,k)) )
         END DO
      END DO
      DO k=k_sf+1-m, k_sf-1 !only for a surface-flux condition and at level k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            invs_fac(i,k)=invs_mom(i,k-1)*impl_mom(i,k)
            invs_mom(i,k)=z1/( disc_mom(i,k) &
                              +impl_mom(i,k)*(z1-invs_fac(i,k)) )
         END DO
      END DO
   END IF

END SUBROUTINE prep_impl_vert_diff

!********************************************************************************

SUBROUTINE calc_impl_vert_diff ( lsflucond, lprecondi, leff_flux, itndcon, &
!
   i_st,i_en, k_tp, k_sf, &
!
   disc_mom, expl_mom, impl_mom, invs_mom, invs_fac, scal_fac, cur_prof, upd_prof, eff_flux)

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for first  horizontal dimension
   k_tp,k_sf, & !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1
!
   itndcon      ! mode index of current tendency  consideration
                ! 0: no consider., (+/-)1: consider curr. tend. in implicit equation only
                !                  (+/-)2: add curr. tend. increment to current profile
                !                  (+/-)3: add related diffusion correction to current profile
                ! for "itndcon>0", 'cur_prof' will be overwritten by input of 'upd_prof'
LOGICAL, INTENT(IN) :: &
!
   lsflucond, & !flux condition at the surface
   lprecondi, & !preconditioning of tridiagonal matrix
   leff_flux    !calculation of effective flux densities required

REAL (KIND=ireals), INTENT(IN) :: &
!
   expl_mom(:,:), & !(explicit part of) diffusion momentum (rho*K/dz)
   impl_mom(:,:), & !(scaled) implicit part of diffusion momentum
   disc_mom(:,:), & !discretis momentum (rho*dz/dt)
!
   invs_mom(:,:), & !inverted momentum
   invs_fac(:,:), & !inversion factor
   scal_fac(:,:)    !scaling factor due to preconditioning

REAL (KIND=ireals), TARGET, INTENT(INOUT) &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
   cur_prof(:,:), & !inp: current vertical variable profile (including gradient corrections)
                    !out: current vertical variable profile (including tendency increment, if "itndcon>=1")
   upd_prof(:,:), & !inp: if "|itndcon|>=1": updated vertical variable profile (includ. tendency increment)
                    !aux: interim solution of diffusion equation
                    !out: updated vertical variable profile by vertical diffusion
   !"itndcon>= 1": 'cur_prof' as output equals 'upd_prof' as input
   !"itndcon<=-1": 'cur_prof' keeps input profile
   !"itndcon = 0": 'upd_prof' is not used as input
!
   eff_flux(:,:)    !out: effective flux density belonging to (semi-)implicit diffusion
                    !     (positiv downward and only, if "leff_flux=.TRUE.")
                    !aux: explicit flux density (positiv upward) and
                    !     full right-hand side of the semi-implicit diffusion equation

INTEGER (KIND=iintegers) :: &
!
   i,k              !horizontal and vertical coordiante indices

REAL (KIND=ireals), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
   old_prof(:,:), & !old variable profile value>
   rhs_prof(:,:)    !profile used at RHS of tridiag. system

!-----------------------------------------------------------------------

!  Preparation:

   SELECT CASE ( ABS(itndcon) ) !discriminate modes of current tendency consideration
   CASE (0) ! no consideration of current tendency
     old_prof => cur_prof
     rhs_prof => cur_prof
   CASE (1) ! consideration of current tendency in implicit part only
     old_prof => upd_prof
     rhs_prof => cur_prof
   CASE (2) ! consideration of current tendency in all parts
     old_prof => upd_prof
     rhs_prof => upd_prof
   CASE (3) ! current profile already contains related gradient correction
     old_prof => cur_prof
     rhs_prof => cur_prof
   END SELECT

   DO k=k_tp+2, k_sf
!DIR$ IVDEP
      DO i=i_st, i_en
         eff_flux(i,k) = expl_mom(i,k) * ( rhs_prof(i,k  ) - rhs_prof(i,k-1) )
      END DO
   END DO
   !Notice that 'expl_mom(i,k_sf)' still contains the whole diffusion momentum!

!Achtung: Korrektur: Richtige Behandlung der unteren Konzentrations-Randbedingung
   IF (.NOT.lsflucond) THEN !only for a surface-concentration condition and
                            !level "k_sf-1" is treated (semi-)implicitly in all
      k=k_sf
!DIR$ IVDEP
      DO i=i_st, i_en
         eff_flux(i,k) = eff_flux(i,k) + impl_mom(i,k  ) * rhs_prof(i,k-1)
      END DO
   !  Note: At level 'k_sf' 'impl_mom' still contains the implicit part without scaling,
   !        and it vanishes at all in case of "lsflucond=T" (surface-flux condition)!
   END IF

!  Resultant right-hand side flux:

   k=k_tp+1
!DIR$ IVDEP
   DO i=i_st, i_en
      eff_flux(i,k  ) = disc_mom(i,k  ) * old_prof(i,k  ) + eff_flux(i,k+1)
   END DO
!  Note: Zero flux condition just below top level.
   DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
      DO i=i_st, i_en
         eff_flux(i,k  ) = disc_mom(i,k  ) * old_prof(i,k  ) + eff_flux(i,k+1) - eff_flux(i,k  )
      END DO
   END DO

   IF (lprecondi) THEN !preconditioning is active
      DO k=k_tp+1, k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            eff_flux(i,k  ) = scal_fac(i,k  ) * eff_flux(i,k  )
         END DO
      END DO
   END IF

!  Save updated profiles (including explicit increments of current tendencies):

   IF (itndcon.GT.0) THEN !consideration of explicit tendencies
      DO k=k_tp+1, k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            cur_prof(i,k) = upd_prof(i,k)
         END DO
      END DO
   END IF

!  Forward substitution:

   k=k_tp+1
!DIR$ IVDEP
   DO i=i_st, i_en
      upd_prof(i,k  ) = eff_flux(i,k  ) * invs_mom(i,k  )
   END DO
   DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
      DO i=i_st, i_en
         upd_prof(i,k  ) = ( eff_flux(i,k  ) + impl_mom(i,k  ) * upd_prof(i,k-1) ) * invs_mom(i,k  )
      END DO
   END DO

!  Backward substitution:

   DO k=k_sf-2, k_tp+1, -1
!DIR$ IVDEP
      DO i=i_st, i_en
         upd_prof(i,k) = upd_prof(i,k  ) + invs_fac(i,k+1) * upd_prof(i,k+1)
      END DO
   END DO

   IF (lprecondi) THEN !preconditioning is active
      DO k=k_tp+1, k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            upd_prof(i,k  ) = scal_fac(i,k  ) * upd_prof(i,k  )
         END DO
      END DO
   END IF

   !Note:
   !'cur_prof(i,k  )' is now the current profile including optional explicit tendencies
   !'upd_prof(i,k  )' contains the new profile value after diffusion

!  Effective flux density by vertical integration of diffusion tendencies:

   IF (leff_flux) THEN
      k=k_tp+1
      DO i=i_st, i_en
        eff_flux(i,k) = z0 !upper zero-flux condition
      END DO

      DO k=k_tp+2, k_sf
!DIR$ IVDEP
         DO i=i_st, i_en
            eff_flux(i,k  ) = eff_flux(i,k-1) + (cur_prof(i,k-1) - upd_prof(i,k-1)) * disc_mom(i,k-1)
         END DO
      END DO
      !Note:
      !'eff_flux' is the vertical flux density of pure diffusion now (positive downward)
   END IF

END SUBROUTINE calc_impl_vert_diff
!++++

!********************************************************************************

SUBROUTINE vert_smooth( i_st, i_en, k_tp, k_sf, &
                        cur_tend, disc_mom, vertsmot, smotfac )

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

REAL (KIND=ireals), INTENT(IN) :: &
!
   disc_mom(:,:), & !discretised momentum (rho*dz/dt) on variable levels
!
   vertsmot         !vertical smoothing factor

REAL (KIND=ireals), INTENT(IN), OPTIONAL :: smotfac(:)

REAL (KIND=ireals), INTENT(INOUT) :: &
!
   cur_tend(:,:)    !inp: current  vertical tendency profile
                    !out: smoothed vertical tendency profile
   !modifications takes place at levels k_tp+1 until k_sf-1.

INTEGER (KIND=iintegers) :: &
!
   i,k, &           !horizontal and vertical coordiante indices
   j0,j1,j2         !altering indices "1" and "2"

REAL (KIND=ireals) :: &
!
   sav_tend(SIZE(cur_tend,1),2), & !saved tendencies of current level
!
   versmot(SIZE(cur_tend,1)), & !smoothing factor of current level
   remfact(SIZE(cur_tend,1))    !remaining factor of current level

!----------------------------------------------------------------------------
   IF (imode_frcsmot.EQ.2 .AND. PRESENT(smotfac)) THEN
     versmot(i_st:i_en) = vertsmot*smotfac(i_st:i_en)
   ELSE
     versmot(:) = vertsmot
   ENDIF

   remfact(i_st:i_en)=z1-versmot(i_st:i_en)
   k=k_tp+1
   j1=1; j2=2
!DIR$ IVDEP
   DO i=i_st,i_en
      sav_tend(i,j1)=cur_tend(i,k)
      cur_tend(i,k) =remfact(i)* cur_tend(i,k)                   &
                    +versmot(i)* cur_tend(i,k+1)*disc_mom(i,k+1) &
                                                /disc_mom(i,k)
   END DO

   remfact(i_st:i_en)=z1-z2*versmot(i_st:i_en)
   DO k=k_tp+2, k_sf-2
      j0=j1; j1=j2; j2=j0
!DIR$ IVDEP
      DO i=i_st,i_en
         sav_tend(i,j1)=cur_tend(i,k)
         cur_tend(i,k) =remfact(i)* cur_tend(i,k)                    &
                       +versmot(i)*(sav_tend(i,j2) *disc_mom(i,k-1)  &
                                   +cur_tend(i,k+1)*disc_mom(i,k+1)) &
                                                   /disc_mom(i,k)
      END DO

   END DO

   remfact(i_st:i_en)=z1-versmot(i_st:i_en)
   k=k_sf-1
   j2=j1
!DIR$ IVDEP
   DO i=i_st,i_en
      cur_tend(i,k) =remfact(i)* cur_tend(i,k)                   &
                    +versmot(i)* sav_tend(i,j2) *disc_mom(i,k-1) &
                                                /disc_mom(i,k)
   END DO

END SUBROUTINE vert_smooth

!********************************************************************************

SUBROUTINE bound_level_interp (i_st,i_en, k_st, k_en, &
                               nvars, pvar, depth, rpdep, auxil)

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_st,k_en, & !start and end index for vertical dimension
!
   nvars        !number of variables to be interpolated onto bound. levels

TYPE (varprf) :: pvar(:) !used variable profiles to be interpolated

REAL (KIND=ireals), TARGET, OPTIONAL, INTENT(IN) :: &
!
   depth(:,:)    !layer depth used for interpolation

REAL (KIND=ireals), TARGET, OPTIONAL, INTENT(INOUT) :: &
!
   rpdep(:,:), & !reciprocal depth of two consecutive layers
   auxil(:,:)    !target for layer depth, when 'depth' contains boundary level height

INTEGER (KIND=iintegers) :: i,k,n

REAL (KIND=ireals), DIMENSION(:,:), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     :: &
!
   usdep     !used  layer depth

LOGICAL :: ldepth, lrpdep, lauxil

   ldepth=PRESENT(depth) !'depth' has to be used
   lrpdep=PRESENT(rpdep) !'rpdep' can be used
   lauxil=PRESENT(auxil) !'depth' contains boundary level height

   IF (ldepth) THEN !depth weighted interpolation
      IF (lauxil) THEN !layer depth needs to be calculated
         usdep => auxil
         DO k=k_en, k_st-1, -1
!DIR$ IVDEP
            DO i=i_st, i_en
               usdep(i,k)=depth(i,k)-depth(i,k+1)
            END DO
         END DO
      ELSE
         usdep => depth
      END IF
      IF (lrpdep) THEN !precalculation of the reciprocal layer depth
         DO k=k_en, k_st, -1
!DIR$ IVDEP
            DO i=i_st, i_en
               rpdep(i,k)=z1/(usdep(i,k-1)+usdep(i,k))
            END DO
         END DO
         DO n=1, nvars
            DO k=k_en, k_st, -1
!DIR$ IVDEP
               DO i=i_st, i_en
                  pvar(n)%bl(i,k)=(pvar(n)%ml(i,k  )*usdep(i,k-1)  &
                                    +pvar(n)%ml(i,k-1)*usdep(i,k))   &
                                    *rpdep(i,k)
               END DO
            END DO
         END DO
      ELSE !no precalculation
         DO n=1, nvars
            DO k=k_en, k_st, -1
!DIR$ IVDEP
               DO i=i_st, i_en
                  pvar(n)%bl(i,k)=(pvar(n)%ml(i,k  )*usdep(i,k-1)  &
                                    +pvar(n)%ml(i,k-1)*usdep(i,k))   &
                                    /(usdep(i,k-1)+usdep(i,k))
               END DO
            END DO
         END DO
      END IF
   ELSE !inverse of main level interpolation
      DO n=1, nvars
         DO k=k_st, k_en
!DIR$ IVDEP
            DO i=i_st, i_en
               pvar(n)%bl(i,k)=z2*pvar(n)%ml(i,k)-pvar(n)%ml(i,k+1)
            END DO
         END DO
      END DO
   END IF

END SUBROUTINE bound_level_interp

!********************************************************************************
!********************************************************************************

#ifdef __COSMO__
 END MODULE turbulence_turbdiff
#endif
#ifdef __ICON__
 END MODULE src_turbdiff
#endif
