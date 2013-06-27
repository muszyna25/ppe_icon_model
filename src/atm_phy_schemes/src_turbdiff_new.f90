!+ Source module for computing diffusion and transfer coefficients 
!+ and implicit vertical diffusion:
!-------------------------------------------------------------------------------

 MODULE src_turbdiff_new

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
! History:
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
!              2010/09/30 Matthias Raschendorfer
!  Substituting control parameter 'itype_diag_t2m' by the already present parameter 'itype_synd'.
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
!  Introduction of some minor formal modification of some parts of the code mainly to achiev better vectorization
!   (results will be modyfied only because of numerical effects).
!  Introduction of OPTIONAL extra output fields 'tket_sso' and 'tket_hshr'.
!              2011/03/23 Matthias Raschendorfer
!  Correction of two bugs introduced together with the last modifications in SUB 'turbtran' 
!   (related to SUB 'diag_level' and the 'tet_l'-gradient used for the flux output).
!  Substitution of run time allocations because of bad performance on some computers.
!              2011/06/17 Matthias Raschendorfer
!  Removing a bug related to the saving of the saturation fraction that had been introduced by the ICON modifications.
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
!  Introduction of SUBs 'prep_impl_vert_diff' and 'calc_impl_vert_diff' substituting the code segment
!   for implicit vertical diffusion calculation for TKE or all other variables, if 'imode_turb' is one of 3 or 4.
!   These SUBs allow the calculation of semi implicit diffusion for arbitrary variables defined on 
!   full- or  half-levels.
!   Diffusion tendency i(rather than flux-density) of pot. temperature is converted in that of ordinary temperatur.
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
!  explicit time tendencies numerically, including smoothing options and a preconditioning the linear system.
!  Parameters 'imode_tran' and 'imode_turb' define whether the TKE-equation is solved in a dianostic (1)
!  or a prognostic (2) form.
!  New parameter 'lsflcnd' is a flag for using a lower flux condition for vertical diffusion.
!  Introduction the flags 'lturatm', 'ltursrf', 'lmomdif', 'lscadif' arranging what tasks of 'organize_turbdiff'
!  are required.
!  Introduction of 'ltkeshs' and removing the case "itype_sher.EQ.3".
!
! Code Description:
! Language: Fortran 90 
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".


! Modules used:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
#ifdef __COSMO__
USE data_parallel,      ONLY :  &
!
    num_compute,     & ! number of compute PEs
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    ncomm_type,      & ! type of communication
    icomm_cart,      & ! communicator for the virtual cartesian topology
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

USE environment, ONLY :  &
!
    exchg_boundaries       ! performs the boundary exchange

USE data_runcontrol , ONLY :   &

    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim           ! 2 dimensional runs
#endif
!__COSMO__----------------------------------------------------------------------

#ifdef __COSMO__
USE data_turbdiff, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
!   get_param,    & !get turbulence parameters being extra converted
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
    rcpv,         & ! cp_d/cp_v - 1
    rcpl,         & ! cp_d/cp_l - 1
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

!-------------------------------------------------------------------------------

#ifdef __COSMO__
USE data_turbdiff, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Parameters for turbulent diffusion and surface-to-atmosphere transfer:
! ------------------------------------------------------------------------
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
    alpha0,       & ! Charnock-parameter
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
    a_stab          ! factor for stability correction of turbulent length scale
!
#ifdef __COSMO__
USE data_turbdiff, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
    tkhmin,       & ! minimal diffusion coefficients for heat
    tkmmin,       & ! minimal diffusion coefficients for momentum
!
    clc_diag,     & ! cloud cover at saturation in statistical cloud diagnostic
    q_crit,       & ! critical value for normalized over-saturation
    c_scld,       & ! factor for liquid water flux density in sub grid scale clouds
!
    impl_s,       & ! implicit weight near the surface (maximal value)
    impl_t,       & ! implicit weight near top of the atmosphere (minimal value)

    epsi,         & ! relative limit of accuracy for comparison of numbers
    frcsmot,      & ! vertical smoothing factor for TKE forcing
    tndsmot,      & ! vertical smoothing factor for diffusion tendencies
    tkesmot,      & ! time smoothing factor for TKE
    stbsmot,      & ! time smoothing factor for stability function
    frcsecu,      & ! security factor for TKE-forcing       (<=1)
    tkesecu,      & ! security factor in  TKE equation      (out of [0; 1])
    stbsecu,      & ! security factor in stability function (out of [0; 1])
!
    it_end,       & ! number of initialization iterations (>=0)
!
! Flake parameters:
!
    h_Ice_min_flk      ! Minimum ice thickness [m]
!
#ifdef __COSMO__
USE data_turbdiff, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
! Switches controlling turbulent diffusion:
! ------------------------------------------
!
    itype_tran,   & ! type of surface-atmosphere transfer
    imode_tran,   & ! mode of TKE-equation in transfer scheme
    icldm_tran,   & ! mode of cloud representation in transfer parametr.
!
    imode_turb,   & ! mode of TKE-equation in turbulence scheme
    icldm_turb,   & ! mode of cloud representation in turbulence parametr.
    itype_sher,   & ! type of shear production for TKE
!
    ltkesso,      & ! consider SSO-wake turbulence production of TKE
    ltkecon,      & ! consider convective buoyancy production of TKE
    lexpcor,      & ! explicit corrections of the implicit calculated
!                   ! turbulent diffusion
    ltmpcor,      & ! consideration of thermal TKE-sources in the enthalpy budget
    lprfcor,      & ! using the profile values of the lowest main level instead of
!                   ! the mean value of the lowest layer for surface flux calulations
    lnonloc,      & ! nonlocal calculation of vertical gradients used
!                   ! for turbulent diffusion
    lcpfluc,      & ! consideration of fluctuations of the heat capacity of air
!
    lsflcnd,      & ! lower flux condition for vertical diffusion calculation
!
! Switches controlling other physical parameterizations:
!
    itype_wcld,   & ! type of water cloud diagnosis 
    itype_synd,   & ! type of diagnostics of synoptical near surface variables
!
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecast with lake model FLake
    lsso,         & ! SSO-Scheme is active
    lconv           ! confection scheme is active

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
!   z2d3=z2/z3     ,&
    z3d2=z3/z2

REAL (KIND=ireals), PARAMETER :: ustmin = 1.e-8_ireals

REAL (KIND=ireals) :: xx

INTEGER (KIND=iintegers) :: &
!
    istat=0, ilocstat=0

LOGICAL :: lerror=.FALSE.

TYPE modvar !model variable
     REAL (KIND=ireals), POINTER :: av(:,:) => NULL() !atmospheric values
     REAL (KIND=ireals), POINTER :: sv(:)   => NULL() !surface     values (concentration of flux density)
     REAL (KIND=ireals), POINTER :: at(:,:) => NULL() !atmospheric time tendencies
     LOGICAL                     :: fc                  !surface values are flux densities
END TYPE modvar

TYPE turvar !turbulence variables
     REAL (KIND=ireals), POINTER :: tkv(:,:) => NULL() !turbulent coefficient for vert. diff.
     REAL (KIND=ireals), POINTER :: dzs(:)   => NULL() !effective surface layer depth
END TYPE turvar

TYPE varprf !variable profile
     REAL (KIND=ireals), POINTER :: bl(:,:) !variable at boundary model levels 
     REAL (KIND=ireals), POINTER :: ml(:,:) !variable at main     model levels
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
           i_stp, i_enp, &
!
           hhl, fr_land, plcov, & 
!
           lai, sai, tai, eai, &
           d_pat, h_can,  &
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
! SHistory:
! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 Matthias Raschendorfer
!  Initial release
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
    ie,  & ! number of grid points in zonal      direction 
    ke,  & ! number of main model levels (start index is 1)
    ke1, & ! number of half model levels (start index is 1) 
    i_stp, i_enp  ! start and end index including the model boundary lines

INTEGER (KIND=iintegers), TARGET, INTENT(INOUT) :: &
!
    kcm             ! index of the lowest model layer higher than the canopy

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
    hhl             ! height of model half levels                   ( m )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
#endif
!
! External parameter fields:
! ----------------------------
    fr_land         ! land portion of a grid point area             ( 1 )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
#endif
!
    plcov,        & ! fraction of plant cover                       ( 1 )
    lai             ! leaf area index                               ( 1 )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: &
#endif
!
    sai,          & ! surface area index                            ( 1 )
    tai,          & ! transpiration area index                      ( 1 )
    eai             ! (evaporative) earth area index                ( 1 )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: &
#endif
!
    h_can,        & ! hight of the vertically resolved canopy
    d_pat           ! horizontal pattern length scale

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,kcm:ke1), OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,kcm:), OPTIONAL, INTENT(INOUT) :: &
#endif
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml           ! effective drag coefficient of canopy elements
                    ! smaller than the turbulent length scale            (1/m)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,kcm-1:ke1), OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,kcm-1:), OPTIONAL, INTENT(INOUT) :: &
#endif
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

! Provisional values for pattern lenth array:
  IF (PRESENT(d_pat)) THEN
     DO i=i_stp, i_enp
        IF (fr_land(i) <= z1d2) THEN
           d_pat(i)=z0
        ELSE
           d_pat(i)=pat_len !should be a 2D external parameter field
        END IF
     END DO
  END IF
! Effective values of the surface area indices:
  IF (PRESENT(sai) .AND. PRESENT(eai) .AND. PRESENT(tai) .AND. &
      PRESENT(lai) .AND. PRESENT(plcov)) THEN

     DO i=i_stp, i_enp
        IF (fr_land(i) <= z1d2) THEN
           sai(i)=c_sea
        ELSE
           tai(i)=MAX( 1.0E-6_ireals, lai(i) )
        END IF
     END DO

     IF (itype_tran.EQ.1) THEN
        DO i=i_stp, i_enp
           IF (fr_land(i) > z1d2) THEN
              sai(i)=tai(i)
              eai(i)=(z1-plcov(i))*sai(i)
              tai(i)=plcov(i)*tai(i)
            END IF
        END DO
     ELSE
        DO i=i_stp, i_enp
           IF (fr_land(i) > z1d2) THEN
              tai(i)=plcov(i)*tai(i)  ! transpiration area index
              eai(i)=c_soil               ! evaporation area index
              sai(i)=c_lnd+tai(i)       ! surface area index

              fakt=EXP( e_surf*LOG( sai(i)) )/sai(i)
            ! effective area indeces by multiplication with the reduction factor fakt:
              sai(i)=fakt*sai(i)
              eai(i)=fakt*eai(i)
              tai(i)=fakt*tai(i)
           END IF
        END DO
     END IF

  END IF

END SUBROUTINE init_canopy

!********************************************************************************
!********************************************************************************

!+ Module procedure organize_turbdiff in "src_turbdiff" for organising the calls 
!+ of surface-to-atmosphere trasnsfer and turbulent diffusion:

SUBROUTINE organize_turbdiff (lstfnct, lsfluse, &
!
          lturatm,ltursrf, iini, &
          lmomdif,lscadif, itnd, &
!
          dt_var,dt_tke, nprv,ntur,ntim, &
!
          ie, ke, ke1, kcm, &
!
          i_st, i_en, i_stp, i_enp, &
!
          l_hori, &
#ifdef __COSMO__
          eddlon, eddlat, edadlat, acrlat, &
#endif
!---------------------------------------------------------------------
          hhl, dp0, &
!
          fr_land, depth_lk, h_ice, gz0, sai, &
!
          d_pat, c_big, c_sml, r_air, &
!
          t_g, qv_s, ps, &
          u, v, w, t, qv, qc, prs, rho, epr, &
!
          ptr, &
!
          tcm, tch, tfm, tfh, tfv, dzm, dzh, &
          tke, tkvm, tkvh, rcld, &
          edr, tket_sso, tket_conv, &
          u_tens, v_tens, t_tens, &
          qv_tens, qc_tens, &
          tketens, &
          qv_conv, ut_sso, vt_sso, &
!
          t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m, &
          shfl_s, lhfl_s, &
!
          ierrstat, errormsg, eroutine)

 
!-------------------------------------------------------------------------------
! Description:
!
! Organizes the CALL of 'turbtran' and 'turbdiff'
! and 'iinit'.

! Method:
!
! All tendency parameters except 'tketens' are OPTIONAL. If they are missing calculated 
! tendencies of SUB 'turbdiff' are automatically added to the related prognostic variables.
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
!
   lstfnct,      & !calculation of stability function required
   lsfluse,      & !use explicit heat flux densities at the suface
   lmomdif,      & !calculation of complete gradient diffusion of horizontal momenum
   lscadif,      & !calculation of complete gradient diffusion of scalar properties
   lturatm,      & !running turbulence model between atmosph. layers (updating diffusion coefficients)
   ltursrf         !running turbulence model at the surface layer (updating transfer coefficients 
                   !                                               and near surface variables) 

REAL (KIND=ireals), INTENT(IN) :: & 
!
   dt_var,       & !time step for ordinary prognostic variables
   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
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
    i_st,    i_en,    & ! start and end index of mass points 
    i_stp,   i_enp      ! start and end index of mass points including model boundary lines

! Time step indices:
! -----------------------------------------------------------------------

REAL (KIND=ireals), INTENT(IN) :: &
!
! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------
!
    l_hori          ! horizontal grid spacing (m)

!---------------------------------------------------------------------

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
#endif
!
    hhl             ! height of model half levels                   ( m )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
    dp0             ! pressure thickness of layer                   (pa )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
#endif
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area             ( 1 )
    depth_lk,     & ! lake depth                                    ( m )
    sai             ! surface area index                            ( 1 )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
!
    d_pat           ! horizontal pattern length scale

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,kcm:ke1), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,kcm:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml           ! effective drag coefficient of canopy elements
                    ! smaller than the turbulent length scale            (1/m)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,kcm-1:ke1), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,kcm-1:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)

! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
#endif
!
    h_ice,        & ! ice thickness                                 (  m  )
    ps              ! surface pressure                              ( pa  )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), TARGET, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), TARGET, INTENT(IN) :: &
#endif
!
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), TARGET, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
#endif
!
! Atmospheric model variables:
! ---------------------------------
!
     u,           & ! zonal wind speed                              ( m/s )
     v,           & ! meridional wind speed                         ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), TARGET, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(IN) :: &
#endif
!
     prs            ! atmospheric pressure                          ( pa  )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
!
     rho,         & ! total density of air                          (kg/m3)
     epr            ! exner pressure                                 (1)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
     w              ! vertical wind speed (defined on half levels)  ( m/s )

TYPE (modvar), OPTIONAL :: ptr(:) !passive tracers

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), TARGET, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), TARGET, INTENT(INOUT) :: &
#endif
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------
!
     gz0,          & ! roughness length * g of the vertically not
                     ! resolved canopy                               (m2/s2)
     tcm,          & ! turbulent transfer coefficients for momentum    --
     tch,          & ! turbulent transfer coefficients for heat        --
!
     tfm,          & ! factor of pure surface layer transfer of momentum          --
     tfh,          & ! factor of pure surface layer transfer of scalars           --
     tfv             ! laminar reduction factor for evaporation        --
 
#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), TARGET, OPTIONAL,  INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL,  INTENT(OUT) :: &
#endif
!
     dzm,          & ! effctive surface layer depth form momentum    ( m  )
     dzh             ! effctive surface layer depth form heat        ( m  )

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=ireals), DIMENSION(ie,ke1,ntim), INTENT(INOUT) :: &
!
     tke             ! q:=SQRT(2*TKE); TKE='turbul. kin. energy'     ( m/s )
                     ! (defined on half levels)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,2:ke1), TARGET, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,2:), TARGET, INTENT(INOUT) :: &
#endif
!
     tkvm,         & ! turbulent diffusion coefficient for momentum  (m/s2 )
     tkvh            ! turbulent diffusion coefficient for heat      (m/s2 )
                     ! (and other scalars)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), TARGET, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
#endif
!
     rcld            ! standard deviation of the saturation deficit        
                     ! (as input and output)                             
                     ! fractional cloud cover (in turbdiff)            --

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, TARGET, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(INOUT) :: &
#endif
!
! Tendency fields for the prognostic variables:
! -----------------------------------------------
!
     u_tens,       & ! u-tendency                                    ( m/s2)
     v_tens,       & ! v-tendency                                    ( m/s2)
     t_tens,       & ! t-tendency                                    ( K/s )
     qv_tens,      & ! qv-tendency                                   ( 1/s )
     qc_tens         ! qc-tendency                                   ( 1/s )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), TARGET, OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
#endif
!
     tketens         ! tendency of q=SQRT(2*TKE)                     ( m/s2)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: &
#endif
!
     qv_conv         ! qv-flux-convergence                            ( 1/s )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
     ut_sso,       & ! u-tendency due to the SSO-Scheme              ( 1/s )
     vt_sso          ! v-tendency due to the SSO-Scheme              ( 1/s )

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), OPTIONAL, TARGET, INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(OUT) :: &
#endif
!
     edr             ! eddy dissipation rate of TKE (EDR)            (m2/s3)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
     tket_conv       ! TKE-tendency due to convective buoyancy       (m2/s3)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: &
#endif
!
     tket_sso        ! TKE-tendency due to SSO wake production       (m2/s3)

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(OUT) :: &
#endif
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

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie), OPTIONAL, TARGET, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), OPTIONAL, TARGET, INTENT(INOUT) :: &
#endif
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2) (positive downward)
     lhfl_s          ! latent   heat flux at the surface             (W/m2) (positive downward)

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

REAL (KIND=ireals), DIMENSION(ie,0:7), TARGET :: &
     dd        !local derived turbulence parameter

INTEGER (KIND=iintegers) :: &
     nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes
     it_start,   & !Startindex der Iterationen
     kem,        & !ke oder ke1
     ntrac,      & !number of passive tracers
     ndiff,      & !number of 1-st order model variables to be diffused
!
! Transformed horizontal indices:
!
! Zonal direction:
     istart,    iend,    & ! start and end index for the inner model domain
     istartpar, iendpar    ! start and end index including model boundary lines

INTEGER (KIND=iintegers) :: &
     n_dec, n_inc,              & !decrease of lowest and increase of largest horiz. index
     iid(nvel,2)                  !horizontal start- and end indices

LOGICAL ::  &
!
     lini,      & !initialization required
     ldovardif, & !berechne (teil-)implizite Vert.diff von Mod.var 1-ter Ordnung
     ldogrdcor, & !mache Gradientkorrektur bei Berechnung der vertikalen Diffusion
     ldomomcor    !mache Gradientkorrektur bei vertikaler Impulsdiffusion 
!
!     Zwischenspeicher fuer rcpv, rcpl und lhocp 
!     und Faktor zur Beruecksichtigung von Fluessigwasser:
!

REAL (KIND=ireals) :: &
!
     zrcpv,zrcpl, & !current value of 'rcpv' and 'rcpl' dependent on 'lcpfluc'
     zlhocp, liqfak !'lhocp' and liquid water factor dependent on mode of cloud representation

REAL (KIND=ireals), TARGET :: &
!
     c_tke,tet_g,c_g,rim, &
     d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
     a_3,a_5,a_6,b_1,b_2, &
     l_scal,fc_min

REAL (KIND=ireals), POINTER :: &
!
!    Pointer fuer Felder der Dichte, des Exnerfaktors und der EDR:
     rhoh(:,:),  exner(:,:), ediss(:,:)

REAL (KIND=ireals), TARGET :: &
!
!    Seicher-Target fuer obige Pointer:
     rho_tar(ie,ke), exner_tar(ie,ke)    

!Declaration of statement functions:

REAL (KIND=ireals) :: &
     zexner, zpres               !Exner factor and its argument

!Exner-factor:
 zexner(zpres)=(zpres/p0ref)**rdocp
 
!-------------------------------------------------------------------------------

 istat=0; ilocstat=0
 errormsg=''; eroutine='organize_turbdiff'; lerror=.FALSE.

!test: always 5 initial iterations
!it_end=5
!test

!Note: 
!Since TKE may not be present at the boundary of the model domain, the turbulence model
!can be applied to the total model domain including the boundary, if 'i_stp', 'i_enp'
!do contain these over all boundary as well.

!Except 'u(_tens)', 'v(_tens)' all variables are defined at horizontal mass positions.
!In the turbulence model all tendencies are calculated at the mass positions.

 ldogrdcor=(lexpcor   .AND. lturatm) !gradient correction has to be done
 ldovardif=(lmomdif   .OR.  lscadif) !any variable has to be diffused
 ldomomcor=(ldogrdcor .AND. lnonloc) !gradient correction for momentum 

 IF (PRESENT(ptr)) THEN !passive tracers are present
    ntrac=UBOUND(ptr,1) !number of tracers
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

 tet_g=grav/cp_d !adiabatic T-gradient

 kem=ke

 istartpar=i_stp; iendpar=i_enp
 
 iid(u_m,1)=i_st; iid(u_m,2)=i_en
 iid(v_m,1)=i_st; iid(v_m,2)=i_en

 IF (iini.GT.0) THEN !an initialization run
    lini=.TRUE.
    IF (iini.EQ.1) THEN !separate initialization before the time loop
       it_start=1 !only 'it_end' iterations for initialization
                  !and an additional one at the first time loop
    ELSE !initialization within the first time step
       it_start=0 !"it_end+1" iterations for initializatzion
    END IF
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

!CALL get_param

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
! SHistory:
! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 
!  Initial release
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

!     Belegung abgeleiteter Konstanten:

      c_tke=d_m**z1d3
      c_g=z1-z1/(a_m*c_tke)-z6*a_m/d_m !=3*0.08

      d_0=d_m

      d_1=z1/a_h
      d_2=z1/a_m
      d_3=z9*a_h
      d_4=z6*a_m
      d_5=z3*(d_h+d_4)
      d_6=d_3+z3*d_4

      rim=z1/(z1+(d_m-d_4)/d_5) !1-Rf_c

      a_3=d_3/(d_2*d_m)
      a_5=d_5/(d_1*d_m)
      a_6=d_6/(d_2*d_m)

      b_1=(z1-d_4/d_m)/d_1
      b_2=(z1-d_4/d_m-c_g)/d_2

      l_scal=MIN( l_hori, tur_len )
!__________________________________________________________________________
!test: frm ohne fc_min-Beschraenkung: Bewirkt Unterschiede!
! JF:       fc_min=(vel_min/MAX( l_hori, tur_len ))**2
      fc_min=z0
!__________________________________________________________________________

      DO i=istartpar,iendpar
         dd(i,0)=d_m

         dd(i,1)=d_1
         dd(i,2)=d_2
         dd(i,3)=d_3
         dd(i,4)=d_4
         dd(i,5)=d_5
         dd(i,6)=d_6

         dd(i,7)=rim
      END DO

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
! History:
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
!  Using the redefined parameter c_g (including the factor 3).
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
!  Optimizations for vectorization: splitted the big loop in Section 4
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
!              2010/09/30 Matthias Raschendorfer
!  Substitution of 'itype_diag_t2m' by 'itype_synd' being already used for that purpose
!   "1": alternative SYNOP-digansostics related to previous Lewis scheme (not included here)
!   "2": SYNOP-diagnostics according to current transfer scheme using SYNOP-z0 'z0d'
!   "3": like "2" but using 'z0d' only for 10m-wind and a specific roughness layer profile for
!        2m-temperature and -humidity.
!  Including the adiabatic lapse rate correction for 't_2m' for "itype_synd.EQ.3" as well.
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
!     Platzh. fuer q=SQRT(2*TKE)-Werte:
!
           q0, &
!
!     Platzh. fuer Hoehendifferenzen und Laengaenskalen:
!   
           dh,l_turb,lh,lm,z0d,z_surf,len1,len2, &
           dz_sg_m, dz_sg_h, dz_g0_m, dz_g0_h, &
           dz_s0_m, dz_sa_m, &
           h_2m, h_10m, a_top, a_atm, a_2m, a_10m, &
!test
q1,q3,&
!
!     sonstiges:
!
           rat_m,rat_h,fac_m,fac_h, &
           fr_sd_h, &
           xf,wf

! Local arrays:

      LOGICAL            ::  &
!
        lo_ice(ie)   ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=ireals), TARGET ::  &
!
!US for a better vectorization (hopefully):
!
        h_top_2d (ie),    &
        h_atm_2d (ie),    &
!
        z0m_2d   (ie),    &
        z0d_2d   (ie),    &
        z2m_2d   (ie),    &
        z10m_2d  (ie),    &
!
        hk_2d    (ie),    &
        hk1_2d   (ie),    &
        h_can_2d (ie),    &
!
        rat_m_2d (ie),    &
        rat_h_2d (ie),    &
        fac_h_2d (ie),    &
        fm2_2d   (ie),    &
        fh2_2d   (ie),    &
        frc_2d   (ie),    &

!
        vel_2d   (ie),    &
        prs_2d   (ie),    &
!
        dz_0a_m  (ie),    &
        dz_0a_h  (ie),    &
        dz_sa_h  (ie),    &
        dz_s0_h  (ie)

     REAL (KIND=ireals), POINTER :: &
!
        vel1_2d  (:),      &
        vel2_2d  (:),      &
        ta_2d    (:),      &
        qda_2d   (:),      &
        ts_2d    (:),      &
        qds_2d   (:),      &
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
           a(ie,ke1:ke1, ndim), & !special thermodynamic variables
        eprs(ie,ke1:ke1),       & !exner pressure
        rclc(ie,ke1:ke1),       & !cloud cover
   len_scale(ie,ke1:ke1),       & !turbulent master length scale
 tketens_tar(ie,ke1:ke1)          !target for total transport of SQRT(TKE)

     INTEGER (KIND=iintegers) ::  &
        k_2d     (ie)

     REAL (KIND=ireals) :: &
!
!     Vertikale Gradienten verschiedener thermodyn. Variablen:
!
        grad(ie,nred) 

     REAL (KIND=ireals) :: &
        rho_2d(ie,ke1:ke1)

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
      ts_2d   => vari(:,ke1,tet_l)
      qds_2d  => vari(:,ke1,h2o_g)

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

      DO i=istartpar, iendpar

         IF (fr_land(i) <= z1d2) THEN
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

      IF (lini) THEN

         DO i=istartpar, iendpar

            IF (fr_land(i) <= z1d2) THEN

!              Ueber Wasserpunkten:

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice

               IF ( lo_ice(i) ) THEN

!                 Bei Eisdecke: 

                  gz0(i)=grav*z0_ice
               ELSE

!                 Bei von Schubspannung abhaengiger Wellenhoehe:

!                 Einfachste Schaetzung der Schubspannung als Impusls-
!                 flussdichte durch die Nebenflaeche ke mit Hilfe
!                 einer diagnostischen TKE ohne Beruecksichtigung von
!                 Feuchte-Effekten und mit neuchtralen Stabilitaets-
!                 funktion und Anwendung der Charnockflormel:

                  l_turb=hhl(i,ke)-hhl(i,ke1)
!test: different turb. length scale
                  l_turb=akt*MAX( len_min, l_turb/(z1+l_turb/l_scal) )
!  l_turb=akt*MIN( l_scal, hhl(i,ke)-hhl(i,ke1) )
!test

                  dh=hhl(i,ke-1)-hhl(i,ke1)

                  vel1=u(i,ke-1)
                  vel2=u(i,ke  )
                  grad(i,u_m)=(vel1-vel2)/dh

                  vel1=v(i,ke-1)
                  vel2=v(i,ke  )
                  grad(i,v_m)=(vel1-vel2)/dh

                  grad(i,tet_l)=z2*(t(i,ke-1)-t(i,ke))/dh &
                                 +tet_g

                  fm2=MAX( grad(i,u_m)**2+grad(i,v_m)**2, fc_min )
                  fh2=grav*grad(i,tet_l)/t(i,ke)

                  ! Vereinfachte Loesung mit Rf=Ri:
                  IF (fh2.GE.(z1-rim)*fm2) THEN
                     ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                     ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                     fakt=z1/rim-z1
                     lm=l_turb*(b_2-(a_6+a_3)*fakt)
                     lh=lm
                  ELSE
                     fakt=fh2/(fm2-fh2)
                     lm=l_turb*(b_2-(a_6+a_3)*fakt)
                     lh=l_turb*(b_1-a_5*fakt)
                  END IF

                  val1=lm*fm2; val2=lh*fh2
                  wert=MAX( val1-val2, rim*val1 )

                  q0=SQRT(d_m*l_turb*wert)

                  tke(i,ke,nvor)=q0

                  wert=MAX(ustmin,lm*q0*SQRT(fm2))
!                 gz0(i)=MAX( grav*len_min, alpha0*wert )
                  gz0(i)=MAX( grav*len_min, &
                                alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF

            tkvm(i,ke)=con_m
            tkvh(i,ke)=con_h
            tkvm(i,ke1)=con_m
            tkvh(i,ke1)=con_h

            rcld(i,ke1)=z0

         END DO
      END IF

! 4)  Berechnung der Transferkoeffizienten:

      DO it_durch=it_start, it_end !Iterationen
!print *,"it_durch=",it_durch

         DO i=istartpar, iendpar

!           Berechnung der benoetigten Laengenskalen:

            ! Dicke der Modell-Prandtl-Schicht
            h_top_2d(i) = hhl(i,ke)-hhl(i,ke1)
            h_atm_2d(i) = h_top_2d(i)*xf        

            ! Hoehe des 2m- und 10m-Niveaus
            h_2m  = z2
            h_10m = z10

            ! Rauhigkeitslaenge
            z0m_2d(i) = gz0(i)/grav
            z_surf      = z0m_2d(i)/sai(i)

            ! turbulente Laengenskala 
            l_tur_z0(i) = akt*z0m_2d(i) 

            ! turbulente Distanz auf der untersten Hauptflaeche
            a_atm = h_atm_2d(i)+z0m_2d(i) 

            ! turbulente Distanz auf der untersten Nebenflaeche
            a_top = h_top_2d(i)+z0m_2d(i) 

!           Laminare Korrektur der Diffusionskoeffizienten:
            tkvm(i,ke1)=MAX( con_m, tkvm(i,ke1) )
            tkvh(i,ke1)=MAX( con_h, tkvh(i,ke1) )

            fakt=z1+(z1-REAL(NINT(fr_land(i)),ireals))*(rat_sea-z1)

            rat_m=tkvm(i,ke1)/con_m
            rat_h=tkvh(i,ke1)/con_h

!           Berechnung der effektiven Widerstandslaenge der L-Schicht:
            dz_sg_m=rlam_mom*z_surf
            dz_sg_h=fakt*rlam_heat*z_surf*(rat_h/rat_m)

!           Berechnung weiterer effektiver Widerstandslaengen fuer Skalare:

!           Bestandesschicht ohne lam. Grenzschicht:
            dz_g0_h=z_surf*LOG(rat_m)

!           Bestandesschicht inclusive lam. Grenzschicht:
            dz_s0_h(i)=dz_sg_h+dz_g0_h

!           Berechnung der effektiven Bestandeshoehe:
            IF (dz_sg_h.eq.z0) THEN
               h_can_2d(i)=rat_can*z0m_2d(i)
            ELSE
               h_can_2d(i)=rat_can*dz_s0_h(i)*LOG(dz_s0_h(i)/dz_sg_h)
            END IF

!           Berechnung weiterer effektiver Widerstandslaengen fuer Impuls:

!           Widerstandslaengen inclusive lam. Grenzschicht:
            wert=z1d2*dz_sg_m 
            dz_s0_m=wert+SQRT(wert**2+h_can_2d(i)*dz_sg_m)

!           Widerstandslaengen ohne lam. Grenzschicht:
            dz_g0_m=dz_s0_m-dz_sg_m

!           der turb. Prandtl-Schicht:

            rat_m=(tkvm(i,ke)*z0m_2d(i))/(tkvm(i,ke1)*a_top)
            rat_m_2d(i)=MIN( z2, MAX( z1d2, rat_m ) )
            
            fac_m=(rat_m_2d(i)-z1)*z0m_2d(i)/h_top_2d(i)
            IF (fac_m.EQ.z1) THEN
               dz_0a_m(i)=z0m_2d(i)*h_atm_2d(i)/a_atm
            ELSE
               dz_0a_m(i)=z0m_2d(i)*LOG(a_atm/(z0m_2d(i)+fac_m*h_atm_2d(i)))/(z1-fac_m)
            END IF

            rat_h=(tkvh(i,ke)*z0m_2d(i))/(tkvh(i,ke1)*a_top)
            rat_h_2d(i)=MIN( z2, MAX( z1d2, rat_h ) )

            fac_h_2d(i)=(rat_h_2d(i)-z1)*z0m_2d(i)/h_top_2d(i)
            IF (fac_h_2d(i).EQ.z1) THEN
               dz_0a_h(i)=z0m_2d(i)*h_atm_2d(i)/a_atm
            ELSE
               dz_0a_h(i)=z0m_2d(i)*LOG(a_atm/(z0m_2d(i)+fac_h_2d(i)*h_atm_2d(i))) &
                                      /(z1-fac_h_2d(i))
            END IF

!           von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):
            dz_sa_m      = dz_s0_m      + dz_0a_m(i)
            dz_sa_h(i) = dz_s0_h(i) + dz_0a_h(i)

!           Reduktionsfaktoren fuer die Bestandesschicht
!           incl. lam. Grenzschicht:

            tfm(i)=dz_0a_m(i)/dz_sa_m      
            tfh(i)=dz_0a_h(i)/dz_sa_h(i)
!if (i.eq.im) then
!  print *,"dz_0a_h,sa=",dz_0a_h(i), dz_sa_h(i)," tfh=",tfh(i)
!endif


!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den
!           Faktor 'rat_lam' gegenueber fuehlbarer Waerme vergroesserten
!           laminaren Transpostwiderstandes:

            tfv(i)=z1/(z1+(rat_lam-z1)*dz_sg_h/dz_sa_h(i))
         END DO

         IF (icldm_tran.EQ.-1) THEN
            !Keine Wolkenberuecksichtigung;
            !Wolkenwasser wird ignoriert (vollkommen trockenes Schema):
   
            zlhocp=z0 !no condensation
            liqfak=z0
         ELSE
            zlhocp=lhocp
            liqfak=z1
         END IF

!        Berechnung der Erhaltungsgroessen in der Prandtl-Schicht:

         IF (lprfcor) THEN
            ks=ke-1
         ELSE
            ks=ke
         END IF

         DO k=ks, ke
            IF (.NOT.PRESENT(epr)) THEN
               DO i=istartpar, iendpar
                  exner(i,k)=zexner(prs(i,k))
               END DO
            END IF   
            DO i=istartpar, iendpar
               vari(i,k,u_m)  = u(i,k)
               vari(i,k,v_m)  = v(i,k)
               vari(i,k,tet_l)=(t(i,k) - zlhocp*qc(i,k))/exner(i,k)
!modif: liqfak
               vari(i,k,h2o_g)=qv(i,k) + liqfak*qc(i,k)
!modif
            END DO
         END DO    

         IF (lprfcor) THEN
            DO i=istartpar, iendpar
               len1=z2*h_top_2d(i)
               len2=(h_top_2d(i)-h_atm_2d(i))**2 &
                   /((hhl(i,ke-1)+hhl(i,ke))*z1d2-hhl(i,ke1)-h_atm_2d(i))
               lm=len1-tfm(i)*h_atm_2d(i)-len2
               lh=len1-tfh(i)*h_atm_2d(i)-len2

               epr_2d(i)=zexner(ps(i))

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
         DO i=istartpar, iendpar
            vel_2d(i)=MAX( vel_min, SQRT(vel1_2d(i)**2+vel2_2d(i)**2) )
         END DO

!        Vorbelegung mit Erhaltungsvariablen auf unterster Hauptflaeche als Referenz-Niveau:

         DO i=istartpar, iendpar
            vari(i,ke1,tet_l)= vari(i,ke,tet_l)
            vari(i,ke1,h2o_g)= vari(i,ke,h2o_g)
         END DO

!        Thermodynamische Hilfsvariablen auf dem Unterrand der Prandtl-Schicht:

         CALL adjust_satur_equil ( ke1,                                    &
!
              istartpar,iendpar, ke1,ke1,                                  &
!  
              lcalrho=.TRUE., lcalepr=.NOT.lprfcor, lcaltdv=.TRUE.,        &
              lpotinp=.FALSE., ladjout=.FALSE.,                            &
!
              icldmod=icldm_tran,                                          &
!
              zrcpv=zrcpv, zrcpl=zrcpl, zlhocp=zlhocp,                     &
!
              prs=RESHAPE(ps  ,(/ie,1/)), t=RESHAPE(t_g,(/ie,1/)),         &
                                         qv=RESHAPE(qv_s,(/ie,1/)),        &
!
              fip=tfh, &
!
              exner=eprs(:,ke1:ke1), rcld=rcld(:,ke1:ke1),                 &
                                     dens=rho_2d(:,ke1:ke1),               &
!
              qst_t=a(:,ke1:ke1,3), g_tet=a(:,ke1:ke1,4),                  &
                                      g_h2o=a(:,ke1:ke1,5),                &
!
              tet_l=vari(:,ke1:ke1,tet_l), q_h2o=vari(:,ke1:ke1,h2o_g),    &
                                             q_liq=vari(:,ke1:ke1,liq)       )

!        Beachte: 'ts_2d' und 'qds_2d' zeigen auf die 'ke1'-Schicht von 'vari(tet_l)' und 'vari(h2o_g)
!                 und sind jetzt die Erhaltungsvariablen am Unterrand der Prandtl-Schicht.
!                 Genauso zeigen auch ta_2d und qda_2d auf die entsprechenden Erhaltungsvariablen.

!        Berechnung der benoetigten Vertikalgradienten und der TKE-Antriebe:

         DO n=1, nvel  
            DO i=istartpar, iendpar
               grad(i,n)=vari(i,ke,n)*tfm(i)/dz_0a_m(i)
            END DO
            !Beachte: Dies ist die Darstellung ohne Nutzung der unteren Randwerte der Prandtl-Schicht
         END DO   
         DO n=nvel+1,nred
            DO i=istartpar, iendpar
               grad(i,n)=(vari(i,ke,n)-vari(i,ke1,n))/dz_0a_h(i)
            END DO
         END DO   

         DO i=istartpar, iendpar
            fm2_2d(i)=MAX( grad(i,u_m)**2+grad(i,v_m)**2, fc_min )
            fh2_2d(i)=g_tet(i)*grad(i,tet_l)+g_vap(i)*grad(i,h2o_g)
         END DO

!        Berechnung der Stabilitaetslaengen:

         IF (it_durch.EQ.it_start .AND. lini) THEN !Startinitialisierung

            DO i=istartpar, iendpar
               IF (fh2_2d(i).GE.(z1-rim)*fm2_2d(i)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  tkvm(i,ke1)=l_tur_z0(i)*(b_2-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=tkvm(i,ke1)
               ELSE
                  fakt=fh2_2d(i)/(fm2_2d(i)-fh2_2d(i))
                  tkvm(i,ke1)=l_tur_z0(i)*(b_2-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=l_tur_z0(i)*(b_1-a_5*fakt)
               END IF

               val1=tkvm(i,ke1)*fm2_2d(i); val2=tkvh(i,ke1)*fh2_2d(i)
               wert=MAX( val1-val2, rim*val1 )
               tke(i,ke1,nvor)=MAX( SQRT(d_m*l_tur_z0(i)*wert), vel_min )

            END DO

         ELSE ! mit Hilfe der vorhergehenden TKE-Werte 

            DO i=istartpar, iendpar
               tkvm(i,ke1)=tkvm(i,ke1)/tke(i,ke1,nvor)
               tkvh(i,ke1)=tkvh(i,ke1)/tke(i,ke1,nvor)
            END DO
         END IF

! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes:

         CALL solve_turb_budgets (ke1,                                                              &
!
                                  istartpar,iendpar, ke1,ke1,                                       &
!
! JF:                             imode_stke=imode_tran,                                            &
                                  imode_stke=0,                                                     &
!
                                  fm2=RESHAPE(fm2_2d,(/ie,1/)), fh2=RESHAPE(fh2_2d,(/ie,1/)),       &
#ifdef SCLM
                                  it_s=it_durch, grd=RESHAPE(grad,(/ie,1,nred/)),                   &
#endif
!-------------------------------------------------------------------------------
                                  tls=len_scale(:,ke1:ke1), tvt=tvt(:,ke1:ke1) )

         DO i=istartpar, iendpar

! 4h)       Bestimmung der durch Wirkung der L-Schicht
!           korrigierten Diffusionskoeffizienten 
!           und der zugehoerigen Transferkoeffizienten:

!           unkorrigierte Diffusionskoeffizienten:

            tkvm(i,ke1)=tke(i,ke1,ntur)*tkvm(i,ke1)
            tkvh(i,ke1)=tke(i,ke1,ntur)*tkvh(i,ke1)

!           Belegung der Felder fuer die Transferkoeffizienten:

            tcm(i)=tkvm(i,ke1)*tfm(i)/(dz_0a_m(i)*vel_2d(i))
            tch(i)=tkvh(i,ke1)*tfh(i)/(dz_0a_h(i)*vel_2d(i))
!if (i.eq.im) then
!  print *,"in turbtran: dz_0a,sa=",dz_0a_h(i),dz_sa_h(i)
!  print *,"vel_2d=",vel_2d(i)
!  print *,"  tkvh=",tkvh(i,ke1)," tch=",tch(i)," tfh=",tfh(i)
!endif

! 4i)       Diagnose von gz0 (fuer den naechsten Zeitschritt) 
!           ueber Wasserflaechen mit der (angepassten) Charnock-Formel
!           und Einschraenkung von z0m_dia ueber Land:
 
            IF (fr_land(i) <= z1d2) THEN

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice
               IF ( lo_ice(i) ) THEN
                  ! Ice-covered grid box
                  gz0(i)=grav*z0_ice
               ELSE
                  velo=(tke(i,ke1,ntur)+tke(i,ke,nvor))*z1d2
                  wert=MAX(ustmin,tcm(i)*vel_2d(i)*SQRT(vel_2d(i)**2+velo**2))
!                 gz0(i)=MAX( grav*len_min, alpha0*wert )
                  gz0(i)=MAX( grav*len_min, alpha0*wert+grav*alpha1*con_m/SQRT(wert) )
               END IF
               !Ueber See gibt es keinen synoptischen Garten
               z0d_2d(i)=z0m_2d(i)
            ELSE
               !Die Rauhigkeitslaenge einer SYNOP Station soll immer
               !kleiner als 10m bleiben:
               z0d_2d(i)=MIN( z10, z0m_dia )
            END IF
         END DO

         IF (it_durch.LT.it_end) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF

      END DO !Iteration

      
! JF:       IF (iini.EQ.1) THEN !only for separate initialization before the time loop
! JF:          RETURN !finish this subroutine
! JF:       END IF

!     Berechnung der Enthalpieflussdichten:
      
      IF (PRESENT(shfl_s)) THEN 
         DO i=istartpar,iendpar
            shfl_s(i)=cp_d*rho_2d(i,ke1)*tkvh(i,ke1)*grad(i,tet_l)
         END DO
      END IF
      IF (PRESENT(lhfl_s)) THEN 
         DO i=istartpar,iendpar
            lhfl_s(i)=lh_v*rho_2d(i,ke1)*tkvh(i,ke1)*grad(i,h2o_g)
         END DO
      END IF   

! 4j) Berechnung der Standardabweichnung des Saettigungsdefizites:

!k->ke1
      DO i=istartpar, iendpar
         rcld(i,ke1)=SQRT(l_tur_z0(i)*tkvh(i,ke1)*d_h)* &
                       ABS(epr_2d(i)*qsat_dT(i)*vari(i,ke1,tet_l)-vari(i,ke1,h2o_g))
      ENDDO

! 4h) Berechnung der Enthalpieflussdichten und der EDR am Unterrand:

      IF (PRESENT(edr)) THEN
         DO i=i_st,i_en
            edr(i,ke1)=tke(i,ke1,ntur)**3/(d_m*l_tur_z0(i))
         END DO
      END IF

      IF (PRESENT(shfl_s)) THEN
         DO i=istartpar,iendpar
            shfl_s(i)=cp_d*rho_2d(i,ke1)*tkvh(i,ke1)*grad(i,tet_l)*epr_2d(i)
            !Note: shfl_s is positive downward and belogns to the T-equation!
         END DO
!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
         IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN 
            !measured SHF has to be used for forcing:
            shfl_s(im)=SHF%mod(0)%val
         ELSEIF (lsurflu) THEN !SHF defined by explicit surface flux density
            SHF%mod(0)%val=shfl_s(im)
            SHF%mod(0)%vst=MAX(i_upd, SHF%mod(0)%vst) !SHF is at least updated
         END IF
#endif
!SCLM --------------------------------------------------------------------------------
      END IF

      IF (PRESENT(lhfl_s)) THEN
         DO i=istartpar,iendpar
            lhfl_s(i)=lh_v*rho_2d(i,ke1)*tkvh(i,ke1)*grad(i,h2o_g)
            !Note: lhfl_s is positive downward!
         END DO
!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
         IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN 
            !measured LHF has to be used for forcing:
            lhfl_s(im)=LHF%mod(0)%val
         ELSEIF (lsurflu) THEN !LHF defined by explicit surface flux density
            LHF%mod(0)%val=lhfl_s(im)
            LHF%mod(0)%vst=MAX(i_upd, LHF%mod(0)%vst) !LHF is at least updated
         END IF
#endif
!SCLM --------------------------------------------------------------------------------
      END IF

!----------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      DO i=istartpar, iendpar

         IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile
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

      IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile

         val2=z1/epsi

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
                  ELSE
                     wert=z0m_2d(i)*LOG((z2m_2d(i)+z0m_2d(i))/            &
                         (z0m_2d(i)+fac_h_2d(i)*z2m_2d(i)))/(z1-fac_h_2d(i))
                  END IF
                  fakt=(dz_s0_h(i)+wert)/dz_sa_h(i)
               END IF
                ts_2d(i) = fakt*ta_2d(i) + (z1-fakt)*t_g(i)/epr_2d(i)
               qds_2d(i) = qv_s(i) + fakt*(qda_2d(i)-qv_s(i))

               prs_2d(i) = ps(i)

             ! t_2m(i) = t_g(i) + (ta_2d(i)-t_g(i))*fakt &
             !           + tet_g*( (h_atm_2d(i)+h_can_2d(i) )*fakt-h_2m )       !Achtung!
             ! qv_2m(i)= qv_s(i) + (qda_2d(i)-qv_s(i))*fakt

            END IF
         END DO   

      ELSE !using only a logarithmic profile above a SYNOP lawn
               
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

               fac_h=(rat_h_2d(i)-z1)*z0d/h_top_2d(i)
               IF (fac_h.EQ.z1) THEN
                  val1=fr_sd_h+h_2m/a_2m
                  val2=fr_sd_h+h_atm_2d(i)/a_atm
               ELSE
                  val1=fr_sd_h+LOG(a_2m/(z0d+fac_h*h_2m))/(z1-fac_h)
                  val2=fr_sd_h+LOG(a_atm/(z0d+fac_h*h_atm_2d(i)))/(z1-fac_h)
               END IF

               fakt=val1/val2

                ts_2d(i) = fakt*ta_2d(i) + (z1-fakt)*t_g(i)/epr_2d(i)
               qds_2d(i) = qv_s(i) + fakt*(qda_2d(i)-qv_s(i))

               prs_2d(i) = ps(i)
            END IF   
         END DO   

      END IF   

      DO i=istartpar, iendpar
         IF (k_2d(i).LT.ke) THEN
!           2m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!           trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommenen
!           linearen Profile der progn. Modellvariablen eine logarith. Interpolation:

            k2=k_2d(i); k1=k2+1
           
            wert=(h_2m    +z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
            fakt=(hk_2d(i)+z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
            fakt=LOG(wert)/LOG(fakt)
!test: different interpol. weight
!  fakt=(z2m_2d(i)-hk1_2d(i))/(hk_2d(i)-hk1_2d(i))
!test
           
            val2=qv(i,k2)+qc(i,k2)
            val1=qv(i,k1)+qc(i,k1)

            qds_2d(i)=val1+fakt*(val2-val1)
          
            IF (.NOT.PRESENT(epr)) THEN
               IF (k2.LT.ks) exner(i,k2)=zexner(prs(i,k2)) 
               IF (k1.LT.ks) exner(i,k1)=zexner(prs(i,k1)) 
            END IF
            val2=(t(i,k2)-zlhocp*qc(i,k2))/exner(i,k2)
            val1=(t(i,k1)-zlhocp*qc(i,k1))/exner(i,k1)
             ts_2d(i)=val1+fakt*(val2-val1)

            prs_2d(i)=prs(i,k1)

            rcl_2d(i)=rcld(i,k1)+fakt*(rcld(i,k2)-rcld(i,k1))
         ELSEIF (z2m_2d(i).LE.z0) THEN
            rcl_2d(i)=rcld(i,ke1)
         ELSE
            fakt=z2m_2d(i)/h_atm_2d(i)
            rcl_2d(i)=rcld(i,ke1)+fakt*(rcld(i,ke)-rcld(i,ke1))
         END IF    
      END DO    

!     Berechnung der zugehoerigen Modell- und Feuchtevariablen im 2m-Niveau:

      DO i=istartpar, iendpar
         wert=ts_2d(i)*(z1+rvd_m_o*qds_2d(i))     !angenaeherte virt. Temp. 
         prs_2d(i)=prs_2d(i)             &        !Druck
                    *EXP(-(z2m_2d(i)-hk1_2d(i)) &
                    *grav/(r_d*wert))
      END DO 
      !Beachte, dass sich die 2m-Werte bislang auf die Erhalturngsvariablen beziehen

      CALL adjust_satur_equil ( ke1,                                     &
!
           istartpar,iendpar, ke1,ke1,                                   &
!  
           lcalrho=.FALSE., lcalepr=.TRUE., lcaltdv=.FALSE.,             &
           lpotinp= .TRUE., ladjout=.TRUE.,                              &
!
           icldmod=icldm_tran,                                           &
!
           zrcpv=zrcpv, zrcpl=zrcpl, zlhocp=zlhocp,                      &
!
           prs=RESHAPE(prs_2d,(/ie,1/)), t=vari(:,ke1:ke1,tet_l),        &
                                        qv=vari(:,ke1:ke1,h2o_g),        &
!
           exner=eprs(:,ke1:ke1), rcld=rclc(:,ke1:ke1),                  &
!
           qst_t=a(:,ke1:ke1,3), g_tet=a(:,ke1:ke1,4),                   &
                                   g_h2o=a(:,ke1:ke1,5),                 &
!
           tet_l=vari(:,ke1:ke1,tet_l), q_h2o=vari(:,ke1:ke1,h2o_g),     &
                                          q_liq=vari(:,ke1:ke1,liq)      )

      DO i=istartpar, iendpar
          t_2m(i)=vari(i,ke1,tet_l)
         qv_2m(i)=vari(i,ke1,h2o_g)
         patm=prs_2d(i)*qv_2m(i) &
             /(rdv+(z1-rdv)*qv_2m(i))               !Wasserdampfdruck

         fakt=patm/zpsat_w( t_2m(i) )
         rh_2m(i)=100.0_ireals*MIN( fakt, z1 )      !relative Feuchte

         wert=LOG(patm/b1)
         td_2m(i)=MIN( (b2w*b3-b4w*wert) &
                   /(b2w-wert), t_2m(i) )           !Taupunktstemperatur
      ENDDO

!     Diagnose der 10m-Groessen:

      CALL diag_level(istartpar,iendpar, z10m_2d, k_2d, hk_2d, hk1_2d)

      DO i=istartpar, iendpar

         IF (k_2d(i).EQ.ke) THEN

!           10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!           in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

            z0d=z0d_2d(i)
            a_atm=h_atm_2d(i)+z0d
            a_10m=h_10m+z0d

            fac_m=(rat_m_2d(i)-z1)*z0d/h_top_2d(i)
            IF (fac_m.EQ.z1) THEN
               val1=h_10m/a_10m
               val2=h_atm_2d(i)/a_atm
            ELSE
               val1=LOG(a_10m/(z0d+fac_m*h_10m))
               val2=LOG(a_atm/(z0d+fac_m*h_atm_2d(i)))
            END IF

            fakt=val1/val2

            u_10m(i)=vel1_2d(i)*fakt
            v_10m(i)=vel2_2d(i)*fakt
              
         ELSE
!           10m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!           trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommen
!           lenearen Profile der progn. Modellvariablen eine logarithm. Interpolation:

            k2=k_2d(i); k1=k2+1

            wert=(h_10m     +z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
            fakt=(hk_2d(i)+z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
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

      !Note: '(u, v)_10m' always belong to mass points! 

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
!  Using the redefined parameter c_g (including the factor 3).
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
!     Platzh. fuer relat. Wolkenanteil u. rezipr. virtueller Faktor:
!
           clcv,virt

      REAL (KIND=ireals) :: &
!
           fr_var, fr_tke    !1/dt_var, 1/dt_tke

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
!     sonstiges:
!
           x1,x2,x3

      REAL (KIND=ireals) :: flukon33, flukon43, flukon53, &
                            flukon34, flukon44, flukon54, flukon55, &
                            flux_3,   flux_4,   flux_5

!---------------------------------------------------------------------
#ifdef __COSMO__
      REAL (KIND=ireals) :: zs11, zs22, zs33, zs12, zs21, zs13, zs23
#endif
!---------------------------------------------------------------------

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
           wind(ie,ke1,nvel),& !horizontale Windkomponenten im
                                  !Massenzentrum, bzw. deren Tendenzen
!
           len_scale(ie,ke1),& !turbulent length-scale
!
           rhon(ie,ke1),& !Luftdichte auf Nebenflaechen (einschl. surface)
!
           frh(ie,ke1),&  !thermal forcing (1/s2) or thermal acceleration (m/s2)
           frm(ie,ke1),&  !mechan. forcing (1/s2) or mechan. accelaration (m/s2)
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
!     Flusskonversionsmatrix, um von den turb. Flussdichten der 2
!     scalaren Erhaltungsgroessen (bzgl. feuchtadiab. vert. Verrueck.)
!     auf diejenigen der nscal scalaren Modellvariablen zu gelangen:
!
!          flukon(nvel+1:nmvar,nvel+1:nred), &
!      
!          wrid aus Effizienzgruenden erstetzt druch Skalare 'flukon33', 'flukon34, etc.!
!
!     Hilfsfelder fuer eine und zwei Variablenschichten:
!
           lay(ie), lays(ie,2), vars(ie,2), &
!
           z0m(ie),  &    !Rauhigkeitslaenge (fuer Impuls)
           vh0(ie),  &    !horizontal velocity scale at lowest model level
           l_pat(ie),&    !Laengenskala der therm. Inhomogenitaeten
                             !der Erdbodenoberflaeche
!
           src(ie),  &   ! arbitrary source term
!
!     Vertikale Gradienten:
!
           grad(nmvar), &
!
!     Time increment and inverse time increment of ordinary prognostic variables:
!
           tinc(ndiff), tinv(ndiff), &
!
           hig(2) ! obere und untere Referenzhoehe bei der Bildung
                  ! nicht-lokaler Gradienten

      REAL (KIND=ireals), POINTER :: &
!
!          Pointer fuer Tendenzfelder:
! 
           utens(:,:), vtens(:,:), &
           ttens(:,:), qvtens(:,:), qctens(:,:), &
!
!          Pointer fuer Felder der Bestandes-Parameter:
!
           dzsm(:), dzsh(:), &
!
           dpat(:),    &
           cbig(:,:), csml(:,:), rair(:,:)

      REAL (KIND=ireals), TARGET :: &
!
!          Seicher-Target fuer obige Pointer:
!
           dzsm_tar(ie), dzsh_tar(ie), &
!
           dpat_tar(ie),   &
           cbig_tar(ie,kcm:ke1), csml_tar(ie,kcm:ke1), rair_tar(ie,kcm-1:ke1)
 
!          Note:
!          The following buffers wouldn't be necessary, if the related pointers above
!          were allowed to be allocated at run time:

      REAL (KIND=ireals), DIMENSION(:,:), POINTER :: &
!
            cur_prof, upd_prof, &
            expl_mom, impl_mom, invs_mom, &
            diff_dep
   
      LOGICAL ::  &
!         
           can_fields,   &  !canopy fields are present
           ltend(ndiff), &  !calculation of tendencies required
           lsfli(ndiff)     !surface value input is a flux density instead of a concentration

      TYPE (modvar) :: dvar(ndiff)  !model variables to be diffused

      TYPE (turvar) :: vtyp(ntyp)   !variable types (momentum and scalars)

      TYPE (varprf) :: pvar(naux+1) !vertical variable profiles

!---- End of header ------------------------------------------------------------

!     nur in LM-Umgebung:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbtran'; lerror=.FALSE.

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

      fr_var=z1/dt_var; fr_tke=z1/dt_tke

      IF (PRESENT(dzm)) THEN
         dzsm => dzm
      ELSE
         dzsm => dzsm_tar
       ! ALLOCATE ( dzsm(ie),         STAT=ilocstat ); istat = istat + ilocstat
      END IF
      IF (PRESENT(dzh)) THEN
         dzsh => dzh
      ELSE
         dzsh => dzsh_tar
       ! ALLOCATE ( dzsh(ie),         STAT=ilocstat ); istat = istat + ilocstat
      END IF

      can_fields=.TRUE.
      IF (PRESENT(d_pat)) THEN
         dpat => d_pat
      ELSE
         dpat => dpat_tar
       ! ALLOCATE ( dpat(ie),         STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF
      IF (PRESENT(c_big)) THEN
         cbig => c_big
      ELSE
         cbig => cbig_tar 
       ! ALLOCATE ( cbig(ie,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF
      IF (PRESENT(c_sml)) THEN
         csml => c_sml
      ELSE
         csml => csml_tar 
       ! ALLOCATE ( csml(ie,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF
      IF (PRESENT(r_air)) THEN
         rair => r_air
      ELSE
         rair => rair_tar 
       ! ALLOCATE ( rair(ie,kcm-1:ke1), STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF

      IF (istat /= 0) THEN
         ierrstat = 1004
         errormsg= &
         'ERROR *** Allocation of space for meteofields failed ***'
         lerror=.TRUE.; RETURN
      ENDIF

      IF (.NOT.can_fields) THEN
!print *,"vor init_canopy kcm=",kcm
         CALL init_canopy(ie=ie, ke=ke, ke1=ke1, kcm=kcm, &
              i_stp=istartpar, i_enp=iendpar, &
              fr_land=fr_land, &
              d_pat=dpat, c_big=cbig, c_sml=csml, r_air=rair) 
!print *,"nach init_canopy kcm=",kcm
      END IF

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

!Achtung: dvar(liq)%sv sollte Null sein!!

      lsfli(:)=.FALSE. !surface values are concentration by default

      dvar(u_m)%av  => u  ; dvar(u_m)%at => utens  ; dvar(u_m)%sv => NULL()
      dvar(v_m)%av  => v  ; dvar(v_m)%at => vtens  ; dvar(v_m)%sv => NULL()
      dvar(tem)%av  => t  ; dvar(tem)%at => ttens  ; dvar(tem)%sv => t_g
      dvar(vap)%av  => qv ; dvar(vap)%at => qvtens ; dvar(vap)%sv => qv_s
      dvar(liq)%av  => qc ; dvar(liq)%at => qctens ; dvar(liq)%sv => NULL()

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
         !measured SHF has to be used for forcing:
         lsfli(tem)=.TRUE.
      END IF
      IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN
         !measured LHF has to be used for forcing:
         lsfli(vap)=.TRUE.
      END IF
      !Note: the measured SHF and LHF must already have been loaded to shfl_s and lhfl_s!
#endif
!SCLM --------------------------------------------------------------------------------
         
      IF (lsfluse) THEN !use explicit heat flux densities at the surface
         lsfli(tem)=.TRUE.; lsfli(vap)=.TRUE.
      END IF

      IF ((lsfli(tem) .AND. .NOT.PRESENT(shfl_s)) .OR. &
          (lsfli(vap) .AND. .NOT.PRESENT(lhfl_s))) THEN
         ierrstat = 1004; lerror=.TRUE.
         errormsg='ERROR *** forcing with not present surface heat flux densities  ***'
         RETURN
      ENDIF

      IF (lsfli(tem)) dvar(tem)%sv => shfl_s
      IF (lsfli(vap)) dvar(vap)%sv => lhfl_s

      IF (PRESENT(ptr)) THEN !passive tracers are present
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
               dvar(n)%sv => NULL(); lsfli(n)=.FALSE.
            END IF
         END DO
      END IF

      vtyp(mom)%tkv => tkvm ; vtyp(mom)%dzs => dzsm
      vtyp(sca)%tkv => tkvh ; vtyp(sca)%dzs => dzsh

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
             
!-----------------------------------------------------------------------

! 1)  Vorbereitungen:


!-----------------------
   IF (lturatm) THEN
!-----------------------

!Achtung bei T-Gleichung in cv-Form
!gesamte Thermodynamik ueberpruefen auch gegenueber satad

      lcircterm=(pat_len.GT.z0)

      ldotkedif=(c_diff .GT.z0)

!     Berechnung abgeleiteter Parameter:

      CALL turb_param

!     Bestimmung der initialen Werte fuer die laminaren Reduktions-
!     koeffizienten und die Standardabw. des Saettigungsdef.:

      IF (lini) THEN 
         DO k=1, ke1
            DO i=istartpar,iendpar    
               rcld(i,k)=z0
            END DO
         END DO
         IF (itype_tran.NE.2) THEN 
            DO i=istartpar,iendpar    
               tfh(i)=z1
               tfm(i)=z1
            END DO
         END IF   
      END IF

      IF (icldm_turb.EQ.-1) THEN
         !Keine Wolkenberuecksichtigung;
         !Wolkenwasser wird ignoriert (vollkommen trockenes Schema):

         zlhocp=z0 !no condensation
         liqfak=z0 !total water without liquid water
      ELSE
         zlhocp=lhocp
         liqfak=z1
      END IF

!     Thermodynamische Hilfsvariablen auf Hauptflaechen:

      CALL adjust_satur_equil (  1,                                                             &
!
           istartpar,iendpar, 1,ke,                                                             &
!  
           lcalrho=.NOT.PRESENT(rho), lcalepr=.NOT.PRESENT(epr), lcaltdv=.TRUE.,                &
           lpotinp=.FALSE., ladjout=.FALSE.,                                                    &
!
           icldmod=icldm_turb,                                                                  &
!
           zrcpv=zrcpv, zrcpl=zrcpl, zlhocp=zlhocp,                                             &
!
           prs=prs, t=t, qv=qv, qc=qc,                                                          &
           exner=exner, rcld=rcld, dens=rhoh,                                                   &
!
           r_cpd=a(:,:,2), qst_t=a(:,:,3), g_tet=a(:,:,4), g_h2o=a(:,:,5),                      &
           tet_l=vari(:,:,tet_l), q_h2o=vari(:,:,h2o_g), q_liq=vari(:,:,liq) )

!     Vorbelegung mit Erhaltungsvariablen auf unterster Hauptflaeche als Referenz-Niveau:

      DO i=istartpar,iendpar
         vari(i,ke1,tet_l)= vari(i,ke,tet_l)
         vari(i,ke1,h2o_g)= vari(i,ke,h2o_g)
      END DO

!     Thermodynamische Hilfsvariablen auf Unterrand der Prandtl-Schicht:

      CALL adjust_satur_equil ( ke1,                                               &
!
           istartpar,iendpar, ke1,ke1,                                             &
!  
           lcalrho=.TRUE. , lcalepr=.TRUE., lcaltdv=.TRUE.,                        &
           lpotinp=.FALSE., ladjout=.FALSE.,                                       &
!
           icldmod=icldm_turb,                                                     &

           zrcpv=zrcpv, zrcpl=zrcpl, zlhocp=zlhocp,                                &
!
           prs=RESHAPE(ps,(/ie,1/)), t=RESHAPE(t_g, (/ie,1/)),                     &
                                    qv=RESHAPE(qv_s,(/ie,1/)),                     &
!
           fip=tfh,                                                                &
!
           exner=a(:,ke1:ke1,1), rcld=rcld(:,ke1:ke1), dens=rhon(:,ke1:ke1),       &
!
           r_cpd=a(:,ke1:ke1,2), qst_t=a(:,ke1:ke1,3), g_tet=a(:,ke1:ke1,4),       &
                                                           g_h2o=a(:,ke1:ke1,5),   &
           tet_l=vari(:,ke1:ke1,tet_l), q_h2o=vari(:,ke1:ke1,h2o_g),               & 
                                          q_liq=vari(:,ke1:ke1,liq)                )

!     Berechnung der horizontalen Windgeschwindigkeiten
!     im Massenzentrum der Gitterbox:

      DO k=1,ke
         DO i=istartpar,iendpar
            vari(i,k,u_m)=u(i,k)
            vari(i,k,v_m)=v(i,k)
         END DO   
      END DO   
      DO i=istartpar,iendpar
         vari(i,ke1,u_m)=vari(i,ke,u_m)*(z1-tfm(i))
         vari(i,ke1,v_m)=vari(i,ke,v_m)*(z1-tfm(i))
      END DO

!     Sichern der auf Massenpunkte interpolierten Windkomponenten:

      DO n=1,nvel
         DO k=1,ke1
            DO i=istartpar,iendpar
               wind(i,k,n)=vari(i,k,n)
            END DO
         END DO
      END DO

!     Berechnung der Schichtdicken und der Dichte auf Nebenflaechen:
        
      DO k=1,ke
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
                               nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dp0)
!___________________________________________________________________________

      pvar(1)%bl => rcld     ; pvar(1)%ml => rcld
      pvar(2)%bl => a(:,:,1) ; pvar(2)%ml => exner
      DO n=2,naux 
         pvar(1+n)%bl => a(:,:,n) ; pvar(1+n)%ml => a(:,:,n)
      END DO
      CALL bound_level_interp( istartpar,iendpar, 2,ke, &
                               nvars=naux+1, pvar=pvar, &
                               depth=dp0, rpdep=hlp)

!print *,"nach interpol"

!     Bestimmung der Rauhigkeitslaenge 
!     und der effektiven Dicke der Modell-Prandtlschicht:

      DO i=istartpar,iendpar    
         z0m(i)=gz0(i)/grav
         vh0(i)=MAX( vel_min, &
                     SQRT(wind(i,ke,u_m)**2+wind(i,ke,v_m)**2) )
      END DO

!bug_mod_2011/09/23: dicke(i,ke1) -> dz0(i,mom), dz0(i,sca) {
      IF (itype_tran.EQ.2) THEN !Transferkoeffizienten mit turbtran
         ! spezifiche effektive Dicke der Prandtlschicht:
         DO i=istartpar,iendpar    
            dzsm(i)=tkvm(i,ke1)*tfm(i)/(vh0(i)*tcm(i))
            dzsh(i)=tkvh(i,ke1)*tfh(i)/(vh0(i)*tch(i))
         END DO 
!print *,"in turbdiff: tkvh=",tkvh(im,ke1)," tch=",tch(im)
!print *,"   tfh=",tfh(im)," dzsh=",dzsh(im),dzsh(im)/tfh(im)
      ELSE 
         ! Unspezifische Hilfskonstruktion fuer neutrale Schichtung:
         DO i=istartpar,iendpar    
            dzsm(i)=z0m(i)*LOG(z1d2*dicke(i,ke)/z0m(i)+z1)
            dzsh(i)=dzsm(i)
         END DO 
         ! Kuerzt sich bei Flussberechnungen wieder heraus!
      END IF
!bug_mod_2011/09/23: dicke(i,ke1) -> dz0(i,mom), dz0(i,sca) }

!     Berechnung der turbulenten Laengenscalen:

      DO i=istartpar,iendpar
         len_scale(i,ke1)=z0m(i)
      END DO
      DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
         DO i=istartpar,iendpar
            IF (cbig(i,k).gt.z0) THEN
!              Die turbulente Laengenskala wird durch die Laengen-
!              skala der lufterfuellten Zwischenraeume limitiert:

               l_turb=z1/(cbig(i,k)*sqrt(z1/exp(rair(i,k))-z1))
               len_scale(i,k)=MIN( dicke(i,k)+len_scale(i,k+1), l_turb )
            ELSE
               len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
            END IF
         END DO
      END DO   
      DO k=kcm-1,1,-1
         DO i=istartpar,iendpar
            len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
         END DO
      END DO   

!     Uebergang von der maximalen turbulenten Laengenskala zur
!     effektiven turbulenten Laengenskala:
!print *,"l_scal=",l_scal

      DO k=ke1,1,-1
!test: different turb. length scale
! IF (k.EQ.ke1) THEN
            DO i=istartpar,iendpar
               lay(i)=l_scal
            END DO
! ELSEIF (k.NE.1) THEN
!    DO i=istartpar,iendpar
!       lay(i)=MAX( l_scal, z1d2*(hhl(i,k-1)-hhl(i,k+1)) )
!    END DO
! END IF    
         DO i=istartpar,iendpar
            len_scale(i,k)=akt*MAX( len_min, &
                                      len_scale(i,k)/(z1+len_scale(i,k)/lay(i)) )
         END DO
      END DO       
!test
!print *,"nach len_scale"

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
               fm2=MAX( grad(u_m)**2+grad(v_m)**2, fc_min )

               ! Vereinfachte Loesung mit Rf=Ri:
               IF (fh2.GE.(z1-rim)*fm2) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=len*(b_2-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2/(fm2-fh2)
                  lm=len*(b_2-(a_6+a_3)*fakt)
                  lh=len*(b_1-a_5*fakt)
               END IF

               val1=lm*fm2; val2=lh*fh2
               wert=MAX( val1-val2, rim*val1 )

               tke(i,k,1)=MAX( vel_min, SQRT(d_m*len*wert) ) !Initialwert fuer SQRT(2TKE)

               tkvm(i,k)=lm*tke(i,k,1)
               tkvh(i,k)=lh*tke(i,k,1)

!              Am Anfang konnte noch keine Advektion oder 
!              Diffusion von q=SQRT(2*TKE) berechnet werden:

               tketens(i,k)=z0
            END DO
 
         END DO   

         DO i=istartpar,iendpar
            tke(i,1,1)=tke(i,2,1)
         END DO

         DO n=2,ntim
            DO k=1,kem
               DO i=istartpar,iendpar
                  tke(i,k,n)=tke(i,k,1)
               END DO
            END DO
         END DO    

      END IF
!print *,"nach tke-init"

! 2)  Berechnung der benoetigten vertikalen Gradienten und
!     Abspeichern auf vari():

!     Am unteren Modellrand:

      DO i=istartpar,iendpar
         lays(i,mom)=z1/dzsm(i)
         lays(i,sca)=z1/dzsh(i)
      END DO
      DO n=1,nmvar
!print *,"n=",n," v(ke)=",vari(im,ke,n)," v(ke1)=",vari(im,ke1,n)
         DO i=istartpar,iendpar
            vari(i,ke1,n)=(vari(i,ke,n)-vari(i,ke1,n))*lays(i,ivtp(n))
         END DO
!print *,"dz=", z1/lays(im,ivtp(n))," grad(ke1)=",vari(im,ke1,n)
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
               DO i=istartpar,iendpar
                  hlp(i,k)=hlp(i,k+1)+vari(i,k,n) &
                                         *(hhl(i,k)-hhl(i,k+1))
               END DO
            END DO

            k1=1
            k2=2
            DO i=istartpar,iendpar
               lays(i,k1)=hhl(i,1)-hhl(i,2)
            END DO

            DO k=2,ke   

!              Berechnung der nicht-lokalen Gradienten als mittlere
!              Differenzenquotienten ueber die stabilitaetsabhaengige
!              Laengenskala in tkvh() bzw tkvm():

!              Speichern der stab.abh. Laengenskala unter lay():

               IF (n.LE.nvel) THEN 
                  DO i=istartpar,iendpar
                     lay(i)=tkvm(i,k)/tke(i,k,nvor)
                  END DO
               ELSE   
                  DO i=istartpar,iendpar
                     lay(i)=tkvh(i,k)/tke(i,k,nvor)
                  END DO
               END IF
                  
!              Bestimmung der nicht-lokalen Gradienten und 
!              Zwischenspeichern derselben auf dem Feld dicke():

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
               DO i=istartpar,iendpar
                  vari(i,k,n)=dicke(i,k)
               END DO
            END DO

         END DO

!        Belegung von dicke() mit den Schichtdicken*rho/dt_tke
!        bzgl. Nebenflaechen:

         DO k=2,ke   
            DO i=istartpar,iendpar
               dicke(i,k)=rhon(i,k)*z1d2*(hhl(i,k-1)-hhl(i,k+1))*fr_tke
            END DO
         END DO      

      ELSE       

!        Berechnung lokaler Gradienten:
         
         DO k=ke,2,-1
            DO i=istartpar,iendpar
               len=(hhl(i,k-1)-hhl(i,k+1))*z1d2
               hlp(i,k)=z1/len
               dicke(i,k)=rhon(i,k)*len*fr_tke
            END DO
         END DO   

         DO n=1,nmvar
            DO k=ke,2,-1
               DO i=istartpar,iendpar
                  vari(i,k,n)=(vari(i,k-1,n)-vari(i,k,n))*hlp(i,k)
               END DO
            END DO   
         END DO      

      END IF
!print *,"nach gradient"

!     Mechanischer Antrieb durch vertikale Scherung:

      DO k=2,kem !von oben nach unten
         DO i=istartpar,iendpar
          ! frm(i,k)=vari(i,k,u_m)**2+vari(i,k,v_m)**2+fc_min
            frm(i,k)=MAX( vari(i,k,u_m)**2+vari(i,k,v_m)**2, fc_min )
         END DO
      END DO    

      !Addition verallgemeinerter Scherterme durch nicht-turbulente subskalige Stroemung:

      DO k=2,kem

         IF (.NOT.lini) THEN !nicht bei der Initialisierung

            IF (lsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
               !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:
                     
!              Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO-Tendenzen:

!achtung
               DO i=istartpar,iendpar
!Achtung: horizontale Pos. beachten!
                 vel1=-(ut_sso(i,k)  *wind(i,k  ,u_m)+vt_sso(i,k)  *wind(i,k  ,v_m))*dp0(i,k-1)
                 vel2=-(ut_sso(i,k-1)*wind(i,k-1,u_m)+vt_sso(i,k-1)*wind(i,k-1,v_m))*dp0(i,k)
   
                 src(i)=MAX( z0, (vel1+vel2)/(dp0(i,k-1)+dp0(i,k)) )
  
!                Beachte: 
!                Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
!                erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
!                Obwohl vel1 und vel2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!
               END DO   

               IF (PRESENT(tket_sso)) THEN
!achtung
                  DO i=istartpar,iendpar
                     tket_sso(i,k)=src(i)
                  END DO
               END IF

!              Addition des Scherterms durch Nachlaufwirbel:

               IF (ltkesso) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
!achtung
                  DO i=istartpar,iendpar
                     frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)
                  END DO
               END IF

            END IF

            IF (lconv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
               !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

!              Addition des Scherterms durch konvektive Zirkulation:

!achtung
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

!Achtung: neue Behandlung der Rauhigkeitsschicht einfuehren
!         wind-tendenzen auf gestaggerte Positionen

!        Vertikalwind und Formreibungskoeff. auf Hauptflaechen:

         DO i=istartpar,iendpar
            lay(i)=z1d2*(w(i,k)+w(i,k+1))                
            src(i)=z1d2*(cbig(i,k)  +csml(i,k) &
                        +cbig(i,k+1)+csml(i,k+1))
         END DO

!        Windbetrag auf der angrenzenden Hauptflaeche:

         IF (k.EQ.kcm) THEN
            DO i=istartpar,iendpar
               velo=z1d2*(w(i,k-1)+w(i,k))
               lays(i,1)=SQRT(wind(i,k-1,u_m)**2+wind(i,k-1,v_m)**2+velo**2)
             END DO 
         END IF

!        Reduzierte Konstanten und implizite Horizontalwind-Tendenzen durch Formreibung:

         DO i=istartpar,iendpar

!           Windbetrag auf der aktuellen Hauptflaeche:
            lays(i,2)=SQRT(wind(i,k,u_m)**2+wind(i,k,v_m)**2+lay(i)**2)

!           Windtendenzen durch Fromreibung:
            wert=src(i)*lays(i,2)
            vars(i,u_m)=-tinc(u_m)*wert*wind(i,k,u_m)/(z1+dt_var*wert)
            vars(i,v_m)=-tinc(v_m)*wert*wind(i,k,v_m)/(z1+dt_var*wert)

!           Windbetrag auf Nebenflaechen:
            can(i,k)=(lays(i,2)*dp0(i,k-1)+lays(i,1)*dp0(i,k))/(dp0(i,k-1)+dp0(i,k))

!           Windbetrag unter Wirkung der Formreibung:
            velo=can(i,k)/(z1+can(i,k)*(cbig(i,k)+csml(i,k))*dt_var)

!           Addition des Scherterms durch Nachlaufwirbel an grossen Rauhigkeitselementen:
            frm(i,k)=frm(i,k)+cbig(i,k)*velo**3/tkvm(i,k)

!           Frequenz der kleinskaligen Rauhigkeitselemente: 
            can(i,k)=csml(i,k)*can(i,k)

         END DO

         DO i=i_st,i_en
            utens(i,k)=utens(i,k)+vars(i,u_m)
            vtens(i,k)=vtens(i,k)+vars(i,v_m)
         END DO

      END DO  !k=kcm,kem !von oben nach unten druch Rauhigkeitsschicht

      !Optionale vertikale Glaettung des mechanischen Antriebs:

      IF (frcsmot.GT.z0) THEN
         DO k=2,kem 
            DO i=istartpar,iendpar
               hlp(i,k)=frm(i,k)
            END DO
         END DO
         CALL vert_smooth ( &
              i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1, k_lw=kem, &
              disc_mom=dicke, cur_tend=hlp ,smo_tend=frm, versmot=frcsmot )
      END IF

!     Thermischer Antrieb:

      DO k=2,ke1 !'frh'(ke1) wird fuer Zirkulationsterm und Temperaturkorrektur benoetigt
         DO i=istartpar,iendpar
            frh(i,k)=a(i,k,4)*vari(i,k,tet_l) &
                      +a(i,k,5)*vari(i,k,h2o_g)
         END DO  
      END DO  

!     Optionale vertikale Glaettung des thermischen Antriebs:

      IF (frcsmot.GT.z0) THEN
         DO k=2,kem
            DO i=istartpar,iendpar
               hlp(i,k)=frh(i,k)
            END DO
         END DO
         CALL vert_smooth ( &
              i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1, k_lw=kem, &
              disc_mom=dicke, cur_tend=hlp ,smo_tend=frh, versmot=frcsmot )
      END IF

!     Belegung von tkvh und tkvm mit den stabilitaetsabhaengigen Laengenmassen:

      DO k=2,kem
         DO i=istartpar,iendpar
            tkvh(i,k)=tkvh(i,k)/tke(i,k,nvor)
            tkvm(i,k)=tkvm(i,k)/tke(i,k,nvor)
         END DO
      END DO

!return

! 3)  Bestimmung von TKE, Stabilitaetsfunktion und Rauhigkeitstermen:

      DO it_durch=it_start, it_end

        !Die Schleife wird nur bei der Initialisierung (d.h. beim ersten Aufruf) wiederholt,
        !um TKE-Gleichgewicht anzunaehern. Die resultierenden TKE-Werte der Zeitstufe 'ntur'
        !gehoeren in diesem Fall dann zur Zeitstufe der uebrigen prognostischen Variablen
        !('nold' bei "leap-frog" oder 'nnow' bei 2-Zeitebenen). 
        !Fuer die folgenden Aufrufe wird die Schleife nur einmal durchlaufen und liefert TKE-Werte
        !die gegenueber den Vorgaengerwerten um einen Zeitschritt weiter in der Zukunft liegen,
        !also wieder zur Zeitstufe der uebrigen prognostischen Variablen gehoeren.

         CALL solve_turb_budgets (1,                                               &
!
                                  istartpar, iendpar, 2, kem,                      &
!
                                  imode_stke=imode_turb,                           &
#ifdef SCLM
                                  it_s=it_durch, grd=vari,                         &
#endif
!-------------------------------------------------------------------------------
                                  fm2=frm, fh2=frh, fcd=can, tls=len_scale, tvt=tketens )

         IF (it_durch.LT.it_end) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF   

      END DO

!return
!     Kein TKE-Gradient am Oberrand:

      DO i=istartpar,iendpar
         tke(i,1,ntur)=tke(i,2,ntur)
      END DO

      IF (iini.EQ.1) THEN !only for separate initialization before the time loop

         DO k=2, kem
            DO i=istartpar,iendpar
               tkvh(i,k)=MAX( con_h, tkhmin, tkvh(i,k)*tke(i,k,ntur) )
               tkvm(i,k)=MAX( con_m, tkmmin, tkvm(i,k)*tke(i,k,ntur) )
! JF:                tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
! JF:                tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
            END DO
         END DO

       ! IF (.NOT.PRESENT(rhoh))  DEALLOCATE ( rhoh , STAT=ilocstat )
       ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
       ! IF (.NOT.PRESENT(dzsm))  DEALLOCATE ( dzsm , STAT=ilocstat )
       ! IF (.NOT.PRESENT(dzsm))  DEALLOCATE ( dzsh , STAT=ilocstat )
       ! IF (.NOT.PRESENT(d_pat)) DEALLOCATE ( dpat,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )
 
         RETURN !finish this subroutine

      END IF   

!  4) Berechnung der effektiven turbulenten Vertikalgradienten
!     und weiterer Turbulenzgroessen:

      DO k=2, ke1

         IF (ltmpcor) THEN

!           Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:
!              
            DO i=istartpar,iendpar
               frm(i,k)=a(i,k,1)*vari(i,k,tet_l)-tet_g &
                         +zlhocp*vari(i,k,liq)               !vert. Temperaturgradient
            END DO   

         !  Dies geschieht schon hier, weil im naechsten Schritt das Feld vari()
         !  durch die effiktiven Gradienten ueberschrieben wird.
         END IF    

         DO i=istartpar,iendpar

!           Berech. der thermodynamischen Flusskonversionsmatrix
!           und der Standardabweichnung des Saettigungsdefizites:

            x1=z1/(z1+zlhocp*a(i,k,3))
            x2=a(i,k,3)*a(i,k,1)
!mod_2011/12/08: temperature conversion of tet-diffusion-tendencies rahter than -fluxesx {
          ! x3=a(i,k,1)*a(i,k,2) !exnr*Cp/Cpd
            x3=a(i,k,2) !Cp/Cpd
!mod_2011/12/08: temperature conversion of tet-diffusion-tendencies rahter than -fluxes }

            clcv = c_scld*rcld(i,k)/(z1+rcld(i,k)*(c_scld-z1))

            rcld(i,k)=SQRT(len_scale(i,k)*tkvh(i,k)*d_h)* &
                        ABS(x2*vari(i,k,tet_l)-vari(i,k,h2o_g))

!           Unter rcld steht jetzt die (geschaetzte) Standardabw.
!           des Saettigungsdefizites!

            flukon33         =z1-clcv*(z1-x1)
            flukon34         =x1*zlhocp/a(i,k,1)*clcv

            flukon53         =-x1*x2*clcv
            flukon54         =x1*clcv
            flukon55         =z1-liqfak

            flukon43         =-flukon53
            flukon44         =z1-flukon54

            flux_3   =(flukon33         *vari(i,k,tet_l) &
                      +flukon34         *vari(i,k,h2o_g)) * x3
            flux_4   =(flukon43         *vari(i,k,tet_l) &
                      +flukon44         *vari(i,k,h2o_g))
            flux_5   =(flukon53         *vari(i,k,tet_l) &
                      +flukon54         *vari(i,k,h2o_g) &
                      +flukon55         *vari(i,k,liq)) 
!if (i.eq.im .and. k.ge.ke) then
!  print *,"k=",k," grad_tet_l=",vari(i,k,tet_l)," grad_tet=",flux_3
!  print *,"k=",k," grad_h2o_g=",vari(i,k,h2o_g)," grad_vap=",flux_4
!  print *,"k=",k," grad_liq_0=",vari(i,k,liq  )," grad_liq=",flux_5
!endif
            vari(i,k,tet)=flux_3 !(Cp/Cpd)*grad(tet)
            vari(i,k,vap)=flux_4
            vari(i,k,liq)=flux_5

            !Achtung: 
            !'vari(i,k,tet)' ist der in der Teta-gleichung benoetigte
            !effective Gradient und entspricht (Cp/Cpd)*tet_flux_dens.
            !Die Indices 'tet', 'tem' und 'tet_l' haben den gleichen Wert
            !so wie auch 'vap' und 'h2o_g' den gleichen Wert haben.
            !Die Unterscheidungen sollen nur den jeweiligen Inhalt verdeutlichen.

            !Im Falle "icldm_turb.NE.-1" ist "liqfak.EQ.z1" und 'flukon55' verschwindet.
            !Im Falle "icldm_turb.EQ.-1" dagegen verschwinden 'flukon53' und 'flukon54'.
            !Dann ist aber "flukon55.EQ.z1". Ferner verschwinden in diesem Fall auch
            !'flukon34' und 'flukon43', wobei die physikalischen Groessen 'tet_l' und 'tet'
            !sowie 'h2o_g' und 'vap' identisch sind.

            a(i,k,3)=clcv !'a(i,k,3)' enthaelt ab jetzt den Wolkenbedeckungsgrad!

         END DO

      END DO   

      DO k=2, ke    
         DO i=istartpar,iendpar
!           Berechn. der Diffusionskoeffizienten:

!Achtung: moeglichst ohne tk?min
!___________________________________________________________________________
!test:
            tkvh(i,k)=MAX( con_h, tkhmin, tkvh(i,k)*tke(i,k,ntur) )
            tkvm(i,k)=MAX( con_m, tkmmin, tkvm(i,k)*tke(i,k,ntur) )

!test: ohne tk?min (wird aber bei Diffusion benutzt!)
! tkvh(i,k)=MAX( con_h, tkvh(i,k)*tke(i,k,ntur) )
! tkvm(i,k)=MAX( con_m, tkvm(i,k)*tke(i,k,ntur) )
            
!test: ohne tk?min und ohen lam. Diffkof. (wird aber bei Diffusion benutzt!)
! tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
! tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)

!test: nur minimaler Diff.koef.
! tkvh(i,k)=tkhmin
! tkvm(i,k)=tkmmin
!___________________________________________________________________________

!           tkvh und tkvm enthalten jetzt nicht mehr Diffusions-
!           laengen, sondern wieder Diffusionskoeffizienten in m^2/s!

         END DO
      END DO

!  5) Untere Randwerte der Turbulenzgroessen:

      IF (itype_tran.NE.2) THEN 
!        Es wurde i.a. noch keine unteren Randwerte fuer 
!        tke(), sowie tkvm(), tkvh() und edr() berechnet:

         DO i=istartpar,iendpar
!Achtung!
!bug_mod_2011/09/23: dicke(i,ke1) -> dz0(i,mom), dz0(i,sca) {
            tkvm(i,ke1)=vh0(i)*dzsm(i)*tcm(i)
            tkvh(i,ke1)=vh0(i)*dzsh(i)*tch(i)
!bug_mod_2011/09/23: dicke(i,ke1) -> dz0(i,mom), dz0(i,sca) }

!           Weil in partur(b,s) keine laminare Grenzschicht
!           beruecksichtigt wird, sind die Reduktionskoeff.
!           immer gleich 1:

            tfh(i)=z1
            tfm(i)=z1

!           Der q=SQRT(2*TKE)-Wert am unteren Modellrand wird mit 
!           Hilfe von c_tke*Ustar bestimmt: 

            tke(i,ke1,ntur)=c_tke*vh0(i)*SQRT(tcm(i))
      
         END DO    

         IF (PRESENT(edr).OR.ltmpcor) THEN
            IF (PRESENT(edr)) THEN
               ediss => edr
            ELSE
               ediss => tketens
            END IF
            DO i=istartpar,iendpar
               ediss(i,ke1)=tke(i,ke1,ntur)**3/(d_m*len_scale(i,ke1))
            END DO
         END IF   

      END IF

      IF (lsfli(tem)) THEN !use explicit shfl_s
         DO i=istartpar,iendpar
            vari(i,ke1,tet)=shfl_s(i)/(cp_d*rhon(i,ke1)*tkvh(i,ke1)*a(i,ke1,1))
         END DO
         !Note:'vari(tet)' belongs to potential temperature,
         !     and 'a(1)' contains Exner factor on half levels.
      END IF
      IF (lsfli(vap)) THEN !use explicit lhfl_s
         DO i=istartpar,iendpar
            vari(i,ke1,vap)=lhfl_s(i)/(lh_v*rhon(i,ke1)*tkvh(i,ke1))
         END DO
      END IF
      !Note: "tkvh(ke1) > 0.0" is required!

! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperaturtendenzen 
!     (ausser der Divergenz des Zirkulationsflusses):

      IF (ltmpcor) THEN

          DO k=2, ke1
            DO i=istartpar,iendpar
               thermik=tkvh(i,k)*frh(i,k)
               phasdif=tkvh(i,k)*frm(i,k) &
                      *(zrcpv*vari(i,k,vap)+zrcpl*vari(i,k,liq))  
  
               tketens(i,k)=len_scale(i,k)/a(i,k,2) &
                       *((ediss(i,k)+thermik)/cp_d+phasdif)
            END DO   

!           Beachte:
!           In 'frm' steht an dieser Stelle der Temp.-Gradient!
!           Wegen der spaeteren Interpolation auf Hauptflaechen,
!           wird mit der turbulenten Laengenskala multipliziert.
         END DO   
   
!achtung
         DO i=i_st,i_en
            ttens(i,1)=ttens(i,1)+tinc(tem)*tketens(i,2) &
                         /(len_scale(i,1)+len_scale(i,2))

         END DO
         DO k=2,ke
!achtung
            DO i=i_st,i_en
               ttens(i,k)=ttens(i,k)+tinc(tem)*(tketens(i,k)+tketens(i,k+1)) &
                            /(len_scale(i,k)+len_scale(i,k+1))
            END DO
         END DO
      END IF   

! 7)  Bestimmung des Zirkulationstermes:
!Achtung: Zirkulationsterm revidieren:

      IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden

         DO i=istartpar,iendpar
           l_pat(i)=l_hori*dpat(i)/(l_hori+dpat(i))
         END DO

         DO k=2,ke1

            IF (k.LT.ke1) THEN
!              Interpolation des Druckes auf die naechst hoehere 
!              Nebenflaeche:

               DO i=istartpar,iendpar
                  lay(i)=(prs(i,k)*dp0(i,k-1)+prs(i,k-1)*dp0(i,k)) &
                          /(dp0(i,k)+dp0(i,k-1))
               END DO   
            ELSE
               DO i=istartpar,iendpar
                  lay(i)=ps(i)
               END DO
            END IF

            DO i=istartpar,iendpar

               fakt=z1-z2*ABS(a(i,k,3)-z1d2)
               len=MAX( l_pat(i), SQRT(fakt*len_scale(i,k)*l_hori) )

!              Berechnung der lokalen Kohaerenzlaenge fuer die Parameterisierung
!              des Zirkulationstermes (kohae_len=l_pat*exnr*grad(tet_v)*r/g):

               fakt=frh(i,k)*lay(i)/(rhon(i,k)*grav**2)
               kohae_len=len*SIGN(z1,fakt)*MIN( ABS(fakt), z1 )

!              Belegung von 'tketens' mit der TKE-Flussdichte des Zirkulationstermes
!              (Zirkulationsfluss oder Drucktransport):

               tketens(i,k)=tkvh(i,k)*kohae_len*frh(i,k)

!              Die Divergenz dieser Flussdichte ist gleichzeitig
!              eine weitere Quelle fuer thermische Energie.

            END DO
!print *,"k=",k
!print *,"len_scale=",MINVAL(len_scale(:,k)),MAXVAL(len_scale(:,k)), &
!         " tketens=",MINVAL(tketens(:,k)),  MAXVAL(tketens(:,k))

!           Addition des Teta-Gradienten, welcher zur Teta-Flussdichte 
!           durch den Zirkulationsterm gehoert:

            IF (ltmpcor) THEN
               DO i=istartpar,iendpar
!bug_2012/01/18: Division durch "a(i,k,1)*a(i,k,2)" fehlte {
                  vari(i,k,tet)=vari(i,k,tet) &
                                 -tketens(i,k)/(tkvh(i,k)*a(i,k,1)*a(i,k,2)*cp_d)
!bug_2012/01/18: Division durch "a(i,k,1)*a(i,k,2)" fehlte }
               END DO
            END IF

            !Die kinetische Energie der Zirkulation wird der thermischen Energie
            !entzogen!

!print *,"k=",k," tketens=",tketens(im,k)
         END DO
!read *,xx
         
      END IF   

! 14) Berechnung der Diffusionstendenz von q=SQRT(2*TKE):

      IF (ldotkedif) THEN 

!        Vorbereitung zur Bestimmung der vertikalen Diffusions-Tendenzen von q**2:

         expl_mom => a(:,:,1)
         impl_mom => a(:,:,2)
         invs_mom => a(:,:,3)
         upd_prof => a(:,:,4)
         diff_dep => a(:,:,5)
         cur_prof => hlp

         DO k=2, ke1
            DO i=istartpar,iendpar
!___________________________________________________________________________
!test: TKE-Diffusion mit Stab.fnkt. fuer Skalare:
               cur_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)
! cur_prof(i,k)=c_diff*tkvh(i,k)
!___________________________________________________________________________
            END DO
         END DO
         !'cur_prof' enthaelt jetzt den Diff.koef. fuer q**2 auf Nebenflaechen

!print *,"k=",2," cur_prof=",cur_prof(im,2)
         DO k=3, ke1
            DO i=istartpar,iendpar
               diff_dep(i,k)=hhl(i,k-1)-hhl(i,k)
               expl_mom(i,k)=rhoh(i,k-1)*z1d2*(cur_prof(i,k-1)+cur_prof(i,k))/diff_dep(i,k)
               impl_mom(i,k)=rhoh(i,k-1)*diff_dep(i,k)*fr_tke
            END DO
!print *,"k=",k," cur_prof=",cur_prof(im,k)," a5=",diff_dep(im,k)
         END DO
!read *,xx

!        Belegung der Eingangsprofile mit (modifizierten) TKE(=z1d2*q**2)-Werten:

         IF (lcircterm) THEN

!           Belegung mit durch den Zirkulationsterm modifizierten (effektiven) TKE-Werten:

            !Beachte:
            !'tketens' ist eine TKE-Flussdichte, deren Vertikalprofil im wesentlichen durch
            !d_z(tet_v)**2 bestimmt ist, was zumindest in der Prandtl-Schicht prop. zu 1/len_scale ist.
            !Die Interpolation auf Hauptflaechen erfolgt daher mit 'tketens*len_scale'.
            !Anderenfalls ist mit grossen Interpolationsfehlern zu rechnen:

            DO i=istartpar,iendpar  
               tketens(i,2)=rhon(i,2)*tketens(i,2)*len_scale(i,2)
               lay(i)=z1d2*tke(i,2,ntur)**2
               cur_prof(i,2)=lay(i)
            END DO
!print *,"k=",2," tketens=",tketens(im,2)," tke=",tke(im,2,ntur)**2,cur_prof(im,2)
            DO k=3,ke1
               DO i=istartpar,iendpar  
                  tketens(i,k)=rhon(i,k)*tketens(i,k)*len_scale(i,k)
                  cur_prof(i,k)=cur_prof(i,k-1)-lay(i) &
                          +(tketens(i,k)+tketens(i,k-1)) &
                           /((len_scale(i,k)+len_scale(i,k-1))*expl_mom(i,k))
                  lay(i)=z1d2*tke(i,k,ntur)**2
                  cur_prof(i,k)=cur_prof(i,k)+lay(i)
               END DO
!print *,"k=",k," tketens=",tketens(im,k)," tke=",tke(im,k,ntur)**2,cur_prof(im,k)
!print *," d_cur_prof=",z1d2*(tketens(im,k)+tketens(im,k-1))*diff_dep(im,k)
            END DO
!read *,xx

         ELSE

!           Belegung mit unveraenderten TKE-Werten:

            DO k=2, ke1
               DO i=istartpar,iendpar  
                  cur_prof(i,k)=z1d2*tke(i,k,ntur)**2
               END DO
            END DO

         END IF     
         !'cur_prof' enthaelt jetzt das (modifizierte) TKE-Profil und
         !'tketens' ist nun wieder frei.


         !In den Diffusionroutinen wird vorausgesetzt, dass ein Flussniveau mit gleichem
         !Vertikalindex wie ein Konzentrationsniveau gerade ueber diesem liegt.
         !Die bisherige Hauptflaechenindizierung musste daher fuer die Uebergabefelder
         !der Routinen 'prep_impl_vert_diff' und 'calc_impl_vert_diff' um "1" verschoben werden.

         CALL prep_impl_vert_diff( lsflucond=.FALSE.,                                  &
              i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1,                          &
              disc_mom=dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
              scal_fac=frm   )

!        Berechnung der vertikalen Diffusionstendenzen von TKE=z1d2*q**2:

         CALL calc_impl_vert_diff ( &
              i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1, itndcon=0,               &
              disc_mom=dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
              scal_fac=frm, cur_prof=cur_prof(:,:), upd_prof=upd_prof, eff_flux=tketens(:,:) )

         !'upd_prof' enthaelt jetzt die mit der Diffusionstendenz aufdatierten (modifizierte) TKE-Werte,
         !'tketens' die effektiven Flussdichten (positiv abwaerts) der (semi-)impliziten Vertikaldiffusion.
         
!        Zuschlag durch Volumenterm innerhalb der Rauhigkeitsschicht:

         DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
            DO i=istartpar,iendpar
               wert=(tketens(i,k)*dp0(i,k)+tketens(i,k+1)*dp0(i,k-1))/(dp0(i,k)+dp0(i,k-1))
                   !effektive TKE-Flussdichte interpoliert auf die k-te Nebenflaeche,
                   !wobei sich 'tketens(:,k)' auf die (k-1)-te Hauptflaeche bezieht!
               upd_prof(i,k)=upd_prof(i,k)+dt_tke*wert*z1d2*(rair(i,k-1)-rair(i,k+1))/dicke(i,k)
            END DO
         END DO   
         !'upd_prof' enthaelt die mit dem Volunmenterm beaufschlagten (modifizierten) TKE-Werte.

!        Speichern der zugehoerigen q-Tendenzen:

         DO k=2,ke
            DO i=istartpar,iendpar
!___________________________________________________________________________
!test: q-Tendenz ueber Zwischenwerte von q:
               tketens(i,k)=(SQRT(z2*upd_prof(i,k))-SQRT(z2*cur_prof(i,k)))*fr_tke
! JF:                tketens(i,k)=(upd_prof(i,k)-cur_prof(i,k))*fr_tke/tke(i,k,ntur)
!modif:Einschraenkung der erweiterten Diffusionstendenzen von TKE
               wert=tke(i,k,ntur)*fr_tke
               tketens(i,k)=MAX( -wert, tketens(i,k) )
! JF:                tketens(i,k)=MAX( -wert, MIN( wert, tketens(i,k) ) )
!modif
!test
!___________________________________________________________________________
            END DO

         END DO

         !Am Unterrand gibt es keine q-Tendenz durch Diffusionsterme:
         DO i=istartpar,iendpar
            tketens(i,ke1)=z0
         END DO

      ELSE !.NOT.ldotkedif      

!        Zuruecksetzen der q-Tendenzen:

         DO k=2,ke1
            DO i=istartpar,iendpar
               tketens(i,k)=z0
            END DO
         END DO

      END IF !ldotkedif

!---------------------------------------------------------------------- 
! 15) Berechnung der turb. horizontalen TKE-Flussdichtedivergenzen
!     und Aufsummieren auf das SQRT(2*TKE)-Tendenzfeld:

!        ** wird erst spaeter eingefuehrt **
!---------------------------------------------------------------------- 

! 13) Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!     des Saettigungsdefizites:

      DO i=istartpar,iendpar  
         rcld(i,1)=rcld(i,2)
      END DO
      DO k=2,kem-1
         DO i=istartpar,iendpar  
            rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
         END DO
      END DO
!     Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
!     der Wert auf der entspr. Nebenflaeche beibehalten.

      DO i=istartpar,iendpar
         dzsm(i)=dzsm(i)/tfm(i)
         dzsh(i)=dzsh(i)/tfh(i)
      END DO

!--------------------     
   ELSE !.NOT.lturatm
!--------------------     

      IF (.NOT.PRESENT(rho)) THEN
         DO k=1,ke
            DO i=istartpar,iendpar    
               virt=z1+rvd_m_o*qv(i,k)-qc(i,k) !virtueller Faktor
               rhoh(i,k)=prs(i,k)/(r_d*virt*t(i,k))
            END DO
         END DO   
      END IF

      CALL bound_level_interp( istartpar,iendpar, 2,ke, &
!___________________________________________________________________________
!test: mass weighted interpolation
!                              nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=hhl, auxil=hlp)
                               nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dp0)
!___________________________________________________________________________

      IF (.NOT.PRESENT(epr)) THEN
         DO k=1, ke
            DO i=istartpar,iendpar
               exner(i,k)=zexner(prs(i,k))
            END DO
         END DO
      END IF   

      DO i=istartpar,iendpar    
         vel1=u(i,ke)
         vel2=v(i,ke)
         velo=MAX( vel_min, SQRT(vel1**2+vel2**2) )

         dzsm(i)=tkvm(i,ke1)/(velo*tcm(i))
         dzsh(i)=tkvh(i,ke1)/(velo*tch(i))
      END DO 

!-----------------------
   END IF !lturatm
!-----------------------

!     Berechnung der Vertikaldiffuion von Modellvariablen auf Hauptflaechen:

      IF (ldovardif .OR. ldogrdcor) THEN !Vertikaldiffusion wird hier berechnet

         !Berechnung der Luftdichte am Unterrand:

         DO i=istartpar,iendpar
            virt=z1+rvd_m_o*qv_s(i) !virtueller Faktor
            rhon(i,ke1)=ps(i)/(r_d*virt*t_g(i))
         END DO   

         !Setzen von Steuerparametern:

         IF (ldogrdcor) THEN
            IF (lnonloc) THEN
               ncorr=1
            ELSE
               ncorr=nvel+1
            END IF
         ELSE
            ncorr=ndiff+1
         END IF

         IF (lmomdif) THEN !wind diffusion required
            nprim=1
         ELSEIF (lscadif) THEN
            nprim=nvel+1
         ELSE
            nprim=ndiff+1
         END IF
         nprim=MIN(nprim, ncorr)

         IF (lscadif) THEN !scalar diffusion required
            nlast=ndiff
         ELSEIF(lmomdif) THEN
            nlast=nvel
         ELSE
            nlast=0
         END IF
         IF (ldogrdcor) nlast=MAX(nlast, nmvar)

         ivtype=0

!-----------------------------------------------------------------        
         DO n=nprim, nlast !loop over all variables to be diffused

            m=MIN(n,nmvar)

            IF (ivtype.EQ.0) THEN
               linisetup=.TRUE.; lnewvtype=.TRUE.
            ELSE
               linisetup=.FALSE.;lnewvtype=.FALSE.
            END IF

            IF (n.LE.nvel) THEN !a wind component

               istart=iid(n,1); iend=iid(n,2)

               IF (ivtype.EQ.sca) lnewvtype=.TRUE.

               lsflucond=.FALSE. !no lower flux condition for momentum!

               ivtype=mom
            ELSE
               istart=i_st; iend=i_en

               IF (ivtype.EQ.mom) THEN
                  lnewvtype=.TRUE.
               END IF

               lsflucond=lsflcnd !use chosen type of lower boundary condition

               ivtype=sca
            END IF

            IF (n.LT.ncorr .OR. n.GT.nmvar) THEN
               igrdcon=0 !keine Gradientkorrektur der Profile
            ELSEIF ((.NOT.lscadif .AND. ivtype.EQ.sca) .OR. (.NOT.lmomdif .AND. ivtype.EQ.mom)) THEN
               igrdcon=1 !nur Profile aus Gradientkorrektur
            ELSE
               igrdcon=2 !korrigierte Profile aus effektiven Gradienten
            END IF

            kgc=2 !uppermost level of gradient correction

            IF (.NOT.ltend(n)) THEN !tendency array not present
               itndcon=0 !no explicit tendency correction
            ELSE
               itndcon=itnd !use chosen mode of tendency correction
            END IF
!test: never expl tendency correction
!itndcon=0
!test

            IF (lsfli(n) .AND. (.NOT.lturatm .OR. n.NE.tem .OR. n.NE.vap)) THEN
               !Load effective surface layer gradients due to given flux values:

               DO i=istart,iend
                  vari(i,ke1,m)=dvar(n)%sv(i)/(rhon(i,ke1)*vtyp(ivtype)%tkv(i,ke1))
               END DO
               IF (n.EQ.tem) THEN !flux density is that of sensible heat
                  DO i=istart,iend
                     vari(i,ke1,m)=vari(i,ke1,m)/(cp_d*zexner(ps(i)))
                  END DO
               ELSEIF (n.eq.vap) THEN !flux density is that of latent heat
                  DO i=istart,iend
                     vari(i,ke1,m)=vari(i,ke1,m)/lh_v
                  END DO
               END IF
               !Note:
               !In this case the surface concentration is not used in 'vert_grad_diff'!
               !'tkv(ke1)' needs to be >0, which is always the case, if calculated by 'turbtran'!
               !Thus "tkvh(ke1)=0.0" should not be forced, if "lsfli=.TRUE"!
               !For tracers it is "m=nmvar"!
               !In case of "lturatm=T" 'vari(ke1)' has already been loaded using shlf_s or lhfl_s!
            END IF      

!           Belegung der Eingangsprofile und -Tendenzen:

            IF (n.EQ.tem) THEN !temperature that needs to be transformed

               cur_prof => hlp

               DO k=1,ke
                  DO i=istart,iend
                     cur_prof(i,k)=t(i,k)/exner(i,k) !potential temperature
                  END DO
               END DO
               DO i=istart,iend
                  cur_prof(i,ke1)=t_g(i)/zexner(ps(i)) 
               END DO
               IF (itndcon.GT.0) THEN !explicit tendencies to be considered
                  DO k=1,ke
                     DO i=istart,iend
                        dicke(i,k)=ttens(i,k)/exner(i,k)
                     END DO
                  END DO
               END IF
            
            ELSE !any other variable that needn't to be transformed

               cur_prof => hlp

               DO k=1,ke
                  DO i=istart,iend
                     cur_prof(i,k)=dvar(n)%av(i,k)
                  END DO
               END DO

               IF (ASSOCIATED(dvar(n)%sv)) THEN !surface variable is present
                  DO i=istart,iend
                     cur_prof(i,ke1)=dvar(n)%sv(i)
                  END DO
               ELSEIF (n.LE.nvel) THEN !a momentum variable
                  DO i=istart,iend
                     cur_prof(i,ke1)=z0 !no slip condition
                  END DO
               ELSE !enforce a zero flux condition as a default
                  DO i=istart,iend
                     cur_prof(i,ke1)=cur_prof(i,ke)
                  END DO
               END IF
               IF (itndcon.GT.0) THEN !explicit tendencies have to be considered
                  DO k=1,ke
                     DO i=istart,iend
                        dicke(i,k)=dvar(n)%at(i,k)
                     END DO
                  END DO
               END IF

            END IF   

!           Berechnung der vertikalen Diffusionstendenzen:

            CALL vert_grad_diff( kcm, kgc,                         &
 
                 istart, iend, 0, ke1,                             &
!
                 dt_var, ivtype, igrdcon, itndcon,                 & 
                 linisetup, lnewvtype, lsflucond, lsfli(n),        & 
!
                 rho=rhoh, rho_n=rhon, hhl=hhl, r_air=rair,        &

                 tkv=vtyp(ivtype)%tkv, dzs=vtyp(ivtype)%dzs,       &
!
                 disc_mom=a(:,:,1), expl_mom=a(:,:,2),             &
                 impl_mom=a(:,:,3), invs_mom=a(:,:,4),             &
                 diff_mom=frh, scal_fac=frm, diff_dep=a(:,:,5),    &
!
                 cur_prof=cur_prof, dif_tend=dicke,                &
                 eff_flux=vari(:,:,m) )

            !Beachte:
            !'frh' und 'frm' sind genauso wie 'a(:,:,1:5)' reine Hilfsspeicher in 'vert_grad_diff'.
            !Weil Fluesse ab "n>=liq=nmvar" nicht mehr benoetigt werden, bleibt
            !'vari' nur bis 'nmvar' dimensioniert und vari(nmvar) wird auch fuer "n>nmvar" benutzt.

!           Sichern der Tendenzen:
       
            IF (n.EQ.tem) THEN
               DO k=1,ke
                  DO i=istart,iend
!mod_2011/12/08: temperature conversion of tet-diffusion-tendencies rahter than -fluxes {
                     dvar(n)%at(i,k)=dvar(n)%at(i,k)+exner(i,k)*dicke(i,k)*tinc(n)
!mod_2011/12/08: temperature conversion of tet-diffusion-tendencies rahter than -fluxes }
                  END DO
               END DO
            ELSE
               DO k=1,ke
                  DO i=istart,iend
                     dvar(n)%at(i,k)=dvar(n)%at(i,k)+dicke(i,k)*tinc(n)
                  END DO
               END DO
            END IF

            IF (n.EQ.vap .AND. PRESENT(qv_conv)) THEN
               DO k=1,ke
                  DO i=istart,iend
                     qv_conv(i,k)=qv_conv(i,k)+dicke(i,k) !always a tendency
                  END DO
               END DO
            END IF

         END DO !nprim, nlast
!-----------------------------------------------------------------        

!DO k=1,ke
! print *,"k=",k," t=",t(im,k),ttens(im,k)
! print *,"k=",k," qv=",qv(im,k),qvtens(im,k)
! print *,"k=",k," qc=",qc(im,k),qctens(im,k)
! print *," u=",u(im,k),utens(im,k)
! print *," v=",v(im,k),vtens(im,k)
! print *,"k=",k," dvar(6)=",dvar(6)%av(im,k),dvar(6)%at(im,k)
! dvar(6)%at(im,k)=z0
!ENDDO
!read *,xx

      ELSE !Vertikaldiffusion wird hier nicht berechnet

         !Berechnung der Oberflaechenflussdichten (poitiv abwaerts):

         DO n=1,nmvar 

            IF (n.LE.nvel) THEN
               ivtype=mom
               istart=iid(n,1); iend=iid(n,2)
            ELSE
               ivtype=sca
               istart=i_st; iend=i_en
            END IF

            DO i=istart,iend
               vari(i,ke1,n)=rhon(i,ke1)*vtyp(ivtype)%tkv(i,ke1)*vari(i,ke1,n)
            END DO

         END DO

      END IF
!return

!     Berechnung der Enthalpieflussdichten:

!Achtung: Ev. in terra benutzen oder implizite Oberlaechentemperatur
!Ist cp-Fluss tatsaechlich der thermische Erdbodenantrieb?
!Was gilt im Falle der T-Gleichung in cv-Form?

      IF (.NOT.lsfluse .OR. .NOT.lsflcnd) THEN !effektive Oberfl.flussdichten wurden neu bestimmt

        IF (PRESENT(shfl_s)) THEN
           DO i=i_st,i_en
              shfl_s(i)=zexner(ps(i))*cp_d*vari(i,ke1,tet)
           END DO
        END IF
        IF (PRESENT(lhfl_s)) THEN
           DO i=i_st,i_en
              lhfl_s(i)=lh_v*vari(i,ke1,vap)
           END DO
        END IF
 
!---------------------------------------------------------------------------------------
#ifdef SCLM 
        IF (lsclm .AND. latmflu) THEN
           SHF%mod(0)%val=zexner(ps(im))*cp_d*vari(im,ke1,tet); SHF%mod(0)%vst=i_cal
           LHF%mod(0)%val=                  lh_v*vari(im,ke1,vap); LHF%mod(0)%vst=i_cal
!print *,"in turbdiff: latmflu=",latmflu," shf=",SHF%mod(0)%val
!print *,"in turbdiff: latmflu=",latmflu," lhf=",LHF%mod(0)%val
!read *,xx
           !Note:
           !IF ".NOT.latmflu", SHF and LHF either are loaded by the fluxes used for the soil budget (lertflu)
           !or they have been loaded above by the explicit SHF and LHF at the surface (lsurflu).
           !SHF and LHF are positive downward and they may have been corrected with 
           !vertical integrated correction tendencies.
           !Thus they always refer to the used flux densities , which are only then equal to the explicit surface
           !flux density, if a lower flux condition is used "lsflcnd=.TRUE.".
        END IF
#endif 
!SCLM-----------------------------------------------------------------------------------

        !Bem: shfl_s und lhfl_s, sowie SHF und LHF sind positiv abwaerts!
      END IF

! 11) Berechnung der Tendenzen infloge horizontaler turb. Diffusion
!     (und aufsummieren auf die Tendenzfelder):

! 16) Deallocierung lokaler dynamischer Felder:
     
    ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
    ! IF (.NOT.PRESENT(d_pat)) DEALLOCATE ( dpat,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )

END SUBROUTINE turbdiff

!********************************************************************************

SUBROUTINE adjust_satur_equil ( ktp,             &
!
   i_st, i_en, k_st, k_en,                       &
!  
   lcalrho, lcalepr, lcaltdv, lpotinp, ladjout,  &
   icldmod,                                      &
   zrcpv, zrcpl, zlhocp,                         &
!
   prs, t, qv, qc, fip,                          & 
   exner, rcld, dens,                            &
   r_cpd, qst_t, g_tet, g_h2o, tet_l, q_h2o, q_liq )

INTEGER (KIND=iintegers), INTENT(IN) :: &
  ktp,                    & !start index of vertical dimension
!
  i_st, i_en,             & !horizontal start- and end-indices
  k_st, k_en                !vertical   start- and end-indices

LOGICAL, INTENT(IN) :: &
  lcalrho, &   !density calculation required
  lcalepr, &   !exner pressure calculation required
  lcaltdv, &   !calculation of special thermodynamic variables required
  lpotinp, &   !input temperature is a potential one
  ladjout      !output of adjusted variables insted of conserved ones

INTEGER (KIND=IINTEGERS), INTENT(IN) :: &
  icldmod      !mode of water cloud representation

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(IN) :: &
  prs, t, qv   !current pressure, temperature and water vapor content

REAL (KIND=ireals), DIMENSION(:,ktp:), OPTIONAL, INTENT(IN) :: &
  qc           !current cloud water content

REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
  fip          !interpolation factor with respect to an reference level

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(INOUT) :: &
  exner, &     !current Exner-factor
  rcld         !inp: standard deviation of oversaturation
               !out: saturation fraction

REAL (KIND=ireals), DIMENSION(:,ktp:), TARGET, INTENT(INOUT) :: &
  tet_l, &     !inp: liquid water potent. temp. (only if 'fip' is present)
               !out: liquid water potent. temp. (or adjust. 't' , if "ladjout")
  q_h2o, &     !inp: total  water content (only if 'fip' is present)
               !out: total  water content       (or adjust. 'qv', if "ladjout") 
  q_liq        !out: liquid water content after adjustment

REAL (KIND=ireals), DIMENSION(:,ktp:), OPTIONAL, INTENT(OUT) :: &
  dens,  &     !current air density 
  r_cpd        !cp/cpd

REAL (KIND=ireals), DIMENSION(:,ktp:), TARGET, INTENT(OUT) :: &
  qst_t, &     !out: d_qsat/d_T
               !aux:                                    TARGET for 'qvap', if ".NOT.ladjout"
  g_tet, &     !out: g/tet-factor of tet_l-gradient in buoyancy term
               !aux: total  water content of ref. lev.; TARGET for 'temp', if ".NOT.ladjout"
  g_h2o        !out: g    -factpr of q_h20-gradient in buoyancy term
               !aux: liq. wat. pot. temp. of ref. lev.; TARGET for 'virt'

REAL (KIND=ireals), INTENT(IN) :: &
  zrcpv,  &    !0-rcpv  swith for cp_d/cp_v - 1
  zrcpl,  &    !0-rcpl  swith for cp_d/cp_l - 1
  zlhocp       !0-lhocp switch for ratio between of latent heat of evap and heat capacity

REAL (KIND=ireals) :: &
  teta, pdry,  &     !corrected pot. temp. and partial pressure of dry air
  clcv, fakt         !cloud cover and auxilary factor

REAL (KIND=ireals), DIMENSION(:,:), POINTER :: &
  temp, &            !corrected temperature
  qvap, &            !corrected water vapour content
  virt               !reciprocal virtual factor

INTEGER (KIND=iintegers) :: &
  i,k

!-------------------------------------------------------------------------------

   !Save conserved variables at reference levels for interpolation:
   IF (PRESENT(fip)) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            g_h2o(i,k)=q_h2o(i,k) !tot. wat. cont.
            g_tet(i,k)=tet_l(i,k) !liq. wat. pot. temp.
         END DO  
      END DO
   END IF     

   !Conserved variables at given levels (with respect to phase change):
   IF (icldmod.EQ.-1 .OR. .NOT.PRESENT(qc)) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            q_h2o(i,k)=qv(i,k) !tot. wat. cont.
            tet_l(i,k)= t(i,k) !liq. wat. temp.
         END DO
      END DO
   ELSE
      DO k=k_st, k_en
         DO i=i_st,i_en
            q_h2o(i,k)=qv(i,k) +       qc(i,k) !tot. wat. cont.
            tet_l(i,k)= t(i,k) - lhocp*qc(i,k) !liq. wat. temp.
         END DO 
      END DO
   END IF

   !Calculation of Exner-pressure:
   IF (lcalepr) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            exner(i,k)=zexner(prs(i,k))
         END DO 
      END DO 
   END IF   

   !Transformation in real liquid water temperature:
   IF (lpotinp) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            tet_l(i,k)=exner(i,k)*tet_l(i,k)
         END DO
      END DO    
   END IF    
      
   !Interpolation of conserved variables (with respect ot phase change)
   !onto levels of interest:
   IF (PRESENT(fip)) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            q_h2o(i,k)=             g_h2o(i,k)*(z1-fip(i))+q_h2o(i,k)*fip(i)
            tet_l(i,k)=exner(i,k)*g_tet(i,k)*(z1-fip(i))+tet_l(i,k)*fip(i)
         END DO

         !The layer between the reference level and the current level is assumed
         !to be a pure CFL layer without any mass. Consequently the Exner pressure
         !values at both levels are treated like being equal.
      END DO
   END IF

   !Berechnung der effektiven Bedeckung mit Wasserwolken:
   IF (icldmod.LE.0) THEN
      !Keine Wolkenberuecksichtigung;
      !alles Wolkenwasser wird ignoriert oder verdunstet:

      DO k=k_st, k_en
         DO i=i_st,i_en
            rcld(i,k)=z0
           q_liq(i,k)=z0
         END DO
      END DO

   ELSEIF (icldmod.EQ.1) THEN
      !Wolkenwasser wird als skalige Wolke interpretiert:

      IF (PRESENT(qc)) THEN
         DO k=k_st, k_en
            DO i=i_st,i_en
               IF ( qc(i,k) .GT. z0 ) THEN
                  rcld(i,k) = z1
               ELSE
                  rcld(i,k) = z0
               END IF
               q_liq(i,k)=qc(i,k)
            END DO
         END DO
      ELSE   
         CALL turb_cloud (ktp,                    &
              i_st,i_en, k_st,k_en,               &
              icldtyp=0,                          &
              prs=prs, ps=ps, t=tet_l, qv=q_h2o,  &
              clcv=rcld, clwc=q_liq                 )
      END IF      

   ELSEIF (icldmod.EQ.2) THEN
      !Spezielle Diagnose von Wasserwolken;
      !unter Breruecksichtigung subskaliger Kondensation:

       CALL turb_cloud (ktp,                   &
            i_st,i_en, k_st,k_en,              &
            icldtyp=itype_wcld,                &
            prs=prs, ps=ps, t=tet_l, qv=q_h2o, &
            clcv=rcld, clwc=q_liq                )
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
      DO k=k_st, k_en
         DO i=i_st,i_en
            temp(i,k)=tet_l(i,k)+zlhocp*q_liq(i,k)         !corrected temperature
            qvap(i,k)=q_h2o(i,k)-       q_liq(i,k)         !corrected water vapor
            virt(i,k)=z1/(z1+rvd_m_o*qvap(i,k)-q_liq(i,k)) !rezipr. virtual factor 
         END DO   
      END DO
   END IF

   IF (lcalrho) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            dens(i,k)=virt(i,k)*prs(i,k)/(r_d*temp(i,k))
         END DO
      END DO
   END IF

   IF (PRESENT(r_cpd)) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            r_cpd(i,k)=z1+zrcpv*qvap(i,k)+zrcpl*q_liq(i,k)  !Cp/Cpd
         END DO
      END DO
   END IF

   IF (.NOT.ladjout) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
            tet_l(i,k)=tet_l(i,k)/exner(i,k)              !liquid water pot. temp.
         END DO
      END DO
   END IF

   IF (lcaltdv) THEN
      DO k=k_st, k_en
         DO i=i_st,i_en
!_________________________________________________________
!test: old qvsat
            pdry=(z1-qvap(i,k))*virt(i,k)*prs(i,k)         !partial pressure of dry air
! JF:             pdry=prs(i,k)
!_________________________________________________________

            qst_t(i,k)=zdqsdt( temp(i,k), zqvap( zpsat_w( temp(i,k) ), pdry ) ) 
                                                                 !d_qsat/d_T
            teta= temp(i,k)/exner(i,k)                       !potential temperature
            fakt=(zlhocp/temp(i,k)-(z1+rvd_m_o)*virt(i,k))/(z1+qst_t(i,k)*zlhocp)
            clcv=c_scld*rcld(i,k)/(z1+rcld(i,k)*(c_scld-z1)) !resulting cloud cover

            g_h2o(i,k)=grav*(rvd_m_o*virt(i,k)+clcv*fakt)    !g    -factor of q_h20-gradient
            g_tet(i,k)=grav*(z1/teta-clcv*fakt*exner(i,k)*qst_t(i,k))
                                                                 !g/tet-factor of tet_l-gradient
         END DO
      END DO
   END IF   

   !Die thermodynamischen Hilfsgroessen wurden hier
   !unter Beruecksichtigung der diagnostizierten
   !Kondensationskorrektur gebildet, indem die entspr.
   !korrigierten Werte fuer t, qv und ql benutzt wurden.
   !Beachte, dass es zu einer gewissen Inkosistenz kommt,
   !wenn die Dichte nicht angepasst wird!

END SUBROUTINE adjust_satur_equil

!********************************************************************************

SUBROUTINE solve_turb_budgets ( ktp, &
!
   i_st, i_en,                       &
   k_st, k_en,                       &
!
   imode_stke,                       &
!
   fm2, fh2,                         &
#ifdef SCLM
   it_s, grd,                        &
#endif
!---------------------------------------------------------------------------------------
   fcd, tls, tvt )

INTEGER (KIND=iintegers), INTENT(IN) :: &  !
!
  ktp,                    & !start index of vertical dimension
  i_st, i_en,             & !horizontal start- and end-indices
  k_st, k_en,             & !vertical   start- and end-indices
!
  imode_stke                !mode of solving the TKE equation

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(IN) :: &
!
  fm2,  &  !squared frequency of mechanical forcing
  fh2      !squared frequency of thermal    forcing

!-------------------------------------------------------------------------------
#ifdef SCLM
INTEGER (KIND=iintegers), INTENT(IN) :: &
!
  it_s     !iteration step

REAL (KIND=ireals), DIMENSION(:,ktp:,:), INTENT(IN) :: &
!
  grd      !vertical gradients (needed only for SUB 'turb_stat' during a SCM-run)
#endif
!-------------------------------------------------------------------------------

REAL (KIND=ireals), DIMENSION(:,kcm:), OPTIONAL, INTENT(IN) :: &
!
  fcd      !frequency of small scale canopy drag 

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(INOUT) :: &
!
  tls      !turbulent master scale

REAL (KIND=ireals), DIMENSION(:,ktp:), TARGET, INTENT(IN) :: &
!
  tvt      !total transport of turbulent velocity scale

INTEGER (KIND=iintegers) :: &
!
  i,k    !running indices for location loops

REAL (KIND=ireals) :: &
!
  val1, val2, wert, fakt, & !auxilary values
  w1,w2, &                  !weighting factors 
  q0,q1,q2,q3, &            !turbulent velocity scales 
!
  d0, d1, d2, d3, d4, d5, d6, &
  a11, a12, a22, a21, &
  a3, a5, a6, be1, be2, &
  gam0, gama, &
  sm, sh, gm, gh

REAL (KIND=ireals), DIMENSION(i_st:i_en)  :: &
! 
  l_dis, & !dissipation length scale
  l_frc, & !forcing     length scale
  frc,   & !effective TKE-forcing
  tvs      !turbulent velocity scale

LOGICAL :: corr

!------------------------------------------------------------------------------------------------------
  
  w1=tkesmot; w2=z1-tkesmot

! Stabilitaetskorrektur der turbulenten Laengenskala bei stabilier Schichtung:

  IF (a_stab.GT.z0) THEN
     DO k=k_st,k_en !von oben nach unten
        DO i=i_st, i_en
           wert=a_stab*SQRT(MAX( z0, fh2(i,k)) )
           tls(i,k)=tke(i,k,nvor)*tls(i,k)/(tke(i,k,nvor)+wert*tls(i,k))
        END DO
     END DO
  END IF     

  IF (kcm.LE.k_en .AND. PRESENT(fcd)) THEN !Es gibt einen Rauhigkeitsbestand

     DO i=i_st, i_en

!       Belegung mit Werten gueltig ausserhalb des Bestandes:

        dd(i,0)=d_m

        dd(i,1)=d_1
        dd(i,2)=d_2
        dd(i,3)=d_3
        dd(i,4)=d_4
        dd(i,5)=d_5
        dd(i,6)=d_6

        dd(i,7)=rim
     END DO

  END IF
     
  DO k=k_st, k_en

!    Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht
!    (ausser Volumenterme, die zur Diffusion gehoeren):

     IF (k.GE.kcm .AND. PRESENT(fcd)) THEN !innerhalb der Rauhigkeitsschicht
!Achtung: neue Behandlung der Rauhigkeitsschicht einfuehren

!       Reduzierte Konstanten:

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

     DO i=i_st, i_en
        l_dis(i)=tls(i,k)*dd(i,0)          !length scale of dissipation
        l_frc(i)=tls(i,k)*dd(i,4)/(z1-c_g) !length scale of forcing

        val1=tkvm(i,k)*fm2(i,k); val2=tkvh(i,k)*fh2(i,k)
!_______________________________________________________________
!test: Keine Beschraenkung von frc nach unten
        frc(i)=MAX( val1-val2, frcsecu*dd(i,7)*val1 )
! frc(i)=val1-val2
!_______________________________________________________________
     END DO
     !Bei "frcsecu=1" wird 'frc' so nach unten beschraenkt, dass die krit. Rf-Zahl (Rf=1-'dd(i,7)')
     ! im TKE-Quellgleichgewicht (ohne Transport und Zeittendenz) nicht ueberschritten wird.
     !Bei "frcsecu<1" wird 'frc' entsrechend weniger eingeschraenkt. 
     !Bei "frcsecu=0" wird Rf nicht > 1 und die Summe der TKE-Quellterme wird nicht negativ.
     !Bei "frcsecu<0" kann "Rf>1" und die Summe der TKE-Quellterme auch negativ werden.

!________________________________________________________________
!test: keine Beschraenkung von frc nach oben (Auskommentierung)
!    IF (frcsecu.GT.0) THEN
!       DO i=i_st, i_en
!          frc(i)=MIN( frc(i), frcsecu*tke(i,k,nvor)**2/l_frc(i)+(z1-frcsecu)*frc(i) )
!       END DO
!    END IF   
!________________________________________________________________
     !'tke(i,k,nvor)**2/l_frc(i)' ist der Maximal-Wert fuer 'frc', mit dem die Abweichung 
     !vom TKE-Gleichgewicht den Wert besitzt, der mit der gegebenen 'tke' bei neutraler Schichtung
     !nicht ueberschritten werden kann.
     !Bei "frcsecu<=0" wird 'frc' nicht nach oben beschraenkt.

!    Berechnung der neuen TKE-Werte:

     IF (imode_stke.EQ.1) THEN !1-st (former) type of prognostic solution
        DO i=i_st, i_en
           q0=tke(i,k,nvor)
           q1=l_dis(i)/dt_tke
           q2=MAX(z0,q0+tvt(i,k)*dt_tke)+frc(i)*dt_tke
           tvs(i)=q1*(SQRT(z1+z4*q2/q1)-z1)*z1d2
        END DO    
     ELSEIF (imode_stke.EQ.2) THEN !2-nd (new) type of prognostic solution
!print *,"tens=",k,q0,tvt(im,k),q0+tvt(im,k)*dt_tke
        DO i=i_st, i_en
           q0=tke(i,k,nvor)
           fakt=z1/(z1+z2*dt_tke*q0/l_dis(i))
           q1=fakt*(tvt(i,k)+frc(i))*dt_tke
           tvs(i)=q1+SQRT(q1**2+fakt*(q0**2+z2*dt_tke*con_m*fm2(i,k)))
        END DO    
     ELSEIF (imode_stke.EQ.-1) THEN !diagn. solution of station. TKE-equation
        DO i=i_st, i_en
           q0=tke(i,k,nvor)
           q2=l_dis(i)*(tvt(i,k)+frc(i))
           q3=l_dis(i)*con_m*fm2(i,k)
           IF (q2.GE.z0) THEN
            ! tvs(i)=SQRT(q2+q3/q0)
              tvs(i)=(q2*q0+q3)**z1d3
           ELSEIF (q2.LT.z0) THEN
              tvs(i)=SQRT(q3/(q0-q2/q0))
           ELSE
              tvs(i)=q3**z1d3
           END IF
        END DO
     ELSE !standard diagnostic solution
        DO i=i_st, i_en
           tvs(i)=SQRT(l_dis(i)*MAX( frc(i), z0 ))
        END DO
     END IF    

     DO i=i_st, i_en
        q2=SQRT(l_frc(i)*MAX( frc(i), z0 ))
!________________________________________________________________
!test: ohne tkesecu*vel_min als untere Schranke
        tke(i,k,ntur)=MAX( tkesecu*vel_min, tkesecu*q2, w1*tke(i,k,nvor)+w2*tvs(i) )
!       tke(i,k,ntur)=MAX(                  tkesecu*q2, w1*tke(i,k,nvor)+w2*tvs(i) )
!________________________________________________________________
     END DO    
     !'q2' ist ein Minimalwert fuer 'tke', mit dem die Abweichung vom TKE-Gleichgewicht
     !den Wert besitzt, der mit dem gegebenen 'frc' bei neutraler Schichtung nicht
     !ueberschritten werden kann.

!    Berechnung der neuen stabilitaetsabhangigen Laengenskalen:

     IF (lstfnct) THEN

        DO i=i_st, i_en

           d0=dd(i,0)
           d1=dd(i,1); d2=dd(i,2); d3=dd(i,3)
           d4=dd(i,4); d5=dd(i,5); d6=dd(i,6)

           w1=stbsmot; w2=z1-stbsmot

           !Loesung bei vorgegebener 'tke':

           fakt=(tls(i,k)/tke(i,k,ntur))**2

           gh=fh2(i,k)*fakt
           gm=fm2(i,k)*fakt

           IF (stbsecu.EQ.z0) THEN !Koeffizientenbelegung ohne Praeconditionierung

              be1=z1
              be2=be1-c_g
      
              a11=d1+(d5-d4)*gh
              a12=d4*gm
              a21=(d6-d4)*gh
              a22=d2+d3*gh+d4*gm

           ELSE !Koeffizientenbelegung mit Praeconditionierung:

              fakt=(tke(i,k,ntur)/tls(i,k))**2

              a11=d1*fakt+(d5-d4)*fh2(i,k)
              a12=d4*fm2(i,k)
              a21=(d6-d4)*fh2(i,k)
              a22=d2*fakt+d3*fh2(i,k)+d4*fm2(i,k)

              be1=z1/MAX( ABS(a11), ABS(a12) )
              be2=z1/MAX( ABS(a21), ABS(a22) )

              a11=be1*a11; a12=be1*a12
              a21=be2*a21; a22=be2*a22

              be1=be1*fakt; be2=be2*fakt*(z1-c_g)

           END IF 

           fakt=a11*a22-a12*a21

           sh=be1*a22-be2*a12
           sm=be2*a11-be1*a21

!          Obere Schranke fuer die Abweichung 'gama' vom TKE-Gleichgewicht:
           gam0=stbsecu/d0+(z1-stbsecu)*(z1-c_g)/d4

           IF (fakt.GT.z0 .AND. sh.GT.z0 .AND. sm.GT.z0) THEN !Loesung moeglich
              fakt=z1/fakt
              sh=sh*fakt
              sm=sm*fakt
              gama=sm*gm-sh*gh
              corr=(gama.GT.gam0)
           ELSE !keine Loesung moeglich
              gama=gam0
              corr=.TRUE.
           END IF
!_________________________________________________
!test: Korrekturbedingung wie zuvor: Bewirkt Unterschiede!
 corr=(fh2(i,k).LT.z0)
!_________________________________________________

!          Korrektur mit beschraenktem 'gama':

           IF (corr) THEN !Korrektur noetig
!             Loesung bei der die TKE-Gleichung in Gleichgewichtsform eingesetzt ist,
!             wobei 'gama' die Abweichung vom Gleichgewicht darstellt,
!             die so beschraenkt wird, dass immer eine Loesung moeglich ist:

!             gama=1/d0       : TKE-Gleichgewicht
!             gama=(z1-c_g)/d4: Maximalwert, damit immmer "be2>=0" und "be1>0"

!_________________________________________________
!test: gama mit vorangegangenem forcing (wie zuvor): Bewirkt Unterschiede!
 gama=frc(i)*tls(i,k)/tke(i,k,ntur)**2
!_________________________________________________
              gama=MIN( gam0, gama )

              be1=z1-d4*gama
              be2=be1-c_g

              a3=d3*gama/d2
              a5=d5*gama/d1
              a6=d6*gama/d2
              be1=be1/d1
              be2=be2/d2

              val1=(fm2(i,k)*be2+(a5-a3+be1)*fh2(i,k))/(z2*be1)
              val2=val1+SQRT(val1**2-(a6+be2)*fh2(i,k)*fm2(i,k)/be1)
              fakt=fh2(i,k)/(val2-fh2(i,k))

              sh=tls(i,k)*(be1-a5*fakt)
              sm=sh*(be2-a6*fakt)/(be1-(a5-a3)*fakt)
           ELSE !Es bleibt nur noch die Multipl. mit 'tls' 
              sh=tls(i,k)*sh; sm=tls(i,k)*sm
           END IF
           
!          Zeitliche Glaettung der Stabilitaetslaenge:
           tkvh(i,k)=sh*w2+tkvh(i,k)*w1
           tkvm(i,k)=sm*w2+tkvm(i,k)*w1

        END DO

     END IF

!------------------------------------------------------------------------------------
#ifdef SCLM
     IF (lsclm .AND. it_s.EQ.it_end) THEN
        CALL turb_stat(lsm=tkvm(im,k),      lsh=tkvh(im,k),        &
                       fm2=fm2(im,k),       fh2=fh2(im,k),         &
                       d_m=dd(im,0), d_h=d_h, a_m=z1/dd(im,2),     &
                       tls=tls(im,k), tvs=tke(im,k,ntur),          &
                       tvt=tvt(im,k),   grd=grd(im,k,:),    k=k)
     END IF
#endif
!SCLM--------------------------------------------------------------------------------

  END DO !k

! Sichern der TKE-Dissipation "edr=q**3/l_dis" u.a. als thermische Quelle:

  IF (PRESENT(edr).OR.ltmpcor) THEN
     IF (PRESENT(edr)) THEN
        ediss => edr
     ELSE
        ediss => tvt
     END IF
     DO k=k_st, k_en
!Achtung: edr vielleicht besser mit effektivem Wert q3*q0**2/l_dis
        DO i=i_st, i_en
           ediss(i,k)=tke(i,k,ntur)**3/(dd(i,0)*tls(i,k))
        END DO
     END DO   
  END IF   

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
    tvt,    & !total transport of turbulent velocity scale [m/s2]
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

SUBROUTINE turb_cloud ( ktp,                 &
!
   istart, iend, kstart, kend,               &
!
   icldtyp,                                  &
!
   prs, ps, t, qv, qc,                       &
!
   rcld, clcv, clwc                             )
   

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
!     icldtyp = 0 : grid scale diagnosis of saturation deficit
!     Cloud water is estimated by grid scale saturation adjustment,
!     which equals the case "icldtyp = 2" in calse of no variance.
!     So cloud cover is either "0" or "1".
!   
!     icldtyp = 1 : empirical diaagnosis of saturation deficit
!     The fractional cloud cover clcv is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.

!     icldtyp = 2: statistical diagnosis of saturation deficit
!     A Gaussion distribution is assumed for the saturation deficit
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
  ktp,               & ! start index of vertical dimension
  istart, iend,      & ! zonal      start and end index
  kstart, kend,      & ! vertical   start and end index
!
  icldtyp              ! type of cloud diagnostics

! Array arguments with intent(in):

REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: & 
  ps     ! surface pressure

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(IN) :: & 
  prs,   & ! atmospheric pressure
  t,     & ! temperature (main levels)
  qv       ! water vapour (")

REAL (KIND=ireals), DIMENSION(:,ktp:), OPTIONAL, INTENT(IN) :: & 
  qc,    & ! cloud water  (")
  rcld     ! standard deviation of saturation deficit

! Array arguments with intent(out):

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(INOUT) :: &
  clcv     ! stratiform subgrid-scale cloud cover

REAL (KIND=ireals), DIMENSION(:,ktp:), INTENT(OUT) :: &
  clwc     ! liquid water content of ""

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,k ! loop indices

REAL (KIND=ireals), PARAMETER :: &
  zsig_max = 1.0E-3_ireals,  & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_ireals,   & ! fraction of saturation specific humidity
  zuc      = 0.95_ireals       ! constant for critical relative humidity

REAL (KIND=ireals), DIMENSION(istart:iend) :: &
  qt, tl, qs, gam 

REAL (KIND=ireals) :: &
  pres, ql, dq, q, sig, uc, & !
  zsigma, zclc1, zq_max       !

!------------------------------------------------------------------
!Festlegung der Formelfunktionen fuer die Turbulenzparametrisierung:

!------------ End of header ----------------------------------------

! Begin Subroutine turb_cloud
! ---------------------------

  zq_max = q_crit*(z1/clc_diag - z1)

  DO k = kstart, kend

     IF (PRESENT(qc)) THEN
        DO i = istart, iend
           qt(i) = qv(i,k) +       qc(i,k) ! total water content
           tl(i) =  t(i,k) - lhocp*qc(i,k) ! liquid water temperature
        END DO 
     ELSE !qv and t already contain the conservation variables
        DO i = istart, iend
           qt(i) = qv(i,k)
           tl(i) =  t(i,k)
        END DO
     END IF

     DO i = istart, iend
!___________________________________________________________
!test: old qsat:
!mod_2011/09/28: zpres=patm -> zpres=pdry {
        pres=(z1-qt(i))/(z1-rvd_m_o*qt(i))*prs(i,k)        ! part. pressure of dry air
! JF:         pres=prs(i,k)
!___________________________________________________________
        qs(i) = zqvap( zpsat_w( tl(i) ), pres )              ! saturation mixing ratio
!mod_2011/09/28: zpres=patm -> zpres=pdry }
        gam(i) = z1 / ( z1 + lhocp*zdqsdt( tl(i), qs(i) ) ) ! slope factor
     END DO
!if (k.eq.2) then
! print *,"qt=",qt(im)," tl=",tl(im)," pres=",prs(im,k)
! print *,"qs=",qs(im)
!end if

     DO i = istart, iend

        pres = prs(i,k)        ! pressure
        dq   = qt(i) - qs(i) ! saturation deficit

        IF ( icldtyp .EQ. 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i)

          ! critical relative humidity
          uc = zuc - uc1 * zsigma * ( z1 - zsigma )  &
                         * ( z1 + uc2*(zsigma-0.5_ireals) )

          ! cloud cover
          clcv(i,k) = MAX( z0,  &
                            MIN( z1, clc_diag * ((qt(i)/qs(i)-uc)/(ucl-uc))**2 ) )

          ! in-cloud water content
          ql = qs(i) * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_ireals ) THEN
            zclc1 = clc_diag * ( (z1-uc)/(ucl-uc) )**2
            ql    = ql + (gam(i)*dq-ql)*(clcv(i,k)-zclc1)/(z1-zclc1)
          END IF
          clwc(i,k) = clcv(i,k) * ql

        ELSE !cloud water diagnosis based on saturation adjustment

          IF ( icldtyp .EQ. 2 ) THEN

            ! Statistical calculation of cloud cover and cloud water content
            ! using the standard deviation of the saturation deficit

            IF (PRESENT(rcld)) THEN !rcld contains standard deviation
              sig = MIN ( zsig_max, rcld(i,k) )
            ELSE !clcv contains standard deviation and will be overwritten by cloud cover
              sig = MIN ( zsig_max, clcv(i,k) )
            END IF

          ELSE !grid scale adjustment wihtout any variance
            
            sig=z0 

          END IF  

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_ireals ) THEN
            clcv(i,k) = ABS ( (SIGN(z1,dq)+z1)*0.5_ireals )
            clwc(i,k) = clcv(i,k) * gam(i) * dq
          ELSE
            q = dq/sig
            clcv(i,k) = MIN ( z1, MAX ( z0, clc_diag * ( z1 + q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              clwc(i,k) = z0
            ELSEIF ( q >= zq_max ) THEN
              clwc(i,k) = gam(i) * dq
            ELSE
              clwc(i,k) = gam(i) * sig * ( q + q_crit ) &
                                           * ( q + zq_max ) / ( z2*( q_crit + zq_max) )
!Achtung: Korrektur:
           !  clwc(i,k) = gam(i) * sig * zq_max *( (q + q_crit)/(zq_max + q_crit) )**2
            ENDIF
          ENDIF

        ENDIF

     ENDDO
!if (k.eq.2) then
! print *,"sig=",MIN ( zsig_max, rcld(im,k) )," clcv=",clcv(im,k)," clwc=",clwc(im,k)
!end if

  ENDDO

END SUBROUTINE turb_cloud

!********************************************************************************

FUNCTION zpsat_w (ztemp)

  REAL (KIND=ireals), INTENT(IN) :: ztemp
  REAL (KIND=ireals) :: zpsat_w

  zpsat_w=b1*exp(b2w*(ztemp-b3)/(ztemp-b4w))

END FUNCTION zpsat_w

FUNCTION zqvap (zpvap, zpres)

  REAL (KIND=ireals), INTENT(IN) :: zpvap, &
                                    zpres !part. pressure of dry air
  REAL (KIND=ireals) :: zqvap

!mod_2011/09/28: zpres=patm -> zpres=pdry {
!________________________________________________
!test: old qsat
! zqvap=rdv*zpvap/(zpres-o_m_rdv*zpvap)
  zqvap=rdv*zpvap
  zqvap=zqvap/(zpres+zqvap)
!________________________________________________
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zqvap

FUNCTION zdqsdt (ztemp, zqsat)

  REAL (KIND=ireals), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=ireals) :: zdqsdt

!mod_2011/09/28: zpres=patm -> zpres=pdry {
!________________________________________________
!test: old qsat
! zdqsdt=b234w*(1.0_ireals+rvd_m_o*zqsat) &
!             *zqsat/(ztemp-b4w)**2
  zdqsdt=b234w*(z1-zqsat)*zqsat/(ztemp-b4w)**2
!________________________________________________
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zdqsdt

!********************************************************************************

SUBROUTINE vert_grad_diff ( kcm, kgc,                  &
!
          i_st, i_en, k_tp, k_sf,                      &
!
          dt_var, ivtype, igrdcon, itndcon,            &
          linisetup, lnewvtype, lsflucond, lsfgrduse,  & 
!
          rho, rho_s, rho_n, hhl, r_air, tkv, dzs,     &
!
          disc_mom, expl_mom, impl_mom, invs_mom,      &
          diff_mom, scal_fac, diff_dep,                &
!
          cur_prof, dif_tend, eff_flux )

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
    lsfgrduse        ! use explicit surface gradient

                  ! DIMENSION(ie,ke) 
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
    rho              ! air density at full levels                   (Kg/m3)

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl              ! half level height                              (m)
 
!Attention: Notice the start index of vertical dimension!
                  ! DIMENSION(ie,2:ke1)
REAL (KIND=ireals), DIMENSION(:,2:), INTENT(IN) :: &
!
    tkv              ! turbulent diffusion coefficient               (m/s2)

                  ! DIMENSION(ie)
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
!
    dzs              ! effective surface layer depth                  (m)

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

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
!
!   Auxilary arrays:
!
    disc_mom, &      ! prep_inp calc_inp: discretis. momentum (rho*dz/dt) on var. levels
    expl_mom, &      ! prep_inp         : diffusion  momentum (rho*K/dz)) 
                     ! prep_out calc_inp: explicit part of diffusion momentum
    impl_mom, &      ! prep_inp         : discretis. momentum (rho*dz/dt) on flux levels
                     ! prep_out calc_inp: implicit  part of diffusion momentum
    invs_mom, &      ! prep_out calc_inp: inversion  momentum
    scal_fac, &      ! prep_out calc_inp: scaling factor due to preconditioning
    diff_dep         ! diffusion depth

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: &
!
    diff_mom         ! aux: saved comlete diffusion momentum (only in case of "itndcon.EQ.3")

                  ! DIMENSION(ie,ke1)
REAL (KIND=ireals), DIMENSION(:,:), INTENT(INOUT) :: &
!
!   Inp-out-variable:
!
    cur_prof, &      ! inp     : current   variable profile (will be overwritten!)
                     ! calc_inp: corrected variable profile
                     ! calc_out: current   variable profile including tendency increment
                     ! smot_inp: current   diffusion tendency 
    eff_flux, &      ! inp     : effective gradient
                     ! out     : effective flux density
    dif_tend         ! inp     : current time tendency of variable
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

REAL (KIND=ireals), DIMENSION(:,:), POINTER :: &
!
    rhon, rhoh

!-------------------------------------------------------------------------

!xx => disc_mom(:,2,:)

  k_lw=k_sf-1; k_hi=k_tp+1

  fr_var=z1/dt_var

  IF (ivtype.EQ.mom) THEN
     tkmin=MAX( con_m, tkmmin )
  ELSE
     tkmin=MAX( con_h, tkhmin )
  END IF

! Initial setup and adoptions for new variable type:

  IF (linisetup .OR. lnewvtype) THEN

     IF (linisetup .OR. .NOT.PRESENT(rho_n)) THEN
        DO k=k_hi,k_lw
           DO i=i_st,i_en
              expl_mom(i,k)=hhl(i,k)-hhl(i,k+1)
           END DO
        END DO
     END IF

     IF (PRESENT(rho_n)) THEN
        rhon => rho_n
     ELSE
        rhon => invs_mom ; rhoh => rho
        CALL bound_level_interp( i_st,i_en, k_hi+1,k_lw, &
                                 nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=expl_mom)

        DO i=i_st,i_en
           rhon(i,k_sf)=rho_s(i)
        END DO
     END IF 

     IF (linisetup) THEN
        DO i=i_st,i_en
           disc_mom(i,k_hi)=rho(i,k_hi)*expl_mom(i,k_hi)*fr_var
        END DO
        DO k=k_hi+1,k_lw
           DO i=i_st,i_en
              disc_mom(i,k)=rho(i,k)*expl_mom(i,k)*fr_var
              diff_dep(i,k)=z1d2*(expl_mom(i,k-1)+expl_mom(i,k))
           END DO
        END DO
     END IF   

     DO k=k_hi+1,k_lw
        DO i=i_st,i_en
           impl_mom(i,k)=rhon(i,k)*diff_dep(i,k)*fr_var
           expl_mom(i,k)=MAX( tkmin, tkv(i,k) ) &
                          *rhon(i,k)/diff_dep(i,k)
        END DO  
     END DO 
     DO i=i_st,i_en
        diff_dep(i,k_sf)=dzs(i)
        impl_mom(i,k_sf)=rhon(i,k_sf)*diff_dep(i,k_sf)*fr_var
        expl_mom(i,k_sf)=rhon(i,k_sf)*tkv(i,k_sf)/diff_dep(i,k_sf)
     END DO
     !Attention: 'tkmin' should be excluded for the surface level 'k_sf'!

     IF (itndcon.EQ.3) THEN
        DO k=k_hi+1,k_sf
           DO i=i_st,i_en
              diff_mom(i,k)=expl_mom(i,k)
           END DO  
        END DO 
     END IF    

     CALL prep_impl_vert_diff( lsflucond, i_st, i_en, k_tp=k_tp, k_sf=k_sf,             &
          disc_mom=disc_mom, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom,   &
          scal_fac=scal_fac )

  END IF

! Optional correction of vertical profiles:

  IF (lsfgrduse .AND. igrdcon.NE.2) THEN !effective surface value from effective surface gradient
     DO i=i_st,i_en
        cur_prof(i,k_sf)=cur_prof(i,k_sf-1)-diff_dep(i,k_sf)*eff_flux(i,k_sf)
     END DO
     !Note that this would be done twice in case of "igrdcon.EQ.2"!
  END IF

  IF (igrdcon.EQ.1) THEN !only correction profiles
     DO k=k_sf,kgc,-1
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-cur_prof(i,k)
        END DO
     END DO
     DO i=i_st,i_en
        cur_prof(i,kgc-1)=z0
     END DO
     DO k=kgc,k_sf
        DO i=i_st,i_en
            cur_prof(i,k)=cur_prof(i,k-1)+(cur_prof(i,k)-diff_dep(i,k)*eff_flux(i,k))
        END DO
     END DO
  ELSEIF (igrdcon.EQ.2) THEN !effektive total profile
     DO k=kgc,k_sf
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-diff_dep(i,k)*eff_flux(i,k)
        END DO
     END DO
  END IF

  IF (itndcon.EQ.3) THEN !add diffusion correction from tendency to current profile
     !Related downward flux densities:
     DO i=i_st,i_en
        eff_flux(i,2)=-dif_tend(i,1)*disc_mom(i,1)*dt_var
     END DO
     DO k=k_hi+1,k_lw
        DO i=i_st,i_en
           eff_flux(i,k+1)=eff_flux(i,k)-dif_tend(i,k)*disc_mom(i,k)*dt_var
        END DO
     END DO
     !Virtual total vertical increment:
     DO k=k_hi+1,k_sf
        DO i=i_st,i_en
           eff_flux(i,k)=eff_flux(i,k)/diff_mom(i,k)+(cur_prof(i,k-1)-cur_prof(i,k))
        END DO
     END DO   
     !Related corrected profile:
     DO k=k_hi+1,k_sf
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-eff_flux(i,k)
        END DO
     END DO    
  END IF    

  IF (itndcon.GE.1) THEN !calculate updated profile by adding tendency increment to current profile 
     DO k=k_hi,k_lw
        DO i=i_st,i_en
            dif_tend(i,k)=cur_prof(i,k)+dif_tend(i,k)*dt_var
        END DO
     END DO
     DO i=i_st,i_en
        dif_tend(i,k_sf)=cur_prof(i,k_sf)
     END DO
  END IF   

  !Final solution of the semi-implicit diffusion equation:

  CALL calc_impl_vert_diff (                                                       &
       i_st, i_en, k_tp=k_tp, k_sf=k_sf, itndcon=itndcon,                          &
       disc_mom=disc_mom, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
       scal_fac=scal_fac, cur_prof=cur_prof, upd_prof=dif_tend, eff_flux=eff_flux    )

  !Calculation of time tendencies:

  DO k=k_hi,k_lw
     DO i=i_st,i_en
        dif_tend(i,k)=(dif_tend(i,k)-cur_prof(i,k))*fr_var
     END DO
  END DO

! Volume correction within the roughness layer:

  IF (PRESENT(r_air)) THEN
     DO k=kcm,k_lw  !r_air-gradient within the roughness layer
        DO i=i_st,i_en
           cur_prof(i,k)=(r_air(i,k)-r_air(i,k+1))/(dt_var*disc_mom(i,k))
        END DO
     END DO
     DO k=kcm,k_lw  !within the roughness layer
        DO i=i_st,i_en
           dif_tend(i,k)=dif_tend(i,k)                          &
                          +z1d2*(eff_flux(i,k+1)+eff_flux(i,k)) & !flux on full level
                               *cur_prof(i,k)                       !r_air-gradient
        END DO
     END DO
  END IF    

! Optional conservative vertical smoothing of tendencies:

  IF (tndsmot.GT.z0) THEN
     DO k=k_hi,k_lw
        DO i=i_st,i_en
           cur_prof(i,k)=dif_tend(i,k)
        END DO
     END DO
     CALL vert_smooth ( &
          i_st, i_en, k_tp=k_tp, k_sf=k_sf, k_lw=k_lw, &
          disc_mom=disc_mom, cur_tend=cur_prof, smo_tend=dif_tend, versmot=tndsmot )
  END IF
    
END SUBROUTINE vert_grad_diff
         
!********************************************************************************

SUBROUTINE prep_impl_vert_diff ( lsflucond, i_st,i_en, k_tp, k_sf, &
!
   scal_fac, disc_mom, expl_mom, impl_mom, invs_mom )

!Achtung: l-Schleifen -> k-Schleifen
!Achtung: Vorzeichenwechsel fur impl. momentum ist uebersichtlicher

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

LOGICAL, INTENT(IN) :: &
!
   lsflucond    !flux condition at the surface

REAL (KIND=ireals), INTENT(IN) :: &
!
   disc_mom(:,:)    !discretis momentum (rho*dz/dt) on variable levels

REAL (KIND=ireals), INTENT(INOUT) :: &
!
   expl_mom(:,:), & !inp: diffusion momentum (rho*K/dz)
                      !out: explicit part of diffusion momentum
   impl_mom(:,:)    !inp: discretis momentum (rho*dz/dt) on flux levels
                      !out: scaled implicit part of diffusion momentum

REAL (KIND=ireals), INTENT(OUT) :: &
!
   scal_fac(:,:), & !scaling factor due to preconditioning
   invs_mom(:,:)    !inversion momentum

INTEGER (KIND=iintegers) :: &
!
   i,k, &             !horizontal and vertical coordiante indices
   m                  !level increment dependent on type of surface condition

!--------------------------------------------------------------------------------------

!  Inverse momentum vector of uppermost level:

   k=k_tp+1

   DO i=i_st, i_en
      impl_mom(i,k+1)= expl_mom(i,k+1) &
!test:
                      *MAX(MIN(expl_mom(i,k+1)/impl_mom(i,k+1), impl_s), impl_t) 
! JF:                       *MAX(impl_s-z1d2*impl_mom(i,k+1)/expl_mom(i,k+1), impl_t)
!test
      scal_fac(i,k  )= z1/SQRT(disc_mom(i,k)+impl_mom(i,k+1))
      invs_mom(i,k  )= z1
      expl_mom(i,k+1)= expl_mom(i,k+1)-impl_mom(i,k+1)
   END DO
 
!  Note: Zero flux condition just below top level.

!  Inverse momentum vector of the other levels:

   IF (lsflucond) THEN
      m=2
   ELSE
      m=1
   END IF

   DO k=k_tp+2, k_sf-m

      DO i=i_st, i_en
         impl_mom(i,k+1)= expl_mom(i,k+1) &
!test:
                         *MAX(MIN(expl_mom(i,k+1)/impl_mom(i,k+1), impl_s), impl_t)
! JF:                          *MAX(impl_s-z1d2*impl_mom(i,k+1)/expl_mom(i,k+1), impl_t)
!test
         scal_fac(i,k  )= z1/SQRT(disc_mom(i,k  )+impl_mom(i,k)+impl_mom(i,k+1))
         impl_mom(i,k  )= scal_fac(i,k-1)*scal_fac(i,k  )*impl_mom(i,k  )
         invs_mom(i,k  )= z1/(z1-invs_mom(i,k-1)*impl_mom(i,k  )**2)
         expl_mom(i,k+1)= expl_mom(i,k+1)-impl_mom(i,k+1)
      END DO

   END DO

   IF (lsflucond) THEN

      k=k_sf-1

      DO i=i_st, i_en
         impl_mom(i,k+1)= z0 !no  implicit surface weight
         scal_fac(i,k  )= z1/SQRT(disc_mom(i,k  )+impl_mom(i,k))
         impl_mom(i,k  )= scal_fac(i,k-1)*scal_fac(i,k  )*impl_mom(i,k  )
         invs_mom(i,k  )= z1/(z1-invs_mom(i,k-1)*impl_mom(i,k  )**2)
      END DO
   END IF
!read *,xx
!print *,"ende prep"

END SUBROUTINE prep_impl_vert_diff

!********************************************************************************

SUBROUTINE calc_impl_vert_diff ( i_st,i_en, k_tp, k_sf, itndcon, &
!
   scal_fac, disc_mom, expl_mom, impl_mom, invs_mom, cur_prof, upd_prof, eff_flux)

!Achtung: l-Schleifen -> k-Schleifen
!Ev. nur fuer gewuenschte levels

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_tp,k_sf, & !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1
!
   itndcon      ! mode index of current tendency  consideration
                ! 0: no consider., 1: consider curr. tend. in implicit equation only
                !                  2: add curr. tend. increment to current profile
                !                  3: add related diffusion correction to current profile

REAL (KIND=ireals), INTENT(IN) :: &
!
   expl_mom(:,:), & !explicit part of diffusion momentum (rho*K/dz)
   impl_mom(:,:), & !scaled implicit part of diffusion momentum
   disc_mom(:,:), & !discretis momentum (rho*dz/dt)
   scal_fac(:,:), & !scaling factor due to preconditioning
!
   invs_mom(:,:)    !inversion momentum

REAL (KIND=ireals), TARGET, INTENT(INOUT) :: &
!
   cur_prof(:,:), & !inp: current vertical variable profile (including gradient corrections)
                      !out: current vertical variable profile (including tendency increment) 
   upd_prof(:,:)    !inp: updated vertical variable profile (including tendency increment)
                      !out: updated vertical variable profile by vertical diffusion
   !'cur_prof' as out equals 'upd_prof' as inp

REAL (KIND=ireals), INTENT(OUT) :: &
!
   eff_flux(:,:)    !effective flux density belonging to (semi-)implicit diffusion (positiv downward)
                      !aux: difference between adjacent profile values of adjacent levels
                      !aux: explicit flux density (positiv upward)
                      !aux: interim solution of diffusion equation
                      !aux: variable differece representing pure diffusion tendency

INTEGER (KIND=iintegers) :: &
!
   i,k              !horizontal and vertical coordiante indices

REAL (KIND=ireals) :: &
!
   res_flux           !resulting right hand side flux density
 
REAL (KIND=ireals), POINTER :: &
!
   old_prof(:,:)      !old variable profile value

!-----------------------------------------------------------------------

!  Preparation:

   IF (itndcon.EQ.0 .OR. itndcon.EQ.3) THEN !old profile = current profile
      old_prof => cur_prof
   ELSE !old profile = updated profile
      old_prof => upd_prof
   END IF

   IF (itndcon.EQ.2) THEN !difference form updated profile
      DO k=k_tp+2, k_sf
         DO i=i_st, i_en
            eff_flux(i,k)=upd_prof(i,k  ) - upd_prof(i,k-1)
         END DO
      END DO
   ELSE !difference from current profile  
      DO k=k_tp+2, k_sf
         DO i=i_st, i_en
            eff_flux(i,k)=cur_prof(i,k  ) - cur_prof(i,k-1)
         END DO
      END DO
   END IF

!  Uppermost level:

   k=k_tp+1

   DO i=i_st, i_en
      eff_flux(i,k+1)=  expl_mom(i,k+1) * eff_flux(i,k+1)
      res_flux         =  disc_mom(i,k  ) * old_prof(i,k  ) &
                                            + eff_flux(i,k+1)
      eff_flux(i,k)  =  scal_fac(i,k  ) * res_flux &
                        * invs_mom(i,k  )            
   END DO
!print *,"kk=",k," disc_mom=",disc_mom(im,k)
!print *," old_prf=",cur_prof(im,k)

!  Note: Zero flux condition just below top level.

!  Intermediate levels:
 
   DO k=k_tp+2, k_sf-2

      DO i=i_st, i_en
         eff_flux(i,k+1)=  expl_mom(i,k+1) * eff_flux(i,k+1)
         res_flux         =  disc_mom(i,k  ) * old_prof(i,k  ) &
                           - eff_flux(i,k  ) + eff_flux(i,k+1)
         eff_flux(i,k)  = (scal_fac(i,k  ) * res_flux          + impl_mom(i,k  ) * eff_flux(i,k-1)) &
                           * invs_mom(i,k  )
      END DO
      !Note:
      !'eff_flux' is an explicit flux density on the right hand side of 'res_flux' 
      !           and a preparation value towards the new profile there after.
!print *,"kk=",k," disc_mom=",disc_mom(im,k)," impl_mom=",impl_mom(im,k)
!print *," old_prf=",cur_prof(im,k)
   END DO

!  Surface level:

   k=k_sf-1

   DO i=i_st, i_en
      eff_flux(i,k+1)=  expl_mom(i,k+1) * eff_flux(i,k+1)
      res_flux         =  disc_mom(i,k  ) * old_prof(i,k  ) &
                        - eff_flux(i,k  ) + eff_flux(i,k+1) &
                        + impl_mom(i,k+1) * cur_prof(i,k+1)
      eff_flux(i,k)  = (scal_fac(i,k  ) * res_flux          + impl_mom(i,k  ) * eff_flux(i,k-1)) &
                        * invs_mom(i,k)
   END DO
      
!print *,"kk=",k," disc_mom=",disc_mom(im,k)," impl_mom=",impl_mom(im,k)
!print *," old_prf=",cur_prof(im,k)
!print *,"kk=",k," impl_mom=",impl_mom(im,k+1)
!print *," old_prf=",cur_prof(im,k+1)
!read *,xx

!  Note: Explicit flux condition at surface level in case of "impl_mom(:,k+1)=0".

!  Backsubstitution:

   IF (itndcon.GT.0) THEN
      DO k=k_sf-1, k_tp+2, -1
         DO i=i_st, i_en
            cur_prof(i,k)=upd_prof(i,k) 
         END DO  
      END DO   
   END IF

   DO k=k_sf-1, k_tp+2, -1
      DO i=i_st, i_en
         eff_flux(i,k-1)=  eff_flux(i,k-1) + impl_mom(i,k  ) * eff_flux(i,k  )  &
                                                                   * invs_mom(i,k-1)
         upd_prof(i,k  )=  scal_fac(i,k  ) * eff_flux(i,k  )
         eff_flux(i,k+1)=  cur_prof(i,k  ) - upd_prof(i,k  )
      END DO
      !Note:
      !'eff_flux(i,k-1)' is the new profile value except the final rescaling
      !'cur_prof(i,k  )' is now the current profile including optional explicit tendencies
      !'eff_flux(i,k+1)' is the difference between the profile values before and after diffusion
      !'upd_prof(i,k  )' contains the new profile value after diffusion
!print *,"k+1=",k+1," cur=",cur_prof(im)," new=",upd_prof(im,k  )
!print *," eff_flux=",eff_flux(im,k+1)
   END DO

   k=k_tp+1

   IF (itndcon.GT.0) THEN
      DO i=i_st, i_en
         cur_prof(i,k)=upd_prof(i,k) 
      END DO  
   END IF

   DO i=i_st, i_en
      upd_prof(i,k  )=  scal_fac(i,k  ) * eff_flux(i,k  )
      eff_flux(i,k+1)=  cur_prof(i,k  ) - upd_prof(i,k  )

      eff_flux(i,k  )=  z0 !upper zero-flux condition
   END DO
!print *,"k+1=",k+1," cur=",cur_prof(im)," new=",upd_prof(im,k  )
!print *," eff_flux=",eff_flux(im,k+1)
!print *

!  Effective flux density by vertical integration of diffusion tendencies:

   DO k=k_tp+2, k_sf
!print *,"k=",k,"oben=",eff_flux(im,k-1)," unten=",eff_flux(im,k  )
      DO i=i_st, i_en
         eff_flux(i,k  ) = eff_flux(i,k-1) + eff_flux(i,k  )*disc_mom(i,k-1) 
      END DO
!print *," eff_flux=",eff_flux(im,k)
   END DO
   !Note:
   !'eff_flux' is the vertical flux density of pure diffusion now (positive downward)

END SUBROUTINE calc_impl_vert_diff

!********************************************************************************

SUBROUTINE vert_smooth( i_st, i_en, k_tp, k_sf, k_lw, &
                        cur_tend, disc_mom, smo_tend, versmot )

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_tp,k_sf, & !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1
   k_lw         !vertical index of the lowest level

REAL (KIND=ireals), INTENT(IN) :: &
!
   cur_tend(:,:), & !current vertical tendency profile
   disc_mom(:,:), & !discretis momentum (rho*dz/dt) on variable levels

   versmot            !vertical smoothing factor

REAL (KIND=ireals), INTENT(OUT) :: &
!
   smo_tend(:,:)    !smoothed vertical tendency profile

INTEGER (KIND=iintegers) :: &
!
   i,k              !horizontal and vertical coordiante indices

REAL (KIND=ireals) :: &
!
   remfact          !remaining factor of current level

!----------------------------------------------------------------------------

   remfact=z1-versmot

   k=k_tp+1

   DO i=i_st,i_en
      smo_tend(i,k)=remfact* cur_tend(i,k)                     &
                     +versmot* cur_tend(i,k+1)*disc_mom(i,k+1) &
                                                /disc_mom(i,k)
   END DO

   remfact=z1-z2*versmot

   DO k=k_tp+2, k_sf-2

      DO i=i_st,i_en
         smo_tend(i,k)=remfact* cur_tend(i,k)                      &
                        +versmot*(cur_tend(i,k-1)*disc_mom(i,k-1)  &
                                 +cur_tend(i,k+1)*disc_mom(i,k+1)) &
                                                   /disc_mom(i,k)
      END DO

   END DO

   remfact=z1-versmot

   k=k_sf-1   

   DO i=i_st,i_en
      smo_tend(i,k)=remfact* cur_tend(i,k)                     &
                     +versmot* cur_tend(i,k-1)*disc_mom(i,k-1) &
                                                /disc_mom(i,k)
   END DO

   DO k=k_sf, k_lw-1

      DO i=i_st,i_en
         smo_tend(i,k)=cur_tend(i,k)
      END DO
 
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

REAL (KIND=ireals), DIMENSION(:,:), POINTER :: &
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
            DO i=i_st, i_en
               usdep(i,k)=depth(i,k)-depth(i,k+1)
            END DO
         END DO
      ELSE
         usdep => depth
      END IF   
      IF (lrpdep) THEN !precalculation of the reciprocal layer depth
         DO k=k_en, k_st, -1
            DO i=i_st, i_en
               rpdep(i,k)=z1/(usdep(i,k-1)+usdep(i,k))
            END DO
         END DO
         DO n=1, nvars
            DO k=k_en, k_st, -1
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
            DO i=i_st, i_en
               pvar(n)%bl(i,k)=z2*pvar(n)%ml(i,k)-pvar(n)%ml(i,k+1)
            END DO
         END DO
      END DO
   END IF   

END SUBROUTINE bound_level_interp

!********************************************************************************
!********************************************************************************

END MODULE src_turbdiff_new
