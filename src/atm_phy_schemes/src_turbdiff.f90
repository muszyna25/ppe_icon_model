!+ Source module for computing diffusion and transfer coefficients 
!-------------------------------------------------------------------------------

 MODULE src_turbdiff  

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
!  'turb_param.incf' contains only statement functions and is renamed into 'stat_funct.incf'
!  Moving content of 'stab_funct.incf' into the SUBROUTINE stab_funct.
! V3_24        2007/04/26 Ulrich Schaettler
!  Introduced call to exchg_boundaries for wind-tendencies again;
!  it is necessary for imode_turb=2/3!
! V4_3         2008/02/25 Matthias Raschendorfer
!  Introduction of a 3D diagnostic field 'edr' for the eddy dissipotion rate.
!  Changing the treatment of parameter field dd(:,:,:) to achieve better vectorization
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
!   by introducing the arrays 'dz0(:,:,mom|sca)' and 'vh0(:,:)
!   removing wrong surface conditions in case "imode_turb.GE.3 .AND. itype_tran.EQ.2".
!              2011/09/28 Matthias Raschendorfer
!  Correction in formula for specific humidity at saturation by using partial pressure of dry air.
!              2011/12/08 Matthias Raschendorfer
!  Introduction of SUBs 'prep_impl_vert_diff' and 'calc_impl_vert_diff' substituting the code segment
!   for implicit vertical diffusion calculation for TKE or all other variables, if 'imode_turb' is one of 3 or 4.
!   These SUBs allow the calculation of semi implicit diffusion for arbitrary variables defined on full- or  half-levels.
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
    a_stab,       & ! factor for stability correction of turbulent length scale
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
    tkesmot,      & ! time smoothing factor for TKE and diffusion coefficients
    wichfakt,     & ! vertical smoothing factor for explicit diffusion tendencies
    securi,       & ! security factor for maximal diffusion coefficients
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
    imode_tran,   & ! mode of surface-atmosphere transfer
    icldm_tran,   & ! mode of cloud representation in transfer parametr.
!
    imode_turb,   & ! mode of turbulent diffusion parametrization
    icldm_turb,   & ! mode of cloud representation in turbulence parametr.
    itype_sher,   & ! type of shear production for TKE
!
    ltkesso,      & ! calculation SSO-wake turbulence production for TKE
    ltkecon,      & ! consider convective buoyancy production for TKE
    lexpcor,      & ! explicit corrections of the implicit calculated
!                   ! turbulent diffusion (only if itype_turb=3)
    ltmpcor,      & ! consideration of thermal TKE-sources in the enthalpy budget
    lprfcor,      & ! using the profile values of the lowest main level instead of
!                   ! the mean value of the lowest layer for surface flux calulations
    lnonloc,      & ! nonlocal calculation of vertical gradients used
!                   ! for turbulent diffusion (only if itype_turb=3)
    lcpfluc,      & ! consideration of fluctuations of the heat capacity of air
    limpltkediff, & ! use semi-implicit TKE diffusion
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
    lsclm, i_cal, im, jm, &
!
    UUA, VVA, UWA, VWA, WWA, UST, TWA, QWA, TTA, TQA, QQA, SHF, LHF, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
#endif
!SCLM---------------------------------------------------------------------------

!-------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE
PUBLIC  :: init_canopy, turbtran, turbdiff, turb_cloud

REAL (KIND=ireals), PARAMETER :: &
    z0=0.0_ireals,&
    z1=1.0_ireals,&
    z2=2.0_ireals,&
    z3=3.0_ireals,&
    z4=4.0_ireals,&
    z5=5.0_ireals,&
    z6=6.0_ireals,&
    z7=7.0_ireals,&
    z8=8.0_ireals,&
    z9=9.0_ireals,&
    z10=10.0_ireals

REAL (KIND=ireals) :: &
    z1d2=z1/z2,&
    z1d3=z1/z3,&  
!   z2d3=z2/z3,&
    z3d2=z3/z2

INTEGER (KIND=iintegers) :: &
    istat=0, ilocstat=0

LOGICAL :: lerror=.FALSE.

INTEGER (KIND=iintegers), PARAMETER :: &
!
!    Indexgrenzen:
!
     nscal=3,     & !Skalare Groessen
     ninv=2,      & !daraus abgeleitete gegenueber vertikalen
                    !(feuchtadiabatischen) Verrueckungen
                    !invarianten Groessen
     nvel=2,      & !Geschwindigkeitskomponenten
     ndiff=nscal+nvel, &
     nred=ninv+nvel,   &
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
     u_s=u_m,     & !zonale Geschw.komp. im gestag. Gitter
     v_s=v_m,     & !meridionale  ,,     ,,   ,,       ,,
     tet=tet_l,   & !pot.Temperatur
     tem=tet,     & !Temperatur
     vap=h2o_g,   & !Wasserdampfmischungsverh.
!
     mom=1,       & !momentum
     sca=2          !scalar

!    Beachte:u_m,v_m, sowie u_s,v_s muessen in [1,nvel] liegen;
!           aber: tet_l,h2o_g in [nvel+1,nred]
!           und   tem (tet),vap,liq in [nvel+1,ndiff]


INTEGER (KIND=iintegers) :: &
     nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes
     it_start      !Startindex der Iterationen

LOGICAL :: lini  !initialization of required

REAL (KIND=ireals), TARGET :: &
     c_tke,tet_g,c_g,rim, &
     d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
     a_3,a_5,a_6,b_1,b_2, &
     l_scal


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
           istartpar, iendpar, &
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
    istartpar, iendpar ! zonal start and end index including the model boundary lines

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
    eai,          & ! (evaporative) earth area index                ( 1 )
!
    h_can,        & ! height of the vertically resolved canopy
    d_pat           ! horizontal pattern length scale

#ifdef __xlC__
REAL (KIND=ireals), DIMENSION(ie,kcm:ke1), OPTIONAL, INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,kcm:), OPTIONAL, INTENT(INOUT) :: &
#endif
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml,        & ! effective drag coefficient of canopy elements
!                   ! smaller than the turbulent length scale            (1/m)
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)

! ----------------
! Local variables:
! ----------------

  INTEGER (KIND=iintegers) ::  &
    i,       & !  loop indices
    kcp        !  buffer for the vertical index of the upper boudary of the canopy

  REAL    (KIND=ireals   ) ::  fakt

!-------------------------------------------------------------------------------
! Begin Subroutine init_canopy
!-------------------------------------------------------------------------------

  kcp=kcm !save current value of 'kcm', that might have been used for allocation before

  IF (PRESENT(h_can) .AND. PRESENT(hhl)) THEN
   ! h_can is a primary external parameter. The initial values of h_can are 0.
   ! If we don't change this, no canopy will be resolved in the vertical direction.

     kcm=ke
     DO WHILE (MAXVAL( h_can - hhl(:,kcm) + hhl(:,ke1) ) .GT. z0)
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
     DO i=istartpar, iendpar
        IF (fr_land(i) <= z1d2) THEN
           d_pat(i)=z0
        ELSE
           d_pat(i)=z1 !should be a 2D external parameter field
        END IF
     END DO
  END IF
! Effective values of the surface area indices:
  IF (PRESENT(sai) .AND. PRESENT(eai) .AND. PRESENT(tai) .AND. &
      PRESENT(lai) .AND. PRESENT(plcov)) THEN

     DO i=istartpar, iendpar
        IF (fr_land(i) <= z1d2) THEN
           sai(i)=c_sea
        ELSE
           tai(i)=MAX( 1.0E-6_ireals, lai(i) )
        END IF
     END DO

     IF (itype_tran.EQ.1) THEN
        DO i=istartpar, iendpar
           IF (fr_land(i) > z1d2) THEN
              sai(i)=tai(i)
              eai(i)=(z1-plcov(i))*sai(i)
              tai(i)=plcov(i)*tai(i)
            END IF
        END DO
     ELSE
        DO i=istartpar, iendpar
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


SUBROUTINE turb_param (l_hori)

   IMPLICIT NONE

   REAL (KIND=ireals), INTENT(IN) :: l_hori

   INTEGER (KIND=iintegers) :: i !loop indices

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

      rim=z1/(z1+(d_m-d_4)/d_5)

      a_3=d_3/(d_2*d_m)
      a_5=d_5/(d_1*d_m)
      a_6=d_6/(d_2*d_m)

      b_1=(z1-d_4/d_m)/d_1
      b_2=(z1-d_4/d_m-c_g)/d_2

      tet_g=grav/cp_d
      l_scal=MIN( l_hori, tur_len )


END SUBROUTINE turb_param

!********************************************************************************

!+ Module procedure turbtran in "src_turbdiff" for computing the coefficients
!+ for turbulent transfer


SUBROUTINE turbtran(iini, dt_tke, nprv, ntur, ntim, &
!
          ie, ke, ke1, &
          istart, iend, istartpar, iendpar, &
          ntstep, &
!    
          l_hori, hhl, &
!    
          fr_land, depth_lk, sai, &
!
          h_ice, ps, t_g, qv_s, &
          u, v, w, t, qv, qc, prs, &
!    
          gz0, tcm, tch, tfm, tfh, tfv, &
          tke, tkvm, tkvh, rcld, edr, &
!    
          t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m, &
!
          ierrstat, errormsg, eroutine)


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
!  Changing the treatment of parameter field dd(:,:,:) to achiev better vectorisation
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

!     Uebergabevariablen:

! 0. Parameters controlling the call of 'turbtran':


REAL (KIND=ireals), INTENT(IN) :: & 
   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
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

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------
!   These variables give the start- and the end-indices of the 
!   forecast for the prognostic variables in a horizontal layer.
!   Note, that the indices for the wind-speeds u and v differ from 
!   the other ones because of the use of the staggered Arakawa-B-grid.
!    
! Zonal direction:
    istart,    iend,    & ! start and end index for the inner model domain
    istartpar, iendpar    ! start and end index including model boundary lines

! Time step indices:
! -----------------------------------------------------------------------

INTEGER (KIND=iintegers), OPTIONAL, INTENT(IN) :: &
!
    ntstep          ! current time step

REAL (KIND=ireals), INTENT(IN) :: &
!
! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------
!
    l_hori          ! horizontal grid spacing (m)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
#endif
!
    hhl             ! height of model half levels                   ( m )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
#endif
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area             ( 1 )
    depth_lk,     & ! lake depth                                    ( m )
    sai,          &  ! surface area index                            ( 1 )
!
! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------
!
    h_ice,        & ! ice thickness                                 (  m  )
    ps,           & ! surface pressure                              ( pa  )
    t_g,          & ! weighted surface temperature                  (  k  )
    qv_s            ! specific water vapor content on the surface   (kg/kg)


#ifdef  __xlC__
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

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), TARGET, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(IN) :: &
#endif
!
     prs            ! atmospheric pressure                          ( pa  )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
#endif
!
     w              ! vertical wind speed (defined on half levels)  ( m/s )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(INOUT) :: &
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
     tfm,          & ! factor of laminar transfer of momentum          --
     tfh,          & ! factor of laminar transfer of scalars           --
     tfv             ! laminar reduction factor for evaporation        --
 
! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=ireals), DIMENSION(ie,ke1,ntim), INTENT(INOUT) :: &
!
     tke             ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                     ! (defined on half levels)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,2:ke1), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,2:), INTENT(INOUT) :: &
#endif
!
     tkvm,         & ! turbulent diffusion coefficients for momentum (m/s2 )
     tkvh            ! turbulent diffusion coefficients for heat     (m/s2 )
                     ! and moisture

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(INOUT) :: &
#endif
!
     rcld            ! standard deviation of the saturation deficit        
                     ! (as input and output)                             
                     ! fractional cloud cover (in turbdiff)            --

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), OPTIONAL, INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: &
#endif
!
     edr             ! eddy dissipation rate of TKE (EDR)            (m2/s3)

#ifdef  __xlC__
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

INTEGER (KIND=iintegers), INTENT(INOUT) :: ierrstat

CHARACTER (LEN=*), INTENT(INOUT) :: eroutine
CHARACTER (LEN=*), INTENT(INOUT) :: errormsg

! Local parameters:

REAL (KIND=ireals), DIMENSION(ie,0:7), TARGET :: &
     dd        !local derived turbulence parameter

! Local scalars:

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &
!
              i_st,i_en,    & !Start- und Endindices in zonale Richtung
              i,k,          & !Diskretis.index fuer lamda, phi, sigma
              it_durch        !Durchgangsindex der Iterationen

!     Lokale real Variablen:

      REAL (KIND=ireals) :: &
!
!          Hilfsvariablen:
!
           wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,             & !  ,,         ,,     ,,      Faktoren
!
!     Platzh. fuer Druck und Ableitung des Wasserdampf-Mischungsverhaeltnisses bei Saettigung
!     nach der Temperatur:
!
           patm, qsat_dT, &
!
!     Platzh. fuer relat. Wolkenanteil, virtuellen Faktor und Exner-Faktor:
!
           virt,exnr, &
!
!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 (fh2,fm2):
!
           fh2,fm2, &
!
!     Platzh. fuer horiz. Geschw.-Komponenten und Geschw.betrag:
!
           vel1,vel2,velh, &
!
!     Platzh. fuer SQRT(2*TKE)-Werte:
!
           q0,q1,q2,q3, &
!
!     Platzh. fuer Hoehendifferenzen und Laengaenskalen:
!   
           dh,l_turb,lh,lm,z0d,z_surf,len1,len2, &
!amk
           z_surf_h, &
!xxx
           dz_sg_m, dz_sg_h, dz_g0_m, dz_g0_h, &
           dz_s0_m, dz_sa_m, &
           h_2m, h_10m, a_top, a_atm, a_2m, a_10m, &
!
!     sonstiges:
!
           a_tet,a_vap,xf,wf, w1, w2, &
           uk1,vk1,uk,vk, &
           rat_m,rat_h,fac_m,fac_h, fr_sd_h

! Local arrays:

      LOGICAL            ::  &
!
           lo_ice(ie)   ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=ireals) ::  &
!
        l_tur_z0 (ie),    &
!
!US for a better vectorization (hopefully):
!
        h_top_2d (ie),    &
        h_atm_2d (ie),    &
!
        z0m_2d   (ie),    &
!amk
        z0h_2d   (ie),    &
!xxx
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
        vel1_2d  (ie),    &
        vel2_2d  (ie),    &
        vel_2d   (ie),    &
        ta_2d    (ie),    &
        qda_2d   (ie),    &
        ts_2d    (ie),    &
        qds_2d   (ie),    &
        rho_2d   (ie),    &
        prs_2d   (ie),    &
!
        dz_0a_m  (ie),    &
        dz_0a_h  (ie),    &
        dz_sa_h  (ie),    &
        dz_s0_h  (ie)

     INTEGER (KIND=iintegers) ::  &
        k_2d     (ie)

     REAL (KIND=ireals) :: &
!
!     Vertikale Gradienten verschiedener thermodyn. Variablen:
!
           grad(ie,nred), &                          
!
!     Zwischenspeicher fuer Wolkenwasser und Bedeckungsgrad:
!
           clc(ie,ke:ke),clcw(ie,ke:ke)

!Declaration of statement functions:

REAL (KIND=ireals) :: &
     zexner, zpres               !Exner factor and its argument

!Exner-factor:
 zexner(zpres)=(zpres/p0ref)**rdocp

!---- End of header ---------------------------------------------------------

! 1)  Vorbereitungen:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbtran'; lerror=.FALSE.


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
 END IF

 nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

!     Unterste Hauptflaeche halbiert den Abstand zur
!     untersten Nebenflaeche:


      xf=z2 
      wf=xf-z1
      xf=z1/xf

!     Festlegung der horiz. Schleifengrenzen: 

      i_st=istartpar 
      i_en=iendpar

!     Berechnung abgeleiteter Parameter:

      CALL turb_param  (l_hori)

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

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren: 

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

      DO i=i_st,i_en

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

         DO i=i_st,i_en

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
                  l_turb=akt*MAX( len_min, l_turb/(z1+l_turb/l_scal) )
!test
!                 l_turb=akt*MIN( l_scal, hhl(i,ke)-hhl(i,ke1) )
!test

                  dh=0.5_ireals*(hhl(i,ke-1)-hhl(i,ke1))

                  grad(i,u_m)=(u(i,ke-1)-u(i,ke))/dh
                  grad(i,v_m)=(v(i,ke-1)-v(i,ke))/dh

                  grad(i,tet_l)=(t(i,ke-1)-t(i,ke))/dh + tet_g

                  fm2=grad(i,u_m)**2+grad(i,v_m)**2
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

                  wert=lm*q0*SQRT(fm2)
                  gz0(i)=MAX( grav*len_min, alpha0*wert )
!                 gz0(i)=MAX( grav*len_min, &
!                               alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF

            tkvm(i,ke)=con_m
            tkvh(i,ke)=con_h
            tkvm(i,ke1)=con_m
            tkvh(i,ke1)=con_h

         END DO
      END IF


! 3)  Diagnose des Fluessigwassergehaltes und des Bedeckungsgrades
!     in der unterstersten Modellschicht:

      CALL turb_cloud (ie, ke,ke1, ke,ke, i_st,i_en, ke,ke, &
           prs, ps, rcld, t, qv, qc, clc, clcw              )

      IF (icldm_tran.EQ.0) THEN
         DO i=i_st,i_en
            clc(i,ke)=z0
         END DO 
       ELSEIF (icldm_tran.EQ.1) THEN
         DO i=i_st,i_en
            IF ( qc(i,ke) > z0 ) THEN
               clc(i,ke)=z1
            ELSE
               clc(i,ke)=z0
            ENDIF
         END DO   
      ENDIF

! 4)  Berechnung der Transferkoeffizienten:

      DO it_durch=it_start, it_end !Iterationen

         DO i=i_st,i_en

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

            z0h_2d(i) = z0m_2d(i)
            z_surf_h    = z0h_2d(i)/sai(i)

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
!xmk
            dz_sg_h=fakt*rlam_heat*z_surf_h*(rat_h/rat_m)
!xxx

!           Berechnung weiterer effektiver Widerstandslaengen fuer Skalare:

!           Bestandesschicht ohne lam. Grenzschicht:
!xmk
            dz_g0_h=z_surf_h*LOG(rat_m)
!xxx

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

!xmk
            rat_h=(tkvh(i,ke)*z0h_2d(i))/(tkvh(i,ke1)*a_top)
            rat_h_2d(i)=MIN( z2, MAX( z1d2, rat_h ) )

            fac_h_2d(i)=(rat_h_2d(i)-z1)*z0h_2d(i)/h_top_2d(i)
            IF (fac_h_2d(i).EQ.z1) THEN
               dz_0a_h(i)=z0h_2d(i)*h_atm_2d(i)/a_atm
            ELSE
               dz_0a_h(i)=z0h_2d(i)*LOG(a_atm/(z0h_2d(i)+fac_h_2d(i)*h_atm_2d(i))) &
                                      /(z1-fac_h_2d(i))
            END IF
!xxx

!           von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):
            dz_sa_m      = dz_s0_m      + dz_0a_m(i)
            dz_sa_h(i) = dz_s0_h(i) + dz_0a_h(i)

!           Reduktionsfaktoren fuer die Bestandesschicht
!           incl. lam. Grenzschicht:

            tfm(i)=dz_0a_m(i)/dz_sa_m      
            tfh(i)=dz_0a_h(i)/dz_sa_h(i)

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den
!           Faktor 'rat_lam' gegenueber fuehlbarer Waerme vergroesserten
!           laminaren Transpostwiderstandes:

            tfv(i)=z1/(z1+(rat_lam-z1)*dz_sg_h/dz_sa_h(i))
         END DO

         DO i=i_st,i_en
!           Berechnung der Windgeschwindigkeiten an den Massepunkten:

            vel1_2d(i)=u(i,ke)
            vel2_2d(i)=v(i,ke)

!           Berechnung der korrigierten Prandtl-Schicht-Werte:

            IF (lprfcor) THEN
               len1=z2*h_top_2d(i)
               len2=(h_top_2d(i)-h_atm_2d(i))**2 &
                 /((hhl(i,ke-1)+hhl(i,ke))*z1d2-hhl(i,ke1)-h_atm_2d(i))
               lm=len1-tfm(i)*h_atm_2d(i)-len2
               lh=len1-tfh(i)*h_atm_2d(i)-len2

               velh=u(i,ke-1)
               vel1_2d(i)=(len1*vel1_2d(i)-len2*velh)/lm
               velh=v(i,ke-1)
               vel2_2d(i)=(len1*vel2_2d(i)-len2*velh)/lm

               ta_2d(i)=(len1*t(i,ke)-h_atm_2d(i)*tfh(i)*t_g(i) &
                                    -len2*t(i,ke-1))/lh
               qda_2d(i)=(len1*qv(i,ke)-h_atm_2d(i)*tfh(i)*qv_s(i) &
                                      -len2*qv(i,ke-1))/lh

               vel_2d(i)=MAX( vel_min, SQRT(vel1_2d(i)**2+vel2_2d(i)**2) )
            ELSE
               ta_2d(i)=t(i,ke)
               qda_2d(i)=qv(i,ke)

               vel_2d(i)=MAX( vel_min, SQRT(vel1_2d(i)**2+vel2_2d(i)**2) )
            END IF

!           Berechnung der unteren Randwerte der Prandtl-Schicht:
!Achtung!
             ts_2d(i)= ta_2d(i)*(z1-tfh(i))+ t_g(i)*tfh(i)
            qds_2d(i)=qda_2d(i)*(z1-tfh(i))+qv_s(i)*tfh(i)
         END DO

         IF (it_durch.LE.it_end) THEN
            !Strange to say this seemingly needless condition makes the NEC-compiler
            !to inline the CALLs of thermodynamic FUNCTIONs and finally to do the
            !vectorization of the follwing loop!

         DO i=i_st,i_en
!           Berechnung von thermodynamischen Groessen:

            exnr   =zexner(ps(i)) !Exner-Faktor
            virt   =z1/(z1+rvd_m_o*qds_2d(i)-clcw(i,ke)) 
                                    !rezip. virtueller Faktor
!mod_2011/09/28: zpres=patm -> zpres=pdry {
            patm   =(z1-qds_2d(i))*virt*ps(i)
         
            qsat_dT=zdqsdt( ts_2d(i), zqvap( zpsat_w( ts_2d(i) ), patm ) )
!mod_2011/09/28: zpres=patm -> zpres=pdry }

            rho_2d(i)=virt*ps(i)/(r_d*ts_2d(i)) !Luftdichte

!           Berechn. der thermodyn. Hilfsfaktoren:
            fakt=lhocp/ts_2d(i)-(z1+rvd_m_o)*virt
            wert=clc(i,ke)*fakt/(z1+qsat_dT*lhocp)

!           Berechn. der thermodyn. Faktoren A_tet und A_q:
            a_tet=z1/ts_2d(i)-wert*qsat_dT
            a_vap=rvd_m_o*virt+wert

!           Berechnung der benoetigten vertikalen Gradienten und
!           abgeleiteter Groessen:

            grad(i,u_m)=tfm(i)*vel1_2d(i)/dz_0a_m(i)
            grad(i,v_m)=tfm(i)*vel2_2d(i)/dz_0a_m(i)

!           Beachte: Fuer die Windkompnenten wurden die unteren
!           Randwerte der Prandtl-Schicht nicht berechnet, daher
!           muss der Faktor tfm hier bei grad(u_m) und grad(v_m)
!           auftauchen!

            grad(i,tet_l)=(ta_2d(i)-ts_2d(i))/dz_0a_h(i) &
                           +tet_g*(z1+clcw(i,ke)*lhocp*virt/ts_2d(i))
            grad(i,h2o_g)=(qda_2d(i)-qds_2d(i))/dz_0a_h(i)

            fm2_2d(i)=grad(i,u_m)**2+grad(i,v_m)**2
            fh2_2d(i)=grav*(a_tet*grad(i,tet_l)+a_vap*grad(i,h2o_g))

            !Beachte:
            !'grad(i,tet_l)' ist bislang der Temperaturgradient und haette zuvor durch
            !den Exnerfaktor dividiert werden muessen, so wie 'a_tet' zuvor damit
            !haette multipliziert werden muessen. Der Exnerfaktor kuerzt sich aber
            !in 'fh2' wieder heraus. Im folgenden wird aber der Teta-Gradient benoetigt:

            grad(i,tet_l)=grad(i,tet_l)/exnr
         END DO      

         END IF

!        Berechnung der atmosphaerischen Forcing-Funktion:

         IF (it_durch.EQ.it_start .AND. lini) THEN !Startinitialisierung

            DO i=i_st,i_en

               IF (fh2_2d(i).GE.(z1-rim)*fm2_2d(i)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=l_tur_z0(i)*(b_2-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2_2d(i)/(fm2_2d(i)-fh2_2d(i))
                  lm=l_tur_z0(i)*(b_2-(a_6+a_3)*fakt)
                  lh=l_tur_z0(i)*(b_1-a_5*fakt)
               END IF

               val1=lm*fm2_2d(i); val2=lh*fh2_2d(i)  
               frc_2d(i)=MAX( val1-val2, rim*val1 )

               tke(i,ke1,nvor)=MAX( SQRT(d_m*l_tur_z0(i)*frc_2d(i)), vel_min )
            END DO       

         ELSE ! mit Hilfe der vorhergehenden TKE-Werte 

            DO i=i_st,i_en
               val1=tkvm(i,ke1)*fm2_2d(i); val2=tkvh(i,ke1)*fh2_2d(i)
               frc_2d(i)=MAX( val1-val2, rim*val1 )/tke(i,ke1,nvor)
            END DO       
         END IF  

! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes:

         fakt=d_4/(z1-c_g)
         w1=tkesmot; w2=z1-tkesmot

         DO i=i_st,i_en

            wert=fakt*l_tur_z0(i)

            q0=tke(i,ke1,nvor)

            frc_2d(i)=MIN( q0**2/wert, frc_2d(i) )

            IF (imode_tran.EQ.1) THEN
               q3=SQRT(d_m*l_tur_z0(i)*frc_2d(i))
            ELSE
               q1=d_m*l_tur_z0(i)/dt_tke
               q2=q0+frc_2d(i)*dt_tke
               q3=q1*(sqrt(z1+z4*q2/q1)-z1)*z1d2
            END IF

            q1=SQRT(wert*frc_2d(i))
            tke(i,ke1,ntur)=MAX( vel_min, q1, w1*q0 + w2*q3 )

!           Was die diversen Beschraenkungen durch MAX- und MIN-
!           Funktionen betrifft, siehe bei den entsprechenden Stellen
!           in 'turbdiff.incf'.

!           frc_2d wird noch zusaetzlich nach oben beschraenkt, damit
!           dieser Wert im Zusammenhang mit q0 keine Singularitaet
!           in der Stabilitaetsfunktion verursachen wuerde. Dies traegt zur
!           numerischen Stabilitaet der iterativ ueber den Zeitindex erfolgenden
!           Bestimmung der Vertikalgradienten bei.

         END DO

         IF (PRESENT(edr)) THEN
            DO i=i_st,i_en
               edr(i,ke1)=tke(i,ke1,ntur)**3/(d_m*l_tur_z0(i))
            END DO
         END IF   


! 4g)    Bestimmung der neuen stabilitaetsabh. Laengenskalen:

         CALL stab_funct(sm=tkvm(:,ke1), sh=tkvh(:,ke1), fm2=fm2_2d, fh2=fh2_2d, &
                         frc=frc_2d, tvs=tke(:,ke1,ntur), tls=l_tur_z0, &
                         i_st=i_st,i_en=i_en)

!-----------------------------------------------------------------------
#ifdef SCLM
         IF (lsclm) THEN
            CALL turb_stat(lsm=tkvm(im,jm,ke1)*l_tur_z0(im,jm), lsh=tkvh(im,jm,ke1)*l_tur_z0(im,jm),   &
                           fm2=fm2_2d(im,jm),                   fh2=fh2_2d(im,jm),                     &
                           d_m=d_m,                             d_h=d_h,                               &
                           tls=l_tur_z0(im,jm),                 tvs=tke(im,jm,ke1,ntur),               &
                           tvt=z0,                              grd=grad(im,jm,:),                   k=ke1 )
         END IF
#endif
!SCLM-------------------------------------------------------------------

         DO i=i_st,i_en

! 4h)       Bestimmung der durch Wirkung der L-Schicht
!           korrigierten Diffusionskoeffizienten 
!           und der zugehoerigen Transferkoeffizienten:

!           unkorrigierte Diffusionskoeffizienten:

            fakt=l_tur_z0(i)*tke(i,ke1,ntur)
            tkvm(i,ke1)=fakt*tkvm(i,ke1)
            tkvh(i,ke1)=fakt*tkvh(i,ke1)

!           Belegung der Felder fuer die Transferkoeffizienten:

            tcm(i)=tkvm(i,ke1)*tfm(i)/(dz_0a_m(i)*vel_2d(i))
            tch(i)=tkvh(i,ke1)*tfh(i)/(dz_0a_h(i)*vel_2d(i))

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
                  velh=(tke(i,ke1,ntur)+tke(i,ke,ntur))*z1d2
                  wert=tcm(i)*vel_2d(i)*SQRT(vel_2d(i)**2+velh**2)
                  gz0(i)=MAX( grav*len_min, alpha0*wert )
!                 gz0(i)=MAX( grav*len_min, alpha0*wert+grav*alpha1*con_m/SQRT(wert) )
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

      IF (iini.EQ.1) THEN !only for separate initialization before the time loop

         RETURN !finish this subroutine

      END IF

!----------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      DO i=i_st,i_en

         IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile
           z2m_2d (i) = z2-h_can_2d(i) !2m ueber dem Bodenniveau des Bestandes
           z10m_2d(i) = z10-z0m_2d(i) !Hoehe in der turbulente Distanz 10m betraegt
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

      CALL diag_level(i_st, i_en, z2m_2d, k_2d, hk_2d, hk1_2d)

      IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile

         DO i=i_st,i_en
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               prs_2d(i)=ps(i)

               IF (z2m_2d(i).LT.z0) THEN
!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

                  fakt=dz_s0_h(i)/dz_sa_h(i)
                  IF (dz_s0_h(i) >= epsi) THEN
                    fakt=MIN( z1, MAX( z0, fakt*EXP(z2m_2d(i)/dz_s0_h(i)) ) )
                  ELSE
! mr: Der Exponent x geht gegen -inf., und 0<=fakt<epsi, so dass fakt*EXP(x) gegen Null geht
                    fakt=z0
                  ENDIF
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
               t_2m(i) = t_g(i) + (ta_2d(i)-t_g(i))*fakt &
                         + tet_g*( (h_atm_2d(i)+h_can_2d(i) )*fakt-h_2m )       !Achtung!
               qv_2m(i)= qv_s(i) + (qda_2d(i)-qv_s(i))*fakt
            END IF
         END DO   

      ELSE !using only a logarithmic profile above a SYNOP lawn
               
         DO i=i_st,i_en
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               prs_2d(i)=ps(i)

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
               t_2m(i) = t_g(i) + (ta_2d(i)-t_g(i))*fakt &
                         + tet_g*(h_atm_2d(i)*fakt-h_2m)
               qv_2m(i)= qv_s(i) + (qda_2d(i)-qv_s(i))*fakt
            END IF   
         END DO   

      END IF   

      DO i=i_st,i_en
         IF (k_2d(i).GT.ke) THEN

!           2m-Niveau liegt oberhalb der untersten Hauptflaeche:

            k = k_2d(i)
            prs_2d(i)=prs(i,k+1)

            IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile
               val1=(z2m_2d(i)+z0m_2d(i))/(hk1_2d(i)+z0m_2d(i))
               val2=(hk_2d(i)+z0m_2d(i))/(hk1_2d(i)+z0m_2d(i))
            ELSE !using only a logarithmic profile above a SYNOP lawn
               val1=(h_2m+z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
               val2=(hk_2d(i)+z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
            END IF

            fakt=LOG(val1)/LOG(val2)

            t_2m(i) =  t(i,k+1) + fakt*( t(i,k)- t(i,k+1)) &
                       + tet_g*( (hk_2d(i)-hk1_2d(i))*fakt - (h_2m-hk1_2d(i)) )
            qv_2m(i) = qv(i,k+1) + fakt*(qv(i,k)-qv(i,k+1))

         END IF
      END DO

      DO i=i_st,i_en

!        Berechnung der zugehoerigen Feuchtewerte:

         wert=grav*(z2m_2d(i)-hk1_2d(i)) &
             /(r_d*t_2m(i)*(z1+rvd_m_o*qv_2m(i))) 
         patm=prs_2d(i)*EXP(-wert)*qv_2m(i) &
             /(rdv+(z1-rdv)*qv_2m(i))               !Wasserdampfdruck

         fakt=patm/zpsat_w( t_2m(i) )
         rh_2m(i)=100.0_ireals*MIN( fakt, z1 )      !relative Feuchte

         wert=LOG(patm/b1)
         td_2m(i)=MIN( (b2w*b3-b4w*wert) &
                   /(b2w-wert), t_2m(i) )           !Taupunktstemperatur

      ENDDO

!     Diagnose der 10m-Groessen:

      CALL diag_level(i_st, i_en, z10m_2d, k_2d, hk_2d, hk1_2d)

      DO i=i_st,i_en

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

            k = k_2d(i)

!           10m-Niveau oberhalb der untersten Modell-Hauptflaeche:

!US this has been forgotten before
            uk =u(i,k)
            vk =v(i,k)
!US

            uk1=u(i,k+1)
            vk1=v(i,k+1)

            val1=(h_10m+z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))
            val2=(hk_2d(i)+z0d_2d(i))/(hk1_2d(i)+z0d_2d(i))

            fakt=LOG(val1)/LOG(val2)

            u_10m(i)=uk1+fakt*(uk-uk1)
            v_10m(i)=vk1+fakt*(vk-vk1)
         END IF

      END DO


CONTAINS



!********************************************************************************

!+ Module procedure stab_funct in "src_turbdiff" for computing the stability
!+ function for vertical diffusion


SUBROUTINE stab_funct (sm, sh, fm2, fh2, frc, tvs, tls, i_st,i_en)

   IMPLICIT NONE

   INTEGER (KIND=iintegers) :: i !loop indices

   REAL (KIND=ireals)        , INTENT(OUT) :: &
        sm(:),  & !stablility function for momentum             [1] 
        sh(:)     !stablility function for scalars (heat)       [1] 

   REAL (KIND=ireals)        , INTENT(IN) :: &
        fm2(:), & !squared forcing frequency for momentum       [1/s2]
        fh2(:), & !squared forcing frequency for scalars (heat) [1/s2]
        frc(:), & !forcing function for TKE                     [m/s2]
        tvs(:), & !turbulent velocity scale SQRT(2*TKE)         [m/s]
        tls(:)    !turbulent length scale

   INTEGER (KIND=iintegers), INTENT(IN) :: &
        i_st, i_en    !start- and end index of horizontal i-loop

   REAL (KIND=ireals)          :: &
!       d0, &
        d1, d2, d3, d4 ,d5 ,d6

   REAL (KIND=ireals) :: &
        a11, a12, a22, a21, &
        a3, a5, a6, be1, be2, &
        gama, fakt, val1, val2, &
        gm, gh

!------------------------------------------------------------------

   DO i=i_st, i_en

      d1=dd(i,1); d2=dd(i,2); d3=dd(i,3)
      d4=dd(i,4); d5=dd(i,5); d6=dd(i,6)

      gama=tls(i)*frc(i)/tvs(i)**2 !entspr. 1/d_m im Gleichgewicht
                                         !und ausserh. des Bestandes
!test gama=1/d0

      fakt=(tls(i)/tvs(i))**2 

!     Folgnde Fallunterscheidung muss gemacht werden,
!     um positiv definete Loesungen fuer die Stabilitaets-
!     funktionen zu ermoeglichen:

      IF (fh2(i).GE.z0) THEN ! stab. Schichtung
!        Allgemeinste der hier verwendeten Loesungen:

         be1=z1
         be2=be1-c_g
       
         gh=fh2(i)*fakt
         gm=fm2(i)*fakt
         
         a11=d1+(d5-d4)*gh
         a12=d4*gm
         a21=(d6-d4)*gh
         a22=d2+d3*gh+d4*gm

         fakt=a11*a22-a12*a21
         sh(i)=(be1*a22-be2*a12)/fakt
         sm(i)=(be2*a11-be1*a21)/fakt
      ELSE ! labile Schichtung
!        Weiter eingeschraenkte Loesungen, bei denen Gm u.
!        Gh unter Einfuehrung der Ri-Zahl eliminiert werden.
!        Dabei wird gama aus der zuvor geloesten TKE-Gleich.
!        genommen. Wegen frc=tls*(sm*fm2-sh*fh2)
!        ist dies aber von den Vorgaengerwerten von sh u. sm
!        aghaengig:

!        Um physikalisch unsinnige Loesungen bez. Singulari-
!        taeten zu vermeiden, muessen be1 u. be2 > 0 sein:

         be1=z1-d4*gama
         be2=be1-c_g

!        Weil tvs schon vorher entsprechend nach unten beschraenkt wurde,
!        ist eine explizite Beschraenkung hier unnoetig!

         a3=d3*gama/d2
         a5=d5*gama/d1
         a6=d6*gama/d2
         be1=be1/d1
         be2=be2/d2

         val1=(fm2(i)*be2+(a5-a3+be1)*fh2(i))/(z2*be1)
         val2=val1+sqrt(val1**2-(a6+be2)*fh2(i)*fm2(i)/be1)
         fakt=fh2(i)/(val2-fh2(i))
         sh(i)=be1-a5*fakt
         sm(i)=sh(i)*(be2-a6*fakt)/(be1-(a5-a3)*fakt)
      END IF   

   END DO

END SUBROUTINE stab_funct


!********************************************************************************

!+ Module procedure diag_level in "src_turbdiff" for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, zdia_2d, k_2d, hk_2d, hk1_2d)

   IMPLICIT NONE

   INTEGER (KIND=iintegers), INTENT(IN) :: &
!
      i_st, i_en !start end end indices of horizontal domain

   REAL (KIND=ireals), INTENT(IN) :: &
!
      zdia_2d(:)   ! diagnostic height

   INTEGER (KIND=iintegers), INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                    !to be used of near surface diagnostics

   REAL (KIND=ireals), INTENT(INOUT) :: &
!
      hk_2d(:), & !mid level height above ground belonging to 'k_2d'
     hk1_2d(:)    !mid level height above ground of the previous layer (below)

   INTEGER (KIND=iintegers) :: i, j

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

END SUBROUTINE turbtran

!********************************************************************************

!+ Module procedure turbdiff in "src_turbdiff" for computing the tendencies
!+ for vertical diffusion


SUBROUTINE turbdiff(iini,lstfnct, dt_var,dt_tke, nprv,ntur,ntim, &
!
          ie, ke, ke1, kcm, &
          istart, iend, istartpar, iendpar, &
          ntstep, &
!    
          l_hori, hhl, dp0, &
!    
          fr_land, depth_lk, sai, &
!
          d_pat, c_big, c_sml, r_air, &
!    
          h_ice, ps, t_g, qv_s, &
          u, v, w, t, qv, qc, prs, rho, epr, &
!    
          gz0, tcm, tch, tfm, tfh, tfv, &
          tke, tkvm, tkvh, rcld, &
          edr, tket_sso, tket_hshr, tket_conv, &
!    
          u_tens, v_tens, t_tens, qv_tens, qc_tens, tketens, &
          qvt_diff, ut_sso, vt_sso, &
!    
          shfl_s, lhfl_s, &
!
          ierrstat, errormsg, eroutine)

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
!     Zirkulationen auf die TKE-Produktion zu beruecksichtigen, was 
!     durch eine Parametrisierung des Drucktransport-Termes erfolgt.
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
!  Changing the treatment of parameter field dd(:,:,:) to achiev better vectorisation
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

! 0. Parameters controlling the call of 'turbdiff':


LOGICAL, INTENT(IN) :: &
!
   lstfnct         !calculation of stability function required

REAL (KIND=ireals), INTENT(IN) :: & 
!
   dt_var,       & !time step for ordinary prognostic variables
   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
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
!   Note, that the indices for the wind-speeds u and v differ from 
!   the other ones because of the use of the staggered Arakawa-B-grid.
!    
! Zonal direction:
    istart,    iend,    & ! start and end index for the inner model domain
    istartpar, iendpar    ! start and end index including model boundary lines

! Time step indices:
! -----------------------------------------------------------------------

INTEGER (KIND=iintegers), OPTIONAL, INTENT(IN) :: &
!
    ntstep          ! current time step

REAL (KIND=ireals), INTENT(IN) :: &
!
! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------
!
    l_hori          ! horizontal grid spacing (m)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
#endif
!
    hhl             ! height of model half levels                   ( m )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
#endif
!
    dp0             ! pressure thickness of layer                   (pa )

#ifdef  __xlC__
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

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
!
    d_pat           ! horizontal pattern length scale

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,kcm:ke1), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,kcm:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml,        & ! effective drag coefficient of canopy elements
!                   ! smaller than the turbulent length scale            (1/m)
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)
#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(IN) :: &
#endif
!
! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------
!
    h_ice,        & ! ice thickness                                 (  m  )
    ps,           & ! surface pressure                              ( pa  )
    t_g,          & ! weighted surface temperature                  (  k  )
    qv_s            ! specific water vapor content on the surface   (kg/kg)

#ifdef  __xlC__
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

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), TARGET, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, INTENT(IN) :: &
#endif
!
     prs            ! atmospheric pressure                          ( pa  )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), TARGET, OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(IN) :: &
#endif
!
     rho,         & ! total density of air                          (kg/m3)
     epr            ! exner pressure                                 (1)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(IN) :: &
#endif
!
     w              ! vertical wind speed (defined on half levels)  ( m/s )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:), INTENT(INOUT) :: &
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
     tfm,          & ! factor of laminar transfer of momentum          --
     tfh,          & ! factor of laminar transfer of scalars           --
     tfv             ! laminar reduction factor for evaporation        --
 
! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=ireals), DIMENSION(ie,ke1,ntim), INTENT(INOUT) :: &
!
     tke             ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                     ! (defined on half levels)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,2:ke1), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,2:), INTENT(INOUT) :: &
#endif
!
     tkvm,         & ! turbulent diffusion coefficients for momentum (m/s2 )
     tkvh            ! turbulent diffusion coefficients for heat     (m/s2 )
                     ! and moisture

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(INOUT) :: &
#endif
!
     rcld            ! standard deviation of the saturation deficit        
                     ! (as input and output)                             
                     ! fractional cloud cover (in turbdiff)            --

#ifdef  __xlC__
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
     qv_tens,      & ! qd-tendency                                   ( 1/s )
     qc_tens         ! qw-tendency                                   ( 1/s )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), INTENT(INOUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), INTENT(INOUT) :: &
#endif
!
     tketens         ! tendency of SQRT(2*TKE)                       ( m/s2)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
     ut_sso,       & ! u-tendency due to the SSO-Scheme              ( 1/s )
     vt_sso          ! v-tendency due to the SSO-Scheme              ( 1/s )

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke1), OPTIONAL, INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: &
#endif
!
     edr             ! eddy dissipation rate of TKE (EDR)            (m2/s3)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(IN) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
#endif
!
     tket_conv       ! TKE-tendency due to convective buoyancy       (m2/s3)

#ifdef  __xlC__
REAL (KIND=ireals), DIMENSION(ie,ke), OPTIONAL, INTENT(OUT) :: &
#else
REAL (KIND=ireals), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: &
#endif
!
     tket_sso,     & ! TKE-tendency due to SSO wake production       (m2/s3)
     tket_hshr,    & ! TKE-tendency due to (sep.) horiz. shear       (m2/s3)
!
     qvt_diff        ! qd-tendency due to diffusion                  ( 1/s )

#ifdef  __xlC__
 REAL (KIND=ireals), DIMENSION(ie), OPTIONAL, INTENT(INOUT) :: &
#else
 REAL (KIND=ireals), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: &
#endif
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2) (positive upward)
     lhfl_s          ! latent   heat flux at the surface             (W/m2) (positive upward)

INTEGER (KIND=iintegers), INTENT(INOUT) :: ierrstat

CHARACTER (LEN=*), INTENT(INOUT) :: eroutine
CHARACTER (LEN=*), INTENT(INOUT) :: errormsg

!-----------------------------------------------------------------------
!
! Local parameters:

REAL (KIND=ireals), DIMENSION(ie,0:7), TARGET :: &
     dd        !local derived turbulence parameter

! Local scalars:

!     Lokale logical Variablen:

      LOGICAL lcalcvdif, & !(teil-)implizite Berechnung der Vertikaldiffusion
              lcircterm    !Zirkulationsterm muss berechnet werden

!     Lokale Integer-Hilfsvariablen:

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &
!
              i,k,          & !Diskretis.index fuer lamda, phi, sigma
              n,ii,kk,      & !Indices fuer diverse Schleifen
              ku,ko,k1,k2,  & !unterer und oberer Schicht-Index
              nkorr,        & !Startindex der Variablen mit Gradientkorrektur
              kem,          & !ke oder ke1
              it_durch        !Durchgangsindex der Iterationen

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
!     Platzh. fuer pot. Temp., relat. Wolkenanteil, rezipr. virtueller Faktor,
!     Temperaturtendenz des Wassersaettigungsmisch.verh. und Druck:
!
           teta,rcl,virt,qsat_dT,patm

      REAL (KIND=ireals) :: &
!
           fr_var, & !1/dt_var
           exnr      !Exner-Faktor

      REAL (KIND=ireals) :: &
!
!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 
!     (fh2,fm2) und die Diffusionskoeff. (bzw. Transferkoeff.) fuer 
!     Skalare und Impuls in m^2/s (kh,km): 
!
           fh2,fm2,kh,km, &
!
!     Platzh. fuer horiz. Geschw.-Komponenten und Geschw.betrag:
!
           vel1,vel2,velh, &
!
!     Platzh. fuer SQRT(2*TKE)-Werte:
!
           q0,q1,q2,q3, &  
!
!     Platzh. fuer die Hoehe ueber Grund, Hoehendifferenzen, obere und
!     untere Hoehe, turbulente Laengaenskalen, Kohaerenzlaenge,
!     Dicke der laminaren Grenzschicht,sowie eine Laenge allgemein:
!   
           h,hu,l_turb,l_diss,lh,lm,kohae_len,len, &
!           
           edh, & ! Kehrwert von Schichtdicken
!
!     Platzhalter fuer Druckdifferenzen
!
           dpo,dpu,dp2, &
!
!     Zwischenspeicher fuer rcpv, rcpl und lhocp
!
           zrcpv,zrcpl,zlhocp, &
!
!     Zwischenspeicher fuer
!
           thermik, & !(negative) Auftriebsproduktion an TKE
           phasdif, & !Temperaturtendenz durch Phasendiffusion
!
           liqfak , & !Faktor zur Beruecksichtigung von Fluessigwasser
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

      LOGICAL            :: lo_ice(ie) ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=ireals), TARGET :: &
!
           vari(ie,ke1,ndiff) ,& !(co-)variances of the velocity-components
                                    !reduced set of variables in the first part of 'turbdiff'
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
           rhoh(ie,ke1),& !Luftdichte auf Hauptflaechen
           rhon(ie,ke1),& !Luftdichte auf Nebenflaechen
                             !einschl. des skin-layers
!
           tfr(ie,ke1),&  !thermisches Forcing (1/s2)
!     
!     Hilfsfeld fuer die korrigierten effektiven Gradienten der
!     ndiff diffundierenden Modellvariablen:
!
           a(ie,ke1,ndiff), &
!
!          a() enthaelt zu beginn die Werte der 5 thermodynamischen
!          Koeffizienten dQs/dT,ex_fakt,cp_fakt,A_tet,A_q 
!          auf den Positionen 1 bis 5 des letzten Index.
!
!          Im Falle der impliziten Berechnung von Flussdivergenzen
!          enthaelt a() hierfuer benoetigte Hilfsgroessen
!
!     3-d Hilfsfelder
!
           hlp(ie,ke1), &
!
!     Flusskonversionsmatrix, um von den turb. Flussdichten der 2
!     scalaren Erhaltungsgroessen (bzgl. feuchtadiab. vert. Verrueck.)
!     auf diejenigen der nscal scalaren Modellvariablen zu gelangen:
!
!          flukon(nvel+1:ndiff,nvel+1:nred), &
!      
!          wrid aus Effizienzgruenden erstetzt druch Skalare 'flukon33', 'flukon34, etc.!
!
!     Hilfsfelder fuer eine und zwei Variablenschichten:
!
           lay(ie), lays(ie,2), &
!
           dz0(ie,2), &   !effective depth of the surface layer
!
           pr(ie),   &    !Druck             (am atm. Modellunterrand)
           tp(ie),   &    !Temperatur              ,,       
           qd(ie),   &    !Wasserdampfgehalt       ,,
           ql(ie),   &    !Fluessigwassergehalt    ,,
           z0m(ie),  &    !Rauhigkeitslaenge (fuer Impuls)
           vh0(ie),  &    !horizontal velocity scale at lowest model level
           l_pat(ie),&    !Laengenskala der therm. Inhomogenitaeten
                             !der Erdbodenoberflaeche
!
           frc(ie),  &   ! forcing
           src(ie),  &   ! source term
!
!     Vertikale Gradienten und turbulente Flussdichten verschiedener 
!     thermodynamischer Variablen:
!
           grad(ndiff), flux(ie,ndiff), &
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
!          Pointer fuer Felder des Exnerfaktors und Bestandes-Parameter:
!
           exner(:,:), &
!
           dpat(:),    &
           cbig(:,:), csml(:,:), rair(:,:)

      REAL (KIND=ireals), TARGET :: &
!
!          Seicher-Target fuer obige Pointer:
!
           exner_tar(ie,ke), &
!
           dpat_tar(ie),   &
           cbig_tar(ie,kcm:ke1), csml_tar(ie,kcm:ke1), rair_tar(ie,kcm:ke1)
 
!          Note:
!          The following buffers wouldn't be necessary, if the related pointers above
!          were allowed to be allocated at run time:

      LOGICAL ::  &
!          
           can_fields, & !canopy fields are present
           ltend(ndiff)  !calculation of tendencies required


      ! Additional local variables for inlining of stab_funct
      REAL (KIND=ireals) :: a11, a12, a22, a21, a3, a5, a6, be1, be2, gama, gm, gh

!---- End of header ------------------------------------------------------------

!Declaration of statement functions:

REAL (KIND=ireals) :: &
     zexner, zpres               !Exner factor and its argument

!Exner-factor:
 zexner(zpres)=EXP(rdocp*LOG(zpres/p0ref))

!     nur in LM-Umgebung:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbtran'; lerror=.FALSE.


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
 END IF

 nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

!     Fuer die Turb.par. benutzter Variablensatz auf Hauptflaechen:
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte) 

!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     und hat einen der Werte u_m,v_m,tet_l,h2o_g,liq;
!                        bzw. u_s,v_s,tem,vap,liq.
!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)
!     Am Ende des Programms wird vari mit den (Co-Varianzen) der
!     Geschwindigkeitskomponenten ueberschrieben, die fuer die
!     Berechung der Horizontaldiffusion benoetigt werden.
!     vari() enthaelt spaeter auch die ndiff (nichtlokalen) vertikalen
!     Gradienten und auch die durch Wirkung der subskaligen Kondensation
!     veraenderten (effektiven) vertikalen Gradienten.
!     Zum Schluss enthaelt vari() fuer die turbulente Horizontaldiff.
!     benoetigte Komponenten des turbulenten Spannungstensors.


      fr_var=z1/dt_var

      IF (PRESENT(epr)) THEN
         exner => epr
      ELSE
         exner => exner_tar
         DO k=1, ke
            DO i=istartpar,iendpar    
               exner(i,k)=zexner(prs(i,k))
            END DO
         END DO
      END IF       

      can_fields=.TRUE.
      IF (PRESENT(d_pat)) THEN
         dpat => d_pat
      ELSE
         dpat => dpat_tar
         can_fields=.FALSE.
      END IF
      IF (PRESENT(c_big)) THEN
         cbig => c_big
      ELSE
         cbig => cbig_tar 
         can_fields=.FALSE.
      END IF
      IF (PRESENT(c_sml)) THEN
         csml => c_sml
      ELSE
         csml => csml_tar 
         can_fields=.FALSE.
      END IF
      IF (PRESENT(r_air)) THEN
         rair => r_air
      ELSE
         rair => rair_tar 
         can_fields=.FALSE.
      END IF

      IF (.NOT.can_fields) THEN
         CALL init_canopy(ie=ie, ke=ke, ke1=ke1, kcm=kcm, &
              istartpar=istartpar, iendpar=iendpar, &
              fr_land=fr_land, &
              d_pat=dpat, c_big=cbig, c_sml=csml, r_air=rair) 
      END IF

      ltend(u_s)=PRESENT(u_tens)
      IF (ltend(u_s)) THEN !calculation of tendencies required
         utens => u_tens !'utens' points to the tendency
      ELSE                 !update of ordinary prognostic variables required
         utens => u      !'utens' points to the prognostic variables
      END IF
      ltend(v_s)=PRESENT(v_tens)
      IF (ltend(v_s)) THEN
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

      DO n=1,ndiff
         IF (ltend(n)) THEN !calculation of tendencies required
            tinc(n)=z1        !no time increment multiplication for tendencies
            tinv(n)=fr_var    !division by time increment for variable increments
         ELSE               !update of prognostic variables required
            tinc(n)=dt_var    !time increment multiplication for tendencies
            tinv(n)=z1        !no division by time increment for variable increments
         END IF
      END DO

      !Note:
      !If a tendency field of an ordinary prognostic variable is not present,
      !the related time step increment due to turbulent diffusion will be
      !added to the prognostic variable directly.
             
!-----------------------------------------------------------------------

! 1)  Vorbereitungen:

!     Setzen von Steuerparametern:

      IF (imode_turb.LT.2) THEN
         lcalcvdif=lexpcor
      ELSE
         lcalcvdif=.TRUE.
      END IF   

      IF (lexpcor) THEN
         IF (lnonloc) THEN
            nkorr=1
         ELSE
            nkorr=nvel+1
         END IF
      ELSE
         nkorr=ndiff+1
      END IF

      IF (itype_tran.EQ.3) THEN
         kem=ke1
      ELSE
         kem=ke
      END IF   

      IF (lcpfluc) THEN
         zrcpv=rcpv
         zrcpl=rcpl
      ELSE
         zrcpv=z0  
         zrcpl=z0  
      END IF

!     Berechnung abgeleiteter Parameter:

      CALL turb_param  (l_hori)

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

!     Bestimmung der initialen Werte fuer die laminaren Reduktions-
!     koeffizienten :

      IF (lini.AND.itype_tran.NE.2) THEN 
         DO i=istartpar,iendpar    
            tfh(i)=z1
            tfm(i)=z1
         END DO
      END IF

!     Berechnung der effektiven Bedeckung mit nichtkonvektiven
!     Wasserwolken:

      IF (icldm_turb.EQ.-1) THEN       
!        Keine Wolkenberuecksichtigung;
!        Wolkenwasser wird ignoriert (vollkommen trockenes Schema):

         DO k=1,ke
            DO i=istartpar,iendpar    
               rcld(i,k)=z0
               hlp(i,k)=z0
            END DO
         END DO

         zlhocp=z0 !no condensation
         liqfak=z0 !total water without liquid water
      
      ELSE

      zlhocp=lhocp; liqfak=z1

      IF (icldm_turb.EQ.0) THEN       

!        Keine Wolkenberuecksichtigung;
!        alles Wolkenwasser wird verdunstet:

         DO k=1,ke
            DO i=istartpar,iendpar    
               rcld(i,k)=z0
               hlp(i,k)=-qc(i,k)
            END DO
         END DO

      ELSEIF (icldm_turb.EQ.1) THEN
!       Wolkenwasser wird als skalige Wolke interpretiert:

        DO k=1,ke
          DO i=istartpar,iendpar
             IF ( qc(i,k) > z0 ) THEN
                rcld(i,k) = z1
             ELSE
                rcld(i,k) = z0
             END IF
             hlp(i,k)=z0
          ENDDO
        ENDDO

      ELSEIF (icldm_turb.EQ.2) THEN
!        spezielle Diagnose von Wasserwolken:

         CALL turb_cloud (ie,  ke,ke1,  1,ke, &
              istartpar,iendpar, 1,ke,        &
              prs, ps, rcld, t, qv, qc,       &
              dicke, hlp                      )

         DO k=1,ke
            DO i=istartpar,iendpar    
               hlp(i,k) =hlp(i,k)-qc(i,k)
               rcld(i,k)=dicke(i,k)
            END DO
         END DO

      END IF
      END IF
      
      DO i=istartpar,iendpar    
         hlp(i,ke1)=hlp(i,ke)
      END DO

      IF (icldm_tran.EQ.0) THEN
         DO i=istartpar,iendpar    
            rcld(i,ke1)=z0
         END DO
      ELSEIF (icldm_tran.EQ.1) THEN
         DO i=istartpar,iendpar    
            IF ( qc(i,ke) > z0 ) THEN
               rcld(i,ke1) = z1
            ELSE
               rcld(i,ke1) = z0
            END IF
         END DO
      ELSE
         DO i=istartpar,iendpar    
            rcld(i,ke1)=rcld(i,ke)
         END DO
      END IF


!     Berechnung der thermodynamischen Hilfsfelder:

!     Bem.:  Am unteren Modellrand gibt es eine skalige Erdboden-
!            oberflaeche mit einer vertikal nicht aufgeloesten 
!            Rauhigkeit der Laenge d+z0m. Diese von Rauhigkeits-
!            elementen durchsetzte 'Mikro'-Rauhigkeitsschicht zaehlt
!            nicht mehr zum atm. Modellgebiet! Sie stellt einen 
!            Mikrobestand der Hoehe d+z0m dar, in dem die konstante
!            Schliessungslaenge akt*z0m gelten soll.

!            Es wird davon ausgegangen, dass das Niveau der Hoehe
!            hhl(i,ke1) (Erdbodenniveau) gerade der Untergrenze der
!            Prandtlschicht entspricht, also um d+z0m hoeher als die
!            Hoehe des eigentlichen (festen) Erdbodens liegt, (wobei
!            d die Verdraengungshoehe der Rauhigkeitselemente und 
!            z0m die effektive Rauhigkeitslaenge ist.

!            Hierbei wird angenommen, dass es am unteren Modellrand
!            keine turb. Fluessigwasserfluesse gibt, dass also der
!            Randwert fuer liq dem Wert der untersten Modellschicht
!            entspricht. Die Flussdichte fuer tet_l entspricht dann
!            der fuer tet und die fuer h2o_g der fuer vap.

      DO k=1,ke1
         IF (k.LT.ke1) THEN
            DO i=istartpar,iendpar    
               tp(i)=t(i,k)+zlhocp*hlp(i,k)
               qd(i)=qv(i,k)-hlp(i,k)
               ql(i)=qc(i,k)+hlp(i,k)
               pr(i)=prs(i,k)
               a(i,k,2)=exner(i,k) !Exner-Faktor
                                       !(wird spaeter auf Nebenflaechen interpoliert)
            END DO   
         ELSE
!           Beachte:
!           tp(), qd(), ql() und pr() sind vom Durchgang k=ke
!           bereits vorhanden.

!           Es werden diese Variablen fuer den atm. Unterrand bestimmt.
!           Dabei wird ueber die Reduktionskoeff. tfh und tfm 
!           die Wirkung der laminaren Grenzschicht zwischen der 
!           Erdbodenumsatzflaeche (im Abstand einer eventuellen
!           Verdraengungshoehe von der festen Erdbodenoberflaeche)
!           und diesem Niveau beruecksichtigt:
 
            DO i=istartpar,iendpar
               tp(i)=(tp(i)+z1d2*tet_g*(hhl(i,ke)-hhl(i,ke1)))*(z1-tfh(i)) &
                      +(t_g(i)+zlhocp*hlp(i,ke1))*tfh(i)
               qd(i)=qd(i)*(z1-tfh(i)) &
                      +(qv_s(i)-hlp(i,ke1))*tfh(i)
               pr(i)=ps(i) 
               a(i,k,2)=zexner(pr(i)) !Exner-Faktor am Boden
            END DO
         END IF

!        Thermdynamische Hilfsvariablen:

         DO i=istartpar,iendpar    
            exnr   =a(i,k,2)    !Exner-Faktor
            teta   =tp(i)/exnr  !potentielle Temperatur

            rcl    =c_scld/(z1+rcld(i,k)*(c_scld-z1))*rcld(i,k) !Bedeckungsgrad

            virt   =z1/(z1+rvd_m_o*qd(i)-ql(i))      !rezipr. virtueller Faktor
            patm   =(z1-qd(i))*virt*pr(i)            !Partialdruck der trockenen Luft
            qsat_dT=zdqsdt( tp(i), zqvap( zpsat_w( tp(i) ), patm ) ) !dQs/dT
            fakt   =(zlhocp/tp(i)-(z1+rvd_m_o)*virt)/(z1+qsat_dT*zlhocp)

            lay(i)=virt

            a(i,k,1)=z1+zrcpv*qd(i)+zrcpl*ql(i) !Cp/Cpd
            a(i,k,3)=qsat_dT

!           Berechn. der thermodyn. Faktoren A_tet und A_q:
            a(i,k,4)=grav*(z1/teta-rcl*fakt*exnr*qsat_dT)
            a(i,k,5)=grav*(rvd_m_o*virt+rcl*fakt)

!           Berechnung der feucht-potentiellen Temp. und des
!           Gesamtwassergehaltes:
            vari(i,k,tet_l)=teta-ql(i)*zlhocp/exnr
            vari(i,k,h2o_g)=qd(i)+liqfak*ql(i)
            vari(i,k,liq  )=ql(i)
         END DO

         IF (PRESENT(rho) .AND. k.LT.ke1) THEN
            DO i=istartpar,iendpar    
               rhoh(i,k)=rho(i,k) !Luftdichte
            END DO
         ELSE    
            DO i=istartpar,iendpar    
               rhoh(i,k)=lay(i)*pr(i)/(r_d*tp(i)) !Luftdichte
            END DO
         END IF   

!        Die thermodynamischen Hilfsgroessen wurden hier
!        unter Beruecksichtigung der diagnostizierten
!        Kondensationskorrektur gebildet, indem die entspr.
!        korrigierten Werte fuer tem, qv und ql benutzt wurden.
!        Beachte, dass es zu einer gewissen Inkosistenz kommt,
!        wenn die Dichte und der Exnerfaktor bereits vorhanden sind!

      END DO

!     Beachte:
!     tp(), qd(), ql() und pr() gehoeren jetzt zu k=ke1 

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

!     Berechnung der Modellschichtdicken:
        
      DO k=1,ke
         DO i=istartpar,iendpar    
            dicke(i,k)=hhl(i,k)-hhl(i,k+1)
         END DO
      END DO   

!     Interpolation der thermodyn. Hilfsgroessen im Feld a(),
!     der Wolkendichte und der Luftdichte auf Nebenflaechen:

      DO k=ke,2,-1
         DO i=istartpar,iendpar
            hlp(i,k)=dp0(i,k)+dp0(i,k-1)
         END DO
      END DO
      DO ii=1,ndiff
         DO k=ke,2,-1
            DO i=istartpar,iendpar
               a(i,k,ii)=(a(i,k,ii)*dp0(i,k-1)+a(i,k-1,ii)*dp0(i,k))/hlp(i,k)
            END DO
         END DO
      END DO
      DO i=istartpar,iendpar    
         rhon(i,ke1)=rhoh(i,ke1)
      END DO
      DO k=ke,2,-1
         DO i=istartpar,iendpar
            rcld(i,k)=(rcld(i,k)*dp0(i,k-1)+rcld(i,k-1)*dp0(i,k))/hlp(i,k)
!test:
            rhon(i,k)=(rhoh(i,k  )*dp0(i,k-1)    &       !mass weighted interpolation
                        +rhoh(i,k-1)*dp0(i,k))/hlp(i,k)
!           rhon(i,k)=(rhoh(i,k  )*dicke(i,k-1)  &       !distance weighted interpolation
!                       +rhoh(i,k-1)*dicke(i,k))   &
!                      /(hhl(i,k-1)-hhl(i,k+1))
!           rhon(i,k)=z2*rhoh(i,k)-rhoh(i,k+1)           !reverse of main level interpolation
                                                               !(leads to a numerical instability)
!test

!           Beachte die andere Behandlung der Dichte!
         END DO
      END DO

!     Bestimmung der initialen Werte fuer die Rauhigkeitslaenge
!     ueber Wasserpunkten:

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

        DO i=istartpar,iendpar

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

      IF (lini.AND.itype_tran.EQ.3) THEN 

         DO i=istartpar,iendpar    
            IF (fr_land(i) <= z1d2) THEN

               ! Use ice surface roughness or open-water surface roughness
               ! according to lo_ice
               IF ( lo_ice(i) ) THEN
                  gz0(i)=grav*z0_ice
               ELSE
!                 Bei von Schubspannung abhaengiger Wellenhoehe:

!                 Einfachste Schaetzung der Schubspannung als Impusls-
!                 flussdichte durch die Nebenflaeche ke mit Hilfe
!                 einer diagnostischen TKE und Bestimmung der Rauhig-
!                 keitslaenge unter Anwendung der Charnockformel:

                  l_turb=dicke(i,ke)
                  l_turb=akt*MAX( len_min, l_turb/(z1+l_turb/l_scal) )
                  edh=z2/(hhl(i,ke-1)-hhl(i,ke1))

                  grad(1)=(vari(i,ke-1,1)-vari(i,ke,1))*edh
                  grad(2)=(vari(i,ke-1,2)-vari(i,ke,2))*edh
                  grad(3)=(vari(i,ke-1,3)-vari(i,ke,3))*edh
                  grad(4)=(vari(i,ke-1,4)-vari(i,ke,4))*edh

                  fh2=a(i,ke,4)*grad(tet_l)+a(i,ke,5)*grad(h2o_g)
                  fm2=grad(u_m)**2+grad(v_m)**2

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

                  q0=MAX( vel_min, SQRT(d_m*l_turb*wert) )

                  wert=lm*q0*SQRT(fm2) 
                  gz0(i)=MAX( grav*len_min, alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF   
         END DO   
      END IF


!     Bestimmung der Rauhigkeitslaenge 
!     und der effektiven Dicke der Modell-Prandtlschicht:

      DO i=istartpar,iendpar    
         z0m(i)=gz0(i)/grav
         vh0(i)=MAX( vel_min, &
                       SQRT(wind(i,ke,u_m)**2+wind(i,ke,v_m)**2) )
      END DO

      IF (itype_tran.EQ.2) THEN !Transferkoeffizienten mit turbtran
         ! spezifiche effektive Dicke der Prandtlschicht:
         DO i=istartpar,iendpar    
            dz0(i,mom)=tkvm(i,ke1)*tfm(i)/(vh0(i)*tcm(i))
            dz0(i,sca)=tkvh(i,ke1)*tfh(i)/(vh0(i)*tch(i))
         END DO 
      ELSE 
         ! Unspezifische Hilfskonsturktion fuer neutrale Schichtung:
         DO i=istartpar,iendpar    
            dz0(i,mom)=z0m(i)*LOG(z1d2*dicke(i,ke)/z0m(i)+z1)
            dz0(i,sca)=dz0(i,mom)
         END DO 
         ! Kuerzt sich bei Flussberechnungen wieder heraus!
      END IF

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

      DO k=1,ke1
        DO i=istartpar,iendpar
          len_scale(i,k)=akt*MAX( len_min,len_scale(i,k)/(z1+len_scale(i,k)/l_scal) )
        END DO
      END DO       

!     Initialisierung der Felder fuer tke,tkvh,tkvm:

      IF (lini) THEN  !nur beim allerersten Durchgang

         DO k=2,kem
            DO i=istartpar,iendpar

!              der Einfachheit halber nur lokale Berechnung der
!              vertikalen Gradienten:

               len=len_scale(i,k)

               IF (k.EQ.ke1) THEN
                  val1=z1/dz0(i,mom)
                  val2=z1/dz0(i,sca)
               ELSE
                  val1=z2/(hhl(i,k+1)-hhl(i,k-1))
                  val2=val1
               END IF

               grad(u_m)=(vari(i,k,1)-vari(i,k-1,1))*val1
               grad(v_m)=(vari(i,k,2)-vari(i,k-1,2))*val1
               grad(tet_l)=(vari(i,k,3)-vari(i,k-1,3))*val2
               grad(h2o_g)=(vari(i,k,4)-vari(i,k-1,4))*val2

               fh2=a(i,k,4)*grad(tet_l)+a(i,k,5)*grad(h2o_g)
               fm2=grad(u_m)**2+grad(v_m)**2

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

               lay(i)=MAX( vel_min, SQRT(d_m*len*wert) ) !Initialwert fuer SQRT(2TKE)

               tkvm(i,k)=lm
               tkvh(i,k)=lh

!              Am Anfang konnte noch keine Advektion oder 
!              Diffusion von SQRT(2*TKE) berechnet werden:

               tketens(i,k)=z0
            END DO        

            DO n=1,ntim
                DO i=istartpar,iendpar
                  tke(i,k,n)=lay(i)
                END DO
            END DO

         END DO    

         DO n=1,ntim
            DO i=istartpar,iendpar
               tke(i,1,n)=tke(i,2,n)
            END DO
         END DO

      ELSE

         DO k=2,kem
            DO i=istartpar,iendpar
               tkvh(i,k)=tkvh(i,k)/tke(i,k,nvor)
               tkvm(i,k)=tkvm(i,k)/tke(i,k,nvor)
            END DO
         END DO    

      END IF

!     tkvh und tkvm enthalten jetzt die stabilitaetsabhaengigen 
!     Laengenmasse, nicht die Diffusionskoeffizienten!


! 2)  Berechnung der benoetigten vertikalen Gradienten und
!     Abspeichern auf vari():

!     Am unteren Modellrand:

!     Impuls-Variablen:
      DO i=istartpar,iendpar
         lay(i)=z1/dz0(i,mom)
      END DO
      DO n=1,nvel
         DO i=istartpar,iendpar
            vari(i,ke1,n)=(vari(i,ke,n)-vari(i,ke1,n))*lay(i)
         END DO
      END DO   

!     Scalar-Variablen:
      DO i=istartpar,iendpar
         lay(i)=z1/dz0(i,sca)
         vari(i,ke1,liq)=z0
      END DO

      DO n=nvel+1,nred
         DO i=istartpar,iendpar
            vari(i,ke1,n)=(vari(i,ke,n)-vari(i,ke1,n))*lay(i)
         END DO
      END DO   
          
!     An den darueberliegenden Nebenflaechen:

      IF (lnonloc) THEN

!        Berechnung nicht-lokaler Gradienten:

         DO n=1,ndiff

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
                     lay(i)=tkvm(i,k)
                  END DO
               ELSE   
                  DO i=istartpar,iendpar
                     lay(i)=tkvh(i,k)
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

!        Belegung von dicke() mit den Schichtdicken
!        bzgl. Nebenflaechen:

         DO k=2,ke   
            DO i=istartpar,iendpar
               dicke(i,k)=(hhl(i,k-1)-hhl(i,k+1))*z1d2
            END DO
         END DO      

      ELSE       

!        Berechnung lokaler Gradienten:
         
         DO k=ke,2,-1
            DO i=istartpar,iendpar
               dicke(i,k)=(hhl(i,k-1)-hhl(i,k+1))*z1d2
               hlp(i,k)=z1/dicke(i,k)
            END DO
         END DO   

         DO n=1,ndiff
            DO k=ke,2,-1
               DO i=istartpar,iendpar
                  vari(i,k,n)=(vari(i,k-1,n)-vari(i,k,n))*hlp(i,k)
               END DO
            END DO   
         END DO      

      END IF

! 3)  Hauptschleife: Bestimmung der TKE und der Stabilitaetsfunkt.:

      DO it_durch=it_start,it_end

         !Die Schleife wird nur bei der Initialisierung (d.h. beim ersten Aufruf) wiederholt,
         !um TKE-Gleichgewicht anzunaehern. Die resultierenden TKE-Werte der Zeitstufe 'ntur'
         !gehoeren in diesem Fall dann zur Zeitstufe 'nvar' ('nold' bei "leap-frog" oder 'nnow'
         !bei 2-Zeitebenen) bzgl. der uebrigen prognostischen Variablen.
         !Fuer die folgenden Aufrufe wird die Schleife nur einmal durchlaufen und liefert TKE-Werte
         !die gegenueber den Vorgaengerwerten um einen Zeitschritt weiter in der Zukunft liegen,
         !also wieder zur Zeitstufe 'nvar' der der uebrigen prognostischen Variablen gehoeren.

         IF (kcm.LE.kem) THEN

!           Es gibt einen Rauhigkeitsbestand:

            DO i=istartpar,iendpar

!              Belegung mit Werten gueltig ausserhalb des Bestandes:

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

         DO k=2,kem

!           Berechnung der atmosphaerischen Vertikal-Antriebe:

            DO i=istartpar,iendpar

               lays(i,1)=vari(i,k,u_m)**2+vari(i,k,v_m)**2
               lays(i,2)=a(i,k,4)*vari(i,k,tet_l) &
                          +a(i,k,5)*vari(i,k,h2o_g)
 
               val1=tkvm(i,k)*lays(i,1); val2=tkvh(i,k)*lays(i,2)

               frc(i)=MAX( dd(i,7)*val1, val1-val2 )

!              Durch die MAX-Funktion wird frc so nach unten 
!              beschraenkt, dass im level-2 Gleichgewicht die
!              krit. Rizahl nie ueberschritten wird. In diesem
!              Gleichgewicht entsrpicht dies dann einem von der
!              Windscherung abhaengigen Minimalwert von TKE.

!              Die Stabilitaetsfunktionen werden aber unabhaengig von
!              dieser Beschraenkung berechnet!

!              Stabilitaetskorrektur der turbulenten Laengenskala
!              bei stabilier Schichtung:

               velh=tke(i,k,nvor)
               wert=a_stab*SQRT(MAX( z0, lays(i,2)) )

               len_scale(i,k)=velh*len_scale(i,k)/(velh+wert*len_scale(i,k))

            END DO

!           Berechnung des Antriebs durch Nachlaufproduktion:

            IF (k.GE.kcm) THEN !Innerhalb des Bestandes:

               IF (k.EQ.ke1) THEN
                  DO i=istartpar,iendpar
                     lay(i)=SQRT(wind(i,ke1,u_m)**2 &
                                +wind(i,ke1,v_m)**2)
                  END DO
               ELSE
                  DO i=istartpar,iendpar
                     dpu=dp0(i,k)
                     dpo=dp0(i,k-1)
                     dp2=dpu+dpo

                     vel1=(wind(i,k,u_m)*dpo+wind(i,k-1,u_m)*dpu)/dp2
                     vel2=(wind(i,k,v_m)*dpo+wind(i,k-1,v_m)*dpu)/dp2
                     lay(i)=SQRT(vel1**2+vel2**2+w(i,k)**2)
                  END DO
               END If

               DO i=istartpar,iendpar

!                 Hilfsgroesse zur Berechn. der reduzierten Konstanten:

                  wert=z3*csml(i,k)*len_scale(i,k) &
                                      *lay(i)/tke(i,k,nvor)

!                 Berechnung der modifizierten Modellparameter:

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

!                 TKE-Forcing durch Addition der Nachlaufproduktion zum vertikalen Forcing:

                  lay(i)=frc(i)+cbig(i,k)*lay(i)**3/tke(i,k,nvor)

               END DO

            ELSE !Ausserhalb des Rauhigkeitsbestandes

!              TKE-Forcing gleich vertikales Forcing:

               DO i=istartpar,iendpar
                  lay(i)=frc(i)
               END DO 

            END IF

            IF (.NOT.(lini.OR.k.eq.ke1)) THEN !nicht bei der Initialisierung oder am Unterrand

               IF (lsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
                  !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:

!                 Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO_Tendenzen:
          
                  DO i=istart,iend
                    dpu=dp0(i,k)
                    dpo=dp0(i,k-1)
                    dp2=dpu+dpo

                    vel1=-(ut_sso(i,k)  *wind(i,k  ,u_m)+vt_sso(i,k)  *wind(i,k  ,v_m))*dpo
                    vel2=-(ut_sso(i,k-1)*wind(i,k-1,u_m)+vt_sso(i,k-1)*wind(i,k-1,v_m))*dpu

                    src(i)=MAX( z0, (vel1+vel2)/dp2 )

!                   Beachte: 
!                   Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
!                   erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
!                   Obwohl vel1 und vel2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!
                  END DO   

                  IF (PRESENT(tket_sso)) THEN
                     DO i=istart,iend
                        tket_sso(i,k)=src(i)
                     END DO
                  END IF

                  IF (ltkesso) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
                     DO i=istart,iend
                        lay(i)=lay(i) + src(i)/tke(i,k,nvor)
                     END DO
                  END IF

               END IF

               IF (lconv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
                  !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

!                 Addition die TKE-Quelle durch konvektive Aktivitaet:

                  DO i=istart,iend
                     lay(i)=lay(i) + MAX( z0, tket_conv(i,k)/tke(i,k,nvor) )
                  END DO

!                 Beachte:  Obwohl tket_conv immer positiv sein muesste, wird zur Sicherheit MAX benutzt!
               END IF
            END IF

!           Beruecksichtigung der 3D_Scherungskorrektur:

            IF (itype_sher.GT.1) THEN !3D-turbulence or separate horizontal shear

!              Addition der 3D-Scherungskorrektur zum bisherigen TKE-Forcing:

               DO i = istart, iend
                  lay(i)=lay(i)+tkvm(i,k)*hlp(i,k)
               END DO

            END IF

!           Berechnung der neuen TKE-Werte:

            DO i=istartpar,iendpar

!              Berechnung einiger Hilfsgroessen:

               q0=tke(i,k,nvor)
               l_turb=len_scale(i,k)
               l_diss=dd(i,0)*l_turb

!              SQRT(2*TKE)-Prognose:

               q1=l_diss/dt_tke
               q2=MAX( z0, q0+tketens(i,k)*dt_tke)+lay(i )*dt_tke

!              In tketens steht die zuvor berechnete Advektions- und 
!              Diffusionstendenz (einschl. Drucktransport).
!              Die MAX-Funkt. verhindert, dass Transportterme die TKE
!              negativ machen koennen.

               q3=q1*(sqrt(z1+z4*q2/q1)-z1)*z1d2

               q1=SQRT(l_turb*frc(i)*dd(i,4)/(z1-c_g))
               tke(i,k,ntur)=MAX( vel_min, q1, q3*(z1-tkesmot)+q0*tkesmot )

!              Die MAX-Funktion und die zeitliche Glaettung
!              sind (vor allem wegen der Wirkung des Zirkulationstermes)
!              nur in Ausnahmefaellen noetig und koennen dann
!              vermeiden, dass die TKE zu klein wird.
!              q1 ist ein Minimalwert, um eine moegliche Singuaritaet bei der
!              Berechnung der Stabilitaetsfunktion fuer labile Schichtung
!              zu vermeiden.

!              Sichern des thermischen TKE-Antriebs:

               tfr(i,k)=lays(i,2)

            END DO    

            IF (PRESENT(edr)) THEN
!              Sichern der TKE-Dissipation "edr=q**3/l_diss" u.a. als thermische Quelle:

               DO i=istartpar,iendpar
                  edr(i,k)=tke(i,k,ntur)**3/(dd(i,0)*len_scale(i,k)) !mit gefiltertem tke-Feld
               END DO   
            END IF   

!           Berechnung der neuen stabilitaetsabhangigen Laengenskalen:

            IF (lstfnct) THEN

!------------------------------------------------------------------

! SR stab_funct inlined in order to avoid excessive call overhead
!    CALL stab_funct(sm=tkvm(:,k), sh=tkvh(:,k), fm2=lays(:,1), fh2=lays(:,2), &
!                    frc=frc, tvs=tke(:,k,ntur), tls=len_scale(:,k), &
!                    i_st=istartpar,i_en=iendpar )


              DO i=istartpar, iendpar

                gama=len_scale(i,k)*frc(i)/tke(i,k,ntur)**2 !entspr. 1/d_m im Gleichgewicht
                                                            !und ausserh. des Bestandes

                fakt=(len_scale(i,k)/tke(i,k,ntur))**2 

          !     Folgende Fallunterscheidung muss gemacht werden,
          !     um positiv definite Loesungen fuer die Stabilitaets-
          !     funktionen zu ermoeglichen:

                 IF (lays(i,2).GE.z0) THEN ! stab. Schichtung
          !        Allgemeinste der hier verwendeten Loesungen:

                   be1=z1
                   be2=be1-c_g
       
                   gh=lays(i,2)*fakt
                   gm=lays(i,1)*fakt
         
                   a11=dd(i,1)+(dd(i,5)-dd(i,4))*gh
                   a12=dd(i,4)*gm
                   a21=(dd(i,6)-dd(i,4))*gh
                   a22=dd(i,2)+dd(i,3)*gh+dd(i,4)*gm

                   fakt=a11*a22-a12*a21
                   tkvh(i,k)=(be1*a22-be2*a12)/fakt
                   tkvm(i,k)=(be2*a11-be1*a21)/fakt
                 ELSE ! labile Schichtung
          !        Weiter eingeschraenkte Loesungen, bei denen Gm u.
          !        Gh unter Einfuehrung der Ri-Zahl eliminiert werden.
          !        Dabei wird gama aus der zuvor geloesten TKE-Gleich.
          !        genommen. Wegen frc=tls*(sm*fm2-sh*fh2)
          !        ist dies aber von den Vorgaengerwerten von sh u. sm
          !        aghaengig:

          !        Um physikalisch unsinnige Loesungen bez. Singulari-
          !        taeten zu vermeiden, muessen be1 u. be2 > 0 sein:

                   be1=z1-dd(i,4)*gama
                   be2=be1-c_g

          !        Weil tvs schon vorher entsprechend nach unten beschraenkt wurde,
          !        ist eine explizite Beschraenkung hier unnoetig!

                   a3=dd(i,3)*gama/dd(i,2)
                   a5=dd(i,5)*gama/dd(i,1)
                   a6=dd(i,6)*gama/dd(i,2)
                   be1=be1/dd(i,1)
                   be2=be2/dd(i,2)

                   val1=(lays(i,1)*be2+(a5-a3+be1)*lays(i,2))/(z2*be1)
                   val2=val1+sqrt(val1**2-(a6+be2)*lays(i,2)*lays(i,1)/be1)
                   fakt=lays(i,2)/(val2-lays(i,2))
                   tkvh(i,k)=be1-a5*fakt
                   tkvm(i,k)=tkvh(i,k)*(be2-a6*fakt)/(be1-(a5-a3)*fakt)
                 END IF   

                tkvm(i,k)=len_scale(i,k)*tkvm(i,k)
                tkvh(i,k)=len_scale(i,k)*tkvh(i,k)
              END DO

            END IF

!------------------------------------------------------------------------------------
#ifdef SCLM
            IF (lsclm .AND. it_durch.EQ.it_end) THEN
               CALL turb_stat(lsm=tkvm(im,jm,k),      lsh=tkvh(im,jm,k),        &
                              fm2=lays(im,jm,1),      fh2=lays(im,jm,2),        &
                              d_m=dd(im,jm,0),        d_h=d_h,                  &
                              tls=len_scale(im,jm,k), tvs=tke(im,jm,k,ntur),    &
                              tvt=tketens(im,jm,k),   grd=vari(im,jm,k,:),   k=k)
            END IF   
#endif
!SCLM--------------------------------------------------------------------------------

            IF (ltmpcor) THEN

!              Sichern der TKE-Dissipation als thermische Quelle:

               DO i=istartpar,iendpar
                  tketens(i,k)=tke(i,k,ntur)**3/(dd(i,0)*len_scale(i,k))
               END DO

            END IF

         END DO !k=2,kem

         IF (it_durch.LT.it_end) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF   


      END DO !Iterationen ueber it_durch

      DO i=istartpar,iendpar
         tke(i,1,ntur)=tke(i,2,ntur)
      END DO

!-----------------------------------------------------------------------
#ifdef SCLM 
      IF (lsclm) THEN
         TKE_SCLM%mod(1)%val=tke(im,jm,1,ntur); TKE_SCLM%mod(1)%vst=i_cal
      END IF
#endif 
!SCLM--------------------------------------------------------------------

      IF (iini.EQ.1) THEN !only for separate initialization before the time loop

         DO k=2, kem
            DO i=istartpar,iendpar
               tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
               tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
            END DO
         END DO
 
         RETURN !finish this subroutine

      END IF   

!  4) Berechnung der effektiven turbulenten Vertikalgradienten
!     und weiterer Turbulenzgroessen:

      DO k=2,kem

         IF (ltmpcor) THEN

!           Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:
!              
            DO i=istartpar,iendpar
               lay(i)=a(i,k,2)*vari(i,k,tet_l)-tet_g &
                       +zlhocp*vari(i,k,liq)               !vert. Temperaturgradient
            END DO   

         !  Dies geschieht schon hier, weil im naechsten Schritt das Feld vari()
         !  durch die effiktiven Gradienten ueberschrieben wird.
         END IF    

         DO i=istartpar,iendpar

!           Berech. der thermodynamischen Flusskonversionsmatrix
!           und der Standardabweichnung der Verteilung des
!           Saettigungsdefizites:

            x1=z1/(z1+zlhocp*a(i,k,3))
            x2=a(i,k,3)*a(i,k,2)
            x3=a(i,k,1) !Cp/Cpd

            rcl = (c_scld * rcld(i,k)) / (z1+rcld(i,k)*(c_scld-z1))

            rcld(i,k)=sqrt(len_scale(i,k)*tkvh(i,k)*d_h)* &
                        abs(x2*vari(i,k,tet_l)-vari(i,k,h2o_g))

!           Unter rcld steht jetzt die (geschaetzte) Standardabw.
!           des Saettigungsdefizites!

            flukon33         =z1-rcl*(z1-x1)
            flukon34         =x1*zlhocp/a(i,k,2)*rcl

            flukon53         =-x1*x2*rcl
            flukon54         =x1*rcl
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

            a(i,k,3)=rcl !'a(i,k,3)' enthaelt ab jetzt den Wolkenbedeckungsgrad!

         END DO

         DO i=istartpar,iendpar
!           Berechn. der Diffusionskoeffizienten:

            tkvh(i,k)=MAX( con_h, tkhmin, tkvh(i,k)*tke(i,k,ntur) )
            tkvm(i,k)=MAX( con_m, tkmmin, tkvm(i,k)*tke(i,k,ntur) )

!           tkvh und tkvm enthalten jetzt nicht mehr Diffusions-
!           laengen, sondern Diffusionskoeffizienten in m^2/s!

         END DO
        
         IF (ltmpcor) THEN

!           Berechnung der Temperaturtendenzen durch TKE-Quellen
!           auf Nebenflaechen (ausser der Divergenz des 
!           Drucktransportes):

            DO i=istartpar,iendpar
               thermik=tkvh(i,k)*tfr(i,k)
               phasdif=tkvh(i,k)*lay(i) &
                      *(zrcpv*vari(i,k,vap)+zrcpl*vari(i,k,liq))  

               tketens(i,k)=len_scale(i,k)/a(i,k,1) &
                       *((tketens(i,k)+thermik)/cp_d+phasdif)
            END DO   

!           Beachte:
!           In tketens() steht bisher die TKE-Dissipation.
!           Wegen der spaeteren Interpolation auf Hauptflaechen,
!           wird mit der turbulenten Laengenskala multipliziert.
         END IF
   
      END DO !k=1, kem

!  5) Untere Randwerte der Turbulenzgroessen:

      IF (kem.EQ.ke1) THEN
!        Untere Randwerte der Turbulenzgroessen und der Transferkoeffizienten
!        mit Hilfe eines vereinfachten Verfahrens zu Testzwecken:

         DO i=istartpar,iendpar

            km=tkvm(i,ke1); kh=tkvh(i,ke1)

!           Dicke der spec. laminaren Grenzschicht:
            val1=z0m(i)*SQRT(con_m/km)
            val2=z0m(i)*SQRT(con_h/kh)

!           Reduktionsfaktoren:
            val1=val1/dz0(i,mom)*(km/con_m)             
            val2=val2/dz0(i,sca)*(kh/con_h)             
            tfm(i)=z1/(z1+rlam_mom *val1)
            tfh(i)=z1/(z1+rlam_heat*val2)
            tfv(i)=(z1+rlam_heat*fakt)/(z1+rat_lam*rlam_heat*val2)
                     !Reduktionsfaktor fuer die Verdunstung aufgrund eines um den
                     !Faktor 'rat_lam' gegenueber fuehlbarer Waerme vergroesserten
                     !laminaren Transpostwiderstandes

            val1=z1/(vh0(i)*dz0(i,mom))
            val2=z1/(vh0(i)*dz0(i,sca))
            tcm(i)=val1*km*tfm(i)
            tch(i)=val2*kh*tfh(i)

!           Berechnung der neuen Rauhigkeitslaenge ueber Meer:

            IF (fr_land(i) <= z1d2) THEN

               ! Use ice surface roughness or open-water surface roughness
               ! according to lo_ice
               IF ( lo_ice(i) ) THEN
                  gz0(i)=grav*z0_ice
               ELSE
!                 Berechnung der Schubspannung:

                  wert=tcm(i)*vh0(i)*SQRT(vh0(i)**2+tke(i,ke1,ntur)**2)
                  
!                 grav*z0 mittels Charnock-Formel: 

                  gz0(i)=MAX( grav*len_min, alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF

         END DO   

      ELSE
!        Untere Randwerte der Turbulenzgroessen mit Hilfe der
!        bereits vorher (in turbtran oder partur(b,s)) berechneten
!        Transferkoeffizienten:

!        Beachte, dass die Konversion der eff. Vertikalgradienten am Unterrand
!        noch nicht stattgefunden hat.

         IF (itype_tran.NE.2) THEN

!           Es wurde i.a. noch keine unteren Randwerte fuer 
!           tke(), sowie tkvm(), tkvh() berechnet:

            DO i=istartpar,iendpar
               tkvm(i,ke1)=vh0(i)*dz0(i,mom)*tcm(i)
               tkvh(i,ke1)=vh0(i)*dz0(i,sca)*tch(i)

!              Weil in partur(b,s) keine laminare Grenzschicht
!              beruecksichtigt wird, sind die Reduktionskoeff.
!              immer gleich 1:

               tfh(i)=z1
               tfm(i)=z1

!              Der SQRT(2*TKE)-Wert am unteren Modellrand wird mit 
!              Hilfe von c_tke*Ustar bestimmt: 

               tke(i,ke1,ntur)=c_tke*vh0(i)*SQRT(tcm(i))
            END DO   

         END IF

         IF (ltmpcor) THEN  
            DO i=istartpar,iendpar
               lay(i)=a(i,ke1,2)*vari(i,ke1,tet_l)-tet_g !vert. Temperaturgradient
               !Beachte, dass hier vari(i,ke1,liq)=0 ist!
            END DO 
         END IF
   

         DO i=istartpar,iendpar

!           Thermischer Antrieb (fh2):

            tfr(i,ke1)=a(i,ke1,4)*vari(i,ke1,tet_l) &
                        +a(i,ke1,5)*vari(i,ke1,h2o_g)

!           Thermodynamische Korrektur des effektiven vertikalen
!           Gradienten der pot. Temperatur am unteren Modellrand:

            vari(i,ke1,tet)=a(i,ke1,1)*vari(i,ke1,tet) !(Cp/Cpd)*grad(tet)

            !Beachte:
            !Am Unterrand wurde grad(liq) = 0 gesetzt. Daher ist dort 
            !grad(tet_l) = grad(tet) und grad(h2o_g)=grad(vap).

            a(i,ke1,3) = (c_scld * rcld(i,ke1)) / (z1+rcld(i,ke1)*(c_scld-z1))
            
            !'rcld(i,ke1)' ist in diesem Fall immer noch der Bedeckungsgrad,
            !dessen durch 'c_scld' abgewandelter Wert nun auch in 'a(i,ke1,3)' steht.
         END DO   

         IF (ltmpcor) THEN  

!           Bestimmung der unteren Randwerte der Enthalpieproduktion
!           durch Wechselwirkung mit der TKE (diss und thermik),
!           sowie durch Phasendiffusion:
          
            DO i=istartpar,iendpar

               thermik=tkvh(i,ke1)*tfr(i,ke1)
               phasdif=tkvh(i,ke1)*lay(i) &
                      *(zrcpv*vari(i,ke1,vap)+zrcpl*vari(i,ke1,liq))
               !Beachte, dass eigentlich vari(i,ke1,liq=0 ist!

               wert=tke(i,ke1,ntur)**3/(d_m*len_scale(i,ke1))
               tketens(i,ke1)=len_scale(i,ke1)/a(i,ke1,1) &
                      *((wert+thermik)/cp_d+phasdif)

            END DO   
                
         END IF    
         
      END IF !untere Randwerte

! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperatur-
!     tendenzen ausser der Divergenz des Drucktransportes:

      IF (ltmpcor) THEN
         DO i=istart,iend
            ttens(i,1)=ttens(i,1)+tinc(tem)*tketens(i,2) &
                         /(len_scale(i,1)+len_scale(i,2))

         END DO
         DO k=2,ke
            DO i=istart,iend
               ttens(i,k)=ttens(i,k)+tinc(tem)*(tketens(i,k)+tketens(i,k+1)) &
                            /(len_scale(i,k)+len_scale(i,k+1))
            END DO
         END DO
      END IF   


! 7)  Bestimmung des Drucktransporttermes (Zirculationstermes):

      lcircterm=(pat_len.GT.z0)

      IF (lcircterm) THEN !Der Zirculationsterm (Drucktransportterm muss berechnet werden)

         DO i=istartpar,iendpar
           IF (dpat(i).GT.z0) THEN
             l_pat(i)=pat_len*dpat(i) 
             l_pat(i)=l_hori*l_pat(i)/(l_hori+l_pat(i))
           ELSE
              l_pat(i)=z0
           END IF
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
                  lay(i)=pr(i)
               END DO
            END IF

            DO i=istartpar,iendpar

               fakt=z1-z2*ABS(a(i,k,3)-z1d2)
               len=MAX( l_pat(i), SQRT(fakt*len_scale(i,k)*l_hori) )

!              Berechnung der lokalen Kohaerenzlaenge fuer die 
!              Drucktransport-Parameterisierung:
!              kohae_len=l_pat*exnr*grad(tet_v)*r/g:

               fakt=tfr(i,k)*lay(i)/(rhon(i,k)*grav**2)
               kohae_len=len*SIGN(z1,fakt)*MIN( ABS(fakt), z1 )

!              Belegung von tketens mit dem zur vert. turb. 
!              TKE-Flussdichte mittels Drucktransport gehoerigen
!              TKE-Gradienten:

               tketens(i,k)=kohae_len*tfr(i,k)

!              Die Divergenz des zugehoerigen Flusses ist gleichzeitig
!              eine weitere Quelle fuer thermische Energie.

            END DO

!           Addition des Gradienten, welcher zur Temperatur-
!           flussdichte durch Drucktransport gehoert:

            IF (ltmpcor) THEN
               DO i=istartpar,iendpar
                  vari(i,k,tet)=vari(i,k,tet)-tketens(i,k)/cp_d
               END DO
            END IF

         END DO   

      ELSE !.NOT.lcircterm

         DO k=2,ke1
            DO i=istartpar,iendpar
               tketens(i,k)=z0
            END DO
         END DO   
         
      END IF   

!     Effektive Schichtdicke multipliziert mit der Dichte:

      IF (lcalcvdif) THEN
         kk=1
      ELSE
         kk=kcm
      END IF

      DO k=kk,ke
         DO i=istartpar,iendpar
            dicke(i,k)=rhoh(i,k)*(hhl(i,k)-hhl(i,k+1))
         END DO
      END DO

! 8)  Explizite Berechnung von Bestandtendenzen:

      IF (kcm.LE.ke) THEN

!        Berechnung des Formwiderstandes im Bestand:

         ku=1
         ko=2
         DO i=istartpar,iendpar
            lays(i,ku)=z0
         END DO
         DO k=ke,kcm,-1   
            DO i=istartpar,iendpar
               lays(i,ko)=w(i,k)
               wert=z1d2*(cbig(i,k)  +csml(i,k) &
                         +cbig(i,k+1)+csml(i,k+1)) &
                   *SQRT(wind(i,k,u_m)**2+wind(i,k,v_m)**2 &
                        +(z1d2*(lays(i,ku)+lays(i,ko)))**2)

!              Durch die min-Funktion, wird verhindert, dass durch
!              die Wirkung der Formreibung eine Richtungsumkehr
!              der Geschwindigkeit erfolgen kann:

               wert=MIN( fr_var,wert )

               wind(i,k,u_m)=-wert*wind(i,k,u_m)*tinc(u_m)
               wind(i,k,v_m)=-wert*wind(i,k,v_m)*tinc(v_m)
            END DO   
            kk=ku
            ku=ko
            ko=kk
         END DO   

!        Berechnung des Volumeneffektes im Bestand:

         DO n=1,ndiff
            ku=1
            ko=2
            IF (n.LE.nvel) THEN
               DO i=istartpar,iendpar
                  lays(i,ku)=-rhon(i,ke1)*tkvm(i,ke1) &
                                             *vari(i,ke1,n)
               END DO
            ELSE
               DO i=istartpar,iendpar
                  lays(i,ku)=-rhon(i,ke1)*tkvh(i,ke1) &
                                             *vari(i,ke1,n)
               END DO
            END IF

            DO k=ke,kcm,-1
               IF (kcm.EQ.1) THEN
                  DO i=istartpar,iendpar
                     lays(i,ku)=lays(i,ku)*z1d2 &
                         *(rair(i,k)-rair(i,k+1))/dicke(i,k)
                  END DO
               ELSE   
                  DO i=istartpar,iendpar
                     lays(i,ko)=-rhon(i,k)*tkvh(i,k)*vari(i,k,n)
                     lays(i,ku)=(lays(i,ko)+lays(i,ku))*z1d2 &
                         *(rair(i,k)-rair(i,k+1))/dicke(i,k)
                  END DO
               END IF
               IF (n.LE.nvel) THEN
                  DO i=istartpar,iendpar
                     wind(i,k,n)=wind(i,k,n)-lays(i,ku)*tinc(n)
                  END DO
               ELSEIF (n.EQ.tem) THEN
                  DO i=istart,iend
                     ttens(i,k)=ttens(i,k)-lays(i,ku)*tinc(tem)
                  END DO
               ELSEIF (n.EQ.vap) THEN
                  DO i=istartpar,iendpar
                     qvtens(i,k)=qvtens(i,k)-lays(i,ku)*tinc(vap)
                  END DO
               ELSEIF (n.EQ.liq) THEN
                  DO i=istart,iend
                     qctens(i,k)=qctens(i,k)-lays(i,ku)*tinc(liq)
                  END DO
               END IF
               kk=ku
               ku=ko
               ko=kk 
            END DO    
         END DO    

      END IF !kcm.LE.ke

! 9) Berechnung des vorgesehenen Anteils der Flussdichtedivergenzen:

      IF (lcalcvdif) THEN
        
         DO k=1,ke
            DO i=istartpar,iendpar
               dicke(i,k)=dicke(i,k)*fr_var
            END DO
         END DO

         DO k=2,ke
            DO i=istartpar,iendpar
               a(i,k,5)=z1d2*(hhl(i,k-1)-hhl(i,k+1))
            END DO
         END DO

         DO i=istartpar,iendpar
            a(i,1,1)=z0
         END DO

         DO n=1, ndiff

!           Vorbereitung zur Bestimmung der vertikalen Diffusions-Tendenzen:

            IF (n.eq.1) THEN

!              Fuer horizontale Windkomponenten:

               DO i=istartpar,iendpar
                  a(i,ke1,5)=dz0(i,mom)/tfm(i)
               END DO
               DO k=2,ke1 
                  DO i=istartpar,iendpar
                     a(i,k,1)=rhon(i,k)*tkvm(i,k)/a(i,k,5)
                     a(i,k,2)=rhon(i,k)*a(i,k,5)*fr_var
                  END DO
               END DO

               CALL prep_impl_vert_diff( &
                 i_st=istartpar, i_en=iendpar, k_tp=0, k_sf=ke1,  &
                 disc_mom=dicke, diff_mom=a(:,:,1), impl_mom=a(:,:,2), invs_mom=a(:,:,3), &
                 lsflucond=(imode_turb.EQ.3) ) 

            ELSEIF (n.eq.nvel+1) THEN  

!              Fuer skalare Variablen:

               DO i=istartpar,iendpar
                  a(i,ke1,5)=dz0(i,sca)/tfh(i)
               END DO
               DO k=2,ke1 
                  DO i=istartpar,iendpar
                     a(i,k,1)=rhon(i,k)*tkvh(i,k)/a(i,k,5)
                     a(i,k,2)=rhon(i,k)*a(i,k,5)*fr_var
                  END DO
               END DO

               CALL prep_impl_vert_diff( &
                 i_st=istartpar, i_en=iendpar, k_tp=0, k_sf=ke1,  &
                 disc_mom=dicke, diff_mom=a(:,:,1), impl_mom=a(:,:,2), invs_mom=a(:,:,3), &
                 lsflucond=(imode_turb.EQ.3) ) 

            END IF

!           Belegung der Eingangsprofile:

            IF (n.EQ.u_m) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     hlp(i,k)=u(i,k)
                  END DO
               END DO   
               DO i=istartpar,iendpar
                  hlp(i,ke1)=z0 !no slip condition
               END DO
            ELSEIF (n.EQ.v_m) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     hlp(i,k)=v(i,k)
                  END DO
               END DO   
               DO i=istartpar,iendpar
                  hlp(i,ke1)=z0 !no slip condition
               END DO
            ELSEIF (n.EQ.tem) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     hlp(i,k)=t(i,k)/exner(i,k)
                  END DO
               END DO   
               DO i=istartpar,iendpar
                  hlp(i,ke1)=t_g(i)/zexner(ps(i))
               END DO
            ELSEIF (n.EQ.vap) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     hlp(i,k)=qv(i,k)
                  END DO
               END DO   
               DO i=istartpar,iendpar
                  hlp(i,ke1)=qv_s(i)
               END DO
            ELSEIF (n.EQ.liq) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     hlp(i,k)=qc(i,k)
                  END DO
               END DO   
               DO i=istartpar,iendpar
                  hlp(i,ke1)=hlp(i,ke) !no cloud water deposition
               END DO
            END IF   

!           Korrigierte Werte durch Vertikalintegration der Gradienten:

            IF (n.GE.nkorr) THEN
               IF (imode_turb.LT.2) THEN !nur Korrekturprofile
                  DO k=ke1,2,-1
                     DO i=istartpar,iendpar
                        hlp(i,k)=hlp(i,k-1)-hlp(i,k)
                     END DO
                  END DO    
                  DO i=istartpar,iendpar
                     hlp(i,1)=z0
                  END DO
                  DO k=2,ke1
                     DO i=istartpar,iendpar
                         hlp(i,k)=hlp(i,k-1)+(hlp(i,k)-a(i,k,5)*vari(i,k,n))
                     END DO
                  END DO
               ELSE !effektives Gesamtprofil
                  DO k=2,ke1
                     DO i=istartpar,iendpar
                        hlp(i,k)=hlp(i,k-1)-a(i,k,5)*vari(i,k,n)
                     END DO
                  END DO
               END IF    
            END IF    

!           Berechnung der vertikalen Diffusionstendenzen der Windkomponenten auf Massenpunkten: 

            CALL calc_impl_vert_diff ( &
                 i_st=istartpar, i_en=iendpar, k_tp=0, k_sf=ke1,  &
                 disc_mom=dicke, diff_mom=a(:,:,1), impl_mom=a(:,:,2), invs_mom=a(:,:,3), &
                 cur_prof=hlp  , upd_prof=a(:,:,4), srf_flux=flux(:,n) )

!           Belegung der Tendenzfelder:

            IF (n.LE.nvel) THEN !Windkomponenten
               DO k=1,kcm-1
                  DO i=istartpar,iendpar
                     wind(i,k,n)=(a(i,k,4)-hlp(i,k))*tinv(n)
                  END DO
               END DO
               DO k=kcm,ke
                  DO i=istartpar,iendpar
                     wind(i,k,n)=wind(i,k,n)+(a(i,k,4)-hlp(i,k))*tinv(n)
                  END DO
               END DO
            ELSEIF (n.EQ.tem) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     ttens(i,k)=ttens(i,k) &
                                  +exner(i,k)*(a(i,k,4)-hlp(i,k))*tinv(n)
                  END DO
               END DO
            ELSEIF (n.EQ.vap) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     qvtens(i,k)=qvtens(i,k) &
                                  +(a(i,k,4)-hlp(i,k))*tinv(n)
                  END DO
               END DO
               IF (PRESENT(qvt_diff)) THEN
                  DO k=1,ke
                     DO i=istartpar,iendpar
                        qvt_diff(i,k)=(a(i,k,4)-hlp(i,k))*fr_var
                     END DO
                  END DO
               END IF
            ELSEIF (n.EQ.liq) THEN
               DO k=1,ke
                  DO i=istartpar,iendpar
                     qctens(i,k)=qctens(i,k) &
                                  +(a(i,k,4)-hlp(i,k))*tinv(n)
                  END DO
               END DO
            END IF

         END DO !n=1,ndiff

      ELSE !.NOT.lcalcvdif
         
         DO n=tet, vap
            DO i=istartpar,iendpar
               flux(i,n)=-rhon(i,ke1)*tkvh(i,ke1)*vari(i,ke1,n)
            END DO
         END DO

      END IF 

! 10) Berechnung der Enthalpieflussdichten:
       
      IF (PRESENT(shfl_s)) THEN
         DO i=istartpar,iendpar
            shfl_s(i)=-zexner(ps(i))*cp_d*flux(i,tet)
         END DO
      END IF
      IF (PRESENT(lhfl_s)) THEN
         DO i=istartpar,iendpar
            lhfl_s(i)=-lh_v*flux(i,vap)
         END DO
      END IF
 
!---------------------------------------------------------------------------------------
#ifdef SCLM 
      IF (lsclm) THEN
         SHF%mod(0)%val=-zexner(ps(im,jm))*cp_d*flux(im,jm,tet); SHF%mod(0)%vst=i_cal
         LHF%mod(0)%val=                  -lh_v*flux(im,jm,vap); LHF%mod(0)%vst=i_cal
      END IF
#endif 
!SCLM-----------------------------------------------------------------------------------

       !Bem: shfl_s und lhfl_s, sowie SHF und LHF sind positiv abwaerts!

! 11) Berechnung der Tendenzen infloge horizontaler turb. Diffusion
!     (und aufsummieren auf die Tendenzfelder):

!----------------------------------------------------------------------
!     IF (lhordiff) THEN

!        Nun werden die noch unter vari(:,:,:,1:2) abgespeicherten effektiven
!        Gradienten nur noch fuer die Berechnung der Komponenten
!        des turbulenten Spannungstensors benoetigt, welche u.U.
!        auch als Eingangsfeld fuer die Bestimmung der turbulenten
!        Horizontaldiffusion gebraucht werden:
         
!        DO k=2,ke1
!           DO i=istartpar,iendpar
!              cvar(u_m,w_m)=-tkvm(i,k)*vari(i,k,u_m)
!              cvar(v_m,w_m)=-tkvm(i,k)*vari(i,k,v_m)

!              q2=z1d3*tke(i,k,ntur)**2
!              fakt=d_m*len_scale(i,k)/tke(i,k,ntur)
!              x1=cvar(u_m,w_m)*vari(i,k,u_m)
!              x2=cvar(v_m,w_m)*vari(i,k,v_m)
!              x3=-tkvh(i,k)*tfr(i,k)

!              vari(i,k,1)=cvar(u_m,w_m)
!              vari(i,k,2)=cvar(v_m,w_m)
!              vari(i,k,3)=q2+fakt*(-z4*x1+z2*x2-z2*x3)
!              vari(i,k,4)=q2+fakt*(+z2*x1-z4*x2-z2*x3)
!           END DO
!        END DO

!        Jetzt werden die in hlp() gespeicherten fh2-Werte nicht mehr
!        benoetigt und hlp() kann anderweitig benutzt werden.

!        Die eigentliche Berechnung der Tendenzen durch horizontaldiffusion
!        fehlt hier noch!
         
!     END IF
!-------------------------------------------------------------------------------------

! 12) Aktualisierung der Tendenzfelder fuer den horizontalen Wind:

      IF (lcalcvdif .OR. kcm.LE.ke) THEN


         IF (lcalcvdif) THEN
            kk=1
         ELSE
            kk=kcm
         END IF

         DO k=kk,ke   
            DO i=istart,iend 
               utens(i,k)=utens(i,k) &
                           +(wind(i,k,u_m)+wind(i,k,u_m))*z1d2
            END DO
            DO i=istart,iend  
               vtens(i,k)=vtens(i,k) &
                           +(wind(i,k,v_m)+wind(i,k,v_m))*z1d2
            END DO
         END DO

      END IF

                  
! 13) Berechnung der Diffusionstendenz von SQRT(2*TKE):

!        Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!        des Saettigungsdefizites und den Drucktransport:

         IF (lcircterm) THEN
            DO k=2,ke1
               DO i=istartpar,iendpar  
                  tketens(i,k)=rhon(i,k)*tkvh(i,k)/tke(i,k,ntur) &
                                            *tketens(i,k)*len_scale(i,k)
               END DO
            END DO
            DO k=2,ke
               DO i=istartpar,iendpar  
                  tketens(i,k)=(tketens(i,k)+tketens(i,k+1)) &
                                /(len_scale(i,k)+len_scale(i,k+1))
               END DO
            END DO
         END IF

         DO i=istartpar,iendpar  
            rcld(i,1)=rcld(i,2)
         END DO
         DO k=2,kem-1
            DO i=istartpar,iendpar  
               rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
            END DO
         END DO
!        Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
!        der Wert auf der entspr. Nebenflaeche beibehalten.

!        tketens wird hiernach als Teilbetrag der SQRT(2*TKE)-Flussdichte
!        behandelt. 

!        Berechnung der vertikalen SQRT(2*TKE)-Flussdichten
!        und Ablegen der selben im Feld tketens:
          
!        Weil die TKE-Gleichung eine Gleichung 2-ter Ordnung ist und
!        die Parametrisierung der Diffusion in dieser Gleichung ohne
!        Zuhilfenahme von Gleichungen hoeherer Ordnung quasi ad hoc 
!        erfolgt, sind die Genauhigkeitsansprueche der numerischen
!        Berechnung ensprechend geringer als bei der Berechnung der
!        Diffusion in den Gleichungen 1-ter Ordnung. Es wird daher
!        hier auf nicht-lokale Gradientbildung und implizite Behandlung
!        generell verzichtet.

!   fuo: For very shallow model levels and TKE values on the order of 1 m2/s2, 
!        the diffusion of TKE is strongly limited by stability constraints,
!        and an implicit solution (by choosing limpltkediff=.true.) should be favored.
!   Ra:  Nevertheless, we need for the present at least an exlicit correction
!        in order to express the circulation term and roughness layer terms.

!        Es soll keinen turb. TKE-FLuss aus der Atm. hinaus geben:

         DO i=istartpar,iendpar
            tketens(i,1)=z0
         END DO

         DO k=2,ke
            DO i=istartpar,iendpar  
               dicke(i,k)=z1d2*(hhl(i,k-1)-hhl(i,k+1))*rhon(i,k)
            END DO
         END DO

         IF (.NOT.limpltkediff) THEN

!           Bei der reinen explizite Diffusion von TKE muss die explizite
!           Flussdichte berechnet werden:

            DO k=ke,2,-1
               DO i=istartpar,iendpar
                  len=hhl(i,k)-hhl(i,k+1)
                  q3=(tke(i,k,ntur)+tke(i,k+1,ntur))*z1d2
                  kh=q3*(len_scale(i,k)+len_scale(i,k+1))*z1d2
                  kh=rhoh(i,k)*MIN( securi*len**2/dt_tke, c_diff*kh )
                     
!-----------------------------------------------------------------------
!                    nur im 1-d-Modell:
   
!                    if (ldruck) then
!                       write (7,*) 'k,tke,p_flux,tke_flux' &
!                      ,k,tke(i,k,ntur),tketens(i,k)    &
!                      ,kh*q3*(tke(i,k+1,ntur)-tke(i,k,ntur))/len
!                    end if
!-----------------------------------------------------------------------

                  tketens(i,k)=tketens(i,k) &
                                +kh*(tke(i,k+1,ntur)-tke(i,k,ntur))/len
               END DO
            END DO

!           Bem.: Der viskosen Zusatzterm kann auch entfallen!
!           fuo:  Note that in order to achieve stability with this
!                 explicit scheme, the following condition 0 <= securi < 0.5
!                 should be met. 

!           Die Diffusionskoeffizienten wurden mit Hilfe der MIN-Funktion
!           auf das numerisch ertraegliche Mass reduziert.

!           Zuvor wurde tketens bereits mit dem Drucktransport belegt,
!           so dass tketens also die Summe aus der eigentlichen TKE-
!           Flussdichte und dem Drucktrainsport enthaelt!

!           Positiv definite Diffusionskorrektur in der Gleichung fuer SQRT(2*TKE)
!           wobei die Laengenskala c_diff*len_scale beschraenkt werden muss:

            DO k=2,ke
               DO i=istartpar,iendpar
                  len=hhl(i,k-1)-hhl(i,k+1)
                  hlp(i,k)=tketens(i,k-1)-tketens(i,k) &
                        +MIN( securi*len**2/(4*tke(i,k,ntur)*dt_tke), c_diff*len_scale(i,k)) &
                             *( (tke(i,k-1,ntur)-tke(i,k+1,ntur))/len )**2 * dicke(i,k)
               END DO
            END DO

         ELSE

!           Bei der impliziten Diffusion von TKE muss nur die Flussdichte bzueglich des
!           Zirkulationstermes explizit behandelt werden:

            DO k=2,ke
               DO i=istartpar,iendpar
                  hlp(i,k)=tketens(i,k-1)-tketens(i,k) 
               END DO
            END DO

         END IF    

!        Volumenkorrektur in der Rauhigkeitsschicht:

         DO k=MAX( 2, kcm ), ke
            DO i=istartpar,iendpar
               hlp(i,k)=hlp(i,k) &
                  +z1d2*(tketens(i,k)*dp0(i,k-1)+tketens(i,k-1)*dp0(i,k)) &
                       /(dp0(i,k)+dp0(i,k-1))*(rair(i,k-1)-rair(i,k+1))
            END DO
         END DO

         DO i=istartpar,iendpar
            hlp(i,1)=hlp(i,2)
            hlp(i,ke1)=hlp(i,ke)
         END DO

!        Die folgende (quellenfreie) vertikale Glaettung (mittels 'wichfakt')
!        ist bei groesseren Zeitschritten (> ca. 30sec) in der expliziten
!        Variante noetig:

         DO k=2,ke
            DO i=istartpar,iendpar
               tketens(i,k)= &
                     -(hlp(i,k)+wichfakt*(hlp(i,k-1)-z2*hlp(i,k) &
                                        +hlp(i,k+1)))/dicke(i,k)
            END DO
         END DO

         IF (limpltkediff) THEN

!           Vorbereitung zur Bestimmung der vertikalen Diffusions-Tendenzen von SQRT(2*TKE):

            DO k=2, ke
               DO i=istartpar,iendpar
                  dicke(i,k)=tke(i,k,ntur)*dicke(i,k)*fr_var
               END DO
            END DO   

            k=2; ku=MOD(k,2)+1
            DO i=istartpar,iendpar
               lays(i,ku)=c_diff*len_scale(i,2)*(tke(i,2,ntur))**2
            END DO

            DO k=3, ke1
               ko=MOD(k-1,2)+1; ku=MOD(k,2)+1
               DO i=istartpar,iendpar
                  lays(i,ku)=c_diff*len_scale(i,k)*(tke(i,k,ntur))**2
                  len=hhl(i,k-1)-hhl(i,k)
                  a(i,k,1)=rhoh(i,k-1)*z1d2*(lays(i,ko)+lays(i,ku))/len
                  a(i,k,2)=rhoh(i,k-1)*len*fr_var
               END DO
            END DO   

            CALL prep_impl_vert_diff( &
                 i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1,  &
                 disc_mom=dicke, diff_mom=a(:,:,1), impl_mom=a(:,:,2), invs_mom=a(:,:,3), &
                 lsflucond=(imode_turb.EQ.3) )

!           Berechnung der vertikalen Diffusionstendenzen von SQRT(2*TKE):

            CALL calc_impl_vert_diff ( &
                 i_st=istartpar, i_en=iendpar, k_tp=1, k_sf=ke1,  &
                 disc_mom=dicke, diff_mom=a(:,:,1), impl_mom=a(:,:,2), invs_mom=a(:,:,3), &
                 cur_prof=tke(:,:,ntur), upd_prof=a(:,:,4) )

            DO k=2,ke
               DO i=istartpar,iendpar
                  tketens(i,k)=tketens(i,k)+ &
                                 (a(i,k,4)-tke(i,k,ntur))/dt_tke
               END DO
            END DO

         END IF

!        Im skin-layer soll es keine TKE-Tendenzen infolge
!        Diffusion oder Drucktransport geben:

         DO i=istartpar,iendpar
            tketens(i,ke1)=z0
         END DO

!        Die folgende MAX-Funktion soll sicherstellen, dass allein durch die Diffusionstendenz
!        von SQRT(2*TKE) (einschl. des Zirkulationstermes) die TKE nicht negativ werden kann.

         DO k=2,ke
            DO i=istartpar,iendpar
               tketens(i,k)=MAX( -tke(i,k,ntur)/dt_tke, tketens(i,k) )
            END DO
         END DO

!        Jetzt enthaelt tketens die SQRT(2*TKE)-Tend. infolge Diffusion
!        und Drucktransport. Zu tketens wird nun an anderer Stelle noch
!        die Advektionstendenz hinzuaddiert und das ganze dann im
!        naechsten Zeitschritt fuer die SQRT(2*TKE)-Prognose benutzt.


! 14) Berechnung der turb. horizontalen TKE-Flussdichtedivergenzen
!     und Aufsummieren auf das SQRT(2*TKE)-Tendenzfeld:

!        ** wird erst spaeter eingefuehrt **


END SUBROUTINE turbdiff


!********************************************************************************

SUBROUTINE prep_impl_vert_diff ( i_st,i_en, k_tp, k_sf, lsflucond, &
                                 disc_mom, diff_mom, impl_mom, invs_mom )


INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for first  horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

LOGICAL, INTENT(IN) :: &
!
   lsflucond    !flux condition at the surface

REAL (KIND=ireals), INTENT(IN) :: &
!
   diff_mom(:,:), & !diffusion momentum (rho*K/dz)
   disc_mom(:,:)    !discretis momentum (rho*dz/dt) on variable levels

REAL (KIND=ireals), INTENT(INOUT) :: &
!
   impl_mom(:,:)    !in : discretis momentum (rho*dz/dt) on flux levels
                      !out: (negative) implicit part of diffusion momentum

REAL (KIND=ireals), INTENT(OUT) :: &
!
   invs_mom(:,:)    !inversion momentum

INTEGER (KIND=iintegers) :: &
!
   i,k, kp1,km1,  & !horizontal and vertical coordiante indices
   l, n,m           !used level index and number of used levels

   n=k_sf-k_tp-1      !counter of the last used level (starting at the column top)

!  Inverse momentum vector of uppermost level:

   l=1; k=k_tp+l; kp1=k+1

   DO i=i_st, i_en
      impl_mom(i,kp1)=-diff_mom(i,kp1) &
                        *MAX(MIN(diff_mom(i,kp1)/impl_mom(i,kp1), impl_s), impl_t) 
      invs_mom(i,k  )= z1/(disc_mom(i,k)-impl_mom(i,kp1))
   END DO
 
!  Note: Zero flux condition just below top level.

!  Inverse momentum vector of the other levels:

   IF (lsflucond) THEN
      m=n-1
   ELSE
      m=n
   END IF

   DO l=2, m
      k=k_tp+l; kp1=k+1; km1=k-1

      DO i=i_st, i_en
         impl_mom(i,kp1)=-diff_mom(i,kp1) &
                           *MAX(MIN(diff_mom(i,kp1)/impl_mom(i,kp1), impl_s), impl_t)
         invs_mom(i,k  )= z1/(disc_mom(i,k  )-impl_mom(i,k)-impl_mom(i,kp1) &
                               -invs_mom(i,km1)*impl_mom(i,k)**2)
      END DO

   END DO

   IF (lsflucond) THEN   
      l=n; k=k_tp+l; kp1=k+1; km1=k-1

      DO i=i_st, i_en
         impl_mom(i,kp1)= z0 !no  implicit surface weight
         invs_mom(i,k  )= z1/(disc_mom(i,k  )-impl_mom(i,k)                   &
                               -invs_mom(i,km1)*impl_mom(i,k)**2)
      END DO
   END IF

END SUBROUTINE prep_impl_vert_diff

!********************************************************************************

SUBROUTINE calc_impl_vert_diff ( i_st,i_en, k_tp, k_sf, &
!
   disc_mom, diff_mom, impl_mom, invs_mom, cur_prof, upd_prof, srf_flux)

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for first  horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

REAL (KIND=ireals), INTENT(IN) :: &
!
   diff_mom(:,:), & !diffusion momentum (rho*K/dz)
   disc_mom(:,:), & !discretis momentum (rho*dz/dt)
!
   invs_mom(:,:), & !inversion momentum
   impl_mom(:,:), & !implicit part of diffusion momentum
!
   cur_prof(:,:)    !current vertical variable profile

REAL (KIND=ireals), INTENT(OUT) :: &
!
   upd_prof(:,:)    !updated vertical variable profile

REAL (KIND=ireals), INTENT(OUT), OPTIONAL :: &
!
   srf_flux(:)      !effective upward vertical surface flux density

INTEGER (KIND=iintegers) :: &
!
   i,k, kp1,km1,  & !horizontal and vertical coordiante indices
   l,n,           & !used level index and number of used levels
   l0,l1            !alternating level indices

REAL (KIND=ireals) :: &
!
   expl_mom(i_en,0:1), & !explicit part of diffusion momentum
   expl_sfl(i_en),     & !explicit part of upward vertical surface flux density
   diff_tnd                             !diffusion tendency (times rho*dz)

   n=k_sf-k_tp-1      !counter of the last used level (starting at the column top)

!  Uppermost level:

   l=1; k=k_tp+l; kp1=k+1; l0=1; l1=0

   DO i=i_st, i_en
      expl_mom(i,l1)= diff_mom(i,kp1)+ impl_mom(i,kp1) 
      diff_tnd        = disc_mom(i,k  )* cur_prof(i,k  ) &
                       +expl_mom(i,l1 )*(cur_prof(i,kp1)-cur_prof(i,k))
      upd_prof(i,k) = diff_tnd         * invs_mom(i,k)            
   END DO

!  Note: Zero flux condition just below top level.

!  Intermediate levels:
 
   DO l=2, n-1
      k=k_tp+l; kp1=k+1; km1=k-1; l0=MOD(l,2); l1=MOD(l+1,2)

      DO i=i_st, i_en
         expl_mom(i,l1)= diff_mom(i,kp1)+ impl_mom(i,kp1)
         diff_tnd        = disc_mom(i,k  )* cur_prof(i,k  )                  &
                          +expl_mom(i,l0) *(cur_prof(i,km1)-cur_prof(i,k)) &
                          +expl_mom(i,l1) *(cur_prof(i,kp1)-cur_prof(i,k))
         upd_prof(i,k) =(diff_tnd         - impl_mom(i,k  )*upd_prof(i,km1))*invs_mom(i,k)
      END DO
   END DO

!  Surface level:

   l=n; k=k_tp+l; kp1=k+1; km1=k-1; l0=MOD(l,2); l1=MOD(l+1,2)

   DO i=i_st, i_en
      expl_mom(i,l1)= diff_mom(i,kp1)+ impl_mom(i,kp1)
      expl_sfl(i)   = expl_mom(i,l1 )*(cur_prof(i,kp1)-cur_prof(i,k))
      diff_tnd        = disc_mom(i,k  )* cur_prof(i,k  )                  &
                       +expl_mom(i,l0) *(cur_prof(i,km1)-cur_prof(i,k)) &
                       +expl_sfl(i)                                         &
                       -impl_mom(i,kp1)* cur_prof(i,kp1)
      upd_prof(i,k) =(diff_tnd          -impl_mom(i,k  )*upd_prof(i,km1))*invs_mom(i,k)
   END DO

   IF (PRESENT(srf_flux)) THEN
      DO i=i_st, i_en
         srf_flux(i)=expl_sfl(i)-impl_mom(i,kp1)*(cur_prof(i,kp1)-upd_prof(i,k))
      END DO
   END IF   
      
!  Note: Explicit flux condition at surface level in case of "impl_mom(:,:,kp1)=0".

!  Backsubstitution:

   DO l=n-1, 1, -1
      k=k_tp+l; kp1=k+1

      DO i=i_st, i_en
         upd_prof(i,k) = upd_prof(i,k)-impl_mom(i,kp1)*upd_prof(i,kp1)*invs_mom(i,k)
      END DO
   END DO

END SUBROUTINE calc_impl_vert_diff

!********************************************************************************

!------------------------------------------------------------------------------------
#ifdef SCLM
SUBROUTINE turb_stat (lsm, lsh, fm2, fh2, d_m, d_h, tls, tvs, tvt, grd, k)

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
    tls,    & !turbulent length scale   [m]
    tvs,    & !turbulent velocity scale [m/s]
    tvt,    & !total transport of turbulent velocity scale [m/s2]
    grd(:)    !vector of vertical gradients of diffused (conserved) variables [{unit}/m]

REAL (KIND=ireals) ::  &
!
   ts_d,                 & !dissipation time scale [s]
   tvs2,                 & !tvs**2 [m2/s2]
   tkm, tkh,             & !turbulent diffusion coefficient for momentum and scalars (heat) [m2/s]
   x1, x2, x3,           & !auxilary TKE source terme values [m2/s3]
   cvar(ndiff+1,ndiff+1)   !covariance matrix  [{unit1}*{unit2}]

   cvar(tet_l,tet_l)=d_h*tls*lsh*grd(tet_l)**2
   cvar(h2o_g,h2o_g)=d_h*tls*lsh*grd(h2o_g)**2
   cvar(tet_l,h2o_g)=d_h*tls*lsh*grd(tet_l)*grd(h2o_g)

   tkm=lsm*tvs; tkh=lsh*tvs

   cvar(u_m  ,w_m)=-tkm*grd(u_m)
   cvar(v_m  ,w_m)=-tkm*grd(v_m)
   cvar(tet_l,w_m)=-tkh*grd(tet_l)
   cvar(h2o_g,w_m)=-tkh*grd(h2o_g)

   TKE_SCLM%mod(k)%val=tvs         ; TKE_SCLM%mod(k)%vst=i_cal

   !Achtung: TKE%mod(k)%val zeigt z.Z. noch auf die alte TKE-Zeitstufe.
   !Somit wird also der alte TKE-Wert mit dem neuen ueberschrieben,
   !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
   !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

   tvs2=tvs**2
   ts_d=d_m*tls/tvs

   x1=cvar(u_m,w_m)*grd(u_m)
   x2=cvar(v_m,w_m)*grd(v_m)
   x3=-tkh*fh2

   BOYPR%mod(k)%val=x3             ; BOYPR%mod(k)%vst=i_cal
   SHRPR%mod(k)%val=-(x1+x2)       ; SHRPR%mod(k)%vst=i_cal
   DISSI%mod(k)%val=tvs2/ts_d      ; DISSI%mod(k)%vst=i_cal 
   TRANP%mod(k)%val=tvs*tvt        ; TRANP%mod(k)%vst=i_cal
                
   cvar(u_m,u_m)=z1d3*tvs2+ts_d*(-z4*x1+z2*x2-z2*x3)
   cvar(v_m,v_m)=z1d3*tvs2+ts_d*(+z2*x1-z4*x2-z2*x3)
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



!********************************************************************************
!********************************************************************************

SUBROUTINE turb_cloud (         &
   ie,  ke, ke1,  kcs,    kce,  &
   istart, iend,  kstart, kend, &
   prs, ps, rcld, t, qv, qc,    &
   clc, clwc                    )
   

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
!     itype_wcld = 1 :
!     The fractional cloud cover clc is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.
!     itype_wcld=2:
!     A Gaussion distribution is assumed for the saturation deficit
!     dq = qt - qs where qt = qv + ql is the total water content and
!     qs is the saturation specific humidity. Using the standard deviation
!     rcld of this distribution (on input) and the conservative grid-scale
!     quantities qt and tl (liquid water temperature), a corrected liquid
!     water content is determined which contains alse the contributions from
!     subgrid-scale clouds. A corresponding cloudiness is also calculated.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT(IN) :: &  ! indices used for allocation of arrays
  ie,                & ! number of grib points of horizontal dimensions
  ke, ke1,           & ! number of full and half level dimension
  kcs, kce,          & ! lower and upper level index bound for output variabls
  istart, iend, & ! zonal      start and end index
  kstart, kend    ! vertical   start and end index

! Array arguments with intent(in):

REAL (KIND=ireals), INTENT(IN) :: & !
  prs (:,:),    &    ! atmospheric pressure
  rcld(:,:),    &    ! standard deviation of saturation deficit
  ps  (:),        &    ! surface pressure
  t   (:,:),    &    ! temperature (main levels)
  qv  (:,:)          ! water vapour (")

REAL (KIND=ireals), OPTIONAL, INTENT(IN) :: & !
  qc  (:,:)          ! cloud water  (")

! Array arguments with intent(out):

REAL (KIND=ireals), INTENT(INOUT) :: &
  clc (:,kcs:),  & ! stratiform subgrid-scale cloud cover
  clwc(:,kcs:)     ! liquid water content of ""

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,k ! loop indices

REAL (KIND=ireals), PARAMETER :: &
  zsig_max = 1.0E-3_ireals,  & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_ireals,   & ! fraction of saturation specific humidity
  zuc      = 0.95_ireals       ! constant for critical relative humidity

REAL (KIND=ireals) :: &
  qt(ie), tl(ie), qs(ie), gam(ie) , &
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
           qt(i) = qc(i,k) +       qv(i,k) ! total water content
           tl(i) =  t(i,k) - lhocp*qc(i,k) ! liquid water temperature
        END DO 
     ELSE !qv and t already contain the conservation variables
        DO i = istart, iend
           qt(i) = qv(i,k)
           tl(i) =  t(i,k)
        END DO
     END IF

     DO i = istart, iend
!mod_2011/09/28: zpres=patm -> zpres=pdry {
        pres=(z1-qt(i))/(z1-rvd_m_o*qt(i))*prs(i,k)        ! part. pressure of dry air
        qs(i) = zqvap( zpsat_w( tl(i) ), pres )              ! saturation mixing ratio
!mod_2011/09/28: zpres=patm -> zpres=pdry }
       gam(i) = z1 / ( z1 + lhocp*zdqsdt( tl(i), qs(i) ) ) ! slope factor
     END DO

     DO i = istart, iend

        pres = prs(i,k)        ! pressure
        dq   = qt(i) - qs(i) ! saturation deficit

        IF ( itype_wcld .EQ. 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i)

          ! critical relative humidity
          uc = zuc - uc1 * zsigma * ( z1 - zsigma )  &
                         * ( z1 + uc2*(zsigma-0.5_ireals) )

          ! cloud cover
          clc(i,k) = MAX( z0,  &
                            MIN( z1, clc_diag * ((qt(i)/qs(i)-uc)/(ucl-uc))**2 ) )

          ! in-cloud water content
          ql = qs(i) * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_ireals ) THEN
            zclc1 = clc_diag * ( (z1-uc)/(ucl-uc) )**2
            ql    = ql + (gam(i)*dq-ql)*(clc(i,k)-zclc1)/(z1-zclc1)
          END IF
          clwc(i,k) = clc(i,k) * ql

        ELSEIF ( itype_wcld .EQ. 2 ) THEN

        ! Statistical calculation of cloud cover and cloud water content
        ! using the standard deviation of the saturation deficit

          sig = MIN ( zsig_max, rcld(i,k) )

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_ireals ) THEN
            clc(i,k)  = ABS ( (SIGN(z1,dq)+z1)*0.5_ireals )
            clwc(i,k) = clc(i,k) * gam(i) * dq
          ELSE
            q = dq/sig
            clc(i,k) = MIN ( z1, MAX ( z0, clc_diag * ( z1 + q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              clwc(i,k) = z0
            ELSEIF ( q >= zq_max ) THEN
              clwc(i,k) = gam(i) * dq
            ELSE
              clwc(i,k) = gam(i) * sig * ( q + q_crit ) &
                                           * ( q + zq_max ) / ( z2*( q_crit + zq_max) )
            ENDIF
          ENDIF

        ENDIF

     ENDDO

  ENDDO

END SUBROUTINE turb_cloud

!********************************************************************************

ELEMENTAL FUNCTION zpsat_w (ztemp)

  REAL (KIND=ireals), INTENT(IN) :: ztemp
  REAL (KIND=ireals) :: zpsat_w

  zpsat_w=b1*exp(b2w*(ztemp-b3)/(ztemp-b4w))

END FUNCTION zpsat_w

ELEMENTAL FUNCTION zqvap (zpvap, zpres)

  REAL (KIND=ireals), INTENT(IN) :: zpvap, &
                                    zpres !part. pressure of dry air
  REAL (KIND=ireals) :: zqvap

!mod_2011/09/28: zpres=patm -> zpres=pdry {
! zqvap=rdv*zpvap/(zpres-o_m_rdv*zpvap)
  zqvap=rdv*zpvap
  zqvap=zqvap/(zpres+zqvap)
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zqvap

ELEMENTAL FUNCTION zdqsdt (ztemp, zqsat)

  REAL (KIND=ireals), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=ireals) :: zdqsdt

!mod_2011/09/28: zpres=patm -> zpres=pdry {
! zdqsdt=b234w*(1.0_ireals+rvd_m_o*zqsat) &
!             *zqsat/(ztemp-b4w)**2
  zdqsdt=b234w*(z1-zqsat)*zqsat/(ztemp-b4w)**2
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zdqsdt

!********************************************************************************
!********************************************************************************

END MODULE src_turbdiff
