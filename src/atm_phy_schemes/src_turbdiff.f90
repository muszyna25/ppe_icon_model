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
!  Reorganization of the code for use in both models COMO and ICON with various formal modifications:
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

!-------------------------------------------------------------------------------

USE environment, ONLY :  &
!
    exchg_boundaries       ! performs the boundary exchange
#endif
!__COSMO__----------------------------------------------------------------------

#ifdef __COSMO__
USE data_turbdiff, ONLY : &
#endif
#ifdef __ICON__
USE mo_data_turbdiff, ONLY : &
#endif
!
!    get_param,    & !get turbulence parameters being extra converted
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
    b2w,          & ! over water (w) and ice (e)
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b234w,        & ! b2w * (b3 - b4w)
    uc1,          & ! variable for computing the rate of cloud cover in 
    uc2,          & ! the unsaturated case
    ucl             !               -- " --

!-------------------------------------------------------------------------------

#ifdef __COSMO__
USE data_turbdiff, ONLY : &
    lsso,         & ! SSO-Scheme is active
    lconv,        &   ! confection scheme is active
    a_hshr,       & ! factor for horizontal shear production of TKE
#endif
#ifdef __ICON__
USE mo_data_turbdiff ,ONLY: &
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
    a_stab,       & ! factor for stability correction of turbulent length scale
!
    tkhmin,       & ! minimal diffusion coefficients for heat
    tkmmin,       & ! minimal diffusion coefficients for momentum
!
    clc_diag,     & ! cloud cover at saturation in statistical cloud diagnostic
    q_crit,       & ! critical value for normalized over-saturation
    c_scld,       & ! factor for liquid water flux density in sub grid scale clouds
!
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
    llake           ! forecast with lake model FLake
#endif
#ifdef __ICON__
USE mo_lnd_nwp_config, ONLY : &
  &  lseaice,      & ! forecast with sea ice model
  &  llake           ! forecast with lake model FLake 
#endif


!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
!
    lsclm, i_cal, im, jm, &
!
    UUA, VVA, UWA, VWA, WWA, UST, TTA, TWA, SHF, LHF, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
#endif
!SCLM---------------------------------------------------------------------------

!-------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE
PUBLIC  :: init_canopy, organize_turbdiff

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
    z1d3=z1/z3 !,&
!    z2d3=z2/z3,&
!    z3d2=z3/z2

INTEGER (KIND=iintegers) :: &
    istat=0, ilocstat=0



LOGICAL :: lerror=.FALSE.

!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
  
!********************************************************************************
!********************************************************************************

!+ Module procedure init_canopy in "src_turbdiff" for initialization and allocation
!+ of special external parameters describing the surface canopy needed for for the
!+ description of surface-to-atmosphere transfer and within canopy diffusion:

SUBROUTINE init_canopy ( ie, je, ke, ke1, kcm, &
#ifdef __ICON__
           itype_tran, &
#endif
!
           istartpar, iendpar, jstartpar, jendpar, &
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
    je,  & ! number of grid points in meridional direction
    ke,  & ! number of main model levels (start index is 1)
    ke1, & ! number of half model levels (start index is 1) 
    istartpar, iendpar, & ! zonal      start and end index including the model boundary lines
    jstartpar, jendpar    ! meridional start and end index including the model boundary lines

INTEGER (KIND=iintegers), TARGET, INTENT(INOUT) :: &
!
    kcm             ! index of the lowest model layer higher than the canopy


#ifdef __ICON__
INTEGER (KIND=iintegers), INTENT(IN) :: itype_tran
#endif

REAL (KIND=ireals), DIMENSION(ie,je,ke1), OPTIONAL, INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=ireals), DIMENSION(ie,je), INTENT(IN) :: &
!
! External parameter fields:
! ----------------------------
    fr_land         ! land portion of a grid point area             ( 1 )

REAL (KIND=ireals), DIMENSION(ie,je), OPTIONAL, INTENT(IN) :: &
!
    plcov,        & ! fraction of plant cover                       ( 1 )
    lai             ! leaf area index                               ( 1 )

REAL (KIND=ireals), DIMENSION(ie,je), OPTIONAL, INTENT(INOUT) :: &
!
    sai,          & ! surface area index                            ( 1 )
    tai,          & ! transpiration area index                      ( 1 )
    eai             ! (evaporative) earth area index                ( 1 )

REAL (KIND=ireals), DIMENSION(ie,je), OPTIONAL, INTENT(INOUT) :: &
!
    h_can,        & ! hight of the vertically resolved canopy
    d_pat           ! horizontal pattern length scale

REAL (KIND=ireals), DIMENSION(ie,je,kcm:ke1), OPTIONAL, INTENT(INOUT) :: &
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
    i,  j,      & !  loop indices
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
     DO WHILE (MAXVAL( h_can - hhl(:,:,kcm) + hhl(:,:,ke1) ) .GT. z0)
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
     DO j=jstartpar, jendpar
     DO i=istartpar, iendpar
        IF (fr_land(i,j).LT.z1d2) THEN
           d_pat(i,j)=z0
        ELSE
           d_pat(i,j)=pat_len !should be a 2D external parameter field
        END IF
     END DO
     END DO
  END IF
! Effective values of the surface area indices:
  IF (PRESENT(sai) .AND. PRESENT(eai) .AND. PRESENT(tai) .AND. &
      PRESENT(lai) .AND. PRESENT(plcov)) THEN

     DO j=jstartpar, jendpar
     DO i=istartpar, iendpar
        IF (fr_land(i,j).LT.z1d2) THEN
           sai(i,j)=c_sea
        ELSE
           tai(i,j)=MAX( 1.0E-6_ireals, lai(i,j) )
        END IF
     END DO
     END DO

     IF (itype_tran.EQ.1) THEN
        DO j=jstartpar, jendpar
        DO i=istartpar, iendpar
           IF (fr_land(i,j).GE.z1d2) THEN
              sai(i,j)=tai(i,j)
              eai(i,j)=(z1-plcov(i,j))*sai(i,j)
              tai(i,j)=plcov(i,j)*tai(i,j)
            END IF
        END DO
        END DO
     ELSE
        DO j=jstartpar, jendpar
        DO i=istartpar, iendpar
           IF (fr_land(i,j).GE.z1d2) THEN
              tai(i,j)=plcov(i,j)*tai(i,j)  ! transpiration area index
              eai(i,j)=c_soil               ! evaporation area index
              sai(i,j)=c_lnd+tai(i,j)       ! surface area index

              fakt=EXP( e_surf*LOG( sai(i,j)) )/sai(i,j)
            ! effective area indeces by multiplication with the reduction factor fakt:
              sai(i,j)=fakt*sai(i,j)
              eai(i,j)=fakt*eai(i,j)
              tai(i,j)=fakt*tai(i,j)
           END IF
        END DO
        END DO
     END IF

  END IF

END SUBROUTINE init_canopy

!********************************************************************************
!********************************************************************************

!+ Module procedure organize_turbdiff in "src_turbdiff" for organising the calls 
!+ of surface-to-atmosphere trasnsfer and turbulent diffusion:

SUBROUTINE organize_turbdiff (action,iini,lstfnct, dt_var,dt_tke, nprv,ntur,ntim, &
!
          ie, je, ke, ke1, kcm, vst, &
          istart, iend, istartu, iendu, istartpar, iendpar, istartv, iendv, &
          jstart, jend, jstartu, jendu, jstartpar, jendpar, jstartv, jendv, &
#ifdef __ICON__
          isso, iconv, itype_tran, imode_tran, icldm_tran, imode_turb,      &
          icldm_turb, itype_sher, ltkesso, ltkecon, lexpcor, ltmpcor,       &
          lprfcor, lnonloc, lcpfluc, limpltkediff, itype_wcld, itype_synd,  &
#endif
          l_hori, &
#ifdef __COSMO__
          ntstep, &
          eddlon, eddlat, edadlat, acrlat, &
#endif
          hhl, dp0, &
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
          edr, tket_sso,             &
#ifdef __COSMO__
          tket_hshr,&
#endif
          tket_conv, &
!    
          u_tens, v_tens, t_tens, qv_tens, qc_tens, tketens, &
          qvt_diff, ut_sso, vt_sso, &
!    
          t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m, shfl_s, lhfl_s, &
!
          ierrstat, errormsg, eroutine)
 
!-------------------------------------------------------------------------------
! Description:
!
! Organizes the CALL of 'turbtran' and 'turbdiff' dependent on the parameters 'action'
! and 'iinit'.

! Method:
!
! All tendency parameters execp 'tketens' are OPTIONAL. If they are missing calculated 
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

CHARACTER (LEN=*), INTENT(IN) :: &
!
   action          !either 'only_tran', 'only_diff' or 'tran_diff'

LOGICAL, INTENT(IN) :: &
!
   lstfnct         !calculation of stability function required

#ifdef __ICON__
INTEGER, INTENT(IN) :: isso, iconv
INTEGER, INTENT(IN) :: itype_tran, imode_tran, icldm_tran, imode_turb
INTEGER, INTENT(IN) :: icldm_turb, itype_sher
LOGICAL, INTENT(IN) :: ltkesso, ltkecon, lexpcor, ltmpcor, lprfcor
LOGICAL, INTENT(IN) :: lnonloc, lcpfluc, limpltkediff
INTEGER, INTENT(IN) :: itype_wcld, itype_synd
#endif

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
    je,           & ! number of grid points in meridional direction
    ke,           & ! index of the lowest main model level
    ke1             ! index of the lowest model half level (=ke+1)

INTEGER (KIND=iintegers), INTENT(INOUT) :: &
!
    kcm             ! level index of the upper canopy bound

INTEGER (KIND=iintegers), INTENT(IN) :: &
!
    vst             ! velocity component staggering

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
    istartu,   iendu,   & ! start and end index for zonal      wind component
    istartv,   iendv,   & ! start and end index for meridional wind component
    istartpar, iendpar, & ! start and end index including model boundary lines
!
! Meridional direction:
    jstart,    jend,    & ! start and end index for the inner model domain
    jstartu,   jendu,   & ! start and end index for zonal      wind component
    jstartv,   jendv,   & ! start and end index for meridional wind component
    jstartpar, jendpar    ! start and end index including model boundary lines

REAL (KIND=ireals), INTENT(IN) :: &
!
! Constants related to the earth, the coordinate system 
! and the reference atmosphere:
! --------------------------------------------------------------------------
!
    l_hori          ! horizontal grid spacing (m)

#ifdef __COSMO__
! Time step indices:
! -----------------------------------------------------------------------

INTEGER (KIND=iintegers), OPTIONAL, INTENT(IN) :: &
!
    ntstep          ! current time step

REAL (KIND=ireals), OPTIONAL, INTENT(IN) :: &
!
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat         ! 1 / (r_earth * dlat)


REAL (KIND=ireals), DIMENSION(je,2), OPTIONAL, INTENT(IN) :: &
!
    acrlat          ! 1 / ( crlat * r_earth)
#endif

REAL (KIND=ireals), DIMENSION(ie,je,ke1), INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=ireals), DIMENSION(ie,je,ke), INTENT(IN) :: &
!
    dp0             ! pressure thickness of layer                   (pa )

REAL (KIND=ireals), DIMENSION(ie,je), INTENT(IN) :: &
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area             ( 1 )
    depth_lk,     & ! lake depth                                    ( m )
    sai             ! surface area index                            ( 1 )

REAL (KIND=ireals), DIMENSION(ie,je), TARGET, OPTIONAL, INTENT(IN) :: &
!
    d_pat           ! horizontal pattern length scale

REAL (KIND=ireals), DIMENSION(ie,je,kcm:ke1), TARGET, OPTIONAL, INTENT(IN) :: &
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml,        & ! effective drag coefficient of canopy elements
!                   ! smaller than the turbulent length scale            (1/m)
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)
REAL (KIND=ireals), DIMENSION(ie,je), INTENT(IN) :: &
!
! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------
!
    h_ice,        & ! ice thickness                                 (  m  )
    ps,           & ! surface pressure                              ( pa  )
    t_g,          & ! weighted surface temperature                  (  k  )
    qv_s            ! specific water vapor content on the surface   (kg/kg)

 REAL (KIND=ireals), DIMENSION(ie,je,ke), TARGET, INTENT(INOUT) :: &
!
! Atmospheric model variables:
! ---------------------------------
!
     u,           & ! zonal wind speed                              ( m/s )
     v,           & ! meridional wind speed                         ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=ireals), DIMENSION(ie,je,ke), TARGET, INTENT(IN) :: &
!
     prs            ! deviation from the reference pressure         ( pa  )

REAL (KIND=ireals), DIMENSION(ie,je,ke), TARGET, OPTIONAL, INTENT(IN) :: &
!
     rho,         & ! total density of air                          (kg/m3)
     epr            ! exner pressure                                 (1)

REAL (KIND=ireals), DIMENSION(ie,je,ke1), INTENT(IN) :: &
!
     w              ! vertical wind speed (defined on half levels)  ( m/s )

REAL (KIND=ireals), DIMENSION(ie,je), INTENT(INOUT) :: &
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

REAL (KIND=ireals), DIMENSION(ie,je,ke1,ntim), INTENT(INOUT) :: &
!
     tke             ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                     ! (defined on half levels)

REAL (KIND=ireals), DIMENSION(ie,je,2:ke1), INTENT(INOUT) :: &
!
     tkvm,         & ! turbulent diffusion coefficients for momentum (m/s2 )
     tkvh            ! turbulent diffusion coefficients for heat     (m/s2 )
                     ! and moisture

REAL (KIND=ireals), DIMENSION(ie,je,ke1), INTENT(INOUT) :: &
!
     rcld            ! standard deviation of the saturation deficit        
                     ! (as input and output)                             
                     ! fractional cloud cover (in turbdiff)            --

REAL (KIND=ireals), DIMENSION(ie,je,ke), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
! Tendency fields for the prognostic variables:
! -----------------------------------------------
!
     u_tens,       & ! u-tendency                                    ( m/s2)
     v_tens,       & ! v-tendency                                    ( m/s2)
     t_tens,       & ! t-tendency                                    ( K/s )
     qv_tens,      & ! qd-tendency                                   ( 1/s )
     qc_tens         ! qw-tendency                                   ( 1/s )

REAL (KIND=ireals), DIMENSION(ie,je,ke1), INTENT(INOUT) :: &
!
     tketens         ! tendency of SQRT(2*TKE)                       ( m/s2)

REAL (KIND=ireals), DIMENSION(ie,je,ke), OPTIONAL, INTENT(IN) :: &
!
     ut_sso,       & ! u-tendency due to the SSO-Scheme              ( 1/s )
     vt_sso          ! v-tendency due to the SSO-Scheme              ( 1/s )

REAL (KIND=ireals), DIMENSION(ie,je,ke1), OPTIONAL, INTENT(OUT) :: &
!
     edr             ! eddy dissipation rate of TKE (EDR)            (m2/s3)

REAL (KIND=ireals), DIMENSION(ie,je,ke), OPTIONAL, INTENT(IN) :: &
!
     tket_conv       ! TKE-tendency due to convective buoyancy       (m2/s3)

REAL (KIND=ireals), DIMENSION(ie,je,ke), OPTIONAL, INTENT(OUT) :: &
!
     tket_sso      ! TKE-tendency due to SSO wake production       (m2/s3)

#ifdef __COSMO__
REAL (KIND=ireals), DIMENSION(ie,je,ke), OPTIONAL, INTENT(OUT) :: &
     tket_hshr       ! TKE-tendency due to (sep.) horiz. shear       (m2/s3)
#endif

REAL (KIND=ireals), DIMENSION(ie,je,ke), OPTIONAL, INTENT(OUT) :: &
!
     qvt_diff        ! qd-tendency due to diffusion                  ( 1/s )

REAL (KIND=ireals), DIMENSION(ie,je), INTENT(OUT) :: &
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

 REAL (KIND=ireals), DIMENSION(ie,je), OPTIONAL, INTENT(INOUT) :: &
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2) (positive upward)
     lhfl_s          ! latent   heat flux at the surface             (W/m2) (positive upward)

INTEGER (KIND=iintegers), INTENT(INOUT) :: ierrstat

CHARACTER (LEN=*), INTENT(INOUT) :: eroutine
CHARACTER (LEN=*), INTENT(INOUT) :: errormsg

!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

REAL (KIND=ireals), DIMENSION(ie,je,0:7), TARGET :: &
     dd        !local derived turbulence parameter

INTEGER (KIND=iintegers) :: &
     nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes
     it_start      !Startindex der Iterationen

LOGICAL :: lini  !initialization of required

LOGICAL :: lsso, lconv


REAL (KIND=ireals), TARGET :: &
     c_tke,tet_g,c_g,rim, &
     d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
     a_3,a_5,a_6,b_1,b_2, &
     l_scal

!Declaration of statement functions:

REAL (KIND=ireals) :: &
     zexner, zpres               !Exner factor and its argument

!-----------------------------------------------------------------------
#ifdef __COSMO__
!INCLUDE 'statement_functs.incf'
 
!REAL (KIND=ireals) :: &
!     zpsat_w ,zqvap, zdqsdt, &   !thermodynamic statement functions and
!     ztemp, zpvap, zqsat         !their missing formal arguments 
 
!Definition of statement functions:
 
!saturation vapour pressure over water:
!zpsat_w(ztemp)=b1*exp(b2w*(ztemp-b3)/(ztemp-b4w))
 
!specific humidity:
!zqvap(zpvap,zpres)=rdv*zpvap/(zpres-o_m_rdv*zpvap)
 
!derivation of zqsat with respect to temperature:
!zdqsdt(ztemp,zqsat)=b234w*(1.0_ireals+rvd_m_o*zqsat) &
!                         *zqsat/(ztemp-b4w)**2
#endif
!__COSMO__----------------------------------------------------------------------

!Exner-factor:
 zexner(zpres)=(zpres/p0ref)**rdocp
 
!-------------------------------------------------------------------------------

 istat=0; ilocstat=0
 errormsg=''; eroutine='organize_turbdiff'; lerror=.FALSE.

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

 !Note:
 !A call with "iini=2" (at the first time step) provides the same result as a
 !  call with "iini=1" (before the time loop) followed by a second 
 !  call with "iini=0" (at the first time step).

 nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

 !Note:
 !It is also possible to use only one time level for TKE ("ntim=1" and thus "nprv=1=ntur").


!KF Translate INTEGER values into LOGIcAL ones
#ifdef __ICON__
    lsso         =(isso  .GT.0)
    lconv        =(iconv .GT.0)
#endif


!print *,"in organize_turbdiff"
 IF (action.EQ.'only_tran') THEN
    IF (.NOT.lerror) CALL turbtran(dt_tke)
 ELSEIF (action.EQ.'only_diff') THEN
    IF (.NOT.lerror) CALL turbdiff(dt_var,dt_tke,lstfnct)
 ELSEIF (action.EQ.'tran_diff') THEN
    IF (.NOT.lerror) CALL turbtran(dt_tke)
    IF (.NOT.lerror) CALL turbdiff(dt_var,dt_tke,lstfnct)
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

   INTEGER (KIND=iintegers) :: i,j !loop indices

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
!print *,"in turb_param l_hori=",l_hori," l_scal=",l_scal

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         dd(i,j,0)=d_m

         dd(i,j,1)=d_1
         dd(i,j,2)=d_2
         dd(i,j,3)=d_3
         dd(i,j,4)=d_4
         dd(i,j,5)=d_5
         dd(i,j,6)=d_6

         dd(i,j,7)=rim
      END DO
      END DO

END SUBROUTINE turb_param

!********************************************************************************

!+ Module procedure turbtran in "src_turbdiff" for computing the coefficients
!+ for turbulent transfer


SUBROUTINE turbtran(dt_tke)

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

! Array arguments with intent(in):

      REAL (KIND=ireals), INTENT (IN) :: &
              dt_tke !Laenge des Zeitschrittes fuer die TKE-Prognose    

! Local parameters:

! Local scalars:

!     Lokale Integer-Parameter:
     
      INTEGER (KIND=iintegers), PARAMETER :: &
!
!             Indexgrenzen:
!
              nred=4,      & !Anzahl der gegenueber vertikalen
                             !feuchtadiabatischen Verrueckungen
                             !invarianten Modellvariablen 
!
!             Zeiger fuer die Variablen :
!
              tet_l=1,     & !feucht-potentielle Temperatur
              h2o_g=2,     & !Gesamtwasseergehalt
              u_m=3,       & !zonale Geschw.komp. im Massenzentrum
              v_m=4          !meridionale  ,,      ,,     ,, 
 
!     Beachte:tet_l,h2o_g muessen in [1,ninv]; u_m,v_m in [ninv+1,nred],
!             liegen!

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &
!
              i_st,i_en,    & !Start- und Endindices in zonale Richtung
              j_st,j_en,    & !Start- und Endindices in merid. Richtung
              i,j,k,        & !Diskretis.index fuer lamda, phi, sigma
              ii,jj,        & !Indices fuer diverse Schleifen
              it_durch        !Durchgangsindex der Iterationen

!     Lokale real Variablen:
#ifdef SCLM
      REAL (KIND=ireals) ::  teta,qvap
#endif
      REAL (KIND=ireals) :: &
!
!          Hilfsvariablen:
!
           wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,             & !  ,,         ,,     ,,      Faktoren
!
!     Platzh. fuer Temperaturen, Druck, Wasserdampf-Mischungsverhaeltniss und
!     dessen Ableitung nach der Temperatur:
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
           lo_ice(ie,je)   ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=ireals) ::  &
!
        l_tur_z0 (ie,je),    &
!
!US for a better vectorization (hopefully):
!
        h_top_2d (ie,je),    &
        h_atm_2d (ie,je),    &
!
        z0m_2d   (ie,je),    &
        z0d_2d   (ie,je),    &
        z2m_2d   (ie,je),    &
        z10m_2d  (ie,je),    &
!
        hk_2d    (ie,je),    &
        hk1_2d   (ie,je),    &
        h_can_2d (ie,je),    &
!
        rat_m_2d (ie,je),    &
        rat_h_2d (ie,je),    &
        fac_h_2d (ie,je),    &
        fm2_2d   (ie,je),    &
        fh2_2d   (ie,je),    &
        frc_2d   (ie,je),    &
!
        vel1_2d  (ie,je),    &
        vel2_2d  (ie,je),    &
        vel_2d   (ie,je),    &
        ta_2d    (ie,je),    &
        qda_2d   (ie,je),    &
        ts_2d    (ie,je),    &
        qds_2d   (ie,je),    &
        rho_2d   (ie,je),    &
        prs_2d   (ie,je),    &
!
        dz_0a_m  (ie,je),    &
        dz_0a_h  (ie,je),    &
        dz_sa_h  (ie,je),    &
        dz_s0_h  (ie,je)

INTEGER (KIND=iintegers) ::  &
        k_2d     (ie,je)

      REAL (KIND=ireals) :: &
!
!     Vertikale Gradienten verschiedener thermodyn. Variablen:
!
           grad(ie,je,nred), &                          
!
!     Zwischenspeicher fuer Wolkenwasser und Bedeckungsgrad:
!
           clc(ie,je,ke:ke),clcw(ie,je,ke:ke)

!---- End of header ---------------------------------------------------------

! 1)  Vorbereitungen:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbtran'; lerror=.FALSE.

!     Unterste Hauptflaeche halbiert den Abstand zur
!     untersten Nebenflaeche:

      xf=z2 
      wf=xf-z1
      xf=z1/xf

!     Festlegung der horiz. Schleifengrenzen: 

      i_st=istartpar 
      i_en=iendpar

      j_st=jstartpar
      j_en=jendpar

      IF (istat /= 0) THEN
         ierrstat = 1004
         errormsg= &
         'ERROR *** Allocation of space for meteofields failed ***'
         lerror=.TRUE.; RETURN
      ENDIF

!     Berechnung abgeleiteter Parameter:

      CALL turb_param

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren: 

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

      DO j=j_st,j_en
      DO i=i_st,i_en

         IF (fr_land(i,j).LT.z1d2) THEN
            ! Water point.
            IF (.NOT. lseaice) THEN
              ! Sea ice model is not used.
              ! Ice surface if SST is less than the salt water freezing temperature.
              lo_ice(i,j) = t_g(i,j) < t0_melt + zt_ice
            ELSE
              ! Sea ice model is used.
              ! Ice surface if ice is present.
              lo_ice(i,j) = h_ice(i,j) > z0
            END IF
            IF (llake) THEN
              ! Lake model is used.
              ! Ice surface if this is a lake point AND ice is present.
              IF ((depth_lk(i,j) > z0) .AND. (h_ice(i,j) >= h_Ice_min_flk)) &
              lo_ice(i,j) = .TRUE.
            END IF
         END IF

      END DO
      END DO

      IF (lini) THEN

         DO j=j_st,j_en
         DO i=i_st,i_en

            IF (fr_land(i,j).LT.z1d2) THEN

!              Ueber Wasserpunkten:

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice

               IF ( lo_ice(i,j) ) THEN

!                 Bei Eisdecke: 

                  gz0(i,j)=grav*z0_ice
               ELSE

!                 Bei von Schubspannung abhaengiger Wellenhoehe:

!                 Einfachste Schaetzung der Schubspannung als Impusls-
!                 flussdichte durch die Nebenflaeche ke mit Hilfe
!                 einer diagnostischen TKE ohne Beruecksichtigung von
!                 Feuchte-Effekten und mit neuchtralen Stabilitaets-
!                 funktion und Anwendung der Charnockflormel:

                  ii=MAX( i-vst, 1 )
                  jj=MAX( j-vst, 1 )

                  l_turb=hhl(i,j,ke)-hhl(i,j,ke1)
                  l_turb=akt*MAX( len_min, l_turb/(z1+l_turb/l_scal) )
!test
!                 l_turb=akt*MIN( l_scal, hhl(i,j,ke)-hhl(i,j,ke1) )
!test

                  dh=hhl(i,j,ke-1)-hhl(i,j,ke1)

                  vel1=(u(i,j,ke-1)+u(ii,j,ke-1))
                  vel2=(u(i,j,ke)+u(ii,j,ke))
                  grad(i,j,u_m)=(vel1-vel2)/dh

                  vel1=(v(i,j,ke-1)+v(i,jj,ke-1))
                  vel2=(v(i,j,ke)+v(i,jj,ke))
                  grad(i,j,v_m)=(vel1-vel2)/dh

                  grad(i,j,tet_l)=z2*(t(i,j,ke-1)-t(i,j,ke))/dh &
                                 +tet_g

                  fm2=grad(i,j,u_m)**2+grad(i,j,v_m)**2
                  fh2=grav*grad(i,j,tet_l)/t(i,j,ke)

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

                  tke(i,j,ke,nvor)=q0

                  wert=lm*q0*SQRT(fm2)
                  gz0(i,j)=MAX( grav*len_min, alpha0*wert )
!                 gz0(i,j)=MAX( grav*len_min, &
!                               alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF

            tkvm(i,j,ke)=con_m
            tkvh(i,j,ke)=con_h
            tkvm(i,j,ke1)=con_m
            tkvh(i,j,ke1)=con_h

         END DO
         END DO
      END IF


! 3)  Diagnose des Fluessigwassergehaltes und des Bedeckungsgrades
!     in der unterstersten Modellschicht:

      CALL turb_cloud (                       &
#ifdef __ICON__
           itype_wcld,                        & 
#endif
           ie, je, ke,ke1, ke,ke,             &
           i_st,i_en, j_st,j_en,       ke,ke, &
           prs, ps, rcld, t, qv, qc,          &
           clc, clcw                            )

      IF (icldm_tran.EQ.0) THEN
         DO j=j_st,j_en
         DO i=i_st,i_en
            clc(i,j,ke)=z0
         END DO 
         END DO 
       ELSEIF (icldm_tran.EQ.1) THEN
         DO j=j_st,j_en
         DO i=i_st,i_en
            IF ( qc(i,j,ke) > z0 ) THEN
               clc(i,j,ke)=z1
            ELSE
               clc(i,j,ke)=z0
            ENDIF
         END DO   
         END DO   
      ENDIF

! 4)  Berechnung der Transferkoeffizienten:

      DO it_durch=it_start, it_end !Iterationen

         DO j=j_st,j_en
         DO i=i_st,i_en

!           Berechnung der benoetigten Laengenskalen:

            ! Dicke der Modell-Prandtl-Schicht
            h_top_2d(i,j) = hhl(i,j,ke)-hhl(i,j,ke1)
            h_atm_2d(i,j) = h_top_2d(i,j)*xf        

            ! Hoehe des 2m- und 10m-Niveaus
            h_2m  = z2
            h_10m = z10

            ! Rauhigkeitslaenge
            z0m_2d(i,j) = gz0(i,j)/grav
            z_surf      = z0m_2d(i,j)/sai(i,j)

            ! turbulente Laengenskala 
            l_tur_z0(i,j) = akt*z0m_2d(i,j) 

            ! turbulente Distanz auf der untersten Hauptflaeche
            a_atm = h_atm_2d(i,j)+z0m_2d(i,j) 

            ! turbulente Distanz auf der untersten Nebenflaeche
            a_top = h_top_2d(i,j)+z0m_2d(i,j) 

!           Laminare Korrektur der Diffusionskoeffizienten:
            tkvm(i,j,ke1)=MAX( con_m, tkvm(i,j,ke1) )
            tkvh(i,j,ke1)=MAX( con_h, tkvh(i,j,ke1) )

            fakt=z1+(z1-REAL(NINT(fr_land(i,j)),ireals))*(rat_sea-z1)

            rat_m=tkvm(i,j,ke1)/con_m
            rat_h=tkvh(i,j,ke1)/con_h

!           Berechnung der effektiven Widerstandslaenge der L-Schicht:
            dz_sg_m=rlam_mom*z_surf
            dz_sg_h=fakt*rlam_heat*z_surf*(rat_h/rat_m)

!           Berechnung weiterer effektiver Widerstandslaengen fuer Skalare:

!           Bestandesschicht ohne lam. Grenzschicht:
            dz_g0_h=z_surf*LOG(rat_m)

!           Bestandesschicht inclusive lam. Grenzschicht:
            dz_s0_h(i,j)=dz_sg_h+dz_g0_h

!           Berechnung der effektiven Bestandeshoehe:
            IF (dz_sg_h.eq.z0) THEN
               h_can_2d(i,j)=rat_can*z0m_2d(i,j)
            ELSE
               h_can_2d(i,j)=rat_can*dz_s0_h(i,j)*LOG(dz_s0_h(i,j)/dz_sg_h)
            END IF

!           Berechnung weiterer effektiver Widerstandslaengen fuer Impuls:

!           Widerstandslaengen inclusive lam. Grenzschicht:
            wert=z1d2*dz_sg_m 
            dz_s0_m=wert+SQRT(wert**2+h_can_2d(i,j)*dz_sg_m)

!           Widerstandslaengen ohne lam. Grenzschicht:
            dz_g0_m=dz_s0_m-dz_sg_m

!           der turb. Prandtl-Schicht:

            rat_m=(tkvm(i,j,ke)*z0m_2d(i,j))/(tkvm(i,j,ke1)*a_top)
            rat_m_2d(i,j)=MIN( z2, MAX( z1d2, rat_m ) )
            
            fac_m=(rat_m_2d(i,j)-z1)*z0m_2d(i,j)/h_top_2d(i,j)
            IF (fac_m.EQ.z1) THEN
               dz_0a_m(i,j)=z0m_2d(i,j)*h_atm_2d(i,j)/a_atm
            ELSE
               dz_0a_m(i,j)=z0m_2d(i,j)*LOG(a_atm/(z0m_2d(i,j)+fac_m*h_atm_2d(i,j)))/(z1-fac_m)
            END IF

            rat_h=(tkvh(i,j,ke)*z0m_2d(i,j))/(tkvh(i,j,ke1)*a_top)
            rat_h_2d(i,j)=MIN( z2, MAX( z1d2, rat_h ) )

            fac_h_2d(i,j)=(rat_h_2d(i,j)-z1)*z0m_2d(i,j)/h_top_2d(i,j)
            IF (fac_h_2d(i,j).EQ.z1) THEN
               dz_0a_h(i,j)=z0m_2d(i,j)*h_atm_2d(i,j)/a_atm
            ELSE
               dz_0a_h(i,j)=z0m_2d(i,j)*LOG(a_atm/(z0m_2d(i,j)+fac_h_2d(i,j)*h_atm_2d(i,j))) &
                                      /(z1-fac_h_2d(i,j))
            END IF

!           von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):
            dz_sa_m      = dz_s0_m      + dz_0a_m(i,j)
            dz_sa_h(i,j) = dz_s0_h(i,j) + dz_0a_h(i,j)

!           Reduktionsfaktoren fuer die Bestandesschicht
!           incl. lam. Grenzschicht:

            tfm(i,j)=dz_0a_m(i,j)/dz_sa_m      
            tfh(i,j)=dz_0a_h(i,j)/dz_sa_h(i,j)

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den
!           Faktor 'rat_lam' gegenueber fuehlbarer Waerme vergroesserten
!           laminaren Transpostwiderstandes:

            tfv(i,j)=z1/(z1+(rat_lam-z1)*dz_sg_h/dz_sa_h(i,j))
         END DO
         END DO      

         DO j=j_st,j_en
         DO i=i_st,i_en
!           Berechnung der Windgeschwindigkeiten an den Massepunkten:

            ii=MAX( i-vst, 1 )
            jj=MAX( j-vst, 1 )

            vel1_2d(i,j)=(u(i,j,ke)+u(ii,j,ke))*z1d2
            vel2_2d(i,j)=(v(i,j,ke)+v(i,jj,ke))*z1d2

!           Berechnung der korrigierten Prandtl-Schicht-Werte:

            IF (lprfcor) THEN
               len1=z2*h_top_2d(i,j)
               len2=(h_top_2d(i,j)-h_atm_2d(i,j))**2 &
                 /((hhl(i,j,ke-1)+hhl(i,j,ke))*z1d2-hhl(i,j,ke1)-h_atm_2d(i,j))
               lm=len1-tfm(i,j)*h_atm_2d(i,j)-len2
               lh=len1-tfh(i,j)*h_atm_2d(i,j)-len2

               velh=(u(i,j,ke-1)+u(ii,j,ke-1))*z1d2
               vel1_2d(i,j)=(len1*vel1_2d(i,j)-len2*velh)/lm
               velh=(v(i,j,ke-1)+v(i,jj,ke-1))*z1d2
               vel2_2d(i,j)=(len1*vel2_2d(i,j)-len2*velh)/lm

               ta_2d(i,j)=(len1*t(i,j,ke)-h_atm_2d(i,j)*tfh(i,j)*t_g(i,j) &
                                    -len2*t(i,j,ke-1))/lh
               qda_2d(i,j)=(len1*qv(i,j,ke)-h_atm_2d(i,j)*tfh(i,j)*qv_s(i,j) &
                                      -len2*qv(i,j,ke-1))/lh

               vel_2d(i,j)=MAX( vel_min, SQRT(vel1_2d(i,j)**2+vel2_2d(i,j)**2) )
            ELSE
               ta_2d(i,j)=t(i,j,ke)
               qda_2d(i,j)=qv(i,j,ke)

               vel_2d(i,j)=MAX( vel_min, SQRT(vel1_2d(i,j)**2+vel2_2d(i,j)**2) )
            END IF

         END DO   
         END DO 

         DO j=j_st,j_en
         DO i=i_st,i_en
!           Berechnung der unteren Randwerte der Prandtl-Schicht:
             ts_2d(i,j)= ta_2d(i,j)*(z1-tfh(i,j))+ t_g(i,j)*tfh(i,j)
            qds_2d(i,j)=qda_2d(i,j)*(z1-tfh(i,j))+qv_s(i,j)*tfh(i,j)
         END DO
         END DO

         IF (it_durch.LE.it_end) THEN
            !Strange to say this seemingly needless condition makes the NEC-compiler
            !to inline the CALLs of thermodynamic FUNCTIONs and finally to do the
            !vectorization of the follwing loop!

         DO j=j_st,j_en
         DO i=i_st,i_en
!           Berechnung von thermodynamischen Groessen:

            exnr   =zexner(ps(i,j)) !Exner-Faktor
            virt   =z1/(z1+rvd_m_o*qds_2d(i,j)-clcw(i,j,ke)) 
                                    !rezip. virtueller Faktor

            qsat_dT=zdqsdt( ts_2d(i,j), zqvap( zpsat_w( ts_2d(i,j) ), ps(i,j) ) )

            rho_2d(i,j)=virt*ps(i,j)/(r_d*ts_2d(i,j)) !Luftdichte

!           Berechn. der thermodyn. Hilfsfaktoren:
            fakt=lhocp/ts_2d(i,j)-(z1+rvd_m_o)*virt
            wert=clc(i,j,ke)*fakt/(z1+qsat_dT*lhocp)

!           Berechn. der thermodyn. Faktoren A_tet und A_q:
            a_tet=z1/ts_2d(i,j)-wert*qsat_dT
            a_vap=rvd_m_o*virt+wert

!           Berechnung der benoetigten vertikalen Gradienten und
!           abgeleiteter Groessen:

            grad(i,j,u_m)=tfm(i,j)*vel1_2d(i,j)/dz_0a_m(i,j)
            grad(i,j,v_m)=tfm(i,j)*vel2_2d(i,j)/dz_0a_m(i,j)

!           Beachte: Fuer die Windkompnenten wurden die unteren
!           Randwerte der Prandtl-Schicht nicht berechnet, daher
!           muss der Faktor tfm hier bei grad(u_m) und grad(v_m)
!           auftauchen!

            grad(i,j,tet_l)=(ta_2d(i,j)-ts_2d(i,j))/dz_0a_h(i,j) &
                           +tet_g*(z1+clcw(i,j,ke)*lhocp*virt/ts_2d(i,j))
            grad(i,j,h2o_g)=(qda_2d(i,j)-qds_2d(i,j))/dz_0a_h(i,j)

            fm2_2d(i,j)=grad(i,j,u_m)**2+grad(i,j,v_m)**2
            fh2_2d(i,j)=grav*(a_tet*grad(i,j,tet_l)+a_vap*grad(i,j,h2o_g))

            !Beachte:
            !'grad(i,j,tet_l)' ist bislang der Temperaturgradient und haette zuvor durch
            !den Exnerfaktor dividiert werden muessen, so wie 'a_tet' zuvor damit
            !haette multipliziert werden muessen. Der Exnerfaktor kuerzt sich aber
            !in 'fh2' wieder heraus. Im folgenden wird aber der Teta-Gradient benoetigt:

            grad(i,j,tet_l)=grad(i,j,tet_l)/exnr
         END DO
         END DO      

         END IF

!        Berechnung der atmosphaerischen Forcing-Funktion:

         IF (it_durch.EQ.it_start .AND. lini) THEN !Startinitialisierung

            DO j=j_st,j_en
            DO i=i_st,i_en

               IF (fh2_2d(i,j).GE.(z1-rim)*fm2_2d(i,j)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=l_tur_z0(i,j)*(b_2-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2_2d(i,j)/(fm2_2d(i,j)-fh2_2d(i,j))
                  lm=l_tur_z0(i,j)*(b_2-(a_6+a_3)*fakt)
                  lh=l_tur_z0(i,j)*(b_1-a_5*fakt)
               END IF

               val1=lm*fm2_2d(i,j); val2=lh*fh2_2d(i,j)  
               frc_2d(i,j)=MAX( val1-val2, rim*val1 )

               tke(i,j,ke1,nvor)=MAX( SQRT(d_m*l_tur_z0(i,j)*frc_2d(i,j)), vel_min )
            END DO
            END DO       

         ELSE ! mit Hilfe der vorhergehenden TKE-Werte 

            DO j=j_st,j_en
            DO i=i_st,i_en
               val1=tkvm(i,j,ke1)*fm2_2d(i,j); val2=tkvh(i,j,ke1)*fh2_2d(i,j)
               frc_2d(i,j)=MAX( val1-val2, rim*val1 )/tke(i,j,ke1,nvor)
            END DO
            END DO       
         END IF  

! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes:

         fakt=d_4/(z1-c_g)
         w1=tkesmot; w2=z1-tkesmot

         DO j=j_st,j_en
         DO i=i_st,i_en

            wert=fakt*l_tur_z0(i,j)

            q0=tke(i,j,ke1,nvor)

            frc_2d(i,j)=MIN( q0**2/wert, frc_2d(i,j) )

            IF (imode_tran.EQ.1) THEN
               q3=SQRT(d_m*l_tur_z0(i,j)*frc_2d(i,j))
            ELSE
               q1=d_m*l_tur_z0(i,j)/dt_tke
               q2=q0+frc_2d(i,j)*dt_tke
               q3=q1*(sqrt(z1+z4*q2/q1)-z1)*z1d2
            END IF

            q1=SQRT(wert*frc_2d(i,j))
            tke(i,j,ke1,ntur)=MAX( vel_min, q1, w1*q0 + w2*q3 )

!           Was die diversen Beschraenkungen durch MAX- und MIN-
!           Funktionen betrifft, siehe bei den entsprechenden Stellen
!           in 'turbdiff.incf'.

!           frc_2d wird noch zusaetzlich nach oben beschraenkt, damit
!           dieser Wert im Zusammenhang mit q0 keine Singularitaet
!           in der Stabilitaetsfunktion verursachen wuerde. Dies traegt zur
!           numerischen Stabilitaet der iterativ ueber den Zeitindex erfolgenden
!           Bestimmung der Vertikalgradienten bei.

         END DO
         END DO

         IF (PRESENT(edr)) THEN
            DO j=j_st,j_en
            DO i=i_st,i_en
               edr(i,j,ke1)=tke(i,j,ke1,ntur)**3/(d_m*l_tur_z0(i,j))
            END DO
            END DO
         END IF   

! 4g)    Bestimmung der neuen stabilitaetsabh. Laengenskalen:

         CALL stab_funct(sm=tkvm(:,:,ke1), sh=tkvh(:,:,ke1), fm2=fm2_2d, fh2=fh2_2d, &
                         frc=frc_2d, tvs=tke(:,:,ke1,ntur), tls=l_tur_z0, &
                         i_st=i_st,i_en=i_en, j_st=j_st,j_en=j_en)

         DO j=j_st,j_en
         DO i=i_st,i_en

! 4h)       Bestimmung der durch Wirkung der L-Schicht
!           korrigierten Diffusionskoeffizienten 
!           und der zugehoerigen Transferkoeffizienten:

!           unkorrigierte Diffusionskoeffizienten:

            fakt=l_tur_z0(i,j)*tke(i,j,ke1,ntur)
            tkvm(i,j,ke1)=fakt*tkvm(i,j,ke1)
            tkvh(i,j,ke1)=fakt*tkvh(i,j,ke1)

!           Belegung der Felder fuer die Transferkoeffizienten:

            tcm(i,j)=tkvm(i,j,ke1)*tfm(i,j)/(dz_0a_m(i,j)*vel_2d(i,j))
            tch(i,j)=tkvh(i,j,ke1)*tfh(i,j)/(dz_0a_h(i,j)*vel_2d(i,j))

! 4i)       Diagnose von gz0 (fuer den naechsten Zeitschritt) 
!           ueber Wasserflaechen mit der (angepassten) Charnock-Formel
!           und Einschraenkung von z0m_dia ueber Land:
 
            IF (fr_land(i,j).LT.z1d2) THEN

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice
               IF ( lo_ice(i,j) ) THEN
                  ! Ice-covered grid box
                  gz0(i,j)=grav*z0_ice
               ELSE
                  velh=(tke(i,j,ke1,ntur)+tke(i,j,ke,ntur))*z1d2
                  wert=tcm(i,j)*vel_2d(i,j)*SQRT(vel_2d(i,j)**2+velh**2)
                  gz0(i,j)=MAX( grav*len_min, alpha0*wert )
!                 gz0(i,j)=MAX( grav*len_min, alpha0*wert+grav*alpha1*con_m/SQRT(wert) )
               END IF
               !Ueber See gibt es keinen synoptischen Garten
               z0d_2d(i,j)=z0m_2d(i,j)
            ELSE
               !Die Rauhigkeitslaenge einer SYNOP Station soll immer
               !kleiner als 10m bleiben:
               z0d_2d(i,j)=MIN( z10, z0m_dia )
            END IF
         END DO
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

      DO j=j_st,j_en
      DO i=i_st,i_en

         IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile
           z2m_2d (i,j) = z2-h_can_2d(i,j) !2m ueber dem Bodenniveau des Bestandes
           z10m_2d(i,j) = z10-z0m_2d(i,j) !Hoehe in der turbulente Distanz 10m betraegt
         ELSE !using only a logarithmic profile above a SYNOP lawn
           z2m_2d (i,j) = h_2m
           z10m_2d(i,j) = h_10m
         END IF   

         !Erste Belegung zweier benachbarter Modellniveaus:            

         hk_2d(i,j)=h_atm_2d(i,j)
         hk1_2d(i,j)=z0
         k_2d(i,j)=ke

      ENDDO
      ENDDO

!     Diagnose der 2m-Groessen:

      CALL diag_level(i_st, i_en, j_st, j_en, z2m_2d, k_2d, hk_2d, hk1_2d)

      IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile

         DO j=j_st,j_en
         DO i=i_st,i_en
            IF (k_2d(i,j).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               prs_2d(i,j)=ps(i,j)

               IF (z2m_2d(i,j).LT.z0) THEN
!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

                  fakt=dz_s0_h(i,j)/dz_sa_h(i,j)
                  IF (dz_s0_h(i,j) >= epsi) THEN
                    fakt=MIN( z1, MAX( z0, fakt*EXP(z2m_2d(i,j)/dz_s0_h(i,j)) ) )
                  ELSE
! mr: Der Exponent x geht gegen -inf., und 0<=fakt<epsi, so dass fakt*EXP(x) gegen Null geht
                    fakt=z0
                  ENDIF
               ELSE
!                 2m-Niveau liegt innerhalb der Modell_Prandtl-Schicht
!                 mit logarithmischen Vertikalprofilen:

                  IF (ABS(z1-fac_h_2d(i,j)) < epsi ) THEN
                     wert=z0m_2d(i,j)*z2m_2d(i,j)/(z2m_2d(i,j)+z0m_2d(i,j))
                  ELSE
                     wert=z0m_2d(i,j)*LOG((z2m_2d(i,j)+z0m_2d(i,j))/            &
                         (z0m_2d(i,j)+fac_h_2d(i,j)*z2m_2d(i,j)))/(z1-fac_h_2d(i,j))
                  END IF
                  fakt=(dz_s0_h(i,j)+wert)/dz_sa_h(i,j)
               END IF
               t_2m(i,j) = t_g(i,j) + (ta_2d(i,j)-t_g(i,j))*fakt &
                         + tet_g*( (h_atm_2d(i,j)+h_can_2d(i,j) )*fakt-h_2m ) 
               qv_2m(i,j)= qv_s(i,j) + (qda_2d(i,j)-qv_s(i,j))*fakt
            END IF
         END DO   
         END DO   

      ELSE !using only a logarithmic profile above a SYNOP lawn
               
         DO j=j_st,j_en
         DO i=i_st,i_en
            IF (k_2d(i,j).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               prs_2d(i,j)=ps(i,j)

               z0d=z0d_2d(i,j)
               a_atm=h_atm_2d(i,j)+z0d
               a_2m=h_2m+z0d

!              Dimensionsloser Widerstand des Rauhigkeitsbestandes der SYNOP-Wiese,
!              wobei die turbulente Geschwindigkeitsskala und der Oberflaechenindex
!              gegenueber dem mittleren Rauhigkeitsbestand gleich bleiben:

               fr_sd_h=MAX( z0, dz_s0_h(i,j)/z0m_2d(i,j)+LOG(z0d/z0m_2d(i,j)) )

               fac_h=(rat_h_2d(i,j)-z1)*z0d/h_top_2d(i,j)
               IF (fac_h.EQ.z1) THEN
                  val1=fr_sd_h+h_2m/a_2m
                  val2=fr_sd_h+h_atm_2d(i,j)/a_atm
               ELSE
                  val1=fr_sd_h+LOG(a_2m/(z0d+fac_h*h_2m))/(z1-fac_h)
                  val2=fr_sd_h+LOG(a_atm/(z0d+fac_h*h_atm_2d(i,j)))/(z1-fac_h)
               END IF

               fakt=val1/val2
               t_2m(i,j) = t_g(i,j) + (ta_2d(i,j)-t_g(i,j))*fakt &
                         + tet_g*(h_atm_2d(i,j)*fakt-h_2m)
               qv_2m(i,j)= qv_s(i,j) + (qda_2d(i,j)-qv_s(i,j))*fakt
            END IF   
         END DO   
         END DO   

      END IF   

      DO j=j_st,j_en
      DO i=i_st,i_en
         IF (k_2d(i,j).GT.ke) THEN

!           2m-Niveau liegt oberhalb der untersten Hauptflaeche:

            k = k_2d(i,j)
            prs_2d(i,j)=prs(i,j,k+1)

            IF (itype_synd.EQ.3) THEN !using an exponetial rougness layer profile
               val1=(z2m_2d(i,j)+z0m_2d(i,j))/(hk1_2d(i,j)+z0m_2d(i,j))
               val2=(hk_2d(i,j)+z0m_2d(i,j))/(hk1_2d(i,j)+z0m_2d(i,j))
            ELSE !using only a logarithmic profile above a SYNOP lawn
               val1=(h_2m+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))
               val2=(hk_2d(i,j)+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))
            END IF

            fakt=LOG(val1)/LOG(val2)

            t_2m(i,j) =  t(i,j,k+1) + fakt*( t(i,j,k)- t(i,j,k+1)) &
                       + tet_g*( (hk_2d(i,j)-hk1_2d(i,j))*fakt - (h_2m-hk1_2d(i,j)) )
            qv_2m(i,j) = qv(i,j,k+1) + fakt*(qv(i,j,k)-qv(i,j,k+1))

         END IF
      END DO    
      END DO

      DO j=j_st,j_en
      DO i=i_st,i_en

!        Berechnung der zugehoerigen Feuchtewerte:

         wert=grav*(z2m_2d(i,j)-hk1_2d(i,j)) &
             /(r_d*t_2m(i,j)*(z1+rvd_m_o*qv_2m(i,j))) 
         patm=prs_2d(i,j)*EXP(-wert)*qv_2m(i,j) &
             /(rdv+(z1-rdv)*qv_2m(i,j))               !Wasserdampfdruck

         fakt=patm/zpsat_w( t_2m(i,j) )
         rh_2m(i,j)=100.0_ireals*MIN( fakt, z1 )      !relative Feuchte

         wert=LOG(patm/b1)
         td_2m(i,j)=MIN( (b2w*b3-b4w*wert) &
                   /(b2w-wert), t_2m(i,j) )           !Taupunktstemperatur

      ENDDO
      ENDDO

!     Diagnose der 10m-Groessen:

      CALL diag_level(i_st, i_en, j_st, j_en, z10m_2d, k_2d, hk_2d, hk1_2d)

      DO j=j_st,j_en
      DO i=i_st,i_en

         ii=MAX( i-vst, 1 )
         jj=MAX( j-vst, 1 )

         IF (k_2d(i,j).EQ.ke) THEN

!           10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!           in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

            z0d=z0d_2d(i,j)
            a_atm=h_atm_2d(i,j)+z0d
            a_10m=h_10m+z0d

            fac_m=(rat_m_2d(i,j)-z1)*z0d/h_top_2d(i,j)
            IF (fac_m.EQ.z1) THEN
               val1=h_10m/a_10m
               val2=h_atm_2d(i,j)/a_atm
            ELSE
               val1=LOG(a_10m/(z0d+fac_m*h_10m))
               val2=LOG(a_atm/(z0d+fac_m*h_atm_2d(i,j)))
            END IF

            fakt=val1/val2

            u_10m(i,j)=vel1_2d(i,j)*fakt
            v_10m(i,j)=vel2_2d(i,j)*fakt
              
         ELSE

            k = k_2d(i,j)

!           10m-Niveau oberhalb der untersten Modell-Hauptflaeche:

!US this has been forgotten before
            uk =(u(i,j,k)+u(ii,j,k))*z1d2
            vk =(v(i,j,k)+v(i,jj,k))*z1d2
!US

            uk1=(u(i,j,k+1)+u(ii,j,k+1))*z1d2
            vk1=(v(i,j,k+1)+v(i,jj,k+1))*z1d2

            val1=(h_10m+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))
            val2=(hk_2d(i,j)+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))

            fakt=LOG(val1)/LOG(val2)

            u_10m(i,j)=uk1+fakt*(uk-uk1)
            v_10m(i,j)=vk1+fakt*(vk-vk1)
         END IF

!-----------------------------------------------------------------------
!        nur im alten 1-d-Modell:

!        Berechnung der der Flussdichten fuer Impuls,
!        fuehlbare und latente Energie und der Ri-Zahl:

!        flimpuls(i,j,ke1)=rho_2d(i,j)*tkvm(i,j,ke1)*SQRT(fm2_2d(i,j))
!        flfbwae(i,j,ke1)=cp_d*rho_2d(i,j)*tkvh(i,j,ke1)*grad(i,j,tet_l)
!        fllawae(i,j,ke1)=lh_v*rho_2d(i,j)*tkvh(i,j,ke1)*grad(i,j,h2o_g)
 
!        IF (fh2_2d(i,j).EQ.z0) THEN
!           rizahl(i,j,ke1)=z0
!        ELSE
!           rizahl(i,j,ke1)=fh2_2d(i,j)/MAX( fm2_2d(i,j), epsi*ABS(fh2_2d(i,j)) )
!        END IF

!        Beachte: flimpuls,flfbwae,fllawae und die Ri-Zahl
!                 werden nur fuer Diagnosezwecke berechnet.
!-----------------------------------------------------------------------
      END DO
      END DO

      IF (PRESENT(shfl_s)) THEN
         DO j=j_st,j_en
         DO i=i_st,i_en
            shfl_s(i,j)=-cp_d*rho_2d(i,j)*tkvh(i,j,ke1)*grad(i,j,tet_l)
         END DO
         END DO
      END IF
      IF (PRESENT(lhfl_s)) THEN
         DO j=j_st,j_en
         DO i=i_st,i_en
            lhfl_s(i,j)=-lh_v*rho_2d(i,j)*tkvh(i,j,ke1)*grad(i,j,h2o_g)
         END DO
         END DO
      END IF
         
!-----------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN

         TKE_SCLM%mod(ke1)%val=tke(im,jm,ke1,ntur); TKE_SCLM%mod(ke1)%vst=i_cal

         !Achtung: TKE%mod(ke1)%val zeigt z.Z. noch auf die alte TKE-Zeitstufe.
         !Somit wird also der alte tke-Wert mit dem neuen ueberschrieben,
         !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
         !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

         velh=tkvm(im,jm,ke1)*SQRT(fm2_2d(im,jm))
         wert=tkvh(im,jm,ke1)*fh2_2d(im,jm)

         teta=-tkvh(im,jm,ke1)*grad(im,jm,tet_l) !teta cov. with w
         qvap=-tkvh(im,jm,ke1)*grad(im,jm,h2o_g) !qvap cov. with w

         UST%mod(ke1)%val=SQRT(velh) ; UST%mod(ke1)%vst=i_cal
         TWA%mod(ke1)%val=teta       ; TWA%mod(ke1)%vst=i_cal

         SHF%mod(0)%val=cp_d*rho_2d(im,jm)*teta ; SHF%mod(0)%vst=i_cal
         LHF%mod(0)%val=lh_v*rho_2d(im,jm)*qvap ; LHF%mod(0)%vst=i_cal
       ! LMO%mod(0)%val=velh**3/wert ; LMO%mod(0)%vst=i_cal
      END IF
#endif
!SCLM-------------------------------------------------------------------

END SUBROUTINE turbtran

!********************************************************************************

!+ Module procedure diag_level in "src_turbdiff" for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, j_st, j_en, zdia_2d, k_2d, hk_2d, hk1_2d)

   IMPLICIT NONE

   INTEGER (KIND=iintegers), INTENT(IN) :: &
!
      i_st, i_en, j_st, j_en !start end end indices of horizontal domain

   REAL (KIND=ireals), INTENT(IN) :: &
!
      zdia_2d(:,:)  !diagnostic height

   INTEGER (KIND=iintegers), INTENT(INOUT) :: &
!
      k_2d(:,:)     !index field of the upper level index
                    !to be used of near surface diagnostics

   REAL (KIND=ireals), INTENT(INOUT) :: &
!
      hk_2d(:,:), & !mid level height above ground belonging to 'k_2d'
     hk1_2d(:,:)    !mid level height above ground of the previous layer (below)

   INTEGER (KIND=iintegers) :: i, j

   LOGICAL :: lcheck

   lcheck=.TRUE. !check whether a diagnostic level is above the current layer

   DO WHILE (lcheck) !loop while previous layer had to be checked
      lcheck=.FALSE. !check next layer ony, if diagnostic level is at least once
                     !above the current layer
      DO j=j_st,j_en
      DO i=i_st,i_en
         IF (hk_2d(i,j)<zdia_2d(i,j) .AND. k_2d(i,j)>1) THEN !diagnostic level is above current layer
            lcheck=lcheck .OR. .TRUE. !for this point or any previous one the
                                      !diagnostic level is above the current layer
            k_2d(i,j)=k_2d(i,j)-1
            hk1_2d(i,j)=hk_2d(i,j)
            hk_2d(i,j)=(hhl(i,j,k_2d(i,j))+hhl(i,j,k_2d(i,j)+1))*z1d2-hhl(i,j,ke1)
          END IF
       END DO   
       END DO   

   END DO   

END SUBROUTINE diag_level

!********************************************************************************

!+ Module procedure turbdiff in "src_turbdiff" for computing the tendencies
!+ for vertical diffusion


SUBROUTINE turbdiff(dt_var,dt_tke,lstfnct)

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

! Array arguments with intent(in):

      REAL (KIND=ireals), INTENT (IN) :: &
             dt_var, &  !Zeitschrittlaenge fuer die turbulente Diffusion
             dt_tke     !Zeitschrittlaenge fuer die TKE-Prognose

      LOGICAL, INTENT (IN) :: lstfnct

!-----------------------------------------------------------------------
!
! Local parameters:

! Local scalars:

!     Lokale logical Variablen:

      LOGICAL limplizit, & !(teil-)implizite Berechnung der Diffus.tend.
              lexplizit, & !(teil-)explizite Berechnung der Diffus.tend.
              lcircterm(ie,je)    !Zirkulationsterm muss berechnet werden

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers), PARAMETER :: &
!
!             Indexgrenzen:
!
              nscal=3,     & !Skalare Groessen
              ninv=2,      & !daraus abgeleitete gegenueber vertikalen
                             !(feuchtadiabatischen) Verrueckungen
                             !invarianten Groessen
              nvel=2,      & !Geschwindigkeitskomponenten
              ndiff=nscal+nvel, &
              nred=ninv+nvel,   &
!              ntmax=3,     & !max. Anzahl der Zeitebenen fuer die TKE
!
!             Zeiger fuer die Variablen :
!
              u_m=1,       & !zonale Geschw.komp. im Massenzentrum
              v_m=2,       & !meridionale  ,,      ,,     ,,
              tet_l=3,     & !feucht-potentielle Temperatur
              tem_l=tet_l, & !Fluessigwasser-Temperatur
              h2o_g=4,     & !Gesamtwasseergehalt
              liq=5,       & !Fluessigwasser  ,,
!
              u_s=u_m,     & !zonale Geschw.komp. im gestag. Gitter
              v_s=v_m,     & !meridionale  ,,     ,,   ,,       ,,
              tet=tet_l,   & !pot.Temperatur
              tem=tet,     & !Temperatur
              vap=h2o_g      !Wasserdampfmischungsverh.

!     Beachte:u_m,v_m, sowie u_s,v_s muessen in [1,nvel] liegen;
!             aber: tet_l,h2o_g in [nvel+1,nred]
!             und   tem (tet),vap,liq in [nvel+1,ndiff]

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &
!
              i,j,k,        & !Diskretis.index fuer lamda, phi, sigma
              n,ii,jj,kk,   & !Indices fuer diverse Schleifen
              ku,ko,k1,k2,  & !unterer und oberer Schicht-Index
              nkorr,nk1,nk2,& !Startindices der explizit zu behandelnden
                              !Variablen
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
!     Platzh. fuer pot. Temp., relat. Wolkenanteil, rezipr. virtueller Faktor
!     und Temperaturtendenz des Wassersaettigungsmisch.verh.:
!
           teta,rcl,virt,qsat_dT

      REAL (KIND=ireals) :: &
!
           exnr !Exner-Faktor

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
           h,hu,l_turb,l_diss,lh,lm,kohae_len,l_lam,len, &
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
           tmp1,tmp2,tmp3,aa,bb, &
           beta,x1,x2,x3

      REAL (KIND=ireals) :: flukon33, flukon43, flukon53, &
                            flukon34, flukon44, flukon54, flukon55, &
                            flux_3,   flux_4,   flux_5

!---------------------------------------------------------------------
#ifdef __COSMO__
      REAL (KIND=ireals) :: zs11, zs22, zs33, zs12, zs21, zs13, zs23
#endif
!---------------------------------------------------------------------

! Local arrays:

      LOGICAL            :: lo_ice(ie,je) ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=ireals), TARGET :: &
!
           vari(ie,je,ke1,ndiff) ,& !(co-)variances of the velocity-components
                                    !reduced set of variables in the first part of 'turbdiff'
!
           dicke(ie,je,ke1),& !(effektive) Dicke einer Modellschicht 
                              !bzw. sonstiges Hilfsfeld, bei k=ke1 steht
                              !eine effekt. Dicke der Modell-Prandtlsch.
!
           wind(ie,je,ke1,nvel),& !horizontale Windkomponenten im
                                  !Massenzentrum, bzw. deren Tendenzen
!
           len_scale(ie,je,ke1),& !turbulent length-scale
!
           rhon(ie,je,ke1),& !Luftdichte auf Nebenflaechen
                             !einschl. des skin-layers
!     
!     Hilfsfeld fuer die korrigierten effektiven Gradienten der
!     ndiff diffundierenden Modellvariablen:
!
           a(ie,je,ke1,ndiff), &
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
           hlp(ie,je,ke1), &
           zaux_rcld(ie,je,ke1),&
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
           lay(ie,je), lays(ie,je,2), &
!
           pr(ie,je),   & !Druck             (am atm. Modellunterrand)
           tp(ie,je),   & !Temperatur              ,,       
           qd(ie,je),   & !Wasserdampfgehalt       ,,
           ql(ie,je),   & !Fluessigwassergehalt    ,,
           z0m(ie,je),  & !Rauhigkeitslaenge (fuer Impuls)
           l_pat(ie,je),& !Laengenskala der therm. Inhomogenitaeten
                          !der Erdbodenoberflaeche
!
           frc(ie,je),  &! forcing
           src(ie,je),  &! source term
!
!     Vertikale Gradienten und turbulente Flussdichten verschiedener 
!     thermodynamischer Variablen:
!
           grad(ndiff), &
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
           utens(:,:,:), vtens(:,:,:), &
           ttens(:,:,:), qvtens(:,:,:), qctens(:,:,:), &
!
!          Pointer fuer Felder des Exnerfaktors und Bestandes-Parameter:
!
           exner(:,:,:), &
!
           dpat(:,:),    &
           cbig(:,:,:), csml(:,:,:), rair(:,:,:)

      REAL (KIND=ireals), TARGET :: &
!
!          Seicher-Target fuer obige Pointer:
!
           exner_tar(ie,je,ke), &
!
           dpat_tar(ie,je),   &
           cbig_tar(ie,je,kcm:ke1), csml_tar(ie,je,kcm:ke1), rair_tar(ie,je,kcm:ke1)
 
!          Note:
!          The following buffers wouldn't be necessary, if the related pointers above
!          were allowed to be allocated at run time:

      LOGICAL ::  &
!          
           can_fields, & !canopy fields are present
           ltend(ndiff)  !calculation of tendencies required

!-------------------------------------------------------------------------------
#ifdef SCLM
      REAL (KIND=iintegers) :: cvar(ndiff,2) !variance and covariance with vertical wind
#endif
!-------------------------------------------------------------------------------

!---- End of header ------------------------------------------------------------

!     nur in LM-Umgebung:

      istat=0; ilocstat=0; ierrstat=0
      errormsg = ''; eroutine='turbtran'; lerror=.FALSE.

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

!print *,"in turbdiff iini=",iini

      IF (PRESENT(epr)) THEN
         exner => epr
      ELSE
         exner => exner_tar
       ! ALLOCATE ( exner(ie,je,ke),     STAT=ilocstat ); istat = istat + ilocstat
         DO k=1, ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               exner(i,j,k)=zexner(prs(i,j,k))
            END DO
            END DO
         END DO
      END IF       

      can_fields=.TRUE.
      IF (PRESENT(d_pat)) THEN
         dpat => d_pat
      ELSE
         dpat => dpat_tar
       ! ALLOCATE ( dpat(ie,je),         STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF
      IF (PRESENT(c_big)) THEN
         cbig => c_big
      ELSE
         cbig => cbig_tar 
       ! ALLOCATE ( cbig(ie,je,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF
      IF (PRESENT(c_sml)) THEN
         csml => c_sml
      ELSE
         csml => csml_tar 
       ! ALLOCATE ( csml(ie,je,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
         can_fields=.FALSE.
      END IF
      IF (PRESENT(r_air)) THEN
         rair => r_air
      ELSE
         rair => rair_tar 
       ! ALLOCATE ( rair(ie,je,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
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
         CALL init_canopy(ie=ie, je=je, ke=ke, ke1=ke1, kcm=kcm, &
#ifdef __ICON__
              itype_tran=itype_tran,                             &
#endif
              istartpar=istartpar, iendpar=iendpar, jstartpar=jstartpar, jendpar=jendpar, &
              fr_land=fr_land, &
              d_pat=dpat, c_big=cbig, c_sml=csml, r_air=rair) 
!print *,"nach init_canopy kcm=",kcm
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
            tinv(n)=z1/dt_var !division by time increment for variable increments
         ELSE               !update of prognostic variables required
            tinc(n)=dt_var    !time increment multiplication for tendencies
            tinv(n)=z1        !no division by time increment for variable increments
         END IF
!print *,"n=",n," tinc=",tinc(n)," tinv=",tinv(n)
      END DO

      !Note:
      !If a tendency field of an ordinary prognostic variable is not present,
      !the related time step increment due to turbulent diffusion will be
      !added to the prognostic variable directly.
             
!-----------------------------------------------------------------------

! 1)  Vorbereitungen:

!     Setzen von Steuerparametern:

      IF (imode_turb.GE.3) THEN
         limplizit=.TRUE.
      ELSE
         limplizit=.FALSE.
      END IF

      IF (imode_turb.EQ.2) THEN
         lexplizit=.TRUE.
         nkorr=1
      ELSEIF (.NOT.lexpcor .OR. imode_turb.EQ.4) THEN
         lexplizit=.FALSE.
         nkorr=ndiff+1
      ELSE
         lexplizit=.TRUE.
         IF (lnonloc) THEN
            nkorr=1
         ELSE
            nkorr=nvel+1
         END IF
      END IF

      nk1=nkorr
      nk2=MAX( nkorr, nvel+1 )

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

      CALL turb_param

!     Bestimmung der initialen Werte fuer die laminaren Reduktions-
!     koeffizienten :

      IF (lini.AND.itype_tran.NE.2) THEN 
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            tfh(i,j)=z1
            tfm(i,j)=z1
         END DO
         END DO
      END IF

!     Berechnung der effektiven Bedeckung mit nichtkonvektiven
!     Wasserwolken:

      IF (icldm_turb.EQ.-1) THEN       
!        Keine Wolkenberuecksichtigung;
!        Wolkenwasser wird ignoriert (vollkommen trockenes Schema):

         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               rcld(i,j,k)=z0
               hlp(i,j,k)=z0
            END DO
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
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               rcld(i,j,k)=z0
               hlp(i,j,k)=-qc(i,j,k)
            END DO
            END DO
         END DO

      ELSEIF (icldm_turb.EQ.1) THEN
!       Wolkenwasser wird als skalige Wolke interpretiert:

        DO k=1,ke
          DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
              IF ( qc(i,j,k) > z0 ) THEN
                rcld(i,j,k) = z1
              ELSE
                rcld(i,j,k) = z0
              END IF
              hlp(i,j,k)=z0
            ENDDO
          ENDDO
        ENDDO

      ELSEIF (icldm_turb.EQ.2) THEN
!        spezielle Diagnose von Wasserwolken:

         zaux_rcld(:,:,:)=0._ireals

         CALL turb_cloud (                                &
#ifdef __ICON__
              itype_wcld,                                 & 
#endif
              ie, je, ke,ke1,           1,ke,             &
              istartpar,iendpar, jstartpar,jendpar, 1,ke, &
              prs, ps, rcld, t, qv, qc,                   &
              zaux_rcld, hlp                                     )

         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               hlp(i,j,k) =hlp(i,j,k)-qc(i,j,k)
               rcld(i,j,k)=zaux_rcld(i,j,k)
            END DO
            END DO
         END DO

      END IF
      END IF
      
      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar    
         hlp(i,j,ke1)=hlp(i,j,ke)
      END DO
      END DO

      IF (icldm_tran.EQ.0) THEN
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            rcld(i,j,ke1)=z0
         END DO
         END DO
      ELSEIF (icldm_tran.EQ.1) THEN
         DO j=jstartpar,jendpar
           DO i=istartpar,iendpar    
             IF ( qc(i,j,ke) > z0 ) THEN
               rcld(i,j,ke1) = z1
             ELSE
               rcld(i,j,ke1) = z0
             END IF
           END DO
         END DO
      ELSE
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            rcld(i,j,ke1)=rcld(i,j,ke)
         END DO
         END DO
      END IF

!print *,"nach rcld"

!     Berechnung der thermodynamischen Hilfsfelder:

!     Bem.:  Am unteren Modellrand gibt es eine skalige Erdboden-
!            oberflaeche mit einer vertikal nicht aufgeloesten 
!            Rauhigkeit der Laenge d+z0m. Diese von Rauhigkeits-
!            elementen durchsetzte 'Mikro'-Rauhigkeitsschicht zaehlt
!            nicht mehr zum atm. Modellgebiet! Sie stellt einen 
!            Mikrobestand der Hoehe d+z0m dar, in dem die konstante
!            Schliessungslaenge akt*z0m gelten soll.

!            Es wird davon ausgegangen, dass das Niveau der Hoehe
!            hhl(i,j,ke1) (Erdbodenniveau) gerade der Untergrenze der
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
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               tp(i,j)=t(i,j,k)+zlhocp*hlp(i,j,k)
               qd(i,j)=qv(i,j,k)-hlp(i,j,k)
               ql(i,j)=qc(i,j,k)+hlp(i,j,k)
               pr(i,j)=prs(i,j,k)
               a(i,j,k,2)=exner(i,j,k) !Exner-Faktor
                                       !(wird spaeter auf Nebenflaechen interpoliert)
            END DO   
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
 
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               tp(i,j)=tp(i,j)*(z1-tfh(i,j)) &
!                     +(t_g(i,j)+zlhocp*hlp(i,j,ke1))*tfh(i,j)
                      +t_g(i,j)*tfh(i,j)
               qd(i,j)=qd(i,j)*(z1-tfh(i,j)) &
!                     +(qv_s(i,j)-hlp(i,j,ke1))*tfh(i,j)
                      +qv_s(i,j)*tfh(i,j)
               pr(i,j)=ps(i,j) 
               a(i,j,k,2)=zexner(pr(i,j)) !Exner-Faktor am Boden
            END DO
            END DO
         END IF

!        Thermdynamische Hilfsvariablen:

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            exnr   =a(i,j,k,2)    !Exner-Faktor
            teta   =tp(i,j)/exnr  !potentielle Temperatur

            rcl    =c_scld/(z1+rcld(i,j,k)*(c_scld-z1))*rcld(i,j,k) !Bedeckungsgrad

            virt   =z1/(z1+rvd_m_o*qd(i,j)-ql(i,j))      !rezipr. virtueller Faktor
            qsat_dT=zdqsdt( tp(i,j), zqvap( zpsat_w( tp(i,j) ), pr(i,j) ) ) !dQs/dT
            fakt   =(zlhocp/tp(i,j)-(z1+rvd_m_o)*virt)/(z1+qsat_dT*zlhocp)

            lay(i,j)=virt

            a(i,j,k,1)=z1+zrcpv*qd(i,j)+zrcpl*ql(i,j) !Cp/Cpd
            a(i,j,k,3)=qsat_dT

!           Berechn. der thermodyn. Faktoren A_tet und A_q:
            a(i,j,k,4)=grav*(z1/teta-rcl*fakt*exnr*qsat_dT)
            a(i,j,k,5)=grav*(rvd_m_o*virt+rcl*fakt)

!           Berechnung der feucht-potentiellen Temp. und des
!           Gesamtwassergehaltes:
            vari(i,j,k,tet_l)=teta-ql(i,j)*zlhocp/exnr
            vari(i,j,k,h2o_g)=qd(i,j)+liqfak*ql(i,j)
            vari(i,j,k,liq  )=ql(i,j)
         END DO
         END DO

         IF (PRESENT(rho) .AND. k.LT.ke1) THEN
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               rhon(i,j,k)=rho(i,j,k) !Luftdichte
            END DO
            END DO
         ELSE    
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               rhon(i,j,k)=lay(i,j)*pr(i,j)/(r_d*tp(i,j)) !Luftdichte
            END DO
            END DO
         END IF   

!        Die thermodynamischen Hilfsgroessen wurden hier
!        unter Beruecksichtigung der diagnostizierten
!        Kondensationskorrektur gebildet, indem die entspr.
!        korrigierten Werte fuer tem, qv und ql benutzt wurden.
!        Beachte, dass es zu einer gewissen Inkosistenz kommt,
!        wenn die Dichte und der Exnerfaktor bereits vorhanden sind!

      END DO
!print *,"nach Erhaltungsvar"

!     Beachte:
!     tp(), qd(), ql() und pr() gehoeren jetzt zu k=ke1 

!     Berechnung der horizontalen Windgeschwindigkeiten
!     im Massenzentrum der Gitterbox:

      DO j=jstartpar,jendpar
         jj=MAX( j-vst, 1 )
         DO k=1,ke
            DO i=istartpar,iendpar
               ii=MAX( i-vst, 1 )
               vari(i,j,k,u_m)=(u(i,j,k)+u(ii,j,k))*z1d2 
               vari(i,j,k,v_m)=(v(i,j,k)+v(i,jj,k))*z1d2
            END DO
         END DO
         DO i=istartpar,iendpar
            vari(i,j,ke1,u_m)=vari(i,j,ke,u_m)*(z1-tfm(i,j))
            vari(i,j,ke1,v_m)=vari(i,j,ke,v_m)*(z1-tfm(i,j))
         END DO
      END DO         

!     Sichern der auf Massenpunkte interpolierten Windkomponenten:

      DO n=1,nvel
         DO k=1,ke1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               wind(i,j,k,n)=vari(i,j,k,n)
            END DO
            END DO
         END DO
      END DO    

!     Berechnung der Modellschichtdicken:
        
      DO k=1,ke
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            dicke(i,j,k)=hhl(i,j,k)-hhl(i,j,k+1)
         END DO
         END DO
      END DO   

!     Interpolation der thermodyn. Hilfsgroessen im Feld a(),
!     der Wolkendichte und der Luftdichte auf Nebenflaechen:

      DO k=ke,2,-1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            hlp(i,j,k)=dp0(i,j,k)+dp0(i,j,k-1)
         END DO
         END DO
      END DO
      DO ii=1,ndiff
         DO k=ke,2,-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,k,ii)=(a(i,j,k,ii)*dp0(i,j,k-1)+a(i,j,k-1,ii)*dp0(i,j,k))/hlp(i,j,k)
            END DO
            END DO
         END DO
      END DO
      DO k=ke,2,-1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            rcld(i,j,k)=(rcld(i,j,k)*dp0(i,j,k-1)+rcld(i,j,k-1)*dp0(i,j,k))/hlp(i,j,k)
            rhon(i,j,k)=(rhon(i,j,k)*dicke(i,j,k-1)  &
                        +rhon(i,j,k-1)*dicke(i,j,k)) &
                       /(hhl(i,j,k-1)-hhl(i,j,k+1))

!           Beachte die andere Behandlung der Dichte!
         END DO
         END DO
      END DO
!print *,"nach interpol"

!     Bestimmun der initialen Werte fuer die Rauhigkeitslaenge
!     ueber Wasserpunkten:

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

      DO j=jstartpar,jendpar
        DO i=istartpar,iendpar

          IF (fr_land(i,j).LT.z1d2) THEN
            ! Water point.
            IF (.NOT. lseaice) THEN
              ! Sea ice model is not used.
              ! Ice surface if SST is less than the salt water freezing temperature.
              lo_ice(i,j) = t_g(i,j) < t0_melt + zt_ice
            ELSE
              ! Sea ice model is used.
              ! Ice surface if ice is present.
              lo_ice(i,j) = h_ice(i,j) > z0
            END IF
            IF (llake) THEN
              ! Lake model is used.
              ! Ice surface if this is a lake point AND ice is present.
              IF ((depth_lk(i,j) > z0) .AND. (h_ice(i,j) >= h_Ice_min_flk)) &
              lo_ice(i,j) = .TRUE.
            END IF
          END IF

        END DO
      END DO
!print *,"nach li_ice"

      IF (lini.AND.itype_tran.EQ.3) THEN 

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            IF (fr_land(i,j).LT.z1d2) THEN

               ! Use ice surface roughness or open-water surface roughness
               ! according to lo_ice
               IF ( lo_ice(i,j) ) THEN
                  gz0(i,j)=grav*z0_ice
               ELSE
!                 Bei von Schubspannung abhaengiger Wellenhoehe:

!                 Einfachste Schaetzung der Schubspannung als Impusls-
!                 flussdichte durch die Nebenflaeche ke mit Hilfe
!                 einer diagnostischen TKE und Bestimmung der Rauhig-
!                 keitslaenge unter Anwendung der Charnockformel:

                  l_turb=dicke(i,j,ke)
                  l_turb=akt*MAX( len_min, l_turb/(z1+l_turb/l_scal) )
                  edh=z2/(hhl(i,j,ke-1)-hhl(i,j,ke1))

               !  DO n=1,nred
               !     grad(n)=(vari(i,j,ke-1,n)-vari(i,j,ke,n))*edh
               !  END DO
               !  nred=4 in all cases,
               !  thus these expressions are vectorized correctly:
                  grad(1)=(vari(i,j,k-1,1)-vari(i,j,k,1))*edh
                  grad(2)=(vari(i,j,k-1,2)-vari(i,j,k,2))*edh
                  grad(3)=(vari(i,j,k-1,3)-vari(i,j,k,3))*edh
                  grad(4)=(vari(i,j,k-1,4)-vari(i,j,k,4))*edh

                  fh2=a(i,j,k,4)*grad(tet_l)+a(i,j,k,5)*grad(h2o_g)
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
                  gz0(i,j)=MAX( grav*len_min, alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF   
         END DO   
         END DO   
      END IF     
!print *,"nach gz0-ini"

!     Bestimmung der Rauhigkeitslaenge 
!     und der effektiven Dicke der Modell-Prandtlschicht:

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar    
         z0m(i,j)=gz0(i,j)/grav
         dicke(i,j,ke1)=z0m(i,j) &
                       *log(z1d2*dicke(i,j,ke)/z0m(i,j)+z1)
      END DO
      END DO

!     Berechnung der tubulenten Laengenscalen:

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         len_scale(i,j,ke1)=z0m(i,j)
      END DO
      END DO
      DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            IF (cbig(i,j,k).gt.z0) THEN

!              Die turbulente Laengenskala wird durch die Laengen-
!              skala der lufterfuellten Zwischenraeume limitiert:

               l_turb=z1/(cbig(i,j,k)*sqrt(z1/exp(rair(i,j,k))-z1))
               len_scale(i,j,k)=MIN( dicke(i,j,k)+len_scale(i,j,k+1), l_turb )
            ELSE
               len_scale(i,j,k)=dicke(i,j,k)+len_scale(i,j,k+1)
            END IF
         END DO
         END DO
      END DO   
      DO k=kcm-1,1,-1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            len_scale(i,j,k)=dicke(i,j,k)+len_scale(i,j,k+1)
         END DO
         END DO
      END DO   

!     Uebergang von der maximalen turbulenten Laengenskala zur
!     effektiven turbulenten Laengenskala:
!print *,"l_scal=",l_scal

      DO k=1,ke1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
               len_scale(i,j,k)=akt*MAX( len_min, &
                                         len_scale(i,j,k)/(z1+len_scale(i,j,k)/l_scal) )
         END DO
         END DO
      END DO       
!print *,"nach len_scale"

!     Initialisierung der Felder fuer tke,tkvh,tkvm:

      IF (lini) THEN  !nur beim allerersten Durchgang

         DO k=2,kem
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              der Einfachheit halber nur lokale Berechnung der
!              vertikalen Gradienten:

               len=len_scale(i,j,k)

               IF (k.EQ.ke1) THEN
                  edh=z1/dicke(i,j,ke1)
               ELSE
                  edh=z2/(hhl(i,j,k+1)-hhl(i,j,k-1))
               END IF

!NEC_CB        DO n=1,nred
!NEC_CB           grad(n)=(vari(i,j,k,n)-vari(i,j,k-1,n))*edh
!NEC_CB        END DO
!              ! nred=4 in all cases,
               ! thus these expressions are vectorized correctly:
               grad(1)=(vari(i,j,k,1)-vari(i,j,k-1,1))*edh
               grad(2)=(vari(i,j,k,2)-vari(i,j,k-1,2))*edh
               grad(3)=(vari(i,j,k,3)-vari(i,j,k-1,3))*edh
               grad(4)=(vari(i,j,k,4)-vari(i,j,k-1,4))*edh

               fh2=a(i,j,k,4)*grad(tet_l)+a(i,j,k,5)*grad(h2o_g)
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

               lay(i,j)=MAX( vel_min, SQRT(d_m*len*wert) ) !Initialwert fuer SQRT(2TKE)

               tkvm(i,j,k)=lm
               tkvh(i,j,k)=lh

!              Am Anfang konnte noch keine Advektion oder 
!              Diffusion von SQRT(2*TKE) berechnet werden:

               tketens(i,j,k)=z0
            END DO 
            END DO        

            DO n=1,ntim
              DO j=jstartpar,jendpar
                DO i=istartpar,iendpar
                  tke(i,j,k,n)=lay(i,j)
                END DO
              END DO
            END DO

         END DO    

         DO n=1,ntim
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tke(i,j,1,n)=tke(i,j,2,n)
            END DO
            END DO
         END DO

      ELSE

         DO k=2,kem
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tkvh(i,j,k)=tkvh(i,j,k)/tke(i,j,k,nvor)
               tkvm(i,j,k)=tkvm(i,j,k)/tke(i,j,k,nvor)
            END DO
            END DO
         END DO    

      END IF
!print *,"nach tke-init"

!     tkvh und tkvm enthalten jetzt die stabilitaetsabhaengigen 
!     Laengenmasse, nicht die Diffusionskoeffizienten!


! 2)  Berechnung der benoetigten vertikalen Gradienten und
!     Abspeichern auf vari():

!     Am unteren Modellrand:

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         hlp(i,j,ke1)=z1/dicke(i,j,ke1)
         vari(i,j,ke1,liq)=z0
      END DO
      END DO
 
      DO n=1,nred 
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            vari(i,j,ke1,n)=(vari(i,j,ke,n)-vari(i,j,ke1,n))*hlp(i,j,ke1)
         END DO
         END DO
      END DO   

!     An den darueberliegenden Nebenflaechen:

      IF (lnonloc) THEN

!        Berechnung nicht-lokaler Gradienten:

         DO n=1,ndiff

!           Berechnung vertikalen Integralfunktionen in hlp():
      
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,ke1)=z0
            END DO 
            END DO 

            DO k=ke,2,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,k)=hlp(i,j,k+1)+vari(i,j,k,n) &
                                         *(hhl(i,j,k)-hhl(i,j,k+1))
               END DO
               END DO
            END DO

            k1=1
            k2=2
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,k1)=hhl(i,j,1)-hhl(i,j,2)
            END DO
            END DO
            DO k=2,ke   

!              Berechnung der nicht-lokalen Gradienten als mittlere
!              Differenzenquotienten ueber die stabilitaetsabhaengige
!              Laengenskala in tkvh() bzw tkvm():

!              Speichern der stab.abh. Laengenskala unter lay():

               IF (n.LE.nvel) THEN 
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lay(i,j)=tkvm(i,j,k)
                  END DO
                  END DO
               ELSE   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lay(i,j)=tkvh(i,j,k)
                  END DO
                  END DO
               END IF
                  
!              Bestimmung der nicht-lokalen Gradienten und 
!              Zwischenspeichern derselben auf dem Feld dicke():

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar


                  lays(i,j,k2)=hhl(i,j,k)-hhl(i,j,k+1)

                  IF (lay(i,j).LE. &
                      z1d2*MIN( lays(i,j,k1), lays(i,j,k2)) ) THEN

!                    Die vertikalen Diffusionswege schneiden weder
!                    eine untere noch eine obere Hauptflaeche. Es 
!                    koennen somit die lokalen Gradienten genommen 
!                    werden. Bei sehr kleinen Diffusionslaengen, muss 
!                    aus num. Gruenden sogar die lokale Berechnung 
!                    gewaehlt werden:

                     dicke(i,j,k)=z2*(vari(i,j,k-1,n)-vari(i,j,k,n)) &
                                    /(hhl(i,j,k-1)-hhl(i,j,k+1))
                  ELSE

!                    Berechn. der benoetigten Referenzhoehen und -level:

                     h=hhl(i,j,k)
                     hu=MAX( h-lay(i,j), hhl(i,j,ke1) )
                     hig(1)=hu+lay(i,j)
                     hig(2)=h+lay(i,j)
 
                     kk=k
111                  IF (hhl(i,j,kk).GT.hu) THEN
                        kk=kk+1
                        GOTO 111
                     END IF
                     ku=kk
                     DO ii=1,2
112                     IF (kk.GT.1) THEN
                           IF (hhl(i,j,kk).LE.hig(ii)) THEN
                              kk=kk-1
                              GOTO 112
                           END IF
                        END IF
                        lev(ii)=kk+1
                     END DO

!                    Berechnung der gemittelten Differenzenquotienten
!                    als Ausdruck fuer die nicht-lokalen Gradienten: 

                     wert=hlp(i,j,ku)-hlp(i,j,k) &
                         +hlp(i,j,lev(2))-hlp(i,j,lev(1)) &
                         +vari(i,j,ku-1,n)*(hu-hhl(i,j,ku)) &
                         -vari(i,j,lev(1)-1,n)*(hig(1)-hhl(i,j,lev(1)))&
                         +vari(i,j,lev(2)-1,n)*(hig(2)-hhl(i,j,lev(2)))

                     dicke(i,j,k)=wert/(lay(i,j)*(h-hu))
                  END IF   
               END DO    
               END DO    
               kk=k1
               k1=k2
               k2=kk

            END DO   

!           Sichern der nicht-lokalen Gradienten im Feld vari():

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vari(i,j,k,n)=dicke(i,j,k)
               END DO
               END DO
            END DO

         END DO

!        Belegung von dicke() mit den Schichtdicken
!        bzgl. Nebenflaechen:

         DO k=2,ke   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=(hhl(i,j,k-1)-hhl(i,j,k+1))*z1d2
            END DO
            END DO
         END DO      

      ELSE       

!        Berechnung lokaler Gradienten:
         
         DO k=ke,2,-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=(hhl(i,j,k-1)-hhl(i,j,k+1))*z1d2
               hlp(i,j,k)=z1/dicke(i,j,k)
            END DO
            END DO
         END DO   

         DO n=1,ndiff
            DO k=ke,2,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vari(i,j,k,n)=(vari(i,j,k-1,n)-vari(i,j,k,n))*hlp(i,j,k)
               END DO
               END DO
            END DO   
         END DO      

      END IF
!print *,"nach gradient"

!------------------------------------------------------------------------------------
#ifdef __COSMO__
      IF (itype_sher.GT.1 .OR. PRESENT(tket_hshr)) THEN
         !3D-turbulent shear or separate horizontal shear mode to be considered as
         !TKE source or at least to be calculated for output:

!        Berechnung der 3D-Korrektur des Scherungsterms der TKE incl. metr. Terme:

         DO k=2,ke

            !Neigung der Modellflaeche:
            DO j = jstart, jend
            DO i = istart, iend
               lay(i,j) = ( hhl(i+1,j,k) - hhl(i-1,j,k) ) * z1d2
               frc(i,j) = ( hhl(i,j+1,k) - hhl(i,j-1,k) ) * z1d2
            END DO
            END DO

            !Scherung durch horizontale Konfluenz:
            DO j = jstart, jend
            DO i = istart, iend
               zs11 = ( u(i,j,k  ) - u(i-1,j,k  )   +       &
                        u(i,j,k-1) - u(i-1,j,k-1) ) * z1d2
               zs22 = ( v(i,j,k  ) - v(i,j-1,k  )   +       &
                        v(i,j,k-1) - v(i,j-1,k-1) ) * z1d2

               zs11 = ( zs11 - lay(i,j)*vari(i,j,k,u_m) ) * eddlon * acrlat(j,1)
               zs22 = ( zs22 - frc(i,j)*vari(i,j,k,v_m) ) * edadlat

               hlp(i,j,k) = z2 * ( zs11**2 + zs22**2 )
            END DO
            END DO

            !Horizontale Scherung:
            DO j = jstart, jend
            DO i = istart-1, iend
               lays(i,j,1) = ( v(i  ,j,k  ) + v(i  ,j-1,k  ) +     &
                               v(i+1,j,k  ) + v(i+1,j-1,k  ) +     &
                               v(i  ,j,k-1) + v(i  ,j-1,k-1) +     &
                               v(i+1,j,k-1) + v(i+1,j-1,k-1) )     &
                               * 0.125_ireals
            END DO
            END DO
            DO j = jstart-1, jend
            DO i = istart, iend
               lays(i,j,2) = ( u(i,j  ,k  ) + u(i-1,j  ,k  ) +     &
                               u(i,j+1,k  ) + u(i-1,j+1,k  ) +     &
                               u(i,j  ,k-1) + u(i-1,j  ,k-1) +     &
                               u(i,j+1,k-1) + u(i-1,j+1,k-1) )     &
                               * 0.125_ireals
            END DO
            END DO
            DO j = jstart, jend
            DO i = istart, iend
               zs12 = lays(i,j,1) - lays(i-1,j,1)
               zs21 = lays(i,j,2) - lays(i,j-1,2)

               zs12 = ( zs12 - lay(i,j)*vari(i,j,k,v_m) ) * acrlat(j,1) * eddlon
               zs21 = ( zs21 - frc(i,j)*vari(i,j,k,u_m) ) * edadlat

               hlp(i,j,k) = hlp(i,j,k) + ( zs12 + zs21 )**2
            END DO
            END DO

            IF (itype_sher.EQ.3 .OR. PRESENT(tket_hshr)) THEN 
               !Separate horizontale Scherungsmode soll berechnet werden:

               DO j = jstart, jend
               DO i = istart, iend
                  src(i,j)=(a_hshr*l_hori)**2 * hlp(i,j,k)**z3d2
               END DO
               END DO

               IF (PRESENT(tket_hshr)) THEN
                  DO j = jstart, jend
                  DO i = istart, iend
                     tket_hshr(i,j,k)=src(i,j)
                  END DO
                  END DO
               END IF    

               IF (itype_sher.EQ.3) THEN
                 !Korrektur durch separate horizontal Scherungsmode:
                  DO j = jstart, jend
                  DO i = istart, iend
                     hlp(i,j,k) = hlp(i,j,k) + src(i,j)/(tke(i,j,k,nvor)*tkvm(i,j,k))
                  END DO
                  END DO
               END IF    
            END IF   

            IF (itype_sher.EQ.2) THEN
               !Korrektur durch 3D Scherung der mittleren Strömung:

               !Vertikale Scherungskorrektur:
               DO j = jstart, jend
               DO i = istart-1, iend
                  lays(i,j,1) = z1d2 * ( w(i,j,k) + w(i+1,j,k) )
               END DO
               END DO
               DO j = jstart-1, jend
               DO i = istart, iend
                  lays(i,j,2) = z1d2 * ( w(i,j,k) + w(i,j+1,k) )
               END DO
               END DO
               DO j = jstart, jend
               DO i = istart, iend
                  zs13 = lays(i,j,1) - lays(i-1,j,1)
                  zs23 = lays(i,j,2) - lays(i,j-1,2)

                  zs33 = ( w(i,j,k+1) - w(i,j,k-1) ) &
                           / ( hhl(i,j,k-1) - hhl(i,j,k+1) )

                  zs13  = ( zs13 - lay(i,j)*zs33 ) * acrlat(j,1) * eddlon
                  zs23  = ( zs23 - frc(i,j)*zs33 ) * edadlat

                  hlp(i,j,k) = hlp(i,j,k) + zs13 * ( z2*vari(i,j,k,u_m) + zs13) &
                                          + zs23 * ( z2*vari(i,j,k,v_m) + zs23) &
                                          + z2 * zs33**2
               END DO
               END DO
            END IF

         END DO

         IF (itype_sher.EQ.2 .AND. kem.EQ.ke1) THEN

            !Neigungskorrektur der vertikalen Scherung am unteren Modellrand
            !bei 3D-Scherung :

            DO j = jstart, jend
            DO i = istart, iend
               !Neigung der Erdoberflaeche:
               lay(i,j) = ( hhl(i+1,j,ke1) - hhl(i-1,j,ke1) )*z1d2
               frc(i,j) = ( hhl(i,j+1,ke1) - hhl(i,j-1,ke1) )*z1d2

               zs33 = w(i,j,ke) / hhl(i,j,ke)
               zs13 = lay(i,j)*zs33
               zs23 = frc(i,j)*zs33

               hlp(i,j,ke1) = ( frc(i,j)*vari(i,j,ke1,u_m) - lay(i,j)*vari(i,j,ke1,v_m) )**2 &
                            + zs13 * ( z2*vari(i,j,ke1,u_m) + zs13 ) &
                            + zs23 * ( z2*vari(i,j,ke1,v_m) + zs23 )
            END DO
            END DO

         END IF

         !Beachte: Von der gesamten 3D-Scherung wurde die reine vertikale Scherung
         !         des Horizontalwindes ausgespart. Diese wird i.f. bestimmt.

      END IF
#endif
!__COSMO__---------------------------------------------------------------------------

! 3)  Hauptschleife: Bestimmung der TKE und der Stabilitaetsfunkt.:
!print *,"vor hauptschleife"

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

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Belegung mit Werten gueltig ausserhalb des Bestandes:

               dd(i,j,0)=d_m

               dd(i,j,1)=d_1
               dd(i,j,2)=d_2
               dd(i,j,3)=d_3
               dd(i,j,4)=d_4
               dd(i,j,5)=d_5
               dd(i,j,6)=d_6

               dd(i,j,7)=rim
            END DO
            END DO

         END IF

         DO k=2,kem

!           Berechnung der atmosphaerischen Vertikal-Antriebe:

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               lays(i,j,1)=vari(i,j,k,u_m)**2+vari(i,j,k,v_m)**2
               lays(i,j,2)=a(i,j,k,4)*vari(i,j,k,tet_l) &
                          +a(i,j,k,5)*vari(i,j,k,h2o_g)
 
               val1=tkvm(i,j,k)*lays(i,j,1); val2=tkvh(i,j,k)*lays(i,j,2)

               frc(i,j)=MAX( dd(i,j,7)*val1, val1-val2 )

!              Durch die MAX-Funktion wird frc so nach unten 
!              beschraenkt, dass im level-2 Gleichgewicht die
!              krit. Rizahl nie ueberschritten wird. In diesem
!              Gleichgewicht entsrpicht dies dann einem von der
!              Windscherung abhaengigen Minimalwert von TKE.

!              Die Stabilitaetsfunktionen werden aber unabhaengig von
!              dieser Beschraenkung berechnet!

!              Stabilitaetskorrektur der turbulenten Laengenskala
!              bei stabilier Schichtung:

               velh=tke(i,j,k,nvor)
               wert=a_stab*SQRT(MAX( z0, lays(i,j,2)) )

               len_scale(i,j,k)=velh*len_scale(i,j,k)/(velh+wert*len_scale(i,j,k))

            END DO
            END DO

!           Berechnung des Antriebs durch Nachlaufproduktion:

            IF (k.GE.kcm) THEN !Innerhalb des Bestandes:

               IF (k.EQ.ke1) THEN
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lay(i,j)=SQRT(wind(i,j,ke1,u_m)**2 &
                                  +wind(i,j,ke1,v_m)**2)
                  END DO
                  END DO
               ELSE
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     dpu=dp0(i,j,k)
                     dpo=dp0(i,j,k-1)
                     dp2=dpu+dpo

                     vel1=(wind(i,j,k,u_m)*dpo+wind(i,j,k-1,u_m)*dpu)/dp2
                     vel2=(wind(i,j,k,v_m)*dpo+wind(i,j,k-1,v_m)*dpu)/dp2
                     lay(i,j)=SQRT(vel1**2+vel2**2+w(i,j,k)**2)
                  END DO
                  END DO
               END If

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

!                 Hilfsgroesse zur Berechn. der reduzierten Konstanten:

                  wert=z3*csml(i,j,k)*len_scale(i,j,k) &
                                      *lay(i,j)/tke(i,j,k,nvor)

!                 Berechnung der modifizierten Modellparameter:

                  dd(i,j,0)=d_m/(z1+d_m*wert)

                  dd(i,j,1)=a_h/(z1+a_h*wert)
                  dd(i,j,2)=a_m/(z1+z2*a_m*wert)

                  dd(i,j,3)=z9*dd(i,j,1)
                  dd(i,j,4)=z6*dd(i,j,2)
                  dd(i,j,5)=z3*(d_h+dd(i,j,4))
                  dd(i,j,6)=dd(i,j,3)+z3*dd(i,j,4)
                  dd(i,j,1)=z1/dd(i,j,1)
                  dd(i,j,2)=z1/dd(i,j,2)

                  dd(i,j,7)=z1/(z1+(dd(i,j,0)-dd(i,j,4))/dd(i,j,5))

!                 TKE-Forcing durch Addition der Nachlaufproduktion zum vertikalen Forcing:

                  lay(i,j)=frc(i,j)+cbig(i,j,k)*lay(i,j)**3/tke(i,j,k,nvor)

               END DO
               END DO

            ELSE !Ausserhalb des Rauhigkeitsbestandes

!              TKE-Forcing gleich vertikales Forcing:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lay(i,j)=frc(i,j)
               END DO 
               END DO 

            END IF

            IF (.NOT.(lini.OR.k.eq.ke1)) THEN !nicht bei der Initialisierung oder am Unterrand

               IF (lsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
                  !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:

!                 Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO_Tendenzen:
          
                  DO j=jstart,jend
                  DO i=istart,iend
                     dpu=dp0(i,j,k)
                     dpo=dp0(i,j,k-1)
                     dp2=dpu+dpo

                     vel1=-(ut_sso(i,j,k)  *wind(i,j,k  ,u_m)      &
                       &  + vt_sso(i,j,k)  *wind(i,j,k  ,v_m))*dpo
                     vel2=-(ut_sso(i,j,k-1)*wind(i,j,k-1,u_m)      &
                       &  + vt_sso(i,j,k-1)*wind(i,j,k-1,v_m))*dpu

                     src(i,j)=MAX( z0, (vel1+vel2)/dp2 )

!                    Beachte: 
!                    Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
!                    erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
!                    Obwohl vel1 und vel2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!
                  END DO   
                  END DO   

                  IF (PRESENT(tket_sso)) THEN
                     DO j=jstart,jend
                     DO i=istart,iend
                        tket_sso(i,j,k)=src(i,j)
                     END DO
                     END DO
                  END IF

                  IF (ltkesso) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
                     DO j=jstart,jend
                     DO i=istart,iend
                        lay(i,j)=lay(i,j) + src(i,j)/tke(i,j,k,nvor)
                     END DO
                     END DO
                  END IF

               END IF

               IF (lconv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
                  !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

!                 Addition die TKE-Quelle durch konvektive Aktivitaet:

                  DO j=jstart,jend
                  DO i=istart,iend
                     lay(i,j)=lay(i,j) + MAX( z0, tket_conv(i,j,k)/tke(i,j,k,nvor) )
                  END DO
                  END DO

!                 Beachte:  Obwohl tket_conv immer positiv sein muesste, wird zur Sicherheit MAX benutzt!
               END IF
            END IF

!           Beruecksichtigung der 3D_Scherungskorrektur:

            IF (itype_sher.GT.1) THEN !3D-turbulence or separate horizontal shear

!              Addition der 3D-Scherungskorrektur zum bisherigen TKE-Forcing:

               DO j = jstart, jend
               DO i = istart, iend
                  lay(i,j)=lay(i,j)+tkvm(i,j,k)*hlp(i,j,k)
               END DO
               END DO

            END IF

!           Berechnung der neuen TKE-Werte:

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Berechnung einiger Hilfsgroessen:

               q0=tke(i,j,k,nvor)
               l_turb=len_scale(i,j,k)
               l_diss=dd(i,j,0)*l_turb

!              SQRT(2*TKE)-Prognose:

               q1=l_diss/dt_tke
               q2=MAX( z0, q0+tketens(i,j,k)*dt_tke)+lay(i,j )*dt_tke

!              In tketens steht die zuvor berechnete Advektions- und 
!              Diffusionstendenz (einschl. Drucktransport).
!              Die MAX-Funkt. verhindert, dass Transportterme die TKE
!              negativ machen koennen.

               q3=q1*(sqrt(z1+z4*q2/q1)-z1)*z1d2

               q1=SQRT(l_turb*frc(i,j)*dd(i,j,4)/(z1-c_g))
               tke(i,j,k,ntur)=MAX( vel_min, q1, q3*(z1-tkesmot)+q0*tkesmot )

!              Die MAX-Funktion und die zeitliche Glaettung
!              sind (vor allem wegen der Wirkung des Zirkulationstermes)
!              nur in Ausnahmefaellen noetig und koennen dann
!              vermeiden, dass die TKE zu klein wird.
!              q1 ist ein Minimalwert, um eine moegliche Singuaritaet bei der
!              Berechnung der Stabilitaetsfunktion fuer labile Schichtung
!              zu vermeiden.

!              Sichern des thermischen TKE-Antriebs:

               hlp(i,j,k)=lays(i,j,2)

            END DO    
            END DO    

            IF (PRESENT(edr)) THEN
!              Sichern der TKE-Dissipation "edr=q**3/l_diss" u.a. als thermische Quelle:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  edr(i,j,k)=tke(i,j,k,ntur)**3/(dd(i,j,0)*len_scale(i,j,k)) !mit gefiltertem tke-Feld
               END DO   
               END DO   
            END IF   

!           Berechnung der neuen stabilitaetsabhangigen Laengenskalen:

            IF (lstfnct) THEN

               CALL stab_funct(sm=tkvm(:,:,k), sh=tkvh(:,:,k), fm2=lays(:,:,1), fh2=lays(:,:,2), &
                               frc=frc, tvs=tke(:,:,k,ntur), tls=len_scale(:,:,k), &
                               i_st=istartpar,i_en=iendpar, j_st=jstartpar,j_en=jendpar)

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  tkvm(i,j,k)=len_scale(i,j,k)*tkvm(i,j,k)
                  tkvh(i,j,k)=len_scale(i,j,k)*tkvh(i,j,k)
               END DO
               END DO

            END IF
!print *,"nach hauptschleife"

!------------------------------------------------------------------------------------
#ifdef SCLM
            IF (lsclm .AND. it_durch.EQ.it_end) THEN

               cvar(tet_l,1)=d_h*tkvh(im,jm,k)*len_scale(im,jm,k) &
                                *vari(im,jm,k,tet_l)**2

               q3=tke(im,jm,k,ntur); q2=q3**2
               l_diss=dd(im,jm,0)*len_scale(im,jm,k); fakt=l_diss/q3

               km=tkvm(im,jm,k)*q3; kh=tkvh(im,jm,k)*q3

               cvar(u_m  ,2)=-km*vari(im,jm,k,u_m)
               cvar(v_m  ,2)=-km*vari(im,jm,k,v_m)
               cvar(tet_l,2)=-kh*vari(im,jm,k,tet_l)

               x1=cvar(u_m,2)*vari(im,jm,k,u_m)
               x2=cvar(v_m,2)*vari(im,jm,k,v_m)
               x3=-kh*lays(im,jm,2)

               !Achtung: TKE%mod(k)%val zeigt z.Z. noch auf die alte TKE-Zeitstufe.
               !Somit wird also der alte tke-Wert mit dem neuen ueberschrieben,
               !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
               !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

               TKE_SCLM%mod(k)%val=q3; TKE_SCLM%mod(k)%vst=i_cal

               BOYPR%mod(k)%val=x3       ; BOYPR%mod(k)%vst=i_cal
               SHRPR%mod(k)%val=-(x1+x2) ; SHRPR%mod(k)%vst=i_cal
               DISSI%mod(k)%val=q2/fakt  ; DISSI%mod(k)%vst=i_cal 
               TRANP%mod(k)%val=q3*tketens(im,jm,k)
                                           TRANP%mod(k)%vst=i_cal
                
               cvar(u_m,1)=z1d3*q2+fakt*(-z4*x1+z2*x2-z2*x3)
               cvar(v_m,1)=z1d3*q2+fakt*(+z2*x1-z4*x2-z2*x3)

               UWA%mod(k)%val=cvar(u_m  ,2); UWA%mod(k)%vst=i_cal
               VWA%mod(k)%val=cvar(v_m  ,2); VWA%mod(k)%vst=i_cal
               TWA%mod(k)%val=cvar(tet_l,2); TWA%mod(k)%vst=i_cal
               UST%mod(k)%val=SQRT(km*SQRT(lays(im,jm,1)))
                                           UST%mod(k)%vst=i_cal

               IF (cvar(u_m,1).GE.z0 .AND. cvar(v_m,1).GE.z0 .AND. &
                   cvar(u_m,1)+cvar(v_m,1).LE.q2) THEN

                  UUA%mod(k)%val=cvar(u_m,1); UUA%mod(k)%vst=i_cal
                  VVA%mod(k)%val=cvar(v_m,1); VVA%mod(k)%vst=i_cal
                  WWA%mod(k)%val=q2-UUA%mod(k)%val-VVA%mod(k)%val
                                            WWA%mod(k)%vst=i_cal
               END IF

               TTA%mod(k)%val=cvar(tet_l,1); TTA%mod(k)%vst=i_cal

            END IF   
#endif
!SCLM--------------------------------------------------------------------------------

            IF (ltmpcor) THEN

!              Sichern der TKE-Dissipation als thermische Quelle:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  tketens(i,j,k)=tke(i,j,k,ntur)**3/(dd(i,j,0)*len_scale(i,j,k))
               END DO
               END DO

            END IF

         END DO !k=2,kem

         IF (it_durch.LT.it_end) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF   

      END DO !Iterationen ueber it_durch

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         tke(i,j,1,ntur)=tke(i,j,2,ntur)
      END DO
      END DO

!-----------------------------------------------------------------------
#ifdef SCLM 
      IF (lsclm) THEN
         TKE_SCLM%mod(1)%val=tke(im,jm,1,ntur); TKE_SCLM%mod(1)%vst=i_cal
      END IF
#endif 
!SCLM--------------------------------------------------------------------
!print *,"nach iterat"

      IF (iini.EQ.1) THEN !only for separate initialization before the time loop

         DO k=2, kem
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tkvh(i,j,k)=tkvh(i,j,k)*tke(i,j,k,ntur)
               tkvm(i,j,k)=tkvm(i,j,k)*tke(i,j,k,ntur)
            END DO
            END DO
         END DO

       ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
       ! IF (.NOT.PRESENT(d_pat)) DEALLOCATE ( dpat,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )
 
         RETURN !finish this subroutine

      END IF   

!  4) Berechnung der effektiven turbulenten Vertikalgradienten
!     und weiterer Turbulenzgroessen:

         DO k=2,kem

            IF (ltmpcor) THEN

!              Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:
!              
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lay(i,j)=a(i,j,k,2)*vari(i,j,k,tet_l)-tet_g &
                          +zlhocp*vari(i,j,k,liq) ! vert. Temperaturgradient
               END DO   
               END DO   

!              Dies geschieht schon hier, weil im naechsten Schritt das Feld vari()
!              durch die effiktiven Gradienten ueberschrieben wird.
            END IF    

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               x2=a(i,j,k,3)*a(i,j,k,2)

               a(i,j,k,3) = (c_scld * rcld(i,j,k)) / (z1+rcld(i,j,k)*(c_scld-z1))

               rcld(i,j,k)=sqrt(len_scale(i,j,k)*tkvh(i,j,k)*d_h)* &
                           abs(x2*vari(i,j,k,tet_l)-vari(i,j,k,h2o_g))

!              Unter rcld steht jetzt die (geschaetzte) Standardabw.
!              des Saettigungsdefizites!

!              Berechn. der Diffusionskoeffizienten:

               tkvh(i,j,k)=MAX( con_h, tkhmin, tkvh(i,j,k)*tke(i,j,k,ntur) )
               tkvm(i,j,k)=MAX( con_m, tkmmin, tkvm(i,j,k)*tke(i,j,k,ntur) )

!              tkvh und tkvm enthalten jetzt nicht mehr Diffusions-
!              laengen, sondern Diffusionskoeffizienten in m^2/s!
            END DO
            END DO

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Berech. der thermodynamischen Flusskonversionsmatrix
!              und der Standardabweichnung der Verteilung des
!              Saettigungsdefizites:

               x1=z1/(z1+zlhocp*a(i,j,k,3))
               x2=a(i,j,k,3)*a(i,j,k,2)
               x3=a(i,j,k,2)*a(i,j,k,1) !exnr*Cp/Cpd

               rcl=a(i,j,k,3)

               flukon33         =z1-rcl*(z1-x1)
               flukon34         =x1*zlhocp/a(i,j,k,2)*rcl

               flukon53         =-x1*x2*rcl
               flukon54         =x1*rcl
               flukon55         =z1-liqfak

               flukon43         =-flukon53
               flukon44         =z1-flukon54

               flux_3   =(flukon33         *vari(i,j,k,tet_l) &
                         +flukon34         *vari(i,j,k,h2o_g)) * x3
               flux_4   =(flukon43         *vari(i,j,k,tet_l) &
                         +flukon44         *vari(i,j,k,h2o_g))
               flux_5   =(flukon53         *vari(i,j,k,tet_l) &
                         +flukon54         *vari(i,j,k,h2o_g) &
                         +flukon55         *vari(i,j,k,liq)) 

               vari(i,j,k,tem)=flux_3 !exnr*(Cp/Cpd)*grad(tet)
               vari(i,j,k,vap)=flux_4
               vari(i,j,k,liq)=flux_5

               !Achtung: 
               !'vari(i,j,k,tem)' ist der in der Temperaturgleichung benoetigte
               !effective Gradient und entspricht exnr*(Cp/Cpd)*tet_flux_dens.
               !Die Indices 'tet', 'tem' und 'tet_l' haben den gleichen Wert
               !und sollen nur den jeweiligen Inhalt verdeutlichen.

               !Im Falle "icldm_turb.NE.-1" ist "liqfak.EQ.z1" und 'flukon55' verschwindet.
               !Im Falle "icldm_turb.EQ.-1" dagegen verschwinden 'flukon53' und 'flukon54'.
               !Dann ist aber "flukon55.EQ.z1". Ferner verschwinden in diesem Fall auch
               !'flukon34' und 'flukon43', wobei die physikalischen Groessen 'tet_l' und 'tet'
               !sowie 'h2o_g' und 'vap' identisch sind.

!-----------------------------------------------------------------------
!     nur im alten 1-d_Modell:

!              Berechnung der der Flussdichten fuer Impuls,
!              fuehlbare und latente Energie (zum Erdboden hin posit.): 

!              km=-rhon(i,j,k)*tkvm(i,j,k)
!              kh=-rhon(i,j,k)*tkvh(i,j,k)

!              flimpuls(i,j,k)=km*sqrt(vari(i,j,k,u_m)**2  &
!                                     +vari(i,j,k,v_m)**2)
!              flfbwae(i,j,k)=cp_d*kh*vari(i,j,k,tem)
!              fllawae(i,j,k)=lh_v*kh*vari(i,j,k,vap)

!        Beachte: flimpuls,flfbwae,fllawael
!                 werden nur fuer Diagnosezwecke berechnet.
!-----------------------------------------------------------------------

            END DO
            END DO
           
            IF (ltmpcor) THEN

!              Berechnung der Temperaturtendenzen durch TKE-Quellen
!              auf Nebenflaechen (ausser der Divergenz des 
!              Drucktransportes):

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  thermik=tkvh(i,j,k)*hlp(i,j,k)
                  phasdif=tkvh(i,j,k)*lay(i,j) &
                         *(zrcpv*vari(i,j,k,vap)+zrcpl*vari(i,j,k,liq))  

                  tketens(i,j,k)=len_scale(i,j,k)/a(i,j,k,1) &
                          *((tketens(i,j,k)+thermik)/cp_d+phasdif)
               END DO   
               END DO   

!              Beachte:
!              In tketens() steht bisher die TKE-Dissipation.
!              Wegen der spaeteren Interpolation auf Hauptflaechen,
!              wird mit der turbulenten Laengenskala multipliziert.
            END IF
   
         END DO !k=1, kem


!  5) Untere Randwerte der Turbulenzgroessen:

         IF (kem.EQ.ke1) THEN

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               km=tkvm(i,j,ke1); kh=tkvh(i,j,ke1)
 
!              Dicke der laminaren Grenzschicht:
               l_lam=z0m(i,j)*SQRT(con_m/km)

!              Reduktionsfaktor:
               fakt=(l_lam/dicke(i,j,ke1))*(kh/con_h)
               tfh(i,j)=z1/(z1+rlam_heat*fakt)
               tfm(i,j)=z1/(z1+rlam_mom *fakt)
               tfv(i,j)=(z1+rlam_heat*fakt)/(z1+rat_lam*rlam_heat*fakt)
                        !Reduktionsfaktor fuer die Verdunstung aufgrund eines um den
                        !Faktor 'rat_lam' gegenueber fuehlbarer Waerme vergroesserten
                        !laminaren Transpostwiderstandes

               velh=MAX( vel_min, &
                         SQRT(wind(i,j,ke,u_m)**2+wind(i,j,ke,v_m)**2) )

!              Neue Transferkoeffizienten:

               fakt=z1/(velh*dicke(i,j,ke1))
               tcm(i,j)=fakt*km*tfm(i,j)
               tch(i,j)=fakt*kh*tfh(i,j)

!              Berechnung der neuen Rauhigkeitslaenge ueber Meer:

               IF (fr_land(i,j).LT.z1d2) THEN

                  ! Use ice surface roughness or open-water surface roughness
                  ! according to lo_ice
                  IF ( lo_ice(i,j) ) THEN
                     gz0(i,j)=grav*z0_ice
                  ELSE
!                    Berechnung der Schubspannung:

                     wert=tcm(i,j)*velh*SQRT(velh**2+tke(i,j,ke1,ntur)**2)
                     
!                    grav*z0 mittels Charnock-Formel: 

                     gz0(i,j)=MAX( grav*len_min, alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
                  END IF
               END IF

            END DO   
            END DO   

            IF (PRESENT(shfl_s)) THEN
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  shfl_s(i,j)=-cp_d*rhon(i,j,ke1)*tkvh(i,j,ke1)*vari(i,j,ke1,tem) &
                                                                  /a(i,j,ke1,2)
               END DO
               END DO
            END IF
            IF (PRESENT(lhfl_s)) THEN
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lhfl_s(i,j)=-lh_v*rhon(i,j,ke1)*tkvh(i,j,ke1)*vari(i,j,ke1,vap)
               END DO
               END DO
            END IF

!---------------------------------------------------------------------------------------
#ifdef SCLM 
            IF (lsclm) THEN

               val1=tkvm(im,jm,ke1)*SQRT(lays(im,jm,1)) !mech_forc(ke1)=lays(im,jm,1)
               val2=tkvh(im,jm,ke1)*lays(im,jm,2)       !ther_forc(ke1)=lays(im,jm,2)
      
               cvar(tet,2)=-tkvh(im,jm,ke1)*vari(im,jm,ke1,tem)/a(im,jm,ke1,2) !tet'w'
               cvar(vap,2)=-tkvh(im,jm,ke1)*vari(im,jm,ke1,vap)                !vap'w'
      
               wert=rhon(im,jm,ke1) !Luftdichte
      
               SHF%mod(0)%val=cp_d*wert*cvar(tet,2); SHF%mod(0)%vst=i_cal
               LHF%mod(0)%val=lh_v*wert*cvar(vap,2); LHF%mod(0)%vst=i_cal
             ! LMO%mod(0)%val=val1**3/val2         ; LMO%mod(0)%vst=i_cal

            END IF
#endif 
!SCLM-----------------------------------------------------------------------------------

         ELSE

!           Untere Randwerte der Turbulenzgroessen mit hilfe der
!           bereits vorher (in turbtran oder partur(b,s) berechneten
!           Transferkoeffizienten:

            IF (itype_tran.NE.2) THEN

!              Es wurde i.a. noch keine unteren Randwerte fuer 
!              tke(), sowie tkvm(), tkvh() berechnet:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  velh=MAX( vel_min, &
                           SQRT(wind(i,j,ke,u_m)**2+wind(i,j,ke,v_m)**2) )
                  wert=velh*dicke(i,j,ke1)
                  tkvh(i,j,ke1)=wert*tch(i,j)
                  tkvm(i,j,ke1)=wert*tcm(i,j)

!                 Weil in partur(b,s) keine laminare Grenzschicht
!                 beruecksichtigt wird, sind die Reduktionskoeff.
!                 immer gleich 1:

!                 tfh(i,j)=z1
!                 tfm(i,j)=z1

!                 Der SQRT(2*TKE)-Wert am unteren Modellrand wird mit 
!                 Hilfe von c_tke*Ustar bestimmt: 

                  tke(i,j,ke1,ntur)=c_tke*velh*SQRT(tcm(i,j))
               END DO   
               END DO   

            END IF
   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Thermischer Antrieb (fh2):

               hlp(i,j,ke1)=a(i,j,ke1,4)*vari(i,j,ke1,tet_l) &
                           +a(i,j,ke1,5)*vari(i,j,ke1,h2o_g)

!              Thermodynamische Korrektur des effektiven vertikalen 
!              Gradienten der pot. Temperatur am unteren Modellrand:

               vari(i,j,ke1,tem)=vari(i,j,ke1,tem) & !grad(tet)
                                *a(i,j,ke1,2)*a(i,j,ke1,1) !exnr*Cp/Cpd

               !Beachte:
               !Am Unterrand soll grad(liq) = 0 sein. Daher ist dort
               !grad(tet_l) = grad(tet).
               !'vari(i,j,ke1,tem)' ist der effective in der Temperaturgleichung
               !benoetigte Gradient und entspricht exnr*(Cp/Cpd)*tet_flux_dens.

               a(i,j,ke1,3)=rcld(i,j,ke1)
            END DO   
            END DO   

            IF (ltmpcor) THEN  

!              Bestimmung der unteren Randwerte der Enthalpieproduktion
!              durch Wechselwirkung mit der TKE (diss und thermik),
!              sowie durch Phasendiffusion:
          
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

                  wert=(t(i,j,ke)-tp(i,j))/dicke(i,j,ke1) !T-Grad.
                  thermik=tkvh(i,j,ke1)*hlp(i,j,ke1)
                  phasdif=tkvh(i,j,ke1)*wert &
                         *(zrcpv*vari(i,j,ke1,vap)+zrcpl*vari(i,j,ke1,liq))

                  wert=tke(i,j,ke1,ntur)**3/(d_m*len_scale(i,j,ke1))
                  tketens(i,j,ke1)=len_scale(i,j,ke1)/a(i,j,ke1,1) &
                         *((wert+thermik)/cp_d+phasdif)

!test             tketens(i,j,ke1)=tketens(i,j,ke)

               END DO   
               END DO   
                
            END IF    
         
         END IF 


! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperatur-
!     tendenzen ausser der Divergenz des Drucktransportes:

      IF (ltmpcor) THEN
!print *,"ltmpcor"
         DO j=jstart,jend
         DO i=istart,iend
            ttens(i,j,1)=ttens(i,j,1)+tinc(tem)*tketens(i,j,2) &
                         /(len_scale(i,j,1)+len_scale(i,j,2))

         END DO
         END DO
         DO k=2,ke
            DO j=jstart,jend
            DO i=istart,iend
               ttens(i,j,k)=ttens(i,j,k)+tinc(tem)*(tketens(i,j,k)+tketens(i,j,k+1)) &
                            /(len_scale(i,j,k)+len_scale(i,j,k+1))
            END DO
            END DO
         END DO
      END IF   


! 7)  Bestimmung des Drucktransporttermes (Zirculationstermes):

      lcircterm(:,:)=.FALSE.

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         IF (dpat(i,j).GT.z0) THEN
            l_pat(i,j)=l_hori*dpat(i,j)/(l_hori+dpat(i,j))
            lcircterm(i,j)=.TRUE.
         ELSE
            l_pat(i,j) = z0
         END IF
      END DO
      END DO


         DO k=2,ke1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tketens(i,j,k)=z0
            END DO
            END DO

            IF (k.LT.ke1) THEN
!              Interpolation des Druckes auf die naechst hoehere 
!              Nebenflaeche:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                 IF (lcircterm(i,j)) &
                  lay(i,j)=(prs(i,j,k)*dp0(i,j,k-1)+prs(i,j,k-1)*dp0(i,j,k)) &
                          /(dp0(i,j,k)+dp0(i,j,k-1))
               END DO   
               END DO   
            ELSE
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                 IF (lcircterm(i,j)) lay(i,j)=pr(i,j)
               END DO
               END DO
            END IF

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               IF (lcircterm(i,j)) THEN
               fakt=z1-z2*ABS(a(i,j,k,3)-z1d2)
               len=MAX( l_pat(i,j), SQRT(fakt*len_scale(i,j,k)*l_hori) )

!              Berechnung der lokalen Kohaerenzlaenge fuer die 
!              Drucktransport-Parameterisierung:
!              kohae_len=l_pat*exnr*grad(tet_v)*r/g:

               fakt=hlp(i,j,k)*lay(i,j)/(rhon(i,j,k)*grav**2)
               kohae_len=len*SIGN(z1,fakt)*MIN( ABS(fakt), z1 )

!              Belegung von tketens mit dem zur vert. turb. 
!              TKE-Flussdichte mittels Drucktransport gehoerigen
!              TKE-Gradienten:

               tketens(i,j,k)=kohae_len*hlp(i,j,k)

!              Die Divergenz des zugehoerigen Flusses ist gleichzeitig
!              eine weitere Quelle fuer thermische Energie.
               ENDIF
            END DO
            END DO
!print *,"k=",k
!print *,"len_scale=",MINVAL(len_scale(:,:,k)),MAXVAL(len_scale(:,:,k)), &
!         " tketens=",MINVAL(tketens(:,:,k)),  MAXVAL(tketens(:,:,k))

         END DO   
  


! 8)  Berechnung der korrigierten Flussdichten, die fuer die 
!     expliziten Berechnungen (Korrekturen) benoetigt werden 
!     und Ablegen derselben im Feld a(), sowie weitere Vorbereitungen:

      IF (lexplizit.or.limplizit) THEN

         wert=securi/dt_var

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            lay(i,j)=(rhon(i,j,ke)+rhon(i,j,ke1)) &
                    *(hhl(i,j,ke)-hhl(i,j,ke1))*z1d2
         END DO
         END DO

         IF (nk2.LE.ndiff) THEN
!           Berechnung der korrigierten Gradienten der thermodyn.
!           (skalaren) Variablen:

!           Beachte:
!           Cp/Cpd bleibt noch unter a(i,j,k,1) gespeichert.

            IF (imode_turb.EQ.2) THEN
!              Die gesamte Rechnung ist explizit, daher werden die
!              unkorrigierten effektiven Gradienten benutzt:

               DO n=nk2,ndiff
                  DO k=2,ke1
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,n)=vari(i,j,k,n)
                     END DO
                     END DO
                  END DO
               END DO
            ELSE   
!              Es wird ein Teil der Flussdichtedivergenzen implizit
!              berechnet. Daher werden nur korrigierte Gradienten
!              benutzt, so dass das explizite Verfahren nur den bisher 
!              nicht berechneten Flussdivergenzanteil uebernimmt:

!              Am unteren Rand ist keine explizite Korrektur noetig,
!              weil grad(liq)=0 ist, und hier kein Unterschied 
!              zwischen lokaler und nicht-lokaler Gradientberechnung
!              besteht:

               DO n=nk2,ndiff
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,ke1,n)=z0
                  END DO
                  END DO
               END DO

!              Fuer die darueberliegenden Schichten werden die lokalen
!              Gradienten der eigentlichen thermodyn. Modellvariablen
!              benoetigt:

               DO k=2,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     edh=z1/dicke(i,j,k)
                     a(i,j,k,tem)=(t(i,j,k-1)-t(i,j,k))*edh+tet_g
                     a(i,j,k,vap)=(qv(i,j,k-1)-qv(i,j,k))*edh
                     a(i,j,k,liq)=(qc(i,j,k-1)-qc(i,j,k))*edh
                  END DO
                  END DO
               END DO

               IF (imode_turb.EQ.3) THEN
!                 Die impliziten Anteile werden mit den thermodyn.
!                 Erhaltungsgroessen bestimmt.
!                 Fuer die explizite Korrektur werden die Differenzen
!                 zwischen den effektiven Gradienten und den lokalen
!                 Gradienten der zugehoerigen Erhaltungsgroessen
!                 benoetigt:

                  DO k=2,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,tem)=vari(i,j,k,tem) & !effectiver Temp.grad.
                         -(a(i,j,k,tem)-zlhocp*a(i,j,k,liq))*a(i,j,k,1) !eff. Grad. der Fluess.wass.temp.
                        a(i,j,k,vap)=vari(i,j,k,vap) &
                         -(a(i,j,k,vap)+liqfak*a(i,j,k,liq))
                        a(i,j,k,liq)=vari(i,j,k,liq)
                     END DO   
                     END DO   
                  END DO   
               ELSE   
!                 Die imliziten Anteile werden mit den eigentlichen
!                 thermodyn. Modellvariablen bestimmt.
!                 Fuer die explizite Korrektur werden die Differenzen
!                 zwischen den effektiven und den lokalen Gradienten
!                 dieser Variablen benoetigt:

                  DO n=nk2,ndiff
                     DO k=2,ke
                        DO j=jstartpar,jendpar
                        DO i=istartpar,iendpar
                           a(i,j,k,n)=vari(i,j,k,n)-a(i,j,k,n)
                        END DO
                        END DO
                     END DO
                  END DO
               END IF   
            END IF

!           Addition des Gradienten, welcher zur Temperatur-
!           flussdichte durch Drucktransport gehoert:

            IF (ltmpcor) THEN
               DO k=2,ke1 
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                    IF (lcircterm(i,j)) &
                     a(i,j,k,tem)=a(i,j,k,tem) &
                                  -tketens(i,j,k)/(cp_d*a(i,j,k,1))
                  END DO
                  END DO
               END DO
            END IF

!           Berechnung der zugehoerigen Flussdichten (mit aus numer.
!           Gruenden reduzierten Diffusionskoeff.):

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,1)=rhon(i,j,k) &
                              *MIN( wert*dicke(i,j,k)**2, tkvh(i,j,k) )
               END DO
               END DO
               DO n=nk2,ndiff
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,n)=-lays(i,j,1)*a(i,j,k,n)
                  END DO 
                  END DO 
               END DO 
            END DO 

!           Quasi impliziete Behandlung der untersten Schicht (Prandtl-
!           schicht):

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,1)=rhon(i,j,ke1)*tkvh(i,j,ke1)
               lays(i,j,2)=lay(i,j) &
                   /(lays(i,j,1)*dt_var/(dicke(i,j,ke1)*securi)+lay(i,j))
            END DO 
            END DO 
            DO n=nk2,ndiff
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

!                 Durch Anpassung der Flussdichte durch die ke-te
!                 Nebenflaeche, kann als untere Randbedingung der
!                 Fluss aus den Boden erhalten bleiben:

                  a(i,j,ke1,n)=-lays(i,j,1)*a(i,j,ke1,n)
                  a(i,j,ke,n)=a(i,j,ke1,n)*(z1-lays(i,j,2)) &
                             +a(i,j,ke,n)*lays(i,j,2)
               END DO
               END DO 
            END DO 

         END IF   

         IF (nk1.LE.nvel) THEN

!           Berechnung der korrigierten Gradienten der 
!           Impulskomponenten:

            IF (imode_turb.EQ.2) THEN
!              Die gesamte Rechnung ist explizit, daher werden die
!              unkorrigierten Gradienten benutzt:

               DO n=nk1,nvel
                  DO k=2,ke1
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,n)=vari(i,j,k,n)
                     END DO
                     END DO
                  END DO
               END DO
            ELSE   
!              Es wird ein Teil der Flussdichtedivergenzen implizit
!              berechnet. Daher werden nur korrigierte Gradienten
!              benutzt, so dass das explizite Verfahren nur den bisher 
!              nicht berechneten Flussdivergenzanteil uebernimmt.
!              Dies sind bei den Impulskomponenten nur die Effekte
!              der nicht-lokalen Gradienten

!              Am unteren Rand ist keine explizite Korrektur noetig,
!              weil hier kein Unterschied zwischen lokaler und 
!              nicht-lokaler Gradientberechnung besteht:

               DO n=1,nvel
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,ke1,n)=z0
                  END DO
                  END DO
               END DO

               DO n=nk1,nvel
                  DO k=2,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,n)=vari(i,j,k,n) &
                        -(wind(i,j,k-1,n)-wind(i,j,k,n))/dicke(i,j,k)
                     END DO
                     END DO
                  END DO
               END DO
            END IF   

!           Berechnung der zugehoerigen Flussdichten (mit aus numer.
!           Gruenden reduzierten Diffusionskoeff.):

            DO k=2,ke 
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,1)=rhon(i,j,k) &
                              *MIN( wert*dicke(i,j,k)**2, tkvm(i,j,k) )
               END DO   
               END DO   
               DO n=nk1,nvel
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,n)=-lays(i,j,1)*a(i,j,k,n)
                  END DO
                  END DO 
               END DO 
            END DO 

!           Quasi impliziete Behandlung der untersten Schicht (Prandtl-
!           schicht):

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,1)=rhon(i,j,ke1)*tkvm(i,j,ke1) 
               lays(i,j,2)=lay(i,j) &
                   /(lays(i,j,1)*dt_var/(dicke(i,j,ke1)*securi)+lay(i,j))
            END DO   
            END DO   
            DO n=nk1,nvel 
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

!                 Durch Anpassung der Flussdichte durch die ke-te
!                 Nebenflaeche, kann als untere Randbedingung der
!                 Fluss aus den Boden erhalten bleiben:

                  a(i,j,ke1,n)=-lays(i,j,1)*a(i,j,ke1,n)
                  a(i,j,ke,n)=a(i,j,ke1,n)*(z1-lays(i,j,2)) &
                             +a(i,j,ke,n)*lays(i,j,2)
               END DO
               END DO 
            END DO 

         END IF   

      END IF

!     Effektive Schichtdicke multipliziert mit der Dichte:

      IF (lexplizit.OR.limplizit) THEN
         kk=1
      ELSE
         kk=kcm
      END IF

      DO k=kk,ke
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            dicke(i,j,k)=(rhon(i,j,k)+rhon(i,j,k+1)) &
                        *(hhl(i,j,k)-hhl(i,j,k+1))*z1d2
         END DO
         END DO
      END DO


! 9)  Explizite Berechnung von Bestandtendenzen:

      IF (kcm.LE.ke) THEN
!print *,"Explizite Berechnung von Bestandtendenzen"

!        Berechnung des Formwiderstandes im Bestand:

         ku=1
         ko=2
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            lays(i,j,ku)=z0
         END DO
         END DO
         DO k=ke,kcm,-1   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,ko)=w(i,j,k)
               wert=z1d2*(cbig(i,j,k)  +csml(i,j,k) &
                         +cbig(i,j,k+1)+csml(i,j,k+1)) &
                   *SQRT(wind(i,j,k,u_m)**2+wind(i,j,k,v_m)**2 &
                        +(z1d2*(lays(i,j,ku)+lays(i,j,ko)))**2)

!              Durch die min-Funktion, wird verhindert, dass durch
!              die Wirkung der Formreibung eine Richtungsumkehr
!              der Geschwindigkeit erfolgen kann:

               wert=MIN( z1/dt_var,wert )

               wind(i,j,k,u_m)=-wert*wind(i,j,k,u_m)*tinc(u_m)
               wind(i,j,k,v_m)=-wert*wind(i,j,k,v_m)*tinc(v_m)
            END DO   
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
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,ku)=-rhon(i,j,ke1)*tkvm(i,j,ke1) &
                                             *vari(i,j,ke1,n)
               END DO
               END DO
            ELSE
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,ku)=-rhon(i,j,ke1)*tkvh(i,j,ke1) &
                                             *vari(i,j,ke1,n)
               END DO
               END DO
            END IF

            DO k=ke,kcm,-1
               IF (kcm.EQ.1) THEN
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lays(i,j,ku)=lays(i,j,ku)*z1d2 &
                         *(rair(i,j,k)-rair(i,j,k+1))/dicke(i,j,k)
                  END DO
                  END DO
               ELSE   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lays(i,j,ko)=-rhon(i,j,k)*tkvh(i,j,k)*vari(i,j,k,n)
                     lays(i,j,ku)=(lays(i,j,ko)+lays(i,j,ku))*z1d2 &
                         *(rair(i,j,k)-rair(i,j,k+1))/dicke(i,j,k)
                  END DO
                  END DO
               END IF
               IF (n.LE.nvel) THEN
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     wind(i,j,k,n)=wind(i,j,k,n)-lays(i,j,ku)*tinc(n)
                  END DO
                  END DO
               ELSEIF (n.EQ.tem) THEN
                  DO j=jstart,jend
                  DO i=istart,iend
                     ttens(i,j,k)=ttens(i,j,k)-lays(i,j,ku)*tinc(tem)
                  END DO
                  END DO
               ELSEIF (n.EQ.vap) THEN
                  DO j=jstart,jend
                  DO i=istartpar,iendpar
                     qvtens(i,j,k)=qvtens(i,j,k)-lays(i,j,ku)*tinc(vap)
                  END DO
                  END DO
               ELSEIF (n.EQ.liq) THEN
                  DO j=jstart,jend
                  DO i=istart,iend
                     qctens(i,j,k)=qctens(i,j,k)-lays(i,j,ku)*tinc(liq)
                  END DO
                  END DO
               END IF
               kk=ku
               ku=ko
               ko=kk 
            END DO    
         END DO    

      END IF   

! 10) Berechnung der Flussdichtedivergenzen:

! 10a) Vorbereitungen: 

      IF (lexplizit.OR.limplizit) THEN

!        Im Feld wind() stehen jetzt innerhalb der Bestandes die dort
!        anfallenden Tendenzen der horizontalen Windkomponenten.
!        I.f. werden weitere Wind-Tendenzen auf dieses Feld addiert.
!        Daher muss es ausserhalb des Bestandes auf 0 gesetzt werden:
 
         DO k=1,kcm-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               wind(i,j,k,u_m)=z0     
               wind(i,j,k,v_m)=z0     
            END DO
            END DO
         END DO

      END IF    

! 10b) Berechnung des expliziten Anteiles der Flussdichtedivergenzen:

       IF (lexplizit) THEN
!print *,"expliziter Anteil"

!         Explizite Berechnung der Flussdichtedivergenzen:

            DO n=nkorr,ndiff
              
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,1)=-a(i,j,2,n)
               END DO
               END DO
               DO k=2,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     hlp(i,j,k)=a(i,j,k,n)-a(i,j,k+1,n)
                  END DO
                  END DO
               END DO

!              Die folgende (quellenfreie) vertikale Glaettung ist
!              bei groesseren Zeitschritten (> ca. 30sec) in der
!              expliziten Variante noetig:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,1,n)=-(hlp(i,j,1)+wichfakt*(-hlp(i,j,1) &
                              +hlp(i,j,2)))/dicke(i,j,1)
               END DO
               END DO
               DO k=2,ke-1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,n)= &
                      -(hlp(i,j,k)+wichfakt*(hlp(i,j,k-1)-2._ireals*hlp(i,j,k) &
                                         +hlp(i,j,k+1)))/dicke(i,j,k)
                  END DO
                  END DO
               END DO
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke,n)=-(hlp(i,j,ke)+wichfakt*(hlp(i,j,ke-1) &
                               -hlp(i,j,ke)))/dicke(i,j,ke)
               END DO
               END DO

!              Beruecksichtigung der vertikalen Cp-Variationen bei
!              der Berechnung der Temperaturtendenzen

               IF (n.eq.tem) THEN
                  DO k=1,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,tem)=a(i,j,k,tem) &
                         /(z1+zrcpv*qv(i,j,k)+zrcpl*qc(i,j,k))
                     END DO
                     END DO
                  END DO
!                 Beachte:
!                 Der Einfachheit halber wurde hier der Cp-Faktor
!                 ohne Berucksichtigung der subskaligen Wolken-
!                 wasserkorrektur berechnet!
               END IF   

            END DO   

!           Abspeichern der Tendenzen in den Tendenzfeldern, damit
!           das Feld a() fuer die ev. folgende implizite Rechnung
!           verfuegbar ist:

            DO n=nk1,nvel
               DO k=1,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     wind(i,j,k,n)=wind(i,j,k,n)+a(i,j,k,n)*tinc(n)
                  END DO
                  END DO
               END DO
            END DO  

            DO k=1,ke
!print *,"k=",k," tt_expl_vor=",MINVAL(ttens(:,:,k)),MAXVAL(ttens(:,:,k))
               DO j=jstart,jend
               DO i=istart,iend
                  ttens(i,j,k) = ttens(i,j,k)+a(i,j,k,tem)*tinc(tem)
                  qvtens(i,j,k)=qvtens(i,j,k)+a(i,j,k,vap)*tinc(vap)
                  qctens(i,j,k)=qctens(i,j,k)+a(i,j,k,liq)*tinc(liq)
               END DO
               END DO
!print *,"k=",k," tt_expl_nac=",MINVAL(ttens(:,:,k)),MAXVAL(ttens(:,:,k))

               IF (PRESENT(qvt_diff)) THEN
                  ! ... expliziten Anteil fuer Berechnung von dqvdt
                  !     nochmals separat sichern
                  DO j=jstart,jend
                  DO i=istart,iend
                     qvt_diff(i,j,k) = a(i,j,k,vap)
                  END DO
                  END DO
               END IF   
            END DO

      END IF   

! 10c) Berechnung des impliziten Anteiles der Flussdichtedivergenzen:

      IF (limplizit) THEN
!print *,"impliziter Anteil"

         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=dicke(i,j,k)/dt_var
            END DO
            END DO
         END DO

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            a(i,j,1,1)=z0
         END DO
         END DO

         DO n=1,nvel

!           Bestimmung der Tendenzen der horiz. Windgeschw.:

            DO k=2,ke 
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,1)=tkvm(i,j,k) &
                         *z2*rhon(i,j,k)/(hhl(i,j,k-1)-hhl(i,j,k+1))
               END DO
               END DO
            END DO
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,ke1,1)=tkvm(i,j,ke1) &
                           *rhon(i,j,ke1)/dicke(i,j,ke1)
            END DO
            END DO

            IF (n.EQ.u_m) THEN
               DO k=1,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     ii=MAX( i-vst ,1 )
                     hlp(i,j,k)=(u(i,j,k)+u(ii,j,k))*z1d2
                  END DO
                  END DO
               END DO   
            ELSEIF (n.EQ.v_m) THEN
               DO k=1,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     jj=MAX( j-vst, 1 )
                     hlp(i,j,k)=(v(i,j,k)+v(i,jj,k))*z1d2
                  END DO
                  END DO
               END DO   
            END IF   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,ke1)=hlp(i,j,ke)*(z1-tfm(i,j))
            END DO
            END DO

            DO k=1,ke-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,3)=a(i,j,k+1,1)
                  a(i,j,k,2)=-(a(i,j,k,1)+a(i,j,k,3)+dicke(i,j,k))
                  a(i,j,k,4)=-dicke(i,j,k)*hlp(i,j,k)
               END DO
               END DO
            END DO
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,ke,2)=-(a(i,j,ke,1)+a(i,j,ke1,1) &
                                        +dicke(i,j,ke))
               a(i,j,ke,4)=-dicke(i,j,ke)*hlp(i,j,ke)
            END DO
            END DO

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,1)=a(i,j,k,1)/a(i,j,k-1,2)
                  a(i,j,k,2)=a(i,j,k,2)-a(i,j,k,1)*a(i,j,k-1,3)
                  a(i,j,k,4)=a(i,j,k,4)-a(i,j,k,1)*a(i,j,k-1,4)
               END DO
               END DO
            END DO

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,ke,4)=a(i,j,ke,4)/a(i,j,ke,2)
            END DO
            END DO
            DO k=ke-1,1,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,4)=(a(i,j,k,4)-a(i,j,k,3)*a(i,j,k+1,4)) &
                                                    /a(i,j,k,2) 
               END DO
               END DO
            END DO

            DO k=1,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  wind(i,j,k,n)=wind(i,j,k,n) &
                               +(a(i,j,k,4)-hlp(i,j,k))*tinv(n)
               END DO
               END DO
            END DO

         END DO   

         DO n=nred,nvel+1,-1

!           Bestimmung der Tendenzen der reduzierten scalaren Var..
!           Dabei soll tet_l als letztes behandelt werden:

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,1)=tkvh(i,j,k) &
                         *z2*rhon(i,j,k)/(hhl(i,j,k-1)-hhl(i,j,k+1))
               END DO
               END DO
            END DO
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,ke1,1)=tkvh(i,j,ke1) &
                           *rhon(i,j,ke1)/dicke(i,j,ke1)
            END DO
            END DO

            IF (n.EQ.tet_l) THEN
               ko=1
               ku=2
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  exnr=exner(i,j,1)
                  wert=z1+zrcpv*qv(i,j,1)+zrcpl*qc(i,j,1)

                  hlp(i,j,1)=(t(i,j,1) &
                             -zlhocp*qc(i,j,1))/exnr
                  dicke(i,j,1)=dicke(i,j,1)*wert
                  lays(i,j,ko)=exnr*wert
               END DO
               END DO

               DO k=2,ke   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     exnr=exner(i,j,k)
                     wert=z1 &
                         +zrcpv*qv(i,j,k)+zrcpl*qc(i,j,k)

                     hlp(i,j,k)=(t(i,j,k) &
                                -zlhocp*qc(i,j,k))/exnr
                     dicke(i,j,k)=dicke(i,j,k)*wert
                     lays(i,j,ku)=exnr*wert

                     wert=(lays(i,j,ku)*dp0(i,j,k-1) &
                          +lays(i,j,ko)*dp0(i,j,k)) &
                         /(dp0(i,j,k-1)+dp0(i,j,k))
                     a(i,j,k,1)=a(i,j,k,1)*wert
                  END DO   
                  END DO   
                  kk=ko
                  ko=ku
                  ku=kk
               END DO   
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
!                 Beachte:
!                 tp(), qd(), ql() und pr() sind fuer die 
!                 Schicht k=ke1 bereits vorhanden.
     
                  exnr=zexner(pr(i,j))
                  wert=z1+zrcpv*qd(i,j)+zrcpl*ql(i,j)

                  hlp(i,j,ke1)=(tp(i,j)-zlhocp*ql(i,j))/exnr

                  a(i,j,ke1,1)=a(i,j,ke1,1)*exnr*wert
!lay(i,j)=exnr
               END DO
               END DO
!do k=1, ke
! print *,"k=",k
! print *,"t_vor=",MINVAL(t(:,:,k)),MAXVAL(t(:,:,k))
! print *,"qv_vor=",MINVAL(qv(:,:,k)),MAXVAL(qv(:,:,k))
! print *,"qc_vor=",MINVAL(qc(:,:,k)),MAXVAL(qc(:,:,k))
! print *,"pr_vor=",MINVAL(prs(:,:,k)),MAXVAL(prs(:,:,k))
! print *,"exnr=",MINVAL(exner(:,:,k)),MAXVAL(exner(:,:,k))
! print *,"hlp=",MINVAL(hlp(:,:,k)),MAXVAL(hlp(:,:,k))
!enddo
! print *,"k=",ke1
! print *,"t_vor=",MINVAL(tp),MAXVAL(tp)
! print *,"qv_vor=",MINVAL(qd),MAXVAL(qd)
! print *,"qc_vor=",MINVAL(ql),MAXVAL(ql)
! print *,"pr_vor=",MINVAL(pr),MAXVAL(pr),"p0ref=",p0ref
! print *,"exnr=",MINVAL(lay),MAXVAL(lay)
! print *,"hlp=",MINVAL(hlp(:,:,ke1)),MAXVAL(hlp(:,:,ke1))

!              Weil tet_l als letztes behandelt wird, konnte
!              dicke() mit der speziellen Korrektur bzgl. 
!              cp_fakt ueberschrieben werden!
            ELSEIF (n.EQ.h2o_g) THEN
               DO k=1,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     hlp(i,j,k)=qv(i,j,k)+liqfak*qc(i,j,k)
                  END DO
                  END DO
               END DO   
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
!                 Beachte:
!                 qd() und ql() sind fuer die Schicht k=ke1
!                 bereits vorhanden.

                  hlp(i,j,ke1)=qd(i,j)+liqfak*ql(i,j)
               END DO
               END DO
            END IF    

            DO k=1,ke-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,3)=a(i,j,k+1,1)
                  a(i,j,k,2)=-(a(i,j,k,1)+a(i,j,k,3)+dicke(i,j,k))
                  a(i,j,k,4)=-dicke(i,j,k)*hlp(i,j,k)
               END DO
               END DO
            END DO

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,ke,2)=-(a(i,j,ke,1)+dicke(i,j,ke))
               a(i,j,ke,4)=-dicke(i,j,ke)*hlp(i,j,ke) &
                          +a(i,j,ke1,1)*(hlp(i,j,ke)-hlp(i,j,ke1))
            END DO
            END DO

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,k,1)=a(i,j,k,1)/a(i,j,k-1,2)
                  a(i,j,k,2)=a(i,j,k,2)-a(i,j,k,1)*a(i,j,k-1,3)
                  a(i,j,k,4)=a(i,j,k,4)-a(i,j,k,1)*a(i,j,k-1,4)
               END DO
               END DO
            END DO

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               vari(i,j,ke,n)=a(i,j,ke,4)/a(i,j,ke,2)
            END DO
            END DO
            DO k=ke-1,1,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vari(i,j,k,n)=(a(i,j,k,4)-a(i,j,k,3)*vari(i,j,k+1,n)) &
                                                      /a(i,j,k,2)
               END DO
               END DO
            END DO

            IF (imode_turb.EQ.3) THEN
!              Es werden lediglich die Diffusionstendenzen der Erhaltungsvariablen
!              zurückberechnet und auf die Tendenzfelder der Modellvariablen gelegt,
!              so dass noch entsprechende explizite Korrekturen nötig sind:

               DO k=1, ke
!if (n.EQ.tet_l) then
! print *,"k=",k," tt_impl_vor=",MINVAL(ttens(:,:,k)),MAXVAL(ttens(:,:,k))
! print *,"hlp=",MINVAL(hlp(:,:,k)),MAXVAL(hlp(:,:,k))," vari=",MINVAL(vari(:,:,k,n)),MAXVAL(vari(:,:,k,n))
!endif
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     vari(i,j,k,n)=(vari(i,j,k,n)-hlp(i,j,k))*tinv(n)
                  END DO
                  END DO
               END DO

            ELSEIF (imode_turb.EQ.4) THEN
!              Für die Konversion in Modellvariablen wird die Fluessigwassertemperatur
!              benötigt:

               IF (n.EQ.tet_l) THEN
                  DO k=1,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        vari(i,j,k,n)=vari(i,j,k,n)*exner(i,j,k)
                     END DO
                     END DO
                  END DO
               END IF
            END IF   

         END DO

         IF (imode_turb.EQ.4) THEN
!           Die um ihre impliziten Diffusionstendenzen erhöhten thermodyn. Erhaltungsgroessen
!           sollen in Modellvariablen konvertiert werden, woraus dann deren Diffusionstendenzen
!           zurueckberechnet werden:

!do k=1, ke
! print *,"vor cloud k=",k
! print *,"exner=",MINVAL(exner(:,:,k)),MAXVAL(exner(:,:,k))
! print *,"prs=",MINVAL(prs(:,:,k)),MAXVAL(prs(:,:,k))
! print *,"rcld=",MINVAL(rcld(:,:,k)),MAXVAL(rcld(:,:,k))
! print *,"tem_l=",MINVAL(vari(:,:,k,tem_l)),MAXVAL(vari(:,:,k,tem_l))
! print *,"h2o_g=",MINVAL(vari(:,:,k,h2o_g)),MAXVAL(vari(:,:,k,h2o_g))
!end do

            CALL turb_cloud (                                         &
#ifdef __ICON__
                 itype_wcld,                                          & 
#endif
                 ie, je, ke,ke1,           1,ke,                      &
                 istartpar,iendpar, jstartpar,jendpar, 1,ke,          &
                 prs, ps, rcld, vari(:,:,:,tem_l), vari(:,:,:,h2o_g), &
                 clc=hlp, clwc=vari(:,:,:,liq))

            DO k=1,ke
!print *,"nach cloud k=",k
!print *,"tem_old=",MINVAL( t(:,:,k)),MAXVAL( t(:,:,k))
!print *,"vap_old=",MINVAL(qv(:,:,k)),MAXVAL(qv(:,:,k))
!print *,"liq_old=",MINVAL(qc(:,:,k)),MAXVAL(qc(:,:,k))
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vari(i,j,k,tem)=vari(i,j,k,tem_l)+zlhocp*vari(i,j,k,liq)
                  vari(i,j,k,vap)=vari(i,j,k,h2o_g)-liqfak*vari(i,j,k,liq)

                   ttens(i,j,k)= ttens(i,j,k)+(vari(i,j,k,tem)- t(i,j,k))*tinv(n)
                  qvtens(i,j,k)=qvtens(i,j,k)+(vari(i,j,k,vap)-qv(i,j,k))*tinv(n)
                  qctens(i,j,k)=qctens(i,j,k)+(vari(i,j,k,liq)-qc(i,j,k))*tinv(n)
               END DO
               END DO
!print *,"tem_new=",MINVAL(vari(:,:,k,tem)),MAXVAL(vari(:,:,k,tem))
!print *,"vap_new=",MINVAL(vari(:,:,k,vap)),MAXVAL(vari(:,:,k,vap))
!print *,"liq_new=",MINVAL(vari(:,:,k,liq)),MAXVAL(vari(:,:,k,liq))
            END DO

         ELSEIF (imode_turb.EQ.3) THEN
            !Es werden lediglich die Tendenzen der Erhaltungsvariablen auf die Temperatur-
            !und Fechte-Tendenzfelder addiert. Zugehörige explizite Kondensationskorrekturen
            !liefern insbesondere auch die Diffusionstendenz des Wolkenwassers:

            DO k=1,ke
               DO j=jstart,jend
               DO i=istart,iend
                  ttens(i,j,k) = ttens(i,j,k)+vari(i,j,k,tet_l)
                  qvtens(i,j,k)=qvtens(i,j,k)+vari(i,j,k,h2o_g)
               END DO
               END DO
!print *,"k=",k," tt_impl_nac=",MINVAL(ttens(:,:,k)),MAXVAL(ttens(:,:,k))
            END DO

         END IF    

         !Beachte, dass 'tet_l', 'tem_l' und 'tem' alle den gleichen Indexwert besitzen!

      END IF   

! 11) Berechnung der Tendenzen infloge horizontaler turb. Diffusion
!     (und aufsummieren auf die Tendenzfelder):

!----------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN

!        Nun werden die noch unter vari(:,:,:,1:2) abgespeicherten effektiven
!        Gradienten nur noch fuer die Berechnung der Komponenten
!        des turbulenten Spannungstensors benoetigt, welche u.U.
!        auch als Eingangsfeld fuer die Bestimmung der turbulenten
!        Horizontaldiffusion gebraucht werden:
         
         DO k=2,ke1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               cvar(u_m,2)=-tkvm(i,j,k)*vari(i,j,k,u_m)
               cvar(v_m,2)=-tkvm(i,j,k)*vari(i,j,k,v_m)

               q2=z1d3*tke(i,j,k,ntur)**2
               fakt=d_m*len_scale(i,j,k)/tke(i,j,k,ntur)
               x1=cvar(u_m,2)*vari(i,j,k,u_m)
               x2=cvar(v_m,2)*vari(i,j,k,v_m)
               x3=-tkvh(i,j,k)*hlp(i,j,k)

               vari(i,j,k,1)=cvar(u_m,2)
               vari(i,j,k,2)=cvar(v_m,2)
               vari(i,j,k,3)=q2+fakt*(-z4*x1+z2*x2-z2*x3)
               vari(i,j,k,4)=q2+fakt*(+z2*x1-z4*x2-z2*x3)
            END DO
            END DO
         END DO

!        Jetzt werden die in hlp() gespeicherten fh2-Werte nicht mehr
!        benoetigt und hlp() kann anderweitig benutzt werden.
         
      END IF
#endif
!SCLM------------------------------------------------------------------

!        ** wird erst spaeter eingefuehrt **

! 12) Aktualisierung der Tendenzfelder fuer den horizontalen Wind:

      IF (limplizit .OR. nk1.EQ.1 .OR. kcm.LE.ke) THEN

!-------------------------------------------------------------------------------------
#ifdef __COSMO__
!        Interpolation der Diffusionstendenzen von horizontalen
!        Impulskomponenten auf die 'gestaggerten' Positionen:

         IF (num_compute.GT.1 .AND. vst.NE.0) THEN

           IF (.NOT.PRESENT(ntstep)) THEN
              ierrstat=1
              errormsg= &
              'ERROR *** ''ntstep'' needed for ''exchg_boundaries'' is not present! ***'
              lerror=.TRUE.; RETURN
           END IF   

           CALL exchg_boundaries                                            &
               ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, ie, je,    &
                (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),      &
                jstartpar, jendpar, 2, nboundlines, my_cart_neigh,          &
                2200+ntstep, .FALSE., ncomm_type, ierrstat, errormsg,       &
                wind(:,:,:,u_m), wind(:,:,:,v_m) )
           IF (ierrstat.NE.0) THEN
              lerror=.TRUE.; RETURN
           END IF
         ENDIF
#endif
!__COSMO__----------------------------------------------------------------------------

         IF (limplizit.OR.(nk1.EQ.1)) THEN
            kk=1
         ELSE
            kk=kcm
         END IF

         DO k=kk,ke   
            DO j=jstartu,jendu     
            DO i=istartu,iendu   
               utens(i,j,k)=utens(i,j,k) &
                           +(wind(i,j,k,u_m)+wind(i+vst,j,k,u_m))*z1d2
            END DO
            END DO
            DO j=jstartv,jendv     
            DO i=istartv,iendv   
               vtens(i,j,k)=vtens(i,j,k) &
                           +(wind(i,j,k,v_m)+wind(i,j+vst,k,v_m))*z1d2
            END DO
            END DO
         END DO

      END IF   
                  
! 13) Berechnung der Diffusionstendenz von SQRT(2*TKE):

!        Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!        des Saettigungsdefizites und den Drucktransport:

            DO k=2,ke1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar  
                 IF (lcircterm(i,j)) &
                  tketens(i,j,k)=rhon(i,j,k)*tkvh(i,j,k)/tke(i,j,k,ntur) &
                                            *tketens(i,j,k)*len_scale(i,j,k)
               END DO
               END DO
            END DO
            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar  
                 IF (lcircterm(i,j)) &
                  tketens(i,j,k)=(tketens(i,j,k)+tketens(i,j,k+1)) &
                                /(len_scale(i,j,k)+len_scale(i,j,k+1))
               END DO
               END DO
            END DO

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar  
            rcld(i,j,1)=rcld(i,j,2)
         END DO
         END DO
         DO k=2,kem-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar  
               rcld(i,j,k)=(rcld(i,j,k)+rcld(i,j,k+1))*z1d2
            END DO
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

!<fuo_impltkediff
!        If limpltkediff=.false. we compute the original explicit diffusion
!        of TKE. Note that in order to achieve stability with this
!        explicit scheme, the following condition 0 <= securi < 0.5
!        should be met. Also, for very shallow model levels and TKE
!        values on the order of 1 m2/s2, the diffusion of TKE is
!        strongly limited by stability constraints, and an implicit
!        solution (by choosing limpltkediff=.true.) should be favored.
!fuo_impltkediff>
!        Nevertheless, we need for the present at least an exlicit correction
!        in order to express the circulation term and roughness layer terms.

!        Es soll keinen turb. TKE-FLuss aus der Atm. hinaus geben:

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            tketens(i,j,1)=z0
         END DO
         END DO

         IF (.not.limpltkediff) THEN

!           Bei der reinen explizite Diffusion von TKE muss die explizite
!           Flussdichte berechnet werden:

            DO k=ke,2,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  dicke(i,j,k)=hhl(i,j,k)-hhl(i,j,k+1)
   
                  q3=(tke(i,j,k,ntur)+tke(i,j,k+1,ntur))*z1d2
                  kh=q3*(len_scale(i,j,k)+len_scale(i,j,k+1))*z1d2
                  kh=(rhon(i,j,k)+rhon(i,j,k+1))*z1d2 &
                    *MIN( securi*dicke(i,j,k)**2/dt_tke, c_diff*kh )
                     
!-----------------------------------------------------------------------
!                    nur im 1-d-Modell:
   
!                    if (ldruck) then
!                       write (7,*) 'k,tke,p_flux,tke_flux' &
!                      ,k,tke(i,j,k,ntur),tketens(i,j,k)    &
!                      ,kh*q3*(tke(i,j,k+1,ntur)-tke(i,j,k,ntur)) &
!                            /dicke(i,j,k)
!                       end if
!-----------------------------------------------------------------------

                  tketens(i,j,k)=tketens(i,j,k) &
                                +kh*(tke(i,j,k+1,ntur)-tke(i,j,k,ntur))&
                                   /dicke(i,j,k)
               END DO
               END DO
            END DO

!           Bem.: Der viskosen Zusatzterm kann auch entfallen!

!           Die Diffusionskoeffizienten wurden mit Hilfe der MIN-Funktion
!           auf das numerisch ertraegliche Mass reduziert.

!           Zuvor wurde tketens bereits mit dem Drucktransport belegt,
!           so dass tketens also die Summe aus der eigentlichen TKE-
!           Flussdichte und dem Drucktrnasport enthaelt!

!           Positiv definite Diffusionskorrektur in der Gleichung fuer SQRT(2*TKE)
!           wobei die Laengenskala c_diff*len_scale beschraenkt werden muss:

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  len=hhl(i,j,k-1)-hhl(i,j,k+1)
                  hlp(i,j,k)=(tketens(i,j,k-1)-tketens(i,j,k))                &
                    & +MIN( securi*len**2/(4._ireals*tke(i,j,k,ntur)*dt_tke), &
                    &       c_diff*len_scale(i,j,k))                          &
                    & *( (tke(i,j,k-1,ntur)-tke(i,j,k+1,ntur))/len )**2
 
               END DO
               END DO
            END DO

         ELSE

!           Bei der impliziten Diffusion von TKE muss nur die Flussdichte bzueglich des
!           Zirkulationstermes explizit behandelt werden:

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,k)=tketens(i,j,k-1)-tketens(i,j,k) 
               END DO
               END DO
            END DO

         END IF    

!        Volumenkorrektur in der Rauhigkeitsschicht:

         DO k=MAX( 2, kcm ), ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,k)=hlp(i,j,k) &
                  +z1d2*(tketens(i,j,k)*dp0(i,j,k-1)+tketens(i,j,k-1)*dp0(i,j,k)) &
                       /(dp0(i,j,k)+dp0(i,j,k-1))*(rair(i,j,k-1)-rair(i,j,k+1))
            END DO
            END DO
         END DO

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            hlp(i,j,1)=hlp(i,j,2)
            hlp(i,j,ke1)=hlp(i,j,ke)
         END DO
         END DO

!        Die folgende (quellenfreie) vertikale Glaettung (mittels 'wichfakt')
!        ist bei groesseren Zeitschritten (> ca. 30sec) in der expliziten
!        Variante noetig:

         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=z1d2*(hhl(i,j,k-1)-hhl(i,j,k+1))*rhon(i,j,k)
               tketens(i,j,k)=                                                                &
                 &  -( hlp(i,j,k)+wichfakt*(hlp(i,j,k-1)-2._ireals*hlp(i,j,k)+hlp(i,j,k+1)) ) &
                 &  /dicke(i,j,k)
            END DO
            END DO
         END DO

         IF (limpltkediff) THEN
!        *************** semi-implicit TKE diffusion tendency *************

         ! define implicit/explicit weight
         !          explicit : beta = 0.0
         !   Crank-Nicholson : beta = 0.5
         !          implicit : beta = 1.0
         beta=0.85_ireals

         ! compute level thickness at full level height
         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,k)=hhl(i,j,k)-hhl(i,j,k+1)
            END DO
            END DO
         END DO

         ! compute level thickness at half level height
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            dicke(i,j,1)=hlp(i,j,1)*rhon(i,j,1)
         END DO
         END DO
         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=z1d2*(hhl(i,j,k-1)-hhl(i,j,k+1))*rhon(i,j,k)
            END DO
            END DO
         END DO
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            dicke(i,j,ke1)=hlp(i,j,ke)*rhon(i,j,ke1)
         END DO
         END DO

         ! compute diffusion constant at full level height
         ! note: kh = c_diff * rhon * len_scale * q**2
         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,k,1)=c_diff * z1d2 * (rhon(i,j,k)+rhon(i,j,k+1)) &
                                 * z1d2 * (len_scale(i,j,k)+len_scale(i,j,k+1)) &
                                 *(z1d2 * (tke(i,j,k,ntur)+tke(i,j,k+1,ntur)) )**2
               ! limiter to mimic old scheme
!               a(i,j,k,1)=MIN( a(i,j,k,1), &
!                    securi*hlp(i,j,k)**2/dt_tke/(z1d2*(rhon(i,j,k)+rhon(i,j,k+1))) )
            END DO
            END DO
         END DO

         ! Top layer (k=1)
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            ! setup matrix (a,b,c) and vector (d)
            tmp1       = dt_tke / tke(i,j,1,ntur) / dicke(i,j,1)
            tmp3       = tmp1 * a(i,j,1,1) / hlp(i,j,1)
            aa         = z0
            bb         = z0 !- tmp3
            a(i,j,1,3) = z0 !+ tmp3
            a(i,j,1,4) = + (z1+(z1-beta)*bb        )*tke(i,j,1,ntur) &
                         + (   (z1-beta)*a(i,j,1,3))*tke(i,j,2,ntur)
            bb         = z1 - beta*bb
            a(i,j,1,3) =    - beta*a(i,j,1,3)
            ! forward elimination
            a(i,j,1,3) = a(i,j,1,3)/bb
            a(i,j,1,4) = a(i,j,1,4)/bb
            a(i,j,1,2) = bb
         END DO
         END DO

         ! The layers from k=2 to k=ke1-1
         DO k=2,ke1-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               ! setup matrix (a,b,c) and vector (d)
               tmp1       = dt_tke / tke(i,j,k,ntur) / dicke(i,j,k)
               tmp2       = tmp1 * a(i,j,k-1,1) / hlp(i,j,k-1)
               tmp3       = tmp1 * a(i,j,k  ,1) / hlp(i,j,k  )
               aa         = + tmp2
               bb         = - tmp2 - tmp3
               a(i,j,k,3) = + tmp3
               a(i,j,k,4) = + (   (z1-beta)*aa        )*tke(i,j,k-1,ntur) &
                            + (z1+(z1-beta)*bb        )*tke(i,j,k  ,ntur) &
                            + (   (z1-beta)*a(i,j,k,3))*tke(i,j,k+1,ntur)
               aa         =    - beta*aa
               bb         = z1 - beta*bb
               a(i,j,k,3) =    - beta*a(i,j,k,3)
               ! forward elimination
               tmp1       = aa/a(i,j,k-1,2)
               a(i,j,k,4) = a(i,j,k,4) - tmp1*a(i,j,k-1,4)*a(i,j,k-1,2)
               a(i,j,k,2) = bb         - tmp1*a(i,j,k-1,3)*a(i,j,k-1,2)
               a(i,j,k,3) = a(i,j,k,3)/a(i,j,k,2)
               a(i,j,k,4) = a(i,j,k,4)/a(i,j,k,2)
            END DO
            END DO
         END DO

         ! The bottom layer (k=ke1)
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            ! setup matrix (a,b,c) and vector (d)
            tmp1         = dt_tke / tke(i,j,ke1,ntur) / dicke(i,j,ke1)
            tmp2         = tmp1 * a(i,j,ke1-1,1) / hlp(i,j,ke1-1)
            aa           = z0 !+ tmp2
            bb           = z0 !- tmp2
            a(i,j,ke1,3) = z0
            a(i,j,ke1,4) = + (   (z1-beta)*aa       )*tke(i,j,ke1-1,ntur) &
                           + (z1+(z1-beta)*bb       )*tke(i,j,ke1  ,ntur)
            aa           =    - beta*aa
            bb           = z1 - beta*bb
            ! forward elimination
            tmp1         = aa/a(i,j,ke1-1,2)
            a(i,j,ke1,4) = a(i,j,ke1,4) - tmp1*a(i,j,ke1-1,4)*a(i,j,ke1-1,2)
            a(i,j,ke1,2) = bb           - tmp1*a(i,j,ke1-1,3)*a(i,j,ke1-1,2)
            a(i,j,ke1,4) = a(i,j,ke1,4)/a(i,j,ke1,2)
         END DO
         END DO

         ! Backsubstitution and storage of new SQRT(2*TKE) field
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            hlp(i,j,ke1) = a(i,j,ke1,4)
         END DO
         END DO
         DO k=ke1-1,1,-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,k) = a(i,j,k,4) - a(i,j,k,3) * hlp(i,j,k+1)
            END DO
            END DO
         END DO

         ! Computation of SQRT(2*TKE) tendencies
         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tketens(i,j,k)=tketens(i,j,k)+ &
                              (hlp(i,j,k)-tke(i,j,k,ntur))/dt_tke
            END DO
            END DO
         END DO

!        *************** end of TKE diffusion tendency *************
         END IF

!        Im skin-layer soll es keine TKE-Tendenzen infolge
!        Diffusion oder Drucktransport geben:

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            tketens(i,j,ke1)=z0
         END DO
         END DO

!        Die folgende MAX-Funktion soll sicherstellen, dass allein durch die Diffusionstendenz
!        von SQRT(2*TKE) (einschl. des Zirkulationstermes) die TKE nicht negativ werden kann.

         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tketens(i,j,k)=MAX( -tke(i,j,k,ntur)/dt_tke, tketens(i,j,k) )
            END DO
            END DO
         END DO

!        Jetzt enthaelt tketens die SQRT(2*TKE)-Tend. infolge Diffusion
!        und Drucktransport. Zu tketens wird nun an anderer Stelle noch
!        die Advektionstendenz hinzuaddiert und das ganze dann im
!        naechsten Zeitschritt fuer die SQRT(2*TKE)-Prognose benutzt.

! 14) Berechnung der turb. horizontalen TKE-Flussdichtedivergenzen
!     und Aufsummieren auf das SQRT(2*TKE)-Tendenzfeld:

!        ** wird erst spaeter eingefuehrt **

! 15) Deallocierung lokaler dynamischer Felder:

    ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
    ! IF (.NOT.PRESENT(d_pat)) DEALLOCATE ( dpat,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )

END SUBROUTINE turbdiff

!********************************************************************************

!+ Module procedure stab_funct in "src_turbdiff" for computing the stability
!+ function for vertical diffusion


SUBROUTINE stab_funct (sm, sh, fm2, fh2, frc, tvs, tls, i_st,i_en, j_st,j_en)

   IMPLICIT NONE

   INTEGER (KIND=iintegers) :: i,j !loop indices

   REAL (KIND=ireals)        , INTENT(INOUT) :: &
        sm(:,:),  & !stablility function for momentum             [1] 
        sh(:,:)     !stablility function for scalars (heat)       [1] 

   REAL (KIND=ireals)        , INTENT(IN) :: &
        fm2(:,:), & !squared forcing frequency for momentum       [1/s2]
        fh2(:,:), & !squared forcing frequency for scalars (heat) [1/s2]
        frc(:,:), & !forcing function for TKE                     [m/s2]
        tvs(:,:), & !turbulent velocity scale SQRT(2*TKE)         [m/s]
        tls(:,:)    !turbulent length scale

   INTEGER (KIND=iintegers), INTENT(IN) :: &
        i_st, i_en, & !start- and end index of horizontal i-loop
        j_st, j_en    !start- and end index of horizontal j-loop

   REAL (KIND=ireals)          :: &
         d1, d2, d3, d4 ,d5 ,d6

   REAL (KIND=ireals) :: &
        a11, a12, a22, a21, &
        a3, a5, a6, be1, be2, &
        gama, fakt, val1, val2, &
        gm, gh

!------------------------------------------------------------------

   DO j=j_st, j_en
   DO i=i_st, i_en

      d1=dd(i,j,1); d2=dd(i,j,2); d3=dd(i,j,3)
      d4=dd(i,j,4); d5=dd(i,j,5); d6=dd(i,j,6)

      gama=tls(i,j)*frc(i,j)/tvs(i,j)**2 !entspr. 1/d_m im Gleichgewicht
                                         !und ausserh. des Bestandes
!test gama=1/d0

      fakt=(tls(i,j)/tvs(i,j))**2 

!     Folgnde Fallunterscheidung muss gemacht werden,
!     um positiv definete Loesungen fuer die Stabilitaets-
!     funktionen zu ermoeglichen:

      IF (fh2(i,j).GE.z0) THEN ! stab. Schichtung
!        Allgemeinste der hier verwendeten Loesungen:

         be1=z1
         be2=be1-c_g
       
         gh=fh2(i,j)*fakt
         gm=fm2(i,j)*fakt
         
         a11=d1+(d5-d4)*gh
         a12=d4*gm
         a21=(d6-d4)*gh
         a22=d2+d3*gh+d4*gm

         fakt=a11*a22-a12*a21
         sh(i,j)=(be1*a22-be2*a12)/fakt
         sm(i,j)=(be2*a11-be1*a21)/fakt
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

         val1=(fm2(i,j)*be2+(a5-a3+be1)*fh2(i,j))/(z2*be1)
         val2=val1+sqrt(val1**2-(a6+be2)*fh2(i,j)*fm2(i,j)/be1)
         fakt=fh2(i,j)/(val2-fh2(i,j))
         sh(i,j)=be1-a5*fakt
         sm(i,j)=sh(i,j)*(be2-a6*fakt)/(be1-(a5-a3)*fakt)
      END IF   

   END DO
   END DO

END SUBROUTINE stab_funct

!********************************************************************************

END SUBROUTINE organize_turbdiff

!********************************************************************************
!********************************************************************************


SUBROUTINE turb_cloud (                      &
#ifdef __ICON__
   itype_wcld,                               &
#endif
   ie, je, ke, ke1,            kcs,    kce,  &
   istart, iend, jstart, jend, kstart, kend, &
!
   prs, ps, rcld, t, qv, qc,                 &
!
   clc, clwc                                   )
   

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
#ifdef __ICON__
INTEGER (KIND=iintegers), INTENT(IN) :: &  ! type of surface-atmosphere transfer
  &  itype_wcld
#endif
INTEGER (KIND=iintegers), INTENT(IN) :: &  ! indices used for allocation of arrays
  ie, je,            & ! number of grib points of horizontal dimensions
  ke, ke1,           & ! number of full and half level dimension
  kcs, kce,          & ! lower and upper level index bound for output variabls
  istart, iend, & ! zonal      start and end index
  jstart, jend, & ! meridional start and end index
  kstart, kend    ! vertical   start and end index

! Array arguments with intent(in):

REAL (KIND=ireals), INTENT(IN) :: & !
  prs (ie,je,ke ),    &    ! base state pressure (")
  rcld(ie,je,ke1),    &    ! standard deviation of saturation deficit
  ps  (ie,je),        &    ! surface pressure
  t   (ie,je,ke ),    &    ! temperature (main levels)
  qv  (ie,je,ke )          ! water vapour (")

REAL (KIND=ireals), OPTIONAL, INTENT(IN) :: & !
  qc  (ie,je,ke )          ! cloud water  (")

! Array arguments with intent(out):

REAL (KIND=ireals), INTENT(INOUT) :: &
  clc (ie,je,kcs:kce),  & ! stratiform subgrid-scale cloud cover
  clwc(ie,je,kcs:kce)     ! liquid water content of ""

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,j,k ! loop indices

REAL (KIND=ireals), PARAMETER :: &
  zsig_max = 1.0E-3_ireals,  & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_ireals,   & ! fraction of saturation specific humidity
  zuc      = 0.95_ireals       ! constant for critical relative humidity

REAL (KIND=ireals) :: &
  qt(ie,je), tl(ie,je), qs(ie,je), gam(ie,je) , &
  pres, ql, dq, q, sig, uc, & !
  zsigma, zclc1, zq_max       !

!------------------------------------------------------------------
!Festlegung der Formelfunktionen fuer die Turbulenzparametrisierung:

!-------------------------------------------------------------------
#ifdef __COSMO__
!INCLUDE 'statement_functs.incf'
!Declaration of statement functions:

!REAL (KIND=ireals) :: &
 
!   zpsat_w ,zqvap, zdqsdt, &  !thermodynamic statement functions and
!   zpres, ztemp, zpvap, zqsat !their formal arguments
 
!Definition of statement functions:
 
!saturation vapour pressure over water:
!zpsat_w(ztemp)=b1*exp(b2w*(ztemp-b3)/(ztemp-b4w))
 
!specific humidity:
!zqvap(zpvap,zpres)=rdv*zpvap/(zpres-o_m_rdv*zpvap)
 
!derivation of zqsat with respect to temperature:
!zdqsdt(ztemp,zqsat)=b234w*(1.0_ireals+rvd_m_o*zqsat) &
!                         *zqsat/(ztemp-b4w)**2
#endif 
!__COSMO------------------------------------------------------------
 
!------------ End of header ----------------------------------------

! Begin Subroutine turb_cloud
! ---------------------------

  zq_max = q_crit*(z1/clc_diag - z1)

  DO k = kstart, kend

     IF (PRESENT(qc)) THEN
        DO j = jstart, jend
        DO i = istart, iend
           qt(i,j) = qc(i,j,k) +       qv(i,j,k) ! total water content
           tl(i,j) =  t(i,j,k) - lhocp*qc(i,j,k) ! liquid water temperature
        END DO 
        END DO 
     ELSE !qv and t already contain the conservation variables
        DO j = jstart, jend
        DO i = istart, iend
           qt(i,j) = qv(i,j,k)
           tl(i,j) =  t(i,j,k)
        END DO
        END DO
     END IF

     DO j = jstart, jend
     DO i = istart, iend
        qs(i,j) = zqvap( zpsat_w( tl(i,j) ), prs(i,j,k) )        ! saturation mixing ratio
       gam(i,j) = z1 / ( z1 + lhocp*zdqsdt( tl(i,j), qs(i,j) ) ) ! slope factor
     END DO
     END DO

     DO j = jstart, jend
     DO i = istart, iend

        pres = prs(i,j,k)        ! pressure
        dq   = qt(i,j) - qs(i,j) ! saturation deficit

        IF ( itype_wcld .EQ. 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i,j)

          ! critical relative humidity
          uc = zuc - uc1 * zsigma * ( z1 - zsigma )  &
                         * ( z1 + uc2*(zsigma-0.5_ireals) )

          ! cloud cover
          clc(i,j,k) = MAX( z0,  &
                            MIN( z1, clc_diag * ((qt(i,j)/qs(i,j)-uc)/(ucl-uc))**2 ) )

          ! in-cloud water content
          ql = qs(i,j) * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_ireals ) THEN
            zclc1 = clc_diag * ( (z1-uc)/(ucl-uc) )**2
            ql    = ql + (gam(i,j)*dq-ql)*(clc(i,j,k)-zclc1)/(z1-zclc1)
          END IF
          clwc(i,j,k) = clc(i,j,k) * ql

        ELSEIF ( itype_wcld .EQ. 2 ) THEN

        ! Statistical calculation of cloud cover and cloud water content
        ! using the standard deviation of the saturation deficit

          sig = MIN ( zsig_max, rcld(i,j,k) )

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_ireals ) THEN
            clc(i,j,k)  = ABS ( (SIGN(z1,dq)+z1)*0.5_ireals )
            clwc(i,j,k) = clc(i,j,k) * gam(i,j) * dq
          ELSE
            q = dq/sig
            clc(i,j,k) = MIN ( z1, MAX ( z0, clc_diag * ( z1 + q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              clwc(i,j,k) = z0
            ELSEIF ( q >= zq_max ) THEN
              clwc(i,j,k) = gam(i,j) * dq
            ELSE
              clwc(i,j,k) = gam(i,j) * sig * ( q + q_crit ) &
                                           * ( q + zq_max ) / ( z2*( q_crit + zq_max) )
            ENDIF
          ENDIF

        ENDIF

     ENDDO
     ENDDO
  ENDDO

END SUBROUTINE turb_cloud

!********************************************************************************

FUNCTION zpsat_w (ztemp)

  REAL (KIND=ireals), INTENT(IN) :: ztemp
  REAL (KIND=ireals) :: zpsat_w

  zpsat_w=b1*exp(b2w*(ztemp-b3)/(ztemp-b4w))

END FUNCTION zpsat_w

FUNCTION zqvap (zpvap, zpres)

  REAL (KIND=ireals), INTENT(IN) :: zpvap, zpres
  REAL (KIND=ireals) :: zqvap

  zqvap=rdv*zpvap/(zpres-o_m_rdv*zpvap)

END FUNCTION zqvap

FUNCTION zdqsdt (ztemp, zqsat)

  REAL (KIND=ireals), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=ireals) :: zdqsdt

  zdqsdt=b234w*(1.0_ireals+rvd_m_o*zqsat) &
              *zqsat/(ztemp-b4w)**2

END FUNCTION zdqsdt

!********************************************************************************
!********************************************************************************

END MODULE src_turbdiff
