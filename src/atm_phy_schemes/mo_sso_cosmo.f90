!>
!! Source module for the sub-grid scale orography
!!
!!-----------------------------------------------------------------------------
!!
!! @par Description:
!!   The module "mo_sso_cosmo" performs calculations related to the
!!   parameterization of sub-grid scale orographic (SSO) effects. The present
!!   SSO scheme is based on Lott and Miller (1997).
!!
!! @par For COSMO:
!!   All global variables of the model that are used by the SSO routines are
!!   imported by USE statements below. The interface of the SSO routines and
!!   the model is provided by the organizational routine "organize_sso".
!! @par For ICON:
!!   The subroutine "sso" is directly called in the interface routine
!!   "mo_nh_interface", the variables are transferred via a variable list.
!!
!! @par reference  The parameterization package has been extracted from the DWD
!!   limited area model COSMO (V4_16 of the DWD lm_f90 library) and is
!!   originally based on an earlier version by M. Miller and F. Lott (ECMWF).
!!
!! @par reference  Lott, F. and M. J. Miller, 1997: A new subgrid-scale
!!   orographic drag parametrization: Its formulation and testing.
!!   Q. J. R. Meteor. Soc., 123, 101-127.
!!
!! @author Francois Lott, LMD
!! @author Jan-Peter Schulz, DWD
!!
!! Current Code Owner: DWD, Jan-Peter Schulz
!!   phone:  +49  69  8062 2751
!!   fax:    +49  69  8062 3721
!!   email:  jan-peter.schulz@dwd.de
!!
!!
!! @par Revision History:
!! Version      Date       Name
!! ------------ ---------- ----
!! V1_0         2011/??/?? Jan-Peter Schulz, DWD
!!  Initial release
!!
!!
!! @par Code Description:
!! Language: Fortran 90.
!! Software Standards: "European Standards for Writing and
!! Documenting Exchangeable Fortran 90 Code".
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
!==============================================================================

MODULE mo_sso_cosmo

! Declarations:
!
! Modules used:

#ifdef __COSMO__

USE data_parameters , ONLY :   &
    ireals     ,    & ! KIND-type parameter for real variables
    iintegers         ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ! number of grid points for this domain
    ie         ,    & ! number of grid points in zonal direction
    je         ,    & ! number of grid points in meridional direction
    ke         ,    & ! number of grid points in vertical direction
    ke1        ,    & ! ke + 1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    
!   zonal direction
    istart     ,    & ! start index for the forecast of w, t, qd, qw and pp
    iend       ,    & ! end index for the forecast of w, t, qd, qw and pp

!   meridional direction
    jstart     ,    & ! start index for the forecast of w, t, qd, qw and pp
    jend       ,    & ! end index for the forecast of w, t, qd, qw and pp

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt         ,    & ! long time-step
    dt2               ! dt*2.

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &
!
! 0. physical constants and related variables
! -------------------------------------------
    r_d        ,    & ! gas constant for dry air
    cp_d       ,    & ! specific heat of dry air at constant pressure
    g                 ! acceleration due to gravity

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :    &
!
! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! reference pressure at full levels             ( Pa  )
    p0hl       ,    & ! reference pressure at half levels             ( Pa  )
    dp0        ,    & ! reference pressure thickness of layer         ( Pa  )
    hhl        ,    & ! geometrical height of half levels             ( m   )

! 2. external parameter fields                                        (unit)
! ----------------------------
    hsurf      ,    & ! height of surface topography                  ( m   )
    sso_stdh   ,    & ! standard deviation of sub-grid scale orography( m   )
    sso_gamma  ,    & ! anisotropy of sub-grid scale orography          --
    sso_theta  ,    & ! angle betw. principal axis of orography and E ( rad )
    sso_sigma  ,    & ! mean slope of sub-grid scale orography          --

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps         ,    & ! surface pressure                              ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields from the sub-grid scale orography scheme
    ut_sso     ,    & ! u-tendency due to SSO                         ( m/s2)
    vt_sso     ,    & ! v-tendency due to SSO                         ( m/s2)
    tt_sso     ,    & ! temperature tendency due to SSO               ( K/s )
    ustr_sso   ,    & ! u-stress (surface momentum flux) due to SSO   ( N/m2)
    vstr_sso   ,    & ! v-stress (surface momentum flux) due to SSO   ( N/m2)
    vdis_sso          ! vert. int. dissipation of kin. en. due to SSO ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :    &
!
! 1. start and end of the forecast
! --------------------------------
                      ! indices for permutation of three time levels
    nold       ,    & ! corresponds to ntstep - 1
    nnow       ,    & ! corresponds to ntstep

! 5. additional control variables
! -------------------------------
    l2tls             ! forecast with 2-TL integration scheme

! end of data_runcontrol

#endif

#ifdef __ICON__

USE mo_kind               , ONLY :   &
    ireals=>wp    , & ! KIND-type parameter for real variables
    iintegers=>i4     ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE mo_physical_constants , ONLY :   &
    r_d   => rd   , & ! gas constant for dry air
    cp_d  => cpd  , & ! specific heat of dry air at constant pressure
    g     => grav     ! acceleration due to gravity

! end of mo_physical_constants

#endif

!==============================================================================

IMPLICIT NONE

PRIVATE

!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: sso

!==============================================================================

! Declarations

! The following variables are tunable parameters for the sub-grid scale
! orography scheme. They are global in this module.

REAL (KIND = ireals) ::      &
!
! Tunable parameters
! ------------------
! Gsigcr  = 0.80_ireals   , &   ! top layer for low level drag
! Gkdrag  = 0.30          , &   ! gw drag constant (Original ECMWF value)
  Gkdrag  = 0.075_ireals  , &   ! gw drag constant
  Gkwake  = 0.50_ireals   , &   ! low level wake drag constant
! Gkdrag  = Gkdrag_read   , &   ! Gkdrag_read read in or set in gme_tuning_constants
! Gkwake  = Gkwake_read   , &   ! Gkwake_read read in or set in gme_tuning_constants
  Grcrit  = 0.25_ireals   , &   ! critical Richardson number
  Gfrcrit = 0.50_ireals   , &   ! critical Froude number

! Security constants
! ------------------
  Gvsec   = 0.10_ireals   , &   ! to secure the projection calculation
  Gssec   = 1.E-12_ireals , &   ! to secure stability
  Gtsec   = 1.E-07_ireals       ! to secure the stress calculation

!==============================================================================
! Module procedures in "mo_sso_cosmo"
!==============================================================================

CONTAINS

#ifdef __COSMO__

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE organize_sso

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure organize_sso is the interface of the model to the
!   parameterization package for sub-grid scale orographic effects.
!
! Externals:
!
!   sso: Lott and Miller SSO scheme
!
!------------------------------------------------------------------------------

! Local scalars and automatic arrays (also to be used in lower level routines):
! -----------------------------------------------------------------------------

  ! Input for the SSO routine "sso"
  ! -------------------------------
  REAL    (KIND=ireals   )  ::  &
    zt   (ie,je,ke),    & ! temperature at full levels
    zu   (ie,je,ke),    & ! zonal wind component
    zv   (ie,je,ke),    & ! meridional wind component
    zfif (ie,je,ke),    & ! geopotential at full levels
    zpf  (ie,je,ke),    & ! pressure at full levels
    zph  (ie,je,ke1),   & ! pressure at half levels
    zdt                   ! actual timestep for 2 or 3 TL-scheme     
  LOGICAL   ::  &
    ldebug                ! debug indicator

  INTEGER (KIND=iintegers) ::  &
    i, j, k, km1, nx

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  &
    izerror

  CHARACTER (LEN=80)       ::  &
    yzerrmsg
!
! End of header
!==============================================================================

  izerror  = 0
  yzerrmsg = '   '
  ldebug   = .FALSE.

  ! Select timelevel and timestep of the computation
  IF ( l2tls ) THEN
    nx  = nnow
    zdt = dt
  ELSE
    nx  = nold
    zdt = dt2
  ENDIF

  ! In order to save some CPU-time, the SSO scheme is called at fixed time
  ! increments nincsso (e.g. every 10 timsteps) and stores the SSO tendencies
  ! and other fields on global arrays which are held fixed in the intermediate
  ! steps. The time increment nincsso can be set on NAMELIST input.

  ! Prepare the input arrays for the SSO scheme.
  DO k = 1, ke 
    km1 = MAX ( 1, k-1 )
    DO j = jstart, jend
      DO i = istart, iend
        zt   (i,j,k) = t(i,j,k,nx)
        zu   (i,j,k) = 0.5_ireals*( u(i,j,k,nx) + u(i-1,j,k,nx) )
        zv   (i,j,k) = 0.5_ireals*( v(i,j,k,nx) + v(i,j-1,k,nx) )
        zfif (i,j,k) = 0.5_ireals*g*( hhl(i,j,k) + hhl(i,j,k+1) )
        zpf  (i,j,k) = p0(i,j,k) + pp(i,j,k,nx)
        zph  (i,j,k) = p0hl(i,j,k)  &
                     + 0.5_ireals*(pp(i,j,k,nx) + pp(i,j,km1,nx))
      ENDDO
    ENDDO
  ENDDO
  DO j = jstart, jend
    DO i = istart, iend
      zph  (i,j,ke1) = ps(i,j,nx)
    ENDDO
  ENDDO

  ! Call to sso.

  CALL sso (                                                 &
       ie    , je    , ke    , ke1     ,                     &
       istart, iend  , jstart, jend    ,                     &
       zpf   , zph   , zfif  , zt      , zu , zv , g*hsurf , &
       sso_stdh, sso_gamma, sso_theta, sso_sigma,            &
       zdt   ,                                               &
       ldebug,                                               &
       ut_sso, vt_sso, tt_sso, ustr_sso, vstr_sso, vdis_sso  )

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE organize_sso

#endif

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE sso (                                                       &
           ie     , je     , ke     , ke1      ,                       &
           istart , iend   , jstart , jend     ,                       &
           ppf    , pph    , pfif   , pt       , pu , pv  , pfis     , &
           psso_stdh, psso_gamma, psso_theta, psso_sigma,              &
           pdt    ,                                                    &
           ldebug ,                                                    &
           pdu_sso, pdv_sso, pdt_sso, pustr_sso, pvstr_sso, pvdis_sso  )

!------------------------------------------------------------------------------
!
! Purpose:   The module procedure sso performs the parameterisation of
!            sub-grid scale orographic effects according to Lott and
!            Miller (1997). It computes the physical tendencies of the
!            prognostic variables u, v and T due to vertical transports
!            by sub-grid scale orographically excited gravity waves.
!
! Method:    The scheme consists of three parts, i.e.
!             - the calculation of lowlevel drag due to wake effects of
!               the sub-grid scale orography
!             - the gravity wave stress and
!             - the stress profile in the vertical
!            The stress is computed based on the low level wind, the 
!            static stability and subgrid orographic parameters (i.e.
!            standard deviation of height, angle of the main orographic
!            axis relative to east, anisotropy and slope of the orography).
!            At each level a wave Richardson number is computed and by
!            requiring that its value is never less than a critical one,
!            a value of stress can be determined for each model level.
!            The critical Froude number of the flow determines the depth
!            of the low level drag.
!
! Reference: Paper on SSO scheme by Lott and Miller (1997).
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input variables
!     ---------------

!     Scalars
!     -------
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      je        ,    & ! number of grid points in second (meridional) direction
      ke        ,    & ! number of grid points in vertical direction
      ke1              ! ke + 1

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend      ,    & ! end index for first (zonal) direction
      jstart    ,    & ! start index for second (meridional) direction
      jend             ! end index for second (meridional) direction

!     Grid scale variables
!     --------------------
      REAL(KIND=ireals), INTENT(IN) :: pt   (ie,je,ke)
!                            temperature at full levels         (K)
      REAL(KIND=ireals), INTENT(IN) :: pu   (ie,je,ke)
!                            zonal wind component               (m/s)
      REAL(KIND=ireals), INTENT(IN) :: pv   (ie,je,ke)
!                            meridional wind component          (m/s)
      REAL(KIND=ireals), INTENT(IN) :: pfif (ie,je,ke)
!                            geopotential at full levels        (m**2/s**2)
      REAL(KIND=ireals), INTENT(IN) :: pfis (ie,je)
!                            geopotential at surface            (m**2/s**2)
      REAL(KIND=ireals), INTENT(IN) :: pph  (ie,je,ke1)
!                            pressure at half levels            (Pa)
      REAL(KIND=ireals), INTENT(IN) :: ppf  (ie,je,ke)
!                            pressure at full levels            (Pa)

!     GWD-Parameters for each model grid point
!     (standard deviation, anisotropy, angle and slope)
!     -------------------------------------------------
      REAL(KIND=ireals), INTENT(IN) :: psso_stdh  (ie,je)
      REAL(KIND=ireals), INTENT(IN) :: psso_gamma (ie,je)
      REAL(KIND=ireals), INTENT(IN) :: psso_theta (ie,je)
      REAL(KIND=ireals), INTENT(IN) :: psso_sigma (ie,je)

      REAL(KIND=ireals) :: pdt           ! time step

      LOGICAL ldebug ! debug control switch

!     Output variables
!     ----------------

!     Tendencies of T, u and v
!
#ifdef __ICON__
      REAL(KIND=ireals), OPTIONAL :: pdt_sso(ie,je,ke)
#else
      REAL(KIND=ireals) :: pdt_sso(ie,je,ke)
#endif
      REAL(KIND=ireals) :: pdv_sso(ie,je,ke)
      REAL(KIND=ireals) :: pdu_sso(ie,je,ke)

!     Surface (u,v) momentum fluxes and vertically integrated dissipation
!
      REAL(KIND=ireals),INTENT(OUT), OPTIONAL :: pustr_sso(ie,je)
      REAL(KIND=ireals),INTENT(OUT), OPTIONAL :: pvstr_sso(ie,je)
      REAL(KIND=ireals),INTENT(OUT), OPTIONAL :: pvdis_sso(ie,je)

!     Arrays and variables local to *sso* or used for communication
!     with higher level subroutines
!     -----------------------------

!     Indicators of significant levels for the operation of the scheme
      INTEGER  mcrit  (ie,je) 
      INTEGER  mkcrith(ie,je) 
      INTEGER  mkenvh (ie,je) 
      INTEGER  mknu   (ie,je) 
      INTEGER  mknu2  (ie,je) 

      LOGICAL  lo_sso (ie,je)

      REAL(KIND=ireals) :: zfi    (ie,je,ke) 
      ! geopotential minus surface geopotential (m**2/s**2)
      REAL(KIND=ireals) :: ztau   (ie,je,ke1) 
      ! gravity wave stress               (Pa)
      REAL(KIND=ireals) :: zstrdu (ie,je,ke1)
      ! flux of u-momentum (GWD and Blocking)
      REAL(KIND=ireals) :: zstrdv (ie,je,ke1)
      ! flux of v-momentum (GWD and Blocking) 
      REAL(KIND=ireals) :: zstab  (ie,je,ke1)
      ! squared Brunt-vaisala frequency   (1/s**2)
      REAL(KIND=ireals) :: zvph   (ie,je,ke1)
      ! wind profile projected onto plane of low level wind
      REAL(KIND=ireals) :: zrho   (ie,je,ke1)
      ! density at half levels            (kg/m**3)
      REAL(KIND=ireals) :: zri    (ie,je,ke1)
      ! mean flow Richardson number       (-) 
      REAL(KIND=ireals) :: zpsi   (ie,je,ke1)
      ! angle between flow and main axis of orography
      !
      REAL(KIND=ireals) :: zzdep  (ie,je,ke)
      !
      REAL(KIND=ireals) :: zdudt  (ie,je)
      ! sso-tendency for zonal wind    (m/s**2)
      REAL(KIND=ireals) :: zdvdt  (ie,je)
      ! sso-tendency for merid.wind    (m/s**2)
      REAL(KIND=ireals) :: zdtdt  (ie,je)
      ! sso-tendency for temperature   (K/s)
      REAL(KIND=ireals) :: zulow  (ie,je)
      ! u-component of low level wind  (m/s)
      REAL(KIND=ireals) :: zvlow  (ie,je)
      ! v-component of low level wind  (m/s)
      ! directional parameters (see *sso_setup*)
      REAL(KIND=ireals) :: zvidis (ie,je)
      REAL(KIND=ireals) :: zd1    (ie,je) 
      REAL(KIND=ireals) :: zd2    (ie,je) 
      REAL(KIND=ireals) :: zdmod  (ie,je)
!
!     Utility variables
!     -----------------
      REAL(KIND=ireals) :: zgdph,zcons1
#ifndef __ICON__
      REAL(KIND=ireals) :: zdedt
!     REAL(KIND=ireals) :: zdis
#endif
      REAL(KIND=ireals) :: zdelp,ztemp,zb,zc,zcs,zss,zconb,zabsv,zzd1
      REAL(KIND=ireals) :: zratio,zbet,zdt2
!     REAL(KIND=ireals) :: zust,zvst

      INTEGER j1,j2,j3      ! loop indices

!     Timestep is already set for 2TL or 3TL scheme, respectively,
!     in calling routine organize_sso.
!     zdt2  = 2._ireals* pdt
      zdt2  = pdt

!     zcons1=1._ireals/(G*pdt*2._ireals)
      zcons1=1._ireals/(G*zdt2)

!     Initialize tendencies and compute geopotential above ground
!     ===========================================================
      DO j3=1,ke
        DO j2=jstart,jend
          DO j1=istart,iend
            pdu_sso(j1,j2,j3) = 0.0_ireals
            pdv_sso(j1,j2,j3) = 0.0_ireals
#ifndef __ICON__
            pdt_sso(j1,j2,j3) = 0.0_ireals
#endif
            zfi    (j1,j2,j3) = pfif(j1,j2,j3)-pfis(j1,j2)
          END DO
        END DO
      END DO

!     Control operation of scheme by selection of points with standard
!     deviation of sub-grid scale orography > 10 m only
!     =================================================
      DO j2=jstart,jend
        DO j1=istart,iend
          IF (psso_stdh(j1,j2).GT.10._ireals) THEN
            lo_sso(j1,j2)=.TRUE.
          ELSE
            lo_sso(j1,j2)=.FALSE.
          ENDIF
        END DO
      END DO

! ========================================================
!     Computation of basic state variables in *sso_setup*
! ========================================================

      CALL sso_setup (                                   &
         ie    , je    , ke    , ke1     ,               &
         istart, iend  , jstart, jend    ,               &
         pph   , ppf   , pu    , pv      , pt  , zfi  ,  &
         psso_stdh, psso_theta, psso_gamma,              &
         lo_sso,                                         &
         zrho  , zri   , zstab, ztau, zvph, zpsi, zzdep, &
         zulow , zvlow , zd1  , zd2 ,zdmod,              &
         mkcrith, mcrit, mkenvh,mknu,mknu2 )

! ========================================================
!     Surface gravity wave stress amplitude
! ========================================================

      CALL gw_stress (                                   &
         ie    , je    , ke1   ,                         &
         istart, iend  , jstart, jend    ,               &
         zrho,zstab,zvph,psso_stdh,psso_sigma,zdmod,     &
         lo_sso,                                         &
         ztau )

! ========================================================
!     Gravity wave stress profile
! ========================================================

      CALL gw_profil(                                    &
         ie    , je    , ke    , ke1     ,               &
         istart, iend  , jstart, jend    ,               &
         pph   , zrho  , zstab , zvph    , zri  ,        &
         ztau  , zdmod , psso_sigma, psso_stdh  ,        &
         mkcrith, mcrit, mkenvh, mknu    , mknu2,        &
         lo_sso )

! ========================================================
!     Computation of SSO effects' tendencies
! ========================================================

!     Initialisation of tendencies for ALL grid points
!     ------------------------------------------------
      DO j2=jstart,jend
        DO j1=istart,iend
        zvidis (j1,j2)=0.0_ireals
        zdudt  (j1,j2)=0.0_ireals
        zdvdt  (j1,j2)=0.0_ireals
        zdtdt  (j1,j2)=0.0_ireals
        END DO
      END DO
 
!     Compute and add low level drag tendencies to the GWD ones
!     ---------------------------------------------------------
      DO j3=1,ke
       DO j2=jstart,jend
        DO j1=istart,iend

        IF (lo_sso(j1,j2)) THEN

!       Gravity wave drag (cf. documentation EQ.4.13)
!       ---------------------------------------------
        zdelp = pph(j1,j2,j3+1)-pph(j1,j2,j3)
        ztemp = -G*(ztau(j1,j2,j3+1)-ztau(j1,j2,j3))                    &
                  /(zvph(j1,j2,ke1)*zdelp)
        zdudt(j1,j2)=(zulow(j1,j2)*zd1(j1,j2)-zvlow(j1,j2)*zd2(j1,j2))  &
                                    *ztemp/zdmod(j1,j2)
        zdvdt(j1,j2)=(zvlow(j1,j2)*zd1(j1,j2)+zulow(j1,j2)*zd2(j1,j2))  &
                                    *ztemp/zdmod(j1,j2)
        IF (j3 < 4) THEN
         zdudt(j1,j2)= SIGN(MIN(ABS(zdudt(j1,j2)),20._ireals/3600._ireals),zdudt(j1,j2))
         zdvdt(j1,j2)= SIGN(MIN(ABS(zdvdt(j1,j2)),20._ireals/3600._ireals),zdvdt(j1,j2))
        ENDIF

!       Low level drag ('blocking') (cf. documentation EQ.4.14 ff.)
!       -----------------------------------------------------------
        IF (j3.GE.mkenvh(j1,j2)) THEN
         zb  = 1.0_ireals-0.18_ireals*psso_gamma(j1,j2)-0.04_ireals*psso_gamma(j1,j2)**2
         zc  = 0.48_ireals*psso_gamma(j1,j2)+0.3_ireals*psso_gamma(j1,j2)**2
         zcs = COS(zpsi(j1,j2,j3))**2
         zss = 1.0_ireals-zcs
         zzd1  =zb*zcs+zc*zss
         zconb =zdt2*Gkwake*psso_sigma(j1,j2)/(2._ireals*psso_stdh(j1,j2))
         zabsv =0.5_ireals*SQRT(pu(j1,j2,j3)**2+pv(j1,j2,j3)**2)
         zratio=(zcs+psso_gamma(j1,j2)*zss)/(psso_gamma(j1,j2)*zcs+zss)
         zbet  =MAX(0._ireals,2._ireals-1._ireals/zratio)*zconb*zzdep(j1,j2,j3)*zzd1*zabsv
!        Partially implicit tendency calculation
!        ---------------------------------------
         zdudt(j1,j2)=-pu(j1,j2,j3)/zdt2*(zbet/(1._ireals+zbet)) 
         zdvdt(j1,j2)=-pv(j1,j2,j3)/zdt2*(zbet/(1._ireals+zbet)) 
        END IF

        pdu_sso(j1,j2,j3)=zdudt(j1,j2)
        pdv_sso(j1,j2,j3)=zdvdt(j1,j2)
! GZ: use differential formulation in order to avoid artifacts (and numerical
!     instabilities) due to possible sign changes at long time steps
!        zust             = pu(j1,j2,j3)+zdt2*zdudt(j1,j2)
!        zvst             = pv(j1,j2,j3)+zdt2*zdvdt(j1,j2)
!        zdis=0.5_ireals*(pu(j1,j2,j3)**2+pv(j1,j2,j3)**2-zust**2-zvst**2)
!        zdedt            = zdis/zdt2
#ifndef __ICON__
        zdedt = -(pu(j1,j2,j3)*zdudt(j1,j2)+pv(j1,j2,j3)*zdvdt(j1,j2))
        zdtdt(j1,j2)     = zdedt       /Cp_d
        pdt_sso(j1,j2,j3)= zdtdt(j1,j2)
!        zvidis(j1,j2)=zvidis(j1,j2)+zdis*zdelp    ! de-activated
        zvidis(j1,j2)=zvidis(j1,j2)+zdedt*zdelp    ! de-activated
#endif
        ENDIF

        END DO
       END DO
      END DO     ! loop over vertical layers

! ======================================================================
!     Flux computations of original code of *GWDRAG* are not required in
!     COSMO model, but may be reactivated by uncommenting the following
!     statements and declaration of the appropriate arrays
! ======================================================================

!     Stress components and dissipation
!     ---------------------------------

      IF(PRESENT(pvdis_sso))THEN
        DO j2=jstart,jend
          DO j1=istart,iend
            IF(lo_sso(j1,j2)) THEN
!              pvdis_sso(j1,j2)=zcons1*zvidis(j1,j2)
              pvdis_sso(j1,j2)=zvidis(j1,j2)/G
            ELSE
              pvdis_sso(j1,j2)=0._ireals
            ENDIF
          END DO
        END DO
      ENDIF

!     Initialize flux at top
!     ----------------------
      DO j2=jstart,jend
        DO j1=istart,iend
          zstrdu(j1,j2,1)=0._ireals
          zstrdv(j1,j2,1)=0._ireals
        END DO
      END DO

!     Increment flux based on tendency in each layer
!     ----------------------------------------------
      DO j3=1,ke
        DO j2=jstart,jend
          DO j1=istart,iend
            zgdph=-G  /(pph(j1,j2,j3+1)-pph(j1,j2,j3))
            zstrdu(j1,j2,j3+1)=pdu_sso(j1,j2,j3)/zgdph + zstrdu(j1,j2,j3)
            zstrdv(j1,j2,j3+1)=pdv_sso(j1,j2,j3)/zgdph + zstrdv(j1,j2,j3)
          END DO
        END DO
      END DO

      IF(PRESENT(pustr_sso))THEN
!     Store flux at surface
!     ---------------------
        DO j2=jstart,jend
          DO j1=istart,iend
            pustr_sso(j1,j2)=zstrdu(j1,j2,ke1)
            pvstr_sso(j1,j2)=zstrdv(j1,j2,ke1)
          END DO
        END DO
      ENDIF
!
!     Control printout
!     ----------------
      IF (ldebug) THEN
        DO j2=jstart,jend
          DO j1=istart,iend
          IF (j1.EQ.55 .and. j2.EQ.64) THEN
            PRINT *, ' '
            PRINT *, ' '
            PRINT *, ' Diagnosis SSO scheme j1=55 j2=64'
            PRINT *, ' '
            PRINT *, ' fis      : ', pfis       (j1,j2)
            PRINT *, ' sso_stdh : ', psso_stdh  (j1,j2)
            PRINT *, ' sso_gamma: ', psso_gamma (j1,j2)
            PRINT *, ' sso_theta: ', psso_theta (j1,j2)
            PRINT *, ' sso_sigma: ', psso_sigma (j1,j2)
            PRINT *, ' '
            PRINT *, '  j3    ph (Pa)     pf (Pa)     u(m/s)    ',  &
           'v (m/s)     T (K)      fif (m^2/s^2)     du_sso     ',  &
           'dv_sss     dt_sso'
          ENDIF
          ENDDO
        ENDDO
  
        DO j2=jstart,jend
          DO j1=istart,iend
          IF (j1.EQ.55 .and. j2.EQ.64) THEN
            DO j3=1,ke
            WRITE (*,'(i3, 9E13.6)') j3, pph(j1,j2,j3+1), ppf(j1,j2,j3),  &
            pu(j1,j2,j3), pv(j1,j2,j3), pt(j1,j2,j3), pfif(j1,j2,j3),     &
            pdu_sso(j1,j2,j3), pdv_sso(j1,j2,j3), pdt_sso(j1,j2,j3)
            ENDDO
          ENDIF
          ENDDO
        ENDDO
      ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sso

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE sso_setup (                                      &
           ie     , je     , ke     , ke1      ,            &
           istart , iend   , jstart , jend     ,            &
           pph ,ppf ,pu ,pv ,pt ,pfi,                       &
           psso_stdh, psso_theta, psso_gamma   ,            &
           lo_sso,                                          &
           prho  , pri   , pstab, ptau, pvph , ppsi, pzdep, &
           pulow , pvlow , pd1  , pd2 , pdmod,              &
           kkcrith, kcrit, kkenvh,kknu,kknu2)

!------------------------------------------------------------------------------
!
! Purpose: The module procedure sso_setup sets up parameters for the
!          parameterization of sub-grid scale orographic effects:
!
!          - definition of various reference model levels
!          - Brunt-Vaisala frequency on half levels
!          - mean wind components in the layer between one and two
!            standard deviations of sso above ground
!          - geometry factors, relating the orientation of sso and wind
!          - Phillips parameters
!          - vertical wind profile in plane of gravity wave stress
!          - basic flow Richardson number
!          - calculation of depth of 'blocked' layer
!          - calculation of layer in which low level wave braking occurs
!          - calculation of assumed vertical profile of sso
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input arrays and variables
!     ==========================
!
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      je        ,    & ! number of grid points in second (meridional) direction
      ke        ,    & ! number of grid points in vertical direction
      ke1              ! ke + 1

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend      ,    & ! end index for first (zonal) direction
      jstart    ,    & ! start index for second (meridional) direction
      jend             ! end index for second (meridional) direction

      REAL(KIND=ireals) :: pph (ie,je,ke1)
      REAL(KIND=ireals) :: ppf (ie,je,ke)
!
      REAL(KIND=ireals) :: pu  (ie,je,ke)
      REAL(KIND=ireals) :: pv  (ie,je,ke)
      REAL(KIND=ireals) :: pt  (ie,je,ke)
      REAL(KIND=ireals) :: pfi (ie,je,ke)

!     subgrid scale orography parameters
      REAL(KIND=ireals) :: psso_stdh (ie,je)
      REAL(KIND=ireals) :: psso_theta(ie,je) 
      REAL(KIND=ireals) :: psso_gamma(ie,je) 

      LOGICAL lo_sso(ie,je)

!     Output arrays
!     =============

      REAL(KIND=ireals) :: prho (ie,je,ke1)
!     density on half levels          (kg/m**3)
      REAL(KIND=ireals) :: pri  (ie,je,ke1)
!     mean flow Richardson number     (-)
      REAL(KIND=ireals) :: pstab(ie,je,ke1)
!     squared Brunt-Vaisala frequency (1/s**2)
      REAL(KIND=ireals) :: ptau (ie,je,ke1)
!     gravity wave stress profile     (Pa)
      REAL(KIND=ireals) :: pvph (ie,je,ke1) 
!     projected flow on half levels   (m/s)
      REAL(KIND=ireals) :: ppsi (ie,je,ke1)
!     angle between orography and blocked flow (1:ke)
!                            or low level flow (ke1)
      REAL(KIND=ireals) :: pzdep(ie,je,ke)
!     height dependency factor for 'blocking' tendency
      REAL(KIND=ireals) :: pulow(ie,je) 
!     low level zonal wind            (m/s)
      REAL(KIND=ireals) :: pvlow(ie,je)
!     low level meridional wind       (m/s)

!     directional parameters
      REAL(KIND=ireals) :: pd1  (ie,je)
      REAL(KIND=ireals) :: pd2  (ie,je)
      REAL(KIND=ireals) :: pdmod(ie,je)

      INTEGER kkcrith(ie,je)
!     maximum level for wave breaking

      INTEGER kcrit  (ie,je) 
!     critical level

      INTEGER kkenvh (ie,je) 
!     index of top of 'envelope' / blocking level

      INTEGER kknu   (ie,je)
!     model level at 4*stdh

      INTEGER kknu2  (ie,je)
!     model level at 3*stdh

!     local arrays and variables
!     ==========================

      REAL(KIND=ireals) :: zvpf   (ie,je,ke)
      ! projected flow on full levels (m/s)
      REAL(KIND=ireals) :: zdp    (ie,je,ke)
      ! pressure difference between layers
      REAL(KIND=ireals) :: zsqst  (ie,je,ke)
      REAL(KIND=ireals) :: znorm  (ie,je)
      REAL(KIND=ireals) :: znup   (ie,je)
      REAL(KIND=ireals) :: znum   (ie,je)
      REAL(KIND=ireals) :: znu    (ie,je) 

      REAL(KIND=ireals) :: zcons1, zcons2 ! utility constants
      REAL(KIND=ireals) :: zu             ! security for low level zonal wind
      REAL(KIND=ireals) :: zb             ! Phillips parameter B
      REAL(KIND=ireals) :: zc             ! Phillips parameter C
      REAL(KIND=ireals) :: zdelp          ! pressure thickness of layers
      REAL(KIND=ireals) :: zvt1,zvt2      ! utility variables for flow projection
      REAL(KIND=ireals) :: zst            ! utility variable for stability calculation
      REAL(KIND=ireals) :: zdwind         ! utility variable for wind shear calculation
      REAL(KIND=ireals) :: zwind          ! utility variable for proj. wind calculation
      REAL(KIND=ireals) :: zggeenv,zggeo,zgvar ! geopotential utility variables
      REAL(KIND=ireals) :: zhcrit 

      INTEGER mknub(ie,je)
      INTEGER mknul(ie,je)

      INTEGER mi3h              ! vertical loop limit
      INTEGER j1,j2,j3          ! loop variables   

      LOGICAL lo1  (ie,je,ke1)
      LOGICAL llo               ! utility switch

!     The following parameter is a tunable constant for the sub-grid scale
!     orography scheme
!     ================

      INTEGER (KIND=iintegers) :: Nktopg                    
                                ! number of topmost layer used to define low level
                                ! flow in case of high vertical model resolution

!-------------------------------------------------------------------------------------

! Begin subroutine

      Nktopg = ke               ! number of topmost layer used to defined low level

!     computational constants
!     =======================

!     mi3h =(ki3e-ki3s+1)/3
      mi3h =ke/3

      zcons1=1._ireals/R_d
      zcons2=G**2/Cp_d
 
!C*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
!C*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
!C*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
!C
!     security on anisotropy factor and presetting of critical levels

      DO j2=jstart,jend
        DO j1=istart,iend
        psso_gamma(j1,j2) = MAX(psso_gamma(j1,j2),Gtsec)
        kknu      (j1,j2) = ke
        kknu2     (j1,j2) = ke
        mknub     (j1,j2) = ke
        mknul     (j1,j2) = ke
        lo1(j1,j2,ke1)    =.FALSE.    ! Initialize control variable
        END DO
      END DO

!
!!!!  define top of low level drag calculation layer (*kkcrit*)
!     and other critical levels
!     ============================================================
!
      DO j3=ke,mi3h,-1     ! vertical loop   
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =4._ireals*psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            kknu(j1,j2)=j3  ! first layer with height > 4*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

      DO j3=ke,mi3h,-1    ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =3._ireals*psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            kknu2(j1,j2)=j3 ! first layer with height > 3*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

      DO j3=ke,mi3h,-1    ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =2._ireals*psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            mknub(j1,j2)=j3  ! first layer with height > 2*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

      DO j3=ke,mi3h,-1    ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            mknul(j1,j2)=j3 ! first layer with height > 1*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Confine critical level indices to be less or equal to Nktopg
!     ============================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        kknu(j1,j2) =MIN(kknu(j1,j2),Nktopg)
        mknub(j1,j2)=MIN(mknub(j1,j2),Nktopg)
        IF(mknub(j1,j2).EQ.Nktopg) mknul(j1,j2)=ke
        IF(mknub(j1,j2).EQ.mknul(j1,j2)) mknub(j1,j2) = mknub(j1,j2) - 1
        END DO
      END DO

!     Initialize various arrays
!     =========================
      DO j2=jstart,jend
        DO j1=istart,iend
        prho (j1,j2,ke1) = 0.0_ireals
        pstab(j1,j2,1  ) = 0.0_ireals
        pstab(j1,j2,ke1) = 0.0_ireals
        pri  (j1,j2,1  ) = 0.0_ireals
        pri  (j1,j2,ke1) = 9999.0_ireals
        ppsi (j1,j2,ke1) = 0.0_ireals
        pvph (j1,j2,1)   = 0.0_ireals
        pulow(j1,j2)     = 0.0_ireals
        pvlow(j1,j2)     = 0.0_ireals
        kkcrith(j1,j2)   = ke
        kkenvh(j1,j2)    = ke ! default for top of envelope layer
        kcrit(j1,j2)     = 1  ! default for critical level
        znu  (j1,j2)     = 0.0_ireals
        znum (j1,j2)     = 0.0_ireals
        lo1  (j1,j2,ke1) = .FALSE.
        END DO
      END DO
 
!     pressure thickness, density and Brunt-Vaisala frequency (squared)
!     =================================================================
      DO j3=ke,2,-1        ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          zdp (j1,j2,j3) = ppf(j1,j2,j3)-ppf(j1,j2,j3-1)
!         density on half levels
          prho (j1,j2,j3) = 2._ireals*pph(j1,j2,j3)*zcons1                 &
     &                     /(pt(j1,j2,j3)+pt(j1,j2,j3-1))
!         squared Brunt-Vaisala frequency on half levels
          pstab(j1,j2,j3)= 2._ireals*zcons2/(pt(j1,j2,j3)+pt(j1,j2,j3-1))  &
     &                    *( 1._ireals-Cp_d*prho(j1,j2,j3)                 &
     &                         *(pt(j1,j2,j3)-pt(j1,j2,j3-1))              &
     &                             / zdp(j1,j2,j3) )
!         security on Brunt-Vaisala frequency
          pstab(j1,j2,j3)=MAX(pstab(j1,j2,j3),gssec)
          zsqst(j1,j2,j3)=SQRT(pstab(j1,j2,j3))
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Definition of blocked flow
!     ==========================
      DO j3=ke,mi3h,-1          ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.mknub(j1,j2).AND.j3.LE.mknul(j1,j2)) THEN
!         only layers with height between one and two stdh contribute 
          pulow(j1,j2) = pulow(j1,j2)                                     &
     &                  +pu   (j1,j2,j3)*(pph(j1,j2,j3+1)-pph(j1,j2,j3))
          pvlow(j1,j2) = pvlow(j1,j2)                                     &
     &                  +pv   (j1,j2,j3)*(pph(j1,j2,j3+1)-pph(j1,j2,j3))
          END IF
          END IF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Division by pressure thickness of contributing layers and
!     determination of wind speed of blocked flow
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        pulow(j1,j2) = pulow(j1,j2) /                                    &
     &              (pph(j1,j2,mknul(j1,j2)+1)-pph(j1,j2,mknub(j1,j2)))
        pvlow(j1,j2) = pvlow(j1,j2) /                                    &
     &              (pph(j1,j2,mknul(j1,j2)+1)-pph(j1,j2,mknub(j1,j2)))
        znorm(j1,j2)=MAX(SQRT(pulow(j1,j2)**2+pvlow(j1,j2)**2),Gvsec)
        pvph(j1,j2,ke1)=znorm(j1,j2)  ! projected flow at lower boundary
        END IF
        END DO
      END DO
 
!     Axes of subgrid scale orography and plane of profiles 
!     =====================================================
      DO j2=jstart,jend      
        DO j1=istart,iend   
        IF (lo_sso(j1,j2)) THEN
        llo=(pulow(j1,j2).LT.Gvsec).AND.(pulow(j1,j2).GE.-Gvsec)
          IF(llo) THEN
          ZU=pulow(j1,j2)+2._ireals*Gvsec
          ELSE
          ZU=pulow(j1,j2)
          ENDIF
!       angle between principal axis of orog. and low-level wind direction
!       ------------------------------------------------------------------
        ppsi (j1,j2,ke1) = psso_theta(j1,j2)-ATAN(pvlow(j1,j2)/ZU)
!       Phillips parameters B and C
!       ---------------------------
        zb           = 1._ireals-0.18_ireals*psso_gamma(j1,j2)                      &
     &                   -0.04_ireals*psso_gamma(j1,j2)**2
        zc           = 0.48_ireals*psso_gamma(j1,j2)+0.3_ireals*psso_gamma(j1,j2)**2
!       projection parameters D1 and D2 (see documentation)
        pd1  (j1,j2) = zb-(zb-zc)*(SIN(ppsi(j1,j2,ke1))**2)
        pd2  (j1,j2) = (zb-zc)*SIN(ppsi(j1,j2,ke1))                   &
     &                        *COS(ppsi(j1,j2,ke1))
        pdmod(j1,j2) = SQRT(pd1(j1,j2)**2+pd2(j1,j2)**2)
        END IF
        END DO
      END DO
 
!     projection of flow into plane of low level stress  (eq.4.7)
!     ===========================================================
      DO j3=1,ke                  ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          zvt1 = pulow(j1,j2)*pu(j1,j2,j3)+pvlow(j1,j2)*pv(j1,j2,j3)
          zvt2 =-pvlow(j1,j2)*pu(j1,j2,j3)+pulow(j1,j2)*pv(j1,j2,j3)
          zvpf(j1,j2,j3)=(zvt1*pd1(j1,j2)+zvt2*pd2(j1,j2))             &
     &                      /(znorm(j1,j2)*pdmod(j1,j2))
          ENDIF

!      initialize stress array *ptau*, depth of blocked layer *pzdep*
!      and angle array *ppsi* (not for j3=ke1) and reset control
!      variable *llo1*
!      -----------------------------------------------------------------
          ptau(j1,j2,j3)  =0.0_ireals
          pzdep(j1,j2,j3) =0.0_ireals
          ppsi(j1,j2,j3)  =0.0_ireals
          lo1(j1,j2,j3)   =.FALSE.
          END DO
        END DO
      END DO                ! end of vertical loop
!C
!     linear interpolation of projected flow to half levels
!     ========================================================
      DO j3=2,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          pvph(j1,j2,j3)=                                        &
     &     ((pph(j1,j2,j3)-ppf(j1,j2,j3-1))*zvpf(j1,j2,j3) +     &
     &      (ppf(j1,j2,j3)-pph(j1,j2,j3  ))*zvpf(j1,j2,j3-1))    &
     &                         /zdp(j1,j2,j3)
            IF(pvph(j1,j2,j3).LT.Gvsec) THEN  ! critical layer
            pvph(j1,j2,j3)=Gvsec
            if (j3.lt.mknub(j1,j2)) kcrit(j1,j2  )=j3
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Brunt-Vaisala frequency and density for lowest level
!     ====================================================
      DO j3=mi3h,ke   ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            IF(j3.GE.(mknub(j1,j2)+1).AND.j3.LE.mknul(j1,j2)) THEN
            zst=zcons2/pt(j1,j2,j3)*(1._ireals-Cp_d*prho(j1,j2,j3)*     &
     &                    (pt(j1,j2,j3)-pt(j1,j2,j3-1))/zdp(j1,j2,j3))
            pstab(j1,j2,ke1)=pstab(j1,j2,ke1)+zst*zdp(j1,j2,j3)
            pstab(j1,j2,ke1)=MAX(pstab(j1,j2,ke1),Gssec)
            prho (j1,j2,ke1)= prho(j1,j2,ke1)                     &
     &                    +pph(j1,j2,j3)*2._ireals*zdp(j1,j2,j3)  &
     &                    *zcons1/(pt(j1,j2,j3)+pt(j1,j2,j3-1))
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     normalization
!     -------------
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        pstab(j1,j2,ke1)=pstab(j1,j2,ke1)                    &
     &      /(ppf(j1,j2,mknul(j1,j2))-ppf(j1,j2,mknub(j1,j2)))
        pstab(j1,j2,ke1)=MAX(pstab(j1,j2,ke1),gssec)
        prho (j1,j2,ke1)=prho(j1,j2,ke1)                     &
     &      /(ppf(j1,j2,mknul(j1,j2))-ppf(j1,j2,mknub(j1,j2)))
        END IF
        END DO
      END DO
 
!     mean flow Richardson number on half levels
!     ==========================================
      DO j3=2,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          zdwind=MAX(ABS(zvpf(j1,j2,j3)-zvpf(j1,j2,j3-1)),Gvsec)
          pri(j1,j2,j3)=pstab(j1,j2,j3)*(zdp(j1,j2,j3)                  &
     &            /(G*prho(j1,j2,j3)*zdwind))**2
          pri(j1,j2,j3)=MAX(pri(j1,j2,j3),Grcrit)
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

!     define top of 'envelope' layer (cf. eq.4.8)
!     ===========================================
      DO j3=2,ke-1     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            IF (j3.GE.kknu2(j1,j2)) THEN
            znum (j1,j2)=znu(j1,j2)
            zwind= (pulow(j1,j2)*pu(j1,j2,j3)                           &
     &             +pvlow(j1,j2)*pv(j1,j2,j3))/                         &
     &         MAX(SQRT(pulow(j1,j2)**2+pvlow(j1,j2)**2),Gvsec)
            zwind=MAX(ABS(zwind),Gvsec)
            zdelp=pph(j1,j2,j3+1)-pph(j1,j2,j3)
!           vertically integrated left side of eq.4.8
            znu(j1,j2) = znu(j1,j2) + (zdelp/G)*                             &
     &              ( (zsqst(j1,j2,j3+1)/prho(j1,j2,j3+1)                    &
     &                +zsqst(j1,j2,j3  )/prho(j1,j2,j3  ) )/2._ireals)/zwind
            IF((znum(j1,j2).LE.gfrcrit).AND.(znu(j1,j2).GT.gfrcrit)          &
     &                          .AND.(kkenvh(j1,j2).EQ.ke))                &
     &      kkenvh(j1,j2)=j3     
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

!     dynamical mixing height for the breaking of gravity waves
!     =========================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        znup(j1,j2)=0.0_ireals
        znum(j1,j2)=0.0_ireals
        END DO
      END DO

      DO j3=ke-1,2,-1  ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            IF (j3.LT.kkenvh(j1,j2)) THEN    ! only above envelope height
            znum(j1,j2)=znup(j1,j2)
            zwind=(pulow(j1,j2)*pu(j1,j2,j3)+pvlow(j1,j2)*pv(j1,j2,j3))/    &
     &            MAX(SQRT(pulow(j1,j2)**2+pvlow(j1,j2)**2),Gvsec)
            zwind=MAX(ABS(zwind),Gvsec)
            zdelp=pph(j1,j2,j3+1)-pph(j1,j2,j3)
            znup(j1,j2) = znup(j1,j2) + (zdelp/G)*                          &
     &            ( (zsqst(j1,j2,j3+1)/prho(j1,j2,j3+1)                     &
     &              +zsqst(j1,j2,j3  )/prho(j1,j2,j3  ) )/2._ireals)/zwind
            IF((znum(j1,j2).LE.1.5_ireals).AND.(znup(j1,j2).GT.1.5_ireals)                &
     &                          .AND.(kkcrith(j1,j2).EQ.ke))              &
     &      kkcrith(j1,j2)=j3
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

!     allow low level wave breaking only above height of 4*stdh 
      DO j2=jstart,jend
        DO j1=istart,iend
        kkcrith(j1,j2)=MIN(kkcrith(j1,j2),kknu(j1,j2))
        END DO
      END DO
 
!     directional information for flow blocking
!     =========================================
      DO j3=mi3h,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.kkenvh(j1,j2)) THEN  ! only within envelope layer
          llo=(pu(j1,j2,j3).LT.Gvsec).AND.(pu(j1,j2,j3).GE.-Gvsec)
            IF(llo) THEN
            ZU=pu(j1,j2,j3)+2._ireals*Gvsec
            ELSE
            ZU=pu(j1,j2,j3)
            ENDIF
          ppsi(j1,j2,j3)=psso_theta(j1,j2)-ATAN(pv(j1,j2,j3)/ZU)
          ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     assumed vertical profile of sso for blocking calculations
!     =========================================================
      DO j3=mi3h,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.kkenvh(j1,j2)) THEN   ! only within envelope layer
          zggeenv= MAX(1._ireals,                                             &
     &     (pfi(j1,j2,kkenvh(j1,j2))+pfi(j1,j2,kkenvh(j1,j2)-1))/2._ireals)
          zggeo  = MAX(pfi(j1,j2,j3),1._ireals)
          zgvar  = MAX(psso_stdh(j1,j2)*G,1._ireals)
          pzdep(j1,j2,j3)=SQRT((zggeenv-zggeo)/(zggeo+zgvar))
          END IF
          END IF
          END DO
        END DO
      END DO                ! end of vertical loop

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sso_setup

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE gw_stress (                                  &
           ie     , je     , ke1    ,                   &
           istart , iend   , jstart , jend     ,        &
           prho,pstab,pvph,psso_stdh,psso_sigma,pdmod,  &
           lo_sso, ptau )

!------------------------------------------------------------------------------
!
! Purpose: The module procedure gw_stress computes the gravity stress
!          amplitude following eq.4.11 of the documentation.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input
!     =====
!
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      je        ,    & ! number of grid points in second (meridional) direction
      ke1              ! ke + 1

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend      ,    & ! end index for first (zonal) direction
      jstart    ,    & ! start index for second (meridional) direction
      jend             ! end index for second (meridional) direction

      REAL(KIND=ireals) :: prho (ie,je,ke1)
      ! density on half levels    (kg/m**3)
      REAL(KIND=ireals) :: pstab(ie,je,ke1)
      ! squared Brunt-Vaisala frequency  (1/s**2)
      REAL(KIND=ireals) :: pvph (ie,je,ke1)
      ! wind on half levels projected on plane of low level wind (m/s)
      REAL(KIND=ireals) :: psso_stdh (ie,je)
      ! standard deviation of sso-height (m)
      REAL(KIND=ireals) :: psso_sigma(ie,je)
      ! mean sso-slope                   (-)
      REAL(KIND=ireals) :: pdmod(ie,je)
      ! projection parameter = SQRT(D1**2+D2**2)    cf. eq.4.7

      LOGICAL lo_sso(ie,je)
      ! 

!     Output
!     ======

      REAL(KIND=ireals) :: ptau(ie,je,ke1)
      ! gravity wave stress amplitude   (Pa)

!     local variables
!     ===============
                                  ! utility variables, which may be used to modify
      REAL(KIND=ireals) :: zblock ! the magnitude of the subgrid scale standard
      REAL(KIND=ireals) :: zeff   ! deviation which enters the stress amplitude
                                  ! calculation (zblock=0.0 in operational version)

      INTEGER j1,j2  ! loop variables

!     gravity wave stress amplitude (eq.4.11)
!     =======================================

      DO j2=jstart,jend
        DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            zblock=0.0_ireals 
            zeff=MAX(0._ireals,2._ireals*psso_stdh(j1,j2)-zblock)
            ptau(j1,j2,ke1)=Gkdrag*prho(j1,j2,ke1)*psso_sigma(j1,j2)  &
     &                     *zeff**2/4._ireals                         &
     &                     /psso_stdh(j1,j2)*pvph(j1,j2,ke1)          &
     &                     *pdmod(j1,j2)*SQRT(pstab(j1,j2,ke1))
          ELSE
            ptau(j1,j2,ke1)=0.0_ireals
          ENDIF
        END DO
      END DO

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gw_stress

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE gw_profil(                                    &
           ie     , je     , ke     , ke1      ,         &
           istart , iend   , jstart , jend     ,         &
           pph    , prho   , pstab  , pvph     , pri ,   &
           ptau   , pdmod  , psso_sigma, psso_stdh   ,   &
           kkcrith, kcrit, kkenvh, kknu, kknu2 ,         &
           lo_sso )

!------------------------------------------------------------------------------
!
! Purpose: The module procedure gw_profil computes the vertical profile of
!          gravity wave stress.
!
! Method:  The stress profile for gravity waves is computed as follows:
!
!          - it is constant (no gravity wave drag) at all levels
!            between the ground and the top of the blocked layer kkenvh
!          - it decreases linearly with height from the top of the
!            blocked layer to 3 * stdh (kknu), to simulate lee waves or
!            non-linear gravity wave breaking
!          - at levels above (kknu) it is constant, except when the
!            wave encounters a critical level (kcrit) or when it breaks
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input
!     =====
!
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      je        ,    & ! number of grid points in second (meridional) direction
      ke        ,    & ! number of grid points in vertical direction
      ke1              ! ke + 1

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend      ,    & ! end index for first (zonal) direction
      jstart    ,    & ! start index for second (meridional) direction
      jend             ! end index for second (meridional) direction

      REAL(KIND=ireals) :: pph  (ie,je,ke1)
      ! half level pressure           (Pa)
      REAL(KIND=ireals) :: prho (ie,je,ke1)
      ! density on half levels    (kg/m**3)
      REAL(KIND=ireals) :: pstab(ie,je,ke1)
      ! squared Brunt-Vaisala frequency  (1/s**2)
      REAL(KIND=ireals) :: pri  (ie,je,ke1)
      ! mean flow Richardson number      ( - )    
      REAL(KIND=ireals) :: pvph (ie,je,ke1)
      ! wind on half levels projected on plane of low level wind (m/s)
      REAL(KIND=ireals) :: ptau (ie,je,ke1)
      !gravity wave stress profile
      REAL(KIND=ireals) :: pdmod(ie,je)
      ! projection parameter = SQRT(D1**2+D2**2)    cf. eq.4.7
      REAL(KIND=ireals) :: psso_stdh (ie,je)
      ! standard deviation of sso-height (m)
      REAL(KIND=ireals) :: psso_sigma(ie,je)
      ! mean sso-slope                   (-)

!     various significant levels
      INTEGER kkcrith(ie,je)
      INTEGER kcrit  (ie,je)
      INTEGER kkenvh (ie,je)
      INTEGER kknu   (ie,je)
      INTEGER kknu2  (ie,je)

      LOGICAL lo_sso (ie,je)
     
!     local arrays and variables
!     ==========================

      REAL(KIND=ireals) :: zdz2   (ie,je,ke)
      REAL(KIND=ireals) :: ztau   (ie,je,ke1)
      REAL(KIND=ireals) :: znorm  (ie,je)
      REAL(KIND=ireals) :: zoro   (ie,je)

      REAL(KIND=ireals) :: zb,zdelp,zdelpt,zalpha,zalfa,zdel,zriw    ! utitility variables
      REAL(KIND=ireals) :: zsqri,zdz2n                               ! utitility variables

      INTEGER j1,j2,j3                 ! loop indices

      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        zoro(j1,j2) = psso_sigma(j1,j2)*pdmod(j1,j2)  &
     &               /(4._ireals*MAX(psso_stdh(j1,j2),1.0_ireals))
!HF  &               /4._ireals/MAX(psso_stdh(j1,j2),1.0_ireals)
        ztau(j1,j2,kknu(j1,j2)+1) = ptau(j1,j2,kknu(j1,j2)+1)
        ztau(j1,j2,ke1          ) = ptau(j1,j2,ke1          )
        ENDIF   
        END DO
      END DO

      DO j3=ke,2,-1     ! vertical loop

!     constant stress up to top of blocking layer
!     ===========================================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.kknu2(j1,j2)) THEN
          ptau(j1,j2,j3)=ztau(j1,j2,ke1)
          ENDIF
        ENDIF   
        END DO
      END DO

!     wave displacement at next level
!     ===============================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
          IF(j3.LT.kknu2(j1,j2)) THEN
          znorm(j1,j2)=Gkdrag*prho(j1,j2,j3)*SQRT(pstab(j1,j2,j3))  &
     &              *pvph(j1,j2,j3)*zoro(j1,j2)
          zdz2(j1,j2,j3)=ptau(j1,j2,j3+1)/MAX(znorm(j1,j2),Gssec)
          ENDIF
        ENDIF   
        END DO
      END DO

!     wave Richardson number, new wave displacement and stress
!     breaking evaluation and critical level
!     ========================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
          IF(j3.LT.kknu2(j1,j2)) THEN    ! only above blocking layer
            IF((ptau(j1,j2,j3+1).LT.Gtsec).OR.(j3.LE.kcrit(j1,j2))) THEN
            ptau(j1,j2,j3)=0.0_ireals    ! above critical levels
            ELSE
            zsqri=SQRT(pri(j1,j2,j3))
            zalfa=SQRT(pstab(j1,j2,j3)*zdz2(j1,j2,j3))/pvph(j1,j2,j3)
!HF         zriw=pri(j1,j2,j3)*(1._ireals-zalfa)/(1+zalfa*zsqri)**2
            zriw=pri(j1,j2,j3)*(1._ireals-zalfa)/(1._ireals+zalfa*zsqri)**2
              IF(zriw.LT.Grcrit) THEN      ! breaking occurs
              zdel=4._ireals/zsqri/Grcrit+1._ireals/Grcrit**2+4._ireals/Grcrit
              zb=1._ireals/Grcrit+2._ireals/zsqri
              zalpha=0.5_ireals*(-zb+SQRT(zdel))
              zdz2n=(pvph(j1,j2,j3)*zalpha)**2/pstab(j1,j2,j3)
              ptau(j1,j2,j3)=znorm(j1,j2)*zdz2n
              ELSE
              ptau(j1,j2,j3)=znorm(j1,j2)*zdz2(j1,j2,j3)
              ENDIF
            ptau(j1,j2,j3)=MIN(ptau(j1,j2,j3),ptau(j1,j2,j3+1))
            ENDIF
          ENDIF
        ENDIF    
        END DO
      END DO

      END DO       ! end of vertical loop

!     reorganisation of stress profile, if breaking occurs at low levels
!     ==================================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        ztau(j1,j2,kkenvh(j1,j2)) =ptau(j1,j2,kkenvh(j1,j2))
        ztau(j1,j2,kkcrith(j1,j2))=ptau(j1,j2,kkcrith(j1,j2))
        ENDIF   
        END DO
      END DO

!     linear decrease between kkenvh and kkcrith
      DO j3=1,ke      ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GT.kkcrith(j1,j2).AND.j3.LT.kkenvh(j1,j2))THEN
          zdelp=pph(j1,j2,j3)-pph(j1,j2,kkenvh(j1,j2))
          zdelpt=pph(j1,j2,kkcrith(j1,j2))-pph(j1,j2,kkenvh(j1,j2))
          ptau(j1,j2,j3)=ztau(j1,j2,kkenvh(j1,j2)) +                    &
     &        (ztau(j1,j2,kkcrith(j1,j2))-ztau(j1,j2,kkenvh(j1,j2)) )*  &
     &            zdelp/zdelpt
          ENDIF
          ENDIF            
          END DO
        END DO
      END DO       ! end of vertical loop

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gw_profil

!------------------------------------------------------------------------------
! End of module mo_sso_cosmo
!------------------------------------------------------------------------------

END MODULE mo_sso_cosmo
