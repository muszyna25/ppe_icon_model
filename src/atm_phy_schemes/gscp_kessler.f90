!>
!! cloud microphysics
!!
!! !------------------------------------------------------------------------------
!!
!! @par Description of *gscp_kessler*:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, water vapor and rain due to cloud microphysical processes related 
!!   to the formation of grid scale precipitation. This includes the sedimentation 
!!   of rain and snow. 
!!   The precipitation fluxes at the surface are also calculated here.
!!
!! Method:
!!   Prognostic one-moment bulk microphysical parameterization.
!!   The sedimentation of rain and snow is computed implicitly.
!!
!! Current Code Owner: DWD, A. Seifert
!!    phone: +49-69-8062-2729,  fax:   +49-69-8062-3721
!!    email: Axel.Seifert@dwd.de
!!
!! @author Guenther Doms
!! @author Axel Seifert
!!
!! @par reference   This is an adaption of subroutine kessler in file src_gscp.f90
!!  of the COSMO-Model. Equation numbers refer to
!!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!!
!------------------------------------------------------------------------------
!!
!! $Id: n/a$
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

MODULE gscp_kessler

!------------------------------------------------------------------------------
!>
!! Description:
!!
!!   The subroutine in this module calculates the rates of change of
!!   temperature, cloud condensate and water vapor due to cloud microphysical
!!   processes related to the formation of grid scale clouds and precipitation.
!!
!==============================================================================
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

#ifdef __COSMO__
USE kind_parameters, ONLY :   &
  wp ,       &    !! KIND-type parameter for real variables
  i4              !! KIND-type parameter for standard integer vaiables

USE data_constants  , ONLY :   &
!   pi,           & !!
!! 2. physical constants and related variables
!! -------------------------------------------
    t0=>t0_melt,  & !! melting temperature of ice
!   r_d,          & !! gas constant for dry air
    r_v,          & !! gas constant for water vapour
    rdv,          & !! r_d / r_v
    o_m_rdv,      & !! 1 - r_d/r_v
    rvd_m_o,      & !! r_v/r_d - 1
    cp_d,         & !! specific heat of dry air
    cpdr,         & !! 1 / cp_d
    lh_v,         & !! latent heat of vapourization
!   lh_f,         & !! latent heat of fusion
    lh_s,         & !! latent heat of sublimation
!   g,            & !! acceleration due to gravity
!   rho_w,        & !! density of liquid water (kg/m^3)
! 3. constants for parametrizations
!! ---------------------------------
    b1,           & !! variables for computing the saturation vapour pressure
    b2w,          & !! over water (w) and ice (i)
    b2i,          & !!               -- " --
    b3,           & !!               -- " --
    b4w,          & !!               -- " --
    b234w,        & !!               -- " --
    b4i,          & !!               -- " --
!> 4. tuning constants for radiation, cloud physics, turbulence
!! ------------------------------------------------------------
    qi0,          & !! cloud ice threshold for autoconversion
    qc0,          & !! cloud water threshold for autoconversion
!> 5. Precision-dependent security parameters (epsilons)
!! ------------------------------------------------------
    repsilon        !! precision of 1.0 in current floating point format

!! end of data_constants

USE data_runcontrol , ONLY :   &
    ldiabf_lh,      & ! include diabatic forcing due to latent heat in RK-scheme
    lsuper_coolw,   & ! switch for supercooled liquid water
    lsppt,          & ! switch, if .true., perturb the physical tendencies
    itype_qxpert_rn,& ! define which hum variables tend. are perturbed
    itype_qxlim_rn    ! type of reduction/removal of the perturbation 
                      ! in case of negative (qv, qc, qi) or 
                      ! supersaturated (qv) values

!------------------------------------------------------------------------------
!! COSMO environment modules
!------------------------------------------------------------------------------

USE utilities,                ONLY : message
USE meteo_utilities,          ONLY : satad       !! saturation adjustment
USE src_stoch_physics,        ONLY : apply_tqx_tend_adj
#endif
! of ifdef __COSMO__

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,           ONLY :  &
    llhn,         & ! main switch for latent heat nudging
    llhnverif       ! main switch for latent heat nudging
#endif

!------------------------------------------------------------------------------

#ifdef __ICON__
USE mo_kind,               ONLY: wp         , &
                                 i4
!USE mo_math_constants    , ONLY: pi
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 lh_v  => alv   , & !! latent heat of vapourization
!                                lh_f  => alf   , & !! latent heat of fusion
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
!                                rho_w => rhoh2o, & !! density of liquid water (kg/m^3)
                                 t0    => tmelt     !! melting temperature of ice/snow

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b4w   => c4les     !!               -- " --
USE mo_satad,              ONLY: satad_v_3d,     &  !! new saturation adjustment
                                 sat_pres_water     !! saturation vapor pressure w.r.t. water
USE mo_exception,          ONLY: message, message_text
#endif

!------------------------------------------------------------------------------

! this can be used by ICON and COSMO
USE gscp_data, ONLY: &          ! all variables are used here

    zconst,   zbev,    zvz0r,                 &
    x7o8,                                     &
    zkcac,   zkphi1,    zkphi2,    zkphi3,    &
    x1o8,      x3o16,   iautocon

#ifdef __ICON__
! this is (at the moment) an ICON part
USE gscp_data, ONLY: &          ! all variables are used here
    kc_c2
#endif

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: kessler

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef __COSMO__
CHARACTER(132) :: message_text = ''
#endif

!==============================================================================

CONTAINS

!==============================================================================
!> Module procedure "kessler" in "gscp_kessler" for computing effects of 
!!  grid scale precipitation including cloud water, cloud ice, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE kessler  (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qr,                  & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_gscp_data.f90
  qc0,                               & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,                           & !! surface precipitation rates
#ifdef __COSMO__
  tinc_lh,                           & !  t-increment due to latent heat 
#endif
#ifdef NUDGING
  tt_lheat,                          & !  t-increments due to latent heating (nud) 
#endif
  l_cv,                              &
  ddt_tend_t     , ddt_tend_qv     , &
  ddt_tend_qc    ,                   & !> ddt_tend_xx are tendencies
  ddt_tend_qr    ,                   & !!    necessary for dynamics
  ddt_diag_au    , ddt_diag_ac     , & !!
  ddt_diag_ev                        ) !! ddt_diag_xxx are optional

!------------------------------------------------------------------------------
!>
!! Description:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. The
!!   variables are updated in this subroutine. Rain and snow are prognostic
!!   variables. The precipitation fluxes at the surface are stored on the
!!   corresponding global fields.
!!
!! Method:
!!
!------------------------------------------------------------------------------
!! Declarations:
!!
!------------------------------------------------------------------------------
!! Modules used: These are declared in the module declaration section
!! -------------

!! Subroutine arguments:
!! --------------------

  INTEGER, INTENT(IN) ::  &
    nvec          ,    & !> number of horizontal points
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart   ,    & !> optional start index for horizontal direction
    ivend     ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=wp), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

#ifdef __ICON__
  REAL(KIND=wp), INTENT(IN) :: &
    qc0              !> cloud ice/water threshold for autoconversion
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(in)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(IN) ::      &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
#endif
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(INOUT) ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
#endif
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qr                     !! specific rain content                         (kg/kg)

#ifdef __COSMO__
  REAL(KIND=wp), INTENT(INOUT) :: &
       tinc_lh(:,:)    ! temperature increments due to heating             ( K/s ) 
#endif  

#ifdef NUDGING
  REAL(KIND=wp), INTENT(INOUT) :: &
       tt_lheat(:,:)       !  t-increments due to latent heating (nudg) ( K/s )
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=wp), DIMENSION(nvec), INTENT(INOUT) ::   &
#else
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
#endif
    prr_gsp                !> precipitation rate of rain, grid-scale        (kg/(m2*s))

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
#endif
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qr         !! tendency qr                                      ( 1/s )

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
#endif
    ddt_diag_au     , & !> optional output autoconversion rate cloud to rain           ( 1/s )
    ddt_diag_ac     , & !! optional output accretion rate cloud to rain                ( 1/s )
    ddt_diag_ev         !! optional output evaporation of rain                         ( 1/s )

  !! Local parameters: None, parameters are in module header, gscp_data or data_constants
  !! ----------------
  
  REAL    (KIND=wp), PARAMETER ::  &
    ! basic constants of the parameterization scheme
    zaau   = 1.0_wp/1000.0_wp,            & ! coef. for autoconversion
    zaac   = 1.72_wp,                     & ! coef. for accretion (neu)
    zbev   = 9.0_wp,                      & ! coef. for drop ventilation,
    ! to avoid illegal operations in power expressions
    znull = 1.E-20_wp

  !> Local scalars:
  !! -------------
  
  INTEGER (KIND=i4) :: &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL (KIND=wp)   ::  &
    ztx   ,            & ! 
    zpx   ,            & !
    fsa3  ,            & !
    zspw  ,            & ! equilibrium vapour pressure over water
    zsa3  ,            & !
    zc3   ,            & !
    zc1c  ,            & !
    zc1   ,            & !
    zx    ,            & !
    zsrmax

  REAL (KIND=wp)   ::  &
    zsqvw ,            & !  specific humidity at water saturation
    zqvts ,            & !  qv-tendency in a layer
    zqcts ,            & !  qc-tendency in a layer
    zqrts ,            & !  qr-tendency in a layer
    ztts  ,            & !  t -tendency in a layer
    zswra ,            & !  autoconversion rate
    zswrk ,            & !  accretion rate
    zsrd                 !  evaporation rate

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in          ,    & !> temperature                                   (  K  )
    qv_in         ,    & !! specific water vapor content                  (kg/kg)
    qc_in         ,    & !! specific cloud water content                  (kg/kg)
    qr_in                !! specific rain content                         (kg/kg)

  REAL (KIND=wp)   ::  &
    fpvsw, fqvs  ,     & !
    zpv          ,     & !
    zdtdh, zphi  ,     & !
    zqrk, lnzqrk ,     & !
    zdtr, ztau   ,     & ! 
    zimr, zzar   ,     & !
    qvg          ,     & !  specific water vapor content: local grid cell value
    qcg          ,     & !  specfic cloud water content:  - "" -
    tg           ,     & !  temperature:                  - "" -
    qrg          ,     & !  specific rain content:        - "" -
    rhog         ,     & !  density                       - "" -
    rhogr        ,     & !  reciprocal density            - "" -
    ppg                  !  full level pressure           - "" -


!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::  &
    zvzr        (nvec),     & !
    zpkr        (nvec),     & !
    zprvr       (nvec),     & !
#ifdef __COSMO__
    zdummy      (nvec,8),   & !
#endif
    zpkm1r      (nvec)        !

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine kessler
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )
  fsa3(ztx) = 3.86E-3_wp - 9.41E-5_wp*(ztx-t0)

  
!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!------------------------------------------------------------------------------

! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)

#ifdef __COSMO__
  z_heat_cap_r = cpdr
#endif

#ifdef __ICON__
  IF (PRESENT(l_cv)) THEN
    IF (l_cv) THEN
      z_heat_cap_r = cvdr
    ELSE
      z_heat_cap_r = cpdr
    ENDIF
  ELSE
    z_heat_cap_r = cpdr
  ENDIF
#endif

  ! Delete precipitation fluxes from previous timestep
  prr_gsp (:) = 0.0_wp
  zpkr    (:) = 0.0_wp
  zprvr   (:) = 0.0_wp
  zvzr    (:) = 0.0_wp
#ifdef __COSMO__
  ! CDIR COLLAPSE
  zdummy(:,:)=0.0_wp
#endif
  ! Optional arguments

  IF (PRESENT(ddt_tend_t)) THEN
    ! save input arrays for final tendency calculation
    t_in(:,:)  = t(:,:)
    qv_in(:,:) = qv(:,:)
    qc_in(:,:) = qc(:,:)
    qr_in(:,:) = qr(:,:)
  END IF

  IF (PRESENT(ivstart)) THEN
    iv_start = ivstart
  ELSE
    iv_start = 1
  END IF
  IF (PRESENT(ivend)) THEN
    iv_end = ivend
  ELSE
    iv_end = nvec
  END IF
  IF (PRESENT(kstart)) THEN
    k_start = kstart
  ELSE
    k_start = 1
  END IF
  IF (PRESENT(idbg)) THEN
    izdebug = idbg
  ELSE
    izdebug = 0
  END IF

! timestep for calculations
  zdtr  = 1.0_wp / zdt

#ifdef NUDGING
  ! add part of latent heating calculated in subroutine kessler to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (llhn) THEN
    ! CALL get_gs_lheating ('add',1,ke) !XL :should not be called from block physics
    tt_lheat(:,:) = tt_lheat(:,:) - t(:,:)
  ENDIF
#endif

! output for various debug levels
  IF (izdebug > 15) CALL message('','gscp_kessler:  Start of kessler')
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   nvec = ',nvec       ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke           ; CALL message('',message_text)
    WRITE (message_text,*) '   ivstart = ',ivstart ; CALL message('',message_text)
    WRITE (message_text,*) '   ivend   = ',ivend   ; CALL message('',message_text)
  END IF

  IF (izdebug > 50) THEN
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN dz  = ',MAXVAL(dz),MINVAL(dz)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN T   = ',MAXVAL(t),MINVAL(t)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN p   = ',MAXVAL(p),MINVAL(p)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN rho = ',MAXVAL(rho),MINVAL(rho)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qv  = ',MAXVAL(qv),MINVAL(qv)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qc  = ',MAXVAL(qc),MINVAL(qc)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qr  = ',MAXVAL(qr),MINVAL(qr)
    CALL message('',message_text)
  ENDIF

! ----------------------------------------------------------------------
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! ----------------------------------------------------------------------

loop_over_levels: DO k = 1, ke

#ifdef __COSMO__
  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat
    tinc_lh(:,k) = tinc_lh(:,k) - t(:,k)
  ENDIF
#endif

  !----------------------------------------------------------------------------
  !  Section 2: Test for clouds and precipitation in present layer.
  !             If no cloud or precipitation points have been found
  !             go to the next layer.
  !----------------------------------------------------------------------------

  loop_over_i: DO iv = iv_start, iv_end

    IF (qr(iv,k) < 1.0E-15_wp) qr(iv,k) = 0.0_wp
    IF (qc(iv,k) < 1.0E-15_wp) qc(iv,k) = 0.0_wp

    qcg   = qc(iv,k)
    qvg   = qv(iv,k)
    qrg   = qr(iv,k)
    tg    = t(iv,k)
    ppg   = p(iv,k)
    rhog  = rho(iv,k)
    rhogr = 1.0_wp / rhog

    zqrk  = qrg * rhog
    zdtdh = 0.5_wp * zdt / dz(iv,k)

    zpkm1r(iv) = zpkr(iv)
    zzar      = zqrk/zdtdh + zprvr(iv) + zpkm1r(iv)
    IF (zqrk > znull) THEN
      zpkr(iv) = zqrk * zvz0r * EXP(x1o8*LOG(zqrk))
    ELSE
      zpkr(iv) = 0.0_wp
    ENDIF
    zpkr(iv) = MIN( zpkr(iv), zzar )
    zzar    = zdtdh * (zzar-zpkr(iv))

    IF (zvzr(iv) == 0.0_wp) zvzr(iv) = zvz0r * EXP(x1o8*LOG(MAX(0.5_wp*zqrk,znull)))

    zimr    = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
    zqrk    = zzar*zimr

    IF (zqrk > znull) THEN
      lnzqrk = LOG(zqrk)
    ELSE
      lnzqrk = 0.0_wp
    ENDIF

    ztts  = 0.0_wp
    zqvts = 0.0_wp
    zqcts = 0.0_wp

    !------------------------------------------------------------------------
    !  Section 5: Calculation of cloud microphysics for cloud case
    !             ( qc > 0)
    !------------------------------------------------------------------------

    IF (qcg > 0.0_wp) THEN

      ! Calculate conversion rates

      IF (iautocon == 0) THEN
        ! Coefficients
        zc1c  = zaac *  EXP(lnzqrk*x7o8)
        zc1   = zaau + zc1c
        zx    = qcg / (1.0_wp + zc1*zdt)

        ! Conversion rates
        zswra  = zaau * zx
        zswrk  = zc1c * zx
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num

        IF (qcg > 1.0E-6_wp) THEN
          ztau  = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
          zphi  = zkphi1 * ztau**zkphi2 * (1.0_wp - ztau**zkphi2)**3
          zswra = zconst * qcg*qcg*qcg*qcg &                      ! S_au
                 * (1.0_wp + zphi/(1.0_wp - ztau)**2)
          zphi  = (ztau/(ztau+zkphi3))**4
          zswrk = zkcac * qcg * qrg * zphi !* zrho1o2(i)          ! S_ac
        ELSE
          zswra = 0.0_wp  ! S_au                    
          zswrk = 0.0_wp  ! S_ac
        ENDIF
      ENDIF

      ! Store tendencies
      zqcts  = - zswra - zswrk
      zqrts  =   zswra + zswrk

      ! Update values
      qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)

    !------------------------------------------------------------------------
    !  Section 7: Calculation of cloud microphysics for
    !             precipitation case without cloud ( qc = 0 )
    !------------------------------------------------------------------------

    ELSEIF ( (zqrk) > 0.0_wp .AND. qcg <= 0.0_wp )    THEN

#ifdef __COSMO__
      zspw  = fpvsw(tg)
      zsqvw = fqvs(zspw, ppg)
#endif

#ifdef __ICON__
      zsqvw = sat_pres_water(tg)/(rhog * r_v *tg)
#endif

      ! Coefficients
      zsrmax = zzar*rhogr * zdtr
      zsa3   = fsa3(tg)
      zc3    = zsa3 * SQRT(zqrk) * (1.0_wp + zbev * EXP(lnzqrk*x3o16))

      ! Conversion rates
      zsrd   = -zc3 * (qvg - zsqvw)
      zsrd   = MIN (zsrmax, zsrd)

      ! Store tendencies
      zqvts =   zsrd
      ztts  = - lh_v*zsrd*z_heat_cap_r
      zqrts = - zsrd

      ! Update values
      qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)

    ENDIF

    !------------------------------------------------------------------------
    ! Section 8: Complete time step
    !------------------------------------------------------------------------

    IF ( k /= ke ) THEN
      ! Store precipitation fluxes and sedimentation velocities for the next level
      zprvr(iv) = qrg*rhog*zvzr(iv)
      zvzr(iv)  = zvz0r * EXP(x1o8 * LOG(MAX((qrg+qr(iv,k+1))*0.5_wp*rhog,znull)))
    ELSE
      ! Precipitation flux at the ground
      prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
    ENDIF

    ! Update of prognostic variables or tendencies
    qr (iv,k) = qrg
!   qrs(iv,k) = qrg
    t  (iv,k) = t (iv,k) + ztts*zdt
    qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvts*zdt )
    qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqcts*zdt )

#ifndef __ICON__
    ! Store optional microphysical rates for diagnostics
    IF (PRESENT(ddt_diag_au)) ddt_diag_au(:,k) = zswra
    IF (PRESENT(ddt_diag_ac)) ddt_diag_ac(:,k) = zswrk
    IF (PRESENT(ddt_diag_ev)) ddt_diag_ev(:,k) = zsrd
#endif

  ENDDO loop_over_i

  IF (izdebug > 25) THEN
    ! Check for negative values
    DO iv = iv_start, iv_end
      IF (qr(iv,k) < 0.0_wp) THEN
        WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qr'
        CALL message('',message_text)
      ENDIF
      IF (qc(iv,k) < 0.0_wp) THEN
        WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qc'
        CALL message('',message_text)
      ENDIF
      IF (qv(iv,k) < 0.0_wp) THEN
        WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qv'
        CALL message('',message_text)
      ENDIF
    ENDDO
  ENDIF

#if defined (__COSMO__)
  ! Do a final saturation adjustment for new values of t, qv and qc
  CALL satad ( 1, t(:,k), qv(:,k),          &
    qc(:,k), t(:,k), p(:,k),                &
    zdummy(:,1),zdummy(:,2),zdummy(:,3),    &
    zdummy(:,4),zdummy(:,5),zdummy(:,6),    &
    zdummy(:,7),zdummy(:,8),                &
    b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,  &
    rvd_m_o, lh_v, cpdr, cp_d,              &
    nvec, 1, iv_start, iv_end, 1 , 1 )

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
    tinc_lh(:,k) = tinc_lh(:,k) + t(:,k)
  ENDIF
#endif

ENDDO loop_over_levels

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydci to model latent
! heating field: add temperature to model latent heating field
IF (llhn .OR. llhnverif) &
! CALL get_gs_lheating ('inc',1,ke)  !XL :this should be called from within the block
     tt_lheat(:,:) = tt_lheat(:,:) + t(:,:)
#endif

#ifdef __ICON__

 CALL satad_v_3d (                             &
               & maxiter  = 10_i4        ,& !> IN
               & tol      = 1.e-3_wp     ,& !> IN
               & te       = t            ,&
               & qve      = qv           ,&
               & qce      = qc           ,&
               & rhotot   = rho          ,&
               & idim     = nvec         ,&
               & kdim     = ke           ,&
               & ilo      = iv_start     ,&
               & iup      = iv_end       ,&
               & klo      = k_start      ,&
               & kup      = ke            &
!              !& count, errstat,
               )
#endif

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

  IF (PRESENT(ddt_tend_t)) THEN

    DO k=k_start,ke
      DO iv=iv_start, iv_end

        ! calculated pseudo-tendencies
        ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
        ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
        ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
        ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)

        ! restore input values
        t (iv,k) = t_in (iv,k)
        qv(iv,k) = qv_in(iv,k)
        qc(iv,k) = qc_in(iv,k)
        qr(iv,k) = qr_in(iv,k)

      END DO
    END DO

  END IF

  IF (izdebug > 25) THEN
    CALL message('mo_gscp', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp T= ' ,&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine kessler
!------------------------------------------------------------------------------

END SUBROUTINE kessler

!==============================================================================

END MODULE gscp_kessler
