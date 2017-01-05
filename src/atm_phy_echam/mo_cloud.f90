#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS STRICT
#endif
#include "fsel.inc"

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module computes large-scale water phase changes
!!
!! @remarks
!!          This routine computes the tendencies of the four prognostic
!!          variables (temperature t, specific humidity q, cloud liquid
!!          water xl, cloud ice xi) due to phase changes (condensation/
!!          deposition, evaporation/sublimation of rain/snow falling
!!          into the unsaturated part of the grid box, melting of snow,
!!          melting/freezing of cloud ice/cloud water, sedimentation of
!!          cloud ice, and precipitation formation in warm, cold and
!!          mixed phase clouds.
!!          The precipitation at the surface (rain and snow) is used in
!!          later for computing the land surface hydrology in *surf*.
!!          The cloud parameters (cloud cover, cloud liquid water and
!!          cloud ice are used for the calculation of radiation at the
!!          next timestep.
!!          Attention:
!!          In the current version the advective tendencies of skewness
!!          and variance are set to zero.
!!
!! @references.
!!     Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
!!     Levkov et al., 1992: Beitr. Phys. Atm. 35-58.          (ice phase)
!!     Beheng, 1994: Atmos. Res. 193-206.                    (warm phase)
!!     Lenderink et al., 1998; KNMI-REPORT NO. 98-13       (condensation)
!!     Tompkins 2002, J. Atmos. Sci.                        (cloud cover)
!!
!! @author U. Lohmann,   MPI-M, Hamburg 1995
!!         G. Lenderink, KNMI,  de Bilt 1998
!!         M. Esch,      MPI-M, Hamburg 1999
!!
!! @par Revision History
!!       - E. Roeckner,  MPI-M, Hamburg  2000
!!       - A. Tompkins,  MPI-M, Hamburg  2000
!!       - U. Schlese,   MPI-M, Hamburg  2003
!!       - P. Stier,     MPI-M, Hamburg  2002-2004: Scavenging parameters
!! - Taken from ECHAM6.2, wrapped in module and modified for ICON
!!   by Monika Esch, MPI-M (2013-11)
!! - updated to ECHAM.6.3 and unified for use in ICON; removed Tompkins scheme
!!   by Monika Esch, MPI-M (2015-05)
!!
!
MODULE mo_cloud

  USE mo_kind,                 ONLY : wp
  USE mo_math_constants,       ONLY : pi
  USE mo_physical_constants,   ONLY : rd, alv, als, rv, vtmpc1, tmelt, rhoh2o
  USE mo_echam_convect_tables, ONLY : prepare_ua_index_spline, lookup_ua_spline      &
                                    , lookup_uaw_spline, lookup_ubc                  &
                                    , lookup_ua_eor_uaw_spline
  USE mo_echam_cloud_config,   ONLY : echam_cloud_config

#ifdef _PROFILE
  USE mo_profile,              ONLY : trace_start, trace_stop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cloud

  ! to simplify access to components of echam_cloud_config
  LOGICAL , POINTER :: lonacc
  INTEGER , POINTER :: jbmin, jbmax
  REAL(wp), POINTER :: cqtmin, cvtfall, crhosno, cn0s   , cthomi , csecfrl, cauloc, &
       &               clmax , clmin  , ccraut , ceffmin, ceffmax, crhoi  ,         &
       &               ccsaut, ccsacl , ccracl , ccwmin , clwprat


CONTAINS
  !>
  !!
  !!
  SUBROUTINE cloud (         kproma,       kbdim,          ktdia                     &
                           , klev,         klevp1                                    &
                           , pdtime                                                  &
    ! - INPUT  1D .
                           , kctop                                                   &
    ! - INPUT  2D .
                           , paphm1                                                  &
                           , papm1                                                   &
                           , pdz                                                     &
                           , pmdry                                                   &
                           , pacdnc                                                  &
                           , ptm1,         prho                                      &
                           , pqm1,         pxlm1,        pxim1                       &
                           , pcair                                                   &
    ! - INPUT/OUTPUT 1D .
                           , paclcov                                                 &
                           , ktype                                                   &
                           , pch_concloud, pcw_concloud                              &
    ! - INPUT/OUTPUT 2D .
                           , pxtecl,       pxteci,       pqtec                       &
                           , ptte                                                    &
                           , pqte,         pxlte,        pxite                       &
                           , paclc                                                   &
    ! - OUTPUT 1D .
                           , pssfl,        prsfl                                     &
    ! - OUTPUT 2D .
                           , prelhum                                                 &
                           , ptte_prc,     pqte_prc                                  &
                           , pxlte_prc,    pxite_prc                                 &
                            )
    !
    !
    !
    INTEGER,  INTENT(IN)    :: kproma, kbdim, ktdia, klev, klevp1
    INTEGER,  INTENT(IN)    :: kctop(kbdim)
    INTEGER,  INTENT(INOUT) :: ktype(kbdim)
    REAL(wp), INTENT(IN)    :: pdtime
    REAL(wp), INTENT(IN)    ::     &
      & paphm1   (kbdim,klevp1)   ,&!< pressure at half levels                   (n-1)
      & papm1    (kbdim,klev)     ,&!< pressure at full levels                   (n-1)
      & pdz      (kbdim,klev)     ,&!< geometric height thickness of layer
      & pmdry    (kbdim,klev)     ,&!< dry air content
      & pacdnc   (kbdim,klev)     ,&!< cloud droplet number concentration (specified)
      & ptm1     (kbdim,klev)     ,&!< temperature                               (n-1)
      & prho     (kbdim,klev)     ,&!< air density
      & pqm1     (kbdim,klev)     ,&!< specific humidity                         (n-1)
      & pxlm1    (kbdim,klev)     ,&!< cloud liquid water                        (n-1)
      & pxim1    (kbdim,klev)     ,&!< cloud ice                                 (n-1)
      & pcair    (kbdim,klev)       !< specific heat of moist air
    REAL(wp), INTENT(INOUT) ::     &
      & paclcov  (kbdim)          ,&!< total cloud cover
      & pch_concloud(kbdim)       ,&!< checking the global heat budget and
      & pcw_concloud(kbdim)         !< checking the local water budget
                                    !< due to convection + cloud processes
    REAL(wp), INTENT(INOUT) ::     &
      & pxtecl   (kbdim,klev)     ,&!< detrained convective cloud liquid water   (n)
      & pxteci   (kbdim,klev)     ,&!< detrained convective cloud ice            (n)
      & pqtec    (kbdim,klev)     ,&!<
      & ptte     (kbdim,klev)     ,&!< tendency of temperature
      & pqte     (kbdim,klev)     ,&!< tendency of specific humidity
      & pxlte    (kbdim,klev)     ,&!< tendency of cloud liquid water
      & pxite    (kbdim,klev)     ,&!< tendency of cloud ice
      & paclc    (kbdim,klev)       !< cloud cover  (now diagnosed in cover)
    REAL(wp),INTENT(INOUT) ::      &
      & ptte_prc(kbdim,klev)      ,&!<
      & pqte_prc(kbdim,klev)        ! OUT
    REAL(wp),INTENT(INOUT) ::      &
      & pxlte_prc(kbdim,klev)     ,&!<
      & pxite_prc(kbdim,klev)       ! OUT
    REAL(wp), INTENT(INOUT)   ::   &! use INOUT to preserve the initialization
      & prsfl    (kbdim)          ,&!< surface rain flux
      & pssfl    (kbdim)            !< surface snow flux
    REAL(wp), INTENT(OUT)   ::     &
      & prelhum  (kbdim,klev)       !< relative humidity
    !
    !   Temporary arrays
    !
    REAL(wp):: zclcpre(kbdim)       ,zclcpre_inv(kbdim)                              &
      &      , zcnd(kbdim)          ,zdep(kbdim)          ,zdp(kbdim)                &
      &      , zevp(kbdim)          ,zxievap(kbdim)       ,zxlevap(kbdim)            &
      &      , zfrl(kbdim)          ,zimlt(kbdim)         ,zsmlt(kbdim)              &
      &      , zrpr(kbdim)          ,zspr(kbdim)          ,zsub(kbdim)               &
      &      , zxiflux(kbdim)       ,zclcauxi(kbdim)      ,zdqsat1(kbdim)            &
      &      , zsacl(kbdim)         ,zqp1(kbdim)                                     &
      &      , zlsdcp(kbdim)        ,zlvdcp(kbdim)        ,zcoeff(kbdim)             &
      &      , ztp1tmp(kbdim)       ,zqp1tmp(kbdim)                                  &
      &      , zrfl(kbdim)          ,zsfl(kbdim)          ,ztp1(kbdim)               &
      &      , zxlb(kbdim)          ,zxib(kbdim)          ,zqrho(kbdim)              &
      &      , zqrho_sqrt(kbdim)    ,zpapm1_inv(kbdim)    ,zpapp1i(kbdim)            &
      &      , zclcov(kbdim)        ,zclcaux(kbdim)       ,zqsed(kbdim)              &
      &      , zxlvi(kbdim)         ,zxlvitop(kbdim)      ,zxlvibot(kbdim)           &
      &      , zxrp1(kbdim)         ,zxsp1(kbdim)         ,zxsp2(kbdim)              &
      &      , zgenti(kbdim)        ,zgentl(kbdim)        ,zauloc(kbdim)             &
      &      , zqsi(kbdim)          ,ztmp1(kbdim)         ,ztmp2(kbdim)              &
      &      , ztmp3(kbdim)         ,ztmp4(kbdim)         ,zxised(kbdim)             &
      &      , zqvdt(kbdim)         ,zqsm1(kbdim)         ,zdtdt(kbdim)              &
      &      , zstar1(kbdim)        ,zlo2(kbdim)          ,za(kbdim)                 &
      &      , ub(kbdim)            ,ua(kbdim)            ,dua(kbdim)                &
      &      , uaw(kbdim)           ,duaw(kbdim)          ,zclten(kbdim)             &
      &      , zcpten(kbdim,klev)   ,zqviten(kbdim)       ,zqten(kbdim,klev)

    INTEGER:: loidx(kbdim), nloidx(kbdim), jjclcpre(kbdim)
    INTEGER:: cond1(kbdim), cond2(kbdim)
    INTEGER:: idx1(kbdim), idx2(kbdim)

    INTEGER:: nclcpre
    INTEGER:: jl, jk, nl, locnt, nlocnt, nphase, i1 , i2, klevtop
    LOGICAL   lo
    !!$ used in Revised Bergeron-Findeisen process only
    !!$  LOGICAL   locc

    REAL(wp):: zdqsat, zqcdif, zfrho, zifrac, zepsec, zxsec                          &
      &      , zqsec, zrc, zcons, ztdif, zsnmlt, zclcstar                            &
      &      , zesi, zsusati, zb1, zb2, zcfac4c, zzeps, zesw, zesat                  &
      &      , zqsw, zsusatw, zdv, zast, zbst, zzepr, zxip1, zxifall, zal1, zal2     &
      &      , zlc, zdqsdt, zlcdqsdt, zdtdtstar, zxilb, zrelhum                      &
      &      , zes, zcor, zqsp1tmp, zoversat, zqcon, zdepos                          &
      &      , zcond, zradl, zf1, zraut, zexm1, zexp, zrac1, zrac2, zrieff           &
      &      , zcolleffi, zc1, zdt2, zsaut, zsaci1, zsaci2, zsacl1, zsacl2, zlamsm   &
      &      , zzdrr, zzdrs, zpretot, zpredel, zpresum, zxlp1                        &
      &      , zxlold, zxiold, zdxicor, zdxlcor, zptm1_inv, zxlp1_d, zxip1_d         &
      &      , zupdate, zlo, zcnt, zclcpre1, zval, zua, zdua, zxitop, zxibot         &
      &      , zqvte, zxlte, zxite, ztte
    !!$ used in Revised Bergeron-Findeisen process only
    !!$  REAL(wp):: zzevp, zeps
    !!$  REAL(wp):: zsupsatw(kbdim)
    !
    ! mpuetz : the following tendencies don't have to be vectors
    !
    REAL(wp):: zmratepr(kbdim,klev), & ! Rain formation rate in cloudy part of the grid box [kg/kg]
      &        zmrateps(kbdim,klev), & ! Ice  formation rate in cloudy part of the grid box  [kg/kg]
      &        zfrain(kbdim,klev),   & ! Rain flux before evaporation [kg/m2/s]
      &        zfsnow(kbdim,klev),   & ! Snow flux before sublimation [kg/m2/s]
      &        zfevapr(kbdim,klev),  & ! Evaporation of rain [kg/m2/s]
      &        zfsubls(kbdim,klev),  & ! Sublimation of snow [kg/m2/s]
      &        zmlwc(kbdim,klev),    & ! In-cloud liquid water mass mixing ratio before rain formation [kg/kg]
      &        zmiwc(kbdim,klev),    & ! In-cloud ice mass mixing ratio before snow formation [kg/kg]
      &        zmsnowacl(kbdim,klev)   ! Accretion rate of snow with cloud droplets
                                       ! in cloudy part of the grid box  [kg/kg]

    ! to simplify access to components of echam_cloud_config
    lonacc   => echam_cloud_config% lonacc
    jbmin    => echam_cloud_config% jbmin
    jbmax    => echam_cloud_config% jbmax
    cqtmin   => echam_cloud_config% cqtmin
    cvtfall  => echam_cloud_config% cvtfall
    crhosno  => echam_cloud_config% crhosno
    cn0s     => echam_cloud_config% cn0s
    cthomi   => echam_cloud_config% cthomi
    csecfrl  => echam_cloud_config% csecfrl
    cauloc   => echam_cloud_config% cauloc
    clmax    => echam_cloud_config% clmax
    clmin    => echam_cloud_config% clmin
    ccraut   => echam_cloud_config% ccraut
    ceffmin  => echam_cloud_config% ceffmin
    ceffmax  => echam_cloud_config% ceffmax
    crhoi    => echam_cloud_config% crhoi
    ccsaut   => echam_cloud_config% ccsaut
    ccsacl   => echam_cloud_config% ccsacl
    ccracl   => echam_cloud_config% ccracl
    ccwmin   => echam_cloud_config% ccwmin
    clwprat  => echam_cloud_config% clwprat

    zmratepr(:,:) = 0._wp
    zmrateps(:,:) = 0._wp
    zfrain(:,:)   = 0._wp
    zfsnow(:,:)   = 0._wp
    zfevapr(:,:)  = 0._wp
    zfsubls(:,:)  = 0._wp
    zmlwc(:,:)    = 0._wp
    zmiwc(:,:)    = 0._wp
    zmsnowacl(:,:)= 0._wp
    ! special treatment for prelhum: actual parameter is stream element and must
    !                                be defined for all array elements.
    prelhum(:,:)  = 0._wp

    ! save the tendencies accumulated before calling this routine

    ptte_prc(1:kproma,:)   =  ptte(1:kproma,:)
    pqte_prc(1:kproma,:)   =  pqte(1:kproma,:)
    pxlte_prc(1:kproma,:)  = pxlte(1:kproma,:)
    pxite_prc(1:kproma,:)  = pxite(1:kproma,:)
    !
    ! Executable statements
    !
#ifdef _PROFILE
    CALL trace_start ('cloud', 10)
#endif

    !
    !   Security parameters
    !
    zepsec = 1.0e-12_wp
    zxsec  = 1.0_wp-zepsec
    zqsec  = 1.0_wp-cqtmin
    !  zeps   = EPSILON(1.0_wp)
    !
    !   Computational constants
    !
    !     ----------------------------------------------------------------------------
    !
    !       1.   Top boundary conditions, air density and geopotential
    !            height at half levels
    !            Set to zero precipitation fluxes etc.
    !
    nclcpre = 0

    DO 111 jl = 1,kproma
       zclcpre(jl)   = 0.0_wp
       zxiflux(jl)   = 0.0_wp
       zrfl(jl)      = 0.0_wp
       zsfl(jl)      = 0.0_wp
111 END DO

    !
    DO 831 jk = ktdia,klev  ! the big jk-loop
    !
    !       1.3   Air density
    !
!IBM* NOVECTOR
      DO jl = 1,kproma
         zqrho(jl)        = 1.3_wp/prho(jl,jk)
         pxtecl(jl,jk)    = MAX(pxtecl(jl,jk),0.0_wp)
         pxteci(jl,jk)    = MAX(pxteci(jl,jk),0.0_wp)
         pqtec(jl,jk)     = MAX(pqtec(jl,jk),0.0_wp)
      END DO

      zqrho_sqrt(1:kproma) = SQRT(zqrho(1:kproma))
      zpapm1_inv(1:kproma) = 1._wp/papm1(1:kproma,jk)

      CALL prepare_ua_index_spline('cloud (1)',kproma,ptm1(1,jk),loidx(1),za(1))
      CALL lookup_ua_spline(kproma,loidx(1),za(1),ua(1),dua(1))
      CALL lookup_uaw_spline(kproma,loidx(1),za(1),uaw(1),duaw(1))
      !
      !     --------------------------------------------------------------------------
      !       2.    Set to zero some local tendencies (increments)
      !     --------------------------------------------------------------------------
      !
!IBM* NOVECTOR
      DO 201 jl = 1,kproma

         zevp(jl)       = 0.0_wp
         zsub(jl)       = 0.0_wp
         zqsed(jl)      = 0.0_wp
         zxlevap(jl)    = 0.0_wp
         zxievap(jl)    = 0.0_wp

         zdp(jl)        = paphm1(jl,jk+1)-paphm1(jl,jk)

         zrc            = 1._wp/pcair(jl,jk)
         zlvdcp(jl)     = alv*zrc
         zlsdcp(jl)     = als*zrc
201   END DO
      !
      !     --------------------------------------------------------------------------
      !
      !       3.   Modification of incoming precipitation fluxes by
      !            melting, sublimation and evaporation
      !
#ifdef _PROFILE
      CALL trace_start ('cloud_loop_3', 13)
#endif
      IF (jk .GT. 1) THEN
        !
        !       3.1   Melting of snow and ice
        !
!IBM* NOVECTOR
        DO 321 jl = 1,kproma

           zcons     = (pmdry(jl,jk)/pdtime)/(zlsdcp(jl)-zlvdcp(jl))
           ztdif     = MAX(0.0_wp,ptm1(jl,jk)-tmelt)
           zsnmlt    = MIN(zxsec*zsfl(jl),zcons*ztdif)
           zrfl(jl)  = zrfl(jl)+zsnmlt
           zsfl(jl)  = zsfl(jl)-zsnmlt
           zsmlt(jl) = zsnmlt/pmdry(jl,jk)*pdtime
           zsnmlt    = MAX(0.0_wp,pxim1(jl,jk)+(pxite(jl,jk)+pxteci(jl,jk))*pdtime)
           zimlt(jl) = FSEL(-ztdif,0.0_wp,zsnmlt)
321     END DO

        IF (nclcpre.GT.0) THEN
        ! equals old zclcpre.gt.0
        !
        !       3.2   Sublimation of snow and ice (Lin et al., 1983)
        !
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
          DO nl = 1,nclcpre
             jl = jjclcpre(nl)
             zesi     = ua(jl)*zpapm1_inv(jl)
             zesi     = MIN(zesi,0.5_wp)
             zqsi(jl) = zesi/(1._wp-vtmpc1*zesi)
             zsusati  = MIN(pqm1(jl,jk)/zqsi(jl)-1.0_wp,0.0_wp)
             zb1      = zlsdcp(jl)**2/(2.43e-2_wp*rv*(ptm1(jl,jk)**2))
             zb2      = 1._wp/(prho(jl,jk)*zqsi(jl)*0.211e-4_wp)
             zcoeff(jl) = 3.e6_wp*2._wp*pi*(zsusati/(prho(jl,jk)*(zb1+zb2)))
             zclcpre_inv(jl) = 1._wp/zclcpre(jl)
          END DO
          !
          ! definition of conditions for the sublimation of snow and ice
!IBM* ASSERT(NODEPS)
          DO nl = 1,nclcpre
             jl = jjclcpre(nl)
             cond1(nl) = INT(FSEL(cqtmin-zsfl(jl)   ,0._wp,1._wp))
             cond2(nl) = INT(FSEL(cqtmin-zrfl(jl)   ,0._wp,1._wp))
          END DO
          i1 = 1
          i2 = 1
          DO nl = 1,nclcpre
             jl = jjclcpre(nl)
             idx1(i1) = jl
             i1 = i1 + cond1(nl)
             idx2(i2) = jl
             i2 = i2 + cond2(nl)
          END DO
          i1 = i1 - 1
          i2 = i2 - 1
          !    old if(zsfl(jl).GT.ctqmin)
          !
          IF (i1.GT.0) THEN
!IBM* ASSERT(NODEPS)
            DO nl = 1,i1
               jl = idx1(nl)
               ztmp1(nl)    = zqrho_sqrt(jl)
               ztmp2(nl)    = zsfl(jl)*zclcpre_inv(jl)/cvtfall
            END DO
            ztmp1(1:i1)   = SQRT(ztmp1(1:i1))
            ztmp2(1:i1)   = ztmp2(1:i1)**(1._wp/1.16_wp)
            ztmp2(1:i1)   = ztmp2(1:i1)/(pi*crhosno*cn0s)
            ztmp2(1:i1)   = SQRT(ztmp2(1:i1))
            ztmp3(1:i1)   = ztmp2(1:i1)**1.3125_wp
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            DO nl = 1,i1
               jl = idx1(nl)
               zclcstar = zclcpre(jl)
               zcfac4c  = 0.78_wp*ztmp2(nl)+232.19_wp*ztmp1(nl)*ztmp3(nl)
               zzeps    = -zxsec*zsfl(jl)*zclcpre_inv(jl)
               zzeps    = MAX(zzeps,zcoeff(jl)*zcfac4c*pmdry(jl,jk))
               zsub(jl) = -(zzeps/pmdry(jl,jk))*pdtime*zclcstar
               zsub(jl) = MIN(zsub(jl),MAX(zxsec*(zqsi(jl)-pqm1(jl,jk)),0.0_wp))
               zsub(jl) = MAX(zsub(jl),0.0_wp)
               zsub(jl) = MIN(zsub(jl),zsfl(jl)/pmdry(jl,jk)*pdtime)
            END DO
          END IF
          !
          !    end if(zsfl(jl).GT.ctqmin)
          !
          !       3.3   Evaporation of rain (Rotstayn, 1997)
          !
          !    old if(zrfl(jl).gt.cqtmin)
          !
          IF (i2.GT.0) THEN
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            DO nl = 1,i2
               jl = idx2(nl)
               zesw      = uaw(jl)*zpapm1_inv(jl)
               zesat     = uaw(jl)/rd
               zesw      = MIN(zesw,0.5_wp)
               zqsw      = zesw/(1._wp-vtmpc1*zesw)
               zsusatw   = MIN(pqm1(jl,jk)/zqsw-1.0_wp,0.0_wp)
               zdv       = 2.21_wp*zpapm1_inv(jl)
               zptm1_inv = 1._wp/ptm1(jl,jk)
               zast      = alv*(alv*zptm1_inv/rv-1.0_wp)*zptm1_inv/0.024_wp
               zbst      = ptm1(jl,jk)/(zdv*zesat)
               ztmp1(nl) = zast+zbst
               ztmp2(nl) = zrfl(jl)*zclcpre_inv(jl)
               ztmp3(nl) = zqsw
               ztmp4(nl) = zsusatw
            END DO
            ztmp2(1:i2) = ztmp2(1:i2)**0.61_wp
!IBM* ASSERT(NODEPS)
            DO nl = 1,i2
               jl = idx2(nl)
               zqsw     = ztmp3(nl)
               zsusatw  = ztmp4(nl)
               zclcstar = zclcpre(jl)
               zzepr    = 870._wp*zsusatw*ztmp2(nl)*zqrho_sqrt(jl)/SQRT(1.3_wp)
               zzepr    = zzepr/ztmp1(nl)
               zzepr    = MAX(-zxsec*zrfl(jl)*zclcpre_inv(jl),zzepr*pmdry(jl,jk))
               zevp(jl) = -(zzepr/pmdry(jl,jk))*pdtime*zclcstar
               zevp(jl) = MIN(zevp(jl),MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_wp))
               zevp(jl) = MAX(zevp(jl),0.0_wp)
               zevp(jl) = MIN(zevp(jl),zrfl(jl)/pmdry(jl,jk)*pdtime)
            END DO
          END IF
          !
          !    end if(zrfl(jl).gt.cqtmin)
          !
        END IF ! nclcpre.GT.0

      ELSE

        DO jl = 1,kproma
           zimlt(jl)  = 0.0_wp
           zsmlt(jl)  = 0.0_wp
        END DO

      END IF ! jk.GT.1
      !
#ifdef _PROFILE
      CALL trace_stop ('cloud_loop_3', 13)
      CALL trace_start ('cloud_loop_4', 14)
#endif
      !
      !     --------------------------------------------------------------------------
      !       4.    Sedimentation of cloud ice from grid-mean values.
      !     --------------------------------------------------------------------------
      !
      !             Updating the tendency 'pxite' to include sedimentation.
      !             At jk=klev, the sedimentation sink is balanced by
      !             precipitation at the surface (through 'zzdrs', see 7.3).
      !             Finally: In-cloud cloud water/ice.
      !

      DO 401 jl = 1,kproma
        zxip1         = pxim1(jl,jk)+(pxite(jl,jk)+pxteci(jl,jk))*pdtime-zimlt(jl)
      !  zxip1         = pxim1(jl,jk)+pxite(jl,jk)*pdtime-zimlt(jl)
        zxip1         = MAX(zxip1,EPSILON(1._wp))
        ztmp1(jl)     = zxip1
        ztmp2(jl)     = prho(jl,jk)*zxip1
401   END DO

      ztmp2(1:kproma) = ztmp2(1:kproma)**0.16_wp

!IBM* NOVECTOR
      DO 402 jl = 1,kproma
        zxifall       = cvtfall*ztmp2(jl)
        zal1          = -zxifall*pdtime/pdz(jl,jk)
        ztmp3(jl)     = zal1
402   END DO

      ztmp3(1:kproma) = EXP(ztmp3(1:kproma))

!IBM* NOVECTOR
      DO 410 jl = 1,kproma
        zxitop        = zxiflux(jl)
        zxip1         = ztmp1(jl)
        zxifall       = cvtfall*ztmp2(jl)
        zal2          = zxitop/(prho(jl,jk)*zxifall)
        zxised(jl)    = MAX(0.0_wp,zxip1*ztmp3(jl)+zal2*(1._wp-ztmp3(jl)))
        zqsed(jl)     = zxised(jl)-zxip1
        zxibot        = MAX(0.0_wp,zxitop-zqsed(jl)*pmdry(jl,jk)/pdtime)
        zqsed(jl)     = (zxitop-zxibot)/pmdry(jl,jk)*pdtime
        zxised(jl)    = zxip1+zqsed(jl)
        zxiflux(jl)   = zxibot
410   END DO

      DO 411 jl = 1,kproma
        zmrateps(jl,jk)=zmrateps(jl,jk)+(ztmp1(jl)-zxised(jl))
411   END DO

      locnt = 0
      nlocnt = 0
      DO 420 jl = 1,kproma
        zclcaux(jl) = paclc(jl,jk)
        IF (zclcaux(jl) .GT. 0.0_wp) THEN    ! locc=T
           locnt = locnt + 1
           loidx(locnt) = jl
        ELSE                                 ! locc=F
           nlocnt = nlocnt + 1
           nloidx(nlocnt) = jl
        END IF
        zesw          = uaw(jl)*zpapm1_inv(jl)
        zesw          = MIN(zesw,0.5_wp)
        zqsw          = zesw/(1._wp-vtmpc1*zesw)
      !  zsupsatw(jl)  = MAX(pqm1(jl,jk)/zqsw-1.0_wp,0.0_wp)
420   END DO

      !
      ! definition of lo2 (new zlo2(jl))
      !
      DO 421 jl = 1,kproma

      !       IF ((ptm1(jl,jk).LT.cthomi).OR.
      !         (ptm1(jl,jk).LT.tmelt.AND.zxised(jl).GT.csecfrli                     &
      !                                         .AND.zsupsatw(jl).LT.zeps)) THEN
      !         cond1(jl) = 1
      !       ELSE
      !         cond1(jl) = 0
      !       END IF

        zlo2(jl)  = FSEL(ptm1(jl,jk)+ptte(jl,jk)*pdtime-tmelt, 0._wp, 1._wp)
        zlo2(jl)  = FSEL(csecfrl-zxised(jl), 0._wp, zlo2(jl))
      !!$  zlo2(jl)  = FSEL(zsupsatw(jl)-zeps, 0._wp, zlo2(jl))
        zlo2(jl)  = FSEL(ptm1(jl,jk)+ptte(jl,jk)*pdtime-cthomi, zlo2(jl), 1._wp)
        cond1(jl) = INT(zlo2(jl))
        zlo2(jl)  = zlo2(jl)-0.5_wp ! zlo2 >= 0  <==> cond1 = 1
421   END DO
      !
      !             In-cloud water/ice calculated from respective grid-means,
      !             partial cloud cover, advective/diffusive tendencies,
      !             detrained cloud water/ice and ice sedimentation.
      !             In-cloud values are required for cloud microphysics.
      !
!IBM* ASSERT(NODEPS)
      DO 430 nl = 1,nlocnt
        jl = nloidx(nl)
        zxib(jl)     = 0.0_wp
        zxlb(jl)     = 0.0_wp
        zclcauxi(jl) = 0.0_wp
        zxip1        = pxim1(jl,jk)+(pxite(jl,jk)+pxteci(jl,jk))*pdtime+zqsed(jl)    &
                                                                      -zimlt(jl)
        zxlp1        = pxlm1(jl,jk)+(pxlte(jl,jk)+pxtecl(jl,jk))*pdtime+zimlt(jl)
        zxievap(jl)  = MAX(0.0_wp,zxip1)
        zxlevap(jl)  = MAX(0.0_wp,zxlp1)
430   END DO

!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO 431 nl = 1,locnt
        jl = loidx(nl)
        zclcauxi(jl) = 1._wp/zclcaux(jl)
        zxip1        = pxim1(jl,jk)+(pxite(jl,jk)+pxteci(jl,jk))*pdtime+zqsed(jl)    &
                                                                      -zimlt(jl)
        zxlp1        = pxlm1(jl,jk)+(pxlte(jl,jk)+pxtecl(jl,jk))*pdtime+zimlt(jl)
        zxib(jl)     = zxip1*zclcauxi(jl)
        zxlb(jl)     = zxlp1*zclcauxi(jl)
      !  zxievap(jl)  = (1.0_wp-zclcaux(jl))*MAX(0.0_wp,zxip1)
      !  zxlevap(jl)  = (1.0_wp-zclcaux(jl))*MAX(0.0_wp,zxlp1)
431   END DO
      !
      !     --------------------------------------------------------------------------
      !       5.    Condensation/deposition and evaporation/sublimation
      !     --------------------------------------------------------------------------
      !
      !             zlc       =  L_{v/s} / c_p
      !             zlcdqsdt  = L dq_sat / c_p dT
      !             zdqsdt    = dq_sat / dT
      !
!IBM* NOVECTOR
      DO 500 jl = 1,kproma
        zlc         = FSEL(zlo2(jl),zlsdcp(jl),zlvdcp(jl))
        zua         = FSEL(zlo2(jl),ua(jl),uaw(jl))
        zdua        = FSEL(zlo2(jl),dua(jl),duaw(jl))
        zqsm1(jl)   = zua*zpapm1_inv(jl)
        zqsm1(jl)   = MIN(zqsm1(jl),0.5_wp)
        zcor        = 1._wp/(1._wp-vtmpc1*zqsm1(jl))
        zqsm1(jl)   = zqsm1(jl)*zcor
        zdqsdt      = zpapm1_inv(jl)*zcor**2*zdua
        zlcdqsdt    = zlc*zdqsdt
        zdtdt(jl)   = ptte(jl,jk)*pdtime-zlvdcp(jl)*(zevp(jl)+zxlevap(jl))           &
                    - zlsdcp(jl)*(zsub(jl)+zxievap(jl))                              &
                    - (zlsdcp(jl)-zlvdcp(jl))*(zsmlt(jl)+zimlt(jl))
        zstar1(jl)  = zclcaux(jl)*zlc*pqte(jl,jk)*pdtime
        zdqsat1(jl) = zdqsdt/(1._wp+zclcaux(jl)*zlcdqsdt)
        zqvdt(jl)   = pqte(jl,jk)*pdtime+zevp(jl)+zsub(jl)+zxievap(jl)+zxlevap(jl)
        zqp1(jl)    = MAX(pqm1(jl,jk)+zqvdt(jl),0.0_wp)
        ztp1(jl)    = ptm1(jl,jk)+zdtdt(jl)
        zgenti(jl)  = 0.0_wp
        zgentl(jl)  = 0.0_wp
        zxib(jl)    = MAX(zxib(jl),0.0_wp)
        zxlb(jl)    = MAX(zxlb(jl),0.0_wp)
      !
      !       Diagnostics: relative humidity
      !
        zrelhum        = pqm1(jl,jk)/zqsm1(jl)
        zrelhum        = MAX(MIN(zrelhum,1._wp),0._wp)
        prelhum(jl,jk) = zrelhum
500   END DO
      !
      DO 520 jl = 1,kproma
        zdtdtstar = zdtdt(jl)+zstar1(jl)
        zdqsat    = zdtdtstar*zdqsat1(jl)
        zxilb     = zxib(jl)+zxlb(jl)
        zqcdif    = (pqte(jl,jk)*pdtime-zdqsat)*zclcaux(jl)
        zqcdif    = MAX(zqcdif,-zxilb*zclcaux(jl))
        zqcdif    = MIN(zqcdif,zqsec*zqp1(jl))
        ztmp1(jl) = zqcdif
520   END DO

      i1 = 0
      DO 530 jl = 1,kproma
        zqcdif = ztmp1(jl)
        IF (zqcdif .LT. 0.0_wp) THEN                 ! cloud dissipation
          zxilb    = zxib(jl)+zxlb(jl)
          zifrac   = zxib(jl)/MAX(zepsec,zxilb)
          zifrac   = MAX(MIN(zifrac,1.0_wp),0.0_wp)
          zdep(jl) = zqcdif*zifrac
          zcnd(jl) = zqcdif*(1.0_wp-zifrac)
        ELSE                                         ! cloud generation
        ! deposition (lo2 = .TRUE.) or condensation (lo2 = .FALSE.)
          zdep(jl) = FSEL(zlo2(jl),zqcdif,0.0_wp)
          zcnd(jl) = FSEL(zlo2(jl),0.0_wp,zqcdif)
        END IF
530   END DO
      !
      !       5.4 Checking for supersaturation of whole grid-box
      !
      DO jl = 1,kproma
        ztp1tmp(jl) = ztp1(jl) + zlvdcp(jl)*zcnd(jl) + zlsdcp(jl)*zdep(jl)
        zqp1tmp(jl) = zqp1(jl) -            zcnd(jl) -            zdep(jl)
        zxip1       = MAX(pxim1(jl,jk)+(pxite(jl,jk)+pxteci(jl,jk))*pdtime            &
                      +zqsed(jl)-zimlt(jl)-zxievap(jl)+zgenti(jl)+zdep(jl),0.0_wp)
        ztmp1(jl)   = zxip1
      END DO

      CALL lookup_ubc(kproma,ztp1tmp(1),ub(1))
      CALL prepare_ua_index_spline('cloud (2)',kproma,ztp1tmp(1),idx1(1),za(1)       &
                                               ,ztmp1(1),nphase,zlo2(1),cond1(1))
      CALL lookup_ua_eor_uaw_spline(kproma,idx1(1),za(1),nphase,cond1(1),ua(1),dua(1))
      zpapp1i(1:kproma) = 1._wp/papm1(1:kproma,jk)

!IBM* NOVECTOR
      DO 540 jl = 1,kproma
        zes         = ua(jl)*zpapp1i(jl)
        zes         = MIN(zes,0.5_wp)
        zcor        = 1._wp/(1._wp-vtmpc1*zes)
        zqsp1tmp    = zes*zcor
        zoversat    = zqsp1tmp*0.01_wp
        zdqsdt      = zpapp1i(jl)*zcor**2*dua(jl)
        zlc         = FSEL(zlo2(jl),zlsdcp(jl),zlvdcp(jl))
        zlcdqsdt    = FSEL(zes-0.4_wp,zqsp1tmp*zcor*ub(jl),zlc*zdqsdt)
        zqcon       = 1._wp/(1._wp+zlcdqsdt)
        zcor        = MAX((zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon,0._wp)
        zupdate     = FSEL(zlo2(jl),zdep(jl),zcnd(jl)) + zcor
        zdep(jl)    = FSEL(zlo2(jl),zupdate,zdep(jl)) ! ice cloud
        zcnd(jl)    = FSEL(zlo2(jl),zcnd(jl),zupdate) ! water cloud
        ztmp2(jl)   = zqsp1tmp
540   END DO

      ! mpuetz: ztmp2 holds inverse of zqsp1tmp
      ztmp2(1:kproma) = 1._wp/ztmp2(1:kproma)
      !
      !       5.5 Change of in-cloud water due to deposition/sublimation and
      !           condensation/evaporation (input for cloud microphysics)
      !
!IBM* ASSERT(NODEPS)
      DO 551 nl = 1,locnt
        jl = loidx(nl)
        zxib(jl) = MAX(zxib(jl)+zdep(jl)*zclcauxi(jl),0.0_wp)
        zxlb(jl) = MAX(zxlb(jl)+zcnd(jl)*zclcauxi(jl),0.0_wp)
 551  END DO

!IBM* ASSERT(NODEPS)
      DO 552 nl = 1,nlocnt
        jl = nloidx(nl)
        zdepos           = MAX(zdep(jl)+zgenti(jl),0.0_wp)
        zcond            = MAX(zcnd(jl)+zgentl(jl),0.0_wp)
        IF (zdepos>0.0_wp .OR. zcond>0.0_wp) THEN
          zclcaux(jl)   = 1.0_wp
          zclcauxi(jl)  = 1.0_wp
          zxib(jl)      = zdepos*zclcauxi(jl)
          zxlb(jl)      = zcond *zclcauxi(jl)
        END IF
 552  END DO

      DO 553 jl = 1,kproma
        ztp1tmp(jl) = ztp1(jl)+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
553   END DO
      !
      !     --------------------------------------------------------------------------
      !       6.    Freezing of cloud water
      !
      !       6.1   Freezing of cloud water for T < 238 K
      !
      DO 610 jl = 1,kproma
        zlo = cthomi - ztp1tmp(jl)
        ! mpuetz: using compute & select is faster than branching
        ! mpuetz: initially zfrl() is zero ?
        zfrl(jl)  = FSEL(zlo,zxlb(jl)*zclcaux(jl),0.0_wp)
        zxib(jl)  = FSEL(zlo,zxib(jl)+zxlb(jl),zxib(jl))
        zxlb(jl)  = FSEL(zlo,0.0_wp,zxlb(jl))
 610  END DO
      !
      !       6.2   Freezing of cloud water between 238 and 273 K
      !
      locnt = 0
      DO 620 jl = 1,kproma
        ! triple floating point compare + high predictability
        ! -> branched logic works best in SMT mode
        lo = ztp1tmp(jl).GT.cthomi.AND. &
             ztp1tmp(jl).LT.tmelt.AND.  &
             zxlb(jl).GT.0._wp

        IF (lo) THEN
          locnt = locnt + 1
          loidx(locnt) = jl
        END IF
620   END DO

      IF (locnt > 0) THEN
!IBM* ASSERT(NODEPS)
        DO 621 nl = 1,locnt
          jl = loidx(nl)
          ztmp1(nl) = 0.66_wp*(tmelt-ztp1tmp(jl))
621     END DO

      ztmp1(1:locnt) = EXP(ztmp1(1:locnt))

!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO 622 nl = 1,locnt
        jl = loidx(nl)
        zfrho    = prho(jl,jk)/(rhoh2o*pacdnc(jl,jk))
        zfrl(jl) = 100._wp*(ztmp1(nl)-1._wp)*zfrho
#if defined (__PGI)
        zfrl(jl) = zxlb(jl)*(1._wp-1._wp/(1._wp+zfrl(jl)*pdtime*zxlb(jl)))
#else
        zfrl(jl) = zxlb(jl)*(1._wp-SWDIV_NOCHK(1._wp,(1._wp+zfrl(jl)*pdtime*zxlb(jl))))
#endif
        ztmp1(nl)= 0.75_wp*zxlb(jl)*zfrho/pi
622   END DO

      ztmp1(1:locnt) = ztmp1(1:locnt)**(1._wp/3._wp)

!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO 623 nl = 1,locnt
        jl = loidx(nl)
        zradl    = ztmp1(nl)
        zval     = 4._wp*pi*zradl*pacdnc(jl,jk)*2.e5_wp*(tmelt-3._wp-ztp1tmp(jl))
        zf1      = zval/prho(jl,jk)
        zf1      = MAX(0.0_wp,zf1)
        zfrl(jl) = zfrl(jl)+pdtime*1.4e-20_wp*zf1
        zfrl(jl) = MAX(0.0_wp,MIN(zfrl(jl),zxlb(jl)))
        zxlb(jl) = zxlb(jl)-zfrl(jl)
        zxib(jl) = zxib(jl)+zfrl(jl)
        zfrl(jl) = zfrl(jl)*zclcaux(jl)
623   END DO
    END IF
    !
    !       6.3 Revised Bergeron-Findeisen process
    !
    !!$        DO 630 jl = 1,kproma
    !!$        locc        = zclcaux(jl) .GT. 0.0_wp
    !!$           IF (locc .AND. zdep(jl)>0._wp .AND. zxlb(jl)>0._wp .AND.           &
    !!$                          zxib(jl)>csecfrl .AND. zsupsatw(jl)<zeps) THEN
    !!$              zzevp        = zxlb(jl)*zclcaux(jl)/pdtime
    !!$              pxlte(jl,jk) = pxlte(jl,jk)-zzevp
    !!$              pxite(jl,jk) = pxite(jl,jk)+zzevp
    !!$              ptte(jl,jk)  = ptte(jl,jk)+(zlsdcp(jl)-zlvdcp(jl))*zzevp
    !!$              zxib(jl)     = zxib(jl)+zxlb(jl)
    !!$              zxlb(jl)     = 0.0_wp
    !!$           END IF
    !!$  630     END DO
    !
#ifdef _PROFILE
    CALL trace_stop ('cloud_loop_4', 14)
    CALL trace_start ('cloud_loop_7', 17)
#endif
    !
    !     ----------------------------------------------------------------------------
    !       7.  Cloud physics and precipitation fluxes at the surface
    !     ----------------------------------------------------------------------------
    !
    zxlb(1:kproma) = MAX(zxlb(1:kproma),1.e-20_wp)
    zxib(1:kproma) = MAX(zxib(1:kproma),1.e-20_wp)
    zmlwc(1:kproma,jk)=zxlb(1:kproma)
    zmiwc(1:kproma,jk)=zxib(1:kproma)

!IBM* NOVECTOR
    DO 701 jl = 1,kproma
      zauloc(jl) = cauloc*pdz(jl,jk)/5000._wp
      zauloc(jl) = MAX(MIN(zauloc(jl),clmax),clmin)
701 END DO

    zxrp1(1:kproma) = 0.0_wp
    zxsp1(1:kproma) = 0.0_wp
    !
!IBM* ASSERT(NODEPS)
    DO 711 nl = 1,nclcpre
      jl = jjclcpre(nl)
      ztmp1(nl) = (zrfl(jl)*zclcpre_inv(jl))/(12.45_wp*zqrho_sqrt(jl))
      ztmp2(nl) = zsfl(jl)*zclcpre_inv(jl)/cvtfall
711 END DO

    ztmp1(1:nclcpre) = ztmp1(1:nclcpre)**(8._wp/9._wp)
    ztmp2(1:nclcpre) = ztmp2(1:nclcpre)**(1._wp/1.16_wp)

!IBM* ASSERT(NODEPS)
    DO 712 nl = 1,nclcpre
      jl = jjclcpre(nl)
      zxrp1(jl) = ztmp1(nl)
      zxsp1(jl) = ztmp2(nl)
712 END DO
    !
    !       7.1   Warm clouds: Coalescence processes after Beheng (1994):
    !             Autoconversion of cloud droplets and collection of cloud
    !             droplets by falling rain. Accretion of cloud droplets by
    !             falling snow (zsacl) is calculated under 7.2
    !
    i1 = 0
    i2 = 0
    locnt = 0
    DO 720 jl = 1,kproma
      zsacl(jl) = 0.0_wp
      zrpr(jl)  = 0.0_wp
      zspr(jl)  = 0.0_wp
      IF (zclcaux(jl).GT.0.0_wp .AND. (zxlb(jl)>cqtmin.OR.zxib(jl)>cqtmin)) THEN
        locnt = locnt + 1
        loidx(locnt) = jl
      END IF
720 END DO

    IF (locnt.GT.0) THEN
      zexm1    = 4.7_wp-1.0_wp
      zexp     = -1._wp/zexm1

!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO nl = 1,locnt
         jl = loidx(nl)
         ztmp1(nl) = (ccraut*1.2e27_wp)/prho(jl,jk)
         ztmp2(nl) = pacdnc(jl,jk)*1.e-6_wp
         ztmp3(nl) = prho(jl,jk)  *1.e-3_wp
         ztmp4(nl) = zxlb(jl)
      END DO
      ztmp2(1:locnt) = ztmp2(1:locnt)**(-3.3_wp)
      ztmp3(1:locnt) = ztmp3(1:locnt)**4.7_wp
      ztmp4(1:locnt) = ztmp4(1:locnt)**zexm1
      DO nl = 1,locnt
         zraut     = ztmp1(nl)*ztmp2(nl)*ztmp3(nl)
         ztmp4(nl) =  1._wp+zraut*pdtime*zexm1*ztmp4(nl)
      END DO
      ztmp4(1:locnt) = ztmp4(1:locnt)**zexp

!IBM* ASSERT(NODEPS)
      DO nl = 1,locnt
         jl = loidx(nl)
         zraut     = zxlb(jl)*(1._wp-ztmp4(nl))
         ztmp1(nl) = -ccracl*zxrp1(jl)*pdtime
         ztmp2(nl) = -ccracl*zauloc(jl)*prho(jl,jk)*zraut*pdtime
         ztmp3(nl) = zxib(jl)*prho(jl,jk)*1000._wp
         ztmp4(nl) = zraut
      END DO
      ztmp1(1:locnt) = EXP(ztmp1(1:locnt))
      ztmp2(1:locnt) = EXP(ztmp2(1:locnt))
      ztmp3(1:locnt) = ztmp3(1:locnt)**0.216_wp

!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      DO nl = 1,locnt
         jl = loidx(nl)
         zraut = ztmp4(nl)
         zxlb(jl) = zxlb(jl)-zraut
         zrac1    = zxlb(jl)*(1._wp-ztmp1(nl))
         zxlb(jl) = zxlb(jl)-zrac1
         zrac2    = zxlb(jl)*(1._wp-ztmp2(nl))
         zxlb(jl) = zxlb(jl)-zrac2
         zclcstar = MIN(zclcaux(jl),zclcpre(jl))
         zrpr(jl) = zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1
         ! zrpr is initialized to zero
         zmratepr(jl,jk)=zraut+zrac1+zrac2
      END DO
      !
      !       7.2  Cold clouds:
      !            Conversion of cloud ice to snow after Levkov et al. 1992:
      !            Aggregation of ice crystals to snow and accretion of ice
      !            by falling snow.
      !            Accrection of cloud droplets by falling snow.
      !            Effective radius of ice crystals after Moss (1995)
      !
!IBM* ASSERT(NODEPS)
      DO nl = 1,locnt
         jl = loidx(nl)
         zrieff    = 83.8_wp*ztmp3(nl)
         zrieff    = MIN(MAX(zrieff,ceffmin),ceffmax)
         ztmp1(nl) = 5113188._wp+2809._wp*zrieff*zrieff*zrieff
         ztmp2(nl) = zqrho(jl)
         ztmp3(nl) = 0.025_wp*(ztp1tmp(jl)-tmelt)
      END DO
      ztmp1(1:locnt) = SQRT(ztmp1(1:locnt))
      ztmp1(1:locnt) = ztmp1(1:locnt)-2261._wp
      ztmp1(1:locnt) = LOG10(ztmp1(1:locnt))
      ztmp2(1:locnt) = ztmp2(1:locnt)**0.33_wp
      ztmp3(1:locnt) = EXP(ztmp3(1:locnt))

!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      DO 721 nl = 1,locnt
         jl = loidx(nl)
         zc1       = 17.5_wp*prho(jl,jk)/crhoi*ztmp2(nl)
         zdt2      = -6._wp/zc1*(ztmp1(nl)/3._wp-2._wp)
         zsaut     = ccsaut/zdt2
         zsaut     = zxib(jl)*(1._wp-1._wp/(1._wp+zsaut*pdtime*zxib(jl)))
         zxib(jl)  = zxib(jl)-zsaut
         zxsp2(nl) = zauloc(jl)*prho(jl,jk)*zsaut
         ztmp1(nl) = zsaut
721   END DO

!IBM* NOVECTOR
      DO 722 nl = 1,locnt
         jl = loidx(nl)
         zsaut     = ztmp1(nl)
         zcolleffi = ztmp3(nl)
         zclcstar = MIN(zclcaux(jl),zclcpre(jl))
         zsaci1    = 0.0_wp
         zsaci2    = 0.0_wp
         zsacl1    = 0.0_wp
         zsacl2    = 0.0_wp
         IF (zxsp1(jl) .GT. cqtmin) THEN
           zlamsm    = (zxsp1(jl)/(pi*crhosno*cn0s))**0.8125_wp
           zsaci1    = pi*cn0s*3.078_wp*zlamsm*zqrho_sqrt(jl)
           zsacl1    = zxlb(jl)*(1._wp-EXP(-zsaci1*ccsacl*pdtime))
           zxlb(jl)  = zxlb(jl)-zsacl1
           zsacl1    = zclcstar*zsacl1
           zsaci1    = zsaci1*zcolleffi*pdtime
           zsaci1    = zxib(jl)*(1._wp-EXP(-zsaci1))
           zxib(jl)  = zxib(jl)-zsaci1
         END IF
         IF (zxsp2(nl) .GT. cqtmin) THEN
           zlamsm    = (zxsp2(nl)/(pi*crhosno*cn0s))**0.8125_wp
           zsaci2    = pi*cn0s*3.078_wp*zlamsm*zqrho_sqrt(jl)
           zsacl2    = zxlb(jl)*(1._wp-EXP(-zsaci2*ccsacl*pdtime))
           zxlb(jl)  = zxlb(jl)-zsacl2
           zsacl2    = zclcaux(jl)*zsacl2
           zsaci2    = zsaci2*zcolleffi*pdtime
           zsaci2    = zxib(jl)*(1._wp-EXP(-zsaci2))
           zxib(jl)  = zxib(jl)-zsaci2
         END IF
         zsacl(jl)    = zsacl1+zsacl2
         zspr(jl)     = zclcaux(jl)*(zsaut+zsaci2) + zclcstar*zsaci1
         ! zspr is initialized to zero

         IF(zclcstar>zepsec .AND. zclcaux(jl)>zepsec) THEN
           zmsnowacl(jl,jk)=zsacl1/zclcstar+zsacl2/zclcaux(jl)
         ELSE IF (zclcstar>zepsec .AND. zclcaux(jl)<=zepsec) THEN
           zmsnowacl(jl,jk)=zsacl1/zclcstar
         ELSE IF (zclcstar<=zepsec .AND. zclcaux(jl)>zepsec) THEN
           zmsnowacl(jl,jk)=zsacl2/zclcaux(jl)
         ELSE
           zmsnowacl(jl,jk)=0._wp
         END IF

         zmrateps(jl,jk)=zmrateps(jl,jk)+zsaut+zsaci1+zsaci2

722   END DO

    END IF ! locnt > 0

    !
    !       7.3 Updating precipitation fluxes. In the lowest layer (klev),
    !           the sedimentation sink of cloud ice is balanced
    !           by precipitation at the surface (through 'zzdrs').
    !           Fraction of precipitating clouds (zclcpre) used for the
    !           calculation of evaporation/sublimation of rain/snow in
    !           the next layer
    !
    zcnt = 0._wp
    nclcpre = 0

    IF (jk.EQ.klev) THEN

!IBM* NOVECTOR
      DO jl = 1,kproma
         zzdrr       =  zrpr(jl)           *pmdry(jl,jk)/pdtime
         zzdrs       = (zspr(jl)+zsacl(jl))*pmdry(jl,jk)/pdtime
         zzdrs       = zzdrs+zxiflux(jl)
         zcons       = (pmdry(jl,jk)/pdtime)/(zlsdcp(jl)-zlvdcp(jl))
         zsnmlt      = MIN(zxsec*zzdrs,zcons*MAX(0._wp,(ztp1tmp(jl)-tmelt)))
         zzdrr       = zzdrr+zsnmlt
         zzdrs       = zzdrs-zsnmlt
         zsmlt(jl)   = zsmlt(jl)+zsnmlt/pmdry(jl,jk)*pdtime
         zpretot     = zrfl(jl)+zsfl(jl)
         zpredel     = zzdrr+zzdrs
         zclcpre(jl) = FSEL(zpredel-zpretot,zclcaux(jl),zclcpre(jl))
         zpresum     = zpretot+zpredel
!#ifdef FAST_AND_DIRTY
!           ! This may trigger a divide by zero. It doesn't harm, because the
!           ! result is adjusted in the next following FSEL/MERGE command.
         zclcpre1       = (zclcaux(jl)*zpredel + zclcpre(jl)*zpretot)/zpresum
!#else
!           ! This ifdef branch does contain clean code for checking
!           ! the correctness of the complete code.
!           IF (zpresum < TINY(zpresum)) THEN
!             zclcpre1     = 0.0_wp
!           ELSE
!             zclcpre1     = (zclcaux(jl)*zpredel + zclcpre(jl)*zpretot)/zpresum
!           ENDIF
!#endif
         zclcpre1       = MAX(zclcpre(jl),zclcpre1)
         zclcpre1       = MAX(0._wp,zclcpre1)
         zclcpre1       = MIN(1._wp,zclcpre1)

         zclcpre(jl)    = FSEL(cqtmin-zpresum,0._wp,zclcpre1)

         zcnt           = zcnt + FSEL(-zclcpre(jl),0.0_wp,1.0_wp)
         cond1(jl)      = INT(FSEL(-zclcpre(jl),0.0_wp,1.0_wp))
         !   Corrected by Junhua Zhang, Philip Stier (01/2004)
         IF (zclcpre(jl) > zepsec) THEN
            zfrain(jl,jk)=(zrfl(jl)+zzdrr)/zclcpre(jl)
            zfsnow(jl,jk)=(zsfl(jl)+zzdrs)/zclcpre(jl)
            zfevapr(jl,jk)=(zevp(jl)*pmdry(jl,jk)/pdtime)/zclcpre(jl)
            zfsubls(jl,jk)=(zsub(jl)*pmdry(jl,jk)/pdtime)/zclcpre(jl)
         ELSE
            zfrain(jl,jk) =0.0_wp
            zfsnow(jl,jk) =0.0_wp
            zfevapr(jl,jk)=0.0_wp
            zfsubls(jl,jk)=0.0_wp
         ENDIF

         zrfl(jl)    = zrfl(jl)+zzdrr-zevp(jl)*pmdry(jl,jk)/pdtime
         zsfl(jl)    = zsfl(jl)+zzdrs-zsub(jl)*pmdry(jl,jk)/pdtime
      END DO

    ELSE

!IBM* NOVECTOR
      DO jl = 1,kproma
         zzdrr          =  zrpr(jl)           *pmdry(jl,jk)/pdtime
         zzdrs          = (zspr(jl)+zsacl(jl))*pmdry(jl,jk)/pdtime
         zpretot        = zrfl(jl)+zsfl(jl)
         zpredel        = zzdrr+zzdrs
         zclcpre(jl)    = FSEL(zpredel-zpretot,zclcaux(jl),zclcpre(jl))
         zpresum        = zpretot+zpredel

#ifdef FAST_AND_DIRTY
         ! This may trigger a divide by zero. It doesn't harm, because the
         ! result is adjusted in the next following FSEL/MERGE command.
         zclcpre1       = (zclcaux(jl)*zpredel + zclcpre(jl)*zpretot)/zpresum
#else
         ! This ifdef branch does contain clean code for checking
         ! the correctness of the complete code.
         IF (zpresum < TINY(zpresum)) THEN
           zclcpre1     = 0.0_wp
         ELSE
           zclcpre1     = (zclcaux(jl)*zpredel + zclcpre(jl)*zpretot)/zpresum
         ENDIF
#endif
         zclcpre1       = MAX(zclcpre(jl),zclcpre1)
         zclcpre1       = MAX(0._wp,zclcpre1)
         zclcpre1       = MIN(1._wp,zclcpre1)
         zclcpre(jl)    = FSEL(cqtmin-zpresum,0._wp,zclcpre1)
         zcnt           = zcnt + FSEL(-zclcpre(jl),0.0_wp,1.0_wp)
         cond1(jl)      = INT(FSEL(-zclcpre(jl),0.0_wp,1.0_wp))
         !   Corrected by Junhua Zhang, Philip Stier (01/2004)
         IF (zclcpre(jl) > zepsec) THEN
           zfrain(jl,jk)=(zrfl(jl)+zzdrr)/zclcpre(jl)
           zfsnow(jl,jk)=(zsfl(jl)+zzdrs)/zclcpre(jl)
           zfevapr(jl,jk)=(zevp(jl)*pmdry(jl,jk)/pdtime)/zclcpre(jl)
           zfsubls(jl,jk)=(zsub(jl)*pmdry(jl,jk)/pdtime)/zclcpre(jl)
         ELSE
           zfrain(jl,jk) =0.0_wp
           zfsnow(jl,jk) =0.0_wp
           zfevapr(jl,jk)=0.0_wp
           zfsubls(jl,jk)=0.0_wp
         ENDIF
         zrfl(jl)       = zrfl(jl)+zzdrr-zevp(jl)*pmdry(jl,jk)/pdtime
         zsfl(jl)       = zsfl(jl)+zzdrs-zsub(jl)*pmdry(jl,jk)/pdtime
      END DO

    END IF

    IF (zcnt > 0._wp) THEN
      nclcpre = 1
      DO jl = 1,kproma
         jjclcpre(nclcpre) = jl
         nclcpre = nclcpre + cond1(jl)
      END DO
      nclcpre = nclcpre - 1
    END IF

#ifdef _PROFILE
    CALL trace_stop ('cloud_loop_7', 17)
    CALL trace_start ('cloud_loop_8', 18)
#endif
    !     ----------------------------------------------------------------------------
    !       8.    Updating tendencies of t, q, xl, xi and final cloud cover
    !     ----------------------------------------------------------------------------
    !

!IBM* NOVECTOR
    DO 820 jl = 1,kproma
    !
    !       8.3   Tendencies of thermodynamic variables
    !
    !       local tendencies due to cloud microphysics
    !
       zqvte  = (-zcnd(jl)-zgentl(jl)+zevp(jl)+zxlevap(jl)                           &
         &       -zdep(jl)-zgenti(jl)+zsub(jl)+zxievap(jl))/pdtime
       zxlte  = (zimlt(jl)-zfrl(jl)-zrpr(jl)-zsacl(jl)                               &
         &                     +zcnd(jl)+zgentl(jl)-zxlevap(jl))/pdtime
       zxite  = (zfrl(jl)-zspr(jl)+zdep(jl)+zgenti(jl)                               &
         &                     -zxievap(jl)-zimlt(jl)+zqsed(jl))/pdtime
       ztte   = (zlvdcp(jl)*(zcnd(jl)+zgentl(jl)-zevp(jl)-zxlevap(jl))               &
         &     + zlsdcp(jl)*(zdep(jl)+zgenti(jl)-zsub(jl)-zxievap(jl))               &
         &     +(zlsdcp(jl)-zlvdcp(jl))                                              &
         &     *(-zsmlt(jl)-zimlt(jl)+zfrl(jl)+zsacl(jl)))/pdtime

       pqte(jl,jk)   = pqte(jl,jk)  + zqvte
       pxlte(jl,jk)  = pxlte(jl,jk) + pxtecl(jl,jk) + zxlte
       pxite(jl,jk)  = pxite(jl,jk) + pxteci(jl,jk) + zxite
       ptte(jl,jk)   = ptte(jl,jk)  + ztte
       zcpten(jl,jk) = ztte                        ! diagnostics [K/s]
       zqten(jl,jk)  = zqvte + zxlte + zxite       ! diagnostics [1/s] ...< 0
820 END DO

!IBM* NOVECTOR
    DO 821 jl = 1,kproma

       zxlp1        = pxlm1(jl,jk) + pxlte(jl,jk)*pdtime
       zxip1        = pxim1(jl,jk) + pxite(jl,jk)*pdtime
    !
    !       8.4   Corrections: Avoid negative cloud water/ice
    !
       zxlold         = zxlp1
       zxiold         = zxip1
       zxlp1_d        = ccwmin - zxlp1
       zxip1_d        = ccwmin - zxip1

       zxlp1          = FSEL(-zxlp1_d,zxlp1,0._wp)
       zxip1          = FSEL(-zxip1_d,zxip1,0._wp)
       zdxlcor        = (zxlp1 - zxlold)/pdtime
       zdxicor        = (zxip1 - zxiold)/pdtime

       zxlp1_d        = MAX(zxlp1_d,0.0_wp)
       paclc(jl,jk)   = FSEL(-(zxlp1_d*zxip1_d),paclc(jl,jk),0._wp)

       pxlte(jl,jk)   = pxlte(jl,jk) + zdxlcor
       pxite(jl,jk)   = pxite(jl,jk) + zdxicor
       pqte(jl,jk)    = pqte(jl,jk) - zdxlcor - zdxicor
       ptte(jl,jk)    = ptte(jl,jk) + zlvdcp(jl)*zdxlcor + zlsdcp(jl)*zdxicor
       ! Here mulitply with the same specific heat as used in the definition
       ! of zlvdcp and zlsdcp ( =Lv/(cp or cv) and Ls/(cp or cv) ) in order
       ! to obtain the specific heating by cloud processes in [W/kg].
       zcpten(jl,jk)  = pcair(jl,jk)                                                 &
                              *(zcpten(jl,jk)+zlvdcp(jl)*zdxlcor+zlsdcp(jl)*zdxicor)
       !
821 END DO

#ifdef _PROFILE
    CALL trace_stop ('cloud_loop_8', 18)
#endif
    !
831 END DO    ! Vertical loop
!
#ifdef _PROFILE
    CALL trace_start ('cloud_loop_9', 19)
#endif
    !
    !     ----------------------------------------------------------------------------
    !       9.    Wet chemistry and in-cloud scavenging
    !     ----------------------------------------------------------------------------
    !     ----------------------------------------------------------------------------
    !       10.    Diagnostics
    !     ----------------------------------------------------------------------------
    !
    !       10.1   Accumulated precipitation at the surface
    !
    DO 911 jl    = 1,kproma
       prsfl(jl) = zrfl(jl)
       pssfl(jl) = zsfl(jl)
911 END DO

    !
    !       10.2   Total cloud cover
    !
    DO 921 jl    = 1,kproma
       zclcov(jl) = 1.0_wp-paclc(jl,1)
921 END DO
    !
    DO 923 jk      = 2,klev
!IBM* NOVECTOR
      DO 922 jl    = 1,kproma
         zclcov(jl) = zclcov(jl) * ((1._wp-MAX(paclc(jl,jk),paclc(jl,jk-1)))         &
                                 / (1._wp-MIN(paclc(jl,jk-1),zxsec)))
922   END DO
923 END DO

    DO 924 jl     = 1,kproma
       zclcov(jl)  = 1.0_wp-zclcov(jl)
       paclcov(jl) = zclcov(jl)
924 END DO
    !
    !      10.3 Derive the tendency increment induced by this routine
    !
    ptte_prc(1:kproma,:)   =  ptte(1:kproma,:)   -  ptte_prc(1:kproma,:)
    pqte_prc(1:kproma,:)   =  pqte(1:kproma,:)   -  pqte_prc(1:kproma,:)
    pxlte_prc(1:kproma,:)   = pxlte(1:kproma,:)   - pxlte_prc(1:kproma,:)
    pxite_prc(1:kproma,:)   = pxite(1:kproma,:)   - pxite_prc(1:kproma,:)

    !
    !      10.4 Vertical integrals of humidity, cloud water and cloud ice
    !
    DO 931 jl   = 1,kproma
       zxlvi(jl) = 0.0_wp
       zclten(jl) = 0.0_wp
       zqviten(jl) = 0.0_wp
931 END DO
    !
    DO 933 jk     = ktdia,klev
       DO 932 jl   = 1,kproma
          zxlvi(jl)  = zxlvi(jl)   + (pxlm1 (jl,jk)+pxlte_prc(jl,jk)*pdtime)*pmdry(jl,jk)
          zclten(jl) = zclten(jl)  + zcpten(jl,jk) *pmdry(jl,jk)     ! [W/m2]
          zqviten(jl)= zqviten(jl) + zqten (jl,jk) *pmdry(jl,jk)     ! [kg/m2s]
932    END DO
933 END DO

!
    DO 934 jl   = 1,kproma
       pch_concloud(jl) = pch_concloud(jl)+zclten(jl)-(alv*prsfl(jl)+als*pssfl(jl)) ! [W/m2]
       pcw_concloud(jl) = pcw_concloud(jl)+zqviten(jl)+prsfl(jl)+pssfl(jl)          ! [kg/m2s]
934 END DO

    ! compare liquid water path below and above convective cloud top
    DO 938 jl = 1,kproma
       zxlvitop(jl) = 0.0_wp
       klevtop = kctop(jl) - 1
       DO 936 jk = ktdia, klevtop
          zxlvitop(jl) = zxlvitop(jl)+pxlm1(jl,jk)*pmdry(jl,jk)
936    END DO
938 END DO

    ! modify ktype where appropriate (to be used in mo_cloud_optics)
    DO 940 jl = 1,kproma
       zxlvibot(jl) = zxlvi(jl) - zxlvitop(jl)
       IF (ktype(jl) .EQ. 2 .AND. zxlvibot(jl) .GT. clwprat * zxlvitop(jl)) THEN
          ktype(jl) = 4
       END IF
940 END DO

#ifdef _PROFILE
    CALL trace_stop ('cloud_loop_9', 19)
    CALL trace_stop ('cloud', 10)
#endif

  END SUBROUTINE cloud

END MODULE mo_cloud
