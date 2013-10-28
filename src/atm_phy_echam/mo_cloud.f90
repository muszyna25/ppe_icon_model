!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
#ifdef __xlC__
@PROCESS HOT
#else
#define SWDIV_NOCHK(a,b) ((a)/(b))
#define FSEL(a,b,c) MERGE(b,c,(a).GE.0._dp)
#endif

MODULE mo_cloud

  USE mo_kind,               ONLY : dp=>wp
  USE mo_math_constants,     ONLY : api=>pi
  USE mo_physical_constants, ONLY : cpd, cvd, cpv, cvv, g=>grav, rd, alv, als, rv   &
                                  , vtmpc1, rhoh2o, tmelt
  USE mo_convect_tables,     ONLY : prepare_ua_index, lookup_ua, lookup_uaw  &
                                  , lookup_ubc, lookup_ua_eor_uaw
  USE mo_echam_cloud_params, ONLY : cqtmin, cvtfall, crhosno, cn0s           &
                                  , cthomi, csecfrl, ncctop, cvarmin         &
                                  , cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs  &
                                  , rbetak, nbetax, tbetai, cauloc           &
                                  , clmax, clmin, jbmin, jbmax, lonacc       &
                                  , ccraut, ceffmin, ceffmax, crhoi, ccsaut  &
                                  , ccsacl, ccracl, cbeta_cs, ccwmin
#ifdef _PROFILE
  USE mo_profile,        ONLY : trace_start, trace_stop
#endif


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cloud

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

SUBROUTINE cloud (         kproma,   kbdim,    ktdia            &
                         , klev,     klevp1                     &
                         , pdelta_time, ptime_step_len          &
! - INPUT  2D .
                         , paphm1                               &
                         , papm1,    papp1                      &
                         , ptm1,     ptvm1,    pgeo,    pvervel &
                         , pacdnc                               &
                         , pqm1,     pxlm1,    pxim1            &
! - INPUT  1D .
                         , knvb                                 &
! - in and out, 2D
                         , pqtec,    pxtec                      &
                         , ptte                                 &
                         , pqte,     pxlte,     pxite           &
                         , paclc                                &
! - OUTPUT 1D .
                         , paclcov,  pqvi                       &
                         , pxlvi,    pxivi                      &
                         , prsfl,    pssfl                      &
! - OUTPUT 2D .
                         , prelhum                              &
                         , ptte_prc, pqte_prc                   &
                         , pxlte_prc,pxite_prc   )
!
!     *Cloud* computes large-scale water phase changes, precipitation,
!             cloud cover, and vertical integrals of specific humidity,
!             cloud liquid water content and cloud ice (diagnostics).
!
!     Subject.
!     --------
!
!          This routine computes the tendencies of the four prognostic
!          variables (temperature t, specific humidity q, cloud liquid
!          water xl, cloud ice xi) due to phase changes (condensation/
!          deposition, evaporation/sublimation of rain/snow falling
!          into the unsaturated part of the grid box, melting of snow,
!          melting/freezing of cloud ice/cloud water, sedimentation of
!          cloud ice, and precipitation formation in warm, cold and
!          mixed phase clouds.
!          The precipitation at the surface (rain and snow) is used in
!          later for computing the land surface hydrology in *surf*.
!          The cloud parameters (cloud cover, cloud liquid water and
!          cloud ice are used for the calculation of radiation at the
!          next timestep.
!          Attention:
!          In the current version the advective tendencies of skewness
!          and variance are set to zero.
!
!     INTERFACE.
!     ----------
!
!     *Call cloud*
!
!     Input arguments.
!     ----- ----------
!  - 2D
!  paphm1   : pressure at half levels                              (n-1)
!  paphp1   : pressure at half levels                              (n+1)
!  papm1    : pressure at full levels                              (n-1)
!  papp1    : pressure at full levels                              (n+1)
!  ptm1     : temperature                                          (n-1)
!  ptvm1    : virtual temperature                                  (n-1)
!  pgeo     : geopotential minus its surface value                 (n-1)
!  pvervel  : vertical velocity in pressure coordinate             (n)
!  pacdnc   : cloud droplet number concentration
!  pqm1     : specific humidity                                    (n-1)
!  pxlm1    : cloud liquid water                                   (n-1)
!  pxim1    : cloud ice                                            (n-1)
!  pxtm1    : tracer (aerosol etc) concentration                   (n-1)
!  - 1D
!  knvb     : provided by subroutine "cover"
!
!     Input, Output arguments.
!     ------------ ----------
!  - 2D
!  pqtec    :
!  pxtec    : detrained convective cloud liquid water or cloud ice (n)
!  ptte     : tendency of temperature accumulated
!  pqte     : tendency of specific humidity accumulated
!  pxlte    : tendency of cloud liquid water accumulated
!  pxite    : tendency of cloud ice
!  paclc    : cloud cover  (now diagnosed in cover)
!
!     Output arguments.
!     ------ ----------
!  - 1D
!  paclcov  : total cloud cover
!  pqvi     : vertically integrated spec. humidity
!  pxlvi    : vertically integrated cloud liquid water
!  pxivi    : vertically integrated cloud ice
!  prsfl    : surface rain flux
!  pssfl    : surface snow flux
!  - 2D
!  prelhum  : relative humidity
!  ptte_prc : tendency of temperature resulting from this subroutine
!  pqte_prc : tendency of specific humidity resulting from this subroutine
!  pxlte_prc: tendency of cloud liquid water resulting from this subroutine
!  pxite_prc: tendency of cloud ice resulting from this subroutine
!
!     Externals.
!     ----------
!  prepare_ua_index_spline, prepare_ua_index
!  lookup_ua_spline, lookup_uaw_spline, lookup_ua_eor_uaw_spline
!  lookup_ua, lookup_uaw, lookup_ua_eor_uaw, lookup_ubc
!  cloud_subm, set_vphysc_var
!
!     Method.
!     -------
!     see References
!
!     References.
!     ----------
!
!     Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
!     Levkov et al., 1992: Beitr. Phys. Atm. 35-58.          (ice phase)
!     Beheng, 1994: Atmos. Res. 193-206.                    (warm phase)
!     Lenderink et al., 1998; KNMI-REPORT NO. 98-13       (condensation)
!     Tompkins 2002, J. Atmos. Sci.                        (cloud cover)
!
!     Authors.
!     -------
!     M.Esch        MPI-Hamburg  1999
!     G.Lenderink   KNMI, de Bilt 1998
!     U.Lohmann     MPI-Hamburg  1995
!
!     Modifications.
!     --------------
!     E.Roeckner    MPI-Hamburg  2000
!     A.Tompkins    MPI-Hamburg  2000
!     U.Schlese     MPI-Hamburg  2003
!     P.Stier       MPI-Hamburg  2002-2004: Scavenging parameters
!
!
  INTEGER, INTENT(IN) :: kbdim, klevp1, klev, kproma, ktdia
  REAL(dp),INTENT(IN) :: pdelta_time, ptime_step_len
  REAL(dp),INTENT(IN) :: paphm1(kbdim,klevp1) ,pvervel(kbdim,klev)  &
                       , papm1(kbdim,klev)    ,pqm1(kbdim,klev)     &
                       , papp1(kbdim,klev)    ,ptm1(kbdim,klev)     &
                       , ptvm1(kbdim,klev)    ,pxlm1(kbdim,klev)    &
                       , pxim1(kbdim,klev)                          &
                       , pgeo(kbdim,klev)

  REAL(dp),INTENT(INOUT) :: pxtec(kbdim,klev)    ,pqtec(kbdim,klev)

  REAL(dp),INTENT(INOUT) :: pxlvi(kbdim)         ,pxivi(kbdim)         ! OUT
  REAL(dp),INTENT(INOUT) :: paclc(kbdim,klev)
  REAL(dp),INTENT(IN)    :: pacdnc(kbdim,klev)
  REAL(dp),INTENT(INOUT) :: prelhum(kbdim,klev)                        ! OUT
  REAL(dp),INTENT(INOUT) :: paclcov(kbdim)       , pqvi(kbdim)         ! OUT

  REAL(dp),INTENT(INOUT) :: ptte(kbdim,klev)     ,pqte(kbdim,klev)
  REAL(dp),INTENT(INOUT) :: pxlte(kbdim,klev)    ,pxite(kbdim,klev)
  REAL(dp),INTENT(INOUT) :: ptte_prc(kbdim,klev) ,pqte_prc(kbdim,klev)  ! OUT
  REAL(dp),INTENT(INOUT) :: pxlte_prc(kbdim,klev),pxite_prc(kbdim,klev) ! OUT
  REAL(dp),INTENT(INOUT) :: pssfl(kbdim)         ,prsfl(kbdim)          ! OUT

!
!   Temporary arrays
!
  REAL(dp):: zclcpre(kbdim)       ,zclcpre_inv(kbdim)                  &
           , zcnd(kbdim)          ,zdep(kbdim)        ,zdp(kbdim)      &
           , zevp(kbdim)          ,zxievap(kbdim)     ,zxlevap(kbdim)  &
           , zfrl(kbdim)          ,zimlt(kbdim)       ,zsmlt(kbdim)    &
           , zrpr(kbdim)          ,zspr(kbdim)        ,zsub(kbdim)     &
           , zxlte(kbdim)         ,zxite(kbdim)       ,zxiflux(kbdim)  &
           , zsacl(kbdim)         ,zdz(kbdim)         ,zqp1(kbdim)     &
           , zlsdcp(kbdim)        ,zlvdcp(kbdim)      ,zximlt(kbdim)   &
           , ztp1tmp(kbdim)       ,zqp1tmp(kbdim)     ,zxisub(kbdim)   &
           , zrfl(kbdim)          ,zsfl(kbdim)        ,ztp1(kbdim)     &
           , zxlb(kbdim)          ,zxib(kbdim)        ,zqrho(kbdim)    &
           , zqrho_sqrt(kbdim)    ,zpapm1_inv(kbdim)  ,zpapp1i(kbdim)  &
           , zclcov(kbdim)        ,zclcaux(kbdim)                      &
           , zqvi(kbdim)          ,zxlvi(kbdim)       ,zxivi(kbdim)    &
           , zclcauxi(kbdim)                                           &
           , zdqsat1(kbdim)                                            &
           , zxrp1(kbdim)         ,zxsp1(kbdim)       ,zxsp2(kbdim)   &
           , zgenti(kbdim)        ,zgentl(kbdim)                       &
           , zcoeff(kbdim)        ,zrhtest(kbdim)                      &
           , zgeoh(kbdim,klevp1)  ,zauloc(kbdim)      ,zqsi(kbdim)     &
           , ztmp1(kbdim)         ,ztmp2(kbdim)       ,ztmp3(kbdim)    &
           , ztmp4(kbdim)         ,zxised(kbdim)      ,zqvdt(kbdim)    &
           , zqsm1(kbdim)         ,ub(kbdim)                           &
           , zdtdt(kbdim)         ,zstar1(kbdim)      ,zlo2(kbdim)     &
           , ua(kbdim)            ,dua(kbdim)                          &
           , uaw(kbdim)           ,duaw(kbdim)        ,zsupsatw(kbdim)


  REAL(dp):: zrho(kbdim,klev)
!
  INTEGER,INTENT(IN) :: knvb(kbdim)
  INTEGER:: loidx(kbdim), nloidx(kbdim), jjclcpre(kbdim)
  INTEGER:: cond1(kbdim), cond2(kbdim), cond3(kbdim)
  INTEGER:: idx1(kbdim), idx2(kbdim), idx3(kbdim)

  INTEGER:: jb, nclcpre
  INTEGER:: jl, jk, nl, locnt, nlocnt, nphase, i1 , i2 , i3, i4
  LOGICAL:: lo, lo1

  REAL(dp):: zdqsat, zqcdif, zfrho                                     &
           , zifrac                                                    &
           , zdtime, zxiupd, zxlupd                                    &
           , zepsec, zxsec, zqsec, ztmst, zcons2, zrcp, zcons          &
           , ztdif, zsnmlt, zximelt, zclcstar, zdpg, zesi              &
           , zsusati, zb1, zb2, zcfac4c, zzeps                         &
           , zsubi, zesw, zesat, zqsw, zsusatw, zdv, zast, zbst        &
           , zzepr, zxip1, zxifall, zal1, zal2, zxim1evp               &
           , zxlm1evp, zxidt, zxldt, zxidtstar, zxldtstar, zlc         &
           , zqst1, zdqsdt, zlcdqsdt, zdtdtstar                        &
           , zxilb, zrelhum, zes                                       &
           , zcor, zqsp1tmp, zoversat, zqcon                           &
           , zdepos, zcond, zradl, zf1, zraut, zexm1, zexp, zrac1      &
           , zrac2, zrieff, zcolleffi, zc1, zdt2, zsaut                &
           , zsaci1, zsaci2, zsacl1, zsacl2, zlamsm, zzdrr             &
           , zzdrs, zpretot, zpredel, zpresum, zxlp1                   &
           , zxlold, zxiold, zdxicor, zdxlcor, zptm1_inv               &
           , zxlp1_d, zxip1_d,zupdate, zdefault, zlo, zcnt, zclcpre1   &
           , zval, zua, zdua, zeps
!
  ! save the tendencies accumulated before calling this routine

     ptte_prc(1:kproma,:)   =  ptte(1:kproma,:)
     pqte_prc(1:kproma,:)   =  pqte(1:kproma,:)
    pxlte_prc(1:kproma,:)   = pxlte(1:kproma,:)
    pxite_prc(1:kproma,:)   = pxite(1:kproma,:)
!
! Executable statements
!
#ifdef _PROFILE
  CALL trace_start ('cloud', 10)
#endif

!
!   Security parameters
!
  zepsec = 1.0e-12_dp
  zxsec  = 1.0_dp-zepsec
  zqsec  = 1.0_dp-cqtmin
  zeps   = EPSILON(1.0_dp)
!
!   Computational constants
!
  zdtime = REAL(pdelta_time,dp)
  ztmst  = REAL(ptime_step_len,dp)
  zcons2 = 1._dp/(ztmst*g)
!
!     ------------------------------------------------------------------
!
!       1.   Top boundary conditions, air density and geopotential
!            height at half levels
!
!       1.1   Set to zero precipitation fluxes etc.
!
  nclcpre = 0

  DO 111 jl = 1,kproma
     zclcpre(jl)   = 0.0_dp
     zxiflux(jl)   = 0.0_dp
     zrfl(jl)      = 0.0_dp
     zsfl(jl)      = 0.0_dp
111 END DO

!
!       1.2   Geopotential at half levels
!
!IBM* UNROLL_AND_FUSE(4)
  DO 132 jk = 2,klev
     DO 131 jl = 1,kproma
        zgeoh(jl,jk)   = 0.5_dp*(pgeo(jl,jk)+pgeo(jl,jk-1))
131  END DO
132 END DO
  DO 133 jl = 1,kproma
     zgeoh(jl,1)      = pgeo(jl,1)+(pgeo(jl,1)-zgeoh(jl,2))
     zgeoh(jl,klevp1) = 0.0_dp
133 END DO
!
  DO 831 jk=ktdia,klev  ! the big jk-loop
!
!       1.3   Air density
!
!IBM* NOVECTOR
     DO jl = 1,kproma
        zrho(jl,jk)      = papm1(jl,jk)/(rd*ptvm1(jl,jk))
        zqrho(jl)     = 1.3_dp/zrho(jl,jk)
        pxtec(jl,jk)  = MAX(pxtec(jl,jk),0.0_dp)
        pqtec(jl,jk)  = MAX(pqtec(jl,jk),0.0_dp)
     END DO

     zqrho_sqrt(1:kproma) = SQRT(zqrho(1:kproma))
     zpapm1_inv(1:kproma) = 1._dp/papm1(1:kproma,jk)

     CALL prepare_ua_index('cloud (1)',kproma,ptm1(1,jk),loidx(1))
     CALL lookup_ua (kproma,loidx(1),ua(1),dua(1))
     CALL lookup_uaw(kproma,loidx(1),uaw(1),duaw(1))
!
!     -------------------------------------------------------------------------
!       2.    Set to zero some local tendencies (increments)
!     -------------------------------------------------------------------------
!
!IBM* NOVECTOR
     DO 201 jl = 1,kproma

        zxisub(jl)     = 0.0_dp
        zevp(jl)       = 0.0_dp
        zsub(jl)       = 0.0_dp

        zdp(jl)        = paphm1(jl,jk+1)-paphm1(jl,jk)
        zdz(jl)        = (zgeoh(jl,jk)-zgeoh(jl,jk+1))/g

        zrcp           = 1._dp/(cpd+(cpv-cpd)*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)     = alv*zrcp
        zlsdcp(jl)     = als*zrcp
201  END DO
!
!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
#ifdef _PROFILE
      CALL trace_start ('cloud_loop_3', 13)
#endif
     IF (jk .GT. 1) THEN
!
!
!       3.1   Melting of snow and ice
!
!IBM* NOVECTOR
        DO 321 jl = 1,kproma

           zcons     = zcons2*(zdp(jl)/(zlsdcp(jl)-zlvdcp(jl)))
           ztdif     = MAX(0.0_dp,ptm1(jl,jk)-tmelt)
           zsnmlt    = MIN(zxsec*zsfl(jl),zcons*ztdif)
           zrfl(jl)  = zrfl(jl)+zsnmlt
           zsfl(jl)  = zsfl(jl)-zsnmlt
           zsmlt(jl) = zsnmlt/(zcons2*zdp(jl))
           zximelt   = MIN(zxsec*zxiflux(jl),zcons*ztdif)  ! mpuetz: zxiflux is always zero here ?
           zxiflux(jl)=zxiflux(jl)-zximelt
           zximlt(jl) = zximelt/(zcons2*zdp(jl))
           zsnmlt    = MAX(0.0_dp,pxim1(jl,jk)+pxite(jl,jk)*ztmst)
           zimlt(jl) = FSEL(-ztdif,0.0_dp,zsnmlt)

321     END DO

        IF (nclcpre.GT.0) THEN
! equals old zclcpre.gt.0
!
!       3.2   Sublimation of snow and ice (Lin et al., 1983)
!
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
           DO nl = 1,nclcpre

              jl = jjclcpre(nl)
              zesi     = ua(jl)*zpapm1_inv(jl)
              zesi     = MIN(zesi,0.5_dp)
              zqsi(jl) = zesi/(1._dp-vtmpc1*zesi)
              zsusati  = MIN(pqm1(jl,jk)/zqsi(jl)-1.0_dp,0.0_dp)
              zb1      = zlsdcp(jl)**2/(2.43e-2_dp*rv*(ptm1(jl,jk)**2))
              zb2      = 1._dp/(zrho(jl,jk)*zqsi(jl)*0.211e-4_dp)
              zcoeff(jl) = 3.e6_dp*2._dp*api*(zsusati/(zrho(jl,jk)*(zb1+zb2)))
              zclcpre_inv(jl) = 1._dp/zclcpre(jl)
           END DO
!
! definition of conditions for the sublimation of snow and ice
!CDIR NODEP
!IBM* ASSERT(NODEPS)
           DO nl = 1,nclcpre

              jl = jjclcpre(nl)
              cond1(nl) = INT(FSEL(cqtmin-zsfl(jl)   ,0._dp,1._dp))
              cond2(nl) = INT(FSEL(cqtmin-zrfl(jl)   ,0._dp,1._dp))
              cond3(nl) = INT(FSEL(cqtmin-zxiflux(jl),0._dp,1._dp))
           END DO

           i1 = 1
           i2 = 1
           i3 = 1

           DO nl = 1,nclcpre
              jl = jjclcpre(nl)

              idx1(i1) = jl
              i1 = i1 + cond1(nl)
              idx2(i2) = jl
              i2 = i2 + cond2(nl)
              idx3(i3) = jl
              i3 = i3 + cond3(nl)
           END DO

           i1 = i1 - 1
           i2 = i2 - 1
           i3 = i3 - 1

!    old if(zsfl(jl).GT.ctqmin)
!
           IF (i1.GT.0) THEN
!CDIR NODEP
!IBM* ASSERT(NODEPS)
              DO nl = 1,i1

                 jl = idx1(nl)

                 ztmp1(nl)    = zqrho_sqrt(jl)
                 ztmp2(nl)    = zsfl(jl)*zclcpre_inv(jl)/cvtfall
              END DO

              ztmp1(1:i1)   = SQRT(ztmp1(1:i1))
              ztmp2(1:i1)   = ztmp2(1:i1)**(1._dp/1.16_dp)
              ztmp2(1:i1)   = ztmp2(1:i1)/(api*crhosno*cn0s)
              ztmp2(1:i1)   = SQRT(ztmp2(1:i1))
              ztmp3(1:i1)   = ztmp2(1:i1)**1.3125_dp

!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
              DO nl = 1,i1

                 jl = idx1(nl)
                 zclcstar = zclcpre(jl)
                 zdpg     = zdp(jl)/g
                 zcfac4c  = 0.78_dp*ztmp2(nl)+232.19_dp*ztmp1(nl)*ztmp3(nl)
                 zzeps    = -zxsec*zsfl(jl)*zclcpre_inv(jl)
                 zzeps    = MAX(zzeps,zcoeff(jl)*zcfac4c*zdpg)
                 zsub(jl) = -(zzeps/zdpg)*ztmst*zclcstar
                 zsub(jl) = MIN(zsub(jl),MAX(zxsec*(zqsi(jl)-pqm1(jl,jk)),0.0_dp))
                 zsub(jl) = MAX(zsub(jl),0.0_dp)

              END DO
           END IF
!
!    end if(zsfl(jl).GT.ctqmin)

!    old if(zxiflux(jl).gt. cqtmin)
!
           IF (i3.GT.0) THEN
!CDIR NODEP
!IBM* ASSERT(NODEPS)
              DO nl = 1,i3

                 jl = idx3(nl)
                 ztmp1(nl) = zqrho_sqrt(jl)
                 ztmp2(nl) = zxiflux(jl)*zclcpre_inv(jl)/cvtfall
              END DO

              ztmp1(1:i3) = SQRT(ztmp1(1:i3))
              ztmp2(1:i3) = ztmp2(1:i3)**(1._dp/1.16_dp)
              ztmp2(1:i3) = ztmp2(1:i3)/(api*crhosno*cn0s)
              ztmp2(1:i3) = SQRT(ztmp2(1:i3))
              ztmp3(1:i3) = ztmp2(1:i3)**1.3125_dp

!CDIR NODEP
!IBM* ASSERT(NODEPS)
              DO nl = 1,i3

                 jl = idx3(nl)
                 zclcstar    = zclcpre(jl)
                 zdpg        = zdp(jl)/g
                 zcfac4c     = 0.78_dp*ztmp2(nl)+232.19_dp*ztmp1(nl)*ztmp3(nl)
                 zzeps       = -zxsec*zxiflux(jl)*zclcpre_inv(jl)
                 zzeps       = MAX(zzeps,zcoeff(jl)*zcfac4c*zdpg)
                 zsubi       = -(zzeps/zdpg)*ztmst*zclcstar
                 zsubi       = MIN(zsubi,MAX(zxsec*(zqsi(jl)-pqm1(jl,jk)),0.0_dp))
                 zsubi       = MAX(zsubi,0.0_dp)
                 zxiflux(jl) = zxiflux(jl)-zsubi*zcons2*zdp(jl)
                 zxisub(jl)  = zsubi
              END DO
           END IF
!
!    end if(zxiflux(jl).gt. cqtmin)

!
!       3.3   Evaporation of rain (Rotstayn, 1997)
!
!    old if(zrfl(jl).gt.cqtmin)
!
           IF (i2.GT.0) THEN
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
              DO nl = 1,i2

                 jl = idx2(nl)
                 zesw      = uaw(jl)*zpapm1_inv(jl)
                 zesat     = uaw(jl)/rd
                 zesw      = MIN(zesw,0.5_dp)
                 zqsw      = zesw/(1._dp-vtmpc1*zesw)
                 zsusatw   = MIN(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
                 zdv       = 2.21_dp*zpapm1_inv(jl)
                 zptm1_inv = 1._dp/ptm1(jl,jk)
                 zast      = alv*(alv*zptm1_inv/rv-1.0_dp)*zptm1_inv/0.024_dp
                 zbst      = ptm1(jl,jk)/(zdv*zesat)
                 ztmp1(nl) = zast+zbst
                 ztmp2(nl) = zrfl(jl)*zclcpre_inv(jl)
                 ztmp3(nl) = zqsw
                 ztmp4(nl) = zsusatw
              END DO

              ztmp2(1:i2) = ztmp2(1:i2)**0.61_dp

!CDIR NODEP
!IBM* ASSERT(NODEPS)
              DO nl = 1,i2

                 jl = idx2(nl)

                 zdpg     = zdp(jl)/g
                 zqsw     = ztmp3(nl)
                 zsusatw  = ztmp4(nl)
                 zclcstar = zclcpre(jl)
                 zzepr    = 870._dp*zsusatw*ztmp2(nl)*zqrho_sqrt(jl)/SQRT(1.3_dp)
                 zzepr    = zzepr/ztmp1(nl)
                 zzepr    = MAX(-zxsec*zrfl(jl)*zclcpre_inv(jl),zzepr*zdpg)
                 zevp(jl) = -(zzepr/zdpg)*ztmst*zclcstar
                 zevp(jl) = MIN(zevp(jl),MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_dp))
                 zevp(jl) = MAX(zevp(jl),0.0_dp)
              END DO
           END IF
!
!    end if(zrfl(jl).gt.cqtmin)
!

     END IF ! nclcpre.GT.0

     ELSE

        DO jl = 1,kproma
           zimlt(jl)  = 0.0_dp
           zsmlt(jl)  = 0.0_dp
           zximlt(jl) = 0.0_dp
        END DO

     END IF ! jk.GT.1
!
#ifdef _PROFILE
     CALL trace_stop ('cloud_loop_3', 13)
     CALL trace_start ('cloud_loop_4', 14)
#endif
!
!     -------------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values.
!     -------------------------------------------------------------------------
!
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!

     DO 401 jl=1,kproma
        zxip1         = pxim1(jl,jk)+pxite(jl,jk)*ztmst-zimlt(jl)
        zxip1         = MAX(zxip1,EPSILON(1._dp))
        ztmp1(jl)     = zxip1
        ztmp2(jl)     = zrho(jl,jk)*zxip1
401  END DO

     ztmp2(1:kproma) = ztmp2(1:kproma)**0.16_dp

!IBM* NOVECTOR
     DO 402 jl=1,kproma
        zxifall       = cvtfall*ztmp2(jl)
        zal1          = -zxifall*g*zrho(jl,jk)*(ztmst/zdp(jl))
        ztmp3(jl)     = zal1
402  END DO

     ztmp3(1:kproma) = EXP(ztmp3(1:kproma))

!IBM* NOVECTOR
     DO 410 jl=1,kproma
        zxip1         = ztmp1(jl)
        zxifall       = cvtfall*ztmp2(jl)
        zal2          = zxiflux(jl)/(zrho(jl,jk)*zxifall)
        zxised(jl)    = zxip1*ztmp3(jl)+zal2*(1._dp-ztmp3(jl))
        zxiflux(jl)   = zxiflux(jl)+(zxip1-zxised(jl))*zcons2*zdp(jl)
        pxite(jl,jk)  = (zxised(jl)-pxim1(jl,jk))/ztmst

        zxievap(jl)    = 0.0_dp
        zxlevap(jl)    = 0.0_dp
410  END DO

     locnt = 0
     nlocnt = 0
     DO 420 jl=1,kproma
        zclcaux(jl) = paclc(jl,jk)
        IF (zclcaux(jl) .GT. 0.0_dp) THEN    ! locc=T
           locnt = locnt + 1
           loidx(locnt) = jl
        ELSE                                 ! locc=F
           nlocnt = nlocnt + 1
           nloidx(nlocnt) = jl
        END IF
        zesw          = uaw(jl)*zpapm1_inv(jl)
        zesw          = MIN(zesw,0.5_dp)
        zqsw          = zesw/(1._dp-vtmpc1*zesw)
        zsupsatw(jl)  = MAX(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
420  END DO

!
! definition of lo2 (new zlo2(jl))
!
     DO 421 jl=1,kproma

!       IF ((ptm1(jl,jk).LT.cthomi).OR.
!         (ptm1(jl,jk).LT.tmelt.AND.zxised(jl).GT.csecfrli.AND.zsupsatw(jl).LT.zeps)) THEN
!         cond1(jl) = 1
!       ELSE
!         cond1(jl) = 0
!       END IF

        zlo2(jl)  = FSEL(ptm1(jl,jk)-tmelt, 0._dp, 1._dp)
        zlo2(jl)  = FSEL(csecfrl-zxised(jl), 0._dp, zlo2(jl))
!!$        zlo2(jl)  = FSEL(zsupsatw(jl)-zeps, 0._dp, zlo2(jl))
        zlo2(jl)  = FSEL(ptm1(jl,jk)-cthomi, zlo2(jl), 1._dp)
        cond1(jl) = INT(zlo2(jl))
        zlo2(jl)  = zlo2(jl)-0.5_dp ! zlo2 >= 0  <==> cond1 = 1
421  END DO

!
!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 430 nl=1,nlocnt
        jl = nloidx(nl)
        zxib(jl)        = 0.0_dp
        zxlb(jl)        = 0.0_dp
        zclcauxi(jl)    = 0.0_dp
430  END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
     DO 431 nl=1,locnt
        jl = loidx(nl)
        zclcauxi(jl)    = 1._dp/zclcaux(jl)
        zxib(jl)        = pxim1(jl,jk)*zclcauxi(jl)
        zxlb(jl)        = pxlm1(jl,jk)*zclcauxi(jl)
431  END DO

     i1 = 0
     i2 = 0

     IF (locnt.GT.0) THEN

        DO 440 nl=1,locnt
           jl = loidx(nl)
           IF (cond1(jl).GT.0) THEN
              i1 = i1 + 1
              idx1(i1) = jl
           ELSE
              i2 = i2 + 1
              idx2(i2) = jl
           END IF
440     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 441 nl=1,i1                       ! ice cloud
           jl = idx1(nl)

           zxite(jl)       = pxtec(jl,jk)
           zxlte(jl)       = 0.0_dp
           zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
           zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
           IF (zxidt .GT. 0.0_dp) THEN
              zxidtstar    = zxidt
              zxib(jl)     = zxib(jl)+zxidt
           ELSE
              zxidtstar    = 0.0_dp
              zxib(jl)     = zxib(jl)+MAX(zxidt*zclcauxi(jl),       &
                   -zxib(jl))
              pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                   -(pxim1(jl,jk)/ztmst+zxite(jl)))
           END IF
           IF (zxldt .GT. 0.0_dp) THEN
              zxldtstar    = zxldt
              zxlb(jl)     = zxlb(jl)+zxldt
           ELSE
              zxldtstar    = 0.0_dp
              zxlb(jl)     = zxlb(jl)+MAX(zxldt*zclcauxi(jl),        &
                   -zxlb(jl))
              pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
           END IF

           zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar
           zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar

441     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 442 nl=1,i2                       !     water cloud
           jl = idx2(nl)

           zxlte(jl)       = pxtec(jl,jk)
           zxite(jl)       = 0.0_dp
           zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                +zximlt(jl)+zimlt(jl)
           zxidt           =  pxite(jl,jk)*ztmst
           IF (zxldt .GT. 0.0_dp) THEN
              zxldtstar    = zxldt
              zxlb(jl)     = zxlb(jl)+zxldt
           ELSE
              zxldtstar    = 0.0_dp
              zxlb(jl)     = zxlb(jl)+MAX(zxldt*zclcauxi(jl),        &
                   -zxlb(jl))
              pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                   -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
           END IF
           IF (zxidt .GT. 0.0_dp) THEN
              zxidtstar    = zxidt
              zxib(jl)     = zxib(jl)+zxidt
           ELSE
              zxidtstar    = 0.0_dp
              zxib(jl)     = zxib(jl)+MAX(zxidt*zclcauxi(jl),        &
                   -zxib(jl))
              pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
           END IF

           zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar
           zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar
442     END DO

     END IF ! locnt > 0

     i3 = 0
     i4 = 0

     IF (nlocnt.GT.0) THEN

        DO 450 nl=1,nlocnt
           jl = nloidx(nl)
           IF (cond1(jl).GT.0) THEN
              i3 = i3 + 1
              idx1(i3) = jl
           ELSE
              i4 = i4 + 1
              idx2(i4) = jl
           END IF
450     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 451 nl=1,i3                       ! ice cloud
           jl = idx1(nl)

           zxlte(jl)    = 0.0_dp
           zxite(jl)    = pxtec(jl,jk)

           zxim1evp     = (pxite(jl,jk) + zxite(jl))*ztmst
           zxlm1evp     = (pxlte(jl,jk)            )*ztmst

           zxidt        = zxim1evp
           zxldt        = zxlm1evp + zximlt(jl) + zimlt(jl)

           zxidtstar    = FSEL(-zxidt,0.0_dp,zxidt)
           zxldtstar    = FSEL(-zxldt,0.0_dp,zxldt)

           zxim1evp     = pxim1(jl,jk) + FSEL(-zxidt,zxim1evp,0._dp)
           zxlm1evp     = pxlm1(jl,jk) + FSEL(-zxldt,zxlm1evp,0._dp)

           zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar+zxim1evp
           zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar+zxlm1evp

           zxiupd       = zxite(jl) - pxim1(jl,jk)/ztmst
           zxlupd       =           - pxlm1(jl,jk)/ztmst

           pxite(jl,jk) = MAX(pxite(jl,jk),FSEL(-zxidt,zxiupd,pxite(jl,jk)))
           pxlte(jl,jk) = MAX(pxlte(jl,jk),FSEL(-zxldt,zxlupd,pxlte(jl,jk)))

451     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 452 nl=1,i4                       ! water cloud
           jl = idx2(nl)

           zxlte(jl)    = pxtec(jl,jk)
           zxite(jl)    = 0.0_dp

           zxim1evp     = (pxite(jl,jk)            )*ztmst
           zxlm1evp     = (pxlte(jl,jk) + zxlte(jl))*ztmst

           zxidt        = zxim1evp
           zxldt        = zxlm1evp + zximlt(jl) + zimlt(jl)

           zxidtstar    = FSEL(-zxidt,0.0_dp,zxidt)
           zxldtstar    = FSEL(-zxldt,0.0_dp,zxldt)

           zxim1evp     = pxim1(jl,jk) + FSEL(-zxidt,zxim1evp,0._dp)
           zxlm1evp     = pxlm1(jl,jk) + FSEL(-zxldt,zxlm1evp,0._dp)

           zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar+zxim1evp
           zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar+zxlm1evp

           zxiupd       =             - pxim1(jl,jk)/ztmst
           zxlupd       = - zxlte(jl) - pxlm1(jl,jk)/ztmst

           pxite(jl,jk) = MAX(pxite(jl,jk),FSEL(-zxidt,zxiupd,pxite(jl,jk)))
           pxlte(jl,jk) = MAX(pxlte(jl,jk),FSEL(-zxldt,zxlupd,pxlte(jl,jk)))

452     END DO

     END IF ! nlocnt > 0

!
!     -------------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!     -------------------------------------------------------------------------
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
        zqsm1(jl)   = MIN(zqsm1(jl),0.5_dp)
        zcor        = 1._dp/(1._dp-vtmpc1*zqsm1(jl))
        zqsm1(jl)   = zqsm1(jl)*zcor
        zqst1       = (zua+0.001_dp*zdua)*zpapm1_inv(jl)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsm1(jl))*1000._dp
        zlcdqsdt    = zlc*zdqsdt

        zdtdt(jl)   = ptte(jl,jk)*ztmst-zlvdcp(jl)*(zevp(jl)             &
                    + zxlevap(jl))                                       &
                    - zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl))       &
                    - (zlsdcp(jl)-zlvdcp(jl))                            &
                    * (zsmlt(jl)+zximlt(jl)+zimlt(jl))

        zstar1(jl)  = zclcaux(jl)*                                      &
                    ( zlc*pqte(jl,jk)*ztmst                             &
                    + zlvdcp(jl)*(zevp(jl)+zxlevap(jl))                 &
                    + zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl))      &
                    )

        zdqsat1(jl) = zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
        zqvdt(jl)   = pqte(jl,jk)*ztmst+zevp(jl)+zsub(jl)              &
                    + zxievap(jl)+zxlevap(jl)+zxisub(jl)
        zqp1(jl)    = MAX(pqm1(jl,jk)+zqvdt(jl),0.0_dp)

        zxib(jl)    = MAX(zxib(jl),0.0_dp)
        zxlb(jl)    = MAX(zxlb(jl),0.0_dp)
!
!       Diagnostics: relative humidity
!
        zrelhum        = pqm1(jl,jk)/zqsm1(jl)
        zrelhum        = MAX(MIN(zrelhum,1._dp),0._dp)
        prelhum(jl,jk) = zrelhum
500  END DO
!
        DO 520 jl = 1,kproma

           zgenti(jl) = 0.0_dp
           zgentl(jl) = 0.0_dp

           ztp1(jl)  = ptm1(jl,jk)+zdtdt(jl)
           zdtdtstar = zdtdt(jl)+zstar1(jl)
           zdqsat    = zdtdtstar*zdqsat1(jl)
           zxilb     = zxib(jl)+zxlb(jl)

           zqcdif    = (zqvdt(jl)-zdqsat)*zclcaux(jl)
           zqcdif    = MAX(zqcdif,-zxilb*zclcaux(jl))
           zqcdif    = MIN(zqcdif,zqsec*zqp1(jl))
           ztmp1(jl) = zqcdif

520     END DO

     i1 = 0
     DO 530 jl = 1,kproma

        zqcdif = ztmp1(jl)

        IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
           zxilb    = zxib(jl)+zxlb(jl)
           zifrac   = zxib(jl)/MAX(zepsec,zxilb)
           zifrac   = MAX(MIN(zifrac,1.0_dp),0.0_dp)
           zdep(jl) = zqcdif*zifrac
           zcnd(jl) = zqcdif*(1.0_dp-zifrac)
        ELSE                                         ! cloud generation

           ! deposition (lo2 = .TRUE.) or condensation (lo2 = .FALSE.)

           zdep(jl) = FSEL(zlo2(jl),zqcdif,0.0_dp)
           zcnd(jl) = FSEL(zlo2(jl),0.0_dp,zqcdif)
        END IF
530  END DO

!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!
     DO jl = 1,kproma

        ztp1tmp(jl)  = ztp1(jl) + zlvdcp(jl)*zcnd(jl) + zlsdcp(jl)*zdep(jl)
        zqp1tmp(jl) = zqp1(jl) -            zcnd(jl) -            zdep(jl)

        zxip1       = MAX(zxised(jl)+zxite(jl)*ztmst-zxievap(jl)           &
                               +zgenti(jl)+zdep(jl),0.0_dp)

        ztmp1(jl)    = zxip1
     END DO

     CALL lookup_ubc(kproma,ztp1tmp(1),ub(1))
     CALL prepare_ua_index('cloud (2)',kproma,ztp1tmp(1),idx1(1),ztmp1(1),nphase,zlo2(1),cond1(1))
     CALL lookup_ua_eor_uaw(kproma,idx1(1),nphase,cond1(1),ua(1),dua(1))

     zpapp1i(1:kproma) = 1._dp/papp1(1:kproma,jk)
     zrhtest(1:kproma) = pqm1(1:kproma,jk)/zqsm1(1:kproma)

!IBM* NOVECTOR
     DO jl = 1,kproma

        zes         = ua(jl)*zpapp1i(jl)
        zes         = MIN(zes,0.5_dp)
        zcor        = 1._dp/(1._dp-vtmpc1*zes)
        zqsp1tmp    = zes*zcor

        zrhtest(jl) = MIN(zrhtest(jl),1._dp)*zqsp1tmp

        zqst1       = (ua(jl)+0.001_dp*dua(jl))*zpapp1i(jl)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsp1tmp)*1000._dp
        zlc         = FSEL(zlo2(jl),zlsdcp(jl),zlvdcp(jl))
        zlcdqsdt    = FSEL(zes-0.4_dp,zqsp1tmp*zcor*ub(jl),zlc*zdqsdt)
        zqcon       = 1._dp/(1._dp+zlcdqsdt)

        ztmp1(jl)   = zqcon
        ztmp2(jl)   = zqsp1tmp
     END DO

     DO 540 jl = 1,kproma
        zqcon    = ztmp1(jl)
        zqsp1tmp = ztmp2(jl)

        zoversat = zqsp1tmp*0.01_dp
        zcor     = MAX((zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon,0._dp)

        zlo      = FSEL(zqsm1(jl)-zqsp1tmp,1._dp,0._dp)
        zlo      = FSEL(zqp1tmp(jl)-zrhtest(jl),0._dp,zlo)
!
        zdefault = MAX(zqp1(jl)-zrhtest(jl),0._dp)
        zupdate  = FSEL(zlo2(jl),zdep(jl),zcnd(jl)) + zcor
        zupdate  = zupdate * (1._dp - zlo) + zlo * FSEL(-zupdate,zupdate,zdefault)
        zdep(jl) = FSEL(zlo2(jl),zupdate,zdep(jl)) ! ice cloud
        zcnd(jl) = FSEL(zlo2(jl),zcnd(jl),zupdate) ! water cloud

540  END DO

     ztmp2(1:kproma) = 1._dp/ztmp2(1:kproma)  ! mpuetz: ztmp2 holds inverse of zqsp1tmp

!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 551 nl = 1,locnt
        jl = loidx(nl)
        zxib(jl) = MAX(zxib(jl)+zdep(jl)*zclcauxi(jl),0.0_dp)
        zxlb(jl) = MAX(zxlb(jl)+zcnd(jl)*zclcauxi(jl),0.0_dp)
 551  END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 552 nl = 1,nlocnt
        jl = nloidx(nl)
        zrelhum   = zqp1tmp(jl)*ztmp2(jl)
        zdepos    = MAX(zdep(jl)+zgenti(jl),0.0_dp)
        zcond     = MAX(zcnd(jl)+zgentl(jl),0.0_dp)
        IF (zdepos>0.0_dp .OR. zcond>0.0_dp) THEN
           zclcaux(jl) =MAX(MIN(zrelhum,1.0_dp),0.01_dp)
           zclcauxi(jl)=1._dp/zclcaux(jl)
           zxib(jl)    = zdepos*zclcauxi(jl)
           zxlb(jl)    = zcond *zclcauxi(jl)
        END IF
 552  END DO

     DO 553 jl = 1,kproma
        ztp1tmp(jl) = ztp1(jl)+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
553  END DO

!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of cloud water for T < 238 K
!
     DO 610 jl = 1,kproma
        zlo = cthomi - ztp1tmp(jl)

        ! mpuetz: using compute & select is faster than branching

        zfrl(jl)  = FSEL(zlo,zxlb(jl)*zclcaux(jl),0.0_dp) ! zfrl(jl)) ! mpuetz: initially zfrl() is zero ?
        zxib(jl)  = FSEL(zlo,zxib(jl)+zxlb(jl),zxib(jl))
        zxlb(jl)  = FSEL(zlo,0.0_dp,zxlb(jl))
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
                zxlb(jl).GT.0._dp

           IF (lo) THEN
              locnt = locnt + 1
              loidx(locnt) = jl
           END IF
620     END DO

     IF (locnt > 0) THEN
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 621 nl = 1,locnt
           jl = loidx(nl)

           ztmp1(nl) = 0.66_dp*(tmelt-ztp1tmp(jl))
621     END DO

        ztmp1(1:locnt) = EXP(ztmp1(1:locnt))

!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
        DO 622 nl = 1,locnt
           jl = loidx(nl)

           zfrho     = zrho(jl,jk)/(rhoh2o*pacdnc(jl,jk))
           zfrl(jl)  = 100._dp*(ztmp1(nl)-1._dp)*zfrho
           zfrl(jl)  = zxlb(jl)*(1._dp-SWDIV_NOCHK(1._dp,(1._dp+zfrl(jl)*ztmst*zxlb(jl))))

           ztmp1(nl)  = 0.75_dp*zxlb(jl)*zfrho/api
622     END DO

        ztmp1(1:locnt) = ztmp1(1:locnt)**(1._dp/3._dp)

!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
        DO 623 nl = 1,locnt
           jl = loidx(nl)

           zradl     = ztmp1(nl)
           zval      = 4._dp*api*zradl*pacdnc(jl,jk)*2.e5_dp* (tmelt-3._dp-ztp1tmp(jl))
           zf1       = zval/zrho(jl,jk)
           zf1       = MAX(0.0_dp,zf1)
           zfrl(jl)  = zfrl(jl)+ztmst*1.4e-20_dp*zf1
           zfrl(jl)  = MAX(0.0_dp,MIN(zfrl(jl),zxlb(jl)))
           zxlb(jl)  = zxlb(jl)-zfrl(jl)
           zxib(jl)  = zxib(jl)+zfrl(jl)
           zfrl(jl)  = zfrl(jl)*zclcaux(jl)

623     END DO
     END IF
!
!       6.3 Revised Bergeron-Findeisen process
!
!!$        DO 630 jl = 1,kproma
!!$        locc        = zclcaux(jl) .GT. 0.0_dp
!!$           IF (locc .AND. zdep(jl)>0._dp .AND. zxlb(jl)>0._dp .AND.    &
!!$                          zxib(jl)>csecfrl .AND. zsupsatw(jl)<zeps) THEN
!!$              zzevp        = zxlb(jl)*zclcaux(jl)/ztmst
!!$              pxlte(jl,jk) = pxlte(jl,jk)-zzevp
!!$              pxite(jl,jk) = pxite(jl,jk)+zzevp
!!$              ptte(jl,jk)  = ptte(jl,jk)+(zlsdcp(jl)-zlvdcp(jl))*zzevp
!!$              zxib(jl)     = zxib(jl)+zxlb(jl)
!!$              zxlb(jl)     = 0.0_dp
!!$           END IF
!!$  630     END DO
!
#ifdef _PROFILE
      CALL trace_stop ('cloud_loop_4', 14)
      CALL trace_start ('cloud_loop_7', 17)
#endif
!
!     -------------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!     -------------------------------------------------------------------------
!

     zxlb(1:kproma) = MAX(zxlb(1:kproma),1.e-20_dp)
     zxib(1:kproma) = MAX(zxib(1:kproma),1.e-20_dp)

!IBM* NOVECTOR
     DO 701 jl = 1,kproma
        zauloc(jl) = cauloc*zdz(jl)/5000._dp
        zauloc(jl) = MAX(MIN(zauloc(jl),clmax),clmin)
!
! mpuetz: the following conditions are constant throughout the module
!         i.e. we should compute an index for all jl where zauloc=0.0

        jb=knvb(jl)
        lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
        lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF(lo .AND. lo1 .AND. lonacc) zauloc(jl)= 0.0_dp
701  END DO

     zxrp1(1:kproma) = 0.0_dp
     zxsp1(1:kproma) = 0.0_dp
!
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 711 nl = 1,nclcpre
        jl = jjclcpre(nl)

        ztmp1(nl) = (zrfl(jl)*zclcpre_inv(jl))/(12.45_dp*zqrho_sqrt(jl))
        ztmp2(nl) = zsfl(jl)*zclcpre_inv(jl)/cvtfall
711  END DO

     ztmp1(1:nclcpre) = ztmp1(1:nclcpre)**(8._dp/9._dp)
     ztmp2(1:nclcpre) = ztmp2(1:nclcpre)**(1._dp/1.16_dp)

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 712 nl = 1,nclcpre
        jl = jjclcpre(nl)

        zxrp1(jl) = ztmp1(nl)
        zxsp1(jl) = ztmp2(nl)
712  END DO
!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2
!
     i1 = 0
     i2 = 0

     locnt = 0
     DO 720 jl=1,kproma
        zsacl(jl) = 0.0_dp
        zrpr(jl)  = 0.0_dp
        zspr(jl)  = 0.0_dp
        IF (zclcaux(jl) .GT. 0.0_dp .AND. (zxlb(jl) > cqtmin .OR. zxib(jl) > cqtmin)) THEN
           locnt = locnt + 1
           loidx(locnt) = jl
        END IF
720  END DO

     IF (locnt.GT.0) THEN

        zexm1    = 4.7_dp-1.0_dp
        zexp     = -1._dp/zexm1

!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
        DO nl = 1,locnt
           jl = loidx(nl)
           ztmp1(nl) = (ccraut*1.2e27_dp)/zrho(jl,jk)
           ztmp2(nl) = pacdnc(jl,jk)*1.e-6_dp
           ztmp3(nl) = zrho(jl,jk)  *1.e-3_dp
           ztmp4(nl) = zxlb(jl)
        END DO

        ztmp2(1:locnt) = ztmp2(1:locnt)**(-3.3_dp)
        ztmp3(1:locnt) = ztmp3(1:locnt)**4.7_dp
        ztmp4(1:locnt) = ztmp4(1:locnt)**zexm1

        DO nl = 1,locnt
           zraut     = ztmp1(nl)*ztmp2(nl)*ztmp3(nl)
           ztmp4(nl) =  1._dp+zraut*ztmst*zexm1*ztmp4(nl)
        END DO

        ztmp4(1:locnt) = ztmp4(1:locnt)**zexp

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO nl = 1,locnt
           jl = loidx(nl)
           zraut     = zxlb(jl)*(1._dp-ztmp4(nl))
           ztmp1(nl) = -ccracl*zxrp1(jl)*ztmst
           ztmp2(nl) = -ccracl*zauloc(jl)*zrho(jl,jk)*zraut*ztmst
           ztmp3(nl) = zxib(jl)*zrho(jl,jk)*1000._dp
           ztmp4(nl) = zraut
        END DO

        ztmp1(1:locnt) = EXP(ztmp1(1:locnt))
        ztmp2(1:locnt) = EXP(ztmp2(1:locnt))
        ztmp3(1:locnt) = ztmp3(1:locnt)**0.216_dp

!CDIR NODEP
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
        DO nl = 1,locnt
           jl = loidx(nl)

           zraut = ztmp4(nl)

           zxlb(jl) = zxlb(jl)-zraut

           zrac1    = zxlb(jl)*(1._dp-ztmp1(nl))
           zxlb(jl) = zxlb(jl)-zrac1
           zrac2    = zxlb(jl)*(1._dp-ztmp2(nl))
           zxlb(jl) = zxlb(jl)-zrac2
           zclcstar = MIN(zclcaux(jl),zclcpre(jl))
           zrpr(jl) = zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1 ! zrpr is initialized to zero

        END DO
!
!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow and accretion of ice
!            by falling snow.
!            Accrection of cloud droplets by falling snow.
!            Effective radius of ice crystals after Moss (1995)
!
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO nl = 1,locnt
           jl = loidx(nl)
           zrieff    = 83.8_dp*ztmp3(nl)
           zrieff    = MIN(MAX(zrieff,ceffmin),ceffmax)

           ztmp1(nl) = 5113188._dp+2809._dp*zrieff*zrieff*zrieff
           ztmp2(nl) = zqrho(jl)
           ztmp3(nl) = 0.025_dp*(ztp1tmp(jl)-tmelt)
        END DO

        ztmp1(1:locnt) = SQRT(ztmp1(1:locnt))
        ztmp1(1:locnt) = ztmp1(1:locnt)-2261._dp
        ztmp1(1:locnt) = LOG10(ztmp1(1:locnt))
        ztmp2(1:locnt) = ztmp2(1:locnt)**0.33_dp
        ztmp3(1:locnt) = EXP(ztmp3(1:locnt))

!CDIR NODEP
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
        DO 721 nl = 1,locnt
           jl = loidx(nl)

           zc1       = 17.5_dp*zrho(jl,jk)/crhoi*ztmp2(nl)
           zdt2      = -6._dp/zc1*(ztmp1(nl)/3._dp-2._dp)
           zsaut     = ccsaut/zdt2
           zsaut     = zxib(jl)*(1._dp-1._dp/(1._dp+zsaut*ztmst*zxib(jl)))

           zxib(jl)  = zxib(jl)-zsaut
           zxsp2(nl) = zauloc(jl)*zrho(jl,jk)*zsaut
           ztmp1(nl) = zsaut

721     END DO


!IBM* NOVECTOR
        DO 722 nl = 1,locnt
           jl = loidx(nl)

           zsaut     = ztmp1(nl)
           zcolleffi = ztmp3(nl)
           zclcstar = MIN(zclcaux(jl),zclcpre(jl))

           zsaci1    = 0.0_dp
           zsaci2    = 0.0_dp
           zsacl1    = 0.0_dp
           zsacl2    = 0.0_dp

           IF (zxsp1(jl) .GT. cqtmin) THEN
              zlamsm    = (zxsp1(jl)/(api*crhosno*cn0s))**0.8125_dp
              zsaci1    = api*cn0s*3.078_dp*zlamsm*zqrho_sqrt(jl)
              zsacl1    = zxlb(jl)*(1._dp-EXP(-zsaci1*ccsacl*ztmst))
              zxlb(jl)  = zxlb(jl)-zsacl1
              zsacl1    = zclcstar*zsacl1
              zsaci1    = zsaci1*zcolleffi*ztmst
              zsaci1    = zxib(jl)*(1._dp-EXP(-zsaci1))
              zxib(jl)  = zxib(jl)-zsaci1
           END IF
           IF (zxsp2(nl) .GT. cqtmin) THEN
              zlamsm    = (zxsp2(nl)/(api*crhosno*cn0s))**0.8125_dp
              zsaci2    = api*cn0s*3.078_dp*zlamsm*zqrho_sqrt(jl)
              zsacl2    = zxlb(jl)*(1._dp-EXP(-zsaci2*ccsacl*ztmst))
              zxlb(jl)  = zxlb(jl)-zsacl2
              zsacl2    = zclcaux(jl)*zsacl2
              zsaci2    = zsaci2*zcolleffi*ztmst
              zsaci2    = zxib(jl)*(1._dp-EXP(-zsaci2))
              zxib(jl)  = zxib(jl)-zsaci2
           END IF
           zsacl(jl)    = zsacl1+zsacl2
           zspr(jl)     = zclcaux(jl)*(zsaut+zsaci2) + zclcstar*zsaci1 ! zspr is initialized to zero

722     END DO

     END IF ! locnt > 0

!
!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer
!
     zcnt = 0._dp
     nclcpre = 0

     IF (jk.EQ.klev) THEN

!IBM* NOVECTOR
        DO jl = 1,kproma

           zzdrr       = zcons2*zdp(jl)*zrpr(jl)
           zzdrs       = zcons2*zdp(jl)*(zspr(jl)+zsacl(jl))

           zzdrs       = zzdrs+zxiflux(jl)
           zcons       = (zcons2*zdp(jl))/(zlsdcp(jl)-zlvdcp(jl))
           zsnmlt      = MIN(zxsec*zzdrs,zcons*MAX(0._dp,(ztp1tmp(jl)-tmelt)))
           zzdrr       = zzdrr+zsnmlt
           zzdrs       = zzdrs-zsnmlt
           zsmlt(jl)   = zsmlt(jl)+zsnmlt/(zcons2*zdp(jl))

           zpretot     = zrfl(jl)+zsfl(jl)
           zpredel     = zzdrr+zzdrs

           zclcpre(jl) = FSEL(zpredel-zpretot,zclcaux(jl),zclcpre(jl))

           zpresum     = zpretot+zpredel

           ! This ifdef branch does contain clean code for checking
           ! the correctness of the complete code.
           IF (zpresum < TINY(zpresum)) THEN
             zclcpre1     = 0.0_dp
           ELSE
             zclcpre1     = (zclcaux(jl)*zpredel + zclcpre(jl)*zpretot)/zpresum
           ENDIF
           zclcpre1       = MAX(zclcpre(jl),zclcpre1)
           zclcpre1       = MAX(0._dp,zclcpre1)
           zclcpre1       = MIN(1._dp,zclcpre1)

           zclcpre(jl)    = FSEL(cqtmin-zpresum,0._dp,zclcpre1)

           zcnt           = zcnt + FSEL(-zclcpre(jl),0.0_dp,1.0_dp)
           cond1(jl)      = INT(FSEL(-zclcpre(jl),0.0_dp,1.0_dp))

           zrfl(jl)    = zrfl(jl)+zzdrr-zcons2*zdp(jl)*zevp(jl)
           zsfl(jl)    = zsfl(jl)+zzdrs-zcons2*zdp(jl)*zsub(jl)
        END DO

     ELSE

!IBM* NOVECTOR
        DO jl = 1,kproma

           zzdrr          = zcons2*zdp(jl)*zrpr(jl)
           zzdrs          = zcons2*zdp(jl)*(zspr(jl)+zsacl(jl))
           zpretot        = zrfl(jl)+zsfl(jl)
           zpredel        = zzdrr+zzdrs

           zclcpre(jl)    = FSEL(zpredel-zpretot,zclcaux(jl),zclcpre(jl))

           zpresum        = zpretot+zpredel

           ! This ifdef branch does contain clean code for checking
           ! the correctness of the complete code.
           IF (zpresum < TINY(zpresum)) THEN
             zclcpre1     = 0.0_dp
           ELSE
             zclcpre1     = (zclcaux(jl)*zpredel + zclcpre(jl)*zpretot)/zpresum
           ENDIF
           zclcpre1       = MAX(zclcpre(jl),zclcpre1)
           zclcpre1       = MAX(0._dp,zclcpre1)
           zclcpre1       = MIN(1._dp,zclcpre1)

           zclcpre(jl)    = FSEL(cqtmin-zpresum,0._dp,zclcpre1)

           zcnt           = zcnt + FSEL(-zclcpre(jl),0.0_dp,1.0_dp)
           cond1(jl)      = INT(FSEL(-zclcpre(jl),0.0_dp,1.0_dp))

           zrfl(jl)       = zrfl(jl)+zzdrr-zcons2*zdp(jl)*zevp(jl)
           zsfl(jl)       = zsfl(jl)+zzdrs-zcons2*zdp(jl)*zsub(jl)
        END DO

     END IF

     IF (zcnt > 0._dp) THEN
        nclcpre = 1
        DO jl = 1,kproma
           jjclcpre(nclcpre) = jl
           nclcpre = nclcpre + cond1(jl)
        END DO
        nclcpre = nclcpre - 1
     END IF

!
#ifdef _PROFILE
      CALL trace_stop ('cloud_loop_7', 17)
      CALL trace_start ('cloud_loop_8', 18)
#endif
!     -------------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!     -------------------------------------------------------------------------
!

!IBM* NOVECTOR
     DO 820 jl = 1,kproma

!       8.3   Tendencies of thermodynamic variables
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
!
        pqte(jl,jk)  = pqte(jl,jk)                                     &
                        +(-zcnd(jl)-zgentl(jl)+zevp(jl)+zxlevap(jl)    &
                          -zdep(jl)-zgenti(jl)+zsub(jl)+zxievap(jl)    &
                                          +zxisub(jl))/ztmst
        ptte(jl,jk)  = ptte(jl,jk)+(zlvdcp(jl)                         &
                        *(zcnd(jl)+zgentl(jl)-zevp(jl)-zxlevap(jl))    &
                                  +zlsdcp(jl)                          &
                        *(zdep(jl)+zgenti(jl)-zsub(jl)-zxievap(jl)     &
                        -zxisub(jl))+(zlsdcp(jl)-zlvdcp(jl))           &
                        *(-zsmlt(jl)-zimlt(jl)-zximlt(jl)+zfrl(jl)     &
                                           +zsacl(jl)))/ztmst
        pxlte(jl,jk) = pxlte(jl,jk)+zxlte(jl)                          &
                        +(zimlt(jl)+zximlt(jl)-zfrl(jl)-zrpr(jl)       &
                        -zsacl(jl)+zcnd(jl)+zgentl(jl)-zxlevap(jl))    &
                                                       /ztmst
        pxite(jl,jk) = pxite(jl,jk)+zxite(jl)+(zfrl(jl)-zspr(jl)       &
                          +zdep(jl)+zgenti(jl)-zxievap(jl))/ztmst
820  END DO

!IBM* NOVECTOR
     DO 821 jl = 1,kproma

        zxlp1        = pxlm1(jl,jk) + pxlte(jl,jk)*ztmst
        zxip1        = pxim1(jl,jk) + pxite(jl,jk)*ztmst
!
!       8.4   Corrections: Avoid negative cloud water/ice
!
        zxlold         = zxlp1
        zxiold         = zxip1
        zxlp1_d        = ccwmin - zxlp1
        zxip1_d        = ccwmin - zxip1

        zxlp1          = FSEL(-zxlp1_d,zxlp1,0._dp)
        zxip1          = FSEL(-zxip1_d,zxip1,0._dp)
        zdxlcor        = (zxlp1 - zxlold)/ztmst
        zdxicor        = (zxip1 - zxiold)/ztmst

        zxlp1_d        = MAX(zxlp1_d,0.0_dp)
        paclc(jl,jk)   = FSEL(-(zxlp1_d*zxip1_d),paclc(jl,jk),0._dp)

        pxlte(jl,jk)   = pxlte(jl,jk) + zdxlcor
        pxite(jl,jk)   = pxite(jl,jk) + zdxicor
        pqte(jl,jk)    = pqte(jl,jk) - zdxlcor - zdxicor
        ptte(jl,jk)    = ptte(jl,jk) + zlvdcp(jl)*zdxlcor + zlsdcp(jl)*zdxicor
!
821  END DO
!
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
!     -------------------------------------------------------------------------
!       10.    Diagnostics
!     -------------------------------------------------------------------------
!
!       10.1   Precipitation at the surface
!
  DO 911 jl    = 1,kproma
     prsfl(jl) = zrfl(jl)
     pssfl(jl) = zsfl(jl)
911 END DO
!
!       10.2   Total cloud cover
!
  DO 921 jl    = 1,kproma
     zclcov(jl) = 1.0_dp-paclc(jl,1)
921 END DO
  !
  DO 923 jk      = 2,klev
!IBM* NOVECTOR
     DO 922 jl    = 1,kproma
        zclcov(jl) = zclcov(jl)                                        &
             * ((1._dp-MAX(paclc(jl,jk),paclc(jl,jk-1)))/(1._dp-MIN(paclc(jl,jk-1),zxsec)))
922  END DO
923 END DO

  DO 924 jl     = 1,kproma
     zclcov(jl)  = 1.0_dp-zclcov(jl)
     paclcov(jl) = zclcov(jl)
924 END DO
!
!       10.3   Vertical integrals of humidity, cloud water and cloud ice
!
  DO 931 jl   = 1,kproma
     zqvi(jl)  = 0.0_dp
     zxlvi(jl) = 0.0_dp
     zxivi(jl) = 0.0_dp
931 END DO
!
  DO 933 jk     = ktdia,klev
     DO 932 jl   = 1,kproma
        zdpg      = (paphm1(jl,jk+1)-paphm1(jl,jk))/g
        zqvi(jl)  = zqvi(jl)+pqm1(jl,jk)*zdpg
        zxlvi(jl) = zxlvi(jl)+pxlm1(jl,jk)*zdpg
        zxivi(jl) = zxivi(jl)+pxim1(jl,jk)*zdpg
932  END DO
933 END DO
!
  DO 934 jl   = 1,kproma
     pqvi(jl)  = zqvi(jl)
     pxlvi(jl) = zxlvi(jl)
     pxivi(jl) = zxivi(jl)
934 END DO
  !
#ifdef _PROFILE
    CALL trace_stop ('cloud_loop_9', 19)
    CALL trace_stop ('cloud', 10)
#endif
!
!      10.4 Derive the tendency increment induced by this routine
!
        ptte_prc(1:kproma,:)   =  ptte(1:kproma,:)   -  ptte_prc(1:kproma,:)
        pqte_prc(1:kproma,:)   =  pqte(1:kproma,:)   -  pqte_prc(1:kproma,:)
       pxlte_prc(1:kproma,:)   = pxlte(1:kproma,:)   - pxlte_prc(1:kproma,:)
       pxite_prc(1:kproma,:)   = pxite(1:kproma,:)   - pxite_prc(1:kproma,:)

  RETURN
END SUBROUTINE cloud

END MODULE mo_cloud
