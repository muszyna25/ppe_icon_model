MODULE mo_cuascent
#ifdef __xlC__
@PROCESS STRICT
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#endif

USE mo_kind,         ONLY : wp
USE mo_physical_constants,ONLY : grav, tmelt, vtmpc1, rv, rd, alv, als
USE mo_echam_conv_constants, ONLY : lmfdudv, lmfmid, cmfcmin, cmfcmax,   &
                            cprcon, entrmid,   &
                            cmfctop, centrmax, cbfac, cminbuoy, cmaxbuoy
USE mo_cuadjust,     ONLY : cuadjtq

#ifdef _PROFILE
USE mo_profile,      ONLY : trace_start, trace_stop
#endif

IMPLICIT NONE

PRIVATE

PUBLIC :: cuasc,cuasct

CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

!this routine contains aerosol mode specific calculation and has to be cleaned later
SUBROUTINE cuasc(    kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxten,    pxtu,     pmfuxt,                       &
           pten,     pqen,     pqsen,                                  &
           pgeo,     pgeoh,    paphp1,   pthvsig,                      &
           pqte,     pverv,    klwmin,                                 &
           ldcum,    ldland,   ktype,    klab,                         &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfu,     pmfub,    pentr,                                  &
           pmfus,    pmfuq,                                            &
           pmful,    plude,    pqude,    pdmfup,                       &
           khmin,    phhatt,   phcbase,  pqsenh,                       &
           pcpen,    pcpcu,                                            &
           kcbot,    kctop,    kctop0,   nmctop,                       &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):-------
           pmwc,     pmrateprecip,       time_step_len                 )
!---End Included for scavenging-----------------------------------------
!
!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac
INTEGER :: jl, jk, jt, ik, icall, ikb, ikt, n, locnt
REAL(wp),INTENT (IN) :: time_step_len
INTEGER, INTENT(IN)  :: nmctop

REAL(wp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pten(kbdim,klev),        pqen(kbdim,klev),                 &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                &
            paphp1(kbdim,klevp1),    pthvsig(kbdim),                   &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                 &
            pverv(kbdim,klev)
!
REAL(wp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfu(kbdim,klev),                                          &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            plu(kbdim,klev),         plude(kbdim,klev),                &
            pqude(kbdim,klev),                                         &
            pmful(kbdim,klev),       pdmfup(kbdim,klev)
REAL(wp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            klab(kbdim,klev),        kcbot(kbdim),                     &
            kctop(kbdim),            kctop0(kbdim)
INTEGER  :: khmin(kbdim),            loidx(kbdim)
REAL(wp) :: phhatt(kbdim,klev)
REAL(wp) :: phcbase(kbdim)
REAL(wp) :: pqsenh(kbdim,klev)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(wp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zmfuu(kbdim),            zmfuv(kbdim),                     &
            zpbase(kbdim),           zqold(kbdim)
REAL(wp) :: zph(kbdim)
REAL(wp) :: zodetr(kbdim,klev)
REAL(wp) :: zoentr(kbdim,klev)
REAL(wp) :: zbuoy(kbdim)
LOGICAL  :: loflag(kbdim)
REAL(wp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
REAL(wp) :: zcons2, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude, zmfuxtk     &
          , zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu, zzdmf    &
          , zdz, zdrodz, zdprho, zalvs, zmse, znevn, zodmax, zga, zdt  &
          , zscod, zqcod, zbuoyz, zscde, zdlev, zlift
!
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):-------
REAL(wp) :: pmwc(kbdim,klev),        pmrateprecip(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!
!      INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN, LOG

#ifdef _PROFILE
CALL trace_start ('cuasc', 40)
#endif

!
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):----
pmwc(1:kproma,:)=0._wp
pmrateprecip(1:kproma,:)=0._wp
!---End Included for scavenging-----------------------------------------
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
  zcons2=1._wp/(grav*time_step_len)
  zqold(1:kproma) = 0.0_wp
  zdlev=3.0E4_wp
!
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
  DO 210 jl=1,kproma
     zmfuu(jl)=0._wp
     zmfuv(jl)=0._wp
     IF(.NOT.ldcum(jl)) ktype(jl)=0
210 END DO
  DO 230 jk=1,klev
     DO 220 jl=1,kproma
        plu(jl,jk)=0._wp
        pmfu(jl,jk)=0._wp
        pmfus(jl,jk)=0._wp
        pmfuq(jl,jk)=0._wp
        pmful(jl,jk)=0._wp
        plude(jl,jk)=0._wp
        pqude(jl,jk)=0._wp
        pdmfup(jl,jk)=0._wp
        IF(.NOT.ldcum(jl).OR.ktype(jl).EQ.3) klab(jl,jk)=0
        IF(.NOT.ldcum(jl).AND.paphp1(jl,jk).LT.4.e4_wp) kctop0(jl)=jk
        IF(jk.LT.kcbot(jl)) klab(jl,jk)=0
220  END DO
     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pmfuxt(jl,jk,jt)=0._wp
2202    END DO
2204 END DO
!
230 END DO
  DO jk=1,klev
     DO jl=1,kproma
        zoentr(jl,jk)=0._wp
        zodetr(jl,jk)=0._wp
     ENDDO
  ENDDO
!
!
!----------------------------------------------------------------------
!
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
  DO 310 jl=1,kproma
     kctop(jl)=klevm1
     IF(.NOT.ldcum(jl)) THEN
        kcbot(jl)=klevm1
        pmfub(jl)=0._wp
        pqu(jl,klev)=0._wp
     END IF
     pmfu(jl,klev)=pmfub(jl)
     pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)             &
                                       +pgeoh(jl,klev))
     pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
     IF(lmfdudv) THEN
        zmfuu(jl)=pmfub(jl)*puu(jl,klev)
        zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
     END IF
310 END DO
!
  DO 3112 jt=1,ktrac
     DO 3110 jl=1,kproma
        IF(.NOT.ldcum(jl)) THEN
           pxtu(jl,klev,jt)=0._wp
        ENDIF
        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
3110 END DO
3112 END DO
!
  DO 320 jl=1,kproma
     ldcum(jl)=.FALSE.
320 END DO
!
!
!
!----------------------------------------------------------------------
!
!     3.5          FIND ORGANIZED ENTRAINMENT AT CLOUD BASE
!                  ----------------------------------------
!
  DO jl=1,kproma
     IF(ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zbuoy(jl)=grav*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb) +        &
                          grav*vtmpc1*(pqu(jl,ikb)-pqenh(jl,ikb))
        IF(zbuoy(jl).GT.0._wp) THEN
           zdz=(pgeo(jl,ikb-1)-pgeo(jl,ikb))/grav
           zdrodz=-LOG(pten(jl,ikb-1)/pten(jl,ikb))/zdz                &
                     -grav/(rd*ptenh(jl,ikb)*(1._wp+vtmpc1*pqenh(jl,ikb)))
! nb zoentr is here a fractional value
           zoentr(jl,ikb-1)=zbuoy(jl)*0.5_wp/(1._wp+zbuoy(jl)*zdz)     &
                                                              + zdrodz
           zoentr(jl,ikb-1)=MIN(zoentr(jl,ikb-1),centrmax)
           zoentr(jl,ikb-1)=MAX(zoentr(jl,ikb-1),0._wp)
        ENDIF
     ENDIF
  ENDDO
!
!
!----------------------------------------------------------------------
!
!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!                  -------------------------------------------------
!
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN
        CALL cubasmc(kproma, kbdim, klev, ik, klab,                    &
                     pten,     pqen,     pqsen,    puen,     pven,     &
                     ktrac,                                            &
                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv                                             )
     ENDIF
!
     locnt = 0
     DO 410 jl=1,kproma
        IF(klab(jl,jk+1).EQ.0) klab(jl,jk)=0
        IF(klab(jl,jk+1).GT.0) THEN
          locnt = locnt + 1
          loidx(locnt) = jl
        END IF
        loflag(jl) = klab(jl,jk+1).GT.0
        zph(jl)=paphp1(jl,jk)
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
              zfac=zmfmax/pmfub(jl)
              pmfu(jl,jk+1)=pmfu(jl,jk+1)*zfac
              pmfus(jl,jk+1)=pmfus(jl,jk+1)*zfac
              pmfuq(jl,jk+1)=pmfuq(jl,jk+1)*zfac
              zmfuu(jl)=zmfuu(jl)*zfac
              zmfuv(jl)=zmfuv(jl)*zfac
           END IF
        END IF
410  END DO
     DO 4102 jt=1,ktrac
        DO 4101 jl=1,kproma
           IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
              zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
              IF(pmfub(jl).GT.zmfmax) THEN
                 zfac=zmfmax/pmfub(jl)
                 pmfuxt(jl,jk+1,jt)=pmfuxt(jl,jk+1,jt)*zfac
              END IF
           END IF
4101    END DO
4102 END DO
!
! RESET PMFUB IF NECESSARY
!
     DO 4103 jl=1,kproma
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           pmfub(jl)=MIN(pmfub(jl),zmfmax)
        END IF
4103 END DO
!
!
!*                 SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
!                  RATES PLUS ORGANIZED DETRAINMENT RATES IN *CUENTR*
!                   -------------------------------------
!
     ik=jk
     CALL cuentr(    kproma, kbdim, klev, klevp1, ik,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           zpbase,   pmfu,     pentr,    zodetr,                       &
           khmin,    pgeoh,                                            &
           zdmfen,   zdmfde)
!
!
!
!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!                  THE CLOUD ENSEMBLE ENTRAINS ENVIRONMENTAL VALUES
!                  IN TURBULENT DETRAINMENT CLOUD ENSEMBLE VALUES
!                  ARE DETRAINED
!                  IN ORGANIZED DETRAINMENT THE DRY STATIC ENERGY AND
!                  MOISTURE THAT ARE NEUTRAL COMPARED TO THE
!                  ENVIRONMENTAL AIR ARE DETRAINED
!                  ---------------------------------------------------
!
     DO 420 n=1,locnt
        jl = loidx(n)

           IF(jk.LT.kcbot(jl)) THEN
              zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zdmfen(jl)=MAX(zdmfen(jl)-MAX(zmftest-zmfmax,0._wp),0._wp)
           END IF
           zdmfde(jl)=MIN(zdmfde(jl),0.75_wp*pmfu(jl,jk+1))
           pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
           IF (ktype(jl).EQ.1 .AND. jk.LT.kcbot(jl)) THEN
              zdprho=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/grav
              zoentr(jl,jk)=zoentr(jl,jk)*zdprho*pmfu(jl,jk+1)
              zmftest=pmfu(jl,jk)+zoentr(jl,jk)-zodetr(jl,jk)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zoentr(jl,jk)=MAX(zoentr(jl,jk)                          &
                                      -MAX(zmftest-zmfmax,0._wp),0._wp)
           ELSE
              zoentr(jl,jk)=0._wp
           ENDIF
           IF(ktype(jl).EQ.1.AND.jk.LT.kcbot(jl).AND.jk.LE.khmin(jl))  &
                                                                   THEN
!          limit organized detrainment to prevent too
!          deep clouds
              zalvs=MERGE(alv,als,ptu(jl,jk+1)>tmelt)
              zmse=pcpcu(jl,jk+1)*ptu(jl,jk+1)+zalvs*pqu(jl,jk+1)      &
                                                 +pgeoh(jl,jk+1)
              ikt=kctop0(jl)
              znevn=(pgeoh(jl,ikt)-pgeoh(jl,jk+1))                     &
                                              *(zmse-phhatt(jl,jk+1))/grav
              IF(znevn.LE.0._wp) znevn=1._wp
              zdprho=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/grav
              zodmax=((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
              zodmax=MAX(zodmax,0._wp)
              zodetr(jl,jk)=MIN(zodetr(jl,jk),zodmax)
           ENDIF
           zodetr(jl,jk)=MIN(zodetr(jl,jk),0.75_wp*pmfu(jl,jk))
           pmfu(jl,jk)=pmfu(jl,jk)+zoentr(jl,jk)-zodetr(jl,jk)
           zqeen=pqenh(jl,jk+1)*zdmfen(jl)
           zqeen=zqeen+pqenh(jl,jk+1)*zoentr(jl,jk)
           zseen=(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))        &
                                             *zdmfen(jl)
           zseen=zseen+(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))  &
                                               *zoentr(jl,jk)
           zscde=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
! find moist static energy that give nonbuoyant air
           zalvs=MERGE(alv,als,ptenh(jl,jk+1)>tmelt)
           zga=zalvs*pqsenh(jl,jk+1)/(rv*(ptenh(jl,jk+1)**2))
           zdt=(plu(jl,jk+1)-vtmpc1*(pqsenh(jl,jk+1)-pqenh(jl,jk+1)))/ &
                       (1._wp/ptenh(jl,jk+1) + vtmpc1*zga)
           zscod=pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1)          &
                                              +pcpcu(jl,jk+1)*zdt
           zscod=MAX(zscod,0._wp)
           zscde=zscde+zodetr(jl,jk)*zscod
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           zqcod=pqsenh(jl,jk+1)+zga*zdt
           zqcod=MAX(zqcod,0._wp)
           zqude=zqude+zodetr(jl,jk)*zqcod
           pqude(jl,jk)=zqude
           plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
           plude(jl,jk)=plude(jl,jk)+plu(jl,jk+1)*zodetr(jl,jk)
           zmfusk=pmfus(jl,jk+1)+zseen-zscde
           zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
           zmfulk=pmful(jl,jk+1)    -plude(jl,jk)
           plu(jl,jk)=zmfulk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))
           pqu(jl,jk)=zmfuqk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))
           ptu(jl,jk)=(zmfusk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))-        &
                                pgeoh(jl,jk))/pcpcu(jl,jk)
           ptu(jl,jk)=MAX(100._wp,ptu(jl,jk))
           ptu(jl,jk)=MIN(400._wp,ptu(jl,jk))
           zqold(jl)=pqu(jl,jk)
420  END DO
!
!
     DO 4204 jt=1,ktrac
!IBM* ASSERT(NODEPS)
        DO 4202 n=1,locnt
           jl = loidx(n)
           zxteen=pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
           zxtude=pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
           zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
           pxtu(jl,jk,jt)=zmfuxtk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))
4202    END DO
4204 END DO
!
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     

     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
          zph,      ptu,      pqu,      loidx, locnt,  icall)

!
!DIR$ IVDEP
!OCL NOVREC
!IBM* ASSERT(NODEPS)
     DO 440 n=1,locnt
        jl = loidx(n)
        IF (pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
           zlift=MIN(zlift,1.0_wp)
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._wp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                    ptenh(jl,jk)*(1._wp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+zlift
           IF(zbuo.GT.0._wp.AND.pmfu(jl,jk).GE.0.01_wp*pmfub(jl).AND.  &
                       jk.GE.kctop0(jl)) THEN
             kctop(jl)=jk
             ldcum(jl)=.TRUE.
             zdnoprc=MERGE(zdlev,1.5e4_wp,ldland(jl))
             zprcon=MERGE(0._wp,cprcon,                                &
                                   zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
             zlnew=plu(jl,jk)/                                         &
                           (1._wp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
             pdmfup(jl,jk)=MAX(0._wp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
             pmrateprecip(jl,jk)=plu(jl,jk)-zlnew
             pmwc(jl,jk)=plu(jl,jk)
!---End Included for scavenging-----------------------------------------
             plu(jl,jk)=zlnew
           ELSE
             klab(jl,jk)=0
             pmfu(jl,jk)=0._wp
           END IF
        END IF
440  END DO

!IBM* ASSERT(NODEPS)
     DO 455 n=1,locnt
        jl = loidx(n)
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)                       &
                                    +pgeoh(jl,jk))*pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
455  END DO
     DO 4554 jt=1,ktrac
!IBM* ASSERT(NODEPS)
        DO 4552 n=1,locnt
           jl = loidx(n)
           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
4552    END DO
4554 END DO
!
     IF(lmfdudv) THEN
        DO jl=1,kproma
           zdmfen(jl)=zdmfen(jl)+zoentr(jl,jk)
           zdmfde(jl)=zdmfde(jl)+zodetr(jl,jk)
        ENDDO
!IBM* ASSERT(NODEPS)
        DO 460 n=1,locnt
           jl = loidx(n)
              IF(ktype(jl).EQ.1.OR.ktype(jl).EQ.3) THEN
                 zz=MERGE(3._wp,2._wp,zdmfen(jl).EQ.0._wp)
              ELSE
                 zz=MERGE(1._wp,0._wp,zdmfen(jl).EQ.0._wp)
              END IF
              zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
              zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
              zdmfdu=MIN(zdmfdu,0.75_wp*pmfu(jl,jk+1))
              zmfuu(jl)=zmfuu(jl)+                                     &
                             zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
              zmfuv(jl)=zmfuv(jl)+                                     &
                             zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
              IF(pmfu(jl,jk).GT.0._wp) THEN
                 puu(jl,jk)=zmfuu(jl)*(1._wp/pmfu(jl,jk))
                 pvu(jl,jk)=zmfuv(jl)*(1._wp/pmfu(jl,jk))
              END IF
460     END DO
     END IF
!
!
!
!                  COMPUTE ORGANIZED ENTRAINMENT
!                  FOR USE AT NEXT LEVEL
!                  ------------------------------
!
!IBM* ASSERT(NODEPS)
     DO n=1,locnt
        jl = loidx(n)
        IF(ktype(jl).EQ.1) THEN
           zbuoyz=grav*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) +           &
                grav*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))-grav*plu(jl,jk)
           zbuoyz=MAX(zbuoyz,0.0_wp)
           zdz=(pgeo(jl,jk-1)-pgeo(jl,jk))/grav
           zdrodz=-LOG(pten(jl,jk-1)/pten(jl,jk))/zdz                  &
                       -grav/(rd*ptenh(jl,jk)*(1._wp+vtmpc1*pqenh(jl,jk)))
           zbuoy(jl)=zbuoy(jl)+zbuoyz*zdz
           zoentr(jl,jk-1)=zbuoyz*0.5_wp/(1._wp+zbuoy(jl)) + zdrodz
           zoentr(jl,jk-1)=MIN(zoentr(jl,jk-1),centrmax)
           zoentr(jl,jk-1)=MAX(zoentr(jl,jk-1),0._wp)
!
        ENDIF
     ENDDO
!
!
480 END DO
!
!
!----------------------------------------------------------------------
!
!     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
!                  ----------------------------------------------------
!                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!                         FROM PREVIOUS CALCULATIONS ABOVE)
!
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
!DIR$ IVDEP
  DO 530 jl=1,kproma
     IF(ldcum(jl)) THEN
        jk=kctop(jl)-1
        zzdmf=cmfctop
        zdmfde(jl)=(1._wp-zzdmf)*pmfu(jl,jk+1)
        plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk)=zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
        pdmfup(jl,jk)=0._wp
        pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
        IF(jk.GE.2) THEN
           plude(jl,jk-1)=pmful(jl,jk)
           pqude(jl,jk-1)=pmfuq(jl,jk)
        ELSE
           plude(jl,jk)=plude(jl,jk)+pmful(jl,jk)
           pqude(jl,jk)=pqude(jl,jk)+pmfuq(jl,jk)
        END IF
     END IF
530 END DO
  DO 5312 jt=1,ktrac
     DO 5310 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
        ENDIF
5310 END DO
5312 END DO
!
  IF(lmfdudv) THEN
!DIR$      IVDEP
     DO 540 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           puu(jl,jk)=puu(jl,jk+1)
           pvu(jl,jk)=pvu(jl,jk+1)
        END IF
540  END DO
  END IF
!
#ifdef _PROFILE
  CALL trace_stop ('cuasc', 40)
#endif

END SUBROUTINE cuasc

SUBROUTINE cuasct(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxten,    pxtu,     pmfuxt,                       &
           pten,     pqen,     pqsen,                                  &
           pgeo,     pgeoh,    paphp1,   pthvsig,                      &
           pqte,     pverv,    klwmin,                                 &
           ldcum,    ldland,   ktype,    klab,                         &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfu,     pmfub,    pentr,                                  &
           pmfus,    pmfuq,                                            &
           pmful,    plude,    pqude,    pdmfup,                       &
           pcpen,    pcpcu,                                            &
           kcbot,    kctop,    kctop0,   nmctop,                       &
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
           pmwc,     pmrateprecip,       time_step_len                 )
!---End Included for scavenging-----------------------------------------
!
!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTRT*.
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac
INTEGER :: jl, jk, jt, ik, icall, is
REAL(wp),INTENT (IN) :: time_step_len
INTEGER, INTENT(IN)  :: nmctop
!
REAL(wp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pten(kbdim,klev),        pqen(kbdim,klev),                 &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                &
            paphp1(kbdim,klevp1),    pthvsig(kbdim),                   &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                 &
            pverv(kbdim,klev)
!
REAL(wp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfu(kbdim,klev),                                          &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            plu(kbdim,klev),         plude(kbdim,klev),                &
            pqude(kbdim,klev),                                         &
            pmful(kbdim,klev),       pdmfup(kbdim,klev)
REAL(wp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            klab(kbdim,klev),        kcbot(kbdim),                     &
            kctop(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(wp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zmfuu(kbdim),            zmfuv(kbdim),                     &
            zpbase(kbdim),           zqold(kbdim)
REAL(wp) :: zph(kbdim)
LOGICAL  :: loflag(kbdim)
INTEGER  :: loidx(kbdim)
REAL(wp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
!
REAL(wp) :: zcons2, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zscde, zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude       &
          , zmfuxtk, zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu  &
          , zzdmf, zdlev, zlift
!
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
REAL(wp) :: pmwc(kbdim,klev),        pmrateprecip(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!
!      INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN
!
#ifdef _PROFILE
CALL trace_start ('cuasct', 40)
#endif
!
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04, UL, 28.3.07):----
pmwc(1:kproma,:)=0._wp
pmrateprecip(1:kproma,:)=0._wp
!---End Included for scavenging-----------------------------------------
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
  zcons2=1._wp/(grav*time_step_len)
  zdlev=3.0E4_wp
!
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
  DO 210 jl=1,kproma
     zmfuu(jl)=0._wp
     zmfuv(jl)=0._wp
     IF(.NOT.ldcum(jl)) ktype(jl)=0
210 END DO
  DO 230 jk=1,klev
     DO 220 jl=1,kproma
        plu(jl,jk)=0._wp
        pmfu(jl,jk)=0._wp
        pmfus(jl,jk)=0._wp
        pmfuq(jl,jk)=0._wp
        pmful(jl,jk)=0._wp
        plude(jl,jk)=0._wp
        pqude(jl,jk)=0._wp
        pdmfup(jl,jk)=0._wp
        IF(.NOT.ldcum(jl).OR.ktype(jl).EQ.3) klab(jl,jk)=0
        IF(.NOT.ldcum(jl).AND.paphp1(jl,jk).LT.4.e4_wp) kctop0(jl)=jk
        IF(jk.LT.kcbot(jl)) klab(jl,jk)=0
220  END DO
     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pmfuxt(jl,jk,jt)=0._wp
2202    END DO
2204 END DO
!
230 END DO
!
!
!----------------------------------------------------------------------
!
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
  DO 310 jl=1,kproma
     kctop(jl)=klevm1
     IF(.NOT.ldcum(jl)) THEN
        kcbot(jl)=klevm1
        pmfub(jl)=0._wp
        pqu(jl,klev)=0._wp
     END IF
     pmfu(jl,klev)=pmfub(jl)
     pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)             &
                                       +pgeoh(jl,klev))
     pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
     IF(lmfdudv) THEN
        zmfuu(jl)=pmfub(jl)*puu(jl,klev)
        zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
     END IF
310 END DO
!
  DO 3112 jt=1,ktrac
     DO 3110 jl=1,kproma
        IF(.NOT.ldcum(jl)) THEN
           pxtu(jl,klev,jt)=0._wp
        ENDIF
        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
3110 END DO
3112 END DO
!
  DO 320 jl=1,kproma
     ldcum(jl)=.FALSE.
320 END DO
!
!
!----------------------------------------------------------------------
!
!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!                  -------------------------------------------------
!
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN
        CALL cubasmc(kproma,   kbdim,    klev,     ik,       klab,     &
                     pten,     pqen,     pqsen,    puen,     pven,     &
                     ktrac,                                            &
                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv                                            )
     ENDIF
!
     is = 0
     DO 409 jl=1,kproma
       IF (klab(jl,jk+1) == 0) klab(jl,jk) = 0
       loflag(jl) = .FALSE.
       IF (klab(jl,jk+1) > 0) THEN
         is = is+1
         loidx(is) = jl 
         loflag(jl) = .TRUE.
       ENDIF
409  ENDDO
     DO 410 jl=1,kproma
        zph(jl)=paphp1(jl,jk)
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
              zfac=zmfmax/pmfub(jl)
              pmfu(jl,jk+1)=pmfu(jl,jk+1)*zfac
              pmfus(jl,jk+1)=pmfus(jl,jk+1)*zfac
              pmfuq(jl,jk+1)=pmfuq(jl,jk+1)*zfac
              zmfuu(jl)=zmfuu(jl)*zfac
              zmfuv(jl)=zmfuv(jl)*zfac
           END IF
        END IF
410  END DO
     DO 4102 jt=1,ktrac
        DO 4101 jl=1,kproma
           IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
              zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
              IF(pmfub(jl).GT.zmfmax) THEN
                 zfac=zmfmax/pmfub(jl)
                 pmfuxt(jl,jk+1,jt)=pmfuxt(jl,jk+1,jt)*zfac
              END IF
           END IF
4101    END DO
4102 END DO
!
! RESET PMFUB IF NECESSARY
!
     DO 4103 jl=1,kproma
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           pmfub(jl)=MIN(pmfub(jl),zmfmax)
        END IF
4103 END DO
!
!
!*                 SPECIFY ENTRAINMENT RATES IN *CUENTRT*
!                  --------------------------------------
!
     ik=jk
     CALL cuentrt(kproma,  kbdim,    klev,     klevp1,   ik,           &
                 ptenh,    pqenh,    pqte,     paphp1,                 &
                 klwmin,   ldcum,    ktype,    kcbot,    kctop0,       &
                 zpbase,   pmfu,     pentr,                            &
                 zdmfen,   zdmfde)
!
!
!
!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!                  ---------------------------------------------------
!
     DO 420 jl=1,kproma
        IF(loflag(jl)) THEN
           IF(jk.LT.kcbot(jl)) THEN
              zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zdmfen(jl)=MAX(zdmfen(jl)-MAX(zmftest-zmfmax,0._wp),0._wp)
           END IF
           zdmfde(jl)=MIN(zdmfde(jl),0.75_wp*pmfu(jl,jk+1))
           pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
           zqeen=pqenh(jl,jk+1)*zdmfen(jl)
           zseen=(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))        &
                                               *zdmfen(jl)
           zscde=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))          &
                                               *zdmfde(jl)
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           pqude(jl,jk)=zqude
           plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
           zmfusk=pmfus(jl,jk+1)+zseen-zscde
           zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
           zmfulk=pmful(jl,jk+1)    -plude(jl,jk)
           plu(jl,jk)=zmfulk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))
           pqu(jl,jk)=zmfuqk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))
           ptu(jl,jk)=(zmfusk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))-        &
                                pgeoh(jl,jk))/pcpcu(jl,jk)
           ptu(jl,jk)=MAX(100._wp,ptu(jl,jk))
           ptu(jl,jk)=MIN(400._wp,ptu(jl,jk))
           zqold(jl)=pqu(jl,jk)
        END IF
420  END DO
!
!
     DO 4204 jt=1,ktrac
        DO 4202 jl=1,kproma
           IF(loflag(jl)) THEN
              zxteen=pxtenh(jl,jk+1,jt)*zdmfen(jl)
              zxtude=pxtu(jl,jk+1,jt)*zdmfde(jl)
              zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
              pxtu(jl,jk,jt)=zmfuxtk*(1._wp/MAX(cmfcmin,pmfu(jl,jk)))
           ENDIF
4202    END DO
4204 END DO
!
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL cuadjtq(kproma, kbdim, klev, ik, zph, ptu, pqu, loidx, is, icall)
!
!DIR$ IVDEP
!OCL NOVREC
     DO 440 jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
           zlift=MIN(zlift,1.0_wp)
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._wp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                    ptenh(jl,jk)*(1._wp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+zlift
           IF(zbuo.GT.0._wp.AND.pmfu(jl,jk).GE.0.1_wp*pmfub(jl)) THEN
              kctop(jl)=jk
              ldcum(jl)=.TRUE.
              zdnoprc=MERGE(zdlev,1.5e4_wp,ldland(jl))
              zprcon=MERGE(0._wp,cprcon,                               &
                               zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
              zlnew=plu(jl,jk)/                                        &
                           (1._wp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk)=MAX(0._wp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
!---Included for scavenging in wetdep_interface (Philip Stier, 19/02/04):-------
             pmrateprecip(jl,jk)=plu(jl,jk)-zlnew
             pmwc(jl,jk)=plu(jl,jk)
!---End Included for scavenging-----------------------------------------
              plu(jl,jk)=zlnew
           ELSE
              klab(jl,jk)=0
              pmfu(jl,jk)=0._wp
           END IF
        END IF
440  END DO
     DO 455 jl=1,kproma
        IF(loflag(jl)) THEN
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))         &
                                 *pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        END IF
455  END DO
     DO 4554 jt=1,ktrac
        DO 4552 jl=1,kproma
           IF(loflag(jl)) THEN
              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
           ENDIF
4552    END DO
4554 END DO
!
     IF(lmfdudv) THEN
        DO 460 jl=1,kproma
           IF(loflag(jl)) THEN
              IF(ktype(jl).EQ.1.OR.ktype(jl).EQ.3) THEN
                 zz=MERGE(3._wp,2._wp,zdmfen(jl).EQ.0._wp)
              ELSE
                 zz=MERGE(1._wp,0._wp,zdmfen(jl).EQ.0._wp)
              END IF
              zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
              zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
              zdmfdu=MIN(zdmfdu,0.75_wp*pmfu(jl,jk+1))
              zmfuu(jl)=zmfuu(jl)+                                     &
                             zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
              zmfuv(jl)=zmfuv(jl)+                                     &
                             zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
              IF(pmfu(jl,jk).GT.0._wp) THEN
                 puu(jl,jk)=zmfuu(jl)*(1._wp/pmfu(jl,jk))
                 pvu(jl,jk)=zmfuv(jl)*(1._wp/pmfu(jl,jk))
              END IF
           END IF
460     END DO
     END IF
!
480 END DO
!
!
!----------------------------------------------------------------------
!
!     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
!                  ----------------------------------------------------
!                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!                         FROM PREVIOUS CALCULATIONS ABOVE)
!
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
!DIR$ IVDEP
  DO 530 jl=1,kproma
     IF(ldcum(jl)) THEN
        jk=kctop(jl)-1
        zzdmf=cmfctop
        zdmfde(jl)=(1._wp-zzdmf)*pmfu(jl,jk+1)
        plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk)=zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
        pdmfup(jl,jk)=0._wp
        pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
        plude(jl,jk-1)=pmful(jl,jk)
        pqude(jl,jk-1)=pmfuq(jl,jk)
     END IF
530 END DO
  DO 5312 jt=1,ktrac
     DO 5310 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
        ENDIF
5310 END DO
5312 END DO
!
  IF(lmfdudv) THEN
!DIR$      IVDEP
     DO 540 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           puu(jl,jk)=puu(jl,jk+1)
           pvu(jl,jk)=pvu(jl,jk+1)
        END IF
540  END DO
  END IF
!
#ifdef _PROFILE
  CALL trace_stop ('cuasct', 40)
#endif
!
END SUBROUTINE cuasct

SUBROUTINE cubasmc(  kproma, kbdim, klev, kk, klab,                    &
           pten,     pqen,     pqsen,    puen,     pven,               &
           ktrac,                                                      &
           pxten,    pxtu,     pmfuxt,                                 &
           pverv,    pgeo,     pgeoh,    ldcum,    ktype,              &
           pmfu,     pmfub,    pentr,    kcbot,                        &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,              &
           pcpen,                                                      &
           pmfuv                                                     )
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, kk
REAL(wp) :: pten(kbdim,klev),        pqen(kbdim,klev),                 &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pqsen(kbdim,klev),       pverv(kbdim,klev),                &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev)
!
REAL(wp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            plu(kbdim,klev),         pmfu(kbdim,klev),                 &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            pmful(kbdim,klev),       pdmfup(kbdim,klev),               &
            pmfuu(kbdim),            pmfuv(kbdim)
REAL(wp) :: pcpen(kbdim,klev)
INTEGER  :: ktype(kbdim),            kcbot(kbdim),                     &
            klab(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(wp) :: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),           &
            pmfuxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jt
REAL(wp) :: zzzmb
!
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!DIR$ IVDEP
!OCL NOVREC
  DO 150 jl=1,kproma
     llo3(jl)=.FALSE.
     IF(.NOT.ldcum(jl) .AND. klab(jl,kk+1) .EQ. 0                       &
              .AND. pqen(jl,kk)    .GT. 0.90_wp*pqsen(jl,kk)            &
              .AND. pverv(jl,kk)   .LT. 0.0_wp                          &
              .AND. pgeoh(jl,kk)/grav .GT. 1500.0_wp) THEN
        llo3(jl)=.TRUE.
        ptu(jl,kk+1)=(pcpen(jl,kk)*pten(jl,kk)                         &
                        +pgeo(jl,kk)-pgeoh(jl,kk+1))/pcpen(jl,kk)
        pqu(jl,kk+1)=pqen(jl,kk)
        plu(jl,kk+1)=0._wp
        zzzmb=MAX(cmfcmin,-pverv(jl,kk)/grav)
        zzzmb=MIN(zzzmb,cmfcmax)
        pmfub(jl)=zzzmb
        pmfu(jl,kk+1)=pmfub(jl)
        pmfus(jl,kk+1)=pmfub(jl)*(pcpen(jl,kk+1)*ptu(jl,kk+1)            &
                                        +pgeoh(jl,kk+1))
        pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
        pmful(jl,kk+1)=0._wp
        pdmfup(jl,kk+1)=0._wp
        kcbot(jl)=kk
        klab(jl,kk+1)=1
        ktype(jl)=3
        pentr(jl)=entrmid
        IF(lmfdudv) THEN
           puu(jl,kk+1)=puen(jl,kk)
           pvu(jl,kk+1)=pven(jl,kk)
           pmfuu(jl)=pmfub(jl)*puu(jl,kk+1)
           pmfuv(jl)=pmfub(jl)*pvu(jl,kk+1)
        END IF
     END IF
150 END DO
!DIR$ IVDEP
!OCL NOVREC
  DO 1504 jt=1,ktrac
     DO 1502 jl=1,kproma
        IF (llo3(jl)) THEN
           pxtu(jl,kk+1,jt)=pxten(jl,kk,jt)
           pmfuxt(jl,kk+1,jt)=pmfub(jl)*pxtu(jl,kk+1,jt)
        ENDIF
1502 END DO
1504 END DO
!
!
  RETURN
END SUBROUTINE cubasmc

SUBROUTINE cuentr(   kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,    podetr,                       &
           khmin,    pgeoh,                                            &
           pdmfen,   pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
!
INTEGER, INTENT (IN) :: kbdim, klev, klevp1, kproma, kk
!
REAL(wp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                 &
            pentr(kbdim),            ppbase(kbdim)
REAL(wp) :: podetr(kbdim,klev)
REAL(wp) :: pgeoh (kbdim,klev)
INTEGER  :: khmin (kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(wp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, ikt, ikh, iklwmin, n, ncnt
REAL(wp) :: zrg, zpmid, zentr, zentest, zzmzk, ztmzk, zorgde, zarg
REAL(wp) :: zrrho(kbdim),zdprho(kbdim)
INTEGER  :: icond1(kbdim),icond2(kbdim),icond3(kbdim),idx(kbdim)
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!

#ifdef _PROFILE
  CALL trace_start ('cuentr', 41)
#endif
!
  zrg=1._wp/grav
!IBM* NOVECTOR
  DO 125 jl=1,kproma
     ppbase(jl) = paphp1(jl,kcbot(jl))
     zrrho(jl)  = (rd*ptenh(jl,kk+1)*(1._wp+vtmpc1*pqenh(jl,kk+1)))    &
                    / paphp1(jl,kk+1)
     zdprho(jl) = (paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid      = 0.5_wp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     icond1(jl) = FSEL(zpmid-paphp1(jl,kk),0._wp,1._wp)
     icond2(jl) = FSEL(0.2e5_wp - (ppbase(jl)-paphp1(jl,kk)),0._wp,1._wp)
     icond3(jl) = FSEL(1.e-5_wp-pqenh(jl,kk+1),0._wp,1._wp)
     pdmfde(jl)=0._wp
     pdmfen(jl)=0._wp
     podetr(jl,kk)=0._wp
125 END DO

  ncnt = 0
  DO jl=1,kproma
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     IF (llo1) THEN
        ncnt = ncnt+1
        idx(ncnt) = jl
     END IF
  END DO

!IBM* ASSERT(NODEPS)
  DO 126 n=1,ncnt
     jl = idx(n)
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho(jl)*zrrho(jl)
     pdmfde(jl)=zentr

     llo2=ktype(jl).EQ.2.AND.(icond2(jl).LT.1.OR.icond1(jl).GT.0)
     IF (llo2) pdmfen(jl)=zentr

     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=ktype(jl).EQ.3 .AND. kk .GE. iklwmin
     IF(llo2) pdmfen(jl)=zentr

     IF(llo2 .AND. icond3(jl).GT.0) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._wp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho(jl)))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho(jl)*zdprho(jl)
     ENDIF

     llo2=ktype(jl).EQ.1 .AND.(kk.GE.iklwmin.OR.icond1(jl).GT.0)
     IF(llo2) pdmfen(jl)=zentr
126 END DO
!
!    organized detrainment, detrainment starts at khmin
!
  IF (kproma > 0) THEN
! kproma>0: mpuetz's workaround to avoid combining of two loops
!IBM* ASSERT(NODEPS)
     DO 127 n=1,ncnt
        jl = idx(n)
        llo2=ktype(jl).EQ.1
        IF(llo2.AND.kk.LE.khmin(jl).AND.kk.GE.kctop0(jl)) THEN
           ikt=kctop0(jl)
           ikh=khmin(jl)
           IF(ikh.GT.ikt) THEN
              zzmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,kk))*zrg
              ztmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,ikt))*zrg
              zarg  =3.1415_wp*(zzmzk/ztmzk)*0.5_wp
              zorgde=TAN(zarg)*3.1415_wp*0.5_wp/ztmzk
              zdprho(jl)=(paphp1(jl,kk+1)-paphp1(jl,kk))*(zrg*zrrho(jl))
              podetr(jl,kk)=MIN(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho(jl)
           ENDIF
        ENDIF
127  END DO
  END IF
!
#ifdef _PROFILE
  CALL trace_stop ('cuentr', 41)
#endif
!
END SUBROUTINE cuentr

SUBROUTINE cuentrt(  kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,                                  &
           pdmfen,   pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
!
INTEGER, INTENT (IN) :: kbdim, klev, klevp1, kproma, kk
!
REAL(wp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                 &
            pentr(kbdim),            ppbase(kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(wp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, iklwmin
REAL(wp) :: zrg, zrrho, zdprho, zpmid, zentr, zentest
          
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
  zrg=1._wp/grav
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._wp+vtmpc1*pqenh(jl,kk+1)))           &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_wp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._wp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_wp.OR.           &
                    paphp1(jl,kk).GT.zpmid)
     pdmfen(jl)=MERGE(zentr,0._wp,llo2)
     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=llo1.AND.(ktype(jl).EQ.1 .OR. ktype(jl).EQ.3) .AND.          &
                   (kk.GE.iklwmin.OR.paphp1(jl,kk).GT.zpmid)
     IF(llo2) pdmfen(jl)=zentr
     IF(llo2 .AND. pqenh(jl,kk+1).GT.1.E-5_wp) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._wp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho*zdprho
     ENDIF
125 END DO
!
  RETURN
END SUBROUTINE cuentrt

END MODULE mo_cuascent
