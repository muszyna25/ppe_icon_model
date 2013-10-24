!>
!! Module contains subroutine cuasc
!! Note: the subroutine contains aerosol mode specific calculation
!! and has to be cleaned later
!!
!! @author M. Tiedtke, ECMWF
!!
!! @par Revision History
!! M. Tiedtke, ECMWF (1986,1989)
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
#endif

MODULE mo_cuasc

  USE mo_kind,               ONLY : dp
  USE mo_cubasmc,            ONLY : cubasmc
  USE mo_cuentr,             ONLY : cuentr
  USE mo_cuadjtq,            ONLY : cuadjtq_idx

  USE mo_physical_constants,  ONLY : g=>grav, tmelt, vtmpc1, rv, rd, alv, als
  USE mo_echam_conv_constants,ONLY : cmfcmin, centrmax, cbfac, cmaxbuoy
  USE mo_echam_cloud_params, ONLY : csecfrl


#ifdef _PROFILE
  USE mo_profile,      ONLY : trace_start, trace_stop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cuasc

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
SUBROUTINE cuasc(  lmfdudv, lmfmid, dlev, cmfctop,                     &
           cprcon,   cminbuoy, nmctop, pdtime,   ptime_step_len,       &
           kproma,   kbdim,    klev,     klevp1, klevm1,               &
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
           kcbot,    kctop,    kctop0,                                 &
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04, UL, 28.3.07):-------
           pmwc,     pmrateprecip,  pmratesnow,                        &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
           plul,     plui,     papp1,                                  &
           pmfull,   pmfuli                                            )
!--- End Included for CDNC/IC ------------------------------------------
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

INTEGER, INTENT (IN) :: nmctop
LOGICAL, INTENT (IN) :: lmfdudv, lmfmid
REAL(dp),INTENT (IN) :: dlev, cmfctop, cprcon, cminbuoy
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac
REAL(dp),INTENT(IN) :: pdtime, ptime_step_len
INTEGER :: jl, jk, jt, ik, icall, ikb, ikt, n, locnt
!!$INTEGER :: jrow, icuasc

REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pten(kbdim,klev),        pqen(kbdim,klev),                 &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                &
            paphp1(kbdim,klevp1),    pthvsig(kbdim),                   &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                 &
            pverv(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfu(kbdim,klev),                                          &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            plu(kbdim,klev),         plude(kbdim,klev),                &
            pqude(kbdim,klev),                                         &
            pmful(kbdim,klev),       pdmfup(kbdim,klev)
REAL(dp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            klab(kbdim,klev),        kcbot(kbdim),                     &
            kctop(kbdim),            kctop0(kbdim)
INTEGER  :: khmin(kbdim),            loidx(kbdim)
REAL(dp) :: phhatt(kbdim,klev)
REAL(dp) :: phcbase(kbdim)
REAL(dp) :: pqsenh(kbdim,klev)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zmfuu(kbdim),            zmfuv(kbdim),                     &
            zpbase(kbdim),           zqold(kbdim)
REAL(dp) :: zph(kbdim)
REAL(dp) :: zodetr(kbdim,klev)
REAL(dp) :: zoentr(kbdim,klev)
REAL(dp) :: zbuoy(kbdim)
LOGICAL  :: loflag(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
REAL(dp) :: zcons2, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude, zmfuxtk     &
          , zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu, zzdmf    &
          , zdz, zdrodz, zdprho, zalvs, zmse, znevn, zodmax, zga, zdt  &
          , zscod, zqcod, zbuoyz, zscde, zdlev
!
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04, UL, 28.3.07):-------
REAL(dp) :: pmwc(kbdim,klev),        pmrateprecip(kbdim,klev),         &
            pmratesnow(kbdim,klev)
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
REAL(dp)    :: plul(kbdim,klev),        plui(kbdim,klev),              &
               pmfull(kbdim,klev),      pmfuli(kbdim,klev),            &
               zrho(kbdim,klev),                                       &
               papp1(kbdim,klev)

REAL(dp)   ::                                                          &
               ztmst,                                                  &
               zeps,                                                   &
               zcdnmin,zicemin,                                        &
               zapmin,                                                 &
               zdtime,zcsaut,                                          &
               zcraut,zsecfrl,zlift

  !--- Assumed updraft velocity in convective clouds [m s-1]:


!--- End Included for CDNC/IC ------------------------------------------
!
!      INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN, LOG

#ifdef _PROFILE
CALL trace_start ('cuasc', 40)
#endif

! IF (ltimer) CALL timer_start(timer_cuasc)

    ztmst  = ptime_step_len
    zdtime = pdtime
!   istep = get_time_step()
!
!--- Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 19/03/06)--
    zcdnmin=40.e6_dp      !
    zicemin=10._dp
    zapmin=1.e7_dp       !
    zeps=EPSILON(1._dp)
    zcraut=1._dp
    zcsaut=1._dp
    zsecfrl=csecfrl

!--- End Included for CDNC/IC ------------------------------------------

!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04, UL, 28.3.07):----
pmwc(1:kproma,:)=0._dp
pmrateprecip(1:kproma,:)=0._dp
pmratesnow(1:kproma,:)=0._dp
!---End Included for scavenging-----------------------------------------
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
!100 CONTINUE
  zcons2=1._dp/(g*ptime_step_len)
  zqold(1:kproma) = 0.0_dp

  zdlev=dlev
!
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
!200 CONTINUE
  DO 210 jl=1,kproma
     zmfuu(jl)=0._dp
     zmfuv(jl)=0._dp
     IF(.NOT.ldcum(jl)) ktype(jl)=0
210 END DO
  DO 230 jk=1,klev
     DO 220 jl=1,kproma
        plu(jl,jk)=0._dp
        pmfu(jl,jk)=0._dp
        pmfus(jl,jk)=0._dp
        pmfuq(jl,jk)=0._dp
        pmful(jl,jk)=0._dp
        plude(jl,jk)=0._dp
        pqude(jl,jk)=0._dp
        pdmfup(jl,jk)=0._dp
!--- Included for prognostic CDNC/IC scheme ----------------------------
        plul(jl,jk)=0._dp
        plui(jl,jk)=0._dp
        pmfull(jl,jk)=0._dp
        pmfuli(jl,jk)=0._dp
!--- End Included for CDNC/IC ------------------------------------------
        IF(.NOT.ldcum(jl).OR.ktype(jl).EQ.3) klab(jl,jk)=0
        IF(.NOT.ldcum(jl).AND.paphp1(jl,jk).LT.4.e4_dp) kctop0(jl)=jk
        IF(jk.LT.kcbot(jl)) klab(jl,jk)=0
220  END DO
     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pmfuxt(jl,jk,jt)=0._dp
2202    END DO
2204 END DO
!
230 END DO
  DO jk=1,klev
     DO jl=1,kproma
        zoentr(jl,jk)=0._dp
        zodetr(jl,jk)=0._dp
        zrho(jl,jk)=papp1(jl,jk)/(rd*ptu(jl,jk))
     ENDDO
  ENDDO
!
!
!----------------------------------------------------------------------
!
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
!300 CONTINUE
  DO 310 jl=1,kproma
     kctop(jl)=klevm1
     IF(.NOT.ldcum(jl)) THEN
        kcbot(jl)=klevm1
        pmfub(jl)=0._dp
        pqu(jl,klev)=0._dp
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
           pxtu(jl,klev,jt)=0._dp
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
!350 CONTINUE
  DO jl=1,kproma
     IF(ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zbuoy(jl)=g*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb) +        &
                          g*vtmpc1*(pqu(jl,ikb)-pqenh(jl,ikb))
        IF(zbuoy(jl).GT.0._dp) THEN
           zdz=(pgeo(jl,ikb-1)-pgeo(jl,ikb))/g
           zdrodz=-LOG(pten(jl,ikb-1)/pten(jl,ikb))/zdz                &
                     -g/(rd*ptenh(jl,ikb)*(1._dp+vtmpc1*pqenh(jl,ikb)))
! nb zoentr is here a fractional value
           zoentr(jl,ikb-1)=zbuoy(jl)*0.5_dp/(1._dp+zbuoy(jl)*zdz)     &
                                                              + zdrodz
           zoentr(jl,ikb-1)=MIN(zoentr(jl,ikb-1),centrmax)
           zoentr(jl,ikb-1)=MAX(zoentr(jl,ikb-1),0._dp)
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
!400 CONTINUE
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN
        CALL cubasmc(lmfdudv,  kproma, kbdim, klev, ik, klab,                    &
                     pten,     pqen,     pqsen,    puen,     pven,     &
                     ktrac,                                            &
                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv                                             &
!------------------------------added by Junhua Zhang--------------------
                    ,plul,plui,pmfull,pmfuli)
!------------------------------------end-------------------------------)
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
!CDIR NODEP
     DO 420 n=1,locnt
        jl = loidx(n)

           IF(jk.LT.kcbot(jl)) THEN
              zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zdmfen(jl)=MAX(zdmfen(jl)-MAX(zmftest-zmfmax,0._dp),0._dp)
           END IF
           zdmfde(jl)=MIN(zdmfde(jl),0.75_dp*pmfu(jl,jk+1))
           pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
           IF (ktype(jl).EQ.1 .AND. jk.LT.kcbot(jl)) THEN
              zdprho=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
              zoentr(jl,jk)=zoentr(jl,jk)*zdprho*pmfu(jl,jk+1)
              zmftest=pmfu(jl,jk)+zoentr(jl,jk)-zodetr(jl,jk)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zoentr(jl,jk)=MAX(zoentr(jl,jk)                          &
                                      -MAX(zmftest-zmfmax,0._dp),0._dp)
           ELSE
              zoentr(jl,jk)=0._dp
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
                                              *(zmse-phhatt(jl,jk+1))/g
              IF(znevn.LE.0._dp) znevn=1._dp
              zdprho=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
              zodmax=((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
              zodmax=MAX(zodmax,0._dp)
              zodetr(jl,jk)=MIN(zodetr(jl,jk),zodmax)
           ENDIF
           zodetr(jl,jk)=MIN(zodetr(jl,jk),0.75_dp*pmfu(jl,jk))
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
                          (1._dp/ptenh(jl,jk+1) + vtmpc1*zga)

           zscod=pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1)          &
                                              +pcpcu(jl,jk+1)*zdt
           zscod=MAX(zscod,0._dp)
           zscde=zscde+zodetr(jl,jk)*zscod
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           zqcod=pqsenh(jl,jk+1)+zga*zdt
           zqcod=MAX(zqcod,0._dp)
           zqude=zqude+zodetr(jl,jk)*zqcod
           pqude(jl,jk)=zqude
              plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
              plude(jl,jk)=plude(jl,jk)+plu(jl,jk+1)*zodetr(jl,jk)

           zmfusk=pmfus(jl,jk+1)+zseen-zscde
           zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
              zmfulk=pmful(jl,jk+1)    -plude(jl,jk)
              plu(jl,jk)=zmfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           pqu(jl,jk)=zmfuqk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           ptu(jl,jk)=(zmfusk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))-        &
                                pgeoh(jl,jk))/pcpcu(jl,jk)
           ptu(jl,jk)=MAX(100._dp,ptu(jl,jk))
           ptu(jl,jk)=MIN(400._dp,ptu(jl,jk))
           zqold(jl)=pqu(jl,jk)
420  END DO
!
!
     DO 4204 jt=1,ktrac
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 4202 n=1,locnt
           jl = loidx(n)
           zxteen=pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
           zxtude=pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
           zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
           pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
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


     CALL cuadjtq_idx(kproma, kbdim, klev, ik,                             &
          zph,      ptu,      pqu,      loidx, locnt,  icall)

!
!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 440 n=1,locnt
        jl = loidx(n)
        IF (pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
           zlift=MIN(zlift,1.0_dp)
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                    ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+zlift
           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.01_dp*pmfub(jl).AND.  &
                       jk.GE.kctop0(jl)) THEN
             kctop(jl)=jk
             ldcum(jl)=.TRUE.
             zdnoprc=MERGE(zdlev,1.5e4_dp,ldland(jl))
             zprcon=MERGE(0._dp,cprcon,                                &
                                   zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
             zlnew=plu(jl,jk)/                                         &
                           (1._dp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
             pdmfup(jl,jk)=MAX(0._dp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
!---Included for scavenging in xtwetdep (Philip Stier, 19/02/04):-------
             pmrateprecip(jl,jk)=plu(jl,jk)-zlnew
             pmwc(jl,jk)=plu(jl,jk)
!---End Included for scavenging-----------------------------------------
             plu(jl,jk)=zlnew
           ELSE
             klab(jl,jk)=0
             pmfu(jl,jk)=0._dp
           END IF
        END IF
440  END DO


!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 455 n=1,locnt
        jl = loidx(n)
              pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)

           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)                       &
                                    +pgeoh(jl,jk))*pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
455  END DO
     DO 4554 jt=1,ktrac
!CDIR NODEP
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
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO 460 n=1,locnt
           jl = loidx(n)
              IF(ktype(jl).EQ.1.OR.ktype(jl).EQ.3) THEN
                 zz=MERGE(3._dp,2._dp,zdmfen(jl).EQ.0._dp)
              ELSE
                 zz=MERGE(1._dp,0._dp,zdmfen(jl).EQ.0._dp)
              END IF
              zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
              zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
              zdmfdu=MIN(zdmfdu,0.75_dp*pmfu(jl,jk+1))
              zmfuu(jl)=zmfuu(jl)+                                     &
                             zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
              zmfuv(jl)=zmfuv(jl)+                                     &
                             zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
              IF(pmfu(jl,jk).GT.0._dp) THEN
                 puu(jl,jk)=zmfuu(jl)*(1._dp/pmfu(jl,jk))
                 pvu(jl,jk)=zmfuv(jl)*(1._dp/pmfu(jl,jk))
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
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO n=1,locnt
        jl = loidx(n)
        IF(ktype(jl).EQ.1) THEN
              zbuoyz=g*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) +           &
                   g*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))-g*plu(jl,jk)
           zbuoyz=MAX(zbuoyz,0.0_dp)
           zdz=(pgeo(jl,jk-1)-pgeo(jl,jk))/g
           zdrodz=-LOG(pten(jl,jk-1)/pten(jl,jk))/zdz                  &
                       -g/(rd*ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk)))
           zbuoy(jl)=zbuoy(jl)+zbuoyz*zdz
           zoentr(jl,jk-1)=zbuoyz*0.5_dp/(1._dp+zbuoy(jl)) + zdrodz
           zoentr(jl,jk-1)=MIN(zoentr(jl,jk-1),centrmax)
           zoentr(jl,jk-1)=MAX(zoentr(jl,jk-1),0._dp)
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
!500 CONTINUE
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
!DIR$ IVDEP
  DO 530 jl=1,kproma
     IF(ldcum(jl)) THEN
        jk=kctop(jl)-1
        zzdmf=cmfctop
        zdmfde(jl)=(1._dp-zzdmf)*pmfu(jl,jk+1)
           plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk)=zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
        pdmfup(jl,jk)=0._dp
        pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
!------------------------Changed by Junhua Zhang----------------
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

! IF (ltimer) CALL timer_stop(timer_cuasc)
!
  END SUBROUTINE cuasc

END MODULE mo_cuasc
