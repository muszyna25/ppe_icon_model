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
! USE mo_cuadjtqi,           ONLY : cuadjtqi

#ifdef __ICON__
  USE mo_physical_constants, ONLY : g=>grav, tmelt, vtmpc1, rv, rd, alv, als
  USE mo_echam_conv_params,  ONLY : lmfdudv, lmfmid, nmctop, cmfcmin, cprcon,      &
                                    cmfctop, centrmax, cbfac, cminbuoy, cmaxbuoy,  &
                                    ncvmicro, dlev
  USE mo_echam_cloud_params, ONLY : csecfrl
!                                   cqtmin, crhosno, cn0s                 &
!                                 , cthomi, ccsacl, clmax, clmin          &
!                                 , ceffmin, ceffmax, crhoi, cauloc
! USE mo_ham_aerosol_params, ONLY : ncdnc
! USE mo_global_variables,   ONLY : ltimer
! USE mo_timer,              ONLY : timer_start, timer_stop, timer_cuasc

#else
  USE mo_control,            ONLY : nn, ltimer
  USE mo_constants,          ONLY : g, tmelt, vtmpc1, rv, rd, alv, als, cpd        &
                                  , vtmpc2, api, ak, rhoh2o
  USE mo_cumulus_flux,       ONLY : lmfdudv, lmfmid, nmctop, cmfcmin, cprcon       &
                                  , cmfctop, centrmax, cbfac, cminbuoy, cmaxbuoy
  USE mo_cloud,              ONLY : cqtmin, crhosno, cn0s                          &
                                  , cthomi, csecfrl, ccsacl, clmax, clmin          &
                                  , ceffmin, ceffmax, crhoi, cauloc
  USE mo_param_switches,     ONLY : nauto, ncvmicro
  USE mo_aero_dummy,         ONLY : ncdnc,              &
                                  & rwet,               &!UL: included for contact freezing
                                  & iaiti,iacci,icoai    !UL: included for contact freezing
  USE mo_conv,               ONLY : na_cv,cdncact_cv,cucov_tm1,ndusol_cv,nduinsolai &
                                  , nduinsolci,nbcsol_cv,nbcinsol,naerinsol         &
                                  , na_cv_diag, cdncact_cv_diag,zrain,zsnow         &
                                  , cucov_tm1_diag,zice
! USE mo_timer,              ONLY : timer_start, timer_stop, timer_cuasc
#endif

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
SUBROUTINE cuasc(    pdtime, ptime_step_len,                           &
           kproma, kbdim, klev, klevp1, klevm1,                        &
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
           pludel,   pludei,   pxtecnl,  pxtecni,                      &
           plul,     plui,     papp1,                                  &
           pmfull,   pmfuli                                            )
!!$           pwcape,   ptkem1,   jrow,     icuasc                        )
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
!!$#ifdef __ICON__
!!$#else
!!$LOGICAL :: lo2,lo
!!$#endif
REAL(dp)    :: zcdnc(kbdim,klev),                                      &
               pludel(kbdim,klev),      pludei(kbdim,klev),            &
               plul(kbdim,klev),        plui(kbdim,klev),              &
               pmfull(kbdim,klev),      pmfuli(kbdim,klev),            &
               zrho(kbdim,klev),        zicnc(kbdim,klev),             &
               zplul(kbdim,klev),       zplui(kbdim,klev),             &
!!$               pwcape(kbdim),           ptkem1(kbdim,klev),            &
               papp1(kbdim,klev),                                      &
               pxtecnl(kbdim,klev),     pxtecni(kbdim,klev)

REAL(dp)   ::                                                          &
!!$               zcucov,                                                 &
               ztmst,                                                  &
!!$               zraut,zrac2,                                            &
!!$               zsaut,zsacl2,zsaci2,zauloc,                             &
               zeps,                                                   &
!!$               zrwetki,zrwetai,zrwetci,                                &
!!$               zccbcki,     zccduai,zccduci,                           &
!!$               zdfarbcki,zdfarduai,zdfarduci,                          &
!!$               zfrzcntdu,zfrzcntbc,zfrzcnt,zfrzimm,                    &
!!$               znaimmdu,znaimmbc,                                      &
!!$               zfrl,zfrln,                                             &
!!$               zfracdusol,zfracduinsolai,zfracduinsolci,               &
!!$               zfracbcsol,zfracbcinsol,                                &
!!$               zradl,zf1,zetaair,                                      &
               zmfulkl,zmfulki,zcdnmin,zicemin,                        &
!!$               zqrho,                                                  &
!!$               zrieff,zri,zrih,                                        &
!!$               zcolleffi,zc1,zdt2,zlamsm,                              &
               zapmin,                                                 &
!!$               zexm1,zxsp2,                                            &
!!$               zexp,zxibold,                                           &
!!$               zsprn1,zself,zsprnself,                                 &
!!$               zdz2,zrprn,zxlbold,                                     &
!!$               zrautn,zrautself,zsacln,                                &
!!$               ztte,zomega,                                            &
               zdtime,zcsaut,                                          &
               zcraut,zsecfrl,zlift

  !--- Assumed updraft velocity in convective clouds [m s-1]:

!!$#ifdef __ICON__
!!$#else
!!$  REAL(dp), PARAMETER :: zwu=2.0_dp
!!$#endif

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
!    zcraut=6.3_dp
    zcraut=1._dp
!    zcsaut=350._dp
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

#ifdef __ICON__
  zdlev=dlev
#else
  IF(klev == 11) THEN
    IF(nn == 21) THEN
      zdlev=1.5E4_dp
    ELSE IF(nn == 31) THEN
      zdlev=2.0E4_dp
    ELSE
      zdlev=3.0E4_dp
    ENDIF
  ELSE
   zdlev=3.0E4_dp
  ENDIF
#endif
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
        pludel(jl,jk)=0._dp
        pludei(jl,jk)=0._dp
        pxtecnl(jl,jk)=0._dp
        pxtecni(jl,jk)=0._dp
        plul(jl,jk)=0._dp
        plui(jl,jk)=0._dp
        pmfull(jl,jk)=0._dp
        pmfuli(jl,jk)=0._dp
        zcdnc(jl,jk)=zcdnmin
        zicnc(jl,jk)=zicemin
        zplul(jl,jk)=0._dp
        zplui(jl,jk)=0._dp
!--- End Included for CDNC/IC ------------------------------------------
        IF(.NOT.ldcum(jl).OR.ktype(jl).EQ.3) klab(jl,jk)=0
        IF(.NOT.ldcum(jl).AND.paphp1(jl,jk).LT.4.e4_dp) kctop0(jl)=jk
        IF(jk.LT.kcbot(jl)) klab(jl,jk)=0
!!$#ifdef __ICON__
!!$#else
!!$        IF (ncvmicro>0 .AND. icuasc.EQ.1) THEN
!!$           zice(jl,jk,jrow)=0._dp
!!$           zrain(jl,jk,jrow)=0._dp
!!$           zsnow(jl,jk,jrow)=0._dp
!!$        ENDIF
!!$#endif
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
        CALL cubasmc(kproma, kbdim, klev, ik, klab,                    &
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
              IF (ncvmicro .EQ. 0) THEN
                 zalvs=MERGE(alv,als,ptu(jl,jk+1)>tmelt)
!!$#ifdef __ICON__
!!$#else
!!$              ELSE
!!$                 zalvs=MERGE(alv,als,ptu(jl,jk+1)>tmelt.OR.(ptu(jl,jk+1)>cthomi &
!!$                      .AND.plui(jl,jk+1)<zsecfrl))
!!$                 IF (icuasc.EQ.2) zalvs=MERGE(alv,als,ptu(jl,jk+1)>tmelt.OR. &
!!$                      (ptu(jl,jk+1)>cthomi.AND.zice(jl,jk+1,jrow)<zsecfrl))
!!$#endif
              ENDIF
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
           IF (ncvmicro .EQ. 0) THEN
              zalvs=MERGE(alv,als,ptenh(jl,jk+1)>tmelt)
!!$#ifdef __ICON__
!!$#else
!!$           ELSE
!!$              zalvs=MERGE(alv,als,ptenh(jl,jk+1)>tmelt &
!!$                   .OR.(ptenh(jl,jk+1)>cthomi.AND.plui(jl,jk+1)<zsecfrl))
!!$              IF (icuasc.EQ.2) zalvs=MERGE(alv,als,ptenh(jl,jk+1)>tmelt &
!!$                   .OR.(ptenh(jl,jk+1)>cthomi.AND.zice(jl,jk+1,jrow)<zsecfrl))
!!$#endif
           ENDIF

           zga=zalvs*pqsenh(jl,jk+1)/(rv*(ptenh(jl,jk+1)**2))

           IF (ncvmicro .EQ. 0) THEN
              zdt=(plu(jl,jk+1)-vtmpc1*(pqsenh(jl,jk+1)-pqenh(jl,jk+1)))/ &
                          (1._dp/ptenh(jl,jk+1) + vtmpc1*zga)
!!$#ifdef __ICON__
!!$#else
!!$           ELSE
!!$              zdt=(plul(jl,jk+1)+plui(jl,jk+1)-vtmpc1*(pqsenh(jl,jk+1) &
!!$                   -pqenh(jl,jk+1)))/(1._dp/ptenh(jl,jk+1) + vtmpc1*zga)
!!$#endif
           ENDIF

           zscod=pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1)          &
                                              +pcpcu(jl,jk+1)*zdt
           zscod=MAX(zscod,0._dp)
           zscde=zscde+zodetr(jl,jk)*zscod
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           zqcod=pqsenh(jl,jk+1)+zga*zdt
           zqcod=MAX(zqcod,0._dp)
           zqude=zqude+zodetr(jl,jk)*zqcod
           pqude(jl,jk)=zqude
!----------added by Junhua Zhang and Ulrike Lohmann for Micro----------
           IF (ncvmicro>0) THEN
              pludel(jl,jk)=plul(jl,jk+1)*zdmfde(jl)
              pludel(jl,jk)=pludel(jl,jk)+plul(jl,jk+1)*zodetr(jl,jk)
              pludei(jl,jk)=plui(jl,jk+1)*zdmfde(jl)
              pludei(jl,jk)=pludei(jl,jk)+plui(jl,jk+1)*zodetr(jl,jk)
           ELSE
              plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
              plude(jl,jk)=plude(jl,jk)+plu(jl,jk+1)*zodetr(jl,jk)
           ENDIF

!---------------------------------end----------------------------------

           zmfusk=pmfus(jl,jk+1)+zseen-zscde
           zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
!----------added by Junhua Zhang and Ulrike Lohmann for Micro----------
           IF (ncvmicro>0) THEN
              zmfulkl=pmfull(jl,jk+1)-pludel(jl,jk)
              zmfulki=pmfuli(jl,jk+1)-pludei(jl,jk)
              plul(jl,jk)=zmfulkl*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
              plui(jl,jk)=zmfulki*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           ELSE
              zmfulk=pmful(jl,jk+1)    -plude(jl,jk)
              plu(jl,jk)=zmfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           ENDIF
!---------------------------------end----------------------------------
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


    IF (ncvmicro .GT. 0) THEN
!!$#ifdef __ICON__
!!$#else
!!$       CALL cuadjtqi(kproma, kbdim, klev, ik,                             &
!!$            zph,      ptu,      pqu,      loflag,   icall,  plui)
!!$       IF (icuasc.EQ.2) CALL cuadjtqi(kproma, kbdim, klev, ik,            &
!!$            zph,      ptu,      pqu,      loflag,   icall,  zice(:,:,jrow))
!!$#endif
    ELSE
     CALL cuadjtq_idx(kproma, kbdim, klev, ik,                             &
          zph,      ptu,      pqu,      loidx, locnt,  icall)
    ENDIF

!
   IF (ncvmicro .EQ. 0) THEN
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
   ENDIF

!!$#ifdef __ICON__
!!$#else
!!$
!!$   IF (ncvmicro > 0) THEN
!!$     DO 441 jl=1,kproma
!!$       IF(loflag(jl)) THEN
!!$         IF (pqu(jl,jk).LT.zqold(jl)) THEN
!!$           klab(jl,jk)=2
!!$           zlift=MAX(cminbuoy,MIN(cmaxbuoy,pthvsig(jl)*cbfac))
!!$           zlift=MIN(zlift,1.0_dp)
!!$           lo2=(ptu(jl,jk).LE.cthomi).OR.                               &
!!$                ((ptu(jl,jk).LT.tmelt).AND.(plui(jl,jk).GT.zsecfrl))
!!$           IF (icuasc.EQ.2) lo2=(ptu(jl,jk).LE.cthomi).OR.              &
!!$                ((ptu(jl,jk).LT.tmelt).AND.(zice(jl,jk,jrow).GT.zsecfrl))
!!$
!!$           IF (lo2) THEN
!!$              plui(jl,jk)=plui(jl,jk)+zqold(jl)-pqu(jl,jk)
!!$           ELSE
!!$              plul(jl,jk)=plul(jl,jk)+zqold(jl)-pqu(jl,jk)
!!$           ENDIF
!!$           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plul(jl,jk)-plui(jl,jk))  &
!!$                   -ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
!!$           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+zlift
!!$           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.01_dp*pmfub(jl).AND.  &
!!$                       jk.GE.kctop0(jl)) THEN
!!$             kctop(jl)=jk
!!$             ldcum(jl)=.TRUE.
!!$!             zcucov=0.035_dp*LOG(1.0_dp+pmfu(jl,jk)*0.0981_dp*3600._dp*24._dp)
!!$             zcucov=pmfu(jl,jk)/(zwu*zrho(jl,jk))
!!$             cucov_tm1(jl,jk,jrow)=MAX(0.005_dp,MIN(zcucov,0.8_dp))
!!$             cucov_tm1_diag(jl,jk,jrow)=cucov_tm1_diag(jl,jk,jrow)+cucov_tm1(jl,jk,jrow)*zdtime
!!$
!!$
!!$!      Calculate in cloud liquid and ice water content
!!$             zplul(jl,jk)=plul(jl,jk)/cucov_tm1(jl,jk,jrow)
!!$             zplui(jl,jk)=plui(jl,jk)/cucov_tm1(jl,jk,jrow)
!!$
!!$!      calculate CDNC
!!$
!!$             IF (zplul(jl,jk)>zeps) THEN
!!$                IF (jk.EQ.kcbot(jl)) THEN
!!$                   zcdnc(jl,jk)=MAX(zcdnmin,cdncact_cv(jl,jk,jrow))
!!$                ELSE IF (jk.LT.kcbot(jl)) THEN
!!$                   zcdnc(jl,jk)=MAX(zcdnmin,cdncact_cv(jl,jk+1,jrow))
!!$                ENDIF
!!$             ENDIF
!!$             cdncact_cv_diag(jl,jk,jrow)=cdncact_cv_diag(jl,jk,jrow)+cdncact_cv(jl,jk,jrow)*zdtime
!!$             IF (ncdnc == 1) THEN
!!$               na_cv_diag(jl,jk,jrow)=na_cv_diag(jl,jk,jrow)+na_cv(jl,jk,jrow)*zdtime
!!$             ENDIF
!!$
!!$!     6.2      FREEZING OF CLOUD WATER
!!$!
!!$             IF (ptu(jl,jk).LE.cthomi) THEN
!!$                zplui(jl,jk)=zplul(jl,jk)+zplui(jl,jk)
!!$                zicnc(jl,jk)=MAX(zcdnc(jl,jk)-zcdnmin,zicemin)
!!$                zcdnc(jl,jk)=0._dp
!!$                zplul(jl,jk)=0._dp
!!$             ENDIF
!!$
!!$!    6.3    FREEZING OF CLOUD WATER BETWEEN 238 AND 273 K
!!$
!!$             lo=zplul(jl,jk).GT.zeps .AND. zcdnc(jl,jk).GE.cqtmin      &
!!$                  .AND.ptu(jl,jk).LT.tmelt.AND.ptu(jl,jk).GT.cthomi
!!$             IF (LO) THEN
!!$                zfracdusol   = MIN(ndusol_cv(jl,jk,jrow)/(cdncact_cv(jl,jk,jrow)+zeps),1._dp)
!!$                zfracduinsolai = MIN(nduinsolai(jl,jk,jrow)/(naerinsol(jl,jk,jrow)+zeps),1._dp)
!!$                zfracduinsolci = MIN(nduinsolci(jl,jk,jrow)/(naerinsol(jl,jk,jrow)+zeps),1._dp)
!!$                zfracbcsol   = MIN(nbcsol_cv(jl,jk,jrow)/(cdncact_cv(jl,jk,jrow)+zeps),1._dp)
!!$                zfracbcinsol = MIN(nbcinsol(jl,jk,jrow)/(naerinsol(jl,jk,jrow)+zeps),1._dp)
!!$                zradl     = (0.75_dp*zplul(jl,jk)*zrho(jl,jk)                      &
!!$                     /(api*rhoh2o*zcdnc(jl,jk)))**(1._dp/3._dp)
!!$                zf1       = 4._dp*api*zradl*zcdnc(jl,jk)/zrho(jl,jk)
!!$                zetaair      = 1.e-5_dp*(1.718_dp+0.0049_dp*(ptu(jl,jk)-tmelt)          &
!!$                     -1.2e-5_dp*(ptu(jl,jk)-tmelt)*(ptu(jl,jk)-tmelt))
!!$                zrwetki      = rwet(iaiti)%ptr(jl,jk,jrow)
!!$                zrwetai      = rwet(iacci)%ptr(jl,jk,jrow)
!!$                zrwetci      = rwet(icoai)%ptr(jl,jk,jrow)
!!$                zccbcki      = 1._dp+1.26_dp*6.6E-8_dp/(zrwetki+zeps) *(101325._dp/papp1(jl,jk))   &
!!$                     *(ptu(jl,jk)/tmelt)
!!$                zccduai      = 1._dp+1.26_dp*6.6E-8_dp/(zrwetai+zeps) *(101325._dp/papp1(jl,jk))   &
!!$                     *(ptu(jl,jk)/tmelt)
!!$                zccduci      = 1._dp+1.26_dp*6.6E-8_dp/(zrwetci+zeps) *(101325._dp/papp1(jl,jk))   &
!!$                     *(ptu(jl,jk)/tmelt)
!!$                IF (zrwetki==0._dp) THEN
!!$                   zdfarbcki=0._dp
!!$                ELSE
!!$                   zdfarbcki    = ak*ptu(jl,jk)*zccbcki/(6._dp*api*zetaair*(zrwetki+zeps))
!!$                ENDIF
!!$                IF (zrwetai==0._dp) THEN
!!$                   zdfarduai=0._dp
!!$                ELSE
!!$                   zdfarduai    = ak*ptu(jl,jk)*zccduai/(6._dp*api*zetaair*(zrwetai+zeps))
!!$                ENDIF
!!$                IF (zrwetci==0._dp) THEN
!!$                   zdfarduci=0._dp
!!$                ELSE
!!$                   zdfarduci    = ak*ptu(jl,jk)*zccduci/(6._dp*api*zetaair*(zrwetci+zeps))
!!$                ENDIF
!!$                zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1014_dp*(ptu(jl,jk)-tmelt)+0.3277_dp)))  ! montmorillonite
!!$                !zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1007_dp*(ptu(jl,jk)-tmelt)+0.6935_dp)))  ! kaolinite
!!$                zfrzcntbc = MIN(1._dp,MAX(0._dp,-(0.0614_dp*(ptu(jl,jk)-tmelt)+0.5730_dp)))
!!$                zfrzcnt   = zplul(jl,jk)/zcdnc(jl,jk)*zrho(jl,jk)*zf1*                     &
!!$                     (zfrzcntdu*(zdfarduai*zfracduinsolai+zdfarduci*zfracduinsolci)   &
!!$                     +zfrzcntbc*zdfarbcki*zfracbcinsol)*(zcdnc(jl,jk)+zicnc(jl,jk))
!!$                zfrzcnt   = zplul(jl,jk)*(1._dp-EXP(-zfrzcnt/zplul(jl,jk)*ztmst))
!!$                znaimmdu  = 32.3_dp*zfracdusol     ! montmorillonite
!!$                !znaimmdu  = 6.15E-2_dp*zfracdusol     ! kaolinite
!!$                znaimmbc  = 2.91E-3_dp*zfracbcsol
!!$                zomega = pverv(jl,jk) - (pwcape(jl) &
!!$                     + 1.33_dp*SQRT(ptkem1(jl,jk)))*zrho(jl,jk)*g
!!$                ztte = zomega/(zrho(jl,jk)*cpd)
!!$                zfrzimm   = -(znaimmdu+znaimmbc)*zrho(jl,jk)/rhoh2o*EXP(tmelt-ptu(jl,jk))*MIN(ztte,0._dp)
!!$                zfrzimm   = zplul(jl,jk)*(1._dp-EXP(-zfrzimm*zplul(jl,jk)/zcdnc(jl,jk)*ztmst))
!!$                zfrl  = zfrzcnt + zfrzimm
!!$                zfrl  = MAX(0.0_dp,MIN(zfrl,zplul(jl,jk)))
!!$                zplul(jl,jk)=zplul(jl,jk)-zfrl
!!$                zplui(jl,jk)=zplui(jl,jk)+zfrl
!!$                !
!!$                ! freezing of cloud droplets
!!$                !
!!$                zfrln=MAX(0.0_dp,MIN(zcdnc(jl,jk)*zfrl/(zplul(jl,jk)+zeps),zcdnc(jl,jk)-zcdnmin))
!!$                zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zfrln,cqtmin)
!!$                zicnc(jl,jk)=zicnc(jl,jk)+zfrln
!!$                !-------Bergeron-Findeisen-Process
!!$                IF (zplui(jl,jk).GT.zsecfrl) THEN
!!$                   zfrl=zfrl+zplul(jl,jk)
!!$                   zcdnc(jl,jk)=cqtmin
!!$                   zplui(jl,jk)=zplul(jl,jk)+zplui(jl,jk)
!!$                   zplul(jl,jk)=0._dp
!!$                END IF
!!$!                zfrl=(als-alv)*zfrl/(cpd+cpd*vtmpc2*MAX(pqu(jl,jk),0.0_dp))*cucov_tm1(jl,jk,jrow)
!!$!                ptu(jl,jk)=ptu(jl,jk)+zfrl
!!$             ENDIF
!!$
!!$!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!!$!             Autoconversion of cloud droplets and collection of cloud
!!$!             droplets by falling rain. Accretion of cloud droplets by
!!$!             falling snow (zsacl) is calculated under 7.2
!!$!
!!$             zraut=0._dp
!!$             zrac2=0._dp
!!$             zsaut=0._dp
!!$             zsacl2=0._dp
!!$             zsaci2=0._dp
!!$             zdz2=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
!!$             zauloc   = cauloc*zdz2/5000._dp
!!$             zauloc   = MAX(MIN(zauloc,clmax),clmin)
!!$             zdnoprc=0._dp
!!$             zdlev=293._dp+2.73_dp*zcdnc(jl,jk)*1.e-6_dp
!!$
!!$!---------Changed by Junhua Zhang for sensitivity study------------------
!!$!             zdnoprc=MERGE(zdlev,1.5e4_dp,ldland(jl))
!!$             zdnoprc=zdlev*g*zrho(jl,jk)
!!$!----------------------------end-----------------------------------------
!!$
!!$             IF(zpbase(jl)-paphp1(jl,jk).GE.zdnoprc) THEN
!!$                IF (zplul(jl,jk).GT.cqtmin.AND.zcdnc(jl,jk).GT.cqtmin) THEN
!!$                   IF (nauto==2) THEN
!!$                      !           Autoconversion rate from Khairoutdinov and Kogan, 2000
!!$                      zraut    = zcraut*1350._dp*(zcdnc(jl,jk)*1.e-6_dp)**(-1.79_dp)
!!$                      zexm1    = 2.47_dp-1.0_dp
!!$                      zexp     = -1._dp/zexm1
!!$                      zraut    = zplul(jl,jk)*(1._dp-(1._dp+zraut*ztmst*zexm1*zplul(jl,jk) &
!!$                           **zexm1)**zexp)
!!$                      zplul(jl,jk) = zplul(jl,jk)-zraut
!!$                      IF (icuasc .EQ. 1) THEN
!!$                         zrac2    =6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
!!$                         zrac2    = zplul(jl,jk)*(1._dp-EXP(-zrac2))
!!$                         zplul(jl,jk) = zplul(jl,jk)-zrac2
!!$                         zrain(jl,jk,jrow)=zrac2+zraut
!!$                      ENDIF
!!$                      IF (icuasc .EQ. 2 .AND. zrain(jl,jk,jrow).GT.cqtmin) THEN
!!$                         zrac2    =6._dp*zrain(jl,jk,jrow)*ztmst
!!$                         zrac2    = zplul(jl,jk)*(1._dp-EXP(-zrac2))
!!$                         zplul(jl,jk) = zplul(jl,jk)-zrac2
!!$                      ENDIF
!!$                      zxlbold=zplul(jl,jk)+zrac2+zraut
!!$                      zrprn=(zraut+zrac2)/(zxlbold+zeps)
!!$                      IF (zplul(jl,jk) .GT. cqtmin) THEN
!!$                         zrprn=MIN(zcdnc(jl,jk)*zrprn,zcdnc(jl,jk)-zcdnmin)
!!$                      ELSE
!!$                         zrprn=MIN(zcdnc(jl,jk)*zrprn,zcdnc(jl,jk))
!!$                      END IF
!!$                      zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrprn,cqtmin)
!!$                   ELSE IF (nauto==1) THEN
!!$                      zraut    = zcraut*1.2e27_dp/zrho(jl,jk)*(zcdnc(jl,jk)*1.e-6_dp)  &
!!$                           **(-3.3_dp)*(zrho(jl,jk)*1.e-3_dp)**4.7_dp
!!$                      zexm1    = 4.7_dp-1.0_dp
!!$                      zexp     = -1._dp/zexm1
!!$                      zraut    = zplul(jl,jk)*(1._dp-(1._dp+zraut*ztmst*zexm1*zplul(jl,jk) &
!!$                           **zexm1)**zexp)
!!$                      zrautn=zraut*7.7e9_dp*zrho(jl,jk)
!!$                      zself=1.289e10_dp*(zrho(jl,jk)*zplul(jl,jk))**2*ztmst
!!$                      zrautself=MIN(zrautn+zself,zcdnc(jl,jk))
!!$                      zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrautself,cqtmin)
!!$                      zplul(jl,jk) = zplul(jl,jk)-zraut
!!$                      IF (icuasc .EQ. 1) THEN
!!$                         zrac2    =6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
!!$                         zrac2    = zplul(jl,jk)*(1._dp-EXP(-zrac2))
!!$                         zplul(jl,jk) = zplul(jl,jk)-zrac2
!!$                         zrain(jl,jk,jrow)=zrac2+zraut
!!$                      ENDIF
!!$                      IF (icuasc .EQ. 2 .AND. zrain(jl,jk,jrow).GT.cqtmin) THEN
!!$                         zrac2    =6._dp*zrain(jl,jk,jrow)*ztmst
!!$                         zrac2    = zplul(jl,jk)*(1._dp-EXP(-zrac2))
!!$                         zplul(jl,jk) = zplul(jl,jk)-zrac2
!!$                      ENDIF
!!$                      zxlbold=zplul(jl,jk)+zrac2+zraut
!!$                      zrprn=zrac2/(zxlbold+zeps)
!!$                      IF (zplul(jl,jk) .GT. cqtmin) THEN
!!$                         zrprn=MIN(zcdnc(jl,jk)*zrprn,zcdnc(jl,jk)-zcdnmin)
!!$                      ELSE
!!$                         zrprn=MIN(zcdnc(jl,jk)*zrprn,zcdnc(jl,jk))
!!$                      END IF
!!$                      zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrprn,cqtmin)
!!$                   ENDIF ! nauto
!!$                ENDIF ! plul > 0
!!$!
!!$!       7.2  Cold clouds:
!!$!            Conversion of cloud ice to snow after Levkov et al. 1992:
!!$!            Aggregation of ice crystals to snow and accretion of ice crystals
!!$!            and cloud water by falling snow.
!!$!            Effective radius of ice crystals after Moss (1995)
!!$!
!!$                IF (zplui(jl,jk).GT.zeps ) THEN
!!$                   zqrho     = 1.3_dp/zrho(jl,jk)
!!$                   zrieff    = 83.8_dp*(zplui(jl,jk)*zrho(jl,jk)*1000._dp)**0.216_dp
!!$                   zrieff    = MIN(MAX(zrieff,ceffmin),ceffmax)
!!$                   zrih      = -2261._dp+SQRT(5113188._dp+2809._dp*zrieff*zrieff*zrieff)
!!$                   zri       = 1.e-6_dp*zrih**(1._dp/3._dp)
!!$                   zcolleffi = EXP(0.025_dp*(ptu(jl,jk)-tmelt))
!!$                   zc1       = 17.5_dp*zrho(jl,jk)/crhoi*zqrho**0.33_dp
!!$                   zdt2      = -6._dp/zc1*LOG10(zri*1.e4_dp)
!!$                   zsaut     = zcsaut/zdt2
!!$                   zsaut     = zplui(jl,jk)*(1._dp-1._dp/(1._dp+zsaut*ztmst*zplui(jl,jk)))
!!$                   zplui(jl,jk)  = zplui(jl,jk)-zsaut
!!$                   zxsp2        = zauloc*zrho(jl,jk)*zsaut
!!$                   IF (icuasc .EQ. 1 .AND. zxsp2 .GT. cqtmin) THEN
!!$                      zlamsm    = (zxsp2/(api*crhosno*cn0s))**0.8125_dp
!!$                      zsaci2    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
!!$                      IF (zplul(jl,jk) .GT. cqtmin) THEN
!!$                         zsacl2    = zplul(jl,jk)*(1._dp-EXP(-zsaci2*ccsacl*ztmst))
!!$                         zplul(jl,jk) = zplul(jl,jk) - zsacl2
!!$                      ENDIF
!!$                      IF (zplui(jl,jk) .GT. cqtmin) THEN
!!$                         zsaci2    = zsaci2*zcolleffi*ztmst
!!$                         zsaci2    = zplui(jl,jk)*(1._dp-EXP(-zsaci2))
!!$                         zplui(jl,jk)      = zplui(jl,jk)-zsaci2
!!$                      ENDIF
!!$                      zsnow(jl,jk,jrow)=(zsaci2+zsaut)*zrho(jl,jk)
!!$                   END IF
!!$                   IF (icuasc .EQ. 2 .AND. zsnow(jl,jk,jrow) .GT. cqtmin) THEN
!!$                      zlamsm    = (zsnow(jl,jk,jrow)/(api*crhosno*cn0s))**0.8125_dp
!!$                      zsaci2    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
!!$                      IF (zplul(jl,jk) .GT. cqtmin) THEN
!!$                         zsacl2    = zplul(jl,jk)*(1._dp-EXP(-zsaci2*ccsacl*ztmst))
!!$                         zplul(jl,jk) = zplul(jl,jk) - zsacl2
!!$                      ENDIF
!!$                      IF (zplui(jl,jk) .GT. cqtmin) THEN
!!$                         zsaci2    = zsaci2*zcolleffi*ztmst
!!$                         zsaci2    = zplui(jl,jk)*(1._dp-EXP(-zsaci2))
!!$                         zplui(jl,jk)      = zplui(jl,jk)-zsaci2
!!$                      ENDIF
!!$                   END IF
!!$                   zxibold=MAX(zplui(jl,jk)+zsaut+zsaci2,0._dp)
!!$                   zsprn1=zicnc(jl,jk)*(zsaci2+zsaut)/(zxibold+zeps)
!!$                   zself=zc1*0.5_dp*zicnc(jl,jk)*ztmst*zplui(jl,jk)
!!$                   zsprnself=MIN(zsprn1+zself,zicnc(jl,jk))
!!$                   zicnc(jl,jk)=MAX(zicnc(jl,jk)-zsprnself,cqtmin)
!!$                   zxlbold=zplul(jl,jk)+zsacl2
!!$                   IF (zplul(jl,jk) .GT. cqtmin) THEN
!!$                      zsacln=MAX( MIN(zcdnc(jl,jk)*(zsacl2)/(zxlbold+zeps),&
!!$                           zcdnc(jl,jk)-zcdnmin ),0._dp )
!!$                   ELSE
!!$                      zsacln=MIN(zcdnc(jl,jk)*(zsacl2)/(zxlbold+zeps),&
!!$                           zcdnc(jl,jk))
!!$                   END IF
!!$                   zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zsacln,cqtmin)
!!$                END IF  ! plui > 0.
!!$
!!$                pdmfup(jl,jk)=MAX(0._dp,(MIN((zraut+zrac2+zsaut+zsacl2+zsaci2)* &
!!$                     cucov_tm1(jl,jk,jrow),plul(jl,jk)+plui(jl,jk)))*pmfu(jl,jk))
!!$                zice(jl,jk,jrow)=zplui(jl,jk)*cucov_tm1(jl,jk,jrow)
!!$                !---Included for scavenging in xtwetdep (UL, 28/05/06):-------
!!$                pmrateprecip(jl,jk)=MAX(0._dp,(zraut+zrac2)*cucov_tm1(jl,jk,jrow))
!!$                pmratesnow(jl,jk)=MAX(0._dp,(zsaut+zsacl2+zsaci2)*cucov_tm1(jl,jk,jrow))
!!$                !---End Included for scavenging-----------------------------------------
!!$             END IF  ! zpbase-zaphp1 > zdnoprc
!!$          ELSE
!!$             klab(jl,jk)=0
!!$             pmfu(jl,jk)=0._dp
!!$          END IF  ! zbuo > 0, +ve mass flux
!!$       END IF  ! pqu > zqold
!!$    END IF  ! loflag
!!$441  END DO
!!$   ENDIF
!!$#endif

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO 455 n=1,locnt
        jl = loidx(n)
!-----------------------Changed by Junhua Zhang----------------------
           IF(ncvmicro>0) THEN
              pmfull(jl,jk)=plul(jl,jk)*pmfu(jl,jk)
              pmfuli(jl,jk)=plui(jl,jk)*pmfu(jl,jk)
           ELSE
              pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           ENDIF
!-------------------------------end----------------------------------

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
           IF (ncvmicro .EQ. 0) THEN
              zbuoyz=g*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) +           &
                   g*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))-g*plu(jl,jk)
           ELSE
              zbuoyz=g*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) +           &
                   g*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))-g*(plul(jl,jk)+plui(jl,jk))
           ENDIF
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
!--- Included for prognostic CDNC/IC scheme ----------------------------
        IF(ncvmicro>0) THEN
           pludel(jl,jk)=zdmfde(jl)*plul(jl,jk+1)
           pludei(jl,jk)=zdmfde(jl)*plui(jl,jk+1)
        ELSE
!--- End Included for CDNC/IC ------------------------------------------
           plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
        ENDIF
        pqude(jl,jk)=zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
        pdmfup(jl,jk)=0._dp
        pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
!------------------------Changed by Junhua Zhang----------------
        IF(ncvmicro>0) THEN
           pmfull(jl,jk)=plul(jl,jk)*pmfu(jl,jk)
           pmfuli(jl,jk)=plui(jl,jk)*pmfu(jl,jk)
           pxtecnl(jl,jk)=zcdnc(jl,jk)
           pxtecni(jl,jk)=zicnc(jl,jk)
           IF(jk.GE.2) THEN
              pludel(jl,jk-1)=pmfull(jl,jk)
              pludei(jl,jk-1)=pmfuli(jl,jk)
              pqude(jl,jk-1)=pmfuq(jl,jk)
           ELSE
              pludel(jl,jk)=pludel(jl,jk)+pmfull(jl,jk)
              pludei(jl,jk)=pludei(jl,jk)+pmfuli(jl,jk)
              pqude(jl,jk)=pqude(jl,jk)+pmfuq(jl,jk)
           END IF
        ELSE
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           IF(jk.GE.2) THEN
              plude(jl,jk-1)=pmful(jl,jk)
              pqude(jl,jk-1)=pmfuq(jl,jk)
           ELSE
              plude(jl,jk)=plude(jl,jk)+pmful(jl,jk)
              pqude(jl,jk)=pqude(jl,jk)+pmfuq(jl,jk)
           END IF
        ENDIF
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
