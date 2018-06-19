! $RCSfile$
! $Revision$ $Date$
!
!>
!! *cuflxn* THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!! *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!!                DOES GLOBAL DIAGNOSTICS
!! *CUDUDV* - UPDATES U AND V TENDENCIES,
!!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!!
!! @author    M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!!
!!
!! *CUCTRACER* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
!!                   IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES
!!
!! *cubidiag* - SOLVES BIDIAGONAL SYSTEM FOR IMPLICIT SOLUTION
!!              OF ADVECTION EQUATION
!!
!! @author    P.BECHTOLD        E.C.M.W.F.              11/02/2004
!!
!! For detailed description refer to subsequent subroutines
!!
!! @par Revision History
!! first implementation into GME/ICON by Kristina Froehlich, DWD (2010-05-27)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#include "consistent_fma.inc"
MODULE mo_cuflxtends

#ifdef __ICON__
  USE mo_kind   ,ONLY: jprb=>wp     , &
    &                  jpim=>i4
#endif

#ifdef __GME__
!  USE parkind1  ,ONLY : jpim     ,jprb
  USE gme_data_parameters, ONLY:  JPRB =>ireals, JPIM => iintegers
#endif

!  USE yomhook   ,ONLY : lhook,   dr_hook
  !KF
  USE mo_cuparameters, ONLY: lphylin  ,rlptrc,  lepcld              ,&
    &                        rcpecons ,rtaumel ,& !rcucov ,rhebc    ,&
    &                        rmfsoltq,  rmfsoluv                    ,&
    &                        rmfsolct, rmfcmin,rg       ,rcpd       ,&
    &                        rlvtt   , rlstt    ,rlmlt    ,rtt      ,&
    &                        lhook, dr_hook, lmfglac, lmfwetb

  USE mo_cufunctions, ONLY: foelhmcu, foeewmcu, foealfcu, &
    & foeewl,   foeewi

  USE mo_fortran_tools, ONLY: t_ptr_tracer

!K.L. for testing:
  USE mo_exception,           ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cuflxn, cudtdqn,cududv,cuctracer

CONTAINS


  !OPTIONS XOPT(HSFUN)
  SUBROUTINE cuflxn &
    & (  kidia,    kfdia,    klon,   ktdia,   klev, rmfcfl, &
    & rhebc_land, rhebc_ocean, rcucov, rhebc_land_trop,     &  
    & rhebc_ocean_trop, rcucov_trop, lmfdsnow, trop_mask, ptsphy, &
    & pten,     pqen,     pqsen,    ptenh,    pqenh,&
    & paph,     pap,      pgeoh,    ldland,   ldlake, ldcum,&
    & kcbot,    kctop,    kdtop,    ktopm2,&
    & ktype,    lddraf,&
    & pmfu,     pmfd,     pmfus,    pmfds,&
    & pmfuq,    pmfdq,    pmful,    plude,    plrain, psnde, &
    & pdmfup,   pdmfdp,   pdpmel,   plglac,&
    & pmflxr,   pmflxs,   prain,    pmfude_rate, pmfdde_rate )
    !
    !!Description:
    !!         M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

    !!         PURPOSE
    !!         -------

    !!         THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
    !!         FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER

    !!         INTERFACE
    !!         ---------
    !!         THIS ROUTINE IS CALLED FROM *CUMASTR*.
    !!
    !!    PARAMETER     DESCRIPTION                                   UNITS
    !!    ---------     -----------                                   -----
    !!    INPUT PARAMETERS (INTEGER):

    !!   *KIDIA*        START POINT
    !!   *KFDIA*        END POINT
    !!   *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KTDIA*        START OF THE VERTICAL LOOP
    !!   *KLEV*         NUMBER OF LEVELS
    !!   *KCBOT*        CLOUD BASE LEVEL
    !!   *KCTOP*        CLOUD TOP LEVEL
    !!   *KDTOP*        TOP LEVEL OF DOWNDRAFTS

    !!   INPUT PARAMETERS (LOGICAL):

    !!   *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
    !!   *LDLAKE*       LAKE MASK (.TRUE. FOR LAKE)
    !!   *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

    !!   INPUT PARAMETERS (REAL):

    !!   *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
    !!   *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !!   *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
    !!   *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
    !!   *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
    !!   *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
    !!   *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
    !!   *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS            PA
    !!   *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2

    !!   UPDATED PARAMETERS (INTEGER):

    !!   *KTYPE*        SET TO ZERO IF LDCUM=.FALSE.

    !!   UPDATED PARAMETERS (LOGICAL):

    !!   *LDDRAF*       SET TO .FALSE. IF LDCUM=.FALSE. OR KDTOP<KCTOP

    !!   UPDATED PARAMETERS (REAL):

    !!   *PMFU*         MASSFLUX IN UPDRAFTS                          KG/(M2*S)
    !!   *PMFD*         MASSFLUX IN DOWNDRAFTS                        KG/(M2*S)
    !!   *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
    !!   *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
    !!   *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
    !!   *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
    !!   *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
    !!   *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
    !    *PSNDE*        DETRAINED SNOW/RAIN                           KG/(M2*S)
    !!   *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
    !!   *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS      KG/(M2*S)
    !    *PMFUDE_RATE*  UPRAFT DETRAINMENT RATE                       KG/(M2*S)
    !!   *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M2*S)

    !!   OUTPUT PARAMETERS (REAL):

    !!   *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)
    !!   *PLGLAC*       FLUX OF FROZEN CLOUD WATER IN UPDRAFTS        KG/(M2*S)
    !!   *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
    !!   *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
    !!   *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
    !!                  (NO EVAPORATION IN DOWNDRAFTS)

    !!         EXTERNALS
    !!         ---------
    !!         NONE

    !!         MODIFICATIONS
    !!         -------------
    !!            99-06-14 : Optimisation        D.SALMOND
    !!            03-08-28 : Clean up LINPHYS    P.BECHTOLD
    !!                       Bugs in Evapor.
    !!            05-02-11 : Extend DDMflux to   P.BECHTOLD
    !!                       surface if zero
    !!       M.Hamrud      01-Oct-2003 CY28 Cleaning

    !----------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
    !USE YOMCST   , ONLY : RG       ,RCPD     ,RLVTT    ,RLSTT    ,RLMLT    ,RTT

    !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
    ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
    ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
    ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
    !USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC
    !USE YOECUMF  , ONLY : RCUCOV   ,RCPECONS ,RTAUMEL  ,RHEBC, RMFCFL

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    REAL(KIND=jprb)   ,INTENT(in)    :: rmfcfl
    REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy
    REAL(KIND=jprb)   ,INTENT(in)    :: rhebc_land, rhebc_ocean
    REAL(KIND=jprb)   ,INTENT(in)    :: rhebc_land_trop, rhebc_ocean_trop
    LOGICAL           ,INTENT(in)    :: lmfdsnow
    REAL(KIND=jprb)   ,INTENT(in)    :: rcucov, rcucov_trop
    REAL(KIND=jprb)   ,INTENT(in)    :: trop_mask(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pqsen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: pap(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    LOGICAL           ,INTENT(in)    :: ldland(klon)
    LOGICAL           ,INTENT(in)    :: ldlake(klon)
    LOGICAL           ,INTENT(in)    :: ldcum(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kcbot(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kctop(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kdtop(klon)
    INTEGER(KIND=jpim),INTENT(out)   :: ktopm2
    INTEGER(KIND=jpim),INTENT(inout) :: ktype(klon)
    LOGICAL ,INTENT(inout) :: lddraf(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfus(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfds(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfuq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfdq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmful(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: plude(klon,klev)
    REAL(KIND=JPRB)   ,INTENT(in)    :: plrain(klon,klev) 
    REAL(KIND=JPRB)   ,INTENT(inout) :: psnde(klon,klev,2)
    REAL(KIND=jprb)   ,INTENT(inout) :: pdmfup(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pdmfdp(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pdpmel(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: plglac(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmflxr(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmflxs(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(out)   :: prain(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfude_rate(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfdde_rate(klon,klev)

    REAL(KIND=jprb) :: zrhebc(klon), zrcucov(klon), zshalfac(klon)
    INTEGER(KIND=jpim) :: ik, ikb, jk, jl
    INTEGER(KIND=jpim) :: idbas(klon)
    LOGICAL :: llddraf

    REAL(KIND=jprb) :: zalfaw, zcons1, zcons1a, zcons2, &
      & zdenom, zdrfl, zdrfl1, zfac, & !, zfoeewi, zfoeewl, &
      & zoealfa, zoeewm, zoelhm, zpdr, zpds, &
      & zrfl, zrfln, zrmin, zrnew, zsnmlt, ztarg, &
      & ztmst, zzp, zten, zwetb , zglac
    REAL(KIND=jprb) :: zhook_handle

    ! Numerical fit to wet bulb temperature
    REAL(KIND=JPRB),PARAMETER :: ZTW1 = 1329.31_JPRB
    REAL(KIND=JPRB),PARAMETER :: ZTW2 = 0.0074615_JPRB
    REAL(KIND=JPRB),PARAMETER :: ZTW3 = 0.85E5_JPRB
    REAL(KIND=JPRB),PARAMETER :: ZTW4 = 40.637_JPRB
    REAL(KIND=JPRB),PARAMETER :: ZTW5 = 275.0_JPRB

    !#include "fcttre.h"

    !--------------------------------------------------------------------
    !*             SPECIFY CONSTANTS

    IF (lhook) CALL dr_hook('CUFLXN',0,zhook_handle)
    ztmst=ptsphy
    zcons1a=rcpd/(rlmlt*rg*rtaumel)
    !ZCONS2=1.0_JPRB/(RG*ZTMST)
    zcons2=rmfcfl/(rg*ztmst)
    IF(lmfwetb) THEN
      zwetb=1.0_JPRB
    ELSE
      zwetb=0.0_JPRB
    ENDIF
    IF(lmfglac) THEN
      zglac=0.5_JPRB
    ELSE
      zglac=1.0_JPRB
    ENDIF

    !*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
    !!                 ---------------------------------

    DO jl=kidia,kfdia
      prain(jl)=0.0_JPRB
      IF(.NOT.ldcum(jl).OR.kdtop(jl) < kctop(jl)) lddraf(jl)=.FALSE.
      IF(.NOT.ldcum(jl)) ktype(jl)=0
      idbas(jl)=klev
      IF(ldland(jl) .OR. ldlake(jl)) THEN
        zrhebc(jl) = rhebc_land*(1._jprb - trop_mask(jl))  + rhebc_land_trop*trop_mask(jl)
      ELSE
        zrhebc(jl) = rhebc_ocean*(1._jprb - trop_mask(jl)) + rhebc_ocean_trop*trop_mask(jl)
      ENDIF
      zrcucov(jl) = rcucov*(1._jprb - trop_mask(jl)) + rcucov_trop*trop_mask(jl)
    ENDDO
    !!TO GET IDENTICAL RESULTS FOR DIFFERENT NPROMA FORCE KTOPM2 TO 2
    ktopm2=2
    DO jk=ktdia-1+ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      ikb=MIN(jk+1,klev)
      DO jl=kidia,kfdia
        pmflxr(jl,jk)=0.0_JPRB
        pmflxs(jl,jk)=0.0_JPRB
        pdpmel(jl,jk)=0.0_JPRB
        psnde(jl,jk,1) = 0.0_JPRB
        psnde(jl,jk,2) = 0.0_JPRB
        IF(ldcum(jl).AND.jk >= kctop(jl)) THEN
          pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*(rcpd*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
          plglac(jl,jk)=pmfu(jl,jk)*plglac(jl,jk)
          llddraf=lddraf(jl).AND.jk >= kdtop(jl)
          IF(llddraf) THEN
            pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*&
              & (rcpd*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
          ELSE
            pmfd(jl,jk)=0.0_JPRB
            pmfds(jl,jk)=0.0_JPRB
            pmfdq(jl,jk)=0.0_JPRB
            pdmfdp(jl,jk-1)=0.0_JPRB
          ENDIF
          IF(llddraf.AND.pmfd(jl,jk)<0._jprb .AND. pmfd(jl,ikb)==0._jprb) THEN
            idbas(jl)=jk
          ENDIF
        ELSE
          pmfu(jl,jk)=0.0_JPRB
          pmfd(jl,jk)=0.0_JPRB
          pmfus(jl,jk)=0.0_JPRB
          pmfds(jl,jk)=0.0_JPRB
          pmfuq(jl,jk)=0.0_JPRB
          pmfdq(jl,jk)=0.0_JPRB
          pmful(jl,jk)=0.0_JPRB
          plglac(jl,jk)=0.0_JPRB
          pdmfup(jl,jk-1)=0.0_JPRB
          pdmfdp(jl,jk-1)=0.0_JPRB
          plude(jl,jk-1)=0.0_JPRB
        ENDIF
      ENDDO
   ENDDO

   pmflxr(:,klev+1)=0.0_JPRB
   pmflxs(:,klev+1)=0.0_JPRB
   psnde(:,ktdia,1)=0.0_JPRB
   psnde(:,ktdia,2)=0.0_JPRB


    !*    1.5          SCALE FLUXES BELOW CLOUD BASE
    !!                 LINEAR DCREASE
    !!                 -----------------------------

!DIR$ IVDEP
!OCL NOVREC
    DO jl=kidia,kfdia
      IF(ldcum(jl)) THEN
        ikb=kcbot(jl)
        ik=ikb+1
        zzp=((paph(jl,klev+1)-paph(jl,ik))/(paph(jl,klev+1)-paph(jl,ikb)))
        IF(ktype(jl) == 3) zzp=zzp*zzp
        pmfu(jl,ik)=pmfu(jl,ikb)*zzp
        IF (lphylin) THEN
          ztarg=ptenh(jl,ikb)
          zoealfa=0.545_JPRB*(TANH(0.17_JPRB*(ztarg-rlptrc))+1.0_JPRB)
          zoelhm=zoealfa*rlvtt+(1.0_JPRB-zoealfa)*rlstt
          pmfus(jl,ik)=(pmfus(jl,ikb)-zoelhm*pmful(jl,ikb))*zzp
        ELSE
          pmfus(jl,ik)=(pmfus(jl,ikb)-foelhmcu(ptenh(jl,ikb))*pmful(jl,ikb))*zzp
        ENDIF
        pmfuq(jl,ik)=(pmfuq(jl,ikb)+pmful(jl,ikb))*zzp
        pmful(jl,ik)=0.0_JPRB
      ENDIF
    ENDDO
    DO jk=ktdia-1+ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=kidia,kfdia
        IF(ldcum(jl).AND.jk > kcbot(jl)+1) THEN
          ikb=kcbot(jl)+1
          zzp=((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
          IF(ktype(jl) == 3) zzp=zzp*zzp
          pmfu(jl,jk)=pmfu(jl,ikb)*zzp
          pmfus(jl,jk)=pmfus(jl,ikb)*zzp
          pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
          pmful(jl,jk)=0.0_JPRB
        ENDIF
        ik=idbas(jl)
        llddraf=lddraf(jl).AND.jk>ik.AND.ik<klev
        IF(llddraf.AND.ik==kcbot(jl)+1) THEN
          zzp=((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ik)))
          IF(ktype(jl) == 3)  zzp=zzp*zzp
          pmfd(jl,jk)=pmfd(jl,ik)*zzp
          pmfds(jl,jk)=pmfds(jl,ik)*zzp
          pmfdq(jl,jk)=pmfdq(jl,ik)*zzp
          pmfdde_rate(jl,jk)=-(pmfd(jl,jk-1)-pmfd(jl,jk))
        ELSEIF(llddraf.AND.ik/=kcbot(jl)+1.AND.jk==ik+1) THEN
          pmfdde_rate(jl,jk)=-(pmfd(jl,jk-1)-pmfd(jl,jk))
        ENDIF
      ENDDO
    ENDDO

    !*    2.            CALCULATE RAIN/SNOW FALL RATES
    !*                  CALCULATE MELTING OF SNOW
    !*                  CALCULATE EVAPORATION OF PRECIP
    !!                  -------------------------------

    DO jk=ktdia-1+ktopm2,klev
      DO jl=kidia,kfdia
        IF(ldcum(jl).AND.jk >= kctop(jl)-1) THEN
          prain(jl)=prain(jl)+pdmfup(jl,jk)
          zten=pten(jl,jk)-zwetb*MAX(0.0_JPRB,pqsen(jl,jk)-pqen(jl,jk))*&
            (ztw1+ztw2*(pap(jl,jk)-ztw3)-ztw4*(pten(jl,jk)-ztw5))
          IF(pmflxs(jl,jk) > 0.0_JPRB.AND.zten > rtt) THEN
            zcons1=zcons1a*(1.0_JPRB+0.5_JPRB*(zten-rtt))
            zfac=zcons1*(paph(jl,jk+1)-paph(jl,jk))
            zsnmlt=MIN(pmflxs(jl,jk),zfac*(zten-rtt))
            pdpmel(jl,jk)=zsnmlt
            IF (lphylin) THEN
              ztarg=zten-zsnmlt/zfac
              zoealfa=0.545_JPRB*(TANH(0.17_JPRB*(ztarg-rlptrc))+1.0_JPRB)
              !>KF use functions!
              zoeewm=zoealfa*foeewl(ztarg) + (1.0_JPRB-zoealfa)*foeewi(ztarg)
              !<KF
              pqsen(jl,jk)=zoeewm/pap(jl,jk)
            ELSE
              pqsen(jl,jk)=foeewmcu(zten-zsnmlt/zfac)/pap(jl,jk)
            ENDIF
          ENDIF
          IF (lphylin) THEN
            zalfaw=0.545_JPRB*(TANH(0.17_JPRB*(zten-rlptrc))+1.0_JPRB)
          ELSE
            zalfaw=foealfcu(zten)
         ENDIF
          !! no liquid precipitation above melting level
          IF(pten(jl,jk) < rtt .AND. zalfaw > 0.0_JPRB) THEN
            plglac(jl,jk)=plglac(jl,jk)+zglac*zalfaw*pdmfup(jl,jk)+zalfaw*pdmfdp(jl,jk)
            zalfaw=0.0_JPRB
         ENDIF
         IF (lmfdsnow .AND. jk<=kcbot(jl)) THEN
           psnde(jl,jk,2)=plrain(jl,jk+1)*pmfude_rate(jl,jk)
           psnde(jl,jk,1)=psnde(jl,jk,2)*(1.0_JPRB-zalfaw)
           psnde(jl,jk,2)=psnde(jl,jk,2)*zalfaw
           psnde(jl,jk,1)=MIN(psnde(jl,jk,1),(1.0_JPRB-zalfaw)*(pdmfup(jl,jk)+pdmfdp(jl,jk))-pdpmel(jl,jk))
           psnde(jl,jk,2)=MIN(psnde(jl,jk,2),zalfaw*(pdmfup(jl,jk)+pdmfdp(jl,jk)))
           psnde(jl,jk,1)=MAX(psnde(jl,jk,1),0.0_JPRB)
           psnde(jl,jk,2)=MAX(psnde(jl,jk,2),0.0_JPRB)
         ENDIF
         pmflxr(jl,jk+1)=pmflxr(jl,jk)+zalfaw*&
            & (pdmfup(jl,jk)+pdmfdp(jl,jk))+pdpmel(jl,jk)-psnde(jl,jk,2)
          pmflxs(jl,jk+1)=pmflxs(jl,jk)+(1.0_JPRB-zalfaw)*&
            & (pdmfup(jl,jk)+pdmfdp(jl,jk))-pdpmel(jl,jk)-psnde(jl,jk,1)
          IF(pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < 0.0_JPRB) THEN
            pdmfdp(jl,jk)=-(pmflxr(jl,jk)+pmflxs(jl,jk)+pdmfup(jl,jk))
            pmflxr(jl,jk+1)=0.0_JPRB
            pmflxs(jl,jk+1)=0.0_JPRB
            pdpmel(jl,jk)  =0.0_JPRB
          ELSEIF(pmflxr(jl,jk+1) < 0.0_JPRB ) THEN
            pmflxs(jl,jk+1)=pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
            pmflxr(jl,jk+1)=0.0_JPRB
          ELSEIF(pmflxs(jl,jk+1) < 0.0_JPRB ) THEN
            pmflxr(jl,jk+1)=pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
            pmflxs(jl,jk+1)=0.0_JPRB
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !!Reminder for conservation:
    !!   pdmfup(jl,jk)+pdmfdp(jl,jk)=pmflxr(jl,jk+1)+pmflxs(jl,jk+1)-pmflxr(jl,jk)-pmflxs(jl,jk)

    DO jl=kidia,kfdia
      IF (ktype(jl)==2) THEN
        zshalfac(jl) = 0.5_jprb
      ELSE
        zshalfac(jl) = 1._jprb
      ENDIF
    ENDDO
    DO jk=ktdia-1+ktopm2,klev
      DO jl=kidia,kfdia
        IF(ldcum(jl).AND.jk >= kcbot(jl)) THEN
          zrfl=pmflxr(jl,jk)+pmflxs(jl,jk)
          IF(zrfl > 1.e-20_JPRB) THEN
            zdrfl1=rcpecons*MAX(0.0_JPRB,pqsen(jl,jk)-pqen(jl,jk))*zrcucov(jl)*&
              & EXP(0.5777_JPRB*LOG(SQRT(paph(jl,jk)/paph(jl,klev+1))/5.09E-3_JPRB*zrfl/zrcucov(jl)))*&
              & (paph(jl,jk+1)-paph(jl,jk))
            zrnew=zrfl-zdrfl1
            zrmin=zrfl-zrcucov(jl)*MAX(0.0_JPRB,zrhebc(jl)*pqsen(jl,jk)-pqen(jl,jk))&
              & *zshalfac(jl)*zcons2*(paph(jl,jk+1)-paph(jl,jk))
            zrnew=MAX(zrnew,zrmin)
            zrfln=MAX(zrnew,0.0_JPRB)
            zdrfl=MIN(0.0_JPRB,zrfln-zrfl)
            IF (lphylin) THEN
              zalfaw=0.545_JPRB*(TANH(0.17_JPRB*(pten(jl,jk)-rlptrc))+1.0_JPRB)
            ELSE
              zalfaw=foealfcu(pten(jl,jk))
            ENDIF
            IF(pten(jl,jk) < rtt) zalfaw=0.0_JPRB
            zpdr=zalfaw*pdmfdp(jl,jk)
            zpds=(1.0_JPRB-zalfaw)*pdmfdp(jl,jk)
            zdenom=1.0_JPRB/MAX(1.e-20_JPRB,pmflxr(jl,jk)+pmflxs(jl,jk))
            pmflxr(jl,jk+1)=pmflxr(jl,jk)+zpdr &
              & +pdpmel(jl,jk)+zdrfl*pmflxr(jl,jk)*zdenom
            pmflxs(jl,jk+1)=pmflxs(jl,jk)+zpds &
              & -pdpmel(jl,jk)+zdrfl*pmflxs(jl,jk)*zdenom
            pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
            IF(pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < 0.0_JPRB) THEN
              pdmfup(jl,jk)=pdmfup(jl,jk)-(pmflxr(jl,jk+1)+pmflxs(jl,jk+1))
              pmflxr(jl,jk+1)=0.0_JPRB
              pmflxs(jl,jk+1)=0.0_JPRB
              pdpmel(jl,jk)  =0.0_JPRB
            ELSEIF(pmflxr(jl,jk+1) < 0.0_JPRB ) THEN
              pmflxs(jl,jk+1)=pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
              pmflxr(jl,jk+1)=0.0_JPRB
            ELSEIF(pmflxs(jl,jk+1) < 0.0_JPRB ) THEN
              pmflxr(jl,jk+1)=pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
              pmflxs(jl,jk+1)=0.0_JPRB
            ENDIF
          ELSE
            pmflxr(jl,jk+1)=0.0_JPRB
            pmflxs(jl,jk+1)=0.0_JPRB
            pdmfdp(jl,jk)  =0.0_JPRB
            pdpmel(jl,jk)  =0.0_JPRB
          ENDIF
        ENDIF
      ENDDO
    ENDDO

     IF (lhook) CALL dr_hook('CUFLXN',1,zhook_handle)
  END SUBROUTINE cuflxn

  !=======================================================================

  SUBROUTINE cudtdqn &
    & (  kidia,    kfdia,    klon,  ktdia,    klev,&
    & ktopm2,   ktype,    kctop,    kdtop,    ldcum,    lddraf,   ptsphy,&
    & paph,     pgeoh,    pgeo,&
    & zdph,                    &
    & pten,     ptenh,    pqen,     pqenh,    pqsen,&
    & plglac,   plude,    psnde,    pmfu,     pmfd,&
    & pmfus,    pmfds,    pmfuq,    pmfdq,&
    & pmful,    pdmfup,   pdpmel,&
    & ptent,    ptenq,    penth )
    !
    !!Description:
    !**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
    !!               DOES GLOBAL DIAGNOSTICS

    !!         M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
    !!         P.BECHTOLD        E.C.M.W.F.     10/05


    !**   INTERFACE.
    !!    ----------

    !!         *CUDTDQ* IS CALLED FROM *CUMASTR*

    !>    PARAMETER     DESCRIPTION                                   UNITS
    !!    ---------     -----------                                   -----
    !!    INPUT PARAMETERS (INTEGER):

    !!   *KIDIA*        START POINT
    !!   *KFDIA*        END POINT
    !!   *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KLEV*         NUMBER OF LEVELS
    !!   *KTYPE*        TYPE OF CONVECTION
    !!                      1 = PENETRATIVE CONVECTION
    !!                      2 = SHALLOW CONVECTION
    !!                      3 = MIDLEVEL CONVECTION
    !!   *KCTOP*        CLOUD TOP LEVEL
    !!   *KDTOP*        TOP LEVEL OF DOWNDRAFTS

    !!   INPUT PARAMETERS (LOGICAL):

    !!   *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
    !!   *LDDRAF*       FLAG: .TRUE. FOR DOWNDRAFT LEVEL

    !!   INPUT PARAMETERS (REAL):

    !!   *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
    !!   *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
    !!   *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
    !!   *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2
    !!    *zdph*         pressure thickness on half levels              Pa
    !!   *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !!   *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
    !!   *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
    !!   *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
    !!   *PQSEN*        SATURATION ENV. SPEC. HUMIDITY (T+1)          KG/KG
    !!   *PLGLAC*       FLUX OF FROZEN CLOUDWATER IN UPDRAFTS         KG/(M2*S)
    !!   *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
    !    *PSNDE*        DETRAINED SNOW/RAIN                           KG/(M2*S)
    !!   *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
    !!   *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
    !!   *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
    !!   *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
    !!   *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
    !!   *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
    !!   *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
    !!   *PDMFUP*       FLUX DIFFERENCE OF PRECIP.                    KG/(M2*S)
    !!   *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)

    !!   UPDATED PARAMETERS (REAL):

    !!   *PTENT*        TEMPERATURE TENDENCY                           K/S
    !!   *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)

    !!   OUTPUT PARAMETERS (REAL):

    !!   *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)

    !----------------------------------------------------------------------

    !!              MODIFICATIONS
    !!       M.Hamrud      01-Oct-2003 CY28 Cleaning

    !!      96-09-20       : changed so can be used with diagnost
    !!                       cloud scheme      D.GREGORY
    !!      99-06-04       : Optimisation      D.SALMOND
    !!      03-08-28       : Clean up LINPHYS  P.BECHTOLD
    !!      05-10-13       : implicit solution P.BECHTOLD

    !=======================================================================

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RG       ,RCPD     ,RLVTT    ,RLSTT    ,RLMLT    ,RTT
    !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
    ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
    ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
    ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
    !USE YOECUMF  , ONLY : RMFSOLTQ
    !USE YOEPHY   , ONLY : LEPCLD
    !USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC
    !USE YOPHNC   , ONLY : LENCLD2

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktopm2
    INTEGER(KIND=jpim),INTENT(in)    :: ktype(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kctop(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kdtop(klon)
    LOGICAL ,INTENT(inout) :: ldcum(klon)
    LOGICAL ,INTENT(in)    :: lddraf(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: zdph(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqsen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: plglac(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: plude(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: psnde(klon,klev,2)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfus(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfds(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfuq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfdq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmful(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pdmfup(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pdpmel(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptent(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptenq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(out)   :: penth(klon,klev)

    LOGICAL :: lltest
    INTEGER(KIND=jpim) :: jk, ik, jl
    REAL(KIND=jprb)    :: ztsphy, zimp, zorcpd, zalv, zoealfa, ztarg,&
      & zzp, zgq, zgs, zgh, zs, zq, zcp_o_cv

    REAL(KIND=jprb) :: zmfus(klon,klev), zmfuq(klon,klev),&
      & zmfds(klon,klev), zmfdq(klon,klev)


    REAL(KIND=jprb),   DIMENSION(klon,klev) :: zdtdt, zdqdt , ZDP
    REAL(KIND=jprb),   DIMENSION(klon,klev) :: zb,    zr1,   zr2
    LOGICAL,           DIMENSION(klon,klev) :: llcumbas
    REAL(KIND=jprb) :: zhook_handle

    !#include "cubidiag.intfb.h"
    !#include "fcttre.h"

    !----------------------------------------------------------------------
    IF (lhook) CALL dr_hook('CUDTDQN',0,zhook_handle)


    !*    1.0          SETUP AND INITIALIZATIONS
    !!                 -------------------------

    zimp=1.0_JPRB-rmfsoltq
    ztsphy=1.0_JPRB/ptsphy
    zorcpd=1.0_JPRB/rcpd

    DO jk=ktdia,klev
      DO jl=kidia,kfdia
        penth(jl,jk)=0.0_JPRB
      ENDDO
    ENDDO

    !!        MASS-FLUX APPROACH SWITCHED ON FOR DEEP CONVECTION ONLY
    !!        IN THE TANGENT-LINEAR AND ADJOINT VERSIONS

    DO jl=kidia,kfdia
      IF (ktype(jl) /= 1.AND.lphylin) ldcum(jl)=.FALSE.
    ENDDO

    !!zero detrained liquid water if diagnostic cloud scheme to be used

    !!this means that detrained liquid water will be evaporated in the
    !!cloud environment and not fed directly into a cloud liquid water
    !!variable

    !LLTEST = (.NOT.LEPCLD.AND..NOT.LENCLD2).OR.(LPHYLIN.AND..NOT.LENCLD2)
    lltest = .NOT.lepcld.OR.lphylin

    IF (lltest) THEN
      DO jk=ktdia,klev
        DO jl=kidia,kfdia
          plude(jl,jk)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF

    DO jk=ktdia,klev
      DO jl=kidia,kfdia
        IF(ldcum(jl)) THEN
!           ZDP(JL,JK)=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
            zdp(jl,jk)=rg/(paph(jl,jk+1)-paph(jl,jk))
          zmfus(jl,jk)=pmfus(jl,jk)
          zmfds(jl,jk)=pmfds(jl,jk)
          zmfuq(jl,jk)=pmfuq(jl,jk)
          zmfdq(jl,jk)=pmfdq(jl,jk)
        ENDIF
      ENDDO
    ENDDO

    !------------------------------------------------------------------------------

    IF ( rmfsoltq>0.0_JPRB ) THEN

      !*    2.0          RECOMPUTE CONVECTIVE FLUXES IF IMPLICIT

      DO jk=ktdia-1+ktopm2,klev
        ik=jk-1
!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          IF(ldcum(jl).AND.jk>=kctop(jl)-1) THEN
            ! compute interpolating coefficients ZGS and ZGQ for half-level values
            zgq=(pqenh(jl,jk)-pqen(jl,ik))/pqsen(jl,jk)
            zgh =rcpd*pten(jl,jk)+pgeo(jl,jk)
            zgs=(rcpd*(ptenh(jl,jk)-pten(jl,ik))+pgeoh(jl,jk)-pgeo(jl,ik))/zgh

            !half-level environmental values for S and Q
            zs =rcpd*(zimp*pten(jl,ik)+zgs*pten(jl,jk))+pgeo(jl,ik)+zgs*pgeo(jl,jk)
            zq =zimp*pqen(jl,ik)+zgq*pqsen(jl,jk)
            zmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*zs
            zmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*zq
            IF(lddraf(jl).AND.jk >= kdtop(jl)) THEN
              zmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*zs
              zmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*zq
            ENDIF
          ENDIF
        ENDDO
      ENDDO

    ENDIF

    !*    3.0          COMPUTE TENDENCIES
    !                  ------------------

    DO jk=ktdia-1+ktopm2,klev

      IF(jk < klev) THEN
        DO jl=kidia,kfdia
          IF(ldcum(jl)) THEN
            IF (lphylin) THEN
              ztarg=pten(jl,jk)
              zoealfa=0.545_JPRB*(TANH(0.17_JPRB*(ztarg-rlptrc))+1.0_JPRB)
              zalv=zoealfa*rlvtt+(1.0_JPRB-zoealfa)*rlstt
            ELSE
              zalv=foelhmcu(pten(jl,jk))
            ENDIF
            !>KF
                    ZDTDT(JL,JK)=ZDP(JL,JK)*ZORCPD*&
            !zdtdt(jl,jk)=rg/zdph(jl,jk)*zorcpd*&
              & (zmfus(jl,jk+1)-zmfus(jl,jk)+&
              & zmfds(jl,jk+1)-zmfds(jl,jk)&
              & +(rlmlt*plglac(jl,jk)&
              & -rlmlt*pdpmel(jl,jk)&
              & -zalv*(pmful(jl,jk+1)-pmful(jl,jk)-&
              & plude(jl,jk)-pdmfup(jl,jk)-psnde(jl,jk,1)-psnde(jl,jk,2))))
            !>KF
                    ZDQDT(JL,JK)=ZDP(JL,JK)*&
            !zdqdt(jl,jk)=rg/zdph(jl,jk)*&
              & (zmfuq(jl,jk+1)-zmfuq(jl,jk)+&
              & zmfdq(jl,jk+1)-zmfdq(jl,jk)+&
              & pmful(jl,jk+1)-pmful(jl,jk)-&
              & plude(jl,jk)-pdmfup(jl,jk)-psnde(jl,jk,1)-psnde(jl,jk,2))
          ENDIF
       ENDDO

      ELSE
        DO jl=kidia,kfdia
          IF(ldcum(jl)) THEN
            IF (lphylin) THEN
              ztarg=pten(jl,jk)
              zoealfa=0.545_JPRB*(TANH(0.17_JPRB*(ztarg-rlptrc))+1.0_JPRB)
              zalv=zoealfa*rlvtt+(1.0_JPRB-zoealfa)*rlstt
            ELSE
              zalv=foelhmcu(pten(jl,jk))
            ENDIF
            !>KF
            ZDTDT(JL,JK)=-ZDP(JL,JK)*ZORCPD*&
            !zdtdt(jl,jk)=-rg/zdph(jl,jk)*zorcpd*&
              & (zmfus(jl,jk)+zmfds(jl,jk)+(rlmlt*pdpmel(jl,jk)-zalv*&
              & (pmful(jl,jk)+pdmfup(jl,jk))))

            !>KF
            ZDQDT(JL,JK)=-ZDP(JL,JK)*&
            !zdqdt(jl,jk)=-rg/zdph(jl,jk)*&
              & (zmfuq(jl,jk)+zmfdq(jl,jk)+&
              & (pmful(jl,jk)+pdmfup(jl,jk)))
          ENDIF
        ENDDO
      ENDIF

    ENDDO

    IF ( rmfsoltq==0.0_JPRB ) THEN

      !*    3.1          UPDATE TENDENCIES
      !!                 -----------------

      DO jk=ktdia-1+ktopm2,klev
        DO jl=kidia,kfdia
          IF(ldcum(jl)) THEN
            ptent(jl,jk)=ptent(jl,jk)+zdtdt(jl,jk)
            ptenq(jl,jk)=ptenq(jl,jk)+zdqdt(jl,jk)
            penth(jl,jk)=zdtdt(jl,jk)*rcpd
          ENDIF
        ENDDO
      ENDDO

    ELSE
      !----------------------------------------------------------------------

      !*    3.2          IMPLICIT SOLUTION
      !!                 -----------------

      !!Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1;
      !!reuse ZMFUS=A
      !!ZDTDT and ZDQDT correspond to the RHS ("constants") of the equation
      !!The solution is in ZR1 and ZR2

      llcumbas(:,:)=.FALSE.
      zb(:,:)=1.0_JPRB
      zmfus(:,:)=0.0_JPRB

      !!Fill vectors A, B and RHS

      DO jk=ktdia-1+ktopm2,klev
        ik=jk+1
        DO jl=kidia,kfdia
          llcumbas(jl,jk)=ldcum(jl).AND.jk>=kctop(jl)-1
          IF(llcumbas(jl,jk)) THEN
            !>KF
             ZZP=RMFSOLTQ*ZDP(JL,JK)*PTSPHY
            !zzp=rmfsoltq*rg/zdph(jl,jk)*ptsphy
            zmfus(jl,jk)=-zzp*(pmfu(jl,jk)+pmfd(jl,jk))
            zdtdt(jl,jk) = zdtdt(jl,jk)*ptsphy+pten(jl,jk)
            zdqdt(jl,jk) = zdqdt(jl,jk)*ptsphy+pqen(jl,jk)
            ! ZDTDT(JL,JK) = (ZDTDT(JL,JK)+PTENT(JL,JK))*PTSPHY+PTEN(JL,JK)
            ! ZDQDT(JL,JK) = (ZDQDT(JL,JK)+PTENQ(JL,JK))*PTSPHY+PQEN(JL,JK)
            IF(jk<klev) THEN
              zb(jl,jk)=1.0_JPRB+zzp*(pmfu(jl,ik)+pmfd(jl,ik))
            ELSE
              zb(jl,jk)=1.0_JPRB
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      CALL cubidiag&
        & ( kidia,    kfdia,   klon,   ktdia, klev, &
        & kctop,    llcumbas, &
        & zmfus,    zb,     zdtdt,   zr1 )

      CALL cubidiag&
        & ( kidia,    kfdia,   klon,   ktdia, klev, &
        & kctop,    llcumbas, &
        & zmfus,    zb,     zdqdt,   zr2 )

      ! Compute tendencies

!PREVENT_INCONSISTENT_IFORT_FMA
      DO jk=ktdia-1+ktopm2,klev
        DO jl=kidia,kfdia
          IF(llcumbas(jl,jk)) THEN
            ptent(jl,jk)=ptent(jl,jk)+(zr1(jl,jk)-pten(jl,jk))*ztsphy
            ptenq(jl,jk)=ptenq(jl,jk)+(zr2(jl,jk)-pqen(jl,jk))*ztsphy
            ! PTENT(JL,JK)=(ZR1(JL,JK)-PTEN(JL,JK))*ZTSPHY
            ! PTENQ(JL,JK)=(ZR2(JL,JK)-PQEN(JL,JK))*ZTSPHY
            penth(jl,jk)=(zr1(jl,jk)-pten(jl,jk))*ztsphy*rcpd
          ENDIF
        ENDDO
      ENDDO

      !----------------------------------------------------------------------
    ENDIF

    IF (lhook) CALL dr_hook('CUDTDQN',1,zhook_handle)
  END SUBROUTINE cudtdqn

  !=======================================================================

  SUBROUTINE cududv &
    & ( kidia,    kfdia,    klon,     ktdia,    klev,&
    & ktopm2,   ktype,    kcbot,    kctop,    ldcum,    ptsphy,&
    & zdph,                                         &
    & paph,     puen,     pven,     pmfu,     pmfd,&
    & puu,      pud,      pvu,      pvd,&
    & ptenu,    ptenv  )
    !
    !!Description:
    !**** *CUDUDV* - UPDATES U AND V TENDENCIES,
    !!               DOES GLOBAL DIAGNOSTIC OF DISSIPATION

    !!         M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
    !!         P.BECHTOLD        E.C.M.W.F.    11/02/05 IMPLICIT SOLVER

    !**   INTERFACE.
    !!    ----------

    !!         *CUDUDV* IS CALLED FROM *CUMASTR*

    !!    PARAMETER     DESCRIPTION                                   UNITS
    !!    ---------     -----------                                   -----
    !!    INPUT PARAMETERS (INTEGER):

    !!   *KIDIA*        START POINT
    !!   *KFDIA*        END POINT
    !!   *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KTDIA*        START OF THE VERTICAL LOOP
    !!   *KLEV*         NUMBER OF LEVELS
    !!   *KTYPE*        TYPE OF CONVECTION
    !!                      1 = PENETRATIVE CONVECTION
    !!                      2 = SHALLOW CONVECTION
    !!                      3 = MIDLEVEL CONVECTION
    !!   *KCBOT*        CLOUD BASE LEVEL
    !!   *KCTOP*        CLOUD TOP LEVEL

    !!   INPUT PARAMETERS (LOGICAL):

    !!   *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

    !!   INPUT PARAMETERS (REAL):

    !!   *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
    !!   *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
    !!   *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
    !!   *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
    !!   *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
    !!   *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
    !!   *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
    !!   *PUD*          U-VELOCITY IN DOWNDRAFTS                       M/S
    !!   *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
    !!   *PVD*          V-VELOCITY IN DOWNDRAFTS                       M/S
    !!   *zdph*         pressure thickness on half levels              Pa

    !!   UPDATED PARAMETERS (REAL):

    !!   *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
    !!   *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2

    !!           METHOD
    !!           -------
    !!      EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
    !!      DEPENDING ON VALUE OF RMFSOLUV:
    !!      0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT
    !
    !!      FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
    !!      FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
    !!                             TO CORRECT TENDENCIES BELOW CLOUD BASE
    !
    !!           EXTERNALS
    !!           ---------
    !!           CUBIDIAG

    !!         MODIFICATIONS
    !!         -------------
    !!            92-09-21 : Update to Cy44      J.-J. MORCRETTE
    !!       M.Hamrud      01-Oct-2003 CY28 Cleaning
    !----------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
    !USE YOMCST        , ONLY : RG

    !USE YOECUMF       , ONLY : RMFSOLUV

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim)               :: ktdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktopm2
    INTEGER(KIND=jpim),INTENT(in)    :: ktype(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kcbot(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kctop(klon)
    LOGICAL ,INTENT(in)    :: ldcum(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: zdph(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: puen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pven(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: puu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pud(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pvu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pvd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptenu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptenv(klon,klev)
    REAL(KIND=jprb) :: zuen(klon,klev),     zven(klon,klev),&
      & zmfuu(klon,klev),    zmfdu(klon,klev),&
      & zmfuv(klon,klev),    zmfdv(klon,klev)

    INTEGER(KIND=jpim) :: ik, ikb, jk, jl

    REAL(KIND=jprb) :: zzp, zimp, ztsphy
    !       ALLOCATABLE ARAYS
    REAL(KIND=jprb),   DIMENSION(klon,klev) :: zdudt, zdvdt, ZDP
    REAL(KIND=jprb),   DIMENSION(klon,klev) :: zb,  zr1,  zr2
    LOGICAL,           DIMENSION(klon,klev) :: llcumbas
    REAL(KIND=jprb) :: zhook_handle

    !#include "cubidiag.intfb.h"
    !----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('CUDUDV',0,zhook_handle)
    zimp=1.0_JPRB-rmfsoluv
    ztsphy=1.0_JPRB/ptsphy

    zmfuu(:,:) = 0._jprb
    zmfdu(:,:) = 0._jprb
    zmfuv(:,:) = 0._jprb
    zmfdv(:,:) = 0._jprb

    DO jk=ktdia,klev
      DO jl=kidia,kfdia
        IF(ldcum(jl)) THEN
          zuen(jl,jk)=puen(jl,jk)
          zven(jl,jk)=pven(jl,jk)
          zdp(jl,jk)=rg/(paph(jl,jk+1)-paph(jl,jk))
        ENDIF
      ENDDO
    ENDDO

    !----------------------------------------------------------------------


    !*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
    !!                 ----------------------------------------------

    DO jk=ktdia-1+ktopm2,klev
      ik=jk-1
      DO jl=kidia,kfdia
        IF(ldcum(jl)) THEN
          zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-zimp*zuen(jl,ik))
          zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-zimp*zven(jl,ik))
          zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-zimp*zuen(jl,ik))
          zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-zimp*zven(jl,ik))
        ENDIF
      ENDDO
    ENDDO

    ! linear fluxes below cloud
    IF(rmfsoluv==0.0_JPRB) THEN
      DO jk=ktdia-1+ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          IF(ldcum(jl).AND.jk > kcbot(jl)) THEN
            ikb=kcbot(jl)
            zzp=((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            IF(ktype(jl) == 3) THEN
              zzp=zzp*zzp
            ENDIF
            zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
            zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
            zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
            zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    !*    1.2          COMPUTE TENDENCIES
    !                  ------------------

    DO jk=ktdia-1+ktopm2,klev

      IF(jk < klev) THEN
        ik=jk+1
        DO jl=kidia,kfdia
          IF(ldcum(jl)) THEN
            !>KF
             ZDUDT(JL,JK)=ZDP(JL,JK)*&
            !zdudt(jl,jk)=rg/zdph(jl,jk)*&
              & (zmfuu(jl,ik)-zmfuu(jl,jk)+zmfdu(jl,ik)-zmfdu(jl,jk))
            !>KF
             ZDVDT(JL,JK)=ZDP(JL,JK)*&
            !zdvdt(jl,jk)=rg/zdph(jl,jk)*&
              & (zmfuv(jl,ik)-zmfuv(jl,jk)+zmfdv(jl,ik)-zmfdv(jl,jk))
          ENDIF
        ENDDO

      ELSE
        DO jl=kidia,kfdia
          IF(ldcum(jl)) THEN
            !>KF
             ZDUDT(JL,JK)=-ZDP(JL,JK)*(ZMFUU(JL,JK)+ZMFDU(JL,JK))
             ZDVDT(JL,JK)=-ZDP(JL,JK)*(ZMFUV(JL,JK)+ZMFDV(JL,JK))
            !zdudt(jl,jk)=-rg/zdph(jl,jk)*(zmfuu(jl,jk)+zmfdu(jl,jk))
            !zdvdt(jl,jk)=-rg/zdph(jl,jk)*(zmfuv(jl,jk)+zmfdv(jl,jk))
          ENDIF
        ENDDO
      ENDIF

    ENDDO

    IF ( rmfsoluv==0.0_JPRB ) THEN

      !*    1.3          UPDATE TENDENCIES
      !                  -----------------

      DO jk=ktdia-1+ktopm2,klev
        DO jl=kidia,kfdia
          IF(ldcum(jl)) THEN
            ptenu(jl,jk)=ptenu(jl,jk)+zdudt(jl,jk)
            ptenv(jl,jk)=ptenv(jl,jk)+zdvdt(jl,jk)
          ENDIF
        ENDDO
      ENDDO

    ELSE
      !----------------------------------------------------------------------

      !*      1.6          IMPLICIT SOLUTION
      !                    -----------------

      !!Fill bi-diagonal Matrix vectors A=k-1, B=k;
      !!reuse ZMFUU=A and ZB=B;
      !!ZDUDT and ZDVDT correspond to the RHS ("constants") of the equation
      !!The solution is in ZR1 and ZR2

      llcumbas(:,:)=.FALSE.
      zb(:,:)=1.0_JPRB
      zmfuu(:,:)=0.0_JPRB

      ! Fill vectors A, B and RHS

      DO jk=ktdia-1+ktopm2,klev
        ik=jk+1
        DO jl=kidia,kfdia
          llcumbas(jl,jk)=ldcum(jl).AND.jk>=kctop(jl)-1
          IF(llcumbas(jl,jk)) THEN
            !>KF
             ZZP=RMFSOLUV*ZDP(JL,JK)*PTSPHY
            !zzp=rmfsoluv*rg/zdph(jl,jk)*ptsphy
            zmfuu(jl,jk)=-zzp*(pmfu(jl,jk)+pmfd(jl,jk))
            zdudt(jl,jk) = zdudt(jl,jk)*ptsphy+zuen(jl,jk)
            zdvdt(jl,jk) = zdvdt(jl,jk)*ptsphy+zven(jl,jk)
            ! ZDUDT(JL,JK) = (PTENU(JL,JK)+ZDUDT(JL,JK))*PTSPHY+ZUEN(JL,JK)
            ! ZDVDT(JL,JK) = (PTENV(JL,JK)+ZDVDT(JL,JK))*PTSPHY+ZVEN(JL,JK)
            IF(jk<klev) THEN
              zb(jl,jk)=1.0_JPRB+zzp*(pmfu(jl,ik)+pmfd(jl,ik))
            ELSE
              zb(jl,jk)=1.0_JPRB
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      CALL cubidiag&
        & ( kidia, kfdia, klon, ktdia, klev, &
        & kctop, llcumbas, &
        & zmfuu,    zb,    zdudt,   zr1 )

      CALL cubidiag&
        & ( kidia, kfdia, klon, ktdia, klev, &
        & kctop, llcumbas, &
        & zmfuu,    zb,    zdvdt,   zr2 )

      ! Compute tendencies

      DO jk=ktdia-1+ktopm2,klev
        DO jl=kidia,kfdia
          IF(llcumbas(jl,jk)) THEN
            ptenu(jl,jk)=ptenu(jl,jk)+(zr1(jl,jk)-zuen(jl,jk))*ztsphy
            ptenv(jl,jk)=ptenv(jl,jk)+(zr2(jl,jk)-zven(jl,jk))*ztsphy
            ! PTENU(JL,JK)=(ZR1(JL,JK)-ZUEN(JL,JK))*ZTSPHY
            ! PTENV(JL,JK)=(ZR2(JL,JK)-ZVEN(JL,JK))*ZTSPHY
          ENDIF
        ENDDO
      ENDDO

    ENDIF
    !----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('CUDUDV',1,zhook_handle)
  END SUBROUTINE cududv

  !=======================================================================

  SUBROUTINE cuctracer &
    & ( kidia,    kfdia,   klon,   ktdia, klev, ktrac,&
    & kctop,     kdtop,   &
    & ldcum,    lddraf,   ptsphy,   &
    & paph,     zdph,                       &
    & pmfu,     pmfd,     pudrate,  pddrate,&
    & pcen,     ptenc  )
    !
    !!Description:
    !**** *CUCTRACER* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
    !!                  IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES

    !!         P.BECHTOLD        E.C.M.W.F.              11/02/2004

    !**   INTERFACE.
    !!    ----------

    !!         *CUTRACER* IS CALLED FROM *CUMASTR*
    !
    !!    PARAMETER     DESCRIPTION                                   UNITS
    !!    ---------     -----------                                   -----
    !!    INPUT PARAMETERS (INTEGER):

    !!   *KIDIA*        START POINT
    !!   *KFDIA*        END POINT
    !!   *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KTDIA*        START OF THE VERTICAL LOOP
    !!   *KLEV*         NUMBER OF LEVELS
    !!   *KTRAC*        NUMBER OF CHEMICAL TRACERS
    !!   *KCTOP*        CLOUD TOP  LEVEL
     !!   *KDTOP*        DOWNDRAFT TOP LEVEL

    !!   INPUT PARAMETERS (LOGICAL):

    !!   *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
    !!   *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST

    !!   INPUT PARAMETERS (REAL):

    !!   *PTSPHY*       PHYSICS TIME-STEP                              S
    !!   *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
    !!   *zdph*         pressure thickness on half levels              Pa
    !!   *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
    !!   *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
    !!   *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
    !!   *PUDRATE       UPDRAFT DETRAINMENT                           KG/(M2*S)
    !!   *PDDRATE       DOWNDRAFT DETRAINMENT                         KG/(M2*S)

    !!   UPDATED PARAMETERS (REAL):

    !!   *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S

    !!         METHOD
    !!         -------
    !!    EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
    !!    DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT

    !!    FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
    !!    FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
    !!                           TO CORRECT TENDENCIES BELOW CLOUD BASE


    !------------------------------------------------------------------------------------
    !!    COMMENTS FOR OFFLINE USERS IN CHEMICAL TRANSPORT MODELS
    !!    (i.e. reading mass fluxes and detrainment rates from ECMWF archive:
    !!     ------------------------------------------------------------------
    !!    KCTOP IS FIRST LEVEL FROM TOP WHERE PMFU>0
    !!    KDTOP IS FIRST LEVEL FROM TOP WHERE PMFD<0
    !!    ATTENTION: ON ARCHIVE DETRAINMENT RATES HAVE UNITS KG/(M3*S), SO FOR USE
    !!               IN CURRENT ROUTINE YOU HAVE TO MULTIPLY ARCHIVED VALUES BY DZ!!
    !!    LDCUM  IS TRUE IF CONVECTION EXISTS, i.e. IF PMFU>0 IN COLUMN OR IF
    !!                      KCTOP>0 AND KCTOP<KLEV
    !!    LDDRAF IS TRUE IF DOWNDRAUGHTS EXIST IF PMFD<0 IN COLUMN OR IF
    !!                      KDTOP>0 AND KDTOP<KLEV
    !!    IF MASSFLUX SATISFIES CFL CRITERIUM M<=DP/Dt IT IS SUFFICIENT TO
    !!    ONLY CONSIDER EXPLICIT SOLUTION (RMFSOLCT=0), IN THIS CASE
    !!    BUT IMPLICIT SOLUTION (RMFSOLCT=1) PART 7.0 IS PREFERED
    !------------------------------------------------------------------------------------

    !!         EXTERNALS
    !!         ---------
    !!         CUBIDIAG

    !!         MODIFICATIONS
    !!         -------------
    !!       M.Hamrud      01-Oct-2003 CY28 Cleaning

    !!Language: Fortran 90.
    !!Software Standards: "European Standards for Writing and
    !!Documenting Exchangeable Fortran 90 Code".
    !=======================================================================

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: ktrac
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim)               :: ktdia
    INTEGER(KIND=jpim),INTENT(in)    :: kctop(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kdtop(klon)
    LOGICAL ,INTENT(in)    :: ldcum(klon)
    LOGICAL ,INTENT(in)    :: lddraf(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: zdph(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pudrate(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pddrate(klon,klev)
!K.L. test with pcen%ptr and ptenc%ptr
!    REAL(KIND=jprb)   ,INTENT(in)    :: pcen(klon,klev,ktrac)
    TYPE(t_ptr_tracer)   ,INTENT(in), OPTIONAL :: pcen(ktrac)
!    REAL(KIND=jprb)   ,INTENT(inout) :: ptenc(klon,klev,ktrac)
    TYPE(t_ptr_tracer)   ,INTENT(inout),OPTIONAL :: ptenc(ktrac)


    ! Set MODULE PARAMETERS for offline
    !REAL(KIND=JPRB) :: RG=9.80665_JPRB, RMFSOLCT=1.0_JPRB, RMFCMIN=1.E-10_JPRB
    !----------------------------------------------------------------------
    INTEGER(KIND=jpim) :: ik, jk, jl, jn

    REAL(KIND=jprb) :: zzp, zmfa, zimp, zerate, zposi, ztsphy

    !     ALLOCATABLE ARAYS
    REAL(KIND=jprb), DIMENSION(klon,klev) :: &
         zcen, & !< Half-level environmental values
         zcu,  & !< Updraft values
         zcd,  & !< Downdraft values
         ztenc, &!< Tendency
         zmfc, & !< Fluxes
         zdp,  & !< Pressure difference
         zb,   &
         zr1
    LOGICAL :: llcumask(klon,klev)
    REAL(KIND=jprb) :: zhook_handle
    REAL(kind=jprb), POINTER :: tenc(:,:), cen(:,:)

    ! Set MODULE PARAMETERS for offline

    !#include "cubidiag.intfb.h"
    !----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('CUCTRACER',0,zhook_handle)
    zimp=1.0_JPRB-rmfsolct
    ztsphy=1.0_JPRB/ptsphy

    !!Initialize Cumulus mask + some setups

!PREVENT_INCONSISTENT_IFORT_FMA
    DO jk=ktdia+1,klev
      DO jl=kidia,kfdia
        llcumask(jl,jk) = ldcum(jl) .AND. jk >= kctop(jl)-1
        IF(ldcum(jl)) THEN
           zdp(jl,jk)=rg/(paph(jl,jk+1)-paph(jl,jk))
        ENDIF
      ENDDO
    ENDDO
    !----------------------------------------------------------------------

    DO jn=1,ktrac

      !*    1.0          DEFINE TRACERS AT HALF LEVELS
      !!                 -----------------------------
      tenc => ptenc(jn)%ptr
      cen => pcen(jn)%ptr
      DO jk=ktdia+1,klev
        ik=jk-1
        DO jl=kidia,kfdia
         zcen(jl,jk)= cen(jl,jk)
         zcd(jl,jk) = cen(jl,ik)
         zcu(jl,jk) = cen(jl,ik)
         zmfc(jl,jk)=0.0_JPRB
         ztenc(jl,jk)=0.0_JPRB
        ENDDO
      ENDDO

      DO jl=kidia,kfdia
        zcu(jl,klev) = cen(jl,klev)
      ENDDO
      !*    2.0          COMPUTE UPDRAFT VALUES
      !!                 ----------------------

      DO jk=klev-1,ktdia+2,-1
        ik=jk+1
        DO jl=kidia,kfdia
          IF ( llcumask(jl,jk) ) THEN
            zerate=pmfu(jl,jk)-pmfu(jl,ik)+pudrate(jl,jk)
            zmfa=1.0_JPRB/MAX(rmfcmin,pmfu(jl,jk))
            IF (jk >=kctop(jl) )  THEN
              zcu(jl,jk)=( pmfu(jl,ik)*zcu(jl,ik)+zerate*cen(jl,jk) &
                & -pudrate(jl,jk)*zcu(jl,ik) )*zmfa
              !!if you have a source term dc/dt=dcdt write
              !!            ZCU(JL,JK)=( PMFU(JL,IK)*ZCU(JL,IK)+ZERATE*PCEN(JL,JK) &
              !!                          -PUDRATE(JL,JK)*ZCU(JL,IK) )*ZMFA
              !!                          +dcdt(jl,ik)*ptsphy
            ENDIF
          ENDIF
        ENDDO
      ENDDO


      !*    3.0          COMPUTE DOWNDRAFT VALUES
      !!                 ------------------------

      DO jk=ktdia+2,klev
        ik=jk-1
        DO jl=kidia,kfdia
          IF ( lddraf(jl).AND.jk==kdtop(jl) ) THEN
            !Nota: in order to avoid final negative Tracer values at LFS the allowed value of ZCD
            !!     depends on the jump in mass flux at the LFS
            !ZCD(JL,JK)=0.5_JPRB*ZCU(JL,JK)+0.5_JPRB*PCEN(JL,IK)
            zcd(jl,jk)=0.1_JPRB*zcu(jl,jk)+0.9_JPRB*cen(jl,ik)
          ELSEIF ( lddraf(jl).AND.jk>kdtop(jl) ) THEN
            zerate=-pmfd(jl,jk)+pmfd(jl,ik)+pddrate(jl,jk)
            zmfa=1._jprb/MIN(-rmfcmin,pmfd(jl,jk))
            zcd(jl,jk)=( pmfd(jl,ik)*zcd(jl,ik)-zerate*cen(jl,ik) &
              & +pddrate(jl,jk)*zcd(jl,ik) )*zmfa
            !!if you have a source term dc/dt=dcdt write
            !!            ZCD(JL,JK)=( PMFD(JL,IK)*ZCD(JL,IK)-ZERATE*PCEN(JL,IK) &
            !!                          &+PDDRATE(JL,JK)*ZCD(JL,IK) &
            !!                          &+dcdt(jl,ik)*ptsphy
          ENDIF
        ENDDO
      ENDDO

      !!In order to avoid negative Tracer at KLEV adjust ZCD
      jk=klev
      ik=jk-1
      DO jl=kidia,kfdia
        IF (lddraf(jl)) THEN
          !KF
           ZPOSI=-ZDP(JL,JK)*(PMFU(JL,JK)*ZCU(JL,JK)+PMFD(JL,JK)*ZCD(JL,JK)&
          !zposi=-rg/zdph(jl,jk)*(pmfu(jl,jk)*zcu(jl,jk)+pmfd(jl,jk)*zcd(jl,jk)&
            & -(pmfu(jl,jk)+pmfd(jl,jk))*cen(jl,ik) )
          IF( cen(jl,jk)+zposi*ptsphy<0.0_JPRB ) THEN
            zmfa=1._jprb/MIN(-rmfcmin,pmfd(jl,jk))
            zcd(jl,jk)=( (pmfu(jl,jk)+pmfd(jl,jk))*cen(jl,ik)-pmfu(jl,jk)*zcu(jl,jk)&
            !  & +pcen(jl,jk)/(ptsphy*rg/zdph(jl,jk)) )*zmfa
              &+cen(jl,jk)/(PTSPHY*ZDP(JL,JK)) )*ZMFA
          ENDIF
        ENDIF
      ENDDO


      !*    4.0          COMPUTE FLUXES
      !!                 --------------

      DO jk=ktdia+1,klev
        ik=jk-1
        DO jl=kidia,kfdia
          IF(llcumask(jl,jk)) THEN
            zmfa=pmfu(jl,jk)+pmfd(jl,jk)
            zmfc(jl,jk)=pmfu(jl,jk)*zcu(jl,jk)+pmfd(jl,jk)*zcd(jl,jk)&
              & -zimp*zmfa*zcen(jl,ik)
          ENDIF
        ENDDO
      ENDDO

      !*    5.0          COMPUTE TENDENCIES = RHS
      !!                 ------------------------

      DO jk=ktdia+1,klev-1
        ik=jk+1
        DO jl=kidia,kfdia
          IF(llcumask(jl,jk)) THEN
                      ZTENC(JL,JK)=ZDP(JL,JK)*(ZMFC(JL,IK)-ZMFC(JL,JK))
                      !ztenc(jl,jk)=rg/zdph(jl,jk)*(zmfc(jl,ik)-zmfc(jl,jk))
          ENDIF
        ENDDO
      ENDDO

      jk=klev
      DO jl=kidia,kfdia
        IF(ldcum(jl)) THEN
           ZTENC(JL,JK)=-ZDP(JL,JK)*ZMFC(JL,JK)
          !ztenc(jl,jk)=-rg/zdph(jl,jk)*zmfc(jl,jk)
        ENDIF
      ENDDO


      IF ( rmfsolct==0.0_JPRB ) THEN


        !*    6.0          UPDATE TENDENCIES
        !!                 -----------------

        DO jk=ktdia+1,klev
          DO jl=kidia,kfdia
            IF(llcumask(jl,jk)) THEN
              tenc(jl,jk) = tenc(jl,jk) + ztenc(jl,jk)
            ENDIF
          ENDDO
        ENDDO

      ELSE

        !---------------------------------------------------------------------------

        !*    7.0          IMPLICIT SOLUTION
        !!                 -----------------

        !!Fill bi-diagonal Matrix vectors A=k-1, B=k;
        !!reuse ZMFC=A and ZB=B;
        !!ZTENC corresponds to the RHS ("constants") of the equation
        !!The solution is in ZR1

        !!Fill vectors A, B and RHS

        DO jk=ktdia+1,klev
          ik=jk+1
          DO jl=kidia,kfdia
            !!LLCUMBAS(JL,JK)=LLCUMASK(JL,JK).AND.JK<=KCBOT(JL)
            IF(llcumask(jl,jk)) THEN
              !>KF
              ZZP=RMFSOLCT*ZDP(JL,JK)*PTSPHY
              !zzp=rmfsolct*rg/zdph(jl,jk)*ptsphy
              zmfc(jl,jk)=-zzp*(pmfu(jl,jk)+pmfd(jl,jk))
              ztenc(jl,jk) = ztenc(jl,jk)*ptsphy+cen(jl,jk)
              !  for implicit solution including tendency source term
              !  ZTENC(JL,JK) = (ZTENC(JL,JK)+PTENC(JL,JK))*PTSPHY+PCEN(JL,JK)
              IF(jk<klev) THEN
                zb(jl,jk)=1.0_JPRB+zzp*(pmfu(jl,ik)+pmfd(jl,ik))
              ELSE
                zb(jl,jk)=1.0_JPRB
              ENDIF
            ELSE
              zb(jl,jk)=1.0_jprb
            ENDIF
          ENDDO
        ENDDO

        CALL cubidiag&
          & ( kidia, kfdia, klon, ktdia, klev,&
          & kctop, llcumask,&
          & zmfc(:,:),  zb,   ztenc(:,:),   zr1 )

        !!Compute tendencies

        DO jk=ktdia+1,klev
          DO jl=kidia,kfdia
            IF(llcumask(jl,jk)) THEN
              tenc(jl,jk) = tenc(jl,jk) + (zr1(jl,jk) - cen(jl,jk)) * ztsphy
              !  for implicit solution including tendency source term
              !    PTENC(JL,JK)=(ZR1(JL,JK)-PCEN(JL,JK))*ZTSPHY
            ENDIF
          ENDDO
        ENDDO

      ENDIF
    END DO
    !---------------------------------------------------------------------------


    IF (lhook) CALL dr_hook('CUCTRACER',1,zhook_handle)
  END SUBROUTINE cuctracer

  !=======================================================================

  SUBROUTINE cubidiag &
    & ( kidia, kfdia, klon, ktdia, klev,&
    & kctop, ld_lcumask,&
    & pa,    pb,   pr,   pu )
    !
    !!Description:
    !!         P. Bechtold         E.C.M.W.F.     07/03

    !!         PURPOSE.
    !!         --------
    !!         SOLVES BIDIAGONAL SYSTEM
    !!         FOR IMPLICIT SOLUTION OF ADVECTION EQUATION
    !
    !!Method:
    !!
    !!         INTERFACE
    !!         ---------

    !!         THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
    !!         IT RETURNS UPDATED VALUE OF QUANTITY

    !!         METHOD.
    !!         --------
    !!         NUMERICAL RECIPES (Cambridge Press)
    !!         DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
    !!         (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
    !!         M  x  U  = R
    !!         ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C

    !!         (  B(kctop-1)    C(kctop-1)    0          0        )
    !!         (  A(kctop)      B(kctop)      C(kctop)   0        )
    !!         (  0             A(jk)         B(jk)      C(jk)    )
    !!         (  0             0             A(klev)    B(klev)  )

    !!    PARAMETER     DESCRIPTION                                   UNITS
    !!    ---------     -----------                                   -----

    !!   INPUT PARAMETERS (INTEGER):

    !!   *KIDIA*        START POINT
    !!   *KFDIA*        END POINT
    !!   *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KLEV*         NUMBER OF LEVELS
    !!   *KCTOP*        CLOUD TOP LEVELS

    !!   INPUT PARAMETERS (REAL):

    !!   *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS
    !!   *PR*           RHS VECTOR CONTAINING "CONSTANTS"

    !!   OUTPUT PARAMETERS (REAL):

    !!   *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

    !!         EXTERNALS
    !!         ---------
    !!         NONE

    !=======================================================================
    !
    ! Declarations:
    !
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    IMPLICIT NONE

    !     DUMMY INTEGER
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kctop(klon)
    LOGICAL ,INTENT(in)    :: ld_lcumask(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pa(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pb(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pr(klon,klev)
    REAL(KIND=jprb)   ,INTENT(out)   :: pu(klon,klev)
    !     DUMMY REALS
    !     DUMMY LOGICALS
    !     LOCALS
    INTEGER(KIND=jpim) :: jk, jl
    REAL(KIND=jprb) :: zbet
    REAL(KIND=jprb) :: zhook_handle

    !----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('CUBIDIAG',0,zhook_handle)

    pu(:,:)=0.0_JPRB

    ! Forward Substitution
    DO jl = kidia,kfdia
      jk = kctop(jl)-1
      IF (jk >= ktdia+1 .AND. jk <= klev) THEN
        IF ( ld_lcumask(jl,jk) ) THEN
          zbet      =1.0_JPRB/(pb(jl,jk)+1.e-35_JPRB)
          pu(jl,jk) = pr(jl,jk) * zbet
        ENDIF
      END IF
    END DO

    DO jk = MAX(ktdia+1, MINVAL(kctop)), klev
      DO jl = kidia,kfdia
        IF ( jk >= kctop(jl) .AND. ld_lcumask(jl,jk) ) THEN
          zbet      = 1.0_JPRB/(pb(jl,jk) + 1.e-35_JPRB)
          pu(jl,jk) =(pr(jl,jk)-pa(jl,jk)*pu(jl,jk-1))*zbet
        ENDIF
      ENDDO
    ENDDO

    IF (lhook) CALL dr_hook('CUBIDIAG',1,zhook_handle)
  END SUBROUTINE cubidiag

END MODULE mo_cuflxtends

