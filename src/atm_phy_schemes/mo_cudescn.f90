! $RCSfile$
! $Revision$ $Date$
!

!>
!! *cudlfsn*  THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!! *cuddrafn* THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!!
!! @author  M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!!
!!
!! @par Revision History
!! first implementation into ICON by Kristina Froehlich (2010-05-27)
!! Added pvbuo by calculation of gusts  by Helmut Frank (2013-03-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#ifndef _OPENMP
#include "consistent_fma.inc"
#endif
MODULE mo_cudescn

#ifdef __ICON__
  USE mo_kind   ,ONLY: jprb=>wp     , &
    &                  jpim=>i4
#endif

#ifdef __GME__
!  USE parkind1  ,ONLY : jpim     ,jprb
  USE gme_data_parameters, ONLY:  JPRB =>ireals, JPIM => iintegers
#endif
  
!  USE yomhook   ,ONLY : lhook,   dr_hook
  
  !KF new - use modules instead of include files
  USE mo_adjust      , ONLY: cuadjtq
  USE mo_cufunctions , ONLY: foelhmcu
  USE mo_cuparameters, ONLY: lphylin  ,rlptrc, rg ,rcpd     ,retv,&
    &                        rlvtt    ,rlstt    , &
    &                        entrdd   ,rmfcmin,       &
    &                        rmfdeps  ,rmfdeps_ocean, lmfdd,   &
    &                        lhook,   dr_hook
  
  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cudlfsn, cuddrafn

CONTAINS

  SUBROUTINE cudlfsn &
    & (kidia,    kfdia,    klon,    ktdia,  klev,&
    & kcbot,    kctop,     ldcum,&
    & ptenh,    pqenh,   &
    & pten,     pqsen,    pgeo,&
    & pgeoh,    paph,     ptu,      pqu,  &
    & pmfub,    prfl,&
    & ptd,      pqd,&
    & pmfd,     pmfds,    pmfdq,    pdmfdp,&
    & kdtop,    lddraf, ldland,   ldlake)
    !!
    !! Description:
    !!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
    !!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES

    !!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89

    !!          PURPOSE.----
    !!          --------
    !!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
    !!          FOR MASSFLUX CUMULUS PARAMETERIZATION

    !!          INTERFACE
    !!          ---------
    !!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
    !!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
    !!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
    ! !         CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
    !!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.

    !!          METHOD.

    !!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
    !!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.

    !!
    !! Code Description:
    !!     PARAMETER     DESCRIPTION                                   UNITS
    !!     ---------     -----------                                   -----
    !!     INPUT PARAMETERS (INTEGER):

    !!    *KIDIA*        START POINT
    !!    *KFDIA*        END POINT
    !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!    *KLEV*         NUMBER OF LEVELS
    !!    *KCBOT*        CLOUD BASE LEVEL
    !!    *KCTOP*        CLOUD TOP LEVEL

    !!    INPUT PARAMETERS (LOGICAL):
    !!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

    !!    INPUT PARAMETERS (REAL):

    !!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
    !!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
    !!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
    !!    *PGEO*         GEOPOTENTIAL                                  M2/S2
    !!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
    !!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
    !!    *PTU*          TEMPERATURE IN UPDRAFTS                        K
    !!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
    !!   *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
    !!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)

    !!    UPDATED PARAMETERS (REAL):

    !!    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)

    !!    OUTPUT PARAMETERS (REAL):

    !!    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
    !!    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
    !!    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
    !!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
    !!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
    !!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)

    !!    OUTPUT PARAMETERS (INTEGER):

    !!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS

    !!    OUTPUT PARAMETERS (LOGICAL):

    !!    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST

    !!          EXTERNALS
    !!          ---------
    !!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS

    !!          MODIFICATIONS
    !!          -------------
    !!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
    !!             99-06-04 : Optimisation        D.SALMOND
    !!        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !!        P. Lopez      20-Jun-2007 CY32R2 Bug correction in latent heat
    !!                                         when LPHYLIN=T.
    !=======================================================================


    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT
    !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
    ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
    ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
    ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
    !USE YOECUMF  , ONLY : RMFDEPS  ,LMFDD
    !USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    INTEGER(KIND=jpim)               :: kcbot(klon) ! Argument NOT used
    INTEGER(KIND=jpim)               :: kctop(klon) ! Argument NOT used
    LOGICAL :: ldcum(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqsen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqu(klon,klev)
!    REAL(KIND=jprb)                  :: plu(klon,klev) ! Argument NOT used
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfub(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: prfl(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pqd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfds(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfdq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pdmfdp(klon,klev)
    INTEGER(KIND=jpim),INTENT(out)   :: kdtop(klon)
    LOGICAL ,INTENT(in)              :: ldland(klon),   ldlake(klon)
    LOGICAL ,INTENT(out)   :: lddraf(klon)
    INTEGER(KIND=jpim) ::            ikhsmin(klon)
    REAL(KIND=jprb) ::     ztenwb(klon,klev),      zqenwb(klon,klev),&
      & zcond(klon),            zph(klon),&
      & zhsmin(klon)
    LOGICAL ::  llo2(klon)

    INTEGER(KIND=jpim) :: icall, ik, ike, is, jk, jl

    REAL(KIND=jprb) :: zbuo, zhsk, zmftop, zoealfa,&
      & zoelhm, zqtest, ztarg, zttest
    REAL(KIND=jprb) :: zhook_handle

    !#include "cuadjtq.intfb.h"
    !#include "fcttre.h"
    !----------------------------------------------------------------------

    !!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
    !!                  ---------------------------------

    IF (lhook) CALL dr_hook('CUDLFSN',0,zhook_handle)
!PREVENT_INCONSISTENT_IFORT_FMA
    DO jl=kidia,kfdia
      lddraf(jl)=.FALSE.
      kdtop(jl)=klev+1
      ikhsmin(jl)=klev+1
      zhsmin(jl)=1.e8_jprb
    ENDDO
    ztenwb(:,:) = 0._jprb
    zqenwb(:,:) = 0._jprb
    zcond (:)   = 0._jprb
    zph   (:)   = 0._jprb


    !orig IF(.NOT.LMFDD) GO TO 300
    IF(lmfdd) THEN

      !----------------------------------------------------------------------

      !!    2.           DETERMINE LEVEL OF FREE SINKING:
      !!                  DOWNDRAFTS SHALL START AT MODEL LEVEL OF MINIMUM
      !!                  OF SATURATION MOIST STATIC ENERGY OR BELOW
      !!                  RESPECTIVELY

      !!                  FOR EVERY POINT AND PROCEED AS FOLLOWS:

      !!                    (1) DETERMINE LEVEL OF MINIMUM OF HS
      !!                    (2) DETERMINE WET BULB ENVIRONMENTAL T AND Q
      !!                    (3) DO MIXING WITH CUMULUS CLOUD AIR
      !!                    (4) CHECK FOR NEGATIVE BUOYANCY
      !!                    (5) IF BUOYANCY>0 REPEAT (2) TO (4) FOR NEXT
      !!                        LEVEL BELOW

      !!                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
      !!                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
      !!                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
      !!                  EVAPORATION OF RAIN AND CLOUD WATER)
      !!                  ----------------------------------------------------

      DO jk=ktdia+2,klev-2

        IF (lphylin) THEN

          DO jl=kidia,kfdia
            ztarg=pten(jl,jk)
            zoealfa=0.545_JPRB*(TANH(0.17_JPRB*(ztarg-rlptrc))+1.0_JPRB)
            zoelhm =zoealfa*rlvtt+(1.0_JPRB-zoealfa)*rlstt
            zhsk=rcpd*pten(jl,jk)+pgeo(jl,jk)+zoelhm*pqsen(jl,jk)
            IF(zhsk < zhsmin(jl)) THEN
              zhsmin(jl)=zhsk
              ikhsmin(jl)=jk
            ENDIF
          ENDDO

        ELSE

          DO jl=kidia,kfdia
            zhsk=rcpd*pten(jl,jk)+pgeo(jl,jk)+foelhmcu(pten(jl,jk))*pqsen(jl,jk)
            IF(zhsk < zhsmin(jl)) THEN
              zhsmin(jl)=zhsk
              ikhsmin(jl)=jk
            ENDIF
          ENDDO

        ENDIF

      ENDDO
      ike=klev-3
      DO jk=ktdia+2,ike

        !!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
        !!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
        !!                  -------------------------------------------

        is=0
        DO jl=kidia,kfdia
          ztenwb(jl,jk)=ptenh(jl,jk)
          zqenwb(jl,jk)=pqenh(jl,jk)
          zph(jl)=paph(jl,jk)
          llo2(jl)=ldcum(jl).AND.prfl(jl) > 0.0_JPRB.AND..NOT.lddraf(jl).AND.&
            & (jk < kcbot(jl).AND.jk > kctop(jl)).AND.&
            & jk >= ikhsmin(jl)
          IF(llo2(jl))THEN
            is=is+1
          ENDIF
        ENDDO
        !orig   IF(IS.EQ.0) GO TO 290
        IF(is == 0) CYCLE

        ik=jk
        icall=2
        CALL cuadjtq &
          & ( kidia,    kfdia,    klon,   klev,&
          & ik,&
          & zph,      ztenwb,   zqenwb,   llo2,     icall)

        !!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
        !!                  AND CHECK FOR NEGATIVE BUOYANCY.
        !!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
        !!                  ----------------------------------------

!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          IF(llo2(jl)) THEN
            zttest=0.5_JPRB*(ptu(jl,jk)+ztenwb(jl,jk))
            zqtest=0.5_JPRB*(pqu(jl,jk)+zqenwb(jl,jk))
            zbuo=zttest*(1.0_JPRB+retv  *zqtest)-&
              & ptenh(jl,jk)*(1.0_JPRB+retv  *pqenh(jl,jk))
            zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)
            zmftop=-MERGE(rmfdeps,rmfdeps_ocean,ldland(jl).OR.ldlake(jl))*pmfub(jl)
            IF(zbuo < 0.0_JPRB.AND.prfl(jl) > 10._jprb*zmftop*zcond(jl)) THEN
              kdtop(jl)=jk
              lddraf(jl)=.TRUE.
              ptd(jl,jk)=zttest
              pqd(jl,jk)=zqtest
              pmfd(jl,jk)=zmftop
              pmfds(jl,jk)=pmfd(jl,jk)*(rcpd*ptd(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1)=-0.5_JPRB*pmfd(jl,jk)*zcond(jl)
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)
            ENDIF
          ENDIF
        ENDDO

        ! 290   continue
      ENDDO

      !300  CONTINUE
    ENDIF

    IF (lhook) CALL dr_hook('CUDLFSN',1,zhook_handle)
  END SUBROUTINE cudlfsn

  !=======================================================================

  SUBROUTINE cuddrafn &
    & ( kidia,    kfdia,    klon,     ktdia,  klev, k950, &
    & lddraf,&
    & ptenh,    pqenh,                   &
    & pgeo,     pgeoh,    paph,     prfl,&
    & zdph,     zdgeoh,         &
    & ptd,      pqd,      pmfu,&
    & pmfd,     pmfds,    pmfdq,    pdmfdp,&
    & pdmfde,   pmfdde_rate,        pkined,&
    & pvbuo )
    !!
    !! Description:

    !!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT

    !!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89

    !!          PURPOSE.
    !!          --------
    !!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
    !!          (I.E. T,Q,U AND V AND FLUXES)

    !!          INTERFACE
    !!          ---------

    !!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
    !!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
    !!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
    !!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS

    !!          METHOD.
    !!          --------
    !!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
    !!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
    !!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.

    !!     PARAMETER     DESCRIPTION                                   UNITS
    !!     ---------     -----------                                   -----
    !!     INPUT PARAMETERS (INTEGER):

    !!    *KIDIA*        START POINT
    !!    *KFDIA*        END POINT
    !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!    *KLEV*         NUMBER OF LEVELS

    !!    INPUT PARAMETERS (LOGICAL):

    !!    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST

    !!    INPUT PARAMETERS (REAL):

    !!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
    !!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
    !!    *PGEO*         GEOPOTENTIAL                                  M2/S2
    !!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
    !!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
    !!    *zdgeoh*       geopot thickness on half levels               m2/s2
    !!    *zdph*         pressure thickness on half levels              Pa
    !!    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)

    !!    UPDATED PARAMETERS (REAL):

    !!   *PRFL*         PRECIPITATION RATE                           KG/(M2*S)

    !!    OUTPUT PARAMETERS (REAL):

    !!    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
    !!    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
    !!    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
    !!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
    !!   *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
    !!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)
    !!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                   KG/(M2*S)
    !!    *PKINED*       DOWNDRAFT KINETIC ENERGY                     M2/S2

    !!          EXTERNALS
    !!          ---------
    !!         *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
    !!          SATURATED DESCENT

    !!          REFERENCE
    !!          ---------
    !!          (TIEDTKE,1989)

    !!          MODIFICATIONS
    !!          -------------
    !!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
    !!             03-08-28 : Clean-up detrainment rates   P. BECHTOLD
    !!        M.Hamrud      01-Oct-2003 CY28 Cleaning

    !=======================================================================
    !
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RG       ,RCPD     ,RETV
    !USE YOECUMF  , ONLY : ENTRDD   ,RMFCMIN  ,NJKT3

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    INTEGER(KIND=jpim),INTENT(in)    :: k950(klon)
    LOGICAL           ,INTENT(in)    :: lddraf(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: zdgeoh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: zdph(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: prfl(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pqd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfds(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfdq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pdmfdp(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pdmfde(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pmfdde_rate(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pkined(klon,klev)
    REAL(KIND=jprb)   ,INTENT(OUT)     :: pvbuo(klon)  ! buoyancy for convective gusts 
    REAL(KIND=jprb) ::     zdmfen(klon),           zdmfde(klon),&
      & zcond(klon),            zoentr(klon),&
      & zbuoy(klon)
    REAL(KIND=jprb)    :: zph(klon)
    LOGICAL            :: llo2(klon)
    INTEGER(KIND=jpim) :: icall, ik, is, jk, jl

    REAL(KIND=jprb) :: zbuo, zbuoyz, zbuoyv, zdmfdp, zdz, zentr, zmfdqk,&
      & zmfdsk, zqdde, zqeen, zrain, &
      & zsdde, zseen, zzentr, zrg, zfacbuo, z_cwdrag , zdkbuo, zdken
    REAL(KIND=jprb) :: zhook_handle

!   calculation of convective gusts
    REAL(KIND=jprb) :: zqprec                    ! rain mixing ratio in downdraft
    REAL(KIND=jprb) :: conv_gust_rain = 0._jprb  ! factor for zqprec

    !#include "cuadjtq.intfb.h"

    IF (lhook) CALL dr_hook('CUDDRAFN',0,zhook_handle)
    zrg=1.0_JPRB/rg
    zfacbuo=0.5_JPRB/(1.0_JPRB+0.5_JPRB)
    z_cwdrag=(3._jprb/8._jprb)*0.506_JPRB/0.2_JPRB/rg
    !!----------------------------------------------------------------------

    !!     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
    !!                     (A) CALCULATING ENTRAINMENT/DETRAINMENT RATES,
    !!                         INCLUDING ORGANIZED ENTRAINMENT DEPENDENT ON
    !!                         NEGATIVE BUOYANCY AND ASSUMING
    !!                        LINEAR DECREASE OF MASSFLUX IN PBL
    !!                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
    !!                         AND MOISTENING IS CALCULATED IN *CUADJTQ*
    !!                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
    !!                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
    !!                    -------------------------------------------------

    DO jl=kidia,kfdia
      zoentr     (jl)  =0.0_JPRB
      zbuoy      (jl)  =0.0_JPRB
      zdmfen     (jl)  =0.0_JPRB
      zdmfde     (jl)  =0.0_JPRB
      zcond      (jl)  =0.0_JPRB
      pmfdde_rate(jl,:)=0.0_JPRB
      pkined     (jl,:)=0.0_JPRB
    ENDDO
    pvbuo        (:)=0.0_JPRB

    DO jk=ktdia+2,klev
      is=0
      DO jl=kidia,kfdia
        zph(jl)=paph(jl,jk)
        llo2(jl)=lddraf(jl).AND.pmfd(jl,jk-1) < 0.0_JPRB
        IF(llo2(jl)) THEN
          is=is+1
        ENDIF
      ENDDO
      IF(is == 0) CYCLE

      DO jl=kidia,kfdia
        IF(llo2(jl)) THEN
          !>KF
          ZENTR=ENTRDD*PMFD(JL,JK-1)*(PGEOH(JL,JK-1)-PGEOH(JL,JK))*ZRG
          ! zentr=entrdd*pmfd(jl,jk-1)*zdgeoh(jl,jk)*zrg
          !<KF
          zdmfen(jl)=zentr
          zdmfde(jl)=zentr
        ENDIF
      ENDDO

      DO jl=kidia,kfdia
        IF (jk > k950(jl)) THEN
          IF(llo2(jl)) THEN
            zdmfen(jl)=0.0_JPRB
            zdmfde(jl)=pmfd(jl,k950(jl))*&
            !  & zdph(jl,jk-1)    /&
                     & (PAPH(JL,JK)-PAPH(JL,JK-1))/&
              & (paph(jl,klev+1)-paph(jl,k950(jl)))
          ENDIF
        ENDIF
      ENDDO


      DO jl=kidia,kfdia
        IF (jk <= k950(jl)) THEN
          IF(llo2(jl)) THEN
            !>KF
            ZDZ=-(PGEOH(JL,JK-1)-PGEOH(JL,JK))*ZRG
            !zdz=-zdgeoh(jl,jk-1)*zrg
            !<KF
            zzentr=zoentr(jl)*zdz*pmfd(jl,jk-1)
            zdmfen(jl)=zdmfen(jl)+zzentr
            zdmfen(jl)=MAX(zdmfen(jl),0.3_JPRB*pmfd(jl,jk-1))
            zdmfen(jl)=MAX(zdmfen(jl),-0.75_JPRB*pmfu(jl,jk)-&
              & (pmfd(jl,jk-1)-zdmfde(jl)))
            zdmfen(jl)=MIN(zdmfen(jl),0.0_JPRB)
          ENDIF

          pdmfde(jl,jk)=zdmfen(jl)-zdmfde(jl)
        ENDIF
      ENDDO

      DO jl=kidia,kfdia
        IF(llo2(jl)) THEN
          pmfd(jl,jk)=pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
          zseen=(rcpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
          zqeen=pqenh(jl,jk-1)*zdmfen(jl)
          zsdde=(rcpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
          zqdde=pqd(jl,jk-1)*zdmfde(jl)
          zmfdsk=pmfds(jl,jk-1)+zseen-zsdde
          zmfdqk=pmfdq(jl,jk-1)+zqeen-zqdde
          pqd(jl,jk)=zmfdqk*(1.0_JPRB/MIN(-rmfcmin,pmfd(jl,jk)))
          ptd(jl,jk)=(zmfdsk*(1.0_JPRB/MIN(-rmfcmin,pmfd(jl,jk)))-pgeoh(jl,jk))/rcpd
          ptd(jl,jk)=MIN(400._jprb,ptd(jl,jk))
          ptd(jl,jk)=MAX(100._jprb,ptd(jl,jk))
          zcond(jl)=pqd(jl,jk)
        ENDIF
      ENDDO

      ik=jk
      icall=2
      CALL cuadjtq &
        & ( kidia,    kfdia,    klon,   klev,&
        & ik,&
        & zph,      ptd,      pqd,      llo2,     icall )

      DO jl=kidia,kfdia
        IF(llo2(jl)) THEN
          zcond(jl)=zcond(jl)-pqd(jl,jk)
          zbuo=ptd(jl,jk)*(1.0_JPRB+retv  *pqd(jl,jk))-&
            & ptenh(jl,jk)*(1.0_JPRB+retv  *pqenh(jl,jk))
          IF(prfl(jl) > 0.0_JPRB.AND.pmfu(jl,jk) > 0.0_JPRB) THEN
            zrain=prfl(jl)/pmfu(jl,jk)
            zbuo=zbuo-ptd(jl,jk)*zrain
          ENDIF
          IF(zbuo >= 0.0_JPRB.OR.prfl(jl) <= (pmfd(jl,jk)*zcond(jl))) THEN
            pmfd(jl,jk)=0.0_JPRB
            zbuo=0.0_JPRB
          ENDIF
          pmfds(jl,jk)=(rcpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
          pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)
          zdmfdp=-pmfd(jl,jk)*zcond(jl)
          pdmfdp(jl,jk-1)=zdmfdp
          prfl(jl)=prfl(jl)+zdmfdp

          !! COMPUTE ORGANIZED ENTRAINMENT FOR USE AT NEXT LEVEL

          zbuoyz=zbuo/ptenh(jl,jk)
          zbuoyv=zbuoyz
          zbuoyz=MIN(zbuoyz,0.0_JPRB)
          !>KF
          zdz=-(pgeo(jl,jk-1)-pgeo(jl,jk))
          !       ZDZ= - zdgeo(jl,jk) !full levels!
          !<KF
          zbuoy(jl)=zbuoy(jl)+zbuoyz*zdz
          zoentr(jl)=rg*zbuoyz*0.5_JPRB/(1.0_JPRB+zbuoy(jl))

          !! STORE DOWNDRAUGHT DETRAINMENT RATES

          pmfdde_rate(jl,jk)=-zdmfde(jl)

          !! COMPUTE KINETIC ENERGY

          zdkbuo=zdz*zbuoyv*zfacbuo
          IF(zdmfen(jl) < 0.0_JPRB)THEN
            zdken=MIN(1.0_JPRB,(1.0_JPRB + rg*z_cwdrag)*&
              & zdmfen(jl)/MIN(-rmfcmin,pmfd(jl,jk-1)))
          ELSE
            zdken=MIN(1.0_JPRB,(1.0_JPRB + rg*z_cwdrag)*&
              & zdmfde(jl)/MIN(-rmfcmin,pmfd(jl,jk-1)))
          ENDIF
          pkined(jl,jk)=MAX(0.0_JPRB,(pkined(jl,jk-1) &
              &  *(1.0_JPRB-zdken)+zdkbuo)/(1.0_JPRB+zdken))

          ! custs generated by buoyancy forces
          zqprec = conv_gust_rain*prfl(jl)/MAX( rmfcmin, -pmfd(jl,jk-1) )
          pvbuo(jl) = pvbuo(jl) + 2._JPRB*( zbuoyv + zqprec )*zdz

        ENDIF
      ENDDO

    ENDDO

    IF (lhook) CALL dr_hook('CUDDRAFN',1,zhook_handle)
  END SUBROUTINE cuddrafn

END MODULE mo_cudescn

