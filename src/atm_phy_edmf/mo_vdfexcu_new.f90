!>
!! Exchange coefficients for EDMF
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2017-8-17)
!!   (IFS cycle CY41r2)
!!
!!-----------------------------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!-----------------------------------------------------------------------------

MODULE mo_vdfexcu_new
 
  PUBLIC :: vdfexcu_new

CONTAINS

!! #ifdef RS6K
!! @PROCESS NOSTRICT
!! #endif
SUBROUTINE VDFEXCU_NEW ( &
                  &KIDIA  , KFDIA  , KLON   , KLEV   , PTMST  , PZ0MM  , &
                  &PHRLW  , PHRSW  , PUM1   , PVM1   , PTM1   , PQM1   , &
!xmk              &PAPHM1 , PGEOM1 , PGEOH  , PGELAT , PCPTGZ , &
                  &PAPHM1 , PGEOM1 , PGEOH  ,          PCPTGZ , &
                  &PKMFL  , PKHFL  , PKQFL  , PCFM   , PCFH   , &
                  &PZINV  , PKH    , PKM    , PZCLDBASE , KPBLTYPE )
!     ------------------------------------------------------------------

!**   *VDFEXCU* - DETERMINES THE EXCHANGE COEFFICIENTS BETWEEN THE
!                 UPPER MODEL LEVELS WITH STABILITY AS A FUNCTION OF
!                 OBUKHOV-L

!     PURPOSE
!     -------

!     DETERMINE EXCHANGE COEFFICIENTS BETWEEN THE UPPER MODEL LEVELS

!     INTERFACE
!     ---------

!     *VDFEXCU* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMUDITY AT T-1
!     *PAPHM1*       PRESSURE AT HALF LEVELS AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PCPTGZ*       DRY STATIC ENERGY
!     *PKMFL*        KINEMATIC MOMENTUM FLUX
!     *PKHFL*        KINEMATIC HEAT FLUX
!     *PKQFL*        KINEMATIC MOISTURE FLUX
!     *PZINV*        INVERSION HEIGHT                  [M]
!     *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!     *PKM*          TURB. DIFF. COEFF. FOR MOM. ABOVE SURF. LAY.  (M2/S)

!     OUTPUT PARAMETERS (REAL):

!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT     (C-STAR IN DOC.)
!                    (ONLY PCFM(*,1:KLEV-1) AND
!                          PCFH(*,1:KLEV-1) ARE COMPUTED)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     AUTHOR.
!     -------
!      A.C.M. BELJAARS  26/03/90.  Original

!     MODIFICATIONS.
!     --------------
!      J.HAGUE          13/01/2003 MASS Vector Functions
!      M. Ko"hler        3/12/2004 Moist Advection-Diffusion incl.
!                                  K,cloud and cloud top entrainment
!      P. Lopez         02/06/2005 Removed option for linearized
!                                  physics (now called separately)
!      M. Ko"hler        1/11/2007 reduced K diffusion above surface or mixed layer
!      M. Ko"hler        2/15/2008 added shear for K and Ri calculation
!      N. Semane+P.Becht 8/06/2012 scaling for small planet
!      A.Beljaars+I.Sandu 1/11/2012 smooth reduction in diffusion above tropopause 
!      +T.Stockdale+P.Bechtold (optional,by default : not active)
!      N.Semane+P.Bechtold 04-10-2012 Add RPLRG/RPLDARE factors for small planet
!      I. Sandu+A.Beljaars 15/3/2013 changed treatment of diffusion in stable conditions,i.e.
!                              LTG functions allover, asymptotic mixing length proportional to PBL
!                              height within stable boundary layers and equal to 30m in free-shear layers,
!                              removed non-resolved shear term
!      F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------------------------

! USE PARKIND1 , ONLY : JPIM, JPRB
! USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
! USE YOMCST   , ONLY : RG, RD, RCPD, RETV
! USE YOEVDF   , ONLY : RKAP, REPDU2, RLAM
! USE YOEVDFS  , ONLY : RCHBA, RCHBB, RCHBD, RCHB23A, RCHBBCD, RCHBCD, &
!  & RCHETA, RCHETB, RCDHALF, RCDHPI2
! USE YOEGWD   , ONLY : YREGWD
! USE YOEPHY   , ONLY : YREPHY
! USE YOMJFH   , ONLY : N_VMASS
! USE YOMDYNCORE,ONLY : RPLRG, RPLDARE

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,vrec     ,&
                & RG       ,RD       ,RCPD     ,RETV     ,&           !yomcst
                & RKAP     ,RVDIFTS  ,REPDU2   ,RLAM     ,&           !yoevdf
                & phihu    ,phimu    ,phims    ,phihs                 !fcvdfs.h
USE mo_edmf_param   ,ONLY : &
                & N_VMASS  ,&                                         !yomjfh
                & REPUST                                              !yos_exc

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
!dmk REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZINV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZCLDBASE(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPBLTYPE(KLON) 

REAL(KIND=JPRB) ::    ZRI(KLON),ZMGEOM(KLON),ZUST(KLON),&
                    & ZDTV(KLON),ZKHVFL(KLON),ZL(KLON),ZPHIM(KLON),&
                    & ZPHIH(KLON)  
REAL(KIND=JPRB) ::    ZDU2(KLON+N_VMASS)

INTEGER(KIND=JPIM) :: JK, JL, JLEN
INTEGER(KIND=JPIM) :: INDX(KFDIA-KIDIA+1), JJ, JJJ
REAL(KIND=JPRB) ::    ZENTRSFC, ZENTRRAD, ZKLEN2, &
                    & ZCB, ZCD, ZCFNC1, ZCONS13, ZRG, &
                    & ZDH, ZDL, ZDRORO, ZEPS, ZETA, &
                    & ZPHIKH, ZPHIKM, ZSCF, &
                    & ZZ, ZWTVENTR, ZKH, ZCFHNEW, &
                    & ZML, ZBASE, ZVSC, ZKCLD, ZDRADFLX(KLON), &
                    & ZREPUST,ZRA, ZLAT
REAL(KIND=JPRB) ::    ZZH, ZIFLTGM, ZIFLTGH, ZIFMOM(KLON), ZIFMOH(KLON), ZDUDZ(KLON)
LOGICAL ::            LLDONE(KLON)
REAL(KIND=JPRB) ::    ZPBLHEIGHT(KLON),ZRIB(KLON),ZSVBOT(KLON),ZRILEV,&
                    & ZRICRI,ZKLENT(KLON,KLEV),ZSV

REAL(KIND=JPRB) ::    ZTMP2(KFDIA-KIDIA+1)
REAL(KIND=JPRB) ::    ZTMP4(KFDIA-KIDIA+1)
REAL(KIND=JPRB) ::    ZTMP5(KFDIA-KIDIA+1)
REAL(KIND=JPRB) ::    ZEPSILON

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

!amk
REAL(KIND=JPRB), PARAMETER :: RPLRG   = 1.0_JPRB
REAL(KIND=JPRB), PARAMETER :: RPLDARE = 1.0_JPRB
!xxx

!dmk #include "surf_inq.h"

! #include "fcvdfs.func.h"


!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFEXCU',0,ZHOOK_HANDLE)

ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

ZENTRSFC  = 0.2_JPRB       ! factor for surface based top entrainment 
ZENTRRAD  = 0.2_JPRB       ! factor for radiative based top entrainment 
ZCD       = 1.0_JPRB
ZCB       = 5.0_JPRB
ZEPS      = 1.E-10_JPRB
ZRICRI   = 0.25_JPRB

!xmk CALL SURF_INQ(YREPHY%YSURF,PREPUST=ZREPUST)
ZREPUST = REPUST

! optimization
ZRG       = 1.0_JPRB/RG
ZCONS13   = 1.0_JPRB/3._JPRB

IF(N_VMASS > 0) THEN
  JLEN=KFDIA-KIDIA+1
ENDIF


!     ------------------------------------------------------------------

!*         2.     PREPARE SCALING COEFFICIENTS FOR MIXED LAYER
!                 --------------------------------------------

DO JL=KIDIA,KFDIA
  ZUST  (JL)=MAX(SQRT(PKMFL(JL)),ZREPUST)
  ZKHVFL(JL)=(PKHFL(JL)+RETV*PTM1(JL,KLEV)*PKQFL(JL))/(RPLRG*RPLDARE)
  IF (ZKHVFL(JL)  <  0.0_JPRB) THEN
    ZL(JL)  = ZUST (JL)**3*PTM1(JL,KLEV)/(RKAP*RG*(ZKHVFL(JL)-ZEPS))
  ENDIF
ENDDO

!          Calculate PBL cloud top radiative flux jump [Km/s] (cooling)
!          for top-driven K and entrainment velocity formulations.

JJJ = 0
DO JL=KIDIA,KFDIA
  IF(KPBLTYPE(JL)==2) THEN
     JJJ=JJJ+1
     INDX(JJJ)=JL
  ENDIF
ENDDO


ZDRADFLX(:) = 0.0_JPRB
!DO JK=KLEV-1,1,-1
 DO JK=1,KLEV-1
! DO JL=KIDIA,KFDIA
!   IF ( KPBLTYPE(JL) == 2 ) THEN
  DO JJ=1,JJJ
    JL=INDX(JJ)
      IF ( PGEOH(JL,JK)*ZRG <= PZINV(JL) .AND. PZINV(JL) < PGEOH(JL,JK-1)*ZRG ) THEN
!       ZDRADFLX(JL) = -  PHRLW(JL,JK+1)                 * (PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZRG
        ZDRADFLX(JL) = - (PHRLW(JL,JK+1)+PHRSW(JL,JK+1)) * (PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZRG
!         ... add solar heating 2nd cld level
        IF ( PZCLDBASE(JL) < PGEOH(JL,JK+1)*ZRG .AND. JK < KLEV-1 ) THEN 
          ZDRADFLX(JL) = ZDRADFLX(JL) - PHRSW(JL,JK+2)   * (PGEOH(JL,JK+1)-PGEOH(JL,JK+2))*ZRG
        ENDIF
        ZDRADFLX(JL) = MAX( ZDRADFLX(JL)/RPLDARE, 0.0_JPRB )    !safety against rad. heating cases
      ENDIF
!   ENDIF
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         3.     VERTICAL LOOP
!                 -------------
!   compute pbl height in stable boundary layers
 DO JL=KIDIA,KFDIA
   LLDONE(JL)=.FALSE.
   ZPBLHEIGHT(JL)=0.0_JPRB
   ZRIB(JL)=0.0_JPRB
   ZSVBOT(JL)=RCPD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV))+PGEOM1(JL,KLEV)+ZEPS
 ENDDO
 DO JK = KLEV-1, 1, -1
    DO JL=KIDIA,KFDIA
      IF (.NOT. LLDONE(JL) .AND. ZKHVFL(JL)  >  0.0_JPRB) THEN
        ZSV=RCPD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK))+PGEOM1(JL,JK)
!   pbl height diag, which considers the winds close to the surf eq to 0
        ZDU2(JL)=MAX(REPDU2, PUM1(JL,JK)**2+PVM1(JL,JK)**2)  
        ZDRORO=(ZSV-ZSVBOT(JL))&
         & /(ZSV-PGEOM1(JL,JK))  
        ZRILEV=(PGEOM1(JL,JK)-PGEOM1(JL,KLEV))*ZDRORO/ZDU2(JL)
!
          IF (ZRILEV  >  ZRICRI) THEN
           ZPBLHEIGHT(JL)=( (ZRILEV-ZRICRI)*PGEOM1(JL,JK+1)&
           & +(ZRICRI-ZRIB(JL))*PGEOM1(JL,JK) )/&
           & ((ZRILEV-ZRIB(JL))*RG)  
           ZRIB(JL)=ZRILEV
           LLDONE(JL)=.TRUE.
          ELSE
           ZRIB(JL)=ZRILEV
          ENDIF
      ENDIF
    ENDDO
 ENDDO

!initialization
DO JK = KLEV, 1, -1
  DO JL=KIDIA,KFDIA
    ZKLENT(JL,JK) =0.0_JPRB
  ENDDO
ENDDO

!***
DO JK = KLEV-1, 1, -1
!***

  DO JL=KIDIA,KFDIA
    PCFM(JL,JK)=0.0_JPRB
    PCFH(JL,JK)=0.0_JPRB
    PKH(JL,JK) =0.0_JPRB
    PKM(JL,JK) =0.0_JPRB
  ENDDO

  IF(N_VMASS <= 0) THEN   ! efficiency of exponentials

!          COMPUTE RI-NUMBER

    DO JL=KIDIA,KFDIA
      ZMGEOM(JL)=PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
      ZDU2(JL)=MAX(REPDU2,(PUM1(JL,JK)-PUM1(JL,JK+1))**2&
                       & +(PVM1(JL,JK)-PVM1(JL,JK+1))**2)
      ZDUDZ(JL)= SQRT(ZDU2(JL)) /ZMGEOM(JL)*RG

      ZDRORO= 2.0_JPRB * (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
       & / ( PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)&
       &   - PGEOM1(JL,JK)-PGEOM1(JL,JK+1))&
       & + RETV*(PQM1(JL,JK)-PQM1(JL,JK+1))
      ZDTV(JL)=( (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
       & + RETV*0.5_JPRB * (PQM1(JL,JK)-PQM1(JL,JK+1))&
       & * (PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)) ) * (1.0_JPRB/RCPD)  
      ZRI(JL)=ZMGEOM(JL)*ZDRORO/ZDU2(JL)
    ENDDO

  ELSE

    DO JL=KIDIA,KFDIA
      ZMGEOM(JL)=PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
      ZDU2(JL)=MAX(REPDU2,(PUM1(JL,JK)-PUM1(JL,JK+1))**2&
                       & +(PVM1(JL,JK)-PVM1(JL,JK+1))**2)
      ZDUDZ(JL)= SQRT(ZDU2(JL)) /ZMGEOM(JL)*RG

      ZTMP2(JL-KIDIA+1)= PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)&
                      & -PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
      ZDTV(JL)=( (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
       & +RETV*0.5_JPRB* (PQM1(JL,JK)-PQM1(JL,JK+1))&
       & *(PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)) ) * (1.0_JPRB/RCPD)
    ENDDO

    CALL VREC(ZTMP4,ZDU2(KIDIA),JLEN)
    CALL VREC(ZTMP5,ZTMP2,JLEN)

    DO JL=KIDIA,KFDIA
      ZDRORO= 2.0_JPRB * (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
       & *ZTMP5(JL-KIDIA+1)&
       & +RETV*(PQM1(JL,JK)-PQM1(JL,JK+1))  
      ZRI(JL)=ZMGEOM(JL)*ZDRORO*ZTMP4(JL-KIDIA+1)
    ENDDO

  ENDIF


  DO JL = KIDIA, KFDIA

!          DIMENSIONLESS COEFFICIENTS MULTIPLIED BY PRESSURE
!          THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE.

    ZZH     = 0.5_JPRB * ZRG * (PGEOM1(JL,JK)+PGEOM1(JL,JK+1)) + PZ0MM(JL)

!          STATICALLY STABLE

    IF ( ZRI(JL) > 0.0_JPRB ) THEN
!     ASYMPTOTIC MIXING LENGTH FOR STABLE SITUATIONS
      IF (PGEOM1(JL,JK)*ZRG <= ZPBLHEIGHT(JL)) THEN
      ZKLENT(JL,JK)=MAX(30.0_JPRB,0.1_JPRB*ZPBLHEIGHT(JL)*RPLRG)
      ZKLENT(JL,JK)=MIN(300.0_JPRB,ZKLENT(JL,JK))
      ELSE
      ZKLENT(JL,JK)=30.0_JPRB
      ENDIF
      ZKLENT(JL,JK)=ZKLENT(JL,JK)/RPLRG

      ZKLEN2  = RKAP * ZZH * ZKLENT(JL,JK) / ( RKAP * ZZH + ZKLENT(JL,JK))
!          COMPUTE STABILITY FUNCTIONS
      ZSCF    = SQRT(1.0_JPRB+ZCD*ZRI(JL))
      ZIFLTGM = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)/ZSCF) !F(LTG),M
      ZIFLTGH = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)*ZSCF) !F(LTG),H

      PCFM(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFLTGM
      PCFH(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFLTGH
      ZPHIM(JL) = 0.0_JPRB
      ZPHIH(JL) = 0.0_JPRB
!          STATICALLY UNSTABLE

    ELSE
!     ASYMPTOTIC MIXING LENGTH FOR UNSTABLE SITUATIONS
      ZKLENT(JL,JK)=RLAM
      ZKLEN2  = RKAP * ZZH * ZKLENT(JL,JK) / ( RKAP * ZZH + ZKLENT(JL,JK) )

!          COMPUTE STABILITY FUNCTIONS
      ZETA  = ZRI(JL)
      ZPHIM(JL) = PHIMU(ZETA)
      ZPHIH(JL) = PHIHU(ZETA)
!JJJ Forcing compiler to create reciprocal may not give most accurate computation, 
!JJJ but gives bit reproducibility (in T159) with previous 38r2 version
      ZIFMOM(JL)  = 1.0_JPRB / (ZPHIM(JL)**2)                              !F(MO),M
      ZIFMOH(JL)  = 1.0_JPRB / (ZPHIM(JL)*ZPHIH(JL))                       !F(MO),H
      PCFM(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFMOM(JL)
      PCFH(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFMOH(JL)
    ENDIF

!          ADD MIXED LAYER PARAMETRIZATION (SURF. AND OUTER LAYER)

!xmkIF (PGEOH(JL,JK-1)*ZRG <= PZINV(JL) ) THEN  ! up to level below entr. level
    IF (PGEOH(JL,JK-1)*ZRG <= PZINV(JL) .AND. ZKHVFL(JL)<0.0_JPRB) THEN  ! up to level below entr. level

      ZZ     = PGEOH(JL,JK)*ZRG
      ZDH    = ZZ/PZINV(JL)
      ZDL    = ZZ/ZL(JL)
      ZETA   = ZZ/ZL(JL)
      ZPHIKH = (1-39._JPRB*ZDL)**(-ZCONS13)
      ZPHIKM = (1-15._JPRB*ZDL)**(-ZCONS13)

!          K,surface

      PCFH(JL,JK) = RKAP / ZPHIKH * ZUST(JL) * ZZ * (1.0_JPRB-ZDH)**2
      PCFM(JL,JK) = RKAP / ZPHIKM * ZUST(JL) * ZZ * (1.0_JPRB-ZDH)**2

!          add cloud-top driven K (Lock et al. 2000, MWR p3187f, equ. 5)
!          (using simplified radiative velocity scale as in Lock, 1998, equ. 12
!          and ignore buoyancy reversal velocity scale)
!          apply K-cloud throughout full PBL (only for stratocumulus)

      ZML   = PZINV(JL)                               ! mixing depth
      ZBASE = 0.0_JPRB                                ! mixing base

      IF ( KPBLTYPE(JL) == 2  .AND.  ZZ >= ZBASE  .AND.  ZZ <= PZINV(JL) ) THEN  
        ZVSC  = ( RG / PTM1(JL,JK) * ZML * ZDRADFLX(JL) ) ** ZCONS13 
        ZKCLD = 0.85_JPRB * RKAP * ZVSC &
            & * (ZZ-ZBASE) ** 2 / ZML &
            & * ( 1 - (ZZ-ZBASE) / ZML ) ** 0.5_JPRB  
        PCFH(JL,JK) = PCFH(JL,JK) + ZKCLD
        PCFM(JL,JK) = PCFM(JL,JK) + ZKCLD * 0.75_JPRB

      ENDIF

    ENDIF

!          ADD ENTRAINMENT
!          entrainment velocity * T,v jump due to lw-radiative cooling & sfc flux
!          (Lock & MacVean, 1999, equ. 11)

    IF ( PGEOH(JL,JK)*ZRG <= PZINV(JL)  .AND.  PZINV(JL) < PGEOH(JL,JK-1)*ZRG ) THEN
      IF ( KPBLTYPE(JL) == 3 ) THEN                     !no entrainment below cumulus
        ZWTVENTR = 0.0_JPRB
      ELSE
        ZWTVENTR = - ZENTRSFC * ZKHVFL(JL)              !sfc flux (little impact)
        IF ( KPBLTYPE(JL) == 2 ) THEN
          ZWTVENTR = ZWTVENTR + ZENTRRAD * ZDRADFLX(JL) !radiation flux jump
        ENDIF
      ENDIF
      ZWTVENTR = MAX(0.0_JPRB,ZWTVENTR)

      ZKH     = ZWTVENTR * ZMGEOM(JL) / SIGN(MAX(ABS(RG*ZDTV(JL)),ZEPSILON),RG*ZDTV(JL))
      ZKH     = MAX(0.0_JPRB,ZKH)
      ZCFHNEW = ZKH

      PCFH(JL,JK) = MAX(PCFH(JL,JK),ZCFHNEW)            !protection against K=0
      PCFM(JL,JK) = MAX(PCFM(JL,JK),ZCFHNEW * 0.75_JPRB)

    ENDIF

!          LIMIT K TO 1000000 M2/S FOR SAFETY

    PCFH(JL,JK) = MIN( PCFH(JL,JK), 1000000.0_JPRB )
    PCFM(JL,JK) = MIN( PCFM(JL,JK), 1000000.0_JPRB )

!! MK: not yet implemented
!!
!!     IF(YREGWD%LRDIFF_STRATO) THEN
!!       IF(YREGWD%NDIFF_STRATO.EQ.1) THEN
!! ! Reduced diffusion in stratosphere just a function of pressure (quadratic from 80-120hPa)
!!         IF (ZRI(JL) > 0.25_JPRB) THEN
!!           IF(PAPHM1(JL,JK)>8000.0_JPRB) THEN
!!             ZRA=(PAPHM1(JL,JK)-8000.0_JPRB)*(PAPHM1(JL,JK)-8000.0_JPRB)/16.0E6_JPRB
!!             ZRA=MAX(MIN(1.0_JPRB,ZRA),0.002_JPRB)
!!           ELSE
!!             ZRA=0.002_JPRB
!!           ENDIF
!!           PCFH(JL,JK) = PCFH(JL,JK)*ZRA
!!           PCFM(JL,JK) = PCFM(JL,JK)*ZRA
!!         ENDIF
!!       ELSEIF(YREGWD%NDIFF_STRATO.EQ.2) THEN
!! ! Reduced diffusion in stratosphere only in tropical lower stratosphere
!!         IF(ABS(PGELAT(JL))<0.5_JPRB) THEN
!!         IF (ZRI(JL) > 0.25_JPRB) THEN
!!           IF(PAPHM1(JL,JK)>8000.0_JPRB) THEN
!!             ZRA=(PAPHM1(JL,JK)-8000.0_JPRB)*(PAPHM1(JL,JK)-8000.0_JPRB)/16.0E6_JPRB
!!           ELSEIF(PAPHM1(JL,JK)>1200.0_JPRB) THEN
!!             ZRA=0.0_JPRB
!!           ELSEIF(PAPHM1(JL,JK)>800.0_JPRB) THEN
!!             ZRA=(PAPHM1(JL,JK)-1200.0_JPRB)*(PAPHM1(JL,JK)-1200.0_JPRB)/16.0E4_JPRB
!!           ELSE
!!             ZRA=1.0_JPRB
!!           ENDIF
!!           ZRA=MAX(MIN(1.0_JPRB,ZRA),(1.0_JPRB+ZRI(JL))**(-4))
!!           PCFH(JL,JK) = PCFH(JL,JK)*ZRA
!!           PCFM(JL,JK) = PCFM(JL,JK)*ZRA
!!         ENDIF
!!         ENDIF
!!       ELSEIF(YREGWD%NDIFF_STRATO.EQ.5) THEN
!! ! Reduced diffusion in tropical lower stratosphere, smooth function of latitude
!!         IF (ZRI(JL) > 0.25_JPRB) THEN
!!           IF(PAPHM1(JL,JK)>12000.0_JPRB) THEN
!!             ZRA=1.0_JPRB
!!           ELSEIF(PAPHM1(JL,JK)>8000.0_JPRB) THEN
!!             ZRA=(PAPHM1(JL,JK)-8000.0_JPRB)*(PAPHM1(JL,JK)-8000.0_JPRB)/16.0E6_JPRB
!!           ELSEIF(PAPHM1(JL,JK)>1200.0_JPRB) THEN
!!             ZRA=0.0_JPRB
!!           ELSEIF(PAPHM1(JL,JK)>800.0_JPRB) THEN
!!             ZRA=(PAPHM1(JL,JK)-1200.0_JPRB)*(PAPHM1(JL,JK)-1200.0_JPRB)/16.0E4_JPRB
!!           ELSE
!!             ZRA=1.0_JPRB
!!           ENDIF
!! ! Restrict reduction to tropics (in extratropics, shear associated with wave activity and noise)
!!           ZLAT=ABS(PGELAT(JL))*57.296_JPRB
!!           IF(ZLAT.LT.20.0_JPRB) THEN
!!             ZRA=ZRA
!!           ELSEIF(ZLAT.LT.30.0_JPRB) THEN
!!             ZRA=1.0_JPRB-(1.0_JPRB-0.01_JPRB*(ZLAT-20.0_JPRB)*(ZLAT-20.0_JPRB))*(1.0_JPRB-ZRA)
!!           ELSE
!!             ZRA=1.0
!!           ENDIF
!!           ZRA=MAX(ZRA,(1.0_JPRB+ZRI(JL))**(-4))
!!           PCFH(JL,JK) = PCFH(JL,JK)*ZRA
!!           PCFM(JL,JK) = PCFM(JL,JK)*ZRA
!!         ENDIF
!!       ENDIF
!!     ENDIF
!!
!!xxx

!          DIFFUSION COEFFICIENT FOR HEAT FOR POSTPROCESSING ONLY IN (M2/S)

    PKH(JL,JK) = PCFH(JL,JK)*RPLDARE
    PKM(JL,JK) = PCFM(JL,JK)*RPLDARE

!          SCALE DIFFUSION COEFFICIENTS FOR USE IN VDFDIFH/M

!xmkZCFNC1 = RG * PAPHM1(JL,JK)&
    ZCFNC1 = RVDIFTS*PTMST*RG**2 * PAPHM1(JL,JK)&
         & /( 0.5_JPRB*RD * ZMGEOM(JL)&
         &    *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  ))&
         &    +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1)))) 
    PCFH(JL,JK) = PCFH(JL,JK) * ZCFNC1 * RPLDARE
    PCFM(JL,JK) = PCFM(JL,JK) * ZCFNC1 * RPLDARE

  ENDDO

!***
ENDDO
!***

IF (LHOOK) CALL DR_HOOK('VDFEXCU',1,ZHOOK_HANDLE)
END SUBROUTINE VDFEXCU_NEW


END MODULE mo_vdfexcu_new
