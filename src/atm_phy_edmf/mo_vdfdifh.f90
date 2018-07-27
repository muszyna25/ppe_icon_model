!>
!! Scalar solver for EDMF DUALM
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
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

MODULE mo_vdfdifh
 
  PUBLIC :: vdfdifh

CONTAINS

SUBROUTINE VDFDIFH(&
 & KIDIA  , KFDIA  , KLON   , KLEV   , KDRAFT , KTOP   , KTILES, &
 & PTMST  , PEXTSHF, PEXTLHF, LDSFCFLX,LDTERRA, LDLAND , &
 & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
 & PHLICE , PTLICE , PTLWML , &
 & PSLM1  , PTM1   , PQM1   , PQTM1  , PAPHM1 , &
 & PCFH   , PCFHTI , PCFQTI , PMFLX  , PSLUH  , PQTUH  , &
 & PTDIF  , PQDIF  , PCPTSTI, PQSTI  , PCAIRTI, PCSATTI, &
 & PDQSTI , PTSKTI , PTSKRAD, PTSM1M , PTSNOW , PTICE  , PSST, &
 & PTSKTIP1,PSLGE  , PTE    , PQTE   , &
 & PJQ    , PSSH   , PSLH   , PSTR   , PG0)  
!     ------------------------------------------------------------------

!**   *VDFDIFH* - DOES THE IMPLICIT CALCULATION FOR DIFFUSION OF S. L.

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.    10-11-89

!     OBUKHOV-L UPDATE      ACMB          26/03/90.
!     SKIN T CLEANING       P. VITERBO    15-11-96.
!     TILE BOUNDARY COND.   ACMB          20-11-98.
!     Surface DDH for TILES P. Viterbo    17-05-2000.
!     New tile coupling, 
!     DDH moved to VDFMAIN  A. Beljaars   2-05-2003.
!     Mass flux terms,
!     Flux b.c. for SCM,
!     Moist generalization  A. Beljaars/M. Ko"hler 3-12-2004.    
!     Multiple mass fluxes  M. Ko"hler               05-2005.
!     Removed option for linearized    P. Lopez   02/06/2005.
!     physics (now called separately) 
!     Upwind discretization M. Ko"hler               10-2008.  
!     Add lake tile                    G. Balsamo 18/04/2008
!     Offline Jacobians EKF P.de Rosnay/G.Balsamo 07/03/2009

!     PURPOSE
!     -------

!     SOLVE TRIDIAGONAL MATRICES FOR DIFFUSION OF DRY STATIC ENERGY
!     AND MOISTURE; IN SO DOING, IT ALSO SOLVES THE SKIN TEMPERATURE
!     EQUATION. UPWIND.

!     INTERFACE
!     ---------

!     *VDFDIFH* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEV*         NUMBER OF LEVELS
!     *KDRAFT*       NUMBER OF UP/DOWN-DRAFTS
!     *KTOP*         INDEX FOR BOUNDARY LAYER TOP

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PFRTI*        FRACTION OF SURFACE AREA COVERED BY TILES
!     *PSLM1*        GENERALIZED LIQUID WATER STATIC ENERGY    AT T-1
!                    (NOTE: In lin/adj physics = Dry static energy)
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMIDITY      AT T-1
!     *PQTM1*        SPECIFIC TOTAL WATER   AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT (C,K-STAR IN DOC.)
!     *PCFHTI*       IDEM FOR HEAT (SURFACE LAYER ONLY)
!     *PCFQTI*       IDEM FOR MOISTURE (SURFACE LAYER ONLY)
!     *PMFLX*        MASSFLUX
!     *PSLUH*        UPDRAFT GENERALIZED LIQUID WATER STATIC ENERGY AT HALF LEVEL
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL
!     *PCPTSTI*      DRY STATIC ENRGY AT SURFACE
!     *PQSTI*        SATURATION Q AT SURFACE
!     *PCAIRTI*      MULTIPLICATION FACTOR FOR Q AT LOWEST MODEL LEVEL
!                    FOR SURFACE FLUX COMPUTATION
!     *PCSATTI*      MULTIPLICATION FACTOR FOR QS AT SURFACE
!                    FOR SURFACE FLUX COMPUTATION
!     *PDQSTI*       D/DT (PQS)
!     *PSSRFLTI*     NET SOLAR RADIATION AT THE SURFACE, FOR EACH TILE
!     *PSLRFL*       NET THERMAL RADIATION AT THE SURFACE
!     *PEMIS*        MODEL SURFACE LONGWAVE EMISSIVITY
!     *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST
!     *PTLICE*       LAKE ICE TEMPERATURE [K]
!     *PHLICE*       LAKE ICE THICKNESS   [m]
!     *PTLWML*       LAKE MEAN WATER TEMPERATURE [K]
!     *PTSKTI*       SKIN TEMPERATURE AT T-1
!     *PTSKRAD*      SKIN TEMPERATURE OF LATEST FULL RADIATION TIMESTEP
!     *PTSM1M*       TOP SOIL LAYER TEMPERATURE
!     *PTSNOW*       SNOW TEMPERATURE 
!     *PTICE*        ICE TEMPERATURE (TOP SLAB)
!     *PSST*         (OPEN) SEA SURFACE TEMPERATURE
!     *PSLGE*        GENERALIZED DRY STATIC ENERGY TENDENCY
!                    (NOTE: In lin/adj physics = Temperature tendency)
!     *PQTE*         TOTAL WATER TENDENCY
!                    (NOTE: In lin/adj physics = Humidity tendency)

!     OUTPUT PARAMETERS (REAL):

!     *PTDIF*        SLG-DOUBLE-TILDE DEVIDED BY ALFA
!     *PQDIF*        QT-DOUBLE-TILDE DEVIDED BY ALFA
!     *PTSKTIP1*     SKIN TEMPERATURE AT T+1
!     *PJQ*          Surface moisture flux                      (kg/m2s)
!     *PSSH*         Surface sensible heat flux                 (W/m2)
!     *PSLH*         Surface latent heat flux                   (W/m2)
!     *PSTR*         Surface net thermal radiation              (W/m2)
!     *PG0*          Surface ground heat flux (solar radiation  (W/m2)
!                    leakage is not included in this term)

!     Additional parameters for flux boundary condtion (in SCM model):

!     *LDSFCFLX*     If .TRUE. flux boundary condtion is used 
!     *LDTERRA*      If .TRUE. flux boundary condition with TERRA fluxes is used 
!     *LDLAND*       Land mask
!     *PEXTSHF*      Specified sensible heat flux (W/m2)
!     *PEXTLHF*      Specified latent heat flux (W/m2)

!     METHOD
!     ------

!     *LU*-DECOMPOSITION (DOWNWARD SCAN), FOLLOWED BY SKIN-TEMPERATURE
!     SOLVER, AND BACK SUBSTITUTION (UPWARD SCAN).

!     EXTERNALS.
!     ----------

!     *VDFDIFH* CALLS:
!         *SURFSEB*

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
! USE YOEVDF   , ONLY : RVDIFTS
! USE YOMCST   , ONLY : RCPD     ,RG     ,RLSTT  ,RLVTT
! USE YOETHF   , ONLY : RVTMP2
! USE YOEPHY   , ONLY : LEOCWA   ,LEOCCO ,LEFLAKE
! USE YOMSEKF  , ONLY : N_SEKF_PT, LUSEKF_REF,LUSE_JATM

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RCPD     ,RG       ,RLVTT  ,&                       !yomcst
                & RVTMP2   ,&                                         !yoethf
                & RVDIFTS                                             !yoevdf
USE mo_edmf_param   ,ONLY : &
                & LEOCWA   ,LEOCCO   ,LEFLAKE  ,&                     !yoephy 
                & N_SEKF_PT          ,LUSEKF_REF         ,LUSE_JATM, &!yomsekf
                & EDMF_CONF
USE mo_surfseb      ,ONLY : surfseb

IMPLICIT NONE

!INTERFACE
! #include "surfseb.h"
!END INTERFACE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTSHF(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTLHF(KLON) 
LOGICAL           ,INTENT(IN)    :: LDSFCFLX 
LOGICAL           ,INTENT(IN)    :: LDTERRA 
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFLTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPSNW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFHTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFQTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFLX(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAIRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSATTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDQSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSM1M(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSKTIP1(KLON,KTILES)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PJQ(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSH(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLH(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTR(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PG0(KLON,KTILES) 

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZMFLX(KLON,0:KLEV), ZMSL(KLON,0:KLEV), ZMQT(KLON,0:KLEV)
REAL(KIND=JPRB) ::    ZAA(KLON,KLEV) ,ZBB(KLON,KLEV) ,ZCC(KLON,KLEV) ,&
                    & ZTYY(KLON,KLEV),ZQYY(KLON,KLEV),ZGAM(KLON,KLEV),&
                    & Z1DP(KLON,KLEV)
REAL(KIND=JPRB) ::    Z1BET(KLON)    ,&
                    & ZAQL(KLON)     ,ZBQL(KLON)     ,ZASL(KLON)     ,&
                    & ZBSL(KLON)     ,ZSL(KLON)      ,ZQL(KLON)
REAL(KIND=JPRB) ::    ZTSRF(KLON,KTILES)  ,ZRHOCHU(KLON,KTILES)  ,& 
                    & ZRHOCQU(KLON,KTILES),ZJS(KLON,KTILES)      ,&
                    & ZSSK(KLON,KTILES)   ,ZTSK(KLON,KTILES)  

INTEGER(KIND=JPIM) :: JK, JL, JT, JD

REAL(KIND=JPRB) ::    ZQDP, ZTPFAC2, ZTPFAC3,&
                    & ZCSNQ, ZCSNS, ZCOEF1,ZCOEF2
REAL(KIND=JPRB) ::    ZHOOK_HANDLE


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFDIFH',0,ZHOOK_HANDLE)

ZTPFAC2=1.0_JPRB/RVDIFTS
ZTPFAC3=1-ZTPFAC2

IF(.NOT.LUSEKF_REF.OR.LUSE_JATM.OR.N_SEKF_PT==0) THEN
! Coupled surface 
!     ------------------------------------------------------------------

!*         1.     FULL MODEL PHYSICS WITH MOIST MASS FLUX PBL
!                 -------------------------------------------

!*         1.0    PRECALCULATION OF MULTIPLE MASS-FLUX TERMS

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    ZMFLX(JL,JK) = 0.0_JPRB
    ZMSL (JL,JK) = 0.0_JPRB
    ZMQT (JL,JK) = 0.0_JPRB
  ENDDO
ENDDO
DO JD=2,KDRAFT  !don't include updraft #1 (test parcel)
  DO JK=KTOP,KLEV-1
    DO JL=KIDIA,KFDIA
!xmk
!     IF ( PMFLX(JL,JK,JD) > 0.0_JPRB .AND. PMFLX(JL,JK-1,JD) > 0.0_JPRB ) THEN
!       ZMFLX(JL,JK) = ZMFLX(JL,JK) + ( PMFLX(JL,JK,JD) + PMFLX(JL,JK-1,JD) ) / 2.0_JPRB
!       ZMSL (JL,JK) = ZMSL (JL,JK) + ( PMFLX(JL,JK,JD) + PMFLX(JL,JK-1,JD) ) &
!                                 & * ( PSLUH(JL,JK,JD) + PSLUH(JL,JK-1,JD) ) / 4.0_JPRB
!       ZMQT (JL,JK) = ZMQT (JL,JK) + ( PMFLX(JL,JK,JD) + PMFLX(JL,JK-1,JD) ) &
!                                 & * ( PQTUH(JL,JK,JD) + PQTUH(JL,JK-1,JD) ) / 4.0_JPRB
!     ENDIF
!---
!new: interpolate flux M*phi from half to full levels NOT M and phi separately as 32r3!
      ZMFLX(JL,JK) = ZMFLX(JL,JK) + ( PMFLX(JL,JK  ,JD) + PMFLX(JL,JK-1,JD) ) / 2.0_JPRB
      ZMSL (JL,JK) = ZMSL (JL,JK) + ( PMFLX(JL,JK  ,JD) * PSLUH(JL,JK  ,JD) &
                                  & + PMFLX(JL,JK-1,JD) * PSLUH(JL,JK-1,JD) ) / 2.0_JPRB
      ZMQT (JL,JK) = ZMQT (JL,JK) + ( PMFLX(JL,JK  ,JD) * PQTUH(JL,JK  ,JD) &
                                  & + PMFLX(JL,JK-1,JD) * PQTUH(JL,JK-1,JD) ) / 2.0_JPRB
!xxx
    ENDDO
  ENDDO
ENDDO

!*         1.1    SETTING OF THE MATRIX A, B AND C.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    Z1DP(JL,JK)=1.0_JPRB/(PAPHM1(JL,JK)-PAPHM1(JL,JK-1))
    ZAA(JL,JK) =(-PCFH(JL,JK-1)-ZMFLX(JL,JK-1))*Z1DP(JL,JK)
    ZCC(JL,JK) =(-PCFH(JL,JK)                 )*Z1DP(JL,JK)
    ZBB(JL,JK) =1.0_JPRB+(PCFH(JL,JK-1)+PCFH(JL,JK)&
     & +ZMFLX(JL,JK))*Z1DP(JL,JK)  
  ENDDO
ENDDO

!          1.1a   THE SURFACE BOUNDARY CONDITION

DO JL=KIDIA,KFDIA
  Z1DP(JL,KLEV)=1.0_JPRB/(PAPHM1(JL,KLEV)-PAPHM1(JL,KLEV-1))
  ZCC(JL,KLEV) =0.0_JPRB
  ZAA(JL,KLEV) =        (-PCFH(JL,KLEV-1)-ZMFLX(JL,KLEV-1))*Z1DP(JL,KLEV)
  ZBB(JL,KLEV) =1.0_JPRB+(PCFH(JL,KLEV-1)                 )*Z1DP(JL,KLEV)  
ENDDO

!          1.1b   THE TOP BOUNDARY CONDITION    

DO JL=KIDIA,KFDIA
  Z1DP(JL,KTOP)=1.0_JPRB/(PAPHM1(JL,KTOP)-PAPHM1(JL,KTOP-1))
  ZAA(JL,KTOP) =0.0_JPRB
  ZCC(JL,KTOP) =         (-PCFH(JL,KTOP)               )*Z1DP(JL,KTOP)
  ZBB(JL,KTOP) =1.0_JPRB+( PCFH(JL,KTOP)+ZMFLX(JL,KTOP))*Z1DP(JL,KTOP)
ENDDO

!*         1.2    SETTING OF RIGHT HAND SIDES.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZTYY(JL,JK) = ZTPFAC2 * PSLM1(JL,JK) &
     & + PTMST * PSLGE(JL,JK) &
     & + ZTPFAC2 * (ZMSL(JL,JK)-ZMSL(JL,JK-1)) * Z1DP(JL,JK)
    ZQYY(JL,JK) = ZTPFAC2 * PQTM1(JL,JK) &
     & + PTMST * PQTE(JL,JK) &
     & + ZTPFAC2 * (ZMQT(JL,JK)-ZMQT(JL,JK-1)) * Z1DP(JL,JK)
  ENDDO
ENDDO

!          1.2a   SURFACE

JK=KLEV
DO JL=KIDIA,KFDIA
  ZTYY(JL,JK) = ZTPFAC2 * PSLM1(JL,JK) &
   & + PTMST * PSLGE(JL,JK) &
   & - ZTPFAC2 * ZMSL(JL,JK-1) * Z1DP(JL,JK)
  ZQYY(JL,JK) = ZTPFAC2 * PQTM1(JL,JK) &
   & + PTMST * PQTE(JL,JK) &
   & - ZTPFAC2 * ZMQT(JL,JK-1) * Z1DP(JL,JK)
ENDDO

!          1.2b   TOP

JK=KTOP
DO JL=KIDIA,KFDIA
  ZTYY(JL,JK) = ZTPFAC2 * PSLM1(JL,JK) &
   & + PTMST * PSLGE(JL,JK) &
   & + ZTPFAC2 *ZMSL(JL,JK) *Z1DP(JL,JK)
  ZQYY(JL,JK) = ZTPFAC2 * PQTM1(JL,JK) &
   & + PTMST * PQTE(JL,JK) &
   & + ZTPFAC2 *ZMQT(JL,JK) *Z1DP(JL,JK)
ENDDO

!*         1.3    ADD MOISTURE FLUX FROM SNOW FROM TILE 7 AS EXPLICIT TERM

IF ( EDMF_CONF == 1 ) THEN  ! This cannot be done if ntiles taken from ICON)

  JK=KLEV
  DO JL=KIDIA,KFDIA
    ZCSNQ=RG*PTMST*PFRTI(JL,7)*PEVAPSNW(JL)*Z1DP(JL,JK)
    ZCSNS=RCPD*RVTMP2*PTSKTI(JL,7)*ZCSNQ
    ZTYY(JL,JK)=ZTYY(JL,JK)-ZCSNS
    ZQYY(JL,JK)=ZQYY(JL,JK)-ZCSNQ
  ENDDO

ENDIF

!*         1.4    TOP LAYER ELIMINATION.

DO JL=KIDIA,KFDIA
  Z1BET(JL)=1.0_JPRB/ZBB(JL,KTOP)
  PTDIF(JL,KTOP)=ZTYY(JL,KTOP)*Z1BET(JL)
  PQDIF(JL,KTOP)=ZQYY(JL,KTOP)*Z1BET(JL)
ENDDO

!*         1.5    ELIMINATION FOR MIDDLE LAYERS.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZGAM(JL,JK)=ZCC(JL,JK-1)*Z1BET(JL)
    Z1BET(JL)=1.0_JPRB/(ZBB(JL,JK)-ZAA(JL,JK)*ZGAM(JL,JK))
    PTDIF(JL,JK)=(ZTYY(JL,JK)-ZAA(JL,JK)*PTDIF(JL,JK-1))*Z1BET(JL)
    PQDIF(JL,JK)=(ZQYY(JL,JK)-ZAA(JL,JK)*PQDIF(JL,JK-1))*Z1BET(JL)
  ENDDO
ENDDO

!*         1.6    BOTTOM LAYER, LINEAR RELATION BETWEEN LOWEST
!                 MODEL LEVEL S AND Q AND FLUXES.

DO JL=KIDIA,KFDIA
  ZGAM(JL,KLEV)=ZCC(JL,KLEV-1)*Z1BET(JL)
  Z1BET(JL)=1.0_JPRB/(ZBB(JL,KLEV)-ZAA(JL,KLEV)*ZGAM(JL,KLEV))
  ZBSL(JL)=(ZTYY(JL,KLEV)-ZAA(JL,KLEV)*PTDIF(JL,KLEV-1))*RVDIFTS*Z1BET(JL)
  ZBQL(JL)=(ZQYY(JL,KLEV)-ZAA(JL,KLEV)*PQDIF(JL,KLEV-1))*RVDIFTS*Z1BET(JL)
  ZQDP=1.0_JPRB/(PAPHM1(JL,KLEV)-PAPHM1(JL,KLEV-1))
  ZASL(JL)=-RG*PTMST*ZQDP*RVDIFTS*Z1BET(JL)
  ZAQL(JL)=ZASL(JL)
ENDDO

ELSE
! Offline Jacobian setup
JK=KLEV
DO JL=KIDIA,KFDIA
  ZBSL(JL)=PSLM1(JL,JK) + PTMST*PSLGE(JL,JK)
  ZBQL(JL)=PQTM1(JL,JK) + PTMST* PQTE(JL,JK)
  ZASL(JL)=0.0_JPRB
  ZAQL(JL)=0.0_JPRB
ENDDO
ENDIF

!*         1.7    PREPARE ARRAY'S FOR CALL TO SURFACE ENERGY
!                 BALANCE ROUTINE

IF (edmf_conf == 1) THEN

  IF (LEOCWA .OR. LEOCCO) THEN
    ZTSRF(KIDIA:KFDIA,1)=PTSKTI(KIDIA:KFDIA,1)
  ELSE
  ZTSRF(KIDIA:KFDIA,1)=PSST(KIDIA:KFDIA)
  ENDIF
  ZTSRF(KIDIA:KFDIA,2)=PTICE(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,3)=PTSM1M(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,4)=PTSM1M(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,5)=PTSNOW(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,6)=PTSM1M(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,7)=PTSNOW(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,8)=PTSM1M(KIDIA:KFDIA)
  IF (LEFLAKE) THEN
    DO JL=KIDIA,KFDIA
      IF(PHLICE(JL) > 1.E-9_JPRB) THEN ! 1.E-9 or H_ICE_MIN_FLK present
        ZTSRF(JL,9)=PTLICE(JL)
      ELSE
        ZTSRF(JL,9)=PTLWML(JL)
      ENDIF
    ENDDO
  ENDIF

ENDIF

ZCOEF1=1.0_JPRB/(RG*RVDIFTS*PTMST)
DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
    ZRHOCHU(JL,JT)=PCFHTI(JL,JT)*ZCOEF1
    ZRHOCQU(JL,JT)=PCFQTI(JL,JT)*ZCOEF1
  ENDDO
ENDDO

!*         1.8    CALL TO SURFACE ENERGY BALANCE ROUTINE
!                 REMEMBER: OUTPUT IS EXTRAPOLATED IN TIME

IF (edmf_conf == 1) THEN

  CALL SURFSEB   (KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,KTILES=KTILES,&
  & PSSKM1M=PCPTSTI,PTSKM1M=PTSKTI,PQSKM1M=PQSTI,&
  & PDQSDT=PDQSTI,PRHOCHU=ZRHOCHU,PRHOCQU=ZRHOCQU,&
  & PALPHAL=PCAIRTI,PALPHAS=PCSATTI,&
  & PSSRFL=PSSRFLTI,PFRTI=PFRTI,PTSRF=ZTSRF,&
  & PHLICE=PHLICE,&
  & PSLRFL=PSLRFL,PTSKRAD=PTSKRAD,PEMIS=PEMIS,&
  & PASL=ZASL,PBSL=ZBSL,PAQL=ZAQL,PBQL=ZBQL,&
  !out
  & PJS=ZJS,PJQ=PJQ,PSSK=ZSSK,PTSK=ZTSK,&
  & PSSH=PSSH,PSLH=PSLH,PSTR=PSTR,PG0=PG0,&
  & PSL=ZSL,PQL=ZQL)

ENDIF

!*         1.9    ADD SNOW EVAPORATION TO FLUXES

!dmk  This needs to be done special !!!
!
!PJQ (KIDIA:KFDIA,7)=PJQ (KIDIA:KFDIA,7)+PEVAPSNW(KIDIA:KFDIA)
!PSLH(KIDIA:KFDIA,7)=PSLH(KIDIA:KFDIA,7)+PEVAPSNW(KIDIA:KFDIA)*RLSTT
!
!xxx

!*         1.10a   Flux boundary condition for 1D model (fluxes in W/m2)
!                 (Over-write output of SURFSEB)

IF (LDSFCFLX) THEN
  DO JT=1,KTILES
    DO JL=KIDIA,KFDIA
!xmk  ZJS(JL,JT)=PEXTSHF(JL)+RCPD*PTSKTI(JL,JT)*RVTMP2*PEXTLHF(JL)/RLVTT
      ZJS(JL,JT)=PEXTSHF(JL)  !no more RVTMP2
      PJQ(JL,JT)=PEXTLHF(JL)/RLVTT
        
      ZSSK(JL,JT)=ZBSL(JL)+ZJS(JL,JT)*(ZASL(JL)-1.0_JPRB/ZRHOCHU(JL,JT)) 
      ZTSK(JL,JT)=ZSSK(JL,JT)/(RCPD*(1.+RVTMP2*PQSTI(JL,JT)))
      PSSH(JL,JT)=PEXTSHF(JL)
      PSLH(JL,JT)=PEXTLHF(JL)
      PSTR(JL,JT)=PSLRFL(JL)
      PG0 (JL,JT)=PEXTSHF(JL)+PEXTLHF(JL)+PSLRFL(JL)+PSSRFLTI(JL,JT)
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    ZSL(JL)=ZJS(JL,1)*ZASL(JL)+ZBSL(JL)
    ZQL(JL)=PJQ(JL,1)*ZAQL(JL)+ZBQL(JL)
  ENDDO
ENDIF

!*         1.10b   Flux boundary condition from TERRA land only (fluxes in W/m2)
!                 (Over-write output of SURFSEB)

IF (edmf_conf == 1) THEN

  IF (LDTERRA) THEN
    DO JT=3,KTILES     ! TERRA goes in tiles 3-8
      DO JL=KIDIA,KFDIA
        IF ( LDLAND(JL) ) THEN
  !xmk    ZJS(JL,JT)=PEXTSHF(JL)+RCPD*PTSKTI(JL,JT)*RVTMP2*PEXTLHF(JL)/RLVTT
          ZJS(JL,JT)=PEXTSHF(JL)  !no more RVTMP2
          PJQ(JL,JT)=PEXTLHF(JL)/RLVTT                !RLVTT or RLSTT ???
          
  !       ZSSK(JL,JT)=ZBSL(JL)+ZJS(JL,JT)*(ZASL(JL)-1.0_JPRB/ZRHOCHU(JL,JT)) 
  !       ZTSK(JL,JT)=ZSSK(JL,JT)/(RCPD*(1.+RVTMP2*PQSTI(JL,JT)))
  
  !----here should be the mean TERRA TSK - maybe separate for snow and soil----
  !----same for fluxes - separate for snow and soil (from vdfmain) ---
  
          PSSH(JL,JT)=PEXTSHF(JL)
          PSLH(JL,JT)=PEXTLHF(JL)
          PSTR(JL,JT)=PSLRFL(JL)
          PG0 (JL,JT)=PEXTSHF(JL)+PEXTLHF(JL)+PSLRFL(JL)+PSSRFLTI(JL,JT)
        ENDIF
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      IF ( LDLAND(JL) ) THEN
        ZSL(JL)=PEXTSHF(JL)      *ZASL(JL)+ZBSL(JL)   ! s,l from land flux - tile average
        ZQL(JL)=PEXTLHF(JL)/RLVTT*ZAQL(JL)+ZBQL(JL)   ! q,l  -"-
      ENDIF
    ENDDO
  ENDIF
  
  !*         1.11   COMPUTE PARAMETERS AT NEW TIME LEVEL 
  
  !dmk  No skin on ICON
  DO JT=1,KTILES
    DO JL=KIDIA,KFDIA
      PTSKTIP1(JL,JT)=ZTPFAC2*ZTSK(JL,JT)+ZTPFAC3*PTSKTI(JL,JT)     !?????
    ENDDO
  ENDDO
  !xxx


ELSE   ! edmf_conf == 2


  DO JT=1,KTILES     ! all tiles
    DO JL=KIDIA,KFDIA
!     ZJS(JL,JT)=PEXTSHF(JL) 
      PJQ(JL,JT)=PEXTLHF(JL)/RLVTT                !RLVTT or RLSTT ???
        
!     ZSSK(JL,JT)=ZBSL(JL)+ZJS(JL,JT)*(ZASL(JL)-1.0_JPRB/ZRHOCHU(JL,JT)) 
!     ZTSK(JL,JT)=ZSSK(JL,JT)/(RCPD*(1.+RVTMP2*PQSTI(JL,JT)))
!     ZTSK(JL,JT)=PTSKTI(JL,JT)

      PSSH(JL,JT)=PEXTSHF(JL)
      PSLH(JL,JT)=PEXTLHF(JL)
!     PSTR(JL,JT)=PSLRFL(JL)                                           ! not needed
!     PG0 (JL,JT)=PEXTSHF(JL)+PEXTLHF(JL)+PSLRFL(JL)+PSSRFLTI(JL,JT)   ! not needed
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    ZSL(JL)=PEXTSHF(JL)      *ZASL(JL)+ZBSL(JL)   ! s,l from land flux - tile average
    ZQL(JL)=PEXTLHF(JL)/RLVTT*ZAQL(JL)+ZBQL(JL)   ! q,l  -"-
  ENDDO

  !*         1.11   COMPUTE PARAMETERS AT NEW TIME LEVEL 
  
  DO JT=1,KTILES
    DO JL=KIDIA,KFDIA
!     PTSKTIP1(JL,JT)=ZTPFAC2*ZTSK(JL,JT)+ZTPFAC3*PTSKTI(JL,JT)
      PTSKTIP1(JL,JT)=PTSKTI(JL,JT)        ! take last time step skin - no skin in ICON
    ENDDO
  ENDDO


ENDIF


!*         1.12   COPY LOWEST MODEL SOLUTION FROM SURFSEB

ZCOEF2=1.0_JPRB/RVDIFTS
DO JL=KIDIA,KFDIA
  PTDIF(JL,KLEV)=ZSL(JL)*ZCOEF2
  PQDIF(JL,KLEV)=ZQL(JL)*ZCOEF2
ENDDO

!*         1.13   BACK-SUBSTITUTION.

DO JK=KLEV-1,KTOP,-1
  DO JL=KIDIA,KFDIA
    PTDIF(JL,JK)=PTDIF(JL,JK)-ZGAM(JL,JK+1)*PTDIF(JL,JK+1)
    PQDIF(JL,JK)=PQDIF(JL,JK)-ZGAM(JL,JK+1)*PQDIF(JL,JK+1)
  ENDDO
ENDDO
  
IF (LHOOK) CALL DR_HOOK('VDFDIFH',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFH


END MODULE mo_vdfdifh
