!>
!! Momentum solver for EDMF DUALM
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

MODULE mo_vdfdifm
 
  PUBLIC :: vdfdifm

CONTAINS

SUBROUTINE VDFDIFM(&
 & KIDIA  , KFDIA  , KLON   , KLEV   , KDRAFT , KTOP  , &
 & PTMST  , PUM1   , PVM1   , PAPHM1 , PCFM   , PMFLX , PUUH  , PVUH  , &
 & PTODC  , PSOTEU , PSOTEV , PSOC   , &
 & PVOM   , PVOL   , PUCURR , PVCURR , UMFL_S , VMFL_S, &
 & PUDIF  , PVDIF  , PTAUX  , PTAUY)  
!     ------------------------------------------------------------------

!**   *VDFDIFM* - DOES THE IMPLICIT CALCULATION FOR DIFFUSION OF MOMENTUM

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.    10-11-89
!     MODIFIED A BELJAARS   E.C.M.W.F.    17-11-02 Turb. Orogr. Drag  
!     MODIFIED A BELJAARS   E.C.M.W.F.    30-09-05 Include Subgr. Oro. in solver  
!     modified A.Beljaars 12/11/02   Ocean current b.c.
!     Introduction of mass flux transport   A. Beljaars 22-04-2005.
!     Multiple mass fluxes                  M. Ko"hler     07-2005.
!     Upwind discretization                 M. Ko"hler     10-2008.
!     modified P.de Rosnay/G.Balsamo      07-03-09 Offline Jacobians EKF 

!     PURPOSE
!     -------

!     SOLVE TRIDIAGONAL MATRICES FOR DIFFUSION OF MOMENTUM.
!     UPWIND.

!     INTERFACE
!     ---------

!     *VDFDIFM* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEV*         NUMBER OF LEVELS
!     *KDRAFT*       NUMBER OF UP/DOWN-DRAFTS
!     *KTOP*         INDEX FOR BOUNDARY LAYER TOP

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!     *PCFM*         PROP.TO EXCH.COEFF. FOR MOMENTUM (C,K-STAR IN DOC.)
!     *PMFLX*        MASSFLUX
!     *PUUH*         UPDRAFT X-MOMENTUM
!     *PVUH*         UPDRAFT Y-MOMENTUM
!     *PTODC*        TURB. OROGR. DRAG (IMPLICIT) COEFFICIENTS 
!     *PSOTEU*       Explicit part of U-tendency from subgrid orography scheme    
!     *PSOTEV*       Explicit part of V-tendency from subgrid orography scheme     
!     *PSOC*         Implicit part of subgrid orography (df/dt=PSOTE-PSOC*f/alpha)
!     *PVOM*         U-TENDENCY
!     *PVOL*         V-TENDENCY
!     *PUCURR*       OCEAN CURRENT X-COMPONENT
!     *PVCURR*       OCEAN CURRENT Y-COMPONENT

!     OUTPUT PARAMETERS (REAL):

!     *PUDIF*        U-DOUBLE-TILDE DEVIDED BY ALFA                [m/s]
!     *PVDIF*        V-DOUBLE-TILDE DEVIDED BY ALFA                [m/s]
!     *PTAUX*        SURFACE STRESS X-COMPONENT                    [N/m2]
!     *PTAUY*        SURFACE STRESS Y-COMPONENT                    [N/m2]

!     METHOD
!     ------

!     *LU*-DECOMPOSITION AND BACK SUBSTITUTION IN ONE DOWNWARD SCAN
!     AND ONE UPWARD SCAN.

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
! USE YOEVDF   , ONLY : RVDIFTS
! USE YOMCST   , ONLY : RG
! USE YOMSEKF  , ONLY : N_SEKF_PT, LUSEKF_REF, LUSE_JATM

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RG       ,&                                         !yomcst
                & RVDIFTS                                             !yoevdf
USE mo_edmf_param   ,ONLY : EDMF_CONF

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFLX(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUUH (KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVUH (KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTODC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: UMFL_S(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: VMFL_S(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUX(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUY(KLON) 

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZMFLX(KLON,0:KLEV), ZMU(KLON,0:KLEV), ZMV(KLON,0:KLEV)
REAL(KIND=JPRB) ::    ZAA(KLON,KLEV) , ZBB(KLON,KLEV)  , ZCC(KLON,KLEV)  ,&
                    & ZUYY(KLON,KLEV), ZVYY(KLON,KLEV) , ZGAM(KLON,KLEV),&
                    & Z1DP(KLON,KLEV), Z1BET(KLON)     ,&
                    & ZBUL(KLON)     , ZBVL(KLON)

INTEGER(KIND=JPIM) :: JK, JL, JD

REAL(KIND=JPRB) ::    ZTPFAC2, ZCOEF2, ZCOEF3, Z1DCOEF, ZAUVL
REAL(KIND=JPRB) ::    ZHOOK_HANDLE


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFDIFM',0,ZHOOK_HANDLE)

ZTPFAC2=1.0_JPRB/RVDIFTS

!*         1.0    PRECALCULATION OF MULTIPLE MASS-FLUX TERMS

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    ZMFLX(JL,JK) = 0.0_JPRB
    ZMU  (JL,JK) = 0.0_JPRB
    ZMV  (JL,JK) = 0.0_JPRB
  ENDDO
ENDDO
DO JD=2,KDRAFT  !don't include updraft #1 (test parcel)
  DO JK=KTOP,KLEV-1
    DO JL=KIDIA,KFDIA
!xmk
!     IF ( PMFLX(JL,JK,JD) > 0.0_JPRB .AND. PMFLX(JL,JK-1,JD) > 0.0_JPRB ) THEN
!       ZMFLX(JL,JK) = ZMFLX(JL,JK) + ( PMFLX(JL,JK,JD) + PMFLX(JL,JK-1,JD) ) / 2.0_JPRB
!       ZMU  (JL,JK) = ZMU  (JL,JK) + ( PMFLX(JL,JK,JD) + PMFLX(JL,JK-1,JD) ) &
!                                 & * ( PUUH(JL,JK,JD)  + PUUH (JL,JK-1,JD) ) / 4.0_JPRB
!       ZMV  (JL,JK) = ZMV  (JL,JK) + ( PMFLX(JL,JK,JD) + PMFLX(JL,JK-1,JD) ) &
!                                 & * ( PVUH(JL,JK,JD)  + PVUH (JL,JK-1,JD) ) / 4.0_JPRB
!     ENDIF
!---
!new: interpolate flux M*phi from half to full levels NOT M and phi separately as 32r3!
      ZMFLX(JL,JK) = ZMFLX(JL,JK) + ( PMFLX(JL,JK  ,JD) + PMFLX(JL,JK-1,JD) ) / 2.0_JPRB
      ZMU  (JL,JK) = ZMU  (JL,JK) + ( PMFLX(JL,JK  ,JD) * PUUH (JL,JK  ,JD) &
                                  & + PMFLX(JL,JK-1,JD) * PUUH (JL,JK-1,JD) ) / 2.0_JPRB
      ZMV  (JL,JK) = ZMV  (JL,JK) + ( PMFLX(JL,JK  ,JD) * PVUH (JL,JK  ,JD) &
                                  & + PMFLX(JL,JK-1,JD) * PVUH (JL,JK-1,JD) ) / 2.0_JPRB
!xxx
    ENDDO
  ENDDO
ENDDO

!*         1.1    SETTING OF THE MATRIX A, B AND C.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    Z1DP(JL,JK)=1.0_JPRB/(PAPHM1(JL,JK)-PAPHM1(JL,JK-1))
    ZAA(JL,JK) =(-PCFM(JL,JK-1)-ZMFLX(JL,JK-1))*Z1DP(JL,JK)
    ZCC(JL,JK) =(-PCFM(JL,JK)                 )*Z1DP(JL,JK)
    ZBB(JL,JK) =1.0_JPRB+(PCFM(JL,JK-1)+PCFM(JL,JK) &
     & +ZMFLX(JL,JK))*Z1DP(JL,JK) &
     & +PTODC(JL,JK)+PSOC(JL,JK)
  ENDDO
ENDDO

!          1.1a   THE SURFACE BOUNDARY CONDITION

DO JL=KIDIA,KFDIA
  Z1DP(JL,KLEV)=1.0_JPRB/(PAPHM1(JL,KLEV)-PAPHM1(JL,KLEV-1))
  ZCC(JL,KLEV) =0.0_JPRB
  ZAA(JL,KLEV) =        (-PCFM(JL,KLEV-1)-ZMFLX(JL,KLEV-1))*Z1DP(JL,KLEV)
  ZBB(JL,KLEV) =1.0_JPRB+(PCFM(JL,KLEV-1)                 )*Z1DP(JL,KLEV) &
               +PTODC(JL,KLEV)+PSOC(JL,KLEV)
ENDDO

!          1.1b   THE TOP BOUNDARY CONDITION    

DO JL=KIDIA,KFDIA
  Z1DP(JL,KTOP)=1.0_JPRB/(PAPHM1(JL,KTOP)-PAPHM1(JL,KTOP-1))
  ZAA(JL,KTOP) =0.0_JPRB
  ZCC(JL,KTOP) =         (-PCFM(JL,KTOP)               )*Z1DP(JL,KTOP)
  ZBB(JL,KTOP) =1.0_JPRB+( PCFM(JL,KTOP)+ZMFLX(JL,KTOP))*Z1DP(JL,KTOP) &
              &+PTODC(JL,KTOP)+PSOC(JL,KTOP)
ENDDO

!*         1.2    SETTING OF RIGHT HAND SIDES.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZUYY(JL,JK) = ZTPFAC2 * PUM1(JL,JK) & 
     & + PTMST * (PVOM(JL,JK)+PSOTEU(JL,JK)) &
     & + ZTPFAC2 * (ZMU(JL,JK)-ZMU(JL,JK-1)) * Z1DP(JL,JK)
    ZVYY(JL,JK) = ZTPFAC2 * PVM1(JL,JK) &
     & + PTMST * (PVOL(JL,JK)+PSOTEV(JL,JK)) &
     & + ZTPFAC2 * (ZMV(JL,JK)-ZMV(JL,JK-1)) * Z1DP(JL,JK)
  ENDDO
ENDDO

!          1.2a   SURFACE

JK=KLEV
DO JL=KIDIA,KFDIA
  ZUYY(JL,JK) = ZTPFAC2 * PUM1(JL,JK) &
   & + PTMST * (PVOM(JL,JK)+PSOTEU(JL,JK)) &
   & - ZTPFAC2 * ZMU(JL,JK-1) * Z1DP(JL,JK)
  ZVYY(JL,JK) = ZTPFAC2 * PVM1(JL,JK) &
   & + PTMST * (PVOL(JL,JK)+PSOTEV(JL,JK)) &
   & - ZTPFAC2 * ZMV(JL,JK-1) * Z1DP(JL,JK)
ENDDO

!          1.2b   TOP

JK=KTOP
DO JL=KIDIA,KFDIA
  ZUYY(JL,JK) = ZTPFAC2 * PUM1(JL,JK) &
   & + PTMST * (PVOM(JL,JK)+PSOTEU(JL,JK)) &
   & + ZTPFAC2 * ZMU(JL,JK) * Z1DP(JL,JK)
  ZVYY(JL,JK) = ZTPFAC2 * PVM1(JL,JK) &
   & + PTMST * (PVOL(JL,JK)+PSOTEV(JL,JK)) &
   & + ZTPFAC2 * ZMV(JL,JK) * Z1DP(JL,JK)
ENDDO

!*         1.3     TOP LAYER ELIMINATION.

DO JL=KIDIA,KFDIA
  Z1BET(JL)=1.0_JPRB/ZBB(JL,KTOP)
  PUDIF(JL,KTOP)=ZUYY(JL,KTOP)*Z1BET(JL)
  PVDIF(JL,KTOP)=ZVYY(JL,KTOP)*Z1BET(JL)
ENDDO

!*         1.4     ELIMINATION FOR MIDDLE LAYERS.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZGAM(JL,JK)=ZCC(JL,JK-1)*Z1BET(JL)
    Z1BET(JL)=1.0_JPRB/(ZBB(JL,JK)-ZAA(JL,JK)*ZGAM(JL,JK))
    PUDIF(JL,JK)=(ZUYY(JL,JK)-ZAA(JL,JK)*PUDIF(JL,JK-1))*Z1BET(JL)
    PVDIF(JL,JK)=(ZVYY(JL,JK)-ZAA(JL,JK)*PVDIF(JL,JK-1))*Z1BET(JL)
  ENDDO
ENDDO

!*         1.5     BOTTOM LAYER ELIMINATION.

ZCOEF2=1.0_JPRB/RVDIFTS
ZCOEF3=RG*PTMST*RVDIFTS

DO JL=KIDIA,KFDIA
  ZGAM(JL,KLEV)=ZCC(JL,KLEV-1)*Z1BET(JL)
  Z1BET(JL)=1.0_JPRB/(ZBB(JL,KLEV)-ZAA(JL,KLEV)*ZGAM(JL,KLEV))

!                  Compute Coefficients A and B (Psi_klev=A*J_Psi + B) 
!                  for U and V (ZAUVL=AU=AV)

  ZBUL(JL)=(ZUYY(JL,KLEV)-ZAA(JL,KLEV)*PUDIF(JL,KLEV-1))*RVDIFTS*Z1BET(JL)
  ZBVL(JL)=(ZVYY(JL,KLEV)-ZAA(JL,KLEV)*PVDIF(JL,KLEV-1))*RVDIFTS*Z1BET(JL)
  ZAUVL=-ZCOEF3*Z1DP(JL,KLEV)*Z1BET(JL)

  Z1DCOEF=1._JPRB/(ZCOEF3-PCFM(JL,KLEV)*ZAUVL)

  IF ( EDMF_CONF == 1 ) THEN
    PTAUX(JL)=PCFM(JL,KLEV)*(ZBUL(JL)-PUCURR(JL))*Z1DCOEF
    PTAUY(JL)=PCFM(JL,KLEV)*(ZBVL(JL)-PVCURR(JL))*Z1DCOEF
  ELSE
    PTAUX(JL)=UMFL_S(JL)
    PTAUY(JL)=VMFL_S(JL)
  ENDIF
  PUDIF(JL,KLEV)=(PTAUX(JL)*ZAUVL+ZBUL(JL))*ZCOEF2
  PVDIF(JL,KLEV)=(PTAUY(JL)*ZAUVL+ZBVL(JL))*ZCOEF2
ENDDO

!*         1.6     BACK-SUBSTITUTION.

DO JK=KLEV-1,KTOP,-1
  DO JL=KIDIA,KFDIA
    PUDIF(JL,JK)=PUDIF(JL,JK)-ZGAM(JL,JK+1)*PUDIF(JL,JK+1)
    PVDIF(JL,JK)=PVDIF(JL,JK)-ZGAM(JL,JK+1)*PVDIF(JL,JK+1)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('VDFDIFM',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFM


END MODULE mo_vdfdifm
