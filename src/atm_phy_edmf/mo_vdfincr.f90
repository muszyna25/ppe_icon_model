!>
!! Diagnostic PBL height for EDMF DUALM
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

MODULE mo_vdfincr
 
  PUBLIC :: vdfincr

CONTAINS

!! #ifdef RS6K
!! @PROCESS HOT(NOVECTOR) 
!! #endif
SUBROUTINE VDFINCR(  KIDIA  , KFDIA  , KLON   , KLEV  , PTMST  , &
 & PUM1   , PVM1   , PSLGM1 , PQTM1  , PAPHM1 , &
 & PTODC  , PSOTEU , PSOTEV , PSOC   , &
 & PUDIF  , PVDIF  , PSLGDIF, PQTDIF , &
 & PVOM   , PVOL   , PSLGE  , PQTE   , PSLGEWODIS, &
 & PVDIS  , PVDISG, PSTRTU  , PSTRTV , PSTRSOU, PSTRSOV, PTOFDU , PTOFDV, &
 ! DIAGNOSTIC OUTPUT
 & PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , JCNT   , LLDIAG , &
 & PDISGW3D)  
!     ------------------------------------------------------------------

!**   *VDFINCR* - INCREMENTS U,V,T AND Q-TENDENCIES; COMPUTE MULTILEVEL
!                 FLUXES AND DISSIPATION.

!     A.C.M. BELJAARS  18/01/90   DERIVED FROM VDIFF (CY34)
!     A.C.M. BELJAARS  26/03/90   OBUKHOV-L UPDATE
!     M. Ko"hler        3/12/2004 Conserved variables (qt and slg)
!     A. Beljaars       4/04/2005 TURBULENT OROGR. DRAG ACMB
!     A  Beljaars      30/09/2005 Include Subgr. Oro. in solver  
!     OCEAN CURRENT B.C.    ACMB          12/11/02.
!     P. Lopez         02/06/2005 Removed option for linearized
!                                 physics (now called separately)
!     P.de Rosnay/G.Balsamo   07/03/2009 Offline Jacobians EKF      

!     PURPOSE
!     -------

!     INCREMENT U,V,T AND Q; COMPUTE MULTILEVEL FLUXES AND DISSIPATION

!     INTERFACE
!     ---------

!     *VDFINCR* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!!!   *KTOP*         FIRST LEVEL INDEX WITHOUT ZERO-DIFFUSION


!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PSLGM1*       GENERALIZED LIQUID WATER STATIC ENERGY (SLG) AT T-1
!!!   *PTM1*         TEMPERATURE AT T-1
!     *PQTM1*        TOTAL WATER AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!!!   *PGEOM1*       GEOPOTENTIAL AT T-1
!!!   *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!     *PTODC*        TURBULENT OROGRAPHIC DRAG COEFFICIENT
!     *PSOTEU*       Explicit part of U-tendency from subgrid orography scheme    
!     *PSOTEV*       Explicit part of V-tendency from subgrid orography scheme     
!     *PSOC*         Implicit part of subgrid orography (df/dt=PSOTE-PSOC*f/alpha)
!     *PUDIF*        U-DOUBLE TILDE DEVIDED BY ALFA
!     *PVDIF*        V-DOUBLE TILDE DEVIDED BY ALFA
!!!   *PUCURR*       U-OCEAN CURRENT
!!!   *PVCURR*       V-OCEAN CURRENT

!     UPDATED PARAMETERS (REAL):

!     *PSLGDIF*      SLG-DOUBLE TILDE DEVIDED BY ALFA (ON ENTRY)
!                    SLG-SINGLE TILDE                 (ON EXIT)
!     *PQTDIF*       QT-DOUBLE TILDE DEVIDED BY ALFA  (ON ENTRY)
!                    QT-SINGLE TILDE                  (ON EXIT)
!     *PVOM*         U-TENDENCY
!     *PVOL*         V-TENDENCY
!     *PSLGE*        SLG-TENDENCY
!     *PQTE*         QT-TENDENCY

!     OUTPUT PARAMETERS (REAL):

!     *PVDIS*        TURBULENT DISSIPATION
!     *PVDISG*       SUBGRID OROGRAPHY DISSIPATION
!     *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM         KG*(M/S)/(M2*S)
!     *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM         KG*(M/S)/(M2*S)
!     *PSTRSOU*      SUBGRID OROGRAPHY FLUX OF U-MOMEMTUM KG*(M/S)/(M2*S)
!     *PSTRSOV*      SUBGRID OROGRAPHY FLUX OF V-MOMEMTUM KG*(M/S)/(M2*S)
!     *PSLGEWODIS*   SLG-TENDENCY MINUS (TOTAL) DISSIPATION
!     *PTOFDU*       TOFD COMP. OF TURBULENT FLUX OF U-MOMEMTUM    KG*(M/S)/(M2*S)
!     *PTOFDV*       TOFD COMP. OF TURBULENT FLUX OF V-MOMEMTUM    KG*(M/S)/(M2*S)
!     *PDISGW3D*     3-D DISSIPATION RATE FOR STOCHASTIC PHYSIC    (M2/S2)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
! USE YOMCST   , ONLY : RG
! USE YOEVDF   , ONLY : RVDIFTS
! USE YOMSEKF  , ONLY : N_SEKF_PT, LUSEKF_REF, LUSE_JATM

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RG       ,&                                         !yomcst
                & RVDIFTS                                             !yoevdf
USE mo_edmf_param   ,ONLY : &
                & N_SEKF_PT          ,LUSEKF_REF         ,LUSE_JATM   !yomsekf

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
!! INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGM1(KLON,KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTODC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIF(KLON,KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLGDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQTDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGEWODIS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDISG(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISGW3D(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDV(KLON) 

!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX2, KLEVX, KFLDX, JCNT
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)
LOGICAL           ,INTENT(IN)    :: LLDIAG

!*         0.2    LOCAL VARIABLES

INTEGER(KIND=JPIM) :: JK, JL, ILEV

REAL(KIND=JPRB) ::    ZCONS1, ZCONS2, &
                    & ZDUDT, ZDVDT, ZVDFDIS, ZSODIS, ZTPFAC2, ZTPFAC3,&
                    & ZDP, ZRG, ZGDPH, ZTETOFDU,ZTETOFDV,&
                    & ZTESOU,ZTESOV,ZHU2,ZHU3,ZU1,ZU2,ZU3,ZV1,ZV2,ZV3
REAL(KIND=JPRB) ::    ZHOOK_HANDLE


!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFINCR',0,ZHOOK_HANDLE)
ZTPFAC2 = 1.0_JPRB/RVDIFTS
ZTPFAC3 = 1.0_JPRB-ZTPFAC2

ZCONS1  = 1.0_JPRB/PTMST
ZCONS2  = 1.0_JPRB/(RG*PTMST)

ZRG     = 1.0_JPRB/RG

IF (.NOT.LUSEKF_REF.OR.LUSE_JATM.OR.N_SEKF_PT==0) THEN
! Coupled surface
  ILEV=1 
ELSEIF (N_SEKF_PT>0) THEN
! Offline EKF surface analysis perturbed runs for the Jacobians
  ILEV=KLEV !loop calculation only on lowest model level 
ENDIF
!     ------------------------------------------------------------------

!*         2.    COMPUTE TENDENCIES AND BUDGETS
!                ------------------------------

DO JL=KIDIA,KFDIA
  PVDIS(JL)    =0.0_JPRB
  PVDISG(JL)   =0.0_JPRB
  PSTRTU(JL,0) =0.0_JPRB
  PSTRTV(JL,0) =0.0_JPRB
  PSTRSOU(JL,0)=0.0_JPRB
  PSTRSOV(JL,0)=0.0_JPRB
  PTOFDU(JL)   =0.0_JPRB
  PTOFDV(JL)   =0.0_JPRB
ENDDO

!*         2.1  VERTICAL LOOP

DO JK=ILEV,KLEV
  DO JL=KIDIA,KFDIA

!   Compute total tendencies (dynamics + vertical diffusion + SO) 
    ZDUDT          = ( PUDIF(JL,JK) - ZTPFAC2 * PUM1(JL,JK) ) * ZCONS1
    ZDVDT          = ( PVDIF(JL,JK) - ZTPFAC2 * PVM1(JL,JK) ) * ZCONS1
     
    ZDP            = PAPHM1(JL,JK)-PAPHM1(JL,JK-1)
    ZGDPH          = -ZDP*ZRG

!   TOFD-tendencies 
    ZHU2=-PTODC(JL,JK)*ZCONS1
    ZTETOFDU=PUDIF(JL,JK)*ZHU2
    ZTETOFDV=PVDIF(JL,JK)*ZHU2
    PTOFDU(JL)=PTOFDU(JL)+ZGDPH*ZTETOFDU
    PTOFDV(JL)=PTOFDV(JL)+ZGDPH*ZTETOFDV

!   Implicit part of SO-tendencies 
    ZHU3=-PSOC(JL,JK)*ZCONS1
    ZTESOU=PUDIF(JL,JK)*ZHU3
    ZTESOV=PVDIF(JL,JK)*ZHU3
    
!   Velocity before VDF 
    ZU1=PUM1(JL,JK)+PVOM(JL,JK)*PTMST
    ZV1=PVM1(JL,JK)+PVOL(JL,JK)*PTMST
!   Velocity after VDF 
    ZU2=PUM1(JL,JK)+(ZDUDT-PSOTEU(JL,JK)-ZTESOU)*PTMST
    ZV2=PVM1(JL,JK)+(ZDVDT-PSOTEV(JL,JK)-ZTESOV)*PTMST
!   Velocity after SO 
    ZU3=PUM1(JL,JK)+ZDUDT*PTMST   
    ZV3=PVM1(JL,JK)+ZDVDT*PTMST   
    ZVDFDIS=0.5_JPRB*(ZU1-ZU2)*(ZU1+ZU2) + 0.5_JPRB*(ZV1-ZV2)*(ZV1+ZV2)
    ZSODIS =0.5_JPRB*(ZU2-ZU3)*(ZU2+ZU3) + 0.5_JPRB*(ZV2-ZV3)*(ZV2+ZV3)
    PDISGW3D(JL,JK) = ZSODIS

!amk------------------------------------------------------------------
    if(lldiag) then
!     Output TOFD tendency
      PEXTRA(JL,JK,42)=PEXTRA(JL,JK,42)+ZTETOFDU*PTMST
      PEXTRA(JL,JK,43)=PEXTRA(JL,JK,43)+ZTETOFDV*PTMST
  
!     Output accumulated gravity wave drag tendency
      PEXTRA(JL,JK,44)=PEXTRA(JL,JK,44)+PSOTEU(JL,JK)*PTMST
      PEXTRA(JL,JK,45)=PEXTRA(JL,JK,45)+PSOTEV(JL,JK)*PTMST
  
!     Output accumulated low level blocking drag tendency
      PEXTRA(JL,JK,46)=PEXTRA(JL,JK,46)+ZTESOU*PTMST
      PEXTRA(JL,JK,47)=PEXTRA(JL,JK,47)+ZTESOV*PTMST
    endif
!xxx------------------------------------------------------------------


!   Integrate VDF-tendencies (including TOFD) to find VDF-stress profile
    PSTRTU(JL,JK)  = (ZDUDT-PVOM(JL,JK)-PSOTEU(JL,JK)-ZTESOU)*ZGDPH+PSTRTU(JL,JK-1)
    PSTRTV(JL,JK)  = (ZDVDT-PVOL(JL,JK)-PSOTEV(JL,JK)-ZTESOV)*ZGDPH+PSTRTV(JL,JK-1)

!   Integrate SO-tendencies to find SO-stress profile
    PSTRSOU(JL,JK) = (PSOTEU(JL,JK)+ZTESOU)*ZGDPH+PSTRSOU(JL,JK-1)
    PSTRSOV(JL,JK) = (PSOTEV(JL,JK)+ZTESOV)*ZGDPH+PSTRSOV(JL,JK-1)

    
    PVOM(JL,JK)    = ZDUDT
    PVOL(JL,JK)    = ZDVDT
    PVDIS(JL)      = PVDIS(JL) +ZVDFDIS*ZDP
    PVDISG(JL)     = PVDISG(JL)+ZSODIS*ZDP
   
    PQTDIF(JL,JK)  =   PQTDIF(JL,JK) + ZTPFAC3 * PQTM1(JL,JK)
    PQTE(JL,JK)    = ( PQTDIF(JL,JK) - PQTM1(JL,JK) ) * ZCONS1

!----------------------------------------------------
!  MPBL - Liquid water static energy
!    Input:  PSLGDIF liquid static energy first guess
!            PQTDIF  total water first guess
!    Output: PSLGE   liquid static energy tendency
!----------------------------------------------------

    PSLGDIF(JL,JK)   =   PSLGDIF(JL,JK) + ZTPFAC3 * PSLGM1(JL,JK)
    PSLGE(JL,JK)     = ( PSLGDIF(JL,JK) + ZVDFDIS+ZSODIS-PSLGM1(JL,JK) ) * ZCONS1
    PSLGEWODIS(JL,JK)= ( PSLGDIF(JL,JK)                - PSLGM1(JL,JK) ) * ZCONS1

  ENDDO
ENDDO


DO JL=KIDIA,KFDIA
  PVDIS(JL) =PVDIS(JL) *ZCONS2
  PVDISG(JL)=PVDISG(JL)*ZCONS2
ENDDO


IF (LHOOK) CALL DR_HOOK('VDFINCR',1,ZHOOK_HANDLE)
END SUBROUTINE VDFINCR


END MODULE mo_vdfincr
