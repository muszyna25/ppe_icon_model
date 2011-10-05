!>
!! Turbulent orographic form drag for EDMF DUALM
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
!!
!!-----------------------------------------------------------------------------
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
!!-----------------------------------------------------------------------------

MODULE mo_vdftofdc
 
  PUBLIC :: vdftofdc

CONTAINS

SUBROUTINE VDFTOFDC(KIDIA,KFDIA,KLON,KLEV,PTMST,&
 & PUM1,PVM1,PGEOM1,PSIGFLT,&
 & PTOFDC)  
!     ------------------------------------------------------------------

!**   *VDFTOFDC* - DETERMINES THE COEFFICIENTS FOR THE 
!                 TURBULENT OROGRAPHIC DRAG PARAMETRIZATION

!     Original  A. BELJAARS   ECMWF    17/11/2002.
!     Modified 

!     PURPOSE
!     -------

!     DETERMINE COEFFICIENTS FOR TURBULENT OROGRAPHIC DRAG

!     INTERFACE
!     ---------

!     *VDFTOFDC* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEV*         NUMBER OF LEVELS

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PSIGFLT*      FILTERED STANDARD DEVIATION OF SUBGRID OROGRAPHY

!     OUTPUT PARAMETERS (REAL):

!     *PTOFDC*        COEFFICIENTS IN DIAGONAL TO BE PASSED ON 
!                    TO IMPLICIT SOLVER. PTOFDC=alpha*DT*stressdiv/PSI

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RG
! USE YOEVDF   , ONLY : RVDIFTS, REPDU2

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                 & RG      ,&                     !yomcst
                 & RVDIFTS ,REPDU2                !yoevdf

IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDC(KLON,KLEV) 
!*            LOCAL STORAGE
!             ----- -------

REAL(KIND=JPRB) ::    ZCOEF(KLON)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) ::    ZZ,ZTOFDALPHA,ZTOFDBETA,ZTOFDCMD,ZTOFDCORR,ZTOFDK1,ZTOFDIH,&
 & ZTOFDKFLT,ZTOFDN1,ZTOFDN2,ZTOFDIC,ZTOFDIZ,ZTOFDIN,ZFACT1,&
 & ZFACT2,ZUABS,ZMAX  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VDFTOFDC',0,ZHOOK_HANDLE)
ZTOFDALPHA=27.
ZTOFDBETA=1.
ZTOFDCMD=0.005
ZTOFDCORR=0.6
ZTOFDK1=0.003
ZTOFDIH=0.00102
ZTOFDKFLT=0.00035
ZTOFDN1=-1.9
ZTOFDN2=-2.8

ZTOFDIC=2.109
ZTOFDIZ=1500.
ZTOFDIN=-1.2

ZFACT1=ZTOFDK1**(ZTOFDN1-ZTOFDN2)/(ZTOFDIH*ZTOFDKFLT**ZTOFDN1)
ZFACT2=ZTOFDALPHA*ZTOFDBETA*ZTOFDCMD*ZTOFDCORR*ZTOFDIC*ZFACT1&
 & *RVDIFTS*PTMST  

ZMAX=5000.

!        1. PREPARE ARRAY INDEPENDENT OF HEIGHT
!           ------- ----- ----------- -- ------
DO JL=KIDIA,KFDIA
  ZCOEF(JL)=ZFACT2*PSIGFLT(JL)**2
ENDDO

!        2. VERTICAL LOOP
!           -------- ----

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZZ=PGEOM1(JL,JK)/RG
    IF (ZZ > ZMAX) THEN
      PTOFDC(JL,JK)=0.
    ELSE
      ZUABS=SQRT(MAX(REPDU2,PUM1(JL,JK)**2+PVM1(JL,JK)**2))
      PTOFDC(JL,JK)=ZCOEF(JL)*ZUABS*EXP(-(ZZ/ZTOFDIZ)**1.5)*ZZ**ZTOFDIN
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('VDFTOFDC',1,ZHOOK_HANDLE)
END SUBROUTINE VDFTOFDC


END MODULE mo_vdftofdc
