!>
!! Wind at 10m for EDMF DUALM
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

MODULE mo_vdffblend
 
  PUBLIC :: vdffblend

CONTAINS

SUBROUTINE VDFFBLEND(KIDIA,KFDIA,KLON,KLEV,&
 & PUM1  ,PVM1  ,PGEOM1, PUCURR, PVCURR, PBLEND,&
 ! OUTPUTS
 & PFBLEND )

!     ------------------------------------------------------------------

!**   *VDFFBLEND* - COMPUTES THE WIND SPEED AT THE BLENDING HEIGHT

!     P. VITERBO         E.C.M.W.F.    10/06/2005. (BASED ON VDFPPCFL)
!     T. Stockdale       E.C.M.W.F.    10/12/2005  Include surface currents

!     PURPOSE
!     -------

!     COMPUTE WIND SPEED AT BLENDING HEIGHT TO BE USE LATER BY CFL INTERPOLATION
!      METHOD

!     INTERFACE
!     ---------

!     *VDFFBLEND* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUM1*         U-COMPONENT WIND AT T-1
!     *PVM1*         V-COMPONENT WIND AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PUCURR*       U-COMPONENT OF OCEAN SFC CURRENT
!     *PVCURR*       V-COMPONENT OF OCEAN SFC CURRENT
!     *PBLEND*       HEIGHT FROM WHICH WIND SPEED IS INTERPOLATED TO 10 M

!     OUTPUT PARAMETERS (REAL):

!     *PFBLEND*      WIND SPEED AT BLENDING HEIGHT

!     METHOD
!     ------

!     LINEAR INTERPOLATION IN PRESSURE FROM THE MODEL LEVEL WIND SPEED TO THE
!      BLENDING HEIGHT

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RG
! USE YOMLUN   , ONLY : NULERR

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RG                              !yomcst
USE mo_edmf_param   ,ONLY : &
                & NULERR                          !yomlun
USE mo_mpi          ,ONLY : abort_mpi

IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLEND(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFBLEND(KLON) 

INTEGER(KIND=JPIM) :: IWIND(KLON)
REAL(KIND=JPRB) ::    ZBLENDR(KLON)

INTEGER(KIND=JPIM) :: IDONE, ITOP, JK, JL

REAL(KIND=JPRB) :: ZUU, ZVV
 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! #include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFFBLEND',0,ZHOOK_HANDLE)

IWIND(:)=0
ITOP=KLEV-4
IDONE=0

DO JL=KIDIA,KFDIA
  IF (PBLEND(JL)*RG  <  PGEOM1(JL,KLEV)) THEN
    ZUU=PUM1(JL,KLEV)-PUCURR(JL)
    ZVV=PVM1(JL,KLEV)-PVCURR(JL)
    PFBLEND(JL)=SQRT(ZUU**2+ZVV**2)
    ZBLENDR(JL)=PGEOM1(JL,KLEV)
    IWIND(JL)=KLEV
    IDONE=IDONE+1
  ELSE
    ZBLENDR(JL)=PBLEND(JL)*RG
  ENDIF
ENDDO

DO JK=KLEV,ITOP,-1
  IF (IDONE  <  KFDIA-KIDIA+1) THEN
    DO JL=KIDIA,KFDIA
      IF (IWIND(JL)  ==  0) THEN
        IF (ZBLENDR(JL)  <  PGEOM1(JL,JK-1) .AND.&
           & ZBLENDR(JL)  >=  PGEOM1(JL,JK)) THEN  
          ZUU=( PUM1(JL,JK-1)*(ZBLENDR(JL)-PGEOM1(JL,JK))&
           & +PUM1(JL,JK)*(PGEOM1(JL,JK-1)-ZBLENDR(JL))&
           & )/(PGEOM1(JL,JK-1)-PGEOM1(JL,JK))  
          ZVV=( PVM1(JL,JK-1)*(ZBLENDR(JL)-PGEOM1(JL,JK))&
           & +PVM1(JL,JK)*(PGEOM1(JL,JK-1)-ZBLENDR(JL))&
           & )/(PGEOM1(JL,JK-1)-PGEOM1(JL,JK))  
          PFBLEND(JL)=SQRT((ZUU-PUCURR(JL))**2+(ZVV-PVCURR(JL))**2)
          IWIND(JL)=JK-1
          IDONE=IDONE+1
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDDO

IF (IDONE  <  KFDIA-KIDIA+1) THEN
  WRITE(NULERR,*) 'ERROR in VDFFBLEND; ITOP too large; '
  WRITE(NULERR,*) 'IDONE/KFDIA-KIDIA+1/ITOP: ', IDONE,KFDIA-KIDIA+1,ITOP
!xmk CALL ABOR1('ERROR in VDFFBLEND; ITOP too large')
  CALL ABORT_MPI
ENDIF

IF (LHOOK) CALL DR_HOOK('VDFFBLEND',1,ZHOOK_HANDLE)
END SUBROUTINE VDFFBLEND


END MODULE mo_vdffblend
