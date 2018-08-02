!>
!! THIS SPECIFIES THE HEIGHT FOR U,V (10 M) 
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY36R1)
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

MODULE mo_sppgust
 
  PUBLIC :: sppgust

CONTAINS

SUBROUTINE SPPGUST(KIDIA, KFDIA, KLON &
 & , PZ0MM, PBUOM, PUSTAR, PU10M, PV10M  &
 ! OUTPUTS
 & , PGUST)  

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! 
! USE YOS_CST   , ONLY : RG       ,&
!  & RLVTT    ,RLSTT    ,RTT
! USE YOS_EXCS  , ONLY : RCHBA    ,RCHBB    ,RCHBD    ,RCHB23A  ,&
!  & RCHBBCD  ,RCHBCD   ,RCHETA   ,RCHETB   ,RCHBHDL  ,&
!  & RCDHALF  ,RCDHPI2
! USE YOS_EXC   , ONLY : RKAP     ,REPDU2

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RKAP                                                !yoevdf  (& yos_exc)

!     ------------------------------------------------------------------

!**   *SPPGUST* - COMPUTES THE area averaged 10 m wind and the gust

!     Author.
!     -------
!     A. Beljaars       E.C.M.W.F.    24/02/2000
!     A. Beljaars       E.C.M.W.F.    15/11/2001 Fix orography problem

!     Modifications.
!     --------------
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     P. Viterbo  ECMWF  12/05/2005 Externalize SURF (based on vdfppgust)
!     A. Beljaars ECMWF  18/02/2006 Revised gust to accomodate stochastic physics

!     PURPOSE
!     -------

!     Compute wind gusts

!     INTERFACE
!     ---------

!     *SPPGUST* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PZ0MM*        AERODYNAMIC ROUGHNESS LENGTH
!     *PBUOM*        BUOYANCY FLUX
!     *PUSTAR*       FRICTION VELOCITY
!      PU10M         U-COMPONENT WIND AT 10 M                         m/s
!      PV10M         V-COMPONENT WIND AT 10 M                         m/s

!     OUTPUT PARAMETERS (REAL):

!     *PGUST*        WIND GUST AT 10 M

!     METHOD
!     ------

!     MO scaling is used to estimate turbulence intensity               

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOM(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTAR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU10M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(:) 

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) ::  &
 & ZUGN, ZCZI, ZUSTAR,&
 & Z1D3, ZIPBL, ZIDL, ZFF10
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! #include "fcsvdfs.h" !  ... not needed???

!     ------------------------------------------------------------------

!*       1.   INITIALIZE CONSTANTS
!             ---------- ----------

!     THIS SPECIFIES THE HEIGHT FOR U,V (10 M) 

IF (LHOOK) CALL DR_HOOK('SPPGUST_MOD:SPP_GUST',0,ZHOOK_HANDLE)

!     For the gust model, a dimensionless number is 
!     used ZUG=(UMAX-U)/u*. It is computed from the gust model 
!     with Kaimal spectrum in (see Beljaars 1987, J. Atmos. Oc. Techn., 4, 
!     613-626). The number is reasonably constant for most conditions, 
!     but dependes on Zi/L in the same way as the standard deviation 
!     of horizontal wind (Panosky et al.1977). ZUG depends also on 
!     the assumed probabilty of exceedence P in the following way for 
!     zi/L large: 

!     P    !  0.10  0.25  0.50  0.75  0.90 
!     ------------------------------------   For an anemometer with 5 s averaging
!     ZUGN !  6.79  6.13  5.49  4.93  4.47

!     P    !  0.10  0.25  0.50  0.75  0.90 
!     ------------------------------------  For an anemometer with 3 s averaging
!     ZUGN !        6.70  6.05  5.49  

!     The stability function is ZUG=ZUGN*( 1+ (0.5/12.)zi/L )^(1./3.)

!     Recimputation with a integration time of 6000 s instead of 600 s
!     The choice of 10 min is for gusts with respect to 10 min averages. 
!     The choice of 100 min is for an area of 60 km at 10 m/s (more model compatible). 

!     P    !  0.10  0.25  0.50  0.75  0.90 
!     ------------------------------------  For an anemometer with 3 s averaging
!     ZUGN !        8.23  7.71  7.26 
!             These numbers apply at 20 m/s

ZUGN=7.71_JPRB
ZCZI=0.5_JPRB/12._JPRB
Z1D3=1._JPRB/3._JPRB
ZIPBL=1000._JPRB


!     ------------------------------------------------------------------

!        2.   COMPUTE HORIZONTAL WIND AND GUTST
!             ---------------------------------

DO JL=KIDIA,KFDIA

!     AREA AVERAGE OF ABSOLUTE 10 M (TO BE USED FOR GUSTS)

  ZUSTAR=PUSTAR(JL)
  ZIDL=-ZIPBL*RKAP*PBUOM(JL)/ZUSTAR**3
  ZFF10=SQRT(PU10M(JL)**2+PV10M(JL)**2)

  IF (ZIDL >= 0.) THEN
    PGUST(JL)=ZFF10+ZUSTAR*ZUGN
  ELSE
    PGUST(JL)=ZFF10+ZUSTAR*ZUGN*(1.0_JPRB-ZCZI*ZIDL)**Z1D3
  ENDIF
  
ENDDO

IF (LHOOK) CALL DR_HOOK('SPPGUST_MOD:SPPGUST',1,ZHOOK_HANDLE)
END SUBROUTINE SPPGUST


END MODULE MO_SPPGUST
