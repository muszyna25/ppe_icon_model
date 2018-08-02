!> 
!! Calculation of Z0 over land and sea
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY37R2)
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

MODULE mo_vupdz0
 
  PUBLIC :: vupdz0

CONTAINS

SUBROUTINE VUPDZ0(KIDIA,KFDIA,KLON,KTILES,KSTEP,CDCONF,&
 & PRVDIFTS,&
 & KTVL,KTVH,PCVL,PCVH,PUMLEV,PVMLEV,&
 & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,&
 & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
 & PHLICE, &
 & PTSKTI,PCHAR,PUCURR,PVCURR,&
 & PZ0MTI,PZ0HTI,PZ0QTI,PBUOMTI,PZDLTI,PRAQTI)  

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! 
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOS_EXC   , ONLY : RKAP     , RZ0ICE   ,REPDU2   ,&
!  & REPUST   ,RNUM     ,RNUH     ,RNUQ     ,RPARZI  
! USE YOS_EXCS  , ONLY : JPRITBL  ,RITBL    ,RCHBA    ,RCHBB    ,&
!  & RCHBD    ,RCHB23A  ,RCHBBCD  ,RCHBCD   ,RCHETA   ,&
!  & RCHETB   ,RCHBHDL  ,RCDHALF  ,RCDHPI2  ,DRITBL  
! USE YOS_CST   , ONLY : RG       ,RD       ,RCPD     ,RETV    
! USE YOS_VEG   , ONLY : RVZ0M    ,RVZ0H
! USE YOS_FLAKE , ONLY : LEFLAKE  ,RH_ICE_MIN_FLK
! USE YOS_EXC   , ONLY : LSCMEC   ,LROUGH   ,REXTZ0M  ,REXTZ0H

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RKAP     ,RZ0ICE   ,REPDU2   ,&                     !yoevdf  (& yos_exc)
      & REPUST   ,RNUM     ,RNUH     ,RNUQ     ,RPARZI   ,& ! -
      & RG       ,RD       ,RCPD     ,RETV                  !yomcst  (& yos_cst)   
USE mo_edmf_param   ,ONLY : &
      & RVZ0M    ,RVZ0H    ,&                               !yos_veg
      & LEFLAKE  ,RH_ICE_MIN_FLK     ,&                     !yoephy  (& yos_flake)
      & LSCMEC   ,LROUGH   ,REXTZ0M  ,REXTZ0H  ,&           !yomct0  (& yos_exc)
      & PSIHU    ,PSIHS                                     !fcsvdfs.h


!! #ifdef DOC
!     ------------------------------------------------------------------

!**   *VUPDZ0* - COMPUTES Z0M,Z0H,Z0Q OVER SEA; SETS Z0H,Z0Q OVER LAND

!     Original   A.C.M. BELJAARS       E.C.M.W.F.    26/03/90.
!     Modified   A.C.M. BELJAARS  26/03/99   Surface tiling
!     Modified   P. Viterbo  ECMWF 12/05/2005 Externalize SURF
!     Modified   A. Beljaars ECMWF 03/12/2005 Roughness tables + TOFD
!     Modified   A. Beljaars ECMWF 17/05/2007 Clean-up of z0 initialization
!     Modified   E. Dutra/G. Balsamo 01/05/2008 Lake tile

!     PURPOSE
!     -------

!     DERIVE Z0M,Z0H AND Z0Q FROM SURFACE FLUXES OVER SEA, SET Z0H AND
!     Z0Q OVER LAND AND DERIVE THE BUOYANCY FLUX.
!     (THE T-1 VALUES ARE UPDATED WITH FLUXES FROM THE PREVIOUS TIME
!      STEP)

!     INTERFACE
!     ---------

!     *VUPDZ0* IS CALLED BY *SURFEXCDRIVER_CTL*

!     Integer (In):
!     *KIDIA*        START OF LOOPS
!     *KFDIA*        END OF LOOPS
!     *KLON*         NUMBER OF POINTS IN PACKET
!     *KTILES*       NUMBER OF TILES
!     *KSTEP*        Time step index

!    Characters (In):
!     *CDCONF*       IFS Configuration

!    Real (In):
!      PRVDIFTS :    Semi-implicit factor for vertical diffusion discretization

!    Integer (in):
!     *KTVL*         LOW VEGETATION TYPE 
!     *KTVH*         HIGH VEGETATION TYPE 

!    Reals (In):
!     *PCVL*         LOW VEGETATION COVER (CLIMATOLOGICAL)
!     *PCVH*         HIGH VEGETATION COVER (CLIMATOLOGICAL)
!     *PUMLEV*       WIND X-COMPONENT AT T-1, lowest model level
!     *PVMLEV*       WIND Y-COMPONENT AT T-1, lowest model level
!     *PTMLEV*       TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*       SPECIFIC HUMUDITY AT T-1, lowest model level
!     *PAPHMS*       PRESSURE AT T-1, surface
!     *PGEOMLEV*     GEOPOTENTIAL T-1, lowest model level
!     *PUSTRTI*      X-STRESS
!     *PVSTRTI*      Y-STRESS
!     *PAHFSTI*      SENSIBLE HEAT FLUX
!     *PEVAPTI*      MOISTURE FLUX
!     *PHLICE*       LAKE ICE THICKNESS
!     *PTSKTI*       SURFACE TEMPERATURE
!     *PCHAR*        "EQUIVALENT" CHARNOCK PARAMETER
!     *PUCURR*       OCEAN CURRENT U-COMPONENT
!     *PVCURR*       OCEAN CURRENT V-COMPONENT

!    Reals (Out):
!     *PZ0MTI*       NEW AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0HTI*       NEW ROUGHNESS LENGTH FOR HEAT
!     *PZ0QTI*       NEW ROUGHNESS LENGTH FOR MOISTURE
!     *PBUOMTI*      BUOYANCY FLUX
!     *PZDLTI*       Z/L AT LOWEST MODEL LEVEL
!     *PRAQTI*       PRELIMINARY AERODYNAMIC RESISTANCE FOR MOISTURE 

!    Additional parameters for boundary condition (in SCM model):

!    *LROUGH*       If .TRUE. surface roughness length is externally specified 
!    *REXTZ0M*      Roughness length for momentum [m]
!    *REXTZ0H*      Roughness length for heat [m]

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------
!! #endif

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP

CHARACTER(LEN=1)  ,INTENT(IN)   ,OPTIONAL :: CDCONF 
 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)   ,OPTIONAL :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)   ,OPTIONAL :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0MTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0HTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0QTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBUOMTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZDLTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAQTI(:,:) 

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JL,JTILE

REAL(KIND=JPRB) :: Z1DZ0Q, ZCON2, ZIPBL, ZNLEV, ZPRH1,&
 & ZPRQ0, ZROWQ, ZROWT, ZTPFAC2, &
 & ZTPFAC3, ZTPFAC4, ZWST2, &
 & ZXLNQ,ZDUA,ZZCDN,ZCDFC  
REAL(KIND=JPRB) :: ZUST(KLON,KTILES),ZUST2(KLON,KTILES)
REAL(KIND=JPRB) :: ZDU2(KLON),ZRHO(KLON)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZLICE(KLON),ZLWAT(KLON) 

LOGICAL :: LLCURR,LLINIT

!             INCLUDE STABILITY FUNCTIONS
!             ------- --------- ---------

! #include "fcsvdfs.h" ! replaced by use statements

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VUPDZ0_MOD:VUPDZ0',0,ZHOOK_HANDLE)

IF (KTILES.LT.8) THEN
  STOP "Wrong number of tiles in VDFUPDZ0"
ENDIF

IF (LEFLAKE .OR. KTILES .GT. 8) THEN
!* FIND LAKE POINTS WITH ICE COVER 
  DO JL=KIDIA,KFDIA
    IF(PHLICE(JL) > RH_ICE_MIN_FLK ) THEN
      ZLICE(JL)=1._JPRB
      ZLWAT(JL)=0._JPRB
    ELSE
      ZLICE(JL)=0._JPRB
      ZLWAT(JL)=1._JPRB
    ENDIF
  ENDDO
ENDIF 

ZTPFAC2=1.0_JPRB/PRVDIFTS
ZTPFAC3=1.0_JPRB-ZTPFAC2
ZTPFAC4=1.0_JPRB+ZTPFAC3

ZCON2  =2.0_JPRB/3._JPRB

!     PBL HEIGHT FOR W* - EFFECT

ZIPBL=RPARZI

IF(PRESENT(PUCURR).AND.PRESENT(PVCURR)) THEN
  LLCURR=.TRUE.
ELSE
  LLCURR=.FALSE.
ENDIF

IF (PRESENT(CDCONF)) THEN
  LLINIT= ( KSTEP == 0 .AND. CDCONF /= 'T' ) 
ELSE
  LLINIT= ( KSTEP == 0)
ENDIF
!     ------------------------------------------------------------------

!*         2.      PRE-COMPUTATION OF TILE INDEPENDENT ARRAYS
!                  
DO JL=KIDIA,KFDIA
  ZRHO(JL)=PAPHMS(JL)/( RD*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )
  IF(LLCURR) THEN
    ZDU2(JL)=MAX(REPDU2,(PUMLEV(JL)-PUCURR(JL))**2+&
     & (PVMLEV(JL)-PVCURR(JL))**2)
  ELSE
    ZDU2(JL)=MAX(REPDU2,PUMLEV(JL)**2+PVMLEV(JL)**2)
  ENDIF
ENDDO
!*         3.   ESTIMATE SURF.FL. FOR STEP 0
!*              (ASSUME NEUTRAL STRATIFICATION)

IF (LLINIT) THEN
  DO JL=KIDIA,KFDIA
!       - Preliminatry value for water
    PZ0MTI(JL,1)=1.E-4_JPRB
!       - Sea ice
    PZ0MTI(JL,2)=RZ0ICE
!       - Wet skin
    PZ0MTI(JL,3)=PCVL(JL)*RVZ0M(KTVL(JL))+PCVH(JL)*RVZ0M(KTVH(JL))&
      & +(1.-PCVL(JL)-PCVH(JL))*RVZ0M(0)  
!       - Low Vegetation
    PZ0MTI(JL,4)=RVZ0M(KTVL(JL))
!       - Exposed snow
    PZ0MTI(JL,5)=RVZ0M(12)
!       - High vegetation
    PZ0MTI(JL,6)=RVZ0M(KTVH(JL))
!       - Sheltered snow
    PZ0MTI(JL,7)=RVZ0M(KTVH(JL))
!       - Bare soil
    PZ0MTI(JL,8)=RVZ0M(0)
    IF (LEFLAKE .OR. KTILES .GT. 8 ) THEN
!       - LAKES   
      PZ0MTI(JL,9)=PZ0MTI(JL,1)*ZLWAT(JL)+PZ0MTI(JL,2)*ZLICE(JL)   
    ENDIF
  ENDDO

  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      ZDUA=SQRT(ZDU2(JL))
      ZZCDN=(RKAP/LOG(1.0_JPRB+PGEOMLEV(JL)/(RG*PZ0MTI(JL,JTILE))))**2
      PUSTRTI(JL,JTILE)=ZRHO(JL)*PUMLEV(JL)*ZDUA*ZZCDN
      PVSTRTI(JL,JTILE)=ZRHO(JL)*PVMLEV(JL)*ZDUA*ZZCDN
      PAHFSTI(JL,JTILE)=0.0_JPRB
      PEVAPTI(JL,JTILE)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

!*         4.      STABILITY PARAMETERS AND FREE CONVECTION
!                  VELOCITY SCALE
!
!                  ZCDFC is a rough estimate of the drag coefficient used  
!                  to convert w*-gustiness at the lowest model level into u*. 
!                  The value is choosen from Fig. 1 on page 37 of the ECMWF 
!                  seminar proceedings on "Atmopshere-surface interaction", 
!                  ie charqacteristic for 1 m/s in unstable situations. 
ZCDFC=2.E-3_JPRB

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZROWQ=PEVAPTI(JL,JTILE)
    ZROWT=PAHFSTI(JL,JTILE)/RCPD
    PBUOMTI(JL,JTILE)=RG*(-RETV*ZROWQ-ZROWT/PTSKTI(JL,JTILE))/ZRHO(JL)
    ZUST2(JL,JTILE)=SQRT(PUSTRTI(JL,JTILE)**2+PVSTRTI(JL,JTILE)**2)/ZRHO(JL)

!     APPLY W* CORRECTION

    IF (PBUOMTI(JL,JTILE)  >  0.0_JPRB) THEN
      ZWST2=(PBUOMTI(JL,JTILE)*ZIPBL)**ZCON2
      ZUST2(JL,JTILE)=ZUST2(JL,JTILE)+ZCDFC*ZWST2
    ENDIF

    ZUST(JL,JTILE)=MAX(SQRT(ZUST2(JL,JTILE)),REPUST)
    PZDLTI(JL,JTILE)=-PGEOMLEV(JL)*RKAP*PBUOMTI(JL,JTILE)/(RG*ZUST(JL,JTILE)**3)
  ENDDO
ENDDO

!*         5.    SETTING OF ROUGHNESS LENGTHS 
!                ----------------------------

JTILE=1
!   - Ocean open water
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RNUM/ZUST(JL,JTILE) + (PCHAR(JL)/RG)*ZUST2(JL,JTILE)
  PZ0HTI(JL,JTILE)=RNUH/ZUST(JL,JTILE)
  PZ0QTI(JL,JTILE)=RNUQ/ZUST(JL,JTILE)
ENDDO

JTILE = 2
!   - Sea ice
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RZ0ICE
  PZ0HTI(JL,JTILE)=RZ0ICE
  PZ0QTI(JL,JTILE)=RZ0ICE
ENDDO

JTILE = 3
!   - Wet skin
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=PCVL(JL)*RVZ0M(KTVL(JL))+PCVH(JL)*RVZ0M(KTVH(JL))+(1.-PCVL(JL)-PCVH(JL))*&
    &RVZ0M(0)
  PZ0HTI(JL,JTILE)=PCVL(JL)*RVZ0H(KTVL(JL))+PCVH(JL)*RVZ0H(KTVH(JL))+(1.-PCVL(JL)-PCVH(JL))*&
    &RVZ0H(0)
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 4
!   - Low Vegetation
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(KTVL(JL))
  PZ0HTI(JL,JTILE)=RVZ0H(KTVL(JL))
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 5
!   - Exposed snow
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(12)
  PZ0HTI(JL,JTILE)=RVZ0H(12)
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 6
!   - High vegetation
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(KTVH(JL))
  PZ0HTI(JL,JTILE)=RVZ0H(KTVH(JL))
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 7
!   - Sheltered snow
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(KTVH(JL))
  PZ0HTI(JL,JTILE)=RVZ0H(KTVH(JL))
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 8
!   - Bare soil
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(0)
  PZ0HTI(JL,JTILE)=RVZ0H(0)
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

IF (LEFLAKE .OR. KTILES .GT. 8) THEN
  JTILE = 9 
  DO JL=KIDIA,KFDIA
    PZ0MTI(JL,JTILE)=PZ0MTI(JL,1)*ZLWAT(JL)+PZ0MTI(JL,2)*ZLICE(JL)
    PZ0HTI(JL,JTILE)=PZ0HTI(JL,1)*ZLWAT(JL)+PZ0HTI(JL,2)*ZLICE(JL)
    PZ0QTI(JL,JTILE)=PZ0QTI(JL,1)*ZLWAT(JL)+PZ0QTI(JL,2)*ZLICE(JL)
  ENDDO 
ENDIF

!*        6.   SCM: Fixed roughness lengths

IF (LSCMEC .AND. LROUGH) THEN
  PZ0MTI(:,:) = REXTZ0M   ! scm namelist parameters
  PZ0HTI(:,:) = REXTZ0H
ENDIF


!*        7.   COMPUTE PRELIMINARY AERODYNAMIC RESISTANCE FOR COMPUTATION 
!*             OF APPARENT SURFACE HUMIDITY IN VDFSURF

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZNLEV=PGEOMLEV(JL)/RG+PZ0MTI(JL,JTILE)
    ZXLNQ=LOG(ZNLEV/PZ0QTI(JL,JTILE))
    Z1DZ0Q=ZNLEV/PZ0QTI(JL,JTILE)
    IF (PZDLTI(JL,JTILE)  >  0.0_JPRB) THEN
      ZPRH1=PSIHS(PZDLTI(JL,JTILE))
      ZPRQ0=PSIHS(PZDLTI(JL,JTILE)/Z1DZ0Q)
    ELSE
      ZPRH1=PSIHU(PZDLTI(JL,JTILE))
      ZPRQ0=PSIHU(PZDLTI(JL,JTILE)/Z1DZ0Q)
    ENDIF
    PRAQTI(JL,JTILE)=(ZXLNQ-ZPRH1+ZPRQ0)/(ZUST(JL,JTILE)*RKAP)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('VUPDZ0_MOD:VUPDZ0',1,ZHOOK_HANDLE)
END SUBROUTINE VUPDZ0


END MODULE MO_VUPDZ0
