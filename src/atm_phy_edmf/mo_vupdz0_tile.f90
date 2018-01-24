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

MODULE mo_vupdz0_tile
 
  PUBLIC :: vupdz0_tile

CONTAINS

SUBROUTINE VUPDZ0_TILE(KIDIA,KFDIA,KLON,JTILE,KSTEP,&
 & PRVDIFTS,&
 & PUMLEV,PVMLEV,&
 & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,&
 & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
 & PHLICE, &
 & PTSKTI,PUCURR,PVCURR,&
 & PZ0MTI,PZ0HTI,PZ0QTI,PBUOMTI,PZDLTI)  

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
      & RG       ,RD       ,RCPD     ,RETV               ,& !yomcst  (& yos_cst)   
      & RCHAR                                               !"EQUIVALENT" CHARNOCK PARAMETER
USE mo_edmf_param   ,ONLY : &
      & RVZ0M    ,RVZ0H    ,&                               !yos_veg
      & LEFLAKE  ,RH_ICE_MIN_FLK     ,&                     !yoephy  (& yos_flake)
      & LSCMEC   ,LROUGH   ,REXTZ0M  ,REXTZ0H  ,&           !yomct0  (& yos_exc)
      & PSIHU    ,PSIHS                                     !fcsvdfs.h
USE mo_lnd_nwp_config,ONLY: isub_water, isub_lake, isub_seaice, llake


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

!    Reals (In):
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
!     *PUCURR*       OCEAN CURRENT U-COMPONENT
!     *PVCURR*       OCEAN CURRENT V-COMPONENT

!    Reals (Out):
!     *PZ0MTI*       NEW AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0HTI*       NEW ROUGHNESS LENGTH FOR HEAT
!     *PZ0QTI*       NEW ROUGHNESS LENGTH FOR MOISTURE
!     *PBUOMTI*      BUOYANCY FLUX
!     *PZDLTI*       Z/L AT LOWEST MODEL LEVEL

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
INTEGER(KIND=JPIM),INTENT(IN)    :: JTILE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKTI(:) 
REAL(KIND=JPRB)   ,INTENT(IN)   ,OPTIONAL :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)   ,OPTIONAL :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0MTI(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0HTI(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0QTI(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBUOMTI(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZDLTI(:) 

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: Z1DZ0Q, ZCON2, ZIPBL, ZNLEV, ZPRH1,&
 & ZPRQ0, ZROWQ, ZROWT, ZTPFAC2, &
 & ZTPFAC3, ZTPFAC4, ZWST2, &
 & ZXLNQ,ZDUA,ZZCDN,ZCDFC  
REAL(KIND=JPRB) :: ZUST(KLON),ZUST2(KLON)
REAL(KIND=JPRB) :: ZDU2(KLON),ZRHO(KLON)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLCURR,LLINIT

!             INCLUDE STABILITY FUNCTIONS
!             ------- --------- ---------

! #include "fcsvdfs.h" ! replaced by use statements

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VUPDZ0_MOD:VUPDZ0',0,ZHOOK_HANDLE)

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

LLINIT= ( KSTEP == 0)


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
    IF ( JTILE == isub_water ) THEN
      PZ0MTI(JL)=1.E-4_JPRB                    ! Charnock?????
    ENDIF
!       - Sea ice
    IF ( JTILE == isub_seaice ) THEN
      PZ0MTI(JL)=RZ0ICE
    ENDIF
!       - LAKES   
    IF ( ( JTILE == isub_lake ) .AND. ( PHLICE(JL) > RH_ICE_MIN_FLK ) ) THEN
        PZ0MTI(JL)=RZ0ICE
!    ELSE
!        PZ0MTI(JL)=1.E-4_JPRB         ! ???????????????
    ENDIF
  ENDDO

  DO JL=KIDIA,KFDIA
    ZDUA=SQRT(ZDU2(JL))
    ZZCDN=(RKAP/LOG(1.0_JPRB+PGEOMLEV(JL)/(RG*PZ0MTI(JL))))**2
    PUSTRTI(JL)=ZRHO(JL)*PUMLEV(JL)*ZDUA*ZZCDN
    PVSTRTI(JL)=ZRHO(JL)*PVMLEV(JL)*ZDUA*ZZCDN
    PAHFSTI(JL)=0.0_JPRB
    PEVAPTI(JL)=0.0_JPRB
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

DO JL=KIDIA,KFDIA
  ZROWQ=PEVAPTI(JL)
  ZROWT=PAHFSTI(JL)/RCPD
  PBUOMTI(JL)=RG*(-RETV*ZROWQ-ZROWT/PTSKTI(JL))/ZRHO(JL)
  ZUST2(JL)=SQRT(PUSTRTI(JL)**2+PVSTRTI(JL)**2)/ZRHO(JL)

!     APPLY W* CORRECTION

  IF (PBUOMTI(JL)  >  0.0_JPRB) THEN
    ZWST2=(PBUOMTI(JL)*ZIPBL)**ZCON2
    ZUST2(JL)=ZUST2(JL)+ZCDFC*ZWST2
  ENDIF

  ZUST(JL)=MAX(SQRT(ZUST2(JL)),REPUST)
  PZDLTI(JL)=-PGEOMLEV(JL)*RKAP*PBUOMTI(JL)/(RG*ZUST(JL)**3)
ENDDO

!*         5.    SETTING OF ROUGHNESS LENGTHS 
!                ----------------------------

IF ( JTILE == isub_water ) THEN             !   - Ocean open water

  DO JL=KIDIA,KFDIA
    PZ0MTI(JL)=RNUM/ZUST(JL) + (RCHAR/RG)*ZUST2(JL)    ! could be updated to tuned ICON value
    PZ0HTI(JL)=RNUH/ZUST(JL)
    PZ0QTI(JL)=RNUQ/ZUST(JL)
  ENDDO

ELSE IF ( JTILE == isub_seaice ) THEN       !   - Sea ice

  DO JL=KIDIA,KFDIA
    PZ0MTI(JL)=RZ0ICE
    PZ0HTI(JL)=RZ0ICE
    PZ0QTI(JL)=RZ0ICE
  ENDDO

ELSE IF ( JTILE == isub_lake ) THEN         !   - Lake

  DO JL=KIDIA,KFDIA
    IF (PHLICE(JL) > RH_ICE_MIN_FLK) THEN
      PZ0MTI(JL)=RZ0ICE
      PZ0HTI(JL)=RZ0ICE
      PZ0QTI(JL)=RZ0ICE
    ELSE
      PZ0MTI(JL)=RNUM/ZUST(JL) + (RCHAR/RG)*ZUST2(JL)
      PZ0HTI(JL)=RNUH/ZUST(JL)
      PZ0QTI(JL)=RNUQ/ZUST(JL)
    ENDIF
  ENDDO

ELSE                                        !   - Land   PZ0M is set in turbtrans interface
                                            !            PZ0H and PZ0Q needs to be calculated here
  DO JL=KIDIA,KFDIA
    PZ0HTI(JL) = PZ0MTI(JL) / 1.0_JPRB      ! approximately as in table (~10, mo_edmf_param.f90) ????
    PZ0QTI(JL) = PZ0MTI(JL) / 1.0_JPRB
  ENDDO

ENDIF


!*        6.   SCM: Fixed roughness lengths

IF (LSCMEC .AND. LROUGH) THEN
  PZ0MTI(:) = REXTZ0M   ! scm namelist parameters
  PZ0HTI(:) = REXTZ0H
ENDIF


IF (LHOOK) CALL DR_HOOK('VUPDZ0_MOD:VUPDZ0',1,ZHOOK_HANDLE)
END SUBROUTINE VUPDZ0_TILE


END MODULE MO_VUPDZ0_TILE
