!>
!! External interface ot computation of quantities at the end of
!! vertical diffusion, includting routines to post-process weather elements
!! and gustiness.
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY31R1)
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

MODULE mo_surfpp
 
  PUBLIC :: surfpp

CONTAINS

SUBROUTINE SURFPP( KIDIA,KFDIA,KLON,KTILES, KDHVTLS, KDHFTLS &
 & , PTSTEP &
! input
 & , PFRTI, PAHFLTI, PG0TI, PSTRTULEV, PSTRTVLEV, PTSKM1M &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTSPP ,PCPTGZLEV &
 & , PAPHMS, PZ0MW, PZ0HW, PZ0QW, PZDL, PQSAPP, PBLEND, PFBLEND, PBUOM &
 & , PZ0M, PEVAPSNW, PSSRFLTI, PSLRFL, PSST &
 & , PUCURR, PVCURR &
! updated
 & , PAHFSTI, PEVAPTI, PTSKE1, PTSKTIP1 &
! output
 & , PDIFTSLEV, PDIFTQLEV, PUSTRTI, PVSTRTI, PTSKTI, PAHFLEV, PAHFLSB, PFWSB  &
 & , PU10M, PV10M, PT2M, PD2M, PQ2M &
 & , PGUST &
 & )

! USE PARKIND1  ,ONLY : JPIM, JPRB
!ifndef INTERFACE
! USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
! USE ABORT_SURF_MOD
! USE SURFPP_CTL_MOD
!endif INTERFACE

!ICON definitions:
USE mo_kind             ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters     ,ONLY : lhook    ,dr_hook           !yomcst  (& yos_exc)

USE mo_edmf_param       ,ONLY : abort_surf
USE mo_surfpp_ctl       ,ONLY : surfpp_ctl

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFPP controls the computation of quantities at the end of
!     vertical diffusion, includting routines to post-process weather elements
!     and gustiness.

!  SURFPP is called by VDFMAIN

!  METHOD:
!    This routine is a shell needed by the surface library  externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005

!  REVISION HISTORY:

!  INTERFACE: 

!  INTERFACE: 

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KTILES   :    Number of files
!      KDHVTLS  :    Number of variables for individual tiles
!      KDHFTLS  :    Number of fluxes for individual tiles


!    Reals (In):
!      PTSTEP    :  Timestep                                          s
!      PFRTI     :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PAHFLTI   :  Surface latent heat flux                         Wm-2
!      PG0TI     :  Surface ground heat flux                         W/m2
!      PSTRTULEV :  TURBULENT FLUX OF U-MOMEMTUM                     kg/(m*s2)
!      PSTRTVLEV :  TURBULENT FLUX OF V-MOMEMTUM                     kg/(m*s2)
!      PTSKM1M   :  Skin temperature, t                              K
!      PUMLEV    :  X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV    :  Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PQMLEV    :  SPECIFIC HUMIDITY                                kg/kg
!      PGEOMLEV  :  Geopotential, lowest atmospehric level           m2/s2
!      PCPTSPP   :  Cp*Ts for post-processing of weather parameters  J/kg
!      PCPTGZLEV :  Geopotential, lowest atmospehric level           J/kg
!      PAPHMS    :  Surface pressure                                 Pa
!      PZ0MW     :  Roughness length for momentum, WMO station       m
!      PZ0HW     :  Roughness length for heat, WMO station           m
!      PZ0QW     :  Roughness length for moisture, WMO station       m
!      PZDL      :  z/L                                              -
!      PQSAPP    :  Apparent surface humidity                        kg/kg
!      PBLEND    :  Blending weight for 10 m wind postprocessing     m
!      PFBLEND   :  Wind speed at blending weight for 10 m wind PP   m/s
!      PBUOM     :  Buoyancy flux, for post-processing of gustiness  ???? 
!      PZ0M     :    AERODYNAMIC ROUGHNESS LENGTH                    m
!      PEVAPSNW :    Evaporation from snow under forest              kgm-2s-1
!      PSSRFLTI  :  NET SOLAR RADIATION AT THE SURFACE, TILED        Wm-2
!      PSLRFL    :  NET THERMAL RADIATION AT THE SURFACE             Wm-2
!      PSST      :  Sea surface temperatute                          K
!      PUCURR    :   U-comp of ocean surface current                 m/s
!      PVCURR    :   V-comp of ocean surface current                 m/s

!    Reals (Updated):
!      PAHFSTI   :  SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI   :  SURFACE MOISTURE FLUX                            kg/m2/s
!      PTSKE1    :  SKIN TEMPERATURE TENDENCY                        K/s
!      PTSKTIP1  :  Tile skin temperature, t+1                       K

!    Reals (Out):
!      PDIFTSLEV :  TURBULENT FLUX OF HEAT                           J/(m2*s)
!      PDIFTQLEV :  TURBULENT FLUX OF SPECIFIC HUMIDITY              kg/(m2*s)
!      PUSTRTI   :  SURFACE U-STRESS                                 N/m2 
!      PVSTRTI   :  SURFACE V-STRESS                                 N/m2 
!      PTSKTI    :  SKIN TEMPERATURE                                 K
!      PAHFLEV   :  LATENT HEAT FLUX  (SNOW/ICE FREE PART)           W/m2
!      PAHFLSB   :  LATENT HEAT FLUX  (SNOW/ICE COVERED PART)        W/m2
!      PFWSB     :  EVAPORATION OF SNOW                              kg/(m**2*s)
!      PU10M     :  U-COMPONENT WIND AT 10 M                         m/s
!      PV10M     :  V-COMPONENT WIND AT 10 M                         m/s
!      PT2M      :  TEMPERATURE AT 2M                                K
!      PD2M      :  DEW POINT TEMPERATURE AT 2M                      K
!      PQ2M      :  SPECIFIC HUMIDITY AT 2M                          kg/kg
!      PGUST     :  GUST AT 10 M                                     m/s
!!!    PDHTLS    :  Diagnostic array for tiles (see module yomcdh)
!!!                    (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)

!     EXTERNALS.
!     ----------

!     ** SURFPP_CTL CALLS SUCCESSIVELY:
!         *SPPCFL*
!         *SPPGUST*
!         *VOSKIN*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFLTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PG0TI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTULEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTVLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTSPP(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0HW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0QW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZDL(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAPP(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLEND(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFBLEND(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOM(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPSNW(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFLTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTIP1(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTSLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVSTRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSKTI(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLSB(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWSB(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(:) 

!ifndef INTERFACE

! Local variables

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFPP',0,ZHOOK_HANDLE)

IF(UBOUND(PFRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PFRTI TOO SHORT!')
ENDIF

IF(UBOUND(PFRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PFRTI TOO SHORT!')
ENDIF

IF(UBOUND(PAHFLTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PAHFLTI TOO SHORT!')
ENDIF

IF(UBOUND(PAHFLTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PAHFLTI TOO SHORT!')
ENDIF

IF(UBOUND(PG0TI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PG0TI TOO SHORT!')
ENDIF

IF(UBOUND(PG0TI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PG0TI TOO SHORT!')
ENDIF

IF(UBOUND(PSTRTULEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PSTRTULEV TOO SHORT!')
ENDIF

IF(UBOUND(PSTRTVLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PSTRTVLEV TOO SHORT!')
ENDIF

IF(UBOUND(PTSKM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PTSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PUMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PUMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PVMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PVMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PQMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PQMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PGEOMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PGEOMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PCPTSPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PCPTSPP TOO SHORT!')
ENDIF

IF(UBOUND(PCPTGZLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PCPTGZLEV TOO SHORT!')
ENDIF

IF(UBOUND(PAPHMS,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PAPHMS TOO SHORT!')
ENDIF

IF(UBOUND(PZ0MW,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PZ0MW TOO SHORT!')
ENDIF

IF(UBOUND(PZ0HW,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PZ0HW TOO SHORT!')
ENDIF

IF(UBOUND(PZ0QW,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PZ0QW TOO SHORT!')
ENDIF

IF(UBOUND(PZDL,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:  PZDL TOO SHORT!')
ENDIF

IF(UBOUND(PQSAPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PQSAPP TOO SHORT!')
ENDIF

IF(UBOUND(PBLEND,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PBLEND TOO SHORT!')
ENDIF

IF(UBOUND(PFBLEND,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PFBLEND TOO SHORT!')
ENDIF

IF(UBOUND(PBUOM,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PBUOM TOO SHORT!')
ENDIF

IF(UBOUND(PZ0M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PZ0M TOO SHORT!')
ENDIF

IF(UBOUND(PEVAPSNW,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PEVAPSNW TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFLTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PSSRFLTI TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFLTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PSSRFLTI TOO SHORT!')
ENDIF

IF(UBOUND(PSLRFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: PSLRFL TOO SHORT!')
ENDIF

IF(UBOUND(PSST,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: PSST TOO SHORT!')
ENDIF

IF(UBOUND(PUCURR,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PUCURR TOO SHORT!')
ENDIF

IF(UBOUND(PVCURR,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PVCURR TOO SHORT!')
ENDIF

IF(UBOUND(PAHFSTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PAHFSTI TOO SHORT!')
ENDIF

IF(UBOUND(PAHFSTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PAHFSTI TOO SHORT!')
ENDIF

IF(UBOUND(PEVAPTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PEVAPTI TOO SHORT!')
ENDIF

IF(UBOUND(PEVAPTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PEVAPTI TOO SHORT!')
ENDIF

IF(UBOUND(PTSKE1,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PTSKE1 TOO SHORT!')
ENDIF

IF(UBOUND(PTSKTIP1,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PTSKTIP1 TOO SHORT!')
ENDIF

IF(UBOUND(PTSKTIP1,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PTSKTIP1 TOO SHORT!')
ENDIF

IF(UBOUND(PDIFTSLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PDIFTSLEV TOO SHORT!')
ENDIF

IF(UBOUND(PDIFTQLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PDIFTQLEV TOO SHORT!')
ENDIF

IF(UBOUND(PUSTRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PUSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PUSTRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PUSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PVSTRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PVSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PVSTRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PVSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PTSKTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP:: FIRST DIMENSION OF PTSKTI TOO SHORT!')
ENDIF

IF(UBOUND(PTSKTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFPP:: SECOND DIMENSION OF PTSKTI TOO SHORT!')
ENDIF

IF(UBOUND(PAHFLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PAHFLEV TOO SHORT!')
ENDIF

IF(UBOUND(PAHFLSB,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PAHFLSB TOO SHORT!')
ENDIF

IF(UBOUND(PFWSB,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PFWSB TOO SHORT!')
ENDIF

IF(UBOUND(PU10M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PU10M TOO SHORT!')
ENDIF

IF(UBOUND(PV10M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PV10M TOO SHORT!')
ENDIF

IF(UBOUND(PT2M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PT2M TOO SHORT!')
ENDIF

IF(UBOUND(PD2M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PD2M TOO SHORT!')
ENDIF

IF(UBOUND(PQ2M,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PQ2M TOO SHORT!')
ENDIF

IF(UBOUND(PGUST,1) < KLON) THEN
  CALL ABORT_SURF('SURFPP: PGUST TOO SHORT!')
ENDIF


CALL SURFPP_CTL( KIDIA,KFDIA,KLON,KTILES, KDHVTLS, KDHFTLS &
 & , PTSTEP &
! input
 & , PFRTI, PAHFLTI, PG0TI, PSTRTULEV, PSTRTVLEV, PTSKM1M &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTSPP ,PCPTGZLEV &
 & , PAPHMS, PZ0MW, PZ0HW, PZ0QW, PZDL, PQSAPP, PBLEND, PFBLEND, PBUOM &
 & , PZ0M, PEVAPSNW, PSSRFLTI, PSLRFL, PSST &
 & , PUCURR, PVCURR &
! updated
 & , PAHFSTI, PEVAPTI, PTSKE1, PTSKTIP1 &
! output
 & , PDIFTSLEV, PDIFTQLEV, PUSTRTI, PVSTRTI, PTSKTI, PAHFLEV, PAHFLSB, PFWSB  &
 & , PU10M, PV10M, PT2M, PD2M, PQ2M &
 & , PGUST &
 & )

IF (LHOOK) CALL DR_HOOK('SURFPP',1,ZHOOK_HANDLE)

!endif INTERFACE

!------------------------------------------------------------------------

END SUBROUTINE SURFPP


END MODULE MO_SURFPP
