!>
!! Top for surface layer parameter calculation
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY37R2)
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

MODULE mo_surfexcdriver_ctl
 
  PUBLIC :: surfexcdriver_ctl

CONTAINS

SUBROUTINE SURFEXCDRIVER_CTL(CDCONF &
 & , KIDIA, KFDIA, KLON, KLEVS, KTILES, KSTEP &
 & , KLEVSN, KLEVI, KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS &
 & , KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, K_VMASS &
 & , PTSTEP, PRVDIFTS &
! input data, non-tiled
 & , KTVL, KTVH, PCVL, PCVH &
 & , PUMLEV, PVMLEV, PTMLEV, PQMLEV, PAPHMS, PGEOMLEV, PCPTGZLEV &
 & , PSST, PTSKM1M, PCHAR, PSSRFL, PSLRFL, PEMIS, PTICE, PTSNOW &
 & , PHLICE,PTLICE,PTLWML &   
 & , PWLMX, PUCURR, PVCURR &
! input data, soil
 & , PTSAM1M, PWSAM1M, KSOTY &
! input data, tiled
 & , PFRTI, PALBTI &
! updated data, tiled
 & , PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI &
! updated data, non-tiled
 & , PZ0M, PZ0H &
! output data, tiled
 & , PSSRFLTI, PQSTI, PDQSTI, PCPTSTI, PCFHTI, PCFQTI, PCSATTI, PCAIRTI &
! output data, non-tiled
 & , PKHLEV, PCFMLEV, PKMFL, PKHFL, PKQFL, PEVAPSNW &
 & , PZ0MW, PZ0HW, PZ0QW, PBLENDPP, PCPTSPP, PQSAPP, PBUOMPP, PZDLPP &
! output data, diagnostics
!dmk & , PDHTLS, PDHTSS, PDHTTS, PDHTIS &   <<< deleted DDH outputs and
 & )                                      ! <<< associated compute_ddh

! USE PARKIND1  ,ONLY : JPIM, JPRB
! 
! USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
! 
! USE YOS_EXC   ,ONLY : RKAP, REPDU2, LEOCWA, LEOCCO, RZ0ICE
! USE YOS_CST   ,ONLY : RG, RD, RSIGMA, RTT, RETV, RCPD
! USE YOS_VEG   ,ONLY : RVTRSR, RVZ0M
! USE YOS_SOIL  ,ONLY : RALFMAXSN
! USE YOS_FLAKE ,ONLY : LEFLAKE    , RH_ICE_MIN_FLK
! USE YOS_THF   ,ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RKAP     ,REPDU2   ,RZ0ICE   ,&                     !yoevdf  (& yos_exc)
      & RG       ,RD       ,RSIGMA   ,RTT      ,RETV     ,& !yomcst  (& yos_cst)
      & RCPD                                                !yomcst  (& yos_cst)
USE mo_edmf_param   ,ONLY : &
      & LEOCWA   ,LEOCCO   ,&                               !yoephy  (& yos_exc)
      & LEFLAKE  ,RH_ICE_MIN_FLK                            !yoephy  (& yos_flake)

USE mo_vupdz0       ,ONLY : vupdz0
USE mo_vexcs        ,ONLY : vexcs

!dmk USE VSURF_MOD
!    USE VEVAP_MOD
!    USE SURFSEB_CTL_MOD
!    USE SRFCOTWO_MOD
!xxx USE VSFLX_MOD


! #ifdef DOC
!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFEXCDRIVER controls the ensemble of routines that prepare
!    the surface exchange coefficients and associated surface quantities
!    needed for the solution of the vertical diffusion equations. 

!  SURFEXCDRIVER is called by VDFMAIN

!  METHOD:
!    This routine is only a shell needed by the surface library
!    externalisation.

!  AUTHOR:
!    P. Viterbo       ECMWF May 2005   

!  REVISION HISTORY:
!    A. Beljaars      10/12/2005  TOFD
!    A. Beljaars      17/05/2007  clean-up of roughness length initialization
!    G. Balsamo       15/11/2007  Use aggregated Z0M for drag and dominant low
!                                 for post-processing of 2m T/TD.
!    E. Dutra/G. Balsamo  01/05/2008  lake tile
!    A. Beljaars/M.Koehler 14/01/2009 Surfcae flux bugfix for stability

!  INTERFACE: 

!    Characters (In):
!      CDCONF   :    IFS Configuration

!    Integers (In):
!      KIDIA    :    Begin point in arrays
!      KFDIA    :    End point in arrays
!      KLON     :    Length of arrays
!      KLEVS    :    Number of soil layers
!      KTILES   :    Number of tiles
!      KSTEP    :    Time step index
!      KLEVSN   :    Number of snow layers (diagnostics) 
!      KLEVI    :    Number of sea ice layers (diagnostics)
!      KDHVTLS  :    Number of variables for individual tiles
!      KDHFTLS  :    Number of fluxes for individual tiles
!      KDHVTSS  :    Number of variables for snow energy budget
!      KDHFTSS  :    Number of fluxes for snow energy budget
!      KDHVTTS  :    Number of variables for soil energy budget
!      KDHFTTS  :    Number of fluxes for soil energy budget
!      KDHVTIS  :    Number of variables for sea ice energy budget
!      KDHFTIS  :    Number of fluxes for sea ice energy budget
!      K_VMASS  :    Controls the use of vector functions in the IBM scientific
!                     library. Set K_VMASS=0 to use standard functions
!      KTVL     :    Dominant low vegetation type
!      KTVH     :    Dominant high vegetation type
!      KSOTY    :    SOIL TYPE                                        (1-6)

!    Reals (In):
!      PTSTEP   :    Timestep
!      PRVDIFTS :    Semi-implicit factor for vertical diffusion discretization

!    Reals with tile index (In): 
!      PFRTI    :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PALBTI   :    Tile albedo                                      (0-1)

!    Reals independent of tiles (In):
!      PCVL     :    LOW VEGETATION COVER (CLIMATOLOGICAL)
!      PCVH     :    HIGH VEGETATION COVER (CLIMATOLOGICAL)
!      PUMLEV   :    X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV   :    Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PTMLEV   :    TEMPERATURE,   lowest atmospheric level          K
!      PQMLEV   :    SPECIFIC HUMIDITY                                kg/kg
!      PAPHMS   :    Surface pressure                                 Pa
!      PGEOMLEV :    Geopotential, lowest atmospehric level           m2/s2
!      PCPTGZLEV:    Geopotential, lowest atmospehric level           J/kg
!      PSST     :    (OPEN) SEA SURFACE TEMPERATURE                   K
!      PUSTRTI  :    X-STRESS                                         N/m2
!      PVSTRTI  :    Y-STRESS                                         N/m2
!      PTSKM1M  :    SKIN TEMPERATURE                                 K
!      PCHAR    :    "EQUIVALENT" CHARNOCK PARAMETER                  -
!      PSSRFL   :    NET SHORTWAVE RADIATION FLUX AT SURFACE          W/m2
!      PSLRFL   :    NET LONGWAVE RADIATION FLUX AT SURFACE           W/m2
!      PEMIS    :    MODEL SURFACE LONGWAVE EMISSIVITY
!      PTSAM1M  :    SURFACE TEMPERATURE                              K
!      PWSAM1M  :    SOIL MOISTURE ALL LAYERS                         m**3/m**3
!      PTICE    :    Ice temperature, top slab                        K
!      PTSNOW   :    Snow temperature                                 K
!      PHLICE   :    Lake ice thickness                               m
!      PTLICE   :    Lake ice temperature                             K
!      PTLWML   :    Lake mean water temperature                      K
!      PWLMX    :    Maximum interception layer capacity              kg/m**2
!      PUCURR   :    u component of ocean surface current             m/s
!      PVCURR   :    v component of ocean surface current             m/s

!    Reals with tile index (In/Out):
!      PUSTRTI  :    SURFACE U-STRESS                                 N/m2 
!      PVSTRTI  :    SURFACE V-STRESS                                 N/m2
!      PAHFSTI  :    SURFACE SENSIBLE HEAT FLUX                       W/m2
!      PEVAPTI  :    SURFACE MOISTURE FLUX                            KG/m2/s
!      PTSKTI   :    SKIN TEMPERATURE                                 K

!    Reals independent of tiles (In/Out):
!      PZ0M     :    AERODYNAMIC ROUGHNESS LENGTH                     m
!      PZ0H     :    ROUGHNESS LENGTH FOR HEAT                        m

!    Reals with tile index (Out):
!      PSSRFLTI :    Tiled NET SHORTWAVE RADIATION FLUX AT SURFACE    W/m2
!      PQSTI    :    Tiled SATURATION Q AT SURFACE                    kg/kg
!      PDQSTI   :    Tiled DERIVATIVE OF SATURATION Q-CURVE           kg/kg/K
!      PCPTSTI  :    Tiled DRY STATIC ENERGY AT SURFACE               J/kg
!      PCFHTI   :    Tiled EXCHANGE COEFFICIENT AT THE SURFACE        ????
!      PCFQTI   :    Tiled EXCHANGE COEFFICIENT AT THE SURFACE        ????
!      PCSATTI  :    MULTIPLICATION FACTOR FOR QS AT SURFACE          -
!                      FOR SURFACE FLUX COMPUTATION
!      PCAIRTI  :    MULTIPLICATION FACTOR FOR Q AT  LOWEST MODEL     - 
!                      LEVEL FOR SURFACE FLUX COMPUTATION

!    Reals independent of tiles (Out):
!      PKHLEV   :    SURFACE LAYER: CH*U                              m/s
!      PCFMLEV  :    PROP. TO EXCH. COEFF. FOR MOMENTUM               ????
!                     (C-STAR IN DOC.) (SURFACE LAYER ONLY)
!      PKMFL    :    Kinematic momentum flux                          ????
!      PKHFL    :    Kinematic heat flux                              ????
!      PKQFL    :    Kinematic moisture flux                          ????
!      PEVAPSNW :    Evaporation from snow under forest               kgm-2s-1
!      PZ0MW    :    Roughness length for momentum, WMO station       m
!      PZ0HW    :    Roughness length for heat, WMO station           m
!      PZ0QW    :    Roughness length for moisture, WMO station       m
!      PBLENDPP :    Blending weight for 10 m wind postprocessing     m
!      PCPTSPP  :    Cp*Ts for post-processing of weather parameters  J/kg
!      PQSAPP   :    Apparent surface humidity for post-processing    kg/kg
!                     of weather parameters
!      PBUOMPP  :    Buoyancy flux, for post-processing of gustiness  ???? 
!      PZDLPP   :    z/L for post-processing of weather parameters    -
!      PDHTLS   :    Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!      PDHTSS   :    Diagnostic array for snow T (see module yomcdh)
!                      (Wm-2 for fluxes)
!      PDHTTS   :    Diagnostic array for soil T (see module yomcdh)
!                      (Wm-2 for fluxes)
!      PDHTIS   :    Diagnostic array for ice T (see module yomcdh)
!                      (Wm-2 for fluxes)


!     EXTERNALS.
!     ----------

!     ** SURFEXCDRIVER_CTL CALLS SUCCESSIVELY:
!         *VUPDZ0*
!         *VSURF*
!         *VEXCS*
!         *VEVAP*
!         *SURFSEB_CTL*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------
! #endif

IMPLICIT NONE

! Declaration of arguments

CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_VMASS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS

INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:)   
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(:)   
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSRFLTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDQSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFHTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFQTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCSATTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAIRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFMLEV(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKMFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKQFL(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0MW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0HW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0QW(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLENDPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCPTSPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSAPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBUOMPP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZDLPP(:)
!dmk 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(:,:,:) 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(:,:,:) 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(:,:,:) 
!xxx REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(:,:,:) 

! Local variables

INTEGER(KIND=JPIM) :: IFRMAX(KLON),IFRLMAX(KLON)

REAL(KIND=JPRB) :: ZZ0MTI(KLON,KTILES) , ZZ0HTI(KLON,KTILES) ,&
                 & ZZ0QTI(KLON,KTILES) , ZBUOMTI(KLON,KTILES),&
                 & ZZDLTI(KLON,KTILES) , ZRAQTI(KLON,KTILES) ,&
                 & ZQSATI(KLON,KTILES) , ZCFMTI(KLON,KTILES) ,&
                 & ZKMFLTI(KLON,KTILES), ZKHFLTI(KLON,KTILES),&
                 & ZKQFLTI(KLON,KTILES), ZZQSATI(KLON,KTILES),&
                 & ZJS(KLON,KTILES)    , ZJQ(KLON,KTILES)    ,&
                 & ZSSK(KLON,KTILES)   , ZTSK(KLON,KTILES)   ,&
                 & ZSSH(KLON,KTILES)   , ZSLH(KLON,KTILES)   ,&
                 & ZSTR(KLON,KTILES)   , ZG0(KLON,KTILES)    ,&
                 & ZRHOCHU(KLON,KTILES), ZRHOCQU(KLON,KTILES),&
                 & ZTSRF(KLON,KTILES)


REAL(KIND=JPRB) :: ZFRMAX(KLON)   , ZFRLMAX(KLON)  , ZALB(KLON)     , &
                 & ZSRFD(KLON)    , ZWETL(KLON)    , ZWETH(KLON)    , &
                 & ZWETHS(KLON)   , ZWETB(KLON)    , ZKHLEV(KLON)   , &
                 & ZTSA(KLON)     , ZCSNW(KLON)    , ZSSRFL1(KLON)  , &
                 & ZCBLENDM(KLON) , ZCBLENDH(KLON) , ZSL(KLON)      , &
                 & ZQL(KLON)      , ZASL(KLON)     , ZBSL(KLON)     , &
                 & ZAQL(KLON)     , ZBQL(KLON)     , ZRHO(KLON)


INTEGER(KIND=JPIM) :: JL, JTILE, IITT
LOGICAL :: LLINIT

REAL(KIND=JPRB) :: ZDUA, ZZCDN, ZQSSN, ZCOR, ZRG, ZRTMST , &
                 & ZZ0MWMO, ZBLENDWMO, ZBLENDZ0, ZCOEF1, ZCONS1
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL         :: LLAND, LLSICE, LLHISSR(KLON)

! #include "fcsttre.h"

!*         1.     Set up of general quantities
!                 ----------------------------

IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER_CTL_MOD:SURFEXCDRIVER_CTL',0,ZHOOK_HANDLE)

ZRTMST      = 1.0_JPRB/PTSTEP    ! optimization
ZRG         = 1.0_JPRB/RG        !     -"-
DO JL=KIDIA,KFDIA
  ZRHO(JL)=PAPHMS(JL)/( RD*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )
ENDDO

!*         1.1  ESTIMATE SURF.FL. FOR STEP 0
!*              (ASSUME NEUTRAL STRATIFICATION)

IF ( KSTEP == 0 .AND. CDCONF /= 'T' ) THEN
  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      PTSKTI(JL,JTILE)=PTSKM1M(JL)
    ENDDO
  ENDDO
  IF ((.NOT. LEOCWA) .AND. (.NOT. LEOCCO)) THEN
    DO JL=KIDIA,KFDIA
      PTSKTI(JL,1)=PSST(JL)
    ENDDO
  ENDIF
ENDIF

!*         1.2  UPDATE Z0

CALL VUPDZ0(KIDIA,KFDIA,KLON,KTILES,KSTEP,CDCONF,&
   & PRVDIFTS,&
   & KTVL,KTVH,PCVL,PCVH,PUMLEV,PVMLEV,&
   & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,&
   & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
   & PHLICE,& 
   & PTSKTI,PCHAR,PUCURR,PVCURR,&
   & ZZ0MTI,ZZ0HTI,ZZ0QTI,ZBUOMTI,ZZDLTI,ZRAQTI)  


!*         1.3  FIND DOMINANT SURFACE TYPE and DOMINANT LOW
!*              parameters for postprocessing

ZFRMAX(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)
ZFRLMAX(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)
IFRMAX(KIDIA:KFDIA)=1
IFRLMAX(KIDIA:KFDIA)=1
DO JTILE=2,KTILES
  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,JTILE)  >  ZFRMAX(JL)) THEN
      ZFRMAX(JL)=PFRTI(JL,JTILE)
      IFRMAX(JL)=JTILE
    ENDIF
    IF (PFRTI(JL,JTILE)  >  ZFRLMAX(JL) .AND. &
      JTILE.NE.6 .AND. JTILE.NE.7) THEN
      ZFRLMAX(JL)=PFRTI(JL,JTILE)
      IFRLMAX(JL)=JTILE
!* for tile wet-skin(s) attribute low-vegetation(4)
      IF (JTILE.EQ.3) IFRLMAX(JL)=JTILE+1
    ENDIF
  ENDDO
ENDDO

!*         Use tile average (log) Z0 for M and H
ZBLENDZ0=10._JPRB
ZCBLENDM(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)&
           &/(LOG(ZBLENDZ0/ZZ0MTI(KIDIA:KFDIA,1)))**2
ZCBLENDH(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,1)&
           &/(LOG(ZBLENDZ0/ZZ0HTI(KIDIA:KFDIA,1)))**2
DO JTILE=2,KTILES
  DO JL=KIDIA,KFDIA
    ZCBLENDM(JL)=ZCBLENDM(JL)&
           &+PFRTI(JL,JTILE)/(LOG(ZBLENDZ0/ZZ0MTI(JL,JTILE)))**2
    ZCBLENDH(JL)=ZCBLENDH(JL)&
           &+PFRTI(JL,JTILE)/(LOG(ZBLENDZ0/ZZ0HTI(JL,JTILE)))**2
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PZ0M(JL)=ZBLENDZ0*EXP(-1._JPRB/SQRT(ZCBLENDM(JL)))
  PZ0H(JL)=ZBLENDZ0*EXP(-1._JPRB/SQRT(ZCBLENDH(JL)))
ENDDO

!     ------------------------------------------------------------------

!*         2.     SURFACE BOUNDARY CONDITIONS FOR T AND Q
!                 ---------------------------------------

!    2.1 Albedo

ZALB(KIDIA:KFDIA)=&
 &  PFRTI(KIDIA:KFDIA,1)*PALBTI(KIDIA:KFDIA,1)&
 & +PFRTI(KIDIA:KFDIA,2)*PALBTI(KIDIA:KFDIA,2)&
 & +PFRTI(KIDIA:KFDIA,3)*PALBTI(KIDIA:KFDIA,3)&
 & +PFRTI(KIDIA:KFDIA,4)*PALBTI(KIDIA:KFDIA,4)&
 & +PFRTI(KIDIA:KFDIA,5)*PALBTI(KIDIA:KFDIA,5)&
 & +PFRTI(KIDIA:KFDIA,6)*PALBTI(KIDIA:KFDIA,6)&
 & +PFRTI(KIDIA:KFDIA,7)*PALBTI(KIDIA:KFDIA,7)&
 & +PFRTI(KIDIA:KFDIA,8)*PALBTI(KIDIA:KFDIA,8) 

IF (LEFLAKE) THEN
  ZALB(KIDIA:KFDIA)=ZALB(KIDIA:KFDIA)&
 & +PFRTI(KIDIA:KFDIA,9)*PALBTI(KIDIA:KFDIA,9)   
ENDIF

DO JL=KIDIA,KFDIA
  IF (ZALB(JL) == 0._JPRB .OR. ZALB(JL) == 1._JPRB) THEN
    WRITE(0,FMT='(''Alb=0. '',17F7.4)') ZALB(JL),(PFRTI(JL,JTILE),PALBTI(JL,JTILE),JTILE=1,KTILES)
  ENDIF
ENDDO
     

ZSSRFL1(KIDIA:KFDIA)=0._JPRB

LLHISSR(KIDIA:KFDIA)=.FALSE.
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
! Disaggregate solar flux but limit to 700 W/m2 (due to inconsistency
!  with albedo)
    PSSRFLTI(JL,JTILE)=((1.0_JPRB-PALBTI(JL,JTILE))/&
   & (1.0_JPRB-ZALB(JL)))*PSSRFL(JL)
    IF (PSSRFLTI(JL,JTILE) > 700._JPRB) THEN
      LLHISSR(JL)=.TRUE.
      PSSRFLTI(JL,JTILE)=700._JPRB
    ENDIF

! Compute averaged net solar flux after limiting to 700 W/m2
    ZSSRFL1(JL)=ZSSRFL1(JL)+PFRTI(JL,JTILE)*PSSRFLTI(JL,JTILE) 
  ENDDO
ENDDO

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    IF (LLHISSR(JL)) THEN
      PSSRFLTI(JL,JTILE)=PSSRFLTI(JL,JTILE)*PSSRFL(JL)/ZSSRFL1(JL)
    ENDIF
    ZSRFD(JL)=PSSRFLTI(JL,JTILE)/(1.0_JPRB-PALBTI(JL,JTILE))  
  ENDDO

!dmk  surface conditions - not needed ???
!  CALL VSURF(KIDIA,KFDIA,KLON,KLEVS,JTILE,&
!   & KTVL,KTVH,&
!   & PTMLEV  ,PQMLEV  ,PAPHMS,&
!   & PTSKTI(:,JTILE),PWSAM1M,PTSAM1M,KSOTY,&
!   & ZSRFD,ZRAQTI(:,JTILE),ZQSATI(:,JTILE),&
!   & PQSTI(:,JTILE)  ,PDQSTI(:,JTILE)  ,&
!   & ZWETB ,PCPTSTI(:,JTILE) ,ZWETL, ZWETH, ZWETHS )
!xxx
  
ENDDO

! DDH diagnostics

!dmk CALL COMPUTE_DDH

!*         3.     EXCHANGE COEFFICIENTS
!                 ---------------------

!*         3.1  SURFACE EXCHANGE COEFFICIENTS

LLINIT=KSTEP == 0
IF (KSTEP <= 3) THEN
  IITT=3
ELSE
  IITT=1
ENDIF
DO JTILE=1,KTILES

  CALL VEXCS(KIDIA,KFDIA,KLON,IITT,K_VMASS,LLINIT,PTSTEP,PRVDIFTS,&
   & PUMLEV,PVMLEV,PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,PCPTGZLEV,&
   & PCPTSTI(:,JTILE),ZQSATI(:,JTILE),&
   & ZZ0MTI(:,JTILE),ZZ0HTI(:,JTILE),&
   & ZZ0QTI(:,JTILE),ZZDLTI(:,JTILE),ZBUOMTI(:,JTILE),&
   & PUCURR,PVCURR,&
   & ZCFMTI(:,JTILE),PCFHTI(:,JTILE),&
   & PCFQTI(:,JTILE),ZKHLEV)

  DO JL=KIDIA,KFDIA
    IF (JTILE == IFRMAX(JL)) THEN 
      PKHLEV(JL)=ZKHLEV(JL)
    ENDIF
  ENDDO
ENDDO

!*         3.2  EQUIVALENT EVAPOTRANSPIRATION EFFICIENCY COEFFICIENT

DO JTILE=1,KTILES
  IF     (JTILE == 1) THEN
    ZTSA(KIDIA:KFDIA)=PSST(KIDIA:KFDIA)
  ELSEIF (JTILE == 2) THEN
    ZTSA(KIDIA:KFDIA)=PTICE(KIDIA:KFDIA)
  ELSEIF (JTILE == 5 .OR. JTILE == 7) THEN
    ZTSA(KIDIA:KFDIA)=PTSNOW(KIDIA:KFDIA)
  ELSEIF (JTILE == 9) THEN
    DO JL=KIDIA,KFDIA
      IF(PHLICE(JL) > RH_ICE_MIN_FLK ) THEN
        ZTSA(JL)=PTLICE(JL)
      ELSE
        ZTSA(JL)=PTLWML(JL)
      ENDIF
    ENDDO
  ELSE
    ZTSA(KIDIA:KFDIA)=PTSAM1M(KIDIA:KFDIA,1)
  ENDIF
!dmk evapotranspiration - not needed ???
!  CALL VEVAP(KIDIA,KFDIA,KLON,PTSTEP,PRVDIFTS,JTILE,&
!   & PWLMX ,PTMLEV  ,PQMLEV  ,PAPHMS,PTSKTI(:,JTILE),ZTSA,&
!   & PQSTI(:,JTILE),PCFQTI(:,JTILE),ZWETB,ZWETL,ZWETH,ZWETHS,&
!   & PCPTSTI(:,JTILE),PCSATTI(:,JTILE),PCAIRTI(:,JTILE),&
!   & ZCSNW)
!xxx
ENDDO

!          COMPUTE SNOW EVAPORATION FROM BELOW TREES i.e. TILE 7

! Note the use of qsat(Tsnow), rather than tile 7 skin. Skin T7 is a
! canopy temperature, definitely not what is desirable. Skin T5 can go
! up (and down ..) freely, not really what we want. The use of
! qsat (Tsnow) is tantamount to neglecting the skin effect there.

!dmk
!DO JL=KIDIA,KFDIA
!  ZQSSN=FOEEW(PTSNOW(JL))/PAPHMS(JL)
!  ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSSN)
!  ZQSSN=ZQSSN*ZCOR
!  PEVAPSNW(JL)=PCFQTI(JL,7)*ZCSNW(JL)*&
!   & (PQMLEV(JL)-ZQSSN)*ZRG*ZRTMST
!ENDDO
!xxx

!-------------------------------------------------------------------------


!*         3.2  CALL TERRA

! CALL nwp_surface_terra_edmf (&   
!   p_patch          =   , & !>in   !!! DELETE
!   ext_data         =   , & !>in   !!! DELETE
!   jb               =   , & !      !!! DELETE
!   nproma           = KLON   , & ! array dimensions
!   i_startidx       = KIDIA  , & ! start index for computations in the parallel program
!   i_endidx         = KFDIA  , & ! end index for computations in the parallel program
!   nsubs0           =   , & ! nsubs0=1 for single tile, nsubs0=2 for multi-tile
!   nsubs1           =   , & ! nsubs1=1 for single tile, nsubs1=#tiles+1 for multi-tile
!   dt               = PTSTEP , &
!  
!   u_ex             = PUMLEV  , & ! zonal wind speed                        ( m/s )
!   v_ex             = PVMLEV  , & ! meridional wind speed                   ( m/s )
!   t_ex             = PTMLEV  , & ! temperature                             (  k  )
!   qv_ex            = PQMLEV  , & ! specific water vapor content            (kg/kg)
!   p0_ex            =   , & ! pressure lowest level                         ( Pa  ) 
!   ps_ex            = PAPHMS  , & ! surface pressure                        ( Pa  )
!  
!   t_snow_ex        =   , & ! temperature of the snow-surface               (  K  )
!   t_snow_mult_ex   =   , & ! temperature of the snow-surface               (  K  )
!   t_s_ex           =   , & ! temperature of the ground surface             (  K  )
!   t_g_ex           =   , & ! surface temperature                           (  K  )
!   qv_s_ex          =   , & ! specific humidity at the surface              (kg/kg)
!   w_snow_ex        =   , & ! water content of snow                         (m H2O)
!   rho_snow_ex      =   , & ! snow density                                  (kg/m**3)
!   rho_snow_mult_ex =   , & ! snow density                                  (kg/m**3)
!   h_snow_ex        =   , & ! snow height                                   (  m  )
!   w_i_ex           =   , & ! water content of interception water           (m H2O)
!   t_so_ex          =   , & ! soil temperature (main level)                 (  K  )
!   w_so_ex          =   , & ! total water conent (ice + liquid water)       (m H20)
!   w_so_ice_ex      =   , & ! ice content                                   (m H20)
!   t_2m_ex          =   , & ! temperature in 2m                             (  K  )
!   u_10m_ex         =   , & ! zonal wind in 10m                             ( m/s )
!   v_10m_ex         =   , & ! meridional wind in 10m                        ( m/s )
!  
!   freshsnow_ex     =   , & ! indicator for age of snow in top of snow layer(  -  )
!   wliq_snow_ex     =   , & ! liquid water content in the snow              (m H2O)
!   wtot_snow_ex     =   , & ! total (liquid + solid) water content of snow  (m H2O)
!   dzh_snow_ex      =   , & ! layer thickness between half levels in snow   (  m  )
!   subsfrac_ex      = PFRTI  , & ! tile fraction                            (  1  )
!   
!   prr_con_ex       =   , & ! precipitation rate of rain, convective        (kg/m2*s)
!   prs_con_ex       =   , & ! precipitation rate of snow, convective        (kg/m2*s)
!   prr_gsp_ex       =   , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
!   prs_gsp_ex       =   , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
!  
!   tch_ex           =   , & ! turbulent transfer coefficient for heat       ( -- )
!   tcm_ex           =   , & ! turbulent transfer coefficient for momentum   ( -- )
!   tfv_ex           =   , & ! laminar reduction factor for evaporation      ( -- )
!  
!   sobs_ex          =   , & ! solar radiation at the ground                 ( W/m2)
!   thbs_ex          =   , & ! thermal radiation at the ground               ( W/m2)
!   pabs_ex          =   , & !!!! photosynthetic active radiation            ( W/m2)
!  
!   runoff_s_ex      =   , & ! surface water runoff; sum over forecast       (kg/m2)
!   runoff_g_ex      =   , & ! soil water runoff; sum over forecast          (kg/m2)
!  
!   t_g              =   , & ! surface temperature (grid mean)               ( K )
!   qv_s             =     & ! surface specific humidity (grid mean)         (kg/kg)
!                    )



!-------------------------------------------------------------------------


!*         3.3  COMPUTE SURFACE FLUXES FOR TILES

ZCONS1 =RG*PTSTEP*PRVDIFTS
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
!   use previous times tep fluxes for heat and moisture
    ZKHFLTI(JL,JTILE)=PAHFSTI(JL,JTILE)/(ZRHO(JL)*RCPD)
    ZKQFLTI(JL,JTILE)=PEVAPTI(JL,JTILE)/ZRHO(JL)

    ZKMFLTI(JL,JTILE)=ZCFMTI(JL,JTILE)*SQRT((PUMLEV(JL)-PUCURR(JL))**2&
   & +(PVMLEV(JL)-PVCURR(JL))**2)/(ZCONS1*ZRHO(JL))
  ENDDO
ENDDO


!*         3.3a   PREPARE ARRAY'S FOR CALL TO SURFACE ENERGY
!                 BALANCE ROUTINE

!xmk: as a first try turn of the first step SEB 
!IF (KSTEP == 0) THEN 
IF (1 == 0) THEN 
!xxx
  IF (LEOCWA .OR. LEOCCO) THEN
    ZTSRF(KIDIA:KFDIA,1)=PTSKTI(KIDIA:KFDIA,1)
  ELSE
    ZTSRF(KIDIA:KFDIA,1)=PSST(KIDIA:KFDIA)
  ENDIF
  ZTSRF(KIDIA:KFDIA,2)=PTICE(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,3)=PTSAM1M(KIDIA:KFDIA,1)
  ZTSRF(KIDIA:KFDIA,4)=PTSAM1M(KIDIA:KFDIA,1)
  ZTSRF(KIDIA:KFDIA,5)=PTSNOW(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,6)=PTSAM1M(KIDIA:KFDIA,1)
  ZTSRF(KIDIA:KFDIA,7)=PTSNOW(KIDIA:KFDIA)
  ZTSRF(KIDIA:KFDIA,8)=PTSAM1M(KIDIA:KFDIA,1)

  ZCOEF1=1.0_JPRB/(RG*PRVDIFTS*PTSTEP)
  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      ZRHOCHU(JL,JTILE)=PCFHTI(JL,JTILE)*ZCOEF1
      ZRHOCQU(JL,JTILE)=PCFQTI(JL,JTILE)*ZCOEF1
    ENDDO
  ENDDO

  ZASL(KIDIA:KFDIA)=0.0_JPRB
  ZBSL(KIDIA:KFDIA)=PCPTGZLEV(KIDIA:KFDIA)
  ZAQL(KIDIA:KFDIA)=0.0_JPRB
  ZBQL(KIDIA:KFDIA)=PQMLEV(KIDIA:KFDIA)


!*         3.3b   CALL TO SURFACE ENERGY BALANCE ROUTINE

!dmk SEB - not needed ???
!  CALL SURFSEB_CTL(KIDIA,KFDIA,KLON,KTILES,&
!   & PCPTSTI,PTSKTI,PQSTI,&
!   & PDQSTI,ZRHOCHU,ZRHOCQU,&
!   & PCAIRTI,PCSATTI,&
!   & PSSRFLTI,PFRTI,ZTSRF,&
!   & PHLICE,&
!   & PSLRFL,PTSKM1M,PEMIS,&
!   & ZASL,ZBSL,ZAQL,ZBQL,&
!   !out
!   & ZJS,ZJQ,ZSSK,ZTSK,&
!   & ZSSH,ZSLH,ZSTR,ZG0,&
!   & ZSL,ZQL)  
!xxx

  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      ZKHFLTI(JL,JTILE)=ZSSH(JL,JTILE)/(ZRHO(JL)*RCPD)
      ZKQFLTI(JL,JTILE)=ZJQ(JL,JTILE)/ZRHO(JL)
    ENDDO
  ENDDO

ENDIF


!          ADD SNOW EVAPORATION FROM BELOW TREES i.e. TILE 7

ZKQFLTI(KIDIA:KFDIA,7)=ZKQFLTI(KIDIA:KFDIA,7)+&
 & ZCSNW(KIDIA:KFDIA)*ZKQFLTI(KIDIA:KFDIA,5)  

!*         3.4  COMPUTE SURFACE FLUXES, WEIGHTED AVERAGE OVER TILES

PKMFL(KIDIA:KFDIA)=0.0_JPRB
PKHFL(KIDIA:KFDIA)=0.0_JPRB
PKQFL(KIDIA:KFDIA)=0.0_JPRB
PCFMLEV(KIDIA:KFDIA)=0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    PKMFL(JL)=PKMFL(JL)+PFRTI(JL,JTILE)*ZKMFLTI(JL,JTILE)
    PKHFL(JL)=PKHFL(JL)+PFRTI(JL,JTILE)*ZKHFLTI(JL,JTILE)
    PKQFL(JL)=PKQFL(JL)+PFRTI(JL,JTILE)*ZKQFLTI(JL,JTILE)
    PCFMLEV(JL)=PCFMLEV(JL)+PFRTI(JL,JTILE)*ZCFMTI(JL,JTILE)
  ENDDO
ENDDO

!*         4.  Preparation for "POST-PROCESSING" of surface weather parameters

!          POST-PROCESSING WITH MINIMUM OF LOCAL AND EFFECTIVE
!          SURFACE ROUGHNESS LENGTH. THE LOCAL ONES ARE FOR
!          WMO-TYPE WIND STATIONS I.E. OPEN TERRAIN WITH GRASS

ZBLENDWMO=75._JPRB
ZZ0MWMO=0.03_JPRB
DO JL=KIDIA,KFDIA
  IF (PZ0M(JL)  >  ZZ0MWMO) THEN
    PZ0MW(JL)=ZZ0MWMO
    PBLENDPP(JL)=ZBLENDWMO
  ELSE
    PZ0MW(JL)=PZ0M(JL)
    PBLENDPP(JL)=PGEOMLEV(JL)*ZRG
  ENDIF
ENDDO

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZZQSATI(JL,JTILE)=PQMLEV(JL)*(1.0_JPRB-PCAIRTI(JL,JTILE))&
     & +PCSATTI(JL,JTILE)*PQSTI(JL,JTILE)  
    ZZQSATI(JL,JTILE)=MAX(1.0E-12_JPRB,ZZQSATI(JL,JTILE))
  ENDDO
ENDDO

!          ROUGHNESS LENGTH FOR HEAT and MOISTURE ARE TAKEN
!          FROM THE DOMINANT LOW-VEG. TYPE
DO JL=KIDIA,KFDIA
  JTILE=IFRLMAX(JL)
  PZ0HW(JL)=ZZ0HTI(JL,JTILE)
  PZ0QW(JL)=ZZ0QTI(JL,JTILE)
  PCPTSPP(JL)=PCPTSTI(JL,JTILE)
  PQSAPP(JL)=ZZQSATI(JL,JTILE)
  PBUOMPP(JL)=ZBUOMTI(JL,JTILE)
  PZDLPP(JL)=ZZDLTI(JL,JTILE)
ENDDO



IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER_CTL_MOD:SURFEXCDRIVER_CTL',1,ZHOOK_HANDLE)

END SUBROUTINE SURFEXCDRIVER_CTL


!dmk SUBROUTINE COMPUTE_DDH
!...


END MODULE MO_SURFEXCDRIVER_CTL
