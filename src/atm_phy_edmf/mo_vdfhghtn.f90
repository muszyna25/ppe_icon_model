!>
!! Top routine for parcel model for EDMF DUALM
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

MODULE mo_vdfhghtn
 
  PUBLIC :: vdfhghtn

CONTAINS

!! !OPTIONS XOPT(HSFUN)
SUBROUTINE VDFHGHTN (KIDIA   , KFDIA   , KLON    , KLEV   , KDRAFT  , PTMST  , &
                   & PUM1    , PVM1    , PTM1    , PQM1   , PLM1    , PIM1   , &
                   & PAPHM1  , PAPM1   , PGEOM1  , PGEOH  , PVERVEL , &
                   & PKMFL   , PKHFL   , PKQFL   , PMFLX  , &
! DIAGNOSTIC OUTPUT
                   & PEXTR2  , KFLDX2  , PEXTRA  , KLEVX  , KFLDX   , JCNT   , LLDIAG   , &
!                   
                   & PUUH    , PVUH    , PSLGUH  , PQTUH  , PFRACB  , PWUH   , &
                   & PZPTOP  , KPTOP   , PZPLCL  , KPLCL  , KPLZB   , &
                   & PWUAVG  , PRICUI  , PMCU    , PDTHV  , &
                   & PFPLVL  , PFPLVN  , PDETR   , &
!amk: for convective preconditioning
                   & PVAR    , &
!xxx
!amk
                   & LDLAND  , &   
!xxx
                   & LDNODECP, KPBLTYPE, PWQT2)
                   
                     
!     ------------------------------------------------------------------

!**   *VDFHGHTN* - DETERMINES THE PBL-HEIGHT AND STRONG UPDRAFT FIELDS
!                  USING A ENTRAINING PARCEL ASCENT METHOD.

!     A.P. SIEBESMA    30/06/99    Original (dry)
!     M. Ko"hler        3/12/2004  Moist Version
!     P. Lopez         02/06/2005 Removed useless option LPHYLIN
!     Roel Neggers     12/04/2005  Multiple updraft extension


!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT AND UPDRAFT FIELDS

!     INTERFACE
!     ---------

!     *VDFHGHTN* IS CALLED BY *VDFMAIN*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KDRAFT*       NUMBER OF EXPLICITLY MODELED DRAFTS - CURRENTLY 3:
!                    1: test parcel
!                    2: rising dry thermals which stop at cloud base or inversion
!                    3: rising dry thermals which become cloudy
!                    (4: downdrafts .. to be done?)

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)        S
!     *PUM1*         X-VELOCITY COMPONENT AT T-1                  M/S
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1                  M/S
!     *PTM1*         TEMPERATURE AT T-1                           K
!     *PQM1*         SPECIFIC HUMUDITY AT T-1                     KG/KG
!     *PLM1*         SPECIFIC CLOUD LIQUID WATER AT T-1           KG/KG
!     *PIM1*         SPECIFIC CLOUD ICE AT T-1                    KG/KG
!!!   *PAM1*         CLOUD FRACTION AT T-1                        KG/KG
!     *PAPHM1*       PRESSURE AT HALF LEVEL AT T-1                PA
!     *PAPM1*        PRESSURE AT FULL LEVEL AT T-1                PA
!     *PGEOM1*       GEOPOTENTIAL AT T-1                          M2/S2
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                   M2/S2
!     *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2  
!     *PKHFL*        SURFACE KINEMATIC HEAT FLUX                  K*M/S
!     *PKQFL*        SURFACE KINEMATIC MOISTURE FLUX              M/S
!!!   *PBIR*         BUOYANCY-FLUX INTEGRAL RATIO (-N/P)
!                    USED FOR DECOUPLING CRITERIA

!     *PVERVEL*      VERTICAL VELOCITY

!     INPUT PARAMETERS (LOGICAL):

!     *LDNODECP*     TRUE:  NEVER DECOUPLE
!                    FALSE: MAYBE DECOUPLE
!!!   *LDRUNDRY*     TRUE:  RUN PARCEL WITHOUT CONDENSATION
!                    FALSE: RUN PARCEL WITH CONDENSATION

!     OUTPUT PARAMETERS (REAL):

!     *PFPLVL*       PBL PRECIPITATION FLUX AS RAIN                KG/(M**2*S)
!     *PFPLVN*       PBL PRECIPITATION FLUX AS SNOW                KG/(M**2*S)

!     *PUUH*         UPDRAFT X-MOMENTUM
!     *PVUH*         UPDRAFT Y-MOMENTUM
!     *PSLGUH*       UPDRAFT GENERALIZED LIQUID STATIC ENERGY (SLG)
!                    AT HALF LEVEL                                M2/S2
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL   KG/KG
!     *PMFLX*        PBL MASS FLUX                                M/S
!     *PZPLCL*       HEIGHT OF LIFTING CONDENSATION LEVEL OF UPDRAFT          M
!     *PZPTOP*       HEIGHT OF LEVEL OF ZERO KINETIC ENERGY (W=0) OF UPDRAFT  M

!     OUTPUT PARAMETERS (INTEGER):

!     *KPLCL*         FIRST HALF LEVEL ABOVE REAL HEIGHT OF UPRAFT LCL
!     *KPTOP*         HIGHEST HALF LEVEL BELOW PZTOP, AND
!                       UPDRAFT TOP FULL LEVEL (PZTOP IS WITHIN THAT LAYER)
!     *KPLZB*         LEVEL OF UPRAFT ZERO BUOYANCY (LAST FULL LEVEL THAT IS POS. BUOYANT)
!     *KPBLTYPE*    -1: not defined yet
!                    0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: stratocumulus
!                    3: shallow cumulus
!                    4: deep cumulus

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RG       ,RD       ,RCPD     ,RETV     ,RLVTT    ,&
!                     & RLSTT    ,RATM     ,RTT      ,RLMLT
! USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!                     & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
!                     & R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
!                     & RTICECU  ,RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
! USE YOEVDF   , ONLY : RKAP     ,RVDIFTS
! USE YOECUMF  , ONLY : RTAUMEL

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                & RG       ,RD       ,RCPD     ,RETV     ,RLVTT    ,& !yomcst
                & RLSTT    ,RATM     ,RTT      ,RLMLT    ,&           ! -
                & RKAP     ,RVDIFTS  ,&                               !yoevdf
                & RTAUMEL                                             !yoecumf
USE mo_edmf_param   ,ONLY : &
                & FOEEWM   ,&                                         !fcttre.f
                & REPUST                                              !yos_exc

USE mo_vdfpdftable  ,ONLY : vdfpdftable
USE mo_vdfparcel    ,ONLY : vdfparcel
USE mo_vdfstcucrit  ,ONLY : vdfstcucrit
USE mo_vdfbuoysort  ,ONLY : vdfbuoysort

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
!! INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPLCL(KLON,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPTOP(KLON,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPLZB(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFLX(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLGUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQTUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRACB(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZPLCL(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZPTOP(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWUAVG(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDETR(KLON,KLEV)
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PBIR(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRICUI(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTHV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMCU(KLON) 
LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON) 
!ldrundry not used now
!! LOGICAL           ,INTENT(IN)    :: LDRUNDRY(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPBLTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWQT2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE(KLON,KLEV) 
!! REAL(KIND=JPRB)   ,INTENT(IN)    :: PQE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVAR(KLON,KLEV)
!amk
LOGICAL                          :: LDLAND(KLON)
!xxx
!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX2, KLEVX, KFLDX, JCNT
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)
LOGICAL           ,INTENT(IN)    :: LLDIAG

!*         0.2    LOCAL VARIABLES

!--- mean & environmental properties ---
REAL(KIND=JPRB) ::    ZUSTAR  (KLON)     , ZWSTAR(KLON)       , ZKHVFL(KLON)       , &
                    & ZUSIGMA(KLON)      , ZWSIGMA(KLON)      , &
                    & ZSLGENH(KLON,0:KLEV),ZQLENH(KLON,0:KLEV), ZQIENH(KLON,0:KLEV), &
                    & ZQTENH(KLON,0:KLEV), ZUENH(KLON,0:KLEV) , ZVENH(KLON,0:KLEV) , &
                    & ZTVEN(KLON,KLEV)   , ZQTM1 (KLON,KLEV)  , &
                    & ZSLGM1(KLON,KLEV)  , ZMGEOM(KLON,0:KLEV), &
                    & ZTENH(KLON,0:KLEV) , ZRHOH (KLON,0:KLEV), ZTHVEN(KLON,KLEV)

REAL(KIND=JPRB) ::    ZWSTARCAPE(KLON), ZMSEFLUX(KLON), ZMSESFC(KLON), ZWSTARMSE(KLON)

!--- updraft parameters ---
REAL(KIND=JPRB) ::    ZWU2H (KLON,0:KLEV,KDRAFT), ZWUH, &
                    & ZQCUH (KLON,0:KLEV,KDRAFT), ZQUH  (KLON,0:KLEV,KDRAFT), &
                    & ZTUH  (KLON,0:KLEV,KDRAFT), ZEPS  (KLON,0:KLEV,KDRAFT), &
                    & ZFRAC (KLON,0:KLEV,KDRAFT), &
                    & ZBUOF (KLON,KLEV,KDRAFT)  , ZMCLD (KLON)              , &
                    & ZABULK(KLON,0:KLEV) , ZWBULK(KLON,0:KLEV)  , &
                    & ZQTBULK(KLON,0:KLEV), ZSLGBULK(KLON,0:KLEV), &
                    & ZUBULK(KLON,0:KLEV) , ZVBULK(KLON,0:KLEV)  , &
                    & ZZPLZB(KLON,KDRAFT) , ZCAPE1(KLON)
                    
REAL(KIND=JPRB) ::    ZQSATM, ZSATDEF, &
                    & ZUPFLXL(KLON,0:KLEV,KDRAFT), ZUPFLXN(KLON,0:KLEV,KDRAFT), &
                    & ZUPGENL(KLON,KLEV,KDRAFT), ZUPGENN(KLON,KLEV,KDRAFT), &
                    & ZDZRHO, ZPFLXTOT, ZPEVAPUP, ZFAC, ZUPMELT, ZAPRECEVAP

REAL(KIND=JPRB) ::    ZFRACB(KLON,KDRAFT), ZMFLXB(KLON,KDRAFT)

REAL(KIND=JPRB) ::    ZFRACMAX , ZFACMAXEXC , ZFRACTEST , ZFACTESTEXC , &
                    & ZFACEXC(KLON,KDRAFT), ZDUMFRAC, ZDUMR, ZMASSCAPDEPTH, &
                    & ZPDFFACPHI(KLON), ZPDFFACW(KLON), ZLOBUKHOV, &
                    & ZSTABILITY(KLON)

LOGICAL ::            LLDONE(KLON,KDRAFT), LLMASSCAP, LLMCIND(KLON), & 
                    & LLWIPE, LLSTCU


INTEGER(KIND=JPIM) :: JK, JL, JD, ITOP

REAL(KIND=JPRB) ::    ZQEXC   , ZTEXC   , ZDZ     ,  ZDB   , &
                    & ZSPEEDENV         , ZSPEEDUP, &
                    & ZCONS10 , ZCFNC1(KLON,0:KLEV)        , ZTVMEAN     , &
                    & ZRG     , ZMFMAX  , ZMFS(KLON,KDRAFT)

!          REMAINING MODEL PARAMETERS

REAL(KIND=JPRB) ::    ZTAUEPS , ZCLDDEPTH     , &
                    & ZW2THRESH               , ZSTABTHRESH    , ZEISTHRESH , ZBIRTHRESH , &
                    & ZTVLIM  , ZCLDDEPTHDP   , ZDZCLOUD(KLON) , ZW2H       , &
                    & ZZFUNC  , ZZFUNC3(KLON) , ZDHRI(KLON)    , ZDHCL      , &
                    & ZDTHVDZ , ZREPUST       , ZGHM1
                    
REAL(KIND=JPRB) ::    ZZI(KLON)

INTEGER(KIND=JPIM) :: IZI(KLON,KDRAFT)

REAL(KIND=JPRB) ::    ZBUOYCU, ZDTHVCUTOP, ZDMDZ

REAL(KIND=JPRB) ::    ZTAUBM

REAL(KIND=JPRB) ::    ZVERVELCRIT

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

INTERFACE
! #include "surf_inq.h"
END INTERFACE


! #include "vdfparcel.intfb.h"
! #include "vdfstcucrit.intfb.h"
! #include "vdfpdftable.intfb.h"
! #include "vdfbuoysort.intfb.h"

! #include "fcttre.h" ! replaced by use statements



!     ------------------------------------------------------------------

!*         1.     INITIALIZATION
!                 --------------

IF (LHOOK) CALL DR_HOOK('VDFHGHTN',0,ZHOOK_HANDLE)

!-- top % of the PDF associated with the test parcel
ZFRACTEST   = 0.002_JPRB    
CALL VDFPDFTABLE (ZFRACTEST, ZFACTESTEXC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor

!-- total convective area fraction that is done with mass flux
!ZFRACMAX    = 0.075_JPRB     
ZFRACMAX    = 0.10_JPRB     
CALL VDFPDFTABLE (ZFRACMAX, ZFACMAXEXC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor

!-- eddy turnover time scale used in parcel entrainment [s]  (Neggers, Siebesma & Jonker, JAS 2002)
ZTAUEPS     = 400._JPRB

!-- threshold parcel vertical velocity squared [m2/s2]
!ZW2THRESH  = -1._JPRB     
ZW2THRESH   = 0._JPRB      

!-- threshold cloud thickness for stcu/cu transition [m]
ZCLDDEPTH   = 2000._JPRB   

!-- threshold cloud thickness used in shallow/deep decision [m]
ZCLDDEPTHDP = 3000._JPRB   
!ZCLDDEPTHDP = 100000._JPRB   

ZSTABTHRESH = 20._JPRB     ! threshold stability (Klein & Hartmann criteria) [K]
!ZEISTHRESH  = 7.0_JPRB    ! threshold stability (Wood & Bretherton) [K]
ZEISTHRESH  = 10.0_JPRB    ! threshold stability (Wood & Bretherton) [K]
ZBIRTHRESH  = 0.1_JPRB     ! threshold BIR (TKE decoupling criteria) [1]
ZTVLIM      = 0.1_JPRB     ! cloud fraction limit in Tv,env calculation

!MK: convert to use statement
!  CALL SURF_INQ(PREPUST=ZREPUST)
ZREPUST     = REPUST
                 
!-- switch for moist mass flux depth limiter - *experimental*
!LLMASSCAP     = .TRUE.
LLMASSCAP     = .FALSE.    
ZMASSCAPDEPTH = 3000._JPRB 
!ZMASSCAPDEPTH = 50000._JPRB 

!-- switch for applying Klein-Hartmann criterion for stratocumulus --
LLSTCU = .TRUE.
!LLSTCU = .FALSE.

!-- factor used in updraft initialization --
DO JL=KIDIA,KFDIA
!xmk
! ZPDFFACW(JL)   = 0.2_JPRB
! ZPDFFACPHI(JL) = 0.2_JPRB
  ZPDFFACW(JL)   = 0.5_JPRB
  ZPDFFACPHI(JL) = 0.5_JPRB
!xxx
  !ZPDFFACW(JL)   = 1.0_JPRB
  !ZPDFFACPHI(JL) = 1.0_JPRB
ENDDO  

!-- updraft precip evaporation constant
ZAPRECEVAP = 0.001_JPRB
!ZAPRECEVAP = 0.000544_JPRB
!ZAPRECEVAP = 0.0001_JPRB

!-- Betts Miller adjustment timescale --
ZTAUBM = 3600._JPRB * 2._JPRB

!-- critical LS ascent for lateral entrainment limiter  --
ZVERVELCRIT = 0.3_JPRB
  
!-- optimization --
ZRG    = 1.0_JPRB/RG


! set some stuff to zero
DO JL=KIDIA,KFDIA
  
  PWUAVG(JL)     = 0.0_JPRB
  KPBLTYPE(JL)   = -1            ! -1 means: yet unknown
  
  ZZI(JL)        = 0._JPRB       !mixed layer scalings
  ZWSTAR(JL)     = 0._JPRB        
  ZWSTARCAPE(JL) = 0._JPRB        
  ZWSTARMSE(JL)  = 0._JPRB        
  
  PRICUI(JL   )  = 1._JPRB       ! 1 / cumulus inversion Richardson number
  PDTHV(JL)      = 0._JPRB
   
  ZCAPE1(JL)     = 0._JPRB

  ZMCLD(JL)      = 0._JPRB
  PMCU(JL)       = 0._JPRB       ! cloud-depth average moist updraft mass flux
  
ENDDO


DO JD=1,KDRAFT
  DO JL=KIDIA,KFDIA
    PZPLCL(JL,JD)  = -100._JPRB  ! default value: -100 (no LCL)
    PZPTOP(JL,JD)  = 0._JPRB     
    KPLCL(JL,JD)   = 0           ! default value: 0 (no PBL cloud)
    KPTOP(JL,JD)   = 0          
    KPLZB(JL,JD)   = 0          
    LLDONE(JL,JD)  = .TRUE.      ! default: TRUE (don't launch the parcel)
    ZFRACB(JL,JD)  = 0._JPRB 
    PFRACB(JL,JD)  = 0._JPRB 
    ZFACEXC(JL,JD) = 0._JPRB 
    ZFACEXC(JL,JD) = 0._JPRB 
    ZMFLXB(JL,JD)  = 0._JPRB 
    IZI(JL,JD)     = 0
  ENDDO
ENDDO


DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    PWQT2(JL,JK) = 0._JPRB  
  ENDDO
ENDDO


DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PDETR(JL,JK) = 0._JPRB
  ENDDO
ENDDO    


!--- parcel half level parameters ---
DO JD=1,KDRAFT
  DO JK=0,KLEV
    DO JL=KIDIA,KFDIA
    PUUH(JL,JK,JD)  = 0.0_JPRB
    PVUH(JL,JK,JD)  = 0.0_JPRB
    PSLGUH(JL,JK,JD)= 0.0_JPRB
    PQTUH(JL,JK,JD) = 0.0_JPRB
    PMFLX(JL,JK,JD) = 0.0_JPRB
    ZTUH(JL,JK,JD)  = 0.0_JPRB
    ZQUH(JL,JK,JD)  = 0.0_JPRB
    ZQCUH(JL,JK,JD) = 0.0_JPRB
    ZEPS(JL,JK,JD)  = 0.0_JPRB
    ZWU2H(JL,JK,JD) = 0.0_JPRB
    ZFRAC(JL,JK,JD) = 0.0_JPRB
    ZUPFLXL(JL,JK,JD)  = 0.0_JPRB
    ZUPFLXN(JL,JK,JD)  = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO


!--- parcel full level parameters ---
DO JD=1,KDRAFT
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZBUOF(JL,JK,JD)    = 0.0_JPRB
      ZUPGENL(JL,JK,JD)  = 0.0_JPRB
      ZUPGENN(JL,JK,JD)  = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO


!     -----------------------------------------------------------------
!
!*         2.     PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!*                OF CONSERVED VARIABLES
!                 -----------------------------------------------------

!*         2.1  FULL LEVEL CPM, SLG, QT AND TV
!*

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZSLGM1(JL,JK) = RCPD * PTM1(JL,JK) + PGEOM1(JL,JK) &
                  & - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)  
      ZQTM1 (JL,JK) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)

!          parcel goes through cloud portion of environment
!          (added ql loading; ql,cld=ql,mean/fc; qv = qsat) 
!          safety: fc>0.1; linear interpolation between overcast 
!                  and cloudy portion for 0<fc<0.1
!                  guaranteed to be < tv from mean conditions

!          grid box mean virtual effect
      ZTVMEAN       = PTM1(JL,JK) * ( 1.0_JPRB + RETV * PQM1(JL,JK) &
                  & - PLM1(JL,JK) - PIM1(JL,JK) )       !qli loading  
      ZTVEN(JL,JK) = ZTVMEAN
      ZTHVEN(JL,JK) = ( PAPM1(JL,JK)/RATM )**(-RD/RCPD) * ZTVEN(JL,JK)
    ENDDO
  ENDDO


!*         2.2  HALF-LEVEL ENVIRONMENT INTERPOLATION (QT, QL, QI, SLG)
!*

  DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA

    IF (JK==1) THEN
      ZGHM1 = PGEOH(JL,JK) + 50000._JPRB*RG   !avoid using top half level (=inf)
    ELSE
      ZGHM1 = PGEOH(JL,JK-1)
    ENDIF  
    
    ZQTENH(JL,JK) = ( ZQTM1(JL,JK+1) *(ZGHM1-PGEOH(JL,JK  )) &
                & +   ZQTM1(JL,JK)   *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZQLENH(JL,JK) = ( PLM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PLM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZQIENH(JL,JK) = ( PIM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PIM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZSLGENH(JL,JK)= ( ZSLGM1(JL,JK+1)*(ZGHM1-PGEOH(JL,JK  )) &
                & +   ZSLGM1(JL,JK)  *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZUENH(JL,JK)  = ( PUM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PUM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZVENH(JL,JK)  = ( PVM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PVM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))

    !calculate T at half levels from sl, for later use in density calculations                
    !ZTENH        = ( ZSLGENH (JL,JK) - PGEOH(JL,JK) &
    !                   & + RLVTT*ZQLENH(JL,JK) + RLSTT*ZQIENH(JL,JK) &
    !                   & ) / RCPD
    
    ZTENH(JL,JK)  =  ( PTM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                 & +   PTM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                 &   )                /(ZGHM1-PGEOH(JL,JK+1))
                
    !air density at half levels
    ZRHOH(JL,JK) = PAPHM1(JL,JK)/(RD*ZTENH(JL,JK))
      
  ENDDO
  ENDDO



!     -----------------------------------------------------------------

!*         3.     RELEASE THE FIRST (TEST) UPDRAFT TO GET PBL HEIGHTS
  
  
  !* set updraft index to 1
  JD = 1   

  DO JL=KIDIA,KFDIA
    
 
    PFRACB(JL,JD) = ZFRACTEST
 
 
    !* 3.1    DETERMINE STABILITY OF BL USING THE SURFACE BUOYANCY FLUX
    !*
    ZKHVFL(JL)  = ( 1.0_JPRB + RETV *  ZQTM1(JL,KLEV) ) * PKHFL(JL) + &
                & ( RETV * ZSLGM1(JL,KLEV) / RCPD )     * PKQFL(JL) 


    IF ( ZKHVFL(JL) >= 0.0_JPRB ) THEN
      
      ! stable BL (no updrafts expected/needed)
      KPBLTYPE(JL)  = 0

    ELSE

      LLDONE(JL,JD) = .FALSE.  !confirm launch
     
     
      !* 3.2    SURFACE LAYER SCALING
      !*
      ZUSTAR(JL)  = MAX( SQRT(PKMFL(JL)), ZREPUST )               !u* (repust=10e-4)
      ZWSTAR(JL)  = (- ZKHVFL(JL) * RG / PTM1(JL,KLEV) * 1000._JPRB ) &   !zi=1000m
                       & ** ( 1._JPRB/3._JPRB) 
      ZWSIGMA(JL)      = 1.2_JPRB &
       & * ( ZUSTAR(JL)**3 &
       & - 1.5_JPRB * RKAP * ZKHVFL(JL) * PGEOH(JL,KLEV-1) / PTM1(JL,KLEV-1) &
       & ) ** ( 1.0_JPRB/3._JPRB )                         ! Kolmogorov 1/3-power
      ZUSIGMA(JL)      = 2.29_JPRB  &
       & * ( ZUSTAR(JL)**3 &
       & + 0.5_JPRB / 12.0_JPRB * RKAP * ZWSTAR(JL)**3 &
       &   ) ** ( 1._JPRB/3._JPRB)

      !  scaling factors between phi* and initial updraft phi excesses
      !ZLOBUKHOV = -ZUSTAR(JL)**3 * PTM1(JL,KLEV) / (RG * RKAP * ZKHVFL(JL))
      !ZPDFFACPHI(JL) = ZRG*PGEOH(JL,KLEV-1)/ZLOBUKHOV
      
      
      !* 3.3    INITIALIZE TEST UPDRAFT
      !*
      
      !get the constant associated with the top ZFRACTEST % of the PDF
      ZFACEXC(JL,1) = ZFACTESTEXC

!amk ... optional 2x ZFACEXC: stronger parcels (1: test and 3: cloudy)
!     ZFACEXC(JL,1) = ZFACEXC(JL,1) * 2.0_JPRB
!xxx      

      !calculate the initial excess values
      ZWU2H(JL,KLEV-1,JD) = ( ZPDFFACW(JL) * ZFACEXC(JL,1) * ZWSIGMA(JL) )**2         
      ZTEXC               = - ZPDFFACPHI(JL) * ZFACEXC(JL,1) * PKHFL(JL) / ZWSIGMA(JL) 
      ZQEXC               = - ZPDFFACPHI(JL) * ZFACEXC(JL,1) * PKQFL(JL) / ZWSIGMA(JL) 
      ZTEXC               = MAX(ZTEXC, 0.0_JPRB)
      ZQEXC               = MAX(ZQEXC, 0.0_JPRB)
      ZWU2H(JL,KLEV-1,JD) = MIN(ZWU2H(JL,KLEV-1,JD),100.0_JPRB) !10 m/s  limit (pure safety)
      ZTEXC               = MIN(ZTEXC              , 10.0_JPRB) !10 K    limit
      ZQEXC               = MIN(ZQEXC              , 0.01_JPRB) !10 g/kg limit
      PQTUH(JL,KLEV-1,JD) = ZQTENH(JL,KLEV-1) + ZQEXC
      ZQCUH(JL,KLEV-1,JD) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
      ZQUH (JL,KLEV-1,JD) = PQTUH(JL,KLEV-1,JD)  - ZQCUH(JL,KLEV-1,JD)
      PSLGUH(JL,KLEV-1,JD)= ZSLGENH(JL,KLEV-1) + RCPD * ZTEXC
      ZTUH (JL,KLEV-1,JD) = ( PSLGUH (JL,KLEV-1,JD) - PGEOH(JL,KLEV-1) &
                       & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1) &
                       & ) / RCPD

! ... u & v: (wind speed assumed to be negatively correlated with T and q excesses)
      ZSPEEDENV        = SQRT( ZUENH(JL,KLEV-1)**2 + ZVENH(JL,KLEV-1)**2 )
!     ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC * ZUSIGMA(JL), 0._JPRB )
!     ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC * ZUSTAR(JL)**2/ ZWSIGMA(JL) , 0._JPRB )

!      PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
!      PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
      PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1)  !mean wind at this half level
      PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1)
!      PUUH(JL,KLEV-1,JD)= PUM1(JL,KLEV)     !10m wind instead of 20m
!      PVUH(JL,KLEV-1,JD)= PVM1(JL,KLEV)
      
    ENDIF

  ENDDO !JL


  !* 3.4   RELEASE THE TEST UPDRAFT #1
  !*          - USED TO MAKE A FIRST GUESS OF THE HEIGHTS OF CLOUD BASE & INVERSION,
  !*            AND TO DETERMINE PBL TYPE.
  !*
  CALL VDFPARCEL (KIDIA   , KFDIA   , KLON    , KLEV    , KDRAFT  , &
                & PGEOH   , PGEOM1  , PAPHM1  , &
                & PUM1    , PVM1    , ZQTM1   , ZSLGM1  , ZTVEN   , &
                & PUUH    , PVUH    , PSLGUH  , PQTUH   , ZWU2H   , ZQCUH  , ZBUOF , & 
                & ZQUH    , ZTUH    , ZEPS    , ZFACEXC , &
                & PZPLCL  , KPLCL   , PZPTOP  , KPTOP   , KPLZB   , &
                & JD      , ZUPGENL , ZUPGENN , &
!amk: for convective preconditioning
                & PVAR    , &
!xxx
!amk
                & LDLAND , &   
!xxx
                & ZTAUEPS , PVERVEL , ZW2THRESH, LLDONE , KPBLTYPE)  



!     -----------------------------------------------------------------

!*         4.     CLASSIFICATION OF THE CONVECTIVE PBL
!                 ------------------------------------


  !* 4.1    CLASSIFY THE CONVECTIVE PBL
  !*
  DO JL=KIDIA,KFDIA
    IF ( KPBLTYPE(JL)/=0 ) THEN
 
!xmk  IF ( PZPLCL(JL,1) > PZPTOP(JL,1) .OR. KPLCL(JL,1) == 0 ) THEN
      IF ( PZPLCL(JL,1) >= PZPTOP(JL,1) .OR. KPLCL(JL,1) == 0 ) THEN
!xxx
      
        !dry convective PBL
        KPBLTYPE(JL)  = 1                   !dry convective PBL
        ZDZCLOUD(JL)  = 0.0_JPRB            !cloud thickness

      ELSE

        !moist convective PBL
        ZDZCLOUD(JL)  = PZPTOP(JL,1) - PZPLCL(JL,1) !cloud thickness
            
        IF (ZDZCLOUD(JL)>ZCLDDEPTHDP .AND. .NOT.LLMASSCAP ) THEN
        
          !deep convection
          KPBLTYPE(JL) = 4
          
        ELSE
        
          IF (LLSTCU) THEN
            KPBLTYPE(JL) = 2   !set the type to stratocumulus for the moment
          ELSE
            KPBLTYPE(JL) = 3   !RN run without Klein-Hartmann criterion!
          ENDIF  
          
        ENDIF
        
      ENDIF

    ENDIF !KPBLTYPE /=0
  ENDDO !JL
  
  
  
  !* 4.2    CHECK THE STRATOCUMULUS/SHALLOW CUMULUS CRITERION (TRIGGER FUNCTION)
  !*        IF SHALLOW CUMULUS IS DIAGNOSED, KPBLTYPE WILL BE SET TO 3
  !*
  CALL VDFSTCUCRIT ( KIDIA   , KFDIA   , KLON  , KLEV  , KDRAFT  , &
                  &  PTM1    , ZSLGM1  , ZQTM1 , PAPM1 , PGEOM1  , &
                  &  ZSTABTHRESH, ZEISTHRESH, ZCLDDEPTH, ZBIRTHRESH, &
                  &  ZDZCLOUD, PZPLCL  , &
                  &  KPTOP   , KPBLTYPE, LDNODECP, &
                  &  ZSTABILITY )
    
  
!     -----------------------------------------------------------------

!*         5.     CLOSURE FOR ORGANIZED UPDRAFTS (JD=2,3)
!                 ---------------------------------------


  !* 5.1    MIXED LAYER SCALINGS
  !*

  DO JL=KIDIA,KFDIA
    
    IF ( KPBLTYPE(JL)/=0 ) THEN    !don't do this for stable PBL
      
      !--- Mixed layer scaling depth ---
      SELECT CASE (KPBLTYPE(JL))
    
        CASE(1)
          !Dry convective PBL - Inversion height
          ZZI(JL)   = PZPTOP(JL,1)
          IZI(JL,1) = KPTOP(JL,1)
                
        CASE(2)
          !Stratocumulus - Inversion height
          !CAUTION: During decoupling in the intermediate regime (e.g. ASTEX/ATEX) the
          !   relevant ML scaling height changes from PBL inversion to level of minimum
          !   buoyancy flux. In the current setup this is not modelled yet!
          ZZI(JL)   = PZPTOP(JL,1)
          IZI(JL,1) = KPTOP(JL,1)
          
        CASE(3)
          !Shallow cumulus - Level of minimum buoyancy flux
          !Assume that the moist updraft LCL is very close to this level
          ZZI(JL)   = PZPLCL(JL,1)
          IZI(JL,1) = KPLCL(JL,1)
          
        CASE(4)
          !Deep cumulus - Only do a dry parcel up to cloud base
          ZZI(JL)   = PZPLCL(JL,1)
          IZI(JL,1) = KPLCL(JL,1)
                
      END SELECT

      !--- Mixed layer convective velocity scale ---
      ZWSTAR(JL) = ( -ZKHVFL(JL) * RG * ZZI(JL) / ZTHVEN(JL,KLEV)  ) ** (1._JPRB/3._JPRB)
      
    ENDIF
    
  ENDDO  

    
    
  !*  5.2    RI NUMBER OF CUMULUS INVERSION
  !* 
  
  !-- Test-updraft cloudy CAPE --
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA
      IF ( ZQCUH(JL,JK,1)>0._JPRB .AND. ZBUOF(JL,JK,1)>0._JPRB .AND. JK<=KPLCL(JL,1) ) THEN  
        ZDZ = ZRG*( PGEOH(JL,JK-1) - PGEOH(JL,JK) )
        ZCAPE1(JL) = ZCAPE1(JL) + ZDZ * ZBUOF(JL,JK,1) 
      ENDIF
    ENDDO
  ENDDO
  
  DO JL=KIDIA,KFDIA
    IF ( KPBLTYPE(JL)==2 .OR. KPBLTYPE(JL)==3 ) THEN
    
      !-- Interpolate LZB height --
      ZZPLZB(JL,1) = PZPTOP(JL,1)
      IF (KPLZB(JL,1)>2) THEN
        ZDZ = (PGEOH(JL,KPLZB(JL,1)-1) - PGEOH(JL,KPLZB(JL,1)))*ZRG
        ZDB = ZBUOF(JL,KPLZB(JL,1),1)  - ZBUOF(JL,KPLZB(JL,1)-1,1)
        IF (ZDB>0._JPRB) THEN
          ZZPLZB(JL,1) = PGEOH(JL,KPLZB(JL,1)) * ZRG + &
                       & ZDZ * ZBUOF(JL,KPLZB(JL,1),1) / ZDB 
          ZZPLZB(JL,1) = MIN( ZZPLZB(JL,1), PGEOH(JL,KPLZB(JL,1)-1)*ZRG )
        ENDIF               
      ENDIF
      
      !-- Cloud layer average test-updraft positive buoyancy --
      ZBUOYCU = 0._JPRB
      IF ( ZZPLZB(JL,1)-PZPLCL(JL,1)>0._JPRB ) THEN  
        ZBUOYCU = ZCAPE1(JL) / ( ZZPLZB(JL,1) - PZPLCL(JL,1) )
      ENDIF  
  
      !-- Inversion theta_v jump --
      !JK = KPLZB(JL,1)    !use level of zero buoyancy (LZB) of test-updraft
      JK = KPTOP(JL,1)    !use top level of test-updraft
      IF (JK>2) THEN
        ZDTHVCUTOP = MAX( ZTHVEN(JL,JK-1)-ZTHVEN(JL,JK), ZTHVEN(JL,JK)-ZTHVEN(JL,JK+1) )
      ENDIF  
      
      !-- Cumulus Ri number - used again in VDFEXCU --
      IF ( ZDTHVCUTOP > 0._JPRB ) THEN 
        PRICUI(JL) = ZBUOYCU * ZRG * ZTHVEN(JL,KLEV) / ZDTHVCUTOP   
      ENDIF  
      
      !-- RN testing: no top entrainment for stcu (yikes) --
      !IF ( ZSTABILITY(JL) > ZSTABTHRESH ) THEN
      !  PRICUI(JL) = 0._JPRB
      !ENDIF

    ENDIF
  ENDDO
  
  

  !* 5.3    CLOSURE OF UPDRAFT AREA FRACTIONS (JD=2,3)
  !*
    
  DO JL=KIDIA,KFDIA
    
    IF ( KPBLTYPE(JL)/=0 ) THEN    !don't do this for stable PBL


      IF (KPBLTYPE(JL)>1) THEN
      
        !--- Transition layer depth, Scale I: Dry entrainment layer depth ---
        !       using thv gradient averaged over a number of layers above h
        !
        ZDTHVDZ  = MAX( 0.01_JPRB, ZTHVEN(JL,IZI(JL,1)-2) - ZTHVEN(JL,IZI(JL,1)) ) * RG / &
               & ( PGEOM1(JL,IZI(JL,1)-2) - PGEOM1(JL,IZI(JL,1)) )

        ZWSTARCAPE(JL) = MAX( 0._JPRB, ZCAPE1(JL) )**0.5_JPRB
        
        ZMSEFLUX(JL)   = MAX( 0._JPRB, - RCPD * PKHFL(JL) - RLVTT * PKQFL(JL) )
        ZMSESFC(JL)    = RCPD * PTM1(JL,KLEV) + RLVTT * PQM1(JL,KLEV)
        ZWSTARMSE(JL)  = ( ZMSEFLUX(JL) * RG / ZMSESFC(JL) * ZZI(JL) ) &
                         & ** ( 1._JPRB/3._JPRB) 
        
        ZW2H      = 0.5_JPRB * ZWSTAR(JL)**2.
!xmk
!       ZW2H      = 1.0_JPRB * ZWSTAR(JL)**2.
!xxx
        !ZW2H      = 0.5_JPRB * MAX( ZWSTAR(JL), ZWSTARCAPE(JL) )**2.
        !ZW2H      = 0.5_JPRB * MAX( ZWSTAR(JL), ZWSTARMSE(JL) )**2.
        
        ZDHRI(JL) = (  ZW2H * ZRG * ZTHVEN(JL,KLEV) / (0.5_JPRB * ZDTHVDZ)  )**0.5
        PDTHV(JL) = ZDHRI(JL) * ZDTHVDZ   !used again in VDFEXCU
      
        ZZFUNC3(JL) = MAX( 0._JPRB, 0.2_JPRB * ZDHRI(JL) / PZPLCL(JL,1) )   !cy32r3
      
      ENDIF
      
      
IF ( LLDIAG ) THEN
    PEXTRA(JL,32,41) = 0.0_JPRB
    PEXTRA(JL,33,41) = 0.0_JPRB
    PEXTRA(JL,34,41) = 0.0_JPRB
    PEXTRA(JL,35,41) = 0.0_JPRB
    PEXTRA(JL,36,41) = 0.0_JPRB
    PEXTRA(JL,37,41) = 0.0_JPRB
    PEXTRA(JL,38,41) = 0.0_JPRB
    PEXTRA(JL,39,41) = 0.0_JPRB
    PEXTRA(JL,40,41) = 0.0_JPRB
    PEXTRA(JL,41,41) = 0.0_JPRB
    PEXTRA(JL,42,41) = 0.0_JPRB
    PEXTRA(JL,43,41) = 0.0_JPRB
    PEXTRA(JL,44,41) = 0.0_JPRB
ENDIF

      !--- Calculation of moist updraft area fraction ---
      SELECT CASE (KPBLTYPE(JL))
    
        CASE(1)
        
          !Dry convective PBL
          !Set area fraction of moist group to zero
          ZFRACB(JL,3) = 0._JPRB

        CASE(2)
        
          !Stratocumulus
          !Set area fraction of moist group to ZFRACMAX
          ZFRACB(JL,3) = ZFRACMAX

        CASE(3)
          
          !Shallow cumulus
          !Flexible updraft area fractions

          !-- Transition layer depth, Scale II: Cumulus condensation depth-scale --
          !
          ZDHCL = MIN(200._JPRB, 0.1_JPRB * ZDZCLOUD(JL))
          ZDHCL = MAX(ZDHCL,0._JPRB)
          
          !Use M/w* ~ Dh / h  (Neggers et al., QJ, 2007)
          IF (PZPLCL(JL,1).GT.0._JPRB) THEN
!xmk
!           ZZFUNC =  0.15_JPRB * ( ZDHCL/PZPLCL(JL,1) )   
!---
            ZZFUNC =  0.25_JPRB * ( ZDHCL/PZPLCL(JL,1) )   
!xxx
          ELSE
            ZZFUNC = 0._JPRB
          ENDIF
                                  
          !-- Choose the minimum of scales I and II --
          ZFRACB(JL,3)  = MIN( ZFRACMAX, ZZFUNC, ZZFUNC3(JL) )
              
          !-- 1st guess for the LCL mass flux --
          ZMFLXB(JL,3)  = ZWSTAR(JL) * ZFRACB(JL,3) * ZRHOH(JL,KPLCL(JL,1))
          
          !-- Switch KPBLTYPE to dry convective if moist updraft is not launched --
          IF (ZFRACB(JL,3).EQ.0._JPRB) THEN
            KPBLTYPE(JL)=1
          ENDIF
          
IF ( LLDIAG ) THEN
    PEXTRA(JL,32,41) = ZZFUNC         ! cum. fraction: 0.1cloud scale<200 / h * 0.25
    PEXTRA(JL,33,41) = ZZFUNC3(JL)    ! cum. fraction: transition layer depth closure
    PEXTRA(JL,34,41) = ZFRACB(JL,3)   ! cum. fraction: final
    PEXTRA(JL,35,41) = ZW2H
    PEXTRA(JL,36,41) = ZDTHVDZ
    ZDTHVDZ  = MAX( 0.01_JPRB, ZTHVEN(JL,IZI(JL,1)-1)-ZTHVEN(JL,IZI(JL,1)) ) * RG / &
           & ( PGEOM1(JL,IZI(JL,1)-1) - PGEOM1(JL,IZI(JL,1)) )
    PEXTRA(JL,37,41) = ZDTHVDZ
    ZDTHVDZ  = MAX( 0.01_JPRB, ZTHVEN(JL,IZI(JL,1)-3)-ZTHVEN(JL,IZI(JL,1)) ) * RG / &
           & ( PGEOM1(JL,IZI(JL,1)-3) - PGEOM1(JL,IZI(JL,1)) )
    PEXTRA(JL,38,41) = ZDTHVDZ
    IF (ZFRACMAX    .EQ. MIN(ZFRACMAX, ZZFUNC, ZZFUNC3(JL)) ) THEN
      PEXTRA(JL,39,41) = 1.0_JPRB
    ENDIF
    IF (ZZFUNC      .EQ. MIN(ZFRACMAX, ZZFUNC, ZZFUNC3(JL)) ) THEN
      PEXTRA(JL,40,41) = 1.0_JPRB
    ENDIF
    IF (ZZFUNC3(JL) .EQ. MIN(ZFRACMAX, ZZFUNC, ZZFUNC3(JL)) ) THEN
      PEXTRA(JL,41,41) = 1.0_JPRB
    ENDIF
    IF (200._JPRB   .EQ. MIN(200._JPRB, 0.1_JPRB * ZDZCLOUD(JL))) THEN
      PEXTRA(JL,42,41) = 1.0_JPRB
    ENDIF
    IF (0.1_JPRB * ZDZCLOUD(JL) .EQ. MIN(200._JPRB, 0.1_JPRB * ZDZCLOUD(JL))) THEN
      PEXTRA(JL,43,41) = 1.0_JPRB
    ENDIF
    PEXTRA(JL,44,41) = ZDZCLOUD(JL)
ENDIF


        CASE(4)
        
          !Deep cumulus
          !Set area fraction of moist group to zero (only allow updraft transport in dry mixed layer)
          ZFRACB(JL,3) = 0._JPRB

      END SELECT !KPBLTYPE
      
      
      !Dry updraft area fraction (JD=2):  a_dry = 0.1 - a_moist
      ZFRACB(JL,2) = MAX( 0._JPRB, ZFRACMAX - ZFRACB(JL,3) )
      
      PFRACB(JL,2) = ZFRACB(JL,2)
      PFRACB(JL,3) = ZFRACB(JL,3)


    ENDIF !KPBLTYPE /=0
    
  ENDDO !JL


      
!     -----------------------------------------------------------------

!*         6.     CALCULATE VERTICAL PROFILES OF ALL UPDRAFTS (JD=2,3)
!                 ----------------------------------------------------


  !*       6.1    CALCULATE THE SCALING FACTORS OF THE UPDRAFT EXCESS WITH THE SURFACE JOINT PDFS
  !*
  DO JD = 2,KDRAFT
    DO JL=KIDIA,KFDIA
      
      IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB ) THEN
        
        !-- Get the PDF scaling factor --
        SELECT CASE (JD)
            
        CASE(2)
            !lower part of top ZFRACMAX %
            ZDUMFRAC = ZFRACMAX - ZFRACB(JL,2)
            CALL VDFPDFTABLE(ZDUMFRAC , ZFACEXC(JL,2), ZDUMR, ZDUMR, 0)
            ZFACEXC(JL,2) = ( ZFRACMAX * ZFACMAXEXC - ZDUMFRAC * ZFACEXC(JL,2) ) / ZFRACB(JL,2)
                    
          CASE(3)
            !upper part of top ZFRACMAX %
            ZDUMFRAC = ZFRACB(JL,JD)
            CALL VDFPDFTABLE(ZDUMFRAC , ZFACEXC(JL,3), ZDUMR, ZDUMR, 0)

!amk ... optional 2x ZFACEXC: stronger parcels (1: test and 3: cloudy)
!           ZFACEXC(JL,3) = ZFACEXC(JL,3) * 2.0_JPRB
!xxx
            
          END SELECT
        
      ENDIF !KPBLTYPE & ZFRACB

    ENDDO !JL
  ENDDO !JD
    
    
  !*       6.2    VERTICAL INTEGRATION OF DRY & MOIST UPDRAFT BUDGETS (JD=2,3)
  !*
  DO JD = 2,KDRAFT
  
    !-- Initialize updraft --
    DO JL=KIDIA,KFDIA
      
      IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB ) THEN
        
        LLDONE(JL,JD) = .FALSE. !confirm launch
          
        ZWU2H(JL,KLEV-1,JD) = ( ZPDFFACW(JL) * ZFACEXC(JL,JD) * ZWSIGMA(JL) )**2 
        ZTEXC            = - ZPDFFACPHI(JL) * ZFACEXC(JL,JD) * PKHFL(JL) / ZWSIGMA(JL) 
        ZQEXC            = - ZPDFFACPHI(JL) * ZFACEXC(JL,JD) * PKQFL(JL) / ZWSIGMA(JL) 
        ZTEXC            = MAX(ZTEXC, 0.0_JPRB)
        ZQEXC            = MAX(ZQEXC, 0.0_JPRB)
        PQTUH(JL,KLEV-1,JD) = ZQTENH(JL,KLEV-1) + ZQEXC
        ZQCUH(JL,KLEV-1,JD) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
        ZQUH (JL,KLEV-1,JD) = PQTUH(JL,KLEV-1,JD)  - ZQCUH(JL,KLEV-1,JD)
        PSLGUH(JL,KLEV-1,JD)= ZSLGENH(JL,KLEV-1) + RCPD * ZTEXC
        ZTUH (JL,KLEV-1,JD) = ( PSLGUH (JL,KLEV-1,JD) - PGEOH(JL,KLEV-1) &
                       & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1) &
                       & ) / RCPD

!   ... u & v: (wind speed assumed to be negatively correlated with T and q excesses)
        ZSPEEDENV        = SQRT( ZUENH(JL,KLEV-1)**2 + ZVENH(JL,KLEV-1)**2 )
!       ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC * ZUSIGMA(JL), 0._JPRB )
!       ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC * ZUSTAR(JL)**2/ ZWSIGMA(JL) , 0._JPRB )

!        PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
!        PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
        PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1)  !mean wind at this half level
        PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1)
!        PUUH(JL,KLEV-1,JD)= PUM1(JL,KLEV)     !10m wind instead of 20m
!        PVUH(JL,KLEV-1,JD)= PVM1(JL,KLEV)
 
      ENDIF !KPBLTYPE & ZFRACB
      
    ENDDO !JL
    
    
    
    
    !-- Release the updraft --
    CALL VDFPARCEL (KIDIA   , KFDIA   , KLON    , KLEV    , KDRAFT  , &
                  & PGEOH   , PGEOM1  , PAPHM1  , &
                  & PUM1    , PVM1    , ZQTM1   , ZSLGM1  , ZTVEN   , &
                  & PUUH    , PVUH    , PSLGUH  , PQTUH   , ZWU2H   , ZQCUH  , ZBUOF , & 
                  & ZQUH    , ZTUH    , ZEPS    , ZFACEXC , &
                  & PZPLCL  , KPLCL   , PZPTOP  , KPTOP   , KPLZB   , &
                  & JD      , ZUPGENL , ZUPGENN , &
                  & PVAR    , &
                  & LDLAND , &   
                  & ZTAUEPS , PVERVEL , ZW2THRESH, LLDONE , KPBLTYPE)  
    
  ENDDO !JD 



  !*        6.3 SOME PROTECTIONAL MEASURES AGAINST UPDRAFT #3 FAILING (ONLY JUST) TO REACH CONDENSATION.
  !*
  DO JL=KIDIA,KFDIA

    IF ( KPBLTYPE(JL)==2 .OR. KPBLTYPE(JL)==3) THEN
      IF ( ZFRACB(JL,3)>0._JPRB .AND. PZPLCL(JL,3)<0._JPRB) THEN
      
        !set cloud base height to updraft #3 top (cloud layer exists but has zero depth)
        KPLCL(JL,3)  = KPTOP(JL,3)
        PZPLCL(JL,3) = PZPTOP(JL,3)
        
      ENDIF
      
    ENDIF
     
  ENDDO !JL 


    
  !*        6.4  LIMITER FOR CLOUDY DEPTH OF MOIST UPDRAFT  *TESTING*
  !* 
  IF (LLMASSCAP) THEN
  
    DO JL=KIDIA,KFDIA
      LLMCIND(JL) = .FALSE.
      IF (KPBLTYPE(JL)==3) THEN
        IF ( PZPLCL(JL,3)>0._JPRB .AND. (PZPTOP(JL,3) - PZPLCL(JL,3))>ZMASSCAPDEPTH ) THEN
          PZPTOP(JL,3) = PZPLCL(JL,3) + ZMASSCAPDEPTH
          LLMCIND(JL) = .TRUE.
        ENDIF
      ENDIF
    ENDDO
    
    DO JK=KLEV-1,1,-1
      DO JL=KIDIA,KFDIA
      
        IF ( LLMCIND(JL) ) THEN
          IF ( PGEOH(JL,JK+1)*ZRG<=PZPTOP(JL,3) .AND. PGEOH(JL,JK)*ZRG>PZPTOP(JL,3) ) THEN
            KPTOP(JL,3) = JK+1
          ENDIF
        ENDIF
        
      ENDDO
    ENDDO
    
  ENDIF


  !*        6.5  UPDRAFT PRECIPITATION FLUXES (RAIN AND SNOW)
  !*             (neglected: rain water loading, see cuascn)
  !     
  DO JD = 3,KDRAFT  !moist updrafts only
    
    DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
      
        ZDZRHO = ZRG * ( PAPHM1(JL,JK)-PAPHM1(JL,JK-1) ) 
        
        !-- Add precip generation to flux [kg/m2/s: tendency * layer depth * air density] --
        ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK-1,JD) + ZFRACB(JL,JD) * ZUPGENL(JL,JK,JD)
        ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK-1,JD) + ZFRACB(JL,JD) * ZUPGENN(JL,JK,JD)
        
        !-- Melting at freezing level (snow->rain) --
        !   - independent of precipitation flux
        !   - cooling by melting is done in vdfmain (chp 6)
        IF (ZUPFLXN(JL,JK,JD)>0._JPRB .AND. PTM1(JL,JK) > RTT) THEN
          ZUPMELT = (1.0_JPRB+0.5_JPRB*(PTM1(JL,JK)-RTT)) * &   ! local rate [kg/kg /s]
                  & (PTM1(JL,JK)-RTT) * RCPD/(RLMLT*RTAUMEL)
          ZUPMELT = ZUPMELT * ZFRACB(JL,JD) * ZDZRHO            ! grid mean flux change [kg/m2/s]
          ZUPMELT = MIN(      ZUPFLXN(JL,JK,JD),  ZUPMELT )
          ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK,JD) - ZUPMELT
          ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK,JD) + ZUPMELT
        ENDIF

        !-- Saturation deficit of mean state T --
         ZQSATM = FOEEWM(PTM1(JL,JK))/PAPM1(JL,JK)
         ZQSATM = MIN(0.5_JPRB,ZQSATM)
         ZQSATM = ZQSATM/(1.0_JPRB-RETV*ZQSATM)
!xmk     ZSATDEF = MAX( 0._JPRB, ZQSATM-PQM1(JL,JK) )

!amk: assume that precipitation falls into environment that is premoistened
!     (on a mixing line just between cloud edge (saturated ) and environment mean) 
         ZSATDEF = MAX( 0._JPRB, (ZQSATM-PQM1(JL,JK))/2.0_JPRB )
!xxx

      !cucalln: (does all levels)
      !  IFLAG=1
      !  CALL SATUR (KIDIA , KFDIA , KLON  , NJKT2 , KLEV,&
      !            & PAP   , ZTP1  , ZQSAT , IFLAG  )  

        ZPFLXTOT = ZUPFLXL(JL,JK,JD) + ZUPFLXN(JL,JK,JD)
!xmk: RH<80%, below cloud base 
        IF ( ZPFLXTOT > 0._JPRB ) THEN
!       IF ( ZPFLXTOT > 0._JPRB .AND. JK > KPLCL(JL,3) .AND. PQM1(JL,JK)/ZQSATM < 0.8_JPRB ) THEN 
!xxx         
          !-- Precip evaporation tendency (Kessler 1969, Tiedtke 1993) --
          ZPEVAPUP = ZAPRECEVAP * ZSATDEF * &                   ! local rate [kg/kg /s]
              & ( ZPFLXTOT / ZFRACB(JL,JD) / 0.00509_JPRB * &
              &   SQRT( PAPM1(JL,JK)/PAPHM1(JL,KLEV) ) &
              & )**0.5777_JPRB
          ZPEVAPUP = ZPEVAPUP * ZFRACB(JL,JD) * ZDZRHO          ! grid mean flux change [kg/m2/s]
          !-- Back-partition evaporation and substract from fluxes --
          ZFAC              = ZUPFLXL(JL,JK,JD) / ZPFLXTOT
          ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK,JD) - ZPEVAPUP * ZFAC
          ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK,JD) - ZPEVAPUP * (1 - ZFAC)
          ZUPFLXL(JL,JK,JD) = MAX(0._JPRB,ZUPFLXL(JL,JK,JD))
          ZUPFLXN(JL,JK,JD) = MAX(0._JPRB,ZUPFLXN(JL,JK,JD))
        ENDIF
        
      ENDDO
    ENDDO
    
    !Add contribution to total flux - weight by updraft area fraction
    DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
        PFPLVL(JL,JK) = PFPLVL(JL,JK) + ZUPFLXL(JL,JK,JD)
        PFPLVN(JL,JK) = PFPLVN(JL,JK) + ZUPFLXN(JL,JK,JD)
      ENDDO
    ENDDO
    
  ENDDO !JD
        
        

!     -----------------------------------------------------------------

!*         7.     CONSTRUCT MASS FLUX PROFILES (JD=2,3)
!                 -------------------------------------


  !*         7.1  DETERMINE THE MIXED LAYER SCALING HEIGHT FOR JD=2,3
  !*
  DO JL=KIDIA,KFDIA
    
    SELECT CASE (KPBLTYPE(JL))
    
      CASE(1)
        !Dry convective PBL - no moist parcel
        IZI(JL,2) = KPTOP(JL,2)      !half level below level of zero kinetic energy
                
      CASE(2)
        !Stratocumulus - no dry parcel
        IZI(JL,3) = KPLCL(JL,3)+1      !half level below lcl
          
      CASE(3)
        !Shallow cumulus - both dry and moist
        IZI(JL,3) = KPLCL(JL,3)+1    !half level below lcl
        IZI(JL,2) = KPTOP(JL,2)      !half level below level of zero kinetic energy

      CASE(4)
        !Deep cumulus - no moist parcel
        IZI(JL,2) = KPTOP(JL,2)  

    END SELECT 

  ENDDO !JL 



  !*         7.2  CONSTRUCT MIXED LAYER MASS FLUXES
  !*               - USE CONSTANT AREA FRACTION, AND MULTIPLY BY PARCEL W
  !*
  DO JD = 2,KDRAFT

    DO JK=KLEV-1,1,-1

      DO JL=KIDIA,KFDIA

        IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB) THEN
          IF (JK>=IZI(JL,JD) ) THEN
          
            ZWUH = MAX( ZWU2H(JL,JK,JD),0._JPRB )
            ZWUH = ZWUH**0.5_JPRB
            PMFLX(JL,JK,JD)  = ZFRACB(JL,JD) * ZWUH * ZRHOH(JL,JK)
          
            IF (ZWU2H(JL,JK,JD)>0._JPRB) THEN
              ZFRAC(JL,JK,JD) = ZFRACB(JL,JD)
            ELSE
              ZFRAC(JL,JK,JD) = 0._JPRB
            ENDIF
            
          ENDIF
        ENDIF
        
      ENDDO !JL
    
    ENDDO !JK
    
  ENDDO !JD
  
    
  
  !*         7.3    BUOYANCY SORTING ON MOIST UPDRAFT IN CLOUD LAYER
  !*
  
  CALL VDFBUOYSORT( KIDIA     , KFDIA   , KLON    , KLEV   , KDRAFT , &
                  & PAPM1     , PGEOM1  , PGEOH   , &
                  & ZQTM1     , ZSLGM1  , &
                  & PFRACB    , KPLCL   , KPTOP   , KPLZB  , &
                  & PQTUH     , PSLGUH  , ZWU2H   , PUUH   , PVUH   , &
                ! DIAGNOSTIC OUTPUT
                  & PEXTR2    , KFLDX2  , PEXTRA  , KLEVX  , KFLDX  , &
                !              
                  & ZABULK    , ZWBULK  , ZQTBULK , ZSLGBULK , ZUBULK , ZVBULK )
  


  !*         7.4    CONSTRUCT CLOUDY MASS FLUX PROFILE (JD=3 ONLY)
  !*
  DO JK=KLEV-2,1,-1

    DO JL=KIDIA,KFDIA
   
      IF ( KPBLTYPE(JL)/=0 .AND. KPBLTYPE(JL)/=4 .AND. ZFRACB(JL,3)>0._JPRB ) THEN
        
        IF (JK>=KPTOP(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN

        ZWUH = MAX(0._JPRB,ZWU2H(JL,JK,3))**0.5_JPRB
        
        IF( JK==KPLCL(JL,3) ) THEN 

          !  Special treatment for cloud base level

          ! -- moist updraft --
          PMFLX(JL,JK,3) = MIN( ZFRACB(JL,3)*ZWUH*ZRHOH(JL,JK), 2._JPRB*PMFLX(JL,JK+1,3))   !cap acceleration through cloud base (can cause instability)
          ZFRAC(JL,JK,3) = ZFRACB(JL,3)

          ! -- dry updraft --
          !    tie dry updraft flux to its value at layer below, to prevent too sharp gradients
          PMFLX(JL,JK,2)  = 1.0_JPRB*PMFLX(JL,JK+1,2)
          ZFRAC(JL,JK,2)  = 1.0_JPRB*ZFRAC(JL,JK+1,2)
          PQTUH(JL,JK,2)  = 0.3_JPRB*( PQTUH(JL,JK+1,2)  - ZQTENH(JL,JK+1)  ) + ZQTENH(JL,JK)
          PSLGUH(JL,JK,2) = 0.3_JPRB*( PSLGUH(JL,JK+1,2) - ZSLGENH(JL,JK+1) ) + ZSLGENH(JL,JK)
          PUUH(JL,JK,2)   = 0.3_JPRB*( PUUH(JL,JK+1,2)   - ZUENH(JL,JK+1)   ) + ZUENH(JL,JK)
          PVUH(JL,JK,2)   = 0.3_JPRB*( PVUH(JL,JK+1,2)   - ZVENH(JL,JK+1)   ) + ZVENH(JL,JK)
          
        ELSEIF ( PMFLX(JL,JK+1,3) > 0._JPRB ) THEN
          
          IF ( JK>=KPTOP(JL,3) ) THEN
            
            !-- convective cloud layer --
            !   buoyancy sorting
            
            ZFRAC(JL,JK,3)  = ZABULK(JL,JK)
            ZWU2H(JL,JK,3)  = ZWBULK(JL,JK) ** 2._JPRB
            PMFLX(JL,JK,3)  = ZABULK(JL,JK) * ZWBULK(JL,JK) * ZRHOH(JL,JK)
            PQTUH(JL,JK,3)  = ZQTBULK(JL,JK)
            PSLGUH(JL,JK,3) = ZSLGBULK(JL,JK)
            PUUH(JL,JK,3)   = ZUBULK(JL,JK)
            PVUH(JL,JK,3)   = ZVBULK(JL,JK)
            
          ELSE 
            
            !-- inversion (layer between tops of test parcel and moist parcel) --
            !   prescribed linear decay of flux
            
           !ZFRAC(JL,JK,3)  = ZFRAC(JL,JK+1,3) * 0.25_JPRB   !fast "exponential" top inversion
           !PMFLX(JL,JK,3)  = PMFLX(JL,JK+1,3) * 0.25_JPRB
            zfac = ( PZPTOP(JL,1) - PGEOH(JL,JK)*ZRG  ) / ( PZPTOP(JL,1) - PZPTOP(JL,3) )
            zfac = MAX( 0._JPRB, MIN(1._JPRB,zfac) )
            ZFRAC(JL,JK,3)  = ZFRAC(JL,KPTOP(JL,3),3) * zfac
            PMFLX(JL,JK,3)  = PMFLX(JL,KPTOP(JL,3),3) * zfac
            
            PQTUH(JL,JK,3)  = ( PQTUH(JL,JK+1,3)  - ZQTENH(JL,JK+1)  ) + ZQTENH(JL,JK)
            PSLGUH(JL,JK,3) = ( PSLGUH(JL,JK+1,3) - ZSLGENH(JL,JK+1) ) + ZSLGENH(JL,JK)
            PUUH(JL,JK,3)   = ( PUUH(JL,JK+1,3)   - ZUENH(JL,JK+1)   ) + ZUENH(JL,JK)
            PVUH(JL,JK,3)   = ( PVUH(JL,JK+1,3)   - ZVENH(JL,JK+1)   ) + ZVENH(JL,JK)
            
          ENDIF
          
          !make sure that updraft #2 does not do any flux here
          PMFLX(JL,JK,2)  = 0._JPRB
          ZFRAC(JL,JK,2)  = 0._JPRB
          PQTUH(JL,JK,2)  = 0._JPRB
          PSLGUH(JL,JK,2) = 0._JPRB
          PUUH(JL,JK,2)   = 0._JPRB
          PVUH(JL,JK,2)   = 0._JPRB
          
        ENDIF
        
        
        !--- limit mass flux covering 50% area (M<rho*w,up*0.5) ---
        !    (detrainment is initiated if strong w,up slowdown)

          PMFLX(JL,JK,3) =   MIN( PMFLX(JL,JK,3) , 0.5_JPRB * ZWUH * ZRHOH(JL,JK) )   

        ENDIF
        
      ENDIF
     
    ENDDO !JL
     
  ENDDO !JK


  !*        7.5  COMPUTE I) VARIANCE TRANSPORT FLUX AND II) BULK UPDRAFT DETRAINMENT, AS USED IN VDFMAIN
  !*
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA
      IF (JK >= KPTOP(JL,1) ) THEN   !highest possible level for any parcel (use test parcel) 
!     IF (JK >  KPTOP(JL,1) ) THEN  
      
        !-- calculate w'qt'qt' at half levels --
!       PWQT2(JL,JK) =  PMFLX(JL,JK,3) * (PQTUH(JL,JK,3) - ZQTENH(JL,JK))**2._JPRB + &
!                    &  PMFLX(JL,JK,2) * (PQTUH(JL,JK,2) - ZQTENH(JL,JK))**2._JPRB
        !-- upwind discretization as in solver
        PWQT2(JL,JK) =  PMFLX(JL,JK,3) * (PQTUH(JL,JK,3) - ZQTM1(JL,JK))**2._JPRB + &
                     &  PMFLX(JL,JK,2) * (PQTUH(JL,JK,2) - ZQTM1(JL,JK))**2._JPRB

      ENDIF 
    ENDDO !JL
  ENDDO !JK

  
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF ( KPBLTYPE(JL)/=0 .AND. KPBLTYPE(JL)/=4 .AND. ZFRACB(JL,3)>0._JPRB ) THEN
        
      IF (JK >= KPTOP(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN  
!     IF (JK >  KPTOP(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN  
        
        ZDZ   =  (   PGEOH(JL,JK-1)   - PGEOH(JL,JK)   ) * ZRG 
        
        ZDMDZ = -(   PMFLX(JL,JK-1,3) - PMFLX(JL,JK,3) &
!                 & + PMFLX(JL,JK-1,2) - PMFLX(JL,JK,2) &  
                 & ) / ZDZ
        
        IF (JK==KPTOP(JL,1)+1) THEN
          ZDMDZ = PMFLX(JL,JK,3) / ZDZ
        ENDIF
        
        ZDMDZ = MAX( 0._JPRB , ZDMDZ )
                        
        PDETR(JL,JK) = ZDMDZ + PMFLX(JL,JK,3) /  &
          &            ( ZTAUEPS * MAX(0.0001_JPRB,ZWU2H(JL,JK,3))**0.5_JPRB )  
        
      ENDIF
              
      ENDIF
    ENDDO !JL
  ENDDO !JK
          


  !*        7.6  Mass flux limit according to CFL criterion
  !*              (preserving the vertical structure of M profile)
  !* 
  ZCONS10 = 1.0_JPRB/(RG*PTMST)

  DO JD = 2,KDRAFT
    DO JL=KIDIA,KFDIA
      ZMFS(JL,JD) = 1.0_JPRB  ! mass flux scaling value (reduction)
    ENDDO
  ENDDO
  
  DO JD = 2,KDRAFT
    DO JK=1,KLEV-1
      DO JL=KIDIA,KFDIA
        IF ( JK >= KPTOP(JL,JD) .AND. KPTOP(JL,JD)>0) THEN
          ZMFMAX = (PAPM1(JL,JK+1)-PAPM1(JL,JK)) * ZCONS10
!test! 4xCFL          
!         ZMFMAX = 4.0_JPRB * ZMFMAX
!         ZMFMAX = 2.0_JPRB * ZMFMAX
          
          IF ( PMFLX(JL,JK,JD) > ZMFMAX ) THEN
            ZMFS(JL,JD) = MIN(ZMFS(JL,JD),ZMFMAX/PMFLX(JL,JK,JD))
          ENDIF
          
!          PMFLX(JL,JK,JD) = MIN( ZMFMAX, PMFLX(JL,JK,JD) )  !correct level by level
          
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  
  DO JD = 2,KDRAFT
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PMFLX(JL,JK,JD) = PMFLX(JL,JK,JD)*ZMFS(JL,JD)
      ENDDO
    ENDDO
  ENDDO



  !*        7.7  CLOUD-LAYER AVERAGE MOIST UPDRAFT MASS FLUX.
  !*               (USED IN VDFEXCU FOR ESTIMATING ENTRAINMENT-K AT CUMULUS PBL TOP)  
  !*
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA
      !-- cloudy mass flux: use moist updraft M (depth average) --
      IF ( ZQCUH(JL,JK,3)>0._JPRB .AND. PMFLX(JL,JK,3)>0._JPRB ) THEN  
        ZMCLD(JL) = ZMCLD(JL) + PMFLX(JL,JK,3) * &
                &     ZRG*( PGEOH(JL,JK-1) - PGEOH(JL,JK) )
      ENDIF
    ENDDO
  ENDDO
  
  DO JL=KIDIA,KFDIA
    IF (PZPTOP(JL,3)-PZPLCL(JL,3)>0._JPRB ) THEN  
!      PMCU(JL)    = ZMCLD(JL) / (PZPTOP(JL,3)-PZPLCL(JL,3))
      PMCU(JL) = ZMFLXB(JL,3)
    ENDIF
  ENDDO


  IF ( LLDIAG ) THEN
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PEXTRA(JL,JK,40) = PEXTRA(JL,JK,40) + PMFLX(JL,JK,2)*PTMST
        PEXTRA(JL,JK,48) = PEXTRA(JL,JK,48) + PMFLX(JL,JK,3)*PTMST
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,22,41) = PMFLX(JL,KPLCL(JL,3),3)
      PEXTRA(JL,27,41) = RCPD*ZKHVFL(JL) &
                     & * PAPHM1(JL,KLEV-1)/RD/PTM1(JL,KLEV)          !rho
    ENDDO
  ENDIF


    
  !*        7.8  SCALING (SEE VDFEXCU)
  !* 
  DO JD = 2,KDRAFT
    DO JK=1,KLEV-1
      DO JL=KIDIA,KFDIA
      
        ZMGEOM(JL,JK)  = PGEOM1(JL,JK)-PGEOM1(JL,JK+1)      
        ZCFNC1(JL,JK)  = RVDIFTS * PTMST * RG**2 * PAPHM1(JL,JK) &
                 & /( ZMGEOM(JL,JK) * RD * 0.5_JPRB &
                 & *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  )) &
                 & +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1))))  
        PMFLX(JL,JK,JD) = ZCFNC1(JL,JK) * ZMGEOM(JL,JK) * ZRG * PMFLX(JL,JK,JD)
        
      ENDDO
    ENDDO
  ENDDO



!     -----------------------------------------------------------------

!*         8.     W-SCALE USED IN CLOUD VARIANCE DISSIPATION (VDFMAIN)
!                 ----------------------------------------------------

  DO JL=KIDIA,KFDIA
    !PWUAVG(JL) = ZWSTAR(JL)
    PWUAVG(JL) = 2.0_JPRB * ZWSTAR(JL)
  ENDDO

  !for use in buoyancy sorting scheme
  DO JD = 1,KDRAFT
    DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
        PWUH(JL,JK,JD) = SQRT( ABS( ZWU2H(JL,JK,JD) ) )
      ENDDO
    ENDDO
  ENDDO



!     -----------------------------------------------------------------

!*         9.     ADVECTIVE FLUX ADJUSTMENTS AT UPDRAFT TOP-LEVEL
!                 -----------------------------------------------
!
!                 Vertical mixing at the top-level is prescribed
!                 and controlled in VDFEXCU.
!

  DO JD = 2,KDRAFT
    DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
        
        !  Remove all mass flux at and above top layer of updraft
        ITOP = MAX( KPTOP(JL,JD), KPTOP(JL,3) )
        
        IF ( JK <= ITOP ) THEN
          
          LLWIPE=.TRUE.    
          
          ! Protect cumulus cloud top: this is treated later in VDFEXCU, dept. on LLRICU=T
          !   Note: do wipe in case of single layer moist convection
          IF ( KPBLTYPE(JL)==3 .AND. KPTOP(JL,3)<KPLCL(JL,3)+1 ) THEN   
            LLWIPE=.FALSE. 
          ENDIF
          
          IF (LLWIPE) THEN
            ZFRAC(JL,JK,JD)  = 0.0_JPRB
            PMFLX(JL,JK,JD)  = 0.0_JPRB
            PWQT2(JL,JK)     = 0.0_JPRB
          ENDIF
          
        ENDIF
        
      ENDDO
    ENDDO
  ENDDO


IF (LHOOK) CALL DR_HOOK('VDFHGHTN',1,ZHOOK_HANDLE)

END SUBROUTINE VDFHGHTN


END MODULE mo_vdfhghtn
