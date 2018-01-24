!>
!! VDF sub-stepping for EDMF DUALM
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2011-9-30)
!!   (IFS cycle CY36R1_DUALM_M8b)
!!
!! Modifications by Dmitrii Mironov, DWD (2016-08-08)
!! - Changes related to the use of prognostic the sea-ice albedo.
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

MODULE mo_vdfouter

  PUBLIC :: vdfouter

CONTAINS

SUBROUTINE VDFOUTER   ( &
 & KIDIA  , KFDIA  , KLON   , KLEV   , KLEVS  , KSTEP  , KTILES , &
 & KTRAC  , KLEVSN , KLEVI  ,  &
 & PTSPHY , PSIGFLT, &
 & PUM1   , PVM1   , PTM1   , PQM1   , PLM1   , PIM1   , PAM1   , &
 & PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , PTSKM1M, &
 & PSSRFL , PSLRFL , PEMIS  , PHRLW  , PHRSW  , &
 & PTSNOW , &
 & PCHAR  , PUCURR , PVCURR , PTSKRAD, &
 & PSOTEU , PSOTEV , PSOBETA, PVERVEL, &
 ! OUTPUT
 & PZ0M   , PZ0H   , &
 & PVAR   , &
 & KHPBLN , KVARTOP, &
 & KPBLTYPE,&
 & PFPLVL , PFPLVN , &
 ! DIAGNOSTIC OUTPUT
 & KFLDX2 , KLEVX  , KFLDX  , LLDIAG,&
 ! OUTPUT TENDENCIES
 & PTE    , PQE    , PLE    , PIE    , PAE    , PVOM   , PVOL   , &
 ! UPDATED FIELDS FOR TILES
 & PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI , &
 ! OUTPUT FLUXES
 & PDIFTS , PDIFTQ , PDIFTL , PDIFTI , PSTRTU , PSTRTV , &
 & PKH    , PKM    , &
 & LDLAND , &
! surface fluxes from TERRA (tq) & TURBTRAN (uv) for EDMF atmospheric transport
 & SHFL_S , LHFL_S , UMFL_S , VMFL_S , &
! TERRA data
 & frac_t , t_g_t  ,qv_s   , t_ice )


!***

!**   *VDFMAIN* - DOES THE VERTICAL EXCHANGE OF U,V,SLG,QT BY TURBULENCE.

!     J.F.GELEYN       20/04/82   Original  
!     C.A.BLONDIN      18/12/86
!     A.C.M. BELJAARS  20/10/89   IFS-VERSION (TECHNICAL REVISION OF CY34)
!     A.C.M. BELJAARS  26/03/90   OBUKHOV-L UPDATE 
!     A.C.M. BELJAARS  30/09/98   SURFACE TILES 
!     P. Viterbo       17/05/2000 Surface DDH for TILES
!     D. Salmond       15/10/2001 FULLIMP mods
!     S. Abdalla       27/11/2001 Passing Zi/L to waves
!     J.Hague          25/06/2003 Tuning for p690
!     P.Viterbo        24/05/2004 Change surface units
!     M. Ko"hler        3/12/2004 Moist Advection-Diffusion
!     A.Beljaars (12-11-02) Introduction of ocean current b.c.
!     A. Beljaars      30/09/2005 Include Subgr. Oro. in solver
!     G. Balsamo       15/01/2007 Include soil type
!     A. Beljaars      27/02/2009 Delete PZIDLWV
!     G. Balsamo       18/04/2008 Include lake tile

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE FOUR
!     PROGNOSTIC VARIABLES U,V,T AND Q DUE TO THE VERTICAL EXCHANGE BY
!     TURBULENT (= NON-MOIST CONVECTIVE) PROCESSES. THESE TENDENCIES ARE
!     OBTAINED AS THE DIFFERENCE BETWEEN THE RESULT OF AN IMPLICIT
!     TIME-STEP STARTING FROM VALUES AT T-1 AND THESE T-1 VALUES. ALL
!     THE DIAGNOSTIC COMPUTATIONS (EXCHANGE COEFFICIENTS, ...) ARE DONE
!      FROM THE T-1 VALUES. AS A BY-PRODUCT THE ROUGHNESS LENGTH OVER SEA
!     IS UPDATED ACCORDINGLY TO THE *CHARNOCK FORMULA. HEAT AND MOISTURE
!     SURFACE FLUXES AND THEIR DERIVATIVES AGAINST TS, WS AND WL
!     (THE LATTER WILL BE LATER WEIGHTED WITH THE SNOW FACTOR IN
!     *VDIFF*), LATER TO BE USED FOR SOIL PROCESSES TREATMENT, ARE ALSO
!     COMPUTED AS WELL AS A STABILITY VALUE TO BE USED AS A DIAGNOSTIC
!     OF THE DEPTH OF THE WELL MIXED LAYER IN CONVECTIVE COMPUTATIONS.

!     INTERFACE.
!     ----------
!          *VDIFF* TAKES THE MODEL VARIABLES AT T-1 AND RETURNS THE VALUES
!     FOR THE PROGNOSTIC TIME T+1 DUE TO VERTICAL DIFFUSION.
!     THE MODEL VARIABLES, THE MODEL DIMENSIONS AND THE DIAGNOSTICS DATA
!     ARE PASSED AS SUBROUTINE ARGUMENTS. CONSTANTS THAT DO NOT CHANGE
!     DURING A MODEL RUN (E.G. PHYSICAL CONSTANTS, SWITCHES ETC.) ARE
!     STORED IN A SINGLE COMMON BLOCK *YOMVDF*, WHICH IS INITIALIZED
!     BY SET-UP ROUTINE *SUVDF*.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*        NUMBER OF SOIL LAYERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                   OF SURFACE BOUNDARY CONDITION)
!    *KTRAC*        Number of tracers
!    *KLEVSN*       Number of snow layers (diagnostics) 
!    *KLEVI*        Number of sea ice layers (diagnostics)
!!!  *KDHVTLS*      Number of variables for individual tiles
!!!  *KDHFTLS*      Number of fluxes for individual tiles
!!!  *KDHVTSS*      Number of variables for snow energy budget
!!!  *KDHFTSS*      Number of fluxes for snow energy budget
!!!  *KDHVTTS*      Number of variables for soil energy budget
!!!  *KDHFTTS*      Number of fluxes for soil energy budget
!!!  *KDHVTIS*      Number of variables for sea ice energy budget
!!!  *KDHFTIS*      Number of fluxes for sea ice energy budget

!!!  *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!!!  *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!!!  *KSOTY*        SOIL TYPE                                    (1-7)

!!   *KCNT*         Index of vdf sub steps.

!     INPUT PARAMETERS (LOGICAL)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):

!    *PCVL*         LOW VEGETATION COVER                          -
!    *PCVH*         HIGH VEGETATION COVER                         -
!    *PSIGFLT*      STANDARD DEVIATION OF FILTERED OROGRAPHY      M
!    *PUM1*         X-VELOCITY COMPONENT                          M/S
!    *PVM1*         Y-VELOCITY COMPONENT                          M/S
!    *PTM1*         TEMPERATURE                                   K
!    *PQM1*         SPECIFIC HUMIDITY                             KG/KG
!    *PLM1*         SPECIFIC CLOUD LIQUID WATER                   KG/KG
!    *PIM1*         SPECIFIC CLOUD ICE                            KG/KG
!    *PAM1*         CLOUD FRACTION                                1
!    *PCM1*         TRACER CONCENTRATION                          KG/KG
!    *PAPM1*        PRESSURE ON FULL LEVELS                       PA
!    *PAPHM1*       PRESSURE ON HALF LEVELS                       PA
!    *PGEOM1*       GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL AT HALF LEVELS                   M2/S2
!    *PTSKM1M*      SKIN TEMPERATURE                              K
!    *PTSAM1M*      SURFACE TEMPERATURE                           K
!    *PWSAM1M*      SOIL MOISTURE ALL LAYERS                      M**3/M**3
!    *PSSRFL*       NET SHORTWAVE RADIATION FLUX AT SURFACE       W/M2
!    *PSLRFL*       NET LONGWAVE RADIATION FLUX AT SURFACE        W/M2
!    *PEMIS*        MODEL SURFACE LONGWAVE EMISSIVITY
!    *PHRLW*        LONGWAVE HEATING RATE                         K/s
!    *PHRSW*        SHORTWAVE HEATING RATE                        K/s
!    *PTSNOW*       SNOW TEMPERATURE                              K
!!   *PTICE*        ICE TEMPERATURE (TOP SLAB)                    K
!    *PTLICE*       LAKE ICE TEMPERATURE                          K
!    *PHLICE*       LAKE ICE THICKNESS                            m
!    *PTLWML*       LAKE MEAN WATER TEMPERATURE                   K
!!!  *PSST*         (OPEN) SEA SURFACE TEMPERATURE                K
!!!  *PFRTI*        TILE FRACTIONS                                (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PALBTI*       BROADBAND ALBEDO FOR TILE FRACTIONS
!    *PWLMX*        MAXIMUM SKIN RESERVOIR CAPACITY               kg/m**2
!    *PCHAR*        "EQUIVALENT" CHARNOCK PARAMETER
!    *PUCURR*       OCEAN CURRENT X-DIRECTION
!    *PVCURR*       OCEAN CURRENT Y-DIRECTION
!    *PTSKRAD*      SKIN TEMPERATURE OF LATEST FULL RADIATION
!                      TIMESTEP                                   K
!    *PCFLX*        TRACER SURFACE FLUX                           kg/(m2 s)
!    *PSOTEU*       Explicit part of U-tendency from subgrid orography scheme    
!    *PSOTEV*       Explicit part of V-tendency from subgrid orography scheme     
!    *PSOBETA*      Implicit part of subgrid orography 

!    *PVERVEL*      VERTICAL VELOCITY

!     INPUT PARAMETERS (LOGICAL):

!     CONTRIBUTIONS TO BUDGETS (OUTPUT,REAL):

!    *PVDIS*        TURBULENT DISSIPATION                         W/M2
!    *PVDISG*       SUBGRID OROGRAPHY DISSIPATION                 W/M2
!    *PAHFLEV*      LATENT HEAT FLUX  (SNOW/ICE FREE PART)        W/M2
!    *PAHFLSB*      LATENT HEAT FLUX  (SNOW/ICE COVERED PART)     W/M2

!     UPDATED PARAMETERS (REAL):

!    *PTE*          TEMPERATURE TENDENCY                          K/S
!    *PQE*          MOISTURE TENDENCY                             KG/(KG S)
!    *PLE*          LIQUID WATER TENDENCY                         KG/(KG S)
!    *PIE*          ICE WATER TENDENCY                            KG/(KG S)
!    *PAE*          CLOUD FRACTION TENDENCY                       1/S)
!    *PVOM*         MERIODINAL VELOCITY TENDENCY (DU/DT)          M/S2
!    *PVOL*         LATITUDE TENDENCY            (DV/DT)          M/S2
!    *PTENC*        TRACER TENDENCY                               KG/(KG S)
!!   *PTSKE1*       SKIN TEMPERATURE TENDENCY                     K/S
!    *PZ0M*         AERODYNAMIC ROUGHNESS LENGTH                  M
!    *PZ0H*         ROUGHNESS LENGTH FOR HEAT                     M

!     UPDATED PARAMETERS FOR TILES (REAL): 

!    *PUSTRTI*      SURFACE U-STRESS                              N/M2 
!    *PVSTRTI*      SURFACE V-STRESS                              N/M2
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2
!    *PEVAPTI*      SURFACE MOISTURE FLUX                         KG/M2/S
!    *PTSKTI*       SKIN TEMPERATURE                              K

!     OUTPUT PARAMETERS (REAL):

!    *PFPLVL*       PBL PRECIPITATION FLUX AS RAIN                KG/(M**2*S)
!    *PFPLVN*       PBL PRECIPITATION FLUX AS SNOW                KG/(M**2*S)
!    *PFHPVL*       ENTHALPY FLUX OF PBL PRECIPITATION AS RAIN    J/(M**2*S)
!    *PFHPVN*       ENTHALPY FLUX OF PBL PRECIPITATION AS SNOW    J/(M**2*S)

!    *PLDIFF*       CONTRIB TO PBL CONDENSATE BY PASSIVE CLOUDS   KG/KG

!    *PFWSB*        EVAPORATION OF SNOW                           KG/(M**2*S)
!    *PU10M*        U-COMPONENT WIND AT 10 M                      M/S
!    *PV10M*        V-COMPONENT WIND AT 10 M                      M/S
!    *PT2M*         TEMPERATURE AT 2M                             K
!    *PD2M*         DEW POINT TEMPERATURE AT 2M                   K
!    *PQ2M*         SPECIFIC HUMIDITY AT 2M                       KG/KG
!    *PGUST*        GUST AT 10 M                                  M/S
!    *PBLH*         PBL HEIGHT (dry diagnostic based on Ri#)      M
!    *PZINV*        PBL HEIGHT (moist parcel, not for stable PBL) M
!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                      EACH TILE                                  W/M2
!    *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST            KG/(M2*S)
!    *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PTOFDU*       TOFD COMP. OF TURBULENT FLUX OF U-MOMEMTUM   KG*(M/S)/(M2*S)
!    *PTOFDV*       TOFD COMP. OF TURBULENT FLUX OF V-MOMEMTUM   KG*(M/S)/(M2*S)
!    *PDIFTS*       TURBULENT FLUX OF HEAT                         J/(M2*S)
!    *PDIFTQ*       TURBULENT FLUX OF SPECIFIC HUMIDITY           KG/(M2*S)
!    *PDIFTL*       TURBULENT FLUX OF LIQUID WATER                KG/(M2*S)
!    *PDIFTI*       TURBULENT FLUX OF ICE WATER                   KG/(M2*S)
!    *PSTRSOU*      SUBGRID OROGRAPHY FLUX OF U-MOMEMTUM    KG*(M/S)/(M2*S)
!    *PSTRSOV*      SUBGRID OROGRAPHY FLUX OF V-MOMEMTUM    KG*(M/S)/(M2*S)

!    *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!                   IN SURFACE LAYER: CH*U                        (M/S)
!    *PKM*          TURB. DIFF. COEFF. FOR MOM. ABOVE SURF. LAY.  (M2/S)
!                   IN SURFACE LAYER: CH*U                        (M/S)
!!!  *PDHTLS*       Diagnostic array for tiles (see module yomcdh)
!!!                    (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!!!  *PDHTSS*       Diagnostic array for snow T (see module yomcdh)
!!!                    (Wm-2 for fluxes)
!!!  *PDHTTS*       Diagnostic array for soil T (see module yomcdh)
!!!                    (Wm-2 for fluxes)
!!!  *PDHTIS*       Diagnostic array for ice T (see module yomcdh)
!                      (Wm-2 for fluxes)

!     METHOD.
!     -------

!          FIRST AN AUXIALIARY VARIABLE CP(Q)T+GZ IS CREATED ON WHICH
!     THE VERTICAL DIFFUSION PROCESS WILL WORK LIKE ON U,V AND Q. THEN
!     ALONG THE VERTICAL AND AT THE SURFACE, EXCHANGE COEFFICIENTS (WITH
!     THE DIMENSION OF A PRESSURE THICKNESS) ARE COMPUTED FOR MOMENTUM
!     AND FOR HEAT (SENSIBLE PLUS LATENT). THE LETTERS M AND H ARE USED
!     TO DISTINGUISH THEM AND THE COMPUTATION IS THE RESULT OF A
!     CONDITIONAL MERGE BETWEEN THE STABLE AND THE UNSTABLE CASE
!     (DEPENDING ON THE SIGN OF THE *RICHARDSON BULK NUMBER).
!          IN THE SECOND PART OF THE ROUTINE THE IMPLICIT LINEAR
!     SYSTEMS FOR U,V FIRST AND T,Q SECOND ARE SOLVED BY A *GAUSSIAN
!     ELIMINATION BACK-SUBSTITUTION METHOD. FOR T AND Q THE LOWER
!     BOUNDARY CONDITION DEPENDS ON THE SURFACE STATE.
!     OVER LAND, TWO DIFFERENT REGIMES OF EVAPORATION PREVAIL:
!     A STOMATAL RESISTANCE DEPENDENT ONE OVER THE VEGETATED PART
!     AND A SOIL RELATIVE HUMIDITY DEPENDENT ONE OVER THE
!     BARE SOIL PART OF THE GRID MESH.
!     POTENTIAL EVAPORATION TAKES PLACE OVER THE SEA, THE SNOW
!     COVERED PART AND THE LIQUID WATER COVERED PART OF THE
!     GRID MESH AS WELL AS IN CASE OF DEW DEPOSITION.
!          FINALLY ONE RETURNS TO THE VARIABLE TEMPERATURE TO COMPUTE
!     ITS TENDENCY AND THE LATER IS MODIFIED BY THE DISSIPATION'S EFFECT
!     (ONE ASSUMES NO STORAGE IN THE TURBULENT KINETIC ENERGY RANGE) AND
!     THE EFFECT OF MOISTURE DIFFUSION ON CP. Z0 IS UPDATED AND THE
!     SURFACE FLUXES OF T AND Q AND THEIR DERIVATIVES ARE PREPARED AND
!     STORED LIKE THE DIFFERENCE BETWEEN THE IMPLICITELY OBTAINED
!     CP(Q)T+GZ AND CP(Q)T AT THE SURFACE.

!     EXTERNALS.
!     ----------

!     *VDFMAIN* CALLS SUCESSIVELY:
!         *VDFSURF*
!         *VDFEXCS*
!         *VDFEVAP*
!         *VDFEXCU*
!         *VDFDIFM*
!         *VDFDIFH*
!         *VDFDIFC*
!         *VDFINCR*
!         *VDFSDRV*
!         *VDFPPCFL*
!         *VDFUPDZ0*

!     REFERENCE.
!     ----------

!          SEE VERTICAL DIFFUSION'S PART OF THE MODEL'S DOCUMENTATION
!     FOR DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
! USE YOEPHY    ,ONLY : LVDFTRAC

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook
USE mo_edmf_param   ,ONLY : LVDFTRAC                              !yoephy

USE mo_vdfmain      ,ONLY : vdfmain 

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PCVL(:)                   !(KLON) 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PCVH(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(:)                !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(:,:)                 !(KLON,KLEV)
! REAL(KIND=JPRB) ,INTENT(IN)    :: PCM1(:,:,:)               !(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(:,:)                !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(:,:)               !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(:,0:)               !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)                !(KLON) 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PTSAM1M(:,:)              !(KLON,KLEVS) 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PWSAM1M(:,:)              !(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(:,:)                !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(:,:)                !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(:)                 !(KLON) 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PHLICE(:)                 !(KLON)
! REAL(KIND=JPRB) ,INTENT(IN)    :: PTLICE(:)                 !(KLON)
! REAL(KIND=JPRB) ,INTENT(IN)    :: PTLWML(:)                 !(KLON)
! REAL(KIND=JPRB) ,INTENT(IN)    :: PALBTI(:,:)               !(KLON,KTILES) 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PWLMX(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(:)                !(KLON) 
! REAL(KIND=JPRB) ,INTENT(IN)    :: PCFLX(:,:)                !(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(:,:)               !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(:,:)               !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOBETA(:,:)              !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(:)                   !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PVDIS(:)                  !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PVDISG(:)                 !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PDISGW3D(:,:)             !(KLON,KLEV) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PAHFLEV(:)                !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PAHFLSB(:)                !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PFWSB(:)                  !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PBIR(:)                   !(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVAR(:,:)                 !(KLON,KLEV)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PU10M(:)                  !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PV10M(:)                  !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PT2M(:)                   !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PD2M(:)                   !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PQ2M(:)                   !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PZINV(:)                  !(KLON)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PBLH(:)                   !(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KHPBLN(:)                 !(KLON) 
! REAL(KIND=JPRB) ,INTENT(INOUT) :: PSSRFLTI(:,:)             !(KLON,KTILES) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PEVAPSNW(:)               !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVL(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVN(:,0:)              !(KLON,0:KLEV)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PFHPVL(:,0:)              !(KLON,0:KLEV)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PFHPVN(:,0:)              !(KLON,0:KLEV)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PGUST(:)                  !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PWUAVG(:)                 !(KLON) 
! LOGICAL         ,INTENT(OUT)   :: LDNODECP(:)               !(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(:)               !(KLON)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PLDIFF(:,:)               !(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KVARTOP(:)                !(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PIE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(:,:)                 !(KLON,KLEV) 
! REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENC(:,:,:)              !(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:,:)              !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:,:)              !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:)              !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:)              !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(:,:)               !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTS(:,0:)              !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQ(:,0:)              !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTL(:,0:)              !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTI(:,0:)              !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(:,0:)              !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(:,0:)              !(KLON,0:KLEV) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PTOFDU(:)                 !(KLON) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PTOFDV(:)                 !(KLON)
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PSTRSOU(:,0:)             !(KLON,0:KLEV) 
! REAL(KIND=JPRB) ,INTENT(OUT)   :: PSTRSOV(:,0:)             !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKM(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(:,:)              !(KLON,KLEV) 
LOGICAL                          :: LDLAND(:)                 !(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SHFL_S(:)                 !(KLON) sensible heat flux from TERRA
REAL(KIND=JPRB)   ,INTENT(IN)    :: LHFL_S(:)                 !(KLON) latent   heat from from TERRA
REAL(KIND=JPRB)   ,INTENT(IN)    :: UMFL_S(:)                 !(KLON) u-flux from TURBTRAN & SFCinterface
REAL(KIND=JPRB)   ,INTENT(IN)    :: VMFL_S(:)                 !(KLON) v-flux from TURBTRAN & SFCinterface
!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX2, KLEVX, KFLDX
! REAL(KIND=JPRB) ,INTENT(INOUT) :: PEXTR2(:,:)               !(KLON,KFLDX2)
! REAL(KIND=JPRB) ,INTENT(INOUT) :: PEXTRA(:,:,:)             !(KLON,KLEVX,KFLDX)
LOGICAL           ,INTENT(IN)    :: LLDIAG

REAL(KIND=JPRB)   ,INTENT(IN)    :: frac_t (:,:)              !(KLON,ntiles_total)  tile fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: t_g_t  (:,:)              !(KLON,ntiles_total)  tiled sfc temperature
REAL(KIND=JPRB)   ,INTENT(IN)    :: qv_s   (:)                !(KLON)               mean  sfc qsat
REAL(KIND=JPRB)   ,INTENT(IN)    :: t_ice  (:)                !(KLON)               ice temperature

!--- dummy variables taken from argument list that are not used

REAL(KIND=JPRB)                  :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)                  :: PTICE(KLON) 
REAL(KIND=JPRB)                  :: PTSKE1(KLON) 
REAL(KIND=JPRB)                  :: PSST(KLON) 
INTEGER(KIND=JPIM)               :: KCNT 
INTEGER(KIND=JPIM)               :: KSOTY(KLON) 
INTEGER(KIND=JPIM)               :: KTVL(KLON) 
INTEGER(KIND=JPIM)               :: KTVH(KLON) 
INTEGER(KIND=JPIM)               :: KDHVTLS 
INTEGER(KIND=JPIM)               :: KDHFTLS 
INTEGER(KIND=JPIM)               :: KDHVTSS 
INTEGER(KIND=JPIM)               :: KDHFTSS 
INTEGER(KIND=JPIM)               :: KDHVTTS 
INTEGER(KIND=JPIM)               :: KDHFTTS 
INTEGER(KIND=JPIM)               :: KDHVTIS 
INTEGER(KIND=JPIM)               :: KDHFTIS 
CHARACTER(LEN=1)                 :: CDCONF 


!*         0.2    LOCAL VARIABLES

! added local variables of deleted arguments that are NOT used

REAL(KIND=JPRB) :: PHLICE  (KLON)          ! lake ice thickness   
REAL(KIND=JPRB) :: PTLICE  (KLON)          ! lake ice temperature
REAL(KIND=JPRB) :: PTLWML  (KLON)          ! lake mean water T   
REAL(KIND=JPRB) :: PCVL    (KLON)          ! low vegetation cover
REAL(KIND=JPRB) :: PCVH    (KLON)          ! high vegetation cover
REAL(KIND=JPRB) :: PWLMX   (KLON)          ! maximum skin reservoir capacity
REAL(KIND=JPRB) :: PVDIS   (KLON)          ! turbulent dissipation
REAL(KIND=JPRB) :: PVDISG  (KLON)          ! SO dissipation
REAL(KIND=JPRB) :: PAHFLEV (KLON)          ! latent heat flux (snow/ice free part)
REAL(KIND=JPRB) :: PAHFLSB (KLON)          ! latent heat flux (snow/ice covered part)
REAL(KIND=JPRB) :: PFWSB   (KLON)          ! evaporation of snow
REAL(KIND=JPRB) :: PBIR    (KLON)          ! BIR buoyancy flux integral ratio
REAL(KIND=JPRB) :: PZINV   (KLON)          ! PBL HEIGHT (moist parcel, not for stable PBL)
REAL(KIND=JPRB) :: PBLH    (KLON)          ! PBL HEIGHT (dry diagnostic based on Ri#)
REAL(KIND=JPRB) :: PEVAPSNW(KLON)          ! evaporation from snow under forest
REAL(KIND=JPRB) :: PGUST   (KLON)          ! 10m gust
REAL(KIND=JPRB) :: PWUAVG  (KLON)          ! w,up averaged
REAL(KIND=JPRB) :: PTOFDU  (KLON)          ! TOFD U flux
REAL(KIND=JPRB) :: PTOFDV  (KLON)          ! TOFD V flux
REAL(KIND=JPRB) :: PU10M   (KLON)          ! 10m U 
REAL(KIND=JPRB) :: PV10M   (KLON)          ! 10m V
REAL(KIND=JPRB) :: PT2M    (KLON)          ! 2m T
REAL(KIND=JPRB) :: PD2M    (KLON)          ! 2m TD
REAL(KIND=JPRB) :: PQ2M    (KLON)          ! 2m Q
LOGICAL         :: LDNODECP(KLON)          ! no decoupling allowed
REAL(KIND=JPRB) :: PDISGW3D(KLON,KLEV)     ! 3D stoch. phys. dissipation
REAL(KIND=JPRB) :: PLDIFF  (KLON,KLEV)     ! contrib to PBL cond. by passive clouds
REAL(KIND=JPRB) :: PFHPVL  (KLON,0:KLEV)   ! PBL rain enthalpy flux
REAL(KIND=JPRB) :: PFHPVN  (KLON,0:KLEV)   ! PBL snow enthalpy flux
REAL(KIND=JPRB) :: PSTRSOU (KLON,0:KLEV)   ! SSO U flux
REAL(KIND=JPRB) :: PSTRSOV (KLON,0:KLEV)   ! SSO V flux
REAL(KIND=JPRB) :: PTSAM1M (KLON,KLEVS)    ! T,soil      
REAL(KIND=JPRB) :: PWSAM1M (KLON,KLEVS)    ! Q,soil
REAL(KIND=JPRB) :: PALBTI  (KLON,KTILES)   ! tile albedo
REAL(KIND=JPRB) :: PSSRFLTI(KLON,KTILES)   ! net SW sfc flux for each tile (use tile ablbedo)
REAL(KIND=JPRB) :: PCM1    (KLON,KLEV,KTRAC) ! tracer - for VDF transport
REAL(KIND=JPRB) :: PTENC   (KLON,KLEV,KTRAC) ! tracer tendency
REAL(KIND=JPRB) :: PCFLX   (KLON,KTRAC)      ! surface tracer flux
REAL(KIND=JPRB) :: PEXTR2  (KLON,KFLDX2)     ! 2D extra variable
REAL(KIND=JPRB) :: PEXTRA  (KLON,KLEVX,KFLDX)! 3D extra variable


REAL(KIND=JPRB) ::    ZUSTRTI(KLON,KTILES),ZVSTRTI(KLON,KTILES),&
 & ZAHFSTI(KLON,KTILES),ZEVAPTI(KLON,KTILES),&
 & ZTSKTI(KLON,KTILES)  
REAL(KIND=JPRB) ::    ZUM1(KLON,KLEV) ,ZVM1(KLON,KLEV) ,&
 & ZTM1(KLON,KLEV) ,&
 & ZQM1(KLON,KLEV) , ZLM1(KLON,KLEV), ZIM1(KLON,KLEV), &
 & ZAM1(KLON,KLEV) , ZCM1(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) ::    ZDIFTQ(KLON,0:KLEV) , ZDIFTL(KLON,0:KLEV),&
 & ZDIFTI(KLON,0:KLEV)  , ZDIFTS(KLON,0:KLEV),&
 & ZSTRTU(KLON,0:KLEV)  , ZSTRTV(KLON,0:KLEV),& 
 & ZSTRSOU(KLON,0:KLEV) , ZSTRSOV(KLON,0:KLEV)  
REAL(KIND=JPRB) ::    ZTSKE1(KLON)
REAL(KIND=JPRB) ::    ZQE(KLON,KLEV,2), ZLE(KLON,KLEV,2), ZIE(KLON,KLEV,2), &
 & ZAE(KLON,KLEV,2), ZTE(KLON,KLEV,2)
REAL(KIND=JPRB) ::    ZVOM(KLON,KLEV,2) ,ZVOL(KLON,KLEV,2)
REAL(KIND=JPRB) ::    ZTENC(KLON,KLEV,MAX(KTRAC,1),2)
REAL(KIND=JPRB) ::    ZTSKE1A(KLON)

REAL(KIND=JPRB) ::    ZTSKM1M(KLON)
REAL(KIND=JPRB) ::    ZT2M(KLON)       ,ZD2M(KLON)       ,ZQ2M(KLON)       ,&
 & ZBLH(KLON)       ,ZU10M(KLON)      ,ZV10M(KLON)      ,&
 & ZGUST(KLON)  
REAL(KIND=JPRB) ::    ZAHFLEV(KLON)    ,ZAHFLSB(KLON)
REAL(KIND=JPRB) ::    ZFWSB(KLON)      ,ZEVAPSNW(KLON)
REAL(KIND=JPRB) ::    ZKH(KLON,KLEV)   ,ZKM(KLON,KLEV)
REAL(KIND=JPRB) ::    ZVDIS(KLON)      ,ZVDISG(KLON)

INTEGER(KIND=JPIM) :: INVDF, JCNT, JROF, JLEV, J1, JTR

REAL(KIND=JPRB) ::    ZTSPHY, ZINVDF
REAL(KIND=JPRB) ::    ZHOOK_HANDLE


! #include "vdfmain.intfb.h"

!     ------------------------------------------------------------------


!*         1.     INITIALIZE CONSTANTS
!                 --------------------

!*         1.1    DEFINE SHORT TIME STEP FOR INNER 
!                 VERTICAL DIFFUSION LOOP 
!                 --------------------------------

IF (LHOOK) CALL DR_HOOK('VDFOUTER',0,ZHOOK_HANDLE)


INVDF = CEILING(PTSPHY/500.0_JPRB) !substeps always smaller than 500s
!INVDF = MAX(INVDF,2)               !at lease 2 sub-steps

!amk 
!INVDF = 1
!xxx

ZINVDF = 1.0_JPRB/INVDF

ZTSPHY = PTSPHY/INVDF


! setup of a few variables that are not used - default values

DO JROF=KIDIA,KFDIA
  PHLICE (JROF) = 0.0_jprb      ! lake ice thickness    (no lakes???)
  PTLICE (JROF) = 273.0_jprb    ! lake ice temperature  (no lakes???)
  PTLWML (JROF) = 273.0_jprb    ! lake mean water T     (no lakes???)
  PCVL   (JROF) = 0.5_jprb      ! low vegetation cover
  PCVH   (JROF) = 0.5_jprb      ! high vegetation cover
  PWLMX  (JROF) = 1.0_jprb      ! maximum skin reservoir capacity (~1mm = 1kg/m2) ??? needs to be done physically by TERRA
  DO J1=1,KLEVS
    PTSAM1M (JROF,J1) = 0.0_jprb! T,soil (only used when TERRA called from EDMF)
    PWSAM1M (JROF,J1) = 0.0_jprb! Q,soil (only used when TERRA called from EDMF)
  ENDDO                         ! for TERRA from EDMF calls take dominant tile lnd_prog_now%t_so_t(jc,jk,jb,1) and w_so_t
ENDDO


!*         1.2    STORE STATE VARIABLES IN LOCAL ARRAYS
!                 -------------------------------------

DO JLEV=1,KLEV
  DO JROF=KIDIA,KFDIA
    ZUM1(JROF,JLEV)=PUM1(JROF,JLEV)
    ZVM1(JROF,JLEV)=PVM1(JROF,JLEV)
    ZTM1(JROF,JLEV)=PTM1(JROF,JLEV)
  ENDDO
ENDDO
DO JLEV=1,KLEV
  DO JROF=KIDIA,KFDIA
    ZQM1(JROF,JLEV)=PQM1(JROF,JLEV)
    ZLM1(JROF,JLEV)=PLM1(JROF,JLEV)
  ENDDO
ENDDO
DO JLEV=1,KLEV
  DO JROF=KIDIA,KFDIA
    ZIM1(JROF,JLEV)=PIM1(JROF,JLEV)
    ZAM1(JROF,JLEV)=PAM1(JROF,JLEV)
  ENDDO
ENDDO

IF (KTRAC > 0 .AND. LVDFTRAC) THEN
  DO JTR=1,KTRAC
    DO JLEV=1,KLEV
      DO JROF=KIDIA,KFDIA
        ZCM1(JROF,JLEV,JTR)=PCM1(JROF,JLEV,JTR)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*         1.3    STORE FLUXES AND STRESSES IN LOCAL ARRAYS
!                 -----------------------------------------

DO JROF=KIDIA,KFDIA
 if ( PTSKM1M(JROF) > 400.0 .or. PTSKM1M(JROF) < 100.0 ) then
  write(*,*) 'vdfouter1: ', PTSKM1M(JROF), PTM1(JROF,KLEV), PTM1(JROF,KLEV-1)
 endif
ENDDO

ZTSKM1M(KIDIA:KFDIA)=PTSKM1M(KIDIA:KFDIA)

ZTSKTI(KIDIA:KFDIA,:)=PTSKTI(KIDIA:KFDIA,:)
ZUSTRTI(KIDIA:KFDIA,:)=PUSTRTI(KIDIA:KFDIA,:)
ZVSTRTI(KIDIA:KFDIA,:)=PVSTRTI(KIDIA:KFDIA,:)
ZAHFSTI(KIDIA:KFDIA,:)=PAHFSTI(KIDIA:KFDIA,:)
ZEVAPTI(KIDIA:KFDIA,:)=PEVAPTI(KIDIA:KFDIA,:)

!*         1.4    INITIALIZE LOCAL ARRAYS FOR ACCUMULATED QUANTITIES
!                 --------------------------------------------------

ZTSKE1A(KIDIA:KFDIA)=0.0_JPRB
PVDIS(KIDIA:KFDIA)=0.0_JPRB
PVDISG(KIDIA:KFDIA)=0.0_JPRB
PAHFLEV(KIDIA:KFDIA)=0.0_JPRB
PAHFLSB(KIDIA:KFDIA)=0.0_JPRB
PFWSB(KIDIA:KFDIA)=0.0_JPRB
PEVAPSNW(KIDIA:KFDIA)=0.0_JPRB

!*         1.5    INITIALIZE FOR TILE FLUXES
!                 --------------------------

PUSTRTI(KIDIA:KFDIA,:)=0.0_JPRB
PVSTRTI(KIDIA:KFDIA,:)=0.0_JPRB
PAHFSTI(KIDIA:KFDIA,:)=0.0_JPRB
PEVAPTI(KIDIA:KFDIA,:)=0.0_JPRB

!*         1.5A   INITIALIZE FOR DDH DIAGNOSTICS
!                 ------------------------------

! unused - deleted

PBIR(KIDIA:KFDIA) = 0.0_JPRB      ! initialize in absence of test vdfmain 
LDNODECP(KIDIA:KFDIA) = .FALSE.   ! allow decoupling and 
                                  ! do not recalculate PBIR

!*         1.5B   INITIALIZE DUMMY TESSEL VARIABLES
!                 ---------------------------------

DO JROF=KIDIA,KFDIA
  KTVL  (JROF) = 0                                   ! KTVL: dummy default, not used
  KTVH  (JROF) = 0                                   ! KTVH: dummy default, not used 
  PTICE (JROF) = 0.0_JPRB                            ! only used in sfcexcdriver (deactivated)
  PTSKE1(JROF) = 0.0_JPRB                            ! only used in sfcexcdriver (deactivated)
  PSST  (JROF) = 0.0_JPRB                            ! SST: not used any more
ENDDO
KDHVTLS = 3                                          !  DDH dimensions
KDHFTLS = 8                                          !   - " -
KDHVTSS = 6                                          !   - " -
KDHFTSS = 9                                          !   - " -
KDHVTTS = 4                                          !   - " -
KDHFTTS = 11                                         !   - " -
KDHVTIS = 4                                          !   - " -
KDHFTIS = 9                                          !   - " -
CDCONF  = 'T'

!  INNER TIME LOOP FOR VERTICAL DIFFUSION

INNER_TIME_LOOP: DO JCNT=1,INVDF

  KCNT=JCNT

  J1=1
  IF(JCNT > 1) J1=2

!*         1.6    STORE TENDENCIES IN LOCAL ARRAYS
!                 --------------------------------

  DO JLEV=1,KLEV
    DO JROF=KIDIA,KFDIA
      ZQE(JROF,JLEV,J1) =PQE(JROF,JLEV)
      ZLE(JROF,JLEV,J1) =PLE(JROF,JLEV)
      ZIE(JROF,JLEV,J1) =PIE(JROF,JLEV)
    ENDDO
  ENDDO
  DO JLEV=1,KLEV
    DO JROF=KIDIA,KFDIA
      ZAE(JROF,JLEV,J1) =PAE(JROF,JLEV)
      ZTE(JROF,JLEV,J1) =PTE(JROF,JLEV)
    ENDDO
  ENDDO
  DO JLEV=1,KLEV
    DO JROF=KIDIA,KFDIA
      ZVOM(JROF,JLEV,J1)=PVOM(JROF,JLEV)
      ZVOL(JROF,JLEV,J1)=PVOL(JROF,JLEV)
    ENDDO
  ENDDO
  ZTSKE1(KIDIA:KFDIA)=PTSKE1(KIDIA:KFDIA)

  IF (KTRAC > 0 .AND. LVDFTRAC) THEN
    DO JTR=1,KTRAC
      DO JLEV=1,KLEV
        DO JROF=KIDIA,KFDIA
          ZTENC(JROF,JLEV,JTR,J1)=PTENC(JROF,JLEV,JTR)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!*         2.2    CALL VDFMAIN EACH SMALL TIME STEP
!                 ---------------------------------

DO JROF=KIDIA,KFDIA
 if ( ZTSKM1M(JROF) > 400.0 .or. ZTSKM1M(JROF) < 100.0  ) then
  write(*,*) 'vdfouter2: ', ZTSKM1M(JROF), PTM1(JROF,KLEV), PTM1(JROF,KLEV-1)
 endif
ENDDO

  CALL VDFMAIN ( CDCONF, &
   & KIDIA  , KFDIA  , KLON   , KLEV   , KLEVS  , KSTEP  , KTILES , &
   & KTRAC  , KLEVSN , KLEVI  , KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS, &
   & KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, &
   & ZTSPHY , KTVL   , KTVH   , KCNT   , PCVL   , PCVH   , PSIGFLT, &
   & ZUM1   , ZVM1   , ZTM1   , ZQM1   , ZLM1   , ZIM1   , ZAM1   , ZCM1   , &
   & PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , ZTSKM1M, PTSAM1M, PWSAM1M, &
   & PSSRFL , PSLRFL , PEMIS  , PHRLW  , PHRSW  , &
   & PTSNOW , PTICE  , &
   & PHLICE , PTLICE , PTLWML , &
   & PSST   , KSOTY  , PFRTI  , PALBTI , PWLMX  , &
   & PCHAR  , PUCURR , PVCURR , PTSKRAD, PCFLX  , &
   & PSOTEU , PSOTEV , PSOBETA, PVERVEL, &
   ! OUTPUT
   & PZ0M   , PZ0H   , &
   & ZVDIS  , ZVDISG , PDISGW3D,ZAHFLEV, ZAHFLSB, ZFWSB  , PBIR   , PVAR   , &
   & ZU10M  , ZV10M  , ZT2M   , ZD2M   , ZQ2M   , PZINV  , ZBLH   , KHPBLN , KVARTOP , &
   & PSSRFLTI,ZEVAPSNW,ZGUST  , PWUAVG , LDNODECP,KPBLTYPE, PLDIFF, &
   & PFPLVL , PFPLVN , PFHPVL , PFHPVN , &
   ! DIAGNOSTIC OUTPUT
   & PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , LLDIAG , &
   ! OUTPUT TENDENDCIES
   & ZTE(:,:,J1)     , ZQE (:,:,J1)    , ZLE (:,:,J1)    , ZIE  (:,:,J1)   , &
   & ZAE(:,:,J1)     , ZVOM(:,:,J1)    , ZVOL(:,:,J1)    , ZTENC(:,:,:,J1) , &
   & ZTSKE1 , &
   ! UPDATED FIELDS FOR TILES
   & ZUSTRTI, ZVSTRTI, ZAHFSTI, ZEVAPTI, ZTSKTI , &
   ! FLUX OUTPUTS
   & ZDIFTS , ZDIFTQ , ZDIFTL , ZDIFTI , ZSTRTU , ZSTRTV , PTOFDU , PTOFDV ,&
   & ZSTRSOU, ZSTRSOV, ZKH    , ZKM    , &
   & LDLAND , &   
! surface fluxes from TERRA for EDMF atmospheric transport
   & SHFL_S , LHFL_S , UMFL_S , VMFL_S , &
   ! TERRA data
   & frac_t , t_g_t  , qv_s   , t_ice )


!*         3.0    UPDATE STATE VARIABLES
!                 ----------------------

  IF(JCNT < INVDF) THEN
    DO JLEV=1,KLEV
      DO JROF=KIDIA,KFDIA
        ZUM1(JROF,JLEV)=ZUM1(JROF,JLEV)+ZVOM(JROF,JLEV,J1)*ZTSPHY
        ZVM1(JROF,JLEV)=ZVM1(JROF,JLEV)+ZVOL(JROF,JLEV,J1)*ZTSPHY
        ZTM1(JROF,JLEV)=ZTM1(JROF,JLEV)+ZTE(JROF,JLEV,J1)*ZTSPHY
        ZQM1(JROF,JLEV)=ZQM1(JROF,JLEV)+ZQE(JROF,JLEV,J1)*ZTSPHY
        ZLM1(JROF,JLEV)=ZLM1(JROF,JLEV)+ZLE(JROF,JLEV,J1)*ZTSPHY
        ZIM1(JROF,JLEV)=ZIM1(JROF,JLEV)+ZIE(JROF,JLEV,J1)*ZTSPHY
        ZAM1(JROF,JLEV)=ZAM1(JROF,JLEV)+ZAE(JROF,JLEV,J1)*ZTSPHY
      ENDDO
    ENDDO
  ENDIF

  IF (JCNT < INVDF .AND. KTRAC > 0 .AND. LVDFTRAC) THEN
    DO JTR=1,KTRAC
      DO JLEV=1,KLEV
        DO JROF=KIDIA,KFDIA
          ZCM1(JROF,JLEV,JTR)=ZCM1(JROF,JLEV,JTR)+ZTENC(JROF,JLEV,JTR,J1)*ZTSPHY
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!*         3.1    UPDATE SKIN TEMPERATURE
!                 -----------------------

DO JROF=KIDIA,KFDIA
 if ( ZTSKM1M(JROF) > 400.0 .or. ZTSKM1M(JROF) < 100.0 ) then
  write(*,*) 'vdfouter3: ', ZTSKM1M(JROF), PTM1(JROF,KLEV), PTM1(JROF,KLEV-1)
 endif
ENDDO

  ZTSKM1M(KIDIA:KFDIA)=ZTSKM1M(KIDIA:KFDIA)+ZTSKE1(KIDIA:KFDIA)*ZTSPHY

DO JROF=KIDIA,KFDIA
 if ( ZTSKM1M(JROF) > 400.0 .or. ZTSKM1M(JROF) < 100.0 ) then
  write(*,*) 'vdfouter4: ', ZTSKM1M(JROF), PTM1(JROF,KLEV), ZTSKE1(JROF)
 endif
ENDDO


!*         4.0    ACCUMULATE AND COMPUTE AVERAGE TENDENCIES
!                 -----------------------------------------

  IF( J1 > 1 )  THEN
    IF( JCNT < INVDF ) THEN
      DO JLEV=1,KLEV
        DO JROF=KIDIA,KFDIA
          ZTE(JROF,JLEV,1) =ZTE(JROF,JLEV,1) +ZTE(JROF,JLEV,2)
          ZQE(JROF,JLEV,1) =ZQE(JROF,JLEV,1) +ZQE(JROF,JLEV,2)
          ZLE(JROF,JLEV,1) =ZLE(JROF,JLEV,1) +ZLE(JROF,JLEV,2)
          ZIE(JROF,JLEV,1) =ZIE(JROF,JLEV,1) +ZIE(JROF,JLEV,2)
          ZAE(JROF,JLEV,1) =ZAE(JROF,JLEV,1) +ZAE(JROF,JLEV,2)
          ZVOM(JROF,JLEV,1)=ZVOM(JROF,JLEV,1)+ZVOM(JROF,JLEV,2)
          ZVOL(JROF,JLEV,1)=ZVOL(JROF,JLEV,1)+ZVOL(JROF,JLEV,2)
        ENDDO
      ENDDO
      IF (KTRAC > 0 .AND. LVDFTRAC) THEN
        DO JTR=1,KTRAC
          DO JLEV=1,KLEV
            DO JROF=KIDIA,KFDIA
              ZTENC(JROF,JLEV,JTR,1)=ZTENC(JROF,JLEV,JTR,1)+ZTENC(JROF,JLEV,JTR,2)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ELSE
      DO JLEV=1,KLEV
        DO JROF=KIDIA,KFDIA
          PTE(JROF,JLEV) =(ZTE(JROF,JLEV,1) +ZTE(JROF,JLEV,2)) *ZINVDF
          PQE(JROF,JLEV) =(ZQE(JROF,JLEV,1) +ZQE(JROF,JLEV,2)) *ZINVDF
          PLE(JROF,JLEV) =(ZLE(JROF,JLEV,1) +ZLE(JROF,JLEV,2)) *ZINVDF
        ENDDO
      ENDDO
      DO JLEV=1,KLEV
        DO JROF=KIDIA,KFDIA
          PIE(JROF,JLEV) =(ZIE(JROF,JLEV,1) +ZIE(JROF,JLEV,2)) *ZINVDF
          PAE(JROF,JLEV) =(ZAE(JROF,JLEV,1) +ZAE(JROF,JLEV,2)) *ZINVDF
        ENDDO
      ENDDO
      DO JLEV=1,KLEV
        DO JROF=KIDIA,KFDIA
          PVOM(JROF,JLEV)=(ZVOM(JROF,JLEV,1)+ZVOM(JROF,JLEV,2))*ZINVDF
          PVOL(JROF,JLEV)=(ZVOL(JROF,JLEV,1)+ZVOL(JROF,JLEV,2))*ZINVDF
        ENDDO
      ENDDO
      IF (KTRAC > 0 .AND. LVDFTRAC) THEN
        DO JTR=1,KTRAC
          DO JLEV=1,KLEV
            DO JROF=KIDIA,KFDIA
              PTENC(JROF,JLEV,JTR)=(ZTENC(JROF,JLEV,JTR,1)+ZTENC(JROF,JLEV,JTR,2))*ZINVDF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDIF

  IF( INVDF == 1 ) THEN
    DO JLEV=1,KLEV
      DO JROF=KIDIA,KFDIA
        PTE(JROF,JLEV) =ZTE(JROF,JLEV,1)
        PQE(JROF,JLEV) =ZQE(JROF,JLEV,1)
        PLE(JROF,JLEV) =ZLE(JROF,JLEV,1)
        PIE(JROF,JLEV) =ZIE(JROF,JLEV,1)
        PAE(JROF,JLEV) =ZAE(JROF,JLEV,1)
        PVOM(JROF,JLEV)=ZVOM(JROF,JLEV,1)
        PVOL(JROF,JLEV)=ZVOL(JROF,JLEV,1)
      ENDDO
    ENDDO
  ENDIF

  ZTSKE1A(KIDIA:KFDIA)=ZTSKE1A(KIDIA:KFDIA)&
   & +ZTSKE1(KIDIA:KFDIA)-PTSKE1(KIDIA:KFDIA)  

!*         4.1    ACCUMULATE AND COMPUTE AVERAGE FLUXES
!                 -------------------------------------

  IF(JCNT == 1) THEN
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PDIFTS(JROF,JLEV)=ZDIFTS(JROF,JLEV)*ZINVDF
        PDIFTQ(JROF,JLEV)=ZDIFTQ(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PDIFTL(JROF,JLEV)=ZDIFTL(JROF,JLEV)*ZINVDF
        PDIFTI(JROF,JLEV)=ZDIFTI(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PSTRTU(JROF,JLEV)=ZSTRTU(JROF,JLEV)*ZINVDF
        PSTRTV(JROF,JLEV)=ZSTRTV(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PSTRSOU(JROF,JLEV)=ZSTRSOU(JROF,JLEV)*ZINVDF
        PSTRSOV(JROF,JLEV)=ZSTRSOV(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
  ELSE
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PDIFTS(JROF,JLEV)=PDIFTS(JROF,JLEV)+ZDIFTS(JROF,JLEV)*ZINVDF
        PDIFTQ(JROF,JLEV)=PDIFTQ(JROF,JLEV)+ZDIFTQ(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PDIFTL(JROF,JLEV)=PDIFTL(JROF,JLEV)+ZDIFTL(JROF,JLEV)*ZINVDF
        PDIFTI(JROF,JLEV)=PDIFTI(JROF,JLEV)+ZDIFTI(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PSTRTU(JROF,JLEV)=PSTRTU(JROF,JLEV)+ZSTRTU(JROF,JLEV)*ZINVDF
        PSTRTV(JROF,JLEV)=PSTRTV(JROF,JLEV)+ZSTRTV(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JROF=KIDIA,KFDIA
        PSTRSOU(JROF,JLEV)=PSTRSOU(JROF,JLEV)+ZSTRSOU(JROF,JLEV)*ZINVDF
        PSTRSOV(JROF,JLEV)=PSTRSOV(JROF,JLEV)+ZSTRSOV(JROF,JLEV)*ZINVDF
      ENDDO
    ENDDO
  ENDIF

!*         4.2    ACCUMULATE FLUXES AND DERIVATIVES OF SURF QUANTITIES
!                 ----------------------------------------------------

  PVDIS(KIDIA:KFDIA)=PVDIS(KIDIA:KFDIA)+ZVDIS(KIDIA:KFDIA)
  PVDISG(KIDIA:KFDIA)=PVDISG(KIDIA:KFDIA)+ZVDISG(KIDIA:KFDIA)
  PAHFLEV(KIDIA:KFDIA)=PAHFLEV(KIDIA:KFDIA)+ZAHFLEV(KIDIA:KFDIA)
  PAHFLSB(KIDIA:KFDIA)=PAHFLSB(KIDIA:KFDIA)+ZAHFLSB(KIDIA:KFDIA)
  PFWSB(KIDIA:KFDIA)=PFWSB(KIDIA:KFDIA)+ZFWSB(KIDIA:KFDIA)
  PEVAPSNW(KIDIA:KFDIA)=PEVAPSNW(KIDIA:KFDIA)+ZEVAPSNW(KIDIA:KFDIA)

!*         4.3    ACCUMULATE TILE FLUXES
!                 ----------------------

  DO J1=1,KTILES
    DO JROF=KIDIA,KFDIA
      PUSTRTI(JROF,J1)=PUSTRTI(JROF,J1)+ZUSTRTI(JROF,J1) 
      PVSTRTI(JROF,J1)=PVSTRTI(JROF,J1)+ZVSTRTI(JROF,J1) 
      PAHFSTI(JROF,J1)=PAHFSTI(JROF,J1)+ZAHFSTI(JROF,J1) 
      PEVAPTI(JROF,J1)=PEVAPTI(JROF,J1)+ZEVAPTI(JROF,J1) 
    ENDDO
  ENDDO

!*         4.4    ACCUMULATE DDH DIAGNOSTICS
!                 --------------------------

! unused - deleted

!*         5.0    VARIABLES CLOSE TO INITIAL TIME
!                 -------------------------------

  IF (JCNT == 1) THEN
    PU10M(KIDIA:KFDIA)=ZU10M(KIDIA:KFDIA)
    PV10M(KIDIA:KFDIA)=ZV10M(KIDIA:KFDIA)
    PT2M(KIDIA:KFDIA)=ZT2M(KIDIA:KFDIA)
    PD2M(KIDIA:KFDIA)=ZD2M(KIDIA:KFDIA)
    PQ2M(KIDIA:KFDIA)=ZQ2M(KIDIA:KFDIA)
    PBLH(KIDIA:KFDIA)=ZBLH(KIDIA:KFDIA)
    PGUST(KIDIA:KFDIA)=ZGUST(KIDIA:KFDIA)
    PKH(KIDIA:KFDIA,1:KLEV)=ZKH(KIDIA:KFDIA,1:KLEV)
    PKM(KIDIA:KFDIA,1:KLEV)=ZKM(KIDIA:KFDIA,1:KLEV)
  ENDIF

ENDDO INNER_TIME_LOOP


!*         6.0    COMPUTE AVERAGE FLUXES
!                 ----------------------

!*         6.1    FOR SURF FLUXES AND DISSIPATION
!                 -------------------------------

PVDIS(KIDIA:KFDIA)   =PVDIS(KIDIA:KFDIA)   *ZINVDF
PVDISG(KIDIA:KFDIA)  =PVDISG(KIDIA:KFDIA)  *ZINVDF
PAHFLEV(KIDIA:KFDIA) =PAHFLEV(KIDIA:KFDIA) *ZINVDF
PAHFLSB(KIDIA:KFDIA) =PAHFLSB(KIDIA:KFDIA) *ZINVDF
PFWSB(KIDIA:KFDIA)   =PFWSB(KIDIA:KFDIA)   *ZINVDF
PEVAPSNW(KIDIA:KFDIA)=PEVAPSNW(KIDIA:KFDIA)*ZINVDF

!*         6.2    FOR TILE SURF FLUXES AND T,SKIN
!                 -------------------------------

PTSKTI(KIDIA:KFDIA,:) =ZTSKTI(KIDIA:KFDIA,:)
PUSTRTI(KIDIA:KFDIA,:)=PUSTRTI(KIDIA:KFDIA,:)*ZINVDF
PVSTRTI(KIDIA:KFDIA,:)=PVSTRTI(KIDIA:KFDIA,:)*ZINVDF
PAHFSTI(KIDIA:KFDIA,:)=PAHFSTI(KIDIA:KFDIA,:)*ZINVDF
PEVAPTI(KIDIA:KFDIA,:)=PEVAPTI(KIDIA:KFDIA,:)*ZINVDF

!*         6.3    COMPUTE T,SKIN AVERAGE TENDENCY
!                 -------------------------------

PTSKE1(KIDIA:KFDIA)=PTSKE1(KIDIA:KFDIA)+ZTSKE1A(KIDIA:KFDIA)*ZINVDF



IF (LHOOK) CALL DR_HOOK('VDFOUTER',1,ZHOOK_HANDLE)
END SUBROUTINE VDFOUTER


END MODULE mo_vdfouter
