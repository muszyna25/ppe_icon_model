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

SUBROUTINE VDFOUTER   ( CDCONF, &
 & KIDIA  , KFDIA  , KLON   , KLEV   , KLEVS  , KSTEP  , KTILES , &
 & KTRAC  , KLEVSN , KLEVI  , KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS, &
 & KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, &
 & PTSPHY , KTVL   , KTVH   , KCNT   , PCVL   , PCVH   , PSIGFLT, &
 & PUM1   , PVM1   , PTM1   , PQM1   , PLM1   , PIM1   , PAM1   , PCM1   , &
 & PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , PTSKM1M, PTSAM1M, PWSAM1M, &
 & PSSRFL , PSLRFL , PEMIS  , PHRLW  , PHRSW  , &
 & PTSNOW , PTICE  , &
 & PHLICE , PTLICE , PTLWML , &
 & PSST   , KSOTY  , PFRTI  , PALBTI , PWLMX  , &
 & PCHAR  , PUCURR , PVCURR , PTSKRAD, PCFLX  , &
 & PSOTEU , PSOTEV , PSOBETA, PVERVEL, &
 ! OUTPUT
 & PZ0M   , PZ0H   , &
 & PVDIS  , PVDISG , PDISGW3D,PAHFLEV, PAHFLSB, PFWSB  , PBIR   , PVAR   , &
 & PU10M  , PV10M  , PT2M   , PD2M   , PQ2M   , PZINV  , PBLH   , KHPBLN , KVARTOP , &
 & PSSRFLTI,PEVAPSNW,PGUST  , PWUAVG , LDNODECP,KPBLTYPE, PLDIFF,&
 & PFPLVL , PFPLVN , PFHPVL , PFHPVN , &
 ! DIAGNOSTIC OUTPUT
 & PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , LLDIAG,&
 ! OUTPUT TENDENCIES
 & PTE    , PQE    , PLE    , PIE    , PAE    , PVOM   , PVOL   , &
 & PTENC  , PTSKE1 , &
 ! UPDATED FIELDS FOR TILES
 & PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI , &
 ! OUTPUT FLUXES
 & PDIFTS , PDIFTQ , PDIFTL , PDIFTI , PSTRTU , PSTRTV , PTOFDU , PTOFDV, &
 & PSTRSOU, PSTRSOV, PKH    , &
!amk
 & LDLAND , &
!xxx
! DDH OUTPUTS
 & PDHTLS , PDHTSS , PDHTTS , PDHTIS &
! TERRA data
 & , ext_data                                                           & !in
 & , jb, jg                                                             & ! -
 & , t_snow_ex, t_snow_mult_ex, t_s_ex, t_g_ex, qv_s_ex                 & !inout
 & , w_snow_ex                                                          & ! -
 & , rho_snow_ex, rho_snow_mult_ex, h_snow_ex, w_i_ex, w_p_ex, w_s_ex   & ! -
 & , t_so_ex, w_so_ex, w_so_ice_ex  &  !, t_2m_ex, u_10m_ex, v_10m_ex   & ! -
 & , freshsnow_ex, snowfrac_lc_ex, snowfrac_ex                          & ! -
 & , wliq_snow_ex, wtot_snow_ex, dzh_snow_ex                            & ! -
 & , prr_con_ex, prs_con_ex, prr_gsp_ex, prs_gsp_ex                     & !in
 & , tch_ex, tcm_ex, tfv_ex                                             & !inout
 & , sobs_ex, thbs_ex, pabs_ex                                          & !in
 & , runoff_s_ex, runoff_g_ex                                           & !inout
 & , t_g, qv_s                                                          & ! -
 & , t_ice, h_ice, t_snow_si, h_snow_si, alb_si                         & ! -
 & , fr_seaice                                                          ) !in
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
!    *KDHVTLS*      Number of variables for individual tiles
!    *KDHFTLS*      Number of fluxes for individual tiles
!    *KDHVTSS*      Number of variables for snow energy budget
!    *KDHFTSS*      Number of fluxes for snow energy budget
!    *KDHVTTS*      Number of variables for soil energy budget
!    *KDHFTTS*      Number of fluxes for soil energy budget
!    *KDHVTIS*      Number of variables for sea ice energy budget
!    *KDHFTIS*      Number of fluxes for sea ice energy budget

!    *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!    *KSOTY*        SOIL TYPE                                    (1-7)

!    *KCNT*         Index of vdf sub steps.

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
!    *PTICE*        ICE TEMPERATURE (TOP SLAB)                    K
!    *PTLICE*       LAKE ICE TEMPERATURE                          K
!    *PHLICE*       LAKE ICE THICKNESS                            m
!    *PTLWML*       LAKE MEAN WATER TEMPERATURE                   K
!    *PSST*         (OPEN) SEA SURFACE TEMPERATURE                K
!    *PFRTI*        TILE FRACTIONS                                (0-1)
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
!    *PTSKE1*       SKIN TEMPERATURE TENDENCY                     K/S
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
!    *PDHTLS*       Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!    *PDHTSS*       Diagnostic array for snow T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTTS*       Diagnostic array for soil T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTIS*       Diagnostic array for ice T (see module yomcdh)
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
USE mo_lnd_nwp_config,ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water
USE mo_ext_data_types,ONLY: t_external_data

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTIS 
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(:)                   !(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(:)                   !(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:)                  !(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCNT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(:)                !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(:,:,:)               !(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(:,:)                !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(:,:)               !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(:,0:)               !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)                !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:)              !(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(:,:)              !(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(:,:)                !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(:,:)                !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:)                 !(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(:)                 !(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(:)                 !(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRTI(:,:)                !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(:,:)               !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(:)                !(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(:,:)                !(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(:,:)               !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(:,:)               !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOBETA(:,:)              !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDISG(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISGW3D(:,:)             !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLEV(:)                !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLSB(:)                !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWSB(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIR(:)                   !(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVAR(:,:)                 !(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(:)                   !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(:)                  !(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLH(:)                   !(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KHPBLN(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSSRFLTI(:,:)             !(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(:)               !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVL(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVN(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPVL(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPVN(:,0:)              !(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(:)                  !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUAVG(:)                 !(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDNODECP(:)               !(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(:)               !(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLDIFF(:,:)               !(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KVARTOP(:)                !(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PIE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAE(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(:,:)                 !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(:,:,:)              !(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(:)                 !(KLON) 
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
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDU(:)                 !(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDV(:)                 !(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOU(:,0:)             !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOV(:,0:)             !(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(:,:)                  !(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   ,optional :: PDHTLS(:,:,:)   !(KLON,KTILES,KDHVTLS+KDHFTLS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   ,optional :: PDHTSS(:,:,:)   !(KLON,KLEVSN,KDHVTSS+KDHFTSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   ,optional :: PDHTTS(:,:,:)   !(KLON,KLEVS,KDHVTTS+KDHFTTS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   ,optional :: PDHTIS(:,:,:)   !(KLON,KLEVI,KDHVTIS+KDHFTIS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(:,:)              !(KLON,KLEV) 
!amk
LOGICAL                          :: LDLAND(:)                 !(KLON)
!xxx
!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(:,:)               !(KLON,KFLDX2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(:,:,:)             !(KLON,KLEVX,KFLDX)
LOGICAL           ,INTENT(IN)    :: LLDIAG

! TERRA data

INTEGER          ,INTENT(IN)                                               :: &
  jb             ,jg                 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,0:nlev_snow,ntiles_total) :: &
  t_snow_mult_ex 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_snow,ntiles_total)   :: &
  rho_snow_mult_ex  
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  t_g_ex         ,qv_s_ex  
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total)             :: &
  t_snow_ex      ,t_s_ex         ,                                            & 
  w_snow_ex      ,rho_snow_ex    ,h_snow_ex       ,                           &
  w_i_ex         ,w_p_ex         ,w_s_ex
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,0:nlev_soil,ntiles_total) :: &
  t_so_ex             
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_soil,ntiles_total)   :: &
  w_so_ex        ,w_so_ice_ex          
!REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                         :: &
! u_10m_ex       ,v_10m_ex             
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total)             :: &
  freshsnow_ex   ,snowfrac_lc_ex ,snowfrac_ex 
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_snow,ntiles_total)   :: &
  wliq_snow_ex   ,wtot_snow_ex   ,dzh_snow_ex          
REAL(KIND=JPRB)  ,INTENT(IN)     ,DIMENSION(KLON)                          :: &
  prr_con_ex     ,prs_con_ex     ,prr_gsp_ex     ,prs_gsp_ex           
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  tch_ex         ,tcm_ex         ,tfv_ex               
REAL(KIND=JPRB)  ,INTENT(IN)     ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  sobs_ex        ,thbs_ex        ,pabs_ex              
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total)             :: &
  runoff_s_ex    ,runoff_g_ex        
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  t_g            ,qv_s
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  t_ice          ,h_ice          ,t_snow_si      ,h_snow_si       ,alb_si
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  fr_seaice
TYPE(t_external_data), INTENT(INOUT)                                       :: &
  ext_data

!*         0.2    LOCAL VARIABLES

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
REAL(KIND=JPRB) ::    ZKH(KLON,KLEV)
REAL(KIND=JPRB) ::    ZVDIS(KLON)      ,ZVDISG(KLON)

REAL(KIND=JPRB) ::    ZDHTLS(KLON,KTILES,KDHVTLS+KDHFTLS)
REAL(KIND=JPRB) ::    ZDHTSS(KLON,KLEVSN,KDHVTSS+KDHFTSS)
REAL(KIND=JPRB) ::    ZDHTTS(KLON,KLEVS,KDHVTTS+KDHFTTS)
REAL(KIND=JPRB) ::    ZDHTIS(KLON,KLEVI,KDHVTIS+KDHFTIS)

INTEGER(KIND=JPIM) :: INVDF, JCNT, JROF, JLEV, J1, J2, JTR

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

ZTSPHY=PTSPHY/INVDF

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

if ( present(PDHTLS) ) then
  PDHTLS(KIDIA:KFDIA,:,KDHVTLS+1:KDHVTLS+KDHFTLS)=0.0_JPRB
  PDHTSS(KIDIA:KFDIA,:,KDHVTSS+1:KDHVTSS+4)=0.0_JPRB
  PDHTTS(KIDIA:KFDIA,:,KDHVTTS+1:KDHVTTS+4)=0.0_JPRB
  PDHTIS(KIDIA:KFDIA,:,KDHVTIS+1:KDHVTIS+4)=0.0_JPRB
endif

PBIR(KIDIA:KFDIA) = 0.0_JPRB      ! initialize in absence of test vdfmain 
LDNODECP(KIDIA:KFDIA) = .FALSE.   ! allow decoupling and 
                                  ! do not recalculate PBIR


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
   & ZSTRSOU, ZSTRSOV, ZKH    , &
!amk
   & LDLAND , &   
!xxx
   ! DDH OUTPUTS
   & ZDHTLS , ZDHTSS , ZDHTTS , ZDHTIS &
   ! TERRA data
   & , ext_data                                                           & !in
   & , jb, jg                                                             & ! -
   & , t_snow_ex, t_snow_mult_ex, t_s_ex, t_g_ex, qv_s_ex                 & !inout
   & , w_snow_ex                                                          & ! -
   & , rho_snow_ex, rho_snow_mult_ex, h_snow_ex, w_i_ex, w_p_ex, w_s_ex   & ! -
   & , t_so_ex, w_so_ex, w_so_ice_ex  &    !, t_2m_ex, u_10m_ex, v_10m_ex & ! -
   & , freshsnow_ex, snowfrac_lc_ex, snowfrac_ex                          & ! -
   & , wliq_snow_ex, wtot_snow_ex, dzh_snow_ex                            & ! -
   & , prr_con_ex, prs_con_ex, prr_gsp_ex, prs_gsp_ex                     & !in
   & , tch_ex, tcm_ex, tfv_ex                                             & !inout
   & , sobs_ex, thbs_ex, pabs_ex                                          & !in
   & , runoff_s_ex, runoff_g_ex                                           & !inout
   & , t_g, qv_s                                                          & ! -
   & , t_ice, h_ice, t_snow_si, h_snow_si, alb_si                         & ! -
   & , fr_seaice                                                          ) !in


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

  if ( present(PDHTLS) ) then
  DO J2=KDHVTLS+1,KDHVTLS+KDHFTLS
    DO J1=1,KTILES
      DO JROF=KIDIA,KFDIA
        PDHTLS(JROF,J1,J2)=PDHTLS(JROF,J1,J2)+ZDHTLS(JROF,J1,J2)*ZINVDF
      ENDDO
    ENDDO
  ENDDO
  DO J2=KDHVTSS+1,KDHVTSS+4
    DO J1=1,KLEVSN
      DO JROF=KIDIA,KFDIA
        PDHTSS(JROF,J1,J2)=PDHTSS(JROF,J1,J2)+ZDHTSS(JROF,J1,J2)*ZINVDF
      ENDDO
    ENDDO
  ENDDO
  DO J2=KDHVTTS+1,KDHVTTS+4
    DO J1=1,KLEVS
      DO JROF=KIDIA,KFDIA
        PDHTTS(JROF,J1,J2)=PDHTTS(JROF,J1,J2)+ZDHTTS(JROF,J1,J2)*ZINVDF
      ENDDO
    ENDDO
  ENDDO
  DO J2=KDHVTIS+1,KDHVTIS+4
    DO J1=1,KLEVI
      DO JROF=KIDIA,KFDIA
        PDHTIS(JROF,J1,J2)=PDHTIS(JROF,J1,J2)+ZDHTIS(JROF,J1,J2)*ZINVDF
      ENDDO
    ENDDO
  ENDDO
  endif

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
! DDH diagnostics
    if ( present(PDHTLS) ) then
    PDHTLS(KIDIA:KFDIA,:,1)=ZDHTLS(KIDIA:KFDIA,:,1)*ZINVDF
    PDHTLS(KIDIA:KFDIA,:,2)=ZDHTLS(KIDIA:KFDIA,:,2)*ZINVDF
    PDHTLS(KIDIA:KFDIA,:,3)=ZDHTLS(KIDIA:KFDIA,:,3)*ZINVDF
    PDHTSS(KIDIA:KFDIA,:,6)=ZDHTSS(KIDIA:KFDIA,:,6)*ZINVDF
    endif
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


! DO JROF=KIDIA,KFDIA
!   IF ( (SUM(PAHFSTI(JROF,:)) == 0.0) .or. (SUM(PEVAPTI(JROF,:)) == 0.0) .or. &
!        (PDIFTS(JROF,KLEV)    == 0.0) .or. (PDIFTQ(JROF,KLEV)    == 0.0) ) THEN
!     write(*,*) 'vdfouter5: ', PDIFTS(JROF,KLEV), PDIFTQ(JROF,KLEV), PAHFSTI(JROF,:), PEVAPTI(JROF,:), &
!       & PFRTI(JROF,:)
!   ENDIF
! ENDDO


IF (LHOOK) CALL DR_HOOK('VDFOUTER',1,ZHOOK_HANDLE)
END SUBROUTINE VDFOUTER


END MODULE mo_vdfouter
