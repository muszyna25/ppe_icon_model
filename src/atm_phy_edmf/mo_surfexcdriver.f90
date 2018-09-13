!>
!! External interface to IFS surface scheme
!!
!! @author Martin Koehler, DWD
!!
!! @par Revision History
!! Imported from IFS by Martin Koehler  (starting 2012-4-30)
!!   (IFS cycle CY36R1)
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

MODULE mo_surfexcdriver
 
  PUBLIC :: surfexcdriver

CONTAINS

SUBROUTINE SURFEXCDRIVER    ( &
! standard input
 &   CDCONF &
 & , KIDIA, KFDIA, KLON, KLEVS, KTILES, KSTEP &
 & , KLEVSN, KLEVI, KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS &
 & , KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, K_VMASS &
 & , PTSTEP, PRVDIFTS &
! input data, non-tiled
 & , KTVL, KTVH, PCVL, PCVH &
 & , PUMLEV, PVMLEV, PTMLEV, PQMLEV, PAPHMS, PAPMS, PGEOMLEV, PCPTGZLEV &
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
! output data, diagnostics                 <<< deleted DDH outputs and
!dmk & , PDHTLS, PDHTSS, PDHTTS, PDHTIS &  <<< associated compute_ddh 
 & ) 

! USE PARKIND1  ,ONLY : JPIM, JPRB
!ifndef INTERFACE
! USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
! USE ABORT_SURF_MOD
!endif INTERFACE

!ICON definitions:
USE mo_kind             ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters     ,ONLY : lhook    ,dr_hook           !yomcst  (& yos_exc)

USE mo_lnd_nwp_config   ,ONLY : nlev_soil, nlev_snow, ntiles_total, ntiles_water
USE mo_ext_data_types   ,ONLY : t_external_data
USE mo_edmf_param       ,ONLY : abort_surf
USE mo_surfexcdriver_ctl,ONLY : surfexcdriver_ctl

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
!    E. Dutra/G. Balsamo May 2008 Add lake tile

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
!      KSOTY    :    SOIL TYPE                                        (1-7)

!    Reals (In):
!      PTSTEP   :    Timestep
!      PRVDIFTS :    Semi-implicit factor for vertical diffusion discretization
!      PCVL     :    Low vegetation fraction
!      PCVH     :    High vegetation fraction

!    Reals with tile index (In): 
!      PFRTI    :    TILE FRACTIONS                                   (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!      PALBTI   :    Tile albedo                                      (0-1)

!    Reals independent of tiles (In):
!      PUMLEV   :    X-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PVMLEV   :    Y-VELOCITY COMPONENT, lowest atmospheric level   m/s
!      PTMLEV   :    TEMPERATURE,   lowest atmospheric level          K
!      PQMLEV   :    SPECIFIC HUMIDITY                                kg/kg
!      PAPHMS   :    Surface pressure                                 Pa
!      PAPMS    :    Pressure, lowest model level                     Pa
!      PGEOMLEV :    Geopotential, lowest atmospehric level           m2/s2
!      PCPTGZLEV:    Geopotential, lowest atmospehric level           J/kg
!      PSST     :    (OPEN) SEA SURFACE TEMPERATURE                   K
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
!         *VSFLX*

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation

!------------------------------------------------------------------------

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPMS(:)
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRTI(:,:) 
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
!dmk REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(:,:,:) 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(:,:,:) 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(:,:,:) 
!xxx REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(:,:,:) 

! TERRA data

INTEGER        ::  jb             ,jg                 
REAL(KIND=JPRB)  ,DIMENSION(KLON,0:nlev_snow,ntiles_total) :: &
  t_snow_mult_ex 
REAL(KIND=JPRB)  ,DIMENSION(KLON,nlev_snow,ntiles_total)   :: &
  rho_snow_mult_ex  
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  t_g_ex         ,qv_s_ex  
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total)             :: &
  t_snow_ex      ,t_s_ex         ,                                            & 
  w_snow_ex      ,rho_snow_ex    ,h_snow_ex       ,                           &
  w_i_ex         ,w_p_ex         ,w_s_ex
REAL(KIND=JPRB)  ,DIMENSION(KLON,0:nlev_soil,ntiles_total) :: &
  t_so_ex             
REAL(KIND=JPRB)  ,DIMENSION(KLON,nlev_soil,ntiles_total)   :: &
  w_so_ex        ,w_so_ice_ex          
REAL(KIND=JPRB)  ,DIMENSION(KLON)                          :: &
  u_10m_ex       ,v_10m_ex             
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total)             :: &
  freshsnow_ex   ,snowfrac_lc_ex ,snowfrac_ex 
REAL(KIND=JPRB)  ,DIMENSION(KLON,nlev_snow,ntiles_total)   :: &
  wliq_snow_ex   ,wtot_snow_ex   ,dzh_snow_ex          
REAL(KIND=JPRB)  ,DIMENSION(KLON)                          :: &
  prr_con_ex     ,prs_con_ex     ,prr_gsp_ex     ,prs_gsp_ex           
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  tch_ex         ,tcm_ex         ,tfv_ex               
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  sobs_ex        ,thbs_ex        ,pabs_ex              
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total)             :: &
  runoff_s_ex    ,runoff_g_ex        
REAL(KIND=JPRB)  ,DIMENSION(KLON)                          :: &
  t_g            ,qv_s
REAL(KIND=JPRB)  ,DIMENSION(KLON)                          :: &
  t_ice          ,h_ice          ,t_snow_si      ,h_snow_si       , alb_si
REAL(KIND=JPRB)  ,DIMENSION(KLON)                          :: &
  fr_seaice
REAL(KIND=JPRB)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  shfl_soil_t    ,lhfl_soil_t    ,shfl_snow_t    ,lhfl_snow_t     ,           &
  shfl_s_t       ,lhfl_s_t       ,qhfl_s_t       ,                            & 
  lhfl_bs_t      ,rstom_t             
REAL(KIND=JPRB)  ,DIMENSION(KLON,nlev_soil,ntiles_total+ntiles_water) :: &
  lhfl_pl_t 
TYPE(t_external_data)                                       :: &
  ext_data

!xxx
  
!ifndef INTERFACE

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFEXCDRIVER',0,ZHOOK_HANDLE)
IF(UBOUND(KTVL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: KTVL TOO SHORT!')
ENDIF

IF(UBOUND(KTVH,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: KTVH TOO SHORT!')
ENDIF

IF(UBOUND(PCVL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PCVL TOO SHORT!')
ENDIF

IF(UBOUND(PCVH,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PCVH TOO SHORT!')
ENDIF

IF(UBOUND(PUMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PUMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PVMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PVMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PTMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PTMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PQMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PQMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PAPHMS,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PAPHMS TOO SHORT!')
ENDIF

IF(UBOUND(PAPMS,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PAPMS TOO SHORT!')
ENDIF

IF(UBOUND(PGEOMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PGEOMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PCPTGZLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PCPTGZLEV TOO SHORT!')
ENDIF

IF(UBOUND(PSST,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PSST TOO SHORT!')
ENDIF

IF(UBOUND(PTSKM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PTSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PCHAR,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PCHAR TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PSSRFL TOO SHORT!')
ENDIF

IF(UBOUND(PSLRFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PSLRFL TOO SHORT!')
ENDIF

IF(UBOUND(PEMIS,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PEMIS TOO SHORT!')
ENDIF

IF(UBOUND(PTICE,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PTICE TOO SHORT!')
ENDIF

IF(UBOUND(PTSNOW,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PTSNOW TOO SHORT!')
ENDIF

IF(UBOUND(PWLMX,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PWLMX TOO SHORT!')
ENDIF

IF(UBOUND(PUCURR,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PUCURR TOO SHORT!')
ENDIF

IF(UBOUND(PVCURR,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PVCURR TOO SHORT!')
ENDIF

IF(UBOUND(PTSAM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PTSAM1M TOO SHORT!')
ENDIF

IF(UBOUND(PTSAM1M,2) < KLEVS) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PTSAM1M TOO SHORT!')
ENDIF

IF(UBOUND(PWSAM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PWSAM1M TOO SHORT!')
ENDIF

IF(UBOUND(PWSAM1M,2) < KLEVS) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PWSAM1M TOO SHORT!')
ENDIF

IF(UBOUND(KSOTY,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF KSOTY TOO SHORT!')
ENDIF

IF(UBOUND(PFRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PFRTI TOO SHORT!')
ENDIF

IF(UBOUND(PFRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PFRTI TOO SHORT!')
ENDIF

IF(UBOUND(PALBTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PALBTI TOO SHORT!')
ENDIF

IF(UBOUND(PALBTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PALBTI TOO SHORT!')
ENDIF

IF(UBOUND(PUSTRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PUSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PUSTRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PUSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PVSTRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PVSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PVSTRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PVSTRTI TOO SHORT!')
ENDIF

IF(UBOUND(PAHFSTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PAHFSTI TOO SHORT!')
ENDIF

IF(UBOUND(PAHFSTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PAHFSTI TOO SHORT!')
ENDIF

IF(UBOUND(PEVAPTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PEVAPTI TOO SHORT!')
ENDIF

IF(UBOUND(PEVAPTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PEVAPTI TOO SHORT!')
ENDIF

IF(UBOUND(PTSKTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PTSKTI TOO SHORT!')
ENDIF

IF(UBOUND(PTSKTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PTSKTI TOO SHORT!')
ENDIF

IF(UBOUND(PZ0M,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PZ0M TOO SHORT!')
ENDIF

IF(UBOUND(PZ0H,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PZ0H TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFLTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PSSRFLTI TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFLTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PSSRFLTI TOO SHORT!')
ENDIF

IF(UBOUND(PQSTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PQSTI TOO SHORT!')
ENDIF

IF(UBOUND(PQSTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PQSTI TOO SHORT!')
ENDIF

IF(UBOUND(PDQSTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PDQSTI TOO SHORT!')
ENDIF

IF(UBOUND(PDQSTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PDQSTI TOO SHORT!')
ENDIF

IF(UBOUND(PCPTSTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PCPTSTI TOO SHORT!')
ENDIF

IF(UBOUND(PCPTSTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PCPTSTI TOO SHORT!')
ENDIF

IF(UBOUND(PCFHTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PCFHTI TOO SHORT!')
ENDIF

IF(UBOUND(PCFHTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PCFHTI TOO SHORT!')
ENDIF

IF(UBOUND(PCFQTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PCFQTI TOO SHORT!')
ENDIF

IF(UBOUND(PCFQTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PCFQTI TOO SHORT!')
ENDIF

IF(UBOUND(PCSATTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PCSATTI TOO SHORT!')
ENDIF

IF(UBOUND(PCSATTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PCSATTI TOO SHORT!')
ENDIF

IF(UBOUND(PCAIRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PCAIRTI TOO SHORT!')
ENDIF

IF(UBOUND(PCAIRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PCAIRTI TOO SHORT!')
ENDIF

IF(UBOUND(PKHLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PKHLEV TOO SHORT!')
ENDIF

IF(UBOUND(PCFMLEV,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PCFMLEV TOO SHORT!')
ENDIF

IF(UBOUND(PKMFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PKMFL TOO SHORT!')
ENDIF

IF(UBOUND(PKHFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PKHFL TOO SHORT!')
ENDIF

IF(UBOUND(PKQFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PKQFL TOO SHORT!')
ENDIF

IF(UBOUND(PEVAPSNW,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PEVAPSNW TOO SHORT!')
ENDIF

IF(UBOUND(PZ0MW,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PZ0MW TOO SHORT!')
ENDIF

IF(UBOUND(PZ0HW,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PZ0HW TOO SHORT!')
ENDIF

IF(UBOUND(PZ0QW,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PZ0QW TOO SHORT!')
ENDIF

IF(UBOUND(PBLENDPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PBLENDPP TOO SHORT!')
ENDIF

IF(UBOUND(PCPTSPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PCPTSPP TOO SHORT!')
ENDIF

IF(UBOUND(PQSAPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PQSAPP TOO SHORT!')
ENDIF

IF(UBOUND(PBUOMPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER: PBUOMPP TOO SHORT!')
ENDIF

IF(UBOUND(PZDLPP,1) < KLON) THEN
  CALL ABORT_SURF('SURFEXCDRIVER:  PZDLPPTOO SHORT!')
ENDIF

!dmk
! IF(UBOUND(PDHTLS,1) < KLON) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PDHTLS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTLS,2) < KTILES) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PDHTLS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTLS,3) < KDHVTLS+KDHFTLS) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: THIRD DIMENSION OF PDHTLS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTSS,1) < KLON) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PDHTSS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTSS,2) < KLEVSN) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PDHTSS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTSS,3) < KDHVTSS+KDHFTSS) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: THIRD DIMENSION OF PDHTSS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTTS,1) < KLON) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PDHTTS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTTS,2) < KLEVS) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PDHTTS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTTS,3) < KDHVTTS+KDHFTTS) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: THIRD DIMENSION OF PDHTTS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTIS,1) < KLON) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PDHTIS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTIS,2) < KLEVI) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: SECOND DIMENSION OF PDHTIS TOO SHORT!')
! ENDIF
! 
! IF(UBOUND(PDHTIS,3) < KDHVTIS+KDHFTIS) THEN
!   CALL ABORT_SURF('SURFEXCDRIVER:: THIRD DIMENSION OF PDHTIS TOO SHORT!')
! ENDIF
!xxx

IF(UBOUND(PHLICE,1) < KLON) THEN 
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PHLICE TOO SHORT!')
ENDIF
IF(UBOUND(PTLICE,1) < KLON) THEN 
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PTLICE TOO SHORT!')
ENDIF
IF(UBOUND(PTLWML,1) < KLON) THEN 
  CALL ABORT_SURF('SURFEXCDRIVER:: FIRST DIMENSION OF PTLWML TOO SHORT!')
ENDIF


CALL SURFEXCDRIVER_CTL(CDCONF &
 & , KIDIA, KFDIA, KLON, KLEVS, KTILES, KSTEP &
 & , KLEVSN, KLEVI, KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS &
 & , KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, K_VMASS &
 & , PTSTEP, PRVDIFTS &
 & , KTVL, KTVH, PCVL, PCVH &
 & , PUMLEV, PVMLEV, PTMLEV, PQMLEV, PAPHMS, PAPMS, PGEOMLEV, PCPTGZLEV &
 & , PSST, PTSKM1M, PCHAR, PSSRFL, PSLRFL, PEMIS, PTICE, PTSNOW &
 & , PHLICE,PTLICE,PTLWML &   
 & , PWLMX, PUCURR, PVCURR &
 & , PTSAM1M, PWSAM1M, KSOTY &
 & , PFRTI, PALBTI &
 & , PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI &
 & , PZ0M, PZ0H &
 & , PSSRFLTI, PQSTI, PDQSTI, PCPTSTI, PCFHTI, PCFQTI, PCSATTI, PCAIRTI &
 & , PKHLEV, PCFMLEV, PKMFL, PKHFL, PKQFL, PEVAPSNW &
 & , PZ0MW, PZ0HW, PZ0QW, PBLENDPP, PCPTSPP, PQSAPP, PBUOMPP, PZDLPP &
!dmk & , PDHTLS, PDHTSS, PDHTTS, PDHTIS )
! TERRA data
 & , ext_data                                                           & !in
 & , jb, jg                                                             & ! -
 & , t_snow_ex, t_snow_mult_ex, t_s_ex, t_g_ex, qv_s_ex                 & !inout
 & , w_snow_ex                                                          & ! -
 & , rho_snow_ex, rho_snow_mult_ex, h_snow_ex, w_i_ex, w_p_ex, w_s_ex   & ! -
 & , t_so_ex, w_so_ex, w_so_ice_ex, u_10m_ex, v_10m_ex                  & ! -
 & , freshsnow_ex, snowfrac_lc_ex, snowfrac_ex                          & ! -
 & , wliq_snow_ex, wtot_snow_ex, dzh_snow_ex                            & ! -
 & , prr_con_ex, prs_con_ex, prr_gsp_ex, prs_gsp_ex                     & !in
 & , tch_ex, tcm_ex, tfv_ex                                             & !inout
 & , sobs_ex, thbs_ex, pabs_ex                                          & !in
 & , runoff_s_ex, runoff_g_ex                                           & !inout
 & , t_g, qv_s                                                          & ! -
 & , t_ice, h_ice, t_snow_si, h_snow_si, alb_si                         & ! -
 & , fr_seaice                                                          & !in
 & , shfl_soil_t, lhfl_soil_t, shfl_snow_t, lhfl_snow_t                 & !out
 & , shfl_s_t   , lhfl_s_t   , qhfl_s_t                                 &
 & , lhfl_bs_t  , lhfl_pl_t  , rstom_t                                  )     

IF (LHOOK) CALL DR_HOOK('SRFEXCDRIVER',1,ZHOOK_HANDLE)

!endif INTERFACE

!------------------------------------------------------------------------

END SUBROUTINE SURFEXCDRIVER


END MODULE MO_SURFEXCDRIVER
