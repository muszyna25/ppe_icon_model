!>
!! Top for surface layer parameter calculation
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
! output data, diagnostics                  <<< deleted DDH outputs and
!dmk & , PDHTLS, PDHTSS, PDHTTS, PDHTIS &   <<< associated compute_ddh
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
 & , shfl_soil_ex, lhfl_soil_ex, shfl_snow_ex, lhfl_snow_ex             & !out
 & , shfl_s_ex   , lhfl_s_ex   , qhfl_s_ex                              &
 & , lhfl_bs_ex  , lhfl_pl_ex  , rstom_ex                               ) !out

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
! USE VSURF_MOD
! USE SURFSEB_CTL_MOD
!dmk USE VEVAP_MOD
!    USE SRFCOTWO_MOD
!xxx USE VSFLX_MOD

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RKAP     ,REPDU2   ,RZ0ICE   ,&                     !yoevdf  (& yos_exc)
      & RG       ,RD       ,RSIGMA   ,RTT      ,RETV     ,& !yomcst  (& yos_cst)
      & RCPD     ,RLVTT                                     !yomcst  (& yos_cst)
USE mo_edmf_param   ,ONLY : &
      & LEOCWA   ,LEOCCO   ,&                               !yoephy  (& yos_exc)
      & LEFLAKE  ,RH_ICE_MIN_FLK     ,&                     !yoephy  (& yos_flake)
      & FOEEW                                               !fcttrm.h (& fcsttre.h)
USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
USE mo_lnd_nwp_config,ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water, &
                                   & isub_water, isub_lake, isub_seaice
USE mo_ext_data_types,ONLY: t_external_data

USE mo_vupdz0       ,ONLY : vupdz0
USE mo_vexcs        ,ONLY : vexcs
USE mo_vsurf        ,ONLY : vsurf
USE mo_vevap        ,ONLY : vevap
USE mo_surfseb_ctl  ,ONLY : surfseb_ctl
USE mo_nwp_sfc_interface_edmf, ONLY : nwp_surface_edmf
USE mo_run_config   ,ONLY : msg_level
USE mo_data_turbdiff,ONLY : t0_melt, zt_ice


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
!      PAPMS    :    Pressure, lowest model level                     Pa
!      PGEOMLEV :    Geopotential, lowest atmospehric level           m2/s2
!      PCPTGZLEV:    DRY STATIC ENERGY, LOWEST MODEL LEVEL            J/kg
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
!dmk 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(:,:,:) 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(:,:,:) 
!    REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(:,:,:) 
!xxx REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(:,:,:) 

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
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  u_10m_ex       ,v_10m_ex             
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
  t_ice          ,h_ice          ,t_snow_si      ,h_snow_si          ,alb_si    
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON)                          :: &
  fr_seaice
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,ntiles_total+ntiles_water):: &
  shfl_soil_ex   ,lhfl_soil_ex   ,shfl_snow_ex   ,lhfl_snow_ex       ,        &
  shfl_s_ex      ,lhfl_s_ex      ,qhfl_s_ex      ,                            & 
  lhfl_bs_ex     ,rstom_ex            
REAL(KIND=JPRB)  ,INTENT(INOUT)  ,DIMENSION(KLON,nlev_soil,ntiles_total+ntiles_water) :: &
  lhfl_pl_ex
TYPE(t_external_data), INTENT(INOUT)                                       :: &
  ext_data

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
                 & ZCM(KLON,KTILES),ZCH(KLON,KTILES)                , &
                 & ZTSA(KLON)     , ZCSNW(KLON)    , ZSSRFL1(KLON)  , &
                 & ZCBLENDM(KLON) , ZCBLENDH(KLON) , ZSL(KLON)      , &
                 & ZQL(KLON)      , ZASL(KLON)     , ZBSL(KLON)     , &
                 & ZAQL(KLON)     , ZBQL(KLON)     , ZRHO(KLON)


INTEGER(KIND=JPIM) :: JL, JTILE, JT, IITT, isubs
LOGICAL :: LLINIT

REAL(KIND=JPRB) :: ZQSSN, ZCOR, ZRG, ZRTMST , &
                 & ZZ0MWMO, ZBLENDWMO, ZBLENDZ0, ZCOEF1, ZCONS1
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL         :: LLHISSR(KLON)

! Globcover 2009: tile index used in TESSEL (IFS) (see also mo_ext_data_state.f90)
!   Number of landcover classes provided by external parameter data
!   Needs to be changed into a variable if landcover classifications 
!   with a different number of classes become available
INTEGER, PARAMETER                 :: num_lcc = 23
REAL(KIND=JPRB), DIMENSION(num_lcc):: jtessel_gcv2009  ! Tessel index table GlobCover2009
                     ! jtessel     GlobCover2009
  DATA jtessel_gcv2009 / 4,   & !  1 irrigated croplands
                     &   4,   & !  2 rainfed croplands
                     &   4,   & !  3 mosaic cropland (50-70%) - vegetation (20-50%)
                     &   4,   & !  4 mosaic vegetation (50-70%) - cropland (20-50%)
                     &   6,   & !  5 closed broadleaved evergreen forest
                     &   6,   & !  6 closed broadleaved deciduous forest
                     &   6,   & !  7 open broadleaved deciduous forest
                     &   6,   & !  8 closed needleleaved evergreen forest
                     &   6,   & !  9 open needleleaved deciduous forest
                     &   6,   & ! 10 mixed broadleaved and needleleaved forest
                     &   4,   & ! 11 mosaic shrubland (50-70%) - grassland (20-50%)
                     &   4,   & ! 12 mosaic grassland (50-70%) - shrubland (20-50%)
                     &   4,   & ! 13 closed to open shrubland
                     &   4,   & ! 14 closed to open herbaceous vegetation
                     &   4,   & ! 15 sparse vegetation
                     &   6,   & ! 16 closed to open forest regulary flooded
                     &   6,   & ! 17 closed forest or shrubland permanently flooded
                     &   4,   & ! 18 closed to open grassland regularly flooded
                     &   8,   & ! 19 artificial surfaces
                     &   8,   & ! 20 bare areas
                     &   1,   & ! 21 water bodies
                     &   5,   & ! 22 permanent snow and ice
                     &   8    / ! 23 undefined

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

!debug
DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
    if ( (PTSKTI(JL,JT) > 400.0)  .or. (PTSKTI(JL,JT) < 100.0  ) ) then
      write(*,*) 'surfexc1: ', JT, PTSKTI(JL,JT), PSST(JL), PTSKM1M(JL)
    endif
  ENDDO
ENDDO
!xxxxx

!*         1.2  UPDATE Z0

CALL VUPDZ0(KIDIA,KFDIA,KLON,KTILES,KSTEP,CDCONF,&
   & PRVDIFTS,&
   & KTVL,KTVH,PCVL,PCVH,PUMLEV,PVMLEV,&
   & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,&
   & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
   & PHLICE,& 
   & PTSKTI,PCHAR,PUCURR,PVCURR,&
   & ZZ0MTI,ZZ0HTI,ZZ0QTI,ZBUOMTI,ZZDLTI,ZRAQTI)

!debug
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    IF ( ZZ0MTI(JL,JTILE) < 0.0_JPRB .OR. ZZ0MTI(JL,JTILE) > 10.0_JPRB ) THEN
      write(*,*) 'surfexc2: ', ZZ0MTI(JL,JTILE)
    ENDIF
  ENDDO
ENDDO
!xxxxx

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

  CALL VSURF(KIDIA,KFDIA,KLON,KLEVS,JTILE,&
   & KTVL,KTVH,&
   & PTMLEV  ,PQMLEV  ,PAPHMS,&
   & PTSKTI(:,JTILE),PWSAM1M,PTSAM1M,KSOTY,&
   & ZSRFD,ZRAQTI(:,JTILE),ZQSATI(:,JTILE),&
   & PQSTI(:,JTILE)  ,PDQSTI(:,JTILE)  ,&
   & ZWETB ,PCPTSTI(:,JTILE) ,ZWETL, ZWETH, ZWETHS )

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

LLINIT=.true.                !????
IITT=3                       !????

DO JTILE=1,KTILES

  CALL VEXCS(KIDIA,KFDIA,KLON,IITT,K_VMASS,LLINIT,PTSTEP,PRVDIFTS,&
   & PUMLEV,PVMLEV,PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,PCPTGZLEV,&
   & PCPTSTI(:,JTILE),ZQSATI(:,JTILE),&
   & ZZ0MTI(:,JTILE),ZZ0HTI(:,JTILE),&
   & ZZ0QTI(:,JTILE),ZZDLTI(:,JTILE),ZBUOMTI(:,JTILE),&
   & PUCURR,PVCURR,&
   & ZCFMTI(:,JTILE),PCFHTI(:,JTILE),&
   & PCFQTI(:,JTILE),ZKHLEV,&
   & ZCM(:,JTILE),   ZCH(:,JTILE))

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
  CALL VEVAP(KIDIA,KFDIA,KLON,PTSTEP,PRVDIFTS,JTILE,&
    & PWLMX ,PTMLEV  ,PQMLEV  ,PAPHMS,PTSKTI(:,JTILE),ZTSA,&
    & PQSTI(:,JTILE),PCFQTI(:,JTILE),ZWETB,ZWETL,ZWETH,ZWETHS,&
    & PCPTSTI(:,JTILE),PCSATTI(:,JTILE),PCAIRTI(:,JTILE),&
    & ZCSNW)
ENDDO

!          COMPUTE SNOW EVAPORATION FROM BELOW TREES i.e. TILE 7

! Note the use of qsat(Tsnow), rather than tile 7 skin. Skin T7 is a
! canopy temperature, definitely not what is desirable. Skin T5 can go
! up (and down ..) freely, not really what we want. The use of
! qsat (Tsnow) is tantamount to neglecting the skin effect there.


DO JL=KIDIA,KFDIA
  ZQSSN=FOEEW(PTSNOW(JL))/PAPHMS(JL)
  ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSSN)
  ZQSSN=ZQSSN*ZCOR
  PEVAPSNW(JL)=PCFQTI(JL,7)*ZCSNW(JL)*&
   & (PQMLEV(JL)-ZQSSN)*ZRG*ZRTMST
ENDDO


!-------------------------------------------------------------------------


!*         3.2a  CALL TERRA

! Assign TERRA exchange coefficients from TESSEL tiled coefficients (including snow!)

DO isubs=1,ntiles_total+ntiles_water
  DO jl=KIDIA,KFDIA
    IF ( isubs <= ntiles_total ) THEN          ! land
      IF ( ext_data%atm%frac_t(jl,jb,isubs) > 0.0_JPRB ) THEN   ! only used tiles 

        JTILE = jtessel_gcv2009(ext_data%atm%lc_class_t(jl,jb,isubs))

! BEST SOLUTION: give snow and no-snow ZCH to TERRA for fractional snow tiles!
!                (pass new argument ZCH_SNOW to terra and use in snow flox calculation)

       !IF ( snowfrac_ex(jl,isubs) > 0.5_jprb ) THEN
        IF ( snowfrac_ex(jl,isubs) > 0.0_jprb ) THEN  !safe: use small snow coefficients!
          SELECT CASE ( JTILE )
            CASE (4)
              JTILE = 5                        ! snow over low vegetation
            CASE (6)
              JTILE = 7                        ! snow over high vegetation
            CASE (8)
              JTILE = 5                        ! snow over bare ground
          END SELECT
        ENDIF

! interception layer (#3) missing ???
! debug: high veg -> low veg (snow or no snow) ???
        IF (JTILE == 6)  JTILE = 4
        IF (JTILE == 7)  JTILE = 5
! this fix is being done because the canopy resistence is not in TERRA and needs to be added
! (chapter 8.2.2 in IFS documentation) ... ra = 1 / (U*Ch), rc = canopy

      ELSE
        JTILE = 8                              ! unused tiles with frac=0, just for safety
      ENDIF

    ELSE

      IF (isubs == isub_water  ) JTILE = 1     ! ocean
      IF (isubs == isub_lake   ) JTILE = 1     ! lake (fake it as ocean????)
      IF (isubs == isub_seaice ) JTILE = 2     ! sea ice

    ENDIF

    tch_ex(jl,isubs) = ZCH(jl,JTILE)
    tcm_ex(jl,isubs) = ZCM(jl,JTILE)
    tfv_ex(jl,isubs) = 1.0_JPRB                ! laminar reduction factor for evaporation (Matthias) ????
  ENDDO
ENDDO

IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
  CALL nwp_surface_edmf (&
    ext_data         = ext_data        , & !>in  
    jb               = jb              , & ! block  
    jg               = jg              , & ! patch
    i_startidx       = KIDIA           , & ! start index for computations in the parallel program
    i_endidx         = KFDIA           , & ! end index for computations in the parallel program
    tcall_sfc_jg     = PTSTEP          , & ! time step
 !                                      
    u_ex             = PUMLEV          , & ! zonal wind speed                              ( m/s )
    v_ex             = PVMLEV          , & ! meridional wind speed                         ( m/s )
    t_ex             = PTMLEV          , & ! temperature                                   (  k  )
    qv_ex            = PQMLEV          , & ! specific water vapor content                  (kg/kg)
    p0_ex            = PAPMS           , & ! pressure lowest level                         ( Pa  ) 
    ps_ex            = PAPHMS          , & ! surface pressure                              ( Pa  )
 !  
    t_snow_ex        = t_snow_ex       , & ! temperature of the snow-surface               (  K  )
    t_snow_mult_ex   = t_snow_mult_ex  , & ! temperature of the snow-surface               (  K  )
    t_s_ex           = t_s_ex          , & ! temperature of the ground surface             (  K  )
    t_g_ex           = t_g_ex          , & ! surface temperature                           (  K  )
    qv_s_ex          = qv_s_ex         , & ! specific humidity at the surface              (kg/kg)
    w_snow_ex        = w_snow_ex       , & ! water content of snow                         (m H2O)
    rho_snow_ex      = rho_snow_ex     , & ! snow density                                  (kg/m**3)
    rho_snow_mult_ex = rho_snow_mult_ex, & ! snow density                                  (kg/m**3)
    h_snow_ex        = h_snow_ex       , & ! snow height                                   (  m  )
    w_i_ex           = w_i_ex          , & ! water content of interception water           (m H2O)
    w_p_ex           = w_p_ex          , & ! water content of pond interception water      (m H2O)
    w_s_ex           = w_s_ex          , & ! water content of interception snow            (m H2O)
    t_so_ex          = t_so_ex         , & ! soil temperature (main level)                 (  K  )
    w_so_ex          = w_so_ex         , & ! total water conent (ice + liquid water)       (m H20)
    w_so_ice_ex      = w_so_ice_ex     , & ! ice content                                   (m H20)
!   t_2m_ex          = t_2m_ex         , & ! temperature in 2m                             (  K  )
    u_10m_ex         = u_10m_ex        , & ! zonal wind in 10m                             ( m/s )
    v_10m_ex         = v_10m_ex        , & ! meridional wind in 10m                        ( m/s )
 !  
    freshsnow_ex     = freshsnow_ex    , & ! indicator for age of snow in top of snow layer(  -  )
    snowfrac_lc_ex   = snowfrac_lc_ex  , & ! snow-cover fraction                           (  -  )
    snowfrac_ex      = snowfrac_ex     , & ! snow-cover fraction                           (  -  )
    wliq_snow_ex     = wliq_snow_ex    , & ! liquid water content in the snow              (m H2O)
    wtot_snow_ex     = wtot_snow_ex    , & ! total (liquid + solid) water content of snow  (m H2O)
    dzh_snow_ex      = dzh_snow_ex     , & ! layer thickness between half levels in snow   (  m  )
 !   
    prr_con_ex       = prr_con_ex      , & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con_ex       = prs_con_ex      , & ! precipitation rate of snow, convective        (kg/m2*s)
    prr_gsp_ex       = prr_gsp_ex      , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp_ex       = prs_gsp_ex      , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
 !                                   
    tch_ex           = tch_ex          , & ! turbulent transfer coefficient for heat       ( -- )
    tcm_ex           = tcm_ex          , & ! turbulent transfer coefficient for momentum   ( -- )
    tfv_ex           = tfv_ex          , & ! laminar reduction factor for evaporation      ( -- )
 !                                   
    sobs_ex          = sobs_ex         , & ! solar radiation at the ground                 ( W/m2)
    thbs_ex          = thbs_ex         , & ! thermal radiation at the ground               ( W/m2)
    pabs_ex          = pabs_ex         , & !!!! photosynthetic active radiation            ( W/m2)
 !                                   
    runoff_s_ex      = runoff_s_ex     , & ! surface water runoff; sum over forecast       (kg/m2)
    runoff_g_ex      = runoff_g_ex     , & ! soil water runoff; sum over forecast          (kg/m2)
 !                                   
    t_g              = t_g             , & ! surface temperature (grid mean)               (  K  )
    qv_s             = qv_s            , & ! surface specific humidity (grid mean)         (kg/kg)
 !
    t_ice            = t_ice           , & ! sea ice temperature                           (  K  )
    h_ice            = h_ice           , & ! sea ice height                                (  m  )
    t_snow_si        = t_snow_si       , & ! sea ice snow temperature                      (  K  )
    h_snow_si        = h_snow_si       , & ! sea ice snow height                           (  m  )
    alb_si           = alb_si          , & ! sea-ice albedo                                (  -  )
    fr_seaice        = fr_seaice       , & ! sea ice fraction                              (  1  )
!
    shfl_soil_ex     = shfl_soil_ex    , & ! sensible heat flux soil/air interface         (W/m2)
    lhfl_soil_ex     = lhfl_soil_ex    , & ! latent   heat flux soil/air interface         (W/m2)
    shfl_snow_ex     = shfl_snow_ex    , & ! sensible heat flux snow/air interface         (W/m2)
    lhfl_snow_ex     = lhfl_snow_ex    , & ! latent   heat flux snow/air interface         (W/m2)
    shfl_s_ex        = shfl_s_ex       , & ! sensible heat flux                            (W/m2)
    lhfl_s_ex        = lhfl_s_ex       , & ! latent heat flux                              (W/m2)
    qhfl_s_ex        = qhfl_s_ex       , & ! moisture flux                                 (W/m2)
    lhfl_bs_ex       = lhfl_bs_ex      , & 
    lhfl_pl_ex       = lhfl_pl_ex      , &
    rstom_ex         = rstom_ex        ) 
ENDIF


IF (msg_level >= 15) THEN
  DO JTILE=1,ntiles_total
    DO JL=KIDIA,KFDIA
      IF ( ABS( shfl_soil_ex(jl,jtile) * (1-snowfrac_ex(jl,jtile)))  >  800.0_JPRB  .OR. & 
           ABS( shfl_snow_ex(jl,jtile) *    snowfrac_ex(jl,jtile) )  >  800.0_JPRB  .OR. & 
           ABS( lhfl_soil_ex(jl,jtile) * (1-snowfrac_ex(jl,jtile)))  > 2000.0_JPRB  .OR. & 
           ABS( lhfl_snow_ex(jl,jtile) *    snowfrac_ex(jl,jtile) )  > 2000.0_JPRB  ) THEN
         write(*,*) 'surfexc3: SHF-soil,-snow,LHF-soil,-snow', &
           jl, jtile, snowfrac_ex(jl,jtile), ext_data%atm%frac_t(jl,jb,jtile), ext_data%atm%lc_class_t(jl,jb,jtile), &
           shfl_soil_ex(jl,jtile), shfl_snow_ex(jl,jtile), &
           lhfl_soil_ex(jl,jtile), lhfl_snow_ex(jl,jtile), &
           tch_ex(jl,jtile), t_g_ex(jl,jtile), PTMLEV(jl), qv_s_ex(jl,jtile), PQMLEV(jl) 
      ENDIF
    ENDDO
  ENDDO
ENDIF

!overwrite fluxes over land from TERRA back to EDMF code
!...this needs to be done by tile properly???????
!...see also mo_vdfmain.f90: ZEXTSHF/ZEXTLHF (done there)
!?? DO JTILE=3,KTILES
!??   DO JL=KIDIA,KFDIA
!??     PAHFSTI(JL,JTILE) = 0.0_JPRB
!??     PEVAPTI(JL,JTILE) = 0.0_JPRB
!??     DO JT=1,ntiles_total+ntiles_water
!??       PAHFSTI(JL,JTILE) = PAHFSTI(JL,JTILE) + subsfrac_ex(JL,JT) * &
!??           ( SHFL_SOIL_EX   (JL,JT) * (1.0_JPRB - SNOWFRAC_EX(JL,JT)) +  &
!??             SHFL_SNOW_EX(JL,JT) *                SNOWFRAC_EX(JL,JT)  ) 
!??       PEVAPTI(JL,JTILE) = PEVAPTI(JL,JTILE) + subsfrac_ex(JL,JT) * &
!??           ( LHFL_SOIL_EX   (JL,JT) * (1.0_JPRB - SNOWFRAC_EX(JL,JT)) +  &
!??             LHFL_SNOW_EX(JL,JT) *                SNOWFRAC_EX(JL,JT)  )/RLVTT
!??     ENDDO
!??   ENDDO
!?? ENDDO

!-------------------------------------------------------------------------


!*         3.3  COMPUTE SURFACE FLUXES FOR TILES

ZCONS1 =RG*PTSTEP*PRVDIFTS
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
!   use previous times step fluxes for heat and moisture
    ZKHFLTI(JL,JTILE)=PAHFSTI(JL,JTILE)/(ZRHO(JL)*RCPD)
    ZKQFLTI(JL,JTILE)=PEVAPTI(JL,JTILE)/ZRHO(JL)

    ZKMFLTI(JL,JTILE)=ZCFMTI(JL,JTILE)*SQRT((PUMLEV(JL)-PUCURR(JL))**2&
   & +(PVMLEV(JL)-PVCURR(JL))**2)/(ZCONS1*ZRHO(JL))
  ENDDO
ENDDO


!*         3.3a   PREPARE ARRAY'S FOR CALL TO SURFACE ENERGY
!                 BALANCE ROUTINE

IF (KSTEP == 0) THEN

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

  CALL SURFSEB_CTL(KIDIA,KFDIA,KLON,KTILES,&
   & PCPTSTI,PTSKTI,PQSTI,&
   & PDQSTI,ZRHOCHU,ZRHOCQU,&
   & PCAIRTI,PCSATTI,&
   & PSSRFLTI,PFRTI,ZTSRF,&
   & PHLICE,&
   & PSLRFL,PTSKM1M,PEMIS,&
   & ZASL,ZBSL,ZAQL,ZBQL,&
   !out
   & ZJS,ZJQ,ZSSK,ZTSK,&
   & ZSSH,ZSLH,ZSTR,ZG0,&
   & ZSL,ZQL)  

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
