!>
!! This module is the interface between the EDMF turbulence to the 
!! surface parameterisations:
!! inwp_sfc  == 1 == surface scheme TERRA
!! inwp_turb == 3 == turbulence EDMF
!!
!! @author Martin Koehler, DWD, Offenbach (2012-05-04)
!!
!! @par Revision History
!! Initial Martin Koehler, DWD, Offenbach (2012-05-04)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_sfc_interface_edmf

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message !, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: zml_soil
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, nsfc_subs, t_tiles,  &
    &                               lseaice, llake, lmulti_snow
  USE mo_satad,               ONLY: sat_pres_water, spec_humi  
  USE mo_soil_ml,             ONLY: terra_multlay
  USE mo_phyparam_soil              ! soil and vegetation parameters for TILES
  
  IMPLICIT NONE 

  PUBLIC  ::  nwp_surface_terra_edmf

  PRIVATE


#ifdef __SX__
! parameters for loop unrolling
INTEGER, PARAMETER :: nlsoil= 7
INTEGER, PARAMETER :: nlsnow= 2
#endif

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!

  SUBROUTINE nwp_surface_terra_edmf (&   
                  p_patch          , & !>in   !!! DELETE
                  ext_data         , & !>in   !!! DELETE
                  jb               , & !      !!! DELETE
                  nproma           , & ! array dimensions
                  i_startidx       , & ! start index for computations in the parallel program
                  i_endidx         , & ! end index for computations in the parallel program
                  nsubs0           , & ! nsubs0=1 for single tile, nsubs0=2 for multi-tile
                  nsubs1           , & ! nsubs1=1 for single tile, nsubs1=#tiles+1 for multi-tile
                  dt               , & ! time step                                     ( s  )
                            
                  u_ex             , & ! zonal wind speed                              ( m/s )
                  v_ex             , & ! meridional wind speed                         ( m/s )
                  t_ex             , & ! temperature                                   (  k  )
                  qv_ex            , & ! specific water vapor content                  (kg/kg)
                  p0_ex            , & !!!! base state pressure                           (Pa) 
                  ps_ex            , & ! surface pressure                              ( pa  )
                                
                  t_snow_ex        , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_ex   , & ! temperature of the snow-surface               (  K  )
                  t_s_ex           , & ! temperature of the ground surface             (  K  )
                  t_g_ex           , & ! weighted surface temperature                  (  K  )
                  qv_s_ex          , & ! specific humidity at the surface              (kg/kg)
                  w_snow_ex        , & ! water content of snow                         (m H2O)
                  rho_snow_ex      , & ! snow density                                  (kg/m**3)
                  rho_snow_mult_ex , & ! snow density                                  (kg/m**3)
                  h_snow_ex        , & ! snow height                                   (  m  )
                  w_i_ex           , & ! water content of interception water           (m H2O)
                  t_so_ex          , & ! soil temperature (main level)                 (  K  )
                  w_so_ex          , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_ex      , & ! ice content                                   (m H20)
                  t_2m_ex          , & ! temperature in 2m                             (  K  )
                  u_10m_ex         , & ! zonal wind in 10m                             ( m/s )
                  v_10m_ex         , & ! meridional wind in 10m                        ( m/s )

                  freshsnow_ex     , & ! indicator for age of snow in top of snow layer(  -  )
                  wliq_snow_ex     , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_ex     , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_ex      , & ! layer thickness between half levels in snow   (  m  )
                  subsfrac_ex      , &
                                
                  prr_con_ex       , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con_ex       , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp_ex       , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp_ex       , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                              
                  tch_ex           , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm_ex           , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv_ex           , & ! laminar reduction factor for evaporation      ( -- )
                         
                  sobs_ex          , & ! solar radiation at the ground                 ( W/m2)
                  thbs_ex          , & ! thermal radiation at the ground               ( W/m2)
                  pabs_ex          , & !!!! photosynthetic active radiation            ( W/m2)
    
                  runoff_s_ex      , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g_ex      , & ! soil water runoff; sum over forecast          (kg/m2)

                  t_g              , & ! surface temperature (grid mean)               ( K )
                  qv_s               & ! surface specific humidity (grid mean)         (kg/kg)
                                     )


  IMPLICIT NONE

  INTEGER,  INTENT(IN)  ::  &
                  nproma,            & ! array dimensions
                  i_startidx,        & ! start index for computations in the parallel program
                  i_endidx,          & ! end index for computations in the parallel program
                  nsubs0,nsubs1    , &
                  jb
  REAL(wp), INTENT(IN)  ::  &
                  dt
  REAL(wp), DIMENSION(nproma), INTENT(IN) :: & 
                  u_ex             , & ! zonal wind speed                              ( m/s )
                  v_ex             , & ! meridional wind speed                         ( m/s )
                  t_ex             , & ! temperature                                   (  k  )
                  qv_ex            , & ! specific water vapor content                  (kg/kg)
                  p0_ex            , & !!!! base state pressure                        ( Pa ) 
                  ps_ex                ! surface pressure                              ( pa  )
  REAL(wp), DIMENSION(nproma,0:nlev_snow,nsfc_subs), INTENT(INOUT) :: &
                  t_snow_mult_ex   , & ! temperature of the snow-surface               (  K  )
                  rho_snow_mult_ex     ! snow density                                  (kg/m**3)
  REAL(wp), DIMENSION(nproma,nsfc_subs), INTENT(INOUT) :: &
                  t_snow_ex        , & ! temperature of the snow-surface (K)
                  t_s_ex           , & ! temperature of the ground surface             (  K  )
                  t_g_ex           , & ! weighted surface temperature                  (  K  )
                  qv_s_ex          , & ! specific humidity at the surface              (kg/kg)
                  w_snow_ex        , & ! water content of snow                         (m H2O)
                  rho_snow_ex      , & ! snow density                                  (kg/m**3)
                  h_snow_ex        , & ! snow height  
                  w_i_ex               ! water content of interception water           (m H2O)
  REAL(wp), DIMENSION(nproma,0:nlev_soil+1,nsfc_subs), INTENT(INOUT) :: &
                  t_so_ex              ! soil temperature (main level)                 (  K  )
  REAL(wp), DIMENSION(nproma,nlev_soil+1,nsfc_subs), INTENT(INOUT) :: &
                  w_so_ex          , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_ex          ! ice content                                   (m H20)
  REAL(wp), DIMENSION(nproma), INTENT(INOUT) :: &
                  t_2m_ex          , & ! temperature in 2m                             (  K  )
                  u_10m_ex         , & ! zonal wind in 10m                             ( m/s )
                  v_10m_ex             ! meridional wind in 10m                        ( m/s )
  REAL(wp), DIMENSION(nproma,nsfc_subs), INTENT(INOUT) :: &
                  freshsnow_ex         ! indicator for age of snow in top of snow layer(  -  )
  REAL(wp), DIMENSION(nproma,nlev_snow,nsfc_subs), INTENT(INOUT) :: &
                  wliq_snow_ex     , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_ex     , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_ex          ! layer thickness between half levels in snow   (  m  )
  REAL(wp), DIMENSION(nproma,nsfc_subs), INTENT(INOUT) :: &
                  subsfrac_ex
  REAL(wp), DIMENSION(nproma), INTENT(IN) ::    &
                  prr_con_ex       , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con_ex       , & ! precipitation rate of snow, convective        (kg/m2*s)
                  prr_gsp_ex       , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp_ex           ! precipitation rate of snow, grid-scale        (kg/m2*s)

  REAL(wp), DIMENSION(nproma,nsfc_subs), INTENT(INOUT) :: &
                  tch_ex           , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm_ex           , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv_ex               ! laminar reduction factor for evaporation      ( -- )
  REAL(wp), DIMENSION(nproma,nsfc_subs), INTENT(IN) :: &
                  sobs_ex          , & ! solar radiation at the ground                 ( W/m2)
                  thbs_ex          , & ! thermal radiation at the ground               ( W/m2)
                  pabs_ex              !!!! photosynthetic active radiation            ( W/m2)
  REAL(wp), DIMENSION(nproma,nsfc_subs), INTENT(INOUT) :: &
                  runoff_s_ex      , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g_ex          ! soil water runoff; sum over forecast          (kg/m2)
  REAL(wp), DIMENSION(nproma), INTENT(INOUT) :: &
                  t_g              , &
                  qv_s

  TYPE(t_patch), TARGET, INTENT(in) :: p_patch         !< grid/patch info
  TYPE(t_external_data), INTENT(in) :: ext_data        !< external data

  TYPE(t_tiles)                     :: p_tiles(nsubs1) !< tiles structure


    ! Local array bounds:
    !
    INTEGER :: isubs           

    ! Local scalars:
    !
    INTEGER :: jc,jg,jk      !loop indices

    REAL(wp) :: ps_t       (nproma)
    REAL(wp) :: prr_con_t  (nproma)
    REAL(wp) :: prs_con_t  (nproma)
    REAL(wp) :: prr_gsp_t  (nproma)
    REAL(wp) :: prs_gsp_t  (nproma)

    REAL(wp) :: u_t (nproma)
    REAL(wp) :: v_t (nproma)
    REAL(wp) :: t_t (nproma)
    REAL(wp) :: qv_t(nproma)
    REAL(wp) :: p0_t(nproma)

    REAL(wp) :: t_snow_now_t (nproma, nsfc_subs)
    REAL(wp) :: t_snow_new_t (nproma, nsfc_subs)

    REAL(wp) :: t_s_now_t  (nproma, nsfc_subs)
    REAL(wp) :: t_s_new_t  (nproma, nsfc_subs)

    REAL(wp) :: t_g_t      (nproma, nsfc_subs)
    REAL(wp) :: qv_s_t     (nproma, nsfc_subs)

    REAL(wp) :: w_snow_now_t (nproma, nsfc_subs)
    REAL(wp) :: w_snow_new_t (nproma, nsfc_subs)
  
    REAL(wp) :: rho_snow_now_t (nproma, nsfc_subs)
    REAL(wp) :: rho_snow_new_t (nproma, nsfc_subs)

    REAL(wp) :: h_snow_t   (nproma, nsfc_subs)

    REAL(wp) :: w_i_now_t  (nproma, nsfc_subs)
    REAL(wp) :: w_i_new_t  (nproma, nsfc_subs)

    REAL(wp) :: t_2m_t     (nproma, nsfc_subs)
    REAL(wp) :: u_10m_t    (nproma, nsfc_subs)
    REAL(wp) :: v_10m_t    (nproma, nsfc_subs)
    REAL(wp) :: freshsnow_t(nproma, nsfc_subs)

    REAL(wp) :: tch_t      (nproma, nsfc_subs)
    REAL(wp) :: tcm_t      (nproma, nsfc_subs)
    REAL(wp) :: tfv_t      (nproma, nsfc_subs)

    REAL(wp) :: sobs_t     (nproma, nsfc_subs)
    REAL(wp) :: thbs_t     (nproma, nsfc_subs)
    REAL(wp) :: pabs_t     (nproma, nsfc_subs)

    REAL(wp) :: runoff_s_t (nproma, nsfc_subs)
    REAL(wp) :: runoff_g_t (nproma, nsfc_subs)

    ! local dummy variable for precipitation rate of graupel, grid-scale
    REAL(wp) :: dummy_prg_gsp(nproma)

    REAL(wp) :: t_snow_mult_now_t(nproma, nlev_snow+1, nsfc_subs)
    REAL(wp) :: t_snow_mult_new_t(nproma, nlev_snow+1, nsfc_subs)

    REAL(wp) :: rho_snow_mult_now_t(nproma, nlev_snow, nsfc_subs)
    REAL(wp) :: rho_snow_mult_new_t(nproma, nlev_snow, nsfc_subs)

    REAL(wp) :: wliq_snow_now_t(nproma, nlev_snow, nsfc_subs)
    REAL(wp) :: wliq_snow_new_t(nproma, nlev_snow, nsfc_subs)

    REAL(wp) :: wtot_snow_now_t(nproma, nlev_snow, nsfc_subs)
    REAL(wp) :: wtot_snow_new_t(nproma, nlev_snow, nsfc_subs)

    REAL(wp) :: dzh_snow_now_t(nproma, nlev_snow, nsfc_subs)
    REAL(wp) :: dzh_snow_new_t(nproma, nlev_snow, nsfc_subs)

    REAL(wp) :: t_so_now_t(nproma, nlev_soil+2, nsfc_subs)
    REAL(wp) :: t_so_new_t(nproma, nlev_soil+2, nsfc_subs)

    REAL(wp) :: w_so_now_t(nproma, nlev_soil+1, nsfc_subs)
    REAL(wp) :: w_so_new_t(nproma, nlev_soil+1, nsfc_subs)

    REAL(wp) :: w_so_ice_now_t(nproma, nlev_soil+1, nsfc_subs)
    REAL(wp) :: w_so_ice_new_t(nproma, nlev_soil+1, nsfc_subs)

    REAL(wp) :: subsfrac_t (nproma, nsfc_subs)
    INTEGER  :: i_count, ic

    REAL(wp) :: t_g_s(nproma), qv_s_s(nproma)
    REAL(wp) :: qv_2m_t     (nproma, nsfc_subs)
    REAL(wp) :: td_2m_t     (nproma, nsfc_subs)
    REAL(wp) :: rh_2m_t     (nproma, nsfc_subs)
    REAL(wp) :: shfl_s_t    (nproma, nsfc_subs)
    REAL(wp) :: lhfl_s_t    (nproma, nsfc_subs)

    REAL(wp) :: dummy1(nproma)

!--------------------------------------------------------------


      ! initialize dummy variable
      dummy_prg_gsp(1:nproma) = 0._wp
  
      ! local variables related to the blocking
  
      jg        = p_patch%id
    
      IF (msg_level >= 15) THEN
        CALL message('mo_nwp_sfc_interface: ', 'call land-surface scheme')
      ENDIF


      IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
        ! check dry case
        IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
          DO jc = i_startidx, i_endidx
            qv_s (jc) = 0._wp
          ENDDO
        ELSE
          ! 
          !> adjust  humidity at water surface because of changed surface pressure
          !
          DO jc = i_startidx, i_endidx
            qv_s (jc) = &
                 &         spec_humi(sat_pres_water(t_g(jc)),&
                 &                                   ps_ex(jc) )
          ENDDO
        ENDIF
      ENDIF

      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN

!---------- Copy input fields for each tile

!----------------------------------
      DO isubs = 1,nsfc_subs
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)

          ps_t(ic)           =  ps_ex     (jc)    
          prr_con_t(ic)      =  prr_con_ex(jc) 
          prs_con_t(ic)      =  prs_con_ex(jc) 
          prr_gsp_t(ic)      =  prr_gsp_ex(jc) 
          prs_gsp_t(ic)      =  prs_gsp_ex(jc) 

          u_t(ic)      =  u_ex (jc)     
          v_t(ic)      =  v_ex (jc)     
          t_t(ic)      =  t_ex (jc)     
          qv_t(ic)     =  qv_ex(jc) 
          p0_t(ic)     =  p0_ex(jc)     

          t_snow_now_t  (ic,isubs)   =  t_snow_ex  (jc,isubs) 
          t_s_now_t     (ic,isubs)   =  t_s_ex     (jc,isubs)   
          t_g_t         (ic,isubs)   =  t_g_ex     (jc,isubs)
          qv_s_t        (ic,isubs)   =  qv_s_ex    (jc,isubs)  
          w_snow_now_t  (ic,isubs)   =  w_snow_ex  (jc,isubs)
          rho_snow_now_t(ic,isubs)   =  rho_snow_ex(jc,isubs)
          w_i_now_t     (ic,isubs)   =  w_i_ex     (jc,isubs)
          freshsnow_t   (ic,isubs)   =  freshsnow_ex(jc,isubs)
          subsfrac_t    (ic,isubs)   =  subsfrac_ex(jc,isubs) 
          runoff_s_t    (ic,isubs)   =  runoff_s_ex(jc,isubs) 
          runoff_g_t    (ic,isubs)   =  runoff_g_ex(jc,isubs)

          t_2m_t(ic,isubs)     =  t_2m_ex (jc) 
          u_10m_t(ic,isubs)    =  u_10m_ex(jc)
          v_10m_t(ic,isubs)    =  v_10m_ex(jc)  
          tch_t(ic,isubs)      =  tch_ex(jc,isubs)
          tcm_t(ic,isubs)      =  tcm_ex(jc,isubs)
          tfv_t(ic,isubs)      =  tfv_ex(jc,isubs)
          sobs_t(ic,isubs)     =  sobs_ex(jc,isubs) 
          thbs_t(ic,isubs)     =  thbs_ex(jc,isubs) 
          pabs_t(ic,isubs)     =  pabs_ex(jc,isubs) 

          t_so_now_t(ic,nlev_soil+2,isubs) = t_so_ex(jc,nlev_soil+2,isubs)

          IF(lmulti_snow) THEN
            t_snow_mult_now_t(ic,nlev_snow+1,isubs) = &
               t_snow_mult_ex(jc,nlev_snow+1,isubs)
            h_snow_t(ic,isubs)  =  h_snow_ex(jc,isubs)
          ENDIF
        ENDDO

       MSNOWI: IF(lmulti_snow) THEN

#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_snow
#else
!CDIR UNROLL=nlsnow
        DO jk=1,nlev_snow
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_snow_mult_now_t  (ic,jk,isubs) = t_snow_mult_ex  (jc,jk,isubs) 
            rho_snow_mult_now_t(ic,jk,isubs) = rho_snow_mult_ex(jc,jk,isubs)
            wliq_snow_now_t    (ic,jk,isubs) = wliq_snow_ex    (jc,jk,isubs) 
            wtot_snow_now_t    (ic,jk,isubs) = wtot_snow_ex    (jc,jk,isubs)
            dzh_snow_now_t     (ic,jk,isubs) = dzh_snow_ex     (jc,jk,isubs) 
          ENDDO
        ENDDO

       END IF MSNOWI

#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count   
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_soil+1
#else
!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_so_now_t(ic,jk,isubs)     = t_so_ex(jc,jk,isubs) 
            w_so_now_t(ic,jk,isubs)     = w_so_ex(jc,jk,isubs) 
            w_so_ice_now_t(ic,jk,isubs) = w_so_ice_ex(jc,jk,isubs)
          ENDDO
        ENDDO
!
!---------- END Copy index list fields


        CALL terra_multlay(                                    &
        &  ie=nproma                                         , & ! array dimensions
        &  istartpar=1,       iendpar=i_count                , & ! optional start/end indicies
        &  nsubs0=1,          nsubs1=nsfc_subs               , & ! nsfc_subs
        &  ke_soil=nlev_soil, ke_snow=nlev_snow              , &
        &  czmls=zml_soil                                    , & ! processing soil level structure 
        &  dt=dt                                             , &

        &  soiltyp_subs = ext_data%atm%soiltyp_t(:,jb,isubs) , & ! type of the soil (keys 0-9)  --
        &  plcov        = ext_data%atm%plcov_t(:,jb,isubs)   , & ! fraction of plant cover      --
        &  rootdp       = ext_data%atm%rootdp_t(:,jb,isubs)  , & ! depth of the roots         ( m  )
        &  sai          = ext_data%atm%sai_t(:,jb,isubs)     , & ! surface area index           --
        &  tai          = ext_data%atm%tai_t(:,jb,isubs)     , & ! surface area index           --
        &  eai          = ext_data%atm%eai_t(:,jb,isubs)     , & ! surface area index           --
        &  rsmin2d      = ext_data%atm%rsmin2d_t(:,jb,isubs) , & ! minimum stomata resistance ( s/m )
!
        &  u  =  u_t(:)                                      , & ! zonal wind speed
        &  v  =  v_t(:)                                      , & ! meridional wind speed 
        &  t  =  t_t(:)                                      , & ! temperature                            (  k  )
        &  qv =  qv_t(:)                                     , & ! specific water vapor content           (kg/kg)
        &  p0 =  p0_t(:)                                     , & ! base state pressure               ( Pa  ) 
        &  ps =  ps_t(:)                                     , & ! surface pressure                       ( pa  )
!                                                            
        &  t_snow_now    = t_snow_now_t(:,isubs)             , & ! temperature of the snow-surface   (  K  )
        &  t_snow_new    = t_snow_new_t(:,isubs)             , & ! temperature of the snow-surface   (  K  )
!                                                            
        &  t_snow_mult_now = t_snow_mult_now_t(:,:,isubs)    , & ! temperature of the snow-surface (  K  )
        &  t_snow_mult_new = t_snow_mult_new_t(:,:,isubs)    , & ! temperature of the snow-surface (  K  )
!                                                            
        &  t_s_now       = t_s_now_t(:,isubs)                , & ! temperature of the ground surface            (  K  )
        &  t_s_new       = t_s_new_t(:,isubs)                , & ! temperature of the ground surface            (  K  )
!                                                            
        &  t_g           =  t_g_t (:,isubs)                  , & ! weighted surface temperature                 (  K  )
        &  qv_s          =  qv_s_t(:,isubs)                  , & ! specific humidity at the surface             (kg/kg)
!                                                            
        &  w_snow_now    = w_snow_now_t(:,isubs)             , & ! water content of snow      (m H2O) 
        &  w_snow_new    = w_snow_new_t(:,isubs)             , & ! water content of snow      (m H2O) 
!                                                            
        &  rho_snow_now  = rho_snow_now_t(:,isubs)           , & ! snow density            (kg/m**3)
        &  rho_snow_new  = rho_snow_new_t(:,isubs)           , & ! snow density            (kg/m**3)
!
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,isubs), & ! snow density       (kg/m**3) 
        &  rho_snow_mult_new = rho_snow_mult_new_t(:,:,isubs), & ! snow density       (kg/m**3) 
!
        &  h_snow        =  h_snow_t(:,isubs)                , & ! snow height
!                                                            
        &  w_i_now       =   w_i_now_t(:,isubs)              , & ! water content of interception water    (m H2O)
        &  w_i_new       =   w_i_new_t(:,isubs)              , & ! water content of interception water    (m H2O)
!                                                            
        &  t_so_now      = t_so_now_t(:,:,isubs)             , & ! soil temperature (main level)    (  K  )
        &  t_so_new      = t_so_new_t(:,:,isubs)             , & ! soil temperature (main level)    (  K  )
!                                                               
        &  w_so_now      = w_so_now_t(:,:,isubs)             , & ! total water content (ice + liquid water) (m H20)
        &  w_so_new      = w_so_new_t(:,:,isubs)             , & ! total water content (ice + liquid water) (m H20)
!                                                               
        &  w_so_ice_now  = w_so_ice_now_t(:,:,isubs)         , & ! ice content   (m H20)
        &  w_so_ice_new  = w_so_ice_new_t(:,:,isubs)         , & ! ice content   (m H20)
!                                                            
        &  t_2m          =  t_2m_t(:,isubs)                  , & ! ,nsfc_subs, temperature in 2m        (  K  )
        &  u_10m         =  u_10m_t(:,isubs)                 , & ! ,nsfc_subs, zonal wind in 10m       ( m/s )
        &  v_10m         =  v_10m_t(:,isubs)                 , & ! ,nsfc_subs,  meridional wind in 10m ( m/s )
        &  freshsnow     =  freshsnow_t(:,isubs)             , & ! indicator for age of snow in top of snow layer       (  -  )
!                                                            
        &  wliq_snow_now = wliq_snow_now_t(:,:,isubs)        , & ! liquid water content in the snow       (m H2O)
        &  wliq_snow_new = wliq_snow_new_t(:,:,isubs)        , & ! liquid water content in the snow       (m H2O)
!                                                              
        &  wtot_snow_now = wtot_snow_now_t(:,:,isubs)        , & ! total (liquid + solid) water content of snow  (m H2O)
        &  wtot_snow_new = wtot_snow_new_t(:,:,isubs)        , & ! total (liquid + solid) water content of snow  (m H2O)
!                                                              
        &  dzh_snow_now  = dzh_snow_now_t(:,:,isubs)         , & ! layer thickness between half levels in snow   (  m  )
        &  dzh_snow_new  = dzh_snow_new_t(:,:,isubs)         , & ! layer thickness between half levels in snow   (  m  )
!                                                            
        &  subsfrac      =  dummy1(:)                        , & ! 
           !                                                 
        &  prr_con       =  prr_con_t(:)                     , & ! precipitation rate of rain, convective        (kg/m2*s)
        &  prs_con       =  prs_con_t(:)                     , & ! precipitation rate of snow, convective        (kg/m2*s)
        &  prr_gsp       =  prr_gsp_t(:)                     , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
        &  prs_gsp       =  prs_gsp_t(:)                     , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
        &  prg_gsp       =  dummy_prg_gsp(:)                 , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
           !                                                 
        &  tch           = tch_t(:,isubs)                    , & ! ,nsfc_subs, turbulent transfer coefficient for heat     ( -- )
        &  tcm           = tcm_t(:,isubs)                    , & ! ,nsfc_subs, turbulent transfer coefficient for momentum ( -- )
        &  tfv           = tfv_t(:,isubs)                    , & ! ,nsfc_subs, laminar reduction factor for evaporation    ( -- )
           !                                                 
        &  sobs          = sobs_t(:,isubs)                   , & ! ,nsfc_subs, solar radiation at the ground   (W/m2)
        &  thbs          = thbs_t(:,isubs)                   , & ! ,nsfc_subs, thermal radiation at the ground (W/m2)
        &  pabs          = pabs_t(:,isubs)                   , & ! ,nsfc_subs, photosynthetic active radiation (W/m2)
           !                                                 
        &  runoff_s      = runoff_s_t(:,isubs)               , & ! surface water runoff; sum over forecast      (kg/m2)
        &  runoff_g      = runoff_g_t(:,isubs)               , & ! soil water runoff; sum over forecast         (kg/m2)
        &  pt_tiles      = p_tiles(:)                          & ! tiles structure
        &                                                    )


!---------- Copy index list fields back to state fields

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          t_snow_ex   (jc,isubs)  = t_snow_new_t  (ic,isubs)         
          t_s_ex      (jc,isubs)  = t_s_new_t     (ic,isubs)              
          t_g_ex      (jc,isubs)  = t_g_t         (ic,isubs)
          qv_s_ex     (jc,isubs)  = qv_s_t        (ic,isubs)                
          w_snow_ex   (jc,isubs)  = w_snow_new_t  (ic,isubs)          
          rho_snow_ex (jc,isubs)  = rho_snow_new_t(ic,isubs)        
          h_snow_ex   (jc,isubs)  = h_snow_t      (ic,isubs)              
          w_i_ex      (jc,isubs)  = w_i_new_t     (ic,isubs)             
          freshsnow_ex(jc,isubs)  = freshsnow_t   (ic,isubs)   
          subsfrac_ex (jc,isubs)  = subsfrac_t    (ic,isubs)
          runoff_s_ex (jc,isubs)  = runoff_s_t    (ic,isubs)  
          runoff_g_ex (jc,isubs)  = runoff_g_t    (ic,isubs)  

          t_so_ex(jc,nlev_soil+2,isubs) = t_so_new_t(ic,nlev_soil+2,isubs)

          IF(lmulti_snow) THEN
            t_snow_mult_ex     (jc,nlev_snow+1,isubs) = &
              t_snow_mult_new_t(ic,nlev_snow+1,isubs)
          ENDIF
        ENDDO

        MSNOWO: IF(lmulti_snow) THEN

#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_snow
#else
!CDIR UNROLL=nlsnow
        DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_snow_mult_ex  (jc,jk,isubs) = t_snow_mult_new_t  (ic,jk,isubs)   
            rho_snow_mult_ex(jc,jk,isubs) = rho_snow_mult_new_t(ic,jk,isubs) 
            wliq_snow_ex    (jc,jk,isubs) = wliq_snow_new_t    (ic,jk,isubs)     
            wtot_snow_ex    (jc,jk,isubs) = wtot_snow_new_t    (ic,jk,isubs)     
            dzh_snow_ex     (jc,jk,isubs) = dzh_snow_new_t     (ic,jk,isubs)      
          ENDDO
        ENDDO
        
        END IF MSNOWO


#ifdef __LOOP_EXCHANGE
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          DO jk=1,nlev_soil+1
#else
!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_so_ex    (jc,jk,isubs) = t_so_new_t    (ic,jk,isubs)          
            w_so_ex    (jc,jk,isubs) = w_so_new_t    (ic,jk,isubs)          
            w_so_ice_ex(jc,jk,isubs) = w_so_ice_new_t(ic,jk,isubs)     
          ENDDO
        ENDDO

      END DO ! isubs - loop over tiles

      i_count = ext_data%atm%lp_count(jb)

      IF (nsfc_subs == 1) THEN 
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp(ic,jb)
          t_g(jc)  = t_g_ex(jc,1)
          qv_s(jc) = qv_s_ex(jc,1) 
        ENDDO
      ELSE ! aggregate fields over tiles
        t_g_s(:)  =  0._wp
        qv_s_s(:) =  0._wp
        DO isubs = 1,nsfc_subs
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            t_g_s(jc) = t_g_s(jc) + ext_data%atm%lc_frac_t(jc,jb,isubs)* &
              t_g_ex(jc,isubs)
            qv_s_s(jc) = qv_s_s(jc) + ext_data%atm%lc_frac_t(jc,jb,isubs)* &
              qv_s_ex(jc,isubs)
          ENDDO
        ENDDO


         ! Apply relaxation if the grid cell contains water - actually, separate
         ! fields carrying SST and/or lake temperature would be needed for such a weighting
         ! to make really sense! What is done here has an unjustified time-step dependence!
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp(ic,jb)
          t_g(jc)  = t_g_s(jc) ! &
    !     ! This does not work in combination with disaggregating the surface radiation flux terms
    !        (1._wp-ext_data%atm%fr_land(jc,jb))*lnd_prog_now%t_g(jc,jb) + &
    !         ext_data%atm%fr_land(jc,jb)*t_g_s(jc)
          qv_s(jc)     = qv_s_s(jc) ! &
    !        (1._wp-ext_data%atm%fr_land(jc,jb))*qv_s_ex(jc) + &
    !         ext_data%atm%fr_land(jc,jb)*qv_s_s(jc)
        ENDDO

      ENDIF
   
    ELSE IF ( atm_phy_nwp_config(jg)%inwp_surface == 2 ) THEN 

          !-------------------------------------------------------------------------
          !> ECHAM version 
          !-------------------------------------------------------------------------
     
    ENDIF !inwp_sfc


  END SUBROUTINE nwp_surface_terra_edmf


END MODULE mo_nwp_sfc_interface_edmf

