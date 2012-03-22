!>
!! This module is the interface between nwp_nh_interface to the 
!! surface parameterisations:
!! inwp_sfc  == 1 == surface scheme TERRA run in COSMO
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
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
MODULE mo_nwp_sfc_interface

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message !, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, zml_soil
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag 
  USE mo_nwp_phy_state,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_state,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv, msg_level
  USe mo_extpar_config,       ONLY: itopo
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, nsfc_subs, t_tiles,  &
    &                               lseaice, llake, lmulti_snow
  USE mo_satad,               ONLY: sat_pres_water, spec_humi  
  USE mo_soil_ml,             ONLY: terra_multlay, terra_multlay_init
  USE mo_phyparam_soil              ! soil and vegetation parameters for TILES
!  USE mo_aggregate_surface,   ONLY: subsmean,subs_disaggregate_radflux,subsmean_albedo
!  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, igbm, iwtr, iice, ilnd
  
  IMPLICIT NONE 

  PUBLIC  ::  nwp_surface, nwp_surface_init

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
  SUBROUTINE nwp_surface( tcall_sfc_jg,                   & !>in
                        & p_patch,                        & !>in
                        & ext_data,                       & !>in
                        & p_prog_rcf,                     & !>in/inout
                        & p_diag ,                        & !>inout
                        & prm_diag,                       & !>inout 
                        & lnd_prog_now, lnd_prog_new,     & !>inout
                        & lnd_diag,                       & !>inout
                        & p_tiles                         ) !>in

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch       !< grid/patch info
    TYPE(t_external_data),       INTENT(in)   :: ext_data      !< external data
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf    !< call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag        !< diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now  !< prog vars for sfc
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new  !< prog vars for sfc
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: tcall_sfc_jg  !< time interval for 
                                                               !< surface
    TYPE(t_tiles),        TARGET,INTENT(in)   :: p_tiles(:)    !< tiles structure

    ! Local array bounds:
    !
    INTEGER :: nblks_c                 !> number of blocks for cells
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: nlev                    !< number of full levels
    INTEGER :: isubs           

    ! Local scalars:
    !
    INTEGER :: jc,jb,jg,jk      !loop indices


    REAL(wp) ::          ps_t        (nproma, p_patch%nblks_c)
    REAL(wp) ::          prr_con_t   (nproma, p_patch%nblks_c)
    REAL(wp) ::          prs_con_t   (nproma, p_patch%nblks_c)
    REAL(wp) ::          prr_gsp_t   (nproma, p_patch%nblks_c)
    REAL(wp) ::          prs_gsp_t   (nproma, p_patch%nblks_c)

    REAL(wp) ::          u_t (nproma,  p_patch%nblks_c)
    REAL(wp) ::          v_t (nproma,  p_patch%nblks_c)
    REAL(wp) ::          t_t (nproma,  p_patch%nblks_c)
    REAL(wp) ::          qv_t(nproma,  p_patch%nblks_c)
    REAL(wp) ::          p0_t(nproma,  p_patch%nblks_c)

    REAL(wp) ::          t_snow_now_t (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          t_snow_new_t (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          t_s_now_t  (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          t_s_new_t  (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          t_g_t      (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          qv_s_t     (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          w_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          w_snow_new_t(nproma, p_patch%nblks_c, nsfc_subs)
  
    REAL(wp) ::          rho_snow_now_t (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          rho_snow_new_t (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          h_snow_t (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          w_i_now_t (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          w_i_new_t (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          t_2m_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          u_10m_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          v_10m_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          freshsnow_t(nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          tch_t      (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          tcm_t      (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          tfv_t      (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          sobs_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          thbs_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          pabs_t     (nproma, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          runoff_s_t (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          runoff_g_t (nproma, p_patch%nblks_c, nsfc_subs)

    ! local dummy variable for precipitation rate of graupel, grid-scale
    REAL(wp) :: dummy_prg_gsp(nproma)

    REAL(wp) :: t_snow_mult_now_t(nproma, nlev_snow+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_snow_mult_new_t(nproma, nlev_snow+1, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: rho_snow_mult_now_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: rho_snow_mult_new_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: wliq_snow_now_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: wliq_snow_new_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: wtot_snow_now_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: wtot_snow_new_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: dzh_snow_now_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: dzh_snow_new_t(nproma, nlev_snow, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: t_so_now_t(nproma, nlev_soil+2, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_so_new_t(nproma, nlev_soil+2, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: w_so_now_t(nproma, nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_new_t(nproma, nlev_soil+1, p_patch%nblks_c, nsfc_subs)

    REAL(wp) :: w_so_ice_now_t(nproma, nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_ice_new_t(nproma, nlev_soil+1, p_patch%nblks_c, nsfc_subs)

!!$    REAL(wp) :: lu_class_frac(nsfc_subs), sum_frac 
!!$    INTEGER  :: i_tile(nproma, nsfc_subs),lu_subs
    REAL(wp) :: subsfrac_t (nproma, p_patch%nblks_c, nsfc_subs)
    INTEGER  :: i_count, ic

    REAL(wp) :: t_g_s(nproma), qv_s_s(nproma)

!--------------------------------------------------------------


    ! initialize dummy variable
    dummy_prg_gsp(1:nproma) = 0._wp

    ! local variables related to the blocking

    nblks_c   = p_patch%nblks_int_c

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev

    ! in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    IF (msg_level >= 15) THEN
      CALL message('mo_nwp_sfc_interface: ', 'call land-surface scheme')
    ENDIF


!$OMP PARALLEL
#ifdef __xlC__
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,i_count,ic, &
!$OMP            t_g_s,qv_s_s,jk), SCHEDULE(guided)
#else
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,i_count,ic, &
!$OMP            t_g_s,qv_s_s,jk), SCHEDULE(guided)
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
        ! check dry case
        IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
          lnd_diag%qv_s (:,jb) = 0._wp
        ELSE
          ! 
          !> adjust  humidity at water surface because of changed surface pressure
          !
          DO jc = i_startidx, i_endidx
            lnd_diag%qv_s (jc,jb) = &
                 &         spec_humi(sat_pres_water(lnd_prog_now%t_g(jc,jb)),&
                 &                                   p_diag%pres_sfc(jc,jb) )
          ENDDO
        ENDIF
      ENDIF
 
      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN

       IF (ext_data%atm%lp_count(jb) == 0) CYCLE ! skip loop if there is no land point

!---------- Copy input fields for each tile

!----------------------------------
       DO isubs = 1,nsfc_subs
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs) 

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)

          ps_t(ic,jb)           =  p_diag%pres_sfc     (jc,jb)    
          prr_con_t(ic,jb)      =  prm_diag%tracer_rate(jc,jb,3) 
          prs_con_t(ic,jb)      =  prm_diag%tracer_rate(jc,jb,4) 
          prr_gsp_t(ic,jb)      =  prm_diag%tracer_rate(jc,jb,1) 
          prs_gsp_t(ic,jb)      =  prm_diag%tracer_rate(jc,jb,2) 

          u_t(ic,jb)      =  p_diag%u         (jc,nlev,jb)     
          v_t(ic,jb)      =  p_diag%v         (jc,nlev,jb)     
          t_t(ic,jb)      =  p_diag%temp      (jc,nlev,jb)     
          qv_t(ic,jb)     =  p_prog_rcf%tracer(jc,nlev,jb,iqv) 
          p0_t(ic,jb)     =  p_diag%pres      (jc,nlev,jb)     

          t_snow_now_t(ic,jb,isubs)          =  lnd_prog_now%t_snow(jc,jb,isubs) 
          t_s_now_t(ic,jb,isubs)             =  lnd_prog_now%t_s(jc,jb,isubs)   
          t_g_t (ic,jb,isubs)                =  lnd_prog_now%t_gt(jc,jb,isubs)
          qv_s_t(ic,jb,isubs)                =  lnd_diag%qv_st(jc,jb,isubs)  
          w_snow_now_t(ic,jb,isubs)          =  lnd_prog_now%w_snow(jc,jb,isubs)
          rho_snow_now_t(ic,jb,isubs)        =  lnd_prog_now%rho_snow(jc,jb,isubs)
          w_i_now_t(ic,jb,isubs)             =  lnd_prog_now%w_i(jc,jb,isubs)
          freshsnow_t(ic,jb,isubs)           =  lnd_diag%freshsnow(jc,jb,isubs)
          subsfrac_t(ic,jb,isubs)            =  lnd_diag%subsfrac(jc,jb,isubs) 
          freshsnow_t(ic,jb,isubs)           =  lnd_diag%freshsnow(jc,jb,isubs)
          subsfrac_t(ic,jb,isubs)            =  lnd_diag%subsfrac(jc,jb,isubs) 
          runoff_s_t(ic,jb,isubs)            =  lnd_diag%runoff_s(jc,jb,isubs) 
          runoff_g_t(ic,jb,isubs)            =  lnd_diag%runoff_g(jc,jb,isubs)

!GZ: Why do these local fields need to be dimensioned with ntiles?
          t_2m_t(ic,jb,isubs)                =  prm_diag%t_2m(jc,jb) 
          u_10m_t(ic,jb,isubs)               =  prm_diag%u_10m(jc,jb)
          v_10m_t(ic,jb,isubs)               =  prm_diag%v_10m(jc,jb)  
          tch_t(ic,jb,isubs)                 =  prm_diag%tch(jc,jb)
          tcm_t(ic,jb,isubs)                 =  prm_diag%tcm(jc,jb)
          tfv_t(ic,jb,isubs)                 =  prm_diag%tfv(jc,jb)
          sobs_t(ic,jb,isubs)                =  prm_diag%swflxsfc(jc,jb) 
          thbs_t(ic,jb,isubs)                =  prm_diag%lwflxsfc(jc,jb) 
          pabs_t(ic,jb,isubs)                =  prm_diag%swflxsfc(jc,jb) 

          t_so_now_t(ic,nlev_soil+2,jb,isubs) = lnd_prog_now%t_so(jc,nlev_soil+2,jb,isubs)

          IF(lmulti_snow) THEN
            t_snow_mult_now_t(ic,nlev_snow+1,jb,isubs) = &
              lnd_prog_now%t_snow_mult(jc,nlev_snow+1,jb,isubs)
            h_snow_t(ic,jb,isubs)  =  lnd_diag%h_snow(jc,jb,isubs)
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
            t_snow_mult_now_t(ic,jk,jb,isubs) = lnd_prog_now%t_snow_mult(jc,jk,jb,isubs) 
            rho_snow_mult_now_t(ic,jk,jb,isubs) = lnd_prog_now%rho_snow_mult(jc,jk,jb,isubs)
            wliq_snow_now_t(ic,jk,jb,isubs) = lnd_prog_now%wliq_snow(jc,jk,jb,isubs) 
            wtot_snow_now_t(ic,jk,jb,isubs) = lnd_prog_now%wtot_snow(jc,jk,jb,isubs)
            dzh_snow_now_t(ic,jk,jb,isubs) = lnd_prog_now%dzh_snow(jc,jk,jb,isubs) 
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
            t_so_now_t(ic,jk,jb,isubs) =  lnd_prog_now%t_so(jc,jk,jb,isubs) 
            w_so_now_t(ic,jk,jb,isubs)     = lnd_prog_now%w_so(jc,jk,jb,isubs) 
            w_so_ice_now_t(ic,jk,jb,isubs) = lnd_prog_now%w_so_ice(jc,jk,jb,isubs)
          ENDDO
        ENDDO
!
!---------- END Copy index list fields


       CALL terra_multlay(                             &
        &  ie=nproma,                                  & ! array dimensions
        &  istartpar=1, iendpar=i_count    , & ! optional start/end indicies
        &  nsubs0=1, nsubs1=nsfc_subs       , & ! nsfc_subs
        &  ke_soil=nlev_soil, ke_snow=nlev_snow      , &
        &  czmls=zml_soil                            , & ! processing soil level structure 
        &  dt=tcall_sfc_jg                           , &
        &  soiltyp_subs  = ext_data%atm%soiltyp_t(:,jb,isubs) , & ! type of the soil (keys 0-9)  --
        &  plcov         = ext_data%atm%plcov_t(:,jb,isubs)   , & ! fraction of plant cover      --
        &  rootdp        = ext_data%atm%rootdp_t(:,jb,isubs)  , & ! depth of the roots         ( m  )
        &  sai           = ext_data%atm%sai_t(:,jb,isubs)     , & ! surface area index           --
        &  tai           = ext_data%atm%tai_t(:,jb,isubs)     , & ! surface area index           --
        &  eai           = ext_data%atm%eai_t(:,jb,isubs)     , & ! surface area index           --
        &  rsmin2d       = ext_data%atm%rsmin2d_t(:,jb,isubs) , & ! minimum stomata resistance ( s/m )
!
        &  u  =  u_t(:,jb)      , & ! zonal wind speed
        &  v  =  v_t(:,jb)      , & ! meridional wind speed 
        &  t  =  t_t(:,jb)      , & ! temperature                            (  k  )
        &  qv =  qv_t(:,jb)     , & ! specific water vapor content           (kg/kg)
        &  p0 =  p0_t(:,jb)     , & ! base state pressure               ( Pa  ) 
        &  ps =  ps_t(:,jb)     , & ! surface pressure                       ( pa  )
!
        &  t_snow_now    = t_snow_now_t(:,jb,isubs)     , & ! temperature of the snow-surface   (  K  )
        &  t_snow_new    = t_snow_new_t(:,jb,isubs)     , & ! temperature of the snow-surface   (  K  )
!
        &  t_snow_mult_now   = t_snow_mult_now_t(:,:,jb,isubs) ,& ! temperature of the snow-surface (  K  )
        &  t_snow_mult_new   = t_snow_mult_new_t(:,:,jb,isubs) ,& ! temperature of the snow-surface (  K  )
!
        &  t_s_now           = t_s_now_t(:,jb,isubs)   , & ! temperature of the ground surface            (  K  )
        &  t_s_new           = t_s_new_t(:,jb,isubs)   , & ! temperature of the ground surface            (  K  )
!
        &  t_g           =  t_g_t (:,jb,isubs)     , & ! weighted surface temperature                 (  K  )
        &  qv_s          =  qv_s_t(:,jb,isubs)     , & ! specific humidity at the surface             (kg/kg)
!
        &  w_snow_now    = w_snow_now_t(:,jb,isubs)  , & ! water content of snow      (m H2O) 
        &  w_snow_new    = w_snow_new_t(:,jb,isubs)  , & ! water content of snow      (m H2O) 
!
        &  rho_snow_now      = rho_snow_now_t(:,jb,isubs)   , & ! snow density            (kg/m**3)
        &  rho_snow_new      = rho_snow_new_t(:,jb,isubs)   , & ! snow density            (kg/m**3)
!
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,jb,isubs), & ! snow density       (kg/m**3) 
        &  rho_snow_mult_new = rho_snow_mult_new_t(:,:,jb,isubs), & ! snow density       (kg/m**3) 
!
        &  h_snow        =  h_snow_t(:,jb,isubs)   , & ! snow height
!
        &  w_i_now       =   w_i_now_t(:,jb,isubs)   , & ! water content of interception water    (m H2O)
        &  w_i_new       =   w_i_new_t(:,jb,isubs)   , & ! water content of interception water    (m H2O)
!
        &  t_so_now          = t_so_now_t(:,:,jb,isubs)  , & ! soil temperature (main level)    (  K  )
        &  t_so_new          = t_so_new_t(:,:,jb,isubs)  , & ! soil temperature (main level)    (  K  )
!
        &  w_so_now          = w_so_now_t(:,:,jb,isubs)  , & ! total water content (ice + liquid water) (m H20)
        &  w_so_new          = w_so_new_t(:,:,jb,isubs)  , & ! total water content (ice + liquid water) (m H20)
!
        &  w_so_ice_now      = w_so_ice_now_t(:,:,jb,isubs)   , & ! ice content   (m H20)
        &  w_so_ice_new      = w_so_ice_new_t(:,:,jb,isubs)   , & ! ice content   (m H20)
!
        &  t_2m          =  t_2m_t(:,jb,isubs),&   ! ,nsfc_subs, temperature in 2m        (  K  )
        &  u_10m         =  u_10m_t(:,jb,isubs),&  ! ,nsfc_subs, zonal wind in 10m       ( m/s )
        &  v_10m         =  v_10m_t(:,jb,isubs),&  ! ,nsfc_subs,  meridional wind in 10m ( m/s )
        &  freshsnow     =  freshsnow_t(:,jb,isubs) , & ! indicator for age of snow in top of snow layer       (  -  )
!
        &  wliq_snow_now     = wliq_snow_now_t(:,:,jb,isubs) , & ! liquid water content in the snow       (m H2O)
        &  wliq_snow_new     = wliq_snow_new_t(:,:,jb,isubs) , & ! liquid water content in the snow       (m H2O)
!
        &  wtot_snow_now     = wtot_snow_now_t(:,:,jb,isubs) , & ! total (liquid + solid) water content of snow  (m H2O)
        &  wtot_snow_new     = wtot_snow_new_t(:,:,jb,isubs) , & ! total (liquid + solid) water content of snow  (m H2O)
!
        &  dzh_snow_now      = dzh_snow_now_t(:,:,jb,isubs) , & ! layer thickness between half levels in snow   (  m  )
        &  dzh_snow_new      = dzh_snow_new_t(:,:,jb,isubs) , & ! layer thickness between half levels in snow   (  m  )
!
        &  subsfrac      =  subsfrac_t(:,jb,isubs)  , & ! 
           !
        &  prr_con       =  prr_con_t(:,jb)  , & ! precipitation rate of rain, convective        (kg/m2*s)
        &  prs_con       =  prs_con_t(:,jb)  , & ! precipitation rate of snow, convective        (kg/m2*s)
        &  prr_gsp       =  prr_gsp_t(:,jb)  , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
        &  prs_gsp       =  prs_gsp_t(:,jb)  , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
        &  prg_gsp       =  dummy_prg_gsp(:)              , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
           !
        &  tch           = tch_t(:,jb,isubs),&  ! ,nsfc_subs, turbulent transfer coefficient for heat     ( -- )
        &  tcm           = tcm_t(:,jb,isubs),&  ! ,nsfc_subs, turbulent transfer coefficient for momentum ( -- )
        &  tfv           = tfv_t(:,jb,isubs),&  ! ,nsfc_subs, laminar reduction factor for evaporation    ( -- )
           ! 
        &  sobs          = sobs_t(:,jb,isubs) ,&  ! ,nsfc_subs, solar radiation at the ground   (W/m2)
        &  thbs          = thbs_t(:,jb,isubs) ,&  ! ,nsfc_subs, thermal radiation at the ground (W/m2)
        &  pabs          = pabs_t(:,jb,isubs) ,&  ! ,nsfc_subs, photosynthetic active radiation (W/m2)
           ! 
        &  runoff_s      = runoff_s_t(:,jb,isubs)   , & ! surface water runoff; sum over forecast      (kg/m2)
        &  runoff_g      = runoff_s_t(:,jb,isubs)   , & ! soil water runoff; sum over forecast         (kg/m2)
        &  pt_tiles      = p_tiles(:)                   & ! tiles structure
        &                                           )


!---------- Copy index list fields back to state fields

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          lnd_prog_new%t_snow(jc,jb,isubs) = t_snow_new_t(ic,jb,isubs)         
          lnd_prog_new%t_s(jc,jb,isubs)  = t_s_new_t(ic,jb,isubs)              
          lnd_prog_new%t_gt(jc,jb,isubs)  = t_g_t (ic,jb,isubs)
          lnd_diag%qv_st(jc,jb,isubs)     = qv_s_t(ic,jb,isubs)                
          lnd_prog_new%w_snow(jc,jb,isubs)  = w_snow_new_t(ic,jb,isubs)          
          lnd_prog_new%rho_snow(jc,jb,isubs)  = rho_snow_new_t(ic,jb,isubs)        
          lnd_diag%h_snow(jc,jb,isubs)       = h_snow_t(ic,jb,isubs)              
          lnd_prog_new%w_i(jc,jb,isubs)  = w_i_new_t(ic,jb,isubs)             
          lnd_diag%freshsnow(jc,jb,isubs) = freshsnow_t(ic,jb,isubs)   
          lnd_diag%subsfrac(jc,jb,isubs)  = subsfrac_t(ic,jb,isubs)
          lnd_diag%runoff_s(jc,jb,isubs) = runoff_s_t(ic,jb,isubs)  
          lnd_diag%runoff_g(jc,jb,isubs) = runoff_g_t(ic,jb,isubs)  

          lnd_prog_new%t_so(jc,nlev_soil+2,jb,isubs) = t_so_new_t(ic,nlev_soil+2,jb,isubs)

          IF(lmulti_snow) THEN
            lnd_prog_new%t_snow_mult(jc,nlev_snow+1,jb,isubs) = &
              t_snow_mult_new_t(ic,nlev_snow+1,jb,isubs)
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
            lnd_prog_new%t_snow_mult(jc,jk,jb,isubs) = t_snow_mult_new_t(ic,jk,jb,isubs)   
            lnd_prog_new%rho_snow_mult(jc,jk,jb,isubs) = rho_snow_mult_new_t(ic,jk,jb,isubs) 
            lnd_prog_new%wliq_snow(jc,jk,jb,isubs) = wliq_snow_new_t(ic,jk,jb,isubs)     
            lnd_prog_new%wtot_snow(jc,jk,jb,isubs) = wtot_snow_new_t(ic,jk,jb,isubs)     
            lnd_prog_new%dzh_snow(jc,jk,jb,isubs)  = dzh_snow_new_t(ic,jk,jb,isubs)      
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
            lnd_prog_new%t_so(jc,jk,jb,isubs) = t_so_new_t(ic,jk,jb,isubs)          
            lnd_prog_new%w_so(jc,jk,jb,isubs) = w_so_new_t(ic,jk,jb,isubs)          
            lnd_prog_new%w_so_ice(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)     
          ENDDO
        ENDDO

       END DO ! isubs - loop over tiles

       i_count = ext_data%atm%lp_count(jb)

       IF (nsfc_subs == 1) THEN 
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp(ic,jb)
           lnd_prog_new%t_g(jc,jb)  = lnd_prog_new%t_gt(jc,jb,1)
           lnd_diag%qv_s(jc,jb)     = lnd_diag%qv_st(jc,jb,1) 
         ENDDO
       ELSE ! aggregate fields over tiles
         t_g_s(:)  =  0._wp
         qv_s_s(:) =  0._wp
         DO isubs = 1,nsfc_subs
!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_lp(ic,jb)
             t_g_s(jc) = t_g_s(jc) + ext_data%atm%lc_frac_t(jc,jb,isubs)* &
               lnd_prog_new%t_gt(jc,jb,isubs)
             qv_s_s(jc) = qv_s_s(jc) + ext_data%atm%lc_frac_t(jc,jb,isubs)* &
               lnd_diag%qv_st(jc,jb,isubs)
           ENDDO
         ENDDO

         ! Apply relaxation if the grid cell contains water - actually, separate
         ! fields carrying SST and/or lake temperature would be needed for such a weighting
         ! to make really sense! What is done here has an unjustified time-step dependence!
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp(ic,jb)
           lnd_prog_new%t_g(jc,jb)  = &
             (1._wp-ext_data%atm%fr_land(jc,jb))*lnd_prog_now%t_g(jc,jb) + &
              ext_data%atm%fr_land(jc,jb)*t_g_s(jc)
           lnd_diag%qv_s(jc,jb)     = &
             (1._wp-ext_data%atm%fr_land(jc,jb))*lnd_diag%qv_s(jc,jb) + &
              ext_data%atm%fr_land(jc,jb)*qv_s_s(jc)
         ENDDO

       ENDIF
   
    
      ELSE IF ( atm_phy_nwp_config(jg)%inwp_surface == 2 ) THEN 

          !-------------------------------------------------------------------------
          !> ECHAM version 
          !-------------------------------------------------------------------------
     

     
    ENDIF !inwp_sfc

  ENDDO    
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE nwp_surface

  !-------------------------------------------------------------------------
  !>
  !! Init surface model TERRA
  !!
  !! Init surface model TERRA.
  !!
  !! @par Revision History
  !! Initial revision by Ekaterina Machulskaya, DWD (2011-07-??)
  !! Modification by Daniel Reienrt, DWD (2011-07-29)
  !! - initialize climatological layer t_so(nlev_soil+2)
  !!
  SUBROUTINE nwp_surface_init( p_patch, ext_data, p_prog_lnd_now, &
    &                          p_prog_lnd_new )
 
                             

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(IN)    :: ext_data
!    TYPE(t_tiles)        , INTENT(INOUT) :: p_tiles(:)
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
!!$    REAL(wp)             , INTENT(IN)   :: subsfrac(nproma,1,nsfc_subs)
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jc,jb,isubs,jk


    REAL(wp) :: t_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_snow_mult_now_t(nproma, 1:nlev_snow+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_s_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_s_new_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: rho_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: rho_snow_mult_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_so_now_t(nproma, 1:nlev_soil+2, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_so_new_t(nproma, 1:nlev_soil+2, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_ice_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_ice_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: wliq_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: wtot_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: dzh_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)

!    INTEGER  :: i_tile(nproma,nsfc_subs),lu_subs
    INTEGER  :: i_count, ic

  !-------------------------------------------------------------------------


    i_nchdom  = MAX(1,p_patch%n_childdom)

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,i_count,ic), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF (itopo == 1) THEN
        DO isubs = 1, nsfc_subs
          DO jc = i_startidx, i_endidx

            ! initialize climatological layer (deepest layer of t_so)
            p_prog_lnd_now%t_so(jc,nlev_soil+2,jb,isubs) = ext_data%atm%t_cl(jc,jb)
            p_prog_lnd_new%t_so(jc,nlev_soil+2,jb,isubs) = ext_data%atm%t_cl(jc,jb)

            p_prog_lnd_now%t_gt(jc,jb,isubs) = p_prog_lnd_now%t_g(jc,jb)
            p_prog_lnd_new%t_gt(jc,jb,isubs) = p_prog_lnd_now%t_g(jc,jb)

          END DO 
        END DO
      ENDIF



       IF (ext_data%atm%lp_count(jb) == 0) CYCLE ! skip loop if there is no land point

!---------- Copy input fields for each tile

!----------------------------------
       DO isubs = 1,nsfc_subs
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs) 

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          t_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%t_snow(jc,jb,isubs) 
          t_s_now_t(ic,jb,isubs)             =  p_prog_lnd_now%t_s(jc,jb,isubs)   
          t_s_new_t(ic,jb,isubs)             =  p_prog_lnd_new%t_s(jc,jb,isubs)   
          w_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%w_snow(jc,jb,isubs)  
          rho_snow_now_t(ic,jb,isubs)        =  p_prog_lnd_now%rho_snow(jc,jb,isubs)
        ENDDO


       IMSNOWI: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
        DO jk=1,nlev_snow+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            t_snow_mult_now_t(ic,jk,jb,isubs)   =  p_prog_lnd_now%t_snow_mult(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

!CDIR UNROLL=nlsnow
        DO jk=1,nlev_snow
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            rho_snow_mult_now_t(ic,jk,jb,isubs) =  p_prog_lnd_now%rho_snow_mult(jc,jk,jb,isubs)
            wliq_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wliq_snow(jc,jk,jb,isubs) 
            wtot_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wtot_snow(jc,jk,jb,isubs) 
            dzh_snow_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%dzh_snow(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

     END IF  IMSNOWI

!CDIR UNROLL=nlsoil+2
        DO jk=1,nlev_soil+2
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            t_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%t_so(jc,jk,jb,isubs) 
            t_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%t_so(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            w_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%w_so(jc,jk,jb,isubs) 
            w_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%w_so(jc,jk,jb,isubs) 
            w_so_ice_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%w_so_ice(jc,jk,jb,isubs) 
            w_so_ice_new_t(ic,jk,jb,isubs)      =  p_prog_lnd_new%w_so_ice(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

         
        CALL terra_multlay_init(                            &
        &  ie=nproma,                                  & ! array dimensions
        &  istartpar=1, iendpar= i_count, & ! optional start/end indicies
!        &  ke=nlev, &! nsubs0=1, nsubs1=nsfc_subs       , & ! nsfc_subs
        &  ke_soil=nlev_soil, ke_snow=nlev_snow      , &
        &  czmls=zml_soil                            , & ! processing soil level structure 
        &  soiltyp_subs  =  ext_data%atm%soiltyp_t(:,jb,isubs)  , & ! type of the soil (keys 0-9)  --
        &  rootdp        =  ext_data%atm%rootdp_t(:,jb,isubs)   , & ! depth of the roots         ( m  )
        &  t_snow_now    =  t_snow_now_t(:,jb,isubs)   , & ! temperature of the snow-surface   (  K  )
        &  t_snow_mult_now   = t_snow_mult_now_t(:,:,jb,isubs) ,& ! temperature of the snow-surface (  K  )
        &  t_s_now           = t_s_now_t(:,jb,isubs)    , & ! temperature of the ground surface            (  K  )
        &  t_s_new           = t_s_new_t(:,jb,isubs)    , & ! temperature of the ground surface            (  K  )
        &  w_snow_now        = w_snow_now_t(:,jb,isubs)    , & ! water content of snow                        (m H2O)
        &  rho_snow_now      = rho_snow_now_t(:,jb,isubs)    , & ! snow density            (kg/m**3)
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,jb,isubs)    , & ! snow density       (kg/m**3)
        &  t_so_now          = t_so_now_t(:,:,jb,isubs)    , & ! soil temperature (main level)    (  K  )
        &  t_so_new          = t_so_new_t(:,:,jb,isubs)    , & ! soil temperature (main level)    (  K  )
        &  w_so_now          = w_so_now_t(:,:,jb,isubs)    , & ! total water content (ice + liquid water) (m H20)
        &  w_so_new          = w_so_new_t(:,:,jb,isubs)    , & ! total water content (ice + liquid water) (m H20)
        &  w_so_ice_now      = w_so_ice_now_t(:,:,jb,isubs)    , & ! ice content   (m H20)
        &  w_so_ice_new      = w_so_ice_new_t(:,:,jb,isubs)    , & ! ice content   (m H20)
        &  wliq_snow_now     = wliq_snow_now_t(:,:,jb,isubs)    , & ! liquid water content in the snow       (m H2O)
        &  wtot_snow_now     = wtot_snow_now_t(:,:,jb,isubs)    , & ! total (liquid + solid) water content of snow  (m H2O)
        &  dzh_snow_now      = dzh_snow_now_t(:,:,jb,isubs)       & ! layer thickness between half levels in snow   (  m  )
                                                      )

!  Recover fields from index list
!
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          p_prog_lnd_now%t_snow(jc,jb,isubs)   = t_snow_now_t(ic,jb,isubs)
          p_prog_lnd_now%t_s(jc,jb,isubs)      = t_s_now_t(ic,jb,isubs)  
          p_prog_lnd_new%t_s(jc,jb,isubs)      = t_s_new_t(ic,jb,isubs) 
          p_prog_lnd_now%w_snow(jc,jb,isubs)   = w_snow_now_t(ic,jb,isubs) 
          p_prog_lnd_now%rho_snow(jc,jb,isubs) = rho_snow_now_t(ic,jb,isubs)
        ENDDO

         IMSNOWO: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
        DO jk=1,nlev_snow+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            p_prog_lnd_now%t_snow_mult(jc,jk,jb,isubs) =  t_snow_mult_now_t(ic,jk,jb,isubs)   
          ENDDO
        ENDDO

!CDIR UNROLL=nlsnow
        DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            p_prog_lnd_now%rho_snow_mult(jc,jk,jb,isubs) = rho_snow_mult_now_t(ic,jk,jb,isubs) 
            p_prog_lnd_now%wliq_snow(jc,jk,jb,isubs) = wliq_snow_now_t(ic,jk,jb,isubs)   
            p_prog_lnd_now%wtot_snow(jc,jk,jb,isubs) = wtot_snow_now_t(ic,jk,jb,isubs)
            p_prog_lnd_now%dzh_snow(jc,jk,jb,isubs)  = dzh_snow_now_t(ic,jk,jb,isubs)    
          ENDDO
        ENDDO

        END IF  IMSNOWO

!CDIR UNROLL=nlsoil+2
        DO jk=1,nlev_soil+2
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            p_prog_lnd_now%t_so(jc,jk,jb,isubs) = t_so_now_t(ic,jk,jb,isubs)          
            p_prog_lnd_new%t_so(jc,jk,jb,isubs) = t_so_new_t(ic,jk,jb,isubs)          
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            p_prog_lnd_now%w_so(jc,jk,jb,isubs) = w_so_now_t(ic,jk,jb,isubs)        
            p_prog_lnd_new%w_so(jc,jk,jb,isubs) = w_so_new_t(ic,jk,jb,isubs)        
            p_prog_lnd_now%w_so_ice(jc,jk,jb,isubs) = w_so_ice_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_ice(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO

      END DO ! isubs

    ENDDO  ! jb loop
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_surface_init

END MODULE mo_nwp_sfc_interface

