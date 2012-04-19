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
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_nwp_phy_state,       ONLY: t_nwp_phy_diag
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
  USE mo_physical_constants,  ONLY: rdocp => rd_o_cpd  ! r_d / cp_d
  
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
    REAL(wp) ::          t_g_now_t  (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          t_g_new_t  (nproma, p_patch%nblks_c, nsfc_subs)
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
    REAL(wp) :: t_s_s(nproma), t_snow_s(nproma), t_snow_mult_s(nproma, nlev_snow),          &
      &         w_snow_s(nproma), wtot_snow_s(nproma, nlev_snow), rho_snow_s(nproma),       &
      &         rho_snow_mult_s(nproma, nlev_snow),  h_snow_s(nproma), freshsnow_s(nproma), &
      &         w_i_s(nproma), t_so_s(nproma, nlev_soil+2), w_so_s(nproma, nlev_soil+1),    &
      &         w_so_ice_s(nproma, nlev_soil+1), runoff_s_s(nproma), runoff_g_s(nproma),    &
      &         tch_s(nproma), tfv_s(nproma), t_2m_s(nproma), qv_2m_s(nproma),              &
      &         td_2m_s(nproma), rh_2m_s(nproma), u_10m_s(nproma), v_10m_s(nproma),         &
      &         shfl_s_s(nproma), lhfl_s_s(nproma)
    REAL(wp) ::          qv_2m_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          td_2m_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          rh_2m_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          shfl_s_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          lhfl_s_t    (nproma, p_patch%nblks_c, nsfc_subs)
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
          t_g_now_t (ic,jb,isubs)            =  lnd_prog_now%t_gt(jc,jb,isubs)
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

          t_2m_t(ic,jb,isubs)                =  prm_diag%t_2m(jc,jb) 
          u_10m_t(ic,jb,isubs)               =  prm_diag%u_10m(jc,jb)
          v_10m_t(ic,jb,isubs)               =  prm_diag%v_10m(jc,jb)  
          tch_t(ic,jb,isubs)                 =  prm_diag%tch(jc,jb)
          tcm_t(ic,jb,isubs)                 =  prm_diag%tcm(jc,jb)
          tfv_t(ic,jb,isubs)                 =  prm_diag%tfv(jc,jb)
          sobs_t(ic,jb,isubs)                =  prm_diag%swflxsfc_t(jc,jb,isubs) 
          thbs_t(ic,jb,isubs)                =  prm_diag%lwflxsfc_t(jc,jb,isubs) 
          pabs_t(ic,jb,isubs)                =  prm_diag%swflxsfc_t(jc,jb,isubs) 

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
        &  t_g           =  t_g_t (:,jb,isubs) , & ! weighted surface temperature                 (  K  )
!EM        &  t_g_now       =  t_g_now_t (:,jb,isubs) , & ! weighted surface temperature                 (  K  )
!EM        &  t_g_new       =  t_g_new_t (:,jb,isubs) , & ! weighted surface temperature                 (  K  )
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
          t_g_new_t (ic,jb,isubs) = t_g_t (ic,jb,isubs)
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

!         CALL subsmean(                                        &
!        &  i_count                                           , & ! 
!        &  t_g_s, qv_s_s, t_s_s, t_snow_s, t_snow_mult_s, w_snow_s, wtot_snow_s, &
!        &  rho_snow_s, rho_snow_mult_s, h_snow_s, freshsnow_s, w_i_s, t_so_s,    &
!        &  w_so_s, w_so_ice_s, runoff_s_s, runoff_g_s, tch_s, tfv_s, t_2m_s,     &
!        &  qv_2m_s, td_2m_s, rh_2m_s, u_10m_s, v_10m_s, shfl_s_s, lhfl_s_s ,     &
!        &  ie=nproma                                         , & ! 
!        &  nsubs0=1, nsubs1=nsfc_subs                        , & ! nsfc_subs
!        &  ke_soil=nlev_soil, ke_snow=nlev_snow              , & ! 
!        &  idx_lst_lndp      = ext_data%atm%idx_lst_lp(:,jb) , & ! index list for land points               (  -  )
!        &  subsfrac          = ext_data%atm%lc_frac_t(:,jb,:), & ! fractions of tiles                       (  -  )
!        &  t                 = t_t(:,jb)                     , & ! temperature                              (  K  )
!        &  qv                = qv_t(:,jb)                    , & ! specific water vapor content             (kg/kg)
!        &  p0                = p0_t(:,jb)                    , & ! base state pressure                      ( Pa  ) 
!        &  ps                = ps_t(:,jb)                    , & ! surface pressure                         ( Pa  )
!        &  t_g_now           = t_g_now_t(:,jb,:)             , & ! surface temperature                      (  K  )
!        &  t_g_new           = t_g_new_t(:,jb,:)             , & ! surface temperature                      (  K  )
!        &  qv_s              = qv_s_t(:,jb,:)                , & ! specific humidity at the surface         (kg/kg)
!        &  t_s               = t_s_new_t(:,jb,:)             , & ! temperature of the ground surface        (  K  )
!        &  t_snow            = t_snow_new_t(:,jb,:)          , & ! temperature of the snow-surface          (  K  )
!        &  t_snow_mult       = t_snow_mult_new_t(:,:,jb,:)   , & ! temperature of the snow-surface          (  K  )
!        &  w_snow            = w_snow_new_t(:,jb,:)          , & ! water content of snow                    (m H2O)
!        &  wtot_snow         = wtot_snow_new_t(:,:,jb,:)     , & ! total (liquid + solid) water content of snow  (m H2O)
!        &  rho_snow          = rho_snow_new_t(:,jb,:)        , & ! snow density                             (kg/m**3)
!        &  rho_snow_mult     = rho_snow_mult_new_t(:,:,jb,:) , & ! snow density                             (kg/m**3)
!        &  h_snow            = h_snow_t(:,jb,:)              , & ! snow height                              (m H2O)
!        &  freshsnow         = freshsnow_t(:,jb,:)           , & ! indicator for age of snow in top of snow layer (  -  )
!        &  w_i               = w_i_new_t(:,jb,:)             , & ! water content of interception water      (m H2O)
!        &  t_so              = t_so_new_t(:,:,jb,:)          , & ! soil temperature (main level)            (  K  )
!        &  w_so              = w_so_new_t(:,:,jb,:)          , & ! total water content (ice + liquid water) (m H20)
!        &  w_so_ice          = w_so_ice_new_t(:,:,jb,:)      , & ! soil ice content                         (m H20)
!        &  runoff_s          = runoff_s_t(:,jb,:)            , & ! surface water runoff; sum over forecast  (kg/m2)
!        &  runoff_g          = runoff_s_t(:,jb,:)            , & ! soil water runoff; sum over forecast     (kg/m2)
!        &  tch               = tch_t(:,jb,:)                 , & ! turbulent transfer coefficient for heat  ( -- )
!        &  tfv               = tfv_t(:,jb,:)                 , & ! laminar reduction factor for evaporation ( -- )
!        &  t_2m              = t_2m_t(:,jb,:)                , & ! temperature in 2m                        (  K  )
!        &  qv_2m             = qv_2m_t(:,jb,:)               , & ! humidity in 2m                           (kg/kg)
!        &  td_2m             = td_2m_t(:,jb,:)               , & ! dew point in 2m                          (  K  )
!        &  rh_2m             = rh_2m_t(:,jb,:)               , & ! relative humidity in 2m                  (  K  )
!        &  u_10m             = u_10m_t(:,jb,:)               , & ! zonal wind in 10m                        ( m/s )
!        &  v_10m             = v_10m_t(:,jb,:)               , & ! meridional wind in 10m                   ( m/s )
!        &  shfl_s            = shfl_s_t(:,jb,:)              , & ! sensible heat flux                       (W/m2 )
!        &  lhfl_s            = lhfl_s_t(:,jb,:)              )   ! latent heat flux                         (W/m2 )

         ! Apply relaxation if the grid cell contains water - actually, separate
         ! fields carrying SST and/or lake temperature would be needed for such a weighting
         ! to make really sense! What is done here has an unjustified time-step dependence!
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp(ic,jb)
           lnd_prog_new%t_g(jc,jb)  = t_g_s(jc) ! &
    !      ! This does not work in combination with disaggregating the surface radiation flux terms
    !         (1._wp-ext_data%atm%fr_land(jc,jb))*lnd_prog_now%t_g(jc,jb) + &
    !          ext_data%atm%fr_land(jc,jb)*t_g_s(jc)
           lnd_diag%qv_s(jc,jb)     = qv_s_s(jc) ! &
    !         (1._wp-ext_data%atm%fr_land(jc,jb))*lnd_diag%qv_s(jc,jb) + &
    !          ext_data%atm%fr_land(jc,jb)*qv_s_s(jc)
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

!-------------------------------------------------------------------------

  SUBROUTINE subsmean(i_count, t_g_s, qv_s_s, t_s_s, t_snow_s, t_snow_mult_s, &
        &  w_snow_s, wtot_snow_s, rho_snow_s, rho_snow_mult_s, h_snow_s, freshsnow_s,   &
        &  w_i_s, t_so_s, w_so_s, w_so_ice_s, runoff_s_s, runoff_g_s, tch_s, tfv_s,     &
        &  t_2m_s, qv_2m_s, td_2m_s, rh_2m_s, u_10m_s, v_10m_s, shfl_s_s, lhfl_s_s,     &
        &  ie, nsubs0, nsubs1, ke_soil, ke_snow, idx_lst_lndp,      &
        &  subsfrac, t, qv, p0, ps, t_g_now, t_g_new, qv_s, t_s, t_snow, t_snow_mult,   &
        &  w_snow, wtot_snow, rho_snow, rho_snow_mult, h_snow, freshsnow, w_i, t_so,    &
        &  w_so, w_so_ice, runoff_s, runoff_g, tch, tfv, t_2m, qv_2m, td_2m, rh_2m,     &
        &  u_10m, v_10m, shfl_s, lhfl_s)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ie, i_count, nsubs0, nsubs1, ke_soil, ke_snow 
  INTEGER, DIMENSION(i_count), INTENT(IN) ::                     &
                  idx_lst_lndp         ! index list for land points                    (  -  )
  REAL(wp), DIMENSION(ie), INTENT(IN) ::                    &
                  t,                 & ! temperature                                   (  K  )
                  qv,                & ! humidity                                      (kg/kg)
                  p0,                & ! pressure                                      ( Pa  )
                  ps                   ! surface pressure                              ( Pa  )
  REAL(wp), DIMENSION(ie,nsubs1), INTENT(IN) ::             &
                  subsfrac,          & ! fractions of tiles                            (  -  )
                  t_g_now,           & ! surface temperature                           (  K  )
                  t_g_new,           & ! surface temperature                           (  K  )
                  qv_s,              & ! specific humidity at the surface              (kg/kg)
                  t_s,               & ! temperature of the ground surface             (  K  )
                  t_snow,            & ! temperature of the snow                       (  K  )
                  w_snow,            & ! water content of snow                         (m H2O)
                  rho_snow,          & ! snow density                                  (kg/m**3)
                  h_snow,            & ! snow height
                  w_i,               & ! water content of interception water           (m H2O)
                  freshsnow,         & ! indicator for age of snow in top of snow layer(  -  )
                  runoff_s,          & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g,          & ! soil water runoff; sum over forecast          (kg/m2)
                  tch,               & ! turbulent transfer coefficient for heat       ( -- )
                  tfv,               & ! laminar reduction factor for evaporation      ( -- )
                  t_2m,              & ! temperature in 2m                             (  K  )
                  qv_2m,             & ! humidity in 2m                                (kg/kg)
                  td_2m,             & ! dew point in 2m                               (  K  )
                  rh_2m,             & ! relative humidity in 2m                       (  K  )
                  u_10m,             & ! zonal wind in 10m                             ( m/s )
                  v_10m,             & ! meridional wind in 10m                        ( m/s )
                  shfl_s,            & ! sensible heat flux                            (W/m2 )
                  lhfl_s               ! latent heat flux                              (W/m2 )
  REAL(wp), DIMENSION(ie,0:ke_snow,nsubs1), INTENT(IN) ::   &
                  t_snow_mult          ! temperature of the snow-surface               (  K  )
  REAL(wp), DIMENSION(ie,1:ke_snow,nsubs1), INTENT(IN) ::   &
                  rho_snow_mult,     & ! snow density                                  (kg/m**3)
                  wtot_snow            ! total water content of snow layers            (m H2O)  
  REAL(wp), DIMENSION(ie,0:ke_soil+1,nsubs1), INTENT(IN) :: &
                  t_so                 ! soil temperature (main level)                 (  K  )
  REAL(wp), DIMENSION(ie,ke_soil+1,nsubs1), INTENT(IN) ::   &
                  w_so,              & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice             ! ice content       
!mean over tiles:
  REAL(wp), DIMENSION(ie), INTENT(OUT) ::                   &
                  t_g_s,             & ! surface temperature                           (  K  )
                  qv_s_s,            & ! specific humidity at the surface              (kg/kg)
                  t_s_s,             & ! temperature of the ground surface             (  K  )
                  t_snow_s,          & ! temperature of the snow                       (  K  )
                  w_snow_s,          & ! water content of snow                         (m H2O)
                  rho_snow_s,        & ! snow density                                  (kg/m**3)
                  h_snow_s,          & ! snow height
                  w_i_s,             & ! water content of interception water           (m H2O)
                  freshsnow_s,       & ! indicator for age of snow in top of snow layer(  -  )
                  runoff_s_s,        & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g_s,        & ! soil water runoff; sum over forecast          (kg/m2)
                  tch_s,             & ! turbulent transfer coefficient for heat       ( -- )
                  tfv_s,             & ! laminar reduction factor for evaporation      ( -- )
                  t_2m_s,            & ! temperature in 2m                             (  K  )
                  qv_2m_s,           & ! humidity in 2m                                (kg/kg)
                  td_2m_s,           & ! dew point in 2m                               (  K  )
                  rh_2m_s,           & ! relative humidity in 2m                       (  K  )
                  u_10m_s,           & ! zonal wind in 10m                             ( m/s )
                  v_10m_s,           & ! meridional wind in 10m                        ( m/s )
                  shfl_s_s,          & ! sensible heat flux                            (W/m2 )
                  lhfl_s_s             ! latent heat flux                              (W/m2 )
  REAL(wp), DIMENSION(ie,0:ke_snow), INTENT(OUT) ::         &
                  t_snow_mult_s        ! temperature of the snow-surface               (  K  )
  REAL(wp), DIMENSION(ie,1:ke_snow), INTENT(OUT) ::         &
                  rho_snow_mult_s,   & ! snow density                                  (kg/m**3)
                  wtot_snow_s          ! total water content of snow layers            (m H2O)  
  REAL(wp), DIMENSION(ie,0:ke_soil+1), INTENT(OUT) ::       &
                  t_so_s               ! soil temperature (main level)                 (  K  )
  REAL(wp), DIMENSION(ie,ke_soil+1), INTENT(OUT) ::         &
                  w_so_s,            & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_s           ! ice content       

  ! Local scalars:

  INTEGER  :: ic,jc,jk,ns,kso,ksn      !loop indices
  INTEGER  :: counter
  REAL(wp), PARAMETER :: small = 1.0E-07_wp

  ! Local arrays:

  REAL(wp) :: tmp(i_count)
  LOGICAL  :: llandmask(i_count,nsubs1)
  REAL(wp) :: t_g_slp(i_count), qv_s_slp(i_count), t_s_slp(i_count), t_snow_slp(i_count),        &
    &         t_snow_mult_slp(i_count, nlev_snow), w_snow_slp(i_count),                          &
    &         wtot_snow_slp(i_count, nlev_snow),                                                 &
    &         rho_snow_slp(i_count), rho_snow_mult_slp(i_count, nlev_snow), h_snow_slp(i_count), &
    &         freshsnow_slp(i_count), w_i_slp(i_count), t_so_slp(i_count, nlev_soil+2),          &
    &         w_so_slp(i_count, nlev_soil+1), w_so_ice_slp(i_count, nlev_soil+1),                &
    &         runoff_s_slp(i_count),                                                             &
    &         runoff_g_slp(i_count), tch_slp(i_count), tfv_slp(i_count), t_2m_slp(i_count),      &
    &         qv_2m_slp(i_count), td_2m_slp(i_count), rh_2m_slp(i_count), u_10m_slp(i_count),    &
    &         v_10m_slp(i_count), shfl_s_slp(i_count), lhfl_s_slp(i_count)
  REAL(wp) :: t_lp(i_count), qv_lp(i_count), p0_lp(i_count), ps_lp(i_count),                     &
    &         t_g_now_tlp(i_count,nsubs1),                                                       &
    &         t_g_new_tlp(i_count,nsubs1), qv_s_tlp(i_count,nsubs1), t_s_tlp(i_count,nsubs1),    &
    &         t_snow_tlp(i_count,nsubs1), t_snow_mult_tlp(i_count,nlev_snow,nsubs1),             &
    &         w_snow_tlp(i_count,nsubs1), wtot_snow_tlp(i_count,nlev_snow,nsubs1),               &
    &         rho_snow_tlp(i_count,nsubs1), rho_snow_mult_tlp(i_count,nlev_snow,nsubs1),         &
    &         h_snow_tlp(i_count,nsubs1), freshsnow_tlp(i_count,nsubs1), w_i_tlp(i_count,nsubs1),&
    &         t_so_tlp(i_count,nlev_soil+2,nsubs1), w_so_tlp(i_count,nlev_soil+1,nsubs1),        &
    &         w_so_ice_tlp(i_count,nlev_soil+1,nsubs1), runoff_s_tlp(i_count,nsubs1),            &
    &         runoff_g_tlp(i_count,nsubs1), tch_tlp(i_count,nsubs1), tfv_tlp(i_count,nsubs1),    &
    &         t_2m_tlp(i_count,nsubs1), qv_2m_tlp(i_count,nsubs1), td_2m_tlp(i_count,nsubs1),    &
    &         rh_2m_tlp(i_count,nsubs1), u_10m_tlp(i_count,nsubs1), v_10m_tlp(i_count,nsubs1),   &
    &         shfl_s_tlp(i_count,nsubs1), lhfl_s_tlp(i_count,nsubs1)
!-------------------------------------------------------------------------

  tmp(:)    = 0._wp
  tch_s(:)  = 0._wp
  qv_s_s(:) = 0._wp

!Set llandmask           !EM: it is a constant field, should be set once, not at each time step!
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
      llandmask(ic,ns) = .TRUE.
    END DO
  END DO
  DO ic = 1, i_count
    llandmask(ic,20) = .FALSE.    !water
  END DO

  DO ic = 1, i_count       
    jc = idx_lst_lndp(ic)
    t_lp(ic) = t(jc)
    qv_lp(ic) = qv(jc)
    p0_lp(ic) = p0(jc)
    ps_lp(ic) = ps(jc)
  END DO

  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count       
      jc = idx_lst_lndp(ic)
      t_g_now_tlp  (ic,ns) = t_g_now  (jc,ns)
      t_g_new_tlp  (ic,ns) = t_g_new  (jc,ns)
      qv_s_tlp     (ic,ns) = qv_s     (jc,ns)
      t_s_tlp      (ic,ns) = t_s      (jc,ns)
      t_snow_tlp   (ic,ns) = t_snow   (jc,ns)
      w_snow_tlp   (ic,ns) = w_snow   (jc,ns)
      rho_snow_tlp (ic,ns) = rho_snow (jc,ns)
      h_snow_tlp   (ic,ns) = h_snow   (jc,ns)
      freshsnow_tlp(ic,ns) = freshsnow(jc,ns)
      w_i_tlp      (ic,ns) = w_i      (jc,ns)
      runoff_s_tlp (ic,ns) = runoff_s (jc,ns)
      runoff_g_tlp (ic,ns) = runoff_g (jc,ns)
      tch_tlp      (ic,ns) = tch      (jc,ns)
      tfv_tlp      (ic,ns) = tfv      (jc,ns)
      t_2m_tlp     (ic,ns) = t_2m     (jc,ns)
      qv_2m_tlp    (ic,ns) = qv_2m    (jc,ns)
      td_2m_tlp    (ic,ns) = td_2m    (jc,ns)
      rh_2m_tlp    (ic,ns) = rh_2m    (jc,ns)
      u_10m_tlp    (ic,ns) = u_10m    (jc,ns)
      v_10m_tlp    (ic,ns) = v_10m    (jc,ns)
      shfl_s_tlp   (ic,ns) = shfl_s   (jc,ns)
      lhfl_s_tlp   (ic,ns) = lhfl_s   (jc,ns)
      DO jk = 0, nlev_soil+1                         !EM order of loops? 
        t_so_tlp    (ic,jk,ns) = t_so    (jc,jk,ns)
      END DO 
      DO jk = 1, nlev_soil+1
        w_so_tlp    (ic,jk,ns) = w_so    (jc,jk,ns)
        w_so_ice_tlp(ic,jk,ns) = w_so_ice(jc,jk,ns)
      END DO 
      DO jk = 0, nlev_snow+1
        t_snow_mult_tlp  (ic,jk,ns) = t_snow_mult  (jc,jk,ns)
      END DO 
      DO jk = 1, nlev_snow+1
        wtot_snow_tlp    (ic,jk,ns) = wtot_snow    (jc,jk,ns)
        rho_snow_mult_tlp(ic,jk,ns) = rho_snow_mult(jc,jk,ns)
      END DO 
    END DO 
  END DO 

!mean (effective) transfer coefficient for scalars
  CALL subsmean_power4(i_count, nsubs0, nsubs1, t_g_now_tlp, t_g_slp, subsfrac)
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
!intermediate storage for mean sensible heat flux
      tmp(ic) = tmp(ic) + tch_tlp(ic,ns)*(t_g_now_tlp(ic,ns) - &
                t_lp(ic) * (ps_lp(ic)/p0_lp(ic))**rdocp)*subsfrac(ic,ns)
    END DO
  END DO
  DO ic = 1, i_count
    tch_slp(ic) = tmp(ic)/SIGN(MAX(ABS(t_g_slp(ic) - t_lp(ic)                 &
      &           *(ps_lp(ic)/p0_lp(ic))**rdocp),small), t_g_slp(ic)-t_lp(ic) &
      &           *(ps_lp(ic)/p0_lp(ic))**rdocp)
    tch_slp(ic) = SIGN(MAX(ABS(tch_slp(ic)),small),tch_slp(ic))
  END DO

!mean (effective) surface humidity
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
!intermediate storage for mean latent heat flux
      tmp(ic) = tmp(ic) + tch_tlp(ic,ns)*(qv_s_tlp(ic,ns)-qv_lp(ic))*subsfrac(ic,ns)
    END DO
  END DO
  DO ic = 1, i_count
    qv_s_slp(ic) = MAX(tmp(ic)/tch_slp(ic) + qv_lp(ic),small)
  END DO

  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_s_tlp, t_s_slp, subsfrac)
  CALL subsmean_power4    (i_count, nsubs0, nsubs1, t_g_new_tlp, t_g_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_snow_tlp, t_snow_slp, subsfrac, &
    &                      w_snow_tlp>0.0_wp)
  DO ksn = 0, ke_snow
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_snow_mult_tlp(:,:,ksn),   &
      &                      t_snow_mult_slp(:,ksn), subsfrac, llandmask)
  END DO
  DO ic = 1, i_count
    counter = COUNT(w_snow_tlp(ic,:) .GT. small)
    IF(counter .EQ. 0) t_snow_slp(ic) = t_s_slp(ic)
  END DO
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_snow_tlp, w_snow_slp, subsfrac, llandmask)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, rho_snow_tlp,  rho_snow_slp, subsfrac,  &
    &                      w_snow_tlp>0.0_wp)
  DO ksn = 1, ke_snow
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, rho_snow_mult_tlp(:,:,ksn), &
                             rho_snow_mult_slp(:,ksn), subsfrac, w_snow_tlp>0.0_wp)
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, wtot_snow_tlp(:,:,ksn),     &
      &                      wtot_snow_slp(:,ksn), subsfrac, llandmask)
  END DO
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, h_snow_tlp, h_snow_slp, subsfrac, llandmask)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, freshsnow_tlp, freshsnow_slp, subsfrac,  &
    &                      w_snow_tlp>0.0_wp)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_i_tlp, w_i_slp, subsfrac, llandmask)
  DO kso = 0, ke_soil+1
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_so_tlp(:,:,kso), t_so_slp(:,kso),    &
      &                      subsfrac, llandmask)
  END DO
  DO kso = 1, ke_soil+1
!    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_tlp(:,:,kso), w_so_slp(:,kso), subsfrac, xlsmmask)
!    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_ice_tlp(:,:,kso), w_so_ice_slp(:,kso), subsfrac, xlsmmask)
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_tlp(:,:,kso), &
      w_so_slp(:,kso), subsfrac, llandmask)           !temporary; should not be averaged over rock
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_ice_tlp(:,:,kso), &
      w_so_ice_slp(:,kso), subsfrac, llandmask)   !temporary; should not be averaged over rock   
  END DO
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, runoff_s_tlp, runoff_s_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, runoff_g_tlp, runoff_g_slp, subsfrac)

  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, tfv_tlp, tfv_slp, subsfrac)

! These variables are averaged inside the radiation module:
! alb_rad, sobs, thbs, pabs

!diagnostics:
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,   t_2m_tlp,   t_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  qv_2m_tlp,  qv_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  td_2m_tlp,  td_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  rh_2m_tlp,  rh_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  u_10m_tlp,  u_10m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  v_10m_tlp,  v_10m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, shfl_s_tlp, shfl_s_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, lhfl_s_tlp, lhfl_s_slp, subsfrac)

  DO ic = 1, i_count                                  
    jc = idx_lst_lndp(ic)
    t_g_s      (jc) = t_g_slp      (ic)               
    qv_s_s     (jc) = qv_s_slp     (ic)
    t_s_s      (jc) = t_s_slp      (ic)
    t_snow_s   (jc) = t_snow_slp   (ic)
    w_snow_s   (jc) = w_snow_slp   (ic)
    rho_snow_s (jc) = rho_snow_slp (ic)
    h_snow_s   (jc) = h_snow_slp   (ic)
    freshsnow_s(jc) = freshsnow_slp(ic)
    w_i_s      (jc) = w_i_slp      (ic)
    runoff_s_s (jc) = runoff_s_slp (ic)
    runoff_g_s (jc) = runoff_g_slp (ic)
    tch_s      (jc) = tch_slp      (ic)
    tfv_s      (jc) = tfv_slp      (ic)
    t_2m_s     (jc) = t_2m_slp     (ic)
    qv_2m_s    (jc) = qv_2m_slp    (ic)
    td_2m_s    (jc) = td_2m_slp    (ic)
    rh_2m_s    (jc) = rh_2m_slp    (ic)
    u_10m_s    (jc) = u_10m_slp    (ic)
    v_10m_s    (jc) = v_10m_slp    (ic)
    shfl_s_s   (jc) = shfl_s_slp   (ic)
    lhfl_s_s   (jc) = lhfl_s_slp   (ic)
    DO jk = 0, nlev_soil+1                         !EM order of loops?
      t_so_s    (jc,jk) = t_so_slp    (ic,jk)
    END DO
    DO jk = 1, nlev_soil+1
      w_so_s    (jc,jk) = w_so_slp    (ic,jk)
      w_so_ice_s(jc,jk) = w_so_ice_slp(ic,jk)
    END DO
    DO jk = 0, nlev_snow+1
      t_snow_mult_s  (jc,jk) = t_snow_mult_slp  (ic,jk)
    END DO
    DO jk = 1, nlev_snow+1
      wtot_snow_s    (jc,jk) = wtot_snow_slp    (ic,jk)
      rho_snow_mult_s(jc,jk) = rho_snow_mult_slp(ic,jk)
    END DO
  ENDDO

  END SUBROUTINE subsmean

!-------------------------------------------------------------------------

  SUBROUTINE subsmean_arithmetic(i_count, nsubs0, nsubs1, field, field_s, subsfrac, mask)

  INTEGER , INTENT(IN) :: i_count, nsubs0, nsubs1
  REAL(wp), DIMENSION(i_count, nsubs1), INTENT(IN) :: field, subsfrac
  LOGICAL , DIMENSION(i_count, nsubs1), INTENT(IN), OPTIONAL :: mask
  REAL(wp), DIMENSION(i_count), INTENT(OUT) :: field_s 

  ! Local scalars:
  INTEGER :: ic, ns

  ! Local arrays:
  REAL(wp):: tmp(i_count)       ! total fraction of tiles with mask=true
  LOGICAL :: mask_tmp(i_count)

  IF(PRESENT(mask)) THEN
    DO ic = 1, i_count
      mask_tmp(ic) = .FALSE.
    END DO 
    DO ns = nsubs0, nsubs1
      DO ic = 1, i_count
        IF(mask(ic,ns)) mask_tmp(ic) = .TRUE.
      END DO
    END DO 
    DO ic = 1, i_count
      IF(mask_tmp(ic)) THEN
        field_s(ic) = 0._wp
        tmp(ic)     = 0._wp
      END IF
    END DO
    DO ns = nsubs0, nsubs1
      DO ic = 1, i_count
        IF(mask(ic,ns)) THEN
          field_s(ic) = field_s(ic) + field(ic,ns)*subsfrac(ic,ns)
          tmp(ic)     = tmp(ic) + subsfrac(ic,ns)
        END IF
      END DO
    END DO
    DO ic = 1, i_count
      IF(mask_tmp(ic)) THEN
        field_s(ic) = field_s(ic)/MAX(tmp(ic),1.E-07_wp)
      END IF
    END DO
  ELSE  ! mask is not present
    DO ic = 1, i_count
      field_s(ic) = 0._wp
    END DO 
    DO ns = nsubs0, nsubs1
      DO ic = 1, i_count
        field_s(ic) = field_s(ic) + field(ic,ns)*subsfrac(ic,ns)
      END DO
    END DO
  END IF

  END SUBROUTINE subsmean_arithmetic

!-------------------------------------------------------------------------

  SUBROUTINE subsmean_power4(i_count, nsubs0, nsubs1, field, field_s, subsfrac)

  INTEGER , INTENT(IN) :: i_count, nsubs0, nsubs1
  REAL(wp), DIMENSION(i_count, nsubs1), INTENT(IN) :: field, subsfrac
  REAL(wp), DIMENSION(i_count), INTENT(OUT) :: field_s 

  ! Local scalars:
  INTEGER :: ic, ns

!Without any mask because this kind of averaging (power4) is being done only for t_g

  DO ic = 1, i_count
    field_s(ic) = 0._wp
  END DO
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
      field_s(ic) = field_s(ic) + field(ic,ns)**4*subsfrac(ic,ns)
    END DO
  END DO
  DO ic = 1, i_count
    field_s(ic) = SQRT(SQRT(field_s(ic)))
  END DO

  END SUBROUTINE subsmean_power4

END MODULE mo_nwp_sfc_interface

