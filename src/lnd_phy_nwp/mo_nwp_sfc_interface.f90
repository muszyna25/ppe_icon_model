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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_sfc_interface

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message !, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, zml_soil
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag 
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv, msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water, &
    &                               lseaice, llake, lmulti_snow, ntiles_lnd, lsnowtile
  USE mo_satad,               ONLY: sat_pres_water, spec_humi  
  USE mo_soil_ml,             ONLY: terra_multlay
  USE mo_nwp_sfc_utils,       ONLY: diag_snowfrac_tg, update_index_lists, update_index_lists_sea
  USE mo_seaice_nwp,          ONLY: seaice_timestep_nwp
  USE mo_phyparam_soil              ! soil and vegetation parameters for TILES
!  USE mo_aggregate_surface,   ONLY: subsmean,subs_disaggregate_radflux,subsmean_albedo
!  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, igbm, iwtr, iice, ilnd
  USE mo_physical_constants,  ONLY: tmelt

  
  IMPLICIT NONE 

  PUBLIC  ::  nwp_surface

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
                        & p_prog_wtr_now, p_prog_wtr_new, & !>inout
                        & lnd_diag                        ) !>inout

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch       !< grid/patch info
    TYPE(t_external_data),       INTENT(inout):: ext_data      !< external data
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf    !< call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag        !< diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now  !< prog vars for sfc
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new  !< prog vars for sfc
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_now !< prog vars for wtr
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_new !< prog vars for wtr
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: tcall_sfc_jg  !< time interval for 
                                                               !< surface

    ! Local array bounds:
    !
    INTEGER :: nblks_c                 !> number of blocks for cells
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: nlev                    !< number of full levels
    INTEGER :: isubs, isubs_snow

    ! Local scalars:
    !
    INTEGER :: jc,jb,jg,jk      !loop indices


    REAL(wp) :: ps_t        (nproma, p_patch%nblks_c)
    REAL(wp) :: prr_con_t   (nproma, p_patch%nblks_c)
    REAL(wp) :: prs_con_t   (nproma, p_patch%nblks_c)
    REAL(wp) :: prr_gsp_t   (nproma, p_patch%nblks_c)
    REAL(wp) :: prs_gsp_t   (nproma, p_patch%nblks_c)

    REAL(wp) :: u_t (nproma,  p_patch%nblks_c)
    REAL(wp) :: v_t (nproma,  p_patch%nblks_c)
    REAL(wp) :: t_t (nproma,  p_patch%nblks_c)
    REAL(wp) :: qv_t(nproma,  p_patch%nblks_c)
    REAL(wp) :: p0_t(nproma,  p_patch%nblks_c)

    REAL(wp) :: sso_sigma_t(nproma,  p_patch%nblks_c)
    INTEGER  :: lc_class_t (nproma,  p_patch%nblks_c, ntiles_total)

    REAL(wp) :: t_snow_now_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_snow_new_t (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: t_s_now_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_new_t  (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: t_g_t      (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_g_now_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_g_new_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: qv_s_t     (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: w_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_snow_new_t(nproma, p_patch%nblks_c, ntiles_total)
  
    REAL(wp) :: rho_snow_now_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_new_t (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: h_snow_t (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: w_i_now_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_i_new_t (nproma, p_patch%nblks_c, ntiles_total)

!   REAL(wp) :: t_2m_t     (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: u_10m_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: v_10m_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: freshsnow_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: snowfrac_t (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: tch_t      (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: tcm_t      (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: tfv_t      (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: sobs_t     (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: thbs_t     (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: pabs_t     (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: runoff_s_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: runoff_g_t (nproma, p_patch%nblks_c, ntiles_total)

    INTEGER  :: soiltyp_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: plcov_t   (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rootdp_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: sai_t     (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: tai_t     (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: eai_t     (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rsmin2d_t (nproma, p_patch%nblks_c, ntiles_total)

    ! local dummy variable for precipitation rate of graupel, grid-scale
    REAL(wp) :: dummy_prg_gsp(nproma)

    REAL(wp) :: t_snow_mult_now_t(nproma, nlev_snow+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_snow_mult_new_t(nproma, nlev_snow+1, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: rho_snow_mult_now_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_mult_new_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: wliq_snow_now_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wliq_snow_new_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: wtot_snow_now_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wtot_snow_new_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: dzh_snow_now_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: dzh_snow_new_t(nproma, nlev_snow, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: t_so_now_t(nproma, nlev_soil+2, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_new_t(nproma, nlev_soil+2, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: w_so_now_t(nproma, nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_new_t(nproma, nlev_soil+1, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: w_so_ice_now_t(nproma, nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_new_t(nproma, nlev_soil+1, p_patch%nblks_c, ntiles_total)

    INTEGER  :: i_count, i_count_snow, ic, icount_init, is1, is2, init_list(2*nproma), it1(nproma), it2(nproma)
    REAL(wp) :: tmp1, tmp2, tmp3
    REAL(wp) :: frac_sv(nproma), frac_snow_sv(nproma), fact1(nproma), fact2(nproma)
    REAL(wp) :: tracer_rate(nproma, p_patch%nblks_c, 4, ntiles_total)
    REAL(wp), PARAMETER :: small = 1.E-06_wp

    REAL(wp) :: t_g_s(nproma), qv_s_s(nproma)
    REAL(wp) :: shfl_s_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: lhfl_s_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: shfl_snow_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: lhfl_snow_t (nproma, p_patch%nblks_c, ntiles_total)
!--------------------------------------------------------------


    ! initialize dummy variable
    dummy_prg_gsp(1:nproma) = 0._wp

    ! local variables related to the blocking

    nblks_c   = p_patch%nblks_int_c

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    IF (msg_level >= 15) THEN
      CALL message('mo_nwp_sfc_interface: ', 'call land-surface scheme')
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,isubs,i_count,ic,t_g_s,qv_s_s,isubs_snow,i_count_snow,&
!$OMP   tmp1,tmp2,tmp3,fact1,fact2,frac_sv,frac_snow_sv,icount_init,init_list,it1,it2,is1,is2) ICON_OMP_GUIDED_SCHEDULE

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
 
      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 .and. &
          & atm_phy_nwp_config(jg)%inwp_turb    /= 3 ) THEN

       IF (ext_data%atm%lp_count(jb) == 0) CYCLE ! skip loop if there is no land point

       ! Copy precipitation fields for subsequent downscaling
       DO isubs = 1,ntiles_total
         i_count = ext_data%atm%gp_count_t(jb,isubs) 
         IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
           tracer_rate(jc,jb,1,isubs) = prm_diag%tracer_rate(jc,jb,1)
           tracer_rate(jc,jb,2,isubs) = prm_diag%tracer_rate(jc,jb,2)
           tracer_rate(jc,jb,3,isubs) = prm_diag%tracer_rate(jc,jb,3)
           tracer_rate(jc,jb,4,isubs) = prm_diag%tracer_rate(jc,jb,4)
         END DO
       END DO

!---------- Preparations for TERRA in the case if snow tiles are considered
       IF(lsnowtile) THEN      ! snow is considered as separate tiles

         DO isubs = 1, ntiles_lnd

           isubs_snow = isubs + ntiles_lnd
           i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow) 

!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
  
             ! Snow and rain fall onto snow-covered tile surface only, 
             ! if 1) the corresponding snow tile already exists and 
             ! 2) the temperature of snow-free tile is below freezing point.
             ! If the temperature of snow-free tile is above freezing point,
             ! precipitation over it will be processed by this tile itself (no snow is created).
             ! If there is no snow tile so far at all, precipitation falls on the snow-free tile,
             ! and the snow tile will be created after TERRA.
             IF(lnd_prog_now%t_snow_t(jc,jb,isubs) > tmelt) THEN
               tracer_rate(jc,jb,1,isubs) = 0._wp
               tracer_rate(jc,jb,2,isubs) = 0._wp
               tracer_rate(jc,jb,3,isubs) = 0._wp
               tracer_rate(jc,jb,4,isubs) = 0._wp
             END IF
           END DO
         END DO
       END IF

!---------- Copy input fields for each tile

!----------------------------------
       DO isubs = 1,ntiles_total
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs) 

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)

          ps_t(ic,jb)           =  p_diag%pres_sfc(jc,jb)    
          prr_con_t(ic,jb)      =  tracer_rate(jc,jb,3,isubs)
          prs_con_t(ic,jb)      =  tracer_rate(jc,jb,4,isubs)
          prr_gsp_t(ic,jb)      =  tracer_rate(jc,jb,1,isubs)
          prs_gsp_t(ic,jb)      =  tracer_rate(jc,jb,2,isubs)

          u_t(ic,jb)      =  p_diag%u         (jc,nlev,jb)     
          v_t(ic,jb)      =  p_diag%v         (jc,nlev,jb)     
          t_t(ic,jb)      =  p_diag%temp      (jc,nlev,jb)     
          qv_t(ic,jb)     =  p_prog_rcf%tracer(jc,nlev,jb,iqv) 
          p0_t(ic,jb)     =  p_diag%pres      (jc,nlev,jb) 
    
          sso_sigma_t(ic,jb)       = ext_data%atm%sso_stdh(jc,jb)
          lc_class_t(ic,jb,isubs)  = ext_data%atm%lc_class_t(jc,jb,isubs)

          t_snow_now_t(ic,jb,isubs)          =  lnd_prog_now%t_snow_t(jc,jb,isubs) 
          t_s_now_t(ic,jb,isubs)             =  lnd_prog_now%t_s_t(jc,jb,isubs)   
          t_g_now_t (ic,jb,isubs)            =  lnd_prog_now%t_g_t(jc,jb,isubs)
          t_g_t (ic,jb,isubs)                =  lnd_prog_now%t_g_t(jc,jb,isubs)
          qv_s_t(ic,jb,isubs)                =  lnd_diag%qv_s_t(jc,jb,isubs)  
          w_snow_now_t(ic,jb,isubs)          =  lnd_prog_now%w_snow_t(jc,jb,isubs)
          rho_snow_now_t(ic,jb,isubs)        =  lnd_prog_now%rho_snow_t(jc,jb,isubs)
          w_i_now_t(ic,jb,isubs)             =  lnd_prog_now%w_i_t(jc,jb,isubs)
          freshsnow_t(ic,jb,isubs)           =  lnd_diag%freshsnow_t(jc,jb,isubs)
          snowfrac_t(ic,jb,isubs)            =  lnd_diag%snowfrac_t(jc,jb,isubs)
          runoff_s_t(ic,jb,isubs)            =  lnd_diag%runoff_s_t(jc,jb,isubs) 
          runoff_g_t(ic,jb,isubs)            =  lnd_diag%runoff_g_t(jc,jb,isubs)
!         t_2m_t(ic,jb,isubs)                =  prm_diag%t_2m(jc,jb) 
          u_10m_t(ic,jb,isubs)               =  prm_diag%u_10m(jc,jb)
          v_10m_t(ic,jb,isubs)               =  prm_diag%v_10m(jc,jb)  
          tch_t(ic,jb,isubs)                 =  prm_diag%tch(jc,jb)
          tcm_t(ic,jb,isubs)                 =  prm_diag%tcm(jc,jb)
          tfv_t(ic,jb,isubs)                 =  prm_diag%tfv(jc,jb)
          sobs_t(ic,jb,isubs)                =  prm_diag%swflxsfc_t(jc,jb,isubs) 
          thbs_t(ic,jb,isubs)                =  prm_diag%lwflxsfc_t(jc,jb,isubs) 
          pabs_t(ic,jb,isubs)                =  prm_diag%swflxsfc_t(jc,jb,isubs) 

          soiltyp_t(ic,jb,isubs)             =  ext_data%atm%soiltyp_t(jc,jb,isubs)
          plcov_t(ic,jb,isubs)               =  ext_data%atm%plcov_t(jc,jb,isubs)
          rootdp_t(ic,jb,isubs)              =  ext_data%atm%rootdp_t(jc,jb,isubs)
          sai_t(ic,jb,isubs)                 =  ext_data%atm%sai_t(jc,jb,isubs)
          tai_t(ic,jb,isubs)                 =  ext_data%atm%tai_t(jc,jb,isubs)
          eai_t(ic,jb,isubs)                 =  ext_data%atm%eai_t(jc,jb,isubs)
          rsmin2d_t(ic,jb,isubs)             =  ext_data%atm%rsmin2d_t(jc,jb,isubs)

          t_so_now_t(ic,nlev_soil+2,jb,isubs) = lnd_prog_now%t_so_t(jc,nlev_soil+2,jb,isubs)

          IF(lmulti_snow) THEN
            t_snow_mult_now_t(ic,nlev_snow+1,jb,isubs) = &
              lnd_prog_now%t_snow_mult_t(jc,nlev_snow+1,jb,isubs)
            h_snow_t(ic,jb,isubs)  =  lnd_diag%h_snow_t(jc,jb,isubs)
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
            t_snow_mult_now_t  (ic,jk,jb,isubs) = lnd_prog_now%t_snow_mult_t  (jc,jk,jb,isubs) 
            rho_snow_mult_now_t(ic,jk,jb,isubs) = lnd_prog_now%rho_snow_mult_t(jc,jk,jb,isubs)
            wliq_snow_now_t    (ic,jk,jb,isubs) = lnd_prog_now%wliq_snow_t    (jc,jk,jb,isubs) 
            wtot_snow_now_t    (ic,jk,jb,isubs) = lnd_prog_now%wtot_snow_t    (jc,jk,jb,isubs)
            dzh_snow_now_t     (ic,jk,jb,isubs) = lnd_prog_now%dzh_snow_t     (jc,jk,jb,isubs) 
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
            t_so_now_t    (ic,jk,jb,isubs) = lnd_prog_now%t_so_t    (jc,jk,jb,isubs) 
            w_so_now_t    (ic,jk,jb,isubs) = lnd_prog_now%w_so_t    (jc,jk,jb,isubs) 
            w_so_ice_now_t(ic,jk,jb,isubs) = lnd_prog_now%w_so_ice_t(jc,jk,jb,isubs)
          ENDDO
        ENDDO
!
!---------- END Copy index list fields


        CALL terra_multlay(                                    &
        &  ie=nproma                                         , & ! array dimensions
        &  istartpar=1,       iendpar=i_count                , & ! optional start/end indicies
        &  nsubs0=jb,          nsubs1=isubs                  , & ! unused except for optional debug output
        &  ke_soil=nlev_soil, ke_snow=nlev_snow              , &
        &  czmls=zml_soil,    ldiag_tg=.FALSE.               , & ! processing soil level structure 
        &  inwp_turb=atm_phy_nwp_config(jg)%inwp_turb        , &
        &  dt=tcall_sfc_jg                                   , &
        &  soiltyp_subs = soiltyp_t(:,jb,isubs)              , & ! type of the soil (keys 0-9)         --    
        &  plcov        = plcov_t(:,jb,isubs)                , & ! fraction of plant cover             --
        &  rootdp       = rootdp_t(:,jb,isubs)               , & ! depth of the roots                ( m  )
        &  sai          = sai_t(:,jb,isubs)                  , & ! surface area index                  --
        &  tai          = tai_t(:,jb,isubs)                  , & ! surface area index                  --
        &  eai          = eai_t(:,jb,isubs)                  , & ! surface area index                  --
        &  rsmin2d      = rsmin2d_t(:,jb,isubs)              , & ! minimum stomata resistance        ( s/m )
!
        &  u  =  u_t(:,jb)                                   , & ! zonal wind speed
        &  v  =  v_t(:,jb)                                   , & ! meridional wind speed 
        &  t  =  t_t(:,jb)                                   , & ! temperature                       (  K  )
        &  qv =  qv_t(:,jb)                                  , & ! specific water vapor content      (kg/kg)
        &  p0 =  p0_t(:,jb)                                  , & ! base state pressure               ( Pa  ) 
        &  ps =  ps_t(:,jb)                                  , & ! surface pressure                  ( pa  )
!
        &  t_snow_now    = t_snow_now_t(:,jb,isubs)          , & ! temperature of the snow-surface   (  K  )
        &  t_snow_new    = t_snow_new_t(:,jb,isubs)          , & ! temperature of the snow-surface   (  K  )
!
        &  t_snow_mult_now = t_snow_mult_now_t(:,:,jb,isubs) , & ! temperature of the snow-surface   (  K  )
        &  t_snow_mult_new = t_snow_mult_new_t(:,:,jb,isubs) , & ! temperature of the snow-surface   (  K  )
!
        &  t_s_now           = t_s_now_t(:,jb,isubs)         , & ! temperature of the ground surface (  K  )
        &  t_s_new           = t_s_new_t(:,jb,isubs)         , & ! temperature of the ground surface (  K  )
!
        &  t_g           =  t_g_t (:,jb,isubs)               , & ! weighted surface temperature      (  K  )
!EM     &  t_g_now       =  t_g_now_t (:,jb,isubs)           , & ! weighted surface temperature      (  K  )
!EM     &  t_g_new       =  t_g_new_t (:,jb,isubs)           , & ! weighted surface temperature      (  K  )
        &  qv_s          =  qv_s_t(:,jb,isubs)               , & ! specific humidity at the surface  (kg/kg)
!
        &  w_snow_now    = w_snow_now_t(:,jb,isubs)          , & ! water content of snow             (m H2O) 
        &  w_snow_new    = w_snow_new_t(:,jb,isubs)          , & ! water content of snow             (m H2O) 
!
        &  rho_snow_now      = rho_snow_now_t(:,jb,isubs)    , & ! snow density                      (kg/m**3)
        &  rho_snow_new      = rho_snow_new_t(:,jb,isubs)    , & ! snow density                      (kg/m**3)
!
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,jb,isubs), & ! snow density                   (kg/m**3) 
        &  rho_snow_mult_new = rho_snow_mult_new_t(:,:,jb,isubs), & ! snow density                   (kg/m**3) 
!
        &  h_snow        =  h_snow_t(:,jb,isubs)             , & ! snow height
!
        &  w_i_now       =   w_i_now_t(:,jb,isubs)           , & ! water content of interception water(m H2O)
        &  w_i_new       =   w_i_new_t(:,jb,isubs)           , & ! water content of interception water(m H2O)
!
        &  t_so_now      = t_so_now_t(:,:,jb,isubs)          , & ! soil temperature (main level)      (  K  )
        &  t_so_new      = t_so_new_t(:,:,jb,isubs)          , & ! soil temperature (main level)      (  K  )
!
        &  w_so_now      = w_so_now_t(:,:,jb,isubs)          , & ! total water content (ice + liquid water) (m H20)
        &  w_so_new      = w_so_new_t(:,:,jb,isubs)          , & ! total water content (ice + liquid water) (m H20)
!
        &  w_so_ice_now  = w_so_ice_now_t(:,:,jb,isubs)      , & ! ice content   (m H20)
        &  w_so_ice_new  = w_so_ice_new_t(:,:,jb,isubs)      , & ! ice content   (m H20)
!
!       &  t_2m          =  t_2m_t(:,jb,isubs)               , & ! temperature in 2m                  (  K  )
        &  u_10m         =  u_10m_t(:,jb,isubs)              , & ! zonal wind in 10m                  ( m/s )
        &  v_10m         =  v_10m_t(:,jb,isubs)              , & ! meridional wind in 10m            ( m/s )
        &  freshsnow     =  freshsnow_t(:,jb,isubs)          , & ! indicator for age of snow in top of snow layer (  -  )
        &  zf_snow       =  snowfrac_t(:,jb,isubs)           , & ! snow-cover fraction                            (  -  )
!
        &  wliq_snow_now = wliq_snow_now_t(:,:,jb,isubs)     , & ! liquid water content in the snow       (m H2O)
        &  wliq_snow_new = wliq_snow_new_t(:,:,jb,isubs)     , & ! liquid water content in the snow       (m H2O)
!                                                            
        &  wtot_snow_now = wtot_snow_now_t(:,:,jb,isubs)     , & ! total (liquid + solid) water content of snow  (m H2O)
        &  wtot_snow_new = wtot_snow_new_t(:,:,jb,isubs)     , & ! total (liquid + solid) water content of snow  (m H2O)
!
        &  dzh_snow_now  = dzh_snow_now_t(:,:,jb,isubs)      , & ! layer thickness between half levels in snow   (  m  )
        &  dzh_snow_new  = dzh_snow_new_t(:,:,jb,isubs)      , & ! layer thickness between half levels in snow   (  m  )
!
        &  prr_con       = prr_con_t(:,jb)                   , & ! precipitation rate of rain, convective        (kg/m2*s)
        &  prs_con       = prs_con_t(:,jb)                   , & ! precipitation rate of snow, convective        (kg/m2*s)
        &  prr_gsp       = prr_gsp_t(:,jb)                   , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
        &  prs_gsp       = prs_gsp_t(:,jb)                   , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
        &  prg_gsp       = dummy_prg_gsp(:)                  , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
!
        &  tch           = tch_t(:,jb,isubs)                 , & ! turbulent transfer coefficient for heat       ( -- )
        &  tcm           = tcm_t(:,jb,isubs)                 , & ! turbulent transfer coefficient for momentum   ( -- )
        &  tfv           = tfv_t(:,jb,isubs)                 , & ! laminar reduction factor for evaporation      ( -- )
!
        &  sobs          = sobs_t(:,jb,isubs)                , & ! solar radiation at the ground                 (W/m2)
        &  thbs          = thbs_t(:,jb,isubs)                , & ! thermal radiation at the ground               (W/m2)
        &  pabs          = pabs_t(:,jb,isubs)                , & ! photosynthetic active radiation               (W/m2)
!
        &  runoff_s      = runoff_s_t(:,jb,isubs)            , & ! surface water runoff; sum over forecast       (kg/m2)
        &  runoff_g      = runoff_g_t(:,jb,isubs)            , & ! soil water runoff; sum over forecast          (kg/m2)
!
        &  zshfl_s       = shfl_s_t   (:,jb,isubs)           , & ! sensible heat flux soil/air interface         (W/m2) 
        &  zlhfl_s       = lhfl_s_t   (:,jb,isubs)           , & ! latent   heat flux soil/air interface         (W/m2) 
        &  zshfl_snow    = shfl_snow_t(:,jb,isubs)           , & ! sensible heat flux snow/air interface         (W/m2) 
        &  zlhfl_snow    = lhfl_snow_t(:,jb,isubs)             & ! latent   heat flux snow/air interface         (W/m2) 
        &                                                    )


        CALL diag_snowfrac_tg(                           &
          &  istart = 1, iend = i_count                , & ! start/end indices
          &  z0_lcc    = ext_data%atm%z0_lcc(:)        , & ! roughness length
          &  lc_class  = lc_class_t        (:,jb,isubs), & ! land-cover class
          &  t_snow    = t_snow_new_t      (:,jb,isubs), & ! snow temp
          &  t_soiltop = t_s_new_t         (:,jb,isubs), & ! soil top temp
          &  w_snow    = w_snow_new_t      (:,jb,isubs), & ! snow WE
          &  rho_snow  = rho_snow_new_t    (:,jb,isubs), & ! snow depth
          &  freshsnow = freshsnow_t       (:,jb,isubs), & ! fresh snow fraction
          &  sso_sigma = sso_sigma_t       (:,jb),       & ! sso stdev
          &  tai       = tai_t             (:,jb,isubs), & ! effective leaf area index
          &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
          &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp


!---------- Copy index list fields back to state fields

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          lnd_prog_new%t_snow_t  (jc,jb,isubs) = t_snow_new_t  (ic,jb,isubs)         
          lnd_prog_new%t_s_t     (jc,jb,isubs) = t_s_new_t     (ic,jb,isubs)              
          lnd_prog_new%t_g_t     (jc,jb,isubs) = t_g_t         (ic,jb,isubs)
          lnd_diag%qv_s_t        (jc,jb,isubs) = qv_s_t        (ic,jb,isubs)                
          lnd_prog_new%w_snow_t  (jc,jb,isubs) = w_snow_new_t  (ic,jb,isubs)          
          lnd_prog_new%rho_snow_t(jc,jb,isubs) = rho_snow_new_t(ic,jb,isubs)        
          lnd_diag%h_snow_t      (jc,jb,isubs) = h_snow_t      (ic,jb,isubs)              
          lnd_prog_new%w_i_t     (jc,jb,isubs) = w_i_new_t     (ic,jb,isubs)             
          lnd_diag%freshsnow_t   (jc,jb,isubs) = freshsnow_t   (ic,jb,isubs) 
          ! Remark: the two snow-cover fraction variables differ only if lsnowtile=true (see below)  
          lnd_diag%snowfrac_lc_t (jc,jb,isubs) = snowfrac_t    (ic,jb,isubs) 
          lnd_diag%snowfrac_t    (jc,jb,isubs) = snowfrac_t    (ic,jb,isubs) 
          lnd_diag%runoff_s_t    (jc,jb,isubs) = runoff_s_t    (ic,jb,isubs)  
          lnd_diag%runoff_g_t    (jc,jb,isubs) = runoff_g_t    (ic,jb,isubs)  

          lnd_prog_new%t_so_t(jc,nlev_soil+2,jb,isubs) = t_so_new_t(ic,nlev_soil+2,jb,isubs)

          IF(lmulti_snow) THEN
            lnd_prog_new%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = &
              t_snow_mult_new_t(ic,nlev_snow+1,jb,isubs)
          ENDIF
          t_g_new_t (ic,jb,isubs) = t_g_t (ic,jb,isubs)
        ENDDO

        IF (lsnowtile .AND. isubs > ntiles_lnd) THEN ! copy snowfrac_t to snow-free tile
!CDIR NODEP,VOVERTAKE,VOB                            ! (needed for index list computation)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            lnd_diag%snowfrac_lc_t(jc,jb,isubs-ntiles_lnd) = lnd_diag%snowfrac_lc_t(jc,jb,isubs)
          ENDDO
        ENDIF

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
            lnd_prog_new%t_snow_mult_t  (jc,jk,jb,isubs) = t_snow_mult_new_t  (ic,jk,jb,isubs)   
            lnd_prog_new%rho_snow_mult_t(jc,jk,jb,isubs) = rho_snow_mult_new_t(ic,jk,jb,isubs) 
            lnd_prog_new%wliq_snow_t    (jc,jk,jb,isubs) = wliq_snow_new_t    (ic,jk,jb,isubs)     
            lnd_prog_new%wtot_snow_t    (jc,jk,jb,isubs) = wtot_snow_new_t    (ic,jk,jb,isubs)     
            lnd_prog_new%dzh_snow_t     (jc,jk,jb,isubs) = dzh_snow_new_t     (ic,jk,jb,isubs)      
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
            lnd_prog_new%t_so_t    (jc,jk,jb,isubs) = t_so_new_t    (ic,jk,jb,isubs)          
            lnd_prog_new%w_so_t    (jc,jk,jb,isubs) = w_so_new_t    (ic,jk,jb,isubs)          
            lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)     
          ENDDO
        ENDDO

       END DO ! isubs - loop over tiles

       IF(lsnowtile) THEN      ! snow is considered as separate tiles
         DO isubs = 1, ntiles_lnd

           isubs_snow = isubs + ntiles_lnd

           ! save previous area fractions for subsequent redistribution computations
           frac_sv(:)      = ext_data%atm%frac_t(:,jb,isubs)
           frac_snow_sv(:) = ext_data%atm%frac_t(:,jb,isubs_snow)

           ! update index lists for snow tiles
           CALL update_index_lists (idx_lst_lp         = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                    lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                    idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                    gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                    idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                    gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                    lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                    partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                    partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                    snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                    snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                    snowfrac           = lnd_diag%snowfrac_lc_t(:,jb,isubs)             )
  
           i_count = ext_data%atm%gp_count_t(jb,isubs) 
           i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)

           ! Check for newly activated grid points that need to be initialized
           icount_init = 0
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs) == 2) THEN
               icount_init = icount_init + 1
               init_list(icount_init) = jc
               it1(icount_init) = isubs      ! target of copy operation
               it2(icount_init) = isubs_snow ! source of copy operation
             ENDIF
           ENDDO
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 2) THEN
               icount_init = icount_init + 1
               init_list(icount_init) = jc
               it1(icount_init) = isubs_snow ! target of copy operation
               it2(icount_init) = isubs      ! source of copy operation
             ENDIF
           ENDDO

           DO ic = 1, icount_init
             jc = init_list(ic)
             is1 = it1(ic)
             is2 = it2(ic)
             lnd_prog_new%t_snow_t  (jc,jb,is1) = lnd_prog_new%t_snow_t  (jc,jb,is2)        
             lnd_prog_new%t_s_t     (jc,jb,is1) = lnd_prog_new%t_s_t     (jc,jb,is2)       
             lnd_prog_new%t_g_t     (jc,jb,is1) = lnd_prog_new%t_g_t     (jc,jb,is2) 
             lnd_diag%qv_s_t        (jc,jb,is1) = lnd_diag%qv_s_t        (jc,jb,is2)             
             lnd_prog_new%w_snow_t  (jc,jb,is1) = lnd_prog_new%w_snow_t  (jc,jb,is2)     
             lnd_prog_new%rho_snow_t(jc,jb,is1) = lnd_prog_new%rho_snow_t(jc,jb,is2)
             lnd_diag%h_snow_t      (jc,jb,is1) = lnd_diag%h_snow_t      (jc,jb,is2)
             lnd_prog_new%w_i_t     (jc,jb,is1) = lnd_prog_new%w_i_t     (jc,jb,is2)        
             lnd_diag%freshsnow_t   (jc,jb,is1) = lnd_diag%freshsnow_t   (jc,jb,is2)
             lnd_diag%snowfrac_lc_t (jc,jb,is1) = lnd_diag%snowfrac_lc_t (jc,jb,is2) 
             lnd_diag%snowfrac_t    (jc,jb,is1) = lnd_diag%snowfrac_t    (jc,jb,is2) 
             lnd_diag%runoff_s_t    (jc,jb,is1) = lnd_diag%runoff_s_t    (jc,jb,is2)
             lnd_diag%runoff_g_t    (jc,jb,is1) = lnd_diag%runoff_g_t    (jc,jb,is2)

             lnd_prog_new%t_so_t    (jc,:,jb,is1) = lnd_prog_new%t_so_t    (jc,:,jb,is2)          
             lnd_prog_new%w_so_t    (jc,:,jb,is1) = lnd_prog_new%w_so_t    (jc,:,jb,is2)        
             lnd_prog_new%w_so_ice_t(jc,:,jb,is1) = lnd_prog_new%w_so_ice_t(jc,:,jb,is2)

             IF (lmulti_snow) THEN
               lnd_prog_new%t_snow_mult_t  (jc,:,jb,is1) = lnd_prog_new%t_snow_mult_t  (jc,:,jb,is2)
               lnd_prog_new%rho_snow_mult_t(jc,:,jb,is1) = lnd_prog_new%rho_snow_mult_t(jc,:,jb,is2)
               lnd_prog_new%wliq_snow_t    (jc,:,jb,is1) = lnd_prog_new%wliq_snow_t    (jc,:,jb,is2)
               lnd_prog_new%wtot_snow_t    (jc,:,jb,is1) = lnd_prog_new%wtot_snow_t    (jc,:,jb,is2)
               lnd_prog_new%dzh_snow_t     (jc,:,jb,is1) = lnd_prog_new%dzh_snow_t     (jc,:,jb,is2)
             ENDIF
           ENDDO

!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 1 .AND. &
                 ext_data%atm%snowtile_flag_t(jc,jb,isubs)      == 1) THEN

               ! compute factors for redistribution of heat and moisture
               fact1(jc) = MIN(1._wp,frac_sv(jc)/     MAX(small,ext_data%atm%frac_t(jc,jb,isubs)     ))
               fact2(jc) = MIN(1._wp,frac_snow_sv(jc)/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow)))
             ENDIF

           END DO

           ! redistribution of heat and moisture between snow-covered and snow-free tiles 
           ! according to their new fractions, in order to keep heat and moisture balances
           DO jk = 1, nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
             DO ic = 1, i_count_snow
               jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

               IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 1 .AND. &
                   ext_data%atm%snowtile_flag_t(jc,jb,isubs)      == 1) THEN

                 tmp1 = lnd_prog_new%t_so_t(jc,jk,jb,isubs) 
                 tmp2 = lnd_prog_new%w_so_t(jc,jk,jb,isubs)
                 tmp3 = lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs)
  
                 lnd_prog_new%t_so_t    (jc,jk,jb,isubs) = lnd_prog_new%t_so_t    (jc,jk,jb,isubs)*fact1(jc) &
                   &                       + lnd_prog_new%t_so_t    (jc,jk,jb,isubs_snow)*(1._wp - fact1(jc))
                 lnd_prog_new%w_so_t    (jc,jk,jb,isubs) = lnd_prog_new%w_so_t    (jc,jk,jb,isubs)*fact1(jc) &
                   &                       + lnd_prog_new%w_so_t    (jc,jk,jb,isubs_snow)*(1._wp - fact1(jc))
                 lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs) = lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs)*fact1(jc) &
                   &                       + lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs_snow)*(1._wp - fact1(jc))
 
                 lnd_prog_new%t_so_t    (jc,jk,jb,isubs_snow) = tmp1*(1._wp - fact2(jc)) &
                   &              + lnd_prog_new%t_so_t    (jc,jk,jb,isubs_snow)*fact2(jc)
                 lnd_prog_new%w_so_t    (jc,jk,jb,isubs_snow) = tmp2*(1._wp - fact2(jc)) &
                   &              + lnd_prog_new%w_so_t    (jc,jk,jb,isubs_snow)*fact2(jc)
                 lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs_snow) = tmp3*(1._wp - fact2(jc)) &
                  &               + lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs_snow)*fact2(jc)

                 IF (jk == 1) THEN
                   lnd_prog_new%t_s_t(jc,jb,isubs)      = lnd_prog_new%t_so_t(jc,jk,jb,isubs)
                   lnd_prog_new%t_s_t(jc,jb,isubs_snow) = lnd_prog_new%t_so_t(jc,jk,jb,isubs_snow)
                 ENDIF
               ENDIF

             END DO
           END DO        ! soil layers
!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             ! snow depth per surface unit -> snow depth per fraction
             lnd_diag%w_snow_eff_t(jc,jb,isubs_snow) = &
               lnd_prog_new%w_snow_t(jc,jb,isubs_snow)/MAX(lnd_diag%snowfrac_lc_t(jc,jb,isubs_snow),small)

             ! reset field for actual snow-cover for grid points / land-cover classes for which there
             ! are seperate snow-free and snow-covered tiles 
             lnd_diag%snowfrac_t(jc,jb,isubs)      = 0._wp
             lnd_prog_new%w_snow_t(jc,jb,isubs)    = 0._wp
             lnd_prog_new%t_snow_t(jc,jb,isubs)    = lnd_prog_new%t_s_t(jc,jb,isubs)
             lnd_prog_new%t_g_t(jc,jb,isubs)       = lnd_prog_new%t_s_t(jc,jb,isubs)

             ! to prevent numerical stability problems, we require at least 1 cm of snow in order to
             ! have a snow-cover fraction of 1 on snow tiles (not critical for the single-layer
             ! snow scheme, but the multi-layer snow model becomes numerically unstable within a few
             ! time steps when associating traces of snow with a snow-cover fraction of 1)
             lnd_diag%snowfrac_t(jc,jb,isubs_snow) = MIN(1._wp,lnd_diag%h_snow_t(jc,jb,isubs_snow)/0.01_wp)

             ! Rediagnose t_g according to the modified snow-cover fraction
             lnd_prog_new%t_g_t(jc,jb,isubs_snow) =  &
               lnd_diag%snowfrac_t(jc,jb,isubs_snow) * lnd_prog_new%t_snow_t(jc,jb,isubs_snow) + &
               (1._wp-lnd_diag%snowfrac_t(jc,jb,isubs_snow))*lnd_prog_new%t_s_t(jc,jb,isubs_snow)

             IF (lmulti_snow) THEN
               lnd_prog_new%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = lnd_prog_new%t_s_t(jc,jb,isubs)
             ENDIF
           END DO

           IF (lmulti_snow) THEN
!CDIR UNROLL=nlsnow
             DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
               DO ic = 1, i_count_snow
                 jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
                 lnd_prog_new%t_snow_mult_t(jc,jk,jb,isubs) = lnd_prog_new%t_s_t(jc,jb,isubs)
                 lnd_prog_new%wliq_snow_t(jc,jk,jb,isubs) = 0._wp
                 lnd_prog_new%wtot_snow_t(jc,jk,jb,isubs) = 0._wp
                 lnd_prog_new%dzh_snow_t (jc,jk,jb,isubs) = 0._wp
               ENDDO
             ENDDO
           ENDIF

         END DO

       ENDIF  !snow tiles

       ! Final step: aggregate t_g and qv_s
       i_count = ext_data%atm%lp_count(jb)

       IF (ntiles_total == 1) THEN 
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp(ic,jb)
           lnd_prog_new%t_g(jc,jb)  = lnd_prog_new%t_g_t(jc,jb,1)
           lnd_diag%qv_s(jc,jb)     = lnd_diag%qv_s_t(jc,jb,1) 
         ENDDO
       ELSE ! aggregate fields over tiles
         t_g_s(:)  =  0._wp
         qv_s_s(:) =  0._wp
         DO isubs = 1,ntiles_total+ntiles_water
!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_lp(ic,jb)
             t_g_s(jc) = t_g_s(jc) + ext_data%atm%frac_t(jc,jb,isubs)* &
               lnd_prog_new%t_g_t(jc,jb,isubs)**4
             qv_s_s(jc) = qv_s_s(jc) + ext_data%atm%frac_t(jc,jb,isubs)* & 
               lnd_diag%qv_s_t(jc,jb,isubs)
           ENDDO
         ENDDO
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp(ic,jb)
           lnd_prog_new%t_g(jc,jb)  = SQRT(SQRT(t_g_s(jc)))
           lnd_diag%qv_s(jc,jb)     = qv_s_s(jc)
         ENDDO

       ENDIF    ! with or without tiles


    
      ELSE IF ( atm_phy_nwp_config(jg)%inwp_surface == 2 ) THEN 

          !-------------------------------------------------------------------------
          !> ECHAM version 
          !-------------------------------------------------------------------------
     

     
      ENDIF !inwp_sfc

    ENDDO    
!$OMP END DO
!$OMP END PARALLEL



!!!!!!!!!!!!! SEAICE deactivated so far (under construction)
!       !
!       ! Call seaice parameterization
!       !
!       IF ( (atm_phy_nwp_config(jg)%inwp_surface == 1) .AND. (lseaice) ) THEN
!         CALL nwp_seaice(p_patch, prm_diag, p_prog_wtr_now, p_prog_wtr_new, &
!    &                    ext_data, tcall_sfc_jg)
!       ENDIF

  END SUBROUTINE nwp_surface



  !>
  !! Interface for seaice parameterization
  !!
  !! Interface for seaice parameterization. Calls seaice time integration scheme 
  !! seaice_timestep_nwp and updates the dynamic seaice index lists.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-08-31)
  !!
  SUBROUTINE nwp_seaice (p_patch, prm_diag, p_prog_wtr_now, p_prog_wtr_new, &
    &                    ext_data, dtime)

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !< grid/patch info
    TYPE(t_nwp_phy_diag),        INTENT(in)   :: prm_diag       !< atm phys vars
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_now !< prog vars for wtr
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_new !< prog vars for wtr
    TYPE(t_external_data),       INTENT(inout):: ext_data       !< external data
    REAL(wp),                    INTENT(in)   :: dtime          !< time interval for 
                                                                !< surface

    ! Local arrays  (local copies)
    !
    REAL(wp) :: shfl_s   (nproma)   ! sensible heat flux at the surface               [W/m^2]
    REAL(wp) :: lhfl_s   (nproma)   ! latent heat flux at the surface                 [W/m^2]
    REAL(wp) :: lwflxsfc (nproma)   ! net long-wave radiation flux at the surface     [W/m^2] 
    REAL(wp) :: swflxsfc (nproma)   ! !< net solar radiation flux at the surface      [W/m^2]
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time  [K]
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level               [m]
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time [K]
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level              [m]
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time       [K]
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level                    [m]
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time      [K]
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level                   [m]

    ! Local array bounds:
    !
    INTEGER :: nblks_c                 !> number of blocks for cells
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:
    !
    INTEGER :: jc, jb, ic              !loop indices
    INTEGER :: i_count

    !-------------------------------------------------------------------------


    ! local variables related to the blocking

    nblks_c   = p_patch%nblks_int_c

    i_nchdom  = MAX(1,p_patch%n_childdom)


    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    IF (msg_level >= 15) THEN
      CALL message('mo_nwp_sfc_interface: ', 'call nwp_seaice scheme')
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,i_count,ic) ICON_OMP_GUIDED_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


      !
      ! Copy input fields
      !
      i_count = ext_data%atm%spi_count(jb) 

      IF (i_count == 0) CYCLE ! skip loop if the index list for the given block is empty

      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_spi(ic,jb)

!!!!! should have the PRIVATE attribute, if not allocated as 2D-arrays!!!!!
!!!!!!!!!!!
        shfl_s   (ic) = prm_diag%shfl_s(jc,jb)     ! sensible heat flux at sfc    [W/m^2]
        lhfl_s   (ic) = prm_diag%lhfl_s(jc,jb)     ! latent heat flux at sfc      [W/m^2]
        lwflxsfc (ic) = prm_diag%lwflxsfc(jc,jb)   ! net lw radiation flux at sfc [W/m^2]
        swflxsfc (ic) = prm_diag%swflxsfc(jc,jb)   ! net solar radiation flux at sfc [W/m^2]
        tice_now (ic) = p_prog_wtr_now%t_ice(jc,jb)
        hice_now (ic) = p_prog_wtr_now%h_ice(jc,jb)
        tsnow_now(ic) = p_prog_wtr_now%t_snow_si(jc,jb)
        hsnow_now(ic) = p_prog_wtr_now%h_snow_si(jc,jb)
          
      ENDDO  ! ic


      ! call seaice time integration scheme
      !
      CALL seaice_timestep_nwp (                               &
                            &   dtime   = dtime,               &
                            &   nsigb   = i_count,             &
                            &   qsen    = shfl_s(:),           & !in 
                            &   qlat    = lhfl_s(:),           & !in
                            &   qlwrnet = lwflxsfc(:),         & !in
                            &   qsolnet = swflxsfc(:),         & !in
                            &   tice_p  = tice_now(:),         & !in
                            &   hice_p  = hice_now(:),         & !in
                            &   tsnow_p = tsnow_now(:),        & !in
                            &   hsnow_p = hsnow_now(:),        & !in
                            &   tice_n  = tice_new(:),         & !out
                            &   hice_n  = hice_new(:),         & !out
                            &   tsnow_n = tsnow_new(:),        & !out
                            &   hsnow_n = hsnow_new(:)         ) !out
! optional arguments dticedt, dhicedt, dtsnowdt, dhsnowdt (tendencies) neglected, so far



      !  Recover fields from index list
      !
      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_spi(ic,jb)

        p_prog_wtr_new%t_ice(jc,jb)     = tice_new(ic)
        p_prog_wtr_new%h_ice(jc,jb)     = hice_new(ic)
        p_prog_wtr_new%t_snow_si(jc,jb) = tsnow_new(ic)
        p_prog_wtr_new%h_snow_si(jc,jb) = hsnow_new(ic)

      ENDDO  ! ic


      ! Update dynamic seaice index list
      !
      CALL update_index_lists_sea (hice_n      = p_prog_wtr_new%h_ice(:,jb),     &
        &                          idx_lst_sp  = ext_data%atm%idx_lst_sp(:,jb),  &
        &                          sp_count    = ext_data%atm%sp_count(jb),      &
        &                          idx_lst_spw = ext_data%atm%idx_lst_spw(:,jb), &
        &                          spw_count   = ext_data%atm%spw_count(jb),     &
        &                          idx_lst_spi = ext_data%atm%idx_lst_spi(:,jb), &
        &                          spi_count   = ext_data%atm%spi_count(jb)      )


    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_seaice


END MODULE mo_nwp_sfc_interface

