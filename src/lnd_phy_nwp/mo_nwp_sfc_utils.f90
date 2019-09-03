!>
!! Utility routines related to the TERRA surface model
!!
!! @author <NAME>, <AFFILIATION>
!!
!!
!! @par Revision History
!!
!! Initial release by <NAME>, <AFFILIATION>  (YYYY-MM-DD)
!!
!! Modifications by Dmitrii Mironov, DWD (2016-08-02)
!! - Changes related to the use of a rate equation 
!!   for the sea-ice albedo.
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
#if defined __xlC__
@PROCESS SPILL(564)
#endif
MODULE mo_nwp_sfc_utils

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text
  USE mo_exception,           ONLY: finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: tmelt, tf_salt, grav, salinity_fac, rhoh2o
  USE mo_math_constants,      ONLY: dbl_eps, rad2deg
  USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, zml_soil, min_rlcell, dzsoil, &
    &                               MODE_IAU, SSTICE_ANA_CLINC, ALB_SI_MISSVAL
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE sfc_flake_data,         ONLY: tpl_T_r, C_T_min, rflk_depth_bs_ref
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ext_data_init,       ONLY: diagnose_ext_aggr, interpol_monthly_mean, vege_clim
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag, t_lnd_state
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_parallel_config,     ONLY: nproma
  USE mo_grid_config,         ONLY: l_limited_area
  USe mo_extpar_config,       ONLY: itopo, itype_vegetation_cycle
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water, &
    &                               lseaice, llake, lmulti_snow, idiag_snowfrac, ntiles_lnd, &
    &                               lsnowtile, isub_water, isub_seaice, isub_lake,    &
    &                               itype_interception, l2lay_rho_snow, lprog_albsi, itype_trvg, &
                                    itype_snowevap
  USE mo_nwp_tuning_config,   ONLY: tune_minsnowfrac
  USE mo_initicon_config,     ONLY: init_mode_soil, ltile_coldstart, init_mode, lanaread_tseasfc, use_lakeiceana
  USE mo_run_config,          ONLY: msg_level
  USE sfc_terra_init,         ONLY: terra_init
  USE sfc_flake,              ONLY: flake_init
  USE sfc_seaice,             ONLY: seaice_init_nwp, alb_seaice_equil, hice_min, frsi_min, &
    &                               hice_ini_min, hice_ini_max, seaice_coldinit_albsi_nwp
  USE sfc_terra_data,         ONLY: cadp, cf_snow     ! soil and vegetation parameters for TILES
  USE turb_data,              ONLY: c_lnd, c_sea
  USE mo_satad,               ONLY: sat_pres_water, sat_pres_ice, spec_humi
  USE mo_sync,                ONLY: global_sum_array, global_max, global_min
  USE mo_nonhydro_types,      ONLY: t_nh_diag, t_nh_state
  USE mo_dynamics_config,     ONLY: nnow_rcf, nnew_rcf
  USE mtime,                  ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToString

  IMPLICIT NONE

  PRIVATE



#ifdef __SX__
! parameters for loop unrolling
INTEGER, PARAMETER :: nlsoil= 8
INTEGER, PARAMETER :: nlsnow= 2
#endif


  PUBLIC :: nwp_surface_init
  PUBLIC :: diag_snowfrac_tg
  PUBLIC :: aggregate_landvars
  PUBLIC :: update_idx_lists_lnd
  PUBLIC :: update_idx_lists_sea
  PUBLIC :: update_sst_and_seaice
  PUBLIC :: update_ndvi_dependent_fields
  PUBLIC :: init_snowtile_lists
  PUBLIC :: init_sea_lists
  PUBLIC :: aggregate_tg_qvs
  PUBLIC :: copy_lnd_prog_now2new
  PUBLIC :: seaice_albedo_coldstart

CONTAINS

  !>
  !! Initialize soil temperature in boundary zone of LAM domain
  !!
  !! Soil tempertures in the boundary zone of the LAM domain are 
  !! filled with meaningful values. This has no immediate impact on the 
  !! prognostic results. It is, however, necessary in order to minimize 
  !! GRIB truncation errors.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2018-02-23)
  !!
  SUBROUTINE init_lamlatbc_phys (p_patch, p_prog_lnd_now, p_prog_lnd_new, p_lnd_diag)

    TYPE(t_patch)     , INTENT(IN)    :: p_patch       ! patch info
    TYPE(t_lnd_prog)  , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
    TYPE(t_lnd_diag)  , INTENT(INOUT) :: p_lnd_diag

    ! local variables
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc

  !-------------------------------------------------------------------------

    ! only boundary cells
    rl_start = 1
    rl_end   = grf_bdywidth_c

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc= i_startidx, i_endidx
        ! tiled
        p_prog_lnd_now%t_g_t(jc,jb,:)    = tmelt
        p_prog_lnd_now%t_s_t(jc,jb,:)    = tmelt
        p_prog_lnd_now%t_so_t(jc,:,jb,:) = tmelt
        !
        p_prog_lnd_new%t_g_t(jc,jb,:)    = tmelt
        p_prog_lnd_new%t_s_t(jc,jb,:)    = tmelt
        p_prog_lnd_new%t_so_t(jc,:,jb,:) = tmelt
        !
        ! agg
        p_prog_lnd_now%t_g(jc,jb)    = tmelt
        p_prog_lnd_new%t_g(jc,jb)    = tmelt
        p_lnd_diag%t_s(jc,jb)        = tmelt
        p_lnd_diag%t_so(jc,:,jb)     = tmelt
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_lamlatbc_phys

  !-------------------------------------------------------------------------
  !>
  !! Init surface model TERRA, lake model Flake and seaice model
  !!
  !! Init surface model TERRA, lake model Flake and seaice model
  !!
  !! @par Revision History
  !! Initial revision by Ekaterina Machulskaya, DWD (2011-07-??)
  !! Modification by Daniel Reinert, DWD (2011-07-29)
  !! - initialize climatological layer t_so(nlev_soil+1)
  !! Modification by Daniel Reinert, DWD (2012-08-31)
  !! - initialize sea-ice model
  !! Modifications by Dmitrii Mironov, DWD (2016-08-04)
  !! - Call to "seaice_init_nwp" is modified with due regard for
  !!   prognostic treatment of the sea-ice albedo.
  !!
  SUBROUTINE nwp_surface_init( p_patch, ext_data, p_prog_lnd_now,           &
    &                          p_prog_lnd_new, p_prog_wtr_now,              &
    &                          p_prog_wtr_new, p_lnd_diag, p_diag, prm_diag )



    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
    TYPE(t_wtr_prog)     , INTENT(INOUT) :: p_prog_wtr_now, p_prog_wtr_new
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag
    TYPE(t_nh_diag), TARGET,INTENT(inout):: p_diag        !< diag vars
    TYPE(t_nwp_phy_diag),  INTENT(in)    :: prm_diag      !< atm phys vars

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jc,jb,isubs,jk
    LOGICAL :: lsnowtile_warmstart


    REAL(wp) :: t_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_snow_mult_now_t(nproma, 1:nlev_snow+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_rhosnowini_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_g_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_new_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: h_snow_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_mult_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_now_t(nproma, 1:nlev_soil  , p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_new_t(nproma, 1:nlev_soil  , p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_now_t(nproma, 1:nlev_soil, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_new_t(nproma, 1:nlev_soil, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wliq_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wtot_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: dzh_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)

    INTEGER  :: soiltyp_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rootdp_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: plcov_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: z0_t     (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: freshsnow_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: snowfrac_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: sso_sigma_t(nproma,  p_patch%nblks_c, ntiles_total)
    INTEGER  :: lc_class_t (nproma,  p_patch%nblks_c, ntiles_total)


    ! local fields for lake model
    !
    REAL(wp) :: fr_lake       (nproma) ! lake fraction
    REAL(wp) :: depth_lk      (nproma) ! lake depth
    REAL(wp) :: fetch_lk      (nproma) ! wind fetch over lake
    REAL(wp) :: dp_bs_lk      (nproma) ! depth of thermally active layer of bot. sediments.
    REAL(wp) :: t_bs_lk       (nproma) ! clim. temp. at bottom of thermally active layer
                                       ! of sediments
    REAL(wp) :: gamso_lk      (nproma) ! attenuation coefficient of lake water with respect
                                       ! to sol. rad.
    REAL(wp) :: t_snow_lk_now (nproma) ! temperature of snow on lake ice
    REAL(wp) :: h_snow_lk_now (nproma) ! depth of snow on lake ice
    REAL(wp) :: t_ice_now     (nproma) ! lake ice temperature
    REAL(wp) :: h_ice_now     (nproma) ! lake ice depth
    REAL(wp) :: t_mnw_lk_now  (nproma) ! mean temperature of the water column
    REAL(wp) :: t_wml_lk_now  (nproma) ! mixed-layer temperature
    REAL(wp) :: t_bot_lk_now  (nproma) ! temperature at the water-bottom sediment interface
    REAL(wp) :: c_t_lk_now    (nproma) ! shape factor (temp. profile in lake thermocline)
    REAL(wp) :: h_ml_lk_now   (nproma) ! mixed-layer thickness
    REAL(wp) :: t_b1_lk_now   (nproma) ! temperature at the bottom of the upper layer
                                       ! of the sediments
    REAL(wp) :: h_b1_lk_now   (nproma) ! thickness of the upper layer of the sediments
    REAL(wp) :: t_scf_lk_now  (nproma) ! lake surface temperature

    LOGICAL  :: lake_mask(nproma)      ! auxiliary field for re-initialization of non-lake points
    LOGICAL  :: iceana_mask(nproma)    ! mask indicating if sea ice fraction analysis should be used in flake_init

    ! local fields for sea ice model
    !
    REAL(wp) :: frsi     (nproma)   ! sea ice fraction
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level
    REAL(wp) :: albsi_now(nproma)   ! sea-ice albedo at previous time level
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level
    REAL(wp) :: albsi_new(nproma)   ! sea-ice albedo at new time level

    INTEGER  :: icount_flk          ! total number of lake points per block
    !
    INTEGER  :: icount_ice          ! total number of sea-ice points per block

    INTEGER  :: i_count, ic, i_count_snow, isubs_snow, jg, kso
    REAL(wp) :: temp, deglat, deglon
    REAL(wp) :: zfrice_thrhld       ! fraction threshold for creating a sea-ice grid point

  !-------------------------------------------------------------------------

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF (lsnowtile .AND. init_mode == MODE_IAU .AND. .NOT. ltile_coldstart) THEN
      lsnowtile_warmstart = .TRUE.
    ELSE
      lsnowtile_warmstart = .FALSE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,isubs,i_count,i_count_snow,icount_ice, &
!$OMP            icount_flk,temp,ic,isubs_snow,frsi,tice_now,hice_now,               &
!$OMP            tsnow_now,hsnow_now,tice_new,hice_new,tsnow_new,hsnow_new,fr_lake,  &
!$OMP            depth_lk,fetch_lk,dp_bs_lk,t_bs_lk,gamso_lk,t_snow_lk_now,          &
!$OMP            h_snow_lk_now,t_ice_now,h_ice_now,t_mnw_lk_now,t_wml_lk_now,        &
!$OMP            t_bot_lk_now,c_t_lk_now,h_ml_lk_now,t_b1_lk_now,h_b1_lk_now,        &
!$OMP            t_scf_lk_now,zfrice_thrhld,lake_mask,albsi_now,albsi_new,           &
!$OMP            iceana_mask,deglat,deglon), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF (itopo == 1) THEN
        DO isubs = 1, ntiles_total
          DO jc = i_startidx, i_endidx

            ! initialize climatological layer (deepest layer of t_so)
            p_prog_lnd_now%t_so_t(jc,nlev_soil+1,jb,isubs) = ext_data%atm%t_cl(jc,jb)
            p_prog_lnd_new%t_so_t(jc,nlev_soil+1,jb,isubs) = ext_data%atm%t_cl(jc,jb)

          END DO
        END DO
      ENDIF

      ! t_s_t: initialization for open water and sea-ice tiles
      ! proper values are needed to perform surface analysis
      ! open water points: set it to SST
      ! sea-ice points   : set it to tf_salt (salt-water freezing point)
      !
      ! Note that after aggregation, t_s is copied to t_so(1)
      !
      DO ic = 1, ext_data%atm%spw_count(jb)
        jc = ext_data%atm%idx_lst_spw(ic,jb)
        p_prog_lnd_now%t_s_t(jc,jb,isub_water) = p_lnd_diag%t_seasfc(jc,jb)
        p_prog_lnd_new%t_s_t(jc,jb,isub_water) = p_lnd_diag%t_seasfc(jc,jb)
      ENDDO

      DO ic = 1, ext_data%atm%spi_count(jb)
        jc = ext_data%atm%idx_lst_spi(ic,jb)
        p_prog_lnd_now%t_s_t(jc,jb,isub_seaice) = tf_salt
        p_prog_lnd_new%t_s_t(jc,jb,isub_seaice) = tf_salt
      ENDDO


      ! Init t_g_t for sea water points
      !
      DO ic = 1, ext_data%atm%spw_count(jb)
        jc = ext_data%atm%idx_lst_spw(ic,jb)
        temp =  p_lnd_diag%t_seasfc(jc,jb)
        p_prog_lnd_now%t_g_t(jc,jb,isub_water) = temp
        p_prog_lnd_new%t_g_t(jc,jb,isub_water) = temp
        ! includes reduction of saturation pressure due to salt content
        p_lnd_diag%qv_s_t(jc,jb,isub_water)    = salinity_fac * &  
          &       spec_humi(sat_pres_water(temp ),p_diag%pres_sfc(jc,jb) )
      END DO


!---------- Copy input fields for each tile

!----------------------------------
      DO isubs = 1,ntiles_total
!----------------------------------

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          t_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%t_snow_t(jc,jb,isubs)
          t_s_now_t(ic,jb,isubs)             =  p_prog_lnd_now%t_s_t(jc,jb,isubs)
          t_s_new_t(ic,jb,isubs)             =  p_prog_lnd_new%t_s_t(jc,jb,isubs)
          w_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%w_snow_t(jc,jb,isubs)
          rho_snow_now_t(ic,jb,isubs)        =  p_prog_lnd_now%rho_snow_t(jc,jb,isubs)

          sso_sigma_t(ic,jb,isubs)           = ext_data%atm%sso_stdh(jc,jb)
          lc_class_t(ic,jb,isubs)            = ext_data%atm%lc_class_t(jc,jb,isubs)
          freshsnow_t(ic,jb,isubs)           = p_lnd_diag%freshsnow_t(jc,jb,isubs)
          h_snow_t(ic,jb,isubs)              = p_lnd_diag%h_snow_t(jc,jb,isubs)

          soiltyp_t(ic,jb,isubs)             =  ext_data%atm%soiltyp_t(jc,jb,isubs)
          rootdp_t(ic,jb,isubs)              =  ext_data%atm%rootdp_t(jc,jb,isubs)
          plcov_t(ic,jb,isubs)               =  ext_data%atm%plcov_t(jc,jb,isubs)

          IF (isubs > ntiles_lnd) THEN
            z0_t(ic,jb,isubs)                =  prm_diag%gz0_t(jc,jb,isubs-ntiles_lnd)/grav
          ELSE
            z0_t(ic,jb,isubs)                =  prm_diag%gz0_t(jc,jb,isubs)/grav
          ENDIF

          IF (itype_vegetation_cycle >= 2) THEN ! use climatological temperature to specify snow density on glaciers
            t_rhosnowini_t(ic,jb,isubs)      = ext_data%atm%t2m_clim_hc(jc,jb)
          ELSE
            t_rhosnowini_t(ic,jb,isubs)      = t_snow_now_t(ic,jb,isubs)
          ENDIF
        ENDDO

        IF(l2lay_rho_snow .OR. lmulti_snow) THEN
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              rho_snow_mult_now_t(ic,jk,jb,isubs) =  p_prog_lnd_now%rho_snow_mult_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

        ENDIF

        IMSNOWI: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
          DO jk=1,nlev_snow+1
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              t_snow_mult_now_t(ic,jk,jb,isubs)   =  p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

!CDIR UNROLL=nlsnow
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              wliq_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wliq_snow_t    (jc,jk,jb,isubs)
              wtot_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wtot_snow_t    (jc,jk,jb,isubs)
              dzh_snow_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%dzh_snow_t     (jc,jk,jb,isubs)
            ENDDO
          ENDDO

        END IF  IMSNOWI

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            t_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%t_so_t(jc,jk,jb,isubs)
            t_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%t_so_t(jc,jk,jb,isubs)
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil
        DO jk=1,nlev_soil
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            w_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%w_so_t(jc,jk,jb,isubs)
            w_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%w_so_t(jc,jk,jb,isubs)
            w_so_ice_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs)
            w_so_ice_new_t(ic,jk,jb,isubs)      =  p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs)
          ENDDO
        ENDDO

        ! Preliminary diagnosis of snow-cover fraction for initialization of split tile index list

        IF (.NOT. lsnowtile_warmstart) THEN

          CALL diag_snowfrac_tg(                           &
            &  istart = 1, iend = i_count                , & ! start/end indices
            &  lc_class  = lc_class_t        (:,jb,isubs), & ! land-cover class
            &  i_lc_urban = ext_data%atm%i_lc_urban      , & ! land-cover class index for urban areas
            &  t_snow    = t_snow_now_t      (:,jb,isubs), & ! snow temp
            &  t_soiltop = t_s_now_t         (:,jb,isubs), & ! soil top temp
            &  w_snow    = w_snow_now_t      (:,jb,isubs), & ! snow WE
            &  rho_snow  = rho_snow_now_t    (:,jb,isubs), & ! snow density
            &  freshsnow = freshsnow_t       (:,jb,isubs), & ! fresh snow fraction
            &  sso_sigma = sso_sigma_t       (:,jb,isubs), & ! sso stdev
            &  z0        = z0_t              (:,jb,isubs), & ! vegetation roughness length
            &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
            &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp

!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)  = snowfrac_t(ic,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
          ENDDO

        ENDIF

      ENDDO

      ! create index lists for snow tiles (first call)
      IF(lsnowtile .AND. .NOT. lsnowtile_warmstart) THEN ! snow is considered as separate tiles
        DO isubs = 1, ntiles_lnd
          isubs_snow = isubs + ntiles_lnd

          CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
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
                                   snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

        END DO
      END IF


      DO isubs = 1,ntiles_total

        i_count = ext_data%atm%lp_count_t(jb,isubs)
        CALL terra_init(                                          &
        &  init_mode         = init_mode_soil                   , & ! coldstart/warmstart/warmstart with snow increments
        &  nvec              = nproma                           , & ! array dimensions
        &  ivstart=1,   ivend= i_count                          , & ! optional start/end indicies
        &  iblock            = jb                               , & ! actual block
        &  ke_soil=nlev_soil-1, ke_snow=nlev_snow               , & ! without lowermost (climat.) soil layer
        &   zmls             = zml_soil                         , & ! processing soil level structure
        &  soiltyp_subs      = soiltyp_t(:,jb,isubs)            , & ! type of the soil (keys 0-9)  --
        &  rootdp            = rootdp_t(:,jb,isubs)             , & ! depth of the roots                ( m  )
        &  plcov             = plcov_t(:,jb,isubs)              , & ! surface fraction covered by plants ( -  )
        &  t_snow_now        = t_snow_now_t(:,jb,isubs)         , & ! temperature of the snow-surface   (  K  )
        &  t_snow_mult_now   = t_snow_mult_now_t(:,:,jb,isubs)  , & ! temperature of the snow-surface   (  K  )
        &  t_rhosnowini      = t_rhosnowini_t(:,jb,isubs)       , & ! temperature used for snow density initialization (  K  )
        &  t_s_now           = t_s_now_t(:,jb,isubs)            , & ! temperature of the ground surface (  K  )
        &  t_s_new           = t_s_new_t(:,jb,isubs)            , & ! temperature of the ground surface (  K  )
        &  w_snow_now        = w_snow_now_t(:,jb,isubs)         , & ! water content of snow             (m H2O)
        &  h_snow            = h_snow_t(:,jb,isubs)             , & ! snow depth                        (m H2O)
        &  rho_snow_now      = rho_snow_now_t(:,jb,isubs)       , & ! snow density                      (kg/m**3)
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,jb,isubs), & ! snow density                      (kg/m**3)
        &  t_so_now          = t_so_now_t(:,:,jb,isubs)         , & ! soil temperature (main level)     (  K  )
        &  t_so_new          = t_so_new_t(:,:,jb,isubs)         , & ! soil temperature (main level)     (  K  )
        &  w_so_now          = w_so_now_t(:,:,jb,isubs)         , & ! total water content (ice + liquid water)     (m H20)
        &  w_so_new          = w_so_new_t(:,:,jb,isubs)         , & ! total water content (ice + liquid water)     (m H20)
        &  w_so_ice_now      = w_so_ice_now_t(:,:,jb,isubs)     , & ! ice content                       (m H20)
        &  w_so_ice_new      = w_so_ice_new_t(:,:,jb,isubs)     , & ! ice content                       (m H20)
        &  wliq_snow_now     = wliq_snow_now_t(:,:,jb,isubs)    , & ! liquid water content in the snow  (m H2O)
        &  wtot_snow_now     = wtot_snow_now_t(:,:,jb,isubs)    , & ! total (liquid + solid) water content of snow (m H2O)
        &  dzh_snow_now      = dzh_snow_now_t(:,:,jb,isubs)       ) ! layer thickness between half levels in snow  (  m  )

        IF (.NOT. lsnowtile_warmstart .OR. isubs > ntiles_lnd) THEN

          CALL diag_snowfrac_tg(                           &
            &  istart = 1, iend = i_count                , & ! start/end indices
            &  lc_class  = lc_class_t        (:,jb,isubs), & ! land-cover class
            &  i_lc_urban = ext_data%atm%i_lc_urban      , & ! land-cover class index for urban areas
            &  t_snow    = t_snow_now_t      (:,jb,isubs), & ! snow temp
            &  t_soiltop = t_s_now_t         (:,jb,isubs), & ! soil top temp
            &  w_snow    = w_snow_now_t      (:,jb,isubs), & ! snow WE
            &  rho_snow  = rho_snow_now_t    (:,jb,isubs), & ! snow depth
            &  freshsnow = freshsnow_t       (:,jb,isubs), & ! fresh snow fraction
            &  sso_sigma = sso_sigma_t       (:,jb,isubs), & ! sso stdev
            &  z0        = z0_t              (:,jb,isubs), & ! vegetation roughness length
            &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
            &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp

        ENDIF

!  Recover fields from index list
!
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          p_prog_lnd_now%t_snow_t(jc,jb,isubs)   = t_snow_now_t(ic,jb,isubs)
          p_prog_lnd_now%t_s_t(jc,jb,isubs)      = t_s_now_t(ic,jb,isubs)
          p_prog_lnd_new%t_s_t(jc,jb,isubs)      = t_s_new_t(ic,jb,isubs)
          p_prog_lnd_now%w_snow_t(jc,jb,isubs)   = w_snow_now_t(ic,jb,isubs)
          p_lnd_diag%h_snow_t(jc,jb,isubs)       = h_snow_t(ic,jb,isubs)
          p_prog_lnd_now%rho_snow_t(jc,jb,isubs) = rho_snow_now_t(ic,jb,isubs)
          IF (.NOT. lsnowtile_warmstart) THEN
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)  = snowfrac_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_t(jc,jb,isubs)     = snowfrac_t(ic,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
            p_prog_lnd_new%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
          ENDIF
          IF (isubs > ntiles_lnd) THEN
            p_lnd_diag%snowfrac_lcu_t(jc,jb,isubs) = MAX(snowfrac_t(ic,jb,isubs),p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
            p_lnd_diag%snowfrac_lcu_t(jc,jb,isubs-ntiles_lnd) = p_lnd_diag%snowfrac_lcu_t(jc,jb,isubs)
          ENDIF
        ENDDO

        IF (lsnowtile .AND. isubs > ntiles_lnd .AND. .NOT. lsnowtile_warmstart) THEN
          ! copy snowfrac_t to snow-free tile (needed for index list computation)
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs-ntiles_lnd) = p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)
          ENDDO

          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = p_prog_lnd_now%w_snow_t(jc,jb,isubs)/            &
                                                   MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
            p_lnd_diag%h_snow_t(jc,jb,isubs)     = p_lnd_diag%h_snow_t(jc,jb,isubs)/                &
                                                   MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
          ENDDO

        ENDIF

        IF (itype_snowevap == 3) THEN ! set snow age to upper limit of 365 days on grid points with glaciers,
          DO ic = 1, i_count          ! and reset hsnow_max to 40 m in agreement with what is done in initicon
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            IF (ext_data%atm%lc_class_t(jc,jb,isubs) == ext_data%atm%i_lc_snow_ice) &
              p_lnd_diag%snow_age(jc,jb)  = 365._wp
              p_lnd_diag%hsnow_max(jc,jb) = MIN(40._wp,p_lnd_diag%hsnow_max(jc,jb))
          ENDDO
        ENDIF

        IF(l2lay_rho_snow .OR. lmulti_snow) THEN
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              p_prog_lnd_now%rho_snow_mult_t(jc,jk,jb,isubs) = rho_snow_mult_now_t(ic,jk,jb,isubs)
            ENDDO
          ENDDO
        ENDIF

        IMSNOWO: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
          DO jk=1,nlev_snow+1
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) =  t_snow_mult_now_t(ic,jk,jb,isubs)
            ENDDO
          ENDDO

!CDIR UNROLL=nlsnow
          DO jk=1,nlev_snow
            IF (lsnowtile .AND. isubs > ntiles_lnd .AND. .NOT. lsnowtile_warmstart) THEN
              DO ic = 1, i_count
                jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = wliq_snow_now_t(ic,jk,jb,isubs)/ &
                                             MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = wtot_snow_now_t(ic,jk,jb,isubs)/ &
                                             MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
                p_prog_lnd_now%dzh_snow_t(jc,jk,jb,isubs)  = dzh_snow_now_t(ic,jk,jb,isubs)/ &
                                            MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
              ENDDO
            ELSE
              DO ic = 1, i_count
                jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = wliq_snow_now_t(ic,jk,jb,isubs)
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = wtot_snow_now_t(ic,jk,jb,isubs)
                p_prog_lnd_now%dzh_snow_t(jc,jk,jb,isubs)  = dzh_snow_now_t(ic,jk,jb,isubs)
              ENDDO
            ENDIF
          ENDDO

        END IF  IMSNOWO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%t_so_t(jc,jk,jb,isubs) = t_so_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%t_so_t(jc,jk,jb,isubs) = t_so_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil
        DO jk=1,nlev_soil
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%w_so_t(jc,jk,jb,isubs) = w_so_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_t(jc,jk,jb,isubs) = w_so_new_t(ic,jk,jb,isubs)
            p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO
      END DO ! isubs


      !===================================================!
      !                                                   !
      !  Warm-start initialization for sea-ice and lake   !
      !                                                   !
      !===================================================!

      !
      ! Warm-start initialization for fresh water lake parameterization
      ! This initialization is performed irrespective of a cold-start initialization.
      !
      IF (llake) THEN

        icount_flk = ext_data%atm%fp_count(jb) ! number of lake points in block jb

        ! Collect data for lake points in 1D-arrays
        DO ic = 1, icount_flk

          jc = ext_data%atm%idx_lst_fp(ic,jb)

          fr_lake      (ic) = ext_data%atm%frac_t(jc,jb,isub_lake)
          depth_lk     (ic) = ext_data%atm%depth_lk   (jc,jb)
          frsi         (ic) = p_lnd_diag%fr_seaice(jc,jb)
          fetch_lk     (ic) = ext_data%atm%fetch_lk   (jc,jb)
          dp_bs_lk     (ic) = ext_data%atm%dp_bs_lk   (jc,jb)
          t_bs_lk      (ic) = ext_data%atm%t_bs_lk    (jc,jb)
          gamso_lk     (ic) = ext_data%atm%gamso_lk   (jc,jb)
          t_snow_lk_now(ic) = p_prog_wtr_now%t_snow_lk(jc,jb)
          h_snow_lk_now(ic) = p_prog_wtr_now%h_snow_lk(jc,jb)
          t_ice_now    (ic) = p_prog_wtr_now%t_ice    (jc,jb)
          h_ice_now    (ic) = p_prog_wtr_now%h_ice    (jc,jb)
          t_mnw_lk_now (ic) = p_prog_wtr_now%t_mnw_lk (jc,jb)
          t_wml_lk_now (ic) = p_prog_wtr_now%t_wml_lk (jc,jb)
          t_bot_lk_now (ic) = p_prog_wtr_now%t_bot_lk (jc,jb)
          c_t_lk_now   (ic) = p_prog_wtr_now%c_t_lk   (jc,jb)
          h_ml_lk_now  (ic) = p_prog_wtr_now%h_ml_lk  (jc,jb)
          t_b1_lk_now  (ic) = p_prog_wtr_now%t_b1_lk  (jc,jb)
          h_b1_lk_now  (ic) = p_prog_wtr_now%h_b1_lk  (jc,jb)
          t_scf_lk_now (ic) = p_prog_lnd_now%t_g_t    (jc,jb,isub_lake)

          ! test implementation: compute approximate mask for Great Lakes
          deglat = p_patch%cells%center(jc,jb)%lat * rad2deg
          deglon = p_patch%cells%center(jc,jb)%lon * rad2deg
          IF (deglat >= 41._wp .AND. deglat <= 49._wp .AND. deglon >= -92._wp .AND. deglon <= -76._wp) THEN
            iceana_mask(ic) = .TRUE.
            IF (deglon > -85._wp .AND. deglat > 46._wp) iceana_mask(ic) = .FALSE.
            IF (deglon < -88._wp .AND. deglat < 46._wp) iceana_mask(ic) = .FALSE.
          ELSE
            iceana_mask(ic) = .FALSE.
          ENDIF
          IF (.NOT. (use_lakeiceana .AND. lanaread_tseasfc(jg)) ) iceana_mask(ic) = .FALSE.
        ENDDO


        CALL flake_init (nflkgb     = icount_flk,             & !in
          &         use_iceanalysis = iceana_mask  (:),       & !in
          &              fr_lake    = fr_lake      (:),       & !in
          &              depth_lk   = depth_lk     (:),       & !in
          &              fr_ice     = frsi         (:),       & !in
          &              fetch_lk   = fetch_lk     (:),       & !inout
          &              dp_bs_lk   = dp_bs_lk     (:),       & !inout
          &              t_bs_lk    = t_bs_lk      (:),       & !inout
          &              gamso_lk   = gamso_lk     (:),       & !inout
          &              t_snow_p   = t_snow_lk_now(:),       & !inout
          &              h_snow_p   = h_snow_lk_now(:),       & !inout
          &              t_ice_p    = t_ice_now    (:),       & !inout
          &              h_ice_p    = h_ice_now    (:),       & !inout
          &              t_mnw_lk_p = t_mnw_lk_now (:),       & !inout
          &              t_wml_lk_p = t_wml_lk_now (:),       & !inout
          &              t_bot_lk_p = t_bot_lk_now (:),       & !inout
          &              c_t_lk_p   = c_t_lk_now   (:),       & !inout
          &              h_ml_lk_p  = h_ml_lk_now  (:),       & !inout
          &              t_b1_lk_p  = t_b1_lk_now  (:),       & !inout
          &              h_b1_lk_p  = h_b1_lk_now  (:),       & !inout
          &              t_scf_lk_p = t_scf_lk_now (:)        )


        !  Recover fields from index list
        !
        DO ic = 1, icount_flk

          jc = ext_data%atm%idx_lst_fp(ic,jb)

          ext_data%atm%fetch_lk(jc,jb) = fetch_lk(ic)
          ext_data%atm%dp_bs_lk(jc,jb) = dp_bs_lk(ic)
          ext_data%atm%t_bs_lk (jc,jb) = t_bs_lk (ic)
          ext_data%atm%gamso_lk(jc,jb) = gamso_lk(ic)

          p_prog_wtr_now%t_snow_lk(jc,jb) = t_snow_lk_now(ic)
          p_prog_wtr_now%h_snow_lk(jc,jb) = h_snow_lk_now(ic)
          p_prog_wtr_now%t_ice    (jc,jb) = t_ice_now    (ic)
          p_prog_wtr_now%h_ice    (jc,jb) = h_ice_now    (ic)
          p_prog_wtr_now%t_mnw_lk (jc,jb) = t_mnw_lk_now (ic)
          p_prog_wtr_now%t_wml_lk (jc,jb) = t_wml_lk_now (ic)
          p_prog_wtr_now%t_bot_lk (jc,jb) = t_bot_lk_now (ic)
          p_prog_wtr_now%c_t_lk   (jc,jb) = c_t_lk_now   (ic)
          p_prog_wtr_now%h_ml_lk  (jc,jb) = h_ml_lk_now  (ic)
          p_prog_wtr_now%t_b1_lk  (jc,jb) = t_b1_lk_now  (ic)
          p_prog_wtr_now%h_b1_lk  (jc,jb) = h_b1_lk_now  (ic)

          p_prog_lnd_now%t_g_t(jc,jb,isub_lake) = t_scf_lk_now(ic)

          ! for consistency, set 
          ! t_so(0) = t_wml_lk       mixed-layer temperature (273.15K if the lake is frozen)
          p_prog_lnd_now%t_s_t(jc,jb,isub_lake) = p_prog_wtr_now%t_wml_lk (jc,jb)
          p_prog_lnd_new%t_s_t(jc,jb,isub_lake) = p_prog_lnd_now%t_s_t(jc,jb,isub_lake)


          ! In addition, initialize prognostic Flake fields at time step 'new'
          p_prog_wtr_new%t_snow_lk(jc,jb) = t_snow_lk_now(ic)
          p_prog_wtr_new%h_snow_lk(jc,jb) = h_snow_lk_now(ic)
          p_prog_wtr_new%t_ice    (jc,jb) = t_ice_now    (ic)
          p_prog_wtr_new%h_ice    (jc,jb) = h_ice_now    (ic)
          p_prog_wtr_new%t_mnw_lk (jc,jb) = t_mnw_lk_now (ic)
          p_prog_wtr_new%t_wml_lk (jc,jb) = t_wml_lk_now (ic)
          p_prog_wtr_new%t_bot_lk (jc,jb) = t_bot_lk_now (ic)
          p_prog_wtr_new%c_t_lk   (jc,jb) = c_t_lk_now   (ic)
          p_prog_wtr_new%h_ml_lk  (jc,jb) = h_ml_lk_now  (ic)
          p_prog_wtr_new%t_b1_lk  (jc,jb) = t_b1_lk_now  (ic)
          p_prog_wtr_new%h_b1_lk  (jc,jb) = h_b1_lk_now  (ic)

          ! keep fr_seaice synchronized with h_ice
          ! i.e. set fr_seaice=1 for frozen lakes
          p_lnd_diag%fr_seaice(jc,jb) = MERGE(1.0_wp, 0.0_wp, &
            &                           p_prog_wtr_now%h_ice(jc,jb)>0._wp)
        ENDDO  ! ic

        ! Re-Initialize lake-specific fields
        ! t_wml_lk
        ! t_mnw_lk
        ! t_bot_lk
        ! c_t_lk
        ! h_ml_lk
        ! h_snow_lk
        ! t_snow_lk
        ! t_bl_lk
        ! h_bl_lk
        ! for non-lake points, because values might be inconsistent due to GRIB-packing
        !
        ! Create lake-mask, in order to re-initialize only non-lake points
        lake_mask(i_startidx:i_endidx) = .FALSE.
        ! set lake-mask to .TRUE. for lake points
        DO ic = 1, icount_flk
          jc = ext_data%atm%idx_lst_fp(ic,jb)
          lake_mask(jc) = .TRUE.
        ENDDO

        DO jc = i_startidx, i_endidx
          IF (.NOT. lake_mask(jc)) THEN
            ! now
            p_prog_wtr_now%t_wml_lk (jc,jb) = tpl_T_r        ! temperature of maximum density
                                                             ! of fresh water
            p_prog_wtr_now%t_mnw_lk (jc,jb) = tpl_T_r
            p_prog_wtr_now%t_bot_lk (jc,jb) = tpl_T_r
            p_prog_wtr_now%c_t_lk   (jc,jb) = C_T_min        ! minimum value
            p_prog_wtr_now%h_ml_lk  (jc,jb) = 0._wp
            p_prog_wtr_now%h_snow_lk(jc,jb) = 0._wp
            p_prog_wtr_now%t_snow_lk(jc,jb) = tmelt          ! fresh water freezing point
            p_prog_wtr_now%t_b1_lk  (jc,jb) = tpl_T_r        ! reference value, bottom-sediment
                                                             ! module is switched off
            p_prog_wtr_now%h_b1_lk  (jc,jb) = rflk_depth_bs_ref ! reference value, bottom-sediment
                                                                ! is switched off
            !
            ! new
            p_prog_wtr_new%t_wml_lk (jc,jb) = tpl_T_r        ! temperature of maximum density
                                                             ! of fresh water
            p_prog_wtr_new%t_mnw_lk (jc,jb) = tpl_T_r
            p_prog_wtr_new%t_bot_lk (jc,jb) = tpl_T_r
            p_prog_wtr_new%c_t_lk   (jc,jb) = C_T_min        ! minimum value
            p_prog_wtr_new%h_ml_lk  (jc,jb) = 0._wp
            p_prog_wtr_new%h_snow_lk(jc,jb) = 0._wp
            p_prog_wtr_new%t_snow_lk(jc,jb) = tmelt          ! fresh water freezing point
            p_prog_wtr_new%t_b1_lk  (jc,jb) = tpl_T_r        ! reference value, bottom-sediment
                                                             ! module is switched off
            p_prog_wtr_new%h_b1_lk  (jc,jb) = rflk_depth_bs_ref ! reference value, bottom-sediment
                                                                ! is switched off
          ENDIF
        ENDDO

      ENDIF  ! llake


      !
      ! Warm-start for sea-ice parameterization scheme
      ! This initialization is performed irrespectively of a coldstart initialization.
      ! At this point it is assumed that both the ice thickness and the ice surface
      ! temperature are available from the previous ICON run or any other educated guess.
      !
      IF ( lseaice ) THEN

        icount_ice = ext_data%atm%spi_count(jb) ! number of sea-ice points in block jb

        DO ic = 1, icount_ice

          jc = ext_data%atm%idx_lst_spi(ic,jb)

          frsi     (ic) = p_lnd_diag%fr_seaice(jc,jb)
          tice_now (ic) = p_prog_wtr_now%t_ice(jc,jb)
          hice_now (ic) = p_prog_wtr_now%h_ice(jc,jb)
          tsnow_now(ic) = p_prog_wtr_now%t_snow_si(jc,jb)
          hsnow_now(ic) = p_prog_wtr_now%h_snow_si(jc,jb)
          albsi_now(ic) = p_prog_wtr_now%alb_si(jc,jb)
        ENDDO  ! jc


        CALL seaice_init_nwp ( icount_ice, frsi,                                    & ! in
          &                    tice_now, hice_now, tsnow_now, hsnow_now, albsi_now, & ! inout
          &                    tice_new, hice_new, tsnow_new, hsnow_new, albsi_new  ) ! inout


        !  Recover fields from index list
        !
        DO ic = 1, icount_ice
          jc = ext_data%atm%idx_lst_spi(ic,jb)

          ! fields at time level now may have changed, potentially!
          p_prog_wtr_now%t_ice(jc,jb)    = tice_now(ic)
          p_prog_wtr_now%h_ice(jc,jb)    = hice_now(ic)
          p_prog_wtr_now%t_snow_si(jc,jb)= tsnow_now(ic)
          p_prog_wtr_now%h_snow_si(jc,jb)= hsnow_now(ic)
          p_prog_wtr_now%alb_si(jc,jb)   = albsi_now(ic)

          p_prog_wtr_new%t_ice(jc,jb)    = tice_new(ic)
          p_prog_wtr_new%h_ice(jc,jb)    = hice_new(ic)
          p_prog_wtr_new%t_snow_si(jc,jb)= tsnow_new(ic)
          p_prog_wtr_new%h_snow_si(jc,jb)= hsnow_new(ic)
          IF (lprog_albsi) THEN
            p_prog_wtr_new%alb_si(jc,jb)   = albsi_new(ic)
          ENDIF

          p_prog_lnd_now%t_g_t(jc,jb,isub_seaice) =  tice_now(ic)
          p_prog_lnd_new%t_g_t(jc,jb,isub_seaice) =  tice_now(ic)
          p_lnd_diag%qv_s_t(jc,jb,isub_seaice)    = spec_humi(sat_pres_ice(tice_now(ic)),&
          &                                   p_diag%pres_sfc(jc,jb) )
        ENDDO  ! ic


        ! Re-initialize sea-ice specific fields h_snow_si, t_snow_si for non sea-ice points,
        ! because values might be inconsistent due to GRIB-packing
        !
        IF ( ntiles_total == 1 ) THEN  ! no tile approach
          zfrice_thrhld = 0.5_wp
        ELSE
          zfrice_thrhld = frsi_min
        ENDIF
        DO jc = i_startidx, i_endidx
          IF (p_lnd_diag%fr_seaice(jc,jb) < zfrice_thrhld) THEN   ! non-seaice point
            ! now
            p_prog_wtr_now%h_snow_si(jc,jb) = 0._wp
            p_prog_wtr_now%t_snow_si(jc,jb) = tmelt   ! snow over ice is not treated explicitly
            !
            ! new
            p_prog_wtr_new%h_snow_si(jc,jb) = 0._wp
            p_prog_wtr_new%t_snow_si(jc,jb) = tmelt   ! snow over ice is not treated explicitly
          ENDIF
        ENDDO  ! jc

      ENDIF  ! lseaice


      ! Re-initialize h_ice and t_ice which are used by both the sea-ice and lake scheme.
      ! 4 cases need to be distinguished in order to make sure, that none of the prognostic
      ! points is overwritten.
      IF      (.NOT. lseaice .AND. .NOT. llake) THEN

        ! Re-initialize all points
        !
        DO jc = i_startidx, i_endidx
          ! now
          p_prog_wtr_now%h_ice  (jc,jb) = 0._wp
          p_prog_wtr_now%t_ice  (jc,jb) = tmelt   ! fresh water freezing point
          !
          ! new
          p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
          p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
        ENDDO

      ELSE IF (      lseaice .AND. .NOT. llake) THEN

        ! Re-initialize all non-sea-ice points
        !
        IF ( ntiles_total == 1 ) THEN  ! no tile approach
          zfrice_thrhld = 0.5_wp
        ELSE
          zfrice_thrhld = frsi_min
        ENDIF
        DO jc = i_startidx, i_endidx
          IF (p_lnd_diag%fr_seaice(jc,jb) < zfrice_thrhld) THEN   ! non-seaice point
            ! now
            p_prog_wtr_now%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_now%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            ! new
            p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            IF (lprog_albsi) THEN
              p_prog_wtr_new%alb_si (jc,jb) = ALB_SI_MISSVAL
              p_prog_wtr_now%alb_si (jc,jb) = ALB_SI_MISSVAL
            ENDIF
          ENDIF
        ENDDO

      ELSE IF (.NOT. lseaice .AND.       llake) THEN

        ! Re-initialize all non-lake points
        !
        DO jc = i_startidx, i_endidx
          IF (.NOT. lake_mask(jc)) THEN   ! non-lake point
            ! now
            p_prog_wtr_now%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_now%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            ! new
            p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
          ENDIF
        ENDDO

      ELSE IF (      lseaice .AND.       llake) THEN

        ! Re-initialize all points which are neither lake nor sea-ice points
        !
        IF ( ntiles_total == 1 ) THEN  ! no tile approach
          zfrice_thrhld = 0.5_wp
        ELSE
          zfrice_thrhld = frsi_min
        ENDIF
        DO jc = i_startidx, i_endidx
          IF ( .NOT. lake_mask(jc) .AND. (p_lnd_diag%fr_seaice(jc,jb) < zfrice_thrhld) ) THEN
            ! now
            p_prog_wtr_now%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_now%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            ! new
            p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            IF (lprog_albsi) THEN
              p_prog_wtr_new%alb_si (jc,jb) = ALB_SI_MISSVAL
              p_prog_wtr_now%alb_si (jc,jb) = ALB_SI_MISSVAL
            ENDIF
          ENDIF
        ENDDO
      ENDIF



      IF(lsnowtile) THEN ! snow is considered as separate tiles
        DO isubs = 1, ntiles_lnd

          isubs_snow = isubs + ntiles_lnd

          CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
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
                                   snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

          i_count = ext_data%atm%gp_count_t(jb,isubs)
          i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count_snow
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

            ! reset field for actual snow-cover for grid points / land-cover classes for which there
            ! are seperate snow-free and snow-covered tiles
            p_lnd_diag%snowfrac_t(jc,jb,isubs)   = 0._wp
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = 0._wp
            p_lnd_diag%h_snow_t(jc,jb,isubs)     = 0._wp
            p_prog_lnd_now%t_snow_t(jc,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)    = p_prog_lnd_now%t_s_t(jc,jb,isubs)

            ! copy rho_snow in order to get the right tile average of snow density
            p_prog_lnd_now%rho_snow_t(jc,jb,isubs) = p_prog_lnd_now%rho_snow_t(jc,jb,isubs_snow)

            ! to prevent numerical stability problems, we require at least 1 cm of snow in order to
            ! have a snow-cover fraction of 1 on snow tiles (not critical for the single-layer
            ! snow scheme, but the multi-layer snow model becomes numerically unstable within a few
            ! time steps when associating traces of snow with a snow-cover fraction of 1)
            p_lnd_diag%snowfrac_t(jc,jb,isubs_snow) = MIN(1._wp,p_lnd_diag%h_snow_t(jc,jb,isubs_snow)*100._wp)

            ! Rediagnose t_g according to the modified snow-cover fraction
            p_prog_lnd_now%t_g_t(jc,jb,isubs_snow) =  &
              p_lnd_diag%snowfrac_t(jc,jb,isubs_snow) * p_prog_lnd_now%t_snow_t(jc,jb,isubs_snow) + &
              (1._wp-p_lnd_diag%snowfrac_t(jc,jb,isubs_snow))*p_prog_lnd_now%t_s_t(jc,jb,isubs_snow)
            ! Rediagnose qv_s_t because of the rediagnosed t_g_t
            p_lnd_diag%qv_s_t(jc,jb,isubs_snow)  =                                                  &
              &                    spec_humi(sat_pres_ice  (p_prog_lnd_now%t_g_t(jc,jb,isubs_snow)),&
              &                    p_diag%pres_sfc(jc,jb) )
             IF (lmulti_snow) THEN
               p_prog_lnd_now%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
             ENDIF

          END DO

          IF (lmulti_snow) THEN
!CDIR UNROLL=nlsnow
            DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
              DO ic = 1, i_count_snow
                jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
                p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = 0._wp
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = 0._wp
                p_prog_lnd_now%dzh_snow_t (jc,jk,jb,isubs) = 0._wp
              ENDDO
            ENDDO
          ENDIF

        END DO
      END IF


      ! Remove snow and w_i on non-existing grid points. This has no impact on the prognostic
      ! results but is needed in order to have meaningful data on the tile-based fields 
      DO isubs = 1, ntiles_total
        DO jc = i_startidx, i_endidx
          IF (ext_data%atm%frac_t(jc,jb,isubs) < 1.e-10_wp) THEN
            p_lnd_diag%h_snow_t(jc,jb,isubs)     = 0._wp
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = 0._wp
            p_lnd_diag%snowfrac_t(jc,jb,isubs)   = 0._wp
            p_prog_lnd_now%w_i_t(jc,jb,isubs)    = 0._wp
          ENDIF
          IF (ext_data%atm%lc_frac_t(jc,jb,isubs) < 1.e-10_wp) THEN
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)   = 0._wp
          ENDIF
        ENDDO
      END DO


    ENDDO  ! jb loop
!$OMP END DO
!$OMP END PARALLEL

    ! Aggregate t_g and qv_s
    ! Loop over all points (land AND water points)
    ! Aggregation has been moved to the end of the subroutine (PR)
    !
    CALL aggregate_tg_qvs( p_patch, ext_data, p_prog_lnd_now , &
     &                           p_lnd_diag )
    p_prog_lnd_new%t_g(:,:)  = p_prog_lnd_now%t_g(:,:)


    ! limited area mode, only
    ! Soil and surface tempertures in the boundary zone of the LAM 
    ! domain are filled with meaningful values
    IF (l_limited_area .AND. jg==1) THEN
      CALL init_lamlatbc_phys(p_patch, p_prog_lnd_now, p_prog_lnd_new, p_lnd_diag)
    ENDIF

  END SUBROUTINE nwp_surface_init


!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  !! Agregate t_g and qv_s_t
  !!
  !! @par Revision History
  !! Initial revision by P Ripodas, DWD (2012-12)
  !!
  !! Segregated from nwp_surface_init, to avoid recoding again
  !-------------------------------------------------------------------------

  SUBROUTINE aggregate_tg_qvs( p_patch, ext_data, p_prog_lnd_now, &
    &                           p_lnd_diag )



    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local :

    INTEGER :: jc,jb,isubs
    REAL(wp) :: t_g_s(nproma), qv_s_s(nproma)
!-------------------------------------------------------------------------

    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,isubs,i_startidx,i_endidx,t_g_s,qv_s_s)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


       IF (ntiles_total == 1) THEN
         DO jc = i_startidx, i_endidx
           p_prog_lnd_now%t_g(jc,jb)  = p_prog_lnd_now%t_g_t(jc,jb,1)
           p_lnd_diag%qv_s(jc,jb)     = p_lnd_diag%qv_s_t(jc,jb,1)
         ENDDO
       ELSE ! aggregate fields over tiles
         t_g_s(:)  =  0._wp
         qv_s_s(:) =  0._wp
         DO isubs = 1,ntiles_total+ntiles_water
           DO jc = i_startidx, i_endidx
             t_g_s(jc) = t_g_s(jc) + ext_data%atm%frac_t(jc,jb,isubs)  &
               &       * p_prog_lnd_now%t_g_t(jc,jb,isubs)**4
             qv_s_s(jc) = qv_s_s(jc) + ext_data%atm%frac_t(jc,jb,isubs) &
               &       * p_lnd_diag%qv_s_t(jc,jb,isubs)
           ENDDO
         ENDDO
         DO jc = i_startidx, i_endidx
           p_prog_lnd_now%t_g(jc,jb)  = SQRT(SQRT(t_g_s(jc)))
           p_lnd_diag%qv_s(jc,jb)     = qv_s_s(jc)
         ENDDO

       ENDIF    ! with or without tiles

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE aggregate_tg_qvs

!-------------------------------------------------------------------------

  SUBROUTINE aggregate_landvars( p_patch, ext_data, lnd_prog, lnd_diag)

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch       !< grid/patch info
    TYPE(t_external_data),       INTENT(in)   :: ext_data      !< external data
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog      !< prog vars for sfc
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc

    ! Local scalars:
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains
    INTEGER :: jc, jb, jk, isubs

    REAL(wp) :: tilefrac ! fractional area covered by tile

    LOGICAL :: lmask(nproma)  ! mask array (TRUE for landpoint)
    INTEGER :: icount         ! index list length per block
    INTEGER :: ic
    INTEGER :: styp           ! soil type index

    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,jk,tilefrac,lmask,icount,styp) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


      IF (ntiles_total == 1) THEN  ! just copy prognostic variables from tile 1
                                   ! to diagnostic aggregated variable
        DO jc = i_startidx, i_endidx
          lnd_diag%t_snow   (jc,jb) = lnd_prog%t_snow_t   (jc,jb,1)
          lnd_diag%t_s      (jc,jb) = lnd_prog%t_s_t      (jc,jb,1)
          lnd_diag%w_snow   (jc,jb) = lnd_prog%w_snow_t   (jc,jb,1)
          lnd_diag%rho_snow (jc,jb) = lnd_prog%rho_snow_t (jc,jb,1)
          lnd_diag%w_i      (jc,jb) = lnd_prog%w_i_t      (jc,jb,1)
          lnd_diag%h_snow   (jc,jb) = lnd_diag%h_snow_t   (jc,jb,1)
          lnd_diag%freshsnow(jc,jb) = lnd_diag%freshsnow_t(jc,jb,1)
          lnd_diag%snowfrac (jc,jb) = lnd_diag%snowfrac_t (jc,jb,1)
          lnd_diag%runoff_s (jc,jb) = lnd_diag%runoff_s_t (jc,jb,1)
          lnd_diag%runoff_g (jc,jb) = lnd_diag%runoff_g_t (jc,jb,1)
          lnd_diag%rstom    (jc,jb) = lnd_diag%rstom_t    (jc,jb,1)

          IF(lmulti_snow) THEN
            lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,1)
          ENDIF
        ENDDO

        IF (itype_interception == 2) THEN
          DO jc = i_startidx, i_endidx
            lnd_diag%w_p    (jc,jb) = lnd_prog%w_p_t      (jc,jb,1)
            lnd_diag%w_s    (jc,jb) = lnd_prog%w_s_t      (jc,jb,1)
          ENDDO
        ENDIF

        IF (itype_trvg == 3) THEN
          DO jc = i_startidx, i_endidx
            lnd_diag%plantevap(jc,jb) = lnd_diag%plantevap_t(jc,jb,1)
          ENDDO
        ENDIF

        DO jk=1,nlev_soil
          DO jc = i_startidx, i_endidx
            lnd_diag%t_so    (jc,jk+1,jb) = lnd_prog%t_so_t    (jc,jk+1,jb,1)
            lnd_diag%w_so    (jc,jk,  jb) = lnd_prog%w_so_t    (jc,jk,  jb,1)
            lnd_diag%w_so_ice(jc,jk,  jb) = lnd_prog%w_so_ice_t(jc,jk,  jb,1)
          ENDDO
        ENDDO

        IF (l2lay_rho_snow .OR. lmulti_snow) THEN
          DO jk=1,nlev_snow
            DO jc = i_startidx, i_endidx
              lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_prog%rho_snow_mult_t(jc,jk,jb,1)
            ENDDO
          ENDDO
        ENDIF

        IF (lmulti_snow) THEN
          DO jk=1,nlev_snow
            DO jc = i_startidx, i_endidx
              lnd_diag%t_snow_mult  (jc,jk,jb) = lnd_prog%t_snow_mult_t(jc,jk,jb,1)
              lnd_diag%wliq_snow    (jc,jk,jb) = lnd_prog%wliq_snow_t(jc,jk,jb,1)
              lnd_diag%wtot_snow    (jc,jk,jb) = lnd_prog%wtot_snow_t(jc,jk,jb,1)
              lnd_diag%dzh_snow     (jc,jk,jb) = lnd_prog%dzh_snow_t(jc,jk,jb,1)
            ENDDO
          ENDDO
        ENDIF

      ELSE ! aggregate fields over tiles


        ! First initialize fields to zero in order to prepare
        ! subsequent summation over the tiles
        !
        lnd_diag%t_snow   (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%t_s      (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%w_snow   (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%w_i      (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%h_snow   (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%freshsnow(i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%snowfrac (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%runoff_s (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%runoff_g (i_startidx:i_endidx,jb)  = 0._wp
        lnd_diag%rstom    (i_startidx:i_endidx,jb)  = 0._wp

        lnd_diag%t_so    (i_startidx:i_endidx,:,jb) = 0._wp
        lnd_diag%w_so    (i_startidx:i_endidx,:,jb) = 0._wp
        lnd_diag%w_so_ice(i_startidx:i_endidx,:,jb) = 0._wp

        IF (itype_interception == 2) THEN
          lnd_diag%w_p    (i_startidx:i_endidx,jb)  = 0._wp
          lnd_diag%w_s    (i_startidx:i_endidx,jb)  = 0._wp
        ENDIF

        IF (itype_trvg == 3) THEN
          lnd_diag%plantevap(i_startidx:i_endidx,jb) = 0._wp
        ENDIF

        IF (l2lay_rho_snow .OR. lmulti_snow) THEN
          lnd_diag%rho_snow_mult(i_startidx:i_endidx,:,jb) = 0._wp
        ENDIF

        IF (lmulti_snow) THEN
          lnd_diag%t_snow_mult  (i_startidx:i_endidx,:,jb) = 0._wp
          lnd_diag%wliq_snow    (i_startidx:i_endidx,:,jb) = 0._wp
          lnd_diag%wtot_snow    (i_startidx:i_endidx,:,jb) = 0._wp
          lnd_diag%dzh_snow     (i_startidx:i_endidx,:,jb) = 0._wp
        ENDIF


        ! Aggregate fields without water tiles
        !
        DO isubs = 1,ntiles_total
          DO jc = i_startidx, i_endidx
            !
            ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
            ! therefore we multiply by inv_frland_from_tiles
            tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
              &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

            lnd_diag%t_snow(jc,jb)    = lnd_diag%t_snow(jc,jb) + tilefrac    &
              &                         * lnd_prog%t_snow_t(jc,jb,isubs)
            lnd_diag%w_snow(jc,jb)    = lnd_diag%w_snow(jc,jb) + tilefrac    &
              &                         * lnd_prog%w_snow_t(jc,jb,isubs)
            lnd_diag%w_i(jc,jb)       = lnd_diag%w_i(jc,jb) + tilefrac       &
              &                         * lnd_prog%w_i_t(jc,jb,isubs)
            lnd_diag%h_snow(jc,jb)    = lnd_diag%h_snow(jc,jb) + tilefrac    &
              &                         * lnd_diag%h_snow_t(jc,jb,isubs)
            lnd_diag%freshsnow(jc,jb) = lnd_diag%freshsnow(jc,jb) + tilefrac &
              &                         * lnd_diag%freshsnow_t(jc,jb,isubs)
            lnd_diag%snowfrac(jc,jb)  = lnd_diag%snowfrac(jc,jb) + tilefrac  &
              &                         * lnd_diag%snowfrac_t(jc,jb,isubs)
            lnd_diag%runoff_s(jc,jb)  = lnd_diag%runoff_s(jc,jb) + tilefrac  &
              &                         * lnd_diag%runoff_s_t(jc,jb,isubs)
            lnd_diag%runoff_g(jc,jb)  = lnd_diag%runoff_g(jc,jb) + tilefrac  &
              &                         * lnd_diag%runoff_g_t(jc,jb,isubs)
            lnd_diag%rstom(jc,jb)     = lnd_diag%rstom(jc,jb) + tilefrac     &
              &                         * lnd_diag%rstom_t(jc,jb,isubs)

            IF(lmulti_snow) THEN
              lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_diag%t_snow_mult(jc,nlev_snow+1,jb)+ &
                                        tilefrac * lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,isubs)
            ENDIF
          ENDDO

          IF (itype_interception == 2) THEN
            DO jc = i_startidx, i_endidx
              !
              ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
              ! therefore we multiply by inv_frland_from_tiles
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%w_p    (jc,jb) = lnd_diag%w_p(jc,jb) + tilefrac       &
                &                       * lnd_prog%w_p_t(jc,jb,isubs)
              lnd_diag%w_s    (jc,jb) = lnd_diag%w_s(jc,jb) + tilefrac       &
                &                       * lnd_prog%w_s_t(jc,jb,isubs)
            ENDDO
          ENDIF

          IF (itype_trvg == 3) THEN
            DO jc = i_startidx, i_endidx
              !
              ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
              ! therefore we multiply by inv_frland_from_tiles
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%plantevap(jc,jb) = lnd_diag%plantevap(jc,jb) + tilefrac       &
                &                       * lnd_diag%plantevap_t(jc,jb,isubs)
            ENDDO
          ENDIF

          DO jk=1,nlev_soil
            DO jc = i_startidx, i_endidx
              !
              ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
              ! therefore we multiply by inv_frland_from_tiles
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)         &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%t_so(jc,jk+1,jb)    = lnd_diag%t_so(jc,jk+1,jb) + tilefrac &
                                             * lnd_prog%t_so_t(jc,jk+1,jb,isubs)
              lnd_diag%w_so(jc,jk,jb)      = lnd_diag%w_so(jc,jk,jb) + tilefrac  &
                                             * lnd_prog%w_so_t(jc,jk,jb,isubs)
              lnd_diag%w_so_ice(jc,jk,jb)  = lnd_diag%w_so_ice(jc,jk,jb) + tilefrac  &
                                             * lnd_prog%w_so_ice_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

          IF (l2lay_rho_snow .OR. lmulti_snow) THEN
            DO jk=1,nlev_snow
              DO jc = i_startidx, i_endidx
                tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                  &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

                lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_diag%rho_snow_mult(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%rho_snow_mult_t(jc,jk,jb,isubs)
              ENDDO
            ENDDO
          ENDIF

          IF (lmulti_snow) THEN
            DO jk=1,nlev_snow
              DO jc = i_startidx, i_endidx
                !
                ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
                ! therefore we multiply by inv_frland_from_tiles
                tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                  &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

                lnd_diag%t_snow_mult(jc,jk,jb)   = lnd_diag%t_snow_mult(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%t_snow_mult_t(jc,jk,jb,isubs)
                lnd_diag%wliq_snow(jc,jk,jb)     = lnd_diag%wliq_snow(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%wliq_snow_t(jc,jk,jb,isubs)
                lnd_diag%wtot_snow(jc,jk,jb)     = lnd_diag%wtot_snow(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%wtot_snow_t(jc,jk,jb,isubs)
                lnd_diag%dzh_snow(jc,jk,jb)      = lnd_diag%dzh_snow(jc,jk,jb) + tilefrac   &
                                                   * lnd_prog%dzh_snow_t(jc,jk,jb,isubs)
              ENDDO
            ENDDO
          ENDIF

        ENDDO  ! isubs


        ! Aggregate fields with water tiles
        DO isubs = 1,ntiles_total + ntiles_water
          DO jc = i_startidx, i_endidx
            tilefrac = ext_data%atm%frac_t(jc,jb,isubs)

            lnd_diag%t_s(jc,jb)       = lnd_diag%t_s(jc,jb) + tilefrac       &
              &                         * lnd_prog%t_s_t(jc,jb,isubs)
          ENDDO
        ENDDO

        ! diagnose rho_snow from aggregated values of w_snow and h_snow; 
        ! by convention, snow density is zero in the absence of snow
        DO jc = i_startidx, i_endidx
          lnd_diag%rho_snow(jc,jb) = rhoh2o*(lnd_diag%w_snow(jc,jb)/MAX(dbl_eps,lnd_diag%h_snow(jc,jb)))
        ENDDO

      ENDIF  ! ntiles_total == 1

      ! copy aggregated snowfrac field to snowfrac_lc
      lnd_diag%snowfrac_lc(i_startidx:i_endidx,jb) = lnd_diag%snowfrac(i_startidx:i_endidx,jb)

      ! fill t_so(1) with t_s
      lnd_diag%t_so(i_startidx:i_endidx,1,jb) = lnd_diag%t_s(i_startidx:i_endidx,jb)

      ! Fill t_so(2:nlev_soil+1) with SST for non-land points (fr_land <= frlnd_thrhld) 
      ! Note that all points with fr_land > frlnd_thrhld are stored in the 
      ! land point index list idx_lst_lp
      !
      ! create mask array
      lmask(i_startidx:i_endidx) = .FALSE.
      !
      icount = ext_data%atm%lp_count(jb)
      DO ic = 1, icount
        jc = ext_data%atm%idx_lst_lp(ic,jb)
        lmask(jc) = .TRUE.
      ENDDO  ! ic
      !
      ! fill non-land points with SST
      ! It is assured that all mixed land/water points only contain pure land temperatures
      ! Note that the climatological layer is explicitly omitted.
      DO jk = 2, nlev_soil
        DO jc = i_startidx, i_endidx
          lnd_diag%t_so(jc,jk,jb) = MERGE(lnd_diag%t_so(jc,jk,jb),lnd_diag%t_s(jc,jb),lmask(jc))
        ENDDO  ! jc
      ENDDO  ! jk
      !
      ! climatological layer: Fill with T_CL over non-land points. 
      ! Land points are already filled with T_CL
      DO jc = i_startidx, i_endidx
        lnd_diag%t_so(jc,nlev_soil+1,jb) = MERGE(lnd_diag%t_so(jc,nlev_soil+1,jb), &
          &                                      ext_data%atm%t_cl(jc,jb),         &
          &                                      lmask(jc))
      ENDDO  ! jc


      ! Make sure that aggregated w_so is always larger than air dryness point 
      ! at points where the soiltype allows infiltration of water.
      DO jk=1,nlev_soil
        DO jc = i_startidx, i_endidx
          styp = ext_data%atm%soiltyp(jc,jb)
          IF ( (styp>=3) .AND. (styp<=8)) THEN   ! 3:sand; 8:peat
            lnd_diag%w_so(jc,jk,jb) = MAX(lnd_diag%w_so(jc,jk,jb),dzsoil(jk)*cadp(styp))
          ELSE
            ! make sure that aggregated w_so=0, w_so_ice=0, where no infiltration is allowed
            lnd_diag%w_so(jc,jk,jb)     = 0._wp
            lnd_diag%w_so_ice(jc,jk,jb) = 0._wp
          ENDIF
        ENDDO
      ENDDO

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE aggregate_landvars



  !! Driver routine to (re-)initialize the snowtile index lists in the case of a restart
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD (2012-08-03)
  !!
  SUBROUTINE init_snowtile_lists(p_patch, ext_data, p_lnd_diag)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jb,isubs,i_count,isubs_snow


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,isubs,i_count,isubs_snow), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO isubs = 1, ntiles_lnd
        isubs_snow = isubs + ntiles_lnd

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
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
                                 snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_snowtile_lists



  !-------------------------------------------------------------------------
  !>
  !! Initialize sea-ice and open water index lists for restart and
  !! non-restart runs. Sea-ice and open water points are distinguished
  !! based on fr_seaice which is provided by analysis.
  !!
  !! Note that fr_seaice is potentially modified.
  !! For fr_seaice in ]0,frsi_min[, it is set to 0
  !! For fr_seaice in ]1-frsi_min,1[, it is set to 1.  
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD (2012-08-03)
  !!
  SUBROUTINE init_sea_lists(p_patch, ext_data, p_lnd_diag, lseaice)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch        !< grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag     !< diag vars for sfc
    LOGICAL              , INTENT(IN)    :: lseaice        !< seaice model on/off

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:
    !
    INTEGER :: jb, ic, jc
    INTEGER :: jg
    INTEGER :: i_count_sea, i_count_ice, i_count_water
    INTEGER :: npoints_ice, npoints_wtr
    REAL(wp):: frac_sea                  ! for sanity check

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_utils:init_sea_lists'
!-------------------------------------------------------------------------


    ! patch ID
    jg = p_patch%id

    i_nchdom = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    IF (lseaice) THEN

    ! generate sea-ice and open-water index list
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,i_count_sea,i_count_ice,i_count_water,frac_sea), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk


      ! This (re)initialization is needed for correct functionality of the IAU iteration
      ext_data%atm%spw_count(jb) = 0
      ext_data%atm%spi_count(jb) = 0

      ! Init sub-index lists for sea points. We distinguish between open-water
      ! (i.e. ice free) points and sea-ice points. diag_lnd%fr_seaice is used
      ! to distinguish between open-water and sea-ice.
      ! These index lists are needed by the sea-ice model
      !
      i_count_sea   = ext_data%atm%sp_count(jb)
      i_count_ice   = 0
      i_count_water = 0

      IF (i_count_sea == 0) CYCLE ! skip loop if the index list for the given block is empty

      ! For fr_seaice in ]0,frsi_min[, set fr_seaice to 0
      ! For fr_seaice in ]1-frsi_min,1[, set fr_seaice to 1. 
      ! This will ensure that sea-ice and water fractions sum up exactly 
      ! to the total sea fraction.
      DO ic = 1, i_count_sea
        jc = ext_data%atm%idx_lst_sp(ic,jb)
        IF (p_lnd_diag%fr_seaice(jc,jb) < frsi_min ) THEN
           p_lnd_diag%fr_seaice(jc,jb) = 0._wp
        ENDIF
        IF (p_lnd_diag%fr_seaice(jc,jb) >= (1._wp-frsi_min) ) THEN
           p_lnd_diag%fr_seaice(jc,jb) = 1._wp
        ENDIF
      ENDDO  ! ic



      IF ( ntiles_total == 1 ) THEN  ! no tile approach

        !
        ! mixed water/ice points are not allowed. A sea point can be either
        ! ice-free or completely ice-covered
        !

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count_sea

          jc = ext_data%atm%idx_lst_sp(ic,jb)


          ! ext_data%atm%lc_frac_t(jc,jb,1) already set to 1 in
          ! init_index_lists for sea points


          IF ( p_lnd_diag%fr_seaice(jc,jb) >= 0.5_wp ) THEN
            !
            ! sea-ice point
            !
            i_count_ice = i_count_ice + 1
            ext_data%atm%idx_lst_spi(i_count_ice,jb) = jc
            ext_data%atm%spi_count(jb)               = i_count_ice
            ! set surface area index (needed by turbtran)
            ext_data%atm%sai_t    (jc,jb,isub_seaice) = c_sea
          ELSE
            !
            ! water point: all sea points with fr_seaice < 0.5
            !
            i_count_water = i_count_water + 1
            ext_data%atm%idx_lst_spw(i_count_water,jb) = jc
            ext_data%atm%spw_count(jb)                 = i_count_water
          ENDIF


          ! Initialize frac_t for sea-ice and water points
          ext_data%atm%frac_t(jc,jb,1) = ext_data%atm%lc_frac_t(jc,jb,1)

        ENDDO  ! ic

      ELSE   ! tile approach


!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count_sea

          jc = ext_data%atm%idx_lst_sp(ic,jb)


          ! set sea-ice area fraction (static)
          ! simply copy from water tile. time-dependent fraction will be set lateron
          ext_data%atm%lc_frac_t(jc,jb,isub_seaice)= ext_data%atm%lc_frac_t(jc,jb,isub_water)

          !
          ! seaice point
          !
          IF ( p_lnd_diag%fr_seaice(jc,jb) >= frsi_min ) THEN
            i_count_ice = i_count_ice + 1
            ext_data%atm%idx_lst_spi(i_count_ice,jb) = jc
            ext_data%atm%spi_count(jb)               = i_count_ice

            ! Initialize frac_t for seaice
            ext_data%atm%frac_t(jc,jb,isub_seaice) = ext_data%atm%lc_frac_t(jc,jb,isub_seaice) &
              &                                    * p_lnd_diag%fr_seaice(jc,jb)

!DR Note that sai at seaice points is initialized with c/=c_sea, a corresponding update
!DR of sai_t needs to be added to the procedure which updates the seaice index list.
            ! set surface area index (needed by turbtran)
            ext_data%atm%sai_t    (jc,jb,isub_seaice)  = c_sea
          ELSE
            ext_data%atm%frac_t(jc,jb,isub_seaice) = 0._wp
          ENDIF

          !
          ! water point: all sea points with fr_seaice < (1-frsi_min)
          !
          IF ( p_lnd_diag%fr_seaice(jc,jb) <= (1._wp-frsi_min) ) THEN
            i_count_water = i_count_water + 1
            ext_data%atm%idx_lst_spw(i_count_water,jb) = jc
            ext_data%atm%spw_count(jb)                 = i_count_water

            ! Update frac_t for water tile
            ext_data%atm%frac_t(jc,jb,isub_water)  = ext_data%atm%lc_frac_t(jc,jb,isub_water)  &
              &                                    * (1._wp - p_lnd_diag%fr_seaice(jc,jb))
          ELSE
            ! necessary, since frac_t(jc,jb,isub_water) has already been initialized
            ! with nonzero values in init_index_lists
            ext_data%atm%frac_t(jc,jb,isub_water)  = 0._wp
          ENDIF

        ENDDO  ! ic


        ! Sanity check
        ! Check whether fractions of seaice and non-seaice covered tiles sum up to total sea fraction. 
        DO ic = 1, i_count_sea
          jc = ext_data%atm%idx_lst_sp(ic,jb)
          frac_sea = ext_data%atm%frac_t(jc,jb,isub_water) + ext_data%atm%frac_t(jc,jb,isub_seaice)
          IF ( ABS(frac_sea - ext_data%atm%lc_frac_t(jc,jb,isub_water)) > dbl_eps ) THEN
            WRITE(message_text,'(a,f12.9)') 'frac_seaice + frac_water: ', frac_sea
            CALL message('', TRIM(message_text))
            WRITE(message_text,'(a,f12.9)') 'tot frac_sea: ',  ext_data%atm%lc_frac_t(jc,jb,isub_water)
            CALL message('', TRIM(message_text))
            CALL finish(routine, 'sea-ice + water fractions do not sum up to total sea fraction')
          END IF
        ENDDO  ! jc

      ENDIF   ! IF (ntiles_total == 1)

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


      ! Some diagnostics: compute total number of sea-ice and open water points
      npoints_ice = SUM(ext_data%atm%spi_count(i_startblk:i_endblk))
      npoints_ice = global_sum_array(npoints_ice)
      npoints_wtr = SUM(ext_data%atm%spw_count(i_startblk:i_endblk))
      npoints_wtr = global_sum_array(npoints_wtr)
      WRITE(message_text,'(a,i3,a,i10)') 'Number of sea-ice points in domain',jg, &
        &  ':',npoints_ice
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,a,i10)') 'Number of water points in domain',jg, &
        &  ':',npoints_wtr
      CALL message('', TRIM(message_text))


    ELSE   ! seaice model switched off

      ! copy sea points index list to open-water index list
      !
      ext_data%atm%spw_count(i_startblk:i_endblk)     =     &
        &                     ext_data%atm%sp_count(i_startblk:i_endblk)
      ext_data%atm%idx_lst_spw(:,i_startblk:i_endblk) =     &
        &                     ext_data%atm%idx_lst_sp(:,i_startblk:i_endblk)

      ext_data%atm%spi_count(i_startblk:i_endblk) = 0

    ENDIF  ! lseaice


  END SUBROUTINE init_sea_lists

  !-------------------------------------------------------------------------

  SUBROUTINE diag_snowfrac_tg(istart, iend, lc_class, i_lc_urban, t_snow, t_soiltop, w_snow, &
    & rho_snow, freshsnow, sso_sigma, z0, snowfrac, t_g, meltrate, snowfrac_u)

    INTEGER, INTENT (IN) :: istart, iend ! start and end-indices of the computation

    INTEGER, INTENT (IN) :: lc_class(:)  ! list of land-cover classes
    INTEGER, INTENT (IN) :: i_lc_urban   ! land-cover class index for urban / artificial surface
    REAL(wp), DIMENSION(:), INTENT(IN) :: t_snow, t_soiltop, w_snow, rho_snow, &
      freshsnow, sso_sigma, z0
    REAL(wp), DIMENSION(:), INTENT(IN), OPTIONAL :: meltrate ! snow melting rate in kg/(m**2*s)

    REAL(wp), DIMENSION(:), INTENT(INOUT) :: snowfrac, t_g
    REAL(wp), DIMENSION(:), INTENT(INOUT), OPTIONAL :: snowfrac_u ! unmodified snow-cover fraction

    INTEGER  :: ic
    REAL(wp) :: h_snow, snowdepth_fac, sso_fac, lc_fac, lc_limit


    SELECT CASE (idiag_snowfrac)
    CASE (1) ! old parameterization depending on SWE only
      DO ic = istart, iend
        snowfrac(ic) = MIN(1.0_wp, w_snow(ic)/cf_snow)
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    CASE (2, 20) ! more advanced parameterization depending on snow depth, accounts also for vegetation and SSO
      DO ic = istart, iend
        IF (w_snow(ic) <= 1.e-6_wp) THEN
          snowfrac(ic) = 0._wp
        ELSE
          h_snow = 1000._wp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
          sso_fac = SQRT(0.025_wp*MAX(25._wp,sso_sigma(ic)*(1._wp-freshsnow(ic))))
          snowdepth_fac = h_snow*(17.5_wp*freshsnow(ic)+5._wp+5._wp/sso_fac*(1._wp-freshsnow(ic)))
          lc_fac   = MAX(1._wp,SQRT(10._wp*z0(ic)))
          IF (lc_class(ic) == i_lc_urban) THEN
            lc_limit = 0.875_wp ! this accounts for the effect of human activities on snow cover
          ELSE
            lc_limit = 1._wp
          ENDIF
          snowfrac(ic) = MAX(0.05_wp,MIN(lc_limit,snowdepth_fac/lc_fac))
        ENDIF
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    CASE (3, 30)  ! similar to option 2, but somewhat less snow cover and limit over high vegetation
      DO ic = istart, iend
        IF (w_snow(ic) <= 1.e-6_wp) THEN
          snowfrac(ic) = 0._wp
        ELSE
          h_snow = 1000._wp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
          sso_fac = SQRT(0.025_wp*MAX(25._wp,sso_sigma(ic)*(1._wp-freshsnow(ic))))
          snowdepth_fac = h_snow*(17.5_wp*freshsnow(ic)+5._wp+5._wp/sso_fac*(1._wp-freshsnow(ic)))
          lc_fac   = MAX(1._wp,SQRT(15.0_wp*z0(ic)))
          IF (lc_class(ic) == i_lc_urban) THEN
            lc_limit = 0.8_wp ! this accounts for the effect of human activities on snow cover
          ELSE
            lc_limit = MAX(0.925_wp,MIN(1._wp,1._wp/MAX(0.1_wp,7.5_wp*z0(ic))**0.125_wp))
          ENDIF
          snowfrac(ic) = MIN(lc_limit,snowdepth_fac/lc_fac)
        ENDIF
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    CASE (4, 40)  ! same as option 3, but even more restrictive snow cover limit over high vegetation
      DO ic = istart, iend
        IF (w_snow(ic) <= 1.e-6_wp) THEN
          snowfrac(ic) = 0._wp
        ELSE
          h_snow = 1000._wp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
          sso_fac = SQRT(0.025_wp*MAX(25._wp,sso_sigma(ic)*(1._wp-freshsnow(ic))))
          snowdepth_fac = h_snow*(17.5_wp*freshsnow(ic)+5._wp+5._wp/sso_fac*(1._wp-freshsnow(ic)))
          lc_fac   = MAX(1._wp,SQRT(15.0_wp*z0(ic)))
          IF (lc_class(ic) == i_lc_urban) THEN
            lc_limit = 0.8_wp ! this accounts for the effect of human activities on snow cover
          ELSE
            lc_limit = MAX(0.85_wp,MIN(1._wp,1._wp/MAX(0.1_wp,7.5_wp*z0(ic))**0.125_wp))
          ENDIF
          snowfrac(ic) = MIN(lc_limit,snowdepth_fac/lc_fac)
        ENDIF
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    END SELECT

    ! For the single-layer scheme, t_soiltop represents the weighted average between the snow-covered
    ! part and the snow-free part in the case of melting snow concurrent with above-freezing soil temperatures
    IF (.NOT. lmulti_snow) THEN
      DO ic = istart, iend
        IF (t_soiltop(ic) > tmelt .AND. t_snow(ic) >= tmelt - dbl_eps) t_g(ic) = t_soiltop(ic)
      ENDDO
    ENDIF

    IF (PRESENT(snowfrac_u)) snowfrac_u(istart:iend) = snowfrac(istart:iend) ! save unmodified snow-cover fraction

    SELECT CASE (idiag_snowfrac)
    CASE (20, 30, 40)
      ! Artificially reduce snow-cover fraction in case of melting snow in order to reduce the ubiquitous
      ! cold bias in such situations
      IF (PRESENT(meltrate)) THEN
        DO ic = istart, iend
          snowfrac(ic) = MIN(snowfrac(ic),MAX(tune_minsnowfrac,1._wp-30000._wp*MAX(0._wp,z0(ic)-0.02_wp)*meltrate(ic)))
        ENDDO
      ENDIF
    END SELECT

  END SUBROUTINE diag_snowfrac_tg




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Updating of dynamic index lists
  !!
  !! Routine updates the following dynamic index lists (if required):
  !!
  !! - dynamic
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD (2012-07-01)
  !!
  SUBROUTINE update_idx_lists_lnd (idx_lst_lp, lp_count, idx_lst, gp_count, idx_lst_snow, &
    &             gp_count_snow, lc_frac, partial_frac, partial_frac_snow, snowtile_flag, &
    &             snowtile_flag_snow, snowfrac)


    INTEGER ,    INTENT(   IN) ::  &   !< static list of all land points of a tile index
      &  idx_lst_lp(:), lp_count       !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic list of all snow-free or mixed land points
      &  idx_lst(:), gp_count          !< of a tile index
                                       !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic list of all snow-covered points for land-cover
      &  idx_lst_snow(:), gp_count_snow !< classes eligible for separate treatment
                                       !< and corresponding grid point counts

    INTEGER ,    INTENT(INOUT) ::  &   !< snowtile flag field for snow-free or mixed points
      &  snowtile_flag(:)              !< -1: no separation between snow tile and snow-free tile
                                       !< 0: inactive
                                       !< 1: active
                                       !< 2: newly activated; initialization from corresponding tile required

    INTEGER ,    INTENT(INOUT) ::  &   !< snowtile flag field for snow-covered points
      &  snowtile_flag_snow(:)         !< -1: no separation between snow tile and snow-free tile
                                       !< 0: inactive
                                       !< 1: active
                                       !< 2: newly activated; initialization from corresponding tile required

    REAL(wp),    INTENT(   IN) ::  &   !< area fraction of land-cover class
      &  lc_frac(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< snow-cover fraction
      &  snowfrac(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< snow-free and snow-covered sub-tile area fraction
      &  partial_frac(:), partial_frac_snow(:)


    ! Local variables
    INTEGER  :: ic, jc, icount, icount_snow
    REAL(wp) :: eps = 1.e-6_wp

    !-------------------------------------------------------------------------

    icount = 0
    icount_snow = 0

!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, lp_count
      jc = idx_lst_lp(ic)

      IF (snowtile_flag(jc) == -1) THEN
        icount = icount + 1
        idx_lst(icount) = jc
        partial_frac(jc) = lc_frac(jc)
      ELSE
        ! Reset snowfrac to 0/1 in case of very small deviations (just to be safe)
        IF (snowfrac(jc) < eps) snowfrac(jc) = 0._wp
        IF (1._wp - snowfrac(jc) < eps) snowfrac(jc) = 1._wp
      ENDIF
    ENDDO

!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, lp_count
      jc = idx_lst_lp(ic)

      IF (snowtile_flag(jc) /= -1) THEN
        IF (snowfrac(jc) > 0._wp) THEN ! snow tile is active
          icount_snow = icount_snow + 1
          idx_lst_snow(icount_snow) = jc
          partial_frac_snow(jc) = lc_frac(jc)*snowfrac(jc)
          IF (snowtile_flag_snow(jc) == 0) THEN
            snowtile_flag_snow(jc) = 2 ! newly activated, initialization needed
          ELSE
            snowtile_flag_snow(jc) = 1
          ENDIF
        ELSE
          snowtile_flag_snow(jc) = 0
          partial_frac_snow(jc) = 0._wp
        ENDIF
        IF (snowfrac(jc) < 1._wp) THEN ! snow-free tile is active
          icount = icount + 1
          idx_lst(icount) = jc
          partial_frac(jc) = lc_frac(jc)*(1._wp-snowfrac(jc))
          IF (snowtile_flag(jc) == 0) THEN
            snowtile_flag(jc) = 2 ! newly activated, initialization needed
          ELSE
            snowtile_flag(jc) = 1
          ENDIF
        ELSE
          snowtile_flag(jc) = 0
          partial_frac(jc) = 0._wp
        ENDIF
      ENDIF
    ENDDO

    gp_count = icount
    gp_count_snow = icount_snow

  END SUBROUTINE update_idx_lists_lnd




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Updating of dynamic index lists for sea points
  !!
  !! Routine updates the following dynamic index lists (if required):
  !!
  !! - dynamic open-water and sea-ice point index list, which are sub-index lists
  !!   of the static sea point index list idx_lst_sp
  !!
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2012-08-31)
  !!
  SUBROUTINE update_idx_lists_sea (hice_n, pres_sfc, idx_lst_spw, spw_count,    &
    &                              idx_lst_spi, spi_count, frac_t_ice,          &
    &                              frac_t_water, lc_frac_t_water, fr_seaice,    &
    &                              hice_old, tice_old, albsi_now, albsi_new,    &
    &                              t_g_t_now, t_g_t_new, t_s_t_now, t_s_t_new,  &
    &                              qv_s_t, t_seasfc )


    REAL(wp),    INTENT(IN)    ::  &   !< sea-ice depth at new time level  [m]
      &  hice_n(:)                     !< dim: (nproma)

    REAL(wp),    INTENT(IN)    ::  &   !< surface pressure                 [Pa]
      &  pres_sfc(:)

    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic open water point index list
      &  idx_lst_spw(:), spw_count     !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic sea-ice point index list
      &  idx_lst_spi(:), spi_count     !< and corresponding grid point counts

    REAL(wp),    INTENT(INOUT) ::  &   !< ice-covered and ice-free ocean area fraction
      &  frac_t_ice(:), frac_t_water(:)

    REAL(wp), INTENT(IN) ::        &   !< total ocean area fraction
      &  lc_frac_t_water(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea-ice fraction
      &  fr_seaice(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea-ice depth at old time level  [m]
      &  hice_old(:)                   !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea-ice temperature at old time level  [K]
      &  tice_old(:)                   !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< prognostic sea ice albedo at old time level  [1]
      &  albsi_now(:)                  !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< prognostic sea ice albedo at new time level  [1]
      &  albsi_new(:)                  !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< temperature of water tile (now)  [K]
      &  t_g_t_now(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< temperature of water tile (new)  [K]
      &  t_g_t_new(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< surface temperature of water tile (now)  [K]
      &  t_s_t_now(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< surface temperature of water tile (new)  [K]
      &  t_s_t_new(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< surface specific humidity        [kg/kg]
      &  qv_s_t(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea surface temperature          [kg/kg]
      &  t_seasfc(:)

    ! Local variables
    INTEGER, DIMENSION(SIZE(idx_lst_spi,1)) :: &
      &   idx_lst_spi_old       !< sea-ice index list local copy
    INTEGER  :: ic, jc          !< loop indices
    INTEGER  :: spi_count_old   !< current seaice grid point count

    !-------------------------------------------------------------------------


    ! save old sea-ice index list and grid point count
    idx_lst_spi_old(:) = idx_lst_spi(:)
    spi_count_old      = spi_count

    ! re-initialize sea-ice index list and grid point count
    idx_lst_spi(:) = 0
    spi_count      = 0



    !
    ! update index list for sea-ice and open water
    !
    ! The current sea-ice model does not allow for new sea-ice points to be
    ! created during model integration. However, the number of sea-ice points 
    ! is allowed to decrease with time, i.e. sea-ice points may be converted 
    ! into water points, but not vice versa.
    !
    ! Loop over old sea-ice points, only

    IF ( ntiles_total == 1 ) THEN  ! no tile approach

!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, spi_count_old
        jc = idx_lst_spi_old(ic)

        IF ( hice_n(jc) >= hice_min )  THEN ! still sea-ice point
          spi_count = spi_count + 1
          idx_lst_spi(spi_count) = jc
          ! sea-ice fraction remains unchanged
        ELSE                        ! sea-ice point has turned into water point
          ! initialize new water tile
          spw_count = spw_count + 1
          idx_lst_spw(spw_count) = jc
          ! Initialize temperature of water tile with salt water freezing point
          t_g_t_new(jc) = tf_salt ! if the SST analysis contains a meaningful water
                                  ! temperature for this point, one may also take
                                  ! the latter
          t_g_t_now(jc) = tf_salt 

          t_seasfc(jc)  = tf_salt

          ! Initialize surface saturation specific humidity for new water tile
          ! includes reduction of saturation pressure due to salt content
          qv_s_t(jc) = salinity_fac * spec_humi(sat_pres_water(t_g_t_new(jc)),pres_sfc(jc))

          ! since sea-ice melted away, the sea-ice fraction is re-set to 0
          fr_seaice(jc)  = 0._wp
          !
          ! resetting of frac_t_water and frac_t_ice is not possible without 
          ! tile approach, since both point to the same array location (frac_t:,:,1). 
          ! frac_t=1 is required without tile approach

          ! reset sea-ice temperature and depth at old time level in order to prevent
          ! other schemes from using them incorrectly
          tice_old(jc) = tmelt
          hice_old(jc) = 0._wp
          !
          ! Reset prognostic sea ice albedo for consistency
          IF (lprog_albsi) THEN
            albsi_now(jc) = ALB_SI_MISSVAL
            albsi_new(jc) = ALB_SI_MISSVAL
          ENDIF

        ENDIF

      ENDDO  ! ic

    ELSE
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, spi_count_old
        jc = idx_lst_spi_old(ic)

        IF ( hice_n(jc) >= hice_min )  THEN ! still sea-ice point
          spi_count = spi_count + 1
          idx_lst_spi(spi_count) = jc
          ! sea-ice fraction remains unchanged
        ELSE                        ! sea-ice point has turned into water point
          ! Check whether we need to initialize a new water tile, or whether a water tile
          ! already exists for the given point:
          IF ( fr_seaice(jc) > (1._wp-frsi_min) ) THEN
            ! water tile does not exist for given point
            ! add new water tile to water-points index list and initialize
            spw_count = spw_count + 1
            idx_lst_spw(spw_count) = jc

            ! Initialize new water tile
            !
            ! Initialize temperature with salt water freezing point
            t_g_t_new(jc) = tf_salt ! if the SST analysis contains a meaningful water
                                    ! temperature for this point, one may also take
                                    ! the latter
            t_g_t_now(jc) = tf_salt 

            t_s_t_new(jc) = tf_salt ! otherwise aggregated t_so and t_s will be 
                                    ! 0 at these points
            t_s_t_now(jc) = tf_salt

            t_seasfc(jc)  = tf_salt

            !
            ! Initialize surface saturation specific humidity
            ! includes reduction of saturation pressure due to salt content
            qv_s_t(jc) = salinity_fac * spec_humi(sat_pres_water(t_g_t_new(jc)),pres_sfc(jc))
          ENDIF

          ! re-set dynamic fractions of water and sea-ice
          !
          ! new sea area fraction is the sum of the current water and sea-ice area fractions
          ! to ensure bit-reproducibility with restart, we use lc_frac_t instead of frac_t_water+frac_t_ice
          frac_t_water(jc) = lc_frac_t_water(jc)
          ! since sea-ice melted away, the sea-ice fraction is re-set to 0
          fr_seaice(jc)  = 0._wp
          frac_t_ice(jc) = 0._wp
          !
          ! reset sea-ice temperature and depth at old time level in order to prevent
          ! other schemes from using them incorrectly
          tice_old(jc) = tmelt
          hice_old(jc) = 0._wp
          !
          ! Reset prognostic sea ice albedo for consistency
          IF (lprog_albsi) THEN
            albsi_now(jc) = ALB_SI_MISSVAL
            albsi_new(jc) = ALB_SI_MISSVAL
          ENDIF
        ENDIF

      ENDDO  ! ic

    ENDIF  ! IF ( ntiles_total == 1 )

  END SUBROUTINE update_idx_lists_sea



  !>
  !! After updating the SST and sea ice fraction (from external parameter files),
  !! the fields t_g_t, t_s_t, h_ice, t_ice are updated accordingly. In addition, 
  !! the dynamic index lists for sea points and some prognostic variables related
  !! to them are updated.
  !!
  !!
  !! @par Revision History
  !! Initial release by Pilar Ripodas (2012-12)
  !! Modification ba Daniel Reinert, DWD (2016-07-22)
  !! - add new mode by which SST is updated with increments from climatology on a daily basis
  !!
  SUBROUTINE update_sst_and_seaice (p_patch, ext_data, p_lnd_state, p_nh_state, &
    &                       sstice_mode, ref_datetime, target_datetime )


    TYPE(t_patch),           INTENT(IN)    :: p_patch
    TYPE(t_external_data),   INTENT(INOUT) :: ext_data
    TYPE(t_lnd_state),       INTENT(INOUT) :: p_lnd_state
    TYPE(t_nh_state),        INTENT(IN)    :: p_nh_state
    INTEGER,                 INTENT(IN)    :: sstice_mode
    TYPE(datetime),          INTENT(IN)    :: ref_datetime
    TYPE(datetime),          INTENT(IN)    :: target_datetime

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks

    ! Local scalars:
    !
    INTEGER :: jb, ic, jc
    INTEGER :: jg
    INTEGER :: count_sea, count_ice, count_water
    INTEGER :: npoints_ice, npoints_wtr, npoints_sea
    INTEGER :: n_now, n_new
    REAL(wp):: t_water
    REAL(wp):: fracwater_old, fracice_old

    ! climatological sst field for the model initialization day
    REAL(wp), ALLOCATABLE :: sst_cl_ini_day(:,:)
    ! climatological sst field for the current day
    REAL(wp), ALLOCATABLE :: sst_cl_cur_day(:,:)
    ! climatological SST increment
    REAL(wp) :: sst_cl_inc
    REAL(wp) :: new_sst             ! updated SST value
    REAL(wp) :: max_inc, min_inc    ! max/min SST increment on given PE
    REAL(wp), ALLOCATABLE :: sst_inc(:,:)

    INTEGER :: ierr

    CHARACTER(len=MAX_DATETIME_STR_LEN) :: target_datetime_PTString
    CHARACTER(len=MAX_DATETIME_STR_LEN) :: ref_datetime_PTString
    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:update_sst_and_seaice'

!_______________

    jg = p_patch%id

    ! SST and sea ice fraction are read from the analysis. 
    ! The SST is updated by climatological increments on a daily basis. 
    ! The sea ice fraction may only change due to the seaice model.
    SELECT CASE(sstice_mode)
    CASE (SSTICE_ANA_CLINC)

      IF (msg_level >= 13) THEN
        WRITE(message_text,'(a)') 'Update SST with climatological increments'
        CALL message('', TRIM(message_text))
        !
        CALL datetimeToString(target_datetime, target_datetime_PTString, ierr)
        CALL datetimeToString(ref_datetime, ref_datetime_PTString, ierr)
        WRITE(message_text,'(a,i2,a,a)') 'Target Date for DOM ',jg,': ',TRIM(target_datetime_PTString) 
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i2,a,a)') 'Reference Date for DOM ',jg,': ',TRIM(ref_datetime_PTString) 
        CALL message('', TRIM(message_text))
      ENDIF

      n_now = nnow_rcf(jg)
      n_new = nnew_rcf(jg)

      ALLOCATE(sst_cl_ini_day(nproma,p_patch%nblks_c), &
        &      sst_cl_cur_day(nproma,p_patch%nblks_c), &
        &      sst_inc(nproma,p_patch%nblks_c), STAT=ierr)
      IF (ierr /= SUCCESS)  CALL finish (routine, 'Allocation of sst_cl_ini_day, sst_cl_cur_day  failed!')

      ! get climatological sst for initial and current day
      !
      CALL interpol_monthly_mean(p_patch, ref_datetime,         &! in
        &                        ext_data%atm_td%sst_m,         &! in
        &                        sst_cl_ini_day                 )! out

      CALL interpol_monthly_mean(p_patch, target_datetime,      &! in
        &                        ext_data%atm_td%sst_m,         &! in
        &                        sst_cl_cur_day                 )! out

      ! exclude nest boundary and halo points
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,sst_cl_inc,new_sst)
      DO jb=i_startblk, i_endblk

        ! loop over all open water points and add climatological increments
        DO ic = 1, ext_data%atm%spw_count(jb)
          jc = ext_data%atm%idx_lst_spw(ic,jb)

          sst_cl_inc = sst_cl_cur_day(jc,jb) - sst_cl_ini_day(jc,jb)
          ! make sure, that the updated SST does not drop below 
          ! the salt-water freezing point
          new_sst = MAX(tf_salt,p_lnd_state%diag_lnd%t_seasfc(jc,jb) + sst_cl_inc)
          !
          p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_water) = new_sst
          p_lnd_state%prog_lnd(n_now)%t_s_t(jc,jb,isub_water) = new_sst
          p_lnd_state%prog_lnd(n_new)%t_g_t(jc,jb,isub_water) = new_sst
          p_lnd_state%prog_lnd(n_new)%t_s_t(jc,jb,isub_water) = new_sst

          ! includes reduction of saturation pressure due to salt content
          p_lnd_state%diag_lnd%qv_s_t(jc,jb,isub_water) = salinity_fac *                       & 
            &   spec_humi( sat_pres_water(p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_water)),&
            &                                  p_nh_state%diag%pres_sfc(jc,jb) )

        ENDDO  ! ic
      ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

      ! aggregate updated t_g_t and qv_s_t
      CALL aggregate_tg_qvs( p_patch, ext_data, p_lnd_state%prog_lnd(n_now) , &
           &                           p_lnd_state%diag_lnd )

      ! debug output
      IF (msg_level >= 13) THEN
        sst_inc(:,:) = 0._wp
        DO jb=i_startblk, i_endblk
          ! loop over all open water points and add climatological increments
          DO ic = 1, ext_data%atm%spw_count(jb)
            jc = ext_data%atm%idx_lst_spw(ic,jb)
            sst_inc(jc,jb) =  sst_cl_cur_day(jc,jb) - sst_cl_ini_day(jc,jb)
          ENDDO
        ENDDO
        max_inc = MAXVAL(sst_inc(:,:))
        min_inc = MINVAL(sst_inc(:,:))
        !
        WRITE(message_text,'(2(a,i2,a,e12.5))') 'max increment DOM', jg, ': ', global_max(max_inc), &
          &                                     ', min increment DOM', jg, ': ', global_min(min_inc)
        CALL message('', TRIM(message_text))
      ENDIF

      DEALLOCATE(sst_cl_ini_day, sst_cl_cur_day, sst_inc)



    CASE DEFAULT

      n_now = nnow_rcf(jg)
      n_new = nnew_rcf(jg)


      ! exclude nest boundary and halo points
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      WRITE(message_text,'(a,3i10)') 'start end  blocks, number of cells', &
           &  i_startblk, i_endblk,  p_patch%n_patch_cells_g
      CALL message('', TRIM(message_text))


      ! re-initialized to zero
      ext_data%atm%spi_count(i_startblk:i_endblk)=0
      ext_data%atm%spw_count(i_startblk:i_endblk)=0


      IF (lseaice) THEN

        ! generate sea-ice and open-water index list
        !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,count_sea,count_ice,count_water, &
!$OMP            fracwater_old,fracice_old, t_water), SCHEDULE(guided)
        DO jb = i_startblk, i_endblk

          ! Init sub-index lists for sea points. We distinguish between open-water
          ! (i.e. ice free) points and sea-ice points. diag_lnd%fr_seaice is used
          ! to distinguish between open-water and sea-ice.
          ! These index lists are needed by the sea-ice model
          !
          count_sea   = ext_data%atm%sp_count(jb)
          count_ice   = 0
          count_water = 0


          IF (count_sea == 0) CYCLE ! skip loop if the index list for the given block is empty

          IF ( ntiles_total == 1 ) THEN  ! no tile approach

            !
            ! mixed water/ice points are not allowed. A sea point can be either
            ! ice-free or completely ice-covered
            !

!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, count_sea

              jc = ext_data%atm%idx_lst_sp(ic,jb)

              IF ( p_lnd_state%diag_lnd%fr_seaice(jc,jb) >= 0.5_wp ) THEN
                !
                ! sea-ice point
                !
                count_ice = count_ice + 1
                ext_data%atm%idx_lst_spi(count_ice,jb) = jc
                ext_data%atm%spi_count(jb)               = count_ice
                !
                fracwater_old = ext_data%atm%frac_t(jc,jb,isub_water)
                IF (fracwater_old  > 0._wp  ) THEN
                  p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_seaice)= tf_salt
                  p_lnd_state%prog_lnd(n_now)%t_s_t(jc,jb,isub_seaice)= tf_salt
                  p_lnd_state%prog_wtr(n_now)%t_ice(jc,jb) = tf_salt
                  p_lnd_state%prog_wtr(n_now)%h_ice(jc,jb) = hice_ini_min +           &
                       p_lnd_state%diag_lnd%fr_seaice(jc,jb) * (hice_ini_max-hice_ini_min)

                  IF (lprog_albsi) THEN
                    ! initialize prognostic seaice albedo with equilibrium value
                    p_lnd_state%prog_wtr(n_now)%alb_si(jc,jb) = &
                      &                             alb_seaice_equil( p_lnd_state%prog_wtr(n_now)%t_ice(jc,jb) )
                    p_lnd_state%prog_wtr(n_new)%alb_si(jc,jb) = p_lnd_state%prog_wtr(n_now)%alb_si(jc,jb)
                  ENDIF
                ELSE
                  ! before was also ice.
                END IF

                ! set sai_t
                ext_data%atm%sai_t    (jc,jb,isub_seaice)  = c_sea

              ELSE
                !
                ! water point: all sea points with fr_seaice < 0.5
                !
                count_water = count_water + 1
                ext_data%atm%idx_lst_spw(count_water,jb) = jc
                ext_data%atm%spw_count(jb)                 = count_water

                p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                p_lnd_state%prog_lnd(n_now)%t_s_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                p_lnd_state%prog_lnd(n_new)%t_g_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                p_lnd_state%prog_lnd(n_new)%t_s_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                t_water = p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                ! includes reduction of saturation pressure due to salt content
                p_lnd_state%diag_lnd%qv_s_t(jc,jb,isub_water)    = salinity_fac *         &
                     &                             spec_humi( sat_pres_water(t_water ),   &
                     &                             p_nh_state%diag%pres_sfc(jc,jb) )

                ! set sai_t
                ext_data%atm%sai_t    (jc,jb,isub_water)  = c_sea

                IF (lprog_albsi) THEN
                  ! re-initialize prognostic seaice albedo with ALB_SI_MISSVAL
                  p_lnd_state%prog_wtr(n_now)%alb_si(jc,jb) = ALB_SI_MISSVAL
                  p_lnd_state%prog_wtr(n_new)%alb_si(jc,jb) = ALB_SI_MISSVAL
                ENDIF

              ENDIF

            ENDDO  ! ic

          ELSE   ! tile approach


!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, count_sea

              jc = ext_data%atm%idx_lst_sp(ic,jb)

              ! NOTE:
              ! Note that in copy_initicon2prog
              ! - for fr_seaice in ]0,frsi_min[, we have set fr_seaice to 0
              ! - for fr_seaice in ]1-frsi_min,1[, we have set fr_seaice to 1
              ! This ensures that frac_t for sea-ice and open water at a given point
              ! sums up to lc_frac_t(isub_water) exactly.

              !
              ! seaice point
              !
              IF ( p_lnd_state%diag_lnd%fr_seaice(jc,jb) >= frsi_min ) THEN
                count_ice = count_ice + 1
                ext_data%atm%idx_lst_spi(count_ice,jb) = jc
                ext_data%atm%spi_count(jb)               = count_ice

                fracice_old = ext_data%atm%frac_t(jc,jb,isub_seaice)              &
                     & / ext_data%atm%lc_frac_t(jc,jb,isub_seaice)
                IF (fracice_old <frsi_min  ) THEN
                  p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_seaice)= tf_salt
                  p_lnd_state%prog_lnd(n_now)%t_s_t(jc,jb,isub_seaice)= tf_salt


                  p_lnd_state%diag_lnd%qv_s_t(jc,jb,isub_seaice)    =                    &
                       &                                     spec_humi( sat_pres_ice(tf_salt ),  &
                       &                                  p_nh_state%diag%pres_sfc(jc,jb) )
                  p_lnd_state%prog_wtr(n_now)%t_ice(jc,jb) = tf_salt
                  p_lnd_state%prog_wtr(n_now)%h_ice(jc,jb) = hice_ini_min +           &
                    &  p_lnd_state%diag_lnd%fr_seaice(jc,jb) * (hice_ini_max-hice_ini_min)

                  IF (lprog_albsi) THEN
                    ! initialize prognostic seaice albedo with equilibrium value
                    p_lnd_state%prog_wtr(n_now)%alb_si(jc,jb) = &
                      &                             alb_seaice_equil( p_lnd_state%prog_wtr(n_now)%t_ice(jc,jb) )
                    p_lnd_state%prog_wtr(n_new)%alb_si(jc,jb) = p_lnd_state%prog_wtr(n_now)%alb_si(jc,jb)
                  ENDIF

                END IF
                ! set sai_t
                ext_data%atm%sai_t    (jc,jb,isub_seaice)  = c_sea
                ! set new frac_t for isub_seaice
                ext_data%atm%frac_t(jc,jb,isub_seaice) = ext_data%atm%lc_frac_t(jc,jb,isub_seaice) &
                     &                                    * p_lnd_state%diag_lnd%fr_seaice(jc,jb)

              ELSE
                ext_data%atm%frac_t(jc,jb,isub_seaice) = 0._wp
              ENDIF

              !
              ! water point: all sea points with fr_seaice < (1-frsi_min)
              !
              IF ( p_lnd_state%diag_lnd%fr_seaice(jc,jb) <= (1._wp-frsi_min) ) THEN
                count_water = count_water + 1
                ext_data%atm%idx_lst_spw(count_water,jb) = jc
                ext_data%atm%spw_count(jb)                 = count_water

                p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                p_lnd_state%prog_lnd(n_now)%t_s_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                p_lnd_state%prog_lnd(n_new)%t_g_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                p_lnd_state%prog_lnd(n_new)%t_s_t(jc,jb,isub_water)= p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                t_water = p_lnd_state%diag_lnd%t_seasfc(jc,jb)
                ! includes reduction of saturation pressure due to salt content
                p_lnd_state%diag_lnd%qv_s_t(jc,jb,isub_water)    =  salinity_fac *         &
                     &                             spec_humi( sat_pres_water(t_water ),    &
                     &                             p_nh_state%diag%pres_sfc(jc,jb) )

                fracwater_old = ext_data%atm%frac_t(jc,jb,isub_water)              &
                     & / ext_data%atm%lc_frac_t(jc,jb,isub_water)

                IF ( p_lnd_state%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
                  ! pure water grid point, no seaice tile
                  IF (lprog_albsi) THEN
                    ! re-initialize prognostic seaice albedo with ALB_SI_MISSVAL
                    p_lnd_state%prog_wtr(n_now)%alb_si(jc,jb) = ALB_SI_MISSVAL
                    p_lnd_state%prog_wtr(n_new)%alb_si(jc,jb) = ALB_SI_MISSVAL
                  ENDIF
                ENDIF

                ! set sai_t
                ext_data%atm%sai_t    (jc,jb,isub_water)  = c_sea
                ! Update frac_t for water tile
                ext_data%atm%frac_t(jc,jb,isub_water)  = ext_data%atm%lc_frac_t(jc,jb,isub_water)  &
                     &                                    * (1._wp - p_lnd_state%diag_lnd%fr_seaice(jc,jb))

              ELSE
                ext_data%atm%frac_t(jc,jb,isub_water)  = 0._wp
              ENDIF

            ENDDO  ! ic

          ENDIF   ! IF (ntiles_total == 1)

        ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL


        ! Some diagnostics: compute total number of sea-ice and open water points
        npoints_ice = SUM(ext_data%atm%spi_count(i_startblk:i_endblk))
        npoints_ice = global_sum_array(npoints_ice)
        npoints_wtr = SUM(ext_data%atm%spw_count(i_startblk:i_endblk))
        npoints_wtr = global_sum_array(npoints_wtr)
        npoints_sea = SUM(ext_data%atm%sp_count(i_startblk:i_endblk))
        npoints_sea = global_sum_array(npoints_sea)
        WRITE(message_text,'(a,i3,a,i10)') 'Number of sea-ice points in domain',jg, &
             &  ':',npoints_ice
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i3,a,i10)') 'Number of water points in domain',jg, &
             &  ':',npoints_wtr
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i3,a,i10)') 'Number of sea points in domain',jg, &
             &  ':',npoints_sea
        CALL message('', TRIM(message_text))


      ELSE   ! seaice model switched off

        ! copy sea points index list to open-water index list
        !
        ext_data%atm%spw_count(i_startblk:i_endblk)     =     &
             &                     ext_data%atm%sp_count(i_startblk:i_endblk)
        ext_data%atm%idx_lst_spw(:,i_startblk:i_endblk) =     &
             &                     ext_data%atm%idx_lst_sp(:,i_startblk:i_endblk)

        ext_data%atm%spi_count(i_startblk:i_endblk) = 0

        DO jb = i_startblk, i_endblk
          count_sea   = ext_data%atm%sp_count(jb)

          !CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, count_sea

            jc = ext_data%atm%idx_lst_sp(ic,jb)
            ! only if the dominant tile is water t_g_t is set to t_seasfc
            IF (p_lnd_state%diag_lnd%fr_seaice(jc,jb) < 0.5_wp ) THEN
              p_lnd_state%prog_lnd(n_now)%t_g_t(jc,jb,isub_water)=   &
                   p_lnd_state%diag_lnd%t_seasfc(jc,jb)
              p_lnd_state%prog_lnd(n_now)%t_s_t(jc,jb,isub_water)=   &
                   p_lnd_state%diag_lnd%t_seasfc(jc,jb)
              p_lnd_state%prog_lnd(n_new)%t_g_t(jc,jb,isub_water)=   &
                   p_lnd_state%diag_lnd%t_seasfc(jc,jb)
              p_lnd_state%prog_lnd(n_new)%t_s_t(jc,jb,isub_water)=   &
                   p_lnd_state%diag_lnd%t_seasfc(jc,jb)
              ! includes reduction of saturation pressure due to salt content
              p_lnd_state%diag_lnd%qv_s_t(jc,jb,isub_water)    =  salinity_fac *        &
                   &   spec_humi( sat_pres_water(p_lnd_state%diag_lnd%t_seasfc(jc,jb) ),&
                   &                                  p_nh_state%diag%pres_sfc(jc,jb) )
            END IF
          ENDDO  ! ic

        ENDDO  ! jb

      ENDIF  ! lseaice


      ! aggregate updated t_g_t and qv_s_t
      CALL aggregate_tg_qvs( p_patch, ext_data, p_lnd_state%prog_lnd(n_now) , &
           &                           p_lnd_state%diag_lnd )

    END SELECT

  END SUBROUTINE update_sst_and_seaice


  !-------------------------------------------------------------------------
  !!
  !! Update fields that depend on the actual ndvi ratio
  !!
  !! @par Revision History
  !! Initial release by Pilar Ripodas (2013-02)
  !!
  !!
  SUBROUTINE update_ndvi_dependent_fields( p_patch, ext_data, lnd_diag )

    TYPE(t_patch)        , INTENT(IN)    :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_diag),      INTENT(IN)    :: lnd_diag

    ! local
    INTEGER  :: jb,jt,ic,jc, jt_in
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk    !> blocks
    INTEGER  :: i_count
    INTEGER  :: lu_subs

    !-------------------------------------------------------------------------


     ! exclude the boundary interpolation zone of nested domains
     rl_start = grf_bdywidth_c+1
     rl_end   = min_rlcell_int

     i_startblk = p_patch%cells%start_block(rl_start)
     i_endblk   = p_patch%cells%end_block(rl_end)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,ic,i_count,jc,lu_subs,jt_in)
     DO jb = i_startblk, i_endblk
       IF (ext_data%atm%lp_count(jb) == 0) CYCLE ! skip loop if there is no land point
       IF (ntiles_lnd == 1) THEN
         i_count = ext_data%atm%lp_count_t(jb,1)
!CDIR NODEP,VOVERTAKE,VOB
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp_t(ic,jb,1)
           ! plant cover
           ext_data%atm%plcov_t  (jc,jb,1)  = ext_data%atm%ndviratio(jc,jb)  &
             &     * MIN(ext_data%atm%ndvi_max(jc,jb),ext_data%atm%plcov_mx(jc,jb))
           ! total area index
           ext_data%atm%tai_t    (jc,jb,1)  = ext_data%atm%plcov_t  (jc,jb,1)  &
             &                                  * ext_data%atm%lai_mx(jc,jb)
           ! surface area index
           ext_data%atm%sai_t    (jc,jb,1)  = c_lnd+ext_data%atm%tai_t(jc,jb,1)

         END DO
       ELSE ! ntiles_lnd > 1
         DO jt=1,ntiles_lnd
           i_count = ext_data%atm%lp_count_t(jb,jt)
!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)

             ! fix for non-dominant land points
             !!! only active for 'old' extpar datasets (20131009 and earlier) !!!
             IF (ext_data%atm%fr_land(jc,jb) < 0.5_wp) THEN
               IF (ext_data%atm%ndviratio(jc,jb) <= 0.0_wp) THEN  ! here: 0=extpar_missval
                 ! ... reset ndviratio to 0.5
                 ext_data%atm%ndviratio(jc,jb) = 0.5_wp
               ENDIF
             ENDIF

             lu_subs = ext_data%atm%lc_class_t(jc,jb,jt)
             IF (lu_subs < 0) CYCLE

             ! plant cover
             ext_data%atm%plcov_t  (jc,jb,jt)  = ext_data%atm%ndviratio(jc,jb)   &
               & * MIN(ext_data%atm%ndvi_max(jc,jb),ext_data%atm%plcovmax_lcc(lu_subs))
             ! total area index
             ext_data%atm%tai_t    (jc,jb,jt)  = ext_data%atm%plcov_t(jc,jb,jt)  &
               & * ext_data%atm%laimax_lcc(lu_subs)
             ! surface area index
             ext_data%atm%sai_t    (jc,jb,jt)  = c_lnd+ ext_data%atm%tai_t (jc,jb,jt)

           END DO !ic
         END DO !jt

       END IF !ntiles
!

       IF (lsnowtile) THEN ! copy static external data fields to snow tile grid points
         DO jt = ntiles_lnd+1, ntiles_total

           jt_in = jt - ntiles_lnd
!CDIR NODEP,VOVERTAKE,VOB
           DO ic = 1, ext_data%atm%lp_count_t(jb,jt)
             jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)
             ext_data%atm%plcov_t(jc,jb,jt)    = ext_data%atm%plcov_t(jc,jb,jt_in)
             ext_data%atm%tai_t(jc,jb,jt)      = ext_data%atm%tai_t(jc,jb,jt_in)
             ext_data%atm%sai_t(jc,jb,jt)      = ext_data%atm%sai_t(jc,jb,jt_in)
           ENDDO !ic

         ENDDO !jt
       ENDIF !lsnowtile

     ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL


    IF (itype_vegetation_cycle >= 2) THEN
      CALL vege_clim (p_patch, ext_data, lnd_diag)
    ENDIF

    CALL diagnose_ext_aggr (p_patch, ext_data)


  END SUBROUTINE update_ndvi_dependent_fields

!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  !! Copies the tile-based prognostic land-state variables from time level now to 
  !! time level new. This has no relevance for the forecast results but
  !! avoids missing values when writing output at an odd integer multiple of the global
  !! physics time step
  !!
  !! @par Revision History
  !! Initial revision by Guenther Zaengl, DWD (2014-06-02)
  !!
  !-------------------------------------------------------------------------
 
  SUBROUTINE copy_lnd_prog_now2new(p_patch, p_prog_lnd_now, p_prog_lnd_new)


    TYPE(t_patch),         INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_new
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: is, ie    !< slices

    ! Local :

    INTEGER :: jb,jt
!-------------------------------------------------------------------------


    ! include nest boundary and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,is,ie)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & is, ie, rl_start, rl_end)

      DO jt = 1, ntiles_total+ntiles_water
        p_prog_lnd_new%t_s_t(is:ie,jb,jt) = p_prog_lnd_now%t_s_t(is:ie,jb,jt)
      ENDDO

      DO jt = 1, ntiles_total
        p_prog_lnd_new%w_i_t(is:ie,jb,jt)        = p_prog_lnd_now%w_i_t(is:ie,jb,jt)
        p_prog_lnd_new%t_so_t(is:ie,:,jb,jt)     = p_prog_lnd_now%t_so_t(is:ie,:,jb,jt)
        p_prog_lnd_new%w_so_t(is:ie,:,jb,jt)     = p_prog_lnd_now%w_so_t(is:ie,:,jb,jt)
        p_prog_lnd_new%w_so_ice_t(is:ie,:,jb,jt) = p_prog_lnd_now%w_so_ice_t(is:ie,:,jb,jt)
        p_prog_lnd_new%t_snow_t(is:ie,jb,jt)     = p_prog_lnd_now%t_snow_t(is:ie,jb,jt)
        p_prog_lnd_new%w_snow_t(is:ie,jb,jt)     = p_prog_lnd_now%w_snow_t(is:ie,jb,jt)
        p_prog_lnd_new%rho_snow_t(is:ie,jb,jt)   = p_prog_lnd_now%rho_snow_t(is:ie,jb,jt)
      ENDDO

      IF (l2lay_rho_snow .OR. lmulti_snow) THEN
        DO jt = 1, ntiles_total
          p_prog_lnd_new%rho_snow_mult_t(is:ie,:,jb,jt) = p_prog_lnd_now%rho_snow_mult_t(is:ie,:,jb,jt)
        ENDDO
      ENDIF

      IF (lmulti_snow) THEN
        DO jt = 1, ntiles_total
          p_prog_lnd_new%t_snow_mult_t(is:ie,:,jb,jt)   = p_prog_lnd_now%t_snow_mult_t(is:ie,:,jb,jt)
          p_prog_lnd_new%wliq_snow_t(is:ie,:,jb,jt)     = p_prog_lnd_now%wliq_snow_t(is:ie,:,jb,jt)
          p_prog_lnd_new%wtot_snow_t(is:ie,:,jb,jt)     = p_prog_lnd_now%wtot_snow_t(is:ie,:,jb,jt)
          p_prog_lnd_new%dzh_snow_t(is:ie,:,jb,jt)      = p_prog_lnd_now%dzh_snow_t(is:ie,:,jb,jt)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE copy_lnd_prog_now2new


  !>
  !! Coldstart initialization of prognostic seaice albedo alb_si
  !!
  !! Coldstart initialization of prognostic seaice albedo alb_si
  !! 
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2016-08-18)
  !!
  SUBROUTINE seaice_albedo_coldstart(p_patch, p_lnd_state, ext_data)

    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data

    ! local
    INTEGER :: jb, jc, ic          ! loop counter
    INTEGER :: jg
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startblk, i_endblk
    !
    REAL(wp):: zfrice_thrhld
    REAL(wp):: frsi(nproma)
    REAL(wp):: t_ice_now(nproma)
    REAL(wp):: alb_si_now(nproma), alb_si_new(nproma)
    !
    TYPE(t_lnd_diag),  POINTER :: lnd_diag
    TYPE(t_wtr_prog),  POINTER :: wtr_prog_now, wtr_prog_new
  !------------------------------------------------------------

    jg = p_patch%id

    WRITE(message_text,'(a,i3)') 'Coldstart initialization of prognostic seaice albedo in domain ',jg
    CALL message('', TRIM(message_text))

    ! set frice_thrhld depending on tile usage
    IF ( ntiles_total == 1 ) THEN  ! no tile approach
      zfrice_thrhld = 0.5_wp
    ELSE
      zfrice_thrhld = frsi_min
    ENDIF

    lnd_diag     => p_lnd_state%diag_lnd
    wtr_prog_now => p_lnd_state%prog_wtr(nnow_rcf(jg))
    wtr_prog_new => p_lnd_state%prog_wtr(nnew_rcf(jg))

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,ic,frsi,t_ice_now,alb_si_now,alb_si_new)
    DO jb = i_startblk, i_endblk
 
      DO ic = 1, ext_data%atm%sp_count(jb)
        jc = ext_data%atm%idx_lst_sp(ic,jb)
        !
        frsi(ic)      = lnd_diag%fr_seaice(jc,jb)
        t_ice_now(ic) = wtr_prog_now%t_ice(jc,jb)
        ! the following 2 lines are required, because 
        ! only a subset of points in alb_si_now(1:sp_count) 
        ! are filled by the following initialization routine 
        ! (i.e. sea-ice points only).
        alb_si_now(ic)= wtr_prog_now%alb_si(jc,jb)
        alb_si_new(ic)= wtr_prog_new%alb_si(jc,jb)
      ENDDO

      CALL  seaice_coldinit_albsi_nwp (                                 & 
        &                    nswgb        = ext_data%atm%sp_count(jb),  & !in
        &                    frice_thrhld = zfrice_thrhld,              & !in
        &                    frsi         = frsi(:),                    & !in
        &                    tice_p       = t_ice_now(:),               & !in
        &                    albsi_p      = alb_si_now(:),              & !inout
        &                    albsi_n      = alb_si_new(:)               & !inout
        &  )

      DO ic = 1, ext_data%atm%sp_count(jb)
        jc = ext_data%atm%idx_lst_sp(ic,jb)
        wtr_prog_now%alb_si(jc,jb) = alb_si_now(ic)
        wtr_prog_new%alb_si(jc,jb) = alb_si_new(ic)
      ENDDO
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE seaice_albedo_coldstart

END MODULE mo_nwp_sfc_utils

