!>
!!  Initialize the physical schemes at start time
!!
!! @author Kristina  Froehlich, DWD
!!
!! @par Revision History
!! First implementations by Kristina Froehlich, DWD, 2010-070-20
!! Include initialitation of SST for APE experiments by P. Ripodas, DWD,2010-11
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

MODULE mo_nwp_phy_init

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_physical_constants,  ONLY: grav, rd_o_cpd, cpd, p0ref, rd, p0sl_bg,         &
    &                               dtdz_standardatm, lh_v=>alv
!   USE mo_math_utilities,      ONLY: sphere_cell_mean_char_length
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag,t_nwp_phy_tend
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ext_data_state,      ONLY: nlev_o3, nmonths
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_exception,           ONLY: message, finish !, message_text
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell, min_rlcell_int, zml_soil, io3_ape,  &
    &                               MODE_COMBINED, MODE_IFSANA, icosmo, ismag,      &
    &                               igme, iedmf, SUCCESS, MAX_CHAR_LENGTH,          &
    &                               MODE_COSMO, iss, iorg, ibc, iso4, idu
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltestcase, iqv, iqc, iqr, iqi, iqs, iqg, iqnc,  &
    &                               iqnr, iqni, iqns, iqng, inccn, ininpot, msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, lrtm_filename,              &
    &                               cldopt_filename, icpl_aero_conv, iprog_aero
  !radiation
  USE mo_newcld_optics,       ONLY: setup_newcld_optics
  USE mo_lrtm_setup,          ONLY: lrtm_setup
  USE mo_radiation_config,    ONLY: ssi_radt, tsi_radt,irad_o3, irad_aero, rad_csalbw
  USE mo_srtm_config,         ONLY: setup_srtm, ssi_amip
  USE mo_radiation_rg_par,    ONLY: rad_aibi
  USE mo_aerosol_util,        ONLY: init_aerosol_dstrb_tanre,                       &
    &                               init_aerosol_props_tanre_rg,                    &
    &                               init_aerosol_props_tanre_rrtm,                  &
    &                               init_aerosol_props_tegen_rg,                    &
    &                               init_aerosol_props_tegen_rrtm,                  &
    &                               zaef_rg, zaea_rg, zaes_rg, zaeg_rg,             &
    &                               zaea_rrtm, zaes_rrtm, zaeg_rrtm
  USE mo_o3_util,             ONLY: o3_pl2ml!, o3_zl2ml
  USE mo_psrad_lrtm_setup,    ONLY: setup_lrtm
  USE mo_psrad_srtm_setup,    ONLY: setup_srtm_psrad => setup_srtm
  USE mo_psrad_cloud_optics,  ONLY: setup_cloud_optics  
  USE mo_psrad_interface,     ONLY: setup_psrad

  ! microphysics
  USE gscp_data,              ONLY: gscp_set_coefficients
  USE mo_mcrph_sb,            ONLY: two_moment_mcrph_init,       &
    &                               set_qnc, set_qnr, set_qni,   &
    &                               set_qns, set_qng
  USE mo_art_clouds_interface,ONLY: art_clouds_interface_2mom_init
  USE mo_cpl_aerosol_microphys, ONLY: lookupcreate_segalkhain, specccn_segalkhain_simple, &
                                      ncn_from_tau_aerosol_speccnconst

  ! convection
  USE mo_cuparameters,        ONLY: sucst,  sucumf,    &
    &                               su_yoethf,         &
    &                               sucldp, suphli,    &
    &                               suvdf , suvdfs
  ! EDMF DUAL turbulence
  USE mo_edmf_param,          ONLY: suct0, su0phy, susekf, susveg, sussoil
  ! turbulence
  USE mo_turbdiff_config,     ONLY: turbdiff_config
  USE mo_data_turbdiff,       ONLY: get_turbdiff_param, lsflcnd, &
                                    impl_s, impl_t,              &
                                    impl_weight
  USE src_turbdiff,           ONLY: organize_turbdiff

  USE mo_nwp_sfc_utils,       ONLY: nwp_surface_init, init_snowtile_lists, init_sea_lists, &
    &                               aggregate_tg_qvs, copy_lnd_prog_now2new
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, lsnowtile, ntiles_water, &
    &                               lseaice
  USE mo_phyparam_soil,       ONLY: csalbw!, z0_lu
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                                sat_pres_ice, &  !! saturation vapor pressure w.r.t. ice
    &                                spec_humi !,qsat_rho !! Specific humidity

  USE data_gwd,               ONLY: sugwwms

  USE mo_nh_testcases_nml,    ONLY: nh_test_name, ape_sst_case, th_cbl, sol_const
  USE mo_ape_params,          ONLY: ape_sst
  USE mo_master_config,       ONLY: isRestart
  USE mo_nwp_parameters,      ONLY: t_phy_params

  USE mo_initicon_config,     ONLY: init_mode

  USE mo_nwp_ww,              ONLY: configure_ww
  USE mo_nwp_tuning_config,   ONLY: tune_gkwake, tune_gkdrag, tune_gfrcrit, tune_grcrit, tune_zceff_min, &
    &                               tune_v0snow, tune_zvz0i
  USE mo_sso_cosmo,           ONLY: sso_cosmo_init_param
  USE mo_cuparameters,        ONLY: sugwd
  USE mo_fortran_tools,       ONLY: init
  USE mtime,                  ONLY: datetime, MAX_DATETIME_STR_LEN, &
    &                               datetimeToString, newDatetime, deallocateDatetime
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,         &
    &                                  calculate_time_interpolation_weights

  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: init_nwp_phy, init_cloud_aero_cpl

CONTAINS


SUBROUTINE init_nwp_phy ( p_patch, p_metrics,                  &
                       &  p_prog_now,  p_diag,                 &
                       &  prm_diag,prm_nwp_tend,               &
                       &  p_prog_lnd_now, p_prog_lnd_new,      &
                       &  p_prog_wtr_now, p_prog_wtr_new,      &
                       &  p_diag_lnd,                          &
                       &  ext_data, phy_params, ini_date, &
                       &  lnest_start)

  TYPE(t_patch),        TARGET,INTENT(in)    :: p_patch
  TYPE(t_nh_metrics),          INTENT(in)    :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout) :: p_prog_now !!the prognostic variables
  TYPE(t_nh_diag),      TARGET,INTENT(inout) :: p_diag  !!the diagostic variables
  TYPE(t_external_data),       INTENT(inout) :: ext_data
  TYPE(t_nwp_phy_diag),        INTENT(inout) :: prm_diag
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout) :: prm_nwp_tend
  TYPE(t_lnd_prog),            INTENT(inout) :: p_prog_lnd_now, p_prog_lnd_new
  TYPE(t_wtr_prog),            INTENT(inout) :: p_prog_wtr_now, p_prog_wtr_new
  TYPE(t_lnd_diag),            INTENT(inout) :: p_diag_lnd
  TYPE(t_phy_params),          INTENT(inout) :: phy_params
  TYPE(datetime),              POINTER       :: ini_date     ! current datetime (mtime)
  LOGICAL, INTENT(IN), OPTIONAL              :: lnest_start

  INTEGER             :: jk, jk1
  REAL(wp)            :: rsltn   ! horizontal resolution
  REAL(wp)            :: pref(p_patch%nlev)
  REAL(wp)            :: zlat, zprat, zn1, zn2, zcdnc
  REAL(wp)            :: zpres, zpres0
  REAL(wp)            :: gz0(nproma), l_hori(nproma)
  REAL(wp)            :: scale_fac ! scale factor used only for RCE cases

  INTEGER             :: icur_date    ! current date converted to integer

  ! Reference atmosphere parameters
  REAL(wp), PARAMETER :: htropo = 11000._wp       ! [m]    tropopause height
  REAL(wp), PARAMETER :: t00    = 288.15_wp       ! [m]    temperature at sea level

  REAL(wp), PARAMETER :: grav_o_rd = grav / rd
  REAL(wp), PARAMETER :: cpd_o_rd  = cpd  / rd

  REAL(wp), PARAMETER :: pr800  = 800._wp / 1013.25_wp
  REAL(wp), PARAMETER :: pr400  = 400._wp / 1013.25_wp

  REAL(wp) :: ttropo, ptropo, temp, zfull

  REAL(wp) :: dz1, dz2, dz3, fact_z0rough
  REAL(wp), ALLOCATABLE :: zrefpres(:,:,:)   ! ref press computed from ref exner
  REAL(wp), ALLOCATABLE :: zreftemp(:,:,:)   ! ref temp computed from ref exner
  REAL(wp), ALLOCATABLE :: zpres_sfc(:,:)    ! ref sfc press
  REAL(wp), ALLOCATABLE :: zpres_ifc(:,:,:)  ! ref press at interfaces

  LOGICAL :: lland, lglac, lshallow, ldetrain_prec
  LOGICAL :: ltkeinp_loc, lgz0inp_loc  !< turbtran switches
  LOGICAL :: linit_mode, lturb_init

  INTEGER :: jb,ic,jc,jt,jg,ist
  INTEGER :: nlev, nlevp1, nlevcm    !< number of full, half and canopy levels
  INTEGER :: nshift                  !< shift with respect to global grid
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !! slices
  INTEGER :: i_nchdom                !! domain index
  INTEGER :: lc_class,i_lc_si

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''

  INTEGER :: nblks_c

  INTEGER :: k1500m                  ! index of first half level above 1500m
  INTEGER :: istatus=0

  REAL(wp) :: hag                    ! height above ground
  REAL(wp) :: h850_standard, h950_standard  ! height of 850hPa and 950hPa level in m

  REAL(wp) :: N_cn0,z0_nccn,z1e_nccn,N_in0,z0_nin,z1e_nin     ! for CCN and IN in case of gscp=5

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
     routine = 'mo_nwp_phy_init:init_nwp_phy'

  CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datetime_string, yyyymmdd

  ! Local control variable for extended turbulence initializations
  IF (ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMO/) == init_mode) ) THEN
    lturb_init = .TRUE.
  ELSE
    lturb_init = .FALSE.
  ENDIF

  ! This is needed for correct flow control when a nested domain is initialized after restarting
  IF (PRESENT(lnest_start)) THEN
    linit_mode = lnest_start
    lturb_init = .TRUE.
  ELSE
    linit_mode = .NOT. isRestart()
  ENDIF

  i_nchdom  = MAX(1,p_patch%n_childdom)

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1
  jg     = p_patch%id

  nshift = p_patch%nshift_total

  i_lc_si= ext_data%atm%i_lc_snow_ice

  nblks_c = p_patch%nblks_c

  dz1 = 0.0_wp
  dz2 = 0.0_wp
  dz3 = 0.0_wp

  IF ( nh_test_name == 'RCE' ) THEN
    ! allocate storage var for press to be used in o3_pl2ml
    ALLOCATE (zrefpres(nproma,nlev,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish (TRIM(routine), &
                 'allocation of zrefpres failed')
    END IF
    ALLOCATE (zreftemp(nproma,nlev,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish (TRIM(routine), &
                 'allocation of zreftemp failed')
    END IF
    ALLOCATE (zpres_sfc(nproma,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish (TRIM(routine), &
                 'allocation of zpres_sfc failed')
    END IF
    ALLOCATE (zpres_ifc(nproma,nlevp1,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish (TRIM(routine), &
                 'allocation of zpres_ifc failed')
    END IF
    zrefpres = 0.0_wp
    zreftemp = 0.0_wp
    zpres_sfc = 0.0_wp
    zpres_ifc = 0.0_wp
  END IF

  ! for both restart and non-restart runs. Could not be included into
  ! mo_ext_data_state/init_index_lists due to its dependence on p_diag_lnd.
  CALL init_sea_lists(p_patch, ext_data, p_diag_lnd, lseaice)


  ! mask field to distinguish between tropics and extratropics (needed for some tuning measures)
  !
  rl_start = 1 ! Initialization should be done for all points
  rl_end   = min_rlcell

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx,i_endidx
      zlat = ABS(p_patch%cells%center(jc,jb)%lat*rad2deg)
      IF (zlat < 25._wp) THEN
        prm_diag%tropics_mask(jc,jb) = 1._wp
      ELSE IF (zlat > 30._wp) THEN
        prm_diag%tropics_mask(jc,jb) = 0._wp
      ELSE
        prm_diag%tropics_mask(jc,jb) = (30._wp-zlat)/5._wp
      ENDIF
      IF (zlat < 12.5_wp) THEN
        prm_diag%innertropics_mask(jc,jb) = 1._wp
      ELSE IF (zlat > 17.5_wp) THEN
        prm_diag%innertropics_mask(jc,jb) = 0._wp
      ELSE
        prm_diag%innertropics_mask(jc,jb) = (17.5_wp-zlat)/5._wp
      ENDIF
    ENDDO
  ENDDO

  IF (linit_mode) THEN

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &  i_startidx, i_endidx, rl_start, rl_end)

      IF (ltestcase .AND. (nh_test_name == 'APE_nwp' .OR. nh_test_name == 'dcmip_tc_52') ) THEN

        ! t_g = ape_sst1

        DO jc = i_startidx, i_endidx
          zlat = p_patch%cells%center(jc,jb)%lat
          p_prog_lnd_now%t_g  (jc,jb)   = ape_sst(ape_sst_case,zlat) ! set SST
          p_prog_lnd_new%t_g  (jc,jb)   = ape_sst(ape_sst_case,zlat)
          p_prog_lnd_now%t_g_t(jc,jb,1) = ape_sst(ape_sst_case,zlat)
          p_prog_lnd_new%t_g_t(jc,jb,1) = ape_sst(ape_sst_case,zlat)
          ! Humidity at water surface = humidity at saturation
          p_diag_lnd%qv_s(jc,jb)     = &
            &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
        END DO


      ELSE IF (ltestcase .AND. nh_test_name == 'wk82' ) THEN

        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)*  &
                    ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
          p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
          p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)
          p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)

          p_diag_lnd%qv_s     (jc,jb) = &
            &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
          p_diag_lnd%qv_s    (jc,jb) = MIN (p_diag_lnd%qv_s(jc,jb) ,   &
            &   p_prog_now%tracer(jc,nlev,jb,iqv))
          p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
        END DO

      ELSE IF (ltestcase .AND. nh_test_name == 'RCE' .AND. atm_phy_nwp_config(jg)%inwp_turb/=ismag) THEN !

        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g  (jc,jb)   = th_cbl(1)
          p_prog_lnd_new%t_g  (jc,jb)   = p_prog_lnd_now%t_g (jc,jb)
          p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)
          p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)
          p_diag_lnd%qv_s(jc,jb)     = &
            &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
          p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
        END DO

      ELSE IF (ltestcase) THEN ! any other testcase

        ! t_g  =  t(nlev)
        ! qv_ s= qv(nlev)
        ! KF increase the surface values to obtain fluxes

        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g  (jc,jb)   = p_diag%temp (jc,nlev,jb)!+0.2_wp
          p_prog_lnd_new%t_g  (jc,jb)   = p_diag%temp (jc,nlev,jb)!+0.2_wp
          p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g  (jc,jb)
          p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g  (jc,jb)
          ! KF NOTE: as long as we have only water as lower boundary
          ! this is the same setting as for APE
          p_diag_lnd%qv_s    (jc,jb) = &
            & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
          p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
        END DO
      ELSE ! For real-case simulations, initialize also qv_s and the tile-based fields

        !t_g_t and qv_s_t are initialized in read_dwdfg_sfc, calculate the aggregated values
        ! needed for example for initializing the turbulence fields
        IF (init_mode /= MODE_IFSANA) THEN
          CALL aggregate_tg_qvs( p_patch, ext_data, p_prog_lnd_now , &
          &                           p_diag_lnd )
          DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g(jc,jb)     =  p_prog_lnd_now%t_g(jc,jb)
          ENDDO

          DO jt = 1, ntiles_total+ntiles_water
            DO jc = i_startidx, i_endidx
              p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g_t(jc,jb,jt)
            END DO
          END DO
        END IF  ! init_mode /= MODE_IFSANA

        ! MODE_IFSANA
        ! t_g:
        ! Note, that in copy_prepicon2prog the entire t_g field is initialized with
        ! t_skin. Lake points are re-initialized with MIN(306.15_wp,tskin).
        !
        ! Here, t_g is re-initialized over sea water points with t_seasfc.
        ! Thus:
        ! t_g = tskin (from IFS), for land, lake and seaice points
        ! t_g = t_seasfc for open water
        !
        ! If l_sst_in==FALSE, then t_seasfc=t_skin (with a limiter), so nothing important happens
        !
        ! qv_s:
        ! Over the sea and over the ice, qv_s is set to the saturated value
        ! Over the land we take the minimum of the saturated value and the value
        ! at the first main level above ground
        !

        ! t_g_t, qv_s and qv_s_t are not initialized in case of MODE_IFSANA
        IF (init_mode == MODE_IFSANA) THEN
          DO ic=1, ext_data%atm%spw_count(jb)
            jc = ext_data%atm%idx_lst_spw(ic,jb)
            IF (lseaice) THEN
              ! all points are open water points
              p_prog_lnd_now%t_g(jc,jb) = p_diag_lnd%t_seasfc(jc,jb)
            ELSE
              ! only points with fr_seaice(jc,jb) <= 0.5_wp are open water points and thus
              ! re-initialized with t_seasfc
              IF (p_diag_lnd%fr_seaice(jc,jb) <= 0.5_wp) THEN   ! water point
                p_prog_lnd_now%t_g(jc,jb) = p_diag_lnd%t_seasfc(jc,jb)
              ENDIF
            ENDIF
            p_diag_lnd%qv_s    (jc,jb)    = &
              & spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

          DO ic=1, ext_data%atm%spi_count(jb)
            jc = ext_data%atm%idx_lst_spi(ic,jb)
            p_diag_lnd%qv_s    (jc,jb)    = &
              & spec_humi(sat_pres_ice(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

          DO ic=1, ext_data%atm%fp_count(jb)
            jc = ext_data%atm%idx_lst_fp(ic,jb)
            ! lake points already initialized in mo_initicon_utils:copy_initicon2prog_sfc
            p_diag_lnd%qv_s    (jc,jb)    = &
              & spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

          DO ic=1, ext_data%atm%lp_count(jb)
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            p_diag_lnd%qv_s(jc,jb) = &
              &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
            p_diag_lnd%qv_s(jc,jb) = MIN (p_diag_lnd%qv_s(jc,jb), &
              &                    p_prog_now%tracer(jc,nlev,jb,iqv))
          END DO

          DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g(jc,jb) = p_prog_lnd_now%t_g(jc,jb)
          ENDDO


          DO jt = 1, ntiles_total+ntiles_water

            DO jc = i_startidx, i_endidx
              p_prog_lnd_now%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
              p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
              p_diag_lnd%qv_s_t(jc,jb,jt) = p_diag_lnd%qv_s(jc,jb)
            ENDDO
          ENDDO
        END IF  ! init_mode == MODE_IFSANA
      ENDIF

      ! Copy t_g to t_seasfc for idealized cases with surface scheme (would be undefined otherwise)
      IF (ltestcase .AND. atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        DO jc = i_startidx, i_endidx
          p_diag_lnd%t_seasfc(jc,jb) = p_prog_lnd_now%t_g(jc,jb)
        ENDDO
      ENDIF

    END DO
    CALL message('mo_nwp_phy_init:', 'initialized surface temp and humidity')

  ELSE  ! in case of restart
    !
    ! necessary, because only t_g(nnow_rcf) is written to the restart file
    ! with the following copy statement the ocean points of t_g(nnew_rcf) are
    ! filled with the correct values.

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &  i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
      ENDDO
      DO jt = 1, ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx
          p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g_t(jc,jb,jt)
        ENDDO
      ENDDO
    ENDDO

    IF (ltestcase .AND. nh_test_name == 'RCE' .AND. atm_phy_nwp_config(jg)%inwp_turb/=ismag) THEN !
     DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &  i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g (jc,jb) = th_cbl(1)
          p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
          p_diag_lnd%qv_s    (jc,jb) = &
          & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
        ENDDO
     ENDDO
    ENDIF

  END IF



  ! index of first half level with height >= 1500 m (above boundary layer)
  k1500m = 1

  !--------------------------------------------------------------
  !>reference pressure according to U.S. standard atmosphere
  ! (with the caveat that the stratosphere is assumed isothermal, which does not hurt
  !  because pref is used for determining model level indices referring to pressures
  !  >= 60 hPa)
  !--------------------------------------------------------------
  ttropo = t00 + dtdz_standardatm*htropo
  ptropo = p0sl_bg*(ttropo/t00)**(-grav/(rd*dtdz_standardatm))
  DO jk = nlev, 1, -1
    jk1 = jk + nshift
    IF ( vct_a(jk1) >= 1500.0_wp .AND. vct_a(jk1+1) < 1500.0_wp ) THEN
      k1500m = jk  ! not jk1 (may result in out-of-bounds)!
    END IF
    zfull = 0.5_wp*(vct_a(jk1) + vct_a(jk1+1))
    IF (zfull < htropo) THEN
      temp = t00 + dtdz_standardatm*zfull
      pref(jk) = p0sl_bg*(temp/t00)**(-grav/(rd*dtdz_standardatm))
    ELSE
      pref(jk) = ptropo*EXP(-grav*(zfull-htropo)/(rd*ttropo))
    ENDIF
  ENDDO


  ! start filling phy_params
  !
  !--------------------------------------------------------------
  !< characteristic gridlength needed by convection and turbulence
  !--------------------------------------------------------------
!   CALL sphere_cell_mean_char_length (p_patch%n_patch_cells_g, phy_params%mean_charlen)
  ! read it directly from the patch%geometry_info
  phy_params%mean_charlen = p_patch%geometry_info%mean_characteristic_length
!   write(0,*) "=============================================="
!   write(0,*) "mean_charlen=", phy_params%mean_charlen, &
!     & p_patch%geometry_info%mean_characteristic_length
!   write(0,*) "=============================================="

  ! compute level index corresponding to the HAG of the 60hPa level 
  ! (currently only needed by mo_nwp_diagnosis:cal_cape_cin) 
  phy_params%k060=1
  DO jk=nlev,1,-1
    IF(pref(jk) >  60.e2_wp) phy_params%k060=jk
  ENDDO


  !------------------------------------------
  !< call for cloud microphysics
  !------------------------------------------

  SELECT CASE ( atm_phy_nwp_config(jg)%inwp_gscp )

  CASE (1,2,3)  ! cloud microphysics from COSMO (V 5.0)
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics')
    CALL gscp_set_coefficients(tune_zceff_min = tune_zceff_min,               &
      &                        tune_v0snow    = tune_v0snow,                  &
      &                        tune_zvz0i     = tune_zvz0i,                   &
      &                        tune_mu_rain   = atm_phy_nwp_config(1)%mu_rain )

  CASE (4) !two moment micrphysics
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics: two-moment')

    IF (jg == 1) CALL two_moment_mcrph_init(igscp=atm_phy_nwp_config(jg)%inwp_gscp, msg_level=msg_level )

    IF (linit_mode) THEN ! Initial condition for number densities
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
       DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
               &                i_startidx, i_endidx, rl_start, rl_end)
          DO jk=1,nlev
             DO jc=i_startidx,i_endidx
                p_prog_now%tracer(jc,jk,jb,iqnc) = set_qnc(p_prog_now%tracer(jc,jk,jb,iqc))
                p_prog_now%tracer(jc,jk,jb,iqnr) = set_qnr(p_prog_now%tracer(jc,jk,jb,iqr))
                p_prog_now%tracer(jc,jk,jb,iqni) = set_qni(p_prog_now%tracer(jc,jk,jb,iqi))
                p_prog_now%tracer(jc,jk,jb,iqns) = set_qns(p_prog_now%tracer(jc,jk,jb,iqs))
                p_prog_now%tracer(jc,jk,jb,iqng) = set_qng(p_prog_now%tracer(jc,jk,jb,iqg))
             END DO
          END DO
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    END IF

  CASE (5) !two moment micrphysics
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics: two-moment')

    IF (jg == 1) CALL two_moment_mcrph_init(atm_phy_nwp_config(jg)%inwp_gscp,&
         &                                  N_cn0,z0_nccn,z1e_nccn,N_in0,z0_nin,z1e_nin,msg_level)

    IF (linit_mode) THEN ! Initial condition for number densities
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
       DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
               &                i_startidx, i_endidx, rl_start, rl_end)
          DO jk=1,nlev
             DO jc=i_startidx,i_endidx
                p_prog_now%tracer(jc,jk,jb,iqnc) = set_qnc(p_prog_now%tracer(jc,jk,jb,iqc))
                p_prog_now%tracer(jc,jk,jb,iqnr) = set_qnr(p_prog_now%tracer(jc,jk,jb,iqr))
                p_prog_now%tracer(jc,jk,jb,iqni) = set_qni(p_prog_now%tracer(jc,jk,jb,iqi))
                p_prog_now%tracer(jc,jk,jb,iqns) = set_qns(p_prog_now%tracer(jc,jk,jb,iqs))
                p_prog_now%tracer(jc,jk,jb,iqng) = set_qng(p_prog_now%tracer(jc,jk,jb,iqg))
             END DO
          END DO
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    END IF

    IF (linit_mode) THEN ! Initial condition for CCN and IN fields
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jk1,jc,i_startidx,i_endidx,zfull) ICON_OMP_GUIDED_SCHEDULE
       DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
               &                i_startidx, i_endidx, rl_start, rl_end)
          DO jk=1,nlev
             DO jc=i_startidx,i_endidx
                jk1 = jk + nshift
                zfull = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
                IF(zfull > z0_nccn) THEN
                   p_prog_now%tracer(jc,jk,jb,inccn) = N_cn0*exp((z0_nccn-zfull)/z1e_nccn)
                ELSE
                   p_prog_now%tracer(jc,jk,jb,inccn) = N_cn0
                END IF
                IF(zfull > z0_nin) THEN
                   p_prog_now%tracer(jc,jk,jb,ininpot)  = N_in0*exp((z0_nin -zfull)/z1e_nin)
                ELSE
                   p_prog_now%tracer(jc,jk,jb,ininpot)  = N_in0
                END IF
             END DO
          END DO
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    END IF
  CASE (6) ! two-moment scheme with prognostic cloud droplet number
           ! and chemical composition taken from the ART extension
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics: ART two-moment')

    IF (jg == 1) CALL art_clouds_interface_2mom_init(msg_level)

    IF (linit_mode) THEN ! Initial condition for number densities
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
       DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
               &                i_startidx, i_endidx, rl_start, rl_end)
          DO jk=1,nlev
             DO jc=i_startidx,i_endidx
                p_prog_now%tracer(jc,jk,jb,iqnc) = set_qnc(p_prog_now%tracer(jc,jk,jb,iqc))
                p_prog_now%tracer(jc,jk,jb,iqnr) = set_qnr(p_prog_now%tracer(jc,jk,jb,iqr))
                p_prog_now%tracer(jc,jk,jb,iqni) = set_qni(p_prog_now%tracer(jc,jk,jb,iqi))
                p_prog_now%tracer(jc,jk,jb,iqns) = set_qns(p_prog_now%tracer(jc,jk,jb,iqs))
                p_prog_now%tracer(jc,jk,jb,iqng) = set_qng(p_prog_now%tracer(jc,jk,jb,iqg))
             END DO
          END DO
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    END IF
  END SELECT

  ! Compute lookup tables for aerosol-microphysics coupling
  IF (jg == 1 .AND. (atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .OR. icpl_aero_conv > 0)) &
    CALL lookupcreate_segalkhain()

  !------------------------------------------
  !< radiation
  !------------------------------------------
  SELECT CASE ( atm_phy_nwp_config(jg)%inwp_radiation )
  CASE (1, 3)

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init RRTM')

    SELECT CASE ( irad_aero )
    ! Note (GZ): irad_aero=2 does no action but is the default in radiation_nml
    ! and therefore should not cause the model to stop
    CASE (0,2,5,6,9)
      !ok
    CASE DEFAULT
      CALL finish('mo_nwp_phy_init: init_nwp_phy',  &
        &      'Wrong irad_aero. For RRTM radiation, this irad_aero is not implemented.')
    END SELECT

!    prm_diag%lfglac (:,:) = ext_data%atm%soiltyp(:,:) == 1  !soiltyp=ice

    ! solar flux (W/m2) in 14 SW bands
    ssi_radt(:) = ssi_amip(:)
    ! solar constant (W/m2)
    tsi_radt    = SUM(ssi_radt(:))


    !--------------------------------------------------
    !< set conditions for Aqua planet or RCE experiment
    !--------------------------------------------------
    IF ( nh_test_name == 'APE_nwp' .OR. nh_test_name == 'dcmip_tc_52' ) THEN
      ssi_radt(:) = ssi_radt(:)*1365._wp/tsi_radt
      tsi_radt = 1365._wp
    ENDIF  ! APE

    IF ( nh_test_name == 'RCE') THEN
      ! solar flux (W/m2) in 14 SW bands
      scale_fac = sol_const/1361.371_wp ! computed relative to amip (1361)
      ssi_radt(:) = scale_fac*ssi_amip(:)
      ! solar constant (W/m2)
      tsi_radt    = SUM(ssi_radt(:))
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_radiation == 1) THEN ! RRTM init
      CALL setup_srtm
      CALL lrtm_setup(lrtm_filename)
      CALL setup_newcld_optics(cldopt_filename)
    ELSE   ! PSRAD init
      CALL setup_psrad
      CALL setup_cloud_optics
      CALL setup_lrtm
      CALL lrtm_setup(lrtm_filename) ! ** necessary because of incorrect USE statements in mo_psrad_lrtm_gas_optics **
      CALL setup_srtm_psrad
    ENDIF

    rl_start = 1  ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk,zprat,zpres,lland,lglac,zn1,&
!$OMP zn2,zcdnc,dz1,dz2,dz3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &  i_startidx, i_endidx, rl_start, rl_end)

      ! Initialize cloud droplet number concentration (acdnc)
      ! like in mo_echam_phy_init.f90
      DO jk = 1,nlev
        ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
        DO jc = i_startidx, i_endidx
          zpres = p0ref * (p_metrics%exner_ref_mc(jc,jk,jb))**(cpd/rd)
          zprat=(MIN(8._wp,80000._wp/zpres))**2

          lland = ext_data%atm%llsm_atm_c(jc,jb)
          lglac = ext_data%atm%soiltyp(jc,jb) == 1
          IF (lland.AND.(.NOT.lglac)) THEN
            zn1= 50._wp
            zn2=220._wp
          ELSE
            zn1= 50._wp
            zn2= 80._wp
          ENDIF
          IF (zpres < 80000._wp) THEN
            zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
          ELSE
            zcdnc=zn2*1.e6_wp
          ENDIF
          prm_diag%acdnc(jc,jk,jb) = zcdnc
          IF ( nh_test_name == 'RCE' ) THEN
            !--- computation of reference pressure field from the reference exner field
            zrefpres(jc,jk,jb) = p0ref * (p_metrics%exner_ref_mc(jc,jk,jb))**(cpd/rd)
            ! here we choose to use temp to compute sfc pres instead of tempv
            zreftemp(jc,jk,jb) = p_metrics%theta_ref_mc(jc,jk,jb)*p_metrics%exner_ref_mc(jc,jk,jb)

          END IF
        END DO !jc
      END DO   !jk
      IF ( nh_test_name == 'RCE') THEN
        ! a ref press field needs to be computed for testcases with a
        ! constant ozone.  the reference field allows the ozone to be
        ! interpolated at a restart without changing due to a changing p field.
        DO jc = i_startidx,i_endidx
          ! we also need pres at interface levels; first we need sfc press...
          ! Height differences between surface and third-lowest main level
          dz1 = p_metrics%z_ifc(jc,nlev,jb)   - p_metrics%z_ifc(jc,nlevp1,jb)
          dz2 = p_metrics%z_ifc(jc,nlev-1,jb) - p_metrics%z_ifc(jc,nlev,jb)
          dz3 = p_metrics%z_mc (jc,nlev-2,jb) - p_metrics%z_ifc(jc,nlev-1,jb)
          ! Compute surface pressure starting from three lowest levels
          zpres_sfc(jc,jb) = p0ref * EXP( cpd_o_rd*LOG(p_metrics%exner_ref_mc(jc,nlev-2,jb))  + &
                             grav_o_rd*(dz1/zreftemp(jc,nlev,jb) + dz2/zreftemp(jc,nlev-1,jb) + &
                             dz3/zreftemp(jc,nlev-2,jb)) )

          zpres_ifc(jc,nlevp1,jb) = zpres_sfc(jc,jb)
        END DO !jc

        ! compute interface from nlev-1 to TOA
        DO jk = nlev,2,-1
          DO jc = 1, i_endidx
            ! pressure at interface levels
            zpres_ifc(jc,jk,jb) = SQRT(zrefpres(jc,jk,jb)*zrefpres(jc,jk-1,jb) )
          END DO
          DO jc = 1, i_endidx !pres at top ifc = pres at top model lev
            zpres_ifc(jc,1,jb) = zrefpres(jc,1,jb)
          END DO
        END DO
      END IF


    !------------------------------------------
    ! APE ozone profile, vertical setting needed only once for NH
    !------------------------------------------
      !IF (irad_o3 == io3_ape .AND. linit_mode) THEN
      IF (irad_o3 == io3_ape ) THEN

!        CALL o3_zl2ml(p_patch%nblks_c,p_patch%npromz_c,        & !
!          &           nlev_o3,      nlev,                      & ! vertical levels in/out
!          &           zf_aux,   p_metrics%z_mc,                & ! vertical in/out
!          &           ext_data%atm_td%o3(:,:,:,nmonths),p_prog%tracer(:,:,:,io3))! o3Field in/out

        IF ( nh_test_name == 'RCE' ) THEN
          CALL o3_pl2ml ( kproma= i_endidx, kbdim=nproma,  &
            & nlev_pres = nlev_o3,klev= nlev ,             &
            & pfoz = ext_data%atm_td%pfoz(:),              &
            & phoz = ext_data%atm_td%phoz(:),              &! in o3-levs
            & ppf = zrefpres (:,:,jb),                  &! in  pres
            & pph = zpres_ifc(:,:,jb),               &! in  pres_halfl
            & o3_time_int = ext_data%atm_td%o3(:,:,jb,nmonths),     &! in
            & o3_clim     = ext_data%atm%o3(:,:,jb) )         ! OUT
        ELSE ! default behaviour
          CALL o3_pl2ml ( kproma= i_endidx, kbdim=nproma,  &
            & nlev_pres = nlev_o3,klev= nlev ,             &
            & pfoz = ext_data%atm_td%pfoz(:),              &
            & phoz = ext_data%atm_td%phoz(:),              &! in o3-levs
            & ppf = p_diag%pres (:,:,jb),                  &! in  pres
            & pph = p_diag%pres_ifc(:,:,jb),               &! in  pres_halfl
            & o3_time_int = ext_data%atm_td%o3(:,:,jb,nmonths),     &! in
            & o3_clim     = ext_data%atm%o3(:,:,jb) )         ! OUT
        ENDIF
      ENDIF

    ENDDO      !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF ( irad_aero == 5 ) THEN

      CALL init_aerosol_props_tanre_rrtm

      CALL init_aerosol_dstrb_tanre (        &
        & kbdim    = nproma,                 & !in
        & pt_patch = p_patch,                & !in
        & aersea   = prm_diag%aersea,        & !out
        & aerlan   = prm_diag%aerlan,        & !out
        & aerurb   = prm_diag%aerurb,        & !out
        & aerdes   = prm_diag%aerdes )         !out

    ELSEIF ( irad_aero == 6 .OR. irad_aero == 9) THEN

      CALL init_aerosol_props_tegen_rrtm

    ELSE

      zaea_rrtm(:,:) = 0.0_wp
      zaes_rrtm(:,:) = 0.0_wp
      zaeg_rrtm(:,:) = 0.0_wp

    ENDIF

    DO ist = 1, UBOUND(csalbw,1)
      rad_csalbw(ist) = csalbw(ist) / (2.0_wp * zml_soil(1))
    ENDDO

  CASE (2)

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init Ritter Geleyn')

    ! Note (GZ): irad_aero=2 does no action but is the default in radiation_nml
    ! and therefore should not cause the model to stop
    SELECT CASE ( irad_aero )
    CASE (0,2,5,6,9)
      !ok
    CASE DEFAULT
      CALL finish('mo_nwp_phy_init: init_nwp_phy',  &
        &      'Wrong irad_aero. For Ritter-Geleyn radiation, this irad_aero is not implemented.')
    END SELECT

    ! solar flux (W/m2) in 14 SW bands
    ssi_radt(:) = ssi_amip(:)
    ! solar constant (W/m2)
    tsi_radt    = SUM(ssi_radt(:))

    IF ( nh_test_name == 'RCE' ) THEN
      tsi_radt = 0._wp
      ! solar flux (W/m2) in 14 SW bands
      scale_fac = sol_const/1361.371_wp ! computed relative to amip (1361)
      ssi_radt(:) = scale_fac*ssi_amip(:)
      ! solar constant (W/m2)
      tsi_radt    = SUM(ssi_radt(:))
    ENDIF

    !------------------------------------------
    !< set conditions for Aqua planet experiment
    !------------------------------------------
    IF ( nh_test_name == 'APE_nwp' .OR. nh_test_name == 'dcmip_tc_52' ) THEN
      ssi_radt(:) = ssi_radt(:)*1365._wp/tsi_radt
      tsi_radt = 1365._wp
    ENDIF

    CALL rad_aibi

    zaef_rg(:,:)= 0.0_wp

    IF ( irad_aero == 5 ) THEN

      CALL init_aerosol_props_tanre_rg

      CALL init_aerosol_dstrb_tanre (        &
        & kbdim    = nproma,                 & !in
        & pt_patch = p_patch,                & !in
        & aersea   = prm_diag%aersea,        & !out
        & aerlan   = prm_diag%aerlan,        & !out
        & aerurb   = prm_diag%aerurb,        & !out
        & aerdes   = prm_diag%aerdes )         !out

    ELSEIF ( irad_aero == 6 .OR. irad_aero == 9) THEN

      CALL init_aerosol_props_tegen_rg

    ELSE

      zaea_rg(:,:) = 0.0_wp
      zaes_rg(:,:) = 0.0_wp
      zaeg_rg(:,:) = 0.0_wp

    ENDIF


    !------------------------------------------
    ! APE ozone profile, vertical setting needed only once for NH
    !------------------------------------------
    IF (irad_o3 == io3_ape .AND. linit_mode ) THEN

      rl_start = 1  ! Initialization should be done for all points
      rl_end   = min_rlcell

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &  i_startidx, i_endidx, rl_start, rl_end)

        CALL o3_pl2ml ( kproma= i_endidx, kbdim=nproma,  &
          & nlev_pres = nlev_o3,klev= nlev ,             &
          & pfoz = ext_data%atm_td%pfoz(:),              &
          & phoz = ext_data%atm_td%phoz(:),              &! in o3-levs
          & ppf = p_diag%pres (:,:,jb),                  &! in  pres
          & pph = p_diag%pres_ifc(:,:,jb),               &! in  pres_halfl
          & o3_time_int = ext_data%atm_td%o3(:,:,jb,nmonths),     &! in
          & o3_clim     = ext_data%atm%o3(:,:,jb) )         ! OUT

      ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF ! (irad_o3 == io3_ape)

    DO ist = 1, UBOUND(csalbw,1)
      rad_csalbw(ist) = csalbw(ist) / (2.0_wp * zml_soil(1))
    ENDDO

  END SELECT !inwp_radiation

  IF ( nh_test_name == 'RCE' ) THEN
    DEALLOCATE (zrefpres)
    DEALLOCATE (zreftemp)
    DEALLOCATE (zpres_sfc)
    DEALLOCATE (zpres_ifc)
  END IF

  !----------------------------------------------------------------
  !< initializations needed both for convection and inwp_cldcover=1
  !----------------------------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 .OR. &
    &  atm_phy_nwp_config(jg)%inwp_cldcover == 1   .OR. &
    &  atm_phy_nwp_config(jg)%inwp_surface == 1    .OR. &
    &  atm_phy_nwp_config(jg)%inwp_turb == iedmf )     THEN

    !This has to be done here because not only convection, but also inwp_cldcover == 1
    !uses mo_cufunctions's foealfa. Therefore, the parameters of the function foealfa
    !have to be initialized by calls of sucst and su_yoethf.

    ! get current date in iso-format "yyyymmddThhmmssZ" (String)
    CALL datetimeToString(ini_date, datetime_string)
    ! convert first 8 characters to integer (yyyy-mm-dd)
    WRITE (yyyymmdd, '(a,a,a)')  datetime_string(1:4), datetime_string(6:7), datetime_string(9:10)
    READ  (yyyymmdd,'(i8)') icur_date

    CALL sucst(54,icur_date,0,0)
    CALL su_yoethf

  ENDIF


  !------------------------------------------
  !< call for convection
  !------------------------------------------
  ! initialization.
  ! k800, k400 will be used for inwp_convection==0 as well. 
  ! Thus we need to make sure that they are initialized.
  prm_diag%k850(:,:) = nlev
  prm_diag%k950(:,:) = nlev
  prm_diag%k800(:,:) = nlev
  prm_diag%k400(:,:) = nlev

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 .OR. &
    &  atm_phy_nwp_config(jg)%inwp_turb == iedmf )     THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init convection')

    ! Please take care for scale-dependent initializations!
    ! rsltn = Average mesh size of ICON grid
    ! needed for RTAU - CAPE calculation
    ! adapted for more general gemoetries
    rsltn = p_patch%geometry_info%mean_characteristic_length


!    WRITE(message_text,'(i3,i10,f20.10)') jg, nsmax, phy_params%mean_charlen
!    CALL message('nwp_phy_init, nsmax=', TRIM(message_text))

    lshallow = atm_phy_nwp_config(jg)%lshallowconv_only
    ldetrain_prec = atm_phy_nwp_config(jg)%ldetrain_conv_prec
    CALL sucumf(rsltn,nlev,pref,phy_params,lshallow,ldetrain_prec)
    CALL suphli
    CALL suvdf
    CALL suvdfs
    CALL sucldp


    ! Initialize fields k850 and k950, which are required for computing the
    ! convective contribution to wind gusts
    rl_start = 1  ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! height of 850 and 950hPa surface for US standard atmosphere in m
    ! For derivation, see documentation of US standard atmosphere
    h850_standard = 1457.235199_wp
    h950_standard = 540.3130233_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,hag,zpres,zpres0) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &  i_startidx, i_endidx, rl_start, rl_end)

      DO jc=i_startidx, i_endidx
        ! initialization
!!$        prm_diag%k850(jc,jb) = nlev
!!$        prm_diag%k950(jc,jb) = nlev
        DO jk=nlev, 1, -1
          ! height above ground
          hag = p_metrics%z_mc(jc,jk,jb)-ext_data%atm%topography_c(jc,jb)

          IF (hag < h850_standard) THEN
            prm_diag%k850(jc,jb) = jk
          ENDIF
          IF (hag < h950_standard) THEN
            prm_diag%k950(jc,jb) = jk
          ENDIF
        ENDDO
        ! security measure
        prm_diag%k950(jc,jb) = MAX(prm_diag%k950(jc,jb),2)
        prm_diag%k850(jc,jb) = MAX(prm_diag%k850(jc,jb),2)

        ! analogous initialization of k800 and k400, based on reference pressure
        ! because this is more meaningful for k400 in the presence of very high orography
        zpres0 = p0ref * (p_metrics%exner_ref_mc(jc,nlev,jb))**(cpd/rd)
!!$        prm_diag%k800(jc,jb) = nlev
!!$        prm_diag%k400(jc,jb) = nlev
        DO jk=nlev-1, 2, -1
          zpres = p0ref * (p_metrics%exner_ref_mc(jc,jk,jb))**(cpd/rd)
          IF (zpres/zpres0 >= pr800) prm_diag%k800(jc,jb) = jk
          IF (zpres/zpres0 >= pr400*SQRT(p0ref/zpres0)) THEN
            prm_diag%k400(jc,jb) = jk
          ELSE
            EXIT
          ENDIF
        ENDDO

      ENDDO  ! jc
    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL message('mo_nwp_phy_init:', 'convection initialized')
  ELSE
    ! initialize parameters that are accessed outside the convection scheme
    phy_params%rcucov           = 0._wp
    phy_params%rcucov_trop      = 0._wp
    phy_params%rhebc_land       = 0._wp
    phy_params%rhebc_ocean      = 0._wp
    phy_params%rhebc_land_trop  = 0._wp
    phy_params%rhebc_ocean_trop = 0._wp
    phy_params%entrorg          = 0._wp
    phy_params%texc             = 0._wp
    phy_params%qexc             = 0._wp
  ENDIF


  !------------------------------------------
  !< surface initialization (including seaice)
  !------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN  ! TERRA
    IF (linit_mode) THEN
      CALL nwp_surface_init(p_patch, ext_data, p_prog_lnd_now, p_prog_lnd_new, &
        &                   p_prog_wtr_now, p_prog_wtr_new, p_diag_lnd, p_diag)
    ELSE
      IF ( lsnowtile ) THEN ! restart mode with snowtiles
        CALL init_snowtile_lists(p_patch, ext_data, p_diag_lnd)
      ENDIF
    ENDIF

    ! Copy timelevel now to timelevel new for land state. This has no impact on the prognostic
    ! results but ensures that the output over non-prognostic grid points (water) is the same
    ! for even and odd multiples of the advection time step
    CALL copy_lnd_prog_now2new(p_patch, p_prog_lnd_now, p_prog_lnd_new)
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init TERRA')
  END IF


  !------------------------------------------
  !< setup for turbulence
  !------------------------------------------

  ! initialize gz0 (roughness length * g)
  !
  IF ( ANY( (/icosmo,igme,ismag/)==atm_phy_nwp_config(jg)%inwp_turb ) .AND. linit_mode ) THEN


    ! gz0 is initialized if we do not start from an own first guess
    IF (lturb_init) THEN

      IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init roughness length')

      IF (turbdiff_config(jg)%lconst_z0) THEN
        ! constant z0 for idealized tests
        prm_diag%gz0(:,:) = grav * turbdiff_config(jg)%const_z0

      ELSE IF (atm_phy_nwp_config(jg)%itype_z0 == 1) THEN
        ! default
        prm_diag%gz0(:,:) = grav * ext_data%atm%z0(:,:)

      ELSE IF (atm_phy_nwp_config(jg)%itype_z0 >= 2) THEN

        rl_start = grf_bdywidth_c + 1 ! land-cover classes are not set for nest-boundary points
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


        ! Scaling factor for SSO contribution to roughness length ("Erdmann Heise formula")
        fact_z0rough = 1.e-5_wp*ATAN(phy_params%mean_charlen/2250._wp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,ic,jt,i_startidx,i_endidx,lc_class,gz0) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, rl_start, rl_end)

          ! specify land-cover-related roughness length over land points
          ! note:  water points are set in turbdiff
          gz0(:) = 0._wp

          DO jt = 1, ntiles_total
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, ext_data%atm%lp_count(jb)
              jc = ext_data%atm%idx_lst_lp(ic,jb)
              lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
              gz0(jc) = gz0(jc) + ext_data%atm%frac_t(jc,jb,jt) * grav * (             &
               (1._wp-p_diag_lnd%snowfrac_t(jc,jb,jt))*ext_data%atm%z0_lcc(lc_class)+  &
                p_diag_lnd%snowfrac_t(jc,jb,jt)*0.5_wp*ext_data%atm%z0_lcc(i_lc_si) ) ! i_lc_si = snow/ice class
            ENDDO
          ENDDO
          IF (atm_phy_nwp_config(jg)%itype_z0 == 3) THEN
            DO ic = 1, ext_data%atm%lp_count(jb)
              jc = ext_data%atm%idx_lst_lp(ic,jb)
              gz0(jc) = gz0(jc) + grav*MIN(fact_z0rough*ext_data%atm%sso_stdh_raw(jc,jb)**2,7.5_wp)
            ENDDO
          ENDIF
          DO jt = ntiles_total+1, ntiles_total+ntiles_water ! required if there are mixed land-water points
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, ext_data%atm%lp_count(jb)
              jc = ext_data%atm%idx_lst_lp(ic,jb)
              lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
              gz0(jc) = gz0(jc) + ext_data%atm%frac_t(jc,jb,jt) * grav*ext_data%atm%z0_lcc(lc_class)
            ENDDO
          ENDDO
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data%atm%lp_count(jb)
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            prm_diag%gz0(jc,jb) = gz0(jc)
          ENDDO
        ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
      ENDIF  !initialize gz0

    END IF

  ENDIF

  IF ( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) THEN

    ! allocate and init implicit weights for tridiagonal solver
    ALLOCATE( turbdiff_config(jg)%impl_weight(nlevp1), &
              STAT=istatus )
    IF(istatus/=SUCCESS)THEN
      CALL finish (TRIM(routine), &
                 'allocation of impl_weight failed')
    ENDIF

    CALL get_turbdiff_param(jg)

    ! using an over implicit value (impl_s) near surface,
    ! reduced to in general slightly off-centered value (impl_t)
    ! in about 1500 m height
    DO jk = 1, k1500m
      impl_weight(jk) = impl_t
    END DO
    DO jk = k1500m+1, nlev
      impl_weight(jk) = impl_t &
                      + (impl_s-impl_t) * (jk-k1500m) / REAL(nlev-k1500m, wp)
    END DO
    impl_weight(nlevp1) = impl_s

  ENDIF

  ! Initialize turbulence models
  !
  IF ( (atm_phy_nwp_config(jg)%inwp_turb == icosmo) .AND. linit_mode ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init COSMO turbulence')

    rl_start = 1 ! Initialization is done also for nest boundary points
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,ic,jc,jt, &
!$OMP            ltkeinp_loc,lgz0inp_loc,nlevcm,l_hori) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      IF (lturb_init) THEN

        ltkeinp_loc = .FALSE.  ! initialize TKE field
        lgz0inp_loc = .FALSE.  ! initialize gz0 field (water points only)

      ELSE
        !
        ! TKE and gz0 are not re-initialized, but re-used from the first guess
        !
        ltkeinp_loc = .TRUE.   ! do NOT re-initialize TKE field (read from FG)
        lgz0inp_loc = .TRUE.   ! do NOT re-initialize gz0 field (read from FG)


        ! Note that TKE in turbtran/turbdiff is defined as the turbulence velocity scale
        ! TVS=SQRT(2*TKE)
        !
        DO jk =1,nlevp1
          p_prog_now%tke(i_startidx:i_endidx,jk,jb)= SQRT(2.0_wp                        &
            &                                * p_prog_now%tke(i_startidx:i_endidx,jk,jb))
        ENDDO
      ENDIF

      l_hori(i_startidx:i_endidx)=phy_params%mean_charlen

      nlevcm = nlevp1

!MR: There should be an initialization for each tile, or the initialization can be
!    executed for 'turbtran' and 'turbdiff' within a single CALL of 'organize_turbdiff'!

      ! turbtran
      CALL organize_turbdiff( &
        &  iini=1, lturatm=.FALSE., ltursrf=.TRUE. , lstfnct=.TRUE. ,         & !only surface-layer turbulence
        &          lnsfdia=.TRUE. , ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, & !including near-surface diagnostics
        &  itnd=0, lum_dif=.FALSE., lvm_dif=.FALSE., lscadif=.FALSE.,         & !and surface-flux calculations
        &          lsrflux=.TRUE. , lsfluse=.FALSE., lqvcrst=.FALSE.,         & !but without vertical diffusion calculation
!
        &  dt_var=atm_phy_nwp_config(jg)%dt_fastphy, &
        &  dt_tke=atm_phy_nwp_config(jg)%dt_fastphy, &
        &  nprv=1, ntur=1, ntim=1, &
        &  ie=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm, &
        &  i_st=i_startidx, i_en=i_endidx, i_stp=i_startidx, i_enp=i_endidx, &
        &  l_hori=l_hori, hhl=p_metrics%z_ifc(:,:,jb), &
        &  dp0=p_diag%dpres_mc(:,:,jb), &
        &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
        &  h_ice=p_prog_wtr_now%h_ice(:,jb), gz0=prm_diag%gz0(:,jb), &
        &  sai=ext_data%atm%sai(:,jb), &
        &  t_g=p_prog_lnd_now%t_g(:,jb), ps=p_diag%pres_sfc(:,jb), &
        &  qv_s=p_diag_lnd%qv_s(:,jb), &
        &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), &
        &  t=p_diag%temp(:,:,jb), prs=p_diag%pres(:,:,jb), &
        &  qv=p_prog_now%tracer(:,:,jb,iqv), qc=p_prog_now%tracer(:,:,jb,iqc), &
        &  tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
        &  tvm=prm_diag%tvm(:,jb), tvh=prm_diag%tvh(:,jb), tkr=prm_diag%tkr(:,jb), &
        &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
        &  tke=p_prog_now%tke(:,:,jb), &
        &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh(:,:,jb), &
        &  rcld=prm_diag%rcld(:,:,jb), &
        &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), &
        &  td_2m=prm_diag%td_2m(:,jb), rh_2m=prm_diag%rh_2m(:,jb), &
        &  u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb), &
        &  shfl_s=prm_diag%shfl_s(:,jb), qvfl_s=prm_diag%qhfl_s(:,jb), &
        &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

      prm_diag%lhfl_s(i_startidx:i_endidx,jb) = &
        &  prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v


      ! turbdiff
      CALL organize_turbdiff( &
        &  iini=1, lturatm=.TRUE. , ltursrf=.FALSE., lstfnct=.TRUE. ,         & !atmosph. turbulence and vertical diffusion
        &          lnsfdia=.FALSE., ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, & !but no surface-layer turbulence (turbtran)
        &  itnd=0, lum_dif=.TRUE. , lvm_dif=.TRUE. , lscadif=.TRUE. ,         & !and thus (implicitly) neither surface-layer diagn.
        &          lsrflux=.FALSE., lsfluse=lsflcnd, lqvcrst=.FALSE.,         & !nor surface-flux calculation (both in turbtran)
!MR: turbulent diffusion can be switched off for initialization!
        &  dt_var=atm_phy_nwp_config(jg)%dt_fastphy, &
        &  dt_tke=atm_phy_nwp_config(jg)%dt_fastphy, &
        &  nprv=1, ntur=1, ntim=1, &
        &  ie=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm, &
        &  i_st=i_startidx, i_en=i_endidx, i_stp=i_startidx, i_enp=i_endidx, &
        &  l_hori=l_hori, hhl=p_metrics%z_ifc(:,:,jb), &
        &  dp0=p_diag%dpres_mc(:,:,jb), &
        &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
        &  h_ice=p_prog_wtr_now%h_ice(:,jb), gz0=prm_diag%gz0(:,jb), &
        &  sai=ext_data%atm%sai(:,jb), &
        &  t_g=p_prog_lnd_now%t_g(:,jb), ps=p_diag%pres_sfc(:,jb), &
        &  qv_s=p_diag_lnd%qv_s(:,jb), &
        &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), &
        &  w=p_prog_now%w(:,:,jb), &
        &  t=p_diag%temp(:,:,jb), prs=p_diag%pres(:,:,jb), &
        &  rho=p_prog_now%rho(:,:,jb), epr=p_prog_now%exner(:,:,jb), &
        &  qv=p_prog_now%tracer(:,:,jb,iqv), qc=p_prog_now%tracer(:,:,jb,iqc), &
!         &  ptr=???, &  ! for the diffusion of additional tracer variables!
        &  tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
        &  tvm=prm_diag%tvm(:,jb), tvh=prm_diag%tvh(:,jb), &
        &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
        &  tke=p_prog_now%tke(:,:,jb), &
        &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh(:,:,jb), &
        &  rcld=prm_diag%rcld(:,:,jb), &
        &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), &
        &  v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb), &
        &  tketens=prm_nwp_tend%ddt_tke(:,:,jb), &
        &  ut_sso=REAL(prm_nwp_tend%ddt_u_sso(:,:,jb),wp), vt_sso=REAL(prm_nwp_tend%ddt_v_sso(:,:,jb),wp), &
        &  shfl_s=prm_diag%shfl_s(:,jb), qvfl_s=prm_diag%qhfl_s(:,jb), &
        &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

      ! preparation for concentration boundary condition. Usually inactive for standard ICON runs.
      IF ( .NOT. lsflcnd ) THEN
        prm_diag%lhfl_s(i_startidx:i_endidx,jb) = &
          &  prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v
      END IF


      ! tile-specific quantities needed by turbtran
      !
      DO jt = 1, ntiles_total+ntiles_water
        prm_diag%gz0_t   (:,jb,jt) = prm_diag%gz0(:,jb)
        prm_diag%tvs_s_t (:,jb,jt) = p_prog_now%tke(:,nlevp1,jb)  !here: SQRT(2*TKE)
        prm_diag%tkvm_s_t(:,jb,jt) = prm_diag%tkvm(:,nlevp1,jb)
        prm_diag%tkvh_s_t(:,jb,jt) = prm_diag%tkvh(:,nlevp1,jb)
        prm_diag%tkr_t   (:,jb,jt) = prm_diag%tkr(:,jb)
      ENDDO


      ! Note that TKE in turbtran/turbdiff is defined as the turbulence velocity scale
      ! TVS=SQRT(2*TKE)
      !
      DO jk =1,nlevp1
        p_prog_now%tke(i_startidx:i_endidx,jk,jb)= 0.5_wp                        &
          &                                * (p_prog_now%tke(i_startidx:i_endidx,jk,jb))**2
      ENDDO
    ENDDO  ! jb
!$OMP END DO

!$OMP END PARALLEL

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'Cosmo turbulence initialized')


  ELSE IF (  atm_phy_nwp_config(jg)%inwp_turb == igme) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init GME turbulence')

    rl_start = grf_bdywidth_c + 1 ! land-cover classes are not set for nest-boundary points
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! paranoia: Make sure that rcld is initialized  (needed by cloud cover scheme)
      prm_diag%rcld(:,:,jb) = 0._wp

      DO jt = 1, ntiles_total+ntiles_water
        prm_diag%gz0_t(i_startidx:i_endidx,jb,jt) = prm_diag%gz0(i_startidx:i_endidx,jb)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == ismag .AND. linit_mode ) THEN

    CALL message('mo_nwp_phy_init:', 'init Smagorinsky turbulence')

    IF(atm_phy_nwp_config(jg)%inwp_surface == 0)THEN
      IF (turbdiff_config(jg)%lconst_z0) THEN
        prm_diag%gz0(:,:) = grav * turbdiff_config(jg)%const_z0
      ELSE
        CALL finish (TRIM(routine), 'Only constant roughness length allowed idealized LES cases!')
      END IF
    ELSE
      !Default: Already set above
    END IF

  ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN  !EDMF DUALM
    CALL suct0
    CALL su0phy
    CALL susekf
    CALL susveg
    CALL sussoil

!$OMP PARALLEL
    ! paranoia: Make sure that rcld is initialized  (needed by cloud cover scheme)
    CALL init(prm_diag%rcld(:,:,:))
!$OMP END PARALLEL

  ENDIF


  ! Gravity wave drag scheme
  !
  IF ( atm_phy_nwp_config(jg)%inwp_gwd == 1 ) THEN  ! IFS gwd scheme

    CALL sugwwms(nflevg=nlev, ppref=pref, klaunch=phy_params%klaunch)
    CALL message('mo_nwp_phy_init:', 'non-orog GWs initialized')

  END IF

  ! SSO scheme
  !
  SELECT CASE ( atm_phy_nwp_config(jg)%inwp_sso )
  CASE ( 1 )                                ! COSMO SSO scheme
    IF (jg == 1) CALL sso_cosmo_init_param(tune_gkwake=tune_gkwake, tune_gkdrag=tune_gkdrag, &
                                           tune_gfrcrit=tune_gfrcrit, tune_grcrit=tune_grcrit)
    IF (linit_mode) prm_diag%ktop_envel(:,:) = nlev
  CASE ( 2 )                                ! IFS SSO scheme
    CALL sugwd(nlev, pref, phy_params )
    IF (linit_mode) prm_diag%ktop_envel(:,:) = nlev
  END SELECT

  !  WW diagnostics
  !
  IF ( atm_phy_nwp_config(jg)%inwp_gscp > 0) THEN
    CALL configure_ww(ini_date, jg, nlev, nshift)
  END IF


END SUBROUTINE init_nwp_phy


  SUBROUTINE init_cloud_aero_cpl ( mtime_date, p_patch, p_metrics, ext_data, prm_diag)

    TYPE(datetime),   POINTER               :: mtime_date
    TYPE(t_patch),               INTENT(in) :: p_patch
    TYPE(t_nh_metrics),          INTENT(in) :: p_metrics
    TYPE(t_external_data),       INTENT(in) :: ext_data

    TYPE(t_nwp_phy_diag),        INTENT(inout) :: prm_diag

    INTEGER          :: imo1, imo2
    INTEGER          :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER          :: jb, jc, jg, nlev

    REAL(wp) :: wgt, zncn(nproma, p_patch%nlev)
    
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    TYPE(datetime), POINTER :: mtime_hour
    
    jg = p_patch%id
    nlev = p_patch%nlev

    IF (irad_aero /= 6) RETURN
    IF (atm_phy_nwp_config(jg)%icpl_aero_gscp /= 1 .AND. icpl_aero_conv /= 1) RETURN

    
    mtime_hour => newDatetime(mtime_date)
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0          
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
    call deallocateDatetime(mtime_hour)
    imo1 = current_time_interpolation_weights%month1
    imo2 = current_time_interpolation_weights%month2
    wgt = current_time_interpolation_weights%weight2

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zncn)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
 
      IF (iprog_aero == 0) THEN
        DO jc = i_startidx, i_endidx

          prm_diag%aerosol(jc,iss,jb) = ext_data%atm_td%aer_ss(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_ss(jc,jb,imo2)   - ext_data%atm_td%aer_ss(jc,jb,imo1)   ) * wgt
          prm_diag%aerosol(jc,iorg,jb) = ext_data%atm_td%aer_org(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_org(jc,jb,imo2)  - ext_data%atm_td%aer_org(jc,jb,imo1)  ) * wgt
          prm_diag%aerosol(jc,ibc,jb) = ext_data%atm_td%aer_bc(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_bc(jc,jb,imo2)   - ext_data%atm_td%aer_bc(jc,jb,imo1)   ) * wgt
          prm_diag%aerosol(jc,iso4,jb) = ext_data%atm_td%aer_so4(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_so4(jc,jb,imo2)  - ext_data%atm_td%aer_so4(jc,jb,imo1)  ) * wgt
          prm_diag%aerosol(jc,idu,jb) = ext_data%atm_td%aer_dust(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_dust(jc,jb,imo2) - ext_data%atm_td%aer_dust(jc,jb,imo1) ) * wgt

        ENDDO
      ENDIF

      CALL ncn_from_tau_aerosol_speccnconst (nproma, nlev, i_startidx, i_endidx, nlev, nlev, &
        p_metrics%z_ifc(:,:,jb), prm_diag%aerosol(:,iss,jb), prm_diag%aerosol(:,iso4,jb),    &
        prm_diag%aerosol(:,iorg,jb), prm_diag%aerosol(:,idu,jb), zncn)

      CALL specccn_segalkhain_simple (nproma, i_startidx, i_endidx, zncn(:,nlev), prm_diag%cloud_num(:,jb))

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_cloud_aero_cpl

END MODULE mo_nwp_phy_init

