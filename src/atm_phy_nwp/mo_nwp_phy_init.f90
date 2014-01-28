!>
!!  Initialize the physical schemes at start time
!!
!! @author Kristina  Froehlich, DWD
!!
!! @par Revision History
!! First implementations by Kristina Froehlich, DWD, 2010-070-20
!! Include initialitation of SST for APE experiments by P. Ripodas, DWD,2010-11
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

MODULE mo_nwp_phy_init

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: grav, rd_o_cpd, cpd, p0ref, rd, p0sl_bg, tmelt, &
    &                               dtdz_standardatm, lh_v=>alv
!   USE mo_math_utilities,      ONLY: sphere_cell_mean_char_length
  USE mo_grid_config,         ONLY: grid_sphere_radius
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag,t_nwp_phy_tend
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ext_data_state,      ONLY: nlev_o3, nmonths
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_vertical_coord_table,ONLY: vct_a, vct
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell, min_rlcell_int, zml_soil, io3_ape,  &
    &                               MODE_COMBINED, MODE_IFSANA, MODE_DWDANA, icosmo,&
    &                               igme, iedmf, SUCCESS, MAX_CHAR_LENGTH,          &
    &                               MODE_COSMODE
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltestcase, iqv, iqc, msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, lrtm_filename,              &
    &                               cldopt_filename
  !radiation
  USE mo_newcld_optics,       ONLY: setup_newcld_optics
  USE mo_lrtm_setup,          ONLY: lrtm_setup
  USE mo_radiation_config,    ONLY: ssi, tsi,irad_o3, irad_aero, rad_csalbw 
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

  ! microphysics
  USE mo_gscp_cosmo,          ONLY: hydci_pp_old_init
  USE gscp_hydci_pp,          ONLY: hydci_pp_init
  USE gscp_hydci_pp_ice,      ONLY: hydci_pp_ice_init
  ! convection
  USE mo_cuparameters,        ONLY: sucst,  sucumf,    &
    &                               su_yoethf,         &
    &                               sucldp, suphli,    &
    &                               suvdf , suvdfs
  USE mo_convect_tables,      ONLY: init_convect_tables
  ! EDMF DUAL turbulence
  USE mo_edmf_param,          ONLY: suct0, su0phy, susekf, susveg, sussoil
  ! turbulence
  USE mo_turbdiff_config,     ONLY: turbdiff_config
  USE mo_data_turbdiff,       ONLY: get_turbdiff_param, lsflcnd, &
                                    impl_s, impl_t,              &
                                    impl_weight
  USE src_turbdiff_new,       ONLY: organize_turbdiff
  USE src_turbdiff,           ONLY: turbtran, turbdiff

  ! air-sea-land interface
  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, iwtr, iice, ilnd !, &
  !   &                                init_sfc_indices
  ! vertical diffusion
  USE mo_echam_vdiff_params,  ONLY: init_vdiff_params, z0m_min, &
    &                                tke_min
  USE mo_vdiff_solver,        ONLY: init_vdiff_solver
  USE mo_nwp_sfc_utils,       ONLY: nwp_surface_init, init_snowtile_lists, init_sea_lists, &
    &                               aggregate_t_g_q_v
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_lnd, lsnowtile, ntiles_water, &
    &                               lseaice, isub_water, isub_lake, isub_seaice
  USE mo_phyparam_soil,       ONLY: csalbw!, z0_lu
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                                sat_pres_ice, &  !! saturation vapor pressure w.r.t. ice
    &                                spec_humi !,qsat_rho !! Specific humidity

  USE data_gwd,               ONLY: sugwwms

  USE mo_nh_testcases_nml,    ONLY: nh_test_name, ape_sst_case, th_cbl
  USE mo_nh_wk_exp,           ONLY: qv_max_wk
  USE mo_ape_params,          ONLY: ape_sst
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_nwp_parameters,      ONLY: t_phy_params

  USE mo_datetime,            ONLY: iso8601
  USE mo_time_config,         ONLY: time_config
  USE mo_initicon_config,     ONLY: init_mode
  USE mo_mcrph_sb,            ONLY: two_moment_mcrph_init

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: init_nwp_phy

CONTAINS

SUBROUTINE init_nwp_phy ( pdtime,                           &
                       &  p_patch, p_metrics,               &
                       &  p_prog_now,  p_diag,              &
                       &  prm_diag,prm_nwp_tend,            &
                       &  p_prog_lnd_now, p_prog_lnd_new,   &
                       &  p_prog_wtr_now, p_prog_wtr_new,   &
                       &  p_diag_lnd,                       &
                       &  ext_data, phy_params)

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

  INTEGER             :: jk, jk1
  REAL(wp)            :: rsltn   ! horizontal resolution
  REAL(wp)            :: pdtime
  REAL(wp)            :: pref(p_patch%nlev)
  REAL(wp)            :: zlat, zprat, zn1, zn2, zcdnc
  REAL(wp)            :: zpres
  REAL(wp)            :: gz0(nproma)

  CHARACTER(len=16)   :: cur_date     ! current date (iso-Format)
  INTEGER             :: icur_date    ! current date converted to integer

  ! Reference atmosphere parameters
  REAL(wp), PARAMETER :: htropo = 11000._wp       ! [m]    tropopause height
  REAL(wp), PARAMETER :: t00    = 288.15_wp       ! [m]    temperature at sea level
  REAL(wp) :: ttropo, ptropo, temp, zfull

  LOGICAL :: lland, lglac
  LOGICAL :: ltkeinp_loc, lgz0inp_loc  !< turbtran switches

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

  INTEGER :: khydromet, ktrac

  INTEGER :: k1500m                  ! index of first half level above 1500m
  INTEGER :: istatus=0

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
     routine = 'mo_nwp_phy_init:init_nwp_phy'

  i_nchdom  = MAX(1,p_patch%n_childdom)

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1
  jg     = p_patch%id

  nshift = p_patch%nshift_total

  i_lc_si= ext_data%atm%i_lc_snow_ice


  ! for both restart and non-restart runs. Could not be included into 
  ! mo_ext_data_state/init_index_lists due to its dependence on p_diag_lnd.
  CALL init_sea_lists(p_patch, ext_data, p_diag_lnd, lseaice)

  IF (.NOT. is_restart_run())THEN

    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &  i_startidx, i_endidx, rl_start, rl_end)

      IF (ltestcase .AND. (nh_test_name == 'APE_nh' .OR. nh_test_name == 'dcmip_tc_52') ) THEN

        ! t_g = ape_sst1
        
        DO jc = i_startidx, i_endidx
          zlat = p_patch%cells%center(jc,jb)%lat
          p_prog_lnd_now%t_g (jc,jb) = ape_sst(ape_sst_case,zlat) ! set SST
          p_prog_lnd_new%t_g (jc,jb) = ape_sst(ape_sst_case,zlat) ! set SST
          ! Humidity at water surface = humidity at saturation
          p_diag_lnd%qv_s    (jc,jb) = &
  !      & qsat_rho(p_prog_lnd_now%t_g (jc,jb),p_prog%rho(jc,nlev,jb))
        & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
        END DO

!        IF( atm_phy_nwp_config(jg)%inwp_radiation > 0 .AND. irad_o3 == io3_ape) THEN
!          DO jc = i_startidx, i_endidx
!            zf_aux( jc,1:nlev_o3,jb) = ext_data%atm_td%zf(1:nlev_o3)
!          ENDDO
!        END IF

      ELSE IF (ltestcase .AND. nh_test_name == 'wk82' ) THEN !
 
        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)*  &
                    ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
          p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb) 
         p_diag_lnd%qv_s     (jc,jb) = &
        & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))  
          p_diag_lnd%qv_s    (jc,jb) = MIN (p_diag_lnd%qv_s(jc,jb) ,   &
                                     &     p_prog_now%tracer(jc,nlev,jb,iqv)) 
        END DO

      ELSE IF (ltestcase .AND. nh_test_name == 'CBL' ) THEN !
 
        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g (jc,jb) = th_cbl(1)
          p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb) 
         p_diag_lnd%qv_s     (jc,jb) = &
        & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))  
        END DO

      ELSE IF (ltestcase) THEN ! any other testcase

        ! t_g  =  t(nlev)
        ! qv_ s= qv(nlev)
        ! KF increase the surface values to obtain fluxes          

        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)!+0.2_wp
          p_prog_lnd_new%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)!+0.2_wp
          ! KF NOTE: as long as we have only water as lower boundary
          ! this is the same setting as for APE
         p_diag_lnd%qv_s    (jc,jb) = &
!                & qsat_rho(p_prog_lnd_now%t_g (jc,jb),p_prog%rho(jc,nlev,jb))
        & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))

        END DO
      ELSE ! For real-case simulations, initialize also qv_s and the tile-based fields

         ! t_g:
         ! Note, that in copy_prepicon2prog the entire t_g field is initialized with 
         ! t_skin.
         ! Here, t_g is re-initialized over open water points with t_seasfc.
         ! Thus:
         ! t_g = tskin (from IFS), for land and seaice points
         ! t_g = t_seasfc for open water and lake points
         !
         ! If l_sst_in==FALSE, then t_seasfc=t_skin (with a limiter), so nothing important happens
         !
         ! qv_s:
         ! Over the sea and over the ice, qv_s is set to the saturated value
         ! Over the land we take the minimum of the saturated value and the value 
         ! at the first main level above ground
         !

         !t_g_t and qv_s_t are initialized in read_dwdfg_sfc, calculate the aggregated values 
         ! needed for example for initializing the turbulence fields
         IF (init_mode /= MODE_IFSANA) THEN
          CALL aggregate_t_g_q_v( p_patch, ext_data, p_prog_lnd_now , &
          &                           p_diag_lnd )    
          DO jc = i_startidx, i_endidx
           p_prog_lnd_new%t_g(jc,jb)     =  p_prog_lnd_now%t_g(jc,jb)
          ENDDO

          DO jt = 1, ntiles_total+ntiles_water          
           DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g_t(jc,jb,jt)
           END DO
          END DO

         END IF

         ! t_g_t  qv_s and qv_s_t are not initialized in case of MODE_IFSANA
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
           p_prog_lnd_now%t_g(jc,jb) = p_diag_lnd%t_seasfc(jc,jb)
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
           p_prog_lnd_new%t_g(jc,jb)     =  p_prog_lnd_now%t_g(jc,jb)
          ENDDO


          DO jt = 1, ntiles_total+ntiles_water
          
           DO jc = i_startidx, i_endidx
          
            p_prog_lnd_now%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
            p_diag_lnd%qv_s_t(jc,jb,jt) = p_diag_lnd%qv_s(jc,jb)
           ENDDO
         ENDDO
        END IF
      ENDIF

      ! Copy t_g to t_seasfc for idealized cases with surface scheme (would be undefined otherwise)
      IF (ltestcase .AND. atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        DO jc = i_startidx, i_endidx
          p_diag_lnd%t_seasfc(jc,jb) = p_prog_lnd_now%t_g(jc,jb)
        ENDDO
      ENDIF
        
    END DO
    CALL message('mo_nwp_phy_init:', 'initialized surface temp and humidity')

  ELSE  ! if is_restart_run()
    !
    ! necessary, because only t_g(nnow_rcf) is written to the restart file
    ! with the following copy statement the ocean points of t_g(nnew_rcf) are 
    ! filled with the correct values.
    !
    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &  i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
      ENDDO
      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN ! the t_g_t does not exist for inwp_surface=0
        DO jt = 1, ntiles_total+ntiles_water
          DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g_t(jc,jb,jt)
          ENDDO            
        ENDDO
      ENDIF

    ENDDO

  END IF


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


  !------------------------------------------
  !< call for cloud microphysics
  !------------------------------------------

  SELECT CASE ( atm_phy_nwp_config(jg)%inwp_gscp )
    
  CASE (1)  ! new hydci_pp from COSMO_V4_23
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics')
    CALL hydci_pp_init
    
  CASE (3)  ! ice nucleation scheme by C.K. based on new hydci_pp from COSMO_V4_23
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics')
    CALL hydci_pp_ice_init

  CASE (4) !two moment micrphysics
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysic:SB')
    CALL two_moment_mcrph_init( )

  CASE (10)  ! old hydci_pp from COSMO_V4_14
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics')
    CALL hydci_pp_old_init          
    
  END SELECT



  !------------------------------------------
  !< radiation
  !------------------------------------------
  IF ( atm_phy_nwp_config(jg)%inwp_radiation == 1 ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init RRTM')

    SELECT CASE ( irad_aero )
    ! Note (GZ): irad_aero=2 does no action but is the default in radiation_nml
    ! and therefore should not cause the model to stop
    CASE (0,2,5,6)
      !ok
    CASE DEFAULT
      CALL finish('mo_nwp_phy_init: init_nwp_phy',  &
        &      'Wrong irad_aero. For RRTM radiation, this irad_aero is not implemented.')
    END SELECT
    
!    prm_diag%lfglac (:,:) = ext_data%atm%soiltyp(:,:) == 1  !soiltyp=ice

    ! solar flux (W/m2) in 14 SW bands
    ssi(:) = ssi_amip(:)
    ! solar constant (W/m2)
    tsi    = SUM(ssi(:))


    !------------------------------------------
    !< set conditions for Aqua planet experiment  
    !------------------------------------------
    IF ( nh_test_name == 'APE_nh' .OR. nh_test_name == 'dcmip_tc_52' ) THEN
      ssi(:) = ssi(:)*1365._wp/tsi
      tsi = 1365._wp
    ENDIF  ! APE

    
    CALL setup_srtm

    CALL lrtm_setup(lrtm_filename)

    CALL setup_newcld_optics(cldopt_filename)
    
    rl_start = 1  ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk,zprat,zpres,lland,lglac,zn1,&
!$OMP zn2,zcdnc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &  i_startidx, i_endidx, rl_start, rl_end)

      ! Initialize cloud droplet number concentration (acdnc)
      ! like in mo_echam_phy_init.f90
      DO jk = 1,nlev
        ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
        DO jc = 1, i_endidx
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
        END DO !jc
      END DO   !jk


    !------------------------------------------
    ! APE ozone profile, vertical setting needed only once for NH
    !------------------------------------------
      IF (irad_o3 == io3_ape .AND. .NOT. is_restart_run()) THEN

!        CALL o3_zl2ml(p_patch%nblks_c,p_patch%npromz_c,        & ! 
!          &           nlev_o3,      nlev,                      & ! vertical levels in/out
!          &           zf_aux,   p_metrics%z_mc,                & ! vertical in/out
!          &           ext_data%atm_td%o3(:,:,:,nmonths),p_prog%tracer(:,:,:,io3))! o3Field in/out
 
        CALL o3_pl2ml ( kproma= i_endidx, kbdim=nproma,  &
          & nlev_pres = nlev_o3,klev= nlev ,             &
          & pfoz = ext_data%atm_td%pfoz(:),              &
          & phoz = ext_data%atm_td%phoz(:),              &! in o3-levs
          & ppf = p_diag%pres (:,:,jb),                  &! in  pres
          & pph = p_diag%pres_ifc(:,:,jb),               &! in  pres_halfl
          & o3_time_int = ext_data%atm_td%o3(:,:,jb,nmonths),     &! in
          & o3_clim     = ext_data%atm%o3(:,:,jb) )         ! OUT 
        
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
      
    ELSEIF ( irad_aero == 6 ) THEN
      
      CALL init_aerosol_props_tegen_rrtm

    ELSE
      
      zaea_rrtm(:,:) = 0.0_wp
      zaes_rrtm(:,:) = 0.0_wp
      zaeg_rrtm(:,:) = 0.0_wp
      
    ENDIF

    DO ist = 1, UBOUND(csalbw,1)
      rad_csalbw(ist) = csalbw(ist) / (2.0_wp * zml_soil(1))
    ENDDO
    
  ELSEIF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init Ritter Geleyn')

    ! Note (GZ): irad_aero=2 does no action but is the default in radiation_nml
    ! and therefore should not cause the model to stop
    SELECT CASE ( irad_aero )
    CASE (0,2,5,6)
      !ok
    CASE DEFAULT
      CALL finish('mo_nwp_phy_init: init_nwp_phy',  &
        &      'Wrong irad_aero. For Ritter-Geleyn radiation, this irad_aero is not implemented.')
    END SELECT

    ! solar flux (W/m2) in 14 SW bands    
    ssi(:) = ssi_amip(:)
    ! solar constant (W/m2)
    tsi    = SUM(ssi(:))


    !------------------------------------------
    !< set conditions for Aqua planet experiment  
    !------------------------------------------
    IF ( nh_test_name == 'APE_nh' .OR. nh_test_name == 'dcmip_tc_52' ) THEN
      ssi(:) = ssi(:)*1365._wp/tsi
      tsi = 1365._wp
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
      
    ELSEIF ( irad_aero == 6 ) THEN

      CALL init_aerosol_props_tegen_rg
        
    ELSE

      zaea_rg(:,:) = 0.0_wp
      zaes_rg(:,:) = 0.0_wp
      zaeg_rg(:,:) = 0.0_wp

    ENDIF 


    !------------------------------------------
    ! APE ozone profile, vertical setting needed only once for NH
    !------------------------------------------
    IF (irad_o3 == io3_ape .AND. .NOT. is_restart_run() ) THEN

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

  ENDIF !inwp_radiation


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
    cur_date = iso8601(time_config%cur_datetime)
    ! convert first 8 characters to integer (yyyymmdd)
    READ(cur_date(1:8),'(i8)') icur_date

    CALL sucst(54,icur_date,0,0)
    CALL su_yoethf

  ENDIF


  !------------------------------------------
  !< call for convection
  !------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 .OR. &
    &  atm_phy_nwp_config(jg)%inwp_turb == iedmf )     THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init convection')

    ! Please take care for scale-dependent initializations!
    ! Spectral resolution corresponding to ICON
    ! needed for RTAU - CAPE calculation
    ! adapted for more general gemoetries
    rsltn = p_patch%geometry_info%mean_characteristic_length 
    

!    WRITE(message_text,'(i3,i10,f20.10)') jg, nsmax, phy_params%mean_charlen
!    CALL message('nwp_phy_init, nsmax=', TRIM(message_text))


    CALL sucumf(rsltn,nlev,pref,phy_params)
    CALL suphli
    CALL suvdf
    CALL suvdfs
    CALL sucldp
    CALL message('mo_nwp_phy_init:', 'convection initialized')
  ENDIF


  !------------------------------------------
  !< surface initialization (including seaice)
  !------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. .NOT. is_restart_run() ) THEN  ! TERRA
    CALL nwp_surface_init(p_patch, ext_data, p_prog_lnd_now, p_prog_lnd_new, &
      &                   p_prog_wtr_now, p_prog_wtr_new, p_diag_lnd, p_diag)
  ELSE IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. is_restart_run()) THEN

    IF ( lsnowtile ) THEN
      CALL init_snowtile_lists(p_patch, ext_data, p_diag_lnd)
    ENDIF
  END IF


  !------------------------------------------
  !< setup for turbulence
  !------------------------------------------

  ! initialize gz0 (roughness length * g)
  !
  IF ( ANY( (/icosmo,igme,10,11,12/)==atm_phy_nwp_config(jg)%inwp_turb ) .AND. .NOT. is_restart_run() ) THEN


    ! gz0 is initialized, if we start from IFS surface (MODE_IFSANA) or 
    ! GME surface (MODE_COMBINED, MODE_COSMODE)
    IF (ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMODE/) == init_mode) ) THEN
 
      IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init roughness length')

      IF (turbdiff_config(jg)%lconst_z0) THEN
        ! constant z0 for idealized tests
        prm_diag%gz0(:,:) = grav * turbdiff_config(jg)%const_z0

      ELSE IF (atm_phy_nwp_config(jg)%itype_z0 == 1) THEN
        ! default
        prm_diag%gz0(:,:) = grav * ext_data%atm%z0(:,:)

      ELSE IF (atm_phy_nwp_config(jg)%itype_z0 == 2) THEN

        rl_start = grf_bdywidth_c + 1 ! land-cover classes are not set for nest-boundary points
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

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

    ENDIF  !operation mode

  ! For 3D Smagorinsky turbulence model 
  ELSE IF (  atm_phy_nwp_config(jg)%is_les_phy .AND. .NOT. is_restart_run() ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init Smagorinsky turbulence')

    IF (turbdiff_config(jg)%lconst_z0) THEN
      ! for idealized tests
      prm_diag%gz0(:,:) = grav * turbdiff_config(jg)%const_z0
      IF(turbdiff_config(jg)%const_z0==0._wp) &
       CALL finish (TRIM(routine), 'roughness length needs to be set for idealized LES cases!')
    ELSE 
      ! default: these are all set in nwp_turbtrans (AD: 11.09.2013)
    ENDIF

  ENDIF

  IF ( ANY( (/icosmo,10,11,12/)==atm_phy_nwp_config(jg)%inwp_turb )) THEN
  
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
  IF ( ANY( (/icosmo,10,11,12/)==atm_phy_nwp_config(jg)%inwp_turb ) .AND. .NOT. is_restart_run() ) THEN
  
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init COSMO turbulence')

    rl_start = 1 ! Initialization is done also for nest boundary points
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,ic,jc,jt, &
!$OMP            ltkeinp_loc,lgz0inp_loc,nlevcm) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      
      IF (ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMODE/) == init_mode) ) THEN

        ltkeinp_loc = .FALSE.  ! do not re-initialize TKE field
        lgz0inp_loc = .FALSE.  ! do re-initialize gz0 field (water points only)

      ELSE  ! init_mode=MODE_DWDANA
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

      IF ( ANY( (/10,11/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN

        nlevcm = nlevp1

        CALL organize_turbdiff( lstfnct=.TRUE., &
          &  lturatm=.FALSE., ltursrf=.TRUE., iini=1, &
          &  ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, &
          &  lmomdif=.FALSE., lscadif=.FALSE., itnd=0, &
          &  dt_var=pdtime, dt_tke=pdtime, &
          &  nprv=1, ntur=1, ntim=1, &
          &  ie=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm, &
          &  i_st=i_startidx, i_en=i_endidx, i_stp=i_startidx, i_enp=i_endidx, &
          &  l_hori=phy_params%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb), &
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
          &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
          &  tke=p_prog_now%tke(:,:,jb), &
          &  tkvm=prm_diag%tkvm(:,2:nlevp1,jb), tkvh=prm_diag%tkvh(:,2:nlevp1,jb), &
          &  rcld=prm_diag%rcld(:,:,jb), &
          &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), &
          &  td_2m=prm_diag%td_2m(:,jb), rh_2m=prm_diag%rh_2m(:,jb), &
          &  u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb), &
          &  shfl_s=prm_diag%shfl_s(:,jb), qvfl_s=prm_diag%qhfl_s(:,jb), &
          &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

        prm_diag%lhfl_s(i_startidx:i_endidx,jb) = &
          &  prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v
        
      ELSE

        CALL turbtran(iini=1, ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, &
!
           &  dt_tke=pdtime, nprv=1, ntur=1, ntim=1,                    &
!
           &  ie=nproma, ke=nlev, ke1=nlevp1,                                           &
           &  istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx, &
!
           &  l_hori=phy_params%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),                &
!
           &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
           &  sai=ext_data%atm%sai(:,jb), h_ice=p_prog_wtr_now%h_ice (:,jb), &
!
           &  ps=p_diag%pres_sfc(:,jb), t_g=p_prog_lnd_now%t_g(:,jb), qv_s=p_diag_lnd%qv_s(:,jb), &
!
           &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), T=p_diag%temp(:,:,jb),   &
           &  qv=p_prog_now%tracer(:,:,jb,iqv), qc=p_prog_now%tracer(:,:,jb,iqc), &
!
           &  prs=p_diag%pres(:,:,jb),  &
!
           &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
           &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
!
           &  tke=p_prog_now%tke(:,:,jb), & !  edr=prm_diag%edr(:,:,jb), &
           &  tkvm=prm_diag%tkvm(:,2:nlevp1,jb), tkvh=prm_diag%tkvh (:,2:nlevp1,jb),  &
           &  rcld=prm_diag%rcld(:,:,jb),                                             &
!
           &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), td_2m=prm_diag%td_2m (:,jb), &
           &  rh_2m=prm_diag%rh_2m(:,jb), u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m (:,jb), &
           &  shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb), qhfl_s=prm_diag%qhfl_s(:,jb),&
           &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

      END IF

      IF ( ANY( (/10,12/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN

        nlevcm = nlevp1

        CALL organize_turbdiff( lstfnct=.TRUE., lsfluse=lsflcnd, &
          &  lturatm=.TRUE., ltursrf=.FALSE., iini=1, &
          &  ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, &
          &  lmomdif=.TRUE., lscadif=.TRUE., itnd=0, &
          &  dt_var=pdtime, dt_tke=pdtime, &
          &  nprv=1, ntur=1, ntim=1, &
          &  ie=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm, &
          &  i_st=i_startidx, i_en=i_endidx, i_stp=i_startidx, i_enp=i_endidx, &
          &  l_hori=phy_params%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb), &
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
          &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
          &  tke=p_prog_now%tke(:,:,jb), &
          &  tkvm=prm_diag%tkvm(:,2:nlevp1,jb), tkvh=prm_diag%tkvh(:,2:nlevp1,jb), &
          &  rcld=prm_diag%rcld(:,:,jb), &
          &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), &
          &  v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb), &
          &  tketens=prm_nwp_tend%ddt_tke(:,:,jb), &
          &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb), vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb), &
          &  shfl_s=prm_diag%shfl_s(:,jb), qvfl_s=prm_diag%qhfl_s(:,jb), &
          &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

        IF ( .NOT. lsflcnd ) THEN
          prm_diag%lhfl_s(i_startidx:i_endidx,jb) = &
            &  prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v
        END IF
        
      ELSE

        CALL turbdiff(iini=1, ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, lstfnct=.TRUE., &
!
           &  dt_var=pdtime, dt_tke=pdtime, nprv=1, ntur=1, ntim=1, &
!
           &  ie=nproma, ke=nlev, ke1=nlevp1, kcm=nlevp1, &
           &  istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx,  &
!
           &  l_hori=phy_params%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),                &
           &  dp0=p_diag%dpres_mc(:,:,jb),                                                &
!
           &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
           &  h_ice=p_prog_wtr_now%h_ice (:,jb),                                        &
!
           &  ps=p_diag%pres_sfc(:,jb), t_g=p_prog_lnd_now%t_g(:,jb), qv_s=p_diag_lnd%qv_s(:,jb), &
!
           &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog_now%w(:,:,jb), T=p_diag%temp(:,:,jb), &
           &  qv=p_prog_now%tracer(:,:,jb,iqv), qc=p_prog_now%tracer(:,:,jb,iqc), &
!
           &  prs=p_diag%pres(:,:,jb), rho=p_prog_now%rho(:,:,jb), epr=p_prog_now%exner(:,:,jb), &
!
           &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
           &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
!
           &  tke=p_prog_now%tke(:,:,jb), & !  edr=prm_diag%edr(:,:,jb), &
           &  tkvm=prm_diag%tkvm(:,2:nlevp1,jb), tkvh=prm_diag%tkvh (:,2:nlevp1,jb),  &
           &  rcld=prm_diag%rcld(:,:,jb),                                             &
!
           &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb), &
           &  tketens=prm_nwp_tend%ddt_tke(:,:,jb), &
           &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb), vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) ,&
!
           &  shfl_s=prm_diag%shfl_s(:,jb), qhfl_s=prm_diag%qhfl_s(:,jb), &
!
           &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

      END IF

      ! tile-specific quantities needed by turbtran
      ! 
      DO jt = 1, ntiles_total+ntiles_water 
        prm_diag%gz0_t   (:,jb,jt) = prm_diag%gz0(:,jb)
        prm_diag%tvs_s_t (:,jb,jt) = p_prog_now%tke(:,nlevp1,jb)  !here: SQRT(2*TKE) 
        prm_diag%tkvm_s_t(:,jb,jt) = prm_diag%tkvm(:,nlevp1,jb)
        prm_diag%tkvh_s_t(:,jb,jt) = prm_diag%tkvh(:,nlevp1,jb)
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


  ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN  !EDMF DUALM
    CALL suct0
    CALL su0phy
    CALL susekf
    CALL susveg
    CALL sussoil

!$OMP PARALLEL WORKSHARE
    ! paranoia: Make sure that rcld is initialized  (needed by cloud cover scheme)
    prm_diag%rcld(:,:,:)    = 0._wp
!$OMP END PARALLEL WORKSHARE

  ENDIF


  ! Gravity wave drag scheme
  !
  IF ( atm_phy_nwp_config(jg)%inwp_gwd == 1 ) THEN  ! IFS gwd scheme

    CALL sugwwms(nflevg=nlev, ppref=pref, klaunch=phy_params%klaunch)
    CALL message('mo_nwp_phy_init:', 'non-orog GWs initialized')

  END IF


END SUBROUTINE init_nwp_phy


END MODULE mo_nwp_phy_init

