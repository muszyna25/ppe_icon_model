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
MODULE mo_nwp_phy_init

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: re, grav, rd_o_cpd, cpd, p0ref, rd, p0sl_bg
  USE mo_math_utilities,      ONLY: mean_domain_values
  USE mo_model_domain_import, ONLY: nroot   
  USE mo_nwp_phy_state,       ONLY: t_nwp_phy_diag,t_nwp_phy_tend
  USE mo_nwp_lnd_state,       ONLY: t_lnd_prog, t_lnd_diag  !, t_tiles
  USE mo_ext_data,            ONLY: t_external_data, nlev_o3, nmonths
  USE mo_nonhydro_state,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_exception,           ONLY: message, finish,message_text
  USE mo_vertical_coord_table,ONLY: vct_a, vct
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell, zml_soil,io3_ape
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltestcase, iqv, iqc, io3, msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  !radiation
  USE mo_newcld_optics,       ONLY: setup_newcld_optics
  USE mo_lrtm_setup,          ONLY: lrtm_setup
  USE mo_radiation_config,    ONLY: ssi, tsi,irad_o3, irad_aero, rad_csalbw 
  USE mo_srtm_config,         ONLY: setup_srtm, ssi_amip
  USE mo_radiation_rg_par,    ONLY: rad_aibi
  USE mo_aerosol_util,        ONLY: init_aerosol_dstrb_tanre,      &
    &                               init_aerosol_props_tanre_rg,   &
    &                               init_aerosol_props_tanre_rrtm, &
    &                               init_aerosol_props_tegen_rg,   &
    &                               init_aerosol_props_tegen_rrtm, &
    &                               zaef_rg, zaea_rg, zaes_rg, zaeg_rg, &
    &                               zaea_rrtm, zaes_rrtm, zaeg_rrtm
  USE mo_o3_util,             ONLY: o3_pl2ml!, o3_zl2ml

  ! microphysics
  USE mo_gscp_cosmo,          ONLY: hydci_pp_init
  ! convection
  USE mo_cuparameters,        ONLY: sucst,  sucumf,    &
    &                               su_yoethf,         &
    &                               sucldp, suphli,    &
    &                               suvdf , suvdfs
  USE mo_convect_tables,      ONLY: init_convect_tables
  ! EDMF DUAL turbulence
  USE mo_edmf_param,          ONLY: suct0, su0phy, susekf
  ! turbulence
  USE mo_turbdiff_config,     ONLY: turbdiff_config
  USE mo_data_turbdiff,       ONLY: get_turbdiff_param
  USE src_turbdiff,           ONLY: init_canopy, organize_turbdiff
  ! for APE_nh experiments

  ! air-sea-land interface
  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, iwtr, iice, ilnd !, &
  !   &                                init_sfc_indices
  ! vertical diffusion
  USE mo_echam_vdiff_params,  ONLY: init_vdiff_params, z0m_min, &
    &                                tke_min
  USE mo_vdiff_solver,        ONLY: init_vdiff_solver
  USE mo_nwp_sfc_interface,   ONLY: nwp_surface_init
  USE mo_phyparam_soil,       ONLY: csalbw
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                                spec_humi !,qsat_rho !! Specific humidity

  USE data_gwd,               ONLY: sugwwms

  USE mo_nh_testcases,        ONLY: nh_test_name, ape_sst_case
  USE mo_nh_wk_exp,           ONLY: qv_max_wk
  USE mo_ape_params,          ONLY: ape_sst
  USE mo_master_control,      ONLY: is_restart_run

  IMPLICIT NONE

  PUBLIC  :: init_nwp_phy

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

SUBROUTINE init_nwp_phy ( pdtime                         , &
                       &  p_patch, p_metrics,              &
                       &  p_prog_now,  p_prog,  p_diag,    &
                       &  prm_diag,prm_nwp_tend,           &
                       &  p_prog_lnd_now, p_prog_lnd_new,  &
                       &  p_diag_lnd,                      &
                       &  ext_data, mean_charlen )

  TYPE(t_patch),        TARGET,INTENT(in)    :: p_patch
  TYPE(t_nh_metrics),          INTENT(in)    :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout) :: p_prog_now !!the prognostic variables
  TYPE(t_nh_prog),      TARGET,INTENT(inout) :: p_prog  !!the prognostic variables
  TYPE(t_nh_diag),      TARGET,INTENT(inout) :: p_diag  !!the diagostic variables
  TYPE(t_external_data),       INTENT(in)    :: ext_data
  TYPE(t_nwp_phy_diag),        INTENT(inout) :: prm_diag
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout) :: prm_nwp_tend
  TYPE(t_lnd_prog),            INTENT(inout) :: p_prog_lnd_now, p_prog_lnd_new
  TYPE(t_lnd_diag),            INTENT(inout) :: p_diag_lnd

  REAL(wp),INTENT(OUT)::  mean_charlen
  INTEGER             :: jk, jk1, nsmax
  REAL(wp)            :: pdtime
  REAL(wp)            :: pref(p_patch%nlevp1)
  REAL(wp), PARAMETER :: h_scal = 8000._wp     ! [m]      scale height
  REAL(wp)            :: zlat, zprat, zn1, zn2, zcdnc
!  REAL(wp)            :: zf_aux(nproma,nlev_o3,p_patch%nblks_c)
  REAL(wp)            :: zpres

  LOGICAL  :: lland, lglac
  
  INTEGER :: jb,jc,jg,ist
  INTEGER :: nlev, nlevp1            !< number of full and half levels
  INTEGER :: nshift                  !< shift with respect to global grid
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !! slices
  INTEGER :: i_nchdom                !! domain index
!  INTEGER :: inwp_turb_init          !< 1: initialize nwp_turb
!                                     !< 0: do not initialize

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''

  INTEGER :: khydromet, ktrac

    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1
    jg     = p_patch%id

    nshift = p_patch%nshift_total

    IF (.NOT. is_restart_run())THEN

    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
             &  i_startidx, i_endidx, rl_start, rl_end)

        IF(ltestcase .AND. nh_test_name == 'APE_nh') THEN

          ! t_g = ape_sst1
          
          DO jc = i_startidx, i_endidx
            zlat = p_patch%cells%center(jc,jb)%lat
            p_prog_lnd_now%t_g (jc,jb) = ape_sst(ape_sst_case,zlat) ! set SST
            p_prog_lnd_new%t_g (jc,jb) = ape_sst(ape_sst_case,zlat) ! set SST
            ! Humidity at water surface = humidity at saturation
            p_diag_lnd%qv_s    (jc,jb) = &
  !        & qsat_rho(p_prog_lnd_now%t_g (jc,jb),p_prog%rho(jc,nlev,jb))
          & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

!!          IF( atm_phy_nwp_config(jg)%inwp_radiation > 0 .AND. irad_o3 == io3_ape) THEN
!            DO jc = i_startidx, i_endidx
!              zf_aux( jc,1:nlev_o3,jb) = ext_data%atm_td%zf(1:nlev_o3)
!            ENDDO
!          END IF

        ELSE IF (ltestcase .AND. nh_test_name == 'wk82' ) THEN !
 
          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)*  &
                      ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
            p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb) 
           p_diag_lnd%qv_s     (jc,jb) = &
          & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))  
            p_diag_lnd%qv_s    (jc,jb) = MIN (p_diag_lnd%qv_s(jc,jb) ,   &
                                       &     p_prog%tracer(jc,nlev,jb,iqv)) 
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
        ELSE ! For real-case simulations, initialize also qv_s
          DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g(jc,jb) =  p_prog_lnd_now%t_g(jc,jb)
            p_diag_lnd%qv_s    (jc,jb) = &
            & spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
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
      ENDDO

    END IF


    !--------------------------------------------------------------
    !< characteristic gridlength needed by convection and turbulence
    !--------------------------------------------------------------
      CALL mean_domain_values (p_patch%level, nroot, mean_charlen)

    !--------------------------------------------------------------
    !>reference pressure
    !--------------------------------------------------------------
      DO jk = nlevp1, 1, -1
        jk1 = jk + nshift
        pref(jk)= p0sl_bg * EXP( -vct_a (jk1)/h_scal)
      ENDDO

    !------------------------------------------
    !< call for cloud microphysics
    !------------------------------------------
  IF (  atm_phy_nwp_config(jg)%inwp_gscp == 1 )  THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics')
    CALL hydci_pp_init

  ENDIF
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

    ssi(:) = ssi_amip(:)
    tsi    = SUM(ssi(:))


    !------------------------------------------
    !< set conditions for Aqua planet experiment  
    !------------------------------------------
    IF ( nh_test_name == 'APE_nh' ) THEN
      ssi(:) = ssi(:)*1365._wp/tsi
      tsi = 1365._wp
    ENDIF  ! APE

    
    CALL setup_srtm

    CALL lrtm_setup

    CALL setup_newcld_optics
    
    rl_start = 1  ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk,zprat,zpres,lland,lglac,zn1,zn2,zcdnc)
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
          & o3_clim     = p_prog%tracer(:,:,jb,io3) )     ! OUT
        
      ENDIF


    ENDDO      !jb
!$OMP END DO
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

    DO ist = 1, 10
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

    
    ssi(:) = ssi_amip(:)
    tsi    = SUM(ssi(:))

    !------------------------------------------
    !< set conditions for Aqua planet experiment  
    !------------------------------------------
    IF ( nh_test_name == 'APE_nh' ) THEN
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
!$OMP DO PRIVATE(jb,i_endidx)
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
          & o3_clim     = p_prog%tracer(:,:,jb,io3) )     ! OUT

      ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    ENDIF ! (irad_o3 == io3_ape)
 
    DO ist = 1, 10
      rad_csalbw(ist) = csalbw(ist) / (2.0_wp * zml_soil(1))
    ENDDO

  ENDIF !inwp_radiation

  !----------------------------------------------------------------
  !< initializations needed both for convection and inwp_cldcover=1 
  !----------------------------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 .OR. &
    &  atm_phy_nwp_config(jg)%inwp_cldcover == 1   .OR. &
    &  atm_phy_nwp_config(jg)%inwp_turb == 3 )     THEN
    
    !This has to be done here because not only convection, but also inwp_cldcover == 1 
    !uses mo_cufunctions's foealfa. Therefore, the parameters of the function foealfa
    !have to be initialized by calls of sucst and su_yoethf.
    
    CALL sucst(54,20020211,0,0)
    CALL su_yoethf

  ENDIF

  !------------------------------------------
  !< call for convection
  !------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 .OR. &
    &  atm_phy_nwp_config(jg)%inwp_turb == 3 )     THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init convection')
! Please take care for scale-dependent initializations!



      nsmax = INT(2._wp*pi*re/mean_charlen)

!        WRITE(message_text,'(i3,i10,f20.10)') jg,nsmax,mean_charlen
!       CALL message('nwp_phy_init, nsmax=', TRIM(message_text))

        CALL sucumf(nsmax,nlevp1,pref)
        CALL suphli
        CALL suvdf
        CALL suvdfs
        CALL sucldp
        CALL message('mo_nwp_phy_init:', 'convection initialized')
  ENDIF


  !------------------------------------------
  !< setup for turbulence
  !------------------------------------------

  ! nsfc_type is used for dimensioning local variables in the NWP interface;
  ! thus, it must be set even if vdiff is not called
 ! IF ( atm_phy_nwp_config(jg)%inwp_turb /= 2) &
    nsfc_type = 1 ! for the time being, nsfc_type must be reset to 1 because land
                  ! is not yet avaliable for vdiff

  IF (  atm_phy_nwp_config(jg)%inwp_turb == 1 .AND. .NOT. is_restart_run() ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init COSMO turbulence')


    ! initialize gz0 (roughness length * g)
    !
    IF (turbdiff_config(jg)%lconst_z0) THEN
      ! for idealized tests
      prm_diag%gz0(:,:) = grav * turbdiff_config(jg)%const_z0
    ELSE 
      ! default
      prm_diag%gz0(:,:) = grav * ext_data%atm%z0(:,:)
    ENDIF

    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    CALL get_turbdiff_param(jg)

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

      CALL init_canopy( ie=nproma, je=1, ke=nlev, ke1=nlevp1, kcm=nlevp1, &
         &  istartpar=i_startidx, iendpar=i_endidx, jstartpar=1, jendpar=1, &
!
!        &  hhl=p_metrics%z_ifc(:,:,jb), &
         &  fr_land=ext_data%atm%fr_land(:,jb), plcov=ext_data%atm%plcov_mx(:,jb), & 
         &  sai=prm_diag%sai(:,jb), lai=ext_data%atm%lai_mx(:,jb), &
         &  tai=prm_diag%tai(:,jb), eai=prm_diag%eai(:,jb) )


      CALL organize_turbdiff(action='tran_diff', iini=1, lstfnct=.TRUE., &
!
         &  dt_var=pdtime, dt_tke=pdtime, nprv=1, ntur=1, ntim=1, &
!
         &  ie=nproma, je=1, ke=nlev, ke1=nlevp1, kcm=nlevp1, vst=0, &
!
         &  istart   =i_startidx, iend   =i_endidx, istartu=i_startidx, iendu=i_endidx, &
         &  istartpar=i_startidx, iendpar=i_endidx, istartv=i_startidx, iendv=i_endidx, &
!
         &  jstart   =1,          jend   =1       , jstartu=1         , jendu=1       , &
         &  jstartpar=1         , jendpar=1       , jstartv=1         , jendv=1       , &
!
         &  l_hori=mean_charlen, hhl=p_metrics%z_ifc(:,:,jb), dp0=p_diag%dpres_mc(:,:,jb), &
!
         &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
         &  sai=prm_diag%sai(:,jb), h_ice=prm_diag%h_ice (:,jb), &
!
         &  ps=p_diag%pres_sfc(:,jb), t_g=p_prog_lnd_now%t_g(:,jb), qv_s=p_diag_lnd%qv_s(:,jb), &
!
         &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb), T=p_diag%temp(:,:,jb), &
         &  qv=p_prog%tracer(:,:,jb,iqv), qc=p_prog%tracer(:,:,jb,iqc), &
!
         &  prs=p_diag%pres(:,:,jb), rho=p_prog%rho(:,:,jb), epr=p_prog%exner(:,:,jb), &
!
         &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
         &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
!
         &  tke=p_prog_now%tke(:,:,jb), & !  edr=prm_diag%edr(:,:,jb), &
         &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh (:,:,jb), rcld=prm_diag%rcld(:,:,jb), &
!
         &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb), &
         &  tketens=prm_nwp_tend%ddt_tke(:,:,jb), &
         &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb), vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) ,&
!
         &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), td_2m=prm_diag%td_2m (:,jb), &
         &  rh_2m=prm_diag%rh_2m(:,jb), u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m (:,jb), &
         &  shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb), &
!
         &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

    ENDDO
!$OMP END DO

!$OMP PARALLEL WORKSHARE
        p_prog %tke (:,:,:) =  p_prog_now%tke (:,:,:)
!$OMP END PARALLEL WORKSHARE

!$OMP END PARALLEL

    CALL message('mo_nwp_phy_init:', 'cosmo turbulence initialized')

  ELSE IF (  atm_phy_nwp_config(jg)%inwp_turb == 2) THEN  !ECHAM vdiff

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init ECHAM turbulence')
      ! Currently the tracer indices are sorted such that we count
      ! the water substances first, and then other species like 
      ! aerosols and their precursors. "ntracer" is the total number 
      ! of tracers (including water substances) handled in the model;
      ! "iqt" is the starting index for non-water species.
      ! Before more sophisticated meta-data structure becomes available, 
      ! it is assumed here that all tracers are subject to turbulent mixing.

    ! For surface processes: 
    ! nsfc_type, iwtr, etc. are set in this subroutine. 
    ! See mo_icoham_sfc_indicies.f90 for further details.

      !<KF temporarly set in, has to moved to general place
    CALL init_convect_tables

!   CALL init_sfc_indices( ltestcase, 'APE' ) !call of a hydrostatic testcase
                                              ! to obtain the demanded parameters

    khydromet = 2 !iqt - 1        ! # of hydrometeors
    ktrac = 1   !ntracer - iqt + 1  ! # of non-water species 

     !IF (p_patch%id == 1) CALL init_vdiff_solver( khydromet, ktrac, nproma, nlev, nsfc_type )
    IF (p_patch%id == 1) CALL init_vdiff_solver( khydromet, ktrac, nlev )

    CALL init_vdiff_params( nlev, nlevp1, nlevp1, vct )

      !KF special setting for ICONAM
    tke_min = 1.e-4_wp
        
!$OMP PARALLEL
!$OMP PARALLEL WORKSHARE
    prm_diag% ustar (:,:)   = 1._wp
    prm_diag% kedisp(:,:)   = 0._wp
    prm_diag% thvvar(:,:,:) = 1.e-4_wp
!$OMP END PARALLEL WORKSHARE
!$OMP END PARALLEL

    IF (iwtr<=nsfc_type) prm_diag%z0m_tile(:,:,iwtr) = 1.e-3_wp !see init_surf in echam (or z0m_oce?)
    IF (iice<=nsfc_type) prm_diag%z0m_tile(:,:,iice) = 1.e-3_wp !see init_surf in echam (or z0m_ice?)
!    IF (ilnd<=nsfc_type) prm_diag%z0m_tile(:,:,ilnd) = z0m_min  ! or maybe a larger value?
    IF (ilnd<=nsfc_type) prm_diag%z0m_tile(:,:,ilnd) = ext_data%atm%z0(:,:)
!    ENDIF
        
    WRITE(message_text,'(a,3I4)') 'init sfc inidces = ',iwtr,iice,nsfc_type
    CALL message('', TRIM(message_text))

    CALL message('mo_nwp_phy_init:', 'echam turbulence initialized')

  ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 3 ) THEN  !EDMF DUALM
    CALL suct0
    CALL su0phy
    CALL susekf
  ENDIF


  IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. .NOT. is_restart_run() ) THEN  ! TERRA
    CALL nwp_surface_init(p_patch, ext_data, p_prog_lnd_now,  &
      &                   p_prog_lnd_new, p_diag_lnd )
  END IF


  IF ( atm_phy_nwp_config(jg)%inwp_gwd == 1 ) THEN  ! IFS gwd scheme

     CALL sugwwms(nflevg= nlevp1, ppref=pref)
     CALL message('mo_nwp_phy_init:', 'non-orog GWs initialized')

  END IF


END SUBROUTINE init_nwp_phy


END MODULE mo_nwp_phy_init

