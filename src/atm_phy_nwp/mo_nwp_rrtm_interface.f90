!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (RRTM or Ritter-Geleyn).
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
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
MODULE mo_nwp_rrtm_interface

!  USE mo_atm_phy_nwp_nml,      ONLY: inwp_radiation, dt_rad, dt_radheat
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_exception,            ONLY: message,  finish !message_tex
  USE mo_ext_data,             ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, &
    &                                io3, ntracer, ntracer_static
  USE mo_grf_interpolation,    ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc!, min_rlcell 
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_ovlparea_start_c
  USE mo_interpolation,        ONLY: t_int_state
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_lnd_state,        ONLY: t_lnd_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch
  USE mo_mpi,                  ONLY: my_process_is_mpi_seq
  USE mo_phyparam_soil,        ONLY: csalb, csalb_snow_min, csalb_snow_max, &
    &                                csalb_snow_fe, csalb_snow_fd, csalb_p, cf_snow
  USE mo_phys_nest_utilities,  ONLY: upscale_rad_input, downscale_rad_output, &
    &                                upscale_rad_input_rg, downscale_rad_output_rg
  USE mo_nonhydro_state,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag !,prm_diag
  USE mo_o3_util,              ONLY: calc_o3_clim
  USE mo_physical_constants,   ONLY: amd, amo3, tmelt
  USE mo_radiation,            ONLY: radiation, pre_radiation_nwp_steps
  USE mo_radiation_config,     ONLY: irad_o3, irad_aero, vmr_co2, rad_csalbw
  USE mo_radiation_rg,         ONLY: fesft
  USE mo_radiation_rg_par,     ONLY: aerdis
  USE mo_satad,                ONLY: qsat_rho
  USE mo_subdivision,          ONLY: p_patch_local_parent
!   USE mo_sync,                 ONLY: SYNC_C, sync_patch_array_mult

  IMPLICIT NONE

  PRIVATE

!!$  PUBLIC :: parameters, &
!!$    &        types,      &
!!$    &        variables,  &
!!$    &        procedures

  PUBLIC ::  nwp_rrtm_radiation, nwp_rrtm_radiation_reduced, nwp_rrtm_ozon
  

  CHARACTER(len=*), PARAMETER:: version = '$Id$'

  REAL(wp), PARAMETER::  &
    & zaeops = 0.05_wp,   &
    & zaeopl = 0.2_wp,    &
    & zaeopu = 0.1_wp,    &
    & zaeopd = 1.9_wp,    &
    & ztrpt  = 30.0_wp,   &
    & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
    & zvobga = 0.007_wp /  19330.0_wp , &
    & zstbga = 0.045_wp /  19330.0_wp!, &
!      & zaeadk(1:3) = (/0.3876E-03_wp,0.6693E-02_wp,0.8563E-03_wp/)

CONTAINS


  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rrtm_ozon ( p_sim_time,pt_patch, &
    & pt_prog_rcf,pt_diag,prm_diag )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog_rcf !<the prognostic variables (with
    !< reduced calling frequency for tracers!
    TYPE(t_nh_diag), TARGET, INTENT(in)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag

    ! for Ritter-Geleyn radiation:
    REAL(wp):: zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    ! for ozone:
    REAL(wp):: &
      & zptop32(nproma,pt_patch%nblks_c), &
      & zo3_hm (nproma,pt_patch%nblks_c), &
      & zo3_top (nproma,pt_patch%nblks_c), &
      & zpbot32(nproma,pt_patch%nblks_c), &
      & zo3_bot (nproma,pt_patch%nblks_c)
    ! for aerosols with Ritter-Geleyn:
    REAL(wp):: &
      & zsign(nproma,pt_patch%nlevp1), &
      & zvdaes(nproma,pt_patch%nlevp1), &
      & zvdael(nproma,pt_patch%nlevp1), &
      & zvdaeu(nproma,pt_patch%nlevp1), &
      & zvdaed(nproma,pt_patch%nlevp1), &
!      & zaeadk(3  ), &
      & zaeqdo   (nproma,pt_patch%nblks_c), zaeqdn,                 &
      & zaequo   (nproma,pt_patch%nblks_c), zaequn,                 &
      & zaeqlo   (nproma,pt_patch%nblks_c), zaeqln,                 &
      & zaeqso   (nproma,pt_patch%nblks_c), zaeqsn


    ! Local scalars:
    INTEGER:: jc,jk,jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------


    !O3

    SELECT CASE (irad_o3)
    CASE (6)
      CALL calc_o3_clim(                             &
        & kbdim      = nproma,                       & ! in
        & p_inc_rad  = atm_phy_nwp_config(jg)%dt_rad,& ! in
        & z_sim_time = p_sim_time,                   & ! in
        & pt_patch   = pt_patch,                     & ! in
        & zvio3      = prm_diag%vio3,                & !inout
        & zhmo3      = prm_diag%hmo3  )                !inout 
    END SELECT

    IF (ntracer + ntracer_static < io3)  RETURN ! no ozon

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx, &
!$OMP       zsign,zvdaes, zvdael, zvdaeu, zvdaed, zaeqsn, zaeqln, zaequn,zaeqdn )
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                        i_startidx, i_endidx, rl_start, rl_end)

      IF ( irad_o3 == 6 .AND. irad_aero == 5 ) THEN

        DO jk = 2, nlevp1
          DO jc = 1,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma, &
          & jcs    = 1, &
          & jce    = i_endidx, &
          & klevp1 = nlevp1, & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) ) !out

        ! 3-dimensional O3
        ! top level
        ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
        DO jc = 1,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zaeqso   (jc,jb) = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,1)
          zaequo   (jc,jb) = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,1)
          zaeqdo   (jc,jb) = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,1)
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
          DO jc = 1,i_endidx
            zaeqsn         = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,jk+1)
            zaeqln         = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,jk+1)
            zaequn         = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,jk+1)
            zaeqdn         = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,jk+1)


            zaeqso(jc,jb)    = zaeqsn
            zaeqlo(jc,jb)    = zaeqln
            zaequo(jc,jb)    = zaequn
            zaeqdo(jc,jb)    = zaeqdn

            zpbot32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,jk+1,jb)))**3
            zo3_bot  (jc,jb) = prm_diag%vio3(jc,jb)* zpbot32(jc,jb)    &
              /( zpbot32(jc,jb) + zo3_hm(jc,jb))
            !O3 content
            zduo3(jc,jk,jb)= zo3_bot (jc,jb)-zo3_top (jc,jb)
            ! store previous bottom values in arrays for top of next layer
            zo3_top (jc,jb) = zo3_bot (jc,jb)
          ENDDO
        ENDDO
        DO jk = 1,nlev
          ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
          DO jc = 1,i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,io3) = &
              & (amo3/amd) * (zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb))
          ENDDO
        ENDDO
      ELSEIF ( irad_o3 == 6 ) THEN !ozone, but no aerosols
         ! 3-dimensional O3
        ! top level
        ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
        DO jc = 1,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO
        ! loop over layers
        DO jk = 1,nlev
          ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
          DO jc = 1,i_endidx
            zpbot32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,jk+1,jb)))**3
            zo3_bot  (jc,jb) = prm_diag%vio3(jc,jb)* zpbot32(jc,jb)    &
              /( zpbot32(jc,jb) + zo3_hm(jc,jb))
            !O3 content
            zduo3(jc,jk,jb)= zo3_bot (jc,jb)-zo3_top (jc,jb)
            ! store previous bottom values in arrays for top of next layer
            zo3_top (jc,jb) = zo3_bot (jc,jb)
          ENDDO
        ENDDO
        DO jk = 1,nlev
          ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
          DO jc = 1,i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,io3) = &
              & (amo3/amd) * (zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb))
          ENDDO
        ENDDO
      ENDIF !irad_o3

    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE nwp_rrtm_ozon
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rrtm_radiation ( p_sim_time,pt_patch, &
    & ext_data,pt_prog_rcf,pt_diag,prm_diag, &
    & lnd_prog_now )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),INTENT(in):: ext_data
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog_rcf !<the prognostic variables (with
    !< reduced calling frequency for tracers!
    TYPE(t_nh_diag), TARGET, INTENT(in)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog_now

    REAL(wp):: albvisdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albvisdif     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdif     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: aclcov        (nproma,pt_patch%nblks_c) !<

    INTEGER :: itype(nproma)   !< type of convection

    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    INTEGER:: jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            &
      & cosmu0_dark  = cosmu0_dark,                       &
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     &
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_radheat ,&
      & p_sim_time   = p_sim_time,                        &
      & pt_patch     = pt_patch,                          &
      & zsmu0        = prm_diag%cosmu0(:,:),              &
      & zsct         = zsct )

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

    IF (msg_level >= 12) &
      &           CALL message('mo_nwp_rad_interface', 'RRTM radiation on full grid')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,itype),SCHEDULE(guided)
    !
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)


      ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
      itype(1:i_endidx) = 0 !INT(field%rtype(1:i_endidx,jb))

      albvisdir(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      albnirdir(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      albvisdif(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      albnirdif(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      prm_diag%tsfctrad(1:i_endidx,jb) = lnd_prog_now%t_g(1:i_endidx,jb)


      CALL radiation(               &
                              !
                              ! input
                              ! -----
                              !
                              ! indices and dimensions
        & jce        =i_endidx             ,&!< in  end   index for loop over block
        & kbdim      =nproma               ,&!< in  dimension of block over cells
        & klev       =nlev                 ,&!< in  number of full levels = number of layers
        & klevp1     =nlevp1               ,&!< in  number of half levels = number of layer ifcs
                              !
        & ktype      =itype                ,&!< in     type of convection
                              !
                              ! surface: albedo + temperature
        & zland      =ext_data%atm%fr_land_smt(:,jb)   ,&!< in     land fraction
        & zglac      =ext_data%atm%fr_glac_smt(:,jb)   ,&!< in     land glacier fraction
                              !
        & cos_mu0    =prm_diag%cosmu0  (:,jb) ,&!< in  cos of zenith angle mu0
        & alb_vis_dir=albvisdir        (:,jb) ,&!< in surface albedo for visible range, direct
        & alb_nir_dir=albnirdir        (:,jb) ,&!< in surface albedo for near IR range, direct
        & alb_vis_dif=albvisdif        (:,jb) ,&!< in surface albedo for visible range, diffuse
        & alb_nir_dif=albnirdif        (:,jb) ,&!< in surface albedo for near IR range, diffuse
        & tk_sfc     =prm_diag%tsfctrad(:,jb) ,&!< in surface temperature
                              !
                              ! atmosphere: pressure, tracer mixing ratios and temperature
        & pp_hl      =pt_diag%pres_ifc  (:,:,jb)     ,&!< in  pres at half levels at t-dt [Pa]
        & pp_fl      =pt_diag%pres      (:,:,jb)     ,&!< in  pres at full levels at t-dt [Pa]
        & tk_fl      =pt_diag%temp      (:,:,jb)     ,&!< in  temperature at full level at t-dt
        & qm_vap     =prm_diag%tot_cld  (:,:,jb,iqv) ,&!< in  water vapor mass mix ratio at t-dt
        & qm_liq     =prm_diag%tot_cld  (:,:,jb,iqc) ,&!< in cloud water mass mix ratio at t-dt
        & qm_ice     =prm_diag%tot_cld  (:,:,jb,iqi) ,&!< in cloud ice mass mixing ratio at t-dt
        & qm_o3      =pt_prog_rcf%tracer(:,:,jb,io3) ,&!< in o3 mass mixing ratio at t-dt
        & cdnc       =prm_diag%acdnc    (:,:,jb)     ,&!< in  cloud droplet numb conc. [1/m**3]
        & cld_frc    =prm_diag%tot_cld  (:,:,jb,icc) ,&!< in  cloud fraction [m2/m2]
                              !
                              ! output
                              ! ------
                              !
        & cld_cvr    =aclcov             (:,jb),&!< out cloud cover in a column [m2/m2]
        & emter_clr  =prm_diag%lwflxclr(:,:,jb),&!< out terrestrial flux, clear sky, net down
        & trsol_clr  =prm_diag%trsolclr(:,:,jb),&!< out sol. transmissivity, clear sky, net down
        & emter_all  =prm_diag%lwflxall(:,:,jb),&!< out terrestrial flux, all sky, net down
        & trsol_all  =prm_diag%trsolall(:,:,jb),&!< out solar transmissivity, all sky, net down
        & opt_halo_cosmu0 = .FALSE. )

      ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_rrtm_radiation
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rrtm_radiation_reduced ( p_sim_time,pt_patch,pt_par_patch, &
    & pt_par_int_state, pt_par_grf_state,ext_data,pt_prog_rcf,pt_diag,prm_diag, &
    & lnd_prog_now )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_int_state),    TARGET,INTENT(in):: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(in) :: pt_par_grf_state  !< grid refinement state
    TYPE(t_external_data),INTENT(in):: ext_data
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog_rcf !<the prognostic variables (with
    !< reduced calling frequency for tracers!
    TYPE(t_nh_diag), TARGET,    INTENT(inout):: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog_now

    REAL(wp):: albvisdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albvisdif     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdif     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: aclcov        (nproma,pt_patch%nblks_c) !<
    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_land  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_glac  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_cosmu0   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albvisdir(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albnirdir(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albvisdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albnirdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tsfc     (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_pres_ifc (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_pres     (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_temp     (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_o3       (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_acdnc    (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tot_cld  (:,:,:,:)
    ! Output fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aclcov   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxclr (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolclr (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolall (:,:,:)
    ! Pointer to parent patach or local parent patch for reduced grid
    TYPE(t_patch), POINTER       :: ptr_pp

    INTEGER :: itype(nproma)   !< type of convection

    REAL(wp):: alb_ther    (nproma,pt_patch%nblks_c) !!


    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    INTEGER:: jk,jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels
    INTEGER:: nblks_par_c       !nblks for reduced grid

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index
    INTEGER:: i_chidx
    LOGICAL:: l_parallel

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            &
      & cosmu0_dark  = cosmu0_dark,                       &
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     &
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_radheat, &
      & p_sim_time   = p_sim_time,                        &
      & pt_patch     = pt_patch,                          &
      & zsmu0        = prm_diag%cosmu0(:,:),              &
      & zsct         = zsct )

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      ! section for computing radiation on reduced grid

      IF (p_test_run) THEN
        prm_diag%lwflxall(:,:,:) = 0._wp
        prm_diag%trsolall(:,:,:) = 0._wp
      ENDIF

      IF (msg_level >= 12) &
        &       CALL message('mo_nwp_rad_interface', 'RRTM radiation on reduced grid')

      IF (my_process_is_mpi_seq()) THEN
        l_parallel = .FALSE.
      ELSE
        l_parallel = .TRUE.
      ENDIF

      i_chidx     =  pt_patch%parent_child_index

      IF (jg == 1 .OR. .NOT. l_parallel) THEN
        ptr_pp => pt_par_patch
        nblks_par_c = pt_par_patch%nblks_c

        ! number of vertical levels
        ! ** for the time being, the radiation grid is assumed to have the same
        !    levels as the main grid **
        ! nlev   = ptr_pp%nlev
        ! nlevp1 = ptr_pp%nlevp1
      ELSE ! Nested domain with MPI parallelization
        ptr_pp      => p_patch_local_parent(jg)
        nblks_par_c =  ptr_pp%nblks_c

        ! number of vertical levels
        ! ** for the time being, the radiation grid is assumed to have the same
        !    levels as the main grid **
        ! nlev   = ptr_pp%nlev
        ! nlevp1 = ptr_pp%nlevp1
      ENDIF

      ALLOCATE (zrg_cosmu0   (nproma,nblks_par_c),          &
        zrg_fr_land  (nproma,nblks_par_c),          &
        zrg_fr_glac  (nproma,nblks_par_c),          &
        zrg_albvisdir(nproma,nblks_par_c),          &
        zrg_albnirdir(nproma,nblks_par_c),          &
        zrg_albvisdif(nproma,nblks_par_c),          &
        zrg_albnirdif(nproma,nblks_par_c),          &
        zrg_tsfc     (nproma,nblks_par_c),          &
        zrg_pres_ifc (nproma,nlevp1,nblks_par_c),   &
        zrg_pres     (nproma,nlev  ,nblks_par_c),   &
        zrg_temp     (nproma,nlev  ,nblks_par_c),   &
        zrg_o3       (nproma,nlev  ,nblks_par_c),   &
        zrg_acdnc    (nproma,nlev  ,nblks_par_c),   &
        zrg_tot_cld  (nproma,nlev  ,nblks_par_c,4), &
        zrg_aclcov   (nproma,       nblks_par_c),   &
        zrg_lwflxclr (nproma,nlevp1,nblks_par_c),   &
        zrg_lwflxall (nproma,nlevp1,nblks_par_c),   &
        zrg_trsolclr (nproma,nlevp1,nblks_par_c),   &
        zrg_trsolall (nproma,nlevp1,nblks_par_c)    )

      rl_start = 1 ! SR radiation is not set up to handle boundaries of nested domains
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


      ! *** this parallel section will be removed later once real data are
      !     are available as input for radiation ***
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
      !
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          &                       i_startidx, i_endidx, rl_start, rl_end)


        ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
        albvisdir(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
        albnirdir(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
        albvisdif(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
        albnirdif(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
        prm_diag%tsfctrad(1:i_endidx,jb) = lnd_prog_now%t_g(1:i_endidx,jb)

      ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

      CALL upscale_rad_input(pt_patch, pt_par_patch, pt_par_grf_state,  &
        & ext_data%atm%fr_land_smt, ext_data%atm%fr_glac_smt,           &
        & prm_diag%cosmu0, albvisdir, albnirdir, albvisdif, albnirdif,  &
        & prm_diag%tsfctrad, pt_diag%pres_ifc,                          &
        & pt_diag%pres, pt_diag%temp,prm_diag%acdnc, prm_diag%tot_cld,  &
        & pt_prog_rcf%tracer(:,:,:,io3),                                &
        & zrg_fr_land, zrg_fr_glac,                                     &
        & zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,      &
        & zrg_albnirdif, zrg_tsfc, zrg_pres_ifc, zrg_pres, zrg_temp,    &
        & zrg_acdnc, zrg_tot_cld, zrg_o3                              )

      rl_start = grf_ovlparea_start_c
      rl_end   = min_rlcell_int

      i_startblk = ptr_pp%cells%start_blk(rl_start,i_chidx)
      i_endblk   = ptr_pp%cells%end_blk(rl_end,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,itype) ,SCHEDULE(guided)
      !
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
          &                         i_startidx, i_endidx, rl_start, rl_end, i_chidx)

        itype(:) = 0

        ! Unfortunately, the coding of SR radiation is not compatible with the presence
        ! of nested domains. Therefore, the normally unused elements of the first block
        ! need to be filled with dummy values
        IF (jg > 1 .AND. jb == i_startblk) THEN
          zrg_fr_land   (1:i_startidx-1,jb) = zrg_fr_land(i_startidx,jb)
          zrg_fr_glac   (1:i_startidx-1,jb) = zrg_fr_glac(i_startidx,jb)
          zrg_cosmu0    (1:i_startidx-1,jb) = zrg_cosmu0 (i_startidx,jb)
          zrg_albvisdir (1:i_startidx-1,jb) = zrg_albvisdir (i_startidx,jb)
          zrg_albnirdir (1:i_startidx-1,jb) = zrg_albnirdir (i_startidx,jb)
          zrg_albvisdif (1:i_startidx-1,jb) = zrg_albvisdif (i_startidx,jb)
          zrg_albnirdif (1:i_startidx-1,jb) = zrg_albnirdif (i_startidx,jb)
          zrg_tsfc      (1:i_startidx-1,jb) = zrg_tsfc      (i_startidx,jb)
          zrg_pres_ifc (1:i_startidx-1,nlevp1,jb) = zrg_pres_ifc (i_startidx,nlevp1,jb)
          DO jk = 1, nlev
            zrg_pres_ifc (1:i_startidx-1,jk,jb) = zrg_pres_ifc (i_startidx,jk,jb)
            zrg_pres     (1:i_startidx-1,jk,jb) = zrg_pres     (i_startidx,jk,jb)
            zrg_temp     (1:i_startidx-1,jk,jb) = zrg_temp     (i_startidx,jk,jb)
            zrg_o3       (1:i_startidx-1,jk,jb) = zrg_o3       (i_startidx,jk,jb)
            zrg_acdnc    (1:i_startidx-1,jk,jb) = zrg_acdnc    (i_startidx,jk,jb)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqv) = zrg_tot_cld(i_startidx,jk,jb,iqv)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqc) = zrg_tot_cld(i_startidx,jk,jb,iqc)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqi) = zrg_tot_cld(i_startidx,jk,jb,iqi)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,icc) = zrg_tot_cld(i_startidx,jk,jb,icc)
          ENDDO
        ENDIF

!#ifdef __BOUNDCHECK
        CALL radiation(               &
                                !
                                ! input
                                ! -----
                                !
                                ! indices and dimensions
          & jce         =i_endidx            ,&!< in  end   index for loop over block
          & kbdim       =nproma              ,&!< in  dimension of block over cells
          & klev        =nlev                ,&!< in  number of full levels = number of layers
          & klevp1      =nlevp1              ,&!< in  number of half levels = number of layer ifcs
                                !
          & ktype       =itype               ,&!< in   type of convection
                                !
          & zland       =zrg_fr_land (:,jb)  ,&!< in land mask,     1. over land
          & zglac       =zrg_fr_glac (:,jb)  ,&!< in glacier mask,  1. over land ice
                                !
          & cos_mu0     =zrg_cosmu0  (:,jb)  ,&!< in    cos of zenith angle mu0
          & alb_vis_dir=zrg_albvisdir(:,jb)  ,&!< in    surface albedo for visible range, direct
          & alb_nir_dir=zrg_albnirdir(:,jb)  ,&!< in    surface albedo for near IR range, direct
          & alb_vis_dif=zrg_albvisdif(:,jb)  ,&!< in    surface albedo for visible range, diffuse
          & alb_nir_dif=zrg_albnirdif(:,jb)  ,&!< in    surface albedo for near IR range, diffuse
          & tk_sfc     =zrg_tsfc     (:,jb)       ,&!< in    surface temperature
                                !
                                ! atmosphere: pressure, tracer mixing ratios and temperature
          & pp_hl      =zrg_pres_ifc(:,:,jb)    ,&!< in    pressure at half levels at t-dt [Pa]
          & pp_fl      =zrg_pres    (:,:,jb)    ,&!< in    pressure at full levels at t-dt [Pa]
          & tk_fl      =zrg_temp    (:,:,jb)    ,&!< in    temperature at full level at t-dt
          & qm_vap     =zrg_tot_cld (:,:,jb,iqv),&!< in    water vapor mass mixing ratio at t-dt
          & qm_liq     =zrg_tot_cld (:,:,jb,iqc),&!< in    cloud water mass mixing ratio at t-dt
          & qm_ice     =zrg_tot_cld (:,:,jb,iqi),&!< in    cloud ice mass mixing ratio at t-dt
          & qm_o3      = zrg_o3     (:,:,jb)    ,&!< in    O3
          & cdnc       =zrg_acdnc   (:,:,jb)    ,&!< in    cloud droplet numb. conc. [1/m**3]
          & cld_frc    =zrg_tot_cld (:,:,jb,icc),&!< in    cld_frac = cloud fraction [m2/m2]
                                !
                                ! output
                                ! ------
                                !
          & cld_cvr    =zrg_aclcov    (:,jb)  ,&!< out   cloud cover in a column [m2/m2]
          & emter_clr  =zrg_lwflxclr(:,:,jb)  ,&!< out   LW (terrestrial) flux,clear sky, net down
          & trsol_clr  =zrg_trsolclr(:,:,jb)  ,&!< out   solar transmissivity, clear sky, net down
          & emter_all  =zrg_lwflxall(:,:,jb)  ,&!< out   terrestrial flux, all sky, net down
          & trsol_all  =zrg_trsolall(:,:,jb)  ,&!< out   solar transmissivity, all sky, net down
          & opt_halo_cosmu0 = .FALSE. )

      ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

      CALL downscale_rad_output(pt_patch, pt_par_patch, pt_par_int_state,         &
        & pt_par_grf_state, zrg_aclcov, zrg_lwflxclr, zrg_lwflxall, zrg_trsolclr, &
        & zrg_trsolall, aclcov, prm_diag%lwflxclr, prm_diag%lwflxall, &
        & prm_diag%trsolclr, prm_diag%trsolall )

      DEALLOCATE (zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif, zrg_albnirdif, &
        zrg_tsfc, zrg_pres_ifc, zrg_pres, zrg_temp, zrg_o3, zrg_acdnc, zrg_tot_cld,       &
        zrg_aclcov, zrg_lwflxclr, zrg_lwflxall, zrg_trsolclr, zrg_trsolall,     &
        zrg_fr_land,zrg_fr_glac)
      
  END SUBROUTINE nwp_rrtm_radiation_reduced


END MODULE mo_nwp_rrtm_interface

