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

MODULE mo_nwp_rad_interface

  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_exception,            ONLY: message,  finish !message_tex
  USE mo_ext_data,             ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run, parallel_radiation_omp

  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, &
    &                                io3, ntracer, ntracer_static
  USE mo_grf_interpolation,    ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc, io3_ape!, min_rlcell 
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
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_radiation_config,     ONLY: irad_o3, irad_aero, vmr_co2, rad_csalbw
  USE mo_radiation_rg,         ONLY: fesft
  USE mo_radiation_rg_par,     ONLY: aerdis
  USE mo_satad,                ONLY: qsat_rho
  USE mo_subdivision,          ONLY: p_patch_local_parent
!   USE mo_sync,                 ONLY: SYNC_C, sync_patch_array_mult

  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation, &
   &  nwp_rrtm_radiation_reduced, nwp_rrtm_ozon_aerosol
  USE mo_nwp_mpiomp_rrtm_interface, ONLY: nwp_omp_rrtm_interface
  
  IMPLICIT NONE

  PRIVATE

!!$  PUBLIC :: parameters, &
!!$    &        types,      &
!!$    &        variables,  &
!!$    &        procedures

  PUBLIC :: nwp_radiation
  

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
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1) or Ritter-Geleyn (2).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time,pt_patch,pt_par_patch, &
    & pt_par_int_state,pt_par_grf_state,ext_data,lnd_diag,pt_prog,pt_prog_rcf,pt_diag,prm_diag, &
    & lnd_prog_now )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'
    
    LOGICAL, INTENT(in)         :: lredgrid        !< use reduced grid for radiation

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_int_state),    TARGET,INTENT(in):: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(in) :: pt_par_grf_state  !< grid refinement state
    TYPE(t_external_data),INTENT(in):: ext_data
    TYPE(t_lnd_diag),     INTENT(in):: lnd_diag      !< diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog_rcf !<the prognostic variables (with
    !< reduced calling frequency for tracers!
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog_now

    REAL(wp) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    
    INTEGER :: jg

    jg = pt_patch%id


    IF (atm_phy_nwp_config(jg)%inwp_radiation == 1 ) THEN
       
      CALL nwp_rrtm_ozon_aerosol ( p_sim_time,pt_patch, &
        & pt_prog_rcf,pt_diag,prm_diag,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )
    
      IF ( .NOT. lredgrid ) THEN

        IF (parallel_radiation_omp) THEN
          CALL nwp_omp_rrtm_interface ( p_sim_time,pt_patch, &
            & ext_data, lnd_diag, pt_prog_rcf, pt_diag, prm_diag, lnd_prog_now )
        ELSE
          CALL nwp_rrtm_radiation ( p_sim_time,pt_patch, &
            & ext_data,lnd_diag,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5, &
            & pt_prog_rcf, pt_diag, prm_diag, lnd_prog_now )
        ENDIF
    
      ELSE 

        CALL nwp_rrtm_radiation_reduced ( p_sim_time,pt_patch,pt_par_patch,  &
          & pt_par_int_state, pt_par_grf_state,ext_data,lnd_diag,            &
          & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,                                   &
          & pt_prog_rcf,pt_diag,prm_diag, lnd_prog_now )
          
      ENDIF

      RETURN
    ENDIF !inwp_radiation = 1


    IF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 .AND. .NOT. lredgrid) THEN
    
      CALL nwp_rg_radiation ( p_sim_time,pt_patch, &
        & ext_data,lnd_diag, &
        & pt_prog,pt_prog_rcf,pt_diag,prm_diag, lnd_prog_now )


    ELSEIF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 .AND. lredgrid) THEN

      CALL nwp_rg_radiation_reduced ( p_sim_time,pt_patch,pt_par_patch, &
        & pt_par_int_state, pt_par_grf_state,ext_data,lnd_diag, &
        & pt_prog,pt_prog_rcf,pt_diag,prm_diag, lnd_prog_now )


    ENDIF !inwp_radiation

  END SUBROUTINE nwp_radiation
  !---------------------------------------------------------------------------------------



  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rg_radiation ( p_sim_time,pt_patch, &
    & ext_data,lnd_diag,pt_prog,pt_prog_rcf,pt_diag,prm_diag, &
    & lnd_prog_now )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

!    REAL(wp), PARAMETER::  &
!      & zqco2 = 0.5014E-03_wp*353.9_wp/330._wp ! CO2 (mixing ratio 353.9 ppm (like vmr_co2))

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),INTENT(in):: ext_data
    TYPE(t_lnd_diag),     INTENT(in):: lnd_diag      !< diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog_rcf !<the prognostic variables (with
    !< reduced calling frequency for tracers!
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog_now

    REAL(wp):: albvisdif     (nproma,pt_patch%nblks_c) !<

    REAL(wp):: zi0        (nproma)  !< solar incoming radiation at TOA   [W/m2]
    ! for Ritter-Geleyn radiation:
    REAL(wp):: zqco2
    REAL(wp):: zsqv (nproma,pt_patch%nlev,pt_patch%nblks_c) !< saturation water vapor mixing ratio
    REAL(wp):: zfls (nproma,pt_patch%nlevp1 ,pt_patch%nblks_c)
    REAL(wp):: zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zduco2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: alb_ther    (nproma,pt_patch%nblks_c) !!
    LOGICAL :: lo_sol (nproma)
    LOGICAL :: losol
    ! For Ritter-Geleyn radiation on reduced grid additionally
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
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
      & zaetr_top(nproma,pt_patch%nblks_c), zaetr_bot, zaetr, &
      &  zaeqdo   (nproma,pt_patch%nblks_c), zaeqdn,                 &
      &  zaequo   (nproma,pt_patch%nblks_c), zaequn,                 &
      & zaeqlo   (nproma,pt_patch%nblks_c), zaeqln,                 &
      & zaeqso   (nproma,pt_patch%nblks_c), zaeqsn


    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    REAL(wp):: zvege, zsnow, zsalb_snow, zsnow_alb
    INTEGER:: jc,jk,jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index
    INTEGER:: ist

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! CO2 help variable for Ritter-Geleyn scheme
    zqco2 = 0.5014E-03_wp*vmr_co2/330.e-6_wp

    ! O3
    SELECT CASE (irad_o3)
    CASE(io3_ape)
      ! APE ozone: do nothing because 
      ! everything is already set in nwp_phy_init
    CASE (6)
      CALL calc_o3_clim(                             &
        & kbdim      = nproma,                       &
        & p_inc_rad  = atm_phy_nwp_config(jg)%dt_rad,&
        & z_sim_time = p_sim_time,                   &
        & pt_patch   = pt_patch,                     &
        & zvio3      = prm_diag%vio3,                &
        & zhmo3      = prm_diag%hmo3  )
    END SELECT

    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            &
      & cosmu0_dark  = cosmu0_dark,                       &
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad     ,&
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_radheat, &
      & p_sim_time   = p_sim_time,                        &
      & pt_patch     = pt_patch,                          &
     !& zsmu0        = prm_diag%cosmu0(1,1),              &
      & zsmu0        = prm_diag%cosmu0(:,:),              &
      & zsct         = zsct )

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx, &
!$OMP       zsign,zvdaes, zvdael, zvdaeu, zvdaed, zaeqsn, zaeqln, zaequn,zaeqdn,zaetr_bot,zaetr )
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                        i_startidx, i_endidx, rl_start, rl_end)

      IF ( irad_o3 == 6 .AND. irad_aero == 5 ) THEN

        DO jk = 2, nlevp1
          DO jc = i_startidx,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = i_startidx,  & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! 3-dimensional O3
        ! top level
        DO jc = i_startidx,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zaeqso   (jc,jb) = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,1)
          zaequo   (jc,jb) = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,1)
          zaeqdo   (jc,jb) = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zaeqsn         = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,jk+1)
            zaeqln         = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,jk+1)
            zaequn         = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,jk+1)
            zaeqdn         = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,jk+1)
            zaetr_bot      = zaetr_top(jc,jb) &
              & * ( MIN (1.0_wp, pt_diag%temp_ifc(jc,jk,jb)/pt_diag%temp_ifc(jc,jk+1,jb)) )**ztrpt

            zaetr          = SQRT(zaetr_bot*zaetr_top(jc,jb))
            zaeq1(jc,jk,jb)= (1._wp-zaetr) &
              & * (ztrbga* pt_diag%dpres_mc(jc,jk,jb)+zaeqln-zaeqlo(jc,jb)+zaeqdn-zaeqdo(jc,jb))
            zaeq2(jc,jk,jb)   = (1._wp-zaetr) * ( zaeqsn-zaeqso(jc,jb) )
            zaeq3(jc,jk,jb)   = (1._wp-zaetr) * ( zaequn-zaequo(jc,jb) )
            zaeq4(jc,jk,jb)   =     zaetr  *   zvobga*pt_diag%dpres_mc(jc,jk,jb)
            zaeq5(jc,jk,jb)   =     zaetr  *   zstbga*pt_diag%dpres_mc(jc,jk,jb)

            zaetr_top(jc,jb) = zaetr_bot
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
        IF (ntracer + ntracer_static >= io3) THEN
          DO jk = 1,nlev
            DO jc = i_startidx,i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,io3) = &
                & (amo3/amd) * (zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF
      ELSEIF ( irad_o3 == 6 ) THEN !ozone, but no aerosols
         ! 3-dimensional O3
        ! top level
        DO jc = i_startidx,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO
        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zpbot32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,jk+1,jb)))**3
            zo3_bot  (jc,jb) = prm_diag%vio3(jc,jb)* zpbot32(jc,jb)    &
              /( zpbot32(jc,jb) + zo3_hm(jc,jb))
            !O3 content
            zduo3(jc,jk,jb)= zo3_bot (jc,jb)-zo3_top (jc,jb)
            ! store previous bottom values in arrays for top of next layer
            zo3_top (jc,jb) = zo3_bot (jc,jb)
          ENDDO
        ENDDO
        IF (ntracer + ntracer_static >= io3) THEN
          DO jk = 1,nlev
            DO jc = i_startidx,i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,io3) = &
                & (amo3/amd) * (zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF
        zaeq1(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq2(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq3(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq4(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq5(i_startidx:i_endidx,:,jb) = 0.0_wp
      ELSEIF ( irad_aero == 5 ) THEN !aerosols, but no ozone:

        DO jk = 2, nlevp1
          DO jc = i_startidx,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = i_startidx,  & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! top level
        DO jc = i_startidx,i_endidx
          zaeqso   (jc,jb) = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,1)
          zaequo   (jc,jb) = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,1)
          zaeqdo   (jc,jb) = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zaeqsn         = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,jk+1)
            zaeqln         = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,jk+1)
            zaequn         = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,jk+1)
            zaeqdn         = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,jk+1)
            zaetr_bot      = zaetr_top(jc,jb) &
              & * ( MIN (1.0_wp, pt_diag%temp_ifc(jc,jk,jb)/pt_diag%temp_ifc(jc,jk+1,jb)) )**ztrpt

            zaetr          = SQRT(zaetr_bot*zaetr_top(jc,jb))
            zaeq1(jc,jk,jb)= (1._wp-zaetr) &
              & * (ztrbga* pt_diag%dpres_mc(jc,jk,jb)+zaeqln-zaeqlo(jc,jb)+zaeqdn-zaeqdo(jc,jb))
            zaeq2(jc,jk,jb)   = (1._wp-zaetr) * ( zaeqsn-zaeqso(jc,jb) )
            zaeq3(jc,jk,jb)   = (1._wp-zaetr) * ( zaequn-zaequo(jc,jb) )
            zaeq4(jc,jk,jb)   =     zaetr  *   zvobga*pt_diag%dpres_mc(jc,jk,jb)
            zaeq5(jc,jk,jb)   =     zaetr  *   zstbga*pt_diag%dpres_mc(jc,jk,jb)

            zaetr_top(jc,jb) = zaetr_bot
            zaeqso(jc,jb)    = zaeqsn
            zaeqlo(jc,jb)    = zaeqln
            zaequo(jc,jb)    = zaequn
            zaeqdo(jc,jb)    = zaeqdn

            zduo3(jc,jk,jb) = &
              & pt_prog_rcf%tracer(jc,jk,jb,io3) * pt_diag%dpres_mc(jc,jk,jb) * amd / amo3
            
          ENDDO
        ENDDO

      ELSEIF (irad_o3 /= 0) THEN !no aerosols and other ozone

        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zduo3(jc,jk,jb) = &
              & pt_prog_rcf%tracer(jc,jk,jb,io3) * pt_diag%dpres_mc(jc,jk,jb) * amd / amo3
            zaeq1(jc,jk,jb) = 0.0_wp
            zaeq2(jc,jk,jb) = 0.0_wp
            zaeq3(jc,jk,jb) = 0.0_wp
            zaeq4(jc,jk,jb) = 0.0_wp
            zaeq5(jc,jk,jb) = 0.0_wp
          ENDDO
        ENDDO

      ELSE !no aerosols and no ozone

        zaeq1(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq2(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq3(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq4(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq5(i_startidx:i_endidx,:,jb) = 0.0_wp
        zduo3(i_startidx:i_endidx,:,jb) = 0.0_wp
        
      ENDIF !irad_o3

    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


    IF (msg_level >= 12) &
      &           CALL message('mo_nwp_rad_interface', 'RG radiation on full grid')

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,zi0,losol,lo_sol, &
!$OMP            ist,zvege,zsnow,zsalb_snow,zsnow_alb),SCHEDULE(guided)
    
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zsqv (jc,jk,jb) = qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDDO

      ! geographical dependent thermal albedo
      alb_ther(i_startidx:i_endidx,jb) = 1._wp-ext_data%atm%emis_rad(i_startidx:i_endidx,jb)

      !------------------------------------------------------------------------------
      ! Calculation of surface albedo taking soil type,              
      ! vegetation and snow/ice conditions into account
      !------------------------------------------------------------------------------
      
      IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        
        DO jc = i_startidx,i_endidx

          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog_now%t_g(jc,jb) >= tmelt-1.7_wp ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF
          albvisdif(jc,jb) = csalb(ist)
          IF ( ext_data%atm%llsm_atm_c(jc,jb)) THEN
            ! ATTENTION: only valid, if nsfc_subs=1
            albvisdif(jc,jb) = csalb(ist) - rad_csalbw(ist)*lnd_prog_now%w_so(jc,1,jb,1)
          ENDIF  ! lsoil, llandmask
          
        ENDDO
        
      ELSE
        
        DO jc = i_startidx,i_endidx
          
          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog_now%t_g(jc,jb) >= tmelt-1.7_wp ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF
          albvisdif(jc,jb) = csalb(ist)
          
        ENDDO
      
      ENDIF

      !sea-ice model not yet implemented
!      IF (atm_phy_nwp_config(jg)%lseaice) THEN
!        DO jc = i_startidx,i_endidx
!          ! In case the sea ice model is used AND water point AND ice is present,
!          ! compute ice albedo for water points with an empirical formula taken from GME.
!          ! The ice albedo is the lower the warmer, and therefore wetter, the ice is.
!          ! Use ice temperature at time level nnow (2-time level scheme in sea ice model).
!!          IF ((.NOT. llandmask(i,j)) .AND. (h_ice(i,j,nnow) > 0.0_ireals))               &
!!            zalso(i,j) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(tmelt-t_ice(i,j,nnow)))) &
!!            * csalb(10)
!          IF (( .NOT. ext_data%atm%llsm_atm_c(jc,jb)) .AND. &
!            & (prm_diag%h_ice(i,j,nnow) > 0.0_wp))          &
!            albvisdif(jc,jb) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(tmelt-t_ice(i,j,nnow)))) &
!            * csalb(10)
!        ENDDO
!      ENDIF

      !lake model not yet implemented
!      IF (atm_phy_nwp_config(jg)%llake) THEN
!        DO jc = i_startidx,i_endidx
!          IF((ext_data%atm%depth_lk(jc,jb)      >  0.0_wp) .AND.    &
!            (prm_diag%h_ice(jc,jb) >= h_Ice_min_flk) ) THEN
!            !  In case the lake model FLake is used AND lake point AND ice is present,
!            !  compute ice albedo for lake points with an empirical formulation 
!            !  proposed by Mironov and Ritter (2004) for use in GME 
!            !  [ice_albedo=function(ice_surface_temperature)].
!            !  Use surface temperature at time level "nnow".
!
!            albvisdif(jc,jb) = EXP(-c_albice_MR*(tpl_T_f-t_s(i,j,nnow))/tpl_T_f)
!            albvisdif(jc,jb) = albedo_whiteice_ref * (1._ireals-zalso(i,j)) +      &
!              albedo_blueice_ref  * albvisdif(jc,jb)
!          ENDIF
!        ENDDO
!      ENDIF

      ! Snow cover and vegetation
      ! -------------------------

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN !soil model switched on
        DO jc = i_startidx,i_endidx
          zvege= 0.0_wp
          zsnow= 0.0_wp
          IF (ext_data%atm%llsm_atm_c(jc,jb)) THEN 
            ! consider effects of aging on solar snow albedo
            ! ATTENTION: only valid, if nsfc_subs=1            
            zsalb_snow = csalb_snow_min + &
              & lnd_diag%freshsnow(jc,jb,1)*(csalb_snow_max-csalb_snow_min)
            zsnow_alb = zsalb_snow*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb)) &
              + csalb_snow_fe * ext_data%atm%for_e(jc,jb)                       &
              + csalb_snow_fd * ext_data%atm%for_d(jc,jb)

          ! account for snow cover and plant cover and compute final solar
          ! snow albedo
            zvege = ext_data%atm%plcov_mx(jc,jb)
            ! ATTENTION: only valid, if nsfc_subs=1
            IF (lnd_prog_now%w_snow(jc,jb,1) > 0.0_wp) THEN
              ! ATTENTION: only valid, if nsfc_subs=1
              zsnow = MIN(1.0_wp, lnd_prog_now%w_snow(jc,jb,1) / cf_snow)
            ENDIF
            albvisdif(jc,jb) = zsnow * zsnow_alb +                               &
              (1.0_wp - zsnow) * (zvege * csalb_p + (1.0_wp - zvege) * albvisdif(jc,jb))
          ENDIF !   llandmask
        ENDDO
      ENDIF !inwp_surface == 1

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog_now%t_g(i_startidx:i_endidx,jb)

      ! CO2 (mixing ratio 353.9 ppm as vmr_co2)
      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zduco2(jc,jk,jb) = zqco2 * pt_diag%dpres_mc (jc,jk,jb)
        ENDDO
      ENDDO

      ! Switch off solar radiation calculations where sun is below horizon:
      WHERE ( prm_diag%cosmu0(i_startidx:i_endidx,jb) > 1.e-8_wp ) !zepmu0 )
        lo_sol(i_startidx:i_endidx) = .TRUE.
      ELSEWHERE
        lo_sol(i_startidx:i_endidx) = .FALSE.
      END WHERE
      losol = ANY(lo_sol(i_startidx:i_endidx))

!#ifdef __BOUNDCHECK
      CALL fesft ( &
                                !  Input:
        & pti = pt_diag%temp_ifc (:,:,jb) , &! Temperature at layer boundaries
        & pdp = pt_diag%dpres_mc (:,:,jb), &! pressure thickness
        & pclc_in= prm_diag%tot_cld  (:,:,jb,icc) , &
        & pqv = prm_diag%tot_cld(:,:,jb,iqv), &! pt_prog_rcf%tracer(:,:,jb,iqv)
        & pqvs = zsqv(:,:,jb), &!saturation water vapor
        & pqcwc = prm_diag%tot_cld    (:,:,jb,iqc) ,&
        & pqiwc = prm_diag%tot_cld    (:,:,jb,iqi) ,&
        & pduco2 = zduco2 (:,:,jb), &! layer CO2 content
        & pduo3  = zduo3(:,:,jb),&! layer O3 content
        & paeq1 = zaeq1(:,:,jb), &
        & paeq2 = zaeq2(:,:,jb),&
        & paeq3 = zaeq3(:,:,jb),&
        & paeq4 = zaeq4(:,:,jb),&
        & paeq5 = zaeq5(:,:,jb),&
        & papre_in =  pt_diag%pres_sfc (:,jb), & ! Surface pressure
        & psmu0 = prm_diag%cosmu0 (:,jb) , & ! Cosine of zenith angle
        & palso = albvisdif(:,jb), & ! solar surface albedo
        & palth = alb_ther(:,jb), & ! thermal surface albedo
        & psct = zsct, &! solar constant (at time of year)
        & kig1s = 1 ,&
        & kig1e = nproma , &
        & ki3s = 1, &
        & ki3e = nlev,&
        & ki1sc= i_startidx, &
        & ki1ec= i_endidx, &
        & lsolar = losol, &! control switch for solar calculations
        !          & lsolar = .TRUE., &! control switch for solar calculations
        & lthermal =.TRUE., &
        & lcrf = .FALSE., &! control switch for cloud-free calcul.
                                ! Output:
        & pflt  = prm_diag%lwflxall(:,:,jb),& !Thermal radiative fluxes at each layer boundary
        & pfls  = zfls  (:,:,jb)  &! solar radiative fluxes at each layer boundary
        & )

      zi0 (i_startidx:i_endidx) = prm_diag%cosmu0(i_startidx:i_endidx,jb) * zsct
      ! compute sw transmissivity trsolall from sw fluxes
      DO jk = 1,nlevp1
        DO jc = i_startidx,i_endidx
          IF (prm_diag%cosmu0(jc,jb) < 1.e-8_wp) THEN
            prm_diag%trsolall(jc,jk,jb) = 0.0_wp
          ELSE
            prm_diag%trsolall(jc,jk,jb) = zfls(jc,jk,jb) / zi0(jc)
          ENDIF
        ENDDO
      ENDDO
      
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE nwp_rg_radiation
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rg_radiation_reduced ( p_sim_time,pt_patch,pt_par_patch, &
    & pt_par_int_state, pt_par_grf_state,ext_data,lnd_diag,pt_prog,pt_prog_rcf,pt_diag,prm_diag, &
    & lnd_prog_now )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

!    REAL(wp), PARAMETER::  &
!      & zqco2 = 0.5014E-03_wp*353.9_wp/330._wp ! CO2 (mixing ratio 353.9 ppm (like vmr_co2))

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_int_state),    TARGET,INTENT(in):: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(in) :: pt_par_grf_state  !< grid refinement state
    TYPE(t_external_data),INTENT(in):: ext_data
    TYPE(t_lnd_diag),     INTENT(in):: lnd_diag      !< diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog_rcf !<the prognostic variables (with
    !< reduced calling frequency for tracers!
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog_now

    REAL(wp):: albvisdif     (nproma,pt_patch%nblks_c) !<
    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_cosmu0   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albvisdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_o3       (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tot_cld  (:,:,:,:)
    ! Output fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fls (:,:,:)
    ! Pointer to parent patach or local parent patch for reduced grid
    TYPE(t_patch), POINTER       :: ptr_pp

    REAL(wp):: zi0        (nproma)  !< solar incoming radiation at TOA   [W/m2]
    ! for Ritter-Geleyn radiation:
    REAL(wp):: zqco2
    REAL(wp):: zsqv (nproma,pt_patch%nlev,pt_patch%nblks_c) !< saturation water vapor mixing ratio
    REAL(wp):: zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zduco2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: alb_ther    (nproma,pt_patch%nblks_c) !!
    LOGICAL :: lo_sol (nproma)
    LOGICAL :: losol
    ! For Ritter-Geleyn radiation on reduced grid additionally
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_alb_ther(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_pres_sfc(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_temp_ifc(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_dpres_mc(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_sqv(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_duco2(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq1(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq2(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq3(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq4(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq5(:,:,:)
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
      & zaetr_top(nproma,pt_patch%nblks_c), zaetr_bot, zaetr, &
      &  zaeqdo   (nproma,pt_patch%nblks_c), zaeqdn,                 &
      &  zaequo   (nproma,pt_patch%nblks_c), zaequn,                 &
      & zaeqlo   (nproma,pt_patch%nblks_c), zaeqln,                 &
      & zaeqso   (nproma,pt_patch%nblks_c), zaeqsn
    REAL(wp), PARAMETER::  &
      & zaeops = 0.05_wp, &
      & zaeopl = 0.2_wp, &
      & zaeopu = 0.1_wp, &
      & zaeopd = 1.9_wp, &
      & ztrpt  = 30.0_wp, &
      & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
      & zvobga = 0.007_wp /  19330.0_wp , &
      & zstbga = 0.045_wp /  19330.0_wp!, &
!      & zaeadk(1:3) = (/0.3876E-03_wp,0.6693E-02_wp,0.8563E-03_wp/)


    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    REAL(wp):: zvege, zsnow, zsalb_snow, zsnow_alb
    INTEGER:: jc,jk,jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels
    INTEGER:: nblks_par_c       !nblks for reduced grid

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index
    INTEGER:: i_chidx
    INTEGER:: ist
    LOGICAL:: l_parallel

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! CO2 help variable for Ritter-Geleyn scheme
    zqco2 = 0.5014E-03_wp*vmr_co2/330.e-6_wp

    ! O3
    SELECT CASE (irad_o3)
    CASE (6)
      CALL calc_o3_clim(                             &
        & kbdim      = nproma,                       &
        & p_inc_rad  = atm_phy_nwp_config(jg)%dt_rad,&
        & z_sim_time = p_sim_time,                   &
        & pt_patch   = pt_patch,                     &
        & zvio3      = prm_diag%vio3,                &
        & zhmo3      = prm_diag%hmo3  )
    END SELECT

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

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx, &
!$OMP       zsign,zvdaes, zvdael, zvdaeu, zvdaed, zaeqsn, zaeqln, zaequn,zaeqdn,zaetr_bot,zaetr )
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                        i_startidx, i_endidx, rl_start, rl_end)

      IF ( irad_o3 == 6 .AND. irad_aero == 5 ) THEN

        DO jk = 2, nlevp1
          DO jc = i_startidx,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = i_startidx,  & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! 3-dimensional O3
        ! top level
        DO jc = i_startidx,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zaeqso   (jc,jb) = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,1)
          zaequo   (jc,jb) = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,1)
          zaeqdo   (jc,jb) = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zaeqsn         = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,jk+1)
            zaeqln         = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,jk+1)
            zaequn         = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,jk+1)
            zaeqdn         = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,jk+1)
            zaetr_bot      = zaetr_top(jc,jb) &
              & * ( MIN (1.0_wp, pt_diag%temp_ifc(jc,jk,jb)/pt_diag%temp_ifc(jc,jk+1,jb)) )**ztrpt

            zaetr          = SQRT(zaetr_bot*zaetr_top(jc,jb))
            zaeq1(jc,jk,jb)= (1._wp-zaetr) &
              & * (ztrbga* pt_diag%dpres_mc(jc,jk,jb)+zaeqln-zaeqlo(jc,jb)+zaeqdn-zaeqdo(jc,jb))
            zaeq2(jc,jk,jb)   = (1._wp-zaetr) * ( zaeqsn-zaeqso(jc,jb) )
            zaeq3(jc,jk,jb)   = (1._wp-zaetr) * ( zaequn-zaequo(jc,jb) )
            zaeq4(jc,jk,jb)   =     zaetr  *   zvobga*pt_diag%dpres_mc(jc,jk,jb)
            zaeq5(jc,jk,jb)   =     zaetr  *   zstbga*pt_diag%dpres_mc(jc,jk,jb)

            zaetr_top(jc,jb) = zaetr_bot
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
        IF (ntracer + ntracer_static >= io3) THEN
          DO jk = 1,nlev
            DO jc = i_startidx,i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,io3) = &
                & (amo3/amd) * (zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF
      ELSEIF ( irad_o3 == 6 ) THEN !ozone, but no aerosols
         ! 3-dimensional O3
        ! top level
        DO jc = i_startidx,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO
        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zpbot32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,jk+1,jb)))**3
            zo3_bot  (jc,jb) = prm_diag%vio3(jc,jb)* zpbot32(jc,jb)    &
              /( zpbot32(jc,jb) + zo3_hm(jc,jb))
            !O3 content
            zduo3(jc,jk,jb)= zo3_bot (jc,jb)-zo3_top (jc,jb)
            ! store previous bottom values in arrays for top of next layer
            zo3_top (jc,jb) = zo3_bot (jc,jb)
          ENDDO
        ENDDO
        IF (ntracer + ntracer_static >= io3) THEN
          DO jk = 1,nlev
            DO jc = i_startidx,i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,io3) = &
                & (amo3/amd) * (zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF
        zaeq1(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq2(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq3(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq4(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq5(i_startidx:i_endidx,:,jb) = 0.0_wp
      ELSEIF ( irad_aero == 5 ) THEN !aerosols, but no ozone:
        zduo3(i_startidx:i_endidx,:,jb) = 0.0_wp

        DO jk = 2, nlevp1
          DO jc = i_startidx,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = i_startidx,  & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! top level
        DO jc = i_startidx,i_endidx
          zaeqso   (jc,jb) = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,1)
          zaequo   (jc,jb) = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,1)
          zaeqdo   (jc,jb) = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zaeqsn         = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,jk+1)
            zaeqln         = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,jk+1)
            zaequn         = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,jk+1)
            zaeqdn         = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,jk+1)
            zaetr_bot      = zaetr_top(jc,jb) &
              & * ( MIN (1.0_wp, pt_diag%temp_ifc(jc,jk,jb)/pt_diag%temp_ifc(jc,jk+1,jb)) )**ztrpt

            zaetr          = SQRT(zaetr_bot*zaetr_top(jc,jb))
            zaeq1(jc,jk,jb)= (1._wp-zaetr) &
              & * (ztrbga* pt_diag%dpres_mc(jc,jk,jb)+zaeqln-zaeqlo(jc,jb)+zaeqdn-zaeqdo(jc,jb))
            zaeq2(jc,jk,jb)   = (1._wp-zaetr) * ( zaeqsn-zaeqso(jc,jb) )
            zaeq3(jc,jk,jb)   = (1._wp-zaetr) * ( zaequn-zaequo(jc,jb) )
            zaeq4(jc,jk,jb)   =     zaetr  *   zvobga*pt_diag%dpres_mc(jc,jk,jb)
            zaeq5(jc,jk,jb)   =     zaetr  *   zstbga*pt_diag%dpres_mc(jc,jk,jb)

            zaetr_top(jc,jb) = zaetr_bot
            zaeqso(jc,jb)    = zaeqsn
            zaeqlo(jc,jb)    = zaeqln
            zaequo(jc,jb)    = zaequn
            zaeqdo(jc,jb)    = zaeqdn
          ENDDO
        ENDDO
      ELSE !no aerosols and no ozone
        zaeq1(i_startidx:i_endidx,:,jb)= 0.0_wp
        zaeq2(i_startidx:i_endidx,:,jb)= 0.0_wp
        zaeq3(i_startidx:i_endidx,:,jb)= 0.0_wp
        zaeq4(i_startidx:i_endidx,:,jb)= 0.0_wp
        zaeq5(i_startidx:i_endidx,:,jb)= 0.0_wp
        zduo3(i_startidx:i_endidx,:,jb)= 0.0_wp
      ENDIF !irad_o3

    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

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
      &  CALL message('mo_nwp_rad_interface', 'Ritter-Geleyn radiation on reduced grid')

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
      zrg_albvisdif(nproma,nblks_par_c),          &
      zrg_alb_ther (nproma,nblks_par_c),          &
      zrg_pres_sfc (nproma,nblks_par_c),          &
      zrg_temp_ifc (nproma,nlevp1,nblks_par_c),   &
      zrg_dpres_mc (nproma,nlev  ,nblks_par_c),   &
      zrg_sqv      (nproma,nlev  ,nblks_par_c),   &
      zrg_duco2    (nproma,nlev  ,nblks_par_c),   &
      zrg_o3       (nproma,nlev  ,nblks_par_c),   &
      zrg_aeq1     (nproma,nlev  ,nblks_par_c),   &
      zrg_aeq2     (nproma,nlev  ,nblks_par_c),   &
      zrg_aeq3     (nproma,nlev  ,nblks_par_c),   &
      zrg_aeq4     (nproma,nlev  ,nblks_par_c),   &
      zrg_aeq5     (nproma,nlev  ,nblks_par_c),   &
      zrg_tot_cld  (nproma,nlev  ,nblks_par_c,4), &
      zrg_fls      (nproma,nlevp1,nblks_par_c),   &
      zrg_lwflxall (nproma,nlevp1,nblks_par_c),   &
      zrg_trsolall (nproma,nlevp1,nblks_par_c)    )


    rl_start = 1 ! SR radiation is not set up to handle boundaries of nested domains
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    
    ! *** this parallel section will be removed later once real data are
    !     are available as input for radiation ***
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,ist,zvege,zsnow,zsalb_snow,zsnow_alb)
    !
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zsqv (jc,jk,jb) = qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDDO

      ! geographical dependent thermal albedo
      alb_ther(i_startidx:i_endidx,jb) = 1._wp-ext_data%atm%emis_rad(i_startidx:i_endidx,jb)
      
      !------------------------------------------------------------------------------
      ! Calculation of surface albedo taking soil type,              
      ! vegetation and snow/ice conditions into account
      !------------------------------------------------------------------------------
      
      IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        
        DO jc = i_startidx,i_endidx

          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog_now%t_g(jc,jb) >= tmelt-1.7_wp ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF
          albvisdif(jc,jb) = csalb(ist)
          IF ( ext_data%atm%llsm_atm_c(jc,jb)) THEN
            ! ATTENTION: only valid, if nsfc_subs=1
            albvisdif(jc,jb) = csalb(ist) - rad_csalbw(ist)*lnd_prog_now%w_so(jc,1,jb,1)
          ENDIF  ! lsoil, llandmask
          
        ENDDO
        
      ELSE
        
        DO jc = i_startidx,i_endidx
          
          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog_now%t_g(jc,jb) >= tmelt-1.7_wp ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF
          albvisdif(jc,jb) = csalb(ist)
          
        ENDDO
      
      ENDIF

      !sea-ice model not yet implemented
!      IF (atm_phy_nwp_config(jg)%lseaice) THEN
!        DO jc = i_startidx,i_endidx
!          ! In case the sea ice model is used AND water point AND ice is present,
!          ! compute ice albedo for water points with an empirical formula taken from GME.
!          ! The ice albedo is the lower the warmer, and therefore wetter, the ice is.
!          ! Use ice temperature at time level nnow (2-time level scheme in sea ice model).
!!          IF ((.NOT. llandmask(i,j)) .AND. (h_ice(i,j,nnow) > 0.0_ireals))               &
!!            zalso(i,j) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(tmelt-t_ice(i,j,nnow)))) &
!!            * csalb(10)
!          IF (( .NOT. ext_data%atm%llsm_atm_c(jc,jb)) .AND. &
!            & (prm_diag%h_ice(i,j,nnow) > 0.0_wp))          &
!            albvisdif(jc,jb) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(tmelt-t_ice(i,j,nnow)))) &
!            * csalb(10)
!        ENDDO
!      ENDIF

      !lake model not yet implemented
!      IF (atm_phy_nwp_config(jg)%llake) THEN
!        DO jc = i_startidx,i_endidx
!          IF((ext_data%atm%depth_lk(jc,jb)      >  0.0_wp) .AND.    &
!            (prm_diag%h_ice(jc,jb) >= h_Ice_min_flk) ) THEN
!            !  In case the lake model FLake is used AND lake point AND ice is present,
!            !  compute ice albedo for lake points with an empirical formulation 
!            !  proposed by Mironov and Ritter (2004) for use in GME 
!            !  [ice_albedo=function(ice_surface_temperature)].
!            !  Use surface temperature at time level "nnow".
!
!            albvisdif(jc,jb) = EXP(-c_albice_MR*(tpl_T_f-t_s(i,j,nnow))/tpl_T_f)
!            albvisdif(jc,jb) = albedo_whiteice_ref * (1._ireals-zalso(i,j)) +      &
!              albedo_blueice_ref  * albvisdif(jc,jb)
!          ENDIF
!        ENDDO
!      ENDIF

      ! Snow cover and vegetation
      ! -------------------------

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN !soil model switched on
        DO jc = i_startidx,i_endidx
          zvege= 0.0_wp
          zsnow= 0.0_wp
          IF (ext_data%atm%llsm_atm_c(jc,jb)) THEN 
            ! consider effects of aging on solar snow albedo
            ! ATTENTION: only valid, if nsfc_subs=1            
            zsalb_snow = csalb_snow_min + &
              & lnd_diag%freshsnow(jc,jb,1)*(csalb_snow_max-csalb_snow_min)
            zsnow_alb = zsalb_snow*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb)) &
              + csalb_snow_fe * ext_data%atm%for_e(jc,jb)                       &
              + csalb_snow_fd * ext_data%atm%for_d(jc,jb)
            ! account for snow cover and plant cover and compute final solar
            ! snow albedo
            zvege = ext_data%atm%plcov_mx(jc,jb)
            ! ATTENTION: only valid, if nsfc_subs=1
            IF (lnd_prog_now%w_snow(jc,jb,1) > 0.0_wp) THEN
              ! ATTENTION: only valid, if nsfc_subs=1
              zsnow = MIN(1.0_wp, lnd_prog_now%w_snow(jc,jb,1) / cf_snow)
            ENDIF
            albvisdif(jc,jb) = zsnow * zsnow_alb +                               &
              (1.0_wp - zsnow) * (zvege * csalb_p + (1.0_wp - zvege) * albvisdif(jc,jb))
          ENDIF !   llandmask
        ENDDO
      ENDIF !inwp_surface == 1

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog_now%t_g(i_startidx:i_endidx,jb)

      ! CO2 (mixing ratio 353.9 ppm as vmr_co2)
      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zduco2(jc,jk,jb) = zqco2 * pt_diag%dpres_mc (jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

    CALL upscale_rad_input_rg(pt_patch, pt_par_patch, pt_par_grf_state,                    &
      & prm_diag%cosmu0, albvisdif, alb_ther, pt_diag%temp_ifc, pt_diag%dpres_mc,          &
      & prm_diag%tot_cld, zsqv ,zduco2, zduo3,                     &
      & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,pt_diag%pres_sfc,                                    &
      & zrg_cosmu0, zrg_albvisdif, zrg_alb_ther, zrg_temp_ifc, zrg_dpres_mc,               &
      & zrg_tot_cld, zrg_sqv ,zrg_duco2, zrg_o3,                                           &
      & zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5,zrg_pres_sfc     )

    rl_start = grf_ovlparea_start_c
    rl_end   = min_rlcell_int

    i_startblk = ptr_pp%cells%start_blk(rl_start,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(rl_end,i_chidx)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,lo_sol,losol,zi0) ,SCHEDULE(guided)
    !
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end, i_chidx)

      ! Switch off solar radiation calculations where sun is below horizon:
      WHERE ( zrg_cosmu0(i_startidx:i_endidx,jb) > 1.e-8_wp ) !zepmu0 )
        lo_sol(i_startidx:i_endidx) = .TRUE.
      ELSEWHERE
        lo_sol(i_startidx:i_endidx) = .FALSE.
      END WHERE
      losol = ANY(lo_sol(i_startidx:i_endidx))

!#ifdef __BOUNDCHECK
      CALL fesft ( &
                                         !  Input:
        & pti = zrg_temp_ifc (:,:,jb) , &! Temperature at layer boundaries
        & pdp = zrg_dpres_mc (:,:,jb), &! pressure thickness
        & pclc_in= zrg_tot_cld  (:,:,jb,icc) , &
        & pqv = zrg_tot_cld(:,:,jb,iqv), &! pt_prog_rcf%tracer(:,:,jb,iqv)
        & pqvs = zrg_sqv(:,:,jb), &!saturation water vapor
        & pqcwc = zrg_tot_cld    (:,:,jb,iqc) ,&
        & pqiwc = zrg_tot_cld    (:,:,jb,iqi) ,&
        & pduco2 = zrg_duco2 (:,:,jb), &! layer CO2 content
        & pduo3  = zrg_o3 (:,:,jb),&! layer O3 content
        & paeq1 = zrg_aeq1(:,:,jb), &
        & paeq2 = zrg_aeq2(:,:,jb),&
        & paeq3 = zrg_aeq3(:,:,jb),&
        & paeq4 = zrg_aeq4(:,:,jb),&
        & paeq5 = zrg_aeq5(:,:,jb),&
        & papre_in = zrg_pres_sfc (:,jb), & ! Surface pressure
        & psmu0 = zrg_cosmu0 (:,jb) , & ! Cosine of zenith angle
        & palso = zrg_albvisdif(:,jb), & ! solar surface albedo
        & palth = zrg_alb_ther(:,jb), & ! thermal surface albedo
        & psct = zsct, &! solar constant (at time of year)
        & kig1s = 1 ,&
        & kig1e = nproma , &
        & ki3s = 1, &
        & ki3e = nlev,&
        & ki1sc= i_startidx, &
        & ki1ec= i_endidx, &
        & lsolar = losol, &! control switch for solar calculations
        !          & lsolar = .TRUE., &! control switch for solar calculations
        & lthermal =.TRUE., &
        & lcrf = .FALSE., &! control switch for cloud-free calcul.
                                ! Output:
        & pflt  = zrg_lwflxall(:,:,jb) ,& ! Thermal radiative fluxes at each layer boundary
        & pfls  = zrg_fls  (:,:,jb)  &! solar radiative fluxes at each layer boundary
        & )
!#else
!        CALL fesft ( &
!                                !  Input:
!          & pti = zrg_temp_ifc (1,1,jb) , &! Temperature at layer boundaries
!          & pdp = zrg_dpres_mc (1,1,jb), &! pressure thickness
!          & pclc_in= zrg_tot_cld  (1,1,jb,icc) , &
!          & pqv = zrg_tot_cld(1,1,jb,iqv), &!pt_prog_rcf%tracer(1,1,jb,iqv)
!          & pqvs = zrg_sqv(1,1,jb), &!saturation water vapor
!          & pqcwc = zrg_tot_cld    (1,1,jb,iqc) ,&
!          & pqiwc = zrg_tot_cld    (1,1,jb,iqi) ,&
!          & pduco2 = zrg_duco2 (1,1,jb), &! layer CO2 content
!          & pduo3 = zrg_o3  (1,1,jb),&! layer O3 content
!          & paeq1 = zrg_aeq1(1,1,jb), &
!          & paeq2 = zrg_aeq2(1,1,jb),&
!          & paeq3 = zrg_aeq3(1,1,jb),&
!          & paeq4 = zrg_aeq4(1,1,jb),&
!          & paeq5 = zrg_aeq5(1,1,jb),&
!          & papre_in = zrg_pres_sfc (1,jb), & ! Surface pressure
!          & psmu0 = zrg_cosmu0 (1,jb) , & ! Cosine of zenith angle
!          & palso = zrg_albvisdif(1,jb), & ! solar surface albedo
!          & palth = zrg_alb_ther(1,jb), & ! thermal surface albedo
!          & psct = zsct, &! solar constant (at time of year)
!          & kig1s = 1 ,&
!          & kig1e = nproma , &
!          & ki3s = 1, &
!          & ki3e = nlev,&
!          & ki1sc= i_startidx, &
!          & ki1ec= i_endidx, &
!          & lsolar = losol, &! control switch for solar calculations
!          !          & lsolar = .TRUE., &! control switch for solar calculations
!          & lthermal =.TRUE., &
!          & lcrf = .FALSE., &! control switch for cloud-free calcul.
!                                ! Output:
!          & pflt  = zrg_lwflxall(1,1,jb) ,& ! Thermal radiative fluxes at each layer boundary
!          & pfls  = zrg_fls  (1,1,jb)  &! solar radiative fluxes at each layer boundary
!          & )
!        !          & pfltf = zfltf (1:i_endidx,1:2,jb),& !cloud_free thermal
!        !          & pflsf = zflsf (1:i_endidx,1:2,jb) ,& !cloud_free solar
!        !          & pflpar= zflpar(1:i_endidx,jb) ,& ! Photosynthetic active radiation
!        !          & pflsp = zflsp (1:i_endidx,jb) ,  & ! direct component of solar radiative flux
!        !          & pflsd = zflsd (1:i_endidx,jb)  ,& ! diffuse downward component of solar flux
!        !          & pflsu = zflsu (1:i_endidx,jb)  ) ! diffuse upward   component of solar flux
!#endif


      zi0 (i_startidx:i_endidx) = zrg_cosmu0(i_startidx:i_endidx,jb) * zsct
      ! compute sw transmissivity trsolall from sw fluxes
      DO jk = 1,nlevp1
        DO jc = i_startidx,i_endidx
          ! This is needed to avoid false synchronization errors
          IF (zrg_cosmu0(jc,jb) < 1.e-8_wp) THEN
            zrg_trsolall(jc,jk,jb) = 0._wp
          ELSE
            zrg_trsolall(jc,jk,jb) = zrg_fls(jc,jk,jb) / zi0(jc)
          ENDIF
        ENDDO
      ENDDO

    ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

    CALL downscale_rad_output_rg(pt_patch, pt_par_patch, pt_par_int_state,         &
      & pt_par_grf_state, zrg_lwflxall, zrg_trsolall, &
      & prm_diag%lwflxall, prm_diag%trsolall )

    DEALLOCATE (zrg_cosmu0,zrg_albvisdif,zrg_alb_ther,                      &
      & zrg_pres_sfc,zrg_temp_ifc,zrg_dpres_mc,zrg_sqv,zrg_duco2,zrg_o3,    &
      & zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5,                       &
      & zrg_tot_cld, zrg_fls,zrg_lwflxall, zrg_trsolall)

  END SUBROUTINE nwp_rg_radiation_reduced

END MODULE mo_nwp_rad_interface

