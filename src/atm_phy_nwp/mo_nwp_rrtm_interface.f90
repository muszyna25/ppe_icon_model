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

  USE mo_aerosol_util,         ONLY: zaea_rrtm,zaes_rrtm,zaeg_rrtm
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_datetime,             ONLY: t_datetime,  month2hour
  USE mo_exception,            ONLY: message,  finish, message_text
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi
  USE mo_grf_intp_data_strc,   ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc, io3_ape!, min_rlcell 
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_ovlparea_start_c
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_lrtm_par,             ONLY: jpband => nbndlw
  USE mo_nwp_lnd_state,        ONLY: t_lnd_prog
  USE mo_model_domain,         ONLY: t_patch, p_patch_local_parent
  USE mo_mpi,                  ONLY: my_process_is_mpi_seq
  USE mo_phys_nest_utilities,  ONLY: upscale_rad_input, downscale_rad_output, &
    &                                upscale_rad_input_rg, downscale_rad_output_rg
  USE mo_nonhydro_types,       ONLY: t_nh_diag
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag
  USE mo_o3_util,              ONLY: calc_o3_clim, calc_o3_gems
  USE mo_radiation,            ONLY: radiation, pre_radiation_nwp_steps
  USE mo_radiation_config,     ONLY: irad_o3, irad_aero, vmr_co2
  USE mo_radiation_rg,         ONLY: fesft
  USE mo_radiation_rg_par,     ONLY: aerdis
  USE mo_satad,                ONLY: qsat_rho
  USE mo_srtm_config,          ONLY: jpsw
  USE mo_sync,                 ONLY: global_max, global_min

  IMPLICIT NONE

  PRIVATE


  PUBLIC ::  nwp_rrtm_radiation, nwp_rrtm_radiation_reduced, nwp_rrtm_ozon_aerosol


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
  SUBROUTINE nwp_rrtm_ozon_aerosol ( p_sim_time, datetime, pt_patch, ext_data, &
    & pt_diag,prm_diag,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_datetime),            INTENT(in) :: datetime
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),       INTENT(inout) :: ext_data
    TYPE(t_nh_diag), TARGET, INTENT(in)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag

    REAL(wp), INTENT(out) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    ! for Ritter-Geleyn radiation:
    REAL(wp):: zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    ! for ozone:
    REAL(wp):: &
      & zptop32(nproma,pt_patch%nblks_c), &
      & zo3_hm (nproma,pt_patch%nblks_c), &
      & zo3_top (nproma,pt_patch%nblks_c), &
      & zpbot32(nproma,pt_patch%nblks_c), &
      & zo3_bot (nproma,pt_patch%nblks_c)
    ! for aerosols:
    REAL(wp):: &
      & zsign(nproma,pt_patch%nlevp1), &
      & zvdaes(nproma,pt_patch%nlevp1), &
      & zvdael(nproma,pt_patch%nlevp1), &
      & zvdaeu(nproma,pt_patch%nlevp1), &
      & zvdaed(nproma,pt_patch%nlevp1), &
      & zaetr_top(nproma,pt_patch%nblks_c), zaetr_bot, zaetr,       &
      & zaeqdo   (nproma,pt_patch%nblks_c), zaeqdn,                 &
      & zaequo   (nproma,pt_patch%nblks_c), zaequn,                 &
      & zaeqlo   (nproma,pt_patch%nblks_c), zaeqln,                 &
      & zaeqso   (nproma,pt_patch%nblks_c), zaeqsn,                 &
      ! for Tegen aerosol
      & z_aer_ss(nproma,pt_patch%nblks_c), &
      & z_aer_or(nproma,pt_patch%nblks_c), &
      & z_aer_bc(nproma,pt_patch%nblks_c), &
      & z_aer_su(nproma,pt_patch%nblks_c), &
      & z_aer_du(nproma,pt_patch%nblks_c), zw



    ! Local scalars:
    INTEGER:: jc,jk,jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index

    INTEGER:: imo1,imo2 !for Tegen aerosol time interpolation

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
    CASE(io3_ape)
      ! APE ozone: do nothing since everything is already
      ! set in mo_nwp_phy_init
    CASE (6)
      CALL calc_o3_clim(                             &
        & kbdim      = nproma,                       & ! in
        & jg         = jg,                           &
        & p_inc_rad  = atm_phy_nwp_config(jg)%dt_rad,& ! in
        & z_sim_time = p_sim_time,                   & ! in
        & pt_patch   = pt_patch,                     & ! in
        & zvio3      = prm_diag%vio3,                & !inout
        & zhmo3      = prm_diag%hmo3  )                !inout
    CASE (7)
      CALL calc_o3_gems(pt_patch,datetime,pt_diag,ext_data)
    END SELECT

    IF ( irad_aero == 6 ) CALL month2hour (datetime, imo1, imo2, zw )
    
    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_endidx, &
!$OMP       zsign,zvdaes, zvdael, zvdaeu, zvdaed, zaeqsn, zaeqln, zaequn,zaeqdn,zaetr_bot,zaetr )
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
          zaetr_top(jc,jb) = 1.0_wp
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
        DO jk = 1,nlev
          ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
          DO jc = 1,i_endidx
            ext_data%atm%o3(jc,jk,jb) = zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb)
          ENDDO
        ENDDO

      ELSEIF ( irad_o3 == 6 .AND. irad_aero == 6 ) THEN

        DO jc = 1,i_endidx

          z_aer_ss(jc,jb) = ext_data%atm_td%aer_ss(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_ss(jc,jb,imo2)   - ext_data%atm_td%aer_ss(jc,jb,imo1)   ) * zw
          z_aer_or(jc,jb) = ext_data%atm_td%aer_org(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_org(jc,jb,imo2)  - ext_data%atm_td%aer_org(jc,jb,imo1)  ) * zw
          z_aer_bc(jc,jb) = ext_data%atm_td%aer_bc(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_bc(jc,jb,imo2)   - ext_data%atm_td%aer_bc(jc,jb,imo1)   ) * zw
          z_aer_su(jc,jb) = ext_data%atm_td%aer_so4(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_so4(jc,jb,imo2)  - ext_data%atm_td%aer_so4(jc,jb,imo1)  ) * zw
          z_aer_du(jc,jb) = ext_data%atm_td%aer_dust(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_dust(jc,jb,imo2) - ext_data%atm_td%aer_dust(jc,jb,imo1) ) * zw

        ENDDO

        DO jk = 2, nlevp1
          DO jc = 1,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO
        
        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = 1,           & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! 3-dimensional O3
        ! top level
        DO jc = 1,i_endidx
          zptop32  (jc,jb) = (SQRT(pt_diag%pres_ifc(jc,1,jb)))**3
          zo3_hm   (jc,jb) = (SQRT(prm_diag%hmo3(jc,jb)))**3
          zaeqso   (jc,jb) = z_aer_ss(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = ( z_aer_or(jc,jb)+z_aer_su(jc,jb) )*zvdael(jc,1)
          zaequo   (jc,jb) = z_aer_bc(jc,jb)              *zvdaeu(jc,1)
          zaeqdo   (jc,jb) =  z_aer_du(jc,jb)              *zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = 1,i_endidx
            zaeqsn         =  z_aer_ss(jc,jb)                  * zvdaes(jc,jk+1)
            zaeqln         = (z_aer_or(jc,jb)+z_aer_su(jc,jb)) * zvdael(jc,jk+1)
            zaequn         = z_aer_bc(jc,jb)                   * zvdaeu(jc,jk+1)
            zaeqdn         =  z_aer_du(jc,jb)                  * zvdaed(jc,jk+1)
            zaetr_bot      = zaetr_top(jc,jb) &
              & * ( MIN (1.0_wp, pt_diag%temp_ifc(jc,jk,jb)/pt_diag%temp_ifc(jc,jk+1,jb)) )**ztrpt

            zaetr          = SQRT(zaetr_bot*zaetr_top(jc,jb))
            zaeq1(jc,jk,jb)= (1.0_wp-zaetr)*( ztrbga*pt_diag%dpres_mc(jc,jk,jb) &
              &            + zaeqln - zaeqlo(jc,jb) )
            zaeq2(jc,jk,jb)   = (1.0_wp-zaetr)*(zaeqsn-zaeqso(jc,jb))
            zaeq3(jc,jk,jb)   = (1.0_wp-zaetr)*(zaeqdn-zaeqdo(jc,jb))
            zaeq4(jc,jk,jb)   = (1.0_wp-zaetr)*(zaequn-zaequo(jc,jb))
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

        DO jk = 1,nlev
          ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
          DO jc = 1,i_endidx
            ext_data%atm%o3(jc,jk,jb) = zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb)
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
            ext_data%atm%o3(jc,jk,jb) = zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb)
          ENDDO
        ENDDO
        zaeq1(1:i_endidx,:,jb) = 0.0_wp
        zaeq2(1:i_endidx,:,jb) = 0.0_wp
        zaeq3(1:i_endidx,:,jb) = 0.0_wp
        zaeq4(1:i_endidx,:,jb) = 0.0_wp
        zaeq5(1:i_endidx,:,jb) = 0.0_wp
        
      ELSEIF ( irad_aero == 5 ) THEN !aerosols, but no ozone:

        DO jk = 2, nlevp1
          DO jc = 1,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = 1,           & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! top level
        DO jc = 1,i_endidx
          zaeqso   (jc,jb) = zaeops*prm_diag%aersea(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = zaeopl*prm_diag%aerlan(jc,jb)*zvdael(jc,1)
          zaequo   (jc,jb) = zaeopu*prm_diag%aerurb(jc,jb)*zvdaeu(jc,1)
          zaeqdo   (jc,jb) = zaeopd*prm_diag%aerdes(jc,jb)*zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = 1,i_endidx
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

      ELSEIF (irad_aero == 6 ) THEN !aerosols, but other ozone:

        DO jc = 1,i_endidx

          z_aer_ss(jc,jb) = ext_data%atm_td%aer_ss(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_ss(jc,jb,imo2)   - ext_data%atm_td%aer_ss(jc,jb,imo1)   ) * zw
          z_aer_or(jc,jb) = ext_data%atm_td%aer_org(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_org(jc,jb,imo2)  - ext_data%atm_td%aer_org(jc,jb,imo1)  ) * zw
          z_aer_bc(jc,jb) = ext_data%atm_td%aer_bc(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_bc(jc,jb,imo2)   - ext_data%atm_td%aer_bc(jc,jb,imo1)   ) * zw
          z_aer_su(jc,jb) = ext_data%atm_td%aer_so4(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_so4(jc,jb,imo2)  - ext_data%atm_td%aer_so4(jc,jb,imo1)  ) * zw
          z_aer_du(jc,jb) = ext_data%atm_td%aer_dust(jc,jb,imo1) + &
            & ( ext_data%atm_td%aer_dust(jc,jb,imo2) - ext_data%atm_td%aer_dust(jc,jb,imo1) ) * zw

        ENDDO

        DO jk = 2, nlevp1
          DO jc = 1,i_endidx
            zsign(jc,jk) = pt_diag%pres_ifc(jc,jk,jb) / 101325._wp
          ENDDO
        ENDDO

        ! The routine aerdis is called to recieve some parameters for the vertical
        ! distribution of background aerosol.
        CALL aerdis ( &
          & kbdim  = nproma,      & !in
          & jcs    = 1,           & !in
          & jce    = i_endidx,    & !in
          & klevp1 = nlevp1,      & !in
          & petah  = zsign(1,1),  & !in
          & pvdaes = zvdaes(1,1), & !out
          & pvdael = zvdael(1,1), & !out
          & pvdaeu = zvdaeu(1,1), & !out
          & pvdaed = zvdaed(1,1) )  !out

        ! top level
        DO jc = 1,i_endidx
          zaeqso   (jc,jb) = z_aer_ss(jc,jb)              *zvdaes(jc,1)
          zaeqlo   (jc,jb) = ( z_aer_or(jc,jb)+z_aer_su(jc,jb) )*zvdael(jc,1)
          zaequo   (jc,jb) =  z_aer_bc(jc,jb)              *zvdaeu(jc,1)
          zaeqdo   (jc,jb) =  z_aer_du(jc,jb)              *zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = 1,i_endidx
            zaeqsn         =  z_aer_ss(jc,jb)                   * zvdaes(jc,jk+1)
            zaeqln         =  (z_aer_or(jc,jb)+z_aer_su(jc,jb)) * zvdael(jc,jk+1)
            zaequn         =  z_aer_bc(jc,jb)                   * zvdaeu(jc,jk+1)
            zaeqdn         =  z_aer_du(jc,jb)                   * zvdaed(jc,jk+1)
            zaetr_bot      = zaetr_top(jc,jb) &
              & * ( MIN (1.0_wp, pt_diag%temp_ifc(jc,jk,jb)/pt_diag%temp_ifc(jc,jk+1,jb)) )**ztrpt

            zaetr          = SQRT(zaetr_bot*zaetr_top(jc,jb))
            zaeq1(jc,jk,jb)= (1.0_wp-zaetr)*( ztrbga*pt_diag%dpres_mc(jc,jk,jb) &
              &            + zaeqln - zaeqlo(jc,jb) )
            zaeq2(jc,jk,jb)   = (1._wp-zaetr) * ( zaeqsn-zaeqso(jc,jb) )
            zaeq3(jc,jk,jb)   = (1.0_wp-zaetr)*(zaeqdn-zaeqdo(jc,jb))
            zaeq4(jc,jk,jb)   = (1.0_wp-zaetr)*(zaequn-zaequo(jc,jb))
            zaeq5(jc,jk,jb)   =     zaetr  *   zstbga*pt_diag%dpres_mc(jc,jk,jb)

            zaetr_top(jc,jb) = zaetr_bot
            zaeqso(jc,jb)    = zaeqsn
            zaeqlo(jc,jb)    = zaeqln
            zaequo(jc,jb)    = zaequn
            zaeqdo(jc,jb)    = zaeqdn

          ENDDO
        ENDDO       
        
      ELSE !no aerosols

        zaeq1(1:i_endidx,:,jb) = 0.0_wp
        zaeq2(1:i_endidx,:,jb) = 0.0_wp
        zaeq3(1:i_endidx,:,jb) = 0.0_wp
        zaeq4(1:i_endidx,:,jb) = 0.0_wp
        zaeq5(1:i_endidx,:,jb) = 0.0_wp

      ENDIF !irad_o3

    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE nwp_rrtm_ozon_aerosol
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rrtm_radiation ( p_sim_time,pt_patch, &
    & ext_data,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,pt_diag,prm_diag,lnd_prog )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),INTENT(in):: ext_data

    REAL(wp), INTENT(in) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)
    
    TYPE(t_nh_diag), TARGET, INTENT(in)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog

    REAL(wp):: albvisdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdif     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: aclcov        (nproma,pt_patch%nblks_c) !<

    INTEGER :: itype(nproma)   !< type of convection

    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    INTEGER:: jc,jb
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
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, &
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
#ifdef __xlC__
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,itype)
#else
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,itype),SCHEDULE(guided)
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)


      ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
      itype(1:i_endidx) = 0 !INT(field%rtype(1:i_endidx,jb))


      ! It may happen that an MPI patch contains only nest boundary points
      ! In this case, no action is needed
      IF (i_startidx > i_endidx) CYCLE


      !Calculate direct albedo from diffuse albedo and solar zenith angle
      !formula as in Ritter-Geleyn's fesft
      DO jc = 1,i_endidx
        albvisdir(jc,jb) =  ( 1.0_wp                                                           &
          &  + 0.5_wp * (prm_diag%cosmu0(jc,jb) * (1.0_wp/prm_diag%albvisdif(jc,jb) - 1.0_wp))) &
          & / (1.0_wp + (prm_diag%cosmu0(jc,jb) * (1.0_wp/prm_diag%albvisdif(jc,jb) - 1.0_wp)))**2
      ENDDO

      ! no distiction between vis and nir albedo
      albnirdir(1:i_endidx,jb) = albvisdir(1:i_endidx,jb)
      albnirdif(1:i_endidx,jb) = prm_diag%albvisdif(1:i_endidx,jb)

      prm_diag%tsfctrad(1:i_endidx,jb) = lnd_prog%t_g(1:i_endidx,jb)


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
        & alb_vis_dif=prm_diag%albvisdif(:,jb),&!< in surface albedo for visible range, diffuse
        & alb_nir_dif=albnirdif        (:,jb) ,&!< in surface albedo for near IR range, diffuse
        & emis_rad=ext_data%atm%emis_rad(:,jb),&!< in longwave surface emissivity
        & tk_sfc     =prm_diag%tsfctrad(:,jb) ,&!< in surface temperature
                              !
                              ! atmosphere: pressure, tracer mixing ratios and temperature
        & pp_hl      =pt_diag%pres_ifc  (:,:,jb)     ,&!< in  pres at half levels at t-dt [Pa]
        & pp_fl      =pt_diag%pres      (:,:,jb)     ,&!< in  pres at full levels at t-dt [Pa]
        & tk_fl      =pt_diag%temp      (:,:,jb)     ,&!< in  temperature at full level at t-dt
        & qm_vap     =prm_diag%tot_cld  (:,:,jb,iqv) ,&!< in  water vapor mass mix ratio at t-dt
        & qm_liq     =prm_diag%tot_cld  (:,:,jb,iqc) ,&!< in cloud water mass mix ratio at t-dt
        & qm_ice     =prm_diag%tot_cld  (:,:,jb,iqi) ,&!< in cloud ice mass mixing ratio at t-dt
        & qm_o3      =ext_data%atm%o3   (:,:,jb)     ,&!< in o3 mass mixing ratio at t-dt
        & cdnc       =prm_diag%acdnc    (:,:,jb)     ,&!< in  cloud droplet numb conc. [1/m**3]
        & cld_frc    =prm_diag%tot_cld  (:,:,jb,icc) ,&!< in  cloud fraction [m2/m2]
        & zaeq1      = zaeq1(:,:,jb)                 ,&!< in aerosol continental
        & zaeq2      = zaeq2(:,:,jb)                 ,&!< in aerosol maritime
        & zaeq3      = zaeq3(:,:,jb)                 ,&!< in aerosol urban
        & zaeq4      = zaeq4(:,:,jb)                 ,&!< in aerosol volcano ashes
        & zaeq5      = zaeq5(:,:,jb)                 ,&!< in aerosol stratospheric background
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
  SUBROUTINE nwp_rrtm_radiation_reduced ( p_sim_time,pt_patch,pt_par_patch,               &
    & pt_par_int_state, pt_par_grf_state,ext_data,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5, &
    & pt_diag,prm_diag,lnd_prog )

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
    REAL(wp),             INTENT(in) ::               &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)    

    TYPE(t_nh_diag), TARGET,    INTENT(inout):: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog

    REAL(wp):: albvisdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdir     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: albnirdif     (nproma,pt_patch%nblks_c) !<
    REAL(wp):: aclcov        (nproma,pt_patch%nblks_c) !<
    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_land  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_glac  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_emis_rad (:,:)
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
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq1(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq2(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq3(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq4(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq5(:,:,:)
    ! Output fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aclcov   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxclr (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolclr (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolall (:,:,:)
    ! Pointer to parent patach or local parent patch for reduced grid
    TYPE(t_patch), POINTER       :: ptr_pp

    INTEGER :: itype(nproma)   !< type of convection

    ! Variables for debug output
    REAL(wp) :: max_albvisdir, min_albvisdir, max_albvisdif, min_albvisdif, &
                max_tsfc, min_tsfc, max_psfc, min_psfc

    REAL(wp), DIMENSION(pt_patch%nlevp1) :: max_pres_ifc, max_pres, max_temp, max_acdnc, &
        max_qv, max_qc, max_qi, max_cc, min_pres_ifc, min_pres, min_temp, min_acdnc, &
        min_qv, min_qc, min_qi, min_cc

    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    INTEGER:: jc,jk,jb
    INTEGER:: jg                     !domain id
    INTEGER:: nlev, nlevp1, nlev_rg  !< number of full and half levels
    INTEGER:: nblks_par_c            !nblks for reduced grid

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
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, &
      & p_sim_time   = p_sim_time,                        &
      & pt_patch     = pt_patch,                          &
      & zsmu0        = prm_diag%cosmu0(:,:),              &
      & zsct         = zsct )

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------


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
      ELSE ! Nested domain with MPI parallelization
        ptr_pp      => p_patch_local_parent(jg)
        nblks_par_c =  ptr_pp%nblks_c
      ENDIF

      ! Add extra layer for atmosphere above model top if requested
      IF (atm_phy_nwp_config(jg)%latm_above_top) THEN
        nlev_rg = nlev + 1
      ELSE
        nlev_rg = nlev
      ENDIF

      ALLOCATE (zrg_cosmu0   (nproma,nblks_par_c),     &
        zrg_fr_land  (nproma,nblks_par_c),             &
        zrg_fr_glac  (nproma,nblks_par_c),             &
        zrg_emis_rad (nproma,nblks_par_c),             &
        zrg_albvisdir(nproma,nblks_par_c),             &
        zrg_albnirdir(nproma,nblks_par_c),             &
        zrg_albvisdif(nproma,nblks_par_c),             &
        zrg_albnirdif(nproma,nblks_par_c),             &
        zrg_tsfc     (nproma,nblks_par_c),             &
        zrg_pres_ifc (nproma,nlev_rg+1,nblks_par_c),   &
        zrg_pres     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_temp     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_o3       (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq1     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq2     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq3     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq4     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq5     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_acdnc    (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_tot_cld  (nproma,nlev_rg  ,nblks_par_c,4), &
        zrg_aclcov   (nproma,          nblks_par_c),   &
        zrg_lwflxclr (nproma,nlev_rg+1,nblks_par_c),   &
        zrg_lwflxall (nproma,nlev_rg+1,nblks_par_c),   &
        zrg_trsolclr (nproma,nlev_rg+1,nblks_par_c),   &
        zrg_trsolall (nproma,nlev_rg+1,nblks_par_c)    )

      rl_start = 1 ! SR radiation is not set up to handle boundaries of nested domains
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


      ! *** this parallel section will be removed later once real data are
      !     are available as input for radiation ***
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          &                       i_startidx, i_endidx, rl_start, rl_end)


      ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM

      !Calculate direct albedo from diffuse albedo and solar zenith angle
      !formula as in Ritter-Geleyn's fesft
      DO jc = 1,i_endidx
        albvisdir(jc,jb) = ( 1.0_wp                                                             &
          &  + 0.5_wp * (prm_diag%cosmu0(jc,jb) * (1.0_wp/prm_diag%albvisdif(jc,jb) - 1.0_wp))) &
          & / (1.0_wp + (prm_diag%cosmu0(jc,jb) * (1.0_wp/prm_diag%albvisdif(jc,jb) - 1.0_wp)))**2
      ENDDO

      ! no distiction between vis and nir albedo
      albnirdir(1:i_endidx,jb) = albvisdir(1:i_endidx,jb)
      albnirdif(1:i_endidx,jb) = prm_diag%albvisdif(1:i_endidx,jb)

      prm_diag%tsfctrad(1:i_endidx,jb) = lnd_prog%t_g(1:i_endidx,jb)

      ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

      CALL upscale_rad_input(pt_patch, pt_par_patch, pt_par_grf_state,  &
        & nlev_rg, ext_data%atm%fr_land_smt, ext_data%atm%fr_glac_smt,  &
        & ext_data%atm%emis_rad,                                        &
        & prm_diag%cosmu0, albvisdir, albnirdir, prm_diag%albvisdif,    &
        & albnirdif, prm_diag%tsfctrad, pt_diag%pres_ifc,               &
        & pt_diag%pres, pt_diag%temp,prm_diag%acdnc, prm_diag%tot_cld,  &
        & ext_data%atm%o3(:,:,:),                                       &
        & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                            &
        & zrg_fr_land, zrg_fr_glac, zrg_emis_rad,                       &
        & zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,      &
        & zrg_albnirdif, zrg_tsfc, zrg_pres_ifc, zrg_pres, zrg_temp,    &
        & zrg_acdnc, zrg_tot_cld, zrg_o3,                               &
        & zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_aeq4, zrg_aeq5 )


      rl_start = grf_ovlparea_start_c
      rl_end   = min_rlcell_int

      i_startblk = ptr_pp%cells%start_blk(rl_start,i_chidx)
      i_endblk   = ptr_pp%cells%end_blk(rl_end,i_chidx)

      ! Debug output of radiation input fields
      IF (msg_level >= 16) THEN
        max_albvisdir = 0._wp
        min_albvisdir = 1.e10_wp
        max_albvisdif = 0._wp
        min_albvisdif = 1.e10_wp
        max_tsfc = 0._wp
        min_tsfc = 1.e10_wp
        max_psfc = 0._wp
        min_psfc = 1.e10_wp
        max_pres_ifc = 0._wp
        max_pres    = 0._wp
        max_temp     = 0._wp
        max_acdnc    = 0._wp
        max_qv = 0._wp
        max_qc = 0._wp
        max_qi = 0._wp
        max_cc  = 0._wp
        min_pres_ifc = 1.e10_wp
        min_pres    = 1.e10_wp
        min_temp     = 1.e10_wp
        min_acdnc    = 1.e10_wp
        min_qv = 1.e10_wp
        min_qc = 1.e10_wp
        min_qi = 1.e10_wp
        min_cc  = 1.e10_wp

        DO jb = i_startblk, i_endblk

         CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end, i_chidx)

         max_albvisdir = MAX(max_albvisdir,MAXVAL(zrg_albvisdir(i_startidx:i_endidx,jb)))
         min_albvisdir = MIN(min_albvisdir,MINVAL(zrg_albvisdir(i_startidx:i_endidx,jb)))
         max_albvisdif = MAX(max_albvisdif,MAXVAL(zrg_albvisdif(i_startidx:i_endidx,jb)))
         min_albvisdif = MIN(min_albvisdif,MINVAL(zrg_albvisdif(i_startidx:i_endidx,jb)))
         max_tsfc = MAX(max_tsfc,MAXVAL(zrg_tsfc(i_startidx:i_endidx,jb)))
         min_tsfc = MIN(min_tsfc,MINVAL(zrg_tsfc(i_startidx:i_endidx,jb)))
         max_psfc = MAX(max_psfc,MAXVAL(zrg_pres_ifc(i_startidx:i_endidx,nlev_rg+1,jb)))
         min_psfc = MIN(min_psfc,MINVAL(zrg_pres_ifc(i_startidx:i_endidx,nlev_rg+1,jb)))
         DO jk = 1, nlev_rg
          max_pres_ifc(jk) = MAX(max_pres_ifc(jk),MAXVAL(zrg_pres_ifc(i_startidx:i_endidx,jk,jb)))
          max_pres(jk)    = MAX(max_pres(jk),MAXVAL(zrg_pres     (i_startidx:i_endidx,jk,jb)))
          max_temp(jk)     = MAX(max_temp(jk),MAXVAL(zrg_temp     (i_startidx:i_endidx,jk,jb)))
          max_acdnc(jk)    = MAX(max_acdnc(jk),MAXVAL(zrg_acdnc    (i_startidx:i_endidx,jk,jb)))
          max_qv(jk) = MAX(max_qv(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqv)))
          max_qc(jk) = MAX(max_qc(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqc)))
          max_qi(jk) = MAX(max_qi(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqi)))
          max_cc(jk)  = MAX(max_cc(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,icc)))
          min_pres_ifc(jk) = MIN(min_pres_ifc(jk),MINVAL(zrg_pres_ifc(i_startidx:i_endidx,jk,jb)))
          min_pres(jk)    = MIN(min_pres(jk),MINVAL(zrg_pres     (i_startidx:i_endidx,jk,jb)))
          min_temp(jk)     = MIN(min_temp(jk),MINVAL(zrg_temp     (i_startidx:i_endidx,jk,jb)))
          min_acdnc(jk)    = MIN(min_acdnc(jk),MINVAL(zrg_acdnc    (i_startidx:i_endidx,jk,jb)))
          min_qv(jk) = MIN(min_qv(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqv)))
          min_qc(jk) = MIN(min_qc(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqc)))
          min_qi(jk) = MIN(min_qi(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqi)))
          min_cc(jk)  = MIN(min_cc(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,icc)))
         ENDDO
        ENDDO ! blocks

        max_albvisdir = global_max(max_albvisdir)
        min_albvisdir = global_min(min_albvisdir)
        max_albvisdif = global_max(max_albvisdif)
        min_albvisdif = global_min(min_albvisdif)
        max_tsfc = global_max(max_tsfc)
        min_tsfc = global_min(min_tsfc)
        max_psfc = global_max(max_psfc)
        min_psfc = global_min(min_psfc)
        max_pres_ifc = global_max(max_pres_ifc)
        max_pres    = global_max(max_pres)
        max_temp     = global_max(max_temp)
        max_acdnc    = global_max(max_acdnc)
        max_qv = global_max(max_qv)
        max_qc = global_max(max_qc)
        max_qi = global_max(max_qi)
        max_cc  = global_max(max_cc)
        min_pres_ifc = global_min(min_pres_ifc)
        min_pres    = global_min(min_pres)
        min_temp     = global_min(min_temp)
        min_acdnc    = global_min(min_acdnc)
        min_qv = global_min(min_qv)
        min_qc = global_min(min_qc)
        min_qi = global_min(min_qi)
        min_cc  = global_min(min_cc)


        WRITE(message_text,'(a,4f12.8)') 'max/min alb = ', max_albvisdir, min_albvisdir, &
          max_albvisdif, min_albvisdif
        CALL message('nwp_nh_interface: ', TRIM(message_text))

        WRITE(message_text,'(a,2f10.3,2f10.2)') 'max/min sfc temp/pres = ', max_tsfc, min_tsfc, &
          max_psfc, min_psfc
        CALL message('nwp_nh_interface: ', TRIM(message_text))

        WRITE(message_text,'(a)') 'max/min pres_ifc, pres, temp, acdnc'
        CALL message('nwp_nh_interface: ', TRIM(message_text))

        DO jk = 1, nlev_rg
          WRITE(message_text,'(i4,4f10.2,2f10.3,2f12.1)') jk,max_pres_ifc(jk), min_pres_ifc(jk), &
            max_pres(jk), min_pres(jk), max_temp(jk), min_temp(jk), max_acdnc(jk), min_acdnc(jk)
          CALL message('nwp_nh_interface: ', TRIM(message_text))
        ENDDO

        WRITE(message_text,'(a)') 'max/min QV, QC, QI, CC'
        CALL message('nwp_nh_interface: ', TRIM(message_text))

        DO jk = 1, nlev_rg
          WRITE(message_text,'(i4,8e13.5)') jk,max_qv(jk), min_qv(jk), max_qc(jk), min_qc(jk), &
             max_qi(jk), min_qi(jk), max_cc(jk), min_cc(jk)
          CALL message('nwp_nh_interface: ', TRIM(message_text))
        ENDDO

      ENDIF ! msg_level >= 16

!$OMP PARALLEL
#ifdef __xlC__
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,itype)
#else
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,itype),SCHEDULE(guided)      
#endif
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
          &                         i_startidx, i_endidx, rl_start, rl_end, i_chidx)

        itype(:) = 0

        ! It may happen that an MPI patch contains only nest boundary points
        ! In this case, no action is needed
        IF (i_startidx > i_endidx) CYCLE

        ! Unfortunately, the coding of SR radiation is not compatible with the presence
        ! of nested domains. Therefore, the normally unused elements of the first block
        ! need to be filled with dummy values
        IF (jg > 1 .AND. jb == i_startblk) THEN
          zrg_fr_land   (1:i_startidx-1,jb) = zrg_fr_land   (i_startidx,jb)
          zrg_fr_glac   (1:i_startidx-1,jb) = zrg_fr_glac   (i_startidx,jb)
          zrg_emis_rad  (1:i_startidx-1,jb) = zrg_emis_rad  (i_startidx,jb)
          zrg_cosmu0    (1:i_startidx-1,jb) = zrg_cosmu0    (i_startidx,jb)
          zrg_albvisdir (1:i_startidx-1,jb) = zrg_albvisdir (i_startidx,jb)
          zrg_albnirdir (1:i_startidx-1,jb) = zrg_albnirdir (i_startidx,jb)
          zrg_albvisdif (1:i_startidx-1,jb) = zrg_albvisdif (i_startidx,jb)
          zrg_albnirdif (1:i_startidx-1,jb) = zrg_albnirdif (i_startidx,jb)
          zrg_tsfc      (1:i_startidx-1,jb) = zrg_tsfc      (i_startidx,jb)
          zrg_pres_ifc (1:i_startidx-1,nlev_rg+1,jb) = zrg_pres_ifc (i_startidx,nlev_rg+1,jb)
          DO jk = 1, nlev_rg
            zrg_pres_ifc (1:i_startidx-1,jk,jb) = zrg_pres_ifc (i_startidx,jk,jb)
            zrg_pres     (1:i_startidx-1,jk,jb) = zrg_pres     (i_startidx,jk,jb)
            zrg_temp     (1:i_startidx-1,jk,jb) = zrg_temp     (i_startidx,jk,jb)
            zrg_o3       (1:i_startidx-1,jk,jb) = zrg_o3       (i_startidx,jk,jb)
            zrg_aeq1     (1:i_startidx-1,jk,jb) = zrg_aeq1     (i_startidx,jk,jb)
            zrg_aeq2     (1:i_startidx-1,jk,jb) = zrg_aeq2     (i_startidx,jk,jb)
            zrg_aeq3     (1:i_startidx-1,jk,jb) = zrg_aeq3     (i_startidx,jk,jb)
            zrg_aeq4     (1:i_startidx-1,jk,jb) = zrg_aeq4     (i_startidx,jk,jb)
            zrg_aeq5     (1:i_startidx-1,jk,jb) = zrg_aeq5     (i_startidx,jk,jb)
            zrg_acdnc    (1:i_startidx-1,jk,jb) = zrg_acdnc    (i_startidx,jk,jb)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqv) = zrg_tot_cld(i_startidx,jk,jb,iqv)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqc) = zrg_tot_cld(i_startidx,jk,jb,iqc)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqi) = zrg_tot_cld(i_startidx,jk,jb,iqi)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,icc) = zrg_tot_cld(i_startidx,jk,jb,icc)
          ENDDO
        ENDIF

        CALL radiation(               &
                                !
                                ! input
                                ! -----
                                !
                                ! indices and dimensions
          & jce         =i_endidx            ,&!< in  end   index for loop over block
          & kbdim       =nproma              ,&!< in  dimension of block over cells
          & klev        =nlev_rg             ,&!< in  number of full levels = number of layers
          & klevp1      =nlev_rg+1           ,&!< in  number of half levels = number of layer ifcs
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
          & emis_rad   =zrg_emis_rad(:,jb)   ,&!< in longwave surface emissivity
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
          & zaeq1      = zrg_aeq1(:,:,jb)       ,&!< in aerosol continental
          & zaeq2      = zrg_aeq2(:,:,jb)       ,&!< in aerosol maritime
          & zaeq3      = zrg_aeq3(:,:,jb)       ,&!< in aerosol urban
          & zaeq4      = zrg_aeq4(:,:,jb)       ,&!< in aerosol volcano ashes
          & zaeq5      = zrg_aeq5(:,:,jb)       ,&!< in aerosol stratospheric background
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

      CALL downscale_rad_output(pt_patch, pt_par_patch, pt_par_int_state,                  &
        & pt_par_grf_state, nlev_rg, zrg_aclcov, zrg_lwflxclr, zrg_lwflxall, zrg_trsolclr, &
        & zrg_trsolall, zrg_tsfc, zrg_albvisdif, zrg_emis_rad, zrg_cosmu0, zrg_tot_cld,    &
        & zrg_pres_ifc, prm_diag%tsfctrad, prm_diag%albvisdif, aclcov, prm_diag%lwflxclr, &
        & prm_diag%lwflxall, prm_diag%trsolclr, prm_diag%trsolall )

      DEALLOCATE (zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif, zrg_albnirdif, &
        zrg_tsfc, zrg_pres_ifc, zrg_pres, zrg_temp, zrg_o3,                               &
        zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5, zrg_acdnc, zrg_tot_cld,             &
        zrg_aclcov, zrg_lwflxclr, zrg_lwflxall, zrg_trsolclr, zrg_trsolall,               &
        zrg_fr_land,zrg_fr_glac,zrg_emis_rad)
      
  END SUBROUTINE nwp_rrtm_radiation_reduced


END MODULE mo_nwp_rrtm_interface


