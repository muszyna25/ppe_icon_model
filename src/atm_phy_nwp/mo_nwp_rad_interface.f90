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
  USE mo_datetime,             ONLY: t_datetime,  month2hour
  USE mo_exception,            ONLY: message,  finish !message_tex
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run, parallel_radiation_mode

  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi
  USE mo_grf_intp_data_strc,   ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc, io3_ape!, min_rlcell 
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_ovlparea_start_c
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch, p_patch_local_parent
  USE mo_mpi,                  ONLY: my_process_is_mpi_seq
  USE mo_phys_nest_utilities,  ONLY: upscale_rad_input_rg, downscale_rad_output_rg
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_o3_util,              ONLY: calc_o3_clim,calc_o3_gems
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_radiation_config,     ONLY: irad_o3, irad_aero, vmr_co2
  USE mo_radiation_rg,         ONLY: fesft
  USE mo_radiation_rg_par,     ONLY: aerdis
  USE mo_satad,                ONLY: qsat_rho
!  USE mo_sync,                 ONLY: SYNC_C, sync_patch_array

  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation, &
   &  nwp_rrtm_radiation_reduced, nwp_rrtm_radiation_repartition, nwp_rrtm_ozon_aerosol
!   USE mo_nwp_mpiomp_rrtm_interface, ONLY: nwp_omp_rrtm_interface
  USE mo_albedo,               ONLY: sfc_albedo

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_radiation
  

  CHARACTER(len=*), PARAMETER:: version = '$Id$'

 CONTAINS
  
  !---------------------------------------------------------------------------------------
  !>
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1) or Ritter-Geleyn (2).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time, datetime, pt_patch,pt_par_patch, &
    & pt_par_int_state,pt_par_grf_state,ext_data,lnd_diag,pt_prog,pt_diag,prm_diag, &
    & lnd_prog, wtr_prog )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'
    
    LOGICAL, INTENT(in)         :: lredgrid        !< use reduced grid for radiation

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_datetime),            INTENT(in) :: datetime
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_int_state),    TARGET,INTENT(in):: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(in) :: pt_par_grf_state  !< grid refinement state
    TYPE(t_external_data),INTENT(inout):: ext_data
    TYPE(t_lnd_diag),     INTENT(in):: lnd_diag      !< diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog   ! time level new
    TYPE(t_wtr_prog),           INTENT(   in):: wtr_prog   ! time level new

    REAL(wp) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    
    INTEGER :: jg

    jg = pt_patch%id


    ! Compute tile-based and aggregated surface-albedo
    CALL sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)

    IF (atm_phy_nwp_config(jg)%inwp_radiation == 1 ) THEN
       
      CALL nwp_rrtm_ozon_aerosol ( p_sim_time, datetime, pt_patch, ext_data, &
        & pt_diag,prm_diag,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )
    
      IF ( .NOT. lredgrid ) THEN

        SELECT CASE(parallel_radiation_mode(pt_patch%id))
        CASE(1) 
          CALL nwp_rrtm_radiation_repartition ( p_sim_time,pt_patch, &
            & ext_data,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,                &
            & pt_diag, prm_diag, lnd_prog   )
!         CASE(2)
!           CALL nwp_omp_rrtm_interface ( p_sim_time,pt_patch, &
!             & ext_data, lnd_diag, pt_diag, prm_diag, lnd_prog )
          
        CASE default
          CALL nwp_rrtm_radiation ( p_sim_time,pt_patch, &
            & ext_data,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,    &
            & pt_diag, prm_diag, lnd_prog   )

       END SELECT
       
      ELSE 

        CALL nwp_rrtm_radiation_reduced ( p_sim_time,pt_patch,pt_par_patch,  &
          & pt_par_int_state, pt_par_grf_state,ext_data,                     &
          & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,                                   &
          & pt_diag,prm_diag, lnd_prog )
          
      ENDIF

      RETURN
    ENDIF !inwp_radiation = 1


    IF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 .AND. .NOT. lredgrid) THEN
    
      CALL nwp_rg_radiation ( p_sim_time, datetime, pt_patch, &
        & ext_data,pt_prog,pt_diag,prm_diag, lnd_prog )


    ELSEIF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 .AND. lredgrid) THEN

      CALL nwp_rg_radiation_reduced ( p_sim_time, datetime, pt_patch,pt_par_patch, &
        & pt_par_int_state, pt_par_grf_state,ext_data,pt_prog, pt_diag, prm_diag,  &
        & lnd_prog )


    ENDIF !inwp_radiation = 2

  END SUBROUTINE nwp_radiation
  !---------------------------------------------------------------------------------------



  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rg_radiation ( p_sim_time, datetime, pt_patch, &
    & ext_data,pt_prog,pt_diag,prm_diag,lnd_prog )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

!    REAL(wp), PARAMETER::  &
!      & zqco2 = 0.5014E-03_wp*353.9_wp/330._wp ! CO2 (mixing ratio 353.9 ppm (like vmr_co2))

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time
    
    TYPE(t_datetime),            INTENT(in) :: datetime
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_external_data)  , INTENT(inout)  :: ext_data
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog

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

    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
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

    ! CO2 help variable for Ritter-Geleyn scheme
    zqco2 = 0.5014E-03_wp*vmr_co2/330.e-6_wp

    CALL nwp_rg_ozon_aerosol ( p_sim_time, datetime, pt_patch, ext_data, &
      & pt_diag,prm_diag,zduo3,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )
    
    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            &
      & cosmu0_dark  = cosmu0_dark,                       &
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad     ,&
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, &
      & p_sim_time   = p_sim_time,                        &
      & pt_patch     = pt_patch,                          &
     !& zsmu0        = prm_diag%cosmu0(1,1),              &
      & zsmu0        = prm_diag%cosmu0(:,:),              &
      & zsct         = zsct )

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF (msg_level >= 15) &
      &           CALL message('mo_nwp_rad_interface', 'RG radiation on full grid')

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
#ifdef __xlC__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,zi0,losol,lo_sol)
#else
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,zi0,losol,lo_sol),SCHEDULE(guided)
#endif
    
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

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)

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

      CALL fesft ( &
                                !  Input:
        & pti = pt_diag%temp_ifc (:,:,jb) , &! Temperature at layer boundaries
        & pdp = pt_diag%dpres_mc (:,:,jb), &! pressure thickness
        & pclc_in= prm_diag%tot_cld  (:,:,jb,icc) , &
        & pqv = prm_diag%tot_cld(:,:,jb,iqv), &
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
        & palso = prm_diag%albvisdif(:,jb), & ! solar surface albedo
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
  SUBROUTINE nwp_rg_radiation_reduced ( p_sim_time, datetime, pt_patch,pt_par_patch, &
    & pt_par_int_state, pt_par_grf_state,ext_data,pt_prog,pt_diag,prm_diag, &
    & lnd_prog )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

!    REAL(wp), PARAMETER::  &
!      & zqco2 = 0.5014E-03_wp*353.9_wp/330._wp ! CO2 (mixing ratio 353.9 ppm (like vmr_co2))

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_datetime),            INTENT(in) :: datetime 
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_int_state),    TARGET,INTENT(in) :: pt_par_int_state  !< " for parent grid
    TYPE(t_gridref_state),TARGET,INTENT(in) :: pt_par_grf_state  !< grid refinement state
    TYPE(t_external_data)       ,INTENT(inout):: ext_data
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog

    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_cosmu0   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albvisdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albeff(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albefffac(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tsfc     (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_o3       (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tot_cld  (:,:,:,:)
    ! Output fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fls (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_flsp (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_flsd (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_flsu (:,:)
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

    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    INTEGER:: jc,jk,jb
    INTEGER:: jg                                !domain id
    INTEGER:: nlev, nlevp1, nlev_rg, nlevp1_rg  !< number of full and half levels
    INTEGER:: nblks_par_c                       !nblks for reduced grid

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

    ! CO2 help variable for Ritter-Geleyn scheme
    zqco2 = 0.5014E-03_wp*vmr_co2/330.e-6_wp

    CALL nwp_rg_ozon_aerosol ( p_sim_time, datetime, pt_patch, ext_data, &
      & pt_diag,prm_diag,zduo3,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )

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

    IF (msg_level >= 15) &
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

    ! Add extra layer for atmosphere above model top if requested
    IF (atm_phy_nwp_config(jg)%latm_above_top) THEN
      nlev_rg   = nlev + 1
      nlevp1_rg = nlevp1 + 1
    ELSE
      nlev_rg   = nlev
      nlevp1_rg = nlevp1
    ENDIF

    ALLOCATE (                                       &
      zrg_cosmu0   (nproma,          nblks_par_c),   &
      zrg_albvisdif(nproma,          nblks_par_c),   &
      zrg_albeff   (nproma,          nblks_par_c),   &
      zrg_albefffac(nproma,          nblks_par_c),   &
      zrg_alb_ther (nproma,          nblks_par_c),   &
      zrg_pres_sfc (nproma,          nblks_par_c),   &
      zrg_tsfc     (nproma,          nblks_par_c),   &
      zrg_temp_ifc (nproma,nlevp1_rg,nblks_par_c),   &
      zrg_dpres_mc (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_sqv      (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_duco2    (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_o3       (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq1     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq2     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq3     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq4     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq5     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_tot_cld  (nproma,nlev_rg  ,nblks_par_c,4), &
      zrg_fls      (nproma,nlevp1_rg,nblks_par_c),   &
      zrg_flsp     (nproma,          nblks_par_c),   &
      zrg_flsd     (nproma,          nblks_par_c),   &
      zrg_flsu     (nproma,          nblks_par_c),   &
      zrg_lwflxall (nproma,nlevp1_rg,nblks_par_c),   &
      zrg_trsolall (nproma,nlevp1_rg,nblks_par_c)  )


    rl_start = 1  !DR 3 seems to be sufficient 
                  ! SR radiation is not set up to handle boundaries of nested domains
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    
    ! *** this parallel section will be removed later once real data are
    !     are available as input for radiation ***
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
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

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)

      ! CO2 (mixing ratio 353.9 ppm as vmr_co2)
      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zduco2(jc,jk,jb) = zqco2 * pt_diag%dpres_mc (jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

    CALL upscale_rad_input_rg( pt_patch%id, pt_par_patch%id,  nlev_rg, nlevp1_rg,            &
      & prm_diag%cosmu0, prm_diag%albvisdif, alb_ther, pt_diag%temp_ifc,                     &
      & pt_diag%dpres_mc, prm_diag%tot_cld, zsqv ,zduco2, zduo3,                             &
      & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,pt_diag%pres_sfc,pt_diag%pres_ifc,                     &
      & zrg_cosmu0, zrg_albvisdif, zrg_alb_ther, zrg_temp_ifc, zrg_dpres_mc,                 &
      & zrg_tot_cld, zrg_sqv ,zrg_duco2, zrg_o3,                                             &
      & zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5,zrg_pres_sfc     )

    rl_start = grf_ovlparea_start_c
    rl_end   = min_rlcell_int

    i_startblk = ptr_pp%cells%start_blk(rl_start,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(rl_end,i_chidx)


!$OMP PARALLEL
#ifdef __xlC__
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,lo_sol,losol,zi0)
#else
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,lo_sol,losol,zi0) ,SCHEDULE(guided)
#endif
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

      CALL fesft ( &
                                         !  Input:
        & pti = zrg_temp_ifc (:,:,jb) , &! Temperature at layer boundaries
        & pdp = zrg_dpres_mc (:,:,jb), &! pressure thickness
        & pclc_in= zrg_tot_cld  (:,:,jb,icc) , &
        & pqv = zrg_tot_cld(:,:,jb,iqv), &
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
        & ki3e = nlev_rg,&
        & ki1sc= i_startidx, &
        & ki1ec= i_endidx, &
        & lsolar = losol, &! control switch for solar calculations
        !          & lsolar = .TRUE., &! control switch for solar calculations
        & lthermal =.TRUE., &
        & lcrf = .FALSE., &! control switch for cloud-free calcul.
                                ! Output:
        & pflt  = zrg_lwflxall(:,:,jb) ,& ! Thermal radiative fluxes at each layer boundary
        & pfls  = zrg_fls  (:,:,jb),  &! solar radiative fluxes at each layer boundary
        & pflsp = zrg_flsp (:,jb), &
        & pflsd = zrg_flsd (:,jb), &
        & pflsu = zrg_flsu (:,jb) &
        & )

      zi0 (i_startidx:i_endidx) = zrg_cosmu0(i_startidx:i_endidx,jb) * zsct
      ! compute sw transmissivity trsolall from sw fluxes
      DO jk = 1,nlevp1_rg
        DO jc = i_startidx,i_endidx
          ! This is needed to avoid false synchronization errors
          IF (zrg_cosmu0(jc,jb) < 1.e-8_wp) THEN
            zrg_trsolall(jc,jk,jb) = 0._wp
          ELSE
            zrg_trsolall(jc,jk,jb) = zrg_fls(jc,jk,jb) / zi0(jc)
          ENDIF
        ENDDO
      ENDDO

      DO jc = i_startidx,i_endidx
        
        IF (zrg_cosmu0(jc,jb) < 1.e-8_wp) THEN
          zrg_albefffac(jc,jb) = 1._wp
          zrg_albeff(jc,jb)    = 0.5_wp
        ELSE
          zrg_albeff(jc,jb) = zrg_flsu(jc,jb) / ( zrg_flsp(jc,jb) + zrg_flsd(jc,jb) )
          zrg_albefffac(jc,jb) = zrg_albeff(jc,jb) / zrg_albvisdif(jc,jb)
        ENDIF
        
        zrg_tsfc(jc,jb) = zrg_temp_ifc(jc,nlevp1_rg,jb)
        
      ENDDO

      ! to avoid division by zero inside downscale_rad_output_rg
      zrg_flsp(i_startidx:i_endidx,jb)=MAX(zrg_flsp(i_startidx:i_endidx,jb),1.e-9_wp)
      zrg_flsd(i_startidx:i_endidx,jb)=MAX(zrg_flsd(i_startidx:i_endidx,jb),1.e-9_wp)

    ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL


    CALL downscale_rad_output_rg(          &
      & jg           = pt_patch%id,        &
      & jgp          = pt_par_patch%id,    &
      & nlev_rg      = nlev_rg,            &
      & rg_lwflxall  = zrg_lwflxall,       &
      & rg_trsolall  = zrg_trsolall,       &
      & tsfc_rg      = zrg_tsfc,           &
      & albeff_rg    = zrg_albeff,         &
      & albefffac_rg = zrg_albefffac,      &
      & flsp_rg      = zrg_flsp,           &
      & flsd_rg      = zrg_flsd,           & 
      & alb_ther_rg  = zrg_alb_ther,       &
      & cosmu0_rg    = zrg_cosmu0,         &
      & tot_cld_rg   = zrg_tot_cld,        &
      & dpres_mc_rg  = zrg_dpres_mc,       &
      & pres_sfc_rg  = zrg_pres_sfc,       &
      & tsfc         = prm_diag%tsfctrad,  &
      & albvisdif    = prm_diag%albvisdif, &
      & zsct         = zsct,               &
      & lwflxall     = prm_diag%lwflxall,  &
      & trsolall     = prm_diag%trsolall )

    DEALLOCATE (zrg_cosmu0,zrg_tsfc,zrg_albvisdif,zrg_alb_ther,             &
      & zrg_pres_sfc,zrg_temp_ifc,zrg_dpres_mc,zrg_sqv,zrg_duco2,zrg_o3,    &
      & zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5,                       &
      & zrg_tot_cld,zrg_fls,zrg_flsp,zrg_flsd,zrg_flsu,zrg_lwflxall,zrg_trsolall)

  END SUBROUTINE nwp_rg_radiation_reduced

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-09-22)
  !!
  SUBROUTINE nwp_rg_ozon_aerosol ( p_sim_time, datetime, pt_patch, ext_data, &
    & pt_diag,prm_diag,zduo3,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_datetime),            INTENT(in) :: datetime
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info
    TYPE(t_external_data),       INTENT(inout) :: ext_data
    TYPE(t_nh_diag), TARGET, INTENT(in)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag

    REAL(wp), INTENT(out) ::                          &
      & zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    REAL(wp)::                                                &
    ! for ozone:    
      & zptop32(nproma,pt_patch%nblks_c),                     &
      & zo3_hm (nproma,pt_patch%nblks_c),                     &
      & zo3_top (nproma,pt_patch%nblks_c),                    &
      & zpbot32(nproma,pt_patch%nblks_c),                     &
      & zo3_bot (nproma,pt_patch%nblks_c),                    &
    ! for aerosols:
      & zsign(nproma,pt_patch%nlevp1),                        &
      & zvdaes(nproma,pt_patch%nlevp1),                       &
      & zvdael(nproma,pt_patch%nlevp1),                       &
      & zvdaeu(nproma,pt_patch%nlevp1),                       &
      & zvdaed(nproma,pt_patch%nlevp1),                       &
      & zaetr_top(nproma,pt_patch%nblks_c), zaetr_bot, zaetr, &
      & zaeqdo   (nproma,pt_patch%nblks_c), zaeqdn,           &
      & zaequo   (nproma,pt_patch%nblks_c), zaequn,           &
      & zaeqlo   (nproma,pt_patch%nblks_c), zaeqln,           &
      & zaeqso   (nproma,pt_patch%nblks_c), zaeqsn,           &
      ! for Tegen aerosol
      & z_aer_ss(nproma,pt_patch%nblks_c), &
      & z_aer_or(nproma,pt_patch%nblks_c), &
      & z_aer_bc(nproma,pt_patch%nblks_c), &
      & z_aer_su(nproma,pt_patch%nblks_c), &
      & z_aer_du(nproma,pt_patch%nblks_c), zw
    ! for aerosols:   
    REAL(wp), PARAMETER::                               &
      & zaeops = 0.05_wp,                               &
      & zaeopl = 0.2_wp,                                &
      & zaeopu = 0.1_wp,                                &
      & zaeopd = 1.9_wp,                                &
      & ztrpt  = 30.0_wp,                               &
      & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
      & zvobga = 0.007_wp /  19330.0_wp ,               &
      & zstbga = 0.045_wp /  19330.0_wp!, &
!      & zaeadk(1:3) = (/0.3876E-03_wp,0.6693E-02_wp,0.8563E-03_wp/)
    
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

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

    ! O3
    SELECT CASE (irad_o3)
    !CASE(io3_ape)
      ! APE ozone: do nothing because 
      ! everything is already set in nwp_phy_init
    CASE (6)
      CALL calc_o3_clim(                             &
        & kbdim      = nproma,                       &
        & jg         = jg,                           &
        & p_inc_rad  = atm_phy_nwp_config(jg)%dt_rad,&
        & z_sim_time = p_sim_time,                   &
        & pt_patch   = pt_patch,                     &
        & zvio3      = prm_diag%vio3,                &
        & zhmo3      = prm_diag%hmo3  )
    CASE (7)
      CALL calc_o3_gems(pt_patch,datetime,pt_diag,ext_data)
    END SELECT

    IF ( irad_aero == 6 ) CALL month2hour (datetime, imo1, imo2, zw )
    
      
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
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            ext_data%atm%o3(jc,jk,jb) = zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb)
          ENDDO
        ENDDO
      ELSEIF ( irad_o3 == 6 .AND. irad_aero == 6 ) THEN
        
        DO jc = i_startidx,i_endidx

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
                    
!          z_aer_ss(jc,jb)= ext_data%atm_td%aer_ss(jc,jb,datetime%month)
!          z_aer_or(jc,jb)= ext_data%atm_td%aer_org(jc,jb,datetime%month)
!          z_aer_bc(jc,jb)= ext_data%atm_td%aer_bc(jc,jb,datetime%month)
!          z_aer_su(jc,jb)= ext_data%atm_td%aer_so4(jc,jb,datetime%month)
!          z_aer_du(jc,jb)= ext_data%atm_td%aer_dust(jc,jb,datetime%month)
          
        ENDDO

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
          zaeqso   (jc,jb) = z_aer_ss(jc,jb)*zvdaes(jc,1)
          zaeqlo   (jc,jb) = ( z_aer_or(jc,jb)+z_aer_su(jc,jb) )*zvdael(jc,1)
          zaequo   (jc,jb) = z_aer_bc(jc,jb)              *zvdaeu(jc,1)
          zaeqdo   (jc,jb) =  z_aer_du(jc,jb)              *zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
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

            ext_data%atm%o3(jc,jk,jb) = zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb)
            
          ENDDO
        ENDDO
                
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
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            ext_data%atm%o3(jc,jk,jb) = zduo3(jc,jk,jb)/pt_diag%dpres_mc(jc,jk,jb)
          ENDDO
        ENDDO

        zaeq1(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq2(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq3(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq4(i_startidx:i_endidx,:,jb) = 0.0_wp
        zaeq5(i_startidx:i_endidx,:,jb) = 0.0_wp
        
      ELSEIF ( irad_aero == 5 ) THEN !aerosols, but other or no ozone:

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

        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zduo3(jc,jk,jb) = ext_data%atm%o3(jc,jk,jb) * pt_diag%dpres_mc(jc,jk,jb)
          ENDDO
        ENDDO
      ELSEIF (irad_aero == 6 ) THEN !aerosols, but other ozone:

        DO jc = i_startidx,i_endidx

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
          
!          z_aer_ss(jc,jb)= ext_data%atm_td%aer_ss(jc,jb,datetime%month)
!          z_aer_or(jc,jb)= ext_data%atm_td%aer_org(jc,jb,datetime%month)
!          z_aer_bc(jc,jb)= ext_data%atm_td%aer_bc(jc,jb,datetime%month)
!          z_aer_su(jc,jb)= ext_data%atm_td%aer_so4(jc,jb,datetime%month)
!          z_aer_du(jc,jb)= ext_data%atm_td%aer_dust(jc,jb,datetime%month)

        ENDDO
        
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
          zaeqso   (jc,jb) = z_aer_ss(jc,jb)              *zvdaes(jc,1)
          zaeqlo   (jc,jb) = ( z_aer_or(jc,jb)+z_aer_su(jc,jb) )*zvdael(jc,1)
          zaequo   (jc,jb) =  z_aer_bc(jc,jb)              *zvdaeu(jc,1)
          zaeqdo   (jc,jb) =  z_aer_du(jc,jb)              *zvdaed(jc,1)
          zaetr_top(jc,jb) = 1.0_wp
        ENDDO

        ! loop over layers
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
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

        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zduo3(jc,jk,jb) = ext_data%atm%o3(jc,jk,jb) * pt_diag%dpres_mc(jc,jk,jb)
          ENDDO
        ENDDO
        
        
      ELSEIF (irad_o3 /= 0) THEN !no aerosols and other ozone
        
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            zduo3(jc,jk,jb) = ext_data%atm%o3(jc,jk,jb) * pt_diag%dpres_mc(jc,jk,jb)
            zaeq1(jc,jk,jb) = 0.0_wp
            zaeq2(jc,jk,jb) = 0.0_wp
            zaeq3(jc,jk,jb) = 0.0_wp
            zaeq4(jc,jk,jb) = 0.0_wp
            zaeq5(jc,jk,jb) = 0.0_wp
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

  END SUBROUTINE nwp_rg_ozon_aerosol
  
END MODULE mo_nwp_rad_interface

