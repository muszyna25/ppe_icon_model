!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (ecRad, RRTM or Ritter-Geleyn).
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_nwp_rad_interface

  USE mo_exception,            ONLY: finish, message
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_nh_testcases_nml,     ONLY: nh_test_name, albedo_set
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma
  USE mo_impl_constants,       ONLY: MODIS, min_rlcell_int
  USE mo_kind,                 ONLY: wp
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation_config,     ONLY: albedo_type, irad_co2, irad_n2o, irad_ch4, irad_cfc11, irad_cfc12
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation,             &
    &                                nwp_rrtm_radiation_reduced,     &
    &                                nwp_ozon_aerosol
  USE mo_nwp_rg_interface,     ONLY: nwp_rg_radiation,               &
    &                                nwp_rg_radiation_reduced
#ifdef __ECRAD
  USE mo_nwp_ecrad_interface,  ONLY: nwp_ecrad_radiation,            &
    &                                nwp_ecrad_radiation_reduced
  USE mo_ecrad,                ONLY: ecrad_conf
#endif
  USE mo_albedo,               ONLY: sfc_albedo, sfc_albedo_modis
  USE mtime,                   ONLY: datetime
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_gpu_util,         ONLY: gpu_d2h_nh_nwp, gpu_h2d_nh_nwp
  USE mo_bc_greenhouse_gases,  ONLY: bc_greenhouse_gases_time_interpolation
#if defined( _OPENACC )
  USE mo_mpi,                  ONLY: i_am_accel_node, my_process_is_work
#endif

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_radiation
  

 CONTAINS
  
  !---------------------------------------------------------------------------------------
  !>
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1), Ritter-Geleyn (2), or 
  !! ecRad(4).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
    & ext_data, lnd_diag, pt_prog, pt_diag, prm_diag, lnd_prog, wtr_prog, linit)

    CHARACTER(len=*), PARAMETER :: &
      &  routine = 'mo_nwp_rad_interface:nwp_radiation'

    LOGICAL,                 INTENT(in)    :: lredgrid        !< use reduced grid for radiation
    LOGICAL, OPTIONAL,       INTENT(in)    :: linit

    REAL(wp),                INTENT(in)    :: p_sim_time

    TYPE(datetime), POINTER, INTENT(in)    :: mtime_datetime
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch     !<grid/patch info.
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_external_data),   INTENT(inout) :: ext_data
    TYPE(t_lnd_diag),        INTENT(in)    :: lnd_diag   !<diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog    !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag    !<the diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog   ! time level new
    TYPE(t_wtr_prog),        INTENT(in)    :: wtr_prog   ! time level new

    REAL(wp) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    INTEGER :: jg, irad
    INTEGER :: jb, jc          !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    LOGICAL :: lacc

    REAL(wp):: zsct        ! solar constant (at time of year)
    REAL(wp):: cosmu0_dark ! minimum cosmu0, for smaller values no shortwave calculations



    ! patch ID
    jg = pt_patch%id
    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

    ! openACC flag (initialization run is left on)
    IF(PRESENT(linit)) THEN
      lacc = .NOT. linit
    ELSE
      lacc = .FALSE.
    ENDIF


    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    IF(ANY((/irad_co2,irad_cfc11,irad_cfc12,irad_n2o,irad_ch4/) == 4)) THEN 
      ! Interpolate greenhouse gas concentrations to the current date and time, 
      !   placing the annual means at the mid points of the current and preceding or following year,
      !   if the current date is in the 1st or 2nd half of the year, respectively.
      ! The data file containing the greenhouse gas concentration is read in the initialisation 
      !   of the NWP physics
      CALL bc_greenhouse_gases_time_interpolation(mtime_datetime)
    END IF

    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation )
    CASE (1)
      ! RRTM
      ! In radiative transfer routine RRTM skips all points with cosmu0<=0. That's why 
      ! points to be skipped need to be marked with a value <=0
      cosmu0_dark = -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    CASE (2)
      ! Ritter-Geleyn
      ! Skipping of points is performed on block- rather than cell-level. I.e. if a block 
      ! contains at least 1 point with cosmu0>1.E-8, radiatve transfer is computed for 
      ! the entire block. Therefore cosmu0_dark = -1.e-9_wp does not work here (crashes).
      ! For all points cosmu0 must be <0.
      cosmu0_dark =  1.e-9_wp   ! minimum cosmu0, for smaller values no shortwave calculations
    CASE (4)
      ! ecRad
      ! ecRad skips cosmu0<=0 as well.
      cosmu0_dark = -1.e-9_wp
    END SELECT


    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            & !in
      & cosmu0_dark  = cosmu0_dark,                       & !in
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     & !in
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, & !in
      & p_sim_time   = p_sim_time,                        & !in
      & pt_patch     = pt_patch,                          & !in
      & zsmu0        = prm_diag%cosmu0(:,:),              & !out
      & zsct         = zsct,                              & !out, optional
      & lacc         = lacc                               ) !in




    ! Compute tile-based and aggregated surface-albedo
    !
    IF ( albedo_type == MODIS ) THEN
      ! MODIS albedo
      CALL sfc_albedo_modis(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag, lacc)
!!
!! in case sfc_albedo_scm has been implemented the following lines should be
!! active instead of the RCEMIP_analytical IF-BLOCK, see below.
!!    ELSE IF ( albedo_type == 3 ) THEN
!!       albedo_value = 0.07_wp                  !should equal albedo_set for RCEMIP_analytical
!!       CALL sfc_albedo_scm(pt_patch, ext_data, albedo_value, prm_diag)
!!
    ELSE
#ifdef _OPENACC
      IF (lacc) CALL finish('nwp_radiation','sfc_albedo not ported to gpu')
#endif
      ! albedo based on tabulated bare soil values
      CALL sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)
    ENDIF

    !FOR RCEMIP, SET ALBEDO TO FIXED VALUE
!! 
!! Use the upper 'IF ( albedo_type == 3 )' as soon as sfc_albedo_scm has been introduced !!!!
!!
    IF( nh_test_name=='RCEMIP_analytical' ) THEN
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          prm_diag%albdif(jc,jb)    = albedo_set
          prm_diag%albvisdif(jc,jb) = albedo_set
          prm_diag%albnirdif(jc,jb) = albedo_set
          prm_diag%albvisdir(jc,jb) = albedo_set
          prm_diag%albnirdir(jc,jb) = albedo_set
        ENDDO
      ENDDO  ! jb
    ENDIF

    
    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------
    !

#ifdef _OPENACC
    IF(lacc) THEN
      CALL message('mo_nh_interface_nwp', 'Device to host copy before Radiation. This needs to be removed once port is finished!')
      CALL gpu_d2h_nh_nwp(pt_patch, prm_diag, ext_data)
      i_am_accel_node = .FALSE.
    ENDIF
#endif
    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation)
    CASE (1) ! RRTM

      CALL nwp_ozon_aerosol ( p_sim_time, mtime_datetime, pt_patch, ext_data, &
        & pt_diag, prm_diag, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 )
    
      IF ( .NOT. lredgrid ) THEN
          
        CALL nwp_rrtm_radiation ( mtime_datetime, pt_patch, ext_data, &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,        &
          & pt_diag, prm_diag, lnd_prog )
       
      ELSE 

        CALL nwp_rrtm_radiation_reduced ( mtime_datetime, pt_patch,pt_par_patch, ext_data, &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                             &
          & pt_diag, prm_diag, lnd_prog )
          
      ENDIF

    CASE (2) ! Ritter-Geleyn

      IF (.NOT. lredgrid) THEN
        CALL nwp_rg_radiation ( p_sim_time, mtime_datetime, pt_patch, &
          & ext_data,pt_prog,pt_diag,prm_diag, lnd_prog, zsct )
      ELSE
        CALL nwp_rg_radiation_reduced ( p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
          & ext_data, pt_prog, pt_diag, prm_diag, lnd_prog, zsct )
      ENDIF

    CASE (4) ! ecRad
#ifdef __ECRAD
      CALL nwp_ozon_aerosol ( p_sim_time, mtime_datetime, pt_patch, ext_data, &
        & pt_diag, prm_diag, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 )

      IF (.NOT. lredgrid) THEN
        CALL nwp_ecRad_radiation ( mtime_datetime, pt_patch, ext_data,      &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                              &
          & pt_diag, prm_diag, lnd_prog, ecrad_conf )
      ELSE
        CALL nwp_ecRad_radiation_reduced ( mtime_datetime, pt_patch,pt_par_patch, &
          & ext_data, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                          &
          & pt_diag, prm_diag, lnd_prog, ecrad_conf )
      ENDIF
#else
      CALL finish(routine,  &
        &      'atm_phy_nwp_config(jg)%inwp_radiation = 4 needs -D__ECRAD.')
#endif
    END SELECT ! inwp_radiation
#ifdef _OPENACC
    IF(lacc) THEN
      CALL message('mo_nh_interface_nwp', 'Host to device copy after Radiation. This needs to be removed once port is finished!')
      CALL gpu_h2d_nh_nwp(pt_patch, prm_diag, ext_data)
      i_am_accel_node = my_process_is_work()
    ENDIF
#endif

  END SUBROUTINE nwp_radiation


END MODULE mo_nwp_rad_interface

