!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (ecRad and RRTM).
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_rad_interface

  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma
  USE mo_impl_constants,       ONLY: MODIS, min_rlcell_int, SUCCESS
  USE mo_kind,                 ONLY: wp
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation_config,     ONLY: albedo_type, albedo_fixed,      &
    &                                irad_co2, irad_n2o, irad_ch4, irad_cfc11, irad_cfc12
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation,             &
    &                                nwp_rrtm_radiation_reduced,     &
    &                                nwp_ozon_aerosol
#ifdef __ECRAD
  USE mo_nwp_ecrad_interface,  ONLY: nwp_ecrad_radiation,            &
    &                                nwp_ecrad_radiation_reduced
  USE mo_ecrad,                ONLY: ecrad_conf
#endif
  USE mo_albedo,               ONLY: sfc_albedo, sfc_albedo_modis, sfc_albedo_scm
  USE mtime,                   ONLY: datetime
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_gpu_util,         ONLY: gpu_d2h_nh_nwp, gpu_h2d_nh_nwp
  USE mo_bc_greenhouse_gases,  ONLY: bc_greenhouse_gases_time_interpolation
#if defined( _OPENACC )
  USE mo_mpi,                  ONLY: i_am_accel_node, my_process_is_work
#endif
  USE mo_bc_aeropt_kinne,      ONLY: set_bc_aeropt_kinne
  USE mo_bc_aeropt_cmip6_volc, ONLY: add_bc_aeropt_cmip6_volc
  USE mo_radiation_config,     ONLY: irad_aero
  USE mo_loopindices,          ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_radiation
  

 CONTAINS
  
  !---------------------------------------------------------------------------------------
  !>
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1) or ecRad(4).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
    & ext_data, lnd_diag, pt_prog, pt_diag, prm_diag, lnd_prog, wtr_prog, zf, dz, linit)

    CHARACTER(len=*), PARAMETER :: &
      &  routine = 'mo_nwp_rad_interface:nwp_radiation'

    LOGICAL,                 INTENT(in)    :: lredgrid        !< use reduced grid for radiation
    LOGICAL, OPTIONAL,       INTENT(in)    :: linit

    REAL(wp),                INTENT(in)    :: p_sim_time   !< simulation time
    REAL(wp),                INTENT(in)    :: zf(:,:,:)    !< model full layer height
    REAL(wp),                INTENT(in)    :: dz(:,:,:)    !< Layer thickness

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

    REAL(wp), ALLOCATABLE :: &
      & od_lw_vr(:,:,:) , & !< LW optical thickness of aerosols    (vertically reversed)
      & od_sw_vr(:,:,:) , & !< SW aerosol optical thickness        (vertically reversed)
      & g_sw_vr (:,:,:) , & !< SW aerosol asymmetry factor         (vertically reversed)
      & ssa_sw_vr(:,:,:), & !< SW aerosol single scattering albedo (vertically reversed)
      & od_lw(:,:,:,:)  , & !< LW optical thickness of aerosols
      & od_sw(:,:,:,:)  , & !< SW aerosol optical thickness
      & g_sw (:,:,:,:)  , & !< SW aerosol asymmetry factor
      & ssa_sw(:,:,:,:)     !< SW aerosol single scattering albedo

    INTEGER :: jg, irad
    INTEGER :: jb, jc, jk              !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: istat
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

#ifdef __ECRAD
    IF (ANY( irad_aero == (/12,13,14,15/) )) THEN

      ALLOCATE(od_lw_vr (nproma,pt_patch%nlev,ecrad_conf%n_bands_lw)                   , &
      &        od_sw_vr (nproma,pt_patch%nlev,ecrad_conf%n_bands_sw)                   , &
      &        g_sw_vr  (nproma,pt_patch%nlev,ecrad_conf%n_bands_sw)                   , &
      &        ssa_sw_vr(nproma,pt_patch%nlev,ecrad_conf%n_bands_sw)                   , &
      &        od_lw    (nproma,pt_patch%nlev,pt_patch%nblks_c,ecrad_conf%n_bands_lw)  , &
      &        od_sw    (nproma,pt_patch%nlev,pt_patch%nblks_c,ecrad_conf%n_bands_sw)  , &
      &        ssa_sw   (nproma,pt_patch%nlev,pt_patch%nblks_c,ecrad_conf%n_bands_sw)  , &
      &        g_sw     (nproma,pt_patch%nlev,pt_patch%nblks_c,ecrad_conf%n_bands_sw)  , &
      &        STAT=istat                                                                )

      IF(istat /= SUCCESS) CALL finish(routine, 'Allocation of od_lw_vr,od_sw_vr, g_sw_vr, &
                                       ssa_sw_vr, od_lw, od_sw, ssa_sw, g_sw failed'       )

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,i_startidx,i_endidx, od_lw_vr, od_sw_vr, ssa_sw_vr, g_sw_vr) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)

        IF (i_startidx>i_endidx) CYCLE

        od_lw_vr(:,:,:)  = 0.0_wp
        od_sw_vr(:,:,:)  = 0.0_wp
        ssa_sw_vr(:,:,:) = 1.0_wp
        g_sw_vr (:,:,:)  = 0.0_wp

        IF (ANY( irad_aero == (/12,13,15/) )) THEN
          CALL set_bc_aeropt_kinne(mtime_datetime, jg, 1, i_endidx, &
            & nproma, pt_patch%nlev, jb, ecrad_conf%n_bands_sw,     &
            & ecrad_conf%n_bands_lw, zf(:,:,jb), dz(:,:,jb),        &
            & od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                    &
            & g_sw_vr (:,:,:), od_lw_vr(:,:,:)                      )
        END IF
        IF (ANY( irad_aero == (/14,15/) )) THEN 
          CALL add_bc_aeropt_cmip6_volc(mtime_datetime, jg, 1,      &
            & i_endidx, nproma, pt_patch%nlev, jb,                  &
            & ecrad_conf%n_bands_sw, ecrad_conf%n_bands_lw,         &
            & zf(:,:,jb), dz(:,:,jb), od_sw_vr(:,:,:),              &
            & ssa_sw_vr(:,:,:), g_sw_vr (:,:,:), od_lw_vr(:,:,:)    )
        END IF
        !
        ! Vertically reverse the fields:
        DO jk = 1, pt_patch%nlev
          od_lw (:,jk,jb,:) = od_lw_vr (:,pt_patch%nlev-jk+1,:)
          od_sw (:,jk,jb,:) = od_sw_vr (:,pt_patch%nlev-jk+1,:)
          ssa_sw(:,jk,jb,:) = ssa_sw_vr(:,pt_patch%nlev-jk+1,:)
          g_sw  (:,jk,jb,:) = g_sw_vr  (:,pt_patch%nlev-jk+1,:)
        ENDDO

      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      DEALLOCATE(od_lw_vr, od_sw_vr, ssa_sw_vr, g_sw_vr, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_lw_vr, &
                                       od_sw_vr, ssa_sw_vr, g_sw_vr failed')

    END IF
#endif

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
    ELSE IF ( albedo_type == 3 ) THEN
      ! globally fixed albedo value for SCM and RCEMIP applications
      CALL sfc_albedo_scm(pt_patch, albedo_fixed, prm_diag)
    ELSE
#ifdef _OPENACC
      IF (lacc) CALL finish('nwp_radiation','sfc_albedo not ported to gpu')
#endif
      ! albedo based on tabulated bare soil values
      CALL sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)
    ENDIF

    
    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------
    !

    !$ACC DATA CREATE(zaeq1, zaeq2, zaeq3, zaeq4, zaeq5) IF(lacc)
    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation)
    CASE (1) ! RRTM

#ifdef _OPENACC
    IF(lacc) THEN
      CALL message('mo_nh_interface_nwp', &
        &  'Device to host copy before nwp_rrtm_radiation. This needs to be removed once port is finished!')
      CALL gpu_d2h_nh_nwp(pt_patch, prm_diag, ext_data)
      i_am_accel_node = .FALSE.
    ENDIF
#endif
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

#ifdef _OPENACC
      IF(lacc) THEN
        CALL message('mo_nh_interface_nwp', &
          &  'Host to device copy after nwp_rrtm_radiation. This needs to be removed once port is finished!')
        CALL gpu_h2d_nh_nwp(pt_patch, prm_diag, ext_data)
        i_am_accel_node = my_process_is_work()
      ENDIF
#endif

    CASE (4) ! ecRad
#ifdef __ECRAD
      !$ACC WAIT
      CALL nwp_ozon_aerosol ( p_sim_time, mtime_datetime, pt_patch, ext_data, &
        & pt_diag, prm_diag, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, use_acc=lacc )

      !$ACC WAIT
      IF (.NOT. lredgrid) THEN
#ifdef _OPENACC
        IF(lacc) THEN
          CALL message('mo_nh_interface_nwp', &
            &  'Device to host copy before nwp_ecRad_radiation (full radiation grid). &
            &  This needs to be removed once port is finished!')
          CALL gpu_d2h_nh_nwp(pt_patch, prm_diag, ext_data)
          !$ACC UPDATE HOST(zaeq1, zaeq2, zaeq3, zaeq4, zaeq5) IF(lacc)
          i_am_accel_node = .FALSE.
        ENDIF
#endif
        CALL nwp_ecRad_radiation ( mtime_datetime, pt_patch, ext_data,      &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                              &
          & od_lw, od_sw, ssa_sw, g_sw,                                     &
          & pt_diag, prm_diag, pt_prog, lnd_prog, ecrad_conf, lacc )
#ifdef _OPENACC
        IF(lacc) THEN
          CALL message('mo_nh_interface_nwp', &
            &  'Host to device copy after nwp_ecRad_radiation (full radiation grid). &
            &  This needs to be removed once port is finished!')
          CALL gpu_h2d_nh_nwp(pt_patch, prm_diag, ext_data)
          i_am_accel_node = my_process_is_work()
        ENDIF
#endif
      ELSE
        !$ACC WAIT
        CALL nwp_ecRad_radiation_reduced ( mtime_datetime, pt_patch,pt_par_patch, &
          & ext_data, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                          &
          & od_lw, od_sw, ssa_sw, g_sw,                                           &
          & pt_diag, prm_diag, pt_prog, lnd_prog, ecrad_conf, lacc )
      ENDIF
#else
      CALL finish(routine,  &
        &      'atm_phy_nwp_config(jg)%inwp_radiation = 4 needs -D__ECRAD.')
#endif

    CASE DEFAULT !Invalid inwp_radiation
      WRITE (message_text, '(a,i2,a)') 'inwp_radiation = ', atm_phy_nwp_config(jg)%inwp_radiation, &
        &                  ' not valid. Valid choices are 0: none, 1:RRTM, 4:ecRad '
      CALL finish(routine,message_text)
    END SELECT ! inwp_radiation

    !$ACC END DATA

    IF( ALLOCATED(od_lw) ) THEN
      DEALLOCATE(od_lw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_lw failed.')
    ENDIF
    IF( ALLOCATED(od_sw) ) THEN
      DEALLOCATE(od_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_sw failed.')
    ENDIF
    IF( ALLOCATED(ssa_sw) ) THEN
      DEALLOCATE(ssa_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of ssa_sw failed.')
    ENDIF    
    IF( ALLOCATED(g_sw) ) THEN
      DEALLOCATE(g_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of g_sw failed.')
    ENDIF

  END SUBROUTINE nwp_radiation


END MODULE mo_nwp_rad_interface

