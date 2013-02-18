!>
!! Provide an implementation of the ocean forcing.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modification by Stephan Lorenz, MPI-M:
!!   - renaming and adjustment to ocean domain and patch_oce (2010-06)
!!   - for parallel ocean: 3-dim ocean grid in v_base        (2011-07)
!!   - adding OMIP fluxes for sea ice                        (2011-09)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_oce_bulk
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
USE mo_kind,                ONLY: wp
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: dtime, ltimer
USE mo_sync,                ONLY: sync_c, sync_patch_array
USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling
USE mo_io_units,            ONLY: filename_max
USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, p_bcast, &
  &                               p_comm_work_test, p_comm_work
USE mo_parallel_config,     ONLY: p_test_run
!USE mo_util_string,         ONLY: t_keyword_list
USE mo_util_netcdf,         ONLY: read_netcdf_data
USE mo_datetime,            ONLY: t_datetime
USE mo_time_config,         ONLY: time_config
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_ext_data_state,      ONLY: ext_data
USE mo_grid_config,         ONLY: nroot
USE mo_ocean_nml,           ONLY: iforc_oce, iforc_type, iforc_len, itestcase_oce,         &
  &                               no_tracer, n_zlev, basin_center_lat,                     &
  &                               basin_center_lon, basin_width_deg, basin_height_deg,     &
  &                               relaxation_param, wstress_coeff, i_apply_bulk,           &
  &                               relax_2d_mon_s, temperature_relaxation, irelax_2d_S,     &
  &                               NO_FORCING, ANALYT_FORC, FORCING_FROM_FILE_FLUX,         &
  &                               FORCING_FROM_FILE_FIELD, FORCING_FROM_COUPLED_FLUX,      &
  &                               FORCING_FROM_COUPLED_FIELD, i_sea_ice
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch, t_patch_3D
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_dbg_nml,             ONLY: idbg_mxmn
USE mo_oce_state,           ONLY: t_hydro_ocean_state
USE mo_exception,           ONLY: finish, message, message_text
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_physical_constants,  ONLY: rho_ref, sal_ref, als, alv, zemiss_def, stbo, tmelt, tf, &
  &                               mu, clw, rho_ref, albedoW
USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, MIN_DOLIC
USE mo_math_utilities,      ONLY: gvec2cvec, cvec2gvec
USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
USE mo_sea_ice,             ONLY: calc_bulk_flux_ice, calc_bulk_flux_oce,     &
  &                               ice_slow, ice_fast
USE mo_coupling_config,     ONLY: is_coupled_run
USE mo_icon_cpl_restart,    ONLY: icon_cpl_write_restart
USE mo_icon_cpl_exchg,      ONLY: ICON_cpl_put, ICON_cpl_get
USE mo_icon_cpl_def_field,  ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

IMPLICIT NONE

! required for reading netcdf files
INCLUDE 'netcdf.inc'

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module    = 'oceBulk     '  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

! Public interface
PUBLIC  :: update_sfcflx

! private implementation
PRIVATE :: update_sfcflx_analytical


CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE update_sfcflx(p_patch_3D, p_os, p_as, p_ice, Qatm, p_sfc_flx, jstep, datetime)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_atmos_fluxes)                        :: Qatm
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_sfc_flx)                             :: p_sfc_flx
    INTEGER, INTENT(IN)                         :: jstep
    TYPE(t_datetime), INTENT(INOUT)             :: datetime
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_sfcflx'
    INTEGER  :: jmon, jdmon, jmon1, jmon2, ylen, yday
    INTEGER  :: iniyear, curyear, offset
    INTEGER  :: jc, jb, i, no_set, k
    INTEGER  :: i_startidx_c, i_endidx_c
    REAL(wp) :: z_tmin, z_relax, rday1, rday2, dtm1, dsec
    !REAL(wp) :: z_c(nproma,n_zlev,p_patch%nblks_c)
    REAL(wp) ::  z_c2(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) ::   Tfw(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)

    ! Local declarations for coupling:
    LOGICAL               :: write_coupler_restart
    INTEGER               :: info, ierror   !< return values form cpl_put/get calls
    INTEGER               :: nbr_hor_points ! = inner and halo points
    INTEGER               :: nbr_points     ! = nproma * nblks
    INTEGER               :: nbr_fields
    INTEGER, ALLOCATABLE  :: field_id(:)
    INTEGER               :: field_shape(3)
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    REAL(wp), PARAMETER   :: seconds_per_month = 2.592e6_wp !TODO: use real month lenght
    TYPE(t_patch), POINTER:: p_patch 
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    all_cells => p_patch%cells%all

    SELECT CASE (iforc_oce)

    CASE (NO_FORCING)                !  10

    ! CALL message(TRIM(routine), 'No  forcing applied' )
      CONTINUE

    CASE (ANALYT_FORC)               !  11

      CALL update_sfcflx_analytical(p_patch_3D, p_os, p_sfc_flx)

    CASE (FORCING_FROM_FILE_FLUX)    !  12

      !-------------------------------------------------------------------------
      ! Applying annual forcing read from file in mo_ext_data:
      !  - stepping daily in monthly data (preliminary solution)

      !  calculate day and month
      jmon  = datetime%month         ! integer current month
      jdmon = datetime%day           ! integer day in month
      yday  = datetime%yeaday        ! integer current day in year
      ylen  = datetime%yealen        ! integer days in year (365 or 366)
      dsec  = datetime%daysec        ! real seconds since begin of day
      !ytim  = datetime%yeatim        ! real time since begin of year

      !jdmon = mod(jdays+1,30)-1     ! no of days in month

      ! To Do: use fraction of month for interpolation
      !frcmon= datetime%monfrc       ! fraction of month
      !rday1 = frcmon+0.5_wp
      !rday2 = 1.0_wp-rday1
      !IF (rday1 > 1.0_wp)  THEN
      !  rday2=rday1
      !  rday1=1.0_wp-rday1
      !END IF

      !njday = int(86400._wp/dtime)  ! no of timesteps per day

      ! Read forcing file in chunks of one year length fixed
      !  - #slo# 2012-02-17: first quick solution for reading NCEP data
      !  - ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean

      ! Check if file should be read:
      !   - for iforc_type=5 only - NCEP type forcing
      !   - read annual data at Jan, 1st: seconds of year are less than a timestep
      !   - or at begin of each run (must not be first of january)
      IF (iforc_type == 5) THEN
        dtm1 = dtime - 1.0_wp

        IF ( (jmon == 1 .AND. jdmon == 1 .AND. dsec < dtm1) .OR. (jstep == 1) ) THEN

          ! use initial date to define correct set (year) of reading NCEP data
          !  - with offset=0 always the first year of NCEP data is used
          iniyear = time_config%ini_datetime%year
          !curyear = time_config%cur_datetime%year  ! not updated each timestep
          curyear = datetime%year
          offset = 0
          no_set = offset + curyear-iniyear + 1 

          idt_src=2  ! output print level (1-5, fix)
       !  IF (idbg_mxmn >= idt_src) THEN
       !    WRITE(message_text,'(a,i2,a,i2,a,e15.5))') 'Read NCEP data: month=', &
       !      &  jmon,' day=',jdmon,' seconds=',dsec
       !    CALL message(TRIM(routine), message_text) 
          WRITE(message_text,'(a,3i5)') 'Read NCEP data: init. year, current year, no. of set:', &
            &                            iniyear, curyear, no_set
          CALL message(TRIM(routine), message_text) 
       !  END IF

          CALL read_forc_data_oce(p_patch, ext_data, no_set)

        END IF

      END IF

      !
      ! use annual forcing-data:
      !
      IF (iforc_len == 1)  THEN

        jmon1=1
        jmon2=1
        rday1=0.5_wp
        rday2=0.5_wp

      !
      ! interpolate monthly forcing-data daily:
      !
      ELSE IF (iforc_len == 12)  THEN

        jmon1=jmon-1
        jmon2=jmon
        rday1=REAL(15-jdmon,wp)/30.0_wp
        rday2=REAL(15+jdmon,wp)/30.0_wp
        IF (jdmon > 15)  THEN
          jmon1=jmon
          jmon2=jmon+1
          rday1=REAL(45-jdmon,wp)/30.0_wp
          rday2=REAL(jdmon-15,wp)/30.0_wp
        END IF

        IF (jmon1 ==  0) jmon1=12
        IF (jmon1 == 13) jmon1=1
        IF (jmon2 ==  0) jmon2=12
        IF (jmon2 == 13) jmon2=1

      !
      ! apply daily forcing-data directly:
      !
      ELSE

        ! - now daily data sets are read in mo_ext_data
        ! - use rday1, rday2, jmon1 = jmon2 = yday for controling correct day in year
        ! - no interpolation applied, 
        jmon1 = yday
        jmon2 = jmon1
        rday1 = 1.0_wp
        rday2 = 0.0_wp

        ! Leap year: read Feb, 28 twice since only 365 data-sets are available
        IF (ylen == 366) then
          IF (yday>59) jmon1=yday-1
          jmon2=jmon1
        ENDIF

      END IF

      !
      ! OMIP data read in mo_ext_data into variable ext_data
      !
      IF (iforc_type >= 1)  THEN

        ! provide OMIP fluxes for wind stress forcing
        ! 1:  wind_u(:,:)   !  'stress_x': zonal wind stress       [m/s]
        ! 2:  wind_v(:,:)   !  'stress_y': meridional wind stress  [m/s]

        ! ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean
        p_sfc_flx%forc_wind_u(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1) + &
          &                          rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)
        p_sfc_flx%forc_wind_v(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2) + &
          &                          rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)

       ! Wind stress boundary condition for vertical diffusion D:
       !   D = d/dz(K_v*du/dz)  where
       ! Boundary condition at surface (upper bound of D at center of first layer)
       !   derived from wind-stress boundary condition Tau read from OMIP data (or elsewhere)
       !   K_v*du/dz(surf) = F_D = Tau/Rho [ m2/s2 ]
       ! discretized:
       !   top_bc_u_c = forc_wind_u / rho_ref
       !
       ! This is equivalent to an additonal forcing term F_u in the velocity equation, i.e. outside
       ! the vertical diffusion, following MITGCM:
       !   F_u = F_D/dz = Tau / (Rho*dz)  [ m/s2 ]

       ! The devision by rho_ref is done in top_bound_cond_horz_veloc (z_scale)

      END IF

      IF (iforc_type == 2 .OR. iforc_type == 5) THEN

      !-------------------------------------------------------------------------
        ! provide OMIP fluxes for sea ice (interface to ocean)
        ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
        ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
        ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
        ! 7:  fclou(:,:),  &  ! Fractional cloud cover
        ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
        ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

        p_as%tafo(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4) + &
          &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
        !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
        p_as%tafo(:,:)  = p_as%tafo(:,:) - tmelt
        p_as%ftdew(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,5) + &
          &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,5)
        p_as%fu10(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,6) + &
          &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,6)
        p_as%fclou(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,7) + &
          &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,7)
        p_as%pao(:,:)   = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,8) + &
          &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,8)
        !  don't - change units to mb/hPa
        !p_as%pao(:,:)   = p_as%pao(:,:) !* 0.01
        p_as%fswr(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,9) + &
          &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,9)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4)
        CALL dbg_print('UpdSfc: Ext data4-ta/mon1' ,z_c2                     ,str_module,idt_src)
        z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
        CALL dbg_print('UpdSfc: Ext data4-ta/mon2' ,z_c2                     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: p_as%tafo'         ,p_as%tafo                ,str_module,idt_src)
        !---------------------------------------------------------------------

      END IF

      IF (iforc_type == 3) THEN

        !-------------------------------------------------------------------------
        ! Apply surface heat and freshwater fluxes (records 4 and 5)
        ! 4:  hflx(:,:)   !  net surface heat flux               [W/m2]
        ! 5:  fwfx(:,:)   !  net freshwater flux                 [m/s]

        p_sfc_flx%forc_hflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4) + &
          &                        rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
        p_sfc_flx%forc_fwfx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,5) + &
          &                        rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,5)

        ! #slo# This is a first try for "simple flux coupling"
        IF (i_sea_ice >= 1) THEN
          Qatm%SWnet  (:,i,:)   = 0.0_wp  ! not available - very hot shot
          DO i = 1, p_ice%kice
            Qatm%LWnet  (:,i,:)   = p_sfc_flx%forc_hflx(:,:)
          ENDDO
          Qatm%sens   (:,:,:) = 0.0_wp
          Qatm%lat    (:,:,:) = 0.0_wp
          Qatm%dsensdT(:,:,:) = 0.0_wp
          Qatm%dlatdT (:,:,:) = 0.0_wp
          Qatm%dLWdT  (:,:,:) = -4.0_wp * zemiss_def*StBo * (p_ice%Tsurf(:,:,:) + tmelt)**3

          ! Fluxes into the water are the same as into the ice
          Qatm%SWnetw (:,:)   = 0.0_wp  ! not available - very hot shot
          Qatm%LWnetw (:,:)   = p_sfc_flx%forc_hflx(:,:)
          Qatm%sensw  (:,:)   = 0.0_wp
          Qatm%latw   (:,:)   = 0.0_wp

          ! p_sfc_flx%forc_hflx is recalculated in upper_ocean_TS in mo_sea_ice.f90, called by
          ! ice_slow

        ENDIF

      END IF
      
      ! this is used for "intermediate complexity flux forcing"
      IF (iforc_type == 4) THEN

        !-------------------------------------------------------------------------
        ! Apply 4 parts of surface heat and 2 parts of freshwater fluxes (records 4 to 9)
        ! 4:  swflx(:,:)   !  surface short wave heat flux        [W/m2]
        ! 5:  lwflx(:,:)   !  surface long  wave heat flux        [W/m2]
        ! 6:  ssflx(:,:)   !  surface sensible   heat flux        [W/m2]
        ! 7:  slflx(:,:)   !  surface latent     heat flux        [W/m2]
        ! 8:  prflx(:,:)   !  total precipitation flux            [m/s]
        ! 9:  evflx(:,:)   !  evaporation flux                    [m/s]

        p_sfc_flx%forc_swflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4) + &
          &                         rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
        p_sfc_flx%forc_lwflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,5) + &
          &                         rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,5)
        p_sfc_flx%forc_ssflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,6) + &
          &                         rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,6)
        p_sfc_flx%forc_slflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,7) + &
          &                         rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,7)
        p_sfc_flx%forc_prflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,8) + &
          &                         rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,8)
        p_sfc_flx%forc_evflx(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,9) + &
          &                         rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,9)

        ! sum of fluxes for ocean boundary condition
        p_sfc_flx%forc_hflx(:,:) = p_sfc_flx%forc_swflx(:,:) + p_sfc_flx%forc_lwflx(:,:) &
          &                      + p_sfc_flx%forc_ssflx(:,:) + p_sfc_flx%forc_slflx(:,:)
        p_sfc_flx%forc_fwfx(:,:) = p_sfc_flx%forc_prflx(:,:) + p_sfc_flx%forc_evflx(:,:)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=2  ! output print level (1-5, fix)
        CALL dbg_print('UpdSfc: Data SW-flux'      ,p_sfc_flx%forc_swflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data LW-flux'      ,p_sfc_flx%forc_lwflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data Sens.  HF'    ,p_sfc_flx%forc_ssflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data Latent HF'    ,p_sfc_flx%forc_slflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data Total  HF'    ,p_sfc_flx%forc_hflx      ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data Precip.'      ,p_sfc_flx%forc_prflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data Evaporation'  ,p_sfc_flx%forc_evflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: Data Freshw. Flux' ,p_sfc_flx%forc_fwfx      ,str_module,idt_src)
        !---------------------------------------------------------------------

        ! call of sea ice model
        IF (i_sea_ice >= 1) THEN

          Qatm%SWnet  (:,1,:) = p_sfc_flx%forc_swflx(:,:)
          Qatm%LWnet  (:,1,:) = p_sfc_flx%forc_lwflx(:,:)
          Qatm%sens   (:,1,:) = p_sfc_flx%forc_ssflx(:,:)
          Qatm%lat    (:,1,:) = p_sfc_flx%forc_slflx(:,:)
          Qatm%SWnetw (:,:)   = p_sfc_flx%forc_swflx(:,:)
          Qatm%LWnetw (:,:)   = p_sfc_flx%forc_lwflx(:,:)
          Qatm%sensw  (:,:)   = p_sfc_flx%forc_ssflx(:,:)
          Qatm%latw   (:,:)   = p_sfc_flx%forc_slflx(:,:)
          Qatm%dsensdT(:,:,:) = 0.0_wp
          Qatm%dlatdT (:,:,:) = 0.0_wp
          Qatm%dLWdT  (:,:,:) = -4.0_wp * zemiss_def*StBo * (p_ice%Tsurf(:,:,:) + tmelt)**3

          ! sum of flux from sea ice to the ocean is stored in p_sfc_flx%forc_hflx
          !  done in mo_sea_ice:upper_ocean_TS

        ENDIF

      END IF

      IF (temperature_relaxation == 2)  THEN

        !-------------------------------------------------------------------------
        ! Apply temperature relaxation data (record 3) from stationary forcing
        !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
        !  - this is not done for temperature_relaxation=3, since init-data is in Celsius

         p_sfc_flx%forc_tracer_relax(:,:,1) = &
           &  rday1*(ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)-tmelt) + &
           &  rday2*(ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)-tmelt)

      END IF

      IF (irelax_2d_S == 2 .AND. no_tracer >1) THEN

        !-------------------------------------------------------------------------
        ! Apply salinity relaxation data (record ??) from stationary forcing

      !  p_sfc_flx%forc_tracer_relax(:,:,2) = &
      !    &  rday1*(ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,x)-tmelt) + &
      !    &  rday2*(ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,x)-tmelt)
        CALL finish(TRIM(ROUTINE),' irelax_2d_S=2 (reading from flux file) not yet implemented')

      END IF

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) THEN
        WRITE(message_text,'(a,i6,2(a,i4),2(a,f12.8))') 'FLUX time interpolation: jt=',jstep, &
          &  ' mon1=',jmon1,' mon2=',jmon2,' day1=',rday1,' day2=',rday2
        CALL message (' ', message_text)
      END IF
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1)
      CALL dbg_print('UpdSfc: Ext data1-u/mon1'  ,z_c2                     ,str_module,idt_src)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)
      CALL dbg_print('UpdSfc: Ext data1-u/mon2'  ,z_c2                     ,str_module,idt_src)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2)
      CALL dbg_print('UpdSfc: Ext data2-v/mon1'  ,z_c2                     ,str_module,idt_src)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)
      CALL dbg_print('UpdSfc: Ext data2-v/mon2'  ,z_c2                     ,str_module,idt_src)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)
      CALL dbg_print('UpdSfc: Ext data3-t/mon1'  ,z_c2                     ,str_module,idt_src)
      z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)
      CALL dbg_print('UpdSfc: Ext data3-t/mon2'  ,z_c2                     ,str_module,idt_src)
      !---------------------------------------------------------------------

      ! bulk formula are calculated globally

      IF (i_sea_ice >= 1) THEN

        IF (iforc_type == 2 .OR. iforc_type == 5) THEN
          CALL calc_bulk_flux_oce(p_patch, p_as, p_os , Qatm)
          CALL calc_bulk_flux_ice(p_patch, p_as, p_ice, Qatm)
        ENDIF
        
        IF ( no_tracer >= 2 ) THEN
          Tfw(:,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
        ELSE
          Tfw = Tf
        ENDIF

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          CALL ice_fast(i_startidx_c, i_endidx_c, nproma, p_ice%kice, dtime, &
            &   p_ice% Tsurf(:,:,jb),   &
            &   p_ice% T1   (:,:,jb),   &
            &   p_ice% T2   (:,:,jb),   &
            &   p_ice% hi   (:,:,jb),   &
            &   p_ice% hs   (:,:,jb),   &
            &   p_ice% Qtop (:,:,jb),   &
            &   p_ice% Qbot (:,:,jb),   & 
            &   Qatm%SWnet  (:,:,jb),   &
            &   Qatm%lat(:,:,jb) + Qatm%sens(:,:,jb) + Qatm%LWnet(:,:,jb),   & 
            &   Qatm%dlatdT(:,:,jb) + Qatm%dsensdT(:,:,jb) + Qatm%dLWdT(:,:,jb),   & 
            &   Tfw         (:,  jb),   &
            &   Qatm%albvisdir(:,:,jb), &
            &   Qatm%albvisdif(:,:,jb), &
            &   Qatm%albnirdir(:,:,jb), &
            &   Qatm%albnirdif(:,:,jb), &
            &   doy=datetime%yeaday)
        ENDDO

        ! Ocean albedo model
        Qatm%albvisdirw = albedoW
        Qatm%albvisdifw = albedoW
        Qatm%albnirdirw = albedoW
        Qatm%albnirdifw = albedoW

        ! (#slo# 2012-12):
        ! sum of flux from sea ice to the ocean is stored in p_sfc_flx%forc_hflx
        ! diagnosis of 4 parts is stored in p_sfc_flx%forc_swflx/lwflx/ssflx/slflx
        ! this diagnosis is done in mo_sea_ice:upper_ocean_TS
        ! 
        ! under ice the conductive heat flux is not yet stored specifically
        ! the sum forc_hflx is aggregated and stored accordingly which cannot be done here

        ! ATTENTION
        !   ice_slow sets the fluxes in Qatm to zero for a new accumulation in ice_fast
        !   this should be done by the coupler if ice_fast is moved to the atmosphere

        CALL ice_slow(p_patch, p_os, p_ice, Qatm, p_sfc_flx)

      ELSE   !  no sea ice

        ! bulk formula applied to boundary forcing for ocean model:
        !  - no sea ice and no temperature relaxation
        !  - apply net surface heat flux in W/m2

        IF (i_apply_bulk == 1) THEN

          IF (iforc_type == 2 .OR. iforc_type == 5) &
            & CALL calc_bulk_flux_oce(p_patch, p_as, p_os, Qatm)

          temperature_relaxation = 0   !  hack

          DO jb = all_cells%start_block, all_cells%end_block
            CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
            DO jc = i_startidx_c, i_endidx_c

              IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
                p_sfc_flx%forc_swflx(jc,jb) = Qatm%SWnetw(jc,jb) ! net SW radiation flux over water
                p_sfc_flx%forc_lwflx(jc,jb) = Qatm%LWnetw(jc,jb) ! net LW radiation flux over water
                p_sfc_flx%forc_ssflx(jc,jb) = Qatm%sensw (jc,jb) ! Sensible heat flux over water
                p_sfc_flx%forc_slflx(jc,jb) = Qatm%latw  (jc,jb) ! Latent heat flux over water
              ELSE
                p_sfc_flx%forc_swflx(jc,jb) = 0.0_wp
                p_sfc_flx%forc_lwflx(jc,jb) = 0.0_wp
                p_sfc_flx%forc_ssflx(jc,jb) = 0.0_wp
                p_sfc_flx%forc_slflx(jc,jb) = 0.0_wp
              END IF
         
       !      p_sfc_flx%forc_hflx(jc,jb)                 &
       !      & =  Qatm%sensw(jc,jb) + Qatm%latw(jc,jb)  & ! Sensible + latent heat flux over water
       !      & +  Qatm%LWnetw(jc,jb)                    & ! net LW radiation flux over water
       !      & +  Qatm%SWin(jc,jb) * (1.0_wp-albedoW)     ! incoming SW radiation flux

            ENDDO
          ENDDO

          ! for the setup with bulk and without sea ice the threshold for temperature is set to tf
          WHERE (p_os%p_prog(nold(1))%tracer(:,1,:,1) .LT. Tf)
            p_os%p_prog(nold(1))%tracer(:,1,:,1) = Tf
          ENDWHERE

          ! sum of fluxes for ocean boundary condition
          p_sfc_flx%forc_hflx(:,:) = p_sfc_flx%forc_swflx(:,:) + p_sfc_flx%forc_lwflx(:,:) &
            &                      + p_sfc_flx%forc_ssflx(:,:) + p_sfc_flx%forc_slflx(:,:)

        ENDIF

      ENDIF  !  sea ice

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=2  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: Bulk SW-flux'      ,p_sfc_flx%forc_swflx     ,str_module,idt_src)
      CALL dbg_print('UpdSfc: Bulk LW-flux'      ,p_sfc_flx%forc_lwflx     ,str_module,idt_src)
      CALL dbg_print('UpdSfc: Bulk Sens.  HF'    ,p_sfc_flx%forc_ssflx     ,str_module,idt_src)
      CALL dbg_print('UpdSfc: Bulk Latent HF'    ,p_sfc_flx%forc_slflx     ,str_module,idt_src)
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: Bulk Total  HF'    ,p_sfc_flx%forc_hflx      ,str_module,idt_src)
      !---------------------------------------------------------------------

    CASE (FORCING_FROM_FILE_FIELD)                                    !  13
      ! 1) Read field data from file
      ! 2) CALL calc_atm_fluxes_from_bulk (p_patch, p_as, p_os, p_ice, Qatm)
      ! 3) CALL update_sfcflx_from_atm_flx(p_patch, p_as, p_os, p_ice, Qatm, p_sfc_flx)

    CASE (FORCING_FROM_COUPLED_FLUX)                                  !  14
      !  use atmospheric fluxes directly, i.e. avoid call to "calc_atm_fluxes_from_bulk"
      !  and do a direct assignment of atmospheric state to surface fluxes.
      !
      IF ( is_coupled_run() ) THEN
        IF (ltimer) CALL timer_start(timer_coupling)

        time_config%cur_datetime = datetime

        nbr_hor_points = p_patch%n_patch_cells
        nbr_points     = nproma * p_patch%nblks_c
        ALLOCATE(buffer(nbr_points,4))
        buffer(:,:) = 0.0_wp

      !
      !  see drivers/mo_ocean_model.f90:
      !
      !   field_id(1) represents "TAUX"   wind stress component
      !   field_id(2) represents "TAUY"   wind stress component
      !   field_id(3) represents "SFWFLX" surface fresh water flux
      !   field_id(4) represents "SFTEMP" surface temperature
      !   field_id(5) represents "THFLX"  total heat flux
      !   field_id(6) represents "ICEATM" ice temperatures and melt potential
      !
      !   field_id(7) represents "SST"    sea surface temperature
      !   field_id(8) represents "OCEANU" u component of ocean surface current
      !   field_id(9) represents "OCEANV" v component of ocean surface current
      !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
      !
      !
        CALL ICON_cpl_get_nbr_fields ( nbr_fields )
        ALLOCATE(field_id(nbr_fields))
        CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
      !
        field_shape(1) = 1
        field_shape(2) = p_patch%n_patch_cells 
        field_shape(3) = 1

      !
      ! buffer is allocated over nproma only

      !
      ! Send fields from ocean to atmosphere
      ! ------------------------------------
      !
        write_coupler_restart = .FALSE.
      !
      ! SST
        buffer(:,1) = RESHAPE(p_os%p_prog(nold(1))%tracer(:,1,:,1), (/nbr_points /) ) + tmelt
        CALL ICON_cpl_put ( field_id(7), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
        IF ( info == 2 ) write_coupler_restart = .TRUE.
      !
      ! zonal velocity
        buffer(:,1) = RESHAPE(p_os%p_diag%u(:,1,:), (/nbr_points /) )
        CALL ICON_cpl_put ( field_id(8), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
        IF ( info == 2 ) write_coupler_restart = .TRUE.
      !
      ! meridional velocity
        buffer(:,1) = RESHAPE(p_os%p_diag%v(:,1,:), (/nbr_points /) )
        CALL ICON_cpl_put ( field_id(9), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
        IF ( info == 2 ) write_coupler_restart = .TRUE.
      !
      ! Ice thickness, concentration, T1 and T2
        buffer(:,1) = RESHAPE(p_ice%hi  (:,1,:), (/nbr_points /) )
        buffer(:,2) = RESHAPE(p_ice%conc(:,1,:), (/nbr_points /) )
        buffer(:,3) = RESHAPE(p_ice%T1  (:,1,:), (/nbr_points /) )
        buffer(:,4) = RESHAPE(p_ice%T2  (:,1,:), (/nbr_points /) )
        field_shape(3) = 4
        CALL ICON_cpl_put ( field_id(10), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
        IF ( info == 2 ) write_coupler_restart = .TRUE.

        IF ( write_coupler_restart ) CALL icon_cpl_write_restart ( 4, field_id(7:10), ierror )
      !
      ! Receive fields from atmosphere
      ! ------------------------------

      !
      ! Apply wind stress - records 0 and 1 of field_id

      ! zonal wind stress
        field_shape(3) = 1
        CALL ICON_cpl_get ( field_id(1), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
        IF (info > 0 ) THEN
            buffer(nbr_hor_points+1:nbr_points,1) = 0.0_wp
            p_sfc_flx%forc_wind_u(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
            CALL sync_patch_array(sync_c, p_patch, p_sfc_flx%forc_wind_u(:,:))
        ENDIF
      !
      ! meridional wind stress
        CALL ICON_cpl_get ( field_id(2), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
        IF (info > 0 ) THEN
            buffer(nbr_hor_points+1:nbr_points,1) = 0.0_wp
            p_sfc_flx%forc_wind_v(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
            CALL sync_patch_array(sync_c, p_patch, p_sfc_flx%forc_wind_v(:,:))
        ENDIF
      !
      ! Apply freshwater flux - 2 parts, precipitation and evaporation - record 3
      ! prflx(:,:)  total precipitation flux            [m/s]
      ! evflx(:,:)  evaporation flux                    [m/s]
        field_shape(3) = 2
        CALL ICON_cpl_get ( field_id(3), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
        IF (info > 0 ) THEN
            buffer(nbr_hor_points+1:nbr_points,1:2) = 0.0_wp
            p_sfc_flx%forc_prflx(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
            p_sfc_flx%forc_evflx(:,:) = RESHAPE(buffer(:,2),(/ nproma, p_patch%nblks_c /) )
            CALL sync_patch_array(sync_c, p_patch, p_sfc_flx%forc_prflx(:,:))
            CALL sync_patch_array(sync_c, p_patch, p_sfc_flx%forc_evflx(:,:))
            ! sum of fluxes for ocean boundary condition
            p_sfc_flx%forc_fwfx(:,:) = p_sfc_flx%forc_prflx(:,:) + p_sfc_flx%forc_evflx(:,:)
        END IF
      !
      ! Apply surface air temperature - not used with temperature_relaxation= -1 - record 4
        field_shape(3) = 1
        CALL ICON_cpl_get ( field_id(4), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
        IF (info > 0 ) THEN
          buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
          p_sfc_flx%forc_tracer_relax(:,:,1) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
        !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
          p_sfc_flx%forc_tracer_relax(:,:,1) = p_sfc_flx%forc_tracer_relax(:,:,1) - tmelt
        END IF
      !
      ! Apply total heat flux - 4 parts - record 5
      ! p_sfc_flx%swflx(:,:)  ocean short wave heat flux                              [W/m2]
      ! p_sfc_flx%lwflx(:,:)  ocean long  wave, latent and sensible heat fluxes (sum) [W/m2]
        field_shape(3) = 2
        CALL ICON_cpl_get ( field_id(5), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
        IF (info > 0 ) THEN
          buffer(nbr_hor_points+1:nbr_points,1:2) = 0.0_wp
          p_sfc_flx%forc_swflx(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
          p_sfc_flx%forc_lwflx(:,:) = RESHAPE(buffer(:,2),(/ nproma, p_patch%nblks_c /) )
          CALL sync_patch_array(sync_c, p_patch, p_sfc_flx%forc_swflx(:,:))
          CALL sync_patch_array(sync_c, p_patch, p_sfc_flx%forc_lwflx(:,:))
          ! sum of fluxes for ocean boundary condition
          p_sfc_flx%forc_hflx(:,:) = p_sfc_flx%forc_swflx(:,:) + p_sfc_flx%forc_lwflx(:,:)
        ENDIF
      ! p_ice%Qtop(:,:)         Surface melt potential of ice                           [W/m2]
      ! p_ice%Qbot(:,:)         Bottom melt potential of ice                            [W/m2]
      ! p_ice%T1  (:,:)         Temperature of the upper ice layer                      [degC]
      ! p_ice%T2  (:,:)         Temperature of the lower ice layer                      [degC]
        field_shape(3) = 4
        CALL ICON_cpl_get ( field_id(6), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
        IF (info > 0 ) THEN
          buffer(nbr_hor_points+1:nbr_points,1:4) = 0.0_wp
          p_ice%Qtop(:,1,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
          p_ice%Qbot(:,1,:) = RESHAPE(buffer(:,2),(/ nproma, p_patch%nblks_c /) )
          p_ice%T1  (:,1,:) = RESHAPE(buffer(:,3),(/ nproma, p_patch%nblks_c /) )
          p_ice%T2  (:,1,:) = RESHAPE(buffer(:,4),(/ nproma, p_patch%nblks_c /) )
          CALL sync_patch_array(sync_c, p_patch, p_ice%Qtop(:,1,:))
          CALL sync_patch_array(sync_c, p_patch, p_ice%Qbot(:,1,:))
          CALL sync_patch_array(sync_c, p_patch, p_ice%T1  (:,1,:))
          CALL sync_patch_array(sync_c, p_patch, p_ice%T2  (:,1,:))
        END IF

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=1  ! output print level (1-5, fix)
        CALL dbg_print('UpdSfc: CPL: SW-flux'       ,p_sfc_flx%forc_swflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: non-solar flux',p_sfc_flx%forc_lwflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: Total  HF'     ,p_sfc_flx%forc_hflx      ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: Melt-pot. top' ,p_ice%Qtop               ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: Melt-pot. bot' ,p_ice%Qbot               ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: Precip.'       ,p_sfc_flx%forc_prflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: Evaporation'   ,p_sfc_flx%forc_evflx     ,str_module,idt_src)
        CALL dbg_print('UpdSfc: CPL: Freshw. Flux'  ,p_sfc_flx%forc_fwfx      ,str_module,idt_src)
        !---------------------------------------------------------------------

        DEALLOCATE(buffer)
        DEALLOCATE(field_id)      

        IF (ltimer) CALL timer_stop(timer_coupling)

        ! call of sea ice model
        IF (i_sea_ice >= 1) THEN

          Qatm%SWnetw (:,:)   = p_sfc_flx%forc_swflx(:,:)
          Qatm%LWnetw (:,:)   = p_sfc_flx%forc_lwflx(:,:)

          CALL ice_slow(p_patch, p_os, p_ice, Qatm, p_sfc_flx)

          ! sum of flux from sea ice to the ocean is stored in p_sfc_flx%forc_hflx
          !  done in mo_sea_ice:upper_ocean_TS

        ENDIF

      ENDIF ! is_coupled

    CASE (FORCING_FROM_COUPLED_FIELD)                                 !  15
      !1) bulk formula to atmospheric state and proceed as above, the only distinction
      !   to OMIP is that atmospheric info is coming from model rather than file

      CALL message(TRIM(routine), 'STOP: Forcing option 15 not implemented yet' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Forcing option not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION DOES NOT EXIST - TERMINATE')

    END SELECT

    !
    ! After final updating of zonal and merdional components (from file, bulk formula, or coupling)
    ! cartesian coordinates are calculated
    !
    IF (iforc_oce > NO_FORCING) THEN
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
                           & p_sfc_flx%forc_wind_v(jc,jb),&
                           & p_patch%cells%center(jc,jb)%lon,&
                           & p_patch%cells%center(jc,jb)%lat,&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb)         = 0.0_wp
            p_sfc_flx%forc_wind_v(jc,jb)         = 0.0_wp
            p_sfc_flx%forc_wind_cc(jc,jb)%x      = 0.0_wp
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: forcing u'       ,p_sfc_flx%forc_wind_u      ,str_module,idt_src)
      CALL dbg_print('UpdSfc: forcing v'       ,p_sfc_flx%forc_wind_v      ,str_module,idt_src)
      idt_src=2  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: forcing cc%x(1)' ,p_sfc_flx%forc_wind_cc%x(1),str_module,idt_src)
      CALL dbg_print('UpdSfc: forcing cc%x(2)' ,p_sfc_flx%forc_wind_cc%x(2),str_module,idt_src)
      !---------------------------------------------------------------------

    END IF

    !-------------------------------------------------------------------------
    ! Apply temperature relaxation to surface boundary condition
    !  - 2011-12: this is alternative to forcing by fluxes, not in addition

    IF (temperature_relaxation >= 1) THEN

      !  - set minimum temperature to tf (-1.9 deg C) for simple temp-relax
      !  - set to zero on land points

      !z_tmin = -1.0_wp
      z_tmin =tf  !  -1.9 deg C

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            p_sfc_flx%forc_tracer_relax(jc,jb,1) &
              & = max(p_sfc_flx%forc_tracer_relax(jc,jb,1), z_tmin)
          ELSE
            p_sfc_flx%forc_tracer_relax(jc,jb,1) = 0.0_wp
          END IF
        END DO
      END DO

      ! Temperature relaxation activated as boundary condition in vertical Diffusion D:
      !   D = d/dz(K_v*dT/dz)  where
      ! Boundary condition at surface (upper bound of D at center of first layer)
      !   is relaxation to temperature (tau = relaxation constant [1/s] ):
      !   K_v*dT/dz(surf) = Q_T = -dz/tau*(T-T*) [ K*m/s ]
      ! discretized (T* = T_data = relaxation-temperature, forc_tracer_relax):
      !   top_bc_tracer = forc_tracer = -(del_zlev_m+h) / relax_param[s] * (tracer - forc_tracer_relax)
      !
      ! This is equivalent to an additonal forcing term in the tracer equation, i.e. outside
      ! the vertical diffusion, following MITGCM:
      !    F_T  = Q_T/dz = -1/tau * (T-T*) [ K/s ]
      ! when using the sign convention
      !   dT/dt = Operators + F_T
      ! i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature if it is warmer than relaxation data) 
      ! 
      ! Extended boundary condition (relaxation term plus heat flux) is not yet implemented

      ! EFFECTIVE RESTORING PARAMETER: 1.0_wp/(relaxation_param*seconds_per_month)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c


          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            z_relax = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb)) / &
              &       (relaxation_param*seconds_per_month)

            p_sfc_flx%forc_tracer(jc,jb, 1) =                             &
              &          - z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,1)  &
              &                    -p_sfc_flx%forc_tracer_relax(jc,jb,1))
          ELSE
            p_sfc_flx%forc_tracer(jc,jb,1) = 0.0_wp
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      z_c2(:,:) = p_sfc_flx%forc_tracer_relax(:,:,1)
      CALL dbg_print('UpdSfc: Temp-relax'        ,z_c2                    ,str_module,idt_src)
      idt_src=2  ! output print level (1-5, fix)
      z_c2(:,:) = p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1)
      CALL dbg_print('UpdSfc: Temp-difference'   ,z_c2                    ,str_module,idt_src)
      z_c2(:,:) = p_sfc_flx%forc_tracer(:,:,1)
      CALL dbg_print('UpdSfc: T-forc-trac [Km/s]',z_c2                    ,str_module,idt_src)
      !---------------------------------------------------------------------

    ENDIF

    !-------------------------------------------------------------------------
    ! Apply salinity relaxation to surface boundary condition

    IF (irelax_2d_S >= 1 .AND. no_tracer >1) THEN

      ! Salinity relaxation activated as boundary condition in vertical Diffusion D:
      !   D = d/dz(K_v*dS/dz)  where
      ! Boundary condition at surface (upper bound of D at center of first layer)
      !   is relaxation to salinity (tau = relaxation constant [1/s] ):
      !   K_v*dS/dz(surf) = Q_S = -dz/tau*(S-S*) [ psu*m/s ]
      ! discretized (S* = S_data = relaxation-salinity, forc_tracer_relax):
      !   top_bc_tracer = forc_tracer = -(del_zlev_m+h) / relax_param[s] * (tracer - forc_tracer_relax)
      !
      ! This is equivalent to an additonal forcing term in the tracer equation, i.e. outside
      ! the vertical diffusion, following MITGCM:
      !    F_S  = Q_S/dz = -1/tau * (S-S*) [ psu/s ]
      ! when using the sign convention
      !   dS/dt = Operators + F_S
      ! i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity if it is saltier than relaxation data) 
      ! note that the freshwater flux is opposite in sign to F_S, see below,
      ! i.e. fwf >0 for  S-S* >0 (i.e. increasing freshwater flux to decrease the salinity)
      ! 
      ! Mixed boundary conditions (relaxation term plus fluxes) can be included accordingly

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

          !z_relax = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)&
          !          &/(relaxation_param*seconds_per_month)


          z_relax = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb)) / &
            &       (relax_2d_mon_S*seconds_per_month)

              p_sfc_flx%forc_tracer(jc,jb,2) =                             &
             &          - z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,2)  &
             &                    -p_sfc_flx%forc_tracer_relax(jc,jb,2))
          ELSE
            p_sfc_flx%forc_tracer(jc,jb,2) = 0.0_wp
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
!Commented out because of NAG-compiler, PK
!       idt_src=1  ! output print level (1-5, fix)
!       z_c2(:,:) = p_sfc_flx%forc_tracer_relax(:,:,2)
!       CALL dbg_print('UpdSfc: Sal.-relax'        ,z_c2                    ,str_module,idt_src)
!       idt_src=2  ! output print level (1-5, fix)
!       z_c2(:,:) = p_sfc_flx%forc_tracer_relax(:,:,2)-p_os%p_prog(nold(1))%tracer(:,1,:,2)
!       CALL dbg_print('UpdSfc: Sal.-difference'   ,z_c2                    ,str_module,idt_src)
!       z_c2(:,:) = p_sfc_flx%forc_tracer(:,:,2)
!       CALL dbg_print('UpdSfc: S-forc-trac [Km/s]',z_c2                    ,str_module,idt_src)
      !---------------------------------------------------------------------

    ENDIF

    ! Heat flux diagnosed for all ocean only relaxation cases
    IF (temperature_relaxation >= 1) THEN

      ! Heat flux diagnosed for relaxation cases, see above
      !   Q_s = Rho*Cp*Q_T  [W/m2]  with density Rho and Cp specific heat capacity
      ! where
      !   Q_T = K_v*dT/dz(surf) = Q_s/Rho/Cp  [K*m/s]

      p_sfc_flx%forc_hflx(:,:) = p_sfc_flx%forc_tracer(:,:,1) * rho_ref * clw

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc:T-relax-hflx [W/m2]',p_sfc_flx%forc_hflx     ,str_module,idt_src)
      !---------------------------------------------------------------------

    END IF

    ! Freshwater flux diagnosed
    IF (irelax_2d_S >= 1 .AND. no_tracer >1) THEN

      ! Freshwater flux at surface W_s diagnosed for relaxation cases (see Griffies)
      !   W_s = -Q_S / S_0   [m/s]  with reference salinity at surface, 
      !   which is set to reference value S_ref=35 psu to avoid instability for low salinity
      ! where
      !   Q_S = K_v*dS/dz(surf) = -W_s*S_ref  [psu*m/s]
      ! from above

      ! p_sfc_flx%forc_fwfx(:,:) = -p_sfc_flx%forc_tracer(:,:,2) / sal_ref
      ! now in m/month for diagnosis
        p_sfc_flx%forc_fwfx(:,:) = -p_sfc_flx%forc_tracer(:,:,2) / sal_ref * 2.592e6_wp

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc:S-relax-fwfx[m/mon]',p_sfc_flx%forc_fwfx     ,str_module,idt_src)
      !---------------------------------------------------------------------

    END IF

    !-------------------------------------------------------------------------
    ! Apply net surface heat flux to boundary condition
    !  - heat flux is applied alternatively to temperature relaxation for coupling
    !  - also done if sea ice model is used since forc_hflx is set in mo_sea_ice
    !  - with OMIP-forcing and sea_ice=0 we need temperature_relaxation=1
    !    since there is no forc_hflx over open water when using OMIP-forcing
    !  - i_apply_bulk=1 provides net surface heat flux globally
    !

    IF (temperature_relaxation == -1 .OR. i_sea_ice >= 1 .OR. i_apply_bulk == 1) THEN

      ! Heat flux boundary condition for diffusion
      !   D = d/dz(K_v*dT/dz)  where
      ! Boundary condition at surface (upper bound of D at center of first layer)
      !   is calculated from net surface heat flux Q_s [W/m2]
      !   which is calculated by the atmosphere (coupled) or read from flux file (see above)
      !   Q_s = Rho*Cp*Q_T  with density Rho and Cp specific heat capacity
      !   K_v*dT/dz(surf) = Q_T = Q_s/Rho/Cp  [K*m/s]
      ! discretized:
      !   top_bc_tracer = forc_tracer = forc_hflx / (rho_ref*clw)

      p_sfc_flx%forc_tracer(:,:,1) = p_sfc_flx%forc_hflx(:,:) / (rho_ref*clw)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: T-forc-hflx[W/m2]' ,p_sfc_flx%forc_hflx     ,str_module,idt_src)
      idt_src=3  ! output print level (1-5, fix)
      z_c2(:,:) = p_sfc_flx%forc_tracer(:,:,1)
      CALL dbg_print('UpdSfc:T-forc-trac[K*m/s]' ,z_c2                    ,str_module,idt_src)
      !---------------------------------------------------------------------

    END IF

    !-------------------------------------------------------------------------
    ! Apply net surface freshwater flux to boundary condition

    IF (irelax_2d_S == -1 .AND. no_tracer >1) THEN

      ! Salinity boundary condition in vertical Diffusion D, see above

      p_sfc_flx%forc_tracer(:,:,2) = -p_sfc_flx%forc_fwfx(:,:) * sal_ref

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('UpdSfc: S-forc-fwfx[m/s]'  ,p_sfc_flx%forc_fwfx     ,str_module,idt_src)
      idt_src=3  ! output print level (1-5, fix)
      z_c2(:,:) = p_sfc_flx%forc_tracer(:,:,2)
      CALL dbg_print('UpdSfc:S-forc-trac[psu*m/s]',z_c2                    ,str_module,idt_src)
      !---------------------------------------------------------------------

    END IF

  END SUBROUTINE update_sfcflx

  !-------------------------------------------------------------------------
  !
  !> Takes thermal calc_atm_fluxes_from_bulk to calculate surface fluxes for ocean forcing:
  !!  heat, freshwater and momentum.
  !!  not active or tested yet (2012/08)
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011). Originally written by D. Notz.
  !
  SUBROUTINE update_sfcflx_from_atm_flx(p_patch_3D, p_as, p_os, p_ice, Qatm, p_sfc_flx)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_atmos_for_ocean),      INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),    INTENT(IN)    :: p_os
    TYPE (t_sea_ice),             INTENT (IN)   :: p_ice
    TYPE (t_atmos_fluxes),        INTENT (INOUT):: Qatm
    TYPE(t_sfc_flx)                             :: p_sfc_flx

    !Local variables 
    REAL(wp) :: z_rho_w = 1.22_wp  !near surface air density [kg/m^3] cf. Large/Yeager, sect 4.1, p.17
    REAL(wp) :: z_C_d0, z_C_d1, z_C_d
    REAL(wp) :: z_norm, z_v, z_relax

    INTEGER :: jc, jb, i
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp):: z_evap        (nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp):: z_Q_freshwater(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_sfcflx_from_atm_flx'
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    all_cells => p_patch%cells%all

    !Relaxation parameter from namelist for salinity.
    z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO i = 1, p_ice%kice
          !surface heat forcing as sum of sensible, latent, longwave and shortwave heat fluxes
          IF (p_ice% hi(jc,jb,i) > 0._wp)THEN

    !  ATTENTION - forc_tracer is INCORRECT here
    !   - forc_tracer is boundary condition in vertical diffusion equation [K*m/s]
    !   - forc_hflx is net surface heat flux [W/m2]
    !       p_sfc_flx%forc_tracer(jc,jb,1)               &
            p_sfc_flx%forc_hflx(jc,jb)                   &
              & =  Qatm%sens(jc,jb,i) + Qatm%lat(jc,jb,i)& ! Sensible + latent heat
              !                                              flux at ice surface
              & +  Qatm%LWnet(jc,jb,i)                   & ! net LW radiation flux over ice surface
              & +  Qatm%bot(jc,jb,i)                       ! Ocean heat flux at ice bottom 
                                                           ! liquid/solid  precipitation rate
            !                                                are zero

            !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
            z_evap(jc,jb) = Qatm%lat(jc,jb,i)/(als*z_rho_w)

          ELSE

    !       p_sfc_flx%forc_tracer(jc,jb,1)             &
            p_sfc_flx%forc_hflx(jc,jb)                 &
            & =  Qatm%sensw(jc,jb) + Qatm%latw(jc,jb)  & ! Sensible + latent heat flux over water
            & +  Qatm%LWnetw(jc,jb)                    & ! net LW radiation flux over water
            & +  Qatm%SWnetw(jc,jb)                      ! net SW radiation flux ove water
                                                         ! liquid/solid  precipitation rate are zero

           !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
            z_evap(jc,jb) = Qatm%latw(jc,jb)/(alv*z_rho_w)
          ENDIF
        END DO

        !calculate surface freshwater flux       
        !following MPI-OM as described in Marsland et al, formula (63)-(65)

        !calculate evaporation from latent heat flux and latent heat of vaporisation
        !This is (63) in Marsland et al.
        !+River runoff +glacial meltwater
        z_Q_freshwater(jc,jb) = (Qatm%rpreci(jc,jb) + Qatm%rprecw(jc,jb)) -  z_evap(jc,jb)  

        !Now the freshwater flux calculation is finished; this is (65) in Marsland et al.
        !Relaxation of top layer salinity to observed salinity
        !
        !  Attention, check consistency in the model:
        !   - salinity relaxation is here in addition to the formulation at the end of update_sfcflx
        !   - also, according to (65) of Marsland, there is a bug below:
        !     multiplication with S1 (tracer(2)) is missing
        !   - has to be checked and merged with salinity boundary condition in update_sfcflx
        !
        p_sfc_flx%forc_tracer(jc,jb,2) =                 &
          & (p_patch_3D%p_patch_1D(1)%del_zlev_m(1)+z_Q_freshwater(jc,jb)) &
          & /p_patch_3D%p_patch_1D(1)%del_zlev_m(1)                        &  !  * tracer(jc,1,jb,2)
          & +z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,2)-p_sfc_flx%forc_tracer_relax(jc,jb,2))


        !calculate wind stress    
        z_norm = sqrt(p_as%u(jc,jb)*p_as%u(jc,jb)+p_as%v(jc,jb)*p_as%v(jc,jb))

        !calculate drag coefficient for wind following 
        ! Kara, Rochford, Hurlburt, Air-Sea Flux Estimates And the 1997-1998 Enso Event
        ! Boundary-Layer Meteorology, 103, 439-458 (2002)
        !
        z_v = MAX(2.5_wp, MIN(p_as%fu10(jc,jb),32.5_wp))

        z_C_d0 = 1.0E-3_wp*(0.692_wp+0.071_wp*z_v-0.00070_wp*z_norm)
        z_C_d1 = 1.0E-3_wp*(0.083_wp-0.0054_wp*z_v-0.000093_wp*z_norm)
        z_C_d  = z_C_d0 + z_C_d1*(p_as%tafo(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1))

        !write(*,*)'final wind stress coeff',z_C_d
        p_sfc_flx%forc_wind_u(jc,jb) = z_rho_w*z_C_d*z_norm &
          &                            *(p_as%u(jc,jb)- p_os%p_diag%u(jc,1,jb))

        p_sfc_flx%forc_wind_v(jc,jb) = z_rho_w*z_C_d*z_norm &
          &                            *(p_as%v(jc,jb) - p_os%p_diag%v(jc,1,jb))
   
      END DO
    END DO

    IF(temperature_relaxation==1)THEN

       p_sfc_flx%forc_tracer(:,:, 1)=  z_relax                                    &
       & *( p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1) )

    ENDIF

  END SUBROUTINE update_sfcflx_from_atm_flx
  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE update_sfcflx_analytical(p_patch_3D, p_os, p_sfc_flx)

  TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
  TYPE(t_hydro_ocean_state)                   :: p_os
  TYPE(t_sfc_flx)                             :: p_sfc_flx
  !
  ! local variables
  INTEGER :: jc, jb
  INTEGER :: i_startidx_c, i_endidx_c
  !INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  !INTEGER :: rl_start_c, rl_end_c

  REAL(wp) :: zonal_str
  REAL(wp) :: z_lat, z_lon, z_lat_deg
  REAL(wp) :: z_forc_period = 1.0_wp !=1.0: single gyre
                                     !=2.0: double gyre
                                     !=n.0: n-gyre
  REAL(wp) :: y_length               !basin extension in y direction in degrees
  REAL(wp) :: z_T_init(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
  REAL(wp) :: z_perlat, z_perlon, z_permax, z_perwid, z_relax, z_dst
  INTEGER  :: z_dolic
  REAL(wp) :: z_temp_max, z_temp_min, z_temp_incr
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_ho_sfcflx'
  !-------------------------------------------------------------------------
  TYPE(t_subset_range), POINTER :: all_cells
  TYPE(t_patch), POINTER :: p_patch
  !-----------------------------------------------------------------------  
  p_patch         => p_patch_3D%p_patch_2D(1)
  !-------------------------------------------------------------------------
  all_cells => p_patch%cells%all


 ! #slo#  Stationary forcing is moved to mo_oce_forcing:init_ho_forcing

    SELECT CASE (itestcase_oce)

    CASE(30,32,27)

      CALL message(TRIM(routine), &
      &  'Testcase (30,32,27) - stationary lat/lon wind forcing &
      &and eventually relax. to T perturbation')
      y_length = basin_height_deg * deg2rad
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

             ! #slo# Warning: s.th. more missing?
             z_lat = p_patch%cells%center(jc,jb)%lat
             z_lon = p_patch%cells%center(jc,jb)%lon

             zonal_str = wstress_coeff*cos(z_forc_period*pi*z_lat-y_length/y_length)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(1) = wstress_coeff*zonal_str*sin(z_lon)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(2) = wstress_coeff*zonal_str*cos(z_lon)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(3) = 0.0_wp
 
             CALL cvec2gvec(p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(3),&
                          & z_lon, z_lat,                      &
                          & p_sfc_flx%forc_wind_u(jc,jb),      &
                          & p_sfc_flx%forc_wind_v(jc,jb))
           ELSE
             p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
             p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
           ENDIF 
        END DO
      END DO
   !  write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)

     IF(no_tracer>=1.AND.temperature_relaxation/=0)THEN

        y_length = basin_height_deg * deg2rad
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c

            IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

              z_T_init(jc,jb) = 20.0_wp- p_patch_3D%p_patch_1D(1)%zlev_m(1)*15.0_wp/4000.0_wp

              z_lat = p_patch%cells%center(jc,jb)%lat
              z_lon = p_patch%cells%center(jc,jb)%lon
 
              ! Add temperature perturbation at new values
              z_perlat = basin_center_lat + 0.1_wp*basin_height_deg
              z_perlon = basin_center_lon + 0.1_wp*basin_width_deg 
              z_permax  = 0.1_wp
              z_perwid  =  10.0_wp

              z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

             z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
             IF (z_dolic > MIN_DOLIC) THEN

               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)

               IF(z_dst<=5.0_wp*deg2rad)THEN
                 z_T_init = z_T_init &
                 &        + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                 &        * sin(pi*p_patch_3D%p_patch_1D(1)%zlev_m(1)/4000.0_wp)
                 !   write(*,*)'z init',jc,jb,p_os%p_prog(nold(1))%tracer(jc,1,jb,1),&
                 !   &z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                 !   & * sin(pi*v_base%zlev_m(1)/4000.0_wp)
               ENDIF
               ! up to here z_init is identically initialized than temperature

               !add local cold perturbation 
               IF(z_dst<=10.5_wp*deg2rad)THEN
                 z_T_init(jc,jb)= z_T_init(jc,jb) - exp(-(z_dst/(z_perwid*deg2rad))**2)
               ENDIF

               p_sfc_flx%forc_tracer_relax(jc,jb,1)=z_T_init(jc,jb)

               p_sfc_flx%forc_tracer(jc,jb, 1)=  z_relax   &          
               & *( p_sfc_flx%forc_tracer_relax(jc,jb,1)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )

               ! write(123,*)'forcing',jc,jb,&
               ! &( p_sfc_flx%forc_tracer_relax(jc,jb,1)    &
               ! & -p_os%p_prog(nold(1))%tracer(jc,1,jb,1)),&
               ! &p_sfc_flx%forc_tracer_relax(jc,jb,1),&
               ! &p_sfc_flx%forc_tracer(jc,jb, 1)
             END IF
           ELSE
             p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
             p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
           ENDIF 
        END DO
      END DO

 !  write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
 !  & minval(p_sfc_flx%forc_tracer_relax)
 !  write(*,*)'max/min-tracer-flux',maxval(p_sfc_flx%forc_tracer),&
 !  & minval(p_sfc_flx%forc_tracer)
 !  write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
 !                                & minval(p_sfc_flx%forc_tracer(:,:,1))
    ENDIF
! ! ! !-----------Old version of Forcing--------------------------------------------------
! ! !!------------Please retain, its also interesting------------------------------------
!!----------------An old version of init corresponds to this forcing--------------------
! !    CASE(32)
! !       CALL message(TRIM(routine), 'Testcase (32): Apply stationary wind forcing' )
! !       y_length = basin_height_deg * deg2rad
! !       DO jb = i_startblk_c, i_endblk_c    
! !         CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
! !          &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
! !         DO jc = i_startidx_c, i_endidx_c
! !           z_lat = p_patch%cells%center(jc,jb)%lat
! !           z_lon = p_patch%cells%center(jc,jb)%lon
! !           IF(v_base%lsm_c(jc,1,jb)<=sea_boundary)THEN
! !             zonal_str = wstress_coeff*cos(z_forc_period*pi*z_lat-y_length/y_length)
! !             p_sfc_flx%forc_wind_cc(jc,jb)%x(1) = wstress_coeff*zonal_str*sin(z_lon)
! !             p_sfc_flx%forc_wind_cc(jc,jb)%x(2) = wstress_coeff*zonal_str*cos(z_lon)
! !             p_sfc_flx%forc_wind_cc(jc,jb)%x(3) = 0.0_wp
! !             CALL cvec2gvec(p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
! !                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
! !                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(3),&
! !                          & z_lon, z_lat,                      &
! !                          & p_sfc_flx%forc_wind_u(jc,jb),      &
! !                          & p_sfc_flx%forc_wind_v(jc,jb))
! !             ! Add temperature perturbation at new values
! !            z_perlat = basin_center_lat + 0.1_wp*basin_height_deg!             !45.5_wp
! !            z_perlon =  0.1_wp*basin_width_deg                                 !4.5_wp
! !            z_permax  = 10.0_wp!20.1_wp
! !            z_perwid  =  5.0_wp!1.5_wp
! !            z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
! ! 
! !             z_dolic = v_base%dolic_c(jc,jb)
! !             IF (z_dolic > 0) THEN
! !               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)
! ! 
! !               !init temperature
! !               z_T_init(jc,jb) = 20.0_wp&
! !               & - v_base%zlev_i(1)*15.0_wp/v_base%zlev_i(z_dolic+1)
! ! 
! !                !add local hot perturbation 
! ! !              IF(z_dst<=3.5_wp*deg2rad)THEN
! !                 z_T_init(jc,jb)= z_T_init(jc,jb)  &
! !                 &   + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
! !                 &   * sin(pi*v_base%zlev_m(1)/v_base%zlev_i(z_dolic+1))
! ! !              ENDIF
! !               !Add local cold perturbation
! !               !IF(z_dst<=5.0_wp*deg2rad)THEN
! !               z_T_init(jc,jb) = z_T_init(jc,jb)     &
! !               &   - z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2)
! !               p_sfc_flx%forc_tracer_relax(jc,jb,1)=z_T_init(jc,jb)
! !               p_sfc_flx%forc_tracer(jc,jb, 1)=z_relax*v_base%del_zlev_i(1)*&
! !               &          (p_sfc_flx%forc_tracer_relax(jc,jb,1)&
! !               &         -p_os%p_prog(nold(1))%tracer(jc,1,jb,1))
! !               !ENDIF 
! !             END IF
! !   ! write(*,*)'Danilovs Wind', jc,jb,p_sfc_flx%forc_wind_cc(jc,jb)%x(1:2), &
! !   ! &p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
! !            ELSE
! !              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
! !              p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
! !              p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
! !            ENDIF 
! !         END DO
! !       END DO
! !       write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)
! !       write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
! !                                   & minval(p_sfc_flx%forc_tracer(:,:,1))
! ! ! !-----------End of Old version of Forcing-------------------------------------------

    CASE (33)
      IF(iforc_oce/=10)THEN 
      y_length = basin_height_deg * deg2rad
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

             ! #slo# Warning: s.th. more missing?
             z_lat = p_patch%cells%center(jc,jb)%lat
             z_lon = p_patch%cells%center(jc,jb)%lon

             zonal_str = wstress_coeff*cos(z_forc_period*pi*z_lat-y_length/y_length)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(1) = wstress_coeff*zonal_str*sin(z_lon)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(2) = wstress_coeff*zonal_str*cos(z_lon)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(3) = 0.0_wp
 
             CALL cvec2gvec(p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(3),&
                          & z_lon, z_lat,                      &
                          & p_sfc_flx%forc_wind_u(jc,jb),      &
                          & p_sfc_flx%forc_wind_v(jc,jb))
           ELSE
             p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
             p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
           ENDIF 
        END DO
      END DO
  !   write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)
      ENDIF
      IF(temperature_relaxation>=1)THEN
      ! CALL message(TRIM(routine), &
      !   &  'Testcase (33): stationary temperature relaxation - latitude dependent')
        z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

        p_sfc_flx%forc_tracer(:,:, 1) = z_relax*( p_sfc_flx%forc_tracer_relax(:,:,1) &
          &                                      -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF

 !  write(*,*)'max/min-tracer-diff',&
 !  &maxval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1)),&
 !  & minval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1))

 !  write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
 !  & minval(p_sfc_flx%forc_tracer_relax)
 !  write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
 !                                & minval(p_sfc_flx%forc_tracer(:,:,1))
    CASE(51)

      CALL message(TRIM(routine), &
      &  'Testcase (51) - stationary lat/lon wind forcing &
      &and eventually relax. to T perturbation')
      y_length = basin_height_deg * deg2rad
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

             ! #slo# Warning: s.th. more missing?
             z_lat = p_patch%cells%center(jc,jb)%lat
             z_lon = p_patch%cells%center(jc,jb)%lon

             p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff * &
             & cos(z_forc_period*pi*(z_lat-y_length)/y_length) 

             p_sfc_flx%forc_wind_v(jc,jb)= 0.0_wp

             CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
                           & p_sfc_flx%forc_wind_v(jc,jb),&
                           & p_patch%cells%center(jc,jb)%lon,&
                           & p_patch%cells%center(jc,jb)%lat,&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
           ELSE
             p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
             p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
           ENDIF 
        END DO
      END DO
  !   write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)

      IF(temperature_relaxation>=1)THEN

        z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

        z_temp_max  = 30.5_wp
        z_temp_min  = 0.5_wp
        z_temp_incr = (z_temp_max-z_temp_min)/(n_zlev-1.0_wp)

      !Add horizontal variation
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lat_deg = z_lat*rad2deg

            IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

              z_temp_max     =0.01_wp*(z_lat_deg-basin_center_lat)*(z_lat_deg-basin_center_lat)
              z_T_init(jc,jb)=30.5_wp

              z_T_init(jc,jb)&
              &=z_T_init(jc,jb)*exp(-z_temp_max/basin_height_deg)
            ELSE
              z_T_init(jc,jb)=0.0_wp
            ENDIF
        END DO
      END DO
      p_sfc_flx%forc_tracer_relax(:,:,1)=z_T_init(:,:)

      p_sfc_flx%forc_tracer(:,:, 1) = z_relax*( p_sfc_flx%forc_tracer_relax(:,:,1) &
          &                                      -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF

  ! CASE (43)
  !   ! no forcing applied
  !   CONTINUE

  ! CASE DEFAULT
  !   CALL message(TRIM(routine), 'STOP: Analytical Forcing for this testcase not implemented' )
  !   CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')
    END SELECT

  END SUBROUTINE update_sfcflx_analytical

  !-------------------------------------------------------------------------
  !>
  !! Read ocean forcing data from netcdf
  !!
  !! Read ocean forcing data for NCEP or other forcing
  !! This routine reads annual data sets of length iforc_len
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2012-02-17)
  !!
  !!
  SUBROUTINE read_forc_data_oce (p_patch, ext_data, no_set)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER,       INTENT(IN)            :: no_set          !  no of set in file to be read

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_oce_bulk:read_forc_data_oce'

    CHARACTER(filename_max) :: ncep_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, i_lev, i_cell_type, no_cells, no_tst, jtime, jt !, jc, jb
    INTEGER :: ncid, dimid,mpi_comm
    INTEGER :: i_start(2),i_count(2), jcells

    REAL(wp):: z_flux(nproma,p_patch%nblks_c,iforc_len)  ! set length is iforc_len, 3rd dimension
    REAL(wp):: z_c   (nproma,iforc_len,p_patch%nblks_c)  ! 2nd dimension is iforc_len
    !TYPE (t_keyword_list), POINTER :: keywords => NULL()

    !-------------------------------------------------------------------------

    !  READ NCEP FORCING

    !-------------------------------------------------------------------------

    !CALL message (TRIM(routine), 'start')

    IF (iforc_oce == 12) THEN

    !DO jg = 1,n_dom
      jg = 1

      i_lev       = p_patch%level
      i_cell_type = p_patch%cell_type

      IF(my_process_is_stdio()) THEN
        !
        WRITE (ncep_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-flux.nc'

        !ncep_file=TRIM('/pool/data/ICON/external/iconR2B04-flux.nc')
        CALL message( TRIM(routine),'Ocean NCEP forcing flux file is: '//TRIM(ncep_file) )
        INQUIRE (FILE=ncep_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'NCEP forcing flux file is not found.')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ncep_file), NF_NOWRITE, ncid))
        !CALL message( TRIM(routine),'Ocean NCEP flux file opened for read' )

        !
        ! get and check number of cells in ncep data
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

        IF(p_patch%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in NCEP flux file do not match.')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid(ncid, 'time', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst))
        !
        ! check - s.b.

      ENDIF
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      CALL p_bcast(no_tst, p_io, mpi_comm)

      !-------------------------------------------------------
      !
      ! Read 12 monthly NCEP data sets for triangle centers using 4-dim routine
      !
      !-------------------------------------------------------

      jcells = p_patch%n_patch_cells  !  global dimension
      jtime  = iforc_len              !  time period to read (not yet)


      ! provide NCEP fluxes for sea ice (interface to ocean)
      ! 1:  'stress_x': zonal wind stress       [m/s]
      ! 2:  'stress_y': meridional wind stress  [m/s]
      ! 3:  'SST"     : sea surface temperature [K]

      ! zonal wind stress
      !write(0,*) ' ncep set 1: dimensions:',p_patch%n_patch_cells_g, p_patch%n_patch_cells, &
      ! &  iforc_len, nproma, p_patch%nblks_c
      !CALL read_netcdf_data (ncid, 'stress_x', p_patch%n_patch_cells_g,      &
      !  &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
      !  &                    iforc_len, z_flx2(:,:,:))
      !write(0,*) ' READ_FORC, READ 1: first data sets: stress-x, block=5, index=1,5:'
      !do jt=1,jtime
      !  write(0,*) 'jt=',jt,' val:',(z_flx2(jc,jt,5),jc=1,5)
      !enddo

      ! start-pointer and length of pointer for reading data:
      ! start: first set (1,1); second year (1,jtime+1)
      i_start(1) = 1
      i_start(2) = jtime*(no_set-1) + 1  ! position pointer to set no_set
      i_count(1) = jcells                ! length of pointer, dim 1 of z_dummy_array
      i_count(2) = jtime                 ! length of pointer, dim 2 of z_dummy_array

      idt_src=2  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) THEN
        !
        WRITE(message_text,'(A,I6,A)')  'Ocean NCEP flux file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )

        WRITE(message_text,'(4(A,I4))')  'NCEP data set: length = ',jtime, &
          &   '; no. of set =',no_set,                                     &
          &   '; pos. of ptr =', i_start(2)
        CALL message( TRIM(routine), TRIM(message_text) )
      END IF

      CALL read_netcdf_data (ncid, 'stress_x', p_patch%n_patch_cells_g,      &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))


      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,1) = z_flux(:,:,jt)
      END DO

      ! meridional wind stress
      CALL read_netcdf_data (ncid, 'stress_y', p_patch%n_patch_cells_g,      &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,2) = z_flux(:,:,jt)
      END DO

      ! SST
      CALL read_netcdf_data (ncid, 'SST', p_patch%n_patch_cells_g,           &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,3) = z_flux(:,:,jt)
      END DO

 !    ! Read complete NCEP data sets for focing ocean model (iforc_type=5)
 !    ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
 !    ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
 !    ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
 !    ! 7:  fclou(:,:),  &  ! Fractional cloud cover
 !    ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
 !    ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

      ! 2m-temperature
      CALL read_netcdf_data (ncid, 'temp_2m', p_patch%n_patch_cells_g,       &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,4) = z_flux(:,:,jt)
      END DO

      ! 2m dewpoint temperature
      CALL read_netcdf_data (ncid, 'dpt_temp_2m', p_patch%n_patch_cells_g,   &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,5) = z_flux(:,:,jt)
      END DO

      ! Scalar wind
      CALL read_netcdf_data (ncid, 'scalar_wind', p_patch%n_patch_cells_g,   &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,6) = z_flux(:,:,jt)
      END DO

      ! cloud cover
      CALL read_netcdf_data (ncid, 'cloud', p_patch%n_patch_cells_g,         &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,7) = z_flux(:,:,jt)
      END DO

      ! sea level pressure
      CALL read_netcdf_data (ncid, 'pressure', p_patch%n_patch_cells_g,      &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,8) = z_flux(:,:,jt)
      END DO

      ! total solar radiation
      CALL read_netcdf_data (ncid, 'tot_solar', p_patch%n_patch_cells_g,     &
        &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
        &                    jtime, i_start, i_count, z_flux(:,:,:))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,9) = z_flux(:,:,jt)
      END DO

      ! precipitation
  !   CALL read_netcdf_data (ncid, 'precip', p_patch%n_patch_cells_g,        &
  !     &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,10) = z_flux(:,:,jt)
  !   END DO

      ! evaporation or downward surface LW flux
  !   CALL read_netcdf_data (ncid, 'evap', p_patch%n_patch_cells_g,          &
  !     &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,11) = z_flux(:,:,jt)
  !   END DO
  !   CALL read_netcdf_data (ncid, 'dlwrf', p_patch%n_patch_cells_g,         &
  !     &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,11) = z_flux(:,:,jt)
  !   END DO

      ! runoff
  !   CALL read_netcdf_data (ncid, 'runoff', p_patch%n_patch_cells_g,        &
  !     &                    p_patch%n_patch_cells, p_patch%cells%glb_index, &
  !     &                    jtime, i_start, i_count, z_flux(:,:,:))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,12) = z_flux(:,:,jt)
  !   END DO


      !
      ! close file
      !
      IF(my_process_is_stdio()) CALL nf(nf_close(ncid))

    !ENDDO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,1)
      CALL dbg_print('ReadFc: NCEP: stress-x'    ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,2)
      CALL dbg_print('ReadFc: NCEP: stress-y'    ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,3)
      CALL dbg_print('ReadFc: NCEP: SST'         ,z_c                     ,str_module,idt_src)
      idt_src=4  ! output print level (1-5, fix)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4)
      CALL dbg_print('ReadFc: NCEP: temp_2m'     ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5)
      CALL dbg_print('ReadFc: NCEP: dpt_temp_2m' ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6)
      CALL dbg_print('ReadFc: NCEP: scalar_wind' ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7)
      CALL dbg_print('ReadFc: NCEP: cloudiness'  ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8)
      CALL dbg_print('ReadFc: NCEP: pressure'    ,z_c                     ,str_module,idt_src)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9)
      CALL dbg_print('ReadFc: NCEP: total solar' ,z_c                     ,str_module,idt_src)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,10)
    ! CALL dbg_print('ReadFc: NCEP: precip.'     ,z_c                     ,str_module,idt_src)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,11)
    ! CALL dbg_print('ReadFc: NCEP: evaporation' ,z_c                     ,str_module,idt_src)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,12)
    ! CALL dbg_print('ReadFc: NCEP: runoff'      ,z_c                     ,str_module,idt_src)
      !---------------------------------------------------------------------

      idt_src=2  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) &
        & CALL message( TRIM(routine),'Ocean NCEP fluxes for external data read' )

    END IF ! iforc_oce=12

  END SUBROUTINE read_forc_data_oce

  !-------------------------------------------------------------------------

  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_oce_bulk netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


END MODULE mo_oce_bulk
