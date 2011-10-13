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
USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling
USE mo_io_units,            ONLY: filename_max
USE mo_mpi,                 ONLY: p_pe, p_io, p_bcast
USE mo_datetime,            ONLY: t_datetime
USE mo_ext_data,            ONLY: ext_data
USE mo_grid_config,         ONLY: nroot
USE mo_ocean_nml,           ONLY: iforc_oce, iforc_omip, iforc_len, itestcase_oce,    &
  &           no_tracer, n_zlev, basin_center_lat, basin_center_lon, basin_width_deg, &
  &                               basin_height_deg, relaxation_param, wstress_coeff,  &
  &                               i_bc_veloc_top, temperature_relaxation,             &
  &                               NO_FORCING, ANALYT_FORC, FORCING_FROM_FILE_FLUX,    &
  &                               FORCING_FROM_FILE_FIELD, FORCING_FROM_COUPLED_FLUX, &
  &                               FORCING_FROM_COUPLED_FIELD, i_dbg_oce, i_sea_ice
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch
USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base
USE mo_exception,           ONLY: finish, message, message_text
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_physical_constants,  ONLY: rho_ref, sfc_press_bar, lsub, lvap, lfreez, cpa, emiss, &
  &                               fr_fac, stefbol, rgas, tmelt, tf, cw, rhoi, rhow, rhos, albedoW
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary,MIN_DOLIC
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec, cvec2gvec
USE mo_sea_ice,             ONLY: t_sea_ice
! USE mo_oce_forcing,         ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
USE mo_sea_ice,             ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
                                  calc_atm_fluxes_from_bulk, set_ice_albedo, set_ice_temp, &
                                  sum_fluxes
USE mo_oce_thermodyn,       ONLY: convert_insitu2pot_temp_func
USE mo_oce_index,           ONLY: print_mxmn, ipl_src
USE mo_master_control,      ONLY: is_coupled_run
USE mo_icon_cpl_exchg,      ONLY: ICON_cpl_put, ICON_cpl_get
USE mo_icon_cpl_def_field,  ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids

IMPLICIT NONE

! required for reading netcdf files
INCLUDE 'netcdf.inc'

CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface
PUBLIC  :: update_sfcflx
!PUBLIC  :: calc_atm_fluxes_from_bulk

! private implementation
PRIVATE :: update_sfcflx_analytical

PRIVATE


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
  SUBROUTINE update_sfcflx(p_patch, p_os, p_as, p_ice, Qatm, QatmAve, p_sfc_flx, jstep, datetime)

  TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
  TYPE(t_hydro_ocean_state)           :: p_os
  TYPE(t_atmos_for_ocean)             :: p_as
  TYPE(t_atmos_fluxes)                :: Qatm, QatmAve
  TYPE(t_sea_ice)                     :: p_ice
  TYPE(t_sfc_flx)                     :: p_sfc_flx
  INTEGER, INTENT(IN)                 :: jstep
  TYPE(t_datetime), INTENT(INOUT)     :: datetime
  !
  ! local variables
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_sfcflx'
  INTEGER  :: jmon, njday, jdmon, jmon1, jmon2!, jdays
  INTEGER  :: jc, jb
  INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER  :: rl_start_c, rl_end_c
  REAL(wp) :: z_tmin, z_relax, rday1, rday2
  REAL(wp) :: z_c(nproma,n_zlev,p_patch%nblks_c)

  ! Local declarations for coupling:
  INTEGER               :: info, ierror !< return values form cpl_put/get calls
  INTEGER               :: nbr_hor_points ! = inner and halo points
  INTEGER               :: nbr_points     ! = nproma * nblks
  INTEGER               :: nbr_fields
  INTEGER, ALLOCATABLE  :: field_id(:)
  INTEGER               :: field_shape(3)
  REAL(wp), ALLOCATABLE :: buffer(:,:)

  !-------------------------------------------------------------------------

  rl_start_c = 1
  rl_end_c   = min_rlcell
  i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

  SELECT CASE (iforc_oce)

  CASE (NO_FORCING)                !  10

  ! CALL message(TRIM(routine), 'No  forcing applied' )
    CONTINUE

  CASE (ANALYT_FORC)

    CALL update_sfcflx_analytical(p_patch, p_os, p_sfc_flx)

  CASE (FORCING_FROM_FILE_FLUX)    !  12

    ! To Do: read forcing file in chunks
    ! Check if file should be read (new chunk if timecriterion is met, or completely)
    ! 

    !-------------------------------------------------------------------------
    ! Applying annual forcing read from file in mo_ext_data:
    !  - stepping daily in monthly data (preliminary solution)

    !  calculate day and month
    jmon  = datetime%month         ! current month
    jdmon = datetime%day           ! day in month

    !jdmon = mod(jdays+1,30)-1     ! no of days in month

    ! To Do: use fraction of month for interpolation
    !frcmon= datetime%monfrc       ! fraction of month
    !rday1 = frcmon+0.5_wp
    !rday2 = 1.0_wp-rday1
    !IF (rday1 > 1.0_wp)  THEN
    !  rday2=rday1
    !  rday1=1.0_wp-rday1
    !END IF

    njday = int(86400._wp/dtime)  ! no of timesteps per day

    !
    ! interpolate monthly forcing-data daily:
    !
    IF (iforc_len == 12)  THEN

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
      ! - use rday1, rday2, jmon1 = jmon2 = yeaday for controling correct day in year
      ! - no interpolation applied, 
      ! - datetime%yeaday may be 366 for a leap year, e.g. year 2000
      jmon1 = datetime%yeaday
      jmon2 = jmon1
      rday1 = 1.0_wp
      rday2 = 0.0_wp
    
    END IF

    !
    ! OMIP data read in mo_ext_data into variable ext_data
    !
    IF (iforc_omip >= 1)  THEN

      ! provide OMIP fluxes for wind stress forcing
      ! 1:  wind_u(:,:)   !  'stress_x': zonal wind stress       [m/s]
      ! 2:  wind_v(:,:)   !  'stress_y': meridional wind stress  [m/s]

      ! ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean
      p_sfc_flx%forc_wind_u(:,:) = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,1) + &
        &                          rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,1)
      p_sfc_flx%forc_wind_v(:,:) = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,2) + &
        &                          rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,2)
     
      ! now done in top_bound_cond_horz_veloc
      !IF(iswm_oce == 1)THEN
      !  z_scale = v_base%del_zlev_m(1)*rho_ref
      !ELSEIF(iswm_oce /= 1)THEN
      !  z_scale = rho_ref
      !ENDIF
      !p_sfc_flx%forc_wind_u=p_sfc_flx%forc_wind_u/z_scale
      !p_sfc_flx%forc_wind_v=p_sfc_flx%forc_wind_v/z_scale

    END IF

    IF (iforc_omip == 2) THEN

    !-------------------------------------------------------------------------
      ! provide OMIP fluxes for sea ice (interface to ocean)
      ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
      ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
      ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
      ! 7:  fclou(:,:),  &  ! Fractional cloud cover
      ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
      ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

      p_as%tafo(:,:)  = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,4) + &
        &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,4)
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      p_as%tafo(:,:)  = p_as%tafo(:,:) - tmelt
      p_as%ftdew(:,:) = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,5) + &
        &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,5)
      p_as%fu10(:,:)  = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,6) + &
        &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,6)
      p_as%fclou(:,:) = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,7) + &
        &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,7)
      p_as%pao(:,:)   = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,8) + &
        &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,8)
      !  don't - change units to mb/hPa
      p_as%pao(:,:)   = p_as%pao(:,:) !* 0.01
      p_as%fswr(:,:)  = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,9) + &
        &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,9)

    END IF

    IF (iforc_omip == 3) THEN

      !-------------------------------------------------------------------------
      ! Apply surface heat and freshwater fluxes

   !  p_as%sfch(:,:)  = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,4) + &
   !    &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,4)
   !  p_as%sfcf(:,:)  = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,4) + &
   !    &               rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,4)

    END IF

    IF (temperature_relaxation == 2 .OR. temperature_relaxation == 3)  THEN

      !-------------------------------------------------------------------------
      ! Apply Temperature relaxation data (array 3) from stationary forcing
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)

       p_sfc_flx%forc_tracer_relax(:,:,1) = &
         &  rday1*(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,3)-tmelt) + &
         &  rday2*(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,3)-tmelt)

    END IF   !  temperature relaxation
     
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
        &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary)THEN
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

    ipl_src=3  ! output print level (1-5, fix)
    IF (i_dbg_oce >= ipl_src) THEN
      WRITE(message_text,'(a,i6,2(a,i2),2(a,f12.8))') 'FLUX time interpolation: jt=',jstep, &
        &  ' mon1=',jmon1,' mon2=',jmon2,' day1=',rday1,' day2=',rday2
      CALL message (' ', message_text)
    END IF
    z_c(:,1,:)=ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,1)
    CALL print_mxmn('Ext data1 (u) mon1',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:)=ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,1)
    CALL print_mxmn('Ext data1 (u) mon2',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:)=ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,2)
    CALL print_mxmn('Ext data2 (v) mon1',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:)=ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,2)
    CALL print_mxmn('Ext data2 (v) mon2',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:)=ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,3)
    CALL print_mxmn('Ext data3 (t) mon1',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:)=ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,3)
    CALL print_mxmn('Ext data3 (t) mon2',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)

    IF (i_sea_ice == 1) THEN
      IF (iforc_omip == 2) &
        CALL calc_atm_fluxes_from_bulk (p_patch, p_as, p_os, p_ice, Qatm)

      CALL update_sfcflx_from_atm_flx(p_patch, p_as, p_os, p_ice, Qatm, p_sfc_flx)
      CALL set_ice_albedo(p_patch,p_ice)
      CALL set_ice_temp(p_patch,p_ice,Qatm)
      CALL sum_fluxes(Qatm, QatmAve)
    ENDIF

  CASE (FORCING_FROM_FILE_FIELD)                                    !  13
    ! 1) Read field data from file
    ! 2) CALL calc_atm_fluxes_from_bulk (p_patch, p_as, p_os, p_ice, Qatm)
    ! 3)CALL update_sfcflx_from_atm_flx(p_patch, p_as, p_os, p_ice, Qatm, p_sfc_flx)

  CASE (FORCING_FROM_COUPLED_FLUX)                                  !  14
    !  use atmospheric fluxes directly, i.e. avoid call to "calc_atm_fluxes_from_bulk"
    !  and do a direct assignment of atmospheric state to surface fluxes.
    !
    IF ( is_coupled_run() ) THEN
      IF (ltimer) CALL timer_start(timer_coupling)

!       CALL message(TRIM(routine), "executing OCEAN coupling")
      nbr_hor_points = p_patch%n_patch_cells
      nbr_points     = nproma * p_patch%nblks_c
      ALLOCATE(buffer(nbr_points,1))
      
    !
    !  see drivers/mo_atmo_model.f90:
    !
    !   field_id(1) represents "TAUX"   wind stress component
    !   field_id(2) represents "TAUY"   wind stress component
    !   field_id(3) represents "SFWFLX" surface fresh water flux
    !   field_id(4) represents "SFTEMP" surface temperature
    !   field_id(5) represents "THFLX"  total heat flux
    !
    !   field_id(6) represents "SST"    sea surface temperature
    !   field_id(7) represents "OCEANU" u component of ocean surface current
    !   field_id(8) represents "OCEANV" v component of ocean surface current
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
    ! SST:
      buffer(:,1) = RESHAPE(p_os%p_prog(nold(1))%tracer(:,1,:,1), (/nbr_points /) )  + 273.15_wp 
      CALL ICON_cpl_put ( field_id(6), field_shape, buffer, ierror )
    !
    ! zonal wind
      buffer(:,1) = RESHAPE(p_os%p_diag%u(:,1,:), (/nbr_points /) )
      CALL ICON_cpl_put ( field_id(7), field_shape, buffer, ierror )
    !
    ! meridional wind
      buffer(:,1) = RESHAPE(p_os%p_diag%v(:,1,:), (/nbr_points /) )
      CALL ICON_cpl_put ( field_id(8), field_shape, buffer, ierror )
    !
    ! Receive fields from atmosphere
    ! ------------------------------
    !
    ! zonal wind stress
      CALL ICON_cpl_get ( field_id(1), field_shape, buffer, info, ierror )
      p_sfc_flx%forc_wind_u(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
    !
    ! meridional wind stress
      CALL ICON_cpl_get ( field_id(2), field_shape, buffer, info, ierror )
      p_sfc_flx%forc_wind_v(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
    !
    ! freshwater flux
      CALL ICON_cpl_get ( field_id(3), field_shape, buffer, info, ierror )
      p_sfc_flx%forc_fwfx(:,:) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
    !
    ! surface temperature
      CALL ICON_cpl_get ( field_id(4), field_shape, buffer, info, ierror )
      p_sfc_flx%forc_tracer_relax(:,:,1) = RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )
    !
    ! total heat flux
      CALL ICON_cpl_get ( field_id(5), field_shape, buffer, info, ierror )
      ! #slo# why accumulated? LL Accumulation will be done by the coupler
      p_sfc_flx%forc_hflx(:,:) = &! p_sfc_flx%forc_hflx(:,:) + &
        & RESHAPE(buffer(:,1),(/ nproma, p_patch%nblks_c /) )

      DEALLOCATE(buffer)
      DEALLOCATE(field_id)      

      IF (ltimer) CALL timer_stop(timer_coupling)
  
    ENDIF ! iforc_omip

  CASE (FORCING_FROM_COUPLED_FIELD)                                 !  15
    !1) bulk formula to atmospheric state and proceed as above, the only distinction
    !   to OMIP is that atmospheric info is coming from model rather than file

    CALL message(TRIM(routine), 'STOP: Forcing option 15 not implemented yet' )
    CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')

  CASE DEFAULT

    CALL message(TRIM(routine), 'STOP: Forcing option not implemented' )
    CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION DOES NOT EXIST - TERMINATE')

  END SELECT

  IF (iforc_oce > NO_FORCING) THEN
    ipl_src=1  ! output print level (1-5, fix)
    z_c(:,1,:)=p_sfc_flx%forc_wind_u(:,:)
    CALL print_mxmn('update forcing u',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:)=p_sfc_flx%forc_wind_v(:,:)
    CALL print_mxmn('update forcing v',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
  END IF

  !-------------------------------------------------------------------------
  ! Apply temperature relaxation to boundary condition

  IF (temperature_relaxation >= 2) THEN

    !  - set minimum temperature to tf (-1.9 deg C) for simple temp-relax
    !  - set to zero on land points

    !z_tmin = -1.0_wp
    z_tmin = tf  !  -1.9 deg C

    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
        &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        IF (v_base%lsm_oce_c(jc,1,jb) <= sea_boundary) THEN
          p_sfc_flx%forc_tracer_relax(jc,jb,1) &
            & = max(p_sfc_flx%forc_tracer_relax(jc,jb,1), z_tmin)
        ELSE
          p_sfc_flx%forc_tracer_relax(jc,jb,1) = 0.0_wp
        END IF
      END DO
    END DO

    ! Temperature relaxation:
    ! #slo# corrected formula: Diffusion
    !   D = d/dz(K_v*dT/dz)  where
    ! Boundary condition at surface (upper bound of D at center of first layer)
    !   is relaxation to temperature:
    !   K_v*dT/dz(surf) = -Q_t = -dz/Tau*(T*-T) = dz/Tau*(T-T*)  [K*m/s]
    ! discretized:
    !   top_bc_tracer = forc_tracer = del_zlev_m+h / relax_param[s] * (tracer - forc_tracer_relax)

    !z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
    !z_relax = v_base%del_zlev_m(1)/(relaxation_param*2.592e6_wp)

    DO jb = i_startblk_c, i_endblk_c    
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
       &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        z_relax = (v_base%del_zlev_m(1)+p_os%p_prog(nold(1))%h(jc,jb)) / &
          &       (relaxation_param*2.592e6_wp)

        IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
          p_sfc_flx%forc_tracer(jc,jb, 1) =                             &
          &          - z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,1)  &
          &                    -p_sfc_flx%forc_tracer_relax(jc,jb,1))
        ELSE
          !p_sfc_flx%forc_tracer_relax(jc,jb,1) = 0.0_wp
          p_sfc_flx%forc_tracer(jc,jb,1) = 0.0_wp
        ENDIF
      END DO
    END DO

    !write(0,*) 'relpar, z_rel:',relaxation_param,z_relax
    ipl_src=1  ! output print level (1-5, fix)
    z_c(:,1,:)=p_sfc_flx%forc_tracer_relax(:,:,1)
    CALL print_mxmn('update temp-relax',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    ipl_src=2  ! output print level (0-5, fix)
    z_c(:,1,:) = p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1)
    CALL print_mxmn('Temp-difference',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    z_c(:,1,:) = p_sfc_flx%forc_tracer(:,:,1)
    CALL print_mxmn('T-forc-tracer-flux',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)

  ENDIF

  !-------------------------------------------------------------------------
  ! Apply net surface heat flux to boundary condition

  IF (temperature_relaxation == -1) THEN

    ! Heat flux boundary condition for diffusion
    !   D = d/dz(K_v*dT/dz)  where
    ! Boundary condition at surface (upper bound of D at center of first layer)
    !   is calculated from net surface heat flux Q_s [W/m2]
    !   which is calculated by the atmosphere (coupled) or read from flux file (see above)
    !   Q_s = Rho*Cp*Q_t  mit Cp specific heat capacity
    !   K_v*dT/dz(surf) = -Q_t = Q_s/Rho/Cp  [K*m/s]
    ! discretized:
    !   top_bc_tracer = forc_tracer = forc_hflx / (rho_ref*cw)

    p_sfc_flx%forc_tracer(:,:,1) = p_sfc_flx%forc_hflx(:,:) / (rho_ref*cw)

    ipl_src=1  ! output print level (1-5, fix)
    z_c(:,1,:) = p_sfc_flx%forc_hflx(:,:)
    CALL print_mxmn('T-forc-nshflx',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)
    ipl_src=2  ! output print level (1-5, fix)
    z_c(:,1,:) = p_sfc_flx%forc_tracer(:,:,1)
    CALL print_mxmn('T-forc-tracer-flux',1,z_c(:,:,:),n_zlev,p_patch%nblks_c,'bul',ipl_src)

  END IF

  END SUBROUTINE update_sfcflx
  !-------------------------------------------------------------------------
  !
  !> Takes thermal calc_atm_fluxes_from_bulk to calculate atmospheric surface fluxes:
  !  heat, freshwater and momentum.
  !! @par Revision History
  !! Rewritten by Einar Olason, MPI-M (2011) based on the previous version and
  !upper_ocean_TS by Dirk Notz
  !
  SUBROUTINE update_sfcflx_from_atm_flx(ppatch, p_as, p_os, ice,QatmAve, p_sfc_flx)
  TYPE(t_patch),            INTENT(IN)     :: ppatch 
  TYPE(t_atmos_for_ocean),  INTENT(IN)     :: p_as
  TYPE(t_hydro_ocean_state),INTENT(INOUT)  :: p_os
  TYPE(t_sea_ice),          INTENT (INOUT) :: ice
  TYPE(t_atmos_fluxes),     INTENT (INOUT) :: QatmAve
  TYPE(t_sfc_flx),          INTENT (INOUT) :: p_sfc_flx

!!Local Variables
  REAL(wp), DIMENSION (nproma, ppatch%nblks_c,ice%kice) ::    &
   draft         ! position of ice-ocean interface below sea level        [m] 
  
  REAL(wp), DIMENSION (nproma, ppatch%nblks_c) ::   & 
   draftAve,    &! average draft of sea ice within a grid cell            [m]
   zUnderIceOld,&! water in upper ocean grid cell below ice (prev. time)  [m]
   heatOceI,    &! heat flux into ocean through formerly ice covered areas[W/m�]
   heatOceW,    &! heat flux into ocean through open water areas          [W/m�]
   delHice,     &! average change in ice thickness within a grid cell     [m]
   snowiceave,  &! average snow to ice conversion within a grid cell      [m]
   evap,        &! evaporated water                                       [m]
   preci,       &! solid precipitation                                    [m]
   precw,       &! liquid precipitation                                   [m]
   z_norm, z_v, z_C_d0, z_C_d1, z_C_d ! drag coefficients for wind

  ! Needs work with FB_BGC_OCE etc.
   REAL(wp)         :: swsum 
   REAL(wp),POINTER :: sao_top(:,:)
   REAL(wp)         :: z_relax
   REAL(wp)         :: z_rho_w = 1.22_wp  !near surface air density [kg/m^3] cf. Large/Yeager, sect 4.1, p.17
!  !-------------------------------------------------------------------------------

  swsum = 0.0_wp
  sao_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

  ! Calculate change in water level 'zo' from liquid and solid precipitation and
  ! evaporation
  precw           (:,:)   = QatmAve% rprecw (:,:) * dtime
  preci           (:,:)   = QatmAve% rpreci (:,:) * dtime
  evap            (:,:)   = (QatmAve% latw(:,:)/ Lvap * dtime * &
                            sum(ice%conc(:,:,:), 3) +           &
                            sum(ice%evapwi(:,:,:) * ice% conc(:,:,:), 3)) /rhow
  p_os%p_prog(nold(1))%h(:,:) = p_os%p_prog(nold(1))%h(:,:) +  precw + preci - evap

  ! Calculate average draft and thickness of water underneath ice in upper ocean
  ! grid box
  zUnderIceOld    (:,:)   = ice%zUnderIce
  draft           (:,:,:) = (rhos * ice%hs + rhoi * ice%hi) / rhow
  draftave        (:,:)   = sum(draft(:,:,:) * ice%conc(:,:,:),3)
  ice%zUnderIce   (:,:)   = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(:,:) - draftave(:,:) 
 
  ! Calculate average change in ice thickness and the snow-to-ice conversion 
  Delhice         (:,:)   = sum((ice% hi(:,:,:) - ice% hiold(:,:,:))*          &
                            ice%conc(:,:,:),3)
  snowiceave      (:,:)   = sum(ice%snow_to_ice(:,:,:) * ice% conc(:,:,:),3)
 

  ! Calculate heat input through formerly ice covered and through open water
  ! areas
  heatOceI        (:,:)   = sum(ice% heatOceI(:,:,:) * ice% conc(:,:,:),3)
  heatOceW        (:,:)   = (QatmAve%SWin(:,:) * (1.0_wp-albedoW) * (1.0_wp-swsum) +    &
                            QatmAve%LWnetw(:,:) + QatmAve%sensw(:,:)+         &
                            QatmAve%latw(:,:))  *  (1.0_wp-sum(ice%conc,3))

  ! Change temperature of upper ocean grid cell according to heat fluxes
!  p_os%p_prog(nold(1))%tracer(:,1,:,1) = p_os%p_prog(nold(1))%tracer(:,1,:,1)&
!                                       & + dtime*(heatOceI + heatOceW) /               &
!                                       & (cw*rhow * ice%zUnderIce)
! TODO: should we also divide with ice%zUnderIce / ( v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) ) ?
  p_sfc_flx%forc_tracer(:,:,1) = (heatOceI + heatOceW) / (cw*rhow)

! TODO: 
!  ! Temperature change of upper ocean grid cell due  to melt-water inflow and
!  ! precipitation
!  p_os%p_prog(nold(1))%tracer(:,1,:,1) = (p_os%p_prog(nold(1))%tracer(:,1,:,1)&
!                          &*zUnderIceOld &
!                          &+ precw*p_as%tafo + preci*0.0_wp + &                             !!!!!!!!!Dirk: times 0.0 ????
!                          &  sum(ice%surfmeltT * ice%surfmelt * ice%conc,3)) / & 
!                          &  (zUnderIceOld + sum(ice%surfmelt*ice%conc,3) +    &
!                          &  precw + preci)
!
!  ! Change salinity of upper ocean grid box from ice growth/melt, snowice
!  ! formation and precipitation
!  p_os%p_prog(nold(1))%tracer(:,1,:,2) = p_os%p_prog(nold(1))%tracer(:,1,:,2)  &
!                                       & + (Delhice(:,:)*rhoi - snowiceave(:,:)*rhos)/rhow *  &
!                                       & MIN(Sice, sao_top(:,:)) / ice%zUnderIce(:,:)
!
  !heatabs         (:,:)   = swsum * QatmAve% SWin * (1 - ice%concsum)

  IF(temperature_relaxation==1)THEN

  !Relaxation parameter from namelist for salinity.
     z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

     p_sfc_flx%forc_tracer(:,:, 1)=  z_relax                                    &
     & *( p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1) )

  ENDIF

  !calculate wind stress    
  z_norm = sqrt(p_as%u*p_as%u+p_as%v*p_as%v)

  !calculate drag coefficient for wind following 
  ! Kara, Rochford, Hurlburt, Air-Sea Flux Estimates And the 1997-1998 Enso Event
  ! Boundary-Layer Meteorology, 103, 439-458 (2002)
  !
  z_v = MAX(2.5_wp, MIN(p_as%fu10,32.5_wp))

  z_C_d0 = 1.0E-3_wp*(0.692_wp+0.071_wp*z_v-0.00070_wp*z_norm)
  z_C_d1 = 1.0E-3_wp*(0.083_wp-0.0054_wp*z_v-0.000093_wp*z_norm)
  z_C_d  = z_C_d0 + z_C_d1*(p_as%tafo-p_os%p_prog(nold(1))%tracer(:,1,:,1))
  !write(*,*)'final wind stress coeff',z_C_d
  p_sfc_flx%forc_wind_u = z_rho_w*z_C_d*z_norm&
                               &*(p_as%u- p_os%p_diag%u(:,1,:))

  p_sfc_flx%forc_wind_v = z_rho_w*z_C_d*z_norm&
                               &*(p_as%v - p_os%p_diag%v(:,1,:))

  END SUBROUTINE update_sfcflx_from_atm_flx  

  !-------------------------------------------------------------------------
  !
  !> Takes thermal calc_atm_fluxes_from_bulk to calculate atmospheric surface fluxes:
  !  heat, freshwater and momentum.
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011). Originally written by D. Notz.
  !
!  SUBROUTINE update_sfcflx_from_atm_flx(ppatch, p_as, p_os, p_ice,Qatm, p_sfc_flx)
!  TYPE(t_patch),                INTENT(in)    :: ppatch
!  TYPE(t_atmos_for_ocean),      INTENT(IN)    :: p_as
!  TYPE(t_hydro_ocean_state),    INTENT(IN)    :: p_os
!  TYPE (t_sea_ice),             INTENT (IN)   :: p_ice
!  TYPE (t_atmos_fluxes),        INTENT (INOUT):: Qatm
!  TYPE(t_sfc_flx)                             :: p_sfc_flx
!
!  !Local variables 
!  REAL(wp) :: z_rho_w = 1.22_wp  !near surface air density [kg/m^3] cf. Large/Yeager, sect 4.1, p.17
!  REAL(wp) :: z_C_d0, z_C_d1, z_C_d
!  REAL(wp) :: z_norm, z_v, z_relax
!
!  INTEGER :: jc, jb, i
!  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
!  INTEGER :: rl_start_c, rl_end_c
!  REAL(wp):: z_evap(nproma,ppatch%nblks_c)
!  REAL(wp):: z_Q_freshwater(nproma,ppatch%nblks_c)
!  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_sfcflx_from_atm_flx'
!  !-------------------------------------------------------------------------
!  CALL message(TRIM(routine), 'start' )
!
!  rl_start_c = 1
!  rl_end_c   = min_rlcell
!
!  i_startblk_c = ppatch%cells%start_blk(rl_start_c,1)
!  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)
!
!  !Relaxation parameter from namelist for salinity.
!  z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
!
!  DO jb = i_startblk_c, i_endblk_c
!    CALL get_indices_c( ppatch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!    &                   rl_start_c, rl_end_c)
!    DO jc = i_startidx_c, i_endidx_c
!      DO i = 1, p_ice%kice
!        !surface heat forcing as sum of sensible, latent, longwave and shortwave heat fluxes
!        IF (p_ice% isice(jc,jb,i))THEN
!
!          p_sfc_flx%forc_tracer(jc,jb,1)             &
!          & =  Qatm%sens(jc,jb,i) + Qatm%lat(jc,jb,i)& ! Sensible + latent heat flux at ice surface
!          & +  Qatm%LWnet(jc,jb,i)                   & ! net LW radiation flux over ice surface
!          & +  Qatm%bot(jc,jb,i)                       ! Ocean heat flux at ice bottom 
!                                                       ! liquid/solid  precipitation rate are zero
!
!          !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
!          z_evap(jc,jb) = Qatm%lat(jc,jb,i)/(Lsub*z_rho_w)
!
!        ELSEIF(.NOT.p_ice% isice(jc,jb,i))THEN
!
!          p_sfc_flx%forc_tracer(jc,jb,1)             &
!          & =  Qatm%sensw(jc,jb) + Qatm%latw(jc,jb)  & ! Sensible + latent heat flux over water
!          & +  Qatm%LWnetw(jc,jb)                    & ! net LW radiation flux over water
!          & +  Qatm%SWin(jc,jb)                        ! incoming SW radiation flux
!                                                       ! liquid/solid  precipitation rate are zero
!
!         !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
!          z_evap(jc,jb) = Qatm%latw(jc,jb)/(Lvap*z_rho_w)
!        ENDIF
!      END DO
!
!      !calculate surface freshwater flux       
!      !following MPI-OM as described in Marsland et al, formula (63)-(65)
!
!      !calculate evaporation from latent heat flux and latent heat of vaporisation
!      !This is (63) in Marsland et al.
!      z_Q_freshwater(jc,jb) = (Qatm%rpreci(jc,jb) + Qatm%rprecw(jc,jb)) -  z_evap(jc,jb) !+River runof +glacial meltwater
!
!      !Now the freshwater flux calculation is finished; this is (65) in Marslkand et al.
!      !Relaxation of top layer salinity to observed salinity
!      !
!      p_sfc_flx%forc_tracer(jc,jb,2) = &
!      &(v_base%del_zlev_m(1)+z_Q_freshwater(jc,jb))&
!      &/v_base%del_zlev_m(1)                       &
!      & + z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,2) - p_sfc_flx%forc_tracer_relax(jc,jb,2))
!
!
!      !calculate wind stress    
!      z_norm = sqrt(p_as%u(jc,jb)*p_as%u(jc,jb)+p_as%v(jc,jb)*p_as%v(jc,jb))
!
!      !calculate drag coefficient for wind following 
!      ! Kara, Rochford, Hurlburt, Air-Sea Flux Estimates And the 1997-1998 Enso Event
!      ! Boundary-Layer Meteorology, 103, 439-458 (2002)
!      !
!      z_v = MAX(2.5_wp, MIN(p_as%fu10(jc,jb),32.5_wp))
!
!      z_C_d0 = 1.0E-3_wp*(0.692_wp+0.071_wp*z_v-0.00070_wp*z_norm)
!      z_C_d1 = 1.0E-3_wp*(0.083_wp-0.0054_wp*z_v-0.000093_wp*z_norm)
!      z_C_d  = z_C_d0 + z_C_d1*(p_as%tafo(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1))
!      !write(*,*)'final wind stress coeff',z_C_d
!      p_sfc_flx%forc_wind_u(jc,jb) = z_rho_w*z_C_d*z_norm&
!                                   &*(p_as%u(jc,jb)- p_os%p_diag%u(jc,1,jb))
!
!      p_sfc_flx%forc_wind_v(jc,jb) = z_rho_w*z_C_d*z_norm&
!                                   &*(p_as%v(jc,jb) - p_os%p_diag%v(jc,1,jb))
! 
!    END DO
!  END DO
!
!  IF(temperature_relaxation==1)THEN
!
!     p_sfc_flx%forc_tracer(:,:, 1)=  z_relax                                    &
!     & *( p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1) )
!
!  ENDIF
!
!  END SUBROUTINE update_sfcflx_from_atm_flx  
  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE update_sfcflx_analytical(p_patch, p_os, p_sfc_flx)

  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch
  TYPE(t_hydro_ocean_state)             :: p_os  
  TYPE(t_sfc_flx)                       :: p_sfc_flx
  !
  ! local variables
  INTEGER :: jc, jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  REAL(wp) :: zonal_str
  REAL(wp) :: z_lat, z_lon
  REAL(wp) :: z_forc_period = 1.0_wp !=1.0: single gyre
                                     !=2.0: double gyre
                                     !=n.0: n-gyre 
  REAL(wp) :: y_length               !basin extension in y direction in degrees
  REAL(wp) :: z_T_init(nproma,p_patch%nblks_c)
  REAL(wp) :: z_perlat, z_perlon, z_permax, z_perwid, z_relax, z_dst
  INTEGER  :: z_dolic
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_ho_sfcflx'
  !-------------------------------------------------------------------------
  rl_start_c   = 1
  rl_end_c     = min_rlcell
  i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)


 ! #slo#  Stationary forcing is moved to mo_oce_forcing:init_ho_forcing

    SELECT CASE (itestcase_oce)

    CASE(30,32,27)

      CALL message(TRIM(routine), &
      &  'Testcase (30,32,27) - stationary lat/lon wind forcing &
      &and eventually relax. to T perturbation')
      y_length = basin_height_deg * deg2rad
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN

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
      write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)

     IF(no_tracer>=1.AND.temperature_relaxation/=0)THEN

        y_length = basin_height_deg * deg2rad
        DO jb = i_startblk_c, i_endblk_c    
          CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
           &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)

          DO jc = i_startidx_c, i_endidx_c

            IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN

              z_T_init(jc,jb) = 20.0_wp- v_base%zlev_m(1)*15.0_wp/4000.0_wp

              z_lat = p_patch%cells%center(jc,jb)%lat
              z_lon = p_patch%cells%center(jc,jb)%lon
 
              ! Add temperature perturbation at new values
              z_perlat = basin_center_lat + 0.1_wp*basin_height_deg
              z_perlon = basin_center_lon + 0.1_wp*basin_width_deg 
              z_permax  = 0.1_wp
              z_perwid  =  10.0_wp

              z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

             z_dolic = v_base%dolic_c(jc,jb)
             IF (z_dolic > MIN_DOLIC) THEN

               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)

               IF(z_dst<=5.0_wp*deg2rad)THEN
                 z_T_init = z_T_init &
                 &        + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                 &        * sin(pi*v_base%zlev_m(1)/4000.0_wp)
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

    write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
    & minval(p_sfc_flx%forc_tracer_relax)
    write(*,*)'max/min-tracer-flux',maxval(p_sfc_flx%forc_tracer),&
    & minval(p_sfc_flx%forc_tracer)
    write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
                                  & minval(p_sfc_flx%forc_tracer(:,:,1))
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
! !           IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
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
      IF(temperature_relaxation>=1)THEN
      ! CALL message(TRIM(routine), &
      !   &  'Testcase (33): stationary temperature relaxation - latitude dependent')
        z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
        
        p_sfc_flx%forc_tracer(:,:, 1) = z_relax*( p_sfc_flx%forc_tracer_relax(:,:,1) &
          &                                      -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF

    write(*,*)'max/min-tracer-diff',&
    &maxval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1)),&
    & minval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1))

    write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
    & minval(p_sfc_flx%forc_tracer_relax)
    write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
                                  & minval(p_sfc_flx%forc_tracer(:,:,1))
  ! CASE (43)
  !   ! no forcing applied
  !   CONTINUE

  ! CASE DEFAULT
  !   CALL message(TRIM(routine), 'STOP: Analytical Forcing for this testcase not implemented' )
  !   CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')
    END SELECT

  END SUBROUTINE update_sfcflx_analytical

  !-------------------------------------------------------------------------

  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


END MODULE mo_oce_bulk
