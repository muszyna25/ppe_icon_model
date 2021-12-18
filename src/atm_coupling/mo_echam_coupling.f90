!>
!! @brief Interface between ECHAM physics and the ocean, through a coupler
!!
!! @author Marco Giorgetta (MPI-M)
!!
!! @par Revision History
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

MODULE mo_echam_coupling

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nonhydro_types      ,ONLY: t_nh_diag
  USE mo_echam_phy_memory    ,ONLY: prm_field
  USE mo_ccycle_config       ,ONLY: ccycle_config
                                
  USE mo_parallel_config     ,ONLY: nproma

  USE mo_run_config          ,ONLY: ltimer, ico2, nlev
  USE mo_timer,               ONLY: timer_start, timer_stop,                &
       &                            timer_coupling_put, timer_coupling_get, &
       &                            timer_coupling_1stget
  USE mo_echam_sfc_indices   ,ONLY: iwtr, iice, ilnd, nsfc_type
  USE mo_echam_phy_config    ,ONLY: echam_phy_config

  USE mo_sync                ,ONLY: sync_c, sync_patch_array

  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif

  USE mo_parallel_config     ,ONLY: nproma

  USE mo_coupling_config     ,ONLY: is_coupled_run
  USE mo_atmo_coupling_frame ,ONLY: lyac_very_1st_get, nbr_inner_cells,     &
    &                               mask_checksum, field_id
  USE mo_exception           ,ONLY: warning, finish, message

  USE mo_yac_finterface      ,ONLY: yac_fput, yac_fget,                     &
    &                               COUPLING, OUT_OF_BOUND

  USE mo_util_dbg_prnt       ,ONLY: dbg_print
  USE mo_dbg_nml             ,ONLY: idbg_mxmn, idbg_val
  USE mo_physical_constants  ,ONLY: amd, amco2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_echam_ocean

  CHARACTER(len=12)           :: str_module = 'InterFaceOce'  ! Output of module for 1 line debug

CONTAINS

  !>
  !! SUBROUTINE interface_echam_ocean -- the interface between
  !! ECHAM physics and the ocean, through a coupler
  !!
  !! This subroutine is called in the time loop of the ICONAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! <li> tendency of the prognostic varibles induced by adiabatic dynamics;
  !! <li> time step;
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients.
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables caused by
  !! the parameterisations.
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE interface_echam_ocean( p_patch , pt_diag)

    ! Arguments

    TYPE(t_patch), TARGET, INTENT(INOUT)    :: p_patch
    TYPE(t_nh_diag), TARGET, INTENT(INOUT)  :: pt_diag

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells  ! = inner and halo points
    INTEGER               :: jg             ! grid index
    INTEGER               :: n              ! nproma loop count
    INTEGER               :: nn             ! block offset
    INTEGER               :: i_blk          ! block loop count
    INTEGER               :: nlen           ! nproma/npromz
    INTEGER               :: info, ierror   !< return values from cpl_put/get calls
    INTEGER               :: no_arr         !  no of arrays in bundle for put/get calls

    REAL(wp), PARAMETER   :: dummy = 0.0_wp

    REAL(wp)              :: scr(nproma,p_patch%alloc_cell_blocks)
    REAL(wp)              :: frac_oce(nproma,p_patch%alloc_cell_blocks)
    REAL(wp)              :: fwf_fac

    REAL(wp), ALLOCATABLE :: buffer(:,:)

    IF ( .NOT. is_coupled_run() ) RETURN

    ! adjust size if larger bundles are used (no_arr > 4 below)

    ALLOCATE(buffer(nproma*p_patch%nblks_c,4))

    ! As YAC does not touch masked data an explicit initialisation
    ! is required as some compilers are asked to initialise with NaN
    ! and as we loop over the full array.

    buffer(:,:) = 0.0_wp

    jg = p_patch%id

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------

    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_field(jg)% u_stress_tile(:,:,iwtr/iice)  and 
    !    prm_field(jg)% v_stress_tile(:,:,iwtr/iice)  which are the wind stress components over water and ice respectively
    !
    ! 2. prm_field(jg)% evap_tile(:,:,iwtr/iice)  evaporation rate over ice-covered and open ocean/lakes, no land;
    !
    ! 3. prm_field(jg)%rsfl + prm_field(jg)%rsfc + prm_field(jg)%ssfl + prm_field(jg)%ssfc
    !    which gives the precipitation rate;
    !
    ! 4. prm_field(jg)% ta(:,nlev,:)  temperature at the lowest model level, or
    !    prm_field(jg)% tas(:,:)      2-m temperature, not available yet, or
    !    prm_field(jg)% shflx_tile(:,:,iwtr) sensible heat flux
    !    ... tbc
    !
    ! 5  prm_field(jg)% lhflx_tile(:,:,iwtr) latent heat flux
    ! 6. shortwave radiation flux at the surface
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. prm_field(jg)% ts_tile(:,:,iwtr)   SST
    ! 2. prm_field(jg)% ocu(:,:) and ocv(:,:) ocean surface current
    ! 3. ... tbc
    !
    ! 
    nbr_hor_cells = p_patch%n_patch_cells
    !
    !  Send fields to ocean:
    !   field_id(1) represents "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   field_id(3) represents "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   field_id(4) represents "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   field_id(5) represents "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   field_id(10) represents "10m_wind_speed"                          - atmospheric wind speed
    !   field_id(11) represents "qtrc(nlev,co2)"                          - co2 mixing ratio
    !   field_id(13) represents "pres_msl"                                - sea level pressure
    !
    !  Receive fields from ocean:
    !   field_id(6) represents "sea_surface_temperature"                  - SST
    !   field_id(7) represents "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !   field_id(8) represents "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !   field_id(9) represents "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !   field_id(12) represents "co2_flux"                                - ocean co2 flux
    !
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to ocean
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    write_coupler_restart = .FALSE.

    ! Calculate fractional ocean mask 
    ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
    !  - lake part is included in land part, must be subtracted as well
    !  - if no lake part is present, subtract land part only
    !  - if no jsbach is present (aquaplanet), frac_oce is 1.
    IF ( mask_checksum > 0 .AND. echam_phy_config(jg)%ljsb ) THEN
      IF ( echam_phy_config(jg)%llake ) THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          frac_oce(n,i_blk) = 1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd) - prm_field(jg)%alake(n,i_blk)
        ENDDO
      ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

      ELSE
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
        DO i_blk = 1, p_patch%nblks_c
          nn = (i_blk-1)*nproma
          IF (i_blk /= p_patch%nblks_c) THEN
            nlen = nproma
          ELSE
            nlen = p_patch%npromz_c
          END IF
          DO n = 1, nlen
            frac_oce(n,i_blk) = 1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd)
          ENDDO
        ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ENDIF
    ELSE
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          frac_oce(:,:) = 1.0
        ENDDO
      ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
    ENDIF
    !
    ! ------------------------------
    !  Send zonal wind stress bundle
    !   field_id(1) represents "surface_downward_eastward_stress" bundle - zonal wind stress component over ice and water
    !
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
         buffer(nn+n,1) = prm_field(jg)%u_stress_tile(n,i_blk,iwtr)
         buffer(nn+n,2) = prm_field(jg)%u_stress_tile(n,i_blk,iice)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(1), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fput called after end of run - id=1, u-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send meridional wind stress bundle
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
         buffer(nn+n,1) = prm_field(jg)%v_stress_tile(n,i_blk,iwtr)
         buffer(nn+n,2) = prm_field(jg)%v_stress_tile(n,i_blk,iice)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(2), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fput called after end of run - id=2, v-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send surface fresh water flux bundle
    !   field_id(3) represents "surface_fresh_water_flux" bundle - liquid rain, snowfall, evaporation
    !
    !   Note: the evap_tile should be properly updated and added;
    !         as long as evaporation over sea-ice is not used in ocean thermodynamics, the evaporation over the
    !         whole ocean part of grid-cell is passed to the ocean
    !         for pre04 a preliminary solution for evaporation in ocean model is to exclude the land fraction
    !         evap.oce = (evap.wtr*frac.wtr + evap.ice*frac.ice)/(1-frac.lnd)
    !
    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
      scr(:,:)      = 0.0_wp
    !
    ! Preliminary: hard-coded correction factor for freshwater imbalance
    fwf_fac = 1.0_wp      ! neutral factor

    ! Aquaplanet coupling: surface types ocean and ice only
    IF (nsfc_type == 2) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
     
          ! total rates of rain and snow over whole cell
          buffer(nn+n,1) = (prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk))*fwf_fac
          buffer(nn+n,2) = (prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk))*fwf_fac
     
          ! evaporation over ice-free and ice-covered water fraction - of whole ocean part
          buffer(nn+n,3) = prm_field(jg)%evap_tile(n,i_blk,iwtr)*prm_field(jg)%frac_tile(n,i_blk,iwtr) + &
            &              prm_field(jg)%evap_tile(n,i_blk,iice)*prm_field(jg)%frac_tile(n,i_blk,iice)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! Full coupling including jsbach: surface types ocean, ice, land
    ELSE IF (nsfc_type == 3) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
    
          ! total rates of rain and snow over whole cell
          buffer(nn+n,1) = (prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk))*fwf_fac
          buffer(nn+n,2) = (prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk))*fwf_fac
    
          ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
          !  - lake part is included in land part, must be subtracted as well
          !    frac_oce(n,i_blk)= 1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd)-prm_field(jg)%alake(n,i_blk)
          !  - sftof = 1-(land+lake) is already available in prm_field
          !frac_oce(n,i_blk)= prm_field(jg)%sftof(n,i_blk) - not yet tested in coupled model

          IF (frac_oce(n,i_blk) <= 0.0_wp) THEN
            ! land part is zero
            buffer(nn+n,3) = 0.0_wp
          ELSE
            buffer(nn+n,3) = (prm_field(jg)%evap_tile(n,i_blk,iwtr)*prm_field(jg)%frac_tile(n,i_blk,iwtr) + &
              &               prm_field(jg)%evap_tile(n,i_blk,iice)*prm_field(jg)%frac_tile(n,i_blk,iice))/frac_oce(n,i_blk)
          ENDIF
          IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) scr(n,i_blk) = buffer(nn+n,3)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
        &  CALL dbg_print('EchOce: evapo-cpl',scr,str_module,3,in_subset=p_patch%cells%owned)
    ELSE
      CALL finish('interface_echam_ocean: coupling only for nsfc_type equals 2 or 3. Check your code/configuration!')
    ENDIF  !  nsfc_type
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 3
    CALL yac_fput ( field_id(3), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND )                  &
         & CALL warning('interface_echam_ocean', &
         &              'YAC says fput called after end of run - id=3, fresh water flux')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send total heat flux bundle
    !   field_id(4) represents "total heat flux" bundle - short wave, long wave, sensible, latent heat flux
    !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = prm_field(jg)%swflxsfc_tile(n,i_blk,iwtr)
        buffer(nn+n,2) = prm_field(jg)%lwflxsfc_tile(n,i_blk,iwtr)
        buffer(nn+n,3) = prm_field(jg)%shflx_tile   (n,i_blk,iwtr)
        buffer(nn+n,4) = prm_field(jg)%lhflx_tile   (n,i_blk,iwtr)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 4
    CALL yac_fput ( field_id(4), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fput called after end of run - id=4, heat flux')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send sea ice flux bundle
    !   field_id(5) represents "atmosphere_sea_ice_bundle" - sea ice surface and bottom melt potentials Qtop, Qbot
    !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = prm_field(jg)%Qtop(n,1,i_blk)
        buffer(nn+n,2) = prm_field(jg)%Qbot(n,1,i_blk)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(5), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND )                  &
         & CALL warning('interface_echam_ocean', &
         &              'YAC says fput called after end of run - id=5, atmos sea ice')

    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
       CALL message('interface_echam_ocean', 'YAC says it is put for restart - ids 1 to 5, atmosphere fields')
    ENDIF

    !
    ! ------------------------------
    !  Send 10m wind speed
    !   field_id(10) represents "10m_wind_speed" - atmospheric wind speed
    !
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        ! as far as no tiles (pre04) are correctly implemented, use the grid-point mean of 10m wind for coupling
        buffer(nn+n,1) = prm_field(jg)%sfcWind(n,i_blk)
      ENDDO
    ENDDO
!!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 1
    CALL yac_fput ( field_id(10), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == OUT_OF_BOUND )                  &
         & CALL warning('interface_echam_ocean', &
         &              'YAC says fput called after end of run - id=10, wind speed')

    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
       CALL message('interface_echam_ocean', 'YAC says it is put for restart - ids 10, wind speed')
    ENDIF

    !
    ! ------------------------------
    !  Send sea level pressure
    !   field_id(13) represents "pres_msl" - atmospheric sea level pressure
    !
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = pt_diag%pres_msl(n,i_blk)
      ENDDO
    ENDDO
!!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 1
    CALL yac_fput ( field_id(13), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == OUT_OF_BOUND )                  &
         & CALL warning('interface_echam_ocean', &
         &              'YAC says fput called after end of run - id=13, sea level pressure')

    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
       CALL message('interface_echam_ocean', 'YAC says it is put for restart - ids 13, sea level pressure')
    ENDIF

#ifndef __NO_ICON_OCEAN__
    IF (ccycle_config(jg)%iccycle /= 0) THEN
       !
       ! ------------------------------
       !  Send co2 mixing ratio
       !   field_id(11) represents "co2_mixing_ratio" - CO2 mixing ratio in ppmv
       !
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
       DO i_blk = 1, p_patch%nblks_c
          nn = (i_blk-1)*nproma
          IF (i_blk /= p_patch%nblks_c) THEN
             nlen = nproma
          ELSE
             nlen = p_patch%npromz_c
          END IF
          SELECT CASE (ccycle_config(jg)%iccycle)
          CASE (1) ! c-cycle with interactive atm. co2 concentration, qtrc in kg/kg
             DO n = 1, nlen
                buffer(nn+n,1)    =  amd/amco2 * 1.0e6_wp * prm_field(jg)%qtrc(n,nlev,i_blk,ico2)
             END DO
          CASE (2) ! c-cycle with prescribed  atm. co2 concentration
             SELECT CASE (ccycle_config(jg)%ico2conc)
             CASE (2) ! constant  co2 concentration, vmr_co2 in m3/m3
                DO n = 1, nlen
                   buffer(nn+n,1) =              1.0e6_wp * ccycle_config(jg)%vmr_co2
                END DO
             CASE (4) ! transient co2 concentration, ghg_co2mmr in kg/kg
                DO n = 1, nlen
                   buffer(nn+n,1) =  amd/amco2 * 1.0e6_wp * ghg_co2mmr
                END DO
             END SELECT
          END SELECT
       ENDDO
!!ICON_OMP_END_PARALLEL_DO
       !
       IF (ltimer) CALL timer_start(timer_coupling_put)

       no_arr = 1
       CALL yac_fput ( field_id(11), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
       IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
          write_coupler_restart = .TRUE.
       ELSE
          write_coupler_restart = .FALSE.
       ENDIF

       IF ( info == OUT_OF_BOUND )                 &
            & CALL warning('interface_echam_ocean', &
            &              'YAC says fput called after end of run - id=11, co2 mr')

       IF (ltimer) CALL timer_stop(timer_coupling_put)
       !
       IF ( write_coupler_restart ) THEN
          CALL message('interface_echam_ocean', 'YAC says it is put for restart - id=11, co2 mr')
       ENDIF

    ENDIF
#endif



    !
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from ocean to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !  Receive fields, only assign values if something was received ( info > 0 )
    !   - ocean fields have undefined values on land, which are not sent to the atmosphere,
    !     therefore buffer is set to zero to avoid unintended usage of ocean values over land
    !
    buffer(:,:) = 0.0_wp
    !
    ! ------------------------------
    !  Receive SST
    !   field_id(6) represents "sea_surface_temperature" - SST
    !
    IF ( .NOT. lyac_very_1st_get ) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ENDIF

    ! Workaround for lake-points on land that are not ocean:
    !  > no preset (via buffer) to avoid update of ts_wtr over land - land points are untouched
    !  > to avoid errors with SST over lakes (not ocean) the untouched land-points should be set to undef 
    !  > the latter would require a discrimination of temperature between ocean water and lake points

    ! This was wrong:
    ! buffer for tsfc in Kelvin (>0)
    !!! buffer(:,:) = 199.99_wp
    ! buffer set to undefined to enforce error on unintended lake grid-points
    !!! buffer(:,:) = -99.999_wp  ! this aborts with lookup table overflow, since lake-points are affected

    CALL yac_fget ( field_id(6), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('interface_echam_ocean', 'YAC says it is get for restart - id=6, SST')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fget called after end of run - id=6, SST')

    IF ( .NOT. lyac_very_1st_get ) THEN
       IF (ltimer) CALL timer_stop(timer_coupling_1stget)
    ENDIF

    lyac_very_1st_get = .FALSE.

    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen

          !  - lake part is included in land part, must be subtracted as well, see frac_oce

          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ts_tile(n,i_blk,iwtr) = dummy
          ELSE
           !! Workaround for missing discrimination between tile_wtr and tile_lake:
           !!   > background: surface temp. ts_tile(lake)=ts_tile(wtr) is calculated in jsbach and was overwritten
           !!     by default values (buffer_wtr) from ocean
           !!   > ts_tile(wtr) is set over ocean (not lake) points only
           !IF ( ext_data(1)%atm%lsm_ctr_c(n,i_blk) < 0 ) prm_field(jg)%ts_tile(n,i_blk,iwtr) = buffer(nn+n,1)
           !!  for dbg_print only
           !IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN
           !  IF ( ext_data(1)%atm%lsm_ctr_c(n,i_blk) < 0 ) THEN
           !    scr(n,i_blk) = buffer(nn+n,1)
           !  ELSE
           !    scr(n,i_blk) = 285.0_wp  !  value over land - for dbg_print
           !  ENDIF
           !ENDIF
           !
           !!   > slo: using frac_oce>0., constant fractional lsm contains ocean part
           IF ( frac_oce(n,i_blk) > EPSILON(1.0_wp) ) prm_field(jg)%ts_tile(n,i_blk,iwtr) = buffer(nn+n,1)
           !  for dbg_print only
           IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN
             IF ( frac_oce(n,i_blk) > 0.0_wp ) THEN
               scr(n,i_blk) = buffer(nn+n,1)
             ELSE
               scr(n,i_blk) = 285.0_wp  !  value over land - for dbg_print
             ENDIF
           ENDIF
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
        &  CALL dbg_print('EchOce: SSToce-cpl',scr,str_module,4,in_subset=p_patch%cells%owned)
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ts_tile(:,:,iwtr))
    END IF
    !
    ! ------------------------------
    !  Receive zonal velocity
    !   field_id(7) represents "eastward_sea_water_velocity" - zonal velocity, u component of ocean surface current
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    CALL yac_fget ( field_id(7), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('interface_echam_ocean', 'YAC says it is get for restart - id=7, u velocity')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fget called after end of run - id=7, u velocity')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ocu(n,i_blk) = dummy
          ELSE
            prm_field(jg)%ocu(n,i_blk) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
    END IF
    !
    !
    ! ------------------------------
    !  Receive meridional velocity
    !   field_id(8) represents "northward_sea_water_velocity" - meridional velocity, v component of ocean surface current
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    CALL yac_fget ( field_id(8), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('interface_echam_ocean', 'YAC says it is get for restart - id=8, v velocity')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fget called after end of run - id=8, v velocity')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ocv(n,i_blk) = dummy
          ELSE
            prm_field(jg)%ocv(n,i_blk) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
    END IF
    !
    ! ------------------------------
    !  Receive sea ice bundle
    !   field_id(9) represents "ocean_sea_ice_bundle" - ice thickness, snow thickness, ice concentration
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 3
    CALL yac_fget ( field_id(9), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('interface_echam_ocean', 'YAC says it is get for restart - id=9, sea ice')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('interface_echam_ocean', 'YAC says fget called after end of run - id=9, sea ice')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%hi  (n,1,i_blk) = dummy
            prm_field(jg)%hs  (n,1,i_blk) = dummy
            prm_field(jg)%conc(n,1,i_blk) = dummy
          ELSE
            prm_field(jg)%hi  (n,1,i_blk) = buffer(nn+n,1)
            prm_field(jg)%hs  (n,1,i_blk) = buffer(nn+n,2)
            prm_field(jg)%conc(n,1,i_blk) = buffer(nn+n,3)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hi  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hs  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%conc(:,1,:))
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          prm_field(jg)%seaice(n,i_blk) = prm_field(jg)%conc(n,1,i_blk)
          prm_field(jg)%siced(n,i_blk)  = prm_field(jg)%hi(n,1,i_blk)
!  snow thickness not yet connected
!!$          prm_field(jg)%...(n,i_blk)    = prm_field(jg)%hs(n,1,i_blk)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
    END IF
    !
    !
    IF (ccycle_config(jg)%iccycle /= 0) THEN
       !
       ! ------------------------------
       !  Receive co2 flux
       !   field_id(12) represents "co2_flux" - ocean co2 flux
       !
       IF (ltimer) CALL timer_start(timer_coupling_get)

       buffer(:,:) = 0.0_wp ! needs to be checked if this is necessary
       CALL yac_fget ( field_id(12), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
       IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
            & CALL message('interface_echam_ocean', 'YAC says it is get for restart - id=12, CO2 flux')
       IF ( info == OUT_OF_BOUND )                      &
            & CALL warning('interface_echam_ocean', 'YAC says fget called after end of run - id=12, CO2 flux')

       IF (ltimer) CALL timer_stop(timer_coupling_get)
       !
       IF ( info > 0 .AND. info < 7 ) THEN
          !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
          DO i_blk = 1, p_patch%nblks_c
             nn = (i_blk-1)*nproma
             IF (i_blk /= p_patch%nblks_c) THEN
                nlen = nproma
             ELSE
                nlen = p_patch%npromz_c
             END IF
             DO n = 1, nlen
                IF ( nn+n > nbr_inner_cells ) THEN
                   prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = dummy
                ELSE
                   prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = buffer(nn+n,1)
                ENDIF
             ENDDO
          ENDDO
!ICON_OMP_END_PARALLEL_DO
          !
          CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%co2_flux_tile(:,:,iwtr))
       ENDIF
       !
    END IF

!---------DEBUG DIAGNOSTICS-------------------------------------------

    ! calculations for debug print output for namelist debug-values >0 only
    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN

      ! u/v-stress on ice and water sent
      scr(:,:) = prm_field(jg)%u_stress_tile(:,:,iwtr)
      CALL dbg_print('EchOce: u_stress.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%u_stress_tile(:,:,iice)
      CALL dbg_print('EchOce: u_stress.ice',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%v_stress_tile(:,:,iwtr)
      CALL dbg_print('EchOce: v_stress.wtr',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%v_stress_tile(:,:,iice)
      CALL dbg_print('EchOce: v_stress.ice',scr,str_module,4,in_subset=p_patch%cells%owned)

      ! rain, snow, evaporation
      scr(:,:) = prm_field(jg)%rsfl(:,:) + prm_field(jg)%rsfc(:,:)
      CALL dbg_print('EchOce: total rain  ',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%ssfl(:,:) + prm_field(jg)%ssfc(:,:)
      CALL dbg_print('EchOce: total snow  ',scr,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: evaporation ',prm_field(jg)%evap   ,str_module,4,in_subset=p_patch%cells%owned)

      ! total: short wave, long wave, sensible, latent heat flux sent
      scr(:,:) = prm_field(jg)%swflxsfc_tile(:,:,iwtr) + &
        &        prm_field(jg)%lwflxsfc_tile(:,:,iwtr) + &
        &        prm_field(jg)%shflx_tile(:,:,iwtr)    + &
        &        prm_field(jg)%lhflx_tile(:,:,iwtr)
      CALL dbg_print('EchOce: totalhfx.wtr',scr,str_module,2,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%swflxsfc_tile(:,:,iwtr)
      CALL dbg_print('EchOce: swflxsfc.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%lwflxsfc_tile(:,:,iwtr)
      CALL dbg_print('EchOce: lwflxsfc.wtr',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%shflx_tile(:,:,iwtr)
      CALL dbg_print('EchOce: shflx.wtr   ',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%lhflx_tile(:,:,iwtr)
      CALL dbg_print('EchOce: lhflx.wtr   ',scr,str_module,4,in_subset=p_patch%cells%owned)

      ! Qtop and Qbot, windspeed sent
      CALL dbg_print('EchOce: ice-Qtop    ',prm_field(jg)%Qtop        ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: ice-Qbot    ',prm_field(jg)%Qbot        ,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: sfcWind     ',prm_field(jg)%sfcWind     ,str_module,3,in_subset=p_patch%cells%owned)

      ! SST, sea ice, ocean velocity received
      scr(:,:) = prm_field(jg)%ts_tile(:,:,iwtr)
      CALL dbg_print('EchOce: ts_tile.iwtr',scr                       ,str_module,2,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: hi(1)       ',prm_field(jg)%hi(:,1,:)   ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: hs(1)       ',prm_field(jg)%hs(:,1,:)   ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: conc(1)     ',prm_field(jg)%conc(:,1,:) ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: siced       ',prm_field(jg)%siced       ,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: seaice      ',prm_field(jg)%seaice      ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: ocu         ',prm_field(jg)%ocu         ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('EchOce: ocv         ',prm_field(jg)%ocv         ,str_module,4,in_subset=p_patch%cells%owned)

      !error?
      !CALL dbg_print('EchOce: ts_tile.iwtr:iwtr',prm_field(jg)%ts_tile(:,:,iwtr:iwtr),str_module,2,in_subset=p_patch%cells%owned)

      ! Fraction of tiles:
      CALL dbg_print('EchOce: frac_oce     ',frac_oce                 ,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%frac_tile(:,:,iwtr)
      CALL dbg_print('EchOce: frac_tile.wtr',scr                      ,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%frac_tile(:,:,iice)
      CALL dbg_print('EchOce: frac_tile.ice',scr                      ,str_module,3,in_subset=p_patch%cells%owned)
      IF ( echam_phy_config(jg)%ljsb ) THEN
      scr(:,:) = prm_field(jg)%frac_tile(:,:,ilnd)
      CALL dbg_print('EchOce: frac_tile.lnd',scr                      ,str_module,4,in_subset=p_patch%cells%owned)
        IF ( echam_phy_config(jg)%llake ) &
          & CALL dbg_print('EchOce: frac_alake   ',prm_field(jg)%alake,str_module,4,in_subset=p_patch%cells%owned)
      ENDIF
    ENDIF

    !---------------------------------------------------------------------

    DEALLOCATE(buffer)

  END SUBROUTINE interface_echam_ocean

  !>
  !! SUBROUTINE destruct_echam_ocean_coupling -- terminates the coupling
  !! between ECHAM physics and the ocean.
  !!
  !! This subroutine is called at the end of the time loop of the ICONAM model.
  
END MODULE mo_echam_coupling

