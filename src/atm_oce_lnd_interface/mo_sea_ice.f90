!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Dirk Notz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
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
MODULE mo_sea_ice
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: success, max_char_length, sea_boundary 
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs, mu, mus, &
    &                               alf, alv, albedoW, clw, cpd, zemiss_def,rd, stbo,tmelt
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_ocean_nml,           ONLY: no_tracer, init_oce_prog, iforc_oce, &
    &                               FORCING_FROM_FILE_FLUX, i_sea_ice
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base, ocean_var_list
  !USE mo_oce_index,           ONLY: print_mxmn, ipl_src
  USE mo_var_list,            ONLY: add_var
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  ! # achim
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean
  USE mo_sea_ice_winton,      ONLY: ice_growth_winton, set_ice_temp_winton
  USE mo_sea_ice_zerolayer,   ONLY: ice_growth_zerolayer, set_ice_temp_zerolayer
  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx, print_maxmin_si, print_cells

  IMPLICIT NONE

  PRIVATE

  ! Public interface

  ! Definition of forcing types
  ! public types
  ! contained in mo_sea_ice_types

  ! public subroutines
  PUBLIC :: construct_sea_ice 
  PUBLIC :: destruct_sea_ice
  PUBLIC :: construct_sfcflx
  PUBLIC :: construct_atmos_for_ocean
  PUBLIC :: construct_atmos_fluxes
  PUBLIC :: destruct_sfcflx
  PUBLIC :: destruct_atmos_for_ocean
  PUBLIC :: destruct_atmos_fluxes

  PUBLIC :: ice_init
  PUBLIC :: set_ice_albedo
  PUBLIC :: sum_fluxes
  PUBLIC :: ave_fluxes
  PUBLIC :: ice_fast
  PUBLIC :: ice_slow
  PUBLIC :: upper_ocean_TS
  PUBLIC :: new_ice_growth
  PUBLIC :: calc_atm_fluxes_from_bulk
  PUBLIC :: prepareAfterRestart

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !to be put into namelist
!  INTEGER :: i_no_ice_thick_class = 1



CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Constructor of sea-ice model, allocates all components and assigns zero. 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE construct_sea_ice(p_patch, p_ice, i_no_ice_thick_class)
    TYPE(t_patch),     INTENT(IN)    :: p_patch
    TYPE (t_sea_ice),  INTENT(INOUT) :: p_ice
    INTEGER,           INTENT(IN)    :: i_no_ice_thick_class

    !Local variables
    !INTEGER i

    INTEGER :: nblks_c, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_sea_ice'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    nblks_c = p_patch%nblks_c

    p_ice%kice = i_no_ice_thick_class

    ALLOCATE(p_ice%isice(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for isice failed')
    END IF

    CALL add_var(ocean_var_list, 'isice', p_ice%restart_isice ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('isice', '', 'ice mask'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'alb', p_ice%alb ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('alb', '', 'albedo of snow-ice system'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'Tsurf', p_ice%Tsurf ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('Tsurf', '', 'surface temperature'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'T1', p_ice%T1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('T1', 'C', 'Temperature upper layer'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'T2', p_ice%T2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('T2', 'C', 'Temperature lower layer'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'E1', p_ice%E1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('E1', 'Jm/kg', 'Energy content upper layer'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'E2', p_ice%E2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('E2', 'Jm/kg', 'Energy content upper layer'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'hi', p_ice%hi ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('hi', 'm', 'ice thickness'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'hs', p_ice%hs ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('hs', 'm', 'snow thickness'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'hiold', p_ice%hiold ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('hiold', 'm', 'ice thickness (last timstep)'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'hsold', p_ice%hsold ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('hsold', 'm', 'snow thickness (last timstep)'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))

    CALL add_var(ocean_var_list, 'Qtop', p_ice%Qtop ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('Qtop', 'W/m^2', 'Energy flux available for surface melting'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'Qbot', p_ice%Qbot ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('Qbot', 'W/m^2', 'Energy flux at ice-ocean interface'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))

    CALL add_var(ocean_var_list, 'heatocei', p_ice%heatocei ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('heatocei', 'J', 'Energy to ocean when all ice is melted'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'snow_to_ice', p_ice%snow_to_ice ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('snow_to_ice', 'm', 'amount of snow that is transformed to ice'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))

    CALL add_var(ocean_var_list, 'surfmelt', p_ice%surfmelt ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('surfmelt', 'm', 'surface melt water running into ocean'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))
    CALL add_var(ocean_var_list, 'surfmeltT', p_ice%surfmeltT ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('surfmeltT', 'C', 'Mean temperature of surface melt water'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))


    CALL add_var(ocean_var_list, 'evapwi', p_ice%evapwi ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('evapwi', 'kg/m^2', 'amount of evaporated water if no ice left'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))

    CALL add_var(ocean_var_list, 'conc', p_ice%conc ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('conc', '', 'ice concentration in each ice class'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,nblks_c/))

    CALL add_var(ocean_var_list, 'ice_u', p_ice%u ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('ice_u', 'm/s', 'zonal velocity'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list, 'ice_v', p_ice%v ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('ice_v', 'm/s', 'meridional velocity'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,nblks_c/))
      
    CALL add_var(ocean_var_list, 'concSum', p_ice%concSum ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('concSum', '', 'total ice concentration'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,nblks_c/))

    CALL add_var(ocean_var_list, 'newice', p_ice%newice ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, &
      &          t_cf_var('newice', 'm', 'new ice groth in open water'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,nblks_c/))

    CALL add_var(ocean_var_list, 'zUnderIce', p_ice%zUnderIce ,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, &
      &          t_cf_var('zUnderIce', 'm', 'water in upper ocean grid cell below ice'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,nblks_c/))

    ALLOCATE(p_ice%hi_lim(i_no_ice_thick_class), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for hi_lim failed')
    END IF


    IF(p_ice%kice==1)THEN
      p_ice%hi_lim = 0.0_wp
    ELSEIF(p_ice%kice==8)THEN
      p_ice%hi_lim(:)=(/ 0.0_wp, 0.1_wp, 0.3_wp, 0.7_wp, 1.1_wp, 1.5_wp, 2.0_wp, 2.5_wp /)
    ENDIF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_sea_ice
  !-------------------------------------------------------------------------
  !
  !> Destructor of sea-ice model, deallocates all components. 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE destruct_sea_ice(p_ice)
    TYPE (t_sea_ice),  INTENT (INOUT) :: p_ice
    !Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_sea_ice'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(p_ice%hi_lim, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for hi_lim failed')
    END IF

    CALL message(TRIM(routine), 'end' )
   
  END SUBROUTINE destruct_sea_ice
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE construct_sfcflx(p_patch, p_sfc_flx)
    !
    TYPE(t_patch),   INTENT(IN)    :: p_patch
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx

    ! Local variables
    INTEGER :: nblks_c, ist

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_sfcflx'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    nblks_c = p_patch%nblks_c

    ALLOCATE(p_sfc_flx%forc_wind_u(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing wind u failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_wind_v(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing wind v failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_hflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing heat flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_fwfx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing freshwater flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_swflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for short wave flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_lwflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for long wave flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_ssflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for surface sensible heat flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_slflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for surface latent heat flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_prflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for precipitation flux failed')
    END IF
    ALLOCATE(p_sfc_flx%forc_evflx(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for evaporation flux failed')
    END IF
    IF(no_tracer>=1)THEN
      ALLOCATE(p_sfc_flx%forc_tracer(nproma,nblks_c, no_tracer), STAT=ist)
      IF (ist/=SUCCESS) THEN
        CALL finish(TRIM(routine),'allocation for tracer forcing failed')
      END IF

      ALLOCATE(p_sfc_flx%forc_tracer_relax(nproma,nblks_c, no_tracer), STAT=ist)
      IF (ist/=SUCCESS) THEN
        CALL finish(TRIM(routine),'allocation for tracer relaxation forcing failed')
      END IF
    ENDIF

    ALLOCATE(p_sfc_flx%forc_wind_cc(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing wind_cc  failed')
    END IF

    p_sfc_flx%forc_wind_u(:,:)       = 0.0_wp
    p_sfc_flx%forc_wind_v(:,:)       = 0.0_wp
    
    ! init of cartesian coordinates:
    p_sfc_flx%forc_wind_cc(:,:)%x(1) = 0.0_wp
    p_sfc_flx%forc_wind_cc(:,:)%x(2) = 0.0_wp
    p_sfc_flx%forc_wind_cc(:,:)%x(3) = 0.0_wp

    IF(no_tracer>=1)THEN
      p_sfc_flx%forc_hflx        (:,:)    = 0.0_wp
      p_sfc_flx%forc_fwfx        (:,:)    = 0.0_wp
      p_sfc_flx%forc_swflx       (:,:)    = 0.0_wp
      p_sfc_flx%forc_lwflx       (:,:)    = 0.0_wp
      p_sfc_flx%forc_ssflx       (:,:)    = 0.0_wp
      p_sfc_flx%forc_slflx       (:,:)    = 0.0_wp
      p_sfc_flx%forc_prflx       (:,:)    = 0.0_wp
      p_sfc_flx%forc_evflx       (:,:)    = 0.0_wp
      p_sfc_flx%forc_tracer      (:,:,:)  = 0.0_wp
      p_sfc_flx%forc_tracer_relax(:,:,:)  = 0.0_wp
    ENDIF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_sfcflx
  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor surface flux forcing for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE destruct_sfcflx(p_sfc_flx)
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx
    !
    ! Local variables

    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_sfcflx'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(p_sfc_flx%forc_wind_u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing wind u failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_wind_v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing wind v failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_hflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing heat flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_fwfx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing freshwater flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_swflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for heat flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_lwflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for heat flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_ssflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for heat flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_slflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for heat flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_prflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for precip flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_evflx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for evap flux failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_tracer, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for tracer forcing failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_tracer_relax, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for tracer relaxation failed')
    END IF
    DEALLOCATE(p_sfc_flx%forc_wind_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing wind cc failed')
    END IF
    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_sfcflx
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmospheric reprsentation  in ocean.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07)
  !
  SUBROUTINE construct_atmos_for_ocean(p_patch, p_as)
    !
    TYPE(t_patch),                INTENT(IN):: p_patch
    TYPE(t_atmos_for_ocean ), INTENT(INOUT) :: p_as

    ! Local variables
    INTEGER :: nblks_c, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_for_ocean'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    nblks_c = p_patch%nblks_c
   
    ALLOCATE(p_as%tafo(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for tafo failed')
    END IF
    ALLOCATE(p_as%ftdew(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ftdew failed')
    END IF
    ALLOCATE(p_as%fclou(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fclou failed')
    END IF

    ALLOCATE(p_as%fu10(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fu10 failed')
    END IF

    ALLOCATE(p_as%fswr(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fswr failed')
    END IF

    ALLOCATE(p_as%pao(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for pao failed')
    END IF

    ALLOCATE(p_as%u(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for u failed')
    END IF
    ALLOCATE(p_as%v(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for v failed')
    END IF


    p_as%tafo (:,:) = 0.0_wp
    p_as%ftdew(:,:) = 0.0_wp
    p_as%fclou(:,:) = 0.0_wp
    p_as%fu10 (:,:) = 0.0_wp
    p_as%fswr (:,:) = 0.0_wp
    p_as%pao  (:,:) = 0.0_wp
    p_as%u    (:,:) = 0.0_wp
    p_as%v    (:,:) = 0.0_wp

    CALL message(TRIM(routine), 'end')

  END SUBROUTINE construct_atmos_for_ocean
  !-------------------------------------------------------------------------
  !
  !>
  !!  Destructor of atmospheric reprsentation  in ocean.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE destruct_atmos_for_ocean(p_as)
    !
    TYPE(t_atmos_for_ocean ), INTENT(INOUT) :: p_as

    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_for_ocean'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

   
    DEALLOCATE(p_as%tafo, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for tafo failed')
    END IF
    DEALLOCATE(p_as%ftdew, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ftdew failed')
    END IF
    DEALLOCATE(p_as%fclou, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fclou failed')
    END IF

    DEALLOCATE(p_as%fu10, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fu10 failed')
    END IF

    DEALLOCATE(p_as%fswr, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fswr failed')
    END IF

    DEALLOCATE(p_as%pao, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for pao failed')
    END IF

    DEALLOCATE(p_as%u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for u failed')
    END IF
    DEALLOCATE(p_as%v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for v failed')
    END IF

    CALL message(TRIM(routine), 'end')

  END SUBROUTINE destruct_atmos_for_ocean
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmos fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE construct_atmos_fluxes(p_patch, p_atm_f, i_no_ice_thick_class)
    !
    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_atmos_fluxes ), INTENT(INOUT) :: p_atm_f
    INTEGER,               INTENT(IN)    :: i_no_ice_thick_class
    ! Local variables
    INTEGER :: nblks_c, ist

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_fluxes'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    nblks_c = p_patch%nblks_c
   
    ALLOCATE(p_atm_f%sens(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sens failed')
    END IF

    ALLOCATE(p_atm_f%lat(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for lat failed')
    END IF

    ALLOCATE(p_atm_f%LWout(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWout failed')
    END IF

    ALLOCATE(p_atm_f%LWnet(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWnet failed')
    END IF

    ALLOCATE(p_atm_f%bot(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sens failed')
    END IF

    ALLOCATE(p_atm_f%dsensdT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for dsensdT failed')
    END IF

    ALLOCATE(p_atm_f%dlatdT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sens failed')
    END IF

    ALLOCATE(p_atm_f%dLWdT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for dLWdT failed')
    END IF

    ALLOCATE(p_atm_f%rprecw(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for rprecw failed')
    END IF

    ALLOCATE(p_atm_f%rpreci(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for rpreci failed')
    END IF

    ALLOCATE(p_atm_f%sensw(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sensw failed')
    END IF

    ALLOCATE(p_atm_f%latw(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for latw failed')
    END IF


    ALLOCATE(p_atm_f%LWoutw(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWoutw failed')
    END IF

    ALLOCATE(p_atm_f%LWnetw(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWnetw failed')
    END IF

    ALLOCATE(p_atm_f%SWin(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for SWin failed')
    END IF

    ALLOCATE(p_atm_f%LWin(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWin failed')
     END IF

    ! Initialise everything with zero
    p_atm_f%sens   (:,:,:) = 0.0_wp
    p_atm_f%lat    (:,:,:) = 0.0_wp
    p_atm_f%LWout  (:,:,:) = 0.0_wp
    p_atm_f%LWnet  (:,:,:) = 0.0_wp
    p_atm_f%bot    (:,:,:) = 0.0_wp
    p_atm_f%dsensdT(:,:,:) = 0.0_wp
    p_atm_f%dlatdT (:,:,:) = 0.0_wp
    p_atm_f%dLWdT  (:,:,:) = 0.0_wp
    p_atm_f%rprecw (:,:)   = 0.0_wp
    p_atm_f%rpreci (:,:)   = 0.0_wp
    p_atm_f%sensw  (:,:)   = 0.0_wp
    p_atm_f%latw   (:,:)   = 0.0_wp
    p_atm_f%LWoutw (:,:)   = 0.0_wp
    p_atm_f%LWnetw (:,:)   = 0.0_wp
    p_atm_f%SWin   (:,:)   = 0.0_wp
    p_atm_f%LWin   (:,:)   = 0.0_wp
    p_atm_f%counter        = 0

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_atmos_fluxes
  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor of atmos fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE destruct_atmos_fluxes(p_atm_f)
    !
    TYPE(t_atmos_fluxes )       :: p_atm_f
    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_fluxes'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

   
    DEALLOCATE(p_atm_f%sens, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sens failed')
    END IF

    DEALLOCATE(p_atm_f%lat, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for lat failed')
    END IF

    DEALLOCATE(p_atm_f%LWout, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWout failed')
    END IF

    DEALLOCATE(p_atm_f%LWnet, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWnet failed')
    END IF

    DEALLOCATE(p_atm_f%bot, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sens failed')
    END IF

    DEALLOCATE(p_atm_f%dsensdT, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for dsensdT failed')
    END IF

    DEALLOCATE(p_atm_f%dlatdT, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sens failed')
    END IF

    DEALLOCATE(p_atm_f%dLWdT, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for dLWdT failed')
    END IF

    DEALLOCATE(p_atm_f%rprecw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for rprecw failed')
    END IF

    DEALLOCATE(p_atm_f%rpreci, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for rpreci failed')
    END IF

    DEALLOCATE(p_atm_f%sensw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sensw failed')
    END IF

    DEALLOCATE(p_atm_f%latw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for latw failed')
    END IF


    DEALLOCATE(p_atm_f%LWoutw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWoutw failed')
    END IF

    DEALLOCATE(p_atm_f%LWnetw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWnetw failed')
    END IF

    DEALLOCATE(p_atm_f%SWin, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for SWin failed')
    END IF

    DEALLOCATE(p_atm_f%LWin, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWin failed')
    END IF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_atmos_fluxes

  !-------------------------------------------------------------------------
  !
  !> ice_init
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_init( p_patch, p_os, ice) !, Qatm, QatmAve)
    TYPE(t_patch), INTENT(in)             :: p_patch 
    TYPE(t_hydro_ocean_state)             :: p_os
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

    !local variables
    REAL(wp), DIMENSION(nproma,ice%kice, p_patch%nblks_c) :: &
      & Tinterface, & ! temperature at snow-ice interface
      & draft,      & ! position of ice-ocean interface below sea level
      & Tfw           ! Ocean freezing temperature [°C]

    !INTEGER i,j,k      ! counter for loops
    INTEGER k      ! counter for loops
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:ice_init'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    !Constructor basic init already done at this point
    !   CALL alloc_mem_commo_ice (ice, Qatm, QatmAve)
    !   CALL ice_zero            (ice, Qatm, QatmAve)

    ! FORALL(i=1:nproma, j=1:p_patch%nblks_c, k=1:ice%kice) 
    !    ice% hi    (i,j,k) = sictho (i,j)
    !    ice% hs    (i,j,k) = sicsno (i,j)
    ! END FORALL

    IF ( no_tracer >= 2 ) THEN
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ELSE
      Tfw(:,:,:) = Tf
    ENDIF
      
    ice% Tsurf(:,:,:)  = Tf
    ice% T1   (:,:,:)  = Tf
    ice% T2   (:,:,:)  = Tf
    ice% conc (:,:,:)  = 0.0_wp
    ice% isice(:,:,:)  = .FALSE.
    draft     (:,:,:)  = 0.0_wp

    ! Stupid initialisation trick for Levitus initialisation
    IF (init_oce_prog == 1) THEN
      WHERE (p_os%p_prog(nold(1))%tracer(:,1,:,1) <= -1.0_wp &
          &     .and. v_base%lsm_oce_c(:,1,:) <= sea_boundary )
        ice%hi(:,1,:) = 2._wp
        ice%conc(:,1,:) = 1._wp
      ENDWHERE
      IF ( no_tracer < 2 ) THEN
        WHERE (p_os%p_prog(nold(1))%tracer(:,:,:,1) <= -1.0_wp    &
          &     .and. v_base%lsm_oce_c(:,:,:) <= sea_boundary )   &
          &             p_os%p_prog(nold(1))%tracer(:,:,:,1) = Tf
      ENDIF
    ENDIF

    WHERE(ice% hi(:,:,:) > 0.0_wp)
      ice% Tsurf (:,:,:) = Tfw(:,:,:)
      ice% T1    (:,:,:) = Tfw(:,:,:)
      ice% T2    (:,:,:) = Tfw(:,:,:)
      Tinterface (:,:,:) = (Tfw(:,:,:) * (ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))+&
        &                    ice%Tsurf(:,:,:)) / (1.0_wp+ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))
      ice% conc  (:,:,:) = 1.0_wp/REAL(ice%kice,wp)
      ice% isice (:,:,:) = .TRUE.
      ice% T1    (:,:,:) = Tfw(:,:,:) + 2._wp/3._wp*(Tinterface(:,:,:)-Tfw(:,:,:))
      ice% T2    (:,:,:) = Tfw(:,:,:) + 1._wp/3._wp*(Tinterface(:,:,:)-Tfw(:,:,:))
      draft      (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
    END WHERE
    
    !ice%zUnderIce (:,:)   = dzw(1) + zo (:,:) &
    !  &                     - sum(draft(:,:,:) * ice%conc(:,:,:),2)
    ice%zUnderIce (:,:)   = v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) &
      &                      - sum(draft(:,:,:) * ice%conc(:,:,:),2)

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE ice_init
  !-------------------------------------------------------------------------  
  !
  !  
  !>
  !! !  ice_fast: Ice routines for atmospheric time step. Sets air-ice fluxes and
  !!    calculates the development of the ice temperature field
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_fast(p_patch, ice,Tfw,Qatm,QatmAve)

    TYPE(t_patch),            INTENT(IN)     :: p_patch 
    !TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    !TYPE(t_atmos_for_ocean),  INTENT(IN)     :: p_as
    REAL(wp),                 INTENT(IN)     :: Tfw(:,:,:)
    TYPE (t_sea_ice),         INTENT (INOUT) :: ice
    TYPE (t_atmos_fluxes),    INTENT (IN)    :: Qatm
    TYPE (t_atmos_fluxes),    INTENT (INOUT) :: QatmAve

    !------------------------------------------------------------------------- 
    CALL prepareAfterRestart(ice)

    !CALL get_atmos_fluxes (p_patch, p_os,p_as,ice, Qatm)
    CALL set_ice_albedo(p_patch,ice)
    
    ! #achim
    IF      ( i_sea_ice == 1 ) THEN
      CALL set_ice_temp_winton  (p_patch,ice, Tfw, Qatm)
    ELSE IF ( i_sea_ice == 2 ) THEN
      CALL set_ice_temp_zerolayer  (p_patch,ice, Tfw, Qatm)
    END IF

    CALL sum_fluxes    (Qatm, QatmAve)

   END SUBROUTINE ice_fast
  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_slow: Ice routines for ocean time step. Calculates average of atmospheric
  ! !           time steps, ice velocity, ice growth rates and updates ice structure
  ! !           accordingly
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_slow(p_patch, p_os,ice, QatmAve, p_sfc_flx)  
    TYPE(t_patch),            INTENT(IN)     :: p_patch 
    TYPE(t_hydro_ocean_state),INTENT(INOUT)  :: p_os
    !TYPE(t_atmos_for_ocean),  INTENT(IN)     :: p_as
    TYPE (t_sea_ice),         INTENT (INOUT) :: ice
    !TYPE (t_atmos_fluxes),    INTENT (INOUT) :: Qatm
    TYPE (t_atmos_fluxes),    INTENT (INOUT) :: QatmAve
    TYPE(t_sfc_flx),          INTENT (INOUT) :: p_sfc_flx
  !-------------------------------------------------------------------------------

    CALL ave_fluxes     (ice, QatmAve)
    !CALL ice_dynamics   (ice, QatmAve)
    
    ! #achim
    IF      ( i_sea_ice == 1 ) THEN
      CALL ice_growth_winton    (p_patch,p_os,ice, QatmAve%rpreci)!, QatmAve%lat)
    ELSE IF ( i_sea_ice == 2 ) THEN
      CALL ice_growth_zerolayer (p_patch,p_os,ice, QatmAve%rpreci)
    END IF

    CALL upper_ocean_TS (p_patch,p_os,ice, QatmAve, p_sfc_flx)
    CALL new_ice_growth (p_patch,ice, p_os,p_sfc_flx)
    !CALL ice_advection  (ice)
    !CALL write_ice      (ice,QatmAve,1,ie,je)
    CALL ice_zero       (ice, QatmAve)
    !sictho = ice%hi   (:,:,1) * ice%conc (:,:,1)
    !sicomo = ice%conc (:,:,1)
    !sicsno = ice%hs   (:,:,1) * ice%conc (:,:,1)

    CALL prepare4restart(ice)
  END SUBROUTINE ice_slow
  !-------------------------------------------------------------------------  
  !
  !  
  !>
  !! !  get_atmos_fluxes: Sets the atmospheric fluxes for the update of the ice
  ! !                 temperature
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE get_atmos_fluxes (p_patch, p_os,p_as,ice, Qatm)
    TYPE(t_patch),            INTENT(IN)    :: p_patch 
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE (t_sea_ice),         INTENT(INOUT) :: ice
    TYPE (t_atmos_fluxes),    INTENT(INOUT) :: Qatm

!#ifdef coupled
    !Qatm% SWin   = 
    !Qatm% LWin   =
    !Qatm% sens   = 
    !Qatm% lat    =
    !Qatm% dsensdT = 
    !Qatm% dlatdT  =
    !Qatm% dLWdT   =
!#elif defined CORE
    !CALL budget_core   (ice, Qatm)
!#else
    CALL calc_atm_fluxes_from_bulk(p_patch, p_as, p_os, ice, Qatm)
!#endif

  END SUBROUTINE get_atmos_fluxes 
  !-------------------------------------------------------------------------  
  !
  !  
  !>
  !! !   sum_fluxes: adds atmospheric fluxes for ocean time stepping. Necessary for
  !!      diagnosis, not for the ice model itself.
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE sum_fluxes        (Qatm, QatmAve)
    TYPE (t_atmos_fluxes), INTENT (IN)    :: Qatm
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

    QatmAve % sens   (:,:,:) = QatmAve % sens   (:,:,:) + Qatm % sens   (:,:,:)
    QatmAve % sensw  (:,:)   = QatmAve % sensw  (:,:)   + Qatm % sensw  (:,:)  
    QatmAve % lat    (:,:,:) = QatmAve % lat    (:,:,:) + Qatm % lat    (:,:,:)
    QatmAve % latw   (:,:)   = QatmAve % latw   (:,:)   + Qatm % latw   (:,:)  
    QatmAve % LWout  (:,:,:) = QatmAve % LWout  (:,:,:) + Qatm % LWout  (:,:,:)
    QatmAve % LWoutw (:,:)   = QatmAve % LWoutw (:,:)   + Qatm % LWoutw (:,:)  
    QatmAve % LWnet  (:,:,:) = QatmAve % LWnet  (:,:,:) + Qatm % LWnet  (:,:,:)
    QatmAve % LWnetw (:,:)   = QatmAve % LWnetw (:,:)   + Qatm % LWnetw (:,:)  
    QatmAve % SWin   (:,:)   = QatmAve % SWin   (:,:)   + Qatm % SWin   (:,:)  
    QatmAve % LWin   (:,:)   = QatmAve % LWin   (:,:)   + Qatm % LWin   (:,:)  
    QatmAve % rprecw (:,:)   = QatmAve % rprecw (:,:)   + Qatm % rprecw (:,:)  
    QatmAve % rpreci (:,:)   = QatmAve % rpreci (:,:)   + Qatm % rpreci (:,:)  
    QatmAve % counter        = QatmAve % counter + 1 

  END SUBROUTINE sum_fluxes  
  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! ave_fluxes: calculates the average of the atmospheric fluxes for ocean time  
  !!   sum_fluxes: adds atmospheric fluxes for ocean time stepping. Necessary for
  !!   diagnosis, not for the ice model itself.
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ave_fluxes (ice, QatmAve)
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve
    !
    !Local variables
    REAL(wp) :: ctr

    !-------------------------------------------------------------------------------

    ctr = REAL(QatmAve% counter,wp)
    QatmAve% sens   (:,:,:) = QatmAve% sens  (:,:,:)  / ctr
    QatmAve% sensw  (:,:)   = QatmAve% sensw (:,:)    / ctr
    QatmAve% lat    (:,:,:) = QatmAve% lat   (:,:,:)  / ctr
    QatmAve% latw   (:,:)   = QatmAve% latw  (:,:)    / ctr
    QatmAve% LWout  (:,:,:) = QatmAve% LWout (:,:,:)  / ctr
    QatmAve% LWoutw (:,:)   = QatmAve% LWoutw(:,:)    / ctr
    QatmAve% LWnet  (:,:,:) = QatmAve% LWnet (:,:,:)  / ctr
    QatmAve% LWnetw (:,:)   = QatmAve% LWnetw(:,:)    / ctr
    QatmAve% SWin   (:,:)   = QatmAve% SWin  (:,:)    / ctr
    QatmAve% LWin   (:,:)   = QatmAve% LWin  (:,:)    / ctr
    QatmAve% rprecw (:,:)   = QatmAve% rprecw(:,:)    / ctr
    QatmAve% rpreci (:,:)   = QatmAve% rpreci(:,:)    / ctr
    ice    % Qbot   (:,:,:) = ice    % Qbot  (:,:,:)  / ctr
    ice    % Qtop   (:,:,:) = ice    % Qtop  (:,:,:)  / ctr

  END SUBROUTINE ave_fluxes
  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! ice_zero: set the avereged fluxes to zero
  !! @par Revision History
  !! Initial release by Einar Olason, MPI-M (2011-09). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_zero (ice,QatmAve)
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

    !Qatm    % sens        (:,:,:) = 0._wp
    !Qatm    % sensw       (:,:)   = 0._wp
    !Qatm    % lat         (:,:,:) = 0._wp
    !Qatm    % latw        (:,:)   = 0._wp
    !Qatm    % LWout       (:,:,:) = 0._wp
    !Qatm    % LWoutw      (:,:)   = 0._wp
    !Qatm    % LWnet       (:,:,:) = 0._wp
    !Qatm    % LWnetw      (:,:)   = 0._wp
    !Qatm    % SWin        (:,:)   = 0._wp
    !Qatm    % LWin        (:,:)   = 0._wp
    !Qatm    % rprecw      (:,:)   = 0._wp
    !Qatm    % rpreci      (:,:)   = 0._wp
                          
    QatmAve % sens        (:,:,:) = 0._wp
    QatmAve % sensw       (:,:)   = 0._wp
    QatmAve % lat         (:,:,:) = 0._wp
    QatmAve % latw        (:,:)   = 0._wp
    QatmAve % LWout       (:,:,:) = 0._wp
    QatmAve % LWoutw      (:,:)   = 0._wp
    QatmAve % LWnet       (:,:,:) = 0._wp
    QatmAve % LWnetw      (:,:)   = 0._wp
    QatmAve % SWin        (:,:)   = 0._wp
    QatmAve % LWin        (:,:)   = 0._wp
    QatmAve % rprecw      (:,:)   = 0._wp
    QatmAve % rpreci      (:,:)   = 0._wp
    QatmAve % counter             = 0 

    ice     % Qbot        (:,:,:) = 0._wp
    ice     % Qtop        (:,:,:) = 0._wp
    ice     % surfmelt    (:,:,:) = 0._wp
    ice     % surfmeltT   (:,:,:) = 0._wp
    ice     % evapwi      (:,:,:) = 0._wp
    ice     % hiold       (:,:,:) = 0._wp
    ice     % snow_to_ice (:,:,:) = 0._wp
    ice     % heatOceI    (:,:,:) = 0._wp

  END SUBROUTINE ice_zero

  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! ice_albedo: set ice albedo 
  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! ice_albedo: set ice albedo 
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE set_ice_albedo(p_patch, ice) 
    TYPE(t_patch),    INTENT(IN)    :: p_patch 
    TYPE (t_sea_ice), INTENT(INOUT) :: ice
    !
    !Local variables
    REAL(wp), PARAMETER :: albtrans   = 0.5_wp
    REAL(wp)            :: albflag(nproma,ice%kice, p_patch%nblks_c)
    !-------------------------------------------------------------------------------

    ! This is Uwe's albedo expression from the old budget function
    albflag (:,:,:) =  1.0_wp/ ( 1.0_wp+albtrans * (ice%tsurf(:,:,:))**2 )
    
    WHERE (ice  % isice(:,:,:))
      WHERE (ice % hs(:,:,:) > 1.e-2_wp)
        ice% alb(:,:,:) =  albflag(:,:,:) * albsm + (1.0_wp-albflag(:,:,:)) * albs
      ELSEWHERE
        ice% alb(:,:,:) =  albflag(:,:,:) * albim + (1.0_wp-albflag(:,:,:)) * albi
      END WHERE
    END WHERE

  END SUBROUTINE set_ice_albedo
  

  
  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! upper_ocean_TS: Adjusts the temperature and salinity of the upper ocean grid
  !!                 cell according to atmospheric heat and fresh-water fluxes,
  !!                 surface melting, ice growth, etc. The upper ocean temperature
  !!                 is also changed in subroutine new_ice_growth and at the
  !!                beginning of subroutine ice_growth
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE upper_ocean_TS(p_patch, p_os,ice, QatmAve, p_sfc_flx)
    TYPE(t_patch),             INTENT(IN)    :: p_patch 
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os
    !TYPE(t_atmos_for_ocean),   INTENT(IN)    :: p_as
    TYPE(t_sea_ice),           INTENT(INOUT) :: ice
    TYPE(t_atmos_fluxes),      INTENT(INOUT) :: QatmAve
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx

    !Local Variables
    ! position of ice-ocean interface below sea level                       [m] 
    REAL(wp) :: draft(nproma,ice%kice, p_patch%nblks_c)
    
    REAL(wp), DIMENSION (nproma, p_patch%nblks_c) ::   & 
      & draftAve,      &! average draft of sea ice within a grid cell             [m]
      & zUnderIceOld,  &! water in upper ocean grid cell below ice (prev. time)   [m]
      & heatOceI,      &! heat flux into ocean through formerly ice covered areas [W/m^2]
      & heatOceW,      &! heat flux into ocean through open water areas           [W/m^2]
      & delHice,       &! average change in ice thickness within a grid cell      [m]
      & snowiceave,    &! average snow to ice conversion within a grid cell       [m]
      & evap,          &! evaporated water                                        [m]
      & preci,         &! solid precipitation                                     [m]
      & precw           ! liquid precipitation                                    [m]

    ! Needs work with FB_BGC_OCE etc.
    !REAL(wp)         :: swsum 
    !REAL(wp),POINTER :: sao_top(:,:)
    !-------------------------------------------------------------------------------

    ! #eoo# What is swsum?
    ! swsum = 0.0_wp
    !sao_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

    ! Ocean points only
    ! Calculate change in water level 'zo' from liquid and solid precipitation and
    ! evaporation
    precw           (:,:)   = QatmAve% rprecw (:,:) * dtime
    preci           (:,:)   = QatmAve% rpreci (:,:) * dtime
    evap            (:,:)   = (QatmAve% latw(:,:)/ alv * dtime * &
      &                       sum(ice%conc(:,:,:), 2) +          &
      &                       sum(ice%evapwi(:,:,:) * ice% conc(:,:,:), 2)) /rho_ref

    ! TODO: This should probably be done via surface fluxes?
    !p_os%p_prog(nold(1))%h(:,:) = p_os%p_prog(nold(1))%h(:,:) +  precw + preci - evap

    ! Calculate average draft and thickness of water underneath ice in upper ocean
    ! grid box
    zUnderIceOld    (:,:)   = ice%zUnderIce(:,:)
    draft           (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
    draftave        (:,:)   = sum(draft(:,:,:) * ice%conc(:,:,:),2)
    ice%zUnderIce   (:,:)   = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(:,:) - draftave(:,:) 
   
    ! Calculate average change in ice thickness and the snow-to-ice conversion 
    Delhice         (:,:)   = sum((ice% hi(:,:,:) - ice% hiold(:,:,:))*          &
      &                       ice%conc(:,:,:),2)
    snowiceave      (:,:)   = sum(ice%snow_to_ice(:,:,:) * ice% conc(:,:,:),2)
   

    ! Calculate heat input through formerly ice covered and through open water areas
    !heatOceW        (:,:)   = (QatmAve%SWin(:,:) * (1.0_wp-albedoW) * (1.0_wp-swsum) +    &
    ! A temporary hack: For iforc_oce == 12 (FORCING_FROM_FILE_FLUX) we have OMIP data (or similar)
    ! and need to apply oceanic albedo to the short-wave flux.  For other cases we assume that the
    ! albedo has already been applied.
    !!! We need a unified albedo calculation !!!
    heatOceI    (:,:)   = sum(ice% heatOceI(:,:,:) * ice% conc(:,:,:),2)
    IF ( iforc_oce == FORCING_FROM_FILE_FLUX ) THEN
      heatOceW  (:,:) = QatmAve%SWin(:,:) * (1.0_wp-albedoW)
    ELSE
      heatOceW  (:,:) = QatmAve%SWin(:,:)
    ENDIF
    heatOceW    (:,:) = ( heatOceW(:,:) + QatmAve%LWnetw(:,:) + QatmAve%sensw(:,:)+ &
      &                 QatmAve%latw(:,:) ) *  (1.0_wp-sum(ice%conc(:,:,:),2))

    ! Change temperature of upper ocean grid cell according to heat fluxes
    !p_os%p_prog(nold(1))%tracer(:,1,:,1) = p_os%p_prog(nold(1))%tracer(:,1,:,1)&
    !  &                                    + dtime*(heatOceI + heatOceW) /               &
    !  &                                    (clw*rho_ref * ice%zUnderIce)
    ! TODO: should we also divide with ice%zUnderIce / ( v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) ) ?
    !p_sfc_flx%forc_tracer(:,:,1) = (heatOceI + heatOceW) / (clw*rho_ref)
    p_sfc_flx%forc_hflx(:,:) = heatOceI(:,:) + heatOceW(:,:)

    ! TODO:
    ! Temperature change of upper ocean grid cell due  to melt-water inflow and
    ! precipitation
    !p_os%p_prog(nold(1))%tracer(:,1,:,1) = (p_os%p_prog(nold(1))%tracer(:,1,:,1) &
    !  &                      *zUnderIceOld                                       &
    !  &                      + precw*p_as%tafo + preci*0.0_wp + &                             !!!!!!!!!Dirk: times 0.0 ????
    !  &                        sum(ice%surfmeltT(:,:,:) * ice%surfmelt * ice%conc(:,:,:),2)) / & 
    !  &                        (zUnderIceOld + sum(ice%surfmelt*ice%conc(:,:,:),2) +    &
    !  &                        precw + preci)
    !
    ! Change salinity of upper ocean grid box from ice growth/melt, snowice
    ! formation and precipitation
    !p_os%p_prog(nold(1))%tracer(:,1,:,2) = p_os%p_prog(nold(1))%tracer(:,1,:,2)  &
    !  &                                    + (Delhice(:,:)*rhoi - snowiceave(:,:)*rhos)/rho_ref *  &
    !  &                                    MIN(Sice, sao_top(:,:)) / ice%zUnderIce(:,:)

    !heatabs         (:,:)   = swsum * QatmAve% SWin(:,:) * (1 - ice%concsum)
    !    CALL print_mxmn('SST',1,p_os%p_prog(nold(1))%tracer(:,1,:,1),1,p_patch%nblks_c,'ice',ipl_src)

  END SUBROUTINE upper_ocean_TS
  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! !! new_ice_growth: Calculates the grid-cell average thickness of new ice 
  !                 forming in open-water areas
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  ! TODO: This needs to be rewritten to take in to account cases where the ice concentration can vary
  ! between 0 and 1
  SUBROUTINE new_ice_growth(p_patch,ice, p_os,p_sfc_flx)
    TYPE(t_patch),             INTENT(IN)    :: p_patch 
    TYPE (t_sea_ice),          INTENT(INOUT) :: ice  
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os
    !TYPE (t_atmos_fluxes),     INTENT(IN)    :: QatmAve
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx

    REAL(wp) :: sst(nproma,p_patch%nblks_c)
    REAL(wp) :: Tfw(nproma,p_patch%nblks_c) ! Ocean freezing temperature [°C]

    if ( no_tracer >= 2 ) then
      Tfw(:,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
    else
      Tfw(:,:) = Tf
    endif
    
    ! Calculate possible super-cooling of the surface layer
    sst = p_os%p_prog(nold(1))%tracer(:,1,:,1) + &
      &      dtime*p_sfc_flx%forc_tracer(:,:,1)/ice%zUnderIce(:,:)

    ice % newice(:,:) = 0.0_wp
    WHERE (sst < Tfw(:,:) .and. v_base%lsm_oce_c(:,1,:) <= sea_boundary )
      ice%newice(:,:) = - (sst - Tfw(:,:)) * ice%zUnderIce(:,:) * clw*rho_ref / (alf*rhoi)
      ! Add energy for new-ice formation due to supercooled ocean to  ocean temperature
      p_sfc_flx%forc_tracer(:,:,1) = &
        &     ice%zUnderIce(:,:) * ( Tfw(:,:) - p_os%p_prog(nold(1))%tracer(:,1,:,1) ) / dtime
    END WHERE

    WHERE(ice%newice(:,:)>0.0_wp)
      WHERE(.NOT.ice%isice(:,1,:))
        ice%Tsurf(:,1,:) = Tfw(:,:)
        ice%T2   (:,1,:) = Tfw(:,:)
        ice%T1   (:,1,:) = Tfw(:,:)
      ENDWHERE
      ice % isice(:,1,:) = .TRUE.
      ice % hi   (:,1,:) = ice%newice(:,:)* (1.0_wp-sum(ice%conc(:,:,:),2))&
                         &+ice%hi(:,1,:)*sum(ice%conc(:,:,:),2)
      !ice % hs   (:,:,1) = 0
      !ice % Tsurf(:,1,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)
  !!!!!!!!!!!DIRK: Where is rhs coming from ???????????????

      !ice % T1   (:,:,1) = T1(:,:,1)
      !ice % T2   (:,:,1) = T2(:,:,1)
      ice % conc (:,1,:) = 1.0_wp
    ENDWHERE
    ice% concSum(:,:)  = SUM(ice% conc(:,:,:),2)

  END SUBROUTINE new_ice_growth


  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of the ice 
  !! temperature and ice growth rates for OMIP forcing
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !
  SUBROUTINE calc_atm_fluxes_from_bulk(p_patch, p_as, p_os, p_ice, Qatm)
    TYPE(t_patch),            INTENT(IN)    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: Qatm


    !Local variables
    REAL(wp), DIMENSION (nproma,p_patch%nblks_c) ::           &
      & Tsurf,          &  ! Surface temperature                             [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & esti,           &  ! water vapor pressure at ice surface             [Pa]
      & estw,           &  ! water vapor pressure at water surface           [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height 
      & sphumidi,       &  ! Specific humididty at ice surface
      & sphumidw,       &  ! Specific humididty at water surface
      & ftdewC,         &  ! Dew point temperature in Celsius                [C]
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl                                   
      & dragl1,         &  ! part of dragl                                   
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fw, fi,     &  ! Enhancment factor for vapor pressure
      & dsphumididesti, &  ! Derivative of sphumidi w.r.t. esti
      & destidT,        &  ! Derivative of esti w.r.t. T
      & dfdT               ! Derivative of f w.r.t. T
    
    INTEGER :: i
    REAL(wp) :: aw,bw,cw,dw,ai,bi,ci,di,AAw,BBw,CCw,AAi,BBi,CCi,alpha,beta

    !CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:calc_atm_fluxes_from_bulk'
    !-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start' )

    Tsurf(:,:)  = p_os%p_prog(nold(1))%tracer(:,1,:,1)  ! set surface temp = mixed layer temp
    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt               ! Change units of tafoK  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt                    ! Change units of ftdewC to C


    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta) 
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/); 
    ! updated from Buck, A. L., New equations for computing vapor pressure and 
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981" 
    !-----------------------------------------------------------------------

    aw=611.21_wp; bw=18.729_wp; cw=257.87_wp; dw=227.3_wp
    ai=611.15_wp; bi=23.036_wp; ci=279.82_wp; di=333.7_wp

    AAw=7.2e-4_wp; BBw=3.20e-6_wp; CCw=5.9e-10_wp
    AAi=2.2e-4_wp; BBi=3.83e-6_wp; CCi=6.4e-10_wp

    alpha=0.62197_wp; beta=0.37803_wp

    fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*ftdewC(:,:)**2)
    esta(:,:) = fa(:,:) * aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
    fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*Tsurf(:,:) **2)
    estw(:,:) = fw(:,:) *aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))
    ! For a given surface salinity we should multiply estw with  1 - 0.000537*S
   
    sphumida(:,:)  = alpha * esta(:,:)/(p_as%pao(:,:)-beta*esta(:,:))
    sphumidw(:,:)  = alpha * estw(:,:)/(p_as%pao(:,:)-beta*estw(:,:))

    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to 
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !-----------------------------------------------------------------------

    humi(:,:)    = 0.39_wp - 0.05_wp*SQRT(esta(:,:)/100._wp)
    fakts(:,:)   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * p_as%fclou(:,:)**2
    Qatm%LWin(:,:) = fakts(:,:) * humi(:,:) * zemiss_def*StBo * tafoK(:,:)**4

    Qatm%LWoutw(:,:) = 4._wp*zemiss_def*StBo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
    Qatm%LWnetw(:,:) = Qatm%LWin(:,:) - Qatm%LWoutw(:,:)

    Qatm%SWin(:,:) = p_as%fswr(:,:)

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to 
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002: 
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------    
    rhoair(:,:)     = p_as%pao(:,:)/(rd*tafoK(:,:)*(1.0_wp+0.61_wp*sphumida(:,:)) )
    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,p_as%fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))
    dragl(:,:)      = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-p_as%tafo(:,:))
    ! Need to keep the drag honest
    dragl(:,:)      = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
    drags(:,:)      = 0.95_wp * dragl(:,:)
    Qatm%sensw(:,:) = drags(:,:)*rhoair(:,:)*cpd*p_as%fu10(:,:) * (p_as%tafo(:,:) -Tsurf(:,:))
    Qatm%latw(:,:)  = dragl(:,:)*rhoair(:,:)*alv*p_as%fu10(:,:) &
      &               * (sphumida(:,:)-sphumidw(:,:))

    DO i = 1, p_ice%kice
      WHERE (p_ice% isice(:,i,:))
        Tsurf(:,:)    = p_ice%Tsurf(:,i,:)
        fi(:,:)       = 1.0_wp+AAi+p_as%pao(:,:)*(BBi+CCi*Tsurf(:,:) **2)
        esti(:,:)     = fi(:,:)*ai*EXP((bi-Tsurf(:,:) /di)*Tsurf(:,:) /(Tsurf(:,:) +ci))
        sphumidi(:,:) = alpha*esti(:,:)/(p_as%pao(:,:)-beta*esti(:,:))
        ! This may not be the best drag parametrisation to use over ice
        dragl(:,:)    = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-p_as%tafo(:,:))
        ! Need to keep the drag honest 
        dragl(:,:)    = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
        drags(:,:)    = 0.95_wp * dragl(:,:)

        ! #achim: ?!
        Qatm%LWout (:,i,:)  = 4._wp*zemiss_def*StBo*tafoK(:,:)**3 * (Tsurf(:,:) &
          &                    - p_as%tafo(:,:))
        Qatm%LWnet (:,i,:)  = Qatm%LWin(:,:) - Qatm%LWout(:,i,:)
        Qatm%dLWdT (:,i,:)  = 4._wp*zemiss_def*StBo*tafoK(:,:)**3
        Qatm%sens  (:,i,:)  = drags(:,:) * rhoair(:,:)*cpd*p_as%fu10(:,:) &
          &                    * (p_as%tafo(:,:) -Tsurf(:,:))
        Qatm%lat   (:,i,:)  = dragl(:,:) * rhoair(:,:)* alf *p_as%fu10(:,:) &
          &                   * (sphumida(:,:)-sphumidi(:,:))

        Qatm%dsensdT(:,i,:) = 0.95_wp*cpd*rhoair(:,:)*p_as%fu10(:,:)&
          &                  *(dragl0(:,:) - 2.0_wp*dragl(:,:))
        dsphumididesti(:,:) = alpha/(p_as%pao(:,:)-beta*esti(:,:)) &
          &                   * (1.0_wp + beta*esti(:,:)/(p_as%pao(:,:)-beta*esti(:,:)))
        destidT(:,:)        = (bi*ci*di-Tsurf(:,:)*(2.0_wp*ci+Tsurf(:,:)))&
          &                   /(di*(ci+Tsurf(:,:))**2) * esti(:,:)
        dfdT(:,:)               = 2.0_wp*CCi*BBi*Tsurf(:,:)
        Qatm%dlatdT(:,i,:)  = alf*rhoair(:,:)*p_as%fu10(:,:)* &
          &                  ( (sphumida(:,:)-sphumidi(:,:))*dragl1(:,:) &
          &                    - dragl(:,:)*dsphumididesti(:,:)*(fi(:,:)*destidT(:,:) &
          &                    + esti(:,:)*dfdT(:,:)) )
      ENDWHERE
    ENDDO

    !Dirk: why zero ?
    Qatm%rpreci(:,:) = 0.0_wp
    Qatm%rprecw(:,:) = 0.0_wp

!!$    !#achim: check stuff
!!$    CALL print_cells(p_as%tafo(:,:),'p_as%tafo')
!!$    CALL print_cells(Qatm%sens(:,1,:),'Qatm%sens')
!!$    CALL print_cells(Qatm%dsensdT(:,1,:),'Qatm%dsensdT')
!!$    CALL print_cells(Qatm%dLWdT(:,1,:),'Qatm%dLWdT')

  END SUBROUTINE calc_atm_fluxes_from_bulk
 
  !-------------------------------------------------------------------------

  SUBROUTINE prepare4restart(p_ice)
    TYPE (t_sea_ice),  INTENT(INOUT) :: p_ice

    WHERE (p_ice%isice(:,:,:))
      p_ice%restart_isice(:,:,:) = 1.0_wp
    ELSEWHERE
      p_ice%restart_isice(:,:,:) = 0.0_wp
    ENDWHERE
  END SUBROUTINE prepare4restart

  SUBROUTINE prepareAfterRestart(p_ice)
    TYPE (t_sea_ice),  INTENT(INOUT) :: p_ice

    IF (is_restart_run()) THEN
      WHERE (p_ice%restart_isice(:,:,:) < 0.5_wp)
        p_ice%isice(:,:,:) = .FALSE.
      ELSEWHERE
        p_ice%isice(:,:,:) = .TRUE.
      ENDWHERE
    END IF
  END SUBROUTINE prepareAfterRestart

END MODULE mo_sea_ice
