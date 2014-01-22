!>
!! Contains the implementation of the initial conditions for the hydrostatic ocean model.
!!
!! Contains the implementation of the initial conditions for the hydrostatic ocean model.
!! This module controls the initial conditions as well as the initialisation of
!! test cases, the top and bottom boundary conditions, and the structure of the
!! forcing quantities.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!
!! @par Revision History
!! Initial version  by Peter Korn (MPI-M)  (2006).
!! Modified by Stephan Lorenz     (MPI-M)  (2010-06).
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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

MODULE mo_ocean_initial_conditions
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_grid_config,        ONLY: nroot, n_dom, grid_sphere_radius, grid_angular_velocity
  USE mo_physical_constants, ONLY: rgrav, sal_ref, sfc_press_bar, tmelt, tf! , SItodBar, rho_ref
  USE mo_math_constants,     ONLY: pi, pi_2, rad2deg, deg2rad
  USE mo_parallel_config,    ONLY: nproma
  USE mo_ocean_nml,          ONLY: iswm_oce, n_zlev, no_tracer, itestcase_oce, i_sea_ice,     &
    & basin_center_lat, basin_center_lon, discretization_scheme,           &
    & basin_height_deg,  basin_width_deg,  use_file_initialConditions,         &
    & initial_temperature_bottom, initial_temperature_top, &
    & initial_salinity_top, initial_salinity_bottom, &
    & use_tracer_x_height, topography_type, topography_height_reference, &
    & sea_surface_height_type, initial_temperature_type, initial_salinity_type, &
    & initial_sst_type, initial_velocity_type, initial_velocity_amplitude

  USE mo_impl_constants,     ONLY: max_char_length, sea, sea_boundary, boundary, land,        &
    & land_boundary,                                             &
    & oce_testcase_zero, oce_testcase_init, oce_testcase_file! , MIN_DOLIC
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, t_geographical_coordinates
  !USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_netcdf_read,        ONLY: read_netcdf_data
  USE mo_sea_ice_types,      ONLY: t_sfc_flx
  USE mo_oce_types,          ONLY: t_hydro_ocean_state
  USE mo_scalar_product,     ONLY: calc_scalar_product_veloc_3d !, map_edges2cell_3D
  USE mo_oce_math_operators, ONLY: grad_fd_norm_oce_3d
  USE mo_oce_thermodyn,      ONLY: convert_insitu2pot_temp_func
  USE mo_oce_ab_timestepping,ONLY: update_time_indices! , calc_vert_velocity
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_ape_params,         ONLY: ape_sst
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  IMPLICIT NONE
  PRIVATE
  INCLUDE 'netcdf.inc'
  !
  PUBLIC :: apply_initial_conditions
  
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=12)           :: str_module    = 'oceInit     '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  
  REAL(wp) :: sphere_radius, u0
  REAL(wp), PARAMETER :: aleph = 0.0_wp
  
  CHARACTER(LEN=*), PARAMETER :: module_name = 'ocean_initial_conditions'

  ! Should be replaced by reading a file
  REAL(wp), PARAMETER :: tprof(20)=&
    & (/ 18.13_wp, 17.80_wp, 17.15_wp, 16.09_wp, 15.04_wp, 13.24_wp, 11.82_wp,  9.902_wp, &
    & 8.484_wp, 7.341_wp, 5.727_wp, 4.589_wp, 3.807_wp, 3.062_wp, 2.481_wp, 2.194_wp, &
    & 1.789_wp, 1.266_wp, 1.070_wp, 0.9211_wp /)

    ! temperature profile for 4-20 layers for testcase 40/45 and similar
    REAL(wp), PARAMETER :: tprof_var(20)= &
      & (/ 25.0_wp, 23.0_wp, 20.0_wp, 15.0_wp, 10.0_wp, 8.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp,&
      & 2.0_wp,  1.0_wp,  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)
    !   & (/ 25.0_wp, 18.0_wp, 12.0_wp, 8.0_wp, 6.0_wp, 4.0_wp, 2.0_wp, 1.0_wp, 0.5_wp, 0.0_wp,&
    !   &     0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)

    REAL(wp) , PARAMETER :: tprof_4layerstommel(4) = (/20.0_wp,10.0_wp,8.0_wp,6.0_wp/)

  
CONTAINS
  

  !-------------------------------------------------------------------------
  SUBROUTINE apply_initial_conditions(patch_2d, patch_3d, ocean_state, external_data, &
    & operators_coeff)
    TYPE(t_patch),TARGET,INTENT(in)   :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_external_data)             :: external_data
    TYPE(t_operator_coeff)            :: operators_coeff
    ! TYPE(t_sfc_flx)                   :: p_sfc_flx

    sphere_radius = grid_sphere_radius
    u0 = (2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

    IF (use_file_initialConditions) THEN
      CALL init_ocean_fromFile(patch_2d, patch_3d, ocean_state)

    ELSE

      CALL init_ocean_bathymetry(patch_3d=patch_3d,  cells_bathymetry=external_data%oce%bathymetry_c(:,:))
      CALL init_ocean_surface_height(patch_3d=patch_3d, ocean_height=ocean_state%p_prog(nold(1))%h(:,:))
      CALL init_ocean_velocity(patch_3d=patch_3d, normal_velocity=ocean_state%p_prog(nold(1))%vn)

      IF (no_tracer > 0) &
        & CALL init_ocean_temperature(patch_3d=patch_3d, ocean_temperature=ocean_state%p_prog(nold(1))%tracer(:,:,:,1))
      IF (no_tracer > 1) &
        & CALL init_ocean_salinity(patch_3d=patch_3d, ocean_salinity=ocean_state%p_prog(nold(1))%tracer(:,:,:,2))

    END IF

    CALL initialize_diagnostic_fields( patch_2d, patch_3d, ocean_state, operators_coeff)
    CALL fill_tracer_x_height(patch_3d, ocean_state)


  END SUBROUTINE apply_initial_conditions
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Initialization of prognostic variables for the hydrostatic ocean model.
  !! Temperature and salinity are read from external data
  !! Finally the prognostic state should be initialized from some restart file.
  !
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M, 2011-09
  !-------------------------------------------------------------------------
  SUBROUTINE init_ocean_fromFile(patch_2d, patch_3d, ocean_state)
    TYPE(t_patch),TARGET, INTENT(in)  :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    
    ! Local Variables
    CHARACTER(LEN=max_char_length), PARAMETER :: method_name = 'mo_ocean_initial_conditions:init_ocean_fromFile'
    CHARACTER(filename_max) :: prog_init_file   !< file name for reading in
    
    LOGICAL :: l_exist
    INTEGER :: i_lev, no_cells, no_levels, jk, jb, jc
    INTEGER :: ncid, dimid
    !INTEGER :: i_startblk_c, i_endblk_c, start_cell_index, end_cell_index, rl_start, rl_end_c
    INTEGER :: start_cell_index, end_cell_index
    
    REAL(wp):: z_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: z_prog(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: all_cells
    !-------------------------------------------------------------------------
    
    CALL message (TRIM(method_name), 'start')
    
    all_cells => patch_2d%cells%ALL
    
    IF(my_process_is_stdio()) THEN
      !
      ! Prognostic variables are read from prog_init_file
      i_lev = patch_2d%level
      WRITE (prog_init_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-prog.nc'
      !prog_init_file="/scratch/local1/m212053/ICON/trunk/icon-dev/grids/&
      !&ts_phc_annual-iconR2B04-L10_50-1000m.nc"
      
      INQUIRE (FILE=prog_init_file, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
        WRITE(message_text,'(3a)') 'netcdf file named ', TRIM(prog_init_file), ' not found!'
        CALL message(TRIM(method_name),TRIM(message_text))
        CALL finish(TRIM(method_name),'File for reading ocean prognostic input not found - ABORT')
      ENDIF
      
      WRITE(message_text,'(3a)') 'netcdf file named ', TRIM(prog_init_file), ' opened for reading'
      CALL message(TRIM(method_name),TRIM(message_text))
      
      ! open file
      CALL nf(nf_open(TRIM(prog_init_file), nf_nowrite, ncid))
      
      ! get number of cells
      CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
      CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
      
      ! check the number of cells
      WRITE(message_text,'(a,i6)') 'No of cells =', no_cells
      CALL message(TRIM(method_name),TRIM(message_text))
      IF(patch_2d%n_patch_cells_g /= no_cells) THEN
        CALL finish(TRIM(method_name),&
          & 'Number of patch cells and cells in ocean prognostic input file do not match - ABORT')
      ENDIF

      ! get number of levels
      CALL nf(nf_inq_dimid(ncid, 'level', dimid))
      CALL nf(nf_inq_dimlen(ncid, dimid, no_levels))
      
      ! check the number of levels
      WRITE(message_text,'(a,i6)') 'No of vertical levels =', no_levels
      CALL message(TRIM(method_name),TRIM(message_text))
      IF(n_zlev /= no_levels) THEN
        CALL finish(TRIM(method_name),&
          & 'Number of vertical levels and &
          & levels in ocean prognostic input file do not match - ABORT')
      ENDIF
      
    ENDIF

    !-------------------------------------------------------
    !
    ! Read ocean init data at cells
    !
    !-------------------------------------------------------
    ! triangle center and edges
    
    ! read temperature
    !  - 2011-11-01, >r7005: read one data set, annual mean only
    !  - "T": annual mean temperature
    CALL read_netcdf_data (ncid, 'T', patch_2d%n_patch_cells_g, patch_2d%n_patch_cells, &
      & patch_2d%cells%decomp_info%glb_index, n_zlev, z_prog)
    
    IF (no_tracer>=1) THEN
      !ocean_state%p_prog(nold(1))%tracer(:,1:n_zlev,:,1) = z_prog(:,1:n_zlev,:)
      ocean_state%p_diag%temp_insitu(:,1:n_zlev,:) = z_prog(:,1:n_zlev,:)
    ELSE
      CALL message( TRIM(method_name),'WARNING: no tracer used, but init temperature attempted')
    END IF

    ! read salinity
    !  - "S": annual mean salinity
    IF (no_tracer > 1) THEN
      CALL read_netcdf_data (ncid, 'S', patch_2d%n_patch_cells_g, patch_2d%n_patch_cells, &
        & patch_2d%cells%decomp_info%glb_index, n_zlev, z_prog)
      ocean_state%p_prog(nold(1))%tracer(:,1:n_zlev,:,2) = z_prog(:,1:n_zlev,:)
    END IF
    
    IF (no_tracer >=2) THEN
      DO jk=1, n_zlev
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            
            ! set values on land to zero/reference
            IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
              
              ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1)=&
                & convert_insitu2pot_temp_func(ocean_state%p_diag%temp_insitu(jc,jk,jb),&
                & ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2),&
                & sfc_press_bar)
            ENDIF
          END DO
        END DO
      END DO
    ELSEIF(no_tracer==1)THEN
      ocean_state%p_prog(nold(1))%tracer(:,1:n_zlev,:,1)=ocean_state%p_diag%temp_insitu(:,1:n_zlev,:)
    ENDIF
    
    ! close file
    IF(my_process_is_stdio()) CALL nf(nf_close(ncid))
    !---------------------------------------------------
    
    DO jk=1, n_zlev
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          
          ! set values on land to zero/reference
          IF ( patch_3d%lsm_c(jc,jk,jb) > sea_boundary ) THEN
            ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = 0.0_wp
            IF (no_tracer>=2) ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2) = 0.0_wp !sal_ref
          ENDIF
          
        END DO
      END DO
    END DO
    
    !---------Debug Diagnostics-------------------------------------------
    idt_src=0  ! output print level - 0: print in any case
    z_c(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    CALL dbg_print('init prognostic - T'       ,z_c                     ,str_module,idt_src)
    IF (no_tracer > 1) THEN
      z_c(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
      CALL dbg_print('init prognostic - S'       ,z_c                   ,str_module,idt_src)
    END IF
    !---------------------------------------------------------------------
    
    CALL message( TRIM(method_name),'Ocean prognostic initialization data read' )
    
  END SUBROUTINE init_ocean_fromFile
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  SUBROUTINE nf(STATUS)
    
    INTEGER, INTENT(in) :: STATUS
    
    IF (STATUS /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(STATUS))
    ENDIF
    
  END SUBROUTINE nf
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE initialize_diagnostic_fields( patch_2d,patch_3d, ocean_state, operators_coeff)
    TYPE(t_patch), TARGET, INTENT(in)             :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff)                        :: operators_coeff

    IF(discretization_scheme==1)THEN

      IF (is_restart_run()) CALL update_time_indices(1)

      CALL calc_scalar_product_veloc_3d( patch_3d,&
        & ocean_state%p_prog(nold(1))%vn,&
        & ocean_state%p_prog(nold(1))%vn,&
        & ocean_state%p_diag,            &
        & operators_coeff)
      
      IF (is_restart_run()) CALL update_time_indices(1)
    ELSE
    ENDIF
    
    !---------Debug Diagnostics-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('recon_fields: p_vn%x(1)'        ,ocean_state%p_diag%p_vn%x(1)  ,str_module,idt_src)
    !---------------------------------------------------------------------
    
    !    IF (.NOT. is_restart_run()) CALL calc_vert_velocity( patch_2D, ocean_state, operators_coeff)
    
  END SUBROUTINE initialize_diagnostic_fields
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! Initialization of test cases for the hydrostatic ocean model.
  !! Currently only some simple test value are set.
  !! Finally the prognostic state should be initialized from some restart file.
  !
  !! @par Revision History
  !! Developed  by Peter Korn, MPI-M, 2006-08
  !
  !-------------------------------------------------------------------------
  SUBROUTINE init_ocean_analytically(patch_3d, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_external_data)             :: external_data
    ! TYPE(t_sfc_flx)                   :: p_sfc_flx
    ! Local Variables
    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: z_dolic
    REAL(wp), ALLOCATABLE :: z_c(:,:,:)
    REAL(wp):: cell_lat, cell_lon
    REAL(wp):: distan, lat_deg, lon_deg, z_tmp
    REAL(wp):: perturbation_lon, perturbation_lat, max_perturbation, perturbation_width !,z_H_0
    REAL(wp):: z_ttrop, z_tpol, z_tpols, z_tdeep, z_tdiff, z_ltrop, z_lpol, z_ldiff
    REAL(wp):: z_temp_max, z_temp_min, z_temp_incr, z_max
    REAL(wp):: t,s
    CHARACTER(LEN=max_char_length) :: initial_sst_type

    REAL(wp) :: salinity_profile(n_zlev)

    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    
    TYPE(t_patch), POINTER   :: patch_2d
    CHARACTER(LEN=max_char_length), PARAMETER :: method_name = 'mo_ocean_initial_conditions:init_ocean_analytically'

    TYPE(t_subset_range), POINTER :: all_cells, owned_cells, all_edges
    !-------------------------------------------------------------------------

    CALL message (TRIM(method_name), 'start')
    
    sphere_radius = grid_sphere_radius
    u0 =(2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

    !-------------------------------------------------------------------------
    
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    owned_cells => patch_2d%cells%owned
    all_edges => patch_2d%edges%ALL
    cell_center => patch_2d%cells%center
    ALLOCATE(z_c(nproma,n_zlev,patch_2d%alloc_cell_blocks))

    ! initialize salinity with reference value rather than with zero
    !  - mainly for plotting purpose
    IF ( no_tracer >= 2) THEN
      ocean_state%p_prog(nold(1))%tracer(:,:,:,2) = sal_ref
    END IF
    
    !IF shallow-water option is NOT selected then)
    IF ( iswm_oce /= 1 )THEN
      
      SELECT CASE (itestcase_oce)
      
      CASE (oce_testcase_zero)
        CALL message(TRIM(method_name), 'you have selected the "no-testcase" option')
      CASE (oce_testcase_init)
        
      CASE (oce_testcase_file)
        CALL finish(TRIM(method_name), 'Initialization from file NOT SUPPORTED YET - TERMINATE')
        !CALL init_from_file(patch_2D)
        
      CASE (30, 31)
        CALL message(TRIM(method_name), 'Simple Initialization of testcases (30, 31)')
        CALL message(TRIM(method_name), ' - here: horizontally homogen, vertical profile for T and S')
        CALL message(TRIM(method_name), ' - Add forcing / restoring / wave for dynamic test')
        
        !Ocean at rest
        ocean_state%p_prog(nold(1))%vn(:,:,:) = 0.0_wp
        
        !init temperature and salinity with vertical profiles
        ! use initial_salinity_type    = 401
        ! use initial_temperature_type = 401

        IF (itestcase_oce == 31) THEN
          CALL message(TRIM(method_name), 'Simple Initialization of testcases (31)')
          CALL message(TRIM(method_name), ' - here: external gravity wave')
          
          ! set sea_surface_height_type = 201

        END IF
        
      CASE (33)
        ! collapsing density front testcase, taken from Stuhne-Peltier (JCP, 2006)
        CALL message(TRIM(method_name), 'Initialization of testcases (33)')
        CALL message(TRIM(method_name), ' - here: Collapsing density front, Stuhne-Peltier')
        ! use initial_salinity_type    = 401
        ! use initial_temperature_type = 201
        

      CASE (34)
        ! Adjusting density front in a basin: vertical wall at basin_center_lon
        CALL message(TRIM(method_name), 'Initialization of testcases (34)')
        CALL message(TRIM(method_name),' - here: Adjusting density front in a basin with vertical wall')
        ! use initial_temperature_type = 202
        
        
      CASE (32) !from Sergy Danilov
        CALL message(TRIM(method_name), 'Simple Initialization of testcases (32)')
        CALL message(TRIM(method_name), ' - here: Danilovs Munk gyre flow')
        ! use initial_temperature_type = 203
        ! use sea_surface_height_type = 202
        
        !p_pos%lon = 0.0_wp
        !p_pos%lat = 0.0_wp
        
        !       DO jb = i_startblk_c, i_endblk_c
        !         CALL get_indices_c(patch_2D, jb, i_startblk_c, i_endblk_c, &
        !          &                start_cell_index, end_cell_index, rl_start, rl_end_c)
        !         DO jc = start_cell_index, end_cell_index
        !
        !           ocean_state%p_prog(nold(1))%tracer(jc,:,jb,1)=20.0_wp
        !
        !           cell_lat = patch_2D%cells%center(jc,jb)%lat
        !           cell_lon = patch_2D%cells%center(jc,jb)%lon
        !
        !           z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
        !           IF (z_dolic > 0) THEN
        !             ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
        !             ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
        !             ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
        !             ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375
        !             DO jk = 1, z_dolic
        !                ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) &
        !               & = 20.0_wp+0.1_wp*v_base%zlev_m(jk)/v_base%zlev_m(z_dolic)
        !               ! ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) &
        !               !   & = 20.0_wp -  v_base%zlev_m(jk)*15.0_wp/4000.0_wp
        !             END DO
        !           END IF
        !         END DO
        !       END DO
        ! DO jk = 1, n_zlev
        ! write(*,*)'Temperature strat',jk,&
        ! &maxval(ocean_state%p_prog(nold(1))%tracer(:,jk,:,1)),&
        ! &minval(ocean_state%p_prog(nold(1))%tracer(:,jk,:,1))
        ! END DO
        

          !         !After hot spot now a cool spot at a slightly different location
          !         ! Add temperature perturbation at new values - 35N; 10W
          !         perturbation_lat = basin_center_lat - 0.1_wp*basin_height_deg!             !45.5_wp
          !         perturbation_lon =  -0.1_wp*basin_width_deg                                 !4.5_wp
          !         max_perturbation  = 10.0_wp!20.1_wp
          !         perturbation_width  =  5.0_wp!1.5_wp
          ! !         perturbation_lat  = 25.0_wp
          ! !         perturbation_lon  = 8.0_wp
          ! !         max_perturbation  = -2.0_wp
          ! !         perturbation_width  =  5.0_wp
          !         DO jb = i_startblk_c, i_endblk_c
          !           CALL get_indices_c(patch_2D, jb, i_startblk_c, i_endblk_c, &
          !            &                start_cell_index, end_cell_index, rl_start, rl_end_c)
          !           DO jc = start_cell_index, end_cell_index
          !             cell_lat = patch_2D%cells%center(jc,jb)%lat
          !             cell_lon = patch_2D%cells%center(jc,jb)%lon
          !             z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
          !             IF (z_dolic > 0) THEN
          !               ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
          !               ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
          !               ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
          !               ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375
          !               distan=sqrt((cell_lat-perturbation_lat*deg2rad)**2+(cell_lon-perturbation_lon*deg2rad)**2)
          !               !write(123,*)'zdist',cell_lat,cell_lon,distan,10.5_wp*deg2rad
          !               !IF(distan<=25.5_wp*deg2rad)cycle
          !               ! at distan > 25.5 degrees:
          !               !  e.g. at 30 deg distance the added perturbation would be ~ exp(-400) ~ 0.0
          !               ! Now without cycle in loop - perturbation is very small at distance>10 deg
          !               !  e.g. at 3 deg distance is
          !               !   T(jk=1)=19.0625+20.1*exp(-4)*sin(pi* 250/4000) = 19.06 + 20.1*0.18*0.06 = 19.28
          !               !   T(jk=4)=13.4375+20.1*exp(-4)*sin(pi*1750/4000) = 13.44 + 20.1*0.18*0.42 = 15.00
          !               DO jk = 1, z_dolic
          !                ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) =          &
          !                 & ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1)          &
          !                 &   - max_perturbation*exp(-(distance/(perturbation_width*deg2rad))**2) &
          !                 &   * sin(pi*v_base%zlev_m(jk)/v_base%zlev_i(z_dolic+1))
          !               ENDDO
          !             END IF
          !           END DO
          !         END DO
          

        !----------------Old version of 32: please retain code, its also interesting
        !----------------An old version of surface forcing corresponds to this
        ! !     CASE (32) !from Sergy Danilov
        ! !       DO jb = i_startblk_c, i_endblk_c
        ! !         CALL get_indices_c(patch_2D, jb, i_startblk_c, i_endblk_c, &
        ! !          &                start_cell_index, end_cell_index, rl_start, rl_end_c)
        ! !         DO jc = start_cell_index, end_cell_index
        ! !           cell_lat = patch_2D%cells%center(jc,jb)%lat
        ! !           cell_lon = patch_2D%cells%center(jc,jb)%lon
        ! !           z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
        ! !           IF (z_dolic > 0) THEN
        ! !             ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
        ! !             ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
        ! !             ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
        ! !             ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375
        ! !             DO jk = 1, z_dolic
        ! !                ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) &
        ! !               & = 20.0_wp-v_base%zlev_m(jk)*15.0_wp/v_base%zlev_i(z_dolic+1)
        ! !             END DO
        ! !           END IF
        ! !         END DO
        ! !       END DO
        ! !       perturbation_lat = basin_center_lat + 0.1_wp*basin_height_deg!             !45.5_wp
        ! !       perturbation_lon =  0.1_wp*basin_width_deg                                 !4.5_wp
        ! !       max_perturbation  = 10.0_wp!20.1_wp
        ! !       perturbation_width  =  5.0_wp!1.5_wp
        ! !       IF (no_tracer > 0 ) THEN
        ! !         DO jb = i_startblk_c, i_endblk_c
        ! !           CALL get_indices_c(patch_2D, jb, i_startblk_c, i_endblk_c, &
        ! !            &                start_cell_index, end_cell_index, rl_start, rl_end_c)
        ! !           DO jc = start_cell_index, end_cell_index
        ! !             cell_lat = patch_2D%cells%center(jc,jb)%lat
        ! !             cell_lon = patch_2D%cells%center(jc,jb)%lon
        ! !             z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
        ! !             IF (z_dolic > 0) THEN
        ! !               distance=sqrt((cell_lat-perturbation_lat*deg2rad)**2+(cell_lon-perturbation_lon*deg2rad)**2)
        ! !               !Local hot perturbation
        ! !               !IF(distance<=5.0_wp*deg2rad)THEN
        ! !               DO jk = 1, z_dolic
        ! !                  ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) =          &
        ! !                  & ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1)          &
        ! !                  &   + max_perturbation*exp(-(distance/(perturbation_width*deg2rad))**2) &
        ! !                  &   * sin(pi*v_base%zlev_m(jk)/v_base%zlev_i(z_dolic+1))
        ! !               END DO
        ! !               !ENDIF
        ! !             END IF
        ! !           END DO
        ! !         END DO
        ! !       END IF ! no_tracer > 0
        !----------------End of old version------------------------------------
        
        
      CASE (40)
        ! Temperature profile depends on latitude and depth
        ! Construct temperature profile
        !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
        !   for maximum tropical temperature see values above
        CALL message(TRIM(method_name), 'Simple Initialization of testcases (40)')
        CALL message(TRIM(method_name), ' - here: simple tropics-pol/vertical temperature profile')
        ! use initial_temperature_type = 204
        
        
      CASE (41)
        ! #slo 2011-10-05#
        !  incorrect (for n_zlev>9) old testcase 40 with z_tpol=0.0 at poles saved for reference
        CALL message(TRIM(method_name), 'Simple Initialization of testcases (41)')
        CALL message(TRIM(method_name), ' - here: old erroneous profile saved for reference')
        
!        ! Temperature profile depends on latitude and depth
!        ! Construct temperature profile
!        !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
!        !   for maximum tropical temperature see values above
!        z_tpol  =  5.0_wp      ! polar temperature - old testcase 40 for n_zlev<6
!        !  z_tpol  =  0.0_wp      ! polar temperature - new #slo# debug
!        z_ltrop = 15.0_wp      ! tropical latitude for temperature gradient
!        z_lpol  = 60.0_wp      ! polar latitude for temperature gradient
!        z_ldiff = z_lpol  - z_ltrop
!
!        DO jb = all_cells%start_block, all_cells%end_block
!          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
!          DO jc = start_cell_index, end_cell_index
!
!            !latitude given in radians
!            cell_lat = patch_2d%cells%center(jc,jb)%lat
!            !transer to latitude in degrees
!            lat_deg = cell_lat*rad2deg
!
!            DO jk=1,n_zlev
!              IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
!
!                ! set maximum tropical temperature from profile above
!                z_ttrop = tprof_var(jk)
!                z_tpol  = MIN(z_tpol,tprof_var(jk))
!                z_tdiff = z_ttrop - z_tpol
!
!                !constant salinity
!                IF(no_tracer==2)THEN
!                  ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2) = sprof(jk) !35.0_wp
!                ENDIF
!
!                IF(ABS(lat_deg)>=z_lpol)THEN
!
!                  ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_tpol
!                  ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = ocean_state%p_diag%temp_insitu(jc,jk,jb)
!
!                ELSEIF(ABS(lat_deg)<=z_ltrop)THEN
!
!                  ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop
!                  ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = ocean_state%p_diag%temp_insitu(jc,jk,jb)
!
!
!                ELSEIF(ABS(lat_deg)<z_lpol .AND. ABS(lat_deg)>z_ltrop)THEN
!                  z_tmp = 0.5_wp*pi*((ABS(lat_deg) - z_ltrop)/z_ldiff)
!                  ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop - z_tdiff*SIN(z_tmp)
!                  ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = ocean_state%p_diag%temp_insitu(jc,jk,jb)
!                ENDIF
!              ENDIF   ! lsm
!            END DO
!          END DO
!        END DO
        
      CASE (43)
        ! #slo# collapsing density front with much weaker gradient than testcase 33
        ! without temperature restoring / relaxation
        CALL message(TRIM(method_name), 'Initialization of testcases (43)')
        CALL message(TRIM(method_name), ' - here: Collapsing density front with weaker gradient')
        
        ! use initial_temperature_type = 205
        ! use initial_salinity_type = 401

      CASE (44)
        ! Temperature is homogeneous in each layer. Varies from 30.5 in top to 0.5 in bottom layer
        CALL message(TRIM(method_name), 'Initialization of testcases (44)')
        CALL message(TRIM(method_name), ' - here: horizontally homogen, stable vertical profile')

        ! use initial_temperature_type = 200
        
      CASE (45)
        ! T and S are horizontally homegeneous. Values are taken from t_prof[_var] and s_prof[_var]
        CALL message(TRIM(method_name), 'Initialization of testcases (45)')
        CALL message(TRIM(method_name), &
          & ' - here: horizontally homogen, use tprof_var and salinity_profile_20levels vertical profiles')
        ! use initial_temperature_type = 402
        ! use initial_salinity_type = 402
        
      CASE (46,461)
        ! T and S are horizontally and vertically homegeneous
        ! Values are taken from namelist and used for comparison with MPIOM; default: t=16 C, s=35 psu
        CALL message(TRIM(method_name), 'Initialization of testcases (46)')
        CALL message(TRIM(method_name), &
          & ' - here: horizontally and vertically homogen')
        ! use initial_temperature_type = 200
        ! use initial_salinity_type = 200

      CASE (47)
        ! T and S are horizontally and vertically homegeneous
        ! include some special init - here Indonesia set to warm/salty surface
        CALL message(TRIM(method_name), 'Initialization of testcases (47)')
        CALL message(TRIM(method_name), &
          & ' - here: horizontally and vertically homogen+warm/salty Indonesia')
        
        ! use initial_temperature_type = 206
        ! use initial_salinity_type = 201

        
      CASE (50)
        ! Testcase for coupled Aquaplanet:
        !  - following APE_ATLAS Equations (2.1) - (2.5)
        !  - use function ape_sst for initializing SST
        !  - decrease maximum temperature vertically by z_temp_incr
        !  - use parameter 'sst_qobs' - maximum temperature = 27, minimum polar temperature = 0 deg C
        CALL message(TRIM(method_name), 'Initialization of testcases (50)')
        CALL message(TRIM(method_name), ' - here: testcase for coupled aquaplanet, using sst_qobs')
        
        ! use initial_salinity_type = 402
        ! use initial_temperature_type = 207
        ! use initial_sst_type='sst1'

        ! Important:
        !   use initial_temperature_top=27.0 initial_temperature_bottom=0.0
        !   to be consistent with the old setup
        
      CASE (1050)
        ! as 50, but salinity is analytically calculated
        CALL message(TRIM(method_name), 'Initialization of testcase (1050)')
        CALL message(TRIM(method_name), ' - here: testcase for coupled aquaplanet, using analytic s')
        
        ! use initial_temperature_type = 200 top_temperature = bottom_temperature = 10
        ! use initial_salinity_type = 202    top_salinity = 34.1, bottom_salinity = 35.0


      CASE (52)
        ! Testcase for coupled Aquaplanet:
        !  - following APE_ATLAS Equations (2.1) - (2.5)
        !  - use function ape_sst for initializing SST
        !  - decrease maximum temperature vertically by z_temp_incr
        !  - now warmer init to avoid growing of sea ice:
        !    maximum temperature = 27, minimum polar temperature = 10 deg C
        CALL message(TRIM(method_name), 'Initialization of testcases (52)')
        CALL message(TRIM(method_name), &
          & ' - here: testcase for coupled aquaplanet, using sst_qobs, min=10 deg C')

        ! use initial_salinity_type = 402
        ! use initial_temperature_type = 207
        ! use initial_sst_type='sst_qobs'

        ! Important:
        !   use initial_temperature_top=27.0 initial_temperature_bottom=10.0
        !   to be consistent with the old setup
        
        
       CASE (51)
        CALL message(TRIM(method_name), 'Simple Initialization of testcases (51)')
        CALL message(TRIM(method_name), &
          & ' - here: horizontally varying T with local perturbation')
        
        ! use initial_temperature_type = 208, top_temperature = 30.5,  bottom_temperature = 0,5
        
      CASE(53)
        CALL message(TRIM(method_name), 'LOCK exchange (53)')
        ocean_state%p_prog(nold(1))%vn = 0.0_wp
        ocean_state%p_prog(nnew(1))%vn = 0.0_wp

        ! use initial_temperature_type = 209
        
      CASE default
        WRITE(0,*)'testcase',itestcase_oce
        CALL finish(TRIM(method_name), 'CHOSEN INITIALIZATION NOT SUPPORTED - TERMINATE')
      END SELECT
      
      !---------Debug Diagnostics-------------------------------------------
      idt_src=0  ! output print level - 0: print in any case
      IF (no_tracer >=1) THEN
        z_c(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
        CALL dbg_print('init testcases  - T'       ,z_c                     ,str_module,idt_src, &
          & in_subset=owned_cells)
      END IF
      IF (no_tracer >= 2) THEN
        z_c(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
        CALL dbg_print('init testcases  - S'       ,z_c                     ,str_module,idt_src, &
          & in_subset=owned_cells)
      END IF
      !---------------------------------------------------------------------
      
      ! Shallow water testcases:
    ELSEIF( iswm_oce == 1 )THEN
      
      SELECT CASE (itestcase_oce)
      
      CASE (oce_testcase_zero)
        
        CALL message(TRIM(method_name), 'you have selected the "no-testcase" option')
        
      CASE (24)
        
        CALL message(TRIM(method_name), 'Shallow-Water-Testcase (24)')
        CALL message(TRIM(method_name), ' - here: h and bathy for solid body rotation (Laeuter Test)')
        ! use topography_type = 200, topography_height_reference = 0
        ! use sea_surface_height_type = 203
        ! use initial_velocity_type = 201

      CASE (25)
        CALL message(TRIM(method_name), 'Shallow-Water-Testcase (25)')
        CALL message(TRIM(method_name), ' - here: h and bathy of Williamson Test 2')
        ! use topography_type = 200, topography_height_reference = 0
        ! use sea_surface_height_type = 204
        ! use initial_velocity_type = 202

      CASE (26)
        CALL message(TRIM(method_name), 'Shallow-Water-Testcase (26)')
        CALL message(TRIM(method_name), ' - here: h and bathy of Williamson Test 5')
        
        ! use topography_type = 201 (test5_oro)
        ! use sea_surface_height_type = 205
        ! use initial_velocity_type = 203, initial_velocity_amplitude = 20.0
        
      CASE(27)!temperature ditribution
        CALL message(TRIM(method_name), 'Shallow-Water-Testcase (27)')
        ocean_state%p_prog(nold(1))%vn = 0.0_wp
        ocean_state%p_prog(nnew(1))%vn = 0.0_wp

        IF(no_tracer>0)THEN
          ocean_state%p_prog(nold(1))%tracer(:,1,:,1) = 0.0_wp
          ocean_state%p_prog(nnew(1))%tracer(:,1,:,1) = 0.0_wp
        ELSE
          CALL finish(TRIM(method_name), 'Number of tracers =0 is inappropriate for this test - TERMINATE')
        ENDIF

        ! use topography_type = 200, topography_height_reference = -200
        ! use initial_temperature_type = 210
        
      CASE(28)
        CALL message(TRIM(method_name), 'Shallow-Water-Testcase (28)')

        ! use initial_velocity_type = 203, initial_velocity_amplitude = 20.0

        perturbation_lat = basin_center_lat! + 0.1_wp*basin_height_deg!             !45.5_wp
        perturbation_lon =  0.0_wp!0.1_wp*basin_width_deg                           !4.5_wp
        !max_perturbation  = 20.0_wp            !20.1_wp
        perturbation_width  =  7.0_wp*pi/64.0_wp !10.0_wp!5.0_wp!1.5_wp
        
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            
            IF(patch_3d%lsm_c(jc,1,jb)<=sea_boundary)THEN
              ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1:no_tracer) = 0.0_wp
              ocean_state%p_prog(nnew(1))%tracer(jc,1,jb,1:no_tracer) = 0.0_wp
              
              distan=SQRT((cell_center(jc,jb)%lat - perturbation_lat * deg2rad)**2 + &
                & (cell_center(jc,jb)%lon - perturbation_lon*deg2rad)**2)
              !Local hot perturbation
              IF(distan<=perturbation_width)THEN
                ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1:no_tracer) =        &
                  & (1.0_wp+COS(pi*distan/perturbation_width))/2.0_wp +2.0_wp
              ENDIF
              ocean_state%p_prog(nnew(1))%tracer(jc,1,jb,1:no_tracer)= ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1:no_tracer)
            ENDIF
          END DO
        END DO
        WRITE(*,*)'max/min tracer at initial time',&
          & MAXVAL( ocean_state%p_prog(nold(1))%tracer(:,1,:,1)),&
          & MINVAL( ocean_state%p_prog(nold(1))%tracer(:,1,:,1))
        
      CASE(29)!State at rest, forced by wind
        CALL message(TRIM(method_name), 'Shallow-Water-Testcase (29)')
        
        ocean_state%p_prog(nold(1))%vn(:,1,:) = 0.0_wp
        
      CASE default
        WRITE(0,*)'testcase',itestcase_oce
        CALL finish(TRIM(method_name), 'CHOSEN INITIALIZATION NOT SUPPORTED in SW MODE - TERMINATE')
      END SELECT
    ENDIF  !  iswm_oce

    !---------------------------------------------------
    DEALLOCATE(z_c )

    CALL message (TRIM(method_name), 'end')
    
  END SUBROUTINE init_ocean_analytically
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_salinity(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    ! replace by file
    REAL(wp), PARAMETER :: sprof(20)=&
      & (/  34.699219_wp, 34.798244_wp, 34.904964_wp, 34.976841_wp, 35.027084_wp, &
      & 35.026825_wp, 34.960835_wp, 34.862324_wp, 34.752468_wp, 34.656761_wp, 34.596603_wp,&
      & 34.594128_wp, 34.628601_wp, 34.678772_wp, 34.717495_wp, 34.738304_wp, 34.741512_wp,&
      & 34.738205_wp, 34.729176_wp, 34.723465_wp /)

    REAL(wp), PARAMETER :: salinity_profile_20levels(20)= &
      & (/ 34.5_wp, 34.6_wp, 34.7_wp, 34.8_wp, 34.9_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp,&
      & 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp/)

    REAL(wp) , PARAMETER :: sprof_4layerstommel(4) = &
      & (/34.699219_wp, 34.798244_wp, 34.904964_wp, 34.976841_wp/)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_salinity'
    !-------------------------------------------------------------------------

    IF (no_tracer < 2) RETURN        ! no salinity
    IF (initial_salinity_type < 200) RETURN ! not analytic salinity

    SELECT CASE (initial_salinity_type)
    !------------------------------
    CASE (200)
      ! uniform salinity or vertically linarly increasing
      CALL message(TRIM(method_name), ': horizontally homogenous, vertically linear')
      CALL tracer_VerticallyLinearlyDecreasing(patch_3d=patch_3d, ocean_tracer=ocean_salinity, &
        & top_value=initial_salinity_top, bottom_value=initial_salinity_bottom)

    !------------------------------
    CASE (201)
      CALL salinity_Uniform_SpecialArea(patch_3d, ocean_salinity)

    !------------------------------
    CASE (202)
      CALL salinity_AnalyticSmoothVerticalProfile(patch_3d, ocean_salinity)

    !------------------------------
    CASE (401)
      ! assign from adhoc array values
      IF (n_zlev==4) THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_salinity, VerticalProfileValue=sprof_4layerstommel)
      ELSEIF  (n_zlev <= 20) THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_salinity, VerticalProfileValue=sprof)
      ELSE
        CALL finish(TRIM(method_name), 'Number of vertical levels to small or to big: >=4 and <=20')
      ENDIF

    !------------------------------
    CASE (402)
      IF  (n_zlev <= 20) THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_salinity, VerticalProfileValue=salinity_profile_20levels)
      ELSE
        CALL finish(TRIM(method_name), 'Number of vertical levels > 20')
      ENDIF

    !------------------------------
     CASE default
      CALL finish(method_name, "unknown initial_salinity_type")

    END SELECT

    CALL dbg_print('init_ocean_salinity', ocean_salinity(:,:,:), &
      & str_module,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE init_ocean_salinity
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_temperature(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    REAL(wp) :: temperature_profile(n_zlev)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_temperature'
    !-------------------------------------------------------------------------

    IF (no_tracer < 1) RETURN        ! no temperature
    IF (initial_temperature_type < 200) RETURN ! not analytic temperature

    SELECT CASE (initial_temperature_type)

    !------------------------------
    CASE (200)
      ! uniform or linearly decreasing temperature
      ! Temperature is homogeneous in each layer.
      CALL message(TRIM(method_name), ': horizontally homogenous, vertically linear')
      CALL tracer_VerticallyLinearlyDecreasing(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & top_value=initial_temperature_top, bottom_value=initial_temperature_bottom)

    !------------------------------
    CASE (201)
      CALL temperature_CollapsingDensityFront_StuhnePeltier(patch_3d, ocean_temperature)

    !------------------------------
    CASE (202)
      CALL temperature_BasinWithVerticalWall(patch_3d, ocean_temperature)

    !------------------------------
    CASE (203)
      CALL temperature_DanilovsMunkGyre(patch_3d, ocean_temperature)

    !------------------------------
    CASE (204)
      CALL temperature_TropicsPolar(patch_3d, ocean_temperature)

    !------------------------------
    CASE (205)
      CALL temperature_CollapsingDensityFront_WeakGrad(patch_3d, ocean_temperature)

    !------------------------------
    CASE (206)
      CALL temperature_Uniform_SpecialArea(patch_3d, ocean_temperature)

    !------------------------------
    CASE (207)
      CALL temperature_APE(patch_3d, ocean_temperature)

    !------------------------------
    CASE (208)
      CALL message(TRIM(method_name), ': horizontally non-homogenous, local pertubation')

      ! first create linearly vertically decreasing temperature, uniform horizontally
      CALL tracer_VerticallyLinearlyDecreasing(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & top_value=initial_temperature_top, bottom_value=initial_temperature_bottom)

      !Add horizontal variation
      CALL temperature_AddHorizontalVariation(patch_3d, ocean_temperature)

      !Add local perturbation
      CALL temperature_AddLocalPerturbation(patch_3d, ocean_temperature)

    !------------------------------
    CASE (209)
      CALL temperature_uniform_SeparationAtLon(patch_3d, ocean_temperature, wallLonDeg=0.0_wp)

    !------------------------------
    CASE (210)
      CALL temperature_uniform_SeparationAtLat(patch_3d, ocean_temperature, wallLatDeg=basin_center_lat)

    !------------------------------
    CASE (401)
      ! assign from adhoc array values
      IF(n_zlev==4)THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_temperature, VerticalProfileValue=tprof_4layerstommel)
      ELSE
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_temperature, VerticalProfileValue=tprof)
      ENDIF

    !------------------------------
    CASE (402)
       CALL fill_FromVerticalArrayProfile(patch_3d, ocean_temperature, VerticalProfileValue=tprof_var)

    !------------------------------
     CASE default
      CALL finish(method_name, "unknown initial_temperature_type")

    END SELECT

    CALL dbg_print('init_ocean_temperature', ocean_temperature(:,:,:), &
      & str_module,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE init_ocean_temperature
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_velocity(patch_3d, normal_velocity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: normal_velocity(:,:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_velocity'
    !-------------------------------------------------------------------------

    IF (initial_velocity_type < 200) RETURN ! not analytic velocity

    SELECT CASE (initial_velocity_type)
    !------------------------------
    CASE (200)
      ! uniform velocity
      CALL message(TRIM(method_name), ': uniform velocity')

    !------------------------------
    CASE (201)
      CALL velocity_usbr(patch_3d, normal_velocity)

    !------------------------------
    CASE (202)
      CALL message(TRIM(method_name), 'Williamson Test 2 ')
      CALL velocity_WilliamsonTest_2_5(patch_3d, normal_velocity, velocity_amplitude=u0)

    !------------------------------
    CASE (203)
      CALL message(TRIM(method_name), 'Williamson Test 5 ')
      CALL velocity_WilliamsonTest_2_5(patch_3d, normal_velocity, velocity_amplitude=initial_velocity_amplitude)

    !------------------------------

    !------------------------------
     CASE default
      CALL finish(method_name, "unknown initial_velocity_type")

    END SELECT

 !   CALL dbg_print('init_ocean_salinity', ocean_salinity(:,:,:), &
 !     & str_module,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE init_ocean_velocity
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !> Initial datum for zonal velocity u, test case unsteady solid body
  ! rotation of L\"auter et al.(2007).
  !
  ! Developed by Th.Heinze, DWD, (2007-03)
  !-------------------------------------------------------------------------
  SUBROUTINE velocity_usbr(patch_3d, vn)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index, level
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: t       ! point of time
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: angle1, angle2, edge_vn, COS_angle1, SIN_angle1

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_usbr_u'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL

    t = 0.0_wp
    angle1 = .25_wp * pi
    COS_angle1 = COS(angle1)
    SIN_angle1 = SIN(angle1)

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat

        angle2 = point_lon + grid_angular_velocity * t
        uu = COS(point_lat) * COS_angle1
        uu = uu + COS(angle2) * SIN(point_lat) * SIN_angle1
        uu = u0 * uu

        vv = SIN(angle2) * SIN_angle1
        vv = -1._wp * u0 * vv

        edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
              & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          vn(edge_index, level, edge_block) = edge_vn
        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE  velocity_usbr
  !-----------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------
  SUBROUTINE velocity_WilliamsonTest_2_5(patch_3d, vn, velocity_amplitude)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)
    REAL(wp), INTENT(in) :: velocity_amplitude

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index, level
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: t       ! point of time
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: angle1, angle2, edge_vn, COS_angle1, SIN_angle1

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_WilliamsonTest_2_5'
    !-------------------------------------------------------------------------

    ! CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat

        uu = COS(point_lat) * COS(aleph)
        uu = uu + COS(point_lon) * SIN(point_lat) * SIN(aleph)
        uu = velocity_amplitude * uu

        vv = SIN(point_lon) * SIN(aleph)
        vv = -1._wp * velocity_amplitude * vv

        edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
              & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          vn(edge_index, level, edge_block) = edge_vn
        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE  velocity_WilliamsonTest_2_5
  !-----------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_uniform_SeparationAtLon(patch_3d, ocean_temperature, wallLonDeg)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp), INTENT(in) :: wallLonDeg

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_uniform_SeparationAtLon'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(jc,jb)%lat * rad2deg
        lon_deg = patch_2d%cells%center(jc,jb)%lon * rad2deg

        IF((lon_deg-basin_center_lon) >= wallLonDeg) THEN
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,jk,jb) = 10.0_wp
          ENDDO
        ELSE
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,jk,jb) = 5.0_wp
          ENDDO
        ENDIF

      END DO
    END DO

  END SUBROUTINE temperature_uniform_SeparationAtLon
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_uniform_SeparationAtLat(patch_3d, ocean_temperature, wallLatDeg)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp), INTENT(in) :: wallLatDeg

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_uniform_SeparationAtLon'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        lat_deg = patch_2d%cells%center(jc,jb)%lat
        lon_deg = patch_2d%cells%center(jc,jb)%lon

        IF((lon_deg-basin_center_lon) >= wallLatDeg)THEN
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,jk,jb) = 10.0_wp
          ENDDO
        ELSE
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,jk,jb) = 5.0_wp
          ENDDO
        ENDIF

      END DO
    END DO

  END SUBROUTINE temperature_uniform_SeparationAtLat
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_AddLocalPerturbation(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_AddLocalPerturbation'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        lat_deg = patch_2d%cells%center(jc,jb)%lat
        lon_deg = patch_2d%cells%center(jc,jb)%lon

        IF(ABS(lon_deg) < 2.5_wp .AND. &
           ABS(lat_deg-basin_center_lat) < 0.25_wp*basin_height_deg) THEN

          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,jk,jb) = ocean_temperature(jc,jk,jb) * 1.1_wp
          END DO

        ENDIF

      END DO
    END DO
  END SUBROUTINE temperature_AddLocalPerturbation
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_AddHorizontalVariation(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, z_temp_max

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_AddHorizontalVariation'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        lat_deg = patch_2d%cells%center(jc,jb)%lat
        z_temp_max=0.01_wp* (lat_deg-basin_center_lat) * (lat_deg-basin_center_lat)

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          ocean_temperature(jc,jk,jb) = ocean_temperature(jc,jk,jb) * &
            & EXP(-z_temp_max/basin_height_deg)!(1.0_wp-exp(-z_temp_max/basin_height_deg))

        END DO
      END DO
    END DO

  END SUBROUTINE temperature_AddHorizontalVariation
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_Uniform_SpecialArea(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! 2012-10-31: Indonesian Archipelago - connected to Atlantic? (4 cpu)
    z_lat1  =  -5.0_wp
    z_lat2  =  10.0_wp
    z_lon1  = 115.0_wp
    z_lon2  = 135.0_wp
    ! 2012-10-31: corresponding NAtl: 1N 20W, 3 levels
    z_lon1  = 145.0_wp
    z_lon2  = 160.0_wp
    ! 2012-11-08: Test homogen with one cell differ at 10N; 70E
    z_lat1  =   0.0_wp
    z_lat2  =  15.0_wp
    z_lon1  =  60.0_wp
    z_lon2  =  80.0_wp

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lat_deg = cell_center(jc, jb)%lat * rad2deg
        lon_deg = cell_center(jc, jb)%lon * rad2deg

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          ocean_temperature(jc,jk,jb) = 5.0_wp

          IF ( (lat_deg >= z_lat1 .AND. lat_deg <= z_lat2) .AND. &
            & (lon_deg >= z_lon1 .AND. lon_deg <= z_lon2) .AND. &
            & ( jk <= 1 ) ) THEN

            ocean_temperature(jc,jk,jb) =  6.0_wp

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE temperature_Uniform_SpecialArea
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE salinity_Uniform_SpecialArea(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! 2012-10-31: Indonesian Archipelago - connected to Atlantic? (4 cpu)
    z_lat1  =  -5.0_wp
    z_lat2  =  10.0_wp
    z_lon1  = 115.0_wp
    z_lon2  = 135.0_wp
    ! 2012-10-31: corresponding NAtl: 1N 20W, 3 levels
    z_lon1  = 145.0_wp
    z_lon2  = 160.0_wp
    ! 2012-11-08: Test homogen with one cell differ at 10N; 70E
    z_lat1  =   0.0_wp
    z_lat2  =  15.0_wp
    z_lon1  =  60.0_wp
    z_lon2  =  80.0_wp

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lat_deg = cell_center(jc, jb)%lat * rad2deg
        lon_deg = cell_center(jc, jb)%lon * rad2deg

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          ocean_salinity(jc,jk,jb) = 35.0_wp

          IF ( (lat_deg >= z_lat1 .AND. lat_deg <= z_lat2) .AND. &
             & (lon_deg >= z_lon1 .AND. lon_deg <= z_lon2) .AND. &
             & ( jk <= 1 ) ) THEN

             ocean_salinity(jc,jk,jb) = 34.8_wp

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE salinity_Uniform_SpecialArea
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_APE(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_APE'
    !-------------------------------------------------------------------------
    ! Important:
    !   use initial_temperature_top=27.0 initial_temperature_bottom=0.0
    !   to be consistent with the old setup
    CALL message(TRIM(method_name), ' using sst1')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
         DO jk=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb))
          ocean_temperature(jc,jk,jb) = MAX( &
            & ape_sst(initial_sst_type, patch_2d%cells%center(jc,jb)%lat) - tmelt,  & ! SST in Celsius
            & initial_temperature_bottom)
        END DO
      END DO
    END DO

    CALL decreaseTracerVerticallyLinearly(patch_3d, ocean_tracer=ocean_temperature, &
      & top_value=initial_temperature_top, bottom_value=initial_temperature_bottom)

  END SUBROUTINE temperature_APE
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_VerticallyLinearlyDecreasing(patch_3d, ocean_tracer, top_value, bottom_value)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in) :: top_value, bottom_value

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_VerticallyLinearlyDecreasing'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        DO jk = 1, MIN(1,patch_3d%p_patch_1d(1)%dolic_c(jc,jb))
          ocean_tracer(jc,jk,jb) = top_value
        ENDDO
      END DO
    END DO

    CALL decreaseTracerVerticallyLinearly(patch_3d, ocean_tracer, top_value, bottom_value)

  END SUBROUTINE tracer_VerticallyLinearlyDecreasing
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! decrease tvertically linerarly the given tracer based on the top level value
  ! of the tracer and using a decres of (top_value - bottom_value) / (n_zlev - 1)
  SUBROUTINE decreaseTracerVerticallyLinearly(patch_3d, ocean_tracer, top_value, bottom_value)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in) :: top_value, bottom_value

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: linear_decrease

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    linear_decrease = (top_value - bottom_value) / (REAL(n_zlev,wp)-1.0_wp)

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          ocean_tracer(jc,jk,jb) &
            & = MAX(ocean_tracer(jc,jk-1,jb) - linear_decrease, bottom_value)
        END DO

      END DO
    END DO

  END SUBROUTINE decreaseTracerVerticallyLinearly
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_CollapsingDensityFront_WeakGrad(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg, z_tmp
    REAL(wp):: z_ldiff, z_ltrop, z_lpol
    REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_CollapsingDensityFront_WeakGrad'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' Collapsing density front with weaker gradient')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Temperature profile in first layer depends on latitude only
    ! Construct temperature profile
    !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
    z_ttrop = 10.0_wp      ! tropical temperature
    z_tpol  =  5.0_wp      ! polar temperature
    z_ttrop =  5.0_wp      ! 2011-09-02: instable stratification
    z_tpol  = 10.0_wp      ! 2011-09-02: instable stratification
    z_lpol  = 70.0_wp      ! polar boundary latitude of transition zone

    z_ttrop = 25.0_wp      ! 2011-09-05: stable stratification
    z_tpol  = 10.0_wp      ! 2011-09-05: stable stratification
    z_tdeep =  5.0_wp      ! 2011-09-05: stable stratification
    z_ltrop = 15.0_wp      ! tropical boundary latitude of transition zone
    z_lpol  = 60.0_wp      ! polar boundary latitude of transition zone
    z_tdiff = z_ttrop - z_tpol
    z_ldiff = z_lpol  - z_ltrop

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(jc,jb)%lat * rad2deg

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

              IF(ABS(lat_deg)>=z_lpol)THEN

                ocean_temperature(jc,jk,jb) =  z_tpol

                ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_tpol
                !            ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1)&
                !             &= 30.0_wp!convert_insitu2pot_temp_func(ocean_state%p_diag%temp_insitu(jc,jk,jb),&
                !                     &ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2),&
                !                     &sfc_press_bar)
                !SItodBar*rho_ref*v_base%zlev_m(jk))!1013.0_wp)SItodBar*101300.0_wp)!

              ELSEIF (ABS(lat_deg)<=z_ltrop) THEN

                ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop
                ocean_temperature(jc,jk,jb) = z_ttrop


              ELSE ! IF(ABS(lat_deg)<z_lpol .AND. ABS(lat_deg)>z_ltrop)THEN
                !   z_tmp = pi*((abs(lat_deg) - z_lpol)/z_ldiff)
                !   ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_tpol + 0.5_wp*z_tdiff*(1.0_wp+cos(z_tmp))
                z_tmp = 0.5_wp*pi*((ABS(lat_deg) - z_ltrop)/z_ldiff)
                ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop - z_tdiff*SIN(z_tmp)
                ocean_temperature(jc,jk,jb) = z_ttrop - z_tdiff*SIN(z_tmp)
                !      if (jk==1) write(*,*) 'zlat,ztmp(deg),temp', &
                !   &  jb,jc,lat_deg,(abs(lat_deg)-z_lpol)/z_ldiff,ocean_state%p_diag%temp_insitu(jc,jk,jb)
              ENDIF
!            ELSE
!              ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_tdeep
!              ocean_temperature(jc,jk,jb) = z_tdeep
!            ENDIF  ! jk=1
        END DO
      END DO
    END DO

  END SUBROUTINE temperature_CollapsingDensityFront_WeakGrad
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! Temperature profile depends on latitude and depth
  ! Construct temperature profile
  !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
  !   for maximum tropical temperature see tprof
  SUBROUTINE temperature_TropicsPolar(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg, z_tmp
    REAL(wp):: z_ldiff, z_ltrop, z_lpol
    REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_TropicsPolar'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Temperature profile depends on latitude and depth
    ! Construct temperature profile
    !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
    !   for maximum tropical temperature see values above
    CALL message(TRIM(method_name), ': simple tropics-pol/vertical temperature profile')

    IF (i_sea_ice == 0) THEN
      z_tpol  =  5.0_wp      ! polar temperature
    ELSE
      z_tpol = tf
    ENDIF
    z_ltrop = 15.0_wp      ! tropical latitude for temperature gradient
    z_lpol  = 60.0_wp      ! polar latitude for temperature gradient
    z_ldiff = z_lpol  - z_ltrop

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(jc,jb)%lat * rad2deg

        ! bugfix: z_tpol was 0 for 10 levels since jk was inner loop
        !         does not effect 4 levels
        z_tpols = z_tpol
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

            ! set maximum tropical temperature from profile above
            !IF(n_zlev<=10)THEN
            !  z_ttrop = tprof_var(jk)
            !  z_tpols = MIN(z_tpols,tprof_var(jk))
            !ELSEIF(n_zlev>10.and.n_zlev<=20)THEN
            z_ttrop = tprof(jk)
            z_tpols = MIN(z_tpols,tprof(jk))
            !ENDIF

            IF(ABS(lat_deg)>=z_lpol)THEN

              ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_tpols
              ocean_temperature(jc,jk,jb) = z_tpols

            ELSEIF(ABS(lat_deg)<=z_ltrop)THEN

              ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop
              ocean_temperature(jc,jk,jb) = z_ttrop


            ELSE ! IF(ABS(lat_deg)<z_lpol .AND. ABS(lat_deg)>z_ltrop)THEN

              z_tdiff = z_ttrop - z_tpols
              z_tmp = 0.5_wp*pi*((ABS(lat_deg) - z_ltrop)/z_ldiff)
              ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop - z_tdiff*SIN(z_tmp)
              ocean_temperature(jc,jk,jb) = z_ttrop - z_tdiff * SIN(z_tmp)

            ENDIF

        END DO
      END DO
    END DO

  END SUBROUTINE temperature_TropicsPolar
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_DanilovsMunkGyre(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: z_dolic
    REAL(wp):: distan, lat_deg, lon_deg, z_tmp
    REAL(wp):: perturbation_lat, perturbation_lon,  max_perturbation, perturbation_width, z_ldiff, z_ltrop, z_lpol
    REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_DanilovsMunkGyre'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    !------------------------------
    CALL message(TRIM(method_name), ': Danilovs Munk gyre flow')

    perturbation_lat = basin_center_lat + 0.1_wp * basin_height_deg
    perturbation_lon = basin_center_lon + 0.1_wp * basin_width_deg
    max_perturbation  = 0.1_wp!20.1_wp
    perturbation_width  = 10.0_wp!1.5_wp

    ! Next update 2011-05-24: due to Danilov the perturbation should be -1 Kelvin, width 3.0
    ! 05-25: max and width larger: -2.0 and 5.0
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        z_dolic = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

        IF (z_dolic > 0) THEN
          ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
          ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
          ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
          ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375
          ocean_temperature(jc,1:z_dolic,jb) = 20.0_wp
          distan = SQRT((cell_center(jc,jb)%lat - perturbation_lat * deg2rad)**2 + &
            & (cell_center(jc,jb)%lon - perturbation_lon * deg2rad)**2)

          !Local hot perturbation
          IF(distan<=5.0_wp*deg2rad)THEN
            DO jk = 1, z_dolic
              ocean_temperature(jc,jk,jb) =          &
                & ocean_temperature(jc,jk,jb)          &
                & + max_perturbation*EXP(-(distan/(perturbation_width*deg2rad))**2) &
              !                &   * sin(pi*v_base%zlev_m(jk)/4000.0_wp)!&
                & * SIN(pi*patch_3d%p_patch_1d(1)%zlev_m(jk) / patch_3d%p_patch_1d(1)%zlev_i(z_dolic+1))
            END DO
          ENDIF !Local hot perturbation

        END IF !(z_dolic > 0)
      END DO
    END DO

  END SUBROUTINE temperature_DanilovsMunkGyre
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_BasinWithVerticalWall(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_BasinWithVerticalWall'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Adjusting density front in a basin: vertical wall at basin_center_lon
    CALL message(TRIM(method_name),' Adjusting density front in a basin with vertical wall')

    !Impose temperature profile. Profile
    !depends on latitude only and is uniform across
    !all vertical layers
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        IF (cell_center(jc,jb)%lon >= basin_center_lon) THEN
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,1:n_zlev,jb) = 30.0_wp
          ENDDO
        ELSE
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            ocean_temperature(jc,1:n_zlev,jb) = 25.0_wp
          ENDDO
        ENDIF
      END DO
    END DO

  END SUBROUTINE temperature_BasinWithVerticalWall
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_CollapsingDensityFront_StuhnePeltier(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: z_dolic
    REAL(wp):: lat_deg, lon_deg, z_tmp
    REAL(wp):: perturbation_lat, perturbation_lon,  max_perturbation, perturbation_width, z_ldiff, z_ltrop, z_lpol
    REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_CollapsingDensityFront_StuhnePeltier'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(TRIM(method_name), ': Collapsing density front, Stuhne-Peltier')

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(jc,jb)%lat * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across
        !all vertical layers
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          IF (ABS(lat_deg) >= 40.0_wp) THEN

            ocean_temperature(jc,jk,jb) = 5.0_wp

            !ocean_state%p_diag%temp_insitu(jc,jk,jb) = 5.0_wp
            !convert_insitu2pot_temp_func(ocean_state%p_diag%temp_insitu(jc,jk,jb),&
            !                      &ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2),&
            !                      &sfc_press_bar)
            !SItodBar*rho_ref*v_base%zlev_m(jk))!1013.0_wp)

          ELSEIF (ABS(lat_deg) <= 20.0_wp) THEN

            ocean_temperature(jc,jk,jb) =  30.0_wp

            ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = 30.0_wp
            !convert_insitu2pot_temp_func(ocean_state%p_diag%temp_insitu(jc,jk,jb),&
            !                     &ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2),&
            !                     &sfc_press_bar)
            !SItodBar*rho_ref*v_base%zlev_m(jk))!1013.0_wp)SItodBar*101300.0_wp)!


          ELSE ! IF (ABS(lat_deg) < 40.0_wp .AND. ABS(lat_deg) > 20.0_wp)THEN

            z_tmp = pi*((ABS(lat_deg) -20.0_wp)/20.0_wp)

            ocean_temperature(jc,jk,jb) = &
              & 5.0_wp + 0.5_wp * 25.0_wp * (1.0_wp + COS(z_tmp))

            ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = 5.0_wp&
            !  & + 0.5_wp*25.0_wp*(1.0_wp+COS(z_tmp))
            !  ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1)= &
            !      convert_insitu2pot_temp_func(ocean_state%p_diag%temp_insitu(jc,jk,jb),&
            !                                                &ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,2),&
            !                                                &sfc_press_bar)
            !SItodBar*rho_ref*v_base%zlev_m(jk))!1013.0_wp)SItodBar*101300.0_wp)!
          ENDIF

        END DO
      END DO
    END DO

   END SUBROUTINE temperature_CollapsingDensityFront_StuhnePeltier
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE salinity_AnalyticSmoothVerticalProfile(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    REAL(wp) :: salinity_profile(n_zlev)
    INTEGER :: jk

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_AnalyticSmoothVerticalProfile'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name),' Creating analytic profile:')

    DO jk=1,n_zlev

      salinity_profile(jk) = &
        & MIN(34.1_wp + LOG(1.3_wp + SQRT(patch_3d%p_patch_1d(1)%zlev_m(jk)) * 0.05), 35.0_wp)

      ! write(0,*) jk, patch_3D%p_patch_1D(1)%zlev_m(jk), " salinity:", salinity_profile(jk)
      WRITE(message_text,*) jk, patch_3d%p_patch_1d(1)%zlev_m(jk), " salinity:", salinity_profile(jk)
      CALL message(TRIM(method_name),TRIM(message_text))
    ENDDO

    CALL fill_FromVerticalArrayProfile(patch_3d=patch_3d, ocean_tracer=ocean_salinity, VerticalProfileValue=salinity_profile)

  END SUBROUTINE salinity_AnalyticSmoothVerticalProfile
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE fill_FromVerticalArrayProfile(patch_3d, ocean_tracer, VerticalProfileValue)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), TARGET :: VerticalProfileValue(:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index

    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !------------------------------
    ! assign from adhoc array values
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          ocean_tracer(jc,jk,jb) = VerticalProfileValue(jk)
        END DO
      END DO
    END DO

  END SUBROUTINE fill_FromVerticalArrayProfile
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
!  FUNCTION tropical_polar_temperature()
!              ! set maximum tropical temperature from profile above
!              !IF(n_zlev<=10)THEN
!              !  z_ttrop = tprof_var(jk)
!              !  z_tpols = MIN(z_tpols,tprof_var(jk))
!              !ELSEIF(n_zlev>10.and.n_zlev<=20)THEN
!              z_ttrop = tprof(jk)
!              z_tpols = MIN(z_tpols,tprof(jk))
!              !ENDIF
!
!              IF(ABS(lat_deg)>=z_lpol)THEN
!
!                ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_tpols
!                ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = z_tpols
!
!              ELSEIF(ABS(lat_deg)<=z_ltrop)THEN
!
!                ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop
!                ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = z_ttrop
!
!
!              ELSE ! IF(ABS(lat_deg)<z_lpol .AND. ABS(lat_deg)>z_ltrop)THEN
!
!                z_tdiff = z_ttrop - z_tpols
!                z_tmp = 0.5_wp*pi*((ABS(lat_deg) - z_ltrop)/z_ldiff)
!                ! ocean_state%p_diag%temp_insitu(jc,jk,jb) = z_ttrop - z_tdiff*SIN(z_tmp)
!                ocean_state%p_prog(nold(1))%tracer(jc,jk,jb,1) = z_ttrop - z_tdiff * SIN(z_tmp)
!
!              ENDIF
!  END FUNCTION tropical_polar_temperature
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_surface_height(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: z_dolic
    REAL(wp):: distan, lat_deg, lon_deg, z_tmp
    REAL(wp):: perturbation_lat, perturbation_lon

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_surface_height'
    !-------------------------------------------------------------------------

    IF (sea_surface_height_type < 200) RETURN ! not analytic sea height

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    SELECT CASE (sea_surface_height_type)
    CASE (200)
      ! 0 height, this is the initialization value,
      ! so no need to explicilty define this case
      ! the whole grid is considered sea
      ocean_height(:,:) = 0.0_wp

    CASE (201)
      ! #slo#: simple elevation between 30W and 30E (pi/3.)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          IF ( patch_3d%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            ocean_height(jc,jb) = 10.0_wp * &
              & SIN(cell_center(jc, jb)%lon * 6.0_wp) * COS(cell_center(jc, jb)%lat * 3.0_wp)

          ENDIF
        END DO
      END DO

    CASE (202)
      ! Add elevation perturbation at new values - 35N; 10W
      ! not clear yet
      perturbation_lat = basin_center_lat + 0.1_wp * basin_height_deg
      perturbation_lon = basin_center_lon + 0.1_wp * basin_width_deg
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index

          IF (patch_3d%p_patch_1d(1)%dolic_c(jc,jb) > 0) THEN

            distan=SQRT((cell_center(jc, jb)%lat - perturbation_lat * deg2rad)**2 + &
              & (cell_center(jc, jb)%lon - perturbation_lon * deg2rad)**2)
            !IF(distan<=15.5_wp*deg2rad) cycle
            IF(distan < 10.0_wp * deg2rad) THEN
              ocean_height(jc,jb) = 0.5_wp & !ocean_state%p_prog(nold(1))%h(jc,jb)&
                & + 0.3_wp * EXP(-(distan/(2.2_wp*deg2rad))**2)
            ENDIF

          ENDIF

        END DO
      END DO

    CASE (203)
      ! test_usbr_h
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          ocean_height(jc,jb)  = test_usbr_h( cell_center(jc, jb)%lon, cell_center(jc, jb)%lat, 0.0_wp)
          ! write(*,*)'h orig, bathy_c:', cell_center(jc, jb)%lon, cell_center(jc, jb)%lat,ocean_state%p_prog(nold(1))%h(jc,jb)!
        END DO
      END DO

    CASE (204)
      ! test2_h
      CALL message(TRIM(method_name), ' h for Williamson Test 2')
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          ocean_height(jc,jb) = test2_h( cell_center(jc, jb)%lon, cell_center(jc, jb)%lat, 0.0_wp)
        END DO
      END DO

    CASE (205)
      ! test5_h
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          ocean_height(jc,jb) = test5_h( cell_center(jc, jb)%lon, cell_center(jc, jb)%lat, 0.0_wp)
        END DO
      END DO


    CASE default
      CALL finish(method_name, "unknown sea_surface_height_type")

    END SELECT

    CALL dbg_print('init_ocean_surface_height', ocean_height, str_module,  1, &
        & in_subset=patch_2d%cells%owned)

  END SUBROUTINE init_ocean_surface_height
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_bathymetry(patch_3d, cells_bathymetry)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET  :: cells_bathymetry(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    ! Local Variables
    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: z_dolic
    REAL(wp):: lat_deg, lon_deg, z_tmp

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_bathymetry'
    !-------------------------------------------------------------------------

    IF (topography_type < 200) RETURN ! not analytic bathymetry

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    SELECT CASE (topography_type)
    CASE (200)
      ! constant depth given by topography_height_reference
      ! the whole grid is considered sea
      patch_3d%lsm_c(:,:,:) = sea
      patch_3d%lsm_e(:,:,:) = sea
      cells_bathymetry(:,:) = topography_height_reference

    CASE (201)
      CALL mountain_orography_Williamson_test5(patch_3d, cells_bathymetry)

    CASE default
      CALL finish(method_name, "unknown topography_type")
    END SELECT

  END SUBROUTINE init_ocean_bathymetry
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  FUNCTION geo_balance_mim(patch_2d, h_e,grad_coeff, rhs_e) result(vn_e)
    !
    TYPE(t_patch) :: patch_2d
    REAL(wp)      :: h_e(:,:)
    REAL(wp)      :: grad_coeff(:,:,:)
    REAL(wp)      :: rhs_e(:,:)!(nproma,n_zlev,patch_2D%nblks_e)
    REAL(wp)      :: vn_e(SIZE(rhs_e,1),SIZE(rhs_e,2))
    !
    !LOCAL VARIABLES
    ! INTEGER,PARAMETER :: nmax_iter= 200 ! maximum number of iterations
    REAL(wp) :: zimpl_coeff = 1.0_wp    !COEFF has to be set appropriately !!!!
    REAL(wp) :: zimpl_prime_coeff
    ! INTEGER  :: n_iter =0               ! number of iterations
    REAL(wp) :: tolerance               ! (relative or absolute) tolerance
    ! REAL(wp) :: z_residual(nmax_iter)
    LOGICAL :: lmax_iter               ! true if reached m iterations
    ! LOGICAL  :: lverbose = .TRUE.
    ! INTEGER  :: jk
    REAL(wp) :: rhstemp(nproma,patch_2d%nblks_e)
    REAL(wp), ALLOCATABLE :: vn_e2(:,:)!(nproma,patch_2D%nblks_e)
    ! INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    ! INTEGER :: rl_start_e, rl_end_e, je,jb
    
    !-----------------------------------------------------------------------
    tolerance         = 1.0e-13_wp  ! solver_tolerance
    zimpl_prime_coeff = (1.0_wp-zimpl_coeff)
    
    ALLOCATE (vn_e2(nproma,patch_2d%nblks_e))
    vn_e2(:,:)   = 0.0_wp
    rhstemp(:,:) = 0.0_wp
    
    
    vn_e2(:,:) = rhs_e(:,:)
    !rhstemp(:,:) = rhs_e(:,:)&
    !  & -zimpl_prime_coeff*lhs_geo_balance_mim(vn_e2,patch_2D,jk,zimpl_coeff,grad_coeff,h_e)
    
    IF (MAXVAL (ABS (rhstemp (:,:))) <= tolerance) THEN
      vn_e(:,:) = vn_e2(:,:)
      PRINT*, "Inv_geo balance GMRES solved by initial guess!",&
        & MAXVAL(ABS(rhstemp(:,:))), MAXVAL(ABS(rhs_e(:,:)))
    ELSE
      vn_e2 = 0.0_wp!rhs_e(:,jk,:)
      
      !  rhstemp(:,:) = rhs_e(:,:)-lhs_geo_balance_mim(vn_e2(:,:),patch_2D, jk,&
      !    &            zimpl_coeff,grad_coeff, h_e)
      ! WRITE(*,*)'max/min residual of inverse primal-flip-flop:',&
      !&jk, maxval(rhstemp),minval(rhstemp)
      
      IF (MAXVAL (ABS (rhstemp (:,:))) >= tolerance) lmax_iter = .TRUE.
      !          IF (lverbose) THEN
      !            IF (lmax_iter) THEN
      !              WRITE (6, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2, 1x, a)') &
      !              &'Inv_flipflop GMRES #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
      !              &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:))), 'GMRES PROBLEM!!!!!!!!!!!!'
      !            ELSE
      !              WRITE (6, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2)') &
      !              &'Inv_flipflop GMRES #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
      !              &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:)))
      !            ENDIF
      !        ENDIF
      vn_e(:,:) = vn_e2(:,:)
    END IF
    
    ! DO jk=1, 1
    !   DO jb = i_startblk_e, i_endblk_e
    !     CALL get_indices_e(patch_2D, jb,&
    !                      & i_startblk_e, i_endblk_e,&
    !                      & i_startidx_e, i_endidx_e,&
    !                      & rl_start_e, rl_end_e)
    !     DO je =  i_startidx_e, i_endidx_e
    !       IF(rhs_e(je,jk,jb)/=0.0_wp)THEN
    !       write(*,*)'RHS:solution:', jk,je,jb,rhs_e(je,jk,jb), inv_flip_flop_e(je,jk,jb)
    !       ENDIF
    !     END DO
    !   END DO
    ! END DO
    
    
    DEALLOCATE (vn_e2)
    
  END FUNCTION geo_balance_mim
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  FUNCTION lhs_geo_balance_mim( x, patch_2d, patch_3d, lev,p_coeff,grad_coeff, h_e) result(llhs)
    TYPE(t_patch),TARGET,INTENT(in) :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    INTEGER :: lev
    REAL(wp),INTENT(inout)          :: x(nproma,patch_2d%nblks_e)!(:,:)
    REAL(wp),INTENT(in)             :: p_coeff
    REAL(wp), INTENT(in)            :: grad_coeff(:,:,:)
    REAL(wp),OPTIONAL,INTENT(in)    :: h_e(nproma,patch_2d%nblks_e)!(SIZE(x,1), SIZE(x,2))!(:,:)
    REAL(wp)                        :: llhs(nproma,patch_2d%nblks_e)!(SIZE(x,1), SIZE(x,2))
    
    !locl variables
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: jc,jb
    REAL(wp) :: z_x_e(nproma,patch_2d%nblks_e)
    REAL(wp) :: z_x_vort(nproma,1,patch_2d%nblks_v)
    REAL(wp) :: z_x_out(SIZE(x,1), 1,SIZE(x,2))!(nproma,patch_2D%nblks_e)
    !REAL(wp) :: z_vt(SIZE(x,1), 1,SIZE(x,2))!(nproma,patch_2D%nblks_e)
    REAL(wp) :: z_grad(SIZE(x,1), 1,SIZE(x,2))!(nproma,patch_2D%nblks_e)
    REAL(wp) :: z_kin(nproma,1,patch_2d%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates)    :: z_pv_cc(nproma,patch_2d%alloc_cell_blocks)
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    all_cells => patch_2d%cells%ALL
    
    z_x_vort(:,:,:)= 0.0_wp
    z_x_out(:,:,:) = 0.0_wp
    z_x_e(:,:)   = x(:,:)
    WRITE(*,*)'warning: edge2cell mapping missing'
    STOP
    !       CALL map_edges2cell_3D( patch_2D, &
    !                             & z_x_e,   &
    !                             &  z_pv_cc,&
    !                             & p_coeff, &
    !                             & level=1)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc =  start_cell_index, end_cell_index
        z_kin(jc,1,jb) = 0.5_wp*DOT_PRODUCT(z_pv_cc(jc,jb)%x,z_pv_cc(jc,jb)%x)
      END DO
    END DO
    
    CALL grad_fd_norm_oce_3d( z_kin,  &
      & patch_3d,  &
      & grad_coeff,&
      & z_grad)
    
    !   z_x_out(:,:,:) = dual_flip_flop(patch_2D, z_x_e, z_x_e, z_x_vort, h_e,&
    !                                  &opt_slev=1, opt_elev=1)
    
    z_x_out=z_x_out!+z_grad
    llhs(1:nproma,1:patch_2d%nblks_e) = p_coeff*z_x_out(1:nproma,1,1:patch_2d%nblks_e)
    !write(*,*)'max/min LHS', maxval(llhs(:,:)),minval(llhs(:,:))
    
  END FUNCTION lhs_geo_balance_mim
  !--------------------------------------------------------------------

  
  !--------------------------------------------------------------------
  !--Below are functions from to implement tests from Williamson shallow-water tests
  !-------------------------------------------------------------------------
  !
  ! !F*UNCTION INTERFACE:
  FUNCTION test0_h( point_lon, point_lat, p_t) result( p_hh)
    !
    ! !DESCRIPTION:
    ! Initial datum for height, test case 0 (conical mountain). \\
    ! Not included in Williamson et al. (1992)
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Revised to programming guide by Th.Heinze, DWD, (2006-12)
    !
    ! !DEFINED PARAMETERS:
    REAL(wp), PARAMETER :: h0=2000._wp  ! basic height level
    REAL(wp), PARAMETER :: h1=1000._wp  ! max height of conical mountain
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! geopotential height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_r1      ! distance point to center of con. mount.
    REAL(wp)             :: z_r2      ! radius of conical mountain
    REAL(wp)             :: z_lon2    ! longitude of center of conical mount.
    REAL(wp)             :: z_lat2    ! latitude of center of conical mount.
    REAL(wp)             :: z_dlon    ! longitudinal distance
    REAL(wp)             :: z_dlat    ! latitudinal distance
    
    ! center and radius of conical mountain
    
    z_lon2  = 1.5_wp * pi
    z_lat2  = pi / 6._wp
    z_r2    = pi / 9._wp
    
    ! distance of point to center (not great arc distance!)
    
    z_dlon = point_lon - z_lon2
    z_dlat = point_lat - z_lat2
    
    z_dlon = z_dlon * z_dlon
    z_dlat = z_dlat * z_dlat
    
    z_r1 = z_dlon + z_dlat
    z_r1 = MIN( z_r2 * z_r2, z_r1)
    z_r1 = SQRT(z_r1)
    
    ! geopotential height
    
    IF( z_r1 < z_r2) THEN              ! point within radius
      
      p_hh = 1._wp - z_r1 / z_r2
      p_hh = h0 + h1 * p_hh
      
    ELSE                               ! point outside of radius
      
      p_hh = h0
      
    ENDIF
    
  END FUNCTION test0_h
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  FUNCTION test2_h( point_lon, point_lat, p_t) result(p_hh)
    !
    ! !DESCRIPTION:
    ! Initial datum for height, test case 2 of Williamson et al.(1992).
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Revised to programming guide by Th.Heinze, DWD, (2006-12)
    !
    ! !DEFINED PARAMETERS:
    REAL(wp), PARAMETER :: h0 = 2.94e4_wp * rgrav  ! maximum height
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_fact1   ! 1st factor
    REAL(wp)             :: z_fact2   ! 2nd factor
    
    ! 1st factor
    
    z_fact1 = sphere_radius * grid_angular_velocity
    z_fact1 = z_fact1 + 0.5_wp * u0
    z_fact1 = z_fact1 * u0 * rgrav
    
    ! 2nd factor
    
    z_fact2 = SIN(point_lat) * COS(aleph)
    z_fact2 = z_fact2 - COS(point_lon) * COS(point_lat) * SIN(aleph)
    z_fact2 = z_fact2 * z_fact2
    
    ! height
    
    p_hh = h0 - z_fact1 * z_fact2
    
  END FUNCTION test2_h
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  FUNCTION test2_vort( point_lon, point_lat, p_t) result(p_vort)
    !
    ! !DESCRIPTION:
    ! Initial datum for relative vorticity, test case 2 of Williamson et al.(1992).
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Revised to programming guide by Th.Heinze, DWD, (2006-12)
    !
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_vort    ! relative vorticity

    p_vort = SIN(point_lat)* COS(aleph)
    p_vort = p_vort - COS(point_lon) * COS(point_lat) * SIN(aleph)
    p_vort = 2._wp * u0 * sphere_radius * p_vort
    
  END FUNCTION test2_vort
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  FUNCTION test5_h( point_lon, point_lat, p_t) result(p_hh)
    !
    ! !DESCRIPTION:
    ! Initial datum for height, test case 5 of Williamson et al.(1992).
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Revised to programming guide by Th.Heinze, DWD, (2007-01)
    !
    ! !DEFINED PARAMETERS:
    REAL(wp), PARAMETER :: h0    = 5960._wp  ! maximum height
    REAL(wp), PARAMETER :: uzero = 20._wp    ! maximum velocity
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_fact1   ! 1st factor
    REAL(wp)             :: z_fact2   ! 2nd factor
    

    ! 1st factor
    
    z_fact1 = sphere_radius * grid_angular_velocity
    z_fact1 = z_fact1 + 0.5_wp * uzero
    z_fact1 = z_fact1 * uzero * rgrav
    
    ! 2nd factor
    
    z_fact2 = SIN(point_lat)
    z_fact2 = z_fact2 * z_fact2
    
    ! height
    
    p_hh = h0 - z_fact1 * z_fact2
    
  END FUNCTION test5_h
  !-------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  !
  ! Initial datum for orography, test case 5 of Williamson et al.(1992).
  !
  ! !REVISION HISTORY:
  ! Developed  by L.Bonaventura  (2002-5).
  ! Revised to programming guide by Th.Heinze, DWD, (2007-02)
  !
  SUBROUTINE mountain_orography_Williamson_test5(patch_3d, cells_bathymetry)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET, INTENT(inout)  :: cells_bathymetry(:,:)

    REAL(wp), PARAMETER :: h_s0  = 2000._wp  ! maximum height of mountain
    
    REAL(wp) :: point_lon     ! longitude of point
    REAL(wp) :: point_lat     ! latitude of point
    REAL(wp) :: p_t           ! point of time
    
    REAL(wp)             :: point_height      ! orography
    
    REAL(wp)             :: z_lon_mc  ! Mountain center, longitude ...
    REAL(wp)             :: z_lat_mc  !          ... and latitude
    REAL(wp)             :: z_rad_mt  ! radius of mountain
    REAL(wp)             :: z_dist_mc ! distance from mountain center
    REAL(wp)             :: z_diff    ! difference of coordinates
    REAL(wp)             :: z_min_dist_sq ! min of square of distances

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    ! Local Variables
    INTEGER :: jb, jc, je, jk
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':mountain_orography_Williamson_test5'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    ! center and radius of mountain
    z_lon_mc = -pi_2
    z_lat_mc = pi / 6._wp
    z_rad_mt = pi / 9._wp
    
    patch_3d%lsm_c(:,:,:) = sea
    patch_3d%lsm_e(:,:,:) = sea

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        point_lon = patch_2d%cells%center(jc, jb)%lon
        point_lat = patch_2d%cells%center(jc, jb)%lat
        p_t       = 0.0_wp

        ! square of distance (in geographical coordinate sense) of point
        ! from mountain center

        z_diff = point_lon - z_lon_mc
        z_diff = z_diff * z_diff
        z_dist_mc = z_diff

        z_diff = point_lat - z_lat_mc
        z_diff = z_diff * z_diff
        z_dist_mc = z_dist_mc + z_diff

        ! if point inside mountain range take its distance, else take mountain radius

        z_diff = z_rad_mt * z_rad_mt
        z_min_dist_sq = MIN ( z_diff, z_dist_mc)
        z_dist_mc = SQRT( z_min_dist_sq)

        ! conical shape of mountain, depending on distance from mountain center

        point_height = z_dist_mc / z_rad_mt
        point_height = 1._wp - point_height

        cells_bathymetry(jc, jb) = h_s0 * point_height

      ENDDO
    ENDDO
    
  END SUBROUTINE mountain_orography_Williamson_test5
  !-----------------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  FUNCTION test6_h(point_lon, point_lat, p_t) result(p_hh)
    !
    ! !DESCRIPTION:
    !
    ! Initial datum for geopotential h, test case 6 of Williamson et al.(1992).
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Modified by Th.Heinze, DWD, (2006-11-02)
    
    ! !DEFINED PARAMETERS:
    !    REAL (wp), PARAMETER  :: h0 = 8000._wp, re_omg_kk = 50._wp
    REAL (wp), PARAMETER :: h0 = 8000._wp, omg_kk = 7.848e-6_wp !(sphere_radius * omg_kk is not 50.)
    ! pripodas
    INTEGER,   PARAMETER :: r = 4
    
    ! !INPUT PARAMETERS:
    REAL(wp) , INTENT(in) :: point_lon, point_lat, p_t
    
    ! !RETURN VALUE:
    REAL(wp)              :: p_hh
    
    ! !LOCAL VARIABLES:
    ! REAL(wp)              :: z_omg, z_phia, z_phib, z_phic, z_r_omega
    REAL(wp)              :: z_phia, z_phib, z_phic, z_r_omega , z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_cosfi2, z_cosfir, z_cosfir2, z_cosfir2m2
    REAL(wp)              :: z_cosdl, z_cosd2l, z_dlon, z_rr1r2
    REAL(wp)              :: z_val
    
    !INTEGER               :: i_r1, i_r1r1, i_r1r2, i_r2, j
    INTEGER :: j
    REAL(wp)              :: z_r, z_r1, z_r1r1, z_r1r2, z_r2
    !-----------------------------------------------------------------------
    
    z_r_omega  = sphere_radius * grid_angular_velocity
    !z_omg     = re_omg_kk / sphere_radius
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
    !    i_r1      = r + 1
    !    i_r2      = r + 2
    !    i_r1r1    = i_r1 * i_r1
    !    i_r1r2    = i_r1 * i_r2
    !    z_rr1r2   = 1._wp / i_r1r2
    
    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r1    = z_r1 * z_r1
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2
    
    z_dlon    = omg_kk * z_r * (3._wp+z_r) - 2.0_wp * grid_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (point_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)
    z_cosd2l  = COS(2._wp * z_dlon)
    
    z_cosfi   = COS(point_lat)
    z_cosfi2  = z_cosfi  * z_cosfi    ! cos^2(lat)
    
    z_cosfir  = z_cosfi

    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO

    z_cosfir2m2 = z_cosfir
    z_cosfir2m2 = z_cosfir2m2 * z_cosfir2m2   ! cos^{2*r1-2}(lat)
    
    z_cosfir  = z_cosfir * z_cosfi    ! cos^{r1}(lat)
    z_cosfir2 = z_cosfir * z_cosfir   ! cos^{2*r1}(lat)
    
    z_val  = -.25_wp + z_r
    z_val  = 2._wp * z_val * z_val - 2.125_wp   ! 2r^2 - r -2
    
    z_phia = z_val * z_cosfi2
    
    z_val  = 2._wp * REAL(r,wp) * z_r
    
    z_phia = z_phia - z_val
    z_val  = z_cosfi2 * z_cosfi2 * z_r1
    z_phia = z_phia + z_val
    
    z_phia = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2m2 * z_phia
    z_val  = .5_wp * z_re_omg_kk * (2._wp * z_r_omega + z_re_omg_kk) * z_cosfi2
    z_phia = z_val + z_phia
    
    z_phib = -1._wp * z_cosfi2 * z_r1r1 + z_r1r1 + 1._wp
    z_phib = z_re_omg_kk * (z_r_omega + z_re_omg_kk) * z_cosfir * z_phib
    z_phib = 2._wp * z_rr1r2 * z_phib
    
    z_phic = z_r1 * z_cosfi2 - 1._wp * z_r2
    z_phic = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2 * z_phic
    
    p_hh   = (z_phia + z_phib * z_cosdl + z_phic * z_cosd2l) * rgrav
    p_hh   = h0 + p_hh
    
  END FUNCTION test6_h
  !-----------------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  FUNCTION test6_u( point_lon, point_lat, p_t) result(p_uu)
    !
    ! !DESCRIPTION:
    !
    ! Initial datum for zonal velocity u, test case 6 of Williamson et al.(1992) .
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Modified by Th.Heinze, DWD, (2006-11-02)
    
    ! !DEFINED PARAMETERS:
    ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER ::  omg_kk = 7.848e-6_wp !(sphere_radius * omg_kk is not 50.)
    ! pripodas
    INTEGER,   PARAMETER :: r = 4
    
    ! !INPUT PARAMETERS:
    REAL(wp) , INTENT(in) :: point_lon, point_lat, p_t
    
    ! !RETURN VALUE:
    REAL(wp)              :: p_uu
    
    ! !LOCAL VARIABLES:
    !REAL(wp)              :: z_omg, z_r_omega
    REAL(wp)              :: z_r_omega, z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_cosfi2, z_sinfi, z_sinfi2
    REAL(wp)              :: z_cosfir, z_cosfirm1
    REAL(wp)              :: z_cosdl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val
    
    !    INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2    !pripodas, better transform to real values
    !-----------------------------------------------------------------------

    z_r_omega  = sphere_radius * grid_angular_velocity
    !z_omg     = re_omg_kk / sphere_radius
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
    
    
    !    i_r1      = r + 1
    !    i_r2      = r + 2
    !    i_r1r2    = i_r1 * i_r2
    !    z_rr1r2   = 1._wp / i_r1r2
    
    z_r       = REAL(r, wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2
    
    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * grid_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (point_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)
    
    z_sinfi   = SIN(point_lat)
    z_sinfi2  = z_sinfi * z_sinfi
    
    z_cosfi   = COS(point_lat)
    z_cosfi2  = z_cosfi * z_cosfi
    
    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfirm1 = z_cosfir
    
    z_val      = z_r * z_sinfi2 - z_cosfi2
    z_val      = z_cosfirm1 * z_val * z_cosdl
    z_val      = z_cosfi + z_val
    p_uu       = z_re_omg_kk * z_val
    
  END FUNCTION test6_u
  !-----------------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  FUNCTION test6_v( point_lon, point_lat, p_t) result(p_vv)
    !
    ! !DESCRIPTION:
    !
    ! Initial datum for meridional velocity v, test case 6 of Williamson
    ! et al.(1992).
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Modified by Th.Heinze, DWD, (2006-11-02)
    
    ! !DEFINED PARAMETERS:
    ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER ::  omg_kk = 7.848e-6_wp !(sphere_radius * omg_kk is not 50.)
    ! pripodas
    INTEGER,   PARAMETER :: r = 4
    
    ! !INPUT PARAMETERS:
    REAL(wp) , INTENT(in) :: point_lon, point_lat, p_t
    
    ! !RETURN VALUE:
    REAL(wp)              :: p_vv
    
    ! !LOCAL VARIABLES:
    !REAL(wp)              :: z_omg, z_r_omega   !pripodas, we use omg_kk and not re_omg_kk
    REAL(wp)              :: z_r_omega, z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_sinfi
    REAL(wp)              :: z_cosfir, z_cosfirm1
    REAL(wp)              :: z_sindl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val
    
    !    INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2    !pripodas, better transform to real values
    !-----------------------------------------------------------------------
    
    z_r_omega = sphere_radius * grid_angular_velocity
    !z_omg     = re_omg_kk / sphere_radius
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
    !    i_r1      = r + 1
    !    i_r2      = r + 2
    !    i_r1r2    = i_r1 * i_r2
    !    z_rr1r2   = 1._wp / i_r1r2
    
    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2
    
    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * grid_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (point_lon - z_dlon) * z_r
    z_sindl   = SIN(z_dlon)
    
    z_sinfi   = SIN(point_lat)
    
    z_cosfi   = COS(point_lat)
    
    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfirm1 = z_cosfir
    
    z_val      = z_cosfirm1 * z_sinfi * z_sindl
    p_vv       = -1._wp * z_re_omg_kk * z_r * z_val
    
  END FUNCTION test6_v
  !-----------------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  FUNCTION test6_vort( point_lon, point_lat, p_t) result(p_vt)
    !
    ! !DESCRIPTION:
    !
    ! Initial datum for relative vorticity, test case 6 of Williamson et al.(1992).
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Modified by Th.Heinze, DWD, (2006-11-02):
    ! - corrected vorticity
    
    ! !DEFINED PARAMETERS:
    ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER ::  omg_kk = 7.848e-6_wp !(sphere_radius * omg_kk is not 50.)
    ! pripodas
    INTEGER,   PARAMETER :: r = 4
    
    ! !INPUT PARAMETERS:
    REAL(wp) , INTENT(in) :: point_lon, point_lat, p_t
    
    ! !RETURN VALUE:
    REAL(wp)              :: p_vt
    
    ! !LOCAL VARIABLES:
    !REAL(wp)              :: z_omg, z_r_omega   !pripodas, we use omg_kk and not re_omg_kk
    REAL(wp)              :: z_r_omega, z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_sinfi
    REAL(wp)              :: z_cosfir
    REAL(wp)              :: z_cosdl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val
    
    !INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2   !pripodas, better transform to real values
    !-----------------------------------------------------------------------
    
    z_r_omega = sphere_radius * grid_angular_velocity
    !z_omg     = re_omg_kk / sphere_radius
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
    ! i_r1      = r + 1
    ! i_r2      = r + 2
    ! i_r1r2    = i_r1 * i_r2
    ! z_rr1r2   = 1._wp / i_r1r2
    
    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2
    
    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * grid_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (point_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)
    
    z_sinfi   = SIN(point_lat)
    
    z_cosfi   = COS(point_lat)
    
    z_cosfir  = z_cosfi
    DO j= 2, r
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    
    z_val     = z_cosfir * z_r1 * z_r2 * z_cosdl
    z_val     = 2._wp - z_val
    p_vt      = omg_kk * z_sinfi * z_val
    
  END FUNCTION test6_vort
  !-----------------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  FUNCTION test_usbr_h( point_lon, point_lat, p_t) result(p_hh)
    !
    ! !DESCRIPTION:
    ! Initial datum for height h, test case unsteady solid body
    ! rotation of L\"auter et al.(2007).
    
    ! !REVISION HISTORY:
    ! Developed by Th.Heinze, DWD, (2007-03)
    !
    ! !DEFINED PARAMETERS:
    REAL(wp), PARAMETER :: d0    = 133681.0_wp  ! additive constant
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_phi_t_k ! 1st summand
    REAL(wp)             :: z_summand ! 2nd summand
    REAL(wp)             :: z_fact    ! factor
    REAL(wp)             :: z_angle1  ! 1st angle
    REAL(wp)             :: z_angle2  ! 2nd angle
    !-----------------------------------------------------------------------
    
    z_angle1 = .25_wp * pi
    z_angle2 = point_lon + grid_angular_velocity * p_t
    
    ! 1st summand: \phi_t(\vec c) \cdot \vec k
    
    z_phi_t_k = SIN(point_lat) * COS(z_angle1)
    z_phi_t_k = z_phi_t_k - COS(z_angle2) * COS(point_lat) * SIN(z_angle1)
    z_phi_t_k = u0 * z_phi_t_k
    
    ! 2nd summand: r_e \grid_angular_velocity \sin \varphi
    
    z_summand = sphere_radius * grid_angular_velocity * SIN(point_lat)
    
    ! one factor
    
    z_fact    = .5_wp *  z_phi_t_k + z_summand
    
    ! height
    
    p_hh      = d0 - z_phi_t_k *  z_fact
    p_hh      = p_hh * rgrav
    !write(*,*)'param:', u0, pi, rgrav,sphere_radius, grid_angular_velocity
    !stop
  END FUNCTION test_usbr_h
  !-------------------------------------------------------------------------
    

  !-------------------------------------------------------------------------
  FUNCTION test_usbr_oro(point_lon,point_lat,p_t) result(point_height)
    !
    ! !DESCRIPTION:
    ! Initial datum for orography, test case unsteady solid body rotation
    ! of L\"auter et al.(2007).
    !
    ! !REVISION HISTORY:
    ! Developed by Th.Heinze, DWD, (2007-03)
    !
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: point_height      ! orography
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_fact    ! factor
    !-----------------------------------------------------------------------

    ! calculate factor
    
    z_fact = sphere_radius * grid_angular_velocity * SIN(point_lat)
    z_fact = z_fact * z_fact
    
    ! height of orography
    
    point_height = .5_wp * z_fact * rgrav
    
  END FUNCTION test_usbr_oro
  !-----------------------------------------------------------------------------------
  
  

  !-------------------------------------------------------------------------
  SUBROUTINE rotate(point_lon, point_lat, p_alpha, p_rotlon, p_rotlat)
    !
    ! !DESCRIPTION:
    ! This subroutine computes the rotated coordinates p\_rotlon, p\_rotlat
    ! for a roatation by angle p\_alpha, given the coordinates p\_lon and p\_lat.
    !
    ! !REVISION HISTORY:
    ! Developed originally by R.Jakob for NCAR shallow water model.
    ! Adapted to ICON code by L.Bonaventura (2002-5).
    ! Adapted to ICON programming guide by Th.Heinze, DWD, (2006-12-12)
    !
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in)  :: point_lon     ! ORIGINAL LONGITUDE
    REAL(wp), INTENT(in)  :: point_lat     ! ORIGINAL LATITUDE
    REAL(wp), INTENT(in)  :: p_alpha   ! ROTATION ANGLE
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(inout) :: p_rotlon  ! ROTATED LONGITUDE
    REAL(wp), INTENT(inout) :: p_rotlat  ! ROTATED LATITUDE
    
    ! !LOCAL VARIABLES:
    REAL(wp)              :: z_test    ! checking value
    !-----------------------------------------------------------------------

    IF (p_alpha == 0.0_wp) THEN       !        NO ROTATION
      
      p_rotlon = point_lon
      p_rotlat = point_lat
      
    ELSE                              !        ROTATION BY ANGLE p_alpha
      
      !     ROTATED LATITUDE
      
      z_test = SIN(point_lat)*COS(p_alpha)- COS(point_lat)*COS(point_lon)*SIN(p_alpha)
      
      IF (z_test > 1.0_wp) THEN
        p_rotlat = pi_2
      ELSEIF (z_test < -1.0_wp) THEN
        p_rotlat = -1.0_wp * pi_2
      ELSE
        p_rotlat = ASIN(z_test)
      ENDIF
      
      !     ROTATED LONGITUDE
      
      z_test = COS(p_rotlat)
      
      IF (z_test == 0.0_wp) THEN
        p_rotlon = 0.0_wp
      ELSE
        z_test = SIN(point_lon)*COS(point_lat)/z_test
        IF (z_test > 1.0_wp) THEN
          p_rotlon = pi_2
        ELSEIF (z_test < -1.0_wp) THEN
          p_rotlon = -1.0_wp * pi_2
        ELSE
          p_rotlon = ASIN(z_test)
        ENDIF
      ENDIF
      
      !        ADJUST FOR CORRECT BRANCH OF INVERSE SINE
      
      z_test = COS(p_alpha)*COS(point_lon)*COS(point_lat) + SIN(p_alpha)*SIN(point_lat)
      
      IF (z_test < 0.0_wp) THEN
        p_rotlon = pi - p_rotlon
      ENDIF
      
    ENDIF
    
  END SUBROUTINE rotate
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !
  ! !IROUTINE:  geostr_balance
  !
  ! !FUNCTION INTERFACE:
  FUNCTION geostr_balance( point_lat, func)  result(p_hh)
    !
    ! !DESCRIPTION:
    ! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
    ! to compute geostrophically balanced initial state used
    ! in test 3.
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Modified by Th.Heinze, DWD, (2006-11-22):
    ! - introduced INTERFACE uu (got an error message with g95 compiler,
    !   scanned the code, this seems to be the correct way, but might be wrong)
    ! Modified by Th.Heinze, DWD, (2006-12-12):
    ! - renamed it to geostr_balance
    !
    ! !REMARKS:
    ! was htmp2 in previous code
    
    ! !INTERFACE:
    INTERFACE                        ! selected function
      
      FUNCTION func(p_t) result(p_vv)
        
        USE mo_kind, ONLY: wp
        
        REAL(wp), INTENT(in) :: p_t
        REAL(wp)             :: p_vv
        
      END FUNCTION func
      
    END INTERFACE
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lat           ! rotated latitude
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh            ! balanced height
    
    ! !LOCAL VARIABLES:
    INTEGER :: j               ! loop index
    
    REAL(wp)             :: z_a             ! left bound
    REAL(wp)             :: z_b             ! right bound
    REAL(wp)             :: cell_lat           ! latitude in loop
    REAL(wp)             :: z_step          ! step
    REAL(wp)             :: z_val, z_val2   ! intermediate values
    !-----------------------------------------------------------------------
    
    z_a = -1._wp * pi_2
    z_b = point_lat
    
    z_step = 0.02_wp * ( z_b - z_a)
    
    p_hh = 0._wp
    
    cell_lat = z_a - 0.5_wp * z_step
    
    DO j = 1, 50
      cell_lat = cell_lat + z_step
      
      z_val = func(cell_lat)
      
      z_val2 = 2._wp * grid_angular_velocity * SIN(cell_lat)
      z_val2 = z_val2 + z_val * TAN(cell_lat)* sphere_radius
      z_val2 = z_val * z_val2
      
      p_hh = p_hh + z_val2 * z_step
      
    ENDDO
    
  END FUNCTION geostr_balance
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  FUNCTION test11_h(lon,lat) result(hh)
    
    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lon,lat
    REAL(wp) :: hh, hdach, alpha, beta, phi2
    REAL(wp)             :: z_rotlon  ! rotated longitude
    REAL(wp)             :: z_rotlat  ! rotated latitude
    
    ! rotate
    CALL rotate( lon, lat, aleph, z_rotlon, z_rotlat)
    ! calculate height
    hh = geostr_balance11( z_rotlat, test11_u2)
    
    
    hdach = 120._wp
    alpha = 1._wp/3._wp
    beta  = 1._wp/15._wp
    phi2  = pi/4._wp
    hh    = hh + hdach*COS(lat)*EXP(-((lon)/alpha)**2)*EXP(-((phi2-lat)/beta)**2)
    
  END FUNCTION test11_h
  !-----------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------
  FUNCTION test11_u(lat) result(uu)
    
    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lat
    REAL(wp) ::  uu, d
    REAL(wp) ::  phi0, phi1, umax, en
    
    phi0=pi/7._wp
    phi1=pi/2._wp - phi0
    en=EXP(-4._wp/(phi0-phi1)**2)
    umax=80._wp
    
    d=.1_wp
    
    IF ((lat.GT.phi0).AND.(lat.LT.phi1))THEN
      uu=umax/en*EXP(1._wp/(lat-phi0)/(lat-phi1))
      IF (uu.LT. 0.001_wp) THEN
        uu  = 0.0_wp
      END IF
      !         print*, "assigning u values", uu
    ELSE
      uu=0._wp
    ENDIF
    
    !    1451 !     For jet on southern hemisphere additionally:
    !    1452 ! if ((lat.lt.-phi0).and.(lat.gt.-phi1)) then
    !    1453 ! uu=+umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For volume tests
    !    1454 ! ! uu=uu-umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For Galewsky tests
    !    1455 ! endif
    
  END FUNCTION test11_u
  !-----------------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------------
  FUNCTION test11_u2(lat) result(uu)
    
    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lat
    REAL(wp) ::  uu, d
    REAL(wp) ::  phi0, phi1, umax, en
    
    phi0=pi/7._wp
    phi1=pi/2._wp - phi0
    en=EXP(-4._wp/(phi0-phi1)**2)
    umax=80._wp
    
    d=.1_wp
    
    IF ((lat.GT.phi0).AND.(lat.LT.phi1))THEN
      uu=umax/en*EXP(1._wp/(lat-phi0)/(lat-phi1))
      IF (uu.LT. 0.001_wp) THEN
        uu  = 0.0_wp
      END IF
      !         print*, "assigning u values", uu
    ELSE
      uu=0._wp
    ENDIF
    
    !     For jet on southern hemisphere additionally:
    ! if ((lat.lt.-phi0).and.(lat.gt.-phi1)) then
    ! uu=+umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For volume tests
    ! ! uu=uu-umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For Galewsky tests
    ! endif
    
  END FUNCTION test11_u2
  
  FUNCTION test11_v(lon,lat) result(vv)
    
    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lon,lat
    REAL(wp) ::  vv
    
    vv = lat
    vv = lon
    vv = 0.0_wp
    
  END FUNCTION test11_v
  !-----------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------
  FUNCTION geostr_balance11( phi, func)  result(p_hh)
    
    ! !DESCRIPTION:
    ! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
    ! to compute geostrophically balanced initial state used
    ! in test 3.
    !
    ! !REVISION HISTORY:
    ! Developed  by L.Bonaventura  (2002-5).
    ! Modified by Th.Heinze, DWD, (2006-11-22):
    ! - introduced INTERFACE uu (got an error message with g95 compiler,
    !   scanned the code, this seems to be the correct way, but might be wrong)
    ! Modified by Th.Heinze, DWD, (2006-12-12):
    ! - renamed it to geostr_balance
    ! Modified by F. Rauser, MPI (2009,10) for testcase 11 galewsky
    !
    ! !REMARKS:
    ! was htmp2 in previous code
    
    ! !INTERFACE:
    INTERFACE                        ! selected function
      
      FUNCTION func(p_t) result(p_vv)
        
        USE mo_kind, ONLY: wp
        
        REAL(wp), INTENT(in) :: p_t
        REAL(wp)             :: p_vv
        
      END FUNCTION func
      
    END INTERFACE
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: phi           ! rotated latitude
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh            ! balanced height
    ! !LOCAL VARIABLES:
    INTEGER :: j               ! loop index
    REAL(wp)             :: phi_a             ! left bound
    REAL(wp)             :: phi_b             ! right bound
    REAL(wp)             :: phidash           ! latitude in loop
    REAL(wp)             :: dphi          ! step
    REAL(wp)             :: u, temp   ! intermediate values

    !-----------------------------------------------------------------------

    phi_a = -0.5_wp * pi
    phi_b = phi
    
    dphi = 0.01_wp * ( phi_b - phi_a)
    
    p_hh = 0._wp
    
    phidash = phi_a - 0.5_wp * dphi
    
    DO j = 1, 100
      phidash = phidash + dphi
      
      u = func(phidash)
      
      temp = 2._wp * grid_angular_velocity * SIN(phidash)
      temp = temp + ( u * TAN(phidash)* sphere_radius)
      temp = sphere_radius *rgrav * u * temp
      
      p_hh = p_hh + temp * dphi
      
    ENDDO
    
    p_hh = 10000._wp - p_hh
    !     print*, "phh", INT(360*phi/pi), INT(p_hh)
    
  END FUNCTION geostr_balance11
  !-----------------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------------
  SUBROUTINE fill_tracer_x_height(patch_3d, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: tracer_idx, jk, jb, jc, start_cell_idx, end_cell_idx
    
    IF (.NOT. use_tracer_x_height) RETURN
    
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    
    DO tracer_idx = 1, no_tracer
      ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration_x_height(:, :, :) = 0.0_wp
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_idx, end_cell_idx)
        DO jc = start_cell_idx, end_cell_idx
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            
            ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration_x_height(jc, jk, jb) = &
              & ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration(jc,jk,jb)   *      &
              & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(jc, jk, jb)
            
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE fill_tracer_x_height
  !-----------------------------------------------------------------------------------
  
END MODULE mo_ocean_initial_conditions
