!>
!!  Contains the data structures for the hydrostatic ocean model.
!!
!!  Contains the data structures to store the hydrostatic & boussinesq ocean model state.
!!  Implementation is based on ICON-Shallow-Water model
!!  to store the shallow water model state and other auxiliary variables.
!!  Constructors and destructors for these data structures are also defined here.
!!
!! @par Revision History
!!  Initial version by Peter Korn (MPI-M), (2006).
!!  Big recoding by P. Korn (MPI-M), (2009/2010)
!!  Modification by Stephan Lorenz, MPI-M, (2010-03-19):
!!   - renaming and adjustment to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2011-07
!!   - 3-dim ocean structures moved from patch_oce to hydro_ocean_base
!
!
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
!!
!----------------------------
MODULE mo_oce_state
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: success, max_char_length
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, use_tracer_x_height, cfl_write
  USE mo_oce_types,           ONLY: t_hydro_ocean_base ,t_hydro_ocean_state ,t_hydro_ocean_prog ,t_hydro_ocean_diag, &
    &                               t_hydro_ocean_aux ,t_hydro_ocean_acc ,t_ptr3d ,t_oce_config ,t_ocean_tracer ,    &
    &                               t_ocean_regions ,t_ocean_region_volumes ,t_ocean_region_areas ,t_ocean_basins 
  USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier,my_process_is_mpi_test
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    &                               success, max_char_length, MIN_DOLIC,               &
    &                               full_coriolis, beta_plane_coriolis,                &
    &                               f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch,t_patch_3d, t_grid_cells, t_grid_edges
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, grid_sphere_radius, grid_angular_velocity, &
    & use_dummy_cell_closure
  USE mo_dynamics_config,     ONLY: nnew,nold
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, t_geographical_coordinates
  
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: add_var,                  &
    &                               new_var_list,             &
    &                               delete_var_list,          &
    &                               default_var_list_settings,&
    &                               add_ref
  USE mo_var_metadata,        ONLY: groups 
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  !  USE mo_ocean_config,        ONLY: ignore_land_points
  
  IMPLICIT NONE
  PRIVATE
  
  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  !public interface
  !
  ! subroutines
  PUBLIC :: construct_hydro_ocean_base
  PUBLIC :: destruct_hydro_ocean_base
  PUBLIC :: construct_hydro_ocean_state
  PUBLIC :: destruct_hydro_ocean_state
  PUBLIC :: construct_patch_3d, destruct_patch_3d
  PUBLIC :: set_oce_tracer_info
  !
  
  PUBLIC :: ocean_restart_list
  PUBLIC :: ocean_default_list
  PUBLIC :: v_base
  PUBLIC :: oce_config
  
  !constructors
  PRIVATE :: construct_hydro_ocean_diag
  PRIVATE :: construct_hydro_ocean_prog
  PRIVATE :: construct_hydro_ocean_aux
  !destructors
  PRIVATE :: destruct_hydro_ocean_diag
  PRIVATE :: destruct_hydro_ocean_aux
  !----------------------------------------------------------------------------
  
  ! variables
  TYPE(t_var_list)                              :: ocean_restart_list
  TYPE(t_var_list)                              :: ocean_default_list
  TYPE(t_hydro_ocean_base) ,TARGET :: v_base
  TYPE(t_oce_config)                            :: oce_config
  
  !-------------------------------------------------------------------------
  
CONTAINS
  
  
  !>
  !! Constructor for hydrostatic ocean state + diagnostic and auxiliary  states.
  !!
  !! Constructor for hydrostatic ocean state
  !! It calls  constructors to single time level
  !! auxiliary and diagnostic states. Then it constructs state array,
  !! whose components (representing multiple time levels).
  !! Initialization of all components with zero.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2007).
  !!  Modification by Stephan Lorenz, MPI-M, (2010-06-01) - no temporary memory array
  !
  !
  SUBROUTINE construct_hydro_ocean_state( patch_2d, p_os )
    
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d(n_dom)
    TYPE(t_hydro_ocean_state), TARGET :: p_os(n_dom)
    
    ! local variables
    INTEGER :: jg
    
    INTEGER :: i_status, jp, prlength ! local prognostic array length
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:construct_hydro_ocean_state'
    
    CALL message(TRIM(routine), 'start to construct hydro_ocean state' )
    
    ! Using Adams-Bashforth semi-implicit timestepping with 3 prognostic time levels:
    prlength = 3
    
    !create state array for each domain
    DO jg = 1, n_dom
      
      ALLOCATE(p_os(jg)%p_prog(1:prlength), stat=i_status)
      IF (i_status/=success) THEN
        CALL finish(TRIM(routine), 'allocation of progn. state array failed')
      END IF
      DO jp = 1, prlength
        CALL construct_hydro_ocean_prog(patch_2d(jg), p_os(jg)%p_prog(jp),jp)
      END DO
      
      CALL construct_hydro_ocean_diag(patch_2d(jg), p_os(jg)%p_diag)
      CALL construct_hydro_ocean_aux(patch_2d(jg),  p_os(jg)%p_aux)
      CALL construct_hydro_ocean_acc(patch_2d(jg),  p_os(jg)%p_acc)
      
      CALL message(TRIM(routine),'construction of hydrostatic ocean state finished')
      
    END DO
    
  END SUBROUTINE construct_hydro_ocean_state
  
  !-------------------------------------------------------------------------
  !>
  !!               Destructor for hydrostatic ocean state.
  !
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
  SUBROUTINE destruct_hydro_ocean_state(p_os)
    TYPE(t_hydro_ocean_state), TARGET,INTENT(inout)   :: p_os(n_dom)
    
    ! local variables
    
    INTEGER :: jg, prlength, ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:destruct_hydro_ocean_state'
    
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start to destruct hydro ocean state ')
    
    prlength = SIZE(p_os(1)%p_prog)
    
    IF (prlength==0) THEN
      CALL finish(TRIM(routine),'prog array has length zero')
    END IF
    
    CALL delete_var_list(ocean_restart_list)
    CALL delete_var_list(ocean_default_list)
    
    DO jg = 1, n_dom
      CALL destruct_hydro_ocean_diag(p_os(jg)%p_diag)
      CALL destruct_hydro_ocean_aux (p_os(jg)%p_aux)
      
      ! destruct state array
      ist = 1
      DEALLOCATE(p_os(jg)%p_prog, stat=ist)
      IF (ist/=success) THEN
        CALL finish(TRIM(routine),'deallocation of state array failed')
      END IF
    END DO
    
    CALL message(TRIM(routine),'destruction of hydrostatic ocean state finished')
    
  END SUBROUTINE destruct_hydro_ocean_state
  
  !-------------------------------------------------------------------------
  !>
  !! Allocation of basic 3-dimensional structure components of hydrostatic ocean state.
  !! Initialization of components with zero.
  !
  !
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2011/06).
  !!
  
  SUBROUTINE construct_hydro_ocean_base(patch_2d, v_base)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hydro_ocean_base), INTENT(inout)    :: v_base
    
    ! local variables
    
    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp, n_zlvm!, ie
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:construct_hydro_ocean_base'
    
    !-------------------------------------------------------------------------
    
    !CALL message(TRIM(routine), 'start to construct basic hydro ocean state')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    n_zlvp = n_zlev + 1
    n_zlvm = n_zlev - 1
    
    ! allocate and set vertical level thickness from the namelist
    ALLOCATE(v_base%del_zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_m failed')
    ENDIF
    ALLOCATE(v_base%zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_m failed')
    ENDIF
    ALLOCATE(v_base%zlev_i(n_zlvp),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_i failed')
    ENDIF
    ALLOCATE(v_base%del_zlev_i(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_i failed')
    ENDIF
    
    !
    !! 3-dim land-sea-mask at cells, edges and vertices
    !
    ! cells
    ALLOCATE(v_base%lsm_c(nproma,n_zlev,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating lsm_c failed')
    ENDIF
    ! edges
    ALLOCATE(v_base%lsm_e(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating lsm_e failed')
    ENDIF
    ! deepest ocean layer in column
    ALLOCATE(v_base%dolic_c(nproma,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating dolic_c failed')
    ENDIF
    ALLOCATE(v_base%dolic_e(nproma,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating dolic_e failed')
    ENDIF
    ! 2-dim basins and areas
    ALLOCATE(v_base%basin_c(nproma,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating basin_c failed')
    ENDIF
    ALLOCATE(v_base%regio_c(nproma,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating regio_c failed')
    ENDIF
    ! 3-dim real land-sea-mask
    ! cells
    ALLOCATE(v_base%wet_c(nproma,n_zlev,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating wet_c failed')
    ENDIF
    ! edges
    ALLOCATE(v_base%wet_e(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating wet_e failed')
    ENDIF
    
    v_base%del_zlev_m = 0._wp
    v_base%del_zlev_i = 0._wp
    v_base%zlev_m     = 0._wp
    v_base%zlev_i     = 0._wp
    
    v_base%lsm_c = 0
    v_base%lsm_e = 0
    v_base%dolic_c = 0
    v_base%dolic_e = 0
    v_base%basin_c = 0
    v_base%regio_c = 0
    
    v_base%wet_c = 0.0_wp
    v_base%wet_e = 0.0_wp
    
  END SUBROUTINE construct_hydro_ocean_base
  
  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2011/06).
  !!
  SUBROUTINE destruct_hydro_ocean_base(v_base)
    
    TYPE(t_hydro_ocean_base), INTENT(inout) :: v_base
    
    ! local variables
    
    INTEGER :: ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:destruct_hydro_ocean_base'
    
    CALL message(TRIM(routine),' start to destruct hydrostatic ocean basic state')
    
    DEALLOCATE(v_base%zlev_m,stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'deallocating zlev_m failed')
    ENDIF
    
  END SUBROUTINE destruct_hydro_ocean_base
  
  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean prognostic state.
  !!               Initialization of components with zero.
  !
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
  SUBROUTINE construct_hydro_ocean_prog(patch_2d, p_os_prog, timelevel)
    
    TYPE(t_patch), INTENT(in), TARGET :: patch_2d
    TYPE(t_hydro_ocean_prog), INTENT(inout)   :: p_os_prog
    INTEGER, INTENT(in)                       :: timelevel
    
    INTEGER :: alloc_cell_blocks, nblks_e !, nblks_v
    INTEGER :: jtrc
    INTEGER, PARAMETER :: max_oce_tracer = 2
    CHARACTER(LEN=max_char_length) :: oce_tracer_names(max_oce_tracer),&
      & oce_tracer_units(max_oce_tracer),&
      & oce_tracer_longnames(max_oce_tracer)
    INTEGER :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(LEN=max_char_length) :: var_suffix
    
    !-------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------
    REAL(wp), POINTER ::        &
      & h(:,:)                ,& ! height of the free surface. Unit: [m]
    ! dimension:(nproma, alloc_cell_blocks)
      & vn(:,:,:)             ,& ! velocity component normal to cell edge. Unit [m/s]
    ! dimension: (nproma, n_zlev, nblks_e)
      & t(:,:,:),s(:,:,:)          ! tracer concentration.
    !-------------------------------------------------------------------------
    WRITE(var_suffix,'(a,i2.2)') '_TL',timelevel
    
    !-------------------------------------------------------------------------
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    
    ! height
    CALL add_var(ocean_restart_list, 'h'//TRIM(var_suffix), p_os_prog%h , &
      & grid_unstructured_cell, za_surface,    &
      & t_cf_var('h'//TRIM(var_suffix), 'm', 'surface elevation at cell center', DATATYPE_FLT64),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/))
    
    !! normal velocity component
    CALL add_var(ocean_restart_list,'vn'//TRIM(var_suffix),p_os_prog%vn,grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('vn'//TRIM(var_suffix), 'm/s', 'normal velocity on edge', DATATYPE_FLT64),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/))
    
    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL set_oce_tracer_info(max_oce_tracer      , &
        & oce_tracer_names    , &
        & oce_tracer_longnames, &
        & oce_tracer_codes    , &
        & oce_tracer_units,     &
        & var_suffix)
      CALL add_var(ocean_restart_list, 'tracers'//TRIM(var_suffix), p_os_prog%tracer , &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity', &
        & DATATYPE_FLT64),&
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
      ! Reference to individual tracer, for I/O
      ALLOCATE(p_os_prog%tracer_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_restart_list, 'tracers'//TRIM(var_suffix),              &
          & oce_tracer_names(jtrc),                 &
          & p_os_prog%tracer_ptr(jtrc)%p,                             &
          & grid_unstructured_cell, za_depth_below_sea,               &
          & t_cf_var(TRIM(oce_tracer_names(jtrc))//TRIM(var_suffix), &
          & oce_tracer_units(jtrc), &
          & oce_tracer_longnames(jtrc), DATATYPE_FLT64), &
          & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/))
      END DO
      
      ! use of the ocean_tracers structure
      ALLOCATE(p_os_prog%ocean_tracers(no_tracer))
      DO jtrc = 1,no_tracer
        ! point the concentration to the 4D tracer
        ! this is a tmeporary solution until the whole code is cleaned
        p_os_prog%ocean_tracers(jtrc)%concentration =>  p_os_prog%tracer(:,:,:,jtrc)
        
        ! allocate a
        IF (use_tracer_x_height) THEN
          !
          CALL add_var(ocean_restart_list, 'ocean_tracers'//TRIM(var_suffix), &
            & p_os_prog%ocean_tracers(jtrc)%concentration_x_height , &
            & grid_unstructured_cell, za_depth_below_sea, &
            & t_cf_var(TRIM(oce_tracer_names(jtrc))//"_x_height",  &
            & oce_tracer_units(jtrc),                            &
            & oce_tracer_longnames(jtrc), DATATYPE_FLT64),       &
            & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
            & ldims=(/nproma,n_zlev,alloc_cell_blocks/))
          
        ENDIF ! use_tracer_x_height
      ENDDO
      
    ENDIF ! no_tracer > 0
    
  END SUBROUTINE construct_hydro_ocean_prog
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean diagnostic state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
  
  SUBROUTINE construct_hydro_ocean_diag(patch_2d,p_os_diag)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hydro_ocean_diag), INTENT(inout)    :: p_os_diag
    
    ! local variables
    
    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    !INTEGER :: jb, jc, jk, je
    !INTEGER :: i_startidx_c, i_endidx_c, i_startidx_e, i_endidx_e
    INTEGER, PARAMETER :: max_oce_tracer = 2
    CHARACTER(LEN=max_char_length) :: oce_tracer_names(max_oce_tracer),&
      & oce_tracer_units(max_oce_tracer),&
      & oce_tracer_longnames(max_oce_tracer)
    INTEGER :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:construct_hydro_ocean_diag'
    
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start to construct diagnostic hydro ocean state')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    
    
    !    CALL add_var(ocean_default_list, 'test', x , GRID_1x1,&
    !    &            no_z_axis, &
    !    &            t_cf_var('monitor_T', 'C', 'monitoring_temperature', DATATYPE_FLT32),&
    !    &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    !    &            ldims=(/1,1,1/),in_group=groups("oce_monitoring"))
    
    
    CALL add_var(ocean_default_list, 'rho', p_os_diag%rho , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rho', 'kg/m^3', 'insitu density', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'rhopot', p_os_diag%rhopot , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rhopot', 'kg/m^3', 'potential density', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'zgrad_rho', p_os_diag%zgrad_rho , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('zgrad_rho', 'kg/m^3', 'vertical density gradiant', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
!   CALL add_var(ocean_default_list, 'vt', p_os_diag%vt, grid_unstructured_edge, &
!     & za_depth_below_sea, &
!     & t_cf_var('vt','m/s','tangential velocity at edges', DATATYPE_FLT32),&
!     & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference,grid_edge),&
!     & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'h_e', p_os_diag%h_e, grid_unstructured_edge,&
      & za_surface, &
      & t_cf_var('h_e','m','surface height ar edges', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference,grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_diag"))
    ! thicknesses
    CALL add_var(ocean_default_list, 'thick_c', p_os_diag%thick_c,  &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('thick_c','m','fluid column thickness at cells', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference,grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag"))
    CALL add_var(ocean_default_list, 'thick_e', p_os_diag%thick_e, &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('thick_e','m','fluid column thickness at edges', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference,grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_restart_list, 'div_mass_flx_c', p_os_diag%div_mass_flx_c,&
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div mass flux','','divergence mass flux at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list, 'mass_flux', p_os_diag%mass_flx_e, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('mass flux','',' mass flux at edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"))
    
    ! velocities
    CALL add_var(ocean_restart_list, 'w', p_os_diag%w, grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('w','m/s','vertical velocity at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'w_old', p_os_diag%w_old, grid_unstructured_cell, &
      & za_depth_below_sea_half,&
      & t_cf_var('w_old','m/s','vertical velocity at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
!   CALL add_var(ocean_restart_list, 'w_e', p_os_diag%w_e, grid_unstructured_cell, &
!     & za_depth_below_sea_half, &
!     & t_cf_var('w_e','m/s','vertical velocity at edges', DATATYPE_FLT32),&
!     & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
!     & ldims=(/nproma,n_zlev+1,nblks_e/),lrestart_cont=.TRUE.)
!   CALL add_var(ocean_default_list, 'w_prev', p_os_diag%w_prev, &
!     & grid_unstructured_edge, za_depth_below_sea_half, &
!     & t_cf_var('w_prev','m/s','vertical velocity at edges', DATATYPE_FLT32),&
!     & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
!     & ldims=(/nproma,n_zlev+1,nblks_e/))
    ! reconstructed u velocity component
    CALL add_var(ocean_default_list, 'u', p_os_diag%u, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('u','m/s','u zonal velocity component', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    ! reconstructed v velocity component
    CALL add_var(ocean_default_list, 'v', p_os_diag%v, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('v','m/s','v meridional velocity component', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    ! reconstrcuted velocity in cartesian coordinates
    !   CALL add_var(ocean_restart_list, 'p_vn', p_os_diag%p_vn, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA, &
    !   &            t_cf_var('p_vn','m/s','normal velocity in cartesian coordinates', DATATYPE_FLT32),&
    !   &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    !   &            ldims=(/nproma,n_zlev,alloc_cell_blocks/))
    ! integrated barotropic stream function
    CALL add_var(ocean_restart_list, 'u_vint', p_os_diag%u_vint, grid_unstructured_cell, &
      & za_surface, &
      & t_cf_var('u_vint','m*m/s','barotropic zonal velocity', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'ptp_vn', p_os_diag%ptp_vn, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('ptp_vn','m/s','normal velocity in cartesian coordinates', &
      & DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),lrestart_cont=.TRUE.)
    ! predicted vn normal velocity component
    CALL add_var(ocean_restart_list, 'vn_pred', p_os_diag%vn_pred, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred','m/s','predicted vn normal velocity component', &
      & DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    ! predicted vn normal velocity component
!    CALL add_var(ocean_restart_list, 'vn_impl_vert_diff', p_os_diag%vn_impl_vert_diff,&
!      & grid_unstructured_edge, za_depth_below_sea, &
!      & t_cf_var('vn_impl_vert_diff','m/s','predicted vn normal velocity component', &
!      & DATATYPE_FLT32),&
!      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
!      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list, 'vn_time_weighted', p_os_diag%vn_time_weighted, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred','m/s','average vn normal velocity component', &
      & DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_default_list, 'w_time_weighted', p_os_diag%w_time_weighted, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('w_prev','m/s','vertical velocity at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    ! vorticity
    CALL add_var(ocean_restart_list, 'vort', p_os_diag%vort, &
      & grid_unstructured_vert, za_depth_below_sea, &
      & t_cf_var('vort','1/s','vorticity', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_vertex),&
      & ldims=(/nproma,n_zlev,nblks_v/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
!   CALL add_var(ocean_restart_list, 'vort_e', p_os_diag%vort_e, &
!     & grid_unstructured_edge, za_depth_below_sea, &
!     & t_cf_var('vort_e','1/s','vorticity at edges', DATATYPE_FLT32),&
!     & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
!     & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! kinetic energy component
    CALL add_var(ocean_default_list, 'kin', p_os_diag%kin, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('kin','J','kinetic energy', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    ! gradient term
    CALL add_var(ocean_restart_list, 'grad', p_os_diag%grad, grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('grad','','gradient', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! divergence component
    CALL add_var(ocean_restart_list, 'div', p_os_diag%div, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div','','divergence', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.TRUE.)
    
    ! pressures
    CALL add_var(ocean_restart_list, 'press_hyd', p_os_diag%press_hyd, grid_unstructured_cell,&
      & za_depth_below_sea, t_cf_var('press_hyd','','hydrostatic pressure', &
      & DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'press_grad', p_os_diag%press_grad, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('press_grad','',' pressure gradient', &
      & DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! horizontal velocity advection
    CALL add_var(ocean_restart_list, 'veloc_adv_horz', p_os_diag%veloc_adv_horz, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('veloc_adv_horz','','horizontal velocity advection', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! vertical velocity advection
    CALL add_var(ocean_restart_list, 'veloc_adv_vert', p_os_diag%veloc_adv_vert, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('veloc_adv_vert','','vertical velocity advection', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! horizontal diffusion
    CALL add_var(ocean_restart_list, 'laplacian_horz', p_os_diag%laplacian_horz, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('laplacian_horz','','horizontal diffusion', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    ! vertical diffusion
    CALL add_var(ocean_restart_list, 'laplacian_vert', p_os_diag%laplacian_vert, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('laplacian_vert','','vertical diffusion', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),lrestart_cont=.TRUE.)
    ! mixed layer depths
    CALL add_var(ocean_default_list, 'mld', p_os_diag%mld , grid_unstructured_cell,za_surface, &
      &          t_cf_var('mld', 'm', 'mixed layer depth', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag","oce_default","oce_essentials"))
    CALL add_var(ocean_default_list, 'condep', p_os_diag%condep , grid_unstructured_cell, za_surface,&
      &         t_cf_var('condep', '', 'convection depth index', DATATYPE_INT8),&
      &         t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      &         ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag","oce_default","oce_essentials"))
    IF (cfl_write) THEN
    CALL add_var(ocean_default_list, 'cfl_vert', p_os_diag%cfl_vert , &
      &          grid_unstructured_cell, za_depth_below_sea_half,&
      &          t_cf_var('cdf_vert', '', 'vertical cfl relation', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      &          ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
      &          in_group=groups("oce_diag"))
    CALL add_var(ocean_default_list, 'cfl_horz', p_os_diag%cfl_horz, &
      &          grid_unstructured_edge, za_depth_below_sea,&
      &          t_cf_var('cfl_horz', '', 'horizontal cfl relation', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      &          ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"))
    ENDIF
    
    !reconstrcuted velocity in cartesian coordinates
    ALLOCATE(p_os_diag%p_vn(nproma,n_zlev,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at cells failed')
    END IF
    
    ALLOCATE(p_os_diag%p_vn_dual(nproma,n_zlev,nblks_v), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at verts failed')
    END IF
    
    !reconstrcuted velocity at edges in cartesian coordinates
    ALLOCATE(p_os_diag%p_vn_mean(nproma,n_zlev,nblks_e), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn_mean at edges failed')
    END IF
    ALLOCATE(p_os_diag%p_mass_flux_sfc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_mass_flux_sfc_cc at cells failed')
    END IF
    ! set all values - incl. last block - of cartesian coordinates to zero (NAG compiler)
    p_os_diag%p_vn     (:,:,:)%x(1)=0.0_wp
    p_os_diag%p_vn     (:,:,:)%x(2)=0.0_wp
    p_os_diag%p_vn     (:,:,:)%x(3)=0.0_wp
    p_os_diag%p_vn_dual(:,:,:)%x(1)=0.0_wp
    p_os_diag%p_vn_dual(:,:,:)%x(2)=0.0_wp
    p_os_diag%p_vn_dual(:,:,:)%x(3)=0.0_wp
    p_os_diag%p_vn_mean(:,:,:)%x(1)=0.0_wp
    p_os_diag%p_vn_mean(:,:,:)%x(2)=0.0_wp
    p_os_diag%p_vn_mean(:,:,:)%x(3)=0.0_wp
    
    p_os_diag%p_mass_flux_sfc_cc(:,:)%x(1)=0.0_wp
    p_os_diag%p_mass_flux_sfc_cc(:,:)%x(2)=0.0_wp
    p_os_diag%p_mass_flux_sfc_cc(:,:)%x(3)=0.0_wp
    
    !remapped velocity at cell edges
    ALLOCATE(p_os_diag%ptp_vn(nproma,n_zlev,nblks_e), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for ptp_vn at edges failed')
    END IF
    ! initialize all components with zero (this is preliminary)
    p_os_diag%ptp_vn    = 0.0_wp
    
    CALL add_var(ocean_restart_list,'temp_insitu',p_os_diag%temp_insitu,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'K', 'in situ temperature', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    !--------------------------------------------------------------------------
    !Add output pf prognostic variables with readable names
    CALL add_var(ocean_default_list,'h',p_os_diag%h , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('h', 'm', 'surface elevation at cell center', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_prog"),&
      & loutput=.TRUE., lrestart=.FALSE.)
    CALL add_var(ocean_default_list,'vn',p_os_diag%vn, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn', 'm/s', 'normal velocity on edge', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_prog"), &
      & loutput=.TRUE.,lrestart=.FALSE.)
    CALL add_var(ocean_default_list, 't',p_os_diag%t,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('t','degC','potential temperature', DATATYPE_FLT32), &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_prog"),&
      & loutput=.TRUE., lrestart=.FALSE.)
    CALL add_var(ocean_default_list, 's',p_os_diag%s,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('s','degC','salinity', DATATYPE_FLT32), &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_prog"),&
      & loutput=.TRUE., lrestart=.FALSE.)
  END SUBROUTINE construct_hydro_ocean_diag
  
  
  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
  SUBROUTINE destruct_hydro_ocean_diag(p_os_diag)
    
    TYPE(t_hydro_ocean_diag), INTENT(inout) :: p_os_diag
    
    ! local variables
    
    INTEGER :: ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:destruct_hydro_ocean_diag'
    
    DEALLOCATE(p_os_diag%p_vn, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn failed')
    END IF
    DEALLOCATE(p_os_diag%p_vn_dual, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn_dual failed')
    END IF
    DEALLOCATE(p_os_diag%ptp_vn, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for ptp_vn failed')
    END IF
    
  END SUBROUTINE destruct_hydro_ocean_diag
  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean auxiliary state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
  SUBROUTINE construct_hydro_ocean_aux(patch_2d, p_os_aux)
    
    TYPE(t_patch),TARGET, INTENT(in)                :: patch_2d
    TYPE(t_hydro_ocean_aux), TARGET,INTENT(inout)   :: p_os_aux
    
    ! local variables
    
    INTEGER ::  ist  !, jtrc
    INTEGER ::  alloc_cell_blocks, nblks_e, nblks_v
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:construct_hydro_ocean_aux'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start to construct hydro ocean auxiliary state')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    
    ! allocation for Adam-Bashford time stepping
    CALL add_var(ocean_restart_list,'g_n',p_os_aux%g_n, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_n','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_aux"),loutput=.TRUE.,lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'g_nm1',p_os_aux%g_nm1, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_nm1','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_aux"),loutput=.TRUE.,lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'g_nimd',p_os_aux%g_nimd, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_nimd','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_aux"),loutput=.TRUE.,lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list,'p_rhs_sfc_eq',p_os_aux%p_rhs_sfc_eq, grid_unstructured_cell,&
      & za_surface, t_cf_var('p_rhs_sfc_eq','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    
    ! allocation for boundary conditions
    CALL add_var(ocean_restart_list,'bc_top_u',p_os_aux%bc_top_u, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_u','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_v',p_os_aux%bc_top_v, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_v','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_vn',p_os_aux%bc_top_vn, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_top_vn','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list,'bc_bot_u',p_os_aux%bc_bot_u, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_u','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_bot_v',p_os_aux%bc_bot_v, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_v','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_bot_vn',p_os_aux%bc_bot_vn, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_bot_vn','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list,'bc_bot_w',p_os_aux%bc_bot_w, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_w','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_w',p_os_aux%bc_top_w, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_w','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list,'bc_bot_tracer',p_os_aux%bc_bot_tracer,grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_tracer','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks,no_tracer/),in_group=groups("oce_aux"))
    CALL add_var(ocean_default_list,'bc_top_tracer',p_os_aux%bc_top_tracer,grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_tracer','','', DATATYPE_FLT32),&
      & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks,no_tracer/),in_group=groups("oce_aux"))
    
    ALLOCATE(p_os_aux%bc_top_veloc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    ALLOCATE(p_os_aux%bc_bot_veloc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    
    p_os_aux%bc_top_veloc_cc(:,:)%x(1) = 0.0_wp
    p_os_aux%bc_top_veloc_cc(:,:)%x(2) = 0.0_wp
    p_os_aux%bc_top_veloc_cc(:,:)%x(3) = 0.0_wp
    p_os_aux%bc_bot_veloc_cc(:,:)%x(1) = 0.0_wp
    p_os_aux%bc_bot_veloc_cc(:,:)%x(2) = 0.0_wp
    p_os_aux%bc_bot_veloc_cc(:,:)%x(3) = 0.0_wp
    
    ! allocation of 3-dim tracer relaxation:
    IF (no_tracer >= 1) THEN
      CALL add_var(ocean_default_list,'data_3dimRelax_Temp',p_os_aux%data_3dimRelax_Temp,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('data_3dimRelax_Temp','','', DATATYPE_FLT32),&
        & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
      CALL add_var(ocean_default_list,'forc_3dimRelax_Temp',p_os_aux%forc_3dimRelax_Temp,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('forc_3dimRelax_Temp','','', DATATYPE_FLT32),&
        & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
    END IF
    IF (no_tracer == 2) THEN
      CALL add_var(ocean_default_list,'data_3dimRelax_Salt',p_os_aux%data_3dimRelax_Salt,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('data_3dimRelax_Salt','','', DATATYPE_FLT32),&
        & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
      CALL add_var(ocean_default_list,'forc_3dimRelax_Salt',p_os_aux%forc_3dimRelax_Salt,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('forc_3dimRelax_Salt','','', DATATYPE_FLT32),&
        & t_grib2_var(255,255,255,DATATYPE_PACK16,grid_reference, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
    END IF
    
  END SUBROUTINE construct_hydro_ocean_aux
  
  SUBROUTINE construct_hydro_ocean_acc(patch_2d, p_os_acc)
    
    TYPE(t_patch),TARGET, INTENT(in)                :: patch_2d
    TYPE(t_hydro_ocean_acc), TARGET,INTENT(inout)   :: p_os_acc
    
    ! local variables
    
    INTEGER ::  ist, jtrc
    INTEGER ::  alloc_cell_blocks, nblks_e, nblks_v
    
    INTEGER, PARAMETER :: max_oce_tracer = 2
    CHARACTER(LEN=max_char_length) :: oce_tracer_names(max_oce_tracer),&
      & oce_tracer_units(max_oce_tracer),&
      & oce_tracer_longnames(max_oce_tracer)
    INTEGER :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(LEN=max_char_length) :: var_suffix
    
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    
    CALL add_var(ocean_default_list, 'h_acc', p_os_acc%h , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('h_acc', 'm', 'surface elevation at cell center', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    CALL add_var(ocean_default_list, 'u_acc', p_os_acc%u, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('u_acc','m/s','zonal velocity component', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    ! reconstructed v velocity component
    CALL add_var(ocean_default_list, 'v_acc', p_os_acc%v, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('v_acc','m/s','meridional velocity component', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    CALL add_var(ocean_default_list, 'rhopot_acc', p_os_acc%rhopot, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('rhopot_acc','kg/m^3','potential density', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default","oce_essentials"))
    CALL add_var(ocean_default_list, 'rho_acc', p_os_acc%rho, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('rho_acc','kg/m^3','insitu density', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    
    CALL add_var(ocean_default_list, 'w_acc', p_os_acc%w, grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('w_acc','m/s','vertical velocity', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    CALL add_var(ocean_default_list, 'mass_flx_e_acc', p_os_acc%mass_flx_e, grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('mass_flx_e_acc','','mass flux at edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'div_mass_flx_c_acc', p_os_acc%div_mass_flx_c, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div_mass_flx_c_acc','','divergence of mass flux', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'u_vint_acc', p_os_acc%u_vint , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('u_vint_acc', 'm*m/s', 'barotropic zonal velocity', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    CALL add_var(ocean_default_list, 'vort_acc', p_os_acc%vort, &
      & grid_unstructured_vert, za_depth_below_sea, &
      & t_cf_var('vort_acc','1/s','vorticity', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_vertex),&
      & ldims=(/nproma,n_zlev,nblks_v/),in_group=groups("oce_default"),lrestart_cont=.TRUE.)
    
    ! kinetic energy component
    CALL add_var(ocean_default_list, 'kin_acc', p_os_acc%kin, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('kin_acc','J','kinetic energy', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    ! & ptp_vn(:,:,:)           ,& ! normal velocity after mapping P^T P
    ! & vn_pred(:,:,:)          ,& ! predicted normal velocity vector at edges.
    ! & vn_impl_vert_diff(:,:,:),& ! predicted normal velocity vector at edges.
    ! & vn_time_weighted(:,:,:) ,&  ! predicted normal velocity vector at edges.
    ! & w_time_weighted(:,:,:)  ,& ! predicted normal velocity vector at edges.
    ! & veloc_adv_horz(:,:,:)   ,& ! horizontal velocity advection
    ! & veloc_adv_vert(:,:,:)   ,& ! vertical velocity advection
    ! & laplacian_horz(:,:,:)   ,& ! horizontal diffusion of horizontal velocity
    ! & laplacian_vert(:,:,:)   ,& ! vertical diffusion of horizontal velocity
    ! & grad(:,:,:)             ,& ! gradient of kinetic energy. Unit [m/s]
    ! & div(:,:,:)              ,& ! divergence. Unit [m/s]
    ! & press_hyd(:,:,:)        ,& ! hydrostatic pressure. Unit [m]
    ! & press_grad(:,:,:)       ,& ! hydrostatic pressure gradient term. Unit [m/s]
    ! & temp_insitu(:,:,:)      ,&
    
    !-------------------------------------------------------------------------
    WRITE(var_suffix,'(a)') '_acc'
    !-------------------------------------------------------------------------
    
    IF ( no_tracer > 0 ) THEN
      CALL set_oce_tracer_info(max_oce_tracer      , &
        & oce_tracer_names    , &
        & oce_tracer_longnames, &
        & oce_tracer_codes    , &
        & oce_tracer_units,var_suffix)
      CALL add_var(ocean_default_list, 'tracers'//TRIM(var_suffix), p_os_acc%tracer , &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity', &
        & DATATYPE_FLT32),&
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
      ! Reference to individual tracer, for I/O
      ALLOCATE(p_os_acc%tracer_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_default_list, 'tracers'//TRIM(var_suffix),&
          & oce_tracer_names(jtrc),                 &
          & p_os_acc%tracer_ptr(jtrc)%p,           &
          & grid_unstructured_cell, za_depth_below_sea,&
          & t_cf_var(oce_tracer_names(jtrc), &
          & oce_tracer_units(jtrc), &
          & oce_tracer_longnames(jtrc), DATATYPE_FLT32), &
          & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
        
      END DO
    ENDIF ! no_tracer > 0
    
  END SUBROUTINE construct_hydro_ocean_acc
  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of auxilliary hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
  SUBROUTINE destruct_hydro_ocean_aux(p_os_aux)
    
    TYPE(t_hydro_ocean_aux), INTENT(inout)      :: p_os_aux
    
    ! local variables
    
    INTEGER :: ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:destruct_hydro_ocean_aux'
    
    DEALLOCATE(p_os_aux%bc_top_veloc_cc, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond cc failed')
    END IF
    
    DEALLOCATE(p_os_aux%bc_bot_veloc_cc, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine),'deallocation of bot boundary cond cc failed')
    END IF
    
  END SUBROUTINE destruct_hydro_ocean_aux
  !-------------------------------------------------------------------------
  
  
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE destruct_patch_3d(patch_3d)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)    :: patch_3d

    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_m)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%inv_del_zlev_m)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_m)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_i)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_i)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_area)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_volume)
    DEALLOCATE(patch_3d%p_patch_1d)
    
  END SUBROUTINE destruct_patch_3d
  
  !-------------------------------------------------------------------------
  !>
  !! Allocation of basic 3-dimensional patch structure. This sbr assumes that
  !! the 2D horizontal patch components is already initialized.
  !
  !
  !! @par Revision History
  !! Developed  by  Peter korn, MPI-M (2012/08).
  !!
  SUBROUTINE construct_patch_3d(patch_3d)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)    :: patch_3d
    
    ! local variables
    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp, n_zlvm!, ie
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_oce_state:construct_patch_3D'
    
    !-------------------------------------------------------------------------
    
    !CALL message(TRIM(routine), 'start to construct basic hydro ocean state')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_3d%p_patch_2d(n_dom)%alloc_cell_blocks
    nblks_e = patch_3d%p_patch_2d(n_dom)%nblks_e
    nblks_v = patch_3d%p_patch_2d(n_dom)%nblks_v
    n_zlvp = n_zlev + 1
    n_zlvm = n_zlev - 1
    
    ALLOCATE(patch_3d%p_patch_1d(n_dom_start:n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating p_patch_1D failed')
    ENDIF
    
    ! allocate and set vertical level thickness from the namelist
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_m failed')
    ENDIF
    ! allocate the inverse of the above
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%inv_del_zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating inv_del_zlev_m failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_m failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_i(n_zlvp),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_i failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_i(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_i failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_area(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating ocean_area failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_volume(n_zlvp),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating ocean_volume failed')
    ENDIF
    
    !
    !! 3-dim land-sea-mask at cells, edges and vertices
    !
    ! cells
    CALL add_var(ocean_default_list, 'lsm_c', patch_3d%lsm_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('lsm_c','','3d lsm on cells', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default"),&
      & isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'lsm_e', patch_3d%lsm_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('lsm_e','','3d lsm on edges', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry","oce_default"),&
      & isteptype=tstep_constant)
    ! surface cells
    CALL add_var(ocean_default_list, 'surface_cell_sea_land_mask', patch_3d%surface_cell_sea_land_mask , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('surface_cell_sea_land_mask', '', 'surface_cell_sea_land_mask', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_vertex),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! surface vertices
    CALL add_var(ocean_default_list, 'surface_edge_sea_land_mask', patch_3d%surface_edge_sea_land_mask , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('surface_edge_sea_land_mask', '', 'surface_edge_sea_land_mask', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_vertex),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! surface vertices
    CALL add_var(ocean_default_list, 'surface_vertex_sea_land_mask', patch_3d%surface_vertex_sea_land_mask , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('surface_vertex_sea_land_mask', '', 'surface_vertex_sea_land_mask', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_vertex),&
      & ldims=(/nproma,nblks_v/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! deepest ocean layer in column
    CALL add_var(ocean_default_list, 'dolic_c', patch_3d%p_patch_1d(n_dom)%dolic_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('dolic_c', '', 'dolic_c', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'dolic_e', patch_3d%p_patch_1d(n_dom)%dolic_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('dolic_e', '', 'dolic_e', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! 2-dim basins and areas
    CALL add_var(ocean_default_list, 'basin_c', patch_3d%basin_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('basin_c', '', 'basin_c', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'regio_c', patch_3d%regio_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('regio_c', '', 'regio_c', DATATYPE_INT8),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    ! 2-dim bottom and column thickness
    CALL add_var(ocean_default_list, 'bottom_thick_c', patch_3d%bottom_thick_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('bottom_thick_c', 'm', 'bottom_thick_c', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'bottom_thick_e', patch_3d%bottom_thick_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('bottom_thick_e', 'm', 'bottom_thick_e', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'column_thick_c', patch_3d%column_thick_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('column_thick_c', 'm', 'column_thick_c', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'column_thick_e', patch_3d%column_thick_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('column_thick_e', 'm', 'column_thick_e', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! 3-dim real land-sea-mask
    ! cells
    CALL add_var(ocean_default_list, 'wet_c', patch_3d%wet_c , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('wet_c', '', '3d lsm on cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'wet_e', patch_3d%wet_e , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('wet_e', '', '3d lsm on edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    ! 3-dim real land-sea-mask with zero on halos
    ! cells
    CALL add_var(ocean_default_list, 'wet_halo_zero_c', patch_3d%wet_halo_zero_c , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('wet_c_halo_zero', '', '3d lsm with halo zero on cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'wet_halo_zero_e', patch_3d%wet_halo_zero_e , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('wet_e_halo_zero', '', '3d lsm with halo zero on edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    
    patch_3d%p_patch_1d(n_dom)%del_zlev_m = 0._wp
    patch_3d%p_patch_1d(n_dom)%del_zlev_i = 0._wp
    patch_3d%p_patch_1d(n_dom)%zlev_m     = 0._wp
    patch_3d%p_patch_1d(n_dom)%zlev_i     = 0._wp
    
    patch_3d%p_patch_1d(n_dom)%ocean_area(:)   = 0._wp
    patch_3d%p_patch_1d(n_dom)%ocean_volume(:) = 0._wp
    
    CALL add_var(ocean_default_list, 'prism_thick_c', patch_3d%p_patch_1d(1)%prism_thick_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism thickness at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_volume', patch_3d%p_patch_1d(1)%prism_volume, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism volume (cells)', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_thick_e', patch_3d%p_patch_1d(n_dom)%prism_thick_e, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism thickness at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_thick_flat_sfc_c', patch_3d%p_patch_1d(n_dom)%prism_thick_flat_sfc_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_c','m','time independent depth at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_thick_flat_sfc_e', patch_3d%p_patch_1d(n_dom)%prism_thick_flat_sfc_e, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_c','m','time independent depth at edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inverse prism_thick_c', patch_3d%p_patch_1d(n_dom)%inv_prism_thick_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inverse prism_thick_c','m','time independent depth at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_center_dist_c', patch_3d%p_patch_1d(n_dom)%prism_center_dist_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_center_dist_c','m','dist between prism centers', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inv_prism_thick_e', patch_3d%p_patch_1d(n_dom)%inv_prism_thick_e, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_c','m','time independent depth at edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inv_prism_center_dist_c', &
      & patch_3d%p_patch_1d(n_dom)%inv_prism_center_dist_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inv_prism_center_dist_c','1/m','inverse of dist between prism centers at cells', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inv_prism_center_dist_e', &
      & patch_3d%p_patch_1d(n_dom)%inv_prism_center_dist_e, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inv_prism_center_dist_e','1/m','inverse of dist betweenprism centers at edges', DATATYPE_FLT32),&
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    
  END SUBROUTINE construct_patch_3d
  
  !-------------------------------------------------------------------------
  !>
  !!
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2011).
  !!
  SUBROUTINE set_oce_tracer_info(max_oce_tracer,&
    & oce_tracer_names,&
    & oce_tracer_longnames,&
    & oce_tracer_codes,&
    & oce_tracer_units,&
    & suffix)
    
    INTEGER, INTENT(in)            :: max_oce_tracer
    CHARACTER(LEN=max_char_length) :: oce_tracer_names(max_oce_tracer),&
      & oce_tracer_units(max_oce_tracer),&
      & oce_tracer_longnames(max_oce_tracer)
    INTEGER :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(LEN=max_char_length), OPTIONAL :: suffix
    
    IF (max_oce_tracer < no_tracer) THEN
      CALL finish('set_oce_tracer_info','Too many tracers! Please provide trace info')
    ENDIF
    IF (PRESENT(suffix)) THEN
      !     write(0,*)'suffix:',suffix
    END IF
    oce_tracer_names(1)     = 't'
    IF (PRESENT(suffix)) THEN
      oce_tracer_names(1) = 't'//TRIM(suffix)
    END IF
    oce_tracer_longnames(1) = 'potential temperature'
    oce_tracer_units(1)     = 'deg C'
    oce_tracer_codes(1)     = 200
    
    oce_tracer_names(2)     = 's'
    IF (PRESENT(suffix)) THEN
      oce_tracer_names(2) = 's'//TRIM(suffix)
    END IF
    oce_tracer_longnames(2) = 'salinity'
    oce_tracer_units(2)     = 'psu'
    oce_tracer_codes(2)     = 201
    
  END SUBROUTINE
  !-------------------------------------------------------------------------
  
  
END MODULE mo_oce_state
