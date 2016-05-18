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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!----------------------------
MODULE mo_ocean_state
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: success, max_char_length
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, use_tracer_x_height, cfl_write,&
    &                               Cartesian_Mixing , &
    &                               k_tracer_dianeutral_parameter,                          &
    &                               k_tracer_isoneutral_parameter, k_tracer_GM_kappa_parameter,    &
    &                               GMRedi_configuration,GMRedi_combined,                   &
    &                               GM_only,Redi_only, type_3dimrelax_salt, type_3dimrelax_temp,  &
    &                               nbgctra, nbgcadv,lhamocc
  USE mo_ocean_types,           ONLY: t_hydro_ocean_base ,t_hydro_ocean_state ,t_hydro_ocean_prog ,t_hydro_ocean_diag, &
    &                               t_hydro_ocean_aux ,t_hydro_ocean_acc, t_oce_config ,t_ocean_tracer
  USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier,my_process_is_mpi_test
  USE mo_parallel_config,     ONLY: nproma
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
  USE mo_grib2,               ONLY: grib2_var, t_grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_INT8, DATATYPE_PACK16, &
    &                               tstep_constant, GRID_LONLAT, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell, grid_unstructured_edge, &
      &                             grid_unstructured_vert, grid_vertex, za_depth_below_sea, za_depth_below_sea_half, za_surface
  !  USE mo_ocean_config,        ONLY: ignore_land_points
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  
  IMPLICIT NONE
  PRIVATE
  
 
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
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_state( patch_3d, ocean_state )
    
    TYPE(t_patch_3D), TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state(n_dom)
    
    ! local variables
    TYPE(t_patch), POINTER :: patch_2d
   
    INTEGER :: i_status, jp, prlength ! local prognostic array length
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_hydro_ocean_state'

    patch_2d => patch_3d%p_patch_2d(1)
    CALL message(TRIM(routine), 'start to construct hydro_ocean state' )
    ocean_state(1)%patch_3d => patch_3d
    ! Using Adams-Bashforth semi-implicit timestepping with 3 prognostic time levels:
    prlength = 3
    
    !create state array for each domain     
    ALLOCATE(ocean_state(1)%p_prog(1:prlength), stat=i_status)
    IF (i_status/=success) THEN
      CALL finish(TRIM(routine), 'allocation of progn. state array failed')
    END IF
    DO jp = 1, prlength
      CALL construct_hydro_ocean_prog(patch_2d, ocean_state(1)%p_prog(jp),jp)
    END DO
    
    CALL construct_hydro_ocean_diag(patch_2d, ocean_state(1)%p_diag)
    CALL construct_hydro_ocean_aux(patch_2d,  ocean_state(1)%p_aux)
    CALL construct_hydro_ocean_acc(patch_2d,  ocean_state(1)%p_acc)
    
    CALL message(TRIM(routine),'construction of hydrostatic ocean state finished')      
    
  END SUBROUTINE construct_hydro_ocean_state
  
  !-------------------------------------------------------------------------
  !>
  !!               Destructor for hydrostatic ocean state.
  !
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_hydro_ocean_state(ocean_state)
    TYPE(t_hydro_ocean_state), TARGET,INTENT(inout)   :: ocean_state(n_dom)
    
    ! local variables
    
    INTEGER :: jg, prlength, ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_state'
    
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start to destruct hydro ocean state ')
    
    prlength = SIZE(ocean_state(1)%p_prog)
    
    IF (prlength==0) THEN
      CALL finish(TRIM(routine),'prog array has length zero')
    END IF
    
    CALL delete_var_list(ocean_restart_list)
    CALL delete_var_list(ocean_default_list)
    
    DO jg = 1, n_dom
      CALL destruct_hydro_ocean_diag(ocean_state(jg)%p_diag)
      CALL destruct_hydro_ocean_aux (ocean_state(jg)%p_aux)
      
      ! destruct state array
      ist = 1
      DEALLOCATE(ocean_state(jg)%p_prog, stat=ist)
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
  
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_base(patch_2d, v_base)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hydro_ocean_base), INTENT(inout)    :: v_base
    
    ! local variables
    
    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_hydro_ocean_base'
    
    !-------------------------------------------------------------------------
    
    !CALL message(TRIM(routine), 'start to construct basic hydro ocean state')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    n_zlvp = n_zlev + 1
    
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
      & routine = 'mo_ocean_state:destruct_hydro_ocean_base'
    
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
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_prog(patch_2d, ocean_state_prog, timelevel)
    
    TYPE(t_patch), INTENT(in), TARGET :: patch_2d
    TYPE(t_hydro_ocean_prog), INTENT(inout)   :: ocean_state_prog
    INTEGER, INTENT(in)                       :: timelevel
    
    INTEGER :: alloc_cell_blocks, nblks_e !, nblks_v
    INTEGER :: jtrc
    INTEGER, PARAMETER :: max_oce_tracer = 2
    CHARACTER(LEN=max_char_length) :: oce_tracer_names(max_oce_tracer+nbgctra),&
      & oce_tracer_units(max_oce_tracer+nbgctra),&
      & oce_tracer_longnames(max_oce_tracer+nbgctra)
    INTEGER :: oce_tracer_codes(max_oce_tracer+nbgctra)
    CHARACTER(LEN=max_char_length) :: var_suffix
    
    !-------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------
    REAL(wp), POINTER ::        &
      & h(:,:)                ,& ! height of the free surface. Unit: [m]
    ! dimension:(nproma, alloc_cell_blocks)
      & vn(:,:,:)             ,& ! velocity component normal to cell edge. Unit [m/s]
    ! dimension: (nproma, n_zlev, nblks_e)
      & t(:,:,:),s(:,:,:)          ! tracer concentration.
    INTEGER                       :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    WRITE(var_suffix,'(a,i2.2)') '_TL',timelevel
    
    !-------------------------------------------------------------------------
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    
    ! height
    CALL add_var(ocean_restart_list, 'h'//TRIM(var_suffix), ocean_state_prog%h , &
      & grid_unstructured_cell, za_surface,    &
      & t_cf_var('h'//TRIM(var_suffix), 'm', 'surface elevation at cell center', DATATYPE_FLT64),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/))
    
    !! normal velocity component
    CALL add_var(ocean_restart_list,'vn'//TRIM(var_suffix),ocean_state_prog%vn,grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('vn'//TRIM(var_suffix), 'm/s', 'normal velocity on edge', DATATYPE_FLT64),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/))
    
    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL set_oce_tracer_info(max_oce_tracer      , &
        & oce_tracer_names    , &
        & oce_tracer_longnames, &
        & oce_tracer_codes    , &
        & oce_tracer_units,     &
        & var_suffix)
      if(lhamocc)then 
      CALL set_bgc_tracer_info(no_tracer,max_oce_tracer+nbgctra      , &
        & oce_tracer_names    , &
        & oce_tracer_longnames, &
        & oce_tracer_codes    , &
        & oce_tracer_units,     &
        & var_suffix)
      endif
      CALL add_var(ocean_restart_list, 'tracers'//TRIM(var_suffix), ocean_state_prog%tracer , &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity', &
        & DATATYPE_FLT64),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer+nbgctra/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
      ! Reference to individual tracer, for I/O
      ALLOCATE(ocean_state_prog%tracer_ptr(no_tracer+nbgctra))
      DO jtrc = 1,no_tracer+nbgctra
        CALL add_ref( ocean_restart_list, 'tracers'//TRIM(var_suffix),              &
          & oce_tracer_names(jtrc),                 &
          & ocean_state_prog%tracer_ptr(jtrc)%p,                             &
          & grid_unstructured_cell, za_depth_below_sea,               &
          & t_cf_var(TRIM(oce_tracer_names(jtrc))//TRIM(var_suffix), &
          & oce_tracer_units(jtrc), &
          & oce_tracer_longnames(jtrc), DATATYPE_FLT64), &
          & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/))
      END DO
      
      ! use of the ocean_tracers structure
      ALLOCATE(ocean_state_prog%ocean_tracers(no_tracer+nbgctra))
      DO jtrc = 1,no_tracer+nbgctra
        ! point the concentration to the 4D tracer
        ! this is a tmeporary solution until the whole code is cleaned
        ocean_state_prog%ocean_tracers(jtrc)%concentration =>  ocean_state_prog%tracer(:,:,:,jtrc)
        
        ! allocate a
!         IF (use_tracer_x_height) THEN
!           !
!           CALL add_var(ocean_restart_list, 'ocean_tracers'//TRIM(var_suffix), &
!             & ocean_state_prog%ocean_tracers(jtrc)%concentration_x_height , &
!             & grid_unstructured_cell, za_depth_below_sea, &
!             & t_cf_var(TRIM(oce_tracer_names(jtrc))//"_x_height",  &
!             & oce_tracer_units(jtrc),                            &
!             & oce_tracer_longnames(jtrc), DATATYPE_FLT64),       &
!             & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
!             & ldims=(/nproma,n_zlev,alloc_cell_blocks/))
!           
!         ENDIF ! use_tracer_x_height
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
  
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_diag(patch_2d,ocean_state_diag)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hydro_ocean_diag), INTENT(inout)    :: ocean_state_diag
    
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
      & routine = 'mo_ocean_state:construct_hydro_ocean_diag'
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF
    
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct diagnostic hydro ocean state...')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    
    ! add monitoring 
    CALL add_var(ocean_default_list, 'volume_Global', ocean_state_diag%monitor%volume , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('volume', 'm^3', 'volume', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'kin_energy_Global', ocean_state_diag%monitor%kin_energy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('kin_energy', 'J', 'kin_energy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'pot_energy_Global', ocean_state_diag%monitor%pot_energy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('pot_energy', 'J', 'pot_energy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'total_energy_Global', ocean_state_diag%monitor%total_energy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('total_energy', 'J', 'total_energy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'total_salt_Global', ocean_state_diag%monitor%total_salt , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('total_salt', 'kg', 'total_salt', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'vorticity_Global', ocean_state_diag%monitor%vorticity , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('vorticity', '', 'vorticity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'enstrophy_Global', ocean_state_diag%monitor%enstrophy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('enstrophy', '', 'enstrophy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'potential_enstrophy_Global', ocean_state_diag%monitor%potential_enstrophy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('potential_enstrophy', '', 'potential_enstrophy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'absolute_vertical_velocity_Global', ocean_state_diag%monitor%absolute_vertical_velocity , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('absolute_vertical_velocity', 'm/s', 'absolute_vertical_velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_ShortWave_Global', ocean_state_diag%monitor%HeatFlux_ShortWave , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_ShortWave', 'W/m2', 'HeatFlux_ShortWave', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_LongWave_Global', ocean_state_diag%monitor%HeatFlux_LongWave , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_LongWave', 'W/m2', 'HeatFlux_LongWave_Global', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Sensible_Global', ocean_state_diag%monitor%HeatFlux_Sensible , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Sensible', 'W/m2', 'HeatFlux_Sensible', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Latent_Global', ocean_state_diag%monitor%HeatFlux_Latent , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Latent', 'W/m2', 'HeatFlux_Latent', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Total_Global', ocean_state_diag%monitor%HeatFlux_Total , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Total', 'W/m2', 'HeatFlux_Total', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Precipitation_Global', ocean_state_diag%monitor%FrshFlux_Precipitation , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Precipitation', 'm/s', 'FrshFlux_Precipitation', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_SnowFall_Global', ocean_state_diag%monitor%FrshFlux_SnowFall , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_SnowFall', 'm/s', 'FrshFlux_SnowFall_Global', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Evaporation_Global', ocean_state_diag%monitor%FrshFlux_Evaporation , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Evaporation', 'm/s', 'FrshFlux_Evaporation_Global', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Runoff_Global', ocean_state_diag%monitor%FrshFlux_Runoff , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Runoff', 'm/s', 'FrshFlux_Runoff', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_TotalSalt_Global', ocean_state_diag%monitor%FrshFlux_TotalSalt , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_TotalSalt', 'm/s', 'FrshFlux_TotalSalt', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_TotalOcean_Global', ocean_state_diag%monitor%FrshFlux_TotalOcean , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_TotalOcean', 'm/s', 'FrshFlux_TotalOcean', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_TotalIce_Global', ocean_state_diag%monitor%FrshFlux_TotalIce , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_TotalIce', 'm/s', 'FrshFlux_TotalIce', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_VolumeIce_Global', ocean_state_diag%monitor%FrshFlux_VolumeIce , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_VolumeIce', 'm/s', 'FrshFlux_VolumeIce', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_VolumeTotal_Global', ocean_state_diag%monitor%FrshFlux_VolumeTotal , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_VolumeTotal', 'm/s', 'FrshFlux_VolumeTotal', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Relax_Global', ocean_state_diag%monitor%HeatFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Relax', 'W/m2', 'HeatFlux_Relax', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Relax_Global', ocean_state_diag%monitor%FrshFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Relax', 'm/m', 'FrshFlux_Relax', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'TempFlux_Relax_Global', ocean_state_diag%monitor%TempFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('TempFlux_Relax', 'T/s', 'TempFlux_Relax', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'SaltFlux_Relax_Global', ocean_state_diag%monitor%SaltFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('SaltFlux_Relax', 'psu/s', 'SaltFlux_Relax', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_volume_nh', ocean_state_diag%monitor%ice_volume_nh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_volume_nh', 'km^3', 'ice_volume_nh', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_volume_sh', ocean_state_diag%monitor%ice_volume_sh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_volume_sh', 'km^3', 'ice_volume_sh', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_extent_nh', ocean_state_diag%monitor%ice_extent_nh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_extent_nh', 'km^2', 'ice_extent_nh', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_extent_sh', ocean_state_diag%monitor%ice_extent_sh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_extent_sh', 'km^2', 'ice_extent_sh', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_framStrait', ocean_state_diag%monitor%ice_framStrait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_framStrait', 'm^3/s', 'ice_framStrait', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'gibraltar', ocean_state_diag%monitor%gibraltar , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('gibraltar', 'Sv', 'gibraltar', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'denmark_strait', ocean_state_diag%monitor%denmark_strait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('denmark_strait', 'Sv', 'denmark_strait', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'drake_passage', ocean_state_diag%monitor%drake_passage , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('drake_passage', 'Sv', 'drake_passage', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'indonesian_throughflow', ocean_state_diag%monitor%indonesian_throughflow , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('indonesian_throughflow', 'Sv', 'indonesian_throughflow', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'scotland_iceland', ocean_state_diag%monitor%scotland_iceland , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('scotland_iceland', 'Sv', 'scotland_iceland', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'mozambique', ocean_state_diag%monitor%mozambique , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('mozambique', 'Sv', 'mozambique', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'framStrait', ocean_state_diag%monitor%framStrait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('framStrait', 'Sv', 'framStrait', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'beringStrait', ocean_state_diag%monitor%beringStrait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('beringStrait', 'Sv', 'beringStrait', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'barentsOpening', ocean_state_diag%monitor%barentsOpening , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('barentsOpening', 'Sv', 'barentsOpening', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'agulhas', ocean_state_diag%monitor%agulhas , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('agulhas', 'Sv', 'agulhas', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'agulhas_long', ocean_state_diag%monitor%agulhas_long , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('agulhas_long', 'Sv', 'agulhas_long', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'agulhas_longer', ocean_state_diag%monitor%agulhas_longer , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('agulhas_longer', 'Sv', 'agulhas_longer', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'florida_strait', ocean_state_diag%monitor%florida_strait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('florida_strait', 'Sv', 'florida_strait', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 't_mean_na_200m', ocean_state_diag%monitor%t_mean_na_200m , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('t_mean_na_200m', 'degC', 't_mean_na_200m', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 't_mean_na_800m', ocean_state_diag%monitor%t_mean_na_800m , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('t_mean_na_800m', 'degC', 't_mean_na_800m', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_ocean_heat_budget', ocean_state_diag%monitor%ice_ocean_heat_budget , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_ocean_heat_budget', '', 'ice_ocean_heat_budget', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_ocean_salinity_budget', ocean_state_diag%monitor%ice_ocean_salinity_budget , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_ocean_salinity_budget', '', 'ice_ocean_salinity_budget', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_ocean_volume_budget', ocean_state_diag%monitor%ice_ocean_volume_budget , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_ocean_volume_budget', '', 'ice_ocean_volume_budget', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat),&
      & in_group=groups("ocean_monitor"),ldims=(/1/))

    
    CALL add_var(ocean_default_list, 'rho', ocean_state_diag%rho , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rho', 'kg/m^3', 'insitu density', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'rhopot', ocean_state_diag%rhopot , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rhopot', 'kg/m^3', 'potential density', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'zgrad_rho', ocean_state_diag%zgrad_rho , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('zgrad_rho', 'kg/m^3', 'vertical density gradiant', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
!   CALL add_var(ocean_default_list, 'vt', ocean_state_diag%vt, grid_unstructured_edge, &
!     & za_depth_below_sea, &
!     & t_cf_var('vt','m/s','tangential velocity at edges', datatype_flt),&
!     & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED,grid_edge),&
!     & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'h_e', ocean_state_diag%h_e, grid_unstructured_edge,&
      & za_surface, &
      & t_cf_var('h_e','m','surface height ar edges', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED,grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_diag"))
    ! thicknesses
    CALL add_var(ocean_default_list, 'thick_c', ocean_state_diag%thick_c,  &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('thick_c','m','fluid column thickness at cells', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED,grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag"))
    CALL add_var(ocean_default_list, 'thick_e', ocean_state_diag%thick_e, &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('thick_e','m','fluid column thickness at edges', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED,grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_diag"))
    
    CALL add_var(ocean_default_list, 'div_mass_flx_c', ocean_state_diag%div_mass_flx_c,&
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div mass flux','','divergence mass flux at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)
    
    CALL add_var(ocean_default_list, 'mass_flux', ocean_state_diag%mass_flx_e, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('mass flux','',' mass flux at edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)
    
    ! velocities
    CALL add_var(ocean_restart_list, 'w', ocean_state_diag%w, grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('w','m/s','vertical velocity at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'w_old', ocean_state_diag%w_old, grid_unstructured_cell, &
      & za_depth_below_sea_half,&
      & t_cf_var('w_old','m/s','vertical velocity at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
!   CALL add_var(ocean_restart_list, 'w_e', ocean_state_diag%w_e, grid_unstructured_cell, &
!     & za_depth_below_sea_half, &
!     & t_cf_var('w_e','m/s','vertical velocity at edges', datatype_flt),&
!     & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
!     & ldims=(/nproma,n_zlev+1,nblks_e/),lrestart_cont=.TRUE.)
!   CALL add_var(ocean_default_list, 'w_prev', ocean_state_diag%w_prev, &
!     & grid_unstructured_edge, za_depth_below_sea_half, &
!     & t_cf_var('w_prev','m/s','vertical velocity at edges', datatype_flt),&
!     & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
!     & ldims=(/nproma,n_zlev+1,nblks_e/))
    ! reconstructed u velocity component
    CALL add_var(ocean_default_list, 'u', ocean_state_diag%u, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('u','m/s','u zonal velocity component', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    ! reconstructed v velocity component
    CALL add_var(ocean_default_list, 'v', ocean_state_diag%v, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('v','m/s','v meridional velocity component', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    ! reconstrcuted velocity in cartesian coordinates
    !   CALL add_var(ocean_restart_list, 'p_vn', ocean_state_diag%p_vn, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA, &
    !   &            t_cf_var('p_vn','m/s','normal velocity in cartesian coordinates', datatype_flt),&
    !   &            grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
    !   &            ldims=(/nproma,n_zlev,alloc_cell_blocks/))
    ! integrated barotropic stream function
    CALL add_var(ocean_default_list, 'u_vint', ocean_state_diag%u_vint, grid_unstructured_cell, &
      & za_surface, &
      & t_cf_var('u_vint','m*m/s','barotropic zonal velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)
    CALL add_var(ocean_default_list, 'v_vint', ocean_state_diag%v_vint, grid_unstructured_cell, &
      & za_surface, &
      & t_cf_var('v_vint','m*m/s','barotropic meridional velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(ocean_restart_list, 'ptp_vn', ocean_state_diag%ptp_vn, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('ptp_vn','m/s','normal velocity in cartesian coordinates', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),loutput=.FALSE., lrestart_cont=.TRUE.)
    ! predicted vn normal velocity component

    CALL add_var(ocean_restart_list, 'zlim', ocean_state_diag%zlim, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('zlim','1','zalesak limiter factor', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),loutput=.true., lrestart_cont=.false.)

    CALL add_var(ocean_restart_list, 'vn_pred', ocean_state_diag%vn_pred, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred','m/s','predicted vn normal velocity component', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list, 'vn_pred_ptp', ocean_state_diag%vn_pred_ptp, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred_ptp','m/s','transformed predicted vn normal velocity component', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)      
      
    ! predicted vn normal velocity component
!    CALL add_var(ocean_restart_list, 'vn_impl_vert_diff', ocean_state_diag%vn_impl_vert_diff,&
!      & grid_unstructured_edge, za_depth_below_sea, &
!      & t_cf_var('vn_impl_vert_diff','m/s','predicted vn normal velocity component', &
!      & datatype_flt),&
!      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
!      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list, 'vn_time_weighted', ocean_state_diag%vn_time_weighted, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred','m/s','average vn normal velocity component', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_default_list, 'w_time_weighted', ocean_state_diag%w_time_weighted, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('w_prev','m/s','vertical velocity at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    ! vorticity
    CALL add_var(ocean_restart_list, 'vort', ocean_state_diag%vort, &
      & grid_unstructured_vert, za_depth_below_sea, &
      & t_cf_var('vort','1/s','vorticity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_vertex),&
      & ldims=(/nproma,n_zlev,nblks_v/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
  
!    CALL add_var(ocean_restart_list, 'potential_vort_e', ocean_state_diag%potential_vort_e, &
!      & grid_unstructured_edge, za_depth_below_sea, &
!      & t_cf_var('vort_e','1/s','potential vorticity at edges', datatype_flt),&
!      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_unstructured, grid_edge),&
!      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_essentials"),lrestart_cont=.TRUE.)

!    CALL add_var(ocean_restart_list, 'potential vort_c', ocean_state_diag%potential_vort_c, &
!      & grid_unstructured_cell, za_depth_below_sea, &
!      & t_cf_var('vort_e','1/s','potential vorticity at cells', datatype_flt),&
!      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
!      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_essentials"),lrestart_cont=.TRUE.)
    
    ! kinetic energy component
    CALL add_var(ocean_default_list, 'kin', ocean_state_diag%kin, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('kin','J','kinetic energy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"))
    
    ! gradient term
    CALL add_var(ocean_default_list, 'grad', ocean_state_diag%grad, grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('grad','','gradient', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)
    
    ! divergence component
    CALL add_var(ocean_default_list, 'div', ocean_state_diag%div, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div','','divergence', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.FALSE.)
    
    ! pressures
    CALL add_var(ocean_restart_list, 'press_hyd', ocean_state_diag%press_hyd, grid_unstructured_cell,&
      & za_depth_below_sea, t_cf_var('press_hyd','','hydrostatic pressure', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'press_grad', ocean_state_diag%press_grad, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('press_grad','',' pressure gradient', &
      & datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! horizontal velocity advection
    CALL add_var(ocean_restart_list, 'veloc_adv_horz', ocean_state_diag%veloc_adv_horz, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('veloc_adv_horz','','horizontal velocity advection', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
    ! vertical velocity advection
    CALL add_var(ocean_default_list, 'veloc_adv_vert', ocean_state_diag%veloc_adv_vert, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('veloc_adv_vert','','vertical velocity advection', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)
    
    ! horizontal diffusion
    CALL add_var(ocean_default_list, 'laplacian_horz', ocean_state_diag%laplacian_horz, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('laplacian_horz','','horizontal diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)
    ! vertical diffusion
    CALL add_var(ocean_default_list, 'laplacian_vert', ocean_state_diag%laplacian_vert, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('laplacian_vert','','vertical diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),lrestart_cont=.FALSE.)
    ! mixed layer depths
    CALL add_var(ocean_default_list, 'mld', ocean_state_diag%mld , grid_unstructured_cell,za_surface, &
      &          t_cf_var('mld', 'm', 'mixed layer depth', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag","oce_default","oce_essentials"))
    ! number of deepest convection layer
    CALL add_var(ocean_default_list, 'condep', ocean_state_diag%condep , grid_unstructured_cell, za_surface,&
      &         t_cf_var('condep', '', 'convection depth index', DATATYPE_INT8),&
      &         grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      &         ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_diag","oce_default","oce_essentials"))
    IF (cfl_write) THEN
    CALL add_var(ocean_default_list, 'cfl_vert', ocean_state_diag%cfl_vert , &
      &          grid_unstructured_cell, za_depth_below_sea_half,&
      &          t_cf_var('cdf_vert', '', 'vertical cfl relation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      &          ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
      &          in_group=groups("oce_diag"))
    CALL add_var(ocean_default_list, 'cfl_horz', ocean_state_diag%cfl_horz, &
      &          grid_unstructured_edge, za_depth_below_sea,&
      &          t_cf_var('cfl_horz', '', 'horizontal cfl relation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      &          ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"))
    ENDIF
    
    !reconstrcuted velocity in cartesian coordinates
    ALLOCATE(ocean_state_diag%p_vn(nproma,n_zlev,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at cells failed')
    END IF
    
    ALLOCATE(ocean_state_diag%p_vn_dual(nproma,n_zlev,nblks_v), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at verts failed')
    END IF
    
    !reconstrcuted velocity at edges in cartesian coordinates
    ALLOCATE(ocean_state_diag%p_vn_mean(nproma,n_zlev,nblks_e), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn_mean at edges failed')
    END IF
    ALLOCATE(ocean_state_diag%p_mass_flux_sfc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for p_mass_flux_sfc_cc at cells failed')
    END IF
    ! set all values - incl. last block - of cartesian coordinates to zero (NAG compiler)
    ocean_state_diag%p_vn     (:,:,:)%x(1)=0.0_wp
    ocean_state_diag%p_vn     (:,:,:)%x(2)=0.0_wp
    ocean_state_diag%p_vn     (:,:,:)%x(3)=0.0_wp
    ocean_state_diag%p_vn_dual(:,:,:)%x(1)=0.0_wp
    ocean_state_diag%p_vn_dual(:,:,:)%x(2)=0.0_wp
    ocean_state_diag%p_vn_dual(:,:,:)%x(3)=0.0_wp
    ocean_state_diag%p_vn_mean(:,:,:)%x(1)=0.0_wp
    ocean_state_diag%p_vn_mean(:,:,:)%x(2)=0.0_wp
    ocean_state_diag%p_vn_mean(:,:,:)%x(3)=0.0_wp
    
    ocean_state_diag%p_mass_flux_sfc_cc(:,:)%x(1)=0.0_wp
    ocean_state_diag%p_mass_flux_sfc_cc(:,:)%x(2)=0.0_wp
    ocean_state_diag%p_mass_flux_sfc_cc(:,:)%x(3)=0.0_wp
    
    !remapped velocity at cell edges
!     ALLOCATE(ocean_state_diag%ptp_vn(nproma,n_zlev,nblks_e), stat=ist)
!     IF (ist/=success) THEN
!       CALL finish(TRIM(routine), 'allocation for ptp_vn at edges failed')
!     END IF
    ! initialize all components with zero (this is preliminary)
    ocean_state_diag%ptp_vn    = 0.0_wp
    
    CALL add_var(ocean_default_list,'temp_insitu',ocean_state_diag%temp_insitu,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'K', 'in situ temperature', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)
    
!     CALL add_var(ocean_restart_list,'temp_horDiffused',ocean_state_diag%temp_horizontally_diffused, grid_unstructured_cell,&
!       & za_depth_below_sea, &
!       & t_cf_var('temp_insitu', 'K', 'horizonatlly diffused temperature', datatype_flt),&
!       & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
!       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.TRUE.)

    !--------------------------------------------------------------------------
    !Add output pf prognostic variables with readable names
    CALL add_var(ocean_default_list,'h',ocean_state_diag%h , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('h', 'm', 'surface elevation at cell center', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_prog"),&
      & loutput=.TRUE., lrestart=.FALSE.)
!     CALL add_var(ocean_default_list,'vn',ocean_state_diag%vn, &
!       & grid_unstructured_edge, za_depth_below_sea, &
!       & t_cf_var('vn', 'm/s', 'normal velocity on edge', datatype_flt),&
!       & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
!       & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_prog"), &
!       & loutput=.TRUE.,lrestart=.FALSE.)
    CALL add_var(ocean_default_list, 't',ocean_state_diag%t,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('t','degC','potential temperature', datatype_flt), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_prog"),&
      & loutput=.TRUE., lrestart=.FALSE.)
    CALL add_var(ocean_default_list, 's',ocean_state_diag%s,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('s','degC','salinity', datatype_flt), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_prog"),&
      & loutput=.TRUE., lrestart=.FALSE.)    
    
    
!       CALL add_var(ocean_restart_list, 'tracers'//TRIM(var_suffix), ocean_state_prog%tracer , &
!         & grid_unstructured_cell, za_depth_below_sea, &
!         & t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity', &
!         & DATATYPE_FLT64),&
!         & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
!         & ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer/), &
!         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
!       
!       ! Reference to individual tracer, for I/O
!       ALLOCATE(ocean_state_prog%tracer_ptr(no_tracer))
!       DO jtrc = 1,no_tracer
!         CALL add_ref( ocean_restart_list, 'tracers'//TRIM(var_suffix),              &
!           & oce_tracer_names(jtrc),                 &
!           & ocean_state_prog%tracer_ptr(jtrc)%p,                             &
!           & grid_unstructured_cell, za_depth_below_sea,               &
!           & t_cf_var(TRIM(oce_tracer_names(jtrc))//TRIM(var_suffix), &
!           & oce_tracer_units(jtrc), &
!           & oce_tracer_longnames(jtrc), DATATYPE_FLT64), &
!           & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
!           & ldims=(/nproma,n_zlev,alloc_cell_blocks/))
!       END DO
!       
!       ! use of the ocean_tracers structure
!       ALLOCATE(ocean_state_prog%ocean_tracers(no_tracer))
!       DO jtrc = 1,no_tracer
!         ocean_state_prog%ocean_tracers(jtrc)%concentration =>  ocean_state_prog%tracer(:,:,:,jtrc)
!       ENDDO
!     

   IF(GMRedi_configuration/=Cartesian_Mixing)THEN
    
      CALL add_var(ocean_restart_list, 'GMRedi_flux_horz',ocean_state_diag%GMRedi_flux_horz, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('GMRedi_flux_horz', '', '1:temperature 2:salinity', &
        & DATATYPE_FLT64),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e,no_tracer+nbgcadv/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(ocean_restart_list, 'GMRedi_flux_vert',ocean_state_diag%GMRedi_flux_vert, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('GMRedi_flux_vert', '', '1:temperature 2:salinity', &
        & DATATYPE_FLT64),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks,no_tracer+nbgcadv/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    
        ocean_state_diag%GMRedi_flux_horz(:,:,:,:)=0.0_wp
        ocean_state_diag%GMRedi_flux_vert(:,:,:,:)=0.0_wp   
    ENDIF  
      
  END SUBROUTINE construct_hydro_ocean_diag
  
  
  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_hydro_ocean_diag(ocean_state_diag)
    
    TYPE(t_hydro_ocean_diag), INTENT(inout) :: ocean_state_diag
    
    ! local variables
    
    INTEGER :: ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_diag'
    
    DEALLOCATE(ocean_state_diag%p_vn, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn failed')
    END IF
    DEALLOCATE(ocean_state_diag%p_vn_dual, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn_dual failed')
    END IF
!     DEALLOCATE(ocean_state_diag%ptp_vn, stat=ist)
!     IF (ist/=success) THEN
!       CALL finish(TRIM(routine), 'deallocation for ptp_vn failed')
!     END IF
    
  END SUBROUTINE destruct_hydro_ocean_diag
  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean auxiliary state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_aux(patch_2d, ocean_state_aux)
    
    TYPE(t_patch),TARGET, INTENT(in)                :: patch_2d
    TYPE(t_hydro_ocean_aux), TARGET,INTENT(inout)   :: ocean_state_aux
    
    ! local variables
    
    INTEGER ::  ist  !, jtrc
    INTEGER ::  alloc_cell_blocks, nblks_e, nblks_v
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_hydro_ocean_aux'
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct hydro ocean auxiliary state...')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    
    ! allocation for Adam-Bashford time stepping
    CALL add_var(ocean_restart_list,'g_n',ocean_state_aux%g_n, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_n','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_aux"),loutput=.TRUE.,lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'g_nm1',ocean_state_aux%g_nm1, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_nm1','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_aux"),loutput=.TRUE.,lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'g_nimd',ocean_state_aux%g_nimd, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_nimd','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_aux"),loutput=.TRUE.,lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list,'p_rhs_sfc_eq',ocean_state_aux%p_rhs_sfc_eq, grid_unstructured_cell,&
      & za_surface, t_cf_var('p_rhs_sfc_eq','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    
    ! allocation for boundary conditions
    CALL add_var(ocean_restart_list,'bc_top_u',ocean_state_aux%bc_top_u, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_u','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_v',ocean_state_aux%bc_top_v, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_v','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_vn',ocean_state_aux%bc_top_vn, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_top_vn','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_WindStress',ocean_state_aux%bc_top_WindStress, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_top_WindStress','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_aux"),lrestart_cont=.FALSE.)
    
    CALL add_var(ocean_restart_list,'bc_bot_vn',ocean_state_aux%bc_bot_vn, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_bot_vn','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    
    CALL add_var(ocean_restart_list,'bc_bot_w',ocean_state_aux%bc_bot_w, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_w','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list,'bc_top_w',ocean_state_aux%bc_top_w, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_w','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_aux"),lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list,'bc_bot_tracer',ocean_state_aux%bc_bot_tracer,grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_tracer','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks,no_tracer+nbgcadv/),in_group=groups("oce_aux"))
    CALL add_var(ocean_default_list,'bc_top_tracer',ocean_state_aux%bc_top_tracer,grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_tracer','','', datatype_flt),&
      & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks,no_tracer+nbgcadv/),in_group=groups("oce_aux"))
    
    ALLOCATE(ocean_state_aux%bc_top_veloc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    
    ocean_state_aux%bc_top_veloc_cc(:,:)%x(1) = 0.0_wp
    ocean_state_aux%bc_top_veloc_cc(:,:)%x(2) = 0.0_wp
    ocean_state_aux%bc_top_veloc_cc(:,:)%x(3) = 0.0_wp
    
    ! allocation of 3-dim tracer relaxation:
    IF (no_tracer>=1 .AND. type_3dimrelax_temp >0) THEN
      CALL add_var(ocean_default_list,'data_3dimRelax_Temp',ocean_state_aux%data_3dimRelax_Temp,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('data_3dimRelax_Temp','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
      CALL add_var(ocean_default_list,'forc_3dimRelax_Temp',ocean_state_aux%forc_3dimRelax_Temp,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('forc_3dimRelax_Temp','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.TRUE.)
      CALL add_var(ocean_default_list,'relax_3dim_coefficient',ocean_state_aux%relax_3dim_coefficient,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('relax_3dim_coefficient','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.TRUE.)
!       ocean_state_aux%relax_3dim_coefficient(:,:,:) = 1.0_wp 
    END IF
    IF (no_tracer==2 .AND. type_3dimrelax_salt >0) THEN
      CALL add_var(ocean_default_list,'data_3dimRelax_Salt',ocean_state_aux%data_3dimRelax_Salt,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('data_3dimRelax_Salt','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
      CALL add_var(ocean_default_list,'forc_3dimRelax_Salt',ocean_state_aux%forc_3dimRelax_Salt,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('forc_3dimRelax_Salt','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
    END IF
    
   IF(GMRedi_configuration/=Cartesian_mixing)THEN
   
     ALLOCATE(ocean_state_aux%slopes(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for slopes at cells failed')
     END IF
     ocean_state_aux%slopes    (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%slopes    (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%slopes    (:,:,:)%x(3)=0.0_wp
     ALLOCATE(ocean_state_aux%PgradTemperature_horz_center(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for PgradTemperature_horz_center at cells failed')
     END IF
     ocean_state_aux%PgradTemperature_horz_center    (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%PgradTemperature_horz_center    (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%PgradTemperature_horz_center    (:,:,:)%x(3)=0.0_wp

     ALLOCATE(ocean_state_aux%PgradSalinity_horz_center(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for PgradSalinity_horz_center at cells failed')
     END IF
     ocean_state_aux%PgradSalinity_horz_center       (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%PgradSalinity_horz_center       (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%PgradSalinity_horz_center       (:,:,:)%x(3)=0.0_wp

     
     CALL add_var(ocean_default_list,'DerivTemperature_vert',ocean_state_aux%DerivTemperature_vert_center,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('DerivTemperature_vert','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
        ! & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
        
     CALL add_var(ocean_default_list,'DerivSalinity_vert',ocean_state_aux%DerivSalinity_vert_center,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('DerivSalinity_vert','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
        
     CALL add_var(ocean_default_list,'slopes_squared',ocean_state_aux%slopes_squared,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('slopes_squared','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
        

     CALL add_var(ocean_default_list,'taper function 1',ocean_state_aux%taper_function_1,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('taper function 1','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)

     CALL add_var(ocean_default_list,'taper function 2',ocean_state_aux%taper_function_2,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('taper function 2','','', datatype_flt),&
        & grib2_var(255,255,255,DATATYPE_PACK16,GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_aux"),loutput=.FALSE.)
        
        
        
     ! set all values - incl. last block - of cartesian coordinates to zero (NAG compiler)
     
!      ocean_state_aux%slopes_squared=0.0_wp
    ENDIF
  END SUBROUTINE construct_hydro_ocean_aux
  
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_acc(patch_2d, ocean_state_acc)
    
    TYPE(t_patch),TARGET, INTENT(in)                :: patch_2d
    TYPE(t_hydro_ocean_acc), TARGET,INTENT(inout)   :: ocean_state_acc
    
    ! local variables
    
    INTEGER ::  ist, jtrc
    INTEGER ::  alloc_cell_blocks, nblks_e, nblks_v
    
    INTEGER, PARAMETER :: max_oce_tracer = 2
    CHARACTER(LEN=max_char_length) :: oce_tracer_names(max_oce_tracer),&
      & oce_tracer_units(max_oce_tracer),&
      & oce_tracer_longnames(max_oce_tracer)
    INTEGER :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(LEN=max_char_length) :: var_suffix
    CHARACTER(LEN=*), PARAMETER :: method_name="construct_hydro_ocean_acc"

    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    
    CALL message(TRIM(method_name), '...')
    
    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    
    CALL add_var(ocean_default_list, 'h_acc', ocean_state_acc%h , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('h_acc', 'm', 'surface elevation at cell center', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    CALL add_var(ocean_default_list, 'h_sqr', ocean_state_acc%h_sqr , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('h_acc', 'm', 'surface elevation squared at cell center', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'u_acc', ocean_state_acc%u, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('u_acc','m/s','zonal velocity component', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    ! reconstructed v velocity component
    CALL add_var(ocean_default_list, 'v_acc', ocean_state_acc%v, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('v_acc','m/s','meridional velocity component', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'rhopot_acc', ocean_state_acc%rhopot, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('rhopot_acc','kg/m^3','potential density', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'rho_acc', ocean_state_acc%rho, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('rho_acc','kg/m^3','insitu density', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    
    CALL add_var(ocean_default_list, 'w_acc', ocean_state_acc%w, grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('w_acc','m/s','vertical velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'mass_flx_e_acc', ocean_state_acc%mass_flx_e, grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('mass_flx_e_acc','','mass flux at edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'div_mass_flx_c_acc', ocean_state_acc%div_mass_flx_c, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div_mass_flx_c_acc','','divergence of mass flux', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'edgeFlux_total_acc', ocean_state_acc%edgeFlux_total, &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('edgeFlux_total_acc', 'm*m/s', 'vertically integrated edge flux', DATATYPE_FLT32),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_default"))
    CALL add_var(ocean_default_list, 'u_vint_acc', ocean_state_acc%u_vint , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('u_vint_acc', 'm*m/s', 'barotropic zonal velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
    CALL add_var(ocean_default_list, 'v_vint_acc', ocean_state_acc%v_vint , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('v_vint_acc', 'm*m/s', 'barotropic meridional velocity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(ocean_default_list, 'vort_acc', ocean_state_acc%vort, &
      & grid_unstructured_vert, za_depth_below_sea, &
      & t_cf_var('vort_acc','1/s','vorticity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_vertex),&
      & ldims=(/nproma,n_zlev,nblks_v/),in_group=groups("oce_default"),lrestart_cont=.TRUE.)
    
    ! kinetic energy component
    CALL add_var(ocean_default_list, 'kin_acc', ocean_state_acc%kin, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('kin_acc','J','kinetic energy', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
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
      CALL add_var(ocean_default_list, 'tracers'//TRIM(var_suffix), ocean_state_acc%tracer , &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
      ! Reference to individual tracer, for I/O
      ALLOCATE(ocean_state_acc%tracer_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_default_list, 'tracers'//TRIM(var_suffix),&
          & oce_tracer_names(jtrc),                 &
          & ocean_state_acc%tracer_ptr(jtrc)%p,           &
          & grid_unstructured_cell, za_depth_below_sea,&
          & t_cf_var(oce_tracer_names(jtrc), &
          & oce_tracer_units(jtrc), &
          & oce_tracer_longnames(jtrc), datatype_flt), &
          & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_default", "oce_essentials"))
        
      END DO

    ENDIF ! no_tracer > 0

    ! Accumulates physics
    CALL add_var(ocean_default_list, 'K_veloc_h_acc', ocean_state_acc%k_veloc_h , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('K_veloc_h_acc', 'kg/kg', 'horizontal velocity diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))
 
    CALL add_var(ocean_default_list, 'A_veloc_v_acc', ocean_state_acc%a_veloc_v , grid_unstructured_edge,&
      & za_depth_below_sea_half, &
      & t_cf_var('A_veloc_v_acc', 'kg/kg', 'vertical velocity diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups("oce_physics","oce_default"))
 
!
!   !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL add_var(ocean_default_list, 'K_tracer_h'//TRIM(var_suffix), ocean_state_acc%k_tracer_h , &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('K_tracer_h'//TRIM(var_suffix), '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e,no_tracer+nbgcadv/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.)
      CALL add_var(ocean_default_list, 'A_tracer_v'//TRIM(var_suffix), ocean_state_acc%a_tracer_v , &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('A_tracer_v'//TRIM(var_suffix), '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks,no_tracer+nbgcadv/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.)
 
      ! Reference to individual tracer, for I/O
 
      ALLOCATE(ocean_state_acc%tracer_horz_physics_ptr(no_tracer+nbgcadv))
      ALLOCATE(ocean_state_acc%tracer_vert_physics_ptr(no_tracer+nbgcadv))
      DO jtrc = 1,no_tracer+nbgcadv
        CALL add_ref( ocean_default_list, 'K_tracer_h'//TRIM(var_suffix),&
          & 'K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc))//TRIM(var_suffix),     &
          & ocean_state_acc%tracer_horz_physics_ptr(jtrc)%p,                             &
          & grid_unstructured_edge, za_depth_below_sea,               &
          & t_cf_var('K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc))//TRIM(var_suffix), &
          & 'kg/kg', &
          & TRIM(oce_config%tracer_longnames(jtrc))//'(K_tracer_h'//TRIM(var_suffix)//')', &
          & datatype_flt), &
          & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
          & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))
        CALL add_ref( ocean_default_list, 'A_tracer_v'//TRIM(var_suffix),&
          & 'A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc))//TRIM(var_suffix),     &
          & ocean_state_acc%tracer_vert_physics_ptr(jtrc)%p,                             &
          & grid_unstructured_cell, za_depth_below_sea_half,            &
          & t_cf_var('A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc))//TRIM(var_suffix), &
          & 'kg/kg', &
          & TRIM(oce_config%tracer_longnames(jtrc))//'(A_tracer_v'//TRIM(var_suffix)//')', &
          & datatype_flt), &
          & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_physics"))
 
      END DO
    ENDIF ! no_tracer > 0

!  IF(GMRedi_configuration==GMRedi_combined&
!  &.OR.GMRedi_configuration==GM_only.OR.GMRedi_configuration==Redi_only)THEN
!
!
!    CALL add_var(ocean_default_list, 'k_tracer_isoneutral_acc', ocean_state_acc%k_tracer_isoneutral, &
!       & grid_unstructured_cell, za_depth_below_sea, &
!       & t_cf_var('k_tracer_isoneutral_acc', '', '1:temperature 2:salinity', datatype_flt),&
!       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
!       & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
!       & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
!
!
!     CALL add_var(ocean_default_list, 'k_tracer_dianeutral_acc', ocean_state_acc%k_tracer_dianeutral, &
!       & grid_unstructured_cell, za_depth_below_sea_half, &
!       & t_cf_var('k_tracer_dianeutral_acc', '', '1:temperature 2:salinity', datatype_flt),&
!       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
!       & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
!       & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
!
!
!    CALL add_var(ocean_default_list, 'k_tracer_GM_kappa_acc', ocean_state_acc%k_tracer_GM_kappa, &
!       & grid_unstructured_cell, za_depth_below_sea, &
!       & t_cf_var('k_tracer_GM_kappa_acc_at_cells', '', '1:temperature 2:salinity', datatype_flt),&
!       & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
!       & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
!       & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
!
!   ENDIF
    
  END SUBROUTINE construct_hydro_ocean_acc
  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of auxilliary hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_hydro_ocean_aux(ocean_state_aux)
    
    TYPE(t_hydro_ocean_aux), INTENT(inout)      :: ocean_state_aux
    
    ! local variables
    
    INTEGER :: ist
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_aux'
    
    DEALLOCATE(ocean_state_aux%bc_top_veloc_cc, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond cc failed')
    END IF
    
  END SUBROUTINE destruct_hydro_ocean_aux
  !-------------------------------------------------------------------------
  
  
  
  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
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
!<Optimize:inUse>
  SUBROUTINE construct_patch_3d(patch_3d)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)    :: patch_3d
    
    ! local variables
    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp, n_zlvm!, ie
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_patch_3D'
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF
    
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
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('lsm_c','','3d lsm on cells', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default"),&
      & isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'lsm_e', patch_3d%lsm_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('lsm_e','','3d lsm on edges', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry","oce_default"),&
      & isteptype=tstep_constant)
    ! surface cells
    CALL add_var(ocean_default_list, 'surface_cell_sea_land_mask', patch_3d%surface_cell_sea_land_mask , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('surface_cell_sea_land_mask', '', 'surface_cell_sea_land_mask', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! surface vertices
    CALL add_var(ocean_default_list, 'surface_edge_sea_land_mask', patch_3d%surface_edge_sea_land_mask , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('surface_edge_sea_land_mask', '', 'surface_edge_sea_land_mask', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! surface vertices
    CALL add_var(ocean_default_list, 'surface_vertex_sea_land_mask', patch_3d%surface_vertex_sea_land_mask , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('surface_vertex_sea_land_mask', '', 'surface_vertex_sea_land_mask', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_vertex),&
      & ldims=(/nproma,nblks_v/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! deepest ocean layer in column
    CALL add_var(ocean_default_list, 'dolic_c', patch_3d%p_patch_1d(n_dom)%dolic_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('dolic_c', '', 'dolic_c', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'dolic_e', patch_3d%p_patch_1d(n_dom)%dolic_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('dolic_e', '', 'dolic_e', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'vertex_bottomLevel', patch_3d%p_patch_1d(n_dom)%vertex_bottomLevel , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('vertex_bottomLevel', '', 'vertex_bottomLevel', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_vertex),&
      & ldims=(/nproma,nblks_v/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! 2-dim basins and areas
    CALL add_var(ocean_default_list, 'basin_c', patch_3d%basin_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('basin_c', '', 'basin_c', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'regio_c', patch_3d%regio_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('regio_c', '', 'regio_c', DATATYPE_INT8),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    ! 2-dim bottom and column thickness
    CALL add_var(ocean_default_list, 'bottom_thick_c', patch_3d%bottom_thick_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('bottom_thick_c', 'm', 'bottom_thick_c', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'bottom_thick_e', patch_3d%bottom_thick_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('bottom_thick_e', 'm', 'bottom_thick_e', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'column_thick_c', patch_3d%column_thick_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('column_thick_c', 'm', 'column_thick_c', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'column_thick_e', patch_3d%column_thick_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('column_thick_e', 'm', 'column_thick_e', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! 3-dim real land-sea-mask
    ! cells
    CALL add_var(ocean_default_list, 'wet_c', patch_3d%wet_c , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('wet_c', '', '3d lsm on cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'wet_e', patch_3d%wet_e , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('wet_e', '', '3d lsm on edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry","oce_default"), &
      & isteptype=tstep_constant)
    ! 3-dim real land-sea-mask with zero on halos
    ! cells
    CALL add_var(ocean_default_list, 'wet_halo_zero_c', patch_3d%wet_halo_zero_c , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('wet_c_halo_zero', '', '3d lsm with halo zero on cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'wet_halo_zero_e', patch_3d%wet_halo_zero_e , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('wet_e_halo_zero', '', '3d lsm with halo zero on edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
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
      & t_cf_var('cons thick','m','prism thickness at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'invConstantPrismThickness', patch_3d%p_patch_1d(1)%invConstantPrismThickness, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inv cons thick','m','inverse prism thickness at cells', DATATYPE_FLT32),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_volume', patch_3d%p_patch_1d(1)%prism_volume, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism volume (cells)', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_thick_e', patch_3d%p_patch_1d(n_dom)%prism_thick_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism thickness at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_thick_flat_sfc_c', patch_3d%p_patch_1d(n_dom)%prism_thick_flat_sfc_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_c','m','time independent depth at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_thick_flat_sfc_e', patch_3d%p_patch_1d(n_dom)%prism_thick_flat_sfc_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_c','m','time independent depth at edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inverse prism_thick_c', patch_3d%p_patch_1d(n_dom)%inv_prism_thick_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inverse prism_thick_c','m','time dependent depth at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_center_dist_c', patch_3d%p_patch_1d(n_dom)%prism_center_dist_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_center_dist_c','m','time dependent distance between prism centers', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'constantPrismCenters_Zdistance', &
      & patch_3d%p_patch_1d(n_dom)%constantPrismCenters_Zdistance, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('constantPrismCenters_Zdistance','m','constant distance between prism centers', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'constantPrismCenters_invZdistance', &
      & patch_3d%p_patch_1d(n_dom)%constantPrismCenters_invZdistance, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('constantPrismCenters_invZdistance','m','inverse constant distance between prism centers', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inv_prism_thick_e', patch_3d%p_patch_1d(n_dom)%inv_prism_thick_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('inv_prism_thick_e','m','time dependent inverse thickeness at edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inv_prism_center_dist_c', &
      & patch_3d%p_patch_1d(n_dom)%inv_prism_center_dist_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inv_prism_center_dist_c','1/m','inverse of dist between prism centers at cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'inv_prism_center_dist_e', &
      & patch_3d%p_patch_1d(n_dom)%inv_prism_center_dist_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('inv_prism_center_dist_e','1/m','inverse of dist between prism centers at edges', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'depth_CellMiddle',   &
      & patch_3d%p_patch_1d(n_dom)%depth_CellMiddle,       &
      & grid_unstructured_cell,                             &
      & za_depth_below_sea,                                 &
      & t_cf_var('depth_CellMiddle','m','depth at the middle of the cells', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'depth_CellInterface',   &
      & patch_3d%p_patch_1d(n_dom)%depth_CellInterface,       &
      & grid_unstructured_cell,                                &
      & za_depth_below_sea_half,                               &
      & t_cf_var('depth_CellInterface','m','depth at cell interfaces', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_geometry"),isteptype=tstep_constant)
     
  END SUBROUTINE construct_patch_3d
  
  !-------------------------------------------------------------------------
  !>
  !!
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2011).
  !!
!<Optimize:inUse>
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
  
  
END MODULE mo_ocean_state
