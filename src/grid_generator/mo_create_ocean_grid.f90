!-------------------------------------------------------------------------------------
! mo_create_ocean_grid
!>
!! A collection of grid tools associated with the ocean grid
!!
!! NOTE : negative elevations are below sea level (sea depth)
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2009-12-4
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

!-------------------------------------------------------------------------------------
MODULE mo_create_ocean_grid
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message_text, message, finish
  USE mo_io_units,             ONLY: filename_max, nnml
  USE mo_namelist,             ONLY: position_nml, open_nml, positioned
  USE mo_local_grid
  USE mo_grid_conditions,      ONLY: read_grid_conditions, get_conditional_cells
  USE mo_io_local_grid,        ONLY: read_new_netcdf_grid, write_netcdf_grid, &
    & read_netcdf_cell_elevation, write_netcdf_vertical_strc
  USE mo_grid_toolbox,         ONLY: get_grid_from_cell_list, &
    & smooth_boundaryfrom_cell_list, &
    & get_conditional_list, set_edge_elev_fromcells, get_grid_conditional_cells,       &
    & until_convergence
  USE mo_local_grid_refinement,ONLY: refine_grid_edgebisection
  USE mo_local_grid_geometry,  ONLY: compute_sphere_grid_geometry
  USE mo_local_grid_hierarchy, ONLY: create_grid_hierarchy
  !   USE mo_grid_checktools
  USE mo_zlayers,              ONLY: read_zlayers, get_grid_zlayers, count_dry_wet_cells
  USE mo_local_grid_optimization, ONLY: read_grid_optimization_param, optimize_grid


  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !-------------------------------------------------------------------------
  PUBLIC :: create_ocean_grid, generate_ocean_grid

  !-------------------------------------------------------------------------
  ! input parameters
  ! NOTE : negative elevations are below sea level (sea depth)
  REAL(wp) :: min_sea_depth     = 0.0_wp ! ocean cells shall have a depth < min_sea_depth
  REAL(wp) :: set_sea_depth     = 0.0_wp ! if negative, depth will be set to this
  REAL(wp) :: set_min_sea_depth = 0.0_wp ! if negative, min depth will be set to this
  LOGICAL  :: only_get_sea_land_mask = .false.
  LOGICAL  :: smooth_ocean_boundary = .true.
  LOGICAL  :: write_all_levels = .false.
  CHARACTER(LEN=filename_max) :: input_file  ! the input grid file
  INTEGER :: refine_iterations
  ! the netcdf file where elevations ares stored, and the related field name
  CHARACTER(LEN=filename_max) :: elevation_file, elevation_field
  CHARACTER(LEN=filename_max) :: output_ocean_file=""  ! the ocean part file of the input grid file
  ! the refined grid file (by edge bisection), including the sea-land mask,
  ! and its ocean part file
!  CHARACTER(LEN=filename_max) :: output_refined_file
  CHARACTER(LEN=filename_max) :: output_atmo_file
  INTEGER :: edge_elev_interp_method=min_interpolation
  CHARACTER(LEN=filename_max) :: ocean3d_file_name = "ocean3Dgrid.nc"
  CHARACTER(LEN=filename_max) :: vertical_out_file_name = "ocean3Dvertical.nc"
  CHARACTER(LEN=filename_max) :: vertical_in_file_name   = "vertical.in"
  CHARACTER(LEN=filename_max) :: parameter_file_name
  ! CHARACTER(LEN=filename_max) :: tmp_file_name   = "tmp.nc"
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !   SUBROUTINE read_param(param_file_name)
  !>
  !! Reads the parameters for the ocean grid generation. Private
  SUBROUTINE read_param(param_file_name)
    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: i_status

!     NAMELIST /create_ocean_grid/ input_file, elevation_file, elevation_field, &
!       & output_ocean_file, output_refined_file, output_atmo_file, &
!       & min_sea_depth, edge_elev_interp_method,ocean3d_file_name,          &
!       & vertical_out_file_name, vertical_in_file_name
    NAMELIST /create_ocean_grid/ input_file, elevation_file, elevation_field, &
      & output_atmo_file, min_sea_depth, set_sea_depth, set_min_sea_depth, &
      & edge_elev_interp_method, only_get_sea_land_mask, output_ocean_file, &
      & smooth_ocean_boundary, refine_iterations, write_all_levels 
    !------------------------------------------------------------------------
    ! set default values
    parameter_file_name = param_file_name
    refine_iterations = 1
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('create_ocean_grid',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,create_ocean_grid)
    ELSE
      WRITE(message_text,'(a,a)') " File", param_file_name, " not POSITIONED"
      CALL finish ('read_param', message_text)
    ENDIF
    CLOSE(nnml)

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "-------- Create Ocean Grid ---------"
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "input_file=",  TRIM(input_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i3.3)') "refine_iterations=", refine_iterations
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) "only_get_sea_land_mask=", only_get_sea_land_mask
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) "smooth_ocean_boundary=", smooth_ocean_boundary
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "output_ocean_file=", TRIM(output_ocean_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "output_atmo_file=", TRIM(output_atmo_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "elevation_file=",TRIM(elevation_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "elevation_field=", TRIM(elevation_field)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "ocean3d_file_name=",TRIM(ocean3d_file_name)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "vertical_in_file_name=",TRIM(vertical_in_file_name)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "vertical_out_file_name=",TRIM(vertical_out_file_name)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,F10.2)') "min_sea_depth=", min_sea_depth
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,F10.2)') "set_sea_depth=", set_sea_depth
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,F10.2)') "set_min_sea_depth=", set_min_sea_depth
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,I2)') "edge_elev_interp_method=", edge_elev_interp_method
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))

    CALL read_grid_optimization_param(param_file_name)

  END SUBROUTINE read_param
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE create_ocean_grid(param_file_name)
  !>
  !! The driving method for generating the ocean grid. Public.
  SUBROUTINE create_ocean_grid(param_file_name)
    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    !------------------------------------------------------------------------
    CALL read_param(param_file_name)
    CALL generate_ocean_grid
    
  END SUBROUTINE create_ocean_grid
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !   SUBROUTINE generate_ocean_grid(param_file_name)
  !>
  !! The main method for generating the ocean grid. Public.
  SUBROUTINE generate_ocean_grid()

    INTEGER :: init_grid_id, ocean_grid_id, refined_grid_id
    TYPE(t_integer_list) :: sea_cell_list

    CHARACTER(LEN=filename_max) :: file_name
    INTEGER :: no_of_conditions, i, level

    ! read initial grid and elevation
    init_grid_id = read_new_netcdf_grid(input_file)
    
    NULLIFY(sea_cell_list%value)
    ocean_grid_id=init_grid_id
    !------------------------------------------------------------------------
    ! If use elevation_file
    IF (elevation_file /= "") THEN
      CALL read_netcdf_cell_elevation(init_grid_id, elevation_file, elevation_field)
      sea_cell_list = get_ocean_cell_list_from_depth(init_grid_id)
    ENDIF
    !--------------------------------------------------------------------
    ! check if we have conditions
    no_of_conditions =  read_grid_conditions(parameter_file_name)
    IF ( no_of_conditions > 0) THEN
      IF (elevation_file /= "") THEN
        ! we cannot have both elevation_file and conditional criteria
        CALL finish('generate_ocean_grid',&
          & 'We cannot have both elevation_file and conditional criteria')
      ENDIF
      CALL get_conditional_cells(init_grid_id, sea_cell_list)
    ENDIF
      
    !------------------------------------------------------------------------
    ! sea_cell_list contains now the sea cells, or is empty in case of aqua planet
   IF (set_sea_depth < 0.0_wp) THEN
      CALL grid_set_sea_depth(init_grid_id, set_sea_depth, sea_cell_list)
    ELSE IF (set_min_sea_depth < 0.0_wp) THEN
      CALL grid_set_min_sea_depth(init_grid_id, set_min_sea_depth, sea_cell_list)
    ENDIF
    
    !------------------------------------------------------------------------
    ! get the new sea_cell_list as defined by the cell%elevation value
    IF (ASSOCIATED(sea_cell_list%value)) DEALLOCATE(sea_cell_list%value)
    sea_cell_list = get_ocean_cell_list_from_depth(init_grid_id)
    ! create sea_land_mask
    CALL fill_sea_land_mask_from_list(init_grid_id, sea_cell_list)
    
    !------------------------------------------------------------------------
    IF (ASSOCIATED(sea_cell_list%value)) DEALLOCATE(sea_cell_list%value)
    ! set boundary indices
    ! CALL create_grid_hierarchy(ocean_grid_id)
    ! write the initial ocean_grid_id grid with the sea-land mask
    ! CALL write_netcdf_grid(ocean_grid_id, output_ocean_file)

    !------------------------------------------------------------------------
    ! refine initial grid
    DO i=1, refine_iterations
      refined_grid_id = refine_grid_edgebisection(init_grid_id)
      CALL set_grid_parent_id(refined_grid_id,0)
      CALL delete_grid(init_grid_id)
      CALL optimize_grid(refined_grid_id)
      
      IF (write_all_levels) THEN
        level = get_grid_level(refined_grid_id)
        WRITE(file_name,'(a,i2.2,a)')  TRIM(output_ocean_file), level, ".nc"
        CALL write_netcdf_grid(refined_grid_id, file_name)
      ENDIF
      
      init_grid_id = refined_grid_id
    ENDDO
    !------------------------------------------------------------------------
    ! smooth the grid
    IF (smooth_ocean_boundary) THEN
      sea_cell_list = get_ocean_cell_list_from_depth(init_grid_id, use_smoothing=.true.)
      CALL fill_sea_land_mask_from_list(init_grid_id, sea_cell_list)
    ENDIF
    !------------------------------------------------------------------------
    CALL compute_sphere_grid_geometry(init_grid_id)
    CALL fill_edges_sea_land_mask(init_grid_id)
    IF (output_atmo_file /= "") THEN
      CALL write_netcdf_grid(ocean_grid_id, output_atmo_file)
    ENDIF
    !------------------------------------------------------------------------

    IF (.NOT. only_get_sea_land_mask) THEN
    
      IF ( .NOT. smooth_ocean_boundary) &   ! if  smooth_ocean_boundary the  sea_cell_list is filled      
        & sea_cell_list = get_ocean_cell_list_from_depth(init_grid_id)
      
      CALL grid_set_parents_from(init_grid_id, parents_from_parentpointers)
      ocean_grid_id = get_grid_from_cell_list(init_grid_id, sea_cell_list)
      CALL delete_grid(init_grid_id)
      
    ELSE
      ocean_grid_id=init_grid_id
    ENDIF
    !------------------------------------------------------------------------
    
   ! CALL create_grid_hierarchy(ocean_grid_id)
    CALL set_edge_elev_fromcells(ocean_grid_id, edge_elev_interp_method)
!     CALL fill_edges_sea_land_mask(ocean_grid_id)
    CALL check_ocean_grid_mask(ocean_grid_id)
    CALL set_boundary_edges_tozero(ocean_grid_id)
    IF (write_all_levels) THEN
      level = get_grid_level(ocean_grid_id)
      WRITE(file_name,'(a,i2.2,a)')  TRIM(output_ocean_file), level, ".nc"
    ELSE
      file_name=output_ocean_file
    ENDIF
    CALL write_netcdf_grid(ocean_grid_id, file_name)
    !CALL message('get_ocean_vertical_strc..','')
    !CALL get_ocean_vertical_strc(refined_grid_id)
    CALL delete_grid(ocean_grid_id)

  END SUBROUTINE generate_ocean_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns the ocean part of the in_grid_id. Private
  !! Ocean cells have elevation < min_sea_depth
  FUNCTION get_ocean_cell_list_from_depth(in_grid_id, use_smoothing) result(smooth_sea_cell_list)
    INTEGER, INTENT(in)  :: in_grid_id
    LOGICAL, INTENT(in), OPTIONAL :: use_smoothing
    
    TYPE(t_integer_list) :: smooth_sea_cell_list

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_integer_list) :: sea_cell_list
    TYPE(t_float_list) :: elevation_list
    INTEGER :: no_of_cells
    LOGICAL :: smooth_thelist

    smooth_thelist = .false.
    IF (PRESENT(use_smoothing)) smooth_thelist = use_smoothing
    

    in_grid => get_grid(in_grid_id)
    in_grid%cells%min_sea_depth = min_sea_depth
    elevation_list%list_size    = in_grid%cells%no_of_allocatedcells
    elevation_list%value       => in_grid%cells%elevation
    CALL get_conditional_list(elevation_list, less_than, min_sea_depth, -1, &
      & sea_cell_list)
    
    ! WRITE(0,*) "sea_cell_list size=", sea_cell_list%list_size
    
    ! smooth the boundary cells until converegnce
    IF (smooth_thelist) THEN
      CALL smooth_boundaryfrom_cell_list(in_grid_id, &
        & sea_cell_list, smooth_sea_cell_list, until_convergence)
      DEALLOCATE(sea_cell_list%value)
    !  WRITE(0,*) "smooth_sea_cell_list size=", smooth_sea_cell_list%list_size
    ELSE
      smooth_sea_cell_list%list_size = sea_cell_list%list_size
      smooth_sea_cell_list%value => sea_cell_list%value
    ENDIF
      
    ! inform what we got
    no_of_cells = in_grid%cells%no_of_existcells
    WRITE(message_text,*) "------------ get_ocean_cell_list_from_depth ---------------"
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i8,i8,i8)') "total, land, sea cells=", no_of_cells, &
      & no_of_cells-smooth_sea_cell_list%list_size, smooth_sea_cell_list%list_size
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) "------------------------------------------------------"
    CALL message ('', TRIM(message_text))

  END FUNCTION get_ocean_cell_list_from_depth
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns the ocean part of the in_grid_id. Private
  !! Ocean cells have elevation < min_sea_depth
  FUNCTION get_land_cell_list_from_mask(in_grid_id, use_smoothing) result(smooth_land_cell_list)
    INTEGER, INTENT(in)  :: in_grid_id
    LOGICAL, INTENT(in), OPTIONAL :: use_smoothing

    TYPE(t_integer_list) :: smooth_land_cell_list

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_integer_list) :: land_cell_list
    TYPE(t_integer_list) :: sea_land_mask_list
    INTEGER :: no_of_cells
    LOGICAL :: smooth_thelist

    smooth_thelist = .false.
    IF (PRESENT(use_smoothing)) smooth_thelist = use_smoothing


    in_grid => get_grid(in_grid_id)
    sea_land_mask_list%list_size  = in_grid%cells%no_of_allocatedcells
    sea_land_mask_list%value      => in_grid%cells%sea_land_mask
    CALL get_conditional_list(sea_land_mask_list, greater_than, land_boundary-1, &
      & -1, land_cell_list)

    ! WRITE(0,*) "sea_cell_list size=", sea_cell_list%list_size

    ! smooth the boundary cells until converegnce
    IF (smooth_thelist) THEN
      CALL smooth_boundaryfrom_cell_list(in_grid_id, &
        & land_cell_list, smooth_land_cell_list, until_convergence)
      DEALLOCATE(land_cell_list%value)
    !  WRITE(0,*) "smooth_sea_cell_list size=", smooth_sea_cell_list%list_size
    ELSE
      smooth_land_cell_list%list_size = land_cell_list%list_size
      smooth_land_cell_list%value     => land_cell_list%value
    ENDIF

    ! inform what we got
    no_of_cells = in_grid%cells%no_of_existcells
    WRITE(message_text,*) "------------- get_land_cell_list_from_mask ----------------"
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i8,i8,i8)') "total, land, sea cells=", no_of_cells, &
      & no_of_cells-smooth_land_cell_list%list_size, smooth_land_cell_list%list_size
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) "------------------------------------------------------"
    CALL message ('', TRIM(message_text))

  END FUNCTION get_land_cell_list_from_mask
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE fill_sea_land_mask_from_list(in_grid_id, sea_cell_list)
  !>
  !! Creates the grid\%cells\%sea_land_mask based on the sea_cell_list. Private
  SUBROUTINE fill_sea_land_mask_from_list(in_grid_id, sea_cell_list)
    INTEGER, INTENT(in) :: in_grid_id
    TYPE(t_integer_list)  :: sea_cell_list

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: in_cells
    TYPE(t_grid_edges), POINTER :: in_edges

    INTEGER :: no_of_sea_cells, i, no_of_cells, no_of_edges
    REAL(wp) :: land_min_elevation

    in_grid => get_grid(in_grid_id)
    in_cells => in_grid%cells
    no_of_cells = in_cells%no_of_existcells
    in_edges => in_grid%edges
    no_of_edges = in_edges%no_of_existedges

    ! fill sea_land_mask
    in_cells%sea_land_mask(1:) = land_inner
    no_of_sea_cells = sea_cell_list%list_size
    DO i=1,no_of_sea_cells
      in_cells%sea_land_mask(sea_cell_list%value(i)) = sea_inner
    ENDDO

    ! make sure that cell elevation is consistent with sea_land_mask
    land_min_elevation = MAX(0.0_wp,min_sea_depth)
    DO i=1,no_of_cells
       IF (in_cells%sea_land_mask(i) == land_inner .AND. &
         & in_cells%elevation(i) < land_min_elevation) THEN
         in_cells%elevation(i) = land_min_elevation
       ENDIF
    ENDDO

    CALL fill_edges_sea_land_mask(in_grid_id)

  END SUBROUTINE fill_sea_land_mask_from_list
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !   SUBROUTINE fill_edges_sea_land_mask(in_grid_id)
  !>
  !! Creates the grid\%cells\%sea_land_mask. Private
  SUBROUTINE fill_edges_sea_land_mask(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: in_cells
    TYPE(t_grid_edges), POINTER :: in_edges

    INTEGER :: edge, no_of_cells, no_of_edges
    INTEGER :: cell1, cell2, product_mask

    in_grid => get_grid(in_grid_id)
    in_cells => in_grid%cells
    no_of_cells = in_cells%no_of_existcells
    in_edges => in_grid%edges
    no_of_edges = in_edges%no_of_existedges
            
    DO cell1=1,no_of_cells
      IF (ABS(in_cells%sea_land_mask(cell1)) == 1) &
        in_cells%sea_land_mask(cell1) = in_cells%sea_land_mask(cell1) * 2
    ENDDO
    
    ! fill edge sea_land_mask
    DO edge=1,no_of_edges
       cell1 = in_edges%get_cell_index(edge,1)
       cell2 = in_edges%get_cell_index(edge,2)

       ! check if either of the boundary cells exist
       IF (cell1 == 0) THEN
         in_edges%sea_land_mask(edge) = land_boundary
         IF (in_cells%sea_land_mask(cell2) > 0) THEN
           in_cells%sea_land_mask(cell2) = land_boundary
         ELSE
           in_cells%sea_land_mask(cell2) = sea_boundary
         ENDIF
         in_edges%elevation(edge) = 0.0_wp
         CYCLE
       ENDIF
       IF (cell2 == 0) THEN
         in_edges%sea_land_mask(edge) = land_boundary
         IF (in_cells%sea_land_mask(cell1) > 0) THEN
           in_cells%sea_land_mask(cell1) = land_boundary
         ELSE
           in_cells%sea_land_mask(cell1) = sea_boundary
         ENDIF
         in_edges%elevation(edge) = 0.0_wp
         CYCLE
       ENDIF

       product_mask = in_cells%sea_land_mask(cell1) * in_cells%sea_land_mask(cell2)
       IF (product_mask < 0) THEN
         ! we have a boundary between sea and land
         in_edges%sea_land_mask(edge) = land_boundary
         IF (in_cells%sea_land_mask(cell1) < 0) THEN
            in_cells%sea_land_mask(cell1) = sea_boundary
            in_cells%sea_land_mask(cell2) = land_boundary
!             in_edges%elevation(edge) = in_cells%elevation(cell2)
            in_edges%elevation(edge) = 0.0_wp
         ELSE
            in_cells%sea_land_mask(cell2) = sea_boundary
            in_cells%sea_land_mask(cell1) = land_boundary
!             in_edges%elevation(edge) = in_cells%elevation(cell1)
            in_edges%elevation(edge) = 0.0_wp
         ENDIF
         CYCLE
       ENDIF

       ! both cells are of the same type
       IF (in_cells%sea_land_mask(cell1) > 0) THEN
         in_edges%sea_land_mask(edge) = land_inner
       ELSE
         in_edges%sea_land_mask(edge) = sea_inner
       ENDIF


    ENDDO !edge=1,no_of_edges

  END SUBROUTINE fill_edges_sea_land_mask
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE set_boundary_edges_tozero(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_edges), POINTER :: in_edges

    INTEGER :: no_of_edges

    in_grid => get_grid(in_grid_id)
    in_edges => in_grid%edges
    no_of_edges = in_edges%no_of_existedges
            
    WHERE(ABS(in_edges%sea_land_mask(:)) == land_boundary)
      in_edges%sea_land_mask(:) = 0
    END WHERE
    
  END SUBROUTINE set_boundary_edges_tozero
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Returns the ocean part of the in_grid_id. Private
  !! Ocean cells have elevation < min_sea_depth
  SUBROUTINE check_ocean_grid_mask(in_grid_id)
    INTEGER, INTENT(in)  :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: in_cells
    TYPE(t_grid_edges), POINTER :: in_edges
    TYPE(t_float_list) :: elevation_list

    INTEGER :: i, cell, no_of_cells, edge, cell1, cell2
    REAL(wp) :: check_min_sea_depth
    LOGICAL :: edge_is_boundary

    ! get sea-land for the original file
    in_grid => get_grid(in_grid_id)
    in_cells => in_grid%cells
    no_of_cells = in_cells%no_of_existcells
    in_edges => in_grid%edges
    
    check_min_sea_depth = in_grid%cells%min_sea_depth 
    elevation_list%list_size    = in_grid%cells%no_of_allocatedcells
    elevation_list%value       => in_grid%cells%elevation

    ! check cell sea land mask
    DO cell=1, no_of_cells
    !  write (*,*) in_cells%elevation(cell), check_min_sea_depth, in_cells%sea_land_mask(cell)
      IF (in_cells%elevation(cell) < check_min_sea_depth) THEN
        IF (in_cells%sea_land_mask(cell) > sea_boundary) THEN
          CALL finish('check_ocean_grid_mask',&
            & 'sea_land_mask(cell) > sea_boundary')
        ELSE
          CYCLE
        ENDIF
      ELSE
        IF (in_cells%sea_land_mask(cell) < land_boundary) THEN
          CALL finish('check_ocean_grid_mask',&
            & 'sea_land_mask(cell)  < land_boundary')
        ENDIF
      ENDIF

      edge_is_boundary = .false.
      DO i=1, in_cells%max_no_of_vertices
        edge=in_cells%get_edge_index(cell,i)
        IF (edge /= 0) THEN
          edge_is_boundary = (ABS(in_edges%sea_land_mask(edge)) == 1)
        ENDIF
        IF (edge_is_boundary) EXIT
      ENDDO

      IF (edge_is_boundary) THEN
        IF (ABS(in_cells%sea_land_mask(cell)) /= 1) THEN
          CALL finish('check_ocean_grid_mask',&
            & 'ABS(in_cells%sea_land_mask(cell)) /= 1')
        ENDIF
      ELSE
        IF (ABS(in_cells%sea_land_mask(cell)) /= 2) THEN
          write(0,*) cell, in_cells%sea_land_mask(cell),&
             in_edges%sea_land_mask(in_cells%get_edge_index(cell,1)), &
             in_edges%sea_land_mask(in_cells%get_edge_index(cell,2)), &
             in_edges%sea_land_mask(in_cells%get_edge_index(cell,3))
          CALL finish('check_ocean_grid_mask',&
            & 'ABS(in_cells%sea_land_mask(cell)) /= 2')
        ENDIF
      ENDIF        
      
    ENDDO

    ! check edge sea land mask
    DO edge=1, in_edges%no_of_existedges
    !  write (*,*) in_cells%elevation(cell), check_min_sea_depth, in_cells%sea_land_mask(cell)
      cell1 = in_edges%get_cell_index(edge, 1)
      cell2 = in_edges%get_cell_index(edge, 2)
      IF (ABS(in_edges%sea_land_mask(edge)) == 2) THEN
        
        IF (in_cells%sea_land_mask(cell1) * in_edges%sea_land_mask(edge) <= 0 .OR.&
            in_cells%sea_land_mask(cell2) * in_edges%sea_land_mask(edge) <= 0) THEN
            WRITE(0,*) in_edges%sea_land_mask(edge), in_cells%sea_land_mask(cell1),&
               in_cells%sea_land_mask(cell2)
            CALL finish('check_ocean_grid_mask',&
            & 'sea_land_mask(cell) * sea_land_mask(edge) <= 0')
        ENDIF
      ELSE
        IF (cell1 == 0 .OR. cell2 == 0) THEN
          cell = MAX(cell1,cell2)
          IF (ABS(in_cells%sea_land_mask(cell) * in_edges%sea_land_mask(edge)) /= 1) THEN
            WRITE(0,*) in_edges%sea_land_mask(edge), in_cells%sea_land_mask(cell)
            CALL finish('check_ocean_grid_mask',&
            & 'sea_land_mask(cell) * sea_land_mask(edge) /= 1')
          ENDIF
        ELSE
          IF (in_cells%sea_land_mask(cell1) * in_cells%sea_land_mask(cell2) >= 0) THEN
            WRITE(0,*) in_edges%sea_land_mask(edge), in_cells%sea_land_mask(cell1), &
              in_cells%sea_land_mask(cell2)
            CALL finish('check_ocean_grid_mask',&
            & 'sea_land_mask(cell1) sea_land_mask(cell2) >= 0')
          ENDIF        
        ENDIF
        
      ENDIF !ABS(in_edges%sea_land_mask(edge)) == 2
        
    ENDDO


  END SUBROUTINE check_ocean_grid_mask
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE get_ocean_vertical_strc(in_grid_id)
  !>
  !! Gets the vertical parameters of in_grid_id. Private
  SUBROUTINE get_ocean_vertical_strc(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_vertical_ocean_structure), TARGET :: vertical

    !INTEGER :: ocean_3Dgrid_id

    ! read vertical structure
    CALL read_zlayers(vertical, vertical_in_file_name)
    CALL get_grid_zlayers(in_grid_id, vertical)

    ! create 3d grid
    !ocean_3Dgrid_id = create_ocean_3D_zgrid(in_grid_id)
    !CALL check_grid_edges(ocean_3Dgrid_id)
    !CALL write_netcdf_grid(ocean_3Dgrid_id, ocean3d_file_name)
    !CALL delete_grid(ocean_3Dgrid_id)

   END SUBROUTINE get_ocean_vertical_strc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! INTEGER FUNCTION create_ocean_3D_zgrid(in_grid_id, vertical)
  !>
  !! Returns the 3D grid from the in_grid_id and vertical. Private
  INTEGER FUNCTION create_ocean_3D_zgrid(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    INTEGER :: out_grid_id
    TYPE(t_grid), POINTER :: in_grid, out_grid
    TYPE(t_grid_cells), POINTER :: cells, out_cells
    ! TYPE(grid_edges), POINTER :: in_edges
    TYPE(t_vertical_ocean_structure), POINTER :: vertical

!    REAL(wp) :: cell_depth
    INTEGER :: no_of_cells, no_of_out_cells
    INTEGER :: no_of_zlayers, no_of_columns
    INTEGER :: i,j,start_cell

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    no_of_cells = cells%no_of_existcells
    ! in_edges => in_grid%edges
    vertical => in_grid%vertical_structure
    no_of_zlayers = vertical%no_of_zlayers

    CALL count_dry_wet_cells(in_grid_id)

    !------------------------------------------
    !  create out_grid
    out_grid_id = new_grid()
    out_grid    => get_grid(out_grid_id)
    out_grid%ncells = in_grid%ncells * no_of_zlayers
    out_grid%nedges = in_grid%nedges * no_of_zlayers
    out_grid%nverts = in_grid%nverts * no_of_zlayers
    out_grid%cells%max_no_of_vertices = in_grid%cells%max_no_of_vertices
    out_grid%verts%max_connectivity   = in_grid%verts%max_connectivity
    CALL allocate_grid_object(out_grid_id)
    out_grid%vertical_structure => vertical
    out_cells => out_grid%cells

    ! make sure we get the right parent pointer
    CALL set_grid_creation(in_grid_id, undefined)
    ! add the cells for each layer to out_grid
    DO i = 0, no_of_zlayers-1
      start_cell = out_cells%no_of_existcells + 1
      out_grid_id = get_grid_conditional_cells(in_grid_id, &
        &  in_grid%cells%no_of_zlayers, greater_than, i, out_grid_id)
      ! set the z layer which the 3D cells belong
      out_cells%no_of_zlayers(start_cell:out_cells%no_of_existcells) = i + 1
    ENDDO

    CALL grid_set_allocated_eq_exist(out_grid_id)
    CALL create_grid_hierarchy(out_grid_id)
    out_grid    => get_grid(out_grid_id)
    out_cells => out_grid%cells

    !------------------------------------------
    ! create the column structure
    no_of_columns = no_of_cells ! = to the surface no_of_cells
    CALL allocate_vertical_columns(vertical, no_of_columns)
    no_of_out_cells = out_cells%no_of_existcells
    vertical%column_cell_index(:,:) = 0
    DO i=1,no_of_out_cells
      vertical%column_cell_index(out_cells%no_of_zlayers(i),out_cells%parent_index(i)) = i
    ENDDO
    ! fill the column size
    vertical%column_size(:) = -1
    DO i=1,no_of_columns
      DO j=1,no_of_zlayers
        IF (vertical%column_cell_index(j,i) == 0) THEN
          vertical%column_size(i) = j-1
          EXIT
        ENDIF
      ENDDO
      IF (vertical%column_size(i) == -1) &
        & vertical%column_size(i) = no_of_zlayers
    ENDDO
    !------------------------------------------

    CALL check_3dzGrid(in_grid_id, out_grid_id) ! run some checks

    CALL write_netcdf_vertical_strc(vertical, vertical_out_file_name)
    CALL deallocate_vertical(vertical)

    !------------------------------------------
    !return
    create_ocean_3D_zgrid = out_grid_id

  END FUNCTION create_ocean_3D_zgrid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE check_3dzGrid(in_grid_id, grid_3D_id)
  !>
  !! Checks the consistency of the ocean grid vertical structure.
  ! Private
  SUBROUTINE check_3dzGrid(in_grid_id, grid_3D_id)
    INTEGER, INTENT(in) :: in_grid_id
    INTEGER, INTENT(in) :: grid_3D_id

    TYPE(t_grid), POINTER :: in_grid, grid_3D
    TYPE(t_grid_cells), POINTER :: in_cells, cells_3D

    TYPE(t_vertical_ocean_structure), POINTER :: vertical
    INTEGER :: no_of_in_cells, no_of_3D_cells
    INTEGER :: no_of_zlayers, no_of_columns
    INTEGER :: i,j,cell_no

    in_grid => get_grid(in_grid_id)
    in_cells => in_grid%cells
    no_of_in_cells = in_cells%no_of_existcells

    grid_3D => get_grid(grid_3D_id)
    cells_3D => grid_3D%cells
    no_of_3D_cells = cells_3D%no_of_existcells
    vertical => grid_3D%vertical_structure
    no_of_zlayers = vertical%no_of_zlayers
    no_of_columns = vertical%no_of_columns

    IF (no_of_columns /= no_of_in_cells) THEN
      CALL finish('check_3dzGrid','no_of_columns /= no_of_in_cells')
    ENDIF
    DO i=1,no_of_columns
     ! PRINT *,i,'column_size=',vertical%column_size(i)

      DO j=1,vertical%column_size(i)
        cell_no =vertical%column_cell_index(j,i)
      !  PRINT *, cell_no, cells_3D%parent_index(cell_no)
        IF (cells_3D%parent_index(cell_no) /= i) THEN
          CALL finish('check_3dzGrid','cells_3D%parent_index(cell_no) /= i')
        ENDIF
        IF (cells_3D%no_of_zlayers(cell_no) /= j) THEN
          CALL finish('check_3dzGrid','cells_3D%no_of_zlayers(cell_no) /= j')
        ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE check_3dzGrid
  !-------------------------------------------------------------------------





END MODULE mo_create_ocean_grid

