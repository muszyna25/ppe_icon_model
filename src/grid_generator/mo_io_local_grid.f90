!>
!!   Contains the IO subroutines for the grid.
!!
!! The subroutines
!!   <i>output_grid</i> and <i>input_grid</i> produce netCDF gridmap files
!!   analogous to those produced by the previous grid generator.
!!
!! @par Revision History
!!  Developed  by Luca Bonaventura and Thomas Heinze (2005).
!!  Modifications by Thomas Heinze (2005-08-29):
!!  - added VERSION CONTROL
!!  - adapted code according to latest programming guide
!!  Modifications by Th.Heinze, DWD (2006-10-25):
!!  - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!    triangle_info, grid_cells, grid_edges and grid_vertices
!!  Modification by P Ripodas, DWD, (2007-02):
!!  - input_grid and output_grid are modified according to the new array
!!  system_orientation at the grid edges.
!!  Modification by Almut Gassmann, MPI-M (2007-03)
!!  - reorganized optimization to include small circle stuff (C-grid option and
!!    equal area subdivision
!!  - sum over dual edge lengths for ave. resolution (great circles) statistics
!!  - directory structure adaptations (mo_io_units now in "shared_general"
!!    and uses different parameters NAMELIST_FILE)
!!  Modification by Almut Gassmann, MPI-M (2007-04)
!!  - modified treatment of area_edge and removed unused edges%neighbor_index
!! Modification by Luis Kornblueh (2008-04-21):
!! - rearrange for use of netCDF and naming scheme for map files to explain
!!   content as are root division and level. The resulting file is architecture
!!   independent and needs to be generated once only.
!! Guenther Zaengl, DWD, 2008-12-08
!!  - add option for stopping grid optimization at a certain level
!!    (needed for combining spring dynamics with nesting)
!! Almut Gassmann, MPI_M, (2009-01-15)
!!  - Cleanings. Adding of mid_edge_len for cells. Adding dummy values for
!!    refinement variables so that the patch generator has not to be run
!!    if the grid is needed on one level. Grid is no longer on the unit sphere
!!    but given already in physical coordinates. Added double periodic option.
!!    Adjust statistics subroutine with respect to the last two points.
!! Almut Gassmann, MPI-M (2009-01-16)
!!  - mid_edge_len seems to be impractical, use rather edge_vert_length and
!!    edge_cell_length
!! Marco Giorgetta, MPI-M (2009-01-27)
!!  - added attribute 'coordinates' for variables 'cell_area' and 'dual_area'
!! Almut Gassmann, MPI-M (2009-03-02)
!!  - beta_spring parameter may be set via namelist input
!! Modification by Marco Giorgetta (2009-05-08)
!! - set default resolution to R2B04
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
!! $Id: n/a$
!!
MODULE mo_io_local_grid
#include "grid_definitions.inc"

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: find_next_free_unit, filename_max
  USE mo_exception,          ONLY: message_text, message, finish, warning
  USE mo_math_constants,     ONLY: pi,pi_2!, eps
!   USE mo_physical_constants, ONLY: re
  USE mo_local_grid
  USE mo_grid_geometry_info, ONLY: sphere_geometry, planar_torus_geometry, t_grid_geometry_info, &
    & read_geometry_info
  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell, &
    & min_rlvert, max_rlvert, min_rledge, max_rledge
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_util_uuid,          ONLY: t_uuid, uuid_generate, &
       &                           uuid_unparse, uuid_string_length

  IMPLICIT NONE

#define lat_is_pole(lat) (ABS(ABS(lat) - pi_2) < 1e-5_wp)
!#define lat_is_pole(lat) (ABS(ABS(lat) - pi_2) < eps )

!   PUBLIC nf

  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_netcdf_grid, write_netcdf_grid
  PUBLIC :: read_netcdf_cell_elevation
  PUBLIC :: write_netcdf_vertical_strc, read_netcdf_vertical_strc
  PUBLIC :: read_no_of_subgrids, read_new_netcdf_grid
  !--------------------------------------------------------------------

CONTAINS

  !--------------------------------------------------------------------
  SUBROUTINE nf(return_status)
    INTEGER, INTENT(in) :: return_status

    IF (return_status /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(return_status))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE read_netcdf_cell_elevation(grid_id, file_name, elevation_field)
    INTEGER, INTENT(in) :: grid_id
    CHARACTER(LEN=filename_max), INTENT(in) :: file_name
    CHARACTER(LEN=filename_max), INTENT(in) :: elevation_field

    TYPE(t_grid), POINTER :: grid_obj
    !   INTEGER:: no_of_cells
    INTEGER :: ncid, varid

    WRITE(message_text,'(a,a,a,a)') 'Read elevation file ', TRIM(file_name), &
      & "  elevation_field=", TRIM(elevation_field)
    CALL message ('', TRIM(message_text))
    !-------------------------------------------------------------------------

    grid_obj => get_grid(grid_id)


    CALL nf(nf_open(TRIM(file_name), nf_nowrite, ncid))


    CALL nf(nf_inq_varid(ncid, TRIM(elevation_field), varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%elevation(1:)))

    CALL nf(nf_close(ncid))

  END SUBROUTINE read_netcdf_cell_elevation
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_netcdf_vertical_strc(vertical_strc, write_file_name)
    TYPE(t_vertical_ocean_structure), INTENT(in) :: vertical_strc
    CHARACTER(LEN=filename_max), INTENT(in) :: write_file_name

    INTEGER :: old_mode
    INTEGER :: ncid
    INTEGER :: dimids(2)

    INTEGER :: dim_zlayers, dim_columns
    INTEGER :: nc_layer_thicknes, nc_layer_bed, nc_layer_middle
    INTEGER :: nc_column_size, nc_column_cell_index
    !----------------------------------------------------------------------

    WRITE(message_text,'(a,a)') 'Write vertical_strc file: ', TRIM(write_file_name)
    CALL message ('', TRIM(message_text))

    !----------------------------------------------------------------------
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(write_file_name), nf_clobber, ncid))
    CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))
    !----------------------------------------------------------------------
    ! Dimensions
    CALL nf(nf_def_dim(ncid, 'no_of_zlayers',   vertical_strc%no_of_zlayers, dim_zlayers))
    CALL nf(nf_def_dim(ncid, 'no_of_columns',   vertical_strc%no_of_columns, dim_columns))
    !----------------------------------------------------------------------
    ! Define variables
    ! the grid column info
 !   INTEGER, POINTER :: column_cell_index(:,:) ! (no_of_zlayers,no_of_columns)
 !   INTEGER, POINTER :: column_size(:) ! (no_of_columns)

    CALL nf(nf_def_var(ncid, 'layer_thickness', nf_double, 1, dim_zlayers,&
      & nc_layer_thicknes))
    CALL nf(nf_put_att_text(ncid, nc_layer_thicknes, 'long_name', 40, &
      & 'layer thicknes'))
    CALL nf(nf_put_att_text(ncid, nc_layer_thicknes, 'units', 2, 'm2'))

    CALL nf(nf_def_var(ncid, 'layer_bed', nf_double, 1, dim_zlayers,&
      & nc_layer_bed))
    CALL nf(nf_put_att_text(ncid, nc_layer_bed, 'long_name', 40, &
      & 'layer bed depth'))
    CALL nf(nf_put_att_text(ncid, nc_layer_bed, 'units', 2, 'm'))

    CALL nf(nf_def_var(ncid, 'layer_middle', nf_double, 1, dim_zlayers,&
      & nc_layer_middle))
    CALL nf(nf_put_att_text(ncid, nc_layer_middle, 'long_name', 40, &
      & 'layer middle depth'))
    CALL nf(nf_put_att_text(ncid, nc_layer_bed, 'units', 2, 'm'))

    CALL nf(nf_def_var(ncid, 'column_size', nf_int, 1, dim_columns, &
      & nc_column_size))
    CALL nf(nf_put_att_text(ncid, nc_column_size, 'long_name', 40, &
      & 'column number of layers '))

    dimids = (/ dim_zlayers, dim_columns /)
    CALL nf(nf_def_var(ncid, 'column_cell_index', nf_int, 2, dimids, nc_column_cell_index))
    CALL nf(nf_put_att_text(ncid, nc_column_cell_index, 'long_name', 40, &
      & 'cell index n (layer,column)'))

    CALL nf(nf_enddef(ncid))
    !----------------------------------------------------------------------
    ! Write variables
    CALL nf(nf_put_var_double(ncid, nc_layer_thicknes, vertical_strc%layer_thicknes(:)))
    CALL nf(nf_put_var_double(ncid, nc_layer_bed, vertical_strc%layer_bed(:)))
    CALL nf(nf_put_var_double(ncid, nc_layer_middle, vertical_strc%layer_middle(:)))
    CALL nf(nf_put_var_int(ncid, nc_column_size, vertical_strc%column_size(:)))
    CALL nf(nf_put_var_int(ncid, nc_column_cell_index, &
      & vertical_strc%column_cell_index(:,:)))
    !------------------------------------------------------------------------
    CALL nf(nf_close(ncid))

  END SUBROUTINE write_netcdf_vertical_strc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE read_netcdf_vertical_strc(vertical_strc, file_name)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_strc
    CHARACTER(LEN=filename_max), INTENT(in) :: file_name

    INTEGER :: ncid, dimid, varid

    !-------------------------------------------------------------------------
    WRITE(message_text,'(a,a)') 'Read vertical_strc file ', TRIM(file_name)
    CALL message ('', TRIM(message_text))

    CALL nf(nf_open(TRIM(file_name), nf_nowrite, ncid))

    ! Dimensions
    CALL nf(nf_inq_dimid(ncid, 'no_of_zlayers', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, vertical_strc%no_of_zlayers))
    CALL nf(nf_inq_dimid(ncid, 'no_of_columns', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, vertical_strc%no_of_columns))

    CALL allocate_vertical(vertical_strc)

    ! read data
    CALL nf(nf_inq_varid(ncid, 'layer_thickness', varid))
    CALL nf(nf_get_var_double   (ncid, varid, vertical_strc%layer_thicknes(:) ))
    CALL nf(nf_inq_varid(ncid, 'layer_bed', varid))
    CALL nf(nf_get_var_double   (ncid, varid, vertical_strc%layer_bed(:) ))
    CALL nf(nf_inq_varid(ncid, 'layer_middle', varid))
    CALL nf(nf_get_var_double   (ncid, varid, vertical_strc%layer_middle(:) ))
    CALL nf(nf_inq_varid(ncid, 'column_size', varid))
    CALL nf(nf_get_var_int   (ncid, varid, vertical_strc%column_size(:) ))
    CALL nf(nf_inq_varid(ncid, 'column_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, vertical_strc%column_cell_index(:,:) ))

    CALL nf(nf_close(ncid))

  END SUBROUTINE read_netcdf_vertical_strc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION read_no_of_subgrids(file_name) result(no_of_subgrids)
    CHARACTER(LEN=filename_max), INTENT(in) :: file_name

    INTEGER :: ncid, varid
    INTEGER :: netcd_status

    !-------------------------------------------------------------------------
    !CALL message('read_no_of_subgrids ', TRIM(file_name))
    CALL nf(nf_open(TRIM(file_name), nf_nowrite, ncid))

    netcd_status = nf_inq_attid(ncid, nf_global, 'no_of_subgrids', varid)
    IF (netcd_status == nf_noerr) THEN
      CALL nf(nf_get_att_int(ncid, nf_global, 'no_of_subgrids', no_of_subgrids))
    ELSE
      CALL message('read_no_of_subgrids ', 'no_of_subgrids not defined!')
       no_of_subgrids = 1
    ENDIF
    
!     CALL nf(nf_close(ncid))
    
  END FUNCTION read_no_of_subgrids
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  INTEGER FUNCTION read_new_netcdf_grid(file_name, read_grid_ids) result(new_grid_id)
    CHARACTER(LEN=filename_max), OPTIONAL, INTENT(in) :: file_name
    LOGICAL, OPTIONAL, INTENT(in) :: read_grid_ids
    
    new_grid_id = new_grid()
    IF (PRESENT(file_name)) &
      CALL set_grid_filename(new_grid_id, file_name)
    
    IF (PRESENT(read_grid_ids)) THEN
      CALL read_netcdf_grid(new_grid_id, read_grid_ids=read_grid_ids)
    ELSE
      CALL read_netcdf_grid(new_grid_id)
    ENDIF
    
  END FUNCTION read_new_netcdf_grid
  !-------------------------------------------------------------------------
          
 
  !-------------------------------------------------------------------------
  SUBROUTINE read_netcdf_grid(grid_id, file_name, read_grid_ids)
    INTEGER, INTENT(in) :: grid_id
    CHARACTER(LEN=filename_max), OPTIONAL, INTENT(in) :: file_name
    LOGICAL, OPTIONAL, INTENT(in) :: read_grid_ids
    
    TYPE(t_grid), POINTER :: grid_obj
    ! INTEGER:: no_of_cells, no_of_edges, no_of_verts
    INTEGER :: ncid, dimid, varid
    INTEGER :: no_of_cells, no_of_edges, no_of_verts
    INTEGER :: max_vert_connect, max_cell_vertices
    INTEGER, POINTER :: tmp_index(:,:)
    REAL(wp), POINTER :: tmp_real(:,:)
    INTEGER :: i,return_status, start_subgrid_id

    CHARACTER(*), PARAMETER :: method_name = "mo_io_local_grid:read_netcdf_grid"
    
    grid_obj => get_grid(grid_id)
    IF (PRESENT(file_name)) &
      CALL set_grid_filename(grid_id, file_name)
    !-------------------------------------------------------------------------
    WRITE(message_text,'(a,a)') 'Read gridmap file ', &
      & TRIM(grid_obj%file_name)
    CALL message ('', TRIM(message_text))

    CALL nf(nf_open(TRIM(grid_obj%file_name), nf_nowrite, ncid))

    IF (grid_obj%level == undefined) THEN
      CALL nf(nf_get_att_int(ncid, nf_global, 'grid_level', grid_obj%level))
    ENDIF

    CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, grid_obj%ncells))
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, grid_obj%nedges))
    CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, grid_obj%nverts))
    CALL nf(nf_inq_dimid(ncid, 'nv', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, grid_obj%cells%max_no_of_vertices))
    CALL nf(nf_inq_dimid(ncid, 'ne', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, grid_obj%verts%max_connectivity))
    
    return_status = nf_inq_attid(ncid, nf_global, 'no_of_subgrids', varid)
    IF (return_status == nf_noerr) THEN
       CALL nf(nf_get_att_int(ncid, nf_global, 'no_of_subgrids', grid_obj%no_of_subgrids))
    ELSE
       grid_obj%no_of_subgrids = 1
    ENDIF
    return_status = nf_inq_attid(ncid, nf_global, 'start_subgrid_id', varid)
    IF (return_status == nf_noerr) THEN
       CALL nf(nf_get_att_int(ncid, nf_global, 'start_subgrid_id', start_subgrid_id))
    ELSE
       start_subgrid_id = 0
    ENDIF
    
    return_status = read_geometry_info(ncid, grid_obj%geometry_info)
    IF (return_status /= 0) &
      CALL warning(method_name, "Cannot not read geometry_info")
        
    print *, 'Read ', grid_obj%ncells, ' cells...'
    CALL allocate_grid_object(grid_id)

    CALL grid_set_exist_eq_allocated(grid_id)
    no_of_cells = grid_obj%ncells
    no_of_edges = grid_obj%nedges
    no_of_verts = grid_obj%nverts
    max_vert_connect  = grid_obj%verts%max_connectivity
    max_cell_vertices = grid_obj%cells%max_no_of_vertices

    ! read cell info
    CALL nf(nf_inq_varid(ncid, 'cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%idx))

    ! read cell geometry
    CALL nf(nf_inq_varid(ncid, 'lon_cell_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%center(:)%lon))
    CALL nf(nf_inq_varid(ncid, 'lat_cell_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%center(:)%lat))
    CALL nf(nf_inq_varid(ncid, 'cell_area_p', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%area))

    return_status = nf_inq_varid(ncid, 'cell_elevation', varid)
    IF (return_status == nf_noerr) THEN
       CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%elevation))
    ENDIF

    CALL nf(nf_inq_varid(ncid, 'cell_sea_land_mask', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%sea_land_mask))
      
    ! read cell connectivity
    ! allocate tmp arrays for reading
    ALLOCATE(tmp_index(no_of_cells,max_cell_vertices),stat=i)
    IF (i > 0) THEN
      CALL finish (method_name, 'failed to ALLOCATE(tmp_index)')
    ENDIF

    CALL nf(nf_inq_varid(ncid, 'edge_of_cell', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_cell_vertices
      grid_obj%cells%get_edge_index(:,i) = tmp_index(:,i)
    ENDDO

    CALL nf(nf_inq_varid(ncid, 'orientation_of_normal', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_cell_vertices
      grid_obj%cells%get_edge_orient(:,i) = tmp_index(:,i)
    ENDDO

    CALL nf(nf_inq_varid(ncid, 'vertex_of_cell', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_cell_vertices
      grid_obj%cells%get_vertex_index(:,i) = tmp_index(:,i)
    ENDDO

    CALL nf(nf_inq_varid(ncid, 'neighbor_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_cell_vertices
      grid_obj%cells%get_neighbor_index(:,i) = tmp_index(:,i)
    ENDDO

    DEALLOCATE(tmp_index)
    ! read cell nest hierarchy info
    CALL nf(nf_inq_varid(ncid, 'parent_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%parent_index))
    CALL nf(nf_inq_varid(ncid, 'parent_cell_type', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%parent_child_type))
    CALL nf(nf_inq_varid(ncid, 'child_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%child_index))
    CALL nf(nf_inq_varid(ncid, 'start_idx_c', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%start_idx))
    CALL nf(nf_inq_varid(ncid, 'end_idx_c', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%end_idx))
    CALL nf(nf_inq_varid(ncid, 'refin_c_ctrl', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%refin_ctrl))

    ! read edge info
    CALL nf(nf_inq_varid(ncid, 'edge_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%idx))
    ! read edge geometry part
    CALL nf(nf_inq_varid(ncid, 'lon_edge_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%center(:)%lon))
    CALL nf(nf_inq_varid(ncid, 'lat_edge_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%center(:)%lat))
    CALL nf(nf_inq_varid(ncid, 'edge_length', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%primal_edge_length))
    CALL nf(nf_inq_varid(ncid, 'dual_edge_length', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%dual_edge_length))
    CALL nf(nf_inq_varid(ncid, 'zonal_normal_primal_edge', varid))
    
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%primal_normal(:)%v1))
    CALL nf(nf_inq_varid(ncid, 'meridional_normal_primal_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%primal_normal(:)%v2))
    CALL nf(nf_inq_varid(ncid, 'zonal_normal_dual_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%dual_normal(:)%v1))
    CALL nf(nf_inq_varid(ncid, 'meridional_normal_dual_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%dual_normal(:)%v2))
    
    CALL nf(nf_inq_varid(ncid, 'edge_system_orientation', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%system_orientation))
    CALL nf(nf_inq_varid(ncid, 'edgequad_area', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%quad_area))

    return_status = nf_inq_varid(ncid, 'edge_elevation', varid)
    IF (return_status == nf_noerr) THEN
       CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%elevation))
    ENDIF

    ! allocate tmp arrays for reading
    ALLOCATE(tmp_real(no_of_edges,2),stat=i)
    IF (i > 0) THEN
      CALL finish (method_name, 'failed to ALLOCATE(tmp_real)')
    ENDIF
    CALL nf(nf_inq_varid(ncid, 'edge_vert_distance', varid))
    CALL nf(nf_get_var_double   (ncid, varid, tmp_real))
    grid_obj%edges%get_edge_vert_length(:,1) =tmp_real(:,1)
    grid_obj%edges%get_edge_vert_length(:,2) =tmp_real(:,2)

    CALL nf(nf_inq_varid(ncid, 'edge_cell_distance', varid))
    CALL nf(nf_get_var_double   (ncid, varid, tmp_real))
    grid_obj%edges%get_edge_cell_length(:,1) =tmp_real(:,1)
    grid_obj%edges%get_edge_cell_length(:,2) =tmp_real(:,2)
    DEALLOCATE(tmp_real)

    ! read edges connectivity
    ! allocate tmp arrays for reading
    ALLOCATE(tmp_index(no_of_edges,2),stat=i)
    IF (i > 0) THEN
      CALL finish (method_name, 'failed to ALLOCATE(tmp_index)')
    ENDIF
    CALL nf(nf_inq_varid(ncid, 'edge_vertices', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    grid_obj%edges%get_vertex_index(:,1) =tmp_index(:,1)
    grid_obj%edges%get_vertex_index(:,2) =tmp_index(:,2)

    CALL nf(nf_inq_varid(ncid, 'adjacent_cell_of_edge', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    grid_obj%edges%get_cell_index(:,1) =tmp_index(:,1)
    grid_obj%edges%get_cell_index(:,2) =tmp_index(:,2)
    DEALLOCATE(tmp_index)

    ! read edges nest hierarchy info
    CALL nf(nf_inq_varid(ncid, 'child_edge_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%child_index))
    CALL nf(nf_inq_varid(ncid, 'parent_edge_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%parent_index))
    CALL nf(nf_inq_varid(ncid, 'edge_parent_type', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%parent_child_type))
    CALL nf(nf_inq_varid(ncid, 'start_idx_e', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%start_idx))
    CALL nf(nf_inq_varid(ncid, 'end_idx_e', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%end_idx))
    CALL nf(nf_inq_varid(ncid, 'refin_e_ctrl', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%refin_ctrl))

    ! read vertex info
    CALL nf(nf_inq_varid(ncid, 'vertex_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%verts%idx))
    CALL nf(nf_inq_varid(ncid, 'longitude_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%verts%vertex(:)%lon))
    CALL nf(nf_inq_varid(ncid, 'latitude_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%verts%vertex(:)%lat))
    
    ! read cartesian coordinates
    CALL nf(nf_inq_varid(ncid, 'cartesian_x_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%verts%cartesian(:)%x(1)))
    CALL nf(nf_inq_varid(ncid, 'cartesian_y_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%verts%cartesian(:)%x(2)))
    CALL nf(nf_inq_varid(ncid, 'cartesian_z_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%verts%cartesian(:)%x(3)))

    !-------------------------------
    ! additional cartesian info
    return_status = nf_inq_varid(ncid, 'edge_middle_cartesian_x', varid)
    IF (return_status /= nf_noerr) THEN
      CALL warning(method_name, "Did not find cartesian grid info")
    ! fill some basic fields
!       CALL geographical_to_cartesian(grid_obj%cells%center, no_of_cells,&
!         & grid_obj%cells%cartesian_center)
    ELSE
      ! read all the cartesian info

      !  edges
      CALL nf(nf_inq_varid(ncid, 'edge_middle_cartesian_x', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_center(:)%x(1)))
      CALL nf(nf_inq_varid(ncid, 'edge_middle_cartesian_y', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_center(:)%x(2)))
      CALL nf(nf_inq_varid(ncid, 'edge_middle_cartesian_z', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_center(:)%x(3)))

      CALL nf(nf_inq_varid(ncid, 'edge_dual_middle_cartesian_x', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%dual_cartesian_center(:)%x(1)))
      CALL nf(nf_inq_varid(ncid, 'edge_dual_middle_cartesian_y', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%dual_cartesian_center(:)%x(2)))
      CALL nf(nf_inq_varid(ncid, 'edge_dual_middle_cartesian_z', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%dual_cartesian_center(:)%x(3)))

      CALL nf(nf_inq_varid(ncid, 'edge_primal_normal_cartesian_x', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_primal_normal(:)%x(1)))
      CALL nf(nf_inq_varid(ncid, 'edge_primal_normal_cartesian_y', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_primal_normal(:)%x(2)))
      CALL nf(nf_inq_varid(ncid, 'edge_primal_normal_cartesian_z', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_primal_normal(:)%x(3)))
      
      CALL nf(nf_inq_varid(ncid, 'edge_dual_normal_cartesian_x', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_dual_normal(:)%x(1)))
      CALL nf(nf_inq_varid(ncid, 'edge_dual_normal_cartesian_y', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_dual_normal(:)%x(2)))
      CALL nf(nf_inq_varid(ncid, 'edge_dual_normal_cartesian_z', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%edges%cartesian_dual_normal(:)%x(3)))

      ! cells
      CALL nf(nf_inq_varid(ncid, 'cell_circumcenter_cartesian_x', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%cartesian_center(:)%x(1)))
      CALL nf(nf_inq_varid(ncid, 'cell_circumcenter_cartesian_y', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%cartesian_center(:)%x(2)))
      CALL nf(nf_inq_varid(ncid, 'cell_circumcenter_cartesian_z', varid))
      CALL nf(nf_get_var_double(ncid, varid, grid_obj%cells%cartesian_center(:)%x(3)))
     
    ENDIF
    !-------------------------------

    CALL nf(nf_inq_varid(ncid, 'dual_area_p', varid))
    CALL nf(nf_get_var_double(ncid, varid, grid_obj%verts%dual_area))
    ! read vertex connectivity
    ! allocate tmp arrays for reading
    ALLOCATE(tmp_index(no_of_verts,max_vert_connect),stat=i)
    IF (i > 0) THEN
      CALL finish (method_name, 'failed to ALLOCATE(tmp_index)')
    ENDIF

    CALL nf(nf_inq_varid(ncid, 'cells_of_vertex', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_vert_connect
      grid_obj%verts%get_cell_index(:,i) = tmp_index(:,i)
    ENDDO

    CALL nf(nf_inq_varid(ncid, 'edges_of_vertex', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_vert_connect
      grid_obj%verts%get_edge_index(:,i) = tmp_index(:,i)
    ENDDO

    CALL nf(nf_inq_varid(ncid, 'edge_orientation', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_vert_connect
      grid_obj%verts%get_edge_orient(:,i) = tmp_index(:,i)
    ENDDO

    CALL nf(nf_inq_varid(ncid, 'vertices_of_vertex', varid))
    CALL nf(nf_get_var_int   (ncid, varid, tmp_index))
    DO i=1,max_vert_connect
      grid_obj%verts%get_neighbor_index(:,i) = tmp_index(:,i)
    ENDDO
    DEALLOCATE(tmp_index)

    ! read verts nest hierarchy info
    CALL nf(nf_inq_varid(ncid, 'parent_vertex_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%verts%parent_index))
    CALL nf(nf_inq_varid(ncid, 'start_idx_v', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%verts%start_idx))
    CALL nf(nf_inq_varid(ncid, 'end_idx_v', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%verts%end_idx))
    CALL nf(nf_inq_varid(ncid, 'refin_v_ctrl', varid))
    CALL nf(nf_get_var_int   (ncid, varid, grid_obj%verts%refin_ctrl))

    return_status = nf_inq_varid(ncid, 'phys_edge_id', varid)
    IF (return_status == nf_noerr) THEN
      CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%subgrid_id))
      CALL nf(nf_inq_varid(ncid, 'phys_cell_id', varid))
      CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%subgrid_id))
    ELSE
      grid_obj%edges%subgrid_id(:) = 0
      grid_obj%edges%subgrid_id(:) = 0
    ENDIF

!    print *, 'grid_obj%start_subgrid_id:', grid_obj%start_subgrid_id
    IF (grid_obj%start_subgrid_id == -1 .AND. start_subgrid_id /= 0) THEN
      ! normalize start_subgrid_id to 0
      grid_obj%edges%subgrid_id(1:no_of_edges) = &
        &  grid_obj%edges%subgrid_id(1:no_of_edges) - start_subgrid_id
      grid_obj%cells%subgrid_id(1:no_of_cells) = &
        & grid_obj%cells%subgrid_id(1:no_of_cells) - start_subgrid_id
      grid_obj%start_subgrid_id = 0

    ELSE IF (grid_obj%start_subgrid_id > 0) THEN
!      print *, 'grid_obj%start_subgrid_id > 0'
      grid_obj%edges%subgrid_id(1:no_of_edges) = grid_obj%edges%subgrid_id(1:no_of_edges) - &
        &   start_subgrid_id + grid_obj%start_subgrid_id
      grid_obj%cells%subgrid_id(1:no_of_cells) = grid_obj%cells%subgrid_id(1:no_of_cells) - &
        & start_subgrid_id + grid_obj%start_subgrid_id
!      print *,  'grid_obj%cells%subgrid_id(1):', grid_obj%cells%subgrid_id(1)

    ELSE
       grid_obj%start_subgrid_id = start_subgrid_id
    ENDIF

    !----------------------------------
    !  read decompositions
    return_status = nf_inq_varid(ncid, 'cell_domain_id', varid)
    IF (return_status == nf_noerr) THEN
      CALL nf(nf_get_var_int(ncid, varid, grid_obj%cells%domain_id(:,:)))
      CALL nf(nf_inq_varid(ncid, 'cell_no_of_domains', varid))
      CALL nf(nf_get_var_int(ncid, varid, grid_obj%cells%no_of_domains))
    ENDIF
    !----------------------------------
        
    !----------------------------------
    IF (PRESENT(read_grid_ids)) THEN
      IF (read_grid_ids) THEN
        CALL nf(nf_get_att_int(ncid, nf_global,'grid_ID', grid_obj%patch_id))
        CALL nf(nf_get_att_int(ncid, nf_global,'parent_grid_ID', &
          & grid_obj%parent_grid_id))
        CALL nf(nf_inq_varid(ncid, 'child_cell_id', varid))
        CALL nf(nf_get_var_int   (ncid, varid, grid_obj%cells%child_id))
        CALL nf(nf_inq_varid(ncid, 'child_edge_id', varid))
        CALL nf(nf_get_var_int   (ncid, varid, grid_obj%edges%child_id))
      ENDIF
    ENDIF

    CALL nf(nf_close(ncid))

    !-------------------------------------------------------------

    grid_obj%is_filled = .true.

  END SUBROUTINE read_netcdf_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_netcdf_grid(grid_id, write_file_name)
    INTEGER, INTENT(in) :: grid_id
    CHARACTER(LEN=filename_max), OPTIONAL, INTENT(in) :: write_file_name

    TYPE(t_grid), POINTER :: grid_obj
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    INTEGER :: old_mode
    INTEGER :: ncid
    INTEGER :: dimids(2), dim_domain_id(2)

    INTEGER :: dim_ncell, dim_nvertex, dim_nedge, dim_two, dim_cell_refine, &
      & dim_edge_refine, dim_vert_refine, dim_nvertex_per_cell,      &
      & dim_ncells_per_edge, dim_nedges_per_vertex, dim_nchilds_per_cell, &
      & dim_list, dim_nchdom, dim_max_decompositions

    INTEGER :: varid_clon, varid_clat, varid_clonv, varid_clatv
    INTEGER :: varid_vlon, varid_vlat, varid_vlonv, varid_vlatv
    INTEGER :: varid_elon, varid_elat, varid_elonv, varid_elatv
    INTEGER :: varid_vx, varid_vy, varid_vz

    INTEGER :: varid_carea, varid_varea

    INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8, &
      & varid9, varid10, varid11, varid12, varid13, varid14, varid15,   &
      & varid16, varid17, varid18, varid19, varid20, varid21, varid22,  &
      & varid23, varid24, varid25, varid26, varid27, varid28, varid29,  &
      & varid30, varid31, varid32, varid33, varid34, varid35, varid36,  &
      & varid37, varid38, varid39, varid40, varid41, varid42, varid43,  &
      & varid44, varid45, varid46, varid47, varid48, varid251,          &
      & varid282, varid403

    INTEGER :: varid_cell_domain_id, varid_cell_no_of_domains    
    INTEGER :: varid_cell_elevation, varid_cell_sea_land_mask
    INTEGER :: varid_edge_elevation, varid_edge_sea_land_mask
    INTEGER :: varid_phys_cell_id, varid_phys_edge_id
    INTEGER :: varid_cell_barycenter_lon, varid_cell_barycenter_lat
    INTEGER :: varid_edge_quad

    INTEGER :: varid_edge_x, varid_edge_y, varid_edge_z
    INTEGER :: varid_dualedge_x, varid_dualedge_y, varid_dualedge_z
    INTEGER :: varid_edgenormal_x, varid_edgenormal_y, varid_edgenormal_z
    INTEGER :: varid_dualedgenormal_x, varid_dualedgenormal_y, varid_dualedgenormal_z
    INTEGER :: varid_circumcenter_x, varid_circumcenter_y, varid_circumcenter_z
    
    INTEGER :: ifs2icon_cell, ifs2icon_edge, ifs2icon_vertex
        
    INTEGER :: no_of_cells, no_of_edges, no_of_verts
    INTEGER :: max_vert_connect, max_cell_vertices
    INTEGER, POINTER :: tmp_index(:,:)
    REAL(wp), POINTER :: tmp_real(:,:)
    INTEGER :: i, j, j1, j2, pole_index

    TYPE(t_uuid) :: uuid
    CHARACTER(len=uuid_string_length) :: uuid_string

    REAL(wp) :: rotation_vector(3)
    REAL(wp), ALLOCATABLE :: zv2dx(:,:), zv2dy(:,:)
    REAL(wp) :: swap(4)

    INTEGER :: ilevel, grid_root
    INTEGER :: itype_optimize=0
    LOGICAL :: l_c_grid = .false.

    INTEGER :: return_status

    !-------------------------------------------------------------------------
    ! get unique grid file identifier for GRIB2 and updated CF-Convention
    CALL uuid_generate(uuid)
    CALL uuid_unparse(uuid, uuid_string)

    grid_obj => get_grid(grid_id)
    cells => grid_obj%cells
    edges => grid_obj%edges
    verts => grid_obj%verts

!     write(0,*) "grid_obj%file_name:", TRIM(grid_obj%file_name)
!     write(0,*) "grid_obj%out_file_name:", TRIM(grid_obj%out_file_name)
    IF (PRESENT(write_file_name)) THEN
!       write(0,*) "write_file_name:", TRIM(write_file_name)
       grid_obj%file_name = write_file_name
    ELSE
       grid_obj%file_name = grid_obj%out_file_name
    ENDIF

    !-------------------------------------------------------------------------
    ! distinguish between gridgeneration for optimization strategies

    ! Dummy settings for special grids
    grid_root = 2
    ! ilevel    = 1level
    ilevel    = grid_obj%level

    !----------------------------------------------------------------------
    no_of_cells = cells%no_of_existcells
    no_of_edges = edges%no_of_existedges
    no_of_verts = verts%no_of_existvertices
    max_vert_connect  = verts%max_connectivity
    max_cell_vertices = cells%max_no_of_vertices

    rotation_vector = (/ 0.0_wp, 0.0_wp, 0.0_wp /)

    WRITE(message_text,'(a,a,a,a,a,i9)')                              &
         &             'Write grid file: ', TRIM(grid_obj%file_name), &
         &             ' uuid ', TRIM(uuid_string),                   &
         &             ' number of cells ', no_of_cells
    CALL message ('', TRIM(message_text))
    !----------------------------------------------------------------------

    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))
    CALL nf(nf_create(TRIM(grid_obj%file_name), nf_clobber, ncid))
    CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))

    !----------------------------------------------------------------------
    ! Global attributes
    CALL nf(nf_put_att_text    (ncid, nf_global, 'title', 21, 'ICON grid description'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'institution', 59, &
      & 'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'source', 8, 'icon-dev'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'uuid' , uuid_string_length, TRIM(uuid_string)))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'grid_mapping_name' , 18, 'lat_long_on_sphere'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_id' , 28, 'urn:ogc:def:cs:EPSG:6.0:6422'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_name',30,'Spherical 2D Coordinate System'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'ellipsoid_name' , 6, 'Sphere'))
!    CALL nf(nf_put_att_double  (ncid, nf_global, 'semi_major_axis' , nf_double, 1, re))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'semi_major_axis' , nf_double, 1, &
      & grid_obj%geometry_info%sphere_radius))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'inverse_flattening' , nf_double, 1, 0.0_wp))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_level', nf_int, 1, ilevel))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_root', nf_int, 1, grid_root))
    
    ! The following three attributes are nontrivial in the presence of grid refinement
    ! and are set here because they are checked in the ICON code
    IF (grid_obj%patch_id < 0) THEN
      CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_ID', nf_int, 1,1))
    ELSE
      CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_ID', nf_int, 1,grid_obj%patch_id))
    ENDIF

    CALL nf(nf_put_att_int     (ncid, nf_global, 'parent_grid_ID', nf_int, 1, &
      & grid_obj%parent_grid_id))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'no_of_subgrids', nf_int, 1, &
      & grid_obj%no_of_subgrids))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'start_subgrid_id', nf_int, 1, &
      & grid_obj%start_subgrid_id))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'max_childdom', nf_int, 1, 1))

!     SELECT CASE (itype_optimize)
!     CASE (0)
!       CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 12, 'natural grid'))
!     CASE (1)
!       IF (l_c_grid) THEN
!         CALL nf(nf_put_att_text(ncid, nf_global, 'grid_optimization', 50, &
!           & 'Heikes-Randall with additional c-grid optimization'))
!       ELSE
!         CALL nf(nf_put_att_text(ncid, nf_global, 'grid_optimization', 27, &
!           & 'Heikes-Randall optimization'))
!       ENDIF
!     CASE (2)
!       CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 22, &
!         & 'equal area subdivision'))
!     CASE (3)
!       CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 30, &
!         & 'c-grid small circle constraint'))
!     END SELECT
    CALL nf(nf_put_att_double  (ncid, nf_global, 'rotation_vector', nf_double, 3, &
      & rotation_vector))


    return_status = write_geometry_info(ncid, grid_obj%geometry_info)
      
    !----------------------------------------------------------------------
    ! Dimensions
    CALL nf(nf_def_dim(ncid, 'cell',   no_of_cells, dim_ncell))
    CALL nf(nf_def_dim(ncid, 'vertex', no_of_verts, dim_nvertex))
    CALL nf(nf_def_dim(ncid, 'edge',   no_of_edges, dim_nedge))
    !
    CALL nf(nf_def_dim(ncid, 'nc',        2, dim_ncells_per_edge))
    CALL nf(nf_def_dim(ncid, 'nv',  cells%max_no_of_vertices, dim_nvertex_per_cell))
    CALL nf(nf_def_dim(ncid, 'ne',  verts%max_connectivity, dim_nedges_per_vertex))
    !
    CALL nf(nf_def_dim(ncid, 'no',        4, dim_nchilds_per_cell))

    ! Dimensions for refinement
    !      write(*,*) "writing Dimensions for refinement..."
    !       call flush(6)
    CALL nf(nf_def_dim(ncid, 'two_grf',     2,    dim_two))
    CALL nf(nf_def_dim(ncid, 'max_chdom',   1, dim_nchdom))
    dim_list = max_rlcell-min_rlcell+1
    CALL nf(nf_def_dim(ncid, 'cell_grf',dim_list, dim_cell_refine))
    dim_list = max_rledge-min_rledge+1
    CALL nf(nf_def_dim(ncid, 'edge_grf',dim_list, dim_edge_refine))
    dim_list = max_rlvert-min_rlvert+1
    CALL nf(nf_def_dim(ncid, 'vert_grf',dim_list, dim_vert_refine))

    CALL nf(nf_def_dim(ncid, 'max_stored_decompositions',   max_decompositions, &
      & dim_max_decompositions))
    !---------------------------------------------------------------------
    !
    ! Grid variables
    !---------------------------------------------------------------------
    !
    ! public grid information
    !
    ! cell part:
    ! write(*,*) "writing public grid information..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'clon', nf_double, 1, dim_ncell, varid_clon))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'long_name', 16, 'center longitude'))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'standard_name', 14, 'grid_longitude'))
    CALL nf(nf_put_att_text(ncid, varid_clon, 'bounds', 13, 'clon_vertices'))
    !
    CALL nf(nf_def_var(ncid, 'clat', nf_double, 1, dim_ncell, varid_clat))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'long_name', 15, 'center latitude'))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'standard_name', 13, 'grid_latitude'))
    CALL nf(nf_put_att_text(ncid, varid_clat, 'bounds', 13, 'clat_vertices'))
    !
    dimids = (/ dim_nvertex_per_cell, dim_ncell /)
    CALL nf(nf_def_var(ncid, 'clon_vertices', nf_double, 2, dimids, varid_clonv))
    CALL nf(nf_put_att_text(ncid, varid_clonv, 'units', 6, 'radian'))
    !
    dimids = (/ dim_nvertex_per_cell, dim_ncell /)
    CALL nf(nf_def_var(ncid, 'clat_vertices', nf_double, 2, dimids, varid_clatv))
    CALL nf(nf_put_att_text(ncid, varid_clatv, 'units', 6, 'radian'))
    !
    ! vertex part:
    !
    ! write(*,*) "writing public grid information vertex part..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'vlon', nf_double, 1, dim_nvertex, varid_vlon))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'long_name', 16, 'vertex longitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'standard_name', 14, 'grid_longitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlon, 'bounds', 13, 'vlon_vertices'))
    !
    CALL nf(nf_def_var(ncid, 'vlat', nf_double, 1, dim_nvertex, varid_vlat))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'long_name', 15, 'vertex latitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'standard_name', 13, 'grid_latitude'))
    CALL nf(nf_put_att_text(ncid, varid_vlat, 'bounds', 13, 'vlat_vertices'))
    !
    dimids = (/ dim_nedges_per_vertex, dim_nvertex /)
    CALL nf(nf_def_var(ncid, 'vlon_vertices', nf_double, 2, dimids, varid_vlonv))
    CALL nf(nf_put_att_text(ncid, varid_vlonv, 'units', 6, 'radian'))
    !
    dimids = (/ dim_nedges_per_vertex, dim_nvertex /)
    CALL nf(nf_def_var(ncid, 'vlat_vertices', nf_double, 2, dimids, varid_vlatv))
    CALL nf(nf_put_att_text(ncid, varid_vlatv, 'units', 6, 'radian'))
    !
    ! edge part:
    !
    ! write(*,*) "writing public grid information edge part..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'elon', nf_double, 1, dim_nedge, varid_elon))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'long_name', 23, 'edge midpoint longitude'))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'standard_name', 14, 'grid_longitude'))
    CALL nf(nf_put_att_text(ncid, varid_elon, 'bounds', 13, 'elon_vertices'))
    !
    CALL nf(nf_def_var(ncid, 'elat', nf_double, 1, dim_nedge, varid_elat))
    CALL nf(nf_put_att_text(ncid, varid_elat, 'long_name', 22, 'edge midpoint latitude') )
    CALL nf(nf_put_att_text(ncid, varid_elat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_elat, 'standard_name', 13, 'grid_latitude'))
    CALL nf(nf_put_att_text(ncid, varid_elat, 'bounds', 13, 'elat_vertices'))
    !
    dimids = (/ dim_nchilds_per_cell, dim_nedge /)
    CALL nf(nf_def_var(ncid, 'elon_vertices', nf_double, 2, dimids, varid_elonv))
    CALL nf(nf_put_att_text(ncid, varid_elonv, 'units', 6, 'radian'))
    !
    dimids = (/ dim_nchilds_per_cell, dim_nedge /)
    CALL nf(nf_def_var(ncid, 'elat_vertices', nf_double, 2, dimids, varid_elatv))
    CALL nf(nf_put_att_text(ncid, varid_elatv, 'units', 6, 'radian'))
    !
    !---------------------------------------------------------------------
    !
    ! ifs2icon variables 
    !
    CALL nf(nf_def_var(ncid, 'ifs2icon_cell_grid', nf_double, 1, dim_ncell, ifs2icon_cell))
    CALL nf(nf_put_att_text(ncid, ifs2icon_cell, 'long_name', 17, 'ifs to icon cells'))
    CALL nf(nf_put_att_text(ncid, ifs2icon_cell, 'coordinates', 9, 'clon clat'))
    !
    CALL nf(nf_def_var(ncid, 'ifs2icon_edge_grid', nf_double, 1, dim_nedge, ifs2icon_edge))
    CALL nf(nf_put_att_text(ncid, ifs2icon_edge, 'long_name', 16, 'ifs to icon edge'))
    CALL nf(nf_put_att_text(ncid, ifs2icon_edge, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'ifs2icon_vertex_grid', nf_double, 1, dim_nvertex, ifs2icon_vertex))
    CALL nf(nf_put_att_text(ncid, ifs2icon_vertex, 'long_name', 18, 'ifs to icon vertex'))
    CALL nf(nf_put_att_text(ncid, ifs2icon_vertex, 'coordinates', 9, 'vlon vlat'))
    !

    !---------------------------------------------------------------------
    !
    ! test variables (areas)
    !
    ! write(*,*) "writing public grid information areas part..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'cell_area', nf_double, 1, dim_ncell, varid_carea))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'long_name', 17, 'area of grid cell'))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'units', 2, 'm2'))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'standard_name', 4, 'area'))
    CALL nf(nf_put_att_text(ncid, varid_carea, 'coordinates', 9, 'clon clat'))
    !
    CALL nf(nf_def_var(ncid, 'dual_area', nf_double, 1, dim_nvertex, varid_varea))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'long_name', 40, &
      & 'areas of dual hexagonal/pentagonal cells'))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'units', 2, 'm2'))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'standard_name', 4, 'area'))
    CALL nf(nf_put_att_text(ncid, varid_varea, 'coordinates', 9, 'vlon vlat'))
    !

    !---------------------------------------------------------------------
    !
    ! private grid information
    !
    ! write(*,*) "writing private grid in(1:no_of_cells)formation ..."
    !       call flush(6)
   !
    CALL nf(nf_def_var(ncid, 'phys_cell_id', nf_int, 1, dim_ncell, varid_phys_cell_id))
    CALL nf(nf_put_att_text(ncid, varid_phys_cell_id, 'long_name', 26, &
      & 'physical domain ID of cell'))
    CALL nf(nf_put_att_text(ncid, varid_phys_cell_id, 'coordinates', 9, 'clon clat'))
!     CALL nf(nf_put_att_text(ncid, varid_phys_cell_id, 'cdi', 6, 'ignore'))

    CALL nf(nf_def_var(ncid, 'phys_edge_id', nf_int, 1, dim_nedge, varid_phys_edge_id))
    CALL nf(nf_put_att_text(ncid, varid_phys_edge_id, 'long_name', 26, &
      & 'physical domain ID of edge'))
    CALL nf(nf_put_att_text(ncid, varid_phys_edge_id, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid_phys_edge_id, 'cdi', 6, 'ignore'))

    CALL nf(nf_def_var(ncid, 'lon_cell_centre', nf_double, 1, dim_ncell, varid1))
    CALL nf(nf_put_att_text(ncid, varid1, 'long_name', 24, 'longitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid1, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid1, 'coordinates', 9, 'clon clat'))
!     CALL nf(nf_put_att_text(ncid, varid1, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_cell_centre', nf_double, 1, dim_ncell, varid2))
    CALL nf(nf_put_att_text(ncid, varid2, 'long_name', 23, 'latitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid2, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid2, 'coordinates', 9, 'clon clat'))
!     CALL nf(nf_put_att_text(ncid, varid2, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_cell_barycenter', nf_double, 1, dim_ncell, &
      varid_cell_barycenter_lat))
    CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lat, 'long_name', 27, &
       'latitude of cell barycenter'))
    CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lat, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lat, 'coordinates', 9, 'clon clat'))
!     CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lat, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lon_cell_barycenter', nf_double, 1, dim_ncell, &
      varid_cell_barycenter_lon))
    CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lon, 'long_name', 28, &
       'longitude of cell barycenter'))
    CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lon, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lon, 'coordinates', 9, 'clon clat'))
!     CALL nf(nf_put_att_text(ncid, varid_cell_barycenter_lon, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'longitude_vertices', nf_double, 1, dim_nvertex, varid3))
    CALL nf(nf_put_att_text(ncid, varid3, 'long_name', 21, 'longitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid3, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid3, 'coordinates', 9, 'vlon vlat'))
!     CALL nf(nf_put_att_text(ncid, varid3, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'latitude_vertices', nf_double, 1, dim_nvertex, varid4))
    CALL nf(nf_put_att_text(ncid, varid4, 'long_name', 20, 'latitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid4, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid4, 'coordinates', 9, 'vlon vlat'))
!     CALL nf(nf_put_att_text(ncid, varid4, 'cdi', 6, 'ignore'))
    !

    
    CALL nf(nf_def_var(ncid, 'lon_edge_centre', nf_double, 1, dim_nedge, varid5))
    CALL nf(nf_put_att_text(ncid, varid5, 'long_name', 28, 'longitudes of edge midpoints'))
    CALL nf(nf_put_att_text(ncid, varid5, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid5, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid5, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_edge_centre', nf_double, 1, dim_nedge, varid6))
    CALL nf(nf_put_att_text(ncid, varid6, 'long_name', 27, 'latitudes of edge midpoints'))
    CALL nf(nf_put_att_text(ncid, varid6, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid6, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid6, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'edge_of_cell', nf_int, 2, dimids, varid7))
    CALL nf(nf_put_att_text(ncid, varid7, 'long_name', 29, 'edges of each cell'))
!     CALL nf(nf_put_att_text(ncid, varid7, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'vertex_of_cell', nf_int, 2, dimids, varid8))
    CALL nf(nf_put_att_text(ncid, varid8, 'long_name', 32, 'vertices of each cell'))
!     CALL nf(nf_put_att_text(ncid, varid8, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'adjacent_cell_of_edge', nf_int, 2, dimids, varid9))
    CALL nf(nf_put_att_text(ncid, varid9, 'long_name', 27, 'cells adjacent to each edge'))
!     CALL nf(nf_put_att_text(ncid, varid9, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vertices', nf_int, 2, dimids, varid10))
    CALL nf(nf_put_att_text(ncid, varid10, 'long_name', 35, &
      & 'vertices at the end of of each edge'))
!     CALL nf(nf_put_att_text(ncid, varid10, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'cells_of_vertex', nf_int, 2, dimids, varid11))
    CALL nf(nf_put_att_text(ncid, varid11, 'long_name', 24, 'cells around each vertex'))
!     CALL nf(nf_put_att_text(ncid, varid11, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'edges_of_vertex', nf_int, 2, dimids, varid12))
    CALL nf(nf_put_att_text(ncid, varid12, 'long_name', 24, 'edges around each vertex'))
!     CALL nf(nf_put_att_text(ncid, varid12, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'vertices_of_vertex', nf_int, 2, dimids, varid13))
    CALL nf(nf_put_att_text(ncid, varid13, 'long_name', 27, 'vertices around each vertex'))
!     CALL nf(nf_put_att_text(ncid, varid13, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_area_p', nf_double, 1, dim_ncell, varid14))
    CALL nf(nf_put_att_text(ncid, varid14, 'long_name', 17, 'area of grid cell'))
    CALL nf(nf_put_att_text(ncid, varid14, 'units', 2, 'm2'))
    CALL nf(nf_put_att_text(ncid, varid14, 'coordinates', 9, 'clon clat'))
!     CALL nf(nf_put_att_text(ncid, varid14, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_elevation', nf_double, 1, dim_ncell, varid_cell_elevation))
    CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'long_name', 29, &
      & 'elevation at the cell centers'))
    CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'units', 1, 'm'))
!     CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'cdi', 6, 'ignore'))
    CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'coordinates', 9, 'clon clat'))
    
    CALL nf(nf_def_var(ncid, 'cell_sea_land_mask', nf_int, 1, dim_ncell, &
      & varid_cell_sea_land_mask))
    CALL nf(nf_put_att_text(ncid, varid_cell_sea_land_mask, 'long_name', 72, &
      & 'sea (-2 inner, -1 boundary) land (2 inner, 1 boundary) mask for the cells'))
    CALL nf(nf_put_att_text(ncid, varid_cell_sea_land_mask, 'units', 8, '2,1,-1,-2'))
    CALL nf(nf_put_att_text(ncid, varid_cell_sea_land_mask, 'coordinates', 9, 'clon clat'))

    dim_domain_id = (/ dim_max_decompositions, dim_ncell /)
    CALL nf(nf_def_var(ncid, 'cell_domain_id', nf_int, 2, dim_domain_id, &
      & varid_cell_domain_id))
    CALL nf(nf_put_att_text(ncid, varid_cell_domain_id, 'long_name', 32, &
      & 'cell domain id for decomposition'))
    CALL nf(nf_put_att_text(ncid, varid_cell_domain_id, 'coordinates', 9, 'clon clat'))
    
    CALL nf(nf_def_var(ncid, 'cell_no_of_domains', nf_int, 1, dim_max_decompositions, &
      & varid_cell_no_of_domains))
    CALL nf(nf_put_att_text(ncid, varid_cell_no_of_domains, 'long_name', 40, &
      & 'number of domains for each decomposition'))
    
    CALL nf(nf_def_var(ncid, 'dual_area_p', nf_double, 1, dim_nvertex, varid15))
    CALL nf(nf_put_att_text(ncid, varid15, 'long_name', 40, &
      & 'areas of dual hexagonal/pentagonal cells'))
    CALL nf(nf_put_att_text(ncid, varid15, 'units', 2, 'm2'))
!     CALL nf(nf_put_att_text(ncid, varid15, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_length', nf_double, 1, dim_nedge, varid16))
    CALL nf(nf_put_att_text(ncid, varid16, 'long_name', 36, &
      & 'lengths of edges of triangular cells'))
    CALL nf(nf_put_att_text(ncid, varid16, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid16, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid16, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_cell_distance', nf_double, 2, dimids, varid40))
    CALL nf(nf_put_att_text(ncid, varid40, 'long_name', 63, &
      & 'distances between edge midpoint and adjacent triangle midpoints'))
    CALL nf(nf_put_att_text(ncid, varid40, 'units', 1, 'm'))
!     CALL nf(nf_put_att_text(ncid, varid40, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'dual_edge_length', nf_double, 1, dim_nedge, varid17))
    CALL nf(nf_put_att_text(ncid, varid17, 'long_name', 71, &
      & 'lengths of dual edges (distances between triangular cell circumcenters)'))
    CALL nf(nf_put_att_text(ncid, varid17, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid17, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid17, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edgequad_area', nf_double, 1, dim_nedge, varid_edge_quad))
    CALL nf(nf_put_att_text(ncid, varid_edge_quad, 'long_name', 57, &
      & 'area around the edge formed by the two adjacent triangles'))
    CALL nf(nf_put_att_text(ncid, varid_edge_quad, 'units', 2, 'm2'))
    CALL nf(nf_put_att_text(ncid, varid_edge_quad, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_elevation', nf_double, 1, dim_nedge, varid_edge_elevation))
    CALL nf(nf_put_att_text(ncid, varid_edge_elevation, 'long_name', 29, &
      & 'elevation at the edge centers'))
    CALL nf(nf_put_att_text(ncid, varid_edge_elevation, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid_edge_elevation, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid_edge_elevation, 'cdi', 6, 'ignore'))

    CALL nf(nf_def_var(ncid, 'edge_sea_land_mask', nf_int, 1, dim_nedge, &
      & varid_edge_sea_land_mask))
    CALL nf(nf_put_att_text(ncid, varid_edge_sea_land_mask, 'long_name', 72, &
      & 'sea (-2 inner, -1 boundary) land (2 inner, 1 boundary) mask for the cells'))
    CALL nf(nf_put_att_text(ncid, varid_edge_sea_land_mask, 'units', 8, '2,1,-1,-2'))
    CALL nf(nf_put_att_text(ncid, varid_edge_sea_land_mask, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid_edge_sea_land_mask, 'cdi', 6, 'ignore'))

    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vert_distance', nf_double, 2, dimids, varid18))
    CALL nf(nf_put_att_text(ncid, varid18, 'long_name', 57, &
      & 'distances between edge midpoint and vertices of that edge'))
    CALL nf(nf_put_att_text(ncid, varid18, 'units', 1, 'm'))
!     CALL nf(nf_put_att_text(ncid, varid18, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'zonal_normal_primal_edge', nf_double, 1, dim_nedge, varid19))
    CALL nf(nf_put_att_text(ncid, varid19, 'long_name', 40, &
      & 'zonal component of normal to primal edge'))
    CALL nf(nf_put_att_text(ncid, varid19, 'units', 6, 'radian'))
!     CALL nf(nf_put_att_text(ncid, varid19, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'meridional_normal_primal_edge', nf_double, 1, dim_nedge, varid20))
    CALL nf(nf_put_att_text(ncid, varid20, 'long_name', 45, &
      & 'meridional component of normal to primal edge'))
    CALL nf(nf_put_att_text(ncid, varid20, 'units', 6, 'radian'))
!     CALL nf(nf_put_att_text(ncid, varid20, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'zonal_normal_dual_edge', nf_double, 1, dim_nedge, varid21))
    CALL nf(nf_put_att_text(ncid, varid21, 'long_name', 38, &
      & 'zonal component of normal to dual edge'))
    CALL nf(nf_put_att_text(ncid, varid21, 'units', 6, 'radian'))
!     CALL nf(nf_put_att_text(ncid, varid21, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'meridional_normal_dual_edge', &
      & nf_double, 1, dim_nedge, varid22))
    CALL nf(nf_put_att_text(ncid, varid22, 'long_name', 43, &
      & 'meridional component of normal to dual edge'))
    CALL nf(nf_put_att_text(ncid, varid22, 'units', 6, 'radian'))
!     CALL nf(nf_put_att_text(ncid, varid22, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'orientation_of_normal', nf_int, 2, dimids, varid23))
    CALL nf(nf_put_att_text(ncid, varid23, 'long_name', 48, &
      & 'orientations of normals to triangular cell edges'))
!     CALL nf(nf_put_att_text(ncid, varid23, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_index', nf_int, 1, dim_ncell, varid24))
    CALL nf(nf_put_att_text(ncid, varid24, 'long_name', 10, 'cell index'))
!     CALL nf(nf_put_att_text(ncid, varid24, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_cell_index', nf_int, 1, dim_ncell, varid25))
    CALL nf(nf_put_att_text(ncid, varid25, 'long_name', 17, 'parent cell index'))
!     CALL nf(nf_put_att_text(ncid, varid25, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_cell_type', nf_int, 1, dim_ncell, varid251))
    CALL nf(nf_put_att_text(ncid, varid251, 'long_name', 16, 'parent cell type'))
!     CALL nf(nf_put_att_text(ncid, varid251, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'neighbor_cell_index', nf_int, 2, dimids, varid26))
    CALL nf(nf_put_att_text(ncid, varid26, 'long_name', 19, 'cell neighbor index'))
!     CALL nf(nf_put_att_text(ncid, varid26, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_cell_index', nf_int, 2, dimids, varid27))
    CALL nf(nf_put_att_text(ncid, varid27, 'long_name', 16, 'child cell index'))
!     CALL nf(nf_put_att_text(ncid, varid27, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'child_cell_id', nf_int, 1, dim_ncell, varid41))
    CALL nf(nf_put_att_text(ncid, varid41, 'long_name', 23, 'domain ID of child cell'))
!     CALL nf(nf_put_att_text(ncid, varid41, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_index', nf_int, 1, dim_nedge, varid28))
    CALL nf(nf_put_att_text(ncid, varid28, 'long_name', 10, 'edge index'))
!     CALL nf(nf_put_att_text(ncid, varid28, 'cdi', 6, 'ignore'))
    !
    !       CALL nf(nf_def_var(ncid, 'edge_parent', NF_INT, 1, dim_nedge, varid281))
    !       CALL nf(nf_put_att_text(ncid, varid281, 'long_name', 10, 'edge parent'))
    !       CALL nf(nf_put_att_text(ncid, varid281, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_parent_type', nf_int, 1, dim_nedge, varid282))
    CALL nf(nf_put_att_text(ncid, varid282, 'long_name', 10, 'edge parent type'))
!     CALL nf(nf_put_att_text(ncid, varid282, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'vertex_index', nf_int, 1, dim_nvertex, varid29))
    CALL nf(nf_put_att_text(ncid, varid29, 'long_name', 14, 'vertices index'))
!     CALL nf(nf_put_att_text(ncid, varid29, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'edge_orientation', nf_int, 2, dimids, varid30))
    CALL nf(nf_put_att_text(ncid, varid30, 'long_name', 16, 'edge orientation'))
!     CALL nf(nf_put_att_text(ncid, varid30, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_system_orientation', nf_int, 1, dim_nedge, varid31))
    CALL nf(nf_put_att_text(ncid, varid31, 'long_name', 23, 'edge system orientation'))
    CALL nf(nf_put_att_text(ncid, varid31, 'coordinates', 9, 'elon elat'))
!     CALL nf(nf_put_att_text(ncid, varid31, 'cdi', 6, 'ignore'))
    !
    ! Variables added for mesh refinement
    !
    ! write(*,*) "writing variables added for mesh refinement ..."
    !       write(*,*) "dim_list=",dim_list
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'refin_c_ctrl', nf_int, 1, dim_ncell, varid32))
    CALL nf(nf_put_att_text(ncid, varid32, 'long_name', 33, 'refinement control flag for cells'))
!     CALL nf(nf_put_att_text(ncid, varid32, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def index_c_list ..."
    !       call flush(6)
    dimids = (/ dim_cell_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_c_list', nf_int, 2, dimids, varid33))
    CALL nf(nf_put_att_text(ncid, varid33, 'long_name', 73, &
      & 'list of start and end indices for each refinement control level for cells'))
!     CALL nf(nf_put_att_text(ncid, varid33, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def start_idx_c..."
    !       call flush(6)
    dimids = (/ dim_cell_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_c', nf_int, 2, dimids, varid43))
    CALL nf(nf_put_att_text(ncid, varid43, 'long_name', 65, &
      & 'list of start indices for each refinement control level for cells'))
!     CALL nf(nf_put_att_text(ncid, varid43, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def end_idx_c..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'end_idx_c', nf_int, 2, dimids, varid44))
    CALL nf(nf_put_att_text(ncid, varid44, 'long_name', 63, &
      & 'list of end indices for each refinement control level for cells'))
!     CALL nf(nf_put_att_text(ncid, varid44, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def refin_e_ctrl..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'refin_e_ctrl', nf_int, 1, dim_nedge, varid34))
    CALL nf(nf_put_att_text(ncid, varid34, 'long_name', 33, 'refinement control flag for edges'))
!     CALL nf(nf_put_att_text(ncid, varid34, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def index_e_list..."
    !       call flush(6)
    dimids = (/ dim_edge_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_e_list', nf_int, 2, dimids, varid35))
    CALL nf(nf_put_att_text(ncid, varid35, 'long_name', 73, &
      & 'list of start and end indices for each refinement control level for edges'))
!     CALL nf(nf_put_att_text(ncid, varid35, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def start_idx_e..."
    !       call flush(6)
    dimids = (/ dim_edge_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_e', nf_int, 2, dimids, varid45))
    CALL nf(nf_put_att_text(ncid, varid45, 'long_name', 65, &
      & 'list of start indices for each refinement control level for edges'))
!     CALL nf(nf_put_att_text(ncid, varid45, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def end_idx_e..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'end_idx_e', nf_int, 2, dimids, varid46))
    CALL nf(nf_put_att_text(ncid, varid46, 'long_name', 63, &
      & 'list of end indices for each refinement control level for edges'))
!     CALL nf(nf_put_att_text(ncid, varid46, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def refin_v_ctrl..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'refin_v_ctrl', nf_int, 1, dim_nvertex, varid36))
    CALL nf(nf_put_att_text(ncid, varid36, 'long_name', 36, &
      & 'refinement control flag for vertices'))
!     CALL nf(nf_put_att_text(ncid, varid36, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def index_v_list..."
    !       call flush(6)
    dimids = (/ dim_vert_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_v_list', nf_int, 2, dimids, varid37))
    CALL nf(nf_put_att_text(ncid, varid37, 'long_name', 76, &
      & 'list of start and end indices for each refinement control level for vertices'))
!     CALL nf(nf_put_att_text(ncid, varid37, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def start_idx_v..."
    !       call flush(6)
    dimids = (/ dim_vert_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_v', nf_int, 2, dimids, varid47))
    CALL nf(nf_put_att_text(ncid, varid47, 'long_name', 68, &
      & 'list of start indices for each refinement control level for vertices'))
!     CALL nf(nf_put_att_text(ncid, varid47, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def end_idx_v..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'end_idx_v', nf_int, 2, dimids, varid48))
    CALL nf(nf_put_att_text(ncid, varid48, 'long_name', 66, &
      & 'list of end indices for each refinement control level for vertices'))
!     CALL nf(nf_put_att_text(ncid, varid48, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def parent_edge_index..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'parent_edge_index', nf_int, 1, dim_nedge, varid38))
    CALL nf(nf_put_att_text(ncid, varid38, 'long_name', 17, 'parent edge index'))
!     CALL nf(nf_put_att_text(ncid, varid38, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def child_edge_index..."
    !       call flush(6)
    dimids = (/ dim_nedge, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_edge_index', nf_int, 2, dimids, varid39))
    CALL nf(nf_put_att_text(ncid, varid39, 'long_name', 16, 'child edge index'))
!     CALL nf(nf_put_att_text(ncid, varid39, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def child_edge_id..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'child_edge_id', nf_int, 1, dim_nedge, varid42))
    CALL nf(nf_put_att_text(ncid, varid42, 'long_name', 23, 'domain ID of child edge'))
!     CALL nf(nf_put_att_text(ncid, varid42, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_vertex_index', nf_int, 1, dim_nvertex, varid403))
    CALL nf(nf_put_att_text(ncid, varid403, 'long_name', 19, 'parent vertex index'))
!     CALL nf(nf_put_att_text(ncid, varid403, 'cdi', 6, 'ignore'))
    !
    

    !-------------------------------------------
    ! define cartesian positions
    
    ! vertexes
    CALL nf(nf_def_var(ncid, 'cartesian_x_vertices', nf_double, 1, dim_nvertex, varid_vx))
    CALL nf(nf_put_att_text(ncid, varid_vx, 'long_name', 40, &
      & 'vertex cartesian coordinate x on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_vx, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_vx, 'coordinates', 9, 'vlon vlat'))
    !
    CALL nf(nf_def_var(ncid, 'cartesian_y_vertices', nf_double, 1, dim_nvertex, varid_vy))
    CALL nf(nf_put_att_text(ncid, varid_vy, 'long_name', 40, &
      & 'vertex cartesian coordinate y on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_vy, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_vy, 'coordinates', 9, 'vlon vlat'))
    !
    CALL nf(nf_def_var(ncid, 'cartesian_z_vertices', nf_double, 1, dim_nvertex, varid_vz))
    CALL nf(nf_put_att_text(ncid, varid_vz, 'long_name', 40, &
      & 'vertex cartesian coordinate x on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_vz, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_vz, 'coordinates', 9, 'vlon vlat'))
    !

    !  edges
    CALL nf(nf_def_var(ncid, 'edge_middle_cartesian_x', nf_double, 1, dim_nedge, varid_edge_x))
    CALL nf(nf_put_att_text(ncid, varid_edge_x, 'long_name', 55, &
      & 'prime edge center cartesian coordinate x on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_edge_x, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_edge_x, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_middle_cartesian_y', nf_double, 1, dim_nedge, varid_edge_y))
    CALL nf(nf_put_att_text(ncid, varid_edge_y, 'long_name', 55, &
      & 'prime edge center cartesian coordinate y on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_edge_y, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_edge_y, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_middle_cartesian_z', nf_double, 1, dim_nedge, varid_edge_z))
    CALL nf(nf_put_att_text(ncid, varid_edge_z, 'long_name', 55, &
      & 'prime edge center cartesian coordinate x on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_edge_z, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_edge_z, 'coordinates', 9, 'elon elat'))
    !
    
    CALL nf(nf_def_var(ncid, 'edge_dual_middle_cartesian_x', nf_double, 1, dim_nedge, &
      & varid_dualedge_x))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_x, 'long_name', 54, &
      & 'dual edge center cartesian coordinate x on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_x, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_x, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_dual_middle_cartesian_y', nf_double, 1, dim_nedge, &
      & varid_dualedge_y))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_y, 'long_name', 54, &
      & 'dual edge center cartesian coordinate y on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_y, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_y, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_dual_middle_cartesian_z', nf_double, 1, dim_nedge, &
      & varid_dualedge_z))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_z, 'long_name', 54, &
      & 'dual edge center cartesian coordinate z on unit sphere'))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_z, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_dualedge_z, 'coordinates', 9, 'elon elat'))
    !
    
    CALL nf(nf_def_var(ncid, 'edge_primal_normal_cartesian_x', nf_double, 1, dim_nedge, &
      & varid_edgenormal_x))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_x, 'long_name', 53, &
      & 'unit normal to the prime edge 3D vector, coordinate x'))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_x, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_x, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_primal_normal_cartesian_y', nf_double, 1, dim_nedge, &
      & varid_edgenormal_y))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_y, 'long_name', 53, &
      & 'unit normal to the prime edge 3D vector, coordinate y'))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_y, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_y, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_primal_normal_cartesian_z', nf_double, 1, dim_nedge, &
      & varid_edgenormal_z))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_z, 'long_name', 53, &
      & 'unit normal to the prime edge 3D vector, coordinate z'))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_z, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_edgenormal_z, 'coordinates', 9, 'elon elat'))
    !
    
    CALL nf(nf_def_var(ncid, 'edge_dual_normal_cartesian_x', nf_double, 1, dim_nedge, &
      & varid_dualedgenormal_x))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_x, 'long_name', 53, &
      & 'unit normal to the dual edge 3D vector, coordinate x'))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_x, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_x, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_dual_normal_cartesian_y', nf_double, 1, dim_nedge, &
      & varid_dualedgenormal_y))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_y, 'long_name', 53, &
      & 'unit normal to the dual edge 3D vector, coordinate y'))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_y, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_y, 'coordinates', 9, 'elon elat'))
    !
    CALL nf(nf_def_var(ncid, 'edge_dual_normal_cartesian_z', nf_double, 1, dim_nedge, &
      & varid_dualedgenormal_z))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_z, 'long_name', 53, &
      & 'unit normal to the dual edge 3D vector, coordinate z'))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_z, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_dualedgenormal_z, 'coordinates', 9, 'elon elat'))
    !
    
    ! cells
    CALL nf(nf_def_var(ncid, 'cell_circumcenter_cartesian_x', nf_double, 1, dim_ncell, &
      & varid_circumcenter_x))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_x, 'long_name', 82, &
      & 'cartesian position of the prime cell circumcenter on the unit sphere, coordinate x'))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_x, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_x, 'coordinates', 9, 'clon clat'))
    !
    CALL nf(nf_def_var(ncid, 'cell_circumcenter_cartesian_y', nf_double, 1, dim_ncell, &
      & varid_circumcenter_y))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_y, 'long_name', 82, &
      & 'cartesian position of the prime cell circumcenter on the unit sphere, coordinate y'))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_y, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_y, 'coordinates', 9, 'clon clat'))
    !
    CALL nf(nf_def_var(ncid, 'cell_circumcenter_cartesian_z', nf_double, 1, dim_ncell, &
      & varid_circumcenter_z))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_z, 'long_name', 82, &
      & 'cartesian position of the prime cell circumcenter on the unit sphere, coordinate z'))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_z, 'units', 6, 'meters'))
    CALL nf(nf_put_att_text(ncid, varid_circumcenter_z, 'coordinates', 9, 'clon clat'))
    !
    !-------------------------------------------
    

    !       write(*,*) "writing def ENDS"
    !       call flush(6)
    CALL nf(nf_enddef(ncid))
    !       write(*,*) "writing def ENDS II"
    !       call flush(6)
    !-------------------------------------------
    
    !-------------------------------------------
    ! put values
    CALL nf(nf_put_var_double(ncid, varid_clon, cells%center(1:no_of_cells)%lon))
    CALL nf(nf_put_var_double(ncid, varid_clat, cells%center(1:no_of_cells)%lat))

    !
    ! vertex part:
    !       write(*,*) "writing verts%vertex..."
    !       call flush(6)
    CALL nf(nf_put_var_double(ncid, varid_vlon,  verts%vertex(1:no_of_verts)%lon))
    CALL nf(nf_put_var_double(ncid, varid_vlat,  verts%vertex(1:no_of_verts)%lat))


    CALL nf(nf_put_var_double(ncid, varid_elon,  edges%center(1:no_of_edges)%lon))
    CALL nf(nf_put_var_double(ncid, varid_elat,  edges%center(1:no_of_edges)%lat))
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !-----------------------------------------------------------------------
    !
    CALL nf(nf_put_var_double(ncid, varid_carea, cells%area(1:no_of_cells)))
    CALL nf(nf_put_var_double(ncid, varid_varea, verts%dual_area(1:no_of_verts)))
    !
    !-----------------------------------------------------------------------
    ! write cell centers
    CALL nf(nf_put_var_double(ncid, varid1,  cells%center(1:no_of_cells)%lon))
    CALL nf(nf_put_var_double(ncid, varid2,  cells%center(1:no_of_cells)%lat))
    
    
    CALL nf(nf_put_var_double(ncid, varid3,  verts%vertex(1:no_of_verts)%lon))
    CALL nf(nf_put_var_double(ncid, varid4,  verts%vertex(1:no_of_verts)%lat))

    CALL nf(nf_put_var_double(ncid, varid5,  edges%center(1:no_of_edges)%lon))
    CALL nf(nf_put_var_double(ncid, varid6,  edges%center(1:no_of_edges)%lat))

    !-----------------------------------------------------------------------
    ! connectivity
    ! allocate tmp arrays for reading
    ALLOCATE(tmp_index(no_of_cells,max_cell_vertices),stat=i)
    IF (i > 0) THEN
      CALL finish ('write_netcdf_grid', 'failed to ALLOCATE(tmp_index)')
    ENDIF
    DO i=1,max_cell_vertices
      tmp_index(1:no_of_cells,i) = cells%get_edge_index(1:no_of_cells,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid7,  tmp_index(1:no_of_cells,:)))
    DO i=1,max_cell_vertices
      tmp_index(1:no_of_cells,i) = cells%get_vertex_index(1:no_of_cells,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid8,  tmp_index(1:no_of_cells,:)))
    DO i=1,max_cell_vertices
      tmp_index(1:no_of_cells,i) = cells%get_neighbor_index(1:no_of_cells,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid26,  tmp_index(1:no_of_cells,:)))
    DO i=1,max_cell_vertices
      tmp_index(1:no_of_cells,i) = cells%get_edge_orient(1:no_of_cells,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid23,  tmp_index(1:no_of_cells,:)))
    DEALLOCATE(tmp_index)

    ALLOCATE(tmp_index(no_of_edges,2),stat=i)
    IF (i > 0) THEN
      CALL finish ('write_netcdf_grid', 'failed to ALLOCATE(tmp_index)')
    ENDIF
    tmp_index(1:no_of_edges,1) = edges%get_cell_index(1:no_of_edges,1)
    tmp_index(1:no_of_edges,2) = edges%get_cell_index(1:no_of_edges,2)
    CALL nf(nf_put_var_int   (ncid, varid9,  tmp_index(1:no_of_edges,:)))
    tmp_index(1:no_of_edges,1) = edges%get_vertex_index(1:no_of_edges,1)
    tmp_index(1:no_of_edges,2) = edges%get_vertex_index(1:no_of_edges,2)
    CALL nf(nf_put_var_int   (ncid, varid10,  tmp_index(1:no_of_edges,:)))
    DEALLOCATE(tmp_index)

    ALLOCATE(tmp_real(no_of_edges,2),stat=i)
    IF (i > 0) THEN
      CALL finish ('write_netcdf_grid', 'failed to ALLOCATE(tmp_real)')
    ENDIF
    tmp_real(1:no_of_edges,1) = edges%get_edge_vert_length(1:no_of_edges,1)
    tmp_real(1:no_of_edges,2) = edges%get_edge_vert_length(1:no_of_edges,2)
    CALL nf(nf_put_var_double(ncid, varid18, tmp_real(1:no_of_edges,:)))
    tmp_real(1:no_of_edges,1) = edges%get_edge_cell_length(1:no_of_edges,1)
    tmp_real(1:no_of_edges,2) = edges%get_edge_cell_length(1:no_of_edges,2)
    CALL nf(nf_put_var_double(ncid, varid40, tmp_real(1:no_of_edges,:)))
    DEALLOCATE(tmp_real)

    ALLOCATE(tmp_index(no_of_verts,max_vert_connect),stat=i)
    IF (i > 0) THEN
      CALL finish ('write_netcdf_grid', 'failed to ALLOCATE(tmp_index)')
    ENDIF
    DO i=1,max_vert_connect
      tmp_index(1:no_of_verts,i) = verts%get_cell_index(1:no_of_verts,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid11, tmp_index(1:no_of_verts,:)))
    DO i=1,max_vert_connect
      tmp_index(1:no_of_verts,i) = verts%get_edge_index(1:no_of_verts,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid12, tmp_index(1:no_of_verts,:)))
    DO i=1,max_vert_connect
      tmp_index(1:no_of_verts,i) = verts%get_neighbor_index(1:no_of_verts,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid13, tmp_index(1:no_of_verts,:)))
    DO i=1,max_vert_connect
      tmp_index(1:no_of_verts,i) = verts%get_edge_orient(1:no_of_verts,i)
    ENDDO
    CALL nf(nf_put_var_int   (ncid, varid30, tmp_index(1:no_of_verts,:)))
    !-----------------------------------------------------------------------
    DEALLOCATE(tmp_index)


    CALL nf(nf_put_var_double(ncid, varid14, cells%area(1:no_of_cells)))
    CALL nf(nf_put_var_double(ncid, varid_cell_elevation, &
       cells%elevation(1:no_of_cells)))
    CALL nf(nf_put_var_int   (ncid, varid_cell_sea_land_mask,&
      cells%sea_land_mask(1:no_of_cells)))
    CALL nf(nf_put_var_int   (ncid, varid_cell_domain_id,&
      cells%domain_id(:,1:no_of_cells)))
    CALL nf(nf_put_var_int   (ncid, varid_cell_no_of_domains,&
      cells%no_of_domains(:)))
    CALL nf(nf_put_var_double(ncid, varid15, verts%dual_area(1:no_of_verts)))
    CALL nf(nf_put_var_double(ncid, varid16, edges%primal_edge_length(1:no_of_edges)))
    CALL nf(nf_put_var_double(ncid, varid17, edges%dual_edge_length(1:no_of_edges)))
    CALL nf(nf_put_var_double(ncid, varid_edge_quad, edges%quad_area(1:no_of_edges)))


    CALL nf(nf_put_var_double(ncid, varid19, edges%primal_normal(1:no_of_edges)%v1))
    CALL nf(nf_put_var_double(ncid, varid20, edges%primal_normal(1:no_of_edges)%v2))
    CALL nf(nf_put_var_double(ncid, varid21, edges%dual_normal(1:no_of_edges)%v1))
    CALL nf(nf_put_var_double(ncid, varid22, edges%dual_normal(1:no_of_edges)%v2) )
    CALL nf(nf_put_var_double(ncid, varid_edge_elevation, edges%elevation(1:no_of_edges)))
    CALL nf(nf_put_var_int   (ncid, varid_edge_sea_land_mask, &
      & edges%sea_land_mask(1:no_of_edges)))

    CALL nf(nf_put_var_int   (ncid, varid24, cells%idx(1:no_of_cells)))
    CALL nf(nf_put_var_int   (ncid, varid25, cells%parent_index(1:no_of_cells)))
    CALL nf(nf_put_var_int   (ncid, varid251, cells%parent_child_type(1:no_of_cells)))

    CALL nf(nf_put_var_int   (ncid, varid27, cells%child_index(1:no_of_cells,:)))
    CALL nf(nf_put_var_int   (ncid, varid41, cells%child_id(1:no_of_cells)))
    CALL nf(nf_put_var_int   (ncid, varid28, edges%idx(1:no_of_edges)))
    !      CALL nf(nf_put_var_int   (ncid, varid281, edges%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid282, edges%parent_child_type(1:no_of_edges)))
    CALL nf(nf_put_var_int   (ncid, varid29, verts%idx(1:no_of_verts)))
    CALL nf(nf_put_var_int   (ncid, varid31, edges%system_orientation(1:no_of_edges)))

    ! Variables added for mesh refinement
    CALL nf(nf_put_var_int   (ncid, varid38, edges%parent_index(1:no_of_edges)))
    CALL nf(nf_put_var_int   (ncid, varid39, edges%child_index(1:no_of_edges,:)))
    CALL nf(nf_put_var_int   (ncid, varid42, edges%child_id(1:no_of_edges)))
    CALL nf(nf_put_var_int   (ncid, varid32, cells%refin_ctrl(1:no_of_cells)))
    ! CALL nf(nf_put_var_int   (ncid, varid33, cells%indlist))
    CALL nf(nf_put_var_int   (ncid, varid43, cells%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid44, cells%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid34, edges%refin_ctrl(1:no_of_edges)))
    ! CALL nf(nf_put_var_int   (ncid, varid35, edges%indlist))
    CALL nf(nf_put_var_int   (ncid, varid45, edges%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid46, edges%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid36, verts%refin_ctrl(1:no_of_verts)))
    ! CALL nf(nf_put_var_int   (ncid, varid37, verts%indlist))
    CALL nf(nf_put_var_int   (ncid, varid47, verts%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid48, verts%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid403,verts%parent_index(1:no_of_verts)))
    !------------------------------------------------------------------------
    CALL nf(nf_put_var_int   (ncid, varid_phys_edge_id, edges%subgrid_id(1:no_of_edges)))
    CALL nf(nf_put_var_int   (ncid, varid_phys_cell_id, cells%subgrid_id(1:no_of_cells)))

    !------------------------------------------------------------------------
    IF (grid_obj%netcdf_flags == netcdf_CF_1_1_convention) THEN
      !------------------------------------------------------------------------
      ! Transpose of index array necessary for CF-1.1 Convention
      !
      ALLOCATE(zv2dx(max_cell_vertices, no_of_cells), zv2dy(max_cell_vertices, no_of_cells))

      IF (max_cell_vertices == 3) THEN
        DO i = 1, no_of_cells
          pole_index=0
          DO j = 1, max_cell_vertices
            zv2dy(j,i) = verts%vertex(cells%get_vertex_index(i,j))%lat
            zv2dx(j,i) = verts%vertex(cells%get_vertex_index(i,j))%lon
            !check if we have a pole
            IF (lat_is_pole(zv2dy(j,i))) &               
               pole_index=j
!                zv2dy(j,i) = SIGN(pi_2,zv2dy(j,i))
          ENDDO

          IF (pole_index > 0) THEN
             SELECT CASE (pole_index)
             CASE (1)
                      j1=2
                      j2=3
             CASE (2)
                      j1=1
                      j2=3
             CASE DEFAULT
                      j1=1
                      j2=2
             END SELECT
             ! wrap around if the center lon is at 180 
             IF (ABS(zv2dx(j1,i) - zv2dx(j2,i)) > pi) THEN
               IF (zv2dx(j1,i) < 0.0_wp) THEN
                 zv2dx(j1,i) = zv2dx(j1,i) + 2.0_wp * pi
               ELSE
                 zv2dx(j2,i) = zv2dx(j2,i) + 2.0_wp * pi
               ENDIF
             ENDIF
             zv2dx(pole_index,i) = (zv2dx(j1,i) + zv2dx(j2,i)) * 0.5_wp
          ENDIF
        ENDDO
      ELSE
        DO i = 1, no_of_cells
          DO j = 1, max_cell_vertices
            IF (cells%get_vertex_index(i,j) > 0) THEN
              zv2dx(j,i) = verts%vertex(cells%get_vertex_index(i,j))%lon
              zv2dy(j,i) = verts%vertex(cells%get_vertex_index(i,j))%lat
            ELSE
              zv2dx(j,i) = fildoub
              zv2dy(j,i) = fildoub
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      WHERE (ABS(zv2dx(:,:)) < EPSILON(0.0_wp))
        zv2dx(:,:) = 0.0_wp
      ENDWHERE
      WHERE (ABS(zv2dy(:,:)) < EPSILON(0.0_wp))
        zv2dy(:,:) = 0.0_wp
      ENDWHERE
!       DO i = 1, no_of_cells
!         DO j = 1, max_cell_vertices
!           IF (ABS(zv2dy(j,i)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
!             zv2dx(j,i) = cells%center(i)%lon
!           ENDIF
!         ENDDO
!       ENDDO
      CALL nf(nf_put_var_double(ncid, varid_clatv, zv2dy))
      CALL nf(nf_put_var_double(ncid, varid_clonv, zv2dx))
      DEALLOCATE(zv2dx, zv2dy)
      !------------------------------------------------------------------------
      ! Transpose of index array necessary for CF-1.1 Convention
      ! Disabled for patches
      ! Has to be rewritten
!       IF (max_vert_connect == 6) THEN ! this works only for triangles
      IF (max_vert_connect == 0) THEN ! this works only for triangles
        ALLOCATE(zv2dx(max_vert_connect, no_of_verts), zv2dy(max_vert_connect, no_of_verts))
        DO i = 1, no_of_verts
          DO j = 1, max_vert_connect
            IF ((verts%get_cell_index(i,j) == 0) .and. (verts%refin_ctrl(i) /= 1)) THEN
              zv2dx(j,i) = cells%center(verts%get_cell_index(i,5))%lon
              zv2dy(j,i) = cells%center(verts%get_cell_index(i,5))%lat
             ELSE IF ((verts%get_cell_index(i,j) < 0) .or. (verts%refin_ctrl(i) == 1)) THEN
              zv2dx(j,i) = 0._wp
              zv2dy(j,i) = 0._wp
            ELSE
              zv2dx(j,i) = cells%center(verts%get_cell_index(i,j))%lon
              zv2dy(j,i) = cells%center(verts%get_cell_index(i,j))%lat
            ENDIF
          ENDDO
        ENDDO
        CALL nf(nf_put_var_double(ncid, varid_vlonv, zv2dx))
        CALL nf(nf_put_var_double(ncid, varid_vlatv, zv2dy))
        DEALLOCATE (zv2dx, zv2dy)
      ENDIF
      !------------------------------------------------------------------------
      ! Transpose of index array necessary for CF-1.1 Convention
      !
      ALLOCATE(zv2dx(4, no_of_edges), zv2dy(4, no_of_edges))
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j)
      DO i = 1, no_of_edges
        ! fill lons
        zv2dx(1,i) = verts%vertex(edges%get_vertex_index(i,1))%lon
        zv2dy(1,i) = verts%vertex(edges%get_vertex_index(i,1))%lat
        zv2dx(3,i) = verts%vertex(edges%get_vertex_index(i,2))%lon
        zv2dy(3,i) = verts%vertex(edges%get_vertex_index(i,2))%lat
        IF (edges%get_cell_index(i,1) > 0) THEN
          zv2dx(4,i) = cells%center(edges%get_cell_index(i,1))%lon
          zv2dy(4,i) = cells%center(edges%get_cell_index(i,1))%lat
        ELSE
          zv2dx(4,i) = 0._wp
          zv2dy(4,i) = 0._wp
        ENDIF
        IF (edges%get_cell_index(i,2) > 0) THEN
          zv2dx(2,i) = cells%center(edges%get_cell_index(i,2))%lon
          zv2dy(2,i) = cells%center(edges%get_cell_index(i,2))%lat
        ELSE
          zv2dx(2,i) = 0._wp
          zv2dy(2,i) = 0._wp
        ENDIF
        DO j = 1, 4
          IF ( lat_is_pole(zv2dy(j,i))) THEN
            zv2dx(j,i) = edges%center(i)%lon
!             zv2dy(j,i) = SIGN(pi_2,zv2dy(j,i))
          ENDIF
        ENDDO

      ENDDO
!$OMP END DO
!       WHERE (ABS(zv2dx(:,:)) < EPSILON(0.0_wp))
!         zv2dx(:,:) = 0.0_wp
!       ENDWHERE
!       WHERE (ABS(zv2dy(:,:)) < EPSILON(0.0_wp))
!         zv2dy(:,:) = 0.0_wp
!       ENDWHERE
!$OMP DO PRIVATE(i,swap)
      DO i = 1, no_of_edges
        IF ( check_orientation(edges%center(i)%lon, &
          & zv2dx(:,i),zv2dy(:,i),4) < 0 ) THEN
          swap(1:4) = zv2dx(4:1:-1,i)
          zv2dx(:,i) = swap(:)
          swap(1:4) = zv2dy(4:1:-1,i)
          zv2dy(:,i) = swap(:)
        ENDIF
!         IF (check_orientation(edges%center(i)%lon, &
!           & zv2dx(:,i),zv2dy(:,i),4) < 0) THEN
!         ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      !
      CALL nf(nf_put_var_double(ncid, varid_elonv, zv2dx))
      CALL nf(nf_put_var_double(ncid, varid_elatv, zv2dy))
      !
      DEALLOCATE (zv2dx, zv2dy)
      !
      !--------------------------------------------------------------------------------------
    ENDIF ! (grid_obj%netcdf_flags == netcdf_CF_1_1_convention)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    ! write cartesian positions

    ! vertexes
    CALL nf(nf_put_var_double(ncid, varid_vx, verts%cartesian(1:no_of_verts)%x(1)))
    CALL nf(nf_put_var_double(ncid, varid_vy, verts%cartesian(1:no_of_verts)%x(2)))
    CALL nf(nf_put_var_double(ncid, varid_vz, verts%cartesian(1:no_of_verts)%x(3)))

    !  edges
    CALL nf(nf_put_var_double(ncid, varid_edge_x, edges%cartesian_center(1:no_of_edges)%x(1)))
    CALL nf(nf_put_var_double(ncid, varid_edge_y, edges%cartesian_center(1:no_of_edges)%x(2)))
    CALL nf(nf_put_var_double(ncid, varid_edge_z, edges%cartesian_center(1:no_of_edges)%x(3)))
    
    CALL nf(nf_put_var_double(ncid, varid_dualedge_x, &
      & edges%dual_cartesian_center(1:no_of_edges)%x(1)))
    CALL nf(nf_put_var_double(ncid, varid_dualedge_y, &
      & edges%dual_cartesian_center(1:no_of_edges)%x(2)))
    CALL nf(nf_put_var_double(ncid, varid_dualedge_z, &
      & edges%dual_cartesian_center(1:no_of_edges)%x(3)))
   
    CALL nf(nf_put_var_double(ncid, varid_edgenormal_x, &
      & edges%cartesian_primal_normal(1:no_of_edges)%x(1)))
    CALL nf(nf_put_var_double(ncid, varid_edgenormal_y, &
      & edges%cartesian_primal_normal(1:no_of_edges)%x(2)))
    CALL nf(nf_put_var_double(ncid, varid_edgenormal_z, &
      & edges%cartesian_primal_normal(1:no_of_edges)%x(3)))
    
    CALL nf(nf_put_var_double(ncid, varid_dualedgenormal_x, &
      & edges%cartesian_dual_normal(1:no_of_edges)%x(1)))
    CALL nf(nf_put_var_double(ncid, varid_dualedgenormal_y, &
      & edges%cartesian_dual_normal(1:no_of_edges)%x(2)))
    CALL nf(nf_put_var_double(ncid, varid_dualedgenormal_z, &
      & edges%cartesian_dual_normal(1:no_of_edges)%x(3)))
    

    ! cells
    CALL nf(nf_put_var_double(ncid, varid_circumcenter_x, &
      & cells%cartesian_center(1:no_of_cells)%x(1)))
    CALL nf(nf_put_var_double(ncid, varid_circumcenter_y, &
      & cells%cartesian_center(1:no_of_cells)%x(2)))
    CALL nf(nf_put_var_double(ncid, varid_circumcenter_z, &
      & cells%cartesian_center(1:no_of_cells)%x(3)))
    !------------------------------------------------------------------------
         
    !------------------------------------------------------------------------
    ! write lon lat of cell barycenters
    CALL nf(nf_put_var_double(ncid, varid_cell_barycenter_lon, &
      & cells%barycenter(1:no_of_cells)%lon))
    CALL nf(nf_put_var_double(ncid, varid_cell_barycenter_lat, &
      & cells%barycenter(1:no_of_cells)%lat))    

    CALL nf(nf_close(ncid))
    !------------------------------------------------------------------------
!     write(*,*) '---',TRIM(grid_obj%file_name), '---'
!     str_idx=LBOUND(verts%start_idx, 1)
!     end_idx=str_idx+SIZE(verts%start_idx, 1)-1
!     DO i=str_idx,end_idx
!       write(*,*) 'verts%start_idx, end:', i, verts%start_idx(i,1), verts%end_idx(i,1)
!     ENDDO
! 
!     str_idx=LBOUND(edges%start_idx, 1)
!     end_idx=str_idx+SIZE(edges%start_idx, 1)-1
!     DO i=str_idx,end_idx
!       write(*,*) 'edges%start_idx, end:', i, edges%start_idx(i,1), edges%end_idx(i,1)
!     ENDDO
! 
!     str_idx=LBOUND(cells%start_idx, 1)
!     end_idx=str_idx+SIZE(cells%start_idx, 1)-1
!     DO i=str_idx,end_idx
!       write(*,*) 'cells%start_idx, end:', i, cells%start_idx(i,1), cells%end_idx(i,1)
!     ENDDO
!     write(*,*) '-------------------'


  END SUBROUTINE write_netcdf_grid
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  INTEGER FUNCTION write_geometry_info(ncid, geometry_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_grid_geometry_info) :: geometry_info

    write_geometry_info = -1
    
    CALL nf(nf_put_att_int      (ncid, nf_global, 'grid_geometry', nf_int, 1,     &
      & geometry_info%geometry_type))
    
    CALL nf(nf_put_att_int   (ncid, nf_global, 'mean_edge_length' , nf_double, 1, &
      & geometry_info%mean_edge_length))

    CALL nf(nf_put_att_double(ncid, nf_global, 'mean_edge_length' , nf_double, 1, &
      & geometry_info%mean_edge_length))
      
    CALL nf(nf_put_att_double  (ncid, nf_global, 'mean_cell_area' , nf_double, 1, &
      & geometry_info%mean_cell_area))
      
    CALL nf(nf_put_att_double   (ncid, nf_global, 'domain_length' , nf_double, 1, &
      & geometry_info%domain_length))
    
    CALL nf(nf_put_att_double   (ncid, nf_global, 'domain_height' , nf_double, 1, &
      & geometry_info%domain_height))
    
    CALL nf(nf_put_att_double   (ncid, nf_global, 'sphere_radius' , nf_double, 1, &
      & geometry_info%sphere_radius))
    
    CALL nf(nf_put_att_double  (ncid, nf_global, 'domain_cartesian_center', nf_double, 3, &
      & geometry_info%center%x))
      
    write_geometry_info = 0

  END FUNCTION write_geometry_info
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION check_orientation (lonc, lon, lat, n)
    INTEGER, INTENT(in) :: n
    REAL(wp), INTENT(in) :: lonc
    REAL(wp), INTENT(in) :: lon(n), lat(n)

    REAL(wp) :: lonl(n), latl(n)

    REAL(wp) :: area

    INTEGER :: i,j

    lonl(:) = lon(:)
    latl(:) = lat(:)

    DO i = 1, n
      lonl(i) = lonl(i) - lonc
      IF (lonl(i) < -pi) THEN
        lonl(i) =  pi+MOD(lonl(i), pi)
      ENDIF
      IF (lonl(i) >  pi) THEN
        lonl(i) = -pi+MOD(lonl(i), pi)
      ENDIF
    ENDDO

    area = 0.0_wp
    DO i = 1, n
      j = MOD(i,n)+1
      area = area+lonl(i)*latl(j)
      area = area-latl(i)*lonl(j)
    ENDDO

    IF (area >= 0.0_wp) THEN
      check_orientation = +1
    ELSE
      check_orientation = -1
    END IF

  END FUNCTION check_orientation
  !-------------------------------------------------------------------------

END MODULE mo_io_local_grid



