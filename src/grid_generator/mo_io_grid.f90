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
MODULE mo_io_grid

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2004-5
  !
  !-------------------------------------------------------------------------
  !
  !
  !

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message_text, message, finish
!  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY: re
  USE mo_grid,               ONLY: t_grid, construct_grid
  USE mo_grid_levels,        ONLY: itype_optimize, l_c_grid
  USE mo_base_geometry,      ONLY: x_rot_angle, y_rot_angle, z_rot_angle
  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell, &
    & min_rlvert, max_rlvert, &
    & min_rledge, max_rledge

  IMPLICIT NONE

  PUBLIC nf

  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: input_grid, write_grid,write_edges_grid


  !  INTEGER, PARAMETER :: unit1=91, unit2=92
  !
  !  INTEGER, PARAMETER, PUBLIC :: binary = 0
  !  INTEGER, PARAMETER, PUBLIC :: netcdf3 = 1


  !--------------------------------------------------------------------
  !BOC

CONTAINS


  !EOC
  !-------------------------------------------------------------------------
  !
  !>
  !!               Reads files with cross reference tables and coordinates of.
  !!
  !!               Reads files with cross reference tables and coordinates of
  !! grid items required by the model. For all grids, the file
  !! <b>GRIDMAP.*</b> is read, where * denotes the grid level.
  !! The files are sequential access, unformatted output files in the directory
  !! <b>input</b>.
  !!
  !! @par Revision History
  !!   Developed and tested  by L.Bonaventura  and T. Heinze (2004-5).
  !!
  SUBROUTINE input_grid(gg, input)
    !
    TYPE(t_grid),                  INTENT(inout) :: gg
    CHARACTER(LEN=filename_max), INTENT(in)    :: input

    INTEGER:: i_nc, i_ne, i_nv

    INTEGER :: ncid, dimid, varid

    !-------------------------------------------------------------------------
    !BOC


    WRITE(message_text,'(a,a)') 'Read gridmap file ', TRIM(input)
    CALL message ('', TRIM(message_text))

    CALL nf(nf_open(TRIM(input), nf_nowrite, ncid))

    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_level', gg%level))

    CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, gg%ncells))
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, gg%nedges))
    CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, gg%nverts))

    i_nc = gg%ncells
    i_ne = gg%nedges
    i_nv = gg%nverts

    CALL construct_grid (gg, i_nc, i_ne, i_nv)

    CALL nf(nf_inq_varid(ncid, 'lon_cell_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%cells%center(:)%lon))
    CALL nf(nf_inq_varid(ncid, 'lat_cell_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%cells%center(:)%lat))
    CALL nf(nf_inq_varid(ncid, 'longitude_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%verts%vertex(:)%lon))
    CALL nf(nf_inq_varid(ncid, 'latitude_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%verts%vertex(:)%lat))
    CALL nf(nf_inq_varid(ncid, 'lon_edge_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%center(:)%lon))
    CALL nf(nf_inq_varid(ncid, 'lat_edge_centre', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%center(:)%lat))
    CALL nf(nf_inq_varid(ncid, 'edge_of_cell', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%edge_index))
    CALL nf(nf_inq_varid(ncid, 'vertex_of_cell', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%vertex_index))
    CALL nf(nf_inq_varid(ncid, 'adjacent_cell_of_edge', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%cell_index))
    CALL nf(nf_inq_varid(ncid, 'edge_vertices', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%vertex_index))
    CALL nf(nf_inq_varid(ncid, 'cells_of_vertex', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%verts%cell_index))
    CALL nf(nf_inq_varid(ncid, 'edges_of_vertex', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%verts%edge_index))
    CALL nf(nf_inq_varid(ncid, 'vertices_of_vertex', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%verts%neighbor_index))
    CALL nf(nf_inq_varid(ncid, 'cell_area_p', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%cells%area))
    CALL nf(nf_inq_varid(ncid, 'dual_area_p', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%verts%dual_area))
    CALL nf(nf_inq_varid(ncid, 'edge_length', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%primal_edge_length))
    CALL nf(nf_inq_varid(ncid, 'dual_edge_length', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%dual_edge_length))
    CALL nf(nf_inq_varid(ncid, 'edge_vert_distance', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%edge_vert_length))
    CALL nf(nf_inq_varid(ncid, 'edge_cell_distance', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%edge_cell_length))
    CALL nf(nf_inq_varid(ncid, 'zonal_normal_primal_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%primal_normal(:)%v1))
    CALL nf(nf_inq_varid(ncid, 'meridional_normal_primal_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%primal_normal(:)%v2))
    CALL nf(nf_inq_varid(ncid, 'zonal_normal_dual_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%dual_normal(:)%v1))
    CALL nf(nf_inq_varid(ncid, 'meridional_normal_dual_edge', varid))
    CALL nf(nf_get_var_double(ncid, varid, gg%edges%dual_normal(:)%v2))
    CALL nf(nf_inq_varid(ncid, 'orientation_of_normal', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%edge_orientation))
    CALL nf(nf_inq_varid(ncid, 'cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%idx))
    CALL nf(nf_inq_varid(ncid, 'parent_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%parent_index))
    CALL nf(nf_inq_varid(ncid, 'parent_cell_type', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%parent_child_type))
    CALL nf(nf_inq_varid(ncid, 'neighbor_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%neighbor_index))
    CALL nf(nf_inq_varid(ncid, 'child_cell_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%cells%child_index))
    CALL nf(nf_inq_varid(ncid, 'edge_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%idx))
    CALL nf(nf_inq_varid(ncid, 'child_edge_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%child_index))
    CALL nf(nf_inq_varid(ncid, 'parent_edge_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%parent_index))
    CALL nf(nf_inq_varid(ncid, 'edge_parent_type', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%parent_child_type))
    CALL nf(nf_inq_varid(ncid, 'vertex_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%verts%idx))
    CALL nf(nf_inq_varid(ncid, 'parent_vertex_index', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%verts%parent_index))

    CALL nf(nf_inq_varid(ncid, 'edge_orientation', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%verts%edge_orientation))
    CALL nf(nf_inq_varid(ncid, 'edge_system_orientation', varid))
    CALL nf(nf_get_var_int   (ncid, varid, gg%edges%system_orientation))

    CALL nf(nf_close(ncid))

    !     DO i_nc=1,gg%nedges
    !        WRITE(*,*) gg%edges%parent_index(i_nc)
    !     ENDDO
    ! these are not filled in the file

  END SUBROUTINE input_grid


  SUBROUTINE nf(STATUS)
    INTEGER, INTENT(in) :: STATUS

    IF (STATUS /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(STATUS))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  ! !IROUTINE: output_grid
  !
  ! !SUBROUTINE INTERFACE:
  SUBROUTINE write_grid(gg, outfile)
    !
    ! !DESCRIPTION:
    ! Creates files  with  cross reference tables and coordinates
    ! of grid items to be read by the model. In most cases for all grid levels,
    ! the file {\bf GRIDMAP.*} is produced, where * denotes the grid level.
    ! The exception is Heikes-Randall + C grid optimization in the last step. Here
    ! only the {\bf GRIDMAP} file of the last level is created. The background
    ! is in this case the {\bf GRIDMAP} files of all the coarser levels would
    ! produce just Heikes-Randall optimized grid without the C grid optimization.
    ! Only in this case gridgen has to be run for all grid levels.\\
    ! The files are produced as sequential access, unformatted output.
    !
    ! !REVISION HISTORY:
    ! Leonidas Lianrdakis, MPI-M (2004-5).
    !
    TYPE(t_grid), INTENT(in) :: gg
    CHARACTER(LEN=filename_max), INTENT(in) :: outfile

    INTEGER :: i_nc, i_ne, i_nv   ! number of cells, edges and vertices

    INTEGER :: old_mode

    INTEGER :: ncid
    INTEGER :: dimids(2)

    INTEGER :: dim_ncell, dim_nvertex, dim_nedge, dim_two, dim_cell_refine, &
      & dim_edge_refine, dim_vert_refine, dim_nvertex_per_cell,      &
      & dim_ncells_per_edge, dim_nedges_per_vertex, dim_nchilds_per_cell, &
      & dim_list, dim_nchdom

    INTEGER :: varid_clon, varid_clat, varid_clonv, varid_clatv
    INTEGER :: varid_vlon, varid_vlat, varid_vlonv, varid_vlatv
    INTEGER :: varid_elon, varid_elat, varid_elonv, varid_elatv

    INTEGER :: varid_carea, varid_varea

    INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8, &
      & varid9, varid10, varid11, varid12, varid13, varid14, varid15,   &
      & varid16, varid17, varid18, varid19, varid20, varid21, varid22,  &
      & varid23, varid24, varid25, varid26, varid27, varid28, varid29,  &
      & varid30, varid31, varid32, varid33, varid34, varid35, varid36,  &
      & varid37, varid38, varid39, varid40, varid41, varid42, varid43,  &
      & varid44, varid45, varid46, varid47, varid48, varid251,           &
      & varid282, varid403, varid49, varid50

    INTEGER :: varid_cell_elevation, varid_cell_sea_land_mask

    REAL(wp), POINTER :: double_pnt_1d(:)
    INTEGER,  POINTER :: int_pnt_1d(:)
    !  INTEGER :: i, j

#ifdef __SX__
    INTEGER :: iargc
    CHARACTER(LEN= 32) :: arg_str
#endif

    !  REAL(wp), ALLOCATABLE :: zv2d(:,:), zv2dx(:,:), zv2dy(:,:)

    REAL(wp) :: rotation_vector(3)

    !  REAL(wp) :: swap(4)

    INTEGER :: ilevel, grid_root

    !EOP
    !-------------------------------------------------------------------------

    ! distinguish between gridgeneration for optimization strategies

    ! Dummy settings for special grids


    grid_root = 0
    ilevel    = 1
    WRITE(message_text,'(a,a)') 'Write gridmap file: ', TRIM(outfile)
    CALL message ('', TRIM(message_text))


    i_nc = gg%ncells
    i_ne = gg%nedges
    i_nv = gg%nverts

    rotation_vector = (/ x_rot_angle, y_rot_angle, z_rot_angle /)

    !----------------------------------------------------------------------
    !
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))

    CALL nf(nf_create(TRIM(outfile), nf_clobber, ncid))
    CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))

    !
    !----------------------------------------------------------------------
    !
    ! Global attributes
    !
    CALL nf(nf_put_att_text    (ncid, nf_global, 'title', 21, 'ICON grid description'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'institution', 59, &
      & 'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'source', 10, 'icon-dev'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'grid_mapping_name' , 18, 'lat_long_on_sphere'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_id' , 28, 'urn:ogc:def:cs:EPSG:6.0:6422'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'crs_name',30,'Spherical 2D Coordinate System'))
    CALL nf(nf_put_att_text    (ncid, nf_global, 'ellipsoid_name' , 6, 'Sphere'))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'semi_major_axis' , nf_double, 1, re))
    CALL nf(nf_put_att_double  (ncid, nf_global, 'inverse_flattening' , nf_double, 1, 0.0_wp))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_level', nf_int, 1, ilevel))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_root', nf_int, 1, grid_root))
    ! The following three attributes are nontrivial in the presence of grid refinement
    ! and are set here because they are checked in the ICON code
    CALL nf(nf_put_att_int     (ncid, nf_global, 'grid_ID', nf_int, 1, 1))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'parent_grid_ID', nf_int, 1, 0))
    CALL nf(nf_put_att_int     (ncid, nf_global, 'max_childdom', nf_int, 1, 1))
    SELECT CASE (itype_optimize)
    CASE (0)
      CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 12, 'natural grid'))
    CASE (1)
      IF (l_c_grid) THEN
        CALL nf(nf_put_att_text(ncid, nf_global, 'grid_optimization', 50, &
          & 'Heikes-Randall with additional c-grid optimization'))
      ELSE
        CALL nf(nf_put_att_text(ncid, nf_global, 'grid_optimization', 27, &
          & 'Heikes-Randall optimization'))
      ENDIF
    CASE (2)
      CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 22, &
        & 'equal area subdivision'))
    CASE (3)
      CALL nf(nf_put_att_text  (ncid, nf_global, 'grid_optimization', 30, &
        & 'c-grid small circle constraint'))
    END SELECT
    CALL nf(nf_put_att_double  (ncid, nf_global, 'rotation_vector',nf_double,3,rotation_vector))
    !
    !      write(*,*) "writing Dimensions..."
    !       call flush(6)
    ! Dimensions
    !
    CALL nf(nf_def_dim(ncid, 'cell',   i_nc, dim_ncell))
    CALL nf(nf_def_dim(ncid, 'vertex', i_nv, dim_nvertex))
    CALL nf(nf_def_dim(ncid, 'edge',   i_ne, dim_nedge))
    !
    CALL nf(nf_def_dim(ncid, 'nc',        2, dim_ncells_per_edge))
    CALL nf(nf_def_dim(ncid, 'nv',  gg%cells%max_no_of_vertices, dim_nvertex_per_cell))
    CALL nf(nf_def_dim(ncid, 'ne',  gg%verts%max_connectivity, dim_nedges_per_vertex))
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

    !
    !---------------------------------------------------------------------
    !
    ! Grid variables
    !
    !---------------------------------------------------------------------
    !
    ! public grid information
    !
    ! cell part:
    !

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
    ! write(*,*) "writing private grid information ..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'lon_cell_centre', nf_double, 1, dim_ncell, varid1))
    CALL nf(nf_put_att_text(ncid, varid1, 'long_name', 24, 'longitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid1, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid1, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_cell_centre', nf_double, 1, dim_ncell, varid2))
    CALL nf(nf_put_att_text(ncid, varid2, 'long_name', 23, 'latitude of cell centre'))
    CALL nf(nf_put_att_text(ncid, varid2, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid2, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'longitude_vertices', nf_double, 1, dim_nvertex, varid3))
    CALL nf(nf_put_att_text(ncid, varid3, 'long_name', 21, 'longitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid3, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid3, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'latitude_vertices', nf_double, 1, dim_nvertex, varid4))
    CALL nf(nf_put_att_text(ncid, varid4, 'long_name', 20, 'latitude of vertices'))
    CALL nf(nf_put_att_text(ncid, varid4, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid4, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lon_edge_centre', nf_double, 1, dim_nedge, varid5))
    CALL nf(nf_put_att_text(ncid, varid5, 'long_name', 28, 'longitudes of edge midpoints'))
    CALL nf(nf_put_att_text(ncid, varid5, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid5, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'lat_edge_centre', nf_double, 1, dim_nedge, varid6))
    CALL nf(nf_put_att_text(ncid, varid6, 'long_name', 27, 'latitudes of edge midpoints'))
    CALL nf(nf_put_att_text(ncid, varid6, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid6, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'edge_of_cell', nf_int, 2, dimids, varid7))
    CALL nf(nf_put_att_text(ncid, varid7, 'long_name', 29, 'edges of each triangular cell'))
    CALL nf(nf_put_att_text(ncid, varid7, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'vertex_of_cell', nf_int, 2, dimids, varid8))
    CALL nf(nf_put_att_text(ncid, varid8, 'long_name', 32, 'vertices of each triangular cell'))
    CALL nf(nf_put_att_text(ncid, varid8, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'adjacent_cell_of_edge', nf_int, 2, dimids, varid9))
    CALL nf(nf_put_att_text(ncid, varid9, 'long_name', 27, 'cells adjacent to each edge'))
    CALL nf(nf_put_att_text(ncid, varid9, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vertices', nf_int, 2, dimids, varid10))
    CALL nf(nf_put_att_text(ncid, varid10, 'long_name', 35, &
      & 'vertices at the end of of each edge'))
    CALL nf(nf_put_att_text(ncid, varid10, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'cells_of_vertex', nf_int, 2, dimids, varid11))
    CALL nf(nf_put_att_text(ncid, varid11, 'long_name', 24, 'cells around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid11, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'edges_of_vertex', nf_int, 2, dimids, varid12))
    CALL nf(nf_put_att_text(ncid, varid12, 'long_name', 24, 'edges around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid12, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'vertices_of_vertex', nf_int, 2, dimids, varid13))
    CALL nf(nf_put_att_text(ncid, varid13, 'long_name', 27, 'vertices around each vertex'))
    CALL nf(nf_put_att_text(ncid, varid13, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_area_p', nf_double, 1, dim_ncell, varid14))
    CALL nf(nf_put_att_text(ncid, varid14, 'long_name', 17, 'area of grid cell'))
    CALL nf(nf_put_att_text(ncid, varid14, 'units', 2, 'm2'))
    CALL nf(nf_put_att_text(ncid, varid14, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_elevation', nf_double, 1, dim_ncell, varid_cell_elevation))
    CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'long_name', 29, &
      & 'elevation at the cell centers'))
    CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid_cell_elevation, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_sea_land_mask', nf_int,1,dim_ncell,varid_cell_sea_land_mask))
    CALL nf(nf_put_att_text(ncid, varid_cell_sea_land_mask, 'long_name', 36, &
      & 'sea (-1) land (1) mask for the cells'))
    CALL nf(nf_put_att_text(ncid, varid_cell_sea_land_mask, 'units', 4, '1,-1'))
    CALL nf(nf_put_att_text(ncid, varid_cell_sea_land_mask, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'dual_area_p', nf_double, 1, dim_nvertex, varid15))
    CALL nf(nf_put_att_text(ncid, varid15, 'long_name', 40, &
      & 'areas of dual hexagonal/pentagonal cells'))
    CALL nf(nf_put_att_text(ncid, varid15, 'units', 2, 'm2'))
    CALL nf(nf_put_att_text(ncid, varid15, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_length', nf_double, 1, dim_nedge, varid16))
    CALL nf(nf_put_att_text(ncid, varid16, 'long_name', 36, &
      & 'lengths of edges of triangular cells'))
    CALL nf(nf_put_att_text(ncid, varid16, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid16, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_cell_distance', nf_double, 2, dimids, varid40))
    CALL nf(nf_put_att_text(ncid, varid40, 'long_name', 63, &
      & 'distances between edge midpoint and adjacent triangle midpoints'))
    CALL nf(nf_put_att_text(ncid, varid40, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid40, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'dual_edge_length', nf_double, 1, dim_nedge, varid17))
    CALL nf(nf_put_att_text(ncid, varid17, 'long_name', 71, &
      & 'lengths of dual edges (distances between triangular cell circumcenters)'))
    CALL nf(nf_put_att_text(ncid, varid17, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid17, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vert_distance', nf_double, 2, dimids, varid18))
    CALL nf(nf_put_att_text(ncid, varid18, 'long_name', 57, &
      & 'distances between edge midpoint and vertices of that edge'))
    CALL nf(nf_put_att_text(ncid, varid18, 'units', 1, 'm'))
    CALL nf(nf_put_att_text(ncid, varid18, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'zonal_normal_primal_edge', nf_double, 1, dim_nedge, varid19))
    CALL nf(nf_put_att_text(ncid, varid19, 'long_name', 40, &
      & 'zonal component of normal to primal edge'))
    CALL nf(nf_put_att_text(ncid, varid19, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid19, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'meridional_normal_primal_edge', nf_double, 1, dim_nedge, varid20))
    CALL nf(nf_put_att_text(ncid, varid20, 'long_name', 45, &
      & 'meridional component of normal to primal edge'))
    CALL nf(nf_put_att_text(ncid, varid20, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid20, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'zonal_normal_dual_edge', nf_double, 1, dim_nedge, varid21))
    CALL nf(nf_put_att_text(ncid, varid21, 'long_name', 38, &
      & 'zonal component of normal to dual edge'))
    CALL nf(nf_put_att_text(ncid, varid21, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid21, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'meridional_normal_dual_edge', nf_double, 1, dim_nedge, varid22))
    CALL nf(nf_put_att_text(ncid, varid22, 'long_name', 43, &
      & 'meridional component of normal to dual edge'))
    CALL nf(nf_put_att_text(ncid, varid22, 'units', 6, 'radian'))
    CALL nf(nf_put_att_text(ncid, varid22, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'orientation_of_normal', nf_int, 2, dimids, varid23))
    CALL nf(nf_put_att_text(ncid, varid23, 'long_name', 48, &
      & 'orientations of normals to triangular cell edges'))
    CALL nf(nf_put_att_text(ncid, varid23, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'cell_index', nf_int, 1, dim_ncell, varid24))
    CALL nf(nf_put_att_text(ncid, varid24, 'long_name', 10, 'cell index'))
    CALL nf(nf_put_att_text(ncid, varid24, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_cell_index', nf_int, 1, dim_ncell, varid25))
    CALL nf(nf_put_att_text(ncid, varid25, 'long_name', 17, 'parent cell index'))
    CALL nf(nf_put_att_text(ncid, varid25, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_cell_type', nf_int, 1, dim_ncell, varid251))
    CALL nf(nf_put_att_text(ncid, varid251, 'long_name', 17, 'parent cell type'))
    CALL nf(nf_put_att_text(ncid, varid251, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nvertex_per_cell /)
    CALL nf(nf_def_var(ncid, 'neighbor_cell_index', nf_int, 2, dimids, varid26))
    CALL nf(nf_put_att_text(ncid, varid26, 'long_name', 19, 'cell neighbor index'))
    CALL nf(nf_put_att_text(ncid, varid26, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_ncell, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_cell_index', nf_int, 2, dimids, varid27))
    CALL nf(nf_put_att_text(ncid, varid27, 'long_name', 16, 'child cell index'))
    CALL nf(nf_put_att_text(ncid, varid27, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'child_cell_id', nf_int, 1, dim_ncell, varid41))
    CALL nf(nf_put_att_text(ncid, varid41, 'long_name', 23, 'domain ID of child cell'))
    CALL nf(nf_put_att_text(ncid, varid41, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_index', nf_int, 1, dim_nedge, varid28))
    CALL nf(nf_put_att_text(ncid, varid28, 'long_name', 10, 'edge index'))
    CALL nf(nf_put_att_text(ncid, varid28, 'cdi', 6, 'ignore'))
    !
    !       CALL nf(nf_def_var(ncid, 'edge_parent', NF_INT, 1, dim_nedge, varid281))
    !       CALL nf(nf_put_att_text(ncid, varid281, 'long_name', 10, 'edge parent'))
    !       CALL nf(nf_put_att_text(ncid, varid281, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_parent_type', nf_int, 1, dim_nedge, varid282))
    CALL nf(nf_put_att_text(ncid, varid282, 'long_name', 10, 'edge parent type'))
    CALL nf(nf_put_att_text(ncid, varid282, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'vertex_index', nf_int, 1, dim_nvertex, varid29))
    CALL nf(nf_put_att_text(ncid, varid29, 'long_name', 14, 'vertices index'))
    CALL nf(nf_put_att_text(ncid, varid29, 'cdi', 6, 'ignore'))
    !
    dimids = (/ dim_nvertex, dim_nedges_per_vertex /)
    CALL nf(nf_def_var(ncid, 'edge_orientation', nf_int, 2, dimids, varid30))
    CALL nf(nf_put_att_text(ncid, varid30, 'long_name', 16, 'edge orientation'))
    CALL nf(nf_put_att_text(ncid, varid30, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'edge_system_orientation', nf_int, 1, dim_nedge, varid31))
    CALL nf(nf_put_att_text(ncid, varid31, 'long_name', 23, 'edge system orientation'))
    CALL nf(nf_put_att_text(ncid, varid31, 'cdi', 6, 'ignore'))
    !
    ! Variables added for mesh refinement
    !
    ! write(*,*) "writing variables added for mesh refinement ..."
    !       write(*,*) "dim_list=",dim_list
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'refin_c_ctrl', nf_int, 1, dim_ncell, varid32))
    CALL nf(nf_put_att_text(ncid, varid32, 'long_name', 33,'refinement control flag for cells'))
    CALL nf(nf_put_att_text(ncid, varid32, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def index_c_list ..."
    !       call flush(6)
    dimids = (/ dim_cell_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_c_list', nf_int, 2, dimids, varid33))
    CALL nf(nf_put_att_text(ncid, varid33, 'long_name', 73, &
      & 'list of start and end indices for each refinement control level for cells'))
    CALL nf(nf_put_att_text(ncid, varid33, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def start_idx_c..."
    !       call flush(6)
    dimids = (/ dim_cell_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_c', nf_int, 2, dimids, varid43))
    CALL nf(nf_put_att_text(ncid, varid43, 'long_name', 65, &
      & 'list of start indices for each refinement control level for cells'))
    CALL nf(nf_put_att_text(ncid, varid43, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def end_idx_c..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'end_idx_c', nf_int, 2, dimids, varid44))
    CALL nf(nf_put_att_text(ncid, varid44, 'long_name', 63, &
      & 'list of end indices for each refinement control level for cells'))
    CALL nf(nf_put_att_text(ncid, varid44, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def refin_e_ctrl..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'refin_e_ctrl', nf_int, 1, dim_nedge, varid34))
    CALL nf(nf_put_att_text(ncid, varid34, 'long_name', 33,'refinement control flag for edges'))
    CALL nf(nf_put_att_text(ncid, varid34, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def index_e_list..."
    !       call flush(6)
    dimids = (/ dim_edge_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_e_list', nf_int, 2, dimids, varid35))
    CALL nf(nf_put_att_text(ncid, varid35, 'long_name', 73, &
      & 'list of start and end indices for each refinement control level for edges'))
    CALL nf(nf_put_att_text(ncid, varid35, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def start_idx_e..."
    !       call flush(6)
    dimids = (/ dim_edge_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_e', nf_int, 2, dimids, varid45))
    CALL nf(nf_put_att_text(ncid, varid45, 'long_name', 65, &
      & 'list of start indices for each refinement control level for edges'))
    CALL nf(nf_put_att_text(ncid, varid45, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def end_idx_e..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'end_idx_e', nf_int, 2, dimids, varid46))
    CALL nf(nf_put_att_text(ncid, varid46, 'long_name', 63, &
      & 'list of end indices for each refinement control level for edges'))
    CALL nf(nf_put_att_text(ncid, varid46, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def refin_v_ctrl..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'refin_v_ctrl', nf_int, 1, dim_nvertex, varid36))
    CALL nf(nf_put_att_text(ncid, varid36, 'long_name', 36, &
      & 'refinement control flag for vertices'))
    CALL nf(nf_put_att_text(ncid, varid36, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def index_v_list..."
    !       call flush(6)
    dimids = (/ dim_vert_refine, dim_two /)
    CALL nf(nf_def_var(ncid, 'index_v_list', nf_int, 2, dimids, varid37))
    CALL nf(nf_put_att_text(ncid, varid37, 'long_name', 76, &
      & 'list of start and end indices for each refinement control level for vertices'))
    CALL nf(nf_put_att_text(ncid, varid37, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def start_idx_v..."
    !       call flush(6)
    dimids = (/ dim_vert_refine, dim_nchdom /)
    CALL nf(nf_def_var(ncid, 'start_idx_v', nf_int, 2, dimids, varid47))
    CALL nf(nf_put_att_text(ncid, varid47, 'long_name', 68, &
      & 'list of start indices for each refinement control level for vertices'))
    CALL nf(nf_put_att_text(ncid, varid47, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def end_idx_v..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'end_idx_v', nf_int, 2, dimids, varid48))
    CALL nf(nf_put_att_text(ncid, varid48, 'long_name', 66, &
      & 'list of end indices for each refinement control level for vertices'))
    CALL nf(nf_put_att_text(ncid, varid48, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def parent_edge_index..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'parent_edge_index', nf_int, 1, dim_nedge, varid38))
    CALL nf(nf_put_att_text(ncid, varid38, 'long_name', 17, 'parent edge index'))
    CALL nf(nf_put_att_text(ncid, varid38, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def child_edge_index..."
    !       call flush(6)
    dimids = (/ dim_nedge, dim_nchilds_per_cell /)
    CALL nf(nf_def_var(ncid, 'child_edge_index', nf_int, 2, dimids, varid39))
    CALL nf(nf_put_att_text(ncid, varid39, 'long_name', 16, 'child edge index'))
    CALL nf(nf_put_att_text(ncid, varid39, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def child_edge_id..."
    !       call flush(6)
    CALL nf(nf_def_var(ncid, 'child_edge_id', nf_int, 1, dim_nedge, varid42))
    CALL nf(nf_put_att_text(ncid, varid42, 'long_name', 23, 'domain ID of child edge'))
    CALL nf(nf_put_att_text(ncid, varid42, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'phys_cell_id', nf_int, 1, dim_ncell, varid49))
    CALL nf(nf_put_att_text(ncid, varid49, 'long_name', 26, 'physical domain ID of cell'))
    CALL nf(nf_put_att_text(ncid, varid49, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'phys_edge_id', nf_int, 1, dim_nedge, varid50))
    CALL nf(nf_put_att_text(ncid, varid50, 'long_name', 26, 'physical domain ID of edge'))
    CALL nf(nf_put_att_text(ncid, varid50, 'cdi', 6, 'ignore'))
    !
    CALL nf(nf_def_var(ncid, 'parent_vertex_index', nf_int, 1, dim_nvertex, varid403))
    CALL nf(nf_put_att_text(ncid, varid403, 'long_name', 19, 'parent vertex index'))
    CALL nf(nf_put_att_text(ncid, varid403, 'cdi', 6, 'ignore'))
    !
    !       write(*,*) "writing def ENDS"
    !       call flush(6)
    CALL nf(nf_enddef(ncid))
    !       write(*,*) "writing def ENDS II"
    !       call flush(6)
    !
    !-------------------------------------------------------------------------
    !
    ! cell part:
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !
    !       ALLOCATE(zv2dx(3, i_nc), zv2dy(3, i_nc))
    !       DO j = 1, 3
    !         DO i = 1, i_nc
    !           zv2dx(j,i) = gg%verts%vertex(gg%cells%vertex_index(i,j))%lon
    !         ENDDO
    !       ENDDO
    !       DO j = 1, 3
    !         DO i = 1, i_nc
    !           zv2dy(j,i) = gg%verts%vertex(gg%cells%vertex_index(i,j))%lat
    !         ENDDO
    !       ENDDO
    !       WHERE (ABS(zv2dx(:,:)) < EPSILON(0.0_wp))
    !         zv2dx(:,:) = 0.0_wp
    !       ENDWHERE
    !       WHERE (ABS(zv2dy(:,:)) < EPSILON(0.0_wp))
    !         zv2dy(:,:) = 0.0_wp
    !       ENDWHERE
    !       DO j = 1, i_nc
    !         DO i = 1, 3
    !           IF (ABS(zv2dy(i,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
    !             zv2dx(i,j) = gg%cells%center(j)%lon
    !           ENDIF
    !         ENDDO
    !       ENDDO
    !       CALL nf(nf_put_var_double(ncid, varid_clatv, zv2dy))
    !       CALL nf(nf_put_var_double(ncid, varid_clonv, zv2dx))
    !       DEALLOCATE(zv2dx, zv2dy)

    !       write(*,*) "writing cells%center..."
    !       call flush(6)
    CALL nf(nf_put_var_double(ncid, varid_clon, gg%cells%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_clat, gg%cells%center(:)%lat))

    !
    ! vertex part:
    !
    !       write(*,*) "writing verts%vertex..."
    !       call flush(6)
    CALL nf(nf_put_var_double(ncid, varid_vlon,  gg%verts%vertex(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_vlat,  gg%verts%vertex(:)%lat))
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !
    !       ALLOCATE(zv2d(6, i_nv))
    !       DO j = 1, 6
    !         DO i = 1, i_nv
    !           IF (gg%verts%cell_index(i,j) == 0) THEN
    !             zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,5))%lon
    !           ELSE
    !             zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,j))%lon
    !           ENDIF
    !         ENDDO
    !       ENDDO
    !       CALL nf(nf_put_var_double(ncid, varid_vlonv, zv2d))
    !       DO j = 1, 6
    !         DO i = 1, i_nv
    !           IF (gg%verts%cell_index(i,j) == 0) THEN
    !             zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,5))%lat
    !           ELSE
    !             zv2d(7-j,i) = gg%cells%center(gg%verts%cell_index(i,j))%lat
    !           ENDIF
    !         ENDDO
    !       ENDDO
    !       CALL nf(nf_put_var_double(ncid, varid_vlatv, zv2d))
    !       DEALLOCATE (zv2d)
    !
    ! edge part:
    !
    !       write(*,*) "writing edges%center..."
    !       call flush(6)
    CALL nf(nf_put_var_double(ncid, varid_elon,  gg%edges%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_elat,  gg%edges%center(:)%lat))
    !
    ! Transpose of index array necessary for CF-1.1 Convention
    !
    !       ALLOCATE(zv2dx(4, i_ne), zv2dy(4, i_ne))
    !       DO i = 1, i_ne
    !         zv2dx(1,i) = gg%verts%vertex(gg%edges%vertex_index(i,1))%lon
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dx(3,i) = gg%verts%vertex(gg%edges%vertex_index(i,2))%lon
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dx(4,i) = gg%cells%center(gg%edges%cell_index(i,1))%lon
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dx(2,i) = gg%cells%center(gg%edges%cell_index(i,2))%lon
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dy(1,i) = gg%verts%vertex(gg%edges%vertex_index(i,1))%lat
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dy(3,i) = gg%verts%vertex(gg%edges%vertex_index(i,2))%lat
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dy(4,i) = gg%cells%center(gg%edges%cell_index(i,1))%lat
    !       ENDDO
    !       DO i = 1, i_ne
    !         zv2dy(2,i) = gg%cells%center(gg%edges%cell_index(i,2))%lat
    !       ENDDO
    !       WHERE (ABS(zv2dx(:,:)) < EPSILON(0.0_wp))
    !         zv2dx(:,:) = 0.0_wp
    !       ENDWHERE
    !       WHERE (ABS(zv2dy(:,:)) < EPSILON(0.0_wp))
    !         zv2dy(:,:) = 0.0_wp
    !       ENDWHERE
    !       DO j = 1, i_ne
    !         DO i = 1, 4
    !           IF ( ABS(zv2dy(i,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
    !             zv2dx(i,j) = gg%edges%center(j)%lon
    !           ENDIF
    !         ENDDO
    !       ENDDO
    !       DO i = 1, i_ne
    !         IF (check_orientation(gg%edges%center(i)%lon, &
    !              zv2dx(:,i),zv2dy(:,i),4) < 0) THEN
    !           swap(1:4) = zv2dx(4:1:-1,i)
    !           zv2dx(:,i) = swap(:)
    !           swap(1:4) = zv2dy(4:1:-1,i)
    !           zv2dy(:,i) = swap(:)
    !         ENDIF
    !         IF (check_orientation(gg%edges%center(i)%lon, &
    !              zv2dx(:,i),zv2dy(:,i),4) < 0) THEN
    !         ENDIF
    !       ENDDO
    !       !
    !       CALL nf(nf_put_var_double(ncid, varid_elonv, zv2dx))
    !       CALL nf(nf_put_var_double(ncid, varid_elatv, zv2dy))
    !       !
    !       DEALLOCATE (zv2dx, zv2dy)
    !
    !-----------------------------------------------------------------------
    !
    CALL nf(nf_put_var_double(ncid, varid_carea, gg%cells%area))
    CALL nf(nf_put_var_double(ncid, varid_varea, gg%verts%dual_area))
    !
    !-----------------------------------------------------------------------
    !
    CALL nf(nf_put_var_double(ncid, varid1,  gg%cells%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid2,  gg%cells%center(:)%lat))
    CALL nf(nf_put_var_double(ncid, varid3,  gg%verts%vertex(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid4,  gg%verts%vertex(:)%lat))
    CALL nf(nf_put_var_double(ncid, varid5,  gg%edges%center(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid6,  gg%edges%center(:)%lat))
    CALL nf(nf_put_var_int   (ncid, varid7,  gg%cells%edge_index))
    CALL nf(nf_put_var_int   (ncid, varid8,  gg%cells%vertex_index))
    CALL nf(nf_put_var_int   (ncid, varid9,  gg%edges%cell_index))
    CALL nf(nf_put_var_int   (ncid, varid10, gg%edges%vertex_index))
    CALL nf(nf_put_var_int   (ncid, varid11, gg%verts%cell_index))
    CALL nf(nf_put_var_int   (ncid, varid12, gg%verts%edge_index))
    CALL nf(nf_put_var_int   (ncid, varid13, gg%verts%neighbor_index))
    CALL nf(nf_put_var_double(ncid, varid14, gg%cells%area))
    double_pnt_1d => gg%cells%elevation(1:)
    CALL nf(nf_put_var_double(ncid, varid_cell_elevation, double_pnt_1d))
    int_pnt_1d   => gg%cells%sea_land_mask(1:)
    CALL nf(nf_put_var_int   (ncid, varid_cell_sea_land_mask, int_pnt_1d))
    CALL nf(nf_put_var_double(ncid, varid15, gg%verts%dual_area))
    CALL nf(nf_put_var_double(ncid, varid16, gg%edges%primal_edge_length))
    CALL nf(nf_put_var_double(ncid, varid17, gg%edges%dual_edge_length))
    CALL nf(nf_put_var_double(ncid, varid18, gg%edges%edge_vert_length))
    CALL nf(nf_put_var_double(ncid, varid40, gg%edges%edge_cell_length))
    CALL nf(nf_put_var_double(ncid, varid19, gg%edges%primal_normal(:)%v1))
    CALL nf(nf_put_var_double(ncid, varid20, gg%edges%primal_normal(:)%v2))
    CALL nf(nf_put_var_double(ncid, varid21, gg%edges%dual_normal(:)%v1))
    CALL nf(nf_put_var_double(ncid, varid22, gg%edges%dual_normal(:)%v2) )
    CALL nf(nf_put_var_int   (ncid, varid23, gg%cells%edge_orientation))
    CALL nf(nf_put_var_int   (ncid, varid24, gg%cells%idx))
    CALL nf(nf_put_var_int   (ncid, varid25, gg%cells%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid251, gg%cells%parent_child_type))
    CALL nf(nf_put_var_int   (ncid, varid26, gg%cells%neighbor_index))
    CALL nf(nf_put_var_int   (ncid, varid27, gg%cells%child_index))
    CALL nf(nf_put_var_int   (ncid, varid41, gg%cells%child_id))
    CALL nf(nf_put_var_int   (ncid, varid28, gg%edges%idx))
    !      CALL nf(nf_put_var_int   (ncid, varid281, gg%edges%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid282, gg%edges%parent_child_type))
    CALL nf(nf_put_var_int   (ncid, varid29, gg%verts%idx))
    CALL nf(nf_put_var_int   (ncid, varid30, gg%verts%edge_orientation))
    CALL nf(nf_put_var_int   (ncid, varid31, gg%edges%system_orientation))

    ! Variables added for mesh refinement
    CALL nf(nf_put_var_int   (ncid, varid38, gg%edges%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid39, gg%edges%child_index))
    CALL nf(nf_put_var_int   (ncid, varid42, gg%edges%child_id))
    CALL nf(nf_put_var_int   (ncid, varid32, gg%cells%refin_ctrl))
    CALL nf(nf_put_var_int   (ncid, varid33, gg%cells%indlist))
    CALL nf(nf_put_var_int   (ncid, varid43, gg%cells%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid44, gg%cells%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid34, gg%edges%refin_ctrl))
    CALL nf(nf_put_var_int   (ncid, varid35, gg%edges%indlist))
    CALL nf(nf_put_var_int   (ncid, varid45, gg%edges%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid46, gg%edges%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid36, gg%verts%refin_ctrl))
    CALL nf(nf_put_var_int   (ncid, varid37, gg%verts%indlist))
    CALL nf(nf_put_var_int   (ncid, varid47, gg%verts%start_idx))
    CALL nf(nf_put_var_int   (ncid, varid48, gg%verts%end_idx))
    CALL nf(nf_put_var_int   (ncid, varid403,gg%verts%parent_index))
    CALL nf(nf_put_var_int   (ncid, varid49, gg%cells%phys_id))
    CALL nf(nf_put_var_int   (ncid, varid50, gg%edges%phys_id))
    !------------------------------------------------------------------------

    CALL nf(nf_close(ncid))

    !   ENDDO

  END SUBROUTINE write_grid
  !
  ! !SUBROUTINE INTERFACE:
  SUBROUTINE write_edges_grid(gg, outfile)
    !
    ! !DESCRIPTION:
    ! Creates files  with  cross reference tables and coordinates
    ! of grid items to be read by the model. In most cases for all grid levels,
    ! the file {\bf GRIDMAP.*} is produced, where * denotes the grid level.
    ! The exception is Heikes-Randall + C grid optimization in the last step. Here
    ! only the {\bf GRIDMAP} file of the last level is created. The background
    ! is in this case the {\bf GRIDMAP} files of all the coarser levels would
    ! produce just Heikes-Randall optimized grid without the C grid optimization.
    ! Only in this case gridgen has to be run for all grid levels.\\
    ! The files are produced as sequential access, unformatted output.
    !
    ! !REVISION HISTORY:
    ! Leonidas Lianrdakis, MPI-M (2004-5).
    !
    TYPE(t_grid), INTENT(in) :: gg
    CHARACTER(LEN=filename_max), INTENT(in) :: outfile

    INTEGER :: i_nc, i_ne, i_nv   ! number of cells, edges and vertices

    INTEGER :: old_mode

    INTEGER :: ncid
    INTEGER :: dimids(2)

    INTEGER :: dim_nvertex, dim_nedge,           &
      & dim_nvertex_per_cell, dim_ncells_per_edge, dim_nedges_per_vertex

    INTEGER :: varid_vlon, varid_vlat
    INTEGER :: varid10

#ifdef __SX__
    INTEGER :: iargc
    CHARACTER(LEN= 32) :: arg_str
#endif

    !  REAL(wp), ALLOCATABLE :: zv2d(:,:), zv2dx(:,:), zv2dy(:,:)


    !  REAL(wp) :: swap(4)

    INTEGER :: ilevel, grid_root

    !EOP
    !-------------------------------------------------------------------------

    ! distinguish between gridgeneration for optimization strategies

    ! Dummy settings for special grids
    grid_root = 0
    ilevel    = 1

    WRITE(message_text,'(a,a)') 'Write edges gridmap file: ', TRIM(outfile)
    CALL message ('', TRIM(message_text))


    i_nc = gg%ncells
    i_ne = gg%nedges
    i_nv = gg%nverts


    !----------------------------------------------------------------------
    !
    CALL nf(nf_set_default_format(nf_format_64bit, old_mode))

    CALL nf(nf_create(TRIM(outfile), nf_clobber, ncid))
    CALL nf(nf_set_fill(ncid, nf_nofill, old_mode))

    !
    !----------------------------------------------------------------------
    !
    ! Global attributes
    !
    !      CALL nf(nf_put_att_text    (ncid, NF_GLOBAL, 'title', 21, 'ICON grid description'))
    !      CALL nf(nf_put_att_text    (ncid, NF_GLOBAL, 'institution', 59, &
    !           'Max Planck Institute for Meteorology/Deutscher Wetterdienst'))
    !
    ! Dimensions
    !
    CALL nf(nf_def_dim(ncid, 'vertex', i_nv, dim_nvertex))
    CALL nf(nf_def_dim(ncid, 'edge',   i_ne, dim_nedge))
    CALL nf(nf_def_dim(ncid, 'nc',        2, dim_ncells_per_edge))
    !
    CALL nf(nf_def_dim(ncid, 'nv',  gg%cells%max_no_of_vertices, dim_nvertex_per_cell))
    CALL nf(nf_def_dim(ncid, 'ne',  gg%verts%max_connectivity, dim_nedges_per_vertex))
    !

    !
    !---------------------------------------------------------------------
    !
    ! Grid variables
    !
    !---------------------------------------------------------------------
    !
    ! public grid information
    !
    ! cell part:
    !
    !
    ! vertex part:
    !
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
    !
    ! edge part:
    !
    !
    dimids = (/ dim_nedge, dim_ncells_per_edge /)
    CALL nf(nf_def_var(ncid, 'edge_vertices', nf_int, 2, dimids, varid10))
    CALL nf(nf_put_att_text(ncid, varid10, 'long_name', 35, &
      & 'vertices at the end of of each edge'))
    CALL nf(nf_put_att_text(ncid, varid10, 'cdi', 6, 'ignore'))
    !
    !
    !
    CALL nf(nf_enddef(ncid))
    !
    !-------------------------------------------------------------------------
    !
    ! vertex part:
    !
    CALL nf(nf_put_var_double(ncid, varid_vlon,  gg%verts%vertex(:)%lon))
    CALL nf(nf_put_var_double(ncid, varid_vlat,  gg%verts%vertex(:)%lat))
    !------------------------------------------------------------------------
    ! edge part:
    CALL nf(nf_put_var_int   (ncid, varid10, gg%edges%vertex_index))

    !------------------------------------------------------------------------
    CALL nf(nf_close(ncid))

    !   ENDDO

  END SUBROUTINE write_edges_grid

END MODULE mo_io_grid



