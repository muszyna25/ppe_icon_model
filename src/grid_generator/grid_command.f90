PROGRAM grid_command

  USE mo_kind
  USE mo_io_units,  ONLY: filename_max
  USE mo_exception, ONLY: finish

  USE mo_create_torus_grid,     ONLY: create_torus_grid
  USE mo_create_ocean_grid,     ONLY: create_ocean_grid
  USE mo_local_patch_hierarchy, ONLY: create_patches
  USE mo_local_grid_refinement, ONLY: grid_refine
  USE mo_grid_toolbox,          ONLY: concatenate_grid_files,create_dual, &
    &   shift_grid_ids
  USE mo_grid_checktools,       ONLY: check_grid_file, grid_statistics_file,&
    & check_inverse_connect_verts, check_compute_sphere_geometry, &
    & check_read_write_grid
  USE mo_local_grid_optimization, ONLY: optimize_grid_file
  USE mo_icosahedron_grid,      ONLY: create_icon_grid
#ifndef __ICON_GRID_GENERATOR__
  USE mo_global_grid_generator, ONLY: global_graph_generator, global_grid_generator, &
    & prepare_gridref
#endif
 
  USE mo_mpi,                   ONLY: p_start, p_stop

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: command_file = 'command.grid'

  CHARACTER(len=32), PARAMETER :: create_dual_c ='create_dual'
  CHARACTER(len=32), PARAMETER :: create_torus_c ='create_torus'
  CHARACTER(len=32), PARAMETER :: create_ocean_c ='create_ocean'
  CHARACTER(len=32), PARAMETER :: refine_grid_c ='refine_grid'
  CHARACTER(len=32), PARAMETER :: create_patch_hierarchy_c ='create_patch_hierarchy'
  CHARACTER(len=32), PARAMETER :: concatenate_grids_c ='concatenate_grids'
  CHARACTER(len=32), PARAMETER :: check_grid_c ='check_grid'
  CHARACTER(len=32), PARAMETER :: check_read_write_c ='read_write_grid'
  CHARACTER(len=32), PARAMETER :: inv_vertex_connect_c ='inv_vertex_connectivity'
  CHARACTER(len=32), PARAMETER :: compute_sphere_geometry_c ='compute_sphere_geometry'
  CHARACTER(len=32), PARAMETER :: grid_statistics_c ='grid_statistics'
  CHARACTER(len=32), PARAMETER :: optimize_grid_c ='optimize_grid'
  CHARACTER(len=32), PARAMETER :: create_icon_grid_c ='create_icon_grid'
  CHARACTER(len=32), PARAMETER :: global_graph_generator_c ='global_graph_generator'
  CHARACTER(len=32), PARAMETER :: global_grid_generator_c ='global_grid_generator'
  CHARACTER(len=32), PARAMETER :: global_grid_refine_c ='global_grid_refine'
  CHARACTER(len=32), PARAMETER :: shift_grid_ids_c ='shift_grid_ids'

  CHARACTER(len=32) :: command =''
  CHARACTER(len=filename_max) :: param_1 = ''
  CHARACTER(len=filename_max) :: param_2 = ''
  CHARACTER(len=filename_max) :: param_3 = ''

  INTEGER :: int_1
!   CALL get_command_argument(1, command)
!   CALL get_command_argument(2, param_1)

  CALL p_start()

  OPEN (500, FILE = command_file,STATUS = 'OLD')
  READ (500, *) command, param_1
  CLOSE(500)
  print *, TRIM(command), ', ', TRIM(param_1)

  SELECT CASE (command)

    CASE (create_ocean_c)
      CALL create_ocean_grid(param_1)

    CASE (refine_grid_c)
      CALL grid_refine(param_1)

    CASE (create_patch_hierarchy_c)
      CALL create_patches(param_1)

    CASE (create_icon_grid_c)
      CALL create_icon_grid(param_1)

    CASE (concatenate_grids_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2, param_3
      CLOSE(500)
!       print *, 'command:', command
!       print *, 'param_1:', param_1
!       print *, 'param_2:', param_2
!       print *, 'param_3:', param_3
!       CALL get_command_argument(3, in_file_2)
!       CALL get_command_argument(4, out_file)
      CALL concatenate_grid_files(param_1, param_2, param_3)

    CASE (create_dual_c)
!       CALL get_command_argument(3, out_file)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL create_dual(param_1, param_2)

    CASE (check_read_write_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL check_read_write_grid(param_1, param_2)

    CASE (shift_grid_ids_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_1
      CLOSE(500)
      CALL shift_grid_ids(param_1, int_1)

    CASE (inv_vertex_connect_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL check_inverse_connect_verts(param_1, param_2)

    CASE (compute_sphere_geometry_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL check_compute_sphere_geometry(param_1, param_2)

    CASE (create_torus_c)
      CALL create_torus_grid(param_1)

    CASE (optimize_grid_c)
      CALL optimize_grid_file(param_1)

    CASE (check_grid_c)
      CALL check_grid_file(param_1)

    CASE (grid_statistics_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL grid_statistics_file(param_1, param_2)

#ifdef __ICON__
    CASE (global_graph_generator_c)
      CALL global_graph_generator()

    CASE (global_grid_generator_c)
      CALL global_grid_generator()
      
    CASE (global_grid_refine_c)
      CALL prepare_gridref()
#endif

    CASE DEFAULT
      CALL finish('grid_command', 'Unkown command')

  END SELECT
  !----------------------------------------------------------
 
  CALL p_stop()

END PROGRAM grid_command
