PROGRAM grid_command

  USE mo_kind
  USE mo_io_units,  ONLY: filename_max
  USE mo_exception, ONLY: finish

  USE mo_create_torus_grid,     ONLY: create_torus_grid
  USE mo_create_ocean_grid,     ONLY: create_ocean_grid, remove_land_points
  USE mo_local_patch_hierarchy, ONLY: create_patches
  USE mo_local_grid_refinement, ONLY: grid_refine, coarsen_grid_file
  USE mo_grid_toolbox,          ONLY: concatenate_grid_files,create_dual, &
    &   shift_grid_ids, grid_file_zero_children
  USE mo_grid_checktools,       ONLY: check_grid_file, grid_statistics_file,&
    & check_inverse_connect_verts, check_read_write_grid, calculate_triangle_properties
  USE mo_local_grid_optimization, ONLY: optimize_grid_file
  USE mo_icosahedron_grid,      ONLY: create_icon_grid
  USE mo_grid_conditions,       ONLY: cut_local_grid, cut_local_grid_ascii, flag_conditional_grid
  USE mo_grid_decomposition,    ONLY: decompose_file_round_robin_opp, &
    & inherit_file_decomposition, print_decomposition_statistics,   &
    & decompose_file_metis, file_pair_opposite_subdomains,          &
    & shrink_file_decomposition, file_reorder_latlon_subdomains,    &
    & file_reorder_lonlat_subdomains, file_cluster_subdomains,      &
    & decompose_file_geometric_medial, decomp_file_geom_medial_cluster
  USE mo_local_grid_geometry,   ONLY:  compute_sphere_geometry, file_rotate_sphere_xaxis
  
#ifndef __ICON_GRID_GENERATOR__
  USE mo_global_grid_generator, ONLY: global_graph_generator, global_grid_generator, &
    & prepare_gridref
#endif
 
  USE mo_mpi,                   ONLY: start_mpi, p_stop

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: command_file = 'command.grid'

  CHARACTER(len=32), PARAMETER :: create_dual_c ='create_dual'
  CHARACTER(len=32), PARAMETER :: create_torus_c ='create_torus'
  CHARACTER(len=32), PARAMETER :: create_ocean_c ='create_ocean'
  CHARACTER(len=32), PARAMETER :: remove_land_points_c ='remove_land_points'
  CHARACTER(len=32), PARAMETER :: refine_grid_c ='refine_grid'
  CHARACTER(len=32), PARAMETER :: cut_local_grid_c ='cut_local_grid'
  CHARACTER(len=32), PARAMETER :: cut_local_grid_ascii_c ='cut_local_grid_ascii'
  CHARACTER(len=32), PARAMETER :: flag_conditional_grid_c ='flag_conditional_grid'
  CHARACTER(len=32), PARAMETER :: zero_chilrden_c ='zero_children'
  CHARACTER(len=32), PARAMETER :: create_patch_hierarchy_c ='create_patch_hierarchy'
  CHARACTER(len=32), PARAMETER :: concatenate_grids_c ='concatenate_grids'
  CHARACTER(len=32), PARAMETER :: check_grid_c ='check_grid'
  CHARACTER(len=32), PARAMETER :: check_read_write_c ='read_write_grid'
  CHARACTER(len=32), PARAMETER :: inv_vertex_connect_c ='inv_vertex_connectivity'
  CHARACTER(len=32), PARAMETER :: compute_sphere_geometry_c ='compute_sphere_geometry'
  CHARACTER(len=32), PARAMETER :: rotate_sphere_xaxis_c ='rotate_sphere_xaxis'
  CHARACTER(len=32), PARAMETER :: grid_statistics_c ='grid_statistics'
  CHARACTER(len=32), PARAMETER :: optimize_grid_c ='optimize_grid'
  CHARACTER(len=32), PARAMETER :: create_icon_grid_c ='create_icon_grid'
  CHARACTER(len=32), PARAMETER :: shift_grid_ids_c ='shift_grid_ids'
  CHARACTER(len=32), PARAMETER :: coarsen_grid_c ='coarsen_grid'
  CHARACTER(len=32), PARAMETER :: decompose_round_robin_opp_c ='decompose_round_robin_opp'
  CHARACTER(len=32), PARAMETER :: decompose_metis_c ='decompose_metis'
  CHARACTER(len=32), PARAMETER :: decompose_geometric_medial_c ='decompose_geometric_medial'
  CHARACTER(len=32), PARAMETER :: decomp_geom_medial_cluster_c = &
    & 'decompose_geometric_medial_cluster'
  CHARACTER(len=32), PARAMETER :: inherit_decomposition_c ='inherit_decomposition'
  CHARACTER(len=32), PARAMETER :: redecompose_pair_opposites_c ='redecompose_pair_opposites'
  CHARACTER(len=32), PARAMETER :: shrink_decomposition_c ='shrink_decomposition'
  CHARACTER(len=32), PARAMETER :: reorder_latlon_subdomains_c ='reorder_latlon_subdomains'
  CHARACTER(len=32), PARAMETER :: reorder_lonlat_subdomains_c ='reorder_lonlat_subdomains'
  CHARACTER(len=32), PARAMETER :: cluster_subdomains_c ='cluster_subdomains'
  CHARACTER(len=32), PARAMETER :: decomposition_statistics_c = &
    & 'decomposition_statistics'
  
  ! drives the old grid generator
  CHARACTER(len=32), PARAMETER :: global_graph_generator_c ='global_graph_generator'
  CHARACTER(len=32), PARAMETER :: global_grid_generator_c ='global_grid_generator'
  CHARACTER(len=32), PARAMETER :: global_grid_refine_c ='global_grid_refine'
  CHARACTER(len=32), PARAMETER :: calculate_triangle_c ='calculate_triangle'

  CHARACTER(len=32) :: command =''
  CHARACTER(len=filename_max) :: param_1 = ''
  CHARACTER(len=filename_max) :: param_2 = ''
  CHARACTER(len=filename_max) :: param_3 = ''

  INTEGER :: int_param_1, int_param_2, int_param_3
  REAL(wp):: real_wp_param_1, real_wp_param_2, real_wp_param_3, real_wp_param_4, &
    & real_wp_param_5, real_wp_param_6
  REAL :: real_param_1
  LOGICAL :: logical_flag
!   CALL get_command_argument(1, command)
!   CALL get_command_argument(2, param_1)

  CALL start_mpi()

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
    
    CASE (cut_local_grid_c)
      CALL cut_local_grid(param_1)
    
    CASE (cut_local_grid_ascii_c)
      CALL cut_local_grid_ascii(param_1)

    CASE (flag_conditional_grid_c)
      CALL flag_conditional_grid(param_1)
    
    CASE (remove_land_points_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL remove_land_points(in_file_name=param_1, out_file_name=param_2)
    
    CASE (concatenate_grids_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2, param_3
      CLOSE(500)
      CALL concatenate_grid_files(param_1, param_2, param_3)

    CASE (coarsen_grid_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2, param_3
      CLOSE(500)
      CALL coarsen_grid_file(param_1, param_2, param_3)

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

    CASE (zero_chilrden_c)
!       CALL get_command_argument(3, out_file)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL grid_file_zero_children(param_1, param_2)

    CASE (shift_grid_ids_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1
      CLOSE(500)
      CALL shift_grid_ids(param_1, int_param_1)

    CASE (inv_vertex_connect_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL check_inverse_connect_verts(param_1, param_2)

    CASE (compute_sphere_geometry_c)
      CALL compute_sphere_geometry(param_file_name=param_1)

    CASE (rotate_sphere_xaxis_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, real_wp_param_1, param_2
      CLOSE(500)
      CALL file_rotate_sphere_xaxis(param_1, real_wp_param_1, param_2)

    CASE (create_torus_c)
      CALL create_torus_grid(param_1)

    CASE (optimize_grid_c)
      CALL optimize_grid_file(param_1)

    CASE (check_grid_c)
      CALL check_grid_file(param_1)

    CASE (decompose_metis_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, param_2, int_param_2, &
        & int_param_3, real_param_1, param_3
      CLOSE(500)
      logical_flag = (param_3 == "use_sea_land_mask")
      CALL decompose_file_metis(grid_file=param_1, dec_size=int_param_1, out_ascii_file=param_2,&
        & decomposition_id=int_param_2, edge_weight=int_param_3, allowed_imbalance=real_param_1, &
        & use_sea_land_mask=logical_flag )
      ! decompose_file_metis(grid_file, dec_size, out_ascii_file, decomposition_id )
    
    CASE (decompose_geometric_medial_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2, param_2
      CLOSE(500)
      CALL decompose_file_geometric_medial(grid_file=param_1, decomposition_size=int_param_1, &
        & out_decomposition_id=int_param_2, out_ascii_file=param_2)
    
    CASE (decomp_geom_medial_cluster_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2,  int_param_3, param_2
      CLOSE(500)
      CALL decomp_file_geom_medial_cluster(grid_file=param_1, decomposition_size=int_param_1, &
        & radiation_split_factor=int_param_2, out_decomposition_id=int_param_3, &
        & out_ascii_file=param_2)
    
    CASE (decompose_round_robin_opp_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2, int_param_3, param_2
      CLOSE(500)
      CALL decompose_file_round_robin_opp(grid_file=param_1,                         &
        &  in_decomposition_id=int_param_1, out_decomposition_id=int_param_2,        &
        &  subdomain_partition=int_param_3, out_ascii_file=param_2  )
    
    
    CASE (inherit_decomposition_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, param_2, int_param_2, param_3
      CLOSE(500)
      CALL inherit_file_decomposition(parent_file=param_1, parent_decomposition_id=int_param_1, &
       & child_file=param_2, child_decomposition_id=int_param_2, out_ascii_file=param_3)

    CASE (redecompose_pair_opposites_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2, param_2
      CLOSE(500)
      CALL file_pair_opposite_subdomains(param_1, int_param_1, int_param_2,&
        & out_ascii_file=param_2)

    CASE (reorder_latlon_subdomains_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2,  param_2
      CLOSE(500)
      CALL file_reorder_latlon_subdomains(param_1, int_param_1, int_param_2, &
        & out_ascii_file=param_2)

    CASE (reorder_lonlat_subdomains_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2,  param_2
      CLOSE(500)
      CALL file_reorder_lonlat_subdomains(param_1, int_param_1, int_param_2, &
        & out_ascii_file=param_2)

    CASE (cluster_subdomains_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2,  param_2
      CLOSE(500)
      CALL file_cluster_subdomains(grid_file=param_1, in_decomposition_id=int_param_1, &
        & out_decomposition_id=int_param_2, &
        & out_ascii_file=param_2)
    
    CASE (shrink_decomposition_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1, int_param_2, param_2
      CLOSE(500)
      CALL shrink_file_decomposition(param_1, int_param_1, int_param_2, out_ascii_file=param_2)
    
    CASE (decomposition_statistics_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, int_param_1
      CLOSE(500)
      CALL print_decomposition_statistics(param_1, int_param_1)

    CASE (grid_statistics_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, param_1, param_2
      CLOSE(500)
      CALL grid_statistics_file(param_1, param_2)
    
    CASE (calculate_triangle_c)
      OPEN (500, FILE = command_file,STATUS = 'OLD')
      READ (500, *) command, real_wp_param_1, real_wp_param_2, &
        & real_wp_param_3, real_wp_param_4, real_wp_param_5, real_wp_param_6
      CLOSE(500)
      CALL calculate_triangle_properties(lon1=real_wp_param_1, lat1=real_wp_param_2,  &
        & lon2=real_wp_param_3, lat2=real_wp_param_4, lon3=real_wp_param_5, lat3=real_wp_param_6 )
      ! decompose_file_metis(grid_file, dec_size, out_ascii_file, decomposition_id )

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
