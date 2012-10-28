!-------------------------------------------------------------------------------------
!>
!! Set of methods for grid decomposition
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 20011-12-6
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
MODULE mo_grid_decomposition
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, finish, warning
  USE mo_io_units,           ONLY: nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_timer,              ONLY: new_timer, timer_start, timer_stop, print_timer, delete_timer
  USE mo_base_geometry
  USE mo_local_grid
  USE mo_grid_toolbox,       ONLY: add_to_list_if_not_exist
  USE mo_io_local_grid,      ONLY: read_new_netcdf_grid, write_netcdf_grid, &
    & write_ascii_decomposition
  USE mo_statistics_tools
  USE mo_decomposition_tools

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: decompose_all_cells, decompose_cluster_all_cells
  PUBLIC :: decompose_file_round_robin_opp
  PUBLIC :: decompose_file_metis, decompose_metis
  PUBLIC :: inherit_file_decomposition, inherit_decomposition
  PUBLIC :: file_pair_opposite_subdomains, pair_opposite_subdomains
  PUBLIC :: print_decomposition_statistics
  PUBLIC :: get_max_domain_id, get_no_of_domains
  PUBLIC :: shrink_file_decomposition, shrink_decomposition
  PUBLIC :: file_reorder_latlon_subdomains, file_reorder_lonlat_subdomains
  PUBLIC :: file_cluster_subdomains, grid_cluster_subdomains
  PUBLIC :: get_max_subdomain_cells
  PUBLIC :: decompose_file_geometric_medial, decomp_file_geom_medial_cluster
  !-------------------------------------------------------------------------

  INTEGER, PARAMETER   :: max_halos = 4096! assume we have max of 32 neighbors, for statistics
  

CONTAINS


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE decompose_file_geometric_medial(grid_file, decomposition_size, &
    & out_decomposition_id, out_ascii_file)
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: decomposition_size, out_decomposition_id

    INTEGER :: grid_id

    TYPE(t_decomposition_structure) :: decomposition_struct

    grid_id = read_new_netcdf_grid(grid_file)
    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)
    
    CALL divide_geometric_medial(decomposition_struct, decomposition_size, &
    &  out_decomposition_id, cluster=.false.)
             
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    
  END SUBROUTINE decompose_file_geometric_medial
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE decomp_file_geom_medial_cluster(grid_file, decomposition_size, &
    & out_decomposition_id, out_ascii_file)
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: decomposition_size, out_decomposition_id

    INTEGER :: grid_id

    TYPE(t_decomposition_structure) :: decomposition_struct

    grid_id = read_new_netcdf_grid(grid_file)
    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)
    
    CALL divide_geometric_medial(decomposition_struct, &
      & decomposition_size = decomposition_size,       &
      & out_decomposition_id = out_decomposition_id,   &
      & cluster=.false.)
             
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    
  END SUBROUTINE decomp_file_geom_medial_cluster
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE decompose_file_round_robin_opp(grid_file,  &
    & in_decomposition_id, out_decomposition_id, subdomain_partition, out_ascii_file)
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: subdomain_partition, in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    TYPE(t_decomposition_structure) :: decomposition_struct

    grid_id = read_new_netcdf_grid(grid_file)
    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)
    
    CALL decompose_round_robin_opp(decomposition_struct,  &
      &  in_decomposition_id  = in_decomposition_id,         &
      &  out_decomposition_id = out_decomposition_id,        &
      &  subdomain_partition  = subdomain_partition)
             
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    

  END SUBROUTINE decompose_file_round_robin_opp
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE grid_round_robin_decomp(grid_id, subdomain_partition, &
    & in_decomposition_id, out_decomposition_id, opposite_subdomain_id)
    INTEGER, INTENT(in)  :: grid_id, subdomain_partition, &
      & in_decomposition_id, out_decomposition_id
    INTEGER, POINTER, OPTIONAL  :: opposite_subdomain_id(:)
    
    TYPE(t_decomposition_structure) :: decomposition_struct

    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)

    CALL decompose_round_robin(decomposition_struct,   &
      & in_decomposition_id = in_decomposition_id,     &
      & out_decomposition_id = out_decomposition_id,   &
      & subdomain_partition = subdomain_partition)
    
  END SUBROUTINE grid_round_robin_decomp
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! decompose a grid file using metis
  SUBROUTINE decompose_file_metis(grid_file, dec_size, out_ascii_file, &
    & decomposition_id, edge_weight, allowed_imbalance, use_sea_land_mask )
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: dec_size, decomposition_id
    LOGICAL, INTENT(in)  :: use_sea_land_mask
    INTEGER, INTENT(in)  :: edge_weight
    REAL, INTENT(in)     :: allowed_imbalance

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    CALL decompose_metis(grid_id, dec_size, decomposition_id, edge_weight, &
      & allowed_imbalance, use_sea_land_mask)
    CALL write_ascii_decomposition(grid_id, decomposition_id, out_ascii_file)
    CALL write_netcdf_grid(grid_id)
    
    CALL delete_grid(grid_id)    

  END SUBROUTINE decompose_file_metis
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! decompose a grid using metis
  SUBROUTINE decompose_metis(grid_id, dec_size, decomposition_id, edge_weight,&
    & allowed_imbalance, use_sea_land_mask)
    !-------------------------------------------------------------------------
#ifdef __HAS_METIS__    
    use iso_c_binding
    INTERFACE F_metis_partgraphrecursive
      SUBROUTINE C_metis_partgraphrecursive(nvtxs, ncon, xadj,&
                    adjncy, vwgt, vsize, adjwgt,              &
                    nparts, tpwgts, ubvec, options,           &
                    edgecut, part) BIND(C, NAME='METIS_PARTGRAPHRECURSIVE')
      use iso_c_binding
      INTEGER(C_INT)  :: nvtxs, nparts, ncon, edgecut
      INTEGER(C_INT), DIMENSION(*)  :: xadj, adjncy, vwgt, adjwgt, &
                    options, part
!      TYPE(C_PTR) :: vsize, tpwgts, ubvec
      INTEGER(C_INT), DIMENSION(*) :: vsize
      REAL(C_FLOAT), DIMENSION(*)  :: tpwgts, ubvec
      END SUBROUTINE C_metis_partgraphrecursive
    END INTERFACE F_metis_partgraphrecursive
    !-------------------------------------------------------------------------
    
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    INTEGER :: no_of_cells, no_of_edges
    INTEGER :: cell_no, neighbor, connect, adjacency_index, adjacency_size
    
    INTEGER, POINTER :: metis_xadj(:), metis_adjncy(:), metis_color(:)
    INTEGER, POINTER :: metis_cell_weights(:), metis_edge_weights(:), metis_options(:)
    INTEGER :: metis_constrains, metis_edgecut
    INTEGER, POINTER :: null_int(:)
    REAL, POINTER :: null_real(:), imbalance_conditions(:)
#endif

    INTEGER, INTENT(in)  :: grid_id, dec_size, decomposition_id
    LOGICAL, INTENT(in)  :: use_sea_land_mask
    INTEGER, INTENT(in)  :: edge_weight
    REAL, INTENT(in)     :: allowed_imbalance
    CHARACTER(*), PARAMETER :: method_name = "decompose_metis"

#ifndef __HAS_METIS__
   CALL finish(method_name, "no metis library is defined")
   RETURN
#else 
    

    cells => get_cells(grid_id)
    edges => get_edges(grid_id)
    no_of_cells=cells%no_of_existcells
    no_of_edges=edges%no_of_existedges

    ! fill metis parameters
    adjacency_size = no_of_edges * 2
    ALLOCATE(metis_xadj(0:no_of_cells))
    ALLOCATE(metis_adjncy(0:adjacency_size))
    ALLOCATE(metis_color(0:no_of_cells - 1))
    ALLOCATE(metis_cell_weights(0:no_of_cells - 1))
    ALLOCATE(metis_edge_weights(0:adjacency_size))
    ALLOCATE(metis_options(0:40))
    ALLOCATE(imbalance_conditions(0:1))
    NULLIFY(null_int)
    NULLIFY(null_real)
    ! ALLOCATE(metis_constrains(0:1))
    
    metis_constrains = 1
    metis_options(:)    =-1

    
    !fill adjacency
    adjacency_index=0
    ! fill edges weights
    DO cell_no = 1, cells%no_of_existcells
      metis_xadj(cell_no - 1) = adjacency_index
      
      ! fill cells weights
      metis_cell_weights(cell_no - 1) = 1
      IF (use_sea_land_mask) THEN
        IF (cells%sea_land_mask(cell_no) > 0) &
          metis_cell_weights(cell_no - 1) = 0
      ENDIF
      
      DO connect=1, cells%max_no_of_vertices
        neighbor = cells%get_neighbor_index(cell_no, connect)
        IF (neighbor <= 0) CYCLE
        metis_adjncy(adjacency_index) = neighbor - 1
        
        ! fill edgesweights
        metis_edge_weights(adjacency_index) = edge_weight
        ! this will be used once the communication has been adapted
!         IF (use_sea_land_mask) THEN
!           IF (cells%sea_land_mask(cell_no) > 0 .OR. &
!               cells%sea_land_mask(neighbor) >0  ) &
!            metis_edge_weights(adjacency_index) = 1
!         ENDIF
        
        adjacency_index = adjacency_index + 1
        IF (adjacency_index > adjacency_size) &
          CALL finish(method_name, "adjacency_index > adjacency_size")
      ENDDO !connect=1, max_no_of_vertices
      
    ENDDO
    metis_xadj(cell_no - 1) = adjacency_index
    imbalance_conditions(0) = allowed_imbalance
    imbalance_conditions(1) = allowed_imbalance
    
    ! call metis
    CALL F_metis_partgraphrecursive(no_of_cells, metis_constrains, &
      & metis_xadj, metis_adjncy,                                  &
      & metis_cell_weights, null_int, metis_edge_weights,          &
      & dec_size, null_real, imbalance_conditions, metis_options, metis_edgecut, metis_color)
! C_NULL_PTR
    cells%no_of_domains(decomposition_id) = dec_size

    !fill the domain ids
!$OMP PARALLEL
!$OMP DO PRIVATE(cell_no)
    DO cell_no = 1, cells%no_of_existcells
      cells%domain_id(decomposition_id, cell_no) = metis_color(cell_no - 1)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    DEALLOCATE(metis_xadj,metis_adjncy,metis_color,   &
      & metis_cell_weights,metis_edge_weights,        &
      & imbalance_conditions, metis_options)
#endif

  END SUBROUTINE decompose_metis
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE inherit_file_decomposition(parent_file, parent_decomposition_id, &
    & child_file, child_decomposition_id, out_ascii_file)
    
    CHARACTER(LEN=*), INTENT(in) :: parent_file, child_file, out_ascii_file
    INTEGER, INTENT(in)  :: parent_decomposition_id, child_decomposition_id

    INTEGER :: parent_grid_id, child_grid_id

    parent_grid_id = read_new_netcdf_grid(parent_file)
    child_grid_id = read_new_netcdf_grid(child_file)
    
    CALL inherit_decomposition(parent_grid_id, parent_decomposition_id, &
      & child_grid_id, child_decomposition_id)
    
    CALL write_netcdf_grid(child_grid_id)
    CALL write_ascii_decomposition(child_grid_id, child_decomposition_id, out_ascii_file)
    
    CALL delete_grid(parent_grid_id)
    CALL delete_grid(child_grid_id)    

  END SUBROUTINE inherit_file_decomposition
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE inherit_decomposition(parent_grid_id, parent_decomposition_id, &
      & child_grid_id, child_decomposition_id)

   INTEGER, INTENT(in) :: parent_grid_id, parent_decomposition_id, &
      & child_grid_id, child_decomposition_id
    
    TYPE(t_grid_cells), POINTER :: parent_cells, child_cells
    INTEGER, POINTER :: cells_per_domain(:), next_new_domain(:)
    
    INTEGER :: child_cell_no, parent_cell_no, domain_id
    INTEGER :: no_of_domains
    
    CHARACTER(*), PARAMETER :: method_name = "inherit_decomposition"
         
    parent_cells => get_cells(parent_grid_id)
    child_cells => get_cells(child_grid_id)
    no_of_domains = parent_cells%no_of_domains(parent_decomposition_id)
    IF (no_of_domains <= 0) &
      & CALL finish (method_name, "not valid parent decomposition")
    
    DO child_cell_no = 1, child_cells%no_of_existcells
      parent_cell_no = child_cells%parent_index(child_cell_no)
      
      IF (parent_cell_no <= 0 .OR. parent_cell_no > parent_cells%no_of_existcells) &
        & CALL finish (method_name, "parent cell number not valid")
      domain_id = parent_cells%domain_id(parent_decomposition_id, parent_cell_no)
      IF (domain_id < 0 .OR. domain_id >= no_of_domains) &
        & CALL finish (method_name, "domain_id not valid")
      
      child_cells%domain_id(child_decomposition_id, child_cell_no) = &
        & domain_id
    ENDDO
    
    child_cells%no_of_domains(child_decomposition_id) = &
      & parent_cells%no_of_domains(parent_decomposition_id)

  END SUBROUTINE inherit_decomposition
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE shrink_file_decomposition(grid_file, in_decomposition_id, &
    & out_decomposition_id, out_ascii_file)
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    
    CALL shrink_decomposition(grid_id, in_decomposition_id, out_decomposition_id)
    
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)

  END SUBROUTINE shrink_file_decomposition
  !-------------------------------------------------------------------------

  !>
  ! Combines two consequetive subdomains into one
  SUBROUTINE shrink_decomposition(grid_id, in_decomposition_id, out_decomposition_id)

   INTEGER, INTENT(in) :: grid_id, in_decomposition_id, &
      & out_decomposition_id
    
    TYPE(t_grid_cells), POINTER :: cells
    INTEGER, POINTER :: cells_per_domain(:), next_new_domain(:)
    
    INTEGER, POINTER :: new_subdomain_id(:)
    INTEGER :: no_of_domains,new_no_of_domains, next_domain_id 
    INTEGER :: return_status 
    
    CHARACTER(*), PARAMETER :: method_name = "shrink_decomposition"
         
    cells => get_cells(grid_id)
    no_of_domains = cells%no_of_domains(in_decomposition_id)
    IF (no_of_domains <= 0) &
      & CALL finish (method_name, "not valid input decomposition")
    
    ALLOCATE(new_subdomain_id(0:no_of_domains-1), stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(new_subdomain_id")
    new_subdomain_id(:) = -1

    ! fill new_subdomain_id
    ! combine consequetive pairs of subdomains
    new_no_of_domains = ((no_of_domains+1)/2)-1
    DO next_domain_id = 0, ((no_of_domains+1)/2)-1
      new_subdomain_id(2*next_domain_id) = next_domain_id
      new_subdomain_id((2*next_domain_id)+1) = next_domain_id
      WRITE(0,*) 2*next_domain_id, (2*next_domain_id)+1, ' -> ', next_domain_id
    ENDDO
    
    CALL fill_redecomposition(grid_id, in_decomposition_id, out_decomposition_id, &
     & new_subdomain_id, new_no_of_domains)

    DEALLOCATE(new_subdomain_id)

  END SUBROUTINE shrink_decomposition
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE file_pair_opposite_subdomains(grid_file, &
    & in_decomposition_id, out_decomposition_id, out_ascii_file )
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    
    CALL grid_pair_opposite_subdomains(grid_id,  &
      & in_decomposition_id, out_decomposition_id)
    
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    
  END SUBROUTINE file_pair_opposite_subdomains
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it reorders the decomposition in pairs of
  !! opposite geographically domains
  SUBROUTINE grid_pair_opposite_subdomains(grid_id,  &
    & in_decomposition_id, out_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, &
      & in_decomposition_id, out_decomposition_id

    TYPE(t_decomposition_structure) :: decomposition_struct

    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)

    CALL pair_opposite_subdomains(decomposition_struct, in_decomposition_id, &
      & out_decomposition_id)
    
  END SUBROUTINE grid_pair_opposite_subdomains
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! count the number of cells in each decomposition
  SUBROUTINE grid_no_of_cells_per_subdomain(grid_id, in_decomposition_id, cells_per_domain)
    INTEGER, INTENT(in)  :: grid_id, in_decomposition_id
    INTEGER, POINTER :: cells_per_domain(:)

    TYPE(t_decomposition_structure) :: decomposition_struct

    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)
    CALL get_no_of_cells_per_subdomain(decomposition_struct, in_decomposition_id, &
      & cells_per_domain)
        
  END SUBROUTINE grid_no_of_cells_per_subdomain
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE file_cluster_subdomains(grid_file, &
    & in_decomposition_id, out_decomposition_id, out_ascii_file )
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    
    CALL grid_cluster_subdomains(grid_id,  &
    & in_decomposition_id, out_decomposition_id)
    
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    
  END SUBROUTINE file_cluster_subdomains
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Cluster the subdomains using the geometrically nearest neigbor (on the sphere)
  !! cluster_size is not used
  SUBROUTINE grid_cluster_subdomains(grid_id,  &
    & in_decomposition_id, out_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, &
      & in_decomposition_id, out_decomposition_id

    TYPE(t_decomposition_structure) :: decomposition_struct

    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)

    CALL cluster_subdomains(decomposition_struct, in_decomposition_id, &
      & out_decomposition_id)
        
  END SUBROUTINE grid_cluster_subdomains
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE file_reorder_latlon_subdomains(grid_file, &
    & in_decomposition_id, out_decomposition_id,  out_ascii_file )
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    
    CALL grid_reorder_latlon_subdomains(grid_id,  &
      & in_decomposition_id, out_decomposition_id)
    
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    
  END SUBROUTINE file_reorder_latlon_subdomains
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE file_reorder_lonlat_subdomains(grid_file, &
    & in_decomposition_id, out_decomposition_id,  out_ascii_file )
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    
    CALL grid_reorder_lonlat_subdomains(grid_id,  &
      & in_decomposition_id, out_decomposition_id)
    
    CALL write_netcdf_grid(grid_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    
  END SUBROUTINE file_reorder_lonlat_subdomains
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE grid_reorder_latlon_subdomains(grid_id,  &
    & in_decomposition_id, out_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, &
      & in_decomposition_id, out_decomposition_id
    
    TYPE(t_decomposition_structure) :: decomposition_struct

    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)

    CALL reorder_latlon_subdomains(decomposition_struct, in_decomposition_id, &
      & out_decomposition_id)
  
  END SUBROUTINE grid_reorder_latlon_subdomains
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE grid_reorder_lonlat_subdomains(grid_id,  &
    & in_decomposition_id, out_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, &
      & in_decomposition_id, out_decomposition_id
    
    TYPE(t_decomposition_structure) :: decomposition_struct

    CALL grid_fill_decomposition_struct(grid_id, decomposition_struct)

    CALL reorder_lonlat_subdomains(decomposition_struct, in_decomposition_id, &
      & out_decomposition_id)
  
  END SUBROUTINE grid_reorder_lonlat_subdomains
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE bubble_sort(in_values, start_index, end_index, sorted_list)
    REAL(wp), POINTER :: in_values(:)
    INTEGER, INTENT(in)  :: start_index, end_index
    INTEGER, POINTER :: sorted_list(:)

    LOGICAL :: keep_ordering 
    INTEGER :: i, exch, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "bubble_sort"

    ALLOCATE( sorted_list(start_index:end_index), stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(sorted_list")
    DO i=start_index, end_index
      sorted_list(i) = i
    ENDDO
                        
    keep_ordering = .true.
    DO WHILE(keep_ordering)
      keep_ordering = .false.
      DO i=start_index, end_index-1
        IF (in_values(sorted_list(i)) > in_values(sorted_list(i+1))) THEN
          exch = sorted_list(i)
          sorted_list(i) = sorted_list(i+1)
          sorted_list(i+1) = exch
          keep_ordering = .true.
        ENDIF
      ENDDO
    ENDDO
        
  END SUBROUTINE bubble_sort
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  ! Under construction
  SUBROUTINE advancing_front_decomposition(grid_id,  &
    & subdomain_size, min_subdomain_size, max_subdomain_size, out_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, subdomain_size, min_subdomain_size, &
      & max_subdomain_size, out_decomposition_id

    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_cartesian_coordinates), POINTER :: subdomain_center(:)
    TYPE(t_cartesian_coordinates) :: v, v2
    INTEGER, POINTER :: opposite_subdomain_id(:), new_subdomain_id(:)

    REAL(wp) :: min_distance, distance
    INTEGER :: cell_no, no_of_domains,in_domain_id, next_domain_id, min_subdomain_id
    INTEGER :: new_no_of_domains, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "advancing_front_decomposition"

    cells => get_cells(grid_id)
    write(0,*) method_name, ":out_decomposition_id=", &
      & out_decomposition_id

    
  END SUBROUTINE advancing_front_decomposition
  !-------------------------------------------------------------------------
  


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE fill_redecomposition(grid_id, in_decomposition_id, out_decomposition_id, &
    & new_subdomain_id, new_no_of_domains)

    INTEGER, INTENT(in)  :: grid_id, &
      & in_decomposition_id, out_decomposition_id, new_no_of_domains
    INTEGER, POINTER ::  new_subdomain_id(:)

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: cell_no, in_domain_id
    CHARACTER(*), PARAMETER :: method_name = "fill_redecomposition"

    cells => get_cells(grid_id)
    DO cell_no = 1, cells%no_of_existcells
      in_domain_id = cells%domain_id(in_decomposition_id, cell_no)
      IF (new_subdomain_id(in_domain_id) < 0 .OR. &
        & new_subdomain_id(in_domain_id) > new_no_of_domains) &
        & CALL finish (method_name, "non-valid new_subdomain_id")
        
      cells%domain_id(out_decomposition_id, cell_no) = &
        & new_subdomain_id(in_domain_id)
    ENDDO
    cells%no_of_domains(out_decomposition_id) = new_no_of_domains

   END SUBROUTINE fill_redecomposition    
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Each cell forms a different domain
  !! plus clustering
  SUBROUTINE decompose_cluster_all_cells(grid_id, decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id
    !-------------------------------------------------------------------------

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: cell_no

    CALL decompose_all_cells(grid_id, decomposition_id)
    CALL grid_cluster_subdomains(grid_id,  &
    & decomposition_id, decomposition_id)    

  END SUBROUTINE decompose_cluster_all_cells
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Each cell forms a different domain
  SUBROUTINE decompose_all_cells(grid_id, decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id
    !-------------------------------------------------------------------------

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: cell_no

    cells => get_cells(grid_id)
    cells%no_of_domains(decomposition_id) = cells%no_of_existcells
!$OMP PARALLEL
!$OMP DO PRIVATE(cell_no)
    DO cell_no = 1, cells%no_of_existcells
      cells%domain_id(decomposition_id, cell_no) = cell_no - 1
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE decompose_all_cells
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  !! returns the maximum number of subdomain cells 
  INTEGER FUNCTION get_max_subdomain_cells(grid_id, in_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, in_decomposition_id

    INTEGER, POINTER :: cells_per_domain(:)
    
    
    ! first count the number of cells in each decomposition
    NULLIFY(cells_per_domain)
    CALL grid_no_of_cells_per_subdomain(grid_id, in_decomposition_id, cells_per_domain)
    get_max_subdomain_cells = MAXVAL(cells_per_domain)
    DEALLOCATE(cells_per_domain)
  
  END FUNCTION get_max_subdomain_cells
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  !! print decompositio statistics for a grid file
  !-------------------------------------------------------------------------
  SUBROUTINE print_decomposition_statistics(grid_file, decomposition_id)
    CHARACTER(LEN=filename_max), INTENT(in) :: grid_file
    INTEGER, INTENT(in) :: decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    CALL grid_decomposition_statistics(grid_id, decomposition_id)
    CALL delete_grid(grid_id)

  END SUBROUTINE print_decomposition_statistics
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! returns the max domain_id for a given decomposition
  INTEGER FUNCTION get_max_domain_id(grid_id, decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id
    
    TYPE(t_grid_cells),    POINTER :: cells
    
    cells => get_cells(grid_id)
    get_max_domain_id = MAXVAL(cells%domain_id(decomposition_id, 1:cells%no_of_existcells))
    
  END FUNCTION get_max_domain_id
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! returns the number of domains for a given decomposition
  INTEGER FUNCTION get_no_of_domains(grid_id, decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id
    
    TYPE(t_grid_cells),    POINTER :: cells
    
    cells => get_cells(grid_id)
    get_no_of_domains = cells%no_of_domains(decomposition_id)
    
  END FUNCTION get_no_of_domains
  !-------------------------------------------------------------------------
   
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE grid_fill_decomposition_struct(grid_id, decomposition_struct)
    INTEGER, INTENT(in)  :: grid_id
    TYPE(t_decomposition_structure) :: decomposition_struct
    
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
        
    cells => get_cells(grid_id)
    edges => get_edges(grid_id)
    verts => get_vertices(grid_id)

    decomposition_struct%no_of_cells = cells%no_of_existcells
    decomposition_struct%no_of_edges = edges%no_of_existedges
    decomposition_struct%no_of_verts = verts%no_of_existvertices
    
    decomposition_struct%cell_cartesian_center => cells%cartesian_center(:)
    decomposition_struct%cell_geo_center       => cells%center(:)
    decomposition_struct%cells_vertex          => cells%get_vertex_index(:,:)

    decomposition_struct%vertex_geo_coord      => verts%vertex(:)
    
    
    decomposition_struct%no_of_decompositions  = max_decompositions
    decomposition_struct%no_of_domains         => cells%no_of_domains(:)
    decomposition_struct%domain_id             => cells%domain_id(:,:)
      
  END SUBROUTINE grid_fill_decomposition_struct
  !-------------------------------------------------------------------------    
   
  !-------------------------------------------------------------------------
  !>
  !! print decompositio statistics for a grid
  SUBROUTINE grid_decomposition_statistics(grid_id, decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id
    
    TYPE(t_grid_cells),    POINTER :: cells
    INTEGER :: max_domain_id
    INTEGER :: max_neighbors, min_neighbors
    INTEGER :: max_halo_cells, min_halo_cells
    REAL(wp) :: cell_ratio, max_halo_cell_ratio

    INTEGER :: no_of_neigbors, halo_cells_to_edge, halo_cells_to_vert
    INTEGER, ALLOCATABLE :: owned_cells(:), land_cells(:), sea_cells(:)

    INTEGER :: cell_no, this_domain_id, return_status
    TYPE(t_statistic) :: neigbors_statistics, halo_cell_to_verts_statistics, &
      & halo_cell_to_edge_statistics, cell_halo_ratios_statistics, owned_cells_statistics

      
    CHARACTER(*), PARAMETER :: method_name = "print_decomposition_statistics"
    
    CALL new(owned_cells_statistics)
    CALL new(neigbors_statistics)
    CALL new(halo_cell_to_verts_statistics)
    CALL new(halo_cell_to_edge_statistics)
    CALL new(cell_halo_ratios_statistics)
    
    cells => get_cells(grid_id)
    max_domain_id = get_max_domain_id(grid_id, decomposition_id)
    IF (cells%no_of_domains(decomposition_id) /= max_domain_id+1) &
      & CALL finish(method_name, "no_of_domains inconsistent to max_domain_id")
    
    ! calculate owned cells
    ALLOCATE(owned_cells(0:max_domain_id),land_cells(0:max_domain_id), &
      & sea_cells(0:max_domain_id), stat=return_status)
    IF (return_status > 0) THEN
        CALL finish (method_name, 'ALLOCATE(owned_cells)')
    ENDIF
    
    owned_cells(:) = 0
    land_cells(:)  = 0
    sea_cells(:)   = 0
    DO cell_no =1, cells%no_of_existcells
      this_domain_id = cells%domain_id(decomposition_id, cell_no)
      owned_cells(this_domain_id) = owned_cells(this_domain_id) + 1
      IF      (cells%sea_land_mask(cell_no) > 0 ) THEN
        land_cells(this_domain_id) = land_cells(this_domain_id)  + 1
      ELSEIF (cells%sea_land_mask(cell_no) < 0 ) THEN
        sea_cells(this_domain_id) = sea_cells(this_domain_id)  + 1
      ENDIF
    ENDDO

    WRITE(*,*) " === Decomposition Statistics ==="
    
    DO this_domain_id = 0, max_domain_id
      CALL get_subdomain_statistics(grid_id, decomposition_id, this_domain_id, &
        & no_of_neigbors, halo_cells_to_edge, halo_cells_to_vert)

      cell_ratio = REAL(halo_cells_to_vert,wp) / REAL(owned_cells(this_domain_id), wp)
      
      WRITE(*,*) "Domain:", this_domain_id, " no_of_neigbors=", no_of_neigbors, &
        " owned_cells=", owned_cells(this_domain_id), &
        " sea_cells=",   sea_cells(this_domain_id), &
        " land_cells=",  land_cells(this_domain_id), &
        " halo_cells=",  halo_cells_to_vert, &
        " ratio=", cell_ratio,  " halo_cells_to_edge=", halo_cells_to_edge
                     
      CALL add_data_to(owned_cells_statistics, owned_cells(this_domain_id))
      CALL add_data_to(neigbors_statistics, no_of_neigbors)
      CALL add_data_to(halo_cell_to_verts_statistics, halo_cells_to_vert)
      CALL add_data_to(halo_cell_to_edge_statistics, halo_cells_to_edge)
      CALL add_data_to(cell_halo_ratios_statistics, cell_ratio)

    ENDDO
    WRITE(*,*) "----"
    WRITE(*,*) "max/min/mean owned cells=", max(owned_cells_statistics), &
      & min(owned_cells_statistics), &
      & mean(owned_cells_statistics)
    WRITE(*,*) "max/min/mean neighbors=", max(neigbors_statistics), min(neigbors_statistics), &
      & mean(neigbors_statistics)
    WRITE(*,*) "max/min/mean halo cells=", max(halo_cell_to_verts_statistics), &
      & min(halo_cell_to_verts_statistics), mean(halo_cell_to_verts_statistics)
    WRITE(*,*) "max/min/mean halo cells to edges=", max(halo_cell_to_edge_statistics), &
      & min(halo_cell_to_edge_statistics), mean(halo_cell_to_edge_statistics)
    WRITE(*,*) "max/min/mean halo cells ratio=", max(cell_halo_ratios_statistics), &
      & min(cell_halo_ratios_statistics), mean(cell_halo_ratios_statistics)

    WRITE(*,*) " === END Decomposition Statistics ==="

    DEALLOCATE(owned_cells, land_cells, sea_cells)
  
    CALL delete(owned_cells_statistics)
    CALL delete(neigbors_statistics)
    CALL delete(halo_cell_to_verts_statistics)
    CALL delete(halo_cell_to_edge_statistics)
    CALL delete(cell_halo_ratios_statistics)
    
  END SUBROUTINE grid_decomposition_statistics
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! get statistics for a domain
  SUBROUTINE get_subdomain_statistics(grid_id, decomposition_id, in_domain_id, &
    & no_of_neigbors, halo_cells_to_edge, halo_cells_to_vert)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id, in_domain_id
    INTEGER, INTENT(out) :: no_of_neigbors, halo_cells_to_edge, halo_cells_to_vert
    
    !-------------------------------------------------------------------------

    TYPE(t_grid),          POINTER :: in_grid
    TYPE(t_grid_cells),    POINTER :: cells
    TYPE(t_grid_edges),    POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    
    INTEGER, PARAMETER   :: max_neighbor_ids = 64! assume we have max of 32 neighbors
    TYPE(t_integer_list) :: neighbor_list
    
    TYPE(t_integer_list) :: halo_cells, boundary_verts

    INTEGER :: edge_no, vert_no, cell_no_1, cell_no_2
    INTEGER :: vert_list_no, cell_in_vertex, bnd_edges
    INTEGER :: return_status
    
    
    CHARACTER(*), PARAMETER :: method_name = "get_subdomain_statistics"

    neighbor_list%allocated_size = max_neighbor_ids
    neighbor_list%list_size      = 0
    ALLOCATE(neighbor_list%value(neighbor_list%allocated_size), stat=return_status)
    IF (return_status > 0) THEN
        CALL finish (method_name, 'ALLOCATE(neighbor_list)')
    ENDIF
    
    halo_cells%allocated_size = max_halos
    boundary_verts%allocated_size = max_halos
    halo_cells%list_size      = 0
    boundary_verts%list_size  = 0
    ALLOCATE(halo_cells%value(halo_cells%allocated_size), &
      & boundary_verts%value(boundary_verts%allocated_size), stat=return_status)
    IF (return_status > 0) THEN
        CALL finish (method_name, 'ALLOCATE(halo_cells, boundary_verts)')
    ENDIF
    

    in_grid  => get_grid(grid_id)
    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts
    bnd_edges = 0
        
    DO edge_no = 1, edges%no_of_existedges
      cell_no_1 = edges%get_cell_index(edge_no, 1)
      cell_no_2 = edges%get_cell_index(edge_no, 2)
      IF (cell_no_1 < 1 .OR. cell_no_2 < 1) CYCLE
      
      IF (cells%domain_id(decomposition_id, cell_no_2) == in_domain_id) THEN
        cell_no_1 = cell_no_2
        cell_no_2 = edges%get_cell_index(edge_no, 1)
      ENDIF
      IF (cells%domain_id(decomposition_id, cell_no_1) /= in_domain_id) CYCLE
      IF (cells%domain_id(decomposition_id, cell_no_2) == in_domain_id) CYCLE
       
      ! here the cell_no_1 belongs to the domain
      ! cell_no_2 is halo that shares the edge_no, which is boundary
      bnd_edges = bnd_edges + 1
      
      ! add the two boundary verts
      return_status = &
        & add_to_list_if_not_exist(boundary_verts, edges%get_vertex_index(edge_no,1))
      return_status = &
        & add_to_list_if_not_exist(boundary_verts, edges%get_vertex_index(edge_no,2))

      ! add the halo cells that share an edge with thid domain
      return_status = add_to_list_if_not_exist(halo_cells, cell_no_2)
        
    ENDDO

!    write(0,*) " bnd_edges=",bnd_edges
    halo_cells_to_edge = halo_cells%list_size
    
    ! now go through the boundary vertices and get all the halo cells
    DO vert_list_no = 1,boundary_verts%list_size
      vert_no = boundary_verts%value(vert_list_no)
      DO cell_in_vertex = 1, verts%max_connectivity
        cell_no_1 = verts%get_cell_index(vert_no,cell_in_vertex)
        IF (cell_no_1 < 1) CYCLE
        IF (cells%domain_id(decomposition_id, cell_no_1) == in_domain_id) CYCLE
        return_status = add_to_list_if_not_exist(halo_cells, cell_no_1)
        ! add (and count) the neighbors
        return_status = add_to_list_if_not_exist(neighbor_list, &
          & cells%domain_id(decomposition_id, cell_no_1))
      ENDDO
    ENDDO
    halo_cells_to_vert = halo_cells%list_size
    no_of_neigbors     = neighbor_list%list_size
    
    ! clean-up
    DEALLOCATE(neighbor_list%value)
    DEALLOCATE(halo_cells%value, boundary_verts%value)

  END SUBROUTINE get_subdomain_statistics
  !-------------------------------------------------------------------------



END MODULE mo_grid_decomposition

