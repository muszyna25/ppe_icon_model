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
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_io_units,           ONLY: nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_timer,              ONLY: new_timer, timer_start, timer_stop, print_timer, delete_timer

  USE mo_local_grid
  USE mo_grid_toolbox,       ONLY: add_to_list_if_not_exist
  USE mo_io_local_grid,      ONLY: read_new_netcdf_grid, write_netcdf_grid, &
    & write_ascii_decomposition

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: decompose_all_cells
  PUBLIC :: redecompose_file_round_robin, redecompose_round_robin
  PUBLIC :: decompose_file_metis, decompose_metis
  PUBLIC :: print_decomposition_statistics
  PUBLIC :: get_max_domain_id, get_no_of_domains
  !-------------------------------------------------------------------------

  INTEGER, PARAMETER   :: max_halos = 4096! assume we have max of 32 neighbors, for statistics
  

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE redecompose_file_round_robin(grid_file, out_ascii_file, &
    & in_decomposition_id, out_decomposition_id)
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    CALL redecompose_round_robin(grid_id, in_decomposition_id, out_decomposition_id)
    CALL write_ascii_decomposition(grid_id, out_decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)
    

  END SUBROUTINE redecompose_file_round_robin
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE redecompose_round_robin(grid_id, in_decomposition_id, out_decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, in_decomposition_id, out_decomposition_id

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER, POINTER :: cells_per_domain(:), next_new_domain(:)
    
    INTEGER :: cell_no, no_of_domains,domain_id
    INTEGER :: return_status
    
    CHARACTER(*), PARAMETER :: method_name = "redecompose_round_robin"

    cells => get_cells(grid_id)
    no_of_domains = cells%no_of_domains(in_decomposition_id)
    write(0,*) "in_decomposition_id, out_decomposition_id:", &
      & in_decomposition_id, out_decomposition_id
    ! first count the number of cells in each decomposition
    ALLOCATE(cells_per_domain(0:no_of_domains-1),&
      & next_new_domain(0:no_of_domains-1), stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(cells_per_domain")
    cells_per_domain(:) = 0  
    DO cell_no = 1, cells%no_of_existcells
      cells_per_domain(cells%get_domain_id(in_decomposition_id, cell_no) ) = &
        & cells_per_domain(cells%get_domain_id(in_decomposition_id, cell_no)) + 1
    ENDDO
    
    ! the next_new_domain indicates what is the next new domain for each cell of
    ! the old domains
    ! initialize next_new_domain
    next_new_domain(0) = 0
    DO domain_id = 0, no_of_domains-2
      next_new_domain(domain_id+1) = &
        & MOD(next_new_domain(domain_id) + cells_per_domain(domain_id), no_of_domains)
    ENDDO
    DO cell_no = 1, cells%no_of_existcells
      domain_id = cells%get_domain_id(in_decomposition_id, cell_no)
      cells%get_domain_id(out_decomposition_id, cell_no) = next_new_domain(domain_id)
      next_new_domain(domain_id) = &
        & MOD(next_new_domain(domain_id) + 1, no_of_domains)
    ENDDO

    ! all cells have been re-distributed
    cells%no_of_domains(out_decomposition_id) = no_of_domains

    DEALLOCATE(cells_per_domain,next_new_domain) 
    
  END SUBROUTINE redecompose_round_robin
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! decompose a grid file using metis
  SUBROUTINE decompose_file_metis(grid_file, dec_size, out_ascii_file, &
    & decomposition_id, use_sea_land_mask )
    
    CHARACTER(LEN=*), INTENT(in) :: grid_file, out_ascii_file
    INTEGER, INTENT(in)  :: dec_size, decomposition_id
    LOGICAL, INTENT(in)  :: use_sea_land_mask

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(grid_file)
    CALL decompose_metis(grid_id, dec_size, decomposition_id, use_sea_land_mask)
    CALL write_ascii_decomposition(grid_id, decomposition_id, out_ascii_file)
    
    CALL delete_grid(grid_id)    

  END SUBROUTINE decompose_file_metis
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! decompose a grid using metis
  SUBROUTINE decompose_metis(grid_id, dec_size, decomposition_id, use_sea_land_mask)
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
    REAL, POINTER :: null_real(:)
#endif

    INTEGER, INTENT(in)  :: grid_id, dec_size, decomposition_id
    LOGICAL, INTENT(in)  :: use_sea_land_mask
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
        metis_edge_weights(adjacency_index) = 1
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
    
    ! call metis
    CALL F_metis_partgraphrecursive(no_of_cells, metis_constrains, &
      & metis_xadj, metis_adjncy,                                &
      & metis_cell_weights, null_int, metis_edge_weights,      &
      & dec_size, null_real, null_real, metis_options, metis_edgecut, metis_color)
! C_NULL_PTR
    cells%no_of_domains(decomposition_id) = dec_size

    !fill the domain ids
!$OMP PARALLEL
!$OMP DO PRIVATE(cell_no)
    DO cell_no = 1, cells%no_of_existcells
      cells%get_domain_id(decomposition_id, cell_no) = metis_color(cell_no - 1)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    DEALLOCATE(metis_xadj,metis_adjncy,metis_color,          &
      & metis_cell_weights,metis_edge_weights,metis_options)

    ! for plotting
!     cells%sea_land_mask(:) = cells%get_domain_id(decomposition_id, :)
!     CALL write_netcdf_grid(grid_id)
#endif

  END SUBROUTINE decompose_metis
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
      cells%get_domain_id(decomposition_id, cell_no) = cell_no - 1
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE decompose_all_cells
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
    get_max_domain_id = MAXVAL(cells%get_domain_id(decomposition_id, 1:cells%no_of_existcells))
    
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
  !! print decompositio statistics for a grid
  SUBROUTINE grid_decomposition_statistics(grid_id, decomposition_id)
    INTEGER, INTENT(in)  :: grid_id, decomposition_id
    
    TYPE(t_grid_cells),    POINTER :: cells
    INTEGER :: max_domain_id
    INTEGER :: max_neighbors, min_neighbors
    INTEGER :: max_halo_cells, min_halo_cells
    REAL(wp) :: cell_ratio, max_halo_cell_ratio

    INTEGER :: no_of_neigbors, halo_cells_to_edge, halo_cells_to_vert
    INTEGER, ALLOCATABLE :: owned_cells(:)

    INTEGER :: cell_no, this_domain_id, return_status

      
    CHARACTER(*), PARAMETER :: method_name = "print_decomposition_statistics"
    
    cells => get_cells(grid_id)
    max_domain_id = get_max_domain_id(grid_id, decomposition_id)
    IF (cells%no_of_domains(decomposition_id) /= max_domain_id+1) &
      & CALL finish(method_name, "no_of_domains incisitent to max_domain_id")
    
    ! calculate owned cells
    ALLOCATE(owned_cells(0:max_domain_id), stat=return_status)
    IF (return_status > 0) THEN
        CALL finish (method_name, 'ALLOCATE(owned_cells)')
    ENDIF
    owned_cells(:) = 0
    DO cell_no =1, cells%no_of_existcells
      owned_cells(cells%get_domain_id(decomposition_id, cell_no)) = &
        & owned_cells(cells%get_domain_id(decomposition_id, cell_no)) + 1
    ENDDO

    WRITE(0,*) " === Decomposition Statistics ==="
    max_neighbors = 0
    min_neighbors = max_domain_id + 2
    max_halo_cells = 0
    min_halo_cells = max_halos + 2
    max_halo_cell_ratio = 0.0_wp
    
    DO this_domain_id = 0, max_domain_id
      CALL get_subdomain_statistics(grid_id, decomposition_id, this_domain_id, &
        & no_of_neigbors, halo_cells_to_edge, halo_cells_to_vert)

      cell_ratio = REAL(halo_cells_to_vert,wp) / REAL(owned_cells(this_domain_id), wp)
      
      WRITE(0,*) "Domain:", this_domain_id, " no_of_neigbors=", no_of_neigbors, &
        " owned_cells=", owned_cells(this_domain_id), " halo_cells=", halo_cells_to_vert, &
        " ratio=", cell_ratio,  " halo_cells_to_edge=", halo_cells_to_edge
               
      max_neighbors = MAX(max_neighbors, no_of_neigbors)
      min_neighbors = MIN(min_neighbors, no_of_neigbors)
      max_halo_cells = MAX(max_halo_cells, halo_cells_to_vert)
      min_halo_cells = MIN(min_halo_cells, halo_cells_to_vert)
      max_halo_cell_ratio = MAX(max_halo_cell_ratio, cell_ratio)

    ENDDO
    WRITE(0,*) "----"
    WRITE(0,*) "max/min neighbors=", max_neighbors, min_neighbors
    WRITE(0,*) "max/min cell halos, max ratio=", max_halo_cells, min_halo_cells, &
      & max_halo_cell_ratio

    WRITE(0,*) " === END Decomposition Statistics ==="

    DEALLOCATE(owned_cells)
  
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
    
    INTEGER, PARAMETER   :: max_neighbor_ids = 32! assume we have max of 32 neighbors
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
      
      IF (cells%get_domain_id(decomposition_id, cell_no_2) == in_domain_id) THEN
        cell_no_1 = cell_no_2
        cell_no_2 = edges%get_cell_index(edge_no, 1)
      ENDIF
      IF (cells%get_domain_id(decomposition_id, cell_no_1) /= in_domain_id) CYCLE
      IF (cells%get_domain_id(decomposition_id, cell_no_2) == in_domain_id) CYCLE
       
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
        IF (cells%get_domain_id(decomposition_id, cell_no_1) == in_domain_id) CYCLE
        return_status = add_to_list_if_not_exist(halo_cells, cell_no_1)
        ! add (and count) the neighbors
        return_status = add_to_list_if_not_exist(neighbor_list, &
          & cells%get_domain_id(decomposition_id, cell_no_1))
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

