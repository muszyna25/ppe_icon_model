!-------------------------------------------------------------------------------------
!>
!! Set of methods for grid decomposition.
!!  It runs only on a single process
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
#define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))
#define d_normalize(v) v%x=v%x/d_norma_3d(v)
#define d_sqrdistance_3d(v1,v2) DOT_PRODUCT((v1%x-v2%x),(v1%x-v2%x))
!-------------------------------------------------------------------------------------
MODULE mo_decomposition_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, finish, warning
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_base_geometry

  IMPLICIT NONE

  PUBLIC :: t_decomposition_structure

  PUBLIC :: cluster_subdomains
  PUBLIC :: reorder_lonlat_subdomains, reorder_latlon_subdomains
  PUBLIC :: pair_opposite_subdomains
  PUBLIC :: decompose_round_robin
  PUBLIC :: decompose_round_robin_opp
  PUBLIC :: get_no_of_cells_per_subdomain
  PUBLIC :: divide_geometric_medial
  PUBLIC :: read_ascii_decomposition

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  !------------------------------
  TYPE t_decomposition_structure
  
    INTEGER :: no_of_cells
    INTEGER :: no_of_verts
    INTEGER :: no_of_edges

    TYPE(t_cartesian_coordinates), POINTER ::  cell_cartesian_center(:) 
    TYPE(t_geographical_coordinates), POINTER :: cell_geo_center(:)
    INTEGER, POINTER :: cells_vertex(:,:)  ! vertices of a cell, DIM(3,no_of_cells)

    TYPE(t_geographical_coordinates), POINTER :: vertex_geo_coord(:)

    !> total number of decomposition hold into this structure
    INTEGER :: no_of_decompositions 
    
    !> The number of domains for each domain decomposition
    INTEGER, POINTER :: no_of_domains(:)
    
    !> Holds the multiple domain ids for domain decompositions
    ! domain_id(no_of_decompositions, no_of_cells)
    INTEGER, POINTER :: domain_id(:,:) ! DIM(no_of_decompositions, no_of_cells)
    
  END TYPE t_decomposition_structure
  !------------------------------

  INTERFACE divide_geometric_medial
    MODULE PROCEDURE decomp_geom_medial_ret
    MODULE PROCEDURE decomp_geometric_medial
  END INTERFACE
  
  !------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_ascii_decomposition(ascii_file_name, cell_owner, no_of_cells)
    INTEGER, POINTER :: cell_owner(:) 
    CHARACTER(LEN=*) :: ascii_file_name 
    INTEGER, INTENT(in) :: no_of_cells
    
    INTEGER :: file_id, return_status, cell_no
    CHARACTER(*), PARAMETER :: method_name = "read_ascii_decomposition"
        
    WRITE(0,*) "Read decomposition from file: ", TRIM(ascii_file_name)
    file_id = find_next_free_unit(10,99)

    OPEN(file_id, FILE=TRIM(ascii_file_name), STATUS='OLD', IOSTAT=return_status)
    IF(return_status /= 0) CALL finish(method_name,&
      & 'Unable to open input file: '//TRIM(ascii_file_name))

    DO cell_no = 1, no_of_cells
      READ(file_id, *, IOSTAT=return_status) cell_owner(cell_no)
      IF(return_status /= 0) CALL finish(method_name,'Error reading: '//TRIM(ascii_file_name))
    ENDDO
    CLOSE(file_id)
    
  END SUBROUTINE read_ascii_decomposition
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_ascii_decomposition(ascii_file_name, cell_owner, no_of_cells)
    INTEGER, POINTER :: cell_owner(:) 
    CHARACTER(LEN=*) :: ascii_file_name 
    INTEGER, INTENT(in) :: no_of_cells


    INTEGER :: file_id, error_status, cell_no
    
    WRITE(message_text,'(a,a)')                          &
      &  'Write decomposition file: ', TRIM(ascii_file_name)
    CALL message ('', TRIM(message_text))
    !----------------------------------------------------------------------
    file_id = find_next_free_unit(100,1000)
    OPEN (file_id, FILE=TRIM(ascii_file_name),IOSTAT=error_status)
    DO cell_no = 1, no_of_cells
      WRITE(file_id,*) cell_owner(cell_no)
    ENDDO
    CLOSE(file_id)    

  END SUBROUTINE write_ascii_decomposition
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE decomp_geom_medial_ret(decomposition_struct, decomposition_size, &
    &  cluster, cells_owner, radiation_onwer, radiation_split_factor)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)    :: decomposition_size   ! Number of processors
    LOGICAL, INTENT(in) :: cluster
    INTEGER, POINTER :: cells_owner(:)
    INTEGER, POINTER, OPTIONAL :: radiation_onwer(:)
    INTEGER, INTENT(in), OPTIONAL :: radiation_split_factor

    INTEGER :: no_of_decompositions, medial_decomposition_id, rad_decomposition_id
    INTEGER :: rad_subdomain_split
    CHARACTER(*), PARAMETER :: method_name = "decompose_geom_medial_ret"

    IF (PRESENT(radiation_onwer)) THEN
      write(0,*) "divide_geometric_medial, redistribute radiation, cluster=", cluster
      no_of_decompositions = 2
    ELSE 
      write(0,*) "divide_geometric_medial, cluster=", cluster
      no_of_decompositions = 1
    ENDIF
    medial_decomposition_id = 1
    rad_decomposition_id = 2
      
      ! temorarily allocate the no_of_domains, domain_id
    CALL allocate_dec(decomposition_struct, no_of_decompositions)
    
    ! the output is in cells_owner, ignore the out_decomposition_id
    CALL decompose_geometric_medial(decomposition_struct, decomposition_size, &
      & medial_decomposition_id)

    IF (cluster) &
      CALL cluster_subdomains(decomposition_struct,  &
        & medial_decomposition_id, medial_decomposition_id)
  
    ! fill the decomposition
    CALL fill_onwers_array(decomposition_struct, medial_decomposition_id, cells_owner)
    
    !-----------------------------
    IF (PRESENT(radiation_onwer)) THEN    
      IF (PRESENT(radiation_split_factor)) THEN
        rad_subdomain_split = radiation_split_factor
      ELSE
        rad_subdomain_split = 6
      ENDIF

      CALL reorder_lonlat_subdomains(decomposition_struct,  &
        & in_decomposition_id = medial_decomposition_id,    &
        & out_decomposition_id = rad_decomposition_id)
          
      CALL decompose_round_robin_opp(decomposition_struct,  &
        & in_decomposition_id = rad_decomposition_id,    &
        & out_decomposition_id = rad_decomposition_id,      &
        & subdomain_partition = rad_subdomain_split)
    
      CALL fill_onwers_array(decomposition_struct, rad_decomposition_id, radiation_onwer)
      
    ENDIF
    
    ! delallocate
    CALL deallocate_dec(decomposition_struct)
      
  END SUBROUTINE decomp_geom_medial_ret
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE decomp_geometric_medial(decomposition_struct, decomposition_size, &
    & out_decomposition_id, cluster )
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)    :: decomposition_size   ! Number of processors
    INTEGER, INTENT(in) ::  out_decomposition_id
    LOGICAL, INTENT(in) :: cluster
    
    CALL decompose_geometric_medial(decomposition_struct, decomposition_size, &
      & out_decomposition_id)
    IF (cluster) &
      CALL cluster_subdomains(decomposition_struct,  &
        & out_decomposition_id, out_decomposition_id)

  END SUBROUTINE decomp_geometric_medial
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE decompose_geometric_medial(decomposition_struct, decomposition_size, &
    & out_decomposition_id)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)    :: decomposition_size   ! Number of processors
    INTEGER, INTENT(in)    ::  out_decomposition_id   
    ! (-1 for cells not in subset)

    INTEGER :: cell, i, j_v, no_of_cells
    REAL(wp), ALLOCATABLE :: cell_desc(:,:), workspace(:,:)

    !-----------------------------------------------------------------------
    no_of_cells = decomposition_struct%no_of_cells

    ! Fill the cell_desc array, it must contain:
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner

    ALLOCATE(cell_desc(4,no_of_cells))


   ! domain decomposition for triangular cells with optional splitting into physical domains
   DO cell = 1, no_of_cells

      cell_desc(1,cell) = decomposition_struct%cell_geo_center(cell)%lat
      cell_desc(2,cell) = decomposition_struct%cell_geo_center(cell)%lon

      ! Using the center of the cells for geometric subdivision leads
      ! to "toothed" edges of the subdivision area
      ! Thus we use the minimum lat/lon as subdision criterion.
      IF (cell_desc(1,cell) >= 0._wp) THEN
        DO i = 1, 3
          j_v = decomposition_struct%cells_vertex(i, cell)
          cell_desc(1,cell) = MAX(cell_desc(1,cell), &
            & decomposition_struct%vertex_geo_coord(j_v)%lat)
          cell_desc(2,cell) = MAX(cell_desc(2,cell), &
            & decomposition_struct%vertex_geo_coord(j_v)%lon)
        ENDDO
      ELSE
        DO i = 1, 3
          j_v = decomposition_struct%cells_vertex(i, cell)
          cell_desc(1,cell) = MIN(cell_desc(1,cell), &
            & decomposition_struct%vertex_geo_coord(j_v)%lat)
          cell_desc(2,cell) = MAX(cell_desc(2,cell), &
            & decomposition_struct%vertex_geo_coord(j_v)%lon)
        ENDDO
      ENDIF

      cell_desc(3,cell) = REAL(cell,wp)
      cell_desc(4,cell) = 0.0_wp

    ENDDO


    ALLOCATE(workspace(4,no_of_cells))
    
    CALL divide_cells_by_location(no_of_cells, &
      & cell_desc(:,:), workspace, 0, decomposition_size-1)

    DEALLOCATE(workspace)

    ! Set owner list,
    ! note that the  cell_desc(3,i) holds the actual cell number
    DO i = 1, no_of_cells
      cell = NINT(cell_desc(3,i))
      decomposition_struct%domain_id(out_decomposition_id, cell) = NINT(cell_desc(4,i))
    ENDDO
    decomposition_struct%no_of_domains(out_decomposition_id) = decomposition_size

    DEALLOCATE(cell_desc)

  END SUBROUTINE decompose_geometric_medial
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Cluster the subdomains using the geometrically nearest neigbor (on the sphere)
  !! cluster_size is not used
  SUBROUTINE cluster_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    TYPE(t_cartesian_coordinates), POINTER :: subdomain_center(:)
    TYPE(t_cartesian_coordinates) :: current_center, v
    
    INTEGER, POINTER :: new_subdomain_id(:)
    
    REAL(wp) :: min_distance, distance, center_weight
    INTEGER :: no_of_domains,in_domain_id,remaining_subdomains, nearest_subdomain
    INTEGER :: next_subdomain
    INTEGER :: return_status
    
    CHARACTER(*), PARAMETER :: method_name = "cluster_subdomains"

    no_of_domains = decomposition_struct%no_of_domains(in_decomposition_id)
    write(0,*) method_name, ":in_decomposition_id, out_decomposition_id=", &
      & in_decomposition_id, out_decomposition_id
      
    ! first calculate the center of each subdomain
    NULLIFY(subdomain_center)
    CALL calclulate_subdomain_centers(decomposition_struct, in_decomposition_id, &
      & subdomain_center)
                             
    ALLOCATE( new_subdomain_id(0:no_of_domains-1), &
      & stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(new_subdomain_id")
    new_subdomain_id(:) = -1

    ! assign the 0 subdomain
    new_subdomain_id(0) = 0
    current_center      = subdomain_center(0)
    center_weight = 1.0_wp
    remaining_subdomains= no_of_domains-1
    
!     DO WHILE(remaining_subdomains > cluster_size)
    DO next_subdomain = 1, no_of_domains-1
    
      ! find the next closest free subdomain
      nearest_subdomain = find_nearest_neigbor(subdomain_center, &
          & current_center, no_of_domains, new_subdomain_id)

      IF (nearest_subdomain < 0) &
        CALL finish(method_name, "cannot find nearest_subdomain")
      
      new_subdomain_id(nearest_subdomain) = next_subdomain
      current_center%x = current_center%x + &
        & center_weight * subdomain_center(nearest_subdomain)%x
      d_normalize(current_center)
      center_weight = center_weight + 1.0_wp
      
      remaining_subdomains = remaining_subdomains - 1
    ENDDO
              
    CALL fill_redecomposition(decomposition_struct, in_decomposition_id, &
      & out_decomposition_id, new_subdomain_id, no_of_domains)
              
    DEALLOCATE(subdomain_center, new_subdomain_id)
    
  END SUBROUTINE cluster_subdomains
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  !! allgigning the oppisite subdomains
  SUBROUTINE decompose_round_robin_opp(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id, subdomain_partition)
    
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)  :: subdomain_partition, in_decomposition_id, out_decomposition_id

    INTEGER, POINTER  :: opposite_subdomain_id(:)

    NULLIFY(opposite_subdomain_id)
    CALL find_opposite_subdomains(decomposition_struct,  &
      & in_decomposition_id = in_decomposition_id, &
      & opposite_subdomain_id = opposite_subdomain_id)    
    CALL decompose_round_robin(decomposition_struct, &
      & in_decomposition_id  = in_decomposition_id,  &
      & out_decomposition_id = out_decomposition_id, &
      & subdomain_partition  = subdomain_partition,  &
      & opposite_subdomain_id = opposite_subdomain_id)
    
  END SUBROUTINE decompose_round_robin_opp
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it re-decomposes in a round_robin way for each subdomain
  SUBROUTINE decompose_round_robin(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id, subdomain_partition, &
    & opposite_subdomain_id)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)  ::  subdomain_partition, &
      & in_decomposition_id, out_decomposition_id

    !> if present allign the oppsite subdomain
    INTEGER, POINTER, OPTIONAL  :: opposite_subdomain_id(:)

    INTEGER, POINTER :: cells_per_domain(:), next_new_domain(:), remain_domain_space(:)
    INTEGER, POINTER :: distributed_cells_per_domain(:)
    
    INTEGER :: cells_group_size
    INTEGER :: cell_no, no_of_domains,in_domain_id
    INTEGER :: max_subdomain_size, next_avail_subdomain, domain_cells
    INTEGER :: check_subdomain_partition, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "decompose_round_robin"

    
    no_of_domains = decomposition_struct%no_of_domains(in_decomposition_id)
    write(0,*) method_name, "in_decomposition_id, out_decomposition_id:", &
      & in_decomposition_id, out_decomposition_id, " subdomain_partition=",subdomain_partition
    
    ! first count the number of cells in each decomposition
    NULLIFY(cells_per_domain)
    CALL get_no_of_cells_per_subdomain(decomposition_struct, in_decomposition_id,&
      & cells_per_domain)
    max_subdomain_size = MAXVAL(cells_per_domain)

    ! the next_new_domain indicates what is the next new domain for each cell of
    ! the old domains
    ALLOCATE(next_new_domain(0:no_of_domains-1), &
      &  remain_domain_space(0:no_of_domains-1), &
      &  distributed_cells_per_domain(0:no_of_domains-1), &
      &  stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(next_new_domain")

    ! find the starting position for reditrubuting,
    !  validate the balance
    check_subdomain_partition = subdomain_partition
    IF (check_subdomain_partition > no_of_domains) &
      & check_subdomain_partition = no_of_domains
      
    cells_group_size = max_subdomain_size / subdomain_partition
    return_status = validate_round_robin(no_of_domains, next_new_domain, remain_domain_space, &
      & cells_per_domain, cells_group_size,max_subdomain_size,  opposite_subdomain_id)
      
    write(0,*) "max_subdomain_size=", max_subdomain_size, &
      & "cells_group_size:",  cells_group_size, "no_of_domains=", no_of_domains, &
      & " new subdomain_partition:", check_subdomain_partition

    
     IF (return_status > 0) THEN
       CALL warning(method_name, "cannot find balanced decomposition")
     ELSEIF (return_status < 0) THEN
       CALL finish(method_name, "cannot find balanced decomposition")
     ENDIF
     
!     IF (return_status /= 0) THEN
!       ! revalidate
!       check_subdomain_partition = MOD(no_of_domains,subdomain_partition)
!       IF (check_subdomain_partition > subdomain_partition / 2) THEN
!         check_subdomain_partition = 2*subdomain_partition - &
!           & MOD(no_of_domains,2*subdomain_partition)
!       ELSE   
!         check_subdomain_partition = subdomain_partition - MOD(no_of_domains,subdomain_partition)
!       ENDIF
!       cells_group_size = max_subdomain_size / check_subdomain_partition
!       write(0,*) "max_subdomain_size=", max_subdomain_size, &
!         & "cells_group_size:",  cells_group_size, "no_of_domains=", no_of_domains, &
!         & " new subdomain_partition:", check_subdomain_partition
!       return_status = validate_round_robin(no_of_domains, next_new_domain, remain_domain_space, &
!         & cells_per_domain, cells_group_size, max_subdomain_size,  opposite_subdomain_id)
!       IF (return_status /= 0) THEN
!         ! we cannot find a suitable decomposition
!         CALL finish(method_name, "cannot find balanced decomposition")
!       ENDIF
!     ENDIF
      
    
    ! re-distribute in groups of cells_group_size
    ! the current group size for each subdomain
    !  will be hold in cells_per_domain
    distributed_cells_per_domain(:) = 0
    remain_domain_space(:) = max_subdomain_size
    DO cell_no = 1, decomposition_struct%no_of_cells
      in_domain_id = decomposition_struct%domain_id(in_decomposition_id, cell_no)
!       write(0,*) cell_no, in_domain_id, next_new_domain(in_domain_id)
      decomposition_struct%domain_id(out_decomposition_id, cell_no) = next_new_domain(in_domain_id)
      remain_domain_space(next_new_domain(in_domain_id)) = &
        remain_domain_space(next_new_domain(in_domain_id)) - 1
!       write(*,*) "cell no=", cell_no, in_domain_id, next_new_domain(in_domain_id), &
!         & remain_domain_space(next_new_domain(in_domain_id))
      IF (remain_domain_space(next_new_domain(in_domain_id)) < 0) THEN
        write(0,*) "in_domain_id=", in_domain_id, " new=", next_new_domain(in_domain_id), &
          & remain_domain_space(next_new_domain(in_domain_id))
        CALL warning(method_name, "remain_domain_space(next_new_domain(in_domain_id)) < 0")
        IF (remain_domain_space(next_new_domain(in_domain_id)) < -4 ) THEN
          CALL finish(method_name, "round robin decomposition cannot balance")
        ENDIF
      ENDIF
      distributed_cells_per_domain(in_domain_id) = distributed_cells_per_domain(in_domain_id) + 1
      IF (distributed_cells_per_domain(in_domain_id) >= cells_group_size) THEN
        next_new_domain(in_domain_id) = &
          & MOD(next_new_domain(in_domain_id) + 1, no_of_domains)
        distributed_cells_per_domain(in_domain_id) = 0
      ENDIF
    ENDDO

    ! all cells have been re-distributed
    decomposition_struct%no_of_domains(out_decomposition_id) = no_of_domains

    DEALLOCATE(cells_per_domain,next_new_domain, remain_domain_space, &
      & distributed_cells_per_domain)
    
  END SUBROUTINE decompose_round_robin
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------    
  INTEGER FUNCTION validate_round_robin(no_of_domains, next_new_domain, remain_domain_space, &
    & cells_per_domain, cells_group_size, max_subdomain_size, opposite_subdomain_id)

    INTEGER , INTENT(in) :: no_of_domains, cells_group_size, max_subdomain_size
    INTEGER, POINTER :: next_new_domain(:), remain_domain_space(:), &
      & cells_per_domain(:)
    INTEGER, POINTER, OPTIONAL :: opposite_subdomain_id(:)

      
    INTEGER :: in_domain_id, next_avail_subdomain, counter, domain_cells
    INTEGER :: group_size

    validate_round_robin = 0
    ! determine the group size to be assigned to ditsributed
    
    next_new_domain(:) = -1
    remain_domain_space(:) = max_subdomain_size
    
    next_avail_subdomain = 0
    DO in_domain_id = 0, no_of_domains-1
      ! if it's alread assigned a new subdomain start, skip it
      IF (next_new_domain(in_domain_id) >= 0) CYCLE
      ! Find the next avaliable subdomain with sufficient size
      IF (PRESENT(opposite_subdomain_id)) THEN
        IF (opposite_subdomain_id(in_domain_id) >= 0) THEN
          group_size = cells_group_size * 2
        ELSE
          group_size = cells_group_size
        ENDIF
      ELSE
        group_size = cells_group_size
      ENDIF
      
      counter = 0
      DO WHILE(remain_domain_space(next_avail_subdomain) < group_size .AND. &
               & counter < no_of_domains)
        next_avail_subdomain = MOD(next_avail_subdomain + 1, no_of_domains)
        counter = counter + 1
      ENDDO
      IF (counter >= no_of_domains) THEN
        validate_round_robin = -1
        RETURN !ERROR
      ENDIF

      ! assign the next avail subdomain to in_domain_id
      next_new_domain(in_domain_id) = next_avail_subdomain
      IF (PRESENT(opposite_subdomain_id)) THEN
        IF (opposite_subdomain_id(in_domain_id) >= 0) THEN
          next_new_domain(opposite_subdomain_id(in_domain_id)) = next_avail_subdomain
          domain_cells = cells_per_domain(in_domain_id) + &
            & cells_per_domain(opposite_subdomain_id(in_domain_id))
        ELSE
          domain_cells = cells_per_domain(in_domain_id)
        ENDIF
      ELSE
        domain_cells = cells_per_domain(in_domain_id)
      ENDIF

      ! compute remaining sizes
      ! this is just for checking, should be removed eventually
      DO WHILE(domain_cells >= (group_size / 2))
!         write(*,*) "in_decomposition_id=", in_domain_id, opposite_subdomain_id(in_domain_id), &
!           &  " goes to ",next_avail_subdomain, &
!           & " remain space=", remain_domain_space(next_avail_subdomain), " - ", &
!           & MIN(group_size, domain_cells)
          
        remain_domain_space(next_avail_subdomain) = &
          & remain_domain_space(next_avail_subdomain) - &
          & MIN(group_size, domain_cells)
        domain_cells = domain_cells - group_size
        IF (remain_domain_space(next_avail_subdomain) < 0) &
          & validate_round_robin = 1 ! ERROR
        next_avail_subdomain = MOD(next_avail_subdomain + 1, no_of_domains)
      ENDDO
    ENDDO ! in_domain_id = 0, no_of_domains-1

    validate_round_robin = 0 ! OK
    
  END FUNCTION validate_round_robin
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it reorders the decomposition in pairs of
  !! opposite geographically domains
  SUBROUTINE pair_opposite_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id

    TYPE(t_cartesian_coordinates), POINTER :: subdomain_center(:)
    TYPE(t_cartesian_coordinates) :: v, v2
    INTEGER, POINTER :: opposite_subdomain_id(:), new_subdomain_id(:)

    REAL(wp) :: min_distance, distance
    INTEGER :: cell_no, no_of_domains,in_domain_id, next_domain_id, min_subdomain_id
    INTEGER :: new_no_of_domains, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "pair_opposite_subdomains"

    no_of_domains = decomposition_struct%no_of_domains(in_decomposition_id)
    write(0,*) method_name, ":in_decomposition_id, out_decomposition_id=", &
      & in_decomposition_id, out_decomposition_id
      
    ! first find opposite subdomains
    CALL find_opposite_subdomains(decomposition_struct,&
      & in_decomposition_id, opposite_subdomain_id)
    
    ALLOCATE(new_subdomain_id(0:no_of_domains-1), stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(new_subdomain_id")
    new_subdomain_id(:)      = -1
         
    ! now combine the opposite submdomains
    next_domain_id = 0
    DO in_domain_id = 0, no_of_domains-1      
      IF (opposite_subdomain_id(in_domain_id) < in_domain_id) CYCLE
      new_subdomain_id(in_domain_id)                        = next_domain_id
      new_subdomain_id(opposite_subdomain_id(in_domain_id)) = next_domain_id
      next_domain_id = next_domain_id + 1
    ENDDO
    new_no_of_domains = next_domain_id
!     ! give a new domain id to domains not paired
!     IF (non_paired_subdomains > 0) THEN
!       DO in_domain_id = 0, no_of_domains-1      
!         IF (opposite_subdomain_id(in_domain_id) < 0) THEN
!           new_subdomain_id(in_domain_id) = next_domain_id
!           next_domain_id = next_domain_id + 1
!         ENDIF
!       ENDDO
!     ENDIF
    
    ! Fill the new domain ids in the cells decomposition
    CALL fill_redecomposition(decomposition_struct, in_decomposition_id, &
      & out_decomposition_id, new_subdomain_id, new_no_of_domains)

    DEALLOCATE(opposite_subdomain_id, new_subdomain_id)
    
  END SUBROUTINE pair_opposite_subdomains
  !-------------------------------------------------------------------------  
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE fill_redecomposition(decomposition_struct, &
    & in_decomposition_id, out_decomposition_id, &
    & new_subdomain_id, new_no_of_domains)

    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in) :: in_decomposition_id, out_decomposition_id, new_no_of_domains
    INTEGER, POINTER ::  new_subdomain_id(:)

    INTEGER :: cell_no, in_domain_id
    CHARACTER(*), PARAMETER :: method_name = "fill_redecomposition"

    DO cell_no = 1, decomposition_struct%no_of_cells
      in_domain_id = decomposition_struct%domain_id(in_decomposition_id, cell_no)
      IF (new_subdomain_id(in_domain_id) < 0 .OR. &
        & new_subdomain_id(in_domain_id) > new_no_of_domains) &
        & CALL finish (method_name, "non-valid new_subdomain_id")
        
      decomposition_struct%domain_id(out_decomposition_id, cell_no) = &
        & new_subdomain_id(in_domain_id)
    ENDDO
    decomposition_struct%no_of_domains(out_decomposition_id) = new_no_of_domains

   END SUBROUTINE fill_redecomposition    
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it reorders the decomposition in pairs of
  !! opposite geographically domains
  SUBROUTINE find_opposite_subdomains(decomposition_struct,  &
    & in_decomposition_id, opposite_subdomain_id)
    TYPE(t_decomposition_structure), INTENT(in)  :: decomposition_struct  
    INTEGER, INTENT(in)  :: in_decomposition_id
    INTEGER, POINTER :: opposite_subdomain_id(:)

    TYPE(t_cartesian_coordinates), POINTER :: subdomain_center(:)
    TYPE(t_cartesian_coordinates) :: v
!     TYPE(t_geographical_coordinates) :: geocoord
    INTEGER, POINTER :: new_subdomain_id(:)

    INTEGER :: cell_no, no_of_domains,in_domain_id, next_domain_id, min_subdomain_id
    INTEGER :: return_status
    
    CHARACTER(*), PARAMETER :: method_name = "find_opposite_subdomains"

    no_of_domains = decomposition_struct%no_of_domains(in_decomposition_id)
    write(0,*) method_name, ":in_decomposition_id=", &
      & in_decomposition_id
      
    ! first calculate the center of each subdomain
    NULLIFY(subdomain_center)
    CALL calclulate_subdomain_centers(decomposition_struct, in_decomposition_id, subdomain_center)

    IF (.NOT. ASSOCIATED(opposite_subdomain_id)) THEN
      ALLOCATE(opposite_subdomain_id(0:no_of_domains-1), stat=return_status)
      IF (return_status > 0) &
        & CALL finish (method_name, "ALLOCATE(opposite_subdomain_id")
    ENDIF
    opposite_subdomain_id(:) = -1
         
    ! find opposite subdomain
    ! brute force should be ok since we traverse only subdomains
!     non_paired_subdomains = 0
    DO in_domain_id = 0, no_of_domains-1
      IF (opposite_subdomain_id(in_domain_id) >= 0) CYCLE

        v%x = -subdomain_center(in_domain_id)%x
        min_subdomain_id = find_nearest_neigbor(subdomain_center, v, &
          & no_of_domains, opposite_subdomain_id)
          
        IF (min_subdomain_id >= 0) THEN
          IF (min_subdomain_id == in_domain_id) THEN
            CALL warning(method_name, "has no opposite")
          ELSE
            opposite_subdomain_id(in_domain_id)     = min_subdomain_id
            opposite_subdomain_id(min_subdomain_id) = in_domain_id
          ENDIF
!           geocoord = cc2gc(subdomain_center(in_domain_id))
!           WRITE(0,*) geocoord%lon, geocoord%lat, "paired: ", in_domain_id, min_subdomain_id,&
!           WRITE(0,*)  "paired: ", in_domain_id, min_subdomain_id,&
!             & v%x, subdomain_center(min_subdomain_id)%x
        ELSE
          WRITE(0,*)method_name, ":domain_id=", in_domain_id, " has no opposite!"
          CALL finish(method_name, "Cannot find opposite")
!           non_paired_subdomains = non_paired_subdomains + 1
        ENDIF
    ENDDO
    
    DEALLOCATE(subdomain_center)
    
  END SUBROUTINE find_opposite_subdomains
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE reorder_lonlat_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id)
    TYPE(t_decomposition_structure), INTENT(in)  :: decomposition_struct  
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id
    
    CALL reorder_geo_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id, 2)
          
  END SUBROUTINE reorder_lonlat_subdomains
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE reorder_latlon_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id)
    TYPE(t_decomposition_structure), INTENT(in)  :: decomposition_struct  
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id
    
    CALL reorder_geo_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id, 1)
          
  END SUBROUTINE reorder_latlon_subdomains
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it reorders the decomposition in geographical coordinates
  SUBROUTINE reorder_geo_subdomains(decomposition_struct,  &
    & in_decomposition_id, out_decomposition_id, latlon_order)
    TYPE(t_decomposition_structure), INTENT(in)  :: decomposition_struct  
    INTEGER, INTENT(in)  :: in_decomposition_id, out_decomposition_id
    
    INTEGER, INTENT(in)  :: latlon_order ! 1=max lat to min lat, max lon to min lon
                                          ! 2=max lon to min lon, max lat to min lat

    TYPE(t_cartesian_coordinates), POINTER :: subdomain_center(:)
    TYPE(t_geographical_coordinates), POINTER :: subdomain_geo_coordinates(:)
    TYPE(t_cartesian_coordinates) :: v, v2
    
    INTEGER, POINTER :: opposite_subdomain_id(:), new_subdomain_id(:), sort_subdomain_id(:)
    
    INTEGER :: cell_no, no_of_domains,in_domain_id, next_domain_id, max_subdomain
    INTEGER :: new_no_of_domains, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "reorder_geo_subdomains"

    no_of_domains = decomposition_struct%no_of_domains(in_decomposition_id)
    write(0,*) method_name, ":in_decomposition_id, out_decomposition_id=", &
      & in_decomposition_id, out_decomposition_id

      
    ! first calculate the center of each subdomain
    NULLIFY(subdomain_center)
    CALL calclulate_subdomain_centers(decomposition_struct, in_decomposition_id, subdomain_center)
    NULLIFY(subdomain_geo_coordinates)
    CALL cartesian_to_geographical(subdomain_center, no_of_domains, subdomain_geo_coordinates, &
      & start_point=0, end_point=no_of_domains-1)
    DEALLOCATE(subdomain_center)

    ! calculate sort value
    SELECT CASE (latlon_order)
    CASE (1)
      CALL bubble_sort_lat_lon(subdomain_geo_coordinates, 0, no_of_domains-1, sort_subdomain_id)      
    CASE (2)
      CALL bubble_sort_lon_lat(subdomain_geo_coordinates, 0, no_of_domains-1, sort_subdomain_id)
    CASE DEFAULT
      CALL finish(method_name, 'Unkown lat_lon_order')
    END SELECT
    
    DEALLOCATE(subdomain_geo_coordinates)
                             
    ALLOCATE( new_subdomain_id(0:no_of_domains-1), &
      & stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(subdomain_sort_value")

    DO in_domain_id=0, no_of_domains-1
      new_subdomain_id(sort_subdomain_id(in_domain_id)) = in_domain_id
    ENDDO
      
    CALL fill_redecomposition(decomposition_struct, in_decomposition_id, out_decomposition_id, &
      & new_subdomain_id, no_of_domains)
              
    DEALLOCATE(sort_subdomain_id, new_subdomain_id)
    
  END SUBROUTINE reorder_geo_subdomains
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! count the number of cells in each decomposition
  SUBROUTINE get_no_of_cells_per_subdomain(decomposition_struct, &
    & in_decomposition_id, cells_per_domain)
    TYPE(t_decomposition_structure), INTENT(in)  :: decomposition_struct  
    INTEGER, INTENT(in)  ::  in_decomposition_id
    INTEGER, POINTER :: cells_per_domain(:)
    
    INTEGER :: cell_no, no_of_domains,domain_id
    INTEGER :: return_status
    
    CHARACTER(*), PARAMETER :: method_name = "get_no_of_cells_per_subdomain"

    no_of_domains = decomposition_struct%no_of_domains(in_decomposition_id)
    
    ! count the number of cells in each decomposition
    IF (.NOT. ASSOCIATED(cells_per_domain)) THEN
      ALLOCATE(cells_per_domain(0:no_of_domains-1), stat=return_status)
      IF (return_status > 0) &
        & CALL finish (method_name, "ALLOCATE(cells_per_domain")
    ENDIF
    cells_per_domain(:) = 0  
    DO cell_no = 1, decomposition_struct%no_of_cells
      cells_per_domain(decomposition_struct%domain_id(in_decomposition_id, cell_no) ) = &
        & cells_per_domain(decomposition_struct%domain_id(in_decomposition_id, cell_no)) + 1
    ENDDO
        
  END SUBROUTINE get_no_of_cells_per_subdomain
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Given a decomposition, it reorders the decomposition in pairs of
  !! opposite geographically domains
  SUBROUTINE calclulate_subdomain_centers(decomposition_struct,  &
    & decomposition_id, subdomain_centers)
    TYPE(t_decomposition_structure), INTENT(in)  :: decomposition_struct  
    INTEGER, INTENT(in)  :: decomposition_id
    TYPE(t_cartesian_coordinates), POINTER :: subdomain_centers(:)
    
    TYPE(t_cartesian_coordinates) :: v
    TYPE(t_geographical_coordinates) :: geocoord
    INTEGER :: cell_no, no_of_domains, domain_id
    INTEGER :: return_status
    
    CHARACTER(*), PARAMETER :: method_name = "calclulate_subdomain_centers"

    no_of_domains = decomposition_struct%no_of_domains(decomposition_id)
    
    ! first calculate the center of each subdomain
    IF (.NOT. ASSOCIATED(subdomain_centers)) THEN
      ALLOCATE(subdomain_centers(0:no_of_domains-1), stat=return_status)
      IF (return_status > 0) &
        & CALL finish (method_name, "ALLOCATE(subdomain_centers")
    ENDIF
    
    subdomain_centers(:)%x(1) = 0.0_wp
    subdomain_centers(:)%x(2) = 0.0_wp
    subdomain_centers(:)%x(3) = 0.0_wp
     
    ! calculate the center of each subdomain
    DO cell_no = 1, decomposition_struct%no_of_cells
      domain_id = decomposition_struct%domain_id(decomposition_id, cell_no)
      subdomain_centers(domain_id)%x = subdomain_centers(domain_id)%x + &
        decomposition_struct%cell_cartesian_center(cell_no)%x
    ENDDO
    ! normalize the center
    DO domain_id = 0, no_of_domains-1
      subdomain_centers(domain_id)%x = subdomain_centers(domain_id)%x / &
        & d_norma_3d(subdomain_centers(domain_id))
    ENDDO
    
  END SUBROUTINE calclulate_subdomain_centers
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE bubble_sort_lon_lat(geo_coord, start_index, end_index, sorted_list)
    TYPE(t_geographical_coordinates), POINTER :: geo_coord(:)
    INTEGER, INTENT(in)  :: start_index, end_index
    INTEGER, POINTER  :: sorted_list(:)

    LOGICAL :: keep_ordering 
    INTEGER :: i, exch, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "bubble_sort"

    write(0,*) " Sorting in min-max lon-lat..."
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
        IF (  geo_coord(sorted_list(i))%lon > geo_coord(sorted_list(i+1))%lon .OR. &
         &  (geo_coord(sorted_list(i))%lon == geo_coord(sorted_list(i+1))%lon .AND. &
         &   geo_coord(sorted_list(i))%lat == geo_coord(sorted_list(i+1))%lat)  ) THEN
          exch = sorted_list(i)
          sorted_list(i) = sorted_list(i+1)
          sorted_list(i+1) = exch
          keep_ordering = .true.
        ENDIF
      ENDDO
    ENDDO
        
  END SUBROUTINE bubble_sort_lon_lat
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE bubble_sort_lat_lon(geo_coord, start_index, end_index, sorted_list)
    TYPE(t_geographical_coordinates), POINTER :: geo_coord(:)
    INTEGER, INTENT(in)  :: start_index, end_index
    INTEGER, POINTER  :: sorted_list(:)

    LOGICAL :: keep_ordering 
    INTEGER :: i, exch, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "bubble_sort"

    write(0,*) " Sorting in min-max lat-lon..."
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
        IF (  geo_coord(sorted_list(i))%lat > geo_coord(sorted_list(i+1))%lat .OR. &
         &  (geo_coord(sorted_list(i))%lat == geo_coord(sorted_list(i+1))%lat .AND. &
         &   geo_coord(sorted_list(i))%lon == geo_coord(sorted_list(i+1))%lon)  ) THEN
          exch = sorted_list(i)
          sorted_list(i) = sorted_list(i+1)
          sorted_list(i+1) = exch
          keep_ordering = .true.
        ENDIF
      ENDDO
    ENDDO
        
  END SUBROUTINE bubble_sort_lat_lon
  !-------------------------------------------------------------------------
   
   
  !-------------------------------------------------------------------------
  !> works on the unit sphere 
  INTEGER FUNCTION find_nearest_neigbor(in_coord, to_pos, size_of_array, array_mask)
    TYPE(t_cartesian_coordinates), POINTER :: in_coord(:)
    TYPE(t_cartesian_coordinates) :: to_pos
    INTEGER :: size_of_array
    INTEGER, POINTER :: array_mask(:)

    TYPE(t_cartesian_coordinates) :: v
    REAL(wp) :: min_distance, distance
    INTEGER :: element   
    
    min_distance = 100.0_wp ! we assume we are on the sphere
    find_nearest_neigbor = -1
    DO element = 0, size_of_array-1
      IF (array_mask(element) >= 0) CYCLE
      v%x = in_coord(element)%x
      distance = d_sqrdistance_3d(to_pos,v)
      IF (distance < min_distance) THEN
        min_distance =  distance
        find_nearest_neigbor = element
      ENDIF
    ENDDO
  
  END FUNCTION find_nearest_neigbor
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE fill_onwers_array(decomposition_struct, decomposition_id, cells_owner)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in) :: decomposition_id   ! Number of processors
    INTEGER, POINTER    :: cells_owner(:)
    
    INTEGER :: cell, return_status
    CHARACTER(*), PARAMETER :: method_name = "fill_onwers_array"

    IF (.NOT. ASSOCIATED(cells_owner)) THEN
      ALLOCATE( cells_owner(decomposition_struct%no_of_cells), &
        & stat=return_status)
      IF (return_status > 0) &
        & CALL finish (method_name, "ALLOCATE(cells_owner")
    ENDIF
     
    DO cell = 1, decomposition_struct%no_of_cells
      cells_owner(cell) = decomposition_struct%domain_id(decomposition_id, cell)
    ENDDO
  
  END SUBROUTINE fill_onwers_array
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  SUBROUTINE allocate_dec(decomposition_struct, no_of_decompositions)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    INTEGER, INTENT(in)    :: no_of_decompositions   ! Number of processors

    INTEGER :: return_status
    CHARACTER(*), PARAMETER :: method_name = "allocate_dec"
    
    ALLOCATE( decomposition_struct%no_of_domains(no_of_decompositions), &
      & decomposition_struct%domain_id(no_of_decompositions, &
      &   decomposition_struct%no_of_cells), &
      & stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(domain_id")

    decomposition_struct%no_of_decompositions = no_of_decompositions
      
  END SUBROUTINE allocate_dec
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_dec(decomposition_struct)
    TYPE(t_decomposition_structure)  :: decomposition_struct
    
    DEALLOCATE( decomposition_struct%no_of_domains, &
      & decomposition_struct%domain_id)
      
  END SUBROUTINE deallocate_dec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Actually divides geometrically by location on cpu_a .. cpu_b
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  RECURSIVE SUBROUTINE divide_cells_by_location(no_of_cells,cell_desc,work,cpu_a,cpu_b)

    INTEGER, INTENT(in) :: no_of_cells, cpu_a, cpu_b

    REAL(wp), INTENT(inout) :: cell_desc(4,no_of_cells),work(4,no_of_cells)
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner
    !
    ! work contains workspace for sort_array_by_row to avoid local allocation there

    INTEGER cpu_m, n_cells_m
    REAL(wp) :: xmax(2), xmin(2), avglat, scalexp
    !-----------------------------------------------------------------------

    ! If there is only 1 CPU for distribution, we are done

    IF(cpu_a==cpu_b) THEN
      cell_desc(4,:) = REAL(cpu_a,wp)
      RETURN
    ENDIF

    ! Get geometric extensions and total number of points of all patches

    xmin(1) = MINVAL(cell_desc(1,:))
    xmin(2) = MINVAL(cell_desc(2,:))
    xmax(1) = MAXVAL(cell_desc(1,:))
    xmax(2) = MAXVAL(cell_desc(2,:))

    ! average latitude in patch
    avglat  = SUM(cell_desc(1,:))/REAL(no_of_cells,wp)

    ! account somehow for convergence of meridians - this formula is just empiric
    scalexp = 1._wp - MAX(0._wp,ABS(xmax(1))-1._wp,ABS(xmin(1))-1._wp)
    xmin(2) = xmin(2)*(COS(avglat))**scalexp
    xmax(2) = xmax(2)*(COS(avglat))**scalexp

    ! Get dimension with biggest distance from min to max
    ! and sort cells in this dimension

    IF(xmax(1)-xmin(1) >= xmax(2)-xmin(2)) THEN
      CALL sort_array_by_row(cell_desc, work, 1)
    ELSE
      CALL sort_array_by_row(cell_desc, work, 2)
    ENDIF

    ! CPU number where to split CPU set

    cpu_m = (cpu_a+cpu_b-1)/2

    ! If the number of CPUs is not even, we have to split the set
    ! of cells accordingly into to differntly sized halfes
    ! in order to end with an equal number of points on every CPU.
    ! Note that DOUBLE arithmetic is used since the integer size
    ! may be exceeded in this calculation!

    n_cells_m = INT(DBLE(no_of_cells)*DBLE(cpu_m-cpu_a+1)/DBLE(cpu_b-cpu_a+1))

    ! If there are only two CPUs, we are done

    IF(cpu_b == cpu_a+1) THEN
      cell_desc(4,1:n_cells_m)         = REAL(cpu_a,wp)
      cell_desc(4,n_cells_m+1:no_of_cells) = REAL(cpu_b,wp)
      RETURN
    ENDIF

    ! Further divide both halves recursively

    CALL divide_cells_by_location(n_cells_m,cell_desc(:,1:n_cells_m),work(:,1:n_cells_m),&
      cpu_a,cpu_m)
    CALL divide_cells_by_location(no_of_cells-n_cells_m,cell_desc(:,n_cells_m+1:no_of_cells),&
      work(:,n_cells_m+1:no_of_cells),cpu_m+1,cpu_b)

  END SUBROUTINE divide_cells_by_location
  !-------------------------------------------------------------------------

 
  !-------------------------------------------------------------------------
  !>
  !! Special quicksort implementation for sorting a 2D array by one selected row.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  RECURSIVE SUBROUTINE sort_array_by_row(x,y,row)

    !

    INTEGER, INTENT(in) :: row ! number of row for sorting

    REAL(wp), INTENT(inout) :: x(:,:) ! array to be sorted
    REAL(wp), INTENT(inout) :: y(:,:) ! workspace

    REAL(wp) :: p
    INTEGER :: n, ipiv, ix, iy, i
    !-----------------------------------------------------------------------

    n = SIZE(x,2)

    IF(n<=1) RETURN

    ipiv = (n+1)/2
    p = x(row,ipiv)
    ix = 0
    iy = 1

#ifdef __SX__
    IF (n >= 12) THEN ! Use vectorized version
      DO i=1,n
        IF(i==ipiv) CYCLE
        IF(x(row,i) < p) THEN
          ix = ix+1
          y(1:4,ix) = x(1:4,i)
        ENDIF
      ENDDO

      y(1:4,ix+1) = x(1:4,ipiv) ! Store pivot

      DO i=1,n
        IF(i==ipiv) CYCLE
        IF(x(row,i) >= p) THEN
          iy = iy+1
          y(1:4,ix+iy) = x(1:4,i)
        ENDIF
      ENDDO
!CDIR COLLAPSE
      x(:,:) = y(:,:)
    ELSE  ! use non-vectorized version
      y(1:4,1) = x(1:4,ipiv) ! Store pivot

      DO i=1,n
        IF(i==ipiv) CYCLE
        IF(x(row,i) < p) THEN
          ix = ix+1
          x(1:4,ix) = x(1:4,i)
        ELSE
          iy = iy+1
          y(1:4,iy) = x(1:4,i)
        ENDIF
      ENDDO

      x(1:4,ix+1:ix+iy) = y(1:4,1:iy)
    ENDIF
#else
    y(1:4,1) = x(1:4,ipiv) ! Store pivot

    DO i=1,n
      IF(i==ipiv) CYCLE
      IF(x(row,i) < p) THEN
        ix = ix+1
        x(1:4,ix) = x(1:4,i)
      ELSE
        iy = iy+1
        y(1:4,iy) = x(1:4,i)
      ENDIF
    ENDDO

    x(1:4,ix+1:ix+iy) = y(1:4,1:iy)
#endif

    ipiv = ix+1 ! New pivot location

    IF(ipiv>2)   CALL sort_array_by_row(x(:,:ipiv-1),y(:,:ipiv-1),row)
    IF(ipiv<n-1) CALL sort_array_by_row(x(:,ipiv+1:),y(:,ipiv+1:),row)

  END SUBROUTINE sort_array_by_row
  !-------------------------------------------------------------------------




END MODULE mo_decomposition_tools

