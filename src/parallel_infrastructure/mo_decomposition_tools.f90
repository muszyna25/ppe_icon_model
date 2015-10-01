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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!-------------------------------------------------------------------------------------
#define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))
#define d_normalize(v) v%x=v%x/d_norma_3d(v)
#define d_sqrdistance_3d(v1,v2) DOT_PRODUCT((v1%x-v2%x),(v1%x-v2%x))
!-------------------------------------------------------------------------------------
MODULE mo_decomposition_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: i8, wp
  USE mo_exception,          ONLY: message_text, message, finish, warning
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_math_utilities
  USE mo_util_sort,          ONLY: quicksort
  USE mo_impl_constants,     ONLY: success
  USE mo_dist_dir,           ONLY: t_dist_dir
  USE ppm_extents,           ONLY: extent, extent_start, extent_size

  IMPLICIT NONE

  PUBLIC :: t_decomposition_structure
  PUBLIC :: t_glb2loc_index_lookup

  PUBLIC :: cluster_subdomains
  PUBLIC :: divide_cells_by_location
  PUBLIC :: t_cell_info, sort_cell_info_by_cell_number
  PUBLIC :: reorder_lonlat_subdomains, reorder_latlon_subdomains
  PUBLIC :: pair_opposite_subdomains
  PUBLIC :: decompose_round_robin
  PUBLIC :: decompose_round_robin_opp
  PUBLIC :: get_no_of_cells_per_subdomain
  PUBLIC :: divide_geometric_medial
  PUBLIC :: read_ascii_decomposition
  PUBLIC :: t_grid_domain_decomp_info
  PUBLIC :: get_local_index
  PUBLIC :: get_valid_local_index
  PUBLIC :: generate_glb_to_loc_lookup
  PUBLIC :: init_glb2loc_index_lookup
  PUBLIC :: set_inner_glb_index
  PUBLIC :: set_outer_glb_index
  PUBLIC :: deallocate_glb2loc_index_lookup
  PUBLIC :: partidx_of_elem_uniform_deco, uniform_partition, &
    &       uniform_partition_start

  PRIVATE


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

  TYPE t_cell_info
    INTEGER :: lat !< latitude coordinate, scaled to integer value range
    INTEGER :: lon !< longitude coordinate, scaled to integer value range
    INTEGER :: cell_number !< cell number (for back-sorting at the end)
    INTEGER :: owner !< will be set to the owner
  END TYPE t_cell_info
  !------------------------------

  TYPE t_glb2loc_index_lookup
    INTEGER, ALLOCATABLE :: inner_glb_index(:)
    INTEGER, ALLOCATABLE :: inner_glb_index_to_loc(:)
    INTEGER, ALLOCATABLE :: outer_glb_index(:)
    INTEGER, ALLOCATABLE :: outer_glb_index_to_loc(:)
    INTEGER :: global_size
  END TYPE  t_glb2loc_index_lookup

  TYPE t_grid_domain_decomp_info

    ! Owner mask:
    ! For cells this is the same as decomp_domain(:,:)==0
    ! index1=nproma, index2=1,nblks_c
    ! For edges, this can not be derived from decomp_domain:
    ! edges at the border are assigned the PE with the bigger number
    ! index1=nproma, index2=1,nblks_e    ! For verts, this can not be derived from decomp_domain:
    ! verts at the border are assigned the PE with the bigger number
    ! index1=nproma, index2=1,nblks_v
    LOGICAL, ALLOCATABLE :: owner_mask(:,:)

    ! The following is only used internally for the coupler
    ! and the icon_comm_lib
    INTEGER, ALLOCATABLE :: owner_local(:)

    ! The following is only used internally for domain decomposition
    INTEGER, ALLOCATABLE :: glb_index(:)
    TYPE (t_glb2loc_index_lookup) :: glb2loc_index

    ! Distributed directory containing owner information
    TYPE(t_dist_dir) :: owner_dist_dir

    ! Domain decomposition flag:
    ! decomp_domain==0: inner domain, decomp_domain>0: boundary, decomp_domain<0: undefined
    ! For cells:
    ! 0=owned, 1=shared edge with owned, 2=shared vertex with owned
    ! index1=nproma, index2=1,nblks_c
    ! For edges:
    ! 0=owned, 1=on owned cell=in domain, 2=exaclty one shared vertex with owned cells
    ! index1=nproma, index2=1,nblks_e
    ! For verts:
    ! 0=owned, 1=on owned cell=in domain, 2=on level 1 cells
    ! index1=nproma, index2=1,nblks_v
    INTEGER, POINTER :: decomp_domain(:,:)

    INTEGER, POINTER :: halo_level(:,:)! just points to the decomp_domain as a more accurate name
  END TYPE

  INTERFACE divide_geometric_medial
    MODULE PROCEDURE decomp_geom_medial_ret
    MODULE PROCEDURE decomp_geometric_medial
  END INTERFACE

  INTERFACE get_valid_local_index
    MODULE PROCEDURE get_valid_local_index_prev
    MODULE PROCEDURE get_valid_local_index_next
  END INTERFACE

  PUBLIC :: OPERATOR(==), OPERATOR(/=)

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE t_cell_info_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE t_cell_info_ne
  END INTERFACE OPERATOR(/=)

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
    TYPE(t_cell_info), ALLOCATABLE :: cell_desc(:)
    REAL(wp) :: lat, lon

    !-----------------------------------------------------------------------
    no_of_cells = decomposition_struct%no_of_cells

     ! Fill the cell_desc array
    ALLOCATE(cell_desc(no_of_cells))


    ! domain decomposition for triangular cells with optional splitting into physical domains
    DO cell = 1, no_of_cells

      lat = decomposition_struct%cell_geo_center(cell)%lat
      lon = decomposition_struct%cell_geo_center(cell)%lon

      ! Using the center of the cells for geometric subdivision leads
      ! to "toothed" edges of the subdivision area
      ! Thus we use the minimum lat/lon as subdision criterion.
      IF (lat >= 0._wp) THEN
        DO i = 1, 3
          j_v = decomposition_struct%cells_vertex(i, cell)
          lat = MAX(lat, decomposition_struct%vertex_geo_coord(j_v)%lat)
          lon = MAX(lon, decomposition_struct%vertex_geo_coord(j_v)%lon)
        ENDDO
      ELSE
        DO i = 1, 3
          j_v = decomposition_struct%cells_vertex(i, cell)
          lat = MIN(lat, decomposition_struct%vertex_geo_coord(j_v)%lat)
          lon = MAX(lon, decomposition_struct%vertex_geo_coord(j_v)%lon)
        ENDDO
      ENDIF

      cell_desc(cell)%lat = fxp_lat(lat)
      cell_desc(cell)%lon = fxp_lon(lon)
      cell_desc(cell)%cell_number = cell
      cell_desc(cell)%owner = 0

    ENDDO


    CALL divide_cells_by_location(no_of_cells, &
      & cell_desc(:), 0, decomposition_size-1)

    ! Set owner list,
    ! note that cell_desc(i)%cell_number holds the cell number, i.e. the
    ! original array index before shuffling due to decomposition
    DO i = 1, no_of_cells
      cell = cell_desc(i)%cell_number
      decomposition_struct%domain_id(out_decomposition_id, cell) = cell_desc(i)%owner
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
    TYPE(t_cartesian_coordinates) :: current_center

    INTEGER, POINTER :: new_subdomain_id(:)

    REAL(wp) :: center_weight
    INTEGER :: no_of_domains, remaining_subdomains, nearest_subdomain
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

    INTEGER, POINTER :: cells_group_size(:)
    INTEGER :: cell_no, no_of_domains,in_domain_id
    INTEGER :: max_subdomain_size
    INTEGER :: check_subdomain_partition, return_status
    REAL    :: check_subdomain_partition_real

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
      &  cells_group_size(0:no_of_domains-1), &
      &  stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(next_new_domain")

    ! find the starting position for reditrubuting,
    !  validate the balance
    check_subdomain_partition = subdomain_partition
    IF (check_subdomain_partition > no_of_domains) &
      & check_subdomain_partition = no_of_domains

 !   cells_group_size = max_subdomain_size / subdomain_partition
    check_subdomain_partition_real = REAL(check_subdomain_partition)
    DO in_domain_id=0,no_of_domains-1
      cells_group_size(in_domain_id) = NINT(REAL(cells_per_domain(in_domain_id)) / check_subdomain_partition_real)
!      write(*,*) "cells_group_size=", cells_group_size(in_domain_id), &
!        & cells_per_domain(in_domain_id), REAL(cells_per_domain(in_domain_id)), &
!        & check_subdomain_partition_real
    ENDDO

    ! gets the first round-robin domain id for each domain
    return_status = validate_round_robin(no_of_domains, next_new_domain, remain_domain_space, &
      & cells_per_domain, cells_group_size, max_subdomain_size,  opposite_subdomain_id)

    write(0,*) "input max_subdomain_size=", max_subdomain_size, &
      &  "no_of_domains=", no_of_domains, " new subdomain_partition:", check_subdomain_partition


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
      IF (distributed_cells_per_domain(in_domain_id) >= cells_group_size(in_domain_id)) THEN
        next_new_domain(in_domain_id) = &
          & MOD(next_new_domain(in_domain_id) + 1, no_of_domains)
        distributed_cells_per_domain(in_domain_id) = 0
      ENDIF
    ENDDO

    ! all cells have been re-distributed
    decomposition_struct%no_of_domains(out_decomposition_id) = no_of_domains

    CALL get_no_of_cells_per_subdomain(decomposition_struct, out_decomposition_id,&
      & cells_per_domain)
    max_subdomain_size = MAXVAL(cells_per_domain)
    WRITE(0,*) " round-robin max_subdomain_size=", max_subdomain_size

    DEALLOCATE(cells_per_domain, next_new_domain, remain_domain_space, &
      & cells_group_size, distributed_cells_per_domain)

  END SUBROUTINE decompose_round_robin
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION validate_round_robin(no_of_domains, next_new_domain, remain_domain_space, &
    & cells_per_domain, cells_group_size, max_subdomain_size, opposite_subdomain_id)

    INTEGER , INTENT(in) :: no_of_domains,  max_subdomain_size
    INTEGER, POINTER :: next_new_domain(:), remain_domain_space(:), &
      & cells_per_domain(:), cells_group_size(:)
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
          group_size = cells_group_size(in_domain_id) + &
            & cells_group_size(opposite_subdomain_id(in_domain_id))
        ELSE
          group_size = cells_group_size(in_domain_id)
        ENDIF
      ELSE
        group_size = cells_group_size(in_domain_id)
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
!      DO WHILE(domain_cells >= (group_size / 2))
      DO WHILE(domain_cells > 1)
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

    INTEGER, POINTER :: opposite_subdomain_id(:), new_subdomain_id(:)

    INTEGER :: no_of_domains,in_domain_id, next_domain_id
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

    INTEGER :: no_of_domains,in_domain_id, min_subdomain_id
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

    INTEGER, POINTER :: new_subdomain_id(:), sort_subdomain_id(:)

    INTEGER :: no_of_domains,in_domain_id
    INTEGER :: return_status

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

    INTEGER :: cell_no, no_of_domains
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
  RECURSIVE SUBROUTINE divide_cells_by_location(n_cells,cell_desc,cpu_a,cpu_b)

    INTEGER, INTENT(in) :: n_cells, cpu_a, cpu_b

    TYPE(t_cell_info), INTENT(inout) :: cell_desc(n_cells)

    INTEGER :: cpu_m, n_cells_m, i, j
    INTEGER(i8) :: lat_sum
    INTEGER :: max_lat_i, max_lon_i, min_lat_i, min_lon_i
    REAL(wp) :: max_lat, max_lon, min_lat, min_lon, avglat, scalexp
    !-----------------------------------------------------------------------

    ! If there is only 1 CPU for distribution, we are done

    IF(cpu_a==cpu_b .OR. n_cells == 0) THEN
      cell_desc(:)%owner = cpu_a
      RETURN
    ENDIF

    ! Get geometric extensions and total number of points of all patches
    min_lat_i = HUGE(min_lat_i)
    max_lat_i = -HUGE(max_lat_i)
    min_lon_i = HUGE(min_lon_i)
    max_lon_i = -HUGE(max_lon_i)
    lat_sum = 0_i8
    DO i = 1, n_cells
      min_lat_i = MIN(min_lat_i, cell_desc(i)%lat)
      max_lat_i = MAX(max_lat_i, cell_desc(i)%lat)
      min_lon_i = MIN(min_lon_i, cell_desc(i)%lon)
      max_lon_i = MAX(max_lon_i, cell_desc(i)%lon)
      lat_sum = lat_sum + INT(cell_desc(i)%lat, i8)
    END DO
    ! min_lat_i = MINVAL(cell_desc(:)%lat)
    ! min_lon_i = MINVAL(cell_desc(:)%lon)
    ! max_lat_i = MAXVAL(cell_desc(:)%lat)
    ! max_lon_i = MAXVAL(cell_desc(:)%lon)

    min_lon = flp_lon(min_lon_i)
    min_lat = flp_lat(min_lat_i)
    max_lon = flp_lon(max_lon_i)
    max_lat = flp_lat(max_lat_i)

    ! average latitude in patch
    avglat  = flp_lat(INT(lat_sum/INT(n_cells,i8)))

    ! account somehow for convergence of meridians - this formula is just empiric
    scalexp = 1._wp - MAX(0._wp,ABS(max_lat)-1._wp,ABS(min_lat)-1._wp)
    min_lon = min_lon*(COS(avglat))**scalexp
    max_lon = max_lon*(COS(avglat))**scalexp

    ! Get dimension with biggest distance from min to max
    ! and sort cells in this dimension

    IF(max_lat-min_lat >= max_lon-min_lon) THEN
      CALL sort_cell_info_by_lat(cell_desc, n_cells)
    ELSE
      CALL sort_cell_info_by_lon(cell_desc, n_cells)
    ENDIF

    ! CPU number where to split CPU set

    cpu_m = (cpu_a+cpu_b-1)/2

    ! If the number of CPUs is not even, we have to split the set
    ! of cells accordingly into to differently sized halfes
    ! in order to end with an equal number of points on every CPU.
    ! Note that DOUBLE arithmetic is used since the integer size
    ! may be exceeded in this calculation!

    n_cells_m = INT(DBLE(n_cells)*DBLE((cpu_b-cpu_a+1)/2)/DBLE(cpu_b-cpu_a+1))

    ! find cells with same lat/lon coordinate to left and right ...
    i = n_cells_m
    j = n_cells_m
    IF (max_lat-min_lat >= max_lon-min_lon) THEN
      DO WHILE (i > 1)
        IF (cell_desc(i - 1)%lat /= cell_desc(n_cells_m)%lat) EXIT
        i = i - 1
      END DO
      DO WHILE (j < n_cells)
        IF (cell_desc(j + 1)%lat /= cell_desc(n_cells_m)%lat) EXIT
        j = j + 1
      END DO
    ELSE
      DO WHILE (i > 1)
        IF (cell_desc(i - 1)%lon /= cell_desc(n_cells_m)%lon) EXIT
        i = i - 1
      END DO
      DO WHILE (j < n_cells)
        IF (cell_desc(j + 1)%lon /= cell_desc(n_cells_m)%lon) EXIT
        j = j + 1
      END DO
    END IF
    ! ... and sort those by cell number to make assignment to partition
    ! repeatable in parallel version
    CALL sort_cell_info_by_cell_number(cell_desc(i:j), j - i + 1)

    ! If there are only two CPUs, we are done

    IF(cpu_b == cpu_a+1) THEN
      cell_desc(1:n_cells_m)%owner = cpu_a
      cell_desc(n_cells_m+1:n_cells)%owner = cpu_b
      RETURN
    ENDIF

    ! Further divide both halves recursively

    CALL divide_cells_by_location(n_cells_m,cell_desc(1:n_cells_m),&
      cpu_a,cpu_m)
    CALL divide_cells_by_location(n_cells-n_cells_m,cell_desc(n_cells_m+1:n_cells),&
      cpu_m+1,cpu_b)

  END SUBROUTINE divide_cells_by_location

  !-------------------------------------------------------------------------
  !> generates the information required to get the local index for a given
  !  global index. The decomp_info%glb_index array needs to be set in order
  !  for this routine to work.
  !  The routine can be used reconstruct the respective information.
  !  Originally this information is generated by
  !  mo_setup_subdivion->divide_patch. The data generated by this routine might
  !  slightly differ from the original data...
  SUBROUTINE generate_glb_to_loc_lookup(n, n_inner, decomp_info)

    INTEGER, INTENT(IN)    :: n        ! Number of local points
    INTEGER, INTENT(IN)    :: n_inner  ! Number of inner global points
    TYPE(t_grid_domain_decomp_info), INTENT(INOUT) :: decomp_info

    INTEGER :: i

    CALL set_inner_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(1:n_inner), &
      &                      (/(i, i = 1, n_inner)/))
    CALL set_outer_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(n_inner+1:), &
      &                      (/(i, i = n_inner+1, n)/))

  END SUBROUTINE generate_glb_to_loc_lookup

  !-------------------------------------------------------------------------

  PURE FUNCTION binary_search(array, key)

    INTEGER, INTENT(IN) :: array(:), key
    INTEGER :: binary_search

    INTEGER :: lb, ub, middle

    lb = 1
    ub = SIZE(array)
    middle = ub / 2

    IF (ub == 0) THEN
      binary_search = 0
      RETURN
    END IF

    DO WHILE (ub >= lb)

      middle = (ub + lb) / 2;

      IF (array(middle) < key) THEN
        lb = middle + 1
      ELSE IF (array(middle) > key) THEN
        ub = middle - 1
      ELSE
        EXIT
      END IF
    END DO

    IF (array(middle) == key) THEN
      binary_search = middle
    ELSE IF (array(middle) > key) THEN
      binary_search = -middle + 1
    ELSE
      binary_search = -middle
    END IF
  END FUNCTION binary_search

  !-------------------------------------------------------------------------
  ! quicksort implementations to sort an array of type(t_cell_info)
  ! according to the lat, lon and cell_number fields

  ! returns the local index for a given global index
  ! in case the global index is not available locally -1 is returned
  ! in case the global index is no valid 0 is returned
  ELEMENTAL FUNCTION get_local_index(glb2loc_index, glb_index)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(in) :: glb_index
    INTEGER :: get_local_index

    INTEGER :: temp

    IF (glb_index > glb2loc_index%global_size .OR. glb_index < 1) THEN
      get_local_index = 0
    ELSE
      ! find in outer indices
      temp = binary_search(glb2loc_index%outer_glb_index(:), glb_index)
      IF (temp > 0) temp = glb2loc_index%outer_glb_index_to_loc(temp)
      ! find in inner indices
      IF (temp <= 0) THEN
        temp = binary_search(glb2loc_index%inner_glb_index(:), glb_index)
        IF (temp > 0) temp = glb2loc_index%inner_glb_index_to_loc(temp)
      END IF

      get_local_index = MERGE(temp, -1, temp > 0)
    END IF

  END FUNCTION get_local_index

  ! returns the local index for a given global index
  ! in case the global index is in the valid range but locally not
  ! available, the index of the next smaller index is return, in
  ! case there is no smaller index or the global index is
  ! invalid 0 is returned
  ELEMENTAL FUNCTION get_valid_local_index_prev(glb2loc_index, glb_index)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(in) :: glb_index
    INTEGER :: get_valid_local_index_prev

    INTEGER :: temp

    IF (glb_index > glb2loc_index%global_size .OR. glb_index < 1) THEN

      get_valid_local_index_prev = 0
    ELSE

      ! find in outer indices
      temp = binary_search(glb2loc_index%outer_glb_index(:), glb_index)
      IF (temp > 0) THEN
        get_valid_local_index_prev = glb2loc_index%outer_glb_index_to_loc(temp)
      ELSE
        ! find in inner indices
        temp = binary_search(glb2loc_index%inner_glb_index(:), glb_index)
        IF (temp > 0) THEN
          get_valid_local_index_prev = glb2loc_index%inner_glb_index_to_loc(temp)
        ELSE
          get_valid_local_index_prev = -temp
        END IF
      END IF
    END IF
  END FUNCTION get_valid_local_index_prev

  ! returns the local index for a given global index
  ! in case the global index is in the valid range but locally not
  ! available, the index of the next greater index is return, in
  ! case there is no greater index the maximum local index + 1, in
  ! case the global index is invalid 0 is returned
  ELEMENTAL FUNCTION get_valid_local_index_next(glb2loc_index, glb_index, use_next)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(in) :: glb_index
    INTEGER :: get_valid_local_index_next
    LOGICAL, INTENT(in) :: use_next

    INTEGER :: temp

    IF (glb_index > glb2loc_index%global_size .OR. &
      & glb_index < 1) THEN

      get_valid_local_index_next = 0
    ELSE

      ! find in outer indices
      temp = binary_search(glb2loc_index%outer_glb_index(:), glb_index)
      IF (temp > 0) THEN
        get_valid_local_index_next = glb2loc_index%outer_glb_index_to_loc(temp)
      ELSE
        ! find in inner indices
        temp = binary_search(glb2loc_index%inner_glb_index(:), glb_index)
        IF (temp > 0) THEN
          get_valid_local_index_next = glb2loc_index%inner_glb_index_to_loc(temp)
        ELSE
          get_valid_local_index_next = - temp + 1
        END IF
      END IF
    END IF
  END FUNCTION get_valid_local_index_next

  !-------------------------------------------------------------------------

  SUBROUTINE init_glb2loc_index_lookup(glb2loc, global_size)

    TYPE (t_glb2loc_index_lookup), INTENT(OUT) :: glb2loc
    INTEGER, INTENT(IN) :: global_size

    INTEGER :: ist

    ALLOCATE(glb2loc%inner_glb_index(0), &
      &      glb2loc%inner_glb_index_to_loc(0), &
      &      glb2loc%outer_glb_index(0), &
      &      glb2loc%outer_glb_index_to_loc(0), stat=ist)
    IF(ist/=success) &
      CALL finish  ("init_glb2loc_index_lookup", &
        &           'allocate in init_glb2loc_index_lookup failed')

    glb2loc%global_size = global_size

  END SUBROUTINE init_glb2loc_index_lookup

  !-------------------------------------------------------------------------

  SUBROUTINE set_inner_glb_index(glb2loc, glb_index, loc_index)

    TYPE (t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc
    INTEGER, INTENT(IN) :: glb_index(:), loc_index(:)

    INTEGER :: ist, num_indices

    num_indices = SIZE(glb_index(:))

    IF (ANY(glb_index(:) < 1 .OR. glb_index(:) > glb2loc%global_size)) &
      CALL finish("set_inner_glb_index", "invalid global index")

    DEALLOCATE(glb2loc%inner_glb_index, &
      &        glb2loc%inner_glb_index_to_loc, stat=ist)
    IF (ist /= success) &
      CALL finish("set_inner_glb_index", "deallocate failed")

    ALLOCATE(glb2loc%inner_glb_index(num_indices), &
      &      glb2loc%inner_glb_index_to_loc(num_indices), stat=ist)
    IF (ist /= success) &
      CALL finish("set_inner_glb_index", "allocate failed")

    glb2loc%inner_glb_index(:) = glb_index(:)
    glb2loc%inner_glb_index_to_loc(:) = loc_index(1:num_indices)

    CALL quicksort(glb2loc%inner_glb_index(:), &
      &            glb2loc%inner_glb_index_to_loc(:))
  END SUBROUTINE set_inner_glb_index

  !-------------------------------------------------------------------------

  SUBROUTINE set_outer_glb_index(glb2loc, glb_index, loc_index)

    TYPE (t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc
    INTEGER, INTENT(IN) :: glb_index(:), loc_index(:)

    INTEGER :: ist, num_indices

    num_indices = SIZE(glb_index(:))

    IF (ANY(glb_index(:) < 1 .OR. glb_index(:) > glb2loc%global_size)) &
      CALL finish("set_outer_glb_index", "invalid global index")

    DEALLOCATE(glb2loc%outer_glb_index, &
      &        glb2loc%outer_glb_index_to_loc, stat=ist)
    IF (ist /= success) &
      CALL finish("set_outer_glb_index", "deallocate failed")

    ALLOCATE(glb2loc%outer_glb_index(num_indices), &
      &      glb2loc%outer_glb_index_to_loc(num_indices), stat=ist)
    IF (ist /= success) &
      CALL finish("set_outer_glb_index", "allocate failed")

    glb2loc%outer_glb_index(:) = glb_index(:)
    glb2loc%outer_glb_index_to_loc(:) = loc_index(1:num_indices)

    CALL quicksort(glb2loc%outer_glb_index(:), &
      &            glb2loc%outer_glb_index_to_loc(:))
  END SUBROUTINE set_outer_glb_index

  !-------------------------------------------------------------------------

  SUBROUTINE deallocate_glb2loc_index_lookup(glb2loc)

    TYPE (t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc

    INTEGER :: ist

    DEALLOCATE(glb2loc%inner_glb_index, &
      &        glb2loc%inner_glb_index_to_loc, &
      &        glb2loc%outer_glb_index, &
      &        glb2loc%outer_glb_index_to_loc, stat=ist)
    IF(ist/=success) &
      CALL finish  ('deallocate_glb2loc_index_lookup', 'deallocate failed')

  END SUBROUTINE deallocate_glb2loc_index_lookup

  !-------------------------------------------------------------------------

  FUNCTION med3_cell_info_lat(a, i, j, k) RESULT(median)
    TYPE(t_cell_info), INTENT(in) :: a(0:)
    INTEGER, INTENT(in) :: i, j, k
    INTEGER :: median
    IF (a(i)%lat < a(j)%lat) THEN
      IF (a(j)%lat < a(k)%lat) THEN
        median = j
      ELSE
        IF (a(i)%lat < a(k)%lat) THEN
          median = k
        ELSE
          median = i
        END IF
      END IF
    ELSE
      IF (a(j)%lat > a(k)%lat) THEN
        median = j
      ELSE
        IF (a(i)%lat < a(k)%lat) THEN
          median = i
        ELSE
          median = k
        END IF
      END IF
    END IF
  END FUNCTION med3_cell_info_lat

  SUBROUTINE cell_info_vec_swap(a, b, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(1:n), b(1:n)
    TYPE(t_cell_info) :: t
    INTEGER :: i
    DO i = 1, n
      t = a(i)
      a(i) = b(i)
      b(i) = t
    END DO
  END SUBROUTINE cell_info_vec_swap

#define SWAP(i,j) temp = a(i) ; a(i) = a(j); a(j) = temp

  RECURSIVE SUBROUTINE sort_cell_info_by_lat(a, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(0:n - 1)
    TYPE(t_cell_info) :: temp

    LOGICAL :: swap_cnt

    INTEGER :: d, pa, pb, pc, pd, pl, pm, pn, pdiff

    IF (n < 7) THEN
      DO pm = 1, n - 1
        pl = pm
        DO WHILE (pl > 0)
          IF (a(pl - 1)%lat <= a(pl)%lat) EXIT
          SWAP(pl, pl - 1)
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pm = n / 2
    IF (n > 7) THEN
      pl = 0
      pn = n - 1
      IF (n > 40) THEN
        d = n / 8
        pl = med3_cell_info_lat(a, pl, pl + d, pl + 2 * d)
        pm = med3_cell_info_lat(a, pm - d, pm, pm + d)
        pn = med3_cell_info_lat(a, pn - 2 * d, pn - d, pn)
      END IF
      pm = med3_cell_info_lat(a, pl, pm, pn)
    END IF
    SWAP(0, pm)
    pb = 1
    pa = pb
    pd = n - 1
    pc = pd
    swap_cnt = .FALSE.
    DO WHILE (.TRUE.)
      DO WHILE (pb <= pc)
        IF (a(pb)%lat > a(0)%lat) EXIT
        IF (a(pb)%lat == a(0)%lat) THEN
          swap_cnt = .TRUE.
          SWAP(pa, pb)
          pa = pa + 1
        END IF
        pb = pb + 1
      END DO
      DO WHILE (pb <= pc)
        IF (a(pc)%lat < a(0)%lat) EXIT
        IF (a(pc)%lat == a(0)%lat) THEN
          swap_cnt = .TRUE.
          SWAP(pc, pd)
          pd = pd - 1
        END IF
        pc = pc - 1
      END DO
      IF (pb > pc) EXIT
      SWAP(pb, pc)
      swap_cnt = .TRUE.
      pb = pb + 1
      pc = pc - 1
    END DO
    IF (.NOT. swap_cnt) THEN  ! Switch to insertion sort
      DO pm = 1, n - 1
        pl = pm
        DO WHILE(pl > 0)
          IF (a(pl - 1)%lat <= a(pl)%lat) EXIT
          SWAP(pl, pl - 1)
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pn =  n
    pdiff = MIN(pa, pb - pa)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(0:pdiff - 1), a(pb - pdiff:pb-1), pdiff)
    pdiff = MIN(pd - pc, pn - pd - 1)
    if (pdiff > 0) &
         CALL cell_info_vec_swap(a(pb:pb + pdiff - 1), a(pn - pdiff:pn - 1), &
         pdiff)
    pdiff = pb - pa
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lat(a, pdiff)
    pdiff = pd - pc
    ! hope the compiler can tail-recurse and save stack space
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lat(a(pn - pdiff:), pdiff)
  END SUBROUTINE sort_cell_info_by_lat

  RECURSIVE SUBROUTINE sort_cell_info_by_lon(a, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(0:n - 1)
    TYPE(t_cell_info) :: temp

    LOGICAL :: swap_cnt

    INTEGER :: d, pa, pb, pc, pd, pl, pm, pn, pdiff

    IF (n < 7) THEN
      DO pm = 1, n - 1
        pl = pm
        DO WHILE (pl > 0)
          IF (a(pl - 1)%lon <= a(pl)%lon) EXIT
          SWAP(pl, pl - 1)
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pm = n / 2
    IF (n > 7) THEN
      pl = 0
      pn = n - 1
      IF (n > 40) THEN
        d = n / 8
        pl = med3_cell_info_lon(a, pl, pl + d, pl + 2 * d)
        pm = med3_cell_info_lon(a, pm - d, pm, pm + d)
        pn = med3_cell_info_lon(a, pn - 2 * d, pn - d, pn)
      END IF
      pm = med3_cell_info_lon(a, pl, pm, pn)
    END IF
    SWAP(0, pm)
    pb = 1
    pa = pb
    pd = n - 1
    pc = pd
    swap_cnt = .FALSE.
    DO WHILE (.TRUE.)
      DO WHILE (pb <= pc)
        IF (a(pb)%lon > a(0)%lon) EXIT
        IF (a(pb)%lon == a(0)%lon) THEN
          swap_cnt = .TRUE.
          SWAP(pa, pb)
          pa = pa + 1
        END IF
        pb = pb + 1
      END DO
      DO WHILE (pb <= pc)
        IF (a(pc)%lon < a(0)%lon) EXIT
        IF (a(pc)%lon == a(0)%lon) THEN
          swap_cnt = .TRUE.
          SWAP(pc, pd)
          pd = pd - 1
        END IF
        pc = pc - 1
      END DO
      IF (pb > pc) EXIT
      SWAP(pb, pc)
      swap_cnt = .TRUE.
      pb = pb + 1
      pc = pc - 1
    END DO
    IF (.NOT. swap_cnt) THEN  ! Switch to insertion sort
      DO pm = 1, n - 1
        pl = pm
        DO WHILE(pl > 0)
          IF (a(pl - 1)%lon <= a(pl)%lon) EXIT
          SWAP(pl, pl - 1)
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pn =  n
    pdiff = MIN(pa, pb - pa)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(0:pdiff - 1), a(pb - pdiff:pb-1), pdiff)
    pdiff = MIN(pd - pc, pn - pd - 1)
    if (pdiff > 0) &
         CALL cell_info_vec_swap(a(pb:pb + pdiff - 1), a(pn - pdiff:pn - 1), &
         pdiff)
    pdiff = pb - pa
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lon(a, pdiff)
    pdiff = pd - pc
    ! hope the compiler can tail-recurse and save stack space
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lon(a(pn - pdiff:), pdiff)
  END SUBROUTINE sort_cell_info_by_lon

  FUNCTION med3_cell_info_lon(a, i, j, k) RESULT(median)
    TYPE(t_cell_info), INTENT(in) :: a(0:)
    INTEGER, INTENT(in) :: i, j, k
    INTEGER :: median
    IF (a(i)%lon < a(j)%lon) THEN
      IF (a(j)%lon < a(k)%lon) THEN
        median = j
      ELSE
        IF (a(i)%lon < a(k)%lon) THEN
          median = k
        ELSE
          median = i
        END IF
      END IF
    ELSE
      IF (a(j)%lon > a(k)%lon) THEN
        median = j
      ELSE
        IF (a(i)%lon < a(k)%lon) THEN
          median = i
        ELSE
          median = k
        END IF
      END IF
    END IF
  END FUNCTION med3_cell_info_lon

  RECURSIVE SUBROUTINE sort_cell_info_by_cell_number(a, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(0:n - 1)
    TYPE(t_cell_info) :: temp

    LOGICAL :: swap_cnt

    INTEGER :: d, pa, pb, pc, pd, pl, pm, pn, pdiff

    IF (n < 7) THEN
      DO pm = 1, n - 1
        pl = pm
        DO WHILE (pl > 0)
          IF (a(pl - 1)%cell_number <= a(pl)%cell_number) EXIT
          SWAP(pl, pl - 1)
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pm = n / 2
    IF (n > 7) THEN
      pl = 0
      pn = n - 1
      IF (n > 40) THEN
        d = n / 8
        pl = med3_cell_info_cell_number(a, pl, pl + d, pl + 2 * d)
        pm = med3_cell_info_cell_number(a, pm - d, pm, pm + d)
        pn = med3_cell_info_cell_number(a, pn - 2 * d, pn - d, pn)
      END IF
      pm = med3_cell_info_cell_number(a, pl, pm, pn)
    END IF
    SWAP(0, pm)
    pb = 1
    pa = pb
    pd = n - 1
    pc = pd
    swap_cnt = .FALSE.
    DO WHILE (.TRUE.)
      DO WHILE (pb <= pc)
        IF (a(pb)%cell_number > a(0)%cell_number) EXIT
        IF (a(pb)%cell_number == a(0)%cell_number) THEN
          swap_cnt = .TRUE.
          SWAP(pa, pb)
          pa = pa + 1
        END IF
        pb = pb + 1
      END DO
      DO WHILE (pb <= pc)
        IF (a(pc)%cell_number < a(0)%cell_number) EXIT
        IF (a(pc)%cell_number == a(0)%cell_number) THEN
          swap_cnt = .TRUE.
          SWAP(pc, pd)
          pd = pd - 1
        END IF
        pc = pc - 1
      END DO
      IF (pb > pc) EXIT
      SWAP(pb, pc)
      swap_cnt = .TRUE.
      pb = pb + 1
      pc = pc - 1
    END DO
    IF (.NOT. swap_cnt) THEN  ! Switch to insertion sort
      DO pm = 1, n - 1
        pl = pm
        DO WHILE(pl > 0)
          IF (a(pl - 1)%cell_number <= a(pl)%cell_number) EXIT
          SWAP(pl, pl - 1)
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pn =  n
    pdiff = MIN(pa, pb - pa)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(0:pdiff - 1), a(pb - pdiff:pb-1), pdiff)
    pdiff = MIN(pd - pc, pn - pd - 1)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(pb:pb + pdiff - 1), a(pn - pdiff:pn - 1), &
         pdiff)
    pdiff = pb - pa
    IF (pdiff > 1) &
         CALL sort_cell_info_by_cell_number(a, pdiff)
    pdiff = pd - pc
    ! hope the compiler can tail-recurse and save stack space
    IF (pdiff > 1) &
         CALL sort_cell_info_by_cell_number(a(pn - pdiff:), pdiff)
  END SUBROUTINE sort_cell_info_by_cell_number

  FUNCTION med3_cell_info_cell_number(a, i, j, k) RESULT(median)
    TYPE(t_cell_info), INTENT(in) :: a(0:)
    INTEGER, INTENT(in) :: i, j, k
    INTEGER :: median
    IF (a(i)%cell_number < a(j)%cell_number) THEN
      IF (a(j)%cell_number < a(k)%cell_number) THEN
        median = j
      ELSE
        IF (a(i)%cell_number < a(k)%cell_number) THEN
          median = k
        ELSE
          median = i
        END IF
      END IF
    ELSE
      IF (a(j)%cell_number > a(k)%cell_number) THEN
        median = j
      ELSE
        IF (a(i)%cell_number < a(k)%cell_number) THEN
          median = i
        ELSE
          median = k
        END IF
      END IF
    END IF
  END FUNCTION med3_cell_info_cell_number

  FUNCTION t_cell_info_eq(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(t_cell_info), INTENT(in) :: a, b
    p = a%lat == b%lat .AND. a%lon == b%lon &
         .AND. a%cell_number == b%cell_number &
         .AND. a%owner == b%owner
  END FUNCTION t_cell_info_eq

  FUNCTION t_cell_info_ne(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(t_cell_info), INTENT(in) :: a, b
    p = a%lat /= b%lat .OR. a%lon /= b%lon &
         .OR. a%cell_number /= b%cell_number &
         .OR. a%owner /= b%owner
  END FUNCTION t_cell_info_ne

  !> compute start integer of uniform interval partition
  ELEMENTAL FUNCTION uniform_partition_start(set_interval, nparts, &
    &                                        part_idx) RESULT(start)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: part_idx
    INTEGER :: start, part_offset, sym_part_idx
    INTEGER(i8) :: sym_size

    part_offset = INT((INT(extent_size(set_interval), i8) &
         &             * INT(part_idx - 1, i8)) / INT(nparts, i8))
    start = extent_start(set_interval) + part_offset

  END FUNCTION uniform_partition_start

  !> compute nth part of integer set interval
  !!
  !! The interval is divided into roughly same sized sub-intervals
  !! forming a uniform partition.
  !! @param set_interval global domain
  !! @param nparts number of parts to decompose into
  !! @param part_idx number of sub-interval to compute
  !! <tt>SIZE(uniform_partition(i)) == SIZE(uniform_partition(nparts - i + 1))</tt>
  !! @return part range corresponding to part_idx
  ELEMENTAL FUNCTION uniform_partition(set_interval, nparts, part_idx) &
    & RESULT(interval)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: part_idx
    TYPE(extent) :: interval
    INTEGER :: start_part, start_next_part

    start_part = uniform_partition_start(set_interval, nparts, part_idx)
    start_next_part = uniform_partition_start(set_interval, nparts, &
         part_idx + 1)
    interval = extent(start_part, start_next_part - start_part)
  END FUNCTION uniform_partition

  ELEMENTAL FUNCTION partidx_of_elem_uniform_deco(set_interval, nparts, &
    &                                             elem_idx) RESULT(part_idx)
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: nparts, elem_idx
    INTEGER :: part_idx

    part_idx = INT((INT(elem_idx - extent_start(set_interval), i8) &
                    * INT(nparts, i8) + INT(nparts, i8) - 1) &
                   / INT(extent_size(set_interval), i8)) + 1

  END FUNCTION partidx_of_elem_uniform_deco

END MODULE mo_decomposition_tools

