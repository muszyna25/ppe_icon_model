#ifdef __PGI
!pgi$g opt=1
#endif
!>
!!               This module provides all routines for dividing patches.
!!
!!               This module provides all routines for dividing patches
!! (including interpolation state) and setting up communication.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_setup_subdivision
  !
  !-------------------------------------------------------------------------
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell,  &
    & min_rledge, max_rledge, min_rlvert, max_rlvert, max_phys_dom,  &
    & min_rlcell_int, min_rledge_int, min_rlvert_int, max_hw
  USE mo_exception,          ONLY: finish, message, get_filename_noext

  USE mo_run_config,         ONLY: msg_level
  USE mo_io_units,           ONLY: filename_max
  USE mo_model_domain,       ONLY: t_patch, p_patch_local_parent, t_pre_patch, &
       c_num_edges, c_parent, c_child, c_phys_id, c_neighbor, c_edge, &
       c_vertex
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info, &
    &                              get_local_index, get_valid_local_index, &
    &                              set_inner_glb_index, set_outer_glb_index, &
    &                              uniform_partition
  USE mo_mpi,                ONLY: p_bcast, proc_split, p_max
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_COMM_NULL, p_int, &
       mpi_in_place, mpi_success, mpi_sum, mpi_lor, p_bool
#endif
  USE mo_mpi,                ONLY: p_comm_work, p_int, &
    & p_pe_work, p_n_work, my_process_is_mpi_parallel, p_alltoall, p_alltoallv

  USE mo_parallel_config,       ONLY:  nproma, ldiv_phys_dom, &
    & division_method, division_file_name, n_ghost_rows, &
    & div_geometric, ext_div_from_file

  USE mo_communication,      ONLY: blk_no, idx_no, idx_1d
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, patch_weight
  USE mo_alloc_patches,ONLY: allocate_basic_patch, allocate_remaining_patch, &
                             deallocate_pre_patch
  USE mo_decomposition_tools, ONLY: divide_cells_by_location, t_cell_info, &
    & sort_cell_info_by_cell_number, read_ascii_decomposition, &
    & partidx_of_elem_uniform_deco
  USE mo_math_utilities,      ONLY: fxp_lat, fxp_lon
  USE mo_master_control,      ONLY: get_my_process_type, ocean_process, &
                                    testbed_process
#ifndef NOMPI
  USE mo_divide_cells_by_location_mpi, ONLY: divide_cells_by_location_mpi, &
       init_divide_cells_by_location_mpi
#endif
  USE mo_grid_config,         ONLY: use_dummy_cell_closure
  USE mo_util_sort,           ONLY: quicksort, insertion_sort
  USE mo_dist_dir,            ONLY: dist_dir_setup
  USE ppm_distributed_array,  ONLY: dist_mult_array, global_array_desc, &
    &                               dist_mult_array_new, &
    &                               dist_mult_array_local_ptr, &
    &                               dist_mult_array_expose, &
    &                               dist_mult_array_delete, dist_mult_array_get, &
    &                               ppm_int
  USE ppm_extents,            ONLY: extent

  IMPLICIT NONE

  PRIVATE

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: decompose_domain


  TYPE nb_flag_list_elem
    INTEGER, ALLOCATABLE :: idx(:)
    INTEGER, ALLOCATABLE :: owner(:)
  END TYPE nb_flag_list_elem

CONTAINS

  !------------------------------------------------------------------
  !>
  !!  Divide patches and interpolation states.
  !!
  !!
  !!  @note Despite its name, this subroutine is also called in
  !!        sequential runs where it simply copies (and initializes)
  !!        the patch data structure.
  !!
  SUBROUTINE decompose_domain( p_patch, p_patch_pre )

    TYPE(t_patch), INTENT(INOUT), TARGET :: p_patch(n_dom_start:)
    TYPE(t_pre_patch), INTENT(INOUT), TARGET :: p_patch_pre(n_dom_start:)

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_subdivision:decompose_domain")

    ! Local variables:
    INTEGER :: jg, jgp, jc, jgc, n, &
      &        n_procs_decomp
    INTEGER :: nprocs(p_patch_pre(1)%n_childdom)
    REAL(wp) :: weight(p_patch_pre(1)%n_childdom)
    TYPE(t_pre_patch), POINTER :: wrk_p_parent_patch_pre
    ! LOGICAL :: l_compute_grid
    INTEGER :: my_process_type, order_type_of_halos
    INTEGER, POINTER :: local_ptr(:)
    TYPE(dist_mult_array) :: dist_cell_owner, dist_cell_owner_p

    CALL message(routine, 'start of domain decomposition')

    ! This subroutine has 2 different operating modes:
    !
    !
    ! - If the subroutine is called from a parallel run,
    !   then the domain is split according to the number of processors
    !
    ! - When called from a single processor run, the domain is
    !   not split but simply copied (and the p_patch_local_parent is
    !   initialized properly).
    !
    n_procs_decomp = p_n_work


    ! -----------------------------------------------------------------------------
    ! Check for processor splitting

    ! Check if the processor set should be split for patches of the 1st generation.
    ! This is the case if patch_weight > 0 for at least one of the root's childs.
    ! For clarity, we require that patch_weight=0 for any other patch

    ! Get weights for 1st generation patches
    proc_split = .FALSE.
    DO jc = 1, p_patch_pre(1)%n_childdom
      jgc = p_patch_pre(1)%child_id(jc)
      weight(jc) = patch_weight(jgc)
      IF(weight(jc) > 0._wp) proc_split = .TRUE.
    ENDDO

    ! Check if weights for other patches are 0 (for clarity only)
    IF(patch_weight(1) /= 0._wp) &
      CALL finish(routine,'Weight for root patch must be 0')
    DO jg = 2, n_dom
      jgp = p_patch_pre(jg)%parent_id
      IF(jgp /= 1 .AND. patch_weight(jg) > 0._wp) &
        CALL finish(routine,'Weight for higher level patch must be 0')
    ENDDO
    ! -----------------------------------------------------------------------------

    IF(proc_split) THEN

      ! -----------------------------------------------------------------------------
      ! CASE 1: DECOMPOSE WITH PROCESSOR SPLITTING
      ! -----------------------------------------------------------------------------

      IF(p_pe_work==0) WRITE(0,*) 'Splitting processor grid for first level patches'
      IF(p_pe_work==0) WRITE(0,'(a,10f12.3)') 'Weights for first level patches:',weight(:)

      ! In this case, the working processor set must be at least as big
      ! as the number of childs of the root patch
      IF(p_patch_pre(1)%n_childdom > n_procs_decomp) &
        CALL finish(routine,'Too few procs for processor grid splitting')

      ! Get the number of procs per patch according to weight(:).
      ! Every patch gets at least 1 proc (of course!):
      nprocs(:) = 1

      ! The remaining procs are divided among patches similar to
      ! the d'Hondt method for elections

      DO n = p_patch_pre(1)%n_childdom+1, n_procs_decomp
        jg = MAXLOC(weight(:)/REAL(nprocs(:)+1,wp),1)
        nprocs(jg) = nprocs(jg)+1
      ENDDO

      IF(p_pe_work==0) THEN
        WRITE(0,*) 'Processor splitting:'
        DO jc = 1, p_patch_pre(1)%n_childdom
          jgc =  p_patch_pre(1)%child_id(jc)
          WRITE(0,'(a,i0,a,f10.3,a,i0,a,i0)')                                              &
            &   'Patch ',jgc,' weight ',weight(jc),' gets ',nprocs(jc),' of ',n_procs_decomp
          IF (nprocs(jc) <= 1) &
            CALL finish(routine,'Processor splitting requires at least 2 PEs per patch')
        ENDDO
      ENDIF

      ! Set proc0, n_proc for all patches ...

      ! ... for the root patch and patch 0 if it exists

      p_patch(n_dom_start:1)%n_proc = n_procs_decomp
      p_patch(n_dom_start:1)%proc0  = 0

      ! ... for 1st generation childs

      n = 0
      DO jc = 1, p_patch_pre(1)%n_childdom
        jgc = p_patch_pre(1)%child_id(jc)
        p_patch(jgc)%proc0  = n
        p_patch(jgc)%n_proc = nprocs(jc)
        n = n + nprocs(jc)
      ENDDO

      ! ... for deeper level descendants

      DO jg = 2, n_dom

        jgp = p_patch_pre(jg)%parent_id

        IF(jgp /= 1) THEN
          p_patch(jg)%n_proc = p_patch(jgp)%n_proc
          p_patch(jg)%proc0  = p_patch(jgp)%proc0
        ENDIF

      ENDDO

    ELSE

      ! -----------------------------------------------------------------------------
      ! CASE 2: STANDARD DECOMPOSITION, NO SPLITTING
      ! -----------------------------------------------------------------------------
      ! No splitting, proc0, n_proc are identical for all patches

      IF(p_pe_work==0) WRITE(0,*) 'No splitting of processor grid'
      p_patch(:)%n_proc = n_procs_decomp
      p_patch(:)%proc0  = 0

    ENDIF

#ifdef NOMPI
    p_patch(:)%comm = 0
#else
    p_patch(:)%comm = MPI_COMM_NULL
#endif
    p_patch(:)%rank = -1

    ! -----------------------------------------------------------------------------
    ! Divide patches
    ! -----------------------------------------------------------------------------

    DO jg = n_dom_start, n_dom

      jgp = p_patch_pre(jg)%parent_id

      IF(jg == n_dom_start) THEN
        NULLIFY(wrk_p_parent_patch_pre) ! Must be NULL for patch
      ELSE
        wrk_p_parent_patch_pre => p_patch_pre(jgp)
      ENDIF

      ! Here comes the actual domain decomposition:
      ! Every cells gets assigned an owner.

      IF (my_process_is_mpi_parallel()) THEN
        ! Set division method
        CALL divide_patch_cells(p_patch_pre(jg), jg, p_patch(jg)%n_proc, &
             p_patch(jg)%proc0, dist_cell_owner, dist_cell_owner_p, &
             wrk_p_parent_patch_pre)
      ELSE
        ! trivial "decomposition"
        CALL create_dist_cell_owner(dist_cell_owner, &
          &                         p_patch_pre(jg)%n_patch_cells_g)
        CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_ptr)
        local_ptr = 0
        IF (jg > n_dom_start) &
          CALL divide_parent_cells(p_patch_pre(jg), p_patch_pre(jgp), &
            &                      dist_cell_owner, dist_cell_owner_p)
      END IF

      ! Please note: Previously, for jg==0 no ghost rows were set.
      ! Currently, we need ghost rows for jg==0 also for dividing the int state and grf state
      ! Have still to check if int state/grf state is needed at all for jg==0,
      ! if this is not the case, the ghost rows can be dropped again.


      ! order_type_of_halos = 0 order for parent (l_compute_grid = false)
      !                       1 order root grid  (l_compute_grid = true)
      !                       2 all halos go to the end, for ocean
      my_process_type = get_my_process_type()
      SELECT CASE (my_process_type)
        CASE (ocean_process)
           order_type_of_halos = 2
        CASE default
        order_type_of_halos = 1
      END SELECT

      CALL divide_patch(p_patch(jg), p_patch_pre(jg), dist_cell_owner, &
        &               n_ghost_rows, order_type_of_halos, p_pe_work)

      IF(jg > n_dom_start) THEN

        order_type_of_halos = 0

        CALL divide_patch(p_patch_local_parent(jg), p_patch_pre(jgp), &
          &               dist_cell_owner_p, 1, order_type_of_halos,  p_pe_work)

        CALL dist_mult_array_delete(dist_cell_owner_p)

        CALL set_pc_idx(p_patch(jg), p_patch_pre(jgp))
        IF (jgp > n_dom_start) THEN

          CALL set_pc_idx(p_patch_local_parent(jg), &
            &             p_patch_pre(p_patch_local_parent(jg)%parent_id))
        END IF
      ENDIF

      CALL dist_mult_array_delete(dist_cell_owner)

    ENDDO

    ! The patches may be discarded now
    DO jg = n_dom_start, n_dom
      CALL deallocate_pre_patch( p_patch_pre(jg) )
    ENDDO

  CONTAINS
!-------------------------------------------------------------------------
  !>
  !! This method_name sets the parent-child-index for cells and edges
  !!
  !! @par Revision History
  !! Developed  by Rainer Johanni, Dec 2011
  !!
    SUBROUTINE set_pc_idx(patch, parent_patch_pre)

      TYPE(t_patch), TARGET, INTENT(inout) :: patch
      TYPE(t_pre_patch), TARGET, INTENT(inout) :: parent_patch_pre

    ! local variables

    INTEGER :: jg, jb, jl, nlen, ip, i

    !-----------------------------------------------------------------------

    patch%cells%pc_idx(:,:) = 0
    patch%edges%pc_idx(:,:) = 0

    DO jb = 1, patch%nblks_c

      IF (jb /= patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch%npromz_c
      ENDIF

      DO jl = 1, nlen

        ip = idx_1d(patch%cells%parent_glb_idx(jl,jb), &
                    patch%cells%parent_glb_blk(jl,jb))
        jg = patch%cells%decomp_info%glb_index(idx_1d(jl, jb))

        DO i = 1, 4
          CALL dist_mult_array_get(parent_patch_pre%cells%dist, c_child, (/ip, i/), jc)
          IF (jc == jg) patch%cells%pc_idx(jl,jb) = i
        END DO
!          IF(patch%cells%pc_idx(jl,jb) == 0) CALL finish('set_pc_idx','cells%pc_idx')

      ENDDO

    ENDDO

    DO jb = 1, patch%nblks_e

      IF (jb /= patch%nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = patch%npromz_e
      ENDIF

      DO jl = 1, nlen

        ip = idx_1d(patch%edges%parent_glb_idx(jl,jb), &
                    patch%edges%parent_glb_blk(jl,jb))
        jg = patch%edges%decomp_info%glb_index(idx_1d(jl, jb))

        DO i = 1, 4
          CALL dist_mult_array_get(parent_patch_pre%edges%child, 1, (/ip, i/), jc)
          IF (jc == jg) patch%edges%pc_idx(jl,jb) = i
        END DO
!          IF(patch%edges%pc_idx(jl,jb) == 0) CALL finish('set_pc_idx','edges%pc_idx')

      ENDDO

    ENDDO

  END SUBROUTINE set_pc_idx

  END SUBROUTINE decompose_domain

  SUBROUTINE init_dist_cell_owner(dist_cell_owner, cell_owner)
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    INTEGER, OPTIONAL, INTENT(IN) :: cell_owner(:)

    INTEGER, POINTER :: local_ptr(:)

    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_ptr)
    IF (PRESENT(cell_owner)) THEN
      local_ptr(:) = cell_owner(LBOUND(local_ptr,1):UBOUND(local_ptr,1))
    ELSE
      local_ptr(:) = 0
    END IF

    CALL dist_mult_array_expose(dist_cell_owner)

  END SUBROUTINE init_dist_cell_owner

  SUBROUTINE create_dist_cell_owner(dist_cell_owner, n_cells_g)

    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner
    INTEGER, INTENT(IN) :: n_cells_g

    TYPE(global_array_desc) :: dist_cell_owner_desc(1)
    TYPE(extent) :: local_chunk(1,1)

    dist_cell_owner_desc(1)%a_rank = 1
    dist_cell_owner_desc(1)%rect(1)%first = 1
    dist_cell_owner_desc(1)%rect(1)%size = n_cells_g
    dist_cell_owner_desc(1)%element_dt = ppm_int

    local_chunk(1,1) = uniform_partition(dist_cell_owner_desc(1)%rect(1), &
      &                                  p_n_work, p_pe_work+1)

    dist_cell_owner = dist_mult_array_new(dist_cell_owner_desc, local_chunk, &
      &                                   p_comm_work)
  END SUBROUTINE create_dist_cell_owner

  !-----------------------------------------------------------------------------
  !>
  !! Divides the cells of a patch (in wrk_p_patch_pre) for parallelization.
  !!
  !! Outputs the subdivsion in cell_owner(:) which must already be allocated
  !! to the correct size (wrk_p_patch_pre%n_patch_cells_g)
  !!
  !! If  wrk_p_parent_patch_pre is associated, this indicates that the patch has a parent
  !! which has consquences for subdivision (cell with the same parent must not
  !! get to different PEs).
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Split out as a separate routine, Rainer Johanni, Oct 2010

  SUBROUTINE divide_patch_cells(wrk_p_patch_pre, patch_no, n_proc, proc0, &
       &                        dist_cell_owner, dist_cell_owner_p, &
       &                        wrk_p_parent_patch_pre)

    TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre
    INTEGER, INTENT(IN)  :: patch_no !> The patch number,
                                     ! used to identify patch specific decomposition
    INTEGER, INTENT(IN)  :: n_proc !> Number of processors for split
    INTEGER, INTENT(IN)  :: proc0  !> First processor of patch
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner !> Cell division
    !> Cell division for parent (only computed when wrk_p_parent_patch_pre
    !! associated)
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner_p
    TYPE(t_pre_patch), POINTER :: wrk_p_parent_patch_pre

    INTEGER, POINTER :: cell_owner(:) !> Cell division
    INTEGER, POINTER :: local_owner_ptr(:)
    INTEGER, ALLOCATABLE :: flag_c(:)
    CHARACTER(LEN=filename_max) :: use_division_file_name ! if ext_div_from_file

    ! Please note: Unfortunately we cannot use p_io for doing I/O,
    ! since this might be the test PE which is never calling this routine
    ! (this is the case in the actual setup).
    ! Thus we use the worker PE 0 for I/O and don't use message() for output.

    IF (division_method(patch_no) == ext_div_from_file) THEN

      ALLOCATE(cell_owner(wrk_p_patch_pre%n_patch_cells_g))

      ! only procs 0 will decompose and then broadcast
      IF (p_pe_work == 0) THEN

        IF (division_method(patch_no) == ext_div_from_file) THEN

          IF (division_file_name(patch_no) == "") THEN
            use_division_file_name = &
              & TRIM(get_filename_noext(wrk_p_patch_pre%grid_filename))//'.cell_domain_ids'
          ELSE
            use_division_file_name = division_file_name(patch_no)
          ENDIF

          CALL read_ascii_decomposition(use_division_file_name, cell_owner, &
            &                           wrk_p_patch_pre%n_patch_cells_g)

        ENDIF ! subddivision_method(patch_no)

        ! Quick check for correct values
        IF(MINVAL(cell_owner(:)) < 0 .or. &
           MAXVAL(cell_owner(:)) >= n_proc) THEN
          WRITE(0,*) "n_porc=",n_proc, " MINAVAL=", MINVAL(cell_owner(:)), &
            " MAXVAL=", MAXVAL(cell_owner(:))
          CALL finish('divide_patch','Illegal subdivision in input file')
        ENDIF

      ENDIF ! IF (p_pe_work == 0)

      ! broadcast cell_owner array to other processes
      CALL p_bcast(cell_owner, 0, comm=p_comm_work)

      ! Set processor offset
      cell_owner(:) = cell_owner(:) + proc0

      CALL create_dist_cell_owner(dist_cell_owner, wrk_p_patch_pre%n_patch_cells_g)
      CALL init_dist_cell_owner(dist_cell_owner, cell_owner)
      DEALLOCATE(cell_owner)

      ! Assign the cell owners of the current patch to the parent cells
      ! for the construction of the local parent, set "-1" elsewhere.
      IF (ASSOCIATED(wrk_p_parent_patch_pre)) &
        CALL divide_parent_cells(wrk_p_patch_pre, wrk_p_parent_patch_pre, &
          &                      dist_cell_owner, dist_cell_owner_p)

    ELSE ! IF (division_method(patch_no) == ext_div_from_file)

      IF(division_method(patch_no) == div_geometric) THEN
        IF(ASSOCIATED(wrk_p_parent_patch_pre)) THEN
          ALLOCATE(flag_c(wrk_p_parent_patch_pre%cells%local_chunk(1,1)%first:&
                          wrk_p_parent_patch_pre%cells%local_chunk(1,1)%first+&
        &                 wrk_p_parent_patch_pre%cells%local_chunk(1,1)%size-1))
        ELSE
          ALLOCATE(flag_c(wrk_p_patch_pre%cells%local_chunk(1,1)%first: &
                          wrk_p_patch_pre%cells%local_chunk(1,1)%first + &
        &                 wrk_p_patch_pre%cells%local_chunk(1,1)%size - 1))
        END IF
      ELSE
        CALL finish('divide_patch','Illegal division_method setting')
      END IF

      ! Built-in subdivison methods
      IF(ASSOCIATED(wrk_p_parent_patch_pre)) THEN

        ! Cells with the same parent must not go to different PEs.
        ! Thus we have to divide in reality the subset of the parent cells
        ! which cover the actual patch and then assign the ownership
        ! of the cells of the actual patch according to the parent cells

        ! Determine the subset of the parent patch covering the actual patch
        ! by flagging the according cells
        CALL compute_flag_c(flag_c, wrk_p_patch_pre, wrk_p_parent_patch_pre)

        ! Divide subset of patch
        ! Receives the PE  numbers for every cell
        CALL divide_subset_geometric(flag_c, n_proc, wrk_p_parent_patch_pre, &
             dist_cell_owner_p, ASSOCIATED(wrk_p_parent_patch_pre))

        CALL dist_mult_array_expose(dist_cell_owner_p)

        ! Owners of the cells of the parent patch are now in dist_cell_owner_p.
        ! Set owners in current patch from this

        CALL compute_dist_owner_from_dist_owner_p( &
          dist_cell_owner, wrk_p_patch_pre, wrk_p_parent_patch_pre, &
          dist_cell_owner_p)

      ELSE

        ! Set subset flags where the "subset" is the whole patch
        flag_c = 1

        ! No parent patch, simply divide current patch

        ! Divide complete patch
        CALL divide_subset_geometric(flag_c, n_proc, wrk_p_patch_pre, &
             dist_cell_owner, ASSOCIATED(wrk_p_parent_patch_pre))
      ENDIF

      ! Set processor offset
      CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_owner_ptr)
      local_owner_ptr = local_owner_ptr + proc0
      CALL dist_mult_array_expose(dist_cell_owner)
    ENDIF ! division_method

  CONTAINS

    SUBROUTINE compute_flag_c(flag_c, wrk_p_patch_pre, wrk_p_parent_patch_pre)

      INTEGER, ALLOCATABLE :: flag_c(:)
      TYPE(t_pre_patch), INTENT(IN) :: wrk_p_patch_pre
      TYPE(t_pre_patch), INTENT(IN) :: wrk_p_parent_patch_pre

      INTEGER, POINTER :: local_parent_ptr(:), local_phys_id_ptr(:)
      INTEGER, ALLOCATABLE :: parent_phys_id_map(:,:), num_send(:), num_recv(:), &
        &                     send_displs(:), recv_displs(:), flag_c_prep(:,:)
      INTEGER :: parent_part_rank, i
      TYPE(extent) :: global_extent

      flag_c = 0

      CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_phys_id, &
        &                            local_phys_id_ptr)
      CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_parent, &
        &                            local_parent_ptr)

      ALLOCATE(parent_phys_id_map(2,LBOUND(local_parent_ptr,1):&
        &                           UBOUND(local_parent_ptr,1)))

      parent_phys_id_map(1, :) = local_parent_ptr
      parent_phys_id_map(2, :) = local_phys_id_ptr

      CALL quicksort(parent_phys_id_map(1, :), parent_phys_id_map(2, :))

      global_extent%first = 1
      global_extent%size = wrk_p_parent_patch_pre%n_patch_cells_g

      ALLOCATE(num_send(p_n_work), num_recv(p_n_work), &
        &      send_displs(p_n_work), recv_displs(p_n_work))

      num_send = 0
      DO i = LBOUND(local_parent_ptr, 1), UBOUND(local_parent_ptr, 1)
        parent_part_rank = &
          partidx_of_elem_uniform_deco(global_extent, p_n_work, &
          &                            local_parent_ptr(i))
        num_send(parent_part_rank) = num_send(parent_part_rank) + 1
      END DO

      CALL p_alltoall(num_send, num_recv, p_comm_work)

      send_displs(1) = 0
      recv_displs(1) = 0
      DO i = 2, p_n_work
        send_displs(i) = send_displs(i-1) + num_send(i-1)
        recv_displs(i) = recv_displs(i-1) + num_recv(i-1)
      END DO
      ALLOCATE(flag_c_prep(2,recv_displs(p_n_work)+num_recv(p_n_work)))
      CALL p_alltoallv(parent_phys_id_map, num_send, send_displs, &
        &              flag_c_prep, num_recv, recv_displs, p_comm_work)

      DO i = 1, SIZE(flag_c_prep, 2)
        flag_c(flag_c_prep(1,i)) = MAX(1, flag_c_prep(2,i))
      END DO
    END SUBROUTINE compute_flag_c

    SUBROUTINE compute_dist_owner_from_dist_owner_p(dist_cell_owner, &
      wrk_p_patch_pre, wrk_p_parent_patch_pre, dist_cell_owner_p)

      TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner
      TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre
      TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_parent_patch_pre
      TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner_p

      INTEGER, POINTER :: local_parent_ptr(:), local_owner_p_ptr(:), &
        &                 local_owner_ptr(:)
      INTEGER :: i, parent_part_rank
      INTEGER, ALLOCATABLE :: parent_cell_idx(:), parent_cell_idx_permutation(:), &
        &                     num_send(:), num_recv(:), send_displs(:), &
        &                     recv_displs(:), inquired_parent_cell_idx(:)
      TYPE(extent) :: global_extent

      CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_parent, &
        &                            local_parent_ptr)
      ALLOCATE(parent_cell_idx(LBOUND(local_parent_ptr,1): &
        &                      UBOUND(local_parent_ptr,1)), &
        &      parent_cell_idx_permutation(LBOUND(local_parent_ptr,1): &
        &                                  UBOUND(local_parent_ptr,1)))
      parent_cell_idx = local_parent_ptr
      parent_cell_idx_permutation = &
        (/(i,i=LBOUND(parent_cell_idx,1),UBOUND(parent_cell_idx,1))/)

      CALL quicksort(parent_cell_idx, parent_cell_idx_permutation)

      global_extent%first = 1
      global_extent%size = wrk_p_parent_patch_pre%n_patch_cells_g

      ALLOCATE(num_send(p_n_work), num_recv(p_n_work), &
        &      send_displs(p_n_work), recv_displs(p_n_work))
      num_send = 0
      DO i = LBOUND(parent_cell_idx,1), UBOUND(parent_cell_idx,1)
        parent_part_rank = &
          partidx_of_elem_uniform_deco(global_extent, p_n_work, &
          &                            parent_cell_idx(i))
        num_send(parent_part_rank) = num_send(parent_part_rank) + 1
      END DO

      CALL p_alltoall(num_send, num_recv, p_comm_work)

      send_displs(1) = 0
      recv_displs(1) = 0
      DO i = 2, p_n_work
        send_displs(i) = send_displs(i-1) + num_send(i-1)
        recv_displs(i) = recv_displs(i-1) + num_recv(i-1)
      END DO
      ALLOCATE(inquired_parent_cell_idx(recv_displs(p_n_work)+num_recv(p_n_work)))
      CALL p_alltoallv(parent_cell_idx, num_send, send_displs, &
        &              inquired_parent_cell_idx, num_recv, recv_displs, p_comm_work)

      CALL dist_mult_array_local_ptr(dist_cell_owner_p, 1, local_owner_p_ptr)

      inquired_parent_cell_idx = local_owner_p_ptr(inquired_parent_cell_idx)

      CALL p_alltoallv(inquired_parent_cell_idx, num_recv, recv_displs, &
        &              parent_cell_idx, num_send, send_displs, p_comm_work)

      CALL create_dist_cell_owner(dist_cell_owner, &
        &                         wrk_p_patch_pre%n_patch_cells_g)
      CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_owner_ptr)

      local_owner_ptr(parent_cell_idx_permutation) = parent_cell_idx

    END SUBROUTINE compute_dist_owner_from_dist_owner_p

  END SUBROUTINE divide_patch_cells

  !-------------------------------------------------------------------------------------------------
  !>
  !! Sets the owner for the division of the cells of the parent patch
  !! with the same subdivision as the child

  SUBROUTINE divide_parent_cells(p_patch_pre, p_parent_patch_pre, &
    &                            dist_cell_owner, dist_cell_owner_p)

    !> Patch for which the parent should be divided
    TYPE(t_pre_patch), INTENT(INOUT) :: p_patch_pre
    !> Parent patch
    TYPE(t_pre_patch), INTENT(INOUT) :: p_parent_patch_pre
    !> Ownership of cells in p_patch_pre
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    !> Output: Cell division for parent. Must be allocated to n_patch_cells of!
    !! the global parent
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner_p

    INTEGER :: child_chunk_first, child_chunk_last, child_chunk_size, i
    INTEGER :: n_cells_parent_g, parent_part_rank
    TYPE(global_array_desc) :: dist_cell_owner_p_desc(1)
    TYPE(extent) :: parent_chunk(1,1)
    INTEGER, POINTER :: parent_chunk_ptr(:), child_chunk_ptr(:), &
      &                 child_chunk_parent_ptr(:)
    INTEGER, ALLOCATABLE :: num_parent_owner_send(:), num_parent_owner_recv(:),&
      &                     parent_owner_send_displ(:), &
      &                     parent_owner_recv_displ(:), &
      &                     parent_cell_owner_map(:,:), parent_chunk_prep(:,:)

    ! get local chunk
    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, child_chunk_ptr)
    child_chunk_first = LBOUND(child_chunk_ptr, 1)
    child_chunk_size = SIZE(child_chunk_ptr)
    child_chunk_last = child_chunk_first + child_chunk_size - 1
    CALL dist_mult_array_local_ptr(p_patch_pre%cells%dist, c_parent, &
      &                            child_chunk_parent_ptr)

    n_cells_parent_g = p_parent_patch_pre%n_patch_cells_g

    dist_cell_owner_p_desc(1)%a_rank = 1
    dist_cell_owner_p_desc(1)%rect(1)%first = 1
    dist_cell_owner_p_desc(1)%rect(1)%size = n_cells_parent_g
    dist_cell_owner_p_desc(1)%element_dt = ppm_int

    ! get parent cells global index for all child cells of child chunk and their
    ! respective owners (parent cells can occur more than once)
    ALLOCATE(parent_cell_owner_map(2, child_chunk_first:child_chunk_last))
    parent_cell_owner_map(1,:) = child_chunk_ptr
    parent_cell_owner_map(2,:) = child_chunk_parent_ptr

    ! sort both arrays according to parent cell global index
    CALL quicksort(parent_cell_owner_map(2,:), parent_cell_owner_map(1,:))

    ! compute number of parent cells per process
    ALLOCATE(num_parent_owner_send(0:p_n_work-1), &
      &      num_parent_owner_recv(0:p_n_work-1), &
      &      parent_owner_send_displ(0:p_n_work-1), &
      &      parent_owner_recv_displ(0:p_n_work-1))
    num_parent_owner_send = 0
    DO i = child_chunk_first, child_chunk_last
      ! compute for each parent cell global index the respective process
      parent_part_rank = &
        partidx_of_elem_uniform_deco(dist_cell_owner_p_desc(1)%rect(1), p_n_work, &
        &                            parent_cell_owner_map(2,i))-1
      num_parent_owner_send(parent_part_rank) = &
        num_parent_owner_send(parent_part_rank) + 1
    END DO

    ! alltoall number of parent cells per process
    CALL p_alltoall(num_parent_owner_send, num_parent_owner_recv, p_comm_work)

    ! alltoallv parent cell global index and owners
    parent_owner_send_displ(0) = 0
    parent_owner_recv_displ(0) = 0
    DO i = 1, p_n_work-1
      parent_owner_send_displ(i) = parent_owner_send_displ(i-1) + &
        &                          num_parent_owner_send(i-1)
      parent_owner_recv_displ(i) = parent_owner_recv_displ(i-1) + &
        &                          num_parent_owner_recv(i-1)
    END DO
    ALLOCATE(parent_chunk_prep(2,parent_owner_recv_displ(p_n_work-1) + &
      &                          num_parent_owner_recv(p_n_work-1)))
    CALL p_alltoallv(parent_cell_owner_map, num_parent_owner_send, &
      &              parent_owner_send_displ, &
      &              parent_chunk_prep, num_parent_owner_recv, &
      &              parent_owner_recv_displ, p_comm_work)

    ! sort according to parent cell global index
    CALL quicksort(parent_chunk_prep(2,:), parent_chunk_prep(1,:))

    ! check that each parent cell global index occurs either 0 or 4 times
    IF (MOD(SIZE(parent_chunk_prep, 2), 4) /= 0) &
      CALL finish('divide_parent_cells','Incomplete parent cell encountered')
    DO i = 1, SIZE(parent_chunk_prep, 2), 4
      IF (ANY((parent_chunk_prep(1,i) /= parent_chunk_prep(1,i+1:i+3)) .OR. &
        &     (parent_chunk_prep(2,i) /= parent_chunk_prep(2,i+1:i+3)))) &
        CALL finish('divide_parent_cells','Incomplete parent cell encountered')
    END DO

    ! set up dist_array
    parent_chunk = p_parent_patch_pre%cells%local_chunk
    dist_cell_owner_p = dist_mult_array_new(dist_cell_owner_p_desc, &
      &                                     parent_chunk, p_comm_work)
    CALL dist_mult_array_local_ptr(dist_cell_owner_p, 1, parent_chunk_ptr)
    parent_chunk_ptr = -1
    DO i = 1, SIZE(parent_chunk_prep, 2), 4
      parent_chunk_ptr(parent_chunk_prep(2,i)) = parent_chunk_prep(1,i)
    END DO

    CALL dist_mult_array_expose(dist_cell_owner_p)

  END SUBROUTINE divide_parent_cells

  !-------------------------------------------------------------------------------------------------
  !>
  !! Divides a patch (in wrk_p_patch_pre) for parallelization.
  !!
  !! Parameters:
  !! cell_owner          owner PE numbers of the global cells
  !! n_boundary_rows     number of boundary rows to be added
  !! order_type_of_halos 1= a "normal" grid for prognostic computations is processed,
  !!                     0= a local parent grid is processed
  !!                       in the first case, a finer distinction between different halo
  !!                       cell/edge/vertex levels is made to optimize communication
  !!                     2=move all halos to the end (for ocean)
  !!
  !! On exit, the entries of wrk_p_patch are set.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Changed for usage for parent patch division, Rainer Johanni, Oct 2010
  !! Major rewrite to reduce memory consumption, Thomas Jahns and Moritz Hanke, Sep 2013

  SUBROUTINE divide_patch(wrk_p_patch, wrk_p_patch_pre, dist_cell_owner, &
    &                     n_boundary_rows, order_type_of_halos, my_proc)

    TYPE(t_patch), INTENT(INOUT) :: wrk_p_patch ! output patch, designated as INOUT because
                                                ! a few attributes are already set
    TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre ! out is required for
                                                        ! the distributed arrays

    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    INTEGER, INTENT(IN) :: n_boundary_rows
    INTEGER, INTENT(IN) :: order_type_of_halos
    INTEGER, INTENT(IN) :: my_proc

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_setup_subdivision::divide_patch'
    INTEGER :: i, j, jl, jb, je, jv, jg, jc, jc_p, jpg, jcg(4), jng, jeg, jvg

#ifdef __PGIC__
    TYPE(nb_flag_list_elem), ALLOCATABLE :: flag2_c_list(:), &
      flag2_v_list(:), flag2_e_list(:)
    INTEGER, ALLOCATABLE :: n2_ilev_c(:), n2_ilev_v(:), n2_ilev_e(:)
#else
    TYPE(nb_flag_list_elem) :: flag2_c_list(0:2*n_boundary_rows), &
         flag2_v_list(0:n_boundary_rows+1), flag2_e_list(0:2*n_boundary_rows+1)
    INTEGER :: n2_ilev_c(0:2*n_boundary_rows), n2_ilev_v(0:n_boundary_rows+1), &
         n2_ilev_e(0:2*n_boundary_rows+1)
#endif
    INTEGER, ALLOCATABLE :: owned_cells(:), owned_edges(:), owned_verts(:)

    IF (msg_level >= 10)  CALL message(routine, 'dividing patch')

    CALL compute_owned_cells()

#ifdef __PGIC__
    ALLOCATE(flag2_c_list(0:2*n_boundary_rows), &
      &      flag2_v_list(0:n_boundary_rows+1), &
      &      flag2_e_list(0:2*n_boundary_rows+1), &
      &      n2_ilev_c(0:2*n_boundary_rows), &
      &      n2_ilev_v(0:n_boundary_rows+1), &
      &      n2_ilev_e(0:2*n_boundary_rows+1))
#endif
    CALL compute_flag_lists(flag2_c_list, flag2_v_list, flag2_e_list, &
      &                     n2_ilev_c, n2_ilev_v, n2_ilev_e, &
      &                     n_boundary_rows, &
    ! removed arguments because NEC has some problems with them...
    !  &                     owned_edges, owned_verts, &
      &                     order_type_of_halos)

    !-----------------------------------------------------------------------------------------------
    ! Get the number of cells/edges/verts and other data for patch allocation
    !-----------------------------------------------------------------------------------------------

    CALL prepare_patch(wrk_p_patch_pre, wrk_p_patch, &
         SUM(n2_ilev_c(:)), SUM(n2_ilev_e(:)), SUM(n2_ilev_v(:)))

    !-----------------------------------------------------------------------------------------------
    ! Set the global ownership for cells, edges and verts (needed for boundary exchange).
    ! Please note that opposed to cells, the global owner for edges/verts is
    ! not unique at the boundaries of the divided patch.
    ! To minimize load imbalance at runtime, we set the owner to the PE with the highest processor number
    ! which is participating at this edge/vert if both PE numbers are even or odd, otherwise
    ! the lowest processor number is chosen
    !-----------------------------------------------------------------------------------------------
    CALL dist_dir_setup(wrk_p_patch%cells%decomp_info%owner_dist_dir, &
      &                 flag2_c_list(0)%idx(1:n2_ilev_c(0)), &
      &                 wrk_p_patch%n_patch_cells_g, p_comm_work, &
      &                 p_pe_work, p_n_work)
    CALL dist_dir_setup(wrk_p_patch%verts%decomp_info%owner_dist_dir, &
      &                 owned_verts, wrk_p_patch%n_patch_verts_g, &
      &                 p_comm_work, p_pe_work, p_n_work)
    CALL dist_dir_setup(wrk_p_patch%edges%decomp_info%owner_dist_dir, &
      &                 owned_edges, wrk_p_patch%n_patch_edges_g, &
      &                 p_comm_work, p_pe_work, p_n_work)

    !-----------------------------------------------------------------------------------------------
    ! Get the indices of local cells/edges/verts within patch and vice versa.

    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_cells, &
         n_patch_cve_g = wrk_p_patch_pre%n_patch_cells_g, &
         patch_id = wrk_p_patch_pre%id, &
         nblks = wrk_p_patch%nblks_c, &
         npromz = wrk_p_patch%npromz_c, &
         cell_type = wrk_p_patch%geometry_info%cell_type, &
         min_rlcve = min_rlcell, &
         min_rlcve_int = min_rlcell_int, &
         max_rlcve = max_rlcell, &
         max_ilev = 2*n_boundary_rows, &
         max_hw_cve = 2*max_hw, &
         flag2_list = flag2_c_list, &
         n2_ilev = n2_ilev_c, &
         decomp_info = wrk_p_patch%cells%decomp_info, &
         start_index = wrk_p_patch%cells%start_index, &
         start_block = wrk_p_patch%cells%start_block, &
         end_index = wrk_p_patch%cells%end_index, &
         end_block = wrk_p_patch%cells%end_block, &
         start_g = wrk_p_patch_pre%cells%start, &
         end_g = wrk_p_patch_pre%cells%end, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .TRUE., &
         refin_ctrl = wrk_p_patch_pre%cells%refin_ctrl, &
         refinement_predicate = refine_cells)
    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_edges, &
         n_patch_cve_g = wrk_p_patch_pre%n_patch_edges_g, &
         patch_id = wrk_p_patch_pre%id, &
         nblks = wrk_p_patch%nblks_e, &
         npromz = wrk_p_patch%npromz_e, &
         cell_type = wrk_p_patch%geometry_info%cell_type, &
         min_rlcve = min_rledge, &
         min_rlcve_int = min_rledge_int, &
         max_rlcve = max_rledge, &
         max_ilev = 2*n_boundary_rows+1, &
         max_hw_cve = 2*max_hw + 1, &
         flag2_list = flag2_e_list, &
         n2_ilev = n2_ilev_e, &
         decomp_info = wrk_p_patch%edges%decomp_info, &
         start_index = wrk_p_patch%edges%start_index, &
         start_block = wrk_p_patch%edges%start_block, &
         end_index = wrk_p_patch%edges%end_index, &
         end_block = wrk_p_patch%edges%end_block, &
         start_g = wrk_p_patch_pre%edges%start, &
         end_g = wrk_p_patch_pre%edges%end, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .FALSE., &
         refin_ctrl = wrk_p_patch_pre%edges%refin_ctrl, &
         refinement_predicate = refine_edges)
    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_verts, &
         n_patch_cve_g = wrk_p_patch_pre%n_patch_verts_g, &
         patch_id = wrk_p_patch_pre%id, &
         nblks = wrk_p_patch%nblks_v, &
         npromz = wrk_p_patch%npromz_v, &
         cell_type = wrk_p_patch%geometry_info%cell_type, &
         min_rlcve = min_rlvert, &
         min_rlcve_int = min_rlvert_int, &
         max_rlcve = max_rlvert, &
         max_ilev = n_boundary_rows+1, &
         max_hw_cve = max_hw + 1, &
         flag2_list = flag2_v_list, &
         n2_ilev = n2_ilev_v, &
         decomp_info = wrk_p_patch%verts%decomp_info, &
         start_index = wrk_p_patch%verts%start_index, &
         start_block = wrk_p_patch%verts%start_block, &
         end_index = wrk_p_patch%verts%end_index, &
         end_block = wrk_p_patch%verts%end_block, &
         start_g = wrk_p_patch_pre%verts%start, &
         end_g = wrk_p_patch_pre%verts%end, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .FALSE., &
         refin_ctrl = wrk_p_patch_pre%verts%refin_ctrl, &
         refinement_predicate = refine_verts)

    ! Sanity checks: are there still elements of the index lists filled with dummy values?
    ! Note: the use of write(0,*) instead of CALL message is intended here
    !       because the debug output is wanted for all PEs where an error occurs
    IF (ANY(wrk_p_patch%cells%start_index(:)<0) .OR. ANY(wrk_p_patch%cells%start_block(:)<0) .OR.&
        ANY(wrk_p_patch%cells%end_index(:)  <0) .OR. ANY(wrk_p_patch%cells%end_block(:)  <0)) THEN
      DO i = min_rlcell, max_rlcell
        WRITE(0,'(a,2i5,2i4,4i7)') 'cells',my_proc,wrk_p_patch%id,i,     &
          wrk_p_patch%cells%start_block(i),wrk_p_patch%cells%start_index(i), &
          wrk_p_patch%cells%end_block(i),  wrk_p_patch%cells%end_index(i)
      ENDDO
      CALL finish('divide_patch','Error in cell start/end indices')
    ENDIF

    IF (ANY(wrk_p_patch%edges%start_index(:)<0) .OR. ANY(wrk_p_patch%edges%start_block(:)<0) .OR.&
        ANY(wrk_p_patch%edges%end_index(:)  <0) .OR. ANY(wrk_p_patch%edges%end_block(:)  <0)) THEN
      DO i = min_rledge, max_rledge
        write(0,'(a,2i5,2i4,4i7)') 'edges',my_proc,wrk_p_patch%id,i,     &
          wrk_p_patch%edges%start_block(i),wrk_p_patch%edges%start_index(i), &
          wrk_p_patch%edges%end_block(i),  wrk_p_patch%edges%end_index(i)
      ENDDO
      CALL finish('divide_patch','Error in edge start/end indices')
    ENDIF

    IF (ANY(wrk_p_patch%verts%start_index(:)<0) .OR. ANY(wrk_p_patch%verts%start_block(:)<0) .OR.&
        ANY(wrk_p_patch%verts%end_index(:)  <0) .OR. ANY(wrk_p_patch%verts%end_block(:)  <0)) THEN
      DO i = min_rlvert, max_rlvert
        write(0,'(a,2i5,2i4,4i7)') 'verts',my_proc,wrk_p_patch%id,i,     &
          wrk_p_patch%verts%start_block(i),wrk_p_patch%verts%start_index(i), &
          wrk_p_patch%verts%end_block(i),  wrk_p_patch%verts%end_index(i)
      ENDDO
    CALL finish('divide_patch','Error in vertex start/end indices')
    ENDIF

    ! Fill 'old' two-dimensional index fields until removing them is completed
    DO j = 1, MAX(1,wrk_p_patch%n_childdom)
      wrk_p_patch%cells%start_blk(:,j) = wrk_p_patch%cells%start_block(:)
      wrk_p_patch%cells%start_idx(:,j) = wrk_p_patch%cells%start_index(:)
      wrk_p_patch%cells%end_blk(:,j)   = wrk_p_patch%cells%end_block(:)
      wrk_p_patch%cells%end_idx(:,j)   = wrk_p_patch%cells%end_index(:)

      wrk_p_patch%edges%start_blk(:,j) = wrk_p_patch%edges%start_block(:)
      wrk_p_patch%edges%start_idx(:,j) = wrk_p_patch%edges%start_index(:)
      wrk_p_patch%edges%end_blk(:,j)   = wrk_p_patch%edges%end_block(:)
      wrk_p_patch%edges%end_idx(:,j)   = wrk_p_patch%edges%end_index(:)

      wrk_p_patch%verts%start_blk(:,j) = wrk_p_patch%verts%start_block(:)
      wrk_p_patch%verts%start_idx(:,j) = wrk_p_patch%verts%start_index(:)
      wrk_p_patch%verts%end_blk(:,j)   = wrk_p_patch%verts%end_block(:)
      wrk_p_patch%verts%end_idx(:,j)   = wrk_p_patch%verts%end_index(:)
    ENDDO

    !-----------------------------------------------------------------------------------------------
    ! Set arrays of divided patch
    !-----------------------------------------------------------------------------------------------

    ! decomp_domain is -1 for invalid locations (at the end of the last strip)

    wrk_p_patch%cells%decomp_info%decomp_domain = -1
    wrk_p_patch%edges%decomp_info%decomp_domain = -1
    wrk_p_patch%verts%decomp_info%decomp_domain = -1

    !---------------------------------------------------------------------------------------

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      DO i=1,wrk_p_patch%geometry_info%cell_type
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_neighbor,(/jg,i/),jng)
!CDIR IEXPAND
        CALL get_local_idx(wrk_p_patch%cells%decomp_info, jng, jc)

        wrk_p_patch%cells%neighbor_idx(jl,jb,i) = idx_no(jc)
        wrk_p_patch%cells%neighbor_blk(jl,jb,i) = blk_no(jc)

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
             (/jg,i/), jeg)
!CDIR IEXPAND
        CALL get_local_idx(wrk_p_patch%edges%decomp_info, jeg, je)

        wrk_p_patch%cells%edge_idx(jl,jb,i) = idx_no(je)
        wrk_p_patch%cells%edge_blk(jl,jb,i) = blk_no(je)

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
             (/jg,i/), jvg)
!CDIR IEXPAND
        CALL get_local_idx(wrk_p_patch%verts%decomp_info, jvg, jv)

        wrk_p_patch%cells%vertex_idx(jl,jb,i) = idx_no(jv)
        wrk_p_patch%cells%vertex_blk(jl,jb,i) = blk_no(jv)

      ENDDO
    ENDDO

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      ! parent and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_parent, (/jg/), jc_p)

      wrk_p_patch%cells%parent_glb_idx(jl,jb)  = idx_no(jc_p)
      wrk_p_patch%cells%parent_glb_blk(jl,jb)  = blk_no(jc_p)
      DO i = 1, 4
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_child, (/jg, i/), jc)
        wrk_p_patch%cells%child_idx(jl,jb,i) = idx_no(jc)
        wrk_p_patch%cells%child_blk(jl,jb,i) = blk_no(jc)
      END DO

      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, (/jg/), &
        &                      wrk_p_patch%cells%num_edges(jl,jb))
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%center, 1, (/jg/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lat)
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%center, 2, (/jg/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lon)
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%refin_ctrl, 1, (/jg/), &
        &                      wrk_p_patch%cells%refin_ctrl(jl,jb))
    ENDDO

    DO j = 0, 2 * n_boundary_rows
      DO i = 1, n2_ilev_c(j)
        jg = flag2_c_list(j)%idx(i)
        jc = get_local_index(wrk_p_patch%cells%decomp_info%glb2loc_index, jg)
        jb = blk_no(jc)
        jl = idx_no(jc)
        wrk_p_patch%cells%decomp_info%decomp_domain(jl,jb) = j
      END DO
    END DO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%edges%decomp_info%glb_index(j)

      ! parent and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

      CALL dist_mult_array_get(wrk_p_patch_pre%edges%parent, 1, (/jg/), jpg)
      DO i = 1, 4
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%child, 1, (/jg, i/), jcg(i))
      END DO

      wrk_p_patch%edges%parent_glb_idx(jl,jb)    = idx_no(jpg)
      wrk_p_patch%edges%parent_glb_blk(jl,jb)    = blk_no(jpg)
      wrk_p_patch%edges%child_idx(jl,jb,1:4) = idx_no(jcg)
      wrk_p_patch%edges%child_blk(jl,jb,1:4) = blk_no(jcg)

      CALL dist_mult_array_get(wrk_p_patch_pre%edges%refin_ctrl, 1, (/jg/), &
        &                      wrk_p_patch%edges%refin_ctrl(jl,jb))
    ENDDO

    DO j = 0, 2 * n_boundary_rows + 1
      DO i = 1, n2_ilev_e(j)
        jg = flag2_e_list(j)%idx(i)
        je = get_local_index(wrk_p_patch%edges%decomp_info%glb2loc_index, jg)
        jb = blk_no(je)
        jl = idx_no(je)
        wrk_p_patch%edges%decomp_info%decomp_domain(jl,jb) = j
      END DO
    END DO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%verts%decomp_info%glb_index(j)

      CALL dist_mult_array_get(wrk_p_patch_pre%verts%vertex, 1, (/jg/), &
        &                      wrk_p_patch%verts%vertex(jl,jb)%lat)
      CALL dist_mult_array_get(wrk_p_patch_pre%verts%vertex, 2, (/jg/), &
        &                      wrk_p_patch%verts%vertex(jl,jb)%lon)
      CALL dist_mult_array_get(wrk_p_patch_pre%verts%refin_ctrl, 1, (/jg/), &
        &                      wrk_p_patch%verts%refin_ctrl(jl,jb))
    ENDDO

    DO j = 0, n_boundary_rows + 1
      DO i = 1, n2_ilev_v(j)
        jg = flag2_v_list(j)%idx(i)
        jv = get_local_index(wrk_p_patch%verts%decomp_info%glb2loc_index, jg)
        jb = blk_no(jv)
        jl = idx_no(jv)
        wrk_p_patch%verts%decomp_info%decomp_domain(jl,jb) = j
      END DO
    END DO

    DEALLOCATE(owned_cells)

  CONTAINS

    SUBROUTINE compute_owned_cells()

      INTEGER, POINTER :: local_chunk(:)
      INTEGER :: local_chunk_first, local_chunk_last, local_chunk_size, i
      INTEGER, ALLOCATABLE :: n_cells_owned_per_proc(:)
      ! n_cells_from_proc(i) = number of cell ownerships transferred
      ! from process i
      INTEGER, ALLOCATABLE :: n_cells_from_proc(:)
      INTEGER, ALLOCATABLE :: glb_index_per_proc(:)
      INTEGER, ALLOCATABLE :: sorted_local_chunk(:)
      INTEGER, ALLOCATABLE :: send_displs(:)
      INTEGER, ALLOCATABLE :: recv_displs(:)

      ! get local chunk
      CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_chunk)
      local_chunk_first = LBOUND(local_chunk, 1)
      local_chunk_size = SIZE(local_chunk)
      local_chunk_last = local_chunk_first + local_chunk_size - 1

      ! for parent patches most of the entries are -1
      ALLOCATE(n_cells_owned_per_proc(-1:p_n_work-1), &
        &      n_cells_from_proc(0:p_n_work-1), &
        &      glb_index_per_proc(local_chunk_first:local_chunk_last), &
        &      sorted_local_chunk(local_chunk_first:local_chunk_last))

      ! compute n_cells_owned_per_proc
      n_cells_owned_per_proc = 0
      DO i = local_chunk_first, local_chunk_last
        n_cells_owned_per_proc(local_chunk(i)) = &
          n_cells_owned_per_proc(local_chunk(i)) + 1

        ! initialise glb_index_per_proc
        glb_index_per_proc(i) = i
      END DO

      ! mpi_alltoall n_cells_owned_per_proc
      CALL p_alltoall(n_cells_owned_per_proc(0:), n_cells_from_proc, &
        &             p_comm_work)

      ! compute glb_index_per_proc sort by owner (see SUBROUTINE quicksort)
      sorted_local_chunk = local_chunk
      CALL quicksort(sorted_local_chunk, glb_index_per_proc)

      ! mpi_alltoallv glb_index_per_proc -> owned_cells
      ALLOCATE(send_displs(0:p_n_work-1), recv_displs(0:p_n_work-1))
      send_displs(0) = n_cells_owned_per_proc(-1)
      recv_displs(0) = 0
      DO i = 1, p_n_work-1
        send_displs(i) = send_displs(i-1) + n_cells_owned_per_proc(i-1)
        recv_displs(i) = recv_displs(i-1) + n_cells_from_proc(i-1)
      END DO
      ALLOCATE(owned_cells(recv_displs(p_n_work-1) + &
        &                  n_cells_from_proc(p_n_work-1)))
      CALL p_alltoallv(glb_index_per_proc, n_cells_owned_per_proc(0:), &
        &              send_displs, owned_cells, n_cells_from_proc, &
        &              recv_displs, p_comm_work)

      ! sort owned_cells
      CALL quicksort(owned_cells)

    END SUBROUTINE compute_owned_cells

    SUBROUTINE prepare_patch(wrk_p_patch_pre, wrk_p_patch, &
         n_patch_cells, n_patch_edges, n_patch_verts)
      !> output patch, designated as INOUT because
      !! a few attributes are already set
      TYPE(t_patch), INTENT(inout) :: wrk_p_patch
      TYPE(t_pre_patch), INTENT(in) :: wrk_p_patch_pre
      INTEGER, INTENT(in) :: n_patch_cells, n_patch_edges, n_patch_verts

      wrk_p_patch%n_patch_cells = n_patch_cells
      wrk_p_patch%n_patch_edges = n_patch_edges
      wrk_p_patch%n_patch_verts = n_patch_verts

      ! save the number of cells/edges/verts of the patch
      wrk_p_patch%n_patch_cells_g = wrk_p_patch_pre%n_patch_cells_g
      wrk_p_patch%n_patch_edges_g = wrk_p_patch_pre%n_patch_edges_g
      wrk_p_patch%n_patch_verts_g = wrk_p_patch_pre%n_patch_verts_g
      !
      ! calculate and save values for the blocking, these are needed for patch allocation.
      ! NB: Avoid the case nblks=0 for empty patches, this might cause troubles
      ! if a empty patch is used somewhere (and npromz gets wrong in the formulas below).
      !
      ! ... for the cells
      wrk_p_patch%nblks_c       = blk_no(wrk_p_patch%n_patch_cells)
      wrk_p_patch%npromz_c      = wrk_p_patch%n_patch_cells - (wrk_p_patch%nblks_c - 1)*nproma
      wrk_p_patch%alloc_cell_blocks = wrk_p_patch%nblks_c
      IF (use_dummy_cell_closure) THEN
        IF (wrk_p_patch%npromz_c == nproma) &
             wrk_p_patch%alloc_cell_blocks = wrk_p_patch%nblks_c + 1
      ENDIF

      ! ... for the edges
      wrk_p_patch%nblks_e       = blk_no(wrk_p_patch%n_patch_edges)
      wrk_p_patch%npromz_e      = wrk_p_patch%n_patch_edges - &
        &                         (wrk_p_patch%nblks_e - 1) * nproma

      ! ... for the vertices
      wrk_p_patch%nblks_v       = blk_no(wrk_p_patch%n_patch_verts)
      wrk_p_patch%npromz_v      = wrk_p_patch%n_patch_verts - (wrk_p_patch%nblks_v - 1)*nproma

      ! Also needed for patch allocation
      wrk_p_patch%max_childdom  = wrk_p_patch_pre%max_childdom

      ! Set other scalar members of patch here too ..
      wrk_p_patch%grid_filename      = wrk_p_patch_pre%grid_filename
      wrk_p_patch%grid_filename_grfinfo = wrk_p_patch_pre%grid_filename_grfinfo
      wrk_p_patch%level              = wrk_p_patch_pre%level
      wrk_p_patch%id                 = wrk_p_patch_pre%id
      wrk_p_patch%cells%max_connectivity = wrk_p_patch_pre%cells%max_connectivity
      wrk_p_patch%verts%max_connectivity = wrk_p_patch_pre%verts%max_connectivity
      wrk_p_patch%parent_id          = wrk_p_patch_pre%parent_id
      wrk_p_patch%parent_child_index = wrk_p_patch_pre%parent_child_index
      wrk_p_patch%child_id(:)        = wrk_p_patch_pre%child_id(:)
      wrk_p_patch%child_id_list(:)   = wrk_p_patch_pre%child_id_list(:)
      wrk_p_patch%n_childdom         = wrk_p_patch_pre%n_childdom
      wrk_p_patch%n_chd_total        = wrk_p_patch_pre%n_chd_total
      wrk_p_patch%nlev               = wrk_p_patch_pre%nlev
      wrk_p_patch%nlevp1             = wrk_p_patch_pre%nlevp1
      wrk_p_patch%nshift             = wrk_p_patch_pre%nshift
      wrk_p_patch%nshift_total       = wrk_p_patch_pre%nshift_total
      wrk_p_patch%nshift_child       = wrk_p_patch_pre%nshift_child
      wrk_p_patch%grid_uuid          = wrk_p_patch_pre%grid_uuid

      !-----------------------------------------------------------------------------------------------
      ! Allocate the required data arrays in patch
      !-----------------------------------------------------------------------------------------------

      CALL allocate_basic_patch(wrk_p_patch)
      CALL allocate_remaining_patch(wrk_p_patch,2) ! 2 = only those needed for parallelization control


    END SUBROUTINE prepare_patch

    ! this routine assumes that both lists are sorted in ascending order
    SUBROUTINE remove_entries_from_ref_list(list, n, ref_list)
      INTEGER, INTENT(inout) :: n
      INTEGER, INTENT(inout) :: list(n)
      INTEGER, INTENT(in) :: ref_list(:)

      INTEGER :: i, j, k, ref_list_size

      ref_list_size = SIZE(ref_list(:))

      IF (n > 0 .AND. ref_list_size > 0) THEN
        ! initialise j
        DO j = 1, ref_list_size
          IF (ref_list(j) >= list(1)) &
            EXIT
        END DO
        IF (j <= ref_list_size) THEN
          k = 0
          DO i = 1, n-1
            IF (list(i) /= ref_list(j)) THEN
              k = k + 1
              list(k) = list(i)
            END IF
            DO j = j, ref_list_size
              IF (ref_list(j) >= list(i+1)) &
                EXIT
            END DO
            IF (j > ref_list_size) EXIT
          END DO
          IF (j > ref_list_size) THEN
            DO i = i+1, n
              k = k + 1
              list(k) = list(i)
            END DO
          ELSE IF (list(i) /= ref_list(j)) THEN
            k = k + 1
            list(k) = list(i)
          END IF
          n = k
        END IF
      END IF
    END SUBROUTINE remove_entries_from_ref_list

    ! this routine assumes that the list is sorted in ascending order
    SUBROUTINE remove_duplicated_entries(list, n)
      INTEGER, INTENT(inout) :: n
      INTEGER, INTENT(inout) :: list(n)

      INTEGER :: i, j

      IF (n == 0) THEN
        RETURN
      END IF

      i = 1

      DO j = 2, n
        IF (list(j-1) /= list(j)) THEN
          i = i + 1
          list(i) = list(j)
        END IF
      END DO

      n = i

    END SUBROUTINE remove_duplicated_entries

    SUBROUTINE compute_flag_lists(flag2_c_list, flag2_v_list, flag2_e_list, &
      &                           n2_ilev_c, n2_ilev_v, n2_ilev_e, &
      &                           n_boundary_rows, &
    ! removed arguments because NEC has some problems with them...
    !  &                           owned_edges, owned_verts, &
      &                           order_type_of_halos)
      INTEGER, INTENT(IN) :: n_boundary_rows
      TYPE(nb_flag_list_elem), INTENT(OUT) :: &
        flag2_c_list(0:2*n_boundary_rows), flag2_v_list(0:n_boundary_rows+1), &
        flag2_e_list(0:2*n_boundary_rows+1)
      INTEGER, INTENT(OUT) :: n2_ilev_c(0:2*n_boundary_rows), n2_ilev_v(0:n_boundary_rows+1), &
        &  n2_ilev_e(0:2*n_boundary_rows+1)
      ! removed arguments because NEC has some problems with them...
      ! INTEGER, ALLOCATABLE, INTENT(OUT) :: owned_edges(:), owned_verts(:)
      INTEGER, INTENT(IN) :: order_type_of_halos

      INTEGER :: n, i, ic, jv, jv_, je, ilev, jc, jc_, k, owner_jc, owner_jc_
      LOGICAL, ALLOCATABLE :: pack_mask(:)
      INTEGER, ALLOCATABLE :: temp_cells(:), temp_vertices(:), &
                              temp_vertices_owner(:), temp_edges(:), &
                              temp_edges_owner(:), inner_edges(:), &
                              edge_cells(:)
      INTEGER :: n_inner_edges, n_temp_edges, n_temp_cells, n_temp_vertices, &
        &        temp_num_edges

      !---------------------------------------------------------------------
      ! flag_c_list(-1)%idx empty dummy list
      ! flag_c_list(0)%idx  all cells jg where cell_owner(jg) == my_proc
      ! flag_c_list(j)%idx j in 1..n_boundary_cells all cells bordering
      !                    on cells in flag_c_list(j - 1)
      ! flag_e_list(-1)%idx empty dummy list
      ! flag_e_list(j)%idx all edges where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the edge
      ! flag_v_list(-1)%idx empty dummy list
      ! flag_v_list(j)%idx all vertices where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the vertex
      !------------------------------------------------------------------------
      ! The purpose of the second set of flags is to control moving the
      ! halo points to the end of the index vector for nearly all cells
      ! even if order_type_of_halos = 1
      !------------------------------------------------------------------------
      ! flag2_e_list(0)%idx all edges which are owned by the local process
      ! flag2_e_list(1)%idx all edges adjacent to cells in flag_c_list(0)
      !                     but owned by other tasks
      ! flag2_e_list(2*i)%idx where i > 0, edges in idx are ONLY adjacent
      !                       to cells in flag_c_list(i)%idx
      ! flag2_e_list(2*i+1)%idx where i > 0, edges adjacent to one cell in
      !                         flag_c_list(i) and either
      !                            one cell in flag_c_list(i+1)
      !                         or edge is an outer edge,
      !                            i.e. adjacent to only one cell
      ! flag2_c_list(0)%idx == flag_c_list(0)%idx
      ! flag2_c_list(2*j-1)%idx where j > 0, contains all cells from
      !                         flag_c_list(j) which are
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_c_list(2*j)%idx   where j > 0, contains all cells from
      !                         flag_c_list(j) which are NOT
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_v_list(0)%idx contains all vertices from flag_v_list(0)%idx
      !                         owned by the current task
      ! flag2_v_list(1)%idx contains all vertices from flag_v_list(0)%idx
      !                     NOT owned by the current task
      ! flag2_v_list(j)%idx for j > 1 == flag_v_list(j-1)%idx
      !------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      ! compute flag2 arrays
      !--------------------------------------------------------------------------

      ALLOCATE(pack_mask(0))

      ! collect cells of level 0
      n2_ilev_c(0) = SIZE(owned_cells(:))

      IF (n2_ilev_c(0) == wrk_p_patch_pre%n_patch_cells_g) THEN

        CALL compute_flag_lists_short(flag2_c_list, flag2_v_list, flag2_e_list, &
          &                           n2_ilev_c, n2_ilev_v, n2_ilev_e, &
          &                           n_boundary_rows &
        ! removed arguments because NEC has some problems with them...
        !  &                           ,owned_edges, owned_verts &
                                      )

        RETURN
      END IF

      ALLOCATE(flag2_c_list(0)%idx(n2_ilev_c(0)), &
        &      flag2_c_list(0)%owner(n2_ilev_c(0)))

      flag2_c_list(0)%idx(:) = owned_cells(:)
      flag2_c_list(0)%owner(:) = my_proc

      n_inner_edges = 0
      n_temp_edges = 0
      n_temp_vertices = 0
      ALLOCATE(inner_edges(n2_ilev_c(0) * &
                           wrk_p_patch_pre%cells%max_connectivity), &
               temp_edges(n2_ilev_c(0) * &
                          wrk_p_patch_pre%cells%max_connectivity), &
               edge_cells(n2_ilev_c(0) * &
                          wrk_p_patch_pre%cells%max_connectivity), &
               temp_vertices(n2_ilev_c(0) * &
                             wrk_p_patch_pre%cells%max_connectivity))

      ! collect inner and outer edges and vertices adjacent to cells of level 0
      DO ic = 1, n2_ilev_c(0)

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
          &                      (/flag2_c_list(0)%idx(ic)/), temp_num_edges)
        DO i = 1, temp_num_edges

          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
            &                      (/flag2_c_list(0)%idx(ic),i/), je)

          CALL dist_mult_array_get(wrk_p_patch_pre%edges%cell, 1, (/je,1/), jc)
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%cell, 1, (/je,2/), jc_)

          IF (jc <= 0 .OR. jc_ <= 0 .OR. jc == jc_) THEN
            n_inner_edges = n_inner_edges + 1
            inner_edges(n_inner_edges) = je
          ELSE

            CALL dist_mult_array_get(dist_cell_owner, 1, (/jc/), owner_jc)
            CALL dist_mult_array_get(dist_cell_owner, 1, (/jc_/), owner_jc_)

            IF (owner_jc == owner_jc_) THEN
              n_inner_edges = n_inner_edges + 1
              inner_edges(n_inner_edges) = je
            ELSE
              n_temp_edges = n_temp_edges + 1
              temp_edges(n_temp_edges) = je
              edge_cells(n_temp_edges) = flag2_c_list(0)%idx(ic)
            END IF
          END IF

          n_temp_vertices = n_temp_vertices + 1
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
            &                      (/flag2_c_list(0)%idx(ic),i/), &
            &                      temp_vertices(n_temp_vertices))
        END DO
      END DO

      ! remove duplicated inner edges
      CALL insertion_sort(inner_edges(1:n_inner_edges))
      CALL remove_duplicated_entries(inner_edges(1:n_inner_edges), n_inner_edges)
      ! remove duplicated vertices
      CALL insertion_sort(temp_vertices(1:n_temp_vertices))
      CALL remove_duplicated_entries(temp_vertices(1:n_temp_vertices), &
                                     n_temp_vertices)

      ! MoHa Note: vertices that are only adjacent to inner edges do not
      ! necessarily also have to be inner vertices...there are special cases...

      IF (SIZE(pack_mask(:)) < n_temp_vertices) THEN
        DEALLOCATE(pack_mask)
        ALLOCATE(pack_mask(n_temp_vertices))
      END IF

      ALLOCATE(temp_vertices_owner(n_temp_vertices))
      CALL compute_vertex_owner(temp_vertices(1:n_temp_vertices), &
        &                       temp_vertices_owner(1:n_temp_vertices))
      pack_mask(1:n_temp_vertices) = temp_vertices_owner(1:n_temp_vertices) &
        &                            == my_proc

      n = COUNT(pack_mask(1:n_temp_vertices))

      ALLOCATE(owned_verts(n))
      owned_verts(:) = PACK(temp_vertices(1:n_temp_vertices), &
                            pack_mask(1:n_temp_vertices))

      ! generate flag2_v_list(0) and flag2_v_list(1)
      IF (order_type_of_halos == 0 .OR. order_type_of_halos == 2 ) THEN
        ALLOCATE(flag2_v_list(0)%idx(n_temp_vertices), &
          &      flag2_v_list(0)%owner(n_temp_vertices), &
          &      flag2_v_list(1)%idx(0), &
          &      flag2_v_list(1)%owner(0))
        n2_ilev_v(0) = n_temp_vertices
        flag2_v_list(0)%idx(:) = temp_vertices(1:n_temp_vertices)
        flag2_v_list(0)%owner(:) = temp_vertices_owner(1:n_temp_vertices)
        n2_ilev_v(1) = 0
      ELSE
        ALLOCATE(flag2_v_list(0)%idx(n), flag2_v_list(0)%owner(n), &
                 flag2_v_list(1)%idx(n_temp_vertices - n), &
                 flag2_v_list(1)%owner(n_temp_vertices - n))
        n2_ilev_v(0) = n
        flag2_v_list(0)%idx(:) = owned_verts(:)
        flag2_v_list(0)%owner(:) = &
          PACK(temp_vertices_owner(1:n_temp_vertices), &
            &  pack_mask(1:n_temp_vertices))
        n2_ilev_v(1) = n_temp_vertices - n
        flag2_v_list(1)%idx(:) = PACK(temp_vertices(1:n_temp_vertices), &
                                      .NOT. pack_mask(1:n_temp_vertices))
        flag2_v_list(1)%owner(:) = &
          PACK(temp_vertices_owner(1:n_temp_vertices), &
            &  .NOT. pack_mask(1:n_temp_vertices))
      END IF

      CALL quicksort(temp_edges(1:n_temp_edges), edge_cells(1:n_temp_edges))
      IF (SIZE(pack_mask(:)) < n_temp_edges) THEN
        DEALLOCATE(pack_mask)
        ALLOCATE(pack_mask(n_temp_edges))
      END IF
      ALLOCATE(temp_edges_owner(n_temp_edges))
      CALL compute_edge_owner(temp_edges(1:n_temp_edges), &
        &                     temp_edges_owner(1:n_temp_edges))
      pack_mask(1:n_temp_edges) = &
        .NOT. temp_edges_owner(1:n_temp_edges) == my_proc
      n = COUNT(pack_mask(1:n_temp_edges))
      ALLOCATE(owned_edges(n_inner_edges + n_temp_edges - n))
      owned_edges(1:n_inner_edges) = inner_edges(1:n_inner_edges)
      owned_edges(n_inner_edges+1:) = PACK(temp_edges(1:n_temp_edges), &
                                             .NOT. pack_mask(1:n_temp_edges))
      CALL insertion_sort(owned_edges(:))

      ! generate flag2_e_list(0) and flag2_e_list(1)
  !    IF (n_boundary_rows > 0 .AND. order_type_of_halos /= 0) THEN
      IF (n_boundary_rows > 0 .AND. order_type_of_halos == 1) THEN

        ALLOCATE(flag2_e_list(0)%idx(n_inner_edges + n_temp_edges - n), &
          &      flag2_e_list(0)%owner(n_inner_edges + n_temp_edges - n), &
          &      flag2_e_list(1)%idx(n), flag2_e_list(1)%owner(n))
        n2_ilev_e(0) = n_inner_edges + n_temp_edges - n
        n2_ilev_e(1) = n
        flag2_e_list(0)%idx(:) = owned_edges(:)
        flag2_e_list(0)%owner(:) = my_proc
        flag2_e_list(1)%idx(:) = PACK(temp_edges(1:n_temp_edges), &
                                        pack_mask(1:n_temp_edges))
        flag2_e_list(1)%owner(:) = PACK(temp_edges_owner(1:n_temp_edges), &
                                        pack_mask(1:n_temp_edges))
      ELSE

        ALLOCATE(flag2_e_list(0)%idx(n_inner_edges + n_temp_edges), &
          &      flag2_e_list(0)%owner(n_inner_edges + n_temp_edges), &
          &      flag2_e_list(1)%idx(0), flag2_e_list(1)%owner(0))
        n2_ilev_e(0) = n_inner_edges + n_temp_edges
        n2_ilev_e(1) = 0
        flag2_e_list(0)%idx(1:n_inner_edges) = inner_edges(1:n_inner_edges)
        flag2_e_list(0)%owner(1:n_inner_edges) = my_proc
        flag2_e_list(0)%idx(n_inner_edges+1:) = temp_edges(1:n_temp_edges)
        flag2_e_list(0)%owner(n_inner_edges+1:) = temp_edges_owner(1:n_temp_edges)
        CALL quicksort(flag2_e_list(0)%idx(:), flag2_e_list(0)%owner(:))
      END IF

      DEALLOCATE(inner_edges)
      ALLOCATE(temp_cells(n_temp_edges))

      DO ilev = 1, n_boundary_rows

        ! collect cells of level 2*ilev-1
        n_temp_cells = 0
        IF (SIZE(temp_cells(:)) < n_temp_edges) THEN

          DEALLOCATE(temp_cells)
          ALLOCATE(temp_cells(n_temp_edges))
        END IF
        DO i = 1, n_temp_edges

          je = temp_edges(i)
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%cell,1, (/je,1/), jc)

          IF (jc == edge_cells(i)) THEN
            CALL dist_mult_array_get(wrk_p_patch_pre%edges%cell,1, (/je,2/),jc)
          END IF

          IF (jc > 0 .AND. jc /= edge_cells(i)) THEN
            n_temp_cells = n_temp_cells + 1
            temp_cells(n_temp_cells) = jc
          END IF
        END DO
        !remove duplicated entries
        CALL insertion_sort(temp_cells(1:n_temp_cells))
        CALL remove_duplicated_entries(temp_cells(1:n_temp_cells), n_temp_cells)

        ! store cells of level 2*ilev-1
        n2_ilev_c(2*ilev-1) = n_temp_cells
        ALLOCATE(flag2_c_list(2*ilev-1)%idx(n_temp_cells), &
          &      flag2_c_list(2*ilev-1)%owner(n_temp_cells))
        flag2_c_list(2*ilev-1)%idx(:) = temp_cells(1:n_temp_cells)
        DO i = 1, n_temp_cells
          CALL dist_mult_array_get(dist_cell_owner, 1, (/temp_cells(i)/), &
            &                      flag2_c_list(2*ilev-1)%owner(i))
        END DO

        ! collect all cells adjacent to the outer vertices
        IF (SIZE(temp_cells(:)) < n_temp_vertices * &
          & wrk_p_patch_pre%verts%max_connectivity) THEN

          DEALLOCATE(temp_cells)
          ALLOCATE(temp_cells(n_temp_vertices * &
            &      wrk_p_patch_pre%verts%max_connectivity))
        END IF
        n_temp_cells = 0
        DO jv = 1, n_temp_vertices
          jv_ = temp_vertices(jv)
          CALL dist_mult_array_get(wrk_p_patch_pre%verts%num_edges, 1, &
            &                      (/jv_/), temp_num_edges)
          DO i = 1, temp_num_edges
            CALL dist_mult_array_get(wrk_p_patch_pre%verts%cell,1,(/jv_,i/),jc)
            IF (jc > 0) THEN
              n_temp_cells = n_temp_cells + 1
              temp_cells(n_temp_cells) = jc
            END IF
          END DO
        END DO
        CALL insertion_sort(temp_cells(1:n_temp_cells))
        CALL remove_duplicated_entries(temp_cells(1:n_temp_cells), n_temp_cells)
        ! remove cells that are on level 2*ilev-2 and level 2*ilev-1
        DO k = -2, -1
          CALL remove_entries_from_ref_list(temp_cells(1:n_temp_cells), &
            n_temp_cells, &
            flag2_c_list(2 * ilev + k)%idx(1:n2_ilev_c(2 * ilev + k)))
        END DO
        ! store cells of level 2*ilev
        n2_ilev_c(2*ilev) = n_temp_cells
        ALLOCATE(flag2_c_list(2*ilev)%idx(n_temp_cells), &
          &      flag2_c_list(2*ilev)%owner(n_temp_cells))
        flag2_c_list(2*ilev)%idx(:) = temp_cells(1:n_temp_cells)
        DO i = 1, n_temp_cells
          CALL dist_mult_array_get(dist_cell_owner, 1, (/temp_cells(i)/), &
            &                      flag2_c_list(2*ilev)%owner(i))
        END DO

        ! get all edges of cells of level 2*ilev and level 2*ilev-1
        IF (SIZE(temp_edges(:)) < (n2_ilev_c(2*ilev) + &
          &                        n2_ilev_c(2*ilev - 1)) * &
          &                       wrk_p_patch_pre%cells%max_connectivity) THEN

          DEALLOCATE(temp_edges)
          ALLOCATE(temp_edges((n2_ilev_c(2*ilev) + &
            &                  n2_ilev_c(2*ilev - 1)) * &
            &                 wrk_p_patch_pre%cells%max_connectivity))
        END IF
        IF (SIZE(edge_cells(:)) < (n2_ilev_c(2*ilev) + &
          &                        n2_ilev_c(2*ilev - 1)) * &
          &                       wrk_p_patch_pre%cells%max_connectivity) THEN

          DEALLOCATE(edge_cells)
          ALLOCATE(edge_cells((n2_ilev_c(2*ilev) + &
            &                  n2_ilev_c(2*ilev - 1)) * &
            &                 wrk_p_patch_pre%cells%max_connectivity))
        END IF
        n_temp_edges = 0
        DO k = -1, 0
          DO ic = 1, n2_ilev_c(2*ilev+k)
            jc = flag2_c_list(2*ilev+k)%idx(ic)
            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges)
            DO i = 1, temp_num_edges
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
                   (/jc,i/), je)
              n_temp_edges = n_temp_edges + 1
              temp_edges(n_temp_edges) = je
              edge_cells(n_temp_edges) = jc
            END DO
          END DO
        END DO
        CALL quicksort(temp_edges(1:n_temp_edges), edge_cells(1:n_temp_edges))
        ! remove all edges of level 2*ilev-2 and 2*ilev-1
        DO k = -2, -1
          CALL remove_entries_from_ref_list(temp_edges(1:n_temp_edges), &
            n_temp_edges, flag2_e_list(2*ilev+k)%idx(1:n2_ilev_e(2*ilev+k)))
        END DO
        ! collect inner edges
        IF(n_temp_edges > 1) THEN
          IF (SIZE(pack_mask(:)) < n_temp_edges) THEN
            DEALLOCATE(pack_mask)
            ALLOCATE(pack_mask(n_temp_edges))
          END IF
          ! mask inner edges (occur twice)
          pack_mask(1:n_temp_edges-1) = (temp_edges(1:n_temp_edges-1) == &
                                         temp_edges(2:n_temp_edges)) .AND. &
                                        (temp_edges(1:n_temp_edges-1) > 0)
          pack_mask(n_temp_edges) = .FALSE.
          n = COUNT(pack_mask(1:n_temp_edges))
          ALLOCATE(flag2_e_list(2*ilev)%idx(n), flag2_e_list(2*ilev)%owner(n))
          flag2_e_list(2*ilev)%idx(:) = PACK(temp_edges(1:n_temp_edges), &
                                               pack_mask(1:n_temp_edges))
          CALL compute_edge_owner(flag2_e_list(2*ilev)%idx(:), &
            &                     flag2_e_list(2*ilev)%owner(:))
          n2_ilev_e(2*ilev) = n
        ELSE
          ALLOCATE(flag2_e_list(2*ilev)%idx(0), flag2_e_list(2*ilev)%owner(0))
          n2_ilev_e(2*ilev) = 0
          pack_mask(1:n_temp_edges) = .FALSE.
        END IF
        ! collect outer edges
        IF (n_temp_edges > 0) THEN
          ! mark inner edges
          pack_mask(2:n_temp_edges) = pack_mask(2:n_temp_edges) .OR. &
                                      pack_mask(1:n_temp_edges-1)
          ! mark invalid edges
          pack_mask(1:n_temp_edges) = pack_mask(1:n_temp_edges) .OR. &
                                      (temp_edges(1:n_temp_edges) <= 0)
          n = n_temp_edges - COUNT(pack_mask(1:n_temp_edges))
        ELSE
          n = 0
        END IF
        ALLOCATE(flag2_e_list(2*ilev+1)%idx(n), flag2_e_list(2*ilev+1)%owner(n))
        flag2_e_list(2*ilev+1)%idx(:) = PACK(temp_edges(1:n_temp_edges), &
                                               .NOT. pack_mask(1:n_temp_edges))
        CALL compute_edge_owner(flag2_e_list(2*ilev+1)%idx(:), &
          &                     flag2_e_list(2*ilev+1)%owner(:))
        n2_ilev_e(2*ilev+1) = n
        temp_edges(1:n) = flag2_e_list(2*ilev+1)%idx(:)
        edge_cells(1:n) = PACK(edge_cells(1:n_temp_edges), &
                               .NOT. pack_mask(1:n_temp_edges))
        n_temp_edges = n

        IF (SIZE(temp_vertices(:)) < SUM(n2_ilev_c(2*ilev-1:2*ilev)) * &
          &                              wrk_p_patch_pre%cells%max_connectivity) THEN

          DEALLOCATE(temp_vertices)
          ALLOCATE(temp_vertices(SUM(n2_ilev_c(2*ilev-1:2*ilev)) * &
            &                        wrk_p_patch_pre%cells%max_connectivity))
        END IF
        ! collect all vertices of level ilev + 1
        n_temp_vertices = 0
        DO k = -1, 0
          DO ic = 1, n2_ilev_c(2*ilev+k)

            jc = flag2_c_list(2*ilev+k)%idx(ic)

            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges)

            DO i = 1, temp_num_edges

              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/jc,i/),jv)

              IF (jv > 0) THEN
                n_temp_vertices = n_temp_vertices + 1
                temp_vertices(n_temp_vertices) = jv
              END IF
            END DO
          END DO
        END DO

        CALL insertion_sort(temp_vertices(1:n_temp_vertices))
        CALL remove_duplicated_entries(temp_vertices(1:n_temp_vertices), &
                                       n_temp_vertices)
        ! remove vertices that are on level ilev and ilev - 1
        DO k = -1, 0
          CALL remove_entries_from_ref_list(temp_vertices(1:n_temp_vertices), &
            n_temp_vertices, flag2_v_list(ilev+k)%idx(1:n2_ilev_v(ilev+k)))
        END DO

        IF (SIZE(temp_vertices_owner(:)) < n_temp_vertices) THEN
          DEALLOCATE(temp_vertices_owner)
          ALLOCATE(temp_vertices_owner(n_temp_vertices))
        END IF

        ALLOCATE(flag2_v_list(ilev+1)%idx(n_temp_vertices), &
          &      flag2_v_list(ilev+1)%owner(n_temp_vertices))
        flag2_v_list(ilev+1)%idx(:) = temp_vertices(1:n_temp_vertices)
        n2_ilev_v(ilev+1) = n_temp_vertices
        CALL compute_vertex_owner(flag2_v_list(ilev+1)%idx(:), &
          &                       flag2_v_list(ilev+1)%owner(:))
      END DO

      DEALLOCATE(temp_cells, temp_vertices, temp_edges, edge_cells, pack_mask)
    END SUBROUTINE compute_flag_lists

    SUBROUTINE compute_flag_lists_short(flag2_c_list, flag2_v_list, &
      &                                 flag2_e_list, n2_ilev_c, n2_ilev_v, &
      &                                 n2_ilev_e, n_boundary_rows &
    ! removed arguments because NEC has some problems with them...
    !  &                                 ,owned_edges, owned_verts &
                                        )
      INTEGER, INTENT(IN) :: n_boundary_rows
      TYPE(nb_flag_list_elem), INTENT(OUT) :: &
        flag2_c_list(0:2*n_boundary_rows), flag2_v_list(0:n_boundary_rows+1), &
        flag2_e_list(0:2*n_boundary_rows+1)
      INTEGER, INTENT(OUT) :: n2_ilev_c(0:2*n_boundary_rows), &
        &                     n2_ilev_v(0:n_boundary_rows+1), &
        &                     n2_ilev_e(0:2*n_boundary_rows+1)
      ! removed arguments because NEC has some problems with them...
      ! INTEGER, ALLOCATABLE, INTENT(OUT) :: owned_edges(:), owned_verts(:)

      INTEGER :: i
      LOGICAL :: position_error

      !---------------------------------------------------------------------
      ! flag_c_list(-1)%idx empty dummy list
      ! flag_c_list(0)%idx  all cells jg where cell_owner(jg) == my_proc
      ! flag_c_list(j)%idx j in 1..n_boundary_cells all cells bordering
      !                    on cells in flag_c_list(j - 1)
      ! flag_e_list(-1)%idx empty dummy list
      ! flag_e_list(j)%idx all edges where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the edge
      ! flag_v_list(-1)%idx empty dummy list
      ! flag_v_list(j)%idx all vertices where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the vertex
      !------------------------------------------------------------------------
      ! The purpose of the second set of flags is to control moving the
      ! halo points to the end of the index vector for nearly all cells
      ! even if order_type_of_halos = 1
      !------------------------------------------------------------------------
      ! flag2_e_list(0)%idx all edges which are owned by the local process
      ! flag2_e_list(1)%idx all edges adjacent to cells in flag_c_list(0)
      !                     but owned by other tasks
      ! flag2_e_list(2*i)%idx where i > 0, edges in idx are ONLY adjacent
      !                       to cells in flag_c_list(i)%idx
      ! flag2_e_list(2*i+1)%idx where i > 0, edges adjacent to one cell in
      !                         flag_c_list(i) and either
      !                            one cell in flag_c_list(i+1)
      !                         or edge is an outer edge,
      !                            i.e. adjacent to only one cell
      ! flag2_c_list(0)%idx == flag_c_list(0)%idx
      ! flag2_c_list(2*j-1)%idx where j > 0, contains all cells from
      !                         flag_c_list(j) which are
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_c_list(2*j)%idx   where j > 0, contains all cells from
      !                         flag_c_list(j) which are NOT
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_v_list(0)%idx contains all vertices from flag_v_list(0)%idx
      !                         owned by the current task
      ! flag2_v_list(1)%idx contains all vertices from flag_v_list(0)%idx
      !                     NOT owned by the current task
      ! flag2_v_list(j)%idx for j > 1 == flag_v_list(j-1)%idx
      !------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      ! compute flag2 arrays
      !--------------------------------------------------------------------------

      position_error = .FALSE.
      DO i = 1, wrk_p_patch_pre%n_patch_cells_g
        position_error = (owned_cells(i) /= i) .OR. position_error
      END DO
      IF (position_error) &
        CALL finish("compute_flag_lists_short", &
          &         "owned_cells content is inconsistent")

      ! set n2_ilev_cve arrays
      n2_ilev_c(:) = 0
      n2_ilev_v(:) = 0
      n2_ilev_e(:) = 0
      n2_ilev_c(0) = wrk_p_patch_pre%n_patch_cells_g
      n2_ilev_v(0) = wrk_p_patch_pre%n_patch_verts_g
      n2_ilev_e(0) = wrk_p_patch_pre%n_patch_edges_g

      ! allocate and set flag2_cve_list arrays
      ALLOCATE(flag2_c_list(0)%idx(wrk_p_patch_pre%n_patch_cells_g), &
        &      flag2_c_list(0)%owner(wrk_p_patch_pre%n_patch_cells_g), &
        &      flag2_v_list(0)%idx(wrk_p_patch_pre%n_patch_verts_g), &
        &      flag2_v_list(0)%owner(wrk_p_patch_pre%n_patch_verts_g), &
        &      flag2_v_list(1)%idx(0), &
        &      flag2_v_list(1)%owner(0), &
        &      flag2_e_list(0)%idx(wrk_p_patch_pre%n_patch_edges_g), &
        &      flag2_e_list(0)%owner(wrk_p_patch_pre%n_patch_edges_g), &
        &      flag2_e_list(1)%idx(0), &
        &      flag2_e_list(1)%owner(0))

      DO i = 1, wrk_p_patch_pre%n_patch_cells_g
        flag2_c_list(0)%idx(i) = i
      END DO
      flag2_c_list(0)%owner(:) = my_proc

      DO i = 1, wrk_p_patch_pre%n_patch_verts_g
        flag2_v_list(0)%idx(i) = i
      END DO
      flag2_v_list(0)%owner(:) = my_proc

      DO i = 1, wrk_p_patch_pre%n_patch_edges_g
        flag2_e_list(0)%idx(i) = i
      END DO
      flag2_e_list(0)%owner(:) = my_proc

      DO i = 1, n_boundary_rows

        ALLOCATE(flag2_c_list(2*i-1)%idx(0), flag2_c_list(2*i)%idx(0), &
          &      flag2_c_list(2*i-1)%owner(0), flag2_c_list(2*i)%owner(0), &
          &      flag2_v_list(i+1)%idx(0), flag2_v_list(i+1)%owner(0), &
          &      flag2_e_list(2*i)%idx(0), flag2_e_list(2*i+1)%idx(0), &
          &      flag2_e_list(2*i)%owner(0), flag2_e_list(2*i+1)%owner(0))
      END DO

      ! allocate and set owner arrays
      ALLOCATE(owned_verts(wrk_p_patch_pre%n_patch_verts_g), &
        &      owned_edges(wrk_p_patch_pre%n_patch_edges_g))
      owned_verts(:) = flag2_v_list(0)%idx(:)
      owned_edges(:) = flag2_e_list(0)%idx(:)
    END SUBROUTINE compute_flag_lists_short

    SUBROUTINE compute_vertex_owner(vertices, owner)
      INTEGER, INTENT(IN) :: vertices(:)
      INTEGER, INTENT(OUT) :: owner(:)

      INTEGER :: i, jv, n, j, jc, &
        &        t_cells(wrk_p_patch_pre%verts%max_connectivity), &
        &        t_cell_owner(wrk_p_patch_pre%verts%max_connectivity), &
        &        a_iown(2), a_mod_iown(2)
      LOGICAL :: swap

      ! compute ownership of vertices
      DO i = 1, SIZE(vertices)
        jv = vertices(i)
        CALL dist_mult_array_get(wrk_p_patch_pre%verts%num_edges, 1, (/jv/), n)
        DO j = 1, n
          CALL dist_mult_array_get(wrk_p_patch_pre%verts%cell, 1, (/jv, j/), &
            &                      t_cells(j))
        END DO
        CALL insertion_sort(t_cells(1:n))
        DO j = 1, n
          CALL dist_mult_array_get(dist_cell_owner, 1, (/t_cells(j)/), &
            &                      t_cell_owner(j))
        END DO
        jc = COUNT(t_cell_owner(1:n) >= 0)
        t_cell_owner(1:jc) = PACK(t_cell_owner(1:n), t_cell_owner(1:n) >= 0)
        n = jc
        a_iown(1) = t_cell_owner(1)
        a_mod_iown(1) = MOD(a_iown(1), 2)
        DO j = 2, n
          a_iown(2) = t_cell_owner(j)
          a_mod_iown(2) = MOD(a_iown(2), 2)
          ! 50% chance of acquiring a vertex
          swap = a_iown(1) > a_iown(2) &
               .NEQV. a_mod_iown(1) == a_mod_iown(2)
          a_iown(1) = MERGE(a_iown(2), a_iown(1), swap)
          a_mod_iown(1) = MERGE(a_mod_iown(2), a_mod_iown(1), swap)
        END DO
        owner(i) = a_iown(1)
      END DO
    END SUBROUTINE compute_vertex_owner

    SUBROUTINE compute_edge_owner(edges, owner)
      INTEGER, INTENT(IN) :: edges(:)
      INTEGER, INTENT(OUT) :: owner(:)

      INTEGER :: i, je, a_idx(2), a_iown(2), a_mod_iown(2)

      DO i = 1, SIZE(edges)

        je = edges(i)
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%cell,1,(/je,1/),a_idx(1))
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%cell,1,(/je,2/),a_idx(2))
        ! outer boundary edges always belong to single adjacent cell owner
        a_idx(:) = MERGE(a_idx(:), MAXVAL(a_idx(:)), a_idx(:) > 0)
        CALL dist_mult_array_get(dist_cell_owner, 1, (/a_idx(1)/), a_iown(1))
        CALL dist_mult_array_get(dist_cell_owner, 1, (/a_idx(2)/), a_iown(2))
        a_iown(:) = MERGE(a_iown(:), MAXVAL(a_iown(:)), a_iown(:) >= 0)
        a_mod_iown(1:2) = MOD(a_iown(1:2), 2)
        ! 50% chance of acquiring an edge
        owner(i) = MERGE(a_iown(1), a_iown(2), &
          &              (a_iown(1) > a_iown(2)) .EQV. &
          &              (a_mod_iown(1) == a_mod_iown(2)))
      END DO
    END SUBROUTINE compute_edge_owner

  END SUBROUTINE divide_patch

  FUNCTION refine_cells(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag == 0 &
         .OR. ( flag >= 1 .AND. irl0 == 1 ) &
         .OR. ( flag == 1 .AND. irl0 == 2 )
  END FUNCTION refine_cells

  FUNCTION refine_edges(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag==0 &
         .OR. ((irl0==1 .OR. irl0==2) .AND. flag>=1) &
         .OR. ((irl0==3 .OR. irl0==4) .AND. (flag==1 .OR. flag==2)) &
         .OR. ((irl0==5 .OR. irl0==6 .OR. irl0<0) .AND. flag==1)
  END FUNCTION refine_edges

  FUNCTION refine_verts(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag == 0 &
         .OR.  (irl0 == 1                .AND. flag  > 0) &
         .OR. ((irl0 == 2 .OR. irl0 < 0) .AND. flag == 1)
  END FUNCTION refine_verts

  SUBROUTINE build_patch_start_end(n_patch_cve, n_patch_cve_g, &
       patch_id, nblks, npromz, cell_type, &
       min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
       flag2_list, n2_ilev, decomp_info, &
       start_index, start_block, end_index, end_block, &
       start_g, end_g, order_type_of_halos, l_cell_correction, refin_ctrl, &
       refinement_predicate)
    INTEGER, INTENT(in) :: n_patch_cve, n_patch_cve_g, &
         patch_id, nblks, npromz, cell_type, &
         min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve
    INTEGER, INTENT(in) :: order_type_of_halos
    LOGICAL, INTENT(in) :: l_cell_correction
    TYPE(dist_mult_array) :: refin_ctrl
    INTEGER, INTENT(in) :: n2_ilev(0:)
    TYPE(nb_flag_list_elem), INTENT(in) :: flag2_list(0:)
    TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info
    INTEGER, DIMENSION(min_rlcve:), INTENT(inout) :: &
         start_index, start_block, end_index, end_block
    INTEGER, DIMENSION(min_rlcve:), INTENT(in) :: start_g, end_g
    INTERFACE
      FUNCTION refinement_predicate(flag, irl0) RESULT(p)
        INTEGER, INTENT(in) :: flag, irl0
        LOGICAL :: p
      END FUNCTION refinement_predicate
    END INTERFACE

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_setup_subdivision::build_patch_start_end'

    INTEGER :: i, ilev, ilev1, ilev_st, irl0, irlev, j, jb, jf, jl, &
         exec_sequence, ref_flag, k, k_start, n_inner, temp

    INTEGER, ALLOCATABLE :: temp_glb_index(:), permutation(:), temp_ilev(:), &
      &                     temp_owner(:)

    ! if all cells/vertices/edges have flag == 0
    IF ((n2_ilev(0) == n_patch_cve) .AND. (n2_ilev(0) == n_patch_cve_g)) THEN

      CALL build_patch_start_end_short(n_patch_cve, n_patch_cve_g, &
        patch_id, nblks, npromz, cell_type, &
        min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
        flag2_list, n2_ilev, decomp_info, &
        start_index, start_block, end_index, end_block, start_g, end_g, &
        order_type_of_halos, l_cell_correction)
      RETURN
    END IF

    ALLOCATE(temp_glb_index(n_patch_cve), temp_owner(n_patch_cve))

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      decomp_info%glb_index(1:n2_ilev(0)) = flag2_list(0)%idx(1:n2_ilev(0))
      decomp_info%owner_local(1:n2_ilev(0)) = flag2_list(0)%owner(1:n2_ilev(0))
      n_inner = n2_ilev(0)

    CASE (1)
      ALLOCATE(temp_ilev(n_patch_cve), permutation(n_patch_cve))
      k = 1
      DO ilev = 0, max_ilev
        temp_glb_index(k:k+n2_ilev(ilev)-1) = flag2_list(ilev)%idx(1:n2_ilev(ilev))
        temp_owner(k:k+n2_ilev(ilev)-1) = flag2_list(ilev)%owner(1:n2_ilev(ilev))
        temp_ilev(k:k+n2_ilev(ilev)-1) = ilev
        k = k + n2_ilev(ilev)
      END DO
      permutation(:) = (/(k, k = 1, n_patch_cve)/)
      CALL quicksort(temp_glb_index(:), permutation(:))
      temp_owner(:) = temp_owner(permutation(:))
      temp_ilev(:) = temp_ilev(permutation(:))
      DEALLOCATE(permutation)
      k = 1
      jf = 1
      j = 0
      DO WHILE (j < n_patch_cve_g .AND. jf <= n_patch_cve)
        j = temp_glb_index(jf)
        CALL dist_mult_array_get(refin_ctrl, 1, (/j/), irl0)
        IF (refinement_predicate(temp_ilev(jf), irl0)) THEN
          decomp_info%glb_index(k) = j
          decomp_info%owner_local(k) = temp_owner(jf)
          k = k + 1
        END IF
        jf = jf + 1
      END DO
      n_inner = k - 1
      DEALLOCATE(temp_ilev, temp_glb_index, temp_owner)

    CASE default
      CALL finish("", "Uknown order_type_of_halos")
    END SELECT

    CALL set_inner_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(1:n_inner), &
      &                      (/(i, i = 1, n_inner)/))

    ! Set start_index/block ... end_index/block for cells/verts/edges.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the inner global indices.

    ! Note: the segments between min_rlcell and min_rlcell_int-1 are reserved for
    ! halo cells; they are set below
    DO i=min_rlcve_int,max_rlcve
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, start_g(i), temp, +1)
      start_index(i) = idx_no(temp)
      start_block(i) = blk_no(temp)
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, end_g(i), temp, -1)
      end_index(i) = idx_no(temp)
      end_block(i) = blk_no(temp)
    ENDDO

    ! Preset remaining index sections with dummy values
    start_index(min_rlcve:min_rlcve_int-1) = -9999
    start_block(min_rlcve:min_rlcve_int-1) = -9999
    end_index(min_rlcve:min_rlcve_int-1)   = -9999
    end_block(min_rlcve:min_rlcve_int-1)   = -9999

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      ! processing local parent grid
      !
      ! Gather all halo cells/edges/verts at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlcve_int-1
      irlev = min_rlcve_int-1
      ref_flag = MERGE(0, 1, l_cell_correction .OR. order_type_of_halos == 2)
      ! write(0,*) "ref_flag=", ref_flag
      j = 1
      DO ilev = ref_flag + 1, max_ilev
        temp_glb_index(j:j+n2_ilev(ilev)-1) &
             = flag2_list(ilev)%idx(1:n2_ilev(ilev))
        temp_owner(j:j+n2_ilev(ilev)-1) &
             = flag2_list(ilev)%owner(1:n2_ilev(ilev))
        j = j + n2_ilev(ilev)
      END DO
      CALL quicksort(temp_glb_index(1:j-1), temp_owner(1:j-1))
      k_start = SUM(n2_ilev(0:ref_flag))
      IF (start_index(irlev)==-9999 .AND. k_start + 1 <= n_patch_cve) THEN
        start_index(irlev) = idx_no(k_start + 1) ! line index
        start_block(irlev) = blk_no(k_start + 1) ! block index
      ENDIF

      DO k = k_start + 1, n_patch_cve
        decomp_info%glb_index(k) = temp_glb_index(k - k_start)
        decomp_info%owner_local(k) = temp_owner(k - k_start)
      END DO
      DEALLOCATE(temp_glb_index, temp_owner)

      ! Set end index when the last point is processed
      end_index(min_rlcve:irlev) = idx_no(n_patch_cve)
      end_block(min_rlcve:irlev) = blk_no(n_patch_cve)
      start_index(min_rlcve:irlev-1) = start_index(irlev)
      start_block(min_rlcve:irlev-1) = start_block(irlev)

    CASE(1)
      k = n_inner
      ! Gather all halo cells/verts/edges except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2 value and will be placed at the
      ! index levels between min_rlcve_int-1 and min_rlcve
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here
      DO ilev = 1, max_ilev
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        DO j = 1, n2_ilev(ilev)
          jf = flag2_list(ilev)%idx(j)
          CALL dist_mult_array_get(refin_ctrl, 1, (/jf/), irl0)
          IF (.NOT. refinement_predicate(ilev, irl0)) THEN
            k = k + 1
            decomp_info%glb_index(k) = jf
            decomp_info%owner_local(k) = flag2_list(ilev)%owner(j)
            jb = blk_no(k) ! block index
            jl = idx_no(k) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (start_index(irlev)==-9999) THEN
              start_index(irlev) = jl
              start_block(irlev) = jb
            ENDIF
            end_index(irlev) = jl
            end_block(irlev) = jb
          END IF
        END DO

        ! Just in case that no grid point is found (may happen for ilev=1)
        IF (.NOT. l_cell_correction .AND. start_index(irlev)==-9999) THEN
          start_index(irlev) = end_index(irlev+1)+1
          start_block(irlev) = end_block(irlev+1)
          end_index(irlev)   = end_index(irlev+1)
          end_block(irlev)   = end_block(irlev+1)
        ENDIF

      ENDDO

      ! Fill start and end indices for remaining index sections
      IF (patch_id == 0) THEN
        ilev1 = min_rlcve_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlcve, min_rlcve_int - max_ilev)
        ilev_st = max_ilev + 1
      ENDIF
      IF (l_cell_correction .AND. cell_type==6) THEN ! for hexagons, there are no even-order halo cells
        DO ilev = 2, max_hw_cve, 2
          irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
          start_index(irlev) = end_index(irlev+1) + 1
          start_block(irlev) = end_block(irlev+1)
          end_index(irlev)   = end_index(irlev+1)
          end_block(irlev)   = end_block(irlev+1)
        ENDDO
      ENDIF
      DO ilev = ilev_st,max_hw_cve
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        start_index(irlev) = end_index(ilev1) + 1
        start_block(irlev) = end_block(ilev1)
        end_index(irlev)   = end_index(ilev1)
        end_block(irlev)   = end_block(ilev1)
      ENDDO
    END SELECT

    ! exec_sequence only serves the purpose to retain the execution
    ! sequence for these three correction as it was in the previous
    ! code version, it might be irrelevant.
    ! see redmine issue #3924
    DO exec_sequence = 0, 1
      IF (.NOT. l_cell_correction .AND. exec_sequence == 0 &
           .OR. l_cell_correction .AND. exec_sequence == 1) THEN
        ! If a PE owns only nest boundary points, it may happen that
        ! one or more index sections of halo cells/verts/edges are
        ! empty. These are filled here
        IF (start_block(0)     > end_block(min_rlcve_int) &
             .OR. start_block(0) == end_block(min_rlcve_int) &
             .AND. start_index(0) > end_index(min_rlcve_int) ) THEN
          DO i = min_rlcve_int-1, min_rlcve, -1
            IF (start_index(i) == -9999 .OR. &
                 start_block(i) == -9999 .OR. &
                 end_index(i)   == -9999 .OR. &
                 end_block(i)   == -9999 ) THEN
              start_index(i) = start_index(i+1)
              start_block(i) = start_block(i+1)
              end_index(i)   = end_index(i+1)
              end_block(i)   = end_block(i+1)
            ENDIF
          ENDDO
        ENDIF

        ! Finally, fill start indices of halo rows with meaningful
        ! values for empty patches (which occur when processor
        ! splitting is applied)
        IF (nblks <= 0 .OR. npromz <= 0) THEN
          DO i=min_rlcve,min_rlcve_int-1
            start_index(i) = start_index(min_rlcve_int)
            start_block(i) = start_block(min_rlcve_int)
            end_index(i)   = end_index(min_rlcve_int)
            end_block(i)   = end_block(min_rlcve_int)
          ENDDO
        ENDIF
      END IF
      ! special treatment for sequential runs with trivial decomposition
      IF (exec_sequence == 0 .AND. SUM(n2_ilev(1:)) == 0) THEN
        end_index(min_rlcve:min_rlcve_int-1) &
             = end_index(min_rlcve_int)
        end_block(min_rlcve:min_rlcve_int-1) &
             = end_block(min_rlcve_int)
        IF (end_index(min_rlcve_int) == nproma) THEN
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int) + 1
          start_index(min_rlcve:min_rlcve_int-1) = 1
        ELSE
          start_index(min_rlcve:min_rlcve_int-1) = &
               &    end_index(min_rlcve_int) + 1
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int)
        END IF
      END IF

    END DO

    CALL set_outer_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(n_inner+1:), &
      &                      (/(i, i = n_inner+1, n_patch_cve)/))

  END SUBROUTINE build_patch_start_end

  SUBROUTINE build_patch_start_end_short(n_patch_cve, n_patch_cve_g, &
       patch_id, nblks, npromz, cell_type, &
       min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
       flag2_list, n2_ilev, decomp_info, &
       start_index, start_block, end_index, end_block, start_g, end_g, &
       order_type_of_halos, l_cell_correction)
    INTEGER, INTENT(in) :: n_patch_cve, n_patch_cve_g, &
         patch_id, nblks, npromz, cell_type, &
         min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve
    INTEGER, INTENT(in) :: order_type_of_halos
    LOGICAL, INTENT(in) :: l_cell_correction
    INTEGER, INTENT(in) :: n2_ilev(0:)
    TYPE(nb_flag_list_elem), INTENT(in) :: flag2_list(0:)
    TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info
    INTEGER, DIMENSION(min_rlcve:), INTENT(inout) :: &
         start_index, start_block, end_index, end_block
    INTEGER, DIMENSION(min_rlcve:), INTENT(in) :: start_g, end_g

    INTEGER :: i, ilev, ilev1, ilev_st, irlev, exec_sequence, temp

    ! if all cells/vertices/edges have flag == 0
    IF ((n2_ilev(0) /= n_patch_cve) .OR. (n2_ilev(0) /= n_patch_cve_g)) &
      CALL finish("build_patch_start_end_short", &
          &       "arguments do not fit for this routine")

    IF (ANY(order_type_of_halos == (/0,1,2/))) THEN
      decomp_info%glb_index(:) = flag2_list(0)%idx(:)
      decomp_info%owner_local(:) = flag2_list(0)%owner(:)
    ELSE
      CALL finish("build_patch_start_end_short", "Uknown order_type_of_halos")
    END IF

    CALL set_inner_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(:), &
      &                      (/(i, i = 1, n_patch_cve)/))

    ! Set start_index/block ... end_index/block for cells/verts/edges.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the inner global indices.

    ! Note: the segments between min_rlcell and min_rlcell_int-1 are reserved for
    ! halo cells; they are set below
    DO i=min_rlcve_int,max_rlcve
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, start_g(i), temp, +1)
      start_index(i) = idx_no(temp)
      start_block(i) = blk_no(temp)
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, end_g(i), temp, -1)
      end_index(i) = idx_no(temp)
      end_block(i) = blk_no(temp)
    ENDDO

    ! Preset remaining index sections with dummy values
    start_index(min_rlcve:min_rlcve_int-1) = -9999
    start_block(min_rlcve:min_rlcve_int-1) = -9999
    end_index(min_rlcve:min_rlcve_int-1)   = -9999
    end_block(min_rlcve:min_rlcve_int-1)   = -9999

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      ! processing local parent grid
      !
      ! Gather all halo cells/edges/verts at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlcve_int-1
      irlev = min_rlcve_int-1

      ! Set end index when the last point is processed
      start_index(min_rlcve:irlev-1) = start_index(irlev)
      start_block(min_rlcve:irlev-1) = start_block(irlev)
      end_index(min_rlcve:irlev) = idx_no(n_patch_cve)
      end_block(min_rlcve:irlev) = blk_no(n_patch_cve)

    CASE(1)
      ! Gather all halo cells/verts/edges except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2 value and will be placed at the
      ! index levels between min_rlcve_int-1 and min_rlcve
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here

      irlev = MAX(min_rlcve, min_rlcve_int - 1)  ! index section into which the halo points are put

      ! Just in case that no grid point is found (may happen for ilev=1)
      IF (.NOT. l_cell_correction) THEN
        start_index(irlev) = end_index(irlev+1)+1
        start_block(irlev) = end_block(irlev+1)
        end_index(irlev)   = end_index(irlev+1)
        end_block(irlev)   = end_block(irlev+1)
      ENDIF

      ! Fill start and end indices for remaining index sections
      IF (patch_id == 0) THEN
        ilev1 = min_rlcve_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlcve, min_rlcve_int - max_ilev)
        ilev_st = max_ilev + 1
      ENDIF
      IF (l_cell_correction .AND. cell_type==6) THEN ! for hexagons, there are no even-order halo cells
        DO ilev = 2, max_hw_cve, 2
          irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
          start_index(irlev) = end_index(irlev+1) + 1
          start_block(irlev) = end_block(irlev+1)
          end_index(irlev)   = end_index(irlev+1)
          end_block(irlev)   = end_block(irlev+1)
        ENDDO
      ENDIF
      DO ilev = ilev_st,max_hw_cve
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        start_index(irlev) = end_index(ilev1) + 1
        start_block(irlev) = end_block(ilev1)
        end_index(irlev)   = end_index(ilev1)
        end_block(irlev)   = end_block(ilev1)
      ENDDO
    END SELECT

    ! exec_sequence only serves the purpose to retain the execution
    ! sequence for these three correction as it was in the previous
    ! code version, it might be irrelevant.
    ! see redmine issue #3924
    DO exec_sequence = 0, 1
      IF (.NOT. l_cell_correction .AND. exec_sequence == 0 &
           .OR. l_cell_correction .AND. exec_sequence == 1) THEN
        ! If a PE owns only nest boundary points, it may happen that
        ! one or more index sections of halo cells/verts/edges are
        ! empty. These are filled here
        IF (start_block(0)     > end_block(min_rlcve_int) &
             .OR. start_block(0) == end_block(min_rlcve_int) &
             .AND. start_index(0) > end_index(min_rlcve_int) ) THEN
          DO i = min_rlcve_int-1, min_rlcve, -1
            IF (start_index(i) == -9999 .OR. &
                 start_block(i) == -9999 .OR. &
                 end_index(i)   == -9999 .OR. &
                 end_block(i)   == -9999 ) THEN
              start_index(i) = start_index(i+1)
              start_block(i) = start_block(i+1)
              end_index(i)   = end_index(i+1)
              end_block(i)   = end_block(i+1)
            ENDIF
          ENDDO
        ENDIF

        ! Finally, fill start indices of halo rows with meaningful
        ! values for empty patches (which occur when processor
        ! splitting is applied)
        IF (nblks <= 0 .OR. npromz <= 0) THEN
          DO i=min_rlcve,min_rlcve_int-1
            start_index(i) = start_index(min_rlcve_int)
            start_block(i) = start_block(min_rlcve_int)
            end_index(i)   = end_index(min_rlcve_int)
            end_block(i)   = end_block(min_rlcve_int)
          ENDDO
        ENDIF
      END IF
      ! special treatment for sequential runs with trivial decomposition
      IF (exec_sequence == 0 .AND. SUM(n2_ilev(1:)) == 0) THEN
        end_index(min_rlcve:min_rlcve_int-1) &
             = end_index(min_rlcve_int)
        end_block(min_rlcve:min_rlcve_int-1) &
             = end_block(min_rlcve_int)
        IF (end_index(min_rlcve_int) == nproma) THEN
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int) + 1
          start_index(min_rlcve:min_rlcve_int-1) = 1
        ELSE
          start_index(min_rlcve:min_rlcve_int-1) = &
               &    end_index(min_rlcve_int) + 1
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int)
        END IF
      END IF
    END DO

  END SUBROUTINE build_patch_start_end_short

  !-------------------------------------------------------------------------
  !>
  !!               Calculates local indices l_i from global indices g_i
  !!               using the mapping in decomp_info
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE get_local_idx(decomp_info, g_i, l_i, opt_mode)

    !
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    INTEGER, INTENT(in) :: g_i
    INTEGER, INTENT(in), OPTIONAL :: opt_mode

    INTEGER, INTENT(out) :: l_i

    INTEGER :: mode

    !-----------------------------------------------------------------------

    ! mode controls the behaviour if the global index is valid, but outside the local domain:
    ! mode > 0 : Map it to the next greater index (even if this is outside the local domain)
    ! mode < 0 : Map it to the next smaller index (even if this is 0)
    ! mode = 0 : Map it to a negative value (negative of global index)

    IF(PRESENT(opt_mode)) THEN
      mode = opt_mode
    ELSE
      mode = 0
    ENDIF

    IF (mode == 0) THEN
      l_i = get_local_index(decomp_info%glb2loc_index, g_i)
      IF (l_i < 0) l_i = -g_i
    ELSE IF (mode > 0) THEN
      l_i = get_valid_local_index(decomp_info%glb2loc_index, g_i, .TRUE.)
    ELSE
      l_i = get_valid_local_index(decomp_info%glb2loc_index, g_i)
    END IF

  END SUBROUTINE get_local_idx

  !-------------------------------------------------------------------------
  !>
  !!               Calculates local line/block indices l_idx, l_blk
  !!               from global line/block indices g_idx, g_blk
  !!               using the mapping in decomp_info
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE get_local_idx_blk(decomp_info, g_idx, g_blk, l_idx, l_blk, opt_mode)

    !
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    INTEGER, INTENT(in) :: g_idx, g_blk
    INTEGER, INTENT(in), OPTIONAL :: opt_mode

    INTEGER, INTENT(out) :: l_idx, l_blk

    INTEGER j_l

    CALL get_local_idx(decomp_info, idx_1d(g_idx, g_blk), j_l, opt_mode)

    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE get_local_idx_blk
  !-------------------------------------------------------------------------
  !>
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE divide_subset_geometric(subset_flag, n_proc, wrk_p_patch_pre, &
                                     dist_cell_owner, lparent_level)

    INTEGER, ALLOCATABLE   :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre ! out is required for
                                                        ! the distributed arrays
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner ! receives the owner PE for every cell
    ! (-1 for cells not in subset)
    LOGICAL, INTENT(IN)    :: lparent_level ! indicates if domain decomposition is executed for
                              ! the parent grid level (true) or for the target grid level (false)

    INTEGER :: j, nc, nc_g, ncs, nce, n_onb_points, jm(1)
    INTEGER :: count_physdom(0:max_phys_dom), count_total, id_physdom(max_phys_dom), &
               num_physdom, proc_count(max_phys_dom), proc_offset(max_phys_dom), checksum, &
               ncell_offset(0:max_phys_dom), ncell_offset_g(0:max_phys_dom)
    TYPE(t_cell_info), ALLOCATABLE :: cell_desc(:)
    REAL(wp) :: corr_ratio(max_phys_dom)
    LOGICAL  :: lsplit_merged_domains, locean
    INTEGER :: my_cell_start, my_cell_end
    INTEGER, POINTER :: owner_chunk_ptr(:)
    INTEGER :: ierror
    CHARACTER(*), PARAMETER :: routine = 'divide_subset_geometric'

    !-----------------------------------------------------------------------

    IF(p_pe_work==0) WRITE(0,*) 'divide_patch: Using geometric area subdivision'

    IF (get_my_process_type() == ocean_process) THEN
      locean = .TRUE.
    ELSE
      locean = .FALSE.
    ENDIF

    lsplit_merged_domains = .FALSE.

    my_cell_start = LBOUND(subset_flag, 1)
    my_cell_end = UBOUND(subset_flag, 1)
    IF (p_n_work > 1) THEN
#ifndef NOMPI
      CALL init_divide_cells_by_location_mpi
#else
      CALL finish("number of processors > 1 but MPI unavailable!")
#endif
    END IF

    ! Check if domain merging has been applied; for large PE numbers, calculating the DD
    ! for each physical domain separately tends to yield a more balanced distribution
    ! of the halo points
    IF (ldiv_phys_dom) THEN
      count_physdom(:) = 0
      count_total      = 0
      num_physdom      = 0
      id_physdom(:)    = 0
      proc_count(:)    = 0
      ncell_offset(:)  = 0
      proc_offset(:)   = 0
      corr_ratio(:)    = 1._wp

      DO j = my_cell_start, my_cell_end
        IF (subset_flag(j) > 0) THEN
          count_physdom(subset_flag(j)) = count_physdom(subset_flag(j)) + 1
          count_total = count_total + 1
        ENDIF
      ENDDO
      IF (p_n_work > 1) THEN
#ifndef NOMPI
        count_physdom(0) = count_total
        CALL mpi_allreduce(mpi_in_place, count_physdom, SIZE(count_physdom), &
             p_int, mpi_sum, p_comm_work, ierror)
        IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
        count_total = count_physdom(0)
#endif
      END IF

      IF (MAXVAL(count_physdom(1:)) < count_total) THEN
        lsplit_merged_domains = .TRUE.
        DO j = 1, max_phys_dom
          IF (count_physdom(j) > 0) THEN
            num_physdom = num_physdom + 1
            id_physdom(num_physdom) = j
            proc_count(num_physdom) = NINT(REAL(count_physdom(j),wp)/REAL(count_total,wp)*REAL(n_proc,wp))
          ENDIF
        ENDDO

        ! Ensure that the sum of partial processor counts matches n_proc
        checksum = SUM(proc_count(1:num_physdom)) - n_proc
        IF (checksum /= 0) THEN
          DO j = 1, num_physdom
            corr_ratio(j) = REAL(proc_count(j),wp) / &
              (REAL(count_physdom(id_physdom(j)),wp)/REAL(count_total,wp)*REAL(n_proc,wp))
          ENDDO
        ENDIF
        IF (checksum > 0) THEN
          DO WHILE (checksum > 0)
            jm = MAXLOC(corr_ratio)
            j = jm(1)
            corr_ratio(j) = 1._wp
            proc_count(j) = proc_count(j) - 1
            checksum = checksum - 1
          ENDDO
        ELSE IF (checksum < 0) THEN
          DO WHILE (checksum < 0)
            jm = MINLOC(corr_ratio)
            j = jm(1)
            corr_ratio(j) = 1._wp
            proc_count(j) = proc_count(j) + 1
            checksum = checksum + 1
          ENDDO
        ENDIF

        ! Compute offset for processor IDs
        DO j = 2, num_physdom
          proc_offset(j) = proc_offset(j-1) + proc_count(j-1)
        ENDDO

        IF(p_pe_work==0) THEN
          WRITE(0,*) 'divide_patch: partial processor counts used for decomposition of merged domain:'
          WRITE(0,*) proc_count(1:num_physdom), 'total: ', n_proc
        ENDIF

      ENDIF

    ENDIF

    IF (.NOT. lsplit_merged_domains) THEN
      num_physdom       = 1
      proc_count(1)     = n_proc
      proc_offset(1)    = 0
      IF (SIZE(subset_flag) > 0) THEN
        id_physdom(1) = p_max(MAXVAL(subset_flag), p_comm_work)
      ELSE
        id_physdom(1) = p_max(0, p_comm_work)
      END IF
    ENDIF


    ! Fill the cell_desc array, it must contain:
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner
    ALLOCATE(cell_desc(my_cell_start:my_cell_end))

    CALL build_cell_info(my_cell_start, my_cell_end, cell_desc, &
         wrk_p_patch_pre, num_physdom, subset_flag, lsplit_merged_domains, &
         lparent_level, locean, id_physdom, ncell_offset, n_onb_points)
    IF (p_n_work > 1) THEN
#ifndef NOMPI
      CALL mpi_allreduce(ncell_offset, ncell_offset_g, SIZE(ncell_offset), &
           p_int, mpi_sum, p_comm_work, ierror)
      IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
#endif
    END IF


    DO j = 1, num_physdom

      ncs = ncell_offset(j-1) + my_cell_start
      nce = ncell_offset(j) + my_cell_start - 1
      nc  = ncell_offset(j) - ncell_offset(j-1)

      IF (p_n_work > 1) THEN
#ifndef NOMPI
        nc_g = ncell_offset_g(j) - ncell_offset_g(j-1)
        CALL divide_cells_by_location_mpi(nc_g, nc, cell_desc(ncs:nce), &
             proc_count(j), p_comm_work)
#endif
      ELSE
        CALL divide_cells_by_location(nc, cell_desc(ncs:nce), &
             0, proc_count(j)-1)
        ! After divide_cells_by_location the cells are sorted by owner,
        ! order them by original cell numbers again
        CALL sort_cell_info_by_cell_number(cell_desc(ncs:nce), nc)
      END IF

      ! Apply shift of processor IDs
      IF (j > 1) THEN
        cell_desc(ncs:nce)%owner = cell_desc(ncs:nce)%owner + proc_offset(j)
      END IF

    ENDDO

    ! set up distributed owner array
    CALL create_dist_cell_owner(dist_cell_owner,wrk_p_patch_pre%n_patch_cells_g)

    ! Initialization of cell owner field
    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, owner_chunk_ptr)
    owner_chunk_ptr = -1

    IF (p_n_work > 1) THEN
#ifndef NOMPI
      CALL set_owners_mpi(my_cell_start, my_cell_end, dist_cell_owner, &
           wrk_p_patch_pre, cell_desc, ncell_offset, num_physdom, n_onb_points)
#endif
    ELSE
      CALL set_owners(1, wrk_p_patch_pre%n_patch_cells_g, &
           owner_chunk_ptr, wrk_p_patch_pre, cell_desc, ncell_offset, &
           num_physdom, n_onb_points)
    END IF

    DEALLOCATE(cell_desc)
  END SUBROUTINE divide_subset_geometric

  SUBROUTINE build_cell_info(range_start, range_end, cell_desc, &
       wrk_p_patch_pre, num_physdom, subset_flag, lsplit_merged_domains, &
       lparent_level, locean, id_physdom, ncell_offset, n_onb_points)
    INTEGER, INTENT(in) :: range_start, range_end
    TYPE(t_cell_info), INTENT(out) :: cell_desc(range_start:range_end)
    TYPE(t_pre_patch), INTENT(inout) :: wrk_p_patch_pre
    LOGICAL, INTENT(in) :: lsplit_merged_domains, lparent_level, locean
    INTEGER, INTENT(in) :: num_physdom
    ! also ends on range_end, but we want implicit shape argument
    INTEGER, INTENT(in) :: subset_flag(range_start:range_end)
    INTEGER, INTENT(inout) :: ncell_offset(0:max_phys_dom)
    INTEGER, INTENT(out) :: n_onb_points
    INTEGER, INTENT(in) :: id_physdom(max_phys_dom)

    INTEGER :: i, j, nc, jd, idp, jv, refin_ctrl
    REAL(wp) :: cclat, cclon, temp_lat, temp_lon

    nc = 0
    n_onb_points = 0

    IF (wrk_p_patch_pre%cell_type == 6) THEN

      DO j = range_start, range_end

        IF (subset_flag(j) <= 0) CYCLE

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%center, 1, (/j/), cclat)
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%center, 2, (/j/), cclon)
        cell_desc(range_start + nc)%lat = fxp_lat(cclat)
        cell_desc(range_start + nc)%lon = fxp_lon(cclon)
        cell_desc(range_start + nc)%cell_number = j
        cell_desc(range_start + nc)%owner = 0

        nc = nc + 1 ! Cell counter
      ENDDO

      ncell_offset(1) = nc

    ELSE ! domain decomposition for triangular cells with optional splitting into physical domains

      DO jd = 1, num_physdom

        idp = id_physdom(jd)

!CDIR NODEP
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          CALL dist_mult_array_get(wrk_p_patch_pre%cells%refin_ctrl, 1, (/j/), &
            &                      refin_ctrl)

          ! Disregard outer nest boundary points for the time being. They do very little
          ! computational work, so they can be added to the closest PEs afterwards
          IF (lparent_level .AND. refin_ctrl == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl >= 1 .AND. refin_ctrl <= 3 .AND. .NOT. locean) THEN
            cell_desc(range_end-n_onb_points)%cell_number = j
            n_onb_points = n_onb_points + 1
            CYCLE
          ENDIF

          CALL dist_mult_array_get(wrk_p_patch_pre%cells%center, 1, (/j/), cclat)
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%center, 2, (/j/), cclon)

          ! Using the center of the cells for geometric subdivision leads
          ! to "toothed" edges of the subdivision area
          ! Thus we use the minimum lat/lon as subdivision criterion.

          IF (cclat >= 0._wp) THEN
            DO i = 1, 3
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/j,i/), jv)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%vertex, 1, (/jv/),&
                &                      temp_lat)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%vertex, 2, (/jv/),&
                &                      temp_lon)
              cclat = MAX(cclat, temp_lat)
              cclon = MAX(cclon, temp_lon)
            ENDDO
          ELSE
            DO i = 1, 3
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/j,i/), jv)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%vertex, 1, (/jv/),&
                &                      temp_lat)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%vertex, 2, (/jv/),&
                &                      temp_lon)
              cclat = MIN(cclat, temp_lat)
              cclon = MAX(cclon, temp_lon)
            ENDDO
          ENDIF

          cell_desc(range_start + nc)%lat = fxp_lat(cclat)
          cell_desc(range_start + nc)%lon = fxp_lon(cclon)
          cell_desc(range_start + nc)%cell_number = j
          cell_desc(range_start + nc)%owner = 0

          nc = nc+1 ! Cell counter

        ENDDO

        ncell_offset(jd) = nc

      ENDDO


    ENDIF

  END SUBROUTINE build_cell_info

  ! Set owner list (of complete patch)
  SUBROUTINE set_owners(range_start, range_end, &
       owner, wrk_p_patch_pre, cell_desc, ncell_offset, &
       num_physdom, n_onb_points)
    INTEGER, INTENT(in) :: range_start, range_end
     !> receives the owner PE for every cell
    INTEGER, INTENT(out) :: owner(range_start:range_end)
    TYPE(t_pre_patch), INTENT(inout) :: wrk_p_patch_pre
    TYPE(t_cell_info), INTENT(inout) :: cell_desc(range_start:range_end)
    !> if > 0 a cell belongs to the subset
    INTEGER, INTENT(in) :: ncell_offset(0:max_phys_dom)
    INTEGER, INTENT(in) :: num_physdom, n_onb_points

    INTEGER :: i, ii, j, jd, jj, jn, &
         first_unassigned_onb, last_unassigned_onb, &
         ncs, nce, num_set_on_this_pass
    TYPE(t_cell_info) :: temp
    INTEGER, POINTER :: cells_num_edges_local_ptr(:), &
      &                 cells_neighbor_local_ptr(:,:)
    INTEGER, ALLOCATABLE :: set_on_this_pass(:)

    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_num_edges, &
      &                            cells_num_edges_local_ptr)
    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_neighbor, &
      &                            cells_neighbor_local_ptr)

    DO jd = 1, num_physdom
      ncs = ncell_offset(jd-1)+1
      nce = ncell_offset(jd)
      DO j = ncs, nce
        owner(cell_desc(j)%cell_number) = cell_desc(j)%owner
      END DO
    ENDDO

    ! Add outer nest boundary points that have been disregarded so far
    IF (n_onb_points > 0) THEN
      first_unassigned_onb = range_end - n_onb_points + 1
      last_unassigned_onb = range_end
      ! Iterations are needed because outer nest boundary row contains
      ! indirect neighbors
      ALLOCATE(set_on_this_pass(n_onb_points))
      DO WHILE (first_unassigned_onb <= last_unassigned_onb)
        num_set_on_this_pass = 0
        DO i = last_unassigned_onb, first_unassigned_onb, -1
          j = cell_desc(i)%cell_number
          DO ii = 1, cells_num_edges_local_ptr(j)
            jn = cells_neighbor_local_ptr(j, ii)
            IF (jn > 0) THEN
              IF (owner(jn) >= 0) THEN
                ! move found cells to end of onb area
                temp = cell_desc(i)
                DO jj = i, last_unassigned_onb - 1
                  cell_desc(jj) = cell_desc(jj + 1)
                END DO
                cell_desc(last_unassigned_onb) = temp
                owner(j) = -owner(jn) - 1
                num_set_on_this_pass = num_set_on_this_pass + 1
                set_on_this_pass(num_set_on_this_pass) = j
                last_unassigned_onb = last_unassigned_onb - 1
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO i = 1, num_set_on_this_pass
          j = set_on_this_pass(i)
          owner(j) = -owner(j) - 1
        END DO
      ENDDO
    ENDIF

  END SUBROUTINE set_owners

#ifndef NOMPI
  ! Set owner list (of complete patch)
  SUBROUTINE set_owners_mpi(range_start, range_end, dist_cell_owner, &
       wrk_p_patch_pre, cell_desc, ncell_offset, num_physdom, n_onb_points)
    USE mpi
    USE ppm_extents
    USE ppm_distributed_array
    INTEGER, INTENT(in) :: range_start, range_end
     !> receives the owner PE for every cell
    TYPE(dist_mult_array), INTENT(inout) :: dist_cell_owner
    TYPE(t_pre_patch), INTENT(inout) :: wrk_p_patch_pre
    TYPE(t_cell_info), INTENT(inout) :: cell_desc(range_start:range_end)
    !> if > 0 a cell belongs to the subset
    INTEGER, INTENT(in) :: ncell_offset(0:max_phys_dom)
    INTEGER, INTENT(in) :: num_physdom, n_onb_points

    INTEGER :: i, ii, j, jd, jj, jn, first_unassigned_onb, &
      &        last_unassigned_onb, ncs, nce, owner_value, owner_idx(1), ierror
    TYPE(t_cell_info) :: temp
    INTEGER, POINTER :: owners_local_ptr(:), cells_num_edges_local_ptr(:), &
      &                 cells_neighbor_local_ptr(:,:)
    INTEGER, ALLOCATABLE :: tmp_owners_local(:)
    LOGICAL :: owner_missing, any_onb_points
    CHARACTER(14), PARAMETER :: routine = 'set_owners_mpi'

    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, owners_local_ptr)
    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_num_edges, &
      &                            cells_num_edges_local_ptr)
    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_neighbor, &
      &                            cells_neighbor_local_ptr)

    DO jd = 1, num_physdom
      ncs = ncell_offset(jd-1) + range_start
      nce = ncell_offset(jd) + range_start - 1
      DO j = ncs, nce
        owners_local_ptr(cell_desc(j)%cell_number) = cell_desc(j)%owner
      END DO
    END DO

    ALLOCATE(tmp_owners_local(range_start:range_end))
    tmp_owners_local = owners_local_ptr

    ! Add outer nest boundary points that have been disregarded so far
    any_onb_points = n_onb_points > 0
    CALL mpi_allreduce(mpi_in_place, any_onb_points, 1, p_bool, &
         mpi_lor, p_comm_work, ierror)
    IF (any_onb_points) THEN

      first_unassigned_onb = range_end - n_onb_points + 1
      last_unassigned_onb = range_end
      ! Iterations are needed because outer nest boundary row contains
      ! indirect neighbours
      owner_missing = first_unassigned_onb <= last_unassigned_onb
      CALL mpi_allreduce(mpi_in_place, owner_missing, 1, p_bool, &
           mpi_lor, p_comm_work, ierror)
      DO WHILE (owner_missing)
        IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
        CALL dist_mult_array_expose(dist_cell_owner)
        DO i = last_unassigned_onb, first_unassigned_onb, -1
          j = cell_desc(i)%cell_number
          DO ii = 1, cells_num_edges_local_ptr(j)
            jn = cells_neighbor_local_ptr(j, ii)
            IF (jn > 0) THEN
              owner_idx(1) = jn
              CALL dist_mult_array_get(dist_cell_owner, 1, owner_idx, &
                &                      owner_value)
              IF (owner_value >= 0) THEN
                ! move found cells to end of onb area
                temp = cell_desc(i)
                DO jj = i, last_unassigned_onb - 1
                  cell_desc(jj) = cell_desc(jj + 1)
                END DO
                cell_desc(last_unassigned_onb) = temp
                tmp_owners_local(j) = owner_value
                last_unassigned_onb = last_unassigned_onb - 1
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        CALL dist_mult_array_unexpose(dist_cell_owner)
        owners_local_ptr = tmp_owners_local
        owner_missing = first_unassigned_onb <= last_unassigned_onb
        CALL mpi_allreduce(mpi_in_place, owner_missing, 1, p_bool, &
             mpi_lor, p_comm_work, ierror)
      ENDDO
    ENDIF

  END SUBROUTINE set_owners_mpi
#endif
  !-------------------------------------------------------------------------

END MODULE mo_setup_subdivision
