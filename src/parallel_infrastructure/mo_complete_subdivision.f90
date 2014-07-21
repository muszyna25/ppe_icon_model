!>
!!               This module provides all routines for dividing patches
!!
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
MODULE mo_complete_subdivision
  !-------------------------------------------------------------------------
  USE mo_impl_constants,     ONLY: min_rlcell, min_rledge, &
    & min_rlcell_int, min_rledge_int, max_phys_dom
  USE mo_exception,          ONLY: finish, message, message_text

  USE mo_model_domain,       ONLY: t_patch, p_patch, p_patch_local_parent, &
    &                              p_phys_patch
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info, get_local_index, &
    &                              get_valid_local_index, t_glb2loc_index_lookup
  USE mo_mpi,                ONLY: p_send, p_recv, p_max, p_min, proc_split, p_sum
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_COMM_NULL
#endif
  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_test, &
    & my_process_is_mpi_seq, process_mpi_all_test_id, process_mpi_all_workroot_id, &
    & my_process_is_mpi_workroot, p_pe_work, p_n_work, my_process_is_mpi_parallel

  USE mo_parallel_config,    ONLY:  p_test_run
  USE mo_communication,      ONLY: setup_comm_pattern, blk_no, idx_no, idx_1d, &
    &                              setup_comm_gather_pattern, t_comm_gather_pattern, &
    &                              ASSIGNMENT(=), delete_comm_gather_pattern
  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
    & grf_bdyintp_end_c, grf_fbk_start_c, grf_fbk_start_e, grf_bdywidth_c, &
    & grf_bdywidth_e
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, n_phys_dom
  USE mo_dist_dir,            ONLY: dist_dir_get_owners
  IMPLICIT NONE

  PRIVATE

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: copy_processor_splitting
  PUBLIC :: set_patch_communicators
  PUBLIC :: finalize_decomposition
  PUBLIC :: setup_phys_patches
  PUBLIC :: set_comm_pat_gather
  PUBLIC :: complete_parallel_setup


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!  In case of a test run: Copies processor splitting to test PE
  SUBROUTINE copy_processor_splitting(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER :: ibuf(n_dom_start:n_dom,2)

    IF(.NOT. p_test_run) RETURN ! Nothing to do

    IF(my_process_is_mpi_test()) THEN
      CALL p_recv(proc_split, process_mpi_all_workroot_id, 1)
      CALL p_recv(ibuf, process_mpi_all_workroot_id, 2)
      p_patch(:)%n_proc = ibuf(:,1)
      p_patch(:)%proc0  = ibuf(:,2)

    ELSEIF(my_process_is_mpi_workroot()) THEN
      CALL p_send(proc_split, process_mpi_all_test_id, 1)
      ibuf(:,1) = p_patch(:)%n_proc
      ibuf(:,2) = p_patch(:)%proc0
      CALL p_send(ibuf, process_mpi_all_test_id, 2)

    ENDIF


  END SUBROUTINE copy_processor_splitting
  !-------------------------------------------------------------------------
  !>
  !!  Sets the communicators in the patches if these have been read.
  SUBROUTINE set_patch_communicators(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER jc, jgc, jg, jgp, n_proc_total, comm, mpierr
    INTEGER, ALLOCATABLE :: patch_no(:)

    IF (my_process_is_mpi_seq()) THEN
      p_patch(:)%comm = p_comm_work
      RETURN
    ENDIF

#ifndef NOMPI
    ! Default if processor set is not split

    p_patch(:)%comm   = p_comm_work
    p_patch(:)%rank   = p_pe_work

    proc_split = .FALSE.

    IF(p_patch(1)%n_childdom <= 1) RETURN ! No splitting for 0 or 1 childs

    ! Check if the processor set is split for childs of root

    ALLOCATE(patch_no(0:p_n_work-1))

    n_proc_total = 0
    patch_no(:) = 0
    DO jc = 1, p_patch(1)%n_childdom
      jgc = p_patch(1)%child_id(jc)
      n_proc_total = n_proc_total + p_patch(jgc)%n_proc
      patch_no(p_patch(jgc)%proc0 : p_patch(jgc)%proc0+p_patch(jgc)%n_proc-1) = jc
    ENDDO

    ! if any processor has no patch assigned, this is an error

    IF(ANY(patch_no == 0)) &
      CALL finish('set_patch_communicators','Unknown patch split mode (1)')

    IF(n_proc_total == p_n_work) THEN

      proc_split = .TRUE.

      ! Split communicator among childs of root patch

      CALL MPI_Comm_split(p_comm_work, patch_no(p_pe_work), p_pe_work, comm, mpierr)

      ! Set comm and rank for childs of root patch

      DO jc = 1, p_patch(1)%n_childdom
        jgc = p_patch(1)%child_id(jc)
        IF(patch_no(p_pe_work) == jc) THEN
          p_patch(jgc)%comm = comm
          CALL MPI_Comm_rank(comm, p_patch(jgc)%rank, mpierr)
        ELSE
          p_patch(jgc)%comm = MPI_COMM_NULL ! We should never use this comm
          p_patch(jgc)%rank = -1
        ENDIF
      ENDDO

      ! and for deeper level descandants

      DO jg = 2, n_dom

        jgp = p_patch(jg)%parent_id

        IF(jgp /= 1) THEN
          p_patch(jg)%comm   = p_patch(jgp)%comm
          p_patch(jg)%rank   = p_patch(jgp)%rank
        ENDIF

      ENDDO

    ELSEIF(n_proc_total /= p_patch(1)%n_childdom * p_n_work) THEN

      ! If the processor is not split and n_proc_totalindicates that not
      ! every processor is working on every patch, this is an error
      CALL finish('set_patch_communicators','Unknown patch split mode (2)')

    ENDIF

    DEALLOCATE(patch_no)
#endif

  END SUBROUTINE set_patch_communicators

  !-----------------------------------------------------------------------------
  ! sets communication patterns and parent-child relationships.
  !
  SUBROUTINE complete_parallel_setup(patch, is_ocean_decomposition)

    TYPE(t_patch), INTENT(INOUT) :: patch(n_dom_start:)
    LOGICAL, INTENT(IN) :: is_ocean_decomposition

    INTEGER :: jg, jgp

    IF (is_ocean_decomposition .AND. (n_dom > n_dom_start)) &
      CALL finish('finalize_decomposition_oce', &
        &         'functionality with local parent patch not implemented')

    DO jg = n_dom_start, n_dom

      jgp = patch(jg)%parent_id

      ! Set communication patterns for boundary exchange
      CALL set_comm_pat_bound_exch(patch(jg))

      ! Set communication patterns for gathering on proc 0
      CALL set_comm_pat_gather(patch(jg))

      IF(jg == n_dom_start) THEN
        ! parent_idx/blk is set to 0 since it just doesn't exist,
        patch(jg)%cells%parent_idx = 0
        patch(jg)%cells%parent_blk = 0
        patch(jg)%edges%parent_idx = 0
        patch(jg)%edges%parent_blk = 0

        ! For parallel runs, child_idx/blk is set to 0 since it makes
        ! sense only on the local parent
        IF (.NOT. my_process_is_mpi_parallel() .OR. &
            is_ocean_decomposition) THEN
          patch(jg)%cells%child_idx  = 0
          patch(jg)%cells%child_blk  = 0
          patch(jg)%edges%child_idx  = 0
          patch(jg)%edges%child_blk  = 0
        END IF

      ELSE

        ! Note: The following call is deprecated and will be removed.
        !
        ! CALL setup_comm_cpy_interpolation(patch(jg), patch(jgp))

        CALL setup_comm_grf_interpolation(patch(jg), patch(jgp))
        CALL setup_comm_ubc_interpolation(patch(jg), patch(jgp))

        CALL set_comm_pat_bound_exch(p_patch_local_parent(jg))
        CALL set_comm_pat_gather(p_patch_local_parent(jg))

        CALL set_parent_child_relations(p_patch_local_parent(jg), patch(jg))

        CALL set_glb_loc_comm(patch(jgp), p_patch_local_parent(jg), &
          &                   patch(jg)%parent_child_index)
      ENDIF

    ENDDO

    CALL set_patch_communicators(patch)

  END SUBROUTINE complete_parallel_setup

  !-----------------------------------------------------------------------------

  SUBROUTINE finalize_decomposition(patch, is_ocean_decomposition)

    TYPE(t_patch), INTENT(INOUT) :: patch(n_dom_start:)
    LOGICAL, INTENT(IN) :: is_ocean_decomposition

    INTEGER :: jg

    IF (is_ocean_decomposition .AND. (n_dom > n_dom_start)) &
      CALL finish('finalize_decomposition_oce', &
        &         'functionality with local parent patch not implemented')

    ! Remap indices in patches and local parents

    DO jg = n_dom_start, n_dom

      CALL remap_patch_indices(patch(jg))

      IF(jg>n_dom_start) THEN
        CALL remap_patch_indices(p_patch_local_parent(jg))
      ENDIF

    ENDDO

  END SUBROUTINE finalize_decomposition

  !-------------------------------------------------------------------------------------------------
  !>
  !! Remaps negative index entries to valid ones
  !!

  SUBROUTINE remap_patch_indices(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch

    INTEGER :: i, j, jb, jl, ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, irl0, irl1, irl2, irl3

    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,p_patch%cell_type

!CDIR IEXPAND
        CALL remap_index(p_patch%cells%decomp_info%glb2loc_index, &
          & p_patch%cells%neighbor_idx(jl,jb,i),    &
          & p_patch%cells%neighbor_blk(jl,jb,i))

        ! edge_idx and vertex_idx should not need a remap !!!
        ! This is only left here if there should change something in the decomposition
!CDIR IEXPAND
        CALL remap_index(p_patch%edges%decomp_info%glb2loc_index, &
          & p_patch%cells%edge_idx(jl,jb,i),        &
          & p_patch%cells%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%verts%decomp_info%glb2loc_index, &
          & p_patch%cells%vertex_idx(jl,jb,i),      &
          & p_patch%cells%vertex_blk(jl,jb,i))

      ENDDO
    ENDDO

    ! ensure that cells%neighbor_idx lies in the correct grid row along the lateral boundary
    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (p_patch%cells%refin_ctrl(jl,jb) > 0 .AND.             &
          p_patch%cells%refin_ctrl(jl,jb) <= grf_bdywidth_c) THEN

        ilc1 = p_patch%cells%neighbor_idx(jl,jb,1)
        ibc1 = p_patch%cells%neighbor_blk(jl,jb,1)
        ilc2 = p_patch%cells%neighbor_idx(jl,jb,2)
        ibc2 = p_patch%cells%neighbor_blk(jl,jb,2)
        ilc3 = p_patch%cells%neighbor_idx(jl,jb,3)
        ibc3 = p_patch%cells%neighbor_blk(jl,jb,3)

        irl0 = p_patch%cells%refin_ctrl(jl,jb)

        IF (ilc1 > 0 .AND. ibc1 > 0) THEN
          irl1 = p_patch%cells%refin_ctrl(ilc1,ibc1)
          IF ( irl1 > 0 .AND. ABS(irl0 - irl1) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,1) = jl
            p_patch%cells%neighbor_blk(jl,jb,1) = jb
          ENDIF
        ENDIF

        IF (ilc2 > 0 .AND. ibc2 > 0) THEN
          irl2 = p_patch%cells%refin_ctrl(ilc2,ibc2)
          IF ( irl2 > 0 .AND. ABS(irl0 - irl2) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,2) = jl
            p_patch%cells%neighbor_blk(jl,jb,2) = jb
          ENDIF
        ENDIF

        IF (ilc3 > 0 .AND. ibc3 > 0) THEN
          irl3 = p_patch%cells%refin_ctrl(ilc3,ibc3)
          IF ( irl3 > 0 .AND. ABS(irl0 - irl3) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,3) = jl
            p_patch%cells%neighbor_blk(jl,jb,3) = jb
          ENDIF
        ENDIF

      ENDIF
    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (p_patch%edges%refin_ctrl(jl,jb) /= 1) THEN
        DO i=1,2
!CDIR IEXPAND
          CALL remap_index(p_patch%cells%decomp_info%glb2loc_index, &
            & p_patch%edges%cell_idx(jl,jb,i),        &
            & p_patch%edges%cell_blk(jl,jb,i))
        ENDDO
      ENDIF

      DO i=1,4
!CDIR IEXPAND
        CALL remap_index(p_patch%verts%decomp_info%glb2loc_index, &
          & p_patch%edges%vertex_idx(jl,jb,i),      &
          & p_patch%edges%vertex_blk(jl,jb,i))
      ENDDO

      DO i=1,4
!CDIR IEXPAND
        CALL remap_index(p_patch%edges%decomp_info%glb2loc_index, &
          & p_patch%edges%quad_idx(jl,jb,i),        &
          & p_patch%edges%quad_blk(jl,jb,i))
      ENDDO

      ! ensure that edges%cell_idx lies in the correct grid row along the lateral boundary
      IF (p_patch%edges%refin_ctrl(jl,jb) >= 2 .AND.            &
          p_patch%edges%refin_ctrl(jl,jb) <= grf_bdywidth_e) THEN

        ilc1 = p_patch%edges%cell_idx(jl,jb,1)
        ibc1 = p_patch%edges%cell_blk(jl,jb,1)
        ilc2 = p_patch%edges%cell_idx(jl,jb,2)
        ibc2 = p_patch%edges%cell_blk(jl,jb,2)
        irl1 = p_patch%cells%refin_ctrl(ilc1,ibc1)
        irl2 = p_patch%cells%refin_ctrl(ilc2,ibc2)

        IF (MOD(p_patch%edges%refin_ctrl(jl,jb),2)==0) THEN
          IF (irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            p_patch%edges%cell_idx(jl,jb,1) = p_patch%edges%cell_idx(jl,jb,2)
            p_patch%edges%cell_blk(jl,jb,1) = p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            p_patch%edges%cell_idx(jl,jb,2) = p_patch%edges%cell_idx(jl,jb,1)
            p_patch%edges%cell_blk(jl,jb,2) = p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ELSE
          IF (irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            p_patch%edges%cell_idx(jl,jb,1) = p_patch%edges%cell_idx(jl,jb,2)
            p_patch%edges%cell_blk(jl,jb,1) = p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            p_patch%edges%cell_idx(jl,jb,2) = p_patch%edges%cell_idx(jl,jb,1)
            p_patch%edges%cell_blk(jl,jb,2) = p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ENDIF

      ENDIF

    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,9-p_patch%cell_type

!CDIR IEXPAND
        CALL remap_index(p_patch%verts%decomp_info%glb2loc_index, &
          & p_patch%verts%neighbor_idx(jl,jb,i),    &
          & p_patch%verts%neighbor_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%edges%decomp_info%glb2loc_index, &
          & p_patch%verts%edge_idx(jl,jb,i),        &
          & p_patch%verts%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%cells%decomp_info%glb2loc_index, &
          & p_patch%verts%cell_idx(jl,jb,i),        &
          & p_patch%verts%cell_blk(jl,jb,i))
      ENDDO

    ENDDO

  END SUBROUTINE remap_patch_indices

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets parent_idx/blk in child and child_idx/blk in parent patches.

  SUBROUTINE set_parent_child_relations(p_pp, p_pc)

    TYPE(t_patch), INTENT(INOUT) :: p_pp   !> divided local parent patch
    TYPE(t_patch), INTENT(INOUT) :: p_pc   !> divided child patch

    INTEGER :: i, j, jl, jb, jc, jc_g, jp, jp_g

    ! Before this call, parent_idx/parent_blk and child_idx/child_blk still point to the global values.
    ! This is changed here.

    ! Attention:
    ! Only inner cells/edges get a valid child index,
    ! indexes for boundary cells/edges are not set.
    ! Therefore when these indexes are used, the code must assure that they are
    ! used only for inner cells/edges!
    ! The main reason for this is that - depending on the number of ghost rows -
    ! there are cells/edges in the parent boundary with missing childs (n_ghost_rows==1)

    ! Set child indices in parent ...

    ! ... cells

    DO j = 1, p_pp%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pp%cells%decomp_info%decomp_domain(jl,jb)>0) THEN
        p_pp%cells%child_idx(jl,jb,:) = 0
        p_pp%cells%child_blk(jl,jb,:) = 0
        CYCLE
      ENDIF

      DO i= 1, 4
        jc_g = idx_1d(p_pp%cells%child_idx(jl,jb,i),p_pp%cells%child_blk(jl,jb,i))
        IF(jc_g<1 .OR. jc_g>p_pc%n_patch_cells_g) &
          & CALL finish('set_parent_child_relations','Invalid cell child index in global parent')
        jc = get_local_index(p_pc%cells%decomp_info%glb2loc_index, jc_g)
        IF(jc <= 0) &
          & CALL finish('set_parent_child_relations','cell child index outside child domain')
        p_pp%cells%child_blk(jl,jb,i) = blk_no(jc)
        p_pp%cells%child_idx(jl,jb,i) = idx_no(jc)
      ENDDO

    ENDDO

    ! ... edges

    DO j = 1, p_pp%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pp%edges%decomp_info%decomp_domain(jl,jb)>1) THEN
        p_pp%edges%child_idx(jl,jb,:) = 0
        p_pp%edges%child_blk(jl,jb,:) = 0
        CYCLE ! only inner edges get a valid parent index
      ENDIF

      DO i= 1, 4

        IF(i==4 .AND. p_pp%edges%refin_ctrl(jl,jb) == -1) THEN
          p_pp%edges%child_blk(jl,jb,i) = blk_no(0)
          p_pp%edges%child_idx(jl,jb,i) = idx_no(0)
          CYCLE
        ENDIF

        jc_g = idx_1d(p_pp%edges%child_idx(jl,jb,i),p_pp%edges%child_blk(jl,jb,i))
        jc = get_valid_local_index(p_pc%edges%decomp_info%glb2loc_index, &
          &                        jc_g, .TRUE.)
        IF(jc == 0) &
          & CALL finish('set_parent_child_relations','edge child index outside child domain')
        p_pp%edges%child_blk(jl,jb,i) = blk_no(jc)
        p_pp%edges%child_idx(jl,jb,i) = SIGN(idx_no(jc),jc_g)
      ENDDO

      p_pp%edges%child_id(jl,jb) = p_pc%id

    ENDDO

    ! Set parent indices in child ...

    ! ... cells

    DO j = 1, p_pc%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jp_g = idx_1d(p_pc%cells%parent_idx(jl,jb),p_pc%cells%parent_blk(jl,jb))
      IF(jp_g<1 .OR. jp_g>p_pp%n_patch_cells_g) &
        & CALL finish('set_parent_child_relations','Inv. cell parent index in global child')

      jp = get_local_index(p_pp%cells%decomp_info%glb2loc_index, jp_g)
      IF(jp <= 0) THEN
        p_pc%cells%parent_blk(jl,jb) = 0
        p_pc%cells%parent_idx(jl,jb) = 0
      ELSE
        p_pc%cells%parent_blk(jl,jb) = blk_no(jp)
        p_pc%cells%parent_idx(jl,jb) = idx_no(jp)
      ENDIF

    ENDDO

    ! ... edges

    DO j = 1, p_pc%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jp_g = idx_1d(p_pc%edges%parent_idx(jl,jb),p_pc%edges%parent_blk(jl,jb))
      IF(jp_g<1 .OR. jp_g>p_pp%n_patch_edges_g) &
        & CALL finish('set_parent_child_relations','Inv. edge parent index in global child')

      jp = get_local_index(p_pp%edges%decomp_info%glb2loc_index, jp_g)
      IF(jp <= 0) THEN
        p_pc%edges%parent_blk(jl,jb) = 0
        p_pc%edges%parent_idx(jl,jb) = 0
      ELSE
        p_pc%edges%parent_blk(jl,jb) = blk_no(jp)
        p_pc%edges%parent_idx(jl,jb) = idx_no(jp)
      ENDIF

    ENDDO

    ! Although this is not really necessary, we set the child index in child
    ! and the parent index in parent to 0 since these have no significance
    ! in the parallel code (and must not be used as they are).

    IF (my_process_is_mpi_parallel()) THEN
      p_pc%cells%child_idx  = 0
      p_pc%cells%child_blk  = 0
      p_pp%cells%parent_idx = 0
      p_pp%cells%parent_blk = 0

      p_pc%edges%child_idx  = 0
      p_pc%edges%child_blk  = 0
      p_pp%edges%parent_idx = 0
      p_pp%edges%parent_blk = 0
    END IF

  END SUBROUTINE set_parent_child_relations


  !-------------------------------------------------------------------------------------------------
  !> Sets the communication patterns for boundary exchange

  SUBROUTINE set_comm_pat_bound_exch(p)

    TYPE(t_patch), INTENT(INOUT):: p

    INTEGER, ALLOCATABLE :: owner_c(:), owner_e(:), owner_v(:)
    INTEGER :: jc

    ALLOCATE(owner_c(p%n_patch_cells))
    ALLOCATE(owner_e(p%n_patch_edges))
    ALLOCATE(owner_v(p%n_patch_verts))

    ! Set the owner arrays for cells/edges/verts which have to be transferred.
    ! The following setting always transfers edges/verts if they are owned
    ! (according to the MIN/MAX PE setting above) by another PE, which implies
    ! that boundary edges/verts (with flag2_e/flag2_v = 1) can be excluded from
    ! prognostic computations if they are followed by a synchronization call.

    owner_c(:) = MERGE(p%cells%decomp_info%owner_local(:), -1, &
      &                p%cells%decomp_info%owner_local(:) /= p_pe_work)

    owner_e(:) = MERGE(p%edges%decomp_info%owner_local(:), -1, &
      &                p%edges%decomp_info%owner_local(:) /= p_pe_work)

    owner_v(:) = MERGE(p%verts%decomp_info%owner_local(:), -1, &
      &                p%verts%decomp_info%owner_local(:) /= p_pe_work)

    ! Set communication patterns for boundary exchange
    CALL setup_comm_pattern(p%n_patch_cells, owner_c, p%cells%decomp_info%glb_index, &
      & p%cells%decomp_info%glb2loc_index, p%comm_pat_c)

    CALL setup_comm_pattern(p%n_patch_edges, owner_e, p%edges%decomp_info%glb_index, &
      & p%edges%decomp_info%glb2loc_index, p%comm_pat_e)

    CALL setup_comm_pattern(p%n_patch_verts, owner_v, p%verts%decomp_info%glb_index, &
      & p%verts%decomp_info%glb2loc_index, p%comm_pat_v)

    DEALLOCATE(owner_e, owner_v)

    ! Set reduced communication pattern containing only level-1 halo cells (immediate neighbors)
    jc = idx_1d(p%cells%end_idx(min_rlcell_int-1,1), &
                p%cells%end_blk(min_rlcell_int-1,1))
    CALL setup_comm_pattern(jc, owner_c(1:jc), p%cells%decomp_info%glb_index, &
      & p%cells%decomp_info%glb2loc_index, p%comm_pat_c1)

    DEALLOCATE(owner_c)

  END SUBROUTINE set_comm_pat_bound_exch

  !-------------------------------------------------------------------------------------------------
  !> Sets the gather communication patterns of a patch

  SUBROUTINE set_comm_pat_gather(p)

    TYPE(t_patch), INTENT(INOUT):: p

    CALL setup_comm_gather_pattern( &
      p%cells%decomp_info%glb2loc_index%global_size, &
      p%cells%decomp_info%owner_local, &
      p%cells%decomp_info%glb_index, p%comm_pat_gather_c)
    CALL setup_comm_gather_pattern( &
      p%verts%decomp_info%glb2loc_index%global_size, &
      p%verts%decomp_info%owner_local, &
      p%verts%decomp_info%glb_index, p%comm_pat_gather_v)
    CALL setup_comm_gather_pattern( &
      p%edges%decomp_info%glb2loc_index%global_size, &
      p%edges%decomp_info%owner_local, &
      p%edges%decomp_info%glb_index, p%comm_pat_gather_e)

  END SUBROUTINE set_comm_pat_gather

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets up communication patterns between global and local parent patches

  SUBROUTINE set_glb_loc_comm(p_pglb, p_ploc, i_chidx)
    TYPE(t_patch), INTENT(IN)    :: p_pglb !> global parent
    TYPE(t_patch), INTENT(INOUT) :: p_ploc !> local parent
    INTEGER, INTENT(IN) :: i_chidx

    INTEGER, ALLOCATABLE :: owner(:)
    LOGICAL, ALLOCATABLE :: mask(:)
    INTEGER :: j, je, icid, jb, jl

    ! Please note:
    ! For creating communication patterns for different amount of data to be transferred
    ! (e.g. only the boundary interpolation zone), create a new communication pattern
    ! by copying the code below and adjusting the limits in the calculation of js/je.

    !-----------------------------------------------------------------------------------------------

    ! child ID
    icid = p_pglb%child_id(i_chidx)

    ! Communication global -> local
    ! Only one pattern is set which can be used everywhere since it doesn't matter
    ! if more cells/edges than needed are set in the local parent.
    ! Please note that only the inner area of the local patch is set and not the boundary cells
    ! (which are needed only for gradient calculation in the moment).

    ! ... cells

    ALLOCATE(owner(MAX(p_ploc%n_patch_cells, p_ploc%n_patch_edges, &
      &                p_pglb%n_patch_cells, p_pglb%n_patch_edges)), &
      &      mask(MAX(p_ploc%n_patch_cells, p_ploc%n_patch_edges, &
      &                p_pglb%n_patch_cells, p_pglb%n_patch_edges)))

    je = idx_1d(p_ploc%cells%end_idx(min_rlcell_int,i_chidx), &
      &         p_ploc%cells%end_blk(min_rlcell_int,i_chidx))

    mask(1:p_ploc%n_patch_cells) = .FALSE.
    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_ploc%cells%child_id(jl,jb)   == icid .AND. &
        &       p_ploc%cells%refin_ctrl(jl,jb) <= grf_bdyintp_start_c
    ENDDO

    owner(1:p_ploc%n_patch_cells) = &
      UNPACK(dist_dir_get_owners(p_pglb%cells%decomp_info%owner_dist_dir, &
        &                        PACK(p_ploc%cells%decomp_info%glb_index(:), &
        &                             mask(1:p_ploc%n_patch_cells))), &
        &    mask(1:p_ploc%n_patch_cells), -1)

    CALL setup_comm_pattern(p_ploc%n_patch_cells, owner(1:p_ploc%n_patch_cells), &
      p_ploc%cells%decomp_info%glb_index, &
      p_pglb%cells%decomp_info%glb2loc_index, p_ploc%comm_pat_glb_to_loc_c)

    ! ... edges

    je = idx_1d(p_ploc%edges%end_idx(min_rledge_int,i_chidx), &
      &         p_ploc%edges%end_blk(min_rledge_int,i_chidx))

    mask(1:p_ploc%n_patch_edges) = .FALSE.
    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_ploc%edges%child_id(jl,jb) == icid .AND. &
        &       p_ploc%edges%refin_ctrl(jl,jb) <= grf_bdyintp_start_e
    ENDDO

    owner(1:p_ploc%n_patch_edges) = &
      UNPACK(dist_dir_get_owners(p_pglb%edges%decomp_info%owner_dist_dir, &
        &                        PACK(p_ploc%edges%decomp_info%glb_index(:), &
        &                             mask(1:p_ploc%n_patch_edges))), &
        &    mask(1:p_ploc%n_patch_edges), -1)

    CALL setup_comm_pattern(p_ploc%n_patch_edges, owner(1:p_ploc%n_patch_edges), &
      p_ploc%edges%decomp_info%glb_index,  &
      p_pglb%edges%decomp_info%glb2loc_index, p_ploc%comm_pat_glb_to_loc_e)

    !-----------------------------------------------------------------------------------------------

    ! Communication local -> global
    ! Here it might get necessary to have different patterns for different start levels
    ! of the copy to the global parent since we may not overwrite arbitrary values there.
    ! Currently only one pattern for feedback is needed (starting at grf_fbk_start_c/e).

    ! ... cells

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      je = idx_1d(p_pglb%cells%end_idx(min_rlcell,p_pglb%n_childdom), &
        &         p_pglb%cells%end_blk(min_rlcell,p_pglb%n_childdom))
    ELSE
      je = idx_1d(p_pglb%cells%end_idx(min_rlcell_int,p_pglb%n_childdom), &
        &         p_pglb%cells%end_blk(min_rlcell_int,p_pglb%n_childdom))
    ENDIF

    mask(1:p_pglb%n_patch_cells) = .FALSE.
    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_pglb%cells%child_id(jl,jb) == icid .AND. &
        &       p_pglb%cells%refin_ctrl(jl,jb) <= grf_fbk_start_c .AND. &
        &       p_pglb%cells%refin_ctrl(jl,jb) >= min_rlcell_int
    ENDDO

    owner(1:p_pglb%n_patch_cells) = &
      UNPACK(dist_dir_get_owners(p_ploc%cells%decomp_info%owner_dist_dir, &
        &                        PACK(p_pglb%cells%decomp_info%glb_index(:), &
        &                             mask(1:p_pglb%n_patch_cells))), &
        &    mask(1:p_pglb%n_patch_cells), -1)

    CALL setup_comm_pattern(p_pglb%n_patch_cells, owner(1:p_pglb%n_patch_cells), &
      p_pglb%cells%decomp_info%glb_index, &
      p_ploc%cells%decomp_info%glb2loc_index, p_ploc%comm_pat_loc_to_glb_c_fbk)

    ! ... edges

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      je = idx_1d(p_pglb%edges%end_idx(min_rledge,i_chidx), &
        &         p_pglb%edges%end_blk(min_rledge,i_chidx))
    ELSE
      je = idx_1d(p_pglb%edges%end_idx(min_rledge_int,i_chidx), &
        &         p_pglb%edges%end_blk(min_rledge_int,i_chidx))
    ENDIF

    mask(1:p_pglb%n_patch_edges) = .FALSE.
    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_pglb%edges%child_id(jl,jb) == icid .AND. &
        &       p_pglb%edges%refin_ctrl(jl,jb) <= grf_fbk_start_e .AND. &
        &       p_pglb%edges%refin_ctrl(jl,jb) >= min_rledge_int
    ENDDO

    owner(1:p_pglb%n_patch_edges) = &
      UNPACK(dist_dir_get_owners(p_ploc%edges%decomp_info%owner_dist_dir, &
        &                        PACK(p_pglb%edges%decomp_info%glb_index(:), &
        &                             mask(1:p_pglb%n_patch_edges))), &
        &    mask(1:p_pglb%n_patch_edges), -1)

    CALL setup_comm_pattern(p_pglb%n_patch_edges, &
      owner(1:p_pglb%n_patch_edges), &
      p_pglb%edges%decomp_info%glb_index, &
      p_ploc%edges%decomp_info%glb2loc_index, p_ploc%comm_pat_loc_to_glb_e_fbk)

    DEALLOCATE(owner, mask)

   END SUBROUTINE set_glb_loc_comm

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for interpolation by direct copying.
  !!
  !! This routine sets up a communication pattern for interpolation by direct copying.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_cpy_interpolation(p_patch, p_parent_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch, p_parent_patch

    INTEGER :: j, jc, jb, jp, p_index_s, p_index_e, i_chidx
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !-----------------------------------------------------------------------

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Get start and end index of the GLOBAL parent cells as used in the interpolation

    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(grf_bdyintp_end_c,i_chidx), &
                       p_parent_patch%cells%end_blk(grf_bdyintp_end_c,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%decomp_info%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%decomp_info%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell

    ALLOCATE(glb_index(p_patch%n_patch_cells), &
      &      owner(p_patch%n_patch_cells))

    glb_index(:) = -1

    DO j = 1, p_patch%n_patch_cells
      jc = idx_no(j)
      jb = blk_no(j)
      jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
      IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
      glb_index(j) = jp
    ENDDO

    owner(:) = &
      UNPACK(dist_dir_get_owners(p_parent_patch%cells%decomp_info%owner_dist_dir, &
        &                        PACK(glb_index(:), glb_index(:) /= -1)), &
        &    glb_index(:) /= -1, -1)

    ! Set up communication pattern

    CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
      & p_parent_patch%cells%decomp_info%glb2loc_index, &
      & p_patch%comm_pat_interpolation_c)

    DEALLOCATE(owner, glb_index)

  END SUBROUTINE setup_comm_cpy_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for grf interpolation.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_grf_interpolation(p_patch, p_parent_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch, p_parent_patch

    INTEGER :: j, n, jc, je, jb, jp
    INTEGER :: num_send, num_recv, np, iss, ise, irs, ire
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !--------------------------------------------------------------------
    ! Cells

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    DO n = 1, 4

      glb_index(:) = -1

      ! Communication to nest boundary points includes halo points in order to save subsequent synchronization
      DO j = 1,p_patch%n_patch_cells
        jc = idx_no(j)
        jb = blk_no(j)
        IF (p_patch%cells%refin_ctrl(jc,jb) > 0 .AND. p_patch%cells%refin_ctrl(jc,jb) <= grf_bdywidth_c &
            .AND. p_patch%cells%pc_idx(jc,jb) == n) THEN
          jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
          glb_index(j) = jp
        ENDIF
      ENDDO

      owner(:) = &
        UNPACK(dist_dir_get_owners(p_parent_patch%cells%decomp_info%owner_dist_dir, &
          &                        PACK(glb_index(:), glb_index(:) /= -1)), &
          &    glb_index(:) /= -1, -1)

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
        & p_parent_patch%cells%decomp_info%glb2loc_index, &
        & p_patch%comm_pat_interpol_scal_grf(n))

    ENDDO

    ! Recompute send/recv processor lists for the cell-based communication patterns
    ! in order to be able to use the lists in the exchange routine
    ! (This is necessary even though the child cells of a given cell are always
    !  owned by the same PE because the halo cells are included in the communication pattern)

    num_send = 0
    num_recv = 0

    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_scal_grf(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_scal_grf(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_scal_grf(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_scal_grf(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DO n = 1, 4
      DEALLOCATE (p_patch%comm_pat_interpol_scal_grf(n)%pelist_send,   &
                  p_patch%comm_pat_interpol_scal_grf(n)%pelist_recv,   &
                  p_patch%comm_pat_interpol_scal_grf(n)%send_startidx, &
                  p_patch%comm_pat_interpol_scal_grf(n)%recv_startidx, &
                  p_patch%comm_pat_interpol_scal_grf(n)%send_count,    &
                  p_patch%comm_pat_interpol_scal_grf(n)%recv_count     )

      p_patch%comm_pat_interpol_scal_grf(n)%np_send = num_send
      p_patch%comm_pat_interpol_scal_grf(n)%np_recv = num_recv

      ALLOCATE (p_patch%comm_pat_interpol_scal_grf(n)%pelist_send(num_send),   &
                p_patch%comm_pat_interpol_scal_grf(n)%pelist_recv(num_recv),   &
                p_patch%comm_pat_interpol_scal_grf(n)%send_startidx(num_send), &
                p_patch%comm_pat_interpol_scal_grf(n)%recv_startidx(num_recv), &
                p_patch%comm_pat_interpol_scal_grf(n)%send_count(num_send),    &
                p_patch%comm_pat_interpol_scal_grf(n)%recv_count(num_recv)     )

      ! The startidx and count fields are not used for the nest interpolation
      ! communication, but the fields have to redimensioned here in order to avoid
      ! segmentation faults in dump/restore mode
      p_patch%comm_pat_interpol_scal_grf(n)%send_startidx(:) = 0
      p_patch%comm_pat_interpol_scal_grf(n)%recv_startidx(:) = 0
      p_patch%comm_pat_interpol_scal_grf(n)%send_count(:) = 0
      p_patch%comm_pat_interpol_scal_grf(n)%recv_count(:) = 0
    ENDDO

    num_send = 0
    num_recv = 0

    ! Now compute "envelope PE lists" for all communication patterns
    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_scal_grf(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_scal_grf(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_scal_grf(j)%pelist_send(num_send) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_scal_grf(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_scal_grf(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_scal_grf(j)%pelist_recv(num_recv) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DEALLOCATE(owner, glb_index)

    !--------------------------------------------------------------------
    ! Edges

    ! For our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(p_patch%n_patch_edges))
    ALLOCATE(owner(p_patch%n_patch_edges))

    DO n = 1, 4

      glb_index(:) = -1

      ! Communication to nest boundary points includes halo points in order to save subsequent synchronization
      DO j = 1,p_patch%n_patch_edges
        je = idx_no(j)
        jb = blk_no(j)
        IF (p_patch%edges%refin_ctrl(je,jb) > 0 .AND. p_patch%edges%refin_ctrl(je,jb) <= grf_bdywidth_e &
            .AND. p_patch%edges%pc_idx(je,jb) == n) THEN
          jp = idx_1d(p_patch%edges%parent_idx(je,jb),p_patch%edges%parent_blk(je,jb))
          glb_index(j) = jp
        ENDIF
      ENDDO

      owner(:) = &
        UNPACK(dist_dir_get_owners(p_parent_patch%edges%decomp_info%owner_dist_dir, &
          &                        PACK(glb_index(:), glb_index(:) /= -1)), &
          &    glb_index(:) /= -1, -1)

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_edges, owner, glb_index,  &
        & p_parent_patch%edges%decomp_info%glb2loc_index, &
        & p_patch%comm_pat_interpol_vec_grf(n))

    ENDDO

    ! Recompute send/recv processor lists for the edge-based communication patterns
    ! in order to be able to use the lists in the exchange routine

    num_send = 0
    num_recv = 0

    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_vec_grf(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_vec_grf(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_vec_grf(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_vec_grf(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DO n = 1, 4
      DEALLOCATE (p_patch%comm_pat_interpol_vec_grf(n)%pelist_send,   &
                  p_patch%comm_pat_interpol_vec_grf(n)%pelist_recv,   &
                  p_patch%comm_pat_interpol_vec_grf(n)%send_startidx, &
                  p_patch%comm_pat_interpol_vec_grf(n)%recv_startidx, &
                  p_patch%comm_pat_interpol_vec_grf(n)%send_count,    &
                  p_patch%comm_pat_interpol_vec_grf(n)%recv_count     )

      p_patch%comm_pat_interpol_vec_grf(n)%np_send = num_send
      p_patch%comm_pat_interpol_vec_grf(n)%np_recv = num_recv

      ALLOCATE (p_patch%comm_pat_interpol_vec_grf(n)%pelist_send(num_send),   &
                p_patch%comm_pat_interpol_vec_grf(n)%pelist_recv(num_recv),   &
                p_patch%comm_pat_interpol_vec_grf(n)%send_startidx(num_send), &
                p_patch%comm_pat_interpol_vec_grf(n)%recv_startidx(num_recv), &
                p_patch%comm_pat_interpol_vec_grf(n)%send_count(num_send),    &
                p_patch%comm_pat_interpol_vec_grf(n)%recv_count(num_recv)     )

      ! The startidx and count fields are not used for the nest interpolation
      ! communication, but the fields have to redimensioned here in order to avoid
      ! segmentation faults in dump/restore mode
      p_patch%comm_pat_interpol_vec_grf(n)%send_startidx(:) = 0
      p_patch%comm_pat_interpol_vec_grf(n)%recv_startidx(:) = 0
      p_patch%comm_pat_interpol_vec_grf(n)%send_count(:) = 0
      p_patch%comm_pat_interpol_vec_grf(n)%recv_count(:) = 0
    ENDDO

    num_send = 0
    num_recv = 0

    ! Now compute "envelope PE lists" for all communication patterns
    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_vec_grf(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_vec_grf(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_vec_grf(j)%pelist_send(num_send) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_vec_grf(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_vec_grf(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_vec_grf(j)%pelist_recv(num_recv) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DEALLOCATE(owner, glb_index)

  END SUBROUTINE setup_comm_grf_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for ubc interpolation.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_ubc_interpolation(p_patch, p_parent_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch, p_parent_patch

    INTEGER :: j, n, jc, je, jb, jp
    INTEGER :: num_send, num_recv, np, iss, ise, irs, ire
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !--------------------------------------------------------------------
    ! Cells

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    DO n = 1, 4

      glb_index(:) = -1

      DO j = 1,p_patch%n_patch_cells
        jc = idx_no(j)
        jb = blk_no(j)
        IF ((p_patch%cells%refin_ctrl(jc,jb) >= grf_bdywidth_c+1 .OR. p_patch%cells%refin_ctrl(jc,jb) <= 0) &
            .AND. p_patch%cells%pc_idx(jc,jb) == n) THEN
          jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
          glb_index(j) = jp
        ENDIF
      ENDDO

      owner(:) = &
        UNPACK(dist_dir_get_owners(p_parent_patch%cells%decomp_info%owner_dist_dir, &
          &                        PACK(glb_index(:), glb_index(:) /= -1)), &
          &    glb_index(:) /= -1, -1)

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
        & p_parent_patch%cells%decomp_info%glb2loc_index, &
        & p_patch%comm_pat_interpol_scal_ubc(n))

    ENDDO

    ! Recompute send/recv processor lists for the cell-based communication patterns
    ! in order to be able to use the lists in the exchange routine
    ! (This is necessary even though the child cells of a given cell are always
    !  owned by the same PE because the halo cells are included in the communication pattern)

    num_send = 0
    num_recv = 0

    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_scal_ubc(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_scal_ubc(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_scal_ubc(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_scal_ubc(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DO n = 1, 4
      DEALLOCATE (p_patch%comm_pat_interpol_scal_ubc(n)%pelist_send,   &
                  p_patch%comm_pat_interpol_scal_ubc(n)%pelist_recv,   &
                  p_patch%comm_pat_interpol_scal_ubc(n)%send_startidx, &
                  p_patch%comm_pat_interpol_scal_ubc(n)%recv_startidx, &
                  p_patch%comm_pat_interpol_scal_ubc(n)%send_count,    &
                  p_patch%comm_pat_interpol_scal_ubc(n)%recv_count     )

      p_patch%comm_pat_interpol_scal_ubc(n)%np_send = num_send
      p_patch%comm_pat_interpol_scal_ubc(n)%np_recv = num_recv

      ALLOCATE (p_patch%comm_pat_interpol_scal_ubc(n)%pelist_send(num_send),   &
                p_patch%comm_pat_interpol_scal_ubc(n)%pelist_recv(num_recv),   &
                p_patch%comm_pat_interpol_scal_ubc(n)%send_startidx(num_send), &
                p_patch%comm_pat_interpol_scal_ubc(n)%recv_startidx(num_recv), &
                p_patch%comm_pat_interpol_scal_ubc(n)%send_count(num_send),    &
                p_patch%comm_pat_interpol_scal_ubc(n)%recv_count(num_recv)     )

      ! The startidx and count fields are not used for the nest interpolation
      ! communication, but the fields have to redimensioned here in order to avoid
      ! segmentation faults in dump/restore mode
      p_patch%comm_pat_interpol_scal_ubc(n)%send_startidx(:) = 0
      p_patch%comm_pat_interpol_scal_ubc(n)%recv_startidx(:) = 0
      p_patch%comm_pat_interpol_scal_ubc(n)%send_count(:) = 0
      p_patch%comm_pat_interpol_scal_ubc(n)%recv_count(:) = 0
    ENDDO

    num_send = 0
    num_recv = 0

    ! Now compute "envelope PE lists" for all communication patterns
    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_scal_ubc(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_scal_ubc(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_scal_ubc(j)%pelist_send(num_send) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_scal_ubc(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_scal_ubc(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_scal_ubc(j)%pelist_recv(num_recv) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DEALLOCATE(owner, glb_index)

    !--------------------------------------------------------------------
    ! Edges

    ! For our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(p_patch%n_patch_edges))
    ALLOCATE(owner(p_patch%n_patch_edges))

    DO n = 1, 4

      glb_index(:) = -1

      DO j = 1,p_patch%n_patch_edges
        je = idx_no(j)
        jb = blk_no(j)
        IF ((p_patch%edges%refin_ctrl(je,jb) >= grf_bdywidth_e+1 .OR. p_patch%edges%refin_ctrl(je,jb) <= 0) &
            .AND. p_patch%edges%pc_idx(je,jb) == n) THEN
          jp = idx_1d(p_patch%edges%parent_idx(je,jb),p_patch%edges%parent_blk(je,jb))
          glb_index(j) = jp
        ENDIF
      ENDDO

      owner(:) = &
        UNPACK(dist_dir_get_owners(p_parent_patch%edges%decomp_info%owner_dist_dir, &
          &                        PACK(glb_index(:), glb_index(:) /= -1)), &
          &    glb_index(:) /= -1, -1)

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_edges, owner, glb_index,  &
        & p_parent_patch%edges%decomp_info%glb2loc_index, &
        & p_patch%comm_pat_interpol_vec_ubc(n))

    ENDDO

    ! Recompute send/recv processor lists for the edge-based communication patterns
    ! in order to be able to use the lists in the exchange routine

    num_send = 0
    num_recv = 0

    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_vec_ubc(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_vec_ubc(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_vec_ubc(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_vec_ubc(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DO n = 1, 4
      DEALLOCATE (p_patch%comm_pat_interpol_vec_ubc(n)%pelist_send,   &
                  p_patch%comm_pat_interpol_vec_ubc(n)%pelist_recv,   &
                  p_patch%comm_pat_interpol_vec_ubc(n)%send_startidx, &
                  p_patch%comm_pat_interpol_vec_ubc(n)%recv_startidx, &
                  p_patch%comm_pat_interpol_vec_ubc(n)%send_count,    &
                  p_patch%comm_pat_interpol_vec_ubc(n)%recv_count     )

      p_patch%comm_pat_interpol_vec_ubc(n)%np_send = num_send
      p_patch%comm_pat_interpol_vec_ubc(n)%np_recv = num_recv

      ALLOCATE (p_patch%comm_pat_interpol_vec_ubc(n)%pelist_send(num_send),   &
                p_patch%comm_pat_interpol_vec_ubc(n)%pelist_recv(num_recv),   &
                p_patch%comm_pat_interpol_vec_ubc(n)%send_startidx(num_send), &
                p_patch%comm_pat_interpol_vec_ubc(n)%recv_startidx(num_recv), &
                p_patch%comm_pat_interpol_vec_ubc(n)%send_count(num_send),    &
                p_patch%comm_pat_interpol_vec_ubc(n)%recv_count(num_recv)     )

      ! The startidx and count fields are not used for the nest interpolation
      ! communication, but the fields have to redimensioned here in order to avoid
      ! segmentation faults in dump/restore mode
      p_patch%comm_pat_interpol_vec_ubc(n)%send_startidx(:) = 0
      p_patch%comm_pat_interpol_vec_ubc(n)%recv_startidx(:) = 0
      p_patch%comm_pat_interpol_vec_ubc(n)%send_count(:) = 0
      p_patch%comm_pat_interpol_vec_ubc(n)%recv_count(:) = 0
    ENDDO

    num_send = 0
    num_recv = 0

    ! Now compute "envelope PE lists" for all communication patterns
    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, 4  ! loop over communication patterns
        iss = p_patch%comm_pat_interpol_vec_ubc(n)%send_limits(np)+1
        ise = p_patch%comm_pat_interpol_vec_ubc(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_vec_ubc(j)%pelist_send(num_send) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, 4  ! loop over communication patterns
        irs = p_patch%comm_pat_interpol_vec_ubc(n)%recv_limits(np)+1
        ire = p_patch%comm_pat_interpol_vec_ubc(n)%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          DO j = 1, 4
            p_patch%comm_pat_interpol_vec_ubc(j)%pelist_recv(num_recv) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    DEALLOCATE(owner, glb_index)

  END SUBROUTINE setup_comm_ubc_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a the physical patches
  !! It works on the divided patch state, so it can be used after a restore also
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2011
  !! Revised version by Moritz Hanke, Nov 2013

  SUBROUTINE setup_phys_patches
    INTEGER :: jp, jg, i
    CHARACTER(LEN=*), PARAMETER :: routine = 'setup_phys_patches'

    INTEGER :: n_patch_cells(max_phys_dom), n_patch_verts(max_phys_dom), &
      &        n_patch_edges(max_phys_dom)
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_cells(max_phys_dom), &
      &                            comm_pat_gather_verts(max_phys_dom), &
      &                            comm_pat_gather_edges(max_phys_dom)

    p_phys_patch(:)%logical_id = -1
    n_patch_cells(:) = 0
    n_patch_verts(:) = 0
    n_patch_edges(:) = 0

    DO jg = 1, n_dom
      CALL setup_phys_patches_cve(p_patch(jg)%n_patch_cells_g, &
        &                         p_patch(jg)%n_patch_cells, &
        &                         p_patch(jg)%cells%decomp_info, &
        &                         p_patch(jg)%cells%phys_id, jg, &
        &                         .TRUE., n_patch_cells(:), &
        &                         comm_pat_gather_cells(:))
      CALL setup_phys_patches_cve(p_patch(jg)%n_patch_verts_g, &
        &                         p_patch(jg)%n_patch_verts, &
        &                         p_patch(jg)%verts%decomp_info, &
        &                         p_patch(jg)%verts%phys_id, jg, &
        &                         .FALSE., n_patch_verts(:), &
        &                         comm_pat_gather_verts(:))
      CALL setup_phys_patches_cve(p_patch(jg)%n_patch_edges_g, &
        &                         p_patch(jg)%n_patch_edges, &
        &                         p_patch(jg)%edges%decomp_info, &
        &                         p_patch(jg)%edges%phys_id, jg, &
        &                         .FALSE., n_patch_edges(:), &
        &                         comm_pat_gather_edges(:))
    ENDDO

    ! MoHa: remark:
    ! We cannot pass p_phys_patch(:)%n_patch_cells and 
    ! p_phys_patch(:)%comm_pat_gather_c directly to setup_phys_patches_cve
    ! because some compiler(intel and pgi) have problems with this type
    ! of array argument... Therefore, the temporary arrays n_patch_*(:) and
    ! comm_pat_gather_*(:) had to be used.
    DO i = 1, max_phys_dom
      p_phys_patch(i)%n_patch_cells = n_patch_cells(i)
      p_phys_patch(i)%comm_pat_gather_c = comm_pat_gather_cells(i)
      CALL delete_comm_gather_pattern(comm_pat_gather_cells(i))
      p_phys_patch(i)%n_patch_verts = n_patch_verts(i)
      p_phys_patch(i)%comm_pat_gather_v = comm_pat_gather_verts(i)
      CALL delete_comm_gather_pattern(comm_pat_gather_verts(i))
      p_phys_patch(i)%n_patch_edges = n_patch_edges(i)
      p_phys_patch(i)%comm_pat_gather_e = comm_pat_gather_edges(i)
      CALL delete_comm_gather_pattern(comm_pat_gather_edges(i))
    ENDDO

    ! Get total number of physical patches

    DO jp = max_phys_dom, 1, -1
      IF(p_phys_patch(jp)%logical_id > 0) EXIT
    ENDDO

    n_phys_dom = jp

    ! Print info and check if there are no unused physical patches

    DO jp = 1, n_phys_dom
      WRITE(message_text,'(a,i4,a,i4)') &
        & 'Physical domain ',jp,' belongs to logical domain ',p_phys_patch(jp)%logical_id
      CALL message (routine, message_text)
      IF(p_phys_patch(jp)%logical_id < 0) THEN
        WRITE(message_text,'(a,i4)') 'Missing physical domain # ',jp
        CALL finish (routine, message_text)
      ENDIF
    ENDDO

  END SUBROUTINE setup_phys_patches

  SUBROUTINE setup_phys_patches_cve(n_g, n, decomp_info, phys_id, &
    &                               curr_patch_idx, set_logical_id, &
    &                               n_patch_cve, comm_pat_gather)
    INTEGER, INTENT(IN) :: n_g, n
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info
    INTEGER, INTENT(IN) :: phys_id(:,:)
    INTEGER, INTENT(IN) :: curr_patch_idx
    LOGICAL, INTENT(IN) :: set_logical_id
    INTEGER, INTENT(INOUT) :: n_patch_cve(:)
    TYPE(t_comm_gather_pattern), INTENT(INOUT) :: comm_pat_gather(:)

    INTEGER, ALLOCATABLE :: owner_local(:)
    INTEGER :: temp_n_patch_cve(max_phys_dom)
    INTEGER :: i, ip

    ALLOCATE(owner_local(n))

    temp_n_patch_cve(:) = 0

    ! Fill with own values of phys_id

    DO i = 1, n
      ip = phys_id(idx_no(i),blk_no(i))
      IF (ip < 1 .OR. ip > max_phys_dom) &
        CALL finish("setup_phys_patches_cve", "invalid phys_id")
      IF (decomp_info%owner_local(i) == p_pe_work) &
        temp_n_patch_cve(ip) = temp_n_patch_cve(ip) + 1
    END DO

    temp_n_patch_cve(:) = p_sum(temp_n_patch_cve(:), p_comm_work)

    IF (set_logical_id) THEN

      ! Check if no other patch uses the same phys_id
      IF (ANY((p_phys_patch(1:max_phys_dom)%logical_id /= -1) .AND. &
        &     (temp_n_patch_cve(:) /= 0))) &
        CALL finish("setup_phys_patches_cve", &
          &         "invalid value in p_phys_patch(:)%logical_id")

      WHERE(temp_n_patch_cve(:) /= 0) &
        p_phys_patch(1:max_phys_dom)%logical_id = curr_patch_idx
    ELSE
      ! Check if no other patch uses the same phys_id
      IF (ANY((p_phys_patch(1:max_phys_dom)%logical_id /= curr_patch_idx) &
        &     .AND. (temp_n_patch_cve(:) /= 0))) &
        CALL finish("setup_phys_patches_cve", &
          &         "invalid value in p_phys_patch(:)%logical_id")
    END IF

    DO ip = 1, max_phys_dom

      IF (p_phys_patch(ip)%logical_id /= curr_patch_idx) CYCLE

      n_patch_cve(ip) = temp_n_patch_cve(ip)

      DO i = 1, n
        owner_local(i) = &
          MERGE(decomp_info%owner_local(i), -1, &
            &   phys_id(idx_no(i),blk_no(i)) == ip)
      END DO

      CALL setup_comm_gather_pattern(n_g, owner_local(:), &
        &                            decomp_info%glb_index(:), &
        &                            comm_pat_gather(ip))
    END DO

    DEALLOCATE(owner_local)

  END SUBROUTINE setup_phys_patches_cve

  !-------------------------------------------------------------------------
  !>
  !!               Calculates local line/block indices l_idx, l_blk
  !!               from global line/block indices g_idx, g_blk
  !!               using the mapping in decomp_info
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!

  !-------------------------------------------------------------------------
  !>
  !! Maps indices which point outside domain (returned as <0 by get_local_idx_blk)
  !! to valid ones which are from the same set of neighbors for the given
  !! cell/edge/vert
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Oct 2011
  !!
  SUBROUTINE remap_index(glb2loc_index, l_idx, l_blk)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(inout) :: l_idx, l_blk

    INTEGER :: j_l

    IF(l_idx>=0) RETURN ! Nothing to do

    j_l = MAX(get_valid_local_index(glb2loc_index, &
      &                             idx_1d(-l_idx, l_blk)), 1)
    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE remap_index

  !-----------------------------------------------------------------------------
END MODULE mo_complete_subdivision

