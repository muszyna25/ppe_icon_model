!>
!!
!! @par Revision History
!!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_build_decomposition

  USE mo_complete_subdivision, ONLY: finalize_decomposition,                  &
    &                                copy_processor_splitting,                &
    &                                complete_parallel_setup
  USE mo_setup_subdivision,    ONLY: decompose_domain
  USE mo_sync,                 ONLY: disable_sync_checks, enable_sync_checks, &
    &                                decomposition_statistics
  USE mo_grid_config,          ONLY: n_dom, n_dom_start
  USE mo_mpi,                  ONLY: my_process_is_mpi_parallel, p_pe_work,   &
    &                                p_comm_work, p_comm_work_test
  USE mo_loopindices,          ONLY: get_indices_e
  USE mo_model_domain,         ONLY: p_patch, t_pre_patch, t_patch_3d,        &
    &                                t_patch, p_patch_local_parent
  USE mo_reshuffle,            ONLY: reshuffle
  USE mo_model_domimp_patches, ONLY: reorder_patch_refin_ctrl,                &
    &                                import_pre_patches, complete_patches, set_parent_loc_idx
  USE mo_parallel_config,      ONLY: p_test_run, l_test_openmp, num_io_procs, division_method
  USE mo_impl_constants,       ONLY: success, max_dom
  USE mo_exception,            ONLY: finish, message, message_text, get_filename_noext
  USE mo_communication,        ONLY: blk_no, idx_no, idx_1d
  USE mo_util_string,          ONLY: int2string

  IMPLICIT NONE
  
  PUBLIC :: build_decomposition

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_build_decomposition'
  
CONTAINS

  !> Main routine for creating the domain decomposition (together with
  !> communication patterns etc.)
  !
  !  @note This routine is called for both: The ocean model and the
  !        atmo_model.
  !
  !  @author F. Prill, DWD (2013-08-06)
  !
  SUBROUTINE build_decomposition(num_lev,nshift,&
    &                            is_ocean_decomposition, patch_3d)
    
    INTEGER, INTENT(in)                 :: num_lev(max_dom),                &
      &                                    nshift(max_dom)
    LOGICAL, INTENT(in)                 :: is_ocean_decomposition
    TYPE(t_patch_3d), POINTER, OPTIONAL :: patch_3d
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'build_decomposition'
    TYPE(t_pre_patch), ALLOCATABLE :: p_patch_pre(:)
    INTEGER                        :: error_status, jg, jgp
    !> If .true., read fields related to grid refinement from separate  grid files:
    LOGICAL                        :: lsep_grfinfo
    
    ! Check p_patch allocation status
    
    IF ( ALLOCATED(p_patch)) THEN
      CALL finish(TRIM(routine), 'p_patch already allocated')
    END IF
    
    ! Allocate p_patch array to start p_patch construction.
    !
    ! At the same time, we allocate the "p_patch_local_parent" which
    ! is the local portion of each p_patch's parent grid.
    ALLOCATE(p_patch             (n_dom_start  :n_dom), &
      & p_patch_local_parent(n_dom_start+1:n_dom), &
      & stat=error_status)
    IF (error_status/=success) CALL finish(TRIM(routine), 'allocation of p_patch failed')
    
    ! --------------------------
    ! Work PEs subdivide patches
    ! --------------------------
                  
    ! compute domain decomposition on-the-fly        
    ALLOCATE(p_patch_pre(n_dom_start:n_dom))
    CALL import_pre_patches(p_patch_pre,num_lev,nshift,lsep_grfinfo)
    ! use internal domain decomposition algorithm
    CALL decompose_domain(p_patch, p_patch_pre)
    DEALLOCATE(p_patch_pre)

    ! set local parent indices 
    DO jg = n_dom_start+1, n_dom
      CALL set_parent_loc_idx(p_patch_local_parent(jg), p_patch(jg))
    END DO

    ! compute the child/edge indices on parent patch from the
    ! corresponding information on parent indices stored in the child
    ! patch
    DO jg = n_dom_start,n_dom
      ! set child indices of patch:
      jgp = p_patch(jg)%parent_id
      IF (jgp >= n_dom_start) THEN
        CALL set_child_indices("patch "//TRIM(int2string(jgp))//" -> "//int2string(jg), &
          &                    p_patch(jg), p_patch(jgp), is_local_parent=.FALSE.)
      END IF
    END DO

    DO jg = n_dom_start, n_dom
      ! set child indices of local parent patch, but don't modify the
      ! "pc_idx" fields of the child patch:
      jgp = p_patch(jg)%parent_id
      IF (jgp >= n_dom_start) THEN
        CALL set_child_indices("loc par patch "//TRIM(int2string(jgp))//" -> "//int2string(jg), &
          &                    p_patch(jg), p_patch_local_parent(jg), is_local_parent=.TRUE.)
      END IF
    END DO

    ! Complete information which is not yet read or calculated
    CALL complete_patches( p_patch, is_ocean_decomposition, lsep_grfinfo)

    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    !--------------------------------------------------------------------------------

    CALL finalize_decomposition(p_patch, is_ocean_decomposition)

    ! reorder patch_local_parents according to their refin_ctrl flags
    DO jg = n_dom_start+1, n_dom
      CALL reorder_patch_refin_ctrl(p_patch_local_parent(jg), p_patch(jg))
    ENDDO

    ! computes communication patterns (done after reordering)
    CALL complete_parallel_setup(p_patch, is_ocean_decomposition)

    IF(.NOT.p_test_run .AND. my_process_is_mpi_parallel()) THEN ! the call below hangs in test mode
      ! Print diagnostic information about domain decomposition
      DO jg = 1, n_dom
        CALL decomposition_statistics(p_patch(jg))
      ENDDO
    ENDIF

    ! set the horizontal attribute pointer of the 3D p_patch
    IF (PRESENT(patch_3d)) THEN
      ALLOCATE(patch_3d, stat=error_status)
      IF (error_status/=success) THEN
        CALL finish(TRIM(routine), 'allocation of patch_3D failed')
      ENDIF
      patch_3d%p_patch_2d => p_patch
    END IF
    
  END SUBROUTINE build_decomposition
  !----------------------------------------------------------------------------


  ! -----------------------------------------------------------------
  !
  ! This subroutine computes the child/edge indices of a given patch
  ! "p_p" from the corresponding information on parent indices stored
  ! in the child patch "p_c". This involves parallel comm., since the
  ! child cells and their parent cells may "live" on different PEs.
  !
  ! We also compute the pc_idx here, i.e. the relative ordering of the
  ! child cells/edges in the parent cell/edge. 
  !
  !          *---------------o--------------*
  !           \__   4     _/  \_    2    __/
  !              \_     _/   3  \_    __/
  !                \__ /__________\ _/
  !                   o            o  
  !                     \__  1  __/
  !                        \_ _/
  !                          *
  !
  ! input:   n_patch_cells, n_patch_cells_g  (on child and parent patch)
  !          n_patch_edges, n_patch_edges_g  (on child and parent patch)
  !          parent_glb_idx/blk  on child patch            (cells and edges)
  !          glb_index           on child and parent patch (cells and edges)
  !          edge_idx/blk        on child patch
  !          decomp_domain, i.e. cell ownership
  !
  ! output:
  !
  ! 06/2016 : F. Prill, DWD
  !
  SUBROUTINE set_child_indices(description, p_c, p_p, is_local_parent)
    CHARACTER(LEN=*), INTENT(IN) :: description !< description string (for debugging purposes)
    TYPE(t_patch), TARGET, INTENT(inout) :: p_c ! child patch
    TYPE(t_patch), TARGET, INTENT(inout) :: p_p ! parent patch
    LOGICAL, INTENT(IN) :: is_local_parent
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":set_child_indices"
    INTEGER              :: i, j, jc_c, jb_c, communicator, i1, i2, iidx,                    &
      &                     tmp, jc_e, jb_e, ordered(4), iglb
    INTEGER, ALLOCATABLE :: in_data(:), dst_idx(:), out_data(:,:), out_pc_data(:,:),         &
      &                     out_data34(:,:), out_data_e(:,:), out_count(:,:),                &
      &                     out_data_e3(:,:), out_data_c3(:,:), out_child_id(:,:)
    LOGICAL              :: lfound, set_pc_idx
    INTEGER, POINTER     :: parent_idx_c(:,:), parent_blk_c(:,:),                            &
      &                     parent_idx_e(:,:), parent_blk_e(:,:),                            &
      &                     lparent_glb_idx_c(:,:), lparent_glb_blk_c(:,:),                  &
      &                     lparent_glb_idx_e(:,:), lparent_glb_blk_e(:,:)


    IF (p_pe_work == 0) THEN
      WRITE (0,*) "set_child_indices: ", TRIM(description), " - Enter."
    END IF

    communicator = p_comm_work
    IF(p_test_run)  communicator = p_comm_work_test

    parent_idx_c => p_c%cells%parent_glb_idx
    parent_blk_c => p_c%cells%parent_glb_blk
    parent_idx_e => p_c%edges%parent_glb_idx
    parent_blk_e => p_c%edges%parent_glb_blk

    IF (is_local_parent) THEN
      ! if "p_p" is a local parent patch, we need to create a mapping
      ! between child patch cells and "p_p" parent cells, but as
      ! global indices (this is why we cannot use
      ! p_c%parent_loc_idx/blk).
      ALLOCATE(lparent_glb_idx_c, source=p_c%cells%parent_loc_idx)
      ALLOCATE(lparent_glb_blk_c, source=p_c%cells%parent_loc_blk)
      DO j = 1, p_c%n_patch_cells
        jc_c = idx_no(j)
        jb_c = blk_no(j)
        ! (we now exploit the fact that p_c and p_p have the same domain decomp.:)
        iglb = p_p%cells%decomp_info%glb_index(idx_1d(p_c%cells%parent_loc_idx(jc_c,jb_c), &
          &                                           p_c%cells%parent_loc_blk(jc_c,jb_c)))
        lparent_glb_idx_c(jc_c,jb_c) = idx_no(iglb)
        lparent_glb_blk_c(jc_c,jb_c) = blk_no(iglb)
      END DO
      parent_idx_c => lparent_glb_idx_c
      parent_blk_c => lparent_glb_blk_c

      ALLOCATE(lparent_glb_idx_e, source=p_c%edges%parent_loc_idx)
      ALLOCATE(lparent_glb_blk_e, source=p_c%edges%parent_loc_blk)
      DO j = 1, p_c%n_patch_edges
        jc_e = idx_no(j)
        jb_e = blk_no(j)
        ! (we now exploit the fact that p_c and p_p have the same domain decomp.:)
        iglb = p_p%edges%decomp_info%glb_index(idx_1d(p_c%edges%parent_loc_idx(jc_e,jb_e), &
          &                                           p_c%edges%parent_loc_blk(jc_e,jb_e)))
        lparent_glb_idx_e(jc_e,jb_e) = idx_no(iglb)
        lparent_glb_blk_e(jc_e,jb_e) = blk_no(iglb)
      END DO
      parent_idx_e => lparent_glb_idx_e
      parent_blk_e => lparent_glb_blk_e
    END IF
    set_pc_idx = .NOT. is_local_parent

    ! -----------------------------------------------------------------
    ! --- create CELL indices -----------------------------------------
    ! -----------------------------------------------------------------

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- send the four child cells to the parent cell ----------------

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))

    ! get the parent cell indices:
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      dst_idx(j) = idx_1d(parent_idx_c(jc_c,jb_c), &
        &                 parent_blk_c(jc_c,jb_c))
    END DO

    in_data(:) = p_c%cells%decomp_info%glb_index(:)

    ALLOCATE(out_data(4,p_p%n_patch_cells))
    out_data(:,:) = 0

    ! communicate child index between processors:
    CALL reshuffle("send cell child indices to parent", dst_idx, in_data, p_p%cells%decomp_info%glb_index, &
      &            p_p%n_patch_cells_g, communicator, out_data)
    DEALLOCATE(in_data, dst_idx)      

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- each child cell sends its three edge indices to the parent
    !     cell. Edges which arrive twice are those of the "central"
    !     child.

    ALLOCATE(dst_idx(3*p_c%n_patch_cells), in_data(3*p_c%n_patch_cells))

    ! get the cell edge indices:
    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      ! loop only over inner domain
      IF (p_c%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE
      DO i=1,3
        jc_e = p_c%cells%edge_idx(jc_c,jb_c,i)
        jb_e = p_c%cells%edge_blk(jc_c,jb_c,i)
        iidx = iidx + 1
        dst_idx(iidx) = idx_1d(parent_idx_c(jc_c,jb_c), &
          &                    parent_blk_c(jc_c,jb_c))
        in_data(iidx) = p_c%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))
      END DO
    END DO

    ALLOCATE(out_data_e(9,p_p%n_patch_cells), & ! each refined cell has 9 sub-edges
      &      out_count(9,p_p%n_patch_cells))
    out_data_e(:,:) = 0
    out_count(:,:)  = 0

    ! communicate between processors:
    CALL reshuffle("send cell-edge indices to parent", dst_idx(1:iidx), in_data(1:iidx), &
      &            p_p%cells%decomp_info%glb_index, p_p%n_patch_cells_g, communicator,   &
      &            out_data_e, out_count)

    DO j = 1, p_p%n_patch_cells
      IF (SUM(out_count(:,j)) == 0)  CYCLE

      ! consistency check: there must be exactly 3 inner edges
      IF (COUNT(out_count(:,j) > 1) /= 3) THEN
        CALL finish(routine, "edge counting went wrong")
      END IF
    END DO

    DEALLOCATE(in_data, dst_idx)      

    ! --- PARENT PATCH -> CHILD PATCH
    ! --- now, send back the "counting" of the edges to the child patch

    ALLOCATE(dst_idx(12*p_p%n_patch_cells), in_data(12*p_p%n_patch_cells))
    dst_idx = 0
    in_data = 0

    iidx = 0
    DO j = 1, p_p%n_patch_cells
      IF (SUM(out_count(:,j)) == 0)  CYCLE

      DO i=1,SIZE(out_data_e,1)
        IF (out_count(i,j) > 1) THEN
          DO i1=1,4
            iidx = iidx + 1
            dst_idx(iidx) = out_data(i1,j)
            in_data(iidx) = out_data_e(i,j)
          END DO
        END IF
      END DO
    END DO

    ALLOCATE(out_data_e3(3,p_c%n_patch_cells))
    out_data_e3(:,:) = 0

    ! communicate between processors:
    CALL reshuffle("send inner child edges", dst_idx(1:iidx), in_data(1:iidx),         &
      &            p_c%cells%decomp_info%glb_index, p_c%n_patch_cells_g, communicator, &
      &            out_data_e3)

    DEALLOCATE(in_data, dst_idx,out_data_e, out_count)      

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- determine the global index of the "inner child cell" no. 3
    !     and send this to parent cell:

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))
    dst_idx = 0
    in_data = 0

    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      ! loop only over inner domain
      IF (p_c%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE

      lfound = .TRUE.
      EDGE_LOOP : DO i=1,3  ! loop over edges
        jc_e = p_c%cells%edge_idx(jc_c,jb_c,i)
        jb_e = p_c%cells%edge_blk(jc_c,jb_c,i)
        iglb = p_c%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))
        IF (.NOT. ANY(out_data_e3(:,idx_1d(jc_c,jb_c)) == iglb)) THEN
          lfound = .FALSE. ; EXIT EDGE_LOOP
        END IF
      END DO EDGE_LOOP

      IF (lfound) THEN
        iidx = iidx + 1
        in_data(iidx) = p_c%cells%decomp_info%glb_index(j)
        dst_idx(iidx) = idx_1d(parent_idx_c(jc_c,jb_c), &
          &                    parent_blk_c(jc_c,jb_c))
      END IF
    END DO

    ALLOCATE(out_data_c3(1,p_p%n_patch_cells))
    out_data_c3(:,:) = 0

    ! communicate between processors:
    CALL reshuffle("send child 3 to parent cell", dst_idx(1:iidx), in_data(1:iidx),    &
      &            p_p%cells%decomp_info%glb_index, p_p%n_patch_cells_g, communicator, &
      &            out_data_c3)

    DEALLOCATE(in_data, dst_idx,out_data_e3)      

    ! Now, we order the (global) child indices, s.t. the child 3 is in
    ! correct position.
    DO j = 1, p_p%n_patch_cells
      ! find position of child no. 3
      i1 = 0
      DO i2=1,4
        IF (out_data(i2,j) == out_data_c3(1,j))  i1=i2
      END DO
      IF (i1 /= 3) THEN
        ! we need to reorder; 
        tmp = out_data(i1,j) ; out_data(i1,j)=out_data(3,j) ; out_data(3,j)=tmp
      END IF
    END DO

    DEALLOCATE(out_data_c3)

    ! --- PARENT PATCH -> CHILD PATCH
    ! --- create CELL pc_idx array ------------------------------------
    !
    ! The "pc_idx" contains the relative ordering of the child cells
    ! in the parent cell.

    ALLOCATE(dst_idx(4*p_p%n_patch_cells), in_data(4*p_p%n_patch_cells))

    iidx = 0
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      IF (out_data(1,j) > 0) THEN
        DO i=1,4
          iidx = iidx + 1
          dst_idx(iidx) = out_data(i,j)
          in_data(iidx) = i
        END DO
      END IF
    END DO

    ALLOCATE(out_pc_data(1,p_c%n_patch_cells))
    out_pc_data(:,:) = 0

    ! communicate child index between processors:
    CALL reshuffle("send pc_idx to child", dst_idx(1:iidx), in_data(1:iidx), p_c%cells%decomp_info%glb_index, &
      &            p_c%n_patch_cells_g, communicator, out_pc_data)

    ! --- as a result, we now have
    !     out_pc_data(1,j) = p_c%cells%pc_idx(j)
    !     out_data(i,j) = cells%child_index(j,i)
    IF (set_pc_idx) THEN
      DO j = 1, p_c%n_patch_cells
        jc_c = idx_no(j)
        jb_c = blk_no(j)
        IF (out_pc_data(1,j) <= 0)  CYCLE
        p_c%cells%pc_idx(jc_c,jb_c) = out_pc_data(1,j)
      END DO
    END IF
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      IF (out_data(1,j) <= 0)  CYCLE
      p_p%cells%child_idx(jc_c,jb_c,:) = idx_no( out_data(:,j) )
      p_p%cells%child_blk(jc_c,jb_c,:) = blk_no( out_data(:,j) )
    END DO

    DEALLOCATE(out_pc_data,out_data, in_data, dst_idx)      

    ! -----------------------------------------------------------------
    ! --- create EDGE indices -----------------------------------------
    ! -----------------------------------------------------------------
    !
    ! The situation is a bit more complicated here: For a given edge,
    ! children 1 and 2 are the sub-triangle edges coinciding with the
    ! parent edge; children 3 and 4 are the edges of the inner
    ! sub-triangles having (roughly) the same orientation as the
    ! parent edge.

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- in a first step, we communicate child edges 3 and 4:
    ALLOCATE(dst_idx(3*p_c%n_patch_cells), in_data(3*p_c%n_patch_cells))

    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      IF (p_c%cells%pc_idx(jc_c,jb_c) == 3) THEN
        DO i=1,3
          jc_e = p_c%cells%edge_idx(jc_c,jb_c,i)
          jb_e = p_c%cells%edge_blk(jc_c,jb_c,i)
          iidx = iidx + 1
          dst_idx(iidx) = idx_1d(parent_idx_e(jc_e,jb_e), &
            &                    parent_blk_e(jc_e,jb_e))
          in_data(iidx) = p_c%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))
        END DO
      END IF
    END DO

    ALLOCATE(out_data34(2,p_p%n_patch_edges))
    out_data34(:,:) = 0

    ! communicate data between processors:
    CALL reshuffle("send child edges 3,4", dst_idx(1:iidx), in_data(1:iidx), p_p%edges%decomp_info%glb_index, &
      &            p_p%n_patch_edges_g, communicator, out_data34)
    DEALLOCATE(in_data, dst_idx)      

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- in a second step, we communicate *all four* edge children

    ALLOCATE(dst_idx(p_c%n_patch_edges), in_data(p_c%n_patch_edges))
    DO j = 1, p_c%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)
      dst_idx(j) = idx_1d(parent_idx_e(jc_e,jb_e), &
        &                 parent_blk_e(jc_e,jb_e))
      in_data(j) = p_c%edges%decomp_info%glb_index(j)
    END DO

    ALLOCATE(out_data(4,p_p%n_patch_edges))
    out_data(:,:) = 0

    ! communicate data between processors:
    CALL reshuffle("send edge indices to parent", dst_idx, in_data, p_p%edges%decomp_info%glb_index, &
      &            p_p%n_patch_edges_g, communicator, out_data)

    ! -----------------------------------------------------------------
    ! - copy results

    DO j = 1, p_p%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)

      ! we need to reorder the edge children 3,4 to the end
      ordered(:) = 0
      iidx = 0
      DO i=1,4 
        IF (.NOT. ANY(out_data34(:,j) == out_data(i,j))) THEN
          iidx = iidx + 1
          ordered(iidx) = out_data(i,j)
        END IF
      END DO
      ordered(3:4)  = out_data34(:,j)
      out_data(:,j) = ordered(:)
    END DO
    DEALLOCATE(dst_idx, in_data)

    ! --- PARENT PATCH -> CHILD PATCH
    ! --- create EDGE pc_idx array ------------------------------------
    !
    ! The "pc_idx" contains the relative ordering of the child edges (1...4).

    ALLOCATE(dst_idx(4*p_p%n_patch_edges), in_data(4*p_p%n_patch_edges))

    iidx = 0
    DO j = 1, p_p%n_patch_edges
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      DO i=1,4
        IF (out_data(i,j) > 0) THEN
          iidx = iidx + 1
          dst_idx(iidx) = out_data(i,j)
          in_data(iidx) = i
        END IF
      END DO
    END DO

    ALLOCATE(out_pc_data(1,p_c%n_patch_edges))
    out_pc_data(:,:) = 0

    ! communicate child index between processors:
    CALL reshuffle("send edge pc_idx", dst_idx(1:iidx), in_data(1:iidx), p_c%edges%decomp_info%glb_index, &
      &            p_c%n_patch_edges_g, communicator, out_pc_data)

    ! --- as a result, we now have
    !     out_pc_data(1,j) = p_c%edges%pc_idx(j)
    !       and 
    !     out_data(:,j) = p_p%edges%child_idx(j,:)
    IF (set_pc_idx) THEN
      DO j = 1, p_c%n_patch_edges
        jc_e = idx_no(j)
        jb_e = blk_no(j)
        p_c%edges%pc_idx(jc_e,jb_e) = out_pc_data(1,j)
      END DO
    END IF
    DO j = 1, p_p%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)
      IF (ALL(out_data(:,j) == 0))  CYCLE
      p_p%edges%child_idx(jc_e,jb_e,:) = idx_no( out_data(:,j) )
      p_p%edges%child_blk(jc_e,jb_e,:) = blk_no( out_data(:,j) )
    END DO

    DEALLOCATE(out_data, out_data34, in_data, dst_idx)      

    ! -----------------------------------------------------------------
    ! --- create CHILD ID's -------------------------------------------
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! --- calculate cells%child_id
    !
    ! To calculate the child_id for each child cell, we loop over
    ! the cells of the child patch and collect the global indices of
    ! their parent cells. Then we send the ID of the child domain to
    ! these parent cells. 

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))

    ! get the parent cell indices:
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      dst_idx(j) = idx_1d(parent_idx_c(jc_c,jb_c), &
        &                 parent_blk_c(jc_c,jb_c))
    END DO
    ! set the cell child_id for the (global) parent grid
    in_data(:) = p_c%id

    ALLOCATE(out_child_id(1,p_p%n_patch_cells))
    out_child_id = -1

    ! communicate ID data between processors:
    CALL reshuffle("send cell id to parent", dst_idx, in_data, p_p%cells%decomp_info%glb_index, &
      &            p_p%n_patch_cells_g, communicator, out_child_id)

    ! communication finished. now copy the result to the local arrays:
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      IF (out_child_id(1,j) < 0)  CYCLE

      ! check for overlap (i.e.: child id already exists)
      IF (p_p%cells%child_id(jc_c,jb_c) == 0) THEN
        p_p%cells%child_id(jc_c,jb_c) = out_child_id(1,j)
      ELSE
        IF (p_p%cells%child_id(jc_c,jb_c) /= out_child_id(1,j)) THEN
          CALL finish(routine, "Cell in DOM "//int2string(p_c%id,'(i0)')//" overlaps with another domain!")
        END IF
      END IF
    END DO
    DEALLOCATE(out_child_id, in_data, dst_idx)

    ! -----------------------------------------------------------------
    ! --- calculate edges%child_id
    !
    ! To calculate the child_id for each child edge, we loop over the
    ! edges of the child domain and collect the global indices of
    ! their parent edges. Then we send the ID of the child domain to
    ! these parent edges. 

    ALLOCATE(dst_idx(p_c%n_patch_edges), in_data(p_c%n_patch_edges))

    ! get the parent cell indices:
    DO j = 1, p_c%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)
      dst_idx(j) = idx_1d(parent_idx_e(jc_e,jb_e), &
        &                 parent_blk_e(jc_e,jb_e))
    END DO
    ! set the edge child_id for the (global) parent grid:
    in_data(:) = p_c%id

    ALLOCATE(out_child_id(1,p_p%n_patch_edges))
    out_child_id = -1

    ! communicate ID data between processors:
    CALL reshuffle("send edge ID to parent", dst_idx, in_data, p_p%edges%decomp_info%glb_index, &
      &            p_p%n_patch_edges_g, communicator, out_child_id)

    ! communication finished. now copy the result to the local arrays:
    DO j = 1, p_p%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)

      IF (out_child_id(1,j) < 0)  CYCLE

      ! check for overlap (i.e.: child id already exists)
      IF (p_p%edges%child_id(jc_e,jb_e) == 0) THEN
        p_p%edges%child_id(jc_e,jb_e) = out_child_id(1,j)
      ELSE
        IF (p_p%edges%child_id(jc_e,jb_e) /= out_child_id(1,j)) THEN
          CALL finish(routine, "Edge in DOM "//int2string(p_c%id,'(i0)')//" overlaps with another domain!")
        END IF
      END IF
    END DO
    DEALLOCATE(out_child_id, in_data, dst_idx)

    IF (is_local_parent) THEN
      DEALLOCATE(lparent_glb_idx_c, lparent_glb_blk_c, lparent_glb_idx_e, lparent_glb_blk_e)
    END IF

    IF (p_pe_work == 0) THEN
      WRITE (0,*) "set_child_indices: ", TRIM(description), " - Done."
    END IF
  END SUBROUTINE set_child_indices

 
END MODULE mo_build_decomposition
