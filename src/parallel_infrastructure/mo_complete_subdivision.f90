!>
!!               This module provides all routines for dividing patches.
!!
!!               This module provides all routines for dividing patches
!! (including interpolation state) and setting up communication.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_complete_subdivision
  ! If METIS is installed, uncomment the following line
  ! (or better adjust configure to recognize that)
  !
  !-------------------------------------------------------------------------
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: success, min_rlcell, max_rlcell,  &
    & min_rledge, max_rledge, min_rlvert, max_rlvert,                &
    & min_rlcell_int, min_rledge_int, min_rlvert_int, max_phys_dom
  USE mo_exception,          ONLY: finish, message, message_text,    &
    &                              get_filename_noext

  USE mo_model_domain,       ONLY: t_patch, p_patch,      &
    &                              p_patch_local_parent,  &
    &                              p_phys_patch
  USE mo_mpi,                ONLY: p_send, p_recv, p_max, p_min, proc_split
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_UNDEFINED, MPI_COMM_NULL
#endif
  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_test, &
    & my_process_is_mpi_seq, process_mpi_all_test_id, process_mpi_all_workroot_id, &
    & my_process_is_mpi_workroot, p_pe_work, p_n_work,                  &
    & get_my_mpi_all_id, my_process_is_mpi_parallel

  USE mo_parallel_config,    ONLY:  nproma, p_test_run    
  USE mo_communication,      ONLY: setup_comm_pattern, blk_no, idx_no, idx_1d
  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
    & grf_bdyintp_end_c, grf_bdyintp_end_e, grf_fbk_start_c, grf_fbk_start_e, &
    & grf_bdywidth_c, grf_bdywidth_e, grf_nudgintp_start_c, grf_nudgintp_start_e
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, n_phys_dom

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: copy_processor_splitting
  PUBLIC :: set_patch_communicators
  PUBLIC :: finalize_decomposition
  PUBLIC :: setup_phys_patches
  PUBLIC :: set_comm_pat_gather
  PUBLIC :: complete_parallel_setup
  PUBLIC :: complete_parallel_setup_oce
  PUBLIC :: finalize_decomposition_oce
  !-------------------------------------------------------------------------
  ! Definition of local parent patches
  ! For any given patch p_patch(jg) and jgp = p_patch(jg)%parent_id,
  ! p_patch_local_parent(jg) has the same resolution as p_patch(jgp)
  ! but it covers only the area of p_patch(jgp) which is covered by its child p_patch(jg)
  ! and it is divided in the same manner as p_patch(jg).
  ! Please note that p_patch_local_parent(1) is undefined if n_dom_start = 1

  ! Please note: The definitions of the local parents are now at the same locations
  ! as the definitions of the respective patch or state
  !-------------------------------------------------------------------------


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!  In case of a test run: Copies processor splitting to test PE
  SUBROUTINE copy_processor_splitting(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER :: ibuf(n_dom_start:n_dom,2)

    IF(.NOT. p_test_run) RETURN ! Nothing to do

    IF(my_process_is_mpi_workroot()) THEN
      CALL p_send(proc_split, process_mpi_all_test_id, 1)
      ibuf(:,1) = p_patch(:)%n_proc
      ibuf(:,2) = p_patch(:)%proc0
      CALL p_send(ibuf, process_mpi_all_test_id, 2)
    ENDIF

    IF(my_process_is_mpi_test()) THEN
      CALL p_recv(proc_split, process_mpi_all_workroot_id, 1)
      CALL p_recv(ibuf, process_mpi_all_workroot_id, 2)
      p_patch(:)%n_proc = ibuf(:,1)
      p_patch(:)%proc0  = ibuf(:,2)
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

  SUBROUTINE complete_parallel_setup

    INTEGER :: jg, jgp

    DO jg = n_dom_start, n_dom

      jgp = p_patch(jg)%parent_id

      ! Set communication patterns for boundary exchange
      CALL set_comm_pat_bound_exch(p_patch(jg))

      ! Set communication patterns for gathering on proc 0
      CALL set_comm_pat_gather(p_patch(jg))

      CALL set_owner_mask(p_patch(jg))
     
     ! Fill the owner_local value
      ! this is done in the set_owner_mask
      ! CALL fill_owner_local(p_patch(jg))

      IF(jg == n_dom_start) THEN

        ! parent_idx/blk is set to 0 since it just doesn't exist,
        ! child_idx/blk is set to 0 since it makes sense only on the local parent
        p_patch(jg)%cells%parent_idx = 0
        p_patch(jg)%cells%parent_blk = 0
        p_patch(jg)%cells%child_idx  = 0
        p_patch(jg)%cells%child_blk  = 0
        p_patch(jg)%edges%parent_idx = 0
        p_patch(jg)%edges%parent_blk = 0
        p_patch(jg)%edges%child_idx  = 0
        p_patch(jg)%edges%child_blk  = 0

      ELSE

        CALL setup_comm_cpy_interpolation(p_patch(jg), p_patch(jgp))
        CALL setup_comm_grf_interpolation(p_patch(jg), p_patch(jgp))
        CALL setup_comm_ubc_interpolation(p_patch(jg), p_patch(jgp))

        CALL set_comm_pat_bound_exch(p_patch_local_parent(jg))
        CALL set_comm_pat_gather(p_patch_local_parent(jg))

        CALL set_parent_child_relations(p_patch_local_parent(jg), p_patch(jg))

        CALL set_glb_loc_comm(p_patch(jgp), p_patch_local_parent(jg), &
          &                   p_patch(jg)%parent_child_index)

        CALL set_owner_mask(p_patch_local_parent(jg))
      ENDIF

    ENDDO

    CALL set_patch_communicators(p_patch)

  END SUBROUTINE complete_parallel_setup

  !-----------------------------------------------------------------------------

  SUBROUTINE finalize_decomposition

    implicit none
    integer jg

    ! Remap indices in patches and local parents

    DO jg = n_dom_start, n_dom

      CALL remap_patch_indices(p_patch(jg))

      IF(jg>n_dom_start) THEN
        CALL remap_patch_indices(p_patch_local_parent(jg))
      ENDIF

    ENDDO
                 
  END SUBROUTINE finalize_decomposition

  !-----------------------------------------------------------------------------
  !>
  ! Fills the in_patch%cells%owner_local using the in_patch%cells%owner_g
  ! Note: At the moment it uses the p_work_pe number which is not the same
  ! as the my_mpi_all_id. It requires 
!   SUBROUTINE fill_owner_local(in_patch)
! 
!     TYPE(t_patch), INTENT(inout) :: in_patch
! 
!     INTEGER :: local_cell_idx, global_cell_idx
!     INTEGER :: i, jb, jl, jb_e, jl_e, jb_v, jl_v, jv, je
!     INTEGER :: owner_id
! 
!     in_patch%edges%owner_local(:) = -1
!     in_patch%verts%owner_local(:) = -1
! 
!     DO local_cell_idx = 1, in_patch%n_patch_cells
!       global_cell_idx = in_patch%cells%glb_index(local_cell_idx)
!       owner_id = in_patch%cells%owner_g(global_cell_idx)
!       in_patch%cells%owner_local(local_cell_idx) = in_patch%cells%owner_g(global_cell_idx)
!       IF (owner_id < 0) CYCLE
! 
!       ! go around the cell edges mark the owner
!       jb = blk_no(local_cell_idx) ! block index
!       jl = idx_no(local_cell_idx) ! line index
!       DO i = 1,in_patch%cells%num_edges(jl,jb)
!         jl_e = in_patch%cells%edge_idx(jl,jb,i)
!         jb_e = in_patch%cells%edge_blk(jl,jb,i)
!         je = idx_1d(jl_e, jb_e)
!         jl_v = in_patch%cells%vertex_idx(jl,jb,i)
!         jb_v = in_patch%cells%vertex_blk(jl,jb,i)
!         jv = idx_1d(jl_v, jb_v)
!         IF (owner_id == p_pe_work) THEN
!           ! Assume we own the edges and vertices of the owned cells
!           ! No process can claim actual ownershipe, as these can be calculated
!           ! only using the halos. No reason to communicate them
!           in_patch%edges%owner_local(je) = p_pe_work
!           in_patch%verts%owner_local(jv) = p_pe_work
!         ELSEIF (in_patch%edges%owner_local(je) < 0) THEN
!           in_patch%edges%owner_local(je) = owner_id
!           in_patch%verts%owner_local(jv) = owner_id
!         ENDIF
! 
!       ENDDO
! 
!     ENDDO
! 
! 
!   END SUBROUTINE fill_owner_local
 
  !-----------------------------------------------------------------------------
  !>
  !! Sets the owner mask
  SUBROUTINE set_owner_mask(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch

    INTEGER :: j, jb, jl, jg

    p_patch%cells%owner_mask = .false.
    p_patch%edges%owner_mask = .false.
    p_patch%verts%owner_mask = .false.
      
    p_patch%cells%owner_local(:) = -1
    p_patch%edges%owner_local(:) = -1
    p_patch%verts%owner_local(:) = -1

    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      jg = p_patch%cells%glb_index(j) ! index in global patch

      p_patch%cells%owner_mask(jl,jb) = p_patch%cells%owner_g(jg)==p_pe_work
      ! fill local owner
      p_patch%cells%owner_local(j) = p_patch%cells%owner_g(jg)

    ENDDO

    DO j = 1, p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      jg = p_patch%edges%glb_index(j) ! index in global patch

      p_patch%edges%owner_mask(jl,jb) = p_patch%edges%owner_g(jg)==p_pe_work
      ! fill local owner
      p_patch%edges%owner_local(j) = p_patch%edges%owner_g(jg)
    
    ENDDO

    DO j = 1, p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      jg = p_patch%verts%glb_index(j) ! index in global patch

      p_patch%verts%owner_mask(jl,jb) = p_patch%verts%owner_g(jg)==p_pe_work
      ! fill local owner
      p_patch%verts%owner_local(j) = p_patch%verts%owner_g(jg)

    ENDDO

  END SUBROUTINE set_owner_mask

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
        CALL remap_index(p_patch%cells%loc_index, &
          & p_patch%cells%neighbor_idx(jl,jb,i),      &
          & p_patch%cells%neighbor_blk(jl,jb,i))

        ! edge_idx and vertex_idx should not need a remap !!!
        ! This is only left here if there should change something in the decomposition
!CDIR IEXPAND
        CALL remap_index(p_patch%edges%loc_index, &
          & p_patch%cells%edge_idx(jl,jb,i),          &
          & p_patch%cells%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%verts%loc_index, &
          & p_patch%cells%vertex_idx(jl,jb,i),        &
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

      DO i=1,2
!CDIR IEXPAND
        CALL remap_index(p_patch%cells%loc_index, &
          & p_patch%edges%cell_idx(jl,jb,i),          &
          & p_patch%edges%cell_blk(jl,jb,i))
      ENDDO

      DO i=1,4
!CDIR IEXPAND
        CALL remap_index(p_patch%verts%loc_index, &
          & p_patch%edges%vertex_idx(jl,jb,i),        &
          & p_patch%edges%vertex_blk(jl,jb,i))
      ENDDO

      DO i=1,4
!CDIR IEXPAND
        CALL remap_index(p_patch%edges%loc_index, &
          & p_patch%edges%quad_idx(jl,jb,i),          &
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
        CALL remap_index(p_patch%verts%loc_index, &
          & p_patch%verts%neighbor_idx(jl,jb,i),      &
          & p_patch%verts%neighbor_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%edges%loc_index, &
          & p_patch%verts%edge_idx(jl,jb,i),          &
          & p_patch%verts%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%cells%loc_index, &
          & p_patch%verts%cell_idx(jl,jb,i),          &
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

      IF(p_pp%cells%decomp_domain(jl,jb)>0) THEN
        p_pp%cells%child_idx(jl,jb,:) = 0
        p_pp%cells%child_blk(jl,jb,:) = 0
        CYCLE
      ENDIF

      DO i= 1, 4
        jc_g = idx_1d(p_pp%cells%child_idx(jl,jb,i),p_pp%cells%child_blk(jl,jb,i))
        IF(jc_g<1 .OR. jc_g>p_pc%n_patch_cells_g) &
          & CALL finish('set_parent_child_relations','Invalid cell child index in global parent')
        jc = p_pc%cells%loc_index(jc_g)
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

      IF(p_pp%edges%decomp_domain(jl,jb)>1) THEN
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
        IF(jc_g<1 .OR. jc_g>p_pc%n_patch_edges_g) &
          & CALL finish('set_parent_child_relations','Inv. edge child index in global parent')
        jc = p_pc%edges%loc_index(ABS(jc_g))
        IF(jc <= 0) &
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

      jp = p_pp%cells%loc_index(jp_g)
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

      jp = p_pp%edges%loc_index(jp_g)
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

    p_pc%cells%child_idx  = 0
    p_pc%cells%child_blk  = 0
    p_pp%cells%parent_idx = 0
    p_pp%cells%parent_blk = 0

    p_pc%edges%child_idx  = 0
    p_pc%edges%child_blk  = 0
    p_pp%edges%parent_idx = 0
    p_pp%edges%parent_blk = 0

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

    owner_c(:) = p%cells%owner_g(p%cells%glb_index(:))
    WHERE(owner_c(:) == p_pe_work) owner_c(:) = -1

    owner_e(:) = p%edges%owner_g(p%edges%glb_index(:))
    WHERE(owner_e(:) == p_pe_work) owner_e(:) = -1

    owner_v(:) = p%verts%owner_g(p%verts%glb_index(:))
    WHERE(owner_v(:) == p_pe_work) owner_v(:) = -1

    ! Set communication patterns for boundary exchange
    CALL setup_comm_pattern(p%n_patch_cells, owner_c, p%cells%glb_index, &
      & p%cells%loc_index, p%comm_pat_c)

    CALL setup_comm_pattern(p%n_patch_edges, owner_e, p%edges%glb_index, &
      & p%edges%loc_index, p%comm_pat_e)

    CALL setup_comm_pattern(p%n_patch_verts, owner_v, p%verts%glb_index, &
      & p%verts%loc_index, p%comm_pat_v)

    DEALLOCATE(owner_c, owner_e, owner_v)

    ! Set reduced communication pattern containing only level-1 halo cells (immediate neighbors)
    jc = idx_1d(p%cells%end_idx(min_rlcell_int-1,1), &
                p%cells%end_blk(min_rlcell_int-1,1))
    ALLOCATE(owner_c(jc))
    owner_c(1:jc) = p%cells%owner_g(p%cells%glb_index(1:jc))
    WHERE(owner_c(:) == p_pe_work) owner_c(:) = -1

    CALL setup_comm_pattern(jc, owner_c, p%cells%glb_index, &
      & p%cells%loc_index, p%comm_pat_c1)

    DEALLOCATE(owner_c)

  END SUBROUTINE set_comm_pat_bound_exch

  !-------------------------------------------------------------------------------------------------
  !> Sets the gather communication patterns of a patch

  SUBROUTINE set_comm_pat_gather(p)

    TYPE(t_patch), INTENT(INOUT):: p

    INTEGER, ALLOCATABLE :: tmp(:)
    INTEGER :: j

    ! For gathering the global fields on p_pe_work==0
    ALLOCATE(tmp(MAX(p%n_patch_cells_g, p%n_patch_edges_g, p%n_patch_verts_g)))

    DO j = 1, SIZE(tmp)
      tmp(j) = j ! Global/local index in global array, i.e. identity!
    ENDDO

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(p%n_patch_cells_g, p%cells%owner_g, tmp, &
        & p%cells%loc_index, p%comm_pat_gather_c)
    ELSE
      ! We don't want to receive any data, i.e. the number of cells is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, p%cells%owner_g, tmp, &
        & p%cells%loc_index, p%comm_pat_gather_c)
    ENDIF

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(p%n_patch_edges_g, p%edges%owner_g, tmp, &
        & p%edges%loc_index, p%comm_pat_gather_e)
    ELSE
      ! We don't want to receive any data, i.e. the number of edges is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, p%edges%owner_g, tmp, &
        & p%edges%loc_index, p%comm_pat_gather_e)
    ENDIF

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(p%n_patch_verts_g, p%verts%owner_g, tmp, &
        & p%verts%loc_index, p%comm_pat_gather_v)
    ELSE
      ! We don't want to receive any data, i.e. the number of edges is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, p%verts%owner_g, tmp, &
        & p%verts%loc_index, p%comm_pat_gather_v)
    ENDIF

    DEALLOCATE(tmp)

  END SUBROUTINE set_comm_pat_gather

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets up communication patterns between global and local parent patches

  SUBROUTINE set_glb_loc_comm(p_pglb, p_ploc, i_chidx)
    TYPE(t_patch), INTENT(IN)    :: p_pglb !> global parent
    TYPE(t_patch), INTENT(INOUT) :: p_ploc !> local parent
    INTEGER, INTENT(IN) :: i_chidx

    INTEGER, ALLOCATABLE :: owner(:)
    INTEGER :: j, js, je, icid, jb, jl

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

    ALLOCATE(owner(p_ploc%n_patch_cells))

    js = idx_1d(p_ploc%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
      &         p_ploc%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    je = idx_1d(p_ploc%cells%end_idx(min_rlcell_int,i_chidx), &
      &         p_ploc%cells%end_blk(min_rlcell_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_pglb%cells%owner_g(p_ploc%cells%glb_index(j))
    ENDDO

    CALL setup_comm_pattern(p_ploc%n_patch_cells, owner, p_ploc%cells%glb_index,  &
      & p_pglb%cells%loc_index, p_ploc%comm_pat_glb_to_loc_c)

    DEALLOCATE(owner)

    ! ... edges

    ALLOCATE(owner(p_ploc%n_patch_edges))

    js = idx_1d(p_ploc%edges%start_idx(grf_bdyintp_start_e,i_chidx), &
      &         p_ploc%edges%start_blk(grf_bdyintp_start_e,i_chidx))
    je = idx_1d(p_ploc%edges%end_idx(min_rledge_int,i_chidx), &
      &         p_ploc%edges%end_blk(min_rledge_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_pglb%edges%owner_g(p_ploc%edges%glb_index(j))
    ENDDO

    CALL setup_comm_pattern(p_ploc%n_patch_edges, owner, p_ploc%edges%glb_index,  &
      & p_pglb%edges%loc_index, p_ploc%comm_pat_glb_to_loc_e)

    DEALLOCATE(owner)

    !-----------------------------------------------------------------------------------------------

    ! Communication local -> global
    ! Here it might get necessary to have different patterns for different start levels
    ! of the copy to the global parent since we may not overwrite arbitrary values there.
    ! Currently only one pattern for feedback is needed (starting at grf_fbk_start_c/e).

    ! ... cells

    ALLOCATE(owner(p_pglb%n_patch_cells))

    js = idx_1d(p_pglb%cells%start_idx(grf_fbk_start_c,i_chidx), &
      &         p_pglb%cells%start_blk(grf_fbk_start_c,i_chidx))
    je = idx_1d(p_pglb%cells%end_idx(min_rlcell_int,i_chidx), &
      &         p_pglb%cells%end_blk(min_rlcell_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_ploc%cells%owner_g(p_pglb%cells%glb_index(j))
    ENDDO

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      js = idx_1d(p_pglb%cells%start_idx(min_rlcell_int-1,1), &
        &         p_pglb%cells%start_blk(min_rlcell_int-1,1))
      je = idx_1d(p_pglb%cells%end_idx(min_rlcell,MAX(1,p_pglb%n_childdom)), &
        &         p_pglb%cells%end_blk(min_rlcell,MAX(1,p_pglb%n_childdom)))

      DO j = js, je

        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index
        IF (p_pglb%cells%child_id(jl,jb)   == icid            .AND. &
            p_pglb%cells%refin_ctrl(jl,jb) <= grf_fbk_start_c .AND. &
            p_pglb%cells%refin_ctrl(jl,jb) >= min_rlcell_int )   THEN

          owner(j) = p_ploc%cells%owner_g(p_pglb%cells%glb_index(j))
        ENDIF
      ENDDO

    ENDIF

    CALL setup_comm_pattern(p_pglb%n_patch_cells, owner, p_pglb%cells%glb_index, &
      & p_ploc%cells%loc_index, p_ploc%comm_pat_loc_to_glb_c_fbk)

    DEALLOCATE(owner)

    ! ... edges

    ALLOCATE(owner(p_pglb%n_patch_edges))

    js = idx_1d(p_pglb%edges%start_idx(grf_fbk_start_e,i_chidx), &
      &         p_pglb%edges%start_blk(grf_fbk_start_e,i_chidx))
    je = idx_1d(p_pglb%edges%end_idx(min_rledge_int,i_chidx), &
      &         p_pglb%edges%end_blk(min_rledge_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_ploc%edges%owner_g(p_pglb%edges%glb_index(j))
    ENDDO

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      js = idx_1d(p_pglb%edges%start_idx(min_rledge_int-1,1), &
        &         p_pglb%edges%start_blk(min_rledge_int-1,1))
      je = idx_1d(p_pglb%edges%end_idx(min_rledge,MAX(1,p_pglb%n_childdom)), &
        &         p_pglb%edges%end_blk(min_rledge,MAX(1,p_pglb%n_childdom)))

      DO j = js, je

        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index
        IF (p_pglb%edges%child_id(jl,jb)   == icid            .AND. &
            p_pglb%edges%refin_ctrl(jl,jb) <= grf_fbk_start_e .AND. &
            p_pglb%edges%refin_ctrl(jl,jb) >= min_rledge_int )   THEN

          owner(j) = p_ploc%edges%owner_g(p_pglb%edges%glb_index(j))
        ENDIF
      ENDDO

    ENDIF

    CALL setup_comm_pattern(p_pglb%n_patch_edges, owner, p_pglb%edges%glb_index, &
      & p_ploc%edges%loc_index, p_ploc%comm_pat_loc_to_glb_e_fbk)

    DEALLOCATE(owner)

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

    ! This routine must not be called in a single CPU run
    IF(my_process_is_mpi_seq()) &
      & CALL finish('setup_comm_cpy_interpolation','must not be called in a single CPU run')

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Get start and end index of the GLOBAL parent cells as used in the interpolation

    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(grf_bdyintp_end_c,i_chidx), &
                       p_parent_patch%cells%end_blk(grf_bdyintp_end_c,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    glb_index(:) = -1
    owner(:)     = -1

    DO j = 1, p_patch%n_patch_cells
      jc = idx_no(j)
      jb = blk_no(j)
      jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
      IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
      glb_index(j) = jp
      owner(j) = p_parent_patch%cells%owner_g(jp)
    ENDDO

    ! Set up communication pattern

    CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
      & p_parent_patch%cells%loc_index, &
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

    INTEGER :: j, n, jc, js, jl, je, jb, jp, p_index_s, p_index_e, i_chidx
    INTEGER :: num_send, num_recv, np, iss, ise, irs, ire
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(my_process_is_mpi_seq()) &
      & CALL finish('setup_comm_grf_interpolation','must not be called in a single CPU run')

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Start and end index of the GLOBAL parent cells as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(grf_bdyintp_end_c,i_chidx), &
                       p_parent_patch%cells%end_blk(grf_bdyintp_end_c,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_cells
        jc = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%cells%pc_idx(jc,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%cells%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
        & p_parent_patch%cells%loc_index, &
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

    ! Start and end index of the GLOBAL parent edges as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%edges%start_idx(grf_bdyintp_start_e,i_chidx), &
                       p_parent_patch%edges%start_blk(grf_bdyintp_start_e,i_chidx))
    p_index_e = idx_1d(p_parent_patch%edges%end_idx(grf_bdyintp_end_e,i_chidx), &
                       p_parent_patch%edges%end_blk(grf_bdyintp_end_e,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%edges%glb_index(p_index_s)
      p_index_e = p_parent_patch%edges%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(p_patch%n_patch_edges))
    ALLOCATE(owner(p_patch%n_patch_edges))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_edges
        je = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%edges%parent_idx(je,jb),p_patch%edges%parent_blk(je,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%edges%pc_idx(je,jb) /= n .OR. p_patch%edges%refin_ctrl(je,jb) > grf_bdywidth_e) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%edges%owner_g(jp)
      ENDDO


      ! include halo edges for nest boundary points in order to avoid extra synchronization
      js = idx_1d(p_patch%edges%start_idx(min_rledge_int-1,1), &
        &         p_patch%edges%start_blk(min_rledge_int-1,1))
      je = idx_1d(p_patch%edges%end_idx(min_rledge_int-2,1), &
        &         p_patch%edges%end_blk(min_rledge_int-2,1))

! GZ: this loop is not vectorized properly, probably because the nested IF clauses are merged
!     in an incorrect way. As the runtime cost of this loop is negligible, we just turn off vectorization
!CDIR NOVECTOR
      DO j = js, je
        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index
        IF (p_patch%edges%refin_ctrl(jl,jb) > 0 .AND. p_patch%edges%refin_ctrl(jl,jb) <= grf_bdywidth_e ) THEN
          jp = idx_1d(p_patch%edges%parent_idx(jl,jb),p_patch%edges%parent_blk(jl,jb))
          IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
          IF(p_patch%edges%pc_idx(jl,jb) /= n) CYCLE
          glb_index(j) = jp
          owner(j) = p_parent_patch%edges%owner_g(jp)
        ENDIF
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_edges, owner, glb_index,  &
        & p_parent_patch%edges%loc_index, &
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

    INTEGER :: j, n, jc, je, jb, jp, p_index_s, p_index_e, i_chidx
    INTEGER :: num_send, num_recv, np, iss, ise, irs, ire
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(my_process_is_mpi_seq()) &
      & CALL finish('setup_comm_ubc_interpolation','must not be called in a single CPU run')

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Start and end index of the GLOBAL parent cells as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_nudgintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_nudgintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(min_rlcell_int,i_chidx), &
                       p_parent_patch%cells%end_blk(min_rlcell_int,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_cells
        jc = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%cells%pc_idx(jc,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%cells%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
        & p_parent_patch%cells%loc_index, &
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

    ! Start and end index of the GLOBAL parent edges as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%edges%start_idx(grf_nudgintp_start_e,i_chidx), &
                       p_parent_patch%edges%start_blk(grf_nudgintp_start_e,i_chidx))
    p_index_e = idx_1d(p_parent_patch%edges%end_idx(min_rledge_int,i_chidx), &
                       p_parent_patch%edges%end_blk(min_rledge_int,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%edges%glb_index(p_index_s)
      p_index_e = p_parent_patch%edges%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(p_patch%n_patch_edges))
    ALLOCATE(owner(p_patch%n_patch_edges))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_edges
        je = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%edges%parent_idx(je,jb),p_patch%edges%parent_blk(je,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%edges%pc_idx(je,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%edges%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_edges, owner, glb_index,  &
        & p_parent_patch%edges%loc_index, &
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

  SUBROUTINE setup_phys_patches

    INTEGER :: jp, jg, n, i, j, jb, jl
    INTEGER, ALLOCATABLE :: glb_phys_id_c(:), glb_phys_id_e(:), glb_phys_id_v(:)
    INTEGER, ALLOCATABLE :: owner(:), glbidx(:)
    CHARACTER(LEN=*), PARAMETER :: routine = 'setup_phys_patches'

    p_phys_patch(:)%logical_id = -1

    DO jg = 1, n_dom

      ! Get global arrays for phys_id

      ! Allocate and set to 0
      ALLOCATE(glb_phys_id_c(p_patch(jg)%n_patch_cells_g))
      ALLOCATE(glb_phys_id_e(p_patch(jg)%n_patch_edges_g))
      ALLOCATE(glb_phys_id_v(p_patch(jg)%n_patch_verts_g))
      glb_phys_id_c(:) = 0
      glb_phys_id_e(:) = 0
      glb_phys_id_v(:) = 0

      ! Fill with own values of phys_id

      DO j = 1, p_patch(jg)%n_patch_cells
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF(.NOT.p_patch(jg)%cells%owner_mask(jl,jb)) CYCLE
        glb_phys_id_c(p_patch(jg)%cells%glb_index(j)) = p_patch(jg)%cells%phys_id(jl,jb)
      ENDDO

      DO j = 1, p_patch(jg)%n_patch_edges
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF(.NOT.p_patch(jg)%edges%owner_mask(jl,jb)) CYCLE
        glb_phys_id_e(p_patch(jg)%edges%glb_index(j)) = p_patch(jg)%edges%phys_id(jl,jb)
      ENDDO

      DO j = 1, p_patch(jg)%n_patch_verts
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF(.NOT.p_patch(jg)%verts%owner_mask(jl,jb)) CYCLE
        glb_phys_id_v(p_patch(jg)%verts%glb_index(j)) = p_patch(jg)%verts%phys_id(jl,jb)
      ENDDO

      ! Get global arrays by obtaining the global maximum
      ! Since p_max works on real valued arrays only, we have to convert to real and back

      glb_phys_id_c = INT( p_max(REAL(glb_phys_id_c,wp), comm=p_comm_work) )
      glb_phys_id_e = INT( p_max(REAL(glb_phys_id_e,wp), comm=p_comm_work) )
      glb_phys_id_v = INT( p_max(REAL(glb_phys_id_v,wp), comm=p_comm_work) )

      ! Get the physical patches contained within current patch

      ! Set logical id of p_phys_patch from cells

      DO j = 1, p_patch(jg)%n_patch_cells_g
        jp = glb_phys_id_c(j)
        IF(jp<1 .OR. jp>max_phys_dom) THEN
          WRITE(message_text,'(a,i4,a,i12)') &
            & 'Patch ',jg,' contains illegal value for cells phys_id: ',jp
          CALL finish  (routine, TRIM(message_text))
        ENDIF
        IF(p_phys_patch(jp)%logical_id < 0) THEN
          p_phys_patch(jp)%logical_id = jg
        ELSE
          ! Check if no other patch uses the same phys_id
          IF(p_phys_patch(jp)%logical_id /= jg) THEN
            WRITE(message_text,'(a,i4,a,i12,a,i12)') &
             & 'Patch ',jg,' contains cells phys_id: ',jp, &
             & ', already used by patch: ',p_phys_patch(jp)%logical_id
            CALL finish  (routine, TRIM(message_text))
          ENDIF
        ENDIF
      ENDDO

      ! Make a check for egdes and verts

      DO j = 1, p_patch(jg)%n_patch_edges_g
        jp = glb_phys_id_e(j)
        IF(jp<1 .OR. jp>max_phys_dom) THEN
          WRITE(message_text,'(a,i4,a,i12)') &
            & 'Patch ',jg,' contains illegal value for edges phys_id: ',jp
          CALL finish  (routine, TRIM(message_text))
        ENDIF
        ! Check if no other patch uses the same phys_id
        IF(p_phys_patch(jp)%logical_id /= jg) THEN
          WRITE(message_text,'(a,i4,a,i12,a,i12)') &
           & 'Patch ',jg,' contains edges phys_id: ',jp, &
           & ', already used by patch: ',p_phys_patch(jp)%logical_id
          CALL finish  (routine, TRIM(message_text))
        ENDIF
      ENDDO

      DO j = 1, p_patch(jg)%n_patch_verts_g
        jp = glb_phys_id_v(j)
        IF(jp<1 .OR. jp>max_phys_dom) THEN
          WRITE(message_text,'(a,i4,a,i12)') &
            & 'Patch ',jg,' contains illegal value for verts phys_id: ',jp
          CALL finish  (routine, TRIM(message_text))
        ENDIF
        ! Check if no other patch uses the same phys_id
        IF(p_phys_patch(jp)%logical_id /= jg) THEN
          WRITE(message_text,'(a,i4,a,i12,a,i12)') &
           & 'Patch ',jg,' contains verts phys_id: ',jp, &
           & ', already used by patch: ',p_phys_patch(jp)%logical_id
          CALL finish  (routine, TRIM(message_text))
        ENDIF
      ENDDO

      ! Set up communication patterns in physical patches

      n = MAX(p_patch(jg)%n_patch_cells_g,p_patch(jg)%n_patch_edges_g,p_patch(jg)%n_patch_verts_g)
      ALLOCATE(owner (n))
      ALLOCATE(glbidx(n))

      DO jp = 1, max_phys_dom ! Loop over physical patches

        IF(p_phys_patch(jp)%logical_id /= jg) CYCLE ! do only for physical patches belonging to jg

        ! cells

        n = 0
        DO i = 1, p_patch(jg)%n_patch_cells_g
          IF(glb_phys_id_c(i) == jp) THEN
            n = n+1
            owner (n) = p_patch(jg)%cells%owner_g(i)
            glbidx(n) = i
          ENDIF
        ENDDO

        p_phys_patch(jp)%n_patch_cells = n

        IF(.NOT.my_process_is_mpi_seq()) THEN
          IF(p_pe_work == 0) THEN
            CALL setup_comm_pattern(n, owner, glbidx, &
              & p_patch(jg)%cells%loc_index, p_phys_patch(jp)%comm_pat_gather_c)
          ELSE
            ! We don't want to receive any data, i.e. the number of cells is 0
            ! and owner/global index are dummies!
            CALL setup_comm_pattern(0, owner, glbidx, &
              & p_patch(jg)%cells%loc_index, p_phys_patch(jp)%comm_pat_gather_c)
          ENDIF
        ENDIF

        ! edges

        n = 0
        DO i = 1, p_patch(jg)%n_patch_edges_g
          IF(glb_phys_id_e(i) == jp) THEN
            n = n+1
            owner (n) = p_patch(jg)%edges%owner_g(i)
            glbidx(n) = i
          ENDIF
        ENDDO

        p_phys_patch(jp)%n_patch_edges = n

        IF(.NOT.my_process_is_mpi_seq()) THEN
          IF(p_pe_work == 0) THEN
            CALL setup_comm_pattern(n, owner, glbidx, &
              & p_patch(jg)%edges%loc_index, p_phys_patch(jp)%comm_pat_gather_e)
          ELSE
            ! We don't want to receive any data, i.e. the number of edges is 0
            ! and owner/global index are dummies!
            CALL setup_comm_pattern(0, owner, glbidx, &
              & p_patch(jg)%edges%loc_index, p_phys_patch(jp)%comm_pat_gather_e)
          ENDIF
        ENDIF

        ! verts

        n = 0
        DO i = 1, p_patch(jg)%n_patch_verts_g
          IF(glb_phys_id_v(i) == jp) THEN
            n = n+1
            owner (n) = p_patch(jg)%verts%owner_g(i)
            glbidx(n) = i
          ENDIF
        ENDDO

        p_phys_patch(jp)%n_patch_verts = n

        IF(.NOT.my_process_is_mpi_seq()) THEN
          IF(p_pe_work == 0) THEN
            CALL setup_comm_pattern(n, owner, glbidx, &
              & p_patch(jg)%verts%loc_index, p_phys_patch(jp)%comm_pat_gather_v)
          ELSE
            ! We don't want to receive any data, i.e. the number of verts is 0
            ! and owner/global index are dummies!
            CALL setup_comm_pattern(0, owner, glbidx, &
              & p_patch(jg)%verts%loc_index, p_phys_patch(jp)%comm_pat_gather_v)
          ENDIF
        ENDIF

      ENDDO

      DEALLOCATE(owner)
      DEALLOCATE(glbidx)
      DEALLOCATE(glb_phys_id_c)
      DEALLOCATE(glb_phys_id_e)
      DEALLOCATE(glb_phys_id_v)

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

  !-------------------------------------------------------------------------
  !>
  !!               Calculates local line/block indices l_idx, l_blk
  !!               from global line/block indices g_idx, g_blk
  !!               using the mapping in loc_index
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!

  !-------------------------------------------------------------------------
  !>
  !! Maps indices which point outside domain (returned as <0 by get_local_index)
  !! to valid ones which are from the same set of neighbors for the given
  !! cell/edge/vert
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Oct 2011
  !!
  SUBROUTINE remap_index(loc_index, l_idx, l_blk)

    INTEGER, INTENT(in) :: loc_index(:)
    INTEGER, INTENT(inout) :: l_idx, l_blk

    INTEGER :: j_l, j_g

    IF(l_idx>=0) RETURN ! Nothing to do

    j_g = idx_1d(-l_idx, l_blk) ! Global index

    ! Safety check only, global index should be in correct range
    IF(j_g < 1 .OR. j_g > UBOUND(loc_index,1)) CALL finish('remap_index','Invalid global index')

    j_l = loc_index(j_g)

    ! Safety check only, local index should be negative
    IF(j_l >= 0) CALL finish('remap_index','Invalid local index')

    ! Remap in the same way as previously in get_local_index (for compatibility only)
    j_l = MAX(ABS(j_l)-1,1)

    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE remap_index

!--------------------------------------------------------------
  SUBROUTINE complete_parallel_setup_oce(p_patch_2D)

   TYPE(t_patch) :: p_patch_2D(:)

    INTEGER :: jg, jgp

    DO jg = n_dom_start, n_dom

      jgp = p_patch_2D(jg)%parent_id

      ! Set communication patterns for boundary exchange
      CALL set_comm_pat_bound_exch(p_patch_2D(jg))

      ! Set communication patterns for gathering on proc 0
      CALL set_comm_pat_gather(p_patch_2D(jg))

      CALL set_owner_mask(p_patch_2D(jg))
     
     ! Fill the owner_local value
      ! this is done in the set_owner_mask
      ! CALL fill_owner_local(p_patch(jg))

      IF(jg == n_dom_start) THEN

        ! parent_idx/blk is set to 0 since it just doesn't exist,
        ! child_idx/blk is set to 0 since it makes sense only on the local parent
        p_patch_2D(jg)%cells%parent_idx = 0
        p_patch_2D(jg)%cells%parent_blk = 0
        p_patch_2D(jg)%cells%child_idx  = 0
        p_patch_2D(jg)%cells%child_blk  = 0
        p_patch_2D(jg)%edges%parent_idx = 0
        p_patch_2D(jg)%edges%parent_blk = 0
        p_patch_2D(jg)%edges%child_idx  = 0
        p_patch_2D(jg)%edges%child_blk  = 0

      ELSE
        CALL finish('complete_parallel_setup_oce','functionality with local parent patch not implemented')
!         CALL setup_comm_cpy_interpolation(p_patch_2D(jg), p_patch_2D(jgp))
!         CALL setup_comm_grf_interpolation(p_patch_2D(jg), p_patch_2D(jgp))
!         CALL setup_comm_ubc_interpolation(p_patch_2D(jg), p_patch_2D(jgp))
! 
!         CALL set_comm_pat_bound_exch(p_patch_local_parent(jg))
!         CALL set_comm_pat_gather(p_patch_local_parent(jg))
! 
!         CALL set_parent_child_relations(p_patch_local_parent(jg), p_patch(jg))
! 
!         CALL set_glb_loc_comm(p_patch(jgp), p_patch_local_parent(jg), &
!           &                   p_patch(jg)%parent_child_index)
! 
!         CALL set_owner_mask(p_patch_local_parent(jg))
      ENDIF

    ENDDO

    CALL set_patch_communicators(p_patch_2D)

  END SUBROUTINE complete_parallel_setup_oce

  !-----------------------------------------------------------------------------

  SUBROUTINE finalize_decomposition_oce(p_patch_2D)

    implicit none

    TYPE(t_patch) :: p_patch_2D(:)
    integer jg

    ! Remap indices in patches and local parents

    DO jg = n_dom_start, n_dom

      CALL remap_patch_indices(p_patch_2D(jg))

      IF(jg>n_dom_start) THEN
        CALL finish('finalize_decomposition_oce','functionality with local parent patch not implemented')
        CALL remap_patch_indices(p_patch_local_parent(jg))
      ENDIF

    ENDDO
                 
  END SUBROUTINE finalize_decomposition_oce

  !-----------------------------------------------------------------------------



END MODULE mo_complete_subdivision

