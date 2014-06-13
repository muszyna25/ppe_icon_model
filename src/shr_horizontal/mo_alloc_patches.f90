!>
!!
!!               The module <i>mo_alloc_domain</i>
!! provides routines to allocate and deallocate patches. The routines
!! included here have been extracted from mo_model_import_domain in order
!! to resolve USE dependencies that led to large compile times on some
!! platforms
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD, 2012-05-14
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


MODULE mo_alloc_patches
  !-------------------------------------------------------------------------

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: success, &
    & max_char_length,  &
    & min_rlcell, max_rlcell, &
    & min_rledge, max_rledge, &
    & min_rlvert, max_rlvert, &
    & max_dom
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_model_domain,       ONLY: t_patch, t_pre_patch
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info, &
    &                              init_glb2loc_index_lookup, &
    &                              t_glb2loc_index_lookup
  USE mo_parallel_config,    ONLY: nproma
  USE mo_grid_config,        ONLY: n_dom, n_dom_start, max_childdom, &
    & dynamics_grid_filename,   dynamics_parent_grid_id,  &
    & radiation_grid_filename, lplane
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords
  USE mo_master_nml,         ONLY: model_base_dir
  USE mo_mpi,                ONLY: my_process_is_mpi_seq

  IMPLICIT NONE

  PRIVATE

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: deallocate_patch
  PUBLIC :: deallocate_glb2loc_index_lookup
  PUBLIC :: destruct_patches
  PUBLIC :: allocate_basic_patch
  PUBLIC :: allocate_pre_patch
  PUBLIC :: deallocate_pre_patch
  PUBLIC :: allocate_remaining_patch
  ! PUBLIC :: allocate_patch
  PUBLIC :: set_patches_grid_filename

  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  SUBROUTINE set_patches_grid_filename( p_patch_pre )

    TYPE(t_pre_patch), TARGET, INTENT(inout) :: p_patch_pre(n_dom_start:)

    INTEGER :: jg
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    !-----------------------------------------------------------------------
    DO jg = n_dom_start, n_dom

      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
      IF (jg==0) THEN
        p_patch_pre(jg)%grid_filename = TRIM(with_keywords(keywords, radiation_grid_filename(1)))
      ELSE
        p_patch_pre(jg)%grid_filename = TRIM(with_keywords(keywords, dynamics_grid_filename(jg)))
      ENDIF
      !     write(0,*) jg, "grid_filename:",TRIM(p_patch_pre(jg)%grid_filename)
    ENDDO

  END SUBROUTINE set_patches_grid_filename
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of all memory that was allocated
  !!               to store the patch information of all levels.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !! Modification by Almut Gassmann, MPI-M (2007-04-04)
  !! - adaptations for new multiple patch structure, cleaning up
  !! Modification by Almut Gassmann, MPI-M (2007-04-19)
  !! - grid information belong here
  !! Modification by Almut Gassmann, MPI-M (2008-10-30)
  !! - add Coriolis destruction
  !!
  SUBROUTINE deallocate_patch( p_patch, lddmode )

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch
    LOGICAL,     OPTIONAL, INTENT(in)    :: lddmode  ! Skip fields not allocated in ldump_dd mode

    !local variables
    CHARACTER(*), PARAMETER :: routine ="mo_alloc_patches:deallocate_patch"
    INTEGER :: ist
    LOGICAL :: l_ddmode
    !-----------------------------------------------------------------------
    !
    IF (PRESENT(lddmode)) THEN
      l_ddmode = lddmode
    ELSE
      l_ddmode = .FALSE.
    ENDIF

    ! Deallocate grid information
    !
    ! Include deallocation of basic patch variables
! GZ: this directive should not be needed here, but the NEC compiler crashes otherwise (internal error in optimization phase)
!CDIR NOIEXPAND
    CALL deallocate_basic_patch(p_patch)

    call deallocate_decomp_info(p_patch%cells%decomp_info)
    call deallocate_decomp_info(p_patch%edges%decomp_info)
    call deallocate_decomp_info(p_patch%verts%decomp_info)

    DEALLOCATE(p_patch%cells%phys_id,    stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell physical ID failed')
    ENDIF

    IF (l_ddmode) RETURN

    DEALLOCATE( p_patch%cells%edge_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell edge_orientation failed')
    ENDIF
    DEALLOCATE( p_patch%cells%area,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell area failed')
    ENDIF
    DEALLOCATE( p_patch%cells%f_c,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell f_c failed')
    ENDIF
    !
    !
    DEALLOCATE(p_patch%edges%phys_id,    stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge physical ID failed')
    ENDIF
    DEALLOCATE( p_patch%edges%cell_idx,  &
      & p_patch%edges%cell_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge cell index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%vertex_idx,  &
      & p_patch%edges%vertex_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge vertex index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%system_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge system orientation failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_idx,  &
      & p_patch%edges%quad_blk,  stat=ist)
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge quad index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge quad orientation failed')
    ENDIF
    DEALLOCATE( p_patch%edges%butterfly_idx,  &
      & p_patch%edges%butterfly_blk,  stat=ist)
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge butterfly index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%center,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge center failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal_normal failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual_normal failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal_cell,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal_normal_cell failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal_cell,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual_normal_cell failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal_vert,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal_normal_vert failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal_vert,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual_normal_vert failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_primal_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for inverse patch edge primal edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_dual_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for inverse patch edge dual edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%edge_vert_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for edge_vert_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_vert_vert_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for inverse vert_vert_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%edge_cell_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for edge_cell_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%area_edge,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge area_edge failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_area,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge quad area failed')
    ENDIF
    DEALLOCATE( p_patch%edges%f_e,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge f_e failed')
    ENDIF
    !
    !
    DEALLOCATE( p_patch%verts%neighbor_idx,  &
      & p_patch%verts%neighbor_blk,  &
      & p_patch%verts%phys_id,       stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex neighbor index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%cell_idx,  &
      & p_patch%verts%cell_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex cell index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%edge_idx,  &
      & p_patch%verts%edge_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex edge index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%edge_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex edge orientation failed')
    ENDIF
    DEALLOCATE( p_patch%verts%num_edges,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex number of edges failed')
    ENDIF
    DEALLOCATE( p_patch%verts%dual_area,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex dual area failed')
    ENDIF
    DEALLOCATE( p_patch%verts%f_v, stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex f_v failed')
    ENDIF


    CALL deallocate_patch_cartesian( p_patch )

  CONTAINS

    SUBROUTINE deallocate_decomp_info( decomp_info )

      TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info

      INTEGER :: ist

      DEALLOCATE(decomp_info%glb_index, &
        &        decomp_info%owner_local, &
        &        decomp_info%owner_mask, &
        &        decomp_info%decomp_domain, stat=ist )
      IF(ist/=success)THEN
        CALL finish  (routine,  'deallocate in deallocate_decomp_info failed')
      ENDIF
      CALL deallocate_glb2loc_index_lookup(decomp_info%glb2loc_index)
    END SUBROUTINE deallocate_decomp_info

  END SUBROUTINE deallocate_patch
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

  SUBROUTINE allocate_patch_cartesian( p_patch )
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch

    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_alloc_patches:allocate_patch_cartesian'

    ALLOCATE( p_patch%edges%primal_cart_normal(nproma,p_patch%nblks_e) )
    p_patch%edges%primal_cart_normal(:,:)%x(1) = 0._wp
    p_patch%edges%primal_cart_normal(:,:)%x(2) = 0._wp
    p_patch%edges%primal_cart_normal(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%edges%dual_cart_normal(nproma,p_patch%nblks_e) )
    p_patch%edges%dual_cart_normal(:,:)%x(1) = 0._wp
    p_patch%edges%dual_cart_normal(:,:)%x(2) = 0._wp
    p_patch%edges%dual_cart_normal(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%cells%cartesian_center(nproma,p_patch%alloc_cell_blocks) )
    p_patch%cells%cartesian_center(:,:)%x(1) = 0._wp
    p_patch%cells%cartesian_center(:,:)%x(2) = 0._wp
    p_patch%cells%cartesian_center(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%edges%cartesian_center(nproma,p_patch%nblks_e) )
    p_patch%edges%cartesian_center(:,:)%x(1) = 0._wp
    p_patch%edges%cartesian_center(:,:)%x(2) = 0._wp
    p_patch%edges%cartesian_center(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%edges%cartesian_dual_middle(nproma,p_patch%nblks_e) )
    p_patch%edges%cartesian_dual_middle(:,:)%x(1) = 0._wp
    p_patch%edges%cartesian_dual_middle(:,:)%x(2) = 0._wp
    p_patch%edges%cartesian_dual_middle(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%verts%cartesian(nproma,p_patch%nblks_v) )
    p_patch%verts%cartesian(:,:)%x(1) = 0._wp
    p_patch%verts%cartesian(:,:)%x(2) = 0._wp
    p_patch%verts%cartesian(:,:)%x(3) = 0._wp

  END SUBROUTINE allocate_patch_cartesian
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_patch_cartesian( p_patch )
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch

    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_alloc_patches:deallocate_patch_cartesian'

    DEALLOCATE( p_patch%edges%primal_cart_normal, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%primal_cart_normal)")
    ENDIF

    DEALLOCATE( p_patch%edges%dual_cart_normal, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%dual_cart_normal)")
    ENDIF

    DEALLOCATE( p_patch%cells%cartesian_center, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%cells%cartesian_center)")
    ENDIF

    DEALLOCATE( p_patch%edges%cartesian_center, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%cartesian_center)")
    ENDIF

    DEALLOCATE( p_patch%edges%cartesian_dual_middle, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%cartesian_dual_middle)")
    ENDIF

    DEALLOCATE( p_patch%verts%cartesian, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%verts%cartesian)")
    ENDIF

  END SUBROUTINE deallocate_patch_cartesian
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> destruct_patches: Calls destruct_patch in a loop
  SUBROUTINE destruct_patches( p_patch )

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch(n_dom_start:)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_alloc_patches:destruct_patches'

    !local variables
    INTEGER :: jg, istart, iend
    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    istart = LBOUND(p_patch, 1)
    iend   = UBOUND(p_patch, 1)

    grid_level_loop: DO jg = istart, iend
      CALL deallocate_patch(p_patch(jg))
    END DO grid_level_loop

    CALL message (routine, 'destruct_patches finished')

  END SUBROUTINE destruct_patches
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !> Allocates all arrays in a basic patch
  !! These are the arrays which have to be read in full size on every PE
  !! since they are needed for domain decomposition.
  !!
  !! @par Revision History
  !! Initial version (split out from read_patch) Rainer Johanni, Oct. 2010
  !! Split into allocate_basic_patch/allocate_remaining_patch, Rainer Johanni, Nov. 2011
  !!
  SUBROUTINE allocate_basic_patch(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch

    INTEGER :: max_childdom

    ! Please note: The following variables in the patch MUST already be set:
    ! - alloc_cell_blocks
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells
    ! - n_patch_edges
    ! - n_patch_verts
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom

    p_patch%cell_type = p_patch%cells%max_connectivity
    max_childdom = p_patch%max_childdom
    !
    ! !grid cells
    !
    p_patch%cells%dummy_cell_block = 0
    p_patch%cells%dummy_cell_index = 0
    ALLOCATE( p_patch%cells%num_edges(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%parent_idx(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%parent_blk(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%pc_idx(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%child_idx(nproma,p_patch%alloc_cell_blocks,4) )
    ALLOCATE( p_patch%cells%child_blk(nproma,p_patch%alloc_cell_blocks,4) )
    ALLOCATE( p_patch%cells%child_id(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%phys_id(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%neighbor_idx(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%neighbor_blk(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%edge_idx(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%edge_blk(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%vertex_idx(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%vertex_blk(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%center(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%refin_ctrl(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%start_idx(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%end_idx(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%start_blk(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%end_blk(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%start_index(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch%cells%end_index(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch%cells%start_block(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch%cells%end_block(min_rlcell:max_rlcell) )

    !
    ! !grid edges
    !
    ALLOCATE( p_patch%edges%parent_idx(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%parent_blk(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%pc_idx(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%child_idx(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%child_blk(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%child_id(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%refin_ctrl(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%start_idx(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%end_idx(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%start_blk(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%end_blk(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%cell_idx(nproma,p_patch%nblks_e,2) )
    ALLOCATE( p_patch%edges%cell_blk(nproma,p_patch%nblks_e,2) )
    ALLOCATE( p_patch%edges%vertex_idx(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%vertex_blk(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%start_index(min_rledge:max_rledge) )
    ALLOCATE( p_patch%edges%end_index(min_rledge:max_rledge) )
    ALLOCATE( p_patch%edges%start_block(min_rledge:max_rledge) )
    ALLOCATE( p_patch%edges%end_block(min_rledge:max_rledge) )

    !
    ! !grid verts
    !
    ALLOCATE( p_patch%verts%vertex(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%refin_ctrl(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%start_idx(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%end_idx(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%start_blk(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%end_blk(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%cell_idx(nproma,p_patch%nblks_v,6) )
    ALLOCATE( p_patch%verts%cell_blk(nproma,p_patch%nblks_v,6) )
    ALLOCATE( p_patch%verts%num_edges(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%start_index(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch%verts%end_index(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch%verts%start_block(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch%verts%end_block(min_rlvert:max_rlvert) )
    ! Set all newly allocated arrays to 0

    p_patch%cells%num_edges = 0
    p_patch%cells%parent_idx = 0
    p_patch%cells%parent_blk = 0
    p_patch%cells%pc_idx = 0
    p_patch%cells%child_idx = 0
    p_patch%cells%child_blk = 0
    p_patch%cells%child_id = 0
    p_patch%cells%phys_id = 0
    p_patch%cells%neighbor_idx = 0
    p_patch%cells%neighbor_blk = 0
    p_patch%cells%edge_idx = 0
    p_patch%cells%edge_blk = 0
    p_patch%cells%vertex_idx = 0
    p_patch%cells%vertex_blk = 0
    p_patch%cells%center(:,:)%lon = 0._wp
    p_patch%cells%center(:,:)%lat = 0._wp
    p_patch%cells%refin_ctrl = 0
    p_patch%cells%start_idx = 0
    p_patch%cells%end_idx = 0
    p_patch%cells%start_blk = 0
    p_patch%cells%end_blk = 0
    p_patch%cells%start_index = 0
    p_patch%cells%end_index = 0
    p_patch%cells%start_block = 0
    p_patch%cells%end_block = 0

    p_patch%edges%parent_idx = 0
    p_patch%edges%parent_blk = 0
    p_patch%edges%pc_idx = 0
    p_patch%edges%child_idx = 0
    p_patch%edges%child_blk = 0
    p_patch%edges%child_id = 0
    p_patch%edges%refin_ctrl = 0
    p_patch%edges%start_idx = 0
    p_patch%edges%end_idx = 0
    p_patch%edges%start_blk = 0
    p_patch%edges%end_blk = 0
    p_patch%edges%start_index = 0
    p_patch%edges%end_index = 0
    p_patch%edges%start_block = 0
    p_patch%edges%end_block = 0

    p_patch%verts%vertex(:,:)%lon = 0._wp
    p_patch%verts%vertex(:,:)%lat = 0._wp
    p_patch%verts%refin_ctrl = 0
    p_patch%verts%start_idx = 0
    p_patch%verts%end_idx = 0
    p_patch%verts%start_blk = 0
    p_patch%verts%end_blk = 0
    p_patch%verts%cell_idx = 0
    p_patch%verts%cell_blk = 0
    p_patch%verts%num_edges = 0
    p_patch%verts%start_index = 0
    p_patch%verts%end_index = 0
    p_patch%verts%start_block = 0
    p_patch%verts%end_block = 0

  END SUBROUTINE allocate_basic_patch
  
  !-------------------------------------------------------------------------
  !> Allocates all arrays in a basic patch
  !! These are the arrays which have to be read in full size on every PE
  !! since they are needed for domain decomposition.
  !!
  !! @par Revision History
  !! Initial version (split out from read_patch) Rainer Johanni, Oct. 2010
  !! Split into allocate_basic_patch/allocate_remaining_patch, Rainer Johanni, Nov. 2011
  !!
  SUBROUTINE allocate_pre_patch(p_patch_pre)

    TYPE(t_pre_patch), INTENT(inout) :: p_patch_pre

    INTEGER :: max_childdom

    ! Please note: The following variables in the patch MUST already be set:
    ! - alloc_cell_blocks
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom

    p_patch_pre%cell_type = p_patch_pre%cells%max_connectivity
    max_childdom = p_patch_pre%max_childdom
    !
    ! !grid cells
    !
    ALLOCATE( p_patch_pre%cells%num_edges(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%parent(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%pc_idx(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%child(p_patch_pre%n_patch_cells_g,4) )
    ALLOCATE( p_patch_pre%cells%child_id(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%phys_id(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%neighbor(p_patch_pre%n_patch_cells_g,p_patch_pre%cell_type) )
    ALLOCATE( p_patch_pre%cells%edge(p_patch_pre%n_patch_cells_g,p_patch_pre%cell_type) )
    ALLOCATE( p_patch_pre%cells%vertex(p_patch_pre%n_patch_cells_g,p_patch_pre%cell_type) )
    ALLOCATE( p_patch_pre%cells%center(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%refin_ctrl(p_patch_pre%n_patch_cells_g) )
    ALLOCATE( p_patch_pre%cells%start_index(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch_pre%cells%end_index(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch_pre%cells%start_block(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch_pre%cells%end_block(min_rlcell:max_rlcell) )

    !
    ! !grid edges
    !
    ALLOCATE( p_patch_pre%edges%parent(p_patch_pre%n_patch_edges_g) )
    ALLOCATE( p_patch_pre%edges%pc_idx(p_patch_pre%n_patch_edges_g) )
    ALLOCATE( p_patch_pre%edges%child(p_patch_pre%n_patch_edges_g,4) )
    ALLOCATE( p_patch_pre%edges%child_id(p_patch_pre%n_patch_edges_g) )
    ALLOCATE( p_patch_pre%edges%refin_ctrl(p_patch_pre%n_patch_edges_g) )
    ALLOCATE( p_patch_pre%edges%cell_idx(nproma,p_patch_pre%nblks_e,2) )
    ALLOCATE( p_patch_pre%edges%cell_blk(nproma,p_patch_pre%nblks_e,2) )
    ALLOCATE( p_patch_pre%edges%vertex_idx(nproma,p_patch_pre%nblks_e,4) )
    ALLOCATE( p_patch_pre%edges%vertex_blk(nproma,p_patch_pre%nblks_e,4) )
    ALLOCATE( p_patch_pre%edges%start_index(min_rledge:max_rledge) )
    ALLOCATE( p_patch_pre%edges%end_index(min_rledge:max_rledge) )
    ALLOCATE( p_patch_pre%edges%start_block(min_rledge:max_rledge) )
    ALLOCATE( p_patch_pre%edges%end_block(min_rledge:max_rledge) )

    !
    ! !grid verts
    !
    ALLOCATE( p_patch_pre%verts%vertex(nproma,p_patch_pre%nblks_v) )
    ALLOCATE( p_patch_pre%verts%refin_ctrl(p_patch_pre%n_patch_verts_g) )
    ALLOCATE( p_patch_pre%verts%cell_idx(nproma,p_patch_pre%nblks_v,6) )
    ALLOCATE( p_patch_pre%verts%cell_blk(nproma,p_patch_pre%nblks_v,6) )
    ALLOCATE( p_patch_pre%verts%num_edges(nproma,p_patch_pre%nblks_v) )
    ALLOCATE( p_patch_pre%verts%start_index(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch_pre%verts%end_index(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch_pre%verts%start_block(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch_pre%verts%end_block(min_rlvert:max_rlvert) )
    ! Set all newly allocated arrays to 0

    p_patch_pre%cells%num_edges = 0
    p_patch_pre%cells%parent = 0
    p_patch_pre%cells%pc_idx = 0
    p_patch_pre%cells%child = 0
    p_patch_pre%cells%child_id = 0
    p_patch_pre%cells%phys_id = 0
    p_patch_pre%cells%neighbor = 0
    p_patch_pre%cells%edge = 0
    p_patch_pre%cells%vertex = 0
    p_patch_pre%cells%center(:)%lon = 0._wp
    p_patch_pre%cells%center(:)%lat = 0._wp
    p_patch_pre%cells%refin_ctrl = 0
    p_patch_pre%cells%start_index = 0
    p_patch_pre%cells%end_index = 0
    p_patch_pre%cells%start_block = 0
    p_patch_pre%cells%end_block = 0

    p_patch_pre%edges%parent = 0
    p_patch_pre%edges%pc_idx = 0
    p_patch_pre%edges%child = 0
    p_patch_pre%edges%child_id = 0
    p_patch_pre%edges%refin_ctrl = 0
    p_patch_pre%edges%start_index = 0
    p_patch_pre%edges%end_index = 0
    p_patch_pre%edges%start_block = 0
    p_patch_pre%edges%end_block = 0

    p_patch_pre%verts%vertex(:,:)%lon = 0._wp
    p_patch_pre%verts%vertex(:,:)%lat = 0._wp
    p_patch_pre%verts%refin_ctrl = 0
    p_patch_pre%verts%cell_idx = 0
    p_patch_pre%verts%cell_blk = 0
    p_patch_pre%verts%num_edges = 0
    p_patch_pre%verts%start_index = 0
    p_patch_pre%verts%end_index = 0
    p_patch_pre%verts%start_block = 0
    p_patch_pre%verts%end_block = 0

  END SUBROUTINE allocate_pre_patch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Deallocates all arrays in a basic patch
  !!
  !! @par Revision History
  !! Initial version, Rainer Johanni, Nov. 2011
  !!
  SUBROUTINE deallocate_basic_patch(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch
    !
    ! !grid cells
    !
    DEALLOCATE( p_patch%cells%num_edges )
    DEALLOCATE( p_patch%cells%parent_idx )
    DEALLOCATE( p_patch%cells%parent_blk )
    DEALLOCATE( p_patch%cells%pc_idx )
    DEALLOCATE( p_patch%cells%child_idx )
    DEALLOCATE( p_patch%cells%child_blk )
    DEALLOCATE( p_patch%cells%child_id )
    DEALLOCATE( p_patch%cells%neighbor_idx )
    DEALLOCATE( p_patch%cells%neighbor_blk )
    DEALLOCATE( p_patch%cells%edge_idx )
    DEALLOCATE( p_patch%cells%edge_blk )
    DEALLOCATE( p_patch%cells%vertex_idx )
    DEALLOCATE( p_patch%cells%vertex_blk )
    DEALLOCATE( p_patch%cells%center )
    DEALLOCATE( p_patch%cells%refin_ctrl )
    DEALLOCATE( p_patch%cells%start_idx )
    DEALLOCATE( p_patch%cells%end_idx )
    DEALLOCATE( p_patch%cells%start_blk )
    DEALLOCATE( p_patch%cells%end_blk )
    DEALLOCATE( p_patch%cells%start_index )
    DEALLOCATE( p_patch%cells%end_index )
    DEALLOCATE( p_patch%cells%start_block )
    DEALLOCATE( p_patch%cells%end_block )
    !
    ! !grid edges
    !
    DEALLOCATE( p_patch%edges%parent_idx )
    DEALLOCATE( p_patch%edges%parent_blk )
    DEALLOCATE( p_patch%edges%pc_idx )
    DEALLOCATE( p_patch%edges%child_idx )
    DEALLOCATE( p_patch%edges%child_blk )
    DEALLOCATE( p_patch%edges%child_id )
    DEALLOCATE( p_patch%edges%refin_ctrl )
    DEALLOCATE( p_patch%edges%start_idx )
    DEALLOCATE( p_patch%edges%end_idx )
    DEALLOCATE( p_patch%edges%start_blk )
    DEALLOCATE( p_patch%edges%end_blk )
    DEALLOCATE( p_patch%edges%start_index )
    DEALLOCATE( p_patch%edges%end_index )
    DEALLOCATE( p_patch%edges%start_block )
    DEALLOCATE( p_patch%edges%end_block )
    ! DEALLOCATE( p_patch%edges%cell_idx )
    ! DEALLOCATE( p_patch%edges%cell_blk )
    ! DEALLOCATE( p_patch%edges%vertex_idx )
    ! DEALLOCATE( p_patch%edges%vertex_blk )
    !
    ! !grid verts
    !
    DEALLOCATE( p_patch%verts%vertex )
    DEALLOCATE( p_patch%verts%refin_ctrl )
    DEALLOCATE( p_patch%verts%start_idx )
    DEALLOCATE( p_patch%verts%end_idx )
    DEALLOCATE( p_patch%verts%start_blk )
    DEALLOCATE( p_patch%verts%end_blk )
    DEALLOCATE( p_patch%verts%start_index )
    DEALLOCATE( p_patch%verts%end_index )
    DEALLOCATE( p_patch%verts%start_block )
    DEALLOCATE( p_patch%verts%end_block )
    ! DEALLOCATE( p_patch%verts%cell_idx )
    ! DEALLOCATE( p_patch%verts%cell_blk )
    ! DEALLOCATE( p_patch%verts%num_edges )

  END SUBROUTINE deallocate_basic_patch

  !-------------------------------------------------------------------------
  !> Deallocates all arrays in a basic patch
  !!
  !! @par Revision History
  !! Initial version, Rainer Johanni, Nov. 2011
  !!
  SUBROUTINE deallocate_pre_patch(p_patch_pre)

    TYPE(t_pre_patch), INTENT(inout) :: p_patch_pre
    !
    ! !grid cells
    !
    DEALLOCATE( p_patch_pre%cells%num_edges )
    DEALLOCATE( p_patch_pre%cells%parent )
    DEALLOCATE( p_patch_pre%cells%pc_idx )
    DEALLOCATE( p_patch_pre%cells%child )
    DEALLOCATE( p_patch_pre%cells%child_id )
    DEALLOCATE( p_patch_pre%cells%neighbor )
    DEALLOCATE( p_patch_pre%cells%edge )
    DEALLOCATE( p_patch_pre%cells%vertex )
    DEALLOCATE( p_patch_pre%cells%center )
    DEALLOCATE( p_patch_pre%cells%refin_ctrl )
    DEALLOCATE( p_patch_pre%cells%start_index )
    DEALLOCATE( p_patch_pre%cells%end_index )
    DEALLOCATE( p_patch_pre%cells%start_block )
    DEALLOCATE( p_patch_pre%cells%end_block )
    !
    ! !grid edges
    !
    DEALLOCATE( p_patch_pre%edges%parent )
    DEALLOCATE( p_patch_pre%edges%pc_idx )
    DEALLOCATE( p_patch_pre%edges%child )
    DEALLOCATE( p_patch_pre%edges%child_id )
    DEALLOCATE( p_patch_pre%edges%refin_ctrl )
    DEALLOCATE( p_patch_pre%edges%start_index )
    DEALLOCATE( p_patch_pre%edges%end_index )
    DEALLOCATE( p_patch_pre%edges%start_block )
    DEALLOCATE( p_patch_pre%edges%end_block )
    !
    ! !grid verts
    !
    DEALLOCATE( p_patch_pre%verts%vertex )
    DEALLOCATE( p_patch_pre%verts%refin_ctrl )
    DEALLOCATE( p_patch_pre%verts%start_index )
    DEALLOCATE( p_patch_pre%verts%end_index )
    DEALLOCATE( p_patch_pre%verts%start_block )
    DEALLOCATE( p_patch_pre%verts%end_block )

  END SUBROUTINE deallocate_pre_patch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Allocates remaining arrays in a full patch.
  !! These are the arrays which have not been allocated in allocate_basic_patch
  !!
  !! @par Revision History
  !! Initial version (split out from read_patch) Rainer Johanni, Oct. 2010
  !! Split into allocate_basic_patch/allocate_remaining_patch, Rainer Johanni, Nov. 2011
  !!
  SUBROUTINE allocate_remaining_patch(p_patch,iopmode)

    TYPE(t_patch), INTENT(inout) :: p_patch
    INTEGER,       INTENT(in)    :: iopmode

    ! operation modes: iopmode=1: allocate all remaining patch arrays
    !                  iopmode=2: allocate only parallelization-related arrays
    !                  iopmode=3: allocate all but the parallelization-related arrays

    ! Please note: The following variables in the patch MUST already be set:
    ! - alloc_cell_blocks
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells
    ! - n_patch_edges
    ! - n_patch_verts
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom

    !
    ! !grid cells
    !
    IF (iopmode /= 2) THEN
      IF (.NOT. ALLOCATED(p_patch%cells%phys_id)) ALLOCATE( p_patch%cells%phys_id(nproma,p_patch%alloc_cell_blocks) )
      ALLOCATE( p_patch%cells%edge_orientation(nproma,p_patch%alloc_cell_blocks,p_patch%cell_type) )
      ALLOCATE( p_patch%cells%area(nproma,p_patch%alloc_cell_blocks) )
      ALLOCATE( p_patch%cells%f_c(nproma,p_patch%alloc_cell_blocks) )
    ENDIF

    !
    ! !grid edges
    !
    IF (iopmode /= 2) THEN
      ALLOCATE( p_patch%edges%phys_id(nproma,p_patch%nblks_e) )
      IF (ALLOCATED(p_patch%edges%cell_idx)) DEALLOCATE(p_patch%edges%cell_idx)
      ALLOCATE( p_patch%edges%cell_idx(nproma,p_patch%nblks_e,2) )
      IF (ALLOCATED(p_patch%edges%cell_blk)) DEALLOCATE(p_patch%edges%cell_blk)
      ALLOCATE( p_patch%edges%cell_blk(nproma,p_patch%nblks_e,2) )
      IF (ALLOCATED(p_patch%edges%vertex_idx)) DEALLOCATE(p_patch%edges%vertex_idx)
      ALLOCATE( p_patch%edges%vertex_idx(nproma,p_patch%nblks_e,4) )
      IF (ALLOCATED(p_patch%edges%vertex_blk)) DEALLOCATE(p_patch%edges%vertex_blk)
      ALLOCATE( p_patch%edges%vertex_blk(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%system_orientation(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%quad_idx(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%quad_blk(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%butterfly_idx(nproma,p_patch%nblks_e,2,2) )
      ALLOCATE( p_patch%edges%butterfly_blk(nproma,p_patch%nblks_e,2,2) )
      ALLOCATE( p_patch%edges%quad_orientation(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%center(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%primal_normal(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%dual_normal(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%primal_normal_cell(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%dual_normal_cell(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%primal_normal_vert(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%dual_normal_vert(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%primal_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%inv_primal_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%dual_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%inv_dual_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%edge_vert_length(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%inv_vert_vert_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%edge_cell_length(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%area_edge(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%quad_area(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%f_e(nproma,p_patch%nblks_e) )
    ENDIF

    !
    ! !grid verts
    !
    IF (iopmode /= 2) THEN
      ALLOCATE( p_patch%verts%phys_id(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%neighbor_idx(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%neighbor_blk(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      IF (ALLOCATED(p_patch%verts%cell_idx)) DEALLOCATE(p_patch%verts%cell_idx)
      ALLOCATE( p_patch%verts%cell_idx(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      IF (ALLOCATED(p_patch%verts%cell_blk)) DEALLOCATE(p_patch%verts%cell_blk)
      ALLOCATE( p_patch%verts%cell_blk(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%edge_idx(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%edge_blk(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%edge_orientation(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      IF (ALLOCATED(p_patch%verts%num_edges)) &
           DEALLOCATE(p_patch%verts%num_edges)
      ALLOCATE( p_patch%verts%num_edges(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%dual_area(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%f_v(nproma,p_patch%nblks_v) )
    ENDIF

    !
    ! ! decomp_info for cells, edges and verts
    !
    IF (iopmode /= 3) THEN

      CALL allocate_decomp_info(p_patch%cells%decomp_info, &
        &                  p_patch%n_patch_cells, p_patch%n_patch_cells_g, &
        &                  p_patch%alloc_cell_blocks)
      CALL allocate_decomp_info(p_patch%edges%decomp_info, &
        &                  p_patch%n_patch_edges, p_patch%n_patch_edges_g, &
        &                  p_patch%nblks_e)
      CALL allocate_decomp_info(p_patch%verts%decomp_info, &
        &                  p_patch%n_patch_verts, p_patch%n_patch_verts_g, &
        &                  p_patch%nblks_v)

      CALL init_decomp_info(p_patch%cells%decomp_info, p_patch%npromz_c, &
        &                   p_patch%alloc_cell_blocks, p_patch%n_patch_cells_g)
      CALL init_decomp_info(p_patch%edges%decomp_info, p_patch%npromz_e, &
        &                   p_patch%nblks_e, p_patch%n_patch_edges_g)
      CALL init_decomp_info(p_patch%verts%decomp_info, p_patch%npromz_v, &
        &                   p_patch%nblks_v, p_patch%n_patch_verts_g)
    END IF

    ! Set all newly allocated arrays to 0

    IF (iopmode /= 2) THEN
      p_patch%cells%phys_id = 0
      p_patch%cells%edge_orientation = 0._wp
      p_patch%cells%area = 0._wp
      p_patch%cells%f_c = 0._wp

      p_patch%edges%phys_id = 0
      p_patch%edges%cell_idx = 0
      p_patch%edges%cell_blk = 0
      p_patch%edges%vertex_idx = 0
      p_patch%edges%vertex_blk = 0
      p_patch%edges%system_orientation = 0._wp
      p_patch%edges%quad_idx = 0
      p_patch%edges%quad_blk = 0
      p_patch%edges%butterfly_idx = 0
      p_patch%edges%butterfly_blk = 0
      p_patch%edges%quad_orientation = 0._wp
      p_patch%edges%center(:,:)%lon = 0._wp
      p_patch%edges%center(:,:)%lat = 0._wp
      p_patch%edges%primal_normal(:,:)%v1 = 0._wp
      p_patch%edges%primal_normal(:,:)%v2 = 0._wp
      p_patch%edges%dual_normal(:,:)%v1 = 0._wp
      p_patch%edges%dual_normal(:,:)%v2 = 0._wp
      p_patch%edges%primal_normal_cell(:,:,:)%v1 = 0._wp
      p_patch%edges%primal_normal_cell(:,:,:)%v2 = 0._wp
      p_patch%edges%dual_normal_cell(:,:,:)%v1 = 0._wp
      p_patch%edges%dual_normal_cell(:,:,:)%v2 = 0._wp
      p_patch%edges%primal_normal_vert(:,:,:)%v1 = 0._wp
      p_patch%edges%primal_normal_vert(:,:,:)%v2 = 0._wp
      p_patch%edges%dual_normal_vert(:,:,:)%v1 = 0._wp
      p_patch%edges%dual_normal_vert(:,:,:)%v2 = 0._wp
      p_patch%edges%primal_edge_length = 0._wp
      p_patch%edges%inv_primal_edge_length = 0._wp
      p_patch%edges%dual_edge_length = 0._wp
      p_patch%edges%inv_dual_edge_length = 0._wp
      p_patch%edges%edge_vert_length = 0._wp
      p_patch%edges%inv_vert_vert_length = 0._wp
      p_patch%edges%edge_cell_length = 0._wp
      p_patch%edges%area_edge = 0._wp
      p_patch%edges%quad_area = 0._wp
      p_patch%edges%f_e = 0._wp

      p_patch%verts%phys_id = 0
      p_patch%verts%neighbor_idx = 0
      p_patch%verts%neighbor_blk = 0
      p_patch%verts%cell_idx = 0
      p_patch%verts%cell_blk = 0
      p_patch%verts%edge_idx = 0
      p_patch%verts%edge_blk = 0
      p_patch%verts%edge_orientation = 0._wp
      p_patch%verts%num_edges = 0
      p_patch%verts%dual_area = 0._wp
      p_patch%verts%f_v = 0._wp
    ENDIF

    IF (iopmode /= 2) CALL allocate_patch_cartesian( p_patch )

  CONTAINS

    SUBROUTINE allocate_decomp_info(decomp_info, n, n_g, n_blk)
      TYPE(t_grid_domain_decomp_info), INTENT(OUT) :: decomp_info
      INTEGER, INTENT(IN) :: n, n_g, n_blk

      ALLOCATE( decomp_info%decomp_domain(nproma,n_blk) )
      ! aliasing the halo_level to decomp_domain
      decomp_info%halo_level => decomp_info%decomp_domain
      ALLOCATE( decomp_info%owner_mask(nproma,n_blk) )
      ALLOCATE( decomp_info%glb_index(n) )
      ALLOCATE( decomp_info%owner_local(n))
    END SUBROUTINE allocate_decomp_info

    SUBROUTINE init_decomp_info(decomp_info, npromz, n_blk, n_g)
      TYPE(t_grid_domain_decomp_info), INTENT(INOUT) :: decomp_info
      INTEGER, INTENT(IN) :: npromz, n_blk, n_g

      INTEGER :: j

      ! Set values which are needed for parallel runs
      ! to the correct values for a single patch owner
      ! decomp_domain/owner mask:
      ! Everywhere 0 or .true. with the exception of unused entries

      decomp_info%decomp_domain = 0
      decomp_info%owner_mask = .TRUE.

      decomp_info%decomp_domain(npromz+1:nproma,n_blk) = -1
      decomp_info%owner_mask(npromz+1:nproma,n_blk) = .FALSE.

      ! The following arrays are currently never needed for non parallel runs,
      ! we set them nonetheless.
      ! MoHa: ...do they really need to be initialised?

      DO j = 1, SIZE(decomp_info%glb_index(:))
        decomp_info%glb_index(j) = j
      ENDDO
      decomp_info%owner_local(:) = 0

      CALL init_glb2loc_index_lookup(decomp_info%glb2loc_index, n_g)

    END SUBROUTINE init_decomp_info

  END SUBROUTINE allocate_remaining_patch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Allocates complete patch by calling allocate_basic_patch/allocate_remaining_patch
  !!
  !! @par Revision History
  !! Initial version, Rainer Johanni, Dec. 2011
  !!
  SUBROUTINE allocate_patch(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch

    CALL allocate_basic_patch(p_patch)
    CALL allocate_remaining_patch(p_patch,1)

  END SUBROUTINE allocate_patch
  !-------------------------------------------------------------------------


END MODULE mo_alloc_patches


