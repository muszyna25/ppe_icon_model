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
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_grid_config,        ONLY: n_dom, n_dom_start, max_childdom, &
    & dynamics_grid_filename,   dynamics_parent_grid_id,  &
    & radiation_grid_filename,  global_cell_type, lplane
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords
  USE mo_master_nml,         ONLY: model_base_dir
  USE mo_mpi,                ONLY: my_process_is_mpi_seq
    
  IMPLICIT NONE
  
  PRIVATE
  
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: deallocate_patch
  PUBLIC :: destruct_patches
  PUBLIC :: allocate_basic_patch
  PUBLIC :: deallocate_basic_patch
  PUBLIC :: allocate_remaining_patch
  PUBLIC :: allocate_patch
  PUBLIC :: set_patches_grid_filename
    
  !-------------------------------------------------------------------------
  
CONTAINS
  
  
  !-------------------------------------------------------------------------
  SUBROUTINE set_patches_grid_filename( p_patch )
    
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch(n_dom_start:)
    
    INTEGER :: jg
    TYPE (t_keyword_list), POINTER :: keywords => NULL()
    
    !-----------------------------------------------------------------------
    DO jg = n_dom_start, n_dom
      
      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
      IF (jg==0) THEN
        p_patch(jg)%grid_filename = TRIM(with_keywords(keywords, radiation_grid_filename(1)))
      ELSE
        p_patch(jg)%grid_filename = TRIM(with_keywords(keywords, dynamics_grid_filename(jg)))
      ENDIF
      !     write(0,*) jg, "grid_filename:",TRIM(p_patch(jg)%grid_filename)
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

    DEALLOCATE( p_patch%cells%decomp_domain,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell decomp_domain failed')
    ENDIF
    DEALLOCATE( p_patch%cells%owner_mask,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell owner_mask failed')
    ENDIF
    DEALLOCATE( p_patch%cells%glb_index,  &
      & p_patch%cells%loc_index,  &
      & p_patch%cells%owner_local,  &
      & p_patch%cells%owner_g,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell data failed')
    ENDIF

    DEALLOCATE( p_patch%edges%decomp_domain,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge decomp_domain failed')
    ENDIF
    DEALLOCATE( p_patch%edges%owner_mask,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge owner_mask failed')
    ENDIF
    DEALLOCATE( p_patch%edges%glb_index,  &
      & p_patch%edges%loc_index,  &
      & p_patch%edges%owner_local,  &
      & p_patch%edges%owner_g,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge data failed')
    ENDIF

    DEALLOCATE( p_patch%verts%decomp_domain,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vert decomp_domain failed')
    ENDIF
    DEALLOCATE( p_patch%verts%owner_mask,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vert owner_mask failed')
    ENDIF
    DEALLOCATE( p_patch%verts%glb_index,  &
      & p_patch%verts%loc_index,  &
      & p_patch%verts%owner_local,  &
      & p_patch%verts%owner_g,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vert data failed')
    ENDIF

    DEALLOCATE(p_patch%cells%phys_id,    &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell physical ID failed')
    ENDIF

    IF (l_ddmode) RETURN

    DEALLOCATE( p_patch%cells%edge_orientation,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell edge_orientation failed')
    ENDIF
    DEALLOCATE( p_patch%cells%area,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell area failed')
    ENDIF
    DEALLOCATE( p_patch%cells%f_c,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch cell f_c failed')
    ENDIF
    !
    !
    DEALLOCATE(p_patch%edges%phys_id,    &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge physical ID failed')
    ENDIF
    DEALLOCATE( p_patch%edges%cell_idx,  &
      & p_patch%edges%cell_blk,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge cell index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%vertex_idx,  &
      & p_patch%edges%vertex_blk,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge vertex index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%system_orientation,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge system orientation failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_idx,  &
      & p_patch%edges%quad_blk,  &
      & stat=ist)
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge quad index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_orientation,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge quad orientation failed')
    ENDIF
    DEALLOCATE( p_patch%edges%butterfly_idx,  &
      & p_patch%edges%butterfly_blk,  &
      & stat=ist)
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge butterfly index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%center,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge center failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge primal_normal failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge dual_normal failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal_cell,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge primal_normal_cell failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal_cell,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge dual_normal_cell failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal_vert,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge primal_normal_vert failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal_vert,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge dual_normal_vert failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_edge_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge primal edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_primal_edge_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for inverse patch edge primal edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_edge_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge dual edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_dual_edge_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for inverse patch edge dual edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%edge_vert_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for edge_vert_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_vert_vert_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for inverse vert_vert_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%edge_cell_length,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for edge_cell_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%area_edge,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge area_edge failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_area,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge quad area failed')
    ENDIF
    DEALLOCATE( p_patch%edges%f_e,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch edge f_e failed')
    ENDIF
    !
    !
    DEALLOCATE( p_patch%verts%neighbor_idx,  &
      & p_patch%verts%neighbor_blk,  &
      & p_patch%verts%phys_id,       &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex neighbor index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%cell_idx,  &
      & p_patch%verts%cell_blk,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex cell index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%edge_idx,  &
      & p_patch%verts%edge_blk,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex edge index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%edge_orientation,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex edge orientation failed')
    ENDIF
    DEALLOCATE( p_patch%verts%num_edges,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex number of edges failed')
    ENDIF
    DEALLOCATE( p_patch%verts%dual_area,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex dual area failed')
    ENDIF
    DEALLOCATE( p_patch%verts%f_v,  &
      & stat=ist )
    IF(ist/=success)THEN
      CALL finish  ('mo_model_domain_import:destruct_patches', &
        & 'deallocate for patch vertex f_v failed')
    ENDIF


    CALL deallocate_patch_cartestian( p_patch )
    
  END SUBROUTINE deallocate_patch
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE allocate_patch_cartestian( p_patch )
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch

    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_alloc_patches:allocate_patch_cartestian'

    ALLOCATE( p_patch%edges%primal_cart_normal(nproma,p_patch%nblks_e) )
    p_patch%edges%primal_cart_normal(:,:)%x(1) = 0._wp
    p_patch%edges%primal_cart_normal(:,:)%x(2) = 0._wp
    p_patch%edges%primal_cart_normal(:,:)%x(3) = 0._wp
    
    ALLOCATE( p_patch%edges%dual_cart_normal(nproma,p_patch%nblks_e) )
    p_patch%edges%dual_cart_normal(:,:)%x(1) = 0._wp
    p_patch%edges%dual_cart_normal(:,:)%x(2) = 0._wp
    p_patch%edges%dual_cart_normal(:,:)%x(3) = 0._wp
    
    ALLOCATE( p_patch%cells%cartesian_center(nproma,p_patch%nblks_c) )
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
  
  END SUBROUTINE allocate_patch_cartestian
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_patch_cartestian( p_patch )
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch

    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_alloc_patches:deallocate_patch_cartestian'

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

  END SUBROUTINE deallocate_patch_cartestian
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !> destruct_patches: Calls destruct_patch in a loop  
  SUBROUTINE destruct_patches( p_patch )
    
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch(n_dom_start:)
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_alloc_patches:destruct_patches'
    
    !local variables
    INTEGER :: jg
    !-----------------------------------------------------------------------
    
    CALL message (TRIM(routine), 'start')
    
    grid_level_loop: DO jg = n_dom_start, n_dom
      
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
    ! - nblks_c
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells
    ! - n_patch_edges
    ! - n_patch_verts
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom
    
    p_patch%cell_type = global_cell_type
    max_childdom = p_patch%max_childdom
    
    !
    ! !grid cells
    !
    ALLOCATE( p_patch%cells%num_edges(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%parent_idx(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%parent_blk(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%pc_idx(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%child_idx(nproma,p_patch%nblks_c,4) )
    ALLOCATE( p_patch%cells%child_blk(nproma,p_patch%nblks_c,4) )
    ALLOCATE( p_patch%cells%child_id(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%phys_id(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%neighbor_idx(nproma,p_patch%nblks_c,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%neighbor_blk(nproma,p_patch%nblks_c,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%edge_idx(nproma,p_patch%nblks_c,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%edge_blk(nproma,p_patch%nblks_c,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%vertex_idx(nproma,p_patch%nblks_c,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%vertex_blk(nproma,p_patch%nblks_c,p_patch%cell_type) )
    ALLOCATE( p_patch%cells%center(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%refin_ctrl(nproma,p_patch%nblks_c) )
    ALLOCATE( p_patch%cells%start_idx(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%end_idx(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%start_blk(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%end_blk(min_rlcell:max_rlcell,max_childdom) )
    
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
    
    !
    ! !grid verts
    !
    ALLOCATE( p_patch%verts%vertex(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%refin_ctrl(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%start_idx(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%end_idx(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%start_blk(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%end_blk(min_rlvert:max_rlvert,max_childdom) )
    
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
    
    p_patch%verts%vertex(:,:)%lon = 0._wp
    p_patch%verts%vertex(:,:)%lat = 0._wp
    p_patch%verts%refin_ctrl = 0
    p_patch%verts%start_idx = 0
    p_patch%verts%end_idx = 0
    p_patch%verts%start_blk = 0
    p_patch%verts%end_blk = 0
    
    
  END SUBROUTINE allocate_basic_patch
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
    !
    ! !grid verts
    !
    DEALLOCATE( p_patch%verts%vertex )
    DEALLOCATE( p_patch%verts%refin_ctrl )
    DEALLOCATE( p_patch%verts%start_idx )
    DEALLOCATE( p_patch%verts%end_idx )
    DEALLOCATE( p_patch%verts%start_blk )
    DEALLOCATE( p_patch%verts%end_blk )
    
  END SUBROUTINE deallocate_basic_patch
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
    
    INTEGER :: jc, je, jv
    
    ! Please note: The following variables in the patch MUST already be set:
    ! - nblks_c
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
      IF (.NOT. ALLOCATED(p_patch%cells%phys_id)) ALLOCATE( p_patch%cells%phys_id(nproma,p_patch%nblks_c) )
      ALLOCATE( p_patch%cells%edge_orientation(nproma,p_patch%nblks_c,p_patch%cell_type) )
      ALLOCATE( p_patch%cells%area(nproma,p_patch%nblks_c) )
      ALLOCATE( p_patch%cells%f_c(nproma,p_patch%nblks_c) )
    ENDIF
    
    IF (iopmode /= 3) THEN
      ALLOCATE( p_patch%cells%decomp_domain(nproma,p_patch%nblks_c) )
      ALLOCATE( p_patch%cells%owner_mask(nproma,p_patch%nblks_c) )
      ALLOCATE( p_patch%cells%glb_index(p_patch%n_patch_cells) )
      ALLOCATE( p_patch%cells%owner_local(p_patch%n_patch_cells))
      ALLOCATE( p_patch%cells%loc_index(p_patch%n_patch_cells_g) )
      ALLOCATE( p_patch%cells%owner_g(p_patch%n_patch_cells_g))

      IF (my_process_is_mpi_seq()) THEN
        p_patch%cells%decomp_domain(:,:) = 0
        p_patch%cells%owner_mask(:,:)    = .true.
        p_patch%cells%owner_local(:)     = 0
        p_patch%cells%owner_g(:)         = 0
      ENDIF
      
    ENDIF
    
    !
    ! !grid edges
    !
    IF (iopmode /= 2) THEN
      ALLOCATE( p_patch%edges%phys_id(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%cell_idx(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%cell_blk(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%vertex_idx(nproma,p_patch%nblks_e,4) )
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

    
    IF (iopmode /= 3) THEN
      ALLOCATE( p_patch%edges%decomp_domain(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%owner_mask(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%glb_index(p_patch%n_patch_edges) )
      ALLOCATE( p_patch%edges%loc_index(p_patch%n_patch_edges_g) )
      ALLOCATE( p_patch%edges%owner_g(p_patch%n_patch_edges_g))
      ALLOCATE( p_patch%edges%owner_local(p_patch%n_patch_edges))

      IF (my_process_is_mpi_seq()) THEN
        p_patch%edges%decomp_domain(:,:) = 0
        p_patch%edges%owner_mask(:,:)    = .true.
        p_patch%edges%owner_local(:)     = 0
        p_patch%edges%owner_g(:)         = 0
      ENDIF

    ENDIF
    
    !
    ! !grid verts
    !
    IF (iopmode /= 2) THEN
      ALLOCATE( p_patch%verts%phys_id(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%neighbor_idx(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%neighbor_blk(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%cell_idx(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%cell_blk(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%edge_idx(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%edge_blk(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%edge_orientation(nproma,p_patch%nblks_v,9-p_patch%cell_type) )
      ALLOCATE( p_patch%verts%num_edges(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%dual_area(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%f_v(nproma,p_patch%nblks_v) )
    ENDIF
    
    IF (iopmode /= 3) THEN
      ALLOCATE( p_patch%verts%decomp_domain(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%owner_mask(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%glb_index(p_patch%n_patch_verts) )
      ALLOCATE( p_patch%verts%loc_index(p_patch%n_patch_verts_g) )
      ALLOCATE( p_patch%verts%owner_g(p_patch%n_patch_verts_g))
      ALLOCATE( p_patch%verts%owner_local(p_patch%n_patch_verts))
      
      IF (my_process_is_mpi_seq()) THEN
        p_patch%verts%decomp_domain(:,:) = 0
        p_patch%verts%owner_mask(:,:)    = .true.
        p_patch%verts%owner_local(:)     = 0
        p_patch%verts%owner_g(:)         = 0
      ENDIF
      
    ENDIF
    
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

    ! Set values which are needed for parallel runs
    ! to the correct values for a single patch owner
    ! decomp_domain/owner mask:
    ! Everywhere 0 or .true. with the exception of unused entries
    
    IF (iopmode /= 3) THEN
      p_patch%cells%decomp_domain = 0
      p_patch%edges%decomp_domain = 0
      p_patch%verts%decomp_domain = 0
      p_patch%cells%owner_mask = .TRUE.
      p_patch%edges%owner_mask = .TRUE.
      p_patch%verts%owner_mask = .TRUE.
    
      p_patch%cells%decomp_domain(p_patch%npromz_c+1:nproma,p_patch%nblks_c) = -1
      p_patch%edges%decomp_domain(p_patch%npromz_e+1:nproma,p_patch%nblks_e) = -1
      p_patch%verts%decomp_domain(p_patch%npromz_v+1:nproma,p_patch%nblks_v) = -1
    
      p_patch%cells%owner_mask(p_patch%npromz_c+1:nproma,p_patch%nblks_c) = .FALSE.
      p_patch%edges%owner_mask(p_patch%npromz_e+1:nproma,p_patch%nblks_e) = .FALSE.
      p_patch%verts%owner_mask(p_patch%npromz_v+1:nproma,p_patch%nblks_v) = .FALSE.
    
      ! The following arrays are currently never needed for non parallel runs,
      ! we set them nontheless.
    
      DO jc = 1, p_patch%n_patch_cells
        p_patch%cells%glb_index(jc) = jc
        p_patch%cells%loc_index(jc) = jc
        p_patch%cells%owner_g(jc) = 0
        p_patch%cells%owner_local(jc) = 0
      ENDDO
    
      DO je = 1, p_patch%n_patch_edges
        p_patch%edges%glb_index(je) = je
        p_patch%edges%loc_index(je) = je
        p_patch%edges%owner_g(je) = 0
        p_patch%edges%owner_local(je) = 0
      ENDDO
    
      DO jv = 1, p_patch%n_patch_verts
        p_patch%verts%glb_index(jv) = jv
        p_patch%verts%loc_index(jv) = jv
        p_patch%verts%owner_g(jv) = 0
        p_patch%verts%owner_local(jv) = 0
      ENDDO
    ENDIF

    IF (iopmode /= 2) CALL allocate_patch_cartestian( p_patch )
    
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


