!-------------------------------------------------------------------------------------
! mo_local_patch_hierarchy
!>
!! Provides the methods for creating the local patch hierarchy
!!
!! $Id: n/a$
!!
!! @par Revision History
!! Current implementation by
!!   Leonidas Linardakis, MPI-M, 2010-01-12
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_local_patch_hierarchy
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  ! USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_io_units,           ONLY: nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_local_grid
  USE mo_io_local_grid,      ONLY: write_netcdf_grid, read_no_of_subgrids
  USE mo_local_grid_hierarchy, ONLY: create_grid_hierarchy, set_bdy_indexing_depth

  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: create_patches
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  TYPE t_patch_node_type
    ! grid info
    INTEGER :: grid_id
    CHARACTER(LEN=filename_max) ::  in_grid_file_name
    CHARACTER(LEN=filename_max) ::  out_grid_file_name
    ! patch hierarchy info
    INTEGER :: patch_id
    INTEGER :: hierarchy_node_type  ! 1= root,  2=inter, 3=leaf
    INTEGER :: tree_level           ! (ie, the depth)
    INTEGER :: no_of_children
    INTEGER :: parent_patch_id
    INTEGER, POINTER :: child_patch_ids(:)
  END TYPE t_patch_node_type
  !-------------------------------------------------------------------------

  INTEGER :: no_of_patches, root_patch_id, root_grid_id
  TYPE(t_patch_node_type), ALLOCATABLE, TARGET :: patch_node(:)
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !   SUBROUTINE create_patches(param_file_name)
  !>
  !! The driving method for creating the patch hierarchy.
  !! Reads the namelists for creating a local patch hierarchy
  !! and call the create_patch_hierarchy
  !! Private
  SUBROUTINE create_patches(param_file_name)
    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER, PARAMETER ::  max_no_of_patches = max_no_of_grid_objects
    CHARACTER(LEN=filename_max) :: in_grid_file_names(max_no_of_patches)
    CHARACTER(LEN=filename_max) :: out_grid_file_names(max_no_of_patches)
    INTEGER :: parent_ids(max_no_of_patches), grid_ids(max_no_of_patches)
    INTEGER :: bdy_indexing_depth
    INTEGER :: i_status
    INTEGER :: total_no_of_patches

    NAMELIST /local_patch_creation/ total_no_of_patches, &
      & in_grid_file_names, out_grid_file_names, &
      & grid_ids, parent_ids, bdy_indexing_depth

    ! set default values
    bdy_indexing_depth = 0
    total_no_of_patches = 0
    
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('local_patch_creation',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,local_patch_creation)
    ELSE
      WRITE(message_text,'(a, a, a)') " File ", TRIM(param_file_name), " not POSITIONED"
      CALL finish('create_patches', message_text)
    ENDIF
    CLOSE(nnml)

    !-------------------------------------------------------------------------
    CALL create_patch_hierarchy(total_no_of_patches, &
      & in_grid_file_names, out_grid_file_names, &
      & grid_ids, parent_ids, bdy_indexing_depth)
    !-------------------------------------------------------------------------

  END SUBROUTINE create_patches
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE init_patch_tree(param_file_name)
  !>
  !! Reads the parameters for creating a local patch hierarchy and sets-up the patch tree.
  !! Private
  SUBROUTINE init_patch_tree(in_grid_file_names, out_grid_file_names, &
        grid_ids, parent_ids)

    CHARACTER(LEN=filename_max), INTENT(in) :: in_grid_file_names(:), &
      & out_grid_file_names(:)
    INTEGER, INTENT(in) :: parent_ids(:), grid_ids(:)

    INTEGER :: parent_id, no_of_children
    INTEGER :: i, i_status

    ! init patches
    ALLOCATE(patch_node(no_of_patches),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('initTree', 'ALLOCATE(patch(noOfPatches))')
    root_patch_id = 0

    DO i = 1, no_of_patches
      patch_node(i)%patch_id = 0
      patch_node(i)%hierarchy_node_type = 0
      patch_node(i)%no_of_children = 0
      patch_node(i)%parent_patch_id = 0
    ENDDO

    ! fill patch info
    DO i = 1, no_of_patches
      ! grid ids should have strictly ascending order
      IF (grid_ids(i) /= i) &
        & CALL finish ('initTree', 'gridIDs(i) /= i')
      IF (parent_ids(i) < 0 .or. parent_ids(i) > no_of_patches) &
        & CALL finish ('initTree', 'parentIDs out of range')

      patch_node(i)%patch_id = i
      patch_node(i)%in_grid_file_name = in_grid_file_names(i)
      patch_node(i)%out_grid_file_name = out_grid_file_names(i)
      patch_node(i)%parent_patch_id = parent_ids(i)
      IF (parent_ids(i) /= 0) &
        & patch_node(parent_ids(i))%no_of_children = patch_node(parent_ids(i))%no_of_children + 1
    ENDDO

    ! fill patch children links
    ! allocate child_ids
    DO i = 1, no_of_patches
      no_of_children = patch_node(i)%no_of_children
      IF (no_of_children > 0) THEN
        ALLOCATE(patch_node(i)%child_patch_ids(no_of_children),stat=i_status)
        IF (i_status > 0) &
          & CALL finish ('initTree', 'ALLOCATE(patch(noOfPatches))')
      ENDIF
      patch_node(i)%no_of_children = 0
    ENDDO
    ! fill child_ids
    DO i = 1, no_of_patches
      parent_id = patch_node(i)%parent_patch_id
      IF (parent_id > 0) THEN
        no_of_children = patch_node(parent_id)%no_of_children + 1
        patch_node(parent_id)%no_of_children = no_of_children
        patch_node(parent_id)%child_patch_ids(no_of_children) = i
      ENDIF
    ENDDO
    ! fill root/leaf properties
    DO i = 1, no_of_patches
      IF (patch_node(i)%parent_patch_id == 0) THEN
        root_patch_id = i
        patch_node(i)%hierarchy_node_type = root_node
      ELSEIF (patch_node(i)%no_of_children == 0) THEN
        patch_node(i)%hierarchy_node_type = leaf_node
      ELSE
        patch_node(i)%hierarchy_node_type = inner_node
      ENDIF
    ENDDO
    IF (root_patch_id == 0) &
      & CALL finish ('initTree', 'no root is defined')

    CALL fill_tree_levels(root_patch_id, 1)

    ! print info
    DO i = 1, no_of_patches
      WRITE(message_text,'(a)') "===================================="
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a, i3)') "patch ID=", patch_node(i)%patch_id
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a, a)') "   inGridFileName=",   TRIM(patch_node(i)%in_grid_file_name)
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a, a)') "   outGridFileName=",  TRIM(patch_node(i)%out_grid_file_name)
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a, i3)') "   nodeType=",        patch_node(i)%hierarchy_node_type
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a, i3)') "   parentPatchID=",   patch_node(i)%parent_patch_id
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a, i3)') "   numOfChildren=",   patch_node(i)%no_of_children
      CALL message ('', TRIM(message_text))
      IF (patch_node(i)%no_of_children > 0) THEN
        WRITE(message_text,'(a, i3)') "   childIDs=",      patch_node(i)%child_patch_ids
        CALL message ('', TRIM(message_text))
      ENDIF
      WRITE(message_text,'(a)') "===================================="
      CALL message ('', TRIM(message_text))
    ENDDO

  END SUBROUTINE init_patch_tree
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   RECURSIVE SUBROUTINE fill_tree_levels(patch_id, level)
  !>
  !! Computes the depth of each node in the patch tree.
  !! Private
  RECURSIVE SUBROUTINE fill_tree_levels(patch_id, level)
    INTEGER, INTENT(in) :: patch_id, level

    INTEGER :: i,next_level

    patch_node(patch_id)%tree_level = level
    next_level = level + 1
    DO i=1, patch_node(patch_id)%no_of_children
      CALL fill_tree_levels(patch_node(patch_id)%child_patch_ids(i), next_level)
    ENDDO

  END SUBROUTINE fill_tree_levels
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE init_patch_grids()
  !>
  !! Fills the hierarchy values in the grid objects.
  !! Private
  SUBROUTINE init_patch_grids()

    TYPE(t_grid),  POINTER :: grid_obj
    TYPE(t_patch_node_type), POINTER :: patch_obj

    INTEGER :: patch_id, grid_id, no_of_children
    INTEGER :: i, start_subgrid_no, no_of_subgrids

    DO patch_id = 1, no_of_patches
      patch_node(patch_id)%grid_id = new_grid()
    ENDDO

    start_subgrid_no = 1
    DO patch_id = 1, no_of_patches
      patch_obj => patch_node(patch_id)
      grid_id = patch_obj%grid_id
      no_of_children =  patch_obj%no_of_children

      grid_obj => get_grid(grid_id)
      grid_obj%file_name           = patch_obj%in_grid_file_name
      grid_obj%out_file_name       = patch_obj%out_grid_file_name
      grid_obj%patch_id            = patch_id
      grid_obj%hierarchy_node_type = patch_obj%hierarchy_node_type

      grid_obj%start_subgrid_id    = start_subgrid_no
      no_of_subgrids = read_no_of_subgrids(grid_obj%file_name)
      ! print *, 'no_of_subgrids:',no_of_subgrids
      ! print *, 'start_subgrid_no:',start_subgrid_no
      start_subgrid_no = start_subgrid_no + no_of_subgrids

      ! grid_obj%level               = patch_obj%tree_level
      IF (patch_obj%parent_patch_id > 0) THEN
        grid_obj%parent_grid_id      = patch_node(patch_obj%parent_patch_id)%grid_id
      ELSE
        grid_obj%parent_grid_id      = 0
      ENDIF

      grid_obj%no_of_children      = no_of_children
      ALLOCATE(grid_obj%child_grid_ids(no_of_children),stat=i)
      IF (i > 0) THEN
        WRITE (message_text, '(a,i8,a)') &
          & 'ALLOCATE grid_obj%child_grid_ids with ', &
          & grid_obj%no_of_children, '  children.'
        CALL finish ('construct_cells', TRIM(message_text))
      ENDIF
      ! fill the children grid id
      DO i=1,no_of_children
        grid_obj%child_grid_ids(i) = patch_node( patch_obj%child_patch_ids(i) )%grid_id
      ENDDO

    ENDDO ! patch_id = 1, no_of_patches

    root_grid_id = patch_node(root_patch_id)%grid_id

  END SUBROUTINE init_patch_grids
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  !   SUBROUTINE create_patch_hierarchy()
  !>
  !! The main method for creating the patch hierarchy.
  !! Private
  SUBROUTINE create_patch_hierarchy(in_no_of_patches, in_grid_file_names, out_grid_file_names, &
        grid_ids, parent_ids, bdy_indexing_depth)

    INTEGER, INTENT(in) :: in_no_of_patches 
    CHARACTER(LEN=filename_max), INTENT(in) :: in_grid_file_names(:), &
      & out_grid_file_names(:)
    INTEGER, INTENT(in) :: parent_ids(:), grid_ids(:), bdy_indexing_depth
    
    no_of_patches = in_no_of_patches
    
    CALL init_patch_tree(in_grid_file_names, out_grid_file_names, &
        grid_ids, parent_ids)

    CALL init_patch_grids()

    CALL set_bdy_indexing_depth(bdy_indexing_depth)
    ! process recursively the patch tree
    CALL create_grid_hierarchy(root_grid_id)
    ! CALL process_patch(root_patch_id)

    ! write and clean the root grid
    CALL write_netcdf_grid(root_grid_id)
    CALL delete_grid(root_grid_id)

    RETURN

  END SUBROUTINE create_patch_hierarchy
  !-------------------------------------------------------------------------


END MODULE mo_local_patch_hierarchy

