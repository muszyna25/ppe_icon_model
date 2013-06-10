!>
!! Basic tools for the patch
!!
!! @par Revision History
!! Initial version  by: Leonidas Linardakis,  MPI-M, Hamburg, June 2013
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_grid_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, warning
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_math_utilities,     ONLY: gvec2cvec, gc2cc
  USE mo_math_types
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: calculate_patch_cartesian_positions
  PUBLIC :: rescale_grid
  PUBLIC :: calculate_edge_area

!  PUBLIC :: create_dummy_cell_closure

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  
  !-------------------------------------------------------------------------
  
CONTAINS

  
  !-------------------------------------------------------------------------
  !> Rescale grids
  ! Note: this does not rescale the cartesian coordinates (used in torus)
  SUBROUTINE rescale_grid( patch, grid_length_rescale_factor)
    
    TYPE(t_patch), INTENT(inout), TARGET ::  patch  ! patch data structure
    REAL(wp), INTENT(in) :: grid_length_rescale_factor
    
    REAL(wp) :: grid_area_rescale_factor
    !-----------------------------------------------------------------------
    grid_area_rescale_factor   = grid_length_rescale_factor * grid_length_rescale_factor
    
    patch%cells%area(:,:)               = &
      & patch%cells%area(:,:)               * grid_area_rescale_factor
    patch%verts%dual_area(:,:)          = &
      & patch%verts%dual_area(:,:)          * grid_area_rescale_factor
    patch%edges%primal_edge_length(:,:) = &
      & patch%edges%primal_edge_length(:,:) * grid_length_rescale_factor
    patch%edges%dual_edge_length(:,:)   = &
      & patch%edges%dual_edge_length(:,:)   * grid_length_rescale_factor
    patch%edges%edge_cell_length(:,:,:) = &
      & patch%edges%edge_cell_length(:,:,:) * grid_length_rescale_factor
    patch%edges%edge_vert_length(:,:,:) = &
      & patch%edges%edge_vert_length(:,:,:) * grid_length_rescale_factor

    ! rescale geometry parameters
    patch%geometry_info%mean_edge_length = &
      & patch%geometry_info%mean_edge_length * grid_length_rescale_factor
    patch%geometry_info%mean_cell_area   = &
      & patch%geometry_info%mean_cell_area   * grid_area_rescale_factor
    patch%geometry_info%domain_length    = &
      & patch%geometry_info%domain_length    * grid_length_rescale_factor
    patch%geometry_info%domain_height    = &
      & patch%geometry_info%domain_height    * grid_length_rescale_factor
    patch%geometry_info%sphere_radius    = &
      & patch%geometry_info%sphere_radius    * grid_length_rescale_factor
    patch%geometry_info%mean_characteristic_length    = &
      & patch%geometry_info%mean_characteristic_length * grid_length_rescale_factor

!     write(0,*) "Rescale grid_length_rescale_factor:", &
!       & grid_length_rescale_factor
!     write(0,*) "Rescale mean_cell_area:", &
!       & patch%geometry_info%mean_cell_area
!     write(0,*) "Rescale mean_characteristic_length:", &
!       & patch%geometry_info%mean_characteristic_length

    IF (patch%geometry_info%mean_characteristic_length == 0.0_wp) &
      & CALL finish("rescale_grid", "mean_characteristic_length=0")
    
  END SUBROUTINE rescale_grid
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! This routine calculates the edge area of the patch
  !!
  SUBROUTINE calculate_edge_area(patch)
    TYPE(t_patch), TARGET, INTENT(inout) :: patch

    INTEGER                       :: je, jb, istart_e, iend_e
    TYPE(t_subset_range), POINTER :: all_edges
    !-------------------------------------------------------------------------
    all_edges => patch%edges%all

    ! a) the control volume associated to each edge is defined as the
    ! quadrilateral whose edges are the primal edge and the associated dual edge
    !----------------------------------------------------------------------------
    ! loop over all blocks and edges

!$OMP PARALLEL DO PRIVATE(jb,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, istart_e, iend_e)
      DO je = istart_e, iend_e
        patch%edges%area_edge(je,jb) =  &
            &    patch%edges%primal_edge_length(je,jb)  &
            &  * patch%edges%dual_edge_length(je,jb)
      END DO
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE calculate_edge_area
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! This method_name calculates the cartesian positions stored in patch.
  !!
  !! These values normally are read from the grid files,
  !! this routine is used only for old grids
  !!
  !! Note: cartesian_dual_middle is used only from the ocean,
  !!   which is using only new grids, therefor it is NOT calculated here
  !!
  !! @par Revision History
  !! Initial release by Jochen Foerstner (2008-05-19)
  !! Modifiaction by A. Gassmann(2010-09-05)
  !! - added also tangential normal, and generalize to lplane
  !!
  SUBROUTINE calculate_patch_cartesian_positions ( patch )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: patch  ! patch on a specific level
    
    TYPE(t_cartesian_coordinates) :: z_vec
    REAL(wp) :: z_lon, z_lat, z_u, z_v  ! location and components of normal
    REAL(wp) :: z_norm                  ! norm of Cartesian normal
    
    INTEGER :: start_index, end_index
    INTEGER :: jb, je                   ! loop indices    
    !-----------------------------------------------------------------------
                
!$OMP PARALLEL    
    ! calculate cells cartesian positions
!$OMP DO PRIVATE(jb,je, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = patch%cells%all%start_block, patch%cells%all%end_block
      CALL get_index_range(patch%cells%all, jb, start_index, end_index)
      DO je = start_index, end_index            
        ! location of cell center
        patch%cells%cartesian_center(je,jb) = gc2cc(patch%cells%center(je,jb))
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
    
    ! calculate verts cartesian positions
!$OMP DO PRIVATE(jb,je, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = patch%verts%all%start_block, patch%verts%all%end_block
      CALL get_index_range(patch%verts%all, jb, start_index, end_index)
      DO je = start_index, end_index            
        ! location of cell center
        patch%verts%cartesian(je,jb) = gc2cc(patch%verts%vertex(je,jb))
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

    ! calculate edges positions
!$OMP DO PRIVATE(jb,je,z_lon,z_lat,z_u,z_v,z_norm,z_vec, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = patch%edges%all%start_block, patch%edges%all%end_block
      CALL get_index_range(patch%edges%all, jb, start_index, end_index)
      DO je = start_index, end_index
            
        ! location of edge midpoint
        patch%edges%cartesian_center(je,jb) = gc2cc(patch%edges%center(je,jb))
        z_lon = patch%edges%center(je,jb)%lon
        z_lat = patch%edges%center(je,jb)%lat
        
        ! zonal and meridional component of primal normal
        z_u = patch%edges%primal_normal(je,jb)%v1
        z_v = patch%edges%primal_normal(je,jb)%v2
        
        ! calculate Cartesian components of primal normal
        CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )
        
        ! compute unit normal to edge je
        z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
        z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)
        
        ! save the values in the according type structure of the patch
!         WRITE(*,*) "----------------------------------"
!         write(*,*) "primal_cart_normal:", patch%edges%primal_cart_normal(je,jb)%x(:)
!         write(*,*) "z_vec:", z_vec%x
!         WRITE(*,*) "----------------------------------"
!         IF ( MAXVAL(ABS(patch%edges%primal_cart_normal(je,jb)%x - z_vec%x)) > 0.0001_wp) &
!           CALL finish("","primal_cart_normal(je,jb)%x /=  z_vec%x")
               
        patch%edges%primal_cart_normal(je,jb)%x(1) = z_vec%x(1)
        patch%edges%primal_cart_normal(je,jb)%x(2) = z_vec%x(2)
        patch%edges%primal_cart_normal(je,jb)%x(3) = z_vec%x(3)
        
        ! zonal and meridional component of dual normal
        z_u = patch%edges%dual_normal(je,jb)%v1
        z_v = patch%edges%dual_normal(je,jb)%v2
        
        ! calculate Cartesian components of dual normal
        CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )
        
        ! compute unit normal to edge je
        z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
        z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)
        
!         WRITE(*,*) "----------------------------------"
!         write(*,*) "dual_cart_normal:", patch%edges%dual_cart_normal(je,jb)%x(:)
!         write(*,*) "z_vec:", z_vec%x
!         WRITE(*,*) "----------------------------------"
!         IF ( MAXVAL(ABS(patch%edges%dual_cart_normal(je,jb)%x - z_vec%x)) > 0.0001_wp) &
!           CALL finish("","dual_cart_normal(je,jb)%x /=  z_vec%x")
        
        ! save the values in the according type structure of the patch
        patch%edges%dual_cart_normal(je,jb)%x(1) = z_vec%x(1)
        patch%edges%dual_cart_normal(je,jb)%x(2) = z_vec%x(2)
        patch%edges%dual_cart_normal(je,jb)%x(3) = z_vec%x(3)
        
      END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
  END SUBROUTINE calculate_patch_cartesian_positions
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  !>
  ! adds a dummy cell at the end of the existing cells
  ! and creates a "dummy" connectivity for edges and cells that have no neighbor
  ! Note that the dummy cell has to be already allocated!
  ! This is indented for the ocean, with no land cells
  !
  ! The subsets must have been filled in order in order to call this routine.
!  SUBROUTINE create_dummy_cell_closure( patch )
!    TYPE(t_patch), INTENT(inout), TARGET ::  patch
!
!    INTEGER :: block, idx, start_idx, end_idx, neighbor
!    INTEGER :: dummy_cell_blk, dummy_cell_idx
!
!    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grid_tools:create_dummy_cell_closure'
!
!    ! the last cell is the dummy cell
!    dummy_cell_blk = patch%nblks_c
!    dummy_cell_idx = patch%npromz_c
!
!
!     write(0,*) "------------ start  create_dummy_cell_closure ---------------"
!
! !$OMP PARALLEL
! !$OMP DO PRIVATE(block, idx, start_idx, end_idx) ICON_OMP_DEFAULT_SCHEDULE
!    DO block = patch%edges%all%start_block, patch%edges%all%end_block
!      CALL get_index_range(patch%edges%all, block, start_idx, end_idx)
!      DO idx = start_idx, end_idx
!        DO neighbor=1,2
!          IF ( patch%edges%cell_idx(idx, block, neighbor) == 0) THEN
!            patch%edges%cell_blk(idx, block, neighbor) = dummy_cell_blk
!            patch%edges%cell_idx(idx, block, neighbor) = dummy_cell_idx
!          ENDIF
!        END DO
!      END DO
!    END DO
! !$OMP END DO
!
! !$OMP PARALLEL DO PRIVATE(block, idx, start_idx, end_idx) ICON_OMP_DEFAULT_SCHEDULE
!    DO block = patch%cells%all%start_block, patch%cells%all%end_block
!      CALL get_index_range(patch%cells%all, block, start_idx, end_idx)
!      DO idx = start_idx, end_idx
!         IF (block /= dummy_cell_blk .OR. idx /= dummy_cell_idx) THEN
!          ! this is not the dummy cell
!          DO neighbor=1, patch%cells%max_connectivity
!            IF ( patch%cells%neighbor_idx(idx, block, neighbor) == 0) THEN
!              patch%cells%neighbor_blk(idx, block, neighbor) = dummy_cell_blk
!              patch%cells%neighbor_idx(idx, block, neighbor) = dummy_cell_idx
!  !             write(0,*) "Replaced neighbor at ", block, idx, neighbor
!            ENDIF
!          END DO
!         ENDIF
!      END DO
!    END DO
! !$OMP END DO
! !$OMP END PARALLEL
!
!  END SUBROUTINE create_dummy_cell_closure
  !-----------------------------------------------------------------------
    
  
END MODULE mo_grid_tools
