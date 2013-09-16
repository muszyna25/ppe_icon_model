!>
!! Contains the implementation of the top and bottom ocean boundary conditions
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2010-04)
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!  methods used are mpi parallelized, LL
!!
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_oce_check_consistency
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: rho_ref
  USE mo_impl_constants,     ONLY: max_char_length, sea_boundary, boundary, sea, min_dolic
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce, i_bc_veloc_top, i_bc_veloc_bot
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_run_config,         ONLY: dtime
  USE mo_exception,          ONLY: message, finish
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  
  PUBLIC :: ocean_check_level_sea_land_mask

CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_check_level_sea_land_mask( patch_3d )
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3d

    INTEGER :: block, idx, start_idx, end_idx, level
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
    TYPE(t_subset_range), POINTER :: all_cells, owned_edges
    TYPE(t_patch), POINTER        :: patch_2d
    CHARACTER(len=*), PARAMETER :: method_name='mo_oce_check_consistency:ocean_check_level_sea_land_mask'
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2D(1)
    all_cells => patch_2d%cells%all
    owned_edges => patch_2d%edges%owned
    !-----------------------------------------------------------------------
    ! check cells
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_idx, end_idx)
      DO idx = start_idx, end_idx

        DO level=1, patch_3d%p_patch_1d(1)%dolic_c(idx, block)
          IF(patch_3d%lsm_c(idx, level, block) > sea_boundary) THEN
            write(0,*) " dolic_c=", patch_3d%p_patch_1d(1)%dolic_c(idx, block), " at level=", level, &
              & "lsm_c= ", patch_3d%lsm_c(idx, level, block)
            CALL finish(method_name,"Inconsistent cell levels vs sea_land_mask")
          ENDIF
        ENDDO

        DO level=patch_3d%p_patch_1d(1)%dolic_c(idx, block)+1, n_zlev
          IF(patch_3d%lsm_c(idx, level, block) < boundary)THEN
            write(0,*) " dolic_c=", patch_3d%p_patch_1d(1)%dolic_c(idx, block), " at level=", level, &
              & "lsm_c= ", patch_3d%lsm_c(idx, level, block)
            CALL finish(method_name,"Inconsistent cell levels vs sea_land_mask")
          ENDIF
        ENDDO

      ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    ! check edges
    DO block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, block, start_idx, end_idx)
      DO idx = start_idx, end_idx

        cell1_idx = patch_2d%edges%cell_idx(idx, block, 1)
        cell1_blk = patch_2d%edges%cell_blk(idx, block, 1)
        cell2_idx = patch_2d%edges%cell_idx(idx, block, 2)
        cell2_blk = patch_2d%edges%cell_blk(idx, block, 2)

        ! check sea
        DO level=1, patch_3d%p_patch_1d(1)%dolic_e(idx, block)
          IF(patch_3d%lsm_e(idx, level, block) > sea_boundary) THEN
            write(0,*) " dolic_e=", patch_3d%p_patch_1d(1)%dolic_e(idx, block), " at level=", level, &
              & "lsm_e= ", patch_3d%lsm_e(idx, level, block)
            CALL finish(method_name,"Inconsistent edge levels vs sea_land_mask")
          ENDIF

          ! check neigboring cells
          IF (cell1_idx < 1 .or. cell2_idx < 1) &
            & CALL finish(method_name,"Sea edge with one cell missing")
          IF (patch_3d%lsm_c(cell1_idx, level, cell1_blk) >= boundary .or. &
              & patch_3d%lsm_c(cell2_idx, level, cell2_blk) >= boundary) THEN
            write(0,*) " at level=", level
            CALL finish(method_name,"Sea edge with one cell land")
          ENDIF

        ENDDO

        ! check land
        DO level=patch_3d%p_patch_1d(1)%dolic_e(idx, block)+1, n_zlev
          IF(patch_3d%lsm_e(idx, level, block) < boundary)THEN
            write(0,*) " dolic_e=", patch_3d%p_patch_1d(1)%dolic_e(idx, block), " at level=", level, &
              & "lsm_e= ", patch_3d%lsm_e(idx, level, block)
            CALL finish(method_name,"Inconsistent edge levels vs sea_land_mask")
          ENDIF

          ! check neigboring cells
          IF (cell1_idx > 0 .and. cell2_idx > 0) THEN
            IF (patch_3d%lsm_c(cell1_idx, level, cell1_blk) < boundary .and. &
                & patch_3d%lsm_c(cell2_idx, level, cell2_blk) < boundary) THEN
              write(0,*) " at level=", level
              CALL finish(method_name,"Land edge with two cells sea")
            ENDIF
          ENDIF

        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE ocean_check_level_sea_land_mask
  !-------------------------------------------------------------------------
    
END MODULE mo_oce_check_consistency
