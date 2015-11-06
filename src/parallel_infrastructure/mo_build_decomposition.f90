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

  USE mo_complete_subdivision, ONLY: finalize_decomposition, &
    &                                copy_processor_splitting, &
    &                                complete_parallel_setup
  USE mo_setup_subdivision,    ONLY: decompose_domain
  USE mo_sync,                 ONLY: disable_sync_checks, enable_sync_checks,                    &
    &                                decomposition_statistics
  USE mo_grid_config,          ONLY: n_dom, n_dom_start
  USE mo_mpi,                  ONLY: my_process_is_mpi_parallel
  USE mo_loopindices,          ONLY: get_indices_e
  USE mo_model_domain,         ONLY: p_patch, t_pre_patch, t_patch_3d, &
    &                                t_patch, p_patch_local_parent
  USE mo_model_domimp_patches, ONLY: reorder_patch_refin_ctrl, &
    & import_pre_patches, complete_patches
  USE mo_parallel_config,      ONLY: p_test_run, l_test_openmp, num_io_procs, division_method
  USE mo_impl_constants,       ONLY: success, max_dom
  USE mo_exception,            ONLY: finish, message, message_text, get_filename_noext

  IMPLICIT NONE
  
  PUBLIC :: build_decomposition
  
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
    INTEGER                        :: error_status,jg
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

  
END MODULE mo_build_decomposition
