!>
!!
!! @par Revision History
!!!
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
MODULE mo_build_decomposition

  USE mo_complete_subdivision
  USE mo_setup_subdivision
  USE mo_sync,                ONLY: disable_sync_checks, enable_sync_checks, decomposition_statistics
  USE mo_grid_config,         ONLY: n_dom, n_dom_start
  USE mo_mpi
  USE mo_kind
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants
  USE mo_model_domain,        ONLY: p_patch,t_patch_3d,t_patch,p_patch_local_parent
  USE mo_model_domimp_patches
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, division_method
  USE mo_impl_constants,      ONLY: success, max_dom
  USE mo_exception,           ONLY: finish, message, message_text, get_filename_noext
  USE mo_run_config,          ONLY: ldump_dd, lread_dd, nproc_dd
  
#ifndef __NO_ICON_ATMO__
  USE mo_dump_restore,        ONLY: restore_patches_netcdf, dump_domain_decomposition
#endif
  
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
  SUBROUTINE build_decomposition(num_lev,num_levp1,nshift,&
    & is_ocean_decomposition,l_restore_states, patch_3d)
    
    INTEGER, INTENT(in) :: num_lev(max_dom),     &
      & num_levp1(max_dom),nshift(max_dom)
    LOGICAL, INTENT(in) :: is_ocean_decomposition,           &
      & l_restore_states
    TYPE(t_patch_3d), POINTER, OPTIONAL :: patch_3d
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'build_decomposition'
    TYPE(t_patch), ALLOCATABLE :: p_patch_global(:)
    INTEGER :: error_status,jg
    LOGICAL :: lrestore
    
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
    
    ! consistency checks:
    
    ! ldump_dd makes not much sense in a test run, we simply stop
    IF (my_process_is_mpi_test() .AND. .NOT. lread_dd .AND. ldump_dd) THEN
      CALL finish(routine, "Flag <dump_dd> was set in test run!")
    ENDIF
    ! Flag <lread_dd> does not make sense in single processor run, we simply stop
    IF (.NOT. l_restore_states .AND. lread_dd .AND. .NOT. my_process_is_mpi_parallel() &
      & .AND. .NOT. my_process_is_mpi_test()) THEN
      CALL finish(routine, "Flag <lread_dd> was set in single processor run!")
    END IF
    
    lrestore = .FALSE.
    
#ifndef __NO_ICON_ATMO__
    lrestore = l_restore_states
    
    IF (lrestore .AND. .NOT. my_process_is_mpi_test()) THEN
      
      ! --------------------------------------------------
      ! CASE 1: Work PEs restore from dump state
      ! --------------------------------------------------
      
      ! Before the restore set p_comm_input_bcast to null
      CALL set_comm_input_bcast(null_comm_type)
      CALL restore_patches_netcdf( p_patch, .TRUE. )
      CALL set_patch_communicators(p_patch)
      ! After the restore is done set p_comm_input_bcast in the
      ! same way as it would be set in parallel_nml_setup when
      ! no restore is wanted:
      CALL set_comm_input_bcast()
      
    ELSE
      
      ! --------------------------------------------------
      ! CASE 2: Work PEs subdivide patches w/out restore
      ! --------------------------------------------------
      
      ! Test run: If lread_dd is not set, import_basic_patches()
      ! goes parallel with the call above, actually the Test PE
      ! reads and broadcasts to the others.  IF lread_dd is set,
      ! the read is in standalone mode.
      IF (my_process_is_mpi_test() .AND. (lread_dd .OR. lrestore)) &
        & CALL set_comm_input_bcast(null_comm_type)
      
      IF (lread_dd .AND. .NOT. my_process_is_mpi_test()) THEN
        
        ! ------------------------------------------------
        ! CASE 2a: Read domain decomposition from file
        ! ------------------------------------------------
        
        CALL restore_patches_netcdf( p_patch, .FALSE. )
        
      ELSE
        
        ! ------------------------------------------------
        ! CASE 2b: compute domain decomposition on-the-fly
        ! ------------------------------------------------
#endif
        
        ALLOCATE(p_patch_global(n_dom_start:n_dom))
        
        CALL import_basic_patches(p_patch_global,num_lev,num_levp1,nshift)
        
        ! use internal domain decomposition algorithm
        
        IF (ldump_dd .AND. .NOT. my_process_is_mpi_parallel()) THEN
          ! If ldump_dd is set in a single processor run, a domain
          ! decomposition for nproc_dd processors is done
          CALL decompose_domain(p_patch, p_patch_global, nproc_dd)
        ELSE
          CALL decompose_domain(p_patch, p_patch_global)
        END IF
        
        DEALLOCATE(p_patch_global)
        
#ifndef __NO_ICON_ATMO__
      ENDIF ! lread_dd
      
      
      IF(ldump_dd) THEN
        IF (my_process_is_mpi_parallel()) THEN
          CALL dump_domain_decomposition(p_patch(n_dom_start))
          DO jg = n_dom_start+1, n_dom
            CALL dump_domain_decomposition(p_patch(jg),p_patch_local_parent(jg))
          ENDDO
        END IF
        ! note: in single processor runs, the dump is done within
        ! decompose_domain
        CALL p_stop
        STOP
      ENDIF
#endif
      
      ! setup communication patterns (also done in sequential runs)
      IF (is_ocean_decomposition) THEN
        CALL complete_parallel_setup_oce(p_patch)
      ELSE
        CALL complete_parallel_setup
      END IF
      
      ! Complete information which is not yet read or calculated
      IF (lrestore .AND. my_process_is_mpi_test()) CALL disable_sync_checks
      CALL complete_patches( p_patch )
      IF (lrestore .AND. my_process_is_mpi_test()) CALL enable_sync_checks
      
      IF (my_process_is_mpi_test() .AND. (lread_dd .OR. lrestore)) &
        & CALL set_comm_input_bcast()
      
#ifndef __NO_ICON_ATMO__
    ENDIF
#endif
    
    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    !--------------------------------------------------------------------------------
    
    IF (is_ocean_decomposition) THEN
      CALL finalize_decomposition_oce(p_patch)
    ELSE
      CALL finalize_decomposition
    END IF
    
    IF(.NOT.p_test_run .AND. my_process_is_mpi_parallel()) THEN ! the call below hangs in test mode
      ! Print diagnostic information about domain decomposition
      DO jg = 1, n_dom
        CALL decomposition_statistics(p_patch(jg))
      ENDDO
    ENDIF
    
    ! reorder patch_local_parents according to their refin_ctrl flags
    DO jg = n_dom_start+1, n_dom
      CALL reorder_patch_refin_ctrl(p_patch_local_parent(jg), p_patch(jg))
    ENDDO
    
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
