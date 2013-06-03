!>
!! A collection of MPI communication tools
!!
!! @par Revision History
!! First version by Leonidas Linardakis,  MPI-M, November 2011.
!!
!! @par
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_icon_comm_interface

  USE mo_kind,            ONLY: wp
  USE mo_io_units,        ONLY: filename_max
  USE mo_exception,       ONLY: message_text, message, finish
  USE mo_parallel_config, ONLY: nproma, icon_comm_debug, p_test_run, use_icon_comm
  USE mo_grid_config,     ONLY: n_dom
  USE mo_model_domain,    ONLY: t_patch
!  USE mo_icoham_dyn_memory,ONLY: p_hydro_state
  USE mo_mpi,             ONLY: my_process_is_mpi_seq, my_process_is_mpi_parallel, &
    & work_mpi_barrier, my_process_is_mpi_test, p_barrier,        &
    & p_comm_work_test, p_comm_work
  USE mo_icon_comm_lib

#ifdef _OPENMP
  USE omp_lib, ONLY: omp_in_parallel
#endif

  IMPLICIT NONE

  PRIVATE

  ! public constants
  PUBLIC :: construct_icon_communication
  PUBLIC :: destruct_icon_communication
  
  PUBLIC :: icon_comm_barrier



  !-------------------------------------------------------------------------
CONTAINS

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE construct_icon_communication(patch, n_dom)
    TYPE(t_patch), TARGET ::  patch(:)
    INTEGER :: n_dom

    INTEGER :: grid_id
    
    CHARACTER(*), PARAMETER :: method_name = "construct_icon_communication"

    !-------------------------------------------------------------------------------------
    ! constrcut patch communicators
    ! this is not the right place, this should be aware of the local decomposition features
    ! as well the coupled setup, parallel io, etc.
    ! Will be moved in the future
    patch(:)%compute_is_parallel = my_process_is_mpi_parallel()
    patch(:)%is_in_parallel_test = p_test_run
    patch(:)%is_test_parallel_process = my_process_is_mpi_test()

    patch(:)%work_communicator = p_comm_work
    patch(:)%parallel_test_communicator = p_comm_work_test
    ! WRITE(0,*) "work_communicator, test=", p_comm_work, p_comm_work_test
    !-------------------------------------------------------------------------------------
    IF (.NOT. use_icon_comm) RETURN

    CALL construct_icon_comm_lib()
    
    DO grid_id = 1, n_dom
      ! create the communication patterns
      CALL init_icon_std_comm_patterns(patch(grid_id))

      ! create the communicators for major variables
!       p_hydro_state(grid_id)%tend_phy%temp_comm = &
!         & new_comm_variable(p_hydro_state(grid_id)%tend_phy%temp, on_cells, &
!         & patch(grid_id))
    
    ENDDO

!     CALL work_mpi_barrier()
!     CALL finish("barrier returns", "")
    
    RETURN

  END SUBROUTINE construct_icon_communication
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE destruct_icon_communication()

    IF (.NOT. use_icon_comm) RETURN

    CALL destruct_icon_comm_lib()
     
  END SUBROUTINE destruct_icon_communication
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  SUBROUTINE icon_comm_barrier(for_patch)
    TYPE(t_patch), TARGET ::  for_patch

    ! write(0,*) "icon_comm_barrier:", for_patch%parallel_test_communicator
#ifndef NOMPI
    IF (for_patch%is_in_parallel_test) THEN
      CALL p_barrier(for_patch%parallel_test_communicator)
    ELSEIF (for_patch%compute_is_parallel) THEN
      CALL p_barrier(for_patch%work_communicator)
    ENDIF
#else
    RETURN
#endif
  END SUBROUTINE icon_comm_barrier
  !-----------------------------------------------------------------------

END MODULE mo_icon_comm_interface
