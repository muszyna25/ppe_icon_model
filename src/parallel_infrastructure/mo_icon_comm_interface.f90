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
  USE mo_parallel_config, ONLY: nproma, icon_comm_debug
  USE mo_grid_config,     ONLY: n_dom
  USE mo_model_domain,    ONLY: p_patch
  USE mo_icoham_dyn_memory,ONLY: p_hydro_state
  USE mo_mpi,             ONLY: my_process_is_mpi_seq
  USE mo_icon_comm_lib

#ifdef _OPENMP
  USE omp_lib, ONLY: omp_in_parallel
#endif

  IMPLICIT NONE

  PRIVATE

  ! public constants
  PUBLIC :: construct_icon_communication
  PUBLIC :: destruct_icon_communication
  
  !-------------------------------------------------------------------------
CONTAINS

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE construct_icon_communication()

    INTEGER :: grid_id
    
    CHARACTER(*), PARAMETER :: method_name = "construct_icon_communication"


    CALL construct_icon_comm_lib()
    
    DO grid_id = 1, n_dom
      ! create the communication patterns
      CALL init_icon_std_comm_patterns(p_patch(grid_id))

      ! create the communicators for major variables
!       p_hydro_state(grid_id)%tend_phy%temp_comm = &
!         & new_comm_variable(p_hydro_state(grid_id)%tend_phy%temp, on_cells, &
!         & p_patch(grid_id))
    
    ENDDO  
    
    RETURN

  END SUBROUTINE construct_icon_communication
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE destruct_icon_communication()

     CALL destruct_icon_comm_lib()
     
  END SUBROUTINE destruct_icon_communication
  

END MODULE mo_icon_comm_interface
