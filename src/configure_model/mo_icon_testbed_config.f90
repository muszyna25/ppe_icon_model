!>
!!        
!! @par Revision History
!!   Created by Leonidas Linardakis, MPI-M, 2011-05-07
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
MODULE mo_icon_testbed_config

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max

  IMPLICIT NONE

  PRIVATE
  
  ! Exported parameters
  PUBLIC :: null_model, test_coupler_model, test_jitter_model, test_halo_communication, &
    & test_radiation_communication, test_netcdf_read_model
    
  ! Exported variables
  PUBLIC :: testbed_model
  PUBLIC :: testbed_iterations
  PUBLIC :: calculate_iterations
  PUBLIC :: no_of_blocks, no_of_layers
  PUBLIC :: o3_filename
  
  ! ---------------
  ! the tesbed modes
  INTEGER, PARAMETER :: null_model          = 0  ! does nothing
  INTEGER, PARAMETER :: test_coupler_model  = 1  ! test the coupler
  INTEGER, PARAMETER :: test_jitter_model   = 3  ! test the jitter
  INTEGER, PARAMETER :: test_halo_communication  = 4  ! test the mpi communication
  INTEGER, PARAMETER :: test_radiation_communication  = 5  ! test the mpi communication
  INTEGER, PARAMETER :: test_netcdf_read_model = 6
  

  INTEGER  :: testbed_model
  ! ---------------
  INTEGER  :: testbed_iterations
  INTEGER  :: calculate_iterations
  INTEGER  :: no_of_blocks, no_of_layers

  CHARACTER(LEN=filename_max) :: o3_filename


CONTAINS
  
  !-------------------------------------------------------------------------
  SUBROUTINE check_icon_testbed_config()
 
  END SUBROUTINE check_icon_testbed_config
  !-------------------------------------------------------------------------
    
END MODULE mo_icon_testbed_config
