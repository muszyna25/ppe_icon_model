!>
!! @par Revision History
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_icon_testbed_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: filename_max
  USE mo_icon_testbed_config, ONLY: &
    & config_testbed_model        => testbed_model,       &
    & config_testbed_iterations   => testbed_iterations,  &
    & config_calculate_iterations => calculate_iterations,&
    & config_no_of_blocks         => no_of_blocks,        &
    & config_no_of_layers         => no_of_layers,        &
    & config_testfile_3D_time          => testfile_3D_time,         &
    & config_testfile_2D_time         => testfile_2D_time,        &
    & null_model
    
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_icon_testbed_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

       
  CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_icon_testbed_namelist( filename )

    ! ------------------------------------------------------------------------
    INTEGER :: testbed_model
    INTEGER :: testbed_iterations
    INTEGER :: calculate_iterations
    INTEGER  :: no_of_blocks, no_of_layers

   CHARACTER(LEN=filename_max) :: testfile_3D_time(2), testfile_2D_time(2)
    
    NAMELIST /testbed_nml/ testbed_model, testbed_iterations, calculate_iterations, &
      & no_of_blocks, no_of_layers, testfile_3D_time, testfile_2D_time

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat
    !0!CHARACTER(len=*), PARAMETER ::   &
    !0!        &  method_name = 'mo_parallel_nml:read_parallel_namelist'

    !--------------------------------------------
    ! set default values
    !--------------------------------------------
    testbed_model         = null_model
    testbed_iterations   = 10
    calculate_iterations = 10
    no_of_blocks         = 16
    no_of_layers         = 80
    testfile_3D_time     = ""
    testfile_2D_time     = ""

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes) 
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('testbed_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, testbed_nml)
    END SELECT
    CALL close_nml
    
    !-----------------------------------------------------
    ! fill_config_testbed
    config_testbed_model        = testbed_model
    config_testbed_iterations   = testbed_iterations
    config_calculate_iterations = calculate_iterations
    config_no_of_blocks         = no_of_blocks
    config_no_of_layers         = no_of_layers
    config_testfile_3D_time     = testfile_3D_time
    config_testfile_2D_time     = testfile_2D_time
    
  END SUBROUTINE read_icon_testbed_namelist
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed_nml
