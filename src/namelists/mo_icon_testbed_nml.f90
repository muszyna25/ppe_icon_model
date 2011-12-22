!>
!!     Contains namelists for parallel run control.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
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
    & config_testbed_mode       => testbed_mode,       &
    & config_testbed_iterations => testbed_iterations, &
    & null_mode
    
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_icon_testbed_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


       
  CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Adapted for I/O PEs, Rainer Johanni, Nov 2010
  !! Leonidas Linardakis, namelist restructuring, Jul 2011
  SUBROUTINE read_icon_testbed_namelist( filename )

    ! ------------------------------------------------------------------------
    INTEGER :: testbed_mode
    INTEGER :: testbed_iterations

    
    NAMELIST /testbed_nml/ testbed_mode, testbed_iterations

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat
    !0!CHARACTER(len=*), PARAMETER ::   &
    !0!        &  method_name = 'mo_parallel_nml:read_parallel_namelist'

    !--------------------------------------------
    ! set default values
    !--------------------------------------------
    testbed_mode       = null_mode
    testbed_iterations = 1
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
    config_testbed_mode       = testbed_mode
    config_testbed_iterations = testbed_iterations
    
  END SUBROUTINE read_icon_testbed_namelist
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed_nml
