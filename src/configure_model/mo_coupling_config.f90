MODULE mo_coupling_config

  USE mo_master_control,  ONLY: are_multiple_models

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'
!>
!! Namelist input to steer the coulpling events
!!
!! Data type to store information about fields comming from the namelist
!!
!! - time_operation (average / accumate / none )
!! - couping_freq coupling frequency in seconds
!!
  TYPE t_cpl_field_nml
     CHARACTER(len=132) :: name
     INTEGER            :: lag
     LOGICAL            :: l_time_accumulation
     LOGICAL            :: l_time_average
     INTEGER            :: dt_coupling
     INTEGER            :: dt_model
     LOGICAL            :: l_diagnostic
     LOGICAL            :: l_activated
  END TYPE t_cpl_field_nml

  TYPE (t_cpl_field_nml), POINTER :: config_cpl_fields(:) => NULL()

  INTEGER :: number_of_coupled_variables
  INTEGER :: config_debug_coupler_level

CONTAINS

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_run()

    is_coupled_run = (number_of_coupled_variables > 0) .AND. are_multiple_models()

  END FUNCTION is_coupled_run
  !------------------------------------------------------------------------


END MODULE mo_coupling_config
