!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_coupling_config

#ifndef YAC_coupling
  USE mo_master_control,  ONLY: are_multiple_models
#endif

  IMPLICIT NONE
  PUBLIC
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

#ifdef YAC_coupling
  LOGICAL :: config_coupled_mode
#endif

CONTAINS

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_run()

#ifdef YAC_coupling
    is_coupled_run = config_coupled_mode
#else
    is_coupled_run = (number_of_coupled_variables > 0) .AND. are_multiple_models()
#endif

  END FUNCTION is_coupled_run
  !------------------------------------------------------------------------


END MODULE mo_coupling_config
