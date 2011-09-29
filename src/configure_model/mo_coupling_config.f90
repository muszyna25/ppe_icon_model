MODULE mo_coupling_config

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
     INTEGER            :: frequency
     INTEGER            :: time_step
     LOGICAL            :: l_diagnostic
     LOGICAL            :: l_activated
  END TYPE t_cpl_field_nml

  TYPE (t_cpl_field_nml), POINTER :: config_cpl_fields(:) => NULL()

END MODULE mo_coupling_config
