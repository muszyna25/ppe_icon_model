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
  TYPE t_field_nml
     CHARACTER(len=132) :: name
     INTEGER            :: lag
     LOGICAL            :: l_time_accumulation
     LOGICAL            :: l_time_average
     INTEGER            :: frequency
     INTEGER            :: time_step
  END TYPE t_field_nml

  TYPE (t_field_nml), POINTER :: config_fields(:) => NULL()

END MODULE mo_coupling_config
