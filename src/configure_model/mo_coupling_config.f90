!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_coupling_config

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

  LOGICAL :: config_coupled_mode

CONTAINS

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_run()

    is_coupled_run = config_coupled_mode

    WRITE ( 6 , * ) " Rene is_coupled_run ", is_coupled_run
    FLUSH ( 6 )

  END FUNCTION is_coupled_run
  !------------------------------------------------------------------------


END MODULE mo_coupling_config
