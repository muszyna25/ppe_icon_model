!>
!! Provides the master control methods and paramaters
!! 
!!        
!! @par Revision History
!!   Created by Leonidas Linardakis, MPI-M, 2011-06-16
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!


!--------------------------------------------------------------------------------------------------------
MODULE mo_master_control
 
  USE mo_io_units,      ONLY: filename_max
  USE mo_master_config
  
  PUBLIC :: master_namelist_filename,                                     &
    & get_my_namelist_filename, get_my_process_type, get_my_process_name, &
    & atmo_process, ocean_process, ps_radiation_process, testbed_process, &
    & my_process_is_ocean, get_my_model_no,                               &
    & are_multiple_models, use_restart_namelists, isRestart,              &
    & process_exists, hamocc_process, my_process_is_hamocc,               &
    & my_process_is_oceanic, icon_output_process, my_process_is_icon_output

CONTAINS

  !------------------------------------------------------------------------
  CHARACTER(len=64) FUNCTION get_my_process_name()

    get_my_process_name = my_model_name

  END FUNCTION get_my_process_name
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  CHARACTER(len=filename_max) FUNCTION get_my_namelist_filename()

    get_my_namelist_filename = my_namelist_filename

  END FUNCTION get_my_namelist_filename
  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  INTEGER FUNCTION get_my_process_type()

    get_my_process_type = my_process_model

  END FUNCTION get_my_process_type
  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  INTEGER FUNCTION get_my_model_no()

    get_my_model_no = my_model_no

  END FUNCTION get_my_model_no
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION are_multiple_models()

    are_multiple_models = multiple_models

  END FUNCTION are_multiple_models
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION process_exists(process_type)
    INTEGER, INTENT(in) :: process_type

    INTEGER :: model_no

    process_exists = .false.
    DO model_no = 1, noOfModels()

      IF (process_type == master_component_models(model_no)%model_type) THEN
        process_exists = .true.
        RETURN
      ENDIF

    ENDDO
    RETURN

  END FUNCTION process_exists
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_ocean()

    my_process_is_ocean = (my_process_model == ocean_process)

  END FUNCTION my_process_is_ocean
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_hamocc()

    my_process_is_hamocc = (my_process_model == hamocc_process)

  END FUNCTION my_process_is_hamocc
  !------------------------------------------------------------------------
 
  !------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_oceanic()

    my_process_is_oceanic = (my_process_model == ocean_process) .or. (my_process_model == hamocc_process) &
      & .or. (my_process_model == icon_output_process) ! FixMe: temporary to make it work with the ocean 

  END FUNCTION my_process_is_oceanic
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_icon_output()

    my_process_is_icon_output = (my_process_model == icon_output_process)

  END FUNCTION my_process_is_icon_output
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION use_restart_namelists()
    use_restart_namelists = (isRestart() .and. read_restart_namelists .and. &
      & (my_process_model /= ocean_process) .and. (my_process_model /= hamocc_process))
  END FUNCTION use_restart_namelists
  !------------------------------------------------------------------------

END MODULE mo_master_control
