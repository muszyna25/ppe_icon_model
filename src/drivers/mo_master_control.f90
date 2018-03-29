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

MODULE mo_master_control

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------

  USE mo_exception,     ONLY: message, finish
  USE mo_mpi,           ONLY: set_process_mpi_name, get_my_global_mpi_id, split_global_mpi_communicator
  USE mo_io_units,      ONLY: filename_max
  USE mo_master_config, ONLY: noOfModels, master_component_models, isRestart, read_restart_namelists
  USE mo_master_nml,    ONLY: read_master_namelist

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  init_master_control, master_namelist_filename,               &
    & get_my_namelist_filename, get_my_process_type, get_my_process_name, &
    & atmo_process, ocean_process, radiation_process, testbed_process,    &
    & my_process_is_ocean, get_my_model_no,                               &
    & are_multiple_models, use_restart_namelists, isRestart
   

  ! ------------------------------------------------------------------------
  INTEGER, PARAMETER :: atmo_process      = 1
  INTEGER, PARAMETER :: ocean_process     = 2
  INTEGER, PARAMETER :: radiation_process = 3
  INTEGER, PARAMETER :: testbed_process   = 99
  ! ------------------------------------------------------------------------

  INTEGER :: my_process_model ! =atmo_process,ocean_process,...
  INTEGER :: my_model_no ! 1,2,3  (id uniquely this process, even if it has the
                         ! same my_process_model with other compnents
                         ! Example: Two different components may run the dummy_process
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=64) :: my_model_name = ""

  CHARACTER(len=filename_max) :: master_namelist_filename = ""

  INTEGER :: my_model_min_rank, my_model_max_rank, my_model_inc_rank 
  LOGICAL :: multiple_models


CONTAINS

  !------------------------------------------------------------------------
  !>
  !!  Initialization of the master control variables
  !!
  INTEGER FUNCTION init_master_control(namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    INTEGER :: master_namelist_status
    INTEGER :: model_no, jg, js, je, jinc

    CHARACTER(LEN=*), PARAMETER :: method_name = "master_control"
    !-----------------------------------------------------------------------

    CALL message(method_name,'start model initialization.')

    master_namelist_status = read_master_namelist(TRIM(namelist_filename))

    !------------------------------------------------------------
    ! some checks

    IF (master_namelist_status == -1) THEN
      CALL finish(method_name,'model identity (atm/oce) can no longer be derived from namelist run_nml!')
    ENDIF

    IF (noOfModels() < 1) THEN
      CALL finish(method_name,'no of models < 1')
    ENDIF

    !------------------------------------------------------------

    master_namelist_filename = TRIM(namelist_filename)

    !------------------------------------------------------------
    ! find what is my process

    multiple_models = noOfModels() > 1    

    IF ( multiple_models ) THEN
      
      CALL set_my_component_null()

      COMPONENT_MODELS: DO model_no = 1, noOfModels()

        !         write(0,*) 'master_component_models:', model_no, trim(master_component_models(model_no)%model_name)

        js   = master_component_models(model_no)%model_min_rank
        je   = master_component_models(model_no)%model_max_rank
        jinc = master_component_models(model_no)%model_inc_rank
        
        DO jg = js, je, jinc

          IF ( get_my_global_mpi_id() == jg ) THEN
            
            CALL set_my_component(model_no,                                      &
                 &                master_component_models(model_no)%model_name,  &
                 &                master_component_models(model_no)%model_type,  &
                 &                master_component_models(model_no)%model_namelist_filename)
            
          ENDIF
          
        ENDDO
        
      ENDDO COMPONENT_MODELS

      CALL split_global_mpi_communicator ( my_model_no )

    ELSE ! only one component    

      model_no = 1
      CALL set_my_component(model_no,                                      &
           &                master_component_models(model_no)%model_name,  &
           &                master_component_models(model_no)%model_type,  &
           &                master_component_models(model_no)%model_namelist_filename)

    ENDIF

    !------------------------------------------------------------
    
    CALL check_my_component()
    
    !------------------------------------------------------------
    
    init_master_control = 0
    
  END FUNCTION init_master_control
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE set_my_component(comp_no, comp_name, comp_id, comp_namelist)

    INTEGER, INTENT(in)          :: comp_no
    CHARACTER(len=*), INTENT(in) :: comp_name
    INTEGER, INTENT(in)          :: comp_id
    CHARACTER(len=*), INTENT(in) :: comp_namelist

    my_model_no          = comp_no
    my_process_model     = comp_id
    my_namelist_filename = TRIM(comp_namelist)
    my_model_name        = TRIM(comp_name)

    my_model_min_rank    = master_component_models(comp_no)%model_min_rank
    my_model_max_rank    = master_component_models(comp_no)%model_max_rank
    my_model_inc_rank    = master_component_models(comp_no)%model_inc_rank

    CALL check_my_component()
    
    CALL set_process_mpi_name(TRIM(my_model_name))

  END SUBROUTINE set_my_component
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE check_my_component()

    CHARACTER(len=*), PARAMETER :: method_name='mo_master_control:check_my_component'

    IF (my_model_no < 1) CALL finish(method_name, 'my_model_no < 1') 
    IF (my_namelist_filename == '') CALL finish(method_name, 'my_namelist_filename = NULL')
    IF (my_model_name == '') CALL finish(method_name, 'my_model_name = NULL')

    SELECT CASE (my_process_model)
      CASE (atmo_process)
      CASE (ocean_process)
      CASE (radiation_process)
      CASE (testbed_process)
      CASE default
        CALL finish("check_my_component","my_process_model is unkown")
    END SELECT

  END SUBROUTINE check_my_component
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE set_my_component_null()

    my_model_no          = 0
    my_process_model     = -1
    my_namelist_filename = ''
    my_model_name        = ''
    my_model_min_rank    = -1
    my_model_max_rank    = -2
    my_model_inc_rank    = -1

  END SUBROUTINE set_my_component_null
  !------------------------------------------------------------------------

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
  LOGICAL FUNCTION my_process_is_ocean()

    my_process_is_ocean = (my_process_model == ocean_process)

  END FUNCTION my_process_is_ocean
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  LOGICAL FUNCTION use_restart_namelists()
    use_restart_namelists = (isRestart() .and. read_restart_namelists .and. (my_process_model /= ocean_process))
  END FUNCTION use_restart_namelists
  !------------------------------------------------------------------------

END MODULE mo_master_control
