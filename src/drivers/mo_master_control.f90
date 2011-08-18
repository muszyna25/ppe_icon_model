!>
!! Provides the master control methods and paramaters
!! 
!!        
!! @par Revision History
!!   Created by Leonidas Linardakis, MPI-M, 2011-06-16
!!
!! @par Copyright
!! 2010-2011 by DWD and MPI-M
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

MODULE mo_master_control

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------

  USE mo_exception,          ONLY: message, finish
  USE mo_mpi,                ONLY: set_process_mpi_name, get_my_global_mpi_id, &
    &                              split_global_mpi_communicator

  USE mo_icon_cpl,           ONLY: get_cpl_local_comm
  USE mo_icon_cpl_init,      ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp, ONLY: icon_cpl_init_comp

  USE mo_io_units,           ONLY: filename_max
  
  USE mo_master_nml,         ONLY: read_master_namelist, lrestart, &
    & no_of_models, master_nml_array
    
  !USE mo_namelist,           ONLY: open_nml,  close_nml

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  init_master_control, get_my_namelist_filename,           &
    & get_my_process_type, get_my_process_name, is_coupled_run,  &
    & atmo_process, ocean_process, radiation_process, dummy_process,  &
    & my_process_is_ocean, is_restart_run, get_my_model_no,           &
    & get_my_couple_id


  ! ------------------------------------------------------------------------
  INTEGER, PARAMETER :: atmo_process  = 1
  INTEGER, PARAMETER :: ocean_process = 2
  INTEGER, PARAMETER :: radiation_process = 3
  INTEGER, PARAMETER :: dummy_process = 99
  ! ------------------------------------------------------------------------
  INTEGER :: my_process_model ! =atmo_process,ocean_process,...
  INTEGER :: my_model_no ! 1,2,3  (id uniquely this process, even if it has the
                         ! same my_process_model with other compnents
                         ! Example: Two different components may run the dummy_process
  INTEGER :: my_coupling_comp_id ! the coupling id for this component 
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=64) :: my_model_name
  
  INTEGER :: my_model_min_rank, my_model_max_rank, my_model_inc_rank 
  LOGICAL :: in_coupled_mode


  CONTAINS

  !------------------------------------------------------------------------
  !>
  !!  Initialization of the master control variables
  !!
  INTEGER FUNCTION init_master_control(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    !
    ! !Local variables
    !
    INTEGER :: master_namelist_status
    INTEGER :: model_no, jg, ierr

    CHARACTER(LEN=*), PARAMETER :: method_name = "master_control"
    !-----------------------------------------------------------------------
    
    CALL message(method_name,'start model initialization.')

    !------------------------------------------------------------
    master_namelist_status = read_master_namelist(TRIM(namelist_filename))
    
    !------------------------------------------------------------
    ! some checks
    IF (master_namelist_status == -1) THEN
      CALL finish(method_name,'model identity (atm/oce) can no longer '//&
                 'be derived from namelist run_nml!')
    ENDIF
    IF (no_of_models < 1) THEN
      CALL finish(method_name,'no_of_models < 1')
    ENDIF
    !------------------------------------------------------------

    !------------------------------------------------------------
    ! find what is my process
    in_coupled_mode = no_of_models > 1

    IF ( in_coupled_mode ) THEN

      CALL icon_cpl_init

      CALL set_my_component_null()
     
      DO model_no =1, no_of_models

!         write(0,*) 'master_nml_array:', model_no, master_nml_array(model_no)%model_name        
        DO jg = master_nml_array(model_no)%model_min_rank,&
          & master_nml_array(model_no)%model_max_rank,&
          & master_nml_array(model_no)%model_inc_rank

          IF ( get_my_global_mpi_id() == jg ) THEN
            
            CALL set_my_component(model_no,             &
               & master_nml_array(model_no)%model_name, &
               & master_nml_array(model_no)%model_type ,&
               & master_nml_array(model_no)%model_namelist_filename)
                 
          ENDIF

        ENDDO
          
      ENDDO !model_no =1, no_of_models

      
    ELSE
      ! only one component    
      model_no=1
      CALL set_my_component(model_no,             &
          & master_nml_array(model_no)%model_name, &
          & master_nml_array(model_no)%model_type ,&
          & master_nml_array(model_no)%model_namelist_filename)
               
    ENDIF
    !------------------------------------------------------------
    ! check if my component is ok
    CALL check_my_component()

    !------------------------------------------------------------
    IF ( in_coupled_mode ) THEN
      ! Inform the coupler about what we are
      CALL icon_cpl_init_comp ( my_model_name, my_model_no, my_process_model, &
          &                     my_coupling_comp_id, ierr )
      ! split the global_mpi_communicator into the components
      CALL split_global_mpi_communicator ( my_model_no )
    ENDIF
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
    my_model_name = TRIM(comp_name)
    
    my_model_min_rank    = master_nml_array(comp_no)%model_min_rank
    my_model_max_rank    = master_nml_array(comp_no)%model_max_rank
    my_model_inc_rank    = master_nml_array(comp_no)%model_inc_rank

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
      CASE (dummy_process)
      CASE default
        CALL finish("set_my_component","my_process_model is unkown")
    END SELECT

  END SUBROUTINE check_my_component
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE set_my_component_null()
  
    my_model_no          = 0
    my_process_model     = 0
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
  INTEGER FUNCTION get_my_couple_id()

    get_my_couple_id = my_coupling_comp_id 
    
  END FUNCTION get_my_couple_id
  !------------------------------------------------------------------------
  
  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_run()

    is_coupled_run = in_coupled_mode
    
  END FUNCTION is_coupled_run
  !------------------------------------------------------------------------
  
  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_restart_run()

    is_restart_run = lrestart
    
  END FUNCTION is_restart_run
  !------------------------------------------------------------------------
  
  !------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_ocean()

    my_process_is_ocean = (my_process_model == ocean_process)
    
  END FUNCTION my_process_is_ocean
  !------------------------------------------------------------------------

END MODULE mo_master_control
