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
    &                              set_process_mpi_communicator

  USE mo_icon_cpl,           ONLY: get_cpl_local_comm!, complist
  USE mo_icon_cpl_init,      ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp, ONLY: icon_cpl_init_comp

  USE mo_io_units,           ONLY: filename_max
  
  USE mo_master_nml,         ONLY: read_master_namelist, lrestart, &
    & no_of_models, master_nml_array
    
  !USE mo_namelist,           ONLY: open_nml,  close_nml

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  init_master_control, get_my_namelist_filename,           &
    & get_my_process_component, get_my_process_name, is_coupled_run,  &
    & atmo_process, ocean_process, radiation_process, dummy_process,  &
    & my_process_is_ocean, is_restart_run


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
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=64) :: my_model_name
  
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
    INTEGER :: model_no, jg, comp_id, str_len, ierr
    INTEGER :: new_comm

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

        write(0,*) 'master_nml_array:', model_no, master_nml_array(model_no)%model_name
        
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
      ! Inform tghe coupler of what we are
      CALL icon_cpl_init_comp ( my_model_name, my_process_model, comp_id, ierr )
      ! make the component communicator available for use within the ICON components
      new_comm = get_cpl_local_comm()
      CALL set_process_mpi_communicator ( new_comm )
    ENDIF   
    !------------------------------------------------------------

    
    !------------------------------------------------------------

!     IF ( in_coupled_mode .AND. atmo_namelist_filename == ocean_namelist_filename ) THEN
!        WRITE ( * , * ) ' Component specific namelists have to be different!'
!     ENDIF

    !------------------------------------------------------------
!     Leonidas: This has to be moved and use the master namelist

!     IF ( in_coupled_mode ) THEN
! 
!        WRITE (complist(comp_id)%nml_name,'(A13)') 'NAMELIST_ICON'
! 
!        SELECT CASE ( my_process_model )
! 
!        CASE ( atmo_process )
! 
!           complist(comp_id)%comp_name      = TRIM(atmo_name)
!           complist(comp_id)%comp_process   = atmo_process
!           complist(comp_id)%min_rank       = atmo_min_rank
!           complist(comp_id)%max_rank       = atmo_max_rank
!           complist(comp_id)%inc_rank       = atmo_inc_rank
! 
!           IF (TRIM(atmo_namelist_filename) /= "") THEN
!              str_len = LEN_TRIM(atmo_namelist_filename)
!              WRITE (complist(comp_id)%nml_name(14:14),'(A1)') '_'
!              WRITE (complist(comp_id)%nml_name(15:15+str_len),'(A)') TRIM(atmo_namelist_filename)
!           ENDIF
! 
!        CASE ( ocean_process )
! 
!           complist(comp_id)%comp_name      = TRIM(ocean_name)
!           complist(comp_id)%l_comp_status  = l_ocean_active
!           complist(comp_id)%min_rank       = ocean_min_rank
!           complist(comp_id)%max_rank       = ocean_max_rank
!           complist(comp_id)%inc_rank       = ocean_inc_rank
! 
!           IF (TRIM(ocean_namelist_filename) /= "") THEN
!              str_len = LEN_TRIM(ocean_namelist_filename)
!              WRITE (complist(comp_id)%nml_name(14:14),'(A1)') '_'
!              WRITE (complist(comp_id)%nml_name(15:15+str_len),'(A)') TRIM(ocean_namelist_filename)
!           ENDIF
! 
!        END SELECT
! 
!     ENDIF
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
  INTEGER FUNCTION get_my_process_component()

    get_my_process_component = my_process_model
    
  END FUNCTION get_my_process_component
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
