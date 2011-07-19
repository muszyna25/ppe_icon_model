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

  USE mo_exception,  ONLY: warning, message, finish
  USE mo_mpi,        ONLY: set_process_mpi_name, get_my_global_mpi_id, &
    &                      set_process_mpi_communicator

  USE mo_icon_cpl,   ONLY: get_cpl_local_comm, complist
  USE mo_icon_cpl_init, ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp, ONLY: icon_cpl_init_comp

  USE mo_io_units,    ONLY: filename_max, nnml
  
  USE mo_master_nml,  ONLY: read_master_namelist,                             &
    &                   l_atmo_active, atmo_name, atmo_namelist_filename,     &
    &                   atmo_min_rank, atmo_max_rank, atmo_inc_rank,          &
    &                   l_ocean_active, ocean_name, ocean_namelist_filename,  &
    &                   ocean_min_rank, ocean_max_rank, ocean_inc_rank
  
  USE mo_coupling_nml, ONLY: read_coupling_namelist

  USE mo_namelist,     ONLY: open_nml,  close_nml

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  init_master_control, get_my_namelist_filename,           &
    & get_my_process_component, get_my_process_name, is_coupled_run,  &
    & atmo_process, ocean_process, radiation_process


  ! ------------------------------------------------------------------------
  INTEGER, PARAMETER :: atmo_process  = 1
  INTEGER, PARAMETER :: ocean_process = 2
  INTEGER, PARAMETER :: radiation_process = 3
  ! ------------------------------------------------------------------------
  INTEGER :: my_process_model
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
    INTEGER :: coupling_namelist_status
    INTEGER :: jg, comp_id, str_len, ierr
    INTEGER :: new_comm
    INTEGER :: nbr_components

    CHARACTER(LEN=*), PARAMETER :: method_name = "master_cotrol"
    !-----------------------------------------------------------------------
    
    CALL message(method_name,'start model initialization.')

    !------------------------------------------------------------
    master_namelist_status = read_master_namelist(TRIM(namelist_filename))
    
    !------------------------------------------------------------

    IF (master_namelist_status == -1) THEN
      ! we are running an old experiment with no master description, no coupling
      in_coupled_mode = .false.
      ! both ocean and atmo read the same namelist
      my_namelist_filename = "NAMELIST_ICON"
     
      CALL finish(method_name,'model identity (atm/oce) can no longer '//&
                 'be derived from namelist run_nml!')
 
     ! read the run_nml in order to figure out which component we run
     !CALL open_nml('NAMELIST_ICON')
     !CALL run_nml_setup
     !CALL close_nml

     !IF (locean) THEN
     !  CALL message(method_name,'ocean_process')
     !  l_ocean_active = .true.
     !  my_process_model = ocean_process
     !ELSE
     !  CALL message(method_name,'atmo_process')
     !  l_atmo_active = .true.
     !  my_process_model = atmo_process
     !ENDIF
     !init_master_control = -1  ! did not find namelist
     !RETURN

   ELSE

      nbr_components = 0

      IF ( l_atmo_active  ) nbr_components = nbr_components + 1
      IF ( l_ocean_active ) nbr_components = nbr_components + 1

      in_coupled_mode = nbr_components > 1

      IF ( in_coupled_mode ) THEN

         CALL icon_cpl_init

         DO jg = ocean_min_rank, ocean_max_rank, ocean_inc_rank

            IF ( get_my_global_mpi_id() == jg ) THEN

               CALL set_my_component("OCEAN", ocean_process , ocean_namelist_filename)
               CALL icon_cpl_init_comp ( 'oce', my_process_model, comp_id, ierr )

!rr  will be moved into mo_ocean_model.f90
!rr               coupling_namelist_status = &
!rr                   read_coupling_namelist(TRIM(ocean_namelist_filename),comp_id)

            ENDIF

         ENDDO

         DO jg = atmo_min_rank, atmo_max_rank, atmo_inc_rank

            IF ( get_my_global_mpi_id() == jg ) THEN

               CALL set_my_component("ATMO", atmo_process ,atmo_namelist_filename)
               CALL icon_cpl_init_comp ( 'atm', my_process_model, comp_id, ierr )

!rr will be moved into mo_atmo_model.f90
!rr               coupling_namelist_status = &
!rr                   read_coupling_namelist(TRIM(atmo_namelist_filename),comp_id)

            ENDIF

         ENDDO

         ! make the component communicator available for use within the ICON components

         new_comm = get_cpl_local_comm()

         CALL set_process_mpi_communicator ( new_comm )

      ELSE

         IF (l_atmo_active.AND.(.NOT.l_ocean_active)) THEN

            CALL set_my_component("ATMO", atmo_process ,atmo_namelist_filename)

         ELSE IF ((.NOT.l_atmo_active).AND.l_ocean_active) THEN

            CALL set_my_component("OCEAN", ocean_process , ocean_namelist_filename)

         ELSE

            CALL finish(method_name,'check master namelist setup!')

         END IF
      ENDIF

    ENDIF
    !------------------------------------------------------------

    
    !------------------------------------------------------------

    IF ( in_coupled_mode .AND. atmo_namelist_filename == ocean_namelist_filename ) THEN
       WRITE ( * , * ) ' Component specific namelists have to be different!'
    ENDIF

    IF ( in_coupled_mode ) THEN

       WRITE (complist(comp_id)%nml_name,'(A13)') 'NAMELIST_ICON'

       SELECT CASE ( my_process_model )

       CASE ( atmo_process )

          complist(comp_id)%comp_name      = TRIM(atmo_name)
          complist(comp_id)%comp_process   = atmo_process
          complist(comp_id)%min_rank       = atmo_min_rank
          complist(comp_id)%max_rank       = atmo_max_rank
          complist(comp_id)%inc_rank       = atmo_inc_rank

          IF (TRIM(atmo_namelist_filename) /= "") THEN
             str_len = LEN_TRIM(atmo_namelist_filename)
             WRITE (complist(comp_id)%nml_name(14:14),'(A1)') '_'
             WRITE (complist(comp_id)%nml_name(15:15+str_len),'(A)') TRIM(atmo_namelist_filename)
          ENDIF

       CASE ( ocean_process )

          complist(comp_id)%comp_name      = TRIM(ocean_name)
          complist(comp_id)%l_comp_status  = l_ocean_active
          complist(comp_id)%min_rank       = ocean_min_rank
          complist(comp_id)%max_rank       = ocean_max_rank
          complist(comp_id)%inc_rank       = ocean_inc_rank

          IF (TRIM(ocean_namelist_filename) /= "") THEN
             str_len = LEN_TRIM(ocean_namelist_filename)
             WRITE (complist(comp_id)%nml_name(14:14),'(A1)') '_'
             WRITE (complist(comp_id)%nml_name(15:15+str_len),'(A)') TRIM(ocean_namelist_filename)
          ENDIF

       END SELECT

    ENDIF

    init_master_control = 0
    
  END FUNCTION init_master_control
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE set_my_component(comp_name, comp_id, comp_namelist)
    
    CHARACTER(len=*), INTENT(in) :: comp_name
    INTEGER, INTENT(in)          :: comp_id
    CHARACTER(len=*), INTENT(in) :: comp_namelist
       
    my_process_model     = comp_id
    my_namelist_filename = TRIM(comp_namelist)
    my_model_name = TRIM(comp_name)
    CALL set_process_mpi_name(TRIM(my_model_name))
    
    SELECT CASE (my_process_model)
      CASE (atmo_process)
      CASE (ocean_process)
      CASE (radiation_process)
      CASE default
        CALL finish("set_my_component","my_process_model is unkown")      
    END SELECT
    
  END SUBROUTINE set_my_component
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

END MODULE mo_master_control
