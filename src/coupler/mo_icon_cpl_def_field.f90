!>
!! Definition of the ICON fields for coupling
!!
!! <Description>
!! 
!! @author Rene Redler, MPI-M
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_icon_cpl_def_field

  USE mo_icon_cpl, ONLY         : ICON_comm,         & ! MPI communicator
     &                            t_cpl_field,       & ! Field type
     &                            t_comp,            & ! Field type
     &                            cpl_fields,        & ! exchange fields
     &                            nbr_ICON_inc,      & ! increment for memory
     &                            nbr_active_fields, & ! total number of coupling fields
     &                            nbr_ICON_fields,   & ! allocated size for fields
     &                            cpl_field_none,    &
     &                            cpl_field_avg,     &
     &                            cpl_field_acc

  USE mo_icon_cpl_restart, ONLY : cpl_init_restart
  USE mo_coupling_config, ONLY  : config_cpl_fields
  USE mo_event_manager, ONLY    : event_add

  IMPLICIT NONE

  PRIVATE

  TYPE(t_cpl_field), POINTER    :: fptr
  TYPE(t_cpl_field), POINTER    :: new_cpl_fields(:)
  INTEGER                       :: global_field_id
  INTEGER                       :: i
  INTEGER                       :: new_dim

  ! Return code for error handling

  INTEGER                       :: ierr

  PUBLIC :: ICON_cpl_def_field, ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids

CONTAINS

  SUBROUTINE ICON_cpl_def_field ( field_name, grid_id, field_id, field_shape, ierror )

    CHARACTER(LEN=*), INTENT (in) :: field_name

    INTEGER, INTENT(in)           :: grid_id          !<  grid id
    INTEGER, INTENT(out)          :: field_id         !<  field id
    INTEGER, INTENT(in)           :: field_shape(3)
    INTEGER, INTENT(out)          :: ierror           !<  returned error code

    ! ----------------
    ! Local variables:
    ! ----------------

    INTEGER                       ::  nbr_max_fields  !< maximum number of fields

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    ierror = 0

    ! -------------------------------------------------------------------
    ! Initialise field
    ! -------------------------------------------------------------------

    DO field_id = 1, nbr_ICON_fields
       IF ( .NOT. cpl_fields(field_id)%l_field_status ) EXIT
    ENDDO

    IF ( field_id > nbr_ICON_fields ) THEN

       ! We need to allocated more memory

       new_dim = nbr_ICON_fields + nbr_ICON_inc

       ALLOCATE ( new_cpl_fields(new_dim), STAT = ierr )
       IF ( ierr > 0 ) THEN
          WRITE ( * , * ) ' Error allocating fields '
#ifndef NOMPI
          CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
       ENDIF

       new_cpl_fields (1:nbr_ICON_fields) = &
           cpl_fields (1:nbr_ICON_fields)

       DO i = nbr_ICON_fields+1, new_dim
          new_cpl_fields(i)%grid_id         = 0
          new_cpl_fields(i)%global_field_id = 0
          new_cpl_fields(i)%field_shape     = 0
          new_cpl_fields(i)%l_field_status  = .FALSE.

          new_cpl_fields(i)%coupling%cdi_varID    = -1
          new_cpl_fields(i)%coupling%cdi_gridID   = -1
          new_cpl_fields(i)%coupling%restart_flag = .FALSE.
          new_cpl_fields(i)%coupling%lag          = 0
          new_cpl_fields(i)%coupling%dt_coupling  = 0
          new_cpl_fields(i)%coupling%dt_model     = 0
          new_cpl_fields(i)%coupling%time_operation = 0

          Nullify ( new_cpl_fields(i)%send_field_acc )

       ENDDO

       !
       !   De-allocate fields vector
       !
       DEALLOCATE ( cpl_fields, STAT = ierr )
       IF (ierr > 0) THEN
          WRITE ( * , * ) ' Error deallocating fields '
#ifndef NOMPI
          CALL MPI_Abort ( ICON_comm, 1, ierr )
#endif
       ENDIF

       ! Assign field ID

       field_id = nbr_ICON_fields + 1

       ! Update field pointer

       cpl_fields => new_cpl_fields

       ! Update size of Field type

       nbr_ICON_fields = new_dim

    ENDIF

    fptr => cpl_fields(field_id)

    ! -------------------------------------------------------------------
    ! Determine global field_id
    !
    ! Currently the global field id is determined from the position in the
    ! namelist file. This requires that all components need to have the
    ! same namelist file for coupling. Possibly this contains coupling
    ! fields which are not handled by an individual component.
    ! On the long term e need to find a different way to get a global
    ! field id, e.g. through the hash value of its character string
    !
    ! Note that the global field id is later used for the MPI message tag
    ! which ensures that a field sent as sea surface temperature will be
    ! received as sea surface temperature.
    !
    ! -------------------------------------------------------------------

    nbr_max_fields = SIZE(config_cpl_fields)

    ! -------------------------------------------------------------------
    ! Rather than doing a string comparision we could used hash values
    ! -------------------------------------------------------------------

    DO global_field_id = 1, nbr_max_fields
       IF ( TRIM(field_name) == TRIM(config_cpl_fields(global_field_id)%name) ) THEN
          WRITE ( * , * ) ' Global Field ID is: ', global_field_id
          EXIT
       ENDIF
    ENDDO

    IF ( global_field_id > nbr_max_fields ) THEN
       WRITE ( * , * ) ' Fieldname not declared in namelist: ', TRIM(field_name)
       WRITE ( * , * ) ' Field will not be exchanged! ', TRIM(field_name)
       field_id = -1
       RETURN
    ENDIF

    ! -------------------------------------------------------------------
    ! Store grid parameters and global index list
    ! -------------------------------------------------------------------

    fptr%l_field_status  = .TRUE.

    fptr%field_name      = TRIM(field_name)

    fptr%grid_id         = grid_id
    fptr%global_field_id = global_field_id

    fptr%field_shape     = field_shape
    nbr_active_fields    = field_id

    ! -------------------------------------------------------------------
    ! Initialize coupling, substitute for the OASIS4 XML reading and storage
    ! -------------------------------------------------------------------

    IF ( config_cpl_fields(global_field_id)%l_time_accumulation ) THEN
       fptr%coupling%time_operation = cpl_field_acc
    ELSE IF ( config_cpl_fields(global_field_id)%l_time_average ) THEN
       fptr%coupling%time_operation = cpl_field_avg
    ELSE
       fptr%coupling%time_operation = cpl_field_none
    ENDIF

    fptr%coupling%cdi_gridID   = -1
    fptr%coupling%cdi_varID    = -1
    fptr%coupling%restart_flag = .FALSE.
    fptr%coupling%lag          = config_cpl_fields(global_field_id)%lag
    fptr%coupling%dt_coupling  = config_cpl_fields(global_field_id)%dt_coupling
    fptr%coupling%dt_model     = config_cpl_fields(global_field_id)%dt_model
    fptr%coupling%l_activated  = config_cpl_fields(global_field_id)%l_activated
    fptr%coupling%l_diagnostic = config_cpl_fields(global_field_id)%l_diagnostic
  
    ! -------------------------------------------------------------------
    ! Prepare the restarting for coupling fields
    ! -------------------------------------------------------------------

    IF ( fptr%coupling%time_operation /= cpl_field_none ) &
       call cpl_init_restart ( field_id, ierror )

    ! -------------------------------------------------------------------
    ! Signal new coupling event and store event_id in  fptr%event_id
    ! -------------------------------------------------------------------

    CALL event_add ( fptr%event_id, fptr%coupling%dt_coupling, &
                     fptr%coupling%dt_model, fptr%coupling%lag )

  END SUBROUTINE ICON_cpl_def_field

  SUBROUTINE  ICON_cpl_get_nbr_fields ( nbr_of_fields )
     INTEGER, INTENT(out) :: nbr_of_fields
     nbr_of_fields = nbr_active_fields
  END SUBROUTINE ICON_cpl_get_nbr_fields

  SUBROUTINE  ICON_cpl_get_field_ids ( nbr_of_fields, field_ids )
     INTEGER, INTENT(IN)  :: nbr_of_fields
     INTEGER, INTENT(out) :: field_ids(nbr_of_fields)
     INTEGER i
     DO i = 1, nbr_active_fields
        field_ids(i) = i
     ENDDO
  END SUBROUTINE ICON_cpl_get_field_ids

END MODULE mo_icon_cpl_def_field
