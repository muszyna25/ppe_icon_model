!>
!! Routine to check whether we need to provide a restart field
!! for the first coupling step
!!
!! <Description>
!! 
!! @author Rene Redler, MPI-M
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2011-12-16)
!!
!! @par Copyright
!! 2011-2012 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and WARRANTY conditions.
!! 
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! &ltol>
!! &ltli> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! &ltli> The code may not be re-distributed without the consent of the authors.
!! &ltli> The copyright notice and statement of authorship must appear in all
!!    copies.
!! &ltli> You accept the warranty conditions (see WARRANTY).
!! &ltli> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_icon_cpl_restart
#ifndef NOMPI
  USE mpi
  USE mo_icon_cpl, ONLY       : ICON_root, ICON_comp_comm, &
       &                        ICON_local_size, ICON_local_rank, &
       &                        cpl_field_none, datatype
#endif
  USE mo_kind, ONLY           : wp
  USE mo_io_units, ONLY       : find_next_free_unit
  USE mo_icon_cpl, ONLY       : t_comp, t_grid, t_cpl_field, &
       &                        comps, grids, cpl_fields,    &
       &                        cplout, debug_coupler_level
  USE mo_time_config, ONLY    : time_config

  IMPLICIT NONE

  PRIVATE
 
  TYPE(t_grid), POINTER      :: gptr ! pointer to grid struct
  TYPE(t_cpl_field), POINTER :: fptr ! pointer to field struct

  INTEGER                :: len
  INTEGER                :: rest_unit
  INTEGER                :: year, month, day, hour, minute, second

  CHARACTER(len=132) :: file_name

  PUBLIC :: cpl_init_restart, cpl_write_restart, cpl_read_restart 

CONTAINS

  ! ---------------------------------------------------------------------

  SUBROUTINE cpl_init_restart ( field_id, ierror )

    INTEGER, INTENT(IN)  :: field_id
    INTEGER, INTENT(out) :: ierror           !<  returned error code

#ifndef NOMPI
    INTEGER :: local_size
    INTEGER :: i

    fptr => cpl_fields(field_id)
    gptr => grids(fptr%grid_id)

    local_size = fptr%field_shape(2) - fptr%field_shape(1) + 1

    ! Collect local sizes of coupling fields

    IF ( ICON_local_rank == ICON_root ) THEN
       ALLOCATE(fptr%local_sizes (ICON_local_size))
       ALLOCATE(fptr%displacement(ICON_local_size))
    ELSE
       ALLOCATE(fptr%local_sizes (1))
       ALLOCATE(fptr%displacement(1))
       ALLOCATE(fptr%global_index(1))
    ENDIF

    CALL MPI_GATHER ( local_size, 1, MPI_INTEGER, fptr%local_sizes, 1, MPI_INTEGER, &
                      ICON_root, ICON_comp_comm, ierror )

    IF ( ICON_local_rank == ICON_root ) THEN

       ! Compute offsets

       fptr%displacement(1) = 0

       DO i = 2, ICON_local_size
          fptr%displacement(i) = fptr%displacement(i-1) + fptr%local_sizes(i-1)
       ENDDO

       fptr%global_size = SUM(fptr%local_sizes)

       ALLOCATE (fptr%global_index(fptr%global_size))

    ENDIF

    ! Collect global indices

    CALL MPI_GATHERV ( gptr%grid_glob_index, local_size, MPI_INTEGER, &
                       fptr%global_index, fptr%local_sizes, fptr%displacement, MPI_INTEGER, &
                       ICON_root, ICON_comp_comm, ierror )
#else
    fptr => cpl_fields(field_id)
    ierror = 0
#endif

  END SUBROUTINE cpl_init_restart

  ! ---------------------------------------------------------------------

  SUBROUTINE cpl_write_restart ( field_id, field_shape, coupling_field, count, ierror )

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of send field
    INTEGER, INTENT(in)    :: count            !<  number of accumulations

    REAL (wp), INTENT(in)  :: coupling_field (field_shape(1):field_shape(2),field_shape(3))

    INTEGER, INTENT(out)   :: ierror           !<  returned error code
#ifndef NOMPI
    REAL(wp), ALLOCATABLE  :: global_field (:,:)

    INTEGER :: i, j, local_size

    ierror = 0

    fptr => cpl_fields(field_id)

    IF ( fptr%coupling%time_operation == cpl_field_none ) RETURN

    local_size = field_shape(2) - field_shape(1) + 1

    IF ( ICON_local_rank == ICON_root ) THEN
       ALLOCATE (global_field (fptr%global_size,field_shape(3)))
    ELSE
       ALLOCATE (global_field (1,field_shape(3)))
    ENDIF

    rest_unit = find_next_free_unit(10,100)

    IF ( debug_coupler_level > 2 ) THEN
       WRITE(file_name(1:22), '(A15,A4,A1,I2.2)') 'debug_wrte_rst_', &
               fptr%field_name, '_', ICON_local_rank
       OPEN ( UNIT = rest_unit, FILE = file_name(1:22), FORM = "FORMATTED", STATUS = "UNKNOWN" )
       DO j = 1, field_shape(3)
       DO i = field_shape(1), field_shape(2)
         WRITE ( rest_unit ,  * ) j, i, coupling_field(i,j)
       ENDDO
       ENDDO
       CLOSE ( UNIT = rest_unit)
    ENDIF

    ! Could we use MPI_Type_vector?

    DO i = 1, field_shape(3)
       CALL MPI_GATHERV ( coupling_field(field_shape(1),i), local_size, datatype, &
                          global_field(:,i), fptr%local_sizes, fptr%displacement, datatype, &
                          ICON_root, ICON_comp_comm, ierror )
    ENDDO

    IF ( ICON_local_rank == ICON_root ) THEN

    second = INT(time_config%end_datetime%second)
    len = 12 + LEN_TRIM(comps(1)%comp_name) + 1 + 16
    WRITE(file_name(1:len), '(A12,A,A1,I4.4,2I2.2,A1,3I2.2,A1)')   &
          "restart_cpl_", TRIM(comps(1)%comp_name), "_", &
          time_config%end_datetime%year,     &
          time_config%end_datetime%month,    &
          time_config%end_datetime%day, "T", &
          time_config%end_datetime%hour,     &
          time_config%end_datetime%minute,   &
          second, "Z"

    OPEN   ( UNIT = rest_unit, FILE = file_name(1:len), FORM = "UNFORMATTED", &
                                      STATUS = "UNKNOWN", POSITION = "APPEND" )

    WRITE  ( unit = rest_unit ) fptr%global_field_id
    WRITE  ( unit = rest_unit ) time_config%end_datetime%year,   &
                                time_config%end_datetime%month,  &
                                time_config%end_datetime%day,    &
                                time_config%end_datetime%hour,   &
                                time_config%end_datetime%minute, &
                                time_config%end_datetime%second
    WRITE  ( unit = rest_unit ) count
    WRITE  ( unit = rest_unit ) fptr%global_size, field_shape(3)

    WRITE  ( unit = rest_unit ) global_field

    CLOSE  ( unit = rest_unit )

    ENDIF
    DEALLOCATE (global_field)
#else
    fptr => cpl_fields(field_id)
    print *, 'Restart is not supported', coupling_field (field_shape(1),field_shape(3)), count
    ierror = 0
#endif

  END SUBROUTINE cpl_write_restart

  ! ---------------------------------------------------------------------

  SUBROUTINE cpl_read_restart ( field_id, field_shape, coupling_field, info, ierror )

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of send field

    REAL (wp), INTENT(out) :: coupling_field (field_shape(1):field_shape(2),field_shape(3))

    INTEGER, INTENT(out)   :: info             !<  returned info code
    INTEGER, INTENT(out)   :: ierror           !<  returned error code

#ifndef NOMPI
    REAL(wp), ALLOCATABLE  :: global_field (:,:)

    INTEGER :: i, j, local_size

    INTEGER :: year, month, day, hour, minute, second
    INTEGER :: nbr_bundles
    INTEGER :: global_field_id
    INTEGER :: count
    INTEGER :: eof
    INTEGER :: bcast_buffer(2)
    LOGICAL :: existing

    fptr => cpl_fields(field_id)

    ierror = 0
    info   = 0
    eof    = 1

    IF ( fptr%coupling%time_operation == cpl_field_none ) RETURN

    second = INT(time_config%cur_datetime%second)
    len = 12 + LEN_TRIM(comps(1)%comp_name) + 1 + 16
    WRITE(file_name(1:len), '(A12,A,A1,I4.4,2I2.2,A1,3I2.2,A1)')   &
         "restart_cpl_", TRIM(comps(1)%comp_name), "_", &
         time_config%cur_datetime%year,     &
         time_config%cur_datetime%month,    &
         time_config%cur_datetime%day, "T", &
         time_config%cur_datetime%hour,     &
         time_config%cur_datetime%minute,   &
         second, "Z"

    INQUIRE ( FILE=TRIM(file_name), EXIST = existing )

    IF ( debug_coupler_level > 0 ) &
       WRITE ( cplout , '(A20,A,L2)' ) 'Trying to open file ', file_name(1:len), existing

    IF ( .NOT. existing ) RETURN

    IF ( ICON_local_rank == ICON_root ) THEN
       rest_unit = find_next_free_unit(10,100)
       OPEN  ( UNIT = rest_unit, FILE = file_name(1:len), FORM = "UNFORMATTED", status = 'OLD' )
    ENDIF

    local_size = field_shape(2) - field_shape(1) + 1

    IF ( ICON_local_rank == ICON_root ) THEN
       ALLOCATE (global_field (fptr%global_size,field_shape(3)))
    ELSE
       ALLOCATE (global_field (               1,field_shape(3)))
    ENDIF

    IF ( ICON_local_rank == ICON_root ) THEN

       global_field_id = -1

       DO WHILE ( eof >=0 ) 
          READ  ( unit = rest_unit, IOSTAT=eof ) global_field_id
          IF ( global_field_id == fptr%global_field_id ) THEN
             READ  ( unit = rest_unit, IOSTAT=eof ) year, month, day, hour, minute, second
             READ  ( unit = rest_unit, IOSTAT=eof ) count
             READ  ( unit = rest_unit ) fptr%global_size, nbr_bundles
             READ  ( unit = rest_unit, IOSTAT=eof ) global_field
             IF ( debug_coupler_level > 0 ) THEN
                 WRITE ( cplout , '(a)' ) ' reading restart ...'
                 WRITE ( cplout , '(a9,2i8)' ) ' with yr ', year,  time_config%cur_datetime%year
                 WRITE ( cplout , '(a9,2i8)' ) ' with mo ', month, time_config%cur_datetime%month
                 WRITE ( cplout , '(a9,2i8)' ) ' with dy ', day,   time_config%cur_datetime%day
                 WRITE ( cplout , '(a9,2i8)' ) ' with hr ', hour,  time_config%cur_datetime%hour
                 WRITE ( cplout , '(a9,2i8)' ) ' with mn', minute, time_config%cur_datetime%minute
                 WRITE ( cplout , '(a9,i8,F8.3)' ) ' with sc', second, &
                   &  time_config%cur_datetime%second
             ENDIF
             EXIT
          ELSE
             READ  ( unit = rest_unit, IOSTAT=eof )
             READ  ( unit = rest_unit, IOSTAT=eof )
             READ  ( unit = rest_unit, IOSTAT=eof )
             READ  ( unit = rest_unit, IOSTAT=eof )
          ENDIF
       ENDDO
       CLOSE ( unit = rest_unit )

    ENDIF

    bcast_buffer(1) = global_field_id
    bcast_buffer(2) = count

    CALL MPI_BCAST ( bcast_buffer, 2, MPI_INTEGER, ICON_root, ICON_comp_comm, ierror )

    global_field_id = bcast_buffer(1)
    count           = bcast_buffer(2)

    ! ids do not match, we have not found the correct field in the restart

    IF ( global_field_id /= fptr%global_field_id ) RETURN

    info = 1

    ! Could we use MPI_Type_vector?
    DO i = 1, field_shape(3)
       CALL MPI_SCATTERV ( global_field (:,i), fptr%local_sizes, fptr%displacement, datatype, &
            coupling_field(:,i), local_size, datatype, &
            ICON_root, ICON_comp_comm, ierror )
    ENDDO

    IF ( debug_coupler_level > 2 ) THEN
       rest_unit = find_next_free_unit(10,100)
       WRITE(file_name(1:22), '(A15,A4,A1,I2.2)') 'debug_read_rst_', &
              fptr%field_name, '_', ICON_local_rank
       OPEN ( UNIT = rest_unit, FILE = file_name(1:22), FORM = "FORMATTED", STATUS = "UNKNOWN" )
       DO j = 1, field_shape(3)
       DO i = field_shape(1), field_shape(2)
         WRITE ( rest_unit ,  * ) j, i, coupling_field(i,j)
       ENDDO
       ENDDO
       CLOSE ( UNIT = rest_unit)
    ENDIF

    IF ( debug_coupler_level > -1 ) THEN
       WRITE ( cplout , '(a,a)' ) 'cpl_read_restart: overwrite coupling field ', &
         &  TRIM(fptr%field_name)
    ENDIF

    DEALLOCATE (global_field)

    CLOSE ( unit = rest_unit )

    IF ( count > 0 ) THEN
       coupling_field(:,:) = coupling_field(:,:) / REAL(count,wp)
    ENDIF
#else
    print *, 'Restart is not supported'
    fptr => cpl_fields(field_id)
    coupling_field(field_shape(1):field_shape(2),field_shape(3)) = 0.0_wp
    info   = 0
    ierror = 0
#endif

  END SUBROUTINE cpl_read_restart

  ! ---------------------------------------------------------------------

END MODULE mo_icon_cpl_restart
