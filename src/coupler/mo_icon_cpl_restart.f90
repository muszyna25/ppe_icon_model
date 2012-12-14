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
  USE mo_util_file, ONLY      : util_symlink, util_rename, util_islink, util_unlink
  USE mo_io_units, ONLY       : find_next_free_unit
  USE mo_icon_cpl, ONLY       : t_comp, t_grid, t_cpl_field, &
       &                        comps, grids, cpl_fields,    &
       &                        cplout, debug_coupler_level
  USE mo_time_config, ONLY    : time_config

  PRIVATE
 
  INCLUDE 'cdi.inc'

  TYPE(t_grid), POINTER      :: gptr ! pointer to grid struct
  TYPE(t_cpl_field), POINTER :: fptr ! pointer to field struct

  INTEGER                :: len
  INTEGER                :: rest_unit
  INTEGER                :: year, month, day, hour, minute, second

  CHARACTER(len=132)     :: file_name
  CHARACTER(len=132)     :: link_name

  PUBLIC :: cpl_init_restart, icon_cpl_write_restart, cpl_read_restart 

CONTAINS

  ! ---------------------------------------------------------------------

  SUBROUTINE cpl_init_restart ( field_id, ierror )

    INTEGER, INTENT(IN)  :: field_id
    INTEGER, INTENT(out) :: ierror           !<  returned error code

#ifndef NOMPI
    INTEGER :: i, ii

    fptr => cpl_fields(field_id)
    gptr => grids(fptr%grid_id)

    ! local size(s) without halo

    fptr%local_size = 0

    ! This is also done in the search routine.
    ! Bad style: provide functions instead
    ! rather than duplicatiing code.

    DO i = fptr%field_shape(1), fptr%field_shape(2)
       IF ( gptr%glob_index_rank(i) == ICON_local_rank ) THEN
          fptr%local_size = fptr%local_size + 1
       ENDIF
    ENDDO

    ! Collect local sizes of coupling fields

    ALLOCATE(fptr%local_index(fptr%local_size))
    ALLOCATE(fptr%local_sizes(ICON_local_size))

    IF ( ICON_local_rank == ICON_root ) THEN
       ALLOCATE(fptr%displacement(ICON_local_size))
    ELSE
       ALLOCATE(fptr%displacement(1))
       ALLOCATE(fptr%global_index(1))
    ENDIF

    CALL MPI_ALLGATHER ( fptr%local_size, 1, MPI_INTEGER, fptr%local_sizes, 1, MPI_INTEGER, &
                         ICON_comp_comm, ierror )

    IF ( ICON_local_rank == ICON_root ) THEN

       ! Compute offsets

       fptr%displacement(1) = 0

       DO i = 2, ICON_local_size
          fptr%displacement(i) = fptr%displacement(i-1) + fptr%local_sizes(i-1)
       ENDDO

       fptr%global_size = SUM(fptr%local_sizes)

       ALLOCATE (fptr%global_index(fptr%global_size))

    ELSE

       fptr%global_size = SUM(fptr%local_sizes)

       ALLOCATE (fptr%global_index(1))

    ENDIF

    ii = 0
    DO i = fptr%field_shape(1), fptr%field_shape(2)
       IF ( gptr%glob_index_rank(i) == ICON_local_rank ) THEN
          ii = ii + 1
          fptr%local_index(ii) = gptr%grid_glob_index(i)
       ENDIF
    ENDDO

    ! Collect global indices without halos

    CALL MPI_GATHERV ( fptr%local_index, fptr%local_size, MPI_INTEGER, &
                       fptr%global_index, fptr%local_sizes, fptr%displacement, MPI_INTEGER, &
                       ICON_root, ICON_comp_comm, ierror )
#else
    fptr => cpl_fields(field_id)
    ierror = 0
#endif

  END SUBROUTINE cpl_init_restart

  ! ---------------------------------------------------------------------

  SUBROUTINE icon_cpl_write_restart ( nfields, field_ids, ierror )

    INTEGER, INTENT(in)    :: nfields            !<  number of field id
    INTEGER, INTENT(in)    :: field_ids(nfields) !<  field id
    INTEGER, INTENT(out)   :: ierror             !<  returned error code

#ifndef NOMPI
    REAL(wp), ALLOCATABLE  :: sorted_field (:,:) !< global sorted field
    REAL(wp), ALLOCATABLE  :: global_field (:,:)
    REAL(wp), ALLOCATABLE  :: buffer_field (:,:)

    INTEGER :: streamID, zaxisID, timeID, vlistID, nmiss
    INTEGER :: count
    INTEGER :: nbr_restart_fields
    INTEGER :: n, i, j, jj

    ierror = 0

    nbr_restart_fields = 0

    DO n = 1, nfields
       IF ( cpl_fields(field_ids(n))%coupling%restart_flag ) nbr_restart_fields = nbr_restart_fields + 1
    ENDDO

    IF ( nbr_restart_fields == 0 ) RETURN

    IF ( ICON_local_rank == ICON_root ) THEN

       len = 12 + LEN_TRIM(comps(1)%comp_name) + 1 + 16

       second = NINT(time_config%cur_datetime%second)

       WRITE(file_name(1:len), '(A12,A,A1,I4.4,2I2.2,A1,3I2.2,A1)')   &
                       "restart_cpl_", TRIM(comps(1)%comp_name), "_", &
                                   time_config%cur_datetime%year,     &
                                   time_config%cur_datetime%month,    &
                                   time_config%cur_datetime%day, "T", &
                                   time_config%cur_datetime%hour,     &
                                   time_config%cur_datetime%minute,   &
                                   second, "Z"

#ifndef CPL_WRITE_RAW_BINARY_RESTART

       WRITE(file_name(len+1:len+3), '(A3)' ) ".nc"

       link_name = "restart_cpl_"//TRIM(comps(1)%comp_name)//".nc"

       IF ( debug_coupler_level > 1 ) &
          WRITE(cplout,*) 'Open file for write ', TRIM(file_name)

       streamID = streamopenwrite(TRIM(file_name),FILETYPE_NC2)

       IF ( streamID < 0 ) &
            WRITE(cplout,*) 'Problem opening file ', TRIM(file_name), ': ', cdiStringError(streamID)

       vlistID = vlistCreate()
       zaxisID = zaxisCreate ( ZAXIS_SURFACE, 1 )

       DO n = 1, nfields

          IF ( .NOT. cpl_fields(field_ids(n))%coupling%restart_flag ) CYCLE 

          fptr => cpl_fields(field_ids(n))

          fptr%coupling%cdi_gridID = gridCreate ( GRID_GENERIC, fptr%global_size*fptr%field_shape(3) )

          CALL gridDefXsize ( fptr%coupling%cdi_gridID, fptr%global_size)
          CALL gridDefXname ( fptr%coupling%cdi_gridID, 'x')
          CALL gridDefYsize ( fptr%coupling%cdi_gridID, fptr%field_shape(3))
          CALL gridDefYname ( fptr%coupling%cdi_gridID, 'y')

          fptr%coupling%cdi_varID = vlistDefVar ( vlistID, fptr%coupling%cdi_gridID, zaxisID, TIME_CONSTANT )

          ! Define the variable names and precision
          CALL vlistDefVarName     ( vlistID, fptr%coupling%cdi_varID, TRIM(fptr%field_name) )
          CALL vlistDefVarDatatype ( vlistID, fptr%coupling%cdi_varID, DATATYPE_FLT64 )

          nmiss = 0

          IF ( debug_coupler_level > 1 ) THEN
             WRITE(cplout,*) "Writing file...", &
               TRIM(file_name),TRIM(fptr%field_name),vlistID, fptr%coupling%cdi_gridID, zaxisID
          ENDIF

       ENDDO

       CALL streamDefVlist ( streamID, vlistID )

       timeID = streamDefTimestep ( streamID, 0 )

#else

       link_name = "restart_cpl_"//TRIM(comps(1)%comp_name)

       OPEN ( UNIT = rest_unit, FILE = file_name(1:len), FORM = "UNFORMATTED", STATUS = "NEW" )

#endif

    ENDIF

    DO n = 1, nfields

       IF ( .NOT. cpl_fields(field_ids(n))%coupling%restart_flag ) CYCLE 

       fptr => cpl_fields(field_ids(n))

       count = fptr%accumulation_count

       ALLOCATE (buffer_field (fptr%local_size,fptr%field_shape(3)) )

       IF ( ICON_local_rank == ICON_root ) THEN
          ALLOCATE (global_field (fptr%global_size,fptr%field_shape(3)))
       ELSE
          ALLOCATE (global_field (1,fptr%field_shape(3)))
       ENDIF

       DO i = 1, fptr%field_shape(3)
          jj = 0
          DO j = fptr%field_shape(1), fptr%field_shape(2)
             IF ( gptr%glob_index_rank(j-fptr%field_shape(1)+1) == ICON_local_rank ) THEN
                jj = jj + 1
                buffer_field(jj,i) = fptr%send_field_acc(j,i)
             ENDIF
          ENDDO

          CALL MPI_GATHERV ( buffer_field(1,i), fptr%local_size, datatype,       &
               global_field(:,i), fptr%local_sizes, fptr%displacement, datatype, &
               ICON_root, ICON_comp_comm, ierror )
       ENDDO

       rest_unit = find_next_free_unit(10,100)

       IF ( debug_coupler_level > 2 ) THEN
          WRITE(file_name(1:22), '(A15,A4,A1,I2.2)') 'debug_wrte_rst_', &
               fptr%field_name, '_', ICON_local_rank
          OPEN ( UNIT = rest_unit, FILE = file_name(1:22), FORM = "FORMATTED", STATUS = "UNKNOWN" )
          DO j = 1, fptr%field_shape(3)
             DO i = 1, fptr%local_size
                WRITE ( rest_unit ,  * ) j, i, buffer_field(i,j)
             ENDDO
          ENDDO
          CLOSE ( UNIT = rest_unit)
       ENDIF

       IF ( ICON_local_rank == ICON_root ) THEN

          ALLOCATE (sorted_field (fptr%global_size,fptr%field_shape(3)))

          DO i = 1, fptr%global_size
             sorted_field(fptr%global_index(i),:) = global_field(i,:)
          ENDDO

#ifndef CPL_WRITE_RAW_BINARY_RESTART

          ! ToDo: count needs to be stored in the file

          sorted_field(:,:) = sorted_field(:,:) / REAL(count,wp)

          CALL streamWriteVar(streamID, fptr%coupling%cdi_varID, sorted_field, nmiss)

#else
          sorted_field(:,:) = sorted_field(:,:)

          WRITE  ( UNIT = rest_unit ) fptr%global_field_id
          WRITE  ( UNIT = rest_unit ) time_config%cur_datetime%year,   &
                                      time_config%cur_datetime%month,  &
                                      time_config%cur_datetime%day,    &
                                      time_config%cur_datetime%hour,   &
                                      time_config%cur_datetime%minute, &
                                      time_config%cur_datetime%second
          WRITE  ( UNIT = rest_unit ) count
          WRITE  ( UNIT = rest_unit ) fptr%global_size, fptr%field_shape(3)
          WRITE  ( UNIT = rest_unit ) sorted_field
#endif
          DEALLOCATE (sorted_field)

       ENDIF

       DEALLOCATE (global_field)
       DEALLOCATE (buffer_field)

    ENDDO

    IF ( ICON_local_rank == ICON_root ) THEN

       IF ( debug_coupler_level > 1 ) &
          WRITE(cplout,*) "Closing file...", TRIM(file_name)

#ifndef CPL_WRITE_RAW_BINARY_RESTART
       CALL streamClose(streamID)

       CALL vlistdestroy(vlistid)
       CALL zaxisdestroy(zaxisid)

       DO n = 1, nfields
          IF ( .NOT. cpl_fields(field_ids(n))%coupling%restart_flag ) CYCLE 
          fptr => cpl_fields(field_ids(n))
          CALL griddestroy(fptr%coupling%cdi_gridID)
       ENDDO
#else
       CLOSE  ( UNIT = rest_unit )
#endif

       IF (util_islink(TRIM(link_name))) ierror = util_unlink(TRIM(link_name))
       ierror = util_symlink(TRIM(file_name),TRIM(link_name))

    ENDIF

#else

    print *, 'Restart is not supported'
    ierror = 0

#endif

  END SUBROUTINE icon_cpl_write_restart

  ! ---------------------------------------------------------------------

  SUBROUTINE cpl_read_restart ( field_id, field_shape, coupling_field, info, ierror )

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of send field

    REAL (wp), INTENT(out) :: coupling_field (field_shape(1):field_shape(2),field_shape(3))

    INTEGER, INTENT(out)   :: info             !<  returned info code
    INTEGER, INTENT(out)   :: ierror           !<  returned error code

    INTEGER                :: vlistID, varID
    INTEGER                :: nvars
    CHARACTER(len=256)     :: name

#ifndef NOMPI
    REAL(wp), ALLOCATABLE  :: global_field (:,:)
    REAL(wp), ALLOCATABLE  :: local_field  (:,:)

    INTEGER :: i, j, jj
    INTEGER :: streamID, nmiss

    INTEGER :: year, month, day, hour, minute, second
    INTEGER :: global_size, nbr_bundles
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

#ifndef CPL_WRITE_RAW_BINARY_RESTART
    file_name = "restart_cpl_"//TRIM(comps(1)%comp_name)//".nc"
#else
    file_name = "restart_cpl_"//TRIM(comps(1)%comp_name)
#endif

    INQUIRE ( FILE=TRIM(file_name), EXIST = existing )

    IF ( debug_coupler_level > 0 ) &
       WRITE ( cplout , '(A20,A,L2)' ) 'Looking for restart file ', TRIM(file_name), existing

    IF ( .NOT. existing ) RETURN

    ALLOCATE (global_field (fptr%global_size,field_shape(3)))
    ALLOCATE (local_field  (fptr%local_size, field_shape(3)))

#ifndef CPL_WRITE_RAW_BINARY_RESTART

    IF ( ICON_local_rank == ICON_root ) THEN

       streamID = streamOpenRead(TRIM(file_name),FILETYPE_NC2)

       IF ( streamID < 0 ) &
          WRITE(cplout,*) 'Problem opening file ', TRIM(file_name), ': ', cdiStringError(streamID)

       IF ( debug_coupler_level > 0 ) &
          WRITE(cplout,*) 'Open restart file ', TRIM(file_name)

       vlistID = streamInqVlist(streamID)

       nvars   = vlistNvars(vlistID)

       ! Reading fields from restart file
       !
       ! This can be done in a more intelligent way rather than comparing strings
       !
       DO varID = 0, nvars - 1
          CALL vlistInqVarName ( vlistID, varID, name )
          IF ( TRIM(fptr%field_name) == TRIM(name) ) THEN
            IF ( debug_coupler_level > 0 ) &
              WRITE(cplout,*) 'Reading ', TRIM(name), ' from restart.'
            CALL streamReadVar(streamID, varID, global_field, nmiss)
          ENDIF
       ENDDO

       count = 1

       CALL streamClose(streamID)

    ENDIF
    
#else

    IF ( ICON_local_rank == ICON_root ) THEN

       rest_unit = find_next_free_unit(10,100)
       OPEN  ( UNIT = rest_unit, FILE = TRIM(file_name), FORM = "UNFORMATTED", status = 'OLD' )

       global_field_id = -1

       DO WHILE ( eof >=0 ) 
          READ ( UNIT = rest_unit, IOSTAT=eof ) global_field_id
          IF ( global_field_id == fptr%global_field_id ) THEN
             READ ( UNIT = rest_unit, IOSTAT=eof ) year, month, day, hour, minute, second
             READ ( UNIT = rest_unit, IOSTAT=eof ) count
             READ ( UNIT = rest_unit ) global_size, nbr_bundles
             READ ( UNIT = rest_unit, IOSTAT=eof ) global_field
             IF ( debug_coupler_level > 0 ) THEN
                 WRITE ( cplout , '(a)' ) ' reading restart ...'
                 WRITE ( cplout , '(a9,2i8)' ) ' with yr ', year,  time_config%cur_datetime%year
                 WRITE ( cplout , '(a9,2i8)' ) ' with mo ', month, time_config%cur_datetime%month
                 WRITE ( cplout , '(a9,2i8)' ) ' with dy ', day,   time_config%cur_datetime%day
                 WRITE ( cplout , '(a9,2i8)' ) ' with hr ', hour,  time_config%cur_datetime%hour
                 WRITE ( cplout , '(a9,2i8)' ) ' with mn ', minute, time_config%cur_datetime%minute
                 WRITE ( cplout , '(a9,i8,F8.3)' ) ' with sc ', second, &
                   &  time_config%cur_datetime%second
             ENDIF
             EXIT
          ELSE
             READ ( UNIT = rest_unit, IOSTAT=eof )
             READ ( UNIT = rest_unit, IOSTAT=eof )
             READ ( UNIT = rest_unit, IOSTAT=eof )
             READ ( UNIT = rest_unit, IOSTAT=eof )
          ENDIF
       ENDDO

       CLOSE ( UNIT = rest_unit )

    ENDIF

    ! ids do not match, we have not found the correct field in the restart

    IF ( global_field_id /= fptr%global_field_id ) THEN
       IF ( debug_coupler_level > 0 ) THEN
          WRITE ( cplout , '(a,2i8)' ) 'cpl_read_restart: ids do not match ', &
               &  global_field_id, fptr%global_field_id 
          RETURN
       ENDIF
    ENDIF

#endif

    bcast_buffer(1) = global_field_id
    bcast_buffer(2) = count

    CALL MPI_BCAST ( bcast_buffer, 2, MPI_INTEGER, &
                     ICON_root, ICON_comp_comm, ierror )

    global_field_id = bcast_buffer(1)
    count           = bcast_buffer(2)

    CALL MPI_BCAST ( global_field, fptr%global_size * field_shape(3), datatype, &
                     ICON_root, ICON_comp_comm, ierror )

    info = 1

    DO i = 1, fptr%local_size
       local_field(i,:) = global_field(fptr%local_index(i),:)
    ENDDO

    IF ( debug_coupler_level > 2 ) THEN
       rest_unit = find_next_free_unit(10,100)
       WRITE (file_name(1:22), '(A15,A4,A1,I2.2)') 'debug_read_rst_', &
              fptr%field_name, '_', ICON_local_rank
       OPEN ( UNIT = rest_unit, FILE = file_name(1:22), FORM = "FORMATTED", STATUS = "UNKNOWN" )
       DO j = 1, field_shape(3)
       DO i = 1, fptr%local_size
         WRITE ( rest_unit ,  * ) j, i, local_field(i,j)
       ENDDO
       ENDDO
       CLOSE ( UNIT = rest_unit)
    ENDIF

    IF ( debug_coupler_level > 0 ) THEN
       WRITE ( cplout , '(a,a)' ) 'cpl_read_restart: overwrite coupling field ', &
         &  TRIM(fptr%field_name)
    ENDIF

    DO i = 1, field_shape(3)
       jj = 0
       DO j = field_shape(1), field_shape(2)
          IF ( gptr%glob_index_rank(j-field_shape(1)+1) == ICON_local_rank ) THEN
             jj = jj + 1
             coupling_field(j,i) = local_field(jj,i)
          ELSE
             coupling_field(j,i) = 0.0_wp
          ENDIF
       ENDDO
    ENDDO

    IF ( count > 0 ) THEN
       coupling_field(:,:) = coupling_field(:,:) / REAL(count,wp)
    ENDIF

    DEALLOCATE (global_field, local_field)

#else

    print *, 'Restart is not supported'
    info   = 0
    ierror = 0

#endif

  END SUBROUTINE cpl_read_restart

  ! ---------------------------------------------------------------------

END MODULE mo_icon_cpl_restart
