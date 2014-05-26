!>
!! This module provides basic methods for reading 
!! a NetCDF file in parallel or sequential in a transparent way.
!!
!! Contains routines for reading data from netcdf-Files of various shape.
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!!
!! @author Daniel Reinert, DWD
!! @author Leonidas Linardakis, MPIM
!!
!!
!! @par Revision History
!! Moved to mo_util_netcd from mo_ext_data by Daniel reinert, DWD (2012-02-01)
!! Moved from mo_util_netcdf, added netcdf_read_oncells, by L. Linardakis, (2013-03-15)
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
#define define_fill_target REAL(wp), TARGET, OPTIONAL
!#define define_fill_target REAL(wp), TARGET, ALLOCATABLE, OPTIONAL
!#define define_fill_target REAL(wp), POINTER, OPTIONAL

MODULE mo_netcdf_read

  USE mo_kind
  USE mo_scatter,            ONLY: scatter_cells_2D, scatter_cells_2D_time, &
    &                              scatter_cells_3D_time, broadcast_array
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message_text, message, warning, finish, em_warn
  USE mo_impl_constants,     ONLY: success, max_char_length
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_units,           ONLY: filename_max

  USE mo_communication,      ONLY: idx_no, blk_no, exchange_data
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work, &
    &                              my_process_is_mpi_workroot
  USE mo_fortran_tools,      ONLY: assign_if_present
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: read_netcdf_data
  PUBLIC :: read_netcdf_data_single
  PUBLIC :: nf
  PUBLIC :: netcdf_open_input, netcdf_close

  PUBLIC :: netcdf_read_0D_real
  PUBLIC :: netcdf_read_1D
  PUBLIC :: netcdf_read_2D_time
  PUBLIC :: netcdf_read_3D_time
  PUBLIC :: netcdf_read_oncells_2D
  PUBLIC :: netcdf_read_oncells_2D_time
  PUBLIC :: netcdf_read_oncells_3D_time
  PUBLIC :: netcdf_read_oncells_2D_extdim
  PUBLIC :: netcdf_read_oncells_3D_extdim

  !--------------------------------------------------------
  !>
  ! use only for testing, these routines are NOT for output
  PUBLIC :: netcdf_write_oncells_2D
  PUBLIC :: netcdf_write_oncells_2D_time
  PUBLIC :: netcdf_write_oncells_3D_time
  !--------------------------------------------------------


  INTERFACE read_netcdf_data
    MODULE PROCEDURE read_netcdf_2d
    MODULE PROCEDURE read_netcdf_2d_int
    MODULE PROCEDURE read_netcdf_3d
    MODULE PROCEDURE read_netcdf_4d
    MODULE PROCEDURE read_netcdf_time
  END INTERFACE read_netcdf_data

  INTERFACE read_netcdf_data_single
    MODULE PROCEDURE read_netcdf_3d_single
  END INTERFACE read_netcdf_data_single

  INTERFACE netcdf_read_0D_real
    MODULE PROCEDURE netcdf_read_REAL_0D_fileid
  END INTERFACE netcdf_read_0D_real

  INTERFACE netcdf_read_1D
    MODULE PROCEDURE netcdf_read_REAL_1D_filename
    MODULE PROCEDURE netcdf_read_REAL_1D_fileid
  END INTERFACE netcdf_read_1D

  INTERFACE netcdf_read_2D_time
    MODULE PROCEDURE netcdf_read_REAL_2D_time_fileid
  END INTERFACE netcdf_read_2D_time

  INTERFACE netcdf_read_3D_time
    MODULE PROCEDURE netcdf_read_REAL_3D_time_fileid
  END INTERFACE netcdf_read_3D_time

  INTERFACE netcdf_read_oncells_2D
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_filename
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_fileid
  END INTERFACE netcdf_read_oncells_2D

  INTERFACE netcdf_read_oncells_2D_time
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_time_filename
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_time_fileid
  END INTERFACE netcdf_read_oncells_2D_time

  INTERFACE netcdf_read_oncells_2D_extdim
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_1extdim_filename
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_1extdim_fileid
  END INTERFACE netcdf_read_oncells_2D_extdim

  INTERFACE netcdf_read_oncells_3D_time
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_3D_time_filename
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_3D_time_fileid
  END INTERFACE netcdf_read_oncells_3D_time

  INTERFACE netcdf_read_oncells_3D_extdim
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_3D_1extdim_filename
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_3D_1extdim_fileid
  END INTERFACE netcdf_read_oncells_3D_extdim

  INTERFACE netcdf_write_oncells_2D
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_2D_filename
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_2D_fileid
  END INTERFACE netcdf_write_oncells_2D

  INTERFACE netcdf_write_oncells_2D_time
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_2D_time_filename
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_2D_time_fileid
  END INTERFACE netcdf_write_oncells_2D_time

  INTERFACE netcdf_write_oncells_3D_time
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_3D_time_filename
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_3D_time_fileid
  END INTERFACE netcdf_write_oncells_3D_time

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS

  !-------------------------------------------------------------------------
  ! used for finding the names of the dimensions in the netcdf files
  CHARACTER(LEN=*), PARAMETER :: std_cells_dim_name_1 = 'cell'
  CHARACTER(LEN=*), PARAMETER :: std_cells_dim_name_2 = 'ncells'
  CHARACTER(LEN=*), PARAMETER :: std_time_dim_name_1  = 'time'

CONTAINS

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_filename(filename, variable_name, fill_array)
    REAL(wp), POINTER :: netcdf_read_REAL_1D_filename(:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_filename'

    file_id = netcdf_open_input(filename)
    netcdf_read_REAL_1D_filename => &
      & netcdf_read_REAL_1D_fileid( &
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array)
    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_read_REAL_1D_filename
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_ONCELLS_2D_filename(filename, variable_name, fill_array, patch)
    REAL(wp), POINTER :: netcdf_read_REAL_ONCELLS_2D_filename(:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_filename'

    file_id = netcdf_open_input(filename)
    netcdf_read_REAL_ONCELLS_2D_filename => &
      & netcdf_read_REAL_ONCELLS_2D_fileid( &
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array, patch=patch)
    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_read_REAL_ONCELLS_2D_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_ONCELLS_2D_time_filename(filename, variable_name, fill_array, patch, &
    & start_timestep, end_timestep )

    REAL(wp), POINTER            :: netcdf_read_REAL_ONCELLS_2D_time_filename(:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_time_filename'

    file_id = netcdf_open_input(filename)
    netcdf_read_REAL_ONCELLS_2D_time_filename =>   &
      & netcdf_read_REAL_ONCELLS_2D_time_fileid(file_id=file_id, variable_name=variable_name, &
      & fill_array=fill_array, patch=patch, start_timestep=start_timestep, end_timestep=end_timestep)
    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_read_REAL_ONCELLS_2D_time_filename
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_ONCELLS_2D_1extdim_filename(filename, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, extdim_name )

    REAL(wp), POINTER            :: netcdf_read_REAL_ONCELLS_2D_1extdim_filename(:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_1extdim_filename'

    file_id = netcdf_open_input(filename)
    netcdf_read_REAL_ONCELLS_2D_1extdim_filename => &
      & netcdf_read_REAL_ONCELLS_2D_1extdim_fileid(file_id=file_id, variable_name=variable_name, &
      & fill_array=fill_array, patch=patch, &
      & start_extdim=start_extdim, end_extdim=end_extdim, extdim_name=extdim_name )
    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_read_REAL_ONCELLS_2D_1extdim_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_ONCELLS_3D_time_filename(filename, variable_name, fill_array, patch, &
    & start_timestep, end_timestep, levelsdim_name )

    REAL(wp), POINTER :: netcdf_read_REAL_ONCELLS_3D_time_filename(:,:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_3D_time_filename'

    file_id = netcdf_open_input(filename)
    netcdf_read_REAL_ONCELLS_3D_time_filename => &
      & netcdf_read_REAL_ONCELLS_3D_time_fileid( &
      & file_id=file_id,                         &
      & variable_name=variable_name,             &
      & fill_array=fill_array,                   &
      & patch=patch,                             &
      & levelsdim_name=levelsdim_name,           &
      & start_timestep=start_timestep,  end_timestep=end_timestep)
    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_read_REAL_ONCELLS_3D_time_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_ONCELLS_3D_1extdim_filename(filename, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, levelsdim_name, extdim_name )

    REAL(wp), POINTER :: netcdf_read_REAL_ONCELLS_3D_1extdim_filename(:,:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target             :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_3D_1extdim_filename'

    file_id = netcdf_open_input(filename)
    netcdf_read_REAL_ONCELLS_3D_1extdim_filename => netcdf_read_REAL_ONCELLS_3D_1extdim_fileid( &
      & file_id=file_id,                          &
      & variable_name=variable_name,              &
      & fill_array=fill_array,                    &
      &  patch=patch,                             &
      & start_extdim=start_extdim,  end_extdim=end_extdim, &
      & levelsdim_name=levelsdim_name,            &
      & extdim_name=extdim_name)
    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_read_REAL_ONCELLS_3D_1extdim_filename
  !-------------------------------------------------------------------------
!  FUNCTION netcdf_read_REAL_0D_filename(filename, variable_name)
!
!    REAL(wp)            :: netcdf_read_REAL_0D_filename
!
!    CHARACTER(LEN=*), INTENT(IN) :: filename
!    CHARACTER(LEN=*), INTENT(IN) :: variable_name
!
!    INTEGER                      :: file_id
!
!    file_id = netcdf_open_input(filename)
!    netcdf_read_REAL_0D_filename = netcdf_read_REAL_0D_fileid(file_id, variable_name)
!  END FUNCTION netcdf_read_REAL_0D_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_0D_fileid(file_id, variable_name)

    REAL(wp)            :: netcdf_read_REAL_0D_fileid

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status
    REAL(wp)     :: zlocal(1)

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_0D_fileid'


    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)
      CALL nf(nf_get_var_double(file_id, varid, zlocal(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(zlocal)

    netcdf_read_REAL_0D_fileid=zlocal(1)

  END FUNCTION netcdf_read_REAL_0D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_fileid(file_id, variable_name, fill_array)

    REAL(wp), POINTER            :: netcdf_read_REAL_1D_fileid(:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_1D_fileid'

    ! trivial return value.
    NULLIFY(netcdf_read_REAL_1D_fileid)

    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 ) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size(1)
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:2))

    IF (PRESENT(fill_array)) THEN
      netcdf_read_REAL_1D_fileid => fill_array
    ELSE
      ALLOCATE( netcdf_read_REAL_1D_fileid(var_size(1)), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_1D_fileid )')
      ENDIF
    ENDIF

    ! check if the size is correct
    IF (SIZE(netcdf_read_REAL_1D_fileid,1) < var_size(1)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(netcdf_read_REAL_1D_fileid,1) > var_size(1)) &
      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(file_id, varid, netcdf_read_REAL_1D_fileid(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(netcdf_read_REAL_1D_fileid)

  END FUNCTION netcdf_read_REAL_1D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_2D_time_fileid(file_id, variable_name, fill_array, dim_names, start_timestep, end_timestep)

    REAL(wp), POINTER            :: netcdf_read_REAL_2D_time_fileid(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(3), count_read_index(3)
    INTEGER :: idim
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_2D_time_fileid'

    ! trivial return value.
    NULLIFY(netcdf_read_REAL_2D_time_fileid)

    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 3 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 2
          IF (TRIM(dim_names(idim)) /= TRIM(var_dim_name(idim))) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),' but dimension name', &
                       TRIM(ADJUSTL(dim_names(idim))),' was expected'
            CALL finish(method_name, 'dimension name mismatch')
          END IF
        END DO
      END IF

      IF (TRIM(var_dim_name(3)) /= 'time' ) THEN
        WRITE(0,*) 'no time dimension found, the last dimension must be time'
        CALL finish(method_name, 'no time dimension found')
      END IF
    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:3))
    file_time_steps=var_size(3)

    ! calculate time range
    IF (PRESENT(start_timestep)) THEN
      start_time = start_timestep
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_timestep)) THEN
      end_time = end_timestep
    ELSE
      end_time = file_time_steps
    ENDIF
!!$    use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "number of time steps < 1")

    IF (PRESENT(fill_array)) THEN
      netcdf_read_REAL_2D_time_fileid => fill_array
    ELSE
      ALLOCATE( netcdf_read_REAL_2D_time_fileid(var_size(1),var_size(2),time_steps), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_2D_time_fileid )')
      ENDIF
      netcdf_read_REAL_2D_time_fileid(:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/1,1,start_time/)
      count_read_index = (/var_size(1),var_size(2),time_steps/)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, count_read_index, &
              netcdf_read_REAL_2D_time_fileid(:,:,:)), variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(netcdf_read_REAL_2D_time_fileid)

  END FUNCTION netcdf_read_REAL_2D_time_fileid
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_3D_time_fileid(file_id, variable_name, fill_array, dim_names, start_timestep, end_timestep)

    REAL(wp), POINTER            :: netcdf_read_REAL_3D_time_fileid(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(4), count_read_index(4)
    INTEGER :: idim
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_3D_time_fileid'

    ! trivial return value.
    NULLIFY(netcdf_read_REAL_3D_time_fileid)

    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 4 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 3
          IF (TRIM(dim_names(idim)) /= TRIM(var_dim_name(idim))) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),' but dimension name ', &
                       TRIM(ADJUSTL(dim_names(idim))),' was expected'
            CALL finish(method_name, 'dimension name mismatch')
          END IF
        END DO
      END IF

      IF (TRIM(var_dim_name(4)) /= 'time' ) THEN
        WRITE(0,*) 'no time dimension found, the last dimension must be time'
        CALL finish(method_name, 'no time dimension found')
      END IF
    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:4))
    file_time_steps=var_size(4)

    ! calculate time range
    IF (PRESENT(start_timestep)) THEN
      start_time = start_timestep
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_timestep)) THEN
      end_time = end_timestep
    ELSE
      end_time = file_time_steps
    ENDIF
!!$    use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "number of time steps < 1")

    IF (PRESENT(fill_array)) THEN
      netcdf_read_REAL_3D_time_fileid => fill_array
    ELSE
      ALLOCATE( netcdf_read_REAL_3D_time_fileid(var_size(1),var_size(2), &
                var_size(3),time_steps), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_3D_time_fileid )')
      ENDIF
      netcdf_read_REAL_3D_time_fileid(:,:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/1,1,1,start_time/)
      count_read_index = (/var_size(1),var_size(2),var_size(3),time_steps/)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, count_read_index, &
              netcdf_read_REAL_3D_time_fileid(:,:,:,:)), variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(netcdf_read_REAL_3D_time_fileid)

  END FUNCTION netcdf_read_REAL_3D_time_fileid
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_ONCELLS_2D_fileid(file_id, variable_name, fill_array, patch)

    REAL(wp), POINTER            :: netcdf_read_REAL_ONCELLS_2D_fileid(:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: total_number_of_cells
    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:)

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_fileid'

    ! trivial return value.
    NULLIFY(netcdf_read_REAL_ONCELLS_2D_fileid)

    total_number_of_cells = patch%n_patch_cells_g

    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 .OR. var_size(1) /= total_number_of_cells) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, &
          & " total_number_of_cells=", total_number_of_cells
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (.NOT. check_is_cell_dim_name(var_dim_name(1))) &
         CALL warning(method_name, "dim_name /= std_cells_dim_name")

    ENDIF

    ALLOCATE( tmp_array(total_number_of_cells), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(file_id, varid, tmp_array(:)), variable_name)
    ENDIF

    IF (PRESENT(fill_array)) THEN
      netcdf_read_REAL_ONCELLS_2D_fileid => fill_array
    ELSE
      ALLOCATE( netcdf_read_REAL_ONCELLS_2D_fileid(nproma, patch%alloc_cell_blocks), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_ONCELLS_2D_fileid )')
      ENDIF
      netcdf_read_REAL_ONCELLS_2D_fileid(:,:) = 0.0_wp
    ENDIF

    CALL scatter_cells_2D(tmp_array, netcdf_read_REAL_ONCELLS_2D_fileid, patch)

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_ONCELLS_2D_fileid
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_REAL_ONCELLS_2D_time_fileid(file_id, variable_name, fill_array, patch, &
    & start_timestep, end_timestep)
    REAL(wp), POINTER            :: netcdf_read_REAL_ONCELLS_2D_time_fileid(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep

    netcdf_read_REAL_ONCELLS_2D_time_fileid => netcdf_read_REAL_ONCELLS_2D_1extdim_fileid(&
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array, patch=patch, &
      & start_extdim=start_timestep, end_extdim=end_timestep, extdim_name="time" )

  END FUNCTION netcdf_read_REAL_ONCELLS_2D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_REAL_ONCELLS_2D_1extdim_fileid(file_id, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, extdim_name )

    REAL(wp), POINTER            :: netcdf_read_REAL_ONCELLS_2D_1extdim_fileid(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    INTEGER :: total_number_of_cells
    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_allocated_step, end_allocated_step
    LOGICAL :: use_time_range
    INTEGER :: start_read_index(2), count_read_index(2)

    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_1extdim_fileid'

    NULLIFY(netcdf_read_REAL_ONCELLS_2D_1extdim_fileid)
    total_number_of_cells = patch%n_patch_cells_g

    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the dimensions have reasonable names
      IF (.NOT. check_is_cell_dim_name(var_dim_name(1))) THEN
        write(0,*) var_dim_name(1)
        WRITE(message_text,*) variable_name, " ", TRIM(var_dim_name(1)), " /= std_cells_dim_name"
        CALL finish(method_name, message_text)
      ENDIF

      IF (PRESENT(extdim_name)) THEN
        IF (TRIM(extdim_name) /= TRIM(var_dim_name(2))) THEN
          WRITE(message_text,*) variable_name, ":", TRIM(extdim_name), "/=", TRIM(var_dim_name(2))
          CALL finish(method_name, message_text)
        ENDIF
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= total_number_of_cells) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", var_size, &
          & " total_number_of_cells=", total_number_of_cells
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:2))
    file_time_steps      = var_size(2)

    ! calculate time range
    IF (PRESENT(start_extdim)) THEN
      start_time = start_extdim
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_extdim)) THEN
      end_time = end_extdim
    ELSE
      end_time = file_time_steps
    ENDIF
    use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "ext dim size < 1")
    !-----------------------
    IF (PRESENT(fill_array)) THEN
      netcdf_read_REAL_ONCELLS_2D_1extdim_fileid => fill_array
    ELSE
      ALLOCATE( netcdf_read_REAL_ONCELLS_2D_1extdim_fileid(nproma, patch%alloc_cell_blocks, time_steps), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_ONCELLS_2D_1extdim_fileid )')
      ENDIF
      netcdf_read_REAL_ONCELLS_2D_1extdim_fileid(:,:,:) = 0.0_wp
    ENDIF
    start_allocated_step = LBOUND(netcdf_read_REAL_ONCELLS_2D_1extdim_fileid, 3)
    end_allocated_step   = UBOUND(netcdf_read_REAL_ONCELLS_2D_1extdim_fileid, 3)
    ALLOCATE( tmp_array(total_number_of_cells, start_allocated_step:end_allocated_step), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF
    IF (SIZE(netcdf_read_REAL_ONCELLS_2D_1extdim_fileid,3) < time_steps) &
      CALL finish(method_name, "allocated size < time_steps")
    !-----------------------

    IF( my_process_is_mpi_workroot()) THEN
  !    CALL nf(nf_get_var_double(file_id, varid, tmp_array(:,:)), variable_name)
      start_read_index = (/ 1, start_time /)
      count_read_index      = (/ total_number_of_cells, time_steps /)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, count_read_index, tmp_array(:,:)), variable_name)
    ENDIF

    CALL scatter_cells_2D_time(tmp_array, netcdf_read_REAL_ONCELLS_2D_1extdim_fileid, patch)

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_ONCELLS_2D_1extdim_fileid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION netcdf_read_REAL_ONCELLS_3D_time_fileid(file_id, variable_name, fill_array, patch, &
    & start_timestep, end_timestep, levelsdim_name)

    REAL(wp), POINTER  :: netcdf_read_REAL_ONCELLS_3D_time_fileid(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    netcdf_read_REAL_ONCELLS_3D_time_fileid => netcdf_read_REAL_ONCELLS_3D_1extdim_fileid(&
      & file_id=file_id,                        &
      & variable_name=variable_name,            &
      & fill_array=fill_array,                  &
      & patch=patch,                            &
      & start_extdim=start_timestep,  end_extdim=end_timestep, &
      & levelsdim_name=levelsdim_name,          &
      & extdim_name="time")

  END FUNCTION netcdf_read_REAL_ONCELLS_3D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION netcdf_read_REAL_ONCELLS_3D_1extdim_fileid(file_id, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, levelsdim_name, extdim_name )

    REAL(wp), POINTER  :: netcdf_read_REAL_ONCELLS_3D_1extdim_fileid(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsdim_name

    INTEGER :: total_number_of_cells
    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_vertical_levels, file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_allocated_step, end_allocated_step
 !   LOGICAL :: use_time_range
    INTEGER :: start_read_index(3), count_read_index(3)

    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:,:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_3D_1extdim_fileid'

    NULLIFY(netcdf_read_REAL_ONCELLS_3D_1extdim_fileid)

    total_number_of_cells = patch%n_patch_cells_g

    IF( my_process_is_mpi_workroot()  ) THEN
      return_status = netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 3 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the dim have reasonable names
      IF (.NOT. check_is_cell_dim_name(var_dim_name(1))) THEN
        write(0,*) var_dim_name(1)
        WRITE(message_text,*) variable_name, " ", TRIM(var_dim_name(3)), " /= std_cells_dim_name"
        CALL finish(method_name, message_text)
      ENDIF

      IF (PRESENT(levelsdim_name)) THEN
        IF (TRIM(levelsdim_name) /= TRIM(var_dim_name(2))) THEN
          WRITE(message_text,*) variable_name, ":", levelsdim_name, "/=",  var_dim_name(2)
          CALL finish(method_name, message_text)
        ENDIF
      ENDIF
      IF (PRESENT(extdim_name)) THEN
        IF (TRIM(extdim_name) /= TRIM(var_dim_name(3))) THEN
          WRITE(message_text,*) variable_name, ":", extdim_name, "/=",  var_dim_name(3)
          CALL finish(method_name, message_text)
        ENDIF
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= total_number_of_cells) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", var_size, &
          & " total_number_of_cells=", total_number_of_cells
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:3))
    file_vertical_levels = var_size(2)
    file_time_steps      = var_size(3)

    ! calculate time range
    IF (PRESENT(start_extdim)) THEN
      start_time = start_extdim
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_extdim)) THEN
      end_time = end_extdim
    ELSE
      end_time = file_time_steps
    ENDIF
  !  use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
  !  write(0,*) "start,end time=", start_time, end_time
    time_steps = end_time - start_time + 1
    IF (time_steps < 1) &
      & CALL finish(method_name, "extdim size < 1")
    !-----------------------
    !-----------------------
    IF (PRESENT(fill_array)) THEN
      netcdf_read_REAL_ONCELLS_3D_1extdim_fileid => fill_array
    ELSE
      ALLOCATE( netcdf_read_REAL_ONCELLS_3D_1extdim_fileid &
        & (nproma, file_vertical_levels, patch%alloc_cell_blocks, time_steps), &
        & stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_ONCELLS_3D_1extdim_fileid )')
      ENDIF
      netcdf_read_REAL_ONCELLS_3D_1extdim_fileid(:,:,:,:) = 0.0_wp
    ENDIF
    start_allocated_step = LBOUND(netcdf_read_REAL_ONCELLS_3D_1extdim_fileid, 4)
    end_allocated_step   = UBOUND(netcdf_read_REAL_ONCELLS_3D_1extdim_fileid, 4)
    ALLOCATE( tmp_array(total_number_of_cells, file_vertical_levels, start_allocated_step:end_allocated_step), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF (SIZE(netcdf_read_REAL_ONCELLS_3D_1extdim_fileid,4) < time_steps) &
      CALL finish(method_name, "allocated size < time_steps")

    IF (file_vertical_levels /= SIZE(netcdf_read_REAL_ONCELLS_3D_1extdim_fileid,2)) &
      CALL finish(method_name, 'file_vertical_levels /= SIZE(fill_array,2)')
    !-----------------------

    IF( my_process_is_mpi_workroot()) THEN
  !    CALL nf(nf_get_var_double(file_id, varid, tmp_array(:,:,:)), variable_name)
      start_read_index = (/ 1, 1, start_time /)
      count_read_index      = (/ total_number_of_cells, file_vertical_levels, time_steps /)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, count_read_index, tmp_array(:,:,:)), variable_name)
    ENDIF


    CALL scatter_cells_3D_time(tmp_array, netcdf_read_REAL_ONCELLS_3D_1extdim_fileid, patch)

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_ONCELLS_3D_1extdim_fileid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_open_output(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: file_id, old_mode

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL nf(nf_set_default_format(nf_format_64bit, old_mode), TRIM(filename))
      CALL nf(nf_create(TRIM(filename), nf_clobber, file_id), TRIM(filename))
      CALL nf(nf_set_fill(file_id, nf_nofill, old_mode), TRIM(filename))
    ELSE
        file_id = -1 ! set it to an invalid value
    ENDIF

    netcdf_open_output = file_id

  END FUNCTION netcdf_open_output
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_open_input(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: file_id

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL nf(nf_open(TRIM(filename), nf_nowrite, file_id), TRIM(filename))
    ELSE
      file_id = -1 ! set it to an invalid value
    ENDIF

    netcdf_open_input = file_id
    
  END FUNCTION netcdf_open_input
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_close(file_id)
    INTEGER, INTENT(IN) :: file_id

    netcdf_close = -1
    IF( my_process_is_mpi_workroot()  ) THEN
        netcdf_close = nf_close(file_id)
    ENDIF

  END FUNCTION netcdf_close
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION check_is_cell_dim_name(dim_name)
    CHARACTER(LEN=*), INTENT(IN) :: dim_name

    SELECT CASE (TRIM(dim_name))
      CASE (std_cells_dim_name_1, std_cells_dim_name_2)
        check_is_cell_dim_name = .TRUE.
      CASE default
        check_is_cell_dim_name = .FALSE.
    END SELECT

  END FUNCTION check_is_cell_dim_name
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION check_is_time_dim_name(dim_name)
    CHARACTER(LEN=*), INTENT(IN) :: dim_name

    SELECT CASE (TRIM(dim_name))
      CASE (std_time_dim_name_1)
        check_is_time_dim_name = .TRUE.
      CASE default
        check_is_time_dim_name = .FALSE.
    END SELECT

  END FUNCTION check_is_time_dim_name
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_inq_var(file_id, name, varid, var_type, var_dims, var_size, var_dim_name)
    INTEGER, INTENT(IN) :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: name
    
    INTEGER, INTENT(OUT) :: varid, var_type, var_dims
    INTEGER, INTENT(OUT) :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max), INTENT(OUT) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER  :: number_of_attributes  
    INTEGER :: var_dims_reference(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: check_var_name
    INTEGER :: i

    netcdf_inq_var = -1
    IF ( .NOT. my_process_is_mpi_workroot() ) RETURN


    netcdf_inq_var = nf_inq_varid(file_id, name, varid)
    CALL nf(netcdf_inq_var, name)
    netcdf_inq_var = nf_inq_var (file_id, varid, check_var_name, var_type, var_dims, &
      & var_dims_reference, number_of_attributes)
    CALL nf(netcdf_inq_var, check_var_name)
    DO i=1, var_dims
      CALL nf(nf_inq_dimlen (file_id, var_dims_reference(i), var_size(i)),     check_var_name)
      CALL nf(nf_inq_dimname(file_id, var_dims_reference(i), var_dim_name(i)), check_var_name)
    ENDDO
!     write(0,*) " Read var_dims, var_size:",  var_dims, var_size
!     write(0,*) " check_var_name:",  check_var_name
!     write(0,*) " name:", name    

  END FUNCTION netcdf_inq_var
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d (file_id, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: file_id          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_2d'

    INTEGER :: varid, mpi_comm, j, jl, jb
    REAL(wp):: z_dummy_array(glb_arr_len)!< local dummy array

    ! Get var ID
    IF( my_process_is_stdio()) CALL nf(nf_inq_varid(file_id, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(file_id, varid, z_dummy_array(:)), routine)
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0._wp

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_2d
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d_int (file_id, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: file_id          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    INTEGER, INTENT(INOUT) :: &          !< output field
      &  var_out(:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_2d_int'

    INTEGER :: varid, mpi_comm, j, jl, jb
    INTEGER :: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF( my_process_is_stdio()) CALL nf(nf_inq_varid(file_id, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF( my_process_is_stdio()) CALL nf(nf_get_var_int(file_id, varid, z_dummy_array(:)), routine)
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_2d_int
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read 3D (inlcuding height) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 3 D by Guenther Zaengl, DWD (2011-07-11)
  !!
  SUBROUTINE read_netcdf_3d (file_id, varname, glb_arr_len, loc_arr_len, glb_index, &
    &                        nlevs, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: file_id          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_3d'

    INTEGER :: varid, mpi_comm, j, jl, jb, jk
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs)!< local dummy array

    !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(file_id, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(file_id, varid, z_dummy_array(:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:,:) = 0._wp

    ! Set var_out from global data
     DO jk = 1, nlevs
       DO j = 1, loc_arr_len

         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch

         var_out(jl,jk,jb) = z_dummy_array(glb_index(j),jk)

       ENDDO
     ENDDO

  END SUBROUTINE read_netcdf_3d
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read 3D dataset from netcdf file in SINGLE PRECISION
  !!
  !! @par Revision History
  !! Initial revision by F. Prill, DWD (2012-02-15)
  !! Optional switch to read 3D field in 2D slices: F. Prill, DWD (2012-12-19)
  !!
  SUBROUTINE read_netcdf_3d_single (file_id, varname, glb_arr_len, loc_arr_len, glb_index, &
    &                               nlevs, var_out, opt_lvalue_add)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: file_id          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global
    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lvalue_add !< If .TRUE., add values to given field

    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_3d_single'
    ! enable this flag to use a 3d buffer (which may be faster)
    LOGICAL, PARAMETER :: luse3dbuffer = .FALSE.
    ! time level (fixed)
    INTEGER, PARAMETER :: itime = 1

    ! local variables:
    INTEGER :: varid, mpi_comm, j, jl, jb, jk, &
      &        istart(3), icount(3), ierrstat, &
      &        dimlen(3), dims(3)
    ! SINGLE PRECISION local array
    REAL(sp), ALLOCATABLE:: tmp_buf(:,:)
    LOGICAL :: lvalue_add

    !-------------------------------------------------------------------------

    lvalue_add = .FALSE.
    CALL assign_if_present(lvalue_add, opt_lvalue_add)

    ! allocate temporary buffer:
    IF (luse3dbuffer) THEN
      ! allocate a buffer for all levels
      ALLOCATE(tmp_buf(glb_arr_len,nlevs), STAT=ierrstat)
    ELSE
      ! allocate a buffer for one vertical level
      ALLOCATE(tmp_buf(glb_arr_len,1), STAT=ierrstat)
    END IF
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! Get var ID
    IF(my_process_is_stdio()) THEN
      CALL nf(nf_inq_varid(file_id, TRIM(varname), varid), routine)

      ! Check variable dimensions:
      CALL nf(NF_INQ_VARDIMID(file_id, varid, dims(:)), routine)
      DO j=1,3
        CALL nf(NF_INQ_DIMLEN  (file_id, dims(j), dimlen(j)), routine)
      END DO
      IF ((dimlen(1) /= glb_arr_len) .OR.  &
        & (dimlen(2) /= nlevs)) THEN
        CALL finish(routine, "Incompatible dimensions!")
      END IF
    END IF

    ! initialize output field:
    var_out(:,:,:) = 0._wp

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data
    IF (luse3dbuffer) THEN
      !-- 3D buffer implementation

      IF(my_process_is_stdio()) THEN
        CALL nf(nf_get_var_real(file_id, varid, tmp_buf(:,:)), routine)
      END IF
      ! broadcast data:
      CALL p_bcast(tmp_buf, p_io, mpi_comm)
      ! Set var_out from global data
      IF (lvalue_add) THEN
        DO jk = 1, nlevs
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = var_out(jl,jk,jb) + REAL(tmp_buf(glb_index(j),jk), wp)
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlevs
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = REAL(tmp_buf(glb_index(j),jk), wp)
          ENDDO
        ENDDO
      END IF

    ELSE
      !-- 2D buffer implementation

      icount = (/ glb_arr_len,1,1 /)
      DO jk=1,nlevs
        istart = (/ 1,jk,itime /)
        IF(my_process_is_stdio()) THEN
          CALL nf(nf_get_vara_real(file_id, varid, &
            &     istart, icount, tmp_buf(:,:)), routine)
        END IF

        ! broadcast data:
        CALL p_bcast(tmp_buf, p_io, mpi_comm)
        ! Set var_out from global data
        IF (lvalue_add) THEN
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = var_out(jl,jk,jb) + REAL(tmp_buf(glb_index(j),1), wp)
          ENDDO
        ELSE
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = REAL(tmp_buf(glb_index(j),1), wp)
          ENDDO
        END IF
      END DO ! jk=1,nlevs

    END IF

    ! clean up
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_netcdf_3d_single


  !-------------------------------------------------------------------------
  !>
  !! Read 4D (inlcuding height and time) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 4 D by Kristina Froehlich, MPI-M (2011-06-16)
  !!
  SUBROUTINE read_netcdf_4d (file_id, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, &
       &                     nlevs, ntime,      var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: file_id          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: ntime         !< time levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:,:)                !< dimensions: nproma, nlevs, nblks, ntime

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_4d'

    INTEGER :: varid, mpi_comm, j, jl, jb, jk, jt
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs,ntime)!< local dummy array

    !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(file_id, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    write(0,*) ' ncep set ',varname,': begin of read - whole time array'
    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(file_id, varid, z_dummy_array(:,:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:,:) = 0._wp

    ! Set var_out from global data
    DO jt = 1, ntime
        DO jk = 1, nlevs
           DO j = 1, loc_arr_len

             jb = blk_no(j) ! Block index in distributed patch
             jl = idx_no(j) ! Line  index in distributed patch

             var_out(jl,jk,jb,jt) = z_dummy_array(glb_index(j),jk,jt)

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_4d


  !-------------------------------------------------------------------------
  !>
  !! Read 3D (inlcuding a period of time) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for time periods by Stephan Lorenz, MPI-M (2012-02-22)
  !!
  SUBROUTINE read_netcdf_time (file_id, varname, glb_arr_len, &
       &                       loc_arr_len, glb_index,     &
       &                       ntime, nstart, ncount, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: file_id          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global
    INTEGER, INTENT(IN) :: ntime         !< number of time steps to read
    INTEGER, INTENT(IN) :: nstart(2)     !< start value for reading in all array dims
    INTEGER, INTENT(IN) :: ncount(2)     !< count value for length of array to read

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)                  !< dimensions: nproma, nblks, ntime

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_time'

    INTEGER :: varid, mpi_comm, j, jl, jb, jt
    REAL(wp):: z_dummy_array(glb_arr_len,ntime)!< local dummy array

    !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(file_id, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data
    z_dummy_array(:,:) = 0.0_wp

    !write(0,*) ' Dimensions: glb, ntime: ',glb_arr_len,ntime
    !write(0,*) ' ncep set ',varname,': begin of read - time period'
    !write(0,*) ' nstart, ncount: ',nstart, ncount

    IF(my_process_is_stdio()) CALL nf(nf_get_vara_double(file_id, varid, &
      &                               nstart(:), ncount(:), z_dummy_array(:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:) = 0.0_wp

    ! Set var_out from global data
    DO jt = 1, ntime
      DO j = 1, loc_arr_len

        jb = blk_no(j) ! Block index in distributed patch
        jl = idx_no(j) ! Line  index in distributed patch

        var_out(jl,jb,jt) = z_dummy_array(glb_index(j),jt)

      ENDDO
    ENDDO

    !write(0,*) ' READ_NETCD_TIME: z_dummy_array stress-x, index 4*64+1,5:'
    !do jt=1,3
    !  write(0,*) 'jt=',jt,' val:',(z_dummy_array(j+4*64,jt),j=1,5)
    !enddo

  END SUBROUTINE read_netcdf_time


  !-------------------------------------------------------------------------
  ! these functions is meant only for checking the read methods
  ! Do not use for regular output
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_2D_filename(filename, variable_name, write_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: write_array(:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_write_REAL_ONCELLS_2D_filename'

    netcdf_write_REAL_ONCELLS_2D_filename = 0

    file_id = netcdf_open_output(filename)

    netcdf_write_REAL_ONCELLS_2D_filename = &
      & netcdf_write_REAL_ONCELLS_2D_fileid(file_id, variable_name, write_array, patch)

    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_write_REAL_ONCELLS_2D_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_2D_time_filename(filename, variable_name, write_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: write_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_write_REAL_ONCELLS_3D_time_filename'

    netcdf_write_REAL_ONCELLS_2D_time_filename = 0

    file_id = netcdf_open_output(filename)

    netcdf_write_REAL_ONCELLS_2D_time_filename = &
      & netcdf_write_REAL_ONCELLS_2D_time_fileid(file_id, variable_name, write_array, patch)

    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_write_REAL_ONCELLS_2D_time_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_3D_time_filename(filename, variable_name, write_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: write_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_write_REAL_ONCELLS_3D_time_filename'

    netcdf_write_REAL_ONCELLS_3D_time_filename = 0

    file_id = netcdf_open_output(filename)

    netcdf_write_REAL_ONCELLS_3D_time_filename = &
      & netcdf_write_REAL_ONCELLS_3D_time_fileid(file_id, variable_name, write_array, patch)

    return_status = netcdf_close(file_id)

  END FUNCTION netcdf_write_REAL_ONCELLS_3D_time_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! this function is meant only for checking the read methods
  ! Do not use for regular output
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_2D_fileid(file_id, variable_name, write_array, patch)
    INTEGER, INTENT(IN)  :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp),      POINTER       :: write_array(:,:)
    TYPE(t_patch), TARGET        :: patch

    REAL(wp), ALLOCATABLE        :: output_array(:)
    INTEGER :: total_number_of_cells
!    INTEGER :: array_time_steps, array_vertical_levels
    INTEGER :: dim_number_of_cells
!    INTEGER :: dim_array_time_steps, dim_vertical_levels
    INTEGER :: output_shape(1), dim_write_shape(1), diff_shape(1)
    INTEGER :: variable_id

!    INTEGER                      :: i,j,t
    CHARACTER(LEN=*), PARAMETER  :: method_name = 'mo_netcdf_read:netcdf_write_REAL_ONCELLS_3D_time_fileid'

    netcdf_write_REAL_ONCELLS_2D_fileid = 0

    ALLOCATE(output_array(MERGE(patch%n_patch_cells_g, 0, &
      &                         my_process_is_mpi_workroot())))
    CALL exchange_data(write_array, output_array, patch%comm_pat_gather_c)

    !----------------------------------------------------------------------
    ! Write only from mpi_workroot
    IF( my_process_is_mpi_workroot()  ) THEN

      ! Write Dimensions
      total_number_of_cells = patch%n_patch_cells_g
      output_shape          = (/ total_number_of_cells /)
      diff_shape            = (SHAPE(output_array) - output_shape)

      IF ( MAXVAL(ABS( diff_shape)) /= 0 ) THEN
        WRITE(0,*) " gather array shape:", SHAPE(output_array)
        WRITE(0,*) " computed array shape:", output_shape
        CALL finish(method_name, "gather array shape is nor correct")
      ENDIF

      CALL nf(nf_def_dim(file_id, 'ncells',  total_number_of_cells,  dim_number_of_cells),  variable_name)
      dim_write_shape = (/ dim_number_of_cells /)

      ! define variable
      ! WRITE(0,*) " define variable:", variable_name
!      CALL nf(nf_def_var(file_id, variable_name, nf_double, 1, dim_write_shape,&
!        & variable_id), variable_name)
      CALL nf(nf_def_var(file_id, variable_name, nf_float, 1, dim_write_shape,&
        & variable_id), variable_name)

      CALL nf(nf_enddef(file_id), variable_name)

!      DO t=1, array_time_steps
      ! write array
      ! WRITE(0,*) " write array..."
      CALL nf(nf_put_var_double(file_id, variable_id, output_array(:)), variable_name)
      ! CALL nf(nf_put_var_real(file_id, variable_id, REAL(output_array(:,:,:))), variable_name)
    ENDIF

    DEALLOCATE(output_array)

  END FUNCTION netcdf_write_REAL_ONCELLS_2D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! this function is meant only for checking the read methods
  ! Do not use for regular output
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_2D_time_fileid(file_id, variable_name, write_array, patch)
    INTEGER, INTENT(IN)  :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp),      POINTER       :: write_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch

    REAL(wp), ALLOCATABLE        :: output_array(:,:)
    INTEGER :: total_number_of_cells, array_time_steps
    INTEGER :: dim_number_of_cells, dim_array_time_steps
    INTEGER :: output_shape(2), dim_write_shape(2), diff_shape(2)
    INTEGER :: variable_id
    INTEGER :: start_timestep, end_timestep, timestep

!    INTEGER                      :: i,j,t
    CHARACTER(LEN=*), PARAMETER  :: method_name = 'mo_netcdf_read:netcdf_write_REAL_ONCELLS_3D_time_fileid'

    netcdf_write_REAL_ONCELLS_2D_time_fileid = 0

    start_timestep = LBOUND(write_array, 3)
    end_timestep   = UBOUND(write_array, 3)

    ALLOCATE(output_array(MERGE(patch%n_patch_cells_g, 0, &
      &                         my_process_is_mpi_workroot()), &
      &                   start_timestep:end_timestep))
    DO timestep = start_timestep, end_timestep
      CALL exchange_data(write_array(:,:,timestep), &
        &                output_array(:,timestep), &
        &                patch%comm_pat_gather_c)
    END DO

    !----------------------------------------------------------------------
    ! Write only from mpi_workroot
    IF( my_process_is_mpi_workroot()  ) THEN

      ! Write Dimensions
      total_number_of_cells = patch%n_patch_cells_g
      array_time_steps      = SIZE(write_array, 3)
      output_shape          = (/ total_number_of_cells, array_time_steps /)
      diff_shape            = (SHAPE(output_array) - output_shape)

      IF ( MAXVAL(ABS( diff_shape)) /= 0 ) THEN
        WRITE(0,*) " gather array shape:", SHAPE(output_array)
        WRITE(0,*) " computed array shape:", output_shape
        CALL finish(method_name, "gather array shape is nor correct")
      ENDIF

      CALL nf(nf_def_dim(file_id, 'ncells',  total_number_of_cells,  dim_number_of_cells),  variable_name)
      CALL nf(nf_def_dim(file_id, 'time',    array_time_steps,       dim_array_time_steps), variable_name)
      dim_write_shape = (/ dim_number_of_cells, dim_array_time_steps /)

      ! define variable
      ! WRITE(0,*) " define variable:", variable_name
!      CALL nf(nf_def_var(file_id, variable_name, nf_double, 2, dim_write_shape,&
!        & variable_id), variable_name)
      CALL nf(nf_def_var(file_id, variable_name, nf_float, 2, dim_write_shape,&
        & variable_id), variable_name)

      CALL nf(nf_enddef(file_id), variable_name)

      ! WRITE(0,*) " write array..."
      CALL nf(nf_put_var_double(file_id, variable_id, output_array(:,:)), variable_name)
      ! CALL nf(nf_put_var_real(file_id, variable_id, REAL(output_array(:,:,:))), variable_name)

    ENDIF

    DEALLOCATE(output_array)

  END FUNCTION netcdf_write_REAL_ONCELLS_2D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! this function is meant only for checking the read methods
  ! Do not use for regular output
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_3D_time_fileid(file_id, variable_name, write_array, patch)
    INTEGER, INTENT(IN)  :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp),      POINTER       :: write_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch

    REAL(wp), ALLOCATABLE        :: output_array(:,:,:)
    INTEGER :: total_number_of_cells, array_vertical_levels, array_time_steps
    INTEGER :: dim_number_of_cells, dim_vertical_levels, dim_array_time_steps
    INTEGER :: output_shape(3), dim_write_shape(3), diff_shape(3)
    INTEGER :: variable_id
    INTEGER :: start_timestep, end_timestep, timestep, nlev

!    INTEGER                      :: i,j,t
    CHARACTER(LEN=*), PARAMETER  :: method_name = 'mo_netcdf_read:netcdf_write_REAL_ONCELLS_3D_time_fileid'

    netcdf_write_REAL_ONCELLS_3D_time_fileid = 0

    nlev           = SIZE(write_array, 2)
    start_timestep = LBOUND(write_array, 4)
    end_timestep   = UBOUND(write_array, 4)

    ALLOCATE(output_array(MERGE(patch%n_patch_cells_g, 0, &
      &                         my_process_is_mpi_workroot()), &
      &                   nlev, start_timestep:end_timestep))
    DO timestep = start_timestep, end_timestep
      CALL exchange_data(write_array(:,:,:, timestep), &
        &                output_array(:, :, timestep), &
        &                patch%comm_pat_gather_c)
    END DO

    !----------------------------------------------------------------------
    ! Write only from mpi_workroot
    IF( my_process_is_mpi_workroot()  ) THEN

      ! Write Dimensions
      total_number_of_cells = patch%n_patch_cells_g
      array_vertical_levels = SIZE(write_array, 2)
      array_time_steps      = SIZE(write_array, 4)
      output_shape          = (/ total_number_of_cells, array_vertical_levels, array_time_steps /)
      diff_shape            = (SHAPE(output_array) - output_shape)

      IF ( MAXVAL(ABS( diff_shape)) /= 0 ) THEN
        WRITE(0,*) " gather array shape:", SHAPE(output_array)
        WRITE(0,*) " computed array shape:", output_shape
        CALL finish(method_name, "gather array shape is nor correct")
      ENDIF

      CALL nf(nf_def_dim(file_id, 'ncells',  total_number_of_cells,  dim_number_of_cells),  variable_name)
      CALL nf(nf_def_dim(file_id, 'plev',    array_vertical_levels,  dim_vertical_levels),  variable_name)
      CALL nf(nf_def_dim(file_id, 'time',    array_time_steps,       dim_array_time_steps), variable_name)
      dim_write_shape = (/ dim_number_of_cells, dim_vertical_levels,  dim_array_time_steps /)

      ! define variable
      ! WRITE(0,*) " define variable:", variable_name
!      CALL nf(nf_def_var(file_id, variable_name, nf_double, 3, dim_write_shape,&
!        & variable_id), variable_name)
      CALL nf(nf_def_var(file_id, variable_name, nf_float, 3, dim_write_shape,&
        & variable_id), variable_name)

      CALL nf(nf_enddef(file_id), variable_name)

!      DO t=1, array_time_steps
!        DO i=1, total_number_of_cells
!          DO j=1, array_vertical_levels
!            write(0,*) t,i,j, ":", output_array(i,j,t)
!          ENDDO
!        ENDDO
!      ENDDO
      ! write array
      ! WRITE(0,*) " write array..."
      CALL nf(nf_put_var_double(file_id, variable_id, output_array(:,:,:)), variable_name)
      ! CALL nf(nf_put_var_real(file_id, variable_id, REAL(output_array(:,:,:))), variable_name)

    ENDIF

    DEALLOCATE(output_array)

  END FUNCTION netcdf_write_REAL_ONCELLS_3D_time_fileid
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  SUBROUTINE nf(STATUS, routine, warnonly, silent)

    INTEGER, INTENT(in)           :: STATUS
    CHARACTER(len=*), INTENT(in) :: routine
    LOGICAL, INTENT(in), OPTIONAL :: warnonly
    LOGICAL, INTENT(in), OPTIONAL :: silent

    LOGICAL :: lwarnonly, lsilent

    lwarnonly = .FALSE.
    lsilent   = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF(PRESENT(silent))   lsilent   = silent

    IF (lsilent) RETURN


    IF (STATUS /= nf_noerr) THEN
      IF (lwarnonly) THEN
        CALL message( TRIM(routine)//' netCDF error', nf_strerror(STATUS), &
          & level=em_warn)
      ELSE
        CALL finish( TRIM(routine)//' netCDF error', nf_strerror(STATUS))
      ENDIF
    ENDIF

  END SUBROUTINE nf

END MODULE mo_netcdf_read
