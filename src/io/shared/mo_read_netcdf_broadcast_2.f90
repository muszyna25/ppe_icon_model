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
#define define_fill_target_int INTEGER, TARGET, OPTIONAL
!#define define_fill_target_int INTEGER, TARGET, ALLOCATABLE, OPTIONAL
!#define define_fill_target_int INTEGER, POINTER, OPTIONAL

MODULE mo_read_netcdf_broadcast_2

  USE mo_kind
  USE mo_scatter,            ONLY: scatter_array, scatter_time_array, &
    &                              broadcast_array
  USE mo_exception,          ONLY: message_text, message, warning, finish, &
    &                              em_warn
  USE mo_impl_constants,     ONLY: success
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_units,           ONLY: filename_max

  USE mo_mpi,                ONLY: my_process_is_mpi_workroot
  USE mo_read_netcdf_distributed, ONLY: var_data_1d_int, &
    &                                   var_data_2d_wp, var_data_2d_int, &
    &                                   var_data_3d_wp, var_data_3d_int
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: nf
  PUBLIC :: netcdf_open_input, netcdf_close

  PUBLIC :: netcdf_read_0D_real
  PUBLIC :: netcdf_read_1D
  PUBLIC :: netcdf_read_1D_extdim_time
  PUBLIC :: netcdf_read_1D_extdim_extdim_time
  PUBLIC :: netcdf_read_2D_int
  PUBLIC :: netcdf_read_2D
  PUBLIC :: netcdf_read_2D_time
  PUBLIC :: netcdf_read_3D
  PUBLIC :: netcdf_read_3D_time
  PUBLIC :: netcdf_read_2D_extdim
  PUBLIC :: netcdf_read_2D_extdim_int
  PUBLIC :: netcdf_read_3D_extdim

  INTERFACE netcdf_read_0D_real
    MODULE PROCEDURE netcdf_read_REAL_0D_fileid
  END INTERFACE netcdf_read_0D_real

  INTERFACE netcdf_read_1D
    MODULE PROCEDURE netcdf_read_REAL_1D_fileid
  END INTERFACE netcdf_read_1D

  INTERFACE netcdf_read_1D_extdim_time
    MODULE PROCEDURE netcdf_read_REAL_1D_extdim_time_fileid
  END INTERFACE netcdf_read_1D_extdim_time

  INTERFACE netcdf_read_1D_extdim_extdim_time
    MODULE PROCEDURE netcdf_read_REAL_1D_extdim_extdim_time_fileid
  END INTERFACE netcdf_read_1D_extdim_extdim_time

  INTERFACE netcdf_read_2D_int
    MODULE PROCEDURE netcdf_read_INT_2D_fileid
    MODULE PROCEDURE netcdf_read_INT_2D_multivar_fileid
  END INTERFACE netcdf_read_2D_int

  INTERFACE netcdf_read_2D
    MODULE PROCEDURE netcdf_read_REAL_2D_fileid
    MODULE PROCEDURE netcdf_read_REAL_2D_multivar_fileid
  END INTERFACE netcdf_read_2D

  INTERFACE netcdf_read_2D_time
    MODULE PROCEDURE netcdf_read_REAL_2D_time_fileid
  END INTERFACE netcdf_read_2D_time

  INTERFACE netcdf_read_2D_extdim
    MODULE PROCEDURE netcdf_read_REAL_2D_extdim_fileid
    MODULE PROCEDURE netcdf_read_REAL_2D_extdim_multivar_fileid
  END INTERFACE netcdf_read_2D_extdim

  INTERFACE netcdf_read_2D_extdim_int
    MODULE PROCEDURE netcdf_read_INT_2D_extdim_fileid
    MODULE PROCEDURE netcdf_read_INT_2D_extdim_multivar_fileid
  END INTERFACE netcdf_read_2D_extdim_int

  INTERFACE netcdf_read_3D
    MODULE PROCEDURE netcdf_read_REAL_3D_fileid
  END INTERFACE netcdf_read_3D

  INTERFACE netcdf_read_3D_time
    MODULE PROCEDURE netcdf_read_REAL_3D_time_fileid
  END INTERFACE netcdf_read_3D_time

  INTERFACE netcdf_read_3D_extdim
    MODULE PROCEDURE netcdf_read_REAL_3D_extdim_fileid
  END INTERFACE netcdf_read_3D_extdim

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS

CONTAINS

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_0D_fileid(file_id, variable_name) result(res)

    REAL(wp)            :: res

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    REAL(wp)     :: zlocal(1)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_0D_fileid'


    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      CALL nf(nf_get_var_double(file_id, varid, zlocal(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(zlocal)

    res=zlocal(1)

  END FUNCTION netcdf_read_REAL_0D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_fileid(file_id, variable_name, fill_array) &
    result(res)

    REAL(wp), POINTER            :: res(:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_1D_fileid'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 ) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size(1)
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:2))

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1)), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_1D_fileid )')
      ENDIF
    ENDIF

    ! check if the size is correct
    IF (SIZE(res,1) < var_size(1)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,1) > var_size(1)) &
      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(file_id, varid, res(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(res)

  END FUNCTION netcdf_read_REAL_1D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_extdim_time_fileid(file_id, variable_name, &
    &                                             fill_array, dim_names, &
    &                                             start_timestep, end_timestep) &
    result(res)

    REAL(wp), POINTER            :: res(:,:,:)

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

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_1D_extdim_time_fileid'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 3 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 2
          IF (TRIM(dim_names(idim)) /= TRIM(var_dim_name(idim))) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),&
                       ' but dimension name', &
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
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1),var_size(2),time_steps), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( res )')
      ENDIF
      res(:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/1,1,start_time/)
      count_read_index = (/var_size(1),var_size(2),time_steps/)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, &
        &                        count_read_index, res(:,:,:)), &
        &     variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(res)

  END FUNCTION netcdf_read_REAL_1D_extdim_time_fileid
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_extdim_extdim_time_fileid( &
    file_id, variable_name, fill_array, dim_names, start_timestep, &
    end_timestep) result(res)

    REAL(wp), POINTER            :: res(:,:,:,:)

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

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_1D_extdim_extdim_time_fileid'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 4 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 3
          IF (TRIM(dim_names(idim)) /= TRIM(var_dim_name(idim))) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),&
                       ' but dimension name ', &
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
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1),var_size(2), &
                var_size(3),time_steps), stat=return_status )
      IF (return_status /= success) &
        CALL finish (method_name, &
          &          'ALLOCATE( netcdf_read_REAL_1D_extdim_extdim_time_fileid )')
      res(:,:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D_fileid,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/1,1,1,start_time/)
      count_read_index = (/var_size(1),var_size(2),var_size(3),time_steps/)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, &
        &                        count_read_index, res(:,:,:,:)), &
        &     variable_name)
    ENDIF

    ! broadcast...
    CALL broadcast_array(res)

  END FUNCTION netcdf_read_REAL_1D_extdim_extdim_time_fileid
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_INT_2D_fileid(file_id, variable_name, fill_array, &
    &                                n_g, glb_index) result(res)
    INTEGER, POINTER             :: res(:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target_int       :: fill_array(:,:)
    INTEGER, INTENT(IN)          :: n_g
    INTEGER, TARGET, INTENT(IN)  :: glb_index(:)

    TYPE(var_data_2d_int) :: fill_arrays(1)
    TYPE(var_data_1d_int) :: glb_index_(1)
    TYPE(var_data_2d_int) :: results(1)

    glb_index_(1)%data => glb_index

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%data => fill_array
      results = netcdf_read_INT_2D_multivar_fileid(file_id=file_id, &
        &                                          variable_name=variable_name,&
        &                                          n_vars=1, &
        &                                          fill_arrays=fill_arrays, &
        &                                          n_g=n_g, &
        &                                          glb_index=glb_index_)
    ELSE
      results = netcdf_read_INT_2D_multivar_fileid(file_id=file_id, &
        &                                          variable_name=variable_name,&
        &                                          n_vars=1, n_g=n_g, &
        &                                          glb_index=glb_index_)
    END IF

    res => results(1)%data

  END FUNCTION netcdf_read_INT_2D_fileid

  FUNCTION netcdf_read_INT_2D_multivar_fileid(file_id, variable_name, n_vars, &
    &                                         fill_arrays, n_g, glb_index) &
    result(res)

    INTEGER, INTENT(IN)               :: n_vars
    INTEGER, INTENT(IN)               :: file_id
    CHARACTER(LEN=*), INTENT(IN)      :: variable_name
    TYPE(var_data_2d_int), OPTIONAL   :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)               :: n_g
    TYPE(var_data_1d_int), INTENT(IN) :: glb_index(n_vars)

    TYPE(var_data_2d_int)             :: res(n_vars)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status, i
    INTEGER, POINTER :: tmp_array(:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_INT_2D_multivar_fileid'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 .OR. var_size(1) /= n_g) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (var_type /= NF_INT) CALL finish(method_name, "invalid var_type")

    ENDIF

    ALLOCATE( tmp_array(n_g), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_int(file_id, varid, tmp_array(:)), variable_name)
    ENDIF

    DO i = 1, n_vars
      IF (PRESENT(fill_arrays)) THEN
        res(i)%data => fill_arrays(i)%data
      ELSE
        ALLOCATE( res(i)%data(nproma, &
          &                   (SIZE(glb_index(i)%data) - 1)/nproma + 1), &
          &       stat=return_status )
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%data(:,:) = 0
      ENDIF

      CALL scatter_array(in_array=tmp_array, out_array=res(i)%data, &
        &                global_index=glb_index(i)%data)
    END DO

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_INT_2D_multivar_fileid
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_2D_fileid(file_id, variable_name, fill_array, &
    &                                 n_g, glb_index) result(res)

    REAL(wp), POINTER            :: res(:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    INTEGER, INTENT(IN)          :: n_g
    INTEGER, TARGET, INTENT(IN)  :: glb_index(:)

    TYPE(var_data_2d_wp) :: fill_arrays(1)
    TYPE(var_data_1d_int) :: glb_index_(1)
    TYPE(var_data_2d_wp) :: results(1)

    glb_index_(1)%data => glb_index

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%data => fill_array
      results = netcdf_read_REAL_2D_multivar_fileid(file_id=file_id, &
        &                                           variable_name=variable_name,&
        &                                           n_vars=1, &
        &                                           fill_arrays=fill_arrays, &
        &                                           n_g=n_g, &
        &                                           glb_index=glb_index_)
    ELSE
      results = netcdf_read_REAL_2D_multivar_fileid(file_id=file_id, &
        &                                           variable_name=variable_name,&
        &                                           n_vars=1, n_g=n_g, &
        &                                           glb_index=glb_index_)
    END IF

    res => results(1)%data

  END FUNCTION netcdf_read_REAL_2D_fileid

  FUNCTION netcdf_read_REAL_2D_multivar_fileid(file_id, variable_name, n_vars, &
    &                                          fill_arrays, n_g, glb_index) &
    result(res)

    INTEGER, INTENT(IN)               :: n_vars
    INTEGER, INTENT(IN)               :: file_id
    CHARACTER(LEN=*), INTENT(IN)      :: variable_name
    TYPE(var_data_2d_wp), OPTIONAL    :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)               :: n_g
    TYPE(var_data_1d_int), INTENT(IN) :: glb_index(n_vars)

    TYPE(var_data_2d_wp)              :: res(n_vars)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status, i
    REAL(wp), POINTER :: tmp_array(:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_2D_multivar_fileid'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 .OR. var_size(1) /= n_g) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ALLOCATE( tmp_array(n_g), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(file_id, varid, tmp_array(:)), &
        &                       variable_name)
    ENDIF

    DO i = 1, n_vars
      IF (PRESENT(fill_arrays)) THEN
        res(i)%data => fill_arrays(i)%data
      ELSE
        ALLOCATE( res(i)%data(nproma, &
          &                   (SIZE(glb_index(i)%data) - 1)/nproma + 1), &
          &       stat=return_status )
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%data(:,:) = 0.0_wp
      ENDIF

      CALL scatter_array(in_array=tmp_array, out_array=res(i)%data, &
        &                global_index=glb_index(i)%data)
    END DO

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_2D_multivar_fileid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n_g) fortran-style: O3(n_g, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_REAL_2D_time_fileid(file_id, variable_name, fill_array, &
    &                                      n_g, glb_index, start_timestep, &
    &                                      end_timestep) result(res)
    REAL(wp), POINTER             :: res(:,:,:)

    INTEGER, INTENT(IN)           :: file_id
    CHARACTER(LEN=*), INTENT(IN)  :: variable_name
    define_fill_target            :: fill_array(:,:,:)
    INTEGER, INTENT(IN)           :: n_g
    INTEGER, INTENT(IN)           :: glb_index(:)
    INTEGER, INTENT(in), OPTIONAL :: start_timestep, end_timestep

    res => netcdf_read_REAL_2D_extdim_fileid( &
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array, &
      & n_g=n_g, glb_index=glb_index, start_extdim=start_timestep, &
      & end_extdim=end_timestep, extdim_name="time" )

  END FUNCTION netcdf_read_REAL_2D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_REAL_2D_extdim_fileid(file_id, variable_name, &
    &                                        fill_array, n_g, glb_index, &
    &                                        start_extdim, end_extdim, &
    &                                        extdim_name ) result(res)

    REAL(wp), POINTER                      :: res(:,:,:)

    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:)
    INTEGER, INTENT(IN)                    :: n_g
    INTEGER, TARGET, INTENT(IN)            :: glb_index(:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(var_data_3d_wp) :: fill_arrays(1)
    TYPE(var_data_1d_int) :: glb_index_(1)
    TYPE(var_data_3d_wp) :: results(1)

    glb_index_(1)%data => glb_index

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%data => fill_array
      results = netcdf_read_REAL_2D_extdim_multivar_fileid( &
        file_id=file_id, variable_name=variable_name, n_vars=1, &
        fill_arrays=fill_arrays,  n_g=n_g, glb_index=glb_index_, &
        start_extdim=start_extdim, end_extdim=end_extdim, &
        extdim_name=extdim_name)
    ELSE
      results = netcdf_read_REAL_2D_extdim_multivar_fileid( &
        file_id=file_id, variable_name=variable_name, n_vars=1, n_g=n_g, &
        glb_index=glb_index_, start_extdim=start_extdim, &
        end_extdim=end_extdim, extdim_name=extdim_name)
    END IF

    res => results(1)%data

  END FUNCTION netcdf_read_REAL_2D_extdim_fileid

  FUNCTION netcdf_read_REAL_2D_extdim_multivar_fileid(file_id, variable_name, &
    &                                                 n_vars, fill_arrays, n_g,&
    &                                                 glb_index, start_extdim, &
    &                                                 end_extdim, extdim_name) &
    result(res)

    INTEGER, INTENT(IN)                    :: n_vars
    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    TYPE(var_data_3d_wp), OPTIONAL         :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)                    :: n_g
    TYPE(var_data_1d_int), INTENT(IN)      :: glb_index(n_vars)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(var_data_3d_wp)          :: res(n_vars)

    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(2), count_read_index(2)

    INTEGER :: return_status, i
    REAL(wp), POINTER :: tmp_array(:,:)
    REAL(wp), POINTER :: tmp_res(:,:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_2D_extdim_multivar_fileid'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
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
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "ext dim size < 1")

    ALLOCATE( tmp_array(n_g, time_steps), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF( my_process_is_mpi_workroot()) THEN

      start_read_index = (/ 1, start_time /)
      count_read_index = (/ n_g, time_steps /)
      CALL nf(nf_get_vara_double(file_id, varid, start_read_index, &
        &                        count_read_index, tmp_array(:,:)), &
        &                        variable_name)
    ENDIF

    DO i = 1, n_vars
      !-----------------------
      IF (PRESENT(fill_arrays)) THEN
        res(i)%data => fill_arrays(i)%data
        IF (SIZE(res(i)%data,3) < time_steps) &
          CALL finish(method_name, "allocated size < time_steps")
      ELSE
        ALLOCATE(res(i)%data(nproma, &
          &                  (SIZE(glb_index(i)%data) - 1) / nproma + 1, &
          &                  time_steps), stat=return_status)
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%data(:,:,:) = 0.0_wp
      ENDIF
      !-----------------------
      tmp_res => res(i)%data(:,:,LBOUND(res(i)%data, 3):UBOUND(res(i)%data, 3))
      CALL scatter_time_array(tmp_array, tmp_res, glb_index(i)%data)
    END DO

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_2D_extdim_multivar_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_INT_2D_extdim_fileid(file_id, variable_name, &
    &                                       fill_array, n_g, glb_index, &
    &                                       start_extdim, end_extdim, &
    &                                       extdim_name ) result(res)

    INTEGER, POINTER                       :: res(:,:,:)

    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target_int                 :: fill_array(:,:,:)
    INTEGER, INTENT(IN)                    :: n_g
    INTEGER, TARGET, INTENT(IN)            :: glb_index(:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(var_data_3d_int) :: fill_arrays(1)
    TYPE(var_data_1d_int) :: glb_index_(1)
    TYPE(var_data_3d_int) :: results(1)

    glb_index_(1)%data => glb_index

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%data => fill_array
      results = netcdf_read_INT_2D_extdim_multivar_fileid( &
        file_id=file_id, variable_name=variable_name, n_vars=1, &
        fill_arrays=fill_arrays,  n_g=n_g, glb_index=glb_index_, &
        start_extdim=start_extdim, end_extdim=end_extdim, &
        extdim_name=extdim_name)
    ELSE
      results = netcdf_read_INT_2D_extdim_multivar_fileid( &
        file_id=file_id, variable_name=variable_name, n_vars=1, n_g=n_g, &
        glb_index=glb_index_, start_extdim=start_extdim, &
        end_extdim=end_extdim, extdim_name=extdim_name)
    END IF

    res => results(1)%data

  END FUNCTION netcdf_read_INT_2D_extdim_fileid

  FUNCTION netcdf_read_INT_2D_extdim_multivar_fileid(file_id, variable_name, &
    &                                                n_vars, fill_arrays, n_g, &
    &                                                glb_index, start_extdim, &
    &                                                end_extdim, extdim_name ) &
    result(res)

    INTEGER, INTENT(IN)                    :: n_vars
    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    TYPE(var_data_3d_int), OPTIONAL        :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)                    :: n_g
    TYPE(var_data_1d_int), INTENT(IN)      :: glb_index(n_vars)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(var_data_3d_int)          :: res(n_vars)

    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_allocated_step, end_allocated_step
    INTEGER :: start_read_index(2), count_read_index(2)

    INTEGER :: return_status, i
    INTEGER, POINTER :: tmp_array(:,:)
    INTEGER, POINTER :: tmp_res(:,:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_INT_2D_extdim_multivar_fileid'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
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
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "ext dim size < 1")

    ALLOCATE( tmp_array(n_g, time_steps), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF( my_process_is_mpi_workroot()) THEN

      start_read_index = (/ 1, start_time /)
      count_read_index = (/ n_g, time_steps /)
      CALL nf(nf_get_vara_int(file_id, varid, start_read_index, &
        &                     count_read_index, tmp_array(:,:)), variable_name)
    ENDIF

    DO i = 1, n_vars
      !-----------------------
      IF (PRESENT(fill_arrays)) THEN
        res(i)%data => fill_arrays(i)%data
        IF (SIZE(res(i)%data,3) < time_steps) &
          CALL finish(method_name, "allocated size < time_steps")
      ELSE
        ALLOCATE(res(i)%data(nproma, &
          &                  (SIZE(glb_index(i)%data) - 1) / nproma + 1, &
          &                  time_steps), stat=return_status)
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%data(:,:,:) = 0
      END IF
      !-----------------------
      tmp_res => res(i)%data(:,:,LBOUND(res(i)%data, 3):UBOUND(res(i)%data, 3))
      CALL scatter_time_array(tmp_array, tmp_res, glb_index(i)%data)
    END DO

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_INT_2D_extdim_multivar_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O2(levels, n_g) fortran-style: O2(n_g, levels)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks)
  FUNCTION netcdf_read_REAL_3D_fileid(file_id, variable_name, fill_array, n_g, &
    &                                 glb_index, levelsdim_name) result(res)

    REAL(wp), POINTER  :: res(:,:,:)

    INTEGER, INTENT(IN)           :: file_id
    CHARACTER(LEN=*), INTENT(IN)  :: variable_name
    define_fill_target            :: fill_array(:,:,:)
    INTEGER, INTENT(IN)           :: n_g
    INTEGER, INTENT(IN)           :: glb_index(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_vertical_levels

    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_3D_fileid'

    NULLIFY(res)

    IF( my_process_is_mpi_workroot() ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 2 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL broadcast_array(var_size(1:2))
    file_vertical_levels = var_size(2)

    !-----------------------
    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res (nproma, file_vertical_levels, &
        &            (SIZE(glb_index) - 1)/nproma + 1), &
        &       stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( res )')
      ENDIF
      res(:,:,:) = 0.0_wp
    ENDIF
    ALLOCATE( tmp_array(n_g, file_vertical_levels), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF (file_vertical_levels /= SIZE(res,2)) &
      CALL finish(method_name, 'file_vertical_levels /= SIZE(fill_array,2)')
    !-----------------------

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(file_id, varid, tmp_array(:,:)), &
        &     variable_name)
    ENDIF

    CALL scatter_array(tmp_array, res, glb_index)

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_3D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION netcdf_read_REAL_3D_time_fileid(file_id, variable_name, fill_array, &
    &                                      n_g, glb_index, start_timestep, &
    &                                      end_timestep, levelsdim_name) &
    & result(res)

    REAL(wp), POINTER  :: res(:,:,:,:)

    INTEGER, INTENT(IN)           :: file_id
    CHARACTER(LEN=*), INTENT(IN)  :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    INTEGER, INTENT(IN)           :: n_g
    INTEGER, INTENT(IN)           :: glb_index(:)
    INTEGER, INTENT(in), OPTIONAL :: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    res => netcdf_read_REAL_3D_extdim_fileid(&
      & file_id=file_id,                        &
      & variable_name=variable_name,            &
      & fill_array=fill_array,                  &
      & n_g=n_g, glb_index=glb_index,           &
      & start_extdim=start_timestep,            &
      & end_extdim=end_timestep,                &
      & levelsdim_name=levelsdim_name,          &
      & extdim_name="time")

  END FUNCTION netcdf_read_REAL_3D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n_g) fortran-style: O3(n_g, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION netcdf_read_REAL_3D_extdim_fileid(file_id, variable_name, &
    &                                        fill_array, n_g, glb_index, &
    &                                        start_extdim, end_extdim, &
    &                                        levelsdim_name, extdim_name ) &
    &  result(res)

    REAL(wp), POINTER  :: res(:,:,:,:)

    INTEGER, INTENT(IN)           :: file_id
    CHARACTER(LEN=*), INTENT(IN)  :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    INTEGER, INTENT(IN)           :: n_g
    INTEGER, INTENT(IN)           :: glb_index(:)
    INTEGER, INTENT(in), OPTIONAL :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsdim_name

    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_vertical_levels, file_time_steps, time_steps, start_time, &
      &        end_time
    INTEGER :: start_read_index(3), count_read_index(3)

    INTEGER :: return_status, i, tt
    REAL(wp), POINTER :: tmp_array(:)
    REAL(wp), POINTER  :: res_level(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_3D_extdim_fileid'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 3 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
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

    time_steps = end_time - start_time + 1
    IF (time_steps < 1) &
      & CALL finish(method_name, "extdim size < 1")

    ALLOCATE( tmp_array(n_g), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    IF (PRESENT(fill_array)) THEN
      res => fill_array(:,:,:,1:time_steps)
    ELSE
      ALLOCATE( res (nproma, file_vertical_levels, &
        &            (SIZE(glb_index) - 1)/nproma + 1, time_steps), &
        &       stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( res )')
      ENDIF
      res(:,:,:,:) = 0.0_wp
    ENDIF

    IF (SIZE(res,4) < time_steps) &
      CALL finish(method_name, "allocated size < time_steps")

    IF (file_vertical_levels /= SIZE(res,2)) &
      CALL finish(method_name, 'file_vertical_levels /= SIZE(fill_array,2)')
    !-----------------------

    DO tt=1, time_steps
      DO i=1, file_vertical_levels
        IF( my_process_is_mpi_workroot()) THEN
          start_read_index = (/ 1, i, tt + start_time - 1 /)
          count_read_index      = (/ n_g, 1, 1 /)
          CALL nf(nf_get_vara_double(file_id, varid, start_read_index, &
            &                        count_read_index, tmp_array(:)),  &
            &     variable_name)
        ENDIF
        
        res_level => res(:,i,:,LBOUND(res, 4)+tt-1)
        CALL scatter_array(tmp_array, res_level, glb_index)
      END DO
    END DO

    DEALLOCATE(tmp_array)

  END FUNCTION netcdf_read_REAL_3D_extdim_fileid
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
  SUBROUTINE netcdf_inq_var(file_id, name, varid, var_type, var_dims, var_size, &
    &                       var_dim_name)
    INTEGER, INTENT(IN) :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: name

    INTEGER, INTENT(OUT) :: varid, var_type, var_dims
    INTEGER, INTENT(OUT) :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max), INTENT(OUT) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER  :: number_of_attributes
    INTEGER :: var_dims_reference(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: check_var_name
    INTEGER :: i

    IF ( .NOT. my_process_is_mpi_workroot() ) RETURN

    CALL nf(nf_inq_varid(file_id, name, varid), name)
    CALL nf(nf_inq_var(file_id, varid, check_var_name, var_type, var_dims, &
      &                var_dims_reference, number_of_attributes), &
      &     check_var_name)
    DO i=1, var_dims
      CALL nf(nf_inq_dimlen (file_id, var_dims_reference(i), var_size(i)), &
        &     check_var_name)
      CALL nf(nf_inq_dimname(file_id, var_dims_reference(i), var_dim_name(i)), &
        &     check_var_name)
    ENDDO
!     write(0,*) " Read var_dims, var_size:",  var_dims, var_size
!     write(0,*) " check_var_name:",  check_var_name
!     write(0,*) " name:", name

  END SUBROUTINE netcdf_inq_var
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

END MODULE mo_read_netcdf_broadcast_2
