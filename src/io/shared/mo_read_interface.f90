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
!! Moved from mo_util_netcdf, added read_onCells, by L. Linardakis, (2013-03-15)
!! Clean up and refactoring of the module, by M. Hanke, DKRZ, (2014-06)
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
#define define_return_pointer REAL(wp), POINTER, OPTIONAL

MODULE mo_read_interface

  USE mo_kind
  USE mo_exception,          ONLY: finish, message_text, message, em_warn
  USE mo_io_config,         ONLY:  read_netcdf_broadcast_method, &
    & read_netcdf_distribute_method,  default_read_method
  USE mo_netcdf_read,      ONLY: netcdf_open_input, netcdf_close, &
    &                            netcdf_read_2D, netcdf_read_2D_extdim, &
    &                            netcdf_read_3D_extdim, netcdf_read_0D_real, &
    &                            netcdf_read_1D, netcdf_read_3D, &
    &                            netcdf_read_1D_extdim_time, &
    &                            netcdf_read_1D_extdim_extdim_time
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data, distrib_nf_open, &
    &                                   distrib_read, distrib_nf_close, &
    &                                   var_data_2d_wp, distrib_inq_var_dims
  USE mo_model_domain, ONLY: t_patch
  USE mo_parallel_config, ONLY: nproma
  USE mo_io_units, ONLY: filename_max


  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: read_netcdf_broadcast_method, read_netcdf_distribute_method
  PUBLIC :: t_stream_id

  PUBLIC :: openInputFile, closeFile

  PUBLIC :: read_0D_real
  PUBLIC :: read_1D
  PUBLIC :: read_1D_extdim_time
  PUBLIC :: read_1D_extdim_extdim_time
  PUBLIC :: read_2D
  PUBLIC :: read_2D_time
  PUBLIC :: read_3D_time
  PUBLIC :: read_2D_extdim
  PUBLIC :: read_3D_extdim

  PUBLIC :: onCells, onVertices, onEdges

  !--------------------------------------------------------

  INTEGER, PARAMETER :: onCells = 1
  INTEGER, PARAMETER :: onVertices = 2
  INTEGER, PARAMETER :: onEdges = 3

  TYPE t_read_info
    ! data required for read_netcdf_distribute_method
    TYPE(t_distrib_read_data), POINTER :: dist_read_info

    ! data required for read_netcdf_broadcast_method
    INTEGER :: n_g
    INTEGER, POINTER :: glb_index(:)
  END TYPE

  TYPE t_stream_id
    INTEGER  :: file_id             ! netcdf file id, or similar
    INTEGER  :: return_status       ! the latest operation return status
    INTEGER  :: input_method        ! read_netcdf_broadcast_method,
                                    ! read_netcdf_distribute_method, etc

    TYPE(t_read_info) :: read_info(3)
  END TYPE t_stream_id
  !--------------------------------------------------------

  INTERFACE openInputFile
    MODULE PROCEDURE openInputFile_dist
    MODULE PROCEDURE openInputFile_bcast
  END INTERFACE openInputFile

  INTERFACE closeFile
    MODULE PROCEDURE closeFile_dist
    MODULE PROCEDURE closeFile_bcast
  END INTERFACE closeFile

  INTERFACE read_0D_real
    MODULE PROCEDURE read_bcast_REAL_0D_streamid
  END INTERFACE read_0D_real

  INTERFACE read_1D
    MODULE PROCEDURE read_bcast_REAL_1D_streamid
  END INTERFACE read_1D

  INTERFACE read_1D_extdim_time
    MODULE PROCEDURE read_bcast_REAL_1D_extdim_time_streamid
  END INTERFACE read_1D_extdim_time

  INTERFACE read_1D_extdim_extdim_time
    MODULE PROCEDURE read_bcast_REAL_1D_extdim_extdim_time_streamid
  END INTERFACE read_1D_extdim_extdim_time

  INTERFACE read_2D
    MODULE PROCEDURE read_dist_REAL_2D_streamid
    MODULE PROCEDURE read_dist_REAL_2D_multivar_streamid
  END INTERFACE read_2D

  INTERFACE read_2D_time
    MODULE PROCEDURE read_dist_REAL_2D_time_streamid
  END INTERFACE read_2D_time

  INTERFACE read_2D_extdim
    MODULE PROCEDURE read_dist_REAL_2D_extdim_streamid
  END INTERFACE read_2D_extdim

  INTERFACE read_3D
    MODULE PROCEDURE read_dist_REAL_3D_streamid
  END INTERFACE read_3D

  INTERFACE read_3D_time
    MODULE PROCEDURE read_dist_REAL_3D_time_streamid
  END INTERFACE read_3D_time

  INTERFACE read_3D_extdim
    MODULE PROCEDURE read_dist_REAL_3D_extdim_streamid
  END INTERFACE read_3D_extdim

CONTAINS

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_bcast_REAL_0D_streamid(file_id, variable_name)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp)    :: read_bcast_REAL_0D_streamid


    read_bcast_REAL_0D_streamid = netcdf_read_0D_real(file_id, variable_name)

  END FUNCTION read_bcast_REAL_0D_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_1D_streamid(file_id, variable_name, fill_array, &
    &                                    return_pointer)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)
    define_return_pointer        :: return_pointer(:)

    REAL(wp), POINTER            :: tmp_pointer(:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_bcast_REAL_1D_streamid'

    tmp_pointer => netcdf_read_1D(file_id, variable_name, fill_array)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_1D_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_1D_extdim_time_streamid(file_id, variable_name, &
    &                                                fill_array, return_pointer, &
    &                                                dim_names, start_timestep, &
    &                                                end_timestep)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_bcast_REAL_1D_extdim_time_streamid'

    tmp_pointer => netcdf_read_1D_extdim_time(file_id, variable_name, &
      &                                       fill_array, dim_names, &
      &                                       start_timestep, end_timestep)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_1D_extdim_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_1D_extdim_extdim_time_streamid( &
    file_id, variable_name, fill_array, return_pointer, dim_names, &
    start_timestep, end_timestep)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_bcast_REAL_1D_extdim_extdim_time_streamid'

    tmp_pointer => netcdf_read_1D_extdim_extdim_time(file_id, variable_name, &
      &                                              fill_array, dim_names, &
      &                                              start_timestep, end_timestep)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_1D_extdim_extdim_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_REAL_2D_multivar_streamid(stream_ids, location, &
    &                                            variable_name, n_var, &
    &                                            fill_array, return_pointer)

    TYPE(t_stream_id), INTENT(INOUT)  :: stream_ids(:)
    INTEGER, INTENT(IN)               :: location
    CHARACTER(LEN=*), INTENT(IN)      :: variable_name
    INTEGER, INTENT(IN)               :: n_var
    TYPE(var_data_2d_wp), OPTIONAL    :: fill_array(:)
    TYPE(var_data_2d_wp), OPTIONAL    :: return_pointer(:)

    REAL(wp), POINTER                 :: tmp_pointer(:,:)
    TYPE(var_data_2d_wp), ALLOCATABLE :: var_data_2d(:)
    CHARACTER(LEN=*), PARAMETER       :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_multivar_streamid'
    TYPE(t_distrib_read_data), ALLOCATABLE :: dist_read_info(:)

    INTEGER :: i

    IF (SIZE(stream_ids) /= n_var) &
      CALL finish(method_name, "number of stream ids does not match number of variables")

    DO i = 2, SIZE(stream_ids)
      IF (stream_ids(i)%file_id /= stream_ids(1)%file_id) &
        CALL finish(method_name, "file ids do not match")
      IF (stream_ids(i)%input_method /= stream_ids(1)%input_method) &
        CALL finish(method_name, "input methods do not match")
    END DO

    CALL check_dimensions(stream_ids(1)%file_id, variable_name, 1, &
      &                   stream_ids(1)%read_info(location)%n_g, location)

    SELECT CASE(stream_ids(1)%input_method)
    CASE (read_netcdf_broadcast_method)
      DO i = 1, n_var

        tmp_pointer => netcdf_read_2D(stream_ids(i)%file_id, variable_name, &
          &                           fill_array(i)%data, &
          &                           stream_ids(i)%read_info(location)%n_g, &
          &                           stream_ids(i)%read_info(location)%glb_index)

        IF (PRESENT(return_pointer)) return_pointer(i)%data => tmp_pointer
      END DO
    CASE (read_netcdf_distribute_method)
    
      ALLOCATE(var_data_2d(n_var), dist_read_info(n_var))

      IF (PRESENT(fill_array)) THEN
        DO i = 1, n_var
          var_data_2d(i)%data => fill_array(i)%data
        END DO
      ELSE
        DO i = 1, n_var
          ALLOCATE(var_data_2d(i)%data(nproma, &
            (SIZE(stream_ids(i)%read_info(location)%glb_index) - 1)/nproma + 1))
          var_data_2d(i)%data(:,:) = 0.0_wp
        END DO
      ENDIF
      IF (PRESENT(return_pointer)) THEN
        DO i = 1, n_var
          return_pointer(i)%data => var_data_2d(i)%data
        END DO
      END IF
      DO i = 1, n_var
        dist_read_info(i) = stream_ids(i)%read_info(location)%dist_read_info
      END DO
      CALL distrib_read(stream_ids(1)%file_id, variable_name, var_data_2d, &
        &               dist_read_info(:))
      DEALLOCATE(var_data_2d)
      DEALLOCATE(dist_read_info)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_multivar_streamid

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_REAL_2D_streamid(stream_id, location, variable_name, &
    &                                   fill_array, return_pointer)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    define_return_pointer        :: return_pointer(:,:)

    REAL(wp), POINTER            :: tmp_pointer(:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_streamid'

    CALL check_dimensions(stream_id%file_id, variable_name, 1, &
      &                   stream_id%read_info(location)%n_g, location)

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_2D(stream_id%file_id, variable_name, &
        &                           fill_array, &
        &                           stream_id%read_info(location)%n_g, &
        &                           stream_id%read_info(location)%glb_index)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        ALLOCATE(tmp_pointer(nproma, &
          (SIZE(stream_id%read_info(location)%glb_index) - 1)/nproma + 1))
        tmp_pointer(:,:) = 0.0_wp
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    
      CALL distrib_read(stream_id%file_id, variable_name, tmp_pointer, &
        &               stream_id%read_info(location)%dist_read_info)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_REAL_2D_time_streamid(stream_id, location, &
    &                                        variable_name, fill_array, &
    &                                        return_pointer, start_timestep, &
    &                                        end_timestep)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep

    CALL check_dimensions(stream_id%file_id, variable_name, 2, &
      &                   stream_id%read_info(location)%n_g, location, &
      &                   (/"time"/))

    CALL read_dist_REAL_2D_extdim_streamid(&
      & stream_id=stream_id, location=location, variable_name=variable_name, &
      & fill_array=fill_array, return_pointer=return_pointer, &
      & start_extdim=start_timestep, end_extdim=end_timestep, &
      & extdim_name="time" )

  END SUBROUTINE read_dist_REAL_2D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_REAL_2D_extdim_streamid(stream_id, location, &
    &                                          variable_name, fill_array, &
    &                                          return_pointer, start_extdim, &
    &                                          end_extdim, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_extdim_streamid'

    IF (PRESENT(extdim_name)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name, 2, &
        &                   stream_id%read_info(location)%n_g, location, &
        &                   (/extdim_name/))
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name, 2, &
        &                   stream_id%read_info(location)%n_g, location)
    END IF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_2D_extdim(stream_id%file_id, variable_name, fill_array, &
         & stream_id%read_info(location)%n_g, &
         & stream_id%read_info(location)%glb_index, start_extdim, end_extdim, &
         & extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        ALLOCATE(tmp_pointer(nproma, &
          (SIZE(stream_id%read_info(location)%glb_index) - 1)/nproma + 1, &
          end_extdim - start_extdim + 1))
        tmp_pointer(:,:,:) = 0.0_wp
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    
      CALL distrib_read(stream_id%file_id, variable_name, tmp_pointer, &
        &               end_extdim - start_extdim + 1, &
        &               stream_id%read_info(location)%dist_read_info)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_extdim_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O2(levels, n) fortran-style: O2(n, levels)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks)
  SUBROUTINE read_dist_REAL_3D_streamid(stream_id, location, &
    &                                   variable_name, fill_array, &
    &                                   return_pointer, levelsDimName)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName

    INTEGER :: var_ndims, var_dimlen(2)
    REAL(wp), POINTER  :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_3D_streamid'

    IF (PRESENT(levelsDimName)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name, 2, &
        &                   stream_id%read_info(location)%n_g, location, &
        &                   (/levelsDimName/))
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name, 2, &
        &                   stream_id%read_info(location)%n_g, location)
    END IF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_3D(stream_id%file_id, variable_name, &
         & fill_array, stream_id%read_info(location)%n_g, &
         & stream_id%read_info(location)%glb_index, levelsDimName)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        CALL distrib_inq_var_dims(stream_id%file_id, variable_name, var_ndims, &
          &                       var_dimlen)
        ALLOCATE(tmp_pointer(nproma, var_dimlen(2), &
          (SIZE(stream_id%read_info(location)%glb_index) - 1)/nproma + 1))
        tmp_pointer(:,:,:) = 0.0_wp
      ENDIF

      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    
      CALL distrib_read(stream_id%file_id, variable_name, tmp_pointer, &
        &               var_dimlen(2), &
        &               stream_id%read_info(location)%dist_read_info)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_3D_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  SUBROUTINE read_dist_REAL_3D_time_streamid(stream_id, location, &
    &                                        variable_name, fill_array, &
    &                                        return_pointer, start_timestep, &
    &                                        end_timestep, levelsDimName)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName

    IF (PRESENT(levelsDimName)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name, 3, &
        &                   stream_id%read_info(location)%n_g, location, &
        &                   (/levelsDimName, "time"/))
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name, 3, &
        &                   stream_id%read_info(location)%n_g, location)
    END IF

    CALL read_dist_REAL_3D_extdim_streamid( &
      & stream_id=stream_id,                 &
      & location=location,                   &
      & variable_name=variable_name,         &
      & fill_array=fill_array,               &
      & return_pointer=return_pointer,       &
      & start_extdim=start_timestep,         &
      & end_extdim=end_timestep,             &
      & levelsDimName=levelsDimName,         &
      & extdim_name="time")

  END SUBROUTINE read_dist_REAL_3D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  SUBROUTINE read_dist_REAL_3D_extdim_streamid(stream_id, location, &
    &                                          variable_name, fill_array, &
    &                                          return_pointer, start_extdim, &
    &                                          end_extdim, levelsDimName, &
    &                                          extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsDimName

    INTEGER :: var_ndims, var_dimlen(3)
    REAL(wp), POINTER  :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_3D_extdim_streamid'

    IF (PRESENT(levelsDimName) .AND. PRESENT(extdim_name)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name, 3, &
        &                   stream_id%read_info(location)%n_g, location, &
        &                   (/levelsDimName, extdim_name/))
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name, 3, &
        &                   stream_id%read_info(location)%n_g, location)
    END IF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_3D_extdim(stream_id%file_id, variable_name, &
         & fill_array, stream_id%read_info(location)%n_g, &
         & stream_id%read_info(location)%glb_index, &
         & start_extdim, end_extdim, levelsDimName, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        CALL distrib_inq_var_dims(stream_id%file_id, variable_name, var_ndims, &
          &                       var_dimlen)
        ALLOCATE(tmp_pointer(nproma, var_dimlen(2), &
          (SIZE(stream_id%read_info(location)%glb_index) - 1)/nproma + 1, &
          end_extdim - start_extdim + 1))
        tmp_pointer(:,:,:,:) = 0.0_wp
      ENDIF

      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    
      CALL distrib_read(stream_id%file_id, variable_name, tmp_pointer, &
        &               var_dimlen(2:3), &
        &               stream_id%read_info(location)%dist_read_info)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT


  END SUBROUTINE read_dist_REAL_3D_extdim_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  TYPE(t_stream_id) FUNCTION openInputFile_dist(filename, patch, input_method)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(t_patch), TARGET, INTENT(IN) :: patch
    INTEGER, OPTIONAL, INTENT(IN) :: input_method

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile_dist'

    IF (PRESENT(input_method)) THEN
      openInputFile_dist%input_method = input_method
    ELSE
      openInputFile_dist%input_method = default_read_method
    END IF

    openInputFile_dist%read_info(onCells)%n_g = patch%n_patch_cells_g
    openInputFile_dist%read_info(onEdges)%n_g = patch%n_patch_edges_g
    openInputFile_dist%read_info(onVertices)%n_g = patch%n_patch_verts_g

    SELECT CASE(openInputFile_dist%input_method)
    CASE (read_netcdf_broadcast_method)

      openInputFile_dist%file_id = netcdf_open_input(filename)

      openInputFile_dist%read_info(onCells)%glb_index => &
        patch%cells%decomp_info%glb_index
      NULLIFY(openInputFile_dist%read_info(onCells)%dist_read_info)

      openInputFile_dist%read_info(onEdges)%glb_index => &
        patch%edges%decomp_info%glb_index
      NULLIFY(openInputFile_dist%read_info(onEdges)%dist_read_info)

      openInputFile_dist%read_info(onVertices)%glb_index => &
        patch%verts%decomp_info%glb_index
      NULLIFY(openInputFile_dist%read_info(onVertices)%dist_read_info)

    CASE (read_netcdf_distribute_method)

      openInputFile_dist%file_id = distrib_nf_open(TRIM(filename))

      openInputFile_dist%read_info(onCells)%dist_read_info => &
        patch%cells%dist_io_data
      openInputFile_dist%read_info(onCells)%glb_index => &
        patch%cells%decomp_info%glb_index

      openInputFile_dist%read_info(onVertices)%dist_read_info => &
        patch%verts%dist_io_data
      openInputFile_dist%read_info(onVertices)%glb_index => &
        patch%verts%decomp_info%glb_index

      openInputFile_dist%read_info(onEdges)%dist_read_info => &
        patch%edges%dist_io_data
      openInputFile_dist%read_info(onEdges)%glb_index => &
        patch%edges%decomp_info%glb_index

    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END FUNCTION openInputFile_dist
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION openInputFile_bcast(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile_bcast'

    openInputFile_bcast = netcdf_open_input(filename)

  END FUNCTION openInputFile_bcast
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE closeFile_dist(stream_id)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:closeFile_dist'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      stream_id%return_status = netcdf_close(stream_id%file_id)
    CASE (read_netcdf_distribute_method)
      CALL distrib_nf_close(stream_id%file_id)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE closeFile_dist
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE closeFile_bcast(file_id, return_status)
    INTEGER, INTENT(IN) :: file_id
    INTEGER, OPTIONAL, INTENT(OUT) :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:closeFile_bcast'

    return_status = netcdf_close(file_id)

  END SUBROUTINE closeFile_bcast

  SUBROUTINE check_dimensions(file_id, variable_name, ref_var_ndims, n_g, &
                              location, extdim_name)

    INTEGER, INTENT(IN) :: file_id
    CHARACTER(LEN=*), INTENT(IN)  :: variable_name
    INTEGER, INTENT(IN) :: ref_var_ndims, n_g
    INTEGER, INTENT(IN) :: location
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name(2:ref_var_ndims)

    INTEGER :: varid, var_ndims, var_dimlen(NF_MAX_VAR_DIMS), &
      &        var_dimids(NF_MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(NF_MAX_VAR_DIMS)
    INTEGER :: i

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:check_dimensions'

    IF (file_id == -1) RETURN

    CALL nf(nf_inq_varid(file_id, variable_name, varid), method_name)
    CALL nf(nf_inq_varndims(file_id, varid, var_ndims), method_name)
    CALL nf(nf_inq_vardimid(file_id, varid, var_dimids), method_name)
    DO i = 1, var_ndims

      CALL nf(nf_inq_dimlen (file_id, var_dimids(i), var_dimlen(i)), method_name)
      CALL nf(nf_inq_dimname(file_id, var_dimids(i), var_dim_name(i)), method_name)
    END DO

    IF (var_ndims /= ref_var_ndims ) THEN
      WRITE(0,*) variable_name, ": var_ndims = ", var_ndims
      CALL finish(method_name, "Dimensions mismatch")
    ENDIF

    ! check if the input has the right shape/size
    IF ( var_dimlen(1) /= n_g) THEN
      WRITE(0,*) variable_name, ": var_ndims = ", var_ndims, " var_dimlen=", &
        &        var_dimlen, " n_g=", n_g
      CALL finish(method_name, "Dimensions mismatch")
    ENDIF

    ! check if the dim have reasonable names

    SELECT CASE(location)
      CASE (onCells)
        IF ((TRIM(var_dim_name(1)) == 'cell') .EQV. &
          & (TRIM(var_dim_name(1)) == 'ncells')) THEN
          write(0,*) var_dim_name(1)
          WRITE(message_text,*) variable_name, " ", TRIM(var_dim_name(1)), &
            &                   " /= std_cells_dim_name"
          CALL finish(method_name, message_text)
        ENDIF
      CASE (onVertices)
        IF ((TRIM(var_dim_name(1)) == 'vertex') .EQV. &
          & (TRIM(var_dim_name(1)) == 'nverts')) THEN
          write(0,*) var_dim_name(1)
          WRITE(message_text,*) variable_name, " ", TRIM(var_dim_name(1)), &
            &                   " /= std_verts_dim_name"
          CALL finish(method_name, message_text)
        ENDIF
      CASE (onEdges)
        IF ((TRIM(var_dim_name(1)) == 'edge') .EQV. &
          & (TRIM(var_dim_name(1)) == 'nedges')) THEN
          write(0,*) var_dim_name(1)
          WRITE(message_text,*) variable_name, " ", TRIM(var_dim_name(1)), &
            &                   " /= std_edge_dim_name"
          CALL finish(method_name, message_text)
        ENDIF
    END SELECT

    IF (PRESENT(extdim_name)) THEN
      DO i = 2, ref_var_ndims
        IF (TRIM(extdim_name(i)) /= TRIM(var_dim_name(i))) THEN
          WRITE(message_text,*) variable_name, ":", TRIM(extdim_name(i)), &
            &                   "/=",  var_dim_name(i)
          CALL finish(method_name, message_text)
        ENDIF
      END DO
    END IF

  END SUBROUTINE check_dimensions

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

END MODULE mo_read_interface
