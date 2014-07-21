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
  USE mo_exception,          ONLY: finish
  USE mo_io_config,         ONLY:  read_netcdf_broadcast_method, &
    & read_netcdf_distribute_method,  default_read_method

  USE mo_netcdf_read,      ONLY: netcdf_open_input, netcdf_close, &
    &                            netcdf_read_2D, netcdf_read_2D_extdim, &
    &                            netcdf_read_3D_extdim, netcdf_read_0D_real, &
    &                            netcdf_read_1D, netcdf_read_2D_time, &
    &                            netcdf_read_3D_time
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data


  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_netcdf_broadcast_method, read_netcdf_distribute_method
  PUBLIC :: t_stream_id

  PUBLIC :: openInputFile, closeFile

  PUBLIC :: read_0D_real
  PUBLIC :: read_1D
  PUBLIC :: read_2D
  PUBLIC :: read_2D_time
  PUBLIC :: read_3D_time
  PUBLIC :: read_2D_extdim
  PUBLIC :: read_3D_extdim

  !--------------------------------------------------------
  TYPE t_stream_id
    INTEGER  :: file_id             ! netcdf file id, or similar
    INTEGER  :: return_status       ! the latest operation return status
    INTEGER  :: input_method        ! read_netcdf_broadcast_method,
                                    ! read_netcdf_distribute_method, etc

    ! data required for read_netcdf_distribute_method
    TYPE(t_distrib_read_data), POINTER :: dist_read_info

    ! data required for read_netcdf_broadcast_method
    INTEGER :: n_g
    INTEGER, POINTER :: glb_index(:)
  END TYPE t_stream_id
  !--------------------------------------------------------

  INTERFACE read_0D_real
    MODULE PROCEDURE read_bcast_REAL_0D_streamid
  END INTERFACE read_0D_real

  INTERFACE read_1D
    MODULE PROCEDURE read_bcast_REAL_1D_streamid
  END INTERFACE read_1D

  INTERFACE read_2D
    MODULE PROCEDURE read_dist_REAL_2D_streamid
  END INTERFACE read_2D

  INTERFACE read_2D_time
    MODULE PROCEDURE read_bcast_REAL_2D_time_streamid
    MODULE PROCEDURE read_dist_REAL_2D_time_streamid
  END INTERFACE read_2D_time

  INTERFACE read_2D_extdim
    MODULE PROCEDURE read_dist_REAL_2D_1extdim_streamid
  END INTERFACE read_2D_extdim

  INTERFACE read_3D_time
    MODULE PROCEDURE read_bcast_REAL_3D_time_streamid
    MODULE PROCEDURE read_dist_REAL_3D_time_streamid
  END INTERFACE read_3D_time

  INTERFACE read_3D_extdim
    MODULE PROCEDURE read_dist_REAL_3D_1extdim_streamid
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
  SUBROUTINE read_bcast_REAL_2D_time_streamid(file_id, variable_name, &
    &                                         fill_array, return_pointer, &
    &                                         dim_names, start_timestep, &
    &                                         end_timestep)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_bcast_REAL_2D_time_streamid'

    tmp_pointer => netcdf_read_2D_time(file_id, variable_name, fill_array, &
      &                                dim_names, start_timestep, &
      &                                end_timestep)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_2D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_3D_time_streamid(file_id, variable_name, &
    &                                         fill_array, return_pointer, &
    &                                         dim_names, start_timestep, &
    &                                         end_timestep)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_bcast_REAL_3D_time_streamid'

    tmp_pointer => netcdf_read_3D_time(file_id, variable_name, fill_array, &
      dim_names, start_timestep, end_timestep)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_3D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_REAL_2D_streamid(stream_id, variable_name, fill_array, &
    &                              return_pointer)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    define_return_pointer        :: return_pointer(:,:)

    REAL(wp), POINTER            :: tmp_pointer(:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_2D(stream_id%file_id, variable_name, &
        &                           fill_array, stream_id%n_g, &
        &                           stream_id%glb_index)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
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
  SUBROUTINE read_dist_REAL_2D_time_streamid(stream_id, variable_name, &
    &                                        fill_array, return_pointer, &
    &                                        start_timestep, end_timestep)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep

    CALL read_dist_REAL_2D_1extdim_streamid(&
      & stream_id=stream_id, variable_name=variable_name, fill_array=fill_array, &
      & return_pointer=return_pointer, start_extdim=start_timestep, &
      & end_extdim=end_timestep, extdim_name="time" )

  END SUBROUTINE read_dist_REAL_2D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_REAL_2D_1extdim_streamid(stream_id, variable_name, &
    &                                           fill_array, return_pointer, &
    &                                           start_extdim, end_extdim, &
    &                                           extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_1extdim_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_2D_extdim(stream_id%file_id, variable_name, fill_array, &
         & stream_id%n_g, stream_id%glb_index, start_extdim, end_extdim, &
         & extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_1extdim_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  SUBROUTINE read_dist_REAL_3D_time_streamid(stream_id, variable_name, &
    &                                        fill_array, return_pointer, &
    &                                        start_timestep, end_timestep, &
    &                                        levelsDimName)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName

    CALL read_dist_REAL_3D_1extdim_streamid( &
      & stream_id=stream_id,                 &
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
  SUBROUTINE read_dist_REAL_3D_1extdim_streamid(stream_id, variable_name, &
    &                                           fill_array, return_pointer, &
    &                                           start_extdim, end_extdim, &
    &                                           levelsDimName, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsDimName

    REAL(wp), POINTER  :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_dist_REAL_3D_1extdim_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_3D_extdim(stream_id%file_id, variable_name, &
         & fill_array, stream_id%n_g, stream_id%glb_index, &
         & start_extdim, end_extdim, levelsDimName, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT


  END SUBROUTINE read_dist_REAL_3D_1extdim_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  TYPE(t_stream_id) FUNCTION openInputFile(filename, input_method, n_g, &
    &                                      glb_index, dist_read_info)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, OPTIONAL, INTENT(IN) :: input_method

    TYPE(t_distrib_read_data), OPTIONAL, POINTER :: dist_read_info
    INTEGER, OPTIONAL, INTENT(IN) :: n_g
    INTEGER, OPTIONAL, POINTER :: glb_index(:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile'

    IF (present(input_method)) THEN
      openInputFile%input_method = input_method
    ELSE
      openInputFile%input_method = default_read_method
    ENDIF

    SELECT CASE(openInputFile%input_method)
    CASE (read_netcdf_broadcast_method)
      IF ((.NOT. PRESENT(n_g)) .OR. (.NOT. PRESENT(glb_index))) &
        CALL finish(method_name, &
          &         "input method read_netcdf_broadcast_method: " // &
          &         "n_g or glb_index is missing")
      openInputFile%file_id = netcdf_open_input(filename)
      openInputFile%n_g = n_g
      openInputFile%glb_index = glb_index
    CASE (read_netcdf_distribute_method)
      IF (.NOT. PRESENT(dist_read_info)) &
        CALL finish(method_name, &
          &         "input method read_netcdf_distribute_method: " // &
          &         "dist_read_info is missing")
      openInputFile%dist_read_info = dist_read_info
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END FUNCTION openInputFile
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE closeFile(stream_id)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:closeFile'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      stream_id%return_status = netcdf_close(stream_id%file_id)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE closeFile
  !-------------------------------------------------------------------------

END MODULE mo_read_interface
