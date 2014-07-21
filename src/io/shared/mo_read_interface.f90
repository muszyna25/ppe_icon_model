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

  USE mo_io_config,         ONLY:  read_netcdf_broadcast_method, &
    & read_netcdf_distribute_method,  default_read_method

  USE mo_netcdf_read,      ONLY: netcdf_open_input, netcdf_close, &
    & netcdf_read_0D_real, netcdf_read_1D, netcdf_read_2D_time, netcdf_read_3D_time, &
    & netcdf_read_2D, netcdf_read_2D_extdim, netcdf_read_3D_extdim
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data


  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_netcdf_broadcast_method, read_netcdf_distribute_method
  PUBLIC :: t_stream_id

  PUBLIC :: openInputFile, closeFile

  PUBLIC :: read_0D_real
  PUBLIC :: read_1D
  PUBLIC :: read_2D_time
  PUBLIC :: read_3D_time
  PUBLIC :: read_onCells_2D
  PUBLIC :: read_onCells_2D_time
  PUBLIC :: read_onCells_3D_time
  PUBLIC :: read_onCells_2D_extdim
  PUBLIC :: read_onCells_3D_extdim

  !--------------------------------------------------------
  TYPE t_stream_id
    INTEGER  :: stream_id           ! future use
    INTEGER  :: file_id             ! netcdf file id, or similar
    INTEGER  :: return_status       ! the latest operation return status
    INTEGER  :: current_state       ! the state of the stream, opened, prefetching, etc. future use
    INTEGER  :: input_method        ! read_netcdf_broadcast_method, read_netcdf_distribute_method, etc
    TYPE(t_patch), POINTER :: patch ! the patch associated with the stream
  END TYPE t_stream_id
  !--------------------------------------------------------

  TYPE t_read_info
    TYPE(t_distrib_read_data) :: dist_read_info
    INTEGER :: n_g
    INTEGER, POINTER :: glb_index(:)
  END TYPE t_read_info

  INTERFACE read_0D_real
    MODULE PROCEDURE read_REAL_0D_streamid
  END INTERFACE read_0D_real

  INTERFACE read_1D
    MODULE PROCEDURE read_REAL_1D_streamid
  END INTERFACE read_1D

  INTERFACE read_2D_time
    MODULE PROCEDURE read_REAL_2D_time_streamid
  END INTERFACE read_2D_time

  INTERFACE read_3D_time
    MODULE PROCEDURE read_REAL_3D_time_streamid
  END INTERFACE read_3D_time

  INTERFACE read_onCells_2D
    MODULE PROCEDURE read_REAL_onCells_2D_streamid
  END INTERFACE read_onCells_2D

  INTERFACE read_onCells_2D_time
    MODULE PROCEDURE read_REAL_onCells_2D_time_streamid
  END INTERFACE read_onCells_2D_time

  INTERFACE read_onCells_2D_extdim
    MODULE PROCEDURE read_REAL_onCells_2D_1extdim_streamid
  END INTERFACE read_onCells_2D_extdim

  INTERFACE read_onCells_3D_time
    MODULE PROCEDURE read_REAL_onCells_3D_time_streamid
  END INTERFACE read_onCells_3D_time

  INTERFACE read_onCells_3D_extdim
    MODULE PROCEDURE read_REAL_onCells_3D_1extdim_streamid
  END INTERFACE read_onCells_3D_extdim

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS

  !-------------------------------------------------------------------------
  ! used for finding the names of the dimensions in the netcdf files
  CHARACTER(LEN=*), PARAMETER :: std_cells_dim_name_1 = 'cell'
  CHARACTER(LEN=*), PARAMETER :: std_cells_dim_name_2 = 'ncells'
  CHARACTER(LEN=*), PARAMETER :: std_time_dim_name_1  = 'time'

CONTAINS

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_0D_streamid(stream_id, variable_name, return_status)

    REAL(wp)    :: read_REAL_0D_streamid
    INTEGER, OPTIONAL :: return_status

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    read_REAL_0D_streamid = netcdf_read_0D_real(stream_id%file_id, variable_name)

  END FUNCTION read_REAL_0D_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_REAL_1D_streamid(stream_id, variable_name, fill_array, return_pointer)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)
    define_return_pointer        :: return_pointer(:)

    REAL(wp), POINTER            :: tmp_pointer(:)
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_1D_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_1D(stream_id%file_id, variable_name, fill_array)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_REAL_1D_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_REAL_2D_time_streamid(stream_id, variable_name, fill_array, return_pointer, &
    & dim_names, start_timestep, end_timestep)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_2D_time_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_2D_time(stream_id%file_id, variable_name, &
        & fill_array, dim_names, start_timestep, end_timestep)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_REAL_2D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_REAL_3D_time_streamid(stream_id, variable_name, fill_array, return_pointer, &
    & dim_names, start_timestep, end_timestep)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_3D_time_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_3D_time(stream_id%file_id, variable_name, &
        & fill_array, dim_names, start_timestep, end_timestep)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_REAL_3D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_REAL_onCells_2D_streamid(stream_id, variable_name, fill_array, return_pointer, patch)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    define_return_pointer        :: return_pointer(:,:)
    TYPE(t_patch), TARGET        :: patch

    REAL(wp), POINTER            :: tmp_pointer(:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_onCells_2D_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_2D(stream_id%file_id, variable_name, &
        &                           fill_array, patch%n_patch_cells_g, &
        &                           patch%cells%decomp_info%glb_index)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_REAL_onCells_2D_streamid
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_REAL_onCells_2D_time_streamid(stream_id, variable_name, fill_array, return_pointer, patch, &
    & start_timestep, end_timestep)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep

    CALL read_REAL_onCells_2D_1extdim_streamid(&
      & stream_id=stream_id, variable_name=variable_name, fill_array=fill_array, &
      & return_pointer=return_pointer, patch=patch, &
      & start_extdim=start_timestep, end_extdim=end_timestep, extdim_name="time" )

  END SUBROUTINE read_REAL_onCells_2D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_REAL_onCells_2D_1extdim_streamid(stream_id, variable_name, fill_array, return_pointer, patch, &
    & start_extdim, end_extdim, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_onCells_2D_1extdim_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_2D_extdim(stream_id%file_id, variable_name, fill_array, &
         & patch%n_patch_cells_g, patch%cells%decomp_info%glb_index, &
         & start_extdim, end_extdim, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_REAL_onCells_2D_1extdim_streamid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  SUBROUTINE read_REAL_onCells_3D_time_streamid(stream_id, variable_name, fill_array, return_pointer, patch, &
    & start_timestep, end_timestep, levelsDimName)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName

    CALL read_REAL_onCells_3D_1extdim_streamid(&
      & stream_id=stream_id,                    &
      & variable_name=variable_name,            &
      & fill_array=fill_array,                  &
      & return_pointer=return_pointer,          &
      & patch=patch,                            &
      & start_extdim=start_timestep,  end_extdim=end_timestep, &
      & levelsDimName=levelsDimName,          &
      & extdim_name="time")

  END SUBROUTINE read_REAL_onCells_3D_time_streamid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  SUBROUTINE read_REAL_onCells_3D_1extdim_streamid(stream_id, variable_name, fill_array, return_pointer, patch, &
    & start_extdim, end_extdim, levelsDimName, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsDimName

    REAL(wp), POINTER  :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_onCells_3D_1extdim_streamid'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_3D_extdim(stream_id%file_id, variable_name, &
         & fill_array, patch%n_patch_cells_g, patch%cells%decomp_info%glb_index, &
         & start_extdim, end_extdim, levelsDimName, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT


  END SUBROUTINE read_REAL_onCells_3D_1extdim_streamid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  TYPE(t_stream_id) FUNCTION openInputFile(filename, input_method)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, OPTIONAL, INTENT(IN) :: input_method
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:openInputFile'

    IF (present(input_method)) THEN
      openInputFile%input_method = input_method
    ELSE
      openInputFile%input_method = default_read_method
    ENDIF

    SELECT CASE(openInputFile%input_method)
    CASE (read_netcdf_broadcast_method)
      openInputFile%file_id = netcdf_open_input(filename)
    CASE (read_netcdf_distribute_method)
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
