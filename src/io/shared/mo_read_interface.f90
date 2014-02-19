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
!! Moved from mo_util_netcdf, added read_oncells, by L. Linardakis, (2013-03-15)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!#define define_fill_target REAL(wp), TARGET, ALLOCATABLE, OPTIONAL
!#define define_fill_target REAL(wp), POINTER, OPTIONAL
#define define_fill_target REAL(wp), TARGET, OPTIONAL

MODULE mo_read_interface

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

  USE mo_io_config,         ONLY:  read_netcdf_broadcast_method, &
    & read_netcdf_distribute_method,  default_read_method

  USE mo_netcdf_read,      ONLY: netcdf_open_input, netcdf_close, &
    & netcdf_read_0D_real, netcdf_read_1D, netcdf_read_2D_time, netcdf_read_3D_time, &
    & netcdf_read_oncells_2D, netcdf_read_oncells_2D_extdim, netcdf_read_oncells_3D_extdim


  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_netcdf_broadcast_method, read_netcdf_distribute_method

  PUBLIC :: openInputFile, closeFile

  PUBLIC :: read_0D_real
  PUBLIC :: read_1D
  PUBLIC :: read_oncells_2D
  PUBLIC :: read_oncells_2D_time
  PUBLIC :: read_oncells_3D_time
  PUBLIC :: read_oncells_2D_extdim
  PUBLIC :: read_oncells_3D_extdim

  !--------------------------------------------------------
  INTERFACE read_0D_real
    MODULE PROCEDURE read_REAL_0D_filename
    MODULE PROCEDURE read_REAL_0D_fileid
  END INTERFACE read_0D_real

  INTERFACE read_1D
    MODULE PROCEDURE read_REAL_1D_filename
    MODULE PROCEDURE read_REAL_1D_fileid
  END INTERFACE read_1D

  INTERFACE read_oncells_2D
    MODULE PROCEDURE read_REAL_ONCELLS_2D_filename
    MODULE PROCEDURE read_REAL_ONCELLS_2D_fileid
  END INTERFACE read_oncells_2D

  INTERFACE read_oncells_2D_time
    MODULE PROCEDURE read_REAL_ONCELLS_2D_time_filename
    MODULE PROCEDURE read_REAL_ONCELLS_2D_time_fileid
  END INTERFACE read_oncells_2D_time

  INTERFACE read_oncells_2D_extdim
    MODULE PROCEDURE read_REAL_ONCELLS_2D_1extdim_filename
    MODULE PROCEDURE read_REAL_ONCELLS_2D_1extdim_fileid
  END INTERFACE read_oncells_2D_extdim

  INTERFACE read_oncells_3D_time
    MODULE PROCEDURE read_REAL_ONCELLS_3D_time_filename
    MODULE PROCEDURE read_REAL_ONCELLS_3D_time_fileid
  END INTERFACE read_oncells_3D_time

  INTERFACE read_oncells_3D_extdim
    MODULE PROCEDURE read_REAL_ONCELLS_3D_1extdim_filename
    MODULE PROCEDURE read_REAL_ONCELLS_3D_1extdim_fileid
  END INTERFACE read_oncells_3D_extdim

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS

  !-------------------------------------------------------------------------
  ! used for finding the names of the dimensions in the netcdf files
  CHARACTER(LEN=*), PARAMETER :: std_cells_dim_name_1 = 'cell'
  CHARACTER(LEN=*), PARAMETER :: std_cells_dim_name_2 = 'ncells'
  CHARACTER(LEN=*), PARAMETER :: std_time_dim_name_1  = 'time'

  INTEGER :: current_read_method

CONTAINS

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_1D_filename(filename, variable_name, fill_array, read_method)
    REAL(wp), POINTER :: read_REAL_1D_filename(:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)
    INTEGER, OPTIONAL :: read_method

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_2D_filename'

    file_id = openInputFile(filename, read_method)
    read_REAL_1D_filename => &
      & read_REAL_1D_fileid( &
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array)
    return_status = closeFile(file_id)

  END FUNCTION read_REAL_1D_filename
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_ONCELLS_2D_filename(filename, variable_name, fill_array, patch, read_method)
    REAL(wp), POINTER :: read_REAL_ONCELLS_2D_filename(:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, OPTIONAL :: read_method

    INTEGER :: file_id
    INTEGER :: return_status    
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_2D_filename'

    file_id = openInputFile(filename, read_method)
    read_REAL_ONCELLS_2D_filename => &
      & read_REAL_ONCELLS_2D_fileid( &
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array, patch=patch)
    return_status = closeFile(file_id)
                              
  END FUNCTION read_REAL_ONCELLS_2D_filename
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_ONCELLS_2D_time_filename(filename, variable_name, fill_array, patch, &
    & start_timestep, end_timestep, read_method )

    REAL(wp), POINTER            :: read_REAL_ONCELLS_2D_time_filename(:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    INTEGER, OPTIONAL :: read_method

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_2D_time_filename'

    file_id = openInputFile(filename, read_method)
    read_REAL_ONCELLS_2D_time_filename =>   &
      & read_REAL_ONCELLS_2D_time_fileid(file_id=file_id, variable_name=variable_name, &
      & fill_array=fill_array, patch=patch, start_timestep=start_timestep, end_timestep=end_timestep)
    return_status = closeFile(file_id)

  END FUNCTION read_REAL_ONCELLS_2D_time_filename
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_ONCELLS_2D_1extdim_filename(filename, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, extdim_name, read_method )

    REAL(wp), POINTER            :: read_REAL_ONCELLS_2D_1extdim_filename(:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name
    INTEGER, OPTIONAL :: read_method

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_2D_1extdim_filename'

    file_id = openInputFile(filename, read_method)
    read_REAL_ONCELLS_2D_1extdim_filename => &
      & read_REAL_ONCELLS_2D_1extdim_fileid(file_id=file_id, variable_name=variable_name, &
      & fill_array=fill_array, patch=patch, &
      & start_extdim=start_extdim, end_extdim=end_extdim, extdim_name=extdim_name )
    return_status = closeFile(file_id)

  END FUNCTION read_REAL_ONCELLS_2D_1extdim_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_ONCELLS_3D_time_filename(filename, variable_name, fill_array, patch, &
    & start_timestep, end_timestep, levelsDimName, read_method )

    REAL(wp), POINTER :: read_REAL_ONCELLS_3D_time_filename(:,:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName
    INTEGER, OPTIONAL :: read_method

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_3D_time_filename'

    file_id = openInputFile(filename, read_method)
    read_REAL_ONCELLS_3D_time_filename => &
      & read_REAL_ONCELLS_3D_time_fileid( &
      & file_id=file_id,                         &
      & variable_name=variable_name,             &
      & fill_array=fill_array,                   &
      & patch=patch,                             &
      & levelsDimName=levelsDimName,           &
      & start_timestep=start_timestep,  end_timestep=end_timestep)
    return_status = closeFile(file_id)

  END FUNCTION read_REAL_ONCELLS_3D_time_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_ONCELLS_3D_1extdim_filename(filename, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, levelsDimName, extdim_name, read_method )

    REAL(wp), POINTER :: read_REAL_ONCELLS_3D_1extdim_filename(:,:,:,:)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target             :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName
    INTEGER, OPTIONAL :: read_method

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_3D_1extdim_filename'

    file_id = openInputFile(filename, read_method)
    read_REAL_ONCELLS_3D_1extdim_filename => read_REAL_ONCELLS_3D_1extdim_fileid( &
      & file_id=file_id,                          &
      & variable_name=variable_name,              &
      & fill_array=fill_array,                    &
      &  patch=patch,                             &
      & start_extdim=start_extdim,  end_extdim=end_extdim, &
      & levelsDimName=levelsDimName,            &
      & extdim_name=extdim_name)
    return_status = closeFile(file_id)

  END FUNCTION read_REAL_ONCELLS_3D_1extdim_filename
  !-------------------------------------------------------------------------
  FUNCTION read_REAL_0D_filename(filename, variable_name)

    REAL(wp), POINTER            :: read_REAL_0D_filename

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    INTEGER                      :: file_id
    INTEGER :: return_status

    file_id = openInputFile(filename, read_method=read_netcdf_broadcast_method)
    read_REAL_0D_filename = read_0D_real(file_id, variable_name)
    return_status = closeFile(file_id)

  END FUNCTION read_REAL_0D_filename
  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_0D_fileid(file_id, variable_name)

    REAL(wp)            :: read_REAL_0D_fileid

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    read_REAL_0D_fileid = netcdf_read_0D_real(file_id, variable_name)

  END FUNCTION read_REAL_0D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_1D_fileid(file_id, variable_name, fill_array)

    REAL(wp), POINTER            :: read_REAL_1D_fileid(:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_1D_fileid'


    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      read_REAL_1D_fileid => netcdf_read_1D(file_id, variable_name, fill_array)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT

  END FUNCTION read_REAL_1D_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_2D_time_fileid(file_id, variable_name, fill_array, dim_names, start_timestep, end_timestep)

    REAL(wp), POINTER            :: read_REAL_2D_time_fileid(:,:,:)

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

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_2D_time_fileid'

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      read_REAL_2D_time_fileid => netcdf_read_2D_time(file_id, variable_name, fill_array, dim_names, start_timestep, end_timestep)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT

  END FUNCTION read_REAL_2D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_3D_time_fileid(file_id, variable_name, fill_array, dim_names, start_timestep, end_timestep)

    REAL(wp), POINTER            :: read_REAL_3D_time_fileid(:,:,:,:)

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

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_3D_time_fileid'

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      read_REAL_3D_time_fileid => netcdf_read_3D_time(file_id, variable_name, fill_array, dim_names, start_timestep, end_timestep)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT

  END FUNCTION read_REAL_3D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION read_REAL_ONCELLS_2D_fileid(file_id, variable_name, fill_array, patch)

    REAL(wp), POINTER            :: read_REAL_ONCELLS_2D_fileid(:,:)

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

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_2D_fileid'

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      read_REAL_ONCELLS_2D_fileid => netcdf_read_oncells_2D(file_id, variable_name, fill_array, patch)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT

  END FUNCTION read_REAL_ONCELLS_2D_fileid
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION read_REAL_ONCELLS_2D_time_fileid(file_id, variable_name, fill_array, patch, &
    & start_timestep, end_timestep)
    REAL(wp), POINTER            :: read_REAL_ONCELLS_2D_time_fileid(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep

    read_REAL_ONCELLS_2D_time_fileid => read_REAL_ONCELLS_2D_1extdim_fileid(&
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array, patch=patch, &
      & start_extdim=start_timestep, end_extdim=end_timestep, extdim_name="time" )

  END FUNCTION read_REAL_ONCELLS_2D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION read_REAL_ONCELLS_2D_1extdim_fileid(file_id, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, extdim_name )

    REAL(wp), POINTER            :: read_REAL_ONCELLS_2D_1extdim_fileid(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target   :: fill_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_2D_1extdim_fileid'

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      read_REAL_ONCELLS_2D_1extdim_fileid => &
         & netcdf_read_oncells_2D_extdim(file_id, variable_name, fill_array, patch, &
         & start_extdim, end_extdim, extdim_name )
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT

  END FUNCTION read_REAL_ONCELLS_2D_1extdim_fileid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION read_REAL_ONCELLS_3D_time_fileid(file_id, variable_name, fill_array, patch, &
    & start_timestep, end_timestep, levelsDimName)

    REAL(wp), POINTER  :: read_REAL_ONCELLS_3D_time_fileid(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target            :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName

    read_REAL_ONCELLS_3D_time_fileid => read_REAL_ONCELLS_3D_1extdim_fileid(&
      & file_id=file_id,                        &
      & variable_name=variable_name,            &
      & fill_array=fill_array,                  &
      & patch=patch,                            &
      & start_extdim=start_timestep,  end_extdim=end_timestep, &
      & levelsDimName=levelsDimName,          &
      & extdim_name="time")

  END FUNCTION read_REAL_ONCELLS_3D_time_fileid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION read_REAL_ONCELLS_3D_1extdim_fileid(file_id, variable_name, fill_array, patch, &
    & start_extdim, end_extdim, levelsDimName, extdim_name )

    REAL(wp), POINTER  :: read_REAL_ONCELLS_3D_1extdim_fileid(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch
    INTEGER, INTENT(in), OPTIONAL:: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsDimName

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

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:read_REAL_ONCELLS_3D_1extdim_fileid'

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      read_REAL_ONCELLS_3D_1extdim_fileid => &
         & netcdf_read_oncells_3D_extdim(file_id, variable_name, fill_array, patch, &
         & start_extdim, end_extdim, levelsDimName, extdim_name )
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT


  END FUNCTION read_REAL_ONCELLS_3D_1extdim_fileid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION openInputFile(filename, read_method)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, OPTIONAL, INTENT(IN) :: read_method
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:openInputFile'

    IF (present(read_method)) THEN
      current_read_method = read_method
    ELSE
      current_read_method = default_read_method
    ENDIF

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      openInputFile = netcdf_open_input(filename)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT
    
  END FUNCTION openInputFile
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION closeFile(file_id)
    INTEGER, INTENT(IN) :: file_id

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_interface:closeFile'

    SELECT CASE(current_read_method)
    CASE (read_netcdf_broadcast_method)
      closeFile = netcdf_close(file_id)
    CASE (read_netcdf_distribute_method)
    CASE default
      CALL finish(method_name, "unknown read_method")
    END SELECT

  END FUNCTION closeFile
  !-------------------------------------------------------------------------
  

END MODULE mo_read_interface
