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
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data, distrib_nf_open
  USE mo_model_domain, ONLY: t_patch


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

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => netcdf_read_2D(stream_id%file_id, variable_name, &
        &                           fill_array, &
        &                           stream_id%read_info(location)%n_g, &
        &                           stream_id%read_info(location)%glb_index)
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

    CALL read_dist_REAL_2D_1extdim_streamid(&
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
  SUBROUTINE read_dist_REAL_2D_1extdim_streamid(stream_id, location, &
    &                                           variable_name, fill_array, &
    &                                           return_pointer, start_extdim, &
    &                                           end_extdim, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
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
         & stream_id%read_info(location)%n_g, &
         & stream_id%read_info(location)%glb_index, start_extdim, end_extdim, &
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

    CALL read_dist_REAL_3D_1extdim_streamid( &
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
  SUBROUTINE read_dist_REAL_3D_1extdim_streamid(stream_id, location, &
    &                                           variable_name, fill_array, &
    &                                           return_pointer, start_extdim, &
    &                                           end_extdim, levelsDimName, &
    &                                           extdim_name )

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
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
         & fill_array, stream_id%read_info(location)%n_g, &
         & stream_id%read_info(location)%glb_index, &
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
  TYPE(t_stream_id) FUNCTION openInputFile_dist(filename, patch, input_method)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(t_patch), INTENT(INOUT) :: patch
    INTEGER, OPTIONAL, INTENT(IN) :: input_method

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile_dist'

    IF (PRESENT(input_method)) THEN
      openInputFile_dist%input_method = input_method
    ELSE
      openInputFile_dist%input_method = default_read_method
    END IF

    SELECT CASE(openInputFile_dist%input_method)
    CASE (read_netcdf_broadcast_method)

      openInputFile_dist%file_id = netcdf_open_input(filename)

      openInputFile_dist%read_info(onCells)%n_g = patch%n_patch_cells_g
      openInputFile_dist%read_info(onCells)%glb_index = &
        patch%cells%decomp_info%glb_index
      NULLIFY(openInputFile_dist%read_info(onCells)%dist_read_info)

      openInputFile_dist%read_info(onEdges)%n_g = patch%n_patch_edges_g
      openInputFile_dist%read_info(onEdges)%glb_index = &
        patch%edges%decomp_info%glb_index
      NULLIFY(openInputFile_dist%read_info(onEdges)%dist_read_info)

      openInputFile_dist%read_info(onVertices)%n_g = patch%n_patch_verts_g
      openInputFile_dist%read_info(onVertices)%glb_index = &
        patch%verts%decomp_info%glb_index
      NULLIFY(openInputFile_dist%read_info(onVertices)%dist_read_info)

    CASE (read_netcdf_distribute_method)

      CALL finish(method_name, "read_netcdf_distribute_method input_method" // &
        &         "not yet supported")

      openInputFile_dist%file_id = distrib_nf_open(TRIM(filename))

      ! this is not yet part of t_patch
      ! openInputFile_dist%read_info(onCells)%dist_read_info = &
      !  patch%cells%dist_read_info
      openInputFile_dist%read_info(onCells)%n_g = -1
      NULLIFY(openInputFile_dist%read_info(onCells)%glb_index)

      ! this is not yet part of t_patch
      ! openInputFile_dist%read_info(onVertices)%dist_read_info = &
      !  patch%verts%dist_read_info
      openInputFile_dist%read_info(onVertices)%n_g = -1
      NULLIFY(openInputFile_dist%read_info(onVertices)%glb_index)

      ! this is not yet part of t_patch
      ! openInputFile_dist%read_info(onEdges)%dist_read_info = &
      !  patch%edges%dist_read_info
      openInputFile_dist%read_info(onEdges)%n_g = -1
      NULLIFY(openInputFile_dist%read_info(onEdges)%glb_index)

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


END MODULE mo_read_interface
