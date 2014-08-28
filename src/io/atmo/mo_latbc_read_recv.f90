!>
!! This module provides basic method to send
!! and receive file data using asynchronous 
!! communication.
!!
!! @author M. Pondkule (DWD)
!!
!! @par Revision History
!! Initial release by M. Pondkule, DWD (2013-04-30)
!! 
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
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
  ! TODO [MP]
  !
  ! NAG Fortran Compiler Release 5.3.2(951)
  ! Warning: ../../../src/io/shared/mo_latbc_read_recv.f90, line 241: Unused local variable I_DOM
  ! Warning: ../../../src/io/shared/mo_latbc_read_recv.f90, line 241: Unused dummy variable VARNAME

MODULE mo_latbc_read_recv

#ifndef NOMPI
  USE mpi
#endif

  USE mo_kind,               ONLY: sp, dp, i8
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_mpi,                ONLY: p_pe_work,  &
    &                              num_work_procs, p_real_sp               
  USE mo_util_cdi,           ONLY: get_cdi_varID
  USE mo_async_latbc_types,  ONLY: t_patch_data, t_reorder_data
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
  
  IMPLICIT NONE
  !-----------------------------------------------------------------
  ! include NetCDF headers (direct NetCDF library calls are required
  ! for output of grid information).
  !-----------------------------------------------------------------
  INCLUDE 'cdi.inc'

  PRIVATE

  PUBLIC  :: prefetch_cdi_2d 
  PUBLIC  :: prefetch_cdi_3d
  PUBLIC  :: compute_data_receive

  CHARACTER(len=*), PARAMETER :: version = &
    &    '$Id: mo_cdi_prefetch.f90 16829 2014-04-30 14:27:20Z mukund.pondkule $'
  CHARACTER(LEN=*), PARAMETER :: modname = '::mo_latbc_read_recv'

CONTAINS

  !-------------------------------------------------------------------------
  !> Read 3D dataset from file.
  ! 
  !  Note: This implementation uses a 2D buffer.
  ! 
  !  @par Revision History
  !  Initial revision by M. Pondkule, DWD (2014-05-19)
  ! 
  SUBROUTINE prefetch_cdi_3d(streamID, varname, patch_data, nlevs, hgrid, ioff)
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: ioff(0:)
#else
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
#endif
    INTEGER,             INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*),    INTENT(IN)    :: varname        !< Var name of field to be read
    TYPE(t_patch_data),  INTENT(IN)    :: patch_data     !< patch data containing information for prefetch 
    INTEGER,             INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    INTEGER,             INTENT(IN)    :: hgrid          !< stored variable location indication
    
    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':prefetch_cdi_3d'
    ! local variables:
    INTEGER                         :: vlistID, varID, zaxisID, gridID,   &
      &                                jk, ierrstat, dimlen(3), nmiss
    REAL(dp), ALLOCATABLE           :: tmp_buf(:) ! temporary local array
    
#ifndef NOMPI
    ! allocate a buffer for one vertical level
    ALLOCATE(tmp_buf(patch_data%n_patch_cells_g), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! initialize temporary buffers
    tmp_buf(:) = 0._dp

    ! get var ID
    vlistID   = streamInqVlist(streamID)
    varID     = get_cdi_varID(streamID, name=TRIM(varname)) 
    zaxisID   = vlistInqVarZaxis(vlistID, varID)
    gridID    = vlistInqVarGrid(vlistID, varID)
    dimlen(1) = gridInqSize(gridID)
    dimlen(2) = zaxisInqSize(zaxisID)

    ! Check variable dimensions:
    IF ((dimlen(1) /= patch_data%n_patch_cells_g) .OR.  &
         & (dimlen(2) /= nlevs)) THEN
       WRITE(message_text,'(s,2i4)') "Horizontal cells: ", dimlen(1), patch_data%n_patch_cells_g, &
         &                           "nlev: ", dimlen(2), nlevs
       CALL message(TRIM(routine), message_text)
       CALL finish(routine, "Incompatible dimensions!")
    END IF
   
    DO jk=1, nlevs
       ! read record as 1D field
       CALL streamReadVarSlice(streamID, varID, jk-1, tmp_buf(:), nmiss)
       ! send 2d buffer using MPI_PUT
       CALL prefetch_proc_send(patch_data, tmp_buf(:), 1, hgrid, ioff)
       tmp_buf(:) = 0._dp 
    ENDDO ! jk=1,nlevs 
    !  ENDIF
  
    ! clean up
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#endif
  
  END SUBROUTINE prefetch_cdi_3d

  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  !  Initial revision by M. Pondkule, DWD (2014-05-15) 
  !
  SUBROUTINE prefetch_cdi_2d (streamID, varname, patch_data, hgrid, ioff)
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: ioff(0:)
#else
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
#endif
    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Varname of field to be read
    TYPE(t_patch_data), INTENT(IN)  :: patch_data  !< patch data containing information for prefetch 
    INTEGER,          INTENT(IN)    :: hgrid          !< stored variable location indication
    
    ! local variables:
    CHARACTER(len=max_char_length), PARAMETER :: &
         routine = modname//':prefetch_cdi_2d_real_id'
    INTEGER       :: varID, nmiss, gridID, vlistID, dimlen(1)
    REAL(dp), ALLOCATABLE :: z_dummy_array(:)       !< local dummy array
  
#ifndef NOMPI
    ! get var ID
    vlistID   = streamInqVlist(streamID) 
    varID     = get_cdi_varID(streamID, name=TRIM(varname))
    gridID    = vlistInqVarGrid(vlistID, varID)
    dimlen(1) = gridInqSize(gridID)
    ALLOCATE(z_dummy_array(patch_data%n_patch_cells_g))

    ! Check variable dimensions:
    IF (dimlen(1) /= patch_data%n_patch_cells_g) THEN
       WRITE(message_text,'(s,2i4)') "Horizontal cells: ", dimlen(1), patch_data%n_patch_cells_g
       CALL message(TRIM(routine), message_text)
       CALL finish(routine, "Incompatible dimensions!")
    END IF
    
    ! prefetch PE reads and puts data in memory window
    ! read record as 1D field
    CALL streamReadVarSlice(streamID, varID, 0, z_dummy_array(:), nmiss)
    CALL prefetch_proc_send(patch_data, z_dummy_array(:), 1, hgrid, ioff)

    DEALLOCATE(z_dummy_array)
#endif

  END SUBROUTINE prefetch_cdi_2d

  !-------------------------------------------------------------------------
  !> compute processor copy 2D or 3D dataset from memory window buffer.
  ! 
  !  @par Revision History
  !  Initial revision by M. Pondkule, DWD (2014-05-15)
  ! 
  SUBROUTINE compute_data_receive (varname, hgrid, nlevs, var_out, eoff, patch_data)
 
    CHARACTER(len=*),    INTENT(IN)    :: varname        !< var name of field to be read
    INTEGER,             INTENT(IN)    :: hgrid          !< stored variable location indication
    INTEGER,             INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(sp),            INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER(i8),         INTENT(INOUT) :: eoff
    TYPE(t_patch_data),  TARGET, INTENT(IN)   :: patch_data
    TYPE(t_reorder_data), POINTER             :: p_ri
    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
         routine = modname//':compute_data_receive'
    ! local variables:
    INTEGER     :: j, jl, jb, jk, i_dom, mpi_error      

#ifndef NOMPI
    CALL MPI_Win_lock(MPI_LOCK_SHARED, p_pe_work, MPI_MODE_NOCHECK, patch_data%mem_win%mpi_win, mpi_error)

    ! initialize output field:
    var_out(:,:,:) = 0._sp

    ! Get pointer to appropriate reorder_info
    SELECT CASE (hgrid)
    CASE(GRID_UNSTRUCTURED_CELL)
       p_ri => patch_data%cells
    CASE(GRID_UNSTRUCTURED_EDGE)
       p_ri => patch_data%edges
    CASE default
       CALL finish(routine,'unknown grid type')
    END SELECT

    DO jk=1, nlevs
       DO j = 1,  p_ri%n_own 
          jb = p_ri%own_blk(j) ! Block index in distributed patch
          jl = p_ri%own_idx(j) ! Line  index in distributed patch
          var_out(jl,jk,jb) = REAL(patch_data%mem_win%mem_ptr_sp(eoff+INT(j,i8)),sp)
    ENDDO
       eoff = eoff + INT(p_ri%n_own,i8) 
    END DO ! jk=1,nlevs

    CALL MPI_Win_unlock(p_pe_work, patch_data%mem_win%mpi_win, mpi_error)
#endif

  END SUBROUTINE compute_data_receive
 
  !-------------------------------------------------------------------------
  !> routine on the prefetch PE to write variable values to memory window
  !
  !  @note This subroutine is called by prefetch PE only.
  !  Initial revision by M. Pondkule, DWD (2014-05-27) 
  !
  SUBROUTINE prefetch_proc_send(patch_data, var1_dp, nlevs, hgrid, ioff)
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: ioff(0:)
#else
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
#endif
    ! local variables  
    TYPE(t_patch_data), TARGET, INTENT(IN)  :: patch_data  !< patch data containing information for prefetch 
    REAL(dp), INTENT(IN) :: var1_dp(:)
    INTEGER, INTENT(IN) :: nlevs
    INTEGER, INTENT(IN) :: hgrid 
    INTEGER :: voff(0:num_work_procs-1) 
    INTEGER :: nv_off_np(0:num_work_procs)
    TYPE(t_reorder_data),  POINTER :: p_ri
    REAL(sp), ALLOCATABLE   :: var3_sp(:)
    INTEGER :: np, nval, nv_off, mpi_error, ierrstat !gridType, gridNumber
    INTEGER :: dst_start, dst_end, src_start, src_end
    CHARACTER(len=max_char_length), PARAMETER :: routine = modname//'::prefetch_proc_send'

    SELECT CASE (hgrid)
       ! 1 is equal to GRID_UNSTRUCTURED_CELL 
    CASE(GRID_UNSTRUCTURED_CELL)
       p_ri => patch_data%cells
       ! 3 is equal to GRID_UNSTRUCTURED_EDGE
    CASE(GRID_UNSTRUCTURED_EDGE)
       p_ri => patch_data%edges
    CASE default
       CALL finish(routine,'unknown grid type')
    END SELECT
    
    ! compute the total offset for each PE
    nv_off       = 0
    nv_off_np(0) = 0
    DO np = 0, num_work_procs-1
       voff(np)        = nv_off
       nval            = p_ri%pe_own(np) * nlevs
       nv_off          = nv_off + nval
       nv_off_np(np+1) = nv_off_np(np) + p_ri%pe_own(np)
    END DO

    ALLOCATE(var3_sp(p_ri%n_glb), STAT=ierrstat) ! Must be allocated to exact size
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      
!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
          DO np = 0, num_work_procs-1
            dst_start = nv_off_np(np)+1
            dst_end   = nv_off_np(np+1)
            src_start = voff(np)+1
            src_end   = voff(np)+p_ri%pe_own(np)
            voff(np)  = src_end
            var3_sp(p_ri%reorder_index(dst_start:dst_end)) = var1_dp(src_start:src_end)
          ENDDO
!$OMP END DO
!$OMP END PARALLEL

    nv_off  = 0_i8
    nval = 0

#ifndef NOMPI
    DO np = 0, num_work_procs-1

       IF(p_ri%pe_own(np) == 0) CYCLE

       nval = p_ri%pe_own(np)*nlevs

       CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, np, MPI_MODE_NOCHECK, patch_data%mem_win%mpi_win, mpi_error)
    
       CALL MPI_PUT(var3_sp(nv_off+1), nval, p_real_sp, np, ioff(np), nval, p_real_sp, &
            & patch_data%mem_win%mpi_win, mpi_error) !MPI_WIN_NULL) 

       CALL MPI_Win_unlock(np, patch_data%mem_win%mpi_win, mpi_error)

       ! Update the offset in var1
       nv_off = nv_off + nval

       ! Update the offset in the memory window on compute PEs
       ioff(np) = ioff(np) + INT(nval,i8) ! replace with mpi_kind

    ENDDO
#endif

    DEALLOCATE(var3_sp, STAT=ierrstat) ! Must be allocated to exact size
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE prefetch_proc_send

  !------------------------------------------------------------------------------------------------

END MODULE mo_latbc_read_recv