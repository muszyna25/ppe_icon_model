!>
!! Module handling the transfer of variable meta info to asynchronous I/O PEs.
!!
!! @author F. Prill
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2014-08-04)
!!
MODULE mo_name_list_output_metadata

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_f_pointer
  USE mo_exception,                         ONLY: finish
  USE mo_var_metadata_types,                ONLY: t_var_metadata
  USE mo_name_list_output_types,            ONLY: t_mem_win
  USE mo_mpi,                               ONLY: p_int, p_comm_work_io,      &
    &                                             p_comm_rank
  USE mo_dynamics_config,                   ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_impl_constants,                    ONLY: TLEV_NNOW, TLEV_NNOW_RCF, TLEV_NNEW, TLEV_NNEW_RCF

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: metainfo_allocate_memory_window
  PUBLIC :: metainfo_write_to_memwin
  PUBLIC :: metainfo_get_from_buffer
  PUBLIC :: metainfo_get_size
  PUBLIC :: metainfo_get_timelevel

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_metadata'

CONTAINS

  !-------------------------------------------------------------------------------------------------
  !> @return size of a single variable's info object
  !
  !  @author F. Prill, DWD
  !
  FUNCTION metainfo_get_size() RESULT(info_size)
    INTEGER :: info_size
    ! local variables
    TYPE(t_var_metadata) :: info  ! dummy meta data object

    info_size = SIZE(TRANSFER(info, (/ 0 /)))
  END FUNCTION metainfo_get_size


  !-------------------------------------------------------------------------------------------------
  !> Allocates an MPI memory window for the meta info of the variables fields.
  !   - allocation for asynchronous I/O mode only.
  !   - allocation on I/O PEs and PE #0 only.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE metainfo_allocate_memory_window(memwin, nvars)

#ifndef NOMPI
#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL
#endif
#endif

    TYPE(t_mem_win),      INTENT(INOUT) :: memwin ! MPI memory window
    INTEGER,              INTENT(IN)    :: nvars  ! total no. of variables

#ifndef NOMPI

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::metainfo_allocate_memory_window"
    INTEGER                         :: ierror, comm_rank, mem_size
    INTEGER, TARGET                 :: dummy(1, 1)
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes, nbytes_int, lb
    TYPE(c_ptr)                     :: c_mem_ptr

    ! total number of required integer variables
    mem_size = metainfo_get_size()
    ! Get amount of bytes per INTEGER variable (in MPI communication)
    CALL MPI_Type_get_extent(p_int, lb, nbytes_int, ierror)
    IF (ierror /= 0) CALL finish(routine, "MPI error!")
    mem_bytes = INT(nvars, mpi_address_kind) * INT(mem_size, mpi_address_kind) &
         * nbytes_int


    comm_rank = p_comm_rank(p_comm_work_io)

    IF (comm_rank == 0) THEN
      CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, ierror)
      IF (ierror /= 0) CALL finish(routine, "MPI error!")
      CALL C_F_POINTER(c_mem_ptr, memwin%mem_ptr_metainfo_pe0, &
           (/ mem_size, nvars /) )
      memwin%mem_ptr_metainfo_pe0 = 0
    ELSE
      mem_bytes = 0_mpi_address_kind
      memwin%mem_ptr_metainfo_pe0 => dummy
    END IF

    ! Create memory window for meta-data communication
    CALL MPI_Win_create(memwin%mem_ptr_metainfo_pe0, mem_bytes, &
      &                 INT(nbytes_int), mpi_info_null, p_comm_work_io, &
      &                 memwin%mpi_win_metainfo, ierror)
    IF (ierror /= 0) CALL finish(routine, "MPI error!")
    IF (comm_rank /= 0) NULLIFY(memwin%mem_ptr_metainfo_pe0)
#endif
  END SUBROUTINE metainfo_allocate_memory_window


  !-------------------------------------------------------------------------------------------------
  !> Store a variable's meta-info to a memory window.
  !  This subroutine does nothing on PEs except compute PE #0 or
  !  if we are not running in asynchronous I/O mode.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE metainfo_write_to_memwin(memwin, ivar, info)
    TYPE(t_mem_win),      INTENT(INOUT) :: memwin ! MPI memory window
    INTEGER,              INTENT(IN)    :: ivar   ! index of variable (corresponds to data memwin)
    TYPE(t_var_metadata), INTENT(IN)    :: info   ! meta data for variable
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::metainfo_write_to_memwin"
    IF (.NOT. ASSOCIATED(memwin%mem_ptr_metainfo_pe0)) THEN
      CALL finish(routine, "Internal error!")
    END IF

    ! copy the info object into the memory window
    memwin%mem_ptr_metainfo_pe0(:, ivar) =  TRANSFER(info, (/ 0 /))
  END SUBROUTINE metainfo_write_to_memwin


  !-------------------------------------------------------------------------------------------------
  !> Retrieve a variable's meta-info from memory window.
  !  This subroutine does nothing on compute PEs or
  !  if we are not running in asynchronous I/O mode.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE metainfo_get_from_buffer(buf, info)
    INTEGER, INTENT(IN)    :: buf(:)
    TYPE(t_var_metadata), INTENT(OUT) :: info   ! meta data for variable
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::metainfo_get_from_memwin"
    ! copy the info object from the memory window
    info = TRANSFER(buf, info)
  END SUBROUTINE metainfo_get_from_buffer

  !-------------------------------------------------------------------------------------------------
  !> Return the output timelevel for given variables info object
  !
  !  @author R. Mueller, MPI
  !
  INTEGER FUNCTION metainfo_get_timelevel(info,domain) RESULT(timelevel)
    TYPE(t_var_metadata), INTENT(IN) :: info
    INTEGER, INTENT(IN)              :: domain

    CHARACTER(LEN=*), PARAMETER       :: routine = modname//"::metainfo_get_timelevel"

    SELECT CASE (info%tlev_source)
    CASE(TLEV_NNOW);     timelevel = nnow(domain)
    CASE(TLEV_NNOW_RCF); timelevel = nnow_rcf(domain)
    CASE(TLEV_NNEW);     timelevel = nnew(domain)
    CASE(TLEV_NNEW_RCF); timelevel = nnew_rcf(domain)
    CASE DEFAULT
      CALL finish(routine,'Unsupported tlev_source')
    END SELECT
  END FUNCTION metainfo_get_timelevel
END MODULE mo_name_list_output_metadata
