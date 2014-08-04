!>
!! Module handling the transfer of variable meta info to asynchronous I/O PEs.
!!
!! @author F. Prill
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2014-08-04)
!!
MODULE mo_name_list_output_metadata

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_f_pointer, c_intptr_t
#endif

  USE mo_exception,                         ONLY: finish
  USE mo_kind,                              ONLY: i8
  USE mo_var_metadata_types,                ONLY: t_var_metadata
  USE mo_name_list_output_types,            ONLY: t_mem_win
  USE mo_mpi,                               ONLY: p_int, p_comm_work_io,      &
    &                                             my_process_is_mpi_test,     &
    &                                             my_process_is_mpi_workroot, &
    &                                             my_process_is_io
  USE mo_name_list_output_config,           ONLY: use_async_name_list_io

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: metainfo_allocate_memory_window
  PUBLIC :: metainfo_write_to_memwin
  PUBLIC :: metainfo_get_from_memwin
  PUBLIC :: metainfo_get_size

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

    TYPE(t_mem_win),      INTENT(INOUT) :: memwin ! MPI memory window
    INTEGER,              INTENT(IN)    :: nvars  ! total no. of variables

#ifndef NOMPI
#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL
#endif

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::metainfo_allocate_memory_window"
    INTEGER                         :: nbytes_int, mpierr
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size, mem_bytes
    TYPE(c_ptr)                     :: c_mem_ptr

#ifdef USE_CRAY_POINTER
    CALL finish(routine, "Cray pointers: not yet implemented!")
#endif
#ifdef __SX__
    CALL finish(routine, "SX no longer supported!")
#endif

    ! total number of required integer variables
    mem_size = nvars * metainfo_get_size()
    ! Get amount of bytes per INTEGER variable (in MPI communication)
    CALL MPI_Type_extent(p_int, nbytes_int, mpierr)
    mem_bytes = MAX(mem_size, 1_i8)*INT(nbytes_int,i8)

    ! TYPE(c_ptr) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
    ! So check if at least c_intptr_t and MPI_ADDRESS_KIND are the same, else we may get
    ! into deep, deep troubles!
    ! There is still a slight probability that TYPE(c_ptr) does not have the size indicated
    ! by c_intptr_t since the standard only requires c_intptr_t is big enough to hold pointers
    ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
    ! in such a way!!!
    ! If c_intptr_t<=0, this type is not defined and we can't do this check, of course.

    IF(c_intptr_t > 0 .AND. c_intptr_t /= MPI_ADDRESS_KIND) &
     & CALL finish(routine,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpierr)

    NULLIFY(memwin%mem_ptr_metainfo_pe0)
    CALL C_F_POINTER(c_mem_ptr, memwin%mem_ptr_metainfo_pe0, (/ mem_size /) )

    ! Create memory window for meta-data communication
    memwin%mem_ptr_metainfo_pe0(:) = 0
    CALL MPI_Win_create( memwin%mem_ptr_metainfo_pe0, mem_bytes, nbytes_int, MPI_INFO_NULL, &
      &                  p_comm_work_io, memwin%mpi_win_metainfo, mpierr )
    IF (mpierr /= 0) CALL finish(routine, "MPI error!")
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
    INTEGER :: info_size, offset

    ! This subroutine does nothing on PEs except compute PE #0 or if
    !  we are not running in asynchronous I/O mode.
    IF (.NOT. use_async_name_list_io .OR. my_process_is_mpi_test()) RETURN
    IF (.NOT. my_process_is_mpi_workroot())                         RETURN
    
    IF (.NOT. ASSOCIATED(memwin%mem_ptr_metainfo_pe0)) THEN
      CALL finish(routine, "Internal error!")
    END IF

    ! get size of a single meta-info field
    info_size = metainfo_get_size()
    offset    = (ivar-1)*info_size
    ! copy the info object into the memory window
    memwin%mem_ptr_metainfo_pe0((offset+1):(offset+info_size)) =  TRANSFER(info, (/ 0 /))
  END SUBROUTINE metainfo_write_to_memwin


  !-------------------------------------------------------------------------------------------------
  !> Retrieve a variable's meta-info from memory window.
  !  This subroutine does nothing on compute PEs or
  !  if we are not running in asynchronous I/O mode.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE metainfo_get_from_memwin(bufr_metainfo, ivar, info)
    INTEGER,              INTENT(IN)    :: bufr_metainfo(:) ! MPI memory window
    INTEGER,              INTENT(IN)    :: ivar   ! index of variable (corresponds to data memwin)
    TYPE(t_var_metadata), INTENT(INOUT) :: info   ! meta data for variable
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::metainfo_get_from_memwin"
    TYPE(t_var_metadata) :: dummy_info  ! dummy meta data object
    INTEGER              :: info_size, offset

    ! This subroutine does nothing on PEs except compute PE #0 or if
    !  we are not running in asynchronous I/O mode.
    IF (.NOT. use_async_name_list_io .OR. my_process_is_mpi_test()) RETURN
    IF (.NOT. my_process_is_io())                                   RETURN

    ! get size of a single meta-info field
    info_size = metainfo_get_size()
    offset    = (ivar-1)*info_size
    ! copy the info object from the memory window
    info = TRANSFER(bufr_metainfo((offset+1):(offset+info_size)), dummy_info)
  END SUBROUTINE metainfo_get_from_memwin

END MODULE mo_name_list_output_metadata

