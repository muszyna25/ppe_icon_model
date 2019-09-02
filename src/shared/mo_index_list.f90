MODULE mo_index_list

  USE mo_kind, ONLY: i1, i2, i4

#ifdef _OPENACC
  USE openacc
  USE, INTRINSIC :: iso_c_binding
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generate_index_list

  INTERFACE generate_index_list
    MODULE PROCEDURE generate_index_list_i1
    MODULE PROCEDURE generate_index_list_batched_i1
  END INTERFACE

#ifdef _OPENACC

  INTERFACE generate_index_list_cuda

    SUBROUTINE generate_index_list_cuda_i1(conditions, startid, endid, indices, nvalid, stream) &
        & BIND(C, name="c_generate_index_list_cuda_i1")

      USE openacc
      USE iso_c_binding

      TYPE(c_devptr),           INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_devptr),           INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      INTEGER(c_int),           INTENT(out)         :: nvalid
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_cuda_i1


    ! NOT IMPLEMENTED (YET)
    SUBROUTINE generate_index_list_cuda_batched_i1                        &
        & (conditions, nbatches, startid, endid, indices, nvalid, stream) &
        & BIND(C, name="c_generate_index_list_cuda_batched_i1")

      USE openacc
      USE iso_c_binding

      TYPE(c_devptr),           INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: nbatches
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_devptr),           INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      INTEGER(c_int),           INTENT(out)         :: nvalid
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_cuda_batched_i1

  END INTERFACE

#endif

  CONTAINS

#ifndef _OPENACC

! Regular CPU implementation with a simple loop

  SUBROUTINE generate_index_list_i1(conditions, indices, startid, endid, nvalid)
    INTEGER(i1), INTENT(in)    :: conditions(:)
    INTEGER,     INTENT(inout) :: indices(:)
    INTEGER,     INTENT(in)    :: startid
    INTEGER,     INTENT(in)    :: endid
    INTEGER,     INTENT(out)   :: nvalid

    INTEGER :: i

    nvalid = 0
    DO i = startid, endid
      IF (conditions(i) /= 0) THEN
        nvalid = nvalid + 1
        indices(nvalid) = i
      END IF
    END DO
  END SUBROUTINE generate_index_list_i1

  SUBROUTINE generate_index_list_batched_i1(conditions, indices, startid, endid, nvalid)
    INTEGER(i1), INTENT(in)    :: conditions(:,:)
    INTEGER,     INTENT(inout) :: indices(:,:)
    INTEGER,     INTENT(in)    :: startid
    INTEGER,     INTENT(in)    :: endid
    INTEGER,     INTENT(inout) :: nvalid(:)

    INTEGER :: i, batch, nbatches(2)
    nbatches = SHAPE(conditions)
    nvalid = 0

    DO batch = 1, nbatches(2)
      CALL generate_index_list_i1(              &
        conditions(:,batch), indices(:, batch), &
        startid, endid, nvalid(batch) )
    END DO

  END SUBROUTINE generate_index_list_i1

#else

! On the GPU call the CUB library through C++

  SUBROUTINE generate_index_list_i1(conditions, indices, startid, endid, nvalid, acc_async_queue)
    INTEGER(i1), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    INTEGER,     INTENT(in), OPTIONAL :: acc_async_queue

    INTEGER(acc_handle_kind) :: stream

    IF ( PRESENT(acc_async_queue) ) THEN
      stream = acc_get_cuda_stream(acc_async_queue)
    ELSE
      stream = acc_get_cuda_stream(acc_async_sync)
    END IF

    CALL generate_index_list_cuda_i1(           &
      & acc_deviceptr(conditions),              &
      & startid, endid,                         &
      & acc_deviceptr(indices), nvalid, stream)

  END SUBROUTINE generate_index_list_i1


  ! This has to be implemented as a proper batch thingy
  SUBROUTINE generate_index_list_batched_i1(conditions, indices, startid, endid, nvalid, acc_async_queue)
    INTEGER(i1), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid(:)
    INTEGER,     INTENT(in), OPTIONAL :: acc_async_queue

    INTEGER(acc_handle_kind) :: stream
    INTEGER :: i, batch, nbatches(2)


    IF ( PRESENT(acc_async_queue) ) THEN
      stream = acc_get_cuda_stream(acc_async_queue)
    ELSE
      stream = acc_get_cuda_stream(acc_async_sync)
    END IF

    nbatches = SHAPE(conditions)

    DO batch = 1, nbatches(2)
      CALL generate_index_list_cuda_i1(           &
        & acc_deviceptr(conditions(:,batch)),     &
        & startid, endid,                         &
        & acc_deviceptr(indices(:,batch)),        &
        & nvalid(batch), stream)
    END DO

  END SUBROUTINE generate_index_list_batched_i1

#endif



END MODULE mo_index_list
