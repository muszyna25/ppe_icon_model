MODULE mo_index_list

  USE mo_kind, ONLY: i1, i4

#ifdef _OPENACC
  USE openacc
  USE, INTRINSIC :: iso_c_binding
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generate_index_list, generate_index_list_batched

  INTERFACE generate_index_list
#ifdef _OPENACC
    MODULE PROCEDURE generate_index_list_i1_gpu
    MODULE PROCEDURE generate_index_list_i4_gpu
#else
    MODULE PROCEDURE generate_index_list_i1
    MODULE PROCEDURE generate_index_list_i4
#endif
  END INTERFACE

! Warning: there will be no GPU -> CPU copy here, array of NUMBER
!  of indices will ONLY be on the GPU!!
  INTERFACE generate_index_list_batched
    MODULE PROCEDURE generate_index_list_batched_i1
    MODULE PROCEDURE generate_index_list_batched_i4
  END INTERFACE

#ifdef _OPENACC

  INTERFACE generate_index_list_cuda

    ! Int8 or char
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

    ! Int32
    SUBROUTINE generate_index_list_cuda_i4(conditions, startid, endid, indices, nvalid, stream) &
        & BIND(C, name="c_generate_index_list_cuda_i4")

      USE openacc
      USE iso_c_binding

      TYPE(c_devptr),           INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_devptr),           INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      INTEGER(c_int),           INTENT(out)         :: nvalid
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_cuda_i4


    SUBROUTINE generate_index_list_cuda_batched_i1               &
        & (batch_size,                                           &
        &  conditions, cond_stride,                              &
        &  startid, endid,                                       &
        &  indices, idx_stride,                                  &
        &  nvalid, stream )                                      &
        & BIND(C, name="c_generate_index_list_cuda_batched_i1")

      USE openacc
      USE iso_c_binding

      INTEGER(c_int),           INTENT(in),  VALUE  :: batch_size
      TYPE(c_devptr),           INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: cond_stride
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_devptr),           INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      INTEGER(c_int),           INTENT(in),  VALUE  :: idx_stride
      TYPE(c_devptr),           INTENT(in),  VALUE  :: nvalid
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_cuda_batched_i1


    SUBROUTINE generate_index_list_cuda_batched_i4               &
        & (batch_size,                                           &
        &  conditions, cond_stride,                              &
        &  startid, endid,                                       &
        &  indices, idx_stride,                                  &
        &  nvalid, stream )                                      &
        & BIND(C, name="c_generate_index_list_cuda_batched_i4")

      USE openacc
      USE iso_c_binding

      INTEGER(c_int),           INTENT(in),  VALUE  :: batch_size
      TYPE(c_devptr),           INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: cond_stride
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_devptr),           INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      INTEGER(c_int),           INTENT(in),  VALUE  :: idx_stride
      TYPE(c_devptr),           INTENT(in),  VALUE  :: nvalid
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_cuda_batched_i4

  END INTERFACE

#endif

  CONTAINS


! Regular CPU implementation with a simple loop

  SUBROUTINE generate_index_list_i1(conditions, indices, startid, endid, nvalid, dummy,lacc)
    INTEGER(i1), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    ! These arguments are used in the OpenACC variant, but not in the CPU one 
    INTEGER,     INTENT(in), OPTIONAL :: dummy
    LOGICAL,     INTENT(IN), OPTIONAL :: lacc

    INTEGER :: i

    nvalid = 0
    DO i = startid, endid
      IF (conditions(i) /= 0) THEN
        nvalid = nvalid + 1
        indices(nvalid) = i
      END IF
    END DO
  END SUBROUTINE generate_index_list_i1

  SUBROUTINE generate_index_list_i4(conditions, indices, startid, endid, nvalid, dummy,lacc)
    INTEGER(i4), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    ! These arguments are used in the OpenACC variant, but not in the CPU one 
    INTEGER,     INTENT(in), OPTIONAL :: dummy
    LOGICAL,     INTENT(IN), OPTIONAL :: lacc

    INTEGER :: i

    nvalid = 0
    DO i = startid, endid
      IF (conditions(i) /= 0) THEN
        nvalid = nvalid + 1
        indices(nvalid) = i
      END IF
    END DO
  END SUBROUTINE generate_index_list_i4

#ifndef _OPENACC

  SUBROUTINE generate_index_list_batched_i1(conditions, indices, startid, endid, nvalid, dummy)
    INTEGER(i1), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(inout)        :: nvalid(:)
    ! This argument is used in the OpenACC variant, but not in the GPU one 
    INTEGER,     INTENT(in), OPTIONAL :: dummy

    INTEGER :: i, batch, batch_size, sh(2)
    sh = SHAPE(conditions)
    batch_size = sh(2)
    nvalid = 0

    DO batch = 1, batch_size
      CALL generate_index_list_i1(              &
        conditions(:,batch), indices(:, batch), &
        startid, endid, nvalid(batch) )
    END DO

  END SUBROUTINE generate_index_list_batched_i1

  SUBROUTINE generate_index_list_batched_i4(conditions, indices, startid, endid, nvalid, dummy)
    INTEGER(i4), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(inout)        :: nvalid(:)
    ! This argument is used in the OpenACC variant, but not in the GPU one
    INTEGER,     INTENT(in), OPTIONAL :: dummy

    INTEGER :: i, batch, batch_size, sh(2)
    sh = SHAPE(conditions)
    batch_size = sh(2)
    nvalid = 0

    DO batch = 1, batch_size
      CALL generate_index_list_i4(              &
        conditions(:,batch), indices(:, batch), &
        startid, endid, nvalid(batch) )
    END DO

  END SUBROUTINE generate_index_list_batched_i4

#else

  ! on the gpu call the cub library through c++

  ! 1 byte
  SUBROUTINE generate_index_list_i1_gpu(conditions, indices, startid, endid, nvalid, acc_async_queue,lacc)
    INTEGER(i1), INTENT(IN)           :: conditions(:)
    INTEGER,     INTENT(INOUT)        :: indices(:)
    INTEGER,     INTENT(IN)           :: startid
    INTEGER,     INTENT(IN)           :: endid
    INTEGER,     INTENT(OUT)          :: nvalid
    INTEGER,     INTENT(IN), OPTIONAL :: acc_async_queue
    LOGICAL,     INTENT(IN), OPTIONAL :: lacc

    INTEGER(acc_handle_kind) :: stream
    LOGICAL :: gen_list_on_gpu

    IF (PRESENT(lacc)) THEN
      gen_list_on_gpu = lacc
    ELSE
      gen_list_on_gpu = .TRUE.
    ENDIF

    ! run on GPU
    IF (gen_list_on_gpu) THEN

      IF ( PRESENT(acc_async_queue) ) THEN
        stream = acc_get_cuda_stream(acc_async_queue)
      ELSE
        stream = acc_get_cuda_stream(acc_async_sync)
      END IF

      CALL generate_index_list_cuda_i1(           &
        & acc_deviceptr(conditions),              &
        & startid, endid,                         &
        & acc_deviceptr(indices), nvalid, stream)

    ! run on CPU
    ELSE
      CALL generate_index_list_i1(conditions, indices, startid, endid, nvalid)
    ENDIF


  END SUBROUTINE generate_index_list_i1_gpu

  ! 4 bytes
  SUBROUTINE generate_index_list_i4_gpu(conditions, indices, startid, endid, nvalid, acc_async_queue, lacc)
    INTEGER(i4), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    INTEGER,     INTENT(in), OPTIONAL :: acc_async_queue
    LOGICAL,     INTENT(in), OPTIONAL :: lacc

    INTEGER(acc_handle_kind) :: stream
    LOGICAL :: gen_list_on_gpu

    IF (PRESENT(lacc)) THEN
      gen_list_on_gpu = lacc
    ELSE
      gen_list_on_gpu = .TRUE.
    ENDIF

    ! run on GPU
    IF (gen_list_on_gpu) THEN

      IF ( PRESENT(acc_async_queue) ) THEN
        stream = acc_get_cuda_stream(acc_async_queue)
      ELSE
        stream = acc_get_cuda_stream(acc_async_sync)
      END IF

      CALL generate_index_list_cuda_i4(           &
        & acc_deviceptr(conditions),              &
        & startid, endid,                         &
        & acc_deviceptr(indices), nvalid, stream)

    ! run on CPU
    ELSE
      CALL generate_index_list_i4(conditions, indices, startid, endid, nvalid)
    ENDIF
  END SUBROUTINE generate_index_list_i4_gpu

  !
  ! Batched
  !

  ! 1 byte
  SUBROUTINE generate_index_list_batched_i1(conditions, indices, startid, endid, nvalid, acc_async_queue)
    ! Attention: all the values, including nvalid (!!)
    !  are updated asynchronously on device ONLY
    INTEGER(i1), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid(:)
    INTEGER,     INTENT(in), OPTIONAL :: acc_async_queue

    INTEGER(acc_handle_kind) :: stream
    INTEGER :: i
    INTEGER :: batch_size, sh(2)
    INTEGER :: cond_stride, idx_stride


    IF ( PRESENT(acc_async_queue) ) THEN
      stream = acc_get_cuda_stream(acc_async_queue)
    ELSE
      stream = acc_get_cuda_stream(acc_async_sync)
    END IF

    sh = SHAPE(conditions)
    batch_size = sh(2)

    ! Hacky way to support non-contiguous slices
    IF ( batch_size > 1 ) THEN
      cond_stride = ( LOC(conditions(1,2)) - LOC(conditions(1,1)) ) / SIZEOF(conditions(1,1))
      idx_stride  = ( LOC(indices   (1,2)) - LOC(indices   (1,1)) ) / SIZEOF(indices   (1,1))
    END IF


    CALL generate_index_list_cuda_batched_i1(      &
        & batch_size,                              &
        & acc_deviceptr(conditions), cond_stride,  &
        & startid, endid,                          &
        & acc_deviceptr(indices), idx_stride,      &
        & acc_deviceptr(nvalid), stream )

  END SUBROUTINE generate_index_list_batched_i1

  ! 4 bytes
  SUBROUTINE generate_index_list_batched_i4(conditions, indices, startid, endid, nvalid, acc_async_queue)
    ! Attention: all the values, including nvalid (!!)
    !  are updated asynchronously on device ONLY
    INTEGER(i4), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid(:)
    INTEGER,     INTENT(in), OPTIONAL :: acc_async_queue

    INTEGER(acc_handle_kind) :: stream
    INTEGER :: i
    INTEGER :: batch_size, sh(2)
    INTEGER :: cond_stride, idx_stride


    IF ( PRESENT(acc_async_queue) ) THEN
      stream = acc_get_cuda_stream(acc_async_queue)
    ELSE
      stream = acc_get_cuda_stream(acc_async_sync)
    END IF

    sh = SHAPE(conditions)
    batch_size = sh(2)

    ! Hacky way to support non-contiguous slices
    IF ( batch_size > 1 ) THEN
      cond_stride = ( LOC(conditions(1,2)) - LOC(conditions(1,1)) ) / SIZEOF(conditions(1,1))
      idx_stride  = ( LOC(indices   (1,2)) - LOC(indices   (1,1)) ) / SIZEOF(indices   (1,1))
    END IF


    CALL generate_index_list_cuda_batched_i4(      &
        & batch_size,                              &
        & acc_deviceptr(conditions), cond_stride,  &
        & startid, endid,                          &
        & acc_deviceptr(indices), idx_stride,      &
        & acc_deviceptr(nvalid), stream )

  END SUBROUTINE generate_index_list_batched_i4

#endif



END MODULE mo_index_list
