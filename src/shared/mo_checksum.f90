!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Utility functions to ensure that DATA does NOT change due to code modifications.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_checksum
    USE ISO_C_BINDING, ONLY: C_INT64_T, C_INT32_T, C_DOUBLE, C_FLOAT
    USE mo_mpi, ONLY: p_comm_size, p_comm_rank, p_comm_work
    USE mo_exception, ONLY: finish
#ifndef NOMPI
#ifndef __SUNPRO_F95
    USE mpi !, ONLY: MPI_INT32_T, MPI_GATHER
#endif
#endif

    IMPLICIT NONE

    PUBLIC printChecksum

    PRIVATE

#ifndef NOMPI
#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#endif
#endif

    INTERFACE printChecksum
        MODULE PROCEDURE printChecksum_1d_int32
        MODULE PROCEDURE printChecksum_1d_int64
        MODULE PROCEDURE printChecksum_1d_float
        MODULE PROCEDURE printChecksum_1d_double

        MODULE PROCEDURE printChecksum_2d_int32
        MODULE PROCEDURE printChecksum_2d_int64
        MODULE PROCEDURE printChecksum_2d_float
        MODULE PROCEDURE printChecksum_2d_double

        MODULE PROCEDURE printChecksum_3d_int32
        MODULE PROCEDURE printChecksum_3d_int64
        MODULE PROCEDURE printChecksum_3d_float
        MODULE PROCEDURE printChecksum_3d_double

        MODULE PROCEDURE printChecksum_4d_int32
        MODULE PROCEDURE printChecksum_4d_int64
        MODULE PROCEDURE printChecksum_4d_float
        MODULE PROCEDURE printChecksum_4d_double

        MODULE PROCEDURE printChecksum_5d_int32
        MODULE PROCEDURE printChecksum_5d_int64
        MODULE PROCEDURE printChecksum_5d_float
        MODULE PROCEDURE printChecksum_5d_double
    END INTERFACE

    INTEGER(KIND = C_INT32_T) :: mold(1)    ! fortran needs a variable of the TARGET TYPE for a TRANSFER(), so this IS it.
    CHARACTER(LEN = 1), PARAMETER :: kNibbles(16) = (/'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'/)
    CHARACTER(LEN = *), PARAMETER :: moduleName = "mo_checksum"

CONTAINS

    CHARACTER(LEN = 8) FUNCTION checksumString(VALUE) RESULT(RESULT)
        INTEGER(KIND = C_INT64_T), VALUE :: VALUE   !ONLY 32 bits are used, but there IS no unsigned IN fortran

        INTEGER :: i
        CHARACTER(LEN = *), PARAMETER :: routine = moduleName//":checksumString"

        IF(VALUE < 0 .OR. VALUE >= 2_C_INT64_T**32) CALL finish(routine, "VALUE range error")
        DO i = 1, 8
            RESULT(i:i) = kNibbles(IAND(15, VALUE) + 1)
            VALUE = ISHFT(VALUE, -4)
        END DO
    END FUNCTION checksumString

    ! This IS the base CASE, all other "implementations" redirect to this functions via a TRANSFER() CALL.
    !
    ! While this algorithm IS NOT a cryptographical hash, it should be reasonably robust:
    ! It IS guaranteed to catch any single bitflip, AND it IS sensitive to the order of the values,
    ! i. e. printChecksum((/ 0, 1 /)) AND printChecksum((/ 1, 0 /)) produce two different results.
    SUBROUTINE printChecksum_1d_int32(prefix, array, opt_comm)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:)
        INTEGER, OPTIONAL, INTENT(IN) :: opt_comm

        INTEGER :: i, communicator, processCount, error
        INTEGER(KIND = C_INT64_T) :: checksum, pseudoRandomBits
        INTEGER(KIND = C_INT64_T), ALLOCATABLE :: processChecksums(:)
        INTEGER(KIND = C_INT64_T), PARAMETER :: prime1 = 2131876679, prime2 = 1665879031    !just two random primes IN the range [2**30, 2**31]
        INTEGER(KIND = C_INT64_T), PARAMETER :: mask = 2_C_INT64_T**32 - 1_C_INT64_T
        CHARACTER(LEN = *), PARAMETER :: routine = moduleName//":printChecksum_1d_int32"

        !compute a process local checksum
        checksum = 0
        pseudoRandomBits = 0
        DO i = 1, SIZE(array, 1)
            checksum = checksum + IEOR(pseudoRandomBits, array(i))  !every entry IS xor'ed with a different bit pattern
            checksum = IAND(mask, checksum + ISHFT(checksum, -32))   !reduce back to 32 bits
            pseudoRandomBits = IAND(mask, pseudoRandomBits + prime1)
        END DO

#ifndef NOMPI
        !gather the process local checksums on process 0
        communicator = p_comm_work
        IF(PRESENT(opt_comm)) communicator = opt_comm
        processCount = p_comm_size(communicator)
        ALLOCATE(processChecksums(processCount))
        CALL MPI_GATHER(checksum, 1, MPI_INT64_T, processChecksums, 1, MPI_INT64_T, 0, communicator, error)
        IF(error /= MPI_SUCCESS) CALL finish(routine, "error in MPI_Gather()")

        !hash the results of the different processes down to a single VALUE AND print that.
        IF(p_comm_rank(communicator) == 0) THEN
            checksum = 0
            pseudoRandomBits = 0
            DO i = 1, processCount
                checksum = checksum + IEOR(pseudoRandomBits, processChecksums(i))  !every entry IS xor'ed with a different bit pattern
                checksum = IAND(mask, checksum + ISHFT(checksum, -32))   !reduce back to 32 bits
                pseudoRandomBits = IAND(mask, pseudoRandomBits + prime2)
            END DO

            !print the RESULT
            print*, prefix//checksumString(checksum)
        END IF
        DEALLOCATE(processChecksums)
#else
        print*, prefix//checksumString(checksum)
#endif
    END SUBROUTINE printChecksum_1d_int32



    ! Derived cases. These all CALL through to printChecksum_1d_int32().

    SUBROUTINE printChecksum_1d_int64(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_1d_int64

    SUBROUTINE printChecksum_1d_float(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_1d_float

    SUBROUTINE printChecksum_1d_double(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_1d_double



    SUBROUTINE printChecksum_2d_int32(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_2d_int32

    SUBROUTINE printChecksum_2d_int64(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_2d_int64

    SUBROUTINE printChecksum_2d_float(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_2d_float

    SUBROUTINE printChecksum_2d_double(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_2d_double



    SUBROUTINE printChecksum_3d_int32(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_3d_int32

    SUBROUTINE printChecksum_3d_int64(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_3d_int64

    SUBROUTINE printChecksum_3d_float(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_3d_float

    SUBROUTINE printChecksum_3d_double(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_3d_double



    SUBROUTINE printChecksum_4d_int32(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_4d_int32

    SUBROUTINE printChecksum_4d_int64(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_4d_int64

    SUBROUTINE printChecksum_4d_float(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_4d_float

    SUBROUTINE printChecksum_4d_double(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_4d_double



    SUBROUTINE printChecksum_5d_int32(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_5d_int32

    SUBROUTINE printChecksum_5d_int64(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_5d_int64

    SUBROUTINE printChecksum_5d_float(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_5d_float

    SUBROUTINE printChecksum_5d_double(prefix, array)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:,:,:)
        CALL printChecksum(prefix, TRANSFER(array, mold))
    END SUBROUTINE printChecksum_5d_double

END MODULE mo_checksum
