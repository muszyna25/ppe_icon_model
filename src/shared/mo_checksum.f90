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
    USE mo_util_string, ONLY: int2string
#ifndef NOMPI
#ifndef __SUNPRO_F95
    USE mpi !, ONLY: MPI_INT32_T, MPI_GATHER
#endif
#endif

    IMPLICIT NONE

    PUBLIC printChecksum, printLocalChecksum

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

    INTERFACE printLocalChecksum
        MODULE PROCEDURE printLocalChecksum_1d_int32
        MODULE PROCEDURE printLocalChecksum_1d_int64
        MODULE PROCEDURE printLocalChecksum_1d_float
        MODULE PROCEDURE printLocalChecksum_1d_double

        MODULE PROCEDURE printLocalChecksum_2d_int32
        MODULE PROCEDURE printLocalChecksum_2d_int64
        MODULE PROCEDURE printLocalChecksum_2d_float
        MODULE PROCEDURE printLocalChecksum_2d_double

        MODULE PROCEDURE printLocalChecksum_3d_int32
        MODULE PROCEDURE printLocalChecksum_3d_int64
        MODULE PROCEDURE printLocalChecksum_3d_float
        MODULE PROCEDURE printLocalChecksum_3d_double

        MODULE PROCEDURE printLocalChecksum_4d_int32
        MODULE PROCEDURE printLocalChecksum_4d_int64
        MODULE PROCEDURE printLocalChecksum_4d_float
        MODULE PROCEDURE printLocalChecksum_4d_double

        MODULE PROCEDURE printLocalChecksum_5d_int32
        MODULE PROCEDURE printLocalChecksum_5d_int64
        MODULE PROCEDURE printLocalChecksum_5d_float
        MODULE PROCEDURE printLocalChecksum_5d_double
    END INTERFACE

    INTERFACE checksum
        MODULE PROCEDURE checksum32
        MODULE PROCEDURE checksum64
    END INTERFACE

    INTEGER(KIND = C_INT32_T) :: mold(1)    ! fortran needs a variable of the TARGET TYPE for a TRANSFER(), so this IS it.
    CHARACTER(LEN = 1), PARAMETER :: kNibbles(16) = (/'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'/)
    INTEGER(KIND = C_INT64_T), PARAMETER :: mask = 2_C_INT64_T**32 - 1_C_INT64_T    ! a bitmask for the 32 low order bits
    CHARACTER(LEN = *), PARAMETER :: moduleName = "mo_checksum"

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! real functionality !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CHARACTER(LEN = 8) FUNCTION checksumString(VALUE) RESULT(RESULT)
        INTEGER(KIND = C_INT64_T), VALUE :: VALUE   !ONLY 32 bits are used, but there IS no unsigned IN fortran

        INTEGER :: i
        CHARACTER(LEN = *), PARAMETER :: routine = moduleName//":checksumString"

        IF(VALUE < 0 .OR. VALUE >= 2_C_INT64_T**32) CALL finish(routine, "VALUE range error")
        DO i = 1, 8
            RESULT(9-i:9-i) = kNibbles(IAND(15_C_INT64_T, VALUE) + 1_C_INT64_T)
            VALUE = ISHFT(VALUE, -4)
        END DO
    END FUNCTION checksumString

    ! This IS the base checksum FUNCTION which IS used to implement both printChecksum_1d_int32() AND printLocalChecksum_1d_int32().
    !
    ! While this algorithm IS NOT a cryptographical hash, it should be reasonably robust:
    ! It IS guaranteed to catch any single bitflip, AND it IS sensitive to the order of the values,
    ! i. e. printChecksum((/ 0, 1 /)) AND printChecksum((/ 1, 0 /)) produce two different results.
    INTEGER(KIND = C_INT64_T) FUNCTION checksum32(array, prime) RESULT(RESULT)
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:)
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: prime

        INTEGER :: i
        INTEGER(KIND = C_INT64_T) :: pseudoRandomBits
        CHARACTER(LEN = *), PARAMETER :: routine = moduleName//":checksum32"

        RESULT = 0
        pseudoRandomBits = 0
        DO i = 1, SIZE(array, 1)
            RESULT = RESULT + IEOR(pseudoRandomBits, INT(array(i), C_INT64_T))  !every entry IS xor'ed with a different bit pattern
            RESULT = IAND(mask, RESULT + ISHFT(RESULT, -32))   !reduce back to 32 bits
            pseudoRandomBits = IAND(mask, pseudoRandomBits + prime)
        END DO
    END FUNCTION checksum32

    ! This IS the base CASE for the global checksums, all other "implementations" redirect to this FUNCTION via a TRANSFER() CALL.
    !
    ! If `opt_lDetails = .TRUE.` IS specified, this also prints a list of all the process local checksums.
    SUBROUTINE printChecksum_1d_int32(prefix, array, opt_comm, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:)
        INTEGER, OPTIONAL, INTENT(IN) :: opt_comm
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails

        INTEGER :: i, communicator, processCount, error
        LOGICAL :: printDetails
        INTEGER(KIND = C_INT64_T) :: hash
        INTEGER(KIND = C_INT64_T), ALLOCATABLE :: processChecksums(:)
        INTEGER(KIND = C_INT64_T), PARAMETER :: prime1 = 2131876679, prime2 = 1665879031    !just two random primes IN the range
                                                                                            ![2**30, 2**31]
        CHARACTER(LEN = *), PARAMETER :: routine = moduleName//":printChecksum_1d_int32"

        !XXX: These two variables are a workaround for the MPI implementation on AIX, which does NOT provide the constants
        !MPI_INT64_T AND MPI_INT32_T.  So, to make this work without those constants, we implicitly reinterprete a C_INT64_T as an
        !array of fortran INTEGERs during the MPI_Gather() CALL.  Warning: This assumes that the SIZE of a fortran INTEGER is a
        !divisor of eight. Should be TRUE on any sane system, but you never know.
        INTEGER :: integerMold(1), integersInInt64
        integersInInt64 = SIZE(TRANSFER(hash, integerMold))

        !compute a process local checksum
        printDetails = .FALSE.
        IF(PRESENT(opt_lDetails)) printDetails = opt_lDetails
        hash = checksum(array, prime1)

#ifndef NOMPI
        !gather the process local checksums on process 0
        communicator = p_comm_work
        IF(PRESENT(opt_comm)) communicator = opt_comm
        processCount = p_comm_size(communicator)
        ALLOCATE(processChecksums(processCount))
        !XXX: Dirty hack ahead. Reinterpreting INTEGER(KIND = C_INT64_T) as array of INTEGER. See comment on integersInInt64.
        CALL MPI_GATHER(hash, integersInInt64, MPI_INTEGER, processChecksums, integersInInt64, MPI_INTEGER, 0, communicator, error)
        IF(error /= MPI_SUCCESS) CALL finish(routine, "error in MPI_Gather()")

        !hash the results of the different processes down to a single VALUE AND print that.
        IF(p_comm_rank(communicator) == 0) THEN
            hash = checksum(processChecksums, prime2)
            IF(printDetails) THEN
                print*, prefix//"details:"
                DO i = 1, processCount
                    print*, "checksum from process "//TRIM(int2string(i - 1))//": "//checksumString(processChecksums(i))
                END DO
            END IF

            !print the RESULT
            print*, prefix//checksumString(hash)
        END IF
        DEALLOCATE(processChecksums)
#else
        print*, prefix//checksumString(hash)
#endif
    END SUBROUTINE printChecksum_1d_int32


    ! This IS the base CASE for the local checksums, all other implementations redirect to this FUNCTION via a TRANSFER() CALL.
    !
    ! IF `opt_lDetails = .TRUE.` IS specified, this also produces a hex dump of the input DATA.
    SUBROUTINE printLocalChecksum_1d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails

        INTEGER :: i
        LOGICAL :: printDetails
        INTEGER(KIND = C_INT64_T) :: hash
        INTEGER(KIND = C_INT64_T), PARAMETER :: prime = 1212113153    !just a random prime IN the range [2**30, 2**31]

        printDetails = .FALSE.
        IF(PRESENT(opt_lDetails)) printDetails = opt_lDetails
        IF(printDetails) THEN
            print*, prefix//"hex dump:"
            DO i = 0, SIZE(array, 1) - 4, 4
                print*, TRIM(int2string(i*4)) &
                    & //": "//checksumString(IAND(mask, INT(array(i + 1), C_INT64_T))) &
                    & // " "//checksumString(IAND(mask, INT(array(i + 2), C_INT64_T))) &
                    & // " "//checksumString(IAND(mask, INT(array(i + 3), C_INT64_T))) &
                    & // " "//checksumString(IAND(mask, INT(array(i + 4), C_INT64_T)))
            END DO
            SELECT CASE(MOD(SIZE(array, 1),4))
                CASE(1)
                    print*, TRIM(int2string(SIZE(array, 1) - 1)) &
                        & //": "//checksumString(IAND(mask, INT(array(SIZE(array, 1) - 0), C_INT64_T)))
                CASE(2)
                    print*, TRIM(int2string(SIZE(array, 1) - 2)) &
                        & //": "//checksumString(IAND(mask, INT(array(SIZE(array, 1) - 1), C_INT64_T))) &
                        & // " "//checksumString(IAND(mask, INT(array(SIZE(array, 1) - 0), C_INT64_T)))
                CASE(3)
                    print*, TRIM(int2string(SIZE(array, 1) - 3)) &
                        & //": "//checksumString(IAND(mask, INT(array(SIZE(array, 1) - 2), C_INT64_T))) &
                        & //": "//checksumString(IAND(mask, INT(array(SIZE(array, 1) - 1), C_INT64_T))) &
                        & // " "//checksumString(IAND(mask, INT(array(SIZE(array, 1) - 0), C_INT64_T)))
            END SELECT
        END IF

        hash = checksum(array, prime)
        print*, prefix//checksumString(hash)
    END SUBROUTINE printLocalChecksum_1d_int32

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! callthroughs to the functions above !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Callthrough to checksum32().

    INTEGER(KIND = C_INT64_T) FUNCTION checksum64(array, prime) RESULT(RESULT)
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:)
        INTEGER(KIND = C_INT64_T), VALUE :: prime

        RESULT = checksum(TRANSFER(array, mold), prime)
    END FUNCTION checksum64



    ! These CALL through to printChecksum_1d_int32().

    SUBROUTINE printChecksum_1d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_1d_int64

    SUBROUTINE printChecksum_1d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_1d_float

    SUBROUTINE printChecksum_1d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_1d_double



    SUBROUTINE printChecksum_2d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_int32

    SUBROUTINE printChecksum_2d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_int64

    SUBROUTINE printChecksum_2d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_float

    SUBROUTINE printChecksum_2d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_double



    SUBROUTINE printChecksum_3d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_int32

    SUBROUTINE printChecksum_3d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_int64

    SUBROUTINE printChecksum_3d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_float

    SUBROUTINE printChecksum_3d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(RESHAPE(array,(/size(array)/)), mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_double



    SUBROUTINE printChecksum_4d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_int32

    SUBROUTINE printChecksum_4d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_int64

    SUBROUTINE printChecksum_4d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_float

    SUBROUTINE printChecksum_4d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_double



    SUBROUTINE printChecksum_5d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_int32

    SUBROUTINE printChecksum_5d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_int64

    SUBROUTINE printChecksum_5d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_float

    SUBROUTINE printChecksum_5d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_double



    ! These CALL through to printLocalChecksum_1d_int32().

    SUBROUTINE printLocalChecksum_1d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_1d_int64

    SUBROUTINE printLocalChecksum_1d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_1d_float

    SUBROUTINE printLocalChecksum_1d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_1d_double



    SUBROUTINE printLocalChecksum_2d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_2d_int32

    SUBROUTINE printLocalChecksum_2d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_2d_int64

    SUBROUTINE printLocalChecksum_2d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_2d_float

    SUBROUTINE printLocalChecksum_2d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_2d_double



    SUBROUTINE printLocalChecksum_3d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_3d_int32

    SUBROUTINE printLocalChecksum_3d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_3d_int64

    SUBROUTINE printLocalChecksum_3d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_3d_float

    SUBROUTINE printLocalChecksum_3d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_3d_double



    SUBROUTINE printLocalChecksum_4d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_4d_int32

    SUBROUTINE printLocalChecksum_4d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_4d_int64

    SUBROUTINE printLocalChecksum_4d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_4d_float

    SUBROUTINE printLocalChecksum_4d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_4d_double



    SUBROUTINE printLocalChecksum_5d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_5d_int32

    SUBROUTINE printLocalChecksum_5d_int64(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        INTEGER(KIND = C_INT64_T), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_5d_int64

    SUBROUTINE printLocalChecksum_5d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_FLOAT), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_5d_float

    SUBROUTINE printLocalChecksum_5d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *), INTENT(IN) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(IN) :: array(:,:,:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDetails
        CALL printLocalChecksum(prefix, TRANSFER(array, mold), opt_lDetails = opt_lDetails)
    END SUBROUTINE printLocalChecksum_5d_double

END MODULE mo_checksum
