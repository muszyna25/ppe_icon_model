!> Interfaces to the functions defined IN support/util_multifile_restart.c
!!
!! Initial implementation: Nathanael Hübbe
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_c_restart_util
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_NULL_CHAR

    IMPLICIT NONE

    PUBLIC :: createEmptyMultifileDir
    PUBLIC :: checkMultifileDir

    PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_c_restart_util"

CONTAINS

    INTEGER FUNCTION createEmptyMultifileDir(path) RESULT(resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: path

        CHARACTER(KIND = C_CHAR) :: pathCopy(LEN(path) + 1)
        INTEGER :: i

        INTERFACE
            INTEGER(C_INT) FUNCTION c_createEmptyMultifileDir(c_path) BIND(C, NAME = "createEmptyMultifileDir")
                IMPORT C_INT, C_CHAR
                CHARACTER(KIND = C_CHAR) :: c_path(*)
            END FUNCTION c_createEmptyMultifileDir
        END INTERFACE

        DO i = 1, LEN(path)
            pathCopy(i) = path(i:i)
        END DO
        pathCopy(LEN(path) + 1) = C_NULL_CHAR
        resultVar = c_createEmptyMultifileDir(pathCopy)
    END FUNCTION createEmptyMultifileDir

    INTEGER FUNCTION checkMultifileDir(path, expectedDomainCount, expectedFileCount) RESULT(resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: path
        INTEGER, VALUE :: expectedDomainCount, expectedFileCount

        CHARACTER(KIND = C_CHAR) :: pathCopy(LEN(path) + 1)
        INTEGER :: i

        INTERFACE
            INTEGER(C_INT) FUNCTION c_checkMultifileDir(c_path, c_expectedDomainCount, c_expectedFileCount) &
              &    BIND(C, NAME = "checkMultifileDir")
                IMPORT C_INT, C_CHAR
                CHARACTER(KIND = C_CHAR) :: c_path(*)
                INTEGER, VALUE :: c_expectedDomainCount, c_expectedFileCount
            END FUNCTION c_checkMultifileDir
        END INTERFACE

        DO i = 1, LEN(path)
            pathCopy(i) = path(i:i)
        END DO
        pathCopy(LEN(path) + 1) = C_NULL_CHAR
        resultVar = c_checkMultifileDir(pathCopy, expectedDomainCount, expectedFileCount)
    END FUNCTION checkMultifileDir

END MODULE mo_c_restart_util
