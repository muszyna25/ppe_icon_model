!> Interfaces to standart C library functions
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

MODULE mo_std_c_lib
    USE ISO_C_BINDING, ONLY: C_PTR, C_CHAR, C_INT, C_SIZE_T, C_ASSOCIATED, C_F_POINTER
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_util_string, ONLY: int2string

    IMPLICIT NONE

    PUBLIC :: strerror

    CHARACTER(*), PARAMETER :: modname = "mo_std_c_lib"

CONTAINS

    FUNCTION strerror(errorNumber) RESULT(resultVar)
        CHARACTER(:), ALLOCATABLE :: resultVar
        INTEGER, VALUE :: errorNumber

        TYPE(C_PTR) :: c_pointer
        INTEGER :: charPointerShape(1), error, i
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER :: f_pointer
        CHARACTER(*), PARAMETER :: routine = modname//":strerror"

        INTERFACE
            FUNCTION c_strerror(c_error) BIND(C, NAME = "strerror") RESULT(c_result)
                IMPORT C_INT, C_PTR
                TYPE(C_PTR) :: c_result
                INTEGER(C_INT), VALUE :: c_error
            END FUNCTION c_strerror

            INTEGER(C_SIZE_T) FUNCTION c_strlen(charPtr) BIND(C, NAME = "strlen")
                IMPORT C_SIZE_T, C_PTR
                TYPE(C_PTR), VALUE :: charPtr
            end function c_strlen
        END INTERFACE

        c_pointer = c_strerror(errorNumber)
        IF(C_ASSOCIATED(c_pointer)) THEN
            charPointerShape(1) = INT(c_strlen(c_pointer))
            CALL C_F_POINTER(c_pointer, f_pointer, charPointerShape)
            ALLOCATE(CHARACTER(LEN = charPointerShape(1)) :: resultVar, STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
            DO i = 1, charPointerShape(1)
                resultVar(i:i) = f_pointer(i)
            END DO
        ELSE
            CALL finish(routine, "strerror("//TRIM(int2string(errorNumber))//") returned NULL")
        END IF
    END FUNCTION strerror

END MODULE mo_std_c_lib
