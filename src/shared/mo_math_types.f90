!>
!!   Contains basic math types
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_math_types
    USE ISO_C_BINDING, ONLY: C_INT64_T
    USE mo_fortran_tools, ONLY: t_Destructible
    USE mo_kind, ONLY: wp, sp, dp
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_cartesian_coordinates
    PUBLIC :: t_geographical_coordinates
    PUBLIC :: t_line
    PUBLIC :: t_tangent_vectors
    PUBLIC :: t_Statistics

    ! cartesian coordinate class
    TYPE t_cartesian_coordinates
        REAL(wp) :: x(3)
    END TYPE t_cartesian_coordinates

    ! geographical coordinate class
    TYPE t_geographical_coordinates
        REAL(wp) :: lon
        REAL(wp) :: lat
    END TYPE t_geographical_coordinates

    ! the two coordinates on the tangent plane
    TYPE t_tangent_vectors
        REAL(wp) :: v1
        REAL(wp) :: v2
    END TYPE t_tangent_vectors

    ! line class
    TYPE t_line
        TYPE(t_geographical_coordinates) :: p1
        TYPE(t_geographical_coordinates) :: p2
    END TYPE t_line

    TYPE, EXTENDS(t_Destructible) :: t_Statistics
        INTEGER(C_INT64_T) :: sampleCount
        REAL(wp) :: MIN, mean, MAX
    CONTAINS
        PROCEDURE :: construct => statistics_construct
        GENERIC :: add => addData_s1d, addData_d1d, addStatistics
        PROCEDURE :: destruct => statistics_destruct

        ! scan a given array AND update the statistics accordingly
        PROCEDURE :: addData_s1d => stastistics_addData_s1d
        PROCEDURE :: addData_d1d => stastistics_addData_d1d
        PROCEDURE :: addStatistics => stastistics_addStatistics  ! update the statistics with the contents of another t_Stastics object
    END TYPE t_Statistics

CONTAINS

    SUBROUTINE statistics_construct(me)
        CLASS(t_Statistics), INTENT(INOUT) :: me

        me%sampleCount = 0_C_INT64_T
        me%MIN = HUGE(me%MIN)
        me%mean = 0.0_wp
        me%MAX = -HUGE(me%MAX)
    END SUBROUTINE statistics_construct

    SUBROUTINE stastistics_addData_s1d(me, DATA)
        CLASS(t_Statistics), INTENT(INOUT) :: me
        REAL(sp), INTENT(IN) :: DATA(:)

        TYPE(t_Statistics) :: newStatistics

        CALL newStatistics%construct()
        newStatistics%sampleCount = SIZE(DATA)
        newStatistics%MIN = REAL(MINVAL(DATA), wp)
        newStatistics%mean = REAL(SUM(DATA), wp)/REAL(newStatistics%sampleCount, wp)
        newStatistics%MAX = REAL(MAXVAL(DATA), wp)
        CALL me%add(newStatistics)
        CALL newStatistics%destruct()
    END SUBROUTINE stastistics_addData_s1d

    SUBROUTINE stastistics_addData_d1d(me, DATA)
        CLASS(t_Statistics), INTENT(INOUT) :: me
        REAL(dp), INTENT(IN) :: DATA(:)

        TYPE(t_Statistics) :: newStatistics

        CALL newStatistics%construct()
        newStatistics%sampleCount = SIZE(DATA)
        newStatistics%MIN = REAL(MINVAL(DATA), wp)
        newStatistics%mean = REAL(SUM(DATA), wp)/REAL(newStatistics%sampleCount, wp)
        newStatistics%MAX = REAL(MAXVAL(DATA), wp)
        CALL me%add(newStatistics)
        CALL newStatistics%destruct()
    END SUBROUTINE stastistics_addData_d1d

    SUBROUTINE stastistics_addStatistics(me, other)
        CLASS(t_Statistics), INTENT(INOUT) :: me
        CLASS(t_Statistics), INTENT(IN) :: other

        INTEGER(C_INT64_T) :: newSampleCount

        newSampleCount = me%sampleCount + other%sampleCount
        me%MIN = MIN(me%MIN, other%MIN)
        me%mean = (me%mean*REAL(me%sampleCount, wp) + other%mean*REAL(other%sampleCount, wp))/REAL(newSampleCount, wp)
        me%MAX = MAX(me%MAX, other%MAX)
        me%sampleCount = newSampleCount
    END SUBROUTINE stastistics_addStatistics

    SUBROUTINE statistics_destruct(me)
        CLASS(t_Statistics), INTENT(INOUT) :: me
        ! empty
    END SUBROUTINE statistics_destruct

END MODULE mo_math_types
