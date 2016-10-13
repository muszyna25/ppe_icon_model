!>
!! A small class to encapsulate the information of how the grid is decomposed,
!! in order to allow quick and easy distribution of input data.
!!
!! All calls within this class are collective.
!!
!! @author N. Hübbe, DWD
!!
!!
!! @par Revision History
!! Initial hack: 2014-08-05 : N. Hübbe, DWD
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_scatter_pattern_base
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: wp, dp, sp, i8
    USE mo_mpi, ONLY: my_process_is_stdio, p_mpi_wtime, p_max
    USE mo_run_config, ONLY: msg_level

    IMPLICIT NONE

PUBLIC :: t_ScatterPattern, t_ScatterPatternPtr, constructScatterPattern, destructScatterPattern, deleteScatterPattern, &
        & lookupScatterPattern

    TYPE, ABSTRACT :: t_ScatterPattern
        INTEGER :: totalPointCount !< number of points needed as input
        INTEGER :: myPointCount !< number of points needed by this PE
        INTEGER :: jg   !< the domain for which this pattern is used
        INTEGER :: communicator !< the communicator to use
        INTEGER(i8) :: distributedData  !< statistic on how much data was distributed
        REAL(dp) :: curStartTime    !< the time when the last distribution call was started
        REAL(dp) :: distributionTime    !< statistic on how long we took to distribute the data

    CONTAINS
        PROCEDURE(interface_distribute_dp), DEFERRED :: distribute_dp   !< distribute double precision data
        PROCEDURE(interface_distribute_sp), DEFERRED :: distribute_sp   !< distribute single precision data
        PROCEDURE(interface_distribute_int), DEFERRED :: distribute_int   !< distribute single precision data


        PROCEDURE :: construct => constructScatterPattern   !< constructor

        PROCEDURE :: globalSize => ScatterPattern_globalSize    !< the count of expected input points for a distribution
        PROCEDURE :: localSize => ScatterPattern_localSize  !< the count of points received by this PE

        PROCEDURE :: resetStatistics => scatterPatternResetStatistics   !< reset the internal statistics
        PROCEDURE :: printStatistics => scatterPatternPrintStatistics   !< print the current statistics and reset them

        PROCEDURE :: startDistribution => scatterPatternStartDistribution   !For use by subclasses only.
        PROCEDURE :: endDistribution => scatterPatternEndDistribution   !For use by subclasses only.

        PROCEDURE :: destruct => destructScatterPattern !< destructor

        GENERIC :: distribute => distribute_dp, distribute_sp, distribute_int
    END TYPE t_ScatterPattern

    TYPE :: t_ScatterPatternPtr
        CLASS(t_ScatterPattern), POINTER :: ptr
    END TYPE t_ScatterPatternPtr

PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_grid_distribution_base"
    LOGICAL, PARAMETER :: debugModule = .false.

    ABSTRACT INTERFACE
        !---------------------------------------------------------------------------------------------------------------------------
        !> do the data distribution for a double precision array
        !---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE interface_distribute_dp(me, globalArray, localArray, ladd_value)
            IMPORT t_ScatterPattern, dp, wp
            CLASS(t_ScatterPattern), INTENT(INOUT) :: me
            REAL(dp), INTENT(INOUT) :: globalArray(:)
            REAL(wp), INTENT(INOUT) :: localArray(:,:)
            LOGICAL, INTENT(IN) :: ladd_value
        END SUBROUTINE interface_distribute_dp
        !---------------------------------------------------------------------------------------------------------------------------
        !> do the data distribution for a double precision array
        !---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE interface_distribute_sp(me, globalArray, localArray, ladd_value)
            IMPORT t_ScatterPattern, sp, wp
            CLASS(t_ScatterPattern), INTENT(INOUT) :: me
            REAL(sp), INTENT(INOUT) :: globalArray(:)
            REAL(wp), INTENT(INOUT) :: localArray(:,:)
            LOGICAL, INTENT(IN) :: ladd_value
        END SUBROUTINE interface_distribute_sp
        !---------------------------------------------------------------------------------------------------------------------------
        !> do the data distribution for a integer
        !---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE interface_distribute_int(me, globalArray, localArray, ladd_value)
            IMPORT t_ScatterPattern
            CLASS(t_ScatterPattern), INTENT(INOUT) :: me
            INTEGER, INTENT(INOUT) :: globalArray(:)
            INTEGER, INTENT(INOUT) :: localArray(:,:)
            LOGICAL, INTENT(IN) :: ladd_value
        END SUBROUTINE interface_distribute_int
    END INTERFACE

    INTERFACE lookupScatterPattern
        MODULE PROCEDURE ScatterPattern_lookupSize
    END INTERFACE lookupScatterPattern

    TYPE(t_ScatterPatternPtr), POINTER :: existingPatterns(:)
    INTEGER :: existingPatternCount = -1

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Query the list of existing scatter patterns for one with the given domain and global size.
    !> This uses the fact that there are generally three different scatter patterns per domain for the edges, cells, and vertices,
    !> all of which have a different global size. Thus the pair (jg, globalSize) uniquely identifies a scatter pattern.
    !-------------------------------------------------------------------------------------------------------------------------------
    FUNCTION ScatterPattern_lookupSize(jg, globalSize) RESULT(result)
        INTEGER, VALUE :: jg, globalSize
        CLASS(t_ScatterPattern), POINTER :: result

        CHARACTER(*), PARAMETER :: routine = modname//":ScatterPattern_lookupSize"
        INTEGER :: i
        CLASS(t_ScatterPatternPtr), POINTER :: temp

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine

        DO i = 1, existingPatternCount
            temp => existingPatterns(i)
            result => temp%ptr
            IF(result%globalSize() == globalSize .and. result%jg == jg) RETURN
        END DO
        result => NULL()

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END FUNCTION ScatterPattern_lookupSize

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destruct and deallocate a t_ScatterPattern
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE deleteScatterPattern(me)
        CLASS(t_ScatterPattern), POINTER, INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":deleteScatterPattern"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        CALL me%destruct()
        DEALLOCATE(me)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE deleteScatterPattern

    !-------------------------------------------------------------------------------------------------------------------------------
    !> constructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE constructScatterPattern(me, jg, loc_arr_len, glb_index, communicator)
        CLASS(t_ScatterPattern), TARGET, INTENT(OUT) :: me
        INTEGER, VALUE :: jg, loc_arr_len, communicator
        INTEGER, INTENT(IN) :: glb_index(:)

        CHARACTER(*), PARAMETER :: routine = modname//":constructScatterPattern"
        CLASS(t_ScatterPatternPtr), POINTER :: temp(:)
        INTEGER :: i, error

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine

        !init the object itself
        me%totalPointCount = p_max(MAXVAL(glb_index), comm = communicator)
        me%myPointCount = loc_arr_len
        me%jg = jg
        me%communicator = communicator
        me%distributedData = 0
        me%distributionTime = 0.0

        !add the new scatter pattern to the list of existing ones
        IF(existingPatternCount == -1) THEN
            ALLOCATE(existingPatterns(8), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            existingPatternCount = 0
        END IF
        IF(SIZE(existingPatterns, 1) == existingPatternCount) THEN
            temp => existingPatterns
            ALLOCATE(existingPatterns(2*existingPatternCount), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            DO i = 1, existingPatternCount
                existingPatterns(i)%ptr => temp(i)%ptr
            END DO
            DEALLOCATE(temp)
        END IF
        existingPatternCount = existingPatternCount + 1
        existingPatterns(existingPatternCount)%ptr => me

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE constructScatterPattern

    !-------------------------------------------------------------------------------------------------------------------------------
    !> accessor
    !-------------------------------------------------------------------------------------------------------------------------------
    INTEGER FUNCTION ScatterPattern_globalSize(me) RESULT(RESULT)
        CLASS(t_ScatterPattern), INTENT(INOUT) :: me
        RESULT = me%totalPointCount
    END FUNCTION ScatterPattern_globalSize

    !-------------------------------------------------------------------------------------------------------------------------------
    !> accessor
    !-------------------------------------------------------------------------------------------------------------------------------
    INTEGER FUNCTION ScatterPattern_localSize(me) RESULT(RESULT)
        CLASS(t_ScatterPattern), INTENT(INOUT) :: me
        RESULT = me%myPointCount
    END FUNCTION ScatterPattern_localSize

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Reset the internal statistics of how long it took to distribute how much data.
    !> This is used in conjunction with printStatistics() to produce statistics output for only a subset of distribute() calls.
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE scatterPatternResetStatistics(me)
        CLASS(t_ScatterPattern), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":scatterPatternResetStatistics"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        me%distributedData = 0_i8
        me%distributionTime = 0.0_dp
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE scatterPatternResetStatistics

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Make some output (if msg_level is high enough) detailing how long it took to distribute data since the last call to
    !> resetStatistics().
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE scatterPatternPrintStatistics(me)
        CLASS(t_ScatterPattern), INTENT(INOUT) :: me
        REAL(dp) :: bandwidth

        CHARACTER(*), PARAMETER :: routine = modname//":scatterPatternPrintStatistics"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        IF(me%distributedData > 0) THEN
            IF(msg_level >= 10 .and. my_process_is_stdio()) THEN
                WRITE(0,*) routine, ": data distribution totals:"
                WRITE(0,'(8X,A,I19,A)')   "amount:    ", me%distributedData, " bytes"
                WRITE(0,'(8X,A,F19.3,A)') "duration:  ", me%distributionTime, " seconds"
                IF (me%distributionTime == 0.0_dp) THEN
                  bandwidth = -1.0_dp
                ELSE
                  bandwidth = REAL(me%distributedData, dp)/(1048576.0_dp*me%distributionTime)
                END IF
                WRITE(0,'(8X,A,F19.3,A)') "bandwidth: ", bandwidth, " MiB/s"
            END IF
        END IF
        call me%resetStatistics()
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE scatterPatternPrintStatistics

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Used by the subclasses to signal when they start distributing data, so that the stastics can be updated.
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE scatterPatternStartDistribution(me)
        CLASS(t_ScatterPattern), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":scatterPatternStartDistribution"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        me%curStartTime = p_mpi_wtime()
        me%distributionTime = me%distributionTime - me%curStartTime
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE scatterPatternStartDistribution

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Used by the subclasses to signal when they are done distributing data, so that the stastics can be updated.
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE scatterPatternEndDistribution(me, bytes)
        CLASS(t_ScatterPattern), INTENT(INOUT) :: me
        INTEGER(i8), INTENT(IN) :: bytes

        CHARACTER(*), PARAMETER :: routine = modname//":scatterPatternEndDistribution"
        REAL(dp) :: curEndTime, bw

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        curEndTime = p_mpi_wtime()
        me%distributionTime = me%distributionTime + curEndTime
        me%distributedData = me%distributedData + bytes
        IF(msg_level >= 20 .and. my_process_is_stdio()) THEN
          IF (curEndTime - me%curStartTime == 0) THEN
            bw = -1.0_dp
          ELSE
            bw = REAL(bytes, dp)/(1048576.0_dp*(curEndTime - me%curStartTime))
          END IF
          WRITE(0,*) routine, ": Distributed ", bytes, " bytes in ", &
               curEndTime - me%curStartTime, " seconds (", bw, " MiB/s"
        END IF
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE scatterPatternEndDistribution

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE destructScatterPattern(me)
        CLASS(t_ScatterPattern), TARGET, INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":destructScatterPattern"
        INTEGER :: i

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine

        call me%printStatistics()

        !remove from existing pattern list
        DO i = 1, existingPatternCount
            IF(ASSOCIATED(existingPatterns(i)%ptr, me)) EXIT
        END DO
        existingPatterns(i)%ptr => existingPatterns(existingPatternCount)%ptr
        existingPatternCount = existingPatternCount - 1

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE destructScatterPattern

END MODULE
