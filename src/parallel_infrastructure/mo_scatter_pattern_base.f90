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
    USE mo_kind, ONLY: wp, dp, sp, i8
    USE mo_mpi, ONLY: my_process_is_stdio, p_mpi_wtime
    USE mo_run_config, ONLY: msg_level

    IMPLICIT NONE

PUBLIC :: t_scatterPattern, constructScatterPattern, destructScatterPattern, deleteScatterPattern

    TYPE, ABSTRACT :: t_scatterPattern
        INTEGER :: myPointCount !< number of points needed by this PE
        INTEGER :: communicator !< the communicator to use
        INTEGER(i8) :: distributedData  !< statistic on how much data was distributed
        REAL(dp) :: curStartTime    !< the time when the last distribution call was started
        REAL(dp) :: distributionTime    !< statistic on how long we took to distribute the data

    CONTAINS
        PROCEDURE(interface_distribute_dp), DEFERRED :: distribute_dp   !< distribute double precision data
        PROCEDURE(interface_distribute_sp), DEFERRED :: distribute_sp   !< distribute single precision data
        PROCEDURE(interface_distribute_int), DEFERRED :: distribute_int   !< distribute single precision data


        PROCEDURE :: construct => constructScatterPattern   !< constructor

        PROCEDURE :: resetStatistics => scatterPatternResetStatistics   !< reset the internal statistics
        PROCEDURE :: printStatistics => scatterPatternPrintStatistics   !< print the current statistics and reset them

        PROCEDURE :: startDistribution => scatterPatternStartDistribution   !For use by subclasses only.
        PROCEDURE :: endDistribution => scatterPatternEndDistribution   !For use by subclasses only.

        PROCEDURE :: destruct => destructScatterPattern !< destructor

        GENERIC :: distribute => distribute_dp, distribute_sp, distribute_int
    END TYPE

PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_grid_distribution_base"
    LOGICAL, PARAMETER :: debugModule = .false.

    ABSTRACT INTERFACE
        !---------------------------------------------------------------------------------------------------------------------------
        !> do the data distribution for a double precision array
        !---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE interface_distribute_dp(me, globalArray, localArray, ladd_value)
            IMPORT t_scatterPattern, dp, wp
            CLASS(t_scatterPattern), INTENT(INOUT) :: me
            REAL(dp), INTENT(INOUT) :: globalArray(:)
            REAL(wp), INTENT(INOUT) :: localArray(:,:)
            LOGICAL, INTENT(IN) :: ladd_value
        END SUBROUTINE interface_distribute_dp
        !---------------------------------------------------------------------------------------------------------------------------
        !> do the data distribution for a double precision array
        !---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE interface_distribute_sp(me, globalArray, localArray, ladd_value)
            IMPORT t_scatterPattern, sp, wp
            CLASS(t_scatterPattern), INTENT(INOUT) :: me
            REAL(sp), INTENT(INOUT) :: globalArray(:)
            REAL(wp), INTENT(INOUT) :: localArray(:,:)
            LOGICAL, INTENT(IN) :: ladd_value
        END SUBROUTINE interface_distribute_sp
        !---------------------------------------------------------------------------------------------------------------------------
        !> do the data distribution for a integer
        !---------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE interface_distribute_int(me, globalArray, localArray, ladd_value)
            IMPORT t_scatterPattern
            CLASS(t_scatterPattern), INTENT(INOUT) :: me
            INTEGER, INTENT(INOUT) :: globalArray(:)
            INTEGER, INTENT(INOUT) :: localArray(:,:)
            LOGICAL, INTENT(IN) :: ladd_value
        END SUBROUTINE interface_distribute_int
    END INTERFACE

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destruct and deallocate a t_scatterPattern
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE deleteScatterPattern(me)
        CLASS(t_scatterPattern), POINTER, INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":deleteScatterPattern"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        CALL me%destruct()
        DEALLOCATE(me)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE deleteScatterPattern

    !-------------------------------------------------------------------------------------------------------------------------------
    !> constructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE constructScatterPattern(me, loc_arr_len, glb_index, communicator)
        CLASS(t_scatterPattern), INTENT(OUT) :: me
        INTEGER, INTENT(IN) :: loc_arr_len, glb_index(:), communicator

        CHARACTER(*), PARAMETER :: routine = modname//":constructScatterPattern"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        me%myPointCount = loc_arr_len
        me%communicator = communicator
        me%distributedData = 0
        me%distributionTime = 0.0
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE constructScatterPattern

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Reset the internal statistics of how long it took to distribute how much data.
    !> This is used in conjunction with printStatistics() to produce statistics output for only a subset of distribute() calls.
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE scatterPatternResetStatistics(me)
        CLASS(t_scatterPattern), INTENT(INOUT) :: me

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
        CLASS(t_scatterPattern), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":scatterPatternPrintStatistics"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        IF(me%distributedData > 0) THEN
            IF(msg_level >= 10 .and. my_process_is_stdio()) THEN
                WRITE(0,*) routine, ": data distribution totals:"
                WRITE(0,'(8X,A,I19,A)')   "amount:    ", me%distributedData, " bytes"
                WRITE(0,'(8X,A,F19.3,A)') "duration:  ", me%distributionTime, " seconds"
                WRITE(0,'(8X,A,F19.3,A)') "bandwidth: ", REAL(me%distributedData, dp)/(1048576.0_dp*me%distributionTime), " MiB/s"
            END IF
        END IF
        call me%resetStatistics()
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE scatterPatternPrintStatistics

    !-------------------------------------------------------------------------------------------------------------------------------
    !> Used by the subclasses to signal when they start distributing data, so that the stastics can be updated.
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE scatterPatternStartDistribution(me)
        CLASS(t_scatterPattern), INTENT(INOUT) :: me

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
        CLASS(t_scatterPattern), INTENT(INOUT) :: me
        INTEGER(i8), INTENT(IN) :: bytes

        CHARACTER(*), PARAMETER :: routine = modname//":scatterPatternEndDistribution"
        REAL(dp) :: curEndTime

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        curEndTime = p_mpi_wtime()
        me%distributionTime = me%distributionTime + curEndTime
        me%distributedData = me%distributedData + bytes
        IF(msg_level >= 20 .and. my_process_is_stdio()) THEN
            WRITE(0,*) routine, ": Distributed ", bytes, " bytes in ", curEndTime - me%curStartTime, " seconds (", &
                &      REAL(bytes, dp)/(1048576.0_dp*(curEndTime - me%curStartTime)), " MiB/s)"
        END IF
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE scatterPatternEndDistribution

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE destructScatterPattern(me)
        CLASS(t_scatterPattern), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":destructScatterPattern"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        call me%printStatistics()
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE destructScatterPattern

END MODULE
