!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Hui Wan, MPI-M (2011-07-15)
!! @author Kristina Froehlich, MPI-M (2011-07-15)
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_time_config

    USE mo_datetime,              ONLY: t_datetime
    USE mo_impl_constants,        ONLY: max_dom
    USE mo_kind,                  ONLY: wp
    USE mtime,                    ONLY: datetime, newDatetime, deallocateDatetime

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: t_time_config, time_config, restart_experiment

    PUBLIC :: getIniTime


    !>
    !! Derived type containing information for time control.
    !!
    TYPE t_time_config

        ! from namelist

        REAL(wp)         :: dt_restart         !< Length of restart cycle in seconds
        INTEGER          :: calendar           !< calendar type

        ! not directly from namelist

        TYPE(t_datetime) :: ini_datetime       !< Starting time of model integration
        TYPE(t_datetime) :: end_datetime       !< Ending   time of model integration
        TYPE(t_datetime) :: cur_datetime       !< Current  time model time

        REAL(wp)         :: sim_time(max_dom)  !< elapsed simulation time (may locally differ between domains!)

        !> LOGICAL is_relative_time: .TRUE., if time loop shall start with
        !> step 0 regardless whether we are in a standard run or in a
        !> restarted run (which means re-initialized run):
        LOGICAL          ::  is_relative_time

    END TYPE t_time_config

    !>
    !! The actual variable
    !!
    TYPE(t_time_config) :: time_config

    LOGICAL :: restart_experiment ! if true, we have to restart the experiment

CONTAINS

    FUNCTION getIniTime() RESULT(RESULT)
        TYPE(datetime) :: RESULT

        TYPE(datetime), POINTER :: tempTime

        ! get ini-datetime in mtime-format
        tempTime => newDatetime(time_config%ini_datetime%year,       &
          &                     time_config%ini_datetime%month,      &
          &                     time_config%ini_datetime%day,        &
          &                     time_config%ini_datetime%hour,       &
          &                     time_config%ini_datetime%minute,     &
          &                     INT(time_config%ini_datetime%second),&
          &                     ms=0)
        RESULT = tempTime
        CALL deallocateDatetime(tempTime)
    END FUNCTION getIniTime

END MODULE mo_time_config
