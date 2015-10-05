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

  USE mo_kind,                  ONLY: wp
  USE mo_datetime,              ONLY: t_datetime
  USE mo_impl_constants,        ONLY: max_dom
  USE mtime,                    ONLY: max_calendar_str_len
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ini_datetime_string, end_datetime_string, calendar
  PUBLIC :: restart_ini_datetime_string, restart_end_datetime_string, restart_calendar
  PUBLIC :: dt_restart, is_relative_time
  PUBLIC :: t_time_config, time_config

  !> namelist parameters (as raw character strings):
  !
  ! these are the namelist settings originating from the restart file:
  CHARACTER(len=32)                   :: restart_ini_datetime_string
  CHARACTER(len=32)                   :: restart_end_datetime_string
  INTEGER                             :: restart_calendar
  !
  !  these are namelist setting which may originate from the restart
  !  file, but with user modifications in the current run:
  CHARACTER(len=32)                   :: ini_datetime_string
  CHARACTER(len=32)                   :: end_datetime_string
  CHARACTER(len=max_calendar_str_len) :: calendar = ''
  REAL(wp)                            :: dt_restart          !< Length of restart cycle in seconds
  LOGICAL                             :: is_relative_time

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
  !! 
  !! The actual variable
  !!
  TYPE(t_time_config) :: time_config
 
END MODULE mo_time_config

