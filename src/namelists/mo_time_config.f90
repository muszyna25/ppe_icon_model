!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_time_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_datetime,              ONLY: t_datetime, date_to_time, add_time, &
    &                                 print_datetime_all
  USE mo_master_nml,            ONLY: lrestart
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
 

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'



 TYPE t_time_config

  ! calendar type
  INTEGER          :: calendar

  TYPE(t_datetime) :: ini_datetime  !< Starting time of model integration
!
  TYPE(t_datetime) :: end_datetime  !< Ending   time of model integration

  TYPE(t_datetime) :: current_datetime  !< Current  time model time 

  ! - data and time structure
  !
  ! current model time, not a namelist variable
  !

  REAL(wp) :: dt_restart    !< Length of restart cycle in seconds

 END TYPE t_time_config

 TYPE(t_time_config) :: time_config


  !>
  !! Derived type containing variables for time control. 
  !!
  !

  !! HW Comment: the character-type variables containing ini_ and end_time
  !! in the format "YYYYMMDDTHHMMSSZ" should be namelist variables
  !! and used for computing ini/end_datetime in this type.



CONTAINS

SUBROUTINE time_setup

    INTEGER  :: calendar      !< calendar type 
  ! time information
  ! ----------------

  ! - data and time structure

  INTEGER:: calendar_old
  INTEGER  :: ini_year_old, ini_month_old, ini_day_old, ini_hour_old, ini_minute_old
  INTEGER  :: restart_year, restart_month, restart_day, restart_hour, restart_minute
  REAL(wp) :: ini_second_old
  REAL(wp) :: restart_second
  REAL(wp) :: cur_datetime_calsec, end_datetime_calsec, length_sec
  
! restart interval
  ! ----------------
    CHARACTER(LEN=132) :: routine = 'time_setup'





END SUBROUTINE time_setup

END MODULE mo_time_config

