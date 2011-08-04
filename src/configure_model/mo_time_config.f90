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
  USE mo_datetime,              ONLY: t_datetime
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_time_config, time_config

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !>
  !! Derived type containing information for time control. 
  !!
  TYPE t_time_config

    ! from namelist 

    REAL(wp) :: dt_restart    !< Length of restart cycle in seconds
    INTEGER  :: calendar      !< calendar type

    ! not directly from namelist  

    TYPE(t_datetime) :: ini_datetime  !< Starting time of model integration
    TYPE(t_datetime) :: end_datetime  !< Ending   time of model integration
    TYPE(t_datetime) :: cur_datetime  !< Current  time model time 
 
  END TYPE t_time_config
  !>
  !! 
  !! The actual variable
  !!
  TYPE(t_time_config) :: time_config

!CONTAINS
!
!  SUBROUTINE configure_time
!  END SUBROUTINE configure_time

END MODULE mo_time_config

