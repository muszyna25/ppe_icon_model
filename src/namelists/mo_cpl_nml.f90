!>
!!        Contains the variables to set up the coupling.
!!
!!        
!! @par Revision History
!!   Created by Rene Redler (2011-03-22)
!!
!! @par Copyright
!! 2010-2011 by DWD and MPI-M
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
!!

MODULE mo_cpl_nml

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------

  USE mo_kind,     ONLY : wp
  USE mo_datetime, ONLY : proleptic_gregorian, &
   &                      date_to_time,        &
   &                      check_date,          &
   &                      print_datetime_all

  USE mo_icon_cpl, ONLY : complist,                           &
   &                      l_debug,                            &
   &                      initial_date, final_date,           &
   &                      nbr_ICON_comps

  USE mo_io_units, ONLY: nnml

  IMPLICIT NONE

  LOGICAL                              :: l_redirect_stdout = .false.

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PRIVATE

  PUBLIC :: cpl_nml_setup, l_redirect_stdout

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! initialization
  ! --------------

  LOGICAL            :: l_time_average       = .false.
  LOGICAL            :: l_time_accumulation  = .false.

  INTEGER            :: coupling_freq = 0
  INTEGER            :: time_step     = 0

  CHARACTER(LEN=22)  :: start_date = ""
  CHARACTER(LEN=22)  :: end_date   = ""


  NAMELIST /icon_cpl/ start_date,          &
                      end_date,            &
                      coupling_freq,       &
                      time_step,           &
                      l_time_average,      &
                      l_time_accumulation, &
                      l_redirect_stdout

  CONTAINS

 !>
 !!  Initialization of variables that contain general information.
 !!
 !!               Initialization of variables that contain general information
 !!               about the coupled model run. The configuration is read from
 !!               namelist 'icon_cpl'.
 !!
 !! @par Revision History
 !!
 SUBROUTINE cpl_nml_setup(index)

   INTEGER, INTENT(in) :: index
                                               
!
! !Local variables
!
   INTEGER          :: istat

   CHARACTER(len=4) :: year_str
   CHARACTER(len=2) :: month_str
   CHARACTER(len=2) :: day_str
   CHARACTER(len=2) :: hour_str
   CHARACTER(len=2) :: minute_str
   CHARACTER(len=2) :: second_str
   CHARACTER(len=2) :: pm_str ! AC | BC [default is AC]

   INTEGER          :: second
   INTEGER          :: i

!rr   CHARACTER(len=max_char_length), PARAMETER :: &
!rr     &  routine = 'mo_cpl_nml/cpl_nml_setup'

!-----------------------------------------------------------------------
   
   !------------------------------------------------------------
   ! 2.0 set up the default values for icon_cpl
   !------------------------------------------------------------

    initial_date%year   = 1850
    initial_date%month  = 1
    initial_date%day    = 1
    initial_date%hour   = 0
    initial_date%minute = 0
    initial_date%second = 0.0_wp

    final_date%calendar = proleptic_gregorian
    final_date%year     = 1950
    final_date%month    = 1
    final_date%day      = 1
    final_date%hour     = 0
    final_date%minute   = 0
    final_date%second   = 0.0_wp

    start_date(1:22)    = '1900-01-01 00:00:00 AC'
    end_date  (1:22)    = '1950-01-01 00:00:00 AC'
    coupling_freq       = 0
    time_step           = 0

    l_time_average      = .FALSE.
    l_time_accumulation = .FALSE.
    l_redirect_stdout   = .FALSE.

   !------------------------------------------------------------
   ! 3.0 Read the namelist to evaluate if a restart file shall
   !     be used to configure and initialize the model, and if
   !     so take read the run_ctl parameters from the restart
   !     file.
   !------------------------------------------------------------
   ! (done so far by all MPI processes)

   OPEN (nnml, file=complist(index)%nml_name, iostat=istat, status='old', &
        action='read', delim='apostrophe')

   READ (nnml, icon_cpl)

    !
    ! set start date of the coupled run
    ! ---------------------------------
    !
    year_str   = start_date(1:4)
    month_str  = start_date(6:7)
    day_str    = start_date(9:10)
    hour_str   = start_date(12:13)
    minute_str = start_date(15:16)
    second_str = start_date(18:19)
    pm_str     = start_date(21:22)

    initial_date%calendar = proleptic_gregorian

    READ(year_str,   '(I4)') initial_date%year
    READ(month_str,  '(I2)') initial_date%month
    READ(day_str,    '(I2)') initial_date%day
    READ(hour_str,   '(I2)') initial_date%hour
    READ(minute_str, '(I2)') initial_date%minute
    READ(second_str, '(I2)') second

    CALL date_to_time ( initial_date )

    !
    ! set end date of the coupled run
    ! -------------------------------
    !
    year_str   = end_date(1:4)
    month_str  = end_date(6:7)
    day_str    = end_date(9:10)
    hour_str   = end_date(12:13)
    minute_str = end_date(15:16)
    second_str = end_date(18:19)
    pm_str     = end_date(21:22)

    final_date%calendar = proleptic_gregorian

    READ(year_str,   '(i4)') final_date%year
    READ(month_str,  '(I2)') final_date%month
    READ(day_str,    '(I2)') final_date%day
    READ(hour_str,   '(I2)') final_date%hour
    READ(minute_str, '(I2)') final_date%minute
    READ(second_str, '(I2)') second

    final_date%second = float(second)

    CALL date_to_time ( final_date )

    IF ( l_debug ) THEN
       CALL check_date         ( initial_date )
       CALL check_date         ( final_date   )
!      CALL print_datetime_all ( initial_date )
!      CALL print_datetime_all ( final_date   )
!      WRITE ( cplout , '(A,A2)' ) 'Before/after Christ: ', pm_str
    ENDIF

    CLOSE (nnml, IOSTAT=istat)

    ! -------------------------------------------------------------------
    ! Assign component namelist input
    ! -------------------------------------------------------------------

    DO i = 1, nbr_ICON_comps

       complist(i)%l_time_average = l_time_average
       complist(i)%l_time_accumulation = l_time_accumulation
       complist(i)%coupling_freq  = coupling_freq
       complist(i)%time_step      = time_step

    ENDDO


 END SUBROUTINE cpl_nml_setup

END MODULE mo_cpl_nml
