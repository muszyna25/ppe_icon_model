!>
!!        
!! @par Revision History
!!       Kristina Froehlich, MPI-M 2011-07-05
!! first implementation for all time control variables
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_time_nml

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_datetime,           ONLY: proleptic_gregorian
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_master_nml,         ONLY: lrestart
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_io_restart_namelist,ONLY: open_and_restore_namelist, close_tmpfile,&
                                  & open_tmpfile, store_and_close_namelist

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  
  PUBLIC :: time_ctl, nml_calendar, nml_ini_datetime,&
    &                  nml_end_datetime, nml_dt_restart

  PRIVATE

  ! time information
  ! ----------------
  !
  ! calendar type
  INTEGER            :: nml_calendar
  !
  ! - data and time of model start and end
  CHARACTER(len=16) :: nml_ini_datetime
  CHARACTER(len=16) :: nml_end_datetime

  ! restart interval
  ! ----------------
  REAL(wp)           :: nml_dt_restart          ! [s] length of a restart cycle 


  NAMELIST /time_ctl/nml_calendar,                  &
    &                nml_ini_datetime, nml_end_datetime,&
    &                nml_dt_restart

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !!  Initialization of variables that contain general information.
  !!
  !!  Initialization of variables that contain general information
  !!  about the model run. The configuration is read from
  !!  namelist 'time_ctl'.
  !!
  !! @par Revision History
  !!  Reading of the namelist and checking of the validity were
  !!  in some other modules in the earlier version of the shallow water model.
  !!  Moved to this module by Hui Wan, MPI-M (2007-02-23)
  !!  The character parameter <i>routine</i> was introduced and used
  !!  for error information by Hui Wan, MPI-M (2007-02-23).
  !!  Modified by Almut Gassmann, MPI-M (2008-09-23)
  !!  - introduced i_cell_type, lidealized and lshallow_water
  !!  Modified by Marco Giorgetta, MPI-M (2009-02-23)
  !!  - lidealized replaced by ltestcase
  !!  Modification by Constantin Junk, MPI-M (2010-02-22)
  !!  - changes to consistency checks
  !!
  SUBROUTINE time_nml_read
                                               
   INTEGER  :: istat, funit

   CHARACTER(len=max_char_length), PARAMETER ::   &
            &  routine = 'mo_time_nml/time_nml_setup'

   ! initial date and time
   nml_calendar       = proleptic_gregorian
   nml_ini_datetime   = "20080901T000000Z"
   !
   ! end date and time
   nml_end_datetime   = "20080901T014000Z"
   !
   ! length of integration = (number of timesteps)*(length of timestep)
   ! - If nsteps is set to a non-zero positive value, then the end date is computed
   !   from the initial date and time, the time step dtime, and nsteps.
   ! - Else if run_day, run_hour, run_minute or run_second is set to a non-zero,
   !   positive value, then the initial date and time and the run_... variables are
   !   used to compute the end date and time and, using dtime, nsteps.
   !   Else nsteps is computed from the initial and end date and time and dtime.
   !
   ! length of restart cycle
   nml_dt_restart     = 86400._wp*30._wp   ! = 30 days
   !
  IF (lrestart) THEN      
 
   ! 2.1 Overwrite the defaults above by values in the restart file

      funit = open_and_restore_namelist('time_ctl')
      READ(funit,NML=run_ctl)
      CALL close_tmpfile(funit) 
  END IF
   !------------------------------------------------------------------------
   !  Read user's (new) specifications. (Done so far by all MPI processes)
   !------------------------------------------------------------------------
    CALL position_nml('time_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, time_ctl)
    END SELECT


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=time_ctl)
    CALL store_and_close_namelist(funit, 'time_ctl')

!    ! write the contents of the namelist to an ASCII file
!    IF(p_pe == p_io) WRITE(nnml_output,nml=time_ctl)


 END SUBROUTINE time_nml_read

END MODULE mo_time_nml
