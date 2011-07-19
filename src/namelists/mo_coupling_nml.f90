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

MODULE mo_coupling_nml

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------

  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH
  USE mo_exception, ONLY: warning, message_text, finish
  USE mo_io_units,  ONLY: filename_max, nnml
  USE mo_namelist,  ONLY: open_nml, position_nml, POSITIONED

  USE mo_icon_cpl,  ONLY : complist,                           &
       &                   l_debug,                            &
       &                   initial_date, final_date,           &
       &                   nbr_ICON_comps


  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  ! ------------------------------------------------------------------------
  ! Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! initialization
  ! --------------

  LOGICAL            :: l_time_average
  LOGICAL            :: l_time_accumulation
  LOGICAL            :: l_redirect_stdout
  INTEGER            :: coupling_freq
  INTEGER            :: time_step

  CHARACTER(LEN=132)  :: start_date = ""
  CHARACTER(LEN=132)  :: end_date   = ""


  NAMELIST /coupling_nml/ coupling_freq,       &
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
  INTEGER FUNCTION read_coupling_namelist (namelist_filename, comp_id)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    INTEGER, INTENT(in)          :: comp_id

    !
    ! Local variables
    !

    INTEGER :: i, istat

    CHARACTER(len=max_char_length), PARAMETER :: &
         &   routine = 'mo_coupling_nml:read_coupling_namelist'

    !-----------------------------------------------------------------------

    !------------------------------------------------------------
    ! 1. Set default values
    !------------------------------------------------------------

    start_date(1:22)    = '1900-01-01T00:00:00Z'
    end_date  (1:22)    = '1950-01-01T00:00:00Z'
    coupling_freq       = 0
    time_step           = 0

    l_time_average      = .FALSE.
    l_time_accumulation = .FALSE.
    l_redirect_stdout   = .FALSE.

    !------------------------------------------------------------------
    ! 2. Read user's specifications (done so far by all MPI processes)
    !------------------------------------------------------------------

    OPEN (nnml, FILE=TRIM(namelist_filename), IOSTAT=istat, &
         & STATUS='old', ACTION='read', DELIM='apostrophe')

    IF (istat/=0) THEN
       CALL warning(namelist_filename,"not found")
       read_coupling_namelist=-1
       RETURN
    ENDIF

    CALL position_nml('coupling_nml',STATUS=istat)

    IF (istat/=POSITIONED) THEN

       CALL finish( TRIM(routine), &
            & 'Namelist coupling_nml not found in file '  &
            & //TRIM(namelist_filename) )

       read_coupling_namelist=-2
       RETURN      
    ENDIF

    READ (nnml, coupling_nml)
    CLOSE (nnml, IOSTAT=istat)

    ! -------------------------------------------------------------------
    ! Assign component namelist input
    ! -------------------------------------------------------------------

    complist(comp_id)%l_time_average      = l_time_average
    complist(comp_id)%l_time_accumulation = l_time_accumulation
    complist(comp_id)%coupling_freq       = coupling_freq
    complist(comp_id)%time_step           = time_step
    complist(comp_id)%l_redirect_stdout   = l_redirect_stdout

  END FUNCTION read_coupling_namelist

END MODULE mo_coupling_nml
