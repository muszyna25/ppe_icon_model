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

  USE mo_impl_constants,  ONLY: max_char_length
  USE mo_exception,       ONLY: finish
  USE mo_io_units,        ONLY: nnml
  USE mo_namelist,        ONLY: open_nml, close_nml, position_nml, POSITIONED

  USE mo_coupling_config, ONLY: t_field_nml, config_fields

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
  INTEGER            :: frequency
  INTEGER            :: time_step
  INTEGER            :: lag
  CHARACTER(len=132) :: name

  NAMELIST /coupling_nml/ name,                &
                          frequency,           &
                          time_step,           &
                          lag,                 &
                          l_time_average,      &
                          l_time_accumulation

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
  SUBROUTINE read_coupling_namelist (namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    !
    ! Local variables
    !

    TYPE(t_field_nml), POINTER :: new_fields(:)

    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: new_dim

    INTEGER :: nbr_fields
    INTEGER :: istat
    LOGICAL :: first
    LOGICAL :: l_redirect_stdout

    CHARACTER(len=max_char_length), PARAMETER :: &
         &   routine = 'mo_coupling_nml:read_coupling_namelist'

    RETURN

    ! -------------------------------------------------------------------
    ! Allocate space for namelist input
    ! -------------------------------------------------------------------

    nbr_fields = 8

    ALLOCATE(config_fields(nbr_fields))

    !--------------------------------------------------------------------
    ! 1. Set default values
    !--------------------------------------------------------------------

    frequency           = 0
    time_step           = 0

    l_time_average      = .FALSE.
    l_time_accumulation = .FALSE.
    l_redirect_stdout   = .FALSE.

    !--------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !--------------------------------------------------------------------

!rr    IF (is_restart_run()) THEN
!rr      funit = open_and_restore_namelist('coupling_nml')
!rr      READ(funit,NML=coupling_nml)
!rr      CALL close_tmpfile(funit)
!rr    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (done so far by all MPI processes)
    !--------------------------------------------------------------------

    first = .TRUE.
    i     = 0

    !--------------------------------------------------------------------
    ! 3.a loop over occurences of namelist group /coupling_nml/
    !--------------------------------------------------------------------

    DO

       CALL open_nml (TRIM(namelist_filename))
       CALL position_nml('coupling_nml', lrewind=first, status=istat)

       frequency = 0
       time_step = 0
       lag       = 0
       l_time_average      = .FALSE.
       l_time_accumulation = .FALSE.

       first = .FALSE.

       !-----------------------------------------------------------------
       ! 3.b if namelist group is present
       !-----------------------------------------------------------------

       SELECT CASE (istat)

       CASE (POSITIONED)

          READ  (nnml, coupling_nml)

          i = i + 1

          !--------------------------------------------------------------
          ! 3.c allocate more memory if needed
          !--------------------------------------------------------------

          IF ( i > nbr_fields ) THEN

             ! we need to allocated more memory

             new_dim = nbr_fields + 8

             ALLOCATE ( new_fields(new_dim), STAT = ierr )
             IF ( ierr > 0 ) &
                  CALL finish ( TRIM(routine), ' Error allocating fields ' )

             new_fields (1:nbr_fields) = config_fields (1:nbr_fields)

             DEALLOCATE ( config_fields, STAT = ierr )
             IF ( ierr > 0 ) &
                  CALL finish ( TRIM(routine), ' Error deallocating fields ' )

             config_fields => new_fields

             ! update size

             nbr_fields = new_dim

          ENDIF

          !--------------------------------------------------------------
          ! 4. Consistency check
          !--------------------------------------------------------------

          if ( frequency < time_step ) &
               CALL finish (TRIM(routine), 'Coupling frequency must be larger that time step' )

          !--------------------------------------------------------------
          ! 5. Fill the configuration state
          !--------------------------------------------------------------

          config_fields(i)%name                = name
          config_fields(i)%frequency           = frequency
          config_fields(i)%time_step           = time_step
          config_fields(i)%lag                 = lag
          config_fields(i)%l_time_average      = l_time_average
          config_fields(i)%l_time_accumulation = l_time_accumulation

       END SELECT

       IF ( istat /= POSITIONED ) EXIT

    END DO

    CALL close_nml

  END SUBROUTINE read_coupling_namelist

END MODULE mo_coupling_nml
