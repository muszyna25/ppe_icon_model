!>
!! @brief namelist setup for JSBACH land scheme
!!
!! Namelist setup for JSBACH land scheme
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Reiner Schnur, MPI-M (2012-04-16)
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2012 by DWD and MPI-M
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
MODULE mo_lnd_jsbach_nml

  USE mo_lnd_jsbach_config,   ONLY: lnd_jsbach_config
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_exception,           ONLY: finish
  USE mo_mpi,                 ONLY: my_process_is_stdio

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_lnd_jsbach_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  INTEGER :: ntiles           !< Number of tiles
  INTEGER :: nsoil            !< Number of soil layers
  INTEGER :: ntsoil           !< Number of soil layers for soil temperature

  NAMELIST /jsbach_ctl_nml/ ntiles, nsoil, ntsoil

CONTAINS
  !>
  !!
  SUBROUTINE read_lnd_jsbach_namelist(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_lnd_jsbach_nml:read_lnd_jsbach_namelist'

    !------------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------------
    ntiles = 12
    nsoil  = 1
    ntsoil = 5

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('jsbach_ctl_nml')
      READ(funit,NML=jsbach_ctl_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------
    ! Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('jsbach_ctl_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, jsbach_ctl_nml)
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! Sanity Check
    !------------------------------------------------------------------
    IF (ntiles /= 12 ) THEN
      CALL finish(TRIM(routine),           &
        &  'Currently, ntiles must be 12.')
    END IF

    IF (nsoil /= 1 .AND. nsoil /= 5) THEN
      CALL finish(TRIM(routine),           &
        &  'nsoil must be 1 or 5.')
    END IF

    IF (ntsoil /= 5) THEN
      CALL finish(TRIM(routine),           &
        &  'ntsoil must be 5.')
    END IF

    !------------------------------------------------------------------
    ! Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=jsbach_ctl_nml)
      CALL store_and_close_namelist(funit, 'jsbach_ctl_nml')
    ENDIF

    !------------------------------------------------------------------
    ! Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=jsbach_ctl_nml)

    !------------------------------------------------------------------
    ! Fill the configuration state
    !------------------------------------------------------------------
    lnd_jsbach_config% ntiles    = ntiles
    lnd_jsbach_config% nsoil     = nsoil
    lnd_jsbach_config% ntsoil    = ntsoil

  END SUBROUTINE read_lnd_jsbach_namelist

END MODULE mo_lnd_jsbach_nml
