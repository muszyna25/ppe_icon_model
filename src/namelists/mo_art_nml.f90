!>
!! Namelist for ART-package
!!
!! Subroutine is called by read_atmo_namelists for setting up the ART-package 
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-12-08)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_art_nml

  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: finish
  USE mo_io_units            ,ONLY: nnml, nnml_output
  USE mo_master_control      ,ONLY: is_restart_run
  USE mo_impl_constants      ,ONLY: max_dom
  USE mo_namelist            ,ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_io_restart_namelist ,ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_art_config          ,ONLY: art_config 

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_art_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! art_nml namelist variables       !
  !----------------------------------!

  LOGICAL :: lart                 !< main switch for using the ART-package
                                  !< .TRUE.: switch ON
                                  !<.FALSE.: switch OFF

  NAMELIST/art_nml/ lart


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for ART-package. 
  !!
  !! This subroutine 
  !! - reads the Namelist for the ART-package
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-12-08)
  !!
  SUBROUTINE read_art_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_art_nml: read_art_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    lart    = .FALSE.           ! ART-package switched off


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('art_nml')
      READ(funit,NML=art_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('art_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, art_nml)
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      art_config(jg)%lart   = lart
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=art_nml)                    
      CALL store_and_close_namelist(funit, 'art_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=art_nml)


  END SUBROUTINE read_art_namelist

END MODULE mo_art_nml
