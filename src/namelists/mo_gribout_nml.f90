!>
!! Namelist for Grib output
!!
!! These subroutines are called by  read_atmo_namelists and do the transport 
!! setup.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-04-20)
!! - moved here from mo_advection_utils
!! Modification by Daniel Reinert, DWD (2011-04-20)
!! - some updates on the structure
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
MODULE mo_gribout_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_gribout_config,      ONLY: gribout_config 

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_gribout_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! gribout_nml namelist variables   !
  !----------------------------------!

  INTEGER :: &                          ! Table 1.2
    & significanceOfReferenceTime       ! 0: Analysis
                                        ! 1: Start of forecast
                                        ! 2: Verifying time of forecast
                                        ! 4: ... 

  INTEGER :: &                          ! Table 1.3
    & productionStatusOfProcessedData   ! 0: Oper. products
                                        ! 1: Oper. test products
                                        ! 2: Research products
                                        ! 3: ...

  INTEGER :: &                          ! Table 1.4
    & typeOfProcessedData               ! 0: Analysis products
                                        ! 1: Forecast products
                                        ! 2: Analysis and forecast products
                                        ! 3: ...

  INTEGER :: &                          ! Table 4.3
    & typeOfGeneratingProcess           ! 0: Analysis
                                        ! 1: Initialization
                                        ! 2: Forecast
                                        ! 3: ...

  INTEGER :: &                          ! Table: backgroundProcess
    & backgroundProcess                 ! 0: main run
                                        ! 1: pre-assimilation
                                        ! 2: assimilation
                                        ! 3: ... 

  INTEGER :: &                          ! Table: generatingProcessIdentifier
    & generatingProcessIdentifier(0:max_dom) ! 1: icogl
                                        ! 2: icrgl
                                        ! 3: icoeu
                                        ! 4: ...



  NAMELIST/gribout_nml/  significanceOfReferenceTime,     &
    &                    productionStatusOfProcessedData, &
    &                    typeOfProcessedData,             &
    &                    typeOfGeneratingProcess,         &
    &                    backgroundProcess,               &
    &                    generatingProcessIdentifier


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for gribout. 
  !!
  !! This subroutine 
  !! - reads the Namelist for gribout
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2013-01-29)
  !!
  SUBROUTINE read_gribout_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_gribout_nml: read_gribout_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    significanceOfReferenceTime     = 1  ! 1: Start of forecast
    productionStatusOfProcessedData = 1  ! 1: Oper. test products
    typeOfProcessedData             = 1  ! 1: Forecast products
    typeOfGeneratingProcess         = 2  ! 2: Forecast
    backgroundProcess               = 0  ! 0: main run
    generatingProcessIdentifier(:)  = 1  ! 1: icogl



    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('gribout_nml')
      READ(funit,NML=gribout_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gribout_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gribout_nml)
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------




    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      gribout_config(jg)%significanceOfReferenceTime     = &
        &                significanceOfReferenceTime
      gribout_config(jg)%productionStatusOfProcessedData = &
        &                productionStatusOfProcessedData
      gribout_config(jg)%typeOfProcessedData             = &
        &                typeOfProcessedData
      gribout_config(jg)%typeOfGeneratingProcess         = &
        &                typeOfGeneratingProcess
      gribout_config(jg)%backgroundProcess               = &
        &                backgroundProcess
      gribout_config(jg)%generatingProcessIdentifier     = &
        &                generatingProcessIdentifier(jg)
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=gribout_nml)                    
      CALL store_and_close_namelist(funit, 'gribout_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=gribout_nml)


  END SUBROUTINE read_gribout_namelist

END MODULE mo_gribout_nml
