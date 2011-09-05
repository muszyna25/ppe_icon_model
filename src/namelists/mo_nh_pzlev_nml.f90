!>
!! Namelist for output on pressure/height levels
!!
!! These subroutines are called by read_atmo_namelists and do the setup 
!! for output on pressure and/or height levels.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-09-05)
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
MODULE mo_nh_pzlev_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nh_pzlev_config,     ONLY: nh_pzlev_config
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_nh_pzlev_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! nh_pzlev_nml namelist variables  !
  !----------------------------------!

  LOGICAL :: lwrite_zlev           !< if .TRUE.: write output on z-levels

  LOGICAL :: lwrite_plev           !< if .TRUE.: write output on p-levels

  INTEGER :: nzlev                 !< number of z-levels

  INTEGER :: nplev                 !< number of p-levels

  INTEGER :: zlevels(100)          !< zlevel heights [m] 

  INTEGER :: plevels(100)          !< plevel heights [m] 


  NAMELIST/nh_pzlev_nml/ lwrite_zlev, lwrite_plev, nzlev, nplev,   &
    &                    zlevels, plevels

CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for output on pressure/height levels 
  !!
  !! This subroutine 
  !! - reads the Namelist for output on p/z-levels
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-09-05)
  !!
  SUBROUTINE read_nh_pzlev_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: jk          !< tracer loop index

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_pzlev_nml: read_nh_pzlev_namelist'

    !--------------------------------------------------------------

    !-----------------------
    ! 1. default settings   
    !-----------------------
    lwrite_zlev = .TRUE.
    lwrite_plev = .TRUE.
    nzlev       = 10
    nplev       = 10

    ! set levels - attention: ordering of the levels must be top-down 
    ! as for the model levels
    DO jk = 1, nzlev
      zlevels(nzlev+1-jk) = REAL(jk-1,wp)*1000._wp ! every 1000 m
    ENDDO
    ! standard pressure levels
    plevels(10) = 100000._wp
    plevels(9)  =  92500._wp
    plevels(8)  =  85000._wp
    plevels(7)  =  70000._wp
    plevels(6)  =  50000._wp
    plevels(5)  =  40000._wp
    plevels(4)  =  30000._wp
    plevels(3)  =  25000._wp
    plevels(2)  =  20000._wp
    plevels(1)  =  10000._wp



    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('nh_pzlev_nml')
      READ(funit,NML=nh_pzlev_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nh_pzlev_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nh_pzlev_nml)
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    ! check whether lwrite_zlev=.TRUE. when lwrite_plev=.TRUE. is chosen
    !
    IF ( lwrite_plev .AND. .NOT. lwrite_zlev)    THEN
      CALL finish( TRIM(routine),                                     &
        &  'lwrite_plev=.TRUE. only possible if lwrite_zlev=.TRUE.')
    ENDIF



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      nh_pzlev_config(jg)%lwrite_zlev = lwrite_zlev
      nh_pzlev_config(jg)%lwrite_plev = lwrite_plev
      nh_pzlev_config(jg)%nzlev       = nzlev
      nh_pzlev_config(jg)%nplev       = nplev
      nh_pzlev_config(jg)%zlevels     = zlevels
      nh_pzlev_config(jg)%plevels     = plevels
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nh_pzlev_nml)                    
      CALL store_and_close_namelist(funit, 'transport_nml')             
    ENDIF


    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nh_pzlev_nml)


  END SUBROUTINE read_nh_pzlev_namelist

END MODULE mo_nh_pzlev_nml
