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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
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
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_nh_pzlev_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! nh_pzlev_nml namelist variables  !
  !----------------------------------!

  INTEGER :: nzlev                 !< number of z-levels

  INTEGER :: nplev                 !< number of p-levels

  INTEGER :: nilev                 !< number of isentropic-levels

  REAL(wp):: zlevels(100)          !< zlevel heights [m] 

  REAL(wp):: plevels(100)          !< plevel heights [Pa] 

  REAL(wp):: ilevels(100)          !< isentropes [K] 


  NAMELIST/nh_pzlev_nml/ nzlev, nplev, nilev,  &
    &                    zlevels, plevels, ilevels

CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for output on pressure/height levels 
  !!
  !! This subroutine 
  !! - reads the Namelist for output on p/i/z-levels
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
    INTEGER :: jk          !< vertical loop index
    INTEGER :: iunit

    REAL(wp) :: p_bot      !< bottom level pressure
    REAL(wp) :: p_top      !< top level pressure
    REAL(wp) :: delp       !< pressure increment based on p_bot, p_top, nplev

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_pzlev_nml: read_nh_pzlev_namelist'

    !--------------------------------------------------------------

    !-----------------------
    ! 1. default settings   
    !-----------------------
    nzlev       = 10
    nplev       = 10
    nilev       = 3

    ! set height and pressure levels - attention: ordering of the levels must be top-down 
    ! as for the model levels
    !
    ! standard height levels (every 1000 m)
    !
    ! Initialize:
    zlevels(:) = 0._wp
    plevels(:) = 0._wp
    DO jk = 1, nzlev
      zlevels(nzlev+1-jk) = REAL(jk-1,wp)*1000._wp
    ENDDO

    ! standard pressure levels (every (p_bot - p_top)/nplev Pa starting from p_bot)
    !
    p_bot = 100000._wp
    p_top = 0._wp
    delp  = (p_bot - p_top)/REAL(nplev,wp)
    
    plevels(nplev) =  p_bot
    DO jk = 2, nplev
      plevels(nplev+1-jk) = plevels(nplev+2-jk) - delp
    ENDDO

    ! standard isentropes (ordering from TOA to bottom)
    !
    ilevels(:)   = 0._wp
    ilevels(1:3) =(/340._wp, 320._wp, 300._wp/) 


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
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nh_pzlev_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nh_pzlev_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nh_pzlev_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      nh_pzlev_config(jg)%nzlev       = nzlev
      nh_pzlev_config(jg)%nplev       = nplev
      nh_pzlev_config(jg)%nilev       = nilev
      nh_pzlev_config(jg)%zlevels     = zlevels
      nh_pzlev_config(jg)%plevels     = plevels
      nh_pzlev_config(jg)%ilevels     = ilevels
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nh_pzlev_nml)                    
      CALL store_and_close_namelist(funit, 'nh_pzlev_nml')             
    ENDIF


    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nh_pzlev_nml)


  END SUBROUTINE read_nh_pzlev_namelist

END MODULE mo_nh_pzlev_nml
