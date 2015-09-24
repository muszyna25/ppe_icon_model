!>
!! Contains the setup of the sleve coordinate
!!
!!        
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_sleve_nml

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_config,       ONLY: isRestart
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_sleve_config       , ONLY: config_min_lay_thckn => min_lay_thckn, &
    &                               config_max_lay_thckn => max_lay_thckn, &
    &                               config_htop_thcknlimit => htop_thcknlimit , &
    &                               config_top_height    => top_height   , &
    &                               config_decay_scale_1 => decay_scale_1, &
    &                               config_decay_scale_2 => decay_scale_2, &
    &                               config_decay_exp     => decay_exp    , &
    &                               config_flat_height   => flat_height  , &
    &                               config_stretch_fac   => stretch_fac  , &
    &                               config_lread_smt     => lread_smt 
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: read_sleve_namelist

  !----------------------------------------------------------------------------
  ! Namelist variables for the SLEVE coordinate
  !---------------------------------------------------------------------------
  !
  ! a) Parameters specifying the distrubution of the coordinate surfaces
  REAL(wp):: min_lay_thckn   ! Layer thickness of lowermost level
  REAL(wp):: max_lay_thckn   ! Maximum layer thickness below htop_thcknlimit
  REAL(wp):: htop_thcknlimit ! Height below which the layer thickness must not exceed max_lay_thckn
  REAL(wp):: stretch_fac     ! Factor for stretching/squeezing the model layer distribution
  REAL(wp):: top_height      ! Height of model top

  ! b) Parameters for SLEVE definition
  REAL(wp):: decay_scale_1    ! Decay scale for large-scale topography component
  REAL(wp):: decay_scale_2    ! Decay scale for small-scale topography component
  REAL(wp):: decay_exp        ! Exponent for decay function
  REAL(wp):: flat_height      ! Height above which the coordinate surfaces are exactly flat
                              ! (not available in the standard SLEVE definition)

  ! c) Parameter for reading in smoothed topography
  LOGICAL :: lread_smt
 
  NAMELIST /sleve_nml/ min_lay_thckn, max_lay_thckn, htop_thcknlimit, top_height,         &
                       decay_scale_1, decay_scale_2, decay_exp, flat_height, stretch_fac, &
                       lread_smt

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for SLEVE coordinate. 
  !!
  !! This subroutine 
  !! - reads the Namelist for SLEVE coordinate
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-07-06)
  !!
  SUBROUTINE read_sleve_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    !0!CHARACTER(len=*), PARAMETER ::  &
    !0!  &  routine = 'mo_sleve_nml:read_sleve_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! a) Parameters determining the distribution of model layers
    !    (if not read in from a table)
    min_lay_thckn   = 50._wp      ! Layer thickness of lowermost layer
    max_lay_thckn   = 25000._wp   ! Maximum layer thickness below htop_thcknlimit
    htop_thcknlimit = 15000._wp   ! Height below which the layer thickness must not exceed max_lay_thckn
    top_height      = 23500._wp   ! Height of model top
    stretch_fac     = 1._wp       ! Scaling factor for stretching/squeezing 
                                  ! the model layer distribution

    ! b) Parameters setting up the decay function of the topographic signal
    decay_scale_1   = 4000._wp    ! Decay scale of large-scale topography component
    decay_scale_2   = 2500._wp    ! Decay scale of small-scale topography component
    decay_exp       = 1.2_wp      ! Exponent for decay function
    flat_height     = 16000._wp   ! Height above which the coordinate surfaces are 
                                      ! flat

    ! c) parameter to switch on/off internal topography smoothing
    lread_smt       = .FALSE.     ! read smoothed topography from file (TRUE/FALSE)

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('sleve_nml')
      READ(funit,NML=sleve_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('sleve_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, sleve_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, sleve_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, sleve_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    config_min_lay_thckn = min_lay_thckn
    config_max_lay_thckn = max_lay_thckn
    config_htop_thcknlimit = htop_thcknlimit
    config_top_height    = top_height
    config_decay_scale_1 = decay_scale_1
    config_decay_scale_2 = decay_scale_2
    config_decay_exp     = decay_exp
    config_flat_height   = flat_height
    config_stretch_fac   = stretch_fac
    config_lread_smt     = lread_smt

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=sleve_nml)                    
      CALL store_and_close_namelist(funit, 'sleve_nml') 
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=sleve_nml)

  END SUBROUTINE read_sleve_namelist

END MODULE mo_sleve_nml
