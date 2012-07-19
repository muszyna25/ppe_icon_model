!>
!! Contains the setup of the sleve coordinate
!!
!!        
!! @par Revision History
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
MODULE mo_sleve_nml

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_sleve_config       , ONLY: config_min_lay_thckn => min_lay_thckn, &
    &                               config_top_height    => top_height   , &
    &                               config_decay_scale_1 => decay_scale_1, &
    &                               config_decay_scale_2 => decay_scale_2, &
    &                               config_decay_exp     => decay_exp    , &
    &                               config_flat_height   => flat_height  , &
    &                               config_stretch_fac   => stretch_fac  , &
    &                               config_lread_smt     => lread_smt 

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: read_sleve_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !----------------------------------------------------------------------------
  ! Namelist variables for the SLEVE coordinate
  !---------------------------------------------------------------------------
  !
  ! a) Parameters specifying the distrubution of the coordinate surfaces
  !     (the initializations are a workaround for a NEC compiler bug)
  REAL(wp):: min_lay_thckn = 1._wp  ! Layer thickness of lowermost level
  REAL(wp):: stretch_fac   = 1._wp  ! Factor for stretching/squeezing the model layer distribution
  REAL(wp):: top_height    = 1._wp  ! Height of model top

  ! b) Parameters for SLEVE definition
  REAL(wp):: decay_scale_1 = 1._wp  ! Decay scale for large-scale topography component
  REAL(wp):: decay_scale_2 = 1._wp  ! Decay scale for small-scale topography component
  REAL(wp):: decay_exp     = 1._wp  ! Exponent for decay function
  REAL(wp):: flat_height   = 1._wp  ! Height above which the coordinate surfaces are exactly flat
                                    ! (not available in the standard SLEVE definition)

  ! c) Parameter for reading in smoothed topography
  LOGICAL :: lread_smt
 
  NAMELIST /sleve_nml/ min_lay_thckn, top_height, decay_scale_1,           &
                       decay_scale_2, decay_exp, flat_height, stretch_fac, &
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
    !0!CHARACTER(len=*), PARAMETER ::  &
    !0!  &  routine = 'mo_sleve_nml:read_sleve_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! a) Parameters determining the distribution of model layers
    !    (if not read in from a table)
    min_lay_thckn   = 50._wp      ! Layer thickness of lowermost layer
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
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('sleve_nml')
      READ(funit,NML=sleve_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('sleve_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, sleve_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    config_min_lay_thckn = min_lay_thckn
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
