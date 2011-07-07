!>
!! Contains the setup of the sleve coordinate
!!
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3919)
!!   Modification by Constantin Junk (2011-03-29)
!!   - added new module mo_sleve_nml.f90 which includes
!!     setting up the SLEVE namelist.
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
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned
  USE mo_sleve_config,        ONLY: sleve_config
  USE mo_master_nml,          ONLY: lrestart
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC


  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the SLEVE coordinate
  ! ----------------------------------------------------------------------------
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
                            ! additional feature not available in the standard SLEVE definition

  NAMELIST /sleve_ctl/ min_lay_thckn, top_height, decay_scale_1,     &
                       decay_scale_2, decay_exp, flat_height, stretch_fac
  !
  !
  !

  PUBLIC:: read_sleve_namelist

CONTAINS

!-------------------------------------------------------------------------
!
!
 !>
 !!  Initialization of the SLEVE coordinate namelist
 !!
 !!
 !! @par Revision History
 !!  Initial version by Guenther Zaengl (2010-07-21)

 SUBROUTINE sleve_nml_setup

!   CHARACTER(len=max_char_length), PARAMETER :: &
!             routine = 'mo_sleve_nml/sleve_nml_setup:'


  !local variable
  INTEGER :: i_status

  !------------------------------------------------------------
  ! 2.0 set up the default values for dynamics_ctl
  !------------------------------------------------------------
  !
  !
  ! a) Parameters determining the distribution of model layers
  !    (if not read in from a table)
  min_lay_thckn   = 50._wp      ! Layer thickness of lowermost layer
  top_height      = 23500._wp   ! Height of model top
  stretch_fac     = 1._wp       ! Scaling factor for stretching/squeezing the model layer distribution

  ! b) Parameters setting up the decay function of the topographic signal
  decay_scale_1   = 4000._wp    ! Decay scale of large-scale topography component
  decay_scale_2   = 2500._wp    ! Decay scale of small-scale topography component
  decay_exp       = 1.2_wp      ! Exponent for decay function
  flat_height     = 16000._wp   ! Height above which the coordinate surfaces are flat
  !
  !
  !------------------------------------------------------------
  ! 3.0 Read the nonhydrostatic namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL position_nml ('sleve_ctl', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, sleve_ctl)
  END SELECT
  !
  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  !currently no consistency check...

  ! write the contents of the namelist to an ASCII file

  IF(p_pe == p_io) WRITE(nnml_output,nml=sleve_ctl)

END SUBROUTINE sleve_nml_setup

  !-------------------------------------------------------------------------
  !
  !
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
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_sleve_namelist
    !
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_sleve_nml: read_sleve_namelist'

    !-----------------------------------------------------------------------

    !-----------------------!
    ! 1. default settings   !
    !-----------------------!

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


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('sleve_ctl')
      READ(funit,NML=sleve_ctl)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('sleve_ctl', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, sleve_ctl)
    END SELECT


    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1,max_dom
      sleve_config(jg)%min_lay_thckn = min_lay_thckn
      sleve_config(jg)%top_height    = top_height
      sleve_config(jg)%decay_scale_1 = decay_scale_1
      sleve_config(jg)%decay_scale_2 = decay_scale_2
      sleve_config(jg)%decay_exp     = decay_exp
      sleve_config(jg)%flat_height   = flat_height
      sleve_config(jg)%stretch_fac   = stretch_fac
    ENDDO


    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=sleve_ctl)                    
    CALL store_and_close_namelist(funit, 'sleve_ctl') 


    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=sleve_ctl)


  END SUBROUTINE read_sleve_namelist


END MODULE mo_sleve_nml
