!>
!!  Namelist for surface physics
!!
!!  these Subroutines are called by control model and construct the
!!  surface scheme composition
!!
!! @author <Kristina Froehlich, DWD>
!!
!!
!! @par Revision History
!! First implementation by Kristina Froehlich, DWD (2010-06-20>)
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
MODULE mo_lnd_nwp_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_master_nml,          ONLY: lrestart
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_nml,          ONLY: lrestart
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_lnd_nwp_config,      ONLY: config_nlev_soil   => nlev_soil  , &
    &                               config_nztlev      => nztlev     , &
    &                               config_nlev_snow   => nlev_snow  , &
    &                               config_nsfc_subs   => nsfc_subs  , &
    &                               config_lseaice     => lseaice    , &
    &                               config_llake       => llake      , &
    &                               config_lmelt       => lmelt      , &
    &                               config_lmelt_var   => lmelt_var  , &
    &                               config_lmulti_snow => lmulti_snow

  IMPLICIT NONE

  PRIVATE

  INTEGER ::  nlev_soil,  nztlev  !! number of soil layers, time integration scheme
  INTEGER ::  nlev_snow           !! number of snow layers
  INTEGER ::  nsfc_subs           !! number of TILES



  LOGICAL ::  lseaice     !> forecast with sea ice model
  LOGICAL ::  llake       !! forecst with lake model FLake
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model


  NAMELIST /lnd_nml/ nlev_soil, nztlev, nlev_snow, nsfc_subs, &
    &                lseaice, llake, lmulti_snow  
   
  PUBLIC :: read_nwp_lnd_namelist, setup_nwp_lnd !(latter to be removed)

 CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Setup NWP physics
  !!
  !! Read namelist for physics. Choose the physical package and subsequent
  !! parameters.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-10-06)
  !!
  SUBROUTINE setup_nwp_lnd

    INTEGER :: i_stat, funit

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------

     nlev_soil       = 7     !> 7 = default value for number of soil layers
     nztlev          = 2     !> 2 = default value for time integration scheme
     nlev_snow       = 1     !> 0 = default value for number of snow layers
     nsfc_subs       = 2     !> 1 = default value for number of TILES



  !> KF  current settings to get NWP turbulence running
     lseaice    = .FALSE.
     llake      = .FALSE.
     lmulti_snow= .FALSE.
    

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('lnd_nml')
      READ(funit,NML=lnd_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processors)
    !------------------------------------------------------------------------

    CALL position_nml ('lnd_nml', status=i_stat)
    IF (i_stat == POSITIONED) THEN
      READ (nnml, lnd_nml)
    ENDIF

  END SUBROUTINE setup_nwp_lnd


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP land physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP land physics
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
  SUBROUTINE read_nwp_lnd_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg            ! loop index

    CHARACTER(len=*), PARAMETER :: routine = 'mo_lnd_nwp_nml:read_nwp_lnd_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    nlev_soil   = 7     !> 7 = default value for number of soil layers
    nztlev      = 2     !> 2 = default value for time integration scheme
    nlev_snow   = 1     !> 0 = default value for number of snow layers
    nsfc_subs   = 2     !> 1 = default value for number of TILES

    !> KF  current settings to get NWP turbulence running
    lseaice     = .FALSE.
    llake       = .FALSE.
    lmulti_snow = .FALSE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('lnd_nml')
      READ(funit,NML=lnd_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processors)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('lnd_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, lnd_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    DO jg = 1,max_dom
      config_nlev_soil   = nlev_soil
      config_nztlev      = nztlev
      config_nlev_snow   = nlev_snow
      config_nsfc_subs   = nsfc_subs
      config_lseaice     = lseaice
      config_llake       = llake
      config_lmulti_snow = lmulti_snow
    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=lnd_nml)                    
    CALL store_and_close_namelist(funit, 'lnd_nml') 

    ! 6. write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=lnd_nml)

  END SUBROUTINE read_nwp_lnd_namelist


END MODULE mo_lnd_nwp_nml

