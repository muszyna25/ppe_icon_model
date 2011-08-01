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

  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_lnd_nwp_config,      ONLY: config_nztlev      => nztlev     , &
    &                               config_nlev_snow   => nlev_snow  , &
    &                               config_nsfc_subs   => nsfc_subs  , &
    &                               config_lseaice     => lseaice    , &
    &                               config_llake       => llake      , &
    &                               config_lmelt       => lmelt      , &
    &                               config_lmelt_var   => lmelt_var  , &
    &                               config_lmulti_snow => lmulti_snow

  IMPLICIT NONE

  PRIVATE

!> Action Variables for physical schemes
! --------------------------------------
  INTEGER ::  nztlev            !< time integration scheme
  INTEGER ::  nlev_snow         !< number of snow layers
  INTEGER ::  nsfc_subs         !< number of TILES
  INTEGER ::  itype_gscp        !< type of grid-scale precipitation physics
  INTEGER ::  itype_trvg        !< type of vegetation transpiration parameterization
  INTEGER ::  itype_evsl        !< type of parameterization of bare soil evaporation
  INTEGER ::  itype_tran        !< type of surface to atmospher transfer
  INTEGER ::  itype_root        !< type of root density distribution
  INTEGER ::  itype_heatcond    !< type of soil heat conductivity
  INTEGER ::  itype_hydbound    !< type of hydraulic lower boundary condition
  INTEGER ::  itype_subs        !< type of subscale surface treatment =1 MOSAIC, =2 TILE       




  LOGICAL ::       &
       lseaice,    & !> forecast with sea ice model
       llake,      & !! forecst with lake model FLake
       lmelt     , & !! soil model with melting process
       lmelt_var , & !! freezing temperature dependent on water content
       lmulti_snow,& !! run the multi-layer snow model
       lstomata   , & ! map of minimum stomata resistance
       l2tls      , & ! forecast with 2-TL integration scheme
       lana_rho_snow                         ! if .TRUE., take rho_snow-values from analysis file 
!--------------------------------------------------------------------
! nwp forcing (right hand side)
!--------------------------------------------------------------------

  NAMELIST/lnd_nml/ nztlev, nlev_snow, nsfc_subs, lseaice, llake, &
    &               lmelt, lmelt_var, lmulti_snow, itype_gscp   , & 
    &               itype_trvg                                  , & 
    &               itype_evsl                                  , & 
    &               itype_tran                                  , & 
    &               itype_root                                  , & 
    &               itype_heatcond                              , & 
    &               itype_hydbound                              , & 
    &               lstomata                                    , & 
    &               l2tls                                       , & 
    &               lana_rho_snow                               , & 
    &               itype_subs            
   
  PUBLIC :: read_nwp_lnd_namelist

 CONTAINS


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

    !0!CHARACTER(len=*), PARAMETER :: routine = 'mo_lnd_nwp_nml:read_nwp_lnd_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    nztlev          = 2       !> 2 = default value for time integration scheme
    nlev_snow       = 1       !> 0 = default value for number of snow layers
    nsfc_subs       = 12      !> 1 = default value for number of TILES
    nztlev          = 2       ! time step scheme 2,3
    lmelt       = .TRUE.      ! soil model with melting process
    lmelt_var   = .TRUE.      ! freezing temperature dependent on water content
    lmulti_snow = .FALSE.     ! run the multi-layer snow model
    !
    itype_gscp=3              ! type of grid-scale precipitation physics
    itype_trvg=2              ! type of vegetation transpiration parameterization
    itype_evsl=2              ! type of parameterization of bare soil evaporation
    itype_tran=2              ! type of surface to atmospher transfer
    itype_root=1              ! type of root density distribution
    itype_heatcond=1          ! type of soil heat conductivity
    itype_hydbound=1          ! type of hydraulic lower boundary condition
    lstomata=.true.           ! map of minimum stomata resistance
    l2tls  =.true.            ! forecast with 2-TL integration scheme
    lana_rho_snow=.false.     ! if .TRUE., take rho_snow-values from analysis file 
    itype_subs=2              ! type of subscale surface treatment =1 MOSAIC, =2 TILE       


    !> KF  current settings to get NWP turbulence running
    lseaice    = .FALSE.
    llake      = .FALSE.
    lmulti_snow= .FALSE.
    


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('lnd_nml')
      READ(funit,NML=lnd_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
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
      config_nztlev      = nztlev
      config_nlev_snow   = nlev_snow
      config_nsfc_subs   = nsfc_subs
      config_lseaice     = lseaice
      config_llake       = llake
      config_lmelt       = lmelt
      config_lmelt_var   = lmelt_var
      config_lmulti_snow = lmulti_snow
    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=lnd_nml)                    
      CALL store_and_close_namelist(funit, 'lnd_nml') 
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=lnd_nml)

  END SUBROUTINE read_nwp_lnd_namelist


END MODULE mo_lnd_nwp_nml

