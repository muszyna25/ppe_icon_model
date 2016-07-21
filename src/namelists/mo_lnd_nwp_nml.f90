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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_lnd_nwp_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_dom, SSTICE_ANA
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  USE mo_lnd_nwp_config,      ONLY: config_nlev_snow   => nlev_snow     , &
    &                               config_ntiles      => ntiles_lnd    , &
    &                               config_frlnd_thrhld => frlnd_thrhld , &
    &                        config_frlndtile_thrhld => frlndtile_thrhld, &
    &                             config_frlake_thrhld => frlake_thrhld , &
    &                               config_frsea_thrhld => frsea_thrhld , &
    &                               config_lseaice     => lseaice       , &
    &                               config_llake       => llake         , &
    &                               config_lmelt       => lmelt         , &
    &                               config_lmelt_var   => lmelt_var     , &
    &                               config_lmulti_snow => lmulti_snow   , &
    &                            config_l2lay_rho_snow => l2lay_rho_snow, &
    &                          config_max_toplaydepth => max_toplaydepth, &
    &                            config_idiag_snowfrac => idiag_snowfrac, &
    &                                 config_cwimax_ml => cwimax_ml     , &
    &                               config_itype_trvg  => itype_trvg    , &
    &                               config_itype_evsl  => itype_evsl    , &
    &                              config_itype_lndtbl => itype_lndtbl  , &
    &                               config_itype_root  => itype_root    , &
    &                               config_lstomata    => lstomata      , &
    &                               config_l2tls       => l2tls         , &
    &                            config_itype_heatcond => itype_heatcond, &
    &                    config_itype_interception => itype_interception, &
    &                            config_itype_hydbound => itype_hydbound, &
    &                            config_lana_rho_snow  => lana_rho_snow , &
    &                            config_lsnowtile      => lsnowtile     , &
    &                            config_sstice_mode  => sstice_mode     , &
    &                           config_sst_td_filename => sst_td_filename,&
    &                            config_ci_td_filename => ci_td_filename

  IMPLICIT NONE

  PRIVATE
! Variable to set the sst and seaice fraction mode      
  INTEGER ::  sstice_mode
!> Action Variables for physical schemes
! --------------------------------------
  INTEGER ::  nlev_snow         !< number of snow layers
  INTEGER ::  ntiles            !< number of static tiles
  REAL(wp)::  frlnd_thrhld      !< fraction threshold for creating a land grid point
  REAL(wp)::  frlndtile_thrhld  !< fraction threshold for retaining the respective 
                                !< tile for a grid point
  REAL(wp)::  frlake_thrhld     !< fraction threshold for creating a lake grid point
  REAL(wp)::  frsea_thrhld      !< fraction threshold for creating a sea grid point
  REAL(wp)::  max_toplaydepth   !< maximum depth of uppermost snow layer for multi-layer snow scheme
  INTEGER ::  itype_trvg        !< type of vegetation transpiration parameterization
  INTEGER ::  itype_evsl        !< type of parameterization of bare soil evaporation
  INTEGER ::  itype_lndtbl      !< choice of table for associating surface parameters to land-cover classes
  INTEGER ::  itype_root        !< type of root density distribution
  INTEGER ::  itype_heatcond    !< type of soil heat conductivity
  INTEGER ::  itype_interception!< type of plant interception
  REAL(wp)::  cwimax_ml         !< scaling parameter for maximum interception storage
  INTEGER ::  itype_hydbound    !< type of hydraulic lower boundary condition
  INTEGER ::  idiag_snowfrac    !< method for diagnosis of snow-cover fraction       

  CHARACTER(LEN=filename_max) :: sst_td_filename, ci_td_filename


  LOGICAL ::       &
       lseaice,    & !> forecast with sea ice model
       llake,      & !! forecast with lake model FLake
       lmelt     , & !! soil model with melting process
       lmelt_var , & !! freezing temperature dependent on water content
       lmulti_snow,& !! run the multi-layer snow model
       l2lay_rho_snow, & ! use two-layer snow density for single-layer snow scheme
       lstomata   , & ! map of minimum stomata resistance
       l2tls      , & ! forecast with 2-TL integration scheme
       lana_rho_snow, &  ! if .TRUE., take rho_snow-values from analysis file 
       lsnowtile      ! if .TRUE., snow is considered as a separate tile
!--------------------------------------------------------------------
! nwp forcing (right hand side)
!--------------------------------------------------------------------

  NAMELIST/lnd_nml/ nlev_snow, ntiles                         , &
    &               frlnd_thrhld, lseaice, llake, lmelt       , &
    &               frlndtile_thrhld, frlake_thrhld           , &
    &               frsea_thrhld, lmelt_var, lmulti_snow      , & 
    &               itype_trvg, idiag_snowfrac, max_toplaydepth, & 
    &               itype_evsl                                , & 
    &               itype_lndtbl                              , & 
    &               itype_root                                , & 
    &               itype_heatcond                            , & 
    &               itype_interception                        , & 
    &               itype_hydbound                            , & 
    &               lstomata                                  , & 
    &               l2tls                                     , & 
    &               lana_rho_snow, l2lay_rho_snow             , & 
    &               lsnowtile                                 , &
    &               sstice_mode                               , &
    &               sst_td_filename                           , &
    &               ci_td_filename, cwimax_ml
   
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
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_lnd_nwp_nml:read_nwp_lnd_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    sstice_mode  = SSTICE_ANA  ! forecast mode, sst and sea ice fraction is read from 
                               ! the analysis, sst ist kept constant, sea ice fraction
                               ! is modified by the sea ice model
                               ! default names for the time dependent SST and CI ext param files
                               ! if sstice=SSTICE_CLIM, <year> is substituted by "CLIM"
    sst_td_filename = "<path>SST_<year>_<month>_<gridfile>"
    ci_td_filename = "<path>CI_<year>_<month>_<gridfile>"


    nlev_snow      = 2       ! 2 = default value for number of snow layers
    ntiles         = 1       ! 1 = default value for number of static surface types
    frlnd_thrhld   = 0.05_wp ! fraction threshold for creating a land grid point

    frlake_thrhld  = 0.05_wp ! fraction threshold for creating a lake grid point

    frsea_thrhld   = 0.05_wp ! fraction threshold for creating a sea grid point
    frlndtile_thrhld = 0.05_wp ! fraction threshold for retaining the respective 
                             ! tile for a grid point
    lmelt          = .TRUE.  ! soil model with melting process
    lmelt_var      = .TRUE.  ! freezing temperature dependent on water content
    lmulti_snow    = .FALSE. ! .TRUE. = run the multi-layer snow model, .FALSE. = use single-layer scheme
    l2lay_rho_snow = .FALSE. ! use two-layer snow density for single-layer snow model
    max_toplaydepth = 0.25_wp ! maximum depth of uppermost snow layer for multi-layer snow scheme (25 cm)
                              ! (also used for simplified two-layer snow density scheme)
    lsnowtile      = .FALSE. ! if .TRUE., snow is considered as a separate tile
    idiag_snowfrac = 1       ! 1: old method based on SWE, 2: more advanced experimental method
    !
    itype_trvg     = 2       ! type of vegetation transpiration parameterization
                             ! Note that this is currently the only available option!
    itype_evsl     = 2       ! type of parameterization of bare soil evaporation
    itype_lndtbl   = 3       ! choice of table for associating surface parameters to land-cover classes
    itype_root     = 2       ! type of root density distribution
                             ! 1: constant
                             ! 2: exponential
    itype_heatcond = 2       ! type of soil heat conductivity
    itype_interception = 1   ! type of plant interception
    cwimax_ml      = 1.e-6_wp ! scaling parameter for maximum interception storage. Almost turned off by default;
                              ! the recommended value to activate interception storage is 5.e-4
    itype_hydbound = 1       ! type of hydraulic lower boundary condition
    lstomata       =.TRUE.   ! map of minimum stomata resistance
    l2tls          =.TRUE.   ! forecast with 2-TL integration scheme
    lana_rho_snow  =.TRUE.   ! if .TRUE., take rho_snow-values from analysis file 


    lseaice    = .TRUE.      ! .TRUE.: sea-ice model is used
    llake      = .TRUE.      ! .TRUE.: lake model is used
    


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('lnd_nml')
      READ(funit,NML=lnd_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('lnd_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, lnd_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, lnd_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, lnd_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml


    !----------------------------------------------------
    ! 4. Sanity check (if necessary)
    !----------------------------------------------------

    !Multi-layer snow model
    !
    IF ( (lmulti_snow) .AND. (nlev_snow <= 1) ) THEN
      CALL finish( TRIM(routine),                                   &
        &  'nlev_snow must be >1 when running the multi-layer snow model')
    ENDIF

    IF ( lmulti_snow .AND. l2lay_rho_snow ) THEN
      CALL finish( TRIM(routine), 'multi-layer snow model cannot be combined with l2lay_rho_snow option')
    ENDIF

    ! For simplicity, in order to avoid further case discriminations
    IF (l2lay_rho_snow) nlev_snow = 2

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg = 1,max_dom
      config_nlev_snow   = nlev_snow
      config_ntiles      = ntiles
      config_frlnd_thrhld = frlnd_thrhld
      config_frlndtile_thrhld = frlndtile_thrhld
      config_frlake_thrhld = frlake_thrhld
      config_frsea_thrhld = frsea_thrhld
      config_lseaice     = lseaice
      config_llake       = llake
      config_lmelt       = lmelt
      config_lmelt_var   = lmelt_var
      config_lmulti_snow = lmulti_snow
      config_max_toplaydepth = max_toplaydepth
      config_idiag_snowfrac = idiag_snowfrac
      config_itype_trvg  = itype_trvg
      config_itype_evsl  = itype_evsl
      config_itype_lndtbl= itype_lndtbl
      config_itype_root  = itype_root
      config_lstomata    = lstomata
      config_l2tls       = l2tls
      config_itype_heatcond = itype_heatcond
      config_itype_interception = itype_interception
      config_cwimax_ml   = cwimax_ml
      config_itype_hydbound = itype_hydbound
      config_lana_rho_snow  = lana_rho_snow
      config_l2lay_rho_snow = l2lay_rho_snow
      config_lsnowtile   = lsnowtile
      config_sstice_mode   = sstice_mode
      config_sst_td_filename = sst_td_filename
      config_ci_td_filename = ci_td_filename

    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=lnd_nml)                    
      CALL store_and_close_namelist(funit, 'lnd_nml') 
    ENDIF


    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=lnd_nml)

  END SUBROUTINE read_nwp_lnd_namelist


END MODULE mo_lnd_nwp_nml

