!>
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
MODULE mo_initicon_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_impl_constants,     ONLY: max_char_length, max_dom, vname_len,      &
    &                              max_var_ml, MODE_IFSANA, MODE_DWDANA,     &
    &                              MODE_DWDANA_INC, MODE_IAU, MODE_IAU_OLD,  &
    &                              MODE_COMBINED,                            &
    &                              MODE_COSMODE, MODE_ICONVREMAP
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_initicon_config,    ONLY: &
    & config_init_mode           => init_mode,           &
    & config_nlevsoil_in         => nlevsoil_in,         &
    & config_zpbl1               => zpbl1,               &
    & config_zpbl2               => zpbl2,               &
    & config_l_sst_in            => l_sst_in,            &
    & config_lread_ana           => lread_ana,           &
    & config_lconsistency_checks => lconsistency_checks, &
    & config_ifs2icon_filename   => ifs2icon_filename,   &
    & config_dwdfg_filename      => dwdfg_filename,      &
    & config_dwdana_filename     => dwdana_filename,     &
    & config_l_coarse2fine_mode  => l_coarse2fine_mode,  &
    & config_lp2cintp_incr       => lp2cintp_incr,       &
    & config_lp2cintp_sfcana     => lp2cintp_sfcana,     &
    & config_ltile_coldstart     => ltile_coldstart,     &
    & config_start_time_avg_fg   => start_time_avg_fg,   &
    & config_end_time_avg_fg     => end_time_avg_fg,     &
    & config_interval_avg_fg     => interval_avg_fg,     &
    & config_filetype            => filetype,            &
    & config_dt_iau              => dt_iau,              &
    & config_timeshift           => timeshift,           &
    & config_type_iau_wgt        => type_iau_wgt,        &
    & config_ana_varlist         => ana_varlist,         &
    & config_rho_incr_filter_wgt => rho_incr_filter_wgt, &
    & config_wgtfac_geobal       => wgtfac_geobal,       &
    & config_ana_varnames_map_file => ana_varnames_map_file, &
    & config_latbc_varnames_map_file => latbc_varnames_map_file

  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings


  IMPLICIT NONE

  PUBLIC :: read_initicon_namelist

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the init_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: init_mode     ! initialization mode
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: l_sst_in      ! logical switch, if sea surface temperature is provided as input
  LOGICAL  :: lread_ana     ! If .TRUE., read analysis fields are read from analysis file
                            ! dwdana_filename. If .FALSE., ICON is soleyly started 
                            ! from first guess fields.   
  LOGICAL  :: lconsistency_checks    ! check validity of input fields (FG and ANA)

  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain
  LOGICAL  :: lp2cintp_incr(max_dom) ! If true, perform parent-to-child interpolation of atmospheric data
                                     ! assimilation increments
  LOGICAL  :: lp2cintp_sfcana(max_dom) ! If true, perform parent-to-child interpolation of
                                       ! surface analysis data
  LOGICAL  :: ltile_coldstart  ! If true, initialize tile-based surface fields from first guess without tiles

  ! Variables controlling computation of temporally averaged first guess fields for DA
  ! The calculation is switched on by setting end_time > start_time
  REAL(wp) :: start_time_avg_fg   ! start time [s]
  REAL(wp) :: end_time_avg_fg     ! end time [s]
  REAL(wp) :: interval_avg_fg     ! averaging interval [s]

  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

  REAL(wp) :: dt_iau        ! Time interval during which incremental analysis update (IAU) is performed [s]. 
                            ! Only required for init_mode=MODE_IAU, MODE_IAU_OLD, MODE_DWDANA_INC
  REAL(wp) :: dt_shift      ! Allows IAU runs to start earlier than the nominal simulation start date without showing up in the output metadata

  INTEGER  :: type_iau_wgt  ! Type of weighting function for IAU.
                            ! 1: Top-hat
                            ! 2: SIN2
                            ! Only required for init_mode=MODE_IAU, MODE_IAU_OLD, MODE_DWDANA_INC
  REAL(wp) :: rho_incr_filter_wgt  ! Vertical filtering weight for density increments 
                                   ! Only applicable for init_mode=MODE_IAU, MODE_IAU_OLD, MODE_DWDANA_INC
  REAL(wp) :: wgtfac_geobal  ! Weighting factor for artificial geostrophic balancing of meridional gradients
                             ! of pressure increments in the tropical stratosphere

  CHARACTER(LEN=vname_len) :: ana_varlist(max_var_ml) ! list of mandatory analysis fields. 
                                                      ! This list can include a subset or the 
                                                      ! entire set of default analysis fields.

  ! IFS2ICON input filename, may contain keywords, by default
  ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: ifs2icon_filename

  ! DWD-FG input filename, may contain keywords, by default
  ! dwdfg_filename = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdfg_filename

  ! DWD-ANA input filename, may contain keywords, by default
  ! dwdana_filename = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdana_filename

  ! analysis file: dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: ana_varnames_map_file      

  ! analysis file: dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names used for lateral boundary nudging.
  CHARACTER(LEN=filename_max) :: latbc_varnames_map_file  

  NAMELIST /initicon_nml/ init_mode, zpbl1, zpbl2, l_coarse2fine_mode,      &
                          nlevsoil_in, l_sst_in, lread_ana,                 &
                          lconsistency_checks, rho_incr_filter_wgt,         &
                          ifs2icon_filename, dwdfg_filename,                &
                          dwdana_filename, filetype, dt_iau, dt_shift,      &
                          type_iau_wgt, ana_varlist, ana_varnames_map_file, &
                          lp2cintp_incr, lp2cintp_sfcana, wgtfac_geobal,    &
                          latbc_varnames_map_file, start_time_avg_fg,       &
                          end_time_avg_fg, interval_avg_fg, ltile_coldstart
                          
CONTAINS

!-------------------------------------------------------------------------
!
!
 !>
 !!  Initialization of the initicon coordinate namelist
 !!
 !!
 !! @par Revision History
 !!  Initial version by Guenther Zaengl (2011-07-11)

 SUBROUTINE read_initicon_namelist( filename )
    
  CHARACTER(LEN=*), INTENT(IN) :: filename

  !local variable
  INTEGER :: i_status
  INTEGER :: z_go_init(8)   ! for consistency check
  INTEGER :: iunit

  CHARACTER(len=*), PARAMETER ::  &
    &  routine = 'mo_initicon_nml: read_initicon_namelist'

  !------------------------------------------------------------
  ! 2.0 set up the default values for initicon
  !------------------------------------------------------------
  !
  !
  init_mode   = MODE_IFSANA    ! Start from IFS analysis
  nlevsoil_in = 4              ! number of soil levels of input data
  zpbl1       = 500._wp        ! AGL heights used for computing vertical 
  zpbl2       = 1000._wp       ! gradients
  l_sst_in    = .TRUE.         ! true: sea surface temperature field provided as input
  lread_ana   = .TRUE.         ! true: read analysis fields from file dwdana_filename
                               ! false: start ICON from first guess file (no analysis)
  lconsistency_checks = .TRUE. ! check validity of input fields  
  filetype    = -1             ! "-1": undefined
  dt_iau      = 10800._wp      ! 3-hour interval for IAU
  dt_shift    = 0._wp          ! do not shift actual simulation start backward
  rho_incr_filter_wgt = 0._wp  ! density increment filtering turned off
  wgtfac_geobal       = 0._wp  ! geostrophic balancing of pressure increments in the tropical stratosphere turned off
  type_iau_wgt= 1              ! Top-hat weighting function
  ana_varlist = ''             ! list of mandatory analysis fields. This list can include a subset 
                               ! or the entire set of default analysis fields. If any of these fields
                               ! is missing in the analysis file, the model aborts. On default 
                               ! this list is empty, meaning that fields which are missing in the 
                               ! analysis file (when compared to the default set), are simply 
                               ! taken from the first guess.
  ana_varnames_map_file = " "
  latbc_varnames_map_file = " "
  ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdfg_filename    = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdana_filename   = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  l_coarse2fine_mode(:) = .FALSE. ! true: apply corrections for coarse-to-fine-mesh interpolation
  lp2cintp_incr(:)      = .FALSE. ! true: perform parent-to-child interpolation of atmospheric data assimilation increments
  lp2cintp_sfcana(:)    = .FALSE. ! true: perform parent-to-child interpolation of surface analysis data
  ltile_coldstart       = .FALSE. ! true: initialize tile-based surface fields from first guess without tiles

  start_time_avg_fg = 0._wp
  end_time_avg_fg   = 0._wp
  interval_avg_fg   = 0._wp

  !------------------------------------------------------------
  ! 3.0 Read the initicon namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL open_nml(TRIM(filename))
  CALL position_nml ('initicon_nml', status=i_status)
  IF (my_process_is_stdio()) THEN
    iunit = temp_defaults()
    WRITE(iunit, initicon_nml)  ! write defaults to temporary text file
  END IF
  SELECT CASE (i_status)
  CASE (positioned)
    READ (nnml, initicon_nml)                                      ! overwrite default settings
    IF (my_process_is_stdio()) THEN
      iunit = temp_settings()
      WRITE(iunit, initicon_nml)  ! write settings to temporary text file
    END IF
  END SELECT
  CALL close_nml


  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  z_go_init = (/MODE_IFSANA,MODE_DWDANA,MODE_DWDANA_INC,MODE_IAU,MODE_IAU_OLD,MODE_COMBINED,MODE_COSMODE,MODE_ICONVREMAP/)
  IF (ALL(z_go_init /= init_mode)) THEN
    CALL finish( TRIM(routine),                         &
      &  'Invalid initialization mode. init_mode must be between 1 and 8')
  ENDIF

  ! Check whether a NetCDF<=>GRIB2 Map File is needed, and if so, whether 
  ! it is provided
  IF (ANY((/MODE_DWDANA,MODE_DWDANA_INC,MODE_IAU,MODE_IAU_OLD,MODE_COMBINED/)==init_mode)) THEN
    ! NetCDF<=>GRIB2 Map File required
    IF(ana_varnames_map_file == ' ') THEN
    CALL finish( TRIM(routine),                         &
      &  'ana_varnames_map_file required, but missing.')
    ENDIF
  ENDIF 

  ! Check whther an analysis file is provided, if lread_ana=.TRUE.
  IF (lread_ana) THEN
    IF (dwdana_filename ==' ') THEN
    CALL finish( TRIM(routine),                         &
      &  'dwdana_filename required, but missing.')
    ENDIF
  ENDIF

  ! Check whether init_mode and lread_ana are consistent
  IF (ANY((/MODE_COMBINED,MODE_COSMODE,MODE_ICONVREMAP/)==init_mode) .AND. lread_ana) THEN
    lread_ana = .FALSE.
    WRITE(message_text,'(a,i2,a)') 'init_mode=', init_mode, &
      '. no analysis required => lread_ana re-set to .FALSE.'
    CALL message(TRIM(routine),message_text)
  ENDIF

  ! Setting the first entry of lp2cintp_incr / lp2cintp_sfcana to true activates parent-to-child interpolation
  ! of DA increments / surface analysis for all domains
  IF (lp2cintp_incr(1)) THEN
    lp2cintp_incr(2:max_dom) = .TRUE.
  ENDIF
  IF (lp2cintp_sfcana(1)) THEN
    lp2cintp_sfcana(2:max_dom) = .TRUE.
  ENDIF
  ! To simplify runtime flow control, set the switches for the global domain to false
  lp2cintp_incr(1)   = .FALSE.
  lp2cintp_sfcana(1) = .FALSE.

  ! Check if settings for temporally averaged first guess output make sense
  IF (end_time_avg_fg > start_time_avg_fg) THEN
    IF (interval_avg_fg < 1.e-10_wp) THEN
      WRITE(message_text,'(a,f8.2,a)') 'averaging interval for first guess output must be positive'
      CALL finish(TRIM(routine),message_text)
    ENDIF
    IF (end_time_avg_fg < start_time_avg_fg + interval_avg_fg) THEN
      WRITE(message_text,'(a,f8.2,a)') &
        'averaging period for first guess output must be larger than averaging interval'
      CALL finish(TRIM(routine),message_text)
    ENDIF
  ENDIF

  ! make sure that dt_shift is negative or 0.
  IF ( dt_shift > 0._wp ) THEN
    WRITE(message_text,'(a,f8.2,a)') 'dt_shift=', dt_shift, &
      ' not allowed. Must be NEGATIVE or 0.'
    CALL finish(TRIM(routine),message_text)
  ENDIF


  !------------------------------------------------------------
  ! 5.0 Fill the configuration state
  !------------------------------------------------------------

  config_init_mode           = init_mode
  config_nlevsoil_in         = nlevsoil_in
  config_zpbl1               = zpbl1
  config_zpbl2               = zpbl2
  config_l_sst_in            = l_sst_in
  config_lread_ana           = lread_ana
  config_lconsistency_checks = lconsistency_checks
  config_ifs2icon_filename   = ifs2icon_filename
  config_dwdfg_filename      = dwdfg_filename
  config_dwdana_filename     = dwdana_filename
  config_l_coarse2fine_mode  = l_coarse2fine_mode
  config_lp2cintp_incr       = lp2cintp_incr
  config_lp2cintp_sfcana     = lp2cintp_sfcana
  config_ltile_coldstart     = ltile_coldstart
  config_start_time_avg_fg   = start_time_avg_fg
  config_end_time_avg_fg     = end_time_avg_fg
  config_interval_avg_fg     = interval_avg_fg
  config_filetype            = filetype
  config_dt_iau              = dt_iau
  config_timeshift%dt_shift  = dt_shift
  config_type_iau_wgt        = type_iau_wgt
  config_ana_varlist         = ana_varlist
  config_ana_varnames_map_file = ana_varnames_map_file
  config_rho_incr_filter_wgt   = rho_incr_filter_wgt
  config_wgtfac_geobal       = wgtfac_geobal
  config_latbc_varnames_map_file = latbc_varnames_map_file

  ! write the contents of the namelist to an ASCII file

  IF(my_process_is_stdio()) WRITE(nnml_output,nml=initicon_nml)

END SUBROUTINE read_initicon_namelist

END MODULE mo_initicon_nml
