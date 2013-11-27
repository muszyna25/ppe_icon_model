!>
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
  USE mo_impl_constants,     ONLY: max_char_length, max_dom, vname_len,  &
    &                              max_var_ml, MODE_IFSANA, MODE_DWDANA, &
    &                              MODE_COMBINED, MODE_COSMODE
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_initicon_config,    ONLY: &
    & config_init_mode          => init_mode,         &
    & config_nlevsoil_in        => nlevsoil_in,       &
    & config_zpbl1              => zpbl1,             &
    & config_zpbl2              => zpbl2,             &
    & config_l_sst_in           => l_sst_in,          &
    & config_lread_ana          => lread_ana,         &
    & config_ifs2icon_filename  => ifs2icon_filename, &
    & config_dwdfg_filename     => dwdfg_filename,    &
    & config_dwdana_filename    => dwdana_filename,   &
    & config_l_coarse2fine_mode => l_coarse2fine_mode,&
    & config_filetype           => filetype,          &
    & config_ana_varlist        => ana_varlist,       &
    & config_ana_varnames_map_file => ana_varnames_map_file
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings


  IMPLICIT NONE

  PUBLIC :: read_initicon_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

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
  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain
  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

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


  NAMELIST /initicon_nml/ init_mode, zpbl1, zpbl2, l_coarse2fine_mode,   &
                          nlevsoil_in, l_sst_in, lread_ana,              &
                          ifs2icon_filename, dwdfg_filename,             &
                          dwdana_filename, filetype, ana_varlist,        &
                          ana_varnames_map_file
  
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
  INTEGER :: i_status, istat
  INTEGER :: z_go_init(4)   ! for consistency check

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_initicon_nml: read_initicon_namelist'

  !------------------------------------------------------------
  ! 2.0 set up the default values for initicon
  !------------------------------------------------------------
  !
  !
  init_mode   = MODE_IFSANA ! Start from IFS analysis
  nlevsoil_in = 4           ! number of soil levels of input data
  zpbl1       = 500._wp     ! AGL heights used for computing vertical 
  zpbl2       = 1000._wp    ! gradients
  l_sst_in    = .TRUE.      ! true: sea surface temperature field provided as input
  lread_ana   = .TRUE.      ! true: read analysis fields from file dwdana_filename
                            ! false: start ICON from first guess file (no analysis) 
  filetype    = -1          ! "-1": undefined
  ana_varlist = ''          ! list of mandatory analysis fields. This list can include a subset 
                            ! or the entire set of default analysis fields. If any of these fields
                            ! is missing in the analysis file, the model aborts. On default 
                            ! this list is empty, meaning that fields which are missing in the 
                            ! analysis file (when compared to the default set), are simply 
                            ! taken from the first guess.
  ana_varnames_map_file = " "
  ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdfg_filename    = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdana_filename   = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  l_coarse2fine_mode(:) = .FALSE. ! true: apply corrections for coarse-to-fine-mesh interpolation

  !------------------------------------------------------------
  ! 3.0 Read the initicon namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL open_nml(TRIM(filename))
  CALL position_nml ('initicon_nml', status=i_status)
  IF (my_process_is_stdio()) WRITE(temp_defaults(), initicon_nml)  ! write defaults to temporary text file
  SELECT CASE (i_status)
  CASE (positioned)
    READ (nnml, initicon_nml)                                      ! overwrite default settings
    IF (my_process_is_stdio()) WRITE(temp_settings(), initicon_nml)  ! write settings to temporary text file
  END SELECT
  CALL close_nml


  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  z_go_init = (/MODE_IFSANA,MODE_DWDANA,MODE_COMBINED,MODE_COSMODE/)
  IF (ALL(z_go_init /= init_mode)) THEN
    CALL finish( TRIM(routine),                         &
      &  'Invalid initialization mode. Must be init_mode=1, 2, 3, or 4')
  ENDIF

  ! Check whether a NetCDF<=>GRIB2 Map File is needed, and if so, whether 
  ! it is provided
  IF (ANY((/MODE_DWDANA,MODE_COMBINED/)==init_mode)) THEN
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
  IF (ANY((/MODE_COMBINED,MODE_COSMODE/)==init_mode) .AND. lread_ana) THEN
    lread_ana = .FALSE.
    WRITE(message_text,'(a,i2,a)') 'init_mode=', init_mode, &
      '. no analysis required => lread_ana re-set to .FALSE.'
    CALL message(TRIM(routine),message_text)
  ENDIF
  


  !------------------------------------------------------------
  ! 5.0 Fill the configuration state
  !------------------------------------------------------------

  config_init_mode          = init_mode
  config_nlevsoil_in        = nlevsoil_in
  config_zpbl1              = zpbl1
  config_zpbl2              = zpbl2
  config_l_sst_in           = l_sst_in
  config_lread_ana          = lread_ana
  config_ifs2icon_filename  = ifs2icon_filename
  config_dwdfg_filename     = dwdfg_filename
  config_dwdana_filename    = dwdana_filename
  config_l_coarse2fine_mode = l_coarse2fine_mode
  config_filetype           = filetype
  config_ana_varlist        = ana_varlist
  config_ana_varnames_map_file = ana_varnames_map_file


  ! write the contents of the namelist to an ASCII file

  IF(my_process_is_stdio()) WRITE(nnml_output,nml=initicon_nml)

END SUBROUTINE read_initicon_namelist

END MODULE mo_initicon_nml
