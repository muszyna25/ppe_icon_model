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
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: max_char_length, max_dom, &
    &                              MODE_IFSANA, MODE_DWDANA
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_initicon_config,    ONLY: &
    & config_init_mode          => init_mode,    &
    & config_nlev_in            => nlev_in,      &
    & config_nlevsoil_in        => nlevsoil_in,  &
    & config_zpbl1              => zpbl1,        &
    & config_zpbl2              => zpbl2,        &
    & config_l_hice_in          => l_hice_in,    &
    & config_l_sst_in           => l_sst_in,     &
    & config_ifs2icon_filename  => ifs2icon_filename, &
    & config_dwdfg_filename     => dwdfg_filename, &
    & config_dwdinc_filename    => dwdinc_filename, &
    & config_l_coarse2fine_mode => l_coarse2fine_mode

  IMPLICIT NONE

  PUBLIC :: read_initicon_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the prep_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: init_mode     ! initialization mode
  INTEGER  :: nlev_in       ! number of model levels of input data
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: l_hice_in     ! Logical switch, if sea-ice thickness field is provided as input
  LOGICAL  :: l_sst_in      ! logical switch, if sea surface temperature is provided as input
  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain

  ! IFS2ICON input filename, may contain keywords, by default
  ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: ifs2icon_filename

  ! DWD-FG input filename, may contain keywords, by default
  ! dwdfg_filename = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdfg_filename

  ! DWD-inc input filename, may contain keywords, by default
  ! dwdinc_filename = "<path>dwdinc_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdinc_filename

  NAMELIST /initicon_nml/ init_mode, nlev_in, zpbl1, zpbl2, l_coarse2fine_mode, &
                          nlevsoil_in, l_hice_in, l_sst_in, ifs2icon_filename, &
                          dwdfg_filename, dwdinc_filename
  
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
  INTEGER :: z_go_init(2)   ! for consistency check

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_initicon_nml: read_initicon_namelist'

  !------------------------------------------------------------
  ! 2.0 set up the default values for initicon
  !------------------------------------------------------------
  !
  !
  init_mode   = MODE_IFSANA ! Start from IFS analysis
  nlev_in     = 91          ! number of levels of input data
  nlevsoil_in = 4           ! number of soil levels of input data
  zpbl1       = 500._wp     ! AGL heights used for computing vertical 
  zpbl2       = 1000._wp    ! gradients
  l_hice_in   = .FALSE.     ! true: sea-ice thickness field provided as input
  l_sst_in    = .FALSE.     ! true: sea surface temperature field provided as input
  ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdfg_filename    = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdinc_filename   = "<path>dwdinc_R<nroot>B<jlev>_DOM<idom>.nc"
  l_coarse2fine_mode(:) = .FALSE. ! true: apply corrections for coarse-to-fine-mesh interpolation

  !------------------------------------------------------------
  ! 3.0 Read the initicon namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL open_nml(TRIM(filename))
  CALL position_nml ('initicon_nml', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, initicon_nml)
  END SELECT
  CALL close_nml

  !------------------------------------------------------------
  ! 4.0 Fill the configuration state
  !------------------------------------------------------------

  config_init_mode         = init_mode
  config_nlev_in           = nlev_in
  config_nlevsoil_in       = nlevsoil_in
  config_zpbl1             = zpbl1
  config_zpbl2             = zpbl2
  config_l_hice_in         = l_hice_in
  config_l_sst_in          = l_sst_in
  config_ifs2icon_filename = ifs2icon_filename
  config_dwdfg_filename    = dwdfg_filename
  config_dwdinc_filename   = dwdinc_filename
  config_l_coarse2fine_mode= l_coarse2fine_mode

  !------------------------------------------------------------
  ! 5.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  z_go_init = (/MODE_IFSANA,MODE_DWDANA/)
  IF (ALL(z_go_init /= init_mode)) THEN
    CALL finish( TRIM(routine),                         &
      &  'Invalid initialization mode. Must be init_mode=1 or 2')
  ENDIF

  ! write the contents of the namelist to an ASCII file

  IF(my_process_is_stdio()) WRITE(nnml_output,nml=initicon_nml)

END SUBROUTINE read_initicon_namelist

END MODULE mo_initicon_nml
