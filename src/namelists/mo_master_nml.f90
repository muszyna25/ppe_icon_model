!>
!! @brief Master namelist.
!!        
!! @par Revision History
!! Created by Rene Redler (2011-03-22)
!!
!! @par Copyright
!! 2010-2011 by DWD and MPI-M
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
MODULE mo_master_nml

  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_exception,      ONLY: warning, message_text, finish
  USE mo_io_units,       ONLY: filename_max, nnml
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED

  IMPLICIT NONE

  PRIVATE
  
  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_master_namelist, lrestart, nml_debug_coupler_level, no_of_models
  PUBLIC :: t_master_nml, master_nml_array


  ! Component models
  !--------------------------------------------------------------
  ! TYPE definitions
  !> Holds a list of initegers
  TYPE t_master_nml
    CHARACTER(len=132) :: model_name
    CHARACTER(len=filename_max) :: model_namelist_filename
    CHARACTER(len=filename_max) :: model_restart_info_filename
    INTEGER :: model_type
    INTEGER :: model_min_rank
    INTEGER :: model_max_rank
    INTEGER :: model_inc_rank
  END TYPE t_master_nml

  INTEGER, PARAMETER :: max_no_of_models=10
  INTEGER :: no_of_models
  TYPE(t_master_nml) :: master_nml_array(max_no_of_models)
  
    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------
    LOGICAL :: lrestart, nml_debug_coupler
    INTEGER :: nml_debug_coupler_level

CONTAINS
  !>
  !! Initialization of variables that contain general information
  !! about the coupled model run. The configuration is read from
  !! namelist 'master_nml'.
  !!
  !! @par Revision History
  !!
  INTEGER FUNCTION read_master_namelist(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    !
    ! Local variables
    !
    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------
    CHARACTER(len=132) :: model_name
    CHARACTER(len=filename_max) :: model_namelist_filename
    CHARACTER(len=filename_max) :: model_restart_info_filename
    INTEGER :: model_type
    INTEGER :: model_min_rank
    INTEGER :: model_max_rank
    INTEGER :: model_inc_rank

    INTEGER :: debug_coupler_level

    NAMELIST /master_model_nml/    &
      model_name,                  &
      model_namelist_filename,     &
      model_restart_info_filename, &
      model_type,                  &
      model_min_rank,              &
      model_max_rank,              &
      model_inc_rank               
    NAMELIST /master_nml/ lrestart, debug_coupler_level

    INTEGER :: istat
    LOGICAL :: rewnd
    CHARACTER(len=*), PARAMETER :: routine = 'mo_master_nml:read_master_namelist'


    !------------------------------------------------------------------
    ! Read  master_nml (done so far by all MPI processes)
    !------------------------------------------------------------------
    lrestart      = .FALSE.
    debug_coupler_level = 0
    OPEN( nnml, FILE=TRIM(namelist_filename), IOSTAT=istat, &
        & STATUS='old', ACTION='read', DELIM='apostrophe')
    IF (istat/=0) THEN
      CALL warning(namelist_filename,"not found")
      read_master_namelist=-1
      RETURN
    ENDIF
    CALL position_nml('master_nml',STATUS=istat)
    IF (istat==POSITIONED) THEN
      READ (nnml, master_nml)
    ENDIF        
    nml_debug_coupler_level = debug_coupler_level
    !------------------------------------------------------------------
    ! Read  master_model_nml (done so far by all MPI processes)
    !------------------------------------------------------------------
    no_of_models = 0
    rewnd = .true.
    DO
      CALL position_nml('master_model_nml', lrewind=rewnd, status=istat)
      IF ( istat /= POSITIONED ) EXIT
      IF (no_of_models >= max_no_of_models) THEN
        CALL finish(routine, 'no_of_models >= max_no_of_models')
      ENDIF
      rewnd=.false.
      
      ! default values
      model_name=''
      model_namelist_filename=''
      model_restart_info_filename=''
      model_type=-1
      model_min_rank=0
      model_max_rank=-1 
      model_inc_rank=1
      
      READ  (nnml, master_model_nml)

      no_of_models=no_of_models+1
      master_nml_array(no_of_models)%model_name              = model_name
      master_nml_array(no_of_models)%model_namelist_filename = model_namelist_filename
      master_nml_array(no_of_models)%model_restart_info_filename=&
        & model_restart_info_filename
      master_nml_array(no_of_models)%model_type              = model_type
      master_nml_array(no_of_models)%model_min_rank          = model_min_rank
      master_nml_array(no_of_models)%model_max_rank          = model_max_rank
      master_nml_array(no_of_models)%model_inc_rank          = model_inc_rank

    ENDDO
      
    CLOSE (nnml, IOSTAT=istat)
   !------------------------------------------------------------------
   read_master_namelist=SUCCESS

  END FUNCTION read_master_namelist

END MODULE mo_master_nml
