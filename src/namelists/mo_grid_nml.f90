!>
!! @brief Contains the variables to set up the grid configuration
!!        
!! @par Revision History
!!   Revision History in mo_model_domimp_setup.f90 (r3965)
!!   Modification by Constantin Junk, MPI-M (2011-04-05)
!!   - moved setup_files to grid_nml_setup and restructured
!!     reading the namelist and setting the default variables
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
MODULE mo_grid_nml

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, MAX_DOM
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_run_nml,            ONLY: lrestart
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, POSITIONED
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC
  PUBLIC :: grid_nml_setup

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  !-------------------------------------------------------------------------
  ! Namelist variables and auxiliary variables
  !-------------------------------------------------------------------------

  INTEGER :: nroot                    ! root division of initial edges
  INTEGER :: start_lev                ! coarsest bisection level
  INTEGER :: n_dom                    ! number of model domains, 1=global domain only 
  INTEGER :: n_dom_start=1 
  INTEGER :: max_childdom
  INTEGER :: parent_id(max_dom-1)  !ID of parent domain

  LOGICAL :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
  LOGICAL :: lredgrid_phys(max_dom)   ! If set to .true. is calculated on a reduced grid
  LOGICAL :: lplane                   ! planar option
  LOGICAL :: l_limited_area            

  REAL(wp) :: corio_lat               ! if lplane = .TRUE., center of the f-plane is located 
                                      ! at this geographical latitude

  REAL(wp) :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                       ! for any of the first level child patches,
                                       ! processor splitting will be performed

  LOGICAL  :: lpatch0                  ! If set to .true. an additional patch one
                                       ! level below the root patch is allocated
                                       ! and read so that physics calculations
                                       ! on a coarser grid are possible

  NAMELIST/grid_ctl/ nroot, start_lev, n_dom, lfeedback, lplane, corio_lat, &
                   & parent_id, l_limited_area, patch_weight, lpatch0,      &
                   & lredgrid_phys

CONTAINS
  !>
  !! @brief Initialization of variables that contain grid configuration
  !!
  !! @par Revision History
  !!  Revision History in mo_model_domimp_setup.f90/setup_files (r3965)
  !!  Modification by Constantin Junk, MPI-M (2011-04-05)
  !!  - moved setup_files to mo_grid_nml
  !!  - renamed setup_files to grid_nml_setup
  !!  - restructured grid_nml_setup
  !!
  SUBROUTINE grid_nml_setup
                                                
    INTEGER  :: ist, i, funit
    CHARACTER(len=max_char_length), PARAMETER ::     &
             &  routine = 'mo_grid_nml/grid_nml_setup'
 
    !------------------------------------------------------------
    ! Set up the default values for grid_ctl
    !------------------------------------------------------------
    nroot       = 2
    start_lev   = 4
    n_dom       = 1
    
    ! Note: the first element of parent_id refers to the first nested domain
    DO i = 1, max_dom-1
      parent_id(i) = i
    ENDDO
  
    lfeedback      = .TRUE.
    lplane         = .FALSE.
    l_limited_area = .FALSE.
    corio_lat      = 0.0_wp
    patch_weight   = 0.0_wp
    lpatch0        = .FALSE.
    lredgrid_phys  = .FALSE.
 
    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values in the previous integration. (Don't get confused by
    ! the name of the subroutine.)
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('grid_ctl')
      READ(funit,NML=grid_ctl) ; write(0,*)
      CALL close_tmpfile(funit); write(0,*)
     !! for testing
     !WRITE (0,*) 'contents of namelist ...'
     !WRITE (0,NML=grid_ctl)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('grid_ctl', status=ist)
    IF (ist == POSITIONED) THEN
      READ (nnml, grid_ctl)
    ENDIF

    !---------------------------------------
    ! Check consistency of the parameters
    !---------------------------------------
    ! Reset lfeedback to false for all model domains if lfeedback(1) = false
    IF (.NOT. lfeedback(1)) lfeedback(2:max_dom) = .FALSE.

    ! convert degrees in radiant for the Coriolis latitude
    corio_lat =  corio_lat/rad2deg

    ! set n_dom_start

    IF(lpatch0) THEN
      n_dom_start = 0
    ELSE
      n_dom_start = 1
      lredgrid_phys = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
    ENDIF

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=grid_ctl)
    CALL store_and_close_namelist(funit, 'grid_ctl')

    ! Write the contents of the namelist to an ASCII file.
    ! Probably will be removed later.

    IF(p_pe == p_io) WRITE(nnml_output,nml=grid_ctl)

 END SUBROUTINE grid_nml_setup

END MODULE mo_grid_nml
