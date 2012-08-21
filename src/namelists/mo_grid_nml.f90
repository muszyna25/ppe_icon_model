!>
!!  Contains the variables to set up the grid configuration
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
MODULE mo_grid_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,               ONLY: wp
! USE mo_exception,          ONLY: message, message_text, finish
  USE mo_io_units,           ONLY: nnml, nnml_output,filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_impl_constants,     ONLY: max_dom, itri
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_physical_constants, ONLY: earth_angular_velocity
!  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist,  &
!                                    & open_and_restore_namelist, close_tmpfile

  USE mo_grid_config,        ONLY:                                         &
    & config_global_cell_type             => global_cell_type,             &
    & config_lfeedback                    => lfeedback,                    &
    & config_ifeedback_type               => ifeedback_type,               &
    & config_start_time                   => start_time,                   &
    & config_end_time                     => end_time,                     &
    & config_lplane                       => lplane,                       &
    & config_corio_lat                    => corio_lat,                    &
    & config_l_limited_area               => l_limited_area,               &
    & config_patch_weight                 => patch_weight,                 &
    & config_lredgrid_phys                => lredgrid_phys,                &
    & config_dynamics_grid_filename       => dynamics_grid_filename,       &
    & config_dynamics_parent_grid_id      => dynamics_parent_grid_id,      &
    & config_radiation_grid_filename      => radiation_grid_filename,      &
    & config_dyn_radiation_grid_link      => dynamics_radiation_grid_link, &
    & config_grid_rescale_factor          => grid_rescale_factor,          &
    & config_grid_angular_velocity        => namelst_grid_angular_velocity,&
!     & config_radiation_grid_distrib       => radiation_grid_distribution,  &
    & check_grid_configuration, max_rad_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_grid_namelist
 !PUBLIC :: fill_grid_nml_configure
  
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary variables
  ! ------------------------------------------------------------------------


  CONTAINS

 !>
 !!  Initialization of grid namelist variables
 !! @par Revision History
 !!  Revision History in mo_model_domimp_setup.f90/setup_files (r3965)
 !!  Modification by Constantin Junk, MPI-M (2011-04-05)
 !!  - moved setup_files to mo_grid_nml
 !!  - renamed setup_files to grid_nml_setup
 !!  - restructured grid_nml_setup
 !!  Leonidas Linardakis, MPI-M, 2011/7/7
 !!  - Restructuring the namelists
 !!
  SUBROUTINE read_grid_namelist( filename )
    
    CHARACTER(LEN=*), INTENT(IN) :: filename                                           
    INTEGER  :: i_status, i

    INTEGER    :: cell_type                ! cell type:

    LOGICAL    :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
    INTEGER    :: ifeedback_type           ! type of feedback (incremental or relaxation)
    REAL(wp)   :: start_time(max_dom)      ! time at which execution of a (nested) model domain starts
    REAL(wp)   :: end_time(max_dom)        ! time at which execution of a (nested) model domain terminates
    LOGICAL    :: lredgrid_phys(max_dom)   ! If set to .true. is calculated on a reduced grid
    LOGICAL    :: l_limited_area            

    LOGICAL    :: lplane                   ! f-plane option
    REAL(wp)   :: corio_lat                ! Latitude, where the f-plane is located if lplane=.true.
  
    REAL(wp)   :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                          ! for any of the first level child patches,
                                          ! processor splitting will be performed

    CHARACTER(LEN=filename_max) :: dynamics_grid_filename(max_dom)
    INTEGER                     :: dynamics_parent_grid_id(max_dom)
    CHARACTER(LEN=filename_max) :: radiation_grid_filename(max_rad_dom)
    INTEGER                     :: dynamics_radiation_grid_link(max_dom)
      
!     INTEGER                     :: radiation_grid_distribution
    
    REAL(wp) :: grid_rescale_factor, grid_angular_velocity


    NAMELIST /grid_nml/ cell_type, lfeedback, ifeedback_type,      &
      &  lplane, corio_lat, l_limited_area, grid_rescale_factor,   &
      &  patch_weight, lredgrid_phys, start_time, end_time,        &
      &  dynamics_grid_filename,  dynamics_parent_grid_id,    &
      &  radiation_grid_filename, dynamics_radiation_grid_link
!       &  radiation_grid_distribution



!    INTEGER  :: funit
!    INTEGER  :: jg, jlev
!    CHARACTER(filename_max) :: patch_file, gridtype
!    INTEGER  ::  patch_level(max_dom)
!    LOGICAL :: l_exist

    !------------------------------------------------------------
    !  set up the default values for grid_nml
    !------------------------------------------------------------
    DO i = 1, max_dom
      dynamics_grid_filename(i)   = ""
      dynamics_parent_grid_id(i)  = i-1
      dynamics_radiation_grid_link(i) = 0
    ENDDO
    dynamics_radiation_grid_link(1) = 1
    DO i = 1, max_rad_dom
      radiation_grid_filename(i)  = ""
    ENDDO

    cell_type   = itri
      
    lfeedback   = .TRUE.
    ifeedback_type = 1
    start_time(:) = 0._wp
    end_time(:)   = 1.e30_wp
    lplane      = .FALSE.
    l_limited_area = .FALSE.
    corio_lat   = 0.0_wp
    patch_weight= 0.0_wp
    lredgrid_phys = .FALSE.
!     radiation_grid_distribution = 0

    !----------------------------------------------------------------
    grid_rescale_factor   = 1.0_wp
    grid_angular_velocity = earth_angular_velocity

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (is_restart_run()) THEN
    ! funit = open_and_restore_namelist('grid_nml')
    ! READ(funit,NML=grid_nml)
    ! CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------
    ! Read the namelist (done so far by all MPI processes)
    !------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('grid_nml', status=i_status)
    IF (i_status == POSITIONED) THEN
      READ (nnml, grid_nml)
    ENDIF
    CALL close_nml

    ! convert degrees in radiant for the Coriolis latitude
    corio_lat =  corio_lat/rad2deg

    ! Reset start and end times for global domain in order to avoid 
    ! interferences with output flow control
    start_time(1) = 0._wp
    end_time(1) = 1.e30_wp

    !-----------------------------------------------------
    !  Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
    ! funit = open_tmpfile()
    ! WRITE(funit,NML=grid_nml)
    ! CALL store_and_close_namelist(funit, 'grid_nml')

      ! write the contents of the namelist to an ASCII file
      WRITE(nnml_output,nml=grid_nml)
    ENDIF


    !-----------------------------------------------------
    config_global_cell_type  = cell_type
    config_lfeedback         = lfeedback
    config_ifeedback_type    = ifeedback_type
    config_start_time        = start_time
    config_end_time          = end_time
    config_lplane            = lplane
    config_corio_lat         = corio_lat
    config_l_limited_area    = l_limited_area
    config_patch_weight      = patch_weight
    config_lredgrid_phys     = lredgrid_phys
!     config_radiation_grid_distrib  = radiation_grid_distribution
    config_dynamics_grid_filename  = dynamics_grid_filename
    config_dynamics_parent_grid_id = dynamics_parent_grid_id
    config_radiation_grid_filename = radiation_grid_filename
    config_dyn_radiation_grid_link = dynamics_radiation_grid_link
    config_grid_rescale_factor     = grid_rescale_factor
    config_grid_angular_velocity   = grid_angular_velocity
    
    ! check the configuration
    CALL check_grid_configuration()
       
  END SUBROUTINE read_grid_namelist
  !-----------------------------------------------------------------------
  
  
END MODULE mo_grid_nml
