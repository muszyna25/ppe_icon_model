!>
!!        Contains the variables to set up the grid configuration
!!
!! @par Revision History
!!  Leonidas Linardakis, MPI-M, 2011/7/7
!!  - Restructuring the namelists
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
MODULE mo_grid_config
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, finish
  USE mo_impl_constants,     ONLY: max_dom, itri, ihex
  USE mo_io_units,           ONLY: filename_max 
  USE mo_physical_constants, ONLY: earth_radius, earth_angular_velocity
  USE mo_parallel_config,    ONLY: division_method, division_file_name

#ifndef NOMPI
! The USE statement below lets this module use the routines from
! mo_read_netcdf_parallel where only 1 processor is reading
! and broadcasting the results
USE mo_read_netcdf_parallel, ONLY:                &
   nf_nowrite, nf_global, nf_noerr, nf_strerror,  &
   nf_open            => p_nf_open,               &
   nf_close           => p_nf_close,              &
   nf_get_att_int     => p_nf_get_att_int ,       &
   nf_get_att_double  => p_nf_get_att_double
#endif

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PRIVATE

  PUBLIC :: init_grid_configuration, get_grid_rescale_factor
  PUBLIC :: max_rad_dom
  
  PUBLIC :: global_cell_type, nroot, start_lev, n_dom, lfeedback,       &
    &       lplane, is_plane_torus, corio_lat, l_limited_area, patch_weight, &
    &       lredgrid_phys, ifeedback_type, start_time, end_time

  PUBLIC :: grid_rescale_factor, grid_length_rescale_factor, &
     & grid_sphere_radius, grid_angular_velocity

  PUBLIC :: namelist_grid_angular_velocity

  PUBLIC :: dynamics_grid_filename,  dynamics_parent_grid_id,     &
    &       radiation_grid_filename, dynamics_radiation_grid_link

! !   PUBLIC :: radiation_grid_distribution
  
  PUBLIC :: n_dom_start, max_childdom     

  PUBLIC :: n_phys_dom

  PUBLIC :: no_of_dynamics_grids, no_of_radiation_grids

  PUBLIC :: use_duplicated_connectivity! , use_dummy_cell_closure
  ! ------------------------------------------------------------------------


#ifdef NOMPI
INCLUDE 'netcdf.inc'
#endif

  INTEGER, PARAMETER  :: max_rad_dom = 3
  ! ------------------------------------------------------------------------
  !Configuration variables
  ! ------------------------------------------------------------------------
  INTEGER  :: global_cell_type
  INTEGER  :: nroot                    ! root division of initial edges
  INTEGER  :: start_lev                ! coarsest bisection level
  INTEGER  :: n_dom                    ! number of model domains, 1=global domain only 
  INTEGER  :: n_dom_start=1 
  INTEGER  :: max_childdom             ! type of feedback (incremental or relaxation)
  INTEGER  :: ifeedback_type
  REAL(wp) :: start_time(max_dom)      ! Time at which execution of a (nested) model domain starts
  REAL(wp) :: end_time(max_dom)        ! Time at which execution of a (nested) model domain terminates
  INTEGER  :: n_phys_dom=1             ! Number of physical domains, computed when reading the patches

  LOGICAL  :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
  LOGICAL  :: lredgrid_phys(max_dom)   ! If set to .true. is calculated on a reduced grid
  LOGICAL  :: l_limited_area

  LOGICAL  :: use_duplicated_connectivity  = .true.  ! if true, the zero connectivity is replaced by the last non-zero value
!  LOGICAL  :: use_dummy_cell_closure = .false.  ! if true then create a dummy cell and connect it to cells and edges with no neigbor
   
!   INTEGER  :: radiation_grid_distribution   ! 0=do nothing
                                       ! 1=redistribute for radiaiton reading from file

  LOGICAL  :: lplane                   ! f-plane option
  LOGICAL  :: is_plane_torus = .false. ! f-plane with doubly periodic boundary==> like a plane torus
  REAL(wp) :: corio_lat                ! Latitude, where the f-plane is located if 
                                       ! lplane or is_plane_torus=.true.

  REAL(wp) :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                       ! for any of the first level child patches,
                                       ! processor splitting will be performed

  REAL(wp) :: grid_rescale_factor = 1.0_wp
  REAL(wp) :: grid_length_rescale_factor = 1.0_wp
!   REAL(wp) :: grid_area_rescale_factor = 1.0_wp
  REAL(wp) :: grid_sphere_radius  = earth_radius
  REAL(wp) :: grid_angular_velocity  = earth_angular_velocity
  REAL(wp) :: namelist_grid_angular_velocity  = earth_angular_velocity

  CHARACTER(LEN=filename_max) :: dynamics_grid_filename(max_dom)
  INTEGER                     :: dynamics_parent_grid_id(max_dom)
  CHARACTER(LEN=filename_max) :: radiation_grid_filename(max_rad_dom)
  INTEGER                     :: dynamics_radiation_grid_link(max_dom)

  INTEGER :: no_of_dynamics_grids  = 0
  INTEGER :: no_of_radiation_grids = 0

  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent variables
  ! -----------------------------------------------------------------------

CONTAINS

 !>
 !!  Initialization of variables that contain grid configuration
 !! @par Revision History
 !!  Leonidas Linardakis, MPI-M, 2011/7/7
 !!  - Restructuring the namelists
  SUBROUTINE init_grid_configuration
                                               
    !local variables
    INTEGER  :: jg, ncid, cell_type
!    INTEGER  :: funit
    LOGICAL  :: file_exists
    CHARACTER(*), PARAMETER :: method_name = "mo_grid_config:init_grid_configuration"

    IF (no_of_dynamics_grids /= 0) &
      CALL finish( method_name, 'should not be called twice')
    
    !-----------------------------------------------------------------------
    ! find out how many grids we have
    ! and check if they exist
    no_of_dynamics_grids  = 0
    no_of_radiation_grids = 0

!     write(0,*) method_name, TRIM(dynamics_grid_filename(1))
    
    jg=1
    DO WHILE (dynamics_grid_filename(jg) /= "")
      INQUIRE (FILE=dynamics_grid_filename(jg), EXIST=file_exists)
      IF (.NOT. file_exists)   THEN
        WRITE (message_text,'(a,a)')  TRIM(dynamics_grid_filename(jg)), &
          " file does not exist"
        CALL finish( TRIM(method_name), TRIM(message_text))
      ENDIF
      jg=jg+1
      IF (jg > max_dom) EXIT
    END DO
    no_of_dynamics_grids  = jg-1
    
    jg=1
    DO WHILE (radiation_grid_filename(jg) /= "")
      INQUIRE (FILE=radiation_grid_filename(jg), EXIST=file_exists)
      IF (.NOT. file_exists)   THEN
        WRITE (message_text,'(a,a)')  TRIM(radiation_grid_filename(jg)), &
          " file does not exist"
        CALL finish( TRIM(method_name), TRIM(message_text))
      ENDIF
      jg=jg+1
      IF (jg > max_rad_dom) EXIT
    END DO
    no_of_radiation_grids = jg-1
    n_dom = no_of_dynamics_grids

    IF (no_of_dynamics_grids < 1) &
      CALL finish( TRIM(method_name), 'no dynamics grid is defined')
      
    ! get here the nroot, eventually it should be moved into the patch info
!     nroot = get_grid_root(dynamics_grid_filename(1))
    CALL nf(nf_open(dynamics_grid_filename(1), nf_nowrite, ncid))
    CALL get_gridfile_root_level(ncid, nroot, start_lev)
    CALL get_gridfile_cell_type(ncid, cell_type)
    CALL nf(nf_close(ncid))

    IF (global_cell_type /= cell_type) &
      CALL finish(method_name, "global_cell_type /= cell_type")

    ! domain geometric properties
!    CALL get_gridfile_rescale_factor(dynamics_grid_filename(1), grid_rescale_factor)
    IF ( grid_rescale_factor <= 0.0_wp ) grid_rescale_factor = 1.0_wp
    CALL get_gridfile_sphere_radius(dynamics_grid_filename(1), grid_sphere_radius)
    grid_sphere_radius = grid_sphere_radius * grid_rescale_factor
    grid_length_rescale_factor = grid_rescale_factor
!     grid_area_rescale_factor   = grid_rescale_factor * grid_rescale_factor
    grid_angular_velocity      = namelist_grid_angular_velocity / grid_rescale_factor
!     write(0,*) "   nroot = ", nroot
    
!     write(0,*) "grid_sphere_radius=",grid_sphere_radius
!     write(0,*) "grid_length_rescale_factor=",grid_length_rescale_factor 
!     write(0,*) "grid_area_rescale_factor=", grid_area_rescale_factor
!     write(0,*) "grid_angular_velocity=", grid_angular_velocity

    IF (no_of_radiation_grids > 0) THEN
      n_dom_start = 0
    ELSE
      n_dom_start = 1
      lredgrid_phys = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
    
      ! the division method starts from 0, shift if there's no 0 grid (ie no reduced radiation)
      DO jg = no_of_dynamics_grids-1, 0, -1
        division_method(jg+1)              = division_method(jg)
        division_file_name(jg+1)           = division_file_name(jg)
      ENDDO

    ENDIF
    
    
    !------------------------------------------------------------
    ! Reset lfeedback to false for all model domains if lfeedback(1) = false
    IF (.NOT. lfeedback(1)) lfeedback(2:max_dom) = .FALSE.
    
    !------------------------------------------------------------
    SELECT CASE (global_cell_type)
    CASE (itri,ihex)
      ! ok
    CASE default
      CALL finish( TRIM(method_name),&
        & 'wrong cell type specifier, "global_cell_type" must be 3 or 6')
    END SELECT
               
!     write(0,*) no_of_dynamics_grids
!     write(0,*) dynamics_grid_filename(1:no_of_dynamics_grids)
!     write(0,*) dynamics_parent_grid_id(1:no_of_dynamics_grids)
!     write(0,*) no_of_radiation_grids
!     IF (no_of_radiation_grids > 0) THEN
!       write(0,*) radiation_grid_filename(1:no_of_radiation_grids)
!       write(0,*) dynamics_radiation_grid_link(1:no_of_dynamics_grids)
!     ENDIF

!     CALL finish("grid_nml_setup","stop")

  END SUBROUTINE init_grid_configuration
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE get_gridfile_root_level( ncid, grid_root, grid_level )
!     CHARACTER(len=*),    INTENT(in)  ::  patch_file   ! name of grid file
    INTEGER,    INTENT(in)     :: ncid
    INTEGER,    INTENT(inout)  :: grid_root, grid_level

    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_root', grid_root))
    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_level', grid_level))

  END SUBROUTINE get_gridfile_root_level
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE get_gridfile_cell_type( ncid, cell_type )
    INTEGER,    INTENT(in)     :: ncid
    INTEGER,    INTENT(inout)  :: cell_type

    INTEGER :: netcd_status
    
    netcd_status = nf_get_att_int(ncid, nf_global,'grid_cell_type', cell_type)
    IF (netcd_status /= nf_noerr) & ! old grid 
      cell_type = global_cell_type

  END SUBROUTINE get_gridfile_cell_type
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
!   SUBROUTINE get_gridfile_rescale_factor( patch_file, rescale_factor )
!     CHARACTER(len=*),    INTENT(in)  ::  patch_file   ! name of grid file
!     REAL(wp),   INTENT(out)          ::  rescale_factor
! 
!     INTEGER :: ncid, netcd_status
! 
!     CALL nf(nf_open(TRIM(patch_file), nf_nowrite, ncid))
!     netcd_status = nf_get_att_double(ncid, nf_global,'earth_rescale_factor', &
!         & rescale_factor)
!     IF (netcd_status /= nf_noerr) rescale_factor = 1.0_wp
!     CALL nf(nf_close(ncid))
! 
!   END SUBROUTINE get_gridfile_rescale_factor
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  SUBROUTINE get_gridfile_sphere_radius( patch_file, sphere_radius )
    CHARACTER(len=*),    INTENT(in)  ::  patch_file   ! name of grid file
    REAL(wp),   INTENT(out)          ::  sphere_radius

    INTEGER :: ncid, netcd_status

    CALL nf(nf_open(TRIM(patch_file), nf_nowrite, ncid))
    netcd_status = nf_get_att_double(ncid, nf_global,'sphere_radius', &
        & sphere_radius )
    IF (netcd_status /= nf_noerr) sphere_radius = earth_radius
    CALL nf(nf_close(ncid))

  END SUBROUTINE get_gridfile_sphere_radius
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  FUNCTION get_grid_rescale_factor( ) result(rescale_factor)
    REAL(wp) ::  rescale_factor

    IF (no_of_dynamics_grids < 1) &
      CALL finish( "get_grid_rescale_factor", 'no dynamics grid is defined')
    rescale_factor = grid_rescale_factor

  END FUNCTION get_grid_rescale_factor
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION get_grid_root( patch_file )
    CHARACTER(len=*),    INTENT(in)  ::  patch_file   ! name of grid file

    INTEGER :: ncid, grid_root

    CALL nf(nf_open(TRIM(patch_file), nf_nowrite, ncid))
    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_root', grid_root))
    CALL nf(nf_close(ncid))

    get_grid_root = grid_root

  END FUNCTION get_grid_root
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE nf(status)
    INTEGER, INTENT(in) :: status
    IF (status /= nf_noerr) THEN
      CALL finish('mo_grid_config netCDF error', nf_strerror(status))
    ENDIF
  END SUBROUTINE nf
  !-------------------------------------------------------------------------

END MODULE mo_grid_config
