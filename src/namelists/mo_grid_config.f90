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
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_dom, max_char_length, itri, ihex, max_char_length
  USE mo_io_units,           ONLY: filename_max 
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_math_constants,     ONLY: rad2deg

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PRIVATE

  PUBLIC :: check_grid_configuration
  
  PUBLIC :: global_cell_type, nroot, start_lev, n_dom, lfeedback,       &
    &       lplane, corio_lat, parent_id, l_limited_area, patch_weight, &
    &       lredgrid_phys

  PUBLIC :: dynamics_grid_filename,  dynamics_parent_grid_id,     &
    &       radiation_grid_filename, dynamics_radiation_grid_link

  PUBLIC :: n_dom_start, max_childdom     

  PUBLIC :: no_of_dynamics_grids, no_of_radiation_grids

  ! ------------------------------------------------------------------------
  !Configuration variables
  ! ------------------------------------------------------------------------
  INTEGER  :: global_cell_type
  INTEGER  :: nroot                    ! root division of initial edges
  INTEGER  :: start_lev                ! coarsest bisection level
  INTEGER  :: n_dom                    ! number of model domains, 1=global domain only 
  INTEGER  :: n_dom_start=1 
  INTEGER  :: max_childdom
  INTEGER  :: parent_id(max_dom-1)     !ID of parent domain

  LOGICAL  :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
  LOGICAL  :: lredgrid_phys(max_dom)   ! If set to .true. is calculated on a reduced grid
  LOGICAL  :: l_limited_area            

  LOGICAL  :: lplane                   ! f-plane option
  REAL(wp) :: corio_lat                ! Latitude at which the f-plane is located 

  REAL(wp) :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                       ! for any of the first level child patches,
                                       ! processor splitting will be performed

!   LOGICAL :: lpatch0                  ! If set to .true. an additional patch one
!                                       ! level below the root patch is allocated
!                                       ! and read so that physics calculations
!                                       ! on a coarser grid are possible

  CHARACTER(LEN=filename_max) :: dynamics_grid_filename(max_dom)
  INTEGER                     :: dynamics_parent_grid_id(max_dom)
  CHARACTER(LEN=filename_max) :: radiation_grid_filename(max_dom)
  INTEGER                     :: dynamics_radiation_grid_link(max_dom)

  INTEGER :: no_of_dynamics_grids, no_of_radiation_grids

  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent variables
  ! -----------------------------------------------------------------------

CONTAINS

 !>
 !!  Initialization of variables that contain grid configuration
 !! @par Revision History
 !!  Leonidas Linardakis, MPI-M, 2011/7/7
 !!  - Restructuring the namelists
  SUBROUTINE check_grid_configuration
                                               
    !local variables
    INTEGER  :: jg
!    INTEGER  :: funit
    LOGICAL  :: file_exists
    CHARACTER(*), PARAMETER :: method_name = "mo_grid_config:check_grid_configuration"

    !-----------------------------------------------------------------------
    ! find out how many grids we have
    ! and check if they exist
    no_of_dynamics_grids  = 0
    no_of_radiation_grids = 0
    
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
      IF (jg > max_dom) EXIT
    END DO
    no_of_radiation_grids = jg-1
    n_dom = no_of_dynamics_grids

    ! some checks
    IF (no_of_dynamics_grids < 1) &
      CALL finish( TRIM(method_name), 'no dynamics grid is defined')
        
    IF (no_of_radiation_grids > 0) THEN
      n_dom_start = 0
    ELSE
      n_dom_start = 1
      lredgrid_phys = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
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
    

    !-----------------------------------------------------
    ! Set dependent variables
    !-----------------------------------------------------

           
!     write(0,*) no_of_dynamics_grids
!     write(0,*) dynamics_grid_filename(1:no_of_dynamics_grids)
!     write(0,*) dynamics_parent_grid_id(1:no_of_dynamics_grids)
!     write(0,*) no_of_radiation_grids
!     IF (no_of_radiation_grids > 0) THEN
!       write(0,*) radiation_grid_filename(1:no_of_radiation_grids)
!       write(0,*) dynamics_radiation_grid_link(1:no_of_dynamics_grids)
!     ENDIF

!     CALL finish("grid_nml_setup","stop")

    !-----------------------------------------------------                                        
    ! Store the namelist for restart                                                              
    !-----------------------------------------------------                                        
   !funit = open_tmpfile()                                                                        
   !WRITE(funit,NML=grid_nml)                                                                     
   !CALL store_and_close_namelist(funit, 'grid_nml')

    ! write the contents of the namelist to an ASCII file
!     IF(p_pe == p_io) WRITE(nnml_output,nml=grid_nml)

    
  END SUBROUTINE check_grid_configuration

END MODULE mo_grid_config
