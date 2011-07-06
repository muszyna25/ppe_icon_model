!>
!!        Contains the variables to set up the grid configuration
!!
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
!!
MODULE mo_grid_configuration
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output,filename_max 
  USE mo_namelist,           ONLY: position_nml, POSITIONED
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PRIVATE

  PUBLIC :: n_dom,n_dom_start, global_cell_type
  
  PUBLIC :: dynamics_grid_filename,  dynamics_parent_grid_id,  &
    & radiation_grid_filename, dynamics_radiation_grid_link,   &
    & no_of_dynamics_grids, no_of_radiation_grids, max_childdom     
  PUBLIC :: l_limited_area, nroot, start_lev, lfeedback, &
    & lplane, corio_lat, parent_id, lredgrid_phys

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary variables
  ! ------------------------------------------------------------------------
  INTEGER    :: global_cell_type
  INTEGER    :: nroot                    ! root division of initial edges
  INTEGER    :: start_lev                ! coarsest bisection level
  INTEGER    :: n_dom                    ! number of model domains, 1=global domain only 
  INTEGER    :: n_dom_start=1 
  INTEGER    :: max_childdom
  INTEGER, DIMENSION (max_dom-1) :: parent_id  !ID of parent domain

  LOGICAL    :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
  LOGICAL    :: lredgrid_phys(max_dom)   ! If set to .true. is calculated on a reduced grid
  LOGICAL    :: lplane                   ! planar option
  LOGICAL    :: l_limited_area            

  ! if lplane: latitude at which tangential plane resides
  REAL(wp)   :: corio_lat                ! Center of the f-plane is located 
                                         ! at this geographical latitude

  REAL(wp)   :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                         ! for any of the first level child patches,
                                         ! processor splitting will be performed

  LOGICAL    :: lpatch0                  ! If set to .true. an additional patch one
                                         ! level below the root patch is allocated
                                         ! and read so that physics calculations
                                         ! on a coarser grid are possible

  CHARACTER(LEN=filename_max) :: dynamics_grid_filename(max_dom)
  INTEGER                     :: dynamics_parent_grid_id(max_dom)
  CHARACTER(LEN=filename_max) :: radiation_grid_filename(max_dom)
  INTEGER                     :: dynamics_radiation_grid_link(max_dom)

  INTEGER :: no_of_dynamics_grids, no_of_radiation_grids

  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent variables
  ! -----------------------------------------------------------------------

  CONTAINS
!
!

 !>
 !!  Initialization of variables that contain grid configuration
 !!
 !!
 !! @par Revision History
 !!  Revision History in mo_model_domimp_setup.f90/setup_files (r3965)
 !!  Modification by Constantin Junk, MPI-M (2011-04-05)
 !!  - moved setup_files to mo_grid_configuration
 !!  - renamed setup_files to grid_nml_setup
 !!  - restructured grid_nml_setup

  SUBROUTINE check_grid_configuration
                                               
    !local variable
    INTEGER  :: i_status, i, jg, jlev, funit
    CHARACTER(filename_max) :: patch_file, gridtype
    INTEGER  ::  patch_level(max_dom)
    LOGICAL :: l_exist


!    CHARACTER(len=max_char_length), PARAMETER :: &
!      &  routine = 'mo_grid_configuration/run_grid_setup'

    !-----------------------------------------------------------------------
    ! clear grid filenames and hierarchy
    no_of_dynamics_grids  = 0
    no_of_radiation_grids = 0
    
    ! Note: the first element of parent_id refers to the first nested domain
    DO i = 1, max_dom-1
      parent_id(i) = i
    ENDDO
  

    !------------------------------------------------------------
    ! 5.0 check the consistency of the parameters
    !------------------------------------------------------------
    ! Reset lfeedback to false for all model domains if lfeedback(1) = false
    IF (.NOT. lfeedback(1)) lfeedback(2:max_dom) = .FALSE.

    !-----------------------------------------------------
    ! Set dependent variables
    !-----------------------------------------------------
    ! convert degrees in radiant for the Coriolis latitude
    corio_lat =  corio_lat/rad2deg

    ! set n_dom_start
    IF(lpatch0) THEN
      n_dom_start = 0
    ELSE
      n_dom_start = 1
      lredgrid_phys = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
    ENDIF


    IF (dynamics_grid_filename(1) == "") THEN
      ! dynamics_grid_filename not filled
      ! we have an old style namelist

      ! fill dynamics_grid_filename
      ! fill level and parent ids
      patch_level(1) = start_lev
      dynamics_parent_grid_id(1) = 0       
      DO jg = 2, n_dom
        dynamics_parent_grid_id(jg) = parent_id(jg-1)
        patch_level(jg) = patch_level(dynamics_parent_grid_id(jg))+1
      ENDDO 
    
      ! fill the grid prefix
      IF (lplane) THEN
        gridtype='plan'
      ELSE
        gridtype='icon'
      END IF
    
      DO jg = 1, n_dom
        jlev = patch_level(jg)
        ! Allow file names without "DOM" specifier if n_dom=1.
        IF (n_dom == 1) THEN
          ! Check if file name without "DOM" specifier exists.
          WRITE (patch_file,'(a,a,i0,a,i2.2,a)') &
              & TRIM(gridtype),'R',nroot,'B',jlev,'-grid.nc'
          INQUIRE (FILE=patch_file, EXIST=l_exist)
          ! Otherwise use file name with "DOM" specifier
          IF (.NOT. l_exist)                                           &
              & WRITE (patch_file,'(a,a,i0,2(a,i2.2),a)')              &
              & TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
        ELSE
          ! n_dom >1 --> "'_DOM',jg" required in file name
          WRITE (patch_file,'(a,a,i0,2(a,i2.2),a)') &
              & TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
        ENDIF
        dynamics_grid_filename(jg) = patch_file
      ENDDO

      IF (n_dom_start == 0) THEN
        ! fill radiation_grid_filename
        jlev = start_lev-1
        jg=0        
        WRITE (patch_file,'(a,a,i0,2(a,i2.2),a)') &
            & TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
        radiation_grid_filename(1) = patch_file
        dynamics_radiation_grid_link(1) = 1
      ENDIF
    
    ENDIF

    ! find out how many grids we have
    jg=1
    DO WHILE (dynamics_grid_filename(jg) /= "")
      jg=jg+1
    END DO
    no_of_dynamics_grids  = jg-1
    jg=1
    DO WHILE (radiation_grid_filename(jg) /= "")
      jg=jg+1
    END DO
    no_of_radiation_grids = jg-1
    n_dom = no_of_dynamics_grids
    IF (no_of_radiation_grids > 0) THEN
      n_dom_start = 0
    ENDIF
           
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
   !WRITE(funit,NML=grid_ctl)                                                                     
   !CALL store_and_close_namelist(funit, 'grid_ctl')

    ! write the contents of the namelist to an ASCII file
!     IF(p_pe == p_io) WRITE(nnml_output,nml=grid_ctl)

    
  END SUBROUTINE check_grid_configuration

END MODULE mo_grid_configuration
