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
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length, itri, ihex
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

  PUBLIC :: read_grid_namelist

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary variables
  ! ------------------------------------------------------------------------

  INTEGER    :: nml_nroot                    ! root division of initial edges
  INTEGER    :: nml_start_lev                ! coarsest bisection level
  INTEGER    :: nml_n_dom                    ! number of model domains, 1=global domain only 
  INTEGER    :: nml_n_dom_start=1 
  INTEGER    :: nml_max_childdom
  INTEGER, DIMENSION (max_dom-1) :: nml_parent_id  !ID of parent domain
  ! cell geometry
  ! -------------
  INTEGER    :: nml_cell_type         ! cell type:

  LOGICAL    :: nml_lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
  LOGICAL    :: nml_lredgrid_phys(max_dom)   ! If set to .true. is calculated on a reduced grid
  LOGICAL    :: nml_lplane                   ! planar option
  LOGICAL    :: nml_l_limited_area            

  ! if lplane: latitude at which tangential plane resides
  REAL(wp)   :: nml_corio_lat                ! Center of the f-plane is located 
                                         ! at this geographical latitude

  REAL(wp)   :: nml_patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                         ! for any of the first level child patches,
                                         ! processor splitting will be performed

  LOGICAL    :: nml_lpatch0                  ! If set to .true. an additional patch one
                                         ! level below the root patch is allocated
                                         ! and read so that physics calculations
                                         ! on a coarser grid are possible

  CHARACTER(LEN=filename_max) :: nml_dynamics_grid_filename(max_dom)
  INTEGER                     :: nml_dynamics_parent_grid_id(max_dom)
  CHARACTER(LEN=filename_max) :: nml_radiation_grid_filename(max_dom)
  INTEGER                     :: nml_dynamics_radiation_gridlink(max_dom)

  INTEGER :: nml_no_of_dynamics_grids, nml_no_of_radiation_grids

  NAMELIST /grid_ctl/ nml_nroot, nml_start_lev, nml_n_dom, nml_lfeedback, nml_lplane, &
    & nml_corio_lat, nml_parent_id, nml_l_limited_area, nml_patch_weight, nml_lpatch0,&
    & nml_lredgrid_phys, nml_cell_type,                                  &
    & nml_dynamics_grid_filename,  nml_dynamics_parent_grid_id,         &
    & nml_radiation_grid_filename, nml_dynamics_radiation_gridlink

  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent variables
  ! -----------------------------------------------------------------------

  INTEGER  :: nroot_d
  INTEGER  :: start_lev_d 
  INTEGER  :: n_dom_d
  INTEGER  :: parent_id_d(max_dom-1)

  LOGICAL  :: lfeedback_d(max_dom)
  LOGICAL  :: lplane_d 
  LOGICAL  :: l_limited_area_d 
  LOGICAL  :: lredgrid_phys_d(max_dom)

  REAL(wp) :: corio_lat_d 
  REAL(wp) :: patch_weight_d(max_dom)

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
 !!  - moved setup_files to mo_grid_nml
 !!  - renamed setup_files to grid_nml_setup
 !!  - restructured grid_nml_setup

  SUBROUTINE read_grid_namelist
                                               
    !local variable
    INTEGER  :: i_status, i, jg, jlev, funit
    CHARACTER(filename_max) :: patch_file, gridtype
    INTEGER  ::  patch_level(max_dom)
    LOGICAL :: l_exist
    
    CHARACTER(*), PARAMETER :: method_name = "check_parallel_configuration"


!    CHARACTER(len=max_char_length), PARAMETER :: &
!      &  routine = 'mo_grid_nml/run_grid_setup'

    !-----------------------------------------------------------------------
    ! clear grid filenames and hierarchy
    nml_no_of_dynamics_grids  = 0
    nml_no_of_radiation_grids = 0
    DO i = 1, max_dom
      nml_dynamics_grid_filename(i)   = ""
      nml_radiation_grid_filename(i)  = ""
      nml_dynamics_parent_grid_id(i)  = 0
      nml_dynamics_radiation_gridlink(i) = 0
    ENDDO
    
    !------------------------------------------------------------
    ! 3.0 set up the default values for grid_ctl
    !------------------------------------------------------------

    nml_nroot       = 2
    nml_start_lev   = 4
    nml_n_dom       = 1
    nml_i_cell_type = itri
    
    ! Note: the first element of parent_id refers to the first nested domain
    DO i = 1, max_dom-1
      nml_parent_id(i) = i
    ENDDO
  
    nml_lfeedback   = .TRUE.
    nml_lplane      = .FALSE.
    nml_l_limited_area = .FALSE.
    nml_corio_lat   = 0.0_wp
    nml_patch_weight= 0.0_wp
    nml_lpatch0     = .FALSE.
    nml_lredgrid_phys = .FALSE.

    ! copy default values to "default" variables
    nroot_d     = nml_nroot
    start_lev_d = nml_start_lev
    n_dom_d     = nml_n_dom
    parent_id_d = nml_parent_id
    lfeedback_d = nml_lfeedback
    lplane_d    = nml_lplane
    l_limited_area_d = nml_l_limited_area
    lredgrid_phys_d  = nml_lredgrid_phys
    corio_lat_d = nml_corio_lat
    patch_weight_d = nml_patch_weight

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
    ! funit = open_and_restore_namelist('grid_ctl')
    ! READ(funit,NML=grid_ctl)
    ! CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------
    ! 4.0 Read the namelist
    !------------------------------------------------------------
    ! (done so far by all MPI processes)

    CALL position_nml ('grid_ctl', status=i_status)
    IF (i_status == POSITIONED) THEN
      READ (nnml, grid_ctl)
    ENDIF

    !------------------------------------------------------------
    ! 5.0 check the consistency of the parameters
    !------------------------------------------------------------
    SELECT CASE (i_cell_type)
    CASE (itri,ihex)
      ! ok
    CASE default
      CALL finish( TRIM(method_name),&
        & 'wrong cell type specifier, "i_cell_type" must be 3 or 6')
    END SELECT


    ! Reset lfeedback to false for all model domains if lfeedback(1) = false
!     IF (.NOT. lfeedback(1)) lfeedback(2:max_dom) = .FALSE.

    !-----------------------------------------------------
    ! Set dependent variables
    !-----------------------------------------------------
    ! convert degrees in radiant for the Coriolis latitude

!     corio_lat =  corio_lat/rad2deg

    ! set n_dom_start

!     IF(lpatch0) THEN
!       n_dom_start = 0
!     ELSE
!       n_dom_start = 1
!       lredgrid_phys = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
!     ENDIF


    IF (dynamics_grid_filename(1) == "") THEN
      ! dynamics_grid_filename not filled
      ! we have an old style namelist

      ! fill dynamics_grid_filename
      ! fill level and parent ids
      nml_patch_level(1) = start_lev
      nml_dynamics_parent_grid_id(1) = 0       
      DO jg = 2, n_dom
        nml_dynamics_parent_grid_id(jg) = parent_id(jg-1)
        nml_patch_level(jg) = patch_level(dynamics_parent_grid_id(jg))+1
      ENDDO 
    
      ! fill the grid prefix
      IF (lplane) THEN
        nml_gridtype='plan'
      ELSE
        nml_gridtype='icon'
      END IF
    
      DO jg = 1, n_dom
        jlev = nml_patch_level(jg)
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
        nml_dynamics_grid_filename(jg) = patch_file
      ENDDO

      IF (n_dom_start == 0) THEN
        ! fill radiation_grid_filename
        jlev = start_lev-1
        jg=0        
        WRITE (patch_file,'(a,a,i0,2(a,i2.2),a)') &
            & TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
        nml_radiation_grid_filename(1) = patch_file
        nml_dynamics_radiation_grid_link(1) = 1
      ENDIF
    
    ENDIF

    ! find out how many grids we have
    jg=1
    DO WHILE (nml_dynamics_grid_filename(jg) /= "")
      jg=jg+1
    END DO
    nml_no_of_dynamics_grids  = jg-1
    jg=1
    DO WHILE (nml_radiation_grid_filename(jg) /= "")
      jg=jg+1
    END DO
    nml_no_of_radiation_grids = jg-1
    nml_n_dom = no_of_dynamics_grids
    IF (nml_no_of_radiation_grids > 0) THEN
      nml_n_dom_start = 0
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
    IF(p_pe == p_io) WRITE(nnml_output,nml=grid_ctl)

    
  END SUBROUTINE read_grid_namelist

END MODULE mo_grid_nml
