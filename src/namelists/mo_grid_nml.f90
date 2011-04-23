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
  USE mo_exception,          ONLY: message, message_text
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, POSITIONED
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_math_constants,     ONLY: rad2deg

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC
  PUBLIC :: grid_nml_setup

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary variables
  ! ------------------------------------------------------------------------

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

  NAMELIST /grid_ctl/ nroot, start_lev, n_dom, lfeedback, lplane, corio_lat,     &
                      parent_id, l_limited_area, patch_weight, lpatch0,          &
                      lredgrid_phys

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

 SUBROUTINE grid_nml_setup
                                               
  !local variable
   INTEGER  :: i_status, i

   CHARACTER(len=max_char_length), PARAMETER :: &
     &  routine = 'mo_grid_nml/run_grid_setup'

!-----------------------------------------------------------------------
   
   !------------------------------------------------------------
   ! 3.0 set up the default values for grid_ctl
   !------------------------------------------------------------

   nroot       = 2
   start_lev   = 4
   n_dom       = 1
   
   ! Note: the first element of parent_id refers to the first nested domain
   DO i = 1, max_dom-1
     parent_id(i) = i
   ENDDO
 
   lfeedback   = .TRUE.
   lplane      = .FALSE.
   l_limited_area = .FALSE.
   corio_lat   = 0.0_wp
   patch_weight= 0.0_wp
   lpatch0     = .FALSE.
   lredgrid_phys = .FALSE.

   ! copy default values to "default" variables
   nroot_d     = nroot
   start_lev_d = start_lev
   n_dom_d     = n_dom
   parent_id_d = parent_id
   lfeedback_d = lfeedback
   lplane_d    = lplane
   l_limited_area_d = l_limited_area
   lredgrid_phys_d  = lredgrid_phys
   corio_lat_d = corio_lat
   patch_weight_d = patch_weight

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

  ! Reset lfeedback to false for all model domains if lfeedback(1) = false
  IF (.NOT. lfeedback(1)) lfeedback(2:max_dom) = .FALSE.

  ! write the contents of the namelist to an ASCII file

  IF(p_pe == p_io) WRITE(nnml_output,nml=grid_ctl)

  ! convert degrees in radiant for the Coriolis latitude
  corio_lat =  corio_lat/rad2deg

  ! set n_dom_start

  IF(lpatch0) THEN
    n_dom_start = 0
  ELSE
    n_dom_start = 1
    lredgrid_phys = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
  ENDIF


 END SUBROUTINE grid_nml_setup

END MODULE mo_grid_nml
