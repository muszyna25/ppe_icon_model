!>
!! @page pagepostprogenf901 Main program to calculate normalized errors as in
!!
!! @author
!!     P. R\'\i{}podas
!!     (DWD)
!!
!!
!! @date 2009-12-14 07:35:01 UTC
!!

!>
!!   This program calculates the l1, l2 and linf normalized errors for the.
!!
!!   This program calculates the l1, l2 and linf normalized errors for the
!!  geopotential, vorticity and/or wind fields produced by ICOSWM.
!!   The errors are calculated respect to a reference. The model and reference
!!  fields should be at the same grid points and the grid points should be in
!!  the same order. (The output of the NCAR STSWM must be interpolated to the
!!  ICOSWM grid for cases 5 and 6)
!!   It reads information from the NAMELIST_ICON file used to run ICOSWM and
!!  reads the patches information from the patch files
!!   It also reads information from the file NAMELIST_POSTPRO.
!!
!! @par Revision History
!! Original version by P. R\\'\\i{}podas (2007-07)
!! Modified by Th. Heinze, DWD (2007-08-09):
!! - adapted to modified module mo_err_norm
!! Modified by P. R\\'\\i{}podas (2009-03):
!! - adapted to the restructured code
!! Modified by P. R\\'\\i{}podas (2009-07):
!! - simplified version to use with the pressure
!!   levels of ICOHDC respect to the T799 ECHAM5
!! Modification by Daniel Reinert, DWD (2010-07-19)
!! - adapted to new type 'external_data'
!!
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
PROGRAM postpro_gen
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!   (Williamson et al. 1992) for test cases without analytical
!    solution
!
!
!
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind
USE mo_model_domain,         ONLY: t_patch
USE mo_ext_data,             ONLY: t_external_data, construct_ext_data
USE mo_model_domain_import,  ONLY: n_dom, nroot, parent_id,      &
                                 & start_lev, import_patches,    &
                                 & grid_nml_setup, destruct_patches
USE mo_hydro_testcases,      ONLY: setup_testcase
USE mo_icoham_dyn_types,     ONLY: t_hydro_atm
USE mo_namelist
USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, SUCCESS
USE mo_io_units,             ONLY: nnml, nout
USE mo_ocean_nml,            ONLY: setup_ocean_nml
USE mo_dynamics_nml,         ONLY: dynamics_nml_setup
USE mo_diffusion_nml,        ONLY: diffusion_nml_setup
USE mo_run_nml,              ONLY: run_nml_setup,ldynamics, ltestcase, latmosphere, &
                                 & lshallow_water, lhydrostatic, locean
USE mo_exception,            ONLY: message,finish
USE mo_err_norm,             ONLY: t_errors, err_norm
USE mo_io_vlist,             ONLY: de_reshape1
USE mo_vertical_coord_table, ONLY: init_vertical_coord_table
USE mo_ha_stepping,          ONLY: prepare_ha_integration
USE mo_intp_state,           ONLY: setup_interpol,               &
                                   construct_2d_interpol_state
USE mo_interpolation,        ONLY: t_int_state
!USE mo_hierarchy_management, ONLY: init_fbk_wgt
USE mo_grf_interpolation,    ONLY: setup_gridref, t_gridref_state, &
                                    construct_2d_gridref_state


IMPLICIT NONE

! CHARACTER(len=*), PARAMETER :: version = '$Id$'

! !LOCAL VARIABLES


CHARACTER (MAX_CHAR_LENGTH) :: varname
CHARACTER (MAX_CHAR_LENGTH) :: icohdc_file, ref_file  !files with the model and
                                                      ! reference "varname" fields


TYPE(t_patch), TARGET, ALLOCATABLE        :: p_patch(:)           ! patch information

TYPE(t_patch), POINTER                    :: p_single_patch    => NULL()
TYPE(t_hydro_atm), TARGET, ALLOCATABLE  :: p_hydro_state(:)
TYPE(t_int_state), TARGET, ALLOCATABLE    :: p_int_state(:)
TYPE(t_gridref_state), TARGET, ALLOCATABLE:: p_grf_state(:)
TYPE(t_external_data), ALLOCATABLE        :: ext_data(:) !< external data


INTEGER                :: istatus, ic, jg, lev
INTEGER                :: ncells, nverts, nblks_c, nblks_v
INTEGER                :: ist

REAL(wp)               :: z_lon, z_lat, z_lon_ref, z_lat_ref

REAL(wp)               :: z_var_diff

REAL(wp), ALLOCATABLE  :: z_var(:), z_var_ref(:)

REAL(wp), ALLOCATABLE  :: parea(:), darea(:)

TYPE(t_errors)           :: z_err_var

!

! Namelist groups to be read from NAMELIST_POSTPRO

NAMELIST /setvar/ varname
NAMELIST /ifiles/ icohdc_file, ref_file

!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! set default values
!--------------------------------------------------------------------

varname     = "VOR"
ref_file = "analytical.txt"

!-----------------------------------------------------------------------
! open namelist file NAMELIST_ICON
!-----------------------------------------------------------------------
CALL open_nml( "NAMELIST_ICON" )

CALL run_nml_setup               !reads namelist group run_nml


  IF (locean) THEN
     CALL setup_ocean_nml
  END IF

  IF (ltestcase) THEN
    IF ((latmosphere .AND. lhydrostatic).OR. lshallow_water) THEN
      CALL setup_testcase
    ENDIF
!    IF (latmosphere .AND. .NOT. lhydrostatic) THEN
!      CALL setup_nh_testcase
!    ENDIF
  ENDIF
CALL grid_nml_setup             !reads namelist group grid_ctl


!-----------------------------------------------------------------------
! import patches
!-----------------------------------------------------------------------
 !check patch allocation status
  IF ( ALLOCATED(p_patch) ) THEN
    CALL finish('postpro', 'patch already allocated')
  END IF
  !
  ! allocate patch array to start patch construction
  ALLOCATE(p_patch(n_dom), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish('postpro', 'allocation of patch failed')
  ENDIF
  CALL import_patches(p_patch)

  IF (ldynamics) THEN
    CALL diffusion_nml_setup(n_dom,parent_id)
    CALL dynamics_nml_setup(n_dom)
  END IF
  IF (n_dom > 1) THEN
    CALL setup_gridref
  ENDIF
  CALL setup_interpol(p_patch)

  !
  ! allocate type for interpolation state
  !
  ALLOCATE (p_int_state(n_dom), &
            p_grf_state(n_dom),STAT=ist)
  IF (ist /= SUCCESS) THEN
    CALL finish('hydro_atmos','allocation for ptr_int_state failed')
  ENDIF

  CALL construct_2d_interpol_state(p_patch, p_int_state)
  IF (n_dom > 1) THEN
    CALL construct_2d_gridref_state (p_patch, p_grf_state)
  ENDIF
!-----------------------
!
!-----------------------
    CALL init_vertical_coord_table(p_patch(1)%nlev)
      ALLOCATE (p_hydro_state(n_dom), STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish('control_model','allocation for p_hydro_state failed')
      ENDIF

    !------------------------------------------------------------------
    ! Create and read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=ist)
    IF (ist /= success) THEN
      CALL finish('postpro','allocation for ext_data failed')
    ENDIF

    CALL construct_ext_data (p_patch, ext_data)

    CALL prepare_ha_integration(p_patch, p_int_state, p_grf_state, p_hydro_state)


CALL close_nml
!-----------------------------------------------------------------------
! open NAMELIST_POSTPRO_GEN
!-----------------------------------------------------------------------

CALL open_nml( "NAMELIST_POSTPRO_GEN" )


!-----------------------------------------------------------------------
! read namelist group setvar
!-----------------------------------------------------------------------
CALL position_nml ('setvar', status=istatus)
SELECT CASE (istatus)
CASE (POSITIONED)
  READ(nnml, setvar)
END SELECT


!-----------------------------------------------------------------------
! read namelist group ifiles
!-----------------------------------------------------------------------
  CALL position_nml ('ifiles', status=istatus)
  SELECT CASE (istatus)
  CASE (POSITIONED)
    READ(nnml, ifiles)
  END SELECT
  icohdc_file =TRIM(ADJUSTL(icohdc_file))
  ref_file    =TRIM(ADJUSTL(ref_file))

CALL close_nml


!-----------------------------------------------------------------------
!Start to calculate errors for each level
!-----------------------------------------------------------------------
! For the moment only for the global domain

!DO jg=1,n_dom

DO jg=1,1

  lev= start_lev + jg - 1
  WRITE(nout,*)

  p_single_patch => p_patch(jg)
!---------------------------------------------------------------------
! Areas of the primal and dual grids (weights for the error calculation)
!---------------------------------------------------------------------
  ncells = p_single_patch%n_patch_cells
  nverts = p_single_patch%n_patch_verts
  nblks_c = p_single_patch%nblks_c
  nblks_v = p_single_patch%nblks_v


  ALLOCATE (parea(ncells))    !primal cell areas
  ALLOCATE (darea(nverts))    !dual cell areas


  CALL de_reshape1( p_single_patch%cells%area, parea )
  CALL de_reshape1( p_single_patch%verts%dual_area, darea )


  ALLOCATE (z_var(ncells), z_var_ref(ncells))


!---------------------------------------------------------------------
! errors for the "varname" field
!---------------------------------------------------------------------

    !open ICOSWM and STSWM files
    OPEN(unit=11, FILE=TRIM(icohdc_file), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file icohdc_file failed')
      CALL finish  ( 'postpro', 'opening file icohdc_file failed')
    END IF

    OPEN(unit=12, FILE=TRIM(ref_file), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file ref_file failed')
      CALL finish  ( 'postpro', 'opening file ref_file failed')
    END IF
    OPEN(unit=13, FILE="var_diff.gmt", status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file var_diff.gmt failed')
      CALL finish  ( 'postpro', 'opening file var_diff.gmt failed')
    END IF



    DO ic=1,ncells

      READ(11,*) z_lon, z_lat, z_var(ic)
      READ(12,*) z_lon_ref, z_lat_ref, z_var_ref(ic)
      z_var_diff= z_var(ic)- z_var_ref(ic)
      WRITE(13,'(2f24.14,g24.14)') z_lon, z_lat, z_var_diff

    ENDDO

    CLOSE(13)
    CLOSE(unit=11)
    CLOSE(unit=12)

    CALL err_norm( ncells, z_var_ref, z_var, parea, z_err_var)

   ! WRITE(nout,'(a,i1,a,i2.2)') TRIM(varname),' normalized errors , R', nroot,'B',lev

    WRITE(nout,'(a,a,a,i1,a,i2.2,g14.4)') &
         & 'Relative l1 error   ', TRIM(varname),'   R', nroot,'B',lev, z_err_var%rel_l1
    WRITE(nout,'(a,a,a,i1,a,i2.2,g14.4)') &
         & 'Relative l2 error   ', TRIM(varname),'   R', nroot,'B',lev, z_err_var%rel_l2
    WRITE(nout,'(a,a,a,i1,a,i2.2,g14.4)') &
         & 'Relative linf error ', TRIM(varname),'   R', nroot,'B',lev, z_err_var%rel_linf

    NULLIFY(p_single_patch)
END DO

!---------------------------------------------------------------------
! destruct patches
!---------------------------------------------------------------------
!

CALL destruct_patches(p_patch)

END PROGRAM postpro_gen
