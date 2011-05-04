!>
!! @page pagepostprof901 Main program to calculate normalized errors as in
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
!! Modified by Daniel Reinert, DWD (2010-07-19)
!! - adapted to new type 'external_data'
!!  Modification by Constantin Junk, MPI-M (2011-22-25)
!!  - renamed setup_dynamics and added call of diffusion_nml_setup
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
PROGRAM postpro
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
USE mo_hydro_testcases,      ONLY:  ctest_name, setup_testcase,  &
                                 & rotate_axis_deg
USE mo_icoham_dyn_types,     ONLY: t_hydro_atm, t_hydro_atm_prog
USE mo_namelist
USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, SUCCESS
USE mo_io_units,             ONLY: nnml, nout
USE mo_ocean_nml,            ONLY: setup_ocean_nml
USE mo_dynamics_nml,         ONLY: dynamics_nml_setup, nnow
USE mo_diffusion_nml,        ONLY: diffusion_nml_setup
USE mo_run_nml,              ONLY: ldynamics, ltestcase, latmosphere, locean, &
                                 & lshallow_water, lhydrostatic, run_nml_setup
USE mo_exception,            ONLY: message,finish
USE mo_err_norm,             ONLY: t_errors, err_norm
USE mo_io_vlist,             ONLY: de_reshape1
USE mo_vertical_coord_table, ONLY: init_vertical_coord_table
USE mo_ha_stepping,          ONLY: prepare_ha_integration
USE mo_intp_state,           ONLY: setup_interpol,               &
                                    construct_2d_interpol_state
USE mo_interpolation,        ONLY: t_int_state
USE mo_grf_interpolation,    ONLY: setup_gridref, t_gridref_state, &
                                    construct_2d_gridref_state
USE mo_math_constants,       ONLY: pi
USE mo_physical_constants,   ONLY: re
USE mo_datetime,             ONLY: rdaylen


IMPLICIT NONE

! CHARACTER(len=*), PARAMETER :: version = '$Id$'

! !LOCAL VARIABLES

LOGICAL   :: lgeo, lvort, lwind   ! if .TRUE. errors are calculated for
                                  ! the geopotential, vorticity and wind field

INTEGER   :: nwfiles              ! number of files with the reference wind field
                                  !  if 1 : one file with the u and v components
                                  !   of the reference wind
                                  !  if 2: one file for the reference u component
                                  !   and one for the v comp.

CHARACTER (MAX_CHAR_LENGTH) :: geofile, ref_geofile   !files with the model and
                                                      ! reference geopotential field
CHARACTER (MAX_CHAR_LENGTH) :: ufile, vfile           !files with the model u and v
                                                      ! components
CHARACTER (MAX_CHAR_LENGTH) :: ref_wfile, ref_ufile, ref_vfile  !files with the reference
                                                                ! u and v, u and v components
CHARACTER (MAX_CHAR_LENGTH) :: vortfile, ref_vortfile           !files with the model and
                                                                ! reference vorticity field


TYPE(t_patch), TARGET, ALLOCATABLE        :: p_patch(:)           ! patch information

TYPE(t_patch), POINTER                    :: p_single_patch    => NULL()
TYPE(t_hydro_atm), TARGET, ALLOCATABLE  :: p_hydro_state(:)
TYPE(t_hydro_atm_prog), POINTER         :: p_prog
TYPE(t_int_state), TARGET, ALLOCATABLE    :: p_int_state(:)
TYPE(t_gridref_state), TARGET, ALLOCATABLE:: p_grf_state(:)
TYPE(t_external_data), ALLOCATABLE        :: ext_data(:) !< external data


INTEGER                :: istatus, ic, jg, iv, lev
INTEGER                :: ncells, nverts, nblks_c, nblks_v
INTEGER                :: ist

REAL(wp)               :: z_lon, z_lat, z_lon_ref, z_lat_ref

REAL(wp)               :: z_h_diff, z_u_diff, z_v_diff, z_vort_diff
REAL(wp)               :: u0, u0dre, z_aleph

REAL(wp), ALLOCATABLE  :: z_h(:), z_h_ref(:),                     &
                        & z_u(:), z_v(:), z_u_ref(:), z_v_ref(:), &
                        & z_vort(:), z_vort_ref(:)
REAL(wp), ALLOCATABLE  :: parea(:), darea(:), topo(:)
REAL(wp), ALLOCATABLE  :: lon(:), lat(:)

TYPE(t_errors)           :: z_err_h, z_err_vort, z_err_wind

!

! Namelist groups to be read from NAMELIST_POSTPRO

NAMELIST /setvar/ lgeo, lvort, lwind, nwfiles
NAMELIST /geofiles/ geofile, ref_geofile
NAMELIST /windfiles/ ufile, vfile , ref_wfile, ref_ufile, ref_vfile
NAMELIST /vortfiles/ vortfile, ref_vortfile

!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! set default values
!--------------------------------------------------------------------

lgeo        = .FALSE.
lvort       = .FALSE.
lwind       = .FALSE.
nwfiles     = 1
ref_geofile = "analytical.txt"
ref_wfile   = "analytical.txt"
ref_ufile   = "analytical.txt"
ref_vfile   = "analytical.txt"
ref_vortfile= "analytical.txt"

!-----------------------------------------------------------------------
! open namelist file NAMELIST_ICON
!-----------------------------------------------------------------------
CALL open_nml( "NAMELIST_ICON")

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
!Here is suposed that lshallowater=.TRUE.
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
! open NAMELIST_POSTPRO
!-----------------------------------------------------------------------

CALL open_nml( "NAMELIST_POSTPRO")


!-----------------------------------------------------------------------
! read namelist group setvar
!-----------------------------------------------------------------------
CALL position_nml ('setvar', status=istatus)
SELECT CASE (istatus)
CASE (POSITIONED)
  READ(nnml, setvar)
END SELECT
!WRITE(6,*) lgeo, lwind, lvort, nwfiles


!-----------------------------------------------------------------------
! read namelist group geofiles
!-----------------------------------------------------------------------
IF (lgeo) THEN
  CALL position_nml ('geofiles', status=istatus)
  SELECT CASE (istatus)
  CASE (POSITIONED)
    READ(nnml, geofiles)
  END SELECT
  geofile=TRIM(ADJUSTL(geofile))
  ref_geofile=TRIM(ADJUSTL(ref_geofile))
  !WRITE(6,*) geofile, ref_geofile
ENDIF


!-----------------------------------------------------------------------
! read namelist group vortfiles
!-----------------------------------------------------------------------
IF (lvort) THEN
  CALL position_nml ('vortfiles', status=istatus)
  SELECT CASE (istatus)
  CASE (POSITIONED)
    READ(nnml, vortfiles)
  END SELECT
  vortfile=TRIM(ADJUSTL(vortfile))
  ref_vortfile=TRIM(ADJUSTL(ref_vortfile))
  !WRITE(6,*) vortfile, ref_vortfile
ENDIF

!-----------------------------------------------------------------------
! read namelist group windfiles
!-----------------------------------------------------------------------
IF (lwind) THEN
  CALL position_nml ('windfiles', status=istatus)
  SELECT CASE (istatus)
  CASE (POSITIONED)
    READ(nnml, windfiles)
  END SELECT
  ufile=TRIM(ADJUSTL(ufile))
  vfile=TRIM(ADJUSTL(vfile))
ENDIF

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


  ALLOCATE (z_h(ncells), z_h_ref(ncells))
  ALLOCATE (z_vort(ncells), z_vort_ref(ncells))
  ALLOCATE (z_u(ncells), z_v(ncells))
  ALLOCATE (z_u_ref(ncells), z_v_ref(ncells))

 IF (ctest_name == 'Will_5') THEN
     ALLOCATE (topo(ncells))
     CALL de_reshape1( ext_data(jg)%atm%topography_c, topo )
 ENDIF
 IF (ctest_name == 'Will_2') THEN
      p_prog => p_hydro_state(jg)%prog(nnow(jg))
    ! ! relative vorticity on main grid
    ! Almut's comment: This is now => p_diag%rel_vort_c
    ! ! get divergence
    ! CALL div (p_prog%vn, p_single_patch, p_diag%div)

    ! take the analytical h

     CALL de_reshape1(p_prog%pres_sfc, z_h_ref)

    ! calculate the analytical u, v and vort

    ALLOCATE (lat(ncells))
    ALLOCATE (lon(ncells))
    CALL de_reshape1( p_single_patch%cells%center(:,:)%lat, lat)
    CALL de_reshape1( p_single_patch%cells%center(:,:)%lon, lon)
    u0    = (2.0_wp*pi*re)/(12.0_wp*rdaylen) ! [m/s]
    u0dre = (2.0_wp*pi)   /(12.0_wp*rdaylen)
    z_aleph = rotate_axis_deg * pi / 180.0_wp
    DO ic=1,ncells
     z_lon = lon(ic)
     z_lat = lat(ic)
     z_vort_ref(ic) =2.0_wp*u0dre*(SIN(z_lat)*COS(z_aleph)&
                                   -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
     z_u_ref(ic)=u0*(COS(z_lat)*COS(z_aleph)+COS(z_lon)*SIN(z_lat)*SIN(z_aleph))
     z_v_ref(ic)=-u0*SIN(z_lon)*SIN(z_aleph)
    END DO

 ENDIF



!---------------------------------------------------------------------
! errors for the geopotential
!---------------------------------------------------------------------
  IF (lgeo) THEN


    !open ICOSWM and STSWM files
    OPEN(unit=11, FILE=TRIM(geofile), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file geofile failed')
      CALL finish  ( 'postpro', 'opening file geofile failed')
    END IF

    IF (ctest_name /= 'Will_2') THEN
     OPEN(unit=12, FILE=TRIM(ref_geofile), status="UNKNOWN", IOSTAT=istatus)
     IF (istatus /=SUCCESS) THEN
       CALL message ( 'postpro', 'opening file ref_geofile failed')
       CALL finish  ( 'postpro', 'opening file ref_geofile failed')
     END IF
    END IF
    OPEN(unit=13, FILe="h_diff.gmt", status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file h_diff.gmt failed')
      CALL finish  ( 'postpro', 'opening file h_diff.gmt failed')
    END IF



    DO ic=1,ncells

 !     READ(11,'(2f24.14,g24.14)') z_lon, z_lat, z_h(ic)
 !     READ(12,'(2f24.14,g24.14)') z_lon_ref, z_lat_ref, z_h_ref(ic)
      READ(11,*) z_lon, z_lat, z_h(ic)
      IF (ctest_name == 'Will_5') THEN
       z_h(ic)=z_h(ic)+topo(ic)
      ENDIF
      IF (ctest_name /= 'Will_2') THEN
       READ(12,*) z_lon_ref, z_lat_ref, z_h_ref(ic)
      END IF
      z_h_diff= z_h(ic)- z_h_ref(ic)
      WRITE(13,'(2f24.14,g24.14)') z_lon, z_lat, z_h_diff

    ENDDO


    CLOSE(13)

    CLOSE(unit=11)
    IF (ctest_name /= 'Will_2') THEN
     CLOSE(unit=12)
    END IF
    CALL err_norm( ncells, z_h_ref, z_h, parea, z_err_h)

    WRITE(nout,'(a,i1,a,i2.2)') 'Geopotential normalized errors , R', nroot,'B',lev

    WRITE(nout,*) 'Relative l1 error h   ', z_err_h%rel_l1
    WRITE(nout,*) 'Relative l2 error h   ', z_err_h%rel_l2
    WRITE(nout,*) 'Relative linf error h ', z_err_h%rel_linf

  ENDIF

!---------------------------------------------------------------------
! errors for the vorticity
! now vorticity output is at cells!!!!!
!---------------------------------------------------------------------
  IF (lvort) THEN

    !open ICOSWM and STSWM files

    OPEN(unit=11, FILE=TRIM(vortfile), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file vortfile failed')
      !CALL finish  ( 'postpro', 'opening file vortfile failed')
    END IF
   IF (ctest_name /= 'Will_2') THEN
    OPEN(unit=12, FILE=TRIM(ref_vortfile), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file ref_vortfile failed')
      !CALL finish  ( 'postpro', 'opening file ref_vortfile failed')
    END IF
   END IF
    OPEN(unit=13, FILe="vor_diff.gmt", status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file vort_diff.gmt failed')
      !CALL finish  ( 'postpro', 'opening file vort_diff.gmt failed')
    END IF

    DO iv=1,ncells

 !     READ(11,'(2f24.14,g24.14)') z_lon, z_lat, z_vort(iv)
 !     READ(12,'(2f24.14,g24.14)') z_lon_ref, z_lat_ref, z_vort_ref(iv)
      READ(11,*) z_lon, z_lat, z_vort(iv)
      IF (ctest_name /= 'Will_2') THEN
       READ(12,*) z_lon_ref, z_lat_ref, z_vort_ref(iv)
      END IF
      z_vort_diff= z_vort(iv)- z_vort_ref(iv)
      WRITE(13,'(2f24.14,g24.14)') z_lon, z_lat, z_vort_diff
    ENDDO

    CLOSE(13)

    CLOSE(unit=11)
   IF (ctest_name /= 'Will_2') THEN
    CLOSE(unit=12)
   ENDIF

    CALL err_norm(ncells, z_vort_ref, z_vort, parea, z_err_vort)

    WRITE(nout,'(a,i1,a,i2.2)') 'Vorticity normalized errors , R', nroot,'B',lev

    WRITE(nout,*) 'Relative l1 error vort   ', z_err_vort%rel_l1
    WRITE(nout,*) 'Relative l2 error vort   ', z_err_vort%rel_l2
    WRITE(nout,*) 'Relative linf error vort ', z_err_vort%rel_linf

  ENDIF

!---------------------------------------------------------------------
! errors for the wind
!---------------------------------------------------------------------
  IF (lwind) THEN

    !open ICOSWM and STSWM files

    OPEN(unit=11, FILE=TRIM(ufile), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file ufile failed')
      !CALL finish  ( 'postpro', 'opening file ufile failed')
    END IF

    OPEN(unit=21, FILE=TRIM(vfile), status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file vfile failed')
      !CALL finish  ( 'postpro', 'opening file vfile failed')
    END IF

   IF (ctest_name /= 'Will_2') THEN
    IF (nwfiles == 1) THEN

      OPEN(unit=12, FILE=TRIM(ref_wfile), status="UNKNOWN", IOSTAT=istatus)
      IF (istatus /=SUCCESS) THEN
        CALL message ( 'postpro', 'opening file ref_wfile failed')
       ! CALL finish  ( 'postpro', 'opening file ref_wfile failed')
      END IF


    ELSEIF (nwfiles == 2) THEN
      WRITE(6,*) nwfiles, TRIM(ref_ufile), TRIM(ref_vfile)
      OPEN(unit=12, FILE=TRIM(ref_ufile), status="UNKNOWN", IOSTAT=istatus)
      IF (istatus /=SUCCESS) THEN
        CALL message ( 'postpro', 'opening file ref_ufile failed')
       ! CALL finish  ( 'postpro', 'opening file ref_ufile failed')
      END IF
      OPEN(unit=22, FILE=TRIM(ref_vfile), status="UNKNOWN", IOSTAT=istatus)
      IF (istatus /=SUCCESS) THEN
        CALL message ( 'postpro', 'opening file ref_vfile failed')
        !CALL finish  ( 'postpro', 'opening file ref_vfile failed')
      END IF

    ENDIF
    ENDIF

    OPEN(unit=13, FILe="u_diff.gmt", status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file u_diff.gmt failed')
      !CALL finish  ( 'postpro', 'opening file u_diff.gmt failed')
    END IF

    OPEN(unit=23, FILe="v_diff.gmt", status="UNKNOWN", IOSTAT=istatus)
    IF (istatus /=SUCCESS) THEN
      CALL message ( 'postpro', 'opening file v_diff.gmt failed')
      !CALL finish  ( 'postpro', 'opening file v_diff.gmt failed')
    END IF

    DO ic=1,ncells

!      READ(11,'(2f24.14,g24.14)') z_lon, z_lat, z_u(ic)
!      READ(21,'(2f24.14,g24.14)') z_lon, z_lat, z_v(ic)
      READ(11,*) z_lon, z_lat, z_u(ic)
      READ(21,*) z_lon, z_lat, z_v(ic)
      IF (ctest_name /= 'Will_2') THEN
       IF (nwfiles == 1) THEN
 !        READ(12,'(2f24.14,2g24.14)') z_lon_ref, z_lat_ref, z_u_ref(ic), &
 !         &                          z_v_ref(ic)
         READ(12,*) z_lon_ref, z_lat_ref, z_u_ref(ic), &
           &                          z_v_ref(ic)
       ELSEIF (nwfiles == 2) THEN
 !       READ(12,'(2f24.14,g24.14)') z_lon_ref, z_lat_ref, z_u_ref(ic)
 !       READ(22,'(2f24.14,g24.14)') z_lon_ref, z_lat_ref, z_v_ref(ic)
         READ(12,*) z_lon_ref, z_lat_ref, z_u_ref(ic)
         READ(22,*) z_lon_ref, z_lat_ref, z_v_ref(ic)
       ENDIF
      ENDIF
      z_u_diff= z_u(ic)- z_u_ref(ic)
      WRITE(13,'(2f24.14,g24.14)') z_lon, z_lat, z_u_diff
      z_v_diff= z_v(ic)- z_v_ref(ic)
      WRITE(23,'(2f24.14,g24.14)') z_lon, z_lat, z_v_diff


    ENDDO

    CLOSE(13)
    CLOSE(23)

    CLOSE(unit=11)
   IF (ctest_name /= 'Will_2') THEN
    IF (nwfiles == 1) THEN
      CLOSE(unit=12)
    ELSEIF (nwfiles == 2) THEN
      CLOSE(unit=12)
      CLOSE(unit=22)
    ENDIF
    ENDIF
    CALL err_norm( ncells, z_u_ref, z_v_ref, z_u, z_v, parea,         &
      &            z_err_wind )

    WRITE(nout,'(a,i1,a,i2.2)') 'Wind normalized errors ,  R', nroot,'B',lev

    WRITE(nout,*) 'Relative l1 error w   ', z_err_wind%rel_l1
    WRITE(nout,*) 'Relative l2 error w   ', z_err_wind%rel_l2
    WRITE(nout,*) 'Relative linf error w ', z_err_wind%rel_linf

  ENDIF
  NULLIFY(p_single_patch)

END DO

!---------------------------------------------------------------------
! destruct patches
!---------------------------------------------------------------------
!

CALL destruct_patches(p_patch)

END PROGRAM postpro
