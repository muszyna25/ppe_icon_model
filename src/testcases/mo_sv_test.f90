!>
!! Module for deformational flow test (stationary vortex)
!!
!! Test case for horizontal advection problems. The deformation
!! test consists of two static vortices located on diametrically opposite
!! sides of the sphere. Because of the evolving sharp gradients, this test
!! is more challenging than the solid body rotation test.
!!
!! Literature:
!! - Nair and Machenauer (2002): The mass-conservative cell-integrated semi-
!! Lagrangian advection scheme on the sphere. Mon. Wea. Rev., 130, 649-667
!! - Nair and Jablonowski (2008): Moving vortices on the sphere: A test
!! case for horizontal advection problems. Mon. Wea. Rev., 136, 699-711
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial version by Daniel Reinert (2009-10-06)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_sv_test

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_physical_constants,  ONLY: earth_radious
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_math_constants,      ONLY: dbl_eps, pi, pi_4
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_eta_coord_diag,      ONLY: half_level_pressure, full_level_pressure
  USE mo_parallel_config,  ONLY: nproma
  USE mo_exception,           ONLY: finish

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

 !DEFINED PARAMETERS:
  REAL(wp), PARAMETER :: zp0     = 100000._wp !< surface pressure
  REAL(wp), PARAMETER :: zt0     = 300._wp    !< temperature (isothermal)

 !VORTEX PARAMETERS: (see Nair and Machenauer (2002))
  REAL(wp), PARAMETER :: rho0    = 3._wp !< constant, regarding radial
                                         !< distance from vortex center
  REAL(wp), PARAMETER :: gamma   = 5._wp
  REAL(wp), PARAMETER :: tottime = 1036800._wp  !< total time 12 days

  ! (npole_lon,npole_lat) new north pole position of the rotated coordinates
  ! relative to the given regular coordinates
  ! Nair et al. (2005), Putman and Lin (2007) setting
  REAL(wp), PARAMETER :: npole_lon = pi - 0.8_wp + pi_4
  REAL(wp), PARAMETER :: npole_lat = pi/4.8_wp

  PUBLIC :: init_hydro_state_prog_svtest, init_sv_wind, get_sv_tracer, &
    &       rotated_sphere


  CONTAINS


  !--------------------------------------------------------------------
  !>
  !! SUBROUTINE init_hydro_state_prog_svtest
  !!
  !! Short description:
  !! Initialization of prognostic state vector for stationary vortex
  !! test case.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert (2009-10-05)
  !!
  SUBROUTINE init_hydro_state_prog_svtest( ptr_patch, ptr_prog, ptr_diag, &
    &                                      ptr_int, ptr_ext_data )

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,INTENT(INOUT) :: ptr_patch
    TYPE(t_int_state),  INTENT(INOUT)  :: ptr_int
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: ptr_prog
    TYPE(t_hydro_atm_diag), INTENT(INOUT) :: ptr_diag
    TYPE(t_external_data), INTENT(INOUT)    :: ptr_ext_data !< external data

    INTEGER  :: jk,jb   !< loop indices
    INTEGER  :: ist     !< status variable
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlev, nlevp1          !< number of full and half levels
    INTEGER  :: ikp1
    INTEGER  :: nlen

    REAL(wp), ALLOCATABLE :: zhelp_c (:,:,:)
  !--------------------------------------------------------------------

    ! allocate memory for temporary field
    nblks_c = ptr_patch%nblks_c
    nlev    = ptr_patch%nlev
    nlevp1  = ptr_patch%nlevp1
    ALLOCATE( zhelp_c(nproma,nlev,nblks_c), STAT=ist)
    IF(ist/=SUCCESS)THEN
       CALL finish('mo_pa_test:init_hydro_state_prog_patest',  &
         &         'allocation of zhelp_c failed')
    ENDIF

    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c


    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_prog%pres_sfc(:,:) = zp0

    ! init temperature (isothermal atmosphere)
    ptr_prog%temp(:,:,:)   = zt0

    ! init vertical velocity field
    ptr_diag%weta(:,:,:)      = 0._wp
    ptr_diag%wpres_mc(:,:,:)  = 0._wp

    ! init surface geopotential
    ptr_diag%geo_ic(:,nlevp1,:) = 0._wp

    ! init Coriolis parameter
    ptr_patch%cells%f_c(:,:) = 0._wp


    ! compute pressure value at half- and full levels as
    ! well as layer thickness at cell- and edge centers.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen,ikp1)
    DO jb = 1, nblks_c
       IF (jb /= nblks_c) THEN
          nlen = nproma
       ELSE
          nlen = npromz_c
       ENDIF
       ! compute half-level pressure
       CALL half_level_pressure( ptr_prog%pres_sfc(:,jb), nproma, nlen, &! in
                                 ptr_diag%pres_ic(:,:,jb)              ) ! out

       ! compute pressure value at full levels
       CALL full_level_pressure( ptr_diag%pres_ic(:,:,jb), nproma, nlen, &! in
                                 ptr_diag%pres_mc(:,:,jb)               ) ! out

       DO jk = 1, nlev
          ikp1 = jk+1
          ! compute layer thickness at cell centers
          zhelp_c(1:nlen,jk,jb) = ptr_diag%pres_ic(1:nlen,ikp1,jb) &
                                - ptr_diag%pres_ic(1:nlen,jk,jb)
          ptr_diag%rdelp_c(1:nlen,jk,jb) = 1.0_wp/zhelp_c(1:nlen,jk,jb)
       ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! compute layer thickness at edge centers
    CALL cells2edges_scalar( zhelp_c, ptr_patch, ptr_int%c_lin_e, &
                             ptr_diag%delp_e )



    ! Initialize wind field
    CALL init_sv_wind ( ptr_patch, ptr_prog )


    ! Initialize tracer field 1 for time t=0
    CALL get_sv_tracer( ptr_patch, ptr_prog%tracer(:,:,:,1), 0._wp )

    ! Initialize tracer field 2 accordingly
    ptr_prog%tracer(:,:,:,2) = ptr_prog%tracer(:,:,:,1)

  END SUBROUTINE init_hydro_state_prog_svtest


  !--------------------------------------------------------------------
  !>
  !! SUBROUTINE init_sv_wind
  !!
  !! Short description:
  !! Initialization of horizontal velocity field for static vortex
  !! test case.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert (2009-10-05)
  !!
  SUBROUTINE init_sv_wind ( ptr_patch, ptr_prog )

    !INPUT PARAMETERS:
    TYPE(t_patch), INTENT(IN)         :: ptr_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: ptr_prog

    INTEGER  :: jb,je,jk        !< loop indices
    INTEGER  :: nblks_e, npromz_e, nlen
    INTEGER  :: nlev            !< number of full levels
    REAL(wp) :: u_wind, v_wind  !< zonal and meridional velocity component
    REAL(wp) :: zomega          !< angular velocity in rad s^-1
    REAL(wp) :: zrho
    REAL(wp) :: zlon, zlat          !< lon/lat in unrotated system
    REAL(wp) :: zlon_rot, zlat_rot  !< lon/lat in rotated system
    REAL(wp) :: zt
    REAL(wp) :: vmax

  !--------------------------------------------------------------------

    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    nlev = ptr_patch%nlev

    ! Since this is a 2D test case, the velocity field
    ! only needs to be calculated for 1 full level. Then the
    ! solution is copied to each vertical level.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nlen,zlon,zlat,zlon_rot,zlat_rot,vmax,  &
!$OMP            zrho,zt,zomega,u_wind,v_wind)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF

      DO je = 1, nlen
        ! location of edge midpoint
        zlon = ptr_patch%edges%center(je,jb)%lon
        zlat = ptr_patch%edges%center(je,jb)%lat

        ! get zlat, zlon in rotated system with npole at
        ! (npole_lon,npole_lat)
        CALL rotated_sphere( npole_lon,npole_lat,  & !<in
         &                   zlon,zlat,            & !<in
         &                   zlon_rot,zlat_rot )     !<inout

        ! get angular velocity depending on zlat_rot
        vmax    = (2._wp*pi*earth_radious/tottime)*1.5_wp*sqrt(3._wp)
        zrho   = rho0 * cos(zlat_rot)
        zt     = tanh(zrho)
        zomega = vmax * (1._wp - zt**2._wp) * zt

        IF (ABS(zrho) < dbl_eps) THEN
          zomega = 0._wp
        ELSE
          zomega = zomega / (earth_radious * zrho)
        ENDIF

        ! calculate zonal and meridional velocity component at edge midpoint
        u_wind = earth_radious * zomega * ( (sin(npole_lat) * cos(zlat))            &
          & - (cos(npole_lat) * cos(zlon - npole_lon) * sin(zlat)) )

        v_wind = earth_radious * zomega * ( cos(npole_lat) * sin(zlon - npole_lon) )

        ! calculate normal wind component
        ptr_prog%vn(je,1,jb) = &
                  u_wind * ptr_patch%edges%primal_normal(je,jb)%v1 &
          &     + v_wind * ptr_patch%edges%primal_normal(je,jb)%v2

      END DO  ! cell loop

    END DO  ! block loop
!$OMP END DO
!$OMP END PARALLEL

    DO jk = 2, nlev
      ptr_prog%vn(:,jk,:) =  ptr_prog%vn(:,jk-1,:)
    END DO

  END SUBROUTINE init_sv_wind


  !--------------------------------------------------------------------
  !!>
  !! SUBROUTINE get_sv_tracer
  !!
  !! Short description:
  !! Analytic tracer field at any given time. This subroutine will
  !! be called for initialization as well as validation purposes
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD (2009-10-06)
  !!
  SUBROUTINE get_sv_tracer( ptr_patch, p_cc, p_sim_time )

    !INPUT PARAMETERS:
    TYPE(t_patch), INTENT(IN) :: ptr_patch

    REAL(wp), INTENT(IN)    :: p_sim_time  !< simulation time in seconds
                                         !< since start
    REAL(wp), INTENT(INOUT) :: p_cc(:,:,:) !< tracer array

    INTEGER   :: jb,jc,jk            !> loop indices
    INTEGER   :: nblks_c, npromz_c, nlen
    INTEGER   :: nlev                !< number of full levels
    REAL(wp)  :: zlon, zlat          !< lon/lat of cell centers (unrotated)
    REAL(wp)  :: zlon_rot, zlat_rot  !< lon/lat in rotated system
    REAL(wp)  :: zinner
    REAL(wp)  :: ztime
    REAL(wp)  :: zomega     !< angular velocity in rad s^-1
    REAL(wp)  :: zrho
    REAL(wp)  :: zt
    REAL(wp)  :: vmax
  !--------------------------------------------------------------------

    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c

    nlev = ptr_patch%nlev

    ! Since this is a 2D test case, the analytic tracer distribution
    ! only needs to be calculated for 1 full level. Then the solution
    ! is copied to each vertical level.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nlen,zlon,zlat,zlon_rot,zlat_rot,vmax,  &
!$OMP            zrho,zt,zomega,ztime,zinner)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jc=1, nlen
        ! location of cell center
        zlon = ptr_patch%cells%center(jc,jb)%lon
        zlat = ptr_patch%cells%center(jc,jb)%lat

        ! get zlat, zlon in rotated system
        CALL rotated_sphere( npole_lon,npole_lat,  & !<in
         &                   zlon,zlat,            & !<in
         &                   zlon_rot,zlat_rot )     !<inout

        ! get angular velocity depending on zlat_rot
        vmax   = (2._wp*pi*earth_radious/tottime)*1.5_wp*sqrt(3._wp)
        zrho   = rho0 * cos(zlat_rot)
        zt     = tanh(zrho)
        zomega = vmax * (1._wp - zt**2._wp) * zt

        IF (ABS(zrho) < dbl_eps) THEN
          zomega = 0._wp
        ELSE
          zomega = zomega / (earth_radious * zrho)
        ENDIF

        ztime   = p_sim_time !!* v0 ??? what is 'scale' ?
        zinner = (zrho/gamma) * sin(zlon_rot - zomega * ztime)

        ! exact solution
        p_cc(jc,1,jb) = 1._wp - tanh(zinner)

      ENDDO  ! cell loop

    ENDDO  ! block loop
!$OMP END DO
!$OMP END PARALLEL

     DO jk = 2, nlev
         p_cc(:,jk,:) =  p_cc(:,jk-1,:)
     ENDDO

  END SUBROUTINE get_sv_tracer


  !--------------------------------------------------------------------
  !!>
  !! SUBROUTINE rotated_sphere
  !!
  !! Short description:
  !! Rotate a sphere so that its NP is at (lc,tc). The rotated coordinates
  !!(rol,rot) corresponding to (la,th) in the regular sph. coordinates
  !!
  !! @par Revision History
  !! Initial version by Ram Nair, NCAR/SCD (2006-08-01)
  !! Adapted to ICON by Daniel Reinert, DWD (2009-10-05)
  !!
  SUBROUTINE rotated_sphere(lc,tc,la,th,rol,rot)

    IMPLICIT NONE
    REAL(wp), INTENT(IN)  :: lc,tc,la,th
    REAL(wp), INTENT(OUT) :: rol,rot

    REAL(wp) :: cost,sint,sinc,cosc
    REAL(wp) :: trm,trm1,trm2,trm3
  !--------------------------------------------------------------------

    sinc = sin(tc)
    cosc = cos(tc)
    cost = cos(th)
    sint = sin(th)

    trm  = cost * cos(la - lc)
    trm1 = cost * sin(la - lc)
    trm2 = sinc * trm  - cosc * sint
    trm3 = sinc * sint + cosc * trm

    rol = atan2(trm1,trm2)
    IF ( rol > 2._wp*pi  ) rol = rol - 2._wp*pi
    IF ( rol < 0._wp     ) rol = rol + 2._wp*pi
    rot = asin(trm3)

  END SUBROUTINE rotated_sphere


END MODULE mo_sv_test

