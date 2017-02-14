!>
!! Initializes shallow water tests for the lshallow_water option in the.
!!
!! Initializes shallow water tests for the lshallow_water option in the
!! hydrostatic code.
!! NOTE: The shallow water model assumes that
!!       "pres_sfc" contains the "thickness" (not the height),
!!       the horizontal normal velocities are as usual,
!!       temperature is not used. For linearization, the variable
!!       sw_ref_height is set here.
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2008-10)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_sw_test
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog
  USE mo_physical_constants,  ONLY: rgrav, rdaylen
  USE mo_math_constants,      ONLY: pi, pi_2
  USE mo_dynamics_config,     ONLY: sw_ref_height
  USE mo_parallel_config,     ONLY: nproma
  USE mo_grid_config,         ONLY: grid_sphere_radius, grid_angular_velocity

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_will2_test, init_will3_test, init_will5_test,&
            init_will6_test, init_usbr_test,  init_swgw_test

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialization of the Williamson (1992) testcase 2.
  !!
  !! Adapted from original shallow water code.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_will2_test (pt_patch, pt_prog, pt_ext_data, p_rotate_axis_deg)

    IMPLICIT NONE

    !INPUT PARAMETER:
    REAL(wp)                      :: p_rotate_axis_deg

    !INPUT/OUTPUT PARAMETER:
    TYPE(t_patch), INTENT(INOUT)    :: pt_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
    TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

    REAL(wp), PARAMETER  :: h0 = 2.94e4_wp * rgrav  ! maximum height
    REAL(wp)  :: u0     ! [m/s]

    REAL(wp)  :: z_fact1, z_fact2, z_lat, z_lon, z_aleph, z_uu, z_vv
    INTEGER   :: jb, jc, jv, je, nblks_e, nblks_c, nblks_v, &
                 npromz_e, npromz_c, npromz_v, nlen

!-----------------------------------------------------------------------
    u0 =(2.0_wp*pi*grid_sphere_radius)/(12.0_wp*rdaylen) ! [m/s]

    ! set reference height
    sw_ref_height = 0.9_wp * h0

    ! rotation in radiant
    z_aleph = p_rotate_axis_deg * pi / 180.0_wp

    ! 1st factor for height (= thickness, because topography is zero)
    z_fact1 = grid_sphere_radius * grid_angular_velocity
    z_fact1 = z_fact1 + 0.5_wp * u0
    z_fact1 = z_fact1 * u0 * rgrav

    ! set topography
    pt_ext_data%atm%topography_c(:,:) = 0.0_wp

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c
    nblks_e   = pt_patch%nblks_e
    npromz_e  = pt_patch%npromz_e
    nblks_v   = pt_patch%nblks_v
    npromz_v  = pt_patch%npromz_v

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1, nlen
        z_lat   = pt_patch%cells%center(jc,jb)%lat
        z_lon   = pt_patch%cells%center(jc,jb)%lon

        z_fact2 = SIN(z_lat) * COS(z_aleph)
        z_fact2 = z_fact2 - COS(z_lon) * COS(z_lat) * SIN(z_aleph)
        z_fact2 = z_fact2 * z_fact2

        ! height
        pt_prog%pres_sfc(jc,jb) = h0 - z_fact1 * z_fact2

        ! Coriolis parameter
        pt_patch%cells%f_c(jc,jb) = 2.0_wp*grid_angular_velocity*(SIN(z_lat)*COS(z_aleph)&
                                  -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
      ENDDO
    ENDDO

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO je = 1, nlen
        z_lat   = pt_patch%edges%center(je,jb)%lat
        z_lon   = pt_patch%edges%center(je,jb)%lon

        z_uu = COS(z_lat) * COS(z_aleph)
        z_uu = z_uu + COS(z_lon) * SIN(z_lat) * SIN(z_aleph)
        z_uu = u0 * z_uu

        z_vv = SIN(z_lon) * SIN(z_aleph)
        z_vv = -1._wp * u0 * z_vv

        ! normal velocity component
        pt_prog%vn(je,1,jb) = &
            z_uu * pt_patch%edges%primal_normal(je,jb)%v1 &
           +z_vv * pt_patch%edges%primal_normal(je,jb)%v2

        ! Coriolis parameter
        pt_patch%edges%f_e(je,jb) = 2.0_wp*grid_angular_velocity*(SIN(z_lat)*COS(z_aleph)&
                                  -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
      ENDDO
    ENDDO

    DO jb = 1, nblks_v
      IF (jb /= nblks_v) THEN
         nlen = nproma
      ELSE
         nlen = npromz_v
      ENDIF
      DO jv = 1, nlen
        z_lat   = pt_patch%verts%vertex(jv,jb)%lat
        z_lon   = pt_patch%verts%vertex(jv,jb)%lon
        ! Coriolis parameter
        pt_patch%verts%f_v(jv,jb) = 2.0_wp*grid_angular_velocity*(SIN(z_lat)*COS(z_aleph)&
                                  -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
      ENDDO
    ENDDO

  END SUBROUTINE init_will2_test


  !-------------------------------------------------------------------------
  !>
  !! Initialization of the Williamson (1992) testcase 3.
  !!
  !! Adapted from original shallow water code.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_will3_test (pt_patch, pt_prog, pt_ext_data, p_rotate_axis_deg)

    IMPLICIT NONE

    !INPUT PARAMETER:
    REAL(wp)                      :: p_rotate_axis_deg

    !INPUT/OUTPUT PARAMETER:
    TYPE(t_patch), INTENT(INOUT)    :: pt_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
    TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

    REAL(wp), PARAMETER  :: h0 = 2.94e4_wp * rgrav  ! maximum height
    REAL(wp)  :: u0

    INTEGER   :: jb, jc, je, jv, nblks_e, nblks_c, nblks_v, &
                 npromz_e, npromz_c, npromz_v, nlen
    REAL(wp)  :: z_aleph, z_lat, z_lon, z_hh, z_rotlon, z_rotlat, z_uu, z_vv
    REAL(wp)  :: inverse_sphere_radius
    
  !-----------------------------------------------------------------------
    ! set reference height
    u0 =(2.0_wp*pi*grid_sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp) ! [m/s]
    inverse_sphere_radius = 1.0_wp / grid_sphere_radius
    sw_ref_height = 0.9_wp * h0

    z_aleph = p_rotate_axis_deg * pi/180.0_wp

    pt_ext_data%atm%topography_c(:,:)=0.0_wp

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c
    nblks_e   = pt_patch%nblks_e
    npromz_e  = pt_patch%npromz_e
    nblks_v   = pt_patch%nblks_v
    npromz_v  = pt_patch%npromz_v

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1, nlen
        z_lat   = pt_patch%cells%center(jc,jb)%lat
        z_lon   = pt_patch%cells%center(jc,jb)%lon

        CALL rotate( z_lon, z_lat, z_aleph, z_rotlon, z_rotlat )
        z_hh = geostr_balance( z_rotlat, symmetric_u_velo, &
          & inverse_sphere_radius, grid_angular_velocity, u0 )
        pt_prog%pres_sfc(jc,jb) = h0 - grid_sphere_radius * rgrav * z_hh

        ! Coriolis parameter
        pt_patch%cells%f_c(jc,jb) = 2.0_wp*grid_angular_velocity*(SIN(z_lat)*COS(z_aleph)&
                                  -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
      ENDDO
    ENDDO

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO je = 1, nlen
        z_lat   = pt_patch%edges%center(je,jb)%lat
        z_lon   = pt_patch%edges%center(je,jb)%lon

        CALL rotate( z_lon, z_lat, z_aleph, z_rotlon, z_rotlat)

        z_uu = COS(z_aleph) * SIN(z_rotlon) * SIN(z_lon)
        z_uu = COS(z_lon)   * COS(z_rotlon) + z_uu
        z_uu = symmetric_u_velo(z_rotlat,u0) * z_uu

        z_vv = COS(z_aleph) * SIN(z_lat) * COS(z_lon)   * SIN(z_rotlon)
        z_vv = z_vv         - SIN(z_lat) * SIN(z_lon)   * COS(z_rotlon)
        z_vv = z_vv         - COS(z_lat) * SIN(z_aleph) * SIN(z_rotlon)
        z_vv = symmetric_u_velo(z_rotlat,u0) * z_vv

        ! normal velocity component
        pt_prog%vn(je,1,jb) = &
            z_uu * pt_patch%edges%primal_normal(je,jb)%v1 &
           +z_vv * pt_patch%edges%primal_normal(je,jb)%v2

        ! Coriolis parameter
        pt_patch%edges%f_e(je,jb) = 2.0_wp*grid_angular_velocity*(SIN(z_lat)*COS(z_aleph)&
                                  -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
      ENDDO
    ENDDO

    DO jb = 1, nblks_v
      IF (jb /= nblks_v) THEN
         nlen = nproma
      ELSE
         nlen = npromz_v
      ENDIF
      DO jv = 1, nlen
        z_lat   = pt_patch%verts%vertex(jv,jb)%lat
        z_lon   = pt_patch%verts%vertex(jv,jb)%lon
        ! Coriolis parameter
        pt_patch%verts%f_v(jv,jb) = 2.0_wp*grid_angular_velocity*(SIN(z_lat)*COS(z_aleph)&
                                  -COS(z_lon)*COS(z_lat)*SIN(z_aleph))
      ENDDO
    ENDDO

  END SUBROUTINE init_will3_test

  !-------------------------------------------------------------------------
  !>
  !! Initialization of the Williamson (1992) testcase 5.
  !!
  !! Adapted from original shallow water code.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_will5_test (pt_patch, pt_prog, pt_ext_data)

    IMPLICIT NONE

    !INPUT/OUTPUT PARAMETER:
    TYPE(t_patch), INTENT(INOUT)    :: pt_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
    TYPE(t_external_data), INTENT(INOUT)    :: pt_ext_data !< external data

    REAL(wp), PARAMETER  :: h0 = 5960.0_wp  ! maximum height [m]
    REAL(wp), PARAMETER  :: u0 = 20.0_wp    ! [m/s]
    REAL(wp), PARAMETER  :: mount_lon = -pi_2
    REAL(wp), PARAMETER  :: mount_lat = pi/6.0_wp
    REAL(wp), PARAMETER  :: mount_rad = pi/9.0_wp
    REAL(wp), PARAMETER  :: mount_hei = 2000.0_wp

    REAL(wp)  :: z_fact1, z_fact2, z_lat, z_lon, z_uu, z_vv, z_diff, &
                 z_dist_mc, z_min_dist_sq
    INTEGER   :: jb, jc, je, nblks_e, nblks_c, &
                 npromz_e, npromz_c, nlen

  !-----------------------------------------------------------------------

    ! set reference height
    sw_ref_height = 0.9_wp * h0

    ! 1st factor for height (= thickness, because topography is zero)
    z_fact1 = grid_sphere_radius * grid_angular_velocity
    z_fact1 = z_fact1 + 0.5_wp * u0
    z_fact1 = z_fact1 * u0 * rgrav

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c
    nblks_e   = pt_patch%nblks_e
    npromz_e  = pt_patch%npromz_e

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1, nlen
        z_lat   = pt_patch%cells%center(jc,jb)%lat
        z_lon   = pt_patch%cells%center(jc,jb)%lon

        ! topography
        z_diff    = z_lon - mount_lon
        z_dist_mc = z_diff * z_diff
        z_diff    = z_lat - mount_lat
        z_dist_mc = z_dist_mc + z_diff * z_diff
        z_diff    = mount_rad * mount_rad
        z_min_dist_sq = MIN ( z_diff, z_dist_mc)
        z_dist_mc = SQRT( z_min_dist_sq)
        pt_ext_data%atm%topography_c(jc,jb) =  mount_hei*&
                               (1.0_wp-z_dist_mc / mount_rad)

        ! second factor for height
        z_fact2 = SIN(z_lat)
        z_fact2 = z_fact2 * z_fact2

        ! thickness= height-surface
        pt_prog%pres_sfc(jc,jb) = h0 - z_fact1 * z_fact2 &
                          - pt_ext_data%atm%topography_c(jc,jb)

      ENDDO
    ENDDO

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO je = 1, nlen
        z_lat   = pt_patch%edges%center(je,jb)%lat
        z_lon   = pt_patch%edges%center(je,jb)%lon

        z_uu = u0 * COS(z_lat)
        z_vv = 0.0_wp

        ! normal velocity component
        pt_prog%vn(je,1,jb) = &
            z_uu * pt_patch%edges%primal_normal(je,jb)%v1 &
           +z_vv * pt_patch%edges%primal_normal(je,jb)%v2

      ENDDO
    ENDDO

  END SUBROUTINE init_will5_test

!-------------------------------------------------------------------------
  !>
  !! Initialization of the Williamson (1992) testcase 6.
  !!
  !! Adapted from original shallow water code.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_will6_test (pt_patch, pt_prog, pt_ext_data)

    IMPLICIT NONE

    !INPUT/OUTPUT PARAMETER:
    TYPE(t_patch), INTENT(INOUT)    :: pt_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
    TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

    REAL(wp), PARAMETER  :: h0     = 8000.0_wp     ! maximum height
    REAL(wp), PARAMETER  :: omg_kk = 7.848e-6_wp  ! grid_angular_velocity = K
    INTEGER,  PARAMETER  :: ir     = 4             ! R

    INTEGER   :: j, jb, jc, je, nblks_e, nblks_c, npromz_e, npromz_c, nlen
    REAL(wp)  :: z_lat, z_lon, z_hh, z_uu, z_vv
    REAL(wp)  :: z_phia, z_phib, z_phic, z_r_omega , z_re_omg_kk
    REAL(wp)  :: z_cosfi, z_cosfi2, z_cosfir, z_cosfir2, z_cosfir2m2, z_cosfirm1
    REAL(wp)  :: z_cosdl, z_cosd2l, z_dlon, z_rr1r2, z_sindl, z_sinfi, z_sinfi2
    REAL(wp)  :: z_val, z_r, z_r1, z_r1r1, z_r1r2, z_r2


  !-----------------------------------------------------------------------
    ! set reference height
    sw_ref_height =  h0

    z_r_omega  = grid_sphere_radius * grid_angular_velocity
    z_re_omg_kk= grid_sphere_radius * omg_kk

    z_r       = REAL(ir,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r1    = z_r1 * z_r1
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    pt_ext_data%atm%topography_c(:,:)=0.0_wp

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c
    nblks_e   = pt_patch%nblks_e
    npromz_e  = pt_patch%npromz_e

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1, nlen
        z_lat   = pt_patch%cells%center(jc,jb)%lat
        z_lon   = pt_patch%cells%center(jc,jb)%lon

        z_dlon    = z_lon * z_r
        z_cosdl   = COS(z_dlon)
        z_cosd2l  = COS(2._wp * z_dlon)

        z_cosfi   = COS(z_lat)
        z_cosfi2  = z_cosfi  * z_cosfi    ! cos^2(lat)

        z_cosfir  = z_cosfi
        DO j= 2, ir-1
          z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
        ENDDO
        z_cosfir2m2 = z_cosfir
        z_cosfir2m2 = z_cosfir2m2 * z_cosfir2m2   ! cos^{2*r1-2}(lat)

        z_cosfir  = z_cosfir * z_cosfi    ! cos^{r1}(lat)
        z_cosfir2 = z_cosfir * z_cosfir   ! cos^{2*r1}(lat)

        z_val  = -.25_wp + z_r
        z_val  = 2._wp * z_val * z_val - 2.125_wp   ! 2r^2 - r -2

        z_phia = z_val * z_cosfi2

        z_val  = 2._wp * REAL(ir,wp) * z_r

        z_phia = z_phia - z_val
        z_val  = z_cosfi2 * z_cosfi2 * z_r1
        z_phia = z_phia + z_val

        z_phia = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2m2 * z_phia
        z_val  = .5_wp * z_re_omg_kk * &
                 (2._wp * z_r_omega + z_re_omg_kk) * z_cosfi2
        z_phia = z_val + z_phia

        z_phib = -1._wp * z_cosfi2 * z_r1r1 + z_r1r1 + 1._wp
        z_phib = z_re_omg_kk * (z_r_omega + z_re_omg_kk) * z_cosfir * z_phib
        z_phib = 2._wp * z_rr1r2 * z_phib

        z_phic = z_r1 * z_cosfi2 - 1._wp * z_r2
        z_phic = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2 * z_phic

        z_hh   = (z_phia + z_phib * z_cosdl + z_phic * z_cosd2l) * rgrav
        z_hh   = h0 + z_hh

        pt_prog%pres_sfc(jc,jb) = z_hh
      ENDDO
    ENDDO

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO je = 1, nlen
        z_lat   = pt_patch%edges%center(je,jb)%lat
        z_lon   = pt_patch%edges%center(je,jb)%lon

        z_dlon    = z_lon * z_r
        z_sinfi   = SIN(z_lat)
        z_cosfi   = COS(z_lat)
        z_cosfir  = z_cosfi
        DO j= 2, ir-1
          z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
        ENDDO
          z_cosfirm1 = z_cosfir

        !u
        z_cosdl   = COS(z_dlon)
        z_sinfi2  = z_sinfi * z_sinfi
        z_cosfi2  = z_cosfi * z_cosfi
        z_val      = z_r * z_sinfi2 - z_cosfi2
        z_val      = z_cosfirm1 * z_val * z_cosdl
        z_val      = z_cosfi + z_val
        z_uu       = z_re_omg_kk * z_val

        !v
        z_sindl   = SIN(z_dlon)
        z_val      = z_cosfirm1 * z_sinfi * z_sindl
        z_vv       = -1._wp * z_re_omg_kk * z_r * z_val

        ! normal velocity component
        pt_prog%vn(je,1,jb) = &
            z_uu * pt_patch%edges%primal_normal(je,jb)%v1 &
           +z_vv * pt_patch%edges%primal_normal(je,jb)%v2

      ENDDO
    ENDDO

  END SUBROUTINE init_will6_test

  !-------------------------------------------------------------------------
  !>
  !! Initialization of a simple gravity test.
  !!
  !! Adapted from original shallow water code.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_swgw_test (pt_patch, pt_prog, pt_ext_data)

    IMPLICIT NONE

    !INPUT/OUTPUT PARAMETER:
    TYPE(t_patch), INTENT(INOUT)    :: pt_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
    TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

    REAL(wp), PARAMETER  :: h0     = 1000.0_wp     ! maximum height
    REAL(wp), PARAMETER  :: &
      c_lon = 0.0_wp,       & ! degrees
      c_lat = 0.0_wp,       &
      r_lon = 1e6_wp,       & ! meters
      r_lat = 1e6_wp,       &
      amplitude = 1.0_wp


    INTEGER   :: jb, jc, je, nblks_e, nblks_c, npromz_e, npromz_c, nlen
    REAL(wp)  :: z_lat, z_lon, z_hh

  !-----------------------------------------------------------------------
    ! set reference height
    sw_ref_height =  h0

    pt_ext_data%atm%topography_c(:,:)=0.0_wp

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c
    nblks_e   = pt_patch%nblks_e
    npromz_e  = pt_patch%npromz_e

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1, nlen
        z_lat   = pt_patch%cells%center(jc,jb)%lat
        z_lon   = pt_patch%cells%center(jc,jb)%lon

        z_hh  = h0 + amplitude*&
                EXP(-( (grid_sphere_radius*(z_lon-pi*c_lon/180.0_wp)/r_lon)**2   &
                     + (grid_sphere_radius*(z_lat-pi*c_lat/180.0_wp)/r_lat)**2 ) )

        pt_prog%pres_sfc(jc,jb) = z_hh
      ENDDO
    ENDDO

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO je = 1, nlen

        pt_prog%vn(je,1,jb) = 0.0_wp

      ENDDO
    ENDDO

  END SUBROUTINE init_swgw_test

  !-------------------------------------------------------------------------
  !>
  !! Initialization of the unsteady solid body rotation test.
  !!
  !! Adapted from original shallow water code.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2008-10)
  !!
  SUBROUTINE init_usbr_test (pt_patch, pt_prog, pt_ext_data)

    IMPLICIT NONE

    !INPUT/OUTPUT PARAMETER:
    TYPE(t_patch), INTENT(INOUT)    :: pt_patch
    TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
    TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

    REAL(wp) :: u0  ! [m/s]
    REAL(wp), PARAMETER :: d0 = 133681.0_wp  ! additive constant

    INTEGER   :: jb, jc, je, nblks_e, nblks_c, npromz_e, npromz_c, nlen
    REAL(wp)  :: z_lat, z_lon, z_hh, z_uu, z_vv, z_or
    REAL(wp)  :: z_phi_t_k
    REAL(wp)  :: z_summand
    REAL(wp)  :: z_fact
    REAL(wp)  :: z_angle1
    REAL(wp)  :: z_angle2
    REAL(wp)  :: z_angle

    u0 = (2.0_wp*pi*grid_sphere_radius)/(12.0_wp*rdaylen) ! [m/s]

  !-----------------------------------------------------------------------
    ! set reference height
    sw_ref_height =  d0 *rgrav

    nblks_c   = pt_patch%nblks_c
    npromz_c  = pt_patch%npromz_c
    nblks_e   = pt_patch%nblks_e
    npromz_e  = pt_patch%npromz_e

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1, nlen
        z_lat   = pt_patch%cells%center(jc,jb)%lat
        z_lon   = pt_patch%cells%center(jc,jb)%lon

        ! height of orography
        z_fact = grid_sphere_radius * grid_angular_velocity * SIN(z_lat)
        z_fact = z_fact * z_fact
        z_or = .5_wp * z_fact * rgrav
        pt_ext_data%atm%topography_c(jc,jb)= z_or

        ! relevant angles
        z_angle1 = .25_wp * pi
        z_angle2 = z_lon

        ! 1st summand: \phi_t(\vec c) \cdot \vec k
        z_phi_t_k = SIN(z_lat) * COS(z_angle1)
        z_phi_t_k = z_phi_t_k - COS(z_angle2) * COS(z_lat) * SIN(z_angle1)
        z_phi_t_k = u0 * z_phi_t_k

        ! 2nd summand: r_e \grid_angular_velocity \sin \varphi
        z_summand = grid_sphere_radius * grid_angular_velocity * SIN(z_lat)

        ! one factor
        z_fact    = .5_wp *  z_phi_t_k + z_summand

        ! height
        z_hh      = d0 - z_phi_t_k *  z_fact
        z_hh      = z_hh * rgrav

        ! thickness = height - oro
        pt_prog%pres_sfc(jc,jb) = z_hh -z_or
      ENDDO
    ENDDO

    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO je = 1, nlen
        z_lat   = pt_patch%edges%center(je,jb)%lat
        z_lon   = pt_patch%edges%center(je,jb)%lon

        z_angle1 = .25_wp * pi
        z_angle2 = z_lon
        z_uu = COS(z_lat) * COS(z_angle1)
        z_uu = z_uu + COS(z_angle2) * SIN(z_lat) * SIN(z_angle1)
        z_uu = u0 * z_uu

        z_angle = z_lon
        z_vv = SIN(z_angle)
        z_angle = .25_wp * pi
        z_vv = z_vv * SIN(z_angle)
        z_vv = -1._wp * u0 * z_vv

        pt_prog%vn(je,1,jb) = &
            z_uu * pt_patch%edges%primal_normal(je,jb)%v1 &
           +z_vv * pt_patch%edges%primal_normal(je,jb)%v2

      ENDDO
    ENDDO

  END SUBROUTINE init_usbr_test

  !-------------------------------------------------------------------------
  !>
  !! Performs  numerical integration between -@f$\frac{\pi}{2}@f$ and @f$\frac{\pi}{2}@f$.
  !!
  !! Performs  numerical integration between -@f$\frac{\pi}{2}@f$ and @f$\frac{\pi}{2}@f$
  !! to compute geostrophically balanced initial state used
  !! in test 3.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2002-5).
  !! Modified by Th.Heinze, DWD, (2006-11-22):
  !! - introduced INTERFACE uu (got an error message with g95 compiler,
  !!   scanned the code, this seems to be the correct way, but might be wrong)
  !! Modified by Th.Heinze, DWD, (2006-12-12):
  !! - renamed it to geostr_balance
  !!
  !! @par Remarks
  !! was htmp2 in previous code
  !!
  FUNCTION geostr_balance( p_lat, func, inverse_sphere_radius, angular_velocity, u0)  RESULT(p_hh)

    INTERFACE                        ! selected function

      FUNCTION func(p_t, u0) RESULT(p_vv)

        USE mo_kind, ONLY: wp
        REAL(wp), INTENT(in) :: p_t, u0
        REAL(wp)             :: p_vv

      END FUNCTION func

    END INTERFACE

    REAL(wp), INTENT(in) :: p_lat           ! rotated latitude
    REAL(wp), INTENT(in) :: inverse_sphere_radius, angular_velocity, u0

    REAL(wp)             :: p_hh            ! balanced height

    INTEGER              :: j               ! loop index

    REAL(wp)             :: z_a             ! left bound
    REAL(wp)             :: z_b             ! right bound
    REAL(wp)             :: z_lat           ! latitude in loop
    REAL(wp)             :: z_step          ! step
    REAL(wp)             :: z_val, z_val2   ! intermediate values


  !-----------------------------------------------------------------------

    z_a = -1._wp * pi_2
    z_b = p_lat

    z_step = 0.02_wp * ( z_b - z_a)

    p_hh = 0._wp

    z_lat = z_a - 0.5_wp * z_step

    DO j = 1, 50
       z_lat = z_lat + z_step

       z_val = func(z_lat, u0)

       z_val2 = 2._wp * angular_velocity * SIN(z_lat)
       z_val2 = z_val2 + z_val * TAN(z_lat)* inverse_sphere_radius
       z_val2 = z_val * z_val2

       p_hh = p_hh + z_val2 * z_step

    ENDDO

  END FUNCTION geostr_balance

  !-------------------------------------------------------------------------
  !>
  !! Specifies zomnally symmetric u field as a function of latitude p_latd.
  !!
  !! Specifies zomnally symmetric u field as a function of latitude p_latd
  !! (bump function). Requires real function bifurc.<br>
  !! Used for test case 3 (see paper by Browning et al.)
  !!
  !! @par Revision History
  !! Developed originally by R.Jakob for NCAR shallow water model.
  !! Adapted to ICON code by L.Bonaventura  (2002-5).
  !! Modified by Th.Heinze, DWD, (2006-12-14):
  !! - renamed it to symmetric_u_velo
  !!
  !! @par Remarks
  !! was us in previous code
  !!
  FUNCTION symmetric_u_velo(p_rlatd, u0) RESULT (p_usres)

    REAL(wp),INTENT(in) :: p_rlatd, u0  ! LATITUDE

    REAL(wp)            :: p_usres  ! MAX AMPLITUDE OF FLOW
                                    ! (12 DAY ROTATION SPEED)
    REAL(wp)            :: z_rlate, z_rlatb  ! NORTH-SOUTH EXTENT OF FLOW FIELD
    REAL(wp)            :: z_xe, z_x         ! FLOW PROFILE TEMPORARIES

  !-----------------------------------------------------------------------

    z_rlatb = pi_2/(-3.0_wp)          ! -pi/6
    z_rlate = pi_2
    z_xe    = 0.3_wp
    z_x     = z_xe * (p_rlatd - z_rlatb) / (z_rlate-z_rlatb)
    p_usres = bifurc(z_x) * bifurc(z_xe - z_x)
    p_usres = u0 * EXP(4.0_wp/z_xe) * p_usres

  END FUNCTION symmetric_u_velo

  !-------------------------------------------------------------------------
  !>
  !! This subroutine computes the rotated coordinates p_rotlon, p_rotlat.
  !!
  !! This subroutine computes the rotated coordinates p_rotlon, p_rotlat
  !! for a roatation by angle p_alpha, given the coordinates p_lon and p_lat.
  !!
  !! @par Revision History
  !! Developed originally by R.Jakob for NCAR shallow water model.
  !! Adapted to ICON code by L.Bonaventura (2002-5).
  !! Adapted to ICON programming guide by Th.Heinze, DWD, (2006-12-12)
  !!
  SUBROUTINE rotate(p_lon, p_lat, p_alpha, p_rotlon, p_rotlat)

    REAL(wp), INTENT(in)  :: p_lon     ! ORIGINAL LONGITUDE
    REAL(wp), INTENT(in)  :: p_lat     ! ORIGINAL LATITUDE
    REAL(wp), INTENT(in)  :: p_alpha   ! ROTATION ANGLE

    REAL(wp), INTENT(out) :: p_rotlon  ! ROTATED LONGITUDE
    REAL(wp), INTENT(out) :: p_rotlat  ! ROTATED LATITUDE

    REAL(wp)              :: z_test    ! checking value

  !-----------------------------------------------------------------------

    IF (p_alpha == 0.0_wp) THEN       !        NO ROTATION

      p_rotlon = p_lon
      p_rotlat = p_lat

    ELSE                              !        ROTATION BY ANGLE p_alpha

!     ROTATED LATITUDE

      z_test = SIN(p_lat)*COS(p_alpha)- COS(p_lat)*COS(p_lon)*SIN(p_alpha)

      IF (z_test > 1.0_wp) THEN
        p_rotlat = pi_2
      ELSEIF (z_test < -1.0_wp) THEN
        p_rotlat = -1.0_wp * pi_2
      ELSE
        p_rotlat = ASIN(z_test)
      ENDIF

!     ROTATED LONGITUDE

      z_test = COS(p_rotlat)

      IF (z_test == 0.0_wp) THEN
        p_rotlon = 0.0_wp
      ELSE
        z_test = SIN(p_lon)*COS(p_lat)/z_test
        IF (z_test > 1.0_wp) THEN
          p_rotlon = pi_2
        ELSEIF (z_test < -1.0_wp) THEN
          p_rotlon = -1.0_wp * pi_2
        ELSE
          p_rotlon = ASIN(z_test)
        ENDIF
      ENDIF

!        ADJUST FOR CORRECT BRANCH OF INVERSE SINE

      z_test = COS(p_alpha)*COS(p_lon)*COS(p_lat) + SIN(p_alpha)*SIN(p_lat)

      IF (z_test < 0.0_wp) THEN
        p_rotlon = pi - p_rotlon
      ENDIF

    ENDIF

  END SUBROUTINE rotate

  !-------------------------------------------------------------------------
  !>
  !! Function .
  !!
  !! Function
  !! @f{align*}{
  !! bifurc (p\_x) &= \begin{cases}
  !! e^{-\frac{1}{p\_x}} &\text{if } p\_x > 0 \\ 0 &\text{if } p\_x \leq 0
  !! \end{cases}
  !! @f}
  !! Used in conjunction with function symmetric_u_velo
  !!
  !! @par Revision History
  !! Developed originally by R.Jakob for NCAR shallow water model.
  !! Adapted to ICON code by L.Bonaventura  (2002-5).
  !! Modified by Th.Heinze, DWD, (2006-12-14):
  !! - renamed it to bifurc
  !!
  !! @par Remarks
  !! was bf2 in previous code
  !!
  FUNCTION bifurc(p_x) RESULT (p_bfres)


    REAL(wp), INTENT(in) :: p_x       ! x value

    REAL(WP)             :: p_bfres   ! result
  !-----------------------------------------------------------------------

    IF (p_x <= 0.0_wp) THEN
      p_bfres = 0.0_wp
    ELSE
      p_bfres = EXP(-1.0_wp/p_x)
    ENDIF

  END FUNCTION bifurc

  !-------------------------------------------------------------------------
END MODULE mo_sw_test


