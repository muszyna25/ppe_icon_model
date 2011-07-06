!>
!! Defines the artificial testcases for the nonhydrostatic atmospheric model.
!! 
!! 
!! @par Revision History
!! Initial release by Almut Gassmann (2008-03-18)
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
!! 
MODULE mo_nh_testcases  
!-------------------------------------------------------------------------  
!  
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006                       
!  
!-------------------------------------------------------------------------  
!  
!  
!  
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_namelist,           ONLY: position_nml, POSITIONED
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_model_domain_import,ONLY: lplane, n_dom
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data,           ONLY: ext_data
  USE mo_math_constants,     ONLY: pi, pi_2
  USE mo_math_utilities,     ONLY: gc2cc, t_cartesian_coordinates, &
                                   t_geographical_coordinates, &
                                   arc_length
  USE mo_parallel_configuration,  ONLY: nproma
  USE mo_run_nml,            ONLY: ltransport, ntracer, iforcing, inwp, &
    &                              iqv, i_cell_type,itopo
  USE mo_grid_configuration, ONLY :  global_cell_type
    
  USE mo_dynamics_nml,       ONLY: nnow,nnew,nnow_rcf, nnew_rcf
  USE mo_atm_phy_nwp_nml,    ONLY: inwp_gscp, inwp_convection
  USE mo_physical_constants, ONLY: grav, cpd, rd, cvd_o_rd, &
   &                               p0ref, re, omega, tmelt, vtmpc1, rv
! USE mo_convect_tables,     ONLY: B1   => c1es,  &
!  &                               B2_w => c3les ,&
!  &                               B2_i => c3ies ,&
!  &                               B4_w => c4les ,&
!  &                               B4_i => c4ies
  USE mo_nonhydro_state,       ONLY: t_nh_state
  
  
  USE mo_interpolation,        ONLY: t_int_state, cells2edges_scalar, edges2cells_scalar
  USE mo_mpi,                  ONLY: p_pe, p_io
  USE mo_vertical_coord_table, ONLY: vct_b
  USE mo_loopindices,          ONLY: get_indices_e, get_indices_c
  USE mo_advection_nml,        ONLY: ctracer_list
  USE mo_ncar_testcases,       ONLY: tracer_q1_q2, tracer_q3
  USE mo_nh_pa_test,           ONLY: init_nh_state_prog_patest
  USE mo_nh_df_test,           ONLY: init_nh_state_prog_dftest
  USE mo_nh_hs_test,           ONLY: init_nh_state_prog_held_suarez
  USE mo_nh_ape_exp,           ONLY: init_nh_state_prog_APE
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_parallel_configuration, ONLY: p_test_run
  USE mo_sync,                 ONLY: SYNC_E, sync_patch_array
  USE mo_satad,                ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                                sat_pres_ice!,  &  !! saturation vapor pressure w.r.t. ice
!   &                                spec_humi!,     &  !! Specific humidity
  USE mo_nh_prog_util,         ONLY: nh_prog_add_random

  
  IMPLICIT NONE  
  
  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  
  CHARACTER(len=MAX_CHAR_LENGTH) :: nh_test_name
  CHARACTER(len=MAX_CHAR_LENGTH) :: ape_sst_case      !SST for APE experiments

  REAL(wp) :: mount_height           ! (m)
  REAL(wp) :: layer_thickness        ! (m)
  REAL(wp) :: nh_brunt_vais          ! (1/s)
  REAL(wp) :: nh_u0                  ! (m/s)
  REAL(wp) :: nh_t0                  ! (K)
  REAL(wp) :: jw_up                  ! amplitude of the u-perturbation (m/s), jabw  
  REAL(wp) :: u0_mrw                 ! (m/s) wind speed for mrw case 
  REAL(wp) :: mount_height_mrw       ! (m) maximum mount height in mrw and mwbr
  REAL(wp) :: mount_half_width       ! (m) half width of mountain in mrw, mwbr and bell
  REAL(wp) :: mount_lonctr_mrw_deg   ! (deg) lon of mountain center in mrw and mwbr
  REAL(wp) :: mount_latctr_mrw_deg   ! (deg) lat of mountain center in mrw and mwbr
  REAL(wp) :: p_int_mwbr_const       ! pressure at interface in mwbr_const test case
  REAL(wp) :: temp_i_mwbr_const      ! temp in isothermal lower layer in mwbr_const
  REAL(wp) :: bruntvais_u_mwbr_const ! brunt vaisala freq in upper layer in mwbr_const
  REAL(wp) :: u0_mwbr_const          ! (m/s) wind speed for mwbr_const case
  REAL(wp) :: rotate_axis_deg        ! (deg) rotation angle
  REAL(wp) :: torus_domain_length    ! (m) length of domain the slice (torus) grid
  LOGICAL  :: lhs_nh_vn_ptb          ! if true, random noise is added to vn in HS_nh test case
  LOGICAL  :: lhs_fric_heat          ! if true, frictional heating is switched on in HS_nh
  REAL(wp) :: hs_nh_vn_ptb_scale     ! amplitude of the random noise
  REAL(wp) :: rh_at_1000hpa          ! relative humidity at 1000 hPa [%]
  REAL(wp) :: qv_max                 ! limit of maximum specific humidity in the tropics [kg/kg]

  LOGICAL  :: linit_tracer_fv  !< finite volume initialization for tracer fields
                               !< if .TRUE.

  INTEGER  :: n_flat_level

  NAMELIST/nh_testcase_ctl/ nh_test_name, mount_height, torus_domain_length, &
                            nh_brunt_vais, nh_u0, nh_t0, layer_thickness,    &
                            n_flat_level, jw_up, u0_mrw, mount_height_mrw,   &
                            mount_half_width, mount_lonctr_mrw_deg,          &
                            mount_latctr_mrw_deg, p_int_mwbr_const,          &
                            temp_i_mwbr_const,  bruntvais_u_mwbr_const,      &
                            u0_mwbr_const, rotate_axis_deg,                  &
                            lhs_nh_vn_ptb, hs_nh_vn_ptb_scale,               &
                            rh_at_1000hpa, qv_max, ape_sst_case,             &
                            linit_tracer_fv, lhs_fric_heat

  PUBLIC :: setup_nh_testcase, layer_thickness, init_nh_testtopo,            &
    &       init_nh_testcase, n_flat_level, nh_test_name, ape_sst_case,      &
    &       mount_height, torus_domain_length, nh_brunt_vais, nh_u0, nh_t0,  &
    &       jw_up, u0_mrw, mount_height_mrw, mount_half_width,               &
    &       mount_lonctr_mrw_deg, mount_latctr_mrw_deg, p_int_mwbr_const,    &
    &       temp_i_mwbr_const,  bruntvais_u_mwbr_const, u0_mwbr_const,       &
    &       rotate_axis_deg, lhs_nh_vn_ptb, hs_nh_vn_ptb_scale,              &
    &       rh_at_1000hpa, qv_max, linit_tracer_fv, lhs_fric_heat

  PRIVATE

! !DEFINED PARAMETERS for jablonowski williamson:  
  REAL(wp), PARAMETER :: eta0  = 0.252_wp ! 
  REAL(wp), PARAMETER :: etat  = 0.2_wp   ! tropopause
  REAL(wp), PARAMETER :: ps0   = 1.e5_wp  ! surface pressure (Pa)
  REAL(wp), PARAMETER :: u0    = 35._wp   ! maximum zonal wind (m/s)
  REAL(wp), PARAMETER :: temp0 = 288._wp  ! horizontal-mean temperature 
                                          ! at surface (K)
  REAL(wp), PARAMETER :: gamma = 0.005_wp ! temperature elapse rate (K/m)
  REAL(wp), PARAMETER :: dtemp = 4.8e5_wp ! empirical temperature difference (K)

  REAL(wp), PARAMETER :: lonC  = pi/9._wp ! longitude of the perturb. centre 
  REAL(wp), PARAMETER :: latC  = 2._wp*lonC !latitude of the perturb. centre
  
! !DEFINED PARAMETERS for mountain induced Rossby wave train:
  REAL(wp), PARAMETER :: pres_sp  = 93000.0_wp  !pressure surface at the south pole
  REAL(wp), PARAMETER :: temp_mrw = 288._wp     !temperature of isothermal atmosphere

  CONTAINS  
  
!-------------------------------------------------------------------------  
!-------------------------------------------------------------------------
!
!
  !>
  !! Defines nonhydrostatic artificial initial conditions.
  !! 
  !! Reads namelist
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-19)
  !! Modification by Daniel Reinert, DWD (2010-07-15)
  !! - moved initialization of topography into new subroutine 
  !!   init_nh_testtopo, which is called after the domain-decomposition.
  !!   (because of possible conflicts with the external-data type)
  !! 
  SUBROUTINE setup_nh_testcase
!
    INTEGER        :: i_status

!-----------------------------------------------------------------------

    ! default values
    nh_test_name           = 'jabw'
    mount_height           = 100.0_wp
    layer_thickness        = -999.0_wp
    n_flat_level           = 2
    nh_u0                  = 0.0_wp
    nh_brunt_vais          = 0.01_wp
    nh_t0                  = 300.0_wp
    jw_up                  = 1.0_wp
    u0_mrw                 = 20.0_wp
    mount_height_mrw       = 2000.0_wp
    mount_half_width       = 1500000._wp
    mount_lonctr_mrw_deg   = 90.0_wp
    mount_latctr_mrw_deg   = 30.0_wp
    u0_mwbr_const          = 20.0_wp
    p_int_mwbr_const       = 70000._wp
    temp_i_mwbr_const      = 288._wp
    bruntvais_u_mwbr_const = 0.025_wp
    rotate_axis_deg        = 0.0_wp
    torus_domain_length    = 100000.0_wp
    lhs_nh_vn_ptb          = .TRUE.
    lhs_fric_heat          = .FALSE.
    hs_nh_vn_ptb_scale     = 1._wp  ! magnitude of the random noise
    rh_at_1000hpa          = 0.7_wp
    qv_max                 = 20.e-3_wp ! 20 g/kg
    ape_sst_case           = 'sst1'
    IF(global_cell_type==3) THEN
      linit_tracer_fv        = .TRUE. ! finite volume initialization for tracer
    ELSE
      linit_tracer_fv        = .FALSE.
    ENDIF


    CALL position_nml ('nh_testcase_ctl', status=i_status)
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nh_testcase_ctl)
    END SELECT

    n_flat_level=MAX(2,n_flat_level)

    IF(p_pe == p_io) WRITE(nnml_output,nml=nh_testcase_ctl)

  END SUBROUTINE setup_nh_testcase


!-------------------------------------------------------------------------
!
!
  !>
  !! Initialize topography for nonhydrostatic artificial testcases.
  !! 
  !! Initialize topography for nonhydrostatic artificial testcases
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-19)
  !! Modification by Daniel Reinert, DWD (2010-07-15)
  !! - moved initialization of topography into new subroutine 
  !!   init_nh_testtopo, which is called after the domain-decomposition.
  !!   (because of possible conflicts with the external-data type)
  !! 
  SUBROUTINE init_nh_testtopo (p_patch)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)

  INTEGER        :: jg, jc, jv, jb, nlen
  REAL(wp)       :: z_lon, z_lat, z_dist
  REAL(wp)       :: zr, zexp
  TYPE(t_geographical_coordinates) :: z_x2_geo
  TYPE(t_cartesian_coordinates)    :: z_x1_cart, z_x2_cart
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = '(mo_nh_testcases) init_nh_testtopo:'
 
  REAL(wp)      :: zsiny, zcosy,tmp1,tmp2,tmp3
  REAL(wp)      :: z_lon_ctr, z_lat_ctr, z_fac1, z_fac2

!-----------------------------------------------------------------------

  ! Initialize topography to zero if idealized topo is used
  IF ( itopo == 0 ) THEN
    DO jg = 1, n_dom
      ext_data(jg)%atm%topography_c(1:nproma,1:p_patch(jg)%nblks_int_c) = 0.0_wp
      ext_data(jg)%atm%topography_v(1:nproma,1:p_patch(jg)%nblks_int_v) = 0.0_wp
    ENDDO
  ENDIF

  SELECT CASE (nh_test_name)

  CASE ('zero', 'HS_jw')
    
    IF(nh_test_name=='HS_jw') CALL message(TRIM(routine),'running the Held-Suarez test')

  CASE ('schaer')
 
    IF(.NOT.lplane) CALL finish(TRIM(routine),'Schaer test case only for lplane=True')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_int_c
        IF (jb /=  p_patch(jg)%nblks_int_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          z_lon  = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
          z_dist = z_lon-z_x2_geo%lon
          ext_data(jg)%atm%topography_c(jc,jb) = 250.0_wp &
          & * EXP(-(z_dist/5000.0_wp)**2)*((COS(pi*z_dist/4000.0_wp))**2)
        ENDDO
      ENDDO 
    ENDDO 
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_int_v
        IF (jb /=  p_patch(jg)%nblks_int_v) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_v
        ENDIF
        DO jv = 1, nlen
          z_lon  = p_patch(jg)%verts%vertex(jv,jb)%lon*torus_domain_length/pi*0.5_wp
          z_dist = z_lon-z_x2_geo%lon
          ext_data(jg)%atm%topography_v(jv,jb) = 250.0_wp &
          & * EXP(-(z_dist/5000.0_wp)**2)*((COS(pi*z_dist/4000.0_wp))**2)
        ENDDO
      ENDDO 
    ENDDO

  CASE ('bell')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_int_c
        IF (jb /=  p_patch(jg)%nblks_int_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (lplane) THEN
            z_lat  = 0.0_wp!p_patch(jg)%cells%center(jc,jb)%lat*torus_domain_length/pi*0.5_wp
            z_lon  = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
            z_dist = SQRT((z_lat-z_x2_geo%lat)**2+(z_lon-z_x2_geo%lon)**2)
          ELSE
            z_x1_cart    = gc2cc(p_patch(jg)%cells%center(jc,jb))
            z_x2_cart    = gc2cc(z_x2_geo)
            z_dist       = arc_length(z_x1_cart,z_x2_cart)
          ENDIF
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height/ &
                 (1.0_wp+ (z_dist/mount_half_width)**2)**1.5_wp
        ENDDO
      ENDDO 
    ENDDO 
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_int_v
        IF (jb /=  p_patch(jg)%nblks_int_v) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_v
        ENDIF
        DO jv = 1, nlen
          IF (lplane) THEN
            z_lat  = 0.0_wp!p_patch(jg)%verts%vertex(jv,jb)%lat*torus_domain_length/pi*0.5_wp
            z_lon  = p_patch(jg)%verts%vertex(jv,jb)%lon*torus_domain_length/pi*0.5_wp
            z_dist = SQRT((z_lat-z_x2_geo%lat)**2+(z_lon-z_x2_geo%lon)**2)
          ELSE
            z_x1_cart    = gc2cc(p_patch(jg)%verts%vertex(jv,jb))
            z_x2_cart    = gc2cc(z_x2_geo)
            z_dist       = arc_length(z_x1_cart,z_x2_cart)
          ENDIF
          ext_data(jg)%atm%topography_v(jv,jb) = mount_height/ &
                 (1.0_wp+ (z_dist/mount_half_width)**2)**1.5_wp
        ENDDO
      ENDDO 
    ENDDO

  CASE ('jabw', 'jabw_s', 'jabw_m')

    DO jg = 1, n_dom 
      DO jb = 1, p_patch(jg)%nblks_int_c
        IF (jb /=  p_patch(jg)%nblks_int_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          z_lat   = p_patch(jg)%cells%center(jc,jb)%lat
          zsiny = SIN(z_lat)
          zcosy = COS(z_lat)
          tmp1  = u0*COS((1._wp-eta0)*pi_2)**1.5_wp
          tmp2  = (-2.0_wp*zsiny**6 * (zcosy*zcosy+1.0_wp/3.0_wp) + &
                  1.0_wp/6.3_wp ) *tmp1
          tmp3  = ( 1.6_wp*zcosy*zcosy*zcosy * (zsiny*zsiny+2.0_wp/3.0_wp)  &
                   - 0.5_wp*pi_2 )*re*omega
          IF ( itopo==0 ) ext_data(jg)%atm%topography_c(jc,jb) = tmp1*(tmp2+tmp3)/grav
          IF (itopo==0 .AND. nh_test_name =='jabw_m' ) THEN
            z_lon = p_patch(jg)%cells%center(jc,jb)%lon
            z_fac1= SIN(latC)*SIN(z_lat)+COS(latC)*COS(z_lat)*COS(z_lon-lonC) 
            z_fac2= re*ACOS(z_fac1)/mount_half_width
            ext_data(jg)%atm%topography_c(jc,jb) = ext_data(jg)%atm%topography_c(jc,jb) &
            & + mount_height*EXP(-z_fac2**2)
          ENDIF 
        ENDDO
      ENDDO
    ENDDO
    DO jg = 1, n_dom 
      DO jb = 1, p_patch(jg)%nblks_int_v
        IF (jb /=  p_patch(jg)%nblks_int_v) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_v
        ENDIF
        DO jv = 1, nlen
          z_lat   = p_patch(jg)%verts%vertex(jv,jb)%lat
          zsiny = SIN(z_lat)
          zcosy = COS(z_lat)
          tmp1  = u0*COS((1._wp-eta0)*pi_2)**1.5_wp
          tmp2  = (-2.0_wp*zsiny**6 * (zcosy*zcosy+1.0_wp/3.0_wp) + &
                  1.0_wp/6.3_wp ) *tmp1
          tmp3  = ( 1.6_wp*zcosy*zcosy*zcosy * (zsiny*zsiny+2.0_wp/3.0_wp)  &
                   - 0.5_wp*pi_2 )*re*omega
          IF ( itopo==0 ) ext_data(jg)%atm%topography_v(jv,jb) = tmp1*(tmp2+tmp3)/grav
          IF (itopo==0 .AND. nh_test_name =='jabw_m' ) THEN
            z_lon = p_patch(jg)%verts%vertex(jv,jb)%lon
            z_fac1= SIN(latC)*SIN(z_lat)+COS(latC)*COS(z_lat)*COS(z_lon-lonC) 
            z_fac2= re*ACOS(z_fac1)/mount_half_width
            ext_data(jg)%atm%topography_v(jv,jb) = ext_data(jg)%atm%topography_v(jv,jb) &
            & + mount_height*EXP(-z_fac2**2)
          ENDIF 
        ENDDO
      ENDDO
    ENDDO

  CASE ('mrw_nh', 'mrw2_nh' , 'mwbr_const')

   CALL message(TRIM(routine),'running mrw, setting topography')
    z_lon_ctr = mount_lonctr_mrw_deg*pi/180.0_wp
    z_lat_ctr = mount_latctr_mrw_deg*pi/180.0_wp

    DO jg = 1, n_dom 
      DO jb = 1, p_patch(jg)%nblks_int_c
        IF (jb /=  p_patch(jg)%nblks_int_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          z_lat   = p_patch(jg)%cells%center(jc,jb)%lat
          z_lon   = p_patch(jg)%cells%center(jc,jb)%lon

          zr = SIN(z_lat_ctr)*SIN(z_lat)+COS(z_lat_ctr)*COS(z_lat)*COS(z_lon-z_lon_ctr) 
          zexp = re*ACOS(zr)/mount_half_width

          IF ( itopo==0 ) THEN
            IF (nh_test_name=='mrw_nh' .OR. nh_test_name=='mwbr_const'  ) THEN
              ext_data(jg)%atm%topography_c(jc,jb) = &
                        mount_height_mrw*EXP( - zexp*zexp )
            ELSEIF (nh_test_name=='mrw2_nh') THEN
              ext_data(jg)%atm%topography_c(jc,jb) = &
                         mount_height_mrw*EXP( - zexp*zexp )*&
                         0.5_wp*(1._wp+COS(pi*zexp*2._wp))
            ENDIF
          ENDIF
          
        ENDDO
      ENDDO
    ENDDO
    DO jg = 1, n_dom 
      DO jb = 1, p_patch(jg)%nblks_int_v
        IF (jb /=  p_patch(jg)%nblks_int_v) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_v
        ENDIF
        DO jv = 1, nlen
          z_lat   = p_patch(jg)%verts%vertex(jv,jb)%lat
          z_lon   = p_patch(jg)%verts%vertex(jv,jb)%lon

          zr = SIN(z_lat_ctr)*SIN(z_lat)+COS(z_lat_ctr)*COS(z_lat)*COS(z_lon-z_lon_ctr) 
          zexp = re*ACOS(zr)/mount_half_width

          IF ( itopo==0 ) THEN
            IF(nh_test_name=='mrw_nh' .OR. nh_test_name=='mwbr_const') THEN
              ext_data(jg)%atm%topography_v(jv,jb) =  &
                        mount_height_mrw*EXP( - zexp*zexp )
            ELSEIF (nh_test_name=='mrw2_nh') THEN
              ext_data(jg)%atm%topography_v(jv,jb) =  &
                        mount_height_mrw*EXP( - zexp*zexp )*&
                        0.5_wp*(1._wp+COS(pi*zexp*2._wp))
            ENDIF
          ENDIF  

        ENDDO
      ENDDO
    ENDDO

   CALL message(TRIM(routine),'topography is initialised ')

  CASE ('PA')
   ! The topography ist initialized in "init_nh_state_prog_patest"
    CALL message(TRIM(routine),'running the Pure 3D-Advection test.')

  CASE ('DF1')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 1.')

  CASE ('DF2')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 2.')

  CASE ('DF3')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 3.')

  CASE ('DF4')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 4.')

  CASE ('HS_nh')
    ! The topography ist initialized in "init_nh_state_prog_held_suarez"
    CALL message(TRIM(routine),'running the Held-Suarez test')
  CASE ('APE_nh')

   ! The topography ist initialized in "init_nh_state_prog_APE"
    CALL message(TRIM(routine),'running Aqua-Planet Experiment')

  CASE DEFAULT

    CALL finish(routine,'wrong input for nh_test_name')

  END SELECT


  END SUBROUTINE init_nh_testtopo



!-------------------------------------------------------------------------
!
!
  !>
  !! Defines nonhydrostatic artificial initial conditions.
  !! 
  !! Initializes meteorological fields
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-04-14)
  !! 
  SUBROUTINE init_nh_testcase (p_patch, p_nh_state, p_int, ntl)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_int_state), INTENT(IN) :: p_int(n_dom)
  INTEGER :: ntl
  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)

  INTEGER        :: jg, je, jc, jb, jk, jt, jjt, jn, ji, niter, &
                    nlen, nblks_e, npromz_e,  nblks_c, npromz_c
  INTEGER        :: i_startidx, i_endidx, i_startblk, icount
  INTEGER        :: nlev, nlevp1        !< number of full and half levels

  REAL(wp), DIMENSION(nproma) ::  &
             z_lat,z_siny,z_cosy, z_fac1, z_fac2, zeta_old, zcoszetav, &
             zsinzetav, z_tavg, z_favg, z_geopot, z_temp, z_fun, z_fund, &
             zeta, zu, zv, z_lon, z_exp, za, zb, zc, z_temp_kp1, z_fac3
  REAL(wp), ALLOCATABLE :: zeta_v(:,:,:)
  REAL(wp), ALLOCATABLE :: zeta_v_e(:,:,:)
  REAL(wp), ALLOCATABLE :: z_wsfc_e(:,:,:), z_wsfc_c(:,:,:)
  REAL(wp), ALLOCATABLE :: z_int_c(:,:)    ! z at interface in mwbr_const
  TYPE(t_nh_state), POINTER       :: p_nhdom
  REAL(wp)              :: zlat, zlon, z_u, bruntvaissq, kappa, zhelp1, &
                           zcoslat, zhelp2, z_pres, z_sfc, z_nlev
  REAL(wp)              :: zhelp1_i, zhelp1_u, zhelp2_i, bruntvaissq_i, &
                           bruntvaissq_u, zhelp3, zhelp4, rkappa
  REAL(wp)              :: theta_v_int !potential temp at the interface in mwbr_const

  REAL(wp) :: z_help

  REAL(wp) :: zsqv
  
  REAL(wp) :: zrhf,z_1_o_rh
  
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                                   '(mo_nh_testcases) init_nh_testcase:' 
! Tracer related variables
  CHARACTER(LEN=1) :: ctracer

!-----------------------------------------------------------------------

  SELECT CASE (nh_test_name)

  CASE ('jabw', 'jabw_s', 'jabw_m', 'APE_nh', 'HS_jw') ! Jablonowski test

  CALL message(TRIM(routine),'Jablonowski test')
  IF ( iforcing == inwp ) THEN
    CALL message(TRIM(routine),' iforcing == inwp')
  ELSE
    CALL message(TRIM(routine),'Attention: iforcing /= inwp')
  ENDIF
  
  DO jg = 1, n_dom

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev

    ALLOCATE (zeta_v(nproma,nlev,p_patch(jg)%nblks_c), &
              zeta_v_e(nproma,nlev,p_patch(jg)%nblks_e) )
    
    zeta_v    = 0._wp
    zeta_v_e  = 0._wp
    nblks_c   = p_patch(jg)%nblks_int_c
    npromz_c  = p_patch(jg)%npromz_int_c
    nblks_e   = p_patch(jg)%nblks_int_e
    npromz_e  = p_patch(jg)%npromz_int_e
    p_nhdom  => p_nh_state(jg)
    p_nhdom%diag%pres_sfc(:,:) = 100000._wp     !set surface pressure to 1000. hPa

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jn,jt,z_lat,z_siny,z_cosy,z_fac1,z_fac2,z_exp,zeta_old,&
!$OMP            zcoszetav,zsinzetav,z_tavg,z_favg,z_geopot,z_temp,z_fun,z_fund,zeta )
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = nlev, 1, -1
        DO jc = 1, nlen
          z_lat(jc) = p_patch(jg)%cells%center(jc,jb)%lat
          z_siny(jc) = SIN(z_lat(jc))
          z_cosy(jc) = COS(z_lat(jc))
          z_fac1(jc) = 1.0_wp/6.3_wp-2.0_wp*(z_siny(jc)**6)*(z_cosy(jc)**2+1.0_wp/3.0_wp)
          z_fac2(jc) = (8.0_wp/5.0_wp*(z_cosy(jc)**3)*(z_siny(jc)**2+2.0_wp/3.0_wp)&
                       -0.25_wp*pi)*re*omega
          z_exp(jc)  = rd*gamma/grav
          zeta_old(jc) = 1.0e-7_wp
        ENDDO
        ! Newton iteration to determine zeta
        DO jn = 1, 100
          DO jc = 1, nlen
            zeta_v(jc,jk,jb) = (zeta_old(jc) - eta0)*pi_2
            zcoszetav(jc)= COS(zeta_v(jc,jk,jb))
            zsinzetav(jc)= SIN(zeta_v(jc,jk,jb))
            z_tavg(jc)   = temp0*(zeta_old(jc)**z_exp(jc))
            z_favg(jc)   = temp0*grav/gamma*(1.0_wp-zeta_old(jc)**z_exp(jc))
            IF (zeta_old(jc) < etat ) THEN
               z_tavg(jc) = z_tavg(jc)+dtemp*((etat-zeta_old(jc))**5)
               z_favg(jc) = z_favg(jc)-rd*dtemp*(                           &
                        (log(zeta_old(jc)/etat)+137.0_wp/60.0_wp)*(etat**5) &
                        -5.0_wp*(etat**4)*zeta_old(jc)                      &
                        +5.0_wp*(etat**3)*(zeta_old(jc)**2)                 &
                        -10.0_wp/3.0_wp*(etat**2)*(zeta_old(jc)**3)         &
                        +1.25_wp*etat*(zeta_old(jc)**4)-0.2_wp*(zeta_old(jc)**5))
            ENDIF
            z_geopot(jc) = z_favg(jc)+u0*(zcoszetav(jc)**1.5_wp)*&
                          (z_fac1(jc)*u0*(zcoszetav(jc)**1.5_wp)+z_fac2(jc)) 
            z_temp(jc)   = z_tavg(jc)+0.75_wp*zeta_old(jc)*pi*u0/rd*zsinzetav(jc)*&
                       SQRT(zcoszetav(jc))*(2.0_wp*u0*z_fac1(jc)*(zcoszetav(jc)**1.5_wp) &
                       + z_fac2(jc))
            z_fun(jc)    = z_geopot (jc)- p_nhdom%metrics%geopot(jc,jk,jb)
            z_fund(jc)   = -rd/zeta_old(jc)*z_temp(jc)
            zeta(jc) = zeta_old(jc) - z_fun(jc)/z_fund(jc)
            zeta_old(jc) = zeta(jc)
          ENDDO ! jc
        ENDDO !jn
        ! Final update for zeta_v
        DO jc = 1, nlen
          zeta_v(jc,jk,jb) = (zeta_old(jc) - eta0)*pi_2
        ENDDO
        ! Save results for all time levels
        ! Use analytic expressions at all model level
        DO jt = 1, ntl
          DO jc = 1, nlen
            p_nhdom%prog(jt)%exner(jc,jk,jb) = zeta_old(jc)**(rd/cpd)
            p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
            &        p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
            p_nhdom%prog(jt)%theta_v(jc,jk,jb) = z_temp(jc) & 
            &        /p_nhdom%prog(jt)%exner(jc,jk,jb)
            p_nhdom%prog(jt)%rho(jc,jk,jb) = &
            &        p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
            &        /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
            IF (jt == 1 ) THEN
              p_nhdom%diag%pres(jc,jk,jb) = p0ref*p_nhdom%prog(jt)%exner(jc,jk,jb)**(cpd/rd) 
              p_nhdom%diag%temp(jc,jk,jb) = z_temp(jc)  
            ENDIF
          ENDDO !jc
        ENDDO !jt
      ENDDO !jk
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    IF (p_test_run) zeta_v_e = 0._wp
    CALL cells2edges_scalar(zeta_v,p_patch(jg),p_int(jg)%c_lin_e,zeta_v_e)
    CALL sync_patch_array(SYNC_E,p_patch(jg),zeta_v_e)

    i_startblk = p_patch(jg)%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,jt,z_lat,z_lon,zu,z_fac1,z_fac2,zv)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

      DO jt = 1, ntl
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_lat(je) = p_patch(jg)%edges%center(je,jb)%lat
            z_lon(je) = p_patch(jg)%edges%center(je,jb)%lon
            zu(je)    = u0*(COS(zeta_v_e(je,jk,jb))**1.5_wp)*(SIN(2.0_wp*z_lat(je))**2)
            IF ((nh_test_name =='jabw').OR.(nh_test_name=='HS_jw')) THEN
             z_fac1(je)= SIN(latC)*SIN(z_lat(je))+COS(latC)*COS(z_lat(je))*COS(z_lon(je)-lonC) 
             z_fac2(je)  = 10._wp*ACOS(z_fac1(je))
             zu(je) = zu(je) + jw_up* EXP(-z_fac2(je)**2)
            END IF
            zv(je) = 0._wp
            p_nhdom%prog(jt)%vn(je,jk,jb) = &
                       zu(je) * p_patch(jg)%edges%primal_normal(je,jb)%v1   &
                   & + zv(je) * p_patch(jg)%edges%primal_normal(je,jb)%v2
          ENDDO
        ENDDO
      ENDDO !jt
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


    IF (iforcing == inwp .AND. (inwp_gscp /= 0 .OR. inwp_convection /= 0)) THEN
      niter = 10 ! experimentation suggest at least 6 iterations
    ELSE
      niter = 1
    ENDIF
    ! Do some iterations to come closer to the moisture and temperature/pressure 
    DO ji = 1, niter

      IF ( ltransport ) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jt,jjt,jc,nlen,zrhf,z_1_o_rh,z_help,zsqv,ctracer,zlat,zlon)
        DO jb = 1, nblks_c
          IF (jb /= nblks_c) THEN 
            nlen = nproma
          ELSE
            nlen = npromz_c
          ENDIF

          DO jt = 1, ntl
            DO jk = 1, nlev
              DO jjt = 1, ntracer
                IF ( iforcing == inwp  ) THEN
                  IF(jjt == iqv .AND. (inwp_gscp /= 0 .OR. inwp_convection /= 0)) THEN
                    DO jc =1, nlen

!                      !KF linear decreasing RH with height like in Hui's testcase
                      zrhf     = rh_at_1000hpa-0.5_wp+p_nhdom%diag%pres(jc,jk,jb)/200000._wp
                      zrhf     = MAX (zrhf,0.0_wp)
                      z_1_o_rh = 1._wp/(zrhf+1.e-6_wp)
                      ! to avoid water vapor pressure > total pressure:
                      z_help = MIN ( sat_pres_water( p_nhdom%diag%temp(jc,jk,jb) ), &
                         & p_nhdom%diag%pres(jc,jk,jb) * z_1_o_rh )
                      IF( p_nhdom%diag%temp(jc,jk,jb) <= tmelt) THEN
                        ! to avoid water vapor pressure > total pressure:
                        z_help = MIN ( sat_pres_ice( p_nhdom%diag%temp(jc,jk,jb) ), &
                         & p_nhdom%diag%pres(jc,jk,jb) * z_1_o_rh )
                      ENDIF
                      ! saturation qv calculated as in mo_satad's qsat_rho
                      zsqv = z_help / ( p_nh_state(jg)%prog(jt)%rho(jc,jk,jb) * rv  &
                         &      * p_nhdom%diag%temp(jc,jk,jb) )
                      p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) = MIN ( zsqv,  zrhf*zsqv )

                      IF (  p_nhdom%diag%pres(jc,jk,jb) <= 10000._wp) &
                         &  p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt)      &
                         &  = MIN ( 5.e-6_wp, p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) )
                    
                      ! Limit QV in the tropics                       
                      p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) = &
                      &   MIN(qv_max,p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt))
  
                    END DO
  
                  ELSE !other tracers than water vapor zero at start
                    p_nhdom%prog(jt)%tracer(:,jk,jb,jjt) = 0._wp
                  ENDIF ! tracer

                ELSE !iforcing /= inwp:

                  ctracer = ctracer_list(jjt:jjt)
  
                  SELECT CASE(ctracer)
  
                  CASE('1')
  
                    DO jc = 1, nlen
                      zlat = p_patch(jg)%cells%center(jc,jb)%lat
                      zlon = p_patch(jg)%cells%center(jc,jb)%lon
                      p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) =  &
                      tracer_q1_q2(zlon, zlat, zeta_v(jc,jk,jb), rotate_axis_deg, 0.6_wp)
                    ENDDO ! cell loop
  
                  CASE('2')
  
                    DO jc =1, nlen
                      zlat = p_patch(jg)%cells%center(jc,jb)%lat
                      zlon = p_patch(jg)%cells%center(jc,jb)%lon
                      p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) =  &
                        tracer_q1_q2(zlon, zlat, zeta_v(jc,jk,jb), rotate_axis_deg, 1.0_wp)
                    ENDDO ! cell loop
  
                  CASE('3')
  
                    DO jc =1, nlen
                      zlat = p_patch(jg)%cells%center(jc,jb)%lat
                      zlon = p_patch(jg)%cells%center(jc,jb)%lon
                      p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) =  &
                        tracer_q3(zlon, zlat, rotate_axis_deg)
                    ENDDO ! cell loop
  
                  CASE('4')
                    p_nhdom%prog(jt)%tracer(:,jk,jb,jjt) = 1._wp
  
                  END SELECT

                ENDIF   ! iforcing                       
               
              ENDDO ! tracer loop
            ENDDO ! vertical level loop
          ENDDO   ! time level loop
        ENDDO ! block loop

!$OMP END DO
!$OMP END PARALLEL
       
      ENDIF !ltransport

      IF ( iforcing == inwp ) THEN 
        
        CALL diagnose_pres_temp ( p_nhdom%metrics, p_nhdom%prog(nnow(jg)),    &
          &                       p_nhdom%prog(nnow_rcf(jg)), p_nhdom%diag,   &
          &                       p_patch(jg),                                &
          &                       opt_calc_temp=.TRUE.,                       &
          &                       opt_calc_pres=.TRUE.                      )
       
        CALL diagnose_pres_temp ( p_nhdom%metrics, p_nhdom%prog(nnew(jg)),    &
          &                       p_nhdom%prog(nnew_rcf(jg)), p_nhdom%diag,   &
          &                       p_patch(jg),                                &
          &                       opt_calc_temp=.TRUE.,                       &
          &                       opt_calc_pres=.TRUE.                       )        
      ENDIF

    ENDDO ! ji


    DEALLOCATE (zeta_v,zeta_v_e)

  ENDDO !jg

  CALL message(TRIM(routine),'End setup Jablonowski test')
  
  CASE ('mrw_nh', 'mrw2_nh')
 
  bruntvaissq = grav*grav/cpd/temp_mrw
  kappa       = rd/cpd
  zhelp1      = bruntvaissq/grav/grav/kappa
  zhelp2      = grav/rd/temp_mrw

  DO jg = 1, n_dom

    ALLOCATE (z_wsfc_c(nproma,1,p_patch(jg)%nblks_c), &
              z_wsfc_e(nproma,1,p_patch(jg)%nblks_e) )

    nblks_c   = p_patch(jg)%nblks_int_c
    npromz_c  = p_patch(jg)%npromz_int_c
    nblks_e   = p_patch(jg)%nblks_int_e
    npromz_e  = p_patch(jg)%npromz_int_e

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    p_nhdom  => p_nh_state(jg)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
!     set the surface pressure
      DO jc = 1, nlen
        zlat= p_patch(jg)%cells%center(jc,jb)%lat
        zcoslat=COS(zlat)
        p_nhdom%diag%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1 * ( u0_mrw *&
                              ( 0.5_wp*u0_mrw + re*omega) * zcoslat*zcoslat - &
                              ext_data(jg)%atm%topography_c(jc,jb)*grav))
      ENDDO !jc
      DO jk = nlev, 1, -1
        IF (jk == nlev) THEN
          ! Use analytic expressions at lowest model level
          DO jt = 1, ntl
            DO jc = 1, nlen
              z_sfc  = ext_data(jg)%atm%topography_c(jc,jb)
              z_nlev = 0.5_wp * ( p_nhdom%metrics%z_ifc(jc,nlev,jb) + &
                       p_nhdom%metrics%z_ifc(jc,nlev+1,jb) )
              z_pres = p_nhdom%diag%pres_sfc(jc,jb) * &
                       EXP(- zhelp2 * ( z_nlev  - z_sfc )  )   !isothermal atm.

              p_nhdom%prog(jt)%exner(jc,jk,jb) = (z_pres/p0ref)**(rd/cpd)
              p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                       p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
              p_nhdom%prog(jt)%theta_v(jc,jk,jb) = temp_mrw & 
                       /p_nhdom%prog(jt)%exner(jc,jk,jb)
              p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                       p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                       /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
              z_temp_kp1(jc) = temp_mrw
              z_temp(jc) = temp_mrw
              IF (jt == 1) THEN
                 p_nhdom%diag%temp(jc,jk,jb) = temp_mrw    
                 p_nhdom%diag%pres(jc,jk,jb) = z_pres
              ENDIF
            ENDDO !jc
          ENDDO !jt
        ELSE
          ! Solve quadratic equation for exner(jk) to obtain exact (discretized)
          ! hydrostatic balance
          IF ( p_patch(jg)%cell_type == 3) THEN
            DO jc = 1, nlen
              z_fac1(jc) = p_nhdom%metrics%wgtfac_c(jc,jk+1,jb)*(z_temp_kp1(jc)              &
                - p_nhdom%metrics%theta_ref_mc(jc,jk+1,jb)*p_nhdom%prog(1)%exner(jc,jk+1,jb))&
                - (1._wp-p_nhdom%metrics%wgtfac_c(jc,jk+1,jb))                               &
                * p_nhdom%metrics%theta_ref_mc(jc,jk,jb)*p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac2(jc) = (1._wp-p_nhdom%metrics%wgtfac_c(jc,jk+1,jb))*z_temp(jc) &
                *p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac3(jc) = p_nhdom%metrics%exner_ref_mc(jc,jk+1,jb)                    &
                -p_nhdom%metrics%exner_ref_mc(jc,jk,jb)-p_nhdom%prog(1)%exner(jc,jk+1,jb)
              za(jc) = ( p_nhdom%metrics%theta_ref_ic(jc,jk+1,jb) &
                *p_nhdom%prog(1)%exner(jc,jk+1,jb)+z_fac1(jc))     &
                /p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)
              zb(jc) = -(za(jc)*z_fac3(jc)+z_fac2(jc)/p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb) &
                + z_fac1(jc)*p_nhdom%metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
              zc(jc) = -(z_fac2(jc)*z_fac3(jc)/p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb) &
                + z_fac2(jc)*p_nhdom%metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
            ENDDO !jc
          ELSEIF ( p_patch(jg)%cell_type == 6) THEN
            DO jc = 1, nlen
              z_fac3(jc) = grav/cpd*p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)
              z_fac1(jc) = p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac2(jc) = p_nhdom%prog(1)%exner(jc,jk+1,jb)/z_temp_kp1(jc)
              za(jc) = 1.0_wp
              zb(jc) = -(z_fac2(jc)*(z_temp(jc)+2.0_wp*z_fac3(jc))-z_fac1(jc))
              zc(jc) = z_temp(jc)*z_fac2(jc)*z_fac1(jc)
            ENDDO !jc
          ENDIF
          DO jt = 1, ntl
            DO jc = 1, nlen
              p_nhdom%prog(jt)%exner(jc,jk,jb) = &
                (zb(jc) + SQRT(zb(jc)**2+4._wp*za(jc)*zc(jc)))/(2._wp*za(jc))
              p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                       p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
              p_nhdom%prog(jt)%theta_v(jc,jk,jb) = z_temp(jc) & 
                       /p_nhdom%prog(jt)%exner(jc,jk,jb)
              p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                       p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                       /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
              z_temp_kp1(jc) = z_temp(jc)
              IF (jt == 1) THEN  
                 p_nhdom%diag%temp(jc,jk,jb) = temp_mrw       
                 p_nhdom%diag%pres(jc,jk,jb) = p0ref*p_nhdom%prog(jt)%exner(jc,jk,jb)**(cpd/rd)  
              ENDIF
            ENDDO !jc
          ENDDO !jt
        ENDIF
      ENDDO !jk
    ENDDO !jb

    i_startblk = p_patch(jg)%edges%start_blk(2,1)
    ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,jt,zlat,z_u)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      DO jt = 1, ntl
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            zlat = p_patch(jg)%edges%center(je,jb)%lat
            z_u = u0_mrw * COS(zlat)  !v component is zero
            p_nhdom%prog(jt)%vn(je,jk,jb) = &
             z_u * p_patch(jg)%edges%primal_normal(je,jb)%v1
          ENDDO
        ENDDO
        DO je = i_startidx, i_endidx
          z_wsfc_e(je,1,jb) = p_nhdom%prog(jt)%vn(je,nlev,jb) * &
             p_nhdom%metrics%ddxn_z_half(je,nlevp1,jb)
        ENDDO !je
      ENDDO !jt
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    CALL edges2cells_scalar(z_wsfc_e,p_patch(jg),p_int(jg)%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)
    i_startblk = p_patch(jg)%cells%start_blk(2,1)

    ! specify a reasonable initial vertical wind speed
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jt,jk)
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch(jg), jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      DO jt = 1, ntl
        DO jc = i_startidx, i_endidx
          p_nhdom%prog(jt)%w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
        ENDDO
        DO jk = nlev, 2, -1
          DO jc = i_startidx, i_endidx
            p_nhdom%prog(jt)%w(jc,jk,jb) = z_wsfc_c(jc,1,jb)*vct_b(jk)**2
          ENDDO !jc
        ENDDO !jk
      ENDDO !jt
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL


!!!!!!! Test !!!!!!!!!

    ! tracer for physics test
    
     IF ( iforcing == inwp ) THEN 
        
      CALL diagnose_pres_temp ( p_nhdom%metrics, p_nhdom%prog(nnow(jg)),       &
        &                       p_nhdom%prog(nnow_rcf(jg)), p_nhdom%diag,      &
        &                       p_patch(jg),                                   &
        &                       opt_calc_temp=.TRUE.,                          &
        &                       opt_calc_pres=.TRUE.                          )
       
      CALL diagnose_pres_temp ( p_nhdom%metrics, p_nhdom%prog(nnew(jg)),       &
        &                       p_nhdom%prog(nnew_rcf(jg)), p_nhdom%diag,      &
        &                       p_patch(jg),                                   &
        &                       opt_calc_temp=.TRUE.,                          &
        &                       opt_calc_pres=.TRUE.                          )
        

!!$OMP DO PRIVATE(jb,jk,jt,jc,nlen,zeta,ctracer,lon,lat,zrhf,z_1_o_rh,zsqv,zhelp)
        DO jb = 1, nblks_c
           IF (jb /= nblks_c) THEN 
              nlen = nproma
           ELSE
              nlen = npromz_c
           ENDIF

        DO jt = 1, ntl

           DO jk = 1, nlev

              DO jjt = 1, ntracer

                  IF(jjt == iqv .AND. (inwp_gscp /= 0 .OR. inwp_convection /= 0) ) THEN
                     
                     DO jc =1, nlen

!KF linear decreasing RH with height like in Hui's testcase
                       zrhf =  rh_at_1000hpa - 0.5_wp +p_nhdom%diag%pres(jc,jk,jb)/200000._wp
                       zrhf    = MAX (zrhf,0.0_wp)
                       z_1_o_rh = 1._wp/(zrhf+1.e-6_wp)
                       ! to avoid water vapor pressure > total pressure:
                       z_help = MIN ( sat_pres_water( p_nhdom%diag%temp(jc,jk,jb) ), &
                         & p_nhdom%diag%pres(jc,jk,jb) * z_1_o_rh )
                       IF( p_nhdom%diag%temp(jc,jk,jb) <= tmelt) THEN
                         ! to avoid water vapor pressure > total pressure:
                         z_help = MIN ( sat_pres_ice( p_nhdom%diag%temp(jc,jk,jb) ), &
                           & p_nhdom%diag%pres(jc,jk,jb) * z_1_o_rh )
                       ENDIF
                       ! saturation qv calculated as in mo_satad's qsat_rho
                       zsqv = z_help / ( p_nh_state(jg)%prog(jt)%rho(jc,jk,jb) * rv  &
                         &      * p_nhdom%diag%temp(jc,jk,jb) )
                       p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) = MIN ( zsqv,  zrhf*zsqv )
                       IF (  p_nhdom%diag%pres(jc,jk,jb) <= 10000._wp) &
                         &  p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt)      &
                         &  = MIN ( 5.e-6_wp, p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) )
                       
                       ! Limit QV in the tropics
                       p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt) = &
                            &   MIN(qv_max,p_nhdom%prog(jt)%tracer(jc,jk,jb,jjt))
                       
                     END DO

                  ELSE !other tracers than water vapor zero at start
                    p_nhdom%prog(jt)%tracer(:,jk,jb,jjt) = 0._wp
                  ENDIF ! tracer
                 
              ENDDO ! tracer loop

           ENDDO ! vertical level loop
           
         ENDDO   ! time level loop

       ENDDO ! block loop
!!$OMP END DO

    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
!!$!     set the surface pressure
!!$      DO jc = 1, nlen
!!$        zlat= p_patch(jg)%cells%center(jc,jb)%lat
!!$        zcoslat=COS(zlat)
!!$        p_nhdom%diag%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1 * ( u0_mrw *&
!!$                              ( 0.5_wp*u0_mrw + re*omega) * zcoslat*zcoslat - &
!!$                              ext_data(jg)%atm%topography_c(jc,jb)*grav))
!!$      ENDDO !jc
      DO jk = nlev, 1, -1
        IF (jk == nlev) THEN
          ! Use analytic expressions at lowest model level
          DO jt = 1, ntl
            DO jc = 1, nlen
              z_sfc  = ext_data(jg)%atm%topography_c(jc,jb)
              z_nlev = 0.5_wp * ( p_nhdom%metrics%z_ifc(jc,nlev,jb) + &
                       p_nhdom%metrics%z_ifc(jc,nlev+1,jb) )
              z_pres = p_nhdom%diag%pres_sfc(jc,jb) * &
                       EXP(- zhelp2 * ( z_nlev  - z_sfc )  )   !isothermal atm.

!              p_nhdom%prog(jt)%exner(jc,jk,jb) = (z_pres/p0ref)**(rd/cpd)
              p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                       p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
              p_nhdom%prog(jt)%theta_v(jc,jk,jb) = temp_mrw &
              &       *( 1._wp +  vtmpc1            &
              &       *p_nhdom%prog(jt)%tracer(jc,jk,jb,iqv) ) &
              &         /p_nhdom%prog(jt)%exner(jc,jk,jb)
              p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                       p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                       /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
              z_temp_kp1(jc) = temp_mrw             &
              &       *( 1._wp +  vtmpc1            &
              &       *p_nhdom%prog(jt)%tracer(jc,jk,jb,iqv) )            

              
              z_temp(jc) = temp_mrw                 &
              &       *( 1._wp +  vtmpc1            &
              &       *p_nhdom%prog(jt)%tracer(jc,jk,jb,iqv) )            

              
              IF (jt == 1) THEN
                 p_nhdom%diag%temp(jc,jk,jb) = temp_mrw    
                 p_nhdom%diag%pres(jc,jk,jb) = z_pres
              ENDIF
            ENDDO !jc
          ENDDO !jt
        ELSE
          DO jc = 1, nlen
            z_temp(jc) = temp_mrw  &
              &       *( 1._wp +  vtmpc1            &
              &       *p_nhdom%prog(nnew(jg))%tracer(jc,jk,jb,iqv) )
          ENDDO
          
          ! Solve quadratic equation for exner(jk) to obtain exact (discretized)
          ! hydrostatic balance
          IF ( p_patch(jg)%cell_type == 3) THEN
            DO jc = 1, nlen
              z_fac1(jc) = p_nhdom%metrics%wgtfac_c(jc,jk+1,jb)*(z_temp_kp1(jc)              &
                - p_nhdom%metrics%theta_ref_mc(jc,jk+1,jb)*p_nhdom%prog(1)%exner(jc,jk+1,jb))&
                - (1._wp-p_nhdom%metrics%wgtfac_c(jc,jk+1,jb))                               &
                * p_nhdom%metrics%theta_ref_mc(jc,jk,jb)*p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac2(jc) = (1._wp-p_nhdom%metrics%wgtfac_c(jc,jk+1,jb))*z_temp(jc) &
                *p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac3(jc) = p_nhdom%metrics%exner_ref_mc(jc,jk+1,jb)                    &
                -p_nhdom%metrics%exner_ref_mc(jc,jk,jb)-p_nhdom%prog(1)%exner(jc,jk+1,jb)
              za(jc) = ( p_nhdom%metrics%theta_ref_ic(jc,jk+1,jb) &
                *p_nhdom%prog(1)%exner(jc,jk+1,jb)+z_fac1(jc))     &
                /p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)
              zb(jc) = -(za(jc)*z_fac3(jc)+z_fac2(jc)/p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb) &
                + z_fac1(jc)*p_nhdom%metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
              zc(jc) = -(z_fac2(jc)*z_fac3(jc)/p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb) &
                + z_fac2(jc)*p_nhdom%metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
            ENDDO !jc
          ELSEIF ( p_patch(jg)%cell_type == 6) THEN
            DO jc = 1, nlen
              z_fac3(jc) = grav/cpd*p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)
              z_fac1(jc) = p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac2(jc) = p_nhdom%prog(1)%exner(jc,jk+1,jb)/z_temp_kp1(jc)
              za(jc) = 1.0_wp
              zb(jc) = -(z_fac2(jc)*(z_temp(jc)+2.0_wp*z_fac3(jc))-z_fac1(jc))
              zc(jc) = z_temp(jc)*z_fac2(jc)*z_fac1(jc)
            ENDDO !jc
          ENDIF
          DO jt = 1, ntl
            DO jc = 1, nlen
              p_nhdom%prog(jt)%exner(jc,jk,jb) = &
                (zb(jc) + SQRT(zb(jc)**2+4._wp*za(jc)*zc(jc)))/(2._wp*za(jc))
              p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                       p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
              p_nhdom%prog(jt)%theta_v(jc,jk,jb) = z_temp(jc) & 
                       /p_nhdom%prog(jt)%exner(jc,jk,jb)
              p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                       p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                       /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
              z_temp_kp1(jc) = z_temp(jc)
              IF (jt == 1) THEN  
                 p_nhdom%diag%temp(jc,jk,jb) = temp_mrw       
                 p_nhdom%diag%pres(jc,jk,jb) = p0ref*p_nhdom%prog(jt)%exner(jc,jk,jb)**(cpd/rd)  
              ENDIF
            ENDDO !jc
          ENDDO !jt
        ENDIF
      ENDDO !jk
    ENDDO !jb

       
    ENDIF !inwp

    DEALLOCATE (z_wsfc_c, z_wsfc_e )

  ENDDO !jg

 CALL message(TRIM(routine),'End setup MRW test')

  CASE ('mwbr_const')
 
  bruntvaissq_i = grav*grav/cpd/temp_i_mwbr_const
  bruntvaissq_u = bruntvais_u_mwbr_const * bruntvais_u_mwbr_const
  kappa       = rd/cpd
  zhelp1_i     = bruntvaissq_i/grav/grav/kappa
  zhelp2_i     = grav/rd/temp_i_mwbr_const
  zhelp1_u     = bruntvaissq_u/grav/grav/kappa
  zhelp3       = u0_mwbr_const/re + 2.0_wp*omega
  theta_v_int  = temp_i_mwbr_const * (p0ref/p_int_mwbr_const)**kappa
  rkappa       = 1.0_wp/kappa

  DO jg = 1, n_dom

    ALLOCATE (z_wsfc_c(nproma,1,p_patch(jg)%nblks_c), &
              z_wsfc_e(nproma,1,p_patch(jg)%nblks_e) )
    ALLOCATE (z_int_c(nproma,p_patch(jg)%nblks_c))
    

    nblks_c   = p_patch(jg)%nblks_int_c
    npromz_c  = p_patch(jg)%npromz_int_c
    nblks_e   = p_patch(jg)%nblks_int_e
    npromz_e  = p_patch(jg)%npromz_int_e

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    p_nhdom  => p_nh_state(jg)


    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
!     set z_int_c, z at the interface
      icount = 0
      DO jc = 1, nlen
        zlat= p_patch(jg)%cells%center(jc,jb)%lat
        zcoslat=COS(zlat)

        z_sfc  = ext_data(jg)%atm%topography_c(jc,jb)
        z_int_c(jc,jb) = LOG(pres_sp/p_int_mwbr_const)/grav/zhelp1_i + &
                         zcoslat*zcoslat*zhelp3*re*u0_mwbr_const/2.0_wp/grav
        IF (z_int_c(jc,jb) <  0._wp ) icount = icount + 1
        IF (z_int_c(jc,jb) >= z_sfc ) THEN

          p_nhdom%diag%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1_i * ( u0_mwbr_const * &
                              zhelp3*re * zcoslat*zcoslat/2.0_wp - &
                              z_sfc*grav))
        ELSE
          zhelp4 = z_sfc - z_int_c(jc,jb)
          p_nhdom%diag%pres_sfc(jc,jb) = p_int_mwbr_const * ( 1.0_wp +       &
                              (EXP(-bruntvaissq_u*zhelp4/grav) - 1.0_wp) *   &
                               bruntvaissq_i/bruntvaissq_u )**rkappa
        ENDIF
      ENDDO !jc
      IF (icount > 0) CALL finish(TRIM(routine), &
        & 'z at interface is negative, needed p_int_mwbr_const < pres_sp')

      DO jk = nlev, 1, -1
        IF (jk == nlev) THEN
          ! Use analytic expressions at lowest model level
          DO jt = 1, ntl
            DO jc = 1, nlen
              z_sfc  = ext_data(jg)%atm%topography_c(jc,jb)
              z_nlev = 0.5_wp * ( p_nhdom%metrics%z_ifc(jc,nlev,jb) + &
                       p_nhdom%metrics%z_ifc(jc,nlev+1,jb) )
              IF (z_nlev < z_int_c(jc,jb) ) THEN

                z_pres = p_nhdom%diag%pres_sfc(jc,jb) * &
                         EXP(- zhelp2_i * ( z_nlev  - z_sfc )  )   !isothermal atm.

                p_nhdom%prog(jt)%exner(jc,jk,jb) = (z_pres/p0ref)**(rd/cpd)
                p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                         p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
                p_nhdom%prog(jt)%theta_v(jc,jk,jb) = temp_mrw & 
                         /p_nhdom%prog(jt)%exner(jc,jk,jb)
                p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                         p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                         /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
                z_temp_kp1(jc) = temp_i_mwbr_const
                z_temp(jc) = temp_i_mwbr_const
                p_nhdom%diag%temp(jc,jk,jb) = temp_i_mwbr_const    
                p_nhdom%diag%pres(jc,jk,jb) = z_pres
               ELSEIF (z_nlev > z_int_c(jc,jb) ) THEN
                zhelp4 = z_nlev - z_int_c(jc,jb)
                z_pres = p_int_mwbr_const * ( 1.0_wp +       &
                              (EXP(-bruntvaissq_u*zhelp4/grav) - 1.0_wp) *   &
                               bruntvaissq_i/bruntvaissq_u )**rkappa
                z_temp(jc) = temp_i_mwbr_const * EXP (bruntvaissq_u*zhelp4/grav) * &
                             (z_pres/p_int_mwbr_const)**kappa
                p_nhdom%prog(jt)%exner(jc,jk,jb) = (z_pres/p0ref)**(rd/cpd)
                p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                         p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
                p_nhdom%prog(jt)%theta_v(jc,jk,jb) =   theta_v_int * & 
                                                       EXP(bruntvaissq_u*zhelp4/grav)
                p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                         p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                         /p_nhdom%prog(jt)%theta_v(jc,jk,jb) 
                p_nhdom%diag%temp(jc,jk,jb) =  z_temp(jc)   
                p_nhdom%diag%pres(jc,jk,jb) = z_pres 
                z_temp_kp1(jc) = z_temp(jc)
              ENDIF

            ENDDO !jc
          ENDDO !jt
        ELSE
          DO jc = 1, nlen
            z_nlev = 0.5_wp * ( p_nhdom%metrics%z_ifc(jc,jk,jb) + &
                       p_nhdom%metrics%z_ifc(jc,jk+1,jb) )
            IF (z_nlev <= z_int_c(jc,jb) ) THEN
                z_temp(jc) = temp_i_mwbr_const
                    
            ELSE
                zhelp4 = z_nlev - z_int_c(jc,jb)
                z_pres = p_int_mwbr_const * ( 1.0_wp +       &
                              (EXP(-bruntvaissq_u*zhelp4/grav) - 1.0_wp) *   &
                               bruntvaissq_i/bruntvaissq_u )**rkappa
                z_temp(jc) = temp_i_mwbr_const * EXP (bruntvaissq_u*zhelp4/grav) * &
                             (z_pres/p_int_mwbr_const)**kappa 

            ENDIF
            p_nhdom%diag%temp(jc,jk,jb) =  z_temp(jc)
          ENDDO
          ! Solve quadratic equation for exner(jk) to obtain exact (discretized)
          ! hydrostatic balance
          IF( p_patch(jg)% == 3) THEN
            DO jc = 1, nlen
              z_fac1(jc) = p_nhdom%metrics%wgtfac_c(jc,jk+1,jb)*(z_temp_kp1(jc)              &
                - p_nhdom%metrics%theta_ref_mc(jc,jk+1,jb)*p_nhdom%prog(1)%exner(jc,jk+1,jb))&
                - (1._wp-p_nhdom%metrics%wgtfac_c(jc,jk+1,jb))                               &
                * p_nhdom%metrics%theta_ref_mc(jc,jk,jb)*p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac2(jc) = (1._wp-p_nhdom%metrics%wgtfac_c(jc,jk+1,jb))*z_temp(jc) &
                *p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac3(jc) = p_nhdom%metrics%exner_ref_mc(jc,jk+1,jb)                    &
                -p_nhdom%metrics%exner_ref_mc(jc,jk,jb)-p_nhdom%prog(1)%exner(jc,jk+1,jb)
              za(jc) = ( p_nhdom%metrics%theta_ref_ic(jc,jk+1,jb) &
                *p_nhdom%prog(1)%exner(jc,jk+1,jb)+z_fac1(jc))     &
                /p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)
              zb(jc) = -(za(jc)*z_fac3(jc)+z_fac2(jc)/p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb) &
                + z_fac1(jc)*p_nhdom%metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
              zc(jc) = -(z_fac2(jc)*z_fac3(jc)/p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb) &
                + z_fac2(jc)*p_nhdom%metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
            ENDDO !jc
          ELSEIF ( p_patch(jg)% == 6) THEN
            DO jc = 1, nlen
              z_fac3(jc) = grav/cpd*p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)
              z_fac1(jc) = p_nhdom%prog(1)%exner(jc,jk+1,jb)
              z_fac2(jc) = p_nhdom%prog(1)%exner(jc,jk+1,jb)/z_temp_kp1(jc)
              za(jc) = 1.0_wp
              zb(jc) = -(z_fac2(jc)*(z_temp(jc)+2.0_wp*z_fac3(jc))-z_fac1(jc))
              zc(jc) = z_temp(jc)*z_fac2(jc)*z_fac1(jc)
            ENDDO !jc
          ENDIF
          DO jt = 1, ntl
            DO jc = 1, nlen
              p_nhdom%prog(jt)%exner(jc,jk,jb) = &
                (zb(jc) + SQRT(zb(jc)**2+4._wp*za(jc)*zc(jc)))/(2._wp*za(jc))
              p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
                       p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
              p_nhdom%prog(jt)%theta_v(jc,jk,jb) = z_temp(jc) & 
                       /p_nhdom%prog(jt)%exner(jc,jk,jb)
              p_nhdom%prog(jt)%rho(jc,jk,jb) = &
                       p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
                       /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
              z_temp_kp1(jc) = z_temp(jc)
              IF (jt == 1) THEN
                  p_nhdom%diag%pres(jc,jk,jb) = p0ref*p_nhdom%prog(jt)%exner(jc,jk,jb)**(cpd/rd)
              ENDIF
            ENDDO !jc
          ENDDO !jt
        ENDIF
      ENDDO !jk
    ENDDO !jb
    i_startblk = p_patch(jg)%edges%start_blk(2,1)
    ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,jt,zlat,z_u)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      DO jt = 1, ntl
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            zlat = p_patch(jg)%edges%center(je,jb)%lat
            z_u = u0_mwbr_const * COS(zlat)  !v component is zero
            p_nhdom%prog(jt)%vn(je,jk,jb) = &
             z_u * p_patch(jg)%edges%primal_normal(je,jb)%v1
          ENDDO
        ENDDO
        DO je = i_startidx, i_endidx
          z_wsfc_e(je,1,jb) = p_nhdom%prog(jt)%vn(je,nlev,jb) * &
             p_nhdom%metrics%ddxn_z_half(je,nlevp1,jb)
        ENDDO !je
      ENDDO !jt
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    CALL edges2cells_scalar(z_wsfc_e,p_patch(jg),p_int(jg)%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)
    i_startblk = p_patch(jg)%cells%start_blk(2,1)

    ! specify a reasonable initial vertical wind speed
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jt,jk)
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch(jg), jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      DO jt = 1, ntl
        DO jc = i_startidx, i_endidx
          p_nhdom%prog(jt)%w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
        ENDDO
        DO jk = nlev, 2, -1
          DO jc = i_startidx, i_endidx
            p_nhdom%prog(jt)%w(jc,jk,jb) = z_wsfc_c(jc,1,jb)*vct_b(jk)**2
          ENDDO !jc
        ENDDO !jk
      ENDDO !jt
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    DEALLOCATE (z_int_c)
    DEALLOCATE (z_wsfc_c, z_wsfc_e )
  ENDDO !jg

  CASE ('zero','bell','schaer')

  ! For the moment we think of a given Brunt Vaisala frequency and a given      
  ! zonal wind. The lplane and the lcorio=F options are assumed

  DO jg = 1, n_dom
    p_nhdom   => p_nh_state(jg)
    nblks_e   = p_patch(jg)%nblks_int_e
    npromz_e  = p_patch(jg)%npromz_int_e

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    DO jt = 1, ntl 
      ! normal wind
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF
        DO jk = 1, nlev
          DO je = 1, nlen
            p_nh_state(jg)%prog(jt)%vn(je,jk,jb) = nh_u0 &
            !(p_nhdom%metrics%geopot(1,jk,1)/grav/1000.0_wp+5.0_wp)& !shear
            *p_patch(jg)%edges%primal_normal(je,jb)%v1
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    nblks_c   = p_patch(jg)%nblks_int_c
    npromz_c  = p_patch(jg)%npromz_int_c
    ! scalars (all is dry!)
    DO jt = 1, ntl 
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF
        DO jk = 1, nlev
          DO jc = 1, nlen
            z_help=(nh_brunt_vais/grav)**2*p_nhdom%metrics%geopot(jc,jk,jb)
            ! profile of theta is explicitly given
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) = nh_t0*EXP(z_help)
          ENDDO
        ENDDO
        DO jk = nlev, 1, -1
          DO jc = 1, nlen
            IF (jk == nlev) THEN
              IF (nh_brunt_vais /= 0.0_wp) THEN
                ! for exner(nlev) (lowermost level): analytical exner
                z_help = (nh_brunt_vais/grav)**2*p_nhdom%metrics%geopot(jc,jk,jb)
                p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) =    &
                  (grav/nh_brunt_vais)**2/nh_t0/cpd*(EXP(-z_help)-1.0_wp)+1.0_wp
              ELSE
                p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) = &
                  1.0_wp-p_nhdom%metrics%geopot(jc,jk,jb)/cpd/nh_t0
              ENDIF
            ELSE ! other levels are hydrostatically balanced with respect to model numerics
              z_help=0.5_wp*(p_nh_state(jg)%prog(jt)%theta_v(jc,jk  ,jb) &
                           + p_nh_state(jg)%prog(jt)%theta_v(jc,jk+1,jb))
              p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) = &
              & p_nh_state(jg)%prog(jt)%exner(jc,jk+1,jb) &
              & -grav/cpd*p_nhdom%metrics%ddqz_z_half(jc,jk+1,jb)/z_help
            ENDIF
          ENDDO
        ENDDO
        DO jk = 1, nlev
          DO jc = 1, nlen
            ! rhotheta has to have the same meaning as exner
            p_nh_state(jg)%prog(jt)%rhotheta_v(jc,jk,jb) = &
                (p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)**cvd_o_rd)*p0ref/rd

       !     ! perturbation in theta_v for gravity test case
       !     p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=&
       !              p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)+ 0.01_wp*sin(&
       !              pi*p_nhdom%metrics%geopot(jc,jk,jb)/grav/10000.0_wp)&
       !            /(1.0_wp+(p_patch(jg)%cells%center(jc,jb)%lon*30.0/pi)**2)
       !            !Das ist fuer dx=500 mit 600 Punkten (auf 2pi verteilt)

       !     ! perturbation in theta_v for Straka test case
       !     z_help = SQRT( ( p_patch(jg)%cells%center(jc,jb)%lon*6.4_wp/pi)**2 &
       !     &      +((p_nhdom%metrics%z_mc(jc,jk,jb)-3000.0_wp)/2000.0_wp)**2)
       !     IF (z_help<=1.0_wp) THEN
       !       p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=&
       !       &   p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) &
       !       &   -15.0_wp*(COS(pi*z_help)+1.0_wp)*0.5_wp&
       !       &   /p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)
       !     ENDIF

            ! exner, rhotheta and theta_v are given, so rho is deduced...
            p_nh_state(jg)%prog(jt)%rho(jc,jk,jb) = &
            &        p_nh_state(jg)%prog(jt)%rhotheta_v(jc,jk,jb) &
            &       /p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)

          ENDDO
        ENDDO
        DO jk = 1, nlevp1
          p_nh_state(jg)%prog(jt)%w(1:nlen,jk,jb) = 0.0_wp
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  CASE ('PA')  ! pure advection test case, no mountain

    DO jg = 1, n_dom
      CALL init_nh_state_prog_patest(p_patch(jg),p_int(jg),             &
        &                       p_nh_state(jg)%prog(nnow(jg)),          &
        &                       p_nh_state(jg)%diag,ext_data(jg),       &
        &                       p_nh_state(jg)%metrics,rotate_axis_deg, &
        &                       linit_tracer_fv )

      CALL init_nh_state_prog_patest(p_patch(jg),p_int(jg),             &
        &                       p_nh_state(jg)%prog(nnew(jg)),          &
        &                       p_nh_state(jg)%diag,ext_data(jg),       &
        &                       p_nh_state(jg)%metrics,rotate_axis_deg, &
        &                       linit_tracer_fv )

    ENDDO !jg

  CASE ('DF1', 'DF2', 'DF3', 'DF4')  ! 2D deformational flow test case, no mountain

    DO jg = 1, n_dom

      CALL init_nh_state_prog_dftest(p_patch(jg),                    &
           &                         p_nh_state(jg)%prog(nnow(jg)),  &
           &                         p_nh_state(jg)%diag,            &
           &                         p_int(jg), ext_data(jg),        &
           &                         p_nh_state(jg)%metrics,         &
           &                         rotate_axis_deg, nh_test_name,  &
           &                         linit_tracer_fv )

      CALL init_nh_state_prog_dftest(p_patch(jg),                    &
           &                         p_nh_state(jg)%prog(nnew(jg)),  &
           &                         p_nh_state(jg)%diag,            &
           &                         p_int(jg), ext_data(jg),        &
           &                         p_nh_state(jg)%metrics,         &
           &                         rotate_axis_deg, nh_test_name,  &
           &                         linit_tracer_fv )

    ENDDO !jg

  CASE ('HS_nh')  ! Held-Suarez test case, no mountain, isothermal atmosphere

    DO jg = 1, n_dom

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev

      CALL init_nh_state_prog_held_suarez(p_patch(jg),p_nh_state(jg)%prog(nnow(jg)),  &
        &                            p_nh_state(jg)%diag,ext_data(jg),           &
        &                            p_nh_state(jg)%metrics)

      CALL init_nh_state_prog_held_suarez(p_patch(jg),p_nh_state(jg)%prog(nnew(jg)),  &
        &                            p_nh_state(jg)%diag,ext_data(jg),           &
        &                            p_nh_state(jg)%metrics)

      IF (lhs_nh_vn_ptb) THEN
          CALL nh_prog_add_random( p_patch(jg), & ! input
               & p_nh_state(jg)%prog(nnow(jg)), & ! in and out
               & hs_nh_vn_ptb_scale, nproma, nlev ) ! input
          !
          CALL message(TRIM(routine),'Initial state used in the &
               & Held-Suarez test: random noised added to the normal wind')
      END IF
      !

    ENDDO !jg

  CASE ('APE_nh_tmp')  ! Aqua-Planet Experiment, no mountain

    DO jg = 1, n_dom
      CALL init_nh_state_prog_APE(p_patch(jg),p_nh_state(jg)%prog(nnow(jg)),    &
        &                            p_nh_state(jg)%diag,ext_data(jg),          &
        &                            p_nh_state(jg)%metrics,rh_at_1000hpa, qv_max)

      CALL init_nh_state_prog_APE(p_patch(jg),p_nh_state(jg)%prog(nnew(jg)),    &
        &                            p_nh_state(jg)%diag,ext_data(jg),          &
        &                            p_nh_state(jg)%metrics,rh_at_1000hpa, qv_max)
      


        CALL diagnose_pres_temp ( p_nh_state(jg)%metrics, p_nh_state(jg)%prog(nnow(jg)),    &
          &                       p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag,   &
          &                       p_patch(jg),                                &
          &                       opt_calc_temp=.TRUE.,                       &
          &                       opt_calc_pres=.TRUE.                      )
       
        CALL diagnose_pres_temp ( p_nh_state(jg)%metrics, p_nh_state(jg)%prog(nnow(jg)),    &
          &                       p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag,   &
          &                       p_patch(jg),                                &
          &                       opt_calc_temp=.TRUE.,                       &
          &                       opt_calc_pres=.TRUE.                       )       

    ENDDO !jg
  END SELECT

 END SUBROUTINE init_nh_testcase


 
END MODULE mo_nh_testcases
  

