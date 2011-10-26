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
  USE mo_exception,          ONLY: message, finish
  USE mo_namelist,           ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, inwp
  USE mo_model_domain_import,ONLY: lplane, n_dom
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data,           ONLY: ext_data
  USE mo_math_constants,     ONLY: pi, pi_2
  USE mo_math_utilities,     ONLY: gc2cc, t_cartesian_coordinates, &
                                   t_geographical_coordinates, &
                                   arc_length
  USE mo_parallel_config,  ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: ltransport, ntracer, iforcing, iqv
  USE mo_extpar_config,         ONLY: itopo
    
  USE mo_dynamics_config,    ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_physical_constants, ONLY: grav, cpd, rd, cvd_o_rd, &
   &                               p0ref, re, omega, tmelt, vtmpc1, rv
! USE mo_convect_tables,     ONLY: B1   => c1es,  &
!  &                               B2_w => c3les ,&
!  &                               B2_i => c3ies ,&
!  &                               B4_w => c4les ,&
!  &                               B4_i => c4ies
  USE mo_nonhydro_state,       ONLY: t_nh_state, duplicate_prog_state
  
  
  USE mo_interpolation,        ONLY: t_int_state, cells2edges_scalar, edges2cells_scalar
  USE mo_mpi,                  ONLY: my_process_is_stdio
  USE mo_vertical_coord_table, ONLY: vct_b
  USE mo_loopindices,          ONLY: get_indices_e, get_indices_c
  USE mo_advection_config,     ONLY: advection_config
  USE mo_ncar_testcases,       ONLY: tracer_q1_q2, tracer_q3
  USE mo_nh_pa_test,           ONLY: init_nh_state_prog_patest
  USE mo_nh_df_test,           ONLY: init_nh_state_prog_dftest
  USE mo_nh_hs_test,           ONLY: init_nh_state_prog_held_suarez
  USE mo_nh_ape_exp,           ONLY: init_nh_state_prog_APE
  USE mo_nh_jabw_exp,          ONLY: init_nh_topo_jabw, init_nh_state_prog_jabw, & 
                                   & init_passive_tracers_nh_jabw, init_nh_inwp_tracers
  USE mo_nh_mrw_exp,           ONLY: init_nh_topo_mrw, init_nh_state_prog_mrw,    &
                                   & init_nh_prog_mwbr_const,                     &
                                   &  mount_lonctr_mrw_deg, mount_latctr_mrw_deg, &
                                   &  u0_mrw,  mount_height_mrw, mount_half_width,&
                                   &  temp_i_mwbr_const,                          &
                                   &  p_int_mwbr_const ,                          &
                                   &  bruntvais_u_mwbr_const
  USE mo_nh_wk_exp,            ONLY: init_nh_topo_wk, init_nh_env_wk,             &
                                   & qv_max_wk, u_infty_wk

  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_sync,                 ONLY: SYNC_E, sync_patch_array
  USE mo_satad,                ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                                sat_pres_ice!,  &  !! saturation vapor pressure w.r.t. ice
!   &                                spec_humi!,     &  !! Specific humidity
  USE mo_nh_prog_util,         ONLY: nh_prog_add_random
  USE mo_nh_init_utils,        ONLY: n_flat_level, layer_thickness

  
  IMPLICIT NONE  
  
  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  
  CHARACTER(len=MAX_CHAR_LENGTH) :: nh_test_name
  CHARACTER(len=MAX_CHAR_LENGTH) :: ape_sst_case      !SST for APE experiments

  REAL(wp) :: mount_height           ! (m)
  REAL(wp) :: nh_brunt_vais          ! (1/s)
  REAL(wp) :: nh_u0                  ! (m/s)
  REAL(wp) :: nh_t0                  ! (K)
  REAL(wp) :: jw_up                  ! amplitude of the u-perturbation (m/s), jabw  
  REAL(wp) :: rotate_axis_deg        ! (deg) rotation angle
  REAL(wp) :: torus_domain_length    ! (m) length of domain the slice (torus) grid
  LOGICAL  :: lhs_nh_vn_ptb          ! if true, random noise is added to vn in HS_nh test case
  LOGICAL  :: lhs_fric_heat          ! if true, frictional heating is switched on in HS_nh
  REAL(wp) :: hs_nh_vn_ptb_scale     ! amplitude of the random noise
  REAL(wp) :: rh_at_1000hpa          ! relative humidity at 1000 hPa [%]
  REAL(wp) :: qv_max                 ! limit of maximum specific humidity in the tropics [kg/kg]

  LOGICAL  :: linit_tracer_fv  !< finite volume initialization for tracer fields
                               !< if .TRUE.



  NAMELIST/nh_testcase_nml/ nh_test_name, mount_height, torus_domain_length, &
                            nh_brunt_vais, nh_u0, nh_t0, layer_thickness,    &
                            n_flat_level, jw_up, u0_mrw, mount_height_mrw,   &
                            mount_half_width, mount_lonctr_mrw_deg,          &
                            mount_latctr_mrw_deg, p_int_mwbr_const,          &
                            temp_i_mwbr_const,  bruntvais_u_mwbr_const,      &
                            rotate_axis_deg,                                 &
                            lhs_nh_vn_ptb, hs_nh_vn_ptb_scale,               &
                            rh_at_1000hpa, qv_max, ape_sst_case,             &
                            linit_tracer_fv, lhs_fric_heat,                  &
                            qv_max_wk, u_infty_wk

  PUBLIC :: read_nh_testcase_namelist, layer_thickness, init_nh_testtopo,    &
    &       init_nh_testcase, n_flat_level, nh_test_name, ape_sst_case,      &
    &       mount_height, torus_domain_length, nh_brunt_vais, nh_u0, nh_t0,  &
    &       jw_up, rh_at_1000hpa,  qv_max,                                   &
    &       rotate_axis_deg, lhs_nh_vn_ptb, hs_nh_vn_ptb_scale,              &
    &       linit_tracer_fv, lhs_fric_heat

  PRIVATE

! !DEFINED PARAMETERS for jablonowski williamson: 
!  The rest of the needed parameters are define in mo_nh_jabw_exp
  REAL(wp), PARAMETER :: ps0   = 1.e5_wp  ! surface pressure (Pa)
  
! !DEFINED PARAMETERS for mountain induced Rossby wave train:
  REAL(wp), PARAMETER :: pres_sp  = 93000.0_wp  !pressure surface at the south pole
  REAL(wp), PARAMETER :: temp_mrw = 288._wp     !temperature of isothermal atmosphere

! !DEFINED PARAMETERS for APE
   REAL(wp), PARAMETER :: zp_ape      = 101325._wp            !< surface pressure
   REAL(wp), PARAMETER :: ztmc_ape    = 25.006_wp             !< total moisture content 

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
  SUBROUTINE read_nh_testcase_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: i_status

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
    ! assuming that default is on triangles the next switch is set
    ! crosscheck follows in the respective module
    linit_tracer_fv        = .TRUE. ! finite volume initialization for tracer
    ! for Weisman-Klemp test case
    qv_max_wk              = 0.014_wp
    u_infty_wk             = 20._wp
    

    CALL open_nml(TRIM(filename))
    CALL position_nml ('nh_testcase_nml', status=i_status)
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nh_testcase_nml)
    END SELECT
    CALL close_nml

    n_flat_level=MAX(2,n_flat_level)

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nh_testcase_nml)

  END SUBROUTINE read_nh_testcase_namelist


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
  INTEGER        :: nblks_c, npromz_c, nblks_v, npromz_v
  REAL(wp)       :: z_lon, z_lat, z_dist
  TYPE(t_geographical_coordinates) :: z_x2_geo
  TYPE(t_cartesian_coordinates)    :: z_x1_cart, z_x2_cart
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = '(mo_nh_testcases) init_nh_testtopo:'
 
  LOGICAL       :: l_modified

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

  CASE ('jabw', 'jabw_s')

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_int_c
     npromz_c  = p_patch(jg)%npromz_c
     nblks_v   = p_patch(jg)%nblks_int_v
     npromz_v  = p_patch(jg)%npromz_v

     CALL init_nh_topo_jabw ( p_patch(jg),ext_data(jg)%atm%topography_c,  &
                          & ext_data(jg)%atm%topography_v, nblks_c, npromz_c, &
                          & nblks_v, npromz_v)
    END DO

  CASE ('jabw_m')  

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_int_c
     npromz_c  = p_patch(jg)%npromz_c
     nblks_v   = p_patch(jg)%nblks_int_v
     npromz_v  = p_patch(jg)%npromz_v

     CALL init_nh_topo_jabw ( p_patch(jg),ext_data(jg)%atm%topography_c,  &
                          & ext_data(jg)%atm%topography_v, nblks_c, npromz_c, &
                          & nblks_v, npromz_v, &
                          & opt_m_height = mount_height, &
                          & opt_m_half_width = mount_half_width )
    END DO

  CASE ('mrw_nh', 'mrw2_nh' , 'mwbr_const')

   CALL message(TRIM(routine),'running mrw, setting topography')

   l_modified = .FALSE.
   IF (nh_test_name=='mrw2_nh'  ) THEN
     l_modified = .TRUE.
   ENDIF

   DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_int_c
     npromz_c  = p_patch(jg)%npromz_c
     nblks_v   = p_patch(jg)%nblks_int_v
     npromz_v  = p_patch(jg)%npromz_v

     CALL init_nh_topo_mrw ( p_patch(jg),ext_data(jg)%atm%topography_c,  &
                          & ext_data(jg)%atm%topography_v, nblks_c, npromz_c, &
                          & nblks_v, npromz_v, l_modified) 
   ENDDO

   CALL message(TRIM(routine),'topography is initialised ')

  CASE ('wk82')  

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_int_c
     npromz_c  = p_patch(jg)%npromz_c
     nblks_v   = p_patch(jg)%nblks_int_v
     npromz_v  = p_patch(jg)%npromz_v

     CALL init_nh_topo_wk ( p_patch(jg),ext_data(jg)%atm%topography_c,  &
                          & ext_data(jg)%atm%topography_v, nblks_c, npromz_c, &
                          & nblks_v, npromz_v )
    END DO

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
    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the Held-Suarez test')
  CASE ('APE_nh')

   ! The topography has been initialized to 0 at the begining of this SUB
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

  INTEGER        :: jg, je, jc, jb, jk, jt,   &
                    nlen, nblks_e, npromz_e,  nblks_c, npromz_c
  INTEGER        :: nlev, nlevp1        !< number of full and half levels

  TYPE(t_nh_state), POINTER       :: p_nhdom
                            
  REAL(wp)              :: p_sfc_jabw  ! surface pressure for the jabw test case, 
                                       ! standard values is 100000 Pa   
  REAL(wp)              :: global_moist
  REAL(wp) :: z_help

  LOGICAL  :: l_hydro_adjust, l_moist
  
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                                   '(mo_nh_testcases) init_nh_testcase:' 
! Tracer related variables
  CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list
!-----------------------------------------------------------------------

  SELECT CASE (nh_test_name)

  CASE ('jabw', 'jabw_s', 'jabw_m', 'HS_jw') ! Jablonowski test

  CALL message(TRIM(routine),'Jablonowski test')
  IF ( iforcing == inwp ) THEN
    CALL message(TRIM(routine),' iforcing == inwp')
  ELSE
    CALL message(TRIM(routine),'Attention: iforcing /= inwp')
  ENDIF

  IF (nh_test_name == "jabw_s" .OR. nh_test_name == "jabw_m") THEN
   jw_up = 0.0_wp
  END IF  
   p_sfc_jabw = ps0

  DO jg = 1, n_dom

    CALL   init_nh_state_prog_jabw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                   & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                   & p_int(jg),                                   &
                                   & p_sfc_jabw,jw_up )


    IF ( ltransport .AND. iforcing /= inwp ) THEN   ! passive tracers
       ! get ctracer_list
       ctracer_list = advection_config(jg)%ctracer_list

       CALL init_passive_tracers_nh_jabw (p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                         & rotate_axis_deg, ctracer_list, p_sfc_jabw) 

    END IF

    IF ( ltransport .AND. iforcing == inwp ) THEN 
     IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                   &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN   !


      CALL init_nh_inwp_tracers ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                & rh_at_1000hpa, qv_max, l_rediag=.TRUE. )




     ELSE

        p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,:) = 0.0_wp   

     END IF

    END IF

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

  ENDDO !jg

  CALL message(TRIM(routine),'End setup Jablonowski test')
  
  CASE ('mrw_nh', 'mrw2_nh')

   CALL message(TRIM(routine),'MRW test')
  
   l_hydro_adjust = .TRUE.
   l_moist = .FALSE.

   DO jg = 1, n_dom

     IF ( iforcing == inwp ) THEN
       CALL message(TRIM(routine),' iforcing == inwp')     
       IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                      &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN 
         l_moist = .TRUE.
       END IF
     ENDIF


     IF (.NOT. l_moist) THEN

       CALL   init_nh_state_prog_mrw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & ext_data(jg)%atm%topography_c,              &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg), l_hydro_adjust, iforcing, l_moist)

    ELSE

       CALL   init_nh_state_prog_mrw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & ext_data(jg)%atm%topography_c,              &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg), l_hydro_adjust, iforcing, l_moist , &
                                     & opt_rh_at_1000hpa= rh_at_1000hpa,           &
                                     & opt_qv_max=qv_max                           )
    END IF
   
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))  
 
  ENDDO !jg

   CALL message(TRIM(routine),'End setup MRW test')

  CASE ('mwbr_const')

   CALL message(TRIM(routine),'mwbr_const test case')
  
   l_hydro_adjust = .TRUE.
   l_moist = .FALSE.

   DO jg = 1, n_dom


     IF ( iforcing == inwp ) THEN
       CALL message(TRIM(routine),' iforcing == inwp')     
       IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                      &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN 
         l_moist = .TRUE.
       END IF
     ENDIF 




     IF (.NOT. l_moist) THEN

       CALL   init_nh_prog_mwbr_const ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                      & p_nh_state(jg)%diag,                        &
                                      & ext_data(jg)%atm%topography_c,              &
                                      & p_nh_state(jg)%metrics,                     &
                                      & p_int(jg), l_hydro_adjust, iforcing, l_moist)

     ELSE
  
       CALL   init_nh_prog_mwbr_const ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                      & p_nh_state(jg)%diag,                        &
                                      & ext_data(jg)%atm%topography_c,              &
                                      & p_nh_state(jg)%metrics,                     &
                                      & p_int(jg), l_hydro_adjust, iforcing, l_moist,&
                                      & opt_rh_at_1000hpa= rh_at_1000hpa,           &
                                      & opt_qv_max=qv_max                           )

    END IF

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg))) 
 
 ENDDO !jg

   CALL message(TRIM(routine),'End setup mwbr_const test')

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



      IF (lhs_nh_vn_ptb) THEN
          CALL nh_prog_add_random( p_patch(jg), & ! input
               & p_nh_state(jg)%prog(nnow(jg)), & ! in and out
               & hs_nh_vn_ptb_scale, nproma, nlev ) ! input
          !
          CALL message(TRIM(routine),'Initial state used in the &
               & Held-Suarez test: random noised added to the normal wind')
      END IF
      !

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

    ENDDO !jg

  CASE ('APE_nh')  ! Aqua-Planet Experiment, no mountain

    p_sfc_jabw   = zp_ape          ! Pa
    global_moist = ztmc_ape        ! kg/m**2 total moisture content
    jw_up = 1._wp



  DO jg = 1, n_dom

    CALL   init_nh_state_prog_jabw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                   & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                   & p_int(jg),                                   &
                                   & p_sfc_jabw,jw_up )



    IF ( ltransport .AND. iforcing == inwp ) THEN   !


      CALL init_nh_inwp_tracers ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                & rh_at_1000hpa, qv_max, l_rediag=.TRUE.,  &
                                & opt_global_moist=global_moist)
 

    END IF

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

  ENDDO !jg

   CALL message(TRIM(routine),'End setup APE_nh test')
  
  CASE ('wk82')

   CALL message(TRIM(routine),'wk82 test')
  
   l_hydro_adjust = .TRUE.

   DO jg = 1, n_dom


  
    CALL   init_nh_env_wk ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg))

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

  ENDDO !jg

   CALL message(TRIM(routine),'End setup wk82 test')

  END SELECT

 END SUBROUTINE init_nh_testcase


 
END MODULE mo_nh_testcases
  

