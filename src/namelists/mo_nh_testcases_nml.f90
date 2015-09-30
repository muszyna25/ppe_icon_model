!>
!! Defines the artificial testcases for the nonhydrostatic atmospheric model.
!! 
!! 
!! @par Revision History
!! Initial release by Almut Gassmann (2008-03-18)
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
MODULE mo_nh_testcases_nml  
!-------------------------------------------------------------------------  
!  
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006                       
!  
!-------------------------------------------------------------------------  
!  
!  
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,                 ONLY: wp
  USE mo_namelist,             ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH 
  USE mo_io_units,             ONLY: nnml
  USE mo_nh_wk_exp,            ONLY: qv_max_wk, u_infty_wk,                          &
                                   & bubctr_lat, bubctr_lon, bubctr_z,               &
                                   & bub_hor_width, bub_ver_width, bub_amp
  USE mo_nh_lim_area_testcases,ONLY: nlayers_nconst,                                 &
                                   & p_base_nconst, theta0_base_nconst, h_nconst,    &
                                   & N_nconst, rh_nconst, rhgr_nconst,               &
                                   & itype_anaprof_uv, nlayers_linwind, h_linwind,   &
                                   & u_linwind, ugr_linwind, vel_const,              &
                                   & itype_topo_ana, schaer_h0,                      &
                                   & schaer_a, schaer_lambda, halfwidth_2dm,         &
                                   & mount_lonc_deg, mount_latc_deg, m_height,       &
                                   & m_width_x, m_width_y, itype_atmo_ana,           &
                                   & nlayers_poly,                                   &
                                   & p_base_poly, h_poly, t_poly,                    &
                                   & tgr_poly, rh_poly, rhgr_poly
  USE mo_nh_init_utils,        ONLY: n_flat_level, layer_thickness
  USE mo_nh_mrw_exp,           ONLY: mount_lonctr_mrw_deg, mount_latctr_mrw_deg,     &
                                   &  u0_mrw,  mount_height_mrw, mount_half_width,   &
                                   &  temp_i_mwbr_const, p_int_mwbr_const,           &
                                   &  bruntvais_u_mwbr_const
  USE mo_nh_dcmip_schaer,      ONLY: lshear_dcmip
  USE mo_nh_dcmip_gw,          ONLY: gw_clat, gw_u0, gw_delta_temp

  IMPLICIT NONE  

  PRIVATE 

  PUBLIC :: read_nh_testcase_namelist,layer_thickness,                       &
    &       n_flat_level, nh_test_name,                                      &
    &       ape_sst_case, ape_sst_val, w_perturb, th_perturb,                &
    &       mount_height, torus_domain_length, nh_brunt_vais, nh_u0, nh_t0,  &
    &       jw_up, rh_at_1000hpa,  qv_max,                                   &
    &       rotate_axis_deg, lhs_nh_vn_ptb, hs_nh_vn_ptb_scale,              & 
    &       linit_tracer_fv, lhs_fric_heat, lcoupled_rho, u_cbl, v_cbl,      &
    &       th_cbl, psfc_cbl, sol_const, zenithang, bubctr_x, bubctr_y

  CHARACTER(len=MAX_CHAR_LENGTH) :: nh_test_name
  CHARACTER(len=MAX_CHAR_LENGTH) :: ape_sst_case      !SST for APE experiments

  LOGICAL  :: lhs_nh_vn_ptb          ! if true, random noise is added to vn in HS_nh test case
  LOGICAL  :: lhs_fric_heat          ! if true, frictional heating is switched on in HS_nh
  LOGICAL  :: linit_tracer_fv  !< finite volume initialization for tracer fields
                               !< if .TRUE.
  LOGICAL  :: lcoupled_rho     !< re-integrate mass equation in PA test cases (TRUE/FALSE)

  REAL(wp) :: mount_height           ! (m)
  REAL(wp) :: nh_brunt_vais          ! (1/s)
  REAL(wp) :: nh_u0                  ! (m/s)
  REAL(wp) :: nh_t0                  ! (K)
  REAL(wp) :: jw_up                  ! amplitude of the u-perturbation (m/s), jabw  
  REAL(wp) :: rotate_axis_deg        ! (deg) rotation angle
  REAL(wp) :: torus_domain_length    ! (m) length of domain the slice (torus) grid
  REAL(wp) :: hs_nh_vn_ptb_scale     ! amplitude of the random noise
  REAL(wp) :: rh_at_1000hpa          ! relative humidity at 1000 hPa [%]
  REAL(wp) :: qv_max                 ! limit of maximum specific humidity in the tropics [kg/kg]
  REAL(wp) :: ape_sst_val            ! (K) value to be used for SST computation for aqua planet
  REAL(wp) :: w_perturb, th_perturb !Random perturbation scale for torus based experiments
  REAL(wp) :: sol_const              ! [W/m2] solar constant
  REAL(wp) :: zenithang              ! [degrees] zenith angle 

  !Linear profiles of variables for LES testcases
  REAL(wp) :: u_cbl(2)   !u_cbl(1) = constant, u_cbl(2) = gradient
  REAL(wp) :: v_cbl(2)   !v_cbl(1) = constant, v_cbl(2) = gradient
  REAL(wp) :: th_cbl(2)  !th_cbl(1) = constant,th_cbl(2) = gradient
  REAL(wp) :: psfc_cbl
  REAL(wp) :: bubctr_x  !X-Center of the warm bubble on torus
  REAL(wp) :: bubctr_y  !Y-Center of the warm bubble on torus

  NAMELIST/nh_testcase_nml/ nh_test_name, mount_height, torus_domain_length, &
                            nh_brunt_vais, nh_u0, nh_t0, layer_thickness,    &
                            n_flat_level, jw_up, u0_mrw, mount_height_mrw,   &
                            mount_half_width, mount_lonctr_mrw_deg,          &
                            mount_latctr_mrw_deg, p_int_mwbr_const,          &
                            temp_i_mwbr_const,  bruntvais_u_mwbr_const,      &
                            rotate_axis_deg,                                 &
                            lhs_nh_vn_ptb, hs_nh_vn_ptb_scale,               &
                            rh_at_1000hpa, qv_max,                           &
                            ape_sst_case, ape_sst_val,                       &
                            linit_tracer_fv, lhs_fric_heat,                  &
                            qv_max_wk, u_infty_wk,                           &
                            bubctr_lat, bubctr_lon, bubctr_z,                &
                            bub_hor_width, bub_ver_width, bub_amp,           &
                            nlayers_nconst,                                  &
                            p_base_nconst, theta0_base_nconst, h_nconst,     &
                            N_nconst, rh_nconst, rhgr_nconst,                &
                            itype_anaprof_uv, nlayers_linwind, h_linwind,    &
                            u_linwind, ugr_linwind, vel_const,               &
                            itype_topo_ana,                                  &
                            schaer_h0, schaer_a, schaer_lambda,              &
                            halfwidth_2dm, mount_lonc_deg, mount_latc_deg,   &
                            m_height, m_width_x, m_width_y, itype_atmo_ana,  &
                            nlayers_poly, p_base_poly, h_poly, t_poly,       &
                            tgr_poly, rh_poly, rhgr_poly, lshear_dcmip,      &
                            lcoupled_rho, gw_clat, gw_u0, gw_delta_temp,     & 
                            u_cbl, v_cbl, th_cbl, w_perturb, th_perturb,    &
                            psfc_cbl, sol_const, zenithang, bubctr_x, bubctr_y
                      

  CONTAINS
!-------------------------------------------------------------------------
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
    ape_sst_val            = 29.0_wp ! 29 degC
    sol_const              = 1361.371_wp ! [W/m2] default value for amip
    zenithang              = 38._wp ! value used for Popke et al. exps with no diurn cycle
    ! assuming that default is on triangles the next switch is set
    ! crosscheck follows in the respective module
    linit_tracer_fv        = .TRUE. ! finite volume initialization for tracer
    ! for Weisman-Klemp test case
    qv_max_wk              = 0.014_wp
    u_infty_wk             = 20._wp
    bub_amp                = 2.0_wp
    bubctr_lon             = 90.0_wp
    bubctr_lat             = 0.0_wp
    bubctr_z               = 1400._wp
    bub_hor_width          = 10000._wp
    bub_ver_width          = 1400._wp
    ! for the limited area test cases
    itype_atmo_ana   = 1
    ! for the piecewise const Brunt-Vaisala freq layers 
    p_base_nconst  = 100000.0_wp
    !theta0_base_nconst  = 300.0_wp
    theta0_base_nconst  = 288.0_wp
    nlayers_nconst  = 1
    h_nconst(:)  = 1.0_wp
      h_nconst(1)  = 0.0_wp
      h_nconst(2)  = 1500.0_wp
      h_nconst(3)  = 12000.0_wp
    N_nconst(:)  = 0.01_wp
      N_nconst(1)  = 0.01_wp
      N_nconst(2)  = 0.001_wp
      N_nconst(3)  = 0.02_wp
    rh_nconst(:)  = 0.5_wp     
   ! rel hum gradients,  positive for decreasing humidity with height
    rhgr_nconst(:)  = 0.0_wp
      rhgr_nconst(1)  = 0.0_wp ! 
      rhgr_nconst(2)  = 0.0_wp ! 
      rhgr_nconst(3)  = 0.0_wp
    ! for the piecewise polytropic layers
    nlayers_poly = 2
    p_base_poly = 100000._wp
    h_poly(:) = 1.0_wp
     h_poly(1) = 0._wp
     h_poly(2) = 12000._wp
    t_poly(:) = 280._wp
     t_poly(1) = 288._wp
     t_poly(2) = 213._wp
    tgr_poly(:) = 0._wp
     tgr_poly(1) = 0.009_wp
     tgr_poly(2) = 0.006_wp
    rh_poly(:) = 0.1_wp
     rh_poly(1)= 0.8_wp
     rh_poly(2)= 0.2_wp
    rhgr_poly(:)= 1.e-6_wp
     rhgr_poly(1)=5.e-5_wp
     rhgr_poly(2)=0._wp   !DR rhgr_poly(1)=0._wp
    ! for the wind profiles
    vel_const = 20.0_wp
    itype_anaprof_uv  = 1
    nlayers_linwind   = 2
    h_linwind(:)      = 1._wp
      h_linwind(1)    = 0._wp
      h_linwind(2)    = 2500._wp
    u_linwind(:)      = 20._wp
      u_linwind(1)    = 5._wp
      u_linwind(2)    = 10._wp
    ugr_linwind(:)    = 0._wp
      ugr_linwind(1)  = 0._wp
      ugr_linwind(2)  = 0._wp
     itype_topo_ana   = 1
    !For the Schaer mountain
    schaer_h0 = 250.0_wp
    schaer_a  = 5000.0_wp
    schaer_lambda = 4000._wp
    halfwidth_2dm = 100000.0_wp
    !For a mountain
    mount_lonc_deg = 90.0_wp
    mount_latc_deg = 0.0_wp
    m_height       = 1000.0_wp
    m_width_x      = 5000.0_wp
    m_width_y      = 5000.0_wp
    lshear_dcmip   = .FALSE.
    ! for PA test cases:
    lcoupled_rho   = .FALSE.

    ! for dcmip_gw_3X test cases
    gw_u0      = 0.0_wp      ! maximum amplitude of the zonal wind [m s^-1]
    gw_clat    = 90._wp      ! center of temperature/density perturbation  [deg]
    gw_delta_temp = 0.01_wp  ! Max amplitude of perturbation [K]


    !For CBL testcases, Anurag Dipankar (MPIM, 2013-04)
    u_cbl(1:2) = 0._wp 
    v_cbl(1:2) = 0._wp 
    th_cbl(1)  = 290._wp
    th_cbl(2)  = 0.006_wp
    psfc_cbl   = 102000._wp
    w_perturb  = 0.05_wp    
    th_perturb = 0.2_wp    

    !For warm bubble experiment on torus
    !Note that (0,0) is the center of the torus
    bubctr_x = 0._wp
    bubctr_y = 0._wp

    CALL open_nml(TRIM(filename))
    CALL position_nml ('nh_testcase_nml', status=i_status)
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nh_testcase_nml)
    END SELECT
    CALL close_nml

    n_flat_level=MAX(2,n_flat_level)

  END SUBROUTINE read_nh_testcase_namelist


!-------------------------------------------------------------------------
END MODULE mo_nh_testcases_nml  
