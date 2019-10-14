!>
!! Initialization routines for lahade-type testcases: 
!! Given an atmosphere at rest in the inertial/absolute frame 
!! and a spherical bottom and rigid lid at model top, rotating 
!! with a constant angular velocity (slip condition at bottom and lid 
!! -> no flux of tangential momentum between boundary and air possible), 
!! analytical solutions of the linearized governing equations are considered  
!! and compared to the numerical solution of ICON relative to 
!! the rotating boundaries: thereby the implementation of the Coriolis acceleration 
!! and the metric terms of advection can be evaluated. 
!! Note: the lahade-testcases are intended for the deep-atmosphere dynamics 
!! (dynamics_config: ldeepatmo = .true.)
!!
!! @Literature
!! - M. LAeuter, D. HAndorf and K. DEthloff (2005) "Unsteady analytical solutions 
!!   of the spherical shallow water equations", J. Comp. Phys., 210, 535-553
!!   (It seems to be custom to name testcases somehow after the original inventors, 
!!   that is why we call it 'LAHADE'.)
!! - A. Staniforth and A. A. White (2008) "Unsteady exact solutions of the flow 
!!   equations for three-dimensional spherical atmospheres", Q. J. R. Meteorol. Soc., 
!!   134, 1615-1626.
!! - M. Baldauf, D. Reinert and G. Zaengl (2014) "An analytical solution for linear 
!!   gravity and sound waves on the sphere as a test for compressible, non-hydrostatic 
!!   numerical models", Q. J. R. Meteorol. Soc., 140, 1974-1985.
!!
!! @par Revision History
!! Initial revision by Sebastian Borchert, DWD (2017-04-21)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nh_lahade

  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: finish, message, message_text
  USE mo_impl_constants,           ONLY: min_rlcell, min_rledge, inoforcing, &
    &                                    MAX_CHAR_LENGTH
  USE mo_loopindices,              ONLY: get_indices_c, get_indices_e
  USE mo_math_constants,           ONLY: pi, deg2rad, dbl_eps
  USE mo_physical_constants,       ONLY: rd, cpd, cvd, rd_o_cpd, p0ref, & 
    &                                    earth_angular_velocity,        &
    &                                    earth_radius, grav
  USE mo_model_domain,             ONLY: t_patch, t_tangent_vectors
  USE mo_nonhydro_types,           ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,           ONLY: t_int_state
  USE mo_intp,                     ONLY: cells2edges_scalar
  USE mo_parallel_config,          ONLY: nproma
  USE mo_sync,                     ONLY: sync_patch_array, SYNC_E
  USE mo_run_config,               ONLY: msg_level
  USE mo_upatmo_config,            ONLY: t_upatmo_dyn_config, t_upatmo_config, &
    &                                    imsg_thr, istatus
  USE mo_name_list_output_types,   ONLY: t_output_name_list
  USE mo_name_list_output_config,  ONLY: is_variable_in_output
  USE mo_name_list_output,         ONLY: istime4name_list_output
  USE mo_util_string,              ONLY: int2string

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_lahade'

  PUBLIC :: lahade, ilahade
  PUBLIC :: init_nh_lahade
  PUBLIC :: nh_lahade_interface
  PUBLIC :: check_nh_lahade

  ! (The following data types are defined here
  ! and not in 'src/namlists/mo_nh_testcases_nml', 
  ! in order to avoid circular dependencies)
  TYPE t_lahade
    INTEGER  :: icase                  !nml-input
    REAL(wp) :: omega                  !nml-input
    REAL(wp) :: bkg_temp               !nml-input
    REAL(wp) :: bkg_pres               !nml-input
    REAL(wp) :: ptb_ctr_lat            !nml-input
    REAL(wp) :: ptb_ctr_lon            !nml-input
    REAL(wp) :: ptb_ctr_hgt            !nml-input
    REAL(wp) :: ptb_rad_min            !nml-input
    REAL(wp) :: ptb_rad_max            !nml-input
    REAL(wp) :: ptb_amp_temp           !nml-input
    REAL(wp) :: ptb_n_rad              !nml-input
    CHARACTER(len=MAX_CHAR_LENGTH) :: output_ptb_var  !nml-input
    LOGICAL  :: lactive                !derived
    LOGICAL  :: lzerograv              !derived
    INTEGER  :: ivarout                !derived
    LOGICAL  :: lupdate                !derived
  END TYPE t_lahade

  ! Identifiers for lahade%icase
  TYPE t_icase
    INTEGER :: ssw       ! Spherical sound wave (currently the only sub-case)
  END TYPE t_icase
  TYPE t_ivarout
    INTEGER :: default   ! No variable subject to output
    INTEGER :: temp      ! Analytical field of perturbation temperature
    INTEGER :: rho       ! --,,-- density
    INTEGER :: pres      ! --,,-- pressure
  END TYPE t_ivarout
  TYPE t_ilahade
    TYPE(t_icase)   :: case
    TYPE(t_ivarout) :: varout
  END TYPE t_ilahade

  TYPE(t_lahade)             :: lahade
  TYPE(t_ilahade), PARAMETER :: ilahade = t_ilahade(    &  !ilahade%
    &                                     t_icase(      &  !ilahade%case%
    &                                                1  &  !ilahade%case%ssw
    &                                            ),     &
    &                                     t_ivarout(    &  !ilahade%varout%
    &                                                0, &  !ilahade%varout%default
    &                                                1, &  !ilahade%varout%temp
    &                                                2, &  !ilahade%varout%rho
    &                                                3  &  !ilahade%varout%pres
    &                                              )    &
    &                                              )      

  ! For convenience:
  TYPE t_id
    INTEGER :: gc  ! Geographical coordinates
    INTEGER :: cc  ! Cartesian coordinates
  END TYPE t_id
  ! Indices for Cartesian components
  TYPE t_icc
    INTEGER :: x
    INTEGER :: y
    INTEGER :: z
  END TYPE t_icc
  ! Indices for geographical coordinates
  TYPE t_igc
    INTEGER :: lon  ! Zonal position
    INTEGER :: lat  ! Meridional position
    INTEGER :: rad  ! Radial position
  END TYPE t_igc
  ! Identifier for direction of coordinate transformation
  TYPE t_itrafo
    INTEGER :: gc2cc  ! Geographical to Cartesian coordinates
    INTEGER :: cc2gc  ! Cartesian to geographical coordinates
  END TYPE t_itrafo
  TYPE t_icoord
    type(t_id)     :: id
    TYPE(t_icc)    :: cc
    TYPE(t_igc)    :: gc
    type(t_itrafo) :: trafo
  END TYPE t_icoord
  TYPE(t_icoord), PARAMETER :: icoord = t_icoord(    &  !icoord%
    &                                   t_id(        &  !icoord%id%
    &                                             1, &  !icoord%id%gc
    &                                             2  &  !icoord%id%cc
    &                                       ),       &                               
    &                                   t_icc(       &  !icoord%cc%
    &                                             1, &  !icoord%cc%x
    &                                             2, &  !icoord%cc%y
    &                                             3  &  !icoord%cc%z
    &                                        ),      &
    &                                   t_igc(       &  !icoord%gc%
    &                                             1, &  !icoord%gc%lon
    &                                             2, &  !icoord%gc%lat
    &                                             3  &  !icoord%gc%rad
    &                                        ),      &
    &                                   t_itrafo(    &  !icoord%trafo%
    &                                             1, &  !icoord%trafo%gc2cc
    &                                             2  &  !icoord%trafo%cc2gc
    &                                           )    &
    &                                           )

  TYPE t_isol
    INTEGER :: num  ! 'extra_3d(..., ...%num)': numerical solution
    INTEGER :: ana  ! 'extra_3d(..., ...%ana)': analytical solution
  END TYPE t_isol
  TYPE(t_isol), PARAMETER :: isol = t_isol(    &  !isol%
    &                                       1, &  !isol%num
    &                                       2  &  !isol%ana
    &                                     )

CONTAINS

  !>
  !! Setup idealized initial conditions for lahade-testcases.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2017-04-21)
  !!
  SUBROUTINE init_nh_lahade ( p_patch,               &  !in
    &                         p_nh_prog,             &  !inout
    &                         p_nh_diag,             &  !inout
    &                         p_int,                 &  !in
    &                         p_metrics,             &  !inout
    &                         grid_sphere_radius,    &  !in
    &                         grid_angular_velocity, &  !in
    &                         grid_rescale_factor,   &  !in
    &                         top_height,            &  !in
    &                         vwind_offctr,          &  !in
    &                         ndyn_substeps,         &  !in
    &                         upatmo_config          )  !in

    ! In/out variables

    TYPE(t_patch),      TARGET, INTENT(INOUT) :: p_patch                !< Patch on which computation is performed
    TYPE(t_nh_prog),            INTENT(INOUT) :: p_nh_prog              !< Prognostic state vector
    TYPE(t_nh_diag),            INTENT(INOUT) :: p_nh_diag              !< Diagnostic state vector
    TYPE(t_int_state),          INTENT(IN)    :: p_int                  !< Interpolation state vector
    TYPE(t_nh_metrics),         INTENT(INOUT) :: p_metrics              !< NH metrics state
    REAL(wp),                   INTENT(IN)    :: grid_sphere_radius     !< (Rescaled) Earth radius
    REAL(wp),                   INTENT(IN)    :: grid_angular_velocity  !< (Rescaled) Earth angular velocity
    REAL(wp),                   INTENT(IN)    :: grid_rescale_factor    !< Rescale factor
    REAL(wp),                   INTENT(IN)    :: top_height             !< Height of model top
    REAL(wp),                   INTENT(IN)    :: vwind_offctr           !< Vertical wind off-centering
    INTEGER,                    INTENT(IN)    :: ndyn_substeps          !< Nominal value of dt_fastphy / dt_dyn
    TYPE(t_upatmo_config),      INTENT(IN)    :: upatmo_config          !< Upper-atmosphere configuration   
 
    ! Local variables

    TYPE(t_tangent_vectors) :: z_en
    REAL(wp), DIMENSION(3)  :: z_ptb_ctr_gc, z_pos_gc
    REAL(wp), DIMENSION(3)  :: z_ptb_ctr_cc, z_pos_cc, z_distvec_cc
    REAL(wp) :: z_r, z_dist, z_dist_rel 
    REAL(wp) :: z_dx_min, z_dz_min, z_vx_max, z_vz_max
    REAL(wp) :: z_en_norm, z_distvec_dot_en, z_distvec_dot_er
    REAL(wp) :: z_sin_lat, z_sin_lon, z_cos_lat, z_cos_lon
    REAL(wp) :: z_bkg_pres, z_bkg_temp, z_bkg_rho, z_bkg_vn, z_bkg_vn_max
    REAL(wp) :: z_ptb_pres, z_ptb_temp, z_ptb_rho, z_ptb_phi, z_ptb_vn, z_ptb_w
    REAL(wp) :: z_ptb_amp_temp, z_ptb_amp_pres, z_ptb_amp_rho, z_ptb_amp_vr
    REAL(wp) :: z_cs, z_dtime_dyn_given, z_dtime_dyn_ubound

    REAL(wp) :: z_me(nproma,p_patch%nlev,p_patch%nblks_e)  ! Height of edge centers

    INTEGER :: jg, jc, je, jk, jb          ! Loop indices for domain, cell, edge, level, block
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER :: nlev                        ! Number of vertical full levels
    INTEGER :: nlevp1                      ! Number of vertical half levels

    LOGICAL :: lmessage
 
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':init_nh_lahade'

    !----------------------------------------------------------------------------

    ! Domain index
    jg = p_patch%id

    ! Number of grid layers
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Print messages?
    lmessage = (msg_level >= imsg_thr%low)

    ! Different variants of lahade-testcases
    SELECT CASE(lahade%icase)
    CASE(ilahade%case%ssw)         
      ! Spherical sound wave (currently the only lahade-testcase): 
      ! gravity is switched off -> air is homogeneous gas trapped 
      ! between concentrical spherical bottom shell and spherical lid shell

      IF (lmessage) CALL message(TRIM(routine), &
        & 'Setup of lahade-testcase: "spherical sound wave" on domain '//TRIM(int2string(jg))//' started')

      !-----------------
      ! 0th Preparation
      !-----------------     

      ! Center of spherical sound wave in geographical coordinates.
      ! (Note: the data types:
      !  - 'src/shared/mo_math_types: t_geographical_coordinates' 
      !  - 'src/shared/mo_math_types: t_cartesian_coordinates' 
      ! are not suitable here, since they are 'normalized' with the radius. 
      ! In addition 'ptb_ctr_lon' and 'ptb_ctr_lat' have already been converted 
      ! to radian in subroutine 'check_nh_lahade' below.)
      z_ptb_ctr_gc(icoord%gc%lon) = lahade%ptb_ctr_lon                       ! Zonal
      z_ptb_ctr_gc(icoord%gc%lat) = lahade%ptb_ctr_lat                       ! Meridional
      z_ptb_ctr_gc(icoord%gc%rad) = grid_sphere_radius + lahade%ptb_ctr_hgt  ! Radial
      
      ! And Cartesian coordinates of center
      z_ptb_ctr_cc = coord_trafo( z_ptb_ctr_gc, icoord%trafo%gc2cc )

      ! In this testcase variant the background state 
      ! of the atmospher does not depend on position, 
      ! so it can be precomputed here
      z_bkg_pres = lahade%bkg_pres
      z_bkg_temp = lahade%bkg_temp
      z_bkg_rho  = z_bkg_pres / ( rd * z_bkg_temp )      

      ! Max. background velocity magnitude at equator and model top
      z_bkg_vn_max = grid_angular_velocity * ( grid_sphere_radius + top_height )

      ! Speed of sound 
      z_cs = SQRT( ( cpd / cvd ) * rd * z_bkg_temp )

      ! Check time step
      IF (jg == 1) THEN
        ! The given time step from namelist input
        z_dtime_dyn_given = upatmo_config%dt_dyn_nom
        ! The upper bound for the time step from the CFL-criterion:
        !
        ! upper-bound-for-time-step ~ min-grid-mesh-size / max-velocity-magnitude
        !
        ! (Both, the horizontal mesh size and the background wind increase 
        ! linearly with height. So checking the criterion at the surface 
        ! is representative for the entire column.)
        z_dx_min = p_patch%geometry_info%mean_characteristic_length
        z_vx_max = MAX( z_cs, ABS(grid_angular_velocity * grid_sphere_radius) )
        ! (The grid layer thickness should be the same for all layers. 
        ! In addition, there is no vertical background velocity.)
        z_dz_min = top_height / REAL(nlev, wp)
        z_vz_max = z_cs
        ! Upper time step bound
        z_dtime_dyn_ubound = MIN( z_dx_min / z_vx_max, &
          &                       z_dz_min / z_vz_max  )
        ! The timestep is a 'critical' quantity, which should rather be reset by hand in the runscript
        IF ( z_dtime_dyn_given > z_dtime_dyn_ubound ) THEN
          WRITE(message_text, '(a,E12.6)') 'The simulation would likely become unstable. ' // &
            & 'Please, choose for dtime a value smaller than ', & 
            & z_dtime_dyn_ubound * REAL( ndyn_substeps, wp ) / grid_rescale_factor
          CALL finish(TRIM(routine),TRIM(message_text))        
        ENDIF
      ENDIF  !IF (jg == 1)

      ! Amplitudes for sound wave
      z_ptb_amp_temp = lahade%ptb_amp_temp                                    ! Temperature (from namelist)
      z_ptb_amp_vr   = ( cvd / rd ) * ( z_ptb_amp_temp / z_bkg_temp ) * z_cs  ! Radial velocity
      z_ptb_amp_pres = ( cpd / cvd ) * ( z_ptb_amp_vr / z_cs ) * z_bkg_pres   ! Pressure
      z_ptb_amp_rho  = ( z_ptb_amp_vr / z_cs ) * z_bkg_rho                    ! Density

      !-------------------------------------------------
      ! 1st Initialize scalar quantities at mass points
      !-------------------------------------------------

      i_rlstart  = 1
      i_rlend    = min_rlcell
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_pos_cc,z_pos_gc,z_dist,z_dist_rel, & 
!$OMP            z_ptb_phi,z_ptb_pres,z_ptb_temp,z_ptb_rho)
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx 
            
            !----------------------------------------
            ! Variables of the background atmosphere
            !----------------------------------------
            
            ! See preparation above
            
            !-------------------------------------
            ! Add initial sound wave perturbation
            !-------------------------------------
            
            ! Geographical coordinates of current position. 
            ! (In principle, 'p_patch%cells%cartesian_center' should 
            ! contain this information, but to make sure that the values 
            ! are really consistent with 'z_ptb_ctr_cc', they are recomputed here)
            z_pos_gc(icoord%gc%lon) = p_patch%cells%center(jc,jb)%lon                ! Zonal
            z_pos_gc(icoord%gc%lat) = p_patch%cells%center(jc,jb)%lat                ! Meridional
            z_pos_gc(icoord%gc%rad) = grid_sphere_radius + p_metrics%z_mc(jc,jk,jb)  ! Radial

            ! Transform geographical to Cartesian coordinates
            z_pos_cc = coord_trafo( z_pos_gc, icoord%trafo%gc2cc )
            
            ! Distance of current position to center of spherical sound wave
            z_dist = SQRT( DOT_PRODUCT( ( z_pos_cc - z_ptb_ctr_cc ), ( z_pos_cc - z_ptb_ctr_cc ) ) )
            
            ! The initial sound wave has a non-vanishing amplitude only for 
            ! 'lahade%ptb_rad_min <= z_dist <= lahade%ptb_rad_max'. 
            ! Non-dimensional distance of current position relative to 'lahade%ptb_rad_min'
            z_dist_rel = ( z_dist - lahade%ptb_rad_min ) / ( lahade%ptb_rad_max - lahade%ptb_rad_min ) 
            
            ! Compute non-dimensional sound wave field, 
            ! which is the basis for all other scalar fields 
            IF ( z_dist_rel >= 0._wp .AND. z_dist_rel <= 1._wp ) THEN
              ! ('lahade%ptb_n_rad' is the number of wave crests in radial direction)
              z_ptb_phi = SIN( pi * z_dist_rel ) * SIN( 2._wp * pi * lahade%ptb_n_rad * z_dist_rel )  & 
                &         + ( lahade%ptb_rad_max - lahade%ptb_rad_min ) / z_dist *                    & 
                &         ( SIN( pi * ( 2._wp * lahade%ptb_n_rad - 1._wp ) * z_dist_rel ) /           & 
                &         ( 2._wp * pi * ( 2._wp * lahade%ptb_n_rad - 1._wp ) )                       & 
                &         - SIN( pi * ( 2._wp * lahade%ptb_n_rad + 1._wp ) * z_dist_rel ) /           & 
                &         ( 2._wp * pi * ( 2._wp * lahade%ptb_n_rad + 1._wp ) ) )     
            ELSE
              z_ptb_phi = 0._wp
            ENDIF
            
            ! Compute pressure, temperature and density perturbation fields of sound wave
            z_ptb_pres = z_ptb_amp_pres * z_ptb_phi
            z_ptb_temp = z_ptb_amp_temp * z_ptb_phi
            z_ptb_rho  = z_ptb_amp_rho * z_ptb_phi
            
            ! Initialize ICON fields (= background state + sound wave perturbation)
            p_nh_diag%pres(jc,jk,jb)    = z_bkg_pres + z_ptb_pres
            p_nh_diag%temp(jc,jk,jb)    = z_bkg_temp + z_ptb_temp
            p_nh_prog%rho(jc,jk,jb)     = z_bkg_rho + z_ptb_rho
            p_nh_prog%exner(jc,jk,jb)   = ( p_nh_diag%pres(jc,jk,jb) / p0ref )**rd_o_cpd
            p_nh_prog%theta_v(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb) / p_nh_prog%exner(jc,jk,jb)
            p_nh_diag%tempv(jc,jk,jb)   = p_nh_diag%temp(jc,jk,jb)
            
          ENDDO  !jc
        ENDDO  !jk
      ENDDO  !jb
!$OMP ENDDO
!$OMP END PARALLEL

      !-----------------------------------------------
      ! 2nd Initialize normal wind component at edges
      !-----------------------------------------------

      ! Interpolate height to edge centers 
      CALL cells2edges_scalar(p_metrics%z_mc, p_patch, p_int%c_lin_e, z_me)
      CALL sync_patch_array(SYNC_E, p_patch, z_me)

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

      i_rlstart  = 1
      i_rlend    = min_rledge
      i_startblk = p_patch%edges%start_block(i_rlstart)
      i_endblk   = p_patch%edges%end_block(i_rlend)
      
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_r,z_pos_cc,z_pos_gc,z_distvec_cc,z_dist,z_en,z_en_norm, & 
!$OMP            z_sin_lat,z_sin_lon,z_cos_lat,z_cos_lon,z_distvec_dot_en,z_dist_rel,z_bkg_vn,z_ptb_vn) 

      DO jb = i_startblk, i_endblk
        
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)
        
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            ! Geographical coordinates of current position 
            z_pos_gc(icoord%gc%lon) = p_patch%edges%center(je,jb)%lon      ! Zonal
            z_pos_gc(icoord%gc%lat) = p_patch%edges%center(je,jb)%lat      ! Meridional
            z_pos_gc(icoord%gc%rad) = grid_sphere_radius + z_me(je,jk,jb)  ! Radial

            ! Transform geographical to Cartesian coordinates.
            ! (The Cartesian coordinates of the edge center are 
            ! stored in 'p_patch%edges%cartesian_center' 
            ! and the Cartesian coordinates of the edge normal unit vector 
            ! are stored in 'p_patch%edges%primal_cart_normal', 
            ! but they are read from file, so to be really sure, 
            ! that the same Cartesian basis is used, 
            ! we recompute these quantities here)
            z_pos_cc = coord_trafo( z_pos_gc, icoord%trafo%gc2cc )

            ! Compute normalized distance vector from the center 
            ! of the spherical sound wave to the current position
            z_distvec_cc = z_pos_cc - z_ptb_ctr_cc 
            z_dist       = SQRT( DOT_PRODUCT( z_distvec_cc, z_distvec_cc ) )
            z_distvec_cc = z_distvec_cc / z_dist

            ! Edge primal normal vector in spherical coordinates
            z_en = p_patch%edges%primal_normal(je,jb)
            ! And just to make sure
            z_en_norm = SQRT( z_en%v1**2 + z_en%v2**2 ) 
            z_en%v1   = z_en%v1 / z_en_norm  ! Zonal component
            z_en%v2   = z_en%v2 / z_en_norm  ! Meridional component   

            ! Scalar product of distance vector and primal normal vector.
            ! (The spherical unit vectors measured in Cartesian coordinates are: 
            ! e_lon =          -sin(lon)*e_x +          cos(lon)*e_y, 
            ! e_lat = -sin(lat)*cos(lon)*e_x - sin(lat)*sin(lon)*e_y + cos(lat)*e_z, 
            ! e_r   =  cos(lat)*cos(lon)*e_x + cos(lat)*sin(lon)*e_y + sin(lat)*e_z, 
            ! The primal normal vector in spherical coordinates is: 
            ! e_n = v1*e_lon + v2*e_lat )
            z_sin_lat = SIN( p_patch%edges%center(je,jb)%lat )
            z_sin_lon = SIN( p_patch%edges%center(je,jb)%lon )
            z_cos_lat = COS( p_patch%edges%center(je,jb)%lat )
            z_cos_lon = COS( p_patch%edges%center(je,jb)%lon )
            z_distvec_dot_en = z_en%v1 *                                           &
              &                (                                                   &
              &                -z_distvec_cc(icoord%cc%x) * z_sin_lon              & 
              &                +z_distvec_cc(icoord%cc%y) * z_cos_lon              & 
              &                )                                                   &
              &                + z_en%v2 *                                         &
              &                (                                                   & 
              &                -z_distvec_cc(icoord%cc%x) * z_sin_lat * z_cos_lon  & 
              &                -z_distvec_cc(icoord%cc%y) * z_sin_lat * z_sin_lon  & 
              &                +z_distvec_cc(icoord%cc%z) * z_cos_lat              &
              &                )

            ! Compute background wind: u0 * e_lon = -Omega * r * cos(lat) * e_lon
            ! (= -Fuehrungsgeschwindigkeit)
            z_r      = z_pos_gc(icoord%gc%rad)
            z_bkg_vn = -grid_angular_velocity * z_r * z_cos_lat * z_en%v1

            ! Compute initial wind perturbation of spherical sound wave
            z_dist_rel = ( z_dist - lahade%ptb_rad_min ) / ( lahade%ptb_rad_max - lahade%ptb_rad_min ) 

            ! (The initial sound wave field is actually based on the specification 
            ! of a radial wind field perturbation. That is the reason, why the following expression 
            ! looks much simpler, than the corresponding expression for the scalar quantities above, 
            ! which are derived from this expression)
            IF ( z_dist_rel >= 0._wp .AND. z_dist_rel <= 1._wp ) THEN
              z_ptb_vn = z_ptb_amp_vr * SIN( pi * z_dist_rel ) *              & 
                &        SIN( 2._wp * pi * lahade%ptb_n_rad * z_dist_rel ) *  & 
                &        z_distvec_dot_en
            ELSE
              z_ptb_vn = 0._wp
            ENDIF

            ! Initialize vn
            p_nh_prog%vn(je,jk,jb) = z_bkg_vn + z_ptb_vn

          ENDDO  ! je
        ENDDO  ! jk
      ENDDO  !jb
!$OMP ENDDO

      !-----------------------------------------------------------
      ! 3rd Initialize vertical wind component at cell interfaces
      !-----------------------------------------------------------

      i_rlstart  = 1
      i_rlend    = min_rlcell
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_pos_cc,z_pos_gc,z_distvec_cc,z_dist,  & 
!$OMP            z_sin_lat,z_sin_lon,z_cos_lat,z_cos_lon,z_distvec_dot_er,z_dist_rel,z_ptb_w) 
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx 

            ! Geographical coordinates of current position
            z_pos_gc(icoord%gc%lon) = p_patch%cells%center(jc,jb)%lon                 ! Zonal
            z_pos_gc(icoord%gc%lat) = p_patch%cells%center(jc,jb)%lat                 ! Meridional
            z_pos_gc(icoord%gc%rad) = grid_sphere_radius + p_metrics%z_ifc(jc,jk,jb)  ! Radial

            ! Transform geographical to Cartesian coordinates
            z_pos_cc = coord_trafo( z_pos_gc, icoord%trafo%gc2cc )

            ! Normalized distance vector to center of sound wave
            z_distvec_cc = z_pos_cc - z_ptb_ctr_cc 
            z_dist       = SQRT( DOT_PRODUCT( z_distvec_cc, z_distvec_cc ) )
            z_distvec_cc = z_distvec_cc / z_dist      

            ! Scalar product of distance vector and radial unit vector e_r.
            ! (The spherical unit vectors measured in Cartesian coordinates are: 
            ! e_lon =          -sin(lon)*e_x +          cos(lon)*e_y, 
            ! e_lat = -sin(lat)*cos(lon)*e_x - sin(lat)*sin(lon)*e_y + cos(lat)*e_z, 
            ! e_r   =  cos(lat)*cos(lon)*e_x + cos(lat)*sin(lon)*e_y + sin(lat)*e_z ) 
            z_sin_lat = SIN( p_patch%cells%center(jc,jb)%lat )
            z_sin_lon = SIN( p_patch%cells%center(jc,jb)%lon )
            z_cos_lat = COS( p_patch%cells%center(jc,jb)%lat )
            z_cos_lon = COS( p_patch%cells%center(jc,jb)%lon )
            z_distvec_dot_er = z_distvec_cc(icoord%cc%x) * z_cos_lat * z_cos_lon  &
              &              + z_distvec_cc(icoord%cc%y) * z_cos_lat * z_sin_lon  & 
              &              + z_distvec_cc(icoord%cc%z) * z_sin_lat

            ! Compute initial wind perturbation of spherical sound wave.
            ! Note: we assume zero topography!
            z_dist_rel = ( z_dist - lahade%ptb_rad_min ) / ( lahade%ptb_rad_max - lahade%ptb_rad_min ) 
            
            IF ( z_dist_rel >= 0._wp .AND. z_dist_rel <= 1._wp ) THEN
              z_ptb_w = z_ptb_amp_vr * SIN( pi * z_dist_rel ) *               &  
                &        SIN( 2._wp * pi * lahade%ptb_n_rad * z_dist_rel ) *  & 
                &        z_distvec_dot_er
            ELSE
              z_ptb_w  = 0._wp
            ENDIF

            ! Initialize w
            p_nh_prog%w(jc,jk,jb) = z_ptb_w
            
          ENDDO  !jc
        ENDDO  !jk

        ! Zero vertical wind at bottom and lid
        DO jc = i_startidx, i_endidx 

          p_nh_prog%w(jc,1,jb)      = 0._wp
          p_nh_prog%w(jc,nlevp1,jb) = 0._wp

        ENDDO  !jc
      ENDDO  !jb
!$OMP ENDDO
!$OMP END PARALLEL

      ! Print some information
      IF (lmessage .AND. jg==1) THEN
        WRITE(message_text, '(a)') 'Some info:'
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f10.2)') 'Radius of Earth [m]: ', earth_radius
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f10.2)') 'Radius of model Earth [m]: ', grid_sphere_radius
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.10)') 'Angular velocity of Earth [rad/s]: ', earth_angular_velocity
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.10)') 'Angular velocity of model Earth [rad/s]: ', grid_angular_velocity
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,E10.4)') 'Gravitational acceleration of model Earth [m/s2]: ', grav
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.7)') 'Background temperature [K]: ', z_bkg_temp
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.5)') 'Background pressure [Pa]: ', z_bkg_pres
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.9)') 'Background density [kg/m3]: ', z_bkg_rho
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.7)') 'Max. magnitude of background wind [m/s]: ', z_bkg_vn_max
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.9)') 'Sound wave temperature amplitude [K]: ', z_ptb_amp_temp
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f16.9)') 'Sound wave pressure amplitude [Pa]: ', z_ptb_amp_pres
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.9)') 'Sound wave density amplitude [kg/m3]: ', z_ptb_amp_rho
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.9)') 'Sound wave velocity amplitude [m/s]: ', z_ptb_amp_vr
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.7)') 'Speed of sound [m/s]: ', z_cs
        CALL message(TRIM(routine),TRIM(message_text))
        
        WRITE(message_text, '(a,f12.7)') 'Model time step [s]: ', upatmo_config%dt_fastphy
        CALL message(TRIM(routine),TRIM(message_text))

        WRITE(message_text, '(a,f12.7)') 'Actual model dynamics time step [s]: ', z_dtime_dyn_given
        CALL message(TRIM(routine),TRIM(message_text))

        WRITE(message_text, '(a,f12.7)') 'Upper bound for model dynamics time step [s]: ', z_dtime_dyn_ubound
        CALL message(TRIM(routine),TRIM(message_text))
      ENDIF

      IF (lmessage) CALL message(TRIM(routine), &
        & 'Setup of lahade-testcase: "spherical sound wave" on domain '//TRIM(int2string(jg))//' finished')
      
    CASE default
      CALL finish(TRIM(routine),'invalid lahade%icase')
    END SELECT

    !-----------------------------------------------------------
    ! 4th Miscellaneous
    !-----------------------------------------------------------

    ! For the sound wave testcases it might be desirable to have "full" control 
    ! over the off-centering parameter for the vertical wind solver: 'nonhydrostatic_nml: vwind_offctr'. 
    ! Unfortunately, the empirical modification of the off-centering parameter 
    ! in 'src/atm_dyn_iconam/mo_vertical_grid' does not allow for a negative off-centering 
    ! (i.e., a stronger weighting of the explicit part as compared to the implicit part). 
    ! For this reason we recompute it here. (The call of the testcase initialization 
    ! takes place after the call of 'set_nh_metrics', within which the off-centering is modified, 
    ! in 'src/atm_dyn_iconam/mo_nh_stepping: prepare_nh_integration'). 
    p_metrics%vwind_impl_wgt(:,:) = 0.5_wp + vwind_offctr   ! ((:,:) -> (jc,jb))
    p_metrics%vwind_expl_wgt(:,:) = 0.5_wp - vwind_offctr

    ! In case the output of the analytical solution of a perturbation quantity is required, 
    ! we should initialize it here for time = 0
    IF (lahade%lupdate) THEN
      CALL nh_lahade_interface ( jstep              = 0,                 &  !in
        &                        sim_time           = 0._wp,             &  !in
        &                        p_patch            = p_patch,           &  !in
        &                        p_metrics          = p_metrics,         &  !in
        &                        p_nh_prog          = p_nh_prog,         &  !in
        &                        p_nh_diag          = p_nh_diag,         &  !inout
        &                        grid_sphere_radius = grid_sphere_radius )  !in
    ENDIF

  END SUBROUTINE init_nh_lahade

  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !>
  !! Updates for the lahade-testcase during integration.
  !! This subroutine is called in 'src/testcases/mo_nh_testcase_interface: nh_testcase_interface'.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2018-06-22)
  !!
  SUBROUTINE nh_lahade_interface ( jstep,             &  !in
    &                              sim_time,          &  !in
    &                              p_patch,           &  !in
    &                              p_metrics,         &  !in
    &                              p_nh_prog,         &  !in
    &                              p_nh_diag,         &  !inout
    &                              grid_sphere_radius )  !in

    ! In/out variables
 
    INTEGER,                    INTENT(IN)    :: jstep               !< Global time step
    REAL(wp),                   INTENT(IN)    :: sim_time            !< Elapsed simulation time on this grid level
    TYPE(t_patch),      TARGET, INTENT(IN)    :: p_patch             !< Grid/patch info
    TYPE(t_nh_metrics),         INTENT(IN)    :: p_metrics           !< NH metrics state
    TYPE(t_nh_prog),            INTENT(IN)    :: p_nh_prog           !< Prognostic state vector
    TYPE(t_nh_diag),            INTENT(INOUT) :: p_nh_diag           !< Diagnostic state vector
    REAL(wp),                   INTENT(IN)    :: grid_sphere_radius  !< (Rescaled) Earth radius

    ! Local variables

    REAL(wp), DIMENSION(3) :: z_ptb_ctr_gc, z_ptb_ctr_cc, z_pos_gc, z_pos_cc
    REAL(wp), DIMENSION(3) :: z_ptb_ctr_rot_cc, z_rot_axis_cc

    REAL(wp) :: z_dist, z_dist_rel, z_rot_angle
    REAL(wp) :: z_bkg_pres, z_bkg_temp, z_bkg_rho, z_cs
    REAL(wp) :: z_ptb_amp_var, z_ptb_amp_temp, z_ptb_amp_vr

    INTEGER :: jg, jc, jk, jb
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER :: nlev

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':nh_lahade_interface'

    !----------------------------------------------------------------------------

    ! So far, the only update is, to get the numerical solution and to compute 
    ! the analytical solution of a selected perturbation quantity 
    ! for the current model time.

    ! Domain index
    jg = p_patch%id

    ! The analytical solution is required only for data output, 
    ! so we limit the following computations to those model time steps 'jstep', 
    ! when some kind of data output is due 
    ! (compare 'src/atm_dyn_iconam/mo_nh_stepping: perform_nh_timeloop')
    IF (jstep == 0 .OR. istime4name_list_output(jstep)) THEN

      ! Different variants of lahade-testcases
      SELECT CASE(lahade%icase)
      CASE(ilahade%case%ssw)         
        ! Spherical sound wave (currently the only lahade-testcase): 
        ! gravity is switched off -> air is homogeneous gas trapped 
        ! between concentrical spherical bottom shell and spherical lid shell

        IF (msg_level >= imsg_thr%high) CALL message(TRIM(routine), &
          & 'Update of lahade-testcase: "spherical sound wave" on domain '//TRIM(int2string(jg))//' started') 
        
        !----------------------------------------------
        !                 Preparation
        !----------------------------------------------
        
        ! Center of spherical sound wave in geographical coordinates at time = 0
        z_ptb_ctr_gc(icoord%gc%lon) = lahade%ptb_ctr_lon                       ! Zonal
        z_ptb_ctr_gc(icoord%gc%lat) = lahade%ptb_ctr_lat                       ! Meridional
        z_ptb_ctr_gc(icoord%gc%rad) = grid_sphere_radius + lahade%ptb_ctr_hgt  ! Radial
        
        ! And Cartesian coordinates of center
        z_ptb_ctr_cc = coord_trafo( z_ptb_ctr_gc, icoord%trafo%gc2cc )
        
        ! Compute rotated center of sound wave at time = 'sim_time':
        ! the angular velocity vector is: omega*e_z, 
        ! the rotation of a radial unit vector e_r around e_z with angular velocity omega 
        ! within the time span t should then follow from: 
        !
        ! e_r(t) = cos(omega*t)*e_r(t=0) - sin(omega*t)*[e_z x e_r(t=0)] 
        !        + {[1 - cos(omega*t)]*[e_z . e_r(t=0)]}*e_z, 
        !
        ! where [a x b] denotes the cross product of a and b, 
        ! and [a . b] is the scalar product of a and b.

        ! Rotation axis in Cartesian coordinates => e_z
        z_rot_axis_cc(icoord%cc%x) = 0._wp
        z_rot_axis_cc(icoord%cc%y) = 0._wp
        z_rot_axis_cc(icoord%cc%z) = 1._wp 

        ! New position of sound wave center at time = 'sim_time' => e_r(t=0) -> e_r(t)
        z_rot_angle      = lahade%omega * sim_time
        z_ptb_ctr_rot_cc = COS( z_rot_angle ) * z_ptb_ctr_cc                             &  
          &              - SIN( z_rot_angle ) *                                          &
          &                cross_product( z_rot_axis_cc, z_ptb_ctr_cc, icoord%id%cc )    &
          &              + ( ( 1._wp - COS( z_rot_angle ) ) *                            & 
          &                DOT_PRODUCT( z_rot_axis_cc, z_ptb_ctr_cc ) ) * z_rot_axis_cc
        
        ! Background state
        z_bkg_pres = lahade%bkg_pres                   ! Pressure
        z_bkg_temp = lahade%bkg_temp                   ! Temperature
        z_bkg_rho  = z_bkg_pres / ( rd * z_bkg_temp )  ! Density
        
        ! Speed of sound 
        z_cs = SQRT( ( cpd / cvd ) * rd * z_bkg_temp )
        
        ! Sound wave amplitude with respect to the perturbation variable selected for output:
        ! auxiliary amplitudes for temperature and radial velocity
        z_ptb_amp_temp = lahade%ptb_amp_temp                                    
        z_ptb_amp_vr   = ( cvd / rd ) * ( z_ptb_amp_temp / z_bkg_temp ) * z_cs  
        SELECT CASE(lahade%ivarout)
        CASE(ilahade%varout%temp)
          ! Temperature
          z_ptb_amp_var = z_ptb_amp_temp
        CASE(ilahade%varout%rho)
          ! Density
          z_ptb_amp_var = ( z_ptb_amp_vr / z_cs ) * z_bkg_rho   
        CASE(ilahade%varout%pres)
          ! Pressure
          z_ptb_amp_var = ( cpd / cvd ) * ( z_ptb_amp_vr / z_cs ) * z_bkg_pres
        CASE default
          CALL finish(TRIM(routine),'invalid lahade%ivarout')
        END SELECT
        
        ! Looping
        nlev       = p_patch%nlev
        i_rlstart  = 1
        i_rlend    = min_rlcell
        i_startblk = p_patch%cells%start_block(i_rlstart)
        i_endblk   = p_patch%cells%end_block(i_rlend)
        
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_pos_cc,z_pos_gc,z_dist,z_dist_rel)
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, i_rlstart, i_rlend)

          !----------------------------------------------
          !         Get numerical solution
          !----------------------------------------------

          ! (Since these computations take place at output times only, 
          ! we assume a 'select-case' within the jb-loop to be computationally bearable)
          IF (lahade%ivarout == ilahade%varout%temp) THEN
            ! Temperature perturbation = theta_v * exner - background temperature
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx 
                p_nh_diag%extra_3d(jc,jk,jb,isol%num) = &
                  & p_nh_prog%theta_v(jc,jk,jb) * p_nh_prog%exner(jc,jk,jb) - z_bkg_temp
              ENDDO  !jc
            ENDDO  !jk
          ELSEIF (lahade%ivarout == ilahade%varout%rho) THEN
            ! Density perturbation = rho - background density
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx 
                p_nh_diag%extra_3d(jc,jk,jb,isol%num) = &
                  & p_nh_prog%rho(jc,jk,jb) - z_bkg_rho
              ENDDO  !jc
            ENDDO  !jk
          ELSEIF (lahade%ivarout == ilahade%varout%pres) THEN
            ! Pressure perturbation = (rho * rd * temp = rho * rd * theta_v * exner) - background pressure
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx 
                p_nh_diag%extra_3d(jc,jk,jb,isol%num) = &
                  & rd * p_nh_prog%rho(jc,jk,jb) * p_nh_prog%theta_v(jc,jk,jb) * &
                  & p_nh_prog%exner(jc,jk,jb) - z_bkg_pres
              ENDDO  !jc
            ENDDO  !jk
          ENDIF

          !----------------------------------------------
          !         Compute analytical solution
          !----------------------------------------------
          
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx 

              ! Geographical coordinates of current position 
              z_pos_gc(icoord%gc%lon) = p_patch%cells%center(jc,jb)%lon                ! Zonal
              z_pos_gc(icoord%gc%lat) = p_patch%cells%center(jc,jb)%lat                ! Meridional
              z_pos_gc(icoord%gc%rad) = grid_sphere_radius + p_metrics%z_mc(jc,jk,jb)  ! Radial

              ! Transform geographical to Cartesian coordinates
              z_pos_cc = coord_trafo( z_pos_gc, icoord%trafo%gc2cc )
              
              ! Distance of current grid point to sound wave center
              z_dist = SQRT( DOT_PRODUCT( ( z_pos_cc - z_ptb_ctr_rot_cc ), ( z_pos_cc - z_ptb_ctr_rot_cc ) ) )
              
              IF ( z_dist >= lahade%ptb_rad_min + z_cs * sim_time .AND. &
                &  z_dist <= lahade%ptb_rad_max + z_cs * sim_time       ) THEN
                
                z_dist_rel = ( z_dist - lahade%ptb_rad_min - z_cs * sim_time ) / &
                  &          ( lahade%ptb_rad_max - lahade%ptb_rad_min ) 
                
                ! Compute sound wave amplitude at current grid point 
                p_nh_diag%extra_3d(jc,jk,jb,isol%ana) = z_ptb_amp_var * (                               &
                  &         ( z_dist - z_cs * sim_time ) / z_dist *                                     &
                  &         SIN( pi * z_dist_rel ) * SIN( 2._wp * pi * lahade%ptb_n_rad * z_dist_rel )  & 
                  &         + ( lahade%ptb_rad_max - lahade%ptb_rad_min ) / z_dist *                    & 
                  &         ( SIN( pi * ( 2._wp * lahade%ptb_n_rad - 1._wp ) * z_dist_rel ) /           & 
                  &         ( 2._wp * pi * ( 2._wp * lahade%ptb_n_rad - 1._wp ) )                       & 
                  &         - SIN( pi * ( 2._wp * lahade%ptb_n_rad + 1._wp ) * z_dist_rel ) /           & 
                  &         ( 2._wp * pi * ( 2._wp * lahade%ptb_n_rad + 1._wp ) ) )                     &
                  &                                                     )
              ELSE
                p_nh_diag%extra_3d(jc,jk,jb,isol%ana) = 0._wp
              ENDIF
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP ENDDO
!$OMP END PARALLEL

        IF (msg_level >= imsg_thr%high) CALL message(TRIM(routine), &
          & 'Update of lahade-testcase: "spherical sound wave" on domain '//TRIM(int2string(jg))//' finished') 

      CASE default
        CALL finish(TRIM(routine),'invalid lahade%icase')
      END SELECT
        
    ENDIF  !Is output due?

  END SUBROUTINE nh_lahade_interface

  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !>
  !! Crosscheck lahade-testcase settings with other namelist settings.
  !! This subroutine is called in 'src/testcases/mo_nh_testcase_check: check_nh_testcase'.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2017-04-21)
  !!
  SUBROUTINE check_nh_lahade ( lvert_nest,             &  !in
    &                          l_nml,                  &  !in
    &                          ivctype,                &  !in
    &                          top_height,             &  !in 
    &                          grid_sphere_radius,     &  !in
    &                          grid_rescale_factor,    &  !in
    &                          first_output_name_list, &  !in
    &                          ldeepatmo,              &  !inout
    &                          lcoriolis,              &  !inout
    &                          l_open_ubc,             &  !inout
    &                          ltransport,             &  !inout
    &                          itopo,                  &  !inout
    &                          iforcing,               &  !inout
    &                          inextra_3d,             &  !inout
    &                          min_lay_thckn,          &  !inout
    &                          grid_angular_velocity,  &  !inout
    &                          upatmo_dyn_config,      &  !inout
    &                          upatmo_config           )  !(opt)in
   
    ! In/out variables

    LOGICAL,                  INTENT(IN)    :: lvert_nest             ! Switch for vertical grid nesting
    LOGICAL,                  INTENT(IN)    :: l_nml                  ! Switch for output
    INTEGER,                  INTENT(IN)    :: ivctype                ! Type of vertical grid (SLEVE etc.)
    REAL(wp),                 INTENT(IN)    :: top_height             ! Height of model top
    REAL(wp),                 INTENT(IN)    :: grid_sphere_radius     ! (Rescaled) Earth radius
    REAL(wp),                 INTENT(IN)    :: grid_rescale_factor    ! Rescale factor
    TYPE(t_output_name_list), POINTER       :: first_output_name_list ! Pointer to a linked list 
                                                                      ! of output name lists    
    LOGICAL,                  INTENT(INOUT) :: ldeepatmo              ! Switch for deep-atmosphere dynamics
    LOGICAL,                  INTENT(INOUT) :: lcoriolis              ! Switch for Coriolis acceleration
    LOGICAL,                  INTENT(INOUT) :: l_open_ubc             ! Switch for open upper boundary condition
    LOGICAL,                  INTENT(INOUT) :: ltransport             ! Switch for tracer transport
    INTEGER,                  INTENT(INOUT) :: itopo                  ! Type of topography
    INTEGER,                  INTENT(INOUT) :: iforcing               ! Switch for physics package 
                                                                      ! (nwp, echam etc.)  
    INTEGER,                  INTENT(INOUT) :: inextra_3d             ! Number of extra fields
    REAL(wp),                 INTENT(INOUT) :: min_lay_thckn          ! Layer thickness of lowermost grid layer
    REAL(wp),                 INTENT(INOUT) :: grid_angular_velocity  ! (Rescaled) Earth angular velocity
    TYPE(t_upatmo_dyn_config),INTENT(INOUT) :: upatmo_dyn_config(:)   ! Deep-atmosphere configuration
    TYPE(t_upatmo_config), OPTIONAL, INTENT(IN) :: upatmo_config(:)   ! Upper-atmosphere configuration

    ! Local variables

    REAL(wp) :: z_ptb_ctr_hgt, z_ptb_rad_max
    REAL(wp) :: z_r, z_cos_lat
    REAL(wp) :: z_small_number

    LOGICAL  :: lupatmo_checked, lcentrifugal

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':check_nh_lahade'

    !----------------------------------------------------------------------------
    
    ! This subroutine is calle only if 
    ! ltestcase = .true. and nh_test_name = 'lahade', 
    ! so we can conclude: 
    lahade%lactive = .TRUE. 

    SELECT CASE(lahade%icase)

    CASE(ilahade%case%ssw)    ! Spherical sound wave (currently the only sub-case) 

      !---------------------------------------
      ! 1st Crosscheck other namelist entries
      !---------------------------------------

      IF (PRESENT(upatmo_config)) THEN
        lupatmo_checked = ANY(upatmo_config(:)%l_status(istatus%checked))
      ELSE
        lupatmo_checked = .FALSE.
      ENDIF
      lcentrifugal = ALL(upatmo_dyn_config(:)%lcentrifugal)
        
      ! This subroutine should have been called before 'check_upatmo'
      IF (lupatmo_checked) THEN
        CALL finish(TRIM(routine),'check_upatmo should be called after this subroutine')
      ENDIF
      ! Deep-atmosphere dynamics have to be switched on
      IF (.NOT. ldeepatmo) THEN
        ldeepatmo = .TRUE.
        CALL message(TRIM(routine),'WARNING! ldeepatmo set to .true.')
      ENDIF
      ! Since the model top height 'top_height' is required in the following, 
      ! only the vertical SLEVE-coordinates are allowed, 
      ! and since this is a more sever namelist switch, 
      ! it seems more secure to stop the program, if necessary
      IF (ivctype/=2) CALL finish(TRIM(routine),'only SLEVE-coordinates (ivctype=2) are allowed for this testcase')
      ! For this testcase a const. vertical layer thickness seems reasonable
      ! (see 'atm_dyn_iconam/mo_init_vgrid: init_sleve_coord')
      IF (min_lay_thckn > 0._wp) THEN
        min_lay_thckn = 0._wp
        CALL message(TRIM(routine),'WARNING! min_lay_thckn set to 0')
      ENDIF
      ! No gravitational acceleration
      lahade%lzerograv = .TRUE. 
      ! Coriolis acceleration has to be switched on
      IF (.NOT. lcoriolis) THEN
        lcoriolis = .TRUE.
        CALL message(TRIM(routine),'WARNING! lcoriolis set to .true.')
      ENDIF      
      ! Centrifugal acceleration has to be switched on 
      IF (.NOT. lcentrifugal) THEN
        upatmo_dyn_config(:)%lcentrifugal = .TRUE.
        CALL message(TRIM(routine),'WARNING! lcentrifugal set to .true.')
      ENDIF
      ! Rigid-lid boundary condition at model top
      IF (l_open_ubc) THEN
        l_open_ubc = .FALSE.
        CALL message(TRIM(routine),'WARNING! l_open_ubc set to .false.')
      ENDIF
      ! Only analytical topography
      IF (itopo /= 0) THEN 
        itopo = 0
        CALL message(TRIM(routine),'WARNING! itopo set to 0')
      ENDIF
      ! No physics forcing
      IF (iforcing /= inoforcing) THEN 
        iforcing = inoforcing
        CALL message(TRIM(routine),'WARNING! iforcing set to 0')
      ENDIF
      ! No tracer transport
      IF (ltransport) THEN 
        ltransport = .FALSE.
        CALL message(TRIM(routine),'WARNING! ltransport set to .false.')
      ENDIF
      ! For this testcase no vertical nesting should take place, 
      ! and since this is a more sever namelist switch, 
      ! it seems more secure to stop the program, if necessary
      IF (lvert_nest) CALL finish(TRIM(routine),'vertical grid nesting (lvert_nest=.true.) not allowed for this testcase')

      !-------------------------------------------------------------------------------------------------
      ! 2nd Unit conversion of namelist entries to which it applies
      ! (To do this here is without alternative, because of the effects on 'grid_angular_velocity' etc.)
      !-------------------------------------------------------------------------------------------------

      ! Note: the units in which the testcase parameters have to be entered in the namelist 
      ! look ugly at first glance, but they were chosen this way, in order to save a lot of 
      ! model-run-prior "pocket calculator" work, if standard-units would have been used instead

      ! Center of spherical sound wave in geographical coordinates (from unit [deg] to [rad])
      lahade%ptb_ctr_lon = lahade%ptb_ctr_lon * deg2rad
      lahade%ptb_ctr_lat = lahade%ptb_ctr_lat * deg2rad      
      ! Height of center of spherical sound wave (from unit [top_height] to [m])
      z_ptb_ctr_hgt      = lahade%ptb_ctr_hgt
      lahade%ptb_ctr_hgt = z_ptb_ctr_hgt * top_height
      ! Inner radius of sound wave (for radii smaller than this radius the initial amplitude 
      ! of the sound wave is zero) (from unit [min{ptb_ctr_hgt,(1-ptb_ctr_hgt)} * top_height] to [m])
      lahade%ptb_rad_min = lahade%ptb_rad_min * MIN( z_ptb_ctr_hgt, 1._wp - z_ptb_ctr_hgt ) * top_height
      ! Outer radius of sound wave (for radii greater than this radius the initial amplitude 
      ! of the sound wave is zero) (from unit [min{ptb_ctr_hgt,(1-ptb_ctr_hgt)} * top_height] to [m])
      z_ptb_rad_max      = lahade%ptb_rad_max
      lahade%ptb_rad_max = z_ptb_rad_max * MIN( z_ptb_ctr_hgt, 1._wp - z_ptb_ctr_hgt ) * top_height
      ! Angular velocity of the testcase Earth:
      ! cosine of latitute of sound wave center
      z_cos_lat = COS( lahade%ptb_ctr_lat )
      IF (ABS(z_cos_lat) < dbl_eps*1000._wp) THEN
        ! The center of the spherically symmetric sound wave is located at one of the poles, 
        ! and, ideally, its evolution should appear the same in the rotational and in the absolute frame, 
        ! so we set the angular velocity to zero
        lahade%omega = 0._wp        
      ELSE
        ! Radial position of center of sound wave 
        ! (This subroutine should be called at a time during the program sequence, 
        ! when 'grid_sphere_radius' has already been rescaled, see below)
        z_r = grid_sphere_radius + lahade%ptb_ctr_hgt
        ! The velocity of the center should be: 'omega * z_r * z_cos_lat', 
        ! and its value should be contained in 'lahade%omega' from the namelist entry, 
        ! so we can transform 'lahade%omega' into the actual angular velocity 'omega' via
        lahade%omega = lahade%omega / (z_r * z_cos_lat)
      ENDIF
      grid_angular_velocity = lahade%omega
      
      !--------------------------------------
      ! 3rd Check gravitational acceleration
      !--------------------------------------

      z_small_number = 0.0001_wp

      IF (lahade%lzerograv .AND. grav > z_small_number) THEN
        WRITE(message_text, '(a)') 'This testcase requires grav=0. ' // &
          & 'Please, set this parameter in src/shared/mo_physical_constants to a very small value and recompile.'
        CALL finish(TRIM(routine),TRIM(message_text))
      ENDIF

      !---------------------------------------------------
      ! 4th Check if grid_sphere_radius has been rescaled
      !---------------------------------------------------

      IF (ABS(grid_sphere_radius/earth_radius - grid_rescale_factor) > z_small_number) THEN
        WRITE(message_text, '(a)') 'Wrong subroutine sequence, ' // &
          & 'grid_sphere_radius has not been rescaled yet.'
        CALL finish(TRIM(routine),TRIM(message_text))
      ENDIF

    CASE default
      CALL finish(TRIM(routine),'invalid lahade%icase')
    END SELECT

    ! Check, if the output of the analytical solution of a perturbation quantity is required
    IF (LEN_TRIM(lahade%output_ptb_var) > 0) THEN
      SELECT CASE(TRIM(lahade%output_ptb_var))
      CASE("temp")
        ! Output analytical solution for temperature perturbation
        lahade%ivarout = ilahade%varout%temp
        CALL message(TRIM(routine), 'Temperature perturbations have been selected for output &
          &in extra_3d1 (numerical solution) and extra_3d2 (analytical solution).')
      CASE("rho")
        ! Output analytical solution for density perturbation
        lahade%ivarout = ilahade%varout%rho
        CALL message(TRIM(routine), 'Density perturbations have been selected for output &
          &in extra_3d1 (numerical solution) and extra_3d2 (analytical solution).')
      CASE("pres")
        ! Output analytical solution for pressure perturbation
        lahade%ivarout = ilahade%varout%pres
        CALL message(TRIM(routine), 'Pressure perturbations have been selected for output &
          &in extra_3d1 (numerical solution) and extra_3d2 (analytical solution).')
      CASE default
        CALL finish(TRIM(routine),'invalid lahade%output_ptb_var')
      END SELECT
      ! The diagnostic field 'p_diag%extra_3d(nproma,nlev,nblks_c,inextra_3d)' has to be allocated: 
      ! - 'extra_3d1': numerical solution for the sound-wave-perturbation-variable
      ! - 'extra_3d2': analytical solution for the sound-wave-perturbation-variable
      IF (inextra_3d > 0) THEN
        CALL message(TRIM(routine),'WARNING! your entry for inextra_3d has been overwritten by 2.')
      ENDIF
      inextra_3d = 2
      ! 'extra_3d1' and 'extra_3d2' should have been added to the 'output_nml: ml_varlist' of choice
      IF (.NOT. l_nml) THEN
        CALL finish(TRIM(routine),'an entry to lahade%output_ptb_var requires &
          &run_nml: output = ..., "nml", ...')
      ELSEIF ((.NOT. is_variable_in_output(first_output_name_list, var_name="extra_3d1")) .OR. &
        &     (.NOT. is_variable_in_output(first_output_name_list, var_name="extra_3d2"))      ) THEN
        CALL finish(TRIM(routine),'an entry to lahade%output_ptb_var requires &
          &output_nml: ml_varlist = ..., "extra_3d1", "extra_3d2",...')
      ENDIF
      ! The computation of the analytical solution requires 'nh_lahade_interface' 
      ! to be called during integration
      lahade%lupdate = lahade%lactive
    ENDIF

  END SUBROUTINE check_nh_lahade

  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !>
  !! Function to transform between coordinates.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2017-04-21)
  !! - Adaption of 'src/shared/mo_math_utilities: gc2cc', 
  !!   which is intended for input/output of data types 
  !!   'src/shared/mo_math_types: t_geographical_coordinates' and 
  !!   'src/shared/mo_math_types: t_cartesian_coordinates'
  !!
  FUNCTION coord_trafo(coords_in, itrafotype) RESULT(coords_out)

    REAL(wp), DIMENSION(3), INTENT(IN) :: coords_in    ! Coordinates before transformation
    INTEGER,                INTENT(IN) :: itrafotype   ! Transformation identifier
    REAL(wp), DIMENSION(3)             :: coords_out   ! Coordinates after transformation
    REAL(wp) :: z_sin_lon, z_sin_lat
    REAL(wp) :: z_cos_lon, z_cos_lat
    REAL(wp) :: z_r, z_r_h
    REAL(wp) :: z_x, z_y, z_z
    REAL(wp), PARAMETER :: z_small_number = EPSILON(1._wp) * 1000._wp
    REAL(wp), PARAMETER :: z_huge_number  = HUGE(1._wp) / 1000._wp

    !--------------------------------------------

    SELECT CASE(itrafotype)
    CASE(icoord%trafo%gc2cc)
      ! Geographical to Cartesian coordinates
      z_sin_lon = SIN(coords_in(icoord%gc%lon))
      z_sin_lat = SIN(coords_in(icoord%gc%lat))
      z_cos_lon = COS(coords_in(icoord%gc%lon))
      z_cos_lat = COS(coords_in(icoord%gc%lat))
      z_r       = coords_in(icoord%gc%rad)
      !
      coords_out(icoord%cc%x) = z_cos_lon * z_cos_lat * z_r 
      coords_out(icoord%cc%y) = z_sin_lon * z_cos_lat * z_r
      coords_out(icoord%cc%z) = z_sin_lat * z_r
    CASE(icoord%trafo%cc2gc)
      ! Cartesian to geographical coordinates
      z_x   = coords_in(icoord%cc%x)
      z_y   = coords_in(icoord%cc%y)
      z_z   = coords_in(icoord%cc%z)
      z_r_h = SQRT( z_x**2 + z_y**2 )
      z_r   = SQRT( z_x**2 + z_y**2 + z_z**2 )
      !
      IF (z_r < z_small_number) THEN
        ! The radius is so small that we assume 
        ! the position to coincide with the center of the Earth
        coords_out(icoord%gc%lon) = 0._wp  ! (Actually not defined)
        coords_out(icoord%gc%lat) = 0._wp  ! (Actually not defined)
        coords_out(icoord%gc%rad) = 0._wp
      ELSEIF (z_r_h < z_small_number) THEN
        ! The length of the projection of the radius vector 
        ! onto the horizontal plane is so small that we assume 
        ! the position to be somewhere on the z-axis
        coords_out(icoord%gc%lon) = 0._wp  ! (Actually not defined)
        coords_out(icoord%gc%lat) = ASIN( z_z / z_r )
        coords_out(icoord%gc%rad) = z_r
      ELSE
        ! (There are likely cleaner ways to compute the longitude)
        coords_out(icoord%gc%lon) = SIGN( ACOS( z_x / z_r_h ), z_y )
        coords_out(icoord%gc%lat) = ASIN( z_z / z_r )
        coords_out(icoord%gc%rad) = z_r
      ENDIF
    CASE default
      ! 'error handling' 
      ! (This function might be called within OMP-threading, 
      ! so that a call of 'finish' might not be recommended, unfortunately. 
      ! Instead we hand over a relatively large number, which hopefully 
      ! points the user to the right direction.)
      coords_out = z_huge_number
    END SELECT

  END FUNCTION coord_trafo

  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------

  !>
  !! Function to compute the cross product of two 'position' vectors 
  !! from coordinate tripls
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2017-04-21)
  !! - Adaption of 'src/shared/mo_math_utilities: vector_product', 
  !!   which is intended for input/output of data type 
  !!   'src/shared/mo_math_types: t_cartesian_coordinates'
  !!
  FUNCTION cross_product(coords_fac1, coords_fac2, icoordtype) RESULT(coords_prod) 
    
    REAL(wp), DIMENSION(3), INTENT(IN) :: coords_fac1, coords_fac2
    INTEGER,                INTENT(IN) :: icoordtype
    REAL(wp), DIMENSION(3)             :: coords_prod
    REAL(wp), DIMENSION(3)             :: coords_aux1, coords_aux2, coords_aux3
    REAL(wp), PARAMETER :: z_huge_number  = HUGE(1._wp) / 1000._wp

    !--------------------------------------------
    
    SELECT CASE(icoordtype)
    CASE(icoord%id%gc)
      ! If we are faced with a coordinate triple in geographical coordinates, 
      ! it might be simpler to first transform it into Cartesian coordinates ...
      coords_aux1 = coord_trafo(coords_fac1, icoord%trafo%gc2cc)
      coords_aux2 = coord_trafo(coords_fac2, icoord%trafo%gc2cc)
      ! ... compute the cross product ...
      coords_aux3(icoord%cc%x) = coords_aux1(icoord%cc%y) * coords_aux2(icoord%cc%z) &
        &                      - coords_aux1(icoord%cc%z) * coords_aux2(icoord%cc%y)
      coords_aux3(icoord%cc%y) = coords_aux1(icoord%cc%z) * coords_aux2(icoord%cc%x) &
        &                      - coords_aux1(icoord%cc%x) * coords_aux2(icoord%cc%z)
      coords_aux3(icoord%cc%z) = coords_aux1(icoord%cc%x) * coords_aux2(icoord%cc%y) &
        &                      - coords_aux1(icoord%cc%y) * coords_aux2(icoord%cc%x) 
      ! ... and transform the result back again
      coords_prod = coord_trafo(coords_aux3, icoord%trafo%cc2gc)
    CASE(icoord%id%cc)
      ! In Cartesian coordinates we might compute the cross product without further ado
      coords_prod(icoord%cc%x) = coords_fac1(icoord%cc%y) * coords_fac2(icoord%cc%z) &
        &                      - coords_fac1(icoord%cc%z) * coords_fac2(icoord%cc%y)
      coords_prod(icoord%cc%y) = coords_fac1(icoord%cc%z) * coords_fac2(icoord%cc%x) &
        &                      - coords_fac1(icoord%cc%x) * coords_fac2(icoord%cc%z)
      coords_prod(icoord%cc%z) = coords_fac1(icoord%cc%x) * coords_fac2(icoord%cc%y) &
        &                      - coords_fac1(icoord%cc%y) * coords_fac2(icoord%cc%x) 
    CASE default
      ! 'error handling' 
      ! (This function might be called within OMP-threading, 
      ! so that a call of 'finish' might not be recommended, unfortunately. 
      ! Instead we hand over a relatively large number, which hopefully 
      ! points the user to the right direction.)
      coords_prod = z_huge_number
    END SELECT      
    
  END FUNCTION cross_product


END MODULE mo_nh_lahade
