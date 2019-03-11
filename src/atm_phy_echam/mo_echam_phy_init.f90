!>
!! @brief Contains subroutines for initializing the ECHAM physics package.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, 2010-07-20
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_echam_phy_init

  ! infrastructure
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mtime,                   ONLY: datetime, OPERATOR(>), OPERATOR(==)
  USE mo_io_config,            ONLY: default_read_method
  USE mo_read_interface,       ONLY: openInputFile, closeFile, read_2D, &
    &                                t_stream_id, on_cells
  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, &
    &                                timer_prep_echam_phy
  USE mo_impl_constants,       ONLY: max_char_length

  ! model configuration
  USE mo_impl_constants,       ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_parallel_config,      ONLY: nproma
  USE mo_master_config,        ONLY: isrestart
  USE mo_run_config,           ONLY: iqt, io3, ntracer, ltestcase

  ! horizontal grid and indices
  USE mo_model_domain,         ONLY: t_patch
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_grid_config,          ONLY: n_dom

  ! vertical grid
  USE mo_vertical_coord_table, ONLY: vct
  
  ! test cases
  USE mo_nh_testcases_nml,     ONLY: nh_test_name, ape_sst_case, th_cbl, tpe_temp
  USE mo_ape_params,           ONLY: ape_sst
  USE mo_physical_constants,   ONLY: tmelt, Tf, albedoW, amd, amo3

  USE mo_sea_ice_nml,          ONLY: albi

  ! echam phyiscs
  USE mo_echam_phy_config,     ONLY: eval_echam_phy_config, eval_echam_phy_tc, print_echam_phy_config, &
    &                                echam_phy_config, echam_phy_tc, dt_zero
  USE mo_echam_phy_memory,     ONLY: prm_field, t_echam_phy_field, &
    &                                prm_tend,  t_echam_phy_tend

  ! radiation
  USE mo_echam_rad_config,     ONLY: eval_echam_rad_config, print_echam_rad_config, echam_rad_config

  ! subgrid scale orographic effects
  USE mo_echam_sso_config,     ONLY: eval_echam_sso_config, print_echam_sso_config

  ! atmospheric gravity wave drag
  USE mo_echam_gwd_config,     ONLY: eval_echam_gwd_config, print_echam_gwd_config

  ! vertical diffusion
  USE mo_echam_vdf_config,     ONLY: eval_echam_vdf_config, print_echam_vdf_config
  USE mo_echam_vdiff_params,   ONLY: init_vdiff_params
  USE mo_vdiff_solver,         ONLY: init_vdiff_solver

#ifndef __NO_JSBACH__
  ! land surface
  USE mo_jsb_model_init,       ONLY: jsbach_init
#endif

  ! carbon cycle
  USE mo_ccycle_config,        ONLY: print_ccycle_config

  ! cumulus convection
  USE mo_echam_cnv_config,     ONLY: alloc_echam_cnv_config, eval_echam_cnv_config, print_echam_cnv_config
  USE mo_convect_tables,       ONLY: init_convect_tables
  USE mo_echam_convect_tables, ONLY: init_echam_convect_tables => init_convect_tables

  ! cloud microphysics
  USE mo_echam_cld_config,     ONLY: print_echam_cld_config, echam_cld_config

  ! Cariolle interactive ozone scheme
  USE mo_lcariolle_externals,  ONLY: read_bcast_real_3d_wrap, &
    &                                read_bcast_real_1d_wrap, &
    &                                closeFile_wrap, openInputFile_wrap, &
    &                                get_constants
  USE mo_lcariolle_types,      ONLY: l_cariolle_initialized_o3, t_avi, t_time_interpolation

  ! water vapour production by methane oxidation
  ! and destruction by photolysis
  USE mo_methox,               ONLY: init_methox

  ! air-sea-land interface
  USE mo_echam_sfc_indices,    ONLY: nsfc_type, iwtr, iice, ilnd, init_sfc_indices

  ! for AMIP boundary conditions
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, calculate_time_interpolation_weights
  USE mo_bc_sst_sic,           ONLY: read_bc_sst_sic, bc_sst_sic_time_interpolation
  USE mo_bc_greenhouse_gases,  ONLY: read_bc_greenhouse_gases, bc_greenhouse_gases_time_interpolation, &
    &                                bc_greenhouse_gases_file_read
  USE mo_bc_aeropt_splumes,    ONLY: setup_bc_aeropt_splumes
  ! psrad
  USE mo_psrad_setup,          ONLY: psrad_basic_setup
  USE mo_psrad_interface,      ONLY: pressure_scale, droplet_scale
  USE mo_atmo_psrad_interface, ONLY: setup_atmo_2_psrad

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: init_echam_phy_params, init_echam_phy_external, init_echam_phy_field
  PUBLIC  :: init_o3_lcariolle

CONTAINS
  !>
  !! Top-level routine for initialization of ECHAM6 physics.
  !! It calls a series of subroutines to initialize tunable parameters,
  !! lookup tables, and the physics state vectors "prm_field" and "prm_tend".
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !! name change to init_echam_phy by Levi Silvers
  !!
  SUBROUTINE init_echam_phy_params( p_patch )

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(:)

    INTEGER :: khydromet, ktrac
    INTEGER :: jg
    INTEGER :: nlev

    LOGICAL :: lany

    ! Shortcuts to components of echam_rad_config
    !
    IF (timers_level > 1) CALL timer_start(timer_prep_echam_phy)

    !-------------------------------------------------------------------
    ! Initialize parameters and lookup tables
    !-------------------------------------------------------------------

    nlev = p_patch(1)%nlev

    ! Diagnostics (all time steps)
    ! ----------------------------
    
    ! For surface processes
    ! nsfc_type, iwtr, etc. are set in this subroutine.
    ! See mo_sfc_indices.f90 for further details.
    !
    CALL init_sfc_indices( nh_test_name )
    
    ! Lookup tables for saturation vapour pressure
    !
    CALL init_convect_tables
    CALL init_echam_convect_tables 


    ! ECHAM physics time control
    ! --------------------------

    ! Evaluate the ECHAM physics configuration variables echam_phy_config(:)
    ! and the derived time control variables echam_phy_tc(:) on all grids
    ! and for all controled processes.
    !
    CALL  eval_echam_phy_config
    CALL  eval_echam_phy_tc
    CALL print_echam_phy_config


    ! Parameterizations (with time control)
    ! -------------------------------------

    ! radiation
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) THEN
      !
      CALL  eval_echam_rad_config
      CALL print_echam_rad_config
      !
      ! Radiation constants for gas and cloud optics
      CALL psrad_basic_setup(.false., nlev, pressure_scale, droplet_scale,               &
        &                    echam_cld_config(1)%cinhoml1 ,echam_cld_config(1)%cinhoml2, &
        &                    echam_cld_config(1)%cinhoml3 ,echam_cld_config(1)%cinhomi)
      !
      ! If there are concurrent psrad processes, set up communication 
      ! between the atmo and psrad processes
      CALL setup_atmo_2_psrad()
      !
    END IF

    ! vertical turbulent mixing and surface
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) THEN
      !
      CALL  eval_echam_vdf_config
      CALL print_echam_vdf_config
      !
      ! Allocate memory for the tri-diagonal solver needed by the implicit
      ! time stepping scheme; Compute time-independent parameters.
      !
      CALL init_vdiff_params( nlev, nlev+1, nlev+1, vct )
      !
      ! Currently the tracer indices are sorted such that we count
      ! the water substances first, and then other species like
      ! aerosols and their precursors. "ntracer" is the total number
      ! of tracers (including water substances) handled in the model;
      ! "iqt" is the starting index for non-water species.
      ! Before more sophisticated meta-data structure becomes available,
      ! it is assumed here that all tracers are subject to turbulent mixing.
      !
      khydromet = iqt - 2        ! # of hydrometeors
      ktrac = ntracer - iqt + 1  ! # of non-water species
      !
      CALL init_vdiff_solver( khydromet, ktrac, nlev )
      !
      ! JSBACH land processes
      !

#ifndef __NO_JSBACH__
      IF (ilnd <= nsfc_type .AND. ANY(echam_phy_config(:)%ljsb)) THEN
        DO jg=1,n_dom
          IF (echam_phy_config(jg)%ljsb) THEN 
            CALL jsbach_init(jg)
          END IF
        END DO ! jg
      END IF ! 
#endif

    ENDIF

    ! carbon cycle
    !
    CALL print_ccycle_config

    ! cumulus convection
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_cnv > dt_zero)
    END DO
    IF (lany) THEN
      CALL alloc_echam_cnv_config 
      CALL  eval_echam_cnv_config
      CALL print_echam_cnv_config
    END IF

    ! cloud processes
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_cld > dt_zero)
    END DO
    IF (lany) THEN
      CALL print_echam_cld_config
    END IF

    ! atmospheric gravity wave drag
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_gwd > dt_zero)
    END DO
    IF (lany) THEN
       CALL  eval_echam_gwd_config
       CALL print_echam_gwd_config
    END IF

    ! subgrid scale orographic effects
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_sso > dt_zero)
    END DO
    IF (lany) THEN
       CALL  eval_echam_sso_config
       CALL print_echam_sso_config
    END IF
 
    ! Cariolle linearized o3 chemistry
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_car > dt_zero)
    END DO
    IF (lany) THEN
      IF (io3 < iqt .OR. ntracer < io3) THEN
        CALL finish('init_echam_phy: mo_echam_phy_init.f90', &
                   &'cannot find an ozone tracer - abort')
      END IF
      IF (n_dom > 1) THEN
        CALL finish('init_echam_phy: mo_echam_phy_init.f90', &
                   &'Cariolle initialization not ready for n_dom>1')
      END IF
      CALL lcariolle_init(                                     &
         & openInputFile_wrap,       closeFile_wrap,           &
         & read_bcast_real_3d_wrap,  read_bcast_real_1d_wrap,  &
         & get_constants                                       )
    END IF

    ! ch4 oxidation and h2o photolysis
    !
    lany=.FALSE.
    DO jg = 1,n_dom
      lany = lany .OR. (echam_phy_tc(jg)%dt_mox > dt_zero)
    END DO
    IF (lany) THEN
      CALL init_methox
    END IF

    IF (timers_level > 1) CALL timer_stop(timer_prep_echam_phy)

  END SUBROUTINE init_echam_phy_params


  SUBROUTINE init_echam_phy_external( p_patch, mtime_current)

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(:)
    TYPE(datetime),  INTENT(in), POINTER    :: mtime_current !< Date and time information

    INTEGER :: jg
    LOGICAL :: lany
    TYPE(t_stream_id) :: stream_id

    CHARACTER(len=max_char_length) :: land_frac_fn
    CHARACTER(len=max_char_length) :: land_phys_fn
    CHARACTER(len=max_char_length) :: land_sso_fn

    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    ! Shortcuts to components of echam_rad_config
    !
    INTEGER, POINTER :: ighg(:), irad_aero(:)
    !
    ighg      => echam_rad_config(1:n_dom)% ighg
    irad_aero => echam_rad_config(1:n_dom)% irad_aero
    
    IF (timers_level > 1) CALL timer_start(timer_prep_echam_phy)

    ! external data on ICON grids:
    
    DO jg= 1,n_dom

      IF (n_dom > 1) THEN
        WRITE(land_frac_fn, '(a,i2.2,a)') 'bc_land_frac_DOM', jg, '.nc'
        WRITE(land_phys_fn, '(a,i2.2,a)') 'bc_land_phys_DOM', jg, '.nc'
        WRITE(land_sso_fn , '(a,i2.2,a)') 'bc_land_sso_DOM' , jg, '.nc'
      ELSE
        land_frac_fn = 'bc_land_frac.nc'
        land_phys_fn = 'bc_land_phys.nc'
        land_sso_fn  = 'bc_land_sso.nc'
      ENDIF

      IF (ilnd <= nsfc_type) THEN

        ! land, glacier and lake masks
        !
        WRITE(message_text,'(2a)') 'Read notsea, glac and lake from file ', TRIM(land_frac_fn)
        CALL message('mo_echam_phy_init:init_echam_phy_external', message_text)
        !
        stream_id = openInputFile(land_frac_fn, p_patch(jg), default_read_method)
        CALL read_2D(stream_id=stream_id, location=on_cells,&
             &          variable_name='notsea',               &
             &          fill_array=prm_field(jg)%lsmask(:,:))
        CALL read_2D(stream_id=stream_id, location=on_cells, &
             &          variable_name='glac',               &
             &          fill_array=prm_field(jg)% glac(:,:))
        IF (echam_phy_config(jg)%llake) THEN
          CALL read_2D(stream_id=stream_id, location=on_cells, &
               &          variable_name='lake',               &
               &          fill_array=prm_field(jg)% alake(:,:))
        ELSE
          ! If running without lakes, set lake fraction to zero.
          prm_field(jg)%alake(:,:) = 0._wp
        END IF
        CALL closeFile(stream_id)

        ! For security
        prm_field(jg)%lsmask(:,:) = MERGE(1._wp, prm_field(jg)%lsmask(:,:), &
                                          prm_field(jg)%lsmask(:,:) > 1._wp - 10._wp*EPSILON(1._wp))
        !
        ! At this point, %lsmask is the fraction of land (incl. glacier and
        ! lakes) in the grid box.
        !
        IF (echam_phy_config(jg)%llake) THEN
          !
          ! Substract lake fraction from %lsmask at inner land points (no ocean fraction, %lsmask==1).
          ! Elsewhere (fractional coastal points and ocean points), set alake to zero.
          WHERE (prm_field(jg)%lsmask(:,:) >= 1._wp) ! Inner land point
            prm_field(jg)%lsmask(:,:) = prm_field(jg)%lsmask(:,:) - prm_field(jg)%alake(:,:)
          ELSE WHERE
            prm_field(jg)%alake(:,:) = 0._wp
          END WHERE
        END IF

        ! roughness length and background albedo
        !
        IF (echam_phy_tc(jg)%dt_vdf > dt_zero .OR. echam_phy_tc(jg)%dt_rad > dt_zero) THEN
          !
          stream_id = openInputFile(land_phys_fn, p_patch(jg), default_read_method)
          !
          IF (echam_phy_tc(jg)%dt_vdf > dt_zero) THEN
            !
            WRITE(message_text,'(2a)') 'Read roughness_length from file: ', TRIM(land_phys_fn)
            CALL message('mo_echam_phy_init:init_echam_phy_external', message_text)
            !
            CALL read_2D(stream_id=stream_id, location=on_cells, &
                  &       variable_name='roughness_length',      &
                  &       fill_array=prm_field(jg)% z0m(:,:))
            !
          END IF
          !
          IF (echam_phy_tc(jg)%dt_rad > dt_zero) THEN
            !
            WRITE(message_text,'(2a)') 'Read albedo           from file: ', TRIM(land_phys_fn)
            CALL message('mo_echam_phy_init:init_echam_phy_external', message_text)
            !
            CALL read_2D(stream_id=stream_id, location=on_cells, &
                 &       variable_name='albedo',                &
                 &       fill_array=prm_field(jg)% alb(:,:))
            !
          END IF
          !
          CALL closeFile(stream_id)
          !
        END IF

        ! orography
        IF (echam_phy_tc(jg)%dt_sso > dt_zero) THEN
          !
          WRITE(message_text,'(2a)') 'Read oroxyz from file: ', TRIM(land_sso_fn)
          CALL message('mo_echam_phy_init:init_echam_phy_external', message_text)
          !
          stream_id = openInputFile(land_sso_fn, p_patch(jg), default_read_method)
          CALL read_2D(stream_id=stream_id, location=on_cells, &
               &       variable_name='oromea',                &
               &       fill_array=prm_field(jg)% oromea(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orostd',                &
             &         fill_array=prm_field(jg)% orostd(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orosig',                &
             &         fill_array=prm_field(jg)% orosig(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orogam',                &
             &         fill_array=prm_field(jg)% orogam(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orothe',                &
             &         fill_array=prm_field(jg)% orothe(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='oropic',                &
             &         fill_array=prm_field(jg)% oropic(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='oroval',                &
             &         fill_array=prm_field(jg)% oroval(:,:))
          CALL closeFile(stream_id)
        END IF

      END IF ! (ilnd <= nsfc_type)

    END DO ! jg

    ! external data:


    ! for radiation
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) THEN
      !
      ! parameterized simple plumes of tropospheric aerosols
      !
      IF (ANY(irad_aero(:) == 18)) THEN
        CALL setup_bc_aeropt_splumes
      END IF
      !
      ! well mixed greenhouse gases, horizontally constant
      !
      IF (ANY(ighg(:) > 0)) THEN
        ! read annual means
        IF (.NOT. bc_greenhouse_gases_file_read) THEN
          CALL read_bc_greenhouse_gases
        END IF
        ! interpolate to the current date and time, placing the annual means at
        ! the mid points of the current and preceding or following year, if the
        ! current date is in the 1st or 2nd half of the year, respectively.
        CALL bc_greenhouse_gases_time_interpolation(mtime_current)
        !
      ENDIF
      !
    END IF

    
    ! for radiation and vertical diffusion
    !
    ! Read sea surface temperature, sea ice concentration and depth
    ! Note: For coupled runs, this is only used for initialization of surface temperatures
    !
    IF (iwtr <= nsfc_type .AND. iice <= nsfc_type) THEN
      !
      IF (.NOT. isrestart()) THEN
        !
        ! interpolation weights for linear interpolation of monthly means to the current time
        current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_current)
        !
        DO jg= 1,n_dom
          !
          IF (echam_phy_tc(jg)%dt_rad > dt_zero .OR. echam_phy_tc(jg)%dt_vdf > dt_zero) THEN
            !
            CALL read_bc_sst_sic(mtime_current%date%year, p_patch(jg))
            !
            CALL bc_sst_sic_time_interpolation(current_time_interpolation_weights                   ,&
                 &                             prm_field(jg)%ts_tile(:,:,iwtr)                      ,&
                 &                             prm_field(jg)%seaice (:,:)                           ,&
                 &                             prm_field(jg)%siced  (:,:)                           ,&
                 &                             p_patch(jg)                                          ,&
                 &                             prm_field(jg)%lsmask(:,:) + prm_field(jg)%alake(:,:) < 1._wp ,&
                 &                             .TRUE. )
            !
          END IF
          !
        END DO
        !
      END IF
      !
    END IF

    IF (timers_level > 1) CALL timer_stop(timer_prep_echam_phy)

  END SUBROUTINE init_echam_phy_external


  !-------------
  !>
  !! Loop over all grid levels and give proper values to some components
  !! of the state vectors "prm_field" and "prm_tend".
  !! This subroutine plays a role similar to "init_g3" in ECHAM6.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !!
  SUBROUTINE init_echam_phy_field( p_patch        ,&
    &                              topography_c   ,&
    &                              z_ifc          ,&
    &                              z_mc           ,&
    &                              ddqz_z_full    ,&
    &                              geopot_agl_ifc ,&
    &                              geopot_agl     ,&
    &                              temp           )

    TYPE(t_patch)    ,INTENT(in) :: p_patch
    REAL(wp)         ,INTENT(in) :: topography_c  (:,  :)
    REAL(wp)         ,INTENT(in) :: z_ifc         (:,:,:)
    REAL(wp)         ,INTENT(in) :: z_mc          (:,:,:)
    REAL(wp)         ,INTENT(in) :: ddqz_z_full   (:,:,:)
    REAL(wp)         ,INTENT(in) :: geopot_agl_ifc(:,:,:)
    REAL(wp)         ,INTENT(in) :: geopot_agl    (:,:,:)
    REAL(wp)         ,INTENT(in) :: temp          (:,:,:)

    ! local variables and pointers

    INTEGER  :: jg, ncd, rls, rle, jb, jbs, jbe, jc, jcs, jce
    REAL(wp) :: zlat

    TYPE(t_echam_phy_field),POINTER :: field => NULL()
    TYPE(t_echam_phy_tend) ,POINTER :: tend  => NULL()

      jg = p_patch%id
    
      field => prm_field(jg)
      tend  => prm_tend (jg)

      ! Inquire current grid level and the total number of grid cells
      ncd = MAX(1,p_patch%n_childdom)
      rls = grf_bdywidth_c+1
      rle = min_rlcell_int
      
      jbs     = p_patch%cells%start_blk(rls,  1)
      jbe     = p_patch%cells%  end_blk(rle,ncd)

      ! Assign initial values for some components of the "field" and
      ! "tend" state vectors.

!$OMP PARALLEL WORKSHARE
      !
      ! constant-in-time fields
      ! initial and re-start
      !
      field%      clon(:,  :) = p_patch% cells% center(:,:)% lon
      field%      clat(:,  :) = p_patch% cells% center(:,:)% lat
      field% areacella(:,  :) = p_patch% cells%   area(:,:)
      field%    coriol(:,  :) = p_patch% cells%    f_c(:,:)
      !
      field%      orog(:,  :) =   topography_c(:,  :)
      field%        zh(:,:,:) =          z_ifc(:,:,:)
      field%        zf(:,:,:) =           z_mc(:,:,:)
      field%        dz(:,:,:) =    ddqz_z_full(:,:,:)
      !
      field%      geoi(:,:,:) = geopot_agl_ifc(:,:,:)
      field%      geom(:,:,:) =     geopot_agl(:,:,:)
 
!$OMP END PARALLEL WORKSHARE

      ! in case of restart, reset output fields of unused parameterizations,
      ! to their intial value
      !
      IF (isrestart()) THEN
         !
         IF ( echam_phy_tc(jg)%dt_rad == dt_zero ) THEN
            field% rld_rt      (:,:,:) = 0.0_wp
            field% rlu_rt      (:,:,:) = 0.0_wp
         END IF
         !
         IF ( echam_phy_tc(jg)%dt_cnv == dt_zero ) THEN
            field% rsfc  (:,:) = 0.0_wp
            field% ssfc  (:,:) = 0.0_wp
            field% rtype (:,:) = 0.0_wp
         END IF
         !
         IF ( echam_phy_tc(jg)%dt_cld == dt_zero ) THEN
            field% rsfl (:,:) = 0.0_wp
            field% ssfl (:,:) = 0.0_wp
         END IF
         !
      END IF

      ! vertical diffusion
      IF (.NOT. isrestart()) THEN
        IF (echam_phy_tc(jg)%dt_vdf > dt_zero) THEN
          IF (iwtr<=nsfc_type) field% z0m_tile(:,:,iwtr) = 1e-3_wp !see init_surf in echam (or z0m_oce?)
          IF (iice<=nsfc_type) field% z0m_tile(:,:,iice) = 1e-3_wp !see init_surf in echam (or z0m_ice?)
          IF (ilnd<=nsfc_type) THEN
            field% z0m_tile(:,:,ilnd) = field%z0m(:,:) ! or maybe a larger value?
            field% z0h_lnd(:,:)       = field%z0m(:,:) ! or maybe a larger value?
          END IF
        END IF
      END IF

      ! surface properties
      !
      ! Initialize some variables for water, ice and land tiles
      ! This can be overridden by the testcases below
      ! initial and re-start
      IF (iwtr <= nsfc_type) THEN
        !
        field% albvisdir_tile(:,:,iwtr) = albedoW ! albedo in the visible range for direct radiation
        field% albnirdir_tile(:,:,iwtr) = albedoW ! albedo in the NIR range for direct radiation
        field% albvisdif_tile(:,:,iwtr) = albedoW ! albedo in the visible range for diffuse radiation
        field% albnirdif_tile(:,:,iwtr) = albedoW ! albedo in the NIR range for diffuse radiation
        field% albedo_tile   (:,:,iwtr) = albedoW
        !
      END IF

      IF (ilnd <= nsfc_type) THEN
        !
        IF (.NOT. isrestart()) THEN
          IF (.NOT. ltestcase) THEN
            field%ts_tile(:,:,ilnd) = field%ts_tile(:,:,iwtr)
          END IF
        END IF
        !
        ! initial and re-start
        field% albvisdir_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the visible range for direct radiation
        field% albnirdir_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the NIR range for direct radiation
        field% albvisdif_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the visible range for diffuse radiation
        field% albnirdif_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the NIR range for diffuse radiation
        field% albedo_tile   (:,:,ilnd) = field%alb(:,:)
        !
      END IF

      IF (iice <= nsfc_type) THEN
        !
        IF (.NOT. isrestart()) THEN
          field%ts_tile(:,:,iice) = field%ts_tile(:,:,iwtr)
        END IF
        !
        ! initial and re-start
        field% albvisdir_tile(:,:,iice) = albi    ! albedo in the visible range for direct radiation
        field% albnirdir_tile(:,:,iice) = albi    ! albedo in the NIR range for direct radiation
        field% albvisdif_tile(:,:,iice) = albi    ! albedo in the visible range for diffuse radiation
        field% albnirdif_tile(:,:,iice) = albi    ! albedo in the NIR range for diffuse radiation
        field% albedo_tile   (:,:,iice) = albi
        !
        IF (.NOT. isrestart()) THEN
          ! The ice model should be able to handle different thickness classes,
          ! but for AMIP we ONLY USE one ice class.
          field% albvisdir_ice(:,:,:) = albi ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi ! albedo in the NIR range for direct radiation
          field% albvisdif_ice(:,:,:) = albi ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi ! albedo in the NIR range for diffuse radiation
          field% Tsurf(:,:,:) = Tf
          field% T1   (:,:,:) = Tf
          field% T2   (:,:,:) = Tf
          WHERE (field%seaice(:,:) > 0.0_wp)
             field% hs   (:,1,:) = 0.1_wp       ! set initial snow depth on sea ice
          ELSEWHERE
             field% hs   (:,1,:) = 0.0_wp
          ENDWHERE
          field% hi   (:,1,:) = field%siced(:,:)
          field% conc (:,1,:) = field%seaice(:,:)
        END IF
        !
      END IF

      ! For idealized test cases

      SELECT CASE (nh_test_name)
      CASE('APE','APE_echam','RCEhydro','RCE_glb') !Note that there is only one surface type in this case
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              field% ts_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
          !
        END DO
!$OMP END PARALLEL DO

      CASE('RCE') !Note that there is only one surface type in this case
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              field% ts_tile(jc,jb,iwtr) = th_cbl(1)
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
          !
        END DO
!$OMP END PARALLEL DO

      CASE('APEi')
        ! The same as APE, except that whenever SST reaches tmelt, we put
        ! 1m-thick ice with a concentration of 0.9 on top
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              ! SST must reach Tf where there's ice. It may be better to modify ape_sst it self.
              field% ts_tile (jc,jb,iwtr) = ape_sst(ape_sst_case,zlat) + Tf
              field% ts_tile (jc,jb,iice) = Tf + tmelt ! K
              field% Tsurf   (jc,1, jb  ) = Tf         ! degC
              field% T1      (jc,1, jb  ) = Tf         ! degC
              field% T2      (jc,1, jb  ) = Tf         ! degC
              field% hs      (jc,1, jb  ) = 0._wp
              IF ( field% ts_tile(jc,jb,iwtr) <= Tf + tmelt ) THEN
                field% Tsurf (jc,1,jb) = field% ts_tile(jc,jb,iice) - tmelt
                field% conc  (jc,1,jb) = 0.9_wp
                field% hi    (jc,1,jb) = 1.0_wp
                field% seaice(jc,  jb) = field%conc(jc,1,jb)
              ELSE
                field% conc  (jc,1,jb) = 0._wp
                field% hi    (jc,1,jb) = 0._wp
                field% seaice(jc,  jb) = field%conc(jc,1,jb)
              ENDIF
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          !
        END DO
!$OMP END PARALLEL DO
        IF (.NOT. isrestart()) THEN
          field% albvisdir_ice(:,:,:) = albi    ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi    ! albedo in the NIR range for direct radiation
          field% albvisdif_ice(:,:,:) = albi    ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi    ! albedo in the NIR range for diffuse radiation
        END IF
       
      CASE('APEc','APEc_nh')
        ! The same as APEi, except we initialize with no ice and don't modify the surface
        ! temperature. This is meant for a coupled run.
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              field% ts_tile (jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
              field% ts_tile (jc,jb,iice) = Tf + tmelt ! K
              field% Tsurf   (jc,1,jb)    = Tf         ! degC
              field% T1      (jc,1,jb)    = Tf         ! degC
              field% T2      (jc,1,jb)    = Tf         ! degC
              field% hs      (jc,1,jb)    = 0._wp
              field% conc    (jc,1,jb)    = 0._wp
              field% hi      (jc,1,jb)    = 0._wp
              field% seaice  (jc,  jb)    = 0._wp
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          !
        END DO
!$OMP END PARALLEL DO
        IF (.NOT. isrestart()) THEN
          field% albvisdir_ice(:,:,:) = albi  ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi  ! albedo in the NIR range for direct radiation
          field% albvisdif_ice(:,:,:) = albi  ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi  ! albedo in the NIR range for diffuse radiation
        END IF
     
      CASE('TPEc', 'TPEo') !Note that there is only one surface type (ilnd) in this case
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 1._wp   ! land fraction = 1
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          !
          IF (.NOT. isrestart()) THEN
            field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
            field% ts_tile(jcs:jce,jb,ilnd) = tpe_temp
          END IF
          !
        END DO
!$OMP END PARALLEL DO

      CASE('JWw-Moist','LDF-Moist','jabw_m')
        !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          ! Set the surface temperature to the same value as the lowest model
          ! level above surface. For this test case, currently we assume
          ! there is no land or sea ice.
          !
          IF (.NOT. isrestart()) THEN
            field% ts_tile(jcs:jce,jb,iwtr) = temp(jcs:jce, p_patch%nlev,jb)
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
        END DO
!$OMP END DO  NOWAIT
!$OMP END PARALLEL

      END SELECT

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      ! initial and re-start
      DO jb = jbs,jbe
        !
        CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
        IF (jcs>jce) CYCLE 
        !
        ! Set surface tiling fractions, wrt. the cell area
        DO jc = jcs,jce
          field% sftlf (jc,jb) = field% lsmask(jc,jb) + field% alake(jc,jb) ! land incl. lakes
          field% sftgif(jc,jb) = field% lsmask(jc,jb) * field% glac (jc,jb) ! land ice
          field% sftof (jc,jb) = 1._wp - field% sftlf(jc,jb)                ! ocean
        END DO
        !
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Settings for total surface
      ! (after tile masks and variables potentially have been overwritten by testcases above)

      ! initial and re-start
      IF (iwtr <= nsfc_type) THEN
        field%albedo   (:,:) = albedoW
      ELSE
        field%albedo   (:,:) = field%alb(:,:)
      END IF

      IF (.NOT. isrestart()) THEN
        !
        IF (iwtr <= nsfc_type) THEN
          field%ts       (:,:) = field%ts_tile(:,:,iwtr)
          field%albvisdir(:,:) = albedoW
          field%albvisdif(:,:) = albedoW
          field%albnirdir(:,:) = albedoW
          field%albnirdif(:,:) = albedoW
        ELSE
          field%ts       (:,:) = field%ts_tile(:,:,ilnd)
          field%albvisdir(:,:) = field%alb(:,:)
          field%albvisdif(:,:) = field%alb(:,:)
          field%albnirdir(:,:) = field%alb(:,:)
          field%albnirdif(:,:) = field%alb(:,:)
        END IF
        !
        field%ts_rad     (:,:) = field%ts(:,:)
        field%ts_rad_rt  (:,:) = field%ts(:,:)
        !
      END IF

      NULLIFY( field,tend )

  END SUBROUTINE init_echam_phy_field


  !-------------
  !>
  !! Initialize the O3 tracer from the Cariolle initial ozone field.
  !! Initialize ozone mass mixing ratios for Cariolle scheme. 
  !! An approximative initialization that considers the atmosphere as being dry is enough.
  !!
  !! @par Revision History
  !! Initial version by Sebastian Rast, MPI-M (2016-11)
  !! Changes:
  !! - separate subroutine, Marco Giorgetta, MPI-M (2018-03)
  !!
  SUBROUTINE init_o3_lcariolle( mtime_current  ,&
    &                           p_patch        ,&
    &                           pres           ,&
    &                           o3              )

    TYPE(datetime)   ,POINTER    :: mtime_current
    TYPE(t_patch)    ,INTENT(in) :: p_patch
    REAL(wp)         ,INTENT(in) :: pres          (:,:,:)
    REAL(wp)         ,INTENT(out):: o3            (:,:,:)

    ! local variables
    INTEGER  :: ncd, rls, rle, jb, jbs, jbe, jcs, jce

    ! Variables for Cariolle ozone scheme
    TYPE(t_avi)                        :: avi
    REAL(wp), TARGET                   :: latc  (nproma)
    REAL(wp), TARGET                   :: pfull (nproma, p_patch%nlev)
    REAL(wp)                           :: vmr_o3(nproma, p_patch%nlev)
    TYPE(t_time_interpolation)         :: time_interpolation
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
    !
    ! External routines of the Cariolle library
    EXTERNAL :: lcariolle_init_o3, lcariolle_lat_intp_li, lcariolle_pres_intp_li

    avi%ldown=.TRUE.
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_current)
    time_interpolation%imonth1=current_time_interpolation_weights%month1_index
    time_interpolation%imonth2=current_time_interpolation_weights%month2_index
    time_interpolation%weight1=current_time_interpolation_weights%weight1
    time_interpolation%weight2=current_time_interpolation_weights%weight2

    ! Inquire current grid level and the total number of grid cells
    ncd = MAX(1,p_patch%n_childdom)
    rls = grf_bdywidth_c+1
    rle = min_rlcell_int
    jbs     = p_patch%cells%start_blk(rls,  1)
    jbe     = p_patch%cells%  end_blk(rle,ncd)
    !
    DO jb = jbs,jbe
      !
      CALL get_indices_c(p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
      IF (jcs>jce) CYCLE 
      !
      pfull(:,:)          =  pres (:,:,jb)
      avi%pres            => pfull
      !
      latc(:)             =  p_patch% cells% center(:,jb)% lat
      avi%cell_center_lat => latc
      !
      CALL lcariolle_init_o3(                                               &
        & jcs,                   jce,                nproma,                &
        & p_patch%nlev,          time_interpolation, lcariolle_lat_intp_li, &
        & lcariolle_pres_intp_li,avi,                vmr_o3                 )
      !
      o3(jcs:jce,:,jb) = vmr_o3(jcs:jce,:)*amo3/amd
      !
    END DO
    !
    l_cariolle_initialized_o3 = .TRUE.

  END SUBROUTINE init_o3_lcariolle

END MODULE mo_echam_phy_init
