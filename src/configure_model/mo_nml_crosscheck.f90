!>
!! This module checks the read-in namelist parameters and, in case of
!! inconsistencies, it tries to correct these.
!!
!!
!! @author Kristina Froehlich, MPI-M (2011-07-12)
!! @author Hui Wan, MPI-M (2011-07-12)
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
MODULE mo_nml_crosscheck

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int32_t
  USE mo_kind,               ONLY: wp, i8
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length, max_dom,                         &
    &                              iecham, ildf_echam, inwp, iheldsuarez,            &
    &                              ildf_dry, inoforcing, ihs_atm_temp,               &
    &                              ihs_atm_theta, tracer_only, inh_atmosphere,       &
    &                              ishallow_water, LEAPFROG_EXPL, LEAPFROG_SI,       &
    &                              NO_HADV, UP, MIURA, MIURA3, FFSL, FFSL_HYB,       &
    &                              MCYCL, MIURA_MCYCL, MIURA3_MCYCL,                 &
    &                              FFSL_MCYCL, FFSL_HYB_MCYCL, ifluxl_sm,            &
    &                              ifluxl_m, ihs_ocean, RAYLEIGH_CLASSIC,            &
    &                              iedmf, icosmo, MODE_IAU, MODE_IAU_OLD,            &
    &                              DEFAULT_RESTART_INTVL, DEFAULT_CHECKPT_INTVL
  USE mo_master_config,      ONLY: tc_exp_stopdate, tc_startdate, tc_stopdate,       &
    &                              experimentReferenceDate,                          &
    &                              experimentStartDate,                              &
    &                              experimentStopDate, isRestart,                    &
    &                              setExpRefdate, setExpStartdate,                   &
    &                              setExpStopdate, setStartdate, setStopdate,        &
    &                              master_nml_calendar => calendar,                  &
    &                              checkpointTimeIntval, restartTimeIntval,          &
    &                              setrestarttimeinterval, setcheckpointtimeinterval
  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN, datetime,                   &
    &                              MAX_CALENDAR_STR_LEN,                             &
    &                              MAX_TIMEDELTA_STR_LEN,                            &
    &                              OPERATOR(>),OPERATOR(/=), newDatetime,            &
    &                              deallocateDatetime, timedelta,                    &
    &                              getPTStringFromMS, newTimedelta, min,             &
    &                              deallocateTimedelta, OPERATOR(+), OPERATOR(==),   &
    &                              timedeltatostring, datetimetostring,              &
    &                              OPERATOR(*), OPERATOR(<), OPERATOR(<=),           &
    &                              setcalendar, calendarToString,                    &
    &                              mtime_proleptic_gregorian => proleptic_gregorian, &
    &                              mtime_year_of_365_days => year_of_365_days,       &
    &                              mtime_year_of_360_days => year_of_360_days,       &
    &                              getTotalMilliSecondsTimeDelta
  USE mo_mtime_extensions,   ONLY: get_datetime_string, datetime_str_equal,          &
    &                              get_timedelta_divide_by_seconds,                  &
    &                              getTimedeltaFromMS, timedelta_str_equal
  USE mo_time_config,        ONLY: t_time_config, dt_restart, time_config,           &
    &                              ini_datetime_string, is_relative_time,            &
    &                              time_nml_calendar => calendar,                    &
    &                              restart_calendar, restart_ini_datetime_string
  USE mo_extpar_config,      ONLY: itopo                                             
  USE mo_io_config,          ONLY: dt_checkpoint, lflux_avg,inextra_2d,              &
    &                              inextra_3d
  USE mo_parallel_config,    ONLY: check_parallel_configuration,                     &
    &                              num_io_procs, itype_comm, num_restart_procs
  USE mo_run_config,         ONLY: nsteps, dtime, iforcing,                          &
    &                              ltransport, ntracer, nlev, ltestcase,             &
    &                              nqtendphy, iqtke, iqv, iqc, iqi,                  &
    &                              iqs, iqr, iqt, iqtvar, ico2, ltimer,              &
    &                              iqni, iqni_nuc, iqg, iqm_max,                     &
    &                              iqh, iqnr, iqns, iqng, iqnh, iqnc,                & 
    &                              inccn, ininact, ininpot,                          &
    &                              activate_sync_timers, timers_level,               &
    &                              output_mode, lart, tc_dt_model,                   &
    &                              mtime_modelTimeStep => modelTimeStep,             &
    &                              setModelTimeStep
  USE mo_dynamics_config,    ONLY: iequations, idiv_method,                          &
    &                              divavg_cntrwgt, sw_ref_height,                    &
    &                              lcoriolis, lshallow_water, ltwotime
  USE mo_advection_config,   ONLY: advection_config

  USE mo_nonhydrostatic_config, ONLY: itime_scheme_nh => itime_scheme,               &
                                      lhdiff_rcf, rayleigh_type, divdamp_order
  USE mo_ha_dyn_config,      ONLY: ha_dyn_config
  USE mo_diffusion_config,   ONLY: diffusion_config
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, icpl_aero_conv
  USE mo_lnd_nwp_config,     ONLY: ntiles_lnd, lsnowtile
  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_radiation_config
  USE mo_echam_conv_config,  ONLY: echam_conv_config
  USE mo_gw_hines_config,    ONLY: gw_hines_config
  USE mo_vdiff_config,       ONLY: vdiff_config
  USE mo_turbdiff_config,    ONLY: turbdiff_config
  USE mo_initicon_config,    ONLY: init_mode, dt_iau, ltile_coldstart
  USE mo_nh_testcases_nml,   ONLY: linit_tracer_fv,nh_test_name
  USE mo_ha_testcases,       ONLY: ctest_name, ape_sst_case

  USE mo_datetime,           ONLY: add_time, print_datetime_all, date_to_time,       &
    &                              string_to_datetime,                               &
    &                              dtime_proleptic_gregorian => proleptic_gregorian, &
    &                              dtime_cly360              => cly360,              &
    &                              dtime_julian_gregorian    => julian_gregorian
  USE mo_meteogram_config,   ONLY: check_meteogram_configuration
  USE mo_master_control,     ONLY: get_my_process_type,      &
    & testbed_process,  atmo_process, ocean_process, radiation_process
  USE mo_time_config,        ONLY: end_datetime_string
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_io_restart,              ONLY: read_restart_header
  USE mo_io_restart_attributes,   ONLY: get_restart_attribute
  USE mo_util_string,             ONLY: tolower, int2string
  USE mo_grid_config,        ONLY: grid_rescale_factor, patch_weight,                &
    &                              start_time, lplane, n_dom, init_grid_configuration

  USE mo_art_config,         ONLY: art_config

  USE mo_gridref_config
  USE mo_interpol_config
  USE mo_sleve_config


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_crosscheck
  PUBLIC :: compute_date_settings


  CHARACTER(LEN = *), PARAMETER :: modname = "mo_nml_crosscheck"


CONTAINS

  SUBROUTINE atm_crosscheck

    INTEGER :: jg
    INTEGER :: jt   ! tracer loop index
    INTEGER :: i_listlen
    INTEGER :: z_go_tri(11)  ! for crosscheck
    CHARACTER(len=*), PARAMETER :: method_name =  'mo_nml_crosscheck:atm_crosscheck'

    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings(dt_restart, nsteps, time_config)


    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()


    !--------------------------------------------------------------------
    ! Grid and dynamics
    !--------------------------------------------------------------------

    ! check the configuration
    CALL init_grid_configuration()
    
    IF (lplane) CALL finish( TRIM(method_name),&
      'Currently a plane version is not available')

    SELECT CASE (iequations)
    CASE(IHS_ATM_TEMP,IHS_ATM_THETA)         ! hydrostatic atm model

      SELECT CASE(iforcing)
      CASE(INOFORCING,IHELDSUAREZ,ILDF_DRY)  ! without moist processes
        ha_dyn_config%ldry_dycore = .TRUE.
      END SELECT

    END SELECT

    lshallow_water = (iequations==ISHALLOW_WATER)

    SELECT CASE (iequations)
    CASE (IHS_ATM_TEMP,IHS_ATM_THETA,ISHALLOW_WATER)

      ltwotime = (ha_dyn_config%itime_scheme/=LEAPFROG_EXPL).AND. &
                 (ha_dyn_config%itime_scheme/=LEAPFROG_SI)

    END SELECT

    !--------------------------------------------------------------------
    ! Testcases (hydrostatic)
    !--------------------------------------------------------------------

    IF ((TRIM(ctest_name)=='GW') .AND. (nlev /= 20)) THEN
      CALL finish(TRIM(method_name),'nlev MUST be 20 for the gravity-wave test case')
    ENDIF

    IF ((TRIM(ctest_name)=='SV') .AND. ntracer /= 2 ) THEN
      CALL finish(TRIM(method_name), &
        & 'ntracer MUST be 2 for the stationary vortex test case')
    ENDIF

    IF ((TRIM(ctest_name)=='DF1') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(method_name), &
        & 'ntracer MUST be >=2 for the deformational flow test case 1')
    ENDIF

    IF ((TRIM(ctest_name)=='DF2') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(method_name), &
        & 'ntracer MUST be >=2 for the deformational flow test case 2')
    ENDIF

    IF ((TRIM(ctest_name)=='DF3') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(method_name), &
        & 'ntracer MUST be >=2 for the deformational flow test case 3')
    ENDIF

    IF ((TRIM(ctest_name)=='DF4') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(method_name), &
        & 'ntracer MUST be >=2 for the deformational flow test case 4')
    ENDIF

    IF ((TRIM(ctest_name)=='APE') .AND. (TRIM(ape_sst_case)=='sst_ice')  ) THEN
      IF (.NOT. lflux_avg)&
      CALL finish(TRIM(method_name), &
        & 'lflux_avg must be set true to run this setup')
    ENDIF

    !--------------------------------------------------------------------
    ! Testcases (nonhydrostatic)
    !--------------------------------------------------------------------
    IF (.NOT. ltestcase .AND. rayleigh_type == RAYLEIGH_CLASSIC) THEN
      CALL finish(TRIM(method_name), &
        & 'rayleigh_type = RAYLEIGH_CLASSIC not applicable to real case runs.')
    ENDIF

    IF ( ( TRIM(nh_test_name)=='APE_nwp'.OR. TRIM(nh_test_name)=='dcmip_tc_52' ) .AND.  &
      &  ( ANY(atm_phy_nwp_config(:)%inwp_surface == 1 ) ) .AND.                       &
      &  ( ANY(atm_phy_nwp_config(:)%inwp_turb    /= iedmf ) ) ) THEN
      CALL finish(TRIM(method_name), &
        & 'surface scheme must be switched off, when running the APE test')
    ENDIF

    !--------------------------------------------------------------------
    ! Shallow water
    !--------------------------------------------------------------------
    IF (iequations==ISHALLOW_WATER.AND.ha_dyn_config%lsi_3d) THEN
      CALL message( TRIM(method_name), 'lsi_3d = .TRUE. not applicable to shallow water model')
    ENDIF

    IF ((iequations==ISHALLOW_WATER).AND.(nlev/=1)) &
    CALL finish(TRIM(method_name),'Multiple vertical level specified for shallow water model')

    !--------------------------------------------------------------------
    ! Hydrostatic atm
    !--------------------------------------------------------------------
    IF (iequations==IHS_ATM_THETA) ha_dyn_config%ltheta_dyn = .TRUE.

    !--------------------------------------------------------------------
    ! Nonhydrostatic atm
    !--------------------------------------------------------------------
    IF (lhdiff_rcf .AND. (itype_comm == 3)) CALL finish(TRIM(method_name), &
      'lhdiff_rcf is available only for idiv_method=1 and itype_comm<=2')

    IF (grf_intmethod_e >= 5 .AND. iequations /= INWP .AND. n_dom > 1) THEN
      grf_intmethod_e = 4
      CALL message( TRIM(method_name), 'grf_intmethod_e has been reset to 4')
    ENDIF

    !--------------------------------------------------------------------
    ! Atmospheric physics, general
    !--------------------------------------------------------------------
    IF ((iforcing==INWP).AND.(iequations/=INH_ATMOSPHERE)) &
    CALL finish( TRIM(method_name), 'NWP physics only implemented in the '//&
               'nonhydrostatic atm model')

    !--------------------------------------------------------------------
    ! NWP physics
    !--------------------------------------------------------------------
    IF (iforcing==inwp) THEN

      DO jg =1,n_dom

        IF( atm_phy_nwp_config(jg)%inwp_satad == 0       .AND. &
          & ((atm_phy_nwp_config(jg)%inwp_convection >0 ) .OR. &
          &  (atm_phy_nwp_config(jg)%inwp_gscp > 0      )    ) ) &
        &  CALL finish( TRIM(method_name),'satad has to be switched on')


        IF( (atm_phy_nwp_config(jg)%inwp_gscp==0) .AND. &
          & (atm_phy_nwp_config(jg)%inwp_convection==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_radiation==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_sso==0)  .AND. &
          & (atm_phy_nwp_config(jg)%inwp_surface == 0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_turb> 0) )   &
        CALL message(TRIM(method_name),' WARNING! NWP forcing set but '//&
                    'only turbulence selected!')


        IF (( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) .AND. &
          & (turbdiff_config(jg)%lconst_z0) ) THEN
          CALL message(TRIM(method_name),' WARNING! NWP forcing set but '//  &
                      'idealized (horizontally homogeneous) roughness '//&
                      'length z0 selected!')
        ENDIF

        IF (.NOT. ltestcase .AND. atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
          CALL finish( TRIM(method_name),'Real-data applications require using a surface scheme!')
        ENDIF

        ! check radiation scheme in relation to chosen ozone and irad_aero=6 to itopo

        IF ( (atm_phy_nwp_config(jg)%inwp_radiation > 0).OR.(echam_phy_config%lrad) )  THEN

          SELECT CASE (irad_o3)
          CASE (0) ! ok
            CALL message(TRIM(method_name),'radiation is used without ozone')
          CASE (2,4,6,7,8,9) ! ok
            CALL message(TRIM(method_name),'radiation is used with ozone')
          CASE default
            CALL finish(TRIM(method_name),'irad_o3 currently has to be 0, 2, 4, 6, 7, 8 or 9.')
          END SELECT

          ! Tegen aerosol and itopo (Tegen aerosol data have to be read from external data file)
          IF ( ( irad_aero == 6 ) .AND. ( itopo /=1 ) ) THEN
            CALL finish(TRIM(method_name),'irad_aero=6 requires itopo=1')
          ENDIF

          IF ( irad_aero /= 6 .AND. (atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .OR. icpl_aero_conv > 0)) THEN
            CALL finish(TRIM(method_name),'aerosol-precipitation coupling requires irad_aero=6')
          ENDIF
        ELSE

          SELECT CASE (irad_o3)
          CASE(0) ! ok
          CASE default
            irad_o3 = 0
            CALL message(TRIM(method_name),'running without radiation => irad_o3 reset to 0')
          END SELECT

        ENDIF !inwp_radiation

        !! check microphysics scheme
        IF (  atm_phy_nwp_config(jg)%mu_rain < 0.0   .OR. &
          &   atm_phy_nwp_config(jg)%mu_rain > 5.0)  THEN
          CALL finish(TRIM(method_name),'mu_rain requires: 0 < mu_rain < 5')
        END IF

        IF (  atm_phy_nwp_config(jg)%mu_snow < 0.0   .OR. &
          &   atm_phy_nwp_config(jg)%mu_snow > 5.0)  THEN
          CALL finish(TRIM(method_name),'mu_snow requires: 0 < mu_snow < 5')
        END IF ! microphysics

        IF (atm_phy_nwp_config(jg)%inwp_surface == 0 .AND. ntiles_lnd > 1) THEN
          ntiles_lnd = 1
          CALL message(TRIM(method_name),'Warning: ntiles reset to 1 because the surface scheme is turned off')
        ENDIF

      ENDDO
    END IF


    !--------------------------------------------------------------------
    ! Tracers and diabatic forcing
    !--------------------------------------------------------------------

    !
    ! Check settings of ntracer
    !
    ! Set tracer indices
    !
    SELECT CASE(iforcing)
    CASE (IECHAM,ILDF_ECHAM)  ! iforcing

      IF (ntracer < 3) CALL finish(TRIM(method_name),'ECHAM physics needs at least 3 tracers')

      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqr    = 0     !! 0: no rain water
      iqs    = 0     !! 0: no snow
      ico2   = 4     !! CO2
      iqm_max= 3     !! end index of water species mixing ratios
      iqt    = 4     !! starting index of non-water species
      nqtendphy = 0  !! number of water species for which convective and turbulent
                     !! tendencies are stored

    CASE (INWP) ! iforcing

      ! ** NWP physics section ** 
      !
      ! IMPORTANT: For NWP physics, five microphysics tracers (QV, QC, QI, QR and QS) must always be
      !            defined because their presence is implicitly assumed in the physics-dynamics interface
      !            when converting between temperature and virtual temperature.
      !
      !            Any additional mass-related microphysics tracers must be numbered in consecutive order
      !            after iqs = 5. The parameter "iqm_max" must signify the highest tracer index carrying a
      !            moisture mixing ratio. Additional tracers for hydrometeor number concentrations (in
      !            case of a two-moment scheme) or other purposes (aerosols, qt variance or anything else)
      !            can be numbered freely after iqm_max. The parameter "iqt", denoting the start index
      !            of tracers not related at all to moisture, is used in configure_advection to specify
      !            the index range of tracers for which advection is turned off in the stratosphere
      !            (i.e. all cloud and precipitation variables including number concentrations)
      !
      !            Note also that the namelist parameter "ntracer" is reset automatically to the correct
      !            value when NWP physics is used in order to avoid multiple namelist changes when playing
      !            around with different physics schemes. 
      !
      ! Default settings valid for all microphysics options
      !
      iqv       = 1     !> water vapour
      iqc       = 2     !! cloud water
      iqi       = 3     !! ice 
      iqr       = 4     !! rain water
      iqs       = 5     !! snow
      nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
      !
      ! The following parameters may be reset depending on the selected physics scheme
      !
      iqm_max   = 5     !! end index of water species mixing ratios
      iqt       = 6     !! start index of other tracers not related at all to moisture
      !
      ntracer   = 5     !! total number of tracers
      !
      ! dummy settings
      iqni     = ntracer+100    !! cloud ice number
      iqni_nuc = ntracer+100    !! activated ice nuclei  
      iqg      = ntracer+100    !! graupel
      iqtvar   = ntracer+100    !! qt variance (for EDMF turbulence)
      iqh      = ntracer+100
      iqnr     = ntracer+100  
      iqns     = ntracer+100
      iqng     = ntracer+100
      iqnh     = ntracer+100
      iqnc     = ntracer+100
      inccn    = ntracer+100
      ininpot  = ntracer+100
      ininact  = ntracer+100

      !
      !
      SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
        

      CASE(2)  ! COSMO-DE (3-cat ice: snow, cloud ice, graupel)

       ! CALL finish('mo_atm_nml_crosscheck', 'Graupel scheme not implemented.')
        
        iqg     = 6       !! graupel
        iqm_max = iqg
        iqt     = iqt + 1

        ntracer = ntracer + 1  !! increase total number of tracers by 1

 
      CASE(3)  ! improved ice nucleation scheme C. Koehler (note: iqm_max does not change!)

        iqni     = 6     !! cloud ice number
        iqni_nuc = 7     !! activated ice nuclei  
        iqt     = iqt + 2

        ntracer = ntracer + 2  !! increase total number of tracers by 2

      CASE(4)  ! two-moment scheme 
      
        iqg  = 6
        iqh  = 7
        iqni = 8        
        iqnr = 9        
        iqns = 10        
        iqng = 11        
        iqnh = 12
        iqnc = 13
        ininact = 14

        nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
        iqm_max   = 7     !! end index of water species mixing ratios
        iqt       = 15    !! start index of other tracers not related at all to moisture
       
        ntracer = 14

      CASE(5)  ! two-moment scheme with CCN and IN budgets
      
        iqg  = 6
        iqh  = 7
        iqni = 8        
        iqnr = 9        
        iqns = 10        
        iqng = 11        
        iqnh = 12
        iqnc = 13
        ininact = 14
        inccn   = 15
        ininpot = 16

        nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
        iqm_max   = 7     !! end index of water species mixing ratios
        iqt       = 17    !! start index of other tracers not related at all to moisture
       
        ntracer = 16
        
      CASE(6)
      
        iqg  = 6
        iqh  = 7
        iqni = 8        
        iqnr = 9        
        iqns = 10        
        iqng = 11        
        iqnh = 12
        iqnc = 13
        
        nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
        iqm_max   = 7     !! end index of water species mixing ratios
        iqt       = 14    !! start index of other tracers not related at all to moisture
        
        ntracer = 13
        
      END SELECT ! microphysics schemes


      IF (atm_phy_nwp_config(jg)%inwp_turb == iedmf) THEN ! EDMF turbulence

        iqtvar = iqt       !! qt variance
        iqt    = iqt + 1   !! start index of other tracers than hydrometeors

        ntracer = ntracer + 1  !! increase total number of tracers by 1

        ntiles_lnd = 5     !! EDMF currently only works with 5 land tiles - consistent with TESSEL
                           !! even if the land model is inactive ntiles_lnd should be 5

      ENDIF

      IF ( (advection_config(jg)%iadv_tke) > 0 ) THEN
        IF ( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) THEN
          iqtke = iqt        !! TKE
 
          ! Note that iqt is not increased, since TKE does not belong to the hydrometeor group.

          ntracer = ntracer + 1  !! increase total number of tracers by 1

          WRITE(message_text,'(a,i3)') 'Attention: TKE is advected, '//&
                                       'ntracer is increased by 1 to ',ntracer
          CALL message(TRIM(method_name),message_text)
        ELSE
          WRITE(message_text,'(a,i2)') 'TKE advection not supported for inwp_turb= ', &
            &                          atm_phy_nwp_config(jg)%inwp_turb
          CALL finish(TRIM(method_name), TRIM(message_text) )
        ENDIF
      ENDIF

      ! Note: Indices for additional tracers are assigned automatically
      ! via add_tracer_ref in mo_nonhydro_state.

      WRITE(message_text,'(a,i3)') 'Attention: NWP physics is used, '//&
                                   'ntracer is automatically reset to ',ntracer
      CALL message(TRIM(method_name),message_text)

      ! take into account additional passive tracers, if present
      ! iqt is not increased, since passive tracers do not belong to the hydrometeor group.
      IF ( advection_config(jg)%npassive_tracer > 0) THEN
        ntracer = ntracer + advection_config(jg)%npassive_tracer
        WRITE(message_text,'(a,i3,a,i3)') 'Attention: passive tracers have been added, '//&
                                     'ntracer is increased by ',advection_config(jg)%npassive_tracer, &
                                     ' to ',ntracer
        CALL message(TRIM(method_name),message_text)
      ENDIF


      IF (lart) THEN
        
        ntracer = ntracer + art_config(jg)%iart_ntracer
        
        WRITE(message_text,'(a,i3,a,i3)') 'Attention: transport of ART tracers is active, '//&
                                     'ntracer is increased by ',art_config(jg)%iart_ntracer, &
                                     ' to ',ntracer
        CALL message(TRIM(method_name),message_text)

      ENDIF

      ! set the nclass_gscp variable for land-surface scheme to number of hydrometeor mixing ratios
      DO jg = 1, n_dom
        atm_phy_nwp_config(jg)%nclass_gscp = iqm_max
      ENDDO

    CASE default ! iforcing

        iqv    = 1     !> water vapour
        iqc    = 2     !! cloud water
        iqi    = 3     !! ice
        iqr    = 0     !! 0: no rain water
        iqs    = 0     !! 0: no snow
        ico2   = 5     !! CO2
        iqm_max= 3     !! end index of water species mixing ratios
        iqt    = 4     !! starting index of non-water species
        nqtendphy = 0  !! number of water species for which convective and turbulent
                       !! tendencies are stored

    END SELECT ! iforcing


    IF (ltransport) THEN
    DO jg = 1,n_dom

      i_listlen = LEN_TRIM(advection_config(jg)%ctracer_list)

      SELECT CASE ( iforcing )
      CASE ( INWP )
      !...........................................................
      ! in NWP physics
      !...........................................................


        ! Force settings for tracer iqtke, if TKE advection is performed
        !
        IF ( advection_config(jg)%iadv_tke > 0 ) THEN

          ! force monotonous slope limiter for vertical advection
          advection_config(jg)%itype_vlimit(iqtke) = 2

          ! force positive definite flux limiter for horizontal advection
          advection_config(jg)%itype_hlimit(iqtke) = 4

          SELECT CASE (advection_config(jg)%iadv_tke)
          CASE (1)
            ! switch off horizontal advection
            advection_config(jg)%ihadv_tracer(iqtke) = 0
          CASE (2)
            ! check whether horizontal substepping is switched on 
            IF (ALL( (/22,32,42,52/) /= advection_config(jg)%ihadv_tracer(iqtke)) ) THEN
              ! choose Miura with substepping
              advection_config(jg)%ihadv_tracer(iqtke) = 22
            ENDIF
          END SELECT

        ENDIF


        IF ( i_listlen /= ntracer ) THEN
          DO jt=1,ntracer
            WRITE(advection_config(jg)%ctracer_list(jt:jt),'(i1.1)')jt
          ENDDO
          WRITE(message_text,'(a,a)') &
            & 'Attention: according to physics, ctracer_list is set to ',&
            & advection_config(jg)%ctracer_list(1:ntracer)
          CALL message(TRIM(method_name),message_text)
        ENDIF


      CASE (inoforcing, iheldsuarez, iecham, ildf_dry, ildf_echam)
      !...........................................................
      ! Other types of adiabatic forcing
      !...........................................................

        IF ( i_listlen < ntracer .AND. i_listlen /= 0 ) THEN
          ntracer = i_listlen
          CALL message(TRIM(method_name),'number of tracers is adjusted according to given list')
        END IF


        IF (echam_phy_config%lrad) THEN
          IF ( izenith > 5)  &
            CALL finish(TRIM(method_name), 'Choose a valid case for rad_nml: izenith.')
        ENDIF
      END SELECT ! iforcing

    END DO ! jg = 1,n_dom
    END IF ! ltransport


    !--------------------------------------------------------------------
    ! Tracer transport
    !--------------------------------------------------------------------
    ! General

    SELECT CASE (iequations)
    CASE (INH_ATMOSPHERE)

      IF ((itime_scheme_nh==tracer_only) .AND. (.NOT.ltransport)) THEN
        WRITE(message_text,'(A,i2,A)') &
          'nonhydrostatic_nml:itime_scheme set to ', tracer_only, &
          '(TRACER_ONLY), but ltransport to .FALSE.'
        CALL finish( TRIM(method_name),TRIM(message_text))
      END IF

    CASE (IHS_ATM_TEMP,IHS_ATM_THETA,ISHALLOW_WATER)

      IF ( (ha_dyn_config%itime_scheme==tracer_only).AND. &
           (.NOT.ltransport)) THEN
        WRITE(message_text,'(A,i2,A)') &
          'ha_dyn_nml:itime_scheme set to ', tracer_only, &
          '(TRACER_ONLY), but ltransport to .FALSE.'
        CALL finish( TRIM(method_name),TRIM(message_text))
      END IF

    END SELECT

#ifndef __ICON_ART
    IF (lart) THEN
      WRITE(message_text,'(A)') &
          'run_nml: lart is set .TRUE. but ICON was compiled without -D__ICON_ART'
        CALL finish( TRIM(method_name),TRIM(message_text))
    ENDIF
#endif

    IF (ltransport) THEN
    DO jg = 1,n_dom


      !----------------------------------------------
      ! Flux computation methods - consistency check

        z_go_tri(1:11)=(/NO_HADV,UP,MIURA,MIURA3,FFSL,FFSL_HYB,MCYCL,       &
          &              MIURA_MCYCL,MIURA3_MCYCL,FFSL_MCYCL,FFSL_HYB_MCYCL/)
        DO jt=1,ntracer
          IF ( ALL(z_go_tri /= advection_config(jg)%ihadv_tracer(jt)) ) THEN
            CALL finish( TRIM(method_name),                                       &
              &  'incorrect settings for TRI-C grid ihadv_tracer. Must be '// &
              &  '0,1,2,3,4,5,6,20,22,32,42 or 52 ')
          ENDIF
        ENDDO

    END DO ! jg = 1,n_dom
    END IF ! ltransport


    !--------------------------------------------------------------------
    ! Horizontal diffusion
    !--------------------------------------------------------------------

    DO jg =1,n_dom

      SELECT CASE( diffusion_config(jg)%hdiff_order )
      CASE(-1)
        WRITE(message_text,'(a,i2.2)') 'Horizontal diffusion '//&
                                       'switched off for domain ', jg
        CALL message(TRIM(method_name),TRIM(message_text))

      CASE(2,3,4,5)
        CONTINUE

      CASE(24,42)
        IF (.NOT.( iequations==IHS_ATM_TEMP)) CALL finish(TRIM(method_name), &
        ' hdiff_order = 24 or 42 only implemented for the hydrostatic atm model')

      CASE DEFAULT
        CALL finish(TRIM(method_name),                       &
          & 'Error: Invalid choice for  hdiff_order. '// &
          & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
      END SELECT

      IF ( diffusion_config(jg)%hdiff_efdt_ratio<=0._wp) THEN
        CALL message(TRIM(method_name),'No horizontal background diffusion is used')
      ENDIF

      IF (lshallow_water)  diffusion_config(jg)%lhdiff_temp=.FALSE.

      IF (itype_comm == 3 .AND. diffusion_config(jg)%hdiff_order /= 5)  &
        CALL finish(TRIM(method_name), 'itype_comm=3 requires hdiff_order = 5')

      IF (itype_comm == 3 .AND. (diffusion_config(jg)%itype_vn_diffu > 1 .OR. &
        diffusion_config(jg)%itype_t_diffu > 1) )                             &
        CALL finish(TRIM(method_name), 'itype_comm=3 requires itype_t/vn_diffu = 1')

    ENDDO

    !--------------------------------------------------------------------
    ! checking the meanings of the io settings
    !--------------------------------------------------------------------

    SELECT CASE(iforcing)
    CASE ( iecham, ildf_echam )
      inextra_2d   = 0
      inextra_3d   = 0

    CASE DEFAULT
    END SELECT


    IF (activate_sync_timers .AND. .NOT. ltimer) THEN
      activate_sync_timers = .FALSE.
      WRITE (message_text,*) &
        & "warning: namelist parameter 'activate_sync_timers' has been set to .FALSE., ", &
        & "because global 'ltimer' flag is disabled."
      CALL message('io_namelist', TRIM(message_text))
    END IF
    IF (timers_level > 9 .AND. .NOT. activate_sync_timers) THEN
      activate_sync_timers = .TRUE.
      WRITE (message_text,*) &
        & "warning: namelist parameter 'activate_sync_timers' has been set to .TRUE., ", &
        & "because global 'timers_level' is > 9."
      CALL message('io_namelist', TRIM(message_text))
    END IF


    !--------------------------------------------------------------------
    ! Realcase runs
    !--------------------------------------------------------------------

    IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN  ! start from dwd analysis with incremental update

      ! check analysis update window
      !
      IF ( (dt_iau > 0._wp) .AND. (dt_iau < dtime)) THEN
        ! If dt_iau is chosen to be larger than 0, it must be >= dtime at least.
        dt_iau = dtime
        WRITE (message_text,'(a,a,f6.2)') "Wrong value for dt_iau. ", &
          &   "If >0 then at least equal to advective/phys tstep ",dtime
        CALL finish('initicon_nml:', TRIM(message_text))
      ENDIF 


      ! IAU modes MODE_IAU_OLD cannot be combined with snowtiles
      ! when performing snowtile warmstart.
      IF ((ntiles_lnd > 1) .AND. (.NOT. ltile_coldstart) .AND. (lsnowtile)) THEN
        IF ( init_mode == MODE_IAU_OLD ) THEN
          WRITE (message_text,'(a,i2)') "lsnowtile=.TRUE. not allowed for IAU-Mode ", init_mode   
          CALL finish(method_name, TRIM(message_text))
        ENDIF
      ENDIF

    ENDIF


    ! check meteogram configuration
    CALL check_meteogram_configuration(num_io_procs)

    CALL land_crosscheck()
  
  END  SUBROUTINE atm_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE land_crosscheck

    CHARACTER(len=*), PARAMETER :: method_name =  'mo_nml_crosscheck:land_crosscheck'

#ifdef __NO_JSBACH__
    IF (echam_phy_config% ljsbach) THEN
      CALL finish(method_name, "This version was compiled without jsbach. Compile with __JSBACH__, or set ljsbach=.FALSE.")
    ENDIF
    echam_phy_config% ljsbach   = .FALSE.     
#else
    IF (echam_phy_config% ljsbach) THEN
      IF (num_restart_procs > 0) THEN
        CALL finish(method_name, "JSBACH currently doesn't work with asynchronous restart. Set num_restart_procs=0 !")
      END IF
      IF (num_io_procs > 0) THEN
        CALL finish(method_name, "JSBACH currently doesn't work with asynchronous IO. Set num_io_procs=0 !")
      END IF
    END IF
#endif

  END SUBROUTINE land_crosscheck
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !> Set time step for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  the date information.
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_timestep_settings()
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_timestep_settings'

    TYPE(timedelta), POINTER             ::  dtime1, dtime2, zero_dt
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) ::  dtime_str2, dtime_string
    INTEGER(i8)                          ::  dtime_ms
    REAL(wp)                             ::  dtime_real
    INTEGER                              ::  jg
    TYPE(datetime), POINTER              ::  reference_dt

    ! --------------------------------------------------------------
    ! PART I: Collect time step as ISO8601 string
    ! --------------------------------------------------------------

    ! The time step is determined either by the "dtime" parameter or
    ! by the "modelTimeStep" parameter from the namelist "run_nml".
    !
    dtime_string = ""
    IF (TRIM(mtime_modelTimeStep) /= "")   dtime_string = TRIM(mtime_modelTimeStep)
    dtime_real = dtime
    IF (dtime_real > 0._wp) THEN
      dtime_ms   = NINT(dtime_real*1000, i8)
      CALL getPTStringFromMS(dtime_ms, dtime_str2)
      
      IF (dtime_string == "") THEN
        dtime_string = TRIM(dtime_str2)
      ELSE
        ! Obviously, both namelist parameters for the time step have
        ! been used. We need to test for equality.
        dtime1 => newTimedelta(dtime_string)
        dtime2 => newTimedelta(dtime_str2)
        IF (dtime1 /= dtime2) THEN
          CALL finish(routine, "Two inconsistent time steps have been specified!")
        END IF
        CALL deallocateTimedelta(dtime1)
        CALL deallocateTimedelta(dtime2)
      END IF
    END IF
    dtime1 => newTimedelta(dtime_string)

    ! Furthermore, the time step may be rescaled by the
    ! "grid_rescale_factor".
    !
    IF (grid_rescale_factor /= 1.0_wp) THEN
      dtime1 = dtime1 * grid_rescale_factor
      CALL timedeltaToString(dtime1, dtime_string)
      IF (dtime_real > 0._wp)  dtime_real = dtime_real * grid_rescale_factor

      IF (get_my_process_type() == atmo_process) THEN
        echam_phy_config%dt_rad = &
          & echam_phy_config%dt_rad * grid_rescale_factor
        
        DO jg=1,max_dom
          atm_phy_nwp_config(jg)%dt_conv = &
            atm_phy_nwp_config(jg)%dt_conv * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_rad  = &
            atm_phy_nwp_config(jg)%dt_rad  * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_sso  = &
            atm_phy_nwp_config(jg)%dt_sso  * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_gwd  = &
            atm_phy_nwp_config(jg)%dt_gwd  * grid_rescale_factor
        END DO
      END IF
    END IF

    ! consistency check
    zero_dt => newTimedelta("PT0S") ! mtime object for zero
    IF (dtime1 <= zero_dt) CALL finish(TRIM(routine),'"dtime" must be positive')

    CALL deallocateTimedelta( dtime1  )
    CALL deallocateTimedelta( zero_dt )


    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 string into "mtime" and old REAL
    ! --------------------------------------------------------------

    CALL setModelTimeStep(dtime_string)
    IF (dtime_real > 0._wp) THEN
      ! In case that we came from the REAL-valued namelist setting of
      ! the time step we try to avoid rounding errors in floating
      ! point arithmetic
      dtime = dtime_real
    ELSE
      ! For conversion to milliseconds, we need an anchor date,
      ! although we know that "dtime" is too small for this to be
      ! relevant.
      reference_dt => newDatetime("1980-06-01T00:00:00.000")
      dtime = REAL(getTotalMilliSecondsTimeDelta(dtime1, reference_dt),wp)/1000._wp
      CALL deallocateDatetime(reference_dt)
    END IF


    ! --------------------------------------------------------------
    ! PART III: Print time step
    ! --------------------------------------------------------------

    CALL message('','')
    WRITE(message_text,'(a,a)') 'Model time step          : ', TRIM(dtime_string)
    CALL message('',message_text)

  END SUBROUTINE compute_timestep_settings


  !---------------------------------------------------------------------------------------
  !> Set restart and checkpoint interval for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  the date information.
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_restart_settings()
    ! local variables:
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_restart_settings'
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: restart_intvl_string, checkpt_intvl_string, &
      &                                     checkpt_intvl2, restart_intvl2, dtime_string
    TYPE(timedelta), POINTER             :: mtime_2_5h, mtime_dt_checkpoint,            &
      &                                     mtime_dt_restart, mtime_dom_start
    TYPE(datetime), POINTER              :: reference_dt
    INTEGER                              :: jg


    ! --------------------------------------------------------------
    ! PART I: Collect the restart and checkpoint intervals 
    !         as ISO8601 strings
    ! --------------------------------------------------------------

    ! --- --- RESTART INTERVAL:
    !
    !         This time interval specifies when the run is supposed to
    !         save its state to a file and stop, later to be resumed.
    !
    !         The restart interval is set by the namelist parameter
    !         "restartTimeIntval" (in the namelist
    !         "master_time_control_nml") - the new way -, and
    !         "dt_restart" (in "time_nml") - the old way.
    !
    !         TODO: The restart interval needs to be multiple of the model
    !               time steps.
    !
    restart_intvl_string = ""
    IF (TRIM(restartTimeIntval) /= "")  restart_intvl_string = TRIM(restartTimeIntval)
    IF (dt_restart > 0._wp) THEN
      restart_intvl2 = "PT"//TRIM(int2string(INT(dt_restart), '(i0)'))//"S"
      IF (TRIM(restart_intvl_string) == "") THEN
        restart_intvl_string = TRIM(restart_intvl2)
      ELSE
        IF (timedelta_str_equal(restart_intvl_string, restart_intvl2)) THEN
          restart_intvl_string = TRIM(restart_intvl2)
        ELSE
          CALL finish(routine, "Inconsistent setting of restart interval: "//&
              &TRIM(restart_intvl_string)//"/"//TRIM(restart_intvl2))
        END IF
      END IF
    END IF
    ! if "restart_intvl_string" still unspecified: set default
    IF (TRIM(restart_intvl_string) == "") THEN
      restart_intvl_string = DEFAULT_RESTART_INTVL
    END IF



    ! --- --- CHECKPOINT INTERVAL:
    !
    !         This time interval specifies when the run is supposed to
    !         save its state to a file (but not to stop afterwards).
    !
    !         The checkpoint interval is set by the namelist parameter
    !         "checkpointTimeIntval" (in the namelist
    !         "master_time_control_nml") - the new way -, and
    !         "dt_checkpoint" (in "io_nml") - the old way.
    !
    !         TODO: The checkpoint interval needs to be multiple of the
    !               model time steps.
    !
    checkpt_intvl_string = ""
    IF (TRIM(checkpointTimeIntval) /= "")  checkpt_intvl_string = TRIM(checkpointTimeIntval)
    IF (dt_checkpoint > 0._wp) THEN
      checkpt_intvl2 = "PT"//TRIM(int2string(INT(dt_checkpoint), '(i0)'))//"S"
      IF (TRIM(checkpt_intvl_string) == "") THEN
        checkpt_intvl_string = TRIM(checkpt_intvl2)
      ELSE
        IF (timedelta_str_equal(checkpt_intvl_string, checkpt_intvl2)) THEN
          checkpt_intvl_string = TRIM(checkpt_intvl2)
        ELSE
          CALL finish(routine, "Inconsistent setting of checkpoint interval: "//&
            &TRIM(checkpt_intvl_string)//"/"//TRIM(checkpt_intvl2))
        END IF
      END IF
    END IF
    ! if "checkpt_intvl_string" still unspecified: set default
    IF (TRIM(checkpt_intvl_string) == "") THEN
      checkpt_intvl_string = DEFAULT_CHECKPT_INTVL
    END IF

    mtime_dt_restart    => newTimedelta(restart_intvl_string)
    mtime_dt_checkpoint => newTimedelta(checkpt_intvl_string)

    ! consistency checks:
    !
    ! When increased sound-wave and gravity-wave damping is chosen
    ! during the spinup phase (i.e. divdamp_order = 24),
    ! checkpointing/restarting is not allowed earlier than three hours
    ! into the integration because the results would not be
    ! bit-identical in this case
    !
    mtime_2_5h          => newTimedelta("PT02H30M")
    IF ((iequations    == inh_atmosphere) .AND. &
      & (divdamp_order == 24)             .AND. &
      & (mtime_dt_checkpoint < mtime_2_5h)) THEN
        WRITE(message_text,'(a)') &
          &  'dt_checkpoint < 2.5 hours not allowed in combination with divdamp_order = 24'
        CALL finish(routine, message_text)
    ENDIF
    CALL deallocateTimedelta(mtime_2_5h)

    ! Writing a checkpoint file exactly at the start time of a nest is
    ! not allowed:
    !
    DO jg =1,n_dom
      dtime_string = "PT"//TRIM(int2string(INT(start_time(jg)), '(i0)'))//"S"
      mtime_dom_start => newTimedelta(dtime_string)
      IF (mtime_dom_start == mtime_dt_checkpoint) THEN
        WRITE(message_text,'(a)') &
          &  'writing a checkpoint file exactly at the start time of a nest is not allowed'
        CALL finish(routine, message_text)
      END IF
      CALL deallocateTimedelta(mtime_dom_start)
    END DO


    ! If undefined, the restart interval is set to the checkpoint
    ! interval. On the other hand, the checkpoint interval is set to
    ! the restart interval, if the first is not specified:
    !
    IF (TRIM(restart_intvl_string) == "")  restart_intvl_string = checkpt_intvl_string
    IF (TRIM(checkpt_intvl_string) == "")  checkpt_intvl_string = restart_intvl_string


    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 string into "mtime" and old REAL
    ! --------------------------------------------------------------

    CALL setRestartTimeInterval   ( restart_intvl_string )
    CALL setCheckpointTimeInterval( checkpt_intvl_string )

    ! For conversion to (milli-)seconds, we need an anchor date:
    !
    reference_dt => newDatetime("1980-06-01T00:00:00.000")
    dt_restart    = REAL(getTotalMilliSecondsTimeDelta(mtime_dt_checkpoint, reference_dt),wp)/1000._wp
    dt_checkpoint = REAL(getTotalMilliSecondsTimeDelta(mtime_dt_restart,    reference_dt),wp)/1000._wp
    CALL deallocateTimedelta(mtime_dt_checkpoint)
    CALL deallocateTimedelta(mtime_dt_restart)
    CALL deallocateDatetime(reference_dt)

    ! --------------------------------------------------------------
    ! PART III: Print restart and checkpoint intervals
    ! --------------------------------------------------------------

    WRITE(message_text,'(a,a)') 'Checkpoint interval      : ', TRIM(checkpt_intvl_string)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Restart interval         : ', TRIM(restart_intvl_string)
    CALL message('',message_text)
    CALL message('','')

  END SUBROUTINE compute_restart_settings


  !---------------------------------------------------------------------------------------
  !> Set start date, end date etc. for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  this date information.
  !
  !  This subroutine checks for contradictory namelist settings and
  !  finally creates 
  !   - data objects of type "t_datetime" 
  !   - data objects from  the mtime library.
  !   - the "nsteps" time loop count
  !
  !  As input this subroutine requires final settings dtime and
  !  dt_restart (these may have been modified in the crosscheck
  !  routine in "mo_nml_crosscheck").
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_date_settings(dt_restart, nsteps, time_config)
    REAL(wp),            INTENT(IN)    :: dt_restart         !< Length of restart cycle in seconds
    INTEGER,             INTENT(INOUT) :: nsteps
    TYPE(t_time_config), INTENT(OUT)   :: time_config
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_date_settings'

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  ini_datetime1, end_datetime1, dstring
    CHARACTER(LEN=32)                     ::  ini_datetime2, end_datetime2
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  start_datetime_string,         & !< run start date
      &                                       stop_datetime_string,          & !< run stop date
      &                                       exp_start_datetime_string,     & !< experiment start date
      &                                       exp_stop_datetime_string,      & !< experiment stop date
      &                                       exp_ref_datetime_string,       & !< experiment reference date
      &                                       cur_datetime_string
    TYPE(datetime),  POINTER              ::  mtime_start, mtime_stop,       &
      &                                       mtime_exp_stop,                &
      &                                       mtime_restart_stop,            &
      &                                       mtime_nsteps_stop
    TYPE(timedelta), POINTER              ::  mtime_dt_restart, mtime_dtime
    INTEGER                               ::  mtime_calendar, dtime_calendar,&
      &                                       errno
    CHARACTER(len=MAX_CALENDAR_STR_LEN)   ::  calendar1, calendar2, calendar


    ! --------------------------------------------------------------
    ! PART I: Collect all the dates as ISO8601 strings
    ! --------------------------------------------------------------

    ! --- Prelimary: set calendars in old and new date module
    !
    !     The calendar may be set by the namelist parameter "calendar"
    !     in "time_nml" or the namelist parameter "calendar" in
    !     "master_time_control_nml".

    calendar1 = TRIM(time_nml_calendar)
    calendar2 = TRIM(master_nml_calendar)
    IF (TRIM(calendar1) /= "")  calendar = calendar1
    IF (TRIM(calendar2) /= "")  calendar = calendar2
    IF ((TRIM(calendar1) /= "") .AND. (TRIM(calendar2) /= "")) THEN
      ! both settings were used; we need to test for equality
      IF (TRIM(tolower(calendar1)) /= TRIM(tolower(calendar2)))  &
        &  CALL finish(routine, "Inconsistent setting of calendar")
    END IF
    SELECT CASE (toLower(calendar))
    CASE ('julian gregorian')
      dtime_calendar  = dtime_julian_gregorian
      mtime_calendar  = -1
      CALL finish(routine, "Julian-Gregorian calendar unsupported by mtime library!")
    CASE ('proleptic gregorian')
      dtime_calendar  = dtime_proleptic_gregorian
      mtime_calendar  = mtime_proleptic_gregorian
    CASE ('365 day year')  
      dtime_calendar  = -1
      mtime_calendar  = mtime_year_of_365_days
      CALL finish(routine, "Year-of-365-days calendar unsupported by datetime module!")
    CASE ('360 day year')  
      dtime_calendar  = dtime_cly360
      mtime_calendar  = mtime_year_of_360_days
    CASE default
      dtime_calendar  = dtime_proleptic_gregorian
      mtime_calendar  = mtime_proleptic_gregorian
      CALL message('','No calendar selected! Use default proleptic Gregorian.')
    END SELECT
    ! setting the calendar for mtime library is needed for subsequent
    ! calculations:
    CALL setCalendar(mtime_calendar) 


    ! --- Now, we merge ISO time stamp strings from the namelists
    !     "time_nml" and "master_time_control_nml". If both namelist
    !     specify the same setting, then we compare the dates and
    !     throw an error in case of inconsistency.  We take also into
    !     account the namelist settings for "nsteps" and "dt_restart".
    !
    !    ***************************  EXPERIMENT RUN   *************
    !    *                                                         *
    !    *                    ***  JOB RUN   ***                   *
    !    *                    *                *                   *
    !  --|--------------------[----------------]-------------------|---------------> (time axis)
    !   exp_start           start            stop                exp_stop
    !

    ! --- --- EXPERIMENT START DATE:
    !
    !         This is the start date for the whole experiment, which
    !         may be much earlier than the start of the current run.
    !
    !         The experiment start date is set by the namelist
    !         parameter "experimentStartDate" (in namelist
    !         "master_time_control_nml") - the new way -, and
    !         "ini_datetime_string" (in "time_nml") - the old way.
    !
    !         Note that if the namelist parameter
    !         "experimentStartDate" (in "master_time_control_nml") is
    !         unspecified, but the namelist parameter
    !         "experimentReferenceDate" (see also below) is set, then
    !         set the experiment start date to the latter.
    !
    exp_start_datetime_string = ""
    ini_datetime1 = TRIM(experimentStartDate)
    ini_datetime2 = TRIM(ini_datetime_string)
    IF (TRIM(ini_datetime1) /= "")  exp_start_datetime_string = ini_datetime1
    IF (TRIM(ini_datetime2) /= "")  exp_start_datetime_string = ini_datetime2
    IF ((TRIM(ini_datetime1) /= "") .AND. (LEN_TRIM(ini_datetime2) > 0)) THEN
      ! both settings were used; we need to test for equality
      IF (.NOT. datetime_str_equal(ini_datetime1, ini_datetime2))  &
        &  CALL finish(routine, "Inconsistent setting of experiment start date: "//&
        &                       TRIM(ini_datetime1)//"/"//TRIM(ini_datetime2))
    END IF
    IF ( (TRIM(ini_datetime1) == "")        .AND.  &
      &  (TRIM(ini_datetime2) == "")        .AND.  &
      &  (TRIM(experimentReferenceDate) /= "")) THEN
      ! set the experiment start date to the reference date
      exp_start_datetime_string = experimentReferenceDate
    END IF
    ! throw an error, if no start date has been specified at all
    IF (TRIM(exp_start_datetime_string) == "") THEN
      CALL finish(routine, "No experiment start date has been set!")
    END IF


    ! --- --- EXPERIMENT STOP DATE:
    !
    !         This is when the whole experiment ends; not that the
    !         current run may stop before this end date has been
    !         reached (and be restarted later).
    !
    !         The experiment stop date is set by the namelist
    !         parameter "experimentStopDate" (in namelist
    !         "master_time_control_nml") - the new way -, and
    !         "end_datetime_string" (in "time_nml") - the old way.
    !
    exp_stop_datetime_string = ""
    end_datetime1 = TRIM(experimentStopDate)
    end_datetime2 = TRIM(end_datetime_string)
    IF (TRIM(end_datetime1) /= "")  exp_stop_datetime_string = end_datetime1
    IF (TRIM(end_datetime2) /= "")  exp_stop_datetime_string = end_datetime2
    IF ((TRIM(end_datetime1) /= "") .AND. (TRIM(end_datetime2) /= "")) THEN
      ! both settings were used; we need to test for equality
      IF (.NOT. datetime_str_equal(end_datetime1, end_datetime2))  &
        &  CALL finish(routine, "Inconsistent setting of experiment stop date: "//&
        &                       TRIM(end_datetime1)//"/"//TRIM(end_datetime2))
    END IF
    ! throw an error, if no start date has been specified at all
    IF (TRIM(exp_stop_datetime_string) == "") THEN
      CALL finish(routine, "No experiment stop date has been set!")
    END IF

    ! --- --- REFERENCE DATE:
    !
    !         This specifies the reference date for the calendar in
    !         use (some kind of "anchor date" on the time line).
    !
    !         If the namelist parameter "experimentReferenceDate (in
    !         "master_time_control_nml") is unspecified, then the
    !         reference date is set to the experiment start date.
    !
    IF (TRIM(experimentReferenceDate) == '') THEN
      exp_ref_datetime_string = exp_start_datetime_string
    ELSE
      exp_ref_datetime_string = experimentReferenceDate
    END IF

    ! --- --- START DATE:
    !
    !         This is the start date of the current run.
    !
    !         a) in case of restart:  start date := restart attribute "tc_startdate"
    !         b) no-restart case:     start date := experiment start date
    !
    !         Special case: In a restart situation, if the calendar or
    !         initial date/time is different from those in the restart
    !         file, we regard this integration as a new one with its
    !         own calendar.  Model time at which the previous run
    !         stopped is thus not relevant.  Simulation will start
    !         from the user-specified initial date/time, which is also
    !         the current model date/time.
    !
    IF (isRestart()) THEN

      IF ((restart_calendar                  /= dtime_calendar) .OR.     &
        & (TRIM(restart_ini_datetime_string) /= TRIM(ini_datetime_string))) THEN

        start_datetime_string = exp_start_datetime_string

      ELSE
        CALL message('','Read restart file meta data ...')
        CALL read_restart_header("atm")      !< TODO: fix this for ocean!!!
        CALL get_restart_attribute('tc_startdate', start_datetime_string)
      END IF

    ELSE
      start_datetime_string = exp_start_datetime_string
    ENDIF

    ! --- --- CURRENT DATE:
    !
    !         This is the current model date, that is stepped forward
    !         during the integration loop. We begin by setting it to
    !         the start date computed above.
    !
    cur_datetime_string = start_datetime_string

    ! --- --- STOP DATE:
    !
    !         This is the date when the current run stops time stepping.
    !
    !         The simulation stops at the next (i.e. earliest) of the
    !         following dates:
    !         a) experiment stop date
    !         b) start date + dt_restart
    !         c) start date + nsteps*dtime
    !
    mtime_dt_restart   => newTimedelta("PT"//TRIM(int2string(INT(dt_restart),'(i0)'))//"S", errno)
    IF (errno /= 0)  CALL finish(routine, "Error in conversion of dt_restart!")
    mtime_dtime        => getTimedeltaFromMS(INT(dtime,i8)*1000)
    IF (.NOT. ASSOCIATED(mtime_dtime))  CALL finish(routine, "Error in conversion of dtime to mtime!")
    mtime_start        => newDatetime(start_datetime_string, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in conversion of start date")
    mtime_exp_stop     => newDatetime(exp_stop_datetime_string, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in conversion of exp stop date")
    mtime_restart_stop => newDatetime(mtime_start, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in initialization of restart date")
    mtime_restart_stop =  mtime_restart_stop + mtime_dt_restart
    IF (nsteps >= 0) THEN   

      ! Special treatment for the hydro atm model
      !
      ! TODO: Is this weird workaround really needed?
      IF ( (iequations == IHS_ATM_TEMP) .OR. &
        &  (iequations == IHS_ATM_THETA)     ) THEN
        
        ! If running the HYDROSTATIC version, let the model integrate
        ! one more step after the desired end of simulation in order
        ! to get the proper output. This additional step is necessary
        ! because the HYDROSTATIC model writes out values of step N
        ! after the integration from N to N+1 is finished. Also note
        ! that this additional step is done only for the regular
        ! output, and is ignored for restart.
        nsteps = nsteps + 1

        ! The additional step is not needed in the NON-hydrostatic
        ! version because in this case the model writes out values of
        ! step N after the integration from N-1 to N is finished.
      END IF
      mtime_nsteps_stop  => newDatetime(mtime_start, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of nsteps  stop date")
      mtime_nsteps_stop = mtime_nsteps_stop + mtime_dtime * INT(nsteps,c_int32_t)
    ELSE
      mtime_nsteps_stop  => newDatetime(mtime_exp_stop, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of nsteps  stop date")
    END IF
    mtime_stop => newDatetime(MIN(MIN(mtime_exp_stop, mtime_restart_stop), mtime_nsteps_stop))

    ! consistency checks:
    !
    CALL datetimeToString(mtime_stop, stop_datetime_string)
    IF (mtime_stop < mtime_start) THEN
      CALL finish(routine, 'The end date and time must not be '// &
        &                  'before the current date and time')
    END IF
    ! If a restart event occurs, check for unsupported combinations of
    ! namelist settings:
    IF (.NOT. (mtime_nsteps_stop < mtime_restart_stop)) THEN
      ! processor splitting cannot be combined with synchronous restart:
      IF ((num_restart_procs == 0) .AND. ANY(patch_weight(1:) > 0._wp)) THEN
        CALL finish(routine, "Processor splitting cannot be combined with synchronous restart!")
      END IF
    END IF
    CALL deallocateDatetime(mtime_start)
    CALL deallocateDatetime(mtime_exp_stop)
    CALL deallocateDatetime(mtime_restart_stop)
    CALL deallocateDatetime(mtime_nsteps_stop)
    CALL deallocateDatetime(mtime_stop)
    CALL deallocateTimedelta(mtime_dt_restart)
    CALL deallocateTimedelta(mtime_dtime)

    ! --- --- NSTEPS
    !
    !         This specifies the number of steps in the time loop
    !         (INTEGER).
    !
    !         If the namelist parameter "nsteps" has not been
    !         specified in the namelist "run_nml", then this value is
    !         computed from the stop date and the start date. On the
    !         other hand, if the namelist parameter "nsteps" is
    !         specified in the namelist "run_nml", then this defines
    !         the end date (see above).
    !
    IF (nsteps < 0) THEN   
      ! User did not specified a value, we need to compute "nsteps" as
      ! (stop date - start date)/dtime:
      nsteps = get_timedelta_divide_by_seconds(start_datetime_string, &
        &                                      stop_datetime_string,  &
        &                                      dtime)
    END IF

    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 strings into "mtime" and "t_datetime" 
    ! --------------------------------------------------------------

    ! --- Second, create date-time objects from the mtime library

    CALL setStartdate   ( start_datetime_string     )
    CALL setStopdate    ( stop_datetime_string      )
    CALL setExpStartdate( exp_start_datetime_string )
    CALL setExpStopdate ( exp_stop_datetime_string  )
    CALL setExpRefdate  ( exp_ref_datetime_string   )
    ! TODO: initialize current date also for mtime here!!!

    ! --- Finally, store the same information in a "t_datetime" data
    !     structure
    time_config%calendar         = dtime_calendar
    time_config%is_relative_time = is_relative_time

    CALL string_to_datetime( start_datetime_string, time_config%ini_datetime )
    time_config%ini_datetime%calendar = dtime_calendar
    CALL date_to_time( time_config%ini_datetime )

    CALL string_to_datetime( end_datetime_string, time_config%end_datetime ) 
    time_config%end_datetime%calendar = dtime_calendar
    CALL date_to_time( time_config%end_datetime )

    CALL string_to_datetime( cur_datetime_string, time_config%cur_datetime )
    time_config%cur_datetime%calendar = dtime_calendar
    CALL date_to_time( time_config%cur_datetime )

    ! --------------------------------------------------------------
    ! PART III: Print all date and time components
    ! --------------------------------------------------------------

    CALL calendarToString(dstring)
    CALL message('','Calendar: '//TRIM(dstring))
    call message('','')

    WRITE(message_text,'(a,a)') 'Experiment reference date: ', TRIM(exp_ref_datetime_string)
    CALL message('',message_text)

    CALL message(' ',' ')
    CALL message(routine,'Initial date and time')
    CALL message(routine,'---------------------')
    WRITE(message_text,'(a,a)') 'Experiment start date    : ', TRIM(start_datetime_string)
    CALL message('',message_text)

    CALL message(' ',' ')
    CALL message(routine,'End date and time')
    CALL message(routine,'-----------------')
    WRITE(message_text,'(a,a)') 'Experiment stop date     : ', TRIM(exp_stop_datetime_string)
    CALL message('',message_text)
  END SUBROUTINE compute_date_settings

END MODULE mo_nml_crosscheck
