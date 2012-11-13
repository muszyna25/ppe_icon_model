!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Kristina Froehlich, MPI-M (2011-07-12)
!! @author Hui Wan, MPI-M (2011-07-12)
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_nml_crosscheck

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish, print_value
  USE mo_impl_constants,     ONLY: max_char_length, max_dom,itconv,itccov,    &
    &                              itrad,itradheat, itsso, itgscp, itsatad,   &
    &                              itupdate, itturb, itsfc, itgwd, iphysproc, &
    &                              iecham, ildf_echam, inwp, iheldsuarez,     &
    &                              ildf_dry, inoforcing, ihs_atm_temp,        &
    &                              ihs_atm_theta, tracer_only, inh_atmosphere,&
    &                              ishallow_water, LEAPFROG_EXPL, LEAPFROG_SI,&
    &                              NO_HADV, UP, MIURA, MIURA3, FFSL, UP3,     &
    &                              MCYCL, MIURA_MCYCL, MIURA3_MCYCL,          &
    &                              ifluxl_sm, ifluxl_m, ihs_ocean,            &
    &                              RAYLEIGH_CLASSIC 
  USE mo_time_config,        ONLY: time_config, restart_experiment
  USE mo_extpar_config,      ONLY: itopo
  USE mo_io_config,          ONLY: dt_checkpoint, lflux_avg,inextra_2d,       &
    &                              inextra_3d, lwrite_cloud, lwrite_extra,    &
    &                              lwrite_omega, lwrite_precip, lwrite_pres,  &
    &                              lwrite_radiation,lwrite_surface,lwrite_tend_phy,& 
    &                              lwrite_tke,lwrite_z3
  USE mo_parallel_config,    ONLY: check_parallel_configuration,              &
    &                              num_io_procs, itype_comm
  USE mo_run_config,         ONLY: lrestore_states, nsteps, dtime, iforcing,  &
    &                              ltransport, ntracer, nlev, ltestcase,      &
    &                              nqtendphy, iqv, iqc, iqi,                  &
    &                              iqs, iqr, iqt, iqtvar, ico2, ltimer,       &
    &                              activate_sync_timers, timers_level,        &
    &                              output_mode, dtime_adv
  USE mo_gridref_config
  USE mo_interpol_config
  USE mo_grid_config
  USE mo_sleve_config

  USE mo_dynamics_config,    ONLY: configure_dynamics,                 &
    &                              iequations, idiv_method,            &
    &                              divavg_cntrwgt, sw_ref_height,      &
    &                              lcoriolis, lshallow_water, ltwotime
  USE mo_advection_config,   ONLY: advection_config, configure_advection

  USE mo_nonhydrostatic_config, ONLY: itime_scheme_nh => itime_scheme, iadv_rcf, &
                                      lhdiff_rcf, rayleigh_type
  USE mo_ha_dyn_config,      ONLY: ha_dyn_config
  USE mo_diffusion_config,   ONLY: diffusion_config, configure_diffusion


  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, configure_atm_phy_nwp
  USE mo_lnd_nwp_config,     ONLY: nlev_soil, nlev_snow ,ntiles_total,  &
    &                              lseaice,  llake, lmelt , lmelt_var, lmulti_snow
  USE mo_lnd_jsbach_config,  ONLY: lnd_jsbach_config, configure_lnd_jsbach
  USE mo_echam_phy_config,   ONLY: echam_phy_config, configure_echam_phy
  USE mo_radiation_config
  USE mo_echam_conv_config,  ONLY: echam_conv_config, configure_echam_convection
  USE mo_gw_hines_config,    ONLY: gw_hines_config
  USE mo_vdiff_config,       ONLY: vdiff_config
  USE mo_turbdiff_config,    ONLY: turbdiff_config
  USE mo_nh_testcases,       ONLY: linit_tracer_fv,nh_test_name
  USE mo_ha_testcases,       ONLY: ctest_name, ape_sst_case

  USE mo_datetime,           ONLY: add_time, print_datetime_all
  USE mo_meteogram_config,   ONLY: check_meteogram_configuration
  USE mo_master_control,     ONLY: is_restart_run, get_my_process_type,      &
    & testbed_process,  atmo_process, ocean_process, radiation_process
  
  USE mo_art_config,         ONLY: art_config

  IMPLICIT NONE

!  PRIVATE

  PUBLIC :: atm_crosscheck, oce_crosscheck

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  !>
  !! Check and, if necessary, adapt simulation lengths
  !!
  !! Check and, if necessary, adapt dt_restart and dt_checkpoint. 
  !! - i.e. makes sure that dt_restart and dt_checkpoint are synchronized 
  !! with a transport event.
  !!
  !!
  !! @par Revision History
  !! Initial revision by Hui Wan, MPI (2011-07)
  !!
  SUBROUTINE resize_simulation_length()

    REAL(wp):: cur_datetime_calsec, end_datetime_calsec, length_sec
    INTEGER :: jg
    CHARACTER(len=*), PARAMETER :: routine =  'resize_simulation_length'
    
    !----------------------------
    ! rescale timestep
    dtime     = dtime     * grid_rescale_factor
    IF (get_my_process_type() == atmo_process) THEN
      dtime_adv = dtime_adv * grid_rescale_factor
      echam_phy_config%dt_rad = &
        & echam_phy_config%dt_rad * grid_rescale_factor
        
      DO jg=1,max_dom
        atm_phy_nwp_config(jg)%dt_conv = &
          atm_phy_nwp_config(jg)%dt_conv * grid_rescale_factor
        atm_phy_nwp_config(jg)%dt_ccov = &
          atm_phy_nwp_config(jg)%dt_ccov * grid_rescale_factor
        atm_phy_nwp_config(jg)%dt_rad  = &
          atm_phy_nwp_config(jg)%dt_rad  * grid_rescale_factor
        atm_phy_nwp_config(jg)%dt_sso  = &
          atm_phy_nwp_config(jg)%dt_sso  * grid_rescale_factor
        atm_phy_nwp_config(jg)%dt_gwd  = &
          atm_phy_nwp_config(jg)%dt_gwd  * grid_rescale_factor
      ENDDO
    ENDIF
    !--------------------------------------------------------------------
    ! Length if this integration
    !--------------------------------------------------------------------
    IF (nsteps/=0) THEN   ! User specified a value

      length_sec = REAL(nsteps,wp)*dtime
      time_config%end_datetime = time_config%cur_datetime
      CALL add_time(length_sec,0,0,0,time_config%end_datetime)

   !HW (2011-07-17): run_day/hour/... not implemented in the restructured version ------
   !ELSE IF (run_day/=0 .OR. run_hour/=0 .OR. run_minute/=0 .OR. run_second/=0.0_wp) THEN
   !  IF (run_day    < 0    ) CALL finish(routine,'"run_day" must not be negative')
   !  IF (run_hour   < 0    ) CALL finish(routine,'"run_hour" must not be negative')
   !  IF (run_minute < 0    ) CALL finish(routine,'"run_minute" must not be negative')
   !  IF (run_second < 0._wp) CALL finish(routine,'"run_second" must not be negative')
   !  !
   !  end_datetime = cur_datetime
   !  CALL add_time(run_second,run_minute,run_hour,run_day,end_datetime)
   !  !
   !  cur_datetime_calsec = (REAL(cur_datetime%calday,wp)+cur_datetime%caltime) &
   !    &                   *REAL(cur_datetime%daylen,wp)
   !  end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
   !    &                   *REAL(end_datetime%daylen,wp)
   !  nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
   !-------------------------

    ELSE
      ! Compute nsteps from cur_datetime, end_datetime and dtime
      !
      cur_datetime_calsec = (REAL(time_config%cur_datetime%calday,wp)  &
                                 +time_config%cur_datetime%caltime   ) &
                           * REAL(time_config%cur_datetime%daylen,wp)
      end_datetime_calsec = (REAL(time_config%end_datetime%calday,wp)  &
                                 +time_config%end_datetime%caltime   ) &
                           * REAL(time_config%end_datetime%daylen,wp)

      IF (end_datetime_calsec < cur_datetime_calsec) &
        & CALL finish(TRIM(routine),'The end date and time must not be '// &
        &            'before the current date and time')

      nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)

    END IF


    IF (iequations/=ihs_ocean) THEN ! atm (ocean does not know iadv_rcf) 
      ! Check whether the end of the restart cycle is synchronized with a transport 
      ! event. If not, adapt dt_restart accordingly.
      ! dtime_adv not available at this point. Thus we need to use iadv_rcf*dtime
      !
      IF (MOD(time_config%dt_restart,REAL(iadv_rcf,wp)*dtime) /= 0) THEN
        time_config%dt_restart =                                              &
          &   REAL(NINT(time_config%dt_restart/(REAL(iadv_rcf,wp)*dtime)),wp) &
          &   * REAL(iadv_rcf,wp)*dtime
        WRITE(message_text,'(a)') &
          &  'length of restart cycle dt_restart synchronized with transport event' 
        CALL message(routine, message_text)
      ENDIF
    ENDIF

 
    ! Length of this integration is limited by length of the restart cycle.
    !
    IF (nsteps > INT(time_config%dt_restart/dtime)) THEN
      nsteps = INT(time_config%dt_restart/dtime)
      restart_experiment = .TRUE.
    ELSE
      restart_experiment = .FALSE.
    ENDIF    
!     nsteps = MIN(nsteps,INT(time_config%dt_restart/dtime))


    CALL message(' ',' ')
    CALL message(routine,'Initial date and time')
    CALL message(routine,'---------------------')
    CALL print_datetime_all(time_config%ini_datetime)  ! print all date and time components

    CALL message(' ',' ')
    CALL message(routine,'End date and time')
    CALL message(routine,'-----------------')
    CALL print_datetime_all(time_config%end_datetime)  ! print all date and time components

    CALL message(' ',' ')
    CALL message(routine,'Length of restart cycle')
    CALL message(routine,'-----------------------')
    WRITE(message_text,'(a,f10.2,a,f16.10,a)') &
         &'dt_restart :',time_config%dt_restart,' seconds =', &
         & time_config%dt_restart/86400._wp, ' days'
    CALL message(routine,message_text)


    ! Reset the value of dt_checkpoint if it is longer than dt_restart
    ! so that at least one restart file is generated at the end of the cycle.
    !
    dt_checkpoint = MIN(dt_checkpoint,time_config%dt_restart)


    IF (iequations/=ihs_ocean) THEN ! atm (ocean does not know iadv_rcf) 
      ! Check whether checkpointing is synchronized with a transport event.
      ! If not, adapt dt_checkpoint accordingly.
      ! dtime_adv not available at this point. Thus we need to use iadv_rcf*dtime
      !
      IF (MOD(dt_checkpoint,REAL(iadv_rcf,wp)*dtime) /= 0) THEN
        dt_checkpoint = REAL(NINT(dt_checkpoint/(REAL(iadv_rcf,wp)*dtime)),wp) &
          &           * REAL(iadv_rcf,wp)*dtime
        WRITE(message_text,'(a)') &
          &  'length of checkpoint cycle dt_checkpoint synchronized with transport event' 
        CALL message(routine, message_text)
      ENDIF
    ENDIF

    WRITE(message_text,'(a,f10.2,a,f16.10,a)')          &
         &'dt_checkpoint :',dt_checkpoint,' seconds =', &
         & dt_checkpoint/86400._wp, ' days'
    CALL message(routine,message_text)

  END SUBROUTINE resize_simulation_length


  SUBROUTINE oce_crosscheck()
    CALL check_parallel_configuration()
    CALL resize_simulation_length()
  END SUBROUTINE oce_crosscheck


  SUBROUTINE atm_crosscheck

    INTEGER :: jg
    INTEGER :: jt   ! tracer loop index
    INTEGER :: i_listlen
    INTEGER :: z_go_hex(3), z_go_tri(8), z_nogo_tri(2)   ! for crosscheck
    REAL(wp):: cur_datetime_calsec, end_datetime_calsec, length_sec
    CHARACTER(len=*), PARAMETER :: routine =  'atm_crosscheck'

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()

    CALL resize_simulation_length()

    ! Special treatment for the hydro atm model

    IF ( (iequations == IHS_ATM_TEMP) .OR. &
         (iequations == IHS_ATM_THETA)     ) THEN

      ! If running the HYDROSTATIC version,
      ! let the model integrate one more step after the desired end of
      ! simulation in order to get the proper output. This additional step is
      ! necessary because the HYDROSTATIC model writes out values of step N
      ! after the integration from N to N+1 is finished. Also note that
      ! this additional step is done only for the regular output, and is
      ! ignored for restart.

      nsteps = nsteps + 1

      ! The additional step is not needed in the NON-hydrostatic version because
      ! in this case the model writes out values of step N
      ! after the integration from N-1 to N is finished.
    ENDIF

    !--------------------------------------------------------------------
    ! Horizontal interpolation
    !--------------------------------------------------------------------
    IF (global_cell_type == 6) THEN
    ! ... check i_cori_method
      IF (i_cori_method <1 .OR. i_cori_method>4) THEN
        CALL finish( TRIM(routine),'value of i_cori_method out of range [1,2,3,4]')
      ENDIF
    ENDIF


    !--------------------------------------------------------------------
    ! Grid and dynamics
    !--------------------------------------------------------------------
    IF (lplane .AND. global_cell_type==3) CALL finish( TRIM(routine),&
      'Currently only the hexagon model can run on a plane')

    IF (global_cell_type==6.AND.idiv_method==2) THEN
      CALL finish( TRIM(ROUTINE),'idiv_method =2 not valid for the hexagonal model')
    ENDIF

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
    IF(global_cell_type==6) THEN
      linit_tracer_fv  = .FALSE.
    ENDIF

    IF ((TRIM(ctest_name)=='GW') .AND. (nlev /= 20)) THEN
      CALL finish(TRIM(routine),'nlev MUST be 20 for the gravity-wave test case')
    ENDIF

    IF ((TRIM(ctest_name)=='SV') .AND. ntracer /= 2 ) THEN
      CALL finish(TRIM(routine), &
        & 'ntracer MUST be 2 for the stationary vortex test case')
    ENDIF 

    IF ((TRIM(ctest_name)=='DF1') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(routine), &
        & 'ntracer MUST be >=2 for the deformational flow test case 1')
    ENDIF 

    IF ((TRIM(ctest_name)=='DF2') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(routine), &
        & 'ntracer MUST be >=2 for the deformational flow test case 2')
    ENDIF 

    IF ((TRIM(ctest_name)=='DF3') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(routine), &
        & 'ntracer MUST be >=2 for the deformational flow test case 3')
    ENDIF 

    IF ((TRIM(ctest_name)=='DF4') .AND. ntracer == 1 ) THEN
      CALL finish(TRIM(routine), &
        & 'ntracer MUST be >=2 for the deformational flow test case 4')
    ENDIF     

    IF ((TRIM(ctest_name)=='APE') .AND. (TRIM(ape_sst_case)=='sst_ice')  ) THEN
      IF (.NOT. lflux_avg)&
      CALL finish(TRIM(routine), &
        & 'lflux_avg must be set true to run this setup')
    ENDIF     

    !--------------------------------------------------------------------
    ! Testcases (nonhydrostatic)
    !--------------------------------------------------------------------
    IF (.NOT. ltestcase .AND. rayleigh_type == RAYLEIGH_CLASSIC) THEN
      CALL finish(TRIM(routine), &
        & 'rayleigh_type = RAYLEIGH_CLASSIC not applicable to real case runs.')
    ENDIF

    IF ((TRIM(nh_test_name)=='APE_nh'.OR. TRIM(nh_test_name)=='dcmip_tc_52') .AND.  &
      & ( ANY(atm_phy_nwp_config(:)%inwp_surface == 1 ) )) THEN
      CALL finish(TRIM(routine), &
        & 'surface scheme must be switched off, when running the APE test')
    ENDIF     


    !--------------------------------------------------------------------
    ! Shallow water
    !--------------------------------------------------------------------
    IF (iequations==ISHALLOW_WATER.AND.ha_dyn_config%lsi_3d) THEN
      CALL message( TRIM(routine), 'lsi_3d = .TRUE. not applicable to shallow water model')
    ENDIF

    IF ((iequations==ISHALLOW_WATER).AND.(nlev/=1)) &
    CALL finish(TRIM(routine),'Multiple vertical level specified for shallow water model')

    !--------------------------------------------------------------------
    ! Hydrostatic atm
    !--------------------------------------------------------------------
    IF (iequations==IHS_ATM_THETA) ha_dyn_config%ltheta_dyn = .TRUE.

    !--------------------------------------------------------------------
    ! Nonhydrostatic atm
    !--------------------------------------------------------------------
    IF (lhdiff_rcf .AND. (itype_comm == 3)) CALL finish(TRIM(routine), &
      'lhdiff_rcf is available only for idiv_method=1 and itype_comm<=2')

    !--------------------------------------------------------------------
    ! Atmospheric physics, general
    !--------------------------------------------------------------------
    IF ((iforcing==INWP).AND.(iequations/=INH_ATMOSPHERE)) &
    CALL finish( TRIM(routine), 'NWP physics only implemented in the '//&
               'nonhydrostatic atm model')

    IF ((iforcing==IECHAM).AND.(iequations==INH_ATMOSPHERE)) &
    CALL finish( TRIM(routine), 'ECHAM physics not implemented in the '//&
               'nonhydrostatic atm model')

    !--------------------------------------------------------------------
    ! NWP physics
    !--------------------------------------------------------------------
    IF (iforcing==inwp) THEN
    
 !     CALL configure_atm_phy_nwp(n_dom,ltestcase)
 
      DO jg =1,n_dom

        IF( atm_phy_nwp_config(jg)%inwp_satad == 0       .AND. &
          & ((atm_phy_nwp_config(jg)%inwp_convection >0 ) .OR. &
          &  (atm_phy_nwp_config(jg)%inwp_gscp > 0      )    ) ) &
        &  CALL finish( TRIM(routine),'satad has to be switched on')


        IF( (atm_phy_nwp_config(jg)%inwp_gscp==0) .AND. &
          & (atm_phy_nwp_config(jg)%inwp_convection==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_radiation==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_sso==0)  .AND. &
          & (atm_phy_nwp_config(jg)%inwp_surface == 0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_turb> 0) )   &
        CALL message(TRIM(routine),' WARNING! NWP forcing set but '//&
                    'only turbulence selected!')


        IF( (atm_phy_nwp_config(jg)%inwp_turb == 1) .AND.                &
          & (turbdiff_config(jg)%lconst_z0) )               THEN
          CALL message(TRIM(routine),' WARNING! NWP forcing set but '//  &
                      'idealized (horizontally homogeneous) roughness '//&
                      'length z0 selected!')
        ENDIF 

        ! check radiation scheme in relation to chosen ozone and irad_aero=6 to itopo

        IF (  atm_phy_nwp_config(jg)%inwp_radiation > 0 )  THEN

          SELECT CASE (irad_o3)
          CASE (0,4,6,7) ! ok
          CASE default
            CALL finish(TRIM(routine),'irad_o3 currently has to be 0 , 4, 6, or 7.')
          END SELECT

          ! Tegen aerosol and itopo (Tegen aerosol data have to be read from external data file)
          IF ( ( irad_aero == 6 ) .AND. ( itopo /=1 ) ) THEN
            CALL finish(TRIM(routine),'irad_aero=6 requires itopo=1')
          ENDIF
          
        ELSE

          SELECT CASE (irad_o3)
          CASE (4,6,7)
            irad_o3 = 0
            CALL message(TRIM(routine),'running without radiation => irad_o3 reset to 0')
          END SELECT
          
        ENDIF !inwp_radiation
        
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
    CASE (IECHAM,ILDF_ECHAM)

      IF (ntracer < 3) &
      CALL finish(TRIM(routine),'ECHAM physics needs at least 3 tracers')

      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      ico2   = 4     !! CO2
      iqt    = 4     !! starting index of non-water species 
      nqtendphy = 0  !! number of water species for which convective and turbulent 
                     !! tendencies are stored

    CASE (INWP)

      ! If iforcing=INWP is chosen, ntracer will automatically be set to 5,
      ! except for the case of ICON-ART. If ICON-ART: then ntracer is adapted from the namelist.
      !
      do jg=1,n_dom
        IF ( (ntracer /= 5) .AND. (.NOT. art_config(jg)%lart) ) THEN
          ntracer = 5
          WRITE(message_text,'(a,i3)') 'Attention: for NWP physics, '//&
                                       'ntracer is set to',ntracer
          CALL message(TRIM(routine),message_text)
        ENDIF
        IF ( (ntracer /= 6) .AND. (atm_phy_nwp_config(jg)%inwp_turb == 3) &
        & .AND. (.NOT. art_config(jg)%lart) ) THEN
          ntracer = 6
          WRITE(message_text,'(a,i3)') 'Attention: for NWP physics, '//&
                                       'ntracer is set to',ntracer
          CALL message(TRIM(routine),message_text)
        ENDIF
      enddo
      
      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqr    = 4     !! rain water
      iqs    = 5     !! snow
      iqtvar = 6     !! qt variance
      
      ! Note: Indices for additional tracers are assigned automatically 
      ! via add_tracer_ref in mo_nonhydro_state.

      iqt    = 6     !! start index of other tracers than hydrometeors
      nqtendphy = 3  !! number of water species for which convective and turbulent 
                     !! tendencies are stored
     
    CASE default

      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      ico2   = 5     !! CO2
      iqt    = 4     !! starting index of non-water species
      nqtendphy = 0  !! number of water species for which convective and turbulent 
                     !! tendencies are stored

    END SELECT


    IF (ltransport) THEN
    DO jg = 1,n_dom

      i_listlen = LEN_TRIM(advection_config(jg)%ctracer_list)

      SELECT CASE ( iforcing )
      CASE ( INWP )
      !...........................................................
      ! in NWP physics
      !...........................................................

        IF ( i_listlen /= ntracer ) THEN
          DO jt=1,ntracer
            WRITE(advection_config(jg)%ctracer_list(jt:jt),'(i1.1)')jt
          ENDDO
          WRITE(message_text,'(a,a)') &
            & 'Attention: according to physics, ctracer_list is set to ',&
            & advection_config(jg)%ctracer_list(1:ntracer)
          CALL message(TRIM(routine),message_text)
        ENDIF


      CASE (inoforcing, iheldsuarez, iecham, ildf_dry, ildf_echam)
      !...........................................................
      ! Other types of adiabatic forcing
      !...........................................................
      
        IF ( i_listlen < ntracer .AND. i_listlen /= 0 ) THEN
          ntracer = i_listlen
          CALL message(TRIM(routine),'number of tracers is adjusted according to given list')
        END IF

      
        IF ((iforcing==IECHAM).AND.(echam_phy_config%lrad)) THEN
          IF ( izenith > 4)  &
            CALL finish(TRIM(routine), 'Coose a valid case for rad_nml: izenith.')
        ENDIF
      END SELECT ! iforcing

    END DO ! jg = 1,n_dom
    END IF ! ltransport


    !--------------------------------------------------------------------
    ! Tracer transport
    !--------------------------------------------------------------------
    ! General

    ! similar check already performed in read_run_namelist.
    !
    !IF(ltransport .AND. ntracer <= 0) THEN
    !  CALL finish( TRIM(routine),'Tracer transport switched on but ntracer <= 0')
    !ENDIF

    !IF (.NOT.ltransport .AND. ntracer > 0) &
    !  CALL finish( TRIM(routine),          &
    !  'either set ltransport = true or ntracer to 0 ')

    SELECT CASE (iequations)
    CASE (INH_ATMOSPHERE)

      IF ((itime_scheme_nh==tracer_only) .AND. (.NOT.ltransport)) THEN
        WRITE(message_text,'(A,i2,A)') &
          'nonhydrostatic_nml:itime_scheme set to ', tracer_only, &
          '(TRACER_ONLY), but ltransport to .FALSE.'
        CALL finish( TRIM(routine),TRIM(message_text))
      END IF

    CASE (IHS_ATM_TEMP,IHS_ATM_THETA,ISHALLOW_WATER)

      IF ( (ha_dyn_config%itime_scheme==tracer_only).AND. &
           (.NOT.ltransport)) THEN
        WRITE(message_text,'(A,i2,A)') &
          'ha_dyn_nml:itime_scheme set to ', tracer_only, &
          '(TRACER_ONLY), but ltransport to .FALSE.'
        CALL finish( TRIM(routine),TRIM(message_text))
      END IF

    END SELECT

    IF (ltransport) THEN
    DO jg = 1,n_dom

      !---------------------------------------
      ! Special check for the hexagonal model

      SELECT CASE (global_cell_type)
      CASE (6)
        WRITE(0,*)'nml_Crosscheck'
        ! 3rd order upwind horizontal advection scheme
        advection_config(jg)%ihadv_tracer(:) = UP3
        ! semi monotonous flux limiter
        advection_config(jg)%itype_hlimit(:) = ifluxl_sm
        WRITE(message_text,'(a,i2)') 'NOTE: For hex grid works only ihadv_tracer =',UP3
        CALL message(TRIM(routine),message_text)
        WRITE(message_text,'(a,i2)')' and itype_hlimit= ', ifluxl_sm
        CALL message(TRIM(routine),message_text)
      END SELECT

      !----------------------------------------------
      ! Flux compuation methods - consistency check

      SELECT CASE (global_cell_type)
      CASE (3)
        z_go_tri(1:8)=(/NO_HADV,UP,MIURA,MIURA3,FFSL,MCYCL,MIURA_MCYCL,MIURA3_MCYCL/)
        DO jt=1,ntracer
          IF ( ALL(z_go_tri /= advection_config(jg)%ihadv_tracer(jt)) ) THEN
            CALL finish( TRIM(routine),                                       &
              &  'incorrect settings for TRI-C grid ihadv_tracer. Must be '// &
              &  '0,1,2,3,4,20,22, or 32 ')
          ENDIF
        ENDDO

        IF ( ntracer > 1 ) THEN
          z_nogo_tri(1:2)=(/MIURA_MCYCL,MIURA3_MCYCL/)
          DO jt=2,ntracer
            IF ( ANY(z_nogo_tri == advection_config(jg)%ihadv_tracer(jt)) ) THEN
              CALL finish( TRIM(routine),                                       &
                &  'TRI-C grid ihadv_tracer: MIURA(3)_MCYCL not allowed for '// &
                &  'any other tracer than qv.')
            ENDIF
          ENDDO
        ENDIF

      CASE (6)
        z_go_hex(1:3) = (/NO_HADV,UP,UP3/)
        DO jt=1,ntracer
          IF ( ALL(z_go_hex /= advection_config(jg)%ihadv_tracer(jt)) ) THEN
            CALL finish( TRIM(routine),                                       &
               'incorrect settings for HEX-C grid ihadv_tracer. Must be 0,1, or 5 ')
          ENDIF
        ENDDO

      END SELECT


      !----------------------------------------------
      ! Limiter - consistency check

      IF (global_cell_type == 6) THEN
        IF ( ANY(advection_config(jg)%itype_hlimit(1:ntracer) == ifluxl_m )) THEN
          CALL finish( TRIM(routine),                                     &
           'incorrect settings for itype_hlimit and hexagonal grid. Must be 0 or 4 ')
        ENDIF
      ENDIF

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
        CALL message(TRIM(routine),TRIM(message_text))

      CASE(2,3,4)
        CONTINUE

      CASE(5)
        IF (global_cell_type==6) CALL finish(TRIM(routine), &
        ' hdiff_order = 5 invalid for hexagonal model.')

      CASE(24,42)
        IF (.NOT.( iequations==IHS_ATM_TEMP)) CALL finish(TRIM(routine), &
        ' hdiff_order = 24 or 42 only implemented for the hydrostatic atm model')

      CASE DEFAULT
        CALL finish(TRIM(routine),                       &
          & 'Error: Invalid choice for  hdiff_order. '// &
          & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
      END SELECT

      IF ( diffusion_config(jg)%hdiff_efdt_ratio<=0._wp) THEN
        CALL message(TRIM(routine),'No horizontal background diffusion is used')
      ENDIF

      IF (lshallow_water)  diffusion_config(jg)%lhdiff_temp=.FALSE.

      IF (itype_comm == 3 .AND. diffusion_config(jg)%hdiff_order /= 5)  &
        CALL finish(TRIM(routine), 'itype_comm=3 requires hdiff_order = 5')
 
      IF (itype_comm == 3 .AND. (diffusion_config(jg)%itype_vn_diffu > 1 .OR. &
        diffusion_config(jg)%itype_t_diffu > 1) )                             &
        CALL finish(TRIM(routine), 'itype_comm=3 requires itype_t/vn_diffu = 1')

    ENDDO

    !--------------------------------------------------------------------
    ! checking the meanings of the io settings
    !--------------------------------------------------------------------
    IF (iequations==ISHALLOW_WATER) THEN
       lwrite_omega     = .FALSE.
       lwrite_pres      = .FALSE.
       lwrite_z3        = .FALSE.
    END IF

    IF (iequations==INH_ATMOSPHERE) THEN
       lwrite_omega     = .FALSE.
    END IF

    SELECT CASE(iforcing)
    CASE ( inwp )
      ! Do nothing. Keep the initial values, if not specified in namelist.
      ! consider special idealized testcase with turbulence only
      IF( .NOT. ltransport  )   THEN
        lwrite_precip    = .FALSE.
        lwrite_cloud     = .FALSE.
        lwrite_radiation = .FALSE.
        lwrite_tke       = .TRUE.
        lwrite_surface   = .FALSE.
        CALL message('io_nml_setup',' ATTENTION! Only TKE output for TURBULENCE ONLY test')
      ENDIF

    CASE ( iecham, ildf_echam )
      lwrite_extra = .FALSE.
      inextra_2d   = 0
      inextra_3d   = 0

    CASE (inoforcing,iheldsuarez,ildf_dry)
       lwrite_tend_phy  = .FALSE.
       lwrite_radiation = .FALSE.                                                         
       lwrite_precip    = .FALSE.                                                         
       lwrite_cloud     = .FALSE.                                                        
       lwrite_tke       = .FALSE.                                                        
       lwrite_surface   = .FALSE.                                                        
    CASE DEFAULT
    END SELECT
  
    IF (( inextra_2D > 0) .OR. (inextra_3D > 0) ) THEN
      lwrite_extra = .TRUE.
      WRITE(message_text,'(a,2I4,a,L4)') &
        &'inextra is',inextra_2d,inextra_3d ,' lwrite_extra has been set', lwrite_extra
      CALL message('io_namelist', TRIM(message_text))
    ENDIF
  
    IF (inextra_2D == 0 .AND. inextra_3D == 0 .AND. lwrite_extra) &
      CALL finish('io_namelist','need to specify number of fields for extra output')

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

    ! check meteogram configuration
    CALL check_meteogram_configuration(num_io_procs)

    !---------------------------------------------------------------
    ! Restart runs are disabled for the "old" asynchronous output
    ! mode implemented in MODULE mo_io_vlist
    ! (potential deadlock observed).
    !---------------------------------------------------------------

    IF (is_restart_run() .AND. output_mode%l_vlist .AND. (num_io_procs>0)) THEN
      CALL finish('atm_crosscheck', &
        &         'Restart runs are disabled for the "old" asynchronous output!')
    END IF

  END  SUBROUTINE atm_crosscheck

!  SUBROUTINE atmospheric_configuration
!  INTEGER :: jg
!  CHARACTER(len=*), PARAMETER :: routine =  'atm_setup'
!  END SUBROUTINE atmospheric_configuration

END MODULE mo_nml_crosscheck

