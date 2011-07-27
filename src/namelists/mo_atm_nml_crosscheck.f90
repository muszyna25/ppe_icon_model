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
MODULE mo_atm_nml_crosscheck

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish, print_value
  USE mo_impl_constants,      ONLY: max_char_length, max_dom,itconv,itccov,    &
    &                               itrad,itradheat, itsso, itgscp, itsatad,   &
    &                               itupdate, itturb, itsfc, itgwd, iphysproc, &
    &                               iecham, ildf_echam, inwp, iheldsuarez,     &
    &                               ildf_dry, inoforcing, ihs_atm_temp,        &
    &                               ihs_atm_theta, tracer_only, inh_atmosphere,&
    &                               ishallow_water, LEAPFROG_EXPL, LEAPFROG_SI,&
    &                               iup3, ifluxl_sm, islopel_m, islopel_sm,    &
    &                               ifluxl_m  
  USE mo_time_config,         ONLY: time_config
  USE mo_parallel_config, ONLY: check_parallel_configuration
  USE mo_run_config,          ONLY: lrestore_states, nsteps, dtime, iforcing,  &
    &                               ltransport, ntracer, nlev, io3, ltestcase, &
    &                               iqcond, ntracer_static,&
    &                               iqv, iqc, iqi, iqs, iqr, iqcond, iqt, io3, &
    &                               ico2
                                  
  USE mo_io_config
  USE mo_gridref_config
  USE mo_interpol_config
  USE mo_grid_config
  USE mo_sleve_config

  USE mo_dynamics_config,     ONLY: configure_dynamics,                 &
    &                               iequations, idiv_method,            &
    &                               divavg_cntrwgt, sw_ref_height,      &
    &                               lcoriolis, lshallow_water, ltwotime
  USE mo_advection_config,    ONLY: advection_config, configure_advection

  USE mo_nonhydrostatic_config, ONLY: itime_scheme_nh => itime_scheme
  USE mo_ha_dyn_config,     ONLY: ha_dyn_config
  USE mo_diffusion_config,  ONLY: diffusion_config, configure_diffusion


  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, tcall_phy, &
    &                              configure_atm_phy_nwp
  USE mo_lnd_nwp_config,     ONLY: nlev_soil, nztlev ,nlev_snow ,nsfc_subs,&
    &                              lseaice,  llake, lmelt , lmelt_var, lmulti_snow
  USE mo_echam_phy_config,   ONLY: echam_phy_config, configure_echam_phy
  USE mo_radiation_config
  USE mo_echam_conv_config,  ONLY: echam_conv_config, configure_echam_convection
  USE mo_gw_hines_config,    ONLY: gw_hines_config
  USE mo_vdiff_config,       ONLY: vdiff_config
  USE mo_nh_testcases,       ONLY: linit_tracer_fv,nh_test_name
  USE mo_ha_testcases,       ONLY: ctest_name

  USE mo_datetime,           ONLY: add_time, print_datetime_all

  IMPLICIT NONE

!  PRIVATE

  PUBLIC :: atm_crosscheck !, atmospheric_configuration

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  SUBROUTINE atm_crosscheck

    INTEGER :: jg
    INTEGER :: jt   ! tracer loop index
    INTEGER :: i_listlen
    REAL(wp):: cur_datetime_calsec, end_datetime_calsec, length_sec
    CHARACTER(len=*), PARAMETER :: routine =  'atm_crosscheck'

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()

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

    ! Length of this integration is limited by length of the restart cycle.
    nsteps = MIN(nsteps,INT(time_config%dt_restart/dtime))

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
    ! Testcases
    !--------------------------------------------------------------------
    IF(global_cell_type==3) THEN
      linit_tracer_fv  = .TRUE. ! like default
    ELSE
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

    !--------------------------------------------------------------------
    ! Shallow water
    !--------------------------------------------------------------------
    IF (iequations==ISHALLOW_WATER.AND.ha_dyn_config%lsi_3d) THEN
      CALL message( TRIM(routine), 'lsi_3d = .TRUE. not appicable to shallow water model')
    ENDIF

    IF ((iequations==ISHALLOW_WATER).AND.(nlev/=1)) &
    CALL finish(TRIM(routine),'Multiple vertical level specified for shallow water model')

    !--------------------------------------------------------------------
    ! Hydrostatic atm
    !--------------------------------------------------------------------
    IF (iequations==IHS_ATM_THETA) ha_dyn_config%ltheta_dyn = .TRUE.

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


        ! check radiation scheme in relation to chosen ozone

        IF (  atm_phy_nwp_config(jg)%inwp_radiation > 0 )  THEN

          SELECT CASE (irad_o3)
          CASE (0,6) ! ok
          CASE default
            CALL finish(TRIM(routine),'irad_o3 currently has to be 0 or 6.')
          END SELECT
        ENDIF

      ENDDO
    END IF



    !--------------------------------------------------------------------
    ! Tracers and diabatic forcing
    !--------------------------------------------------------------------
    SELECT CASE (iforcing)
    CASE(IECHAM,ILDF_ECHAM)

      IF (ntracer < 3) &
      CALL finish(TRIM(routine),'ECHAM physics needs at least 3 tracers')

    CASE(INWP)

      ! If iforcing=INWP is chosen, we do not want to set ntracer explicitly.
      ! Instead, ntracer and ntracer_static will be set automatically, 
      ! according to the selected radiation scheme.
      !IF (ntracer < 5) &
      !CALL finish(TRIM(routine),'NWP physics needs at least 5 tracers')

    END SELECT

    !
    ! Tracer indices need to be set before further checking
    !
    SELECT CASE(iforcing)
    CASE (IECHAM,ILDF_ECHAM)

      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqcond = iqi   !! index of last hydrometeor to ease summation over all of them
      iqt    = 4     !! starting index of non-water species 
      io3    = 5     !! O3
      ico2   = 6     !! CO2

    CASE (INWP)

      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqr    = 4     !! rain water
      iqs    = 5     !! snow
      iqcond = iqs   !! index of last hydrometeor to ease summation over all of them
      io3    = 6     !! O3
      ico2   = 7     !! CO2
      iqt    = 6     !! start index of other tracers than hydrometeors

    CASE default

      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqcond = iqi   !! index of last hydrometeor to ease summation over all of them
      iqt    = 4     !! starting index of non-water species
      io3    = 5     !! O3
      ico2   = 6     !! CO2

    END SELECT

    ntracer_static = 0

    IF (ltransport) THEN
    DO jg = 1,n_dom

      i_listlen = LEN_TRIM(advection_config(jg)%ctracer_list)

      SELECT CASE ( iforcing )
      CASE ( INWP )
      !...........................................................
      ! in NWP physics
      !...........................................................
        
        SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation)
        CASE (0)
          IF ( ntracer /= iqcond ) THEN
            ntracer = iqcond
            WRITE(message_text,'(a,i3)') 'Attention: for NWP physics, '//&
                                         'ntracer is set to',iqcond
            CALL message(TRIM(routine),message_text)
          ENDIF
          IF ( i_listlen /= iqcond ) THEN
            DO jt=1,ntracer
              WRITE(advection_config(jg)%ctracer_list(jt:jt),'(i1.1)')jt
            ENDDO
            WRITE(message_text,'(a)') &
              & 'Attention: according to physics, ctracer_list is set to ',&
              & advection_config(jg)%ctracer_list(1:ntracer)
            CALL message(TRIM(routine),message_text)
          ENDIF
        CASE (1)
          ntracer_static = 1
          IF ( ntracer /= iqcond ) THEN
            ntracer = iqcond
            WRITE(message_text,'(a,i3)') &
              &  'Attention: according to physics, ntracer is set to', iqcond   
            CALL message(TRIM(routine),message_text)
            WRITE(message_text,'(a)') &
              &  'In addition, there is one static tracer for O3'
            CALL message(TRIM(routine),message_text)
          ENDIF
          IF ( i_listlen /= ntracer ) THEN
            DO jt=1,ntracer
              WRITE(advection_config(jg)%ctracer_list(jt:jt),'(i1.1)')jt
            ENDDO
            WRITE(message_text,'(a)') &
              & 'Attention: according to physics, ctracer_list is set to ',&
              &   advection_config(jg)%ctracer_list(1:ntracer)
            CALL message(TRIM(routine),message_text)
          ENDIF
        CASE (2)
          SELECT CASE (irad_o3)
          CASE (0)
            IF ( ntracer /= iqcond  ) THEN
              ntracer = iqcond
              WRITE(message_text,'(a,i3)') &
                &  'Attention: according to physics, ntracer is set to', iqcond      
              CALL message(TRIM(routine),message_text)
            ENDIF
            IF ( i_listlen /= ntracer ) THEN
              DO jt=1,ntracer
                WRITE(advection_config(jg)%ctracer_list(jt:jt),'(i1.1)')jt
              ENDDO
              WRITE(message_text,'(a)') &
                & 'Attention: according to physics, ctracer_list is set to ', &
                &  advection_config(jg)%ctracer_list(1:ntracer)
              CALL message(TRIM(routine),message_text)
            ENDIF
          CASE (6)
            ntracer_static = 1
            IF ( ntracer /= iqcond  ) THEN
              ntracer = iqcond
              WRITE(message_text,'(a,i3)') &
                &  'Attention: according to physics, ntracer is set to', iqcond
              CALL message(TRIM(routine),message_text)           
              WRITE(message_text,'(a)') &
                &  'In addition, there is one static tracer for O3'
              CALL message(TRIM(routine),message_text)           
            ENDIF
            IF ( i_listlen /= ntracer ) THEN
              DO jt=1,ntracer
                WRITE(advection_config(jg)%ctracer_list(jt:jt),'(i1.1)')jt
              ENDDO
              WRITE(message_text,'(a)') &
                & 'Attention: according to physics with radiation and O3 ', &
                &  'ctracer_list is set to ', &
                &  advection_config(jg)%ctracer_list(1:ntracer)
              CALL message(TRIM(routine),message_text)
            ENDIF
          END SELECT
        END SELECT


        IF ( ( atm_phy_nwp_config(jg)%inwp_radiation > 0 )      &
          &  .AND. (irad_o3==0 .OR. irad_o3==6) )         THEN
          IF ( advection_config(jg)%ihadv_tracer(io3) /= 0 ) THEN
            advection_config(jg)%ihadv_tracer(io3) = 0
            WRITE(message_text,'(a,i1,a)') &
              & 'Attention: Since irad_o3 is set to ',irad_o3,', ihadv_tracer(io3) is set to 0.'
            CALL message(TRIM(routine),message_text)
          ENDIF
          IF ( advection_config(jg)%ivadv_tracer(io3) /= 0 ) THEN
            advection_config(jg)%ivadv_tracer(io3) = 0
            WRITE(message_text,'(a,i1,a)') &
              & 'Attention: Since irad_o3 is set to ',irad_o3,', ivadv_tracer(io3) is set to 0.'
            CALL message(TRIM(routine),message_text)
          ENDIF
        ENDIF


      CASE (inoforcing, iheldsuarez, iecham, ildf_dry, ildf_echam)
      !...........................................................
      ! Other types of adiabatic forcing
      !...........................................................
      
        IF ( i_listlen < ntracer .AND. i_listlen /= 0 ) THEN
          ntracer = i_listlen
          CALL message(TRIM(routine),'number of tracers is adjusted according to given list')
        END IF

        IF ((iforcing==IECHAM).AND.(echam_phy_config%lrad).AND. &
            (irad_o3 > 0) .AND. (io3 > ntracer) ) THEN

          CALL print_value('irad_o3' ,irad_o3)
          CALL print_value('io3    ' ,io3)
          CALL print_value('ntracer' ,ntracer)
          CALL finish(TRIM(routine), 'Not enough tracers for ECHAM physics with RRTM.')
        END IF
      
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
        ! 3rd order upwind horizontal advection scheme
        advection_config(jg)%ihadv_tracer(:) = iup3
        ! semi monotonous flux limiter
        advection_config(jg)%itype_hlimit(:) = ifluxl_sm
        WRITE(message_text,'(a)') 'Attention: on the hexagonal grid, ',              &
          &  'ihadv_tracer(:) = iup3 and itype_hlimit(:) = ifluxl_sm are the only ', &
          &  'possible settings. Please adjust your namelist settings accordingly'
        CALL message(TRIM(routine),message_text)
      END SELECT

      !----------------------------------------------
      ! Flux compuation methods - consistency check

      SELECT CASE (global_cell_type)
      CASE (3)
        IF ( ANY(advection_config(jg)%ihadv_tracer(1:ntracer) > 3))   THEN
          CALL finish( TRIM(routine),                                       &
               'incorrect settings for TRI-C grid ihadv_tracer. Must be 0,1,2, or 3 ')
        ENDIF
      CASE (6)
        IF (ANY(advection_config(jg)%ihadv_tracer(1:ntracer) == 2) .OR.     &
         &  ANY(advection_config(jg)%ihadv_tracer(1:ntracer) == 3))   THEN 
          CALL finish( TRIM(routine),                                       &
               'incorrect settings for HEX-C grid ihadv_tracer. Must be 0,1, or 4 ')
        ENDIF
      END SELECT


      !----------------------------------------------
      ! Limiter - consistency check

      IF (global_cell_type == 6) THEN
        IF ( ANY(advection_config(jg)%itype_hlimit(1:ntracer) == islopel_sm ) .OR.   &
          &  ANY(advection_config(jg)%itype_hlimit(1:ntracer) == islopel_m  ) .OR.   &
          &  ANY(advection_config(jg)%itype_hlimit(1:ntracer) == ifluxl_m   )) THEN
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

      CASE(2,4)
        CONTINUE

      CASE(3)
        IF (global_cell_type==3) CALL finish(TRIM(routine), &
        ' hdiff_order = 3 invalid for triangular model.')

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


  END  SUBROUTINE atm_crosscheck

!  SUBROUTINE atmospheric_configuration
!  INTEGER :: jg
!  CHARACTER(len=*), PARAMETER :: routine =  'atm_setup'
!  END SUBROUTINE atmospheric_configuration

END MODULE mo_atm_nml_crosscheck

