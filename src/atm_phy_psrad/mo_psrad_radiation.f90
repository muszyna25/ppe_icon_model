!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide interface to radiation routines. 
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the radiation code.  Mostly it organizes and calculates the 
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19): 
!!
!!         Hauke Schmidt, MPI-M, Hamburg (2009-12-18): Few modifications to
!!              allow specific solar irradiance for AMIP-type and preindustrial 
!!              simulations.
!!         Luis Kornblueh, MPI-M, Hamburg (2010-04-06): Never ever use write 
!!              directly 
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              Extracted public parameters into new module mo_radiation_parameters
!!              to avoid circular dependencies in submodels
!!                                      (2010-06-03):
!!              Added submodel calls, decl_sun_cur
!!         Robert Pincus, U. Colorado, while visiting MPI-M (2011-08-16) 
!!              Replaced underlying SW and LW schemes
!!         Dagmar Popke, MPI-M, Hamburg (2013-11-15):
!!              Implementation of RCE
!!         Sebastian Rast, MPI-M, Hamburg (2015-01-06):
!!              Adaption to ICON from echam-6.3 (revision 3735)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously contained in the ECHAM5 routines rad_int.f90, 
!!   radiation.f90 and prerad.f90.  Modifications were also made to provide
!!   a cleaner interface to the aerosol and cloud properties. Contributors to
!!   the code from which the present routines were derived include:  M. Jarraud,
!!   ECMWF (1983-06); M.A. Giorgetta, MPI-M (2002-05); U. Schulzweida,  MPI-M
!!   (2002-05); P. Stier MPI-M \& Caltech (2004-04, 2006-07), M. Thomas MPI-M 
!!   (2007-06); U. Schlese, MPI-M (2007-06); M. Esch, MPI-M (2007-06); S.J. 
!!   Lorenz, MPI-M (2007-11); T. Raddatz, MPI-M (2006-05); I. Kirchner.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
MODULE mo_psrad_radiation

  USE mo_kind,            ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_physical_constants,       ONLY: vtmpc1, rae,           &
       &                        amco2, amch4, amn2o, amo2, amd
!  USE mo_control,         ONLY: lcouple, lmidatm
!  USE mo_time_base,       ONLY: get_calendar_type, JULIAN
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_mpi,             ONLY: my_process_is_stdio
  USE mo_namelist,        ONLY: open_nml, position_nml, close_nml, POSITIONED
  USE mo_io_units,        ONLY: nnml, nnml_output
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist
!  USE mo_param_switches,  ONLY: lrad
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
!  USE mo_time_control,    ONLY: l_orbvsop87, get_orbit_times,                 &
!       &                        p_bcast_event, current_date, next_date,       &
!       &                        previous_date, radiation_date,                &
!       &                        prev_radiation_date,get_date_components,      &
!       &                        lresume, lstart, get_month_len
  USE mo_echam_convect_tables,  ONLY : prepare_ua_index_spline, lookup_ua_spline
  USE mo_datetime,        ONLY: t_datetime, rdaylen
  USE mo_get_utc_date_tr, ONLY: get_utc_date_tr
! amu0_x must now be taken from prm_field (has to be passed to the respect. routines
!  USE mo_geoloc,          ONLY: coslon_2d, &
!       &                        sinlon_2d, sinlat_2d, coslat_2d
  USE mo_model_domain,    ONLY: p_patch
! the present mo_orbit is different from the one in echam6
  USE mo_psrad_orbit,     ONLY: cecc, cobld, clonp, &
                              & orbit_kepler, orbit_vsop87, &
                              & get_orbit_times
  USE mo_psrad_radiation_parameters, ONLY: solar_parameters, &
                              & l_interp_rad_in_time,        &
                              & zepzen
! the present mo_bc_solar_irradiance is "old" and not the one used in echam-6.3.
!  USE mo_solar_irradiance,ONLY: get_solar_irradiance, set_solar_irradiance, &
!                                get_solar_irradiance_m, set_solar_irradiance_m
! cloud optics does not exist in icon
  USE mo_psrad_cloud_optics,    ONLY: setup_cloud_optics  
! greenhouse gases in mo_radiation_config, but probably in different "units" (ppm...)
!  USE mo_greenhouse_gases, ONLY: mmr_co2, mmr_ch4, mmr_n2o, ghg_cfcvmr
! ozone: read by mo_bc_ozone in icon, does not contain full functionality like here.
!  USE mo_o3clim,          ONLY: pre_o3clim_4, pre_o3clim_3, o3clim,            &
!       &                        read_o3clim_3
!  USE mo_o3_lwb,          ONLY: o3_lwb
! used for interactive CO2, can be switched off for the moment.
!  USE mo_submodel,        ONLY: lco2
! the following module is for diagnostic purposes only
!  USE mo_memory_cfdiag,   ONLY: locfdiag, &
!                              irlu, srsu, irld, srsd, irlucs, srsucs, irldcs, srsdcs
! namelist parameters. Reorganize reading
  USE mo_psrad_radiation_parameters, ONLY: nb_sw 
  USE mo_radiation_config,           ONLY: ih2o=>irad_h2o,     &
                                           ico2=>irad_co2,     &
                                           ich4=>irad_ch4,     &
                                           io3=>irad_o3,       &
                                           io2=>irad_o2,       &
                                           in2o=>irad_n2o,     &
                                           icfc11=>irad_cfc11, &
                                           icfc12=>irad_cfc12, &
                                           iaero=>irad_aero,   &
                                           mmr_co2,            &
                                           mmr_ch4,            &
                                           mmr_n2o,            &
                                           mmr_o2,             &
                                           vmr_cfc11,          &
                                           vmr_cfc12,          &
                                           nmonth,             &
                                           isolrad,            &
                                           ldiur,              &
                                           lyr_perp,           &
                                           yr_perp
!!$  USE mo_radiation_parameters, ONLY: ldiur, lradforcing,                          &
!!$                                     l_interp_rad_in_time, zepzen,                &
!!$                                     lyr_perp, yr_perp, nmonth, isolrad, nb_sw,   &
!!$                                     lw_spec_samp, sw_spec_samp,                  & 
!!$                                     lw_gpts_ts,   sw_gpts_ts,   rad_perm,        &
!!$                                     l_do_sep_clear_sky, i_overlap,               &
!!$                                     ih2o, ico2, ich4, io3, io2, in2o, icfc,      &
!!$                                     ighg, fco2, nmonth, iaero,                   &
!!$                                     co2vmr, ch4vmr, o2vmr, n2ovmr, cfcvmr,       &
!!$                                     co2mmr, ch4mmr, o2mmr, n2ommr,               &
!!$                                     ch4_v, n2o_v, cemiss, solc,                  &
!!$                                     psct, psctm, ssi_factor,                     &
!!$                                     flx_ratio_cur, flx_ratio_rad,                & 
!!$                                     decl_sun_cur,solar_parameters

! following module for diagnostic of radiative forcing only
!  USE mo_radiation_forcing,ONLY: prepare_forcing

  USE mo_rrtm_params,   ONLY : nbndsw
! new to icon
!  USE mo_psrad_srtm_setup,ONLY : ssi_default, ssi_preind, ssi_amip,           &
!                             & ssi_RCEdiurnOn, ssi_RCEdiurnOff
  USE mo_psrad_interface,ONLY : setup_psrad, &
!, psrad_interface, &
                                lw_strat, sw_strat
  USE mo_psrad_spec_sampling, ONLY : spec_sampling_strategy, &
                             & set_spec_sampling_lw, set_spec_sampling_sw, get_num_gpoints

  IMPLICIT NONE
  
  PRIVATE

    LOGICAL            :: lradforcing(2)=(/.FALSE.,.FALSE./)
    INTEGER            :: lw_gpts_ts=1,  &
                        & lw_spec_samp=1,&
                        & rad_perm=0,    &
                        & sw_gpts_ts=1,  &
                        & sw_spec_samp=1

  
  PUBLIC :: pre_psrad_radiation, setup_psrad_radiation, psrad_radiation

  CONTAINS

  SUBROUTINE pre_psrad_radiation(p_patch,datetime_radiation,ltrig_rad,amu0_x,rdayl_x)
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !

    TYPE(t_patch), INTENT(IN)        :: p_patch
    TYPE(t_datetime), INTENT(IN)     :: datetime_radiation !< date and time of radiative transfer calculation
    LOGICAL         , INTENT(IN)     :: ltrig_rad !< .true. if radiative transfer calculation has to be done at current time step
    REAL(wp), INTENT(OUT)            :: amu0_x(:,:), rdayl_x(:,:)

    LOGICAL  :: l_rad_call, l_write_solar
    INTEGER  :: icurrentyear, icurrentmonth, iprevmonth, i
    REAL(wp) :: rasc_sun, decl_sun, dist_sun, time_of_day, zrae
    REAL(wp) :: solcm, orbit_date, flx_ratio_cur
    LOGICAL  :: l_orbvsop87
    !
    ! 1.0 Compute orbital parameters for current time step
    ! --------------------------------
!!$    l_rad_call = .FALSE.
    CALL get_orbit_times(datetime_radiation, time_of_day, &
         &               orbit_date)

    IF (l_orbvsop87) THEN 
      CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
    ELSE
      CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
    END IF
!!$    decl_sun_cur = decl_sun       ! save for aerosol and chemistry submodels
    CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
!!$         &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
         &                p_patch,                         &
         &                flx_ratio_cur, amu0_x, rdayl_x)

!!$    IF (lrad) THEN
!!$
!!$      SELECT CASE (isolrad)
!!$      CASE (0)
!!$        solc = SUM(ssi_default)
!!$      CASE (1)
!!$        CALL get_solar_irradiance(current_date, next_date)
!!$        CALL set_solar_irradiance(solc)
!!$      CASE (2)
!!$        solc = SUM(ssi_preind)
!!$      CASE (3)
!!$        solc = SUM(ssi_amip)
!!$      CASE (4)
!!$        solc = SUM(ssi_RCEdiurnOn)
!!$      CASE (5)
!!$        solc = SUM(ssi_RCEdiurnOff)
!!$      CASE default
!!$        WRITE (message_text, '(a,i2,a)') &
!!$             'isolrad = ', isolrad, ' in radctl namelist is not supported'
!!$        CALL message('pre_radiation', message_text)
!!$      END SELECT
!!$      psct = flx_ratio_cur*solc
!!$
!!$    END IF ! lrad

    !
    ! 2.0 Prepare time dependent quantities for rad (on radiation timestep)
    ! --------------------------------
!!$    IF (phy_config%lrad .AND. ltrig_rad) THEN
!!$      l_rad_call = .TRUE.
!!$      CALL get_orbit_times(datetime_radiation, time_of_day , &
!!$           &               orbit_date)
!!$
!!$      IF ( l_orbvsop87 ) THEN 
!!$        CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
!!$      ELSE
!!$        CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
!!$      END IF
!!$      CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
!!$           &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
!!$           &                p_patch,                         &
!!$           &                flx_ratio_rad ,amu0m_x, rdaylm_x)
!!$      !
!!$      ! consider curvature of the atmosphere for high zenith angles
!!$      !
!!$      zrae = rae*(rae+2.0_wp)
!!$      amu0m_x(:,:)  = rae/(SQRT(amu0m_x(:,:)**2+zrae)-amu0m_x(:,:))
!!$      !
!!$      ! For the calculation of radiative transfer, a maximum zenith angle
!!$      ! of about 84 degrees is applied in order to avoid to much overshooting
!!$      ! when the extrapolation of the radiative fluxes from night time
!!$      ! regions to daytime regions is done for time steps at which no
!!$      ! radiation calculation is performed. This translates into cosines
!!$      ! of the zenith angle > 0.1.  This approach limits the calculation of the 
!!$      ! curvature effect above, and should be reconsidered in the future
!!$      ! 
!!$      amu0m_x(:,:) = MAX(amu0m_x(:,:),0.1_wp)
!!$      !
!!$      ! --- Prepare Ozone climatology
!!$      !
!!$      SELECT CASE (io3)
!!$      CASE (3) 
!!$        CALL pre_o3clim_3(nmonth)
!!$      CASE (4) 
!!$        CALL pre_o3clim_4
!!$      END SELECT
!!$      !
!!$
!!$      !++jsr&hs
!!$      ! 3.0 Prepare possibly time dependent total solar and spectral irradiance
!!$      ! --------------------------------
!!$      ! ATTENTION: 
!!$      ! This part requires some further work. Currently, a solar constant of
!!$      ! 1361.371 is used as default. This is the TSI averaged over the
!!$      ! years 1979 to 1988, and should be used for AMIP type runs. If lcouple is
!!$      ! true, a solar constant of 1360.875 is used, the average for the years 1844
!!$      ! to 1856. This should be used for a preindustrial control run.
!!$      ! The spectral distribution of this TSI is currently also prescribed for
!!$      ! these two cases depending on the lcouple switch.
!!$      ! For transient CMIP5 simulations, the available time
!!$      ! varying TSI and SSI has to be read in and used here.
!!$
!!$      SELECT CASE (isolrad)
!!$      CASE (0)
!!$        solcm = SUM(ssi_default)
!!$        ssi_factor = ssi_default
!!$      CASE (1)
!!$        CALL get_solar_irradiance_m(prev_radiation_date, radiation_date, nb_sw)
!!$        CALL set_solar_irradiance_m(solcm, ssi_factor, nb_sw)
!!$      CASE (2)
!!$        solcm = SUM(ssi_preind)
!!$        ssi_factor = ssi_preind
!!$      CASE (3)
!!$        solcm = SUM(ssi_amip)
!!$        ssi_factor = ssi_amip
!!$      CASE (4)
!!$        solcm = SUM(ssi_RCEdiurnOn)
!!$        ssi_factor = ssi_RCEdiurnOn
!!$      CASE (5)
!!$        solcm = SUM(ssi_RCEdiurnOff)
!!$        ssi_factor = ssi_RCEdiurnOff
!!$      CASE default
!!$        WRITE (message_text, '(a,i2,a)') &
!!$             'isolrad = ', isolrad, ' in radctl namelist is not supported'
!!$        CALL message('pre_radiation', message_text)
!!$      END SELECT
!!$      psctm = flx_ratio_rad*solcm
!!$      ssi_factor(:) = ssi_factor(:)/solcm
!!$
!!$      ! output of solar constant every month
!!$
!!$      CALL get_date_components(current_date, month=icurrentmonth, &
!!$           year=icurrentyear)
!!$      CALL get_date_components(previous_date, month=iprevmonth)
!!$      l_write_solar = icurrentmonth/=iprevmonth
!!$      IF (l_write_solar .OR. lresume .OR. lstart) THEN
!!$        CALL message('','')
!!$        WRITE (message_text,'(a,i0,a,i2.2,a,f6.1)') &
!!$             'Total solar constant [W/m^2] for ',      &
!!$             icurrentyear, '-',                        &
!!$             icurrentmonth, ' = ', solc
!!$        CALL message('',message_text)
!!$        CALL message('','')
!!$        DO i = 1, nb_sw
!!$          WRITE (message_text,'(a,i2,a,f7.5)') &
!!$               '   solar constant fraction: band ', i, &
!!$               ' = ', ssi_factor(i)
!!$          CALL message('',message_text)
!!$        END DO
!!$      END IF
!!$      !--jsr&hs
!!$
!!$    END IF ! lrad .AND. l_trigrad

  END SUBROUTINE pre_psrad_radiation

  SUBROUTINE setup_psrad_radiation(file_name)
!!$    USE mo_aero_kinne,       ONLY: su_aero_kinne
!!$    USE mo_aero_volc,        ONLY: su_aero_volc
!!$    USE mo_aero_volc_tab,    ONLY: su_aero_prop_ham, su_aero_prop_crow, &
!!$                                   read_aero_volc_tables
!!$    USE mo_solar_irradiance, ONLY: init_solar_irradiance

    CHARACTER(len=*), INTENT(IN)      :: file_name
    INTEGER :: istat, funit

    NAMELIST /psrad_nml/ lradforcing, &
                       & lw_gpts_ts,  &
                       & lw_spec_samp,&
                       & rad_perm,    &
                       & sw_gpts_ts,  &
                       & sw_spec_samp
    !
    ! 1.0 Read psrad_nml namelist 
    ! --------------------------------
    CALL open_nml(TRIM(file_name))
    CALL position_nml ('psrad_nml', status=istat)
    SELECT CASE (istat)
      CASE (POSITIONED) 
        READ (nnml, psrad_nml)
    END SELECT
    CALL close_nml
    ! store namelist for restart
    IF(my_process_is_stdio()) THEN
      funit = open_tmpfile()
      WRITE(funit,NML=psrad_nml)
      CALL store_and_close_namelist(funit, 'psrad_nml')
    END IF
    ! write the contents of the namelist to an ASCII file
    IF (my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=psrad_nml)
    END IF

    !
    ! 3.0 If radiation is active check NAMELIST variable conformance
    ! --------------------------------
    IF (phy_config%lrad) THEN

      CALL setup_psrad
      nb_sw = nbndsw
      !
      ! --- Spectral sampling strategy
      !
      lw_strat = set_spec_sampling_lw(lw_spec_samp, num_gpts_ts=lw_gpts_ts) 
      sw_strat = set_spec_sampling_sw(sw_spec_samp, num_gpts_ts=sw_gpts_ts) 
      WRITE (message_text, '("LW sampling strategy: ", i2, " using ", i3, " g-points per time step")') &
                 lw_spec_samp, get_num_gpoints(lw_strat)
      CALL message('',message_text)
      WRITE (message_text, '("SW sampling strategy: ", i2, " using ", i3, " g-points per time step")') &
                 sw_spec_samp, get_num_gpoints(sw_strat)
      CALL message('',message_text)

      !
      CALL message('','lrad = .TRUE.  --> Doing radiation radiation')
      !
      ! --- Check  H2O
      !
      SELECT CASE (ih2o)
      CASE(0)
        CALL message('','irad_h2o = 0 --> no H2O(gas,liquid,ice) in radiation')
      CASE(1)
        CALL message('','irad_h2o = 1 --> prognostic H2O(gas,liquid,ice)')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_h2o =', ih2o, ' in radiation_nml namelist is not supported'
        CALL message('', message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_h2o')
      END SELECT
      !
      ! --- Check  CO2
      ! 
      SELECT CASE (ico2)
!!$      CASE(0)
!!$        CALL message('','ico2 = 0 --> no CO2 in radiation')
!!$        co2mmr=co2vmr*amco2/amd   ! Necessary for use with lco2=.TRUE.
!!$      CASE(1)  
!!$        IF (lco2) THEN
!!$          WRITE (message_text, '(a,e16.8)') &
!!$               'irad_co2 = 1 --> Initial CO2 mass mixing ratio=', mmr_co2
!!$          CALL message('',message_text)
!!$          co2mmr = co2vmr*amco2/amd
!!$        ELSE
!!$          CALL finish('setup_psrad_radiation','irad_co2=1 (interactive CO2) not '// &
!!$               &      'a valid choice for lco2=.false.')
!!$        END IF
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_co2 = 2 --> CO2 mass mixing ratio=', mmr_co2
        CALL message('',message_text)
!!$        co2mmr = co2vmr*amco2/amd
!!$      CASE(4)
!!$        CALL message('','irad_co2 = 4 --> CO2 mass mixing ratio from scenario')
!!$        IF (ABS(fco2-1._wp) > EPSILON(1._wp)) THEN
!!$           WRITE (message_text, '(a,e16.8,a)') &
!!$                'fco2 = ', fco2, ' --> Factor for CO2 scenario'
!!$           CALL message('',message_text)
!!$        END IF
!!$        co2mmr = co2vmr*fco2*amco2/amd    ! This is only a dummy value for the first
!!$                                          ! initialization of co2m1 in the co2-module. 
!!$                                          ! co2m1 will be overwritten with the correct
!!$                                          ! values as soon as the ghg-data are
!!$                                          ! interpolated to the right date.
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_co2 = ', ico2, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_co2')
      END SELECT
!!$
!!$      IF (ico2 /= 4 .AND. ABS(fco2-1._wp) > EPSILON(1._wp)) THEN
!!$         WRITE (message_text, '(a,i2,a)') &
!!$              'ico2 = ', ico2, ' and fco2 != 1. --> Ignoring fco2'
!!$         CALL message('',message_text)
!!$      END IF
!!$      !
!!$      ! --- Check CH4
!!$      ! 
      SELECT CASE (ich4)
      CASE(0)
        CALL message('','irad_ch4 = 0 --> no CH4 in radiation')
      CASE(1)
        CALL message('','irad_ch4 = 1 --> transported CH4 is not yet implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_ch4')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_ch4 = 2 --> CH4 mass mixing ratio=', mmr_ch4
        CALL message('',message_text)
!!$        ch4mmr = ch4vmr*amch4/amd
      CASE(3)
        WRITE (message_text, '(a,e16.8)') &
             'irad_ch4 = 3 --> CH4 (trop) mass mixing ratio =', mmr_ch4
        CALL message('',message_text)
!!$        ch4mmr = ch4vmr*amch4/amd
!!$      CASE(4)
!!$        CALL message('','irad_ch4 = 4 --> CH4 volume mixing ratio from scenario')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_ch4 =', ich4, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_ch4')
      END SELECT
      !
      ! --- Check O3
      ! 
      SELECT CASE (io3)
      CASE(0)
        CALL message('','irad_o3  = 0 --> no O3 in radiation')
      CASE(1)
        CALL message('','irad_o3  = 1 --> transported O3 is not yet implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_o3')
      CASE(2)
        CALL message('','irad_o3  = 2 --> spectral O3 climatology (ECHAM4), no implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_o3')
      CASE(8)
        CALL message('','irad_o3  = 8 --> gridpoint O3 climatology from NetCDF file, not implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_o3')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_o3  =', io3, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_o3')
      END SELECT
      !
      ! --- Check N2O
      ! 
      SELECT CASE (in2o)
      CASE(0)
        CALL message('','irad_n2o = 0 --> no N2O in radiation')
      CASE(1)
        CALL message('','irad_n2o = 1 --> transported N2O is not yet implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_n2o')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_n2o = 2 --> N2O mass mixing ratio=', mmr_n2o
        CALL message('',message_text)
!!$        n2ommr = n2ovmr*amn2o/amd
      CASE(3)
        WRITE (message_text, '(a,e16.8)') &
             'irad_n2o = 3 --> N2O (trop) mass mixing ratio=', mmr_n2o
        CALL message('',message_text)
!!$        n2ommr = n2ovmr*amn2o/amd
!!$      CASE(4)
!!$        CALL message('','irad_n2o = 4 --> N2O volume mixing ratio from scenario')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_n2o =',in2o,' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated in2o')
      END SELECT
      !
      ! --- Check CFCs
      ! 
      SELECT CASE (icfc11)
      CASE(0)
        CALL message('','irad_cfc11 = 0 --> no CFCs in radiation')
      CASE(1)
        CALL message('','irad_cfc11 = 1 --> transported CFCs not yet implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_cfc11')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_cfc11 = 2 --> CFC11    volume mixing ratio=', vmr_cfc11
        CALL message('',message_text)
!!$      CASE(4)
!!$        CALL message('','irad_cfc11 = 4 --> CFC11 volume mixing ratio from scenario')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_cfc11=', icfc11, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_cfc11')
      END SELECT

      SELECT CASE (icfc12)
      CASE(0)
        CALL message('','irad_cfc12 = 0 --> no CFCs in radiation')
      CASE(1)
        CALL message('','irad_cfc12 = 1 --> transported CFCs not yet implemented')
        CALL finish('setup_psrad_radiation','Run terminated irad_cfc12')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_cfc12 = 2 --> CFC12    volume mixing ratio=', vmr_cfc12
        CALL message('',message_text)
!!$      CASE(4)
!!$        CALL message('','irad_cfc12 = 4 --> CFC12 volume mixing ratio from scenario')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_cfc12=', icfc12, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_cfc12')
      END SELECT
!!$      !
!!$      ! --- Check Scenario
!!$      ! 
!!$      SELECT CASE (ighg)
!!$      CASE(0)
!!$        CALL message('','ighg = 0 --> no scenario, fixed greenhouse gases and/or cfc')
!!$      CASE(1)
!!$        CALL message('','ighg = 1 --> greenhouse gases from scenario, check setting of switches')
!!$      END SELECT
!!$      !
!!$      ! --- Check O2
!!$      ! 
!!$      SELECT CASE (io2)
!!$      CASE(0)
!!$        CALL message('','io2  = 0 --> no O2  in radiation')
!!$      CASE(2)
!!$        WRITE (message_text, '(a,e16.8)') &
!!$             'io2  = 2 --> O2    volume mixing ratio=', o2vmr
!!$        CALL message('',message_text)
!!$        o2mmr = o2vmr*amo2/amd
!!$      CASE default
!!$        WRITE (message_text, '(a,i2,a)') &
!!$             'io2 =', io2, ' in radctl namelist is not supported'
!!$        CALL message('',message_text)
!!$        CALL finish('setup_psrad_radiation','Run terminated io2')
!!$      END SELECT
!!$      !
!!$      ! --- Check aerosol
!!$      ! 
      SELECT CASE (iaero)
      CASE(0)
        CALL message('','irad_aero= 0 --> no aerosol in radiation')
!!$      CASE(1)
!!$        CALL message('','iaero= 1 --> prognostic aerosol (sub model)')
!!$      CASE(3)
!!$        CALL message('','iaero= 3 --> Kinne climatology')
!!$        CALL su_aero_kinne(nb_sw)
!!$      CASE(5)
!!$        CALL message('','iaero= 5 --> Kinne climatology + Stenchikov volcanic aerosol')
!!$        CALL su_aero_kinne(nb_sw)
!!$        CALL su_aero_volc(nb_sw)
!!$      CASE(6)
!!$        CALL message('','iaero= 6 --> Kinne climatology + Stenchikov volcanic aerosols + HAM volcanic aerosol')
!!$        CALL su_aero_kinne(nb_sw)
!!$        CALL su_aero_volc(nb_sw)
!!$        CALL su_aero_prop_ham
!!$        CALL read_aero_volc_tables
!!$      CASE(7)
!!$        CALL message('','iaero= 7 --> Kinne climatology + Crowley volcanic aerosol')
!!$        CALL su_aero_kinne(nb_sw)
!!$        CALL su_aero_prop_crow
!!$        CALL read_aero_volc_tables
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_aero=', iaero, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_aero')
      END SELECT
      !
      ! --- Check annual cycle
      ! 
      SELECT CASE (nmonth)
      CASE(0)
        CALL message('','nmonth=0 --> annual cycle on')
      CASE(1:12)
        WRITE (message_text, '(a,i2.2,a)') &
             'nmonth = ', nmonth, ' --> perpetual month'
        CALL message('',message_text)
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'nmonth=', nmonth, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated nmonth')
      END SELECT
      !
      ! --- Check Shortwave Model
      ! 
      CALL message('','  --> USE AER RRTM Shortwave Model')
      !
      ! --- Check Longwave Model
      ! 
      CALL message('','  --> USE New (V4) LRTM Model')
      !
      ! --- Check solar constant
      !
      SELECT CASE (isolrad)
      CASE (0) 
        CALL message('','isolrad = 0 --> standard rrtm solar constant')
      CASE (1) 
        CALL message('','isolrad = 1 --> time dependent spectrally resolved solar constant read from file')
!!$        CALL init_solar_irradiance(nb_sw)
      CASE (2) 
        CALL message('','isolrad = 2 --> preindustrial solar constant')
      CASE (3) 
        CALL message('','isolrad = 3 --> solar constant for amip runs')
!!$      CASE (4)
!!$        CALL message('','isolrad = 4 --> solar constant for rad.-convective eq. runs with diurnal cycle ON')
!!$      CASE (5)
!!$        CALL message('','isolrad = 5 --> solar constant for rad.-convective eq. runs with diurnal cycle OFF')
      CASE default 
        WRITE (message_text, '(a,i3,a)') &
             'Run terminated isolrad = ', isolrad, ' not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation', message_text)
      END SELECT
      !
      ! --- Check diurnal cycle
      ! 
      IF (ldiur) THEN
        CALL message('','ldiur =.TRUE.  --> diurnal cycle on')
      ELSE
        CALL message('','ldiur =.FALSE. --> diurnal cycle off')
      ENDIF
!!$      !
!!$      ! --- Check for diagnosis of instantaneous aerosol radiative forcing
!!$      ! 
!!$      CALL message('','instantaneous forcing diagnostic:')
!!$      WRITE (message_text,'(a16,L3,a18,L3)')       &
!!$           ' solar radiation: ',   lradforcing(1), &
!!$           ' thermal radiation: ', lradforcing(2)
!!$      CALL message('',message_text)
      !
      ! --- Check perpetual orbit
      ! 
      IF (yr_perp.NE.-99999)  lyr_perp = .TRUE.
!!$      CALL p_bcast (lyr_perp, p_io)
!!$
      IF (lyr_perp) THEN
!!$        IF (l_orbvsop87) THEN
          WRITE (message_text, '(a,i0,a)') &
               'yr_perp=', yr_perp, ' --> perpetual year for orbit'
          CALL message('',message_text)
!!$        ELSE
!!$          WRITE (message_text, '(a,i0,a,l1,a)') &
!!$               'yr_perp = ', yr_perp, ' l_orbvsop87 = ',l_orbvsop87,' not allowed!'
!!$          CALL message('',message_text)
!!$          CALL finish('setup_psrad_radiation', &
!!$               ' yr_perp.ne.-99999 cannot run  PCMDI-orbit (l_orbvsop87=.F.).')
!!$        END IF
      END IF
      !
      ! 4.0 Initialization for radiation
      ! -------------------------------
      !
      ! --- resolution/run dependent cloud optical parameters (tuning)
      !
      CALL setup_cloud_optics
      !
!!$      ! --- Ozone climatology
!!$      ! 
!!$      IF (io3==3) CALL read_o3clim_3
!!$      !
    ELSE
      CALL message('','phy_config%lrad = .FALSE. --> no radiation')
!!$      co2mmr = co2vmr*amco2/amd
    ENDIF
  END SUBROUTINE setup_psrad_radiation

  SUBROUTINE psrad_radiation ( &
    & jg         ,&!< in  domain index
    & jb         ,&!< in  block index
    & kproma     ,&!< in  end index for loop over block
    & kbdim      ,&!< in  dimension of block over cells
    & klev       ,&!< in  number of full levels = number of layers
    & klevp1     ,&!< in  number of half levels = number of layer interfaces
    & ktrac      ,&!< in number of non-water tracers
    & ktype      ,&!< in  type of convection
    & loland     ,&!< in  land-sea mask. (1. = land, 0. = sea/lakes)
    & loglac     ,&!< in  fraction of land covered by glaciers
    & pcos_mu0   ,&!< in  cosine of solar zenith angle
    & alb_vis_dir,&!< in  surface albedo for visible range, direct
    & alb_nir_dir,&!< in  surface albedo for near IR range, direct
    & alb_vis_dif,&!< in  surface albedo for visible range, diffuse
    & alb_nir_dif,&!< in  surface albedo for near IR range, diffuse
    & tk_sfc     ,&!< in  grid box mean surface temperature
    & pp_hl      ,&!< in  pressure at half levels at t-dt [Pa]
    & pp_fl      ,&!< in  pressure at full levels at t-dt [Pa]
    & tk_fl      ,&!< in  tk_fl  = temperature at full level at t-dt
    & qm_vap     ,&!< in  qm_vap = water vapor mass mixing ratio at t-dt
    & qm_liq     ,&!< in  qm_liq = cloud water mass mixing ratio at t-dt
    & qm_ice     ,&!< in  qm_ice = cloud ice mass mixing ratio at t-dt
    & geom1      ,&!< in  pgeom1 = geopotential above ground at t-dt [m2/s2]
    & cdnc       ,&!< in  cloud droplet number concentration
    & cld_frc    ,&!< in  cloud fraction
    & pxtm1       &!< tracer concentration
    &              )
    INTEGER, INTENT(in)  :: &
    & jg,             & !< domain index
    & jb,             & !< block index
    & kproma,         & !< end   index for loop over block
    & kbdim,          & !< dimension of block over cells
    & klev,           & !< number of full levels = number of layers
    & klevp1,         & !< number of half levels = number of layer interfaces 
    & ktrac,          & !< number of non-water tracers
    & ktype(kbdim)      !< convection type

    LOGICAL, INTENT(IN)  :: &
    & loland(kbdim),     & !< land mask
    & loglac(kbdim)        !< glacier mask

    REAL(wp), INTENT(IN) :: &
    & pcos_mu0(kbdim),   & !< cosine of solar zenith angle
    & alb_vis_dir(kbdim),& !< surface albedo for visible range and direct light
    & alb_nir_dir(kbdim),& !< surface albedo for NIR range and direct light
    & alb_vis_dif(kbdim),& !< surface albedo for visible range and diffuse light
    & alb_nir_dif(kbdim),& !< surface albedo for NIR range and diffuse light
    & tk_sfc(kbdim),     & !< Surface temperature
    & pp_hl(kbdim,klevp1),& !< pressure at half levels [Pa]
    & pp_fl(kbdim,klev),  & !< Pressure at full levels [Pa]
    & tk_fl(kbdim,klev),  & !< Temperature on full levels [K]
    & qm_vap(kbdim,klev), & !< Water vapor mixing ratio
    & qm_liq(kbdim,klev), & !< Liquid water mixing ratio
    & qm_ice(kbdim,klev), & !< Ice water mixing ratio
    & geom1(kbdim,klev),  & !< Geopotential height at t-dt
    & cdnc(kbdim,klev),   & !< Cloud drop number concentration
    & cld_frc(kbdim,klev),& !< Cloud fraction
    & pxtm1(kbdim,klev,ktrac)!< non-water tracers
    INTEGER              :: jk, jl, idx(kbdim)
    REAL(wp)             ::         &
    & cos_mu0(kbdim),               &
    & pp_sfc(kbdim),                &
    & ppd_hl(kbdim,klev),           &
    & tk_hl(kbdim,klevp1),          &
    & xq_vap(kbdim,klev),           &
    & za(kbdim),                    & !< Spline interpolation arrays for qsat
    & ua(kbdim),                    &
    & xq_sat(kbdim,klev)              !< Saturation mixing ratio 

    !
    ! 1.0 calculate variable input parameters (location and state variables)
    ! --------------------------------
    ! 
    ! --- solar zenith angle
    !
    ! Get the local cosine of the solar zenith angle, setting a minimum positive definite 
    ! value for smooth interpretation in time (if desired(
    !
    cos_mu0(1:kproma) = pcos_mu0(1:kproma)
    IF (l_interp_rad_in_time) cos_mu0(1:kproma) = MAX(cos_mu0(1:kproma), zepzen) 
    !
    ! --- Pressure (surface and distance between half levels)
    !
    pp_sfc(1:kproma)   = pp_hl(1:kproma,klevp1)
    ppd_hl(1:kproma,:) = pp_hl(1:kproma,2:klev+1)-pp_hl(1:kproma,1:klev)
    !
    ! --- temperature at half levels
    !
    DO jk=2,klev
      DO jl = 1, kproma
        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
             & - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
             & - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
      END DO
    END DO
    DO jl = 1, kproma
      tk_hl(jl,klevp1) = tk_sfc(jl)
      tk_hl(jl,1)      = tk_fl(jl,1)-pp_fl(jl,1)*(tk_fl(jl,1) - tk_hl(jl,2))  &
           &             / (pp_fl(jl,1)-pp_hl(jl,2))
    END DO
    !
    ! --- phases of water substance
    !
    xq_vap(1:kproma,:) = MAX(qm_vap(1:kproma,:),EPSILON(1.0_wp))
    DO jk = 1, klev
      CALL prepare_ua_index_spline('radiation', kproma, tk_fl(:,jk), idx(:),za(:))
      CALL lookup_ua_spline(kproma, idx(:), za(:), ua(:))
      xq_sat(1:kproma,jk) = ua(1:kproma)/pp_fl(1:kproma,jk)
      xq_sat(1:kproma,jk) = MIN(xq_sat(1:kproma,jk),0.5_wp)
      xq_sat(1:kproma,jk) = xq_sat(1:kproma,jk)/(1.0_wp-vtmpc1                &
           &                * xq_sat(1:kproma,jk))
      xq_sat(1:kproma,jk) = MAX(2.0_wp*EPSILON(1.0_wp),xq_sat(1:kproma,jk))
    END DO



  END SUBROUTINE psrad_radiation

  END MODULE mo_psrad_radiation
