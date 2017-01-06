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

  USE mo_kind,            ONLY: wp, i8
  USE mo_model_domain,    ONLY: t_patch
  USE mo_physical_constants,       ONLY: rae
  USE mo_exception,       ONLY: finish, message, message_text, print_value
  USE mo_mpi,             ONLY: my_process_is_stdio
  USE mo_namelist,        ONLY: open_nml, position_nml, close_nml, POSITIONED
  USE mo_io_units,        ONLY: nnml, nnml_output
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist
  USE mo_impl_constants,      ONLY: io3_clim, io3_ape, io3_amip
  USE mo_ext_data_types,      ONLY: t_external_atmos_td
  USE mo_ext_data_state,      ONLY: ext_data, nlev_o3
  USE mo_bc_ozone,            ONLY: o3_plev, nplev_o3, plev_full_o3, plev_half_o3
  USE mo_o3_util,             ONLY: o3_pl2ml, o3_timeint
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
!  USE mo_time_control,    ONLY: l_orbvsop87, get_orbit_times,                 &
!       &                        p_bcast_event, current_date, next_date,       &
!       &                        previous_date, radiation_date,                &
!       &                        prev_radiation_date,get_date_components,      &
!       &                        lresume, lstart, get_month_len
  USE mo_psrad_orbit,     ONLY: orbit_kepler, orbit_vsop87, &
                              & get_orbit_times
  USE mo_psrad_orbit_nml, ONLY: read_psrad_orbit_namelist
  USE mo_psrad_radiation_parameters, ONLY: solar_parameters, &
                              & l_interp_rad_in_time,        &
                              & zepzen
! the present mo_bc_solar_irradiance is "old" and not the one used in echam-6.3.
!  USE mo_solar_irradiance,ONLY: get_solar_irradiance, set_solar_irradiance, &
!                                get_solar_irradiance_m, set_solar_irradiance_m
! cloud optics does not exist in icon
  USE mo_psrad_cloud_optics,    ONLY: setup_cloud_optics  
  USE mo_bc_greenhouse_gases,   ONLY: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcvmr
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqt, ico2, ntracer
! ozone: read by mo_bc_ozone in icon, does not contain full functionality like here.
!  USE mo_o3clim,          ONLY: pre_o3clim_4, pre_o3clim_3, o3clim,            &
!       &                        read_o3clim_3
!  USE mo_o3_lwb,          ONLY: o3_lwb
! used for interactive CO2, can be switched off for the moment.
  USE mo_radiation_config,           ONLY: irad_h2o,           &
                                           irad_co2,           &
                                           irad_ch4,           &
                                           irad_o3,            &
                                           irad_o2,            &
                                           irad_n2o,           &
                                           irad_cfc11,         &
                                           irad_cfc12,         &
                                           irad_aero,          &
                                           ighg,               &
                                           vmr_co2, mmr_co2,   &
                                           vmr_ch4, mmr_ch4,   &
                                           vmr_n2o, mmr_n2o,   &
                                           vmr_o2,  mmr_o2,    &
                                           vmr_cfc11,          &
                                           vmr_cfc12,          &
                                           fh2o, fco2, fch4,   &
                                           fn2o, fo3, fo2,     &
                                           fcfc,               &
                                           ch4_v=>vpp_ch4,     &
                                           n2o_v=>vpp_n2o,     &
                                           vmr_o2,             &
                                           nmonth,             &
                                           isolrad,            &
                                           ldiur,              &
                                           lyr_perp,           &
                                           yr_perp,            &
                                           lradforcing,        &
                                           tsi,                &
                                           tsi_radt,           &
                                           ssi_radt
  USE mo_psrad_radiation_parameters, ONLY: nb_sw,              & 
                                     irad_aero_forcing,        &
                                     l_interp_rad_in_time,     &
                                     zepzen,                   &
                                     lw_spec_samp,             &
                                     sw_spec_samp,             &
                                     lw_gpts_ts,               &
                                     sw_gpts_ts,               &
                                     rad_perm,                 &
                                     cemiss,                   &
                                     solc,                     &
                                     psct,                     &
                                     psctm,                    &
                                     ssi_factor,               &
                                     flx_ratio_cur,            &
                                     flx_ratio_rad,            &
                                     solar_parameters

! following module for diagnostic of radiative forcing only
  USE mo_psrad_radiation_forcing,ONLY: prepare_psrad_radiation_forcing

  USE mo_rrtm_params,   ONLY : nbndsw
! new to icon
  USE mo_psrad_srtm_setup,ONLY : ssi_default, ssi_preind, ssi_amip,           &
                             & ssi_RCEdiurnOn, ssi_RCEdiurnOff
  USE mo_psrad_interface,ONLY : setup_psrad, psrad_interface, &
                                lw_strat, sw_strat
  USE mo_psrad_spec_sampling, ONLY : set_spec_sampling_lw, set_spec_sampling_sw, get_num_gpoints
  USE mo_psrad_orbit_config,  ONLY : psrad_orbit_config

  USE mtime, ONLY: datetime
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: pre_psrad_radiation, setup_psrad_radiation, psrad_radiation

  CONTAINS

  SUBROUTINE pre_psrad_radiation( p_patch,          datetime_radiation, &
                                & current_datetime, ltrig_rad,          &
                                & amu0_x,           rdayl_x,            &
                                & amu0m_x,          rdaylm_x            )
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !

    TYPE(t_patch),           INTENT(in) :: p_patch
    TYPE(datetime), POINTER, INTENT(in) :: datetime_radiation, & !< date and time of radiative transfer calculation
         &                                 current_datetime       !< current time step
    LOGICAL,                 INTENT(in) :: ltrig_rad !< .true. if radiative transfer calculation has to be done at current time step
    REAL(wp),                INTENT(out) :: amu0_x(:,:), rdayl_x(:,:), &
         &                                  amu0m_x(:,:), rdaylm_x(:,:)

    LOGICAL  :: l_write_solar
    INTEGER(i8) :: icurrentyear
    INTEGER  :: icurrentmonth, i
    INTEGER, SAVE :: iprevmonth=-9999
    REAL(wp) :: rasc_sun, decl_sun, dist_sun, time_of_day, zrae
    REAL(wp) :: orbit_date
    REAL(wp) :: solcm
    LOGICAL  :: l_orbvsop87, l_sph_symm_irr

    l_orbvsop87 = psrad_orbit_config%l_orbvsop87
    l_sph_symm_irr = psrad_orbit_config%l_sph_symm_irr

    !
    ! 1.0 Compute orbital parameters for current time step
    ! --------------------------------
    CALL get_orbit_times(current_datetime, time_of_day, orbit_date)

    IF (l_orbvsop87) THEN 
      CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
    ELSE
      CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
    END IF

!!$    decl_sun_cur = decl_sun       ! save for aerosol and chemistry submodels
    CALL solar_parameters(decl_sun,       dist_sun,         time_of_day,       &
         &                ldiur,          l_sph_symm_irr,   p_patch,           &
         &                flx_ratio_cur,  amu0_x,           rdayl_x            )

    IF (phy_config%lrad) THEN
!!$
      SELECT CASE (isolrad)
      CASE (0)
        solc = SUM(ssi_default)
      CASE (1)
!!$        CALL get_solar_irradiance(current_date, next_date)
!!$        CALL set_solar_irradiance(solc)
        solc = tsi
        continue ! solar irradiance was read in echam_phy_bcs_global
      CASE (2)
        solc = SUM(ssi_preind)
      CASE (3)
        solc = SUM(ssi_amip)
      CASE (4)
        solc = SUM(ssi_RCEdiurnOn)
      CASE (5)
        solc = SUM(ssi_RCEdiurnOff)
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'isolrad = ',isolrad, ' in radiation_nml namelist is not supported'
        CALL finish('pre_psrad_radiation', message_text)
      END SELECT
      psct = flx_ratio_cur*solc

    END IF ! lrad

    !
    ! 2.0 Prepare time dependent quantities for rad (on radiation timestep)
    ! --------------------------------
    IF (phy_config%lrad .AND. ltrig_rad) THEN

      CALL get_orbit_times(datetime_radiation, time_of_day, orbit_date)

      IF ( l_orbvsop87 ) THEN 
        CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
      ELSE
        CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
      END IF
      CALL solar_parameters(decl_sun,        dist_sun,          time_of_day,        &
           &                ldiur,           l_sph_symm_irr,    p_patch,            &
           &                flx_ratio_rad,   amu0m_x,           rdaylm_x            )
      !
      ! consider curvature of the atmosphere for high zenith angles
      !
      zrae = rae*(rae+2.0_wp)
      amu0m_x(:,:)  = rae/(SQRT(amu0m_x(:,:)**2+zrae)-amu0m_x(:,:))
      !
      ! For the calculation of radiative transfer, a maximum zenith angle
      ! of about 84 degrees is applied in order to avoid to much overshooting
      ! when the extrapolation of the radiative fluxes from night time
      ! regions to daytime regions is done for time steps at which no
      ! radiation calculation is performed. This translates into cosines
      ! of the zenith angle > 0.1.  This approach limits the calculation of the 
      ! curvature effect above, and should be reconsidered in the future
      ! 
      amu0m_x(:,:) = MAX(amu0m_x(:,:),0.1_wp)
!!$      !
!!$      ! --- Prepare Ozone climatology
!!$      !
!!$      SELECT CASE (irad_o3)
!!$      CASE (3) 
!!$        CALL pre_o3clim_3(nmonth)
!!$      CASE (4) 
!!$        CALL pre_o3clim_4
!!$      END SELECT
!!$      !
!!$
      !++jsr&hs
      ! 3.0 Prepare possibly time dependent total solar and spectral irradiance
      ! --------------------------------
      ! ATTENTION: 
      ! This part requires some further work. Currently, a solar constant of
      ! 1361.371 is used as default. This is the TSI averaged over the
      ! years 1979 to 1988, and should be used for AMIP type runs. If lcouple is
      ! true, a solar constant of 1360.875 is used, the average for the years 1844
      ! to 1856. This should be used for a preindustrial control run.
      ! The spectral distribution of this TSI is currently also prescribed for
      ! these two cases depending on the lcouple switch.
      ! For transient CMIP5 simulations, the available time
      ! varying TSI and SSI has to be read in and used here.

      SELECT CASE (isolrad)
      CASE (0)
        solcm = SUM(ssi_default)
        ssi_factor = ssi_default
      CASE (1)
!!$        CALL get_solar_irradiance_m(prev_radiation_date, radiation_date, nb_sw)
!!$        CALL set_solar_irradiance_m(solcm, ssi_factor, nb_sw)
        solcm=tsi_radt
        ssi_factor=ssi_radt
        continue ! solar irradiance was read in echam_phy_bcs_global
      CASE (2)
        solcm = SUM(ssi_preind)
        ssi_factor = ssi_preind
      CASE (3)
        solcm = SUM(ssi_amip)
        ssi_factor = ssi_amip
      CASE (4)
        solcm = SUM(ssi_RCEdiurnOn)
        ssi_factor = ssi_RCEdiurnOn
      CASE (5)
        solcm = SUM(ssi_RCEdiurnOff)
        ssi_factor = ssi_RCEdiurnOff
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'isolrad = ', isolrad, ' in radctl namelist is not supported'
        CALL message('pre_radiation', message_text)
      END SELECT
      psctm = flx_ratio_rad*solcm
      ssi_factor(:) = ssi_factor(:)/solcm

      ! output of solar constant every month

!!$      CALL get_date_components(current_date, month=icurrentmonth, &
!!$           year=icurrentyear)
!!$      CALL get_date_components(previous_date, month=iprevmonth)

      icurrentmonth = datetime_radiation%date%month
      icurrentyear = datetime_radiation%date%year
      l_write_solar = icurrentmonth/=iprevmonth
      iprevmonth = icurrentmonth
      IF (l_write_solar) THEN
        CALL message('','')
        WRITE (message_text,'(a,i0,a,i2.2,a,f6.1)') &
             'Total solar constant [W/m^2] for ',      &
             icurrentyear, '-',                        &
             icurrentmonth, ' = ', solcm
        CALL message('',message_text)
        CALL message('','')
        DO i = 1, nb_sw
          WRITE (message_text,'(a,i2,a,f7.5)') &
               '   solar constant fraction: band ', i, &
               ' = ', ssi_factor(i)
          CALL message('',message_text)
        END DO
      END IF
      !--jsr&hs
    ELSE
      amu0m_x(:,:) = 0.0_wp
      rdaylm_x(:,:) = 0.0_wp
    END IF ! lrad .AND. l_trigrad

  END SUBROUTINE pre_psrad_radiation

  SUBROUTINE setup_psrad_radiation(file_name)
!!$    USE mo_aero_kinne,       ONLY: su_aero_kinne
!!$    USE mo_aero_volc,        ONLY: su_aero_volc
!!$    USE mo_aero_volc_tab,    ONLY: su_aero_prop_ham, su_aero_prop_crow, &
!!$                                   read_aero_volc_tables
!!$    USE mo_solar_irradiance, ONLY: init_solar_irradiance

    CHARACTER(len=*), INTENT(IN)      :: file_name
    INTEGER :: istat, funit

    NAMELIST /psrad_nml/ lradforcing,       & ! switch for short and longwave
     ! radiative forcing calculation by double call to radiation (default:
     ! (/.FALSE.,.FALSE./)) 
                       & irad_aero_forcing, & ! key number of aerosols 
     ! in reference radiation computation for radiative forcing calculation
     ! by double call to radiation
                       & lw_gpts_ts,        &
                       & lw_spec_samp,      &
                       & rad_perm,          &
                       & sw_gpts_ts,        &
                       & sw_spec_samp

    ! 0.9 Read psrad_orbit namelist
    CALL read_psrad_orbit_namelist(file_name)
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

      CALL message('','')
      CALL message('','PSrad setup')
      CALL message('','===========')
      CALL message('','- New (V4) LRTM Model')
      CALL message('','- AER RRTM Shortwave Model')
      CALL message('','')
      CALL print_value('radiation time step in [s]',phy_config%dt_rad)
      CALL message('','')
      !
      CALL setup_psrad
      nb_sw = nbndsw
      !
      ! --- Spectral sampling strategy
      !
      lw_strat = set_spec_sampling_lw(lw_spec_samp, num_gpts_ts=lw_gpts_ts) 
      sw_strat = set_spec_sampling_sw(sw_spec_samp, num_gpts_ts=sw_gpts_ts) 
      WRITE (message_text, '("LW sampling strategy", i2, ", using ", i3, " g-points per rad. time step")') &
                 lw_spec_samp, get_num_gpoints(lw_strat)
      CALL message('',message_text)
      WRITE (message_text, '("SW sampling strategy", i2, ", using ", i3, " g-points per rad. time step")') &
                 sw_spec_samp, get_num_gpoints(sw_strat)
      CALL message('',message_text)


      CALL message('','')
      CALL message('','Sources of volume/mass mixing ratios used in radiation')
      CALL message('','------------------------------------------------------')
      !
      ! --- Check  H2O
      !
      SELECT CASE (irad_h2o)
      CASE(0)
        CALL message('','irad_h2o   = 0 --> no H2O(gas,liq,ice) in radiation')
      CASE(1)
        CALL message('','irad_h2o   = 1 --> H2O   (gas,liq,ice) mass mixing ratios from tracer fields')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_h2o   =', irad_h2o, ' in radiation_nml namelist is not supported'
        CALL message('', message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_h2o')
      END SELECT
      !
      ! --- Check  CO2
      ! 
      SELECT CASE (irad_co2)
      CASE(0)
        CALL message('','irad_co2   = 0 --> no CO2 in radiation')
      CASE(1)
        IF ( iqt <= ico2 .AND. ico2 <= ntracer) THEN
          CALL message('','irad_co2   = 1 --> CO2   mass mixing ratio from tracer field')
        ELSE
          CALL finish('setup_psrad_radiation','irad_co2 = 1 (CO2 tracer in radiation) is not '// &
               &      'a valid choice because no CO2 tracer is available')
        END IF
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_co2   = 2 --> CO2   volume mixing ratio from radiation_nml namelist =', vmr_co2
        CALL message('',message_text)
      CASE(4)
        CALL message('','irad_co2   = 4 --> CO2   volume mixing ratio from ghg scenario file')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_co2   = ', irad_co2, ' in radiation_nml namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_co2')
      END SELECT
      !
      ! --- Check CH4
      ! 
      SELECT CASE (irad_ch4)
      CASE(0)
        CALL message('','irad_ch4   = 0 --> no CH4 in radiation')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_ch4   = 2 --> CH4   volume mixing ratio from radiation_nml namelist =', vmr_ch4
        CALL message('',message_text)
      CASE(3)
        WRITE (message_text, '(a,e16.8)') &
             'irad_ch4   = 3 --> CH4   tanh-profile with surface volume mixing ratio from radiation_nml namelist =', vmr_ch4
        CALL message('',message_text)
      CASE(4)
        CALL message('','irad_ch4   = 4 --> CH4   tanh-profile with surface volume mixing ratio from ghg scenario file')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_ch4   =', irad_ch4, ' in radiation_nml namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_ch4')
      END SELECT
      !
      ! --- Check N2O
      ! 
      SELECT CASE (irad_n2o)
      CASE(0)
        CALL message('','irad_n2o   = 0 --> no N2O in radiation')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_n2o   = 2 --> N2O   volume mixing ratio from radiation_nml namelist =', vmr_n2o
        CALL message('',message_text)
      CASE(3)
        WRITE (message_text, '(a,e16.8)') &
             'irad_n2o   = 3 --> N2O   tanh-profile with surface volume mixing ratio from radiation_nml namelist =', vmr_n2o
        CALL message('',message_text)
      CASE(4)
        CALL message('','irad_n2o   = 4 --> N2O   tanh-profile with surface volume mixing ratio from ghg scenario file')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_n2o   =',irad_n2o,' in radiation_nml namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_n2o')
      END SELECT
      !
      ! --- Check CFCs
      ! 
      SELECT CASE (irad_cfc11)
      CASE(0)
        CALL message('','irad_cfc11 = 0 --> no CFC11 in radiation')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_cfc11 = 2 --> CFC11 volume mixing ratio from radiation_nml namelist =', vmr_cfc11
        CALL message('',message_text)
      CASE(4)
        CALL message('','irad_cfc11 = 4 --> CFC11 volume mixing ratio from ghg scenario file')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_cfc11 =', irad_cfc11, ' in radiation_nml namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_cfc11')
      END SELECT

      SELECT CASE (irad_cfc12)
      CASE(0)
        CALL message('','irad_cfc12 = 0 --> no CFC12 in radiation')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_cfc12 = 2 --> CFC12 volume mixing ratio from radiation_nml namelist =', vmr_cfc12
        CALL message('',message_text)
      CASE(4)
        CALL message('','irad_cfc12 = 4 --> CFC12 volume mixing ratio from ghg scenario file')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_cfc12 =', irad_cfc12, ' in radiation_nml namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_cfc12')
      END SELECT
      !
      ! --- Check O3
      ! 
      SELECT CASE (irad_o3)
      CASE(0)
        CALL message('','irad_o3    = 0 --> no O3 in radiation')
!!$      CASE(2)
!!$        CALL message('','irad_o3    = 2 --> O3    periodic-in-time 3-dim. volume mixing ratio from file')
      CASE(4)
        CALL message('','irad_o3    = 4 --> O3    constant-in-time 3-dim. volume mixing ratio from file')
      CASE(8)
        CALL message('','irad_o3    = 8 --> O3    transient 3-dim. volume mixing ratio from file')
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_o3    =', irad_o3, ' in radiation_nml namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_o3')
      END SELECT
      !
      ! --- Check O2
      ! 
      SELECT CASE (irad_o2)
      CASE(0)
        CALL message('','irad_o2    = 0 --> no O2  in radiation')
      CASE(2)
        WRITE (message_text, '(a,e16.8)') &
             'irad_o2    = 2 --> O2    volume mixing ratio from radiation_nml namelist =', vmr_o2
        CALL message('',message_text)
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'irad_o2    =', irad_o2, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_psrad_radiation','Run terminated irad_o2')
      END SELECT
!!$      !
!!$      ! --- Check aerosol
!!$      ! 
!!$      SELECT CASE (irad_aero)
!!$      CASE(0)
!!$        CALL message('','irad_aero= 0 --> no aerosol in radiation')
!!$      CASE(1)
!!$        CALL message('','irad_aero= 1 --> prognostic aerosol (sub model)')
!!$      CASE(3)
!!$        CALL message('','irad_aero= 3 --> Kinne climatology')
!!$        CALL su_aero_kinne(nb_sw)
!!$      CASE(5)
!!$        CALL message('','irad_aero= 5 --> Kinne climatology + Stenchikov volcanic aerosol')
!!$        CALL su_aero_kinne(nb_sw)
!!$        CALL su_aero_volc(nb_sw)
!!$      CASE(6)
!!$        CALL message('','irad_aero= 6 --> Kinne climatology + Stenchikov volcanic aerosols + HAM volcanic aerosol')
!!$        CALL su_aero_kinne(nb_sw)
!!$        CALL su_aero_volc(nb_sw)
!!$        CALL su_aero_prop_ham
!!$        CALL read_aero_volc_tables
!!$      CASE(7)
!!$        CALL message('','irad_aero= 7 --> Kinne climatology + Crowley volcanic aerosol')
!!$        CALL su_aero_kinne(nb_sw)
!!$        CALL su_aero_prop_crow
!!$        CALL read_aero_volc_tables
!!$      CASE default
!!$        WRITE (message_text, '(a,i2,a)') &
!!$             'irad_aero=', irad_aero, ' in radctl namelist is not supported'
!!$        CALL message('',message_text)
!!$        CALL finish('setup_psrad_radiation','Run terminated irad_aero')
!!$      END SELECT
!
      !
      ! --- Check scaling factors
      !
      CALL message('','')
      CALL message('','Multiplication factors applied in radiation to vol./mass mixing ratio sources')
      CALL message('','-----------------------------------------------------------------------------')
      CALL print_value('H2O(gas,liq,ice): fh2o =',fh2o)
      CALL print_value('CO2             : fco2 =',fco2)
      CALL print_value('CH4             : fch4 =',fch4)
      CALL print_value('N2O             : fn2o =',fn2o)
      CALL print_value('O3              : fo3  =',fo3 )
      CALL print_value('O2              : fo2  =',fo2 )
      CALL print_value('CFC11 and CFC12 : fcfc =',fcfc)
      CALL message('','')
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
      CASE (4)
        CALL message('','isolrad = 4 --> solar constant for rad.-convective eq. runs with diurnal cycle ON')
      CASE (5)
        CALL message('','isolrad = 5 --> solar constant for rad.-convective eq. runs with diurnal cycle OFF')
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
      !
      ! --- Check for diagnosis of instantaneous aerosol radiative forcing
      ! 
      CALL message('','instantaneous forcing diagnostic:')
      WRITE (message_text,'(a18,L3,a20,L3,a39,I3)')       &
           ' solar radiation: ',   lradforcing(1), &
           ' thermal radiation: ', lradforcing(2), &
           ' irad_aero_forcing for reference aerosols: ', irad_aero_forcing
      CALL message('',message_text)
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
!!$      IF (irad_o3==3) CALL read_o3clim_3
!!$      !
    ENDIF
  END SUBROUTINE setup_psrad_radiation

  SUBROUTINE psrad_radiation ( &
    & current_date, &!< in current date
    & jg         ,&!< in  domain index
    & jb         ,&!< in  block index
    & kproma     ,&!< in  end index for loop over block
    & kbdim      ,&!< in  dimension of block over cells
    & klev       ,&!< in  number of full levels = number of layers
    & klevp1     ,&!< in  number of half levels = number of layer interfaces
    & ktype      ,&!< in  type of convection
    & loland     ,&!< in  land-sea mask. (1. = land, 0. = sea/lakes)
    & loglac     ,&!< in  fraction of land covered by glaciers
    & this_datetime,&!< in  actual time step
    & pcos_mu0   ,&!< in  cosine of solar zenith angle
    & geoi       ,&!< in  geopotential wrt surface at layer interfaces
    & geom       ,&!< in  geopotential wrt surface at layer centres
    & oromea     ,&!< in  orography in m
    & alb_vis_dir,&!< in  surface albedo for visible range, direct
    & alb_nir_dir,&!< in  surface albedo for near IR range, direct
    & alb_vis_dif,&!< in  surface albedo for visible range, diffuse
    & alb_nir_dif,&!< in  surface albedo for near IR range, diffuse
    & tk_sfc     ,&!< in  grid box mean surface temperature
    & pp_hl      ,&!< in  pressure at half levels at t-dt [Pa]
    & pp_fl      ,&!< in  pressure at full levels at t-dt [Pa]
    & tk_fl      ,&!< in  tk_fl  = temperature at full level at t-dt
    & xm_trc     ,&!< in  tracer mass mixing ratio
    & xm_ozn     ,&!< inout  ozone mixing ratio
    & cdnc       ,&!< in  cloud droplet number concentration
    & cld_frc    ,&!< in  cloud fraction
    & cld_cvr    ,&!< cloud cover in a column
    & vis_frc_sfc,&!< visible (250-680) fraction of net surface radiation
    & par_dn_sfc ,&!< downward Photosynth. Active Radiation (PAR) at surface
    & nir_dff_frc,&!< diffuse fraction of downw. surf. near-infrared radiation
    & vis_dff_frc,&!< diffuse fraction of downward surface visible radiation
    & par_dff_frc,&!< iffuse fraction of downward surface PAR
    & lw_flx_up_sfc,&!< longwave upward surface radiation
    & lw_net_clr_bnd,&!< clear-sky net longwave  at TOA (:,1) and surface (:,2)
    & sw_net_clr_bnd,&!< clear-sky net shortwave at TOA (:,1) and surface (:,2)
    & lw_net_clr ,&!< clear-sky net longwave  at all levels
    & sw_net_clr ,&!< clear-sky net shortwave at all levels
    & lw_net     ,&!< all-sky net longwave  at all levels
    & sw_net      &!< all-sky net shortwave at all levels
    &              )

    TYPE(datetime), POINTER, INTENT(in) :: current_date
    
    INTEGER, INTENT(in)  :: &
    & jg,             & !< domain index
    & jb,             & !< block index
    & kproma,         & !< end   index for loop over block
    & kbdim,          & !< dimension of block over cells
    & klev,           & !< number of full levels = number of layers
    & klevp1,         & !< number of half levels = number of layer interfaces
    & ktype(kbdim)      !< convection type

    LOGICAL, INTENT(IN)  :: &
    & loland(kbdim),     & !< land mask
    & loglac(kbdim)        !< glacier mask

    TYPE(datetime), POINTER :: this_datetime !< actual time step

    REAL(wp), INTENT(IN) :: &
    & pcos_mu0(kbdim),   & !< cosine of solar zenith angle
    & geoi(kbdim,klevp1),& !< geopotential wrt surface at layer interfaces
    & geom(kbdim,klev),  & !< geopotential wrt surface at layer centres
    & oromea(kbdim),     & !< orography in m
    & alb_vis_dir(kbdim),& !< surface albedo for visible range and direct light
    & alb_nir_dir(kbdim),& !< surface albedo for NIR range and direct light
    & alb_vis_dif(kbdim),& !< surface albedo for visible range and diffuse light
    & alb_nir_dif(kbdim),& !< surface albedo for NIR range and diffuse light
    & tk_sfc(kbdim),     & !< Surface temperature
    & pp_hl(kbdim,klevp1),& !< pressure at half levels [Pa]
    & pp_fl(kbdim,klev),  & !< Pressure at full levels [Pa]
    & tk_fl(kbdim,klev),  & !< Temperature on full levels [K]
    & xm_trc(kbdim,klev,ntracer), & !< tracer mixing ratio
    & cdnc(kbdim,klev),   & !< Cloud drop number concentration
    & cld_frc(kbdim,klev)   !< Cloud fraction
    REAL(wp), INTENT(INOUT) :: &
    & xm_ozn(kbdim,klev)    !< ozone
    REAL(wp), INTENT(OUT) ::      &
    & cld_cvr(:),              & !< Cloud cover in a column
    & vis_frc_sfc(kbdim),      & !< Visible (250-680) fraction of net surface radiation
    & par_dn_sfc(kbdim),       & !< Downward Photosynthetically Active Radiation (PAR) at surface
    & nir_dff_frc(kbdim),      & !< Diffuse fraction of downward surface near-infrared radiation
    & vis_dff_frc(kbdim),      & !< Diffuse fraction of downward surface visible radiation
    & par_dff_frc(kbdim),      & !< Diffuse fraction of downward surface PAR
    & lw_flx_up_sfc(kbdim),    & !< longwave upward surface radiation
    & lw_net_clr_bnd(kbdim,2), & !< Clear-sky net longwave  at TOA (:,1) and surface (:,2) 
    & sw_net_clr_bnd(kbdim,2), & !< Clear-sky net shortwave at TOA (:,1) and surface (:,2) 
    & lw_net_clr(kbdim,klevp1),& !< Clear-sky net longwave  at all levels
    & sw_net_clr(kbdim,klevp1),& !< Clear-sky net shortwave at all levels
    & lw_net(kbdim,klevp1),    & !< All-sky net longwave  at all levels
    & sw_net(kbdim,klevp1)       !< All-sky net shortwave at all levels

    INTEGER              :: jk, jl, idx(kbdim), iaero_call, number_rad_call, i_rad_call
    INTEGER              :: knwtrc  !< number of non-water tracers
    INTEGER              :: selmon  !< index to select a calendar month

    REAL(wp)             ::         &
    & cos_mu0(kbdim),               &
    & pp_sfc(kbdim),                &
!!$    & ppd_hl(kbdim,klev),           &
    & tk_hl(kbdim,klevp1),          &
    & xm_vap(kbdim,klev),           &
    & za(kbdim),                    & !< Spline interpolation arrays for qsat
    & ua(kbdim),                    &
    & xm_liq(kbdim,klev),           & !< cloud water
    & xm_ice(kbdim,klev),           & !< cloud ice
    & xc_frc(kbdim,klev),           & !< cloud fraction
    & xm_co2(kbdim,klev),           & !< CO2 mixing ratio
    & zo3_timint(kbdim,nplev_o3),   & !< intermediate value of ozon
    & xm_o3(kbdim,klev),            & !< O3 mixing ratio
    & xm_o2(kbdim,klev),            & !< O2 mixing ratio
    & xm_ch4(kbdim,klev),           & !< Methane mixing ratio
    & xm_n2o(kbdim,klev),           & !< Nitrous Oxide mixing ratio
    & xm_cfc(kbdim,klev,2),         & !< CFC mixing ratio
    & flux_factor(kbdim),           & !< 1D Scratch Array for diagnostics
    & flx_uplw    (kbdim,klevp1),   & !<   All-sky   upward longwave  flux [Wm2]
    & flx_uplw_clr(kbdim,klevp1),   & !< Clear-sky   upward longwave  flux [Wm2]
    & flx_dnlw    (kbdim,klevp1),   & !<   All-sky downward longwave  flux [Wm2]
    & flx_dnlw_clr(kbdim,klevp1),   & !< Clear-sky downward longwave  flux [Wm2]
    & flx_upsw    (kbdim,klevp1),   & !<   All-sky   upward shortwave flux [Wm2]
    & flx_upsw_clr(kbdim,klevp1),   & !< Clear-sky   upward shortwave flux [Wm2]
    & flx_dnsw    (kbdim,klevp1),   & !<   All-sky downward shortwave flux [Wm2]
    & flx_dnsw_clr(kbdim,klevp1)      !< Clear-sky downward shortwave flux [Wm2]

    TYPE(t_external_atmos_td) ,POINTER :: atm_td

    knwtrc = ntracer-iqt+1 ! tracers iqt:ntracer are non-water tracers
    
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
!!$    ppd_hl(1:kproma,:) = pp_hl(1:kproma,2:klev+1)-pp_hl(1:kproma,1:klev)
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
    !     vapor
    xm_vap(1:kproma,:) = gas_profile(kproma, klev, irad_h2o,                 &
         &                           gas_val      = xm_trc(1:kproma,:,iqv),  &
         &                           gas_factor   = fh2o)
    !     cloud water
    xm_liq(1:kproma,:) = gas_profile(kproma, klev, irad_h2o,                 &
         &                           gas_val      = xm_trc(1:kproma,:,iqc),  &
         &                           gas_epsilon  = 0.0_wp,                  &
         &                           gas_factor   = fh2o)
    !     cloud ice
    xm_ice(1:kproma,:) = gas_profile(kproma, klev, irad_h2o,                 &
         &                           gas_val      = xm_trc(1:kproma,:,iqi),  &
         &                           gas_epsilon  = 0.0_wp,                  &
         &                           gas_factor   = fh2o)
    !
    ! --- cloud cover
    ! 
    xc_frc(1:kproma,1:klev) = MERGE(cld_frc(1:kproma,1:klev), 0._wp, &
         xm_liq(1:kproma,1:klev) > 0.0_wp .OR. xm_ice(1:kproma,1:klev) > 0.0_wp)
    !
    cld_cvr(1:kproma) = 1.0_wp - xc_frc(1:kproma,1)
    DO jk = 2, klev
      cld_cvr(1:kproma) = cld_cvr(1:kproma)                                    &
           &        *(1.0_wp-MAX(xc_frc(1:kproma,jk),xc_frc(1:kproma,jk-1))) &
           &        /(1.0_wp-MIN(xc_frc(1:kproma,jk-1),1.0_wp-EPSILON(1.0_wp)))
    END DO
    cld_cvr(1:kproma) = 1.0_wp-cld_cvr(1:kproma)   
    !
    ! --- gases
    !
    ! CO2: use CO2 tracer only if the CO2 index is in the correct range
    IF ( iqt <= ico2 .AND. ico2 <= ntracer ) THEN
      xm_co2(1:kproma,:) = gas_profile(kproma, klev, irad_co2,                 &
           &                           gas_mmr      = mmr_co2,                 &
           &                           gas_scenario = ghg_co2mmr,              &
           &                           gas_val      = xm_trc(1:kproma,:,ico2), &
           &                           gas_factor   = fco2)
    ELSE
      xm_co2(1:kproma,:) = gas_profile(kproma, klev, irad_co2,       &
           &                           gas_mmr      = mmr_co2,       &
           &                           gas_scenario = ghg_co2mmr,    &
           &                           gas_factor   = fco2)
    END IF

    xm_ch4(1:kproma,:)   = gas_profile(kproma, klev, irad_ch4,       &
         &                             gas_mmr      = mmr_ch4,       &
         &                             gas_scenario = ghg_ch4mmr,    &
         &                             pressure = pp_fl, xp = ch4_v, &
         &                             gas_factor   = fch4)

    xm_n2o(1:kproma,:)   = gas_profile(kproma, klev, irad_n2o,       &
         &                             gas_mmr      = mmr_n2o,       &
         &                             gas_scenario = ghg_n2ommr,    &
         &                             pressure = pp_fl, xp = n2o_v, &
         &                             gas_factor   = fn2o)

    xm_cfc(1:kproma,:,1) = gas_profile(kproma, klev, irad_cfc11,     &
         &                             gas_mmr      = vmr_cfc11,     &
         &                             gas_scenario = ghg_cfcvmr(1), &
         &                             gas_factor   = fcfc)

    xm_cfc(1:kproma,:,2) = gas_profile(kproma, klev, irad_cfc12,     &
         &                             gas_mmr      = vmr_cfc12,     &
         &                             gas_scenario = ghg_cfcvmr(2), &
         &                             gas_factor   = fcfc)

    ! O3: provisionally construct here the ozone profiles
    atm_td => ext_data(jg)%atm_td
    SELECT CASE(irad_o3)
    CASE default
      CALL finish('radiation','o3: this "irad_o3" is not supported')
    CASE(0)
      xm_ozn(:,:) = 0.0_wp
    CASE(io3_clim, io3_ape)

      IF(irad_o3 == io3_ape) THEN
        selmon=1 ! select 1st month of file
      ELSE
        selmon=9 ! select 9th month of file
      ENDIF

      CALL o3_pl2ml ( kproma = kproma, kbdim = kbdim,        &
           &          nlev_pres = nlev_o3, klev = klev,      &
           &          pfoz = atm_td%pfoz(:),                 &
           &          phoz = atm_td%phoz(:),                 &! in o3-levs
           &          ppf  = pp_fl(:,:),                     &! in  app1
           &          pph  = pp_hl(:,:),                     &! in  aphp1
           &          o3_time_int = atm_td%o3(:,:,jb,selmon),&! in
           &          o3_clim     = xm_ozn(:,:)              )! OUT

    CASE(io3_amip)
      CALL o3_timeint(kproma = kproma, kbdim = kbdim,        &
           &          nlev_pres=nplev_o3,                    &
           &          ext_o3=o3_plev(:,:,jb,:),              &
           &          current_date=current_date,             &
           &          o3_time_int=zo3_timint                 )
      CALL o3_pl2ml ( kproma = kproma, kbdim = kbdim,         &
           &          nlev_pres = nplev_o3, klev = klev,      &
           &          pfoz = plev_full_o3,                   &
           &          phoz = plev_half_o3,                   &
           &          ppf  = pp_fl(:,:),                     &
           &          pph  = pp_hl(:,:),                     &
           &          o3_time_int = zo3_timint,              &
           &          o3_clim     = xm_ozn(:,:)              )
    END SELECT

    xm_o3(1:kproma,:)    = gas_profile(kproma, klev, irad_o3,             &
         &                             gas_scenario_v = xm_ozn(1:kproma,:), &
         &                             gas_factor     = fo3)

    xm_o2(1:kproma,:)    = gas_profile(kproma, klev, irad_o2,        &
         &                             gas_mmr      = mmr_o2,        &
         &                             gas_factor   = fo2)

    ! 2.0 Radiation used to advance model, provide standard diagnostics, and radiative forcing if desired
    !
    ! --------------------------------
    ! 2.1 Radiation call (number of calls depends on whether forcing is desired)
    ! --------------------------------
    number_rad_call = 1
    IF (lradforcing(1).OR.lradforcing(2)) number_rad_call = 2 

    DO i_rad_call = 1,number_rad_call
      iaero_call = irad_aero
      IF (i_rad_call < number_rad_call) iaero_call = irad_aero_forcing

      CALL psrad_interface(   current_date    ,jg                               ,&
           & iaero_call      ,kproma          ,kbdim           ,klev            ,& 
!!$           & jb              ,knwtrc          ,ktype           ,nb_sw           ,&
           & jb                               ,ktype           ,nb_sw           ,&
           & loland          ,loglac          ,cemiss          ,this_datetime   ,&
           & cos_mu0         ,geoi            ,geom            ,oromea          ,&
           & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
           & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
           & tk_hl           ,tk_sfc          ,xm_vap          ,xm_liq          ,&
           & xm_ice          ,cdnc            ,xc_frc          ,xm_o3           ,&
           & xm_co2          ,xm_ch4          ,xm_n2o          ,xm_cfc          ,&
!!$           & xm_o2           ,xm_trc(:,:,iqt:ntracer)                           ,&
           & xm_o2                                                              ,&
           & flx_uplw        ,flx_uplw_clr    ,flx_dnlw        ,flx_dnlw_clr    ,&
           & flx_upsw        ,flx_upsw_clr    ,flx_dnsw        ,flx_dnsw_clr    ,&
           & vis_frc_sfc     ,par_dn_sfc      ,nir_dff_frc     ,vis_dff_frc     ,&
           & par_dff_frc                                                         )
      !
      ! Compute net fluxes from up/down fluxes, and normalize solar fluxes by current value of solar constant
      ! and zenith angle as they are renormalized when heating rates are calculated.
      !
      flux_factor(1:kproma) = 1._wp / (psctm*cos_mu0(1:kproma))
      lw_net    (1:kproma,1:klevp1) =  flx_dnlw    (1:kproma, 1:klevp1) - flx_uplw    (1:kproma, 1:klevp1)       
      lw_flx_up_sfc (1:kproma)      =  flx_uplw(1:kproma,klevp1)
      sw_net    (1:kproma,1:klevp1) = (flx_dnsw    (1:kproma, 1:klevp1) - flx_upsw    (1:kproma, 1:klevp1)) * &
           &                      SPREAD(flux_factor(1:kproma),2,klevp1)       
      lw_net_clr(1:kproma,1:klevp1) =  flx_dnlw_clr(1:kproma, 1:klevp1) - flx_uplw_clr(1:kproma, 1:klevp1)       
      sw_net_clr(1:kproma,1:klevp1) = (flx_dnsw_clr(1:kproma, 1:klevp1) - flx_upsw_clr(1:kproma, 1:klevp1)) * &
           &                      SPREAD(flux_factor(1:kproma),2,klevp1)
      par_dn_sfc(1:kproma) = par_dn_sfc(1:kproma) * flux_factor(1:kproma)

      IF (i_rad_call < number_rad_call) CALL prepare_psrad_radiation_forcing(jg, & 
           & kproma          ,kbdim           ,klevp1          ,jb          ,&
           & lw_net          ,sw_net          ,lw_net_clr      ,sw_net_clr   )
    END DO
    !
    ! 2.1 Fluxes to advance to the model, compute cloud radiative effect 
    ! --------------------------------
    !
    ! --- Total (net) fluxes, used to advance the model 
    !
    !
    ! --- Clear sky fluxes
    lw_net_clr_bnd(1:kproma,1)    = lw_net_clr(1:kproma,1)
    lw_net_clr_bnd(1:kproma,2)    = lw_net_clr(1:kproma,klevp1)
    sw_net_clr_bnd(1:kproma,1)    = sw_net_clr(1:kproma,1)        
    sw_net_clr_bnd(1:kproma,2)    = sw_net_clr(1:kproma,klevp1) 

  END SUBROUTINE psrad_radiation

  !---------------------------------------------------------------------------
  !>
  !! GAS_PROFILE:  Determines Gas distributions based on case specification
  !! 
  !! @par Revsision History 
  !! B. Stevens (2009-08). 
  !! H. Schmidt (2010-08): profile calculation added for scenario case.
  !!
  !! Description: This routine calculates the gas distributions for one of
  !! five cases:  (0) no gas present; (1) prognostic gas; (2) specified 
  !! mixing ratio; (3) mixing ratio decaying with height given profile;
  !! (4) scenario run with different mixing ratio, if profile parameters are
  !! given a vertical profile is calculated as in (3).
  !
  FUNCTION gas_profile (kproma, klev, igas, gas_mmr, gas_scenario, gas_mmr_v, &
       &                gas_scenario_v, gas_val, xp, pressure,                &
       &                gas_epsilon, gas_factor)

    INTEGER, INTENT (IN) :: kproma, klev, igas
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_mmr              ! for igas = 2 and 3
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_scenario         ! for igas = 4
    REAL (wp), OPTIONAL, INTENT (IN) :: pressure(:,:), xp(3) ! for igas = 3 and 4
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_mmr_v(:,:)       ! for igas = 2
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_scenario_v(:,:)  ! for igas = 4
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_val(:,:)         ! for igas = 1
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_epsilon
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_factor

    REAL (wp) :: gas_profile(kproma,klev), zx_d, zx_m, eps, fgas
    LOGICAL :: gas_initialized

    gas_initialized = .FALSE.

    IF (PRESENT(gas_epsilon)) THEN
       eps = gas_epsilon
    ELSE
       eps = EPSILON(1.0_wp)
    END IF

    IF (PRESENT(gas_factor)) THEN
       fgas = gas_factor
    ELSE
       fgas = 1.0_wp
    END IF

    SELECT CASE (igas)

    CASE (0)                             ! 0: set concentration to zero
      gas_profile(1:kproma,:) = eps
      gas_initialized = .TRUE.

    CASE (1)                             ! 1: horizontally and vertically variable
      IF (PRESENT(gas_val)) THEN
        gas_profile(1:kproma,:) = MAX(gas_val(1:kproma,:)*fgas, eps)
        gas_initialized = .TRUE.
      END IF

    CASE (2)
      IF (PRESENT(gas_mmr)) THEN         ! 2a: horizontally and vertically constant
        gas_profile(1:kproma,:) = MAX(gas_mmr*fgas, eps)
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_mmr_v)) THEN  ! 2b: = (1)
        gas_profile(1:kproma,:) = MAX(gas_mmr_v(1:kproma,:)*fgas, eps)
        gas_initialized = .TRUE.
      END IF

    CASE (3)                             ! 3: horizontally constant and tanh-profile in the vertical
      IF (PRESENT(gas_mmr) .AND. PRESENT(xp) .AND. PRESENT(pressure)) THEN
        zx_m = (gas_mmr+xp(1)*gas_mmr)*0.5_wp
        zx_d = (gas_mmr-xp(1)*gas_mmr)*0.5_wp
        gas_profile(1:kproma,:)=MAX((1-(zx_d/zx_m)*TANH(LOG(pressure(1:kproma,:)   &
             &                      /xp(2)) /xp(3))) * zx_m * fgas, eps)
        gas_initialized = .TRUE.
      END IF

    CASE (4,8)
      IF (PRESENT(gas_scenario)) THEN
        IF (PRESENT(xp) .AND. PRESENT(pressure)) THEN ! 4a: = (3)
          ! comment H. Schmidt: If the respective parameters are present, a vertical 
          ! profile is calculated as in option (3). This allows a seamless
          ! continuation of preindustrial control with scenarios. The treatment here is
          ! inconsistent with having two different options for the constant
          ! concentration cases (2 without and 3 with profile). However, instead
          ! of adding a fifth option, it seems more advisable to clean up the 
          ! complete handling of radiation switches (including ighg), later.
          zx_m = (gas_scenario+xp(1)*gas_scenario)*0.5_wp
          zx_d = (gas_scenario-xp(1)*gas_scenario)*0.5_wp
          gas_profile(1:kproma,:)=MAX((1-(zx_d/zx_m)*TANH(LOG(pressure(1:kproma,:)   &
             &                        /xp(2)) /xp(3))) * zx_m * fgas, eps)
        ELSE                                          ! 4b: = (2a)
          gas_profile(1:kproma,:)=MAX(gas_scenario*fgas, eps)
        ENDIF
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_scenario_v)) THEN          ! 4c: = (1)
        gas_profile(1:kproma,:) = MAX(gas_scenario_v(1:kproma,:)*fgas, eps)
        gas_initialized = .TRUE.
      END IF

    END SELECT

    IF (.NOT. gas_initialized) &
         CALL finish('radiation','gas_profile options not supported')

  END FUNCTION gas_profile
  !---------------------------------------------------------------------------


  END MODULE mo_psrad_radiation
