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
  USE mo_physical_constants,       ONLY: vtmpc1, rae,           &
       &                        amco2, amch4, amn2o, amo2, amd
!  USE mo_control,         ONLY: lcouple, lmidatm
!  USE mo_time_base,       ONLY: get_calendar_type, JULIAN
  USE mo_exception,       ONLY: finish, message, message_text
!  USE mo_mpi,             ONLY: p_parallel, p_parallel_io, p_bcast, p_io
!  USE mo_namelist,        ONLY: open_nml, position_nml, POSITIONED
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

! the present mo_bc_solar_irradiance is "old" and not the one used in echam-6.3.
!  USE mo_solar_irradiance,ONLY: get_solar_irradiance, set_solar_irradiance, &
!                                get_solar_irradiance_m, set_solar_irradiance_m
! cloud optics does not exist in icon
!  USE mo_cloud_optics,    ONLY: setup_cloud_optics  
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

! must be introduced to icon
!  USE mo_rrtm_params,   ONLY : nbndsw
! new to icon
!  USE mo_srtm_setup,    ONLY : ssi_default, ssi_preind, ssi_amip,           &
!                             & ssi_RCEdiurnOn, ssi_RCEdiurnOff
!  USE mo_psrad_interface,ONLY : setup_psrad, psrad_interface, &
!                                lw_strat, sw_strat
!  USE mo_spec_sampling, ONLY : spec_sampling_strategy, &
!                             & set_spec_sampling_lw, set_spec_sampling_sw, get_num_gpoints

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: pre_psrad_radiation, setup_psrad_radiation, psrad_radiation

  CONTAINS

  SUBROUTINE pre_psrad_radiation(datetime_radiation,ltrig_rad)
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !

    TYPE(t_datetime), INTENT(IN)     :: datetime_radiation !< date and time of radiative transfer calculation
    LOGICAL         , INTENT(IN)     :: ltrig_rad !< .true. if radiative transfer calculation has to be done at current time step

    LOGICAL  :: l_rad_call, l_write_solar
    INTEGER  :: icurrentyear, icurrentmonth, iprevmonth, i
    REAL(wp) :: rasc_sun, decl_sun, dist_sun, time_of_day, zrae
    REAL(wp) :: solcm, orbit_date
    LOGICAL  :: l_orbvsop87
    !
    ! 1.0 Compute orbital parameters for current time step
    ! --------------------------------
    l_rad_call = .FALSE.
    CALL get_orbit_times(datetime_radiation, time_of_day, &
         &               orbit_date)

    IF (l_orbvsop87) THEN 
      CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
    ELSE
      CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
    END IF
!!$    decl_sun_cur = decl_sun       ! save for aerosol and chemistry submodels
!!$    CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
!!$         &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
!!$         &                flx_ratio_cur, amu0_x, rdayl_x)
!!$
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
    IF (phy_config%lrad .AND. ltrig_rad) THEN
      l_rad_call = .TRUE.
      CALL get_orbit_times(datetime_radiation, time_of_day , &
           &               orbit_date)

      IF ( l_orbvsop87 ) THEN 
        CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
      ELSE
        CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
      END IF
!!$      CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
!!$           &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
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

    END IF ! lrad .AND. l_trigrad

  END SUBROUTINE pre_psrad_radiation

  SUBROUTINE setup_psrad_radiation
  END SUBROUTINE setup_psrad_radiation

  SUBROUTINE psrad_radiation
  END SUBROUTINE psrad_radiation

  END MODULE mo_psrad_radiation
