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

  USE mo_kind,                ONLY: wp, i8
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants      ,ONLY: min_rlcell_int, grf_bdywidth_c

  USE mo_physical_constants,  ONLY: rae, amd, amco2, amch4, amn2o, amo2, amc11, amc12
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_run_config,          ONLY: iqv, iqc, iqi, ico2, io3, &
                                    ntracer, lart
  USE mo_echam_phy_config,    ONLY: echam_phy_tc
  USE mo_echam_rad_config,    ONLY: echam_rad_config
  USE mo_bc_greenhouse_gases, ONLY: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr
  USE mo_bc_ozone,            ONLY: ext_ozone
  USE mo_o3_util,             ONLY: o3_pl2ml, o3_timeint

  USE mo_psrad_orbit,         ONLY: orbit_kepler, orbit_vsop87, get_orbit_times
  USE mo_psrad_solar_data,    ONLY : ssi_default, ssi_amip,                    &
                                     ssi_cmip5_picontrol, ssi_cmip6_picontrol, &
                                     ssi_RCEdiurnOn, ssi_RCEdiurnOff,          &
                                     ssi_radt, tsi_radt
  USE mo_psrad_solar_parameters, ONLY:                         &
                                     psctm,                    &
                                     ssi_factor,               &
                                     solar_parameters

  USE mo_psrad_general,       ONLY: nbndsw

  USE mo_atmo_psrad_interface,ONLY : atmo_psrad_interface

  USE mtime, ONLY: datetime, getTotalSecondsTimeDelta !, datetimeToString

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: pre_psrad_radiation, psrad_radiation

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
    LOGICAL,                 INTENT(in) :: ltrig_rad !< .true. if SW radiative transfer calculation has to be done at current time step
    REAL(wp),                INTENT(out) :: amu0_x(:,:), rdayl_x(:,:), &
         &                                  amu0m_x(:,:), rdaylm_x(:,:)

    LOGICAL  :: l_write_solar
    INTEGER(i8) :: icurrentyear
    INTEGER  :: icurrentmonth, i
    INTEGER, SAVE :: iprevmonth=-9999
    REAL(wp) :: rasc_sun, decl_sun, dist_sun
    REAL(wp) :: time_of_day, orbit_date
    REAL(wp) :: time_of_day_rt, orbit_date_rt
    REAL(wp) :: dt_ext
    REAL(wp) :: tsi

    ! Shortcuts to components of echam_rad_config
    !
    LOGICAL , POINTER :: l_orbvsop87, ldiur, l_sph_symm_irr, lyr_perp
    INTEGER , POINTER :: isolrad, icosmu0, yr_perp
    REAL(wp), POINTER :: fsolrad, cecc, cobld, clonp
    INTEGER  :: jg
!     CHARACTER(LEN=32)               :: datestring, datestring2
!     CHARACTER(LEN=*), PARAMETER     :: method_name="pre_psrad_radiation"

    !
    jg = p_patch%id
    isolrad        => echam_rad_config(jg)% isolrad
    fsolrad        => echam_rad_config(jg)% fsolrad
    l_orbvsop87    => echam_rad_config(jg)% l_orbvsop87
    cecc           => echam_rad_config(jg)% cecc
    cobld          => echam_rad_config(jg)% cobld
    clonp          => echam_rad_config(jg)% clonp
    lyr_perp       => echam_rad_config(jg)% lyr_perp
    yr_perp        => echam_rad_config(jg)% yr_perp
    ldiur          => echam_rad_config(jg)% ldiur
    l_sph_symm_irr => echam_rad_config(jg)% l_sph_symm_irr
    icosmu0        => echam_rad_config(jg)% icosmu0

    !
    ! 1.0 Compute orbital parameters for current time step
    ! --------------------------------
    !
    ! "time_of_day"  defines the local noon longitude for the time of this time step.
    ! "orbit_date" is not used. Instead always "orbit_date_rt" is used.
    CALL get_orbit_times(current_datetime, lyr_perp, yr_perp, time_of_day, orbit_date)
    !
    ! "time_of_day_rt" defines the local noon longitude for the time of the
    ! radiative transfer calculation.
    ! "orbit_date_rt" defines the orbit position at the radiation time.
    ! This orbit position is kept constant through the radiation interval.
    CALL get_orbit_times(datetime_radiation, lyr_perp, yr_perp, time_of_day_rt, orbit_date_rt)


    ! Compute the orbital parameters of Earth for "orbit_date_rt".
    IF (l_orbvsop87) THEN 
      CALL orbit_vsop87 (             orbit_date_rt, rasc_sun, decl_sun, dist_sun)
    ELSE
      CALL orbit_kepler (cecc, cobld, clonp, orbit_date_rt, rasc_sun, decl_sun, dist_sun)
    END IF

!     CALL datetimeToString(current_datetime,   datestring)
!     CALL datetimeToString(datetime_radiation, datestring2)
!     WRITE(message_text,'(4a)') ' current datetime: ', datestring, '; radiation datetime: ', datestring2
!     CALL message (method_name, message_text)


    dt_ext = 0.0_wp
    ! Compute cos(zenith angle) "amu0_x" for the current time "time_of_day", if the
    ! SW fluxes should be adjusted to the current sun (icosmu0=1:4), or the radiation
    ! time "time_of_day_rt", if the heating shall be computed for the sun position
    ! used for the radiative transfer (icosmu0=0).
    ! In both cases use the orbit parameters "decl_sun" and "dist_sun" valid for
    ! the "orbit_date_rt".
    SELECT CASE (icosmu0)
    CASE (0)
       CALL solar_parameters(decl_sun,        time_of_day_rt,  &
            &                icosmu0,         dt_ext,          &
            &                ldiur,           l_sph_symm_irr,  &
            &                p_patch,                          &
            &                amu0_x,          rdayl_x          )
    CASE (1:4)
       CALL solar_parameters(decl_sun,        time_of_day,     &
            &                icosmu0,         dt_ext,          &
            &                ldiur,           l_sph_symm_irr,  &
            &                p_patch,                          &
            &                amu0_x,          rdayl_x          )
    CASE DEFAULT
       CALL finish('mo_psrad_radiation/pre_psrad_radiation','invalid icosmu0, must be in 0:4')
    END SELECT
    !
    !
    ! 2.0 Prepare time dependent quantities for rad (on radiation timestep)
    ! --------------------------------
    IF (ltrig_rad) THEN

      ! Compute cos(zenith angle) "amu0m_x" for the radiation time "time_of_day_rt"
      ! and the orbit parameters "decl_sun" and "dist_sun" valid for "orbit_date_rt".
      !
      ! "amu0m_x" is needed for the incoming SW flux field at the top of the atmosphere
      ! and for the optical paths.
      !
      SELECT CASE (icosmu0)
      CASE (0)
         dt_ext = 0.0_wp
      CASE (1:4)
         ! Extend the sunlit area for the radiative transfer calculations over an extended area
         ! including a rim of width dt_rad/2/86400*2pi (in radian) around the sunlit hemisphere.
         dt_ext = getTotalSecondsTimeDelta(echam_phy_tc(p_patch%id)%dt_rad,current_datetime)
      END SELECT
      !
      CALL solar_parameters(decl_sun,        time_of_day_rt,   &
           &                icosmu0,         dt_ext,           &
           &                ldiur,           l_sph_symm_irr,   &
           &                p_patch,                           &
           &                amu0m_x,         rdaylm_x          )
      !
      ! Consider curvature of the atmosphere for high zenith angles:
      ! The atmospheric path for a zenith angle mu0 through a spherical shell of
      ! thickness H and an inner radius R (and ratio rae=H/R) is shorter than
      ! a path at equal zenith angle through a plane parallel layer of thickness H.
      !
      WHERE (rdaylm_x == 1.0_wp) 
         amu0m_x(:,:)  = rae/(SQRT(amu0m_x(:,:)**2+rae*(rae+2.0_wp))-amu0m_x(:,:))
      END WHERE

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
        tsi = SUM(ssi_default)
        ssi_factor = ssi_default
      CASE (1)
        tsi=tsi_radt
        ssi_factor=ssi_radt
        continue ! solar irradiance was read in echam_phy_bcs_global
      CASE (2)
        tsi = SUM(ssi_cmip5_picontrol)
        ssi_factor = ssi_cmip5_picontrol
      CASE (3)
        tsi = SUM(ssi_amip)
        ssi_factor = ssi_amip
      CASE (4)
        tsi = SUM(ssi_RCEdiurnOn)
        ssi_factor = ssi_RCEdiurnOn
      CASE (5)
        tsi = SUM(ssi_RCEdiurnOff)
        ssi_factor = ssi_RCEdiurnOff
      CASE (6)
        tsi = SUM(ssi_cmip6_picontrol)
        ssi_factor = ssi_cmip6_picontrol
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'isolrad = ', isolrad, ' in radctl namelist is not supported'
        CALL message('pre_radiation', message_text)
      END SELECT
      psctm = tsi/dist_sun**2 * fsolrad
      ssi_factor(:) = ssi_factor(:)/tsi

      ! output of solar constant every month

      icurrentmonth = datetime_radiation%date%month
      icurrentyear = datetime_radiation%date%year
      l_write_solar = icurrentmonth/=iprevmonth
      iprevmonth = icurrentmonth
      IF (l_write_solar) THEN
        CALL message('','')
        WRITE (message_text,'(a,i0,a,i2.2,a,f6.1)') &
             'Total solar constant [W/m^2] for ',      &
             icurrentyear, '-',                        &
             icurrentmonth, ' = ', tsi
        CALL message('',message_text)
        CALL message('','')
        DO i = 1, nbndsw
          WRITE (message_text,'(a,i2,a,f7.5)') &
               '   solar constant fraction: band ', i, &
               ' = ', ssi_factor(i)
          CALL message('',message_text)
        END DO
      END IF
    END IF ! ltrig_rad

  END SUBROUTINE pre_psrad_radiation
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE psrad_radiation( &
    & patch          ,&!< in  domain index
    & klev           ,&!< in  number of full levels = number of layers
    & klevp1         ,&!< in  number of half levels = number of layer interfaces
    & ktype          ,&!< in  type of convection
    & loland         ,&!< in  land-sea mask. (1. = land, 0. = sea/lakes)
    & loglac         ,&!< in  fraction of land covered by glaciers
    & this_datetime  ,&!< in  actual time step
    & pcos_mu0       ,&!< in  cosine of solar zenith angle
    & daylght_frc    ,&!< in  daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
    & alb_vis_dir    ,&!< in  surface albedo for visible range, direct
    & alb_nir_dir    ,&!< in  surface albedo for near IR range, direct
    & alb_vis_dif    ,&!< in  surface albedo for visible range, diffuse
    & alb_nir_dif    ,&!< in  surface albedo for near IR range, diffuse
    & tk_sfc         ,&!< in  grid box mean surface temperature
    & zf             ,&!< in  geometric height at full level      [m]
    & zh             ,&!< in  geometric height at half level      [m]
    & dz             ,&!< in  geometric height thickness of layer [m]
    & pp_hl          ,&!< in  pressure at half levels at t-dt [Pa]
    & pp_fl          ,&!< in  pressure at full levels at t-dt [Pa]
    & tk_fl          ,&!< in  tk_fl  = temperature at full level at t-dt
    & xm_dry         ,&!< in  dry air mass in layer [kg/m2]
    & xm_trc         ,&!< in  tracer  mass in layer [kg/m2]
    & xm_ozn         ,&!< inout ozone mass mixing ratio [kg/kg]
    !
    & cdnc           ,&!< in  cloud droplet number concentration
    & cld_frc        ,&!< in  cloud fraction
    & cld_cvr        ,&!< out cloud cover in a column
    !
    & lw_dnw_clr     ,&!< out clear-sky downward longwave  at all levels
    & lw_upw_clr     ,&!< out clear-sky upward   longwave  at all levels
    & sw_dnw_clr     ,&!< out clear-sky downward shortwave at all levels
    & sw_upw_clr     ,&!< out clear-sky upward   shortwave at all levels
    & lw_dnw         ,&!< out all-sky   downward longwave  at all levels
    & lw_upw         ,&!< out all-sky   upward   longwave  at all levels
    & sw_dnw         ,&!< out all-sky   downward shortwave at all levels
    & sw_upw         ,&!< out all-sky   upward   shortwave at all levels
    !
    & vis_dn_dir_sfc ,&!< out all-sky downward direct visible radiation at surface
    & par_dn_dir_sfc ,&!< out all-sky downward direct PAR     radiation at surface
    & nir_dn_dir_sfc ,&!< out all-sky downward direct near-IR radiation at surface
    & vis_dn_dff_sfc ,&!< out all-sky downward diffuse visible radiation at surface
    & par_dn_dff_sfc ,&!< out all-sky downward diffuse PAR     radiation at surface
    & nir_dn_dff_sfc ,&!< out all-sky downward diffuse near-IR radiation at surface
    & vis_up_sfc     ,&!< out all-sky upward visible radiation at surface
    & par_up_sfc     ,&!< out all-sky upward PAR     radiation at surfac
    & nir_up_sfc     ) !< out all-sky upward near-IR radiation at surface

    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch

    INTEGER, INTENT(in)     :: &
    & klev,                    & !< number of full levels = number of layers
    & klevp1,                  & !< number of half levels = number of layer interfaces
    & ktype(:,:)               !< convection type

    LOGICAL, INTENT(IN)     :: &
    & loland(:,:),           & !< land mask
    & loglac(:,:)              !< glacier mask

    TYPE(datetime), POINTER :: this_datetime !< actual time step

    REAL(wp), INTENT(IN)    :: &
    & pcos_mu0(:,:),         & !< cosine of solar zenith angle
    & daylght_frc(:,:),      & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
    & alb_vis_dir(:,:),      & !< surface albedo for visible range and direct light
    & alb_nir_dir(:,:),      & !< surface albedo for NIR range and direct light
    & alb_vis_dif(:,:),      & !< surface albedo for visible range and diffuse light
    & alb_nir_dif(:,:),      & !< surface albedo for NIR range and diffuse light
    & tk_sfc(:,:),           & !< Surface temperature
    & zf(:,:,:),          & !< geometric height at full level      [m]
    & zh(:,:,:),          & !< geometric height at half level      [m]
    & dz(:,:,:),          & !< geometric height thickness of layer [m]
    & xm_dry(:,:,:),      & !< dry air mass in layer [kg/m2]
    & pp_hl(:,:,:),       & !< pressure at half levels [Pa]
    & pp_fl(:,:,:),       & !< Pressure at full levels [Pa]
    & tk_fl(:,:,:),       & !< Temperature on full levels [K]
    & xm_trc(:,:,:,:),    & !< tracer mass in layer [kg/m2]
    & cdnc(:,:,:),        & !< Cloud drop number concentration
    & cld_frc(:,:,:)        !< Cloud fraction
    REAL(wp), INTENT(INOUT) :: &
    & xm_ozn(:,:,:)         !< ozone mixing ratio  [kg/kg]

    ! OUT
    REAL(wp), INTENT(INOUT)   :: &
    & cld_cvr(:,:),            & !< Cloud cover in a column
    & lw_dnw_clr(:,:,:),& !< Clear-sky downward longwave  at all levels
    & lw_upw_clr(:,:,:),& !< Clear-sky upward   longwave  at all levels
    & sw_dnw_clr(:,:,:),& !< Clear-sky downward shortwave at all levels
    & sw_upw_clr(:,:,:),& !< Clear-sky upward   shortwave at all levels
    & lw_dnw(:,:,:),    & !< All-sky   downward longwave  at all levels
    & lw_upw(:,:,:),    & !< All-sky   upward   longwave  at all levels
    & sw_dnw(:,:,:),    & !< All-sky   downward shortwave at all levels
    & sw_upw(:,:,:)       !< All-sky   upward   shortwave at all levels

    REAL (wp), INTENT (INOUT) :: &
    & vis_dn_dir_sfc(:,:)  , & !< Diffuse downward flux surface visible radiation 
    & par_dn_dir_sfc(:,:)  , & !< Diffuse downward flux surface PAR
    & nir_dn_dir_sfc(:,:)  , & !< Diffuse downward flux surface near-infrared radiation
    & vis_dn_dff_sfc(:,:)  , & !< Direct  downward flux surface visible radiation 
    & par_dn_dff_sfc(:,:)  , & !< Direct  downward flux surface PAR
    & nir_dn_dff_sfc(:,:)  , & !< Direct  downward flux surface near-infrared radiation
    & vis_up_sfc    (:,:)  , & !< Upward  flux surface visible radiation 
    & par_up_sfc    (:,:)  , & !< Upward  flux surface PAR
    & nir_up_sfc    (:,:)      !< Upward  flux surface near-infrared radiation


    REAL (wp) ::      &
    & xm_vap(nproma,klev, patch%nblks_c),           & !< water vapor mass in layer [kg/m2]
    & xm_liq(nproma,klev, patch%nblks_c),           & !< cloud water mass in layer [kg/m2]
    & xm_ice(nproma,klev, patch%nblks_c),           & !< cloud ice   mass in layer [kg/m2]
    & xm_co2(nproma,klev, patch%nblks_c),           & !< CO2 mass in layer [kg/m2]
    & xm_o3(nproma,klev, patch%nblks_c),            & !< O3  mass in layer [kg/m2]
    & xm_o2(nproma,klev, patch%nblks_c),            & !< O2  mass in layer [kg/m2]
    & xm_ch4(nproma,klev, patch%nblks_c),           & !< CH4 mass in layer [kg/m2]
    & xm_n2o(nproma,klev, patch%nblks_c),           & !< N2O mass in layer [kg/m2]
    & xm_cfc(nproma,klev,2, patch%nblks_c),         &!< CFC mass in layer [kg/m2]
    & xc_frc(nproma,klev, patch%nblks_c)              !< cloud fraction

    REAL (wp) ::      &
    & tk_hl(nproma,klevp1, patch%nblks_c)

    REAL (wp) :: &
    & pp_sfc(nproma,patch%nblks_c)

    INTEGER  :: jg             
    INTEGER  :: i_nchdom, rl_start, rl_end
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    ! Shortcuts to components of echam_rad_config
    !
    INTEGER , POINTER :: irad_aero
    !
    jg         =  patch%id
    irad_aero  => echam_rad_config(jg)% irad_aero
    
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
 
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      CALL calculate_temperature_pressure (        &
        & jcs            = jcs                    ,&!< in  start index for loop over block
        & jce            = jce                    ,&!< in  end   index for loop over block
        & kbdim          = nproma                 ,&!< in  dimension of block over cells
        & klev           = klev                   ,&!< in  number of full levels = number of layers
        & klevp1         = klevp1                 ,&!< in  number of half levels = number of layer interfaces
        ! in 
        & pp_hl          = pp_hl(:,:,jb)          ,&
        & pp_fl          = pp_fl(:,:,jb)          ,&
        & tk_fl          = tk_fl(:,:,jb)          ,&
        & tk_sfc         = tk_sfc(:,jb)           ,&
        !out      
        & pp_sfc         = pp_sfc(:,jb)           ,&
        & tk_hl          = tk_hl(:,:,jb))         

   END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
  
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       
      CALL psrad_get_gas_profiles (                &
        & jg             = jg                     ,&!< in  domain index
        & jb             = jb                     ,&!< in  block index
        & jcs            = jcs                    ,&!< in  start index for loop over block
        & jce            = jce                    ,&!< in  end   index for loop over block
        & kbdim          = nproma                 ,&!< in  dimension of block over cells
        & klev           = klev                   ,&!< in  number of full levels = number of layers
        & klevp1         = klevp1                 ,&!< in  number of half levels = number of layer interfaces
        & this_datetime  = this_datetime          ,&!< in  actual time step
        & pp_hl          = pp_hl  (:,:,jb)        ,&!< in  pressure at half levels at t-dt [Pa]
        & pp_fl          = pp_fl  (:,:,jb)        ,&!< in  pressure at full levels at t-dt [Pa]
        & xm_dry         = xm_dry (:,:,jb)        ,&!< in  dry air mass in layer [kg/m2]
        & xm_trc         = xm_trc (:,:,jb,:)      ,&!< in  tracer  mass in layer [kg/m2]
        & cld_frc        = cld_frc(:,:,jb)        ,&!< in   cloud fraction [m2/m2]
        & xm_ozn         = xm_ozn (:,:,jb)        ,&!< inout  ozone  mass mixing ratio [kg/kg]
        & xm_vap         = xm_vap (:,:,jb)        ,& !< water vapor mass in layer [kg/m2]
        & xm_liq         = xm_liq (:,:,jb)        ,& !< cloud water mass in layer [kg/m2]
        & xm_ice         = xm_ice (:,:,jb)        ,& !< cloud ice   mass in layer [kg/m2]
        & xm_co2         = xm_co2 (:,:,jb)        ,& !< CO2 mass in layer [kg/m2]
        & xm_o3          = xm_o3  (:,:,jb)        ,& !< O3  mass in layer [kg/m2]
        & xm_o2          = xm_o2  (:,:,jb)        ,& !< O2  mass in layer [kg/m2]
        & xm_ch4         = xm_ch4 (:,:,jb)        ,& !< CH4 mass in layer [kg/m2]
        & xm_n2o         = xm_n2o (:,:,jb)        ,& !< N2O mass in layer [kg/m2]
        & xm_cfc         = xm_cfc (:,:,:,jb)      ,& !< CFC mass in layer [kg/m2]
        & xc_frc         = xc_frc (:,:,jb)        ,&
        & cld_cvr        = cld_cvr(:,jb)          )

    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
!     CALL psrad_interface(                                                   &
    CALL atmo_psrad_interface(                                                   &
      & patch,                                                              &
      & irad_aero     ,klev                                                ,& 
      & ktype                                                              ,&
      & psctm, ssi_factor,                                                  &
      & loland          ,loglac          ,this_datetime                    ,&
      & pcos_mu0        ,daylght_frc                                       ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & zf              ,zh              ,dz                               ,&
      & pp_sfc          ,pp_fl                                             ,&
      & tk_sfc          ,tk_fl           ,tk_hl                            ,&
      & xm_dry          ,xm_vap          ,xm_liq          ,xm_ice          ,&
      & cdnc            ,xc_frc                                            ,&
      & xm_co2          ,xm_ch4          ,xm_n2o          ,xm_cfc          ,&
      & xm_o3           ,xm_o2                                             ,&
      & lw_upw          ,lw_upw_clr      ,lw_dnw          ,lw_dnw_clr      ,&
      & sw_upw          ,sw_upw_clr      ,sw_dnw          ,sw_dnw_clr      ,&
      & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
      & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
      & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       )     
     !-------------------------------------------------------------------

  END SUBROUTINE psrad_radiation
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE calculate_temperature_pressure (               &
      & jcs           ,&!< in  start index for loop over block
      & jce           ,&!< in  end   index for loop over block
      & kbdim         ,&!< in  dimension of block over cells
      & klev          ,&!< in  number of full levels = number of layers
      & klevp1        ,&!< in  number of half levels = number of layer interfaces
      & pp_hl         ,&! in
      & pp_fl         ,&
      & tk_fl         ,&
      & tk_sfc        ,&
      & pp_sfc        ,& ! out
      & tk_hl   )

  
    INTEGER, INTENT(in) :: &
      & jcs               ,&
      & jce               ,&
      & kbdim             ,&
      & klev              ,&
      & klevp1
    ! in 
    REAL(wp), INTENT(in) :: &
      & pp_hl(kbdim,klevp1),&
      & pp_fl(kbdim,klev)  ,&
      & tk_fl(kbdim,klev)  ,&
      & tk_sfc(kbdim)

    !out      
    REAL(wp), INTENT(inout) :: &
      & pp_sfc (kbdim)        ,&
      & tk_hl  (kbdim,klevp1)  
    
    INTEGER              :: jk, jl
   !
    ! 1.0 calculate variable input parameters (location and state variables)
    ! --------------------------------
    ! 
    ! --- Pressure (surface and distance between half levels)
    !
    pp_sfc(jcs:jce)   = pp_hl(jcs:jce,klevp1)
    !
    ! --- temperature at half levels
    !
    DO jk=2,klev
      DO jl = jcs, jce
        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
             & - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
             & - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
      END DO
    END DO
    DO jl = jcs, jce
      tk_hl(jl,klevp1) = tk_sfc(jl)
      tk_hl(jl,1)      = tk_fl(jl,1)-pp_fl(jl,1)*(tk_fl(jl,1) - tk_hl(jl,2))  &
           &             / (pp_fl(jl,1)-pp_hl(jl,2))
    END DO


  END SUBROUTINE calculate_temperature_pressure 
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE psrad_get_gas_profiles ( &
    & jg             ,&!< in  domain index
    & jb             ,&!< in  block index
    & jcs            ,&!< in  start index for loop over block
    & jce            ,&!< in  end   index for loop over block
    & kbdim          ,&!< in  dimension of block over cells
    & klev           ,&!< in  number of full levels = number of layers
    & klevp1         ,&!< in  number of half levels = number of layer interfaces
    & this_datetime  ,&!< in  actual time step
    & pp_hl          ,&!< in  pressure at half levels at t-dt [Pa]
    & pp_fl          ,&!< in  pressure at full levels at t-dt [Pa]
    & xm_dry         ,&!< in  dry air mass in layer [kg/m2]
    & xm_trc         ,&!< in  tracer  mass in layer [kg/m2]
    & cld_frc        ,&!< in   cloud fraction [m2/m2]
    & xm_ozn         ,&  !< inout ozone mass mixing ratio [kg/kg]
    & xm_vap,         & !< water vapor mass in layer [kg/m2]
    & xm_liq,         & !< cloud water mass in layer [kg/m2]
    & xm_ice,         & !< cloud ice   mass in layer [kg/m2]
    & xm_co2,         & !< CO2 mass in layer [kg/m2]
    & xm_o3,          & !< O3  mass in layer [kg/m2]
    & xm_o2,          & !< O2  mass in layer [kg/m2]
    & xm_ch4,         & !< CH4 mass in layer [kg/m2]
    & xm_n2o,         & !< N2O mass in layer [kg/m2]
    & xm_cfc,         & !< CFC mass in layer [kg/m2]
    & xc_frc,         &
    & cld_cvr  )
     


    INTEGER, INTENT(in)     :: &
    & jg,                      & !< domain index
    & jb,                      & !< block index
    & jcs,                     & !< start index for loop over block
    & jce,                     & !< end   index for loop over block
    & kbdim,                   & !< dimension of block over cells
    & klev,                    & !< number of full levels = number of layers
    & klevp1                     !< number of half levels = number of layer interfaces

    TYPE(datetime), POINTER :: this_datetime !< actual time step

    REAL(wp), INTENT(IN)    :: &
    & pp_hl(kbdim,klevp1),     & !< pressure at half levels [Pa]
    & pp_fl(kbdim,klev),       & !< Pressure at full levels [Pa]
    & xm_dry(kbdim,klev),      & !< dry air mass in layer [kg/m2]
    & xm_trc(kbdim,klev,ntracer),&  !< tracer mass in layer [kg/m2]
    & cld_frc(kbdim,klev)        !< cloud fraction in layer [m2/m2]

    REAL(wp), INTENT(INOUT) :: &
    & xm_ozn(kbdim,klev)         !< ozone mixing ratio  [kg/kg]

    REAL(wp), INTENT(INOUT) ::      &
    & xm_vap(kbdim,klev),           & !< water vapor mass in layer [kg/m2]
    & xm_liq(kbdim,klev),           & !< cloud water mass in layer [kg/m2]
    & xm_ice(kbdim,klev),           & !< cloud ice   mass in layer [kg/m2]
    & xm_co2(kbdim,klev),           & !< CO2 mass in layer [kg/m2]
    & xm_o3(kbdim,klev),            & !< O3  mass in layer [kg/m2]
    & xm_o2(kbdim,klev),            & !< O2  mass in layer [kg/m2]
    & xm_ch4(kbdim,klev),           & !< CH4 mass in layer [kg/m2]
    & xm_n2o(kbdim,klev),           & !< N2O mass in layer [kg/m2]
    & xm_cfc(kbdim,klev,2),         & !< CFC mass in layer [kg/m2]
    & xc_frc (:,:)                 ,&
    & cld_cvr(:)    

    INTEGER             :: jk
    INTEGER             :: jtrc      !< tracer index

    REAL(wp),ALLOCATABLE:: zo3_timint(:,:) !< intermediate value of ozon
    REAL(wp)            :: mmr       !< local mass mixing ratio

    ! vertical profile parameters (vpp) of CH4 and N2O
    REAL(wp), PARAMETER :: vpp_ch4(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp/)
    REAL(wp), PARAMETER :: vpp_n2o(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)

    ! Shortcuts to components of echam_rad_config
    !
    INTEGER , POINTER   :: irad_h2o, irad_co2, irad_ch4, irad_n2o, irad_o3, irad_o2, irad_cfc11, irad_cfc12
    REAL(wp), POINTER   ::            vmr_co2,  vmr_ch4,  vmr_n2o,           vmr_o2,  vmr_cfc11,  vmr_cfc12
    REAL(wp), POINTER   :: frad_h2o, frad_co2, frad_ch4, frad_n2o, frad_o3, frad_o2, frad_cfc

    irad_h2o   => echam_rad_config(jg)% irad_h2o
    irad_co2   => echam_rad_config(jg)% irad_co2
    irad_ch4   => echam_rad_config(jg)% irad_ch4
    irad_n2o   => echam_rad_config(jg)% irad_n2o
    irad_o3    => echam_rad_config(jg)% irad_o3
    irad_o2    => echam_rad_config(jg)% irad_o2
    irad_cfc11 => echam_rad_config(jg)% irad_cfc11
    irad_cfc12 => echam_rad_config(jg)% irad_cfc12
    !
    vmr_co2    => echam_rad_config(jg)% vmr_co2
    vmr_ch4    => echam_rad_config(jg)% vmr_ch4
    vmr_n2o    => echam_rad_config(jg)% vmr_n2o
    vmr_o2     => echam_rad_config(jg)% vmr_o2
    vmr_cfc11  => echam_rad_config(jg)% vmr_cfc11
    vmr_cfc12  => echam_rad_config(jg)% vmr_cfc12
    !
    frad_h2o   => echam_rad_config(jg)% frad_h2o
    frad_co2   => echam_rad_config(jg)% frad_co2
    frad_ch4   => echam_rad_config(jg)% frad_ch4
    frad_n2o   => echam_rad_config(jg)% frad_n2o
    frad_o3    => echam_rad_config(jg)% frad_o3
    frad_o2    => echam_rad_config(jg)% frad_o2
    frad_cfc   => echam_rad_config(jg)% frad_cfc
    !
    ! --- phases of water substance
    !
    !     vapor
    xm_vap(jcs:jce,:) = gas_profile(jcs, jce, klev, irad_h2o, xm_dry,     &
         &                          gas_val      = xm_trc(:,:,iqv),       &
         &                          gas_factor   = frad_h2o)
    !     cloud water
    xm_liq(jcs:jce,:) = gas_profile(jcs, jce, klev, irad_h2o, xm_dry,     &
         &                          gas_val      = xm_trc(:,:,iqc),       &
         &                          gas_epsilon  = 0.0_wp,                &
         &                          gas_factor   = frad_h2o)
    !     cloud ice
    xm_ice(jcs:jce,:) = gas_profile(jcs, jce, klev, irad_h2o, xm_dry,     &
         &                          gas_val      = xm_trc(:,:,iqi),       &
         &                          gas_epsilon  = 0.0_wp,                &
         &                          gas_factor   = frad_h2o)
    !
    ! --- cloud cover
    ! 
    xc_frc(jcs:jce,1:klev) = MERGE(cld_frc(jcs:jce,1:klev), 0._wp, &
         xm_liq(jcs:jce,1:klev) > 0.0_wp .OR. xm_ice(jcs:jce,1:klev) > 0.0_wp)
    !
    cld_cvr(jcs:jce) = 1.0_wp - xc_frc(jcs:jce,1)
    DO jk = 2, klev
      cld_cvr(jcs:jce) = cld_cvr(jcs:jce)                                          &
           &            *(1.0_wp-MAX(xc_frc(jcs:jce,jk),xc_frc(jcs:jce,jk-1)))     &
           &            /(1.0_wp-MIN(xc_frc(jcs:jce,jk-1),1.0_wp-EPSILON(1.0_wp)))
    END DO
    cld_cvr(jcs:jce) = 1.0_wp-cld_cvr(jcs:jce)   
    !
    ! --- gases
    !
    ! CO2: use CO2 tracer only if the CO2 index is in the correct range
    jtrc=MIN(ico2,ntracer)
    mmr = vmr_co2 * amco2/amd
    xm_co2(jcs:jce,:)   = gas_profile(jcs, jce, klev, irad_co2, xm_dry,   &
         &                            gas_mmr      = mmr,                 &
         &                            gas_scenario = ghg_co2mmr,          &
         &                            gas_val      = xm_trc(:,:,jtrc),    &
         &                            gas_factor   = frad_co2)

    mmr = vmr_ch4 * amch4/amd
    xm_ch4(jcs:jce,:)   = gas_profile(jcs, jce, klev, irad_ch4, xm_dry,   &
         &                            gas_mmr      = mmr,                 &
         &                            gas_scenario = ghg_ch4mmr,          &
         &                            pressure = pp_fl, xp = vpp_ch4,     &
         &                            gas_factor   = frad_ch4)

    mmr = vmr_n2o * amn2o/amd
    xm_n2o(jcs:jce,:)   = gas_profile(jcs, jce, klev, irad_n2o, xm_dry,   &
         &                            gas_mmr      = mmr,                 &
         &                            gas_scenario = ghg_n2ommr,          &
         &                            pressure = pp_fl, xp = vpp_n2o,     &
         &                            gas_factor   = frad_n2o)

    mmr = vmr_cfc11 * amc11/amd
    xm_cfc(jcs:jce,:,1) = gas_profile(jcs, jce, klev, irad_cfc11, xm_dry, &
         &                            gas_mmr      = mmr,                 &
         &                            gas_scenario = ghg_cfcmmr(1),       &
         &                            gas_factor   = frad_cfc)

    mmr = vmr_cfc12 * amc12/amd
    xm_cfc(jcs:jce,:,2) = gas_profile(jcs, jce, klev, irad_cfc12, xm_dry, &
         &                            gas_mmr      = mmr,                 &
         &                            gas_scenario = ghg_cfcmmr(2),       &
         &                            gas_factor   = frad_cfc)

    ! O3
    SELECT CASE(irad_o3)
      !
    CASE default
      CALL finish('radiation','o3: this "irad_o3" is not supported')
      !
    CASE(0) ! no ozone 
      xm_o3(jcs:jce,:)    = gas_profile(jcs, jce, klev, irad_o3, xm_dry,  &
           &                            gas_scenario_v = xm_ozn(:,:),     &
           &                            gas_factor     = frad_o3)
      !
    CASE(1) ! ozone tracer
      jtrc=MIN(io3,ntracer)
      xm_o3 (jcs:jce,:)   = gas_profile(jcs, jce, klev, irad_o3, xm_dry,  &
           &                            gas_val        = xm_trc(:,:,jtrc),&
           &                            gas_factor     = frad_o3)
      !
    CASE(4) ! ozone is constant in time
      CALL o3_pl2ml ( jcs=jcs, jce=jce, kbdim = kbdim,       &
           &          nlev_pres = ext_ozone(jg)%nplev_o3,    &
           &          klev = klev,                           &
           &          pfoz = ext_ozone(jg)%plev_full_o3,     &
           &          phoz = ext_ozone(jg)%plev_half_o3,     &
           &          ppf  = pp_fl(:,:),                     &! in  app1
           &          pph  = pp_hl(:,:),                     &! in  aphp1
           &          o3_time_int = ext_ozone(jg)%o3_plev(:,:,jb,1),&! in
           &          o3_clim     = xm_ozn(:,:)              )! OUT
      xm_o3(jcs:jce,:)    = gas_profile(jcs, jce, klev, irad_o3, xm_dry,  &
           &                            gas_scenario_v = xm_ozn(:,:),     &
           &                            gas_factor     = frad_o3)

    CASE(2,8)
      ALLOCATE(zo3_timint(kbdim,ext_ozone(jg)%nplev_o3))
      CALL o3_timeint(jcs=jcs, jce=jce, kbdim = kbdim,       &
           &          nlev_pres=ext_ozone(jg)%nplev_o3,      &
           &          ext_o3=ext_ozone(jg)%o3_plev(:,:,jb,:),&
           &          current_date=this_datetime,            &
           &          o3_time_int=zo3_timint                 )
      CALL o3_pl2ml ( jcs=jcs, jce=jce, kbdim = kbdim,       &
           &          nlev_pres = ext_ozone(jg)%nplev_o3,    &
           &          klev = klev,                           &
           &          pfoz = ext_ozone(jg)%plev_full_o3,     &
           &          phoz = ext_ozone(jg)%plev_half_o3,     &
           &          ppf  = pp_fl(:,:),                     &
           &          pph  = pp_hl(:,:),                     &
           &          o3_time_int = zo3_timint,              &
           &          o3_clim     = xm_ozn(:,:)              )
      DEALLOCATE(zo3_timint)
      xm_o3(jcs:jce,:)    = gas_profile(jcs, jce, klev, irad_o3, xm_dry,  &
           &                             gas_scenario_v = xm_ozn(:,:),    &
           &                             gas_factor     = frad_o3)
      !
    CASE(10) ! ozone from ART
      IF(.NOT. lart) CALL finish('psrad:mo_psrad_radiation', &
        & 'irad_o3=10 not supported without lart = .True.'   )
      !
      xm_o3(jcs:jce,:)    = gas_profile(jcs, jce, klev, 8, xm_dry,        &
        &                                gas_scenario_v = xm_ozn(:,:),    &
        &                                gas_factor     = frad_o3)

    END SELECT

    mmr = vmr_o2 * amo2/amd
    xm_o2(jcs:jce,:)    = gas_profile(jcs, jce, klev, irad_o2, xm_dry,    &
         &                            gas_mmr      = mmr,                 &
         &                            gas_factor   = frad_o2)

  END SUBROUTINE psrad_get_gas_profiles
  !-------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !>
  !! GAS_PROFILE:  Determines Gas distributions based on case specification
  !!
  !!               Input  units: mass mixing ratio with respect to dry air [kg/kg]
  !!               Output units: mass content in layer [kg/m2]
  !! 
  !! @par Revsision History 
  !! B. Stevens   (2009-08). 
  !! H. Schmidt   (2010-08): profile calculation added for scenario case.
  !! M. Giorgetta (2016-09): change output units from kg/kg to kg/m2
  !!
  !! Description:
  !! This routine calculates the gas distributions for one of five cases:
  !! - igas=0: no gas present;
  !! - igas=1: prognostic gas;
  !! - igas=2: specified mixing ratio;
  !! - igas=3: mixing ratio decaying with height given profile;
  !! - igas=4: scenario run with different mixing ratio, if profile parameters
  !!           are given a vertical profile is calculated as for igas=3.
  !
  FUNCTION gas_profile (jcs, jce, klev, igas, xm_dry,           &
       &                gas_mmr, gas_scenario, gas_mmr_v,       &
       &                gas_scenario_v, gas_val, xp, pressure,  &
       &                gas_epsilon, gas_factor)

    INTEGER,             INTENT (IN) :: jcs                  ! start index
    INTEGER,             INTENT (IN) :: jce                  ! end   index and horizontal dimension
    INTEGER,             INTENT (IN) :: klev                 ! vertical dimensions
    INTEGER,             INTENT (IN) :: igas                 ! gas case
    REAL (wp),           INTENT (IN) :: xm_dry(:,:)          ! dry air content    [kg/m2]
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_mmr              ! for igas = 2 and 3 [kg/kg]
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_scenario         ! for igas = 4       [kg/kg]
    REAL (wp), OPTIONAL, INTENT (IN) :: pressure(:,:), xp(3) ! for igas = 3 and 4 [kg/kg]
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_mmr_v(:,:)       ! for igas = 2       [kg/kg]
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_scenario_v(:,:)  ! for igas = 4       [kg/kg]
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_val(:,:)         ! for igas = 1       [kg/m2]
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_epsilon
    REAL (wp), OPTIONAL, INTENT (IN) :: gas_factor

    REAL (wp) :: gas_profile(jcs:jce,klev), zx_d, zx_m, eps, fgas
    LOGICAL   :: gas_initialized

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
      gas_profile(jcs:jce,:) = 0.0_wp
      gas_initialized = .TRUE.

    CASE (1)                             ! 1: horizontally and vertically variable
      IF (PRESENT(gas_val)) THEN
        gas_profile(jcs:jce,:) = gas_val(jcs:jce,:)
        gas_initialized = .TRUE.
      END IF

    CASE (2)
      IF (PRESENT(gas_mmr)) THEN         ! 2a: horizontally and vertically constant
        gas_profile(jcs:jce,:) = gas_mmr  * xm_dry(jcs:jce,:)
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_mmr_v)) THEN  ! 2b: = (1)
        gas_profile(jcs:jce,:) = gas_mmr_v(jcs:jce,:)  * xm_dry(jcs:jce,:)
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_scenario_v)) THEN ! 2c: = (1)
        gas_profile(jcs:jce,:) = gas_scenario_v(jcs:jce,:) * xm_dry(jcs:jce,:)
        gas_initialized = .TRUE.
      END IF

    CASE (3)                             ! 3: horizontally constant and tanh-profile in the vertical
      IF (PRESENT(gas_mmr) .AND. PRESENT(xp) .AND. PRESENT(pressure)) THEN
        zx_m = (gas_mmr+xp(1)*gas_mmr)*0.5_wp
        zx_d = (gas_mmr-xp(1)*gas_mmr)*0.5_wp
        gas_profile(jcs:jce,:)=(1-(zx_d/zx_m)*TANH(LOG(pressure(jcs:jce,:)   &
             &                  /xp(2)) /xp(3))) * zx_m * xm_dry(jcs:jce,:)
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
          gas_profile(jcs:jce,:)=(1-(zx_d/zx_m)*TANH(LOG(pressure(jcs:jce,:)   &
             &                    /xp(2)) /xp(3))) * zx_m * xm_dry(jcs:jce,:)
        ELSE                                          ! 4b: = (2a)
          gas_profile(jcs:jce,:)=gas_scenario * xm_dry(jcs:jce,:)
        ENDIF
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_scenario_v)) THEN          ! 4c: = (1)
        gas_profile(jcs:jce,:) = gas_scenario_v(jcs:jce,:) * xm_dry(jcs:jce,:)
        gas_initialized = .TRUE.
      END IF

    END SELECT

    IF (.NOT. gas_initialized) &
         CALL finish('radiation','gas_profile options not supported')

    gas_profile(jcs:jce,:) = MAX(fgas * gas_profile(jcs:jce,:),eps)
    
  END FUNCTION gas_profile
  !---------------------------------------------------------------------------


  END MODULE mo_psrad_radiation
