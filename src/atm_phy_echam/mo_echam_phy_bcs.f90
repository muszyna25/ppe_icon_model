!>
!! @brief Prepares boundary conditions needed for ECHAM physics
!!
!! @author Marco Giorgetta (MPI-M)
!!
!! @par Revision History
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

MODULE mo_echam_phy_bcs

  USE mo_kind                       ,ONLY: wp, i8
  USE mtime,                         ONLY: datetime, newDatetime, operator(+), &
       &                                   newTimedelta, timedelta, max_timedelta_str_len, &
       &                                   getPTStringFromSeconds, getNoOfSecondsElapsedInDayDateTime
  USE mo_model_domain               ,ONLY: t_patch

  USE mo_master_config              ,ONLY: isRestart
  USE mo_echam_phy_config           ,ONLY: echam_phy_config
  USE mo_radiation_config           ,ONLY: ighg, isolrad, tsi, tsi_radt, ssi_radt, irad_o3, irad_aero

  USE mo_echam_sfc_indices          ,ONLY: nsfc_type, iwtr, iice

  USE mo_impl_constants             ,ONLY: io3_amip

  USE mo_bcs_time_interpolation,     ONLY: t_time_interpolation_weights,         &
       &                                   calculate_time_interpolation_weights
  
  USE mo_bc_greenhouse_gases        ,ONLY: bc_greenhouse_gases_time_interpolation
  USE mo_bc_sst_sic                 ,ONLY: get_current_bc_sst_sic_year, read_bc_sst_sic, &
    &                                      bc_sst_sic_time_interpolation
  USE mo_bc_solar_irradiance        ,ONLY: read_bc_solar_irradiance, ssi_time_interpolation
  USE mo_psrad_radiation            ,ONLY: pre_psrad_radiation
  USE mo_bc_ozone                   ,ONLY: read_bc_ozone
  USE mo_bc_aeropt_kinne            ,ONLY: read_bc_aeropt_kinne
  USE mo_bc_aeropt_stenchikov       ,ONLY: read_bc_aeropt_stenchikov

  USE mo_echam_phy_memory           ,ONLY: prm_field

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_bcs_global
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_bcs'

CONTAINS
  !>
  !! SUBROUTINE echam_phy_bcs_global
  !!
  !! This subroutine is called in the time loop of the model and serves
  !! to prepare the boundary conditions for the ECHAM physics.
  !!
  !! This routine contains all time dependent preparations that are valid
  !! for the whole globe and need to be done outside of the loop over
  !! rows, in which the echam_phy_main routine is called for each row
  !! of the block.
  !!
  !! Note that each call of this subroutine deals with a single grid
  !! with index jg rather than the entire grid tree.

  SUBROUTINE echam_phy_bcs_global( mtime_old,    &
    &                              jg           ,&! in
    &                              patch        ,&! in
    &                              dtadv_loc    ,&! in
    &                              ltrig_rad    ) ! out

    ! Arguments

    TYPE(datetime) , POINTER  ,INTENT(in)    :: mtime_old
    INTEGER                   ,INTENT(in)    :: jg             !< grid index
    TYPE(t_patch)  , TARGET   ,INTENT(in)    :: patch          !< description of grid jg
    REAL(wp)                  ,INTENT(in)    :: dtadv_loc      !< timestep of advection and physics on grid jg
    LOGICAL                   ,INTENT(out)   :: ltrig_rad      !< trigger for radiation transfer computation

    ! Local variables

    TYPE(datetime) , POINTER                 :: radiation_time !< date and time for radiative transfer
    TYPE(timedelta), POINTER                 :: td_radiation_offset
    CHARACTER(len=max_timedelta_str_len)     :: dstring

    REAL(wp)                                 :: dsec           !< [s] time increment of datetime_radtran wrt. datetime

    LOGICAL                                  :: is_1st_call = .TRUE.

    TYPE(t_time_interpolation_weights), SAVE :: current_time_interpolation_weights
    TYPE(t_time_interpolation_weights), SAVE :: radiation_time_interpolation_weights 

    ! Local parameters

!!$    CHARACTER(*), PARAMETER :: method_name = "echam_phy_bcs_global"

    !-------------------------------------------------------------------------
    ! Prepare some global parameters or parameter arrays
    !-------------------------------------------------------------------------

    ! Check whether the radiative transfer needs to be calculated in this
    ! timestep. Then boundary conditions must be prepared for this purpose.
    !
    IF (echam_phy_config%lrad) THEN
      ltrig_rad   = ( is_1st_call .AND. (.NOT.isRestart()) ) .OR. &
        &           ( MOD(getNoOfSecondsElapsedInDayDateTime(mtime_old),NINT(echam_phy_config%dt_rad)) == 0 )
    ELSE
      ltrig_rad = .FALSE.
    END IF

    ! Set the time instance datetime_radtran for the zenith angle to be used
    ! in the radiative transfer. All other input for the radiative transfer
    ! is for datetime, i.e. the start date and time of the current timestep.
    !

    radiation_time => newDatetime(mtime_old)
    IF (ltrig_rad) THEN
      dsec = 0.5_wp*(echam_phy_config%dt_rad - dtadv_loc) ! [s] time increment for zenith angle
      CALL getPTStringFromSeconds(dsec, dstring)
      td_radiation_offset => newTimedelta(dstring)
      radiation_time = radiation_time + td_radiation_offset
    END IF

    ! interpolation weights for linear interpolation
    ! of monthly means onto the actual integration time step
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_old)

    ! Read and interpolate in time monthly mean SST for AMIP simulations
    ! SST is needed for turbulent vertical fluxes and for radiation.
    !
    IF (echam_phy_config%lamip) THEN
      IF (iwtr <= nsfc_type .OR. iice <= nsfc_type) THEN
        IF (mtime_old%date%year /= get_current_bc_sst_sic_year()) THEN
          CALL read_bc_sst_sic(mtime_old%date%year, patch)
        END IF
        CALL bc_sst_sic_time_interpolation(current_time_interpolation_weights, &
          &                                prm_field(jg)%lsmask(:,:)         , &
          &                                prm_field(jg)%ts_tile(:,:,iwtr)   , &
          &                                prm_field(jg)%seaice(:,:)         , &
          &                                prm_field(jg)%siced(:,:)          , &
          &                                patch                              )

        ! The ice model should be able to handle different thickness classes, 
        ! but for AMIP we ONLY USE one ice class.
        IF (iice <= nsfc_type) THEN
          prm_field(jg)%conc(:,1,:) = prm_field(jg)%seaice(:,:)
          prm_field(jg)%hi  (:,1,:) = prm_field(jg)%siced (:,:)
        END IF
      END IF
    END IF

    ! total solar irradiation at the mean sun earth distance
    IF (isolrad==1) THEN
      CALL read_bc_solar_irradiance(mtime_old%date%year, .FALSE.)
      CALL ssi_time_interpolation(current_time_interpolation_weights, .FALSE., tsi)
    END IF

    ! quantities needed for the radiative transfer only
    !
    IF (ltrig_rad) THEN
      !
      ! interpolation weights for linear interpolation
      ! of monthly means onto the radiation time step
      radiation_time_interpolation_weights = calculate_time_interpolation_weights(radiation_time)
      !
      ! total and spectral solar irradiation at the mean sun earth distance
      IF (isolrad==1) THEN
        CALL read_bc_solar_irradiance(mtime_old%date%year,.TRUE.)
        CALL ssi_time_interpolation(radiation_time_interpolation_weights,.TRUE.,tsi_radt,ssi_radt)
      END IF
      !
      ! greenhouse gas concentrations, assumed constant in horizontal dimensions
      IF (ighg > 0) THEN
        CALL bc_greenhouse_gases_time_interpolation(mtime_old)
      END IF
      !
      ! ozone concentration
      IF (irad_o3 == io3_amip) THEN
        CALL read_bc_ozone(mtime_old%date%year, patch)
      END IF
      !
      ! tropospheric aerosol optical properties
      IF (irad_aero == 13) THEN
        CALL read_bc_aeropt_kinne(mtime_old%date%year, patch)
      END IF
      !
      ! stratospheric aerosol optical properties
      IF (irad_aero == 14) THEN
        CALL read_bc_aeropt_stenchikov(mtime_old)
      END IF
      !
      ! tropospheric and stratospheric aerosol optical properties
      IF (irad_aero == 15) THEN
        CALL read_bc_aeropt_kinne     (mtime_old%date%year, patch)
        CALL read_bc_aeropt_stenchikov(mtime_old)
      END IF
      ! tropospheric background aerosols (Kinne) and stratospheric
      ! aerosols (Stenchikov) + simple plumes (analytical, nothing to be read
      ! here, initialization see init_echam_phy (mo_echam_phy_init)) 
      IF (irad_aero == 18) THEN
        CALL read_bc_aeropt_kinne     (1850_i8, patch)
        CALL read_bc_aeropt_stenchikov(mtime_old)
      END IF
      !
    END IF ! ltrig_rad

    CALL pre_psrad_radiation( &
            & patch,                           radiation_time,               &
            & mtime_old,                       ltrig_rad,                    &
            & prm_field(jg)%cosmu0,            prm_field(jg)%daylght_frc,    &
            & prm_field(jg)%cosmu0_rt,         prm_field(jg)%daylght_frc_rt )

    is_1st_call = .FALSE.

  END SUBROUTINE echam_phy_bcs_global

END MODULE mo_echam_phy_bcs
