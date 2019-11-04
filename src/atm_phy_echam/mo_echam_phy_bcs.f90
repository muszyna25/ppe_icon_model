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

  USE mo_kind                       ,ONLY: wp
  USE mo_exception                  ,ONLY: warning
  USE mtime                         ,ONLY: datetime , newDatetime ,                        &
       &                                   timedelta, newTimedelta, max_timedelta_str_len, &
       &                                   operator(+), operator(-), operator(*),          &
       &                                   operator(<=),operator(>),                       &   
       &                                   getPTStringFromSeconds,                         &
       &                                   getTotalSecondsTimeDelta,                       &
       &                                   isCurrentEventActive, deallocateDatetime
  USE mo_model_domain               ,ONLY: t_patch, p_patch
  USE mo_impl_constants             ,ONLY: max_dom

  USE mo_echam_phy_memory           ,ONLY: t_echam_phy_field, prm_field
  USE mo_echam_phy_config           ,ONLY: echam_phy_config, echam_phy_tc, dt_zero
  USE mo_echam_rad_config           ,ONLY: echam_rad_config
  USE mo_ccycle_config              ,ONLY: ccycle_config
  USE mo_psrad_solar_data           ,ONLY: ssi_radt, tsi_radt, tsi
  USE mo_psrad_radiation            ,ONLY: pre_psrad_radiation

  USE mo_echam_sfc_indices          ,ONLY: nsfc_type, iwtr, iice

  USE mo_bcs_time_interpolation     ,ONLY: t_time_interpolation_weights,         &
       &                                   calculate_time_interpolation_weights
  
  USE mo_bc_greenhouse_gases        ,ONLY: bc_greenhouse_gases_time_interpolation
  USE mo_bc_sst_sic                 ,ONLY: get_current_bc_sst_sic_year, read_bc_sst_sic, &
    &                                      bc_sst_sic_time_interpolation
  USE mo_bc_solar_irradiance        ,ONLY: read_bc_solar_irradiance, ssi_time_interpolation
  USE mo_bc_ozone                   ,ONLY: read_bc_ozone
  USE mo_bc_aeropt_kinne            ,ONLY: read_bc_aeropt_kinne
  USE mo_bc_aeropt_stenchikov       ,ONLY: read_bc_aeropt_stenchikov
  USE mo_atmo_psrad_interface       ,ONLY: dtrad_shift
#if defined( _OPENACC )
  USE mo_mpi                   ,ONLY: i_am_accel_node, my_process_is_work
#endif

  ! for 6hourly sst and ice data
  USE mo_time_config,          ONLY: time_config
  USE mo_reader_sst_sic,       ONLY: t_sst_sic_reader
  USE mo_interpolate_time,     ONLY: t_time_intp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_bcs
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_bcs'

  TYPE(t_sst_sic_reader), TARGET :: sst_sic_reader
  TYPE(t_time_intp)      :: sst_intp
  TYPE(t_time_intp)      :: sic_intp
  REAL(wp), ALLOCATABLE  :: sst_dat(:,:,:,:)
  REAL(wp), ALLOCATABLE  :: sic_dat(:,:,:,:)

CONTAINS
  !>
  !! SUBROUTINE echam_phy_bcs
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

  SUBROUTINE echam_phy_bcs( patch        ,&! in
    &                       mtime_old,    &
    &                       dtadv_loc    ) ! out

    ! Arguments

    TYPE(t_patch)  , TARGET   ,INTENT(in)    :: patch          !< description of this grid
    TYPE(datetime) , POINTER  ,INTENT(in)    :: mtime_old
    REAL(wp)                  ,INTENT(in)    :: dtadv_loc      !< timestep of advection and physics on this grid

    ! Local variables

    ! mtime currently does not work with arrays of datetime pointers
    ! therefore a type is constructed around the mtime pointer
    TYPE t_radtime_domains
      TYPE(datetime) , POINTER               :: radiation_time => NULL() !< date and time for radiative transfer
    END TYPE t_radtime_domains
    !
    TYPE(t_radtime_domains), SAVE            :: radtime_domains(max_dom)

    TYPE(timedelta), POINTER                 :: td_radiation_offset
    CHARACTER(len=max_timedelta_str_len)     :: dstring

    REAL(wp)                                 :: dsec           !< [s] time increment of datetime_radtran wrt. datetime
    REAL(wp)                                 :: dtrad_loc      !< [s] time local radiation time step

    LOGICAL                                  :: luse_rad       !< use LW radiation
    LOGICAL                                  :: ltrig_rad      !< trigger for LW radiative transfer computation

    LOGICAL                                  :: ghg_time_interpol_already_done

    TYPE(t_time_interpolation_weights), SAVE :: current_time_interpolation_weights
    TYPE(t_time_interpolation_weights), SAVE :: radiation_time_interpolation_weights 

    LOGICAL, ALLOCATABLE                     :: mask_sftof(:,:)

!!$    CHARACTER(*), PARAMETER :: method_name = "echam_phy_bcs"

    ! Shortcuts to components of echam_cld_config
    !
    INTEGER          :: jc, jb, jg
    INTEGER, POINTER :: ighg, isolrad, irad_o3, irad_aero
    TYPE(t_echam_phy_field) , POINTER    :: field
    !
#ifdef _OPENACC
    i_am_accel_node = .FALSE.    ! Deactivate GPUs
#endif
    !
    jg        =  patch%id ! grid index
    
    !-------------------------------------------------------------------------
    ! Prepare some global parameters or parameter arrays
    !-------------------------------------------------------------------------

    ! interpolation weights for linear interpolation
    ! of monthly means onto the actual integration time step
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_old)

    ! Read and interpolate in time monthly mean SST for AMIP simulations
    ! SST is needed for turbulent vertical fluxes and for radiation.
    !
    IF (echam_phy_tc(jg)%dt_rad > dt_zero .OR. echam_phy_tc(jg)%dt_vdf > dt_zero) THEN
      !
      IF (echam_phy_config(jg)%lamip) THEN
       field => prm_field(jg)
       IF (iwtr <= nsfc_type .OR. iice <= nsfc_type) THEN
        !
        IF (.NOT.  ANY(echam_phy_config(:)%lsstice)) THEN
          IF (mtime_old%date%year /= get_current_bc_sst_sic_year()) THEN
            CALL read_bc_sst_sic(mtime_old%date%year, patch)
          END IF
          ALLOCATE(mask_sftof(SIZE(field%sftof,1), SIZE(field%sftof,2)))
          !$ACC DATA PRESENT( field%sftof ) &
          !$ACC       CREATE( mask_sftof )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jb = 1, SIZE(field%sftof,2)
            !$ACC LOOP GANG VECTOR
            DO jc = 1, SIZE(field%sftof,1)
              mask_sftof(jc, jb) = field%sftof(jc,jb) > 0._wp
            END DO
          END DO
          !$ACC END PARALLEL
          CALL bc_sst_sic_time_interpolation(current_time_interpolation_weights , &
            &                                field%ts_tile(:,:,iwtr)            , &
            &                                field%seaice (:,:)                 , &
            &                                field%siced  (:,:)                 , &
            &                                patch                              , &
            &                                mask_sftof                         , &
            &                                .FALSE., lopenacc = .TRUE. )
          !$ACC END DATA
          DEALLOCATE(mask_sftof)

        ELSE
          !
          ! READ 6-hourly sst values (dyamond+- setup, preliminary)
          CALL sst_sic_reader%init(p_patch(jg), 'sst-sic-runmean_G.nc')
          CALL sst_intp%init(sst_sic_reader, time_config%tc_current_date, "SST")
          CALL sst_intp%intp(time_config%tc_current_date, sst_dat)
          !$ACC DATA PRESENT( field%ts_tile ) &
          !$ACC      CREATE( sst_dat )
          !$ACC UPDATE DEVICE( sst_dat )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jb = LBOUND(field%ts_tile, 2), UBOUND(field%ts_tile,2)
            !$ACC LOOP GANG VECTOR
            DO jc = LBOUND(field%ts_tile, 1), UBOUND(field%ts_tile,1)
              IF (sst_dat(jc,1,jb,1) > 0.0_wp) THEN
                prm_field(jg)%ts_tile(jc,jb,iwtr) = sst_dat(jc,1,jb,1)
              END IF
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
          !
          CALL sic_intp%init(sst_sic_reader, time_config%tc_current_date, "SIC")
          CALL sic_intp%intp(time_config%tc_current_date, sic_dat)
          !$ACC DATA PRESENT( field%seaice, sic_dat )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jb = LBOUND(field%seaice, 2), UBOUND(field%seaice,2)
            !$ACC LOOP GANG VECTOR
            DO jc = LBOUND(field%seaice, 1), UBOUND(field%seaice,1)
              field%seaice(jc,jb) = sic_dat(jc,1,jb,1)   !160:164
              field%seaice(jc,jb) = MERGE(0.99_wp, field%seaice(jc,jb),  &
                                              & field%seaice(jc,jb) > 0.99_wp)
              field%seaice(jc,jb) = MERGE(0.0_wp, field%seaice(jc,jb),  &
                                              & field%seaice(jc,jb) <= 0.01_wp)
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA

          ! set ice thickness
          !$ACC DATA PRESENT( field%seaice, field%siced, sic_dat, p_patch(jg)%cells%center )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jb = LBOUND(field%siced, 2), UBOUND(field%siced,2)
            !$ACC LOOP GANG VECTOR
            DO jc = LBOUND(field%siced, 1), UBOUND(field%siced,1)
              IF (field%seaice(jc,jb) > 0.0_wp) THEN
                field%siced(jc,jb) = MERGE(2.0_wp, 1.0_wp, p_patch(jg)%cells%center(jc,jb)%lat > 0.0_wp)
              ELSE
                field%siced(jc,jb) = 0.0_wp
              END IF
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
        !
        END IF
       END IF
        ! The ice model should be able to handle different thickness classes, 
        ! but for AMIP we only use one ice class.
        IF (iice <= nsfc_type) THEN
          !$ACC DATA PRESENT( field%conc, field%seaice, field%hi, field%siced )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jb = LBOUND(field%conc,3), UBOUND(field%conc,3)
            !$ACC LOOP GANG VECTOR
            DO jc = LBOUND(field%conc,1), UBOUND(field%conc,1)
              field%conc(jc,1,jb) = field%seaice(jc,jb)
              field%hi  (jc,1,jb) = field%siced (jc,jb)
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
        END IF
      END IF
      !
    END IF

    ! IF radiation is used in this experiment (dt_rad>0):
    ! - luse_rad : Check whether radiative heating is used in this timestep
    ! - ltrig_rad: Check whether radiative transfer needs to be calculated
    IF (echam_phy_tc(jg)%dt_rad >  dt_zero) THEN
       luse_rad  = (echam_phy_tc(jg)%sd_rad <= mtime_old) .AND. &
            &      (echam_phy_tc(jg)%ed_rad >  mtime_old)
       ltrig_rad = isCurrentEventActive(echam_phy_tc(jg)%ev_rad,mtime_old)
    ELSE
       luse_rad  = .FALSE.
       ltrig_rad = .FALSE.
    END IF

    IF (luse_rad) THEN

      ! total solar irradiation at the mean sun earth distance
      IF (echam_rad_config(jg)% isolrad == 1) THEN
        CALL read_bc_solar_irradiance(mtime_old%date%year, .FALSE.)
        CALL ssi_time_interpolation(current_time_interpolation_weights, .FALSE., tsi)
      END IF

    !
    ! quantities needed for the radiative transfer only
    !
    IF (ltrig_rad) THEN
      !
      ! Set the time instance datetime_radtran for the zenith angle to be used
      ! in the radiative transfer. All other input for the radiative transfer
      ! is for datetime, i.e. the start date and time of the current timestep.
      !
      IF (ASSOCIATED(radtime_domains(jg)%radiation_time)) &
        & CALL deallocateDatetime(radtime_domains(jg)%radiation_time) 
      radtime_domains(jg)%radiation_time => newDatetime(mtime_old)
      dtrad_loc = getTotalSecondsTimeDelta(echam_phy_tc(jg)%dt_rad,mtime_old) ! [s] local time step of radiation
      dsec = 0.5_wp*(dtrad_loc - dtadv_loc) + dtrad_shift                     ! [s] time increment for zenith angle
      CALL getPTStringFromSeconds(dsec, dstring)
      td_radiation_offset => newTimedelta(dstring)
      radtime_domains(jg)%radiation_time = radtime_domains(jg)%radiation_time + td_radiation_offset
      !
      ! interpolation weights for linear interpolation
      ! of monthly means onto the radiation time step
      radiation_time_interpolation_weights = calculate_time_interpolation_weights(radtime_domains(jg)%radiation_time)
      !
      ! total and spectral solar irradiation at the mean sun earth distance
      IF (echam_rad_config(jg)% isolrad == 1) THEN
        CALL read_bc_solar_irradiance(mtime_old%date%year,.TRUE.)
        CALL ssi_time_interpolation(radiation_time_interpolation_weights,.TRUE.,tsi_radt,ssi_radt)
      END IF
      !
      ! ozone concentration
      IF   (      echam_rad_config(jg)% irad_o3 ==  2 &       ! climatological annual cycle defined by monthly data
           & .OR. echam_rad_config(jg)% irad_o3 ==  4 &       ! constant in time
           & .OR. echam_rad_config(jg)% irad_o3 ==  8 &       ! transient monthly means
           & .OR. echam_rad_config(jg)% irad_o3 == 10 ) THEN  ! coupled to ART
        CALL read_bc_ozone(mtime_old%date%year, patch)
      END IF
      !
      ! tropospheric aerosol optical properties
      IF (echam_rad_config(jg)% irad_aero == 13) THEN
        CALL read_bc_aeropt_kinne(mtime_old, patch)
      END IF
      !
      ! stratospheric aerosol optical properties
      IF (echam_rad_config(jg)% irad_aero == 14) THEN
        CALL read_bc_aeropt_stenchikov(mtime_old, patch)
      END IF
      !
      ! tropospheric and stratospheric aerosol optical properties
      IF (echam_rad_config(jg)% irad_aero == 15) THEN
        CALL read_bc_aeropt_kinne     (mtime_old, patch)
        CALL read_bc_aeropt_stenchikov(mtime_old, patch)
      END IF
      !
      ! tropospheric background aerosols (Kinne) and stratospheric
      ! aerosols (Stenchikov) + simple plumes (analytical, nothing to be read
      ! here, initialization see init_echam_phy (mo_echam_phy_init)) 
      IF (echam_rad_config(jg)% irad_aero == 18) THEN
        CALL read_bc_aeropt_kinne     (mtime_old, patch)
        CALL read_bc_aeropt_stenchikov(mtime_old, patch)
      END IF
      !
      ! greenhouse gas concentrations, assumed constant in horizontal dimensions
      ghg_time_interpol_already_done = .FALSE.
      IF  ( echam_rad_config(jg)%irad_co2   == 4 .OR. &
          & echam_rad_config(jg)%irad_ch4   == 4 .OR. &
          & echam_rad_config(jg)%irad_n2o   == 4 .OR. &
          & echam_rad_config(jg)%irad_cfc11 == 4 .OR. &
          & echam_rad_config(jg)%irad_cfc12 == 4      ) THEN
        CALL bc_greenhouse_gases_time_interpolation(mtime_old)
        ghg_time_interpol_already_done = .TRUE.
      END IF
      !
    END IF ! ltrig_rad
    !
    ! co2 concentration for carbon cycle
    IF  ( ccycle_config(jg)%iccycle  == 2 .AND. &       ! c-cycle is used with prescribed co2 conc.
        & ccycle_config(jg)%ico2conc == 4 .AND. &       ! co2 conc. is read from scenario file
        & .NOT. ghg_time_interpol_already_done  ) THEN  ! time interpolation still to be done
      CALL bc_greenhouse_gases_time_interpolation(mtime_old)
    END IF

    IF ( luse_rad ) THEN
       CALL pre_psrad_radiation(                                             &
            & patch,                     radtime_domains(jg)%radiation_time, &
            & mtime_old,                 ltrig_rad,                          &
            & prm_field(jg)%cosmu0,      prm_field(jg)%daylght_frc,          &
            & prm_field(jg)%cosmu0_rt,   prm_field(jg)%daylght_frc_rt )
#ifdef _OPENACC
       CALL warning("GPU:echam_phy_bcs", "GPU device synchronization")
#endif
       !$ACC UPDATE DEVICE( prm_field(jg)%cosmu0, prm_field(jg)%cosmu0_rt,         &
       !$ACC                prm_field(jg)%daylght_frc, prm_field(jg)%daylght_frc_rt )
    END IF

    END IF ! luse_rad

#ifdef _OPENACC
    i_am_accel_node = my_process_is_work()    ! Activate GPUs
#endif

  END SUBROUTINE echam_phy_bcs

END MODULE mo_echam_phy_bcs
