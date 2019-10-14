!>
!! The radiation scheme ecRad expects information from the host model (i.e., ICON) to be
!!   copied to different ecRad data structures: t_ecrad_single_level_type, t_ecrad_gas_type,
!!   t_ecrad_thermodynamics_type and t_ecrad_cloud_type. Similarly, the output of ecRad, i.e.
!!   the radiative fluxes, are stored in a data structure t_ecrad_flux_type. 
!! This module offers subroutines that transfer the data from ICON variables into the
!!   correct ecRad data structure. The subroutines are written in a way that they can be used
!!   on the reduced radiation grid as well as on the full radiation grid. This ensures 
!!   consistency between the two modes.
!!
!! @author Daniel Rieger, Deutscher Wetterdienst, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-10)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_nwp_ecrad_utilities

  USE mo_kind,                   ONLY: wp
  USE mo_math_constants,         ONLY: rad2deg
  USE mo_exception,              ONLY: finish
  USE mo_impl_constants,         ONLY: MAX_CHAR_LENGTH
  USE mo_math_types,             ONLY: t_geographical_coordinates
  USE mo_physical_constants,     ONLY: rd
  USE mo_radiation_config,       ONLY: vmr_co2, vmr_n2o, vmr_o2, vmr_ch4,        &
                                   &   vmr_cfc11, vmr_cfc12,                     &
                                   &   irad_h2o, irad_o3, irad_co2,              &
                                   &   irad_n2o, irad_ch4,                       &
                                   &   irad_o2, irad_cfc11, irad_cfc12,          &
                                   &   vpp_ch4, vpp_n2o
  USE mtime,                     ONLY: datetime
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: ecrad_set_gas_units,                      &
                                   &   t_ecrad_conf,                             &
                                   &   t_ecrad_single_level_type,                &
                                   &   t_ecrad_thermodynamics_type,              &
                                   &   t_ecrad_gas_type, t_ecrad_flux_type,      &
                                   &   t_ecrad_cloud_type,                       &
                                   &   IMassMixingRatio, IVolumeMixingRatio,     &
                                   &   ecRad_IH2O, ecRad_ICO2, ecRad_IO3,        &
                                   &   ecRad_IN2O, ecRad_ICH4,                   &
                                   &   ecRad_IO2, ecRad_ICFC11, ecRad_ICFC12,    &
                                   &   nweight_par_ecrad, iband_par_ecrad,       &
                                   &   weight_par_ecrad
#endif


  IMPLICIT NONE

  PRIVATE
#ifdef __ECRAD
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_utilities'


  PUBLIC :: ecrad_set_single_level
  PUBLIC :: ecrad_set_thermodynamics
  PUBLIC :: ecrad_set_clouds
  PUBLIC :: ecrad_set_gas
  PUBLIC :: ecrad_store_fluxes


CONTAINS


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_single_level:
  !! Set ecRad single level information, i.e. fill t_ecrad_single_level_type with information
  !! from ICON.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-10)
  !!
  SUBROUTINE ecrad_set_single_level(ecrad_single_level, current_datetime, cell_center, cosmu0, tsfc, &
    &                               albvisdif, albnirdif, albvisdir, albnirdir, emis_rad, i_startidx, i_endidx)

    TYPE(t_ecrad_single_level_type), INTENT(inout) :: &
      &  ecrad_single_level       !< ecRad single level information
    TYPE(datetime), INTENT(in) :: &
      &  current_datetime         !< Current date and time
    TYPE(t_geographical_coordinates), INTENT(in) :: &
      &  cell_center(:)           !< lon/lat information of cell centers
    REAL(wp), INTENT(in)     :: &
      &  cosmu0(:),             & !< Cosine of solar zenith angle
      &  tsfc(:),               & !< Surface temperature
      &  albvisdif(:),          & !< Surface albedo for visible range (diffuse)
      &  albnirdif(:),          & !< Surface albedo for near IR range (diffuse)
      &  albvisdir(:),          & !< Surface albedo for visible range (direct)
      &  albnirdir(:),          & !< Surface albedo for near IR range (direct)
      &  emis_rad(:)              !< Longwave surface emissivity
    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
! Local variables
    INTEGER                  :: &
      &  jc                       !< loop index

    DO jc = i_startidx, i_endidx
        ecrad_single_level%cos_sza(jc)            = cosmu0(jc)
        ecrad_single_level%skin_temperature(jc)   = tsfc(jc)
        ecrad_single_level%sw_albedo(jc,1)        = albvisdif(jc)
        ecrad_single_level%sw_albedo(jc,2)        = albnirdif(jc)
        ecrad_single_level%sw_albedo_direct(jc,1) = albvisdir(jc)
        ecrad_single_level%sw_albedo_direct(jc,2) = albnirdir(jc)
        ecrad_single_level%lw_emissivity(jc,1)    = emis_rad(jc)
        ecrad_single_level%iseed(jc)              = create_rdm_seed(cell_center(jc)%lon, &
          &                                                         cell_center(jc)%lat, &
          &                                                         current_datetime)
    ENDDO !jc

  END SUBROUTINE ecrad_set_single_level
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_thermodynamics:
  !! Set ecRad thermodynamics information, i.e. fill t_ecrad_thermodynamics_type with information
  !! from ICON.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-10)
  !!
  SUBROUTINE ecrad_set_thermodynamics(ecrad_thermodynamics, temp, pres, pres_ifc, tsfc, &
    &                                 nlev, nlevp1, i_startidx, i_endidx)

    TYPE(t_ecrad_thermodynamics_type), INTENT(inout) :: &
      &  ecrad_thermodynamics     !< ecRad thermodynamics information
    REAL(wp), INTENT(in)     :: &
      &  temp(:,:),             & !< Full level temperature field
      &  pres(:,:),             & !< Full level pressure field
      &  pres_ifc(:,:),         & !< Half level pressure field
      &  tsfc(:)                  !< Surface temperature
    INTEGER, INTENT(in)      :: &
      &  nlev, nlevp1,          & !< Number of vertical full and half levels
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
! Local variables
    INTEGER                  :: &
      &  jc, jk                   !< loop indices

      DO jk=1,nlevp1
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%pressure_hl(jc,jk)    = pres_ifc(jc,jk)
        ENDDO !jc
      ENDDO !jk
      ! Temperature at half levels is interpolated in the same way as in rrtm so far.
      DO jk=2,nlev
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%temperature_hl(jc,jk) =                                             &
            &                  (temp(jc,jk-1) * pres(jc,jk-1)  * ( pres(jc,jk) - pres_ifc(jc,jk) ) &
            &                + temp(jc,jk) * pres(jc,jk) * ( pres_ifc(jc,jk) - pres(jc,jk-1)))     &
            &                / ( pres_ifc(jc,jk) * (pres(jc,jk) - pres(jc,jk-1) ) )
        ENDDO !jc
      ENDDO !jk
      DO jc = i_startidx, i_endidx
        ecrad_thermodynamics%temperature_hl(jc,nlevp1)    = tsfc(jc)
        ecrad_thermodynamics%temperature_hl(jc,1)         = temp(jc,1)                        &
          &                   + ( pres_ifc(jc,1) - pres(jc,1) )                               &
          &                   * (temp(jc,1)      - ecrad_thermodynamics%temperature_hl(jc,2)) &
          &                   / (pres(jc,1)      - pres_ifc(jc,2) )
      ENDDO !jc

      ! Directly provide full level temperature and pressure to rrtm gas_optics in ecrad (see rrtm_pass_temppres_fl).
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%pressure_fl(jc,jk)    = pres(jc,jk)
          ecrad_thermodynamics%temperature_fl(jc,jk) = temp(jc,jk)
        ENDDO !jc
      ENDDO !jk
      
      CALL ecrad_thermodynamics%calc_saturation_wrt_liquid(istartcol=i_startidx, iendcol=i_endidx)

  END SUBROUTINE ecrad_set_thermodynamics
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_clouds:
  !! Set ecRad clouds information, i.e. fill t_ecrad_cloud_type with information
  !! from ICON.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-10)
  !!
  SUBROUTINE ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, qc, qi, clc, temp, pres, acdnc, fr_glac, fr_land, fact_reffc, &
    &                         clc_min, nlev, i_startidx, i_endidx)

    TYPE(t_ecrad_cloud_type), INTENT(inout) :: &
      &  ecrad_cloud              !< ecRad cloud information
    TYPE(t_ecrad_thermodynamics_type), INTENT(inout) :: &
      &  ecrad_thermodynamics     !< ecRad thermodynamics information
    REAL(wp), INTENT(in)     :: &
      &  qc(:,:),               & !< Total cloud water (gridscale + subgridscale)
      &  qi(:,:),               & !< Total cloud ice   (gridscale + subgridscale)
      &  clc(:,:),              & !< Cloud cover
      &  temp(:,:),             & !< Full level temperature field
      &  pres(:,:),             & !< Full level pressure field
      &  acdnc(:,:),            & !< Cloud droplet numb. conc. (m-3)
      &  fr_glac(:),            & !< fraction of land covered by glaciers
      &  fr_land(:),            & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  fact_reffc,            & !< Factor in the calculation of cloud droplet effective radius
      &  clc_min                  !< Minimum cloud cover value to be considered as partly cloudy
    INTEGER, INTENT(in)      :: &
      &  nlev,                  & !< Number of vertical full levels
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
! Local variables
    REAL(wp)                 :: &
      &  lwc, iwc,              & !< Cloud liquid and ice water content
      &  liwcfac                  !< Factor to calculate cloud liquid and ice water content
    INTEGER                  :: &
      &  jc, jk                   !< loop indices

    ! Currently hardcoded values decorrelation length need adaption
    CALL ecrad_cloud%set_overlap_param(ecrad_thermodynamics, 2000._wp, istartcol=i_startidx, iendcol=i_endidx)

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
        ecrad_cloud%q_liq(jc,jk)    = qc(jc,jk)
        ecrad_cloud%q_ice(jc,jk)    = qi(jc,jk)
        IF ( clc(jc,jk) > clc_min ) THEN
          liwcfac = 1000.0_wp / clc(jc,jk) * pres(jc,jk) / temp(jc,jk) / rd
          ecrad_cloud%fraction(jc,jk) = clc(jc,jk)
        ELSE
          liwcfac = 0._wp
          ecrad_cloud%fraction(jc,jk) = 0._wp
        ENDIF
        lwc                         = qc(jc,jk) * liwcfac
        iwc                         = qi(jc,jk) * liwcfac
        ! Careful with acdnc input: A division is performed and it is not checked for 0 as the function used
        ! to create acdnc returns always positive values
        ecrad_cloud%re_liq(jc,jk)   = reff_droplet(lwc, acdnc(jc,jk), fr_land(jc), fr_glac(jc), fact_reffc)
        ecrad_cloud%re_ice(jc,jk)   = reff_crystal(iwc)
      ENDDO
    ENDDO

  END SUBROUTINE ecrad_set_clouds
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_gas:
  !! Set ecRad gas information, i.e. fill t_ecrad_gas_type with information
  !! from ICON namelist parameters. The parameters irad_xyz and vmr_xyz are used
  !! (xyz=h2o, o3, co2, o2, cfc11, cfc12, n20 and ch4). Not all options are available for
  !! each of the gases. irad_xyz can have the following values:
  !!   0         : Set the gas globally constant to a concentration of zero
  !!   1         : Use information from a tracer variable
  !!   2         : Set the concentration to a globally constant value taken from vmr_xyz
  !!   3         : Use vmr_xyz at the surface, tanh-decay with height
  !!   7/9/79/97 : Use climatologies (only implemented for ozone)
  !! The finish calls in case default should never trigger as the values for irad_xyz
  !! were already checked in mo_nml_crosscheck.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-13)
  !!
  SUBROUTINE ecrad_set_gas(ecrad_gas, ecrad_conf, o3, qv, pres, i_startidx, i_endidx, nlev)

    TYPE(t_ecrad_gas_type), INTENT(inout) :: &
      &  ecrad_gas                !< ecRad cloud information
    TYPE(t_ecrad_conf), INTENT(in) :: &
      &  ecrad_conf               !< ecRad configuration type
    REAL(wp), INTENT(in)     :: &
      &  o3(:,:),               & !< Ozone concentration
      &  qv(:,:),               & !< Water vapor
      &  pres(:,:)                !< Full level pressure
    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlev                     !< Number of vertical full levels
    ! Local variables
    REAL(wp), ALLOCATABLE    :: &
      &  ch4(:,:),              & !< Methane volume mixing ratio
      &  n2o(:,:)                 !< N2O volume mixing ratio
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &  routine = modname//'::ecrad_set_gas'

    ! Water Vapor
    SELECT CASE(irad_h2o)
      CASE(0) ! No water vapor
        CALL ecrad_gas%put_well_mixed(ecRad_IH2O, IVolumeMixingRatio, 0._wp,  istartcol=i_startidx, iendcol=i_endidx)
      CASE(1) ! Use values from diagnosed water vapor content
        CALL ecrad_gas%put(ecRad_IH2O, IMassMixingRatio, qv(:,:))
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_h2o = 0, 1')
    END SELECT

    ! Ozone
    SELECT CASE(irad_o3)
      CASE(0) ! No Ozone
        CALL ecrad_gas%put_well_mixed(ecRad_IO3,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(7,9,79,97) ! Use values from GEMS/MACC (different profiles)
        CALL ecrad_gas%put(ecRad_IO3,  IMassMixingRatio, o3(:,:))
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_o3 = 0, 7, 9, 79, 97')
    END SELECT

    !CO2
    SELECT CASE(irad_co2)
      CASE(0) ! No CO2
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_co2
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IVolumeMixingRatio, vmr_co2,  istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_co2 = 0, 2')
    END SELECT

    !O2
    SELECT CASE(irad_o2)
      CASE(0) ! No O2
        CALL ecrad_gas%put_well_mixed(ecRad_IO2,IVolumeMixingRatio, 0._wp,   istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_o2
        ! O2 is hardcoded within rrtm gas solver of ecRad, option for ecRad_IO2 only in case of psrad gas solver
        ! We still put it in ecRad, because this bug should be fixed with the next ecRad release
        CALL ecrad_gas%put_well_mixed(ecRad_IO2,IVolumeMixingRatio, vmr_o2,  istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_o2 = 0, 2')
    END SELECT

    !CFC11
    SELECT CASE(irad_cfc11)
      CASE(0) ! No CFC11
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_cfc11
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IVolumeMixingRatio, vmr_cfc11,istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_cfc11 = 0, 2')
    END SELECT

    !CFC12
    SELECT CASE(irad_cfc12)
      CASE(0) ! No CFC12
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_cfc12
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IVolumeMixingRatio, vmr_cfc12,istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_cfc12 = 0, 2')
    END SELECT

    !N2O
    SELECT CASE(irad_n2o)
      CASE(0) ! No N2O
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_n2o
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IVolumeMixingRatio, vmr_n2o,istartcol=i_startidx,iendcol=i_endidx)
      CASE(3) ! Tanh profile
        ALLOCATE(n2o(i_startidx:i_endidx,nlev))
        n2o(:,:)=gas_profile(vmr_n2o, pres, vpp_n2o, i_startidx, i_endidx, nlev)
        CALL ecrad_gas%put(ecRad_IN2O,  IVolumeMixingRatio, n2o(:,:))
        DEALLOCATE(n2o)
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_n2o = 0, 2, 3')
    END SELECT

    !CH4
    SELECT CASE(irad_ch4)
      CASE(0) ! No CH4
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_ch4
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IVolumeMixingRatio, vmr_ch4,istartcol=i_startidx,iendcol=i_endidx)
      CASE(3) ! Tanh profile
        ALLOCATE(ch4(i_startidx:i_endidx,nlev))
        ch4(:,:)=gas_profile(vmr_ch4, pres, vpp_ch4, i_startidx, i_endidx, nlev)
        CALL ecrad_gas%put(ecRad_ICH4,  IVolumeMixingRatio, ch4(:,:))
        DEALLOCATE(ch4)
      CASE DEFAULT
        CALL finish(TRIM(routine),'Current implementation only supports irad_ch4 = 0, 2, 3')
    END SELECT

    CALL ecrad_set_gas_units(ecrad_conf, ecrad_gas)
    ! Possible further gases to be added: CO, HCFC22, CCl4, NO2

  END SUBROUTINE ecrad_set_gas
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_store_fluxes:
  !! Stores radiative fluxes calculated by ecRad from the ecRad type t_ecrad_flux_type 
  !! in the corresponding ICON data structure.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-13)
  !!
  SUBROUTINE ecrad_store_fluxes(ecrad_flux, cosmu0, trsolall, trsol_up_toa, trsol_up_sfc, trsol_par_sfc,  &
    &                           trsol_dn_sfc_diff, trsolclr_sfc, lwflxall, lwflx_up_sfc_rs, lwflxclr_sfc, &
    &                           cosmu0mask, i_startidx, i_endidx, nlevp1)

    TYPE(t_ecrad_flux_type), INTENT(inout) :: &
      &  ecrad_flux               !< ecRad cloud information
    REAL(wp), INTENT(inout)  :: &
      &  cosmu0(:),             & !< Cosine of solar zenith angle
      &  trsolall(:,:),         & !< solar transmissivity, all sky, net down
      &  trsol_up_toa(:),       & !< upward solar transmissivity at TOA
      &  trsol_up_sfc(:),       & !< upward solar transmissivity at surface
      &  trsol_par_sfc(:),      & !< downward transmissivity for photosynthetically active rad. at surface
      &  trsol_dn_sfc_diff(:),  & !< downward diffuse solar transmissivity at surface
      &  trsolclr_sfc(:),       & !< clear-sky net transmissivity at surface
      &  lwflxall(:,:),         & !< terrestrial flux, all sky, net down
      &  lwflx_up_sfc_rs(:),    & !< longwave upward flux at surface
      &  lwflxclr_sfc(:)          !< longwave clear-sky flux at surface
    LOGICAL, INTENT(in)      :: &
      &  cosmu0mask(:)            !< Mask if cosmu0 > 0
    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlevp1                   !< Number of vertical half levels
    ! Local Variables
    INTEGER                  :: &
      &  jband, jc, jk            !< Loop indices

      ! Initialize output fields
      trsolall(:,:)        = 0._wp
      trsol_up_toa(:)      = 0._wp
      trsol_up_sfc(:)      = 0._wp
      trsol_par_sfc(:)     = 0._wp
      trsol_dn_sfc_diff(:) = 0._wp
      trsolclr_sfc(:)      = 0._wp

      ! Store output of 3-D Fluxes
      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx
          IF ( cosmu0mask(jc) ) THEN 
            ! solar transmissivity, all sky, net down
            trsolall(jc,jk)     = (ecrad_flux%sw_dn(jc,jk)-ecrad_flux%sw_up(jc,jk))/cosmu0(jc)
          ENDIF
          ! terrestrial flux, all sky, net down
          lwflxall(jc,jk)       = ecrad_flux%lw_dn(jc,jk)-ecrad_flux%lw_up(jc,jk)
        ENDDO
      ENDDO

      ! Store output of 2-D Fluxes
      DO jc = i_startidx, i_endidx
        lwflx_up_sfc_rs(jc) = ecrad_flux%lw_up(jc,nlevp1)
        lwflxclr_sfc(jc)    = ecrad_flux%lw_dn_clear(jc,nlevp1) - ecrad_flux%lw_up_clear(jc,nlevp1)

        IF ( cosmu0mask(jc) ) THEN
          trsol_up_toa(jc)  = ecrad_flux%sw_up(jc,1)/cosmu0(jc)
          trsol_up_sfc(jc)  = ecrad_flux%sw_up(jc,nlevp1)/cosmu0(jc)
          DO jband = 1, nweight_par_ecrad
            trsol_par_sfc(jc)   = trsol_par_sfc(jc)                                                       &
              &  + (weight_par_ecrad(jband) * ecrad_flux%sw_dn_surf_band(iband_par_ecrad(jband),jc))      &
              &  / cosmu0(jc)
          ENDDO
          trsol_dn_sfc_diff(jc) = (ecrad_flux%sw_dn(jc,nlevp1) - ecrad_flux%sw_dn_direct(jc,nlevp1))      &
            &                     / cosmu0(jc)
          trsolclr_sfc(jc)      = (ecrad_flux%sw_dn_clear(jc,nlevp1) - ecrad_flux%sw_up_clear(jc,nlevp1)) &
            &                     / cosmu0(jc)
        ENDIF
      ENDDO

  END SUBROUTINE ecrad_store_fluxes
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function create_rdm_seed:
  !! Create a unique but reproducable random seed for the McICA solver
  !! 
  !! -------Algorithm taken from IFS, courtesy of R.J. Hogan-------
  !! This method gives a unique value for roughly every 1-km square
  !! on the globe and every minute.  (lat * rad2deg) gives rough
  !! latitude in degrees, which we multiply by 100 to give a unique
  !! value for roughly every km. lon*60*100 gives a unique number
  !! for roughly every km of longitude around the equator, which we
  !! multiply by 180*100 so there is no overlap with the latitude
  !! values.  The result can be contained in a 32-byte integer (but
  !! since random numbers are generated with the help of integer
  !! overflow, it should not matter if the number did overflow).
  !! 
  !! A more simple algorithm using the sum of int(lat), int(lon) and
  !! int(simtime) creates stripe patterns in the instantaneous fluxes.
  !! --------------------------------------------------------------
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !!
  FUNCTION create_rdm_seed(lon,lat,current_datetime)
    ! In:
    REAL(wp),       INTENT(in) :: &
      &  lon, lat             !< Longitude and Latitude value (radian)
    TYPE(datetime), INTENT(in) :: &
      &  current_datetime     !< Current date and time
    ! Out:
    INTEGER              :: &
      &  create_rdm_seed
    ! Local:
    INTEGER              :: &
      &  time, day            !< time of the day in minutes, the day of the month

    day  = INT(current_datetime%date%day)
    time = (INT(current_datetime%time%hour) * 60) + INT(current_datetime%time%minute)

    create_rdm_seed = time + day +  NINT(lon*108000000._wp + (lat * rad2deg * 100._wp) )

  END FUNCTION create_rdm_seed
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function reff_crystal:
  !! Function to calculate effective radius of ice crystals (extracted for the use in ecRad
  !! from mo_newcld_optics.f90). Author of the original code: Bjorn Stevens, MPI-M, Hamburg
  !! see ECHAM5 documentation (Roeckner et al, MPI report 349)
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !!
  REAL(wp) FUNCTION reff_crystal(ziwc)  ![m]
    REAL(wp), INTENT (IN)  :: &
      &  ziwc                    !< ice water content (g/m3)
    REAL(wp)                :: &
      &  reff_crystal_min,     & !< Minimum value ice crystal effective radius
      &  reff_crystal_max        !< Maximum value ice crystal effective radius
    
    ! Minimum and maximum value (derived from file ECHAM6_CldOptProps.nc)
    reff_crystal_min = 4._wp
    reff_crystal_max = 99._wp ! 124._wp < modified as values > 100 mu m lead to crashes in ecRad
                              ! (if delta_eddington_scat_od = .false.)

    reff_crystal = MAX(reff_crystal_min ,MIN(reff_crystal_max  ,83.8_wp*ziwc**0.216_wp)) * 1.e-6_wp
  END FUNCTION
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function reff_droplet:
  !! Function to calculate effective radius of water droplets (extracted for the use in ecRad
  !! from mo_newcld_optics.f90). Author of the original code: Bjorn Stevens, MPI-M, Hamburg
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !!
  REAL(wp) FUNCTION reff_droplet(zlwc,zcdnc,zland,zglac,zfact)  ![m]
    REAL (wp), INTENT (IN)  :: &
      &  zlwc,                 & !< liquid water content (g/m3)
      &  zcdnc,                & !< cloud drop number concentration
      &  zglac,                & !< fraction of land covered by glaciers
      &  zland,                & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zfact                   !< factor
    REAL(wp)                :: &
      &  zkap,                 & !< Factor
      &  reff_droplet_min,     & !< Minimum value ice crystal effective radius
      &  reff_droplet_max        !< Maximum value ice crystal effective radius
    REAL (wp), PARAMETER ::    &
      &  zkap_cont = 1.143_wp, & !< continental (Martin et al. ) breadth param
      &  zkap_mrtm = 1.077_wp    !< maritime (Martin et al.) breadth parameter

    ! Minimum and maximum value (derived from file ECHAM6_CldOptProps.nc)
    reff_droplet_min = 2.e-6_wp
    reff_droplet_max = 32.e-6_wp

    zkap         = zkap_cont*(zland-zglac) + zkap_mrtm*(1.0_wp-zland+zglac)
    reff_droplet = MAX(reff_droplet_min,MIN(reff_droplet_max,zfact*zkap*(zlwc / (zcdnc * 1.e-6_wp))**(1.0_wp/3.0_wp)))
  END FUNCTION
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function gas_profile:
  !! Function to calculate gas profile decaying with height with a tanh function.
  !! Extracted for the use in ecRad from mo_radiation.f90.
  !! Author of the original code: Bjorn Stevens, MPI-M, Hamburg
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-14)
  !!
  FUNCTION gas_profile(vmr_gas, pres, xp, i_startidx, i_endidx, nlev)

    REAL(wp), INTENT (in)   :: &
      &  vmr_gas,              & !< Constant volume mixing ratio of gas specified via namelist
      &  pres(:,:),            & !< Full level pressure
      &  xp(3)                   !< Gas-specific coefficient
    INTEGER,  INTENT (in)   :: &
      &  i_startidx, i_endidx, & !< Start and end index of nproma loop in current block
      &  nlev                    !< Number of vertical full levels
    REAL(wp)                :: &
      &  zx_d, zx_m,           &
      &  gas_profile(i_startidx:i_endidx,nlev) !< Profile to be calculated

    zx_m = (vmr_gas+xp(1)*vmr_gas)*0.5_wp
    zx_d = (vmr_gas-xp(1)*vmr_gas)*0.5_wp

    gas_profile(i_startidx:i_endidx,:) = (1._wp-(zx_d/zx_m)*TANH(LOG(pres(i_startidx:i_endidx,:)   &
          &                             /xp(2)) /xp(3))) * zx_m
  END FUNCTION
#endif
END MODULE mo_nwp_ecrad_utilities
