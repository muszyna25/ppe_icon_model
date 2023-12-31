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
  USE mo_math_constants,         ONLY: rad2deg, pi
  USE mo_exception,              ONLY: finish
  USE mo_math_types,             ONLY: t_geographical_coordinates
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_physical_constants,     ONLY: rd, grav
  USE mo_radiation_config,       ONLY: vmr_co2, vmr_n2o, vmr_o2, vmr_ch4,        &
                                   &   vmr_cfc11, vmr_cfc12,                     &
                                   &   irad_h2o, irad_o3, irad_co2,              &
                                   &   irad_n2o, irad_ch4,                       &
                                   &   irad_o2, irad_cfc11, irad_cfc12,          &
                                   &   vpp_ch4, vpp_n2o, tsi_radt
  USE mo_nwp_tuning_config,      ONLY: tune_difrad_3dcont
  USE mtime,                     ONLY: datetime
  USE mo_bc_greenhouse_gases,    ONLY: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr
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

  USE mo_exception,              ONLY: message
  USE mo_grid_config,            ONLY: l_scm_mode
  USE mo_scm_nml,                ONLY: lon_scm, lat_scm 

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
  PUBLIC :: add_3D_diffuse_rad

  ! helper functions to be removed once acc is merged into libecrad
  PUBLIC :: ecrad_acc_allocation
  PUBLIC :: ecrad_acc_deallocation
  PUBLIC :: update_host_pre_ecrad
  PUBLIC :: update_device_post_ecrad

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
    &                               albvisdif, albnirdif, albvisdir, albnirdir, emis_rad, i_startidx, i_endidx, use_acc)

    TYPE(t_ecrad_single_level_type), INTENT(inout) :: &
      &  ecrad_single_level       !< ecRad single level information
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  current_datetime         !< Current date and time
    TYPE(t_geographical_coordinates), INTENT(in), TARGET :: &
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
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
! Local variables
    INTEGER                  :: &
      &  jc,                    & !< loop index
      &  seed_in_time             !< temporal contribution to the seed

    TYPE(t_geographical_coordinates), TARGET, ALLOCATABLE :: scm_center(:)
    TYPE(t_geographical_coordinates), POINTER             :: ptr_center(:)

    LOGICAL                  :: lacc


    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if

    ! SCM: read lat/lon for horizontally uniform zenith angle
    IF ( l_scm_mode ) THEN
#ifdef _OPENACC
      IF (lacc) CALL finish('ecrad_set_single_level','l_scm_mode not ported to gpu')
#endif
      ALLOCATE(scm_center(SIZE(cell_center)))
      DO jc = i_startidx, i_endidx
        scm_center(jc)%lat = lat_scm * pi/180.
        scm_center(jc)%lon = lon_scm * pi/180.
      ENDDO
      ptr_center => scm_center
    ELSE
      ptr_center => cell_center
    ENDIF

    seed_in_time = create_rdm_seed_in_time(current_datetime)

    !$ACC PARALLEL DEFAULT(NONE) PRESENT(cosmu0, tsfc, albvisdif, albnirdif, albvisdir,  &
    !$ACC                               albnirdir, emis_rad, ecrad_single_level) &
    !$ACC                        COPYIN(ptr_center) IF(lacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
        ecrad_single_level%cos_sza(jc)            = cosmu0(jc)
        ecrad_single_level%skin_temperature(jc)   = tsfc(jc)
        ecrad_single_level%sw_albedo(jc,1)        = albvisdif(jc)
        ecrad_single_level%sw_albedo(jc,2)        = albnirdif(jc)
        ecrad_single_level%sw_albedo_direct(jc,1) = albvisdir(jc)
        ecrad_single_level%sw_albedo_direct(jc,2) = albnirdir(jc)
        ecrad_single_level%lw_emissivity(jc,1)    = emis_rad(jc)
        ecrad_single_level%iseed(jc)              = create_rdm_seed(ptr_center(jc)%lon, &
          &                                                         ptr_center(jc)%lat, &
          &                                                         seed_in_time)
    ENDDO !jc
    !$ACC END PARALLEL

    IF ( l_scm_mode ) THEN
      DEALLOCATE(scm_center)
    ENDIF

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
    &                                 nlev, nlevp1, i_startidx, i_endidx, use_acc)

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
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
! Local variables
    INTEGER                  :: &
      &  jc, jk                   !< loop indices
    LOGICAL                  :: lacc

      if (present(use_acc)) then
        lacc = use_acc
      else
        lacc = .false.
      end if

      !$ACC DATA PRESENT(ecrad_thermodynamics, temp, pres, pres_ifc) IF(lacc)

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(NONE) ASYNC(1) IF(lacc)
      DO jk=1,nlevp1
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%pressure_hl(jc,jk)    = pres_ifc(jc,jk)
        ENDDO !jc
      ENDDO !jk
      !$ACC END PARALLEL

      ! Temperature at half levels is interpolated in the same way as in rrtm so far.
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP SEQ
      DO jk=2,nlev
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%temperature_hl(jc,jk) =                                             &
            &                  (temp(jc,jk-1) * pres(jc,jk-1)  * ( pres(jc,jk) - pres_ifc(jc,jk) ) &
            &                + temp(jc,jk) * pres(jc,jk) * ( pres_ifc(jc,jk) - pres(jc,jk-1)))     &
            &                / ( pres_ifc(jc,jk) * (pres(jc,jk) - pres(jc,jk-1) ) )
        ENDDO !jc
      ENDDO !jk
      !$ACC END PARALLEL

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) ASYNC(1) IF(lacc)
      DO jc = i_startidx, i_endidx
        ecrad_thermodynamics%temperature_hl(jc,nlevp1) = temp(jc,nlev) + (pres_ifc(jc,nlevp1) - pres(jc,nlev)) * &
                               (temp(jc,nlev-1) - temp(jc,nlev))/(pres(jc,nlev-1) - pres(jc,nlev))
        ecrad_thermodynamics%temperature_hl(jc,1)         = temp(jc,1)                        &
          &                   + ( pres_ifc(jc,1) - pres(jc,1) )                               &
          &                   * (temp(jc,1)      - ecrad_thermodynamics%temperature_hl(jc,2)) &
          &                   / (pres(jc,1)      - pres_ifc(jc,2) )
      ENDDO !jc
      !$ACC END PARALLEL

      ! Directly provide full level temperature and pressure to rrtm gas_optics in ecrad (see rrtm_pass_temppres_fl).
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(NONE) ASYNC(1) IF(lacc)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%pressure_fl(jc,jk)    = pres(jc,jk)
          ecrad_thermodynamics%temperature_fl(jc,jk) = temp(jc,jk)
        ENDDO !jc
      ENDDO !jk
      !$ACC END PARALLEL

#ifndef __ECRAD_ACC
      !$ACC UPDATE HOST( ecrad_thermodynamics%pressure_hl, ecrad_thermodynamics%temperature_hl ) IF( lacc )
#endif
      CALL ecrad_thermodynamics%calc_saturation_wrt_liquid(istartcol=i_startidx, iendcol=i_endidx, use_acc=lacc)
#ifndef __ECRAD_ACC
      !$ACC ENTER DATA CREATE(ecrad_thermodynamics%h2o_sat_liq) IF(lacc)
      !$ACC UPDATE DEVICE(ecrad_thermodynamics%h2o_sat_liq) IF(lacc)
#endif

      !$ACC END DATA

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
  SUBROUTINE ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, qc, qi, clc, temp, pres, acdnc, fr_glac, fr_land, &
    &                         qr,qs,qg,reff_liq, reff_frz, reff_rain, reff_snow, reff_graupel,                     & 
    &                         icpl_reff, fact_reffc, clc_min, nlev, i_startidx, i_endidx, use_acc)

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
      &  fact_reffc,            & !< Factor in the calculation of cloud droplet effective radius
      &  clc_min                  !< Minimum cloud cover value to be considered as partly cloudy
    REAL(wp), POINTER, INTENT(in)     :: &
      &  acdnc(:,:),            & !< Cloud droplet numb. conc. (m-3)
      &  fr_glac(:),            & !< fraction of land covered by glaciers
      &  fr_land(:),            & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  qr(:,:),               & !< rain
      &  qs(:,:),               & !< snow
      &  qg(:,:),               & !< graupel
      &  reff_liq(:,:),         & !< effective radius of the liquid phase (external)
      &  reff_frz(:,:),         & !< effective radius of the frozen phase (external)
      &  reff_rain(:,:),        & !< effective radius of the rain phase (external)
      &  reff_snow(:,:),        & !< effective radius of the snow phase (external)
      &  reff_graupel(:,:)        !< effective radius of the graupel phase (external)

    INTEGER, INTENT(in)      :: &
      &  icpl_reff,             & !< Option for effective radius
      &  nlev,                  & !< Number of vertical full levels
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
    
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
! Local variables
    REAL(wp)                 :: &
      &  lwc, iwc,              & !< Cloud liquid and ice water content
      &  liwcfac                  !< Factor to calculate cloud liquid and ice water content
    INTEGER                  :: &
      &  jc, jk                   !< loop indices
    LOGICAL                  :: lacc

    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if

    ! Currently hardcoded values decorrelation length need adaption
#ifndef __ECRAD_ACC
    !$ACC UPDATE HOST(ecrad_thermodynamics%pressure_hl, ecrad_thermodynamics%temperature_hl) IF(lacc)
#endif
    CALL ecrad_cloud%set_overlap_param(ecrad_thermodynamics, 2000._wp, istartcol=i_startidx, iendcol=i_endidx, &
      &                                use_acc=lacc)
#ifndef __ECRAD_ACC
    !$ACC UPDATE DEVICE(ecrad_cloud%overlap_param) IF(lacc)
#endif

    !$ACC DATA PRESENT(ecrad_cloud, qc, qi, clc, temp, pres, &
    !$ACC              acdnc, fr_land, fr_glac, reff_frz, reff_liq) IF(lacc)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(liwcfac, lwc, iwc)
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
        IF ( icpl_reff == 0 ) THEN ! No external calculationcof reff.
          lwc                         = qc(jc,jk) * liwcfac
          iwc                         = qi(jc,jk) * liwcfac
          ! Careful with acdnc input: A division is performed and it is not checked for 0 as the function used
          ! to create acdnc returns always positive values
          ecrad_cloud%re_liq(jc,jk)   = reff_droplet(lwc, acdnc(jc,jk), fr_land(jc), fr_glac(jc), fact_reffc)
          ecrad_cloud%re_ice(jc,jk)   = reff_crystal(iwc)
        END IF
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    IF ( icpl_reff > 0 ) THEN
      IF (.NOT. ASSOCIATED(reff_liq) .OR. .NOT. ASSOCIATED(reff_frz)) THEN
        CALL finish('ecrad_set_clouds','effective radius fields not associated')
      ENDIF
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ecrad_cloud%re_liq(jc,jk) = MAX(MIN(reff_liq(jc,jk),32.0e-6_wp),2.0e-6_wp)  
          ecrad_cloud%re_ice(jc,jk) = MAX(MIN(reff_frz(jc,jk),99.0e-6_wp),5.0e-6_wp) 
        ENDDO
      ENDDO
    !$ACC END PARALLEL
    ENDIF

    !$ACC END DATA
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
  !!   11        : Read ozone from SCM input instead of here.
  !! The finish calls in case default should never trigger as the values for irad_xyz
  !! were already checked in mo_nml_crosscheck.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-13)
  !!
  SUBROUTINE ecrad_set_gas(ecrad_gas, ecrad_conf, o3, qv, pres, i_startidx, i_endidx, nlev, use_acc)

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
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
    ! Local variables
    REAL(wp), ALLOCATABLE    :: &
      &  ch4(:,:),              & !< Methane volume mixing ratio
      &  n2o(:,:)                 !< N2O volume mixing ratio
    LOGICAL                  :: lacc
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//'::ecrad_set_gas'

    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if

#ifndef __ECRAD_ACC
    !$ACC UPDATE HOST(o3, qv, pres) IF(lacc)
    lacc = .false.
#endif

    ! Water Vapor
    SELECT CASE(irad_h2o)
      CASE(0) ! No water vapor
        CALL ecrad_gas%put_well_mixed(ecRad_IH2O, IVolumeMixingRatio, 0._wp,  istartcol=i_startidx, iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(1) ! Use values from diagnosed water vapor content
        CALL ecrad_gas%put(ecRad_IH2O, IMassMixingRatio, qv(:,:), use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_h2o = 0, 1')
    END SELECT

    ! Ozone
    SELECT CASE(irad_o3)
      CASE(0) ! No Ozone
        CALL ecrad_gas%put_well_mixed(ecRad_IO3,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(7,9,79,97) ! Use values from GEMS/MACC (different profiles)
        CALL ecrad_gas%put(ecRad_IO3,  IMassMixingRatio, o3(:,:), use_acc=lacc)
      CASE(11) ! Ozone is read from SCM input file
        CALL message('mo_nwp_ecrad_utilities: irad_o3=11', &
          &          'Ozone used for radiation is read from SCM input file')
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_o3 = 0, 7, 9, 79, 97, 11')
    END SELECT

    !CO2
    SELECT CASE(irad_co2)
      CASE(0) ! No CO2
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(2) ! Constant value derived from namelist parameter vmr_co2
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IVolumeMixingRatio, vmr_co2,  istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IMassMixingRatio, ghg_co2mmr,  istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_co2 = 0, 2, 4')
    END SELECT

    !O2
    SELECT CASE(irad_o2)
      CASE(0) ! No O2
        CALL ecrad_gas%put_well_mixed(ecRad_IO2,IVolumeMixingRatio, 0._wp,   istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(2) ! Constant value derived from namelist parameter vmr_o2
        ! O2 is hardcoded within rrtm gas solver of ecRad, option for ecRad_IO2 only in case of psrad gas solver
        ! We still put it in ecRad, because this bug should be fixed with the next ecRad release
        CALL ecrad_gas%put_well_mixed(ecRad_IO2,IVolumeMixingRatio, vmr_o2,  istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_o2 = 0, 2')
    END SELECT

    !CFC11
    SELECT CASE(irad_cfc11)
      CASE(0) ! No CFC11
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx,&
          &                           use_acc=lacc)
      CASE(2) ! Constant value derived from namelist parameter vmr_cfc11
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IVolumeMixingRatio, vmr_cfc11,istartcol=i_startidx,iendcol=i_endidx,&
          &                           use_acc=lacc)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IMassMixingRatio, ghg_cfcmmr(1),istartcol=i_startidx, &
          &                           iendcol=i_endidx, use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_cfc11 = 0, 2, 4')
    END SELECT

    !CFC12
    SELECT CASE(irad_cfc12)
      CASE(0) ! No CFC12
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx,&
          &                           use_acc=lacc)
      CASE(2) ! Constant value derived from namelist parameter vmr_cfc12
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IVolumeMixingRatio, vmr_cfc12,istartcol=i_startidx,iendcol=i_endidx,&
          &                           use_acc=lacc)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IMassMixingRatio, ghg_cfcmmr(2),istartcol=i_startidx, &
          &                           iendcol=i_endidx, use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_cfc12 = 0, 2, 4')
    END SELECT

    !N2O
    SELECT CASE(irad_n2o)
      CASE(0) ! No N2O
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(2) ! Constant value derived fromecrad_set_gas namelist parameter vmr_n2o
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IVolumeMixingRatio, vmr_n2o,istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(3) ! Tanh profile
        ALLOCATE(n2o(i_startidx:i_endidx,nlev))
        !$ACC ENTER DATA CREATE(n2o) ASYNC(1) IF(lacc)
        CALL gas_profile(vmr_n2o, pres, vpp_n2o, i_startidx, i_endidx, nlev, n2o(:,:), use_acc=lacc)
        CALL ecrad_gas%put(ecRad_IN2O,  IVolumeMixingRatio, n2o(:,:), use_acc=lacc)
        !$ACC WAIT
        !$ACC EXIT DATA DELETE(n2o) IF(lacc)
        DEALLOCATE(n2o)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IMassMixingRatio, ghg_n2ommr,istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_n2o = 0, 2, 3, 4')
    END SELECT

    !CH4
    SELECT CASE(irad_ch4)
      CASE(0) ! No CH4
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(2) ! Constant value derived from namelist parameter vmr_ch4
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IVolumeMixingRatio, vmr_ch4,istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE(3) ! Tanh profile
        ALLOCATE(ch4(i_startidx:i_endidx,nlev))
        !$ACC ENTER DATA CREATE(ch4) ASYNC(1) IF(lacc)
        CALL gas_profile(vmr_ch4, pres, vpp_ch4, i_startidx, i_endidx, nlev, ch4(:,:), use_acc=lacc)
        CALL ecrad_gas%put(ecRad_ICH4,  IVolumeMixingRatio, ch4(:,:), use_acc=lacc)
        !$ACC WAIT
        !$ACC EXIT DATA DELETE(ch4) IF(lacc)
        DEALLOCATE(ch4)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IMassMixingRatio, ghg_ch4mmr,istartcol=i_startidx,iendcol=i_endidx, &
          &                           use_acc=lacc)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_ch4 = 0, 2, 3, 4')
    END SELECT

    CALL ecrad_set_gas_units(ecrad_conf, ecrad_gas, use_acc=lacc)
    ! Possible further gases to be added: CO, HCFC22, CCl4, NO2
    !$ACC UPDATE DEVICE(ecrad_gas%mixing_ratio) IF(lacc)

#ifndef __ECRAD_ACC
    !$ACC UPDATE DEVICE(ecrad_gas%mixing_ratio) IF(lacc)
#endif

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
  SUBROUTINE ecrad_store_fluxes(jg, ecrad_flux, cosmu0, trsolall, trsol_up_toa, trsol_up_sfc, trsol_par_sfc,  &
    &                           trsol_dn_sfc_diff, trsolclr_sfc, lwflxall, lwflx_up_sfc_rs, lwflxclr_sfc, &
    &                           lwflx_up    , lwflx_dn    , swflx_up    , swflx_dn,                       &
    &                           lwflx_up_clr, lwflx_dn_clr, swflx_up_clr, swflx_dn_clr,                   &
    &                           cosmu0mask, i_startidx, i_endidx, nlevp1, use_acc)

    INTEGER, INTENT(in)   :: &
      &  jg                       !< domain index
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
      &  lwflxclr_sfc(:),       & !< longwave clear-sky flux at surface
      &  lwflx_up(:,:),         & !< longwave  3D upward   flux            [W/m2]
      &  lwflx_dn(:,:),         & !< longwave  3D downward flux            [W/m2]
      &  swflx_up(:,:),         & !< shortwave 3D upward   flux            [W/m2]
      &  swflx_dn(:,:),         & !< shortwave 3D downward flux            [W/m2]
      &  lwflx_up_clr(:,:),     & !< longwave  3D upward   flux clear-sky  [W/m2]
      &  lwflx_dn_clr(:,:),     & !< longwave  3D downward flux clear-sky  [W/m2]
      &  swflx_up_clr(:,:),     & !< shortwave 3D upward   flux clear-sky  [W/m2]
      &  swflx_dn_clr(:,:)        !< shortwave 3D downward flux clear-sky  [W/m2]

    LOGICAL, INTENT(in)      :: &
      &  cosmu0mask(:)            !< Mask if cosmu0 > 0
    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlevp1                   !< Number of vertical half levels
    LOGICAL, INTENT(in), OPTIONAL :: use_acc

    ! Local Variables
    INTEGER                  :: &
      &  jband, jc, jk            !< Loop indices

    LOGICAL                  :: lacc

    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if
    
      !$ACC DATA PRESENT(ecrad_flux , cosmu0, trsolall, trsol_up_toa, trsol_up_sfc, trsol_par_sfc,      &
      !$ACC             trsol_dn_sfc_diff, trsolclr_sfc, lwflxall, lwflx_up_sfc_rs, lwflxclr_sfc,       &
      !$ACC             lwflx_up, lwflx_dn, swflx_up, swflx_dn, lwflx_up_clr, lwflx_dn_clr,             &
      !$ACC             swflx_up_clr, swflx_dn_clr, cosmu0mask) IF(lacc)

      ! Initialize output fields
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP SEQ
      DO jk = 1, nlevp1
        !$ACC LOOP GANG(STATIC:1) VECTOR
        DO jc = 1,SIZE(trsolall,1)
          trsolall(jc,jk)      = 0._wp
        ENDDO
      ENDDO
      !$ACC LOOP GANG(STATIC:1) VECTOR
      DO jc = 1,SIZE(trsolall,1)
        trsol_up_toa(jc)      = 0._wp
        trsol_up_sfc(jc)      = 0._wp
        trsol_par_sfc(jc)     = 0._wp
        trsol_dn_sfc_diff(jc) = 0._wp
        trsolclr_sfc(jc)      = 0._wp
      ENDDO
      !$ACC END PARALLEL

      ! Store output of 3-D Fluxes
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
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
      !$ACC END PARALLEL

      IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN    
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlevp1
          DO jc = i_startidx, i_endidx
            ! LW/SW, up/down, all/clear 3D fluxes
            lwflx_up    (jc,jk)   = ecrad_flux%lw_up(jc,jk)
            lwflx_dn    (jc,jk)   = ecrad_flux%lw_dn(jc,jk)
  
            swflx_up    (jc,jk)   = ecrad_flux%sw_up(jc,jk)       * tsi_radt
            swflx_dn    (jc,jk)   = ecrad_flux%sw_dn(jc,jk)       * tsi_radt
            lwflx_up_clr(jc,jk)   = ecrad_flux%lw_up_clear(jc,jk)
            lwflx_dn_clr(jc,jk)   = ecrad_flux%lw_dn_clear(jc,jk)
            swflx_up_clr(jc,jk)   = ecrad_flux%sw_up_clear(jc,jk) * tsi_radt
            swflx_dn_clr(jc,jk)   = ecrad_flux%sw_dn_clear(jc,jk) * tsi_radt   
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      END IF

      ! Store output of 2-D Fluxes
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        lwflx_up_sfc_rs(jc) = ecrad_flux%lw_up(jc,nlevp1)
        lwflxclr_sfc(jc)    = ecrad_flux%lw_dn_clear(jc,nlevp1) - ecrad_flux%lw_up_clear(jc,nlevp1)

        IF ( cosmu0mask(jc) ) THEN
          trsol_up_toa(jc)  = ecrad_flux%sw_up(jc,1)/cosmu0(jc)
          trsol_up_sfc(jc)  = ecrad_flux%sw_up(jc,nlevp1)/cosmu0(jc)

          trsol_dn_sfc_diff(jc) = (ecrad_flux%sw_dn(jc,nlevp1) - ecrad_flux%sw_dn_direct(jc,nlevp1))      &
            &                     / cosmu0(jc)
          trsolclr_sfc(jc)      = (ecrad_flux%sw_dn_clear(jc,nlevp1) - ecrad_flux%sw_up_clear(jc,nlevp1)) &
            &                     / cosmu0(jc)
        ENDIF
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP SEQ
      DO jband = 1, nweight_par_ecrad
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          IF ( cosmu0mask(jc) ) THEN
            trsol_par_sfc(jc)   = trsol_par_sfc(jc)                                                       &
              &  + (weight_par_ecrad(jband) * ecrad_flux%sw_dn_surf_band(iband_par_ecrad(jband),jc))      &
              &  / cosmu0(jc)
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC END DATA

  END SUBROUTINE ecrad_store_fluxes
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE add_3D_diffuse_rad:
  !! Adds 3D contribution to diffuse radiation by reflection of direct solar radiation on scattered low clouds
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, Deutscher Wetterdienst, Offenbach (2019-12-06)
  !!
  SUBROUTINE add_3D_diffuse_rad(ecrad_flux, clc, pres, temp, cosmu0, trsol_dn_sfc_diff, i_startidx, i_endidx, nlev, &
      &                         use_acc)

    TYPE(t_ecrad_flux_type), INTENT(inout) :: ecrad_flux !< ecRad cloud information

    REAL(wp), INTENT(in)  :: &
      &  cosmu0(:),             & !< Cosine of solar zenith angle
      &  clc(:,:),              & !< cloud cover fraction
      &  pres(:,:),             & !< pressure
      &  temp(:,:)                !< temperature

    REAL(wp), INTENT(inout)  :: trsol_dn_sfc_diff(:) !< downward diffuse solar transmissivity at surface

    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlev                     !< Number of vertical levels
    LOGICAL, INTENT(in), OPTIONAL :: use_acc

    ! Local Variables
    INTEGER                   ::  jc, jk  !< Loop indices

    REAL(wp), PARAMETER :: zdecorr = 2000.0_wp, & ! decorrelation length scale for cloud overlap scheme
                           epsi    = 1.e-20_wp

    REAL(wp) :: zcloud(i_endidx), ccmax, ccran, deltaz, alpha
    LOGICAL  :: lacc

    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if

    !$ACC DATA PRESENT(ecrad_flux, pres, clc, temp, cosmu0, trsol_dn_sfc_diff) CREATE(zcloud) IF(lacc)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      zcloud(jc)     = 0.0_wp
    ENDDO
    !$ACC END PARALLEL

    ! Calculate low-level cloud cover fraction
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
    !$ACC LOOP SEQ
    DO jk = 2, nlev
      !$ACC LOOP GANG VECTOR PRIVATE(ccmax, ccran, deltaz, alpha)
      DO jc = i_startidx, i_endidx
        IF (pres(jc,jk)/pres(jc,nlev) > 0.75_wp) THEN
          ccmax = MAX(clc(jc,jk),  zcloud(jc))
          ccran = clc(jc,jk) + zcloud(jc) - clc(jc,jk)*zcloud(jc)

          ! layer thickness [m] between level jk and next upper level jk-1
          deltaz = (pres(jc,jk)-pres(jc,jk-1))/(pres(jc,jk-1)+pres(jc,jk)) * &
                   (temp(jc,jk-1)+temp(jc,jk))*rd/grav

          alpha  = MIN(EXP(-deltaz/zdecorr), clc(jc,jk-1)/MAX(epsi,clc(jc,jk)) )

          zcloud(jc) = alpha * ccmax + (1-alpha) * ccran
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      IF (cosmu0(jc) > 0.05_wp) THEN
        trsol_dn_sfc_diff(jc) = MIN(ecrad_flux%sw_dn(jc,nlev+1)/cosmu0(jc), trsol_dn_sfc_diff(jc) + &
          tune_difrad_3dcont*ecrad_flux%sw_dn(jc,nlev+1)/cosmu0(jc)*zcloud(jc)*(1._wp-zcloud(jc))**2)
      ENDIF
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE add_3D_diffuse_rad
  !---------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------
  !>
  !! Function create_rdm_seed_in_time:
  !! Create a unique but reproducable random seed in time for the McICA solver
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !!
  FUNCTION create_rdm_seed_in_time(current_datetime)
    !$ACC ROUTINE SEQ
    ! In:
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  current_datetime     !< Current date and time
    ! Out:
    INTEGER              :: &
      &  create_rdm_seed_in_time
    ! Local:
    INTEGER              :: &
      &  time, day            !< time of the day in minutes, the day of the month

    day  = INT(current_datetime%date%day)
    time = (INT(current_datetime%time%hour) * 60) + INT(current_datetime%time%minute)

    create_rdm_seed_in_time = time + day

  END FUNCTION create_rdm_seed_in_time
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
  FUNCTION create_rdm_seed(lon,lat,seed_in_time)
    !$ACC ROUTINE SEQ
    ! In:
    REAL(wp),       INTENT(in) :: &
      &  lon, lat             !< Longitude and Latitude value (radian)
    INTEGER, INTENT(in) :: &
      &  seed_in_time     !< contribution to the seed in time
    ! Out:
    INTEGER              :: &
      &  create_rdm_seed

    create_rdm_seed = seed_in_time +  NINT(lon*108000000._wp + (lat * rad2deg * 100._wp) )

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
    !$ACC ROUTINE SEQ
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
    !$ACC ROUTINE SEQ
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
  SUBROUTINE gas_profile(vmr_gas, pres, xp, i_startidx, i_endidx, nlev, profile, use_acc)

    REAL(wp), INTENT (IN)         :: &
      &  vmr_gas,                    & !< Constant volume mixing ratio of gas specified via namelist
      &  pres(:,:),                  & !< Full level pressure
      &  xp(3)                         !< Gas-specific coefficient
    INTEGER,  INTENT (IN)         :: &
      &  i_startidx, i_endidx,       & !< Start and end index of nproma loop in current block
      &  nlev                          !< Number of vertical full levels
    REAL(wp), INTENT (OUT)        :: &
      &  profile(i_startidx:i_endidx,nlev) !< Profile to be calculated
    LOGICAL, INTENT(IN), OPTIONAL :: &
      &  use_acc
    REAL(wp)                      :: &
      &  zx_d, zx_m
    INTEGER                       :: &
      &  jc, jk                        !< loop indices
    LOGICAL                       :: &
      &  lacc

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    zx_m = (vmr_gas+xp(1)*vmr_gas)*0.5_wp
    zx_d = (vmr_gas-xp(1)*vmr_gas)*0.5_wp

    !$ACC PARALLEL DEFAULT(NONE) PRESENT(profile, pres, xp) ASYNC(1) IF(lacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk=1,nlev
      DO jc = i_startidx, i_endidx
        profile(jc,jk) = (1._wp-(zx_d/zx_m)*TANH(LOG(pres(jc,jk) /xp(2)) /xp(3))) * zx_m
      ENDDO !jc
    ENDDO !jk
    !$ACC END PARALLEL

  END SUBROUTINE

#ifndef __ECRAD_ACC
  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_acc_allocation:
  !! The enter data create of the derived ecrad types that happen in the
  !! subruntines:
  !!   - ecrad_single_level%allocate
  !!   - ecrad_thermodynamics%allocate
  !!   - ecrad_gas%allocate
  !!   - ecrad_cloud%allocate
  !!   - ecrad_cloud%create_fractional_std
  !!   - ecrad_aerosol%allocate_direct
  !!   - ecrad_flux%allocate
  !! are extracted here. Also the arrays that are intitialized there are updated on the GPU.
  !! @par Revision History
  !! Initial release by Daniel Hupp, MeteoSwiss, Winterthur (2022-01-11)
  !!
  SUBROUTINE ecrad_acc_allocation(ecrad_conf, ecrad_single_level, &
      &  ecrad_thermodynamics, ecrad_gas, ecrad_cloud, ecrad_aerosol, &
      &  ecrad_flux, use_acc)

    USE mo_ecrad,                  ONLY: t_ecrad_aerosol_type

    TYPE(t_ecrad_conf),               INTENT(IN)   :: ecrad_conf
    TYPE(t_ecrad_single_level_type),  INTENT(IN)   :: ecrad_single_level
    TYPE(t_ecrad_thermodynamics_type),INTENT(IN)   :: ecrad_thermodynamics
    TYPE(t_ecrad_gas_type),           INTENT(IN)   :: ecrad_gas
    TYPE(t_ecrad_cloud_type),         INTENT(IN)   :: ecrad_cloud
    TYPE(t_ecrad_aerosol_type),       INTENT(IN)   :: ecrad_aerosol
    TYPE(t_ecrad_flux_type),          INTENT(IN)   :: ecrad_flux

    LOGICAL, INTENT(IN), OPTIONAL :: use_acc

    LOGICAL                  :: lacc

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    ! CALL ecrad_single_level%allocate(nproma_rad, 2, 1, .true., use_acc=lacc) !< use_sw_albedo_direct, 2 bands
    !$ACC ENTER DATA CREATE(ecrad_single_level%cos_sza, ecrad_single_level%lw_emissivity, &
    !$ACC   ecrad_single_level%sw_albedo, ecrad_single_level%sw_albedo_direct, ecrad_single_level%iseed) IF(lacc)
    !$ACC ENTER DATA CREATE(ecrad_single_level%skin_temperature) &
    !$ACC   IF(lacc .AND.        ecrad_single_level%is_simple_surface)
    !$ACC ENTER DATA CREATE(ecrad_single_level%lw_emission) &
    !$ACC   IF(lacc .AND. (.NOT. ecrad_single_level%is_simple_surface))

    ! CALL ecrad_thermodynamics%allocate(nproma_rad, nlev_rg, use_h2o_sat=.false., rrtm_pass_temppres_fl=.true., &
    !                                    use_acc=lacc)
    !$ACC ENTER DATA CREATE(ecrad_thermodynamics%pressure_hl, ecrad_thermodynamics%temperature_hl) IF(lacc)
    !$ACC ENTER DATA CREATE(ecrad_thermodynamics%pressure_fl, ecrad_thermodynamics%temperature_fl) &
    !$ACC   IF(lacc .AND. ecrad_thermodynamics%rrtm_pass_temppres_fl)

    ! CALL ecrad_gas%allocate(nproma_rad, nlev_rg, use_acc=lacc)
    !$ACC ENTER DATA COPYIN(ecrad_gas%mixing_ratio) IF(lacc)

    ! CALL ecrad_cloud%allocate(nproma_rad, nlev_rg, use_acc=lacc)
    !$ACC ENTER DATA CREATE(ecrad_cloud%q_liq, ecrad_cloud%re_liq, ecrad_cloud%q_ice, ecrad_cloud%re_ice, &
    !$ACC   ecrad_cloud%fraction, ecrad_cloud%overlap_param, ecrad_cloud%fractional_std, &
    !$ACC   ecrad_cloud%inv_cloud_effective_size) IF(lacc)
    ! CALL ecrad_cloud%create_fractional_std(nproma_rad, nlev_rg, 1._wp, use_acc=lacc)
    !$ACC UPDATE DEVICE(ecrad_cloud%fractional_std) IF(lacc)

    ! CALL ecrad_aerosol%allocate_direct(ecrad_conf, nproma_rad, 1, nlev_rg, use_acc=lacc)
    !$ACC ENTER DATA CREATE(ecrad_aerosol%od_sw, ecrad_aerosol%ssa_sw, ecrad_aerosol%g_sw) &
    !$ACC   IF(lacc .AND. ecrad_conf%use_aerosols .AND. ecrad_conf%do_sw)
    !$ACC ENTER DATA CREATE(ecrad_aerosol%od_lw, ecrad_aerosol%ssa_lw, ecrad_aerosol%g_lw) &
    !$ACC   IF(lacc .AND. ecrad_conf%use_aerosols .AND. ecrad_conf%do_lw)
    !$ACC UPDATE DEVICE(ecrad_aerosol%ssa_lw, ecrad_aerosol%g_lw) &
    !$ACC   IF(lacc .AND. ecrad_conf%use_aerosols .AND. ecrad_conf%do_lw)

    ! CALL ecrad_flux%allocate(ecrad_conf, 1, nproma_rad, nlev_rg, use_acc=lacc)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_up, ecrad_flux%lw_dn) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_up_clear, ecrad_flux%lw_dn_clear) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw .AND. ecrad_conf%do_clear )
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_up_band, ecrad_flux%lw_dn_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw .AND. ecrad_conf%do_save_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_up_clear_band, ecrad_flux%lw_dn_clear_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw .AND. ecrad_conf%do_clear .AND. ecrad_conf%do_save_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_derivatives) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw .AND. ecrad_conf%do_lw_derivatives)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_dn_surf_g) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_dn_surf_clear_g) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw .AND. ecrad_conf%do_clear)
    !$ACC ENTER DATA CREATE(ecrad_flux%lw_dn_surf_canopy) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_lw .AND. ecrad_conf%do_canopy_fluxes_lw)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_up, ecrad_flux%sw_dn) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_direct) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_sw_direct)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_up_clear, ecrad_flux%sw_dn_clear) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_clear)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_direct_clear) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_clear .AND. ecrad_conf%do_sw_direct)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_up_band, ecrad_flux%sw_dn_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_save_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_direct_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_sw_direct .AND. ecrad_conf%do_save_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_up_clear_band, ecrad_flux%sw_dn_clear_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_clear .AND. ecrad_conf%do_save_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_direct_clear_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_clear .AND. ecrad_conf%do_sw_direct .AND. &
    !$ACC      ecrad_conf%do_save_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_surf_band, ecrad_flux%sw_dn_direct_surf_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_surface_sw_spectral_flux)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_surf_clear_band, ecrad_flux%sw_dn_direct_surf_clear_band) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_clear)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_diffuse_surf_g,ecrad_flux%sw_dn_direct_surf_g) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_diffuse_surf_clear_g, ecrad_flux%sw_dn_direct_surf_clear_g) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_clear)
    !$ACC ENTER DATA CREATE(ecrad_flux%sw_dn_diffuse_surf_canopy, ecrad_flux%sw_dn_direct_surf_canopy) &
    !$ACC   IF(lacc .AND. ecrad_conf%do_sw .AND. ecrad_conf%do_canopy_fluxes_sw)
    !$ACC ENTER DATA CREATE(ecrad_flux%cloud_cover_lw, ecrad_flux%cloud_cover_sw) &
    !$ACC   IF(lacc)

  END SUBROUTINE

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_acc_deallocation:
  !! The exit data delete of the derived ecrad types that happen in the
  !! subruntines:
  !!   - ecrad_single_level%allocate
  !!   - ecrad_thermodynamics%allocate
  !!   - ecrad_gas%allocate
  !!   - ecrad_cloud%allocate
  !!   - ecrad_cloud%create_fractional_std
  !!   - ecrad_aerosol%allocate_direct
  !!   - ecrad_flux%allocate
  !! are extracted here.
  !! @par Revision History
  !! Initial release by Daniel Hupp, MeteoSwiss, Winterthur (2022-01-11)
  !!
  SUBROUTINE ecrad_acc_deallocation(ecrad_conf, ecrad_single_level, &
      &  ecrad_thermodynamics, ecrad_gas, ecrad_cloud, ecrad_aerosol, &
      &  ecrad_flux, use_acc)

    USE mo_ecrad,                  ONLY: t_ecrad_aerosol_type

    TYPE(t_ecrad_conf),               INTENT(IN)   :: ecrad_conf
    TYPE(t_ecrad_single_level_type),  INTENT(IN)   :: ecrad_single_level
    TYPE(t_ecrad_thermodynamics_type),INTENT(IN)   :: ecrad_thermodynamics
    TYPE(t_ecrad_gas_type),           INTENT(IN)   :: ecrad_gas
    TYPE(t_ecrad_cloud_type),         INTENT(IN)   :: ecrad_cloud
    TYPE(t_ecrad_aerosol_type),       INTENT(IN)   :: ecrad_aerosol
    TYPE(t_ecrad_flux_type),          INTENT(IN)   :: ecrad_flux

    LOGICAL, INTENT(IN), OPTIONAL :: use_acc

    LOGICAL                  :: lacc

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    ! CALL ecrad_single_level%deallocate(use_acc=lacc)
    !$ACC EXIT DATA DELETE(ecrad_single_level%cos_sza) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%cos_sza))
    !$ACC EXIT DATA DELETE(ecrad_single_level%skin_temperature) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%skin_temperature))
    !$ACC EXIT DATA DELETE(ecrad_single_level%sw_albedo) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%sw_albedo))
    !$ACC EXIT DATA DELETE(ecrad_single_level%sw_albedo_direct) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%sw_albedo_direct))
    !$ACC EXIT DATA DELETE(ecrad_single_level%lw_emissivity) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%lw_emissivity))
    !$ACC EXIT DATA DELETE(ecrad_single_level%lw_emission) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%lw_emission))
    !$ACC EXIT DATA DELETE(ecrad_single_level%spectral_solar_scaling) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%spectral_solar_scaling))
    !$ACC EXIT DATA DELETE(ecrad_single_level%iseed) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%iseed))

    ! CALL ecrad_thermodynamics%deallocate(use_acc=lacc)
    !$ACC EXIT DATA DELETE(ecrad_thermodynamics%pressure_hl) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%pressure_hl))
    !$ACC EXIT DATA DELETE(ecrad_thermodynamics%temperature_hl) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%temperature_hl))
    !$ACC EXIT DATA DELETE(ecrad_thermodynamics%pressure_fl) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%pressure_fl))
    !$ACC EXIT DATA DELETE(ecrad_thermodynamics%temperature_fl) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%temperature_fl))
    !$ACC EXIT DATA DELETE(ecrad_thermodynamics%h2o_sat_liq) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%h2o_sat_liq))

    ! CALL ecrad_gas%deallocate(use_acc=lacc)
    !$ACC EXIT DATA DELETE(ecrad_gas%mixing_ratio) IF(lacc .AND. ALLOCATED(ecrad_gas%mixing_ratio))
 
    ! CALL ecrad_cloud%deallocate(use_acc=lacc)
    !$ACC EXIT DATA DELETE(ecrad_cloud%q_liq)                    IF(lacc .AND. ALLOCATED(ecrad_cloud%q_liq))
    !$ACC EXIT DATA DELETE(ecrad_cloud%re_liq)                   IF(lacc .AND. ALLOCATED(ecrad_cloud%re_liq))
    !$ACC EXIT DATA DELETE(ecrad_cloud%q_ice)                    IF(lacc .AND. ALLOCATED(ecrad_cloud%q_ice))
    !$ACC EXIT DATA DELETE(ecrad_cloud%re_ice)                   IF(lacc .AND. ALLOCATED(ecrad_cloud%re_ice))
    !$ACC EXIT DATA DELETE(ecrad_cloud%fraction)                 IF(lacc .AND. ALLOCATED(ecrad_cloud%fraction))
    !$ACC EXIT DATA DELETE(ecrad_cloud%overlap_param)            IF(lacc .AND. ALLOCATED(ecrad_cloud%overlap_param))
    !$ACC EXIT DATA DELETE(ecrad_cloud%fractional_std)           IF(lacc .AND. ALLOCATED(ecrad_cloud%fractional_std))
    !$ACC EXIT DATA DELETE(ecrad_cloud%inv_cloud_effective_size) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_cloud%inv_cloud_effective_size))
    !$ACC EXIT DATA DELETE(ecrad_cloud%inv_inhom_effective_size) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_cloud%inv_inhom_effective_size))

    IF ( ecrad_conf%use_aerosols ) THEN
      ! CALL ecrad_aerosol%deallocate(use_acc=lacc)
      !$ACC EXIT DATA DELETE(ecrad_aerosol%mixing_ratio) IF(lacc .AND. ALLOCATED(ecrad_aerosol%mixing_ratio))
      !$ACC EXIT DATA DELETE(ecrad_aerosol%od_sw)        IF(lacc .AND. ALLOCATED(ecrad_aerosol%od_sw))
      !$ACC EXIT DATA DELETE(ecrad_aerosol%ssa_sw)       IF(lacc .AND. ALLOCATED(ecrad_aerosol%ssa_sw))
      !$ACC EXIT DATA DELETE(ecrad_aerosol%g_sw)         IF(lacc .AND. ALLOCATED(ecrad_aerosol%g_sw))
      !$ACC EXIT DATA DELETE(ecrad_aerosol%od_lw)        IF(lacc .AND. ALLOCATED(ecrad_aerosol%od_lw))
      !$ACC EXIT DATA DELETE(ecrad_aerosol%ssa_lw)       IF(lacc .AND. ALLOCATED(ecrad_aerosol%ssa_lw))
      !$ACC EXIT DATA DELETE(ecrad_aerosol%g_lw)         IF(lacc .AND. ALLOCATED(ecrad_aerosol%g_lw))
    ENDIF

    ! CALL ecrad_flux%deallocate(use_acc=lacc)
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_up)                        IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn)                        IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_up_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up_clear))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_clear))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_up)                        IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn)                        IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_up_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up_clear))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_clear))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct)                 IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_clear) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_clear))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_up_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_up_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_up_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_up_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_surf_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_surf_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_surf_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_surf_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_surf_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_surf_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_clear_band))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn_surf_canopy) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_surf_canopy))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_diffuse_surf_canopy) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_diffuse_surf_canopy))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_surf_canopy) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_canopy))
    !$ACC EXIT DATA DELETE(ecrad_flux%cloud_cover_sw)               IF(lacc .AND. ALLOCATED(ecrad_flux%cloud_cover_sw))
    !$ACC EXIT DATA DELETE(ecrad_flux%cloud_cover_lw)               IF(lacc .AND. ALLOCATED(ecrad_flux%cloud_cover_lw))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_derivatives)               IF(lacc .AND. ALLOCATED(ecrad_flux%lw_derivatives))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn_surf_g)                 IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_surf_g))
    !$ACC EXIT DATA DELETE(ecrad_flux%lw_dn_surf_clear_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_surf_clear_g))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_diffuse_surf_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_diffuse_surf_g))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_surf_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_g))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_diffuse_surf_clear_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_diffuse_surf_clear_g))
    !$ACC EXIT DATA DELETE(ecrad_flux%sw_dn_direct_surf_clear_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_clear_g))
  END SUBROUTINE

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE update_host_pre_ecrad:
  !! The values of the arrays that are need on the CPU for the ecrad calculatione are updated here.
  !! @par Revision History
  !! Initial release by Daniel Hupp, MeteoSwiss, Winterthur (2022-01-11)
  !! 
  SUBROUTINE update_host_pre_ecrad(ecrad_conf, ecrad_single_level, &
      &  ecrad_thermodynamics, ecrad_gas, ecrad_cloud, ecrad_aerosol, &
      &  ecrad_flux, use_acc)

    USE mo_ecrad,                  ONLY: t_ecrad_aerosol_type

    TYPE(t_ecrad_conf),               INTENT(IN)   :: ecrad_conf
    TYPE(t_ecrad_single_level_type),  INTENT(IN)   :: ecrad_single_level
    TYPE(t_ecrad_thermodynamics_type),INTENT(IN)   :: ecrad_thermodynamics
    TYPE(t_ecrad_gas_type),           INTENT(IN)   :: ecrad_gas
    TYPE(t_ecrad_cloud_type),         INTENT(IN)   :: ecrad_cloud
    TYPE(t_ecrad_aerosol_type),       INTENT(IN)   :: ecrad_aerosol
    TYPE(t_ecrad_flux_type),          INTENT(IN)   :: ecrad_flux

    LOGICAL, INTENT(IN), OPTIONAL :: use_acc

    LOGICAL                  :: lacc

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !$ACC UPDATE HOST(ecrad_single_level%cos_sza)                IF(lacc .AND. ALLOCATED(ecrad_single_level%cos_sza))
    !$ACC UPDATE HOST(ecrad_single_level%skin_temperature) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%skin_temperature))
    !$ACC UPDATE HOST(ecrad_single_level%sw_albedo)              IF(lacc .AND. ALLOCATED(ecrad_single_level%sw_albedo))
    !$ACC UPDATE HOST(ecrad_single_level%sw_albedo_direct) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%sw_albedo_direct))
    !$ACC UPDATE HOST(ecrad_single_level%lw_emissivity) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%lw_emissivity))
    !$ACC UPDATE HOST(ecrad_single_level%lw_emission) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%lw_emission))
    !$ACC UPDATE HOST(ecrad_single_level%spectral_solar_scaling) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_single_level%spectral_solar_scaling))
    !$ACC UPDATE HOST(ecrad_single_level%iseed)                  IF(lacc .AND. ALLOCATED(ecrad_single_level%iseed))

    !$ACC UPDATE HOST(ecrad_thermodynamics%pressure_hl)    IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%pressure_hl))
    !$ACC UPDATE HOST(ecrad_thermodynamics%temperature_hl) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%temperature_hl))
    !$ACC UPDATE HOST(ecrad_thermodynamics%pressure_fl)    IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%pressure_fl))
    !$ACC UPDATE HOST(ecrad_thermodynamics%temperature_fl) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%temperature_fl))
    !$ACC UPDATE HOST(ecrad_thermodynamics%h2o_sat_liq)    IF(lacc .AND. ALLOCATED(ecrad_thermodynamics%h2o_sat_liq))

    ! $ACC UPDATE HOST(ecrad_gas%mixing_ratio) IF(lacc .AND. ALLOCATED(ecrad_gas%mixing_ratio))

    !$ACC UPDATE HOST(ecrad_cloud%q_liq)                    IF(lacc .AND. ALLOCATED(ecrad_cloud%q_liq))
    !$ACC UPDATE HOST(ecrad_cloud%re_liq)                   IF(lacc .AND. ALLOCATED(ecrad_cloud%re_liq))
    !$ACC UPDATE HOST(ecrad_cloud%q_ice)                    IF(lacc .AND. ALLOCATED(ecrad_cloud%q_ice))
    !$ACC UPDATE HOST(ecrad_cloud%re_ice)                   IF(lacc .AND. ALLOCATED(ecrad_cloud%re_ice))
    !$ACC UPDATE HOST(ecrad_cloud%fraction)                 IF(lacc .AND. ALLOCATED(ecrad_cloud%fraction))
    !$ACC UPDATE HOST(ecrad_cloud%overlap_param)            IF(lacc .AND. ALLOCATED(ecrad_cloud%overlap_param))
    !$ACC UPDATE HOST(ecrad_cloud%fractional_std)           IF(lacc .AND. ALLOCATED(ecrad_cloud%fractional_std))
    !$ACC UPDATE HOST(ecrad_cloud%inv_cloud_effective_size) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_cloud%inv_cloud_effective_size))
    !$ACC UPDATE HOST(ecrad_cloud%inv_inhom_effective_size) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_cloud%inv_inhom_effective_size))

    !$ACC UPDATE HOST(ecrad_aerosol%mixing_ratio) IF(lacc .AND. ALLOCATED(ecrad_aerosol%mixing_ratio))
    !$ACC UPDATE HOST(ecrad_aerosol%od_sw)        IF(lacc .AND. ALLOCATED(ecrad_aerosol%od_sw))
    !$ACC UPDATE HOST(ecrad_aerosol%ssa_sw)       IF(lacc .AND. ALLOCATED(ecrad_aerosol%ssa_sw))
    !$ACC UPDATE HOST(ecrad_aerosol%g_sw)         IF(lacc .AND. ALLOCATED(ecrad_aerosol%g_sw))
    !$ACC UPDATE HOST(ecrad_aerosol%od_lw)        IF(lacc .AND. ALLOCATED(ecrad_aerosol%od_lw))
    !$ACC UPDATE HOST(ecrad_aerosol%ssa_lw)       IF(lacc .AND. ALLOCATED(ecrad_aerosol%ssa_lw))
    !$ACC UPDATE HOST(ecrad_aerosol%g_lw)         IF(lacc .AND. ALLOCATED(ecrad_aerosol%g_lw))

    !$ACC UPDATE HOST(ecrad_flux%cloud_cover_sw) IF(lacc)
    !$ACC UPDATE HOST(ecrad_flux%cloud_cover_lw) IF(lacc)

  END SUBROUTINE

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE update_device_post_ecrad:
  !! The results of the ecrad computaiton are updated on the GPU.
  !! @par Revision History
  !! Initial release by Daniel Hupp, MeteoSwiss, Winterthur (2022-01-11)
  !! 
  SUBROUTINE update_device_post_ecrad(ecrad_conf, ecrad_flux, use_acc)

    TYPE(t_ecrad_conf),               INTENT(IN)   :: ecrad_conf
    TYPE(t_ecrad_flux_type),          INTENT(IN)   :: ecrad_flux

    LOGICAL, INTENT(IN), OPTIONAL :: use_acc

    LOGICAL                  :: lacc

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF 

    !$ACC UPDATE DEVICE(ecrad_flux%lw_up)                        IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn)                        IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_up_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up_clear))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_clear))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_up)                        IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn)                        IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_up_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up_clear))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_clear)                  IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_clear))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct)                 IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_clear)           IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_clear))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_up_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up_band))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_band))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_up_clear_band)             IF(lacc .AND. ALLOCATED(ecrad_flux%lw_up_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn_clear_band)             IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_up_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_band)                   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_up_clear_band)             IF(lacc .AND. ALLOCATED(ecrad_flux%sw_up_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_clear_band)             IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_band)            IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_surf_band)              IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_surf_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_surf_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_surf_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_surf_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_surf_clear_band) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_clear_band))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn_surf_canopy)            IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_surf_canopy))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_diffuse_surf_canopy) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_diffuse_surf_canopy))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_surf_canopy) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_canopy))
    !$ACC UPDATE DEVICE(ecrad_flux%cloud_cover_sw)               IF(lacc .AND. ALLOCATED(ecrad_flux%cloud_cover_sw))
    !$ACC UPDATE DEVICE(ecrad_flux%cloud_cover_lw)               IF(lacc .AND. ALLOCATED(ecrad_flux%cloud_cover_lw))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_derivatives)               IF(lacc .AND. ALLOCATED(ecrad_flux%lw_derivatives))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn_surf_g)                 IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_surf_g))
    !$ACC UPDATE DEVICE(ecrad_flux%lw_dn_surf_clear_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%lw_dn_surf_clear_g))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_diffuse_surf_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_diffuse_surf_g))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_surf_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_g))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_diffuse_surf_clear_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_diffuse_surf_clear_g))
    !$ACC UPDATE DEVICE(ecrad_flux%sw_dn_direct_surf_clear_g) &
    !$ACC   IF(lacc .AND. ALLOCATED(ecrad_flux%sw_dn_direct_surf_clear_g))


  END SUBROUTINE
#endif

#endif
END MODULE mo_nwp_ecrad_utilities
