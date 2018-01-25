!--------------------------------------------------------------------
!
! Serialization routine for ECHAM surface
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_surface

  USE mo_kind,               ONLY: vp, wp
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memor,    ONLY: t_echam_phy_fiel
  IMPLICIT NONE

  LOGICAL :: writeIn = .TRUE.
  LOGICAL :: writeOut = .TRUE.

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jg, kproma, kbdim, klev, ksfc_type, idx_wtr, &
                             idx_ice, idx_lnd, pdtime, field, pfac_sfc,   &
                             aa, aa_btm, bb, bb_btm, pcpt_tile, nblock,   &
                             pch_tile)
    INTEGER, INTENT(IN)              :: jg
    INTEGER, INTENT(IN)              :: kproma, kbdim
    INTEGER, INTENT(IN)              :: klev, ksfc_type
    INTEGER, INTENT(IN)              :: idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN)              :: pdtime
    TYPE(t_echam_phy_field) ,POINTER :: field
    REAL(wp),INTENT(IN)              :: pfac_sfc(:)
    REAL(wp),INTENT(INOUT)           :: aa(:,:,:,:)
    REAL(wp),INTENT(INOUT)           :: aa_btm(:,:,:,:)
    REAL(wp),INTENT(INOUT)           :: bb(:,:,:)
    REAL(wp),INTENT(INOUT)           :: bb_btm(:,:,:)
    REAL(wp),INTENT(INOUT)           :: pcpt_tile(:,:)
    INTEGER, OPTIONAL,INTENT(IN)     :: nblock
    REAL(wp),OPTIONAL,INTENT(IN)     :: pch_tile(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (writeIn) THEN
    !$ser verbatim call datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim call init()
    !$ser savepoint echam_surface-input jg=jg nblock=nblock date=date
#if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif 
    !$ser data jg=jg                                &
    !$ser&     kproma=kproma                        &
    !$ser&     kbdim=kbdim                          &
    !$ser&     kice=field%kice                      &
    !$ser&     klev=klev                            &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     idx_ice=idx_ice                      &
    !$ser&     idx_lnd=idx_lnd                      &
    !$ser&     pdtime=pdtime                        &
    !$ser&     pfrc=field%frac_tile                 &
    !$ser&     pcfh_tile=field%cfh_tile             &
    !$ser&     pcfm_tile=field%cfm_tile             &
    !$ser&     pfac_sfc=pfac_sfc                    &
    !$ser&     pocu=field%ocu                       &
    !$ser&     pocv=field%ocv                       &
    !$ser&     aa=aa                                &
    !$ser&     aa_btm=aa_btm                        &
    !$ser&     bb=bb                                &
    !$ser&     bb_btm=bb_btm                        &
    !$ser&     pcpt_tile=pcpt_tile                  &
    !$ser&     pqsat_tile=field%qs_sfc_tile         &
    !$ser&     ptsfc_tile=field%ts_tile             &
    !$ser&     plhflx_tile=field%lhflx_tile         &
    !$ser&     pshflx_tile=field%shflx_tile         &
    !$ser&     nblock=nblock                        &
    !$ser&     lsm=field%lsmask                     &
    !$ser&     alake%field%alake                    &
    !$ser&     pu=field%ua                          &
    !$ser&     pv=field%va                          &
    !$ser&     ptemp=field%ta                       &
    !$ser&     pq=field%qtrc                        &
    !$ser&     prsfl=field%rsfl                     &
    !$ser&     prsfc=field%rsfc                     &
    !$ser&     pssfl=field%ssfl                     &
    !$ser&     pssfc=field%ssfc                     &
    !$ser&     rlds=field%rlds                      &
    !$ser&     rlus=field%rlus                      &
    !$ser&     rsds=field%rsds                      &
    !$ser&     rsus=field%rsus                      &
    !$ser&     rvds_dir=field%rvds_dir              &
    !$ser&     rpds_dir=field%rpds_dir              &
    !$ser&     rnds_dir=field%rnds_dir              &
    !$ser&     rvds_dif=field%rvds_dif              &
    !$ser&     rpds_dif=field%rpds_dif              &
    !$ser&     rnds_dif=field%rnds_dif              &
    !$ser&     ps=field%presi_old                   &
    !$ser&     pcosmu0=field%cosmu0                 &
    !$ser&     pch_tile=pch_tile                    &
    !$ser&     pcsat=field%csat                     &
    !$ser&     pcair=field%pcair                    &
    !$ser&     z0m_tile=field%z0m_tile              &
    !$ser&     z0h_lnd=field%z0h_lnd                &
    !$ser&     albvisdir=field%albvisdir            &
    !$ser&     albnirdir=field%albnirdir            &
    !$ser&     albvisdif=field%albvisdif            &
    !$ser&     albnirdif=field%albnirdif            &
    !$ser&     albvisdir_tile=field%albvisdir_tile  &
    !$ser&     albnirdir_tile=field%albnirdir_tile  &
    !$ser&     albvisdif_tile=field%albvisdif_tile  &
    !$ser&     albnirdif_tile=field%albnirdif_tile  &
    !$ser&     albedo=field%albedo                  &
    !$ser&     albedo_tile=field%albedo_tile        &
    !$ser&     rsns_tile=field%swflxfc_tile         &
    !$ser&     rlns_tile=field%lwflxsfc_tile        &
    !$ser&     Tsurf=field%Tsurf                    &
    !$ser&     T1=field%T1                          &
    !$ser&     T2=field%T2                          &
    !$ser&     hi=field%hi                          &
    !$ser&     hs=field%hs                          &
    !$ser&     Qtop=field%Qtop                      &
    !$ser&     Qbot=field%Qbot                      &
    !$ser&     conc=field%conc                      &
    !$ser&     albvisdir_ice=field%albvisdir_ice    &
    !$ser&     albvisdif_ice=field%albnirdir_ice    &
    !$ser&     albnirdir_ice=field%albvisdif_ice    &
    !$ser&     albnirdif_ice=field%albnirdif_ice
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jg, kproma, kbdim, klev, ksfc_type, idx_wtr,     &
                              idx_ice, idx_lnd, pdtime, field, aa, aa_btm, bb, &
                              bb_btm, pcpt_tile, nblock, q_snocpymlt)
    INTEGER, INTENT(IN)              :: jg
    INTEGER, INTENT(IN)              :: kproma, kbdim
    INTEGER, INTENT(IN)              :: klev, ksfc_type
    INTEGER, INTENT(IN)              :: idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN)              :: pdtime
    TYPE(t_echam_phy_field) ,POINTER :: field
    REAL(wp),INTENT(INOUT)           :: aa(:,:,:,:)
    REAL(wp),INTENT(INOUT)           :: aa_btm(:,:,:,:)
    REAL(wp),INTENT(INOUT)           :: bb(:,:,:)
    REAL(wp),INTENT(INOUT)           :: bb_btm(:,:,:)
    REAL(wp),INTENT(INOUT)           :: pcpt_tile(:,:)
    INTEGER, OPTIONAL,INTENT(IN)     :: nblock
    REAL(wp),OPTIONAL,INTENT(INOUT)  :: q_snocpymlt(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim if (writeOut) then
    !$ser verbatim call datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim call init()
    !$ser savepoint echam_surface-output jg=jg nblock=nblock date=date
    !$ser mode write
    !$ser data jg=jg                                &
    !$ser&     kproma=kproma                        &
    !$ser&     kbdim=kbdim                          &
    !$ser&     kice=field%kice                      &
    !$ser&     klev=klev                            &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     idx_ice=idx_ice                      &
    !$ser&     idx_lnd=idx_lnd                      &
    !$ser&     pdtime=pdtime                        &
    !$ser&     aa=aa                                &
    !$ser&     aa_btm=aa_btm                        &
    !$ser&     bb=bb                                &
    !$ser&     bb_btm=bb_btm                        &
    !$ser&     pcpt_tile=pcpt_tile                  &
    !$ser&     pqsat_tile=field%qs_sfc_tile         &
    !$ser&     ptsfc_tile=field%ts_tile             &
    !$ser&     pu_stress_gbm=field%u_stress         &
    !$ser&     pv_stress_gbm=field%v_stress         &
    !$ser&     plhflx_gbm=field%lhflx               &
    !$ser&     pshflx_gbm=field%shflx               &
    !$ser&     pevap_gbm=field%evap                 &
    !$ser&     pu_stress_tile=field%u_stress_tile   &
    !$ser&     pv_stress_tile=field%v_stress_tile   &
    !$ser&     plhflx_tile=field%lhflx_tile         &
    !$ser&     pshflx_tile=field%shflx_tile         &
    !$ser&     pevap_tile=field%evap_tile           &
    !$ser&     nblock=nblock                        &
    !$ser&     rlus=field%rlus                      &
    !$ser&     pcsat=field%csat                     &
    !$ser&     pcair=field%pcair                    &
    !$ser&     q_snocpymlt=q_snocpymlt              &
    !$ser&     z0m_tile=field%z0m_tile              &
    !$ser&     z0h_lnd=field%z0h_lnd                &
    !$ser&     albvisdir=field%albvisdir            &
    !$ser&     albnirdir=field%albnirdir            &
    !$ser&     albvisdif=field%albvisdif            &
    !$ser&     albnirdif=field%albnirdif            &
    !$ser&     albvisdir_tile=field%albvisdir_tile  &
    !$ser&     albnirdir_tile=field%albnirdir_tile  &
    !$ser&     albvisdif_tile=field%albvisdif_tile  &
    !$ser&     albnirdif_tile=field%albnirdif_tile  &
    !$ser&     albedo=field%albedo                  &
    !$ser&     albedo_tile=field%albedo_tile        &
    !$ser&     ptsfc=field%ts                       &
    !$ser&     ptsfc_rad=field%ts_rad               &
    !$ser&     rsns_tile=field%swflxfc_tile         &
    !$ser&     rlns_tile=field%lwflxsfc_tile        &
    !$ser&     lake_ice_frc=field%lake_ice_frc      &
    !$ser&     Tsurf=field%Tsurf                    &
    !$ser&     T1=field%T1                          &
    !$ser&     T2=field%T2                          &
    !$ser&     hi=field%hi                          &
    !$ser&     hs=field%hs                          &
    !$ser&     Qtop=field%Qtop                      &
    !$ser&     Qbot=field%Qbot                      &
    !$ser&     albvisdir_ice=field%albvisdir_ice    &
    !$ser&     albvisdif_ice=field%albnirdir_ice    &
    !$ser&     albnirdir_ice=field%albvisdif_ice    &
    !$ser&     albnirdif_ice=field%albnirdif_ice
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim endif

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_surface
