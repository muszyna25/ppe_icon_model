!--------------------------------------------------------------------
!
! Serialization routine for ECHAM surface
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_update_surface

  USE mo_kind,               ONLY: vp, wp
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field
  USE mo_run_config,         ONLY: iqv
  IMPLICIT NONE

  LOGICAL :: writeIn = .FALSE.
  LOGICAL :: writeOut = .FALSE.
  LOGICAL :: serializeStepIn = .TRUE.
  LOGICAL :: serializeStepOut = .TRUE.
  LOGICAL, PARAMETER :: singleStepIn = .TRUE.
  LOGICAL, PARAMETER :: singleStepOut = .TRUE.
  INTEGER, PARAMETER :: singleBlock = 7

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jb, jg, jcs, kproma, kbdim, klev, klevp1, ksfc_type, &
                             idx_wtr, idx_ice, idx_lnd, pdtime, field,        &
                             pcfh_tile, pcfm_tile, pfac_sfc, aa, aa_btm, bb,  &
                             bb_btm, pcpt_tile, pqsat_tile, nblock,      &
                             pco2, pch_tile)
    INTEGER, INTENT(IN)              :: jb, jg, jcs, kproma, kbdim, klev, klevp1, ksfc_type
    INTEGER, INTENT(IN)              :: idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN)              :: pdtime
    TYPE(t_echam_phy_field),POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT)              :: pcfh_tile(:,:)
    REAL(wp),INTENT(INOUT)              :: pcfm_tile(:,:)
    REAL(wp),INTENT(INOUT)              :: pfac_sfc(:)
    REAL(wp),INTENT(INOUT)              :: aa(:,:,:,:)
    REAL(wp),INTENT(INOUT)              :: aa_btm(:,:,:,:)
    REAL(wp),INTENT(INOUT)              :: bb(:,:,:)
    REAL(wp),INTENT(INOUT)              :: bb_btm(:,:,:)
    REAL(wp),INTENT(INOUT)              :: pcpt_tile(:,:)
    REAL(wp),INTENT(INOUT)              :: pqsat_tile(:,:)
    INTEGER,INTENT(IN)                  :: nblock
    REAL(wp),INTENT(INOUT)              :: pco2(:)
    REAL(wp),INTENT(INOUT)              :: pch_tile(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date


    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeIn = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepIn .and. writeIn) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_update_surface')
    !$ser savepoint echam_update_surface-input jb=jb jg=jg jcs=jcs kproma=kproma kbdim=kbdim &
    !$ser&          kice=field%kice klev=klev klevp1=klevp1 ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd pdtime=pdtime iqv=iqv date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data pfrc=field%frac_tile(:,jb,:)                 &
    !$ser&     pcfh_tile=pcfh_tile                          &
    !$ser&     pcfm_tile=pcfm_tile                          &
    !$ser&     pfac_sfc=pfac_sfc                            &
    !$ser&     pocu=field%ocu(:,jb)                         &
    !$ser&     pocv=field%ocv(:,jb)                         &
    !$ser&     aa=aa                                        &
    !$ser&     aa_btm=aa_btm                                &
    !$ser&     bb=bb                                        &
    !$ser&     bb_btm=bb_btm                                &
    !$ser&     pcpt_tile=pcpt_tile                          &
    !$ser&     pqsat_tile=pqsat_tile                        &
    !$ser&     ptsfc_tile=field%ts_tile(:,jb,:)             &
    !$ser&     plhflx_tile=field%lhflx_tile(:,jb,:)         &
    !$ser&     pshflx_tile=field%shflx_tile(:,jb,:)         &
    !$ser&     lsm=field%lsmask(:,jb)                       &
    !$ser&     alake%field%alake(:,jb)                      &
    !$ser&     pu=field%ua(:,klev,jb)                       &
    !$ser&     pv=field%va(:,klev,jb)                       &
    !$ser&     ptemp=field%ta(:,klev,jb)                    &
    !$ser&     pq=field%qtrc(:,klev,jb,iqv)                 &
    !$ser&     pco2=pco2                                    &
    !$ser&     prsfl=field%rsfl(:,jb)                       &
    !$ser&     prsfc=field%rsfc(:,jb)                       &
    !$ser&     pssfl=field%ssfl(:,jb)                       &
    !$ser&     pssfc=field%ssfc(:,jb)                       &
    !$ser&     rlds=field%rlds(:,jb)                        &
    !$ser&     rlus=field%rlus(:,jb)                        &
    !$ser&     rsds=field%rsds(:,jb)                        &
    !$ser&     rsus=field%rsus(:,jb)                        &
    !$ser&     rvds_dir=field%rvds_dir(:,jb)                &
    !$ser&     rpds_dir=field%rpds_dir(:,jb)                &
    !$ser&     rnds_dir=field%rnds_dir(:,jb)                &
    !$ser&     rvds_dif=field%rvds_dif(:,jb)                &
    !$ser&     rpds_dif=field%rpds_dif(:,jb)                &
    !$ser&     rnds_dif=field%rnds_dif(:,jb)                &
    !$ser&     ps=field%presi_old(:,klevp1,jb)              &
    !$ser&     pcosmu0=field%cosmu0(:,jb)                   &
    !$ser&     pch_tile=pch_tile                            &
    !$ser&     pcsat=field%csat(:,jb)                       &
    !$ser&     pcair=field%cair(:,jb)                       &
    !$ser&     z0m_tile=field%z0m_tile(:,jb,:)              &
    !$ser&     z0h_lnd=field%z0h_lnd(:,jb)                  &
    !$ser&     albvisdir=field%albvisdir(:,jb)              &
    !$ser&     albnirdir=field%albnirdir(:,jb)              &
    !$ser&     albvisdif=field%albvisdif(:,jb)              &
    !$ser&     albnirdif=field%albnirdif(:,jb)              &
    !$ser&     albvisdir_tile=field%albvisdir_tile(:,jb,:)  &
    !$ser&     albnirdir_tile=field%albnirdir_tile(:,jb,:)  &
    !$ser&     albvisdif_tile=field%albvisdif_tile(:,jb,:)  &
    !$ser&     albnirdif_tile=field%albnirdif_tile(:,jb,:)  &
    !$ser&     albedo=field%albedo(:,jb)                    &
    !$ser&     albedo_tile=field%albedo_tile(:,jb,:)        &
    !$ser&     pco2_flux_tile=field%co2_flux_tile(:,jb,:)   &
    !$ser&     rsns_tile=field%swflxsfc_tile(:,jb,:)        &
    !$ser&     rlns_tile=field%lwflxsfc_tile(:,jb,:)        &
    !$ser&     Tsurf=field%Tsurf(:,:,jb)                    &
    !$ser&     T1=field%T1(:,:,jb)                          &
    !$ser&     T2=field%T2(:,:,jb)                          &
    !$ser&     hi=field%hi(:,:,jb)                          &
    !$ser&     hs=field%hs(:,:,jb)                          &
    !$ser&     Qtop=field%Qtop(:,:,jb)                      &
    !$ser&     Qbot=field%Qbot(:,:,jb)                      &
    !$ser&     conc=field%conc(:,:,jb)                      &
    !$ser&     albvisdir_ice=field%albvisdir_ice(:,:,jb)    &
    !$ser&     albnirdir_ice=field%albnirdir_ice(:,:,jb)    &
    !$ser&     albvisdif_ice=field%albvisdif_ice(:,:,jb)    &
    !$ser&     albnirdif_ice=field%albnirdif_ice(:,:,jb)
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim IF (singleStepIn) THEN
    !$ser verbatim   serializeStepIn = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jb, jg, jcs, kproma, kbdim, klev, ksfc_type, idx_wtr,     &
                              idx_ice, idx_lnd, pdtime, field, aa, aa_btm, bb, &
                              bb_btm, pcpt_tile, pqsat_tile, q_snocpymlt)
    INTEGER, INTENT(IN)       :: jb, jg, jcs, kproma, kbdim, klev, ksfc_type
    INTEGER, INTENT(IN)       :: idx_wtr, idx_ice, idx_lnd
    REAL(wp), INTENT(IN)      :: pdtime
    TYPE(t_echam_phy_field) ,POINTER, INTENT(INOUT) :: field
    REAL(wp), INTENT(INOUT)  :: aa(:,:,:,:)
    REAL(wp), INTENT(INOUT)  :: aa_btm(:,:,:,:)
    REAL(wp), INTENT(INOUT)  :: bb(:,:,:)
    REAL(wp), INTENT(INOUT)  :: bb_btm(:,:,:)
    REAL(wp), INTENT(INOUT)  :: pcpt_tile(:,:)
    REAL(wp), INTENT(INOUT)  :: pqsat_tile(:,:)
    REAL(wp), INTENT(INOUT)  :: q_snocpymlt(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeOut = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepOut .and. writeOut) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_update_surface')
    !$ser savepoint echam_update_surface-output jb=jb jg=jg jcs=jcs kproma=kproma kbdim=kbdim &
    !$ser&          kice=field%kice klev=klev ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd pdtime=pdtime iqv=iqv date=TRIM(date)
    !$ser mode write
    !$ser data aa=aa                                        &
    !$ser&     aa_btm=aa_btm                                &
    !$ser&     bb=bb                                        &
    !$ser&     bb_btm=bb_btm                                &
    !$ser&     pcpt_tile=pcpt_tile                          &
    !$ser&     pqsat_tile=pqsat_tile                        &
    !$ser&     ptsfc_tile=field%ts_tile(:,jb,:)             &
    !$ser&     pu_stress_gbm=field%u_stress(:,jb)           &
    !$ser&     pv_stress_gbm=field%v_stress(:,jb)           &
    !$ser&     plhflx_gbm=field%lhflx(:,jb)                 &
    !$ser&     pshflx_gbm=field%shflx(:,jb)                 &
    !$ser&     pevap_gbm=field%evap(:,jb)                   &
    !$ser&     pu_stress_tile=field%u_stress_tile(:,jb,:)   &
    !$ser&     pv_stress_tile=field%v_stress_tile(:,jb,:)   &
    !$ser&     plhflx_tile=field%lhflx_tile(:,jb,:)         &
    !$ser&     pshflx_tile=field%shflx_tile(:,jb,:)         &
    !$ser&     pevap_tile=field%evap_tile(:,jb,:)           &
    !$ser&     pco2nat=field%fco2nat(:,jb)                  &
    !$ser&     rlus=field%rlus(:,jb)                        &
    !$ser&     pcsat=field%csat(:,jb)                       &
    !$ser&     pcair=field%cair(:,jb)                       &
    !$ser&     q_snocpymlt=q_snocpymlt                      &
    !$ser&     z0m_tile=field%z0m_tile(:,jb,:)              &
    !$ser&     z0h_lnd=field%z0h_lnd(:,jb)                  &
    !$ser&     albvisdir=field%albvisdir(:,jb)              &
    !$ser&     albnirdir=field%albnirdir(:,jb)              &
    !$ser&     albvisdif=field%albvisdif(:,jb)              &
    !$ser&     albnirdif=field%albnirdif(:,jb)              &
    !$ser&     albvisdir_tile=field%albvisdir_tile(:,jb,:)  &
    !$ser&     albnirdir_tile=field%albnirdir_tile(:,jb,:)  &
    !$ser&     albvisdif_tile=field%albvisdif_tile(:,jb,:)  &
    !$ser&     albnirdif_tile=field%albnirdif_tile(:,jb,:)  &
    !$ser&     albedo=field%albedo(:,jb)                    &
    !$ser&     albedo_tile=field%albedo_tile(:,jb,:)        &
    !$ser&     co2_flux_tile=field%co2_flux_tile(:,jb,:)    &
    !$ser&     ptsfc=field%ts(:,jb)                         &
    !$ser&     ptsfc_rad=field%ts_rad(:,jb)                 &
    !$ser&     rsns_tile=field%swflxsfc_tile(:,jb,:)        &
    !$ser&     rlns_tile=field%lwflxsfc_tile(:,jb,:)        &
    !$ser&     lake_ice_frc=field%lake_ice_frc(:,jb)        &
    !$ser&     Tsurf=field%Tsurf(:,:,jb)                    &
    !$ser&     T1=field%T1(:,:,jb)                          &
    !$ser&     T2=field%T2(:,:,jb)                          &
    !$ser&     hi=field%hi(:,:,jb)                          &
    !$ser&     hs=field%hs(:,:,jb)                          &
    !$ser&     Qtop=field%Qtop(:,:,jb)                      &
    !$ser&     Qbot=field%Qbot(:,:,jb)                      &
    !$ser&     albvisdir_ice=field%albvisdir_ice(:,:,jb)    &
    !$ser&     albnirdir_ice=field%albnirdir_ice(:,:,jb)    &
    !$ser&     albvisdif_ice=field%albvisdif_ice(:,:,jb)    &
    !$ser&     albnirdif_ice=field%albnirdif_ice(:,:,jb)
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_update_surface
