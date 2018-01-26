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
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field
  IMPLICIT NONE

  LOGICAL :: writeIn = .FALSE.
  LOGICAL :: writeOut = .FALSE.

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jb, nlev, nlevp1, iqv, jg, kproma, kbdim, klev, ksfc_type, idx_wtr, &
                             idx_ice, idx_lnd, pdtime, field, pfac_sfc,   &
                             aa, aa_btm, bb, bb_btm, pcpt_tile, nblock,   &
                             pch_tile)
    INTEGER, INTENT(IN)              :: jb
    INTEGER, INTENT(IN)              :: nlev
    INTEGER, INTENT(IN)              :: nlevp1
    INTEGER, INTENT(IN)              :: iqv
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
    !$ser verbatim call init('echam_surface')
    !$ser savepoint echam_surface-input jb=jb nlev=nlev nlevp1=nlevp1 iqv=iqv &
    !$ser&          date=TRIM(date)
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
    !$ser&     pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     pcfh_tile=field%cfh_tile(:,jb,:)     &
    !$ser&     pcfm_tile=field%cfm_tile(:,jb,:)     &
    !$ser&     pfac_sfc=pfac_sfc                    &
    !$ser&     pocu=field%ocu(:,jb)                 &
    !$ser&     pocv=field%ocv(:,jb)                 &
    !$ser&     aa=aa                                &
    !$ser&     aa_btm=aa_btm                        &
    !$ser&     bb=bb                                &
    !$ser&     bb_btm=bb_btm                        &
    !$ser&     pcpt_tile=pcpt_tile                  &
    !$ser&     pqsat_tile=field%qs_sfc_tile(:,jb,:) &
    !$ser&     ptsfc_tile=field%ts_tile(:,jb,:)     &
    !$ser&     plhflx_tile=field%lhflx_tile(:,jb,:) &
    !$ser&     pshflx_tile=field%shflx_tile(:,jb,:) &
    !$ser&     nblock=nblock                        &
    !$ser&     lsm=field%lsmask(:,jb)               &
    !$ser&     alake%field%alake(:,jb)              &
    !$ser&     pu=field%ua(:,nlev,jb)               &
    !$ser&     pv=field%va(:,nlev,jb)               &
    !$ser&     ptemp=field%ta(:,nlev,jb)            &
    !$ser&     pq=field%qtrc(:,nlev,jb,iqv)         &
    !$ser&     prsfl=field%rsfl(:,jb)               &
    !$ser&     prsfc=field%rsfc(:,jb)               &
    !$ser&     pssfl=field%ssfl(:,jb)               &
    !$ser&     pssfc=field%ssfc(:,jb)               &
    !$ser&     rlds=field%rlds(:,jb)                &
    !$ser&     rlus=field%rlus(:,jb)                &
    !$ser&     rsds=field%rsds(:,jb)                &
    !$ser&     rsus=field%rsus(:,jb)                &
    !$ser&     rvds_dir=field%rvds_dir(:,jb)              &
    !$ser&     rpds_dir=field%rpds_dir(:,jb)              &
    !$ser&     rnds_dir=field%rnds_dir(:,jb)              &
    !$ser&     rvds_dif=field%rvds_dif(:,jb)              &
    !$ser&     rpds_dif=field%rpds_dif(:,jb)              &
    !$ser&     rnds_dif=field%rnds_dif(:,jb)              &
    !$ser&     ps=field%presi_old(:,nlevp1,jb)                  &
    !$ser&     pcosmu0=field%cosmu0(:,jb)                 &
    !$ser&     pch_tile=pch_tile                    &
    !$ser&     pcsat=field%csat(:,jb)                     &
    !$ser&     pcair=field%cair(:,jb)                    &
    !$ser&     z0m_tile=field%z0m_tile(:,jb,:)              &
    !$ser&     z0h_lnd=field%z0h_lnd(:,jb)                &
    !$ser&     albvisdir=field%albvisdir(:,jb)            &
    !$ser&     albnirdir=field%albnirdir(:,jb)            &
    !$ser&     albvisdif=field%albvisdif(:,jb)            &
    !$ser&     albnirdif=field%albnirdif(:,jb)            &
    !$ser&     albvisdir_tile=field%albvisdir_tile(:,jb,:)  &
    !$ser&     albnirdir_tile=field%albnirdir_tile(:,jb,:)  &
    !$ser&     albvisdif_tile=field%albvisdif_tile(:,jb,:)  &
    !$ser&     albnirdif_tile=field%albnirdif_tile(:,jb,:)  &
    !$ser&     albedo=field%albedo(:,jb)                  &
    !$ser&     albedo_tile=field%albedo_tile(:,jb,:)        &
    !$ser&     rsns_tile=field%swflxsfc_tile(:,jb,:)         &
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
    !$ser&     albvisdif_ice=field%albnirdir_ice(:,:,jb)    &
    !$ser&     albnirdir_ice=field%albvisdif_ice(:,:,jb)   &
    !$ser&     albnirdif_ice=field%albnirdif_ice(:,:,jb)
    !NOser verbatim writeIn = .FALSE.
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jb, jg, kproma, kbdim, klev, ksfc_type, idx_wtr,     &
                              idx_ice, idx_lnd, pdtime, field, aa, aa_btm, bb, &
                              bb_btm, pcpt_tile, nblock, q_snocpymlt)
    INTEGER, INTENT(IN)              :: jb
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
    !$ser verbatim call init('echam_surface')
    !$ser savepoint echam_surface-output jb=jb date=TRIM(date)
    !$ser mode write
    !$ser data jg=jg                                        &
    !$ser&     kproma=kproma                                &
    !$ser&     kbdim=kbdim                                  &
    !$ser&     kice=field%kice                              &
    !$ser&     klev=klev                                    &
    !$ser&     ksfc_type=ksfc_type                          &
    !$ser&     idx_wtr=idx_wtr                              &
    !$ser&     idx_ice=idx_ice                              &
    !$ser&     idx_lnd=idx_lnd                              &
    !$ser&     pdtime=pdtime                                &
    !$ser&     aa=aa                                        &
    !$ser&     aa_btm=aa_btm                                &
    !$ser&     bb=bb                                        &
    !$ser&     bb_btm=bb_btm                                &
    !$ser&     pcpt_tile=pcpt_tile                          &
    !$ser&     pqsat_tile=field%qs_sfc_tile(:,jb,:)         &
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
    !$ser&     nblock=nblock                                &
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
    !$ser&     albvisdif_ice=field%albnirdir_ice(:,:,jb)    &
    !$ser&     albnirdir_ice=field%albvisdif_ice(:,:,jb)    &
    !$ser&     albnirdif_ice=field%albnirdif_ice(:,:,jb)
    !NOser verbatim writeOut = .FALSE.
    !$ser verbatim endif

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_surface
