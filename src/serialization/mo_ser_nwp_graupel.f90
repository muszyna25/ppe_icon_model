!--------------------------------------------------------------------
!
! Serialization routine for NWP GRAUPEL
!
!--------------------------------------------------------------------

MODULE mo_ser_nwp_graupel

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning, finish
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqr, iqs,       &
                                   iqni, iqni_nuc, iqg, iqh, iqnr, iqns,     &
                                   iqng, iqnh, iqnc
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_ser_nml,            ONLY: ser_graupel

  IMPLICIT NONE

  PUBLIC :: serialize_graupel_input
  PUBLIC :: serialize_graupel_output

  CONTAINS

  SUBROUTINE serialize_graupel_input(jg, nproma, nlev, p_metrics, p_prog, &
                                     p_prog_rcf, p_diag, prm_diag, prm_nwp_tend)
    INTEGER, INTENT(IN)    :: jg, nproma, nlev
    TYPE(t_nh_metrics)     , INTENT(in)   :: p_metrics
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog          !<the dyn prog vars
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag)        , INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag        !<the atm phys vars
    TYPE(t_nwp_phy_tend)   , TARGET, INTENT(inout):: prm_nwp_tend    !< atm tend vars

    LOGICAL :: lnumber
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim lnumber = atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .and. iprog_aero > 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_graupel) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_graupel:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_graupel:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( p_metrics%z_ifc ) IF ( ASSOCIATED(p_metrics%z_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( p_metrics%ddqz_z_full ) IF ( ASSOCIATED(p_metrics%ddqz_z_full) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog%w ) IF ( ASSOCIATED(p_prog%w) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog%rho ) IF ( ASSOCIATED(p_prog%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog_rcf%tracer ) IF ( ASSOCIATED(p_prog_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog_rcf%tke ) IF ( ASSOCIATED(p_prog_rcf%tke) )
    !$ser verbatim   !$ACC UPDATE HOST( p_diag%temp ) IF ( ASSOCIATED(p_diag%temp) )
    !$ser verbatim   !$ACC UPDATE HOST( p_diag%pres ) IF ( ASSOCIATED(p_diag%pres) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%cloud_num ) IF ( ASSOCIATED(prm_diag%cloud_num) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%tt_lheat ) IF ( ASSOCIATED(prm_diag%tt_lheat) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%rain_gsp_rate ) IF ( ASSOCIATED(prm_diag%rain_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%snow_gsp_rate ) IF ( ASSOCIATED(prm_diag%snow_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%ice_gsp_rate ) IF ( ASSOCIATED(prm_diag%ice_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%graupel_gsp_rate  ) IF ( ASSOCIATED(prm_diag%graupel_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%qrs_flux ) IF ( ASSOCIATED(prm_diag%qrs_flux) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%rain_gsp ) IF ( ASSOCIATED(prm_diag%rain_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%ice_gsp ) IF ( ASSOCIATED(prm_diag%ice_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%snow_gsp ) IF ( ASSOCIATED(prm_diag%snow_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%hail_gsp ) IF ( ASSOCIATED(prm_diag%hail_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%graupel_gsp ) IF ( ASSOCIATED(prm_diag%graupel_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%prec_gsp ) IF ( ASSOCIATED(prm_diag%prec_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_nwp_tend%ddt_temp_gscp ) IF ( ASSOCIATED(prm_nwp_tend%ddt_temp_gscp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_nwp_tend%ddt_tracer_gscp ) IF ( ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_graupel-input jg=jg nproma=nproma nlev=nlev date=TRIM(date)
#if defined(SERIALIZE_CREATE_REFERENCE)
    !$ser mode       write
#elif defined(SERIALIZE_READ_REFERENCE)
    !$ser mode       read
#elif defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser mode       read-perturb
#else
    !$ser verbatim   CALL finish('nwp_graupel-input', 'SERIALIZATION MODE IS NOT SET')
#endif
    !$ser data nwp_graupel_p_metrics_z_ifc=p_metrics%z_ifc(:,:,:) IF (ASSOCIATED(p_metrics%z_ifc))
    !$ser data nwp_graupel_p_metrics_ddqz_z_full=p_metrics%ddqz_z_full(:,:,:) IF (ASSOCIATED(p_metrics%ddqz_z_full))
    !$ser data nwp_graupel_p_prog_w=p_prog%w(:,:,:) IF (ASSOCIATED(p_prog%w))
    !$ser data nwp_graupel_p_prog_rho=p_prog%rho(:,:,:) IF (ASSOCIATED(p_prog%rho))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqnc=p_prog_rcf%tracer(:,:,:,iqnc) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqnr=p_prog_rcf%tracer(:,:,:,iqnr) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqni=p_prog_rcf%tracer(:,:,:,iqni) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqns=p_prog_rcf%tracer(:,:,:,iqns) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqng=p_prog_rcf%tracer(:,:,:,iqng) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqc=p_prog_rcf%tracer(:,:,:,iqc) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqv=p_prog_rcf%tracer(:,:,:,iqv) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqi=p_prog_rcf%tracer(:,:,:,iqi) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqr=p_prog_rcf%tracer(:,:,:,iqr) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqs=p_prog_rcf%tracer(:,:,:,iqs) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqg=p_prog_rcf%tracer(:,:,:,iqg) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_tke=p_prog_rcf%tke(:,:,:) IF (ASSOCIATED(p_prog_rcf%tke))
    !$ser data nwp_graupel_p_diag_temp=p_diag%temp(:,:,:) IF (ASSOCIATED(p_diag%temp))
    !$ser data nwp_graupel_p_diag_pres=p_diag%pres(:,:,:) IF (ASSOCIATED(p_diag%pres))
    !$ser data nwp_graupel_prm_diag_cloud_num=prm_diag%cloud_num(:,:) IF (ASSOCIATED(prm_diag%cloud_num))
    !$ser data nwp_graupel_prm_diag_tt_lheat=prm_diag%tt_lheat(:,:,:) IF (ASSOCIATED(prm_diag%tt_lheat))
    !$ser data nwp_graupel_prm_diag_rain_gsp_rate=prm_diag%rain_gsp_rate(:,:) IF (ASSOCIATED(prm_diag%rain_gsp_rate))
    !$ser data nwp_graupel_prm_diag_snow_gsp_rate=prm_diag%snow_gsp_rate(:,:) IF (ASSOCIATED(prm_diag%snow_gsp_rate))
    !$ser data nwp_graupel_prm_diag_ice_gsp_rate=prm_diag%ice_gsp_rate(:,:) IF (ASSOCIATED(prm_diag%ice_gsp_rate))
    !$ser data nwp_graupel_prm_diag_graupel_gsp_rate=prm_diag%graupel_gsp_rate (:,:) IF (ASSOCIATED(prm_diag%graupel_gsp_rate))
    !$ser data nwp_graupel_prm_diag_qrs_flux=prm_diag%qrs_flux(:,:,:) IF (ASSOCIATED(prm_diag%qrs_flux))
    !$ser data nwp_graupel_prm_diag_rain_gsp=prm_diag%rain_gsp(:,:) IF (ASSOCIATED(prm_diag%rain_gsp))
    !$ser data nwp_graupel_prm_diag_ice_gsp=prm_diag%ice_gsp(:,:) IF (ASSOCIATED(prm_diag%ice_gsp))
    !$ser data nwp_graupel_prm_diag_snow_gsp=prm_diag%snow_gsp(:,:) IF (ASSOCIATED(prm_diag%snow_gsp))
    !$ser data nwp_graupel_prm_diag_hail_gsp=prm_diag%hail_gsp(:,:) IF (ASSOCIATED(prm_diag%hail_gsp))
    !$ser data nwp_graupel_prm_diag_graupel_gsp=prm_diag%graupel_gsp(:,:) IF (ASSOCIATED(prm_diag%graupel_gsp))
    !$ser data nwp_graupel_prm_diag_prec_gsp=prm_diag%prec_gsp(:,:) IF (ASSOCIATED(prm_diag%prec_gsp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_temp_gscp=prm_nwp_tend%ddt_temp_gscp(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqv=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqc=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqi=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqi) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqr=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqr) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqs=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqs) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
#if defined( _OPENACC )
#if defined(SERIALIZE_READ_REFERENCE) || defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser verbatim CALL warning('GPU:mo_ser_nwp_graupel:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( p_metrics%z_ifc ) IF ( ASSOCIATED(p_metrics%z_ifc) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_metrics%ddqz_z_full ) IF ( ASSOCIATED(p_metrics%ddqz_z_full) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_prog%w ) IF ( ASSOCIATED(p_prog%w) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_prog%rho ) IF ( ASSOCIATED(p_prog%rho) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_prog_rcf%tracer ) IF ( ASSOCIATED(p_prog_rcf%tracer) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_prog_rcf%tke ) IF ( ASSOCIATED(p_prog_rcf%tke) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_diag%temp ) IF ( ASSOCIATED(p_diag%temp) )
    !$ser verbatim !$ACC UPDATE DEVICE( p_diag%pres ) IF ( ASSOCIATED(p_diag%pres) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%cloud_num(:,:) ) IF ( ASSOCIATED(prm_diag%cloud_num) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%tt_lheat ) IF ( ASSOCIATED(prm_diag%tt_lheat) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%rain_gsp_rate ) IF ( ASSOCIATED(prm_diag%rain_gsp_rate) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%snow_gsp_rate ) IF ( ASSOCIATED(prm_diag%snow_gsp_rate) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%ice_gsp_rate ) IF ( ASSOCIATED(prm_diag%ice_gsp_rate) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%graupel_gsp_rate  ) IF ( ASSOCIATED(prm_diag%graupel_gsp_rate) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%qrs_flux ) IF ( ASSOCIATED(prm_diag%qrs_flux) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%rain_gsp ) IF ( ASSOCIATED(prm_diag%rain_gsp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%ice_gsp ) IF ( ASSOCIATED(prm_diag%ice_gsp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%snow_gsp ) IF ( ASSOCIATED(prm_diag%snow_gsp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%hail_gsp ) IF ( ASSOCIATED(prm_diag%hail_gsp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%graupel_gsp ) IF ( ASSOCIATED(prm_diag%graupel_gsp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_diag%prec_gsp ) IF ( ASSOCIATED(prm_diag%prec_gsp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_nwp_tend%ddt_temp_gscp ) IF ( ASSOCIATED(prm_nwp_tend%ddt_temp_gscp) )
    !$ser verbatim !$ACC UPDATE DEVICE( prm_nwp_tend%ddt_tracer_gscp ) IF ( ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_graupel_input

  SUBROUTINE serialize_graupel_output(jg, nproma, nlev, p_metrics, p_prog, &
                                      p_prog_rcf, p_diag, prm_diag, prm_nwp_tend)
    INTEGER, INTENT(IN)    :: jg, nproma, nlev
    TYPE(t_nh_metrics)     , INTENT(in)   :: p_metrics
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog          !<the dyn prog vars
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag)        , INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag        !<the atm phys vars
    TYPE(t_nwp_phy_tend)   , TARGET, INTENT(inout):: prm_nwp_tend    !< atm tend vars

    LOGICAL :: lnumber
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim lnumber = atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .and. iprog_aero > 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_graupel) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_graupel:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_graupel:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( p_metrics%z_ifc ) IF ( ASSOCIATED(p_metrics%z_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( p_metrics%ddqz_z_full ) IF ( ASSOCIATED(p_metrics%ddqz_z_full) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog%w ) IF ( ASSOCIATED(p_prog%w) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog%rho ) IF ( ASSOCIATED(p_prog%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog_rcf%tracer ) IF ( ASSOCIATED(p_prog_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE HOST( p_prog_rcf%tke ) IF ( ASSOCIATED(p_prog_rcf%tke) )
    !$ser verbatim   !$ACC UPDATE HOST( p_diag%temp ) IF ( ASSOCIATED(p_diag%temp) )
    !$ser verbatim   !$ACC UPDATE HOST( p_diag%pres ) IF ( ASSOCIATED(p_diag%pres) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%cloud_num ) IF ( ASSOCIATED(prm_diag%cloud_num) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%tt_lheat ) IF ( ASSOCIATED(prm_diag%tt_lheat) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%rain_gsp_rate ) IF ( ASSOCIATED(prm_diag%rain_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%snow_gsp_rate ) IF ( ASSOCIATED(prm_diag%snow_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%ice_gsp_rate ) IF ( ASSOCIATED(prm_diag%ice_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%graupel_gsp_rate  ) IF ( ASSOCIATED(prm_diag%graupel_gsp_rate) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%qrs_flux ) IF ( ASSOCIATED(prm_diag%qrs_flux) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%rain_gsp ) IF ( ASSOCIATED(prm_diag%rain_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%ice_gsp ) IF ( ASSOCIATED(prm_diag%ice_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%snow_gsp ) IF ( ASSOCIATED(prm_diag%snow_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%hail_gsp ) IF ( ASSOCIATED(prm_diag%hail_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%graupel_gsp ) IF ( ASSOCIATED(prm_diag%graupel_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_diag%prec_gsp ) IF ( ASSOCIATED(prm_diag%prec_gsp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_nwp_tend%ddt_temp_gscp ) IF ( ASSOCIATED(prm_nwp_tend%ddt_temp_gscp) )
    !$ser verbatim   !$ACC UPDATE HOST( prm_nwp_tend%ddt_tracer_gscp ) IF ( ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_graupel-output jg=jg nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data nwp_graupel_p_metrics_z_ifc=p_metrics%z_ifc(:,:,:) IF (ASSOCIATED(p_metrics%z_ifc))
    !$ser data nwp_graupel_p_metrics_ddqz_z_full=p_metrics%ddqz_z_full(:,:,:) IF (ASSOCIATED(p_metrics%ddqz_z_full))
    !$ser data nwp_graupel_p_prog_w=p_prog%w(:,:,:) IF (ASSOCIATED(p_prog%w))
    !$ser data nwp_graupel_p_prog_rho=p_prog%rho(:,:,:) IF (ASSOCIATED(p_prog%rho))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqnc=p_prog_rcf%tracer(:,:,:,iqnc) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqnr=p_prog_rcf%tracer(:,:,:,iqnr) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqni=p_prog_rcf%tracer(:,:,:,iqni) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqns=p_prog_rcf%tracer(:,:,:,iqns) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqng=p_prog_rcf%tracer(:,:,:,iqng) IF (lnumber)
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqc=p_prog_rcf%tracer(:,:,:,iqc) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqv=p_prog_rcf%tracer(:,:,:,iqv) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqi=p_prog_rcf%tracer(:,:,:,iqi) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqr=p_prog_rcf%tracer(:,:,:,iqr) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqs=p_prog_rcf%tracer(:,:,:,iqs) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_iqg=p_prog_rcf%tracer(:,:,:,iqg) IF (ASSOCIATED(p_prog_rcf%tracer))
    !$ser data nwp_graupel_p_prog_rcf_tracer_tke=p_prog_rcf%tke(:,:,:) IF (ASSOCIATED(p_prog_rcf%tke))
    !$ser data nwp_graupel_p_diag_temp=p_diag%temp(:,:,:) IF (ASSOCIATED(p_diag%temp))
    !$ser data nwp_graupel_p_diag_pres=p_diag%pres(:,:,:) IF (ASSOCIATED(p_diag%pres))
    !$ser data nwp_graupel_prm_diag_cloud_num=prm_diag%cloud_num(:,:) IF (ASSOCIATED(prm_diag%cloud_num))
    !$ser data nwp_graupel_prm_diag_tt_lheat=prm_diag%tt_lheat(:,:,:) IF (ASSOCIATED(prm_diag%tt_lheat))
    !$ser data nwp_graupel_prm_diag_rain_gsp_rate=prm_diag%rain_gsp_rate(:,:) IF (ASSOCIATED(prm_diag%rain_gsp_rate))
    !$ser data nwp_graupel_prm_diag_snow_gsp_rate=prm_diag%snow_gsp_rate(:,:) IF (ASSOCIATED(prm_diag%snow_gsp_rate))
    !$ser data nwp_graupel_prm_diag_ice_gsp_rate=prm_diag%ice_gsp_rate(:,:) IF (ASSOCIATED(prm_diag%ice_gsp_rate))
    !$ser data nwp_graupel_prm_diag_graupel_gsp_rate=prm_diag%graupel_gsp_rate (:,:) IF (ASSOCIATED(prm_diag%graupel_gsp_rate))
    !$ser data nwp_graupel_prm_diag_qrs_flux=prm_diag%qrs_flux(:,:,:) IF (ASSOCIATED(prm_diag%qrs_flux))
    !$ser data nwp_graupel_prm_diag_rain_gsp=prm_diag%rain_gsp(:,:) IF (ASSOCIATED(prm_diag%rain_gsp))
    !$ser data nwp_graupel_prm_diag_ice_gsp=prm_diag%ice_gsp(:,:) IF (ASSOCIATED(prm_diag%ice_gsp))
    !$ser data nwp_graupel_prm_diag_snow_gsp=prm_diag%snow_gsp(:,:) IF (ASSOCIATED(prm_diag%snow_gsp))
    !$ser data nwp_graupel_prm_diag_hail_gsp=prm_diag%hail_gsp(:,:) IF (ASSOCIATED(prm_diag%hail_gsp))
    !$ser data nwp_graupel_prm_diag_graupel_gsp=prm_diag%graupel_gsp(:,:) IF (ASSOCIATED(prm_diag%graupel_gsp))
    !$ser data nwp_graupel_prm_diag_prec_gsp=prm_diag%prec_gsp(:,:) IF (ASSOCIATED(prm_diag%prec_gsp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_temp_gscp=prm_nwp_tend%ddt_temp_gscp(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqv=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqc=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqi=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqi) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqr=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqr) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser data nwp_graupel_prm_nwp_tend_ddt_tracer_gscp_iqs=prm_nwp_tend%ddt_tracer_gscp(:,:,:,iqs) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp))
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_graupel_output

END MODULE mo_ser_nwp_graupel

