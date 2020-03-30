!--------------------------------------------------------------------
!
! Serialization routine for NWP nwp_turbdiff
!
!--------------------------------------------------------------------

MODULE mo_ser_nwp_tudif

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning, finish
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_lnd_types,      ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqnc, iqni, iqns, iqs, &
                                   iqtke
  USE mo_advection_config,   ONLY: advection_config
  USE mo_ser_nml,            ONLY: ser_turbdiff_interface, ser_turbdiff, ser_vertdiff

  IMPLICIT NONE

  PUBLIC :: serialize_turbdiff_interface_input
  PUBLIC :: serialize_turbdiff_interface_output
  PUBLIC :: serialize_turbdiff_input
  PUBLIC :: serialize_turbdiff_output
  PUBLIC :: serialize_vertdiff_input
  PUBLIC :: serialize_vertdiff_output

  CONTAINS

  SUBROUTINE serialize_turbdiff_interface_input(jg, nproma, nlev, prog, &
                                                prog_rcf, prog_now_rcf, diag, metrics, &
                                                prm_diag, prm_nwp_tend, wtr_prog_now, &
                                                lnd_prog_now, lnd_diag)
    INTEGER, INTENT(IN)    :: jg, nproma, nlev
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog 
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog_rcf
    TYPE(t_nh_prog), TARGET, INTENT(IN) :: prog_now_rcf
    TYPE(t_nh_diag), TARGET, INTENT(INOUT) :: diag
    TYPE(t_nh_metrics), INTENT(IN) :: metrics
    TYPE(t_nwp_phy_diag), INTENT(INOUT) :: prm_diag
    TYPE(t_nwp_phy_tend), TARGET, INTENT(INOUT) :: prm_nwp_tend
    TYPE(t_wtr_prog), INTENT(IN) :: wtr_prog_now
    TYPE(t_lnd_prog), INTENT(IN) :: lnd_prog_now
    TYPE(t_lnd_diag), INTENT(INOUT) :: lnd_diag

    LOGICAL :: ltke_adv
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim ltke_adv = advection_config(jg)%iadv_tke > 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_turbdiff_interface) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_tudif:turbdiff_interface-input','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_tudif:turbdiff_interface-input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(prog%w) IF(ASSOCIATED(prog%w))
    !$ser verbatim   !$ACC UPDATE HOST(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE HOST(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dpres_mc) IF(ASSOCIATED(diag%dpres_mc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE HOST(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE HOST(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE HOST(diag%hdef_ic) IF(ASSOCIATED(diag%hdef_ic))
    !$ser verbatim   !$ACC UPDATE HOST(diag%div_ic) IF(ASSOCIATED(diag%div_ic))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dwdx) IF(ASSOCIATED(diag%dwdx))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dwdy) IF(ASSOCIATED(diag%dwdy))
    !$ser verbatim   !$ACC UPDATE HOST(diag%ddt_tracer_adv) IF(ltke_adv)
    !$ser verbatim   !$ACC UPDATE HOST(metrics%wgtfac_c) IF(ASSOCIATED(metrics%wgtfac_c))
    !$ser verbatim   !$ACC UPDATE HOST(metrics%z_ifc) IF(ASSOCIATED(metrics%z_ifc))
    !$ser verbatim   !$ACC UPDATE HOST(metrics%geopot_agl) IF(ASSOCIATED(metrics%geopot_agl))
    !$ser verbatim   !$ACC UPDATE HOST(metrics%geopot_agl_ifc) IF(ASSOCIATED(metrics%geopot_agl_ifc))
    !$ser verbatim   !$ACC UPDATE HOST(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE HOST(prog_now_rcf%tke) IF(ASSOCIATED(prog_now_rcf%tke))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke) IF(ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_pconv) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_hsh) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser verbatim   !$ACC UPDATE HOST(wtr_prog_now%t_ice) IF(ASSOCIATED(wtr_prog_now%t_ice))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_prog_now%t_g) IF(ASSOCIATED(lnd_prog_now%t_g))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_prog_now%t_snow_t) IF(ASSOCIATED(lnd_prog_now%t_snow_t))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))

#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_turbdiff_interface-input jg=jg nproma=nproma nlev=nlev date=TRIM(date)
#if defined(SERIALIZE_CREATE_REFERENCE)
    !$ser mode       write
#elif defined(SERIALIZE_READ_REFERENCE)
    !$ser mode       read
#elif defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser mode       read-perturb
#else
    !$ser verbatim   CALL finish('nwp_turbdiff_interface-input', 'SERIALIZATION MODE IS NOT SET')
#endif
    !$ser data nwp_turbdiff_interface_prog_w=prog%w(:,:,:) IF (ASSOCIATED(prog%w))
    !$ser data nwp_turbdiff_interface_prog_rho=prog%rho(:,:,:) IF (ASSOCIATED(prog%rho))
    !$ser data nwp_turbdiff_interface_prog_exner=prog%exner(:,:,:) IF (ASSOCIATED(prog%exner))
    !$ser data nwp_turbdiff_interface_prog_dpres_mc=diag%dpres_mc(:,:,:) IF (ASSOCIATED(diag%dpres_mc))
    !$ser data nwp_turbdiff_interface_prog_pres_sfc=diag%pres_sfc(:,:) IF (ASSOCIATED(diag%pres_sfc))
    !$ser data nwp_turbdiff_interface_diag_u=diag%u(:,:,:) IF (ASSOCIATED(diag%u))
    !$ser data nwp_turbdiff_interface_diag_v=diag%v(:,:,:) IF (ASSOCIATED(diag%v))
    !$ser data nwp_turbdiff_interface_diag_temp=diag%temp(:,:,:) IF (ASSOCIATED(diag%temp))
    !$ser data nwp_turbdiff_interface_diag_pres=diag%pres(:,:,:) IF (ASSOCIATED(diag%pres))
    !$ser data nwp_turbdiff_interface_diag_hdef_ic=diag%hdef_ic(:,:,:) IF (ASSOCIATED(diag%hdef_ic))
    !$ser data nwp_turbdiff_interface_diag_div_ic=diag%div_ic(:,:,:) IF (ASSOCIATED(diag%div_ic))
    !$ser data nwp_turbdiff_interface_diag_dwdx=diag%dwdx(:,:,:) IF (ASSOCIATED(diag%dwdx))
    !$ser data nwp_turbdiff_interface_diag_dwdy=diag%dwdy(:,:,:) IF (ASSOCIATED(diag%dwdy))
    !$ser data nwp_turbdiff_interface_diag_ddt_tracer_adv_iqtke=diag%ddt_tracer_adv(:,:,:,iqtke) IF (ltke_adv)
    !$ser data nwp_turbdiff_interface_metrics_wgtfac_c=metrics%wgtfac_c(:,:,:) IF (ASSOCIATED(metrics%wgtfac_c))
    !$ser data nwp_turbdiff_interface_metrics_z_ifc=metrics%z_ifc(:,:,:) IF (ASSOCIATED(metrics%z_ifc))
    !$ser data nwp_turbdiff_interface_metrics_geopot_agl=metrics%geopot_agl(:,:,:) IF (ASSOCIATED(metrics%geopot_agl))
    !$ser data nwp_turbdiff_interface_metrics_geopot_agl_ifc=metrics%geopot_agl_ifc(:,:,:) IF (ASSOCIATED(metrics%geopot_agl_ifc))
    !$ser data nwp_turbdiff_interface_metrics_z_mc=metrics%z_mc(:,:,:) IF (ASSOCIATED(metrics%z_mc))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqv=prog_rcf%tracer(:,:,:,iqv) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqc=prog_rcf%tracer(:,:,:,iqc) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqs=prog_rcf%tracer(:,:,:,iqs) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqi=prog_rcf%tracer(:,:,:,iqi) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_now_rcf=prog_now_rcf%tke(:,:,:) IF (ASSOCIATED(prog_now_rcf%tke))
    !$ser data nwp_turbdiff_interface_prm_diag_gz0=prm_diag%gz0(:,:) IF (ASSOCIATED(prm_diag%gz0))
    !$ser data nwp_turbdiff_interface_prm_diag_tvm=prm_diag%tvm(:,:) IF (ASSOCIATED(prm_diag%tvm))
    !$ser data nwp_turbdiff_interface_prm_diag_tvh=prm_diag%tvh(:,:) IF (ASSOCIATED(prm_diag%tvh))
    !$ser data nwp_turbdiff_interface_prm_diag_tkvm=prm_diag%tkvm(:,:,:) IF (ASSOCIATED(prm_diag%tkvm))
    !$ser data nwp_turbdiff_interface_prm_diag_tkvh=prm_diag%tkvh(:,:,:) IF (ASSOCIATED(prm_diag%tkvh))
    !$ser data nwp_turbdiff_interface_prm_diag_shfl_s=prm_diag%shfl_s(:,:) IF (ASSOCIATED(prm_diag%shfl_s))
    !$ser data nwp_turbdiff_interface_prm_diag_qhfl_s=prm_diag%qhfl_s(:,:) IF (ASSOCIATED(prm_diag%qhfl_s))
    !$ser data nwp_turbdiff_interface_prm_diag_qcfl_s=prm_diag%qcfl_s(:,:) IF (ASSOCIATED(prm_diag%qcfl_s))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tracer_turb_iqv=prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tracer_turb_iqc=prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tracer_turb_iqi=prm_nwp_tend%ddt_tracer_turb(:,:,:,iqi) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tke=prm_nwp_tend%ddt_tke(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tke_pconv=prm_nwp_tend%ddt_tke_pconv(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tke_hsh=prm_nwp_tend%ddt_tke_hsh(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_u_turb=prm_nwp_tend%ddt_u_turb(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_v_turb=prm_nwp_tend%ddt_v_turb(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_temp_turb=prm_nwp_tend%ddt_temp_turb(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_u_sso=prm_nwp_tend%ddt_u_sso(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_v_sso=prm_nwp_tend%ddt_v_sso(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser data nwp_turbdiff_interface_prm_wtr_prog_now_t_ice=wtr_prog_now%t_ice(:,:) IF (ASSOCIATED(wtr_prog_now%t_ice))
    !$ser data nwp_turbdiff_interface_lnd_prg_now_t_g=lnd_prog_now%t_g(:,:) IF (ASSOCIATED(lnd_prog_now%t_g))
    !$ser data nwp_turbdiff_interface_lnd_prg_now_t_snow_t=lnd_prog_now%t_snow_t(:,:,:) IF (ASSOCIATED(lnd_prog_now%t_snow_t))
    !$ser data nwp_turbdiff_interface_lnd_diag_qv_s=lnd_diag%qv_s(:,:) IF (ASSOCIATED(lnd_diag%qv_s))

#if defined(_OPENACC)
#if defined(SERIALIZE_READ_REFERENCE) || defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser verbatim CALL warning('GPU:mo_ser_nwp_tudif:turbdiff_interface-input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%w) IF(ASSOCIATED(prog%w))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%dpres_mc) IF(ASSOCIATED(diag%dpres_mc))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%hdef_ic) IF(ASSOCIATED(diag%hdef_ic))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%div_ic) IF(ASSOCIATED(diag%div_ic))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%dwdx) IF(ASSOCIATED(diag%dwdx))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%dwdy) IF(ASSOCIATED(diag%dwdy))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%ddt_tracer_adv) IF(ltke_adv)
    !$ser verbatim   !$ACC UPDATE DEVICE(metrics%wgtfac_c) IF(ASSOCIATED(metrics%wgtfac_c))
    !$ser verbatim   !$ACC UPDATE DEVICE(metrics%z_ifc) IF(ASSOCIATED(metrics%z_ifc))
    !$ser verbatim   !$ACC UPDATE DEVICE(metrics%geopot_agl) IF(ASSOCIATED(metrics%geopot_agl))
    !$ser verbatim   !$ACC UPDATE DEVICE(metrics%geopot_agl_ifc) IF(ASSOCIATED(metrics%geopot_agl_ifc))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog_now_rcf%tke) IF(ASSOCIATED(prog_now_rcf%tke))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tke) IF(ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tke_pconv) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tke_hsh) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_u_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_v_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser verbatim   !$ACC UPDATE DEVICE(wtr_prog_now%t_ice) IF(ASSOCIATED(wtr_prog_now%t_ice))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_prog_now%t_g) IF(ASSOCIATED(lnd_prog_now%t_g))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_prog_now%t_snow_t) IF(ASSOCIATED(lnd_prog_now%t_snow_t))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))
#endif
#endif

    !$ser verbatim END IF

  END SUBROUTINE serialize_turbdiff_interface_input

  SUBROUTINE serialize_turbdiff_input(jg, jb, nproma, nlev, prog, &
                                      prog_rcf, diag, metrics, &
                                      prm_diag, prm_nwp_tend, &
                                      lnd_prog_now, lnd_diag, z_tvs)
    INTEGER, INTENT(IN)    :: jg, jb, nproma, nlev
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog 
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog_rcf
    TYPE(t_nh_diag), TARGET, INTENT(INOUT) :: diag
    TYPE(t_nh_metrics), INTENT(IN) :: metrics
    TYPE(t_nwp_phy_diag), INTENT(INOUT) :: prm_diag
    TYPE(t_nwp_phy_tend), TARGET, INTENT(INOUT) :: prm_nwp_tend
    TYPE(t_lnd_prog), INTENT(IN) :: lnd_prog_now
    TYPE(t_lnd_diag), INTENT(INOUT) :: lnd_diag
    REAL(wp) :: z_tvs(:,:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_turbdiff) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_tudif:turbdiff-input','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_tudif:turbdiff-input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(prog%w) IF(ASSOCIATED(prog%w))
    !$ser verbatim   !$ACC UPDATE HOST(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE HOST(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dpres_mc) IF(ASSOCIATED(diag%dpres_mc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE HOST(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE HOST(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE HOST(diag%hdef_ic) IF(ASSOCIATED(diag%hdef_ic))
    !$ser verbatim   !$ACC UPDATE HOST(diag%div_ic) IF(ASSOCIATED(diag%div_ic))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dwdx) IF(ASSOCIATED(diag%dwdx))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dwdy) IF(ASSOCIATED(diag%dwdy))
    !$ser verbatim   !$ACC UPDATE HOST(metrics%z_ifc) IF(ASSOCIATED(metrics%z_ifc))
    !$ser verbatim   !$ACC UPDATE HOST(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke) IF(ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_pconv) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_hsh) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_prog_now%t_g) IF(ASSOCIATED(lnd_prog_now%t_g))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))
    !$ser verbatim   !$ACC UPDATE HOST(z_tvs)

#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_turbdiff-input jg=jg jb=jb nproma=nproma nlev=nlev date=TRIM(date)
#if defined(SERIALIZE_CREATE_REFERENCE)
    !$ser mode       write
#elif defined(SERIALIZE_READ_REFERENCE)
    !$ser mode       read
#elif defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser mode       read-perturb
#else
    !$ser verbatim   CALL finish('nwp_turbdiff-input', 'SERIALIZATION MODE IS NOT SET')
#endif
    !$ser data nwp_turbdiff_prog_w=prog%w(:,:,jb) IF (ASSOCIATED(prog%w))
    !$ser data nwp_turbdiff_prog_rho=prog%rho(:,:,jb) IF (ASSOCIATED(prog%rho))
    !$ser data nwp_turbdiff_prog_exner=prog%exner(:,:,jb) IF (ASSOCIATED(prog%exner))
    !$ser data nwp_turbdiff_prog_dpres_mc=diag%dpres_mc(:,:,jb) IF (ASSOCIATED(diag%dpres_mc))
    !$ser data nwp_turbdiff_prog_pres_sfc=diag%pres_sfc(:,jb) IF (ASSOCIATED(diag%pres_sfc))
    !$ser data nwp_turbdiff_diag_u=diag%u(:,:,jb) IF (ASSOCIATED(diag%u))
    !$ser data nwp_turbdiff_diag_v=diag%v(:,:,jb) IF (ASSOCIATED(diag%v))
    !$ser data nwp_turbdiff_diag_temp=diag%temp(:,:,jb) IF (ASSOCIATED(diag%temp))
    !$ser data nwp_turbdiff_diag_pres=diag%pres(:,:,jb) IF (ASSOCIATED(diag%pres))
    !$ser data nwp_turbdiff_diag_hdef_ic=diag%hdef_ic(:,:,jb) IF (ASSOCIATED(diag%hdef_ic))
    !$ser data nwp_turbdiff_diag_div_ic=diag%div_ic(:,:,jb) IF (ASSOCIATED(diag%div_ic))
    !$ser data nwp_turbdiff_diag_dwdx=diag%dwdx(:,:,jb) IF (ASSOCIATED(diag%dwdx))
    !$ser data nwp_turbdiff_diag_dwdy=diag%dwdy(:,:,jb) IF (ASSOCIATED(diag%dwdy))
    !$ser data nwp_turbdiff_metrics_z_ifc=metrics%z_ifc(:,:,jb) IF (ASSOCIATED(metrics%z_ifc))
    !$ser data nwp_turbdiff_prog_rcf_tracer_iqv=prog_rcf%tracer(:,:,jb,iqv) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_prog_rcf_tracer_iqc=prog_rcf%tracer(:,:,jb,iqc) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_prm_diag_gz0=prm_diag%gz0(:,jb) IF (ASSOCIATED(prm_diag%gz0))
    !$ser data nwp_turbdiff_prm_diag_tvm=prm_diag%tvm(:,jb) IF (ASSOCIATED(prm_diag%tvm))
    !$ser data nwp_turbdiff_prm_diag_tvh=prm_diag%tvh(:,jb) IF (ASSOCIATED(prm_diag%tvh))
    !$ser data nwp_turbdiff_prm_diag_tkvm=prm_diag%tkvm(:,:,jb) IF (ASSOCIATED(prm_diag%tkvm))
    !$ser data nwp_turbdiff_prm_diag_tkvh=prm_diag%tkvh(:,:,jb) IF (ASSOCIATED(prm_diag%tkvh))
    !$ser data nwp_turbdiff_prm_diag_shfl_s=prm_diag%shfl_s(:,jb) IF (ASSOCIATED(prm_diag%shfl_s))
    !$ser data nwp_turbdiff_prm_diag_qhfl_s=prm_diag%qhfl_s(:,jb) IF (ASSOCIATED(prm_diag%qhfl_s))
    !$ser data nwp_turbdiff_prm_diag_qcfl_s=prm_diag%qcfl_s(:,jb) IF (ASSOCIATED(prm_diag%qcfl_s))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tracer_turb_iqv=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tracer_turb_iqc=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tke=prm_nwp_tend%ddt_tke(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tke_pconv=prm_nwp_tend%ddt_tke_pconv(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tke_hsh=prm_nwp_tend%ddt_tke_hsh(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_u_turb=prm_nwp_tend%ddt_u_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_v_turb=prm_nwp_tend%ddt_v_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_temp_turb=prm_nwp_tend%ddt_temp_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_u_sso=prm_nwp_tend%ddt_u_sso(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_v_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser data nwp_turbdiff_lnd_prg_now_t_g=lnd_prog_now%t_g(:,jb) IF (ASSOCIATED(lnd_prog_now%t_g))
    !$ser data nwp_turbdiff_lnd_diag_qv_s=lnd_diag%qv_s(:,jb) IF (ASSOCIATED(lnd_diag%qv_s))
    !$ser data nwp_turbdiff_z_tvs=z_tvs(:,:,:)

#ifdef _OPENACC
#if defined(SERIALIZE_READ_REFERENCE) || defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser verbatim CALL warning('GPU:mo_ser_nwp_tudif:turbdiff-input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%w) IF(ASSOCIATED(prog%w))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%dpres_mc) IF(ASSOCIATED(diag%dpres_mc))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%hdef_ic) IF(ASSOCIATED(diag%hdef_ic))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%div_ic) IF(ASSOCIATED(diag%div_ic))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%dwdx) IF(ASSOCIATED(diag%dwdx))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%dwdy) IF(ASSOCIATED(diag%dwdy))
    !$ser verbatim   !$ACC UPDATE DEVICE(metrics%z_ifc) IF(ASSOCIATED(metrics%z_ifc))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tke) IF(ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tke_pconv) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tke_hsh) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_u_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_v_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_prog_now%t_g) IF(ASSOCIATED(lnd_prog_now%t_g))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(z_tvs)
#endif
#endif

    !$ser verbatim END IF

  END SUBROUTINE serialize_turbdiff_input

  SUBROUTINE serialize_vertdiff_input(jg, jb, nproma, nlev, prog, &
                                      prog_rcf, diag, metrics, prm_diag, &
                                      prm_nwp_tend, lnd_prog_now, lnd_diag, &
                                      zvari, zrhon)
    INTEGER, INTENT(IN)    :: jg, jb, nproma, nlev
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog 
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog_rcf
    TYPE(t_nh_diag), TARGET, INTENT(INOUT) :: diag
    TYPE(t_nwp_phy_diag), TARGET, INTENT(INOUT) :: prm_diag
    TYPE(t_nh_metrics), INTENT(IN) :: metrics
    TYPE(t_nwp_phy_tend), TARGET, INTENT(INOUT) :: prm_nwp_tend
    TYPE(t_lnd_prog), INTENT(IN) :: lnd_prog_now
    TYPE(t_lnd_diag), INTENT(INOUT) :: lnd_diag
    REAL(wp) :: zvari(:,:,:)
    REAL(wp) :: zrhon(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_vertdiff) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_tudif:vertdiff-input','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_tudif:vertdiff-input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE HOST(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE HOST(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE HOST(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(metrics%z_ifc) IF(ASSOCIATED(metrics%z_ifc))
    !$ser verbatim   !$ACC UPDATE HOST(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_prog_now%t_g) IF(ASSOCIATED(lnd_prog_now%t_g))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))
    !$ser verbatim   !$ACC UPDATE HOST(zvari)
    !$ser verbatim   !$ACC UPDATE HOST(zrhon)

#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_vertdiff-input jg=jg jb=jb nproma=nproma nlev=nlev date=TRIM(date)
#if defined(SERIALIZE_CREATE_REFERENCE)
    !$ser mode       write
#elif defined(SERIALIZE_READ_REFERENCE)
    !$ser mode       read
#elif defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser mode       read-perturb
#else
    !$ser verbatim   CALL finish('nwp_vertdiff-input', 'SERIALIZATION MODE IS NOT SET')
#endif
    !$ser data nwp_vertdiff_prog_rho=prog%rho(:,:,jb) IF (ASSOCIATED(prog%rho))
    !$ser data nwp_vertdiff_prog_exner=prog%exner(:,:,jb) IF (ASSOCIATED(prog%exner))
    !$ser data nwp_vertdiff_diag_pres_sfc=diag%pres_sfc(:,jb) IF (ASSOCIATED(diag%pres_sfc))
    !$ser data nwp_vertdiff_diag_u=diag%u(:,:,jb) IF (ASSOCIATED(diag%u))
    !$ser data nwp_vertdiff_diag_v=diag%v(:,:,jb) IF (ASSOCIATED(diag%v))
    !$ser data nwp_vertdiff_diag_temp=diag%temp(:,:,jb) IF (ASSOCIATED(diag%temp))
    !$ser data nwp_vertdiff_diag_pres=diag%pres(:,:,jb) IF (ASSOCIATED(diag%pres))
    !$ser data nwp_vertdiff_metrics_z_ifc=metrics%z_ifc(:,:,jb) IF (ASSOCIATED(metrics%z_ifc))
    !$ser data nwp_vertdiff_prog_rcf_iqv=prog_rcf%tracer(:,:,jb,iqv) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prog_rcf_iqc=prog_rcf%tracer(:,:,jb,iqc) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prog_rcf_iqi=prog_rcf%tracer(:,:,jb,iqi) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prog_rcf_iqs=prog_rcf%tracer(:,:,jb,iqs) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_tracer_turb_iqv=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_tracer_turb_iqc=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_tracer_turb_iqi=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_u_turb=prm_nwp_tend%ddt_u_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_v_turb=prm_nwp_tend%ddt_v_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_temp_turb=prm_nwp_tend%ddt_temp_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser data nwp_vertdiff_lnd_prog_now_t_g=lnd_prog_now%t_g(:,jb) IF (ASSOCIATED(lnd_prog_now%t_g))
    !$ser data nwp_vertdiff_lnd_diag_qv_s=lnd_diag%qv_s(:,jb) IF (ASSOCIATED(lnd_diag%qv_s))
    !$ser data nwp_vertdiff_prm_diag_tvm=prm_diag%tvm(:,jb) IF (ASSOCIATED(prm_diag%tvm))
    !$ser data nwp_vertdiff_prm_diag_tvh=prm_diag%tvh(:,jb) IF (ASSOCIATED(prm_diag%tvh))
    !$ser data nwp_vertdiff_prm_diag_tkvm=prm_diag%tkvm(:,:,jb) IF (ASSOCIATED(prm_diag%tkvm))
    !$ser data nwp_vertdiff_prm_diag_tkvh=prm_diag%tkvh(:,:,jb) IF (ASSOCIATED(prm_diag%tkvh))
    !$ser data nwp_vertdiff_prm_diag_shfl_s=prm_diag%shfl_s(:,jb) IF (ASSOCIATED(prm_diag%shfl_s))
    !$ser data nwp_vertdiff_prm_diag_qhfl_s=prm_diag%qhfl_s(:,jb) IF (ASSOCIATED(prm_diag%qhfl_s))
    !$ser data nwp_vertdiff_prm_diag_qcfl_s=prm_diag%qcfl_s(:,jb) IF (ASSOCIATED(prm_diag%qcfl_s))
    !$ser data nwp_vertdiff_zvari_1=zvari(:,:,1)
    !$ser data nwp_vertdiff_zvari_2=zvari(:,:,2)
    !$ser data nwp_vertdiff_zvari_3=zvari(:,:,3)
    !$ser data nwp_vertdiff_zvari_4=zvari(:,:,4)
    !$ser data nwp_vertdiff_zvari_5=zvari(:,:,5)
    !$ser data nwp_vertdiff_zrhon=zrhon(:,:)

#ifdef _OPENACC
#if defined(SERIALIZE_READ_REFERENCE) || defined(SERIALIZE_PERTURB_REFERENCE)
    !$ser verbatim CALL warning('GPU:mo_ser_nwp_tudif:vertdiff-input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE DEVICE(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(metrics%z_ifc) IF(ASSOCIATED(metrics%z_ifc))
    !$ser verbatim   !$ACC UPDATE DEVICE(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_prog_now%t_g) IF(ASSOCIATED(lnd_prog_now%t_g))
    !$ser verbatim   !$ACC UPDATE DEVICE(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))
    !$ser verbatim   !$ACC UPDATE DEVICE(zvari)
    !$ser verbatim   !$ACC UPDATE DEVICE(zrhon)
#endif
#endif

    !$ser verbatim END IF

  END SUBROUTINE serialize_vertdiff_input


  SUBROUTINE serialize_turbdiff_interface_output(jg, nproma, nlev, prog, &
                                                 prog_rcf, diag, prm_diag, prm_nwp_tend, lnd_diag)
    INTEGER, INTENT(IN)    :: jg, nproma, nlev
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog 
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog_rcf
    TYPE(t_nh_diag), TARGET, INTENT(INOUT) :: diag
    TYPE(t_nwp_phy_diag), INTENT(INOUT) :: prm_diag
    TYPE(t_nwp_phy_tend), TARGET, INTENT(INOUT) :: prm_nwp_tend
    TYPE(t_lnd_diag), TARGET, INTENT(INOUT) :: lnd_diag

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_turbdiff_interface) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_tudif:turbdiff_interface-output','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_tudif:turbdiff_interface-output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(prog%w) IF(ASSOCIATED(prog%w))
    !$ser verbatim   !$ACC UPDATE HOST(prog%rho) IF(ASSOCIATED(prog%rho))
    !$ser verbatim   !$ACC UPDATE HOST(prog%exner) IF(ASSOCIATED(prog%exner))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dpres_mc) IF(ASSOCIATED(diag%dpres_mc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres_sfc) IF(ASSOCIATED(diag%pres_sfc))
    !$ser verbatim   !$ACC UPDATE HOST(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE HOST(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE HOST(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE HOST(diag%pres) IF(ASSOCIATED(diag%pres))
    !$ser verbatim   !$ACC UPDATE HOST(diag%hdef_ic) IF(ASSOCIATED(diag%hdef_ic))
    !$ser verbatim   !$ACC UPDATE HOST(diag%div_ic) IF(ASSOCIATED(diag%div_ic))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dwdx) IF(ASSOCIATED(diag%dwdx))
    !$ser verbatim   !$ACC UPDATE HOST(diag%dwdy) IF(ASSOCIATED(diag%dwdy))
    !$ser verbatim   !$ACC UPDATE HOST(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke) IF(ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_pconv) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_hsh) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_sso) IF(ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser verbatim   !$ACC UPDATE HOST(lnd_diag%qv_s) IF(ASSOCIATED(lnd_diag%qv_s))
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_turbdiff_interface-output jg=jg nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data nwp_turbdiff_interface_prog_w=prog%w(:,:,:) IF (ASSOCIATED(prog%w))
    !$ser data nwp_turbdiff_interface_prog_rho=prog%rho(:,:,:) IF (ASSOCIATED(prog%rho))
    !$ser data nwp_turbdiff_interface_prog_exner=prog%exner(:,:,:) IF (ASSOCIATED(prog%exner))
    !$ser data nwp_turbdiff_interface_prog_dpres_mc=diag%dpres_mc(:,:,:) IF (ASSOCIATED(diag%dpres_mc))
    !$ser data nwp_turbdiff_interface_prog_pres_sfc=diag%pres_sfc(:,:) IF (ASSOCIATED(diag%pres_sfc))
    !$ser data nwp_turbdiff_interface_diag_u=diag%u(:,:,:) IF (ASSOCIATED(diag%u))
    !$ser data nwp_turbdiff_interface_diag_v=diag%v(:,:,:) IF (ASSOCIATED(diag%v))
    !$ser data nwp_turbdiff_interface_diag_temp=diag%temp(:,:,:) IF (ASSOCIATED(diag%temp))
    !$ser data nwp_turbdiff_interface_diag_pres=diag%pres(:,:,:) IF (ASSOCIATED(diag%pres))
    !$ser data nwp_turbdiff_interface_diag_hdef_ic=diag%hdef_ic(:,:,:) IF (ASSOCIATED(diag%hdef_ic))
    !$ser data nwp_turbdiff_interface_diag_div_ic=diag%div_ic(:,:,:) IF (ASSOCIATED(diag%div_ic))
    !$ser data nwp_turbdiff_interface_diag_dwdx=diag%dwdx(:,:,:) IF (ASSOCIATED(diag%dwdx))
    !$ser data nwp_turbdiff_interface_diag_dwdy=diag%dwdy(:,:,:) IF (ASSOCIATED(diag%dwdy))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqv=prog_rcf%tracer(:,:,:,iqv) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqc=prog_rcf%tracer(:,:,:,iqc) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqs=prog_rcf%tracer(:,:,:,iqs) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prog_rcf_tracer_iqi=prog_rcf%tracer(:,:,:,iqi) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_turbdiff_interface_prm_diag_gz0=prm_diag%gz0(:,:) IF (ASSOCIATED(prm_diag%gz0))
    !$ser data nwp_turbdiff_interface_prm_diag_tvm=prm_diag%tvm(:,:) IF (ASSOCIATED(prm_diag%tvm))
    !$ser data nwp_turbdiff_interface_prm_diag_tvh=prm_diag%tvh(:,:) IF (ASSOCIATED(prm_diag%tvh))
    !$ser data nwp_turbdiff_interface_prm_diag_tkvm=prm_diag%tkvm(:,:,:) IF (ASSOCIATED(prm_diag%tkvm))
    !$ser data nwp_turbdiff_interface_prm_diag_tkvh=prm_diag%tkvh(:,:,:) IF (ASSOCIATED(prm_diag%tkvh))
    !$ser data nwp_turbdiff_interface_prm_diag_shfl_s=prm_diag%shfl_s(:,:) IF (ASSOCIATED(prm_diag%shfl_s))
    !$ser data nwp_turbdiff_interface_prm_diag_qhfl_s=prm_diag%qhfl_s(:,:) IF (ASSOCIATED(prm_diag%qhfl_s))
    !$ser data nwp_turbdiff_interface_prm_diag_qcfl_s=prm_diag%qcfl_s(:,:) IF (ASSOCIATED(prm_diag%qcfl_s))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tracer_turb_iqv=prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tracer_turb_iqc=prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tracer_turb_iqi=prm_nwp_tend%ddt_tracer_turb(:,:,:,iqi) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tke=prm_nwp_tend%ddt_tke(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tke_pconv=prm_nwp_tend%ddt_tke_pconv(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_pconv))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_tke_hsh=prm_nwp_tend%ddt_tke_hsh(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_u_turb=prm_nwp_tend%ddt_u_turb(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_v_turb=prm_nwp_tend%ddt_v_turb(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_temp_turb=prm_nwp_tend%ddt_temp_turb(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_u_sso=prm_nwp_tend%ddt_u_sso(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_u_sso))
    !$ser data nwp_turbdiff_interface_prm_nwp_tend_ddt_v_sso=prm_nwp_tend%ddt_v_sso(:,:,:) IF (ASSOCIATED(prm_nwp_tend%ddt_v_sso))
    !$ser data nwp_turbdiff_interface_lnd_diag_qv_s=lnd_diag%qv_s(:,:) IF (ASSOCIATED(lnd_diag%qv_s))

    !$ser verbatim END IF

  END SUBROUTINE serialize_turbdiff_interface_output


  SUBROUTINE serialize_turbdiff_output(jg, jb, nproma, nlev, prm_diag, prm_nwp_tend, &
                                       z_tvs, zrhon, zvari)
    INTEGER, INTENT(IN)    :: jg, jb, nproma, nlev
    TYPE(t_nwp_phy_diag), INTENT(INOUT) :: prm_diag
    TYPE(t_nwp_phy_tend), TARGET, INTENT(INOUT) :: prm_nwp_tend
    REAL(wp) :: z_tvs(:,:,:)
    REAL(wp) :: zrhon(:,:)
    REAL(wp) :: zvari(:,:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_turbdiff) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_tudif:turbdiff-output','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_tudif:turbdiff-output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%gz0) IF(ASSOCIATED(prm_diag%gz0))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke) IF(ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tke_hsh) IF(ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE HOST(z_tvs)
    !$ser verbatim   !$ACC UPDATE HOST(zrhon)
    !$ser verbatim   !$ACC UPDATE HOST(zvari)
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_turbdiff-output jg=jg jb=jb nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data nwp_turbdiff_prm_diag_gz0=prm_diag%gz0(:,jb) IF (ASSOCIATED(prm_diag%gz0))
    !$ser data nwp_turbdiff_prm_diag_tvm=prm_diag%tvm(:,jb) IF (ASSOCIATED(prm_diag%tvm))
    !$ser data nwp_turbdiff_prm_diag_tvh=prm_diag%tvh(:,jb) IF (ASSOCIATED(prm_diag%tvh))
    !$ser data nwp_turbdiff_prm_diag_tkvm=prm_diag%tkvm(:,:,jb) IF (ASSOCIATED(prm_diag%tkvm))
    !$ser data nwp_turbdiff_prm_diag_tkvh=prm_diag%tkvh(:,:,jb) IF (ASSOCIATED(prm_diag%tkvh))
    !$ser data nwp_turbdiff_prm_diag_shfl_s=prm_diag%shfl_s(:,jb) IF (ASSOCIATED(prm_diag%shfl_s))
    !$ser data nwp_turbdiff_prm_diag_qhfl_s=prm_diag%qhfl_s(:,jb) IF (ASSOCIATED(prm_diag%qhfl_s))
    !$ser data nwp_turbdiff_prm_diag_qcfl_s=prm_diag%qcfl_s(:,jb) IF (ASSOCIATED(prm_diag%qcfl_s))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tracer_turb_iqv=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tracer_turb_iqc=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tke=prm_nwp_tend%ddt_tke(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_tke))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_tke_hsh=prm_nwp_tend%ddt_tke_hsh(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_tke_hsh))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_u_turb=prm_nwp_tend%ddt_u_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_v_turb=prm_nwp_tend%ddt_v_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser data nwp_turbdiff_prm_nwp_tend_ddt_temp_turb=prm_nwp_tend%ddt_temp_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser data nwp_turbdiff_z_tvs=z_tvs(:,:,:)
    !$ser data nwp_turbdiff_zrhon=zrhon(:,:)
    !$ser data nwp_turbdiff_zvari_1=zvari(:,:,1)
    !$ser data nwp_turbdiff_zvari_2=zvari(:,:,2)
    !$ser data nwp_turbdiff_zvari_3=zvari(:,:,3)
    !$ser data nwp_turbdiff_zvari_4=zvari(:,:,4)
    !$ser data nwp_turbdiff_zvari_5=zvari(:,:,5)

    !$ser verbatim END IF

  END SUBROUTINE serialize_turbdiff_output

  SUBROUTINE serialize_vertdiff_output(jg, jb, nproma, nlev, &
                                       prog_rcf, diag, prm_diag, &
                                       prm_nwp_tend, &
                                       zvari)
    INTEGER, INTENT(IN)    :: jg, jb, nproma, nlev
    TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: prog_rcf
    TYPE(t_nh_diag), TARGET, INTENT(INOUT) :: diag
    TYPE(t_nwp_phy_diag), TARGET, INTENT(INOUT) :: prm_diag
    TYPE(t_nwp_phy_tend), TARGET, INTENT(INOUT) :: prm_nwp_tend
    REAL(wp) :: zvari(:,:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    INTEGER, SAVE :: ser_count = 0

    !$ser verbatim ser_count = ser_count + 1
    !$ser verbatim IF(ser_count <= ser_vertdiff) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_nwp_tudif:vertdiff-output','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim   CALL warning('GPU:mo_ser_nwp_tudif:vertdiff-output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(diag%u) IF(ASSOCIATED(diag%u))
    !$ser verbatim   !$ACC UPDATE HOST(diag%v) IF(ASSOCIATED(diag%v))
    !$ser verbatim   !$ACC UPDATE HOST(diag%temp) IF(ASSOCIATED(diag%temp))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvm) IF(ASSOCIATED(prm_diag%tvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tvh) IF(ASSOCIATED(prm_diag%tvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvm) IF(ASSOCIATED(prm_diag%tkvm))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%tkvh) IF(ASSOCIATED(prm_diag%tkvh))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%shfl_s) IF(ASSOCIATED(prm_diag%shfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qhfl_s) IF(ASSOCIATED(prm_diag%qhfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prm_diag%qcfl_s) IF(ASSOCIATED(prm_diag%qcfl_s))
    !$ser verbatim   !$ACC UPDATE HOST(prog_rcf%tracer) IF(ASSOCIATED(prog_rcf%tracer))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_u_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_v_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_temp_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser verbatim   !$ACC UPDATE HOST(prm_nwp_tend%ddt_tracer_turb) IF(ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser verbatim   !$ACC UPDATE HOST(zvari)

#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_vertdiff-output jg=jg jb=jb nproma=nproma nlev=nlev date=TRIM(date)
    !$ser mode write
    !$ser data nwp_vertdiff_diag_u=diag%u(:,:,jb) IF (ASSOCIATED(diag%u))
    !$ser data nwp_vertdiff_diag_v=diag%v(:,:,jb) IF (ASSOCIATED(diag%v))
    !$ser data nwp_vertdiff_diag_temp=diag%temp(:,:,jb) IF (ASSOCIATED(diag%temp))
    !$ser data nwp_vertdiff_prog_rcf_iqv=prog_rcf%tracer(:,:,jb,iqv) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prog_rcf_iqc=prog_rcf%tracer(:,:,jb,iqc) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prog_rcf_iqi=prog_rcf%tracer(:,:,jb,iqi) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prog_rcf_iqs=prog_rcf%tracer(:,:,jb,iqs) IF (ASSOCIATED(prog_rcf%tracer))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_tracer_turb_iqv=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_tracer_turb_iqc=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_tracer_turb_iqi=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi) IF (ASSOCIATED(prm_nwp_tend%ddt_tracer_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_u_turb=prm_nwp_tend%ddt_u_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_u_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_v_turb=prm_nwp_tend%ddt_v_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_v_turb))
    !$ser data nwp_vertdiff_prm_nwp_tend_ddt_temp_turb=prm_nwp_tend%ddt_temp_turb(:,:,jb) IF (ASSOCIATED(prm_nwp_tend%ddt_temp_turb))
    !$ser data nwp_vertdiff_prm_diag_tvm=prm_diag%tvm(:,jb) IF (ASSOCIATED(prm_diag%tvm))
    !$ser data nwp_vertdiff_prm_diag_tvh=prm_diag%tvh(:,jb) IF (ASSOCIATED(prm_diag%tvh))
    !$ser data nwp_vertdiff_prm_diag_tkvm=prm_diag%tkvm(:,:,jb) IF (ASSOCIATED(prm_diag%tkvm))
    !$ser data nwp_vertdiff_prm_diag_tkvh=prm_diag%tkvh(:,:,jb) IF (ASSOCIATED(prm_diag%tkvh))
    !$ser data nwp_vertdiff_prm_diag_shfl_s=prm_diag%shfl_s(:,jb) IF (ASSOCIATED(prm_diag%shfl_s))
    !$ser data nwp_vertdiff_prm_diag_qhfl_s=prm_diag%qhfl_s(:,jb) IF (ASSOCIATED(prm_diag%qhfl_s))
    !$ser data nwp_vertdiff_prm_diag_qcfl_s=prm_diag%qcfl_s(:,jb) IF (ASSOCIATED(prm_diag%qcfl_s))
    !$ser data nwp_vertdiff_zvari_1=zvari(:,:,1)
    !$ser data nwp_vertdiff_zvari_2=zvari(:,:,2)
    !$ser data nwp_vertdiff_zvari_3=zvari(:,:,3)
    !$ser data nwp_vertdiff_zvari_4=zvari(:,:,4)
    !$ser data nwp_vertdiff_zvari_5=zvari(:,:,5)

    !$ser verbatim END IF

  END SUBROUTINE serialize_vertdiff_output

END MODULE mo_ser_nwp_tudif
