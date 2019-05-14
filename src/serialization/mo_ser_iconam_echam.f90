!--------------------------------------------------------------------
!
! Serialization routine for ECHAM-ICONAM interface
!
!--------------------------------------------------------------------

MODULE mo_ser_iconam_echam

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  IMPLICIT NONE

  PUBLIC :: serialize_iconam_input
  PUBLIC :: serialize_iconam_output

  CONTAINS

  SUBROUTINE serialize_iconam_input(jg, field, tend,                                       &
                                    pt_int_state, p_metrics, pt_prog_old, pt_prog_old_rcf, &
                                    pt_prog_new, pt_prog_new_rcf, pt_diag)
    INTEGER, INTENT(IN)    :: jg
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend
    TYPE(t_int_state)     , INTENT(IN), TARGET      :: pt_int_state
    TYPE(t_nh_metrics)    , INTENT(IN), TARGET      :: p_metrics
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_old
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_old_rcf
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_new
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_new_rcf
    TYPE(t_nh_diag)       , INTENT(INOUT), TARGET   :: pt_diag

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .TRUE.
    REAL(wp), POINTER :: c_lin_e(:,:,:), rbf_vec_coeff_c(:,:,:,:)
    INTEGER, POINTER  :: rbf_vec_idx_c(:,:,:), rbf_vec_blk_c(:,:,:)

    c_lin_e => pt_int_state%c_lin_e
    rbf_vec_idx_c => pt_int_state%rbf_vec_idx_c
    rbf_vec_blk_c => pt_int_state%rbf_vec_blk_c
    rbf_vec_coeff_c => pt_int_state%rbf_vec_coeff_c

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_iconam:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_iconam:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%vor ) IF( ASSOCIATED(field%vor) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%tv ) IF( ASSOCIATED(field%tv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_new ) IF( ASSOCIATED(field%presm_new) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%dz ) IF( ASSOCIATED(field%dz) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mh2o ) IF( ASSOCIATED(field%mh2o) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mdry ) IF( ASSOCIATED(field%mdry) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%xref ) IF( ASSOCIATED(field%xref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%omega ) IF( ASSOCIATED(field%omega) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_new ) IF( ASSOCIATED(field%presi_new) )
    !$ser verbatim   !$ACC UPDATE HOST( field%clon ) IF( ASSOCIATED(field%clon) )
    !$ser verbatim   !$ACC UPDATE HOST( field%clat ) IF( ASSOCIATED(field%clat) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mtrc ) IF( ASSOCIATED(field%mtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mtrcvi ) IF( ASSOCIATED(field%mtrcvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mh2ovi ) IF( ASSOCIATED(field%mh2ovi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mairvi ) IF( ASSOCIATED(field%mairvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mdryvi ) IF( ASSOCIATED(field%mdryvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mrefvi ) IF( ASSOCIATED(field%mrefvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%seaice ) IF( ASSOCIATED(field%seaice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%siced ) IF( ASSOCIATED(field%siced) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftof ) IF( ASSOCIATED(field%sftof) )
    !$ser verbatim   !$ACC UPDATE HOST( field%conc ) IF( ASSOCIATED(field%conc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%hi ) IF( ASSOCIATED(field%hi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0_rt ) IF( ASSOCIATED(field%cosmu0_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%daylght_frc ) IF( ASSOCIATED(field%daylght_frc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%daylght_frc_rt ) IF( ASSOCIATED(field%daylght_frc_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua ) IF( ASSOCIATED(tend%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va ) IF( ASSOCIATED(tend%va) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta ) IF( ASSOCIATED(tend%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_dyn ) IF( ASSOCIATED(tend%ua_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_dyn ) IF( ASSOCIATED(tend%va_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_dyn ) IF( ASSOCIATED(tend%ta_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc ) IF( ASSOCIATED(tend%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_dyn ) IF( ASSOCIATED(tend%qtrc_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%mtrcvi_phy ) IF( ASSOCIATED(tend%mtrcvi_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%mtrc_phy ) IF( ASSOCIATED(tend%mtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%c_lin_e )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%rbf_vec_idx_c )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%rbf_vec_blk_c )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%rbf_vec_coeff_c )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old%vn ) IF( ASSOCIATED(pt_prog_old%vn) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old%theta_v ) IF( ASSOCIATED(pt_prog_old%theta_v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old%exner ) IF( ASSOCIATED(pt_prog_old%exner) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old_rcf%tracer ) IF( ASSOCIATED(pt_prog_old_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%vn ) IF( ASSOCIATED(pt_prog_new%vn) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%rho ) IF( ASSOCIATED(pt_prog_new%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%w ) IF( ASSOCIATED(pt_prog_new%w) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%exner ) IF( ASSOCIATED(pt_prog_new%exner) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%theta_v ) IF( ASSOCIATED(pt_prog_new%theta_v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new_rcf%tracer ) IF( ASSOCIATED(pt_prog_new_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%u ) IF( ASSOCIATED(pt_diag%u) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%v ) IF( ASSOCIATED(pt_diag%v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%vor ) IF( ASSOCIATED(pt_diag%vor) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp ) IF( ASSOCIATED(pt_diag%temp) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%tempv ) IF( ASSOCIATED(pt_diag%tempv) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp_ifc ) IF( ASSOCIATED(pt_diag%temp_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres ) IF( ASSOCIATED(pt_diag%pres) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres_ifc ) IF( ASSOCIATED(pt_diag%pres_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%ddt_tracer_adv ) IF( ASSOCIATED(pt_diag%ddt_tracer_adv) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%ddt_vn_phy ) IF( ASSOCIATED(pt_diag%ddt_vn_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%exner_pr ) IF( ASSOCIATED(pt_diag%exner_pr) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%ddt_exner_phy ) IF( ASSOCIATED(pt_diag%ddt_exner_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%exner_dyn_incr ) IF( ASSOCIATED(pt_diag%exner_dyn_incr) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_iconam-input jg=jg date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_iconam_field_ua=field%ua IF (ASSOCIATED(field%ua))
    !$ser data echam_iconam_field_va=field%va IF (ASSOCIATED(field%va))
    !$ser data echam_iconam_field_vor=field%vor IF (ASSOCIATED(field%vor))
    !$ser data echam_iconam_field_ta=field%ta IF (ASSOCIATED(field%ta))
    !$ser data echam_iconam_field_tv=field%tv IF (ASSOCIATED(field%tv))
    !$ser data echam_iconam_field_presm_old=field%presm_old IF (ASSOCIATED(field%presm_old))
    !$ser data echam_iconam_field_presm_new=field%presm_new IF (ASSOCIATED(field%presm_new))
    !$ser data echam_iconam_field_rho=field%rho IF (ASSOCIATED(field%rho))
    !$ser data echam_iconam_field_mair=field%mair IF (ASSOCIATED(field%mair))
    !$ser data echam_iconam_field_dz=field%dz IF (ASSOCIATED(field%dz))
    !$ser data echam_iconam_field_mh2o=field%mh2o IF (ASSOCIATED(field%mh2o))
    !$ser data echam_iconam_field_mdry=field%mdry IF (ASSOCIATED(field%mdry))
    !$ser data echam_iconam_field_mref=field%mref IF (ASSOCIATED(field%mref))
    !$ser data echam_iconam_field_xref=field%xref IF (ASSOCIATED(field%xref))
    !$ser data echam_iconam_field_omega=field%omega IF (ASSOCIATED(field%omega))
    !$ser data echam_iconam_field_presi_old=field%presi_old IF (ASSOCIATED(field%presi_old))
    !$ser data echam_iconam_field_presi_new=field%presi_new IF (ASSOCIATED(field%presi_new))
    !$ser data echam_iconam_field_clon=field%clon IF (ASSOCIATED(field%clon))
    !$ser data echam_iconam_field_clat=field%clat IF (ASSOCIATED(field%clat))
    !$ser data echam_iconam_field_mtrc=field%mtrc IF (ASSOCIATED(field%mtrc))
    !$ser data echam_iconam_field_qtrc=field%qtrc IF (ASSOCIATED(field%qtrc))
    !$ser data echam_iconam_field_mtrcvi=field%mtrcvi IF (ASSOCIATED(field%mtrcvi))
    !$ser data echam_iconam_field_mh2ovi=field%mh2ovi IF (ASSOCIATED(field%mh2ovi))
    !$ser data echam_iconam_field_mairvi=field%mairvi IF (ASSOCIATED(field%mairvi))
    !$ser data echam_iconam_field_mdryvi=field%mdryvi IF (ASSOCIATED(field%mdryvi))
    !$ser data echam_iconam_field_mrefvi=field%mrefvi IF (ASSOCIATED(field%mrefvi))
    !$ser data echam_iconam_field_ts_tile=field%ts_tile IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_iconam_field_seaice=field%seaice IF (ASSOCIATED(field%seaice))
    !$ser data echam_iconam_field_siced=field%siced IF (ASSOCIATED(field%siced))
    !$ser data echam_iconam_field_sftof=field%sftof IF (ASSOCIATED(field%sftof))
    !$ser data echam_iconam_field_conc=field%conc IF (ASSOCIATED(field%conc))
    !$ser data echam_iconam_field_hi=field%hi IF (ASSOCIATED(field%hi))
    !$ser data echam_iconam_field_cosmu0=field%cosmu0 IF (ASSOCIATED(field%cosmu0))
    !$ser data echam_iconam_field_cosmu0_rt=field%cosmu0_rt IF (ASSOCIATED(field%cosmu0_rt))
    !$ser data echam_iconam_field_daylght_frc=field%daylght_frc IF (ASSOCIATED(field%daylght_frc))
    !$ser data echam_iconam_field_daylght_frc_rt=field%daylght_frc_rt IF (ASSOCIATED(field%daylght_frc_rt))
    !$ser data echam_iconam_tend_qtrc_phy=tend%qtrc_phy IF (ASSOCIATED(tend%qtrc_phy))
    !$ser data echam_iconam_tend_ua=tend%ua IF (ASSOCIATED(tend%ua))
    !$ser data echam_iconam_tend_va=tend%va IF (ASSOCIATED(tend%va))
    !$ser data echam_iconam_tend_ta=tend%ta IF (ASSOCIATED(tend%ta))
    !$ser data echam_iconam_tend_ua_phy=tend%ua_phy IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_iconam_tend_va_phy=tend%va_phy IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_iconam_tend_ta_phy=tend%ta_phy IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_iconam_tend_ua_dyn=tend%ua_dyn IF (ASSOCIATED(tend%ua_dyn))
    !$ser data echam_iconam_tend_va_dyn=tend%va_dyn IF (ASSOCIATED(tend%va_dyn))
    !$ser data echam_iconam_tend_ta_dyn=tend%ta_dyn IF (ASSOCIATED(tend%ta_dyn))
    !$ser data echam_iconam_tend_qtrc=tend%qtrc IF (ASSOCIATED(tend%qtrc))
    !$ser data echam_iconam_tend_qtrc_dyn=tend%qtrc_dyn IF (ASSOCIATED(tend%qtrc_dyn))
    !$ser data echam_iconam_tend_mtrcvi_phy=tend%mtrcvi_phy IF (ASSOCIATED(tend%mtrcvi_phy))
    !$ser data echam_iconam_tend_mtrc_phy=tend%mtrc_phy IF (ASSOCIATED(tend%mtrc_phy))
    !$ser data echam_iconam_pt_int_state_c_lin_e=c_lin_e
    !$ser data echam_iconam_pt_int_state_rbf_vec_idx_c=rbf_vec_idx_c
    !$ser data echam_iconam_pt_int_state_rbf_vec_blk_c=rbf_vec_blk_c
    !$ser data echam_iconam_pt_int_state_rbf_vec_coeff_c=rbf_vec_coeff_c
    !$ser data echam_iconam_pt_prog_old_vn=pt_prog_old%vn IF (ASSOCIATED(pt_prog_old%vn))
    !$ser data echam_iconam_pt_prog_old_theta_v=pt_prog_old%theta_v IF (ASSOCIATED(pt_prog_old%theta_v))
    !$ser data echam_iconam_pt_prog_old_exner=pt_prog_old%exner IF (ASSOCIATED(pt_prog_old%exner))
    !$ser data echam_iconam_pt_prog_old_rcf_tracer=pt_prog_old_rcf%tracer IF (ASSOCIATED(pt_prog_old_rcf%tracer))
    !$ser data echam_iconam_pt_prog_new_vn=pt_prog_new%vn IF (ASSOCIATED(pt_prog_new%vn))
    !$ser data echam_iconam_pt_prog_new_rho=pt_prog_new%rho IF (ASSOCIATED(pt_prog_new%rho))
    !$ser data echam_iconam_pt_prog_new_w=pt_prog_new%w IF (ASSOCIATED(pt_prog_new%w))
    !$ser data echam_iconam_pt_prog_new_exner=pt_prog_new%exner IF (ASSOCIATED(pt_prog_new%exner))
    !$ser data echam_iconam_pt_prog_new_theta_v=pt_prog_new%theta_v IF (ASSOCIATED(pt_prog_new%theta_v))
    !$ser data echam_iconam_pt_prog_new_rcf%tracer=pt_prog_new_rcf%tracer IF (ASSOCIATED(pt_prog_new_rcf%tracer))
    !$ser data echam_iconam_pt_diag_u=pt_diag%u IF (ASSOCIATED(pt_diag%u))
    !$ser data echam_iconam_pt_diag_v=pt_diag%v IF (ASSOCIATED(pt_diag%v))
    !$ser data echam_iconam_pt_diag_vor=pt_diag%vor IF (ASSOCIATED(pt_diag%vor))
    !$ser data echam_iconam_pt_diag_temp=pt_diag%temp IF (ASSOCIATED(pt_diag%temp))
    !$ser data echam_iconam_pt_diag_tempv=pt_diag%tempv IF (ASSOCIATED(pt_diag%tempv))
    !$ser data echam_iconam_pt_diag_temp_ifc=pt_diag%temp_ifc IF (ASSOCIATED(pt_diag%temp_ifc))
    !$ser data echam_iconam_pt_diag_pres=pt_diag%pres IF (ASSOCIATED(pt_diag%pres))
    !$ser data echam_iconam_pt_diag_pres_ifc=pt_diag%pres_ifc IF (ASSOCIATED(pt_diag%pres_ifc))
    !$ser data echam_iconam_pt_diag_ddt_tracer_adv=pt_diag%ddt_tracer_adv IF (ASSOCIATED(pt_diag%ddt_tracer_adv))
    !$ser data echam_iconam_pt_diag_ddt_vn_phy=pt_diag%ddt_vn_phy IF (ASSOCIATED(pt_diag%ddt_vn_phy))
    !$ser data echam_iconam_pt_diag_ddt_exner_phy=pt_diag%ddt_exner_phy IF (ASSOCIATED(pt_diag%ddt_exner_phy))
    !$ser data echam_iconam_pt_diag_exner_dyn_incr=pt_diag%exner_dyn_incr IF (ASSOCIATED(pt_diag%exner_dyn_incr))
    !$ser data echam_iconam_pt_diag_exner_pr=pt_diag%exner_pr IF (ASSOCIATED(pt_diag%exner_pr))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_iconam:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%vor ) IF( ASSOCIATED(field%vor) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tv ) IF( ASSOCIATED(field%tv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_new ) IF( ASSOCIATED(field%presm_new) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%dz ) IF( ASSOCIATED(field%dz) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mh2o ) IF( ASSOCIATED(field%mh2o) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mdry ) IF( ASSOCIATED(field%mdry) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%xref ) IF( ASSOCIATED(field%xref) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%omega ) IF( ASSOCIATED(field%omega) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_new ) IF( ASSOCIATED(field%presi_new) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%clon ) IF( ASSOCIATED(field%clon) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%clat ) IF( ASSOCIATED(field%clat) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mtrc ) IF( ASSOCIATED(field%mtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mtrcvi ) IF( ASSOCIATED(field%mtrcvi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mh2ovi ) IF( ASSOCIATED(field%mh2ovi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mairvi ) IF( ASSOCIATED(field%mairvi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mdryvi ) IF( ASSOCIATED(field%mdryvi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mrefvi ) IF( ASSOCIATED(field%mrefvi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%seaice ) IF( ASSOCIATED(field%seaice) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%siced ) IF( ASSOCIATED(field%siced) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%sftof ) IF( ASSOCIATED(field%sftof) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%conc ) IF( ASSOCIATED(field%conc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%hi ) IF( ASSOCIATED(field%hi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%cosmu0_rt ) IF( ASSOCIATED(field%cosmu0_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%daylght_frc ) IF( ASSOCIATED(field%daylght_frc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%daylght_frc_rt ) IF( ASSOCIATED(field%daylght_frc_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua ) IF( ASSOCIATED(tend%ua) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va ) IF( ASSOCIATED(tend%va) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta ) IF( ASSOCIATED(tend%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_dyn ) IF( ASSOCIATED(tend%ua_dyn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_dyn ) IF( ASSOCIATED(tend%va_dyn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_dyn ) IF( ASSOCIATED(tend%ta_dyn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc ) IF( ASSOCIATED(tend%qtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_dyn ) IF( ASSOCIATED(tend%qtrc_dyn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%mtrcvi_phy ) IF( ASSOCIATED(tend%mtrcvi_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%mtrc_phy ) IF( ASSOCIATED(tend%mtrc_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_int_state%c_lin_e )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_int_state%rbf_vec_idx_c )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_int_state%rbf_vec_blk_c )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_int_state%rbf_vec_coeff_c )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_old%vn ) IF( ASSOCIATED(pt_prog_old%vn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_old%theta_v ) IF( ASSOCIATED(pt_prog_old%theta_v) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_old%exner ) IF( ASSOCIATED(pt_prog_old%exner) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_old_rcf%tracer ) IF( ASSOCIATED(pt_prog_old_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_new%vn ) IF( ASSOCIATED(pt_prog_new%vn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_new%rho ) IF( ASSOCIATED(pt_prog_new%rho) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_new%w ) IF( ASSOCIATED(pt_prog_new%w) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_new%exner ) IF( ASSOCIATED(pt_prog_new%exner) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_new%theta_v ) IF( ASSOCIATED(pt_prog_new%theta_v) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_new_rcf%tracer ) IF( ASSOCIATED(pt_prog_new_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%u ) IF( ASSOCIATED(pt_diag%u) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%v ) IF( ASSOCIATED(pt_diag%v) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%vor ) IF( ASSOCIATED(pt_diag%vor) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%temp ) IF( ASSOCIATED(pt_diag%temp) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%tempv ) IF( ASSOCIATED(pt_diag%tempv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%temp_ifc ) IF( ASSOCIATED(pt_diag%temp_ifc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%pres ) IF( ASSOCIATED(pt_diag%pres) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%pres_ifc ) IF( ASSOCIATED(pt_diag%pres_ifc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%ddt_tracer_adv ) IF( ASSOCIATED(pt_diag%ddt_tracer_adv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%ddt_vn_phy ) IF( ASSOCIATED(pt_diag%ddt_vn_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%exner_pr ) IF( ASSOCIATED(pt_diag%exner_pr) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%ddt_exner_phy ) IF( ASSOCIATED(pt_diag%ddt_exner_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%exner_dyn_incr ) IF( ASSOCIATED(pt_diag%exner_dyn_incr) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_iconam_input

  SUBROUTINE serialize_iconam_output(jg, field, tend,           &
                                     pt_int_state, p_metrics, pt_prog_old, pt_prog_old_rcf, &
                                     pt_prog_new, pt_prog_new_rcf, pt_diag)
    INTEGER, INTENT(IN)    :: jg
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend
    TYPE(t_int_state)     , INTENT(IN), TARGET      :: pt_int_state
    TYPE(t_nh_metrics)    , INTENT(IN), TARGET      :: p_metrics
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_old
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_old_rcf
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_new
    TYPE(t_nh_prog)       , INTENT(INOUT), TARGET   :: pt_prog_new_rcf
    TYPE(t_nh_diag)       , INTENT(INOUT), TARGET   :: pt_diag

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .TRUE.
    REAL(wp), POINTER :: c_lin_e(:,:,:), rbf_vec_coeff_c(:,:,:,:)
    INTEGER, POINTER  :: rbf_vec_idx_c(:,:,:), rbf_vec_blk_c(:,:,:)

    c_lin_e => pt_int_state%c_lin_e
    rbf_vec_idx_c => pt_int_state%rbf_vec_idx_c
    rbf_vec_blk_c => pt_int_state%rbf_vec_blk_c
    rbf_vec_coeff_c => pt_int_state%rbf_vec_coeff_c

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_iconam:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_iconam:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%vor ) IF( ASSOCIATED(field%vor) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%tv ) IF( ASSOCIATED(field%tv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_new ) IF( ASSOCIATED(field%presm_new) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%dz ) IF( ASSOCIATED(field%dz) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mh2o ) IF( ASSOCIATED(field%mh2o) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mdry ) IF( ASSOCIATED(field%mdry) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%xref ) IF( ASSOCIATED(field%xref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%omega ) IF( ASSOCIATED(field%omega) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_new ) IF( ASSOCIATED(field%presi_new) )
    !$ser verbatim   !$ACC UPDATE HOST( field%clon ) IF( ASSOCIATED(field%clon) )
    !$ser verbatim   !$ACC UPDATE HOST( field%clat ) IF( ASSOCIATED(field%clat) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mtrc ) IF( ASSOCIATED(field%mtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mtrcvi ) IF( ASSOCIATED(field%mtrcvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mh2ovi ) IF( ASSOCIATED(field%mh2ovi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mairvi ) IF( ASSOCIATED(field%mairvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mdryvi ) IF( ASSOCIATED(field%mdryvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mrefvi ) IF( ASSOCIATED(field%mrefvi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%seaice ) IF( ASSOCIATED(field%seaice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%siced ) IF( ASSOCIATED(field%siced) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftof ) IF( ASSOCIATED(field%sftof) )
    !$ser verbatim   !$ACC UPDATE HOST( field%conc ) IF( ASSOCIATED(field%conc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%hi ) IF( ASSOCIATED(field%hi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0_rt ) IF( ASSOCIATED(field%cosmu0_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%daylght_frc ) IF( ASSOCIATED(field%daylght_frc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%daylght_frc_rt ) IF( ASSOCIATED(field%daylght_frc_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua ) IF( ASSOCIATED(tend%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va ) IF( ASSOCIATED(tend%va) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta ) IF( ASSOCIATED(tend%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_dyn ) IF( ASSOCIATED(tend%ua_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_dyn ) IF( ASSOCIATED(tend%va_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_dyn ) IF( ASSOCIATED(tend%ta_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc ) IF( ASSOCIATED(tend%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_dyn ) IF( ASSOCIATED(tend%qtrc_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%mtrcvi_phy ) IF( ASSOCIATED(tend%mtrcvi_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%mtrc_phy ) IF( ASSOCIATED(tend%mtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%c_lin_e )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%rbf_vec_idx_c )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%rbf_vec_blk_c )
    !$ser verbatim   !$ACC UPDATE HOST( pt_int_state%rbf_vec_coeff_c )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old%vn ) IF( ASSOCIATED(pt_prog_old%vn) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old%theta_v ) IF( ASSOCIATED(pt_prog_old%theta_v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old%exner ) IF( ASSOCIATED(pt_prog_old%exner) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_old_rcf%tracer ) IF( ASSOCIATED(pt_prog_old_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%vn ) IF( ASSOCIATED(pt_prog_new%vn) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%rho ) IF( ASSOCIATED(pt_prog_new%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%w ) IF( ASSOCIATED(pt_prog_new%w) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%exner ) IF( ASSOCIATED(pt_prog_new%exner) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new%theta_v ) IF( ASSOCIATED(pt_prog_new%theta_v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_new_rcf%tracer ) IF( ASSOCIATED(pt_prog_new_rcf%tracer) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%u ) IF( ASSOCIATED(pt_diag%u) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%v ) IF( ASSOCIATED(pt_diag%v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%vor ) IF( ASSOCIATED(pt_diag%vor) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp ) IF( ASSOCIATED(pt_diag%temp) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%tempv ) IF( ASSOCIATED(pt_diag%tempv) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp_ifc ) IF( ASSOCIATED(pt_diag%temp_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres ) IF( ASSOCIATED(pt_diag%pres) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres_ifc ) IF( ASSOCIATED(pt_diag%pres_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%ddt_tracer_adv ) IF( ASSOCIATED(pt_diag%ddt_tracer_adv) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%ddt_vn_phy ) IF( ASSOCIATED(pt_diag%ddt_vn_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%exner_pr ) IF( ASSOCIATED(pt_diag%exner_pr) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%ddt_exner_phy ) IF( ASSOCIATED(pt_diag%ddt_exner_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%exner_dyn_incr ) IF( ASSOCIATED(pt_diag%exner_dyn_incr) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_iconam-output jg=jg date=TRIM(date) 
    !$ser mode write
    !$ser data echam_iconam_field_ua=field%ua IF (ASSOCIATED(field%ua))
    !$ser data echam_iconam_field_va=field%va IF (ASSOCIATED(field%va))
    !$ser data echam_iconam_field_vor=field%vor IF (ASSOCIATED(field%vor))
    !$ser data echam_iconam_field_ta=field%ta IF (ASSOCIATED(field%ta))
    !$ser data echam_iconam_field_tv=field%tv IF (ASSOCIATED(field%tv))
    !$ser data echam_iconam_field_presm_old=field%presm_old IF (ASSOCIATED(field%presm_old))
    !$ser data echam_iconam_field_presm_new=field%presm_new IF (ASSOCIATED(field%presm_new))
    !$ser data echam_iconam_field_rho=field%rho IF (ASSOCIATED(field%rho))
    !$ser data echam_iconam_field_mair=field%mair IF (ASSOCIATED(field%mair))
    !$ser data echam_iconam_field_dz=field%dz IF (ASSOCIATED(field%dz))
    !$ser data echam_iconam_field_mh2o=field%mh2o IF (ASSOCIATED(field%mh2o))
    !$ser data echam_iconam_field_mdry=field%mdry IF (ASSOCIATED(field%mdry))
    !$ser data echam_iconam_field_mref=field%mref IF (ASSOCIATED(field%mref))
    !$ser data echam_iconam_field_xref=field%xref IF (ASSOCIATED(field%xref))
    !$ser data echam_iconam_field_omega=field%omega IF (ASSOCIATED(field%omega))
    !$ser data echam_iconam_field_presi_old=field%presi_old IF (ASSOCIATED(field%presi_old))
    !$ser data echam_iconam_field_presi_new=field%presi_new IF (ASSOCIATED(field%presi_new))
    !$ser data echam_iconam_field_clon=field%clon IF (ASSOCIATED(field%clon))
    !$ser data echam_iconam_field_clat=field%clat IF (ASSOCIATED(field%clat))
    !$ser data echam_iconam_field_mtrc=field%mtrc IF (ASSOCIATED(field%mtrc))
    !$ser data echam_iconam_field_qtrc=field%qtrc IF (ASSOCIATED(field%qtrc))
    !$ser data echam_iconam_field_mtrcvi=field%mtrcvi IF (ASSOCIATED(field%mtrcvi))
    !$ser data echam_iconam_field_mh2ovi=field%mh2ovi IF (ASSOCIATED(field%mh2ovi))
    !$ser data echam_iconam_field_mairvi=field%mairvi IF (ASSOCIATED(field%mairvi))
    !$ser data echam_iconam_field_mdryvi=field%mdryvi IF (ASSOCIATED(field%mdryvi))
    !$ser data echam_iconam_field_mrefvi=field%mrefvi IF (ASSOCIATED(field%mrefvi))
    !$ser data echam_iconam_field_ts_tile=field%ts_tile IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_iconam_field_seaice=field%seaice IF (ASSOCIATED(field%seaice))
    !$ser data echam_iconam_field_siced=field%siced IF (ASSOCIATED(field%siced))
    !$ser data echam_iconam_field_sftof=field%sftof IF (ASSOCIATED(field%sftof))
    !$ser data echam_iconam_field_conc=field%conc IF (ASSOCIATED(field%conc))
    !$ser data echam_iconam_field_hi=field%hi IF (ASSOCIATED(field%hi))
    !$ser data echam_iconam_field_cosmu0=field%cosmu0 IF (ASSOCIATED(field%cosmu0))
    !$ser data echam_iconam_field_cosmu0_rt=field%cosmu0_rt IF (ASSOCIATED(field%cosmu0_rt))
    !$ser data echam_iconam_field_daylght_frc=field%daylght_frc IF (ASSOCIATED(field%daylght_frc))
    !$ser data echam_iconam_field_daylght_frc_rt=field%daylght_frc_rt IF (ASSOCIATED(field%daylght_frc_rt))
    !$ser data echam_iconam_tend_qtrc_phy=tend%qtrc_phy IF (ASSOCIATED(tend%qtrc_phy))
    !$ser data echam_iconam_tend_ua=tend%ua IF (ASSOCIATED(tend%ua))
    !$ser data echam_iconam_tend_va=tend%va IF (ASSOCIATED(tend%va))
    !$ser data echam_iconam_tend_ta=tend%ta IF (ASSOCIATED(tend%ta))
    !$ser data echam_iconam_tend_ua_phy=tend%ua_phy IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_iconam_tend_va_phy=tend%va_phy IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_iconam_tend_ta_phy=tend%ta_phy IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_iconam_tend_ua_dyn=tend%ua_dyn IF (ASSOCIATED(tend%ua_dyn))
    !$ser data echam_iconam_tend_va_dyn=tend%va_dyn IF (ASSOCIATED(tend%va_dyn))
    !$ser data echam_iconam_tend_ta_dyn=tend%ta_dyn IF (ASSOCIATED(tend%ta_dyn))
    !$ser data echam_iconam_tend_qtrc=tend%qtrc IF (ASSOCIATED(tend%qtrc))
    !$ser data echam_iconam_tend_qtrc_dyn=tend%qtrc_dyn IF (ASSOCIATED(tend%qtrc_dyn))
    !$ser data echam_iconam_tend_mtrcvi_phy=tend%mtrcvi_phy IF (ASSOCIATED(tend%mtrcvi_phy))
    !$ser data echam_iconam_tend_mtrc_phy=tend%mtrc_phy IF (ASSOCIATED(tend%mtrc_phy))
    !$ser data echam_iconam_pt_int_state_c_lin_e=c_lin_e
    !$ser data echam_iconam_pt_int_state_rbf_vec_idx_c=rbf_vec_idx_c
    !$ser data echam_iconam_pt_int_state_rbf_vec_blk_c=rbf_vec_blk_c
    !$ser data echam_iconam_pt_int_state_rbf_vec_coeff_c=rbf_vec_coeff_c
    !$ser data echam_iconam_pt_prog_old_vn=pt_prog_old%vn IF (ASSOCIATED(pt_prog_old%vn))
    !$ser data echam_iconam_pt_prog_old_rcf_tracer=pt_prog_old_rcf%tracer IF (ASSOCIATED(pt_prog_old_rcf%tracer))
    !$ser data echam_iconam_pt_prog_new_vn=pt_prog_new%vn IF (ASSOCIATED(pt_prog_new%vn))
    !$ser data echam_iconam_pt_prog_new_rho=pt_prog_new%rho IF (ASSOCIATED(pt_prog_new%rho))
    !$ser data echam_iconam_pt_prog_new_w=pt_prog_new%w IF (ASSOCIATED(pt_prog_new%w))
    !$ser data echam_iconam_pt_prog_new_exner=pt_prog_new%exner IF (ASSOCIATED(pt_prog_new%exner))
    !$ser data echam_iconam_pt_prog_new_theta_v=pt_prog_new%theta_v IF (ASSOCIATED(pt_prog_new%theta_v))
    !$ser data echam_iconam_pt_prog_new_rcf%tracer=pt_prog_new_rcf%tracer IF (ASSOCIATED(pt_prog_new_rcf%tracer))
    !$ser data echam_iconam_pt_diag_u=pt_diag%u IF (ASSOCIATED(pt_diag%u))
    !$ser data echam_iconam_pt_diag_v=pt_diag%v IF (ASSOCIATED(pt_diag%v))
    !$ser data echam_iconam_pt_diag_vor=pt_diag%vor IF (ASSOCIATED(pt_diag%vor))
    !$ser data echam_iconam_pt_diag_temp=pt_diag%temp IF (ASSOCIATED(pt_diag%temp))
    !$ser data echam_iconam_pt_diag_tempv=pt_diag%tempv IF (ASSOCIATED(pt_diag%tempv))
    !$ser data echam_iconam_pt_diag_temp_ifc=pt_diag%temp_ifc IF (ASSOCIATED(pt_diag%temp_ifc))
    !$ser data echam_iconam_pt_diag_pres=pt_diag%pres IF (ASSOCIATED(pt_diag%pres))
    !$ser data echam_iconam_pt_diag_pres_ifc=pt_diag%pres_ifc IF (ASSOCIATED(pt_diag%pres_ifc))
    !$ser data echam_iconam_pt_diag_ddt_tracer_adv=pt_diag%ddt_tracer_adv IF (ASSOCIATED(pt_diag%ddt_tracer_adv))
    !$ser data echam_iconam_pt_diag_ddt_vn_phy=pt_diag%ddt_vn_phy IF (ASSOCIATED(pt_diag%ddt_vn_phy))
    !$ser data echam_iconam_pt_diag_ddt_exner_phy=pt_diag%ddt_exner_phy IF (ASSOCIATED(pt_diag%ddt_exner_phy))
    !$ser data echam_iconam_pt_diag_exner_dyn_incr=pt_diag%exner_dyn_incr IF (ASSOCIATED(pt_diag%exner_dyn_incr))
    !$ser data echam_iconam_pt_diag_exner_pr=pt_diag%exner_pr IF (ASSOCIATED(pt_diag%exner_pr))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_iconam_output

END MODULE mo_ser_iconam_echam
