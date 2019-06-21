!--------------------------------------------------------------------
!
! Serialization routine for ECHAM VDF
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_vdf

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqt
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  ! Full interface serialization
  PUBLIC :: serialize_vdf_input
  PUBLIC :: serialize_vdf_output
  ! Additional serialization check points
  PUBLIC :: serialize_vdf_chk_A_output
  PUBLIC :: serialize_vdf_chk_B_output
  PUBLIC :: serialize_vdf_chk_C_output
  PUBLIC :: serialize_vdf_chk_D_output
  PUBLIC :: serialize_vdf_chk_E_output
  PUBLIC :: serialize_vdf_chk_F_output
  PUBLIC :: serialize_vdf_chk_G_output
  PUBLIC :: serialize_vdf_chk_H_output
  ! Serialization for vdiff_down routine
  PUBLIC :: serialize_vdf_vd_input
  PUBLIC :: serialize_vdf_vd_output
  ! Serialization for update_surface routine
  PUBLIC :: serialize_vdf_us_input
  PUBLIC :: serialize_vdf_us_output
  ! Serialization for vdiff_up routine
  PUBLIC :: serialize_vdf_vu_input
  PUBLIC :: serialize_vdf_vu_output
  ! Serialization for nsurf_diag routine
  PUBLIC :: serialize_vdf_nd_input
  PUBLIC :: serialize_vdf_nd_output

  CONTAINS

  SUBROUTINE serialize_vdf_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile ) IF( ASSOCIATED(field%frac_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc))
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta))
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old))
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old))
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua))
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va))
    !$ser verbatim   !$ACC UPDATE HOST( field%ocu ) IF( ASSOCIATED(field%ocu))
    !$ser verbatim   !$ACC UPDATE HOST( field%ocv ) IF( ASSOCIATED(field%ocv))
    !$ser verbatim   !$ACC UPDATE HOST( field%zf ) IF( ASSOCIATED(field%zf))
    !$ser verbatim   !$ACC UPDATE HOST( field%zh ) IF( ASSOCIATED(field%zh))
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmax ) IF( ASSOCIATED(field%tasmax))
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmin ) IF( ASSOCIATED(field%tasmin))
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair))
    !$ser verbatim   !$ACC UPDATE HOST( field%mref ) IF( ASSOCIATED(field%mref))
    !$ser verbatim   !$ACC UPDATE HOST( field%geom ) IF( ASSOCIATED(field%geom))
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile ) IF( ASSOCIATED(field%z0m_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx_tile ) IF( ASSOCIATED(field%lhflx_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx_tile ) IF( ASSOCIATED(field%shflx_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%lsmask ) IF( ASSOCIATED(field%lsmask))
    !$ser verbatim   !$ACC UPDATE HOST( field%alake ) IF( ASSOCIATED(field%alake))
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfl ) IF( ASSOCIATED(field%rsfl))
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfc ) IF( ASSOCIATED(field%rsfc))
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfl ) IF( ASSOCIATED(field%ssfl))
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfc ) IF( ASSOCIATED(field%ssfc))
    !$ser verbatim   !$ACC UPDATE HOST( field%rlds ) IF( ASSOCIATED(field%rlds))
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus ) IF( ASSOCIATED(field%rlus))
    !$ser verbatim   !$ACC UPDATE HOST( field%rsds ) IF( ASSOCIATED(field%rsds))
    !$ser verbatim   !$ACC UPDATE HOST( field%rsus ) IF( ASSOCIATED(field%rsus))
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dir ) IF( ASSOCIATED(field%rvds_dir))
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dir ) IF( ASSOCIATED(field%rpds_dir))
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dir ) IF( ASSOCIATED(field%rnds_dir))
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dif ) IF( ASSOCIATED(field%rvds_dif))
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dif ) IF( ASSOCIATED(field%rpds_dif))
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dif ) IF( ASSOCIATED(field%rnds_dif))
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0))
    !$ser verbatim   !$ACC UPDATE HOST( field%csat ) IF( ASSOCIATED(field%csat))
    !$ser verbatim   !$ACC UPDATE HOST( field%cair ) IF( ASSOCIATED(field%cair))
    !$ser verbatim   !$ACC UPDATE HOST( field%z0h_lnd ) IF( ASSOCIATED(field%z0h_lnd))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir ) IF( ASSOCIATED(field%albvisdir))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir ) IF( ASSOCIATED(field%albnirdir))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif ) IF( ASSOCIATED(field%albvisdif))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif ) IF( ASSOCIATED(field%albnirdif))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_tile ) IF( ASSOCIATED(field%albvisdir_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_tile ) IF( ASSOCIATED(field%albnirdir_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_tile ) IF( ASSOCIATED(field%albvisdif_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_tile ) IF( ASSOCIATED(field%albnirdif_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo ) IF( ASSOCIATED(field%albedo))
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo_tile ) IF( ASSOCIATED(field%albedo_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%co2_flux_tile ) IF( ASSOCIATED(field%co2_flux_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%swflxsfc_tile ) IF( ASSOCIATED(field%swflxsfc_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%lwflxsfc_tile ) IF( ASSOCIATED(field%lwflxsfc_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%Tsurf ) IF( ASSOCIATED(field%Tsurf))
    !$ser verbatim   !$ACC UPDATE HOST( field%T1 ) IF( ASSOCIATED(field%T1))
    !$ser verbatim   !$ACC UPDATE HOST( field%T2 ) IF( ASSOCIATED(field%T2))
    !$ser verbatim   !$ACC UPDATE HOST( field%hi ) IF( ASSOCIATED(field%hi))
    !$ser verbatim   !$ACC UPDATE HOST( field%hs ) IF( ASSOCIATED(field%hs))
    !$ser verbatim   !$ACC UPDATE HOST( field%Qtop ) IF( ASSOCIATED(field%Qtop))
    !$ser verbatim   !$ACC UPDATE HOST( field%Qbot ) IF( ASSOCIATED(field%Qbot))
    !$ser verbatim   !$ACC UPDATE HOST( field%conc ) IF( ASSOCIATED(field%conc))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_ice ) IF( ASSOCIATED(field%albvisdir_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_ice ) IF( ASSOCIATED(field%albnirdir_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_ice ) IF( ASSOCIATED(field%albvisdif_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_ice ) IF( ASSOCIATED(field%albnirdif_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%coriol ) IF( ASSOCIATED(field%coriol))
    !$ser verbatim   !$ACC UPDATE HOST( field%tv ) IF( ASSOCIATED(field%tv))
    !$ser verbatim   !$ACC UPDATE HOST( field%aclc ) IF( ASSOCIATED(field%aclc))
    !$ser verbatim   !$ACC UPDATE HOST( field%tottem1 ) IF( ASSOCIATED(field%tottem1))
    !$ser verbatim   !$ACC UPDATE HOST( field%ustar ) IF( ASSOCIATED(field%ustar))
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar_tile ) IF( ASSOCIATED(field%wstar_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%thvsig ) IF( ASSOCIATED(field%thvsig))
    !$ser verbatim   !$ACC UPDATE HOST( field%rld_rt ) IF( ASSOCIATED(field%rld_rt))
    !$ser verbatim   !$ACC UPDATE HOST( field%rlu_rt ) IF( ASSOCIATED(field%rlu_rt))
    !$ser verbatim   !$ACC UPDATE HOST( field%totte ) IF( ASSOCIATED(field%totte))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev))
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_vdf ) IF( ASSOCIATED(tend%ua_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_vdf ) IF( ASSOCIATED(tend%va_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_vdf ) IF( ASSOCIATED(tend%qtrc_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy))
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy))
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy))
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_vdf_frac_tile=field%frac_tile(:,jb,:) IF (ASSOCIATED(field%frac_tile))
    !$ser data echam_vdf_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_vdf_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_vdf_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_vdf_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_vdf_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_vdf_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_vdf_ocu=field%ocu(:,jb) IF (ASSOCIATED(field%ocu))
    !$ser data echam_vdf_ocv=field%ocv(:,jb) IF (ASSOCIATED(field%ocv))
    !$ser data echam_vdf_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_vdf_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_vdf_tasmax=field%tasmax(:,jb) IF (ASSOCIATED(field%tasmax))
    !$ser data echam_vdf_tasmin=field%tasmin(:,jb) IF (ASSOCIATED(field%tasmin))
    !$ser data echam_vdf_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser data echam_vdf_mref=field%mref(:,:,jb) IF (ASSOCIATED(field%mref))
    !$ser data echam_vdf_geom=field%geom(:,:,jb) IF (ASSOCIATED(field%geom))
    !$ser data echam_vdf_z0m_tile=field%z0m_tile(:,jb,:) IF (ASSOCIATED(field%z0m_tile))
    !$ser data echam_vdf_ts_tile=field%ts_tile(:,jb,:) IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_vdf_lhflx_tile=field%lhflx_tile(:,jb,:) IF (ASSOCIATED(field%lhflx_tile))
    !$ser data echam_vdf_shflx_tile=field%shflx_tile(:,jb,:) IF (ASSOCIATED(field%shflx_tile))
    !$ser data echam_vdf_lsmask=field%lsmask(:,jb) IF (ASSOCIATED(field%lsmask))
    !$ser data echam_vdf_alake=field%alake(:,jb) IF (ASSOCIATED(field%alake))
    !$ser data echam_vdf_rsfl=field%rsfl(:,jb) IF (ASSOCIATED(field%rsfl))
    !$ser data echam_vdf_rsfc=field%rsfc(:,jb) IF (ASSOCIATED(field%rsfc))
    !$ser data echam_vdf_ssfl=field%ssfl(:,jb) IF (ASSOCIATED(field%ssfl))
    !$ser data echam_vdf_ssfc=field%ssfc(:,jb) IF (ASSOCIATED(field%ssfc))
    !$ser data echam_vdf_rlds=field%rlds(:,jb) IF (ASSOCIATED(field%rlds))
    !$ser data echam_vdf_rlus=field%rlus(:,jb) IF (ASSOCIATED(field%rlus))
    !$ser data echam_vdf_rsds=field%rsds(:,jb) IF (ASSOCIATED(field%rsds))
    !$ser data echam_vdf_rsus=field%rsus(:,jb) IF (ASSOCIATED(field%rsus))
    !$ser data echam_vdf_rvds_dir=field%rvds_dir(:,jb) IF (ASSOCIATED(field%rvds_dir))
    !$ser data echam_vdf_rpds_dir=field%rpds_dir(:,jb) IF (ASSOCIATED(field%rpds_dir))
    !$ser data echam_vdf_rnds_dir=field%rnds_dir(:,jb) IF (ASSOCIATED(field%rnds_dir))
    !$ser data echam_vdf_rvds_dif=field%rvds_dif(:,jb) IF (ASSOCIATED(field%rvds_dif))
    !$ser data echam_vdf_rpds_dif=field%rpds_dif(:,jb) IF (ASSOCIATED(field%rpds_dif))
    !$ser data echam_vdf_rnds_dif=field%rnds_dif(:,jb) IF (ASSOCIATED(field%rnds_dif))
    !$ser data echam_vdf_cosmu0=field%cosmu0(:,jb) IF (ASSOCIATED(field%cosmu0))
    !$ser data echam_vdf_csat=field%csat(:,jb) IF (ASSOCIATED(field%csat))
    !$ser data echam_vdf_cair=field%cair(:,jb) IF (ASSOCIATED(field%cair))
    !$ser data echam_vdf_z0h_lnd=field%z0h_lnd(:,jb) IF (ASSOCIATED(field%z0h_lnd))
    !$ser data echam_vdf_albvisdir=field%albvisdir(:,jb) IF (ASSOCIATED(field%albvisdir))
    !$ser data echam_vdf_albnirdir=field%albnirdir(:,jb) IF (ASSOCIATED(field%albnirdir))
    !$ser data echam_vdf_albvisdif=field%albvisdif(:,jb) IF (ASSOCIATED(field%albvisdif))
    !$ser data echam_vdf_albnirdif=field%albnirdif(:,jb) IF (ASSOCIATED(field%albnirdif))
    !$ser data echam_vdf_albvisdir_tile=field%albvisdir_tile(:,jb,:) IF (ASSOCIATED(field%albvisdir_tile))
    !$ser data echam_vdf_albnirdir_tile=field%albnirdir_tile(:,jb,:) IF (ASSOCIATED(field%albnirdir_tile))
    !$ser data echam_vdf_albvisdif_tile=field%albvisdif_tile(:,jb,:) IF (ASSOCIATED(field%albvisdif_tile))
    !$ser data echam_vdf_albnirdif_tile=field%albnirdif_tile(:,jb,:) IF (ASSOCIATED(field%albnirdif_tile))
    !$ser data echam_vdf_albedo=field%albedo(:,jb) IF (ASSOCIATED(field%albedo))
    !$ser data echam_vdf_albedo_tile=field%albedo_tile(:,jb,:) IF (ASSOCIATED(field%albedo_tile))
    !$ser data echam_vdf_co2_flux_tile=field%co2_flux_tile(:,jb,:) IF (ASSOCIATED(field%co2_flux_tile))
    !$ser data echam_vdf_swflxsfc_tile=field%swflxsfc_tile(:,jb,:) IF (ASSOCIATED(field%swflxsfc_tile))
    !$ser data echam_vdf_lwflxsfc_tile=field%lwflxsfc_tile(:,jb,:) IF (ASSOCIATED(field%lwflxsfc_tile))
    !$ser data echam_vdf_Tsurf=field%Tsurf(:,:,jb) IF (ASSOCIATED(field%Tsurf))
    !$ser data echam_vdf_T1=field%T1(:,:,jb) IF (ASSOCIATED(field%T1))
    !$ser data echam_vdf_T2=field%T2(:,:,jb) IF (ASSOCIATED(field%T2))
    !$ser data echam_vdf_hi=field%hi(:,:,jb) IF (ASSOCIATED(field%hi))
    !$ser data echam_vdf_hs=field%hs(:,:,jb) IF (ASSOCIATED(field%hs))
    !$ser data echam_vdf_Qtop=field%Qtop(:,:,jb) IF (ASSOCIATED(field%Qtop))
    !$ser data echam_vdf_Qbot=field%Qbot(:,:,jb) IF (ASSOCIATED(field%Qbot))
    !$ser data echam_vdf_conc=field%conc(:,:,jb) IF (ASSOCIATED(field%conc))
    !$ser data echam_vdf_albvisdir_ice=field%albvisdir_ice(:,:,jb) IF (ASSOCIATED(field%albvisdir_ice))
    !$ser data echam_vdf_albnirdir_ice=field%albnirdir_ice(:,:,jb) IF (ASSOCIATED(field%albnirdir_ice))
    !$ser data echam_vdf_albvisdif_ice=field%albvisdif_ice(:,:,jb) IF (ASSOCIATED(field%albvisdif_ice))
    !$ser data echam_vdf_albnirdif_ice=field%albnirdif_ice(:,:,jb) IF (ASSOCIATED(field%albnirdif_ice))
    !$ser data echam_vdf_coriol=field%coriol(:,jb) IF (ASSOCIATED(field%coriol))
    !$ser data echam_vdf_tv=field%tv(:,:,jb) IF (ASSOCIATED(field%tv))
    !$ser data echam_vdf_aclc=field%aclc(:,:,jb) IF (ASSOCIATED(field%aclc))
    !$ser data echam_vdf_tottem1=field%tottem1(:,:,jb) IF (ASSOCIATED(field%tottem1))
    !$ser data echam_vdf_ustar=field%ustar(:,jb) IF (ASSOCIATED(field%ustar))
    !$ser data echam_vdf_wstar_tile=field%wstar_tile(:,jb,:) IF (ASSOCIATED(field%wstar_tile))
    !$ser data echam_vdf_thvsig=field%thvsig(:,jb) IF (ASSOCIATED(field%thvsig))
    !$ser data echam_vdf_rld_rt=field%rld_rt(:,:,jb) IF (ASSOCIATED(field%rld_rt))
    !$ser data echam_vdf_rlu_rt=field%rlu_rt(:,:,jb) IF (ASSOCIATED(field%rlu_rt))
    !$ser data echam_vdf_totte=field%totte(:,:,jb) IF (ASSOCIATED(field%totte))
    !$ser data echam_vdf_q_rlw_nlev=field%q_rlw_nlev(:,jb) IF (ASSOCIATED(field%q_rlw_nlev))
    !$ser data echam_vdf_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_vdf_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_vdf_q_phy=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_vdf_ua_vdf=tend%ua_vdf(:,:,jb) IF (ASSOCIATED(tend%ua_vdf))
    !$ser data echam_vdf_va_vdf=tend%va_vdf(:,:,jb) IF (ASSOCIATED(tend%va_vdf))
    !$ser data echam_vdf_qtrc_vdf=tend%qtrc_vdf(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_vdf))
    !$ser data echam_vdf_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_vdf_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_vdf_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_vdf_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_vdf:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%frac_tile ) IF( ASSOCIATED(field%frac_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_old ) IF( ASSOCIATED(field%presi_old))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ua ) IF( ASSOCIATED(field%ua))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%va ) IF( ASSOCIATED(field%va))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ocu ) IF( ASSOCIATED(field%ocu))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ocv ) IF( ASSOCIATED(field%ocv))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zf ) IF( ASSOCIATED(field%zf))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zh ) IF( ASSOCIATED(field%zh))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tasmax ) IF( ASSOCIATED(field%tasmax))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tasmin ) IF( ASSOCIATED(field%tasmin))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mair ) IF( ASSOCIATED(field%mair))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mref ) IF( ASSOCIATED(field%mref))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%geom ) IF( ASSOCIATED(field%geom))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%z0m_tile ) IF( ASSOCIATED(field%z0m_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ts_tile ) IF( ASSOCIATED(field%ts_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%lhflx_tile ) IF( ASSOCIATED(field%lhflx_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%shflx_tile ) IF( ASSOCIATED(field%shflx_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%lsmask ) IF( ASSOCIATED(field%lsmask))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%alake ) IF( ASSOCIATED(field%alake))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsfl ) IF( ASSOCIATED(field%rsfl))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsfc ) IF( ASSOCIATED(field%rsfc))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ssfl ) IF( ASSOCIATED(field%ssfl))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ssfc ) IF( ASSOCIATED(field%ssfc))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlds ) IF( ASSOCIATED(field%rlds))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlus ) IF( ASSOCIATED(field%rlus))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsds ) IF( ASSOCIATED(field%rsds))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsus ) IF( ASSOCIATED(field%rsus))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvds_dir ) IF( ASSOCIATED(field%rvds_dir))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpds_dir ) IF( ASSOCIATED(field%rpds_dir))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnds_dir ) IF( ASSOCIATED(field%rnds_dir))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvds_dif ) IF( ASSOCIATED(field%rvds_dif))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpds_dif ) IF( ASSOCIATED(field%rpds_dif))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnds_dif ) IF( ASSOCIATED(field%rnds_dif))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%csat ) IF( ASSOCIATED(field%csat))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%cair ) IF( ASSOCIATED(field%cair))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%z0h_lnd ) IF( ASSOCIATED(field%z0h_lnd))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albvisdir ) IF( ASSOCIATED(field%albvisdir))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albnirdir ) IF( ASSOCIATED(field%albnirdir))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albvisdif ) IF( ASSOCIATED(field%albvisdif))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albnirdif ) IF( ASSOCIATED(field%albnirdif))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albvisdir_tile ) IF( ASSOCIATED(field%albvisdir_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albnirdir_tile ) IF( ASSOCIATED(field%albnirdir_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albvisdif_tile ) IF( ASSOCIATED(field%albvisdif_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albnirdif_tile ) IF( ASSOCIATED(field%albnirdif_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albedo ) IF( ASSOCIATED(field%albedo))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albedo_tile ) IF( ASSOCIATED(field%albedo_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%co2_flux_tile ) IF( ASSOCIATED(field%co2_flux_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%swflxsfc_tile ) IF( ASSOCIATED(field%swflxsfc_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%lwflxsfc_tile ) IF( ASSOCIATED(field%lwflxsfc_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%Tsurf ) IF( ASSOCIATED(field%Tsurf))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%T1 ) IF( ASSOCIATED(field%T1))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%T2 ) IF( ASSOCIATED(field%T2))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%hi ) IF( ASSOCIATED(field%hi))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%hs ) IF( ASSOCIATED(field%hs))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%Qtop ) IF( ASSOCIATED(field%Qtop))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%Qbot ) IF( ASSOCIATED(field%Qbot))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%conc ) IF( ASSOCIATED(field%conc))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albvisdir_ice ) IF( ASSOCIATED(field%albvisdir_ice))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albnirdir_ice ) IF( ASSOCIATED(field%albnirdir_ice))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albvisdif_ice ) IF( ASSOCIATED(field%albvisdif_ice))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%albnirdif_ice ) IF( ASSOCIATED(field%albnirdif_ice))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%coriol ) IF( ASSOCIATED(field%coriol))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tv ) IF( ASSOCIATED(field%tv))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%aclc ) IF( ASSOCIATED(field%aclc))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tottem1 ) IF( ASSOCIATED(field%tottem1))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ustar ) IF( ASSOCIATED(field%ustar))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%wstar_tile ) IF( ASSOCIATED(field%wstar_tile))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%thvsig ) IF( ASSOCIATED(field%thvsig))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rld_rt ) IF( ASSOCIATED(field%rld_rt))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlu_rt ) IF( ASSOCIATED(field%rlu_rt))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%totte ) IF( ASSOCIATED(field%totte))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy))
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_vdf ) IF( ASSOCIATED(tend%ua_vdf))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_vdf ) IF( ASSOCIATED(tend%va_vdf))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_vdf ) IF( ASSOCIATED(tend%qtrc_vdf))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_phy ) IF( ASSOCIATED(tend%va_phy))
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy))
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_input

  SUBROUTINE serialize_vdf_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%sfcWind ) IF( ASSOCIATED(field%sfcWind))
    !$ser verbatim   !$ACC UPDATE HOST( field%tas ) IF( ASSOCIATED(field%tas))
    !$ser verbatim   !$ACC UPDATE HOST( field%dew2 ) IF( ASSOCIATED(field%dew2))
    !$ser verbatim   !$ACC UPDATE HOST( field%uas ) IF( ASSOCIATED(field%uas))
    !$ser verbatim   !$ACC UPDATE HOST( field%vas ) IF( ASSOCIATED(field%vas))
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmax ) IF( ASSOCIATED(field%tasmax))
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmin ) IF( ASSOCIATED(field%tasmin))
    !$ser verbatim   !$ACC UPDATE HOST( field%sfcWind_tile ) IF( ASSOCIATED(field%sfcWind_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%tas_tile ) IF( ASSOCIATED(field%tas_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%dew2_tile ) IF( ASSOCIATED(field%dew2_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%uas_tile ) IF( ASSOCIATED(field%uas_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%vas_tile ) IF( ASSOCIATED(field%vas_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile ) IF( ASSOCIATED(field%z0m_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m ) IF( ASSOCIATED(field%z0m))
    !$ser verbatim   !$ACC UPDATE HOST( field%totte ) IF( ASSOCIATED(field%totte))
    !$ser verbatim   !$ACC UPDATE HOST( field%sh_vdiff ) IF( ASSOCIATED(field%sh_vdiff))
    !$ser verbatim   !$ACC UPDATE HOST( field%qv_vdiff ) IF( ASSOCIATED(field%qv_vdiff))
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress ) IF( ASSOCIATED(field%u_stress))
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress ) IF( ASSOCIATED(field%v_stress))
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx ) IF( ASSOCIATED(field%lhflx))
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx ) IF( ASSOCIATED(field%shflx))
    !$ser verbatim   !$ACC UPDATE HOST( field%evap ) IF( ASSOCIATED(field%evap))
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress_tile ) IF( ASSOCIATED(field%u_stress_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress_tile ) IF( ASSOCIATED(field%v_stress_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx_tile ) IF( ASSOCIATED(field%lhflx_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx_tile ) IF( ASSOCIATED(field%shflx_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%evap_tile ) IF( ASSOCIATED(field%evap_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%fco2nat ) IF( ASSOCIATED(field%fco2nat))
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus ) IF( ASSOCIATED(field%rlus))
    !$ser verbatim   !$ACC UPDATE HOST( field%csat ) IF( ASSOCIATED(field%csat))
    !$ser verbatim   !$ACC UPDATE HOST( field%cair ) IF( ASSOCIATED(field%cair))
    !$ser verbatim   !$ACC UPDATE HOST( field%z0h_lnd ) IF( ASSOCIATED(field%z0h_lnd))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir ) IF( ASSOCIATED(field%albvisdir))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir ) IF( ASSOCIATED(field%albnirdir))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif ) IF( ASSOCIATED(field%albvisdif))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif ) IF( ASSOCIATED(field%albnirdif))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_tile ) IF( ASSOCIATED(field%albvisdir_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_tile ) IF( ASSOCIATED(field%albnirdir_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_tile ) IF( ASSOCIATED(field%albvisdif_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_tile ) IF( ASSOCIATED(field%albnirdif_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo ) IF( ASSOCIATED(field%albedo))
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo_tile ) IF( ASSOCIATED(field%albedo_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%co2_flux_tile ) IF( ASSOCIATED(field%co2_flux_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%ts ) IF( ASSOCIATED(field%ts))
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad ) IF( ASSOCIATED(field%ts_rad))
    !$ser verbatim   !$ACC UPDATE HOST( field%swflxsfc_tile ) IF( ASSOCIATED(field%swflxsfc_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%lwflxsfc_tile ) IF( ASSOCIATED(field%lwflxsfc_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%lake_ice_frc ) IF( ASSOCIATED(field%lake_ice_frc))
    !$ser verbatim   !$ACC UPDATE HOST( field%Tsurf ) IF( ASSOCIATED(field%Tsurf))
    !$ser verbatim   !$ACC UPDATE HOST( field%T1 ) IF( ASSOCIATED(field%T1))
    !$ser verbatim   !$ACC UPDATE HOST( field%T2 ) IF( ASSOCIATED(field%T2))
    !$ser verbatim   !$ACC UPDATE HOST( field%hi ) IF( ASSOCIATED(field%hi))
    !$ser verbatim   !$ACC UPDATE HOST( field%hs ) IF( ASSOCIATED(field%hs))
    !$ser verbatim   !$ACC UPDATE HOST( field%Qtop ) IF( ASSOCIATED(field%Qtop))
    !$ser verbatim   !$ACC UPDATE HOST( field%Qbot ) IF( ASSOCIATED(field%Qbot))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_ice ) IF( ASSOCIATED(field%albvisdir_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_ice ) IF( ASSOCIATED(field%albnirdir_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_ice ) IF( ASSOCIATED(field%albvisdif_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_ice ) IF( ASSOCIATED(field%albnirdif_ice))
    !$ser verbatim   !$ACC UPDATE HOST( field%ustar ) IF( ASSOCIATED(field%ustar))
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar_tile ) IF( ASSOCIATED(field%wstar_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%thvsig ) IF( ASSOCIATED(field%thvsig))
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar ) IF( ASSOCIATED(field%wstar))
    !$ser verbatim   !$ACC UPDATE HOST( field%qs_sfc_tile ) IF( ASSOCIATED(field%qs_sfc_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%hdtcbl ) IF( ASSOCIATED(field%hdtcbl))
    !$ser verbatim   !$ACC UPDATE HOST( field%ri_atm ) IF( ASSOCIATED(field%ri_atm))
    !$ser verbatim   !$ACC UPDATE HOST( field%mixlen ) IF( ASSOCIATED(field%mixlen))
    !$ser verbatim   !$ACC UPDATE HOST( field%cfm ) IF( ASSOCIATED(field%cfm))
    !$ser verbatim   !$ACC UPDATE HOST( field%cfm_tile ) IF( ASSOCIATED(field%cfm_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%cfh ) IF( ASSOCIATED(field%cfh))
    !$ser verbatim   !$ACC UPDATE HOST( field%cfh_tile ) IF( ASSOCIATED(field%cfh_tile))
    !$ser verbatim   !$ACC UPDATE HOST( field%cfv ) IF( ASSOCIATED(field%cfv))
    !$ser verbatim   !$ACC UPDATE HOST( field%cftotte ) IF( ASSOCIATED(field%cftotte))
    !$ser verbatim   !$ACC UPDATE HOST( field%cfthv ) IF( ASSOCIATED(field%cfthv))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_snocpymlt ) IF( ASSOCIATED(field%q_snocpymlt))
    !$ser verbatim   !$ACC UPDATE HOST( field%kedisp ) IF( ASSOCIATED(field%kedisp))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_vdf ) IF( ASSOCIATED(field%q_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_vdf_vi ) IF( ASSOCIATED(field%q_vdf_vi))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_vdf ) IF( ASSOCIATED(tend%ua_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_vdf ) IF( ASSOCIATED(tend%va_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_vdf ) IF( ASSOCIATED(tend%qtrc_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_vdf ) IF( ASSOCIATED(tend%ta_vdf))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy))
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy))
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy))
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_impl ) IF( ASSOCIATED(field%q_rlw_impl))
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rlw_impl ) IF( ASSOCIATED(tend%ta_rlw_impl))
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_sfcWind=field%sfcWind(:,jb) IF (ASSOCIATED(field%sfcWind))
    !$ser data echam_vdf_tas=field%tas(:,jb) IF (ASSOCIATED(field%tas))
    !$ser data echam_vdf_dew2=field%dew2(:,jb) IF (ASSOCIATED(field%dew2))
    !$ser data echam_vdf_uas=field%uas(:,jb) IF (ASSOCIATED(field%uas))
    !$ser data echam_vdf_vas=field%vas(:,jb) IF (ASSOCIATED(field%vas))
    !$ser data echam_vdf_tasmax=field%tasmax(:,jb) IF (ASSOCIATED(field%tasmax))
    !$ser data echam_vdf_tasmin=field%tasmin(:,jb) IF (ASSOCIATED(field%tasmin))
    !$ser data echam_vdf_sfcWind_tile=field%sfcWind_tile(:,jb,:) IF (ASSOCIATED(field%sfcWind_tile))
    !$ser data echam_vdf_tas_tile=field%tas_tile(:,jb,:) IF (ASSOCIATED(field%tas_tile))
    !$ser data echam_vdf_dew2_tile=field%dew2_tile(:,jb,:) IF (ASSOCIATED(field%dew2_tile))
    !$ser data echam_vdf_uas_tile=field%uas_tile(:,jb,:) IF (ASSOCIATED(field%uas_tile))
    !$ser data echam_vdf_vas_tile=field%vas_tile(:,jb,:) IF (ASSOCIATED(field%vas_tile))
    !$ser data echam_vdf_z0m_tile=field%z0m_tile(:,jb,:) IF (ASSOCIATED(field%z0m_tile))
    !$ser data echam_vdf_z0m=field%z0m(:,jb) IF (ASSOCIATED(field%z0m))
    !$ser data echam_vdf_totte=field%totte(:,:,jb) IF (ASSOCIATED(field%totte))
    !$ser data echam_vdf_sh_vdiff=field%sh_vdiff(:,jb) IF (ASSOCIATED(field%sh_vdiff))
    !$ser data echam_vdf_qv_vdiff=field%qv_vdiff(:,jb) IF (ASSOCIATED(field%qv_vdiff))
    !$ser data echam_vdf_ts_tile=field%ts_tile(:,jb,:) IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_vdf_u_stress=field%u_stress(:,jb) IF (ASSOCIATED(field%u_stress))
    !$ser data echam_vdf_v_stress=field%v_stress(:,jb) IF (ASSOCIATED(field%v_stress))
    !$ser data echam_vdf_lhflx=field%lhflx(:,jb) IF (ASSOCIATED(field%lhflx))
    !$ser data echam_vdf_shflx=field%shflx(:,jb) IF (ASSOCIATED(field%shflx))
    !$ser data echam_vdf_evap=field%evap(:,jb) IF (ASSOCIATED(field%evap))
    !$ser data echam_vdf_u_stress_tile=field%u_stress_tile(:,jb,:) IF (ASSOCIATED(field%u_stress_tile))
    !$ser data echam_vdf_v_stress_tile=field%v_stress_tile(:,jb,:) IF (ASSOCIATED(field%v_stress_tile))
    !$ser data echam_vdf_lhflx_tile=field%lhflx_tile(:,jb,:) IF (ASSOCIATED(field%lhflx_tile))
    !$ser data echam_vdf_shflx_tile=field%shflx_tile(:,jb,:) IF (ASSOCIATED(field%shflx_tile))
    !$ser data echam_vdf_evap_tile=field%evap_tile(:,jb,:) IF (ASSOCIATED(field%evap_tile))
    !$ser data echam_vdf_fco2nat=field%fco2nat(:,jb) IF (ASSOCIATED(field%fco2nat))
    !$ser data echam_vdf_rlus=field%rlus(:,jb) IF (ASSOCIATED(field%rlus))
    !$ser data echam_vdf_csat=field%csat(:,jb) IF (ASSOCIATED(field%csat))
    !$ser data echam_vdf_cair=field%cair(:,jb) IF (ASSOCIATED(field%cair))
    !$ser data echam_vdf_z0h_lnd=field%z0h_lnd(:,jb) IF (ASSOCIATED(field%z0h_lnd))
    !$ser data echam_vdf_albvisdir=field%albvisdir(:,jb) IF (ASSOCIATED(field%albvisdir))
    !$ser data echam_vdf_albnirdir=field%albnirdir(:,jb) IF (ASSOCIATED(field%albnirdir))
    !$ser data echam_vdf_albvisdif=field%albvisdif(:,jb) IF (ASSOCIATED(field%albvisdif))
    !$ser data echam_vdf_albnirdif=field%albnirdif(:,jb) IF (ASSOCIATED(field%albnirdif))
    !$ser data echam_vdf_albvisdir_tile=field%albvisdir_tile(:,jb,:) IF (ASSOCIATED(field%albvisdir_tile))
    !$ser data echam_vdf_albnirdir_tile=field%albnirdir_tile(:,jb,:) IF (ASSOCIATED(field%albnirdir_tile))
    !$ser data echam_vdf_albvisdif_tile=field%albvisdif_tile(:,jb,:) IF (ASSOCIATED(field%albvisdif_tile))
    !$ser data echam_vdf_albnirdif_tile=field%albnirdif_tile(:,jb,:) IF (ASSOCIATED(field%albnirdif_tile))
    !$ser data echam_vdf_albedo=field%albedo(:,jb) IF (ASSOCIATED(field%albedo))
    !$ser data echam_vdf_albedo_tile=field%albedo_tile(:,jb,:) IF (ASSOCIATED(field%albedo_tile))
    !$ser data echam_vdf_co2_flux_tile=field%co2_flux_tile(:,jb,:) IF (ASSOCIATED(field%co2_flux_tile))
    !$ser data echam_vdf_ts=field%ts(:,jb) IF (ASSOCIATED(field%ts))
    !$ser data echam_vdf_ts_rad=field%ts_rad(:,jb) IF (ASSOCIATED(field%ts_rad))
    !$ser data echam_vdf_swflxsfc_tile=field%swflxsfc_tile(:,jb,:) IF (ASSOCIATED(field%swflxsfc_tile))
    !$ser data echam_vdf_lwflxsfc_tile=field%lwflxsfc_tile(:,jb,:) IF (ASSOCIATED(field%lwflxsfc_tile))
    !$ser data echam_vdf_lake_ice_frc=field%lake_ice_frc(:,jb) IF (ASSOCIATED(field%lake_ice_frc))
    !$ser data echam_vdf_Tsurf=field%Tsurf(:,:,jb) IF (ASSOCIATED(field%Tsurf))
    !$ser data echam_vdf_T1=field%T1(:,:,jb) IF (ASSOCIATED(field%T1))
    !$ser data echam_vdf_T2=field%T2(:,:,jb) IF (ASSOCIATED(field%T2))
    !$ser data echam_vdf_hi=field%hi(:,:,jb) IF (ASSOCIATED(field%hi))
    !$ser data echam_vdf_hs=field%hs(:,:,jb) IF (ASSOCIATED(field%hs))
    !$ser data echam_vdf_Qtop=field%Qtop(:,:,jb) IF (ASSOCIATED(field%Qtop))
    !$ser data echam_vdf_Qbot=field%Qbot(:,:,jb) IF (ASSOCIATED(field%Qbot))
    !$ser data echam_vdf_albvisdir_ice=field%albvisdir_ice(:,:,jb) IF (ASSOCIATED(field%albvisdir_ice))
    !$ser data echam_vdf_albnirdir_ice=field%albnirdir_ice(:,:,jb) IF (ASSOCIATED(field%albnirdir_ice))
    !$ser data echam_vdf_albvisdif_ice=field%albvisdif_ice(:,:,jb) IF (ASSOCIATED(field%albvisdif_ice))
    !$ser data echam_vdf_albnirdif_ice=field%albnirdif_ice(:,:,jb) IF (ASSOCIATED(field%albnirdif_ice))
    !$ser data echam_vdf_ustar=field%ustar(:,jb) IF (ASSOCIATED(field%ustar))
    !$ser data echam_vdf_wstar_tile=field%wstar_tile(:,jb,:) IF (ASSOCIATED(field%wstar_tile))
    !$ser data echam_vdf_thvsig=field%thvsig(:,jb) IF (ASSOCIATED(field%thvsig))
    !$ser data echam_vdf_wstar=field%wstar(:,jb) IF (ASSOCIATED(field%wstar))
    !$ser data echam_vdf_qs_sfc_tile=field%qs_sfc_tile(:,jb,:) IF (ASSOCIATED(field%qs_sfc_tile))
    !$ser data echam_vdf_hdtcbl=field%hdtcbl(:,jb) IF (ASSOCIATED(field%hdtcbl))
    !$ser data echam_vdf_ri_atm=field%ri_atm(:,:,jb) IF (ASSOCIATED(field%ri_atm))
    !$ser data echam_vdf_mixlen=field%mixlen(:,:,jb) IF (ASSOCIATED(field%mixlen))
    !$ser data echam_vdf_cfm=field%cfm(:,:,jb) IF (ASSOCIATED(field%cfm))
    !$ser data echam_vdf_cfm_tile=field%cfm_tile(:,jb,:) IF (ASSOCIATED(field%cfm_tile))
    !$ser data echam_vdf_cfh=field%cfh(:,:,jb) IF (ASSOCIATED(field%cfh))
    !$ser data echam_vdf_cfh_tile=field%cfh_tile(:,jb,:) IF (ASSOCIATED(field%cfh_tile))
    !$ser data echam_vdf_cfv=field%cfv(:,:,jb) IF (ASSOCIATED(field%cfv))
    !$ser data echam_vdf_cftotte=field%cftotte(:,:,jb) IF (ASSOCIATED(field%cftotte))
    !$ser data echam_vdf_cfthv=field%cfthv(:,:,jb) IF (ASSOCIATED(field%cfthv))
    !$ser data echam_vdf_q_snocpymlt=field%q_snocpymlt IF (ASSOCIATED(field%q_snocpymlt))
    !$ser data echam_vdf_kedisp=field%kedisp(:,jb) IF (ASSOCIATED(field%kedisp))
    !$ser data echam_vdf_q_vdf=field%q_vdf(:,:,jb) IF (ASSOCIATED(field%q_vdf))
    !$ser data echam_vdf_q_vdf_vi=field%q_vdf_vi(:,jb) IF (ASSOCIATED(field%q_vdf_vi))
    !$ser data echam_vdf_ua_vdf=tend%ua_vdf(:,:,jb) IF (ASSOCIATED(tend%ua_vdf))
    !$ser data echam_vdf_va_vdf=tend%va_vdf(:,:,jb) IF (ASSOCIATED(tend%va_vdf))
    !$ser data echam_vdf_qtrc_vdf=tend%qtrc_vdf(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_vdf))
    !$ser data echam_vdf_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_vdf_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_vdf_ta_vdf=tend%ta_vdf(:,:,jb) IF (ASSOCIATED(tend%ta_vdf))
    !$ser data echam_vdf_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_vdf_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_vdf_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_vdf_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser data echam_vdf_q_rlw_impl=field%q_rlw_impl(:,jb) IF (ASSOCIATED(field%q_rlw_impl))
    !$ser data echam_vdf_ta_rlw_impl=tend%ta_rlw_impl(:,jb) IF (ASSOCIATED(tend%ta_rlw_impl))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_output

  SUBROUTINE serialize_vdf_chk_A_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field, zxt_emis, zco2, dummy, dummyx, zqx)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
        zxt_emis(:,:),        &
        zco2(:),              &
        dummy(:,:),           &
        dummyx(:,:),          &
        zqx(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_A:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_A:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( zxt_emis ) IF( SIZE(zxt_emis) > 0 )
    !$ser verbatim   !$ACC UPDATE HOST( zco2 )
    !$ser verbatim   !$ACC UPDATE HOST( dummy )
    !$ser verbatim   !$ACC UPDATE HOST( dummyx )
    !$ser verbatim   !$ACC UPDATE HOST( zqx )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_a-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_a_zco2=zco2                              &
    !$ser&     echam_vdf_chk_a_dummy=dummy                         &
    !$ser&     echam_vdf_chk_a_dummyx=dummyx                          &
    !$ser&     echam_vdf_chk_a_zqx=zqx
    !$ser data echam_vdf_chk_a_zxt_emis=zxt_emis IF (SIZE(zxt_emis)>0)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_A_output

  SUBROUTINE serialize_vdf_chk_B_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_B:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_B:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar ) IF( ASSOCIATED(field%wstar) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qs_sfc_tile ) IF( ASSOCIATED(field%qs_sfc_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%hdtcbl ) IF( ASSOCIATED(field%hdtcbl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ri_atm ) IF( ASSOCIATED(field%ri_atm) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mixlen ) IF( ASSOCIATED(field%mixlen) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfm ) IF( ASSOCIATED(field%cfm) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfm_tile ) IF( ASSOCIATED(field%cfm_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfh ) IF( ASSOCIATED(field%cfh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfh_tile ) IF( ASSOCIATED(field%cfh_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfv ) IF( ASSOCIATED(field%cfv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cftotte ) IF( ASSOCIATED(field%cftotte) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfthv ) IF( ASSOCIATED(field%cfthv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap_tile )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_b-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_b_wstar=field%wstar IF (ASSOCIATED(field%wstar))
    !$ser data echam_vdf_chk_b_qs_sfc_tile=field%qs_sfc_tile IF (ASSOCIATED(field%qs_sfc_tile))
    !$ser data echam_vdf_chk_b_hdtcbl=field%hdtcbl IF (ASSOCIATED(field%hdtcbl))
    !$ser data echam_vdf_chk_b_ri_atm=field%ri_atm IF (ASSOCIATED(field%ri_atm))
    !$ser data echam_vdf_chk_b_mixlen=field%mixlen IF (ASSOCIATED(field%mixlen))
    !$ser data echam_vdf_chk_b_cfm=field%cfm IF (ASSOCIATED(field%cfm))
    !$ser data echam_vdf_chk_b_cfm_tile=field%cfm_tile IF (ASSOCIATED(field%cfm_tile))
    !$ser data echam_vdf_chk_b_cfh=field%cfh IF (ASSOCIATED(field%cfh))
    !$ser data echam_vdf_chk_b_cfh_tile=field%cfh_tile IF (ASSOCIATED(field%cfh_tile))
    !$ser data echam_vdf_chk_b_cfv=field%cfv IF (ASSOCIATED(field%cfv))
    !$ser data echam_vdf_chk_b_cftotte=field%cftotte IF (ASSOCIATED(field%cftotte))
    !$ser data echam_vdf_chk_b_cfthv=field%cfthv IF (ASSOCIATED(field%cfthv))
    !$ser data echam_vdf_chk_b_lhflx_tile=field%lhflx_tile             &
    !$ser&     echam_vdf_chk_b_shflx_tile=field%shflx_tile             &
    !$ser&     echam_vdf_chk_b_evap_tile=field%evap_tile
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_B_output

  SUBROUTINE serialize_vdf_chk_C_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_C:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_C:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%q_snocpymlt ) IF( ASSOCIATED(field%q_snocpymlt) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_c-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_c_q_snocpymlt=field%q_snocpymlt IF (ASSOCIATED(field%q_snocpymlt))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_C_output

  SUBROUTINE serialize_vdf_chk_D_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), INTENT(INOUT) :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_D:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_D:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%kedisp ) IF( ASSOCIATED(field%kedisp) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_vdf ) IF( ASSOCIATED(field%q_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_vdf_vi ) IF( ASSOCIATED(field%q_vdf_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_vdf ) IF( ASSOCIATED(tend%ua_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_vdf ) IF( ASSOCIATED(tend%va_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_vdf ) IF( ASSOCIATED(tend%qtrc_vdf) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_d-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_d_kedisp=field%kedisp IF (ASSOCIATED(field%kedisp))
    !$ser data echam_vdf_chk_d_q_vdf=field%q_vdf IF (ASSOCIATED(field%q_vdf))
    !$ser data echam_vdf_chk_d_q_vdf_vi=field%q_vdf_vi IF (ASSOCIATED(field%q_vdf_vi))
    !$ser data echam_vdf_chk_d_ua_vdf=tend%ua_vdf IF (ASSOCIATED(tend%ua_vdf))
    !$ser data echam_vdf_chk_d_va_vdf=tend%va_vdf IF (ASSOCIATED(tend%va_vdf))
    !$ser data echam_vdf_chk_d_qtrc_vdf=tend%qtrc_vdf IF (ASSOCIATED(tend%qtrc_vdf))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_D_output

  SUBROUTINE serialize_vdf_chk_E_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field, tend, q_snocpymlt, q_vdf, tend_ua_vdf, tend_va_vdf,&
                                tend_qtrc_vdf)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), INTENT(INOUT) :: tend
    REAL(wp), INTENT(IN)  ::  &
      q_snocpymlt(:),         &
      q_vdf(:,:),             &
      tend_ua_vdf(:,:),       &
      tend_va_vdf(:,:),       &
      tend_qtrc_vdf(:,:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_E:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_E:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( q_snocpymlt )
    !$ser verbatim   !$ACC UPDATE HOST( q_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( tend_ua_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( tend_va_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( tend_qtrc_vdf )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_e-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_e_q_snocpymlt=q_snocpymlt
    !$ser data echam_vdf_chk_e_q_vdf=q_vdf
    !$ser data echam_vdf_chk_e_tend_ua_vdf=tend_ua_vdf
    !$ser data echam_vdf_chk_e_tend_va_vdf=tend_va_vdf
    !$ser data echam_vdf_chk_e_tend_qtrc_vdf=tend_qtrc_vdf
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_E_output

  SUBROUTINE serialize_vdf_chk_F_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field, tend, tend_ta_sfc)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), INTENT(INOUT) :: tend
    REAL(wp), INTENT(IN)  ::  &
      tend_ta_sfc(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_F:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_F:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( tend_ta_sfc )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_sfc ) IF( ASSOCIATED(tend%ta_sfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_f-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_f_tend_ta_sfc=tend_ta_sfc
    !$ser data echam_vdf_chk_f_ta_sfc=tend%ta_sfc IF (ASSOCIATED(tend%ta_sfc))
    !$ser data echam_vdf_chk_f_q_phy=field%q_phy IF (ASSOCIATED(field%q_phy))
    !$ser data echam_vdf_chk_f_q_phy_vi=field%q_phy_vi IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_vdf_chk_f_ta_phy=tend%ta_phy IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_vdf_chk_f_ta=field%ta IF (ASSOCIATED(field%ta))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_F_output

  SUBROUTINE serialize_vdf_chk_G_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field, tend, tend_ta_vdf, q_rlw_impl, tend_ta_rlw_impl)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), INTENT(INOUT) :: tend
    REAL(wp), INTENT(IN)  ::  &
      tend_ta_vdf(:,:),       &
      q_rlw_impl(:),       &
      tend_ta_rlw_impl(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_G:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_G:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( tend_ta_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_vdf ) IF( ASSOCIATED(tend%ta_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( q_rlw_impl )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_impl ) IF( ASSOCIATED(field%q_rlw_impl) )
    !$ser verbatim   !$ACC UPDATE HOST( tend_ta_rlw_impl )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rlw_impl ) IF( ASSOCIATED(tend%ta_rlw_impl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%tottem1 ) IF( ASSOCIATED(field%tottem1) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_g-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_g_tend_ta_vdf=tend_ta_vdf
    !$ser data echam_vdf_chk_g_ta_vdf=tend%ta_vdf IF (ASSOCIATED(tend%ta_vdf))
    !$ser data echam_vdf_chk_g_q_phy=field%q_phy IF (ASSOCIATED(field%q_phy))
    !$ser data echam_vdf_chk_g_q_phy_vi=field%q_phy_vi IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_vdf_chk_g_ua_phy=tend%ua_phy IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_vdf_chk_g_va_phy=tend%va_phy IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_vdf_chk_g_ta_phy=tend%ta_phy IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_vdf_chk_g_qtrc_phy=tend%qtrc_phy IF (ASSOCIATED(tend%qtrc_phy))
    !$ser data echam_vdf_chk_g_ua=field%ua IF (ASSOCIATED(field%ua))
    !$ser data echam_vdf_chk_g_va=field%va IF (ASSOCIATED(field%va))
    !$ser data echam_vdf_chk_g_ta=field%ta IF (ASSOCIATED(field%ta))
    !$ser data echam_vdf_chk_g_qtrc=field%qtrc IF (ASSOCIATED(field%qtrc))
    !$ser data echam_vdf_chk_g_q_rlw_impl=q_rlw_impl
    !$ser data echam_vdf_chk_g_field_q_rlw_impl=field%q_rlw_impl IF (ASSOCIATED(field%q_rlw_impl))
    !$ser data echam_vdf_chk_g_tend_ta_rlw_impl=tend_ta_rlw_impl
    !$ser data echam_vdf_chk_g_ta_rlw_impl=tend%ta_rlw_impl IF (ASSOCIATED(tend%ta_rlw_impl))
    !$ser data echam_vdf_chk_g_tottem1=field%tottem1 IF (ASSOCIATED(field%tottem1))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_G_output

  SUBROUTINE serialize_vdf_chk_H_output(jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type, pdtime, &
                                field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev, ntrac, nsfc_type
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), INTENT(INOUT) :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_chk_H:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_chk_H:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ustar ) IF( ASSOCIATED(field%ustar) )
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar ) IF( ASSOCIATED(field%wstar) )
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar_tile ) IF( ASSOCIATED(field%wstar_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qs_sfc_tile ) IF( ASSOCIATED(field%qs_sfc_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%hdtcbl ) IF( ASSOCIATED(field%hdtcbl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ri_atm ) IF (ASSOCIATED(field%ri_atm) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mixlen ) IF( ASSOCIATED(field%mixlen) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfm ) IF( ASSOCIATED(field%cfm) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfm_tile ) IF( ASSOCIATED(field%cfm_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfh ) IF( ASSOCIATED(field%cfh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfh_tile ) IF( ASSOCIATED(field%cfh_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfv ) IF( ASSOCIATED(field%cfv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cftotte ) IF( ASSOCIATED(field%cftotte) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cfthv ) IF( ASSOCIATED(field%cfthv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%thvsig ) IF( ASSOCIATED(field%thvsig) )
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress ) IF( ASSOCIATED(field%u_stress) )
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress ) IF( ASSOCIATED(field%v_stress) )
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx ) IF( ASSOCIATED(field%lhflx) )
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx ) IF( ASSOCIATED(field%shflx) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus ) IF( ASSOCIATED(field%rlus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsus ) IF( ASSOCIATED(field%rsus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%csat ) IF( ASSOCIATED(field%csat) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cair ) IF( ASSOCIATED(field%cair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0h_lnd ) IF( ASSOCIATED(field%z0h_lnd) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir ) IF( ASSOCIATED(field%albvisdir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir ) IF( ASSOCIATED(field%albnirdir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif ) IF( ASSOCIATED(field%albvisdif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif ) IF( ASSOCIATED(field%albnirdif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo ) IF( ASSOCIATED(field%albedo) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts ) IF( ASSOCIATED(field%ts) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad ) IF( ASSOCIATED(field%ts_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress_tile ) IF( ASSOCIATED(field%u_stress_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress_tile ) IF( ASSOCIATED(field%v_stress_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx_tile ) IF( ASSOCIATED(field%lhflx_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx_tile ) IF( ASSOCIATED(field%shflx_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap_tile ) IF( ASSOCIATED(field%evap_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%lwflxsfc_tile ) IF( ASSOCIATED(field%lwflxsfc_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%swflxsfc_tile ) IF( ASSOCIATED(field%swflxsfc_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_tile ) IF( ASSOCIATED(field%albvisdir_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_tile ) IF( ASSOCIATED(field%albnirdir_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_tile ) IF( ASSOCIATED(field%albvisdif_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_tile ) IF( ASSOCIATED(field%albnirdif_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo_tile ) IF( ASSOCIATED(field%albedo_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%co2_flux_tile ) IF( ASSOCIATED(field%co2_flux_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_snocpymlt ) IF( ASSOCIATED(field%q_snocpymlt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_impl ) IF( ASSOCIATED(field%q_rlw_impl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%Tsurf ) IF( ASSOCIATED(field%Tsurf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%T1 ) IF( ASSOCIATED(field%T1) )
    !$ser verbatim   !$ACC UPDATE HOST( field%T2 ) IF( ASSOCIATED(field%T2) )
    !$ser verbatim   !$ACC UPDATE HOST( field%Qtop ) IF( ASSOCIATED(field%QTop) )
    !$ser verbatim   !$ACC UPDATE HOST( field%Qbot ) IF( ASSOCIATED(field%Qbot) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_ice ) IF( ASSOCIATED(field%albvisdir_ice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_ice ) IF( ASSOCIATED(field%albnirdir_ice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_ice ) IF( ASSOCIATED(field%albvisdif_ice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_ice ) IF( ASSOCIATED(field%albnirdif_ice) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_sfc ) IF( ASSOCIATED(tend%ta_sfc) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rlw_impl ) IF( ASSOCIATED(tend%ta_rlw_impl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%totte ) IF( ASSOCIATED(field%totte) )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m ) IF( ASSOCIATED(field%z0m) )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile ) IF( ASSOCIATED(field%z0m_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%kedisp ) IF( ASSOCIATED(field%kedisp) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sh_vdiff ) IF( ASSOCIATED(field%sh_vdiff) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qv_vdiff ) IF( ASSOCIATED(field%qv_vdiff) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_vdf ) IF( ASSOCIATED(field%q_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_vdf_vi ) IF( ASSOCIATED(field%q_vdf_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_vdf ) IF( ASSOCIATED(tend%ta_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_vdf ) IF( ASSOCIATED(tend%va_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_vdf ) IF( ASSOCIATED(tend%ua_vdf) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_vdf ) IF( ASSOCIATED(tend%qtrc_vdf) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_chk_h-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev &
    !$ser&          ntrac=ntrac nsfc_type=nsfc_type &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_chk_h_ustar=field%ustar IF (ASSOCIATED(field%ustar))
    !$ser data echam_vdf_chk_h_wstar=field%wstar IF (ASSOCIATED(field%wstar))
    !$ser data echam_vdf_chk_h_wstar_tile=field%wstar_tile IF (ASSOCIATED(field%wstar_tile))
    !$ser data echam_vdf_chk_h_qs_sfc_tile=field%qs_sfc_tile IF (ASSOCIATED(field%qs_sfc_tile))
    !$ser data echam_vdf_chk_h_hdtcbl=field%hdtcbl IF (ASSOCIATED(field%hdtcbl))
    !$ser data echam_vdf_chk_h_ri_atm=field%ri_atm IF (ASSOCIATED(field%ri_atm))
    !$ser data echam_vdf_chk_h_mixlen=field%mixlen IF (ASSOCIATED(field%mixlen))
    !$ser data echam_vdf_chk_h_cfm=field%cfm IF (ASSOCIATED(field%cfm))
    !$ser data echam_vdf_chk_h_cfm_tile=field%cfm_tile IF (ASSOCIATED(field%cfm_tile))
    !$ser data echam_vdf_chk_h_cfh=field%cfh IF (ASSOCIATED(field%cfh))
    !$ser data echam_vdf_chk_h_cfh_tile=field%cfh_tile IF (ASSOCIATED(field%cfh_tile))
    !$ser data echam_vdf_chk_h_cfv=field%cfv IF (ASSOCIATED(field%cfv))
    !$ser data echam_vdf_chk_h_cftotte=field%cftotte IF (ASSOCIATED(field%cftotte))
    !$ser data echam_vdf_chk_h_cfthv=field%cfthv IF (ASSOCIATED(field%cfthv))
    !$ser data echam_vdf_chk_h_thvsig=field%thvsig IF (ASSOCIATED(field%thvsig))
    !$ser data echam_vdf_chk_h_u_stress=field%u_stress IF (ASSOCIATED(field%u_stress))
    !$ser data echam_vdf_chk_h_v_stress=field%v_stress IF (ASSOCIATED(field%v_stress))
    !$ser data echam_vdf_chk_h_lhflx=field%lhflx IF (ASSOCIATED(field%lhflx))
    !$ser data echam_vdf_chk_h_shflx=field%shflx IF (ASSOCIATED(field%shflx))
    !$ser data echam_vdf_chk_h_evap=field%evap IF (ASSOCIATED(field%evap))
    !$ser data echam_vdf_chk_h_rlus=field%rlus IF (ASSOCIATED(field%rlus))
    !$ser data echam_vdf_chk_h_rsus=field%rsus IF (ASSOCIATED(field%rsus))
    !$ser data echam_vdf_chk_h_csat=field%csat IF (ASSOCIATED(field%csat))
    !$ser data echam_vdf_chk_h_cair=field%cair IF (ASSOCIATED(field%cair))
    !$ser data echam_vdf_chk_h_z0h_lnd=field%z0h_lnd IF (ASSOCIATED(field%z0h_lnd))
    !$ser data echam_vdf_chk_h_albvisdir=field%albvisdir IF (ASSOCIATED(field%albvisdir))
    !$ser data echam_vdf_chk_h_albnirdir=field%albnirdir IF (ASSOCIATED(field%albnirdir))
    !$ser data echam_vdf_chk_h_albvisdif=field%albvisdif IF (ASSOCIATED(field%albvisdif))
    !$ser data echam_vdf_chk_h_albnirdif=field%albnirdif IF (ASSOCIATED(field%albnirdif))
    !$ser data echam_vdf_chk_h_albedo=field%albedo IF (ASSOCIATED(field%albedo))
    !$ser data echam_vdf_chk_h_ts=field%ts IF (ASSOCIATED(field%ts))
    !$ser data echam_vdf_chk_h_ts_rad=field%ts_rad IF (ASSOCIATED(field%ts_rad))
    !$ser data echam_vdf_chk_h_ts_tile=field%ts_tile IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_vdf_chk_h_u_stress_tile=field%u_stress_tile IF (ASSOCIATED(field%u_stress_tile))
    !$ser data echam_vdf_chk_h_v_stress_tile=field%v_stress_tile IF (ASSOCIATED(field%v_stress_tile))
    !$ser data echam_vdf_chk_h_lhflx_tile=field%lhflx_tile IF (ASSOCIATED(field%lhflx_tile))
    !$ser data echam_vdf_chk_h_shflx_tile=field%shflx_tile IF (ASSOCIATED(field%shflx_tile))
    !$ser data echam_vdf_chk_h_evap_tile=field%evap_tile IF (ASSOCIATED(field%evap_tile))
    !$ser data echam_vdf_chk_h_lwflxsfc_tile=field%lwflxsfc_tile IF (ASSOCIATED(field%lwflxsfc_tile))
    !$ser data echam_vdf_chk_h_swflxsfc_tile=field%swflxsfc_tile IF (ASSOCIATED(field%swflxsfc_tile))
    !$ser data echam_vdf_chk_h_albvisdir_tile=field%albvisdir_tile IF (ASSOCIATED(field%albvisdir_tile))
    !$ser data echam_vdf_chk_h_albnirdir_tile=field%albnirdir_tile IF (ASSOCIATED(field%albnirdir_tile))
    !$ser data echam_vdf_chk_h_albvisdif_tile=field%albvisdif_tile IF (ASSOCIATED(field%albvisdif_tile))
    !$ser data echam_vdf_chk_h_albnirdif_tile=field%albnirdif_tile IF (ASSOCIATED(field%albnirdif_tile))
    !$ser data echam_vdf_chk_h_albedo_tile=field%albedo_tile IF (ASSOCIATED(field%albedo_tile))
    !$ser data echam_vdf_chk_h_co2_flux_tile=field%co2_flux_tile IF (ASSOCIATED(field%co2_flux_tile))
    !$ser data echam_vdf_chk_h_q_snocpymlt=field%q_snocpymlt IF (ASSOCIATED(field%q_snocpymlt))
    !$ser data echam_vdf_chk_h_q_rlw_impl=field%q_rlw_impl IF (ASSOCIATED(field%q_rlw_impl))
    !$ser data echam_vdf_chk_h_Tsurf=field%Tsurf IF (ASSOCIATED(field%Tsurf))
    !$ser data echam_vdf_chk_h_T1=field%T1 IF (ASSOCIATED(field%T1))
    !$ser data echam_vdf_chk_h_T2=field%T2 IF (ASSOCIATED(field%T2))
    !$ser data echam_vdf_chk_h_Qtop=field%Qtop IF (ASSOCIATED(field%Qtop))
    !$ser data echam_vdf_chk_h_Qbot=field%Qbot IF (ASSOCIATED(field%Qbot))
    !$ser data echam_vdf_chk_h_albvisdir_ice=field%albvisdir_ice IF (ASSOCIATED(field%albvisdir_ice))
    !$ser data echam_vdf_chk_h_albnirdir_ice=field%albnirdir_ice IF (ASSOCIATED(field%albnirdir_ice))
    !$ser data echam_vdf_chk_h_albvisdif_ice=field%albvisdif_ice IF (ASSOCIATED(field%albvisdif_ice))
    !$ser data echam_vdf_chk_h_albnirdif_ice=field%albnirdif_ice IF (ASSOCIATED(field%albnirdif_ice))
    !$ser data echam_vdf_chk_h_totte=field%totte IF (ASSOCIATED(field%totte))
    !$ser data echam_vdf_chk_h_ta_sfc=tend%ta_sfc IF (ASSOCIATED(tend%ta_sfc))
    !$ser data echam_vdf_chk_h_ta_rlw_impl=tend%ta_rlw_impl IF (ASSOCIATED(tend%ta_rlw_impl))
    !$ser data echam_vdf_chk_h_totte=field%totte IF (ASSOCIATED(field%totte))
    !$ser data echam_vdf_chk_h_z0m=field%z0m IF (ASSOCIATED(field%z0m))
    !$ser data echam_vdf_chk_h_z0m_tile=field%z0m_tile IF (ASSOCIATED(field%z0m_tile))
    !$ser data echam_vdf_chk_h_kedisp=field%kedisp IF (ASSOCIATED(field%kedisp))
    !$ser data echam_vdf_chk_h_sh_vdiff=field%sh_vdiff IF (ASSOCIATED(field%sh_vdiff))
    !$ser data echam_vdf_chk_h_qv_vdiff=field%qv_vdiff IF (ASSOCIATED(field%qv_vdiff))
    !$ser data echam_vdf_chk_h_q_vdf=field%q_vdf IF (ASSOCIATED(field%q_vdf))
    !$ser data echam_vdf_chk_h_q_vdf_vi=field%q_vdf_vi IF (ASSOCIATED(field%q_vdf_vi))
    !$ser data echam_vdf_chk_h_ta_vdf=tend%ta_vdf IF (ASSOCIATED(tend%ta_vdf))
    !$ser data echam_vdf_chk_h_ua_vdf=tend%ua_vdf IF (ASSOCIATED(tend%ua_vdf))
    !$ser data echam_vdf_chk_h_va_vdf=tend%va_vdf IF (ASSOCIATED(tend%va_vdf))
    !$ser data echam_vdf_chk_h_qtrc_vdf=tend%qtrc_vdf IF (ASSOCIATED(tend%qtrc_vdf))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_chk_H_output

  SUBROUTINE serialize_vdf_vd_input(jg, jb, jcs, kproma, kbdim, klev, klevm1, klevp1, ktrac, &
                             ksfc_type, idx_wtr, idx_ice, idx_lnd, pdtime,            &
                             field, pxm1, pxt_emis, pthvvar, pxvar, pwstar,           &
                             pqsat_tile, phdtcbl, pri, pri_tile, pmixlen, pcfm,       &
                             pcfm_tile, pcfh, pcfh_tile, pcfv, pcftotte, pcfthv,      &
                             aa, aa_btm, bb, bb_btm, pfactor_sfc, pcpt_tile,          &
                             pcptgz, pzthvvar, pztottevn)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, kproma
    INTEGER, INTENT(IN)    :: kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN)    :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      & pxm1(:,:),         &
      & pxt_emis(:,:),     &
      & pthvvar(:,:),      &
      & pxvar(:,:),        &
      & pwstar(:),         &
      & pqsat_tile(:,:),   &
      & phdtcbl(:),        &
      & pri(:,:),          &
      & pri_tile(:,:),     &
      & pmixlen(:,:),      &
      & pcfm(:,:),         &
      & pcfm_tile(:,:),    &
      & pcfh(:,:),         &
      & pcfh_tile(:,:),    &
      & pcfv(:,:),         &
      & pcftotte (:,:),    &
      & pcfthv (:,:),      &
      & aa(:,:,:,:),       &
      & aa_btm(:,:,:,:),   &
      & bb(:,:,:),         &
      & bb_btm(:,:,:),     &
      & pfactor_sfc(:),    &
      & pcpt_tile(:,:),    &
      & pcptgz(:,:),       &
      & pzthvvar(:,:),     &
      & pztottevn(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_vd:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_vd:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%coriol )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh )
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%ocu )
    !$ser verbatim   !$ACC UPDATE HOST( field%ocv )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua )
    !$ser verbatim   !$ACC UPDATE HOST( field%va )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( pxm1 )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old )
    !$ser verbatim   !$ACC UPDATE HOST( field%tv )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclc )
    !$ser verbatim   !$ACC UPDATE HOST( pxt_emis ) IF (SIZE(pxt_emis)>0)
    !$ser verbatim   !$ACC UPDATE HOST( pthvvar )
    !$ser verbatim   !$ACC UPDATE HOST( pxvar )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%tottem1 )
    !$ser verbatim   !$ACC UPDATE HOST( field%ustar )
    !$ser verbatim   !$ACC UPDATE HOST( pwstar )
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pqsat_tile )
    !$ser verbatim   !$ACC UPDATE HOST( phdtcbl )
    !$ser verbatim   !$ACC UPDATE HOST( pri )
    !$ser verbatim   !$ACC UPDATE HOST( pri_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pmixlen )
    !$ser verbatim   !$ACC UPDATE HOST( pcfm )
    !$ser verbatim   !$ACC UPDATE HOST( pcfm_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfh )
    !$ser verbatim   !$ACC UPDATE HOST( pcfh_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfv )
    !$ser verbatim   !$ACC UPDATE HOST( pcftotte )
    !$ser verbatim   !$ACC UPDATE HOST( pcfthv )
    !$ser verbatim   !$ACC UPDATE HOST( aa )
    !$ser verbatim   !$ACC UPDATE HOST( aa_btm )
    !$ser verbatim   !$ACC UPDATE HOST( bb )
    !$ser verbatim   !$ACC UPDATE HOST( bb_btm )
    !$ser verbatim   !$ACC UPDATE HOST( pfactor_sfc )
    !$ser verbatim   !$ACC UPDATE HOST( pcpt_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcptgz )
    !$ser verbatim   !$ACC UPDATE HOST( pzthvvar )
    !$ser verbatim   !$ACC UPDATE HOST( field%thvsig )
    !$ser verbatim   !$ACC UPDATE HOST( pztottevn )
    !$ser verbatim   !$ACC UPDATE HOST( field%csat )
    !$ser verbatim   !$ACC UPDATE HOST( field%cair )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0h_lnd )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_vd-input jg=jg jb=jb jcs=jcs kproma=kproma &
    !$ser&          kbdim=kbdim klev=klev klevm1=klevm1 klevp1=klevp1 &
    !$ser&          ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd iqv=iqv iqc=iqc iqi=iqi &
    !$ser&         iqt=iqt pdtime=pdtime date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_vdf_vd_pcoriol=field%coriol(:,jb)           &
    !$ser&     echam_vdf_vd_pzf=field%zf(:,:,jb)                 &
    !$ser&     echam_vdf_vd_pzh=field%zh(:,:,jb)                 &
    !$ser&     echam_vdf_vd_pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     echam_vdf_vd_ptsfc_tile=field%ts_tile(:,jb,:)     &
    !$ser&     echam_vdf_vd_pocu=field%ocu(:,jb)                 &
    !$ser&     echam_vdf_vd_pocv=field%ocv(:,jb)                 &
    !$ser&     echam_vdf_vd_ppsfc=field%presi_old(:,klevp1,jb)   &
    !$ser&     echam_vdf_vd_pum1=field%ua(:,:,jb)                &
    !$ser&     echam_vdf_vd_pvm1=field%va(:,:,jb)                &
    !$ser&     echam_vdf_vd_ptm1=field%ta(:,:,jb)                &
    !$ser&     echam_vdf_vd_pqm1=field%qtrc(:,:,jb,iqv)          &
    !$ser&     echam_vdf_vd_pxlm1=field%qtrc(:,:,jb,iqc)         &
    !$ser&     echam_vdf_vd_pxim1=field%qtrc(:,:,jb,iqi)         &
    !$ser&     echam_vdf_vd_pxm1=pxm1                            &
    !$ser&     echam_vdf_vd_pmair=field%mair(:,:,jb)             &
    !$ser&     echam_vdf_vd_pmref=field%mref(:,:,jb)             &
    !$ser&     echam_vdf_vd_paphm1=field%presi_old(:,:,jb)       &
    !$ser&     echam_vdf_vd_papm1=field%presm_old(:,:,jb)        &
    !$ser&     echam_vdf_vd_ptvm1=field%tv(:,:,jb)               &
    !$ser&     echam_vdf_vd_paclc=field%aclc(:,:,jb)             &
    !$ser&     echam_vdf_vd_pthvvar=pthvvar                      &
    !$ser&     echam_vdf_vd_pxvar=pxvar                          &
    !$ser&     echam_vdf_vd_pz0m_tile=field%z0m_tile(:,jb,:)     &
    !$ser&     echam_vdf_vd_ptottem1=field%tottem1(:,:,jb)       &
    !$ser&     echam_vdf_vd_pustar=field%ustar(:,jb)             &
    !$ser&     echam_vdf_vd_pwstar=pwstar                        &
    !$ser&     echam_vdf_vd_pwstar_tile=field%wstar_tile(:,jb,:) &
    !$ser&     echam_vdf_vd_pqsat_tile=pqsat_tile                &
    !$ser&     echam_vdf_vd_phdtcbl=phdtcbl                      &
    !$ser&     echam_vdf_vd_pri=pri                              &
    !$ser&     echam_vdf_vd_pri_tile=pri_tile                    &
    !$ser&     echam_vdf_vd_pmixlen=pmixlen                      &
    !$ser&     echam_vdf_vd_pcfm=pcfm                            &
    !$ser&     echam_vdf_vd_pcfm_tile=pcfm_tile                  &
    !$ser&     echam_vdf_vd_pcfh=pcfh                            &
    !$ser&     echam_vdf_vd_pcfh_tile=pcfh_tile                  &
    !$ser&     echam_vdf_vd_pcfv=pcfv                            &
    !$ser&     echam_vdf_vd_pcftotte=pcftotte                    &
    !$ser&     echam_vdf_vd_pcfthv=pcfthv                        &
    !$ser&     echam_vdf_vd_aa=aa                                &
    !$ser&     echam_vdf_vd_aa_btm=aa_btm                        &
    !$ser&     echam_vdf_vd_bb=bb                                &
    !$ser&     echam_vdf_vd_bb_btm=bb_btm                        &
    !$ser&     echam_vdf_vd_pfactor_sfc=pfactor_sfc              &
    !$ser&     echam_vdf_vd_pcpt_tile=pcpt_tile                  &
    !$ser&     echam_vdf_vd_pcptgz=pcptgz                        &
    !$ser&     echam_vdf_vd_pzthvvar=pzthvvar                    &
    !$ser&     echam_vdf_vd_pthvsig=field%thvsig(:,jb)           &
    !$ser&     echam_vdf_vd_pztottevn=pztottevn                  &
    !$ser&     echam_vdf_vd_pcsat=field%csat(:,jb)               &
    !$ser&     echam_vdf_vd_pcair=field%cair(:,jb)               &
    !$ser&     echam_vdf_vd_paz0lh=field%z0h_lnd(:,jb)
    !$ser data echam_vdf_vd_pxt_emis=pxt_emis IF (SIZE(pxt_emis)>0)
    !$ser data echam_vdf_vd_pxtm1=field%qtrc(:,:,jb,iqt:) IF (SIZE(field%qtrc(:,:,jb,iqt:))>0)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_vdf_vd:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%coriol )
    !$ser verbatim !$ACC UPDATE DEVICE( field%zf )
    !$ser verbatim !$ACC UPDATE DEVICE( field%zh )
    !$ser verbatim !$ACC UPDATE DEVICE( field%frac_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ts_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ocu )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ocv )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presi_old )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ua )
    !$ser verbatim !$ACC UPDATE DEVICE( field%va )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ta )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( pxm1 )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%mair )
    !$ser verbatim !$ACC UPDATE DEVICE( field%mref )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presi_old )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presm_old )
    !$ser verbatim !$ACC UPDATE DEVICE( field%tv )
    !$ser verbatim !$ACC UPDATE DEVICE( field%aclc )
    !$ser verbatim !$ACC UPDATE DEVICE( pxt_emis ) IF (SIZE(pxt_emis)>0)
    !$ser verbatim !$ACC UPDATE DEVICE( pthvvar )
    !$ser verbatim !$ACC UPDATE DEVICE( pxvar )
    !$ser verbatim !$ACC UPDATE DEVICE( field%z0m_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%tottem1 )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ustar )
    !$ser verbatim !$ACC UPDATE DEVICE( pwstar )
    !$ser verbatim !$ACC UPDATE DEVICE( field%wstar_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pqsat_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( phdtcbl )
    !$ser verbatim !$ACC UPDATE DEVICE( pri )
    !$ser verbatim !$ACC UPDATE DEVICE( pri_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pmixlen )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfm )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfm_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfh )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfh_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfv )
    !$ser verbatim !$ACC UPDATE DEVICE( pcftotte )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfthv )
    !$ser verbatim !$ACC UPDATE DEVICE( aa )
    !$ser verbatim !$ACC UPDATE DEVICE( aa_btm )
    !$ser verbatim !$ACC UPDATE DEVICE( bb )
    !$ser verbatim !$ACC UPDATE DEVICE( bb_btm )
    !$ser verbatim !$ACC UPDATE DEVICE( pfactor_sfc )
    !$ser verbatim !$ACC UPDATE DEVICE( pcpt_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pcptgz )
    !$ser verbatim !$ACC UPDATE DEVICE( pzthvvar )
    !$ser verbatim !$ACC UPDATE DEVICE( field%thvsig )
    !$ser verbatim !$ACC UPDATE DEVICE( pztottevn )
    !$ser verbatim !$ACC UPDATE DEVICE( field%csat )
    !$ser verbatim !$ACC UPDATE DEVICE( field%cair )
    !$ser verbatim !$ACC UPDATE DEVICE( field%z0h_lnd )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_vd_input

  SUBROUTINE serialize_vdf_vd_output(jg, jb, jcs, kproma, kbdim, klev, klevm1, klevp1, ktrac, &
                             ksfc_type, idx_wtr, idx_ice, idx_lnd, pdtime,             &
                             field, pwstar, pqsat_tile, phdtcbl, pri, pri_tile,        &
                             pmixlen, pcfm, pcfm_tile, pcfh, pcfh_tile, pcfv,          &
                             pcftotte, pcfthv, aa, aa_btm, bb, bb_btm,                 &
                             pfactor_sfc, pcpt_tile, pcptgz, pzthvvar,                 &
                             pztottevn, pch_tile, pbn_tile, pbhn_tile,                 &
                             pbm_tile, pbh_tile)
    INTEGER, INTENT(IN)     :: jg, jb, jcs, kproma
    INTEGER, INTENT(IN)     :: kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN)     :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp), INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp), INTENT(INOUT) :: &
      pwstar(:),            &
      pqsat_tile(:,:),      &
      phdtcbl(:),           &
      pri(:,:),             &
      pri_tile(:,:),        &
      pmixlen(:,:),         &
      pcfm(:,:),            &
      pcfm_tile(:,:),       &
      pcfh(:,:),            &
      pcfh_tile(:,:),       &
      pcfv(:,:),            &
      pcftotte (:,:),       &
      pcfthv (:,:),         &
      aa(:,:,:,:),          &
      aa_btm(:,:,:,:),      &
      bb(:,:,:),            &
      bb_btm(:,:,:),        &
      pfactor_sfc(:),       &
      pcpt_tile(:,:),       &
      pcptgz(:,:),          &
      pzthvvar(:,:),        &
      pztottevn(:,:),       &
      pch_tile(:,:),        &
      pbn_tile(:,:),        &
      pbhn_tile(:,:),       &
      pbm_tile(:,:),        &
      pbh_tile(:,:)


    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_vd:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_vd:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ustar )
    !$ser verbatim   !$ACC UPDATE HOST( pwstar )
    !$ser verbatim   !$ACC UPDATE HOST( field%wstar_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pqsat_tile )
    !$ser verbatim   !$ACC UPDATE HOST( phdtcbl )
    !$ser verbatim   !$ACC UPDATE HOST( pri )
    !$ser verbatim   !$ACC UPDATE HOST( pri_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pmixlen )
    !$ser verbatim   !$ACC UPDATE HOST( pcfm )
    !$ser verbatim   !$ACC UPDATE HOST( pcfm_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfh )
    !$ser verbatim   !$ACC UPDATE HOST( pcfh_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfv )
    !$ser verbatim   !$ACC UPDATE HOST( pcftotte )
    !$ser verbatim   !$ACC UPDATE HOST( pcfthv )
    !$ser verbatim   !$ACC UPDATE HOST( aa )
    !$ser verbatim   !$ACC UPDATE HOST( aa_btm )
    !$ser verbatim   !$ACC UPDATE HOST( bb )
    !$ser verbatim   !$ACC UPDATE HOST( bb_btm )
    !$ser verbatim   !$ACC UPDATE HOST( pfactor_sfc )
    !$ser verbatim   !$ACC UPDATE HOST( pcpt_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcptgz )
    !$ser verbatim   !$ACC UPDATE HOST( pzthvvar )
    !$ser verbatim   !$ACC UPDATE HOST( field%thvsig )
    !$ser verbatim   !$ACC UPDATE HOST( pztottevn )
    !$ser verbatim   !$ACC UPDATE HOST( pch_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbn_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbhn_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbm_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbh_tile )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_vd-output jg=jg jb=jb jcs=jcs kproma=kproma &
    !$ser&          kbdim=kbdim klev=klev klevm1=klevm1 klevp1=klevp1 &
    !$ser&          ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd iqv=iqv iqc=iqc iqi=iqi &
    !$ser&          iqt=iqt pdtime=pdtime date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_vd_pustar=field%ustar(:,jb)             &
    !$ser&     echam_vdf_vd_pwstar=pwstar                        &
    !$ser&     echam_vdf_vd_pwstar_tile=field%wstar_tile(:,jb,:) &
    !$ser&     echam_vdf_vd_pqsat_tile=pqsat_tile                &
    !$ser&     echam_vdf_vd_phdtcbl=phdtcbl                      &
    !$ser&     echam_vdf_vd_pri=pri                              &
    !$ser&     echam_vdf_vd_pri_tile=pri_tile                    &
    !$ser&     echam_vdf_vd_pmixlen=pmixlen                      &
    !$ser&     echam_vdf_vd_pcfm=pcfm                            &
    !$ser&     echam_vdf_vd_pcfm_tile=pcfm_tile                  &
    !$ser&     echam_vdf_vd_pcfh=pcfh                            &
    !$ser&     echam_vdf_vd_pcfh_tile=pcfh_tile                  &
    !$ser&     echam_vdf_vd_pcfv=pcfv                            &
    !$ser&     echam_vdf_vd_pcftotte=pcftotte                    &
    !$ser&     echam_vdf_vd_pcfthv=pcfthv                        &
    !$ser&     echam_vdf_vd_aa=aa                                &
    !$ser&     echam_vdf_vd_aa_btm=aa_btm                        &
    !$ser&     echam_vdf_vd_bb=bb                                &
    !$ser&     echam_vdf_vd_bb_btm=bb_btm                        &
    !$ser&     echam_vdf_vd_pfactor_sfc=pfactor_sfc              &
    !$ser&     echam_vdf_vd_pcpt_tile=pcpt_tile                  &
    !$ser&     echam_vdf_vd_pcptgz=pcptgz                        &
    !$ser&     echam_vdf_vd_pzthvvar=pzthvvar                    &
    !$ser&     echam_vdf_vd_pthvsig=field%thvsig(:,jb)           &
    !$ser&     echam_vdf_vd_pztottevn=pztottevn                  &
    !$ser&     echam_vdf_vd_pch_tile=pch_tile                    &
    !$ser&     echam_vdf_vd_pbn_tile=pbn_tile                    &
    !$ser&     echam_vdf_vd_pbhn_tile=pbhn_tile                  &
    !$ser&     echam_vdf_vd_pbm_tile=pbm_tile                    &
    !$ser&     echam_vdf_vd_pbh_tile=pbh_tile
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_vd_output

  SUBROUTINE serialize_vdf_us_input(jb, jg, jcs, kproma, kbdim, klev, klevp1, ksfc_type, &
                             idx_wtr, idx_ice, idx_lnd, pdtime, field,        &
                             pcfh_tile, pcfm_tile, pfac_sfc, aa, aa_btm, bb,  &
                             bb_btm, pcpt_tile, pqsat_tile, nblock,      &
                             pco2, pch_tile)
    INTEGER, INTENT(IN)              :: jb, jg, jcs, kproma, kbdim, klev, klevp1, ksfc_type
    INTEGER, INTENT(IN)              :: idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN)              :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
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
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_update_surface:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_update_surface:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfh_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfm_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pfac_sfc )
    !$ser verbatim   !$ACC UPDATE HOST( field%ocu )
    !$ser verbatim   !$ACC UPDATE HOST( field%ocv )
    !$ser verbatim   !$ACC UPDATE HOST( aa )
    !$ser verbatim   !$ACC UPDATE HOST( aa_btm )
    !$ser verbatim   !$ACC UPDATE HOST( bb )
    !$ser verbatim   !$ACC UPDATE HOST( bb_btm )
    !$ser verbatim   !$ACC UPDATE HOST( pcpt_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pqsat_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%lsmask )
    !$ser verbatim   !$ACC UPDATE HOST( field%alake )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua )
    !$ser verbatim   !$ACC UPDATE HOST( field%va )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( pco2 )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfl )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfc )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfl )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfc )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlds )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsds )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsus )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dir )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dir )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dir )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dif )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dif )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dif )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old )
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0 )
    !$ser verbatim   !$ACC UPDATE HOST( pch_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%csat )
    !$ser verbatim   !$ACC UPDATE HOST( field%cair )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0h_lnd )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo )
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%co2_flux_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%swflxsfc_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%lwflxsfc_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%Tsurf )
    !$ser verbatim   !$ACC UPDATE HOST( field%T1 )
    !$ser verbatim   !$ACC UPDATE HOST( field%T2 )
    !$ser verbatim   !$ACC UPDATE HOST( field%hi )
    !$ser verbatim   !$ACC UPDATE HOST( field%hs )
    !$ser verbatim   !$ACC UPDATE HOST( field%Qtop )
    !$ser verbatim   !$ACC UPDATE HOST( field%Qbot )
    !$ser verbatim   !$ACC UPDATE HOST( field%conc )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_ice )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_ice )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_ice )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_ice )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('init')
    !$ser savepoint echam_vdf_us-input jb=jb jg=jg jcs=jcs kproma=kproma kbdim=kbdim &
    !$ser&          kice=field%kice klev=klev klevp1=klevp1 ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd pdtime=pdtime iqv=iqv date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_vdf_us_pfrc=field%frac_tile(:,jb,:)                 &
    !$ser&     echam_vdf_us_pcfh_tile=pcfh_tile                          &
    !$ser&     echam_vdf_us_pcfm_tile=pcfm_tile                          &
    !$ser&     echam_vdf_us_pfac_sfc=pfac_sfc                            &
    !$ser&     echam_vdf_us_pocu=field%ocu(:,jb)                         &
    !$ser&     echam_vdf_us_pocv=field%ocv(:,jb)                         &
    !$ser&     echam_vdf_us_aa=aa                                        &
    !$ser&     echam_vdf_us_aa_btm=aa_btm                                &
    !$ser&     echam_vdf_us_bb=bb                                        &
    !$ser&     echam_vdf_us_bb_btm=bb_btm                                &
    !$ser&     echam_vdf_us_pcpt_tile=pcpt_tile                          &
    !$ser&     echam_vdf_us_pqsat_tile=pqsat_tile                        &
    !$ser&     echam_vdf_us_ptsfc_tile=field%ts_tile(:,jb,:)             &
    !$ser&     echam_vdf_us_plhflx_tile=field%lhflx_tile(:,jb,:)         &
    !$ser&     echam_vdf_us_pshflx_tile=field%shflx_tile(:,jb,:)         &
    !$ser&     echam_vdf_us_lsm=field%lsmask(:,jb)                       &
    !$ser&     echam_vdf_us_alake=field%alake(:,jb)                      &
    !$ser&     echam_vdf_us_pu=field%ua(:,klev,jb)                       &
    !$ser&     echam_vdf_us_pv=field%va(:,klev,jb)                       &
    !$ser&     echam_vdf_us_ptemp=field%ta(:,klev,jb)                    &
    !$ser&     echam_vdf_us_pq=field%qtrc(:,klev,jb,iqv)                 &
    !$ser&     echam_vdf_us_pco2=pco2                                    &
    !$ser&     echam_vdf_us_prsfl=field%rsfl(:,jb)                       &
    !$ser&     echam_vdf_us_prsfc=field%rsfc(:,jb)                       &
    !$ser&     echam_vdf_us_pssfl=field%ssfl(:,jb)                       &
    !$ser&     echam_vdf_us_pssfc=field%ssfc(:,jb)                       &
    !$ser&     echam_vdf_us_rlds=field%rlds(:,jb)                        &
    !$ser&     echam_vdf_us_rlus=field%rlus(:,jb)                        &
    !$ser&     echam_vdf_us_rsds=field%rsds(:,jb)                        &
    !$ser&     echam_vdf_us_rsus=field%rsus(:,jb)                        &
    !$ser&     echam_vdf_us_rvds_dir=field%rvds_dir(:,jb)                &
    !$ser&     echam_vdf_us_rpds_dir=field%rpds_dir(:,jb)                &
    !$ser&     echam_vdf_us_rnds_dir=field%rnds_dir(:,jb)                &
    !$ser&     echam_vdf_us_rvds_dif=field%rvds_dif(:,jb)                &
    !$ser&     echam_vdf_us_rpds_dif=field%rpds_dif(:,jb)                &
    !$ser&     echam_vdf_us_rnds_dif=field%rnds_dif(:,jb)                &
    !$ser&     echam_vdf_us_ps=field%presi_old(:,klevp1,jb)              &
    !$ser&     echam_vdf_us_pcosmu0=field%cosmu0(:,jb)                   &
    !$ser&     echam_vdf_us_pch_tile=pch_tile                            &
    !$ser&     echam_vdf_us_pcsat=field%csat(:,jb)                       &
    !$ser&     echam_vdf_us_pcair=field%cair(:,jb)                       &
    !$ser&     echam_vdf_us_z0m_tile=field%z0m_tile(:,jb,:)              &
    !$ser&     echam_vdf_us_z0h_lnd=field%z0h_lnd(:,jb)                  &
    !$ser&     echam_vdf_us_albvisdir=field%albvisdir(:,jb)              &
    !$ser&     echam_vdf_us_albnirdir=field%albnirdir(:,jb)              &
    !$ser&     echam_vdf_us_albvisdif=field%albvisdif(:,jb)              &
    !$ser&     echam_vdf_us_albnirdif=field%albnirdif(:,jb)              &
    !$ser&     echam_vdf_us_albvisdir_tile=field%albvisdir_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albnirdir_tile=field%albnirdir_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albvisdif_tile=field%albvisdif_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albnirdif_tile=field%albnirdif_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albedo=field%albedo(:,jb)                    &
    !$ser&     echam_vdf_us_albedo_tile=field%albedo_tile(:,jb,:)        &
    !$ser&     echam_vdf_us_pco2_flux_tile=field%co2_flux_tile(:,jb,:)   &
    !$ser&     echam_vdf_us_rsns_tile=field%swflxsfc_tile(:,jb,:)        &
    !$ser&     echam_vdf_us_rlns_tile=field%lwflxsfc_tile(:,jb,:)        &
    !$ser&     echam_vdf_us_Tsurf=field%Tsurf(:,:,jb)                    &
    !$ser&     echam_vdf_us_T1=field%T1(:,:,jb)                          &
    !$ser&     echam_vdf_us_T2=field%T2(:,:,jb)                          &
    !$ser&     echam_vdf_us_hi=field%hi(:,:,jb)                          &
    !$ser&     echam_vdf_us_hs=field%hs(:,:,jb)                          &
    !$ser&     echam_vdf_us_Qtop=field%Qtop(:,:,jb)                      &
    !$ser&     echam_vdf_us_Qbot=field%Qbot(:,:,jb)                      &
    !$ser&     echam_vdf_us_conc=field%conc(:,:,jb)                      &
    !$ser&     echam_vdf_us_albvisdir_ice=field%albvisdir_ice(:,:,jb)    &
    !$ser&     echam_vdf_us_albnirdir_ice=field%albnirdir_ice(:,:,jb)    &
    !$ser&     echam_vdf_us_albvisdif_ice=field%albvisdif_ice(:,:,jb)    &
    !$ser&     echam_vdf_us_albnirdif_ice=field%albnirdif_ice(:,:,jb)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined ( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_vdf_us:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%frac_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfh_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfm_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pfac_sfc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ocu )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ocv )
    !$ser verbatim !$ACC UPDATE DEVICE( aa )
    !$ser verbatim !$ACC UPDATE DEVICE( aa_btm )
    !$ser verbatim !$ACC UPDATE DEVICE( bb )
    !$ser verbatim !$ACC UPDATE DEVICE( bb_btm )
    !$ser verbatim !$ACC UPDATE DEVICE( pcpt_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pqsat_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ts_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%lhflx_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%shflx_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%lsmask )
    !$ser verbatim !$ACC UPDATE DEVICE( field%alake )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ua )
    !$ser verbatim !$ACC UPDATE DEVICE( field%va )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ta )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( pco2 )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rsfl )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rsfc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ssfl )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ssfc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rlds )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rlus )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rsds )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rsus )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rvds_dir )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rpds_dir )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rnds_dir )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rvds_dif )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rpds_dif )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rnds_dif )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presi_old )
    !$ser verbatim !$ACC UPDATE DEVICE( field%cosmu0 )
    !$ser verbatim !$ACC UPDATE DEVICE( pch_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%csat )
    !$ser verbatim !$ACC UPDATE DEVICE( field%cair )
    !$ser verbatim !$ACC UPDATE DEVICE( field%z0m_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%z0h_lnd )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albvisdir )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albnirdir )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albvisdif )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albnirdif )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albvisdir_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albnirdir_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albvisdif_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albnirdif_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albedo )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albedo_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%co2_flux_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%swflxsfc_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%lwflxsfc_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%Tsurf )
    !$ser verbatim !$ACC UPDATE DEVICE( field%T1 )
    !$ser verbatim !$ACC UPDATE DEVICE( field%T2 )
    !$ser verbatim !$ACC UPDATE DEVICE( field%hi )
    !$ser verbatim !$ACC UPDATE DEVICE( field%hs )
    !$ser verbatim !$ACC UPDATE DEVICE( field%Qtop )
    !$ser verbatim !$ACC UPDATE DEVICE( field%Qbot )
    !$ser verbatim !$ACC UPDATE DEVICE( field%conc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albvisdir_ice )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albnirdir_ice )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albvisdif_ice )
    !$ser verbatim !$ACC UPDATE DEVICE( field%albnirdif_ice )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_us_input

  SUBROUTINE serialize_vdf_us_output(jb, jg, jcs, kproma, kbdim, klev, ksfc_type, idx_wtr,     &
                              idx_ice, idx_lnd, pdtime, field, aa, aa_btm, bb, &
                              bb_btm, pcpt_tile, pqsat_tile, q_snocpymlt)
    INTEGER, INTENT(IN)       :: jb, jg, jcs, kproma, kbdim, klev, ksfc_type
    INTEGER, INTENT(IN)       :: idx_wtr, idx_ice, idx_lnd
    REAL(wp), INTENT(IN)      :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp), INTENT(INOUT)  :: aa(:,:,:,:)
    REAL(wp), INTENT(INOUT)  :: aa_btm(:,:,:,:)
    REAL(wp), INTENT(INOUT)  :: bb(:,:,:)
    REAL(wp), INTENT(INOUT)  :: bb_btm(:,:,:)
    REAL(wp), INTENT(INOUT)  :: pcpt_tile(:,:)
    REAL(wp), INTENT(INOUT)  :: pqsat_tile(:,:)
    REAL(wp), INTENT(INOUT)  :: q_snocpymlt(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_us:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_us:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( aa )
    !$ser verbatim   !$ACC UPDATE HOST( aa_btm )
    !$ser verbatim   !$ACC UPDATE HOST( bb )
    !$ser verbatim   !$ACC UPDATE HOST( bb_btm )
    !$ser verbatim   !$ACC UPDATE HOST( pcpt_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pqsat_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress )
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress )
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx )
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap )
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%lhflx_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%shflx_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%fco2nat )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus )
    !$ser verbatim   !$ACC UPDATE HOST( field%csat )
    !$ser verbatim   !$ACC UPDATE HOST( field%cair )
    !$ser verbatim   !$ACC UPDATE HOST( q_snocpymlt )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0h_lnd )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo )
    !$ser verbatim   !$ACC UPDATE HOST( field%albedo_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%co2_flux_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad )
    !$ser verbatim   !$ACC UPDATE HOST( field%swflxsfc_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%lwflxsfc_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%lake_ice_frc )
    !$ser verbatim   !$ACC UPDATE HOST( field%Tsurf )
    !$ser verbatim   !$ACC UPDATE HOST( field%T1 )
    !$ser verbatim   !$ACC UPDATE HOST( field%T2 )
    !$ser verbatim   !$ACC UPDATE HOST( field%hi )
    !$ser verbatim   !$ACC UPDATE HOST( field%hs )
    !$ser verbatim   !$ACC UPDATE HOST( field%Qtop )
    !$ser verbatim   !$ACC UPDATE HOST( field%Qbot )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdir_ice )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdir_ice )
    !$ser verbatim   !$ACC UPDATE HOST( field%albvisdif_ice )
    !$ser verbatim   !$ACC UPDATE HOST( field%albnirdif_ice )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('init')
    !$ser savepoint echam_vdf_us-output jb=jb jg=jg jcs=jcs kproma=kproma kbdim=kbdim &
    !$ser&          kice=field%kice klev=klev ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd pdtime=pdtime iqv=iqv date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_us_aa=aa                                        &
    !$ser&     echam_vdf_us_aa_btm=aa_btm                                &
    !$ser&     echam_vdf_us_bb=bb                                        &
    !$ser&     echam_vdf_us_bb_btm=bb_btm                                &
    !$ser&     echam_vdf_us_pcpt_tile=pcpt_tile                          &
    !$ser&     echam_vdf_us_pqsat_tile=pqsat_tile                        &
    !$ser&     echam_vdf_us_ptsfc_tile=field%ts_tile(:,jb,:)             &
    !$ser&     echam_vdf_us_pu_stress_gbm=field%u_stress(:,jb)           &
    !$ser&     echam_vdf_us_pv_stress_gbm=field%v_stress(:,jb)           &
    !$ser&     echam_vdf_us_plhflx_gbm=field%lhflx(:,jb)                 &
    !$ser&     echam_vdf_us_pshflx_gbm=field%shflx(:,jb)                 &
    !$ser&     echam_vdf_us_pevap_gbm=field%evap(:,jb)                   &
    !$ser&     echam_vdf_us_pu_stress_tile=field%u_stress_tile(:,jb,:)   &
    !$ser&     echam_vdf_us_pv_stress_tile=field%v_stress_tile(:,jb,:)   &
    !$ser&     echam_vdf_us_plhflx_tile=field%lhflx_tile(:,jb,:)         &
    !$ser&     echam_vdf_us_pshflx_tile=field%shflx_tile(:,jb,:)         &
    !$ser&     echam_vdf_us_pevap_tile=field%evap_tile(:,jb,:)           &
    !$ser&     echam_vdf_us_pco2nat=field%fco2nat(:,jb)                  &
    !$ser&     echam_vdf_us_rlus=field%rlus(:,jb)                        &
    !$ser&     echam_vdf_us_pcsat=field%csat(:,jb)                       &
    !$ser&     echam_vdf_us_pcair=field%cair(:,jb)                       &
    !$ser&     echam_vdf_us_q_snocpymlt=q_snocpymlt                      &
    !$ser&     echam_vdf_us_z0m_tile=field%z0m_tile(:,jb,:)              &
    !$ser&     echam_vdf_us_z0h_lnd=field%z0h_lnd(:,jb)                  &
    !$ser&     echam_vdf_us_albvisdir=field%albvisdir(:,jb)              &
    !$ser&     echam_vdf_us_albnirdir=field%albnirdir(:,jb)              &
    !$ser&     echam_vdf_us_albvisdif=field%albvisdif(:,jb)              &
    !$ser&     echam_vdf_us_albnirdif=field%albnirdif(:,jb)              &
    !$ser&     echam_vdf_us_albvisdir_tile=field%albvisdir_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albnirdir_tile=field%albnirdir_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albvisdif_tile=field%albvisdif_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albnirdif_tile=field%albnirdif_tile(:,jb,:)  &
    !$ser&     echam_vdf_us_albedo=field%albedo(:,jb)                    &
    !$ser&     echam_vdf_us_albedo_tile=field%albedo_tile(:,jb,:)        &
    !$ser&     echam_vdf_us_co2_flux_tile=field%co2_flux_tile(:,jb,:)    &
    !$ser&     echam_vdf_us_ptsfc=field%ts(:,jb)                         &
    !$ser&     echam_vdf_us_ptsfc_rad=field%ts_rad(:,jb)                 &
    !$ser&     echam_vdf_us_rsns_tile=field%swflxsfc_tile(:,jb,:)        &
    !$ser&     echam_vdf_us_rlns_tile=field%lwflxsfc_tile(:,jb,:)        &
    !$ser&     echam_vdf_us_lake_ice_frc=field%lake_ice_frc(:,jb)        &
    !$ser&     echam_vdf_us_Tsurf=field%Tsurf(:,:,jb)                    &
    !$ser&     echam_vdf_us_T1=field%T1(:,:,jb)                          &
    !$ser&     echam_vdf_us_T2=field%T2(:,:,jb)                          &
    !$ser&     echam_vdf_us_hi=field%hi(:,:,jb)                          &
    !$ser&     echam_vdf_us_hs=field%hs(:,:,jb)                          &
    !$ser&     echam_vdf_us_Qtop=field%Qtop(:,:,jb)                      &
    !$ser&     echam_vdf_us_Qbot=field%Qbot(:,:,jb)                      &
    !$ser&     echam_vdf_us_albvisdir_ice=field%albvisdir_ice(:,:,jb)    &
    !$ser&     echam_vdf_us_albnirdir_ice=field%albnirdir_ice(:,:,jb)    &
    !$ser&     echam_vdf_us_albvisdif_ice=field%albvisdif_ice(:,:,jb)    &
    !$ser&     echam_vdf_us_albnirdif_ice=field%albnirdif_ice(:,:,jb)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_us_output

  SUBROUTINE serialize_vdf_vu_input(jb, jcs, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, &
                             idx_wtr, pdtime, field, pcfm_tile, aa, pcptgz, pztottevn, bb, &
                             pzthvvar, pxvar, pkedisp)
    INTEGER, INTENT(IN)    :: jb, jcs, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, idx_wtr
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      pcfm_tile(:,:),       &
      aa(:,:,:,:),          &
      pcptgz(:,:),          &
      pztottevn(:,:),       &
      bb(:,:,:),            &
      pzthvvar(:,:),        &
      pxvar(:,:),           &
      pkedisp(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_vu:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_vu:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pcfm_tile )
    !$ser verbatim   !$ACC UPDATE HOST( aa )
    !$ser verbatim   !$ACC UPDATE HOST( pcptgz )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua )
    !$ser verbatim   !$ACC UPDATE HOST( field%va )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%geom )
    !$ser verbatim   !$ACC UPDATE HOST( pztottevn )
    !$ser verbatim   !$ACC UPDATE HOST( bb )
    !$ser verbatim   !$ACC UPDATE HOST( pzthvvar )
    !$ser verbatim   !$ACC UPDATE HOST( pxvar )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_echa_vdf_vu-input jb=jb jcs=jcs kproma=kproma kbdim=kbdim klev=klev &
    !$ser&          klevm1=klevm1 ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_vdf_vu_pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     echam_vdf_vu_pcfm_tile=pcfm_tile                  &
    !$ser&     echam_vdf_vu_aa=aa                                &
    !$ser&     echam_vdf_vu_pcptgz=pcptgz                        &
    !$ser&     echam_vdf_vu_pum1=field%ua(:,:,jb)                &
    !$ser&     echam_vdf_vu_pvm1=field%va(:,:,jb)                &
    !$ser&     echam_vdf_vu_ptm1=field%ta(:,:,jb)                &
    !$ser&     echam_vdf_vu_pmair=field%mair(:,:,jb)             &
    !$ser&     echam_vdf_vu_pmref=field%mref(:,:,jb)             &
    !$ser&     echam_vdf_vu_pqm1=field%qtrc(:,:,jb,iqv)          &
    !$ser&     echam_vdf_vu_pxlm1=field%qtrc(:,:,jb,iqc)         &
    !$ser&     echam_vdf_vu_pxim1=field%qtrc(:,:,jb,iqi)         &
    !$ser&     echam_vdf_vu_pgeom1=field%geom(:,:,jb)            &
    !$ser&     echam_vdf_vu_pztottevn=pztottevn                  &
    !$ser&     echam_vdf_vu_bb=bb                                &
    !$ser&     echam_vdf_vu_pzthvvar=pzthvvar                    &
    !$ser&     echam_vdf_vu_pxvar=pxvar                          &
    !$ser&     echam_vdf_vu_pz0m_tile=field%z0m_tile(:,jb,:)
    !$ser data echam_vdf_vu_pxtm1=field%qtrc(:,:,jb,iqt:) IF (SIZE(field%qtrc(:,:,jb,iqt:))>0)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_vdf_vu:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%frac_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pcfm_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( aa )
    !$ser verbatim !$ACC UPDATE DEVICE( pcptgz )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ua )
    !$ser verbatim !$ACC UPDATE DEVICE( field%va )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ta )
    !$ser verbatim !$ACC UPDATE DEVICE( field%mair )
    !$ser verbatim !$ACC UPDATE DEVICE( field%mref )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%geom )
    !$ser verbatim !$ACC UPDATE DEVICE( pztottevn )
    !$ser verbatim !$ACC UPDATE DEVICE( bb )
    !$ser verbatim !$ACC UPDATE DEVICE( pzthvvar )
    !$ser verbatim !$ACC UPDATE DEVICE( pxvar )
    !$ser verbatim !$ACC UPDATE DEVICE( field%z0m_tile )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_vu_input

  SUBROUTINE serialize_vdf_vu_output(jb, jcs, kproma, kbdim, klev, klevm1, ktrac, &
                              ksfc_type, idx_wtr, pdtime, field, bb, pxvar, &
                              pkedisp, pute_vdf, pvte_vdf, pq_vdf, tend_qtrc_vdf, &
                              pthvvar)
    INTEGER, INTENT(IN)    :: jb, jcs, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, idx_wtr
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      bb(:,:,:),            &
      pxvar(:,:),           &
      pkedisp(:),           &
      pute_vdf(:,:),        &
      pvte_vdf(:,:),        &
      pq_vdf(:,:),          &
      tend_qtrc_vdf(:,:,:), &
      pthvvar(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_vu:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_vu:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( bb )
    !$ser verbatim   !$ACC UPDATE HOST( pxvar )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pkedisp )
    !$ser verbatim   !$ACC UPDATE HOST( pute_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( pvte_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( pq_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( tend_qtrc_vdf )
    !$ser verbatim   !$ACC UPDATE HOST( field%z0m )
    !$ser verbatim   !$ACC UPDATE HOST( pthvvar )
    !$ser verbatim   !$ACC UPDATE HOST( field%totte )
    !$ser verbatim   !$ACC UPDATE HOST( field%sh_vdiff )
    !$ser verbatim   !$ACC UPDATE HOST( field%qv_vdiff )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_vu-output jb=jb jcs=jcs kproma=kproma kbdim=kbdim klev=klev &
    !$ser&          klevm1=klevm1 ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_vu_bb=bb                                  &
    !$ser&     echam_vdf_vu_pxvar=pxvar                            &
    !$ser&     echam_vdf_vu_pz0m_tile=field%z0m_tile(:,jb,:)       &
    !$ser&     echam_vdf_vu_pkedisp=pkedisp                        &
    !$ser&     echam_vdf_vu_pute_vdf=pute_vdf                      &
    !$ser&     echam_vdf_vu_pvte_vdf=pvte_vdf                      &
    !$ser&     echam_vdf_vu_pq_vdf=pq_vdf                          &
    !$ser&     echam_vdf_vu_pqte_vdf=tend_qtrc_vdf(:,:,iqv)        &
    !$ser&     echam_vdf_vu_pxlte_vdf=tend_qtrc_vdf(:,:,iqc)       &
    !$ser&     echam_vdf_vu_pxite_vdf=tend_qtrc_vdf(:,:,iqi)       &
    !$ser&     echam_vdf_vu_pz0m=field%z0m(:,jb)                   &
    !$ser&     echam_vdf_vu_pthvvar=pthvvar                        &
    !$ser&     echam_vdf_vu_ptotte=field%totte(:,:,jb)             &
    !$ser&     echam_vdf_vu_psh_vdiff=field%sh_vdiff(:,jb)         &
    !$ser&     echam_vdf_vu_pqv_vdiff=field%qv_vdiff(:,jb)
    !$ser data echam_vdf_vu_pxtte_vdf=tend_qtrc_vdf(:,:,iqt:) IF (SIZE(tend_qtrc_vdf(:,:,iqt:))>0)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_vu_output

  SUBROUTINE serialize_vdf_nd_input(jb, jcs, kproma, kbdim, klev, klevp1, ksfc_type, idx_lnd, field, &
                             pxm1, pcptgz, pcpt_tile, pbn_tile, pbhn_tile, pbh_tile, pbm_tile, pri_tile)
    INTEGER, INTENT(IN)    :: jb, jcs, kproma, kbdim, klev, klevp1, ksfc_type, idx_lnd
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      pxm1(:,:),        &
      pcptgz(:,:),      &
      pcpt_tile(:,:),   &
      pbn_tile(:,:),    &
      pbhn_tile(:,:),   &
      pbh_tile(:,:),    &
      pbm_tile(:,:),    &
      pri_tile(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_nd:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_nd:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old )
    !$ser verbatim   !$ACC UPDATE HOST( pxm1 )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua )
    !$ser verbatim   !$ACC UPDATE HOST( field%va )
    !$ser verbatim   !$ACC UPDATE HOST( field%ocu )
    !$ser verbatim   !$ACC UPDATE HOST( field%ocv )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh )
    !$ser verbatim   !$ACC UPDATE HOST( pcptgz )
    !$ser verbatim   !$ACC UPDATE HOST( pcpt_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbn_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbhn_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbh_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pbm_tile )
    !$ser verbatim   !$ACC UPDATE HOST( pri_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmax )
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmin )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_nd-input jb=jb jcs=jcs kproma=kproma kbdim=kbdim klev=klev klevp1=klevp1 &
    !$ser&          ksfc_type=ksfc_type idx_lnd=idx_lnd iqv=iqv date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_vdf_nd_pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     echam_vdf_nd_pqm1=field%qtrc(:,klev,jb,iqv)       &
    !$ser&     echam_vdf_nd_ptm1=field%ta(:,klev,jb)             &
    !$ser&     echam_vdf_nd_papm1=field%presm_old(:,klev,jb)     &
    !$ser&     echam_vdf_nd_paphm1=field%presi_old(:,klevp1,jb)  &
    !$ser&     echam_vdf_nd_pxm1=pxm1(:,klev)                    &
    !$ser&     echam_vdf_nd_pum1=field%ua(:,klev,jb)             &
    !$ser&     echam_vdf_nd_pvm1=field%va(:,klev,jb)             &
    !$ser&     echam_vdf_nd_pocu=field%ocu(:,jb)                 &
    !$ser&     echam_vdf_nd_pocv=field%ocv(:,jb)                 &
    !$ser&     echam_vdf_nd_pzf=field%zf(:,klev,jb)              &
    !$ser&     echam_vdf_nd_pzs=field%zh(:,klev+1,jb)            &
    !$ser&     echam_vdf_nd_pcptgz=pcptgz(:,klev)                &
    !$ser&     echam_vdf_nd_pcpt_tile=pcpt_tile                  &
    !$ser&     echam_vdf_nd_pbn_tile=pbn_tile                    &
    !$ser&     echam_vdf_nd_pbhn_tile=pbhn_tile                  &
    !$ser&     echam_vdf_nd_pbh_tile=pbh_tile                    &
    !$ser&     echam_vdf_nd_pbm_tile=pbm_tile                    &
    !$ser&     echam_vdf_nd_pri_tile=pri_tile                    &
    !$ser&     echam_vdf_nd_ptasmax=field%tasmax(:,jb)           &
    !$ser&     echam_vdf_nd_ptasmin=field%tasmin(:,jb)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined ( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_vdf_nd:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%frac_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ta )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presm_old )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presi_old )
    !$ser verbatim !$ACC UPDATE DEVICE( pxm1 )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ua )
    !$ser verbatim !$ACC UPDATE DEVICE( field%va )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ocu )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ocv )
    !$ser verbatim !$ACC UPDATE DEVICE( field%zf )
    !$ser verbatim !$ACC UPDATE DEVICE( field%zh )
    !$ser verbatim !$ACC UPDATE DEVICE( pcptgz )
    !$ser verbatim !$ACC UPDATE DEVICE( pcpt_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pbn_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pbhn_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pbh_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pbm_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( pri_tile )
    !$ser verbatim !$ACC UPDATE DEVICE( field%tasmax )
    !$ser verbatim !$ACC UPDATE DEVICE( field%tasmin )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_nd_input

  SUBROUTINE serialize_vdf_nd_output(jb, jcs, kproma, kbdim, ksfc_type, idx_lnd, field)
    INTEGER, INTENT(IN)    :: jb, jcs, kproma, kbdim, ksfc_type, idx_lnd
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .FALSE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_vdf_nd:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_vdf_nd:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%sfcWind )
    !$ser verbatim   !$ACC UPDATE HOST( field%tas )
    !$ser verbatim   !$ACC UPDATE HOST( field%dew2 )
    !$ser verbatim   !$ACC UPDATE HOST( field%uas )
    !$ser verbatim   !$ACC UPDATE HOST( field%vas )
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmax )
    !$ser verbatim   !$ACC UPDATE HOST( field%tasmin )
    !$ser verbatim   !$ACC UPDATE HOST( field%sfcWind_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%tas_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%dew2_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%uas_tile )
    !$ser verbatim   !$ACC UPDATE HOST( field%vas_tile )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdf_nd-output jb=jb jcs=jcs kproma=kproma kbdim=kbdim &
    !$ser&          ksfc_type=ksfc_type idx_lnd=idx_lnd date=TRIM(date)
    !$ser mode write
    !$ser data echam_vdf_nd_psfcWind_gbm=field%sfcWind(:,jb)             &
    !$ser&     echam_vdf_nd_ptas_gbm=field%tas(:,jb)                     &
    !$ser&     echam_vdf_nd_pdew2_gbm=field%dew2(:,jb)                   &
    !$ser&     echam_vdf_nd_puas_gbm=field%uas(:,jb)                     &
    !$ser&     echam_vdf_nd_pvas_gbm=field%vas(:,jb)                     &
    !$ser&     echam_vdf_nd_ptasmax=field%tasmax(:,jb)                   &
    !$ser&     echam_vdf_nd_ptasmin=field%tasmin(:,jb)                   &
    !$ser&     echam_vdf_nd_psfcWind_tile=field%sfcWind_tile(:,jb,:)     &
    !$ser&     echam_vdf_nd_ptas_tile=field%tas_tile(:,jb,:)             &
    !$ser&     echam_vdf_nd_pdew2_tile=field%dew2_tile(:,jb,:)           &
    !$ser&     echam_vdf_nd_puas_tile=field%uas_tile(:,jb,:)             &
    !$ser&     echam_vdf_nd_pvas_tile=field%vas_tile(:,jb,:)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_vdf_nd_output

END MODULE mo_ser_echam_vdf
