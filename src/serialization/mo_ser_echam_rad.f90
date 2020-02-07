!--------------------------------------------------------------------
!
! Serialization routine for ECHAM RAD
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_rad

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  PUBLIC :: serialize_rad_input
  PUBLIC :: serialize_rad_output

  CONTAINS

  SUBROUTINE serialize_rad_input(jg, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
#if defined(SERIALIZE_ECHAM_RAD) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
#if defined( SERIALIZE_CREATE_REFERENCE )
    LOGICAL, SAVE :: lenabled = .TRUE.
#else
    LOGICAL, SAVE :: lenabled = .FALSE.
#endif
    LOGICAL, SAVE :: lactive = .TRUE.

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_rad:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_rad:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !NOACC UPDATE HOST( field%ts_rad ) IF( ASSOCIATED( field%ts_rad )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%ts_rad_rt ) IF( ASSOCIATED( field%ts_rad_rt )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%sftlf ) IF( ASSOCIATED( field%sftlf )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%sftgif ) IF( ASSOCIATED( field%sftgif )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%cosmu0_rt ) IF( ASSOCIATED( field%cosmu0_rt )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%daylght_frc_rt ) IF( ASSOCIATED( field%daylght_frc_rt )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%albvisdir ) IF( ASSOCIATED( field%albvisdir )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%albnirdir ) IF( ASSOCIATED( field%albnirdir )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%albvisdif ) IF( ASSOCIATED( field%albvisdif )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%albnirdif ) IF( ASSOCIATED( field%albnirdif )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%emissivity ) IF( ASSOCIATED( field%emissivity )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%zf ) IF( ASSOCIATED( field%zf )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%zh ) IF( ASSOCIATED( field%zh )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%dz ) IF( ASSOCIATED( field%dz )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED( field%presi_old )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED( field%presm_old )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%ta ) IF( ASSOCIATED( field%ta )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%mdry ) IF( ASSOCIATED( field%mdry )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%mtrc ) IF( ASSOCIATED( field%mtrc ))  
    !$ser verbatim   !NOACC UPDATE HOST( field%o3 ) IF( ASSOCIATED( field%o3 )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%acdnc ) IF( ASSOCIATED( field%acdnc )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%aclc ) IF( ASSOCIATED( field%aclc )) 
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_rad-input jg=jg nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_rad_ts_rad=field%ts_rad(:,:) IF (ASSOCIATED(field%ts_rad)) 
    !$ser data echam_rad_ta_rad_rt=field%ts_rad_rt(:,:) IF (ASSOCIATED(field%ts_rad_rt)) 
    !$ser data echam_rad_sftlf=field%sftlf(:,:) IF (ASSOCIATED(field%sftlf)) 
    !$ser data echam_rad_sftgif=field%sftgif(:,:) IF (ASSOCIATED(field%sftgif)) 
    !$ser data echam_rad_cosmu0_rt=field%cosmu0_rt(:,:) IF (ASSOCIATED(field%cosmu0_rt)) 
    !$ser data echam_rad_daylght_frc_rt=field%daylght_frc_rt(:,:) IF (ASSOCIATED(field%daylght_frc_rt)) 
    !$ser data echam_rad_albvisdir=field%albvisdir(:,:) IF (ASSOCIATED(field%albvisdir)) 
    !$ser data echam_rad_albnirdir=field%albnirdir(:,:) IF (ASSOCIATED(field%albnirdir)) 
    !$ser data echam_rad_albvisdif=field%albvisdif(:,:) IF (ASSOCIATED(field%albvisdif)) 
    !$ser data echam_rad_albnirdif=field%albnirdif(:,:) IF (ASSOCIATED(field%albnirdif)) 
    !$ser data echam_rad_emissivity=field%emissivity(:,:) IF (ASSOCIATED(field%emissivity)) 
    !$ser data echam_rad_zf=field%zf(:,:,:) IF (ASSOCIATED(field%zf)) 
    !$ser data echam_rad_zh=field%zh(:,:,:) IF (ASSOCIATED(field%zh)) 
    !$ser data echam_rad_dz=field%dz(:,:,:) IF (ASSOCIATED(field%dz)) 
    !$ser data echam_rad_presi_old=field%presi_old(:,:,:) IF (ASSOCIATED(field%presi_old)) 
    !$ser data echam_rad_presm_old=field%presm_old(:,:,:) IF (ASSOCIATED(field%presm_old)) 
    !$ser data echam_rad_ta=field%ta(:,:,:) IF (ASSOCIATED(field%ta)) 
    !$ser data echam_rad_mdry=field%mdry(:,:,:) IF (ASSOCIATED(field%mdry)) 
    !$ser data echam_rad_mtrc=field%mtrc(:,:,:,:) IF (ASSOCIATED(field%mtrc))
    !$ser data echam_rad_o3=field%o3(:,:,:) IF (ASSOCIATED(field%o3)) 
    !$ser data echam_rad_acdnc=field%acdnc(:,:,:) IF (ASSOCIATED(field%acdnc)) 
    !$ser data echam_rad_aclc=field%aclc(:,:,:) IF (ASSOCIATED(field%aclc)) 
    !$ser data echam_rad_aclcov=field%aclcov(:,:) IF (ASSOCIATED(field%aclcov)) 
    !$ser data echam_rad_rldcs_rt=field%rldcs_rt(:,:,:) IF (ASSOCIATED(field%rldcs_rt))
    !$ser data echam_rad_rlucs_rt=field%rlucs_rt(:,:,:) IF (ASSOCIATED(field%rlucs_rt))
    !$ser data echam_rad_rsdcs_rt=field%rsdcs_rt(:,:,:) IF (ASSOCIATED(field%rsdcs_rt))
    !$ser data echam_rad_rsucs_rt=field%rsucs_rt(:,:,:) IF (ASSOCIATED(field%rsucs_rt))
    !$ser data echam_rad_rld_rt=field%rld_rt(:,:,:) IF (ASSOCIATED(field%rld_rt))
    !$ser data echam_rad_rlu_rt=field%rlu_rt(:,:,:) IF (ASSOCIATED(field%rlu_rt))
    !$ser data echam_rad_rsd_rt=field%rsd_rt(:,:,:) IF (ASSOCIATED(field%rsd_rt))
    !$ser data echam_rad_rsu_rt=field%rsu_rt(:,:,:) IF (ASSOCIATED(field%rsu_rt))
    !$ser data echam_rad_rvds_dir_rt=field%rvds_dir_rt(:,:) IF (ASSOCIATED(field%rvds_dir_rt))
    !$ser data echam_rad_rpds_dir_rt=field%rpds_dir_rt(:,:) IF (ASSOCIATED(field%rpds_dir_rt))
    !$ser data echam_rad_rnds_dir_rt=field%rnds_dir_rt(:,:) IF (ASSOCIATED(field%rnds_dir_rt))
    !$ser data echam_rad_rvds_dif_rt=field%rvds_dif_rt(:,:) IF (ASSOCIATED(field%rvds_dif_rt))
    !$ser data echam_rad_rpds_dif_rt=field%rpds_dif_rt(:,:) IF (ASSOCIATED(field%rpds_dif_rt))
    !$ser data echam_rad_rnds_dif_rt=field%rnds_dif_rt(:,:) IF (ASSOCIATED(field%rnds_dif_rt))
    !$ser data echam_rad_rvus_rt=field%rvus_rt(:,:) IF (ASSOCIATED(field%rvus_rt))
    !$ser data echam_rad_rpus_rt=field%rpus_rt(:,:) IF (ASSOCIATED(field%rpus_rt))
    !$ser data echam_rad_rnus_rt=field%rnus_rt(:,:) IF (ASSOCIATED(field%rnus_rt))
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lactive = .FALSE.
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_rad:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !NOACC UPDATE DEVICE( field%ts_rad ) IF( ASSOCIATED( field%ts_rad )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%ts_rad_rt ) IF( ASSOCIATED( field%ts_rad_rt )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%sftlf ) IF( ASSOCIATED( field%sftlf )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%sftgif ) IF( ASSOCIATED( field%sftgif )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%cosmu0_rt ) IF( ASSOCIATED( field%cosmu0_rt )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%daylght_frc_rt ) IF( ASSOCIATED( field%daylght_frc_rt )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%albvisdir ) IF( ASSOCIATED( field%albvisdir )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%albnirdir ) IF( ASSOCIATED( field%albnirdir )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%albvisdif ) IF( ASSOCIATED( field%albvisdif )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%albnirdif ) IF( ASSOCIATED( field%albnirdif )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%emissivity ) IF( ASSOCIATED( field%emissivity )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%zf ) IF( ASSOCIATED( field%zf )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%zh ) IF( ASSOCIATED( field%zh )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%dz ) IF( ASSOCIATED( field%dz )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%presi_old ) IF( ASSOCIATED( field%presi_old )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED( field%presm_old )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED( field%ta )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%mdry ) IF( ASSOCIATED( field%mdry )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%mtrc ) IF( ASSOCIATED( field%mtrc ))
    !$ser verbatim   !NOACC UPDATE DEVICE( field%o3 ) IF( ASSOCIATED( field%o3 )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%acdnc ) IF( ASSOCIATED( field%acdnc )) 
    !$ser verbatim   !NOACC UPDATE DEVICE( field%aclc ) IF( ASSOCIATED( field%aclc )) 
#endif
#endif
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_rad_input

  SUBROUTINE serialize_rad_output(jg, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
#if defined(SERIALIZE_ECHAM_RAD) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .TRUE.

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_rad:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_rad:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !NOACC UPDATE HOST( field%o3 ) IF( ASSOCIATED( field%o3 )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%aclcov ) IF( ASSOCIATED( field%aclcov )) 
    !$ser verbatim   !NOACC UPDATE HOST( field%rldcs_rt ) IF( ASSOCIATED( field%rldcs_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rlucs_rt ) IF( ASSOCIATED( field%rlucs_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rsdcs_rt ) IF( ASSOCIATED( field%rsdcs_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rsucs_rt ) IF( ASSOCIATED( field%rsucs_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rld_rt ) IF( ASSOCIATED( field%rld_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rlu_rt ) IF( ASSOCIATED( field%rlu_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rsd_rt ) IF( ASSOCIATED( field%rsd_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rsu_rt ) IF( ASSOCIATED( field%rsu_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rvds_dir_rt ) IF( ASSOCIATED( field%rvds_dir_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rpds_dir_rt ) IF( ASSOCIATED( field%rpds_dir_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rnds_dir_rt ) IF( ASSOCIATED( field%rnds_dir_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rvds_dif_rt ) IF( ASSOCIATED( field%rvds_dif_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rpds_dif_rt ) IF( ASSOCIATED( field%rpds_dif_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rnds_dif_rt ) IF( ASSOCIATED( field%rnds_dif_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rvus_rt ) IF( ASSOCIATED( field%rvus_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rpus_rt ) IF( ASSOCIATED( field%rpus_rt ))
    !$ser verbatim   !NOACC UPDATE HOST( field%rnus_rt ) IF( ASSOCIATED( field%rnus_rt ))
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_rad-output jg=jg nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_rad_o3=field%o3(:,:,:) IF (ASSOCIATED(field%o3)) 
    !$ser data echam_rad_aclcov=field%aclcov(:,:) IF (ASSOCIATED(field%aclcov)) 
    !$ser data echam_rad_rldcs_rt=field%rldcs_rt(:,:,:) IF (ASSOCIATED(field%rldcs_rt))
    !$ser data echam_rad_rlucs_rt=field%rlucs_rt(:,:,:) IF (ASSOCIATED(field%rlucs_rt))
    !$ser data echam_rad_rsdcs_rt=field%rsdcs_rt(:,:,:) IF (ASSOCIATED(field%rsdcs_rt))
    !$ser data echam_rad_rsucs_rt=field%rsucs_rt(:,:,:) IF (ASSOCIATED(field%rsucs_rt))
    !$ser data echam_rad_rld_rt=field%rld_rt(:,:,:) IF (ASSOCIATED(field%rld_rt))
    !$ser data echam_rad_rlu_rt=field%rlu_rt(:,:,:) IF (ASSOCIATED(field%rlu_rt))
    !$ser data echam_rad_rsd_rt=field%rsd_rt(:,:,:) IF (ASSOCIATED(field%rsd_rt))
    !$ser data echam_rad_rsu_rt=field%rsu_rt(:,:,:) IF (ASSOCIATED(field%rsu_rt))
    !$ser data echam_rad_rvds_dir_rt=field%rvds_dir_rt(:,:) IF (ASSOCIATED(field%rvds_dir_rt))
    !$ser data echam_rad_rpds_dir_rt=field%rpds_dir_rt(:,:) IF (ASSOCIATED(field%rpds_dir_rt))
    !$ser data echam_rad_rnds_dir_rt=field%rnds_dir_rt(:,:) IF (ASSOCIATED(field%rnds_dir_rt))
    !$ser data echam_rad_rvds_dif_rt=field%rvds_dif_rt(:,:) IF (ASSOCIATED(field%rvds_dif_rt))
    !$ser data echam_rad_rpds_dif_rt=field%rpds_dif_rt(:,:) IF (ASSOCIATED(field%rpds_dif_rt))
    !$ser data echam_rad_rnds_dif_rt=field%rnds_dif_rt(:,:) IF (ASSOCIATED(field%rnds_dif_rt))
    !$ser data echam_rad_rvus_rt=field%rvus_rt(:,:) IF (ASSOCIATED(field%rvus_rt))
    !$ser data echam_rad_rpus_rt=field%rpus_rt(:,:) IF (ASSOCIATED(field%rpus_rt))
    !$ser data echam_rad_rnus_rt=field%rnus_rt(:,:) IF (ASSOCIATED(field%rnus_rt))
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lactive = .FALSE.
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_rad_output

END MODULE mo_ser_echam_rad
