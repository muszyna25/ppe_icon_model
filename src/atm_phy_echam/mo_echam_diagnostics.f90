MODULE mo_echam_diagnostics
  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_statistics          ,ONLY: levels_horizontal_mean
  USE mo_name_list_output_init, ONLY: isRegistered
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field

  PUBLIC echam_global_diagnostics

CONTAINS
  SUBROUTINE echam_global_diagnostics(patch)
    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch

    REAL(wp) :: t2m_global

    ! global mean t2m if requested for output
    t2m_global = 0.0_wp
    IF ( isRegistered("t2m_global") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%tas(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & t2m_global)
    END IF
    prm_field(patch%id)%t2m_global = t2m_global
  END SUBROUTINE echam_global_diagnostics
END MODULE mo_echam_diagnostics
