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

    REAL(wp) :: tas_gmean

    ! global mean t2m if requested for output
    tas_gmean = 0.0_wp
    IF ( isRegistered("tas_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%tas(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & tas_gmean)
    END IF
    prm_field(patch%id)%tas_gmean = tas_gmean
  END SUBROUTINE echam_global_diagnostics
END MODULE mo_echam_diagnostics
