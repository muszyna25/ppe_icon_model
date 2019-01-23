MODULE mo_echam_diagnostics
  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_statistics          ,ONLY: levels_horizontal_mean
  USE mo_name_list_output_init, ONLY: isRegistered
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field

  PUBLIC echam_global_diagnostics

CONTAINS
  SUBROUTINE echam_global_diagnostics(patch)
    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    REAL(wp)                           :: scr(nproma,patch%alloc_cell_blocks)

    REAL(wp) :: tas_gmean, rsdt_gmean, rsut_gmean, rlut_gmean, prec_gmean, evap_gmean, radtop_gmean, fwfoce_gmean

    ! global mean t2m, tas_gmean, if requested for output
    tas_gmean = 0.0_wp
    IF ( isRegistered("tas_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%tas(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & tas_gmean)
    END IF
    prm_field(patch%id)%tas_gmean = tas_gmean

    ! global mean toa incident shortwave radiation, rsdt
    rsdt_gmean = 0.0_wp
    IF ( isRegistered("rsdt_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rsdt(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rsdt_gmean)
    END IF
    prm_field(patch%id)%rsdt_gmean = rsdt_gmean

    ! global mean toa outgoing shortwave radiation, rsut
    rsut_gmean = 0.0_wp
    IF ( isRegistered("rsut_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rsut(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rsut_gmean)
    END IF
    prm_field(patch%id)%rsut_gmean = rsut_gmean

    ! global mean toa outgoing longwave radiation, rlut
    rlut_gmean = 0.0_wp
    IF ( isRegistered("rlut_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rlut(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rlut_gmean)
    END IF
    prm_field(patch%id)%rlut_gmean = rlut_gmean

    ! global mean precipitation flux, prec
    prec_gmean = 0.0_wp
    IF ( isRegistered("prec_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%pr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & prec_gmean)
    END IF
    prm_field(patch%id)%prec_gmean = prec_gmean

    ! global mean evaporation flux, evap
    evap_gmean = 0.0_wp
    IF ( isRegistered("evap_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%evap(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & evap_gmean)
    END IF
    prm_field(patch%id)%evap_gmean = evap_gmean

    ! global mean toa total radiation, radtop, derived variable
    radtop_gmean = 0.0_wp
    IF ( isRegistered("radtop_gmean") ) THEN
      scr(:,:) = 0.0_wp
      scr(:,:) = prm_field(patch%id)%rsdt(:,:) - prm_field(patch%id)%rsut(:,:) - prm_field(patch%id)%rlut(:,:)
      call levels_horizontal_mean( scr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & radtop_gmean)
    END IF
    prm_field(patch%id)%radtop_gmean = radtop_gmean

    ! global mean freshwater flux over ocean area, fwfoce, derived variable
    fwfoce_gmean = 0.0_wp
    IF ( isRegistered("fwfoce_gmean") ) THEN
      scr(:,:) = 0.0_wp
      scr(:,:) = (prm_field(patch%id)%pr(:,:) + prm_field(patch%id)%evap(:,:))*prm_field(patch%id)%sftof(:,:)
      call levels_horizontal_mean( scr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & fwfoce_gmean)
    END IF
    prm_field(patch%id)%fwfoce_gmean = fwfoce_gmean

    ! global mean ice cover fraction, icefrc - not set in atmosphere
 !  icefrc_gmean = 0.0_wp
 !  IF ( isRegistered("icefrc_gmean") ) THEN
 !    call levels_horizontal_mean( prm_field(patch%id)%icefrc(:,:), &
 !        & patch%cells%area(:,:), &
 !        & patch%cells%owned, &
 !        & icefrc_gmean)
 !  END IF
 !  prm_field(patch%id)%icefrc_gmean = icefrc_gmean

  END SUBROUTINE echam_global_diagnostics
END MODULE mo_echam_diagnostics
