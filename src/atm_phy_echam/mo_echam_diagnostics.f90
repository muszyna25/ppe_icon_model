MODULE mo_echam_diagnostics
  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_statistics          ,ONLY: levels_horizontal_mean
  USE mo_name_list_output_init, ONLY: isRegistered
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field
  !$ser verbatim USE mo_ser_echam_global_diag, ONLY: serialize_input,&
  !$ser verbatim                                     serialize_output

  PUBLIC echam_global_diagnostics

CONTAINS
  SUBROUTINE echam_global_diagnostics(patch)
    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    REAL(wp)                           :: scr(nproma,patch%alloc_cell_blocks)

    REAL(wp) :: tas_gmean, rsdt_gmean, rsut_gmean, rlut_gmean, prec_gmean, evap_gmean, radtop_gmean, fwfoce_gmean
    TYPE(t_echam_phy_field), POINTER    :: field
    INTEGER  :: jc, jk

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_input(patch%id, prm_field(patch%id))

    ! global mean t2m, tas_gmean, if requested for output
    tas_gmean = 0.0_wp
    IF ( isRegistered("tas_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%tas(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & tas_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%tas_gmean = tas_gmean

    ! global mean toa incident shortwave radiation, rsdt
    rsdt_gmean = 0.0_wp
    IF ( isRegistered("rsdt_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rsdt(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rsdt_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%rsdt_gmean = rsdt_gmean

    ! global mean toa outgoing shortwave radiation, rsut
    rsut_gmean = 0.0_wp
    IF ( isRegistered("rsut_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rsut(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rsut_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%rsut_gmean = rsut_gmean

    ! global mean toa outgoing longwave radiation, rlut
    rlut_gmean = 0.0_wp
    IF ( isRegistered("rlut_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rlut(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rlut_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%rlut_gmean = rlut_gmean

    ! global mean precipitation flux, prec
    prec_gmean = 0.0_wp
    IF ( isRegistered("prec_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%pr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & prec_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%prec_gmean = prec_gmean

    ! global mean evaporation flux, evap
    evap_gmean = 0.0_wp
    IF ( isRegistered("evap_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%evap(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & evap_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%evap_gmean = evap_gmean

    ! global mean toa total radiation, radtop, derived variable
    radtop_gmean = 0.0_wp
    IF ( isRegistered("radtop_gmean") ) THEN

      field => prm_field(patch%id)

      !$ACC DATA PRESENT( field%rsdt, field%rsut, field%rlut ) &
      !$ACC       CREATE( scr )

      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP SEQ
      DO jk = 1, patch%alloc_cell_blocks
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nproma
          scr(jc,jk) = 0.0_wp
          scr(jc,jk) = field%rsdt(jc,jk) - field%rsut(jc,jk) - field%rlut(jc,jk)
        END DO
      END DO
      !$ACC END PARALLEL
      call levels_horizontal_mean( scr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & radtop_gmean, lopenacc=.TRUE.)

      !$ACC END DATA

      NULLIFY(field)
    END IF
    prm_field(patch%id)%radtop_gmean = radtop_gmean

    ! global mean freshwater flux over ocean area, fwfoce, derived variable
    fwfoce_gmean = 0.0_wp
    IF ( isRegistered("fwfoce_gmean") ) THEN

      field => prm_field(patch%id)

      !$ACC DATA PRESENT( field%pr, field%evap, field%sftof ) &
      !$ACC       CREATE( scr )

      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP SEQ
      DO jk = 1, patch%alloc_cell_blocks
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nproma
          scr(jc,jk) = 0.0_wp
          scr(jc,jk) = (field%pr(jc,jk) + field%evap(jc,jk))*field%sftof(jc,jk)
        END DO
      END DO
      !$ACC END PARALLEL
      call levels_horizontal_mean( scr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & fwfoce_gmean, lopenacc=.TRUE.)

      !$ACC END DATA

      NULLIFY(field)
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

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_output(patch%id, prm_field(patch%id))

  END SUBROUTINE echam_global_diagnostics
END MODULE mo_echam_diagnostics
