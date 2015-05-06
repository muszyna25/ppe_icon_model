!>
!!  This module contains the interface routines for calling the RTTOV library
!!  for generating synthetic satellite images
!!
!! @par Revision History
!!  Developed by Guenther Zaengl, DWD (2015-04-20)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

MODULE mo_rttov_interface
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message_text, message, finish
USE mo_model_domain,        ONLY: t_patch, t_grid_cells, p_patch_local_parent, p_patch
USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state, p_int_state_local_parent
USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, p_grf_state, p_grf_state_local_parent
USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi, iqs, ltimer
USE mo_nwp_phy_state,       ONLY: prm_diag
USE mo_nonhydro_state,      ONLY: p_nh_state
USE mo_ext_data_state,      ONLY: ext_data
USE mo_ext_data_types,      ONLY: t_external_atmos
USE mo_nwp_lnd_state,       ONLY: p_lnd_state
USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
USE mo_impl_constants,      ONLY: min_rlcell_int
USE mo_loopindices,         ONLY: get_indices_c
USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c
USE mo_communication,       ONLY: exchange_data, exchange_data_mult
USE mo_nh_vert_interp,      ONLY: prepare_lin_intp, prepare_extrap, lin_intp, z_at_plevels
USE mo_satad,               ONLY: qsat_rho
USE mo_physical_constants,  ONLY: rd, vtmpc1, earth_radius
USE mo_cufunctions,         ONLY: foealfa
USE mo_math_constants,      ONLY: rad2deg
USE mo_timer,               ONLY: timer_start, timer_stop, timer_synsat
USE mo_synsat_config
#ifdef __USE_RTTOV
USE mo_rttov_ifc,           ONLY: rttov_fill_input, rttov_direct_ifc
#endif

IMPLICIT NONE

PRIVATE


PUBLIC :: copy_rttov_ubc, rttov_driver


CONTAINS


!>
!! Driver routine to call the RTTOV library for computing synthetic satellite images
!! For the time being, using a reduced radiation grid is mandatory
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2015-04-30
!!
SUBROUTINE rttov_driver (jg, jgp, nnow, nimg)

  INTEGER, INTENT(IN) :: jg, jgp ! grid ID and parent grid ID
  INTEGER, INTENT(IN) :: nnow    ! time level nnow valid for long time step
  INTEGER, INTENT(IN) :: nimg    ! number of synthetic satellite images contained in synsat_arr

  TYPE(t_grid_cells),     POINTER :: p_gcp
  TYPE(t_patch),          POINTER :: p_pp
  TYPE(t_nh_prog),        POINTER :: p_nh_prog
  TYPE(t_nh_diag),        POINTER :: p_nh_diag
  TYPE(t_nh_metrics),     POINTER :: p_nh_metrics
  TYPE(t_external_atmos), POINTER :: p_extdata
  TYPE(t_lnd_prog),       POINTER :: p_lnd_prog
  TYPE(t_lnd_diag),       POINTER :: p_lnd_diag

  REAL(wp), DIMENSION(nproma, nlev_rttov, p_patch_local_parent(jg)%nblks_c) :: &
    temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov

  REAL(wp) :: rg_synsat(nproma, nimg, p_patch_local_parent(jg)%nblks_c)

  REAL(wp), DIMENSION(nproma, p_patch_local_parent(jg)%nblks_c) :: &
    rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m

  INTEGER, DIMENSION(nproma, p_patch_local_parent(jg)%nblks_c) :: rg_stype, rg_wtype


  ! Local variables for RTTOV calls (with RTTOV-specific memory layout)
  REAL(wp), DIMENSION(nlev_rttov,nproma) :: temp, pres, qv
  REAL(wp), DIMENSION(6,nlev_rttov,nproma) :: clc, cld

  INTEGER,  DIMENSION(nproma*MAXVAL(numchans(:))) :: iprof, ichan
  REAL(wp), DIMENSION(MAXVAL(numchans(:)),nproma) :: emiss, t_b, t_b_clear, rad, rad_clear

  REAL(wp) :: pres_rttov(nlev_rttov), r_sat, sat_a(nproma), sat_z(nproma), alpha_e, r_atm, lon

  INTEGER :: idg(nproma), ish(nproma)

  INTEGER :: jb, jc, jk, i_startblk, i_endblk, is, ie, j, k
  INTEGER :: nlev_rg, isens, n_profs, ncalc, iprint, istatus

  IF (ltimer) CALL timer_start(timer_synsat)

  nlev_rg = p_patch(jgp)%nlev

  p_gcp        => p_patch_local_parent(jg)%cells
  p_pp         => p_patch_local_parent(jg)
  p_nh_prog    => p_nh_state(jg)%prog(nnow)
  p_nh_diag    => p_nh_state(jg)%diag
  p_nh_metrics => p_nh_state(jg)%metrics
  p_extdata    => ext_data(jg)%atm
  p_lnd_prog   => p_lnd_state(jg)%prog_lnd(nnow)
  p_lnd_diag   => p_lnd_state(jg)%diag_lnd


  ! distance of satellite from middle of the earth
  r_sat       = 35880.e3_wp + earth_radius



  CALL prepare_rttov_input(jg, jgp, nlev_rg, p_nh_metrics%z_ifc, p_nh_diag%pres,               &
    p_nh_diag%dpres_mc, p_nh_diag%temp, prm_diag(jg)%tot_cld, prm_diag(jg)%clc,                &
    p_nh_prog%tracer(:,:,:,iqs), prm_diag(jg)%con_udd(:,:,:,3), prm_diag(jg)%con_udd(:,:,:,7), &
    p_extdata%fr_land, p_extdata%fr_lake, p_lnd_diag%fr_seaice, prm_diag(jg)%cosmu0,           &
    p_nh_diag%pres_sfc, p_lnd_prog%t_g, prm_diag(jg)%t_2m, prm_diag(jg)%qv_2m,                 &
    prm_diag(jg)%u_10m, prm_diag(jg)%v_10m, prm_diag(jg)%buffer_rttov,                         &
    temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov,                  &
    rg_stype, rg_wtype, rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m)

  ! Define RTTOV levels
  CALL define_rttov_levels (nlev_rttov, pres_rttov)

  ! Choices for RTTOV ice cloud scheme and ice crystal shape
  idg(:) = iwc2effdiam ! McFarquhar et al. (2003); seems to the only numerically stable option
  ish(:) = iceshape    ! hexagonal crystals

  ! Call RTTOV library - hopefully it is threadsafe
  !
  i_startblk = p_gcp%start_block(grf_bdyintp_start_c)
  i_endblk   = p_gcp%end_block(min_rlcell_int)

#ifdef __USE_RTTOV
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jc,jk,j,k,pres,temp,qv,clc,cld,isens,n_profs,ncalc,iprof,ichan,emiss,lon, &
!$OMP           alpha_e,r_atm,sat_z,sat_a,t_b,t_b_clear,rad,rad_clear,iprint,istatus)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, is, ie, grf_bdyintp_start_c, min_rlcell_int)


    ! Copy input variables into RTTOV buffer
    DO jk = 1, nlev_rttov
      DO jc = is, ie
        pres(jk,jc) = pres_rttov(jk)
        temp(jk,jc) = temp_rttov(jc,jk,jb)
        qv(jk,jc)   = qv_rttov(jc,jk,jb)
        ! cld(1) = stratiform cloud water, computed as total - convective
        cld(1,jk,jc) = MAX(0._wp, qc_rttov(jc,jk,jb) - qcc_rttov(jc,jk,jb))
        ! cld(3) = convective cloud water
        cld(3,jk,jc) = qcc_rttov(jc,jk,jb)
        ! cld(6) = cloud ice
        cld(6,jk,jc) = qi_rttov(jc,jk,jb)
        ! clc(1) = cloud fraction
        clc(1,jk,jc) = clc_rttov(jc,jk,jb)
      ENDDO
    ENDDO


    istatus = rttov_fill_input(                            &
          press      = pres(:,is:ie) ,                     &
          temp       = temp(:,is:ie),                      &
          humi       = qv(:,is:ie),                        &
          t2m        = rg_t2m(is:ie,jb),                   &
          q2m        = rg_qv2m(is:ie,jb),                  &
          psurf      = rg_psfc(is:ie,jb),                  &
          hsurf      = rg_hsfc(is:ie,jb),                  &
          u10m       = rg_u10m(is:ie,jb),                  &
          v10m       = rg_v10m(is:ie,jb),                  &
          stemp      = rg_tsfc(is:ie,jb),                  &
          stype      = rg_stype(is:ie,jb),                 &
          watertype  = rg_wtype(is:ie,jb),                 &
          latgroundp = p_gcp%center(is:ie,jb)%lat*rad2deg, &
          satzenith  = (/(0.0_wp, jc=is,ie)/),             &
          sunZenith  = rg_cosmu0(is:ie,jb),                & ! actually unused for addsolar=.false.
          cloud      = cld(:,:,is:ie),                     &
          cfrac      = clc(:,:,is:ie),                     &
          idg        = idg(is:ie),                         &
          ish        = ish(is:ie),                         &
          addsolar   = .false.,                            &
          addrefrac  = .true.,                             &
          rttov9_compat = .false.,                         &
          addinterp  = .false.,                            &
          ivect      = 1)


    sensor_loop: DO isens = 1, num_sensors

      ! Set up instrument, channel and profile numbers for RTTOV

      n_profs = ie-is+1
      ncalc = n_profs * numchans(isens)
      DO j = 1, numchans(isens)
        iprof(j:ncalc:numchans(isens)) = (/ (k, k=1,n_profs) /)
        ichan(j:ncalc:numchans(isens)) = &
           (/ (sat_compute(isens)%nchan_list(j), k=1,n_profs) /)
      ENDDO


      ! Set/compute some sensor dependent quantities

      DO jc = is, ie 
        ! Since the emissitivy is intent(inout) in RTTOV, we have to 
        ! reinitialize it
        DO k = 1,  numchans(isens)
          emiss(k, jc-is+1) = sat_compute(isens)%emissivity(k)
        ENDDO
        lon = p_gcp%center(jc,jb)%lat - sat_compute(isens)%longitude
        ! Calculate the satellite zenith angle 
        alpha_e  = ACOS(COS(p_gcp%center(jc,jb)%lat) * COS(lon))
        r_atm    = SQRT(r_sat**2 + earth_radius**2 -2*r_sat*earth_radius*COS(alpha_e))
        sat_z(jc) = ASIN(SIN(alpha_e)*r_sat/r_atm) * rad2deg
        sat_z(jc) = MIN(ABS(sat_z(jc)), zenmax10)

        ! Calculate the satellite azimuth angle
        sat_a(jc) = rad2deg * (1. + ATAN2(TAN(lon), SIN(p_gcp%center(jc,jb)%lat)))
      ENDDO

      iprint = 0

      istatus = rttov_direct_ifc(                                 &
             isens,                                               &
             iprof(1:ncalc),                                      &
             ichan(1:ncalc),                                      &
             emiss(1:numchans(isens), 1:n_profs),                 &
             satAzim   = sat_a(is:ie),                            &
             satZenith = sat_z(is:ie),                            &
             T_b       = T_b(1:numchans(isens), 1:n_profs),       &
             T_b_clear = T_b_clear(1:numchans(isens), 1:n_profs), &
             rad       = rad      (1:numchans(isens), 1:n_profs), &
             radClear  = rad_clear(1:numchans(isens), 1:n_profs), &
             iprint    = iprint)

 !     IF (istatus /= NO_ERROR) THEN
 !       CALL print_profile(idims,j)
 !       WRITE(*,*) '*** RTTOV ERROR (rttov_direct_ifc, proc ', &
 !            my_cart_id,') ',TRIM(rttov_ifc_errMsg(istatus))
 !       WRITE(yerror,'(" ERROR *** while computing synthetic satellite &
 !            &images: rttov_fill_input ",I3)') istatus
 !       ierror = 9006
 !       RETURN
 !     ENDIF

 !     IF ((sat_compute(isens)%satellite == 'METEOSAT').and.&
 !         (sat_compute(isens)%nsat_id == 7)) THEN
 !       syn => synme7
 !     ELSEIF ((sat_compute(isens)%satellite == 'MSG').and.&
 !             (sat_compute(isens)%nsat_id == 2)) THEN
 !       syn => synmsg
 !     ELSE
 !       PRINT *,'*** ERROR in src_sat_tbs: satellite '//&
 !               TRIM(sat_compute(isens)%satellite), &
 !               sat_compute(isens)%nsat_id,' not implemented.'
 !       CYCLE
 !     ENDIF

      IF (sat_compute(isens)%lcloud_tem) THEN
        DO k = 1, numchans(isens)
          rg_synsat(is:ie,(k-1)*4+1, jb) = T_b(k, 1:n_profs) 
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_tem) THEN
        DO k = 1, numchans(isens)
          rg_synsat(is:ie, (k-1)*4+2, jb) = T_b_clear(k, 1:n_profs) 
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lcloud_rad) THEN
        DO k = 1, numchans(isens)
          rg_synsat(is:ie, (k-1)*4+3, jb) = Rad(k, 1:n_profs) 
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_rad) THEN
        DO k = 1, numchans(isens)
          rg_synsat(is:ie, (k-1)*4+4, jb) = Rad_clear(k, 1:n_profs) 
        ENDDO
      ENDIF

    END DO sensor_loop

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

  CALL downscale_rttov_output(jg, jgp, nimg, rg_synsat, prm_diag(jg)%synsat_arr)

  IF (ltimer) CALL timer_stop(timer_synsat)

END SUBROUTINE rttov_driver


!>
!! This routine interpolates the input fields for RTTOV to the reduced radiation grid,
!! combined with a vertical interpolation to the predefined RTTOV model levels
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2015-04-30
!!
SUBROUTINE prepare_rttov_input(jg, jgp, nlev_rg, z_ifc, pres, dpres, temp, tot_cld, clc, qs, &
  det_rate, con_upd, fr_land, fr_lake, fr_seaice, cosmu0, psfc, tsfc, t2m, qv2m, u10m, v10m, &
  buffer_rttov, temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov,    &
  rg_stype, rg_wtype, rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m)


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                                  &
    z_ifc(:,:,:), pres(:,:,:), dpres(:,:,:), temp(:,:,:), tot_cld(:,:,:,:), clc(:,:,:),    &
    qs(:,:,:), det_rate(:,:,:), con_upd(:,:,:), fr_land(:,:), fr_lake(:,:), fr_seaice(:,:),&
    cosmu0(:,:), tsfc(:,:), psfc(:,:), t2m(:,:), qv2m(:,:), u10m(:,:), v10m(:,:)

  ! Buffer for auxiliary temperature and pressure levels above the vertical nest interface
  REAL(wp), INTENT(IN) :: buffer_rttov(:,:,:)

  ! Corresponding output fields (on reduced grid)
  REAL(wp), INTENT(OUT) ::                                                                  &
    temp_rttov(:,:,:), qv_rttov(:,:,:), qc_rttov(:,:,:), qcc_rttov(:,:,:), qi_rttov(:,:,:), &
    qs_rttov(:,:,:), clc_rttov(:,:,:), rg_cosmu0(:,:), rg_tsfc(:,:), rg_psfc(:,:),          &
    rg_hsfc(:,:), rg_t2m(:,:), rg_qv2m(:,:), rg_u10m(:,:), rg_v10m(:,:)

  INTEGER, INTENT(OUT) :: rg_stype(:,:), rg_wtype(:,:)

  ! Local auxiliary fields
  REAL(wp), DIMENSION(nproma,nlev_rg,p_patch_local_parent(jg)%nblks_c)  ::   &
    rg_pres, rg_temp, rg_qv, rg_z_mc

  REAL(wp), DIMENSION(nproma,nlev_rg+1,p_patch_local_parent(jg)%nblks_c) ::  &
    rg_qc_vi, rg_qcc_vi, rg_qi_vi, rg_qs_vi, rg_clc_vi, rg_z_ifc

  REAL(wp), DIMENSION(nproma,nlev_rg) :: zdpres, zclc, zqc, zqcc, zqi, zqs

  REAL(wp), DIMENSION(nproma,nlev_rttov,p_patch_local_parent(jg)%nblks_c) ::  &
    wfac_lin, wfac_lin_i, z3d_rttov, pres3d_rttov, qc_vi_rttov, qcc_vi_rttov, &
    qi_vi_rttov, qs_vi_rttov, clc_vi_rttov

  REAL(wp), DIMENSION(nproma,p_patch_local_parent(jg)%nblks_c) :: wfacpbl1, wfacpbl2

  INTEGER,  DIMENSION(nproma,nlev_rttov,p_patch_local_parent(jg)%nblks_c) :: &
    idx0_lin, idx0_lin_i

  INTEGER,  DIMENSION(nproma,p_patch_local_parent(jg)%nblks_c) :: &
    kpbl1, kpbl2, bot_idx_lin, bot_idx_lin_i


  ! Convenience pointers
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_patch),      POINTER     :: p_pp

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  ! Indices
  INTEGER :: jb, jc, jk, jk1, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c, npromz_c

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg

  REAL(wp) :: qv_aux, pres_rttov(nlev_rttov), zfr_land, zfr_lake, zfr_si, zrho, &
              zdet_rate, zcon_upd

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'vertical interpolation and upscaling of RTTOV input',&
      jg,' =>',jgp
    CALL message('prepare_rttov_input',message_text)
  ENDIF

  ! Number of levels of the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter


  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf_state_local_parent(jg)%fbk_wgt_bln


  ! Define RTTOV levels
  CALL define_rttov_levels (nlev_rttov, pres_rttov)


  ! Upscaling of RTTOV input variables from compute grid to reduced grid
  i_startblk = p_gcp%start_block(grf_bdyintp_start_c)
  i_endblk   = p_gcp%end_block(min_rlcell_int)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,zdpres,zrho,qv_aux,zclc,zqc,zqcc,zqi,zqs, &
!$OMP            zfr_land,zfr_lake,zfr_si,zdet_rate,zcon_upd)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,          &
      i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

    DO jk = 1, nlev_rttov
      pres3d_rttov(i_startidx:i_endidx,jk,jb) = pres_rttov(jk)*100._wp  ! convert in Pa
    ENDDO

    DO jc = i_startidx, i_endidx
      rg_clc_vi(jc,nlevp1_rg,jb) = 0._wp
      rg_qc_vi(jc,nlevp1_rg,jb)  = 0._wp
      rg_qcc_vi(jc,nlevp1_rg,jb) = 0._wp
      rg_qi_vi(jc,nlevp1_rg,jb)  = 0._wp
      rg_qs_vi(jc,nlevp1_rg,jb)  = 0._wp

      rg_z_ifc(jc,nlevp1_rg,jb) =                                     &
        z_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        z_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        z_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        z_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
    ENDDO

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      DO jk = 1, nlev
        jk1 = jk + nshift
#else
    DO jk = 1, nlev
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
#endif

        rg_z_ifc(jc,jk1,jb) =                                       &
          z_ifc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          z_ifc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          z_ifc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          z_ifc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        rg_pres(jc,jk1,jb) =                                       &
          pres(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          pres(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          pres(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          pres(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        zdpres(jc,jk1) =                                            &
          dpres(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          dpres(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          dpres(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          dpres(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        rg_temp(jc,jk1,jb) =                                       &
          temp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          temp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          temp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          temp(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        qv_aux =                                                          &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),iqv)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),iqv)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),iqv)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),iqv)*p_fbkwgt(jc,jb,4)

        ! limit qv to saturation and convert into volume mixing ratio
        zrho = rg_pres(jc,jk1,jb)/(rd*(rg_temp(jc,jk1,jb)+vtmpc1*qv_aux))
        qv_aux = MIN(qv_aux, qsat_rho(rg_temp(jc,jk1,jb),zrho))
        rg_qv(jc,jk1,jb) = qv_aux / (1._wp - qv_aux) * 1.6077_wp

        zclc(jc,jk1) =                                            &
          clc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          clc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          clc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          clc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        zqc(jc,jk1) =                                                     &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),iqc)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),iqc)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),iqc)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),iqc)*p_fbkwgt(jc,jb,4)

        ! convective cloud water (computation taken over from cloud cover scheme)
        !
        ! averaged convective updraft detrainment rate
        zdet_rate =                                                    &
          det_rate(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          det_rate(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          det_rate(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          det_rate(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        ! averaged updraft condensate
        zcon_upd =                                                    &
          con_upd(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          con_upd(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          con_upd(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          con_upd(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        ! diagnosed averaged convective cloud water, must not exceed total QC
        zqcc(jc,jk1) = MIN((zdet_rate/zrho) / (zdet_rate/zrho + 1.0_wp/1500._wp) * &
                       zcon_upd * foealfa(rg_temp(jc,jk1,jb)), zqc(jc,jk1) )

        zqi(jc,jk1) =                                                     &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),iqi)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),iqi)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),iqi)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),iqi)*p_fbkwgt(jc,jb,4)

        zqs(jc,jk1) =                                            &
          qs(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          qs(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          qs(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          qs(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO


    ! Fill levels above model top
    DO jk = 1, nshift
      DO jc = i_startidx, i_endidx
        rg_z_ifc(jc,jk,jb) = buffer_rttov(jc,jk,jb)
        rg_pres (jc,jk,jb) = buffer_rttov(jc,nshift+jk,jb)
        zdpres(jc,jk)      = buffer_rttov(jc,2*nshift+jk,jb)
        rg_temp (jc,jk,jb) = buffer_rttov(jc,3*nshift+jk,jb)
        qv_aux             = buffer_rttov(jc,4*nshift+jk,jb)
        !
        ! convert qv into volume mixing ratio
        rg_qv(jc,jk,jb)    = qv_aux / (1._wp - qv_aux) * 1.6077_wp
        !
        ! for the cloud variables, we assume that the top of the nested model domain is not
        ! lower than the upper boundary of moisture physics calculation
        zclc(jc,jk) = 0._wp
        zqc(jc,jk)  = 0._wp
        zqcc(jc,jk) = 0._wp
        zqi(jc,jk)  = 0._wp
        zqs(jc,jk)  = 0._wp
      ENDDO
    ENDDO

    ! Compute rg_z_mc
    DO jk = 1, nlev_rg
      DO jc = i_startidx, i_endidx
        rg_z_mc(jc,jk,jb) = 0.5_wp*(rg_z_ifc(jc,jk,jb)+rg_z_ifc(jc,jk+1,jb))
      ENDDO
    ENDDO

    ! Vertically integrate the cloud quantities
    DO jk = nlev_rg, 1, -1
      DO jc = i_startidx, i_endidx
        rg_clc_vi(jc,jk,jb) = rg_clc_vi(jc,jk+1,jb) + zdpres(jc,jk)*zclc(jc,jk)
        rg_qc_vi(jc,jk,jb)  = rg_qc_vi(jc,jk+1,jb)  + zdpres(jc,jk)*zqc(jc,jk)
        rg_qcc_vi(jc,jk,jb) = rg_qcc_vi(jc,jk+1,jb) + zdpres(jc,jk)*zqcc(jc,jk)
        rg_qi_vi(jc,jk,jb)  = rg_qi_vi(jc,jk+1,jb)  + zdpres(jc,jk)*zqi(jc,jk)
        rg_qs_vi(jc,jk,jb)  = rg_qs_vi(jc,jk+1,jb)  + zdpres(jc,jk)*zqs(jc,jk)
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx

      rg_hsfc(jc,jb) = rg_z_ifc(jc,nlevp1_rg,jb) * 0.001_wp ! convert in km

      zfr_land =                                                 &
        fr_land(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_land(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_land(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_land(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      zfr_lake =                                                 &
        fr_lake(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_lake(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_lake(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_lake(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      zfr_si =                                                     &
        fr_seaice(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_seaice(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_seaice(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_seaice(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      IF (zfr_land > 0.5_wp) THEN
        rg_stype(jc,jb) = 0
        rg_wtype(jc,jb) = 0
      ELSE IF (zfr_si > 0.5_wp) THEN
        rg_stype(jc,jb) = 2
        rg_wtype(jc,jb) = 1
      ELSE 
        rg_stype(jc,jb) = 1
        IF (zfr_lake > 0.5_wp) THEN
          rg_wtype(jc,jb) = 0
        ELSE
          rg_wtype(jc,jb) = 1
        ENDIF
      ENDIF
      
      rg_cosmu0(jc,jb) =                                        &
        cosmu0(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        cosmu0(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        cosmu0(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        cosmu0(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      rg_tsfc(jc,jb) =                                        &
        tsfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        tsfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        tsfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        tsfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      rg_psfc(jc,jb) =  MIN(1099.9_wp, 0.01_wp* (             & ! convert in hPa
        psfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        psfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        psfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        psfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4) ) )

      rg_t2m(jc,jb) =                                        &
        t2m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        t2m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        t2m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        t2m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      qv_aux =                                                &
        qv2m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        qv2m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        qv2m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        qv2m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ! convert qv into volume mixing ratio
      rg_qv2m(jc,jb) = qv_aux / (1._wp - qv_aux) * 1.6077_wp

      rg_u10m(jc,jb) =                                        &
        u10m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        u10m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        u10m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        u10m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      rg_v10m(jc,jb) =                                        &
        v10m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        v10m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        v10m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        v10m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL



  ! Compute interpolation coefficients

  nblks_c  = p_gcp%end_block(min_rlcell_int)
  npromz_c = p_gcp%end_index(min_rlcell_int)

  CALL prepare_extrap(rg_z_mc, nblks_c, npromz_c, nlev_rg,    &
                      kpbl1, wfacpbl1, kpbl2, wfacpbl2     )

  CALL z_at_plevels(rg_pres, rg_temp, rg_z_mc, pres3d_rttov, z3d_rttov, &
                    nblks_c, npromz_c, nlev_rg, nlev_rttov,             &
                    kpbl1, wfacpbl1, kpbl2, wfacpbl2                    )

  CALL prepare_lin_intp(rg_z_mc, z3d_rttov, nblks_c, npromz_c, nlev_rg, &
                        nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin)                       

  CALL prepare_lin_intp(rg_z_ifc, z3d_rttov, nblks_c, npromz_c, nlevp1_rg, &
                        nlev_rttov, wfac_lin_i, idx0_lin_i, bot_idx_lin_i)                       

  ! Execute interpolation

  CALL lin_intp(rg_temp, temp_rttov,  nblks_c, npromz_c, nlev_rg,   &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,        &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.FALSE., &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.               )

  CALL lin_intp(rg_qv, qv_rttov,  nblks_c, npromz_c, nlev_rg,       &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,        &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,  &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.               )

  CALL lin_intp(rg_qc_vi, qc_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,         &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,   &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.                )

  CALL lin_intp(rg_qcc_vi, qcc_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,           &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,     &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.                  )

  CALL lin_intp(rg_qi_vi, qi_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,         &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,   &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.                )

  CALL lin_intp(rg_qs_vi, qs_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,         &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,   &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.                )

  CALL lin_intp(rg_clc_vi, clc_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,           &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,     &
                l_pd_limit=.TRUE., l_extrapol=.FALSE.                  )


  ! Calculate layer averages of cloud variables from vertically integrated fields

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,          &
      i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

    DO jk = 1, nlev_rttov-1
      DO jc = i_startidx, i_endidx
        qc_rttov(jc,jk,jb)  = (qc_vi_rttov(jc,jk,jb)-qc_vi_rttov(jc,jk+1,jb)) / &
                              (100._wp*(pres_rttov(jk+1)-pres_rttov(jk)))
        IF (qc_rttov(jc,jk,jb) < 1.e-8_wp) qc_rttov(jc,jk,jb) = 0._wp
        qcc_rttov(jc,jk,jb) = (qcc_vi_rttov(jc,jk,jb)-qcc_vi_rttov(jc,jk+1,jb)) / &
                              (100._wp*(pres_rttov(jk+1)-pres_rttov(jk)))
        IF (qcc_rttov(jc,jk,jb) < 1.e-8_wp) qcc_rttov(jc,jk,jb) = 0._wp
        qi_rttov(jc,jk,jb)  = (qi_vi_rttov(jc,jk,jb)-qi_vi_rttov(jc,jk+1,jb)) / &
                              (100._wp*(pres_rttov(jk+1)-pres_rttov(jk)))
        IF (qi_rttov(jc,jk,jb) < 1.e-8_wp) qi_rttov(jc,jk,jb) = 0._wp
        qs_rttov(jc,jk,jb)  = (qs_vi_rttov(jc,jk,jb)-qs_vi_rttov(jc,jk+1,jb)) / &
                              (100._wp*(pres_rttov(jk+1)-pres_rttov(jk)))
        IF (qs_rttov(jc,jk,jb) < 1.e-8_wp) qs_rttov(jc,jk,jb) = 0._wp
        clc_rttov(jc,jk,jb) = MIN(1._wp,(clc_vi_rttov(jc,jk,jb)-clc_vi_rttov(jc,jk+1,jb)) / &
                              (100._wp*(pres_rttov(jk+1)-pres_rttov(jk))) )
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE prepare_rttov_input



!>
!! Back-interpolation of synthetic satellite images from reduced grid to compute grid
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2015-04-30
!!
SUBROUTINE downscale_rttov_output(jg, jgp, nimg, rg_satimg, satimg)

  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nimg     ! number of sat images (= number of "levels" in satimg field)

  ! Array with sat images on reduced grid
  REAL(wp), INTENT(INOUT) :: rg_satimg(:,:,:)

  ! Array with sat images on full grid
  REAL(wp), INTENT(OUT) :: satimg(:,:,:)

  ! Indices
  INTEGER :: i_chidx, jb, jk, jc, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4
  INTEGER :: rl_start, rl_end

  ! Convenience pointers
  TYPE(t_patch),         POINTER :: p_pp
  TYPE(t_int_state),     POINTER :: p_int
  TYPE(t_gridref_state), POINTER :: p_grf

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

  p_pp  => p_patch_local_parent(jg)
  p_int => p_int_state_local_parent(jg)
  p_grf => p_grf_state_local_parent(jg)

  ! pointers to child index/block
  iidx => p_pp%cells%child_idx
  iblk => p_pp%cells%child_blk

  i_chidx  = p_patch(jg)%parent_child_index

  ! Synchronize sat image array before interpolation
  CALL exchange_data(p_pp%comm_pat_c, rg_satimg)

  ! Execute interpolation from reduced grid to full grid
  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 1, 1, &
    &                         rg_satimg, satimg, overshoot_fac=1.0_wp)

  ! Fill nest boundary points by copying the values from the reduced grid to the full grid
  rl_start = -1
  rl_end   = -2
  i_startblk = p_pp%cells%start_block(rl_start)
  i_endblk   = p_pp%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jc1,jc2,jc3,jc4,jb1,jb2,jb3,jb4)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,       &
                       i_startidx, i_endidx, rl_start, rl_end)

 
    DO jc = i_startidx, i_endidx
      jc1 = iidx(jc,jb,1)
      jc2 = iidx(jc,jb,2)
      jc3 = iidx(jc,jb,3)
      jc4 = iidx(jc,jb,4)
      jb1 = iblk(jc,jb,1)
      jb2 = iblk(jc,jb,2)
      jb3 = iblk(jc,jb,3)
      jb4 = iblk(jc,jb,4)

      DO jk = 1,nimg
        satimg(jc1,jk,jb1) = rg_satimg(jc,jk,jb)
        satimg(jc2,jk,jb2) = rg_satimg(jc,jk,jb)
        satimg(jc3,jk,jb3) = rg_satimg(jc,jk,jb)
        satimg(jc4,jk,jb4) = rg_satimg(jc,jk,jb)
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE downscale_rttov_output


!>
!! Define the pressure levels (interface/half levels, to be precise) used by RTTOV
!! Unit is hPa!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2015-04-20
!!
SUBROUTINE define_rttov_levels (nlev_rt, pres_rttov)

  INTEGER, INTENT(IN) :: nlev_rt
  REAL(wp), INTENT(OUT) :: pres_rttov(:)


      SELECT CASE (nlev_rt)
      CASE (43) ! RTTOV9 levels and RTTOV10 levels (in RTTOV-9 compat mode,
                ! using RTTOV9 coefficients)
        pres_rttov =                                                        &
         (/           0.10,   0.29,   0.69,   1.42,  2.611,  4.407,   6.95, &
             10.37,  14.81,  20.40,  27.26,  35.51,  45.29,  56.73,  69.97, &
             85.18, 102.05, 122.04, 143.84, 167.95, 194.36, 222.94, 253.71, &
            286.60, 321.50, 358.28, 396.81, 436.95, 478.54, 521.46, 565.54, &
            610.60, 656.43, 702.73, 749.12, 795.09, 839.95, 882.80, 922.46, &
            957.44, 985.88,1005.43,1013.25  /)
      CASE (44) ! RTTOV10 levels (using RTTOV9 coefficients)
        pres_rttov =                                                        &
         (/  0.005,   0.10,   0.29,   0.69,   1.42,  2.611,  4.407,   6.95, &
             10.37,  14.81,  20.40,  27.26,  35.51,  45.29,  56.73,  69.97, &
             85.18, 102.05, 122.04, 143.84, 167.95, 194.36, 222.94, 253.71, &
            286.60, 321.50, 358.28, 396.81, 436.95, 478.54, 521.46, 565.54, &
            610.60, 656.43, 702.73, 749.12, 795.09, 839.95, 882.80, 922.46, &
            957.44, 985.88,1005.43,1013.25  /)
      CASE (50) ! RTTOV-10 levels (in RTTOV-9 compat mode)
        pres_rttov =                                                         &
         (/             0.01,    0.10,    0.20,     0.50,    0.80,    1.20,  &
               1.60,    2.20,    2.70,    3.50,     4.20,    5.00,    6.95,  &
              10.37,   14.81,   20.40,   27.26,    35.51,   45.29,   56.73,  &
              69.97,   85.18,  102.05,  122.04,   143.84,  167.95,  194.36,  &
             222.94,  253.71,  286.60,  321.50,   358.38,  396.81,  436.95,  &
             478.54,  521.46,  565.54,  610.60,   656.43,  702.73,  749.12,  &
             795.09,  839.95,  882.80,  922.46,   957.44,  985.88, 1005.43,  &
            1025.00, 1050.00 /)
      CASE (51) ! RTTOV-10 levels
        pres_rttov =                                                         &
         (/    0.005,   0.01,    0.10,    0.20,     0.50,    0.80,    1.20,  &
               1.60,    2.20,    2.70,    3.50,     4.20,    5.00,    6.95,  &
              10.37,   14.81,   20.40,   27.26,    35.51,   45.29,   56.73,  &
              69.97,   85.18,  102.05,  122.04,   143.84,  167.95,  194.36,  &
             222.94,  253.71,  286.60,  321.50,   358.38,  396.81,  436.95,  &
             478.54,  521.46,  565.54,  610.60,   656.43,  702.73,  749.12,  &
             795.09,  839.95,  882.80,  922.46,   957.44,  985.88, 1005.43,  &
            1025.00, 1050.00 /)
      CASE (53) ! RTTOV-10 levels (in RTTOV-9 compat mode)
        pres_rttov =                                                         &
         (/             0.0131,    0.0304,    0.0644,    0.1263,    0.2324,  &
             0.4052,    0.6749,    1.0801,    1.6691,    2.5011,    3.6462,  &
             5.1864,    7.2150,    9.8368,   13.1672,   17.3308,   22.4601,  &
            28.6937,   36.1735,   45.0430,   55.4433,   67.5109,   81.3744,  &
            97.1505,  114.9420,  134.8320,  156.8850,  181.1390,  207.6090,  &
           236.2780,  267.1010,  300.0000,  334.8650,  371.5530,  409.8890,  &
           449.6680,  490.6520,  532.5770,  575.1540,  618.0710,  660.9960,  &
           703.5860,  745.4840,  786.3280,  825.7550,  863.4050,  898.9280,  &
           931.9850,  962.2590,  989.4510, 1013.2900, 1033.5400, 1050.0000  /)
      CASE (54) ! RTTOV-10 levels
        pres_rttov =                                                         &
         (/  0.0050,    0.0131,    0.0304,    0.0644,    0.1263,    0.2324,  &
             0.4052,    0.6749,    1.0801,    1.6691,    2.5011,    3.6462,  &
             5.1864,    7.2150,    9.8368,   13.1672,   17.3308,   22.4601,  &
            28.6937,   36.1735,   45.0430,   55.4433,   67.5109,   81.3744,  &
            97.1505,  114.9420,  134.8320,  156.8850,  181.1390,  207.6090,  &
           236.2780,  267.1010,  300.0000,  334.8650,  371.5530,  409.8890,  &
           449.6680,  490.6520,  532.5770,  575.1540,  618.0710,  660.9960,  &
           703.5860,  745.4840,  786.3280,  825.7550,  863.4050,  898.9280,  &
           931.9850,  962.2590,  989.4510, 1013.2900, 1033.5400, 1050.0000  /)
      ! case (100) hires coefficients. only available for iasi, airs, cris
      ! case (101)
      CASE default
        call finish('define_rttov_levels', 'failed to determine pressure levels')
      END SELECT

END SUBROUTINE define_rttov_levels


!>
!! This routine copies additional model levels to the local parent grid if the RTTOV library
!! is called on a vertically nested grid
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2015-04-26
!!
SUBROUTINE copy_rttov_ubc (jg, jgc)

  ! Input grid parameters
  INTEGER, INTENT(in) :: jg, jgc

  ! Local fields

  INTEGER :: nshift

  nshift = p_patch(jgc)%nshift
  CALL exchange_data_mult(p_patch_local_parent(jgc)%comm_pat_glb_to_loc_c, 5, 5*nshift,                            &
       RECV1=prm_diag(jgc)%buffer_rttov(:,         1:  nshift,:), SEND1=p_nh_state(jg)%metrics%z_ifc(:,1:nshift,:),&
       RECV2=prm_diag(jgc)%buffer_rttov(:,  nshift+1:2*nshift,:), SEND2=p_nh_state(jg)%diag%pres(:,1:nshift,:),    &
       RECV3=prm_diag(jgc)%buffer_rttov(:,2*nshift+1:3*nshift,:), SEND3=p_nh_state(jg)%diag%dpres_mc(:,1:nshift,:),&
       RECV4=prm_diag(jgc)%buffer_rttov(:,3*nshift+1:4*nshift,:), SEND4=p_nh_state(jg)%diag%temp(:,1:nshift,:),    &
       RECV5=prm_diag(jgc)%buffer_rttov(:,4*nshift+1:5*nshift,:), SEND5=prm_diag(jg)%tot_cld(:,1:nshift,:,iqv)     )
END SUBROUTINE copy_rttov_ubc


END MODULE mo_rttov_interface



