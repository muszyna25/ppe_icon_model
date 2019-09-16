!! @brief Subroutine interface_echam_vdf calls the vertical diffusion and the surface schemes.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_vdf

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_exception           ,ONLY: finish, warning

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_vdf, &
    &                               timer_vdf_dn, timer_vdf_sf, timer_vdf_up

  USE mo_echam_rad_config    ,ONLY: echam_rad_config
  USE mo_ccycle_config       ,ONLY: ccycle_config
  USE mo_physical_constants  ,ONLY: amco2, amd
  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqt, ico2
  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep  ,ONLY: vdiff_up
  USE mo_vdiff_solver        ,ONLY: nvar_vdiff, nmatrix, imh, imqv, ih_vdiff=>ih, iqv_vdiff=>iqv
  
  USE mo_echam_sfc_indices   ,ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface             ,ONLY: update_surface
  USE mo_surface_diag        ,ONLY: nsurf_diag
  USE mo_run_config          ,ONLY: lart
  !$ser verbatim USE mo_ser_echam_vdf, ONLY: serialize_vdf_input, serialize_vdf_output,&
  !$ser verbatim                             serialize_vdf_vd_input, serialize_vdf_vd_output,&
  !$ser verbatim                             serialize_vdf_us_input, serialize_vdf_us_output,&
  !$ser verbatim                             serialize_vdf_vu_input, serialize_vdf_vu_output,&
  !$ser verbatim                             serialize_vdf_nd_input, serialize_vdf_nd_output,&
  !$ser verbatim                             serialize_vdf_chk_A_output,&
  !$ser verbatim                             serialize_vdf_chk_B_output,&
  !$ser verbatim                             serialize_vdf_chk_C_output,&
  !$ser verbatim                             serialize_vdf_chk_D_output,&
  !$ser verbatim                             serialize_vdf_chk_E_output,&
  !$ser verbatim                             serialize_vdf_chk_F_output,&
  !$ser verbatim                             serialize_vdf_chk_G_output,&
  !$ser verbatim                             serialize_vdf_chk_H_output
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_vdf

CONTAINS

  SUBROUTINE interface_echam_vdf(jg, jb, jcs, jce     ,&
       &                         nproma,nlev,ntracer  ,& 
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         pdtime               )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev,ntracer
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Pointers
    !
    LOGICAL                 ,POINTER    :: lparamcpl
    LOGICAL                 ,POINTER    :: lmig
    INTEGER                 ,POINTER    :: fc_vdf
    REAL(wp)                ,POINTER    :: vmr_co2
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER  :: jl, jk, jsfc, jt
    INTEGER  :: nlevm1, nlevp1
    INTEGER  :: ntrac
    !
    REAL(wp) :: zxt_emis(nproma,ntracer-iqt+1)   !< tracer tendency due to surface emission

    REAL(wp) :: zcpt_sfc_tile(nproma,nsfc_type)  !< dry static energy at surface
   
    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)
    REAL(wp) :: ri_tile(nproma,nsfc_type)           !< Richardson number
    REAL(wp) :: zaa    (nproma,nlev,3,nmatrix)      !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nproma,3,nsfc_type,imh:imqv)!< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nproma,nlev,nvar_vdiff)             !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nproma,nsfc_type,ih_vdiff:iqv_vdiff)!< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF, JSBACH
    REAL(wp) :: zfactor_sfc(nproma)
    REAL(wp) :: zco2       (nproma)          !< co2 value passed on to jsbach

    REAL(wp) :: zqx      (nproma,nlev)       !< total cloud condensate
    REAL(wp) :: zcptgz   (nproma,nlev)       !< dry static energy
    REAL(wp) :: zthvvar  (nproma,nlev)       !< intermediate value of thvvar
    REAL(wp) :: dummy    (nproma,nlev)       !< to replace thvvar
    REAL(wp) :: dummyx   (nproma,nlev)       !< to replace xvar
    REAL(wp) :: ztottevn (nproma,nlev)       !< intermediate value of TTE
    REAL(wp) :: zch_tile (nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zchn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zcdn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zcfnc_tile(nproma,nsfc_type) !< for "nsurf_diag"
    REAL(wp) :: zbn_tile (nproma,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbhn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbm_tile (nproma,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbh_tile (nproma,nsfc_type)  !< for "nsurf_diag"

    ! output variables of vdiff_down
    !
    REAL(wp) :: wstar      (nproma)
    REAL(wp) :: qs_sfc_tile(nproma,nsfc_type)
    REAL(wp) :: hdtcbl     (nproma)
    REAL(wp) :: ri_atm     (nproma,nlev)
    REAL(wp) :: mixlen     (nproma,nlev)
    REAL(wp) :: cfm        (nproma,nlev)
    REAL(wp) :: cfm_tile   (nproma,nsfc_type)
    REAL(wp) :: cfh        (nproma,nlev)
    REAL(wp) :: cfh_tile   (nproma,nsfc_type)
    REAL(wp) :: cfv        (nproma,nlev)
    REAL(wp) :: cftotte    (nproma,nlev)
    REAL(wp) :: cfthv      (nproma,nlev)

    ! output variables of update_surface
    !
    REAL(wp) :: q_snocpymlt     (nproma)
    REAL(wp) :: tend_ta_sfc     (nproma)
    REAL(wp) :: q_rlw_impl      (nproma)
    REAL(wp) :: tend_ta_rlw_impl(nproma)

    ! output variables of vdiff_up
    !
    REAL(wp) :: kedisp       (nproma)
    REAL(wp) :: q_vdf        (nproma,nlev)
    REAL(wp) :: tend_ta_vdf  (nproma,nlev)
    REAL(wp) :: tend_ua_vdf  (nproma,nlev)
    REAL(wp) :: tend_va_vdf  (nproma,nlev)
    REAL(wp), TARGET :: tend_qtrc_vdf(nproma,nlev,ntracer)
    REAL(wp), POINTER, CONTIGUOUS :: tend_qtrc_vdf_iqt(:,:,:)
    REAL(wp), TARGET :: tend_qtrc_vdf_dummy(nproma,nlev,0)
    REAL(wp) :: mmr_co2

    IF (ltimer) CALL timer_start(timer_vdf)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    lmig      => echam_phy_config(jg)%lmig
    fc_vdf    => echam_phy_config(jg)%fc_vdf
    vmr_co2   => echam_rad_config(jg)%vmr_co2
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_vdf_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    nlevm1 = nlev-1
    nlevp1 = nlev+1
    ntrac  = ntracer-iqt+1  ! number of tracers excluding water vapour and hydrometeors

    !$ACC DATA PCREATE( zxt_emis ) IF( ntrac > 0 )
    !$ACC DATA PRESENT( field ),                                                &
    !$ACC      PCREATE( zcpt_sfc_tile, ri_tile, zqx, zcptgz, zbn_tile,          &
    !$ACC               zbhn_tile, zbm_tile, zbh_tile, dummy, dummyx,           &
    !$ACC               wstar, qs_sfc_tile, hdtcbl, ri_atm, mixlen, cfm,        &
    !$ACC               cfm_tile, cfh, cfh_tile, cfv, cftotte, cfthv, zaa,      &
    !$ACC               zaa_btm, zbb, zbb_btm, zfactor_sfc, zcpt_sfc_tile,      &
    !$ACC               zthvvar, ztottevn, zch_tile, kedisp, tend_ua_vdf,       &
    !$ACC               tend_va_vdf, q_vdf, tend_qtrc_vdf, q_snocpymlt, zco2,   &
    !$ACC               tend_ta_sfc, q_rlw_impl, tend_ta_rlw_impl, tend_ta_vdf )

    !$ser verbatim zaa = 0._wp
    !$ser verbatim zbb = 0._wp
    !$ser verbatim !$ACC UPDATE DEVICE( zaa, zbb )

!!$    ! Emission of aerosols or other tracers (not implemented yet)
!!$    IF (ntrac>0) THEN
!!$       CALL tracer_emission()
!!$    ENDIF
    !
    ! default setting for all tracers
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( ntrac > 0 )
    !$ACC LOOP SEQ
    DO jt = 1,ntrac
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,jce
        zxt_emis(jl,jt) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL
    !

!    ! set emissions of co2, if any (hardcoded, co2-tracer 5 eq. zxt_emis(:,2))
!    !  should never be set for lmig or lart set to true

!    IF (ntrac>=2 .AND. .NOT. lmig .AND. .NOT. lart) THEN
!      IF (ccycle_config(jg)%iccy_co2conc .EQ. 1 .AND. ccycle_config(jg)%iccy_co2flux .EQ. 2) THEN
!        !$ACC DATA PRESENT( field%fco2nat )
!        !$ACC PARALLEL DEFAULT(PRESENT)
!        !$ACC LOOP GANG VECTOR
!        DO jl = jcs,jce
!          zxt_emis(jl,2) = field%fco2nat(jl,jb)
!        END DO
!        !$ACC END PARALLEL
!        !$ACC END DATA
!      ELSE
!        CALL warning('echam_vdf','co2 emissions not possible in this setup. Please check your settings!')
!      END IF
!    END IF
!!$    ! Dry deposition of aerosols or other tracers (not implemented yet)
!!$    CALL dry_deposition()

    !
    ! setting of co2 for jsbach
    !
    mmr_co2 = vmr_co2   * amco2/amd

    IF(ccycle_config(jg)%iccy_co2conc .EQ. 0) THEN   !  default
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = 1,nproma
        zco2(jl) = mmr_co2
      END DO
      !$ACC END PARALLEL
    ELSE IF (ccycle_config(jg)%iccy_co2conc .EQ. 2 .OR. ccycle_config(jg)%iccy_co2conc .EQ. 4) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = 1,nproma
        zco2(jl) = ghg_co2mmr
      END DO
      !$ACC END PARALLEL
    ELSE IF (ccycle_config(jg)%iccy_co2conc .EQ. 1) THEN
      !$ACC DATA PRESENT( field%qtrc )
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = 1,nproma
        zco2(jl) = field% qtrc(jl,nlev,jb,ico2)
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    ELSE
      CALL finish('echam_vdf: setting of co2 not recommended. Please check your setting!')
    END IF

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          ! Set dummy values to zero to prevent invalid floating point operations:
          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1,nlev
            DO jl = 1,nproma
              dummy (jl,jk) = 0._wp
              dummyx(jl,jk) = 0._wp 
            END DO
          END DO
          !$ACC END PARALLEL
          !
          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1,nlev
            DO jl = jcs,jce
              zqx(jl,jk) =  field%qtrc(jl,jk,jb,iqc) + field%qtrc(jl,jk,jb,iqi)
            END DO
          END DO
          !$ACC END PARALLEL
          !
          ! Serialbox2 intermediate output serialization
          !$ser verbatim call serialize_vdf_chk_A_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
          !$ser verbatim   field, zxt_emis, zco2, dummy, dummyx, zqx)
          !
          IF (ltimer) CALL timer_start(timer_vdf_dn)
          !
          ! Turbulent mixing, part I:
          ! - computation of exchange coefficients in the atmosphere and at the surface;
          ! - build up the tridiagonal linear algebraic system;
          ! - downward sweep (Gaussian elimination from top till level nlev-1)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdf_vd_input(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, nlevm1, nlevp1, ntrac, nsfc_type, iwtr, iice,&
          !$ser verbatim   ilnd, pdtime, field, zqx, zxt_emis, dummy, dummyx,&
          !$ser verbatim   wstar, qs_sfc_tile, hdtcbl, ri_atm, ri_tile, mixlen, cfm,&
          !$ser verbatim   cfm_tile, cfh, cfh_tile, cfv, cftotte, cfthv, zaa, zaa_btm,&
          !$ser verbatim   zbb, zbb_btm, zfactor_sfc, zcpt_sfc_tile, zcptgz, zthvvar, ztottevn)
          !
          ! DA: this routine is async aware, so it's safe not not wait here
          !
          CALL vdiff_down(jg,                              &! in
               &          jb,                              &! in  used for debugging only
               &          jcs, jce, nproma,                &! in
               &          nlev, nlevm1, nlevp1,            &! in
               &          ntrac, nsfc_type,                &! in
               &          iwtr, iice, ilnd,                &! in, indices of different surface types
               &          pdtime,                          &! in, time step
               &          field%coriol(:,jb),              &! in, Coriolis parameter
               &          field%   zf(:,:,jb),             &! in, geopot. height above sea level, full level
               &          field%   zh(:,:,jb),             &! in, geopot. height above sea level, half level
               &          field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &          field% ts_tile(:,jb,:),          &! in, surface temperature
               &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
               &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
               &          field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
               &          field%   ua(:,:,jb),             &! in, um1
               &          field%   va(:,:,jb),             &! in, vm1
               &          field%   ta(:,:,jb),             &! in, tm1
               &          field% qtrc(:,:,jb,iqv),         &! in, qm1
               &          field% qtrc(:,:,jb,iqc),         &! in, xlm1
               &          field% qtrc(:,:,jb,iqi),         &! in, xim1
               &                  zqx(:,:),                &! in, xlm1 + xim1
               &          field% qtrc(:,:,jb,iqt:),        &! in, xtm1
               &          field% mair(:,:,jb),             &! in, moist     air mass
               &          field% mref(:,:,jb),             &! in, reference air mass
               &          field% presi_old(:,:,jb),        &! in, aphm1
               &          field% presm_old(:,:,jb),        &! in, apm1
               &          field%   tv(:,:,jb),             &! in, virtual temperaturea
               &          field% aclc(:,:,jb),             &! in, cloud fraction
               &          zxt_emis,                        &! in, zxtems
               &          dummy(:,:),                      &! in, variance of theta_v at step t-dt
               &          dummyx(:,:),                     &! in
               &          field% z0m_tile(:,jb,:),         &! in
               &          field%  tottem1(:,:,jb),         &! in, TTE at step t-dt
               &          field%  ustar  (:,  jb),         &! inout
               &                  wstar  (:),              &! out, convective velocity scale
               &          field%  wstar_tile(:,jb,:),      &! inout, convective velocity scale (each sfc type)
               &                 qs_sfc_tile(:,   :),      &! out, sfc specific humidity at saturation
               &                 hdtcbl  (:),              &! out, for output
               &                 ri_atm  (:,:),            &! out, for output
               &                 ri_tile (:,:),            &! out, for nsurf_diag
               &                 mixlen  (:,:),            &! out, for output
               &                 cfm     (:,:),            &! out, for output
               &                 cfm_tile(:,   :),         &! out, for output and "vdiff_up"
               &                 cfh     (:,:),            &! out, for output
               &                 cfh_tile(:,   :),         &! out, for output and "vdiff_up"
               &                 cfv     (:,:),            &! out, for output
               &                 cftotte (:,:),            &! out, for output
               &                 cfthv   (:,:),            &! out, for output
               &          zaa, zaa_btm, zbb, zbb_btm,      &! out, for "vdiff_up"
               &          zfactor_sfc(:),                  &! out, for "vdiff_up"
               &          zcpt_sfc_tile(:,:),              &! out, for "vdiff_up"
               &          zcptgz(:,:),                     &! out, for "vdiff_up"
               &          zthvvar(:,:),                    &! out, for "vdiff_up"
               &          field%   thvsig(:,  jb),         &! out, for "cucall"
               &          ztottevn (:,:),                  &! out, for "vdiff_up"
               &          zch_tile(:,:),                   &! out, for "nsurf_diag"
!!$               &          zchn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcdn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcfnc_tile(:,:),                 &! out, for "nsurf_diag"
               &          zbn_tile(:,:),                   &! out, for "nsurf_diag"
               &          zbhn_tile(:,:),                  &! out, for "nsurf_diag"
               &          zbm_tile(:,:),                   &! out, for "nsurf_diag"
               &          zbh_tile(:,:),                   &! out, for "nsurf_diag"
               &          pcsat = field% csat(:,jb),       &! in, optional, area fraction with wet land surface
               &          pcair = field% cair(:,jb),       &! in, optional, area fraction with wet land surface (air)
               &          paz0lh = field% z0h_lnd(:,jb))    ! in, optional, roughness length for heat over land
          !
          ! DA: before the rest of the kernels here are async, we need to wait
          !$ACC WAIT
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_vdf_vd_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, nlevm1, nlevp1, ntrac, nsfc_type, iwtr, iice,&
          !$ser verbatim   ilnd, pdtime, field, wstar, qs_sfc_tile, hdtcbl, ri_atm,&
          !$ser verbatim   ri_tile, mixlen, cfm, cfm_tile, cfh, cfh_tile, cfv,& 
          !$ser verbatim   cftotte, cfthv, zaa, zaa_btm, zbb, zbb_btm, zfactor_sfc,&
          !$ser verbatim   zcpt_sfc_tile, zcptgz, zthvvar, ztottevn, zch_tile,&
          !$ser verbatim   zbn_tile, zbhn_tile, zbm_tile, zbh_tile)
          !
          ! store in memory for output, if requested
          !
          IF (ASSOCIATED(field% wstar)) THEN
            !$ACC DATA PRESENT( field%wstar, wstar )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl=jcs,jce
              field% wstar(jl,jb) = wstar(jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% qs_sfc_tile)) THEN
            !$ACC DATA PRESENT( field%qs_sfc_tile, qs_sfc_tile )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jsfc=1,nsfc_type
              !$ACC LOOP GANG VECTOR
              DO jl=jcs,jce
                field% qs_sfc_tile(jl,jb,jsfc) = qs_sfc_tile(jl,jsfc)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% hdtcbl)) THEN
            !$ACC DATA PRESENT( field%hdtcbl, hdtcbl )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl=jcs,jce
              field% hdtcbl(jl,jb) = hdtcbl(jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% ri_atm)) THEN
            !$ACC DATA PRESENT( field%ri_atm, ri_atm )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% ri_atm(jl,jk,jb) = ri_atm(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% mixlen)) THEN
            !$ACC DATA PRESENT( field%mixlen, mixlen )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% mixlen(jl,jk,jb) = mixlen(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cfm)) THEN
            !$ACC DATA PRESENT( field%cfm, cfm )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% cfm(jl,jk,jb) = cfm(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cfm_tile)) THEN
            !$ACC DATA PRESENT( field%cfm_tile, cfm_tile )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jsfc=1,nsfc_type
              !$ACC LOOP GANG VECTOR
              DO jl=jcs,jce
                field% cfm_tile(jl,jb,jsfc) = cfm_tile(jl,jsfc)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cfh)) THEN
            !$ACC DATA PRESENT( field%cfh, cfh )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% cfh(jl,jk,jb) = cfh(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cfh_tile)) THEN
            !$ACC DATA PRESENT( field%cfh_tile, cfh_tile )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jsfc=1,nsfc_type
              !$ACC LOOP GANG VECTOR
              DO jl=jcs,jce
                field% cfh_tile(jl,jb,jsfc) = cfh_tile(jl,jsfc)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cfv)) THEN
            !$ACC DATA PRESENT( field%cfv, cfv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% cfv(jl,jk,jb) = cfv(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cftotte)) THEN
            !$ACC DATA PRESENT( field%cftotte, cftotte )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% cftotte(jl,jk,jb) = cftotte(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% cfthv)) THEN
            !$ACC DATA PRESENT( field%cfthv, cfthv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk=1,nlev
              !$ACC LOOP VECTOR
              DO jl=jcs,jce
                field% cfthv(jl,jk,jb) = cfthv(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ltimer) CALL timer_stop(timer_vdf_dn)
          !
          !
          ! Surface processes that provide time-dependent lower boundary
          ! condition for wind, temperature, tracer concentration, etc.
          !
          !$ACC DATA PRESENT( field%lhflx_tile, field%shflx_tile, field%evap_tile )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jsfc=1,nsfc_type
            !$ACC LOOP GANG VECTOR
            DO jl=jcs,jce
              field% lhflx_tile(jl,jb,jsfc) = 0._wp
              field% shflx_tile(jl,jb,jsfc) = 0._wp
              field% evap_tile (jl,jb,jsfc) = 0._wp
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA

          !
          ! Serialbox2 intermediate output serialization
          !$ser verbatim call serialize_vdf_chk_B_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
          !$ser verbatim   field)
          !
          IF (ltimer) CALL timer_start(timer_vdf_sf)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdf_us_input(jb, jg, jcs, jce, nproma,&
          !$ser verbatim   nlev, nlevp1, nsfc_type, iwtr, iice, ilnd, pdtime,&
          !$ser verbatim   field, cfh_tile, cfm_tile, zfactor_sfc, zaa, zaa_btm, zbb,&
          !$ser verbatim   zbb_btm, zcpt_sfc_tile, qs_sfc_tile, jb, zco2, zch_tile)
          !
          CALL update_surface(jg, jcs, jce, nproma, field%kice,               &! in
               &              nlev, nsfc_type,                                &! in
               &              iwtr, iice, ilnd,                               &! in, indices of surface types
               &              pdtime,                                         &! in, time step
               &              field%frac_tile(:,jb,:),                        &! in, area fraction
               &                     cfh_tile(:,   :),                        &! in, from "vdiff_down"
               &                     cfm_tile(:,   :),                        &! in, from "vdiff_down"
               &              zfactor_sfc(:),                                 &! in, from "vdiff_down"
               &              field% ocu (:,jb),                              &! in, ocean sfc velocity, u-component
               &              field% ocv (:,jb),                              &! in, ocean sfc velocity, v-component
               &              zaa, zaa_btm, zbb, zbb_btm,                     &! inout
               &               zcpt_sfc_tile(:,:),                            &! inout, from "vdiff_down", for "vdiff_up"
               &                 qs_sfc_tile(:,   :),                         &! inout, from "vdiff_down", for "vdiff_up"
               &              field% ts_tile(:,jb,:),                         &! inout
               &              field%u_stress    (:,  jb),                     &! out
               &              field%v_stress    (:,  jb),                     &! out
               &              field% lhflx      (:,  jb),                     &! out
               &              field% shflx      (:,  jb),                     &! out
               &              field%  evap      (:,  jb),                     &! out, for "cucall"
               &              field%u_stress_tile  (:,jb,:),                  &! out
               &              field%v_stress_tile  (:,jb,:),                  &! out
               &              field% lhflx_tile    (:,jb,:),                  &! out
               &              field% shflx_tile    (:,jb,:),                  &! out
               &              field%  evap_tile    (:,jb,:),                  &! out
               &              field%  fco2nat      (:,  jb),                  &! inout
               &              nblock = jb,                                    &! in
               &              lsm = field%lsmask(:,jb),                       &!< in, land-sea mask
               &              alake = field%alake(:,jb),                      &! in, lake fraction
               &              pu    = field% ua(:,nlev,jb),                   &! in, um1
               &              pv    = field% va(:,nlev,jb),                   &! in, vm1
               &              ptemp = field% ta(:,nlev,jb),                   &! in, tm1
               &              pq = field% qtrc(:,nlev,jb,iqv),                &! in, qm1
               &              pco2 = zco2,                                    &! in, co2, lowest level or fixed value
               &              prsfl = field% rsfl(:,jb),                      &! in, rain surface large scale (from cloud)
               &              prsfc = field% rsfc(:,jb),                      &! in, rain surface concective (from cucall)
               &              pssfl = field% ssfl(:,jb),                      &! in, snow surface large scale (from cloud)
               &              pssfc = field% ssfc(:,jb),                      &! in, snow surface concective (from cucall)
               &              rlds        = field% rlds (:,jb),               &! in,  downward surface  longwave flux [W/m2]
               &              rlus        = field% rlus (:,jb),               &! inout, upward surface  longwave flux [W/m2]
               &              rsds        = field% rsds (:,jb),               &! in,  downward surface shortwave flux [W/m2]
               &              rsus        = field% rsus (:,jb),               &! in,  upward surface shortwave flux [W/m2]
               !
               &              rvds_dir   = field%rvds_dir   (:,jb),           &! in, all-sky downward direct visible radiation at surface
               &              rpds_dir   = field%rpds_dir   (:,jb),           &! in, all-sky downward direct PAR     radiation at surface
               &              rnds_dir   = field%rnds_dir   (:,jb),           &! in, all-sky downward direct near-IR radiation at surface
               &              rvds_dif   = field%rvds_dif   (:,jb),           &! in, all-sky downward diffuse visible radiation at surface
               &              rpds_dif   = field%rpds_dif   (:,jb),           &! in, all-sky downward diffuse PAR     radiation at surface
               &              rnds_dif   = field%rnds_dif   (:,jb),           &! in, all-sky downward diffuse near-IR radiation at surface
               !
               &              ps = field% presi_old(:,nlevp1,jb),             &! in, paphm1, half level pressure
               &              pcosmu0 = field% cosmu0(:,jb),                  &! in, amu0_x, cos of zenith angle
               &              pch_tile = zch_tile(:,:),                       &! in, from "vdiff_down" for JSBACH
               &              pcsat = field%csat(:,jb),                       &! inout, area fraction with wet land surface
               &              pcair = field%cair(:,jb),                       &! inout, area fraction with wet land surface (air)
               &              q_snocpymlt    = q_snocpymlt(:),                &! out, heating  by melting snow on the canopy [W/m2]
               &              z0m_tile = field% z0m_tile(:,jb,:),             &! inout, roughness length for momentum over tiles
               &              z0h_lnd  = field% z0h_lnd (:,jb),               &! out, roughness length for heat over land
               &              albvisdir      = field% albvisdir     (:,jb)  , &! inout
               &              albnirdir      = field% albnirdir     (:,jb)  , &! inout
               &              albvisdif      = field% albvisdif     (:,jb)  , &! inout
               &              albnirdif      = field% albnirdif     (:,jb)  , &! inout
               &              albvisdir_tile = field% albvisdir_tile(:,jb,:), &! inout
               &              albnirdir_tile = field% albnirdir_tile(:,jb,:), &! inout
               &              albvisdif_tile = field% albvisdif_tile(:,jb,:), &! inout
               &              albnirdif_tile = field% albnirdif_tile(:,jb,:), &! inout
               &              albedo         = field% albedo        (:,jb)  , &! inout
               &              albedo_tile    = field% albedo_tile(:,jb,:),    &! inout
               &              emissivity     = field% emissivity    (:,jb)  , &! inout
               &              pco2_flux_tile = field% co2_flux_tile(:,jb,:),  &! inout
               &              ptsfc     = field%ts    (:,jb),                 &! out
               &              ptsfc_rad = field%ts_rad(:,jb),                 &! out
               &              rlns_tile = field%lwflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              rsns_tile = field%swflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              lake_ice_frc = field%lake_ice_frc(:,jb),        &! out
               &              Tsurf = field% Tsurf(:,:,jb),                   &! inout, for sea ice
               &              T1    = field% T1   (:,:,jb),                   &! inout, for sea ice
               &              T2    = field% T2   (:,:,jb),                   &! inout, for sea ice
               &              hi    = field% hi   (:,:,jb),                   &! in, for sea ice
               &              hs    = field% hs   (:,:,jb),                   &! in, for sea ice
               &              conc  = field% conc (:,:,jb),                   &! in, for sea ice
               &              Qtop  = field% Qtop (:,:,jb),                   &! out, for sea ice
               &              Qbot  = field% Qbot (:,:,jb),                   &! out, for sea ice
               &              albvisdir_ice = field% albvisdir_ice(:,:,jb),   &! inout ice albedos
               &              albnirdir_ice = field% albnirdir_ice(:,:,jb),   &! inout
               &              albvisdif_ice = field% albvisdif_ice(:,:,jb),   &! inout
               &              albnirdif_ice = field% albnirdif_ice(:,:,jb))    ! inout
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_vdf_us_output(jb, jg, jcs, jce, nproma,&
          !$ser verbatim   nlev, nsfc_type, iwtr, iice, ilnd,&
          !$ser verbatim   pdtime, field, zaa, zaa_btm, zbb, zbb_btm,&
          !$ser verbatim   zcpt_sfc_tile, qs_sfc_tile, q_snocpymlt)
          !
          IF (ltimer) CALL timer_stop(timer_vdf_sf)
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% q_snocpymlt)) THEN
            !$ACC DATA PRESENT( field%q_snocpymlt, q_snocpymlt )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              field% q_snocpymlt(jl, jb)   = q_snocpymlt(jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF

          !
          !
          !
          ! Turbulent mixing, part II:
          ! - Elimination for the lowest model level using boundary conditions
          !   provided by the surface model(s);
          ! - Back substitution to get solution of the tridiagonal system;
          ! - Compute tendencies and additional diagnostics.
          !
          ! Serialbox2 intermediate output serialization
          !$ser verbatim call serialize_vdf_chk_C_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
          !$ser verbatim   field )
          !
          IF (ltimer) CALL timer_start(timer_vdf_up)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdf_vu_input(jb, jcs, jce, nproma, nlev,&
          !$ser verbatim   nlevm1, ntrac, nsfc_type, iwtr, pdtime, field,&
          !$ser verbatim   cfm_tile, zaa, zcptgz, ztottevn, zbb, zthvvar, dummyx, kedisp)
          !
          IF(ntracer .GT. iqt) THEN
            tend_qtrc_vdf_iqt => tend_qtrc_vdf(:,:,iqt:)
          ELSE
            tend_qtrc_vdf_iqt => tend_qtrc_vdf_dummy
          ENDIF
 
          CALL vdiff_up(jcs, jce, nproma, nlev, nlevm1,  &! in
               &        ntrac, nsfc_type,                &! in
               &        iwtr,                            &! in, indices of different sfc types
               &        pdtime,                          &! in, time steps
               &        field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &               cfm_tile(:,   :),         &! in
               &        zaa,                             &! in, from "vdiff_down"
               &         zcptgz(:,:),                    &! in, from "vdiff_down"
               &        field%   ua(:,:,jb),             &! in, um1
               &        field%   va(:,:,jb),             &! in, vm1
               &        field%   ta(:,:,jb),             &! in, tm1
               &        field% mair(:,:,jb),             &! in, moist     air mass [kg/m2]
               &        field% mref(:,:,jb),             &! in, reference air mass [kg/m2]
               &        field% qtrc(:,:,jb,iqv),         &! in, qm1
               &        field% qtrc(:,:,jb,iqc),         &! in, xlm1
               &        field% qtrc(:,:,jb,iqi),         &! in, xim1
               &        field% qtrc(:,:,jb,iqt:),        &! in, xtm1
               &        field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
               &             ztottevn(:,:),              &! in, TTE at intermediate time step
               &        zbb,                             &! in
               &        zthvvar(:,:),                    &! inout
               &        dummyx(:,:),                     &! inout
               &        field% z0m_tile(:,jb,:),         &! inout
               &                 kedisp(:),              &! out, vert. integr. diss. kin. energy [W/m2]
               &          tend_ua_vdf(:,:),              &! out
               &          tend_va_vdf(:,:),              &! out
               &                q_vdf(:,:),              &! out   heating W/m2
               &        tend_qtrc_vdf(:,:,iqv),          &! out
               &        tend_qtrc_vdf(:,:,iqc),          &! out
               &        tend_qtrc_vdf(:,:,iqi),          &! out
!               &        tend_qtrc_vdf(:,:,iqt:),         &! out
               &        tend_qtrc_vdf_iqt, & ! out
               &        field%   z0m   (:,  jb),         &! out, for the next step
               &        dummy(:,:),                      &! 
               &        field%      totte(:,:,jb),       &! out
               &        field%   sh_vdiff(:,  jb),       &! out, for energy diagnostic
               &        field%   qv_vdiff(:,  jb)        )! out, for energy diagnostic
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_vdf_vu_output(jb, jcs, jce, nproma, nlev,&
          !$ser verbatim   nlevm1, ntrac, nsfc_type, iwtr, pdtime, field, zbb,&
          !$ser verbatim   dummyx, kedisp, tend_ua_vdf, tend_va_vdf, q_vdf,&
          !$ser verbatim   tend_qtrc_vdf, dummy)
          !
          IF (ltimer) CALL timer_stop(timer_vdf_up)
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% kedisp  )) THEN
            !$ACC DATA PRESENT( field%kedisp, kedisp )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              field% kedisp  (jl, jb)   = kedisp (jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(field% q_vdf)) THEN
            !$ACC DATA PRESENT( field%q_vdf, q_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                field% q_vdf(jl,jk,jb) = q_vdf(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_vdf_vi)) THEN
            !$ACC DATA PRESENT( field%q_vdf_vi, q_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              field% q_vdf_vi(jl,jb) = SUM(q_vdf(jl,:))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_vdf)) THEN
            !$ACC DATA PRESENT( tend%ua_vdf, tend_ua_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                tend% ua_vdf(jl,jk,jb) = tend_ua_vdf(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_vdf)) THEN
            !$ACC DATA PRESENT( tend%va_vdf, tend_va_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                tend% va_vdf(jl,jk,jb) = tend_va_vdf(jl,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_vdf )) THEN
            !$ACC DATA PRESENT( tend%qtrc_vdf, tend_qtrc_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                tend% qtrc_vdf(jl,jk,jb,iqv) = tend_qtrc_vdf(jl,jk,iqv)
                tend% qtrc_vdf(jl,jk,jb,iqc) = tend_qtrc_vdf(jl,jk,iqc)
                tend% qtrc_vdf(jl,jk,jb,iqi) = tend_qtrc_vdf(jl,jk,iqi)
                !$ACC LOOP SEQ
                DO jt = iqt,ntracer
                  tend% qtrc_vdf(jl,jk,jb,jt) = tend_qtrc_vdf(jl,jk,jt)
                END DO
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          ! Serialbox2 intermediate output serialization
          !$ser verbatim call serialize_vdf_chk_D_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
          !$ser verbatim   field, tend)
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(field% q_snocpymlt)) THEN
            !$ACC DATA PRESENT( q_snocpymlt, field%q_snocpymlt )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              q_snocpymlt(jl) = field% q_snocpymlt(jl,  jb)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(field% q_vdf)) THEN
            !$ACC DATA PRESENT( q_vdf, field%q_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                q_vdf(jl,jk) = field% q_vdf(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_vdf)) THEN
            !$ACC DATA PRESENT( tend_ua_vdf, tend%ua_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                tend_ua_vdf(jl,jk) = tend% ua_vdf(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_vdf)) THEN
            !$ACC DATA PRESENT( tend_va_vdf, tend%va_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                tend_va_vdf(jl,jk) = tend% va_vdf(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_vdf )) THEN
            !$ACC DATA PRESENT( tend_qtrc_vdf, tend%qtrc_vdf )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1,nlev
              !$ACC LOOP VECTOR
              DO jl = jcs,jce
                tend_qtrc_vdf(jl,jk,iqv) = tend% qtrc_vdf(jl,jk,jb,iqv)
                tend_qtrc_vdf(jl,jk,iqc) = tend% qtrc_vdf(jl,jk,jb,iqc)
                tend_qtrc_vdf(jl,jk,iqi) = tend% qtrc_vdf(jl,jk,jb,iqi)
                !$ACC LOOP SEQ
                DO jt = iqt,ntracer
                  tend_qtrc_vdf(jl,jk,jt) = tend% qtrc_vdf(jl,jk,jb,jt)
                END DO
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          ! Serialbox2 intermediate output serialization
          !$ser verbatim call serialize_vdf_chk_E_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
          !$ser verbatim   field, tend, q_snocpymlt, q_vdf, tend_ua_vdf,&
          !$ser verbatim   tend_va_vdf, tend_qtrc_vdf)
          !
       END IF
       !
       !
       ! Surface effect on the lowermost layer
       !
       IF (echam_phy_config(jg)%ljsb) THEN
          !
          ! convert    heating
          ! q_snocpymlt = heating for melting of snow on canopy
          !             = cooling of atmosphere --> negative sign
          !$ACC DATA PRESENT( tend_ta_sfc, q_snocpymlt, field%qconv )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            tend_ta_sfc(jl) = -q_snocpymlt(jl) * field% qconv(jl,nlev,jb)
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
          !
          IF (ASSOCIATED(tend% ta_sfc)) THEN
            !$ACC DATA PRESENT( tend%ta_sfc, tend_ta_sfc )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs, jce
              tend% ta_sfc(jl,jb) = tend_ta_sfc(jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF

          ! for output: accumulate heating
          IF (ASSOCIATED(field% q_phy)) THEN
            !$ACC DATA PRESENT( field%q_phy, q_snocpymlt )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs, jce
              field% q_phy(jl,nlev,jb) = field% q_phy(jl,nlev,jb) - q_snocpymlt(jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_phy_vi)) THEN
            !$ACC DATA PRESENT( field%q_phy_vi, q_snocpymlt )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jl = jcs, jce
              field% q_phy_vi(jl, jb) = field% q_phy_vi(jl, jb) - q_snocpymlt(jl)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          ! accumulate tendencies for later updating the model state
          SELECT CASE(fc_vdf)
          CASE(0)
             ! diagnostic, do not use tendency
          CASE(1)
             ! use tendency to update the model state
             !$ACC DATA PRESENT( tend%ta_phy, tend_ta_sfc )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG VECTOR
             DO jl = jcs, jce
               tend% ta_phy(jl,nlev,jb) = tend% ta_phy(jl,nlev,jb) + tend_ta_sfc(jl)
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
!!$          CASE(2)
!!$             ! use tendency as forcing in the dynamics
!!$             ...
          END SELECT
          !
          ! update physics state for input to the next physics process
          IF (lparamcpl) THEN
             !$ACC DATA PRESENT( field%ta, tend_ta_sfc )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG VECTOR
             DO jl = jcs, jce
                field% ta(jl,nlev,jb) = field% ta(jl,nlev,jb) + tend_ta_sfc(jl)*pdtime
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
          END IF
          !
          ! Serialbox2 intermediate output serialization
          !$ser verbatim call serialize_vdf_chk_F_output(jg, jb, jcs, jce, nproma,&
          !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
          !$ser verbatim   field, tend, tend_ta_sfc)
          !
       END IF
       !
       !
       ! Vertical diffusion effect on the atmospheric column
       !
       ! convert    heating
       !$ACC DATA PRESENT( q_vdf, field%qconv )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1,nlev
         !$ACC LOOP VECTOR
         DO jl = jcs, jce
           tend_ta_vdf(jl,jk) = q_vdf(jl,jk) * field% qconv(jl,jk,jb)
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       IF (ASSOCIATED(tend% ta_vdf)) THEN
         !$ACC DATA PRESENT( tend%ta_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs, jce
             tend% ta_vdf(jl,jk,jb) = tend_ta_vdf(jl,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF

       ! for output: accumulate heating
       IF (ASSOCIATED(field% q_phy)) THEN
         !$ACC DATA PRESENT( field%q_phy, q_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs, jce
             field% q_phy(jl,jk,jb) = field% q_phy(jl,jk,jb) + q_vdf(jl,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_phy_vi)) THEN 
         !$ACC DATA PRESENT( field%q_phy_vi, q_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs, jce
           field% q_phy_vi(jl, jb) = field% q_phy_vi(jl, jb) + SUM(q_vdf(jl,:))
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_vdf)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          !$ACC DATA PRESENT( tend%ua_phy, tend_ua_vdf, tend%va_phy, tend_va_vdf, tend%ta_phy, &
          !$ACC               tend%qtrc_phy, tend_qtrc_vdf )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1,nlev
            DO jl = jcs, jce
              tend%   ua_phy(jl,jk,jb)      = tend%   ua_phy(jl,jk,jb)      + tend_ua_vdf  (jl,jk)
              tend%   va_phy(jl,jk,jb)      = tend%   va_phy(jl,jk,jb)      + tend_va_vdf  (jl,jk)
              tend%   ta_phy(jl,jk,jb)      = tend%   ta_phy(jl,jk,jb)      + tend_ta_vdf  (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqv)  = tend% qtrc_phy(jl,jk,jb,iqv)  + tend_qtrc_vdf(jl,jk,iqv)
              tend% qtrc_phy(jl,jk,jb,iqc)  = tend% qtrc_phy(jl,jk,jb,iqc)  + tend_qtrc_vdf(jl,jk,iqc)
              tend% qtrc_phy(jl,jk,jb,iqi)  = tend% qtrc_phy(jl,jk,jb,iqi)  + tend_qtrc_vdf(jl,jk,iqi)
              !$ACC LOOP SEQ
              DO jt = iqt,ntracer
                tend% qtrc_phy(jl,jk,jb,jt) = tend% qtrc_phy(jl,jk,jb,jt) + tend_qtrc_vdf(jl,jk,jt)
              END DO
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          !$ACC DATA PRESENT( field%ua, tend_ua_vdf, field%va, tend_va_vdf, field%ta,  &
          !$ACC               field%qtrc, tend_qtrc_vdf )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1,nlev
            DO jl = jcs, jce
              field%   ua(jl,jk,jb)      = field%   ua(jl,jk,jb)      + tend_ua_vdf  (jl,jk)     *pdtime
              field%   va(jl,jk,jb)      = field%   va(jl,jk,jb)      + tend_va_vdf  (jl,jk)     *pdtime
              field%   ta(jl,jk,jb)      = field%   ta(jl,jk,jb)      + tend_ta_vdf  (jl,jk)     *pdtime
              field% qtrc(jl,jk,jb,iqv)  = field% qtrc(jl,jk,jb,iqv)  + tend_qtrc_vdf(jl,jk,iqv) *pdtime
              field% qtrc(jl,jk,jb,iqc)  = field% qtrc(jl,jk,jb,iqc)  + tend_qtrc_vdf(jl,jk,iqc) *pdtime
              field% qtrc(jl,jk,jb,iqi)  = field% qtrc(jl,jk,jb,iqi)  + tend_qtrc_vdf(jl,jk,iqi) *pdtime
              !$ACC LOOP SEQ
              DO jt = iqt,ntracer
                field% qtrc(jl,jk,jb,jt) = field% qtrc(jl,jk,jb,jt) + tend_qtrc_vdf(jl,jk,jt)*pdtime
              END DO
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       !
       ! Correction related to implicitness, due to the fact that surface model only used
       ! part of longwave radiation to compute new surface temperature
       ! 
       !$ACC DATA PRESENT( q_rlw_impl, field%rld_rt, field%rlu_rt, field%rlds, field%rlus, field%q_rlw_nlev )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         q_rlw_impl(jl) =                                               &
              &  ( (field%rld_rt(jl,nlev,jb)-field%rlu_rt(jl,nlev,jb))  & ! ( rln  from "radiation", at top of layer nlev
              &   -(field%rlds  (jl,jb)     -field%rlus  (jl,jb)     )) & !  -rlns from "radheating" and "update_surface")
              & -field%q_rlw_nlev(jl,jb)                                       ! -old heating in layer nlev from "radheating"
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       IF (ASSOCIATED(field%q_rlw_impl)) THEN
         !$ACC DATA PRESENT( field%q_rlw_impl, q_rlw_impl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
          field%q_rlw_impl(jl,jb) = q_rlw_impl(jl)
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF

       ! convert    heating
       !$ACC DATA PRESENT( tend_ta_rlw_impl, q_rlw_impl, field%qconv )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         tend_ta_rlw_impl(jl) = q_rlw_impl(jl) * field% qconv(jl,nlev,jb)
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       IF (ASSOCIATED(tend%ta_rlw_impl)) THEN
         !$ACC DATA PRESENT( tend%ta_rlw_impl, tend_ta_rlw_impl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           tend%ta_rlw_impl(jl,jb) = tend_ta_rlw_impl(jl)
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF

       ! for output: accumulate heating
       IF (ASSOCIATED(field% q_phy)) THEN
         !$ACC DATA PRESENT( field%q_phy, q_rlw_impl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% q_phy(jl,nlev,jb) = field% q_phy(jl,nlev,jb) + q_rlw_impl(jl)
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_phy_vi)) THEN
         !$ACC DATA PRESENT( field%q_phy_vi, q_rlw_impl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% q_phy_vi(jl,jb) = field% q_phy_vi(jl,jb) + q_rlw_impl(jl)
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_vdf)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          !$ACC DATA PRESENT( tend%ta_phy, tend_ta_rlw_impl )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            tend%ta_phy(jl,nlev,jb) = tend%ta_phy(jl,nlev,jb) + tend_ta_rlw_impl(jl)
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          !$ACC DATA PRESENT( field%ta, tend_ta_rlw_impl )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% ta(jl,nlev,jb) = field% ta(jl,nlev,jb) + tend_ta_rlw_impl(jl)*pdtime
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF

       ! 2-tl-scheme
       !$ACC DATA PRESENT( field%tottem1, field%totte )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1,nlev
         !$ACC LOOP VECTOR
         DO jl = jcs,jce
           field% tottem1(jl,jk,jb) = field% totte (jl,jk,jb)
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       ! Serialbox2 intermediate output serialization
       !$ser verbatim call serialize_vdf_chk_G_output(jg, jb, jcs, jce, nproma,&
       !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
       !$ser verbatim   field, tend, tend_ta_vdf, q_rlw_impl, tend_ta_rlw_impl)


       ! Turbulent mixing, part III:
       ! - Further diagnostics.
       !----------------------------------------------------------------------------------------
       ! Serialbox2 input fields serialization
       !$ser verbatim call serialize_vdf_nd_input(jb, jcs, jce, nproma, nlev, nlevp1, nsfc_type,&
       !$ser verbatim   ilnd, field, zqx, zcptgz, zcpt_sfc_tile, zbn_tile,&
       !$ser verbatim   zbhn_tile, zbh_tile, zbm_tile, ri_tile)
       !
       CALL nsurf_diag(jcs, jce, nproma, nsfc_type,     &! in
            &          ilnd,                            &! in
            &          field%frac_tile(:,jb,:),         &! in
            &          field%  qtrc(:,nlev,jb,iqv),     &! in humidity qm1
            &          field%    ta(:,nlev,jb),         &! in tm1
            &          field% presm_old(:,nlev,jb),     &! in, apm1
            &          field% presi_old(:,nlevp1,jb),   &! in, aphm1
            &                  zqx(:,nlev),             &! in, xlm1 + xim1
            &          field%   ua(:,nlev,jb),          &! in, um1
            &          field%   va(:,nlev,jb),          &! in, vm1
            &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
            &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
            &          field% zf  (:,nlev  ,jb),        &! in, height of lowermost full level (m)
            &          field% zh  (:,nlev+1,jb),        &! in, surface height    (m)
            &          zcptgz(:,nlev),                  &! in dry static energy
            &          zcpt_sfc_tile(:,:),              &! in dry static energy
            &          zbn_tile(:,:),                   &! in for diagnostic
            &          zbhn_tile(:,:),                  &! in for diagnostic
            &          zbh_tile(:,:),                   &! in for diagnostic
            &          zbm_tile(:,:),                   &! in for diagnostic
            &          ri_tile(:,:),                    &! in 
            &          field%sfcWind(:,  jb),           &! out 10m windspeed
            &          field%    tas(:,  jb),           &! out temperature in 2m
            &          field%   dew2(:,  jb),           &! out dew point temperature in 2m
            &          field%    uas(:,  jb),           &! out zonal wind in 10m
            &          field%    vas(:,  jb),           &! out meridional wind in 10m
            &          field%tasmax (:,  jb),           &! out max 2m temperature
            &          field%tasmin (:,  jb),           &! out min 2m temperature
            &          field%sfcWind_tile(:,jb,:),      &! out 10m windspeed on tiles
            &          field%    tas_tile(:,jb,:),      &! out temperature in 2m on tiles
            &          field%   dew2_tile(:,jb,:),      &! out dew point temperature in 2m on tiles
            &          field%    uas_tile(:,jb,:),      &! out zonal wind in 10m on tiles
            &          field%    vas_tile(:,jb,:)       )! out meridional wind in 10m on tiles
       !
       !----------------------------------------------------------------------------------------
       ! Serialbox2 output fields serialization
       !$ser verbatim call serialize_vdf_nd_output(jb, jcs, jce, nproma,&
       !$ser verbatim   nsfc_type, ilnd, field)

       !
       !
    ELSE
       !
       !
       ! vdiff_down
       !$ACC DATA PRESENT( field%ustar )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         field% ustar(jl,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       IF (ASSOCIATED(field% wstar)) THEN
         !$ACC DATA PRESENT( field%wstar )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% wstar(jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !$ACC DATA PRESENT( field%wstar_tile )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP SEQ
       DO jsfc = 1,nsfc_type
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% wstar_tile(jl,jb,jsfc) = 0.0_wp
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       IF (ASSOCIATED(field% qs_sfc_tile)) THEN
         !$ACC DATA PRESENT( field%qs_sfc_tile )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jsfc = 1,nsfc_type
           !$ACC LOOP GANG VECTOR
           DO jl = jcs,jce
             field% qs_sfc_tile (jl,jb,jsfc) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% hdtcbl  )) THEN
         !$ACC DATA PRESENT( field%hdtcbl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% hdtcbl(jl, jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% ri_atm  )) THEN
         !$ACC DATA PRESENT( field%ri_atm )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% ri_atm(jk,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% mixlen)) THEN
         !$ACC DATA PRESENT( field%mixlen )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% mixlen(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cfm)) THEN
         !$ACC DATA PRESENT( field%cfm )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% cfm(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cfm_tile)) THEN
         !$ACC DATA PRESENT( field%cfm_tile )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jsfc = 1,nsfc_type
           !$ACC LOOP GANG VECTOR
           DO jl = jcs,jce
             field% cfm_tile(jl,jb,jsfc) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cfh)) THEN
         !$ACC DATA PRESENT( field%cfh )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% cfh(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cfh_tile)) THEN
         !$ACC DATA PRESENT( field%cfh_tile )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jsfc = 1,nsfc_type
           !$ACC LOOP GANG VECTOR
           DO jl = jcs,jce
             field% cfh_tile(jl,jb,jsfc) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cfv)) THEN
         !$ACC DATA PRESENT( field%cfv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% cfv(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cftotte )) THEN
         !$ACC DATA PRESENT( field%cftotte )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% cftotte(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% cfthv)) THEN
         !$ACC DATA PRESENT( field%cfthv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% cfthv(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !$ACC DATA PRESENT( field%thvsig )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         field% thvsig(jl,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       ! update_surface
       !$ACC DATA PRESENT( field%ts_tile, field%u_stress, field%v_stress, field%lhflx, field%shflx, field%evap, &
       !$ACC               field%u_stress_tile, field%v_stress_tile, field%lhflx_tile, field%shflx_tile,        &
       !$ACC               field%evap_tile, field%lwflxsfc_tile, field%swflxsfc_tile, field%rlus, field%rsus,   &
       !$ACC               field%csat, field%cair, field%z0h_lnd, field%albvisdir, field%albnirdir,             &
       !$ACC               field%albvisdif, field%albnirdif, field%albvisdir_tile, field%albnirdir_tile,        &
       !$ACC               field%albvisdif_tile, field%albnirdif_tile, field%albedo, field%albedo_tile,         &
       !$ACC               field%co2_flux_tile, field%ts, field%ts_rad, field%lwflxsfc_tile, field%swflxsfc_tile )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         field% u_stress  (jl,jb) = 0.0_wp
         field% v_stress  (jl,jb) = 0.0_wp
         field% lhflx     (jl,jb) = 0.0_wp
         field% shflx     (jl,jb) = 0.0_wp
         field%  evap     (jl,jb) = 0.0_wp
         field% rlus      (jl,jb) = 0.0_wp
         field% rsus      (jl,jb) = 0.0_wp
         field% csat      (jl,jb) = 0.0_wp
         field% cair      (jl,jb) = 0.0_wp
         field% z0h_lnd   (jl,jb) = 0.0_wp
         field% albvisdir (jl,jb) = 0.0_wp
         field% albnirdir (jl,jb) = 0.0_wp
         field% albvisdif (jl,jb) = 0.0_wp
         field% albnirdif (jl,jb) = 0.0_wp
         field% albedo    (jl,jb) = 0.0_wp
         field% ts        (jl,jb) = 0.0_wp
         field% ts_rad    (jl,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP SEQ
       DO jsfc = 1,nsfc_type
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% ts_tile        (jl,jb,jsfc) = 0.0_wp
           field% u_stress_tile  (jl,jb,jsfc) = 0.0_wp
           field% v_stress_tile  (jl,jb,jsfc) = 0.0_wp
           field% lhflx_tile     (jl,jb,jsfc) = 0.0_wp
           field% shflx_tile     (jl,jb,jsfc) = 0.0_wp
           field%  evap_tile     (jl,jb,jsfc) = 0.0_wp
           field% lwflxsfc_tile  (jl,jb,jsfc) = 0.0_wp
           field% swflxsfc_tile  (jl,jb,jsfc) = 0.0_wp
           field% albvisdir_tile (jl,jb,jsfc) = 0.0_wp
           field% albnirdir_tile (jl,jb,jsfc) = 0.0_wp
           field% albvisdif_tile (jl,jb,jsfc) = 0.0_wp
           field% albnirdif_tile (jl,jb,jsfc) = 0.0_wp
           field% albedo_tile    (jl,jb,jsfc) = 0.0_wp
           field% co2_flux_tile  (jl,jb,jsfc) = 0.0_wp
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       IF (ASSOCIATED(field% q_snocpymlt)) THEN
         !$ACC DATA PRESENT( field%q_snocpymlt )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% q_snocpymlt (jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_rlw_impl)) THEN
         !$ACC DATA PRESENT( field%q_rlw_impl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% q_rlw_impl(jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       !$ACC DATA PRESENT( field%Tsurf, field%T1, field%T2, field%Qtop, field%Qbot, field%albvisdir_ice, &
       !$ACC               field%albnirdir_ice, field%albvisdif_ice, field%albnirdif_ice )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1,nlev
         !$ACC LOOP VECTOR
         DO jl = jcs,jce
           field% Tsurf          (jl,jk,jb) = 0.0_wp
           field% T1             (jl,jk,jb) = 0.0_wp
           field% T2             (jl,jk,jb) = 0.0_wp
           field% Qtop           (jl,jk,jb) = 0.0_wp
           field% Qbot           (jl,jk,jb) = 0.0_wp
           field% albvisdir_ice  (jl,jk,jb) = 0.0_wp
           field% albnirdir_ice  (jl,jk,jb) = 0.0_wp
           field% albvisdif_ice  (jl,jk,jb) = 0.0_wp
           field% albnirdif_ice  (jl,jk,jb) = 0.0_wp
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       IF (ASSOCIATED(tend% ta_sfc     )) THEN
         !$ACC DATA PRESENT( tend%ta_sfc )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           tend% ta_sfc(jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ta_rlw_impl)) THEN
         !$ACC DATA PRESENT( tend%ta_rlw_impl )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           tend% ta_rlw_impl(jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! vdiff_up
       !$ACC DATA PRESENT( field%totte, field%z0m, field%z0m_tile )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1,nlev
         !$ACC LOOP VECTOR
         DO jl = jcs,jce
           field% totte(jl,jk,jb) = 0.0_wp
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         field% z0m(jl,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP SEQ
       DO jsfc = 1,nsfc_type
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% z0m_tile(jl,jb,jsfc) = 0.0_wp
         END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       IF (ASSOCIATED(field% kedisp)) THEN
         !$ACC DATA PRESENT( field%kedisp )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% kedisp(jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !$ACC DATA PRESENT( field%sh_vdiff, field%qv_vdiff )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,jce
         field% sh_vdiff(jl,jb) = 0.0_wp
         field% qv_vdiff(jl,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       IF (ASSOCIATED(field% q_vdf)) THEN
         !$ACC DATA PRESENT( field%q_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
         !$ACC LOOP VECTOR
           DO jl = jcs,jce
             field% q_vdf(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_vdf_vi)) THEN
         !$ACC DATA PRESENT( field%q_vdf_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jl = jcs,jce
           field% q_vdf_vi(jl,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% ta_vdf)) THEN
         !$ACC DATA PRESENT( tend%ta_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
         !$ACC LOOP VECTOR
           DO jl = jcs,jce
             tend% ta_vdf(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ua_vdf)) THEN
         !$ACC DATA PRESENT( tend%ua_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
         !$ACC LOOP VECTOR
           DO jl = jcs,jce
             tend% ua_vdf(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% va_vdf)) THEN
         !$ACC DATA PRESENT( tend%va_vdf )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
         !$ACC LOOP VECTOR
           DO jl = jcs,jce
             tend% va_vdf(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% qtrc_vdf)) THEN
          !$ACC DATA PRESENT( tend%qtrc_vdf )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1,nlev
          !$ACC LOOP VECTOR
            DO jl = jcs,jce
              tend% qtrc_vdf(jl,jk,jb,iqv)  = 0.0_wp
              tend% qtrc_vdf(jl,jk,jb,iqc)  = 0.0_wp
              tend% qtrc_vdf(jl,jk,jb,iqi)  = 0.0_wp
              !$ACC LOOP SEQ
              DO jt = iqt,ntracer
                tend% qtrc_vdf(jl,jk,jb,jt) = 0.0_wp
              END DO
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF
       !
       ! Serialbox2 intermediate output serialization
       !$ser verbatim call serialize_vdf_chk_H_output(jg, jb, jcs, jce, nproma,&
       !$ser verbatim   nlev, ntrac, nsfc_type, pdtime,&
       !$ser verbatim   field, tend)
       !
    END IF

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_vdf_output(jg, jb, jcs, jce, nproma, nlev, field, tend)

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_vdf)
    NULLIFY(vmr_co2)
    NULLIFY(field)
    NULLIFY(tend)

    !$ACC END DATA
    !$ACC END DATA

    IF (ltimer) CALL timer_stop(timer_vdf)

  END SUBROUTINE interface_echam_vdf

END MODULE mo_interface_echam_vdf
