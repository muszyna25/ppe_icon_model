!>
!! @brief Subroutine interface_echam_cnv calls the Tiedtke/Nordeng convection scheme.
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

MODULE mo_interface_echam_cnv

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cnv

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqt
  USE mo_cumastr             ,ONLY: cumastr
  !$ser verbatim USE mo_ser_echam_cnv, ONLY: serialize_cnv_input,&
  !$ser verbatim                             serialize_cnv_output

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_cnv

CONTAINS

  SUBROUTINE interface_echam_cnv(jg,jb,jcs,jce        ,&
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
    LOGICAL                             :: lparamcpl
    INTEGER                             :: fc_cnv
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER                             :: itype(nproma)    !< type of convection
    !
    REAL(wp)                            :: q_cnv(nproma,nlev)
    !
    REAL(wp)                            :: tend_ta_cnv  (nproma,nlev)
    REAL(wp)                            :: tend_ua_cnv  (nproma,nlev)
    REAL(wp)                            :: tend_va_cnv  (nproma,nlev)
    REAL(wp), TARGET                    :: tend_qtrc_cnv(nproma,nlev,ntracer)
    REAL(wp), POINTER, CONTIGUOUS :: tend_qtrc_cnv_iqt(:,:,:)
    REAL(wp), TARGET :: tend_qtrc_cnv_dummy(nproma,nlev,0)
    !
    INTEGER  :: nlevm1, nlevp1
    INTEGER  :: ntrac                        !< # of tracers excluding water vapour and hydrometeors
                                             !< which are convectively transported

    LOGICAL  :: ldland(nproma)               !< land sea mask for using dlev_land or dlev_ocean

    INTEGER  :: ictop(nproma)                !< level index of cnovective cloud top
    REAL(wp) :: ztop(nproma)                 !< convective cloud top pressure   [Pa]
    REAL(wp) :: zta    (nproma,nlev)         !< provisional temperature         [K]
    REAL(wp), TARGET :: zqtrc  (nproma,nlev,ntracer) !< provisional mass mixing ratios  [kg/kg]
    REAL(wp) :: zua    (nproma,nlev)         !< provisional zonal      wind     [m/s]
    REAL(wp) :: zva    (nproma,nlev)         !< provisional meridional wind     [m/s]
    !
    REAL(wp) :: zqtrc_cnd(nproma,nlev)       !< cloud condensate mixing ratio   [kg/kg]
    REAL(wp) :: ztend_qv(nproma,nlev)        !< moisture tendency from dynamics and physics before convection
    REAL(wp), POINTER, CONTIGUOUS :: zqtrc_iqt(:,:,:)
    REAL(wp), TARGET :: zqtrc_dummy(nproma,nlev,0)
    INTEGER  :: jc, jk, jt

    IF (ltimer) CALL timer_start(timer_cnv)

    ! associate pointers
    lparamcpl = echam_phy_config(jg)%lparamcpl
    fc_cnv    = echam_phy_config(jg)%fc_cnv
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_cnv_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    nlevm1 = nlev-1
    nlevp1 = nlev+1

    !$ACC DATA CREATE( itype, q_cnv, tend_ta_cnv, tend_ua_cnv, tend_va_cnv, tend_qtrc_cnv, ldland, &
    !$ACC              ztop, zta, zqtrc, zua, zva, zqtrc_cnd, ztend_qv )
    
    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !$ACC DATA PRESENT( field%ta, field%qtrc, field%ua, field%va, tend%qtrc_dyn, tend%qtrc_phy, field%sftlf, field%rtype )
          !
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = jcs, jce
              zta  (jc,jk) = field% ta(jc,jk,jb)
              !$ACC LOOP SEQ
              DO jt = 1, ntracer
                zqtrc(jc,jk,jt) = MAX(field% qtrc(jc,jk,jb,jt), 0.0_wp)
              END DO
              zua  (jc,jk) =     field% ua  (jc,jk,jb)
              zva  (jc,jk) =     field% va  (jc,jk,jb)
            END DO
          END DO
          !$ACC END PARALLEL
          !
          ! Prepare input
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            DO jc = jcs, jce
              zqtrc_cnd(jc,jk)   = zqtrc(jc,jk,iqc) + zqtrc(jc,jk,iqi)
              ztend_qv (jc,jk)   = tend%qtrc_dyn(jc,jk,jb,iqv) + tend%qtrc_phy(jc,jk,jb,iqv)
            END DO
          END DO
          !$ACC END PARALLEL
          !
          ! number of tracers excluding water vapour and hydrometeors
          ntrac = ntracer-iqt+1
          !
          ! land sea mask for using dlev_land or dlev_ocean
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            ldland(jc)    = field% sftlf(jce,jb) > 0._wp
          END DO
          !$ACC END PARALLEL
          !
          IF(ntracer .GT. iqt) THEN
            zqtrc_iqt => zqtrc(:,:,iqt:)
            tend_qtrc_cnv_iqt => tend_qtrc_cnv  (:,:,iqt:)
          ELSE
            zqtrc_iqt => zqtrc_dummy
            tend_qtrc_cnv_iqt => tend_qtrc_cnv_dummy
          END IF
          CALL cumastr(jg, jb,                       &! in
               &       jcs, jce, nproma,             &! in
               &       nlev, nlevp1, nlevm1,         &! in
               &       pdtime,                       &! in
               &       field% zf       (:,:,jb),     &! in
               &       field% zh       (:,:,jb),     &! in
               &       field% mref     (:,:,jb),     &! in
               &             zta       (:,:),        &! in
               &             zqtrc     (:,:,   iqv), &! in
               &             zqtrc_cnd (:,:),        &! in
               &             zua       (:,:),        &! in
               &             zva       (:,:),        &! in
               &       ntrac,                        &! in
               &             ldland    (:),          &! in
!               &             zqtrc     (:,:,   iqt:),&! in
               &             zqtrc_iqt (:,:,:),&! in
               &       field% omega    (:,:,jb),     &! in
               &       field% evap     (:,  jb),     &! in
               &       field% presm_new(:,:,jb),     &! in
               &       field% presi_new(:,:,jb),     &! in
               &       field% geom     (:,:,jb),     &! in
               &       field% geoi     (:,:,jb),     &! in
               &             ztend_qv  (:,:),        &! in
               &       field% thvsig   (:,  jb),     &! in
               &       itype           (:),          &! out
               &       ictop           (:),          &! out
               &       field% rsfc     (:,  jb),     &! out
               &       field% ssfc     (:,  jb),     &! out
               &       field% con_dtrl (:,jb),       &! out
               &       field% con_dtri (:,jb),       &! out
               &       field% con_iteqv(:,jb),       &! out
               &              q_cnv    (:,:),        &! out
               &        tend_ua_cnv    (:,:),        &! out
               &        tend_va_cnv    (:,:),        &! out
               &        tend_qtrc_cnv  (:,:,iqv),    &! out
!               &        tend_qtrc_cnv  (:,:,iqt:),   &! out
               &        tend_qtrc_cnv_iqt(:,:,:),   &! out
               &        tend_qtrc_cnv  (:,:,iqc),    &! out
               &        tend_qtrc_cnv  (:,:,iqi),    &! out
               &             ztop      (:)           )! out
          !
          ! keep minimum conv. cloud top pressure (= max. conv. cloud top height) of this output interval
          IF (ASSOCIATED(field% topmax)) THEN
            !$ACC DATA PRESENT( field%topmax )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% topmax(jc,jb) = MIN(field% topmax(jc,jb),ztop(jc))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% q_cnv)) THEN
            !$ACC DATA PRESENT( field% q_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field% q_cnv(jc,jk,jb) = q_cnv(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_cnv_vi)) THEN
            !$ACC DATA PRESENT( field% q_cnv_vi )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_cnv_vi(jc,jb) = SUM(q_cnv(jc,:))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_cnv)) THEN
            !$ACC DATA PRESENT( tend% ua_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend% ua_cnv(jc,jk,jb) = tend_ua_cnv(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_cnv)) THEN
            !$ACC DATA PRESENT( tend% va_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend% va_cnv(jc,jk,jb) = tend_va_cnv(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_cnv )) THEN
            !$ACC DATA PRESENT( tend% qtrc_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend% qtrc_cnv(jc,jk,jb,iqv)  = tend_qtrc_cnv(jc,jk,iqv)
                tend% qtrc_cnv(jc,jk,jb,iqc)  = tend_qtrc_cnv(jc,jk,iqc)
                tend% qtrc_cnv(jc,jk,jb,iqi)  = tend_qtrc_cnv(jc,jk,iqi)
                !$ACC LOOP SEQ
                DO jt = iqt, ntracer
                  tend% qtrc_cnv(jc,jk,jb,jt) = tend_qtrc_cnv(jc,jk,jt)
                END DO
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF

          !$ACC END DATA
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(field% q_cnv)) THEN
            !$ACC DATA PRESENT( field% q_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                q_cnv(jc,jk) = field% q_cnv(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_cnv)) THEN
            !$ACC DATA PRESENT( tend% ua_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend_ua_cnv(jc,jk) = tend% ua_cnv(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_cnv)) THEN
            !$ACC DATA PRESENT( tend% va_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend_va_cnv(jc,jk) = tend% va_cnv(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_cnv )) THEN
            !$ACC DATA PRESENT( tend% qtrc_cnv )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend_qtrc_cnv(jc,jk,iqv)  = tend% qtrc_cnv(jc,jk,jb,iqv)
                tend_qtrc_cnv(jc,jk,iqc)  = tend% qtrc_cnv(jc,jk,jb,iqc)
                tend_qtrc_cnv(jc,jk,iqi)  = tend% qtrc_cnv(jc,jk,jb,iqi)
                !$ACC LOOP SEQ
                DO jt = iqt, ntracer
                  tend_qtrc_cnv(jc,jk,jt) = tend% qtrc_cnv(jc,jk,jb,jt)
                END DO
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
       END IF
       !
       ! convert    heating
       !$ACC DATA PRESENT( field%qconv )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1, nlev
          !$ACC LOOP VECTOR
          DO jc = jcs, jce
            tend_ta_cnv(jc,jk) = q_cnv(jc,jk) * field% qconv(jc,jk,jb)
          END DO
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       IF (ASSOCIATED(tend% ta_cnv)) THEN
         !$ACC DATA PRESENT( tend%ta_cnv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = jcs, jce
              tend% ta_cnv(jc,jk,jb) = tend_ta_cnv(jc,jk)
            END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF

       ! for output: accumulate heating
       IF (ASSOCIATED(field% q_phy)) THEN
         !$ACC DATA PRESENT( field%q_phy )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_phy(jc,jk,jb) = field% q_phy(jc,jk,jb) +     q_cnv(jc,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_phy_vi)) THEN
         !$ACC DATA PRESENT( field%q_phy_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_phy_vi(jc, jb) = field% q_phy_vi(jc, jb) + SUM(q_cnv(jc,:))
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_cnv)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          !$ACC DATA PRESENT( tend%ua_phy, tend%va_phy, tend%ta_phy, tend%qtrc_phy )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = jcs, jce
              ! use tendency to update the model state
              tend%   ua_phy(jc,jk,jb)      = tend%   ua_phy(jc,jk,jb)      + tend_ua_cnv  (jc,jk)
              tend%   va_phy(jc,jk,jb)      = tend%   va_phy(jc,jk,jb)      + tend_va_cnv  (jc,jk)
              tend%   ta_phy(jc,jk,jb)      = tend%   ta_phy(jc,jk,jb)      + tend_ta_cnv  (jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqv)  = tend% qtrc_phy(jc,jk,jb,iqv)  + tend_qtrc_cnv(jc,jk,iqv)
              tend% qtrc_phy(jc,jk,jb,iqc)  = tend% qtrc_phy(jc,jk,jb,iqc)  + tend_qtrc_cnv(jc,jk,iqc)
              tend% qtrc_phy(jc,jk,jb,iqi)  = tend% qtrc_phy(jc,jk,jb,iqi)  + tend_qtrc_cnv(jc,jk,iqi)
              !$ACC LOOP SEQ
              DO jt = iqt, ntracer
                tend% qtrc_phy(jc,jk,jb,jt) = tend% qtrc_phy(jc,jk,jb,jt) + tend_qtrc_cnv(jc,jk,jt)
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
       SELECT CASE(fc_cnv)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1,2)
          ! use tendency to update the physics state
          IF (lparamcpl) THEN
             ! prognostic
             !$ACC DATA PRESENT( field%ua, field%va, field%ta, field%qtrc )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = 1, nlev
               !$ACC LOOP VECTOR
               DO jc = jcs, jce
                 field%   ua(jc,jk,jb)      = field%   ua(jc,jk,jb)      + tend_ua_cnv  (jc,jk)     *pdtime
                 field%   va(jc,jk,jb)      = field%   va(jc,jk,jb)      + tend_va_cnv  (jc,jk)     *pdtime
                 field%   ta(jc,jk,jb)      = field%   ta(jc,jk,jb)      + tend_ta_cnv  (jc,jk)     *pdtime
                 field% qtrc(jc,jk,jb,iqv)  = field% qtrc(jc,jk,jb,iqv)  + tend_qtrc_cnv(jc,jk,iqv) *pdtime
                 field% qtrc(jc,jk,jb,iqc)  = field% qtrc(jc,jk,jb,iqc)  + tend_qtrc_cnv(jc,jk,iqc) *pdtime
                 field% qtrc(jc,jk,jb,iqi)  = field% qtrc(jc,jk,jb,iqi)  + tend_qtrc_cnv(jc,jk,iqi) *pdtime
                 !$ACC LOOP SEQ
                 DO jt = iqt, ntracer
                   field% qtrc(jc,jk,jb,jt) = field% qtrc(jc,jk,jb,jt) + tend_qtrc_cnv(jc,jk,jt)*pdtime
                 END DO
               END DO
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
             !
             ! diagnostic
             ! store convection type as real value
             !$ACC DATA PRESENT( field%rtype, field%ictop )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG VECTOR
             DO jc = jcs, jce
               field% rtype(jc,jb) = REAL(itype(jc),wp)
               field% ictop(jc,jb) = ictop(jc)
             END DO
             !
             !$ACC END PARALLEL
             !$ACC END DATA
          END IF
       END SELECT
       !
    ELSE
       !$ACC DATA PRESENT( field%rtype, field%ictop, field%rsfc, field%ssfc, field%con_dtrl, field%con_dtri, &
       !$ACC               field%con_iteqv )
       !
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jc = jcs, jce
         field% rtype    (jc,jb) = 0.0_wp
         field% ictop    (jc,jb) = nlevm1
         field% rsfc     (jc,jb) = 0.0_wp
         field% ssfc     (jc,jb) = 0.0_wp
         field% con_dtrl (jc,jb) = 0.0_wp
         field% con_dtri (jc,jb) = 0.0_wp
         field% con_iteqv(jc,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !
       IF (ASSOCIATED(field% q_cnv   )) THEN
         !$ACC DATA PRESENT( field%q_cnv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1,nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_cnv(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_cnv_vi)) THEN
         !$ACC DATA PRESENT( field%q_cnv_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_cnv_vi(jc,  jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% ta_cnv)) THEN
         !$ACC DATA PRESENT( tend%ta_cnv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_cnv(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ua_cnv)) THEN
         !$ACC DATA PRESENT( tend%ua_cnv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ua_cnv(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% va_cnv)) THEN
         !$ACC DATA PRESENT( tend%va_cnv )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% va_cnv(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% qtrc_cnv)) THEN
          !$ACC DATA PRESENT( tend%qtrc_cnv )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = jcs, jce
              tend% qtrc_cnv(jc,jk,jb,iqv)  = 0.0_wp
              tend% qtrc_cnv(jc,jk,jb,iqc)  = 0.0_wp
              tend% qtrc_cnv(jc,jk,jb,iqi)  = 0.0_wp
              DO jt = iqt, ntracer
                tend% qtrc_cnv(jc,jk,jb,jt) = 0.0_wp
              END DO
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF
       !
       !$ACC END DATA
       !
    END IF

    !$ACC END DATA

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_cnv_output(jg, jb, jcs, jce, nproma, nlev, field, tend)

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_cnv)

  END SUBROUTINE interface_echam_cnv

END MODULE mo_interface_echam_cnv
