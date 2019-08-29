!>
!! @brief Subroutine interface_echam_cld calls the Lohmann&Roeckner cloud scheme.
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

MODULE mo_interface_echam_cld

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
       &                            t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cld

  USE mo_run_config          ,ONLY: iqv, iqc, iqi
  USE mo_cloud               ,ONLY: cloud
  !$ser verbatim USE mo_ser_echam_cld, ONLY: serialize_cld_input,&
  !$ser verbatim                             serialize_cld_output

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_cld

CONTAINS

  SUBROUTINE interface_echam_cld(jg,jb,jcs,jce        ,&
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
    INTEGER                 ,POINTER    :: fc_cld
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER                             :: itype(nproma)    !< type of convection
    INTEGER                             :: jc,jk
    !
    REAL(wp)                            :: hur  (nproma,nlev)
    REAL(wp)                            :: q_cld(nproma,nlev)
    !
    REAL(wp)                            :: tend_ta_cld  (nproma,nlev)
    REAL(wp)                            :: tend_qtrc_cld(nproma,nlev,ntracer)

    IF (ltimer) call timer_start(timer_cld)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_cld    => echam_phy_config(jg)%fc_cld
    field     => prm_field(jg)
    tend      => prm_tend(jg)

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_cld_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    IF ( is_in_sd_ed_interval ) THEN
       !$ACC DATA PRESENT( field%rtype, field% qconv ) &
       !$ACC       CREATE( itype, hur, q_cld, tend_ta_cld, tend_qtrc_cld )
       !
       IF ( is_active ) THEN
          !
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jc = 1, nproma
            itype(jc) = NINT(field%rtype(jc,jb))
          END DO
          !$ACC END PARALLEL
          !
          CALL cloud(jg,                           &! in
               &     jb,                           &! in
               &     jcs, jce, nproma, nlev,       &! in
               &     pdtime,                       &! in
               &     field% ictop    (:,  jb),     &! in (from "cucall")
               &     field% presm_old(:,:,jb),     &! in
               &     field% dz       (:,:,jb),     &! in
               &     field% mref     (:,:,jb),     &! in
               &     field% rho      (:,:,jb),     &! in
               &     field% cpair    (:,:,jb),     &! in
               &     field% acdnc    (:,:,jb),     &! in  acdnc
               &     field% ta       (:,:,jb),     &! in  tm1
               &     field% qtrc     (:,:,jb,iqv), &! in  qm1
               &     field% qtrc     (:,:,jb,iqc), &! in  xlm1
               &     field% qtrc     (:,:,jb,iqi), &! in  xim1
               !
               &     itype,                        &! inout
               &     field% aclc     (:,:,jb),     &! inout
               !
               &     field% aclcov   (:,  jb),     &! out
               &     field% rsfl     (:,  jb),     &! out
               &     field% ssfl     (:,  jb),     &! out
               &            hur      (:,:),        &! out
               &            q_cld    (:,:),        &! out
               &      tend_qtrc_cld  (:,:,iqv),        &! out
               &      tend_qtrc_cld  (:,:,iqc),        &! out
               &      tend_qtrc_cld  (:,:,iqi)         )! out
          !
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jc = 1, nproma
            field% rtype(jc,jb) = REAL(itype(jc),wp)
          END DO
          !$ACC END PARALLEL
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% hur)) THEN
            !$ACC DATA PRESENT( field%hur )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field%hur(jc,jk,jb) = hur(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(field% q_cld)) THEN
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field% q_cld(jc,jk,jb) = q_cld(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
          END IF
          IF (ASSOCIATED(field% q_cld_vi)) THEN
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_cld_vi(jc,  jb) = SUM(q_cld(jc,:))
            END DO
            !$ACC END PARALLEL
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_cld )) THEN
            !$ACC DATA PRESENT( tend, tend%qtrc_cld )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
               tend% qtrc_cld(jc,jk,jb,iqv) = tend_qtrc_cld(jc,jk,iqv)
               tend% qtrc_cld(jc,jk,jb,iqc) = tend_qtrc_cld(jc,jk,iqc)
               tend% qtrc_cld(jc,jk,jb,iqi) = tend_qtrc_cld(jc,jk,iqi)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(field% q_cld)) THEN
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                q_cld(jc,jk) = field% q_cld(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_cld )) THEN
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                tend_qtrc_cld(jc,jk,iqv) = tend%qtrc_cld(jc,jk,jb,iqv)
                tend_qtrc_cld(jc,jk,iqc) = tend%qtrc_cld(jc,jk,jb,iqc)
                tend_qtrc_cld(jc,jk,iqi) = tend%qtrc_cld(jc,jk,jb,iqi)
              END DO
            END DO
            !$ACC END PARALLEL
          END IF
          !
       END IF
       !
       ! convert    heating
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1, nlev
         !$ACC LOOP VECTOR
         DO jc = jcs, jce
           tend_ta_cld(jc,jk) = q_cld(jc,jk) * field% qconv(jc,jk,jb)
         END DO
       END DO
       !$ACC END PARALLEL
       !
       IF (ASSOCIATED(tend% ta_cld)) THEN
         !$ACC DATA PRESENT( tend, tend%ta_cld )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_cld(jc,jk,jb) = tend_ta_cld(jc,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       
       ! for output: accumulate heating
       IF (ASSOCIATED(field% q_phy   )) THEN
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_phy(jc,jk,jb) = field% q_phy(jc,jk,jb) + q_cld(jc,jk)
           END DO
         END DO
         !$ACC END PARALLEL
       END IF
       IF (ASSOCIATED(field% q_phy_vi)) THEN
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_phy_vi(jc,jb) = field% q_phy_vi(jc, jb) + SUM(q_cld(jc,:))
         END DO
         !$ACC END PARALLEL
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_cld)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = jcs, jce
              tend%   ta_phy(jc,jk,jb)      = tend%   ta_phy(jc,jk,jb)     + tend_ta_cld  (jc,jk)
              tend% qtrc_phy(jc,jk,jb,iqv)  = tend% qtrc_phy(jc,jk,jb,iqv) + tend_qtrc_cld(jc,jk,iqv)
              tend% qtrc_phy(jc,jk,jb,iqc)  = tend% qtrc_phy(jc,jk,jb,iqc) + tend_qtrc_cld(jc,jk,iqc)
              tend% qtrc_phy(jc,jk,jb,iqi)  = tend% qtrc_phy(jc,jk,jb,iqi) + tend_qtrc_cld(jc,jk,iqi)
            END DO
          END DO
          !$ACC END PARALLEL
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
         SELECT CASE(fc_cld)
         CASE(0)
            ! diagnostic, do not use tendency
         CASE(1)
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field%   ta(jc,jk,jb)      = field%   ta(jc,jk,jb)      + tend_ta_cld(jc,jk)*pdtime
                field% qtrc(jc,jk,jb,iqv)  = field% qtrc(jc,jk,jb,iqv)  + tend_qtrc_cld(jc,jk,iqv)*pdtime
                field% qtrc(jc,jk,jb,iqc)  = field% qtrc(jc,jk,jb,iqc)  + tend_qtrc_cld(jc,jk,iqc)*pdtime
                field% qtrc(jc,jk,jb,iqi)  = field% qtrc(jc,jk,jb,iqi)  + tend_qtrc_cld(jc,jk,iqi)*pdtime
              END DO
            END DO
            !$ACC END PARALLEL
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
         END SELECT
       END IF
       !
       !$ACC END DATA
    ELSE
       !
       IF (ASSOCIATED(field% q_cld   )) THEN
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_cld (jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
       END IF
       IF (ASSOCIATED(field% q_cld_vi)) THEN
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_cld_vi(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
       END IF
       !
       IF (ASSOCIATED(tend% ta_cld)) THEN
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_cld(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
       END IF
       !
       IF (ASSOCIATED(tend% qtrc_cld)) THEN
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% qtrc_cld(jc,jk,jb,iqv) = 0.0_wp
             tend% qtrc_cld(jc,jk,jb,iqc) = 0.0_wp
             tend% qtrc_cld(jc,jk,jb,iqi) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
       END IF
       !
    END IF

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_cld_output(jg, jb, jcs, jce, nproma, nlev, field, tend)

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_cld)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_cld)

  END SUBROUTINE interface_echam_cld

END MODULE mo_interface_echam_cld
