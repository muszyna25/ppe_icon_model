!>
!! @brief Subroutine interface_echam_gwd calls the Hines gravity wave scheme.
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

MODULE mo_interface_echam_gwd

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
       &                            t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_gwd

  USE mo_gw_hines            ,ONLY: gw_hines
!!$  USE mo_math_constants      ,ONLY: pi
  !$ser verbatim USE mo_ser_echam_gwd, ONLY: serialize_gwd_input,&
  !$ser verbatim                             serialize_gwd_output

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_gwd

CONTAINS

  SUBROUTINE interface_echam_gwd(jg, jb,jcs,jce       ,&
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
    INTEGER                             :: fc_gwd
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp)                            :: q_gwd(nproma,nlev)
    !
    REAL(wp)                            :: tend_ta_gwd(nproma,nlev)
    REAL(wp)                            :: tend_ua_gwd(nproma,nlev)
    REAL(wp)                            :: tend_va_gwd(nproma,nlev)
    !
    REAL(wp)                            :: zdis_gwd(nproma,nlev) !< out, energy dissipation rate [J/s/kg]
    INTEGER                             :: nc, jc, jk            !< number of cells/columns from jcs to jce
!!$    REAL(wp)                            :: zlat_deg(nproma)       !< latitude in deg N

    IF (ltimer) call timer_start(timer_gwd)    

    ! associate pointers
    lparamcpl = echam_phy_config(jg)%lparamcpl
    fc_gwd    = echam_phy_config(jg)%fc_gwd
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_gwd_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    IF ( is_in_sd_ed_interval ) THEN

       !$ACC DATA CREATE( q_gwd, tend_ta_gwd, tend_ua_gwd, tend_va_gwd, zdis_gwd, field%mair, field%qconv )
       !
       IF ( is_active ) THEN
          !
          ! number of cells/columns from index jcs to jce
          nc = jce-jcs+1
          !
!!$          ! latitude [degN]
!!$          zlat_deg(jcs:jce) = field% clat(jcs:jce,jb) * 180._wp/pi
          !
          CALL gw_hines (jg                       ,&
               &         nproma                   ,&
               &         jcs                      ,&
               &         jce                      ,&
               &         nc                       ,&
               &         nlev                     ,&
               &         field% presi_old(:,:,jb) ,&
               &         field% presm_old(:,:,jb) ,&
               &         field%   zh(:,:,jb)      ,&
               &         field%  rho(:,:,jb)      ,&
               &         field% mair(:,:,jb)      ,&
               &         field%   ta(:,:,jb)      ,&
               &         field%   ua(:,:,jb)      ,&
               &         field%   va(:,:,jb)      ,&
!!$               &         zlat_deg(:)              ,&
!!$               &         aprflux(:,krow)          ,&
               &         zdis_gwd(:,:)            ,&
               &         tend_ua_gwd(:,:)         ,&
               &         tend_va_gwd(:,:)          )
          !
          ! heating
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              q_gwd(jc,jk) = zdis_gwd(jc,jk) * field%mair(jc,jk,jb)
            END DO
          END DO
          !$ACC END PARALLEL
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% q_gwd)) THEN
            !$ACC DATA PRESENT( field%q_gwd )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                field% q_gwd(jc,jk,jb) = q_gwd(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_gwd_vi)) THEN
            !$ACC DATA PRESENT( field%q_gwd_vi )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_gwd_vi(jc,jb) = SUM(q_gwd(jc,:))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_gwd)) THEN
            !$ACC DATA PRESENT( tend%ua_gwd )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend% ua_gwd(jc,jk,jb) = tend_ua_gwd(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_gwd)) THEN
            !$ACC DATA PRESENT( tend%va_gwd )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend% va_gwd(jc,jk,jb) = tend_va_gwd(jc,jk)
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
          IF (ASSOCIATED(field% q_gwd)) THEN
            !$ACC DATA PRESENT( field%q_gwd )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                q_gwd(jc,jk) = field% q_gwd(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_gwd)) THEN
            !$ACC DATA PRESENT( tend%ua_gwd )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend_ua_gwd(jc,jk) = tend% ua_gwd(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_gwd)) THEN
            !$ACC DATA PRESENT( tend%va_gwd )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend_va_gwd(jc,jk) = tend% va_gwd(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
       END IF
       !
       ! convert    heating
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP SEQ
       DO jk = 1, nlev
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           tend_ta_gwd(jc,jk) = q_gwd(jc,jk) * field% qconv(jc,jk,jb)
         END DO
       END DO
       !$ACC END PARALLEL
       !
       IF (ASSOCIATED(tend% ta_gwd)) THEN
         !$ACC DATA PRESENT( tend%ta_gwd )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% ta_gwd(jc,jk,jb) = tend_ta_gwd(jc,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF

       ! for output: accumulate heating
       IF (ASSOCIATED(field% q_phy   )) THEN
         !$ACC DATA PRESENT( field%q_phy )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             field% q_phy(jc,jk,jb) = field% q_phy(jc,jk,jb) + q_gwd(jc,jk)
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
           field% q_phy_vi(jc,jb) = field% q_phy_vi(jc,jb) + SUM(q_gwd(jc,:))
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_gwd)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          !$ACC DATA PRESENT( tend%ta_phy, tend%ua_phy, tend%va_phy )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              ! use tendency to update the model state
              tend% ta_phy(jc,jk,jb) = tend% ta_phy(jc,jk,jb) + tend_ta_gwd(jc,jk)
              tend% ua_phy(jc,jk,jb) = tend% ua_phy(jc,jk,jb) + tend_ua_gwd(jc,jk)
              tend% va_phy(jc,jk,jb) = tend% va_phy(jc,jk,jb) + tend_va_gwd(jc,jk)
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
       SELECT CASE(fc_gwd)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1,2)
          ! use tendency to update the physics state
          IF (lparamcpl) THEN
             !$ACC DATA PRESENT( field%ta, field%ua, field%va )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP SEQ
             DO jk = 1, nlev
               !$ACC LOOP GANG VECTOR
               DO jc = jcs, jce
                 field% ta(jc,jk,jb) = field% ta(jc,jk,jb) + tend_ta_gwd(jc,jk)*pdtime
                 field% ua(jc,jk,jb) = field% ua(jc,jk,jb) + tend_ua_gwd(jc,jk)*pdtime
                 field% va(jc,jk,jb) = field% va(jc,jk,jb) + tend_va_gwd(jc,jk)*pdtime
               END DO
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
          END IF
       END SELECT
       !
       !$ACC END DATA
       !
    ELSE
       !
       IF (ASSOCIATED(field% q_gwd)) THEN
         !$ACC DATA PRESENT( field%q_gwd )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             field% q_gwd(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_gwd_vi)) THEN
         !$ACC DATA PRESENT( field%q_gwd_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_gwd_vi(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% ta_gwd)) THEN
         !$ACC DATA PRESENT( tend%ta_gwd )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% ta_gwd(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ua_gwd)) THEN
         !$ACC DATA PRESENT( tend%ua_gwd )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
         tend% ua_gwd(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% va_gwd)) THEN
         !$ACC DATA PRESENT( tend%va_gwd )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% va_gwd(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
    END IF

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_gwd_output(jg, jb, jcs, jce, nproma, nlev, field, tend)

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend )

    IF (ltimer) call timer_stop(timer_gwd)

  END SUBROUTINE interface_echam_gwd

END MODULE mo_interface_echam_gwd
