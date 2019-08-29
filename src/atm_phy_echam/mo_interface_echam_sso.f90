!>
!! @brief Subroutine interface_echam_sso calls the SSO parameterization scheme.
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

MODULE mo_interface_echam_sso

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_sso

  USE mo_echam_sso_config    ,ONLY: echam_sso_config
  USE mo_ssodrag             ,ONLY: ssodrag
  !$ser verbatim USE mo_ser_echam_sso, ONLY: serialize_sso_input,&
  !$ser verbatim                             serialize_sso_output
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_sso

CONTAINS

  SUBROUTINE interface_echam_sso(jg, jb,jcs,jce       ,&
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
    INTEGER                             :: fc_sso
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    LOGICAL                             :: lsftlf   ! <-- provisional

    ! Local variables
    !
    REAL(wp)                            ::    u_stress_sso(nproma)
    REAL(wp)                            ::    v_stress_sso(nproma)
    REAL(wp)                            :: dissipation_sso(nproma)
    !
    REAL(wp)                            :: q_sso(nproma,nlev)
    !
    REAL(wp)                            :: tend_ta_sso(nproma,nlev)
    REAL(wp)                            :: tend_ua_sso(nproma,nlev)
    REAL(wp)                            :: tend_va_sso(nproma,nlev)
    !
    REAL(wp)                            :: zdis_sso(nproma,nlev)  !<  out, energy dissipation rate [J/s/kg]
    REAL(wp)                            :: zscale(nproma)         !< area scaling factor
    INTEGER                             :: jc, jk

    IF (ltimer) call timer_start(timer_sso)

    ! associate pointers
    lparamcpl = echam_phy_config(jg)%lparamcpl
    fc_sso    = echam_phy_config(jg)%fc_sso
    field     => prm_field(jg)
    tend      => prm_tend (jg)
    lsftlf    = echam_sso_config(jg)%lsftlf        ! <-- provisional

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_sso_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    IF ( is_in_sd_ed_interval ) THEN
       !
       !$ACC DATA PRESENT( field%qconv ) &
       !$ACC       CREATE( u_stress_sso, v_stress_sso, dissipation_sso, q_sso, tend_ta_sso,  &
       !$ACC               tend_ua_sso, tend_va_sso, zdis_sso, zscale )
       !

       IF ( is_active ) THEN
          !
          ! area scaling factor                     ! <-- provisional
          IF (lsftlf) THEN                          ! <-- provisional
             !$ACC DATA PRESENT( field%sftlf )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG VECTOR
             DO jc = 1, nproma
               zscale(jc) = field% sftlf (jc,jb)
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
          ELSE                                      ! <-- provisional
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG VECTOR
             DO jc = 1, nproma
               zscale(jc) = 1.0_wp                     ! <-- provisional
             END DO
             !$ACC END PARALLEL
          END IF                                    ! <-- provisional
          !
          CALL ssodrag(jg                           ,& ! in,  grid index
               &       jcs, jce                     ,& ! in,  start and end index
               &       nproma                       ,& ! in,  dimension of block of cells/columns
               &       nlev                         ,& ! in,  number of levels
               !
               &       pdtime                       ,& ! in,  time step length
               &       field% coriol(:,jb)          ,& ! in,  Coriolis parameter (1/s)
               &       field% zf  (:,:,jb)          ,& ! in,  full level height (m)
               &       field% zh  (:,nlev+1,jb)     ,& ! in,  surface height    (m)
               !
               &       field% presi_old(:,:,jb)     ,& ! in,  p at half levels
               &       field% presm_old(:,:,jb)     ,& ! in,  p at full levels
               &       field% mair(:,:,jb)          ,& ! in,  air mass
               &       field%   ta(:,:,jb)          ,& ! in,  T
               &       field%   ua(:,:,jb)          ,& ! in,  u
               &       field%   va(:,:,jb)          ,& ! in,  v
               !
               &       field% oromea(:,jb)          ,& ! in,  Mean Orography (m)
               &       field% orostd(:,jb)          ,& ! in,  SSO standard deviation (m)
               &       field% orosig(:,jb)          ,& ! in,  SSO slope
               &       field% orogam(:,jb)          ,& ! in,  SSO Anisotropy
               &       field% orothe(:,jb)          ,& ! in,  SSO Angle
               &       field% oropic(:,jb)          ,& ! in,  SSO Peaks elevation (m)
               &       field% oroval(:,jb)          ,& ! in,  SSO Valleys elevation (m)
               !
               &              zscale(:)             ,& ! in,  area fraction of land incl. lakes
               !                                              where the SSO params are valid
               !
               &          u_stress_sso(:)           ,& ! out, u-gravity wave stress
               &          v_stress_sso(:)           ,& ! out, v-gravity wave stress
               &       dissipation_sso(:)           ,& ! out, dissipation by gravity wave drag
               !
               &          zdis_sso(:,:)             ,& ! out, energy dissipation rate
               &       tend_ua_sso(:,:)             ,& ! out, tendency of zonal wind
               &       tend_va_sso(:,:)              ) ! out, tendency of meridional wind
          !
          ! heating
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              q_sso(jc,jk) = zdis_sso(jc,jk) * field%mair(jc,jk,jb)
            END DO
          END DO
          !$ACC END PARALLEL
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% u_stress_sso)) THEN
            !$ACC DATA PRESENT( field%u_stress_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% u_stress_sso(jc,jb) = u_stress_sso(jc)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% v_stress_sso)) THEN
            !$ACC DATA PRESENT( field%v_stress_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% v_stress_sso(jc,jb) = v_stress_sso(jc)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% dissipation_sso)) THEN
            !$ACC DATA PRESENT( field%dissipation_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% dissipation_sso(jc,jb) = dissipation_sso(jc)
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(field% q_sso)) THEN
            !$ACC DATA PRESENT( field%q_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                field% q_sso(jc,jk,jb) = q_sso(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_sso_vi)) THEN
            !$ACC DATA PRESENT( field%q_sso_vi )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_sso_vi(jc,jb) = SUM(q_sso(jc,:))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_sso)) THEN
            !$ACC DATA PRESENT( tend%ua_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend% ua_sso(jc,jk,jb) = tend_ua_sso(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_sso)) THEN
            !$ACC DATA PRESENT( tend%va_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend% va_sso(jc,jk,jb) = tend_va_sso(jc,jk)
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
          IF (ASSOCIATED(field% q_sso)) THEN
            !$ACC DATA PRESENT( field%q_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                q_sso(jc,jk) = field% q_sso(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          IF (ASSOCIATED(tend% ua_sso)) THEN
            !$ACC DATA PRESENT( tend%ua_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend_ua_sso(jc,jk) = tend% ua_sso(jc,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(tend% va_sso)) THEN
            !$ACC DATA PRESENT( tend%va_sso )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              !$ACC LOOP GANG VECTOR
              DO jc = jcs, jce
                tend_va_sso(jc,jk) = tend% va_sso(jc,jk,jb)
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
           tend_ta_sso(jc,jk) = q_sso(jc,jk) * field% qconv(jc,jk,jb)
         END DO
       END DO
       !$ACC END PARALLEL
       !
       IF (ASSOCIATED(tend% ta_sso)) THEN
         !$ACC DATA PRESENT( tend%ta_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% ta_sso(jc,jk,jb) = tend_ta_sso(jc,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF

       ! for output: accumulate heating
       IF (ASSOCIATED(field% q_phy)) THEN
         !$ACC DATA PRESENT( field%q_phy )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             field% q_phy (jc,jk,jb) = field% q_phy(jc,jk,jb) + q_sso(jc,jk)
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
         field% q_phy_vi(jc,jb) = field% q_phy_vi(jc,jb) + SUM(q_sso(jc,:))
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_sso)
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
              tend% ta_phy(jc,jk,jb) = tend% ta_phy(jc,jk,jb) + tend_ta_sso(jc,jk)
              tend% ua_phy(jc,jk,jb) = tend% ua_phy(jc,jk,jb) + tend_ua_sso(jc,jk)
              tend% va_phy(jc,jk,jb) = tend% va_phy(jc,jk,jb) + tend_va_sso(jc,jk)
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
          !$ACC DATA PRESENT( field%ta, field%ua, field%va )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% ta(jc,jk,jb) = field% ta(jc,jk,jb) + tend_ta_sso(jc,jk)*pdtime
              field% ua(jc,jk,jb) = field% ua(jc,jk,jb) + tend_ua_sso(jc,jk)*pdtime
              field% va(jc,jk,jb) = field% va(jc,jk,jb) + tend_va_sso(jc,jk)*pdtime
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF
       !
       !$ACC END DATA
       !
    ELSE
       !
       IF (ASSOCIATED(field% u_stress_sso)) THEN
         !$ACC DATA PRESENT( field%u_stress_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% u_stress_sso(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% v_stress_sso)) THEN
         !$ACC DATA PRESENT( field%v_stress_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% v_stress_sso(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% dissipation_sso)) THEN
         !$ACC DATA PRESENT( field%dissipation_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% dissipation_sso(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(field% q_sso)) THEN
         !$ACC DATA PRESENT( field%dissipation_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             field% q_sso(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_sso_vi)) THEN
         !$ACC DATA PRESENT( field%q_sso_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_sso_vi(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% ta_sso)) THEN
         !$ACC DATA PRESENT( tend%ta_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% ta_sso(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ua_sso)) THEN
         !$ACC DATA PRESENT( tend%ua_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% ua_sso(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% va_sso)) THEN
         !$ACC DATA PRESENT( tend%va_sso )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP SEQ
         DO jk = 1, nlev
           !$ACC LOOP GANG VECTOR
           DO jc = jcs, jce
             tend% va_sso(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
    END IF

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_sso_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
    
    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_sso)

  END SUBROUTINE interface_echam_sso

END MODULE mo_interface_echam_sso
