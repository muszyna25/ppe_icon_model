!>
!! @brief Subroutine interface_echam_rht calls the radiative heating scheme.
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

MODULE mo_interface_echam_rht

  USE mo_kind                   ,ONLY: wp
  USE mo_exception              ,ONLY: finish
  USE mtime                     ,ONLY: datetime

  USE mo_echam_phy_config       ,ONLY: echam_phy_config
  USE mo_echam_phy_memory       ,ONLY: t_echam_phy_field, prm_field, &
    &                                  t_echam_phy_tend,  prm_tend

  USE mo_timer                  ,ONLY: ltimer, timer_start, timer_stop, timer_rht

  USE mo_radheating             ,ONLY: radheating
  USE mo_psrad_solar_parameters ,ONLY: psctm
  !$ser verbatim USE mo_ser_echam_rht, ONLY: serialize_rht_input,&
  !$ser verbatim                             serialize_rht_output

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_rht

CONTAINS

  SUBROUTINE interface_echam_rht(jg, jb,jcs,jce       ,&
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
    INTEGER                 ,POINTER    :: fc_rht
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER                             :: nlevp1, jc, jk
    !
    REAL(wp)                            :: q_rad(nproma,nlev)
    REAL(wp)                            :: q_rlw(nproma,nlev)
    REAL(wp)                            :: q_rsw(nproma,nlev)
    !
    REAL(wp)                            :: tend_ta_rad(nproma,nlev)

    IF (ltimer) CALL timer_start(timer_rht)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_rht    => echam_phy_config(jg)%fc_rht
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_rht_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    IF ( is_in_sd_ed_interval ) THEN
       !$ACC DATA PRESENT( field%cosmu0, field%daylght_frc, field%ts_rad, field%ts_rad_rt, &
       !$ACC               field%rsd_rt, field%rsu_rt, field%rsdcs_rt, field%rsucs_rt,     &
       !$ACC               field%rld_rt, field%rlu_rt, field%rldcs_rt, field%rlucs_rt,     &
       !$ACC               field%rvds_dir_rt, field%rpds_dir_rt, field%rnds_dir_rt,        &
       !$ACC               field%rvds_dif_rt, field%rpds_dif_rt, field%rnds_dif_rt,        &
       !$ACC               field%rvus_rt, field%rpus_rt, field%rnus_rt, field%rsdt,        &
       !$ACC               field%rsut, field%rsds, field%rsus, field%rsutcs, field%rsdscs, &
       !$ACC               field%rsuscs, field%rvds_dir, field%rpds_dir, field%rnds_dir,   &
       !$ACC               field%rvds_dif, field%rpds_dif, field%rnds_dif, field%rvus,     &
       !$ACC               field%rpus, field%rnus, field%rlut, field%rlds, field%rlus,     &
       !$ACC               field%rlutcs, field%rldscs, field%q_rlw_nlev, field%qconv,      &
       !$ACC               field%emissivity )                                              &
       !$ACC       CREATE( q_rad, q_rsw, q_rlw, tend_ta_rad )
       !
       IF (is_active) THEN
          !
          nlevp1 = nlev+1
          !
          CALL radheating (                                   &
               !
               ! input
               ! -----
               !
               & jcs        = jcs                            ,&! loop start index
               & jce        = jce                            ,&! loop end index
               & kbdim      = nproma                         ,&! dimension size
               & klev       = nlev                           ,&! vertical dimension size
               & klevp1     = nlevp1                         ,&! vertical dimension size
               !
               & rsdt0      = psctm(jg)                      ,&! toa incident shortwave radiation for sun in zenith
               & cosmu0     = field%cosmu0    (:,jb)         ,&! solar zenith angle at current time
               & daylght_frc= field%daylght_frc(:,jb)        ,&! daylight fraction
               !
               & emiss      = field%emissivity (:,jb)        ,&! lw sfc emissivity
               & tsr        = field%ts_rad (:,jb)            ,&! radiative surface temperature at current   time [K]
               & tsr_rt     = field%ts_rad_rt(:,jb)          ,&! radiative surface temperature at radiation time [K]
               !
               & rsd_rt     = field%rsd_rt           (:,:,jb),&! all-sky   shortwave downward flux at radiation time [W/m2]
               & rsu_rt     = field%rsu_rt           (:,:,jb),&! all-sky   shortwave upward   flux at radiation time [W/m2]
               !
               & rsdcs_rt   = field%rsdcs_rt         (:,:,jb),&! clear-sky shortwave downward flux at radiation time [W/m2]
               & rsucs_rt   = field%rsucs_rt         (:,:,jb),&! clear-sky shortwave upward   flux at radiation time [W/m2]
               !
               & rld_rt     = field%rld_rt           (:,:,jb),&! all-sky   longwave  downward flux at radiation time [W/m2]
               & rlu_rt     = field%rlu_rt           (:,:,jb),&! all-sky   longwave  upward   flux at radiation time [W/m2]
               !
               & rldcs_rt   = field%rldcs_rt         (:,:,jb),&! clear-sky longwave  downward flux at radiation time [W/m2]
               & rlucs_rt   = field%rlucs_rt         (:,:,jb),&! clear-sky longwave  upward   flux at radiation time [W/m2]
               !
               & rvds_dir_rt= field%rvds_dir_rt        (:,jb),&!< out  all-sky downward direct visible radiation at surface
               & rpds_dir_rt= field%rpds_dir_rt        (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
               & rnds_dir_rt= field%rnds_dir_rt        (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
               & rvds_dif_rt= field%rvds_dif_rt        (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
               & rpds_dif_rt= field%rpds_dif_rt        (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
               & rnds_dif_rt= field%rnds_dif_rt        (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
               & rvus_rt    = field%rvus_rt            (:,jb),&!< out  all-sky upward visible radiation at surface
               & rpus_rt    = field%rpus_rt            (:,jb),&!< out  all-sky upward PAR     radiation at surfac
               & rnus_rt    = field%rnus_rt            (:,jb),&!< out  all-sky upward near-IR radiation at surface
               !
               ! output
               ! ------
               !
               & rsdt       = field%rsdt               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
               & rsut       = field%rsut               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
               & rsds       = field%rsds               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
               & rsus       = field%rsus               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
               !
               & rsutcs     = field%rsutcs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
               & rsdscs     = field%rsdscs             (:,jb),&! clear-sky shortwave downward flux at current   time [W/m2]
               & rsuscs     = field%rsuscs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
               !
               & rvds_dir   = field%rvds_dir           (:,jb),&!< out  all-sky downward direct visible radiation at surface
               & rpds_dir   = field%rpds_dir           (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
               & rnds_dir   = field%rnds_dir           (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
               & rvds_dif   = field%rvds_dif           (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
               & rpds_dif   = field%rpds_dif           (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
               & rnds_dif   = field%rnds_dif           (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
               & rvus       = field%rvus               (:,jb),&!< out  all-sky upward visible radiation at surface
               & rpus       = field%rpus               (:,jb),&!< out  all-sky upward PAR     radiation at surfac
               & rnus       = field%rnus               (:,jb),&!< out  all-sky upward near-IR radiation at surface
               !
               & rlut       = field%rlut               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
               & rlds       = field%rlds               (:,jb),&! all-sky   longwave  downward flux at current   time [W/m2]
               & rlus       = field%rlus               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
               !
               & rlutcs     = field%rlutcs             (:,jb),&! clear-sky longwave  upward   flux at current   time [W/m2]
               & rldscs     = field%rldscs             (:,jb),&! clear-sky longwave  downward flux at current   time [W/m2]
               !
               & q_rsw      = q_rsw                    (:,:) ,&! rad. heating by SW           [W/m2]
               & q_rlw      = q_rlw                    (:,:)  )! rad. heating by LW           [W/m2]
          !
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = jcs, jce
              q_rad(jc,jk) = q_rsw(jc,jk)+q_rlw(jc,jk) ! rad. heating by SW+LW        [W/m2]
            END DO
          END DO
          !$ACC END PARALLEL
          !
          ! for output: SW+LW heating
          IF (ASSOCIATED(field% q_rad)) THEN
            !$ACC DATA PRESENT( field%q_rad )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field% q_rad(jc,jk,jb) = q_rad(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_rad_vi)) THEN
            !$ACC DATA PRESENT( field%q_rad_vi )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_rad_vi(jc, jb) = SUM(q_rad(jc,:))
            END DO
            !$acc END PARALLEL
            !$acc END DATA
          END IF
          !
          ! for output: SW heating
          IF (ASSOCIATED(field% q_rsw)) THEN
            !$ACC DATA PRESENT( field%q_rsw )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field% q_rsw(jc,jk,jb) = q_rsw(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_rsw_vi)) THEN
            !$ACC DATA PRESENT( field%q_rsw_vi )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_rsw_vi(jc,jb) = SUM(q_rsw(jc,:))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          !
          ! for output: LW heating
          IF (ASSOCIATED(field% q_rlw)) THEN 
            !$ACC DATA PRESENT( field%q_rlw )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field% q_rlw(jc,jk,jb) = q_rlw(jc,jk)
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          IF (ASSOCIATED(field% q_rlw_vi)) THEN
            !$ACC DATA PRESENT( field%q_rlw_vi )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% q_rlw_vi(jc,jb) = SUM(q_rlw(jc,:))
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
          
          ! store LW heating in lowermost layer separately,
          ! which is needed for computing q_rlw_impl
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            field% q_rlw_nlev(jc,jb) = q_rlw(jc,nlev)
          END DO
          !$ACC END PARALLEL
          !
       ELSE
          CALL finish('mo_interface_echam_rht','interface_echam_rht must not be called with is_active=.FALSE.')
       END IF
       !
       ! convert    heating
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = 1, nlev
         !$ACC LOOP VECTOR
         DO jc = jcs, jce
           tend_ta_rad(jc,jk) = q_rad(jc,jk) * field% qconv(jc,jk,jb)
         END DO
       END DO
       !$ACC END PARALLEL
       !
       IF (ASSOCIATED(tend% ta_rad)) THEN
         !$ACC DATA PRESENT( tend%ta_rad )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_rad(jc,jk,jb) = tend_ta_rad(jc,jk)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ta_rsw)) THEN
         !$ACC DATA PRESENT( tend%ta_rsw )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_rsw(jc,jk,jb) = q_rsw(jc,jk) * field% qconv(jc,jk,jb)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ta_rlw)) THEN
         !$ACC DATA PRESENT( tend%ta_rlw )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_rlw(jc,jk,jb) = q_rlw(jc,jk) * field% qconv(jc,jk,jb)
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate heating
       IF (ASSOCIATED(field% q_phy   )) THEN
         !$ACC DATA PRESENT( field%q_phy )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_phy(jc,jk,jb) = field% q_phy(jc,jk,jb) + q_rad(jc,jk)
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
           field% q_phy_vi(jc,jb) = field% q_phy_vi(jc,jb) + SUM(q_rad(jc,:))
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_rht)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
         !$ACC DATA PRESENT( tend%ta_phy )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_phy(jc,jk,jb) = tend% ta_phy(jc,jk,jb) + tend_ta_rad (jc,jk)
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
       SELECT CASE(fc_rht)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1,2)
          ! use tendency to update the physics state
          IF (lparamcpl) THEN
            !$ACC DATA PRESENT( tend%ta )
            !$ACC PARALLEL DEFAULT(PRESENT)
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = jcs, jce
                field% ta(jc,jk,jb) = field% ta(jc,jk,jb) + tend_ta_rad(jc,jk)*pdtime
              END DO
            END DO
            !$ACC END PARALLEL
            !$ACC END DATA
          END IF
       END SELECT
       !$ACC END DATA
       !
    ELSE
       !
       IF (ASSOCIATED(field% q_rad)) THEN
         !$ACC DATA PRESENT( field%q_rad )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_rad(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_rad_vi)) THEN
         !$ACC DATA PRESENT( field%q_rad_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_rad_vi(jc, jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(field% q_rsw)) THEN
         !$ACC DATA PRESENT( field%q_rsw )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_rsw(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_rsw_vi)) THEN
         !$ACC DATA PRESENT( field%q_rsw_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_rsw_vi(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(field% q_rlw)) THEN
         !$ACC DATA PRESENT( field%q_rlw )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             field% q_rlw(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(field% q_rlw_vi)) THEN
         !$ACC DATA PRESENT( field%q_rlw_vi )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG VECTOR
         DO jc = jcs, jce
           field% q_rlw_vi(jc,jb) = 0.0_wp
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% ta_rad)) THEN
         !$ACC DATA PRESENT( tend%ta_rad )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_rad(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ta_rsw)) THEN
         !$ACC DATA PRESENT( tend%ta_rsw )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_rsw(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       IF (ASSOCIATED(tend% ta_rlw)) THEN
         !$ACC DATA PRESENT( tend%ta_rlw )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jc = jcs, jce
             tend% ta_rlw(jc,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       !$ACC DATA PRESENT( field%q_rlw_nlev )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP VECTOR
       DO jc = jcs, jce
         field% q_rlw_nlev(jc,jb) = 0.0_wp
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
    END IF

    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_rht_output(jg, jb, jcs, jce, nproma, nlev, field, tend)

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_rht)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_rht)

  END SUBROUTINE interface_echam_rht

END MODULE mo_interface_echam_rht
