!>
!! @brief Subroutine interface_echam_mox calls the CH4 oxidation and H2O photolysis schemes.
!!
!! @author S. Rast, MPI-M
!!
!! @par Revision History
!!   original source S. Rast, MPI-M (03-24-2017)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_mox
  
  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_mox

  USE mo_run_config          ,ONLY: iqv
  USE mo_methox              ,ONLY: methox
  !$ser verbatim USE mo_ser_echam_mox, ONLY: serialize_mox_input,&
  !$ser verbatim                             serialize_mox_output

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_mox

CONTAINS
  
  SUBROUTINE interface_echam_mox(jg, jb,jcs,jce       ,&
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
    INTEGER                             :: fc_mox
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp)                            :: tend_qtrc_mox(nproma,nlev,ntracer)
    INTEGER                             :: jc, jk

    IF (ltimer) call timer_start(timer_mox)

    ! associate pointers
    lparamcpl = echam_phy_config(jg)%lparamcpl
    fc_mox    = echam_phy_config(jg)%fc_mox
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    ! Serialbox2 input fields serialization
    !$ser verbatim call serialize_mox_input(jg, jb, jcs, jce, nproma, nlev, field, tend)

    !$ACC DATA CREATE( tend_qtrc_mox )

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          CALL methox( jcs, jce,                  &
               &       nproma, nlev,              &
               &       field% presm_old(:,:,jb),  &
               &       field% qtrc(:,:,jb,iqv),   &
               &       tend_qtrc_mox(:,:,iqv)     )
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(tend% qtrc_mox)) THEN
             !$ACC DATA PRESENT( tend%qtrc_mox)
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP SEQ
             DO jk = 1, nlev
               !$ACC LOOP GANG VECTOR
               DO jc = jcs, jce
                 tend% qtrc_mox(jc,jk,jb,iqv) = tend_qtrc_mox(jc,jk,iqv)
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
          IF (ASSOCIATED(tend% qtrc_mox)) THEN
             !$ACC DATA PRESENT( tend%qtrc_mox )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP SEQ
             DO jk = 1, nlev
             !$ACC LOOP GANG VECTOR
               DO jc= jcs, jce
                 tend_qtrc_mox(jc,jk,iqv) = tend% qtrc_mox(jc,jk,jb,iqv)
               END DO
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
          END IF
          !
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_mox)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          !$ACC DATA PRESENT( tend%qtrc_phy )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
            ! use tendency to update the model state
              tend% qtrc_phy(jc,jk,jb,iqv) = tend% qtrc_phy(jc,jk,jb,iqv) + tend_qtrc_mox(jc,jk,iqv)
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
          !$ACC DATA PRESENT( field%qtrc )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              field% qtrc(jc,jk,jb,iqv)  = field% qtrc(jc,jk,jb,iqv) + tend_qtrc_mox(jc,jk,iqv)*pdtime
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF
       !
    ELSE
       !
       IF (ASSOCIATED(tend% qtrc_mox)) THEN
          !$ACC DATA PRESENT( tend%qtrc_mox )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR
            DO jc = jcs, jce
              tend% qtrc_mox(jc,jk,jb,iqv) = 0.0_wp
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF
       !
    END IF

    !$ACC END DATA
    
    ! Serialbox2 output fields serialization
    !$ser verbatim call serialize_mox_output(jg, jb, jcs, jce, nproma, nlev, field, tend)

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_mox)

  END SUBROUTINE interface_echam_mox

END MODULE mo_interface_echam_mox
