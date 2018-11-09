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
    LOGICAL                 ,POINTER    :: lparamcpl
    INTEGER                 ,POINTER    :: fc_mox
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp)                            :: tend_qtrc_mox(nproma,nlev,ntracer)

    IF (ltimer) call timer_start(timer_mox)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_mox    => echam_phy_config(jg)%fc_mox
    field     => prm_field(jg)
    tend      => prm_tend (jg)

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
             tend% qtrc_mox(jcs:jce,:,jb,iqv) = tend_qtrc_mox(jcs:jce,:,iqv)
          END IF
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(tend% qtrc_mox)) THEN
             tend_qtrc_mox(jcs:jce,:,iqv) = tend% qtrc_mox(jcs:jce,:,jb,iqv)
          END IF
          !
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_mox)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend% qtrc_phy(jcs:jce,:,jb,iqv) = tend% qtrc_phy(jcs:jce,:,jb,iqv) + tend_qtrc_mox(jcs:jce,:,iqv)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field% qtrc(jcs:jce,:,jb,iqv)  = field% qtrc(jcs:jce,:,jb,iqv) + tend_qtrc_mox(jcs:jce,:,iqv)*pdtime
       END IF
       !
    ELSE
       !
       IF (ASSOCIATED(tend% qtrc_mox)) THEN
          tend% qtrc_mox(jcs:jce,:,jb,iqv) = 0.0_wp
       END IF
       !
    END IF

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_mox)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_mox)

  END SUBROUTINE interface_echam_mox

END MODULE mo_interface_echam_mox
