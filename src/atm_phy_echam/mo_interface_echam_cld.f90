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

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_cld

CONTAINS

  SUBROUTINE interface_echam_cld(jg,jb,jcs,jce        ,&
       &                         nproma,nlev          ,& 
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         pdtime               )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev
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
    INTEGER  :: itype(nproma)    !< type of convection

    IF (ltimer) call timer_start(timer_cld)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_cld    => echam_phy_config(jg)%fc_cld
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          itype(:) = NINT(field%rtype(:,jb))
          !
          CALL cloud(jg,                           &! in
               &     jce, nproma, nlev,            &! in
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
               &     field% relhum   (:,:,jb),     &! out
               &     field% q_cld    (:,:,jb),     &! out
               &      tend% qtrc_cld (:,:,jb,iqv), &! out
               &      tend% qtrc_cld (:,:,jb,iqc), &! out
               &      tend% qtrc_cld (:,:,jb,iqi)  )! out
          !
          field% rtype(:,jb) = REAL(itype(:),wp)
          !
          ! vertical integral
          field% q_cld_vi(jcs:jce,jb) = SUM(field% q_cld(jcs:jce,:,jb),DIM=2)
          !
       END IF
       !
       ! convert    heating
       tend% ta_cld(jcs:jce,:,jb) = field% q_cld(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field% q_cld(jcs:jce,:,jb)
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_cld)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)     + tend%   ta_cld(jcs:jce,:,jb)
          tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv) + tend% qtrc_cld(jcs:jce,:,jb,iqv)
          tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc) + tend% qtrc_cld(jcs:jce,:,jb,iqc)
          tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi) + tend% qtrc_cld(jcs:jce,:,jb,iqi)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field%   ta(jcs:jce,:,jb)      = field%   ta(jcs:jce,:,jb)      + tend%   ta_cld(jcs:jce,:,jb)    *pdtime
          field% qtrc(jcs:jce,:,jb,iqv)  = field% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_cld(jcs:jce,:,jb,iqv)*pdtime
          field% qtrc(jcs:jce,:,jb,iqc)  = field% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_cld(jcs:jce,:,jb,iqc)*pdtime
          field% qtrc(jcs:jce,:,jb,iqi)  = field% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_cld(jcs:jce,:,jb,iqi)*pdtime
       END IF
       !
    ELSE
       !
       tend% ta_cld(jcs:jce,:,jb) = 0.0_wp
       !
       tend% qtrc_cld(jcs:jce,:,jb,iqv) = 0.0_wp
       tend% qtrc_cld(jcs:jce,:,jb,iqc) = 0.0_wp
       tend% qtrc_cld(jcs:jce,:,jb,iqi) = 0.0_wp
       !
    END IF

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_cld)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_cld)

  END SUBROUTINE interface_echam_cld

END MODULE mo_interface_echam_cld
