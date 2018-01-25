!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
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

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev, iqv, iqc, iqi

  USE mtime                  ,ONLY: datetime
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cloud

  USE mo_cloud               ,ONLY: cloud

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_cld

CONTAINS

  SUBROUTINE echam_cld(is_in_sd_ed_interval, &
       &               is_active,            &
       &               jg,jb,jcs,jce,        &
       &               datetime_old,         &
       &               pdtime                )

    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    INTEGER                 ,INTENT(in) :: jg                  !< grid  index
    INTEGER                 ,INTENT(in) :: jb                  !< block index
    INTEGER                 ,INTENT(in) :: jcs, jce            !< start/end column index within this block
    TYPE(datetime)          ,POINTER    :: datetime_old        !< generic input, not used in echam_sso
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    !
    INTEGER  :: itype(nproma)    !< type of convection

    IF (ltimer) call timer_start(timer_cloud)

    ! associate pointers
    field => prm_field(jg)
    tend  => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          itype(:) = NINT(field%rtype(:,jb))
          !
          CALL cloud(jg,                           &! in
               &     jcs, jce, nproma, nlev,       &! in
               &     pdtime,                       &! in
               &     field% ictop    (:,  jb),     &! in (from "cucall")
               &     field% presm_old(:,:,jb),     &! in
               &     field% dz       (:,:,jb),     &! in
               &     field% mdry     (:,:,jb),     &! in
               &     field% rho      (:,:,jb),     &! in
               &     field% cpair    (:,:,jb),     &! in
               &     field% acdnc    (:,:,jb),     &! in  acdnc
               &     field% ta       (:,:,jb),     &! in  tm1
               &     field% qtrc     (:,:,jb,iqv), &! in  qm1
               &     field% qtrc     (:,:,jb,iqc), &! in  xlm1
               &     field% qtrc     (:,:,jb,iqi), &! in  xim1
               &      tend%   ta_phy (:,:,jb),     &! in  tte
               &      tend% qtrc_phy (:,:,jb,iqv), &! in  qte
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
       ! accumulate tendencies
       tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)     + tend%   ta_cld(jcs:jce,:,jb)
       tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv) + tend% qtrc_cld(jcs:jce,:,jb,iqv)
       tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc) + tend% qtrc_cld(jcs:jce,:,jb,iqc)
       tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi) + tend% qtrc_cld(jcs:jce,:,jb,iqi)
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

    IF (ltimer) call timer_stop(timer_cloud)

  END SUBROUTINE echam_cld
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_cld
