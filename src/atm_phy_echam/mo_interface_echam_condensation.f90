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

MODULE mo_interface_echam_condensation

  USE mo_kind                ,ONLY: wp

  USE mo_model_domain        ,ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev, iqv, iqc, iqi

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, t_echam_phy_tend

  USE mo_cloud               ,ONLY: cloud

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cloud

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_condensation

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_condensation(is_in_sd_ed_interval,    &
       &                                  is_active,               &
       &                                  patch, rl_start, rl_end, &
       &                                  jks,                     &
       &                                  field, tend,             &
       &                                  ictop,                   &
       &                                  pdtime                   )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch
    INTEGER                 ,INTENT(in)    :: rl_start, rl_end
    INTEGER                 ,INTENT(in)    :: jks                                !< start vertical level
    TYPE(t_echam_phy_field) ,POINTER       :: field    
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    INTEGER                 ,INTENT(in)    :: ictop (:,:)                        !< from massflux
    REAL(wp)                ,INTENT(in)    :: pdtime

    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (ltimer) CALL timer_start(timer_cloud)
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       !
       CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       !
       CALL echam_condensation(is_in_sd_ed_interval,          &
            &                  is_active,                     &
            &                  jb,jcs,jce, nproma,            &
            &                  jks,                           &
            &                  field, tend,                   &
            &                  ictop(:,jb),                   &
            &                  pdtime                         )
       !
    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_cloud)

  END SUBROUTINE interface_echam_condensation
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE echam_condensation(is_in_sd_ed_interval,  &
       &                        is_active,             &
       &                        jb,jcs,jce, nbdim,     &
       &                        jks,                   &
       &                        field, tend,           &
       &                        ictop,                 &
       &                        pdtime                 )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    INTEGER                 ,INTENT(in)    :: jb                  !< block index
    INTEGER                 ,INTENT(in)    :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in)    :: nbdim               !< size of this block 
    INTEGER                 ,INTENT(in)    :: jks                 !< start vertical level
    TYPE(t_echam_phy_field) ,POINTER       :: field
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    INTEGER                 ,INTENT(in)    :: ictop  (nbdim)      !< from massflux
    REAL(wp)                ,INTENT(in)    :: pdtime

    ! local
    INTEGER  :: itype(nbdim)              !< type of convection

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          itype(:) = NINT(field%rtype(:,jb))
          !
          CALL cloud(jce, nbdim, jks, nlev,           &! in
               &     pdtime,                       &! in
               &     ictop,                        &! in (from "cucall")
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
               &      tend% ta       (:,:,jb),     &! in  tte
               &      tend% qtrc     (:,:,jb,iqv), &! in  qte
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
       tend% ta  (jcs:jce,:,jb)      = tend% ta  (jcs:jce,:,jb)     + tend% ta_cld  (jcs:jce,:,jb)
       tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv) + tend% qtrc_cld(jcs:jce,:,jb,iqv)
       tend% qtrc(jcs:jce,:,jb,iqc)  = tend% qtrc(jcs:jce,:,jb,iqc) + tend% qtrc_cld(jcs:jce,:,jb,iqc)
       tend% qtrc(jcs:jce,:,jb,iqi)  = tend% qtrc(jcs:jce,:,jb,iqi) + tend% qtrc_cld(jcs:jce,:,jb,iqi)
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

  END SUBROUTINE echam_condensation
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_condensation
