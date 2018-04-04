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

#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_interface_echam_cnv

  USE mo_kind                ,ONLY: wp

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev, nlevm1, nlevp1, iqv, iqc, iqi, iqt, ntracer

  USE mtime                  ,ONLY: datetime
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_convection

  USE mo_cumastr             ,ONLY: cumastr

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_cnv

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE echam_cnv(is_in_sd_ed_interval, &
       &               is_active,            &
       &               jg,jb,jcs,jce,        &
       &               datetime_old,         &
       &               pdtime                )

    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    INTEGER                 ,INTENT(in) :: jg                  !< grid  index
    INTEGER                 ,INTENT(in) :: jb                  !< block index
    INTEGER                 ,INTENT(in) :: jcs, jce            !< start/end column index within this block
    TYPE(datetime)          ,POINTER    :: datetime_old        !< generic input, not used in echam_cnv
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    ! 
    INTEGER  :: itype(nproma)                !< type of convection
    INTEGER  :: ntrac                        !< # of tracers excluding water vapour and hydrometeors
                                             !< which are convectively transported


    REAL(wp) :: ztop(nproma)                 !< convective cloud top pressure   [Pa]
    REAL(wp) :: zta    (nproma,nlev)         !< provisional temperature         [K]
    REAL(wp) :: zqtrc  (nproma,nlev,ntracer) !< provisional mass mixing ratios  [kg/kg]
    REAL(wp) :: zua    (nproma,nlev)         !< provisional zonal      wind     [m/s]
    REAL(wp) :: zva    (nproma,nlev)         !< provisional meridional wind     [m/s]
    REAL(wp) :: zqtrc_cnd(nproma,nlev)       !< cloud condensate mixing ratio   [kg/kg]
    REAL(wp) :: ztend_qv(nproma,nlev)        !< moisture tendency from dynamics and physics before convection

    IF (ltimer) CALL timer_start(timer_convection)

    ! associate pointers
    field => prm_field(jg)
    tend  => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! Update physics state for input to convection
          zta  (jcs:jce,:)   =     field% ta  (jcs:jce,:,jb)   + pdtime*tend%   ta_phy(jcs:jce,:,jb)
          zqtrc(jcs:jce,:,:) = MAX(field% qtrc(jcs:jce,:,jb,:) + pdtime*tend% qtrc_phy(jcs:jce,:,jb,:), 0.0_wp)
          zua  (jcs:jce,:)   =     field% ua  (jcs:jce,:,jb)   + pdtime*tend%   ua_phy(jcs:jce,:,jb)
          zva  (jcs:jce,:)   =     field% va  (jcs:jce,:,jb)   + pdtime*tend%   va_phy(jcs:jce,:,jb)
          !
          zqtrc_cnd(jcs:jce,:) = zqtrc(jcs:jce,:,iqc) + zqtrc(jcs:jce,:,iqi)
          ztend_qv (jcs:jce,:) = tend%qtrc_dyn(jcs:jce,:,jb,iqv) + tend%qtrc_phy(jcs:jce,:,jb,iqv)
          !
          ! number of tracers excluding water vapour and hydrometeors
          ntrac = ntracer-iqt+1
          !
          CALL cumastr(jg,                           &! in
               &       jce, nproma,                  &! in
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
               &       field% lfland   (:,  jb),     &! in
               &             zqtrc     (:,:,   iqt:),&! in
               &       field% omega    (:,:,jb),     &! in
               &       field% evap     (:,  jb),     &! in
               &       field% presm_new(:,:,jb),     &! in
               &       field% presi_new(:,:,jb),     &! in
               &       field% geom     (:,:,jb),     &! in
               &       field% geoi     (:,:,jb),     &! in
               &             ztend_qv  (:,:),        &! in
               &       field% thvsig   (:,  jb),     &! in
               &       itype           (:),          &! out
               &       field% ictop    (:,  jb),     &! out
               &       field% rsfc     (:,  jb),     &! out
               &       field% ssfc     (:,  jb),     &! out
               &       field% con_dtrl (:,jb),       &! out
               &       field% con_dtri (:,jb),       &! out
               &       field% con_iteqv(:,jb),       &! out
               &       field% q_cnv    (:,:,jb),     &! out
               &        tend%   ua_cnv (:,:,jb),     &! out
               &        tend%   va_cnv (:,:,jb),     &! out
               &        tend% qtrc_cnv (:,:,jb,iqv), &! out
               &        tend% qtrc_cnv (:,:,jb,iqt:),&! out
               &        tend% qtrc_cnv (:,:,jb,iqc), &! out
               &        tend% qtrc_cnv (:,:,jb,iqi), &! out
               &             ztop      (:)           )! out
          !
          ! store convection type as real value
          field% rtype(jcs:jce,jb) = REAL(itype(jcs:jce),wp)
          !
          ! keep minimum conv. cloud top pressure (= max. conv. cloud top height) of this output interval
          field% topmax(jcs:jce,jb) = MIN(field% topmax(jcs:jce,jb),ztop(jcs:jce))
          !
          ! vertical integral
          field% q_cnv_vi(jcs:jce,jb) = SUM(field% q_cnv(jcs:jce,:,jb),DIM=2)
          !
       END IF
       !
       ! convert    heating
       tend% ta_cnv(jcs:jce,:,jb) = field% q_cnv(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field% q_cnv(jcs:jce,:,jb)
       !
       ! accumulate tendencies
       tend%   ua_phy(jcs:jce,:,jb)      = tend%   ua_phy(jcs:jce,:,jb)      + tend%   ua_cnv(jcs:jce,:,jb)
       tend%   va_phy(jcs:jce,:,jb)      = tend%   va_phy(jcs:jce,:,jb)      + tend%   va_cnv(jcs:jce,:,jb)
       tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)      + tend%   ta_cnv(jcs:jce,:,jb)
       !
       tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv)
       tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc)  + tend% qtrc_cnv(jcs:jce,:,jb,iqc)
       tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi)  + tend% qtrc_cnv(jcs:jce,:,jb,iqi)
       tend% qtrc_phy(jcs:jce,:,jb,iqt:) = tend% qtrc_phy(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)
       !
    ELSE
       !
       field% rtype    (jcs:jce,jb) = 0.0_wp
       field% ictop    (jcs:jce,jb) = nlevm1
       field% rsfc     (jcs:jce,jb) = 0.0_wp
       field% ssfc     (jcs:jce,jb) = 0.0_wp
       field% con_dtrl (jcs:jce,jb) = 0.0_wp
       field% con_dtri (jcs:jce,jb) = 0.0_wp
       field% con_iteqv(jcs:jce,jb) = 0.0_wp
       !
       tend% ta_cnv(jcs:jce,:,jb) = 0.0_wp
       tend% ua_cnv(jcs:jce,:,jb) = 0.0_wp
       tend% va_cnv(jcs:jce,:,jb) = 0.0_wp
       !
       tend% qtrc_cnv(jcs:jce,:,jb,iqv ) = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqc ) = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqi ) = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqt:) = 0.0_wp
       !
    END IF

    IF (ltimer) CALL timer_stop(timer_convection)

  END SUBROUTINE echam_cnv
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_cnv
