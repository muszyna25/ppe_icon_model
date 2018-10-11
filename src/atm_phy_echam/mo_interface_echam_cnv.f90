!>
!! @brief Subroutine interface_echam_cnv calls the Tiedtke/Nordeng convection scheme.
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

MODULE mo_interface_echam_cnv

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cnv

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqt, ntracer
  USE mo_cumastr             ,ONLY: cumastr

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_cnv

CONTAINS

  SUBROUTINE interface_echam_cnv(jg,jb,jcs,jce        ,&
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
    INTEGER                 ,POINTER    :: fc_cnv
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER  :: itype(nproma)                !< type of convection
    INTEGER  :: nlevm1, nlevp1
    INTEGER  :: ntrac                        !< # of tracers excluding water vapour and hydrometeors
                                             !< which are convectively transported

    LOGICAL  :: ldland(nproma)               !< land sea mask for using dlev_land or dlev_ocean

    REAL(wp) :: ztop(nproma)                 !< convective cloud top pressure   [Pa]
    REAL(wp) :: zta    (nproma,nlev)         !< provisional temperature         [K]
    REAL(wp) :: zqtrc  (nproma,nlev,ntracer) !< provisional mass mixing ratios  [kg/kg]
    REAL(wp) :: zua    (nproma,nlev)         !< provisional zonal      wind     [m/s]
    REAL(wp) :: zva    (nproma,nlev)         !< provisional meridional wind     [m/s]
    !
    REAL(wp) :: zqtrc_cnd(nproma,nlev)       !< cloud condensate mixing ratio   [kg/kg]
    REAL(wp) :: ztend_qv(nproma,nlev)        !< moisture tendency from dynamics and physics before convection

    IF (ltimer) CALL timer_start(timer_cnv)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_cnv    => echam_phy_config(jg)%fc_cnv
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    nlevm1 = nlev-1
    nlevp1 = nlev+1
    
    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          zta  (jcs:jce,:)   =     field% ta  (jcs:jce,:,jb)
          zqtrc(jcs:jce,:,:) = MAX(field% qtrc(jcs:jce,:,jb,:), 0.0_wp)
          zua  (jcs:jce,:)   =     field% ua  (jcs:jce,:,jb)
          zva  (jcs:jce,:)   =     field% va  (jcs:jce,:,jb)
          !
          ! Prepare input
          zqtrc_cnd(jcs:jce,:)   = zqtrc(jcs:jce,:,iqc) + zqtrc(jcs:jce,:,iqi)
          ztend_qv (jcs:jce,:)   = tend%qtrc_dyn(jcs:jce,:,jb,iqv) + tend%qtrc_phy(jcs:jce,:,jb,iqv)
          !
          ! number of tracers excluding water vapour and hydrometeors
          ntrac = ntracer-iqt+1
          !
          ! land sea mask for using dlev_land or dlev_ocean
          ldland(jcs:jce)    = field% sftlf(jcs:jce,jb) > 0._wp
          !
          CALL cumastr(jg, jb,                       &! in
               &       jcs, jce, nproma,             &! in
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
               &             ldland    (:),          &! in
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
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_cnv)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend%   ua_phy(jcs:jce,:,jb)      = tend%   ua_phy(jcs:jce,:,jb)      + tend%   ua_cnv(jcs:jce,:,jb)
          tend%   va_phy(jcs:jce,:,jb)      = tend%   va_phy(jcs:jce,:,jb)      + tend%   va_cnv(jcs:jce,:,jb)
          tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)      + tend%   ta_cnv(jcs:jce,:,jb)
          tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv)
          tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc)  + tend% qtrc_cnv(jcs:jce,:,jb,iqc)
          tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi)  + tend% qtrc_cnv(jcs:jce,:,jb,iqi)
          tend% qtrc_phy(jcs:jce,:,jb,iqt:) = tend% qtrc_phy(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field%   ua(jcs:jce,:,jb)      = field%   ua(jcs:jce,:,jb)      + tend%   ua_cnv(jcs:jce,:,jb)     *pdtime
          field%   va(jcs:jce,:,jb)      = field%   va(jcs:jce,:,jb)      + tend%   va_cnv(jcs:jce,:,jb)     *pdtime
          field%   ta(jcs:jce,:,jb)      = field%   ta(jcs:jce,:,jb)      + tend%   ta_cnv(jcs:jce,:,jb)     *pdtime
          field% qtrc(jcs:jce,:,jb,iqv)  = field% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv) *pdtime
          field% qtrc(jcs:jce,:,jb,iqc)  = field% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_cnv(jcs:jce,:,jb,iqc) *pdtime
          field% qtrc(jcs:jce,:,jb,iqi)  = field% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_cnv(jcs:jce,:,jb,iqi) *pdtime
          field% qtrc(jcs:jce,:,jb,iqt:) = field% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)*pdtime
       END IF
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
       tend%   ta_cnv(jcs:jce,:,jb)      = 0.0_wp
       tend%   ua_cnv(jcs:jce,:,jb)      = 0.0_wp
       tend%   va_cnv(jcs:jce,:,jb)      = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqv ) = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqc ) = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqi ) = 0.0_wp
       tend% qtrc_cnv(jcs:jce,:,jb,iqt:) = 0.0_wp
       !
    END IF

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_cnv)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_cnv)

  END SUBROUTINE interface_echam_cnv

END MODULE mo_interface_echam_cnv
