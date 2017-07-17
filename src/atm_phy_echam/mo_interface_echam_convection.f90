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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_interface_echam_convection

  USE mo_kind                ,ONLY: wp

  USE mo_model_domain        ,ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev, nlevm1, nlevp1, iqv, iqc, iqi, iqt, ntracer

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, t_echam_phy_tend

  USE mo_cumastr             ,ONLY: cumastr

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_convection

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_convection

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_convection(is_in_sd_ed_interval,    &
       &                                is_active,               &
       &                                patch, rl_start, rl_end, &
       &                                field, tend,             &
       &                                ntrac, ictop,            &
       &                                pdtime                   )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch
    INTEGER                 ,INTENT(in)    :: rl_start, rl_end
    TYPE(t_echam_phy_field) ,POINTER       :: field    
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    INTEGER                 ,INTENT(in)    :: ntrac                              !< # of tracers excl. water vapour and hydrometeors
    INTEGER                 ,INTENT(inout) :: ictop (:,:)                        !< from massflux
    REAL(wp)                ,INTENT(in)    :: pdtime

    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (ltimer) CALL timer_start(timer_convection)
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       !
       CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       !
       CALL echam_convection(is_in_sd_ed_interval,           &
            &                is_active,                      &
            &                jb, jcs,jce, nproma, ntrac,     &
            &                field, tend,                    &
            &                ictop(:,jb),                    &
            &                pdtime                          )
       !
    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_convection)

  END SUBROUTINE interface_echam_convection
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE echam_convection(is_in_sd_ed_interval,         &
       &                      is_active,                    &
       &                      jb,jcs,jce, nbdim, ntrac,     &
       &                      field, tend,                  &
       &                      ictop,                        &
       &                      pdtime                        )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    INTEGER                 ,INTENT(in)    :: jb                  !< block index
    INTEGER                 ,INTENT(in)    :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in)    :: nbdim               !< size of this block 
    INTEGER                 ,INTENT(in)    :: ntrac               !< non-water tracers
    TYPE(t_echam_phy_field) ,POINTER       :: field
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    INTEGER                 ,INTENT(inout) :: ictop  (nbdim)      !< from massflux
    REAL(wp)                ,INTENT(in)    :: pdtime

    ! local 
    INTEGER  :: itype(nbdim)                !< type of convection

    REAL(wp) :: ztop(nbdim)                 !< convective cloud top pressure   [Pa]
    REAL(wp) :: zta    (nbdim,nlev)         !< provisional temperature         [K]
    REAL(wp) :: zqtrc  (nbdim,nlev,ntracer) !< provisional mass mixing ratios  [kg/kg]
    REAL(wp) :: zua    (nbdim,nlev)         !< provisional zonal      wind     [m/s]
    REAL(wp) :: zva    (nbdim,nlev)         !< provisional meridional wind     [m/s]
    REAL(wp) :: zqtrc_cnd(nbdim,nlev)       !< cloud condensate mixing ratio   [kg/kg]
    REAL(wp) :: ztend_qv(nbdim,nlev)        !< moisture tendency from dynamics and physics before convection

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! Update physics state for input to convection
          zta  (jcs:jce,:)   =     field% ta  (jcs:jce,:,jb)   + pdtime*tend% ta  (jcs:jce,:,jb)
          zqtrc(jcs:jce,:,:) = MAX(field% qtrc(jcs:jce,:,jb,:) + pdtime*tend% qtrc(jcs:jce,:,jb,:), 0.0_wp)
          zua  (jcs:jce,:)   =     field% ua  (jcs:jce,:,jb)   + pdtime*tend% ua  (jcs:jce,:,jb)
          zva  (jcs:jce,:)   =     field% va  (jcs:jce,:,jb)   + pdtime*tend% va  (jcs:jce,:,jb)
          !
          zqtrc_cnd(jcs:jce,:) = zqtrc(jcs:jce,:,iqc) + zqtrc(jcs:jce,:,iqi)
          ztend_qv (jcs:jce,:) = tend%qtrc_dyn(jcs:jce,:,jb,iqv) + tend%qtrc_phy(jcs:jce,:,jb,iqv)
          !
          CALL cumastr(jce, nbdim,                      &! in
               &       nlev, nlevp1, nlevm1,         &! in
               &       pdtime,                       &! in
               &       field% zf       (:,:,jb),     &! in
               &       field% zh       (:,:,jb),     &! in
               &       field% mdry     (:,:,jb),     &! in
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
               &       ictop           (:),          &! out
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
       tend% ua(jcs:jce,:,jb) = tend% ua(jcs:jce,:,jb) + tend% ua_cnv(jcs:jce,:,jb)
       tend% va(jcs:jce,:,jb) = tend% va(jcs:jce,:,jb) + tend% va_cnv(jcs:jce,:,jb)
       tend% ta(jcs:jce,:,jb) = tend% ta(jcs:jce,:,jb) + tend% ta_cnv(jcs:jce,:,jb)
       !
       tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv)
       tend% qtrc(jcs:jce,:,jb,iqc)  = tend% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_cnv(jcs:jce,:,jb,iqc)
       tend% qtrc(jcs:jce,:,jb,iqi)  = tend% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_cnv(jcs:jce,:,jb,iqi)
       tend% qtrc(jcs:jce,:,jb,iqt:) = tend% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)
       !
    ELSE
       !
       field% rtype    (jcs:jce,jb) = 0.0_wp
       ictop           (jcs:jce)    = nlevm1
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

  END SUBROUTINE echam_convection
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_convection
