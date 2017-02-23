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

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_run_config,          ONLY: nlev, nlevm1, nlevp1,    &
    &                               iqv, iqc, iqi, iqt, ntracer
  USE mo_echam_phy_config,    ONLY: echam_phy_config
  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_ham_aerosol_params,  ONLY: ncdnc, nicnc
  USE mo_cumastr,             ONLY: cumastr
  USE mo_cloud,               ONLY: cloud

  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_convection
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_convection

  INTEGER  :: ntrac   !< # of tracers excluding water vapour and hydrometeors
                      !< (handled by sub-models, e.g., chemical species)
  REAL(wp) :: pdtime  !< time step

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_convection(patch, rl_start, rl_end, field, tend, &
    & ictop, zconv, zq_phy, ntracer_wout_vapour, in_pdtime)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER :: field    
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    INTEGER  :: ictop (:,:)             !< from massflux
    REAL(wp) :: zconv  (:,:,:)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zq_phy (:,:,:)       !< heating by whole ECHAM physics    [W/m2]
    INTEGER         ,INTENT(IN)  :: ntracer_wout_vapour   !< # of tracers excluding water vapour and hydrometeors
    REAL(wp)        ,INTENT(IN)  :: in_pdtime         !< time step

!     INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    
    
    pdtime = in_pdtime           ! the module timestep length
    ntrac = ntracer_wout_vapour  !# of tracers excluding water vapour and hydrometeors

!     jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (ltimer) CALL timer_start(timer_convection)
    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      CALL echam_convection(jb,jcs,jce, nproma, field, tend, ictop(:,jb), zconv(:,:,jb), zq_phy(:,:,jb))

    ENDDO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_convection)
 
  END SUBROUTINE interface_echam_convection
  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_convection(jb,jcs,jce, nbdim, field, tend, ictop, zconv, zq_phy)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block 
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    INTEGER  :: ictop (nbdim)             !< from massflux
    REAL(wp) :: zconv  (nbdim,nlev)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]

    ! local 
    INTEGER  :: itype(nbdim)              !< type of convection

    REAL(wp) :: ztop(nbdim)           !<  convective cloud top pressure [Pa]
    REAL(wp) :: zta    (nbdim,nlev)         !< provisional temperature         [K]
    REAL(wp) :: zqtrc  (nbdim,nlev,ntracer) !< provisional mass mixing ratios  [kg/kg]
    REAL(wp) :: zua    (nbdim,nlev)         !< provisional zonal      wind     [m/s]
    REAL(wp) :: zva    (nbdim,nlev)         !< provisional meridional wind     [m/s]
    REAL(wp) :: zqtrc_cnd(nbdim,nlev) !<  cloud condensate mixing ratio [kg/kg]
    REAL(wp) :: ztend_qv(nbdim,nlev)  !<  moisture tendency from dynamics and physics before convection
    REAL(wp) :: zq_cnv (nbdim,nlev)       !< heating by convection             [W/m2]

    !-------------------------------------------------------------------
    ! Update physics state for input to next parameterization
    zta  (:,:)   =     field% ta  (:,:,jb)   + pdtime*tend% ta  (:,:,jb)
    zqtrc(:,:,:) = MAX(field% qtrc(:,:,jb,:) + pdtime*tend% qtrc(:,:,jb,:), 0.0_wp)
    zua  (:,:)   =     field% ua  (:,:,jb)   + pdtime*tend% ua  (:,:,jb)
    zva  (:,:)   =     field% va  (:,:,jb)   + pdtime*tend% va  (:,:,jb)

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------

    zqtrc_cnd(:,:) = zqtrc(:,:,iqc) + zqtrc(:,:,iqi)
    ztend_qv(:,:)  = tend%qtrc_dyn(:,:,jb,iqv) + tend%qtrc_phy(:,:,jb,iqv)

    IF (echam_phy_config%lconv) THEN

      CALL cumastr(jce, nbdim,                   &! in
        &          nlev, nlevp1, nlevm1,         &! in
        &          pdtime,                       &! in
        &          field% zf       (:,:,jb),     &! in
        &          field% zh       (:,:,jb),     &! in
        &          field% mdry     (:,:,jb),     &! in
        &                zta       (:,:),        &! in
        &                zqtrc     (:,:,   iqv), &! in
        &                zqtrc_cnd (:,:),        &! in
        &                zua       (:,:),        &! in
        &                zva       (:,:),        &! in
        &          ntrac,                        &! in
        &          field% lfland   (:,  jb),     &! in
        &                zqtrc     (:,:,   iqt:),&! in
        &          field% omega    (:,:,jb),     &! in
        &          field% evap     (:,  jb),     &! in
        &          field% presm_new(:,:,jb),     &! in
        &          field% presi_new(:,:,jb),     &! in
        &          field% geom     (:,:,jb),     &! in
        &          field% geoi     (:,:,jb),     &! in
        &                ztend_qv  (:,:),        &! in
        &          field% thvsig   (:,  jb),     &! in
        &          itype,                        &! out
        &          ictop,                        &! out
        &          field% rsfc     (:,  jb),     &! out
        &          field% ssfc     (:,  jb),     &! out
        &          field% con_dtrl (:,jb),       &! out
        &          field% con_dtri (:,jb),       &! out
        &          field% con_iteqv(:,jb),       &! out
        &                   zq_cnv (:,:),        &! out
        &           tend%   ua_cnv (:,:,jb),     &! out
        &           tend%   va_cnv (:,:,jb),     &! out
        &           tend% qtrc_cnv (:,:,jb,iqv), &! out
        &           tend% qtrc_cnv (:,:,jb,iqt:),&! out
        &           tend% qtrc_cnv (:,:,jb,iqc), &! out
        &           tend% qtrc_cnv (:,:,jb,iqi), &! out
        &                ztop      (:)           )! out

      ! store convection type as real value
      field% rtype(:,jb) = REAL(itype(:),wp)

      ! keep minimum conv. cloud top pressure (= max. conv. cloud top height) of this output interval
      field% topmax(:,jb) = MIN(field% topmax(:,jb),ztop(:))

      ! heating accumulated
      zq_phy(:,:) = zq_phy(:,:) + zq_cnv(:,:)

      ! tendency
      tend% ta_cnv(:,:,jb) = zq_cnv(:,:)*zconv(:,:)

      ! tendencies accumulated
      tend%   ua(:,:,jb)      = tend%   ua(:,:,jb)      + tend%   ua_cnv(:,:,jb)
      tend%   va(:,:,jb)      = tend%   va(:,:,jb)      + tend%   va_cnv(:,:,jb)
      tend%   ta(:,:,jb)      = tend%   ta(:,:,jb)      + tend%   ta_cnv(:,:,jb)
      tend% qtrc(:,:,jb,iqv)  = tend% qtrc(:,:,jb,iqv)  + tend% qtrc_cnv(:,:,jb,iqv)
      tend% qtrc(:,:,jb,iqc)  = tend% qtrc(:,:,jb,iqc)  + tend% qtrc_cnv(:,:,jb,iqc)
      tend% qtrc(:,:,jb,iqi)  = tend% qtrc(:,:,jb,iqi)  + tend% qtrc_cnv(:,:,jb,iqi)
      tend% qtrc(:,:,jb,iqt:) = tend% qtrc(:,:,jb,iqt:) + tend% qtrc_cnv(:,:,jb,iqt:)

    ELSE ! NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

      ictop(:)   = nlev-1
      field% rtype(:,jb) = 0.0_wp

    ENDIF !lconv

    !-------------------------------------------------------------
    ! Update provisional physics state
    !
!!$    field% ta  (:,:,jb)     = field% ta  (:,:,jb)     + tend% ta  (:,:,jb)    *pdtime
!!$    field% qtrc(:,:,jb,iqv) = field% qtrc(:,:,jb,iqv) + tend% qtrc(:,:,jb,iqv)*pdtime
    field% qtrc(:,:,jb,iqc) = field% qtrc(:,:,jb,iqc) + tend% qtrc(:,:,jb,iqc)*pdtime
    field% qtrc(:,:,jb,iqi) = field% qtrc(:,:,jb,iqi) + tend% qtrc(:,:,jb,iqi)*pdtime
    !
    !-------------------------------------------------------------

  END SUBROUTINE echam_convection
  !---------------------------------------------------------------------


END MODULE mo_interface_echam_convection
