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
    &                               iqv, iqc, iqi, iqt
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_ham_aerosol_params,  ONLY: ncdnc, nicnc
  USE mo_cumastr,             ONLY: cucall
  USE mo_cloud,               ONLY: cloud

  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_cucall, timer_cloud
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_convection

  INTEGER  :: jks     !< start index for vertical loops
  INTEGER  :: ntrac   !< # of tracers excluding water vapour and hydrometeors
                      !< (handled by sub-models, e.g., chemical species)
  REAL(wp) :: pdtime  !< time step
  REAL(wp) :: zcd     !< specific heat of dry air          [J/K/kg]
  REAL(wp) :: zcv     !< specific heat of water vapor      [J/K/kg]

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_convection(patch, rl_start, rl_end, start_level, field, tend, &
    & zcair, ntracer_wout_vapour, in_zcd, in_zcv, in_pdtime)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end
    INTEGER         ,INTENT(IN)  :: start_level            !< start vertical level
    TYPE(t_echam_phy_field),   POINTER :: field    
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp) :: zcair(:,:,:) !  (nproma,nlev,patch%nblks_c)       !< specific heat of moist air        [J/K/kg]
    INTEGER         ,INTENT(IN)  :: ntracer_wout_vapour   !< # of tracers excluding water vapour and hydrometeors
    REAL(wp) :: in_zcd                       !< specific heat of dry air          [J/K/kg]
    REAL(wp) :: in_zcv                       !< specific heat of water vapor      [J/K/kg]
    REAL(wp)        ,INTENT(IN)  :: in_pdtime         !< time step

    INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    
    
    pdtime = in_pdtime           ! the module timestep length
    jks   = start_level          ! start level
    ntrac = ntracer_wout_vapour  !# of tracers excluding water vapour and hydrometeors
    zcd = in_zcd
    zcv = in_zcv


    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      CALL echam_cumulus_condensation(jb,jcs,jce, nproma, field, tend, zcair(:,:,jb))

    ENDDO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
 
  END SUBROUTINE interface_echam_convection
  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_cumulus_condensation(jb,jcs,jce, nbdim, field, tend, zcair)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block 
    REAL(wp)        ,INTENT(IN) :: zcair  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    ! local 
    REAL(wp) :: zqtec  (nbdim,nlev)       !< tracer tendency due to entrainment/detrainment
    INTEGER  :: itype(nbdim)              !< type of convection
    INTEGER  :: ictop (nbdim)             !< from massflux
    INTEGER  :: ilab   (nproma,nlev)

    itype(jcs:jce) = 0

    ! 7.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
    !       AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS

    tend% xl_dtr(jcs:jce,:,jb) = 0._wp
    tend% xi_dtr(jcs:jce,:,jb) = 0._wp
    zqtec  (jcs:jce,:) = 0._wp

    field% rsfc(:,jb) = 0._wp
    field% ssfc(:,jb) = 0._wp

    ! 7.2   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION

    IF (phy_config%lconv) THEN

      IF (ltimer) call timer_start(timer_cucall)

      CALL cucall( jce, nbdim, nlev,          &! in
        &          nlevp1, nlevm1,            &! in
        &          ntrac,                     &! in     tracers
!        &          jb,                        &! in     row index
        &          pdtime,                    &! in
        &          field% lfland(:,jb),       &! in     loland
        &          field% ta(:,:,jb),         &! in     tm1
        &          field% ua(:,:,jb),         &! in     um1
        &          field% va(:,:,jb),         &! in     vm1
        &          field% qtrc(:,:,jb,iqv),   &! in     qm1
        &          field% qtrc(:,:,jb,iqc),   &! in     xlm1
        &          field% qtrc(:,:,jb,iqi),   &! in     xim1
        &          field% qtrc(:,:,jb,iqt:),  &! in     xtm1
        &          tend% qtrc(:,:,jb,iqv),    &! in     qte  for internal updating
        &          tend% qtrc(:,:,jb,iqc),    &! in     xlte
        &          tend% qtrc(:,:,jb,iqi),    &! in     xite
        &          field% omega(:,:,jb),      &! in     vervel
        &          field% evap(:,jb),         &! in     qhfla (from "vdiff")
        &          field% geom(:,:,jb),       &! in     geom1
        &          field% presm_new(:,:,jb),  &! in     app1
        &          field% presi_new(:,:,jb),  &! in     aphp1
        &          field% thvsig(:,jb),       &! in           (from "vdiff")
        &          tend% ta(:,:,jb),          &! in     tte  for internal updating
        &          tend% ua(:,:,jb),          &! in     vom  for internal updating
        &          tend% va(:,:,jb),          &! in     vol  for internal updating
        &          tend% qtrc(:,:,jb,iqt:),   &! in     xtte for internal updating
        &          zqtec,                     &! inout
        &          field% ch_concloud(:,jb),  &! inout condensational heat
        &          field% cw_concloud(:,jb),  &! inout condensational heat
        &          field% rsfc(:,jb),         &! out
        &          field% ssfc(:,jb),         &! out
        &          tend% xl_dtr(:,:,jb),      &! inout  xtecl
        &          tend% xi_dtr(:,:,jb),      &! inout  xteci
        &          itype,                     &! inout
        &          ictop,                     &! out
        &          ilab,                      &! out
        &          field% topmax(:,jb),       &! inout
        &          echam_conv_config%cevapcu, &! in
        &          zcd, zcv,                  &! in
        &          tend% qtrc_dyn(:,:,jb,iqv),&! in     qte by transport
        &          tend% qtrc_phy(:,:,jb,iqv),&! in     qte by physics
        &          field% con_dtrl(:,jb),     &! inout detrained liquid
        &          field% con_dtri(:,jb),     &! inout detrained ice
        &          field% con_iteqv(:,jb),    &! inout v. int. tend of water vapor within conv
        &          tend%  ta_cnv(:,:,jb),     &! out
        &          tend%  ua_cnv(:,:,jb),     &! out
        &          tend%  va_cnv(:,:,jb),     &! out
        &          tend%qtrc_cnv(:,:,jb,iqv), &! out
        &          tend%qtrc_cnv(:,:,jb,iqt:) )! out

      IF (ltimer) CALL timer_stop(timer_cucall)

      field% rtype(jcs:jce,jb) = REAL(itype(jcs:jce),wp)

!!$      ! heating accumulated
!!$      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_cnv(jcs:jce,:)
!!$
!!$      ! tendency
!!$      tend% temp_cnv(jcs:jce,:,jb) = zq_cnv(jcs:jce,:)*zconv(jcs:jce,:)
!!$
      ! tendencies accumulated
      tend%   ua(jcs:jce,:,jb)      = tend%   ua(jcs:jce,:,jb)      + tend%   ua_cnv(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb)      = tend%   va(jcs:jce,:,jb)      + tend%   va_cnv(jcs:jce,:,jb)
      tend%   ta(jcs:jce,:,jb)      = tend%   ta(jcs:jce,:,jb)      + tend%   ta_cnv(jcs:jce,:,jb)
      tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv)
      tend% qtrc(jcs:jce,:,jb,iqt:) = tend% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)


    ELSE ! NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

      ilab(jcs:jce,1:nlev) = 0
      ictop(jcs:jce)       = nlev-1

      tend%   ua_cnv(jcs:jce,:,jb)      = 0._wp
      tend%   va_cnv(jcs:jce,:,jb)      = 0._wp
      tend%   ta_cnv(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_cnv(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_cnv(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lconv

    !-------------------------------------------------------------
    ! 7. LARGE SCALE CONDENSATION.
    !-------------------------------------------------------------
    IF(phy_config%lcond) THEN

      !IF (lcotra) CALL get_col_pol( tend%ta(:,:,jb),tend%qtrc(:,:,jb,iqv),jb )

      IF (ncdnc==0 .AND. nicnc==0) THEN

        field% rsfl(:,jb) = 0._wp
        field% ssfl(:,jb) = 0._wp

        IF (ltimer) CALL timer_start(timer_cloud)

        CALL cloud(jce, nproma, jks, nlev, nlevp1, &! in
          &        pdtime,                    &! in
          &        ictop,                     &! in (from "cucall")
          &        field% presi_old(:,:,jb),  &! in
          &        field% presm_old(:,:,jb),  &! in
!          &        field% presm_new(:,:,jb), &! in
          &        field% acdnc (:,:,jb),     &! in. acdnc
          &        field%   ta  (:,:,jb),     &! in. tm1
          &        field%   tv  (:,:,jb),     &! in. ztvm1
          &        field% qtrc  (:,:,jb,iqv), &! in.  qm1
          &        field% qtrc  (:,:,jb,iqc), &! in. xlm1
          &        field% qtrc  (:,:,jb,iqi), &! in. xim1
          &        zcair(:,:),                &! in
          &        field% geom  (:,:,jb),     &! in. geom1
          &        field% aclcov(:,  jb),     &! out
          &        field%  qvi  (:,  jb),     &! out
          &        field% xlvi  (:,  jb),     &! out
          &        field% xivi  (:,  jb),     &! out
          &        itype,                     &!
          &        field% ch_concloud(:,jb),  &! inout condens. heat
          &        field% cw_concloud(:,jb),  &! inout condens. heat
          &         tend% xl_dtr(:,:,jb),     &! inout  xtecl
          &         tend% xi_dtr(:,:,jb),     &! inout  xteci
          &        zqtec,                     &! inout (there is a clip inside)
          &         tend% ta  (:,:,jb),     &! inout.  tte
          &         tend% qtrc  (:,:,jb,iqv), &! inout.  qte
          &         tend% qtrc  (:,:,jb,iqc), &! inout. xlte
          &         tend% qtrc  (:,:,jb,iqi), &! inout. xite
          &        field% cld_dtrl(:,jb),     &! inout detrained liquid
          &        field% cld_dtri(:,jb),     &! inout detrained ice
          &        field% cld_iteq(:,jb),     &! inout v. int. tend of qv,qc, and qi within cloud
!          &         tend% x_dtr(:,:,jb),      &! inout (there is a clip inside)
          &        field% aclc  (:,:,jb),     &! inout
          &        field% ssfl  (:,  jb),     &! out
          &        field% rsfl  (:,  jb),     &! out
          &        field% relhum(:,:,jb),     &! out
          &        tend%  ta_cld(:,:,jb),     &! out
          &        tend%qtrc_cld(:,:,jb,iqv), &! out
          &        tend%qtrc_cld(:,:,jb,iqc), &! out
          &        tend%qtrc_cld(:,:,jb,iqi)  )! out

        IF (ltimer) CALL timer_stop(timer_cloud)

      ELSE IF (ncdnc>0 .AND. nicnc>0) THEN
!0      CALL cloud_cdnc_icnc(...) !!skipped in ICON
      ELSE
        CALL finish('echam_phy_main', ' check setting of ncdnc and nicnc.')
      END IF

    ELSE ! NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.

      field% rsfl (jcs:jce,  jb) = 0._wp
      field% ssfl (jcs:jce,  jb) = 0._wp
      field% aclc (jcs:jce,:,jb) = 0._wp

      tend%   ta_cld(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqc)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqi)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lcond

  END SUBROUTINE echam_cumulus_condensation
  !---------------------------------------------------------------------


END MODULE mo_interface_echam_convection
