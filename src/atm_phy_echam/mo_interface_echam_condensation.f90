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
MODULE mo_interface_echam_condensation

  USE mo_kind,                ONLY: wp
  USE mo_run_config,          ONLY: nlev, iqv, iqc, iqi
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_cloud,               ONLY: cloud

  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_cloud
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_condensation

  INTEGER  :: jks     !< start index for vertical loops
                      !< (handled by sub-models, e.g., chemical species)
  REAL(wp) :: pdtime  !< time step

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_condensation(patch, rl_start, rl_end, start_level, field, tend, &
    & zcpair, zconv, zq_phy, ictop, in_pdtime)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end
    INTEGER         ,INTENT(IN)  :: start_level            !< start vertical level
    TYPE(t_echam_phy_field),   POINTER :: field    
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp) :: zcpair (:,:,:)       !< specific heat of moist air at const. pressure [J/K/kg]
    REAL(wp) :: zconv  (:,:,:)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zq_phy (:,:,:)       !< heating by whole ECHAM physics    [W/m2]
    INTEGER  :: ictop (:,:)             !< from massflux
    REAL(wp)        ,INTENT(IN)  :: in_pdtime         !< time step

!     INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    
    
    pdtime = in_pdtime           ! the module timestep length
    jks   = start_level          ! start level


!     jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (ltimer) CALL timer_start(timer_cloud)
    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      CALL echam_condensation(jb,jcs,jce, nproma, field, tend, zcpair(:,:,jb), zconv(:,:,jb), zq_phy(:,:,jb), ictop(:,jb))

    ENDDO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_cloud)
 
  END SUBROUTINE interface_echam_condensation
  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_condensation(jb,jcs,jce, nbdim, field, tend, zcpair, zconv, zq_phy, ictop)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block 
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(IN) :: zcpair  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    REAL(wp) :: zconv  (nbdim,nlev)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]
    INTEGER  :: ictop (nbdim)             !< from massflux

    ! local 
    INTEGER  :: itype(nbdim)              !< type of convection

    REAL(wp) :: zq_cld (nbdim,nlev)       !< heating by stratiform clouds      [W/m2]


    !-------------------------------------------------------------
    ! 7. LARGE SCALE CONDENSATION.
    !-------------------------------------------------------------
    !IF (lcotra) CALL get_col_pol( tend%ta(:,:,jb),tend%qtrc(:,:,jb,iqv),jb )
    itype(:) = NINT(field%rtype(:,jb))

    CALL cloud(jce, nbdim, jks, nlev,        &! in
      &        pdtime,                       &! in
      &        ictop,                        &! in (from "cucall")
      &        field% presm_old(:,:,jb),     &! in
      &        field% dz       (:,:,jb),     &! in
      &        field% mdry     (:,:,jb),     &! in
      &        field% rho      (:,:,jb),     &! in
      &               zcpair   (:,:),        &! in
      &        field% acdnc    (:,:,jb),     &! in  acdnc
      &        field% ta       (:,:,jb),     &! in  tm1
      &        field% qtrc     (:,:,jb,iqv), &! in  qm1
      &        field% qtrc     (:,:,jb,iqc), &! in  xlm1
      &        field% qtrc     (:,:,jb,iqi), &! in  xim1
      &         tend% ta       (:,:,jb),     &! in  tte
      &         tend% qtrc     (:,:,jb,iqv), &! in  qte
      !
      &        itype,                        &! inout
      &        field% aclc     (:,:,jb),     &! inout
      !
      &        field% aclcov   (:,  jb),     &! out
      &        field% rsfl     (:,  jb),     &! out
      &        field% ssfl     (:,  jb),     &! out
      &        field% relhum   (:,:,jb),     &! out
      &               zq_cld   (:,:),        &! out
      &         tend% qtrc_cld (:,:,jb,iqv), &! out
      &         tend% qtrc_cld (:,:,jb,iqc), &! out
      &         tend% qtrc_cld (:,:,jb,iqi)  )! out


    field% rtype(:,jb) = REAL(itype(:),wp)

    ! heating accumulated
    zq_phy(:,:) = zq_phy(:,:) + zq_cld(:,:)

    ! tendency
    tend% ta_cld(:,:,jb) = zq_cld(:,:)*zconv(:,:)

    ! tendencies accumulated
    tend%   ta(:,:,jb)      = tend%   ta(:,:,jb)      + tend%   ta_cld(:,:,jb)
    tend% qtrc(:,:,jb,iqv)  = tend% qtrc(:,:,jb,iqv)  + tend% qtrc_cld(:,:,jb,iqv)
    tend% qtrc(:,:,jb,iqc)  = tend% qtrc(:,:,jb,iqc)  + tend% qtrc_cld(:,:,jb,iqc)
    tend% qtrc(:,:,jb,iqi)  = tend% qtrc(:,:,jb,iqi)  + tend% qtrc_cld(:,:,jb,iqi)

  END SUBROUTINE echam_condensation
  !---------------------------------------------------------------------


END MODULE mo_interface_echam_condensation
