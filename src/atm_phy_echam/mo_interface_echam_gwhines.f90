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

MODULE mo_interface_echam_gwhines

  USE mo_kind                ,ONLY: wp
  USE mo_math_constants      ,ONLY: pi

  USE mo_model_domain        ,ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, t_echam_phy_tend

  USE mo_gw_hines            ,ONLY: gw_hines

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_gw_hines

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_gwhines

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_gwhines(is_in_sd_ed_interval,    &
       &                             is_active,               &
       &                             patch, rl_start, rl_end, &
       &                             field, tend              )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    TYPE(t_patch)  , TARGET ,INTENT(in)    :: patch
    INTEGER                 ,INTENT(in)    :: rl_start, rl_end
    TYPE(t_echam_phy_field) ,POINTER       :: field    
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend

    INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (ltimer) call timer_start(timer_gw_hines)    
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       !
       CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       !
       CALL echam_gw_hines(is_in_sd_ed_interval,          &
            &              is_active,                     &
            &              jg, jb,jcs,jce, nproma,        &
            &              field, tend                    )
    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
    IF (ltimer) call timer_stop(timer_gw_hines)

  END SUBROUTINE interface_echam_gwhines
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE echam_gw_hines(is_in_sd_ed_interval,  &
       &                    is_active,             &
       &                    jg, jb,jcs,jce, nbdim, &
       &                    field, tend            )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    INTEGER                 ,INTENT(in)    :: jg
    INTEGER                 ,INTENT(in)    :: jb                  !< block index
    INTEGER                 ,INTENT(in)    :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in)    :: nbdim               !< size of this block  
    TYPE(t_echam_phy_field) ,POINTER       :: field
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend

    ! local
    REAL(wp) :: zdis_gwd(nbdim,nlev)  !< out, energy dissipation rate [J/s/kg]
    INTEGER  :: nc                    !< number of cells/columns from jcs to jce
    REAL(wp) :: zlat_deg(nbdim)       !< latitude in deg N

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! number of cells/columns from index jcs to jce
          nc = jce-jcs+1
          !
          ! latitude [degN]
          zlat_deg(jcs:jce) = field% clat(jcs:jce,jb) * 180._wp/pi
          !
          CALL gw_hines (jg                       ,&
               &         nbdim                    ,&
               &         jcs                      ,&
               &         jce                      ,&
               &         nc                       ,&
               &         nlev                     ,&
               &         field% presi_old(:,:,jb) ,&
               &         field% presm_old(:,:,jb) ,&
               &         field%   zh(:,:,jb)      ,&
               &         field%  rho(:,:,jb)      ,&
               &         field% mair(:,:,jb)      ,&
               &         field%   ta(:,:,jb)      ,&
               &         field%   ua(:,:,jb)      ,&
               &         field%   va(:,:,jb)      ,&
               &         zlat_deg(:)              ,&
!!$            &         aprflux(:,krow)          ,&
               &         zdis_gwd(:,:)            ,&
               &         tend%   ua_gwd(:,:,jb)   ,&
               &         tend%   va_gwd(:,:,jb)    )
          !
          ! heating
          field% q_gwd(jcs:jce,:,jb) = zdis_gwd(jcs:jce,:) * field%mair(jcs:jce,:,jb)
          !
          ! vertical integral
          field% q_gwd_vi(jcs:jce,jb) = SUM(field% q_gwd(jcs:jce,:,jb),DIM=2)
          !
       END IF
       !
       ! convert    heating
       tend% ta_gwd(jcs:jce,:,jb) = field% q_gwd(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field% q_gwd(jcs:jce,:,jb)
       !
       ! tendencies accumulated
       tend% ta_phy(jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb) + tend% ta_gwd(jcs:jce,:,jb)
       tend% ua_phy(jcs:jce,:,jb) = tend% ua_phy(jcs:jce,:,jb) + tend% ua_gwd(jcs:jce,:,jb)
       tend% va_phy(jcs:jce,:,jb) = tend% va_phy(jcs:jce,:,jb) + tend% va_gwd(jcs:jce,:,jb)
       !
    ELSE
       !
       tend% ta_gwd(jcs:jce,:,jb) = 0.0_wp
       tend% ua_gwd(jcs:jce,:,jb) = 0.0_wp
       tend% va_gwd(jcs:jce,:,jb) = 0.0_wp
       !
    END IF

  END SUBROUTINE echam_gw_hines
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_gwhines
