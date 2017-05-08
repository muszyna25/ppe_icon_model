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

MODULE mo_interface_echam_sso

  USE mo_kind                ,ONLY: wp

  USE mo_model_domain        ,ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, t_echam_phy_tend
  
  USE mo_ssodrag             ,ONLY: ssodrag
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_ssodrag

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_sso

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_sso(is_in_sd_ed_interval,    &
       &                         is_active,               &
       &                         patch, rl_start, rl_end, &
       &                         field, tend,             &
       &                         pdtime                   )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch
    INTEGER                 ,INTENT(in)    :: rl_start, rl_end
    TYPE(t_echam_phy_field) ,POINTER       :: field    
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    REAL(wp)                ,INTENT(in)    :: pdtime

    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jg             !< grid index
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    jg = patch%id
    
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)
 
    IF (ltimer) call timer_start(timer_ssodrag)
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       !
       CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       !
       CALL echam_ssodrag(is_in_sd_ed_interval,          &
            &             is_active,                     &
            &             jg, jb,jcs,jce, nproma,        &
            &             field, tend,                   &
            &             pdtime                         )
    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------

    IF (ltimer) call timer_stop(timer_ssodrag)

  END SUBROUTINE interface_echam_sso
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE echam_ssodrag(is_in_sd_ed_interval,  &
       &                   is_active,             &
       &                   jg, jb,jcs,jce, nbdim, &
       &                   field, tend,           &
       &                   pdtime                 )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in)    :: is_active
    INTEGER                 ,INTENT(in)    :: jg                  !< grid  index
    INTEGER                 ,INTENT(in)    :: jb                  !< block index
    INTEGER                 ,INTENT(in)    :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in)    :: nbdim               !< size of this block 
    TYPE(t_echam_phy_field) ,POINTER       :: field
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend
    REAL(wp)                ,INTENT(in)    :: pdtime

    ! local
    REAL(wp) :: zdis_sso(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]
    INTEGER  :: nc

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! number of cells/columns from index jcs to jce
          nc = jce-jcs+1
          !
          CALL ssodrag(jg                           ,& ! in,  grid index
               &       nc                           ,& ! in,  number of cells/columns in loop (jce-jcs+1)
               &       nbdim                        ,& ! in,  dimension of block of cells/columns
               &       nlev                         ,& ! in,  number of levels
               !
               &       pdtime                       ,& ! in,  time step length
               &       field% coriol(:,jb)          ,& ! in,  Coriolis parameter (1/s)
               &       field% zf  (:,:,jb)          ,& ! in,  full level height (m)
               &       field% zh  (:,nlev+1,jb)     ,& ! in,  surface height    (m)
               !
               &       field% presi_old(:,:,jb)     ,& ! in,  p at half levels
               &       field% presm_old(:,:,jb)     ,& ! in,  p at full levels
               &       field% mair(:,:,jb)          ,& ! in,  air mass
               &       field%   ta(:,:,jb)          ,& ! in,  T
               &       field%   ua(:,:,jb)          ,& ! in,  u
               &       field%   va(:,:,jb)          ,& ! in,  v
               !
               &       field% oromea(:,jb)          ,& ! in,  Mean Orography (m)
               &       field% orostd(:,jb)          ,& ! in,  SSO standard deviation (m)
               &       field% orosig(:,jb)          ,& ! in,  SSO slope
               &       field% orogam(:,jb)          ,& ! in,  SSO Anisotropy
               &       field% orothe(:,jb)          ,& ! in,  SSO Angle
               &       field% oropic(:,jb)          ,& ! in,  SSO Peaks elevation (m)
               &       field% oroval(:,jb)          ,& ! in,  SSO Valleys elevation (m)
               !
               &       field% u_stress_sso(:,jb)    ,& ! out, u-gravity wave stress
               &       field% v_stress_sso(:,jb)    ,& ! out, v-gravity wave stress
               &       field% dissipation_sso(:,jb) ,& ! out, dissipation by gravity wave drag
               !
               &       zdis_sso(:,:)                ,& ! out, energy dissipation rate
               &       tend%   ua_sso(:,:,jb)       ,& ! out, tendency of zonal wind
               &       tend%   va_sso(:,:,jb)        ) ! out, tendency of meridional wind
          !
          ! heating
          field% q_sso(jcs:jce,:,jb) = zdis_sso(jcs:jce,:) * field%mair(jcs:jce,:,jb)
          !
          ! vertical integral
          field% q_sso_vi(jcs:jce,jb) = SUM(field% q_sso(jcs:jce,:,jb),DIM=2)
          !
       END IF
       !
       ! convert    heating
       tend% ta_sso(jcs:jce,:,jb) = field% q_sso(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field% q_sso(jcs:jce,:,jb)

       ! accumulate tendencies
       tend% ta(jcs:jce,:,jb) = tend% ta(jcs:jce,:,jb) + tend% ta_sso(jcs:jce,:,jb)
       tend% ua(jcs:jce,:,jb) = tend% ua(jcs:jce,:,jb) + tend% ua_sso(jcs:jce,:,jb)
       tend% va(jcs:jce,:,jb) = tend% va(jcs:jce,:,jb) + tend% va_sso(jcs:jce,:,jb)

    ELSE
       !
       field% u_stress_sso   (jcs:jce,jb) = 0.0_wp
       field% v_stress_sso   (jcs:jce,jb) = 0.0_wp
       field% dissipation_sso(jcs:jce,jb) = 0.0_wp
       !
       tend% ta_sso(jcs:jce,:,jb) = 0.0_wp
       tend% ua_sso(jcs:jce,:,jb) = 0.0_wp
       tend% va_sso(jcs:jce,:,jb) = 0.0_wp
       !
    END IF
    
  END SUBROUTINE echam_ssodrag
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_sso
