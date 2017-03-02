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

  USE mo_kind,                ONLY: wp
  USE mo_run_config,          ONLY: nlev
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_ssodrag
  USE mo_ssortns,             ONLY: ssodrag

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_sso

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_sso(patch, rl_start, rl_end, field, tend, zconv, zq_phy, pdtime)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER :: field    
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(IN)    :: zconv  (nproma,nlev,patch%nblks_c)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp)        ,INTENT(INOUT) :: zq_phy (nproma,nlev,patch%nblks_c)       !< heating by whole ECHAM physics    [W/m2]
    REAL(wp)        ,INTENT(IN) :: pdtime         !< time step

    INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)
 
   ! 6.2   CALL SUBROUTINE SSODRAG
    IF (ltimer) call timer_start(timer_ssodrag)

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      CALL echam_ssodrag(jb,jcs,jce, nproma, field, tend, zconv(:,:,jb), zq_phy(:,:,jb), pdtime)
    ENDDO
!$OMP END PARALLEL DO 

    IF (ltimer) call timer_stop(timer_ssodrag)

  END SUBROUTINE interface_echam_sso
   !-------------------------------------------------------------------

   !---------------------------------------------------------------------
  SUBROUTINE echam_ssodrag(jb,jcs,jce, nbdim, field, tend, zconv, zq_phy, pdtime)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block 
    REAL(wp)        ,INTENT(IN) :: zconv  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(INOUT) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]
    REAL(wp)        ,INTENT(IN) :: pdtime         !< time step

    ! Temporary array used by SSODRAG
    REAL(wp) :: zdis_sso(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]
    REAL(wp) :: zq_sso (nbdim,nlev)       !< heating by subgrid scale orogr.   [W/m2]
    INTEGER :: nc

    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    CALL ssodrag( nc                                        ,& ! in,  number of cells/columns in loop (jce-jcs+1)
                  nbdim                                     ,& ! in,  dimension of block of cells/columns
                  nlev                                      ,& ! in,  number of levels
                  !
                  pdtime                                    ,& ! in,  time step length
                  field% coriol(:,jb)                       ,& ! in,  Coriolis parameter (1/s)
                  field% zf  (:,:,jb)                       ,& ! in,  full level height (m)
                  field% zh  (:,nlev+1,jb)                  ,& ! in,  surface height    (m)
                  !
                  field% presi_old(:,:,jb)                  ,& ! in,  p at half levels
                  field% presm_old(:,:,jb)                  ,& ! in,  p at full levels
                  field% mair(:,:,jb)                       ,& ! in,  air mass
                  field%   ta(:,:,jb)                       ,& ! in,  T
                  field%   ua(:,:,jb)                       ,& ! in,  u
                  field%   va(:,:,jb)                       ,& ! in,  v
                  !
                  field% oromea(:,jb)                       ,& ! in,  Mean Orography (m)
                  field% orostd(:,jb)                       ,& ! in,  SSO standard deviation (m)
                  field% orosig(:,jb)                       ,& ! in,  SSO slope
                  field% orogam(:,jb)                       ,& ! in,  SSO Anisotropy
                  field% orothe(:,jb)                       ,& ! in,  SSO Angle
                  field% oropic(:,jb)                       ,& ! in,  SSO Peaks elevation (m)
                  field% oroval(:,jb)                       ,& ! in,  SSO Valleys elevation (m)
                  !
                  field% u_stress_sso(:,jb)                 ,& ! out, u-gravity wave stress
                  field% v_stress_sso(:,jb)                 ,& ! out, v-gravity wave stress
                  field% dissipation_sso(:,jb)              ,& ! out, dissipation by gravity wave drag
                  !
                  zdis_sso(:,:)                             ,& ! out, energy dissipation rate
                  tend%   ua_sso(:,:,jb)                    ,& ! out, tendency of zonal wind
                  tend%   va_sso(:,:,jb)                     ) ! out, tendency of meridional wind


      ! heating
      zq_sso(jcs:jce,:) = zdis_sso(jcs:jce,:) * field%mair(jcs:jce,:,jb)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_sso(jcs:jce,:)

      ! tendency
      tend% ta_sso(jcs:jce,:,jb) = zq_sso(jcs:jce,:)*zconv(jcs:jce,:)

      ! tendencies accumulated
      tend%   ta(jcs:jce,:,jb) = tend%   ta(jcs:jce,:,jb) + tend%   ta_sso(jcs:jce,:,jb)
      tend%   ua(jcs:jce,:,jb) = tend%   ua(jcs:jce,:,jb) + tend%   ua_sso(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb) = tend%   va(jcs:jce,:,jb) + tend%   va_sso(jcs:jce,:,jb)

  END SUBROUTINE echam_ssodrag
  !---------------------------------------------------------------------

END MODULE mo_interface_echam_sso
