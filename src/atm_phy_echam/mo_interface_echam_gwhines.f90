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

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi
  USE mo_run_config,          ONLY: nlev
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_gw_hines
  USE mo_gw_hines,            ONLY: gw_hines
  USE mo_ssortns,             ONLY: ssodrag

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_gwhines

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_gwhines(patch, rl_start, rl_end, field, tend, zconv, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER :: field    
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(IN)    :: zconv  (nproma,nlev,patch%nblks_c)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp)        ,INTENT(INOUT) :: zq_phy (nproma,nlev,patch%nblks_c)       !< heating by whole ECHAM physics    [W/m2]

    INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    ! 6.1   CALL SUBROUTINE GW_HINES
    IF (phy_config%lgw_hines) THEN

      IF (ltimer) call timer_start(timer_gw_hines)
!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        CALL echam_gw_hines(patch, jg, jb,jcs,jce, nproma, field, tend, zconv(:,:,jb), zq_phy(:,:,jb))
      ENDDO
!$OMP END PARALLEL DO 

      IF (ltimer) call timer_stop(timer_gw_hines)

    ENDIF

  END SUBROUTINE interface_echam_gwhines
   !-------------------------------------------------------------------

 
  !---------------------------------------------------------------------
  SUBROUTINE echam_gw_hines(patch, jg, jb,jcs,jce, nbdim, field, tend, zconv, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN) :: jg
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block  
    REAL(wp)        ,INTENT(IN) :: zconv  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(INOUT) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]

    ! Temporary array used by GW_HINES
    REAL(wp) :: zdis_gwh(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]

    INTEGER  :: nc    !< number of cells/columns from (jce-jcs+1)
    REAL(wp) :: zlat_deg(nbdim)           !< latitude in deg N
    REAL(wp) :: zq_gwh (nbdim,nlev)       !< heating by atm. gravity waves     [W/m2]
 
    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    zlat_deg(jcs:jce) = field% clat(jcs:jce,jb) * 180._wp/pi

    CALL gw_hines ( jg                       ,&
      &             nbdim                    ,&
      &             jcs                      ,&
      &             jce                      ,&
      &             nc                       ,&
      &             nlev                     ,&
      &             field% presi_old(:,:,jb) ,&
      &             field% presm_old(:,:,jb) ,&
      &             field%   zh(:,:,jb)      ,&
      &             field%  rho(:,:,jb)      ,&
      &             field% mair(:,:,jb)      ,&
      &             field%   ta(:,:,jb)      ,&
      &             field%   ua(:,:,jb)      ,&
      &             field%   va(:,:,jb)      ,&
      &             zlat_deg(:)              ,&
!!$        &             aprflux(:,krow)          ,&
      &             zdis_gwh(:,:)            ,&
      &             tend%   ua_gwh(:,:,jb)   ,&
      &             tend%   va_gwh(:,:,jb) )


    ! heating
    zq_gwh(jcs:jce,:) = zdis_gwh(jcs:jce,:) * field%mair(jcs:jce,:,jb)

    ! heating accumulated
    zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_gwh(jcs:jce,:)

    ! tendency
    tend% ta_gwh(jcs:jce,:,jb) = zq_gwh(jcs:jce,:)*zconv(jcs:jce,:)

    ! tendencies accumulated
    tend%   ta(jcs:jce,:,jb) = tend%   ta(jcs:jce,:,jb) + tend%   ta_gwh(jcs:jce,:,jb)
    tend%   ua(jcs:jce,:,jb) = tend%   ua(jcs:jce,:,jb) + tend%   ua_gwh(jcs:jce,:,jb)
    tend%   va(jcs:jce,:,jb) = tend%   va(jcs:jce,:,jb) + tend%   va_gwh(jcs:jce,:,jb)

  END SUBROUTINE echam_gw_hines
  !---------------------------------------------------------------------


END MODULE mo_interface_echam_gwhines
