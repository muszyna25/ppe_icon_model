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

MODULE mo_interface_echam_sso

  USE mo_kind                ,ONLY: wp

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev

  USE mtime                  ,ONLY: datetime
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_ssodrag

  USE mo_ssodrag             ,ONLY: ssodrag
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_sso

CONTAINS

  SUBROUTINE echam_sso(is_in_sd_ed_interval, &
       &               is_active,            &
       &               jg, jb,jcs,jce,       &
       &               datetime_old,         &
       &               pdtime                )

    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    INTEGER                 ,INTENT(in) :: jg                  !< grid  index
    INTEGER                 ,INTENT(in) :: jb                  !< block index
    INTEGER                 ,INTENT(in) :: jcs, jce            !< start/end column index within this block
    TYPE(datetime)          ,POINTER    :: datetime_old        !< generic input, not used in echam_sso
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    !
    REAL(wp) :: zdis_sso(nproma,nlev)  !<  out, energy dissipation rate [J/s/kg]
    INTEGER  :: nc

    IF (ltimer) call timer_start(timer_ssodrag)

    ! associate pointers
    field => prm_field(jg)
    tend  => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! number of cells/columns from index jcs to jce
          nc = jce-jcs+1
          !
          CALL ssodrag(jg                           ,& ! in,  grid index
               &       nc                           ,& ! in,  number of cells/columns in loop (jce-jcs+1)
               &       nproma                       ,& ! in,  dimension of block of cells/columns
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
       tend% ta_phy(jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb) + tend% ta_sso(jcs:jce,:,jb)
       tend% ua_phy(jcs:jce,:,jb) = tend% ua_phy(jcs:jce,:,jb) + tend% ua_sso(jcs:jce,:,jb)
       tend% va_phy(jcs:jce,:,jb) = tend% va_phy(jcs:jce,:,jb) + tend% va_sso(jcs:jce,:,jb)

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
    
    IF (ltimer) call timer_stop(timer_ssodrag)

  END SUBROUTINE echam_sso

END MODULE mo_interface_echam_sso
