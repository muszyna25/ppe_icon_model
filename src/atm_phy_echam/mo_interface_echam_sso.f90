!>
!! @brief Subroutine interface_echam_sso calls the SSO parameterization scheme.
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
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_sso

  USE mo_echam_sso_config    ,ONLY: echam_sso_config
  USE mo_ssodrag             ,ONLY: ssodrag
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_sso

CONTAINS

  SUBROUTINE interface_echam_sso(jg, jb,jcs,jce       ,&
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
    INTEGER                 ,POINTER    :: fc_sso
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    LOGICAL                 ,POINTER    :: lsftlf   ! <-- provisional

    ! Local variables
    !
    REAL(wp) :: zdis_sso(nproma,nlev)  !<  out, energy dissipation rate [J/s/kg]
    INTEGER  :: nc
    REAL(wp) :: zscale(nproma)         !< area scaling factor

    IF (ltimer) call timer_start(timer_sso)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_sso    => echam_phy_config(jg)%fc_sso
    field     => prm_field(jg)
    tend      => prm_tend (jg)
    lsftlf    => echam_sso_config(jg)%lsftlf        ! <-- provisional

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! number of cells/columns from index jcs to jce
          nc = jce-jcs+1
          !
          ! area scaling factor                     ! <-- provisional
          IF (lsftlf) THEN                          ! <-- provisional
             zscale(:) = field% sftlf (:,jb)
          ELSE                                      ! <-- provisional
             zscale(:) = 1.0_wp                     ! <-- provisional
          END IF                                    ! <-- provisional
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
               &              zscale(:)             ,& ! in,  area fraction of land incl. lakes
               !                                              where the SSO params are valid
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
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_sso)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend% ta_phy(jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb) + tend% ta_sso(jcs:jce,:,jb)
          tend% ua_phy(jcs:jce,:,jb) = tend% ua_phy(jcs:jce,:,jb) + tend% ua_sso(jcs:jce,:,jb)
          tend% va_phy(jcs:jce,:,jb) = tend% va_phy(jcs:jce,:,jb) + tend% va_sso(jcs:jce,:,jb)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field% ta(jcs:jce,:,jb) = field% ta(jcs:jce,:,jb) + tend% ta_sso(jcs:jce,:,jb)*pdtime
          field% ua(jcs:jce,:,jb) = field% ua(jcs:jce,:,jb) + tend% ua_sso(jcs:jce,:,jb)*pdtime
          field% va(jcs:jce,:,jb) = field% va(jcs:jce,:,jb) + tend% va_sso(jcs:jce,:,jb)*pdtime
       END IF
       !
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
    
    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_sso)
    NULLIFY(field)
    NULLIFY(tend)
    NULLIFY(lsftlf)                                 ! <-- provisional

    IF (ltimer) call timer_stop(timer_sso)

  END SUBROUTINE interface_echam_sso

END MODULE mo_interface_echam_sso
