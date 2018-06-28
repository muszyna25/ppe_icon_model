!>
!! @brief Subroutine interface_echam_gwd calls the Hines gravity wave scheme.
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

MODULE mo_interface_echam_gwd

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
       &                            t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_gwd

  USE mo_gw_hines            ,ONLY: gw_hines
!!$  USE mo_math_constants      ,ONLY: pi

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_gwd

CONTAINS

  SUBROUTINE interface_echam_gwd(jg, jb,jcs,jce       ,&
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
    INTEGER                 ,POINTER    :: fc_gwd
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp) :: zdis_gwd(nproma,nlev) !< out, energy dissipation rate [J/s/kg]
    INTEGER  :: nc                    !< number of cells/columns from jcs to jce
!!$    REAL(wp) :: zlat_deg(nproma)       !< latitude in deg N

    IF (ltimer) call timer_start(timer_gwd)    

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_gwd    => echam_phy_config(jg)%fc_gwd
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          ! number of cells/columns from index jcs to jce
          nc = jce-jcs+1
          !
!!$          ! latitude [degN]
!!$          zlat_deg(jcs:jce) = field% clat(jcs:jce,jb) * 180._wp/pi
          !
          CALL gw_hines (jg                       ,&
               &         nproma                   ,&
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
!!$               &         zlat_deg(:)              ,&
!!$               &         aprflux(:,krow)          ,&
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
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_gwd)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend% ta_phy(jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb) + tend% ta_gwd(jcs:jce,:,jb)
          tend% ua_phy(jcs:jce,:,jb) = tend% ua_phy(jcs:jce,:,jb) + tend% ua_gwd(jcs:jce,:,jb)
          tend% va_phy(jcs:jce,:,jb) = tend% va_phy(jcs:jce,:,jb) + tend% va_gwd(jcs:jce,:,jb)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field% ta(jcs:jce,:,jb) = field% ta(jcs:jce,:,jb) + tend% ta_gwd(jcs:jce,:,jb)*pdtime
          field% ua(jcs:jce,:,jb) = field% ua(jcs:jce,:,jb) + tend% ua_gwd(jcs:jce,:,jb)*pdtime
          field% va(jcs:jce,:,jb) = field% va(jcs:jce,:,jb) + tend% va_gwd(jcs:jce,:,jb)*pdtime
       END IF
       !
    ELSE
       !
       tend% ta_gwd(jcs:jce,:,jb) = 0.0_wp
       tend% ua_gwd(jcs:jce,:,jb) = 0.0_wp
       tend% va_gwd(jcs:jce,:,jb) = 0.0_wp
       !
    END IF

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_gwd)
    NULLIFY(field)
    NULLIFY(tend )

    IF (ltimer) call timer_stop(timer_gwd)

  END SUBROUTINE interface_echam_gwd

END MODULE mo_interface_echam_gwd