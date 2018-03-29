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

MODULE mo_interface_echam_gwd

  USE mo_kind                ,ONLY: wp

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev

  USE mtime                  ,ONLY: datetime
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_gw_hines

  USE mo_gw_hines            ,ONLY: gw_hines
!!$  USE mo_math_constants      ,ONLY: pi

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_gwd

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE echam_gwd(is_in_sd_ed_interval, &
       &               is_active,            &
       &               jg, jb,jcs,jce,       &
       &               datetime_old,         &
       &               pdtime                )

    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    INTEGER                 ,INTENT(in) :: jg
    INTEGER                 ,INTENT(in) :: jb                  !< block index
    INTEGER                 ,INTENT(in) :: jcs, jce            !< start/end column index within this block
    TYPE(datetime)          ,POINTER    :: datetime_old        !< generic input, not used in echam_gwd
    REAL(wp)                ,INTENT(in) :: pdtime              !< generic input, not used in echam_gwd

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    !
    REAL(wp) :: zdis_gwd(nproma,nlev) !< out, energy dissipation rate [J/s/kg]
    INTEGER  :: nc                    !< number of cells/columns from jcs to jce
!!$    REAL(wp) :: zlat_deg(nproma)       !< latitude in deg N

    IF (ltimer) call timer_start(timer_gw_hines)    

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

    IF (ltimer) call timer_stop(timer_gw_hines)

  END SUBROUTINE echam_gwd
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_gwd
