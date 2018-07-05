!>
!! @brief Subroutine interface_echam_car calls Cariolle's linearized ozone scheme.
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

MODULE mo_interface_echam_car

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_car

  USE mo_run_config          ,ONLY: io3
  USE mo_physical_constants  ,ONLY: amd, amo3
  USE mo_bcs_time_interpolation ,ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  USE mo_lcariolle_types     ,ONLY: t_avi, t_time_interpolation
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_car

CONTAINS

  SUBROUTINE interface_echam_car(jg,jb,jcs,jce        ,&
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
    INTEGER                 ,POINTER    :: fc_car
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp)                            :: do3dt(nproma,nlev)
    TYPE(t_time_interpolation)          :: time_interpolation
    TYPE(t_time_interpolation_weights)  :: current_time_interpolation_weights
    TYPE(t_avi)                         :: avi

    EXTERNAL :: lcariolle_do3dt, lcariolle_lat_intp_li, lcariolle_pres_intp_li

    IF (ltimer) call timer_start(timer_car)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_car    => echam_phy_config(jg)%fc_car
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          current_time_interpolation_weights = calculate_time_interpolation_weights(datetime_old)
          time_interpolation% imonth1 = current_time_interpolation_weights% month1_index
          time_interpolation% imonth2 = current_time_interpolation_weights% month2_index
          time_interpolation% weight1 = current_time_interpolation_weights% weight1
          time_interpolation% weight2 = current_time_interpolation_weights% weight2
          !
          ALLOCATE(avi%o3_vmr(nproma,nlev), avi%vmr2molm2(nproma,nlev), avi%cell_center_lat(nproma), avi%lday(nproma))
          !
          avi%ldown=.TRUE.
          avi%o3_vmr(jcs:jce,:)        =  field% qtrc(jcs:jce,:,jb,io3)*amd/amo3
          avi%tmprt                    => field% ta  (:,:,jb)
          avi%vmr2molm2(jcs:jce,:)     =  field% mdry(jcs:jce,:,jb) / amd * 1.e3_wp
          avi%pres                     => field% presm_old(:,:,jb)
          avi%cell_center_lat(jcs:jce) =  field% clat(jcs:jce,jb)
          avi%lday(jcs:jce)            =  field% cosmu0(jcs:jce,jb) > 1.e-3_wp
          !
          CALL lcariolle_do3dt(jcs,                   jce,                    &
               &               nproma,                nlev,                   &
               &               time_interpolation,                            &
               &               lcariolle_lat_intp_li, lcariolle_pres_intp_li, &
               &               avi,                   do3dt                   )
          !
          DEALLOCATE(avi%o3_vmr, avi%vmr2molm2, avi%cell_center_lat, avi%lday)
          !
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_car)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend% qtrc_phy(jcs:jce,:,jb,io3) = tend% qtrc_phy(jcs:jce,:,jb,io3) + do3dt(jcs:jce,:)*amo3/amd
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field% qtrc(jcs:jce,:,jb,io3)  = field% qtrc(jcs:jce,:,jb,io3)  +  do3dt(jcs:jce,:)*amo3/amd*pdtime
       END IF
       !
    END IF
       
    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_car)
    NULLIFY(field)
    NULLIFY(tend )
    
    IF (ltimer) call timer_stop(timer_car)

  END SUBROUTINE interface_echam_car

END MODULE mo_interface_echam_car
