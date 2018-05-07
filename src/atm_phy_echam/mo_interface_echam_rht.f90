!>
!! @brief Subroutine interface_echam_rht calls the radiative heating scheme.
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

MODULE mo_interface_echam_rht

  USE mo_kind                   ,ONLY: wp
  USE mtime                     ,ONLY: datetime

  USE mo_echam_phy_config       ,ONLY: echam_phy_config
  USE mo_echam_phy_memory       ,ONLY: t_echam_phy_field, prm_field, &
    &                                  t_echam_phy_tend,  prm_tend

  USE mo_timer                  ,ONLY: ltimer, timer_start, timer_stop, timer_rht

  USE mo_radheating             ,ONLY: radheating
  USE mo_psrad_solar_parameters ,ONLY: psctm
  USE mo_ext_data_state         ,ONLY: ext_data


  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_rht

CONTAINS

  SUBROUTINE interface_echam_rht(jg, jb,jcs,jce       ,&
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
    INTEGER                 ,POINTER    :: fc_rht
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER                             :: nlevp1
    
    IF (ltimer) CALL timer_start(timer_rht)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_rht    => echam_phy_config(jg)%fc_rht
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       nlevp1 = nlev+1
       !
       CALL radheating (                                   &
            !
            ! input
            ! -----
            !
            & jcs        = jcs                            ,&! loop start index
            & jce        = jce                            ,&! loop end index
            & kbdim      = nproma                         ,&! dimension size
            & klev       = nlev                           ,&! vertical dimension size
            & klevp1     = nlevp1                         ,&! vertical dimension size
            !
            & rsdt0      = psctm                          ,&! toa incident shortwave radiation for sun in zenith
            & cosmu0     = field%cosmu0    (:,jb)         ,&! solar zenith angle at current time
            & daylght_frc= field%daylght_frc(:,jb)        ,&! daylight fraction
            !
            & emiss      = ext_data(jg)%atm%emis_rad(:,jb),&! lw sfc emissivity
            & tsr        = field%ts_rad (:,jb)            ,&! radiative surface temperature at current   time [K]
            & tsr_rt     = field%ts_rad_rt(:,jb)          ,&! radiative surface temperature at radiation time [K]
            !
            & rsd_rt     = field%rsd_rt           (:,:,jb),&! all-sky   shortwave downward flux at radiation time [W/m2]
            & rsu_rt     = field%rsu_rt           (:,:,jb),&! all-sky   shortwave upward   flux at radiation time [W/m2]
            !
            & rsdcs_rt   = field%rsdcs_rt         (:,:,jb),&! clear-sky shortwave downward flux at radiation time [W/m2]
            & rsucs_rt   = field%rsucs_rt         (:,:,jb),&! clear-sky shortwave upward   flux at radiation time [W/m2]
            !
            & rld_rt     = field%rld_rt           (:,:,jb),&! all-sky   longwave  downward flux at radiation time [W/m2]
            & rlu_rt     = field%rlu_rt           (:,:,jb),&! all-sky   longwave  upward   flux at radiation time [W/m2]
            !
            & rldcs_rt   = field%rldcs_rt         (:,:,jb),&! clear-sky longwave  downward flux at radiation time [W/m2]
            & rlucs_rt   = field%rlucs_rt         (:,:,jb),&! clear-sky longwave  upward   flux at radiation time [W/m2]
            !
            & rvds_dir_rt= field%rvds_dir_rt        (:,jb),&!< out  all-sky downward direct visible radiation at surface
            & rpds_dir_rt= field%rpds_dir_rt        (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
            & rnds_dir_rt= field%rnds_dir_rt        (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
            & rvds_dif_rt= field%rvds_dif_rt        (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
            & rpds_dif_rt= field%rpds_dif_rt        (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
            & rnds_dif_rt= field%rnds_dif_rt        (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
            & rvus_rt    = field%rvus_rt            (:,jb),&!< out  all-sky upward visible radiation at surface
            & rpus_rt    = field%rpus_rt            (:,jb),&!< out  all-sky upward PAR     radiation at surfac
            & rnus_rt    = field%rnus_rt            (:,jb),&!< out  all-sky upward near-IR radiation at surface
            !
            ! output
            ! ------
            !
            & rsdt       = field%rsdt               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
            & rsut       = field%rsut               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
            & rsds       = field%rsds               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
            & rsus       = field%rsus               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
            !
            & rsutcs     = field%rsutcs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
            & rsdscs     = field%rsdscs             (:,jb),&! clear-sky shortwave downward flux at current   time [W/m2]
            & rsuscs     = field%rsuscs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
            !
            & rvds_dir   = field%rvds_dir           (:,jb),&!< out  all-sky downward direct visible radiation at surface
            & rpds_dir   = field%rpds_dir           (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
            & rnds_dir   = field%rnds_dir           (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
            & rvds_dif   = field%rvds_dif           (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
            & rpds_dif   = field%rpds_dif           (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
            & rnds_dif   = field%rnds_dif           (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
            & rvus       = field%rvus               (:,jb),&!< out  all-sky upward visible radiation at surface
            & rpus       = field%rpus               (:,jb),&!< out  all-sky upward PAR     radiation at surfac
            & rnus       = field%rnus               (:,jb),&!< out  all-sky upward near-IR radiation at surface
            !
            & rlut       = field%rlut               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
            & rlds       = field%rlds               (:,jb),&! all-sky   longwave  downward flux at current   time [W/m2]
            & rlus       = field%rlus               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
            !
            & rlutcs     = field%rlutcs             (:,jb),&! clear-sky longwave  upward   flux at current   time [W/m2]
            & rldscs     = field%rldscs             (:,jb),&! clear-sky longwave  downward flux at current   time [W/m2]
            !
            & q_rsw      = field%q_rsw            (:,:,jb),&! rad. heating by SW           [W/m2]
            & q_rlw      = field%q_rlw            (:,:,jb) )! rad. heating by LW           [W/m2]
       !
       ! vertical integral
       field% q_rsw_vi(jcs:jce,jb) = SUM(field% q_rsw(jcs:jce,:,jb),DIM=2)
       field% q_rlw_vi(jcs:jce,jb) = SUM(field% q_rlw(jcs:jce,:,jb),DIM=2)
       !
       ! convert    heating
       tend%ta_rsw(jcs:jce,:,jb) = field%q_rsw(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       tend%ta_rlw(jcs:jce,:,jb) = field%q_rlw(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field%q_rsw(jcs:jce,:,jb) + field%q_rlw(jcs:jce,:,jb)
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_rht)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend% ta_phy(jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb) + tend% ta_rsw (jcs:jce,:,jb) + tend% ta_rlw (jcs:jce,:,jb)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field% ta(jcs:jce,:,jb) = field% ta(jcs:jce,:,jb) + (tend% ta_rsw(jcs:jce,:,jb) + tend% ta_rlw (jcs:jce,:,jb))*pdtime
       END IF
       !
    ELSE
       !
       tend%ta_rsw(jcs:jce,:,jb) = 0.0_wp
       tend%ta_rlw(jcs:jce,:,jb) = 0.0_wp
       !
    END IF

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_rht)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_rht)

  END SUBROUTINE interface_echam_rht

END MODULE mo_interface_echam_rht
