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
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_interface_echam_radheating

  USE mo_kind,                ONLY: wp

  USE mo_model_domain        ,ONLY: t_patch
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config,          ONLY: nlev, nlevp1
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, t_echam_phy_tend
  USE mo_ext_data_state,      ONLY: ext_data

  USE mo_psrad_solar_parameters, ONLY: psctm
  USE mo_radheating,          ONLY: radheating

  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_radheat

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_radheating

CONTAINS

  !-------------------------------------------------------------------
  SUBROUTINE interface_echam_radheating(is_in_sd_ed_interval,    &
       &                                patch, rl_start, rl_end, &
       &                                field, tend              )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    TYPE(t_patch)   ,TARGET ,INTENT(in)    :: patch           !< grid/patch info
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

      
    IF (ltimer) CALL timer_start(timer_radheat)
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
       !
       CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
       !
       CALL echam_radheating(is_in_sd_ed_interval,          &
            &                jg, jb,jcs,jce, nproma,        &
            &                field, tend                    )

    END DO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_radheat)
 
  END SUBROUTINE interface_echam_radheating
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE echam_radheating(is_in_sd_ed_interval,  &
       &                      jg, jb,jcs,jce, nbdim, &
       &                      field, tend            )

    LOGICAL                 ,INTENT(in)    :: is_in_sd_ed_interval
    INTEGER                 ,INTENT(in)    :: jg
    INTEGER                 ,INTENT(in)    :: jb                  !< block index
    INTEGER                 ,INTENT(in)    :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in)    :: nbdim               !< size of this block
    TYPE(t_echam_phy_field) ,POINTER       :: field
    TYPE(t_echam_phy_tend)  ,POINTER       :: tend

    IF ( is_in_sd_ed_interval ) THEN
       !
       CALL radheating (                                   &
            !
            ! input
            ! -----
            !
            & jcs        = jcs                            ,&! loop start index
            & jce        = jce                            ,&! loop end index
            & kbdim      = nbdim                          ,&! dimension size
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
       ! accumulate tendencies
       tend% ta_phy(jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb) + tend% ta_rsw (jcs:jce,:,jb) + tend% ta_rlw (jcs:jce,:,jb)
       !
    ELSE
       !
       tend%ta_rsw(jcs:jce,:,jb) = 0.0_wp
       tend%ta_rlw(jcs:jce,:,jb) = 0.0_wp
       !
    END IF
    !
  END SUBROUTINE echam_radheating
  !---------------------------------------------------------------------

END MODULE mo_interface_echam_radheating
