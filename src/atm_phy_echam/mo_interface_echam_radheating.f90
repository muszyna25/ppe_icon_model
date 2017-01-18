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
  USE mo_run_config,          ONLY: nlev, nlevp1
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_radheat
  USE mo_radheating,          ONLY: radheating
  USE mo_psrad_radiation_parameters, ONLY: psctm

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_radheating

CONTAINS

  !----------------------------------------------------------------
  SUBROUTINE interface_echam_radheating(patch, rl_start, rl_end, field, tend, zconv, zq_rlw, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER :: field    
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp) :: zconv  (nproma,nlev,patch%nblks_c)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zq_rlw (nproma,nlev,patch%nblks_c)       !< heating by long  wave radiation   [W/m2]
    REAL(wp) :: zq_phy (nproma,nlev,patch%nblks_c)       !< heating by whole ECHAM physics    [W/m2]

    REAL(wp) :: zq_rsw (nproma,nlev)       !< heating by short wave radiation   [W/m2]
    INTEGER  :: jg             
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
    

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    IF (phy_config%lrad) THEN
      
      IF (ltimer) CALL timer_start(timer_radheat)
      ! 4.2 RADIATIVE HEATING
      !----------------------
      ! radheat first computes the shortwave and longwave radiation for the current time step from transmissivity and
      ! the longwave flux at the radiation time step and, from there, the radiative heating due to sw and lw radiation.
      ! If radiation is called every time step, the longwave flux is not changed.

!$OMP PARALLEL DO PRIVATE(jcs,jce, zq_rsw)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        CALL echam_radheating( jg, jb,jcs,jce, nproma, field, zq_rsw, zq_rlw(:,:,jb))

        ! heating accumulated
        zq_phy(jcs:jce,:,jb) = zq_phy(jcs:jce,:,jb) + zq_rsw(jcs:jce,:) + zq_rlw(jcs:jce,:,jb)
        
        ! tendencies
        tend%ta_rsw(jcs:jce,:,jb) = zq_rsw(jcs:jce,:)    * zconv(jcs:jce,:,jb)
        tend%ta_rlw(jcs:jce,:,jb) = zq_rlw(jcs:jce,:,jb) * zconv(jcs:jce,:,jb)

        ! tendencies accumulated
        tend% ta(jcs:jce,:,jb) = tend% ta     (jcs:jce,:,jb) &
          &                    + tend% ta_rsw (jcs:jce,:,jb) &
          &                    + tend% ta_rlw (jcs:jce,:,jb)


      ENDDO
!$OMP END PARALLEL DO 
      IF (ltimer) CALL timer_stop(timer_radheat)

    ELSE   
      ! If computation of radiative heating is by-passed
      ! this shound not be needed, as these field are zeroed in the initialization
!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        tend%ta_rsw(jcs:jce,:,jb) = 0.0_wp
        tend%ta_rlw(jcs:jce,:,jb) = 0.0_wp

        field%rsdt(jcs:jce,jb)= 0.0_wp
      ENDDO
!$OMP END PARALLEL DO 

    END IF ! lrad
 
  END SUBROUTINE interface_echam_radheating
   !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_radheating( jg, jb,jcs,jce, nbdim, field, zq_rsw, zq_rlw)
    INTEGER         ,INTENT(IN) :: jg
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field    ! in

    REAL(wp)        ,INTENT(INOUT) :: zq_rsw (nbdim,nlev)       !< out, heating by short wave radiation   [W/m2]
    REAL(wp)        ,INTENT(INOUT) :: zq_rlw (nbdim,nlev)      !<  out, heating by long  wave radiation   [W/m2]

    CALL radheating (                                &
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
      & q_rsw      = zq_rsw                   (:,:) ,&! rad. heating by SW           [W/m2]
      & q_rlw      = zq_rlw                   (:,:) ) ! rad. heating by LW           [W/m2]

  END SUBROUTINE echam_radheating
  !---------------------------------------------------------------------

END MODULE mo_interface_echam_radheating
