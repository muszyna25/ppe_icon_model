!>
!! Provides utilities for abstracting the NWP and ECHAM physics differenes
!!
!! This module provides an interface to ART routines.  The interface
!! is written in such a way, that ICON will compile and run properly,
!! even if the ART-routines are not available at compile time.
!!
!!
!! @author Luis Kornblueh, MPIM
!!
!! @par Revision History
!! Initial revision by Luis Kornblueh, MPIM (2019-10-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_util

  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: iecham, inwp
  USE mo_echam_phy_memory,              ONLY: prm_field
  USE mo_nwp_phy_state,                 ONLY: prm_diag  
  USE mo_physical_constants,            ONLY: grav
  
  IMPLICIT NONE

  PRIVATE

  TYPE t_art_phys_container
    REAL(wp), POINTER :: swflx_par_sfc(:,:) => NULL()
    REAL(wp), POINTER :: rpds_dir(:,:)      => NULL()
    REAL(wp), POINTER :: rpds_dif(:,:)      => NULL()    
    REAL(wp), POINTER :: u_10m(:,:)         => NULL()
    REAL(wp), POINTER :: v_10m(:,:)         => NULL()
    REAL(wp), POINTER :: rain_gsp_rate(:,:) => NULL()
    REAL(wp), POINTER :: rain_con_rate(:,:) => NULL()
    REAL(wp), POINTER :: rh_2m(:,:)         => NULL()
    REAL(wp), POINTER :: t_2m(:,:)          => NULL()
    REAL(wp), POINTER :: td_2m(:,:)         => NULL()    
    REAL(wp), POINTER :: tcm(:,:)           => NULL()
    REAL(wp), POINTER :: tch(:,:)           => NULL()
    REAL(wp), POINTER :: gz0(:,:)           => NULL()
    REAL(wp), POINTER :: z0m(:,:)           => NULL()            
  CONTAINS
    PROCEDURE :: fill =>  t_art_phys_container_fill
  END TYPE t_art_phys_container

  PUBLIC :: t_art_phys_container
  PUBLIC :: td_and_t2rh
  
CONTAINS

  ! Variable mapping
  !
  ! echam prm_fields    | nwp prm_diag  | description [unit]
  !--------------------------------------------------------------------------------------------------------------                      
  ! rpds_dir + rpds_dif | swflx_par_sfc | shortwave downward photosynthetically active flux at the surface [W/m2]
  ! uas                 | u_10m         | zonal wind in 10m [m s-1]
  ! vas                 | v_10m         | meridional wind in 10m [m s-1]
  ! rsfl                | rain_gsp_rate | gridscale rain rate [kg m-2 s-1]
  ! rsfc                | rain_con_rate | convective rain rate [kg m-2 s-1]
  ! grav*z0m            | gz0           | grid cell average z0 * gravity
  ! cfh                 | tch           | turbulent transfer coefficients for heat [?]
  ! cfm                 | tcm           | turbulent transfer coefficients for momentum [?]
  ! dew2, tas -> Magnus | rh_2m         | relative humidity in 2m [%]
  !
  ! Magnus formula (coefficients from Sonntag, 1990)
  
  ELEMENTAL FUNCTION td_and_t2rh(td, t) RESULT (rh)
    REAL(wp) :: rh                ! result in [%]  
    REAL(wp), INTENT(in) :: td    ! input expected in [K]
    REAL(wp), INTENT(in) :: t     ! input expected in [K]
    REAL(wp), PARAMETER :: b =  17.62_wp
    REAL(wp), PARAMETER :: c = -30.03_wp
    REAL(wp), PARAMETER :: t0 = 273.15_wp

    rh = 100.0_wp*(exp((b*(td-t0))/(c+td))/exp((b*(t-t0))/(c+t)))

  END FUNCTION td_and_t2rh

  SUBROUTINE t_art_phys_container_fill(this, physic, jd)
    CLASS(t_art_phys_container) :: this
    INTEGER, INTENT(in) :: physic
    INTEGER, INTENT(in), OPTIONAL :: jd

    INTEGER :: jdl
    
    jdl = 1
    IF (PRESENT(jd)) jdl = jd

    SELECT CASE (physic)
    CASE (inwp)
      this%swflx_par_sfc => prm_diag(jdl)%swflx_par_sfc
      this%u_10m         => prm_diag(jdl)%u_10m           
      this%v_10m         => prm_diag(jdl)%v_10m           
      this%rain_gsp_rate => prm_diag(jdl)%rain_gsp_rate   
      this%rain_con_rate => prm_diag(jdl)%rain_con_rate   
      this%rh_2m         => prm_diag(jdl)%rh_2m           
      this%tch           => prm_diag(jdl)%tch
      this%tcm           => prm_diag(jdl)%tcm
      this%gz0           => prm_diag(jdl)%gz0                 
    CASE (iecham)
      this%rpds_dir      => prm_field(jdl)%rpds_dir
      this%rpds_dif      => prm_field(jdl)%rpds_dif      
      this%u_10m         => prm_field(jdl)%uas
      this%v_10m         => prm_field(jdl)%vas           
      this%rain_gsp_rate => prm_field(jdl)%rsfl
      this%rain_con_rate => prm_field(jdl)%rsfc
      this%t_2m          => prm_field(jdl)%tas
      this%td_2m         => prm_field(jdl)%dew2
      this%tcm           => prm_field(jdl)%cfm(:,ubound(prm_field(jdl)%cfm,2),:)                 
      this%tch           => prm_field(jdl)%cfh(:,ubound(prm_field(jdl)%cfh,2),:)                 
      this%z0m           => prm_field(jdl)%z0m                
    END SELECT

  END SUBROUTINE t_art_phys_container_fill

END MODULE mo_art_util


