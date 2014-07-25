#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!!  !MODULE:  mo_nh_initicon_types\\
!!
!! Description:  Data type definition for initicon
!
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-13)
!! - encapsulated type definitions by Daniel Reinert, DWD (2012-12-17)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nh_initicon_types

  USE mo_kind,                 ONLY: wp
  USE mo_var_metadata_types,   ONLY: VARNAME_LEN
  USE mtime,                   ONLY: datetime, timedelta

  IMPLICIT NONE
  PRIVATE


  !
  !variables
  PUBLIC :: t_initicon_state  !> state vector for initicon
  PUBLIC :: t_pi_atm_in
  PUBLIC :: t_pi_sfc_in
  PUBLIC :: t_pi_atm
  PUBLIC :: t_pi_sfc


  ! atmospheric input variables
  TYPE :: t_pi_atm_in ! surface geopotential is regarded as
                      ! atmospheric variable here because the atmospheric fields cannot be processed without it

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: psfc, phi_sfc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: temp, pres, z3d_ifc, w_ifc, z3d, u, v, omega, &
      &                                         w, vn, qv, qc, qi, qr, qs

  END TYPE t_pi_atm_in


  ! surface input variables
  TYPE :: t_pi_sfc_in

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, sst, snowalb,snowweq, snowdens, &
                                              skinres, ls_mask, seaice, phi
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, wsoil

  END TYPE t_pi_sfc_in


  ! 
  TYPE :: t_pi_atm

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: vn, u, v, w, temp, theta_v, exner, rho, &
                                               pres, qv, qc, qi, qr, qs

  END TYPE t_pi_atm


  !
  TYPE :: t_pi_sfc

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, sst,  snowalb, snowweq, snowdens, &
                                              skinres, ls_mask, seaice
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, wsoil

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: w_so

  END TYPE t_pi_sfc



  ! complete state vector type
  !
  TYPE :: t_initicon_state
    
    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: topography_c

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: z_ifc, z_mc

    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE, DIMENSION(:)    :: grp_vars_fg, grp_vars_ana
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE, DIMENSION(:)    :: grp_vars_fg_default, grp_vars_ana_default

    INTEGER :: ngrp_vars_fg, ngrp_vars_ana, ngrp_vars_fg_default, ngrp_vars_ana_default

    TYPE (t_pi_atm_in)     :: atm_in
    TYPE (t_pi_sfc_in)     :: sfc_in
    TYPE (t_pi_atm)        :: atm
    TYPE (t_pi_atm)        :: atm_inc
    TYPE (t_pi_sfc)        :: sfc
    TYPE (t_pi_sfc)        :: sfc_inc

  END TYPE t_initicon_state
 

END MODULE mo_nh_initicon_types
