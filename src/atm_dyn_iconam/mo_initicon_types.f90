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
MODULE mo_initicon_types

  USE mo_kind,                 ONLY: wp
  USE mo_var_metadata_types,   ONLY: VARNAME_LEN
  USE mo_dictionary,           ONLY: t_dictionary

  IMPLICIT NONE
  PRIVATE


  !
  !variables
  PUBLIC :: t_initicon_state  !> state vector for initicon
  PUBLIC :: t_pi_atm_in
  PUBLIC :: t_pi_sfc_in
  PUBLIC :: t_pi_atm
  PUBLIC :: t_pi_sfc
  PUBLIC :: t_sfc_inc
  PUBLIC :: geop_ml_var, alb_snow_var
  PUBLIC :: ana_varnames_dict

  ! atmospheric input variables
  TYPE :: t_pi_atm_in ! surface geopotential is regarded as
                      ! atmospheric variable here because the atmospheric fields cannot be processed without it

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), POINTER, DIMENSION(:,:)   :: psfc, phi_sfc
    REAL(wp), POINTER, DIMENSION(:,:,:) :: temp, pres, z3d_ifc, w_ifc, z3d, u, v, omega, &
      &                                         w, vn, qv, qc, qi, qr, qs, rho, theta_v, tke, tke_ifc

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
                                               pres, qv, qc, qi, qr, qs, tke

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


  ! surface field increments
  TYPE :: t_sfc_inc
    !
    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: w_so
    REAL(wp), ALLOCATABLE, DIMENSION (:,:)   :: h_snow
    REAL(wp), ALLOCATABLE, DIMENSION (:,:)   :: freshsnow

  END TYPE t_sfc_inc


  ! complete state vector type
  !
  TYPE :: t_initicon_state
    
    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: topography_c

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: z_ifc, z_mc

    TYPE (t_pi_atm_in)     :: atm_in
    TYPE (t_pi_sfc_in)     :: sfc_in
    TYPE (t_pi_atm)        :: atm
    TYPE (t_pi_atm)        :: atm_inc
    TYPE (t_pi_sfc)        :: sfc
    TYPE (t_sfc_inc)       :: sfc_inc

  END TYPE t_initicon_state
 

  CHARACTER(LEN=VARNAME_LEN) :: geop_ml_var  ! model level surface geopotential
  CHARACTER(LEN=VARNAME_LEN) :: alb_snow_var ! snow albedo


  ! dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  TYPE (t_dictionary) :: ana_varnames_dict

END MODULE mo_initicon_types
