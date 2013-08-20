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
!! $Id: n/a$
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-13)
!! - encapsulated type definitions by Daniel Reinert, DWD (2012-12-17)
!!
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_nh_initicon_types

  USE mo_kind,                 ONLY: wp
  USE mo_var_metadata,         ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: temp, pres, z3d, u, v, omega, &
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

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: vn, u, v, w, temp, theta_v, exner, rho, &
                                               pres, qv, qc, qi, qr, qs

  END TYPE t_pi_atm


  !
  TYPE :: t_pi_sfc

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, sst,  snowalb, snowweq, snowdens, &
                                              skinres, ls_mask, seaice
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, wsoil

    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE, DIMENSION(:)    :: grp_vars_fg, grp_vars_ana
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE, DIMENSION(:)    :: grp_vars_fg_default, grp_vars_ana_default

    INTEGER :: ngrp_vars_fg, ngrp_vars_ana, ngrp_vars_fg_default, ngrp_vars_ana_default

  END TYPE t_pi_sfc



  ! complete state vector type
  !
  TYPE :: t_initicon_state
    
    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: topography_c, topography_v

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: z_ifc, z_mc

    TYPE (t_pi_atm_in) :: atm_in
    TYPE (t_pi_sfc_in) :: sfc_in
    TYPE (t_pi_atm)    :: atm
    TYPE (t_pi_sfc)    :: sfc

  END TYPE t_initicon_state
 

END MODULE mo_nh_initicon_types
