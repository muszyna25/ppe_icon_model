!>
!> Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------
!!
! Note: This file shall be called "mo_phyparam_turb.f90" once changed in COSMO.
!
!> Description:
!!  This module contains variables that are used in the turbulence
!!  parameterizations. With these variables a tuning of the scheme is
!!  possible.
!!
!! @author DWD, Ulrich.Schaettler
!!  phone:  +49  69  8062 2739
!!  fax:    +49  69  8062 3721
!!  email:  ulrich.schaettler@dwd.de
!!
!!
!! @par Revision History
!> ---------- ---------- ----
!! 3.21       2006/12/04 Ulrich Schaettler
!!  Initial release
!! V3_23        2007/03/30 Matthias Raschendorfer
!!  Importing 'rat_lam' from data_soil
!!  and 'clc_diag', 'q_crit', 'akt' from data_constants.
!!  Introduction of some parameters from turb_param.incf.
!!  Initialization of all parameters with default values.
!! V4_10        2009/09/11 Matthias Raschendorfer
!!  Introduction of 'a_hshr' and 'a_stab'.
!! V4_13        2010/05/11 Michael Gertz
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

MODULE data_turbulence

! Modules used:

#ifdef __COSMO__

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

#endif

#ifdef __ICON__

USE mo_kind, ONLY:     &
    ireals   =>wp    , &
    iintegers=>i4

#endif

!==============================================================================

IMPLICIT NONE
PUBLIC           ! All constants and variables in this module are public

!==============================================================================


! Variables for tuning the turbulence parameterizations
! -----------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs!

!==============================================================================

! 1. Parameters describing physical properties of the lower boundary
!    of the atmosphere:
!------------------------------------------

REAL (KIND=ireals) ::          &
  rlam_mom   =  0.0_ireals,    & ! scaling factor of the laminar boundary layer for momentum
  rlam_heat  =  1.0_ireals,    & ! scaling factor of the laminar boundary layer for heat

  rat_lam    =  1.0_ireals,    & ! ratio of laminar scaling factors for vapour and heat
  rat_can    =  1.0_ireals,    & ! ratio of canopy height over z0m
  rat_sea    =  1.0_ireals,    & ! ratio of laminar scaling factors for heat over sea and land

  z0m_dia    =  0.2_ireals,    & ! roughness length of a typical synoptic station [m]

  alpha0     =  0.0123_ireals, & ! Charnock-parameter
  alpha1     =  0.0000_ireals    ! parameter scaling the molek. roughness of water waves


! 2. Parameters that should be external parameter fields being not yet
!    available:
!------------------------------------------

REAL (KIND=ireals) ::          &
  c_lnd      = 2.0_ireals,     & ! surface area density of the roughness elements over land
  c_sea      = 1.5_ireals,     & ! surface area density of the waves over sea
  c_soil     = 1.0_ireals,     & ! surface area density of the (evaporative) soil surface
  e_surf     = 1.0_ireals        ! exponent to get the effective surface area


! 3. Parameters that should be dynamical fields being not yet available:
!------------------------------------------

REAL (KIND=ireals) ::          &
  zt_ice     = -1.7_ireals,    & !freezing temperature of sea ice
  z0_ice     =  0.001_ireals     !roughness length of sea ice


! 4. Parameters for modelling turbulent diffusion:
!------------------------------------------

REAL (KIND=ireals) ::         &
  tur_len    = 500.0_ireals,  & ! assymtotic maximal turbulent length scale [m]
  pat_len    = 500.0_ireals,  & ! lenth scale of subscale surface patterns over land [m]
                                ! (should be dependent on location)
  len_min    =  1.0E-6_ireals,& ! minimal turbulent length scale [m]

  vel_min    =  0.01_ireals,  & ! minimal velocity scale [m/s]

  akt        =  0.4_ireals,   & ! von Karman-constant

  ! Length scale factors for pressure destruction of turbulent
  a_heat     =  0.74_ireals,  & ! scalar (heat) transport
  a_mom      =  0.92_ireals,  & ! momentum transport

  ! Length scale factors for dissipation of
  d_heat     =  10.1_ireals,  & ! scalar (temperature) variance
  d_mom      =  16.6_ireals,  & ! momentum variance

  ! Length scale factors for turbulent transport (vertical diffusion)
  c_diff     =  0.20_ireals,  & ! of TKE

  ! Length scale factor for separate horizontal shear production
  a_hshr     =  0.20_ireals,  & ! of TKE

  ! Length scale factor for the stability correction
  a_stab     =  0.00_ireals,  & ! no stability correction so far

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):
  clc_diag   =  0.5_ireals,   & !cloud cover at saturation
  q_crit     =  4.0_ireals,   & !critical value for normalized over-saturation
  c_scld     =  1.0_ireals,   & !factor for liquid water flux density in sub grid scale clouds

  ! Minimal diffusion coefficients in [m^2/s] for vertical
  tkhmin     =  0.1_ireals,   & ! scalar (heat) transport
  tkmmin     =  0.1_ireals      ! momentum transport


! 5. Numerical parameters:
!-------------------------

REAL (KIND=ireals) ::         &
  epsi       =  1.0E-6_ireals,& ! relative limit of accuracy for comparison of numbers
  tkesmot    =  0.15_ireals,  & ! time smoothing factor for TKE and diffusion coefficients
  wichfakt   =  0.15_ireals,  & ! vertical smoothing factor for explicit diffusion tendencies
  securi     =  0.50_ireals     ! security factor for maximal diffusion coefficients

INTEGER (KIND=iintegers) :: &
  it_end     =  1                ! number of initialization iterations (>=0)

!==============================================================================

END MODULE data_turbulence
