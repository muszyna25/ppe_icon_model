!>
!!  Data module for variables concerned with the FLake Model
!!
!!
!! Description:
!!  This module contains common definitions, variables and constants for the
!!  FLake Model. These are:
!!
!!  - flake_albedo_ref:
!!   "reference" values of albedo for the lake water, lake ice and snow.
!!    As in "flake_paramoptic_ref", two ice categories, viz. white ice and blue
!!    ice, and two snow categories, viz. dry snow and melting snow, are used.
!!
!!  - flake_configure:
!!    Switches and reference values of parameters
!!    that configure the lake model FLake are set.
!!
!!  - flake_derivedtypes:
!!    Derived type(s) is(are) defined.
!!
!!  - flake_parameters:
!!    Values of empirical constants of the lake model FLake
!!    and of several thermodynamic parameters are set.
!!
!!  - flake_paramoptic_ref:
!!    This module contains "reference" values of the optical characteristics
!!    of the lake water, lake ice and snow.
!!
!!
!! @author Dmitrii Mironov, DWD
!!
!!
!! @par Revision History
!! ---------- ---------- ----
!! 3.18       2006/03/03 Dmitrii Mironov
!!  Initial release
!! V4_13        2010/05/11 Michael Gertz
!!  Adaptions to SVN
!!
!! Implementation into ICON Kristina Froehlich
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
!!

MODULE data_flake

! Modules used:

#ifdef __COSMO__
USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

#endif

#ifdef __ICON__

USE mo_kind,               ONLY: ireals=>wp     , &
                                 iintegers=>i4
#endif
!==============================================================================

IMPLICIT NONE
PUBLIC           ! All constants and variables in this module are public

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
!==============================================================================
!
! Declarations

!  flake_albedo_ref:
!  -----------------

!  Albedo for water, ice and snow.
REAL (KIND = ireals), PARAMETER ::        &
  albedo_water_ref       = 0.07_ireals  , & ! Water
  albedo_whiteice_ref    = 0.60_ireals  , & ! White ice
  albedo_blueice_ref     = 0.10_ireals  , & ! Blue ice
  albedo_drysnow_ref     = 0.60_ireals  , & ! Dry snow
  albedo_meltingsnow_ref = 0.10_ireals      ! Melting snow

!  Empirical parameters.
REAL (KIND = ireals), PARAMETER :: &
  c_albice_MR = 95.6_ireals        ! Constant in the interpolation formula for
                                   ! the ice albedo (Mironov and Ritter 2004)

!------------------------------------------------------------------------------

!  flake_configure:
!  ----------------

! Switches and reference values of parameters
LOGICAL, PARAMETER :: &
  lflk_botsed_use   = .FALSE.

  ! .TRUE. indicates that the bottom-sediment scheme is used to compute the
  ! depth penetrated by the thermal wave, the temperature at this depth and
  ! the bottom heat flux.  Otherwise, the heat flux at the water-bottom
  ! sediment interface is set to zero, the depth penetrated by the thermal
  ! wave is set to a reference value defined below, and the temperature at
  ! this depth is set to the temperature of maximum density of the fresh water.

REAL (KIND = ireals), PARAMETER :: &
  rflk_depth_bs_ref = 10.0_ireals

  ! Reference value of the depth of the thermally active layer of bottom
  ! sediments [m].  This value is used to (formally) define the depth
  ! penetrated by the thermal wave in case the bottom-sediment scheme
  ! is not used.

!------------------------------------------------------------------------------

!  flake_derived_types:
!  --------------------

!  Derived type(s) is(are) defined.

!  Maximum value of the wave-length bands
!  in the exponential decay law for the radiation flux.
!  A storage for a ten-band approximation is allocated,
!  although a smaller number of bands is actually used.

INTEGER (KIND = iintegers), PARAMETER :: &
  nband_optic_max = 10_iintegers

!  Define TYPE "opticpar_medium"
TYPE t_opticpar_medium

  INTEGER (KIND = iintegers)          :: &
    nband_optic                            ! Number of wave-length bands

  REAL (KIND = ireals)                :: &
    frac_optic       (nband_optic_max),  & ! Fractions of total radiation flux
    extincoef_optic  (nband_optic_max)     ! Extinction coefficients

END TYPE t_opticpar_medium

!------------------------------------------------------------------------------

! flake_parameters:
! -----------------

!  Values of empirical constants of the lake model FLake
!  and of several thermodynamic parameters are set.

!  Dimensionless constants in the equations for the mixed-layer depth and for
!  the shape factor with respect to the temperature profile in the thermocline

REAL (KIND = ireals), PARAMETER ::      &
  c_cbl_1    = 0.17_ireals , & ! Constant in the CBL entrainment equation
  c_cbl_2    = 1.0_ireals  , & ! Constant in the CBL entrainment equation
  c_sbl_ZM_n = 0.5_ireals  , & ! Constant in the ZM1996 equation for the
                               !       equilibrium SBL depth
  c_sbl_ZM_s = 10.0_ireals , & ! Constant in the ZM1996 equation for the
                               !       equilibrium SBL depth
  c_sbl_ZM_i = 20.0_ireals , & ! Constant in the ZM1996 equation for the
                               !       equilibrium SBL depth
  c_relax_h  = 0.010_ireals, & ! Constant in the relaxation equation for the
                               !       SBL depth
  c_relax_C  = 0.0030_ireals   ! Constant in the relaxation equation for the
                               !       shape factor with respect to the
                               !       temperature profile in the thermocline

!  Parameters of the shape functions
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments,
!  B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function
!  at "zeta=0" ad "zeta=1", respectively.

REAL (KIND = ireals), PARAMETER ::         &
  C_T_min       = 0.5_ireals           , & ! Minimum value of the shape factor
                                           !    C_T (thermocline)
  C_T_max       = 0.8_ireals           , & ! Maximum value of the shape factor
                                           !    C_T (thermocline)
  Phi_T_pr0_1   = 40._ireals/3._ireals , & ! Constant in the expression for the
                                           !    T shape-function derivative
  Phi_T_pr0_2   = 20._ireals/3._ireals , & ! Constant in the expression for the
                                           !    T shape-function derivative
  C_TT_1        = 11._ireals/18._ireals, & ! Constant in the expression for
                                           !    C_TT (thermocline)
  C_TT_2        = 7._ireals/45._ireals , & ! Constant in the expression for
                                           !    C_TT (thermocline)
  C_B1          = 2._ireals/3._ireals  , & ! Shape factor (upper layer of
                                           !    bottom sediments)
  C_B2          = 3._ireals/5._ireals  , & ! Shape factor (lower layer of
                                           !    bottom sediments)
  Phi_B1_pr0    = 2._ireals            , & ! B1 shape-function derivative
                                           !
  C_S_lin       = 0.5_ireals           , & ! Shape factor (linear temperature
                                           !    profile in the snow layer)
  Phi_S_pr0_lin = 1._ireals            , & ! S shape-function derivative
                                           !    (linear profile)
  C_I_lin       = 0.5_ireals           , & ! Shape factor (linear temperature
                                           !    profile in the ice layer)
  Phi_I_pr0_lin = 1._ireals            , & ! I shape-function derivative
                                           !    (linear profile)
  Phi_I_pr1_lin = 1._ireals            , & ! I shape-function derivative
                                           !    (linear profile)
  Phi_I_ast_MR  = 2._ireals            , & ! Constant in the MR2004 expression
                                           !    for I shape factor
  C_I_MR        = 1._ireals/12._ireals , & ! Constant in the MR2004 expression
                                           !    for I shape factor
  H_Ice_max     = 3._ireals                ! Maximum ice tickness in the
                                           !    Mironov and Ritter ice model [m]

!  Security constants
REAL (KIND = ireals), PARAMETER ::       &
  h_Snow_min_flk = 1.0E-5_ireals       , & ! Minimum snow thickness [m]
  h_Ice_min_flk  = 1.0E-9_ireals       , & ! Minimum ice thickness [m]
  h_ML_min_flk   = 1.0E-2_ireals       , & ! Minimum mixed-layer depth [m]
  h_ML_max_flk   = 1.0E+3_ireals       , & ! Maximum mixed-layer depth [m]
  H_B1_min_flk   = 1.0E-3_ireals       , & ! Minimum thickness of the upper
                                           !     layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6_ireals           ! Minimum value of the surface
                                           !     friction velocity [m s^{-1}]

!  Security constant(s)
REAL (KIND = ireals), PARAMETER ::         &
  c_small_flk    = 1.0E-10_ireals            ! A small number

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::        &
  tpl_grav          = 9.81_ireals      , & ! Acceleration due to gravity
                                           !        [m s^{-2}]
  tpl_T_r           = 277.13_ireals    , & ! Temperature of maximum density
                                           !        of fresh water [K]
  tpl_T_f           = 273.15_ireals    , & ! Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05_ireals, & ! Constant in the fresh-water
                                           !        equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03_ireals   , & ! Maximum density of fresh water
                                           !        [kg m^{-3}]
  tpl_rho_I         = 9.1E+02_ireals   , & ! Density of ice [kg m^{-3}]
  tpl_rho_S_min     = 1.0E+02_ireals   , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max     = 4.0E+02_ireals   , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S   = 2.0E+02_ireals   , & ! Empirical parameter [kg m^{-4}] in
                                           ! the expression for the snow density
  tpl_L_f           = 3.3E+05_ireals   , & ! Latent heat of fusion
                                           !   [J kg^{-1}]
  tpl_c_w           = 4.2E+03_ireals   , & ! Specific heat of water
                                           !   [J kg^{-1} K^{-1}]
  tpl_c_I           = 2.1E+03_ireals   , & ! Specific heat of ice
                                           !   [J kg^{-1} K^{-1}]
  tpl_c_S           = 2.1E+03_ireals   , & ! Specific heat of snow
                                           !   [J kg^{-1} K^{-1}]
  tpl_kappa_w       = 5.46E-01_ireals  , & ! Molecular heat conductivity of water
                                           !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_I       = 2.29_ireals      , & ! Molecular heat conductivity of ice
                                           !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_min   = 0.2_ireals       , & ! Minimum molecular heat conductivity
                                           !   of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max   = 1.5_ireals       , & ! Maximum molecular heat conductivity
                                           !   of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S = 1.3_ireals           ! Empirical parameter in expression
                                           !   for the snow heat conductivity
                                           !   [J m^{-2} s^{-1} K^{-1}]

!------------------------------------------------------------------------------

! flake_paramoptic_ref:
! ---------------------

!  This module contains "reference" values of the optical characteristics
!  of the lake water, lake ice and snow. These reference values may be used
!  if no information about the optical characteristics of the lake in question
!  is available. An exponential decay law for the solar radiation flux is assumed.
!  In the simplest one-band approximation,
!  the extinction coefficient for water is set to a large value,
!  leading to the absorption of 95% of the incoming radiation
!  within the uppermost 1 m of the lake water.
!  The extinction coefficients for ice and snow are taken from
!  Launiainen and Cheng (1998). The estimates for the ice correspond
!  to the uppermost 0.1 m of the ice layer and to the clear sky conditions
!  (see Table 2 in op. cit.).
!  Very large values of the extinction coefficients for ice and snow ("opaque")
!  can be used to prevent penetration of the solar radiation
!  through the snow-ice cover.

#ifdef __COSMO__
INTEGER (KIND = iintegers), PRIVATE :: & ! Help variable(s)
  i                                      ! DO loop index
#endif

#ifdef __ICON__
INTEGER (KIND = iintegers) :: i
#endif
!  Optical characteristics for water, ice and snow.
!  The simplest one-band approximation is used as a reference.
TYPE (t_opticpar_medium), PARAMETER ::                           &
  ! Water (reference)
  opticpar_water_ref = t_opticpar_medium(1,                      &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/3._ireals, (1.E+10_ireals,i=2,nband_optic_max)/))      , &
  ! Transparent Water (two-band)
  opticpar_water_trans = t_opticpar_medium(2,                             &
    (/0.10_ireals, 0.90_ireals, (0._ireals,i=3,nband_optic_max)/),      &
    (/2.0_ireals, 0.20_ireals, (1.E+10_ireals,i=3,nband_optic_max)/)) , &
  ! Transparent Water (one-band)
!_nu  opticpar_water_trans = opticpar_medium(1,                    &
!_nu    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
!_nu    (/0.30_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  ! White ice
  opticpar_whiteice_ref = t_opticpar_medium(1,                   &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/17.1_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  ! Blue ice
  opticpar_blueice_ref = t_opticpar_medium(1,                    &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/8.4_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))     , &
  ! Dry snow
  opticpar_drysnow_ref = t_opticpar_medium(1,                    &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/25.0_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  ! Melting snow
  opticpar_meltingsnow_ref = t_opticpar_medium(1,                &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/15.0_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  ! Opaque ice
  opticpar_ice_opaque = t_opticpar_medium(1,                     &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/1.0E+07_ireals, (1.E+10_ireals,i=2,nband_optic_max)/)) , &
  ! Opaque snow
  opticpar_snow_opaque = t_opticpar_medium(1,                    &
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/1.0E+07_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))

!==============================================================================

END MODULE data_flake
