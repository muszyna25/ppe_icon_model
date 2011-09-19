!+ Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------

MODULE mo_data_turbdiff

!------------------------------------------------------------------------------
!
! Description:
!  This module matches all parameters needed for module src_turbdiff:
!
! Current Code Owner: DWD, Matthias.Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  Matthias.Raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
!              2010/12/30 Matthias Raschendorfer
! Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
USE mo_kind, ONLY : &
!
    ireals   =>wp, & ! KIND-type parameter for real variables
    iintegers=>i4    ! KIND-type parameter for standard integer variables

USE mo_physical_constants, ONLY : &
!
! Physical constants and related variables:
! -------------------------------------------
!
    r_d      => rd,       & ! gas constant for dry air
    rvd_m_o  => vtmpc1,   & ! r_v/r_d - 1
    cp_d     => cpd,      & ! specific heat for dry air
    lh_v     => alv,      & ! evaporation heat
    rdocp    => rd_o_cpd, & ! r_d / cp_d
    lhocp    => alvdcp,   & ! lh_v / cp_d
    t0_melt  => tmelt,    & ! absolute zero for temperature (K)
    b3       => tmelt,    & !          -- " --

    rdv, o_m_rdv, rcpv, rcpl, con_m, con_h, &
    grav, p0ref, &
    uc1, ucl


USE mo_convect_tables, ONLY : &
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    b1       => c1es,     & ! variables for computing the saturation steam pressure
    b2w      => c3les,    & ! over water (w) and ice (e)
    b2i      => c3ies,    & !               -- " --
    b4w      => c4les,    & !               -- " --
    b4i      => c4ies,    & !               -- " --
    b234w    => c5les       ! b2w * (b3 - b4w)

USE mo_math_constants, ONLY : &
!
    uc2      => sqrt3       ! used in cloud cover diagnostics based on rel. humid.

! Parameters for tuning the turbulence parameterizations
! -----------------------------------------------------

USE data_turbulence, ONLY : rlam_mom, & ! scaling factor of the laminar boudary layer for momentum
     &                      rlam_heat,& ! scaling factor of the laminar boudary layer for heat
     &                      rat_can,  & ! factor for the canopy height
     &                      rat_sea,  & ! ratio of laminar scaling factors for heat
                                        ! over sea and land
     &                      rat_lam,  & ! ratio of laminar scaling factors for vapour and heat
     &                      z0m_dia,  & ! roughness length of a typical synoptic station
     &                      alpha0,   & ! Charnock-parameter
     &                      alpha1,   & ! parameter scaling the molek. roughness of water waves
!
     &                      c_lnd,    & ! surface area index of the land exept the leaves
     &                      c_soil,   & ! surface area index of the (evaporative) soil
     &                      c_sea,    & ! surface area index of the waves over sea
     &                      e_surf,   & ! exponent to get the effective surface area
     &                      zt_ice,   & ! freezing temperature of sea ice
     &                      z0_ice,   & ! roughness length of sea ice
     &                      tur_len,  & ! maximal turbulent length scale [m]
     &                      pat_len,  & ! lenth scale of subscale patterns over land [m]
     &                      len_min,  & ! minimal turbulent length scale [m]
     &                      vel_min,  & ! minimal velocity scale [m/s]
     &                      akt,      & ! von Karman-constant
     &                      a_heat,   & ! factor for turbulent heat transport
     &                      a_mom,    & ! factor for turbulent momentum transport
     &                      d_heat,   & ! factor for turbulent heat dissipation
     &                      d_mom,    & ! factor for turbulent momentum dissipation
     &                      c_diff,   & ! factor for turbulent diffusion of TKE
     &                      a_hshr,   & ! factor for horizontal shear production of TKE
     &                      a_stab,   & ! factor for stability correction of 
                                        !turbulent length scale
     &                      clc_diag, & ! cloud cover at saturation in statistical 
                                        !cloud diagnostic
     &                      q_crit,   & ! critical value for normalized over-saturation
     &                      c_scld,   & ! factor for liquid water flux density in
                                        ! sub grid scale clouds
     &                      tkhmin,   & ! minimal diffusion coefficients for heat
     &                      tkmmin,   & ! minimal diffusion coefficients for momentum
     &                      epsi,     & ! relative limit of accuracy for 
                                        ! comparison of numbers
     &                      tkesmot,  & ! time smoothing factor for TKE and 
                                        ! diffusion coefficients
     &                      wichfakt, & ! vertical smoothing factor for 
                                        ! explicit diffusion tendencies
     &                      securi,   & ! security factor for maximal diffusion coefficients
     &                      it_end      ! number of initialization iterations (>=0)


! Flake parameters:
! -----------------------------------------------------

USE data_flake, ONLY : &
!
  h_Ice_min_flk

! Switches controlling other physical parameterizations:
! -----------------------------------------------------



!==============================================================================

IMPLICIT NONE

#ifdef __COSMO__
PUBLIC          ! All constants and variables in this module are public
!==============================================================================

INTEGER (KIND=iintegers) :: &
!
! Switches controlling turbulent diffusion:
! ------------------------------------------
!
    itype_tran   =2,       & ! type of surface-atmosphere transfer
    imode_tran   =1,       & ! mode of surface-atmosphere transfer
    icldm_tran   =0,       & ! mode of cloud representation in transfer parametr.
!
    imode_turb   =3,       & ! mode of turbulent diffusion parametrization
    icldm_turb   =2,       & ! mode of cloud representation in turbulence parametr.
    itype_sher   =1          ! type of shear production for TKE

LOGICAL :: &
!
    ltkesso      =.FALSE., & ! calculation SSO-wake turbulence production for TKE
    ltkecon      =.FALSE., & ! consider convective buoyancy production for TKE
    lexpcor      =.FALSE., & ! explicit corrections of the implicit calculated turbul. diff.
    ltmpcor      =.FALSE., & ! consideration of thermal TKE-sources in the enthalpy budget
    lprfcor      =.FALSE., & ! using the profile values of the lowest main level instead of
                             ! the mean value of the lowest layer for surface flux calulations
    lnonloc      =.FALSE., & ! nonlocal calculation of vertical gradients used for turbul. diff.
    lcpfluc      =.FALSE., & ! consideration of fluctuations of the heat capacity of air
    limpltkediff =.TRUE.     ! use semi-implicit TKE diffusion

INTEGER (KIND=iintegers) :: &
!
! Switches controlling other physical parameterizations:
!
    itype_wcld   =2,       & ! type of water cloud diagnosis
    itype_synd   =2          ! type of diagnostics of synoptical near surface variables


LOGICAL :: &
!
    lseaice      =.FALSE., & ! sea ice parameterization active
    llake        =.FALSE. !, & ! lake model active
!
! Attention:
! The given initializations are default settings of the boundary layer
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs!
!

!-----------------------------------------------------------------------------
!CONTAINS
!-----------------------------------------------------------------------------
#endif
!==============================================================================
!==============================================================================

END MODULE mo_data_turbdiff
