!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_atm_phy_nwp_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: max_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_phy_nwp_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for atm dynamics
  !!--------------------------------------------------------------------------
  TYPE :: t_atm_phy_nwp_config

    ! namelist variables

  INTEGER ::  inwp_gscp        !> microphysics
  INTEGER ::  inwp_satad       !! saturation adjustment
  INTEGER ::  inwp_convection  !! convection
  INTEGER ::  inwp_radiation   !! radiation
  INTEGER ::  inwp_sso         !! sso
  INTEGER ::  inwp_gwd         !! non-orographic gravity wave drag
  INTEGER ::  inwp_cldcover    !! cloud cover
  INTEGER ::  inwp_turb        !! turbulence
  INTEGER ::  inwp_surface     !! surface including soil, ocean, ice,lake

  REAL(wp) :: dt_conv   !> field element for convection
  REAL(wp) :: dt_ccov   !! field element for subscale cloud cover
  REAL(wp) :: dt_rad    !! "-"                     radiation
  REAL(wp) :: dt_radheat!! "-" rad. heating from radiative fluxes with updated cosmu0 
  REAL(wp) :: dt_sso    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gwd    !! "-"  for subscale gravity waves
  REAL(wp) :: dt_gscp   !! field element for microphysics
  REAL(wp) :: dt_turb   !! field element for turbulence
  REAL(wp) :: dt_sfc    !! field element for surface
  REAL(wp) :: dt_satad  !! field element for sat. adjustment
  REAL(wp) :: dt_update !! field element for tracer phys update

  LOGICAL ::  lseaice     !> forecast with sea ice model
  LOGICAL ::  llake       !! forecst with lake model FLake
  LOGICAL ::  l3dturb     !! 3D-turbulence (additional horizontal diffusion)
  LOGICAL ::  lprog_tke   !! prognostic treatment of TKE (for itype_turb= 5/7)
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model


  END TYPE t_atm_phy_nwp_config

  !>
  !!
  TYPE(t_atm_phy_nwp_config) :: atm_phy_nwp_config(max_dom) !< shape: (n_dom)

END MODULE mo_atm_phy_nwp_config
