!>
!! Data types and variables used by the ECHAM6 physics package implemented
!! in the ICOHAM model.
!!
!! This module contains
!! <ol>
!! <li> definition of data types for organising the physical quantities in the
!!      ECHAM physics package,
!! <li> the actual variables that are declared of these types, and
!! <li> subroutines for (de-)allocating memory for the variables.
!! </ol>
!! This module is similar to mo_memory_g3b in ECHAM.
!!
!! @author Kristina Froehlich (DWD)
!! @author Marco Giorgetta (MPI-M)
!! @author Hui Wan (MPI-M)
!!
!! @par Revision History
!!  First version by Hui Wan, Marco Giorgetta and Kristina Froehlich, 2010-10-28.
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
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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
MODULE mo_echam_phy_memory

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_exception,           ONLY: message, finish
  USE mo_run_nml,             ONLY: nlev, nlevp1, nproma, ntracer
  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, igbm
  USE mo_echam_phy_nml,       ONLY: lvdiff
  USE mo_model_domain,        ONLY: t_patch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_field, prm_tend                                 !< variables
  PUBLIC :: construct_echam_phy_state, destruct_echam_phy_state !< subroutines
  PUBLIC :: t_echam_phy_field, t_echam_phy_tend

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_memory'

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------
  !>
  !! Derived data type: t_echam_phy_field
  !!
  !! This structure contains two kinds of components:
  !! <ol>
  !! <li> quantities involved in the parameterisation scheme but not in the
  !!      dynamical core, e.g., cloud cover, 10-m wind speed;
  !! <li> atmospheric state variables shared by dynamics and physics, e.g.
  !!      wind, temperature, pressure and tracer concentrations.
  !!      At each time step, the dynamical core provides initial values of
  !!      these quantites on the dynamics grid, which are then interpolated
  !!      to the physics grid and passed on to the physics package.
  !!      In the physics package, these variables may be updated once or
  !!      even more times, depending on the actual numerical schemes chosen
  !!      for the physics-dynamics coupling and the coupling between
  !!      individual parameterisation schemes;
  !! </ol>
  !!
  !! All components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma,           nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy, ntracers)
  !! </ol>
  !! Currently the physics grid has the same spatial resolution as the
  !! dynamics grid, but is unstaggered. This means
  !!
  !!    nlev_phy = nlev
  !!   nblks_phy = patch%nblks_int_c
  !!
  !! In the long run, the physics grid and dynamics grid may differ in
  !! horizontal and/or vertical resolution, or even in shape.

  TYPE t_echam_phy_field

    ! Meteorology and tracers
    REAL(wp),ALLOCATABLE :: &
      & u         (:,:,:),  &!< [m/s]   zonal wind
      & v         (:,:,:),  &!< [m/s]   meridional wind
      & vor       (:,:,:),  &!< [1/s]   relative vorticity
      & temp      (:,:,:),  &!< [K]     temperature          (tm1  of memory_g1a in ECHAM)
      & tv        (:,:,:),  &!< [K]     virtual temperature  (tvm1 of memory_g1a in ECHAM)
      & q         (:,:,:,:),&!< [kg/kg] tracer concentration (qm1, xlm1, xim1 of memory_g1a in ECHAM)
      & qx        (:,:,:),  &!< [kg/kg] total concentration of hydrometeors
      & omega     (:,:,:),  &!< [Pa/s]  vertical velocity in pressure coord. ("vervel" in ECHAM)
      & geoi      (:,:,:),  &!< [m2/s2] geopotential at half levels (vertical interfaces)
      & geom      (:,:,:),  &!< [m2/s2] geopotential at full levels (layer ave. or mid-point value)
      & presi_old (:,:,:),  &!< [Pa]    pressure at half levels at time step "old"
      & presm_old (:,:,:),  &!< [Pa]    pressure at full levels at time step "old"
      & presi_new (:,:,:),  &!< [Pa]    pressure at half levels at time step "new"
      & presm_new (:,:,:)    !< [Pa]    pressure at full levels at time step "new"


    ! Radiation
    REAL(wp),ALLOCATABLE :: &
      !
      ! insolation at TOA
      & cosmu0      (:,  :),  &!< [ ]    cos of zenith angle mu0 for radiative heating  calculation
      & flxdwswtoa  (:,  :),  &!< [W/m2] downward shortwave flux at TOA
      !
      ! shortwave surface albedo
      & albvisdir   (:,  :),  &!< [ ]    surface albedo for visible range, direct
      & albvisdif   (:,  :),  &!< [ ]    surface albedo for visible range, diffuse
      & albnirdir   (:,  :),  &!< [ ]    surface albedo for near IR range, direct
      & albnirdif   (:,  :),  &!< [ ]    surface albedo for near IR range, diffuse
      !
      ! shortwave net surface fluxes
      & vissfc      (:,  :),  &!< [ ]    solar transmissivity in VIS, net downward
      & visdffsfc   (:,  :),  &!< [ ]    diffuse fraction in VIS net downw. flux (?)
      & nirsfc      (:,  :),  &!< [ ]    solar transmissivity in NIR, net downward
      & nirdffsfc   (:,  :),  &!< [ ]    diffuse fraction in NIR net downw. flux (?)
      & parsfc      (:,  :),  &!< [ ]    solar transmissivity in PAR, downward
      & pardffsfc   (:,  :),  &!< [ ]    diffuse fraction in PAR net downw. flux (?)
      !
      ! shortwave net transmissivity
      & trsolclr    (:,:,:),  &!< [ ]    solar transmissivity  , clear sky, net downward
      & trsolall    (:,:,:),  &!< [ ]    solar transmissivity  , all   sky, net downward
      !
      ! longwave net fluxes
      & emterclr    (:,:,:),  &!< [W/m2] terrestrial emissivity, clear sky, net downward
      & emterall    (:,:,:)    !< [W/m2] terrestrial emissivity, all   sky, net downward


    ! Cloud and precipitation
    REAL(wp),ALLOCATABLE :: &
      & aclc      (:,:,:),  &!< [m2/m2] fractional cloud cover   (was aclc in echam)
      & aclcac    (:,:,:),  &!< accumulated cloud cover
      & aclcov    (:,  :),  &!< (vertically integrated) total cloud cover
      & acdnc     (:,:,:),  &!< cloud droplet number concentration [1/m**3]
      & xvar      (:,:,:),  &!< variance of total water amount qv+qi+ql [kg/kg] (memory_g3b)
      & xskew     (:,:,:),  &!< skewness of total water amount qv+qi+ql [kg/kg]
      & relhum    (:,:,:),  &!< relative humidity (relhum of memory_g3b in ECHAM)
      & aprl      (:,  :),  &!< sfc precip amount, rain+snow, large scale   [kg/m**2]
      & aprc      (:,  :),  &!< sfc precip amount, rain+snow, convective    [kg/m**2]
      & aprs      (:,  :),  &!< sfc snow   amount, large scale + convective [kg/m**2]
      & rsfl      (:,  :),  &!< sfc rain flux, large scale [kg/m**2/s]
      & rsfc      (:,  :),  &!< sfc rain flux, convective  [kg/m**2/s]
      & ssfl      (:,  :),  &!< sfc snow flux, large scale [kg/m**2/s]
      & ssfc      (:,  :),  &!< sfc snow flux, convective  [kg/m**2/s]
      & qvi       (:,  :),  &!< vertically integrated water vapor [kg/m**2s]
      & xlvi      (:,  :),  &!< vertically integrated cloud water [kg/m**2s]
      & xivi      (:,  :)    !< vertically integrated cloud ice [kg/m**2s]

    REAL(wp),ALLOCATABLE :: &
      & rintop (:,  :),     &!< low lever inversion, computed by "cover" (memory_g3b)
      & rtype  (:,  :),     &!< type of convection 0...3. (in memory_g3b in ECHAM)
      & topmax (:,  :),     &!< maximum height of convective cloud tops [Pa] (memory_g3b)
      & thvsig (:,  :)       !< Std. dev. of virtual potential temperature at the upper
                             !< interface of the lowest model layer.
                             !< Computed in "vdiff" by getting the square root of
                             !< thvvar(:,nlev-1,:). Used by "cucall".

    ! Turbulence

    REAL(wp),ALLOCATABLE :: &
      & tke       (:,:,:),  &!< turbulent kinetik energy at step n+1
      & tkem0     (:,:,:),  &!< turbulent kinetik energy at step n
      & tkem1     (:,:,:)    !< turbulent kinetik energy at step n-1

    ! need only for vdiff ++++
    REAL(wp),ALLOCATABLE :: &
      & ri        (:,:,:),  &!< moist Richardson number at layer interfaces
      & mixlen    (:,:,:),  &!< mixing length at layer interfaces
      & thvvar    (:,:,:)    !< variance of virtual potential temperature at layer interfaces.
                             !< Computed in "vdiff" by solving a prognostic equation of
                             !< the variance. Used for getting "thvsig".

    REAL(wp),ALLOCATABLE :: &
      & cfm    (:,:,:),     &!< turbulent exchange coefficient
      & cfm_sfc(:,:,:),     &!< turbulent exchange coefficient
      & cfh    (:,:,:),     &!< turbulent exchange coefficient
      & cfh_sfc(:,:,:),     &!< turbulent exchange coefficient
      & cfv    (:,:,:),     &!< turbulent exchange coefficient
      & cftke  (:,:,:),     &!< turbulent exchange coefficient
      & cfthv  (:,:,:)       !< turbulent exchange coefficient


    REAL(wp),ALLOCATABLE :: &
      & coriol(:,:),        &!< Coriolis parameter, needed for diagnosing PBL height.
      & ghpbl (:,:),        &!< geopotential of the top of the atmospheric boundary layer
      & z0m   (:,:,:),      &!< aerodynamic roughness length (grid-box mean and over each surface type)
      & ustar (:,:),        &!<
      & kedisp(:,:),        &!< time-mean (or integrated?) vertically integrated dissipation of kinetic energy
      & ocu   (:,:),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)          !< northward velocity of ocean surface current

    ! need only for vdiff ----

    ! Surface parameters

    LOGICAL, ALLOCATABLE :: &
      & lfland(:,:),        &!< .TRUE. when fraction of land > 0.
      & lfglac(:,:)          !< .TRUE. when fraction of glaciated land > 0.

    REAL(wp),ALLOCATABLE :: &
      & lsmask(:,:),        &!< land-sea mask. (1. = land, 0. = sea/lakes) (slm in memory_g3b)
      & glac  (:,:),        &!< fraction of land covered by glaciers (glac in memory_g3b)
      & seaice(:,:),        &!< ice cover given as the fraction of (1- slm) (seaice in memory_g3b)
      & icefrc(:,:),        &!< ice cover given as the fraction of grid box (friac  in memory_g3b)
      & tsfc  (:,:,:),      &!< surface temperature over land/water/ice (tsw/l/i in memory_g3b)
      & qs_sfc(:,:,:)        !< saturation specitifc humidity at surface

!!$    ! Variables for debugging
!!$    REAL(wp),ALLOCATABLE :: &
!!$      ! 2d arrays
!!$      & debug_2d_1(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_2(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_3(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_4(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_5(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_6(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_7(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_8(:,:),    &!< 2d variable for debug purposes
!!$      ! 3d arrays
!!$      & debug_3d_1(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_2(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_3(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_4(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_5(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_6(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_7(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_8(:,:,:)    !< 3d variable for debug purposes

  END TYPE t_echam_phy_field

  !>
  !! Data type containing the tendencies returned by the individual parameterizations
  !!
  !! The components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma, nlev_phy, nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy, ntracers)
  !! </ol>
  !!
  TYPE t_echam_phy_tend

    REAL(wp), ALLOCATABLE ::   &
      !
      ! all processes
      !
      &    u        (:,:,:)  , & !< accumulated tendency
      &    v        (:,:,:)  , & !< accumulated tendency
      & temp        (:,:,:)  , & !< accumulated tendency
      &    q        (:,:,:,:), & !< accumulated tendency
      !
      ! cloud microphysics
      !
      & temp_cld    (:,:,:)  , & !< T tendency from cloud microphysical processes
      &    q_cld    (:,:,:,:), & !< tracer tendency from cloud microphysical process
      !
      ! cumulus convection
      !
      & temp_cnv    (:,:,:),   & !< temperature tendency from cumulus convection
      &    u_cnv    (:,:,:),   & !< Zonal      wind tendency from cumulus convection
      &    v_cnv    (:,:,:),   & !< Meridional wind tendency from cumulus convection
      &    q_cnv    (:,:,:,:), & !< tracer tendency due to cumulus convection
      &    x_dtr    (:,:,:),   & !< cloud water/ice tendency due to detrainment (memory_g3b:xtec)
      !
      ! vertical turbulent mixing ("vdiff")
      !
      & temp_vdf  (:,:,:),     & !< temperature tendency due to turbulent mixing
      &    u_vdf  (:,:,:),     & !< u-wind tendency due to turbulent mixing
      &    v_vdf  (:,:,:),     & !< v-wind tendency due to turbulent mixing
      &    q_vdf  (:,:,:,:),   & !< tracer tendency due to turbulent mixing
      !
      ! atmospheric gravity wave drag
      !
!!$      & u_gwd       (:,:,:)  , & !< ZonalW-tendency from gravity wave drag
!!$      & v_gwd       (:,:,:)  , & !< MeridW-tendency from gravity wave drag
!!$      & temp_gwd    (:,:,:)  , & !< Temp-tendency from gravity wave drag
      !
      ! subgrid scale orographic (sso) blocking and gravity wave drag
      !
!!$      & u_sso       (:,:,:)  , & !< ZonalW-tendency from sso drag
!!$      & v_sso       (:,:,:)  , & !< MeridW-tendency from sso drag
!!$      & temp_sso    (:,:,:)  , & !< Temp-tendency from sso drag
      !
      ! radiation
      !
      & temp_radsw  (:,:,:)  , & !< Temp-tendency from radiation
      & temp_radlw  (:,:,:)      !< Temp-tendency from radiation

  END TYPE t_echam_phy_tend

  !!--------------------------------------------------------------------------
  !!                            VARIABLES
  !!--------------------------------------------------------------------------
  !! The variable names have the prefix "prm_" in order to emphasize that they
  !! are defined for and used in parameterisations.

  TYPE(t_echam_phy_field),TARGET,ALLOCATABLE :: prm_field(:)  !< shape: (n_dom)
  TYPE(t_echam_phy_tend ),TARGET,ALLOCATABLE :: prm_tend (:)  !< shape: (n_dom)

CONTAINS
  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR (DE-)ALLOCATING MEMORY
  !!--------------------------------------------------------------------------
  !>
  !! Allocate memory for the physics state
  !!
  SUBROUTINE construct_echam_phy_state( p_patch )

    TYPE(t_patch),INTENT(IN) :: p_patch(:)
    INTEGER :: ndomain, jg, ist
    !---

    CALL message(TRIM(thismodule),'Construction of ECHAM physics state started.')

    ! Allocate prm_field and prm_tend as arrays

    ndomain = SIZE(p_patch)
    ALLOCATE( prm_field(ndomain), prm_tend(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of prm_field/tend array failed')

    ! Allocate memory for the omponents of prm_field and prm_tend on all grid levels

    DO jg = 1,ndomain
      CALL construct_echam_phy_field( p_patch(jg), prm_field(jg) )
      CALL construct_echam_phy_tend ( p_patch(jg), prm_tend (jg) )
    ENDDO
    CALL message(TRIM(thismodule),'Construction of ECHAM physics state finished.')

  END SUBROUTINE construct_echam_phy_state
  !-------------
  !>
  !! Release memory used by the physics state
  !!
  SUBROUTINE destruct_echam_phy_state

    INTEGER :: ndomain  !< total # of grid levels/domains
    INTEGER :: jg       !< grid level/domain index
    INTEGER :: ist      !< system status code
    !---
    CALL message(TRIM(thismodule),'Destruction of ECHAM physics state started.')

    ndomain = SIZE(prm_field)

    DO jg = 1,ndomain
      CALL destruct_echam_phy_field( prm_field(jg) )
      CALL destruct_echam_phy_tend ( prm_tend (jg) )
    ENDDO

    DEALLOCATE( prm_field, prm_tend, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'deallocation of prm_field/tend array failed')

    CALL message(TRIM(thismodule),'Destruction of ECHAM physics state finished.')

  END SUBROUTINE destruct_echam_phy_state
  !-------------
  !>
  !!
  !!
  SUBROUTINE construct_echam_phy_field( p_patch, field )

    TYPE(t_patch),INTENT(IN)                :: p_patch
    TYPE(t_echam_phy_field),INTENT(INOUT) :: field
    INTEGER :: nblks, ist

    nblks = p_patch%nblks_c

    ! Meteorology
    ALLOCATE( field% u         (nproma,nlev  ,nblks),          &
      &       field% v         (nproma,nlev  ,nblks),          &
      &       field% vor       (nproma,nlev  ,nblks),          &
      &       field% temp      (nproma,nlev  ,nblks),          &
      &       field% tv        (nproma,nlev  ,nblks),          &
      &       field% q         (nproma,nlev  ,nblks,ntracer),  &
      &       field% qx        (nproma,nlev  ,nblks),          &
      &       field% omega     (nproma,nlev  ,nblks),          &
      &       field% geoi      (nproma,nlevp1,nblks),          &
      &       field% geom      (nproma,nlev  ,nblks),          &
      &       field% presi_old (nproma,nlevp1,nblks),          &
      &       field% presm_old (nproma,nlev  ,nblks),          &
      &       field% presi_new (nproma,nlevp1,nblks),          &
      &       field% presm_new (nproma,nlev  ,nblks),          &
      &       STAT=ist                                         )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of ECHAM meteorol. fields failed')

    ! Radiation
    ALLOCATE( field% cosmu0    (nproma,       nblks),          &
      &       field% flxdwswtoa(nproma,       nblks),          &
      !
      &       field% albvisdir (nproma,       nblks),          &
      &       field% albvisdif (nproma,       nblks),          &
      &       field% albnirdir (nproma,       nblks),          &
      &       field% albnirdif (nproma,       nblks),          &
      !
      &       field% vissfc    (nproma,       nblks),          &
      &       field% visdffsfc (nproma,       nblks),          &
      &       field% nirsfc    (nproma,       nblks),          &
      &       field% nirdffsfc (nproma,       nblks),          &
      &       field% parsfc    (nproma,       nblks),          &
      &       field% pardffsfc (nproma,       nblks),          &
      !
      &       field% emterclr  (nproma,nlevp1,nblks),          &
      &       field% emterall  (nproma,nlevp1,nblks),          &
      &       field% trsolclr  (nproma,nlevp1,nblks),          &
      &       field% trsolall  (nproma,nlevp1,nblks),          &
      &       STAT=ist                                         )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of radiation field state failed')

    ! Cloud and precipitation
    ALLOCATE( field% aclc   (nproma,nlev  ,nblks), &
      &       field% aclcac (nproma,nlev  ,nblks), &
      &       field% aclcov (nproma,       nblks), &
      &       field% acdnc  (nproma,nlev  ,nblks), &
      &       field% xvar   (nproma,nlev  ,nblks), &
      &       field% xskew  (nproma,nlev  ,nblks), &
      &       field% relhum (nproma,nlev  ,nblks), &
      &       field% aprl   (nproma,       nblks), &
      &       field% aprc   (nproma,       nblks), &
      &       field% aprs   (nproma,       nblks), &
      &       field% rsfl   (nproma,       nblks), &
      &       field% rsfc   (nproma,       nblks), &
      &       field% ssfl   (nproma,       nblks), &
      &       field% ssfc   (nproma,       nblks), &
      &       field%  qvi   (nproma,       nblks), &
      &       field% xlvi   (nproma,       nblks), &
      &       field% xivi   (nproma,       nblks), &
      &       field% rintop (nproma,       nblks), &
      &       field% rtype  (nproma,       nblks), &
      &       field% topmax (nproma,       nblks), &
      &       field% tke    (nproma,nlev  ,nblks), &
      &       field% thvsig (nproma,       nblks), &
      &       STAT=ist                            )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of cloud fields failed')

    ! Turbulence
    IF (lvdiff) THEN

      ALLOCATE( field% ri     (nproma,nlev,nblks), &
        &       field% mixlen (nproma,nlev,nblks), &
        &       field% thvvar (nproma,nlev,nblks), &
        &       field% tkem0  (nproma,nlev,nblks), &
        &       field% tkem1  (nproma,nlev,nblks), &
        &       STAT=ist                           )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),   &
        & 'allocation of turbulence fields failed')

      ALLOCATE( field% cfm    (nproma,nlev,     nblks), &
        &       field% cfm_sfc(nproma,nsfc_type,nblks), &
        &       field% cfh    (nproma,nlev,     nblks), &
        &       field% cfh_sfc(nproma,nsfc_type,nblks), &
        &       field% cfv    (nproma,nlev,     nblks), &
        &       field% cftke  (nproma,nlev,     nblks), &
        &       field% cfthv  (nproma,nlev,     nblks)  )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),   &
        & 'allocation of turbulent exchange coefficients failed')

      ALLOCATE( field% coriol (nproma,nblks),                &
        &       field% ghpbl  (nproma,nblks),                &
        &       field% z0m    (nproma,igbm:nsfc_type,nblks), &
        &       field% ustar  (nproma,nblks),                &
        &       field% kedisp (nproma,nblks),                &
        &       field% ocu    (nproma,nblks),                &
        &       field% ocv    (nproma,nblks),                &
        &       STAT=ist                                     )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),   &
        & 'allocation of failed')

    ENDIF ! lvdiff

    ! Surface
    ALLOCATE( field% lfland (nproma, nblks),                 &
      &       field% lfglac (nproma, nblks),                 &
      &       field% lsmask (nproma, nblks),                 &
      &       field% glac   (nproma, nblks),                 &
      &       field% seaice (nproma, nblks),                 &
      &       field% icefrc (nproma, nblks),                 &
      &       field% tsfc   (nproma, igbm:nsfc_type, nblks), &
      &       field% qs_sfc (nproma,      nsfc_type, nblks), &
      &       STAT=ist                                       )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),  &
      & 'allocation of surface fields failed')

!!$    ! Variables for debugging
!!$    ALLOCATE( field% debug_2d_1(nproma,         nblks), &
!!$      &       field% debug_2d_2(nproma,         nblks), &
!!$      &       field% debug_2d_3(nproma,         nblks), &
!!$      &       field% debug_2d_4(nproma,         nblks), &
!!$      &       field% debug_2d_5(nproma,         nblks), &
!!$      &       field% debug_2d_6(nproma,         nblks), &
!!$      &       field% debug_2d_7(nproma,         nblks), &
!!$      &       field% debug_2d_8(nproma,         nblks), &
!!$      &       field% debug_3d_1(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_2(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_3(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_4(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_5(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_6(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_7(nproma, nlev,   nblks), &
!!$      &       field% debug_3d_8(nproma, nlev,   nblks), &
!!$      &       STAT=ist                          )
!!$
!!$    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of debugging fields failed')

  END SUBROUTINE construct_echam_phy_field
  !-------------
  !>
  !!
  !!
  SUBROUTINE destruct_echam_phy_field( field )

    TYPE(t_echam_phy_field),INTENT(INOUT) :: field
    INTEGER :: ist

    DEALLOCATE( field% u,         &
      &         field% v,         &
      &         field% temp,      &
      &         field% tv,        &
      &         field% q,         &
      &         field% qx,        &
      &         field% omega,     &
      &         field% geoi,      &
      &         field% geom,      &
      &         field% presi_old, &
      &         field% presm_old, &
      &         field% presi_new, &
      &         field% presm_new, &
      &         field% vor,       &
      &         STAT=ist          )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of ECHAM meteorology fields failed')

!!$    DEALLOCATE( field% debug_2d_1, &
!!$      &         field% debug_2d_2, &
!!$      &         field% debug_2d_3, &
!!$      &         field% debug_2d_4, &
!!$      &         field% debug_2d_5, &
!!$      &         field% debug_2d_6, &
!!$      &         field% debug_2d_7, &
!!$      &         field% debug_2d_8, &
!!$      &         field% debug_3d_1, &
!!$      &         field% debug_3d_2, &
!!$      &         field% debug_3d_3, &
!!$      &         field% debug_3d_4, &
!!$      &         field% debug_3d_5, &
!!$      &         field% debug_3d_6, &
!!$      &         field% debug_3d_7, &
!!$      &         field% debug_3d_8, &
!!$      &         STAT=ist          )
!!$
!!$    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
!!$      & 'deallocation of debugging fields failed')

    ! Radiation
    DEALLOCATE( field% cosmu0,    &
      &         field% flxdwswtoa,&
      !
      &         field% albvisdir, &
      &         field% albvisdif, &
      &         field% albnirdir, &
      &         field% albnirdif, &
      !
      &         field% vissfc,    &
      &         field% visdffsfc, &
      &         field% nirsfc,    &
      &         field% nirdffsfc, &
      &         field% parsfc,    &
      &         field% pardffsfc, &
      !
      &         field% emterclr,  &
      &         field% emterall,  &
      &         field% trsolclr,  &
      &         field% trsolall,  &
      &         STAT=ist          )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of radiation field state failed')

    ! Cloud and precipitation
    DEALLOCATE( field% aclc,      &
      &         field% aclcac,    &
      &         field% aclcov,    &
      &         field% acdnc,     &
      &         field% xvar,      &
      &         field% xskew,     &
      &         field% relhum,    &
      &         field% aprl,      &
      &         field% aprc,      &
      &         field% aprs,      &
      &         field% rsfl,      &
      &         field% rsfc,      &
      &         field% ssfl,      &
      &         field% ssfc,      &
      &         field%  qvi,      &
      &         field% xlvi,      &
      &         field% xivi,      &
      &         field% rintop,    &
      &         field% rtype,     &
      &         field% topmax,    &
      &         field% tke,       &
      &         field% thvsig,    &
      &         STAT=ist          )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of cloud fields failed')

    ! Turbulence
    IF (lvdiff) THEN

      DEALLOCATE( field% ri,        &
        &         field% mixlen,    &
        &         field% thvvar,    &
        &         field% tkem0,     &
        &         field% tkem1,     &
        &         STAT=ist          )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
        & 'deallocation of turbulence fields failed')

      DEALLOCATE( field% cfm    , &
        &         field% cfm_sfc, &
        &         field% cfh    , &
        &         field% cfh_sfc, &
        &         field% cfv    , &
        &         field% cftke  , &
        &         field% cfthv    )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
        & 'deallocation of turbulent exchange coefficients failed')

      DEALLOCATE( field% coriol,    &
        &         field% ghpbl,     &
        &         field% z0m,       &
        &         field% ustar,     &
        &         field% kedisp,    &
        &         field% ocu,       &
        &         field% ocv,       &
        &         STAT=ist          )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
        & 'deallocation of surface fields failed')

    ENDIF ! lvdiff

    ! Surface
    DEALLOCATE( field% lfland,    &
      &         field% lfglac,    &
      &         field% lsmask,    &
      &         field% glac,      &
      &         field% seaice,    &
      &         field% icefrc,    &
      &         field% tsfc,      &
      &         field% qs_sfc     )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of surface fields failed')



  END SUBROUTINE destruct_echam_phy_field
  !-------------
  !>
  !!
  !!
  SUBROUTINE construct_echam_phy_tend( p_patch, tend )

    TYPE(t_patch),INTENT(IN)               :: p_patch
    TYPE(t_echam_phy_tend),INTENT(INOUT) :: tend
    INTEGER :: nblks, ist

    nblks = p_patch%nblks_c

    ALLOCATE( tend%    u      (nproma,nlev,nblks),          &
      &       tend%    v      (nproma,nlev,nblks),          &
      &       tend% temp      (nproma,nlev,nblks),          &
      &       tend%    q      (nproma,nlev,nblks,ntracer),  &
      &       tend% temp_cld  (nproma,nlev,nblks),          &
      &       tend%    q_cld  (nproma,nlev,nblks,ntracer),  &
      &       tend% temp_cnv  (nproma,nlev,nblks),          &
      &       tend%    q_cnv  (nproma,nlev,nblks,ntracer),  &
      &       tend%    u_cnv  (nproma,nlev,nblks),          &
      &       tend%    v_cnv  (nproma,nlev,nblks),          &
      &       tend%    x_dtr  (nproma,nlev,nblks),          &
      &       tend% temp_vdf  (nproma,nlev,nblks),          &
      &       tend%    u_vdf  (nproma,nlev,nblks),          &
      &       tend%    v_vdf  (nproma,nlev,nblks),          &
      &       tend%    q_vdf  (nproma,nlev,nblks,ntracer),  &
      &       tend% temp_radsw(nproma,nlev,nblks),          &
      &       tend% temp_radlw(nproma,nlev,nblks),          &
      &       STAT=ist                                      )

    IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'allocation of ECHAM tendency state failed')

  END SUBROUTINE construct_echam_phy_tend
  !-------------
  !>
  !!
  !!
  SUBROUTINE destruct_echam_phy_tend( tend )

    TYPE(t_echam_phy_tend),INTENT(INOUT) :: tend
    INTEGER :: ist

    DEALLOCATE( tend%    u,      &
      &         tend%    v,      &
      &         tend% temp,      &
      &         tend%    q,      &
      &         tend% temp_cld,  &
      &         tend%    q_cld,  &
      &         tend% temp_cnv,  &! cumulus convection
      &         tend%    q_cnv,  &
      &         tend%    u_cnv,  &
      &         tend%    v_cnv,  &
      &         tend%    x_dtr,  &
      &         tend% temp_vdf,  &! turbulent mixing
      &         tend%    q_vdf,  &
      &         tend%    u_vdf,  &
      &         tend%    v_vdf,  &
      &         tend% temp_radsw,&! radiation, shortwave
      &         tend% temp_radlw,&! radiation, longwave
      &         STAT=ist         )

    IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'deallocation of ECHAM tendency state failed')

  END SUBROUTINE destruct_echam_phy_tend
  !-------------

END MODULE mo_echam_phy_memory
