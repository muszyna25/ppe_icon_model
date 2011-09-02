#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
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
!! This module has a functionality similar to mo_memory_g3b in ECHAM,
!! but uses derived data types in order to allow for local refinement.
!!
!! @author Hui Wan (MPI-M)
!! @author Marco Giorgetta (MPI-M)
!! @author Kristina Froehlich (DWD, MPI-M)
!! @author Luis Kornblueh (MPI-M)
!!
!! @par Revision History
!! First version by Hui Wan, Marco Giorgetta and Kristina Froehlich, 2010-10-28.
!! Memory allocation method changed from explicit allocation to Luis' 
!! infrastructure by Hui Wan (MPI-M, 2011-04-24)
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
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,  ONLY: nproma
  USE mo_advection_config,    ONLY: advection_config
  USE mo_icoham_sfc_indices,  ONLY: nsfc_type
!   USE mo_echam_phy_nml,       ONLY: lvdiff
   USE mo_echam_phy_config,    ONLY: get_lvdiff
!   USE mo_echam_phy_config,    ONLY: echam_phy_config
  USE mo_model_domain,        ONLY: t_patch

  USE mo_linked_list,  ONLY: t_var_list
  USE mo_var_list,     ONLY: default_var_list_settings, &
                           & add_var, add_ref,          &
                           & new_var_list,              &
                           & delete_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants 


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_field, prm_tend                         !< variables
  PUBLIC :: prm_field_list, prm_tend_list               !< variable lists
  PUBLIC :: construct_echam_phy_state                   !< subroutine
  PUBLIC :: destruct_echam_phy_state                    !< subroutines
  PUBLIC :: t_echam_phy_field, t_echam_phy_tend         !< derived types

#ifdef HAVE_F95
  PUBLIC :: t_ptr2d, t_ptr3d
#endif

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_memory'

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------
  !>
  !! Derived data types for building pointer arrays
  !!
  TYPE t_ptr2d
    REAL(wp),POINTER :: p(:,:)    ! pointer to 2D (spatial) array
  END TYPE t_ptr2d

  TYPE t_ptr3d
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr3d

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
    REAL(wp),POINTER ::     &
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

    TYPE(t_ptr3d),ALLOCATABLE :: q_ptr(:)


    ! Radiation
    REAL(wp),POINTER ::       &
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
      & emterall    (:,:,:),  & !< [W/m2] terrestrial emissivity, all   sky, net downward
      & o3          (:,:,:)     !< temporary set ozone mass mixing ratio  

    ! Cloud and precipitation
    REAL(wp),POINTER ::     &
      & aclc      (:,:,:),  &!< [m2/m2] cloud area fractional
      & aclcac    (:,:,:),  &!< area fractional accumulated over output interval
      & aclcov    (:,  :),  &!< total cloud cover accumulated over output interval
      & acdnc     (:,:,:),  &!< cloud droplet number concentration [1/m**3]
      & xvar      (:,:,:),  &!< variance of total water amount qv+qi+ql [kg/kg] (memory_g3b)
      & xskew     (:,:,:),  &!< skewness of total water amount qv+qi+ql [kg/kg]
      & relhum    (:,:,:),  &!< relative humidity (relhum of memory_g3b in ECHAM)
      & aprl      (:,  :),  &!< (time accum) sfc precip amount, rain+snow, large-scale     [kg m-2]
      & aprc      (:,  :),  &!< (time accum) sfc precip amount, rain+snow, convective      [kg m-2]
      & aprs      (:,  :),  &!< (time accum) sfc snowfall amount, large scale + convective [kg m-2]
      & rsfl      (:,  :),  &!< sfc rain flux, large scale [kg m-2 s-1]
      & rsfc      (:,  :),  &!< sfc rain flux, convective  [kg m-2 s-1]
      & ssfl      (:,  :),  &!< sfc snow flux, large scale [kg m-2 s-1]
      & ssfc      (:,  :),  &!< sfc snow flux, convective  [kg m-2 s-1]
      & qvi       (:,  :),  &!< (time accum) vertically integrated water vapor [kg/m**2s]
      & xlvi      (:,  :),  &!< (time accum) vertically integrated cloud water [kg/m**2s]
      & xivi      (:,  :)    !< (time accum) vertically integrated cloud ice [kg/m**2s]

    REAL(wp),POINTER :: &
      & rintop (:,  :),     &!< low lever inversion, computed by "cover" (memory_g3b)
      & rtype  (:,  :),     &!< type of convection 0...3. (in memory_g3b in ECHAM)
      & topmax (:,  :),     &!< maximum height of convective cloud tops [Pa] (memory_g3b)
      & thvsig (:,  :)       !< Std. dev. of virtual potential temperature at the upper
                             !< interface of the lowest model layer.
                             !< Computed in "vdiff" by getting the square root of
                             !< thvvar(:,nlev-1,:). Used by "cucall".

    ! Turbulence

    REAL(wp),POINTER ::     &
      & tke       (:,:,:),  &!< turbulent kinetik energy at step n+1
      & tkem0     (:,:,:),  &!< turbulent kinetik energy at step n
      & tkem1     (:,:,:)    !< turbulent kinetik energy at step n-1

    ! need only for vdiff ++++
    REAL(wp),POINTER ::     &
      & ri        (:,:,:),  &!< moist Richardson number at layer interfaces
      & mixlen    (:,:,:),  &!< mixing length at layer interfaces
      & thvvar    (:,:,:)    !< variance of virtual potential temperature at layer interfaces.
                             !< Computed in "vdiff" by solving a prognostic equation of
                             !< the variance. Used for getting "thvsig".

    REAL(wp),POINTER ::      &
      & cfm     (:,:,:),     &!< turbulent exchange coefficient
      & cfm_tile(:,:,:),     &!< turbulent exchange coefficient
      & cfh     (:,:,:),     &!< turbulent exchange coefficient
      & cfh_tile(:,:,:),     &!< turbulent exchange coefficient
      & cfv     (:,:,:),     &!< turbulent exchange coefficient
      & cftke   (:,:,:),     &!< turbulent exchange coefficient
      & cfthv   (:,:,:)       !< turbulent exchange coefficient

    TYPE(t_ptr2d),ALLOCATABLE :: cfm_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: cfh_tile_ptr(:)

    REAL(wp),POINTER ::     &
      & coriol(:,:),        &!< Coriolis parameter, needed for diagnosing PBL height.
      & ghpbl (:,:),        &!< geopotential of the top of the atmospheric boundary layer
      & z0m_tile(:,:,:),    &!< aerodynamic roughness length (over each surface type)
      & z0m   (:,:),        &!< aerodynamic roughness length (grid box mean)
      & ustar (:,:),        &!<
      & kedisp(:,:),        &!< time-mean (or integrated?) vertically integrated dissipation of kinetic energy
      & ocu   (:,:),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)          !< northward velocity of ocean surface current

    TYPE(t_ptr2d),ALLOCATABLE :: z0m_tile_ptr(:)

    ! need only for vdiff ----

    ! Surface variables

    LOGICAL, POINTER :: &
      & lfland(:,:),        &!< .TRUE. when fraction of land > 0.
      & lfglac(:,:)          !< .TRUE. when fraction of glaciated land > 0.

    REAL(wp),POINTER :: &
      & lsmask(:,:),        &!< land-sea mask. (1. = land, 0. = sea/lakes) (slm in memory_g3b)
      & glac  (:,:),        &!< fraction of land covered by glaciers (glac in memory_g3b)
      & seaice(:,:),        &!< ice cover given as the fraction of (1- slm) (seaice in memory_g3b)
      & icefrc(:,:),        &!< ice cover given as the fraction of grid box (friac  in memory_g3b)
      & tsfc_tile (:,:,:),  &!< surface temperature over land/water/ice (tsw/l/i in memory_g3b)
      & tsfc      (:,  :),  &!< surface temperature, grid box mean
      & qs_sfc_tile(:,:,:)   !< saturation specitifc humidity at surface 

    TYPE(t_ptr2d),ALLOCATABLE ::   tsfc_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: qs_sfc_tile_ptr(:)

    REAL(wp),POINTER :: &
      & lhflx_ac  (:,  :), &!< (time accum) grid box mean latent   heat flux at surface 
      & shflx_ac  (:,  :), &!< (time accum) grid box mean sensible heat flux at surface 
      &  evap_ac  (:,  :), &!< (time accum) grid box mean evaporation at surface 
      & lhflx_tile(:,:,:), &!< (instantaneous) latent   heat flux at surface 
      & shflx_tile(:,:,:), &!< (instantaneous) sensible heat flux at surface 
      &  evap_tile(:,:,:)   !< (instantaneous) evaporation at surface 

    TYPE(t_ptr2d),ALLOCATABLE :: lhflx_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: shflx_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE ::  evap_tile_ptr(:)

    REAL(wp),POINTER :: &
      & u_stress_ac  (:,  :), &!< (time accum) grid box mean wind stress 
      & v_stress_ac  (:,  :), &!< (time accum) grid box mean wind stress 
      & u_stress_tile(:,:,:), &!< (instantaneous) wind stress 
      & v_stress_tile(:,:,:)   !< (instantaneous) wind stress 

    TYPE(t_ptr2d),ALLOCATABLE :: u_stress_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: v_stress_tile_ptr(:)

!!$    ! Variables for debugging
!!$    REAL(wp),POINTER :: &
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

    REAL(wp), POINTER ::   &
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
      & temp_cld    (:,:,:)  , & !< temperature tendency from cloud microphysical processes
      &    q_cld    (:,:,:,:), & !< tracer tendency from cloud microphysical process
      !
      ! cumulus convection
      !
      & temp_cnv    (:,:,:),   & !< temperature tendency from cumulus convection
      &    u_cnv    (:,:,:),   & !< u-wind tendency from cumulus convection
      &    v_cnv    (:,:,:),   & !< v-wind tendency from cumulus convection
      &    q_cnv    (:,:,:,:), & !< tracer tendency from cumulus convection
      &    x_dtr    (:,:,:),   & !< cloud water/ice tendency due to detrainment (memory_g3b:xtec)
      !
      ! vertical turbulent mixing ("vdiff")
      !
      & temp_vdf  (:,:,:),     & !< temperature tendency due to turbulent mixing
      &    u_vdf  (:,:,:),     & !< u-wind tendency due to turbulent mixing
      &    v_vdf  (:,:,:),     & !< v-wind tendency due to turbulent mixing
      &    q_vdf  (:,:,:,:),   & !< tracer tendency due to turbulent mixing
      !
      ! Hines param. for atmospheric gravity waves
      !
      & u_gwh       (:,:,:)  , & !< u-wind tendency from Hines gravity wave param.
      & v_gwh       (:,:,:)  , & !< v-wind tendency from Hines gravity wave param.
      & temp_gwh    (:,:,:)  , & !< temperature tendency from Hines gravity wave param.
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

    TYPE(t_ptr3d),ALLOCATABLE ::     q_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_cld_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_cnv_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_vdf_ptr(:)

  END TYPE t_echam_phy_tend

  !!--------------------------------------------------------------------------
  !!                          STATE VARIABLES 
  !!--------------------------------------------------------------------------
  !! The variable names have the prefix "prm_" in order to emphasize that they
  !! are defined for and used in parameterisations.

  TYPE(t_echam_phy_field),ALLOCATABLE,TARGET :: prm_field(:)  !< shape: (n_dom)
  TYPE(t_echam_phy_tend ),ALLOCATABLE,TARGET :: prm_tend (:)  !< shape: (n_dom)

  !!--------------------------------------------------------------------------
  !!                          VARIABLE LISTS
  !!--------------------------------------------------------------------------
  TYPE(t_var_list),ALLOCATABLE :: prm_field_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),ALLOCATABLE :: prm_tend_list (:)  !< shape: (n_dom)
 
CONTAINS
  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_echam_phy_state( ntracer, patch_array )

    INTEGER,INTENT(IN) :: ntracer
    TYPE(t_patch),INTENT(IN) :: patch_array(:)
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list
    INTEGER :: ndomain, jg, ist, nblks, nlev

    !---

    CALL message(TRIM(thismodule),'Construction of ECHAM physics state started.')

    ! Allocate pointer arrays prm_field and prm_tend, 
    ! as well as the corresponding list arrays.

    ndomain = SIZE(patch_array)
    ALLOCATE( prm_field(ndomain), prm_tend(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      &'allocation of prm_field/tend array failed')

    ALLOCATE( prm_field_list(ndomain), prm_tend_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      &'allocation of prm_field/tend list array failed')

    ! Build a field list and a tendency list for each grid level.
    ! This includes memory allocation. 

    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c
      nlev  = patch_array(jg)%nlev

      ! get ctracer_list
      ctracer_list = advection_config(jg)%ctracer_list

      WRITE(listname,'(a,i2.2)') 'prm_field_D',jg
      CALL new_echam_phy_field_list(nproma, nlev, nblks, ntracer, ctracer_list, &
                                   & nsfc_type, TRIM(listname), 'prm_',         &
                                   & prm_field_list(jg), prm_field(jg)          )

      WRITE(listname,'(a,i2.2)') 'prm_tend_D',jg
      CALL new_echam_phy_tend_list( nproma, nlev, nblks, ntracer, ctracer_list, &
                                  & TRIM(listname), 'prm_tend_',                &
                                  & prm_tend_list(jg), prm_tend(jg)             )
    ENDDO
    CALL message(TRIM(thismodule),'Construction of ECHAM physics state finished.')

  END SUBROUTINE construct_echam_phy_state
  !-------------
  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_echam_phy_state

    INTEGER :: ndomain  !< total # of grid levels/domains
    INTEGER :: jg       !< grid level/domain index
    INTEGER :: ist      !< system status code

    !---
    CALL message(TRIM(thismodule),'Destruction of ECHAM physics state started.')

    ndomain = SIZE(prm_field)

    DO jg = 1,ndomain
      CALL delete_var_list( prm_field_list(jg) )
      CALL delete_var_list( prm_tend_list (jg) )
    ENDDO

    DEALLOCATE( prm_field_list, prm_tend_list, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of prm_field/tend list array failed')

    DEALLOCATE( prm_field, prm_tend, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of prm_field/tend array failed')

    CALL message(TRIM(thismodule),'Destruction of ECHAM physics state finished.')

  END SUBROUTINE destruct_echam_phy_state
  !-------------
  !>
  !!
  !!
  SUBROUTINE new_echam_phy_field_list( kproma, klev, kblks, ktracer,      &
                                     & ctracer_list, ksfc_type, listname, &
                                     & prefix, field_list, field          )

    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer, ksfc_type  !< dimension sizes

    CHARACTER(len=*),INTENT(IN)    :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list

    TYPE(t_var_list),       INTENT(INOUT) :: field_list
    TYPE(t_echam_phy_field),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shapesfc(3)
!0!    INTEGER :: shape4d(4)
    INTEGER :: nbits, morebits, jsfc, jtrc

    CHARACTER(LEN=1) :: csfc

    nbits = 16
    morebits = 24 

    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)

    ! Register a field list and apply default settings

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! &       field% u         (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind')
    grib2_desc = t_grib2_var(0, 2, 2, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'u', field%u,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% v         (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind')
    grib2_desc = t_grib2_var(0, 2, 3, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'v', field%v,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% vor       (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vorticity', 's-1', 'relative vorticity')
    grib2_desc = t_grib2_var(0, 2, 12, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vo', field%vor,                          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% temp      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('temperature', 'K', 'temperature')
    grib2_desc = t_grib2_var(0, 0, 0, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'t', field%temp,                          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% tv        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vtemp', field%tv,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! OZONE 
    ! &       field% o3        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('ozone', 'g/g', 'ozone mixing ratio')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'o3', field%o3,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% q         (nproma,nlev  ,nblks,ntracer),  &
    CALL add_var( field_list, prefix//'q', field%q,                                    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                        &
                & t_cf_var('q', 'kg kg-1', ''),                                &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims = (/kproma,klev,kblks,ktracer/),                       &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    ALLOCATE(field%q_ptr(ktracer))                                                   
    DO jtrc = 1,ktracer                                                                   
      CALL add_ref( field_list, prefix//'q',                                     &
                  & prefix//'q_'//ctracer_list(jtrc:jtrc), field%q_ptr(jtrc)%p,  &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                        &
                  & t_cf_var('q_'//ctracer_list(jtrc:jtrc), 'kg kg-1', ''),      &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=(/kproma,klev,kblks/))                                        
    END DO                                                                                

    ! &       field% qx        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('condensated_water', 'kg kg-1', 'cloud water + cloud ice')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'qx', field%qx,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% omega     (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vertical_velocity', 'Pa s-1', 'vertical velocity')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'omega', field%omega,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% geom      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface')
    grib2_desc = t_grib2_var(0, 3, 5, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ghm', field%geom,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% presm_old (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step')
    grib2_desc = t_grib2_var(0, 3, 0, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presm_old', field%presm_old,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% presm_new (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step')
    grib2_desc = t_grib2_var(0, 3, 0, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presm_new', field%presm_new,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )


    !-- Variables defined at layer interfaces --

    shape3d = (/kproma,klev+1,kblks/)

    ! &       field% geoi      (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface')
    grib2_desc = t_grib2_var(0, 3, 5, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ghi', field%geoi,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% presi_old (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step')
    grib2_desc = t_grib2_var(0, 3, 0, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presi_old', field%presi_old,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% presi_new (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step')
    grib2_desc = t_grib2_var(0, 3, 0, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presi_new', field%presi_new,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )

    !------------------
    ! Radiation
    !------------------
    ! 2D variables

   !ALLOCATE( field% cosmu0    (nproma,       nblks),          &
    cf_desc    = t_cf_var('cosmu0', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cosmu0', field%cosmu0,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% flxdwswtoa(nproma,       nblks),          &
    cf_desc    = t_cf_var('flxdwswtoa', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'flxdwswtoa', field%flxdwswtoa,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albvisdir (nproma,       nblks),          &
    cf_desc    = t_cf_var('albvisdir', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir', field%albvisdir,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albvisdif (nproma,       nblks),          &
    cf_desc    = t_cf_var('albvisdif', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif', field%albvisdif,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albnirdir (nproma,       nblks),          &
    cf_desc    = t_cf_var('albnirdir', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir', field%albnirdir,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albnirdif (nproma,       nblks),          &
    cf_desc    = t_cf_var('albnirdif', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif', field%albnirdif,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% vissfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('vissfc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vissfc', field%vissfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% visdffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('visdffsfc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'visdffsfc', field%visdffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% nirsfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('nirsfc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'nirsfc', field%nirsfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% nirdffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('nirdffsfc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'nirdffsfc', field%nirdffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% parsfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('parsfc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'parsfc', field%parsfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% pardffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('pardffsfc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pardffsfc', field%pardffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )


    !---- 3D variables defined at layer interfaces ----

    shape3d = (/kproma,klev+1,kblks/)

    ! &       field% emterclr  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('emterclr', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'emterclr', field%emterclr,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% emterall  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('emterall', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'emterall', field%emterall,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% trsolclr  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('trsolclr', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'trsolclr', field%trsolclr,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% trsolall  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('trsolall', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'trsolall', field%trsolall,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d )


    !-------------------------
    ! Cloud and precipitation
    !-------------------------
    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)

   !ALLOCATE( field% aclc   (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('ACLC', 'm2 m-2', 'cloud area fraction, instantaneous')
    grib2_desc = t_grib2_var(0,6,1, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aclc', field%aclc,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% aclcac (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('ACLCAC', 's m2 m-2', &
               & 'cloud area fraction accumulated over output interval')
    grib2_desc = t_grib2_var(0,6,255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aclcac', field%aclcac,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% aclcov (nproma,       nblks), &
    cf_desc    = t_cf_var('ACLCOV', 's m2 m-2', &
               & 'total cloud cover accumulated over output interval')
    grib2_desc = t_grib2_var(0,6, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aclcov', field%aclcov,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% acdnc  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('acdnc', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'acdnc', field%acdnc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% xvar   (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('variance_of_total_water', '', 'subgrid variance of total water')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xvar', field%xvar,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% xskew  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('skewness_of_total_water', '', 'skewness of total water')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xskew', field%xskew,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% relhum (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('relative_humidity', '', 'relative humidity')
    grib2_desc = t_grib2_var(0, 1, 1, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'r', field%relhum,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% aprl   (nproma,       nblks), &
    cf_desc    = t_cf_var('APRL', 'kg m-2', &
               & 'large-scale precipitation amount accumulated over output interval')
    grib2_desc = t_grib2_var(0, 1, 9, morebits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ncpcp', field%aprl,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% aprc   (nproma,       nblks), &
    cf_desc    = t_cf_var('APRC', 'kg m-2', &
               & 'convective precipitation amount accumulated over output interval')
    grib2_desc = t_grib2_var(0, 1, 37, morebits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cprat', field%aprc,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% aprs   (nproma,       nblks), &
    cf_desc    = t_cf_var('APRS', 'kg m-2', &
               & 'snowfall (large-scale + convective) amount accumulated over output interval')
    grib2_desc = t_grib2_var(0, 1, 66, morebits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'sprate', field%aprs,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% rsfl   (nproma,       nblks), &
    cf_desc    = t_cf_var('RSFL', 'kg m-2 s-1',    &
               & 'instantaneous large-scale precipitation flux (water)')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsfl', field%rsfl,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% rsfc   (nproma,       nblks), &
    cf_desc    = t_cf_var('RSFC', 'kg m-2 s-1',    &
               & 'instantaneous convective precipitation flux (water)')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsfc', field%rsfc,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% ssfl   (nproma,       nblks), &
    cf_desc    = t_cf_var('SSFL', 'kg m-2 s-1',    &
               & 'instantaneous large-scale precipitation flux (snow)')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ssfl', field%ssfl,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% ssfc   (nproma,       nblks), &
    cf_desc    = t_cf_var('SSFC', 'kg m-2 s-1',    &
               & 'instantaneous convective precipitation flux (snow)')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ssfc', field%ssfc,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field%  qvi   (nproma,       nblks), &
    cf_desc    = t_cf_var('total_vapour', 'kg m-2', 'vertically integrated water vapour')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'qvi', field%qvi,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% xlvi   (nproma,       nblks), &
    cf_desc    = t_cf_var('total_cloud_water', 'kg m-2',&
               & 'vertically integrated cloud water'    )
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xlvi', field%xlvi,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% xivi   (nproma,       nblks), &
    cf_desc    = t_cf_var('total_cloud_ice', 'kg m-2',&
               & 'vertically integrated cloud ice'    )
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xivi', field%xivi,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% rintop (nproma,       nblks), &
    cf_desc    = t_cf_var('rintop', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rintop', field%rintop,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% rtype  (nproma,       nblks), &
    cf_desc    = t_cf_var('convection_type', '', 'convection_type (0...3)')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rtype', field%rtype,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% topmax (nproma,       nblks), &
    cf_desc    = t_cf_var('topmax', 'Pa', 'maximum height of convective cloud tops')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'topmax', field%topmax,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% tke    (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('turbulent_kinetic_energy', 'm2 s-2', 'turbulent kinetic energy')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tke', field%tke,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% thvsig (nproma,       nblks), &
    cf_desc    = t_cf_var('thvsig', 'K', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'thvsig', field%thvsig,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !--------------------
    ! Turbulence
    !--------------------
    IF (get_lvdiff()) THEN

      shape2d  = (/kproma,            kblks/)
      shape3d  = (/kproma, klev,      kblks/)
     !shapesfc = (/kproma, ksfc_type, kblks/)
      shapesfc = (/kproma, kblks, ksfc_type/)

     !ALLOCATE( field% ri     (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('richardson_number', ' ', 'moist Richardson number')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ri', field%ri,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% mixlen (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('mixing_length', 'm', 'mixing_length')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'mixlen', field%mixlen,                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% thvvar (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('thvvar', 'K2',                           &
                 & 'subgrid variance of virtual potential temperature')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'thvvar', field%thvvar,                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% tkem0  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('tke', 'm2 s-2', 'TKE at step t')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'tkem0', field%tkem0,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% tkem1  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('tke', 'm2 s-2', 'TKE at step t-dt')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'tkem1', field%tkem1,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

     !---------
     !ALLOCATE( field% cfm    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_momentum', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfm', field%cfm,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfm_tile(nproma,nsfc_type,nblks), &
      CALL add_var( field_list, prefix//'cfm_tile', field%cfm_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('turb_exchng_coeff_momentum', '', ''),              &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

      ALLOCATE(field%cfm_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( field_list, prefix//'cfm_tile',                              &
                    & prefix//'cfm_tile_'//csfc, field%cfm_tile_ptr(jsfc)%p,       &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_momentum_'//csfc, '', ''),       &
                    & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

      !---------
      ! &       field% cfh    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_heat', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfh', field%cfh,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfh_tile(nproma,nsfc_type,nblks), &
      CALL add_var( field_list, prefix//'cfh_tile', field%cfh_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('turb_exchng_coeff_heat', '', ''),                  &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

      ALLOCATE(field%cfh_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( field_list, prefix//'cfh_tile',                              &
                    & prefix//'cfh_tile_'//csfc, field%cfh_tile_ptr(jsfc)%p,       &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_heat_'//csfc, '', ''),           &
                    & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

      !---------
      ! &       field% cfv    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_water_var', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfv', field%cfv,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cftke  (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_tke', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cftke', field%cftke,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfthv  (nproma,nlev,     nblks)  )
      cf_desc    = t_cf_var('turb_exchng_coeff_thv', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfthv', field%cfthv,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

     !ALLOCATE( field% coriol (nproma,nblks),                &
      cf_desc    = t_cf_var('Coriolis_param', 's-1', 'Coriolis parameter')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'coriol', field%coriol,                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ghpbl  (nproma,nblks),                &
      cf_desc    = t_cf_var('geopot_pbl_top', '', 'geopotential of PBL top')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ghpbl', field%ghpbl,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      !-----------------------------------
      ! &       field% z0m(nproma,nblks), &
      cf_desc    = t_cf_var('z0m', '', 'aerodynamic roughness length, grid box mean')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'z0m', field%z0m,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% z0m_tile(nproma,nsfc_type,nblks), &
      CALL add_var( field_list, prefix//'z0m_tile', field%z0m_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('z0m_tile', '', 'aerodynamic roughness length'),    &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

      ALLOCATE(field%z0m_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( field_list, prefix//'z0m_tile',                              &
                    & prefix//'z0m_tile_'//csfc, field%z0m_tile_ptr(jsfc)%p,       &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('z0m_tile_'//csfc, '',''),                          &
                    & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO
      !-----------------------------------

      ! &       field% ustar  (nproma,nblks),                &
      cf_desc    = t_cf_var('fricktion_velocity', 'm s-1', 'friction velocity')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ustar', field%ustar,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% kedisp (nproma,nblks),                &
      cf_desc    = t_cf_var('KE dissipation rate', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'kedisp', field%kedisp,                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ocu    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_u', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, morebits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ocu', field%ocu,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ocv    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_v', '', '')
      grib2_desc = t_grib2_var(255, 255, 255, morebits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ocv', field%ocv,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ENDIF ! lvdiff

    !-----------------------
    ! Surface
    !-----------------------
   !ALLOCATE( field% lfland (nproma, nblks),                 &
    cf_desc    = t_cf_var('lfland', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lfland', field%lfland,        &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,               &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

    ! &       field% lfglac (nproma, nblks),                 &
    cf_desc    = t_cf_var('lfglac', '', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lfglac', field%lfglac,         &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

  ! ALLOCATE( field% lfland (kproma, kblks), &
  !         & field% lfglac (kproma, kblks)  ) 

    ! &       field% lsmask (nproma, nblks),                 &
    cf_desc    = t_cf_var('land_cover', '', 'land cover')
    grib2_desc = t_grib2_var(2, 0, 0, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'land', field%lsmask,                   &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% glac   (nproma, nblks),                 &
    cf_desc    = t_cf_var('glacier_cover', '', 'fraction of land covered by glaciers')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'glac', field%glac,                     &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% seaice (nproma, nblks),                 &
    cf_desc    = t_cf_var('sea_ice_cover', '', 'fraction of ocean covered by sea ice')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'seaice', field%seaice,                 &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% icefrc (nproma, nblks),                 &
    cf_desc    = t_cf_var('ice_cover', '', 'ice cover given as fraction of grid box')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'icefrc', field%icefrc,                 &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !-----------------------------------
    ! &       field% tsfc(nproma,nblks), &
    cf_desc    = t_cf_var('surface_temperature', '', 'surface temperature')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc', field%tsfc,                     &
              & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% tsfc_tile(nproma,nsfc_type,nblks), &
    CALL add_var( field_list, prefix//'tsfc_tile', field%tsfc_tile,            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var('tsfc_tile', '', 'skin temperature'),               &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    ALLOCATE(field%tsfc_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 
      CALL add_ref( field_list, prefix//'tsfc_tile',                             &
                  & prefix//'tsfc_tile_'//csfc, field%tsfc_tile_ptr(jsfc)%p,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('tsfc_tile_'//csfc, '', ''),                        &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shape2d )
    END DO
    !-----------------------------------

    ! &       field% qs_sfc_tile (nproma,nsfc_type, nblks), &
    CALL add_var( field_list, prefix//'qs_sfc_tile', field%qs_sfc_tile,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var('qs_sfc_tile', '', ''),                             &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    ALLOCATE(field%qs_sfc_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 
      CALL add_ref( field_list, prefix//'qs_sfc_tile',                           &
                  & prefix//'qs_sfc_tile_'//csfc, field%qs_sfc_tile_ptr(jsfc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('qs_sfc_tile_'//csfc, '', ''),                      &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shape2d )
    END DO

    !---------------------------
    ! Surface fluxes
    !---------------------------
    ! Accumulated gridbox mean

    CALL add_var( field_list, prefix//'evap_ac', field%evap_ac,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                  &
                & t_cf_var('evap_ac', 'kg m-2', 'evaporation'//           &
                & ' accumulated over output interval'),                   &
                & t_grib2_var(2,0,6,morebits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d                                           )

    CALL add_var( field_list, prefix//'lhflx_ac', field%lhflx_ac,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                  &
                & t_cf_var('lhflx_ac', 'W m-2 s', 'latent heat flux'//    &
                & 'accumulated over output interval'),                    &
                & t_grib2_var(2, 0, 6, nbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d                                           )

    CALL add_var( field_list, prefix//'shflx_ac', field%shflx_ac,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                  &
                & t_cf_var('shflx_ac', 'W m-2 s', 'sensible heat flux'//  &
                & 'accumulated over output interval'),                    &
                & t_grib2_var(2, 0, 6, nbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d                                           )

    !---------------------------------
    ! Instantaneous values over tiles

    CALL add_var( field_list, prefix//'evap_tile', field%evap_tile,       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                  &
                & t_cf_var('evap_tile', '', ''),                          &
                & t_grib2_var(2, 0, 6, nbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.      )

    CALL add_var( field_list, prefix//'lhflx_tile', field%lhflx_tile,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                  &
                & t_cf_var('lhflx_tile', 'W m-2', 'latent heat flux'),    &
                & t_grib2_var(2, 0, 6, nbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.      )

    CALL add_var( field_list, prefix//'shflx_tile', field%shflx_tile,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                  &
                & t_cf_var('shflx_tile', 'W m-2', 'sensible heat flux'),  &
                & t_grib2_var(2, 0, 6, nbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.      )

    ALLOCATE(field%evap_tile_ptr(ksfc_type))
    ALLOCATE(field%lhflx_tile_ptr(ksfc_type))
    ALLOCATE(field%shflx_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 

      CALL add_ref( field_list, prefix//'evap_tile',                             &
                  & prefix//'evap_tile_'//csfc, field%evap_tile_ptr(jsfc)%p,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('evap_tile_'//csfc, '', ''),                        &
                  & t_grib2_var(2,0,6, nbits, GRID_REFERENCE, GRID_CELL),        &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'lhflx_tile',                            &
                  & prefix//'lhflx_tile_'//csfc, field%lhflx_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('lhflx_tile_'//csfc, '', ''),                       &
                  & t_grib2_var(2,0,6, nbits, GRID_REFERENCE, GRID_CELL),        &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'shflx_tile',                            &
                  & prefix//'shflx_tile_'//csfc, field%shflx_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                  & t_cf_var('shflx_tile_'//csfc, '', ''),                       &
                  & t_grib2_var(2,0,6, nbits, GRID_REFERENCE, GRID_CELL),        &
                  & ldims=shape2d )
    END DO

    !-----------------------------------------
    ! wind stress, accumulated grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'u_stress_ac', field%u_stress_ac,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var('u_stress_ac', 'N m-2 s', 'surface wind stress'//   &
                & 'accumulated over output interval'),                         &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shape2d                                                )

    CALL add_var( field_list, prefix//'v_stress_ac', field%v_stress_ac,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var('v_stress_ac', 'N m-2 s', 'surface wind stress'//   &
                & 'accumulated over output interval'),                         &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shape2d                                                )

    ! wind stress, instantaneous tile values 

    CALL add_var( field_list, prefix//'u_stress_tile', field%u_stress_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var('u_stress_tile', '', ''),                           &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    CALL add_var( field_list, prefix//'v_stress_tile', field%v_stress_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var('v_stress_tile', '', ''),                           &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    ALLOCATE(field%u_stress_tile_ptr(ksfc_type))
    ALLOCATE(field%v_stress_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 

      CALL add_ref( field_list, prefix//'u_stress_tile',                            &
                  & prefix//'u_stress_tile_'//csfc, field%u_stress_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                          &
                  & t_cf_var('u_stress_tile_'//csfc, '', ''),                       &
                  & t_grib2_var(255,255,255, nbits, GRID_REFERENCE, GRID_CELL),     &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'v_stress_tile',                            &
                  & prefix//'v_stress_tile_'//csfc, field%v_stress_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                          &
                  & t_cf_var('v_stress_tile_'//csfc, '', ''),                       &
                  & t_grib2_var(255,255,255, nbits, GRID_REFERENCE, GRID_CELL),     &
                  & ldims=shape2d )
    END DO

  END SUBROUTINE new_echam_phy_field_list
  !-------------
  !>
  !!
  !!
  SUBROUTINE new_echam_phy_tend_list( kproma, klev, kblks, ktracer,   &
                                    & ctracer_list, listname, prefix, &
                                    & tend_list, tend )

    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer  !< dimension sizes

    CHARACTER(len=*),INTENT(IN)    :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list

    TYPE(t_var_list)      ,INTENT(INOUT) :: tend_list
    TYPE(t_echam_phy_tend),INTENT(INOUT) :: tend

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d(3), shape_trc(4)
    INTEGER :: nbits, jtrc
    !------------------------------

    nbits = 16 ! "entropy" of horizontal slice

    shape3d   = (/kproma, klev, kblks/)
    shape_trc = (/kproma, klev, kblks, ktracer/)

    CALL new_var_list( tend_list, TRIM(listname) )
    CALL default_var_list_settings( tend_list, lrestart=.FALSE. )

    !------------------------------
    ! Temperature tendencies
    !------------------------------
    ! &       tend% temp      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp', tend%temp,                    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_radsw(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_radsw', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_radsw', tend%temp_radsw,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_radlw(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_radlw', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_radlw', tend%temp_radlw,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_cld  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_cloud', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_cld', tend%temp_cld,            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_convective', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_cnv', tend%temp_cnv,            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_turbulent', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_vdf', tend%temp_vdf,            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_Hines_gw', 'K s-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_gwh', tend%temp_gwh,            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! U-wind tendencies
    !------------------------------
   !ALLOCATE( tend%    u      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u', tend%u,                          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_convective', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_cnv', tend%u_cnv,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_turbulent', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_vdf', tend%u_vdf,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_Hines_gw', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_gwh', tend%u_gwh,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! V-wind tendencies
    !------------------------------
    ! &       tend%    v      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v', tend%v,                          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_cnv', tend%v_cnv,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_turbulent', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_vdf', tend%v_vdf,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_Hines_gw', 'm s-2', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_gwh', tend%v_gwh,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! Detrainment
    !------------------------------
    ! &       tend%    x_dtr  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('detrain_rate', 's-1', '')
    grib2_desc = t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'x_dtr', tend%x_dtr,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !-------------------
    ! Tracer tendencies
    !-------------------
    ! Tracer arrays for (model) internal use                                               
                                                                                            
    CALL add_var( tend_list, prefix//'q', tend%q,                              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                        &
                & t_cf_var('tend_q', 's-1', 'tracer tendency'),                &
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    CALL add_var( tend_list, prefix//'q_cld', tend%q_cld,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                           &
                & t_cf_var('tend_q_cld', 's-1', 'tracer tendency condensational'),&           
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),   &
                & ldims = shape_trc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.              )

    CALL add_var( tend_list, prefix//'q_cnv', tend%q_cnv,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                        &
                & t_cf_var('tend_q_cnv', 's-1', 'tracer tendency convective'), &           
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    CALL add_var( tend_list, prefix//'q_vdf', tend%q_vdf,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                        &
                & t_cf_var('tend_q_vdf', 's-1', 'tracer tendency turbulent'),  &           
                & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    ! Referrence to individual tracer, for I/O                                            
                                                                                            
    ALLOCATE(tend%     q_ptr(ktracer))
    ALLOCATE(tend% q_cld_ptr(ktracer))
    ALLOCATE(tend% q_cnv_ptr(ktracer))
    ALLOCATE(tend% q_vdf_ptr(ktracer))

    DO jtrc = 1,ktracer                                                                   
                                                                                          
      CALL add_ref( tend_list, prefix//'q',                                        &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc), tend%q_ptr(jtrc)%p,      &       
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                          &       
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc), 's-1', ''),        &       
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),  &       
                  & ldims=(/kproma,klev,kblks/))                                        
                                                                                          
      CALL add_ref( tend_list, prefix//'q_cld',                                           &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_cld', tend%q_cld_ptr(jtrc)%p, &       
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                                 &       
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_cld', 's-1', ''),       &       
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),         &       
                  & ldims=(/kproma,klev,kblks/))                                        

      CALL add_ref( tend_list, prefix//'q_cnv',                                           &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_cnv', tend%q_cnv_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                                 &
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_cnv', 's-1', ''),       &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/))                                        

      CALL add_ref( tend_list, prefix//'q_vdf',                                           &
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_vdf', tend%q_vdf_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                                 &
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_vdf', 's-1', ''),       &
                  & t_grib2_var(255, 255, 255, nbits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/))                                        
    END DO                                                                                

  END SUBROUTINE new_echam_phy_tend_list
  !-------------

END MODULE mo_echam_phy_memory
