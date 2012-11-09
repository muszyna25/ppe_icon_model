#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif

!>
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the ECHAM physics.
!!  This module should be an analogon to 'mo_nonhydro_types.f90'

!!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps
!!
!! @author Kristina Froehlich, DWD
!! @author Marco Giorgetta, MPI-M
!!
!! $Id: n/a$
!!
!! @par Revision History
!! Initial  by Kristina Froehlich (2009-06-10)
!! Memory allocation method changed from explicit allocation to Luis' 
!! infrastructure by Kristina Froehlich (MPI-M, 2011-04-27)
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
MODULE mo_nwp_phy_types

  USE mo_kind,                ONLY: wp
  USE mo_fortran_tools,       ONLY: t_ptr_2d3d
  USE mo_model_domain,        ONLY: t_patch
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_nwp_parameters,      ONLY: t_phy_params
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var

  IMPLICIT NONE
  PRIVATE


  ! !VERSION CONTROL:
  CHARACTER(len=*), PARAMETER :: version = &
    & '$Id$'

  !public interface
  !
  !types
  PUBLIC :: t_nwp_phy_diag
  PUBLIC :: t_nwp_phy_tend


  !
  !!data structure defining model states
  !
  !!diagnostic variables
  !

  TYPE t_nwp_phy_diag

    TYPE(t_ptr_2d3d),ALLOCATABLE :: tot_ptr(:)  !< pointer array: one pointer for each tot var (grid+subgrid)
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tci_ptr(:)  !< pointer array: total column-integrated values
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tav_ptr(:)  !< pointer array: average of tci

    TYPE(t_ptr_2d3d),ALLOCATABLE :: cfm_ptr(:)  !< pointer array: average of cfm
    TYPE(t_ptr_2d3d),ALLOCATABLE :: cfh_ptr(:)  !< pointer array: average of cfh
    TYPE(t_ptr_2d3d),ALLOCATABLE :: z0m_ptr(:)  !< pointer array: average of z0m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: albvisdif_t_ptr(:)!< pointer array: tile-specific albedo
    TYPE(t_ptr_2d3d),ALLOCATABLE :: swflxsfc_t_ptr(:) !< pointer array: shortwave net flux at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: lwflxsfc_t_ptr(:) !< pointer array: longwave net flux at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tcm_t_ptr(:) !< pointer array: turbulent transfer coefficients for momentum
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tch_t_ptr(:) !< pointer array: turbulent transfer coefficients for heat
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tfm_t_ptr(:) !< pointer array: factor of laminar transfer of momentum
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tfh_t_ptr(:) !< pointer array: factor of laminar transfer of heat 
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tfv_t_ptr(:) !< pointer array: laminar reduction factor for evaporation
    TYPE(t_ptr_2d3d),ALLOCATABLE :: gz0_t_ptr(:) !< pointer array: roughness length * gravity
    TYPE(t_ptr_2d3d),ALLOCATABLE :: t_2m_t_ptr(:)  !< pointer array: temperature at 2m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: qv_2m_t_ptr(:) !< pointer array: specific water vapor content at 2m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: rh_2m_t_ptr(:) !< pointer array: relative humidity at 2m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: td_2m_t_ptr(:) !< pointer array: dew point at 2m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: u_10m_t_ptr(:) !< pointer array: zonal wind at 2am
    TYPE(t_ptr_2d3d),ALLOCATABLE :: v_10m_t_ptr(:) !< pointer array: meridional wind at 2m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: shfl_s_t_ptr(:) !< pointer array: surface sensible heat flux 
    TYPE(t_ptr_2d3d),ALLOCATABLE :: lhfl_s_t_ptr(:) !< pointer array: surface latent heat flux

    REAL(wp), POINTER ::  &
      &   rain_gsp_rate(:,:),  & !! grid-scale surface rain rate                         [kg/m2/s]
      &   snow_gsp_rate(:,:),  & !! grid_scale surface snow rate                         [kg/m2/s]
      &   rain_con_rate(:,:),  & !! convective surface rain rate                         [kg/m2/s]
      &   snow_con_rate(:,:),  & !! convective surface snow_rate                         [kg/m2/s]
      &   rain_con_rate_3d(:,:,:),  & !! 3d convective rain rate                         [kg/m2/s]
      &   snow_con_rate_3d(:,:,:),  & !! 3d convective snow_rate                         [kg/m2/s]
      &   rain_gsp(:,:),       & !! accumulated grid-scale surface rain                  [kg/m2]
      &   snow_gsp(:,:),       & !! accumulated grid_scale surface snow                  [kg/m2]
      &   rain_con(:,:),       & !! accumulated convective surface rain                  [kg/m2]
      &   snow_con(:,:),       & !! accumulated convective surface snow                  [kg/m2]
      &   tot_prec(:,:),       & !! accumulated grid-scale plus convective surface       [kg/m2]
                                 !! total precipitation
      &   tot_prec_rate_avg(:,:),   & !! average since model start of                    [kg/m2/s]
                                !! grid-scale plus convective surface 
                                !! total precipitation rate
      &   con_prec_rate_avg(:,:),   & !! average since model start of                    [kg/m2/s]
                                !! convective surface precipitation rate
      &   gsp_prec_rate_avg(:,:),   & !! average since model start of                    [kg/m2/s]
                                !! grid-scale surface precipitation rate
      &   cape    (:,:),       & !! convective available energy
      &   con_gust(:,:),       & !! convective gusts near surface
      &   con_udd(:,:,:,:),    & !!(nproma,nlev,nblks,8) convective up/downdraft fields
                                !! 1= convective updraft mass flux (pmfu)
                                !! 2= convective downdraft mass flux (pmfd)
                                !! 3= updraft   detrainment rate  (pmfude_rate)
                                !! 4= downdraft   detrainment rate (pmfdde_rate)
                                !! 5= temperature in updraft region (ptu)
                                !! 6= humidity in updraft region (pqu)
                                !! 7= condensate in updraft region (plu)
      &  rain_upd(:,:),        & !! total precipitation produced in updrafts [kg/m2/s]
      &  shfl_s(:,:),          & !! sensible heat flux (surface) ( W/m2)
      &  shfl_s_t(:,:,:),      & !! sensible heat flux (surface) ( W/m2)
      &  lhfl_s(:,:),          & !! latent   heat flux (surface) ( W/m2)
      &  lhfl_s_t(:,:,:),      & !! latent   heat flux (surface) ( W/m2)
      &  qhfl_s(:,:),          & !!      moisture flux (surface) ( Kg/m2/s)
                                !!      = evaporation rate at surface
      &  shfl_s_a(:,:),        & !! average or accumulated since model start of shfl_s [W/m2]
      &  lhfl_s_a(:,:),        & !! average or accumulated since model start of lhfl_s [W/m2]
      &  qhfl_s_avg(:,:),      & !! average since model start of qhfl_s ( Kg/m2/s) 
                                !! = average of evaporation rate at surface
      &  tot_cld(:,:,:,:),     & !! total cloud variables (cc,qv,qc,qi)
      &  tot_cld_vi(:,:,:),    & !! vertically integrated tot_cld (cc,qv,qc,qi) 
                                !! for cc, instead of the vertically integrated value, 
                                !! this is the cloud cover assuming maximum-random overlap
      &  tot_cld_vi_avg(:,:,:),& !! average since model start of the 
                                !! vertically integrated tot_cld (cc,qv,qc,qi)  
      &  cosmu0(:,:),          & !! cosine of solar zenith angle
      &  albvisdif(:,:),       & !! surface albedo for visible range, diffuse
      &  albvisdif_t(:,:,:),   & !! tile-based surface albedo for visible range, diffuse
      &  vio3(:,:),            & !! vertically integrated ozone amount (Pa O3)
      &  hmo3(:,:),            & !! height of O3 maximum (Pa)
      &  flxdwswtoa(:,:),      & !! downward shortwave flux at TOA [W/m2]
      &  tsfctrad(:,:),        & !! surface temperature at trad [K]
      &  lwflxclr(:,:,:),      & !! longwave clear-sky net flux [W/m2]
      &  lwflxall(:,:,:),      & !! longwave net flux           [W/m2]
      &  lwflxsfc(:,:),        & !! longwave net flux at surface [W/m2]
      &  lwflxsfc_t(:,:,:),    & !! tile-based longwave net flux at surface [W/m2]
      &  trsolclr(:,:,:),      & !! shortwave clear-sky net tranmissivity []
      &  trsolall(:,:,:),      & !! shortwave net tranmissivity []
      &  swflxsfc(:,:),        & !! shortwave net flux at surface [W/m2]
      &  swflxsfc_t(:,:,:),    & !! tile-based shortwave net flux at surface [W/m2]
      &  swflxtoa(:,:),        & !! shortwave net flux at toa [W/m2]
      &  lwflxsfc_a(:,:),      & !! longwave net flux at surface [W/m2], accumulated or mean since last output
      &  swflxsfc_a(:,:),      & !! shortwave net flux at surface [W/m2], accumulated or mean since last output
      &  lwflxtoa_a(:,:),      & !! longwave net flux at toa [W/m2], accumulated or mean since last output
      &  swflxtoa_a(:,:),      & !! shortwave net flux at toa [W/m2], accumulated or mean since last output
                                !! _a means average values if lflux_avg=.TRUE.
                                !! and accumulated values if lflux_avg=.FALSE., default is .FALSE.
      &  acdnc(:,:,:)            !! cloud droplet number concentration [1/m**3]



    !> Parameter fields for turbulence
    REAL(wp), POINTER ::  &
      rcld(:,:,:)      ,    & !> standard deviation of the saturation deficit    --
      tcm(:,:)        ,    & !! turbulent transfer coefficients for momentum    --
      tch(:,:)        ,    & !! turbulent transfer coefficients for heat        --
      tfm(:,:)        ,    & !! factor of laminar transfer of momentum          --
      tfh(:,:)        ,    & !! factor of laminar transfer of scalars           --
      tfv(:,:)        ,    & !! laminar reduction factor for evaporation        --
      gz0(:,:),            & !! roughness length * g of the vertically not
                                !! resolved canopy                               (m2/s2)
      sai(:,:),            & !! surface area index                            ( 1 )
      tai(:,:),            & !! transpiration area index                      ( 1 )
      eai(:,:),            & !! (evaporative) earth area index                ( 1 )
      tkvm(:,:,:),         & !! turbulent diffusion coefficients for momentum (m/s2 )
      tkvh(:,:,:),         & !! turbulent diffusion coefficients for heat     (m/s2 )
      t_2m(:,:)       ,    & !! temperature in 2m                             (  K  )
      t_2m_s6avg(:,:),     & !! 6 hourly sample 2 m temperature average       (  K  )
      qv_2m (:,:)     ,    & !! specific water vapor content in 2m            (kg/kg)
      qv_2m_s6avg(:,:),    & !! 6 hourly sample 2 m specific water vapor content average   (kg/kg)
      td_2m (:,:)     ,    & !! dew-point in 2m                               (  K  )
      rh_2m (:,:)     ,    & !! relative humidity in 2m                       (  %  )
      u_10m (:,:)     ,    & !! zonal wind in 10m                             ( m/s )
      v_10m (:,:)     ,    & !! meridional wind in 10m                        ( m/s )
      u_10m_s6avg (:,:),   & !! 6 hourly sample 10m zonal wind  average       ( m/s )
      v_10m_s6avg (:,:),   & !! 6 hourly sample 10m  meridional wind average  ( m/s )
      edr   (:,:,:)    ,   & !! eddy dissipation rate
      tcm_t(:,:,:)     ,   & !! turbulent transfer coefficients for momentum    --
      tch_t(:,:,:)     ,   & !! turbulent transfer coefficients for heat        --
      tfm_t(:,:,:)     ,   & !! factor of laminar transfer of momentum          --
      tfh_t(:,:,:)     ,   & !! factor of laminar transfer of scalars           --
      tfv_t(:,:,:)     ,   & !! laminar reduction factor for evaporation        --
      gz0_t(:,:,:)     ,   & !! roughness length * g                          (m2/s2)
      t_2m_t (:,:,:)   ,   & !! temperature at 2m                             (  K  )
      qv_2m_t(:,:,:)   ,   & !! specific water vapor content at 2m            (kg/kg)
      td_2m_t(:,:,:)   ,   & !! dew-point at 2m                               (  K  )
      rh_2m_t(:,:,:)   ,   & !! relative humidity at 2m                       (  %  )
      u_10m_t(:,:,:)   ,   & !! zonal wind at 10m                             ( m/s )
      v_10m_t(:,:,:)         !! meridional wind at 10m                        ( m/s )

    ! need only for vdiff (and some for EDMF)
    REAL(wp),POINTER :: &
      & ri      (:,:,:),    &!< moist Richardson number at layer interfaces
      & mixlen  (:,:,:),    &!< mixing length at layer interfaces
      & thvvar  (:,:,:),    &!< variance of virtual potential temperature at layer interfaces.
                                !< Computed in "vdiff" by solving a prognostic equation of
                                !< the variance. Used for getting "thvsig".
      & z0m_tile(:,:,:),    &!< aerodynamic roughness length
                                !< (grid-box mean and over each surface type)
      & z0m     (:,:)  ,    &!< aerodynamic roughness length
      & ustar   (:,:)  ,    &!<
      & kedisp  (:,:)  ,    &!< time-mean (or integrated?)
                                !< vertically integrated dissipation of kinetic energy
      & ocu     (:,:)  ,    &!< eastward  velocity of ocean surface current
      & ocv     (:,:)        !< northward velocity of ocean surface current


    REAL(wp),POINTER :: &
      & cfm    (:,:,:),     &!< turbulent exchange coefficient
      & cfm_tile(:,:,:),    &!< turbulent exchange coefficient
      & cfh    (:,:,:),     &!< turbulent exchange coefficient
      & cfh_tile(:,:,:),    &!< turbulent exchange coefficient
      & cfv    (:,:,:),     &!< turbulent exchange coefficient
      & cftke  (:,:,:),     &!< turbulent exchange coefficient
      & cfthv  (:,:,:),     &!< turbulent exchange coefficient
      & ghpbl  (:,:)         !< geopotential of the top of the atmospheric boundary layer


    ! for old aerosol climatology from COSMO (to be used with inwp_radiation==2)
    REAL(wp),POINTER :: &
      & aersea  (:,:),      &
      & aerlan  (:,:),      &
      & aerurb  (:,:),      &
      & aerdes  (:,:)

    INTEGER, POINTER :: &
      &  mbas_con(:,:),     & !< cloud base level index
      &  mtop_con(:,:),     & !< cloud top  level index
      &  ktype   (:,:)        !< Type of convection

    LOGICAL, POINTER :: &
      & locum     (:,:)       !< convective  activity indicator

    !> (Optional:) Additional diagnostic fields:
    REAL(wp), POINTER ::  &
      rh(:,:,:)               !> relative humidity

  END TYPE t_nwp_phy_diag
  !
  ! !---tendencies of type global!
  !
  TYPE t_nwp_phy_tend

    REAL(wp),POINTER ::  &
      ddt_temp_radsw  (:,:,:)  ,& !! Temp-tendency from shortwave radiation
      ddt_temp_radlw  (:,:,:)  ,& !! Temp-tendency from longwave radiation
      ddt_temp_turb   (:,:,:)  ,& !! Temp-tendency from turbulence
      ddt_temp_drag   (:,:,:)  ,& !! Temp-tendency from sso + gravity-wave drag + Rayleigh friction
      ddt_temp_pconv  (:,:,:)  ,& !! Temp-tendency from convective prec
      ddt_u_turb      (:,:,:)  ,& !! ZonalW-tendency from turbulence
      ddt_u_gwd       (:,:,:)  ,& !! ZonalW-tendency from gravity wave drag
      ddt_u_raylfric  (:,:,:)  ,& !! ZonalW-tendency from artificial Rayleigh friction
      ddt_u_sso       (:,:,:)  ,& !! ZonalW-tendency from sso drag
      ddt_u_pconv     (:,:,:)  ,& !! ZonalW-tendency from convective prec
      ddt_v_turb      (:,:,:)  ,& !! MeridW-tendency from turbulence
      ddt_v_gwd       (:,:,:)  ,& !! MeridW-tendency from gravity wave drag
      ddt_v_raylfric  (:,:,:)  ,& !! MeridW-tendency from artificial Rayleigh friction
      ddt_v_sso       (:,:,:)  ,& !! MeridW-tendency from sso drag
      ddt_v_pconv     (:,:,:)  ,& !! MeridW-tendency from convective prec
      ddt_tracer_turb (:,:,:,:),& !! Hydromet-tendency from turbulence
      ddt_tracer_pconv(:,:,:,:),& !! Hydromet-tendency from convective prec
      ddt_tke         (:,:,:)     !! tendency for turbulent kinetic energy [m^2/s^3]

    TYPE(t_ptr_2d3d),ALLOCATABLE ::  &
      &  tracer_turb_ptr(:)    ,& !< pointer array: one pointer for each component
      &  tracer_conv_ptr(:)       !< pointer array: one pointer for each component
  END TYPE t_nwp_phy_tend

END MODULE mo_nwp_phy_types
