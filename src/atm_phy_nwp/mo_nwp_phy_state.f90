#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
MODULE mo_nwp_phy_state
!>
!!  !MODULE:  mo_nwp_phy_state\\
!!
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the ECHAM physics.
!!  Constructors and destructors for these data structures as well as
!!  initialization of fields are also defined here.
!!  This module should be an analogon to 'mo_hydro_state.f90'

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

! !USES:

USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: icc, success, max_char_length,      &
  &                               VINTP_METHOD_UV, VINTP_TYPE_P_OR_Z, &
&                                 VINTP_METHOD_LIN,VINTP_METHOD_QV
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: nqtendphy, iqv, iqc, iqi
USE mo_exception,           ONLY: message, finish !,message_text
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: n_dom
USE mo_icoham_sfc_indices,  ONLY: nsfc_type
USE mo_linked_list,         ONLY: t_var_list
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
USE mo_lnd_nwp_config,      ONLY: nsfc_subs
USE mo_var_list,            ONLY: default_var_list_settings, &
  &                               add_var, add_ref, new_var_list, delete_var_list, &
  &                               add_var_list_reference, create_vert_interp_metadata, &
  &                               groups
USE mo_nwp_parameters,      ONLY: t_phy_params
USE mo_cf_convention,       ONLY: t_cf_var
USE mo_grib2,               ONLY: t_grib2_var
USE mo_io_config,           ONLY: lflux_avg
USE mo_cdi_constants 

IMPLICIT NONE
PRIVATE

INCLUDE 'netcdf.inc'

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!public interface
!
! subroutines
PUBLIC :: construct_nwp_phy_state
PUBLIC :: destruct_nwp_phy_state
!
!types
PUBLIC :: t_nwp_phy_diag, t_nwp_phy_tend
#ifdef HAVE_F95
PUBLIC :: t_ptr_phy
#endif
!
!variables
PUBLIC :: prm_diag 
PUBLIC :: prm_nwp_tend
PUBLIC :: phy_params
PUBLIC :: prm_nwp_diag_list  !< variable lists
PUBLIC :: prm_nwp_tend_list  !< variable lists

!
!!data structure defining model states
!
!!diagnostic variables
!

  TYPE t_ptr_phy
    REAL(wp),POINTER :: p_3d(:,:,:) ! pointer to 3D (spatial) array
    REAL(wp),POINTER :: p_2d(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_phy


TYPE t_nwp_phy_diag

   TYPE(t_ptr_phy),ALLOCATABLE :: tot_ptr(:)  !< pointer array: one pointer for each tot var (grid+subgrid)
   TYPE(t_ptr_phy),ALLOCATABLE :: tci_ptr(:)  !< pointer array: total column-integrated values
   TYPE(t_ptr_phy),ALLOCATABLE :: tav_ptr(:)  !< pointer array: average of tci

   TYPE(t_ptr_phy),ALLOCATABLE :: cfm_ptr(:)  !< pointer array: average of cfm
   TYPE(t_ptr_phy),ALLOCATABLE :: cfh_ptr(:)  !< pointer array: average of cfh
   TYPE(t_ptr_phy),ALLOCATABLE :: z0m_ptr(:)  !< pointer array: average of z0m
   TYPE(t_ptr_phy),ALLOCATABLE :: tra_rate_ptr(:)   !< pointer array: tracer rates
   TYPE(t_ptr_phy),ALLOCATABLE :: swflxsfc_t_ptr(:) !< pointer array: shortwave net flux at surface
   TYPE(t_ptr_phy),ALLOCATABLE :: lwflxsfc_t_ptr(:) !< pointer array: longwave net flux at surface
   TYPE(t_ptr_phy),ALLOCATABLE :: tcm_t_ptr(:) !< pointer array: turbulent transfer coefficients for momentum
   TYPE(t_ptr_phy),ALLOCATABLE :: tch_t_ptr(:) !< pointer array: turbulent transfer coefficients for heat
   TYPE(t_ptr_phy),ALLOCATABLE :: tfm_t_ptr(:) !< pointer array: factor of laminar transfer of momentum
   TYPE(t_ptr_phy),ALLOCATABLE :: tfh_t_ptr(:) !< pointer array: factor of laminar transfer of heat 
   TYPE(t_ptr_phy),ALLOCATABLE :: tfv_t_ptr(:) !< pointer array: laminar reduction factor for evaporation
   TYPE(t_ptr_phy),ALLOCATABLE :: gz0_t_ptr(:) !< pointer array: roughness length * gravity
   TYPE(t_ptr_phy),ALLOCATABLE :: sai_t_ptr(:) !< pointer array: surface area index
   TYPE(t_ptr_phy),ALLOCATABLE :: tai_t_ptr(:) !< pointer array: transpiration area index
   TYPE(t_ptr_phy),ALLOCATABLE :: eai_t_ptr(:) !< pointer array: evaporation area index
   TYPE(t_ptr_phy),ALLOCATABLE :: t_2m_t_ptr(:)  !< pointer array: temperature at 2m
   TYPE(t_ptr_phy),ALLOCATABLE :: qv_2m_t_ptr(:) !< pointer array: specific water vapor content at 2m
   TYPE(t_ptr_phy),ALLOCATABLE :: rh_2m_t_ptr(:) !< pointer array: relative humidity at 2m
   TYPE(t_ptr_phy),ALLOCATABLE :: td_2m_t_ptr(:) !< pointer array: dew point at 2m
   TYPE(t_ptr_phy),ALLOCATABLE :: u_10m_t_ptr(:) !< pointer array: zonal wind at 2am
   TYPE(t_ptr_phy),ALLOCATABLE :: v_10m_t_ptr(:) !< pointer array: meridional wind at 2m
   TYPE(t_ptr_phy),ALLOCATABLE :: shfl_s_t_ptr(:) !< pointer array: surface sensible heat flux 
   TYPE(t_ptr_phy),ALLOCATABLE :: lhfl_s_t_ptr(:) !< pointer array: surface latent heat flux

  REAL(wp), POINTER ::  &
       &   tracer_rate(:,:,:) , & !> (nproma,nblks,4) precipitation rate of rain and snow [kg/m2/s]
       &   rain_gsp(:,:),       & !! accumulated grid-scale surface rain                  [kg/m2]
       &   snow_gsp(:,:),       & !! accumulated grid_scale surface snow                  [kg/m2]
       &   rain_con(:,:),       & !! accumulated convective surface rain                  [kg/m2]
       &   snow_con(:,:),       & !! accumulated convective surface snow                  [kg/m2]
       &   tot_prec(:,:),       & !! accumulated grid-scale plus convective surface total precipitation
       &   tot_prec_rate_avg(:,:),   & !! average since model start of 
                                  !! grid-scale plus convective surface 
                                  !! total precipitation rate
       &   con_prec_rate_avg(:,:),   & !! average since model start of 
                                  !! convective surface precipitation rate
       &   gsp_prec_rate_avg(:,:),   & !! average since model start of 
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
       &  rain_upd(:,:),        & !! total precipitation produced in updrafts (prain)
       &  shfl_s(:,:),          & !! sensible heat flux (surface) ( W/m2)
       &  shfl_s_t(:,:,:),      & !! sensible heat flux (surface) ( W/m2)
       &  lhfl_s(:,:),          & !! latent   heat flux (surface) ( W/m2)
       &  lhfl_s_t(:,:,:),      & !! latent   heat flux (surface) ( W/m2)
       &  qhfl_s(:,:),          & !!      moisture flux (surface) ( Kg/m2/s)
                                  !!      = evaporation rate at surface
       &  shfl_s_a(:,:),        & !! average or accumulated since model start of shfl_s ( W/m2)
       &  lhfl_s_a(:,:),        & !! average or accumulated since model start of lhfl_s ( W/m2)
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
       h_ice(:,:),          & !! ice thickness                                 (  m  )
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
       sai_t(:,:,:)     ,   & !! surface area index                            ( 1 )
       tai_t(:,:,:)     ,   & !! transpiration area index                      ( 1 )
       eai_t(:,:,:)     ,   & !! (evaporative) earth area index                ( 1 )
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
   TYPE(t_ptr_phy),ALLOCATABLE :: tracer_turb_ptr(:)  !< pointer array: one pointer for each component
   TYPE(t_ptr_phy),ALLOCATABLE :: tracer_conv_ptr(:)  !< pointer array: one pointer for each component
END TYPE t_nwp_phy_tend

!!--------------------------------------------------------------------------
!!                          STATE VARIABLES 
!!--------------------------------------------------------------------------
  TYPE(t_nwp_phy_diag), ALLOCATABLE, TARGET :: prm_diag(:) !< shape: (n_dom)
  TYPE(t_nwp_phy_tend), ALLOCATABLE :: prm_nwp_tend(:)     !< shape: (n_dom)
!-------------------------------------------------------------------------
  
!!--------------------------------------------------------------------------
!!                          VARIABLE LISTS
!!--------------------------------------------------------------------------
  TYPE(t_var_list),ALLOCATABLE :: prm_nwp_diag_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),ALLOCATABLE :: prm_nwp_tend_list(:)  !< shape: (n_dom)

!!-------------------------------------------------------------------------
!! Parameters of various physics parameterizations that have to be 
!! domain-dependent (computed during physics initialization phase)
!!-------------------------------------------------------------------------
  TYPE (t_phy_params), ALLOCATABLE :: phy_params(:)  !< shape: (n_dom)


CONTAINS

!-------------------------------------------------------------------------

SUBROUTINE construct_nwp_phy_state( p_patch)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)

CHARACTER(len=max_char_length) :: listname
INTEGER ::  jg,ist, nblks_c, nlev, nlevp1

!-------------------------------------------------------------------------

CALL message('mo_nwp_phy_state:construct_nwp_state', &
  'start to construct 3D state vector')


  ! Allocate pointer arrays prm_diag_nwp and prm_nwp_tend, 
  ! as well as the corresponding list arrays.

  ALLOCATE(prm_diag(n_dom), prm_nwp_diag_list(n_dom),STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of diagnostic physical array and list failed')
  ENDIF

  ALLOCATE(prm_nwp_tend(n_dom), prm_nwp_tend_list(n_dom), STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of tendency physical array and list failed')
  ENDIF

  DO jg = 1, n_dom

     !determine size of arrays
     nblks_c = p_patch(jg)%nblks_c
     
     ! number of vertical levels
     nlev   = p_patch(jg)%nlev
     nlevp1 = p_patch(jg)%nlevp1
     
     WRITE(listname,'(a,i2.2)') 'prm_diag_of_domain_',jg
     
     CALL new_nwp_phy_diag_list( jg, nlev, nlevp1, nblks_c,&
                               & TRIM(listname), prm_nwp_diag_list(jg), prm_diag(jg))
     !
     WRITE(listname,'(a,i2.2)') 'prm_tend_of_domain_',jg
     CALL new_nwp_phy_tend_list ( jg, nlev, nblks_c,&
                                & TRIM(listname), prm_nwp_tend_list(jg), prm_nwp_tend(jg))
  ENDDO


  ! Allocate variable of type t_phy_params containing domain-dependent parameters
  !
  ALLOCATE(phy_params(n_dom), STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of phy_params array failed')
  ENDIF

  
  CALL message('mo_nwp_phy_state:construct_nwp_state', &
    'construction of state vector finished')

END SUBROUTINE construct_nwp_phy_state

!
SUBROUTINE destruct_nwp_phy_state

  INTEGER :: jg, ist  !< grid level/domain index

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
  'start to destruct 3D state vector')

  DO jg = 1,n_dom
    CALL delete_var_list( prm_nwp_diag_list(jg) )
    CALL delete_var_list( prm_nwp_tend_list (jg) )
  ENDDO


  DEALLOCATE(prm_diag, prm_nwp_diag_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:destruct_nwp_phy_state', &
       &  'deallocation of NWP physics diagnostic array and list failed')
  ENDIF
 
  DEALLOCATE(prm_nwp_tend, prm_nwp_tend_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:destruct_nwp_phy_state', &
         &' deallocation of NWP physics tendencies array and list failed') 
  ENDIF

  DEALLOCATE(phy_params, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:destruct_nwp_phy_state', &
         &' deallocation of phy_params array failed') 
  ENDIF

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
    'destruction of 3D state vector finished')

END SUBROUTINE destruct_nwp_phy_state

     !
SUBROUTINE new_nwp_phy_diag_list( k_jg, klev, klevp1, kblks,   &
                     & listname, diag_list, diag)

    INTEGER,INTENT(IN) :: klev, klevp1, kblks,k_jg !< dimension sizes

    CHARACTER(len=*),INTENT(IN)     :: listname
    CHARACTER(len=max_char_length)  :: vname_prefix
    CHARACTER(LEN=1)                :: csfc

    TYPE(t_var_list)    ,INTENT(INOUT) :: diag_list
    TYPE(t_nwp_phy_diag),INTENT(INOUT) :: diag

    ! Local variables

    INTEGER :: n_updown = 7 !> number of up/downdrafts variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shapesfc(3), shape3dsubs(3)
    INTEGER :: shape3dkp1(3)
    INTEGER :: ientr,  kcloud
    INTEGER :: jsfc 
    CHARACTER(len=NF_MAX_NAME) :: long_name
    CHARACTER(len=21) :: name
    CHARACTER(len=3)  :: prefix
    CHARACTER(len=8)  :: meaning
    CHARACTER(len=10) :: varunits  ! variable units, depending on "lflux_avg"
 
    ientr = 16 ! bits "entropy" of horizontal slice

    shape2d    = (/nproma,            kblks            /)
    shape3d    = (/nproma, klev,      kblks            /)
    shape3dkp1 = (/nproma, klevp1,    kblks            /)
    shape3dsubs = (/nproma, kblks,    nsfc_subs        /)

    IF( atm_phy_nwp_config(k_jg)%inwp_turb == 2) THEN
      shapesfc   = (/nproma,          kblks, nsfc_type /)
    ENDIF

    ! Register a field list and apply default settings

    CALL new_var_list( diag_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( diag_list,                 &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    !-------------------
    ! Clouds and precip
    !------------------
    ! 2D and 3D variables

    kcloud= 4


    ! &      diag%tracer_rate(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tracer_rate','kg m-2 ','precipitation rate of rain and snow')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tracer_rate', diag%tracer_rate,          &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,cf_desc, grib2_desc, &
          &                                   ldims=(/nproma,kblks,4/),&
          &         lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ALLOCATE( diag%tra_rate_ptr(kcloud))
    vname_prefix='tra_rate_'

           !qv
        CALL add_ref( diag_list, 'tracer_rate',                                     &
                    & TRIM(vname_prefix)//'qv', diag%tra_rate_ptr(1)%p_2d,          &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qv', '','unknown'),             &
                    & t_grib2_var(0, 1, 201, ientr, GRID_REFERENCE, GRID_CELL),     &
                    & ldims=shape2d)
           !qc
        CALL add_ref( diag_list, 'tracer_rate',                                     &
                    & TRIM(vname_prefix)//'qc', diag%tra_rate_ptr(2)%p_2d,          &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qc', '','unknown'),             &
                    & t_grib2_var(0, 1, 202, ientr, GRID_REFERENCE, GRID_CELL),     &
                    & ldims=shape2d)
!          !qi
!       CALL add_ref( diag_list, 'tracer_rate',                                     &
!                   & TRIM(vname_prefix)//'qi', diag%tra_rate_ptr(3)%p_2d,          &
!                   & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
!                   & t_cf_var(TRIM(vname_prefix)//'qi', '','unknown'),             &
!                   & t_grib2_var(0, 1, 203, ientr, GRID_REFERENCE, GRID_CELL),     &
!                   & ldims=shape2d)
           !qr
        CALL add_ref( diag_list, 'tracer_rate',                                     &
                    & TRIM(vname_prefix)//'qr', diag%tra_rate_ptr(3)%p_2d,          &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qr', '','precipitation rate'),  &
                    & t_grib2_var(0, 1, 24, ientr, GRID_REFERENCE, GRID_CELL),      &
                    & ldims=shape2d)
           !qs
        CALL add_ref( diag_list, 'tracer_rate',                                     &
                    & TRIM(vname_prefix)//'qs', diag%tra_rate_ptr(4)%p_2d,          &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qs', '','snowfall rate'),       &
                    & t_grib2_var(0, 1, 25, ientr, GRID_REFERENCE, GRID_CELL),      &
                    & ldims=shape2d)



    ! &      diag%rain_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp ', 'kg m-2 ', 'gridscale rain ')
    grib2_desc = t_grib2_var(0, 1, 77, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp', diag%rain_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
                & ldims=shape2d, in_group=groups("precip_vars") )

    ! &      diag%snow_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp', 'kg m-2 ', 'gridscale snow')
    grib2_desc = t_grib2_var(0, 1, 15, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp', diag%snow_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
                & ldims=shape2d, istatproc=TSTEP_ACCUM,                      &
                & in_group=groups("precip_vars") )

    ! &      diag%rain_con(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con', 'kg m-2 ', 'convective rain')
    grib2_desc = t_grib2_var(0, 1, 76, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_con', diag%rain_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, in_group=groups("precip_vars") )

    ! &      diag%snow_con(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con', 'kg m-2', 'convective snow')
    grib2_desc = t_grib2_var(0, 1, 14, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_con', diag%snow_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, istatproc=TSTEP_ACCUM,                       &
                & in_group=groups("precip_vars") )

    ! &      diag%tot_prec(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec', '', 'total precip')
    grib2_desc = t_grib2_var(0, 1, 52, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec', diag%tot_prec,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, istatproc=TSTEP_ACCUM,                       &
                & in_group=groups("precip_vars") )

    ! &      diag%tot_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_rate_avg', '', 'total precip, time average')
    grib2_desc = t_grib2_var(0, 1, 52, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec_rate_avg', diag%tot_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%con_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('con_prec_rate_avg', '', 'convective precip, time average')
    grib2_desc = t_grib2_var(0, 1, 10, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_prec_rate_avg', diag%con_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d,  lrestart=.FALSE.,                           &
                & in_group=groups("additional_precip_vars") )

    ! &      diag%gsp_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('gsp_prec_rate_avg', '', 'gridscale precip, time average')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'gsp_prec_rate_avg', diag%gsp_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars") )

    ! &      diag%cape(nproma,nblks_c)
    cf_desc    = t_cf_var('cape', 'J kg-1 ', 'conv avail pot energy')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cape', diag%cape,                               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars") )

    ! &      diag%con_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('con_gust', 'm s-1 ', 'convective gusts')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_gust', diag%con_gust,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, lrestart=.FALSE. )
   
    ! &      diag%rain_upd(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_upd', 'unit ', 'rain in updroughts')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_upd', diag%rain_upd,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%con_udd(nproma,nlev,nblks,8)
    cf_desc    = t_cf_var('con_udd', 'unit ', 'convective up/downdraft fields')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_udd', diag%con_udd,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                & ldims=(/nproma,klev,kblks,n_updown/),                       &
                & lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%mbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mbas_con', '', 'cloud base level index')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mbas_con', diag%mbas_con,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%mtop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mtop_con', '', 'cloud top level index')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mtop_con', diag%mtop_con,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%locum(nproma,nblks_c)
    cf_desc    = t_cf_var('locum', '', 'convective  activity indicator')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'locum', diag%locum,                              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%ktype(nproma,nblks_c)
    cf_desc    = t_cf_var('ktype', '', 'convective up/downdraft fields')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ktype', diag%ktype,                              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,              &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

   ! &      diag%tot_cld_vi(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tot_cld_vi', 'unit ','vertical integr total cloud variables')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi', diag%tot_cld_vi,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                &                                 ldims=(/nproma,kblks,4/),   &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the seperate variables belonging to the container tot_cld_vi
    ALLOCATE( diag%tci_ptr(kcloud))
    vname_prefix='T'
       
    !QV
    CALL add_ref( diag_list, 'tot_cld_vi',                        &
      & TRIM(vname_prefix)//'QV', diag%tci_ptr(1)%p_2d,           &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                    &
      & t_cf_var(TRIM(vname_prefix)//'QV', 'kg m**-2','column_integrated_water_vapour'),&
      & t_grib2_var( 0, 1, 64, ientr, GRID_REFERENCE, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !qc
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & TRIM(vname_prefix)//'QC', diag%tci_ptr(2)%p_2d,            &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
      & t_cf_var(TRIM(vname_prefix)//'QC', 'kg m**-2',             &
      & 'total_column-integrated_cloud_water'),                    &
      & t_grib2_var( 0, 1, 69, ientr, GRID_REFERENCE, GRID_CELL),  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !qi
    CALL add_ref( diag_list, 'tot_cld_vi',&
      & TRIM(vname_prefix)//'QI', diag%tci_ptr(3)%p_2d,            &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
      & t_cf_var(TRIM(vname_prefix)//'qi', 'kg m**-2',             &
      & 'total_column-integrated_cloud_water'),                    &
      & t_grib2_var(0, 1, 70, ientr, GRID_REFERENCE, GRID_CELL),   &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))
    
    !CC
    CALL add_ref( diag_list, 'tot_cld_vi',            &
      & TRIM(vname_prefix)//'CC', diag%tci_ptr(4)%p_2d,                  &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                           &
      & t_cf_var(TRIM(vname_prefix)//'CC', '%','total_column_integrated_cloud_cover'), &
      & t_grib2_var(0, 6, 1, ientr, GRID_REFERENCE, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars") )


   ! &      diag%tot_cld_vi_avg(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tot_cld_vi_avg', 'unit ','vertical integr total cloud variables')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi_avg', diag%tot_cld_vi_avg,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                &                                 ldims=(/nproma,kblks,4/)  , &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! fill the seperate variables belonging to the container tot_cld_vi_avg
    ALLOCATE( diag%tav_ptr(kcloud))
    vname_prefix='avg_'

           !CC
    CALL add_ref( diag_list, 'tot_cld_vi_avg',                                 &
                & TRIM(vname_prefix)//'cc', diag%tav_ptr(1)%p_2d,              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var(TRIM(vname_prefix)//'cc', '','tci_cloud_cover_avg'), &
                & t_grib2_var(0, 6, 1, ientr, GRID_REFERENCE, GRID_CELL),      &
                & ldims=shape2d, lrestart=.FALSE. )

       !QV
    CALL add_ref( diag_list, 'tot_cld_vi_avg',               &
                & TRIM(vname_prefix)//'qv', diag%tav_ptr(2)%p_2d,              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var(TRIM(vname_prefix)//'qv', '','tci_specific_humidity_avg'), &
                & t_grib2_var( 0, 1, 64, ientr, GRID_REFERENCE, GRID_CELL),    &
                & ldims=shape2d, lrestart=.FALSE. )

       !qc
    CALL add_ref( diag_list, 'tot_cld_vi_avg',      &
                & TRIM(vname_prefix)//'qc', diag%tav_ptr(3)%p_2d,              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var(TRIM(vname_prefix)//'qc', '',                       &
                & 'tci_specific_cloud_water_content_avg'),                     &
                & t_grib2_var(0, 6, 18, ientr, GRID_REFERENCE, GRID_CELL),     &
                & ldims=shape2d, lrestart=.FALSE. )

       !qi
    CALL add_ref( diag_list, 'tot_cld_vi_avg',          &
                & TRIM(vname_prefix)//'qi', diag%tav_ptr(4)%p_2d,              & 
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                & t_cf_var(TRIM(vname_prefix)//'qi', '',                       &
                & 'tci_specific_cloud_ice_content_avg'),                       &
                & t_grib2_var(0, 6, 19, ientr, GRID_REFERENCE, GRID_CELL),     &
                & ldims=shape2d, lrestart=.FALSE. )

   ! &      diag%tot_cld(nproma,nlev,nblks_c,4)
    cf_desc    = t_cf_var('tot_cld', ' ','total cloud variables (cc,qv,qc,qi)')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld', diag%tot_cld,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=(/nproma,klev,kblks,4/) ,                               &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                & initval_r=0.0_wp )

    ALLOCATE( diag%tot_ptr(kcloud))
    vname_prefix='tot_'

           !QV
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'qv', diag%tot_ptr(iqv)%p_3d,                &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var(TRIM(vname_prefix)//'qv', '','total_specific_humidity'),&
                    & t_grib2_var( 0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL),         &
                    & ldims=shape3d,                                                   &
                    & vert_interp=create_vert_interp_metadata(                         &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                    &
                    &             vert_intp_method=VINTP_METHOD_QV,                    &
                    &             l_satlimit=.FALSE.,                                  & 
                    &             lower_limit=2.5e-6_wp, l_restore_pbldev=.FALSE. ),   &
                    & in_group=groups("cloud_diag") )

           !QC
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'qc', diag%tot_ptr(iqc)%p_3d,                &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var(TRIM(vname_prefix)//'qc', '',                           &
                    & 'total_specific_cloud_water_content'),                           &
                    & t_grib2_var(0, 1, 83, ientr, GRID_REFERENCE, GRID_CELL),         &
                    & ldims=shape3d,                                                   &
                    & vert_interp=create_vert_interp_metadata(                         &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                    &
                    &             vert_intp_method=VINTP_METHOD_LIN,                   &
                    &             l_loglin=.FALSE.,                                    &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,               &
                    &             lower_limit=0._wp ),                                 &
                    & in_group=groups("cloud_diag") )

           !QI
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'qi', diag%tot_ptr(iqi)%p_3d,                &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var(TRIM(vname_prefix)//'qi', '',                           &
                    & 'total_specific_cloud_ice_content'),                             &
                    & t_grib2_var(0, 1, 84, ientr, GRID_REFERENCE, GRID_CELL),         &
                    & ldims=shape3d,                                                   &
                    & vert_interp=create_vert_interp_metadata(                         &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                    &
                    &             vert_intp_method=VINTP_METHOD_LIN,                   &
                    &             l_loglin=.FALSE.,                                    &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,               &
                    &             lower_limit=0._wp ),                                 &
                    & in_group=groups("cloud_diag") )

           !CC
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'cc', diag%tot_ptr(icc)%p_3d,                &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var(TRIM(vname_prefix)//'cc', '','total_cloud_cover'),      &
                    & t_grib2_var(0, 6, 22, ientr, GRID_REFERENCE, GRID_CELL),         &
                    & ldims=shape3d,                                                   &
                    & vert_interp=create_vert_interp_metadata(                         &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                    &
                    &             vert_intp_method=VINTP_METHOD_LIN,                   &
                    &             l_loglin=.FALSE.,                                    &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,               &
                    &             lower_limit=0._wp ),                                 &
                    & in_group=groups("cloud_diag") )


        ! &      diag%acdnc(nproma,nlsacev,nblks_c)
        cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'acdnc', diag%acdnc,                               &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,              &
          & ldims=shape3d, lrestart=.FALSE. )
        diag%acdnc(:,:,:)=220.0e6_wp

    !------------------
    ! Radiation
    !------------------
    ! 2D variables

        !        diag%albvisdif    (nproma,       nblks),          &
        cf_desc    = t_cf_var('albvisdif', '', '')
        grib2_desc = t_grib2_var(192, 128, 243, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'albvisdif', diag%albvisdif,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d, lrestart=.FALSE.,                                    &
          & in_group=groups("rad_vars") )


        ! These variables only make sense, if the land-surface scheme is switched on.
        IF ( atm_phy_nwp_config(k_jg)%inwp_surface == 1 ) THEN

          !        diag%albvisdif_t    (nproma, nblks, nsfc_subs),          &
          cf_desc    = t_cf_var('albvisdif_t', '', '')
          grib2_desc = t_grib2_var(192, 128, 243, ientr, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'albvisdif_t', diag%albvisdif_t,               &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
            & ldims=(/nproma,kblks,nsfc_subs/), lcontainer=.TRUE., lrestart=.FALSE. )


          ! &      diag%swflxsfc_t(nproma,nblks_c,nsfc_subs)
          cf_desc    = t_cf_var('SOB_S_T', 'W m-2', 'tile-based shortwave net flux at surface')
          grib2_desc = t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'SOB_S_T', diag%swflxsfc_t,                       &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,            &
            & ldims=(/nproma,kblks,nsfc_subs/), lcontainer=.TRUE., lrestart=.FALSE., &
            & loutput=.FALSE.)

          ! fill the seperate variables belonging to the container swflxsfc_t
          ALLOCATE(diag%swflxsfc_t_ptr(nsfc_subs))
          DO jsfc = 1,nsfc_subs
            WRITE(csfc,'(i1)') jsfc 
            CALL add_ref( diag_list, 'SOB_S_T',                                &
               & 'SOB_S_T_'//TRIM(ADJUSTL(csfc)),                              &
               & diag%swflxsfc_t_ptr(jsfc)%p_2d,                               &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
               & t_cf_var('swflxsfc_t_'//TRIM(csfc), '', ''),                  &
               & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),       &
               & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,               &
               & in_group=groups("rad_vars"))
          ENDDO

          ! &      diag%lwflxsfc_t(nproma,nblks_c,nsfc_subs)
          cf_desc    = t_cf_var('THB_S_T', 'W m-2', 'tile_based longwave net flux at surface')
          grib2_desc = t_grib2_var(0, 5, 0, ientr, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'THB_S_T', diag%lwflxsfc_t,                        &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,             &
            & ldims=(/nproma,kblks,nsfc_subs/),lcontainer=.TRUE., lrestart=.FALSE.,   &
            & loutput=.FALSE.)

          ! fill the seperate variables belonging to the container lwflxsfc_t
          ALLOCATE(diag%lwflxsfc_t_ptr(nsfc_subs))
          DO jsfc = 1,nsfc_subs
            WRITE(csfc,'(i1)') jsfc 
            CALL add_ref( diag_list, 'THB_S_T',                                &
               & 'THB_S_T_'//TRIM(ADJUSTL(csfc)),                              &
               & diag%lwflxsfc_t_ptr(jsfc)%p_2d,                               &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
               & t_cf_var('lwflxsfc_t_'//TRIM(csfc), '', ''),                  &
               & t_grib2_var(0, 5, 0, ientr, GRID_REFERENCE, GRID_CELL),       &
               & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,               &
               & in_group=groups("rad_vars"))
          ENDDO
        ENDIF


        !        diag%cosmu0    (nproma,       nblks),          &
        cf_desc    = t_cf_var('cosmu0', '', '')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'cosmu0', diag%cosmu0,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%tsfctrad(nproma,nblks_c)
        cf_desc    = t_cf_var('tsfctrad', 'K', 'surface temperature at trad')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tsfctrad', diag%tsfctrad,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &       diag% flxdwswtoa(nproma,       nblks),          &
        cf_desc    = t_cf_var('flxdwswtoa', 'W m-2', 'downward shortwave flux at TOA')
        grib2_desc = t_grib2_var(0, 4, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'flxdwswtoa', diag%flxdwswtoa,                 &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d, lrestart=.FALSE.)


        ! &      diag%swflxsfc(nproma,nblks_c)
        cf_desc    = t_cf_var('SOB_S', 'W m-2', 'shortwave net flux at surface')
        grib2_desc = t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'SOB_S', diag%swflxsfc,                        &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d)



        ! &      diag%swflxtoa(nproma,nblks_c)
        cf_desc    = t_cf_var('SOB_T', 'W m-2', 'shortwave net flux at TOA')
        grib2_desc = t_grib2_var(0, 4, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'SOB_T', diag%swflxtoa,                        &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d, lrestart=.FALSE.,                                    &
          & in_group=groups("rad_vars"))

        ! &      diag%lwflxsfc(nproma,nblks_c)
        cf_desc    = t_cf_var('THB_S', 'W m-2', 'longwave net flux at surface')
        grib2_desc = t_grib2_var(0, 5, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'THB_S', diag%lwflxsfc,                        &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                      &
          & in_group=groups("rad_vars"))


        IF (lflux_avg ) THEN
            prefix = "A"
            meaning = "mean"
            varunits= "W/m**2"
        ELSE
            prefix = "ACC"
            meaning = "acc." 
            varunits= "J/m**2"    
        END IF
        WRITE(name,'(A,A5)') TRIM(prefix),"THB_T"
        WRITE(long_name,'(A26,A4,A18)') "longwave  net flux at TOA ", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                          TRIM(long_name))
        grib2_desc = t_grib2_var(0, 5, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%lwflxtoa_a,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d, istatproc=TSTEP_ACCUM)


        ! &      diag%swflxtoa_a(nproma,nblks_c)
        WRITE(name,'(A,A5)') TRIM(prefix),"SOB_T"
        WRITE(long_name,'(A26,A4,A18)') "shortwave net flux at TOA ", meaning,  &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                         TRIM(long_name))
        grib2_desc = t_grib2_var(0, 4, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name) , diag%swflxtoa_a,                  &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

        
        ! &      diag%lwflxsfc_a(nproma,nblks_c)
        WRITE(name,'(A,A5)') TRIM(prefix),"THB_S"
        WRITE(long_name,'(A30,A4,A18)') "longwave  net flux at surface ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                      TRIM(long_name))
        grib2_desc = t_grib2_var(0, 5, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%lwflxsfc_a,                 &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%swflxsfc_a(nproma,nblks_c)
        WRITE(name,'(A,A5)') TRIM(prefix),"SOB_S"
        WRITE(long_name,'(A30,A4,A18)') "shortwave net flux at surface ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                    TRIM(long_name) )
        grib2_desc = t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%swflxsfc_a,                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%vio3(nproma,nblks_c)
        cf_desc    = t_cf_var('vio3', '', 'vertically integrated ozone amount')
        grib2_desc = t_grib2_var(0, 14, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'vio3', diag%vio3,                           &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%hmo3(nproma,nblks_c)
        cf_desc    = t_cf_var('hmo3', 'Pa', 'height of O3 maximum (Pa)')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'hmo3', diag%hmo3,                           &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )



   ! for old aerosol climatology from COSMO (to be used with inwp_radiation==2)

        ! &      diag%aersea(nproma,nblks_c)
        cf_desc    = t_cf_var('aersea', '', '')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aersea', diag%aersea,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. ) 

        
        ! &      diag%aerlan(nproma,nblks_c)
        cf_desc    = t_cf_var('aerlan', '', '')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aerlan', diag%aerlan,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. ) 

        
        ! &      diag%aerurb(nproma,nblks_c)
        cf_desc    = t_cf_var('aerurb', '', '')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aerurb', diag%aerurb,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. ) 

        
        ! &      diag%aerdes(nproma,nblks_c)
        cf_desc    = t_cf_var('aerdes', '', '')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aerdes', diag%aerdes,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. ) 


        !------------------
        !Radiation 3D variables

        ! &      diag%lwflxclr(nproma,nlev,nblks_c)
        cf_desc    = t_cf_var('lwflxclr', 'W m-2', 'longwave clear-sky net flux')
        grib2_desc = t_grib2_var(0, 5, 6, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxclr', diag%lwflxclr,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
          & ldims=shape3dkp1, lrestart=.FALSE. )

        ! &      diag% lwflxall(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('lwflxall', 'W m-2 ', 'longwave net flux')
        grib2_desc = t_grib2_var(0, 5, 5, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxall', diag%lwflxall,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
          & ldims=shape3dkp1 )

        ! &      diag% trsolclr(nproma,nlev,nblks_c)
        cf_desc    = t_cf_var('trsolclr', '', 'shortwave clear-sky net tranmissivity')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsolclr', diag%trsolclr,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
          & ldims=shape3dkp1, lrestart=.FALSE. )

        ! &      diag%  trsolall(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('trsolall', '', 'shortwave net tranmissivity')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsolall', diag%trsolall,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1)


        !------------------
        !Turbulence 2D variables
        
        ! &      diag%shfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('SHFL_S', 'W m-2 ', 'surface sensible heat flux')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'SHFL_S', diag%shfl_s,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars") )


        WRITE(name,'(A,A6)') TRIM(prefix),"SHFL_S"
        WRITE(long_name,'(A27,A4,A18)') "surface sensible heat flux ", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name),  TRIM(varunits), TRIM(long_name))
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%shfl_s_a,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d )


        ! &      diag%lhfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('LHFL_S', 'W m-2 ', 'surface latent heat flux')
        grib2_desc = t_grib2_var(0, 0, 10, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'LHFL_S', diag%lhfl_s,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars"))

                    
        WRITE(name,'(A,A6)') TRIM(prefix),"LHFL_S"
        WRITE(long_name,'(A27,A4,A18)') "surface latent   heat flux ", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name))
        grib2_desc = t_grib2_var(0, 0, 10, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%lhfl_s_a,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d )


        ! &      diag%qhfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('qhfl_s', 'Kg m-2 s-1', 'surface moisture flux')
        grib2_desc = t_grib2_var(2, 0, 6, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qhfl_s', diag%qhfl_s,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars"))


        cf_desc    = t_cf_var('qhfl_s_avg', 'Kg m-2 s-1', 'surface moisture flux')
        grib2_desc = t_grib2_var(2, 0, 6, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qhfl_s_avg', diag%qhfl_s_avg,               &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d )


        ! &      diag%tcm(nproma,nblks_c)
        cf_desc    = t_cf_var('tcm', ' ','turbulent transfer coefficients for momentum')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tcm', diag%tcm,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE.,                                  &
          & in_group=groups("pbl_vars") )


        ! &      diag%tch(nproma,nblks_c)
        cf_desc    = t_cf_var('tch', ' ','turbulent transfer coefficients for heat')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tch', diag%tch,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE.,                                  &
          & in_group=groups("pbl_vars") )

        
        ! &      diag%tfm(nproma,nblks_c)
        cf_desc    = t_cf_var('tfm', ' ','factor of laminar transfer of momentum')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfm', diag%tfm,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%tfh(nproma,nblks_c)
        cf_desc    = t_cf_var('tfh', ' ','factor of laminar transfer of scalars')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfh', diag%tfh,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%tfv(nproma,nblks_c)
        cf_desc    = t_cf_var('tfv', ' ','laminar reduction factor for evaporation')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfv', diag%tfv,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%gz0(nproma,nblks_c)
        cf_desc    = t_cf_var('gz0', 'm2 s-2 ','roughness length times gravity')
        grib2_desc = t_grib2_var(2, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'gz0', diag%gz0,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d )
        diag%gz0(:,:)=0.01_wp


        ! &      diag%sai(nproma,nblks_c)
        cf_desc    = t_cf_var('sai', ' ','surface area index')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sai', diag%sai,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d)


        ! &      diag%tai(nproma,nblks_c)
        cf_desc    = t_cf_var('tai', ' ','transpiration area index')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tai', diag%tai,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%eai(nproma,nblks_c)
        cf_desc    = t_cf_var('eai', ' ','(evaporative) earth area index')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'eai', diag%eai,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%t_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('T_2M', 'K ','temperature in 2m')
        grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'T_2M', diag%t_2m,                           &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )


        ! &      diag%t_2m_s6avg(nproma,nblks_c)
        cf_desc    = t_cf_var('T_2M_s6avg', 'K ','temperature in 2m')
        grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'T_2M_s6avg', diag%t_2m_s6avg,               &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%qv_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('QV_2M', 'kg kg-1 ','specific water vapor content in 2m')
        grib2_desc = t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'QV_2M', diag%qv_2m,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )


        ! &      diag%qv_2m_s6avg(nproma,nblks_c)
        cf_desc    = t_cf_var('QV_2M_s6avg', 'kg kg-1 ','specific water vapor content in 2m')
        grib2_desc = t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'QV_2M_s6avg', diag%qv_2m_s6avg,             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%rh_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('rh_2m', '%','relative humidity in 2m')
        grib2_desc = t_grib2_var(0, 1, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'rh_2m', diag%rh_2m,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%td_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('TD_2M', 'K ','dew-point in 2m')
        grib2_desc = t_grib2_var(0, 0, 6, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TD_2M', diag%td_2m,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )


        ! &      diag%u_10m(nproma,nblks_c)
        cf_desc    = t_cf_var('U_10M', 'm s-1 ','zonal wind in 10m')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'U_10M', diag%u_10m,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )

        ! &      diag%u_10m_s6avg(nproma,nblks_c)
        cf_desc    = t_cf_var('U_10M_s6avg', 'm s-1 ','zonal wind in 10m')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'U_10M_s6avg', diag%u_10m_s6avg,             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%v_10m(nproma,nblks_c)
        cf_desc    = t_cf_var('V_10M', 'm s-1 ','meridional wind in 10m')
        grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'V_10M', diag%v_10m,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )


        ! &      diag%v_10m_s6avg(nproma,nblks_c)
        cf_desc    = t_cf_var('V_10M_s6avg', 'm s-1 ','meridional wind in 10m')
        grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'V_10M_s6avg', diag%v_10m_s6avg,             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%h_ice(nproma,nblks_c)
        cf_desc    = t_cf_var('h_ice', 'm','ice thickness')
        grib2_desc = t_grib2_var(10, 2, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'h_ice', diag%h_ice,                         &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

!tiled quantities
        ! &      diag%shfl_s_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('SHFL_S_T', 'W m-2 ', 'tile-based surface sensible heat flux')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'SHFL_S_T', diag%shfl_s_t,                                &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container shfl_s_t
        ALLOCATE(diag%shfl_s_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'SHFL_S_T',                            &
             & 'SHFL_S_T_'//TRIM(ADJUSTL(csfc)),                          &
             & diag%shfl_s_t_ptr(jsfc)%p_2d,                              &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('shfl_s_t_'//TRIM(csfc), '', ''),                 &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%lhfl_s_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('LHFL_S_T', 'W m-2 ', 'tile-based surface latent heat flux')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'LHFL_S_T', diag%lhfl_s_t,                                &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container lhfl_s_t
        ALLOCATE(diag%lhfl_s_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'LHFL_S_T',                            &
             & 'LHFL_S_T_'//TRIM(ADJUSTL(csfc)),                          &
             & diag%lhfl_s_t_ptr(jsfc)%p_2d,                              &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     & 
             & t_cf_var('lhfl_s_t_'//TRIM(csfc), '', ''),                 &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    & 
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tcm_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TCM_T', ' ', &
          & 'tile-based turbulent transfer coefficients for momentum')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TCM_T', diag%tcm_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container tcm_t
        ALLOCATE(diag%tcm_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TCM_T',                               &
             & 'TCM_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%tcm_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     & 
             & t_cf_var('tcm_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    & 
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tch_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TCH_T', ' ', 'tile-based turbulent transfer coefficients for heat')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TCH_T', diag%tch_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)                       
          
        ! fill the separate variables belonging to the container tch_t
        ALLOCATE(diag%tch_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TCH_T',                               &
             & 'TCH_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%tch_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('tch_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tfm_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TFM_T', ' ', 'tile-based factor of laminar transfer of momentum')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TFM_T', diag%tfm_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container tfm_t
        ALLOCATE(diag%tfm_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TFM_T',                               &
             & 'TFM_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%tfm_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('tfm_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tfh_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TFH_T', ' ', 'tile-based factor of laminar transfer of heat')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TFH_T', diag%tfh_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container tfh_t
        ALLOCATE(diag%tfh_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TFH_T',                               &
             & 'TFH_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%tfh_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('tfh_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tfv_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TFV_T', ' ', 'tile-based laminar reduction factor for evaporation')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TFV_T', diag%tfv_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container tfv_t
        ALLOCATE(diag%tfv_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TFV_T',                               &
             & 'TFV_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%tfv_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('tfv_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%gz0_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('GZ0_T', 'm2 s-2 ', 'tile-based roughness length times gravity')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'GZ0_T', diag%gz0_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container gz0_t
        ALLOCATE(diag%gz0_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'GZ0_T',                               &
             & 'GZ0_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%gz0_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('gz0_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%sai_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('SAI_T', ' ', 'tile-based surface area index')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'SAI_T', diag%sai_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container sai_t
        ALLOCATE(diag%sai_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'SAI_T',                               &
             & 'SAI_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%sai_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('sai_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tai_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TAI_T', ' ', 'tile-based transpiration area index')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TAI_T', diag%tai_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container tai_t
        ALLOCATE(diag%tai_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TAI_T',                               &
             & 'TAI_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%tai_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('tai_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%eai_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('EAI_T', ' ', 'tile-based (evaporative) earth area index')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'EAI_T', diag%eai_t,                                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container eai_t
        ALLOCATE(diag%eai_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'EAI_T',                               &
             & 'EAI_T_'//TRIM(ADJUSTL(csfc)),                             &
             & diag%eai_t_ptr(jsfc)%p_2d,                                 &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('eai_t_'//TRIM(csfc), '', ''),                    &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%t_2m_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('T_2M_T', 'K ', 'tile-based temperature in 2m')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'T_2M_T', diag%t_2m_t,                                    &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container t_2m_t
        ALLOCATE(diag%t_2m_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'T_2M_T',                              &
             & 'T_2M_T_'//TRIM(ADJUSTL(csfc)),                            &
             & diag%t_2m_t_ptr(jsfc)%p_2d,                                &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('t_2m_t_'//TRIM(csfc), '', ''),                   &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%qv_2m_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('QV_2M_T', 'kg kg-1 ', 'tile-based water vapor content in 2m')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'QV_2M_T', diag%qv_2m_t,                                  &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
        ! fill the separate variables belonging to the container qv_2m_t
        ALLOCATE(diag%qv_2m_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'QV_2M_T',                             &
             & 'QV_2M_T_'//TRIM(ADJUSTL(csfc)),                           &
             & diag%qv_2m_t_ptr(jsfc)%p_2d,                               &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('qv_2m_t_'//TRIM(csfc), '', ''),                  &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%td_2m_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('TD_2M_T', 'K ', 'tile-based dew point in 2m')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'TD_2M_T', diag%td_2m_t,                                  &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
        ! fill the separate variables belonging to the container td_2m_t
        ALLOCATE(diag%td_2m_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'TD_2M_T',                             &
             & 'TD_2M_T_'//TRIM(ADJUSTL(csfc)),                           &
             & diag%td_2m_t_ptr(jsfc)%p_2d,                               &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('td_2m_t_'//TRIM(csfc), '', ''),                  &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%rh_2m_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('RH_2M_T', '% ', 'tile-based relative humidity in 2m')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'RH_2M_T', diag%rh_2m_t,                                  &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
        ! fill the separate variables belonging to the container rh_2m_t
        ALLOCATE(diag%rh_2m_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'RH_2M_T',                             &
             & 'RH_2M_T_'//TRIM(ADJUSTL(csfc)),                           &
             & diag%rh_2m_t_ptr(jsfc)%p_2d,                               &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('rh_2m_t_'//TRIM(csfc), '', ''),                  &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%u_10m_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('U_10M_T', 'm s-1 ', 'tile-based zonal wind in 2m')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'U_10M_T', diag%u_10m_t,                                  &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      
        ! fill the separate variables belonging to the container u_10m_t
        ALLOCATE(diag%u_10m_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'U_10M_T',                             &
             & 'U_10M_T_'//TRIM(ADJUSTL(csfc)),                           &
             & diag%u_10m_t_ptr(jsfc)%p_2d,                               &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('u_10m_t_'//TRIM(csfc), '', ''),                  &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%v_10m_t(nproma,nblks_c,nsfc_subs)
        cf_desc    = t_cf_var('V_10M_T', 'm s-1 ', 'tile-based meridional wind in 2m')
        grib2_desc = t_grib2_var(0, 0, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'V_10M_T', diag%v_10m_t,                                  &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs, &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

        ! fill the separate variables belonging to the container v_10m_t
        ALLOCATE(diag%v_10m_t_ptr(nsfc_subs))
        DO jsfc = 1,nsfc_subs
          WRITE(csfc,'(i1)') jsfc
          CALL add_ref( diag_list, 'V_10M_T',                             &
             & 'V_10M_T_'//TRIM(ADJUSTL(csfc)),                           &
             & diag%v_10m_t_ptr(jsfc)%p_2d,                               &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
             & t_cf_var('v_10m_t_'//TRIM(csfc), '', ''),                  &
             & t_grib2_var(0, 4, 0, ientr, GRID_REFERENCE, GRID_CELL),    &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO
  !
  ! vdiff
  !

    IF( atm_phy_nwp_config(k_jg)%inwp_turb == 2) THEN

      ! &      diag%cfm_tile(nproma,nblks_c)
      cf_desc    = t_cf_var('cfm_tile','',&
         & 'turbulent exchange coefficient of momentum at surface')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'cfm_tile', diag%cfm_tile,                     &
        & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shapesfc ,&
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.           )

      ALLOCATE(diag%cfm_ptr(nsfc_type))
      DO jsfc = 1,nsfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'cfm_tile',                                       &
                    & 'cfm_tile_'//csfc, diag%cfm_ptr(jsfc)%p_2d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_momentum_'//csfc, '', ''),       &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d, lrestart=.FALSE. )
      END DO

      ! &      diag%cfh_tile(nproma,nblks_c)
      cf_desc    = t_cf_var('cfh_tile', '','turbulent exchange coefficient of heat at surface')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'cfh_tile', diag%cfh_tile,                             &
        & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.           )

       ALLOCATE(diag%cfh_ptr(nsfc_type))
       DO jsfc = 1,nsfc_type
         WRITE(csfc,'(i1)') jsfc 
         CALL add_ref( diag_list, 'cfh_tile',                                      &
                    & 'cfh_tile_'//csfc, diag%cfh_ptr(jsfc)%p_2d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_heat_'//csfc, '', ''),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d, lrestart=.FALSE. )
       END DO

       ! &      diag%ghpbl(nproma,nblks_c)
       cf_desc    = t_cf_var('gh_pbl','','turbulent exchange coefficient of momentum at surface')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'ghpbl', diag%ghpbl,                          &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

       ! &      diag%z0m(nproma,nblks_c)
       cf_desc    = t_cf_var('z0m', '','geopotential of the top of the atmospheric boundary layer')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'z0m', diag%z0m,                              &
         & GRID_UNSTRUCTURED_CELL,ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

       ! &      diag%z0m_tile(nproma,nblks_c,nsfc_type)
       cf_desc    = t_cf_var('z0m_tile', '',&
         &'geopotential of the top of the atmospheric boundary layer')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'z0m_tile', diag%z0m_tile,                    &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shapesfc ,&
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.           )
       diag%z0m_tile(:,:,:)=1.e-3_wp

       ALLOCATE(diag%z0m_ptr(nsfc_type))
       DO jsfc = 1,nsfc_type
         WRITE(csfc,'(i1)') jsfc 
         CALL add_ref( diag_list, 'z0m_tile',                                       &
                    & 'z0m_tile_'//csfc, diag%z0m_ptr(jsfc)%p_2d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_heat_'//csfc, '', ''),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
       END DO

       ! &      diag%ustar(nproma,nblks_c)
       cf_desc    = t_cf_var('ustar', 'm s-1','friction velocity')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'ustar', diag%ustar,                         &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

       ! &      diag%kedisp(nproma,nblks_c)
       cf_desc    = t_cf_var('kedisp','','KE dissipation rate')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'kedisp', diag%kedisp,                       &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

       ! &      diag%ocu(nproma,nblks_c)
       cf_desc    = t_cf_var('ocu', 'm s-1','eastward  velocity of ocean surface current')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'ocu', diag%ocu,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

       ! &      diag%ocv(nproma,nblks_c)
       cf_desc    = t_cf_var('ocv','','northward velocity of ocean surface current')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'ocv', diag%ocv,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)



       !------------------
       !Turbulence 3D variables

       ! &      diag%ri(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('ri', '', '  moist Richardson number at layer interfaces')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'ri', diag%ri,                               &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1, lrestart=.FALSE. ) 

       ! &      diag%mixlen(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('mixlen', 'm s-2', 'mixing length at layer interfaces')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'mixlen', diag%mixlen,                       &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
          & ldims=shape3dkp1, lrestart=.FALSE. )

       ! &      diag%thvvar(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('thvvar', 'm s-2', &
         &'variance of virtual potential temperature at layer interfaces')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'thvvar', diag%thvvar,                       &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1 ) 

       ! &      diag%cfm(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('cfm', '', 'turbulent exchange coefficient')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'cfm', diag%cfm,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1, lrestart=.FALSE. ) 

       ! &      diag%cfh(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('cfh', '', 'turbulent exchange coefficient')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'cfh', diag%cfh,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1, lrestart=.FALSE. ) 

       ! &      diag%cfv(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('cfv', '','turbulent exchange coefficient')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'cfv', diag%cfv,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1, lrestart=.FALSE. )

       ! &      diag%cftke(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('cftke', '', 'turbulent exchange coefficient')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'cftke', diag%cftke,                         &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1, lrestart=.FALSE. )

       ! &      diag%cfthv(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('cfthv', '','turbulent exchange coefficient')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'cfthv', diag%cfthv,                         &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,        &
         & ldims=shape3dkp1, lrestart=.FALSE. ) 

    ENDIF  !inwp_turb == 2 (vdiff)


  !
  ! EDMF
  !

    IF( atm_phy_nwp_config(k_jg)%inwp_turb == 3) THEN

       ! &      diag%z0m(nproma,nblks_c)
       cf_desc    = t_cf_var('z0m', '','geopotential of the top of the atmospheric boundary layer')
       grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'z0m', diag%z0m,                             &
         & GRID_UNSTRUCTURED_CELL,ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ENDIF  !inwp_turb == 3 (EDMF)


    !------------------
    !Turbulence 3D variables

   ! &      diag%tkvm(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('tkvm', 'm**2/s', ' turbulent diffusion coefficients for momentum')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
      & GRID_UNSTRUCTURED_CELL,  ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
      & ldims=shape3d, in_group=groups("pbl_vars") )

   ! &      diag%tkvh(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('tkvh', 'm**2/s', ' turbulent diffusion coefficients for heat')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,          &
      & ldims=shape3d, in_group=groups("pbl_vars") ) 

   ! &      diag%rcld(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('rcld', '', 'standard deviation of the saturation deficit')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rcld', diag%rcld,                             &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,          &
      & ldims=shape3dkp1 )

   ! &      diag%edr(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('edr', '', 'eddy dissipation rate')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'edr', diag%edr,                               &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,          &
      & ldims=shape3dkp1, lrestart=.FALSE. ) 


    CALL message('mo_nwp_phy_state:construct_nwp_phy_diag', &
                 'construction of NWP physical fields finished')  

END SUBROUTINE new_nwp_phy_diag_list



SUBROUTINE new_nwp_phy_tend_list( k_jg, klev,  kblks,   &
                     & listname, phy_tend_list, phy_tend)

    INTEGER,INTENT(IN) :: k_jg, klev, kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list)    ,INTENT(INOUT) :: phy_tend_list
    TYPE(t_nwp_phy_tend),INTENT(INOUT) :: phy_tend

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d(3), shape3dkp1(3), shape4d(4)
    INTEGER :: ientr, ktracer

    ientr = 16 ! "entropy" of horizontal slice

    shape3d    = (/nproma, klev  , kblks            /)
    shape3dkp1 = (/nproma, klev+1, kblks            /)
    shape4d    = (/nproma, klev  , kblks, nqtendphy /)


    CALL new_var_list( phy_tend_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( phy_tend_list,             &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    
    !------------------------------
    ! Temperature tendencies
    !------------------------------

   ! &      phy_tend%ddt_temp_radsw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radsw', 'K s-1', &
         &                            'short wave radiative temperature tendency')
    grib2_desc = t_grib2_var(192, 162, 100, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radsw', phy_tend%ddt_temp_radsw,      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_temp_radlw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('temp_tend_radlw', 'K s-1', &
         &                            'long wave radiative temperature tendency')
    grib2_desc = t_grib2_var(192, 162, 101, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radlw', phy_tend%ddt_temp_radlw,      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_temp_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('temp_tend_turb', 'K s-1', &
         &                            'turbulence temperature tendency')
    grib2_desc = t_grib2_var(192, 162, 121, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_turb', phy_tend%ddt_temp_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_temp_drag(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temp_tend_drag', 'K s-1', &
         &                'sso + gwdrag temperature tendency')
    grib2_desc = t_grib2_var(192, 162, 125, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_drag', phy_tend%ddt_temp_drag,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_temp_pconv(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_temp_pconv', 'K s-1', &
         &                            'convection temperature tendency')
    grib2_desc = t_grib2_var(192, 162, 128, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_pconv', phy_tend%ddt_temp_pconv,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )


    !------------------------------
    ! Zonal Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_u_turb(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_turb', 'm s-2', &
         &                            'turbulence tendency of zonal wind')
    grib2_desc = t_grib2_var(192, 162, 119, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_turb', phy_tend%ddt_u_turb,          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_u_sso(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_sso', 'm s-2', &
         &                            'sso tendency of zonal wind')
    grib2_desc = t_grib2_var(192, 162, 123, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_sso', phy_tend%ddt_u_sso,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d, &
                & in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_u_gwd(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_gwd', 'm s-2', &
         &                            'GWD tendency of zonal wind')
    grib2_desc = t_grib2_var(192, 128, 220, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_gwd', phy_tend%ddt_u_gwd,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d, &
                & in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_u_raylfric(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_raylfric', 'm s-2', &
         &                'Rayleigh friction tendency of zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_raylfric', phy_tend%ddt_u_raylfric,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                & ldims=shape3d,  lrestart=.FALSE.)

   ! &      phy_tend%ddt_u_pconv(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_pconv', 'm s-2', &
         &                            'convection tendency of zonal wind')
    grib2_desc = t_grib2_var(192, 162, 126, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_pconv', phy_tend%ddt_u_pconv,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d)


    !------------------------------
    ! Meridional Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_v_turb(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_turb', 'm s-2', &
         &                            'turbulence tendency of meridional wind')
    grib2_desc = t_grib2_var(192, 162, 120, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_turb', phy_tend%ddt_v_turb,          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_v_sso(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_sso', 'm s-2', &
         &                            'sso tendency of meridional wind')
    grib2_desc = t_grib2_var(192, 162, 124, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_sso', phy_tend%ddt_v_sso,                    &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d, &
                & in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_v_gwd(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_gwd', 'm s-2', &
         &                            'GWD tendency of meridional wind')
    grib2_desc = t_grib2_var(192, 128, 221, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_gwd', phy_tend%ddt_v_gwd,                    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d, &
                & in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_v_raylfric(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_raylfric', 'm s-2', &
         &                'Rayleigh friction tendency of meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_raylfric', phy_tend%ddt_v_raylfric,          &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d, lrestart=.FALSE.)

   ! &      phy_tend%ddt_v_pconv(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_pconv', 'm s-2', &
         &                            'convection tendency of meridional wind')
    grib2_desc = t_grib2_var(192, 162, 127, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_pconv', phy_tend%ddt_v_pconv,                &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! Moist tracer tendencies
    !------------------------------

   ! &      phy_tend%ddt_tracer_turb(nproma,nlev,nblks,nqtendphy),          &
    cf_desc    = t_cf_var('ddt_tracer_turb', 's-1', &
         &                            'turbulence tendency of tracers')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_turb', phy_tend%ddt_tracer_turb,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape4d,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ktracer=nqtendphy
    ALLOCATE( phy_tend%tracer_turb_ptr(nqtendphy) )

         !qv
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qv_turb', phy_tend%tracer_turb_ptr(1)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qv_turb', 'kg kg**-1 s**-1',                     &
                    & 'tendency_of_specific_humidity_due_to_turbulence'),            &
                    & t_grib2_var(192, 162, 122, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d, lrestart=.FALSE.)

         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qc_turb', phy_tend%tracer_turb_ptr(2)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qc_turb', 'kg kg**-1 s**-1',                     &
                    & 'tendency_of_specific_cloud_water_due_to_turbulence'),         &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d, lrestart=.FALSE.)

         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qi_turb', phy_tend%tracer_turb_ptr(3)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qi_turb', 'kg kg**-1 s**-1',                     &
                    & 'tendency_of_specific_cloud_ice_due_to_turbulence'),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d, lrestart=.FALSE.)


   ! &      phy_tend%ddt_tracer_pconv(nproma,nlev,nblks,nqtendphy),          &
    cf_desc    = t_cf_var('ddt_tracer_pconv', 's-1', &
         &                            'convective tendency y of tracers')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_pconv', phy_tend%ddt_tracer_pconv,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape4d,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ktracer=nqtendphy
    ALLOCATE( phy_tend%tracer_conv_ptr(nqtendphy) )

         !qv
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qv_conv', phy_tend%tracer_conv_ptr(1)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qv_conv', 'kg kg**-1 s**-1',                     &
                    & 'tendency_of_specific_humidity_due_to_convection'),            &
                    & t_grib2_var(192, 162, 129, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)
         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qc_conv', phy_tend%tracer_conv_ptr(2)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qc_conv', 'kg kg**-1 s**-1',                     &
                    & 'tendency_of_specific_cloud_water_due_to_convection'),         &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)
         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qi_conv', phy_tend%tracer_conv_ptr(3)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qi_conv', 'kg kg**-1 s**-1',                     &
                    & 'tendency_of_specific_cloud_ice_due_to_convection'),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)


    !------------------------------
    ! TKE tendency
    !------------------------------

   ! &      phy_tend%ddt_tke(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_tke', 'm 2 s-3'          , &
         &                            'tendency of turbulent kinetic energy ')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke', phy_tend%ddt_tke,             &
                GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,&
              & ldims=shape3dkp1 )

!
    CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'construction of NWP physical tendency fields finished')

END SUBROUTINE new_nwp_phy_tend_list



!
!-------------------------------------------------------------------------

END MODULE mo_nwp_phy_state
!<
