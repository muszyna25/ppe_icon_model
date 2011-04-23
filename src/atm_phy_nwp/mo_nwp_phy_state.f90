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
USE mo_impl_constants,      ONLY: SUCCESS
USE mo_run_nml,             ONLY: nproma, ntracer, iqcond
USE mo_exception,           ONLY: message, finish,message_text
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: n_dom
USE mo_atm_phy_nwp_nml,     ONLY: inwp_turb
!USE mo_icoham_sfc_indices,  ONLY: nsfc_type, igbm

IMPLICIT NONE
PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!public interface
!
! subroutines
PUBLIC :: construct_nwp_phy_state
PUBLIC :: destruct_nwp_phy_state
!
!variables
PUBLIC :: t_nwp_phy_diag, t_nwp_phy_tend
PUBLIC :: prm_diag, prm_nwp_tend
PUBLIC :: mean_charlen
!
!!data structure defining model states
!
!!diagnostic variables
!

INTEGER :: n_updown = 7 !> number of up/downdrafts variables

REAL(wp), ALLOCATABLE :: mean_charlen(:)

TYPE t_nwp_phy_diag
  REAL(wp), ALLOCATABLE ::  &
       &   tracer_rate(:,:,:) , & !> (nproma,nblks,4) precipitation rate of rain and snow
       &   rain_gsp(:,:),       & !! accumulated grid-scale surface rain
       &   snow_gsp(:,:),       & !! accumulated grid_scale surface snow
       &   rain_con(:,:),       & !! accumulated convective surface rain
       &   snow_con(:,:),       & !! accumulated convective surface snow
       &   tot_prec(:,:),       & !! accumulated grid-scale surface total precipitation
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
       &  lhfl_s(:,:),          & !! latent   heat flux (surface) ( W/m2)
       &  qhfl_s(:,:),          & !!      moisture flux (surface) ( W/m2)
       &  tot_cld(:,:,:,:),     & !! total cloud variables (cc,qv,qc,qi)
       &  tot_cld_vi(:,:,:),    & !! total cloud variables (cc,qv,qc,qi) vertical integrated
       &  cosmu0(:,:),          & !! cosine of solar zenith angle
       &  vio3(:,:),            & !! vertically integrated ozone amount (Pa O3)
       &  hmo3(:,:),            & !! height of O3 maximum (Pa)
       &  flxdwswtoa(:,:),      & !! downward shortwave flux at TOA [W/m2]
       &  tsfctrad(:,:),        & !! surface temperature at trad [K]
       &  lwflxclr(:,:,:),      & !! longwave clear-sky net flux [W/m2]
       &  lwflxall(:,:,:),      & !! longwave net flux           [W/m2]
       &  lwflxsfc(:,:),        & !! longwave net flux at surface [W/m2]
       &  trsolclr(:,:,:),      & !! shortwave clear-sky net tranmissivity []
       &  trsolall(:,:,:),      & !! shortwave net tranmissivity []
       &  swflxsfc(:,:),        & !! shortwave net flux at surface [W/m2]
       &  swflxtoa(:,:),        & !! shortwave net flux at toa [W/m2]
       &  lwflxsfc_avg(:,:),    & !! longwave net flux at surface [W/m2], mean since model start
       &  swflxsfc_avg(:,:),    & !! shortwave net flux at surface [W/m2], mean since model start
       &  lwflxtoa_avg(:,:),    & !! longwave net flux at toa [W/m2], mean since model start
       &  swflxtoa_avg(:,:),    & !! shortwave net flux at toa [W/m2], mean since model start
       &  acdnc(:,:,:)            !! cloud droplet number concentration [1/m**3]
  

 !> Parameter fields for turbulence
  REAL(wp), ALLOCATABLE ::  &
      rcld(:,:,:)      ,    & !> standard deviation of the saturation deficit    --
       tcm(:,:)        ,    & !! turbulent transfer coefficients for momentum    --
       tch(:,:)        ,    & !! turbulent transfer coefficients for heat        --
       tfm(:,:)        ,    & !! factor of laminar transfer of momentum          --
       tfh(:,:)        ,    & !! factor of laminar transfer of scalars           --
       tfv(:,:)         ,   & !! laminar reduction factor for evaporation        --
       gz0(:,:),            & !! roughness length * g of the vertically not
                              !! resolved canopy   !! surface area index       ( 1 )
       sai(:,:),            & !! surface area index                            ( 1 )
       tai(:,:),            & !! transpiration area index                      ( 1 )
       eai(:,:),            & !! (evaporative) earth area index                ( 1 )
       tkvm(:,:,:),         & !! turbulent diffusion coefficients for momentum (m/s2 )
       tkvh(:,:,:),         & !! turbulent diffusion coefficients for heat     (m/s2 )
       h_ice(:,:),          & !! ice thickness                                 (  m  )
       t_2m(:,:)       ,    & !! temperature in 2m                             (  K  )
       qv_2m (:,:)     ,    & !! specific water vapor content in 2m            (kg/kg)
       td_2m (:,:)     ,    & !! dew-point in 2m                               (  K  )
       rh_2m (:,:)     ,    & !! relative humidity in 2m                       (  %  )
       u_10m (:,:)     ,    & !! zonal wind in 10m                             ( m/s )
       v_10m (:,:)     ,    & !! meridional wind in 10m                        ( m/s )
       edr   (:,:,:)          !! eddy dissipation rate
!

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
      & cfthv  (:,:,:),     &!< turbulent exchange coefficient
      & ghpbl (:,:),        &!< geopotential of the top of the atmospheric boundary layer
      & z0m   (:,:,:),      &!< aerodynamic roughness length
                             !< (grid-box mean and over each surface type)
      & ustar (:,:),        &!<
      & kedisp(:,:),        &!< time-mean (or integrated?)
                             !< vertically integrated dissipation of kinetic energy
      & ocu   (:,:),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)          !< northward velocity of ocean surface current


    ! for old aerosol climatology from COSMO (to be used with inwp_radiation==2)
    REAL(wp),ALLOCATABLE :: &
      & aersea  (:,:),      &
      & aerlan  (:,:),      &
      & aerurb  (:,:),      &
      & aerdes  (:,:)
    
!    IF( nextra > 0) THEN
!   REAL(wp),ALLOCATABLE :: &
!      & extra   (:,:,:,nextra) !< experimental diagnostic fields (e.g. cc_conv)
!    ENDIF


INTEGER, ALLOCATABLE ::         &
       &  mbas_con(:,:),        & !!cloud base level index
       &  mtop_con(:,:),        & !! cloud top  level index
       &  ktype   (:,:)           !!  Type of convection

LOGICAL, ALLOCATABLE ::         & !!
       & locum     (:,:)          !! convective  activity indicator
END TYPE t_nwp_phy_diag
!
! !---tendencies of type global!
!
TYPE t_nwp_phy_tend
    REAL(wp),ALLOCATABLE ::  &
             ddt_temp_radsw  (:,:,:)  ,& !! Temp-tendency from shortwave radiation
             ddt_temp_radlw  (:,:,:)  ,& !! Temp-tendency from longwave radiation
             ddt_temp_turb   (:,:,:)  ,& !! Temp-tendency from turbulence
             ddt_temp_gwd    (:,:,:)  ,& !! Temp-tendency from gravity wave drag
             ddt_temp_sso    (:,:,:)  ,& !! Temp-tendency from sso drag
             ddt_temp_pconv  (:,:,:)  ,& !! Temp-tendency from convective prec
             ddt_temp_pscl   (:,:,:)  ,& !! Temp-tendency from grid scale prec
             ddt_u_turb      (:,:,:)  ,& !! ZonalW-tendency from turbulence
             ddt_u_gwd       (:,:,:)  ,& !! ZonalW-tendency from gravity wave drag
             ddt_u_sso       (:,:,:)  ,& !! ZonalW-tendency from sso drag
             ddt_u_pconv     (:,:,:)  ,& !! ZonalW-tendency from convective prec
             ddt_v_turb      (:,:,:)  ,& !! MeridW-tendency from turbulence
             ddt_v_gwd       (:,:,:)  ,& !! MeridW-tendency from gravity wave drag
             ddt_v_sso       (:,:,:)  ,& !! MeridW-tendency from sso drag
             ddt_v_pconv     (:,:,:)  ,& !! MeridW-tendency from convective prec
             ddt_tracer_turb (:,:,:,:),& !! Hydromet-tendency from turbulence
             ddt_tracer_pconv(:,:,:,:),& !! Hydromet-tendency from convective prec
             ddt_tracer_pscl (:,:,:,:),& !! Hydromet-tendency from grid scale prec
             ddt_tke         (:,:,:)     !! tendency for turbulent knetic energy
END TYPE t_nwp_phy_tend

  TYPE(t_nwp_phy_diag), ALLOCATABLE, TARGET :: prm_diag(:)     !< shape: (n_dom)
  TYPE(t_nwp_phy_tend), ALLOCATABLE, TARGET :: prm_nwp_tend(:) !< shape: (n_dom)
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------

SUBROUTINE construct_nwp_phy_state( p_patch)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)

!INTEGER :: jg
INTEGER :: ist
!-----------------------------------------------------------------------

CALL message('mo_nwp_phy_state:construct_nwp_state', &
  'sta rt to construct 3D state vector')

  !This array is only defined via its patches
  ALLOCATE(mean_charlen(n_dom), STAT=ist)
  IF (ist/=success) THEN
  CALL finish('mo_nwp_phys_state:', 'allocation of mean_charlen failed')
  ENDIF

  !create state arrays
  ALLOCATE(prm_diag(n_dom), STAT=ist)
  IF(ist/=SUCCESS)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of diagnostic physical array failed')
  ENDIF

  ALLOCATE(prm_nwp_tend(n_dom), STAT=ist)
  IF(ist/=SUCCESS)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of diagnostic physical array failed')
  ENDIF

  CALL construct_nwp_phy_diag(p_patch)
!
  CALL construct_nwp_phy_tend(p_patch )

  CALL message('mo_nwp_phy_state:construct_nwp_state', &
    'construction of state vector finished')

END SUBROUTINE construct_nwp_phy_state

!
SUBROUTINE destruct_nwp_phy_state

  INTEGER :: ist !> status

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
  'start to destruct 3D state vector')

  CALL destruct_nwp_phy_diag

  CALL destruct_nwp_phy_tend


  !This array is only defined via its patches
  DEALLOCATE(mean_charlen, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_phy_state', &
         &       'deallocation of mean_charlen failed')
  ENDIF

  DEALLOCATE(prm_diag, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_phy_state', &
       &       'deallocation of NWP physics diagnostic array failed')
  ENDIF
 
  DEALLOCATE(prm_nwp_tend, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_phy_state', &
         &       'deallocation of NWP physics tendencies array failed') 
  ENDIF




  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
    'destruction of 3D state vector finished')

END SUBROUTINE destruct_nwp_phy_state


SUBROUTINE construct_nwp_phy_diag(p_patch)
!
!!DESCRIPTION: Allocation of fields neede by the NWP physics
!               Initialization
!
!!INPUT PARAMETERS:
TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)

INTEGER :: jg

!!LOCAL VARIABLES:
INTEGER :: ist, nblks_c, nblks_e, nblks_v
INTEGER :: nlev, nlevp1           !< number of full and half levels
!-----------------------------------------------------------------------

CALL message('mo_nwp_phy_state:construct_nwp_phy_diag', &
             'start to construct NWP physical fields')

DO jg = 1, n_dom

  !determine size of arrays
  nblks_c = p_patch(jg)%nblks_c
  nblks_e = p_patch(jg)%nblks_e
  nblks_v = p_patch(jg)%nblks_v

  ! number of vertical levels
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1


!allocate fields
!
  ALLOCATE(prm_diag(jg)%tracer_rate(nproma,nblks_c,4), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for precipitation rate of rain failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%rain_gsp(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for accumulated grdscl rain failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%snow_gsp(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for accumulated grdscl snow failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%rain_con(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for accumulated rain failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%snow_con(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for accumulated convective snow failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%tot_prec(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for accumulated total precipitaion failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%cape(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for cape failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%con_gust(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for convective gusts failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%rain_upd(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'allocation for rain in updrafts failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%con_udd(nproma,nlev,nblks_c,n_updown), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for convective up/downdraft fields failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%shfl_s(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for surface sensible heat flux failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%lhfl_s(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for surface latent heat flux failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%qhfl_s(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for surface moisture flux failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%tot_cld(nproma,nlev,nblks_c,4), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for tot cloud components failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%tot_cld_vi(nproma,nblks_c,4), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for tot cloud components vert. integ.  failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%cosmu0(nproma,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for cosine of solar zenith angle failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%vio3(nproma,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for vio3 failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%hmo3(nproma,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for hmo3 failed')    
  ENDIF  

  ALLOCATE(prm_diag(jg)%flxdwswtoa(nproma,nblks_c), STAT=ist ) !radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation of downward shortwave flux at TOA failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%tsfctrad(nproma,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for surface temperature at trad failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%lwflxclr(nproma,nlevp1,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for longwave clear-sky net flux failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%lwflxall(nproma,nlevp1,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for longwave net flux failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%lwflxsfc(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for longwave surface net flux failed')    
  ENDIF
  
  ALLOCATE(prm_diag(jg)%trsolclr(nproma,nlevp1,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for shortwave clear-sky net tranmissivity failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%trsolall(nproma,nlevp1,nblks_c), STAT=ist ) !needed for radiation
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for shortwave net tranmissivity failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%swflxsfc(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for shortwave surface net flux failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%swflxtoa(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for shortwave surface net flux failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%lwflxsfc_avg(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for time-avg longwave surface net flux failed')
  ENDIF
    
  ALLOCATE(prm_diag(jg)%swflxsfc_avg(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for time-avg shortwave surface net flux failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%lwflxtoa_avg(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for time-avg longwave toa net flux failed')
  ENDIF
    
  ALLOCATE(prm_diag(jg)%swflxtoa_avg(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for time-avg shortwave toa net flux failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%acdnc(nproma,nlev,nblks_c), STAT=ist ) !needed for inwp_radiation=1
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for cloud droplet number concentration failed')    
  ENDIF

  ALLOCATE(prm_diag(jg)%aersea(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for aersea failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%aerlan(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for aerlan failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%aerurb(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for aerurb failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%aerdes(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for aerdes failed')
  ENDIF
  
  ALLOCATE(prm_diag(jg)%mbas_con(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for level index of con-cloud base failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%mtop_con(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for level index of con-cloud top failed')
  ENDIF
  ALLOCATE(prm_diag(jg)%ktype(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for convective type failed')
  ENDIF

  ALLOCATE(prm_diag(jg)%locum(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'allocation for convective activity field failed')
  ENDIF

!  IF(inwp_turb == 1) THEN
    ALLOCATE(prm_diag(jg)%gz0(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%tcm(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%rcld(nproma,nlevp1,nblks_c),STAT=ist )
    ALLOCATE(prm_diag(jg)%tch(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%tfm(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%tfh(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%tfv(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%sai(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%tai(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%eai(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%tkvm(nproma,nlevp1,nblks_c), STAT=ist ) !remove startlevel2 KF
    ALLOCATE(prm_diag(jg)%tkvh(nproma,nlevp1,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%h_ice(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%t_2m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%td_2m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%qv_2m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%rh_2m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%u_10m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%v_10m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%v_10m(nproma,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)%edr(nproma,nlevp1,nblks_c), STAT=ist )

    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
           'allocation for turbulencs related fields failed')
    ENDIF

!  IF(inwp_turb == 2) THEN

!    WRITE(message_text,'(a,2I4)') 'init value = ',nsfc_type,igbm
!      CALL message('', TRIM(message_text))

    ALLOCATE( prm_diag(jg)% ri     (nproma,nlevp1,nblks_c), STAT=ist )
    ALLOCATE( prm_diag(jg)% mixlen (nproma,nlevp1,nblks_c), STAT=ist )
    ALLOCATE( prm_diag(jg)% thvvar (nproma,nlevp1,nblks_c), STAT=ist )

      IF (ist/=SUCCESS) CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag',   &
        & 'allocation of vdiff  turbulence fieldss failed')

    ALLOCATE(prm_diag(jg)% cfm    (nproma,nlevp1,     nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)% cfm_sfc(nproma,3,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)% cfh    (nproma,nlevp1,     nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)% cfh_sfc(nproma,3,nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)% cfv    (nproma,nlevp1,     nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)% cftke  (nproma,nlevp1,     nblks_c), STAT=ist )
    ALLOCATE(prm_diag(jg)% cfthv  (nproma,nlevp1,     nblks_c), STAT=ist )

    IF (ist/=SUCCESS) CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag',   &
        & 'allocation of turbulence fields failed')


      ALLOCATE( prm_diag(jg)% ghpbl  (nproma,nblks_c),                &
        &       prm_diag(jg)% z0m    (nproma,0:3,nblks_c),              & !KF set a fix value
        &       prm_diag(jg)% ustar  (nproma,nblks_c),                &
        &       prm_diag(jg)% kedisp (nproma,nblks_c),                &
        &       prm_diag(jg)% ocu    (nproma,nblks_c),                &
        &       prm_diag(jg)% ocv    (nproma,nblks_c),                &
        &       STAT=ist                                     )
  IF (ist/=SUCCESS) CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag',   &
        & 'allocation of turbulence fields failed')

!    ENDIF
!  ENDIF

!
!!---------------------------------------------------------------------------------
!! Initialize variables
!!---------------------------------------------------------------------------------
!
  prm_diag(jg)%tracer_rate = 0._wp !>  precipitation rate of rain and snow
  prm_diag(jg)%rain_gsp    = 0._wp !!  accumulated gridscale rain
  prm_diag(jg)%snow_gsp    = 0._wp !!  accumulated gridscale snow
  prm_diag(jg)%rain_con    = 0._wp !!  accumulated convective rain
  prm_diag(jg)%snow_con    = 0._wp !!  accumulated convective snow
  prm_diag(jg)%tot_prec    = 0._wp !!  accumulated total precipitation
  prm_diag(jg)%cape        = 0._wp !!
  prm_diag(jg)%con_gust    = 0._wp !!
  prm_diag(jg)%con_udd     = 0._wp !!
  prm_diag(jg)%rain_upd    = 0._wp !!
  prm_diag(jg)%shfl_s      = 0._wp !!
  prm_diag(jg)%lhfl_s      = 0._wp !!
  prm_diag(jg)%qhfl_s      = 0._wp !!
  prm_diag(jg)%tot_cld     = 0._wp !!
  prm_diag(jg)%tot_cld_vi  = 0._wp !!
  prm_diag(jg)%flxdwswtoa  = 0._wp
  prm_diag(jg)%cosmu0      = 0._wp
  prm_diag(jg)%vio3        = 0._wp
  prm_diag(jg)%hmo3        = 0._wp
  prm_diag(jg)%lwflxclr    = 0._wp
  prm_diag(jg)%lwflxall    = 0._wp
  prm_diag(jg)%lwflxsfc    = 0._wp
  prm_diag(jg)%trsolclr    = 0._wp
  prm_diag(jg)%trsolall    = 0._wp
  prm_diag(jg)%swflxsfc    = 0._wp
  prm_diag(jg)%swflxtoa    = 0._wp
  prm_diag(jg)%lwflxsfc_avg= 0._wp
  prm_diag(jg)%swflxsfc_avg= 0._wp
  prm_diag(jg)%lwflxtoa_avg= 0._wp
  prm_diag(jg)%swflxtoa_avg= 0._wp  
  !
  prm_diag(jg)%acdnc       = 220._wp*1.e6_wp

  prm_diag(jg)%mbas_con     =  0
  prm_diag(jg)%mtop_con     =  0
  prm_diag(jg)%ktype        =  0

  prm_diag(jg)%locum        = .FALSE.

!  IF(inwp_turb== 1) THEN
    prm_diag(jg)%gz0      = 0.01_wp
    prm_diag(jg)%rcld     = 0._wp
    prm_diag(jg)%tcm      = 0._wp
    prm_diag(jg)%tch      = 0._wp
    prm_diag(jg)%tfm      = 0._wp
    prm_diag(jg)%tfh      = 0._wp
    prm_diag(jg)%tfv      = 0._wp
    prm_diag(jg)%sai      = 0._wp
    prm_diag(jg)%tai      = 0._wp
    prm_diag(jg)%eai      = 0._wp
    prm_diag(jg)%tkvm     = 0._wp
    prm_diag(jg)%tkvh     = 0._wp
    prm_diag(jg)%h_ice    = 0._wp
    prm_diag(jg)%t_2m     = 0._wp
    prm_diag(jg)%qv_2m    = 0._wp
    prm_diag(jg)%td_2m    = 0._wp
    prm_diag(jg)%rh_2m    = 0._wp
    prm_diag(jg)%u_10m    = 0._wp
    prm_diag(jg)%v_10m    = 0._wp
    prm_diag(jg)%edr      = 0._wp
 
    IF(inwp_turb == 2) THEN
      prm_diag(jg)% ri      = 0._wp 
      prm_diag(jg)% mixlen  = 0._wp 
      prm_diag(jg)% thvvar  = 0._wp 
      prm_diag(jg)% cfm     = 0._wp 
      prm_diag(jg)% cfm_sfc = 0._wp 
      prm_diag(jg)% cfh     = 0._wp 
      prm_diag(jg)% cfh_sfc = 0._wp 
      prm_diag(jg)% cfv     = 0._wp 
      prm_diag(jg)% cftke   = 0._wp 
      prm_diag(jg)% cfthv   = 0._wp 
      prm_diag(jg)% ghpbl   = 0._wp 
      prm_diag(jg)% z0m     = 0._wp 
      prm_diag(jg)% ustar   = 0._wp 
      prm_diag(jg)% kedisp  = 0._wp 
      prm_diag(jg)% ocu     = 0._wp 
      prm_diag(jg)% ocv     = 0._wp 
    ENDIF

ENDDO

CALL message('mo_nwp_phy_state:construct_nwp_phy_diag', &
  'construction of NWP physical fields finished')

END SUBROUTINE construct_nwp_phy_diag


SUBROUTINE construct_nwp_phy_tend(p_patch)
!
!!DESCRIPTION: Allocation of fields neede by the NWP physics
!               Initialization
!
!!INPUT PARAMETERS:
TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)
!
INTEGER :: jg
!
!!LOCAL VARIABLES:
INTEGER :: ist, nblks_c, nblks_e, nblks_v
INTEGER :: nlev, nlevp1           !< number of full and half levels
!-----------------------------------------------------------------------

CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
             'start to construct NWP physical tendency fields')

DO jg = 1, n_dom

  !determine size of arrays
  nblks_c = p_patch(jg)%nblks_c
  nblks_e = p_patch(jg)%nblks_e
  nblks_v = p_patch(jg)%nblks_v

  ! number of vertical levels
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

!allocate fields
  !

  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_radsw(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for shortwave radiative tendency for temperature failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_radlw(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for longwave radiative tendency for temperature failed')
  ENDIF
  !  
  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_turb(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for turbulence tendency for temperature failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_gwd(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for gwd tendency for temperature failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_sso(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for sso tendency for temperature failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_pconv(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for convective tendency for temperature failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_temp_pscl(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for microphysical tendency for temperature failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_u_turb(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for turbulence tendency for zonal wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_u_gwd(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for gwd tendency for zonal wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_u_sso(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for sso tendency for zonal wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_u_pconv(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for convective tendency for zonal wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_v_turb(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for turbulence tendency for meridional wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_v_gwd(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for gwd tendency for meridional wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_v_sso(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for sso tendency for meridional wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_v_pconv(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for convective tendency for meridional wind failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_tracer_turb(nproma,nlev,nblks_c,iqcond), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for turbulent tendency for tracer failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_tracer_pconv(nproma,nlev,nblks_c,iqcond), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for convective tendency for tracer failed')
  ENDIF
  !
  ALLOCATE(prm_nwp_tend(jg)%ddt_tracer_pscl(nproma,nlev,nblks_c,iqcond), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for microphysical tendency for tracer failed')
  ENDIF

  ALLOCATE(prm_nwp_tend(jg)%ddt_tke(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'allocation for tke tendency for tracer failed')
  ENDIF
!
!! Initialize fields
!
  prm_nwp_tend(jg)%ddt_temp_radsw  (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_temp_radlw  (:,:,:)   = 0._wp  
  prm_nwp_tend(jg)%ddt_temp_turb   (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_temp_gwd    (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_temp_sso    (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_temp_pconv  (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_temp_pscl   (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_u_turb      (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_u_gwd       (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_u_sso       (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_u_pconv     (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_v_turb      (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_v_gwd       (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_v_sso       (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_v_pconv     (:,:,:)   = 0._wp
  prm_nwp_tend(jg)%ddt_tracer_turb (:,:,:,:) = 0._wp
  prm_nwp_tend(jg)%ddt_tracer_pconv(:,:,:,:) = 0._wp
  prm_nwp_tend(jg)%ddt_tracer_pscl (:,:,:,:) = 0._wp
  prm_nwp_tend(jg)%ddt_tke         (:,:,:)   = 0._wp

ENDDO

CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
  'construction of NWP physical tendency fields finished')

END SUBROUTINE construct_nwp_phy_tend


SUBROUTINE destruct_nwp_phy_diag
!
! !DESCRIPTION:  Deallocation
!
! !LOCAL VARIABLES:
INTEGER :: ist, jg

CALL message('mo_nwp_phy_state:destruct_nwp_phy_diag', &
  'start to destruct physical state')
!
DO jg = 1, n_dom
!
  DEALLOCATE(prm_diag(jg)%tracer_rate, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:deconstruct_nwp_phy_diag', &
      'deallocation for precipitation rate of rain failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%rain_gsp, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:deconstruct_nwp_phy_diag', &
      'deallocation for accumulated grdscl rain failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%snow_gsp, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:deconstruct_nwp_phy_diag', &
      'deallocation for accumulated grdscl snow failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%rain_con, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:deconstruct_nwp_phy_diag', &
      'deallocation for accumulated conv. rain failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%snow_con, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:deconstruct_nwp_phy_diag', &
      'deallocation for accumulated conv. snow failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%tot_prec, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:deconstruct_nwp_phy_diag', &
      'deallocation for accumulated total precipitation failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%cape, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'deallocation for cape failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%con_gust, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'deallocation for convective gusts failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%rain_upd, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for rain in updrafts failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%con_udd, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for convective up/downdraft fields failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%shfl_s, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for surface sensible heat flux failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%lhfl_s, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for surface latent heat flux failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%qhfl_s, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for surface moisture flux failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%tot_cld, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for total cloud components failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%tot_cld_vi, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
    'deallocation for total cloud vert. integ. components failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%cosmu0, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for cosine of solar zenith angle failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%vio3, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for vio3 failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%hmo3, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for hmo3 failed')
  ENDIF
  
  DEALLOCATE(prm_diag(jg)%flxdwswtoa, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation of downward shortwave flux at TOA failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%tsfctrad, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for surface temperature at trad failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%lwflxclr, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for longwave clear-sky net flux failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%lwflxall, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for longwave net flux failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%lwflxsfc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for longwave net surface flux failed')
  ENDIF
  
  DEALLOCATE(prm_diag(jg)%trsolclr, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for shortwave clear-sky net tranmissivity failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%trsolall, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for shortwave net tranmissivity failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%swflxsfc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for shortwave net surface flux failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%swflxtoa, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for shortwave net toa flux failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%lwflxsfc_avg, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for averaged longwave net surface flux failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%swflxsfc_avg, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for averaged shortwave net surface flux failed')
  ENDIF
  
  DEALLOCATE(prm_diag(jg)%lwflxtoa_avg, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for averaged longwave net toa flux failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%swflxtoa_avg, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for averaged shortwave net toa flux failed')
  ENDIF
  
  DEALLOCATE(prm_diag(jg)%acdnc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for cloud droplet number concentration failed')
  ENDIF
  
  DEALLOCATE(prm_diag(jg)%mbas_con, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for level index of con-cloud base failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%mtop_con, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
         'deallocation for level index of con-cloud top failed')
  ENDIF
  DEALLOCATE(prm_diag(jg)%ktype, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'deallocation for convective type failed')
  ENDIF

  DEALLOCATE(prm_diag(jg)%locum, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
      'deallocation for convective activity field failed')
  ENDIF

  IF(inwp_turb == 1) THEN
    DEALLOCATE(prm_diag(jg)%gz0,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%rcld,  STAT=ist )
    DEALLOCATE(prm_diag(jg)%tcm,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%tch,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%tfm,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%tfh,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%tfv,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%sai,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%tai,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%eai,   STAT=ist )
    DEALLOCATE(prm_diag(jg)%tkvm,  STAT=ist )
    DEALLOCATE(prm_diag(jg)%tkvh,  STAT=ist )
    DEALLOCATE(prm_diag(jg)%h_ice, STAT=ist )
    DEALLOCATE(prm_diag(jg)%t_2m,  STAT=ist )
    DEALLOCATE(prm_diag(jg)%td_2m, STAT=ist )
    DEALLOCATE(prm_diag(jg)%qv_2m, STAT=ist )
    DEALLOCATE(prm_diag(jg)%rh_2m, STAT=ist )
    DEALLOCATE(prm_diag(jg)%u_10m, STAT=ist )
    DEALLOCATE(prm_diag(jg)%v_10m, STAT=ist )
    DEALLOCATE(prm_diag(jg)%edr ,  STAT=ist )

    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nwp_phy_state:construct_nwp_phy_diag', &
           'deallocation for turbulencs related fields failed')
    ENDIF

    ELSE IF(inwp_turb == 2 )THEN
      DEALLOCATE( prm_diag(jg)% ri,        &
        &         prm_diag(jg)% mixlen,    &
        &         prm_diag(jg)% thvvar,    &
        &         STAT=ist                 )
     DEALLOCATE( prm_diag(jg)% cfm    , &
        &        prm_diag(jg)% cfm_sfc, &
        &        prm_diag(jg)% cfh    , &
        &        prm_diag(jg)% cfh_sfc, &
        &        prm_diag(jg)% cfv    , &
        &        prm_diag(jg)% cftke  , &
        &        prm_diag(jg)% cfthv    )

     DEALLOCATE( prm_diag(jg)% ghpbl,     &
        &        prm_diag(jg)% z0m,       &
        &        prm_diag(jg)% ustar,     &
        &        prm_diag(jg)% kedisp,    &
        &        prm_diag(jg)% ocu,       &
        &        prm_diag(jg)% ocv,       &
        &        STAT=ist          )
   CALL message('mo_nwp_phy_state:destruct_nwp_phy_diag', &
        'destruction of destruct physical state finished')


    ENDIF

!
ENDDO

CALL message('mo_nwp_phy_state:destruct_nwp_phy_diag', &
  'destruction of destruct physical state finished')

END SUBROUTINE destruct_nwp_phy_diag


SUBROUTINE destruct_nwp_phy_tend
!
! !DESCRIPTION:  Deallocation

! !LOCAL VARIABLES:
INTEGER :: ist,jg

CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
  'start to destruct physical state')
!
!
DO jg = 1, n_dom

!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_radsw, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for shortwave radiative tendency for temperature failed')
!  ENDIF
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_radlw, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for longwave radiative tendency for temperature failed')
!  ENDIF
!!
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_turb, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for turbulence tendency for temperature failed')
!  ENDIF
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier2')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_gwd, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for gwd tendency for temperature failed')
!  ENDIF
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier3')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_sso, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for sso tendency for temperature failed')
!  ENDIF
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier4')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_pconv, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for convective tendency for temperature failed')
!  ENDIF
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier5')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_temp_pscl, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for microphysical tendency for temperature failed')
!  ENDIF
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier6')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_u_turb, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for turbulence tendency for zonal wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier7')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_u_gwd, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for gwd tendency for zonal wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier8')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_u_sso, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for sso tendency for zonal wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier9')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_u_pconv, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for convective tendency for zonal wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier10')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_v_turb, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for turbulence tendency for meridional wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier11')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_v_gwd, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for gwd tendency for meridional wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier12')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_v_sso, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for sso tendency for meridional wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier13')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_v_pconv, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for convective tendency for meridional wind failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier14')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_tracer_turb, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for turbulent tendency for tracer failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier15')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_tracer_pconv, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for convective tendency for tracer failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier16')
!!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_tracer_pscl, STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for microphysical tendency for tracer failed')
!  ENDIF
! CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  'hier17')
!
!  CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!  ' destruct physical state tendencies finished')
!
!  DEALLOCATE(prm_nwp_tend(jg)%ddt_tke   ,  STAT=ist )
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for microphysical tendency for tke failes')
!  ENDIF
!
!  IF (ist/=SUCCESS)THEN
!    CALL finish('mo_nwp_phy_state:destruct_nwp_phy_tend', &
!      'deallocation for turbulent tendency for tracer failed')
!  ENDIF
ENDDO

CALL message('mo_nwp_phy_state:destruct_nwp_phy_tend', &
  'destruction of destruct physical state finished')

END SUBROUTINE destruct_nwp_phy_tend

!
!-------------------------------------------------------------------------

END MODULE mo_nwp_phy_state
!<
