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
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
USE mo_run_nml,             ONLY: nproma, ntracer, iqcond
USE mo_exception,           ONLY: message, finish !,message_text
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: n_dom
USE mo_icoham_sfc_indices,  ONLY: nsfc_type
USE mo_linked_list,         ONLY: t_var_list
USE mo_atm_phy_nwp_nml,     ONLY: inwp_turb
USE mo_var_list,            ONLY: default_var_list_settings, &
                                & add_var, add_ref,          &
                                & new_var_list,              &
                                & delete_var_list
USE mo_cf_convention
USE mo_grib2
USE mo_cdi_constants 


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
PUBLIC :: prm_nwp_diag_list, prm_nwp_tend_list  !< variable lists
PUBLIC :: t_ptr_phy
!
!!data structure defining model states
!
!!diagnostic variables
!
INTEGER :: n_updown = 7 !> number of up/downdrafts variables

REAL(wp), POINTER :: mean_charlen(:)


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

  REAL(wp), POINTER ::  &
       &   tracer_rate(:,:,:) , & !> (nproma,nblks,4) precipitation rate of rain and snow
       &   rain_gsp(:,:),       & !! accumulated grid-scale surface rain
       &   snow_gsp(:,:),       & !! accumulated grid_scale surface snow
       &   rain_con(:,:),       & !! accumulated convective surface rain
       &   snow_con(:,:),       & !! accumulated convective surface snow
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
       &  lhfl_s(:,:),          & !! latent   heat flux (surface) ( W/m2)
       &  qhfl_s(:,:),          & !!      moisture flux (surface) ( Kg/m2/s)
                                  !!      = evaporation rate at surface
       &  shfl_s_avg(:,:),      & !! average since model start of shfl_s ( W/m2)
       &  lhfl_s_avg(:,:),      & !! average since model start of lhfl_s ( W/m2)
       &  qhfl_s_avg(:,:),      & !! average since model start of qhfl_s ( Kg/m2/s) 
                                  !! = average of evaporation rate at surface
       &  tot_cld(:,:,:,:),     & !! total cloud variables (cc,qv,qc,qi)
       &  tot_cld_vi(:,:,:),    & !! vertically integrated tot_cld (cc,qv,qc,qi) 
                                  !! for cc, instead of the vertically integrated value, 
                                  !! this is the cloud cover assuming maximum-random overlap
       &  tot_cld_vi_avg(:,:,:),& !! average since model start of the 
                                  !! vertically integrated tot_cld (cc,qv,qc,qi)  
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
       &  lwflxsfc_avg(:,:),    & !! longwave net flux at surface [W/m2], mean since last output
       &  swflxsfc_avg(:,:),    & !! shortwave net flux at surface [W/m2], mean since last output
       &  lwflxtoa_avg(:,:),    & !! longwave net flux at toa [W/m2], mean since last output
       &  swflxtoa_avg(:,:),    & !! shortwave net flux at toa [W/m2], mean since last output
       &  acdnc(:,:,:)            !! cloud droplet number concentration [1/m**3]
  

 !> Parameter fields for turbulence
  REAL(wp), POINTER ::  &
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
       t_2m_s6avg(:,:),     & !! 6 hourly sample 2 m temperature average       (  K  )
       qv_2m (:,:)     ,    & !! specific water vapor content in 2m            (kg/kg)
       qv_2m_s6avg(:,:),     & !! 6 hourly sample 2 m specific water vapor content average   (kg/kg)
       td_2m (:,:)     ,    & !! dew-point in 2m                               (  K  )
       rh_2m (:,:)     ,    & !! relative humidity in 2m                       (  %  )
       u_10m (:,:)     ,    & !! zonal wind in 10m                             ( m/s )
       v_10m (:,:)     ,    & !! meridional wind in 10m                        ( m/s )
       u_10m_s6avg (:,:),    & !! 6 hourly sample 10m zonal wind  average      ( m/s )
       v_10m_s6avg (:,:),    & !! 6 hourly sample 10m  meridional wind average ( m/s )
       edr   (:,:,:)          !! eddy dissipation rate
!

    ! need only for vdiff ++++
    REAL(wp),POINTER :: &
      & ri        (:,:,:),  &!< moist Richardson number at layer interfaces
      & mixlen    (:,:,:),  &!< mixing length at layer interfaces
      & thvvar    (:,:,:)    !< variance of virtual potential temperature at layer interfaces.
                             !< Computed in "vdiff" by solving a prognostic equation of
                             !< the variance. Used for getting "thvsig".

    REAL(wp),POINTER :: &
      & cfm    (:,:,:),     &!< turbulent exchange coefficient
      & cfm_tile(:,:,:),     &!< turbulent exchange coefficient
      & cfh    (:,:,:),     &!< turbulent exchange coefficient
      & cfh_tile(:,:,:),     &!< turbulent exchange coefficient
      & cfv    (:,:,:),     &!< turbulent exchange coefficient
      & cftke  (:,:,:),     &!< turbulent exchange coefficient
      & cfthv  (:,:,:),     &!< turbulent exchange coefficient
      & ghpbl  (:,:),       &!< geopotential of the top of the atmospheric boundary layer
      & z0m_tile(:,:,:),    &!< aerodynamic roughness length
                             !< (grid-box mean and over each surface type)
      & z0m    (:,:),       &!< aerodynamic roughness length
      !< (grid-box mean and over each surface type)
      & ustar (:,:),        &!<
      & kedisp(:,:),        &!< time-mean (or integrated?)
                             !< vertically integrated dissipation of kinetic energy
      & ocu   (:,:),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)          !< northward velocity of ocean surface current


    ! for old aerosol climatology from COSMO (to be used with inwp_radiation==2)
    REAL(wp),POINTER :: &
      & aersea  (:,:),      &
      & aerlan  (:,:),      &
      & aerurb  (:,:),      &
      & aerdes  (:,:)

INTEGER, POINTER ::         &
       &  mbas_con(:,:),        & !!cloud base level index
       &  mtop_con(:,:),        & !! cloud top  level index
       &  ktype   (:,:)           !!  Type of convection

LOGICAL, POINTER ::         & !!
       & locum     (:,:)          !! convective  activity indicator
END TYPE t_nwp_phy_diag
!
! !---tendencies of type global!
!
TYPE t_nwp_phy_tend

    REAL(wp),POINTER ::  &
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
   TYPE(t_ptr_phy),ALLOCATABLE :: tracer_turb_ptr(:)  !< pointer array: one pointer for each component
   TYPE(t_ptr_phy),ALLOCATABLE :: tracer_conv_ptr(:)  !< pointer array: one pointer for each component
   TYPE(t_ptr_phy),ALLOCATABLE :: tracer_pscl_ptr(:)  !< pointer array: one pointer for each component
END TYPE t_nwp_phy_tend

!!--------------------------------------------------------------------------
!!                          STATE VARIABLES 
!!--------------------------------------------------------------------------
  TYPE(t_nwp_phy_diag), ALLOCATABLE :: prm_diag(:)     !< shape: (n_dom)
  TYPE(t_nwp_phy_tend), ALLOCATABLE :: prm_nwp_tend(:) !< shape: (n_dom)
!-------------------------------------------------------------------------
  
!!--------------------------------------------------------------------------
!!                          VARIABLE LISTS
!!--------------------------------------------------------------------------
  TYPE(t_var_list),ALLOCATABLE :: prm_nwp_diag_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),ALLOCATABLE :: prm_nwp_tend_list(:)  !< shape: (n_dom)
 

CONTAINS

!-------------------------------------------------------------------------

SUBROUTINE construct_nwp_phy_state( p_patch)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)

CHARACTER(len=MAX_CHAR_LENGTH) :: listname
INTEGER ::  jg,ist, nblks_c, nlev, nlevp1

!-----------------------------------------------------------------------

CALL message('mo_nwp_phy_state:construct_nwp_state', &
  'start to construct 3D state vector')

  !This array is only defined via its patches
  ALLOCATE(mean_charlen(n_dom), STAT=ist)
  IF (ist/=success) THEN
  CALL finish('mo_nwp_phys_state:', 'allocation of mean_charlen failed')
  ENDIF

  ! Allocate pointer arrays prm_diag_nwp and prm_nwp_tend, 
  ! as well as the corresponding list arrays.

  ALLOCATE(prm_diag(n_dom), prm_nwp_diag_list(n_dom),STAT=ist)
  IF(ist/=SUCCESS)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of diagnostic physical array and list failed')
  ENDIF

  ALLOCATE(prm_nwp_tend(n_dom), prm_nwp_tend_list(n_dom), STAT=ist)
  IF(ist/=SUCCESS)THEN
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
     
     CALL new_nwp_phy_diag_list( nlev, nlevp1, nblks_c,&
                               & TRIM(listname), prm_nwp_diag_list(jg), prm_diag(jg))
     !
     WRITE(listname,'(a,i2.2)') 'prm_tend_of_domain_',jg
     CALL new_nwp_phy_tend_list ( nlev, nblks_c,&
                                & TRIM(listname), prm_nwp_tend_list(jg), prm_nwp_tend(jg))
  ENDDO
  
  CALL message('mo_nwp_phy_state:construct_nwp_state', &
    'construction of state vector finished')

END SUBROUTINE construct_nwp_phy_state

!
SUBROUTINE destruct_nwp_phy_state

  INTEGER :: jg,ist  !< grid level/domain index

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
  'start to destruct 3D state vector')

  DO jg = 1,n_dom
    CALL delete_var_list( prm_nwp_diag_list(jg) )
    CALL delete_var_list( prm_nwp_tend_list (jg) )
  ENDDO

  !This array is only defined via its patches
  DEALLOCATE(mean_charlen, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_phy_state', &
         &       'deallocation of mean_charlen failed')
  ENDIF

  DEALLOCATE(prm_diag, prm_nwp_diag_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_phy_state', &
       &  'deallocation of NWP physics diagnostic array and list failed')
  ENDIF
 
  DEALLOCATE(prm_nwp_tend, prm_nwp_tend_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_phy_state', &
         &' deallocation of NWP physics tendencies array and list failed') 
  ENDIF

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
    'destruction of 3D state vector finished')

END SUBROUTINE destruct_nwp_phy_state

     !
SUBROUTINE new_nwp_phy_diag_list( klev, klevp1, kblks,   &
                     & listname, diag_list, diag)

    INTEGER,INTENT(IN) :: klev, klevp1, kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN)     :: listname
    CHARACTER(len= MAX_CHAR_LENGTH) :: vname_prefix
    CHARACTER(LEN=1)                :: csfc

    TYPE(t_var_list)    ,INTENT(INOUT) :: diag_list
    TYPE(t_nwp_phy_diag),INTENT(INOUT) :: diag

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape4d(4), shapesfc(3)
    INTEGER :: shape3dkp1(3)
    INTEGER :: ientr,  kcloud
    INTEGER :: jsfc 


    ientr = 16 ! bits "entropy" of horizontal slice

    shape2d    = (/nproma,            kblks            /)
    shape3d    = (/nproma, klev,      kblks            /)
    shape3dkp1 = (/nproma, klevp1,    kblks            /)
    shape4d    = (/nproma, klev,      kblks, iqcond    /)

    IF(inwp_turb == 2) THEN
      shapesfc   = (/nproma,          kblks, nsfc_type /)
    ENDIF

    ! Register a field list and apply default settings

    CALL new_var_list( diag_list, TRIM(listname) )
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



    ! &      diag%tracer_rate(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tracer_rate','kg m-2 ','precipitation rate of rain and snow')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tracer_rate', diag%tracer_rate,          &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,cf_desc, grib2_desc, &
          &                                   ldims=(/nproma,kblks,4/),&
          &                                           lrestart=.FALSE. )

    ! &      diag%rain_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp ', 'kg m-2 ', 'gridscale rain ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp', diag%rain_gsp,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
                &                                              ldims=shape2d,&
          &                                           lrestart=.FALSE. )

    ! &      diag%snow_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp', 'kg m-2 ', 'gridscale snow')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp', diag%snow_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
                & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%rain_con(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con', 'kg m-2 ', 'convective rain')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_con', diag%rain_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                &                              ldims=shape2d,lrestart=.FALSE. )

    ! &      diag%snow_con(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con', 'kg m-2', 'convective snow')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_con', diag%snow_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                &   ldims=shape2d ,lrestart=.FALSE. )

    ! &      diag%tot_prec(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec', '', 'total precip')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec', diag%tot_prec,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d,  lrestart=.FALSE. )

    ! &      diag%tot_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_rate_avg', '', 'total precip, time average')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec_rate_avg', diag%tot_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%con_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('con_prec_rate_avg', '', 'convective precip, time average')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_prec_rate_avg', diag%con_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d,  lrestart=.FALSE. )

    ! &      diag%gsp_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('gsp_prec_rate_avg', '', 'gridscale precip, time average')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'gsp_prec_rate_avg', diag%gsp_prec_rate_avg,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
                &  lrestart=.FALSE. )

 !               &  lmiss=.true.,     missval=0._wp                      )
    ! &      diag%cape(nproma,nblks_c)
    cf_desc    = t_cf_var('cape', 'J kg-1 ', 'conv avail pot energy')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cape', diag%cape,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)
 !               &  lmiss=.true.,     missval=0._wp                      )

    ! &      diag%con_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('con_gust', 'm s-1 ', 'convective gusts')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_gust', diag%con_gust,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)
 !               &  lmiss=.true.,     missval=0._wp                      )
   
    ! &      diag%rain_upd(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_upd', 'unit ', 'rain in updroughts')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_upd', diag%rain_upd,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)
 !               &  lmiss=.true.,     missval=0._wp                      )

    ! &      diag%con_udd(nproma,nblks_c)
    cf_desc    = t_cf_var('con_udd', 'unit ', 'convective up/downdraft fields')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_udd', diag%con_udd,                             &
       &                     GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,&
       &                     ldims=(/nproma,klev,kblks,n_updown/),  lrestart=.FALSE. )

    ! &      diag%mbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mbas_con', '', 'cloud base level index')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mbas_con', diag%mbas_con,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,ldims=shape2d ,&
                &  lrestart=.FALSE. )

    ! &      diag%mtop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mtop_con', '', 'cloud top level index')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mtop_con', diag%mtop_con,                  &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,ldims=shape2d ,&
                & lrestart=.FALSE.)

    ! &      diag%locum(nproma,nblks_c)
    cf_desc    = t_cf_var('locum', '', 'convective  activity indicator')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'locum', diag%locum,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,ldims=shape2d,&
                &   lrestart=.FALSE.  )

    ! &      diag%ktype(nproma,nblks_c)
    cf_desc    = t_cf_var('ktype', '', 'convective up/downdraft fields')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ktype', diag%ktype,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc,       &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE. )

   ! &      diag%tot_cld_vi(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tot_cld_vi', 'unit ','vertical integr total cloud variables')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi', diag%tot_cld_vi,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                &                                 ldims=(/nproma,kblks,4/),   &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

    ! fill the seperate variables belonging to the container tot_cld_vi
    kcloud= 4
    ALLOCATE( diag%tci_ptr(kcloud))
    vname_prefix='T'
       
    !CC
    CALL add_ref( diag_list, 'tot_cld_vi',            &
      & TRIM(vname_prefix)//'CC', diag%tci_ptr(1)%p_2d,                  &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                           &
      & t_cf_var(TRIM(vname_prefix)//'CC', '%','total_column_integrated_cloud_cover'), &
      & t_grib2_var(6, 1, 1, ientr, GRID_REFERENCE, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE. )

    !QV
    CALL add_ref( diag_list, 'tot_cld_vi',                        &
      & TRIM(vname_prefix)//'QV', diag%tci_ptr(2)%p_2d,           &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                    &
      & t_cf_var(TRIM(vname_prefix)//'QV', 'kg m**-2','column_integrated_water_vapour'),&
      & t_grib2_var( 1, 64, 1, ientr, GRID_REFERENCE, GRID_CELL), &
      & ldims=shape2d)

    !qc
    CALL add_ref( diag_list, 'tot_cld_vi',              &
      & TRIM(vname_prefix)//'QC', diag%tci_ptr(3)%p_2d,            &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
      & t_cf_var(TRIM(vname_prefix)//'QC', 'kg m**-2',             &
      & 'total_column-integrated_cloud_water'),                    &
      & t_grib2_var( 1, 69, 1, ientr, GRID_REFERENCE, GRID_CELL),  &
      & ldims=shape2d, lrestart=.FALSE.)

    !qi
    CALL add_ref( diag_list, 'tot_cld_vi',&
      & TRIM(vname_prefix)//'QI', diag%tci_ptr(4)%p_2d,            &
      & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
      & t_cf_var(TRIM(vname_prefix)//'qi', 'kg m**-2',             &
      & 'total_column-integrated_cloud_water'),                    &
      & t_grib2_var(1, 70, 1, ientr, GRID_REFERENCE, GRID_CELL),  &
      & ldims=shape2d, lrestart=.FALSE.)
    

   ! &      diag%tot_cld_vi_avg(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tot_cld_vi_avg', 'unit ','vertical integr total cloud variables')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi_avg', diag%tot_cld_vi_avg,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                &                                 ldims=(/nproma,kblks,4/)  , &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

  ! fill the seperate variables belonging to the container tot_cld_vi_avg
    ALLOCATE( diag%tav_ptr(kcloud))
    vname_prefix='avg_'

           !CC
        CALL add_ref( diag_list, 'tot_cld_vi_avg',                        &
                    & TRIM(vname_prefix)//'cc', diag%tav_ptr(1)%p_2d,     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,              &
                    & t_cf_var(TRIM(vname_prefix)//'cc', '','tci_cloud_cover_avg'),   &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d, lrestart=.FALSE. )

           !QV
        CALL add_ref( diag_list, 'tot_cld_vi_avg',               &
                    & TRIM(vname_prefix)//'qv', diag%tav_ptr(2)%p_2d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                              &
                    & t_cf_var(TRIM(vname_prefix)//'qv', '','tci_specific_humidity_avg'), &
                    & t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d, lrestart=.FALSE. )

           !qc
        CALL add_ref( diag_list, 'tot_cld_vi_avg',      &
                    & TRIM(vname_prefix)//'qc', diag%tav_ptr(3)%p_2d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qc', '',                        &
                    & 'tci_specific_cloud_water_content_avg'),                      &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d, lrestart=.FALSE. )

           !qi
        CALL add_ref( diag_list, 'tot_cld_vi_avg',          &
                    & TRIM(vname_prefix)//'qi', diag%tav_ptr(4)%p_2d,              & 
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'qi', '',                       &
                    & 'tci_specific_cloud_ice_content_avg'),                       &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d, lrestart=.FALSE. )

   ! &      diag%tot_cld(nproma,nblks_c,4)
    cf_desc    = t_cf_var('tot_cld', ' ','total cloud variables (cc,qv,qc,qi)')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld', diag%tot_cld,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,&
                &                         ldims=(/nproma,klev,kblks,4/)    ,&
                 & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)


    ALLOCATE( diag%tot_ptr(kcloud))
    vname_prefix='tot_'

           !CC
        CALL add_ref( diag_list, 'tot_cld',                                         &
                    & TRIM(vname_prefix)//'cc', diag%tot_ptr(1)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var(TRIM(vname_prefix)//'cc', '','total_cloud_cover'),   &
                    & t_grib2_var(192, 128, 164, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d)

           !QV
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'qv', diag%tot_ptr(2)%p_3d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var(TRIM(vname_prefix)//'qv', '','total_specific_humidity'),&
                    & t_grib2_var( 0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d)

           !qc
        CALL add_ref( diag_list, 'tot_cld',                                         &
                    & TRIM(vname_prefix)//'qc', diag%tot_ptr(3)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var(TRIM(vname_prefix)//'qc', '',                        &
                    & 'total_specific_cloud_water_content'),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d)

           !qi
        CALL add_ref( diag_list, 'tot_cld',                                         &
                    & TRIM(vname_prefix)//'qi', diag%tot_ptr(4)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var(TRIM(vname_prefix)//'qi', '',                        &
                    & 'total_specific_cloud_ice_content'),                          &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d)


        ! &      diag%acdnc(nproma,nlsacev,nblks_c)
        cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'acdnc', diag%acdnc,                                 &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d, &
          &  initval_r=220.0e6_wp )
        diag%acdnc(:,:,:)=220.0e6_wp

    !------------------
    ! Radiation
    !------------------
    ! 2D variables

        !        diag%cosmu0    (nproma,       nblks),          &
        cf_desc    = t_cf_var('cosmu0', '', '')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'cosmu0', diag%cosmu0,                   &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tsfctrad(nproma,nblks_c)
        cf_desc    = t_cf_var('tsfctrad', 'W m-2 ', 'surface temperature at trad')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tsfctrad', diag%tsfctrad,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
!                &  lmiss=.true.,     missval=0._wp                      )

        ! &       diag% flxdwswtoa(nproma,       nblks),          &
        cf_desc    = t_cf_var('flxdwswtoa', 'W m-2', 'downward shortwave flux at TOA')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'flxdwswtoa', diag%flxdwswtoa,           &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%swflxsfc(nproma,nblks_c)
        cf_desc    = t_cf_var('swflxsfc', 'W m-2', ' shortwave net flux at surface')
        grib2_desc = t_grib2_var(0, 4, 9, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'swflxsfc', diag%swflxsfc,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%swflxtoa(nproma,nblks_c)
        cf_desc    = t_cf_var('swflxtoa', 'W m-2', ' shortwave net flux at TOA')
        grib2_desc = t_grib2_var(0,4, 9, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'swflxtoa', diag%swflxtoa,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%lwflxsfc(nproma,nblks_c)
        cf_desc    = t_cf_var('lwflxsfc', 'W m-2', 'longwave net flux at surface')
        grib2_desc = t_grib2_var(0, 5, 5, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxsfc', diag%lwflxsfc,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%lwflxtoa_avg(nproma,nblks_c)
        cf_desc    = t_cf_var('lwflxtoa', 'W m-2', &
          &                          'longwave net flux at TOA mean since model start')
        grib2_desc = t_grib2_var(0, 5, 5, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxtoa_avg', diag%lwflxtoa_avg,                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%swflxtoa_avg(nproma,nblks_c)
        cf_desc    = t_cf_var('swflxtoa', 'W m-2', &
          &                         'shortwave net flux at TOA mean since model start')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'swflxtoa_avg', diag%swflxtoa_avg,                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%lwflxsfc_avg(nproma,nblks_c)
        cf_desc    = t_cf_var('lwflxsfc', 'W m-2', &
          &                      'longwave net flux at surface mean since model start')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxsfc_avg', diag%lwflxsfc_avg,                      &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%swflxsfc_avg(nproma,nblks_c)
        cf_desc    = t_cf_var('swflxsfc', 'W m-2', &
          &                     'shortwave net flux at surface mean since model start')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'swflxsfc_avg', diag%swflxsfc_avg,                      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)
 !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%vio3(nproma,nblks_c)
        cf_desc    = t_cf_var('vio3', '', 'vertically integrated ozone amount')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'vio3', diag%vio3,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%hmo3(nproma,nblks_c)
        cf_desc    = t_cf_var('hmo3', 'Pa', 'height of O3 maximum (Pa)')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'hmo3', diag%hmo3,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
!                &  lmiss=.true.,     missval=0._wp                      )


   ! for old aerosol climatology from COSMO (to be used with inwp_radiation==2)

        ! &      diag%aersea(nproma,nblks_c)
        cf_desc    = t_cf_var('aersea', '', '')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aersea', diag%aersea,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) 
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%aerlan(nproma,nblks_c)
        cf_desc    = t_cf_var('aerlan', '', '')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aerlan', diag%aerlan,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) 
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%aerurb(nproma,nblks_c)
        cf_desc    = t_cf_var('aerurb', '', '')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aerurb', diag%aerurb,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) 
        !                &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%aerdes(nproma,nblks_c)
        cf_desc    = t_cf_var('aerdes', '', '')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'aerdes', diag%aerdes,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) 
        !                &  lmiss=.true.,     missval=0._wp                      )

        !------------------
        !Radiation 3D variables

        ! &      diag%lwflxclr(nproma,nlev,nblks_c)
        cf_desc    = t_cf_var('lwflxclr', 'W m-2', 'longwave clear-sky net flux')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxclr', diag%lwflxclr,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )
        !                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag% lwflxall(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('lwflxall', 'W m-2 ', 'longwave net flux')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxall', diag%lwflxall,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 )
        !                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag% trsolclr(nproma,nlev,nblks_c)
        cf_desc    = t_cf_var('trsolclr', '', 'shortwave clear-sky net tranmissivity')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsolclr', diag%trsolclr,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d)
!                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%  trsolall(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('trsolall', '', 'shortwave net tranmissivity')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsolall', diag%trsolall,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1)
        !                &  lmiss=.true.,     missval=0._wp                      )

        !------------------
        !Turbulence 2D variables
        
        ! &      diag%shfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('shfl_s', 'W m-2 ', 'surface sensible heat flux')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'shfl_s', diag%shfl_s,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        cf_desc    = t_cf_var('shfl_s_avg', 'W m-2 ', 'surface sensible heat flux, time avg')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'shfl_s_avg', diag%shfl_s_avg,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%lhfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('lhfl_s', 'W m-2 ', 'surface latent heat flux')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lhfl_s', diag%lhfl_s,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

        cf_desc    = t_cf_var('lhfl_s_avg', 'W m-2 ', 'surface latent heat flux, time avg')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lhfl_s_avg', diag%lhfl_s_avg,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%qhfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('qhfl_s', 'Kg m-2 s-1', 'surface moisture flux')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qhfl_s', diag%qhfl_s,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

        cf_desc    = t_cf_var('qhfl_s_avg', 'Kg m-2 s-1', 'surface moisture flux')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qhfl_s_avg', diag%qhfl_s_avg,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tcm(nproma,nblks_c)
        cf_desc    = t_cf_var('tcm', ' ','turbulent transfer coefficients for momentum')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tcm', diag%tcm,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
!                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tch(nproma,nblks_c)
        cf_desc    = t_cf_var('tch', ' ','turbulent transfer coefficients for heat')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tch', diag%tch,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )
        
        ! &      diag%tfm(nproma,nblks_c)
        cf_desc    = t_cf_var('tfm', ' ','factor of laminar transfer of momentum')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfm', diag%tfm,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tfh(nproma,nblks_c)
        cf_desc    = t_cf_var('tfh', ' ','factor of laminar transfer of scalars')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfh', diag%tfh,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tfv(nproma,nblks_c)
        cf_desc    = t_cf_var('tfv', ' ','laminar reduction factor for evaporation')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfv', diag%tfv,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tgz0(nproma,nblks_c)
        cf_desc    = t_cf_var('gz0', ' ','roughness length time gravity')
        grib2_desc = t_grib2_var(2, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'gz0', diag%gz0,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d ,&
          &  initval_r=0.01_wp                      )
        diag%gz0(:,:)=0.01_wp


        ! &      diag%sai(nproma,nblks_c)
        cf_desc    = t_cf_var('sai', ' ','surface area index')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sai', diag%sai,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%tai(nproma,nblks_c)
        cf_desc    = t_cf_var('tai', ' ','transpiration area index')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tai', diag%tai,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%eai(nproma,nblks_c)
        cf_desc    = t_cf_var('eai', ' ','(evaporative) earth area index')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'eai', diag%eai,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d)! ,  &
 !               &  lmiss=.true.,     missval=0.01_wp                      )

        ! &      diag%t_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('t_2m', 'K ','temperature in 2m')
        grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 't_2m', diag%t_2m,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%t_2m_s6avg(nproma,nblks_c)
        cf_desc    = t_cf_var('t_2m_s6avg', 'K ','temperature in 2m')
        grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 't_2m_s6avg', diag%t_2m_s6avg,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !               &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%qv_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('qv_2m', 'kg kg-1 ','specific water vapor content in 2m')
        grib2_desc = t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qv_2m', diag%qv_2m,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
!                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%qv_2m_s6avg(nproma,nblks_c)
        cf_desc    = t_cf_var('qv_2m_s6avg', 'kg kg-1 ','specific water vapor content in 2m')
        grib2_desc = t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qv_2m_s6avg', diag%qv_2m_s6avg,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%rh_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('rh_2m', '%','relative humidity in 2m')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'rh_2m', diag%rh_2m,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
        !                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%td_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('td_2m', ' ','dew-point in 2m')
        grib2_desc = t_grib2_var(0, 0, 6, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'td_2m', diag%td_2m,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
!                &  lmiss=.true.,     missval=0._wp                      )

        ! &      diag%u_10m(nproma,nblks_c)
        cf_desc    = t_cf_var('u_10m', ' ','zonal wind in 10m')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'u_10m', diag%u_10m,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&

  ! &      diag%u_10m_s6avg(nproma,nblks_c)
    cf_desc    = t_cf_var('u_10m_s6avg', ' ','zonal wind in 10m')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'u_10m_s6avg', diag%u_10m_s6avg,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
!                &  lmiss=.true.,     missval=0._wp                      )

  ! &      diag%v_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('v_10m', ' ','meridional wind in 10m')
    grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'v_10m', diag%v_10m,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

  ! &      diag%v_10m_s6avg(nproma,nblks_c)
    cf_desc    = t_cf_var('v_10m_s6avg', ' ','meridional wind in 10m')
    grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'v_10m_s6avg', diag%v_10m_s6avg,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

  ! &      diag%h_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('h_ice', 'm','ice thickness')
    grib2_desc = t_grib2_var(10, 2, 1, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'h_ice', diag%h_ice,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )


  ! +++vdiff
  !
  !
    IF(inwp_turb == 2) THEN

  ! &      diag%cfm_tile(nproma,nblks_c)
    cf_desc    = t_cf_var('cfm_tile','',&
               & 'turbulent exchange coefficient of momentum at surface')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cfm_tile', diag%cfm_tile,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shapesfc ,&
               & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

    ALLOCATE(diag%cfm_ptr(nsfc_type))
      DO jsfc = 1,nsfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'cfm_tile',                                       &
                    & 'cfm_tile_'//csfc, diag%cfm_ptr(jsfc)%p_2d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_momentum_'//csfc, '', ''),       &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

  ! &      diag%cfh_tile(nproma,nblks_c)
    cf_desc    = t_cf_var('cfh_tile', '','turbulent exchange coefficient of heat at surface')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cfh_tile', diag%cfh_tile,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
               & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )

     ALLOCATE(diag%cfh_ptr(nsfc_type))
      DO jsfc = 1,nsfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'cfh_tile',                                       &
                    & 'cfh_tile_'//csfc, diag%cfh_ptr(jsfc)%p_2d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                       &
                    & t_cf_var('turb_exchng_coeff_heat_'//csfc, '', ''),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

  ! &      diag%ghpbl(nproma,nblks_c)
    cf_desc    = t_cf_var('gh_pbl','','turbulent exchange coefficient of momentum at surface')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ghpbl', diag%ghpbl,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

  ! &      diag%z0m_tile(nproma,nblks_c,nsfc_type)
    cf_desc    = t_cf_var('z0m_tile', '',&
                         &'geopotential of the top of the atmospheric boundary layer')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'z0m_tile', diag%z0m_tile,                    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shapesfc ,&
                &  initval_r=1.e-3_wp ,&
               & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.           )
!         diag%z0m_tile(:,:,:)=1.e-3_wp

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

    ENDIF  !inwp_turb == 2 (vdiff)



  ! &      diag%z0m(nproma,nblks_c)
    cf_desc    = t_cf_var('z0m', '','geopotential of the top of the atmospheric boundary layer')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'z0m', diag%z0m,                             &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d )!,&
!                &  lmiss=.true.,     missval=0._wp                      )

    ! &      diag%ustar(nproma,nblks_c)
    cf_desc    = t_cf_var('ustar', 'm s-1','friction velocity')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ustar', diag%ustar,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

  ! &      diag%kedisp(nproma,nblks_c)
    cf_desc    = t_cf_var('kedisp','','KE dissipation rate')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'kedisp', diag%kedisp,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

    ! &      diag%ocu(nproma,nblks_c)
    cf_desc    = t_cf_var('ocu', 'm s-1','eastward  velocity of ocean surface current')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ocu', diag%ocu,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
 !               &  lmiss=.true.,     missval=0._wp                      )

  ! &      diag%ocv(nproma,nblks_c)
    cf_desc    = t_cf_var('ocv','','northward velocity of ocean surface current')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ocv', diag%ocv,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d) !,&
  !              &  lmiss=.true.,     missval=0._wp                      )

    !------------------
    !Turbulence 3D variables

   ! &      diag%tkvm(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('tkvm', 'm s-2', ' turbulent diffusion coefficients for momentum')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
                & GRID_UNSTRUCTURED_CELL,  ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d) !,&
  !              &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%tkvh(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('tkvh', 'm s-2', ' turbulent diffusion coefficients for heat')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d) !,&
  !              &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%rcld(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('rcld', '', 'standard deviation of the saturation deficit')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rcld', diag%rcld,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d) !,&
   !             &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%edr(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('edr', '', 'eddy dissipation rate')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'edr', diag%edr,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) !,&
    !            &  lmiss=.true.,     missval=0._wp                      )

    ! need only for vdiff ++++
    ! &      diag%ri(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('ri', '', '  moist Richardson number at layer interfaces')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ri', diag%ri,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) !,&
   !             &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%mixlen(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('mixlen', 'm s-2', 'mixing length at layer interfaces')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mixlen', diag%mixlen,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) !,&
   !             &  lmiss=.true.,     missval=0._wp                      )

    ! &      diag%thhvvar(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('thvvar', 'm s-2', &
                     &'variance of virtual potential temperature at layer interfaces')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'thvvar', diag%thvvar,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) !,&
    !            &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%cfm(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('cfm', '', 'turbulent exchange coefficient')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cfm', diag%cfm,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) !,&
    !            &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%cfh(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('cfh', '', 'turbulent exchange coefficient')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cfh', diag%cfh,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) !,&
   !             &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%cfv(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('cfv', '','turbulent exchange coefficient')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cfv', diag%cfv,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) 
    !            &  lmiss=.true.,     missval=0._wp                      )

   ! &      diag%cftke(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('cftke', '', 'turbulent exchange coefficient')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cftke', diag%cftke,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 )

   ! &      diag%cfthv(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('cfthv', '','turbulent exchange coefficient')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cfthv', diag%cfthv,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3dkp1 ) 
 !               &  lmiss=.true.,     missval=0._wp                      )

    CALL message('mo_nwp_phy_state:construct_nwp_phy_diag', &
                 'construction of NWP physical fields finished')  

END SUBROUTINE new_nwp_phy_diag_list


SUBROUTINE new_nwp_phy_tend_list( klev,  kblks,   &
                     & listname, phy_tend_list, phy_tend)

    INTEGER,INTENT(IN) :: klev, kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list)    ,INTENT(INOUT) :: phy_tend_list
    TYPE(t_nwp_phy_tend),INTENT(INOUT) :: phy_tend

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d(3), shape4d(4)
    INTEGER :: ientr, ktracer, jtrc

    ientr = 16 ! "entropy" of horizontal slice

    shape3d    = (/nproma, klev,  kblks        /)
    shape4d    = (/nproma, klev,  kblks, iqcond/)


    CALL new_var_list( phy_tend_list, TRIM(listname) )
    CALL default_var_list_settings( phy_tend_list,             &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    
    !------------------------------
    ! Temperature tendencies
    !------------------------------

   ! &      phy_tend%ddt_temp_radsw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radsw', 'K s-1', &
         &                            'short wave radiative temperature tendency')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radsw', phy_tend%ddt_temp_radsw,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_temp_radlw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('temp_tend_radlw', 'K s-1', &
         &                            'long wave radiative temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radlw', phy_tend%ddt_temp_radlw,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d  )

   ! &      phy_tend%ddt_temp_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('temp_tend_turb', 'K s-1', &
         &                            'turbulence temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_turb', phy_tend%ddt_temp_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_temp_sso(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temp_tend_sso', 'K s-1', &
         &                            'sso temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_sso', phy_tend%ddt_temp_sso,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_temp_gwd(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_temp_gwd', 'K s-1', &
         &                            'GWD temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_gwd', phy_tend%ddt_temp_gwd,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d ) !,&
!                &  lmiss=.true.,     missval=0._wp                      )

   ! &      phy_tend%ddt_temp_pconv(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_temp_pconv', 'K s-1', &
         &                            'convection temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_pconv', phy_tend%ddt_temp_pconv,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_temp_pscl(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_temp_pscl', 'K s-1', &
         &                            'cloud microphysical temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_pscl', phy_tend%ddt_temp_pscl,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )


    !------------------------------
    ! Zonal Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_u_turb(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_turb', 'm s-2', &
         &                            'turbulence tendency of zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_turb', phy_tend%ddt_u_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d)

   ! &      phy_tend%ddt_u_sso(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_sso', 'm s-2', &
         &                            'sso tendency of zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_sso', phy_tend%ddt_u_sso,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_u_gwd(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_gwd', 'm s-2', &
         &                            'GWD tendency of zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_gwd', phy_tend%ddt_u_gwd,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_u_pconv(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_u_pconv', 'm s-2', &
         &                            'convection tendency of zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_pconv', phy_tend%ddt_u_pconv,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d)


    !------------------------------
    ! Meridional Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_v_turb(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_turb', 'm s-2', &
         &                            'turbulence tendency of meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_turb', phy_tend%ddt_v_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d)

   ! &      phy_tend%ddt_v_sso(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_sso', 'm s-2', &
         &                            'sso tendency of meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_sso', phy_tend%ddt_v_sso,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_v_gwd(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_gwd', 'm s-2', &
         &                            'GWD tendency of meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_gwd', phy_tend%ddt_v_gwd,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

   ! &      phy_tend%ddt_v_pconv(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('ddt_v_pconv', 'm s-2', &
         &                            'convection tendency of meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_pconv', phy_tend%ddt_v_pconv,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! Moist tracer tendencies
    !------------------------------

   ! &      phy_tend%ddt_tracer_turb(nproma,nlev,nblks,iqcond),          &
    cf_desc    = t_cf_var('ddt_tracer_turb', 's-1', &
         &                            'turbulence tendency of tracers')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_turb', phy_tend%ddt_tracer_turb,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape4d,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)


    ktracer=iqcond
    ALLOCATE( phy_tend%tracer_turb_ptr(ktracer) )

         !qv
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qv_turb', phy_tend%tracer_turb_ptr(1)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var('ddt_qv_turb', 'kg kg**-1 s**-1',                               &
                    & 'tendency_of_specific_humidity_due_to_turbulence'),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d)

         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qc_turb', phy_tend%tracer_turb_ptr(2)%p_3d,              &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var('ddt_qc_turb', 'kg kg**-1 s**-1',                                  &
                    & 'tendency_of_specific_cloud_water_due_to_turbulence'),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),    &
                    & ldims=shape3d)

         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qi_turb', phy_tend%tracer_turb_ptr(3)%p_3d,            &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qi_turb', 'kg kg**-1 s**-1',                                &
                    & 'tendency_of_specific_cloud_ice_due_to_turbulence'),  &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)


   ! &      phy_tend%ddt_tracer_pconv(nproma,nlev,nblks,iqcond),          &
    cf_desc    = t_cf_var('ddt_tracer_pconv', 's-1', &
         &                            'convective tendency y of tracers')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_pconv', phy_tend%ddt_tracer_pconv,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape4d,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

    ktracer=iqcond
    ALLOCATE( phy_tend%tracer_conv_ptr(ktracer) )

         !qv
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qv_conv', phy_tend%tracer_conv_ptr(1)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var('ddt_qv_conv', 'kg kg**-1 s**-1',                    &
                    & 'tendency_of_specific_humidity_due_to_convection'),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d)
         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qc_conv', phy_tend%tracer_conv_ptr(2)%p_3d,              &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                            &
                    & t_cf_var('ddt_qc_conv', 'kg kg**-1 s**-1',                       &
                    & 'tendency_of_specific_cloud_water_due_to_convection'),           &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),    &
                    & ldims=shape3d)
         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qi_conv', phy_tend%tracer_conv_ptr(3)%p_3d,            &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('ddt_qi_conv', 'kg kg**-1 s**-1',                      &
                    & 'tendency_of_specific_cloud_ice_due_to_convection'),            &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)


   ! &      phy_tend%ddt_tracer_pscl(nproma,nlev,nblks,iqcond)
    cf_desc    = t_cf_var('ddt_tracer_pscl', 's-1', &
         &                            'cloud microphysical tendency y of tracers')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_pscl', phy_tend%ddt_tracer_pscl,        &
                & GRID_UNSTRUCTURED_CELL,ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape4d ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

    ktracer=ntracer
    ALLOCATE( phy_tend%tracer_pscl_ptr(ktracer) )

          !qv
        CALL add_ref( phy_tend_list, 'ddt_tracer_pscl',                                 &
                    & 'ddt_qv_pscl', phy_tend%tracer_pscl_ptr(1)%p_3d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                             &
                    & t_cf_var('ddt_qv_pscl', 'kg kg**-1 s**-1',                        &
                    & 'tendency_of_specific_humidity_due_to_grid_scale_precipitation'), &
                    & t_grib2_var(1, 200, 105, ientr, GRID_REFERENCE, GRID_CELL),       &
                    & ldims=shape3d)
         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_pscl',                                   &
                    & 'ddt_qc_pscl', phy_tend%tracer_pscl_ptr(2)%p_3d,                    &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                               &
                    & t_cf_var('ddt_qc_pscl', 'kg kg**-1 s**-1',                          &
                    & 'tendency_of_specific_cloud_water_due_to_grid_scale_precipitation'),&
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),       &
                    & ldims=shape3d)
         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_pscl',                                 &
                    & 'ddt_qi_pscl', phy_tend%tracer_pscl_ptr(3)%p_3d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                             &
                    & t_cf_var('ddt_qi_pscl', 'kg kg**-1 s**-1',                        &
                    & 'tendency_of_specific_cloud_ice_due_to_grid_scale_precipitation'),&
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),     &
                    & ldims=shape3d)
         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_pscl',                                  &
                    & 'ddt_qr_pscl', phy_tend%tracer_pscl_ptr(4)%p_3d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                              &
                    & t_cf_var('ddt_qr_pscl', 'kg kg**-1 s**-1',                         &
                    & 'tendency_of_rain_due_to_grid_scale_precipitation'),               &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),      &
                    & ldims=shape3d)
         !qs
        CALL add_ref( phy_tend_list, 'ddt_tracer_pscl',                                  &
                    & 'ddt_qs_pscl', phy_tend%tracer_pscl_ptr(5)%p_3d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                              &
                    & t_cf_var('ddt_qs_pscl', 'kg kg**-1 s**-1',                         &
                    & 'tendency_of_snow_due_to_grid_scale_precipitation'),               &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),      &
                    & ldims=shape3d)

    !------------------------------
    ! TKE tendency
    !------------------------------

   ! &      phy_tend%ddt_tke(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_tke', 'm 2 s-3'          , &
         &                            'tendency of turbulent kinetic energy ')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke', phy_tend%ddt_tke,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, ldims=shape3d )

!
    CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'construction of NWP physical tendency fields finished')

END SUBROUTINE new_nwp_phy_tend_list



!
!-------------------------------------------------------------------------

END MODULE mo_nwp_phy_state
!<
