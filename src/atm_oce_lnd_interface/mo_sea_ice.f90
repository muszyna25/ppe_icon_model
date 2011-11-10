!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Dirk Notz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_sea_ice
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
USE mo_kind,                ONLY: wp
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: dtime
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch
USE mo_exception,           ONLY: finish, message
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell 
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_utilities,      ONLY: t_cartesian_coordinates
USE mo_physical_constants,  ONLY: rhoi, rhos, rhow,ki,ks,tf,albi,albim,albsm,albs,&
 &                                mus,ci, Lfreez, I_0, Lsub, Lvap, albedoW, clw,Sice,&
 &                                cpa, emiss,fr_fac,rgas, stefbol,tmelt   
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_ocean_nml,           ONLY: no_tracer, init_oce_prog
USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base
USE mo_oce_index,           ONLY: print_mxmn, ipl_src
IMPLICIT NONE

PRIVATE


CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface


! Definition of forcing types
! public types
PUBLIC  :: t_sea_ice
PUBLIC  :: t_sfc_flx
PUBLIC  :: t_atmos_fluxes
PUBLIC  :: t_atmos_for_ocean

! public subroutines
PUBLIC :: construct_sea_ice 
PUBLIC :: destruct_sea_ice
PUBLIC :: construct_sfcflx
PUBLIC :: construct_atmos_for_ocean
PUBLIC :: construct_atmos_fluxes
PUBLIC :: destruct_sfcflx
PUBLIC :: destruct_atmos_for_ocean
PUBLIC :: destruct_atmos_fluxes

PUBLIC :: ice_init
PUBLIC :: ice_growth
PUBLIC :: set_ice_temp
PUBLIC :: set_ice_albedo
PUBLIC :: sum_fluxes
PUBLIC :: ave_fluxes
PUBLIC :: ice_fast
PUBLIC :: ice_slow
PUBLIC :: upper_ocean_TS
PUBLIC :: new_ice_growth
PUBLIC :: calc_atm_fluxes_from_bulk

!to be put into namelist
INTEGER :: i_no_ice_thick_class = 1

!------  Definition of surface flux type---------------------
TYPE t_sfc_flx

! The forcing is specified as fluxes at the air-sea interface defined on cell-centers
! dimension: (nproma, nblks_c)
  REAL(wp), POINTER ::           &
    &  forc_wind_u(:,:),         & ! forcing of zonal component of velocity equation,
    &  forc_wind_v(:,:),         & ! forcing of meridional component of velocity equation,
    &  forc_hflx(:,:),           & ! forcing of temperature tracer with surface heat flux [W/m2]
    &  forc_fwfx(:,:),           & ! forcing of salinity tracer with surface freshwater flux [m/s]
    &  forc_swflx(:,:),          & ! surface short wave heat flux [W/m2]
    &  forc_lwflx(:,:),          & ! surface long wave heat flux [W/m2]
    &  forc_ssflx(:,:),          & ! surface sensible heat flux [W/m2]
    &  forc_slflx(:,:),          & ! surface latent heat flux [W/m2]
    &  forc_prflx(:,:),          & ! total precipitation flux [m/s]
    &  forc_evflx(:,:),          & ! evaporation flux [m/s]
    &  forc_tracer(:,:,:),       & ! tracer flux. Last index refers to tracer id
    &  forc_tracer_relax(:,:,:)    ! tracer relaxation: contains data to which is relaxated.
                                   ! Last index refers to tracer id (1=temperature, 2=salinity)
  TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
  & ALLOCATABLE :: forc_wind_cc(:,:) 

END TYPE t_sfc_flx

! global type variables
TYPE(t_sfc_flx), PUBLIC, TARGET :: v_sfc_flx


!------  Definition of representation of atm state in ocean model---
!
!representation of atmosphere in ocean model. Data are coming either from
!atmosphere model or from file. These fields are transformed via bulk fomulas
!into atmospheric fluxes, the fluxes are then used to set the oceans surface 
!boundary conditions
TYPE t_atmos_for_ocean

  REAL(wp), ALLOCATABLE :: &
      tafo(:,:),   &  ! 2 m air temperature                              [C]
      ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
      fclou(:,:),  &  ! Fractional cloud cover
      fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
      fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]
      pao(:,:),    &  !Surface atmospheric pressure                      [hPa]
      u(:,:),      &  !wind in reference height                          [m/s]
      v(:,:)       

END TYPE t_atmos_for_ocean



!------  Definition of forcing---------------------
TYPE t_atmos_fluxes

  REAL(wp), ALLOCATABLE :: &
    sens(:,:,:),           & ! Sensible heat flux at ice surface        [W/m2]
    lat(:,:,:),            & ! Latent heat flux at ice surface          [W/m2]
    LWout(:,:,:),          & ! outgoing LW radiation flux at ice surface[W/m2]
    LWnet(:,:,:),          & ! net LW radiation flux at ice surface     [W/m2]
    bot(:,:,:),            & ! Ocean heat flux at ice bottom            [W/m2]
    dsensdT(:,:,:),        & ! d sensible Flux / d T_surf             [W/m2/K]
    dlatdT(:,:,:),         & ! d latent Flux / d T_surf               [W/m2/K]
    dLWdT(:,:,:)             ! d radiation Flux / d T_surf            [W/m2/K]

  REAL(wp), ALLOCATABLE :: &
    rprecw(:,:),           & ! liquid precipitation rate                [m/s]
    rpreci(:,:),           & ! solid  precipitation rate                [m/s]
    sensw(:,:),            & ! Sensible heat flux over water            [W/m2]
    latw(:,:),             & ! Latent heat flux over water              [W/m2]
    LWoutw(:,:),           & ! outgoing LW radiation flux over water    [W/m2]
    LWnetw(:,:),           & ! net LW radiation flux over water         [W/m2]
    SWin(:,:),             & ! incoming SW radiation flux               [W/m2]
    LWin(:,:)                ! incoming LW radiation flux               [W/m2]

  INTEGER ::     counter

  REAL(wp), ALLOCATABLE ::   &
    &  forc_wind_u(:,:),     & ! forcing of zonal component of velocity equation,
    &  forc_wind_v(:,:),     & ! forcing of meridional component of velocity equation,
    &  forc_swflx(:,:),      & ! surface short wave heat flux [W/m2]
    &  forc_lwflx(:,:),      & ! surface long wave heat flux [W/m2]
    &  forc_ssflx(:,:),      & ! surface sensible heat flux [W/m2]
    &  forc_slflx(:,:),      & ! surface latent heat flux [W/m2]
    &  forc_prflx(:,:),      & ! total precipitation flux [m/s]
    &  forc_evflx(:,:),      & ! evaporation flux [m/s]
    &  forc_hflx(:,:),       & ! forcing of temperature tracer with surface heat flux [W/m2]
    &  forc_fwfx(:,:),       & ! forcing of salinity tracer with surface freshwater flux [m/s]
    &  forc_tracer(:,:,:),   & ! tracer flux. Last index refers to tracer id
    &  forc_tracer_relax(:,:,:) ! tracer relaxation: contains data to which is relaxated. 
                                ! Last index refers to tracer id (1=temperature, 2=salinity)
  TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
  & ALLOCATABLE :: forc_wind_cc(:,:) 

END TYPE t_atmos_fluxes

TYPE t_sea_ice

! The description of the sea-ice state, defined on cell-centers
! dimension: (nproma, nblks_c)

  LOGICAL, ALLOCATABLE :: &
    isice(:,:,:)    ! Logical field that marks ice-covered grid cells
  
  REAL(wp), ALLOCATABLE :: &
    alb(:,:,:)         ,  & ! Albedo of snow-ice system
    Tsurf(:,:,:)       ,  & ! Surface temperature                      [C]
    T1 (:,:,:)         ,  & ! Temperature upper layer                  [C]
    T2 (:,:,:)         ,  & ! Temperature lower layer                  [C]
    E1(:,:,:)          ,  & ! Energy content upper layer               [Jm/kg]
    E2(:,:,:)          ,  & ! Energy content lower layer               [Jm/kg]
    hi(:,:,:)          ,  & ! Ice thickness                            [m]
    hs(:,:,:)          ,  & ! Snow thickness                           [m]
    hiold(:,:,:)       ,  & ! Ice thickness at previous time step      [m]
    hsold(:,:,:)       ,  & ! Snow thickness at previous time step     [m]
    Qtop(:,:,:)        ,  & ! Energy flux available for surface melting[W/m2]
    Qbot(:,:,:)        ,  & ! Energy flux at ice-ocean interface       [W/m2]
    heatocei(:,:,:)    ,  & ! Energy to ocean when all ice is melted   [J]
    snow_to_ice(:,:,:) ,  & ! amount of snow that is transformed to ice[m]
    surfmelt(:,:,:)    ,  & ! surface melt water running into ocean    [m]
    surfmeltT(:,:,:)   ,  & ! Mean temperature of surface melt water   [C]
    evapwi(:,:,:)      ,  & ! amount of evaporated water if no ice left[kg/m2]
    conc(:,:,:)             ! ice concentration in each ice class
    
  REAL(wp), ALLOCATABLE :: &
    u(:,:)          ,  & ! Zonal velocity                            [m/s]
    v(:,:)          ,  & ! Meridional velocity                       [m/s]
    concSum(:,:)    ,  & ! Total ice concentration within a grid cell    
    newice(:,:)     ,  & ! New ice growth in open water              [m]
    zUnderIce(:,:)       ! water in upper ocean grid cell below ice  [m]

   INTEGER ::  kice = 1   ! Number of ice-thickness classes

  REAL(wp), ALLOCATABLE ::  hi_lim(:)   ! Thickness limits 

END TYPE t_sea_ice

! global type variables
TYPE(t_sea_ice),PUBLIC, SAVE, TARGET :: v_sea_ice


CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Constructor of sea-ice model, allocates all components and assigns zero. 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE construct_sea_ice(ppatch, p_ice, i_no_ice_thick_class)
  TYPE(t_patch),     INTENT(in)     :: ppatch
  TYPE (t_sea_ice),  INTENT (INOUT) :: p_ice
  INTEGER, INTENT(IN)               :: i_no_ice_thick_class

  !Local variables
  !INTEGER i

  INTEGER :: nblks_c, ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_sea_ice'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  nblks_c = ppatch%nblks_c

  p_ice%kice = i_no_ice_thick_class

  ALLOCATE(p_ice%isice(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for isice failed')
  END IF

  ALLOCATE(p_ice%alb(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for alb failed')
  END IF

 ALLOCATE(p_ice%Tsurf(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for Tsurf failed')
  END IF

 ALLOCATE(p_ice%T1(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for T1 failed')
  END IF

  ALLOCATE(p_ice%T2(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for T2 failed')
  END IF

 ALLOCATE(p_ice%E1(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for E1 failed')
  END IF

  ALLOCATE(p_ice%E2(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for E2 failed')
  END IF

  ALLOCATE(p_ice%hi(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hi failed')
  END IF

  ALLOCATE(p_ice%hs(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hs failed')
  END IF

  ALLOCATE(p_ice%hiold(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hiold failed')
  END IF

  ALLOCATE(p_ice%hsold(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hsold failed')
  END IF

  ALLOCATE(p_ice%Qtop(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for Qtop failed')
  END IF

  ALLOCATE(p_ice%Qbot(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for Qbot failed')
  END IF

  ALLOCATE(p_ice%heatocei(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for heatocei failed')
  END IF

  ALLOCATE(p_ice%snow_to_ice(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for snow_to_ice failed')
  END IF

  ALLOCATE(p_ice%surfmelt(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for surfmelt failed')
  END IF

  ALLOCATE(p_ice%surfmeltT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for surfmeltT failed')
  END IF

  ALLOCATE(p_ice%evapwi(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for evapwi failed')
  END IF

  ALLOCATE(p_ice%conc(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for conc failed')
  END IF

  ALLOCATE(p_ice%u(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for u failed')
  END IF
    
  ALLOCATE(p_ice%v(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for v failed')
  END IF

  ALLOCATE(p_ice%concSum(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for concSum failed')
  END IF

  ALLOCATE(p_ice%newice(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for newice failed')
  END IF

  ALLOCATE(p_ice%zUnderIce(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for zUnderIce failed')
  END IF

  ALLOCATE(p_ice%hi_lim(i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hi_lim failed')
  END IF


  p_ice%alb        = 0.0_wp       
  p_ice%Tsurf      = 0.0_wp       
  p_ice%T1         = 0.0_wp        
  p_ice%T2         = 0.0_wp      
  p_ice%E1         = 0.0_wp        
  p_ice%E2         = 0.0_wp          
  p_ice%hi         = 0.0_wp          
  p_ice%hs         = 0.0_wp        
  p_ice%hiold      = 0.0_wp     
  p_ice%hsold      = 0.0_wp       
  p_ice%Qtop       = 0.0_wp       
  p_ice%Qbot       = 0.0_wp       
  p_ice%heatocei   = 0.0_wp    
  p_ice%snow_to_ice= 0.0_wp  
  p_ice%surfmelt   = 0.0_wp    
  p_ice%surfmeltT  = 0.0_wp   
  p_ice%evapwi     = 0.0_wp      
  p_ice%conc       = 0.0_wp         
    
  
 
  p_ice%snow_to_ice= 0.0_wp 
  p_ice%surfmelt   = 0.0_wp 
  p_ice%surfmeltT  = 0.0_wp 
  p_ice%evapwi     = 0.0_wp 
  p_ice%conc       = 0.0_wp 
  p_ice%u          = 0.0_wp 
  p_ice%v          = 0.0_wp 
  p_ice%concSum    = 0.0_wp 
  p_ice%newice     = 0.0_wp 
  p_ice%zUnderIce  = 0.0_wp 

  IF(p_ice%kice==1)THEN
    p_ice%hi_lim = 0.0_wp
  ELSEIF(p_ice%kice==8)THEN
    p_ice%hi_lim(:)=(/ 0.0_wp, 0.1_wp, 0.3_wp, 0.7_wp, 1.1_wp, 1.5_wp, 2.0_wp, 2.5_wp /)
  ENDIF

  CALL message(TRIM(routine), 'end' )

END SUBROUTINE construct_sea_ice
!-------------------------------------------------------------------------
  !
  !> Destructor of sea-ice model, deallocates all components. 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE destruct_sea_ice(p_ice)
  TYPE (t_sea_ice),  INTENT (INOUT)    :: p_ice
  !Local variables
  INTEGER :: ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_sea_ice'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )
  DEALLOCATE(p_ice%alb, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for alb failed')
  END IF

 DEALLOCATE(p_ice%Tsurf, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for Tsurf failed')
  END IF

 DEALLOCATE(p_ice%T1, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for T1 failed')
  END IF

  DEALLOCATE(p_ice%T2, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for T2 failed')
  END IF

 DEALLOCATE(p_ice%E1, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for E1 failed')
  END IF

  DEALLOCATE(p_ice%E2, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for E2 failed')
  END IF

  DEALLOCATE(p_ice%hi, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hi failed')
  END IF

  DEALLOCATE(p_ice%hs, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hs failed')
  END IF

  DEALLOCATE(p_ice%hiold, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hiold failed')
  END IF

  DEALLOCATE(p_ice%hsold, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hsold failed')
  END IF

  DEALLOCATE(p_ice%Qtop, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for Qtop failed')
  END IF

  DEALLOCATE(p_ice%Qbot, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for Qbot failed')
  END IF

  DEALLOCATE(p_ice%heatocei, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for heatocei failed')
  END IF

  DEALLOCATE(p_ice%snow_to_ice, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for snow_to_ice failed')
  END IF

  DEALLOCATE(p_ice%surfmelt, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for surfmelt failed')
  END IF

  DEALLOCATE(p_ice%surfmeltT, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for surfmeltT failed')
  END IF

  DEALLOCATE(p_ice%evapwi, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for evapwi failed')
  END IF

  DEALLOCATE(p_ice%conc, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for conc failed')
  END IF

  DEALLOCATE(p_ice%u, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for u failed')
  END IF
    
  DEALLOCATE(p_ice%v, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'adellocation for v failed')
  END IF

  DEALLOCATE(p_ice%concSum, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for concSum failed')
  END IF

  DEALLOCATE(p_ice%newice, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for newice failed')
  END IF

  DEALLOCATE(p_ice%zUnderIce, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for zUnderIce failed')
  END IF

  CALL message(TRIM(routine), 'end' )
 
END SUBROUTINE destruct_sea_ice
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE construct_sfcflx(ppatch, p_sfc_flx)
  !
  TYPE(t_patch), INTENT(in) :: ppatch
  TYPE(t_sfc_flx)           :: p_sfc_flx

  ! Local variables
  INTEGER :: nblks_c, ist, jc, jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_sfcflx'

  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  rl_start_c = 1
  rl_end_c = min_rlcell

  i_startblk_c = ppatch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)


  nblks_c = ppatch%nblks_c

  ALLOCATE(p_sfc_flx%forc_wind_u(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing wind u failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_wind_v(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing wind v failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_hflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing heat flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_fwfx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing freshwater flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_swflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for short wave flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_lwflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for long wave flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_ssflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for surface sensible heat flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_slflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for surface latent heat flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_prflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for precipitation flux failed')
  END IF
  ALLOCATE(p_sfc_flx%forc_evflx(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for evaporation flux failed')
  END IF
  IF(no_tracer>=1)THEN
    ALLOCATE(p_sfc_flx%forc_tracer(nproma,nblks_c, no_tracer), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for tracer forcing failed')
    END IF

    ALLOCATE(p_sfc_flx%forc_tracer_relax(nproma,nblks_c, no_tracer), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for tracer relaxation forcing failed')
    END IF
  ENDIF

  ALLOCATE(p_sfc_flx%forc_wind_cc(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing wind_cc  failed')
  END IF

  p_sfc_flx%forc_wind_u   = 0.0_wp
  p_sfc_flx%forc_wind_v   = 0.0_wp
 
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
    ENDDO
  END DO
  IF(no_tracer>=1)THEN
    p_sfc_flx%forc_hflx         = 0.0_wp
    p_sfc_flx%forc_fwfx         = 0.0_wp
    p_sfc_flx%forc_swflx        = 0.0_wp
    p_sfc_flx%forc_lwflx        = 0.0_wp
    p_sfc_flx%forc_ssflx        = 0.0_wp
    p_sfc_flx%forc_slflx        = 0.0_wp
    p_sfc_flx%forc_prflx        = 0.0_wp
    p_sfc_flx%forc_evflx        = 0.0_wp
    p_sfc_flx%forc_tracer       = 0.0_wp
    p_sfc_flx%forc_tracer_relax = 0.0_wp
  ENDIF

  CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_sfcflx
  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor surface flux forcing for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE destruct_sfcflx(p_sfc_flx)
  TYPE(t_sfc_flx) :: p_sfc_flx
  !
  ! Local variables

  INTEGER :: ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_sfcflx'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  DEALLOCATE(p_sfc_flx%forc_wind_u, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for forcing wind u failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_wind_v, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for forcing wind v failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_hflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for forcing heat flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_fwfx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for forcing freshwater flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_swflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for heat flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_lwflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for heat flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_ssflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for heat flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_slflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for heat flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_prflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for precip flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_evflx, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for evap flux failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_tracer, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for tracer forcing failed')
  END IF
  DEALLOCATE(p_sfc_flx%forc_tracer_relax, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for tracer relaxation failed')
  END IF
   DEALLOCATE(p_sfc_flx%forc_wind_cc, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for forcing wind cc failed')
   END IF
  CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_sfcflx
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmospheric reprsentation  in ocean.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07)
  !
  SUBROUTINE construct_atmos_for_ocean(ppatch, p_as)
  !
  TYPE(t_patch), INTENT(in):: ppatch
  TYPE(t_atmos_for_ocean ) :: p_as

  ! Local variables
  INTEGER :: nblks_c, ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_for_ocean'

  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  nblks_c = ppatch%nblks_c
 
   ALLOCATE(p_as%tafo(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for tafo failed')
   END IF
   ALLOCATE(p_as%ftdew(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for ftdew failed')
   END IF
    ALLOCATE(p_as%fclou(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fclou failed')
    END IF

    ALLOCATE(p_as%fu10(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fu10 failed')
    END IF

    ALLOCATE(p_as%fswr(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fswr failed')
    END IF

    ALLOCATE(p_as%pao(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for pao failed')
    END IF

    ALLOCATE(p_as%u(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for u failed')
    END IF
    ALLOCATE(p_as%v(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for v failed')
    END IF


    p_as%tafo  = 0.0_wp
    p_as%ftdew = 0.0_wp
    p_as%fclou = 0.0_wp
    p_as%fu10  = 0.0_wp
    p_as%fswr  = 0.0_wp
    p_as%pao   = 0.0_wp
    p_as%u     = 0.0_wp
    p_as%v     = 0.0_wp

  CALL message(TRIM(routine), 'end')

  END SUBROUTINE construct_atmos_for_ocean
  !-------------------------------------------------------------------------
  !
  !>
  !!  Destructor of atmospheric reprsentation  in ocean.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE destruct_atmos_for_ocean(p_as)
  !
  TYPE(t_atmos_for_ocean ) :: p_as

  ! Local variables
  INTEGER :: ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_for_ocean'

  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

 
   DEALLOCATE(p_as%tafo, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for tafo failed')
   END IF
   DEALLOCATE(p_as%ftdew, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for ftdew failed')
   END IF
    DEALLOCATE(p_as%fclou, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fclou failed')
    END IF

    DEALLOCATE(p_as%fu10, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fu10 failed')
    END IF

    DEALLOCATE(p_as%fswr, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fswr failed')
    END IF

    DEALLOCATE(p_as%pao, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for pao failed')
    END IF

    DEALLOCATE(p_as%u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for u failed')
    END IF
    DEALLOCATE(p_as%v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for v failed')
    END IF

  CALL message(TRIM(routine), 'end')

  END SUBROUTINE destruct_atmos_for_ocean
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmos fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE construct_atmos_fluxes(ppatch, p_atm_f, i_no_ice_thick_class)
  !
  TYPE(t_patch), INTENT(in)   :: ppatch
  TYPE(t_atmos_fluxes )       :: p_atm_f
  INTEGER, INTENT(IN)         :: i_no_ice_thick_class
  ! Local variables
  INTEGER :: nblks_c, ist

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_fluxes'

  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  nblks_c = ppatch%nblks_c
 
   ALLOCATE(p_atm_f%sens(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sens failed')
   END IF

   ALLOCATE(p_atm_f%lat(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for lat failed')
   END IF

   ALLOCATE(p_atm_f%LWout(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWout failed')
   END IF

   ALLOCATE(p_atm_f%LWnet(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWnet failed')
   END IF

   ALLOCATE(p_atm_f%bot(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sens failed')
   END IF

   ALLOCATE(p_atm_f%dsensdT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for dsensdT failed')
   END IF

   ALLOCATE(p_atm_f%dlatdT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sens failed')
   END IF

   ALLOCATE(p_atm_f%dLWdT(nproma,i_no_ice_thick_class,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for dLWdT failed')
   END IF

   ALLOCATE(p_atm_f%rprecw(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for rprecw failed')
   END IF

   ALLOCATE(p_atm_f%rpreci(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for rpreci failed')
   END IF

   ALLOCATE(p_atm_f%sensw(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sensw failed')
   END IF

   ALLOCATE(p_atm_f%latw(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for latw failed')
   END IF


   ALLOCATE(p_atm_f%LWoutw(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWoutw failed')
   END IF

   ALLOCATE(p_atm_f%LWnetw(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWnetw failed')
   END IF

   ALLOCATE(p_atm_f%SWin(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for SWin failed')
   END IF

   ALLOCATE(p_atm_f%LWin(nproma,nblks_c), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWin failed')
   END IF

  CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_atmos_fluxes
  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor of atmos fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE destruct_atmos_fluxes(p_atm_f)
  !
  TYPE(t_atmos_fluxes )       :: p_atm_f
  ! Local variables
  INTEGER :: ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_fluxes'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

 
   DEALLOCATE(p_atm_f%sens, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for sens failed')
   END IF

   DEALLOCATE(p_atm_f%lat, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for lat failed')
   END IF

   DEALLOCATE(p_atm_f%LWout, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for LWout failed')
   END IF

   DEALLOCATE(p_atm_f%LWnet, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for LWnet failed')
   END IF

   DEALLOCATE(p_atm_f%bot, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for sens failed')
   END IF

   DEALLOCATE(p_atm_f%dsensdT, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for dsensdT failed')
   END IF

   DEALLOCATE(p_atm_f%dlatdT, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for sens failed')
   END IF

   DEALLOCATE(p_atm_f%dLWdT, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for dLWdT failed')
   END IF

   DEALLOCATE(p_atm_f%rprecw, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for rprecw failed')
   END IF

   DEALLOCATE(p_atm_f%rpreci, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for rpreci failed')
   END IF

   DEALLOCATE(p_atm_f%sensw, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for sensw failed')
   END IF

   DEALLOCATE(p_atm_f%latw, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for latw failed')
   END IF


   DEALLOCATE(p_atm_f%LWoutw, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for LWoutw failed')
   END IF

   DEALLOCATE(p_atm_f%LWnetw, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for LWnetw failed')
   END IF

   DEALLOCATE(p_atm_f%SWin, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for SWin failed')
   END IF

   DEALLOCATE(p_atm_f%LWin, STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for LWin failed')
   END IF

  CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_atmos_fluxes

!-------------------------------------------------------------------------
!
!> ice_init
!!
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE ice_init( ppatch, p_os, ice) !, Qatm, QatmAve)
  TYPE(t_patch), INTENT(in)             :: ppatch 
  TYPE(t_hydro_ocean_state)             :: p_os
  TYPE (t_sea_ice),      INTENT (INOUT) :: ice
!  TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
!  TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve
  
!local variables
  REAL(wp):: Tinterface(nproma,i_no_ice_thick_class, ppatch%nblks_c) ! temperature at snow-ice interface
  REAL(wp):: draft(nproma,i_no_ice_thick_class, ppatch%nblks_c)     ! position of ice-ocean interface below sea level
  
!  INTEGER i,j,k      ! counter for loops
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:ice_init'
  !-------------------------------------------------------------------------

  CALL message(TRIM(routine), 'start' )

!Constructor basic init already done at this point
!   CALL alloc_mem_commo_ice (ice, Qatm, QatmAve)
!   CALL ice_zero            (ice, Qatm, QatmAve)

! FORALL(i=1:nproma, j=1:ppatch%nblks_c, k=1:i_no_ice_thick_class) 
!    ice% hi    (i,j,k) = sictho (i,j)
!    ice% hs    (i,j,k) = sicsno (i,j)
! END FORALL
  
  ice% Tsurf  = Tf
  ice% T1     = Tf
  ice% T2     = Tf
  ice% conc   = 0.0_wp
  ice% isice  = .FALSE.
  Tinterface  = Tf
  draft       = 0.0_wp

! Stupid initialisation trick for Levitus initialisation
  IF (init_oce_prog == 1) THEN
    WHERE (p_os%p_prog(nold(1))%tracer(:,1,:,1) <= -1.0_wp )
      ice%hi(:,1,:) = 2._wp
      ice%conc(:,1,:) = 1._wp
    ENDWHERE
    WHERE (p_os%p_prog(nold(1))%tracer(:,:,:,1) <= -1.0_wp ) &
      &   p_os%p_prog(nold(1))%tracer(:,:,:,1) = Tf
  ENDIF

  WHERE(ice% hi(:,:,:) > 0.0_wp)
    Tinterface (:,:,:) = (Tf * (ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))+&
                        ice%Tsurf(:,:,:)) / (1.0_wp+ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))
    ice% conc  (:,:,:) = 1.0_wp/REAL(i_no_ice_thick_class,wp)
    ice% isice (:,:,:) = .TRUE.
    ice% T1    (:,:,:) = Tf + 2._wp/3._wp*(Tinterface(:,:,:)-Tf)
    ice% T2    (:,:,:) = Tf + 1._wp/3._wp*(Tinterface(:,:,:)-Tf)
    draft      (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rhow
  END WHERE
  
!   ice%zUnderIce (:,:)   = dzw(1) + zo (:,:) &
!                           - sum(draft(:,:,:) * ice%conc(:,:,:),2)
   ice%zUnderIce (:,:)   = v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) &
                           - sum(draft(:,:,:) * ice%conc(:,:,:),2)

  CALL message(TRIM(routine), 'end' )

END SUBROUTINE ice_init
!-------------------------------------------------------------------------  
!
!  
!>
!! !  ice_fast: Ice routines for atmospheric time step. Sets air-ice fluxes and
!!    calculates the development of the ice temperature field
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
 SUBROUTINE ice_fast(ppatch, ice,Qatm,QatmAve)
   TYPE(t_patch),            INTENT(IN)     :: ppatch 
!   TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
!   TYPE(t_atmos_for_ocean),  INTENT(IN)     :: p_as
   TYPE (t_sea_ice),         INTENT (INOUT) :: ice
   TYPE (t_atmos_fluxes),    INTENT (INOUT) :: Qatm
   TYPE (t_atmos_fluxes),    INTENT (INOUT) :: QatmAve
!------------------------------------------------------------------------- 
!    CALL get_atmos_fluxes (ppatch, p_os,p_as,ice, Qatm)
    CALL set_ice_albedo   (ppatch,ice)
    CALL set_ice_temp     (ppatch,ice, Qatm)
    CALL sum_fluxes       (Qatm, QatmAve)
 END SUBROUTINE ice_fast
!-------------------------------------------------------------------------------
!
!  
!>
!! !  ice_slow: Ice routines for oceand time step. Calculates average of atmospheric
! !           time steps, ice velocity, ice growth rates and updates ice structure
! !           accordingly
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
 SUBROUTINE ice_slow(ppatch, p_os,ice, QatmAve, p_sfc_flx)  
   TYPE(t_patch),         INTENT(IN)      :: ppatch 
   TYPE(t_hydro_ocean_state),INTENT(INOUT):: p_os
!   TYPE(t_atmos_for_ocean),INTENT(IN)     :: p_as
   TYPE (t_sea_ice),      INTENT (INOUT)  :: ice
!   TYPE (t_atmos_fluxes), INTENT (INOUT)  :: Qatm
   TYPE (t_atmos_fluxes), INTENT (INOUT)  :: QatmAve
   TYPE(t_sfc_flx),        INTENT (INOUT) :: p_sfc_flx
!-------------------------------------------------------------------------------

   CALL ave_fluxes     (ice, QatmAve)
!   !CALL ice_dynamics   (ice, QatmAve)
   CALL ice_growth     (ppatch,p_os,ice, QatmAve%rpreci)!, QatmAve%lat)
   CALL upper_ocean_TS (ppatch,p_os,ice, QatmAve, p_sfc_flx)
   CALL new_ice_growth (ppatch,ice, p_os,p_sfc_flx)
!   CALL ice_advection  (ice)
!   CALL write_ice(ice,QatmAve,1,ie,je)
   CALL ice_zero       (ice, QatmAve)
!   sictho = ice%hi   (:,:,1) * ice%conc (:,:,1)
!   sicomo = ice%conc (:,:,1)
!   sicsno = ice%hs   (:,:,1) * ice%conc (:,:,1)
 END SUBROUTINE ice_slow
!-------------------------------------------------------------------------  
!
!  
!>
!! !  get_atmos_fluxes: Sets the atmospheric fluxes for the update of the ice
! !                 temperature
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE get_atmos_fluxes (ppatch, p_os,p_as,ice, Qatm)
  TYPE(t_patch),            INTENT(IN)     :: ppatch 
  TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
  TYPE(t_atmos_for_ocean),  INTENT(IN)     :: p_as
   TYPE (t_sea_ice),        INTENT (INOUT) :: ice
   TYPE (t_atmos_fluxes),   INTENT (INOUT) :: Qatm
 

#ifdef coupled
   Qatm% SWin   = 
   Qatm% LWin   =
   Qatm% sens   = 
   Qatm% lat    =
   Qatm% dsensdT = 
   Qatm% dlatdT  =
   Qatm% dLWdT   =
#elif defined CORE
   CALL budget_core   (ice, Qatm)
#else
   CALL calc_atm_fluxes_from_bulk(ppatch, p_as, p_os, ice, Qatm)
#endif

END SUBROUTINE get_atmos_fluxes 
!-------------------------------------------------------------------------  
!
!  
!>
!! !   sum_fluxes: adds atmospheric fluxes for ocean time stepping. Necessary for
!!      diagnosis, not for the ice model itself.
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE sum_fluxes        (Qatm, QatmAve)
  TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
  TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve
  QatmAve % sens   (:,:,:) = QatmAve % sens    + Qatm % sens
  QatmAve % sensw  (:,:)   = QatmAve % sensw   + Qatm % sensw
  QatmAve % lat    (:,:,:) = QatmAve % lat     + Qatm % lat
  QatmAve % latw   (:,:)   = QatmAve % latw    + Qatm % latw
  QatmAve % LWout  (:,:,:) = QatmAve % LWout   + Qatm % LWout
  QatmAve % LWoutw (:,:)   = QatmAve % LWoutw  + Qatm % LWoutw
  QatmAve % LWnet  (:,:,:) = QatmAve % LWnet   + Qatm % LWnet
  QatmAve % LWnetw (:,:)   = QatmAve % LWnetw  + Qatm % LWnetw
  QatmAve % SWin   (:,:)   = QatmAve % SWin    + Qatm % SWin
  QatmAve % LWin   (:,:)   = QatmAve % LWin    + Qatm % LWin
  QatmAve % rprecw (:,:)   = QatmAve % rprecw  + Qatm % rprecw
  QatmAve % rpreci (:,:)   = QatmAve % rpreci  + Qatm % rpreci
  QatmAve % counter        = QatmAve % counter + 1 
END SUBROUTINE sum_fluxes  
!-------------------------------------------------------------------------------
!
!  
!>
!! ! ave_fluxes: calculates the average of the atmospheric fluxes for ocean time  
!!   sum_fluxes: adds atmospheric fluxes for ocean time stepping. Necessary for
!!   diagnosis, not for the ice model itself.
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE ave_fluxes (ice, QatmAve)
  TYPE (t_sea_ice),      INTENT (INOUT) :: ice
  TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve
  !
  !Local variables
  REAL(wp) :: ctr
!-------------------------------------------------------------------------------
  ctr = REAL(QatmAve% counter,wp)
  QatmAve% sens   (:,:,:) = QatmAve% sens   / ctr
  QatmAve% sensw  (:,:)   = QatmAve% sensw  / ctr
  QatmAve% lat    (:,:,:) = QatmAve% lat    / ctr
  QatmAve% latw   (:,:)   = QatmAve% latw   / ctr
  QatmAve% LWout  (:,:,:) = QatmAve% LWout  / ctr
  QatmAve% LWoutw (:,:)   = QatmAve% LWoutw / ctr
  QatmAve% LWnet  (:,:,:) = QatmAve% LWnet  / ctr
  QatmAve% LWnetw (:,:)   = QatmAve% LWnetw / ctr
  QatmAve% SWin   (:,:)   = QatmAve% SWin   / ctr
  QatmAve% LWin   (:,:)   = QatmAve% LWin   / ctr
  QatmAve% rprecw (:,:)   = QatmAve% rprecw / ctr
  QatmAve% rpreci (:,:)   = QatmAve% rpreci / ctr
  ice    % Qbot   (:,:,:) = ice    % Qbot   / ctr
  ice    % Qtop   (:,:,:) = ice    % Qtop   / ctr

END SUBROUTINE ave_fluxes
!-------------------------------------------------------------------------------
!
!  
!>
!! ! ice_zero: set the avereged fluxes to zero
!! @par Revision History
!! Initial release by Einar Olason, MPI-M (2011-09). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE ice_zero (ice,QatmAve)
  TYPE (t_sea_ice),      INTENT (INOUT) :: ice
!  TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
  TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

!  Qatm    % sens        (:,:,:) = 0._wp
!  Qatm    % sensw       (:,:)   = 0._wp
!  Qatm    % lat         (:,:,:) = 0._wp
!  Qatm    % latw        (:,:)   = 0._wp
!  Qatm    % LWout       (:,:,:) = 0._wp
!  Qatm    % LWoutw      (:,:)   = 0._wp
!  Qatm    % LWnet       (:,:,:) = 0._wp
!  Qatm    % LWnetw      (:,:)   = 0._wp
!  Qatm    % SWin        (:,:)   = 0._wp
!  Qatm    % LWin        (:,:)   = 0._wp
!  Qatm    % rprecw      (:,:)   = 0._wp
!  Qatm    % rpreci      (:,:)   = 0._wp
                        
  QatmAve % sens        (:,:,:) = 0._wp
  QatmAve % sensw       (:,:)   = 0._wp
  QatmAve % lat         (:,:,:) = 0._wp
  QatmAve % latw        (:,:)   = 0._wp
  QatmAve % LWout       (:,:,:) = 0._wp
  QatmAve % LWoutw      (:,:)   = 0._wp
  QatmAve % LWnet       (:,:,:) = 0._wp
  QatmAve % LWnetw      (:,:)   = 0._wp
  QatmAve % SWin        (:,:)   = 0._wp
  QatmAve % LWin        (:,:)   = 0._wp
  QatmAve % rprecw      (:,:)   = 0._wp
  QatmAve % rpreci      (:,:)   = 0._wp
  QatmAve % counter             = 0 

  ice     % Qbot        (:,:,:) = 0._wp
  ice     % Qtop        (:,:,:) = 0._wp
  ice     % surfmelt    (:,:,:) = 0._wp
  ice     % surfmeltT   (:,:,:) = 0._wp
  ice     % evapwi      (:,:,:) = 0._wp
  ice     % hiold       (:,:,:) = 0._wp
  ice     % snow_to_ice (:,:,:) = 0._wp
  ice     % heatOceI    (:,:,:) = 0._wp

  END SUBROUTINE ice_zero

!-------------------------------------------------------------------------------
!
!  
!>
!! ! ice_albedo: set ice albedo 
!-------------------------------------------------------------------------------
!
!  
!>
!! ! ice_albedo: set ice albedo 
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE set_ice_albedo(ppatch, ice) 
  TYPE(t_patch),    INTENT(IN)     :: ppatch 
  TYPE (t_sea_ice), INTENT (INOUT) :: ice
  !
  !Local variables
  REAL(wp), PARAMETER :: albtrans   = 0.5_wp
  REAL(wp)            :: albflag(nproma,i_no_ice_thick_class, ppatch%nblks_c)
  !-------------------------------------------------------------------------------

  ! This is Uwe's albedo expression from the old budget function
  albflag (:,:,:) =  1.0_wp/ ( 1.0_wp+albtrans * (ice%tsurf(:,:,:))**2 )
  
  WHERE (ice  % isice)
    WHERE (ice % hs > 1.e-2_wp)
      ice% alb(:,:,:) =  albflag(:,:,:) * albsm + (1.0_wp-albflag(:,:,:)) * albs
    ELSEWHERE
      ice% alb(:,:,:) =  albflag(:,:,:) * albim + (1.0_wp-albflag(:,:,:)) * albi
    END WHERE
  END WHERE
END SUBROUTINE set_ice_albedo
!-------------------------------------------------------------------------------
!
!  
!>
!! ! set_ice_temp:: calculate new ice + snow temperatures according to sec.2a from
!!           Winton, M., 2000: A Reformulated Three-Layer Sea Ice Model,   
!!           J. Atmos. Oce. Tech., 17, 525-531. 
!!
!!           doi: 10.1175/1520-0426(2000)017<0525:ARTLSI> (put into google)
!!
!! This function changes:
!! ice % Ts       the new surface temperature   for each ice category     [�C]
!! ice % T1       the new upper ice+snow temp.  for each ice category     [�C]
!! ice % T2       the new lower ice temperature for each ice category     [�C]
!! ice % Qbot     Heat flux available for freezing/melting at ice bottom  [W/m�]
!! ice % Qtop     Heat flux available for melting at ice surface          [W/m�]
!!
!!           all "dtime" in this function are atmospheric time step
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE set_ice_temp(ppatch,ice, Qatm) 
  TYPE(t_patch),         INTENT(IN)     :: ppatch 
  TYPE (t_sea_ice),      INTENT (INOUT) :: ice
  TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm

  !!Local variables
  REAL(wp), DIMENSION (nproma,i_no_ice_thick_class, ppatch%nblks_c) ::           &
    A,      & ! Eq. 7
    A1,     & ! Eq. 16
    A1a,    & ! First two terms of Eq. 16 and 19
    B,      & ! Eq. 8
    B1,     & ! Eq. 17
    B1a,    & ! First three terms of Eq. 17 and 20
    C1,     & ! Eq. 18
    D,      & ! 1./(6*dT*K2) + rhoi*hi*C for Eq. 16, 17, 19, 20
    iK1B,   & ! 1./(K1 + B) (used in eq. 16 and 17)
    K1,     & ! Winton's K 1/2 (eq. 5)
    K2,     & ! Winton's K 3/2 (eq. 10)
    SWin3D, & ! Short-wave radiation field splitted into ice categories
    Tsurfm, & ! Surface melting temperature
  
  ! The following fields are copied from 'ice' for easier readability
    hi,     & ! Ice thickness                             [m] 
    hs,     & ! Snow thickness                            [m]
    Tsurf,  & ! Surface temperature                       [�C]
    T1,     & ! Temperature of the upper ice+snow layer   [�C]
    T2        ! Temperature of the lower layer            [�C]

  REAL(wp) ::   &
    idt2      ! 1 / (2*dt)
  
  INTEGER   &
    i,j,k         ! counter for loops
 !-------------------------------------------------------------------------------
  ! Copy fields from 'ice' structure into local variables for easy readability
  hi   (:,:,:) = ice% hi   (:,:,:)
  hs   (:,:,:) = ice% hs   (:,:,:)
  Tsurf(:,:,:) = ice% Tsurf(:,:,:)
  T1   (:,:,:) = ice% T1   (:,:,:)
  T2   (:,:,:) = ice% T2   (:,:,:)

  idt2   =  1.0_wp / (2.0_wp*dtime)

  ! Create array of shortwave radiation split up into ice categories
  ! (purely for computational reasons)
  FORALL(i=1:nproma, j=1:ppatch%nblks_c, k=1:i_no_ice_thick_class, ice % isice (i,k,j)) &
      SWin3d(i,k,j) = Qatm% SWin(i,j)

  ! Calculate new ice temperature wherever there is ice 
  ! lat > 0, sens > 0, LWnet >0 , SWin > 0  for downward flux 
  ! dlatdT, dsensdT, dLWdT >0 for downward flux increasing with increasing Tsurf
  !

  WHERE (ice % isice (:,:,:) )
    B   (:,:,:) = -Qatm% dlatdT - Qatm% dsensdT - Qatm% dLWdT             ! Eq.  8
    A   (:,:,:) = -Qatm% lat - Qatm% sens - Qatm% LWnet -                     &
                  (1.0_wp - ice% alb) * I_0 *  SWin3d  - Tsurf* B        ! Eq.  7
    K1  (:,:,:)  =  4.0_wp * ki * ks / (ks * hi + 4.0_wp * ki * hs)          ! Eq.  5
    K2  (:,:,:)  =  2.0_wp * ki / hi                                     ! Eq. 10
    D   (:,:,:)  =  1.0_wp / (6.0_wp * dtime * K2 + rhoi*hi*ci)                 
    iK1B(:,:,:)  =  1.0_wp / (K1 + B)

   ! Set temperature at which surface is fully liquid
    WHERE (hs(:,:,:) > 1e-6_wp) 
      Tsurfm(:,:,:)  =  0.0_wp
    ELSEWHERE
      Tsurfm(:,:,:)  =  - muS
    END WHERE

    
    A1a   (:,:,:)  =  rhoi*hi * idt2 * ci + K2* (4.0_wp * dtime * K2 + rhoi*hi*ci)*D 
    A1    (:,:,:)  =  A1a + K1*B * iK1B                              ! Eq. 16
    B1a   (:,:,:)  =  -rhoi*hi* (ci*T1 - Lfreez*muS/T1) * idt2 - I_0              & !originally L instead of Lfreez
                       - K2*(4.0_wp*dtime*K2*Tf+rhoi*hi*ci*T2)*D
    B1    (:,:,:)  =  B1a + A*K1*iK1B                                ! Eq. 17
    C1    (:,:,:)  =  - rhoi*hi * Lfreez * muS * idt2                     ! Eq. 18  !originally L instead of Lfreez
    T1    (:,:,:)  =  -(B1 + SQRT(B1*B1-4.0_wp*A1*C1)) / (2.0_wp*A1)           ! Eq. 21
    Tsurf (:,:,:)  =  (K1*T1-A) * iK1B                               ! Eq.  6


    WHERE ( Tsurf(:,:,:) > Tsurfm(:,:,:) ) 
      A1           (:,:,:)  =  A1a + K1                              ! Eq. 19
      B1           (:,:,:)  =  B1a - K1*Tsurfm                       ! Eq. 20
      T1           (:,:,:)  =  -(B1 + SQRT(B1*B1-4.0_wp*A1*C1)) / (2.0_wp*A1)  ! Eq. 21
      Tsurf        (:,:,:)  =  Tsurfm                               
      ! Sum up heatfluxes available for melting at ice surface for each
      ! atmopheric time step. ice%Qtop will be averaged in ave_fluxes
      ice%Qtop     (:,:,:)  =  ice% Qtop + K1*(T1-Tsurf) - (A + B*Tsurf)! Eq. 22
    END WHERE
   
   
    T2         (:,:,:)  =  ( 2.0_wp*dtime*K2*(T1+2.0_wp*Tf) + rhoi*hi*ci*T2) * D  ! Eq. 15
    ! Sum up conductive heatflux at ice-ocean interface for each atmospheric
    ! time step. ice%Qtop will be averaged in ave_fluxes
    ! The ocean heat flux is calculated at the beginning of ice_growth
    ice% Qbot  (:,:,:)  =  ice% Qbot - 4.0_wp*Ki*(Tf-T2)/hi               ! Eq. 23

    ! Save new temperatures in 'ice' structure
    ice% Tsurf (:,:,:)  =  Tsurf
    ice% T1    (:,:,:)  =  T1
    ice% T2    (:,:,:)  =  T2
  END WHERE

  ipl_src=1  ! output print level (1-5, fix)
  CALL print_mxmn('ice%Tsurf',1,ice%Tsurf(:,1,:),1,ppatch%nblks_c,'ice',ipl_src)

END SUBROUTINE set_ice_temp
!-------------------------------------------------------------------------------
!
!  
!>
!! ! set_ice_temp:: ice_growth - change ice and snow thickness (Winton 2000, section 2b)
!! This function changes:
!! ice % hs       new snow thickness for each ice category                [m]
!! ice % hi       new ice  thickness for each ice category                [m]
!! ice % hsold    old snow thickness for each ice category                [m]
!! ice % hiold    old ice  thickness for each ice category                [m]
!! ice % T1       the new upper ice+snow temp.  for each ice category     [�C]
!! ice % T2       the new lower ice temperature for each ice category     [�C]
!! ice % evapwi   amount of evaporated water from the mixed layer
!!                in previously ice covered areas if all ice is gone      [kg/m�]
!! ice % heatOceI to contain the energy that is available to the mixed layer
!!                in previously ice covered areas if all ice is gone      [J]
!!
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE ice_growth(ppatch, p_os, ice, rpreci)!, lat)
 TYPE(t_patch),    INTENT(IN)     :: ppatch 
  TYPE(t_hydro_ocean_state),INTENT(IN):: p_os
  TYPE (t_sea_ice),INTENT (INOUT) :: ice
  REAL(wp),        INTENT (IN )   :: rpreci(:,:) ! water equiv. solid 
                                                 ! precipitation rate [m/s] DIMENSION (ie,je),
!  REAL(wp),        INTENT (IN )   :: lat(:,:,:)  ! lat. heat flux  [W/m�] DIMENSION (ie,je,kice),

  !!Local variables
!!Local Variables
  REAL(wp), DIMENSION (nproma,i_no_ice_thick_class, ppatch%nblks_c) ::         &
    below_water, & ! Thickness of snow layer below water line           [m]
    C1,          & ! L  * rhos * hs                                     [J/m�]
    C2,          & ! E1 * rhoi * h1                                     [J/m�]
    C3,          & ! E2 * rhoi * h2                                     [J/m�]
    delh2,       & ! increase of bottom layer thickness (Eq. 24)        [m]
    draft,       & ! depth of ice-ocean interface below sea level       [m]
    E1,          & ! Energy content of upper ice+snow layer             [J/kg]
    E2,          & ! Energy content of lower ice      layer             [J/kg]
    f1,          & ! Fraction of upper ice in new ice layer (Eq. 37)
    h1,          & ! Thickness of upper ice layer                       [m]
    h2,          & ! Thickness of lower ice layer                       [m] 
    new_snow3d,  & ! New snow fallen onto each ice category             [m]
!    subli,       & ! Amount of ice+snow that is sublimated away         [kg/m�]
    Tbar,        & ! Dummy temperature for new temperature calculation  [�C]
    surfmeltsn,  & ! Surface melt water from snow melt with T=0�C       [m]
    surfmelti1,  & ! Surface melt water from upper ice with T=-muS      [m]
    surfmelti2,  & ! Surface melt water from lower ice with T=-muS      [m]
    heatocei,    & ! Oceanic heat flux                                  [W/m^2]

    ! The following fields are copied from 'ice' for easier readability
    Qbot,        & ! Energy flux available for bottom melting/freezing [W/m�]
    Qtop,        & ! Energy flux available for surface melting         [W/m�]
    hi,          & ! ice thickness                                     [m]
    hs,          & ! snow thickness                                    [m]
    Tsurf,       & ! surface temperature                               [m]
    T1,          & ! temperature of upper ice+snow layer               [�C]
    T2             ! temperature of lower ice      layer               [�C]
 
  INTEGER k
  REAL(wp),POINTER :: sst(:,:)

  delh2=0._wp
  sst =>p_os%p_prog(nold(1))%tracer(:,1,:,1)
!  !-------------------------------------------------------------------------------
  ! Calculate snow fall and create array split into ice categories
  new_snow3d (:,1,:)   = rpreci (:,:) * dtime * rhow / rhos 
  FORALL(k=2:i_no_ice_thick_class)  new_snow3d(:,k,:)  = new_snow3d(:,1,:)

  ! Add oceanic heat flux to energy available at the bottom of the ice.
  ! Currently (as in growth.f90): all energy available in upper ocean grid cell 
  ! is supplied to the ice and the upper ocean temperature is held at the 
  ! freezing point. This is not very physical.
  DO k=1,i_no_ice_thick_class
    WHERE (ice%isice(:,k,:)) 
      heatOceI(:,k,:) = ( sst - Tf ) * ice%zUnderIce * clw*rhow/dtime
      ice%Qbot(:,k,:) = ice%Qbot(:,k,:) + heatOceI(:,k,:)
    ENDWHERE
  END DO
 
  ! Do the following wherever there is ice
  !isice: &  
  WHERE (ice% isice)
    ! Copy fields from 'ice' structure for better readability of the code
    hi                (:,:,:) = ice%   hi
    hs                (:,:,:) = ice%   hs
    Qbot              (:,:,:) = ice%   Qbot
    Qtop              (:,:,:) = ice%   Qtop
    Tsurf             (:,:,:) = ice%   Tsurf
    T1                (:,:,:) = ice%   T1
    T2                (:,:,:) = ice%   T2
  
   
    ! Save ice thickness at previous time step for calculation of heat and salt
    ! flux into ocean in subroutine upper_ocean_TS
    ice % hiold (:,:,:) = hi
    ice % hsold (:,:,:) = hs

    h1(:,:,:) = hi(:,:,:) / 2.0_wp
    h2(:,:,:) = h1(:,:,:)


    ! Apply mass increasing changes first. 
    ! 1. Snow fall

    hs(:,:,:) = hs(:,:,:) + new_snow3d(:,:,:)
    
    ! 2. Bottom ice-growth  (maybe add frazil ice?)

    ! #eoo# Eqns. 24, 27--29 and 31--36 appear to be missing rhoi or rhos to get the proper units
    ! for Delta h. But these are included in this program
    WHERE (Qbot < 0.0_wp) 
      delh2(:,:,:)  = Qbot * dtime / (rhoi * (ci * (Tf + muS) - Lfreez))  ! Eq. 24 & 25
      T2   (:,:,:)  = (delh2*Tf + h2 * T2) / (delh2 + h2)                 ! Eq. 26
      h2   (:,:,:)  = h2 + delh2
    END WHERE

    ! Now mass decreasing changes. 
    ! 1. Evaporation
    ! #eoo# Not in Winton - does this count the latent fluxes twice?

!    subli(:,:,:) = lat  / Lsub * dtime;    ![kg/m�]
!    WHERE     (subli <= hs*rhos )         
!      hs(:,:,:) = hs - subli / rhos
!    ELSEWHERE (subli <= hs*rhos + h1*rhoi )          ! if all snow is gone
!      hs(:,:,:) = 0.0_wp
!      h1(:,:,:) = h1 - (subli - hs*rhos) / rhoi
!    ELSEWHERE (subli <= hs*rhos + (h1+h2)*rhoi )     ! if upper ice is gone
!      hs(:,:,:) = 0.0_wp
!      h1(:,:,:) = 0.0_wp
!      h2(:,:,:) = h2 - (subli - hs*rhos - h1*rhoi) / rhoi
!    ELSEWHERE                                         ! if all ice is gone
!      hs(:,:,:) = 0.0_wp
!      h1(:,:,:) = 0.0_wp
!      h2(:,:,:) = 0.0_wp
!      ice% evapwi(:,:,:) = (subli - hs*rhos - (h1+h2)*rhoi) * Lsub / Lvap
!    END WHERE
   
 
   ! 2. surface ablation (if any) 

   E1(:,:,:) = ci * ( T1+muS ) - Lfreez*(1.0_wp+muS/T1)   ! Eq.  1 (energy upper layer) 
   E2(:,:,:) = ci * ( T2+muS ) - Lfreez                   ! Eq. 25 (energy lower layer), originally L instead of Lfreez
   C1(:,:,:) = Lfreez  * rhos * hs
   C2(:,:,:) = E1 * rhoi * h1
   C3(:,:,:) = E2 * rhoi * h2
  
   WHERE ( Qtop(:,:,:) > 0.0_wp ) 
      surfmeltsn   (:,:,:) = MIN(Qtop*dtime / (Lfreez * rhos), hs)
      hs           (:,:,:) = hs - surfmeltsn                            ! Eq. 27
      ice%surfmelt (:,:,:) = surfmeltsn * rhos/rhow
      WHERE (hs(:,:,:) <= 0.0_wp) 
        surfmelti1   (:,:,:) = MIN((Qtop*dtime-C1) / (-E1*rhoi), h1)
        h1           (:,:,:) = h1 - surfmelti1                          ! Eq. 28
        ice%surfmelt (:,:,:) = ice%surfmelt + surfmelti1 * rhoi/rhow
        WHERE (h1(:,:,:) <= 0.0_wp) 
          surfmelti2   (:,:,:) = MIN((Qtop*dtime-C1+C2) / (-E2*rhoi), h2)
          h2           (:,:,:) = h2 - surfmelti2                        ! Eq. 29
          ice%surfmelt (:,:,:) = ice%surfmelt + surfmelti2 * rhoi/rhow
          WHERE (h2(:,:,:) <= 0.0_wp) 
           ice% heatOceI(:,:,:) = Qtop + (-C1 + C2 + C3)/dtime               ! Eq. 30
! Flux - not heat
!           ice% heatOceI(:,:,:) = Qtop*dtime - C1 + C2 + C3                ! Eq. 30
          END WHERE
        END WHERE
      END WHERE
      ! Calculate average temperature of surface melt water 
      ! T(snow) = 0�C, T(ice) = -muS �C
      ice%surfmeltT = (surfmelti1+surfmelti2) * (-muS) /  ice%surfmelt
    END WHERE
   
    C1(:,:,:) = Lfreez    * rhos * hs
    C2(:,:,:) = E1 * rhoi * h1
    C3(:,:,:) = E2 * rhoi * h2
 
   ! 3. bottom ablation (if any)

    WHERE ( Qbot(:,:,:) > 0.0_wp ) 
      h2 (:,:,:) = h2 - MIN(Qbot * dtime/ (-E2*rhoi), h2)                  ! Eq. 31
      WHERE (h2(:,:,:) <= 0.0_wp) 
        h1 (:,:,:) = h1 - MIN((Qbot * dtime  + C3) / (-E1*rhoi), h1)       ! Eq. 32
        WHERE (h1(:,:,:) <= 0.0_wp) 
          hs (:,:,:) = hs(:,:,:) - MIN((Qbot * dtime+C3+C2)/(Lfreez*rhos), hs(:,:,:))!33
          WHERE (hs (:,:,:) <= 0.0_wp) 
           ice% heatOceI(:,:,:) = ice% heatOceI + Qbot + (-C1 + C2 + C3)/dtime    !34
! Flux - not heat
!           ice% heatOceI(:,:,:) = ice% heatOceI + Qbot * dtime - C1 + C2 + C3   !34
          END WHERE
        END WHERE
      END WHERE
    END WHERE

    ! Calculate ice thickness and draft (ice+snow depth below water line)
    hi          (:,:,:) = h1 + h2
    draft       (:,:,:) = (rhoi*hi+rhos*hs) / rhow
    below_water (:,:,:) = draft - hi

    
    ! snow -> ice conversion for snow below waterlevel
    ! Currently not quite physical: Snow is pushed together to form new ice,
    ! hence snow thickness decreases more than ice thickness by rhoi/rhos
    ! ( analogue to old growth.f90 sea-ice model )
    ! Salt content of snow ice is equal to that of normal ice, salt is 
    ! removed from the ocean
    ! Temperature of new upper ice is calculated as described in the paragraph 
    ! below Eq. 36
    WHERE (below_water (:,:,:) > 0.0_wp) 
      ice% snow_to_ice  (:,:,:) = below_water * rhoi / rhos
      hs                (:,:,:) = hs - ice% snow_to_ice
      f1                (:,:,:) = h1 / (h1+below_water)
      Tbar              (:,:,:) = f1  * ( T1 - Lfreez* muS/(ci*T1) ) - (1.0_wp-f1)*muS 
      !originally L instead of Lfreez. see also below
      T1                (:,:,:) = 0.5_wp * ( Tbar - SQRT(Tbar*Tbar + 4.0_wp*muS*Lfreez/ci) )
      h1                (:,:,:) = h1 + below_water
      hi                (:,:,:) = h1 + h2
    END WHERE
 
    ! Even up upper and lower layer
    WHERE ( h1(:,:,:) < h2(:,:,:)  ) 
      f1   (:,:,:)  =  h1 / (0.5_wp*hi)                                
      Tbar (:,:,:)  =  f1 * ( T1 - Lfreez*muS/(ci*T1) ) + (1.0_wp-f1)*T2        ! Eq. 39
      T1   (:,:,:)  =  0.5_wp * ( Tbar - SQRT(Tbar*Tbar + 4.0_wp*muS*Lfreez/ci) )  ! Eq. 38
    ELSEWHERE ( h1(:,:,:) > h2(:,:,:) ) 
      f1   (:,:,:)  =  h1 / (0.5_wp*hi) - 1.0_wp
      T2   (:,:,:)  =  f1 * ( T1 - Lfreez*muS/(ci*T1) ) + (1.0_wp-f1)*T2        ! Eq. 40
    END WHERE
  
    ! T2 can get above bulk melting temperature. If this happens, use additional
    ! energy to melt equal thickness of upper and lower layer (last para.
    ! section 2)
    ! Energy available for melting: -h2 * ci * (T2+muS)
    ! Energy needed for melting lower layer: L
    ! Energy needed for melting upper layer: -(ci*(T1+muS)-L*(1+muS/T1)) (Eq. 1)
    WHERE (t2 (:,:,:) > -muS)                  
      hi (:,:,:)  =  hi - h2*ci*(T2+muS) /                          &
                     ( 0.5_wp*Lfreez - 0.5_wp*(ci*(T1+muS) - Lfreez*(1.0_wp+muS/T1)) )
      T2 (:,:,:)  =  -muS
    END WHERE

    ! Is this necessary?
    WHERE (hi(:,:,:) <= 0.0_wp) 
      T1 (:,:,:)   =  Tf
      T2 (:,:,:)   =  Tf
      ice%isice    =  .FALSE.
      ice%conc     = 0.0_wp
      hi           = 0.0_wp
    ELSEWHERE
      ice%heatOceI =  ice%heatOceI - heatOceI
    END WHERE
  
    ! Save new values in 'ice' structure
    ice%hi (:,:,:) = hi
    ice%hs (:,:,:) = hs
    ice%T1 (:,:,:) = T1
    ice%T2 (:,:,:) = T2
    
  END WHERE !isice
   
  RETURN

END SUBROUTINE ice_growth
!-------------------------------------------------------------------------------
!
!  
!>
!! ! upper_ocean_TS: Adjusts the temperature and salinity of the upper ocean grid
!!                 cell according to atmospheric heat and fresh-water fluxes,
!!                 surface melting, ice growth, etc. The upper ocean temperature
!!                 is also changed in subroutine new_ice_growth and at the
!!                beginning of subroutine ice_growth
!!
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
SUBROUTINE upper_ocean_TS(ppatch, p_os,ice, QatmAve, p_sfc_flx)
  TYPE(t_patch),       INTENT(IN)     :: ppatch 
  TYPE(t_hydro_ocean_state),INTENT(INOUT):: p_os
!  TYPE(t_atmos_for_ocean),INTENT(IN)  :: p_as
  TYPE(t_sea_ice),     INTENT (INOUT) :: ice
  TYPE(t_atmos_fluxes),INTENT (INOUT) :: QatmAve
  TYPE(t_sfc_flx),     INTENT (INOUT) :: p_sfc_flx

!!Local Variables
  REAL(wp), DIMENSION (nproma,i_no_ice_thick_class, ppatch%nblks_c) ::    &
   draft         ! position of ice-ocean interface below sea level        [m] 
  
  REAL(wp), DIMENSION (nproma, ppatch%nblks_c) ::   & 
   draftAve,    &! average draft of sea ice within a grid cell            [m]
   zUnderIceOld,&! water in upper ocean grid cell below ice (prev. time)  [m]
   heatOceI,    &! heat flux into ocean through formerly ice covered areas[W/m�]
   heatOceW,    &! heat flux into ocean through open water areas          [W/m�]
   delHice,     &! average change in ice thickness within a grid cell     [m]
   snowiceave,  &! average snow to ice conversion within a grid cell      [m]
   evap,        &! evaporated water                                       [m]
   preci,       &! solid precipitation                                    [m]
   precw         ! liquid precipitation                                   [m]

  ! Needs work with FB_BGC_OCE etc.
   REAL(wp)         :: swsum 
   REAL(wp),POINTER :: sao_top(:,:)
!  !-------------------------------------------------------------------------------

! #eoo# What is swsum?
  swsum = 0.0_wp
  sao_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

  ! Calculate change in water level 'zo' from liquid and solid precipitation and
  ! evaporation
  precw           (:,:)   = QatmAve% rprecw (:,:) * dtime
  preci           (:,:)   = QatmAve% rpreci (:,:) * dtime
  evap            (:,:)   = (QatmAve% latw(:,:)/ Lvap * dtime * &
                            sum(ice%conc(:,:,:), 2) +           &
                            sum(ice%evapwi(:,:,:) * ice% conc(:,:,:), 2)) /rhow

  ! TODO: This should probably be done via surface fluxes?
  !p_os%p_prog(nold(1))%h(:,:) = p_os%p_prog(nold(1))%h(:,:) +  precw + preci - evap

  ! Calculate average draft and thickness of water underneath ice in upper ocean
  ! grid box
  zUnderIceOld    (:,:)   = ice%zUnderIce
  draft           (:,:,:) = (rhos * ice%hs + rhoi * ice%hi) / rhow
  draftave        (:,:)   = sum(draft(:,:,:) * ice%conc(:,:,:),2)
  ice%zUnderIce   (:,:)   = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(:,:) - draftave(:,:) 
 
  ! Calculate average change in ice thickness and the snow-to-ice conversion 
  Delhice         (:,:)   = sum((ice% hi(:,:,:) - ice% hiold(:,:,:))*          &
                            ice%conc(:,:,:),2)
  snowiceave      (:,:)   = sum(ice%snow_to_ice(:,:,:) * ice% conc(:,:,:),2)
 

  ! Calculate heat input through formerly ice covered and through open water
  ! areas
  heatOceI        (:,:)   = sum(ice% heatOceI(:,:,:) * ice% conc(:,:,:),2)
  heatOceW        (:,:)   = (QatmAve%SWin(:,:) * (1.0_wp-albedoW) * (1.0_wp-swsum) +    &
                            QatmAve%LWnetw(:,:) + QatmAve%sensw(:,:)+         &
                            QatmAve%latw(:,:))  *  (1.0_wp-sum(ice%conc,2))

  ! Change temperature of upper ocean grid cell according to heat fluxes
!  p_os%p_prog(nold(1))%tracer(:,1,:,1) = p_os%p_prog(nold(1))%tracer(:,1,:,1)&
!                                       & + dtime*(heatOceI + heatOceW) /               &
!                                       & (clw*rhow * ice%zUnderIce)
! TODO: should we also divide with ice%zUnderIce / ( v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) ) ?
  p_sfc_flx%forc_tracer(:,:,1) = (heatOceI + heatOceW) / (clw*rhow)

! TODO:
!  ! Temperature change of upper ocean grid cell due  to melt-water inflow and
!  ! precipitation
!  p_os%p_prog(nold(1))%tracer(:,1,:,1) = (p_os%p_prog(nold(1))%tracer(:,1,:,1)&
!                          &*zUnderIceOld &
!                          &+ precw*p_as%tafo + preci*0.0_wp + &                             !!!!!!!!!Dirk: times 0.0 ????
!                          &  sum(ice%surfmeltT * ice%surfmelt * ice%conc,2)) / & 
!                          &  (zUnderIceOld + sum(ice%surfmelt*ice%conc,2) +    &
!                          &  precw + preci)
!
!  ! Change salinity of upper ocean grid box from ice growth/melt, snowice
!  ! formation and precipitation
!  p_os%p_prog(nold(1))%tracer(:,1,:,2) = p_os%p_prog(nold(1))%tracer(:,1,:,2)  &
!                                       & + (Delhice(:,:)*rhoi - snowiceave(:,:)*rhos)/rhow *  &
!                                       & MIN(Sice, sao_top(:,:)) / ice%zUnderIce(:,:)

  !heatabs         (:,:)   = swsum * QatmAve% SWin * (1 - ice%concsum)

END SUBROUTINE upper_ocean_TS
!-------------------------------------------------------------------------------
!
!  
!>
!! !! new_ice_growth: Calculates the grid-cell average thickness of new ice 
!                 forming in open-water areas
!!
!! @par Revision History
!! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
!! Dirk Notz, following MPI-OM. Code transfered to ICON.
!!
! TODO: This needs to be rewritten to take in to account cases where the ice concentration can vary
! between 0 and 1
SUBROUTINE new_ice_growth(ppatch,ice, p_os,p_sfc_flx)
  TYPE(t_patch),                INTENT(IN) :: ppatch 
  TYPE (t_sea_ice),         INTENT (INOUT) :: ice  
  TYPE(t_hydro_ocean_state),INTENT(INOUT)  :: p_os
!  TYPE (t_atmos_fluxes),    INTENT (IN   ) :: QatmAve
  TYPE(t_sfc_flx),     INTENT (INOUT) :: p_sfc_flx

  REAL(wp), DIMENSION (nproma,ppatch%nblks_c) :: sst

! Calculate possible super-cooling of the surface layer
  sst = p_os%p_prog(nold(1))%tracer(:,1,:,1) + &
   &      dtime*p_sfc_flx%forc_tracer(:,:,1)/ice%zUnderIce

  ice % newice = 0.0_wp
  WHERE (sst <= Tf)
    ice%newice(:,:) = - (sst - Tf) * ice%zUnderIce * clw*rhow / (Lfreez*rhoi)
    ! Add energy for new-ice formation due to supercooled ocean to  ocean temperature
    p_sfc_flx%forc_tracer(:,:,1) = &
      &     ice%zUnderIce * ( Tf - p_os%p_prog(nold(1))%tracer(:,1,:,1) ) / dtime
  END WHERE

  WHERE(ice%newice>0.0_wp)
    WHERE(.NOT.ice%isice(:,1,:))
      ice%Tsurf(:,1,:) = Tf
      ice%T2   (:,1,:) = Tf
      ice%T1   (:,1,:) = Tf
    ENDWHERE
    ice % isice(:,1,:) = .TRUE.
    ice % hi   (:,1,:) = ice%newice* (1.0_wp-sum(ice%conc,2))&
                       &+ice%hi(:,1,:)*sum(ice%conc,2)
     !ice % hs   (:,:,1) = 0
!     ice % Tsurf(:,1,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)
!!!!!!!!!!!DIRK: Where is rhs coming from ???????????????

!     ice % T1   (:,:,1) = T1(:,:,1)
!     ice % T2   (:,:,1) = T2(:,:,1)
    ice % conc (:,1,:) = 1.0_wp
  ENDWHERE
  ice% concSum(:,:)  = SUM(ice% conc(:,:,:),2)
END SUBROUTINE new_ice_growth


  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of the ice 
  !! temperature and ice growth rates for OMIP forcing
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !
  SUBROUTINE calc_atm_fluxes_from_bulk(ppatch, p_as, p_os, p_ice, Qatm)
  TYPE(t_patch),            INTENT(IN)    :: ppatch
  TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
  TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
  TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
  TYPE(t_atmos_fluxes),     INTENT(INOUT) :: Qatm


  !Local variables
  REAL(wp), DIMENSION (nproma,ppatch%nblks_c) ::           &
    z_Tsurf,      &  ! Surface temperature                             [C]
    z_tafoK,      &  ! Air temperature at 2 m in Kelvin                [K]
    z_fu10lim,    &  ! wind speed at 10 m height in range 2.5...32     [m/s]
    z_esta,       &  ! water vapor pressure at 2 m height              [Pa]
    z_esti,       &  ! water vapor pressure at ice surface             [Pa]
    z_estw,       &  ! water vapor pressure at water surface           [Pa]
    z_sphumida ,  &  ! Specific humididty at 2 m height 
    z_sphumidi ,  &  ! Specific humididty at ice surface
    z_sphumidw ,  &  ! Specific humididty at water surface
    z_ftdewC,     &  ! Dew point temperature in Celsius                [C]
    z_rhoair ,    &  ! air density                                     [kg/m³]
    z_dragl1,     &  ! part of z_dragl                                   
    z_dragl ,     &  ! Drag coefficient for latent   heat flux
    z_drags ,     &  ! Drag coefficient for sensible heat flux (=0.95 z_dragl)
    z_xlat ,      &  ! latitude limited to 60S...60N
    z_fakts ,     &  ! Effect of cloudiness on LW radiation
    z_humi           ! Effect of air humidity on LW radiation
  
  INTEGER :: jc, jb, i
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:calc_atm_fluxes_from_bulk'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  rl_start_c = 1
  rl_end_c   = min_rlcell
  i_startblk_c = ppatch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)

  !TODO: Can't we get rid of the do-loops?
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( ppatch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
 
      z_Tsurf (jc,jb) = p_os%p_prog(nold(1))%tracer(jc,1,jb,1)        ! set surface temp = mixed layer temp
      z_tafoK (jc,jb) = p_as%tafo  (jc,jb) + tmelt                    ! Change units of z_tafoK  to Kelvin
      z_ftdewC(jc,jb) = p_as%ftdew (jc,jb) - tmelt                    ! Change units of z_ftdewC to C

      z_xlat   (jc,jb) = MIN(ABS(ppatch%cells%center(jc,jb)%lat*rad2deg),60.0_wp) 

      !-----------------------------------------------------------------------
      ! Compute water vapor pressure and specific humididty in 2m height (z_esta) 
      ! and at water surface (z_estw) according to "Buck Research Manual (1996)
      ! (see manuals for instruments at http://www.buck-research.com/); 
      ! updated from Buck, A. L., New equations for computing vapor pressure and 
      ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981" 
      !-----------------------------------------------------------------------

      ! Buck 1981
      z_esta(jc,jb)  = 611.21_wp * EXP( (18.729_wp-z_ftdewC(jc,jb)/227.3_wp)*z_ftdewC(jc,jb)&
                                    &/ (z_ftdewC(jc,jb)+257.87_wp) )
      ! Buck 1996
      !z_esta(:,:) = 611.21 * EXP( (18.678-z_ftdewC/234.5)*z_ftdewC/ (z_ftdewC+257.14) )
      ! Buck 1981
      z_estw(jc,jb)  = 611.21_wp*EXP( (18.729_wp-z_Tsurf(jc,jb)/227.3_wp)&
                     & * z_Tsurf(jc,jb) /  (z_Tsurf(jc,jb) +257.87_wp) )
      ! Buck 1996
      !z_estw(:,:) = 611.21 * EXP( (18.678-z_Tsurf /234.5)*z_Tsurf/  (z_Tsurf +257.14) )
      z_estw(jc,jb)  = 0.9815_wp * z_estw(jc,jb)
      !or more accurate: (1-5.27e-4 * mixed layer salinity) * z_estw (Kraus and
      ! Businger, 1994)

      z_sphumida(jc,jb)  = 0.62197_wp * z_esta(jc,jb)/(p_as%pao(jc,jb)-0.37803_wp*z_esta(jc,jb))
      z_sphumidw (jc,jb) = 0.62197_wp * z_estw(jc,jb)/(p_as%pao(jc,jb)-0.37803_wp*z_estw(jc,jb))

      !-----------------------------------------------------------------------
      !  Compute longwave radiation according to 
      !         Koch 1988: A coupled Sea Ice - Atmospheric Boundary Layer Model,
      !                    Beitr.Phys.Atmosph., 61(4), 344-354.
      !  or (ifdef QLOBERL)
      !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
      !         long-wave radiation of the Earth with consideration of the effect
      !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
      !         cited by: Budyko, Climate and Life, 1974.
      !         Note that for z_humi, z_esta is given in [mmHg] in the original
      !         publication. Therefore, 0.05*sqrt(z_esta/100) is used rather than
      !         0.058*sqrt(z_esta)
      !-----------------------------------------------------------------------

      z_humi   (jc,jb) = 0.601_wp+ 5.95_wp*1.0e-7_wp*z_esta(jc,jb)*EXP(1500.0_wp/z_tafoK(jc,jb))
      z_fakts  (jc,jb) =  1.0_wp + 0.3_wp*p_as%fclou(jc,jb)**2
      Qatm%LWin(jc,jb) = z_fakts(jc,jb) * z_humi(jc,jb) * emiss*StefBol * z_tafoK(jc,jb)**4

      Qatm%LWoutw(jc,jb) = emiss*StefBol * (z_Tsurf(jc,jb)+273.15_wp)**4
      Qatm%LWnetw(jc,jb) = Qatm%LWin(jc,jb) - Qatm%LWoutw(jc,jb)

      Qatm%SWin(jc,jb) = p_as%fswr(jc,jb)

      !-----------------------------------------------------------------------
      !  Calculate bulk equations according to 
      !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002: 
      !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
      !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
      !-----------------------------------------------------------------------    
      z_rhoair  (jc,jb) = p_as%pao(jc,jb)/(rgas*z_tafoK(jc,jb)*(1.0_wp+0.61_wp*z_sphumida(jc,jb)) )
      z_fu10lim (jc,jb) = MAX (2.5_wp, MIN(32.5_wp,p_as%fu10(jc,jb)) )
      z_dragl1  (jc,jb) = 1e-3_wp*(-0.0154_wp + 0.5698_wp/z_fu10lim(jc,jb)                 &
                        & -0.6743_wp/(z_fu10lim(jc,jb) * z_fu10lim(jc,jb)))
      z_dragl   (jc,jb) = 1.0e-3_wp*(0.8195_wp+0.0506_wp*z_fu10lim(jc,jb)  &
                        &-0.0009_wp*z_fu10lim(jc,jb)*z_fu10lim(jc,jb)) + z_dragl1(jc,jb)   &
                        &* (z_Tsurf(jc,jb)-p_as%tafo(jc,jb))
      z_dragl   (jc,jb) = MIN (z_dragl(jc,jb), 3.0E-3_wp)
      z_drags   (jc,jb) = 0.96_wp * z_dragl(jc,jb)
      Qatm%sensw(jc,jb) = z_drags(jc,jb)*z_rhoair(jc,jb)*cpa*p_as%fu10(jc,jb)             &
                        & * (p_as%tafo(jc,jb) -z_Tsurf(jc,jb))  *fr_fac
      Qatm%latw (jc,jb) = z_dragl(jc,jb)*z_rhoair(jc,jb)*Lvap*p_as%fu10(jc,jb)          &
                        & * (z_sphumida(jc,jb)-z_sphumidw(jc,jb))*fr_fac

      DO i = 1, p_ice%kice
        IF (p_ice% isice(jc,i,jb))THEN
          z_Tsurf(jc,jb) = p_ice%Tsurf(jc,i,jb)
          ! z_esti is calculated according to Buck Research Manuals, 1996 (see z_esta)
          z_esti     (jc,jb) = 611.15_wp*EXP( (23.036_wp-z_Tsurf(jc,jb)/333.7_wp) &
                             & *z_Tsurf(jc,jb)/(z_Tsurf(jc,jb) + 279.82_wp) )
          z_sphumidi (jc,jb) = 0.62197_wp*z_esti(jc,jb)/(p_as%pao(jc,jb)-0.37803_wp*z_esti(jc,jb))
          z_dragl    (jc,jb) = MAX (0.5e-3_wp, 1.0e-3_wp * (0.8195_wp+0.0506_wp*z_fu10lim(jc,jb) &
                             & -0.0009_wp*z_fu10lim(jc,jb) * z_fu10lim(jc,jb)) + z_dragl1(jc,jb) &
                             & * (z_Tsurf(jc,jb)-p_as%tafo(jc,jb)) )
          z_drags    (jc,jb) = 0.95_wp * z_dragl(jc,jb)

          Qatm%LWout (jc,i,jb) = emiss*StefBol * (z_Tsurf(jc,jb)+273.15_wp)**4
          Qatm%LWnet (jc,i,jb) = Qatm%LWin(jc,jb) - Qatm%LWout(jc,i,jb)
          Qatm%dLWdT (jc,i,jb) = - 4.0_wp * emiss*StefBol * (z_Tsurf(jc,jb) + 273.15_wp)**3
          Qatm%sens  (jc,i,jb) = z_drags(jc,jb) * z_rhoair(jc,jb)*cpa*p_as%fu10(jc,jb)&
                               & * (p_as%tafo(jc,jb) -z_Tsurf(jc,jb))   *fr_fac
          Qatm%lat   (jc,i,jb) = z_dragl(jc,jb) * z_rhoair(jc,jb)* Lfreez *p_as%fu10(jc,jb)&
                               & * (z_sphumida(jc,jb)-z_sphumidi(jc,jb))*fr_fac
    
          Qatm%dsensdT(jc,i,jb)= 0.95_wp*z_dragl1(jc,jb)*z_drags(jc,jb)*z_rhoair(jc,jb)&
                               & *cpa * p_as%fu10(jc,jb)       &
                               & * (p_as%tafo(jc,jb) - z_Tsurf(jc,jb)) *  fr_fac &
                               & -z_drags(jc,jb)*z_rhoair(jc,jb) *cpa*p_as%fu10(jc,jb)
          Qatm%dlatdT(jc,i,jb) = z_dragl1(jc,jb) * z_rhoair(jc,jb)*Lfreez * p_as%fu10(jc,jb)&
                               & *(z_sphumida(jc,jb)-z_sphumidi(jc,jb))*fr_fac
        ENDIF
      ENDDO

      !Dirk: why zero ?
      Qatm%rpreci(jc,jb) = 0.0_wp
      Qatm%rprecw(jc,jb) = 0.0_wp

    END DO
  END DO
 
  END SUBROUTINE calc_atm_fluxes_from_bulk
 
  !-------------------------------------------------------------------------


END MODULE mo_sea_ice
