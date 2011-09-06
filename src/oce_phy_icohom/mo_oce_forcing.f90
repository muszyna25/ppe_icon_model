!>
!! Provide an implementation of the ocean forcing.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modification by Stephan Lorenz, MPI-M (2010-06):
!!   - renaming and adjustment to ocean domain and patch_oce
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!    adapted to structures discussed in 2010-01.
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
MODULE mo_oce_forcing
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
USE mo_kind,                ONLY: wp
USE mo_parallel_config,     ONLY: nproma
USE mo_ocean_nml,           ONLY: no_tracer
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch
USE mo_oce_state,           ONLY: t_hydro_ocean_state
USE mo_exception,           ONLY: finish, message !, message_text
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary!,    &
! &                               land, sea, boundary,                  &
! &                               min_rlcell, min_rledge, min_rlvert
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec, cvec2gvec
!USE mo_param_ice !,           ONLY: kice
USE mo_sea_ice,             ONLY: t_sea_ice
IMPLICIT NONE

PRIVATE


CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface

! public subroutines
PUBLIC  :: construct_sfcflx
PUBLIC  :: construct_atmos_for_ocean
PUBLIC  :: construct_atmos_fluxes
PUBLIC  :: destruct_sfcflx
PUBLIC  :: destruct_atmos_for_ocean
PUBLIC  :: destruct_atmos_fluxes

! public types
PUBLIC  :: t_sfc_flx
PUBLIC  :: t_atmos_fluxes
PUBLIC  :: t_atmos_for_ocean


!------  Definition of surface flux type---------------------

TYPE t_sfc_flx

! The forcing is specified as fluxes at the air-sea interface defined on cell-centers
! dimension: (nproma, nblks_c)
  REAL(wp), POINTER ::          &
    &  forc_wind_u(:,:),            & !forcing of zonal component of velocity equation,
    &  forc_wind_v(:,:),            & !forcing of meridional component of velocity equation,
    &  forc_tracer(:,:,:),          & !tracer flux. Last index refers to tracer id (1=heat, 2=fresh-water flux)
    &  forc_tracer_relax(:,:,:)       !tracer relaxation: contains data to which is relaxated. e.g. climatology. 
                                      !Last index refers to tracer id (1=temperature, 2=salinity)
  TYPE(t_cartesian_coordinates),    & !wind forcing with cartesian vector, located at cell centers
  & ALLOCATABLE :: forc_wind_cc(:,:) 

END TYPE t_sfc_flx

! publice variables
TYPE(t_sfc_flx),PUBLIC,TARGET :: v_sfc_flx

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

  REAL(wp), ALLOCATABLE ::          &
    &  forc_wind_u(:,:),            & !forcing of zonal component of velocity equation,
    &  forc_wind_v(:,:),            & !forcing of meridional component of velocity equation,
    &  forc_tracer(:,:,:),          & !tracer flux. Last index refers to tracer id (1=heat, 2=fresh-water flux)
    &  forc_tracer_relax(:,:,:)       !tracer relaxation: contains data to which is relaxated. 
                                      !Last index refers to tracer id (1=temperature, 2=salinity)
  TYPE(t_cartesian_coordinates),    & !wind forcing with cartesian vector, located at cell centers
  & ALLOCATABLE :: forc_wind_cc(:,:) 

END TYPE t_atmos_fluxes


CONTAINS


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
  INTEGER :: nblks_c, ist, jc,jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:construct_ho_sfcflx'

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
    p_sfc_flx%forc_tracer       = 0.0_wp
    p_sfc_flx%forc_tracer_relax = 0.0_wp
  ENDIF
  CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_sfcflx
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
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:construct_atmos_for_ocean'

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
  INTEGER :: nblks_c, ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:destruct_atmos_for_ocean'

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
  INTEGER :: nblks_c, ist, jc,jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:construct_ho_sfcflx'

  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  nblks_c = ppatch%nblks_c
 
   ALLOCATE(p_atm_f%sens(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sens failed')
   END IF

   ALLOCATE(p_atm_f%lat(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for lat failed')
   END IF

   ALLOCATE(p_atm_f%LWout(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWout failed')
   END IF

   ALLOCATE(p_atm_f%LWnet(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for LWnet failed')
   END IF

   ALLOCATE(p_atm_f%bot(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sens failed')
   END IF

   ALLOCATE(p_atm_f%dsensdT(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for dsensdT failed')
   END IF

   ALLOCATE(p_atm_f%dlatdT(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
   IF (ist/=SUCCESS) THEN
     CALL finish(TRIM(routine),'allocation for sens failed')
   END IF

   ALLOCATE(p_atm_f%dLWdT(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
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
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:construct_ho_sfcflx'
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
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:destruct_sfcflx'
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
!   DEALLOCATE(p_sfc_flx%forc_freshw, STAT=ist)
!   IF (ist/=SUCCESS) THEN
!     CALL finish(TRIM(routine),'deallocation for forcing freshwater failed')
!   END IF
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

END MODULE mo_oce_forcing
