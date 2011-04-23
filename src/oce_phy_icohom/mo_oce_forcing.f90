!>
!! Provide an implementation of the ocean forcing.
!!
!! Provide an implementation of the parameters used for top forcing
!! of the hydrostatic ocean model.
!! Although forcing is part of the time-variant diagnostics, its different
!! types are stored in this module.
!! Routine init_ho_forcing controls which structure of surface forcing is
!! constructed, initialised, updated, and destructed.
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
USE mo_run_nml,             ONLY: nproma
USE mo_ocean_nml,           ONLY: iforc_oce, analyt_stat, no_tracer!, wstress_coeff
USE mo_model_domain,        ONLY: t_patch
USE mo_exception,           ONLY: finish, message !, message_text
USE mo_math_constants,      ONLY: pi
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary!,    &
! &                               land, sea, boundary,                     &
! &                               min_rlcell, min_rledge, min_rlvert,      &
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec
IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface

! public subroutines
PUBLIC :: construct_ho_sfcflx, destruct_ho_sfcflx, update_ho_sfcflx
!PUBLIC :: construct_ho_coreforc, destruct_ho_coreforc, update_ho_coreforc
!PUBLIC :: construct_ho_fullforc, destruct_ho_fullforc, update_ho_fullforc

! public types
PUBLIC :: t_ho_sfc_flx, t_ho_core_forc, t_ho_full_forc

! public variables
!PUBLIC :: p_sfc_flx!, f_coreforc, f_fullforc

! Definition of forcing types

TYPE t_ho_sfc_flx

! The forcing is specified as fluxes at the air-sea interface defined on cell-centers
! dimension: (nproma, nblks_c)

  REAL(wp), ALLOCATABLE :: &
    &  forc_wind_u(:,:),   &    !forcing of zonal component of velocity equation,
    &  forc_wind_v(:,:),   &    !forcing of meridional component of velocity equation,
    &  forc_freshw(:,:),   &    !(mass-)forcing of height equation (fresh-water flux)
    &  forc_tracer(:,:,:), &    !forcing of tracers. Last index refers to tracer id (1=temperature, 2=salinity)
    &  forc_tracer_relax(:,:,:) !tracer relaxation. Last index refers to tracer id (1=temperature, 2=salinity)
!     &  forc_heat(:,:),     & ! forcing of temperature equation (heat-flux)
!     &  forc_temp(:,:),     & ! forcing of temperature equation (restoring)
!     &  forc_sal(:,:)         ! forcing of salinity equation (restoring)

    TYPE(t_cartesian_coordinates), ALLOCATABLE :: forc_wind_cc(:,:) !wind forcing with cartesian vector, located at cell centers

END TYPE t_ho_sfc_flx

TYPE t_ho_full_forc

! Components related to complete set of forcing functions (cleanup if overlap with other arrays)
! Nomenclature adheres to mo_ncar_fluxes from MPIOM.
! Arrays used in bulk formula

  REAL(wp), ALLOCATABLE ::  &
    &  wind_stress_u(:,:),  &  ! windstress coefficients normal to edge,     dim: (nproma, nblks_e)
    &  wind_stress_v(:,:),  &  ! windstress coefficients tangential to edge, dim: (nproma, nblks_e)
    &  du           (:,:),  &  !                                             dim: (nproma, nblks_e)
    &  dv           (:,:),  &  !                                             dim: (nproma, nblks_e)
    &  uvdel        (:,:),  &  !                                             dim: (nproma, nblks_e)
    &  uvdelc       (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  qo           (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  z            (:,:),  &  !                                             dim: (nproma, nblks_e)
    &  qpre         (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  qla          (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  qse          (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  qnsw         (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  qnlw         (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  qnet         (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  fw           (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  evap         (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  ce           (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  cd           (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  cde          (:,:),  &  !                                             dim: (nproma, nblks_e)
    &  ch           (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  ustar        (:,:),  &  !                                             dim: (nproma, nblks_c)
    &  bstar        (:,:)      !                                             dim: (nproma, nblks_c)
END TYPE t_ho_full_forc

TYPE t_ho_core_forc
! external CORE forcing functions themselves

  REAL(wp), ALLOCATABLE ::  &
    &  corv10n  (:,:),  &  ! 10m normal wind speed [m s-1]        dim: (nproma, nblks_e)
    &  corv10t  (:,:),  &  ! 10m tangential wind speed [m s-1]    dim: (nproma, nblks_e)
    &  cort10   (:,:),  &  ! 10m air temperature [C]              dim: (nproma, nblks_c)
    &  corq10   (:,:),  &  ! 10m humidity [kg/kg]                 dim: (nproma, nblks_c)
    &  corqdlw  (:,:),  &  ! downward longwave radiation [W m-2]  dim: (nproma, nblks_c)
    &  corqdsw  (:,:),  &  ! downward shortwave radiation [W m-2] dim: (nproma, nblks_c)
    &  corprec  (:,:),  &  ! precipitation [kg m-2]               dim: (nproma, nblks_c)
    &  corrunoff(:,:)      ! runoff [kg m-2]                      dim: (nproma, nblks_c)
END TYPE t_ho_core_forc

! SAVE p_sfc_flx, f_coreforc, f_fullforc
! TYPE(t_ho_sfc_flx)    :: p_sfc_flx
! TYPE(t_ho_core_forc)  :: f_coreforc
! TYPE(t_ho_full_forc)  :: f_fullforc


CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Construct surface flux forcing for hydrostatic ocean
  !!
  !! Construct surface flux forcing for hydrostatic ocean consisting of
  !! fluxes at the air-sea interface defined on cell-centers
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE construct_ho_sfcflx(ppatch, p_sfc_flx)
  !
  TYPE(t_patch), INTENT(in) :: ppatch
  TYPE(t_ho_sfc_flx)        :: p_sfc_flx

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
  ALLOCATE(p_sfc_flx%forc_freshw(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing freshwater failed')
  END IF


  ALLOCATE(p_sfc_flx%forc_tracer(nproma,nblks_c, no_tracer), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for tracer forcing failed')
  END IF

  ALLOCATE(p_sfc_flx%forc_tracer_relax(nproma,nblks_c, no_tracer), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for tracer relaxation forcing failed')
  END IF

  ALLOCATE(p_sfc_flx%forc_wind_cc(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for forcing wind_cc  failed')
  END IF

  p_sfc_flx%forc_wind_u   = 0.0_wp
  p_sfc_flx%forc_wind_v   = 0.0_wp
  p_sfc_flx%forc_freshw   = 0.0_wp

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
    ENDDO
  END DO


  p_sfc_flx%forc_tracer       = 0.0_wp
  p_sfc_flx%forc_tracer_relax = 0.0_wp

  CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_ho_sfcflx
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Destruct surface flux forcing for hydrostatic ocean
  !!
  !! Destruct surface flux forcing for hydrostatic ocean ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE destruct_ho_sfcflx(p_sfc_flx)
  TYPE(t_ho_sfc_flx) :: p_sfc_flx
  !

  ! Local variables

  INTEGER :: ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:destruct_ho_sfcflx'

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
  DEALLOCATE(p_sfc_flx%forc_freshw, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for forcing freshwater failed')
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

  END SUBROUTINE destruct_ho_sfcflx
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !! Update surface flux forcing for hydrostatic ocean ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE update_ho_sfcflx(p_patch, p_sfc_flx)

  TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
  TYPE(t_ho_sfc_flx)                :: p_sfc_flx
  !
  ! local variables

!  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: jc, jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start, rl_end


  !REAL(wp), PARAMETER :: alpha = 0.0_wp 
  !REAL(wp), PARAMETER :: u_0 = 1.0_wp
  !REAL(wp) :: u_0, rforc, hforc
  !REAL(wp) :: tau_0
  REAL(wp) :: z_lat, z_lon  !, z_L_phi
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:update_ho_sfcflx'
  !-------------------------------------------------------------------------
    rl_start = 1
    rl_end = min_rlcell
    i_startblk_c = p_patch%cells%start_blk(rl_start,1)
    i_endblk_c   = p_patch%cells%end_blk(rl_end,1)


  !IF(FORCING_TYPE==ANALYTICAL_STATIONARY)THEN
  SELECT CASE (iforc_oce)

  CASE (0)
    CALL message(TRIM(routine), 'No wind forcing applied' )

  CASE (analyt_stat)

    CALL message(TRIM(routine), 'Apply stationary wind forcing' )
    !This implements formula (75), (76) in the Williamson JCP paper from 1992

    ! #slo# 2010-12-08: Stommel gyre test (see MITGCM and Comblen et al. 2008)
    ! Tau = tau_0 * sin(Pi*phi/L_phi)
    ! A 60x60 deg ocean basin is used here - 30W-30E, 30S-30N (=-pi/3,+pi/3), L_phi=60 deg = pi/6
    !Latitudes in ICON vary from -pi/2 to pi/2, if the latitude equals 0, then the
    !forcing should be maximal, i.e. the cos has to evaluate to 1 and its argument should equal 0.
    !tau_0   = wstress_coeff       !  see MITGCM
    !z_L_phi = pi/6.0_wp

    ! #slo# 2011-03-10: Munk gyre flow (see Danilov, Oc.Dyn. 2010)
    !   - now a basin in hexagonal form along the original icosahedral triangles is used
    !     10 triangles in r2b0 are water: 50W-50E, 30S-30N
    !   - in order to avoid negative forcing the scaling factor in cos(z_lat*3) is reduced
    !   - tau_0*rho_0 is set via wstress_coeff=0.2
 
    DO jb = i_startblk_c, i_endblk_c

      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, &
      &                                i_startidx_c, i_endidx_c,rl_start,rl_end)
      DO jc = i_startidx_c, i_endidx_c

        z_lat = p_patch%cells%center(jc,jb)%lat
        z_lon = p_patch%cells%center(jc,jb)%lon

        IF(p_patch%patch_oce%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
!         p_sfc_flx%forc_wind_u(jc,jb) = cos(z_lat*3.0_wp)
          p_sfc_flx%forc_wind_u(jc,jb) = cos(z_lat*2.8_wp)
        ELSE
          p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
        ENDIF !write(*,*)'forc',jc,jb, z_lat, z_lat*180.0_wp/pi, p_sfc_flx%forc_wind_u(jc,jb) 
        p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp

        !Init cartesian wind
        CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
                       & p_sfc_flx%forc_wind_v(jc,jb),&
                       & p_patch%cells%center(jc,jb)%lon,&
                       & p_patch%cells%center(jc,jb)%lat,&
                       & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                       & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                       & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
                       !write(*,*)'sfc forcing', jc,jb,p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
      END DO
   END DO
   write(*,*)'max/min-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)

! CASE(NCEP)
!   CALL message(TRIM(routine), 'STOP: Option not implemented yet' )
! CASE(CORE)
!   CALL message(TRIM(routine), 'STOP: Option not implemented yet' )

  CASE DEFAULT
    CALL message(TRIM(routine), 'STOP: Forcing option not implemented yet' )
!   CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')
  END SELECT

  END SUBROUTINE update_ho_sfcflx

  !-------------------------------------------------------------------------
END MODULE mo_oce_forcing
