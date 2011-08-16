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
!!  Modified by Stephan Lorenz,     MPI-M (2011-07)
!!   - for parallel ocean: 3-dim ocean grid in v_base
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
MODULE mo_oce_bulk
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
USE mo_io_units,            ONLY: filename_max
USE mo_mpi,                 ONLY: p_pe, p_io, p_bcast
USE mo_ext_data,            ONLY: ext_data
USE mo_grid_config,         ONLY: nroot
USE mo_ocean_nml,           ONLY: iforc_oce, itestcase_oce, no_tracer,                          &
  &                               basin_center_lat, basin_center_lon, basin_width_deg,&
  &                               basin_height_deg, relaxation_param, wstress_coeff,  &
  &                               i_bc_veloc_top, temperature_relaxation,             &
  &                               NO_FORCING, ANALYT_FORC, FORCING_FROM_FILE_FLUX,    &
  &                               FORCING_FROM_FILE_FIELD, FORCING_FROM_COUPLED_FLUX, &
  &                               FORCING_FROM_COUPLED_FIELD
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch
USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base
USE mo_exception,           ONLY: finish, message, message_text
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_physical_constants,  ONLY: rho_ref, sfc_press_bar, lsub, lvap, cpa, emiss, &
  &                               fr_fac, stefbol, rgas, L, tmelt
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec, cvec2gvec
USE mo_sea_ice,             ONLY: t_sea_ice
USE mo_oce_forcing,         ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
USE mo_oce_thermodyn,       ONLY:convert_insitu2pot_temp_func
IMPLICIT NONE

! required for reading netcdf files
INCLUDE 'netcdf.inc'

PRIVATE


CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface

! public subroutines
PUBLIC  :: update_sfcflx

! private implementation
PRIVATE :: calc_atm_fluxes_from_bulk
PRIVATE :: update_sfcflx_analytical

CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  SUBROUTINE update_sfcflx(p_patch, p_os, p_as, p_ice, Qatm, p_sfc_flx, jstep)

  TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
  TYPE(t_hydro_ocean_state)           :: p_os
  TYPE(t_atmos_for_ocean), INTENT(IN) :: p_as
  TYPE (t_atmos_fluxes)               :: Qatm
  TYPE (t_sea_ice)                    :: p_ice
  TYPE(t_sfc_flx)                     :: p_sfc_flx
  INTEGER, INTENT(IN)                 :: jstep
  !
  ! local variables
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_sfcflx'
  INTEGER :: jmon, jmonst, njday, jdays, jdmon, jmon1, jmon2
  INTEGER :: jc, jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c
  REAL(wp):: z_relax, rday1, rday2
  !REAL(wp):: z_omip_data(nproma,12,p_patch%nblks_c,3)  ! 3 arrays to be read
  !-------------------------------------------------------------------------
  rl_start_c = 1
  rl_end_c   = min_rlcell
  i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)


  SELECT CASE (iforc_oce)

  CASE (NO_FORCING)

    CALL message(TRIM(routine), 'No  forcing applied' )

  CASE (ANALYT_FORC)

    CALL update_sfcflx_analytical(p_patch, p_os, p_sfc_flx)

  CASE (FORCING_FROM_FILE_FLUX)

    ! To Do:
    ! 2 possibilities: read forcing file either completely or in chunks
    !
    ! a) Check if file should be read (new chunk if timecriterion is met, or completely)
    ! 
    ! Quick and dirty:
    !  - model currently starts at first of March
    !  - later: use datetime for correct evaluation of month
    !jmonst = ini_datetime%month
    jmonst = 3
    !  calculate day and month
    njday = int(86400._wp/dtime)  ! no of timesteps per day
    jdays = (jstep-1)/njday + 1   ! no of days from start
    jmon  = (jdays-1)/30 + jmonst ! month 
    if (jmon > 12) jmon=1
    jdmon = mod(jdays+1,30)-1     ! no of days in month
    !jdstp = mod(jstep-1,njday)    ! 

    write(*,*) ' jstep, njday, jdays, jmon, jdmon : ',jstep, njday, jdays, jmon, jdmon

    !
    ! interpolate omip-data daily:
    !
    !IF ( MOD(jstep-1,njday)==0) THEN

      jmon1=jmon-1
      jmon2=jmon
      rday1=REAL(15-jdmon,wp)/30.0_wp
      rday2=REAL(15+jdmon,wp)/30.0_wp
      if (jdmon > 15)  then
        jmon1=jmon
        jmon2=jmon+1
        rday1=REAL(15+jdmon,wp)/30.0_wp
        rday2=REAL(15-jdmon,wp)/30.0_wp
      end if

      if (jmon1 ==  0) jmon1=12
      if (jmon1 == 13) jmon1=1
      if (jmon2 ==  0) jmon2=12
      if (jmon2 == 13) jmon2=1

      write(*,*) ' jstep, jmon1, jmon2, rday1, rday2: ',jstep, jmon1, jmon2, rday1, rday2

      !
      ! OMIP data read in mo_ext_data into variable ext_data
      ! ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean
      !

      p_sfc_flx%forc_wind_u(:,:) = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,1) + &
        &                          rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,1)
      p_sfc_flx%forc_wind_v(:,:) = rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,2) + &
        &                          rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,2)

      IF (temperature_relaxation == 2)  THEN
         p_sfc_flx%forc_tracer_relax(:,:,1) = &
           &  rday1*ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,3) + &
           &  rday2*ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,3)

        ! subtract tmelt and set to zero
        DO jb = i_startblk_c, i_endblk_c
          CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
            &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
          DO jc = i_startidx_c, i_endidx_c
            IF (v_base%lsm_oce_c(jc,1,jb) <= sea_boundary) THEN
              p_sfc_flx%forc_tracer_relax(jc,jb,1) &
                & = p_sfc_flx%forc_tracer_relax(jc,jb,1) - tmelt
            ELSE
              p_sfc_flx%forc_tracer_relax(jc,jb,1) = 0.0_wp
            END IF
          END DO
        END DO

      END IF
       

      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
          &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
                           & p_sfc_flx%forc_wind_v(jc,jb),&
                           & p_patch%cells%center(jc,jb)%lon,&
                           & p_patch%cells%center(jc,jb)%lat,&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb)         = 0.0_wp
            p_sfc_flx%forc_wind_v(jc,jb)         = 0.0_wp
            p_sfc_flx%forc_wind_cc(jc,jb)%x      = 0.0_wp
          ENDIF
        END DO
      END DO

    !END IF  ! interpolate daily

    WRITE(*,*)'MAX/MIN omip_data-u mon=',jmon1,&
      &        MAXVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,1)),&
      &        MINVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,1)) 
    WRITE(*,*)'MAX/MIN omip_data-v mon=',jmon1,&
      &        MAXVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,2)),&
      &        MINVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,2)) 
    WRITE(*,*)'MAX/MIN omip_data-t mon=',jmon1,&
      &        MAXVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,3)),&
      &        MINVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon1,:,3)) 
    WRITE(*,*)'MAX/MIN omip_data-u mon=',jmon2,&
      &        MAXVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,1)),&
      &        MINVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,1)) 
    WRITE(*,*)'MAX/MIN omip_data-v mon=',jmon2,&
      &        MAXVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,2)),&
      &        MINVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,2)) 
    WRITE(*,*)'MAX/MIN omip_data-t mon=',jmon2,&
      &        MAXVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,3)),&
      &        MINVAL(ext_data(1)%oce%omip_forc_mon_c(:,jmon2,:,3)) 

    WRITE(*,*)'MAX/MIN forc_wind_u:',       1, &
      &        MAXVAL(p_sfc_flx%forc_wind_u(:,:)),&
      &        MINVAL(p_sfc_flx%forc_wind_u(:,:)) 
    WRITE(*,*)'MAX/MIN forc_wind_v:',       1, &
      &        MAXVAL(p_sfc_flx%forc_wind_v(:,:)),&
      &        MINVAL(p_sfc_flx%forc_wind_v(:,:))
    WRITE(*,*)'MAX/MIN forc_tracer_relax:', 1, &
      &        MAXVAL(p_sfc_flx%forc_tracer_relax(:,:,1)),&
      &        MINVAL(p_sfc_flx%forc_tracer_relax(:,:,1))

  CASE (FORCING_FROM_FILE_FIELD)
    ! 1) Read field data from file
    ! 2) CALL calc_atm_fluxes_from_bulk (p_patch, p_as, p_os, p_ice, Qatm)
    ! 3)CALL update_sfcflx_from_atm_flx(p_patch, p_as, p_os, p_ice, Qatm, p_sfc_flx)

  CASE (FORCING_FROM_COUPLED_FLUX,FORCING_FROM_COUPLED_FIELD)
    !Depending on coupling type, apply one of two following ways:
    !1) bulk formula to atmospheric state and proceed as above, the only distinction
    !   to OMIP is that atmospheric info is coming from model rather than file
    !2) use atmospheric fluxes directly, i.e. avoid call to "calc_atm_fluxes_from_bulk"
    !    and do a direct assignment of atmospheric state to surface fluxes.

  CASE DEFAULT

    CALL message(TRIM(routine), 'STOP: Forcing option not implemented yet' )
    CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')

  END SELECT

  ! Temperature relaxation: This is a provisional solution
  IF(temperature_relaxation==1)THEN
    z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
    DO jb = i_startblk_c, i_endblk_c    
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
       &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c

        IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
          p_sfc_flx%forc_tracer(jc,jb, 1)=                       &!Compare form of relaxation with MPI-OM
          &          z_relax&!*v_base%del_zlev_m(1)     &
          &          *(p_sfc_flx%forc_tracer_relax(jc,jb,1)-     &
          &            p_os%p_prog(nold(1))%tracer(jc,1,jb,1))

        ENDIF
      END DO
    END DO

    write(*,*)'max/min-tracer-diff',&
    & maxval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1)),&
    & minval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1))
    
    write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
    & minval(p_sfc_flx%forc_tracer_relax)
    write(*,*)'max/min-tracer-flux',maxval(p_sfc_flx%forc_tracer),&
    & minval(p_sfc_flx%forc_tracer)
    write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
                                & minval(p_sfc_flx%forc_tracer(:,:,1))
  ENDIF


  END SUBROUTINE update_sfcflx
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
  TYPE (t_sea_ice),         INTENT(IN)    :: p_ice
  TYPE (t_atmos_fluxes),    INTENT(INOUT) :: Qatm


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
  
  INTEGER i

  INTEGER :: jc,jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:calc_atm_fluxes_from_bulk'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  rl_start_c = 1
  rl_end_c   = min_rlcell
  i_startblk_c = ppatch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)

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
      z_dragl   (jc,jb) = MAX ( 0.5e-3_wp,1.0e-3_wp*(0.8195_wp+0.0506_wp*z_fu10lim(jc,jb)  &
                        &-0.0009_wp*z_fu10lim(jc,jb)*z_fu10lim(jc,jb)) + z_dragl1(jc,jb)   &
                        &* (z_Tsurf(jc,jb)-p_as%tafo(jc,jb)) )
      z_dragl   (jc,jb) = MIN (z_dragl(jc,jb), 3.0E-3_wp)
      z_drags   (jc,jb) = 0.96_wp * z_dragl(jc,jb)
      Qatm%sensw(jc,jb) = z_drags(jc,jb)*z_rhoair(jc,jb)*cpa*p_as%fu10(jc,jb)             &
                        & * (p_as%tafo(jc,jb) -z_Tsurf(jc,jb))  *fr_fac
      Qatm%latw (jc,jb) = z_dragl(jc,jb)*z_rhoair(jc,jb)*L  *p_as%fu10(jc,jb)             &
                        & * (z_sphumida(jc,jb)-z_sphumidw(jc,jb))*fr_fac

      DO i = 1, p_ice%kice
        IF (p_ice% isice(jc,jb,i))THEN
          z_Tsurf(jc,jb) = p_ice%Tsurf(jc,jb,i)
          ! z_esti is calculated according to Buck Research Manuals, 1996 (see z_esta)
          z_esti     (jc,jb) = 611.15_wp*EXP( (23.036_wp-z_Tsurf(jc,jb)/333.7_wp) &
                             & *z_Tsurf(jc,jb)/(z_Tsurf(jc,jb) + 279.82_wp) )
          z_sphumidi (jc,jb) = 0.62197_wp*z_esti(jc,jb)/(p_as%pao(jc,jb)-0.37803_wp*z_esti(jc,jb))
          z_dragl    (jc,jb) = MAX (0.5e-3_wp, 1.0e-3_wp * (0.8195_wp+0.0506_wp*z_fu10lim(jc,jb) &
                             & -0.0009_wp*z_fu10lim(jc,jb) * z_fu10lim(jc,jb)) + z_dragl1(jc,jb) &
                             & * (z_Tsurf(jc,jb)-p_as%tafo(jc,jb)) )
          z_drags    (jc,jb) = 0.96_wp * z_dragl(jc,jb)

          Qatm%LWout (jc,jb,i) = emiss*StefBol * (z_Tsurf(jc,jb)+273.15_wp)**4
          Qatm%LWnet (jc,jb,i) = Qatm%LWin(jc,jb) - Qatm%LWout(jc,jb,i)
          Qatm%dLWdT (jc,jb,i) = - 4.0_wp * emiss*StefBol * (z_Tsurf(jc,jb) + 273.15_wp)**3
          Qatm%sens  (jc,jb,i) = z_drags(jc,jb) * z_rhoair(jc,jb)*cpa*p_as%fu10(jc,jb)&
                               & * (p_as%tafo(jc,jb) -z_Tsurf(jc,jb))   *fr_fac
          Qatm%lat   (jc,jb,i) = z_dragl(jc,jb) * z_rhoair(jc,jb)* L *p_as%fu10(jc,jb)&
                               & * (z_sphumida(jc,jb)-z_sphumidi(jc,jb))*fr_fac
    
          Qatm%dsensdT(jc,jb,i)= 0.96_wp*z_dragl1(jc,jb)*z_drags(jc,jb)*z_rhoair(jc,jb)&
                               & *cpa * p_as%fu10(jc,jb)       &
                               & * (p_as%tafo(jc,jb) - z_Tsurf(jc,jb)) *  fr_fac &
                               & -z_drags(jc,jb)*z_rhoair(jc,jb) *cpa*p_as%fu10(jc,jb)
          Qatm%dlatdT(jc,jb,i) = z_dragl1(jc,jb) * z_rhoair(jc,jb)*L * p_as%fu10(jc,jb)&
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
  !
  !> Takes thermal calc_atm_fluxes_from_bulk to calculate atmospheric surface fluxes:
  !  heat, freshwater and momentum.
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011).
  !
  SUBROUTINE update_sfcflx_from_atm_flx(ppatch, p_as, p_os, p_ice,Qatm, p_sfc_flx)
  TYPE(t_patch),                INTENT(in)    :: ppatch
  TYPE(t_atmos_for_ocean),      INTENT(IN)    :: p_as
  TYPE(t_hydro_ocean_state),    INTENT(IN)    :: p_os
  TYPE (t_sea_ice),             INTENT (IN)   :: p_ice
  TYPE (t_atmos_fluxes),        INTENT (INOUT):: Qatm
  TYPE(t_sfc_flx)                             :: p_sfc_flx

  !Local variables 
  REAL(wp) :: z_rho_w = 1.22_wp  !near surface air density [kg/m^3] cf. Large/Yeager, sect 4.1, p.17
  REAL(wp) :: z_C_d0, z_C_d1, z_C_d
  REAL(wp) :: z_norm, z_v, z_relax

  INTEGER :: jc,jb, i
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c
  REAL(wp):: z_evap(nproma,ppatch%nblks_c)
  REAL(wp):: z_Q_freshwater(nproma,ppatch%nblks_c)
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_sfcflx_from_atm_flx'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  rl_start_c = 1
  rl_end_c   = min_rlcell

  i_startblk_c = ppatch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)

  !Relaxation parameter from namelist for salinity.
  z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( ppatch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      DO i = 1, p_ice%kice
        !surface heat forcing as sum of sensible, latent, longwave and shortwave heat fluxes
        IF (p_ice% isice(jc,jb,i))THEN

          p_sfc_flx%forc_tracer(jc,jb,1)             &
          & =  Qatm%sens(jc,jb,i) + Qatm%lat(jc,jb,i)& ! Sensible + latent heat flux at ice surface
          & +  Qatm%LWnet(jc,jb,i)                   & ! net LW radiation flux over ice surface
          & +  Qatm%bot(jc,jb,i)                       ! Ocean heat flux at ice bottom 
                                                       ! liquid/solid  precipitation rate are zero

          !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
          z_evap(jc,jb) = Qatm%lat(jc,jb,i)/(Lsub*z_rho_w)

        ELSEIF(.NOT.p_ice% isice(jc,jb,i))THEN

          p_sfc_flx%forc_tracer(jc,jb,1)             &
          & =  Qatm%sensw(jc,jb) + Qatm%latw(jc,jb)  & ! Sensible + latent heat flux over water
          & +  Qatm%LWnetw(jc,jb)                    & ! net LW radiation flux over water
          & +  Qatm%SWin(jc,jb)                        ! incoming SW radiation flux
                                                       ! liquid/solid  precipitation rate are zero

         !This prepares freshwater flux calculation below; eq. (64) in Marsland et al.
          z_evap(jc,jb) = Qatm%latw(jc,jb)/(Lvap*z_rho_w)
        ENDIF
      END DO

      !calculate surface freshwater flux       
      !following MPI-OM as described in Marsland et al, formula (63)-(65)

      !calculate evaporation from latent heat flux and latent heat of vaporisation
      !This is (63) in Marsland et al.
      z_Q_freshwater(jc,jb) = (Qatm%rpreci(jc,jb) + Qatm%rprecw(jc,jb)) -  z_evap(jc,jb) !+River runof +glacial meltwater

      !Now the freshwater flux calculation is finished; this is (65) in Marslkand et al.
      !Relaxation of top layer salinity to observed salinity
      !
      p_sfc_flx%forc_tracer(jc,jb,2) = &
      &(v_base%del_zlev_m(1)+z_Q_freshwater(jc,jb))&
      &/v_base%del_zlev_m(1)                       &
      & + z_relax*(p_os%p_prog(nold(1))%tracer(jc,1,jb,2) - p_sfc_flx%forc_tracer_relax(jc,jb,2))


      !calculate wind stress    
      z_norm = sqrt(p_as%u(jc,jb)*p_as%u(jc,jb)+p_as%v(jc,jb)*p_as%v(jc,jb))

      !calculate drag coefficient for wind following 
      ! Kara, Rochford, Hurlburt, Air-Sea Flux Estimates And the 1997-1998 Enso Event
      ! Boundary-Layer Meteorology, 103, 439-458 (2002)
      !
      z_v = MAX(2.5_wp, MIN(p_as%fu10(jc,jb),32.5_wp))

      z_C_d0 = 1.0E-3_wp*(0.692_wp+0.071_wp*z_v-0.00070_wp*z_norm)
      z_C_d1 = 1.0E-3_wp*(0.083_wp-0.0054_wp*z_v-0.000093_wp*z_norm)
      z_C_d  = z_C_d0 + z_C_d1*(p_as%tafo(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1))
write(*,*)'final wind stress coeff',z_C_d
      p_sfc_flx%forc_wind_u(jc,jb) = z_rho_w*z_C_d*z_norm&
                                   &*(p_as%u(jc,jb)- p_os%p_diag%u(jc,1,jb))

      p_sfc_flx%forc_wind_v(jc,jb) = z_rho_w*z_C_d*z_norm&
                                   &*(p_as%v(jc,jb) - p_os%p_diag%v(jc,1,jb))
 
    END DO
  END DO

  IF(temperature_relaxation==1)THEN

     p_sfc_flx%forc_tracer(:,:, 1)=  z_relax                                    &
     & *( p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1) )

  ENDIF

  END SUBROUTINE update_sfcflx_from_atm_flx  
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
  SUBROUTINE update_sfcflx_analytical(p_patch, p_os, p_sfc_flx)

  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch
  TYPE(t_hydro_ocean_state)             :: p_os  
  TYPE(t_sfc_flx)                       :: p_sfc_flx
  !
  ! local variables
  INTEGER :: jc, jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  REAL(wp) :: zonal_str, z_tmp
  REAL(wp) :: z_thick_forc
  REAL(wp) :: z_lat, z_lon, z_lat_deg
  REAL(wp) :: z_forc_period = 1.0_wp !=1.0: single gyre
                                     !=2.0: double gyre
                                     !=n.0: n-gyre 
  REAL(wp) :: y_length               !basin extension in y direction in degrees
  REAL(wp) :: z_T_init(nproma,p_patch%nblks_c)
  REAL(wp) :: z_T_init_insitu(nproma,1,p_patch%nblks_c)
  REAL(wp) :: z_perlat, z_perlon, z_permax, z_perwid, z_relax, z_dst
  INTEGER  :: z_dolic
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:update_ho_sfcflx'
  !-------------------------------------------------------------------------
  rl_start_c   = 1
  rl_end_c     = min_rlcell
  i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

    SELECT CASE (itestcase_oce)

    CASE(27)
      CALL message(TRIM(routine), 'Testcase (27): Apply stationary wind forcing' )
      y_length = basin_height_deg * deg2rad
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN

            z_lat = p_patch%cells%center(jc,jb)%lat
            z_lon = p_patch%cells%center(jc,jb)%lon

            p_sfc_flx%forc_wind_u(jc,jb) =&
            & cos(z_forc_period*pi*(z_lat-y_length)/y_length) 
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


  !           IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
  !             !domain is from 30 deg north to -30 deg south, length=60 degrees 
  !             !zonal_str = cos(pi*(z_lat-60.0_wp*deg2rad)/(120.0_wp*deg2rad)) 
  !             zonal_str = cos(pi*(z_lat-60*deg2rad)/(60.0_wp*deg2rad))
  !             p_sfc_flx%forc_wind_cc(jc,jb)%x(1) = zonal_str*sin(z_lon)
  !             p_sfc_flx%forc_wind_cc(jc,jb)%x(2) = zonal_str*cos(z_lon)
  !             p_sfc_flx%forc_wind_cc(jc,jb)%x(3) = 0.0_wp
  !             CALL cvec2gvec(p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
  !                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
  !                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(3),&
  !                          & z_lon, z_lat,                      &
  !                          & p_sfc_flx%forc_wind_u(jc,jb),      &
  !                          & p_sfc_flx%forc_wind_v(jc,jb))
  ! ! write(*,*)'Danilovs Wind', jc,jb,p_sfc_flx%forc_wind_cc(jc,jb)%x(1:2), &
  ! ! &p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
  !            ELSE
  !              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
  !              p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
  !              p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
  !            ENDIF 
        END DO
      END DO

    CASE(30,34,35,36)

      !analytical forcing 
      IF(iforc_oce==11)THEN
        CALL message(TRIM(routine), &
          &  'Testcase (30,34,35; iforc=11): stationary wind forcing - u=cos(n*(lat-lat_0)/lat_0')

        !Latitudes in ICON vary from -pi/2 to pi/2
        y_length = basin_height_deg * deg2rad
        DO jb = i_startblk_c, i_endblk_c
          CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, &
          &                                i_startidx_c, i_endidx_c,rl_start_c,rl_end_c)
          DO jc = i_startidx_c, i_endidx_c
            z_lat = p_patch%cells%center(jc,jb)%lat
            z_lon = p_patch%cells%center(jc,jb)%lon
          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
              p_sfc_flx%forc_wind_u(jc,jb) = cos(z_forc_period*pi*(z_lat-y_length)&
                                           &/y_length) 
            ELSE
              p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
            ENDIF !write(*,*)'forc',jc,jb, z_lat, z_lat*180.0_wp/pi, p_sfc_flx%forc_wind_u(jc,jb) 
            p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp

            !Init cartesian wind
          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
              CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                             & p_sfc_flx%forc_wind_v(jc,jb),      &
                             & p_patch%cells%center(jc,jb)%lon,   &
                             & p_patch%cells%center(jc,jb)%lat,   &
                             & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                             & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                             & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
                           ! write(*,*)'sfc forcing', jc,jb,p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
            ELSE
              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
            ENDIF
          END DO
        END DO
        write(*,*)'max/min-analytical Forcing',&
        &maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)
      ENDIF
    CASE(31)
      ! 3d-gravity wave needs initialization only
      CONTINUE

    CASE(32)
      CALL message(TRIM(routine), &
      &  'Testcase (32) - Munk Gyre: stationary lat/lon wind forcing and relax. to T perturbation')
      y_length = basin_height_deg * deg2rad
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          z_T_init(jc,jb) = 20.0_wp- v_base%zlev_m(1)*15.0_wp/4000.0_wp

          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lon = p_patch%cells%center(jc,jb)%lon

          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN

             zonal_str = cos(z_forc_period*pi*z_lat-y_length/y_length)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(1) = zonal_str*sin(z_lon)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(2) = zonal_str*cos(z_lon)
             p_sfc_flx%forc_wind_cc(jc,jb)%x(3) = 0.0_wp
 
             CALL cvec2gvec(p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(3),&
                          & z_lon, z_lat,                      &
                          & p_sfc_flx%forc_wind_u(jc,jb),      &
                          & p_sfc_flx%forc_wind_v(jc,jb))
 
             ! Add temperature perturbation at new values
             z_perlat = basin_center_lat + 0.1_wp*basin_height_deg
             z_perlon = basin_center_lon + 0.1_wp*basin_width_deg 
             z_permax  = 0.1_wp
             z_perwid  =  10.0_wp

             z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)

             z_lat = p_patch%cells%center(jc,jb)%lat
             z_lon = p_patch%cells%center(jc,jb)%lon

             z_dolic = v_base%dolic_c(jc,jb)
             IF (z_dolic > 0) THEN

               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)

                IF(z_dst<=5.0_wp*deg2rad)THEN
                   z_T_init = z_T_init &
                   &        + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                   &        * sin(pi*v_base%zlev_m(1)/4000.0_wp)

                   !   write(*,*)'z init',jc,jb,p_os%p_prog(nold(1))%tracer(jc,1,jb,1),&
                   !   &z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                   !   & * sin(pi*v_base%zlev_m(1)/4000.0_wp)
                ENDIF
               ! up to here z_init is identically initialized than temperature

              !add local cold perturbation 
                IF(z_dst<=10.5_wp*deg2rad)THEN
                  z_T_init(jc,jb)= z_T_init(jc,jb) - exp(-(z_dst/(z_perwid*deg2rad))**2)
                ENDIF

                p_sfc_flx%forc_tracer_relax(jc,jb,1)=z_T_init(jc,jb)

                p_sfc_flx%forc_tracer(jc,jb, 1)=  z_relax              &!*v_base%del_zlev_m(1)     &
                & *( p_sfc_flx%forc_tracer_relax(jc,jb,1)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )

                ! write(123,*)'forcing',jc,jb,&
                ! &( p_sfc_flx%forc_tracer_relax(jc,jb,1)    &
                ! & -p_os%p_prog(nold(1))%tracer(jc,1,jb,1)),&
                ! &p_sfc_flx%forc_tracer_relax(jc,jb,1),&
                ! &p_sfc_flx%forc_tracer(jc,jb, 1)
             END IF
           ELSE
             p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
             p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
             p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
           ENDIF 
        END DO
      END DO
    write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
    & minval(p_sfc_flx%forc_tracer_relax)
    write(*,*)'max/min-tracer-flux',maxval(p_sfc_flx%forc_tracer),&
    & minval(p_sfc_flx%forc_tracer)
    write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)
    write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
                                  & minval(p_sfc_flx%forc_tracer(:,:,1))

! ! ! !-----------Old version of Forcing--------------------------------------------------
! ! !!------------Please retain, its also interesting------------------------------------
!!----------------An old version of init corresponds to this forcing--------------------
! !    CASE(32)
! !       CALL message(TRIM(routine), 'Testcase (32): Apply stationary wind forcing' )
! !       y_length = basin_height_deg * deg2rad
! !       DO jb = i_startblk_c, i_endblk_c    
! !         CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
! !          &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
! !         DO jc = i_startidx_c, i_endidx_c
! !           z_lat = p_patch%cells%center(jc,jb)%lat
! !           z_lon = p_patch%cells%center(jc,jb)%lon
! !           IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
! !             zonal_str = cos(z_forc_period*pi*z_lat-y_length/y_length)
! !             p_sfc_flx%forc_wind_cc(jc,jb)%x(1) = zonal_str*sin(z_lon)
! !             p_sfc_flx%forc_wind_cc(jc,jb)%x(2) = zonal_str*cos(z_lon)
! !             p_sfc_flx%forc_wind_cc(jc,jb)%x(3) = 0.0_wp
! !             CALL cvec2gvec(p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
! !                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
! !                          & p_sfc_flx%forc_wind_cc(jc,jb)%x(3),&
! !                          & z_lon, z_lat,                      &
! !                          & p_sfc_flx%forc_wind_u(jc,jb),      &
! !                          & p_sfc_flx%forc_wind_v(jc,jb))
! !             ! Add temperature perturbation at new values
! !            z_perlat = basin_center_lat + 0.1_wp*basin_height_deg!             !45.5_wp
! !            z_perlon =  0.1_wp*basin_width_deg                                 !4.5_wp
! !            z_permax  = 10.0_wp!20.1_wp
! !            z_perwid  =  5.0_wp!1.5_wp
! !            z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
! ! 
! !             z_dolic = v_base%dolic_c(jc,jb)
! !             IF (z_dolic > 0) THEN
! !               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)
! ! 
! !               !init temperature
! !               z_T_init(jc,jb) = 20.0_wp&
! !               & - v_base%zlev_i(1)*15.0_wp/v_base%zlev_i(z_dolic+1)
! ! 
! !                !add local hot perturbation 
! ! !              IF(z_dst<=3.5_wp*deg2rad)THEN
! !                 z_T_init(jc,jb)= z_T_init(jc,jb)  &
! !                 &   + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
! !                 &   * sin(pi*v_base%zlev_m(1)/v_base%zlev_i(z_dolic+1))
! ! !              ENDIF
! !               !Add local cold perturbation
! !               !IF(z_dst<=5.0_wp*deg2rad)THEN
! !               z_T_init(jc,jb) = z_T_init(jc,jb)     &
! !               &   - z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2)
! !               p_sfc_flx%forc_tracer_relax(jc,jb,1)=z_T_init(jc,jb)
! !               p_sfc_flx%forc_tracer(jc,jb, 1)=z_relax*v_base%del_zlev_i(1)*&
! !               &          (p_sfc_flx%forc_tracer_relax(jc,jb,1)&
! !               &         -p_os%p_prog(nold(1))%tracer(jc,1,jb,1))
! !               !ENDIF 
! !             END IF
! !   ! write(*,*)'Danilovs Wind', jc,jb,p_sfc_flx%forc_wind_cc(jc,jb)%x(1:2), &
! !   ! &p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
! !            ELSE
! !              p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
! !              p_sfc_flx%forc_wind_u(jc,jb)       = 0.0_wp
! !              p_sfc_flx%forc_wind_v(jc,jb)       = 0.0_wp
! !            ENDIF 
! !         END DO
! !       END DO
! !       write(*,*)'max/min-Wind-Forcing',maxval(p_sfc_flx%forc_wind_u), minval(p_sfc_flx%forc_wind_u)
! !       write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
! !                                   & minval(p_sfc_flx%forc_tracer(:,:,1))
! ! ! !-----------End of Old version of Forcing-------------------------------------------

    CASE (33)
      CALL message(TRIM(routine), &
        &  'Testcase (33): stationary temperature relaxation - latitude dependent')
      z_relax = relaxation_param/(30.0_wp*24.0_wp*3600.0_wp)
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          !transer to latitude in degrees
          z_lat_deg = p_patch%cells%center(jc,jb)%lat*rad2deg
          IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN

              !constant salinity
               IF(no_tracer==2)THEN
                 p_os%p_prog(nold(1))%tracer(jc,1,jb,2) = 35.0_wp
               ENDIF
              IF(abs(z_lat_deg)>=40.0_wp)THEN

                z_T_init_insitu(jc,1,jb) = 5.0_wp
                z_T_init(jc,jb)&
                &= convert_insitu2pot_temp_func(z_T_init_insitu(jc,1,jb),&
                                                &p_os%p_prog(nold(1))%tracer(jc,1,jb,2),&
                                                &sfc_press_bar)
              ELSEIF(abs(z_lat_deg)<=20.0_wp)THEN

                z_T_init_insitu(jc,1,jb) = 30.0_wp
                z_T_init(jc,jb)&
                &= convert_insitu2pot_temp_func(z_T_init_insitu(jc,1,jb),&
                                                &p_os%p_prog(nold(1))%tracer(jc,1,jb,2),&
                                                &sfc_press_bar)

              ELSEIF(abs(z_lat_deg)<40.0_wp .AND. abs(z_lat_deg)>20.0_wp)THEN
                z_tmp = pi*((abs(z_lat_deg) -20.0_wp)/20.0_wp)

                z_T_init_insitu(jc,1,jb) =  5.0_wp&
                                        & + 0.5_wp*25.0_wp*(1.0_wp+cos(z_tmp))
                z_T_init(jc,jb)&
                &= convert_insitu2pot_temp_func(z_T_init_insitu(jc,1,jb),&
                                               &p_os%p_prog(nold(1))%tracer(jc,1,jb,2),&
                                               &sfc_press_bar)
              ENDIF
          ENDIF

          p_sfc_flx%forc_tracer_relax(jc,jb,1)=z_T_init(jc,jb)

          p_sfc_flx%forc_tracer(jc,jb, 1)=  &  !0.0_wp!               &
          &          -z_relax               &  !*v_base%del_zlev_m(1) &
          &          *( p_sfc_flx%forc_tracer_relax(jc,jb,1)     &
          &           -p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )
        END DO
      END DO  
   write(*,*)'max/min-tracer-diff',&
   &maxval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1)),&
    & minval(p_sfc_flx%forc_tracer_relax(:,:,1)-p_os%p_prog(nold(1))%tracer(:,1,:,1))

    write(*,*)'max/min-tracer-relaxation',maxval(p_sfc_flx%forc_tracer_relax),&
    & minval(p_sfc_flx%forc_tracer_relax)
    write(*,*)'max/min-tracer-flux',maxval(p_sfc_flx%forc_tracer),&
    & minval(p_sfc_flx%forc_tracer)
    write(*,*)'max/min-Temp-Flux',maxval(p_sfc_flx%forc_tracer(:,:,1)),&
                                  & minval(p_sfc_flx%forc_tracer(:,:,1))

    CASE DEFAULT
      CALL message(TRIM(routine), 'STOP: Analytical Forcing for this testcase not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION NOT SUPPORTED - TERMINATE')
    END SELECT


    !Final modification of surface wind forcing according to surface boundary condition 
    SELECT CASE (i_bc_veloc_top) !The value of 'top_boundary_condition' was set in namelist

    CASE (0)

      CALL message (TRIM(routine),'ZERO top velocity boundary conditions chosen')

      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
          &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
            p_sfc_flx%forc_wind_u(jc,jb)      = 0.0_wp
            p_sfc_flx%forc_wind_v(jc,jb)      = 0.0_wp
            p_sfc_flx%forc_wind_cc(jc,jb)%x(:)= 0.0_wp
        END DO
      END DO

    CASE (1)!Forced by wind stored 

      z_thick_forc = v_base%del_zlev_m(1)      !z_thick_forc= 1.0_wp
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
          &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary)THEN

            p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff*p_sfc_flx%forc_wind_u(jc,jb)&
              &/(rho_ref*z_thick_forc)

            p_sfc_flx%forc_wind_v(jc,jb) = wstress_coeff*p_sfc_flx%forc_wind_v(jc,jb)&
              &/(rho_ref*z_thick_forc)
           !write(*,*)'top bc gc:', jc,jb, top_bc_u_c(jc,jb), top_bc_v_c(jc,jb)

            p_sfc_flx%forc_wind_cc(jc,jb)%x =&
              & wstress_coeff*p_sfc_flx%forc_wind_cc(jc,jb)%x&
              &/(rho_ref*z_thick_forc)

!             p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff*( p_sfc_flx%forc_wind_u(jc,jb)    &
!               &               - p_os%p_diag%u(jc,1,jb) )&
!               &               /(rho_ref*z_thick_forc)
!             p_sfc_flx%forc_wind_v(jc,jb) = wstress_coeff*( p_sfc_flx%forc_wind_v(jc,jb)    &
!               &               - p_os%p_diag%v(jc,1,jb) )&
!               &               /(rho_ref*z_thick_forc)
! 
!             p_sfc_flx%forc_wind_cc(jc,jb)%x =&
!             & wstress_coeff*p_sfc_flx%forc_wind_cc(jc,jb)%x&
!             &/(rho_ref*z_thick_forc)
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb)    =0.0_wp
            p_sfc_flx%forc_wind_v(jc,jb)    =0.0_wp
            p_sfc_flx%forc_wind_cc(jc,jb)%x =0.0_wp
          ENDIF
        END DO
      END DO

    CASE (2)!forced by difference between wind field in p_os%p_aux%bc_top_veloc and ocean velocity at top layer

      z_thick_forc = v_base%del_zlev_m(1)      !z_thick_forc= 1.0_wp
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,  &
          &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary)THEN

            p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff*( p_sfc_flx%forc_wind_u(jc,jb)    &
              &               - p_os%p_diag%u(jc,1,jb) )&
              &               /(rho_ref*z_thick_forc)
            p_sfc_flx%forc_wind_v(jc,jb) = wstress_coeff*( p_sfc_flx%forc_wind_v(jc,jb)    &
              &               - p_os%p_diag%v(jc,1,jb) )&
              &               /(rho_ref*z_thick_forc)

            p_sfc_flx%forc_wind_cc(jc,jb)%x =&
            & wstress_coeff*p_sfc_flx%forc_wind_cc(jc,jb)%x&
            &/(rho_ref*z_thick_forc)
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb)    =0.0_wp
            p_sfc_flx%forc_wind_v(jc,jb)    =0.0_wp
            p_sfc_flx%forc_wind_cc(jc,jb)%x =0.0_wp
          ENDIF
        END DO
      END DO
    END SELECT

  END SUBROUTINE update_sfcflx_analytical

  !-------------------------------------------------------------------------

  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


END MODULE mo_oce_bulk
