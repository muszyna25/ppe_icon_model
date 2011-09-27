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
USE mo_ocean_nml,           ONLY: no_tracer, itestcase_oce, iforc_oce, analyt_forc, &
  &                               wstress_coeff, iforc_stat_oce, basin_height_deg
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch
USE mo_oce_index,           ONLY: print_mxmn, jkc, jkdim, ipl_src
USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base
USE mo_exception,           ONLY: finish, message !, message_text
USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary!,    &
! &                               land, sea, boundary,                  &
! &                               min_rlcell, min_rledge, min_rlvert
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec, cvec2gvec
!USE mo_param_ice !,           ONLY: kice
USE mo_sea_ice,             ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
IMPLICIT NONE

PRIVATE


CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface

! public subroutines
PUBLIC  :: init_sfcflx


CONTAINS


  !-------------------------------------------------------------------------
  !
  !>
  !! Initialization of stationary surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2011-09)
  !
  SUBROUTINE init_sfcflx(ppatch, p_sfc_flx)
  !
  TYPE(t_patch), INTENT(in) :: ppatch
  TYPE(t_sfc_flx)           :: p_sfc_flx

  ! Local variables
  INTEGER :: nblks_c, jc, jb
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: rl_start_c, rl_end_c

  REAL(wp) :: z_lat, z_lon, y_length, y_center
  REAL(wp) :: z_forc_period = 1.0_wp !=1.0: single gyre
                                     !=2.0: double gyre
                                     !=n.0: n-gyre 
  REAL(wp) :: z_c(nproma,1,ppatch%nblks_c)

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:init_ho_sfcflx'

  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  rl_start_c = 1
  rl_end_c = min_rlcell

  i_startblk_c = ppatch%cells%start_blk(rl_start_c,1)
  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)

  nblks_c = ppatch%nblks_c

  ! analytical forcing 

  IF (iforc_oce == ANALYT_FORC) THEN

    IF (itestcase_oce == 27 .OR. itestcase_oce == 29) iforc_stat_oce = 1

    SELECT CASE (iforc_stat_oce)

    CASE (0)

      CALL message(TRIM(routine), &
        &  'iforc_stat_oce=0: no stationary wind forcing applied' )

    CASE (1)

      CALL message(TRIM(routine), 'Testcase (27,29): Apply stationary wind forcing' )
      y_length = basin_height_deg * deg2rad
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN

            z_lat = ppatch%cells%center(jc,jb)%lat
            z_lon = ppatch%cells%center(jc,jb)%lon

            p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff * &
            & cos(z_forc_period*pi*(z_lat-y_length)/y_length) 
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
          ENDIF !write(*,*)'forc',jc,jb, z_lat, z_lat*180.0_wp/pi, p_sfc_flx%forc_wind_u(jc,jb) 
          p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp

          !Init cartesian wind
          CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
                         & p_sfc_flx%forc_wind_v(jc,jb),&
                         & ppatch%cells%center(jc,jb)%lon,&
                         & ppatch%cells%center(jc,jb)%lat,&
                         & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                         & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                         & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
  !write(*,*)'sfc forcing', jc,jb,p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)

        END DO
      END DO

    CASE (2)

      CALL message(TRIM(routine), &
        &  'iforc_stat_oce=2: stationary wind forcing - u=cos(n*(lat-lat_0)/lat_0)')
      
      ! Latitudes vary from -pi/2 to pi/2
      y_length = basin_height_deg * deg2rad
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c( ppatch, jb, i_startblk_c, i_endblk_c, &
        &                                i_startidx_c, i_endidx_c,rl_start_c,rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon
          IF (v_base%lsm_oce_c(jc,1,jb)<=sea_boundary) THEN
            p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * cos(z_forc_period*pi*(z_lat-y_length)&
              &                                / y_length) 
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
          ENDIF
          !write(*,*)'forc',jc,jb, z_lat, z_lat*180.0_wp/pi, p_sfc_flx%forc_wind_u(jc,jb) 
          p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
      
          !Init cartesian wind
          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                           & p_sfc_flx%forc_wind_v(jc,jb),      &
                           & ppatch%cells%center(jc,jb)%lon,   &
                           & ppatch%cells%center(jc,jb)%lat,   &
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
          ! write(*,*)'sfc forcing', jc,jb,p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
          ELSE
            p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
          ENDIF
        END DO
      END DO

    CASE (3)

      CALL message(TRIM(routine), &
        &  'iforc_stat_oce=3: stationary wind forcing - u=cos(n*lat/lat_0)')
      
      ! Use here global scale:
      !y_length = basin_height_deg * deg2rad
      !y_center = basin_center_lat * deg2rad
      y_length = 180.0_wp * deg2rad
      y_center = -60.0_wp * deg2rad
      z_forc_period = 3.0_wp
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c( ppatch, jb, i_startblk_c, i_endblk_c, &
          &                             i_startidx_c, i_endidx_c,rl_start_c,rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon
          IF (v_base%lsm_oce_c(jc,1,jb)<=sea_boundary) THEN
            p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * cos(z_forc_period*pi*(z_lat-y_center)&
              &                                / y_length) 
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
          ENDIF
          p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp
      
          !Init cartesian wind
          IF(v_base%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                           & p_sfc_flx%forc_wind_v(jc,jb),      &
                           & ppatch%cells%center(jc,jb)%lon,   &
                           & ppatch%cells%center(jc,jb)%lat,   &
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
          ELSE
            p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
          ENDIF
        END DO
      END DO

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Stationary Analytical Forcing not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN STATIONARY FORCING OPTION NOT SUPPORTED - TERMINATE')

    END SELECT

    ipl_src=0  ! output print level (1-5, fix)
    z_c(:,1,:)=p_sfc_flx%forc_wind_u(:,:)
    CALL print_mxmn('analytical forcing u',1,z_c(:,:,:),1,ppatch%nblks_c,'per',ipl_src)
    z_c(:,1,:)=p_sfc_flx%forc_wind_v(:,:)
    CALL print_mxmn('analytical forcing v',1,z_c(:,:,:),1,ppatch%nblks_c,'per',ipl_src)

  END IF

  END SUBROUTINE init_sfcflx


END MODULE mo_oce_forcing
