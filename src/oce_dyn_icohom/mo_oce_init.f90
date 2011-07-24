!>
!! Contains the implementation of the initial conditions for the hydrostatic ocean model.
!!
!! Contains the implementation of the initial conditions for the hydrostatic ocean model.
!! This module controls the initial conditions as well as the initialisation of
!! test cases, the top and bottom boundary conditions, and the structure of the
!! forcing quantities.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!
!! @par Revision History
!! Initial version  by Peter Korn (MPI-M)  (2006).
!! Modified by Stephan Lorenz     (MPI-M)  (2010-06).
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_oce_init
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
USE mo_kind,               ONLY: wp
USE mo_physical_constants, ONLY: re, rre, omega, rgrav
USE mo_math_constants
USE mo_ocean_nml,          ONLY: iswm_oce, n_zlev, no_tracer , &
   &                             itestcase_oce,  &
   &                             basin_center_lat, basin_center_lon,idisc_scheme,&
   &                             basin_height_deg,  basin_width_deg
!  &                             iforc_oce, inoforcing, analyt_forc, core_forc, full_forc, &
!  &                             core_annwind, wstress_coeff
USE mo_impl_constants,     ONLY: max_char_length, sea, sea_boundary, & !success,               &
  & min_rlcell, min_rledge,  &
  & oce_testcase_zero, oce_testcase_init, oce_testcase_file                             
!  &                              land , boundary
USE mo_dynamics_config,    ONLY: nold,nnew
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e!, get_indices_v
USE mo_exception,          ONLY: finish, message!, message_text
USE mo_model_domain,       ONLY: t_patch
USE mo_math_utilities,     ONLY: gc2cc,cvec2gvec, t_cartesian_coordinates,&
  &                              t_geographical_coordinates
USE mo_oce_forcing,         ONLY: t_ho_sfc_flx!, t_ho_core_forc, t_ho_full_forc,  &
!  &                              t_ho_sfc_flx, f_coreforc, f_fullforc,          &
!  &                              construct_ho_sfcflx, destruct_ho_sfcflx
USE mo_oce_state,         ONLY: t_hydro_ocean_state
USE mo_scalar_product,            ONLY: map_cell2edges, map_edges2cell, map_edges2edges, &
  &                                     calc_scalar_product_for_veloc
IMPLICIT NONE
PRIVATE

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

! public interface
!
! public subroutines
!PUBLIC :: init_ho_forcing_old
!PUBLIC :: destruct_ho_forcing
PUBLIC :: init_ho_testcases

!private implementation
!PRIVATE :: init_ho_forcing_analytical
!PRIVATE :: init_ho_core_forcing
!PRIVATE :: init_ho_full_forcing
!PRIVATE :: init_ho_forcing_U_10_CORE
!PRIVATE :: init_testcase
!PRIVATE :: init_four_layer_Stommel
!PRIVATE :: init_twenty_layer_box
!PRIVATE :: init_twenty_layer_ocean

REAL(wp), PARAMETER :: aleph = 0.0_wp 
REAL (wp), PARAMETER :: u0 =(2.0_wp*pi*re)/(12.0_wp*24.0_wp*3600.0_wp)
CONTAINS

! 
! !-------------------------------------------------------------------------
  !>
  !! Initialization of test cases for the hydrostatic ocean model.
  !! Currently only some simple test value are set.
  !! Finally the prognostic state should be initialized
  !! from some restart file.
  !
  !! @par Revision History
  !! Developed  by Peter Korn, MPI-M, 2006-08
  !
  SUBROUTINE init_ho_testcases(ppatch, p_os, p_sfc_flx)
  TYPE(t_patch)                     :: ppatch
  TYPE(t_hydro_ocean_state), TARGET :: p_os 
  TYPE(t_ho_sfc_flx)                :: p_sfc_flx
  ! Local Variables
    INTEGER :: jb, jc, je, jk
    INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start, rl_end_e,rl_end_c
    INTEGER :: z_dolic
    REAL(wp)::  z_lat, z_lon 
    REAL(wp):: z_dst, z_lat_deg, z_lon_deg, z_tmp
    REAL(wp):: z_perlon, z_perlat, z_permax, z_perwid
    !TYPE(t_cartesian_coordinates)    :: p_x 
    !TYPE(t_geographical_coordinates) :: p_pos
    !REAL(wp), PARAMETER :: tprof(20)=(/ 18.13, 17.80, 17.15, 16.09, 15.04, 13.24, 11.82, 9.902, &
    !  &   8.484, 7.341, 5.727, 4.589, 3.807, 3.062, 2.481, 2.194, 1.789, 1.266, 1.070, 0.9211 /)

    !REAL(wp), PARAMETER :: sprof(20)=(/  34.699219, 34.798244, 34.904964, 34.976841, 35.027084, &
    !  & 35.026825, 34.960835, 34.862324, 34.752468, 34.656761, 34.596603, 34.594128, 34.628601, &
    !  & 34.678772, 34.717495, 34.738304, 34.741512, 34.738205, 34.729176, 34.723465 /)

    REAL(wp) , PARAMETER :: tprof_4layerStommel(4) = (/20.0_wp,10.0_wp,8.0_wp,6.0_wp/)
    REAL(wp) , PARAMETER :: sprof_4layerStommel(4) = &
    &(/34.699219_wp, 34.798244_wp, 34.904964_wp, 34.976841_wp/)
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_init:init_ho_testcases'
  !-------------------------------------------------------------------------
  CALL message (TRIM(routine), 'start')

  rl_start = 1
  rl_end_c = min_rlcell
  rl_end_e = min_rledge
  i_startblk_e = ppatch%edges%start_blk(rl_start,1)
  i_endblk_e   = ppatch%edges%end_blk(rl_end_e,1)
  i_startblk_c = ppatch%cells%start_blk(rl_start,1)
  i_endblk_c   = ppatch%cells%end_blk(rl_end_c,1)    

  !IF shallow-water option is NOT selected then)
  IF ( iswm_oce /= 1 )THEN

    SELECT CASE (itestcase_oce)

    CASE (oce_testcase_zero)
      CALL message(TRIM(routine), 'you have selected the "no-testcase" option')
    CASE (oce_testcase_init)

    CASE (oce_testcase_file)
      CALL finish(TRIM(routine), 'Initialization from file NOT SUPPORTED YET - TERMINATE')
      !CALL init_from_file(ppatch)

    CASE (30)
        CALL message(TRIM(routine), 'Simple Initialization of testcases (30)')
        CALL message(TRIM(routine), ' - here: 4-layer Stommel Testcase for T and S')
        !CALL init_testcase(ppatch, itestcase)
        ! #slo# 2011-05-11 - dolic set to 4 no longer nec.
        !ppatch%patch_oce%dolic_c(:,:)      = 4
        !init temperature and salinity with vertical profiles
        DO jk=1,n_zlev
          DO jb = i_startblk_c, i_endblk_c    
            CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                      & i_startidx_c, i_endidx_c, rl_start, rl_end_c)

            DO jc = i_startidx_c, i_endidx_c
               !z_lat = ppatch%cells%center(jc,jb)%lat
               !z_lon = ppatch%cells%center(jc,jb)%lon

              IF ( ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
                !Temperature
                p_os%p_prog(nold(1))%tracer(jc,jk,jb,1) = tprof_4layerStommel(jk)
                !Salinity
                p_os%p_prog(nold(1))%tracer(jc,jk,jb,2) = sprof_4layerStommel(jk)

              !Temperature
                p_os%p_prog(nnew(1))%tracer(jc,jk,jb,1) = tprof_4layerStommel(jk)
                !Salinity
                p_os%p_prog(nnew(1))%tracer(jc,jk,jb,2) = sprof_4layerStommel(jk)
                !p_os%p_prog(nold(1))%h(jc,jb) = 1.0E-7*test5_h( z_lon, z_lat, 0.0_wp)
              ENDIF
            END DO
          END DO
        END DO

    CASE (31)
      CALL message(TRIM(routine), 'Simple Initialization of testcases (31)')
      CALL message(TRIM(routine), ' - here: external gravity wave')
      !Flat surface of the ocean
      p_os%p_prog(nold(1))%h(:,:) = 0.0_wp

      !Ocean at rest 
      p_os%p_prog(nold(1))%vn(:,:,:) = 0.0_wp

      ! #slo# 2011-01-07: init elevation for simple 3-d SW-like test mode
      !all other prognostic variables: vn, s, t are initialized identical zero
      !CALL message(TRIM(routine), 'Simple Initialization of h')

      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon

          ! #slo#: simple elevation between 30W and 30E (pi/3.)
          IF ( ppatch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
            p_os%p_prog(nold(1))%h(jc,jb) = 10.0_wp * &
              &    sin(z_lon*6.0_wp) * cos(z_lat*3.0_wp)
          ELSE
            p_os%p_prog(nold(1))%h(jc,jb) = 0.0_wp
          ENDIF
        END DO
      END DO

      WRITE(*,*)'MIN/MAX:h-init:',&
      &        minval(p_os%p_prog(nold(1))%h(:,:)), &
      &        maxval(p_os%p_prog(nold(1))%h(:,:))

      !init temperature and salinity with vertical profiles
      DO jk=1,n_zlev
        DO jb = i_startblk_c, i_endblk_c    
          CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                          & i_startidx_c, i_endidx_c, rl_start, rl_end_c)
  
          DO jc = i_startidx_c, i_endidx_c
            IF ( ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
              !Temperature
              p_os%p_prog(nold(1))%tracer(jc,jk,jb,1) = tprof_4layerStommel(jk)
              !Salinity
              p_os%p_prog(nold(1))%tracer(jc,jk,jb,2) = sprof_4layerStommel(jk)
            ENDIF
          END DO
        END DO
      END DO 

    CASE (32) !from Sergy Danilov
      CALL message(TRIM(routine), 'Simple Initialization of testcases (32)')
      CALL message(TRIM(routine), ' - here: Danilovs Munk gyre flow')

      !p_pos%lon = 0.0_wp
      !p_pos%lat = 0.0_wp

      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c
          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon

          z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
          IF (z_dolic > 0) THEN
            ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
            ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
            ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
            ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375
            DO jk = 1, z_dolic
               p_os%p_prog(nold(1))%tracer(jc,jk,jb,1) &
              & = 20.0_wp-ppatch%patch_oce%zlev_m(jk)*15.0_wp/ppatch%patch_oce%zlev_i(z_dolic+1)
              ! p_os%p_prog(nold(1))%tracer(jc,jk,jb,1) &
              !   & = 20.0_wp -  ppatch%patch_oce%zlev_m(jk)*15.0_wp/4000.0_wp
            END DO
          END IF
        END DO
      END DO

      ! #slo# 2011-05-23 - new formulation and suggested values for perturbation radius
      !
      ! Add temperature perturbation centered at 45.5 N and 4.5 W - old values
      z_perlat = basin_center_lat + 0.1_wp*basin_height_deg!             !45.5_wp
      z_perlon =  0.1_wp*basin_width_deg                                 !4.5_wp
      z_permax  = 10.0_wp!20.1_wp
      z_perwid  =  5.0_wp!1.5_wp

      ! Add temperature perturbation at new values - 35N; 10W
      !z_perlat  = 35.0_wp
      !z_perlon  = 10.0_wp
      ! Next update 2011-05-24: due to Danilov the perturbation should be -1 Kelvin, width 3.0
      ! 05-25: max and width larger: -2.0 and 5.0
      !z_permax  = -2.0_wp
      !z_perwid  =  5.0_wp

      IF (no_tracer > 0 ) THEN

        DO jb = i_startblk_c, i_endblk_c    
          CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
           &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)

          DO jc = i_startidx_c, i_endidx_c
            z_lat = ppatch%cells%center(jc,jb)%lat
            z_lon = ppatch%cells%center(jc,jb)%lon

            z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
            IF (z_dolic > 0) THEN
              ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
              ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
              ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
              ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375

              z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)
              !write(123,*)'zdist',z_lat,z_lon,z_dst,10.5_wp*deg2rad
              !IF(z_dst<=25.5_wp*deg2rad)cycle
              ! at distance > 25.5 degrees: 
              !  e.g. at 30 deg distance the added perturbation would be ~ exp(-400) ~ 0.0
              ! Now without cycle in loop - perturbation is very small at z_dst>10 deg
              !  e.g. at 3 deg distance is
              !   T(jk=1)=19.0625+20.1*exp(-4)*sin(pi* 250/4000) = 19.06 + 20.1*0.18*0.06 = 19.28
              !   T(jk=4)=13.4375+20.1*exp(-4)*sin(pi*1750/4000) = 13.44 + 20.1*0.18*0.42 = 15.00
              !Local hot perturbation
              !IF(z_dst<=5.0_wp*deg2rad)THEN
              DO jk = 1, z_dolic
                p_os%p_prog(nold(1))%tracer(jc,jk,jb,1) =          &
                & p_os%p_prog(nold(1))%tracer(jc,jk,jb,1)          &
                &   + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                &   * sin(pi*ppatch%patch_oce%zlev_m(jk)/ppatch%patch_oce%zlev_i(z_dolic+1))
              END DO
              !ENDIF
            END IF
          END DO
        END DO

!         !After hot spot now a cool spot at a slightly different location
!         ! Add temperature perturbation at new values - 35N; 10W
!         z_perlat = basin_center_lat - 0.1_wp*basin_height_deg!             !45.5_wp
!         z_perlon =  -0.1_wp*basin_width_deg                                 !4.5_wp
!         z_permax  = 10.0_wp!20.1_wp
!         z_perwid  =  5.0_wp!1.5_wp
! !         z_perlat  = 25.0_wp
! !         z_perlon  = 8.0_wp
! !         z_permax  = -2.0_wp
! !         z_perwid  =  5.0_wp
!         DO jb = i_startblk_c, i_endblk_c
!           CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
!            &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)
!           DO jc = i_startidx_c, i_endidx_c
!             z_lat = ppatch%cells%center(jc,jb)%lat
!             z_lon = ppatch%cells%center(jc,jb)%lon
!             z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
!             IF (z_dolic > 0) THEN
!               ! jk=1:  250m  T= 20 - 0.9375 = 19.0625
!               ! jk=2:  750m  T= 20 - 2.8125 = 17.1875
!               ! jk=3: 1250m  T= 20 - 4.6875 = 15.3125
!               ! jk=4: 1750m  T= 20 - 6.5625 = 13.4375
!               z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)
!               !write(123,*)'zdist',z_lat,z_lon,z_dst,10.5_wp*deg2rad
!               !IF(z_dst<=25.5_wp*deg2rad)cycle
!               ! at distance > 25.5 degrees: 
!               !  e.g. at 30 deg distance the added perturbation would be ~ exp(-400) ~ 0.0
!               ! Now without cycle in loop - perturbation is very small at z_dst>10 deg
!               !  e.g. at 3 deg distance is
!               !   T(jk=1)=19.0625+20.1*exp(-4)*sin(pi* 250/4000) = 19.06 + 20.1*0.18*0.06 = 19.28
!               !   T(jk=4)=13.4375+20.1*exp(-4)*sin(pi*1750/4000) = 13.44 + 20.1*0.18*0.42 = 15.00   
!               DO jk = 1, z_dolic
!                p_os%p_prog(nold(1))%tracer(jc,jk,jb,1) =          &
!                 & p_os%p_prog(nold(1))%tracer(jc,jk,jb,1)          &
!                 &   - z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
!                 &   * sin(pi*ppatch%patch_oce%zlev_m(jk)/ppatch%patch_oce%zlev_i(z_dolic+1))
!               ENDDO
!             END IF
!           END DO
!         END DO

        ! Add elevation perturbation at new values - 35N; 10W
        ! not clear yet
  !     DO jb = i_startblk_c, i_endblk_c    
  !       CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
  !        &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)
  !       DO jc = i_startidx_c, i_endidx_c
       
  !         z_lat = ppatch%cells%center(jc,jb)%lat
  !         z_lon = ppatch%cells%center(jc,jb)%lon
  !         z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
  !         IF (z_dolic > 0) THEN
  !           z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)
  !           !IF(z_dst<=15.5_wp*deg2rad) cycle
  !           IF(z_dst>=15.5_wp*deg2rad) cycle
  !           p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)&
  !                                        &-0.3_wp*exp(-(z_dst/(2.2_wp*deg2rad))**2)
  !         ENDIF
  !       END DO
  !     END DO

      END IF ! no_tracer > 0

    CASE (33) !collpasing density front testcase, taken from Stuhne-Peltier (JCP, 2006)

      p_os%p_prog(nold(1))%tracer(:,:,:,1) = 20.0_wp

      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c

          !latitude given in radians
          z_lat = ppatch%cells%center(jc,jb)%lat
          !transer to latitude in degrees
          z_lat_deg = z_lat*rad2deg
          !Impose emperature profile. Profile
          !depends on latitude only and is uniform across
          !all vertical layers
          IF(abs(z_lat_deg)>=40.0_wp)THEN

            p_os%p_prog(nold(1))%tracer(jc,1:n_zlev,jb,1) = 5.0_wp

          ELSEIF(abs(z_lat_deg)<=20.0_wp)THEN

            p_os%p_prog(nold(1))%tracer(jc,1:n_zlev,jb,1) = 30.0_wp

          ELSEIF(abs(z_lat_deg)<40.0_wp .AND. abs(z_lat_deg)>20.0_wp)THEN
            z_tmp = (abs(z_lat_deg) -20.0_wp)/20.0_wp
            p_os%p_prog(nold(1))%tracer(jc,1:n_zlev,jb,1) = 5.0_wp&
                                    & + 0.5_wp*25.0_wp*(1.0_wp+cos(z_tmp))
          ENDIF
          !write(90,*)'lat-degrees', jc,jb,z_lat, z_lat_deg, p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
        END DO
      END DO

   CASE (34)

      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
         &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c
          IF(ppatch%patch_oce%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
          !latitude given in radians
          z_lon = ppatch%cells%center(jc,jb)%lon
          !transer to latitude in degrees
          z_lon_deg = z_lon*rad2deg
          !Impose emperature profile. Profile
          !depends on latitude only and is uniform across
          !all vertical layers
          IF(z_lon_deg>=basin_center_lon*rad2deg)THEN

            p_os%p_prog(nold(1))%tracer(jc,1:n_zlev,jb,1) = 5.0_wp
          ELSE
            p_os%p_prog(nold(1))%tracer(jc,1:n_zlev,jb,1) = 30.0_wp
          ENDIF
          ENDIF
        END DO
      END DO


    CASE DEFAULT
     CALL finish(TRIM(routine), 'CHOSEN INITIALIZATION NOT SUPPORTED - TERMINATE')
    END SELECT


!IF shallow-water option is selected then chosse between 2 options)
ELSEIF( iswm_oce == 1 )THEN

  SELECT CASE (itestcase_oce)

    CASE (oce_testcase_zero)

      CALL message(TRIM(routine), 'you have selected the "no-testcase" option')

    CASE (25)

      CALL message(TRIM(routine), 'Shallow-Water-Testcase (25)')
      CALL message(TRIM(routine), ' - here: h and bathy for solid body rotation (Laeuter Test)')

      ppatch%patch_oce%lsm_oce_c(:,:,:) = sea
      ppatch%patch_oce%lsm_oce_e(:,:,:) = sea
      !init height 
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                    & i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c
          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon
  
          p_os%p_prog(nold(1))%h(jc,jb) = test_usbr_h( z_lon, z_lat, 0.0_wp)

          ! #slo# - bathymetry not taken into account due to constant bottom cell depth
          !       - use ppatch%patch_oce%zlev_i(dolic(jc,jb)) instead
          ppatch%patch_oce%bathymetry_c(jc,jb) = test_usbr_oro( z_lon, z_lat, 0.0_wp )

          ! write(*,*)'h orig, bathy_c:', z_lon, z_lat,p_os%p_prog(nold(1))%h(jc,jb)!, &
          !                                            ppatch%patch_oce%bathymetry_c(jc,jb)
        END DO
      END DO

      !init normal velocity
      DO jb = i_startblk_e, i_endblk_e    
        CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,&
                    & i_startidx_e, i_endidx_e, rl_start, rl_end_e)

        DO je = i_startidx_e, i_endidx_e
          z_lat = ppatch%edges%center(je,jb)%lat
          z_lon = ppatch%edges%center(je,jb)%lon

          p_os%p_prog(nold(1))%vn(je,:,jb) = test_usbr_u(z_lon, z_lat,0.0_wp)* &
            &                                ppatch%edges%primal_normal(je,jb)%v1&
            &                                + test_usbr_v(z_lon, z_lat,0.0_wp)* &
            &                                ppatch%edges%primal_normal(je,jb)%v2
          ! write(*,*)'vn', je,jb,p_os%p_prog(nold(1))%vn(je,1,jb),z_lon, z_lat 

        END DO
      END DO
    CASE (26)
      CALL message(TRIM(routine), 'Shallow-Water-Testcase (26)')
      CALL message(TRIM(routine), ' - here: h and bathy of Williamson Test 5')

      ppatch%patch_oce%lsm_oce_c(:,:,:) = sea
      ppatch%patch_oce%lsm_oce_e(:,:,:) = sea
      !init height 
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                    & i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c
          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon

         p_os%p_prog(nold(1))%h(jc,jb) = test5_h( z_lon, z_lat, 0.0_wp)
         ! #slo# - bathymetry not taken into account due to constant bottom cell depth
         !       - use ppatch%patch_oce%zlev_i(dolic(jc,jb)) instead

        ppatch%patch_oce%bathymetry_c(jc,jb) = test5_oro( z_lon, z_lat, 0.0_wp )
        ! write(*,*)'h orig, bathy_c:', z_lon, z_lat,p_os%p_prog(nold(1))%h(jc,jb)!, &
        !                               ppatch%patch_oce%bathymetry_c(jc,jb)
        END DO
      END DO

      !init normal velocity
      DO jb = i_startblk_e, i_endblk_e    
        CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,&
                    & i_startidx_e, i_endidx_e, rl_start, rl_end_e)

        DO je = i_startidx_e, i_endidx_e
          z_lat = ppatch%edges%center(je,jb)%lat
          z_lon = ppatch%edges%center(je,jb)%lon

          p_os%p_prog(nold(1))%vn(je,1,jb) = &
          &   test5_u(z_lon, z_lat,0.0_wp)*ppatch%edges%primal_normal(je,jb)%v1  &
          & + test5_v(z_lon, z_lat,0.0_wp)*ppatch%edges%primal_normal(je,jb)%v2
          ! write(*,*)'vn', je,jb,p_os%p_prog(nold(1))%vn(je,1,jb),z_lon, z_lat 
        END DO
      END DO

    CASE(27)
      p_os%p_prog(nold(1))%h  = 0.0_wp
      p_os%p_prog(nold(1))%vn = 0.0_wp
      p_os%p_prog(nnew(1))%vn = 0.0_wp
      p_os%p_prog(nold(1))%tracer(:,1,:,1) = 0.0_wp
      p_os%p_prog(nnew(1))%tracer(:,1,:,1) = 0.0_wp
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c,&
                    & i_startidx_c, i_endidx_c, rl_start, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
         ppatch%patch_oce%bathymetry_c(jc,jb) = -200._wp
         ppatch%patch_oce%dolic_c(jc,jb)      = 1 

          !latitude given in radians
           z_lat = ppatch%cells%center(jc,jb)%lat
           z_lat_deg = z_lat*rad2deg
           !Impose emperature profile. Profile
           !depends on latitude only
           IF(abs(z_lat_deg-basin_center_lat)>=0.0_wp*basin_height_deg)THEN 
             p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = 5.0_wp 
             p_os%p_prog(nnew(1))%tracer(jc,1,jb,1) = 5.0_wp
           ELSEIF(abs(z_lat_deg-basin_center_lat)<0.0_wp*basin_height_deg)THEN 
             p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = 10.0_wp
             p_os%p_prog(nnew(1))%tracer(jc,1,jb,1) = 10.0_wp
           ENDIF
!           !write(90,*)'lat-degrees', jc,jb,z_lat, z_lat_deg, p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
        END DO
      END DO

     CASE(28)
      !init normal velocity
      DO jb = i_startblk_e, i_endblk_e    
        CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,&
                    & i_startidx_e, i_endidx_e, rl_start, rl_end_e)

        DO je = i_startidx_e, i_endidx_e
          z_lat = ppatch%edges%center(je,jb)%lat
          z_lon = ppatch%edges%center(je,jb)%lon
          IF(ppatch%patch_oce%lsm_oce_e(je,1,jb)<=sea_boundary)THEN
            p_os%p_prog(nold(1))%vn(je,1,jb) = &
            &   (test5_u(z_lon, z_lat,0.0_wp)*ppatch%edges%primal_normal(je,jb)%v1  &
            & + test5_v(z_lon, z_lat,0.0_wp)*ppatch%edges%primal_normal(je,jb)%v2)/40.0_wp
            ! write(*,*)'vn', je,jb,p_os%p_prog(nold(1))%vn(je,1,jb),z_lon, z_lat 
            p_os%p_prog(nnew(1))%vn(je,1,jb) = p_os%p_prog(nold(1))%vn(je,1,jb)
            p_os%p_diag%h_e(je,jb) = 1.0_wp
          ENDIF
        END DO
      END DO
      z_perlat = basin_center_lat! + 0.1_wp*basin_height_deg!             !45.5_wp
      z_perlon =  0.0_wp!0.1_wp*basin_width_deg                                 !4.5_wp
      z_permax  = 20.0_wp!20.1_wp
      z_perwid  =  10.0_wp!5.0_wp!1.5_wp
      
      DO jb = i_startblk_c, i_endblk_c    
        CALL get_indices_c(ppatch, jb, i_startblk_c, i_endblk_c, &
          &                i_startidx_c, i_endidx_c, rl_start, rl_end_c)

        DO jc = i_startidx_c, i_endidx_c
          z_lat = ppatch%cells%center(jc,jb)%lat
          z_lon = ppatch%cells%center(jc,jb)%lon

          IF(ppatch%patch_oce%lsm_oce_c(jc,1,jb)<=sea_boundary)THEN
            p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = 0.0_wp
            p_os%p_prog(nnew(1))%tracer(jc,1,jb,1) = 0.0_wp

            p_os%p_prog(nold(1))%h(jc,jb) = 1.0_wp!test5_h( z_lon, z_lat, 0.0_wp)

            z_dst=sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)
              !Local hot perturbation
            IF(z_dst<=5.0_wp*deg2rad)THEN
            p_os%p_prog(nold(1))%tracer(jc,1,jb,1) =          &
            !& 20.0_wp &!p_os%p_prog(nold(1))%tracer(jc,1,jb,1)          &
            &   z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) 

write(*,*)'init temp',p_os%p_prog(nold(1))%tracer(jc,1,jb,1), z_dst,&
&z_dst/(z_perwid*deg2rad),&
&z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2)
            ENDIF
            p_os%p_prog(nnew(1))%tracer(jc,1,jb,1)= p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
            p_os%p_prog(nnew(1))%h(jc,jb)         = p_os%p_prog(nold(1))%h(jc,jb)
          ENDIF
        END DO
      END DO
IF(idisc_scheme==1)THEN
CALL calc_scalar_product_for_veloc( ppatch,                &
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_diag%h_e,        &
                                    & p_os%p_diag)
ENDIF
! CALL rbf_vec_interpol_edge( p_os%p_prog(nold(1))%vn,&
!                           & ppatch,                &
!                           & p_int,                  &
!                           & p_os%p_diag%vt,         &
!                           & opt_slev=1, opt_elev=n_zlev)
! CALL rbf_vec_interpol_cell( p_os%p_prog(nold(1))%vn,&
!                           & ppatch,&
!                           & p_int,&
!                           & p_os%p_diag%u,  &
!                           & p_os%p_diag%v, &
!                           & opt_slev=1, opt_elev=n_zlev)
    CASE DEFAULT
     CALL finish(TRIM(routine), 'CHOSEN INITIALIZATION NOT SUPPORTED - TERMINATE')
  END SELECT
ENDIF

END SUBROUTINE init_ho_testcases


!------------Below are functions from to implement tests from Williamson shallow-water tests
!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !*IROUTINE:  test0_h
!  
! !F*UNCTION INTERFACE: 
  FUNCTION test0_h( p_lon, p_lat, p_t) RESULT( p_hh)
!
! !DESCRIPTION:
! Initial datum for height, test case 0 (conical mountain). \\
! Not included in Williamson et al. (1992)
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2006-12)
!  
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: h0=2000._wp  ! basic height level
    REAL(wp), PARAMETER  :: h1=1000._wp  ! max height of conical mountain
 
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_hh      ! geopotential height

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_r1      ! distance point to center of con. mount.
    REAL(wp)             :: z_r2      ! radius of conical mountain
    REAL(wp)             :: z_lon2    ! longitude of center of conical mount.
    REAL(wp)             :: z_lat2    ! latitude of center of conical mount.
    REAL(wp)             :: z_dlon    ! longitudinal distance 
    REAL(wp)             :: z_dlat    ! latitudinal distance 

!EOP  
!-----------------------------------------------------------------------  
!BOC

! center and radius of conical mountain

    z_lon2  = 1.5_wp * pi 
    z_lat2  = pi / 6._wp
    z_r2    = pi / 9._wp

! distance of point to center (not great arc distance!)

    z_dlon = p_lon - z_lon2
    z_dlat = p_lat - z_lat2

    z_dlon = z_dlon * z_dlon 
    z_dlat = z_dlat * z_dlat 

    z_r1 = z_dlon + z_dlat
    z_r1 = MIN( z_r2 * z_r2, z_r1)
    z_r1 = SQRT(z_r1)

! geopotential height

    IF( z_r1 < z_r2) THEN              ! point within radius

       p_hh = 1._wp - z_r1 / z_r2  
       p_hh = h0 + h1 * p_hh  

    ELSE                               ! point outside of radius

       p_hh = h0

    ENDIF

  END FUNCTION test0_h

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test2_h
!  
! !FUNCTION INTERFACE: 
  FUNCTION test2_h( p_lon, p_lat, p_t) RESULT(p_hh)
!
! !DESCRIPTION:
! Initial datum for height, test case 2 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2006-12)
!  
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: h0 = 2.94e4_wp * rgrav  ! maximum height
 
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_hh      ! height

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_fact1   ! 1st factor
    REAL(wp)             :: z_fact2   ! 2nd factor

!EOP  
!-----------------------------------------------------------------------  
!BOC

! 1st factor

    z_fact1 = re * omega 
    z_fact1 = z_fact1 + 0.5_wp * u0 
    z_fact1 = z_fact1 * u0 * rgrav

! 2nd factor

    z_fact2 = SIN(p_lat) * COS(aleph)
    z_fact2 = z_fact2 - COS(p_lon) * COS(p_lat) * SIN(aleph)
    z_fact2 = z_fact2 * z_fact2
    
! height

    p_hh = h0 - z_fact1 * z_fact2

  END FUNCTION test2_h

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test2_u
!  
! !FUNCTION INTERFACE: 
  FUNCTION test2_u( p_lon, p_lat, p_t) RESULT( p_uu)
!
! !DESCRIPTION:
! Initial datum for zonal velocity u, test case 2 of Williamson et al.(1992).
! Revised to programming guide by Th.Heinze, DWD, (2006-12)
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
!
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_uu      ! zonal velocity

!EOP  
!-----------------------------------------------------------------------  
!BOC

    p_uu = COS(p_lat) * COS(aleph)
    p_uu = p_uu + COS(p_lon) * SIN(p_lat) * SIN(aleph)
    p_uu = u0 * p_uu

  END FUNCTION test2_u

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test2_v
!  
! !FUNCTION INTERFACE: 
  FUNCTION test2_v( p_lon, p_lat, p_t) RESULT(p_vv)
!
! !DESCRIPTION:
! Initial datum for meridional velocity v, test case 2 of Williamson 
! et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2006-12)
!
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_vv      ! meridional velocity

!EOP  
!-----------------------------------------------------------------------  
!BOC

    p_vv = SIN(p_lon) * SIN(aleph)
    p_vv = -1._wp * u0 * p_vv

  END FUNCTION test2_v

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test2_vort
!  
! !FUNCTION INTERFACE: 
  FUNCTION test2_vort( p_lon, p_lat, p_t) RESULT(p_vort)
!
! !DESCRIPTION:
! Initial datum for relative vorticity, test case 2 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2006-12)
!
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_vort    ! relative vorticity

!EOP  
!-----------------------------------------------------------------------  
!BOC

    p_vort = SIN(p_lat)* COS(aleph)
    p_vort = p_vort - COS(p_lon) * COS(p_lat) * SIN(aleph)
    p_vort = 2._wp * u0 * rre * p_vort

  END FUNCTION test2_vort

!EOC  
!-------------------------------------------------------------------------  
!BOP
!  
! !IROUTINE:  test5_h
!  
! !FUNCTION INTERFACE: 
  FUNCTION test5_h( p_lon, p_lat, p_t) RESULT(p_hh)
!
! !DESCRIPTION:
! Initial datum for height, test case 5 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2007-01)
!  
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: h0    = 5960._wp  ! maximum height
    REAL(wp), PARAMETER  :: uzero = 20._wp    ! maximum velocity
 
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_hh      ! height

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_fact1   ! 1st factor
    REAL(wp)             :: z_fact2   ! 2nd factor

!EOP  
!-----------------------------------------------------------------------  
!BOC

! 1st factor

    z_fact1 = re * omega 
    z_fact1 = z_fact1 + 0.5_wp * uzero 
    z_fact1 = z_fact1 * uzero * rgrav

! 2nd factor

    z_fact2 = SIN(p_lat) 
    z_fact2 = z_fact2 * z_fact2
    
! height

    p_hh = h0 - z_fact1 * z_fact2

  END FUNCTION test5_h

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test5_u
!  
! !FUNCTION INTERFACE: 
  FUNCTION test5_u( p_lon, p_lat, p_t) RESULT( p_uu)
!
! !DESCRIPTION:
! Initial datum for zonal velocity u, test case 5 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2007-02)
!
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: uzero = 20._wp    ! maximum velocity

! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_uu      ! zonal velocity

!EOP  
!-----------------------------------------------------------------------  
!BOC

    p_uu = COS(p_lat) * COS(aleph)
    p_uu = p_uu + COS(p_lon) * SIN(p_lat) * SIN(aleph)
    p_uu = uzero * p_uu

  END FUNCTION test5_u

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test5_v
!  
! !FUNCTION INTERFACE: 
  FUNCTION test5_v( p_lon, p_lat, p_t) RESULT(p_vv)
!
! !DESCRIPTION:
! Initial datum for meridional velocity v, test case 5 of Williamson 
! et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2007-02)
!
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: uzero = 20._wp    ! maximum velocity

! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_vv      ! meridional velocity

!EOP  
!-----------------------------------------------------------------------  
!BOC

    p_vv = SIN(p_lon) * SIN(aleph)
    p_vv = -1._wp * uzero * p_vv

  END FUNCTION test5_v

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test5_oro
!  
! !FUNCTION INTERFACE:   
  FUNCTION test5_oro(p_lon, p_lat, p_t) RESULT(p_or)
!
! !DESCRIPTION:
! Initial datum for orography, test case 5 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Revised to programming guide by Th.Heinze, DWD, (2007-02)
!
!  
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: h_s0  = 2000._wp  ! maximum height of mountain
 
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_or      ! orography

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_lon_mc  ! Mountain center, longitude ...
    REAL(wp)             :: z_lat_mc  !          ... and latitude
    REAL(wp)             :: z_rad_mt  ! radius of mountain
    REAL(wp)             :: z_dist_mc ! distance from mountain center
    REAL(wp)             :: z_diff    ! difference of coordinates
    REAL(wp)             :: z_min_dist_sq ! min of square of distances

!EOP  
!-----------------------------------------------------------------------  
!BOC

! center and radius of mountain

    z_lon_mc = -pi_2
    z_lat_mc = pi / 6._wp    
    z_rad_mt = pi / 9._wp

! square of distance (in geographical coordinate sense) of point 
! from mountain center
 
    z_diff = p_lon - z_lon_mc
    z_diff = z_diff * z_diff
    z_dist_mc = z_diff

    z_diff = p_lat - z_lat_mc
    z_diff = z_diff * z_diff
    z_dist_mc = z_dist_mc + z_diff

! if point inside mountain range take its distance, else take mountain radius

    z_diff = z_rad_mt * z_rad_mt   
    z_min_dist_sq = MIN ( z_diff, z_dist_mc)
    z_dist_mc = SQRT( z_min_dist_sq)

! conical shape of mountain, depending on distance from mountain center

    p_or = z_dist_mc / z_rad_mt
    p_or = 1._wp - p_or
    p_or = h_s0 * p_or

  END FUNCTION test5_oro

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test6_h
!  
! !FUNCTION INTERFACE: 
  FUNCTION test6_h(p_lon, p_lat, p_t) RESULT(p_hh)
!
! !DESCRIPTION:
!
! Initial datum for geopotential h, test case 6 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Modified by Th.Heinze, DWD, (2006-11-02)

! !DEFINED PARAMETERS:  
!    REAL (wp), PARAMETER  :: h0 = 8000._wp, re_omg_kk = 50._wp
    REAL (wp), PARAMETER  :: h0 = 8000._wp, omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_hh

! !LOCAL VARIABLES:  
   ! REAL(wp)              :: z_omg, z_phia, z_phib, z_phic, z_r_omega 
    REAL(wp)              :: z_phia, z_phib, z_phic, z_r_omega , z_re_omg_kk 
    REAL(wp)              :: z_cosfi, z_cosfi2, z_cosfir, z_cosfir2, z_cosfir2m2
    REAL(wp)              :: z_cosdl, z_cosd2l, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

    !INTEGER               :: i_r1, i_r1r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r1, z_r1r2, z_r2

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_r_omega  = re * omega
    !z_omg     = re_omg_kk / re
    z_re_omg_kk= re * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
!    i_r1      = r + 1
!    i_r2      = r + 2
!    i_r1r1    = i_r1 * i_r1
!    i_r1r2    = i_r1 * i_r2
!    z_rr1r2   = 1._wp / i_r1r2 

    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r1    = z_r1 * z_r1
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    z_dlon    = omg_kk * z_r * (3._wp+z_r) - 2.0_wp * omega
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)  
    z_cosd2l  = COS(2._wp * z_dlon)  

    z_cosfi   = COS(p_lat)
    z_cosfi2  = z_cosfi  * z_cosfi    ! cos^2(lat)

    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfir2m2 = z_cosfir 
    z_cosfir2m2 = z_cosfir2m2 * z_cosfir2m2   ! cos^{2*r1-2}(lat)

    z_cosfir  = z_cosfir * z_cosfi    ! cos^{r1}(lat)
    z_cosfir2 = z_cosfir * z_cosfir   ! cos^{2*r1}(lat)

    z_val  = -.25_wp + z_r
    z_val  = 2._wp * z_val * z_val - 2.125_wp   ! 2r^2 - r -2

    z_phia = z_val * z_cosfi2
    
    z_val  = 2._wp * REAL(r,wp) * z_r

    z_phia = z_phia - z_val 
    z_val  = z_cosfi2 * z_cosfi2 * z_r1
    z_phia = z_phia + z_val 

    z_phia = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2m2 * z_phia 
    z_val  = .5_wp * z_re_omg_kk * (2._wp * z_r_omega + z_re_omg_kk) * z_cosfi2
    z_phia = z_val + z_phia 
    
    z_phib = -1._wp * z_cosfi2 * z_r1r1 + z_r1r1 + 1._wp 
    z_phib = z_re_omg_kk * (z_r_omega + z_re_omg_kk) * z_cosfir * z_phib
    z_phib = 2._wp * z_rr1r2 * z_phib
 
    z_phic = z_r1 * z_cosfi2 - 1._wp * z_r2
    z_phic = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2 * z_phic

    p_hh   = (z_phia + z_phib * z_cosdl + z_phic * z_cosd2l) * rgrav
    p_hh   = h0 + p_hh

  END FUNCTION test6_h

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test6_u
!  
! !FUNCTION INTERFACE: 
  FUNCTION test6_u( p_lon, p_lat, p_t) RESULT(p_uu)
!
! !DESCRIPTION:
!
! Initial datum for zonal velocity u, test case 6 of Williamson et al.(1992) .
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).    
! Modified by Th.Heinze, DWD, (2006-11-02)

! !DEFINED PARAMETERS:  
   ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER  ::  omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_uu

! !LOCAL VARIABLES:  
    !REAL(wp)              :: z_omg, z_r_omega 
    REAL(wp)              :: z_r_omega, z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_cosfi2, z_sinfi, z_sinfi2
    REAL(wp)              :: z_cosfir, z_cosfirm1
    REAL(wp)              :: z_cosdl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

!    INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2    !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_r_omega  = re * omega
    !z_omg     = re_omg_kk / re
    z_re_omg_kk= re * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    

    
!    i_r1      = r + 1
!    i_r2      = r + 2
!    i_r1r2    = i_r1 * i_r2
!    z_rr1r2   = 1._wp / i_r1r2 

    z_r       = REAL(r, wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * omega
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)  

    z_sinfi   = SIN(p_lat)
    z_sinfi2  = z_sinfi * z_sinfi

    z_cosfi   = COS(p_lat)
    z_cosfi2  = z_cosfi * z_cosfi

    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfirm1 = z_cosfir 

    z_val      = z_r * z_sinfi2 - z_cosfi2
    z_val      = z_cosfirm1 * z_val * z_cosdl
    z_val      = z_cosfi + z_val
    p_uu       = z_re_omg_kk * z_val 

  END FUNCTION test6_u

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test6_v
!  
! !FUNCTION INTERFACE: 
  FUNCTION test6_v( p_lon, p_lat, p_t) RESULT(p_vv)
!
! !DESCRIPTION:
!
! Initial datum for meridional velocity v, test case 6 of Williamson 
! et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).    
! Modified by Th.Heinze, DWD, (2006-11-02)

! !DEFINED PARAMETERS:  
   ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER  ::  omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_vv

! !LOCAL VARIABLES:  
    !REAL(wp)              :: z_omg, z_r_omega   !pripodas, we use omg_kk and not re_omg_kk
    REAL(wp)              :: z_r_omega, z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_sinfi
    REAL(wp)              :: z_cosfir, z_cosfirm1
    REAL(wp)              :: z_sindl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

!    INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2    !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_r_omega = re * omega
    !z_omg     = re_omg_kk / re
    z_re_omg_kk= re * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
!    i_r1      = r + 1
!    i_r2      = r + 2
!    i_r1r2    = i_r1 * i_r2
!    z_rr1r2   = 1._wp / i_r1r2 

    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2 

    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * omega
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_sindl   = SIN(z_dlon)  

    z_sinfi   = SIN(p_lat)

    z_cosfi   = COS(p_lat)

    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfirm1 = z_cosfir 

    z_val      = z_cosfirm1 * z_sinfi * z_sindl
    p_vv       = -1._wp * z_re_omg_kk * z_r * z_val 

  END FUNCTION test6_v

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test6_vort
!  
! !FUNCTION INTERFACE: 
   FUNCTION test6_vort( p_lon, p_lat, p_t) RESULT(p_vt)
!
! !DESCRIPTION:
!
! Initial datum for relative vorticity, test case 6 of Williamson et al.(1992).
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).    
! Modified by Th.Heinze, DWD, (2006-11-02):
! - corrected vorticity

! !DEFINED PARAMETERS:  
   ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER  ::  omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_vt

! !LOCAL VARIABLES:  
    !REAL(wp)              :: z_omg, z_r_omega   !pripodas, we use omg_kk and not re_omg_kk
    REAL(wp)              :: z_r_omega, z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_sinfi
    REAL(wp)              :: z_cosfir
    REAL(wp)              :: z_cosdl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

    !INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2   !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_r_omega = re * omega
    !z_omg     = re_omg_kk / re
    z_re_omg_kk= re * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
   ! i_r1      = r + 1
   ! i_r2      = r + 2
   ! i_r1r2    = i_r1 * i_r2
   ! z_rr1r2   = 1._wp / i_r1r2 

    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * omega
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)  

    z_sinfi   = SIN(p_lat)

    z_cosfi   = COS(p_lat)

    z_cosfir  = z_cosfi
    DO j= 2, r
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO

    z_val     = z_cosfir * z_r1 * z_r2 * z_cosdl
    z_val     = 2._wp - z_val
    p_vt      = omg_kk * z_sinfi * z_val

  END FUNCTION test6_vort

!EOC  
!-------------------------------------------------------------------------  
!BOP
!  
! !IROUTINE:  test_usbr_h
!  
! !FUNCTION INTERFACE: 
  FUNCTION test_usbr_h( p_lon, p_lat, p_t) RESULT(p_hh)
!
! !DESCRIPTION:
! Initial datum for height h, test case unsteady solid body 
! rotation of L\"auter et al.(2007).

! !REVISION HISTORY:  
! Developed by Th.Heinze, DWD, (2007-03)
!  
! !DEFINED PARAMETERS:  
    REAL(wp), PARAMETER  :: d0    = 133681.0_wp  ! additive constant
 
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_hh      ! height

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_phi_t_k ! 1st summand
    REAL(wp)             :: z_summand ! 2nd summand
    REAL(wp)             :: z_fact    ! factor
    REAL(wp)             :: z_angle1  ! 1st angle
    REAL(wp)             :: z_angle2  ! 2nd angle

!EOP  
!-----------------------------------------------------------------------  
!BOC

! relevant angles

    z_angle1 = .25_wp * pi
    z_angle2 = p_lon + omega * p_t

! 1st summand: \phi_t(\vec c) \cdot \vec k

    z_phi_t_k = SIN(p_lat) * COS(z_angle1)
    z_phi_t_k = z_phi_t_k - COS(z_angle2) * COS(p_lat) * SIN(z_angle1)
    z_phi_t_k = u0 * z_phi_t_k

! 2nd summand: r_e \Omega \sin \varphi

    z_summand = re * omega * SIN(p_lat)

! one factor

    z_fact    = .5_wp *  z_phi_t_k + z_summand
    
! height

    p_hh      = d0 - z_phi_t_k *  z_fact
    p_hh      = p_hh * rgrav
!write(*,*)'param:', u0, pi, rgrav,re, omega
!stop
  END FUNCTION test_usbr_h

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test_usbr_u
!  
! !FUNCTION INTERFACE: 
  FUNCTION test_usbr_u( p_lon, p_lat, p_t) RESULT( p_uu)
!
! !DESCRIPTION:
! Initial datum for zonal velocity u, test case unsteady solid body 
! rotation of L\"auter et al.(2007).

! !REVISION HISTORY:  
! Developed by Th.Heinze, DWD, (2007-03)
!
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_uu      ! zonal velocity

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_angle1  ! 1st angle
    REAL(wp)             :: z_angle2  ! 2nd angle

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_angle1 = .25_wp * pi
    z_angle2 = p_lon + omega * p_t
    p_uu = COS(p_lat) * COS(z_angle1)
    p_uu = p_uu + COS(z_angle2) * SIN(p_lat) * SIN(z_angle1)
    p_uu = u0 * p_uu

  END FUNCTION test_usbr_u

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test_usbr_v
!  
! !FUNCTION INTERFACE: 
  FUNCTION test_usbr_v( p_lon, p_lat, p_t) RESULT(p_vv)
!
! !DESCRIPTION:
! Initial datum for meridional velocity v, test case unsteady solid body 
! rotation of L\"auter et al.(2007).

! !REVISION HISTORY:  
! Developed by Th.Heinze, DWD, (2007-03)
!
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_vv      ! meridional velocity

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_angle   ! angle

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_angle = p_lon + omega * p_t
    p_vv = SIN(z_angle) 
    z_angle = .25_wp * pi
    p_vv = p_vv * SIN(z_angle) 
    p_vv = -1._wp * u0 * p_vv

  END FUNCTION test_usbr_v

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  test_usbr_oro
!  
! !FUNCTION INTERFACE:   
  FUNCTION test_usbr_oro(p_lon, p_lat, p_t) RESULT(p_or)
!
! !DESCRIPTION:
! Initial datum for orography, test case unsteady solid body rotation
! of L\"auter et al.(2007).
!
! !REVISION HISTORY:  
! Developed by Th.Heinze, DWD, (2007-03)
! 
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time

! !RETURN VALUE:  
    REAL(wp)             :: p_or      ! orography

! !LOCAL VARIABLES:  
    REAL(wp)             :: z_fact    ! factor

!EOP  
!-----------------------------------------------------------------------  
!BOC

! calculate factor

    z_fact = re * omega * SIN(p_lat)
    z_fact = z_fact * z_fact

! height of orography

    p_or = .5_wp * z_fact * rgrav

  END FUNCTION test_usbr_oro

!EOC  

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  rotate
!  
! !SUBROUTINE INTERFACE: 
  SUBROUTINE rotate(p_lon, p_lat, p_alpha, p_rotlon, p_rotlat)
!
! !DESCRIPTION:
! This subroutine computes the rotated coordinates p\_rotlon, p\_rotlat
! for a roatation by angle p\_alpha, given the coordinates p\_lon and p\_lat.
!
! !REVISION HISTORY:  
! Developed originally by R.Jakob for NCAR shallow water model.  
! Adapted to ICON code by L.Bonaventura (2002-5).
! Adapted to ICON programming guide by Th.Heinze, DWD, (2006-12-12)
!
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in)  :: p_lon     ! ORIGINAL LONGITUDE
    REAL(wp), INTENT(in)  :: p_lat     ! ORIGINAL LATITUDE
    REAL(wp), INTENT(in)  :: p_alpha   ! ROTATION ANGLE
   
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(out) :: p_rotlon  ! ROTATED LONGITUDE
    REAL(wp), INTENT(out) :: p_rotlat  ! ROTATED LATITUDE
    
! !LOCAL VARIABLES:  
    REAL(wp)              :: z_test    ! checking value
    
!EOP  
!-----------------------------------------------------------------------  
!BOC

    IF (p_alpha == 0.0_wp) THEN       !        NO ROTATION

      p_rotlon = p_lon
      p_rotlat = p_lat

    ELSE                              !        ROTATION BY ANGLE p_alpha

!     ROTATED LATITUDE

      z_test = SIN(p_lat)*COS(p_alpha)- COS(p_lat)*COS(p_lon)*SIN(p_alpha)

      IF (z_test > 1.0_wp) THEN
        p_rotlat = pi_2
      ELSEIF (z_test < -1.0_wp) THEN
        p_rotlat = -1.0_wp * pi_2
      ELSE
        p_rotlat = ASIN(z_test)
      ENDIF

!     ROTATED LONGITUDE
       
      z_test = COS(p_rotlat)
       
      IF (z_test == 0.0_wp) THEN
        p_rotlon = 0.0_wp
      ELSE
        z_test = SIN(p_lon)*COS(p_lat)/z_test
        IF (z_test > 1.0_wp) THEN
          p_rotlon = pi_2
        ELSEIF (z_test < -1.0_wp) THEN
          p_rotlon = -1.0_wp * pi_2
        ELSE
          p_rotlon = ASIN(z_test)
        ENDIF
      ENDIF

!        ADJUST FOR CORRECT BRANCH OF INVERSE SINE
       
      z_test = COS(p_alpha)*COS(p_lon)*COS(p_lat) + SIN(p_alpha)*SIN(p_lat)
 
      IF (z_test < 0.0_wp) THEN
        p_rotlon = pi - p_rotlon
      ENDIF
 
    ENDIF
    
  END SUBROUTINE rotate

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  geostr_balance
!
! !FUNCTION INTERFACE: 
  FUNCTION geostr_balance( p_lat, func)  RESULT(p_hh)
!
! !DESCRIPTION:
! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
! to compute geostrophically balanced initial state used
! in test 3.
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
! Modified by Th.Heinze, DWD, (2006-11-22): 
! - introduced INTERFACE uu (got an error message with g95 compiler,
!   scanned the code, this seems to be the correct way, but might be wrong)
! Modified by Th.Heinze, DWD, (2006-12-12): 
! - renamed it to geostr_balance
!
! !REMARKS:
! was htmp2 in previous code

! !INTERFACE:
    INTERFACE                        ! selected function

      FUNCTION func(p_t) RESULT(p_vv)  

        USE mo_kind, ONLY: wp
	       
        REAL(wp), INTENT(in) :: p_t
        REAL(wp)             :: p_vv

      END FUNCTION func
       
    END INTERFACE
   
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lat           ! rotated latitude

! !RETURN VALUE:  
    REAL(wp)             :: p_hh            ! balanced height

! !LOCAL VARIABLES:  
    INTEGER              :: j               ! loop index

    REAL(wp)             :: z_a             ! left bound
    REAL(wp)             :: z_b             ! right bound
    REAL(wp)             :: z_lat           ! latitude in loop
    REAL(wp)             :: z_step          ! step
    REAL(wp)             :: z_val, z_val2   ! intermediate values


!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_a = -1._wp * pi_2
    z_b = p_lat
    
    z_step = 0.02_wp * ( z_b - z_a)

    p_hh = 0._wp

    z_lat = z_a - 0.5_wp * z_step

    DO j = 1, 50
       z_lat = z_lat + z_step
       
       z_val = func(z_lat)

       z_val2 = 2._wp * omega * SIN(z_lat)
       z_val2 = z_val2 + z_val * TAN(z_lat)* rre
       z_val2 = z_val * z_val2

       p_hh = p_hh + z_val2 * z_step

    ENDDO

  END FUNCTION geostr_balance

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  zero
!  
! !FUNCTION INTERFACE: 
  
  FUNCTION zero(lon,lat,t) RESULT(uu)
!
! !DESCRIPTION:
!
! Dummy constant function zero, used in the
! initialization of some test cases.
!
! !REVISION HISTORY:  
! Developed  by L.Bonaventura  (2002-5).
    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lon,lat,t
    REAL(wp) :: uu

!EOP  
!-----------------------------------------------------------------------  
!BOC

    uu=0._wp
  END FUNCTION zero



 FUNCTION test11_h(lon,lat,t) RESULT(hh)

    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lon,lat,t
    REAL(wp) :: hh, hdach, alpha, beta, phi2
    REAL(wp)             :: z_rotlon  ! rotated longitude
    REAL(wp)             :: z_rotlat  ! rotated latitude

     ! rotate
     CALL rotate( lon, lat, aleph, z_rotlon, z_rotlat)
     ! calculate height
     hh = geostr_balance11( z_rotlat, test11_u2)


     hdach = 120._wp
     alpha = 1._wp/3._wp
     beta  = 1._wp/15._wp
     phi2  = pi/4._wp
     hh    = hh + hdach*cos(lat)*exp(-((lon)/alpha)**2)*exp(-((phi2-lat)/beta)**2)

  END FUNCTION test11_h


 FUNCTION test11_u(lon,lat,t) RESULT(uu)

    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lon,lat,t
    REAL(wp) ::  uu, d
     REAL(wp) ::  phi0, phi1, umax, en

    phi0=pi/7._wp
    phi1=pi/2._wp - phi0
    en=exp(-4._wp/(phi0-phi1)**2)
    umax=80._wp

    d=.1_wp

    if ((lat.gt.phi0).and.(lat.lt.phi1))then
         uu=umax/en*exp(1._wp/(lat-phi0)/(lat-phi1))
         if (uu.lt. 0.001_wp) then
               uu  = 0.0_wp
         end if
 !         print*, "assigning u values", uu
    else
         uu=0._wp
    endif

!    1451 !     For jet on southern hemisphere additionally:
!    1452 ! if ((lat.lt.-phi0).and.(lat.gt.-phi1)) then
!    1453 ! uu=+umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For volume tests
!    1454 ! ! uu=uu-umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For Galewsky tests
!    1455 ! endif

  END FUNCTION test11_u

  FUNCTION test11_u2(lat) RESULT(uu)

    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lat
    REAL(wp) ::  uu, d
    REAL(wp) ::  phi0, phi1, umax, en

    phi0=pi/7._wp
    phi1=pi/2._wp - phi0
    en=exp(-4._wp/(phi0-phi1)**2)
    umax=80._wp

    d=.1_wp

    if ((lat.gt.phi0).and.(lat.lt.phi1))then
         uu=umax/en*exp(1._wp/(lat-phi0)/(lat-phi1))
         if (uu.lt. 0.001_wp) then
               uu  = 0.0_wp
         end if
 !         print*, "assigning u values", uu
    else
         uu=0._wp
    endif

 !     For jet on southern hemisphere additionally:
 ! if ((lat.lt.-phi0).and.(lat.gt.-phi1)) then
 ! uu=+umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For volume tests
 ! ! uu=uu-umax/en*exp(1._wp/(lat+phi0)/(lat+phi1))!!! For Galewsky tests
 ! endif

  END FUNCTION test11_u2

  FUNCTION test11_v(lon,lat,t) RESULT(vv)

    IMPLICIT NONE
    REAL(wp) , INTENT(in):: lon,lat,t
    REAL(wp) ::  vv

    vv = 0.0_wp

  END FUNCTION test11_v


 FUNCTION geostr_balance11( phi, func)  RESULT(p_hh)

! !DESCRIPTION:
 ! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
 ! to compute geostrophically balanced initial state used
 ! in test 3.
 !
 ! !REVISION HISTORY:
 ! Developed  by L.Bonaventura  (2002-5).
 ! Modified by Th.Heinze, DWD, (2006-11-22):
 ! - introduced INTERFACE uu (got an error message with g95 compiler,
 !   scanned the code, this seems to be the correct way, but might be wrong)
 ! Modified by Th.Heinze, DWD, (2006-12-12):
 ! - renamed it to geostr_balance
 ! Modified by F. Rauser, MPI (2009,10) for testcase 11 galewsky
 !
 ! !REMARKS:
 ! was htmp2 in previous code

 ! !INTERFACE:
     INTERFACE                        ! selected function

       FUNCTION func(p_t) RESULT(p_vv)

         USE mo_kind, ONLY: wp

         REAL(wp), INTENT(in) :: p_t
         REAL(wp)             :: p_vv

       END FUNCTION func

     END INTERFACE

 ! !INPUT PARAMETERS:
     REAL(wp), INTENT(in) :: phi           ! rotated latitude
 ! !RETURN VALUE:
     REAL(wp)             :: p_hh            ! balanced height
 ! !LOCAL VARIABLES:
     INTEGER              :: j               ! loop index
     REAL(wp)             :: phi_a             ! left bound
     REAL(wp)             :: phi_b             ! right bound
     REAL(wp)             :: phidash           ! latitude in loop
     REAL(wp)             :: dphi          ! step
     REAL(wp)             :: u, temp   ! intermediate values
 !EOP
 !-----------------------------------------------------------------------
 !BOC

     phi_a = -0.5_wp * pi
     phi_b = phi

     dphi = 0.01_wp * ( phi_b - phi_a)

     p_hh = 0._wp

     phidash = phi_a - 0.5_wp * dphi

     DO j = 1, 100
        phidash = phidash + dphi

       u = func(phidash)

        temp = 2._wp * omega * SIN(phidash)
        temp = temp + ( u * TAN(phidash)* rre)
        temp = re *rgrav * u * temp

        p_hh = p_hh + temp * dphi

     ENDDO

     p_hh = 10000._wp - p_hh
 !     print*, "phh", INT(360*phi/pi), INT(p_hh)

   END FUNCTION geostr_balance11


END MODULE mo_oce_init
