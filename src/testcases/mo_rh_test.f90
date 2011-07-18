!>
!!  This module contains parameters, initialization subroutines and.
!!
!!  This module contains parameters, initialization subroutines and
!!  functions to be used in the 3D Rossby-Haurwitz wave test case of the
!!  hydrostatic dynamical core.
!!
!! @par Revision History
!!  Original version by Hui Wan, MPI-M (2006-09-10)
!!  Modification by Hui Wan, MPI-M (2007-07-26):
!!  - removed the *grid* structure according to recent changes by Almut.
!!  Modification by Almut Gassmann, MPI-M (2008-04-24)
!!  - minor changes for NCAR workshop
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
!!
MODULE mo_rh_test
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2005
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: re, omega, rd, grav
  USE mo_vertical_coord_table,ONLY: vct_a,vct_b
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data,            ONLY: t_external_data
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog
  USE mo_parallel_config,  ONLY: nproma

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: init_hydro_state_prog_rhtest


  REAL(wp), PARAMETER :: ps0   = 95500.0_wp        ! surface pressure (Pa)
  REAL(wp), PARAMETER :: lr0   = 0.0065_wp*rd/grav ! lapse rate
  REAL(wp), PARAMETER :: temp0 = 288._wp           ! horizontal-mean
                                                   ! temperature at surface (K)
  REAL(wp), PARAMETER :: rgeo0 = lr0/rd/temp0

  INTEGER  :: rh_wavenum   ! zonal wavenumber, used in namelist
  REAL(wp) :: A0  ! = 50._wp/re/m, amplitude of zonal flow
  REAL(wp) :: Am  ! = 50._wp/re/m, amplitude of wave
  REAL(wp) :: rh_init_shift_deg  ! shift of the initial state

  PUBLIC :: rh_wavenum, rh_init_shift_deg

 CONTAINS

!--------------------------------------------------------------------

!-------------------------------------------------------------------------
!

 !>
 !!               Initialization of prognostic state vector.
 !!
 !!
 !! @par Revision History
 !!  Original version  by Hui Wan, MPI-M (2006-09)
 !!  Modification by Hui Wan, MPI-M (2007-07-26):
 !!  - removed the *grid* structure according to recent changes by Almut.
 !! @par
 !!  arguments
 !!
 SUBROUTINE init_hydro_state_prog_rhtest(pt_patch, pt_prog, pt_ext_data)


  IMPLICIT NONE

  TYPE(t_patch), INTENT(INOUT)            :: pt_patch
  TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
  TYPE(t_external_data), INTENT(INOUT)    :: pt_ext_data  !< external data

  REAL(wp) :: tmp1,tmp2,tmp3, zgeo, zu,zv
  REAL(wp) :: lat,lon
  REAL(wp) :: zsinlat,zcoslat,zsinmlon,zcosmlon
  REAL(wp) :: zpk,zpk1,zpres
  REAL(wp) :: zshift

  INTEGER :: m, jb, je, jc, jk, nlen, nblks_c, nblks_e, npromz_e, npromz_c
  INTEGER :: nlev              !< number of full levels

!--------------------------------------------------------------------

  nlev   = pt_patch%nlev

  m   = rh_wavenum
  A0  = 50._wp/re/REAL(m,wp)
  Am  = 50._wp/re/REAL(m,wp)

  zshift = rh_init_shift_deg*pi/180._wp
!
!---initialize the prognostic variables
!
     nblks_c   = pt_patch%nblks_int_c
     npromz_c  = pt_patch%npromz_int_c
     nblks_e   = pt_patch%nblks_int_e
     npromz_e  = pt_patch%npromz_int_e

     ! topography
     pt_ext_data%atm%topography_c(:,:) = 0.0_wp

     ! surface pressure
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nlen,lat,lon,zsinlat,zcoslat,zcosmlon,tmp1,tmp2,tmp3,zgeo)
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        DO jc =1,nlen

           lat      = pt_patch%cells%center(jc,jb)%lat
           lon      = pt_patch%cells%center(jc,jb)%lon + zshift
           zsinlat  = SIN(lat)
           zcoslat  = COS(lat)
           zcosmlon = COS(REAL(m,wp)*lon)

           tmp1 =   0.5_wp*A0*( 2._wp*omega + A0 )*zcoslat*zcoslat       &
                & + 0.25_wp*Am*Am*zcoslat**(2*m)                         &
                & *( (REAL(m,wp)+1._wp)*zcoslat*zcoslat                  &
                &   +(2._wp*REAL(m,wp)*REAL(m,wp) - REAL(m,wp) -2._wp)   &
                &   - 2._wp*REAL(m,wp)*REAL(m,wp)/zcoslat/zcoslat )

           tmp2 =  2._wp*(omega+A0)*Am/(REAL(m,wp)+1._wp)/(REAL(m,wp)+2._wp) &
                & * zcoslat**m *( REAL(m,wp)*REAL(m,wp)                      &
                & + 2._wp*REAL(m,wp) + 2._wp                                 &
                & -(REAL(m,wp)+1._wp)*(REAL(m,wp)+1._wp)*zcoslat*zcoslat )

           tmp3 =  0.25_wp*Am*Am * zcoslat**(2*m)  &
                & *( (REAL(m,wp)+1._wp)*zcoslat*zcoslat -REAL(m,wp) -2._wp )

           zgeo = tmp1 + tmp2*zcosmlon + tmp3*( 2._wp*zcosmlon*zcosmlon-1._wp )

           pt_prog%pres_sfc(jc,jb) = &
                  ps0*(1._wp+rgeo0*re*re*zgeo)**(1._wp/lr0)
        ENDDO
     ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,je,nlen,lat,lon,zsinlat,zcoslat,zsinmlon,zcosmlon,zu,zv)
     ! normal wind
     DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF

        DO je =1,nlen
           lat = pt_patch%edges%center(je,jb)%lat
           lon = pt_patch%edges%center(je,jb)%lon + zshift

           zsinlat  = SIN(lat)
           zcoslat  = COS(lat)
           zsinmlon = SIN(REAL(m,wp)*lon)
           zcosmlon = COS(REAL(m,wp)*lon)

           zu =  re*A0*zcoslat                                          &
              & +re*Am*( REAL(m,wp)*zsinlat*zsinlat - zcoslat*zcoslat ) &
              &    *zcoslat**(m-1) *zcosmlon

           zv = -re*REAL(m,wp)*Am* zcoslat**(m-1) *zsinlat*zsinmlon

           pt_prog%vn(je,:,jb) =                                     &
             &          zu * pt_patch%edges%primal_normal(je,jb)%v1  &
             &        + zv * pt_patch%edges%primal_normal(je,jb)%v2

        ENDDO
     ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,jc,jk,nlen,zpk,zpk1,zpres)
     ! temperature
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! the uppermost layer
        DO jc = 1, nlen

           zpk  = vct_a(1) + vct_b(1)*pt_prog%pres_sfc(jc,jb)
           zpk1 = vct_a(2) + vct_b(2)*pt_prog%pres_sfc(jc,jb)
           zpres = 0.5_wp*(zpk1-zpk)

           pt_prog%temp(jc,1,jb) = temp0*(zpres/ps0)**lr0
        END DO

        ! the other vertical layers

        DO jk=2,nlev
           DO jc=1,nlen

              zpk  = vct_a(jk)   + vct_b(jk)  *pt_prog%pres_sfc(jc,jb)
              zpk1 = vct_a(jk+1) + vct_b(jk+1)*pt_prog%pres_sfc(jc,jb)
              zpres = EXP( (zpk1*LOG(zpk1)-zpk*LOG(zpk))/(zpk1-zpk) - 1._wp)

              pt_prog%temp(jc,jk,jb) = temp0*(zpres/ps0)**lr0
           ENDDO
        ENDDO
     ENDDO  ! block loop
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_hydro_state_prog_rhtest

END MODULE mo_rh_test
