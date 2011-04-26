!>
!! This module contains a subroutine used to initialize the
!! hydrostatic prognostic state for the local diabatic
!! forcing test case
!!
!! mo_ldf_init_prog_state is called from mo_hydro_testcases
!! when switching on local diabatic forcing
!!
!! @author Constantin Junk, MPI-M
!!
!! @par Revision History
!! Original implementation by Constantin Junk, MPI-M (2011-03-28)
!! Modification by Constantin Junk, MPI-M (2011-04-05)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_ldf_init

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi, pi_2
  USE mo_physical_constants,  ONLY: re, rgrav, omega, rd, tmelt
  USE mo_advection_nml,       ONLY: ctracer_list
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data,            ONLY: t_external_data
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_run_nml,             ONLY: ltransport, ntracer, nproma, iforcing, &
    &                               ildf_echam, iqv, iqt
  USE mo_vertical_coord_table,ONLY: ceta
  USE mo_ncar_testcases,      ONLY: tracer_q1_q2, tracer_q3,regrot, turnwi
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                               sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
    &                               spec_humi          !! Specific humidity
  USE mo_eta_coord_diag,      ONLY: half_level_pressure,full_level_pressure
  USE mo_exception,           ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ldf_init_prog_state

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! Variables for initialisation of hydrostatic prognostic state

  REAL(wp), PARAMETER :: eta0  = 0.252_wp !
  REAL(wp), PARAMETER :: etat  = 0.2_wp   ! tropopause
  REAL(wp), PARAMETER :: ps0   = 1.e5_wp  ! surface pressure (Pa)
  REAL(wp), PARAMETER :: u0    = 25._wp
  REAL(wp), PARAMETER :: temp0 = 288._wp  ! horizontal-mean temperature
                                          ! at surface (K)
  REAL(wp), PARAMETER :: gamma = 0.005_wp ! temperature elapse rate (K/m)
  REAL(wp), PARAMETER :: dtemp = 4.8e5_wp ! empirical temperature difference (K)


CONTAINS


 !>
 !! Initialization of prognostic state vector for the 
 !! Local Diabatic Forcing testcase. The initialisation
 !! of the prognostic variables is similar to the JWs
 !! testcase. However, to account for easterlies in the
 !! tropics, the zonal wind profile was changed and, thus,
 !! the surface geopotential and zonal temperature profile
 !! as well. The meridional wind is set to zero by design.
 !!
 !! @par Revision History
 !!  Original version for JWw by Hui Wan, MPI-M (2006-02)
 !!  Modified by Constantin Junk, MPI-M (2011-28-03)
 !!   - removed perturbation part and optional parameters
 !!   - added different zonal wind field, temperature
 !!     and surface geopotential initialisation 
 !! @par
 !!  arguments
 !!
 SUBROUTINE ldf_init_prog_state(pt_patch, pt_hydro_prog,      &
   &                            pt_hydro_diag, pt_ext_data,   &
   &                            p_rotate_axis_deg,            &
   &                            opt_lrh_linear_pres,          &
   &                            opt_rh_at_1000hpa             )

  IMPLICIT NONE

  TYPE(t_patch), INTENT(INOUT)            :: pt_patch
  TYPE(t_hydro_atm_prog), INTENT(INOUT)   :: pt_hydro_prog
  TYPE(t_hydro_atm_diag), INTENT(INOUT)   :: pt_hydro_diag
  TYPE(t_external_data), INTENT(INOUT)    :: pt_ext_data   !< external data

  REAL(wp),INTENT(IN) :: p_rotate_axis_deg

  LOGICAL, INTENT(IN),OPTIONAL :: opt_lrh_linear_pres
  REAL(wp),INTENT(IN),OPTIONAL :: opt_rh_at_1000hpa

  !local variables

  CHARACTER(LEN=1) :: ctracer

  REAL(wp) :: lon,lat,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5, rot_lon, rot_lat
  REAL(wp) :: zeta,zcos12z,zcos32z,zsinz,zcosysq,zsin2ysq,zsiny,zcosy,ztemp
  REAL(wp) :: zrotate_axis_rad, zu, zv

  REAL(wp) :: etav(pt_patch%nlev)        !
  REAL(wp) :: sin_etav(pt_patch%nlev)    !
  REAL(wp) :: cos_12_etav(pt_patch%nlev) !
  REAL(wp) :: cos_32_etav(pt_patch%nlev) !
  REAL(wp) :: temp_avg(pt_patch%nlev)    ! horizontally averaged temperature

  REAL(wp) :: zsqv, zrhf

  INTEGER :: nblks_c, nblks_e, nblks_v, npromz_e, npromz_c, npromz_v, &
             nlen, jt, jb, je, jc, jk, jv
  INTEGER :: nlev                        !< number of full levels

  LOGICAL  :: lrh_linear_pres
  REAL(wp) :: rh_at_1000hpa


!--------------------------------------------------------------------
! First check optional arguments

  IF (PRESENT(opt_lrh_linear_pres)) THEN
    lrh_linear_pres = opt_lrh_linear_pres
  ELSE
    lrh_linear_pres = .FALSE.
  ENDIF

  IF (lrh_linear_pres) THEN
     IF (PRESENT(opt_rh_at_1000hpa)) THEN
        rh_at_1000hpa = opt_rh_at_1000hpa
     ELSE
        rh_at_1000hpa = 0.75_wp
     ENDIF
  ENDIF

!-----------------------

  ! number of vertical levels
  nlev = pt_patch%nlev

  tmp0 = rd*gamma*rgrav

  zrotate_axis_rad = p_rotate_axis_deg * pi/180._wp
  DO jk = 1, nlev
     zeta            = ceta(jk)
     etav(jk)        = pi_2*(zeta-eta0)
     sin_etav(jk)    = SIN(etav(jk))
     cos_12_etav(jk) = COS(etav(jk))**0.5_wp
     cos_32_etav(jk) = COS(etav(jk))**1.5_wp
     temp_avg(jk)    = zeta**tmp0 * temp0
     IF ( zeta < etat ) temp_avg(jk) = temp_avg(jk) + (etat-zeta)**5*dtemp
  END DO

!
!---initialize the prognostic variables
!
     nblks_c   = pt_patch%nblks_int_c
     npromz_c  = pt_patch%npromz_int_c
     nblks_e   = pt_patch%nblks_int_e
     npromz_e  = pt_patch%npromz_int_e
     nblks_v   = pt_patch%nblks_int_v
     npromz_v  = pt_patch%npromz_int_v

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nlen,lon,lat,rot_lon,rot_lat,zsiny,zcosy,tmp1,tmp2,tmp3,tmp4,tmp5)
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! set the topography
        !

        DO jc = 1, nlen

           lon   = pt_patch%cells%center(jc,jb)%lon
           lat   = pt_patch%cells%center(jc,jb)%lat

           IF ( ABS(p_rotate_axis_deg) < 1.0e-8_wp ) THEN
              rot_lon = lon
              rot_lat = lat
           ELSE
              CALL regrot(lon,lat,rot_lon,rot_lat,0._wp,&
                          -pi_2+zrotate_axis_rad,1)
           ENDIF

           zsiny = SIN(rot_lat)
           zcosy = COS(rot_lat)

           tmp1  = u0*COS((1._wp-eta0)*pi_2)**1.5_wp
           tmp2  = tmp1*(-2._wp*zsiny**6*(1._wp/3._wp+zcosy**2) + 1.0_wp/6.3_wp)
           tmp3  = omega*re*(8._wp/5._wp*zcosy**3*(2._wp/3._wp+zsiny**2) - pi/4_wp)
           tmp4  = 8._wp*(1._wp-2._wp/3._wp*zsiny**2)*zsiny**4 - 88._wp/105._wp
           tmp5  = zcosy**3*(4._wp*zcosy - 8._wp/3._wp) +pi/2._wp -3.2_wp/1.5_wp


           pt_ext_data%atm%topography_c(jc,jb) = (tmp1*(tmp2+tmp3+tmp4)+tmp5)*rgrav
           ! Coriolis parameter
           pt_patch%cells%f_c(jc,jb) = 2.0_wp*omega*(SIN(lat)*COS(zrotate_axis_rad)&
                                     -COS(lon)*COS(lat)*SIN(zrotate_axis_rad))

        ENDDO
     ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,jv,nlen,lon,lat)

     DO jb = 1, nblks_v
        IF (jb /= nblks_v) THEN
           nlen = nproma
        ELSE
           nlen = npromz_v
        ENDIF
        DO jv = 1, nlen
           lon= pt_patch%verts%vertex(jv,jb)%lon
           lat= pt_patch%verts%vertex(jv,jb)%lat
           ! Coriolis parameter
           pt_patch%verts%f_v(jv,jb) = 2.0_wp*omega*(SIN(lat)*COS(zrotate_axis_rad)&
                                     -COS(lon)*COS(lat)*SIN(zrotate_axis_rad))
        ENDDO
     ENDDO
!$OMP END DO

     ! set the surface pressure to constant value
     !

!$OMP WORKSHARE
     pt_hydro_prog%pres_sfc = ps0
!$OMP END WORKSHARE

     ! normal wind (The meridional wind is zero by design.)
     !
!$OMP DO PRIVATE(jb,jk,je,nlen,zcos32z,lon,lat,rot_lon,rot_lat,&
!$OMP            zcosysq,zsin2ysq,tmp1,tmp2,tmp3,zu,zv)
     DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF
        DO jk = 1, nlev
           zcos32z = cos_32_etav(jk)

           DO je = 1, nlen
              lon= pt_patch%edges%center(je,jb)%lon
              lat= pt_patch%edges%center(je,jb)%lat

              IF ( ABS(p_rotate_axis_deg) < 1.0e-8_wp ) THEN
                 rot_lon = lon
                 rot_lat = lat
              ELSE
                 CALL regrot(lon,lat,rot_lon,rot_lat,0._wp,&
                             -pi_2+zrotate_axis_rad,1)
              ENDIF

              zcosysq = COS(rot_lat)*COS(rot_lat)
              zsin2ysq = SIN(2.0_wp*rot_lat)*SIN(2.0_wp*rot_lat)

              zu = u0*zcos32z*zsin2ysq - 4._wp*zcosysq
              zv = 0._wp

              IF ( ABS(p_rotate_axis_deg) >= 1.0e-8_wp ) THEN
                 tmp1 = zu
                 tmp2 = zv
                 ! rotate wind components
                 CALL turnwi(tmp1, tmp2, zu, zv,lon,lat,rot_lon,rot_lat,&
                             0.0d0,-pi_2+zrotate_axis_rad,-1)
                 IF ( ABS(zu)<1.0e-10_wp ) zu = 0._wp
                 IF ( ABS(pi_2-lat)<1.0e-8_wp .OR. ABS(pi_2+lat)<1.0e-8_wp ) &
                                           zv = 0._wp
              ENDIF

              pt_hydro_prog%vn(je,jk,jb) = &
                       zu * pt_patch%edges%primal_normal(je,jb)%v1   &
                   & + zv * pt_patch%edges%primal_normal(je,jb)%v2

              ! Coriolis parameter
              pt_patch%edges%f_e(je,jb) = 2.0_wp*omega*(SIN(lat)*COS(zrotate_axis_rad)&
                                        -COS(lon)*COS(lat)*SIN(zrotate_axis_rad))
           ENDDO ! edge loop
        ENDDO ! vertical level loop
     ENDDO ! block loop
!$OMP END DO

   ! temperature

!$OMP DO PRIVATE(jb,jk,jc,nlen,zeta,zsinz,zcos12z,zcos32z,ztemp,lon,lat,rot_lon,&
!$OMP            rot_lat,zsiny,zcosy,tmp1,tmp2,tmp3,tmp4)
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF
        DO jk = 1, nlev

           zeta    = ceta(jk)
           zsinz   = sin_etav(jk)
           zcos12z = cos_12_etav(jk)
           zcos32z = cos_32_etav(jk)
           ztemp   = temp_avg(jk)

           DO jc = 1, nlen
              lon   = pt_patch%cells%center(jc,jb)%lon
              lat   = pt_patch%cells%center(jc,jb)%lat

              IF ( ABS(p_rotate_axis_deg) < 1.0e-8_wp ) THEN
                 rot_lon = lon
                 rot_lat = lat
              ELSE
                 CALL regrot(lon,lat,rot_lon,rot_lat,0._wp,&
                             -pi_2+zrotate_axis_rad,1)
              ENDIF

              zsiny = SIN(rot_lat)
              zcosy = COS(rot_lat)

              tmp1  = 1.5_wp*zeta*pi_2*u0/rd * zsinz*zcos12z
              tmp2  = 2.0_wp*u0*zcos32z*(-2._wp*zsiny**6*(1._wp/3._wp+zcosy**2) +10._wp/63._wp)
              tmp3  = omega*re*(8._wp/5._wp*zcosy**3*(2._wp/3._wp+zsiny**2) - pi/4._wp)
              tmp4  = 8._wp*(1._wp-2._wp/3._wp*zsiny**2)*zsiny**4 - 8.8_wp/10.5_wp

              pt_hydro_prog%temp(jc,jk,jb) = ztemp + tmp1*(tmp2+tmp3+tmp4)

           ENDDO ! cell loop
        ENDDO ! vertical level loop
     ENDDO ! block loop
!$OMP END DO

     ! tracer
     IF ( ltransport ) THEN

!$OMP DO PRIVATE(jb,jk,jt,jc,nlen,zeta,ctracer,lon,lat,zrhf,zsqv)
        DO jb = 1, nblks_c
           IF (jb /= nblks_c) THEN
              nlen = nproma
           ELSE
              nlen = npromz_c
           ENDIF
           ! compute half-level pressure
           CALL half_level_pressure( pt_hydro_prog%pres_sfc(:,jb), &! input
&                                                     nproma,nlen, &! input
&                                    pt_hydro_diag%pres_ic(:,:,jb)) ! output
           ! compute pressure value at full levels
           CALL full_level_pressure( pt_hydro_diag%pres_ic(:,:,jb),&! input
&                                                       nproma, nlen,  &
&                                   pt_hydro_diag%pres_mc    (:,:,jb)) ! output

           DO jk = 1, nlev

              zeta = ceta(jk)

              DO jt = 1, ntracer

                 ctracer = ctracer_list(jt:jt)

                 SELECT CASE(ctracer)

                 CASE('1')

                    DO jc = 1, nlen
                       lat = pt_patch%cells%center(jc,jb)%lat
                       lon = pt_patch%cells%center(jc,jb)%lon
                       pt_hydro_prog%tracer(jc,jk,jb,jt) =  &
                         tracer_q1_q2(lon, lat, zeta, p_rotate_axis_deg, 0.6_wp)
                    ENDDO ! cell loop

                 CASE('2')

                    DO jc =1, nlen
                       lat = pt_patch%cells%center(jc,jb)%lat
                       lon = pt_patch%cells%center(jc,jb)%lon
                       pt_hydro_prog%tracer(jc,jk,jb,jt) =  &
                         tracer_q1_q2(lon, lat, zeta, p_rotate_axis_deg, 1.0_wp)
                    ENDDO ! cell loop

                 CASE('3')

                    DO jc =1, nlen
                       lat = pt_patch%cells%center(jc,jb)%lat
                       lon = pt_patch%cells%center(jc,jb)%lon
                       pt_hydro_prog%tracer(jc,jk,jb,jt) =  &
                         tracer_q3(lon, lat, p_rotate_axis_deg)
                    ENDDO ! cell loop

                 CASE('4')
                    pt_hydro_prog%tracer(:,jk,jb,jt) = 1._wp

                 END SELECT

                 IF (iforcing==ildf_echam) THEN
                   IF(jt == iqv ) THEN

                     DO jc =1, nlen

                       IF (lrh_linear_pres) THEN  ! rel. humidity is a linear func. of pressure
                         zrhf =  rh_at_1000hpa - 0.5_wp                      &
                              & +pt_hydro_diag%pres_mc(jc,jk,jb)/200000._wp
                         zrhf = MAX(0._wp,zrhf)
                       ELSE
                        !IF( pt_hydro_diag%pres_mc(jc,jk,jb) <= 50000._wp) zrhf = 0.2_wp
                         IF( pt_hydro_diag%pres_mc(jc,jk,jb) <= 70000._wp) zrhf = 0.55_wp
                         IF( pt_hydro_diag%pres_mc(jc,jk,jb) >  70000._wp) zrhf = 0.9_wp
                       ENDIF

                       IF ( pt_hydro_prog%temp(jc,jk,jb) <= tmelt) THEN
                          zsqv = spec_humi( sat_pres_ice( pt_hydro_prog%temp(jc,jk,jb) ), &
                                        & pt_hydro_diag%pres_mc(jc,jk,jb)            )
                       ELSE
                          zsqv = spec_humi( sat_pres_water( pt_hydro_prog%temp(jc,jk,jb) ),  &
                                        & pt_hydro_diag%pres_mc(jc,jk,jb)             )
                       END IF

                       pt_hydro_prog%tracer(jc,jk,jb,jt) = zrhf*zsqv

                       IF( pt_hydro_diag%pres_mc(jc,jk,jb) <= 10000._wp) &
                           pt_hydro_prog%tracer(jc,jk,jb,jt) &
                           & = MIN ( 5.e-6_wp, pt_hydro_prog%tracer(jc,jk,jb,jt) )
                     END DO
!                    flev= FLOAT(nlev)
!                    flev2 = flev*flev
!                    fheight=(flev-(jk-1))**1.5
!                    pt_hydro_prog%tracer(:,jk,jb,jt) = 3.e-3_wp*EXP(-fheight/flev)
!
                     IF(jb == 1) THEN
                       WRITE(message_text,'(a,i4,f15.4,f15.10,f10.5)') &
                         & 'jwtest: qv',jk,                            &
                         & pt_hydro_diag%pres_mc(nlen,jk,jb)*0.01_wp,  &
                         & pt_hydro_prog%tracer(nlen,jk,jb,jt),        &
                         & pt_hydro_prog%temp(nlen,jk,jb)
                       CALL message('', TRIM(message_text))
                     ENDIF

                   ELSE IF (jt==iqt) THEN
                     pt_hydro_prog%tracer(:,jk,jb,jt) = 1._wp
                   ELSE ! other hydrometeors
                     pt_hydro_prog%tracer(:,jk,jb,jt) = 0._wp
                   ENDIF !tracer
                 ENDIF   !iforcing

              ENDDO ! tracer loop
           ENDDO ! vertical level loop
        ENDDO ! block loop
!$OMP END DO
     ENDIF

!$OMP END PARALLEL


 END SUBROUTINE ldf_init_prog_state

!-------------------------------------------------------------------------
END MODULE mo_ldf_init



