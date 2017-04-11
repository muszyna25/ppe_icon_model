MODULE mo_nwp_ww

!>
!!  Module to calculate weather interpretation, WW.
!!
!! @author Helmut Frank, DWD
!!
!! $Id: n/a$
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

#ifdef ONLYWW

  USE mo_wwonly, ONLY: wp, vct_a, &
                       rdv, O_m_rdv, rd_o_cpd, tmelt, alvdcp, b3, b1, b2w, b4w, &
                       max_dom, kstart_moist

#else

  USE mo_kind,               ONLY: wp

  USE mo_physical_constants, ONLY: rdv, O_m_rdv, tmelt,    &
                                   rd_o_cpd, alvdcp,          &
                                   b3    => tmelt
  USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                   b2w   => c3les , & !! pressure over water (l) and ice (i)
                                   b4w   => c4les     !!               -- " --
  USE mo_impl_constants,     ONLY: max_dom

  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_vertical_coord_table,  ONLY: vct_a

  USE mo_exception,             ONLY: message, message_text
  USE mtime,                    ONLY: datetime, newDatetime, timeDelta

#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: configure_ww, ww_diagnostics
#ifdef ONLYWW
  PUBLIC :: rkgrenz1, rkgrenz2
#else
  PUBLIC :: ww_datetime

  TYPE t_datetime_ptr
    TYPE(datetime), POINTER :: ptr 
  END TYPE t_datetime_ptr

  TYPE(t_datetime_ptr) :: ww_datetime(max_dom)
#endif

!  The following indices are used in the diagnostics of WW
  INTEGER :: ih_500hPa(max_dom),  &    ! index of layer next to 500 hPa for US standard atmosphere
             ih_700hPa(max_dom),  &    ! index of layer next to 700 hPa for US standard atmosphere
             ih_850hPa(max_dom),  &    ! index of layer next to 850 hPa for US standard atmosphere
             ih_950hPa(max_dom),  &    ! index of layer next to 950 hPa for US standard atmosphere
             ihb500hPa(max_dom),  &    ! index of next layer below 500 hPa
             ihb950hPa(max_dom)        ! index of next layer below 950 hPa

! The following indices are used to calculate a fog index
  INTEGER :: ifog_temp(max_dom),  &    ! calculates minimum spread in this layer
             ifog_wind(max_dom)        ! calculate mean wind in this layer

  REAL(wp) :: rkgrenz1,   &    ! limit for thunderstorm: convective precipitation in 3 hours
              rkgrenz2         ! limit for strong thunderstorm: convective precip. in 3 hours

CONTAINS

  SUBROUTINE configure_ww( ini_datetime, jg, nlev, nshift_total)
!
!   Set some constants for the calculation of WW
!
    CHARACTER(LEN=24), PARAMETER :: routine = 'atm_phy_nwp:configure_ww'

    TYPE(datetime),   INTENT(IN), POINTER     :: ini_datetime  ! init datetime (mtime)
    INTEGER,  INTENT(IN) :: jg           !< patch 
    INTEGER,  INTENT(IN) :: nlev         !< number of full vertical levels 
    INTEGER,  INTENT(IN) :: nshift_total 

!   The following values are calculated as height difference between the height of the pressure level
!   and the height of 1000 hPa for a US standard atmosphere; h_1000hPa=110.88 m
    REAL(wp), PARAMETER :: h_500hPa = 5463.56_wp    ! height of 500 hPa level
    REAL(wp), PARAMETER :: h_700hPa = 2091.30_wp    ! height of 700 hPa level
    REAL(wp), PARAMETER :: h_850hPa = 1346.42_wp    ! height of 850 hPa level
    REAL(wp), PARAMETER :: h_950hPa =  429.46_wp    ! height of 950 hPa level
!   The followin heights are used to calculate a fox index
    REAL(wp), PARAMETER :: hfog_temp = 877.62_wp   ! h(900 hPa)-h(1000 hPa)
    REAL(wp), PARAMETER :: hfog_wind = 298.74_wp   ! h(965 hPa)-h(1000 hPa)

    INTEGER             :: jk, jk1

#ifndef ONLYWW
    ww_datetime(jg)%ptr => newDatetime(ini_datetime)
#endif

!   Find next layer below height of 500 hPa, and next layer to height of 500 hPa
    ihb500hPa(jg) = kstart_moist(jg)
    DO jk = MAX( kstart_moist(jg),2) , nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_500hPa) THEN
        ihb500hPa(jg) = jk
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_500hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_500hPa) ) THEN
          ih_500hPa(jg) = jk
        ELSE
          ih_500hPa(jg) = jk-1
        END IF
        EXIT
      END IF
    ENDDO

    ih_700hPa(jg) = ihb500hPa(jg)
    DO jk = ihb500hPa(jg), nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_700hPa) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_700hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_700hPa) ) THEN
          ih_700hPa(jg) = jk
        ELSE
          ih_700hPa(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO

    ih_850hPa(jg) = ih_700hPa(jg)+1
    DO jk = ih_700hPa(jg)+1, nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_850hPa) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_850hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_850hPa) ) THEN
          ih_850hPa(jg) = jk
        ELSE
          ih_850hPa(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO

    ifog_temp(jg) = ih_850hPa(jg)+1
    DO jk = ih_850hPa(jg)+1, nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < hfog_temp) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - hfog_temp) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - hfog_temp) ) THEN
          ifog_temp(jg) = jk
        ELSE
          ifog_temp(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO

    ihb950hPa(jg) = ifog_temp(jg)
    DO jk = ifog_temp(jg), nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_950hPa) THEN
        ihb950hPa(jg) = jk
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_950hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_950hPa) ) THEN
          ih_950hPa(jg) = jk
        ELSE
          ih_950hPa(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO

    ifog_wind(jg) = ih_950hPa(jg)
    DO jk = ih_950hPa(jg), nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < hfog_wind) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - hfog_wind) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - hfog_wind) ) THEN
          ifog_wind(jg) = jk
        ELSE
          ifog_wind(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO

#ifdef ONLYWW
    WRITE( *,'(a,i4,a,6i4)') 'Domain', jg,                                    &
  &    '; ih_500hPa, ihb500hPa, ih_700hPa, ih_850hPa, ih_950hPa, ihb950hPa ', &
  &       ih_500hPa(jg), ihb500hPa(jg), ih_700hPa(jg), ih_850hPa(jg),         &
  &       ih_950hPa(jg), ihb950hPa(jg)
    WRITE( *,'(2(a,i4),i4)') 'Domain', jg,                                    &
  &    '; ifog_temp, ifog_wind ', ifog_temp(jg), ifog_wind(jg)
#else
    WRITE( message_text,'(a,i4,a,6i4)') 'Domain', jg,                         &
  &    '; ih_500hPa, ihb500hPa, ih_700hPa, ih_850hPa, ih_950hPa, ihb950hPa ', &
  &       ih_500hPa(jg), ihb500hPa(jg), ih_700hPa(jg), ih_850hPa(jg),         &
  &       ih_950hPa(jg), ihb950hPa(jg)
    CALL message( TRIM(routine), message_text)
    WRITE( message_text,'(2(a,i4),i4)') 'Domain', jg,                         &
  &    '; ifog_temp, ifog_wind ', ifog_temp(jg), ifog_wind(jg)
    CALL message( TRIM(routine), message_text)
#endif

  END SUBROUTINE configure_ww


  SUBROUTINE ww_diagnostics( ie, ke, ke1, i_startidx, i_endidx, jg,       &
                             t   , qv   , qc , u   , v   , pf   , ph,     &
                             t_2m, td_2m, t_g, clct, clcm, u_10m, v_10m,  &
                             rain_gsp0, rain_gsp, rain_con0, rain_con,    &
                             snow_gsp0, snow_gsp, snow_con0, snow_con,    &
                             bas_con, top_con, time_diff, ymodel, iww)
!
! Description:
! The subroutime evaluates for each grid point of the model a "ww-number";
! that means, a number corresponding to those numbers that
! would be used for weather observation (e.g. thunderstorm=95, 
! heavy rain=65 etc.)
! This subroutine is based on auswert_lm.f of the ww-library of DWD.


  INTEGER, INTENT(IN)  :: ie, ke, ke1, &        ! horizontal and vertical dimensions
                          i_startidx, i_endidx  ! start and end indices of loops in horizontal
  INTEGER, INTENT(IN)  :: jg                    !  patch 

  REAL(wp), INTENT(IN) :: t (ie,ke),  &  ! temperature (K)
                          qv(ie,ke),  &  ! specific humidity (on model levels) (kg/kg)
                          qc(ie,ke),  &  ! specific cloud liquid water (on model levels) (kg/kg)
                          u (ie,ke),  &  ! zonal comp. of wind vel. (m/s)
                          v (ie,ke),  &  ! meridional comp. of wind vel. (m/s)
                          pf(ie,ke),  &  ! pressure on full levels (Pa)
                          ph(ie,ke1)     ! pressure on half levels (Pa)

  REAL(wp), INTENT(IN) :: t_2m (ie),  &  ! temperature at 2m above the ground (K)
                          td_2m(ie),  &  ! dewpoint temperature at 2m above the ground (K)
                          t_g  (ie),  &  ! surface temperature (weighted from t_s and t_snow) (K)
                          clct (ie),  &  ! total cloud cover
                          clcm (ie),  &  ! medium cloud cover
                          u_10m(ie),  &  ! zonal wind component at 10m above the ground (m/s)
                          v_10m(ie),  &  ! meridional wind component at 10m above the ground (m/s)
                          rain_gsp0(ie), & ! grid scale rain, valid at previous call of ww_diagnostics
                          rain_gsp (ie), & ! grid scale rain, valid now
                          rain_con0(ie), & ! convective rain, valid at previous call of ww_diagnostics
                          rain_con (ie), & ! convective rain, valid now
                          snow_gsp0(ie), & ! grid scale snow, valid at previous call of ww_diagnostics
                          snow_gsp (ie), & ! grid scale snow, valid now
                          snow_con0(ie), & ! convective snow, valid at previous call of ww_diagnostics
                          snow_con (ie)    ! convective snow, valid now

  INTEGER,  INTENT(IN) :: bas_con(ie),   & ! base index of main convective cloud
                          top_con(ie)      ! top index of main convective cloud

  TYPE(timeDelta),        POINTER :: time_diff ! time since previous call
  CHARACTER(LEN=4), INTENT(IN) :: ymodel   ! Use scalings factors for this NWP model

  INTEGER, INTENT(OUT) :: iww(ie)          ! significant weather code

    REAL(wp) :: rgdiff         !  Difference total precip. in mm/hour
    REAL(wp) :: rkdiff         !  Difference convective precip. in mm/hour
    REAL(wp) :: rldiff         !  Difference grid scale precip. in mm/hour
    REAL(wp) :: rrdiff         !  Difference rain in mm/hour
    REAL(wp) :: rsdiff         !  Difference snow in mm/hour

    INTEGER :: i, k
    INTEGER :: iwolk, iwolkc

    REAL(wp) :: rko            ! KO-Index
    REAL(wp) :: rkgr_1, rkgr_2 ! limits for rkdiff
    REAL(wp) :: vbetr          ! wind speed in 10 m a. g.
    REAL(wp) :: tblmax         ! max. value of t in lowest 100 hPa
    REAL(wp) :: tdbl, tdblmax  ! value/max. value of td in lowest 100 hPa
    REAL(wp) :: dt_ke          ! spread in lowest model level
    REAL(wp) :: vbl            ! mean value of wind speed in lowest 35 hPa
    REAL(wp) :: neb_i          ! Nebulae-Index
    REAL(wp) :: zfrac, qvmin   ! variables used in determination of dew point temperature
    REAL(wp) :: sf1            ! scale factor for showers
    REAL(wp) :: sf2            ! time scaling factor
    REAL(wp) :: dp             ! pressure difference
    INTEGER  :: igfb, irrb, isprb
    REAL(wp) :: dhour

    REAL(wp), PARAMETER :: rkogrenz = 1._wp          ! limit for rko
    REAL(wp), PARAMETER :: ms_kn = 1._wp/0.51444_wp  ! factor to convert numbers in m/s to knots

! preset some constants for calculation of TD and for conversi
    sf1 = 1._wp

    dhour = time_diff%day*24. + time_diff%hour + REAL(time_diff%minute,wp)/60.
    sf2 = dhour/3._wp

    SELECT CASE (ymodel)
      CASE( 'ICON')  ! these values are the same as for COSMO-EU (LME)
        rkgrenz1 = 2.25_wp
        rkgrenz2 = 6.0_wp
      CASE( 'GME ')
        rkgrenz1 = 0.6_wp
        rkgrenz2 = 3.0_wp
      CASE( 'GME6')         ! old GME with mesh width <= 60 km (NI<=128)
        rkgrenz1 = 0.75_wp
        rkgrenz2 = 3.0_wp
      CASE( 'LME ')
        rkgrenz1 = 2.25_wp
        rkgrenz2 = 6.0_wp
      CASE( 'LMK ')
        rkgrenz1 = 1.80_wp
        rkgrenz2 = 6.0_wp
      CASE( 'TEST')
        PRINT *, 'Test ww_diagnostics with rkgrenz1=', rkgrenz1, '  rkgrenz2=',rkgrenz2
      CASE DEFAULT
        PRINT *, 'Unknown model ', ymodel,' in ww_diagnostics'
        iww(:) = -1
        RETURN
    END SELECT

!     Limits for convective precipitation

    rkgr_1 = sf2 * rkgrenz1
    rkgr_2 = sf2 * rkgrenz2

    DO i = i_startidx, i_endidx
      iww(i) = -9
      rgdiff =  rain_gsp (i) + snow_gsp (i) + rain_con (i) + snow_con (i)     &
   &          - rain_gsp0(i) - snow_gsp0(i) - rain_con0(i) - snow_con0(i)

!.... Total precipitation < 0.045 mm (pro 3h) => weather without precipitation
!        IF ( rgdiff < (sf2*0.045_wp) ) GO TO 400
WW_PRECIP: IF ( rgdiff >= (sf2*0.045_wp) ) THEN

         rkdiff = rain_con(i) + snow_con(i) - rain_con0(i) - snow_con0(i)
         rldiff = rain_gsp(i) + snow_gsp(i) - rain_gsp0(i) - snow_gsp0(i)
         rrdiff = rain_gsp(i) + rain_con(i) - rain_gsp0(i) - rain_con0(i)
         rsdiff = snow_gsp(i) + snow_con(i) - snow_gsp0(i) - snow_con0(i)

!.... Thunderstorm

         IF (top_con(i) > 0 .AND. top_con(i) <= ke1) THEN
            dp = ph(i,bas_con(i)) - pf(i,top_con(i))
         ELSE
            dp   = 0.0_wp
         ENDIF

!     Determination of KO-Index (rko)

         CALL ko_index( pf(i,ih_950hPa(jg)), pf(i,ih_850hPa(jg)), pf(i,ih_700hPa(jg)), pf(i,ih_500hPa(jg)),  &
     &                   t(i,ih_950hPa(jg)),  t(i,ih_850hPa(jg)),  t(i,ih_700hPa(jg)),  t(i,ih_500hPa(jg)),  &
     &                  qv(i,ih_950hPa(jg)), qv(i,ih_850hPa(jg)), qv(i,ih_700hPa(jg)), qv(i,ih_500hPa(jg)),  &
     &                  rko)

         IF ( rkdiff > rkgr_1 .AND. (top_con(i) > 0 .AND. top_con(i) <= ke1)) THEN
            IF ( dp > 400.e2_wp .AND. t(i, top_con(i)) < tmelt-25._wp) THEN
               IF (rko < rkogrenz) THEN
                  IF (rko < -6._wp .AND. t(i, top_con(i)) < tmelt-45..AND.  rkdiff > rkgr_2) THEN
                    iww(i) = 96
                  ELSE
                    iww(i) = 95
                  ENDIF
                  CYCLE
               ENDIF
            ENDIF
         ENDIF

         CALL gefr( jg, ke, ke1, ph(i,:), qc(i,:), t(i,:), t_2m(i), t_g(i), &
                    igfb, irrb, isprb )

!.... Convective precipitation

CPRECIP: IF ( rkdiff > rldiff) THEN
            IF ( ( (rsdiff > rrdiff .AND. t_2m(i) <= tmelt+5._wp) .OR.  &
       &         t_2m(i) < tmelt-10._wp)  .OR. irrb == 0) THEN
!....          Snow shower
               IF (     rgdiff >= (sf1 * sf2 * 0.75_wp)) THEN
                 iww(i) = 86
               ELSE IF (rgdiff >= (sf1 * sf2 * 0.10_wp)) THEN
                 iww(i) = 85
               ELSE
!              Clouds
                 iww(i) = clct2ww( clct(i) )
               END IF
            ELSE
               IF (igfb == 1) THEN
!....             Freezing rain (without drizzle)
                  IF (rgdiff > (sf2 * 1.5_wp)) THEN
                    iww(i) = 67
                  ELSE
                    iww(i) = 66
                  END IF
               ELSE
!....             Rain shower 
                  IF (     rgdiff > (sf1 * sf2 * 24._wp)) THEN
                    iww(i) = 82
                  ELSE IF (rgdiff > (sf1 * sf2 * 1.5_wp)) THEN
                    iww(i) = 81
                  ELSE IF (rgdiff > (sf1 * sf2 * 0.6_wp)) THEN
                    iww(i) = 80
                  ELSE
!                   Clouds
                    iww(i) = clct2ww( clct(i) )
                  END IF
               ENDIF
            ENDIF

         ELSE CPRECIP

!.... Large-scale precipitaion

            iwolk = 0
            iwolkc = 0
            DO k = ihb950hPa(jg)+1, ke
               IF ( qc(i,k) > 1.e-8_wp) THEN
                  iwolkc = 1
                  EXIT
               ENDIF
            ENDDO
            IF (iwolkc == 1.AND. clcm(i) < 0.95_wp) iwolk = 1

            IF ( ((rsdiff > rrdiff .AND. t_2m(i) <= tmelt+5._wp) .OR.t_2m(i) < tmelt-10._wp)  &
     &       .OR.( irrb == 0 .AND. rgdiff >  (sf2*0.6_wp) )                                   &
     &       .OR.( irrb == 0 .AND. rgdiff <= (sf2*0.6_wp) .AND.iwolk == 0)                    &
     &       .OR.( igfb  > 0 .AND. rgdiff <= (sf2*0.6_wp) .AND.iwolk == 1.AND.isprb == 0) ) THEN

!.... Snow
               IF (rgdiff > (sf2 * 6._wp)) THEN
                 iww(i) = 75
               ELSE IF (rgdiff > (sf2 * 0.75_wp)) THEN
                 iww(i) = 73
               ELSE
                 iww(i) = 70
               END IF
            ELSE
               IF (igfb > 0) THEN

!.... Freezing rain (with drizzle)

                  IF (rgdiff <= (sf2 * 0.6_wp).AND.iwolk == 1) THEN
                   iww(i) = 56
                  ELSE IF (rgdiff <= (sf2 * 1.5_wp)) THEN
                    iww(i) = 66
                  ELSE
                    iww(i) = 67
                  END IF
               ELSE

!.... Rain
                  IF (rgdiff <= (sf2 * 0.6_wp).AND.iwolk == 1) THEN
                    iww(i) = 50
                  ELSE IF (rgdiff <= (sf2 * 1.5_wp)) THEN
                    iww(i) = 60
                  ELSE IF (rgdiff <= (sf2 * 12._wp)) THEN
                    iww(i) = 63
                  ELSE
                    iww(i) = 65
                  END IF
               ENDIF
            ENDIF
         ENDIF CPRECIP

      ELSE WW_PRECIP

!.... Weather without precipitation

!        calculate fog index
         vbetr   = SQRT (u_10m(i)**2 + v_10m(i)**2)
         tblmax  = t_2m (i)
         tdblmax = td_2m(i)
         vbl = 0._wp
         DO k = ifog_wind(jg), ke
           vbl = vbl + SQRT(u(i,k)**2 + v(i,k)**2)*(ph(i,k+1)-ph(i,k))
         ENDDO
         vbl = vbl / ( ph(i,ke1)-ph(i,ifog_wind(jg)) )
         DO k = ifog_temp(jg), ke
           tblmax = MAX( t(i,k), tblmax)
           qvmin = MAX( qv(i,k), 1.e-12_wp)
           zfrac = LOG( pf(i,k)*qvmin/(b1*(rdv+O_m_rdv*qvmin)))
           tdbl = (b2w*b3-b4w*zfrac)/(b2w-zfrac)
           tdbl = MIN( tdbl, t(i,k))
           tdblmax = MAX( tdbl, tdblmax)
         ENDDO
         dt_ke = t(i,ke) - tdbl
         neb_i = ( t_2m(i)-tblmax) + ( t_2m(i)-tdblmax) + ms_kn*vbl

         IF ( neb_i < 5._wp .AND. dt_ke < 2._wp) THEN
!           Fog
            IF ( (t_2m(i) < tmelt-1.0_wp .OR. t_g(i) < tmelt-1.0_wp)   &
        &       .AND. vbetr > 1.5_wp) THEN
              iww(i) = 48
            ELSE
              iww(i) = 45
            END IF

         ELSE
!          Clouds
           iww(i) = clct2ww( clct(i) )
         END IF

      END IF WW_PRECIP
   
    ENDDO

  END SUBROUTINE ww_diagnostics


  SUBROUTINE ko_index ( p950,  p850,  p700,  p500,      &
     &                  t950,  t850,  t700,  t500,      &
     &                 qv950, qv850, qv700, qv500,      &
     &                 rko)
!
! Description:
! Calculation of thunderstorm indices KO as for COSMO-EU as in subroutine index.f in ww-library.
!
! Method:
! Calculation of the equivalent potential temperature: then,  
! using the appropriate formulas, calculation of KO index.
!

!=======================================================================
!
      REAL(wp), INTENT(IN) :: p950, p850, p700, p500       ! single values of pressure field p in hPa
      REAL(wp), INTENT(IN) :: t950, t850, t700, t500       ! single values of temperature field t in degree C
      REAL(wp), INTENT(IN) :: qv950, qv850, qv700, qv500   ! single values of specific humidity field qv

      REAL(wp), INTENT(OUT)::  rko                         ! KO-Index
!-----------------------------------------------------------------------
!
      REAL(wp), PARAMETER :: qvmin = 1.e-12_wp            ! minimum value for specific humidity
      REAL(wp)            :: ta950, ta850, ta700, ta500   ! equivalent potential temperature
!
  
!     950 hPa
      ta950 = ( t950 + alvdcp*MAX( qvmin, qv950) )*(1.e5/p950)**rd_o_cpd

!     850 hPa
      ta850 = ( t850 + alvdcp*MAX( qvmin, qv850) )*(1.e5/p850)**rd_o_cpd

!     700 hPa
      ta700 = ( t700 + alvdcp*MAX( qvmin, qv700) )*(1.e5/p700)**rd_o_cpd

!     500 hPa
      ta500 = ( t500 + alvdcp*MAX( qvmin, qv500) )*(1.e5/p500)**rd_o_cpd

      rko = 0.5_wp *( ta700 + ta500 - ta950 - ta850 )

    END SUBROUTINE ko_index


    SUBROUTINE gefr ( jg, ke, ke1, p_gefr, qc_gefr, t_gefr, t2m, tg, &
                      igfb, irrb, isprb )
!
! Description:
! Determination of igfb, irrb and isprb (indicators for freezing)
!
! Method:
! Main part:
!    Set igfb, irrb and isprb to 0. If the main freezing condition is fulfilled,
!    set igfb to 1. From level ihb500hPa(jg) to ke (bottom level of the model) calculate
!    the products t*p, if t>0; build the sum of those products. Dependent on
!    this sum resp. the value of t_2m  and t_g set the value for irrb resp.
!    igfb.
!
! Additional conditions for drizzle:
!   Determine the index of the bottom layer of low clouds (kwug); 
!   dependent on the value of t at this level continue (or return to the
!   calling program):
!   If kwug equals ke set isprb to 1 and return to the calling program);
!   otherwise build the sum of the products t*p from kwug+1 to ke; dependent on
!   the value of this sum resp. the value of the temperature at level ke set
!   isprb and return to the calling program.
!
! Current Code Owner: DWD, Ulrich Pflueger
!   Deutscher Wetterdienst, Abt. FE1
!   Frankfurter Str. 135
!   63067 Offenbach/Main
!   Telefon: 069/8062-2753
!   FAX: 069/8062-3721
!   email: ulrich.pflueger@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V1_0         2014/02/27 Ulrich Pflueger
!  Initial release - based on V1_31 gefr.f (library wwlm)
!  and               based on V1_9  gefr.f (library wwlmk)
!
! Code Description:
! Language: Fortran 90.
!=======================================================================
!
      INTEGER, INTENT(IN)  :: jg                ! patch index
      INTEGER, INTENT(IN)  :: ke, ke1           ! number of levels
      REAL(wp), INTENT(IN) :: p_gefr(ke1)       ! pressure at half levels
      REAL(wp), INTENT(IN) :: qc_gefr(ke)       ! profile of cloud water qc
      REAL(wp), INTENT(IN) :: t_gefr(ke)        ! profile of temperature t
      REAL(wp), INTENT(IN) :: t2m, tg           ! 2m temp. t_2m, surface temp. t_g

      INTEGER, INTENT(OUT) :: igfb,         &   ! indicator for freezing
                              irrb,         &   ! indicator for rain
                              isprb             ! indicator for drizzle
!
! Local scalars:

      INTEGER kk                         ! loop variable
      INTEGER :: kwug                    ! index of bottom level of low clouds

      REAL(wp) :: sumtpr                 ! sum of products t*p for t>0 from level ihb500hPa(jg) to ke
      REAL(wp) :: sumtpsp                ! sum of products t*p from level kwug+1 to ke
      REAL(wp) :: dp                     ! help variables to calculate
                                         ! sumtpr resp. sumtpsp
      REAL(wp), PARAMETER :: xtpsp = -20.e2_wp          ! limit for sumtpsp
      REAL(wp), PARAMETER :: xtsp  = tmelt - 2.0_wp     ! limit for t at level kwug
      REAL(wp), PARAMETER :: xtpr  = 100.e2_wp          ! limit for sumtpr
!
!- End of header
!

!   **********************************************************************
!   * Limits for t und t*p (variable name ending with "r" for rain,      *
!   * ending with "sp" for drizzle)                                      *
!   **********************************************************************

      igfb = 0
      irrb = 0
      isprb = 0

!.... Freezing condition fulfilled

      IF ( t2m < tmelt-0.1_wp .OR. tg < tmelt-0.1_wp) igfb = 1

!.... Additional condition for rain

      sumtpr = 0._wp

      DO kk = ihb500hPa(jg), ke
         IF ( t_gefr(kk) > tmelt) THEN
            dp = p_gefr(kk+1)  - p_gefr(kk)
            sumtpr = sumtpr + (t_gefr(kk)-tmelt) * dp
         ENDIF
      ENDDO

      IF (sumtpr > xtpr .OR. (t2m > tmelt+0.5_wp .AND. tg >= tmelt-0.1_wp)) THEN
         irrb = 1
      ELSE
         IF (igfb == 1) igfb = 2
      ENDIF

      kwug = ihb950hPa(jg)

      DO kk = ihb950hPa(jg)+1, ke
         IF ( qc_gefr(kk) > 1.e-8) kwug = kk
      ENDDO   

!.... First additional condition

      IF ( t_gefr(kwug) <= xtsp) RETURN

!.... Return, if kwug = ke

      IF (kwug == ke) THEN
         isprb = 1
         RETURN
      ENDIF

!.... Second additional

      sumtpsp = 0._wp
      DO  kk = kwug+1, ke
         dp = p_gefr(kk+1)  - p_gefr(kk)
         sumtpsp = sumtpsp + (t_gefr(kk)-tmelt) * dp
      ENDDO

      IF (sumtpsp > xtpsp .OR. t_gefr(ke) > tmelt) THEN
         isprb = 1
      ENDIF

    END SUBROUTINE gefr

    INTEGER FUNCTION clct2ww( clct)
!
!   WW key (significant weather) for different cloud covers
!
      REAL(wp), INTENT(IN) :: clct
      IF (      clct > 0.8125_wp) THEN
        clct2ww = 3
      ELSE IF ( clct > 0.4375_wp) THEN
        clct2ww = 2
      ELSE IF ( clct > 0.0625_wp) THEN
        clct2ww = 1
      ELSE
        clct2ww = 0
      END IF
    END  FUNCTION clct2ww
    

END MODULE mo_nwp_ww
