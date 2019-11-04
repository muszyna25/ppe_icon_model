!! ===========================================================================================================================
!! Implementation of tides by computation of the Sun's and Moon's full tidal potential
!! This will be used in the pressure gradient calculation
!!
!! Authors: Kai Logemann, Helmholtz-Zentrum Geesthacht, and Leonidas Linardakis, Max-Planck-Institute for Meteorology, Hamburg
!!
!! ===========================================================================================================================

!----------------------------
#include "icon_definitions.inc"
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_tides
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, dp
   USE mtime,                    ONLY: datetime, datetimeToString, deallocateDatetime,              &
       &                               timedelta, newTimedelta, deallocateTimedelta,                &
       &                               MAX_DATETIME_STR_LEN, newDatetime,                           &
       &                               MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,           &
       &                               OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),          &
       &                               ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),     &
       &                               event, eventGroup, newEvent,                                 &
       &                               addEventToEventGroup, isCurrentEventActive
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_ocean_nml,              ONLY: n_zlev, tide_startdate, tides_esl_damping_coeff
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state, &
    & t_operator_coeff
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,        ONLY: nproma
  USE mo_impl_constants,         ONLY: sea_boundary
  USE mo_ocean_math_operators,   ONLY: grad_fd_norm_oce_2d_3d

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: tide
  
  CHARACTER(LEN=12)  :: str_module = 'tides'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
  !-------------------------------------------------------------------------

  
CONTAINS

  
  SUBROUTINE tide(patch_3d,mtime_current,tides_potential)
  
  IMPLICIT NONE
  
  TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
  TYPE(datetime), POINTER              :: mtime_current
  REAL(wp), INTENT(inout)              :: tides_potential(:,:)
  
  TYPE(t_patch), POINTER :: patch_2D
  TYPE(t_subset_range), POINTER :: all_cells
  INTEGER :: jb,jc,jk,je,start_index,end_index,start_edge_index,end_edge_index,icel1,icel2
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: datestri
  INTEGER  :: jahr,monat,tag,stunde,minute
  DOUBLE PRECISION :: sekunde,gezhochfahr,rs,rm
  DOUBLE PRECISION :: dtim,gst,dcl_m,alp_m,h_m,dcl_s,alp_s,h_s,mpot,spot,longitude,latitude
  double precision tide_start
  
  
  
  patch_2D        => patch_3d%p_patch_2d(1)
  all_cells       => patch_2D%cells%ALL
  
  read(tide_startdate(1:4),'(i4)') jahr
  read(tide_startdate(6:7),'(i2.2)') monat
  read(tide_startdate(9:10),'(i2.2)') tag
  read(tide_startdate(12:13),'(i2.2)') stunde
  read(tide_startdate(15:16),'(i2.2)') minute
  sekunde = 0.d0
  
  call timing(jahr,monat,tag,stunde,minute,sekunde,tide_start)  ! start time of tidal spin-up (number of days since 2010 January 0.0)
  
  CALL datetimeToString(mtime_current, datestri)
  
  READ(datestri(1:4),'(i4)') jahr
  READ(datestri(6:7),'(i2.2)') monat
  READ(datestri(9:10),'(i2.2)') tag
  READ(datestri(12:13),'(i2.2)') stunde
  READ(datestri(15:16),'(i2.2)') minute
  READ(datestri(18:23),'(f6.3)') sekunde
  
  call timing(jahr,monat,tag,stunde,minute,sekunde,dtim)  ! compute the number of days since 2010 January 0.0
  call get_gst(jahr,monat,tag,stunde,minute,sekunde,gst) ! compute the Greenwich siderial time gst
  call moon_declination(dtim,gst,dcl_m,alp_m,h_m) ! declination (deg) and right ascension (hours) and Greenwich hour angle h (hours) of the Moon after Duffet and Zwart (1979)
  call sun_declination(dtim,gst,dcl_s,alp_s,h_s) ! declination (deg) dcl_s, right ascension (hours) alp_s and Greenwich hour angle h (hours) of the Sun after Duffet and Zwart (1979)
  
  call dist_earth_sun_moon(jahr,monat,tag,stunde+minute*1./60.,rs,rm)  ! distance between rs = Earth - Sun and rm = Earth - Moon  (m)
  
  
   gezhochfahr = min(1.d0,(dtim-tide_start)/30.d0)  ! spin-up over the first 30 days, when starting at "tide_start"
!  gezhochfahr = min(1.d0,(dtim+3286.d0)/30.d0)  ! spin-up over the first 30 days, when starting at 2001-01-01 00:00
!  gezhochfahr = min(1.d0,(dtim+2160.d0)/30.d0)  ! spin-up over the first 30 days, when starting at 2004-02-01 00:00
  
! Compute tidal potential ==================================================
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, start_index,end_index,longitude,latitude,mpot,spot) ICON_OMP_DEFAULT_SCHEDULE

  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, start_index, end_index)

    DO jc = start_index, end_index

      longitude = patch_2d%cells%center(jc,jb)%lon
      latitude = patch_2d%cells%center(jc,jb)%lat           
      call Moon_pot(dcl_m,h_m,latitude,longitude,rm,mpot) ! computation of the Moon's tidal potential
      call Sun_pot(dcl_s,h_s,latitude,longitude,rs,spot)  ! computation of the Sun's tidal potential

      tides_potential(jc,jb) = tides_esl_damping_coeff * gezhochfahr * (mpot + spot) ! combined tidal potential 
  
    END DO
  END DO
!ICON_OMP_END_PARALLEL_DO
  
  
  RETURN
END SUBROUTINE tide
  
!-------------------------------------------------------------------------
! find the number of days since 2010 January 0.0
subroutine timing(year,mon,day,hour,minute,sekunde,d)
implicit none

double precision d,d0,sekunde
integer year,mon,day,hour,minute

integer molen(12),ylen(1904:2100),y,m

data molen /31,28,31,30,31,30,31,31,30,31,30,31/



do y = 1904,2100,1
  ylen(y) = 365
enddo
do y = 1904,2100,4
  ylen(y) = 366
  if(year.eq.y) molen(2) = 29
enddo


d0 = 0.
do y = 1904,2009,1
  d0 = d0 + ylen(y)
enddo

d = 0.
do y = 1904,year-1,1
  d = d + ylen(y)
enddo


do m = 1,mon-1,1
  d = d + molen(m)
enddo

d = d + day
d = d + hour*1.d0/24.
d = d + minute*1.d0/(24.*60)
d = d + sekunde*1.d0/(24.*60*60)

d = d - d0

return
end subroutine timing

! ==========================================================================
! declination (deg) and right ascension (hours) of the Moon after Duffet and Zwart (1979)

subroutine moon_declination(d,gst,dcl,alp,h)
implicit none

double precision d,dcl,alp,gst,h
double precision ma,lama
double precision epsg,omg,e,pi,rad,l,l0,mm,n,p0,n0,ev,ae,a3,c
double precision mms,ec,a4,ls,v,lss,ns,i,lamm,betm,y,x,at,t,eps

pi = 3.141592653589793d0
rad = pi/180.d0


epsg = 279.557208d0  ! ecliptic longitude at epoch 2010.0
omg = 283.112438d0   ! ecliptic longitude of perigee at epoch 2010.0
e = 0.016705d0       ! eccentricity of orbit at epoch 2010.0

ma = 360.d0/365.242191d0*d + epsg - omg

ma = mod(ma,360.d0)
if(ma.lt.0.) ma = ma + 360.d0   ! Sun's mean anomaly


lama = 360.d0/365.242191d0*d + 360.d0/pi*e*dsin((360.d0/365.242191d0*d + epsg - omg)*rad) + epsg

lama = mod(lama,360.d0)
if(lama.lt.0.) lama = lama + 360.d0   ! Sun’s ecliptic longitude


! ----------------------------------------------------------------------------------

l0 = 91.929336d0  ! Moon’s mean longitude at the epoch
p0 = 130.143076d0 ! mean longitude of the perigee at the epoch
n0 = 291.682547d0 ! mean longitude of the node at the epoch

l = 13.1763966d0*d + l0  ! Moon’s mean longitude
l = mod(l,360.d0)
if(l.lt.0.) l = l + 360.d0

mm = l - 0.1114041d0*d - p0  ! Moon’s mean anomaly
mm = mod(mm,360.d0)
if(mm.lt.0.) mm = mm + 360.d0

n = n0 - 0.0529539*d  ! ascending node’s mean longitude
n = mod(n,360.d0)
if(n.lt.0.) n = n + 360.d0


! --------------------------------------------------------------

C = l-lama

Ev = 1.2739*dsin(rad*(2*C-Mm))  ! corrections for evection

Ae = 0.1858*dsin(rad*Ma) ! annual equation

A3 = 0.37*dsin(rad*Ma)  ! third correction


! -----------------------------------------------------------

Mms = Mm + Ev - Ae -A3   ! Moon’s corrected anomaly

Ec = 6.2886*dsin(rad*Mms)  ! equation of the centre

A4 = 0.214d0*dsin(rad*2*Mms) ! another correction term

ls = l + Ev + Ec - Ae + A4 ! Moon’s corrected longitude

V = 0.6583d0*dsin(rad*2*(ls-lama)) ! variation

lss = ls + v  ! Moon’s true orbital longitude


! ---------------------------------------------------------

Ns = N - 0.16d0*dsin(rad*Ma) ! corrected longitude of the node
i = 5.145396d0    ! inclination of Moon’s orbit

y = dsin(rad*(lss-Ns))*cos(rad*i)
x = dcos(rad*(lss-Ns)) + 1.e-8


at = datan(y/x)*1.d0/rad

call conv(at,x,y)


lamm = at + Ns  ! ecliptic longitude

if(lamm.gt.360.) lamm = lamm-360.d0


betm = dasin( dsin(rad*(lss-Ns))*dsin(rad*i))*1.d0/rad  ! ecliptic latitude


! -------------------------------------------------------------

T = (d + 3651.5d0)/36525.0d0 ! number of Julian centuries since epoch 2000 January 1.5

eps = 23.43929166663333336d0 - 0.013004166666666666d0*T - 1.6666666666666665d-7*T**2 + 5.027777777777778d-7*T**3 ! mean obliquity of the ecliptic

y = dsin(rad*lamm)*dcos(rad*eps) - dtan(rad*betm)*dsin(rad*eps)
x = dcos(rad*lamm) + 1.d-8


alp = datan(y/x)*1.d0/rad

call conv(alp,x,y)

alp = alp*0.06666666666666667d0 ! right ascension (hours)

call adj2(alp)

dcl = dasin(dsin(rad*betm)*dcos(rad*eps) + dcos(rad*betm)*dsin(rad*eps)*dsin(rad*lamm))*1.d0/rad ! declination

h = gst - alp  ! hour angle (hours)

call adj2(h)


return
end subroutine moon_declination

! ============================================================
subroutine conv(at,x,y)
implicit none

double precision at,x,y

if(at.lt.0.) at = at + 360.d0

if((x.ge.0).and.(y.ge.0.)) then
  if(at.gt.90.) at = at - 180.d0
endif

if((x.lt.0).and.(y.ge.0.)) then
  if(at.gt.270.) at = at - 180.d0
endif

if((x.lt.0).and.(y.lt.0.)) then
  if(at.lt.90.) at = at + 180.d0
endif

if((x.ge.0).and.(y.lt.0.)) then
  if(at.lt.180.) at = at + 180.d0
endif

return
end subroutine conv

! ==========================================================================
! declination (deg) and right ascension (hours) of the Sun after Duffet and Zwart (1979)

subroutine sun_declination(d,gst,dcl,alp,h)
implicit none

double precision d,dcl,alp,gst,h
double precision eps,t,ma,pi,rad,lam,n,ec,x,y,bet,azi

pi = 3.141592653589793d0
rad = pi/180.d0


n = 360.d0/365.242191d0*d  

call adj(n)

Ma = n + 279.557208d0-283.112438d0  !mean anomaly

call adj(ma)

ec = 360.d0/pi*0.016705d0*dsin(rad*Ma)

lam = n + ec + 279.557208  ! longitude of the Sun
call adj(lam)


! ===================================================================

bet = 0.d0

T = (d + 3651.5d0)/36525.0d0 ! number of Julian centuries since epoch 2000 January 1.5

eps = 23.43929166663333336d0 - 0.013004166666666666d0*T - 1.6666666666666665d-7*T**2 +5.027777777777778d-7*T**3 ! mean obliquity of the ecliptic

y = dsin(rad*lam)*dcos(rad*eps) - dtan(rad*bet)*dsin(rad*eps)
x = dcos(rad*lam) + 1.d-8


alp = datan(y/x)*1.d0/rad

call conv(alp,x,y)

alp = alp*0.06666666666666667d0 ! right ascension (hours)

dcl = dasin(dsin(rad*bet)*dcos(rad*eps) + dcos(rad*bet)*dsin(rad*eps)*dsin(rad*lam))*1.d0/rad ! declination

h = gst - alp

call adj2(h)



return
end subroutine sun_declination

! =======================================================================
! Umformung zu postivem Wert zwischen 0 und 360.

subroutine adj(x)
implicit none

double precision x

1 continue
if(x.lt.0) then
  x = x + 360.d0
  goto 1
endif

2 continue
if(x.gt.360.) then
  x = x - 360.d0
  goto 2
endif

return
end subroutine  adj

! =======================================================================
! Umformung zu positivem Wert zwischen 0. und 24.

subroutine adj2(x)
implicit none

double precision x

1 continue
if(x.lt.0) then
  x = x + 24.d0
  goto 1
endif

2 continue
if(x.gt.24.) then
  x = x - 24.d0
  goto 2
endif

return
end subroutine adj2

! =========================================================================
subroutine get_gst(year,mon,day,hour,minute,sekunde,gst) ! computation of Greenwich siderial time
implicit none

integer year,mon,day,hour,minute
double precision sekunde,gst,d,s,t0,t,ut

call timing(year,mon,day,0,0,0.d0,d)

s = d - 7193.5d0 + 10845.d0
t = s/36525.0d0
t0 = 6.697374558d0 + 2400.051336d0*t + 0.000025862d0*t**2

call adj2(t0)

ut = hour*1.d0 + minute*1.d0/60 + sekunde*1.d0/3600 

ut = ut*1.002737909d0

gst = t0 + ut

call adj2(gst)

return
end subroutine get_gst

! ==========================================================================

subroutine Moon_pot(dcl,h,lat,lon,r,pot) ! computation of the Moon's tidal potential
implicit none
     
double precision pot
double precision dcl,h

double precision pi,rad,gam,m,r,pot0,lat,lon,er,tet,phi,codec,latd,lond
double precision cosg,wrz


pi = 3.141592653589793d0
rad = pi/180.d0


gam = 6.67408d-11 ! gravitational constant [m^3/kg/s]
m = 7.349d22 ! Moon's mass [kg]
er = 6371.d3 ! earth radius [m]

pot0 = gam*m/r

latd = lat/rad
lond = lon/rad
  
tet = (90.d0 - latd)*rad
codec = (90.d0 - dcl)*rad
phi = (h*360.d0/24.d0 + lond)*rad

  
cosg = dcos(tet)*dcos(codec) + dsin(tet)*dsin(codec)*dcos(phi)  
wrz = dsqrt(1.d0 - 2.d0*er/r*cosg + (er/r)**2)
  
pot = -pot0*(1.d0 + er/r*cosg - 1.d0/wrz)
  

return
end subroutine Moon_pot


! ==========================================================================

subroutine Sun_pot(dcl,h,lat,lon,r,pot) ! computation of the Sun's tidal potential
implicit none
     
double precision pot
double precision dcl,h

double precision pi,rad,gam,m,r,pot0,lat,lon,er,tet,phi,codec,latd,lond
double precision cosg,wrz


pi = 3.141592653589793d0
rad = pi/180.d0


gam = 6.67408d-11 ! gravitaional constant [m^3/kg/s]
m = 1.9884d30 ! Sun's mass [kg]
er = 6371.d3 ! earth radius [m]

pot0 = gam*m/r

latd = lat/rad
lond = lon/rad
  
tet = (90.d0 - latd)*rad
codec = (90.d0 - dcl)*rad
phi = (h*360.d0/24.d0 + lond)*rad

  
cosg = dcos(tet)*dcos(codec) + dsin(tet)*dsin(codec)*dcos(phi)  
wrz = dsqrt(1.d0 - 2.d0*er/r*cosg + (er/r)**2)
  
pot = -pot0*(1.d0 + er/r*cosg - 1.d0/wrz)
  

return
end subroutine Sun_pot

! ===============================================================================================================
subroutine dist_earth_sun_moon(y,mon,day,h,rs,rm) ! distance between rs = Earth - Sun and rm = Earth - Moon  (m)
implicit none

integer y,mon,day,ys,mons,a,b,c,d
real h
double precision jd,t,epsg,omg,e,ms,ee,ny,pi,rad,rs,r0
double precision lams,p0,l0,l,dd,cc,A3,Ae,Ev,Mm,Mms,ec,rm,n
integer ih,im


pi = 4.d0*datan(1.d0)
rad = pi/180.d0

r0 = 1.495985d11   ! semi-major axis (m)

if(mon.lt.3) then
  ys = y-1
  mons = mon+12
else
  ys = y
  mons = mon
endif

a = dint(0.01d0*ys)
b = 2 - a + dint(0.25d0*a)
c = dint(365.25d0*ys)
d = dint(30.6001d0*(mons +1))

jd = b + c + d + day + h/24.d0 + 1720994.5d0
t = (jd - 2415020.0d0)/36525.d0

epsg = 279.6966778d0 + 36000.76892d0*t + 0.0003025d0*t**2 ! ecliptic longitute (degrees)
call adj(epsg)

omg = 281.2208444d0 + 1.719175d0*t + 0.000452778d0*t**2 ! ecliptic longitude of perigee (degrees)
call adj(omg)

e = 0.01675104d0 - 0.0000418d0*t - 0.000000126d0*t**2 ! eccentricity of orbit

ms = epsg - omg
call adj(ms)

ms = pi/180.d0*ms

call findEE(e,ms,ee)

ny = 2.d0*datan( dsqrt((1.d0+e)/(1.d0-e))*dtan(0.5*ee))
ny = 180.d0/pi*ny
call adj(ny)

rs = r0 *(1.d0 - e**2)/(1.d0+e*dcos(ny*rad))

! =========================================================================
! distance Earth - Moon (m)

ih = int(h)
im = int((h-ih)*1./60 + 0.5)

call timing(y,mon,day,ih,im,0.d0,DD)  ! DD = number of days since 2010 January 0.0

a = 384401.d3 ! semi-major axis (m)
e = 0.054900 ! eccentricity

n = 360.d0/365.242191d0*DD
call adj(n)

ms = ms*180./pi  ! Sun’s mean anomaly (deg)
call adj(ms)

lams = ny + omg   ! Sun’s ecliptic longitude (deg)
call adj(lams)

l0 = 91.929336d0   ! Moon’s mean longitude at the epoch (deg)
p0 = 130.143076d0  ! mean longitude of the perigee at the epoch (deg)

l = 13.1763966d0*DD + l0  ! Moon's mean longitude (deg)
call adj(l)

cc = l - lams
call adj(cc)

A3 = 0.37d0*dsin(ms*rad)
Ae = 0.1858d0*dsin(ms*rad)
Mm = l - 0.111404d0*DD-p0  ! Moon’s mean anomaly (deg)
call adj(Mm)

Ev = 1.2739*dsin((2.d0*cc - mm)*rad) ! evection

Mms = Mm + Ev - Ae - A3  ! Moon’s corrected anomaly (deg)

Ec = 6.2886*dsin(Mms*rad) 

rm = a*(1.d0 - e**2)/(1.d0 + e*dcos((Mms + Ec)*rad))


return
end subroutine dist_earth_sun_moon

! ================================================

subroutine findee(e,ms,ee)
implicit none

integer it
double precision e,ms,ee,ee1,ee2,incr,pi,v1,v2,fkt

INTEGER ie, ie1,ie2,iincr

pi = 4.d0*datan(1.d0)
it = 0


ie1 = 1
ie2 = 10
iincr = 1

ee1  = 0.d0
ee2  = 2.d0*pi
incr = 0.2d0*pi

2 continue

it = it + 1
ee=ee1
! do ee = ee1,ee2,incr
do ie = ie1,ie2,iincr
  ee = ee + incr
  fkt = ee - e*dsin(ee)
  if(fkt.gt.ms) then
    v1 = ee-incr
    v2 = ee
    goto 1
  endif
enddo

1 continue
ee1 = v1
ee2 = v2
incr = 0.1d0*incr

if(it.lt.8) goto 2

ee = 0.5d0*(v1+v2)

return
end subroutine findee  
  

END MODULE mo_ocean_tides
!=============================================================================
