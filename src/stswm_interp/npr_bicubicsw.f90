      SUBROUTINE npr_bicubicsw(nlon,nlat,glat,wiy,sg1,sg2,gifile)
!
!@(#) Library gmtri: Global gridpoint model, triangular grid
!@(#) Module pr_bicubics.f90, V1.26 from 3/1/04, extracted: 6/30/04
!
!
!=======================================================================
!
!     This routine interpolates scalar field sg defined on a latitude-
!     longitude grid with nlon longitudes, nlev vertical levels, and
!     nlat latitudes to a scalar field s defined on a icosahedral mesh
!     with ni intervals along each icosahedral edge.  The array glat
!     contains the latitudes on the latitude-longitude grid.
!     The bicubic interpolation uses a 12-point stencil.
!     This routine is based on a program by J. Baumgardner
!
!=======================================================================
!
! Current Code Owner: DWD, D. Liermann
!    phone: +49-69-8062-2732, fax: +49-69-8062-3721
!    email: doerte.liermann@dwd.de
!
!=======================================================================
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2000/01/28 D. Liermann
!  Initial Release
! 1.14       2002/01/16 Helmut P. Frank
!  Introduce KIND-notation for all reals with real2kind.pl
! 1.15       2004/08/04 P. Ripodas
! Adapted to be use within the ICON project
!
!=======================================================================
!
!------------------------------------------------------------------------------
!
! Modules used:
!
!------------------------------------------------------------------------------

      USE mo_kind, ONLY : dp

!==============================================================================

      IMPLICIT NONE

!     INCLUDE "gme_param.h"
!      INCLUDE "gme_horgrid.h"
!
      INTEGER nlon, nlat
      INTEGER indexs(2)
      REAL(KIND=dp) :: sg1(nlon,nlat),sg2(nlon,nlat),sn1,sn2
      REAL(KIND=dp) :: glat(-1:nlat+2)
      REAL(KIND=dp) :: wiy(4,2,0:nlat)
!
      REAL(KIND=dp) :: dlon, dloni, sixth, s1, s2, s3, s4
      REAL(KIND=dp) :: xlat, xlon, x1, x2, x3, x4, y1, y2, y3, y4
      REAL(KIND=dp) :: pi
      INTEGER ilat, ilon, ilonm1
      INTEGER ilonp1, ilonp2, jlon, jlonm1, jlonp1, jlonp2, mlat
      CHARACTER(len=80) :: gifile
!
      pi    = 2._dp*asin(1._dp)
      dlon  = 2._dp*pi/REAL(nlon,dp)
      dloni = 1._dp/dlon
      sixth = 1._dp/6._dp
      mlat  = nlat + nlat - 1

!
      OPEN (12,FILE="outinterp.dat")
      OPEN (11,FILE=TRIM(gifile))
!
      reading_loop: DO
!
           READ (11,*,end=100) xlon, xlat
           !WRITE(6,*)xlon,xlat
           xlat   = xlat*pi/180._dp
           xlon   = xlon*pi/180._dp
           IF(xlon .lt. 0._dp) xlon = xlon + 2._dp*pi
!
           CALL npr_getindex(nlon, nlat, glat, xlat, xlon, indexs)
           ilat = indexs(1)
           ilon = indexs(2)
!
           ilonp1 = mod(ilon, nlon) + 1
           ilonp2 = mod(ilonp1, nlon) + 1
           ilonm1 = mod(ilon+nlon-2, nlon) + 1
           jlon   = mod(ilon+nlon/2-1, nlon) + 1
           jlonp1 = mod(jlon, nlon) + 1
           jlonp2 = mod(jlonp1, nlon) + 1
           jlonm1 = mod(jlon+nlon-2, nlon) + 1

!              Compute the x-interpolants for the four latitudes.

           x2 =  (xlon - dlon*(REAL(ilon,dp) - 1._dp))*dloni
           x1 = -x2 - 1._dp
           x3 = -x2 + 1._dp
           x4 =  x2 - 2._dp

           IF(ilat .le. 1) THEN
              s1 = sg1(jlonp1,2-ilat)*x2  &
     &           + sg1(jlon  ,2-ilat)*x3
           ELSE
              s1 = sg1(ilonp1,ilat-1)*x2  &
     &           + sg1(ilon  ,ilat-1)*x3
           ENDIF

           IF(ilat .EQ. 0) THEN
              s2 = sg1(jlonm1,1     )*x2*x3*x4*sixth         &
     &           + sg1(jlon  ,1     )*x1*x3*x4*0.5_dp    &
     &           + sg1(jlonp1,1     )*x1*x2*x4*0.5_dp    &
     &           + sg1(jlonp2,1     )*x1*x2*x3*sixth
           ELSE
              s2 = sg1(ilonm1,ilat  )*x2*x3*x4*sixth         &
     &           + sg1(ilon  ,ilat  )*x1*x3*x4*0.5_dp    &
     &           + sg1(ilonp1,ilat  )*x1*x2*x4*0.5_dp    &
     &           + sg1(ilonp2,ilat  )*x1*x2*x3*sixth
           ENDIF

           IF(ilat .EQ. nlat) THEN
              s3 = sg1(jlonm1,nlat  )*x2*x3*x4*sixth         &
     &           + sg1(jlon  ,nlat  )*x1*x3*x4*0.5_dp    &
     &           + sg1(jlonp1,nlat  )*x1*x2*x4*0.5_dp    &
     &           + sg1(jlonp2,nlat  )*x1*x2*x3*sixth
           ELSE
              s3 = sg1(ilonm1,ilat+1)*x2*x3*x4*sixth         &
     &           + sg1(ilon  ,ilat+1)*x1*x3*x4*0.5_dp    &
     &           + sg1(ilonp1,ilat+1)*x1*x2*x4*0.5_dp    &
     &           + sg1(ilonp2,ilat+1)*x1*x2*x3*sixth
           ENDIF

           IF(ilat .ge. nlat-1) THEN
              s4 = sg1(jlonp1,mlat-ilat)*x2  &
     &           + sg1(jlon  ,mlat-ilat)*x3
           ELSE
              s4 = sg1(ilonp1,ilat+2)*x2     &
     &           + sg1(ilon  ,ilat+2)*x3
           ENDIF

!          Compute the y-interpolant.

           y1 = xlat - wiy(1,1,ilat)
           y2 = xlat - wiy(2,1,ilat)
           y3 = xlat - wiy(3,1,ilat)
           y4 = xlat - wiy(4,1,ilat)

           sn1 = s1*y2*y3*y4*wiy(1,2,ilat)  &
     &                + s2*y1*y3*y4*wiy(2,2,ilat)      &
     &                + s3*y1*y2*y4*wiy(3,2,ilat)      &
     &                + s4*y1*y2*y3*wiy(4,2,ilat)

           IF(ilat .le. 1) THEN
              s1 = sg2(jlonp1,2-ilat)*x2  &
     &           + sg2(jlon  ,2-ilat)*x3
           ELSE
              s1 = sg2(ilonp1,ilat-1)*x2  &
     &           + sg2(ilon  ,ilat-1)*x3
           ENDIF

           IF(ilat .EQ. 0) THEN
              s2 = sg2(jlonm1,1     )*x2*x3*x4*sixth         &
     &           + sg2(jlon  ,1     )*x1*x3*x4*0.5_dp    &
     &           + sg2(jlonp1,1     )*x1*x2*x4*0.5_dp    &
     &           + sg2(jlonp2,1     )*x1*x2*x3*sixth
           ELSE
              s2 = sg2(ilonm1,ilat  )*x2*x3*x4*sixth         &
     &           + sg2(ilon  ,ilat  )*x1*x3*x4*0.5_dp    &
     &           + sg2(ilonp1,ilat  )*x1*x2*x4*0.5_dp    &
     &           + sg2(ilonp2,ilat  )*x1*x2*x3*sixth
           ENDIF

           IF(ilat .EQ. nlat) THEN
              s3 = sg2(jlonm1,nlat  )*x2*x3*x4*sixth         &
     &           + sg2(jlon  ,nlat  )*x1*x3*x4*0.5_dp    &
     &           + sg2(jlonp1,nlat  )*x1*x2*x4*0.5_dp    &
     &           + sg2(jlonp2,nlat  )*x1*x2*x3*sixth
           ELSE
              s3 = sg2(ilonm1,ilat+1)*x2*x3*x4*sixth         &
     &           + sg2(ilon  ,ilat+1)*x1*x3*x4*0.5_dp    &
     &           + sg2(ilonp1,ilat+1)*x1*x2*x4*0.5_dp    &
     &           + sg2(ilonp2,ilat+1)*x1*x2*x3*sixth
           ENDIF

           IF(ilat .ge. nlat-1) THEN
              s4 = sg2(jlonp1,mlat-ilat)*x2  &
     &           + sg2(jlon  ,mlat-ilat)*x3
           ELSE
              s4 = sg2(ilonp1,ilat+2)*x2     &
     &           + sg2(ilon  ,ilat+2)*x3
           ENDIF

!          Compute the y-interpolant.

           y1 = xlat - wiy(1,1,ilat)
           y2 = xlat - wiy(2,1,ilat)
           y3 = xlat - wiy(3,1,ilat)
           y4 = xlat - wiy(4,1,ilat)

           sn2 = s1*y2*y3*y4*wiy(1,2,ilat)  &
     &                + s2*y1*y3*y4*wiy(2,2,ilat)      &
     &                + s3*y1*y2*y4*wiy(3,2,ilat)      &
     &                + s4*y1*y2*y3*wiy(4,2,ilat)


           WRITE(12,*)xlon*180._dp/pi, xlat*180._dp/pi, sn1,sn2
           !WRITE(6,*)xlon*180._dp/pi, xlat*180._dp/pi, s

      ENDDO reading_loop
100   CLOSE(11)
      CLOSE(12)
!

      END
