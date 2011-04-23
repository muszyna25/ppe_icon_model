      SUBROUTINE npr_getindex(nlon, nlat, glat, gilat, gilon, indexs)
!
!@(#) Library gmtri: Global gridpoint model, triangular grid
!@(#) Module pr_getindex.f90, V1.26 from 3/1/04, extracted: 6/30/04
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
! 1.14       2002/01/16 Helmut P. Frank; A.Mueller
!  Introduce KIND-notation for all reals with real2kind.pl
!  Replace amax1 by MAX
! 1.17       2002/05/03 A. Mueller
!  Code cleanup
! 1.18       2004/08/04 P. Ripodas
!  Adapted to use within the ICON project
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!
!=======================================================================
!
!
!------------------------------------------------------------------------------
!
! Modules used:
!
!------------------------------------------------------------------------------

      USE mo_kind, ONLY : dp

!==============================================================================

      IMPLICIT NONE

      INTEGER nlon, nlat
      REAL(KIND=dp) :: glat(-1:nlat+2)
      REAL(KIND=dp) :: gilon, gilat
      INTEGER indexs(2)
      REAL(KIND=dp) :: pi, dlati,  dlon, dloni,  xlat,  xlon
      INTEGER ilat  ,ilon
!
!
      pi    = 2._dp*asin(1._dp)
      dlati = 1._dp/(glat(nlat/2+1) - glat(nlat/2))
      dlon  = 2._dp*pi/REAL(nlon,dp)
      dloni = 1._dp/dlon
!
!
          xlat   = gilat
          ilat   = INT((xlat - glat(1))*dlati + 1._dp)
          IF(xlat .ge. glat(ilat+1)) ilat = ilat + 1
!
          indexs(1) = ilat

          xlon   = gilon
          IF (xlon.lt.0._dp) xlon = xlon + 2._dp*pi
          ilon   = INT(xlon*dloni + 1._dp)
          IF (ilon.EQ.nlon+1) ilon = nlon
!
          indexs(2) = ilon
!
!=======================================================================
!
      RETURN
      END
