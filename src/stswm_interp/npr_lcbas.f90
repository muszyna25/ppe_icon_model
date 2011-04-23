      SUBROUTINE npr_lcbas(grd,bas1,bas2)
!
!@(#) Library gmtri: Global gridpoint model, triangular grid
!@(#) Module pr_lcbas.f90, V1.26 from 3/1/04, extracted: 6/30/04
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
! 1.1        2000/01/28 J. Baumgardner
!  Initial Release
! 1.14       2002/01/16 Helmut P. Frank
!  Introduce KIND-notation for all reals with real2kind.pl
! 1.15       2004/08/04 P.Ripodas
!  adapted to use within the ICON project
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!
!=======================================================================
!
! Evaluate the partial Lagrangian cubic basis functions (denominator
! only ) for the grid points and gather grid values
!
!------------------------------------------------------------------------------
!
! Modules used:
!
!------------------------------------------------------------------------------

      USE mo_kind, ONLY : dp

!==============================================================================

      IMPLICIT NONE
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
      REAL(KIND=dp) :: grd(4)               ! grid stencil
!
! Output arguments
!
      REAL(KIND=dp) :: bas1(4),            &! grid values on stencil
     &     bas2(4)              ! lagrangian basis functions
!
!---------------------------Local variables-----------------------------
!
      REAL(KIND=dp) :: x0mx1,              &! |
     &     x0mx2,              &! |
     &     x0mx3,              &! |- grid value differences used in
     &     x1mx2,              &! |  weights
     &     x1mx3,              &! |
     &     x2mx3                ! |
!
!-----------------------------------------------------------------------
!
      x0mx1   = grd(1) - grd(2)
      x0mx2   = grd(1) - grd(3)
      x0mx3   = grd(1) - grd(4)
      x1mx2   = grd(2) - grd(3)
      x1mx3   = grd(2) - grd(4)
      x2mx3   = grd(3) - grd(4)
!
      bas1(1) = grd(1)
      bas1(2) = grd(2)
      bas1(3) = grd(3)
      bas1(4) = grd(4)
!
      bas2(1) =  1._dp/ ( x0mx1 * x0mx2 * x0mx3 )
      bas2(2) = -1._dp/ ( x0mx1 * x1mx2 * x1mx3 )
      bas2(3) =  1._dp/ ( x0mx2 * x1mx2 * x2mx3 )
      bas2(4) = -1._dp/ ( x0mx3 * x1mx3 * x2mx3 )
!
      END
