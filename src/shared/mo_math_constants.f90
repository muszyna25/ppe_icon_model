!>
!!  Module determines main mathematical constants to be used by shallow water.
!!
!!  Module determines main mathematical constants to be used by shallow water
!!  prototype. Values are taken from glibc 2.2.5: /usr/include/math.h
!!  These constants
!!  are provided to more significant digits than is necessary for a 64-bit
!!  double precision number; they may be used for other purposes where the
!!  extra precision is necessary or useful.
!!
!! @par Revision History
!!  Developed  by Luis Kornblueh (2004)
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Modified according to style guide by Thomas Heinze (2005-06-24):
!!   - module renamed from mo_math to mo_math_constants
!!   - eps moved from mo_physical_constants
!!   - pid180i renamed to rad2deg
!!  Including some more constants from math.h by Thomas Heinze (2005-07-18)
!!  pid5 renamed to pi_5 by Thomas Heinze (2005-07-26)
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Thomas Heinze (2006-05-18):
!!  - introduced dbl_eps
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
!! <li> In case you intend to Modified according to style guideuse the code commercially, we oblige you to sign
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
MODULE mo_math_constants

USE mo_kind,      ONLY:  wp

IMPLICIT NONE

PUBLIC

CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

!
! ! Mathematical constants (in brackets original C names)
!
! ! euler     (mo_E       )  -- e
! ! log2e     (mo_LOG2E   )  -- log2(e)
! ! log10e    (mo_LOG10E  )  -- log10(e)
! ! ln2       (mo_LN2     )  -- ln(2)
! ! ln10      (mo_LN10    )  -- ln(10)
! ! pi        (mo_PI      )  -- pi
! ! pi_2      (mo_PI_2    )  -- pi/2
! ! pi_4      (mo_PI_4    )  -- pi/4
! ! rpi       (mo_1_PI    )  -- 1/pi
! ! rpi_2     (mo_2_PI    )  -- 2/pi
! ! rsqrtpi_2 (mo_2_SQRTPI)  -- 2/(sqrt(pi))
! ! sqrt2     (mo_SQRT2   )  -- sqrt(2)
! ! sqrt1_2   (mo_SQRT1_2 )  -- 1/sqrt(2)
! ! sqrt3                    -- sqrt(3)
! ! sqrt1_3                  -- 1/sqrt(3)
! !
!
!

REAL (wp), PARAMETER ::  euler     = 2.71828182845904523536028747135266250_wp
REAL (wp), PARAMETER ::  log2e     = 1.44269504088896340735992468100189214_wp
REAL (wp), PARAMETER ::  log10e    = 0.434294481903251827651128918916605082_wp
REAL (wp), PARAMETER ::  ln2       = 0.693147180559945309417232121458176568_wp
REAL (wp), PARAMETER ::  ln10      = 2.30258509299404568401799145468436421_wp
REAL (wp), PARAMETER ::  pi        = 3.14159265358979323846264338327950288_wp
REAL (wp), PARAMETER ::  pi_2      = 1.57079632679489661923132169163975144_wp
REAL (wp), PARAMETER ::  pi_4      = 0.785398163397448309615660845819875721_wp
REAL (wp), PARAMETER ::  rpi       = 0.318309886183790671537767526745028724_wp
REAL (wp), PARAMETER ::  rpi_2     = 0.636619772367581343075535053490057448_wp
REAL (wp), PARAMETER ::  rsqrtpi_2 = 1.12837916709551257389615890312154517_wp
REAL (wp), PARAMETER ::  sqrt2     = 1.41421356237309504880168872420969808_wp
REAL (wp), PARAMETER ::  sqrt1_2   = 0.707106781186547524400844362104849039_wp
REAL (wp), PARAMETER ::  sqrt3     = 1.7320508075688772935274463415058723_wp
REAL (wp), PARAMETER ::  sqrt1_3   = 0.5773502691896257645091487805019575_wp

!
! ! some more useful constants
! ! pi_5    -- half angle of pentagon
! ! rad2deg -- conversion factor from radians to degree
! ! deg2rad -- conversion factor from degree to radians
! ! eps     -- residual bound for solvers
!

REAL (wp), PARAMETER ::  pi_5      = pi*0.2_wp
REAL (wp), PARAMETER ::  pi2       = pi*2.0_wp
REAL (wp), PARAMETER ::  rad2deg   = 180.0_wp/pi
REAL (wp), PARAMETER ::  deg2rad   = pi/180.0_wp
REAL (wp), PARAMETER ::  eps       = 1.e-8_wp
REAL (wp), PARAMETER ::  dbl_eps   = EPSILON(1._wp)

!
! ! phi0 is  the latitude of the lowest major triangle corner
! ! and the latitude of the major hexagonal faces centers
! ! phi0 = 0.5_wp*pi - 2._wp*acos(1.0_wp/(2._wp*sin(pi/5._wp)))
!

REAL (wp), PARAMETER ::  phi0  = 0.46364760900080614903_wp


!--------------------------------------------------------------------

END MODULE mo_math_constants

