!>
!!  Module determines kinds for different precisions.
!!
!!  Module determines kinds for different precisions
!!  Number model from which the SELECTED_*\\_KIND are requested: <br>
!!
!! @f{tabular}{{r@{\hspace*{3em}}c@{\hspace*{3em}}c}
!!                     &4 byte REAL     &8 byte REAL        \\\
!!        CRAY:        &-               &precision =   13   \\\
!!                     &                &exponent  = 2465   \\\
!!        IEEE:        &precision = 6   &precision =   15   \\\
!!                     &exponent  = 37  &exponent  =  307
!! @f}
!! \\medskip
!!
!!  Most likely this are the only possible models.
!!
!! @par Revision History
!!  Working precision and comments are added by Luis Kornblueh (2001)
!!  Modification by Luis Kornblueh (2010-02-16):
!!  Working precision selection is parameterized by values, which can be as
!!  well used for selection of MPI data types associated with the respective
!!  Fortran datatypes.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_kind

  IMPLICIT NONE

  PRIVATE

  !--------------------------------------------------------------------
  !
  ! Floating point section
  ! ----------------------
  !
  INTEGER, PARAMETER :: ps =   6
  INTEGER, PARAMETER :: rs =  37
  !
  INTEGER, PARAMETER :: pd =  12
  INTEGER, PARAMETER :: rd = 307
  !
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs) !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
  !
  INTEGER, PARAMETER :: wp = dp                        !< selected working precision
  !
#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif

#ifdef __MIXED_PRECISION_2
  INTEGER, PARAMETER :: vp2 = sp
#else
  INTEGER, PARAMETER :: vp2 = wp
#endif
  !
  ! Integer section
  ! ---------------
  !
  INTEGER, PARAMETER :: pi2 =  4
  INTEGER, PARAMETER :: pi4 =  9
  INTEGER, PARAMETER :: pi8 = 14  ! could be larger, but SX cannot do some operations otherwise
  !
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(pi2)   !< at least 2 byte integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)   !< at least 4 byte integer
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)   !< at least 8 byte integer
  !
  !
  ! The following variable is made available internally only. configure needs to detect
  ! the addressing mode and according to this mo_kind has to be updated by an preprocessor
  ! directive and #include <config.h>. This needs some changes.
  !
  INTEGER, PARAMETER :: wi = i4                       !< selected working precission
  !
  PUBLIC :: sp, dp, wp, vp, vp2, i4, i8, i2
  !
  !--------------------------------------------------------------------

END MODULE mo_kind

