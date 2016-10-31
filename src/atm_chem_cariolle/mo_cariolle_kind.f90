MODULE mo_cariolle_kind
  USE mo_kind, ONLY: dp, i4
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wp, wi
!!$  INTEGER, PARAMETER :: pd =  12
!!$  INTEGER, PARAMETER :: rd = 307
!!$  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)
!!$  INTEGER, PARAMETER :: pi8= 8
!!$  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)  
!!$  !
!!$  INTEGER, PARAMETER :: wp = dp
!!$  INTEGER, PARAMETER :: wi = i8

INTEGER, PARAMETER :: wp=dp
INTEGER, PARAMETER :: wi=i4
END MODULE mo_cariolle_kind
