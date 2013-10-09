!> Module whose sole purpose is to provide the current SVN revision
!  number.
!
!  @author 10/2013 : F. Prill, DWD
MODULE mo_util_svn

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    SUBROUTINE private_util_print_svn_version() BIND(C,NAME='print_svn_version') 
    END SUBROUTINE private_util_print_svn_version
  END INTERFACE

  PUBLIC :: printSVNVersion

CONTAINS

  SUBROUTINE printSVNVersion()
    CALL private_util_print_svn_version
  END SUBROUTINE printSVNVersion

END MODULE mo_util_svn
