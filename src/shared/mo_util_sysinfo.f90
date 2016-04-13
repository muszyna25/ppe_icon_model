!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_sysinfo

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    SUBROUTINE private_util_user_name(name, name_len) BIND(C,NAME='util_user_name') 
#if defined (__SUNPRO_F95) 
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(out) :: name
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
#endif
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_user_name
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_os_system(name, name_len) BIND(C,NAME='util_os_system') 
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(out) :: name
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
#endif
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_os_system
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_node_name(name, name_len) BIND(C,NAME='util_node_name') 
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(out) :: name
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
#endif
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_node_name
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_get_maxrss(maxrss) BIND(C,NAME='util_get_maxrss') 
#if defined (__SUNPRO_F95) 
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
#else
      IMPORT :: C_INT
#endif
      INTEGER(C_INT), INTENT(out) :: maxrss
    END SUBROUTINE private_util_get_maxrss
  END INTERFACE


  INTERFACE
    SUBROUTINE private_util_compiler_release(release_str, rstr_len) BIND(C,NAME='util_compiler_release') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: release_str
      INTEGER(C_INT), INTENT(out) :: rstr_len
    END SUBROUTINE private_util_compiler_release
  END INTERFACE


  PUBLIC :: util_user_name
  PUBLIC :: util_os_system
  PUBLIC :: util_node_name
  PUBLIC :: util_get_maxrss
  PUBLIC :: util_compiler_release

CONTAINS

  SUBROUTINE util_user_name(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_user_name(name, name_len)
  END SUBROUTINE util_user_name

  SUBROUTINE util_os_system(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_os_system(name, name_len)
  END SUBROUTINE util_os_system

  SUBROUTINE util_node_name(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_node_name(name, name_len)
  END SUBROUTINE util_node_name

  SUBROUTINE util_get_maxrss(maxrss)
    INTEGER, INTENT(out) :: maxrss
    CALL private_util_get_maxrss(maxrss)
  END SUBROUTINE util_get_maxrss

  SUBROUTINE util_compiler_release(release_str, rstr_len)
    CHARACTER(len=*), INTENT(out) :: release_str
    INTEGER, INTENT(out) :: rstr_len
    CALL private_util_compiler_release(release_str, rstr_len)
  END SUBROUTINE util_compiler_release

END MODULE mo_util_sysinfo
