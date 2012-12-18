!>
!! This module contains debugging utilities, especially subroutines
!! for writing REAL arrays to NetCDF files (for debugging purposes).
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation,            F. Prill, DWD (2011-11-14)
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
!!
MODULE mo_util_debug

! enable the following directive for disabling
! NetCDF dumps at compile time:
!define DISABLE_DUMP 1

  !
  ! debugging utilities
  !
  USE mo_kind,                  ONLY: wp
  USE mo_util_string,           ONLY: int2string
  USE mo_util_netcdf,           ONLY: nf
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH
  USE mo_exception,             ONLY: finish

  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: dump_array_to_netcdf
  PUBLIC :: debug_step
  PUBLIC :: ldebug_enable

  INCLUDE 'netcdf.inc'

  ! The global variables "debug_step", "ldebug_enable" are useful when debugging output is
  ! desired only for certain steps inside a loop. For example, one may set "debug_step" to 
  ! the current iteration and use this value at some other place (where the
  ! original counter is not available).

  INTEGER :: debug_step    = 0       !< global counter
  LOGICAL :: ldebug_enable = .TRUE.  !< enabling/disabling dumps (during runtime)

  INTERFACE dump_array_to_netcdf
    MODULE PROCEDURE dump_array_to_netcdf_1d
    MODULE PROCEDURE dump_array_to_netcdf_2d
    MODULE PROCEDURE dump_array_to_netcdf_3d
  END INTERFACE
  !
CONTAINS

  !> Dumps real array to NetCDF file
  !
  !  Note: When using this SR, take care of OpenMP regions!

  SUBROUTINE dump_array_to_netcdf_1d(zfilename, p_array)
    CHARACTER (len=*), INTENT(IN) :: zfilename
    REAL(wp)         , INTENT(IN) :: p_array(:)
    ! local variables
    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: idim, ncfile, ncid_var, &
         &     ncid_dim(ndims), icount(ndims)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_util_debug:dump_array_to_netcdf_1d'


    IF (.NOT. ldebug_enable) RETURN

#ifndef DISABLE_DUMP
    WRITE (*,*) "Dumping ", zfilename   
    ! create NetCDF file:
    CALL nf(nf_create("00_"//TRIM(zfilename)//"_"//TRIM(int2string(debug_step))//".nc", &
      &               nf_clobber, ncfile), routine)
    ! create dimensions:
    DO idim=1,ndims
      CALL nf(nf_def_dim(ncfile, 'dim'//int2string(idim), SIZE(p_array,idim), &
        &     ncid_dim(idim)), routine)
      icount(idim) = SIZE(p_array,idim)
    END DO
    ! create variable:
    CALL nf(nf_def_var(ncfile, "var", NF_DOUBLE, ndims, ncid_dim(:), ncid_var), routine)
    ! End of definition mode
    CALL nf(nf_enddef(ncfile), routine)
    ! put data:
    CALL nf(nf_put_vara_double(ncfile, ncid_var, (/ (1, idim=1,ndims) /), &
      &                        icount, p_array), routine)
    ! close file
    CALL nf(nf_close(ncfile), routine)
#endif

  END SUBROUTINE dump_array_to_netcdf_1d


  !> Dumps real array to NetCDF file
  !
  !  Note: When using this SR, take care of OpenMP regions!

  SUBROUTINE dump_array_to_netcdf_2d(zfilename, p_array)
    CHARACTER (len=*), INTENT(IN) :: zfilename
    REAL(wp)         , INTENT(IN) :: p_array(:,:)
    ! local variables
    INTEGER, PARAMETER :: ndims = 2
    INTEGER :: idim, ncfile, ncid_var, &
         &     ncid_dim(ndims), icount(ndims)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_util_debug:dump_array_to_netcdf_2d'


    IF (.NOT. ldebug_enable) RETURN

#ifndef DISABLE_DUMP
    WRITE (*,*) "Dumping ", zfilename
    ! create NetCDF file:
    CALL nf(nf_create("00_"//TRIM(zfilename)//"_"//TRIM(int2string(debug_step))//".nc", &
      &               nf_clobber, ncfile), routine)
    ! create dimensions:
    DO idim=1,ndims
      CALL nf(nf_def_dim(ncfile, 'dim'//int2string(idim), SIZE(p_array,idim), &
        &     ncid_dim(idim)), routine)
      icount(idim) = SIZE(p_array,idim)
    END DO
    ! create variable:
    CALL nf(nf_def_var(ncfile, "var", NF_DOUBLE, ndims, ncid_dim(:), ncid_var), &
      &     routine)
    ! End of definition mode
    CALL nf(nf_enddef(ncfile), routine)
    ! put data:
    CALL nf(nf_put_vara_double(ncfile, ncid_var, (/ (1, idim=1,ndims) /), &
      &                        icount, p_array), routine)
    ! close file
    CALL nf(nf_close(ncfile), routine)
#endif

  END SUBROUTINE dump_array_to_netcdf_2d


  !> Dumps real array to NetCDF file
  !
  !  Note: When using this SR, take care of OpenMP regions!

  SUBROUTINE dump_array_to_netcdf_3d(zfilename, p_array)
    CHARACTER (len=*), INTENT(IN) :: zfilename
    REAL(wp)         , INTENT(IN) :: p_array(:,:,:)
    ! local variables
    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: idim, ncfile, ncid_var, &
         &     ncid_dim(ndims), icount(ndims)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_util_debug:dump_array_to_netcdf_3d'


    IF (.NOT. ldebug_enable) RETURN

#ifndef DISABLE_DUMP
    WRITE (*,*) "Dumping ", zfilename   
    ! create NetCDF file:
    CALL nf(nf_create("00_"//TRIM(zfilename)//"_"//TRIM(int2string(debug_step))//".nc", &
      &               nf_clobber, ncfile), routine)
    ! create dimensions:
    DO idim=1,ndims
      CALL nf(nf_def_dim(ncfile, 'dim'//int2string(idim), SIZE(p_array,idim), &
        &     ncid_dim(idim)), routine)
      icount(idim) = SIZE(p_array,idim)
    END DO
    ! create variable:
    CALL nf(nf_def_var(ncfile, "var", NF_DOUBLE, ndims, ncid_dim(:), ncid_var), &
      &     routine)
    ! End of definition mode
    CALL nf(nf_enddef(ncfile), routine)
    ! put data:
    CALL nf(nf_put_vara_double(ncfile, ncid_var, (/ (1, idim=1,ndims) /), &
      &                        icount, p_array), routine)
    ! close file
    CALL nf(nf_close(ncfile), routine)
#endif

  END SUBROUTINE dump_array_to_netcdf_3d

END MODULE mo_util_debug
