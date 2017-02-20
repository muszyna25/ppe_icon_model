!>
!! @brief This subroutine contains external subprograms referring to
!! subprograms provided by icon for the Cariolle scheme.
!! documentation: cr2016_10_22_rjs
!!
!! @author Sebastian Rast, MPI-M
!!
!! @par Revision History
!!  Original version Sebastian Rast (2016)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_lcariolle_externals
USE mo_kind,               ONLY: wp
USE mo_read_interface,     ONLY: read_bcast_REAL_3D, read_1D,  &
                               & closeFile, openInputFile
USE mo_physical_constants, ONLY: avo
USE mo_exception,          ONLY: finish
IMPLICIT NONE
PRIVATE
PUBLIC    :: read_bcast_real_3d_wrap, read_bcast_real_1d_wrap, & 
           & closeFile_wrap, openInputFile_wrap, get_constants
CONTAINS
  SUBROUTINE read_bcast_real_3d_wrap( &
       & file_id, variable_name, &
       & n1,      n2,            &
       & n3,      a_temp         )
    ! read a 3d-field from a netcdf file using the read_bcast_REAL_3D subroutine
    ! of ICON. This subroutine reads a 3d-field and sends it to all processors
  INTEGER, INTENT(IN)               :: file_id !< of necdf file
  CHARACTER(LEN=*), INTENT(IN)      :: variable_name !< in netcdf file
  INTEGER, INTENT(IN)               :: n1,n2,n3 !< dimensions in netcdf file
  REAL(wp), INTENT(INOUT)           :: a_temp(n1,n2,n3)
  REAL(wp), POINTER                 :: return_pointer(:,:,:)
  CALL read_bcast_REAL_3D(file_id=file_id, variable_name=variable_name, &
                          return_pointer=return_pointer)
  a_temp=return_pointer
END SUBROUTINE read_bcast_real_3d_wrap
SUBROUTINE read_bcast_real_1d_wrap( &
     & file_id, variable_name, &
     & n1, a_temp              )
    ! read a 1d-field from a netcdf file using the read_1D subroutine
    ! of ICON. This subroutine reads a 1d-field and sends it to all processors
  INTEGER, INTENT(IN)               :: file_id !< of netcdf file
  CHARACTER(LEN=*), INTENT(IN)      :: variable_name !< in netcdf file
  INTEGER, INTENT(IN)               :: n1 !< dimension in netcdf file
  REAL(wp), INTENT(INOUT)           :: a_temp(n1)
  REAL(wp), POINTER                 :: return_pointer(:)
  CALL read_1D(file_id, variable_name, return_pointer=return_pointer)
  a_temp=return_pointer
END SUBROUTINE read_bcast_real_1d_wrap
SUBROUTINE closeFile_wrap(file_id)
  ! close a netcdf file. Since this subroutine is overloaded, we need a wrapper
  INTEGER, INTENT(IN)               :: file_id
  CALL closeFile(file_id)
END SUBROUTINE closeFile_wrap
INTEGER FUNCTION openInputFile_wrap(filename)
  ! open a netcdf file. Since this subroutine is overloaded, we need a wrapper
  CHARACTER(LEN=*), INTENT(IN)      :: filename
  openInputFile_wrap=openInputFile(filename)
END FUNCTION openInputFile_wrap
REAL(wp) FUNCTION get_constants(constant_name)
  ! This functions is meant to create the constants needed inside the
  ! Cariolle scheme. When callint this function, you should make sure that
  ! all values are the same as those of the host model ICON
  CHARACTER(LEN=*), INTENT(IN)      :: constant_name
  SELECT CASE (TRIM(constant_name)) 
    CASE ('avogadro')
      get_constants=avo
    CASE DEFAULT
      CALL finish('get_constants: mo_lcariolle_externals', &
         & 'no method to get constant '//TRIM(constant_name))
  END SELECT
END FUNCTION get_constants
END MODULE mo_lcariolle_externals
