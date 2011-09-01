!>
!! Provides interfaces for vendor accelarated math functions
!!
!! @author Leonidas Linardakis, MPI-M
!!
!!
!! @par Revision History
!! - Initial version by Leonidas Linardakis (2011-09-01).
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_fast_math_functions

  USE mo_io_units, ONLY: nerr, nlog, filename_max
  USE mo_mpi,      ONLY: run_is_global_mpi_parallel, p_abort, my_process_is_stdio, &
    & get_my_global_mpi_id
  USE mo_kind,     ONLY: wp

#ifdef __USE_MATH_LIB__
#ifdef  __xlC__
!  INCLUDE 'mass.include'
#endif
#endif

  IMPLICIT NONE

  PRIVATE

 ! PUBLIC :: cube_root_fc
  PUBLIC :: cube_root_rt

  INTERFACE cube_root_rt
    MODULE PROCEDURE vector_cube_root_rt_vsize
!     MODULE PROCEDURE vector_cube_root_rt
!    MODULE PROCEDURE scalar_cube_root_rt
  END INTERFACE
  
!  INTERFACE cube_root_fc
!    MODULE PROCEDURE scalar_cube_root_fc
!  END INTERFACE
 

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  
  !-------------------------------
  !>
  SUBROUTINE vector_cube_root_rt_vsize(in_vector, out_vector, vector_size)
    REAL(wp), INTENT(in)  :: in_vector(:)
    REAL(wp), INTENT(inout) :: out_vector(:)
    INTEGER, INTENT(in), OPTIONAL  :: vector_size

    INTEGER :: in_vector_size

    IF (PRESENT(vector_size)) THEN
      in_vector_size = vector_size
    ELSE
      in_vector_size = SIZE(in_vector)
    ENDIF
    IF (in_vector_size < 1) RETURN
    
#ifdef __USE_MATH_LIB__
#ifdef  __xlC__
     CALL vcbrt(out_vector,in_vector,in_vector_size)
#else
     ! no compiler specific optimization
     out_vector(1:in_vector_size) = in_vector(1:in_vector_size)**0.33333333_wp     
#endif

#else
     ! no __USE_MATH_LIB__ 
     out_vector(1:in_vector_size) = in_vector(1:in_vector_size)**0.33333333_wp     

#endif

  END SUBROUTINE vector_cube_root_rt_vsize
  !-------------------------------
  
!   !-------------------------------
!   !>
!   SUBROUTINE vector_cube_root_rt(in_vector, out_vector)
!     REAL(wp), INTENT(in)  :: in_vector(:)
!     REAL(wp), INTENT(out) :: out_vector(:)
! 
!     INTEGER :: in_vector_size
! 
!     in_vector_size = SIZE(in_vector)
!     
! #ifdef __USE_MATH_LIB__
! #ifdef  __xlC__
!      CALL vcbrt(out_vector,in_vector,in_vector_size)
! #else
!      ! no compiler specific optimization
!      out_vector(1:in_vector_size) = in_vector(1:in_vector_size)**0.33333333_wp     
! #endif
! 
! #else
!      ! no __USE_MATH_LIB__ 
!      out_vector(1:in_vector_size) = in_vector(1:in_vector_size)**0.33333333_wp     
! 
! #endif
! 
!   END SUBROUTINE vector_cube_root_rt
!   !-------------------------------
!   
  !-------------------------------
  !>
!  SUBROUTINE scalar_cube_root_rt(in_scalar, out_scalar)
!    REAL(wp), INTENT(in)  :: in_scalar
!    REAL(wp), INTENT(out) :: out_scalar
     
!#ifdef __USE_MATH_LIB__
!#ifdef  __xlC__
!     out_scalar = cbrt(in_scalar)
!#else
!     ! no compiler specific optimization
!     out_scalar = in_scalar**0.33333333_wp
!#endif
!
!#else
!     ! no __USE_MATH_LIB__ 
!     out_scalar = in_scalar**0.33333333_wp
!#endif
!    
!  END SUBROUTINE scalar_cube_root_rt
  !-------------

  
  !-------------------------------
  !>
!  FUNCTION scalar_cube_root_fc(in_scalar) result(out_scalar)
!    REAL(wp), INTENT(in)  :: in_scalar
!    REAL(wp) :: out_scalar
!     
!#ifdef __USE_MATH_LIB__
!#ifdef  __xlC__
!     out_scalar = cbrt(in_scalar)
!#else
!     ! no compiler specific optimization
!     out_scalar = in_scalar**0.33333333_wp
!#endif
!
!#else
!     ! no __USE_MATH_LIB__ 
!     out_scalar = in_scalar**0.33333333_wp
!#endif
!    
!  END FUNCTION scalar_cube_root_fc
  !-------------

END MODULE mo_fast_math_functions
