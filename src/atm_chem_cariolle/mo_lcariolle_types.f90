!>
!! @brief This module contains the definition of all derived types used
!! in the Cariolle scheme.
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
MODULE mo_lcariolle_types
USE mo_lcariolle_kind, ONLY: wp,wi
IMPLICIT NONE
PRIVATE
PUBLIC :: t_avi, t_pvi, nlatx, nlevx, nmonthx, pvi, avi, t_time_interpolation
! number of latitudes, pressure layers and months used in the Cariolle
! climatology
INTEGER(wi), PARAMETER          :: nlatx=65, nlevx=91, nmonthx=12
! The derived type t_avi is meant to create variables containing all
! variables that have to be filled in the host model and being passed to
! the Cariolle submodel.
TYPE t_avi
REAL(wp), POINTER :: tmprt(:,:)
REAL(wp), POINTER :: vmr2molm2(:,:) 
REAL(wp), POINTER :: pres(:,:)
REAL(wp), POINTER :: o3_vmr(:,:)
LOGICAL, POINTER  :: lday(:)
REAL(wp), POINTER :: cell_center_lat(:)
LOGICAL           :: ldown
LOGICAL           :: l_initialized_o3=.FALSE.
END type t_avi
TYPE(t_avi)       :: avi
! The derived type t_pvi is meant to host all variables connected to
! the Cariolle climatology of coefficients determining the ozone tendencies.
TYPE t_pvi
REAL(wp), DIMENSION(0:nlatx+1,nlevx,0:nmonthx+1) :: a1, a2, a3, &
                                                  & a4, a5, a6, &
                                                  & a7, a8
REAL(wp)                                         :: plev(nlevx),   &
                                                  & rlat(0:nlatx+1)
REAL(wp)          :: delta_lat
LOGICAL           :: l_lat_sn
REAL(wp)          :: avogadro
END type t_pvi
TYPE(t_pvi)       :: pvi
! The variables of type t_time_interpolation contain the weights and indices
! for a linear time interpolation of the monthly given climatology of
! coefficients of the Cariolle scheme.
TYPE t_time_interpolation
INTEGER(wi)          :: imonth1,imonth2
REAL(wp)             :: weight1,weight2
END TYPE t_time_interpolation
END MODULE mo_lcariolle_types
