MODULE mo_cariolle_types
USE mo_cariolle_kind, ONLY: wp,wi
IMPLICIT NONE
PRIVATE
PUBLIC :: t_avi, t_pvi, nlatx, nlevx, nmonthx, pvi, avi
INTEGER(wi), PARAMETER          :: nlatx=65, nlevx=91, nmonthx=12
TYPE t_avi
REAL(wp), POINTER :: tmprt(:,:)
REAL(wp), POINTER :: vmr2molm2(:,:) 
REAL(wp), POINTER :: pres(:,:)
REAL(wp), POINTER :: o3_vmr(:,:)
REAL(wp), POINTER :: cell_center_lat(:)
LOGICAL           :: ldown
END type t_avi
TYPE(t_avi)       :: avi
TYPE t_pvi
REAL(wp), DIMENSION(0:nlatx+1,nlevx,nmonthx) :: a1, a2, a3, &
                                          & a4, a5, a6, &
                                          & a7, a8
REAL(wp)                                 :: plev(nlevx),   &
                                          & rlat(0:nlatx+1)
REAL(wp)          :: delta_lat
LOGICAL           :: l_lat_sn
END type t_pvi
TYPE(t_pvi)       :: pvi
TYPE t_lat_interpolation
INTEGER(wi), POINTER :: ilat1(:),ilat2(:)
REAL(wp), POINTER :: weight1(:),weight2(:)
END type t_lat_interpolation
TYPE t_time_interpolation
INTEGER(wi), POINTER :: imonth1,imonth2
REAL(wp), POINTER :: weight1,weight2
END TYPE t_time_interpolation

END MODULE mo_cariolle_types
