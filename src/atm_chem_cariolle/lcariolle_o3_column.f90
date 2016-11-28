!>
!! @brief This subroutine calculates the overhead ozone column in mole/m^2
!! The subroutine can take any number of columns jcb,...,jce, but the columns must
!! comprise the full ozone column.
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
SUBROUTINE lcariolle_o3_column(                    &
         & jcb,          jce,           NCX,       &
         & nlev,         o3_vmr,        vmr2molm2, &
         & ldown,        o3_column                 )
USE mo_lcariolle_kind,   ONLY: wp,wi 
IMPLICIT NONE
INTEGER(wi),INTENT(IN)  :: &
     & jcb,jce,  & !< begin, end index of column
     & NCX,      & !< first dim of fields as in calling subprogram
     & nlev        !< number of levels in column
REAL(wp),INTENT(IN)     :: o3_vmr(NCX,nlev)    !< ozone VMR 
REAL(wp),INTENT(IN)     :: vmr2molm2(NCX,nlev) !< conversion factor from VMR to mole/m^2
LOGICAL,INTENT(IN)      :: ldown               !< .true. if layers are counted from top
                                               !< to bottom, .false. otherwise
REAL(wp),INTENT(INOUT)  :: o3_column(NCX,nlev) !< overhead ozone column (mole/m^2)

INTEGER(wi)             :: ii,ilev,incr

IF (ldown) THEN
  ilev=1
  incr=1
ELSE
  ilev=nlev
  incr=-1
END IF
o3_column(jcb:jce,ilev) = 0.5_wp * o3_vmr(jcb:jce,ilev) * vmr2molm2(jcb:jce,ilev)
DO ii=1,nlev-1
   ilev=ilev+incr
   o3_column(jcb:jce,ilev) = o3_column(jcb:jce,ilev-incr) +                   &
         & 0.5_wp * o3_vmr(jcb:jce,ilev-incr) * vmr2molm2(jcb:jce,ilev-incr) +&
         & 0.5_wp * o3_vmr(jcb:jce,ilev) * vmr2molm2(jcb:jce,ilev)
END DO

END SUBROUTINE lcariolle_o3_column
