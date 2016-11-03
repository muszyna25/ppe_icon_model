SUBROUTINE lcariolle_o3_column(                    &
         & jcb,          jce,           NCX,       &
         & nlev,         o3_vmr,        vmr2molm2, &
         & ldown,        o3_column                 )
USE mo_lcariolle_kind,   ONLY: wp,wi 
IMPLICIT NONE
INTEGER(wi),INTENT(IN)  :: jcb,jce,NCX,nlev
REAL(wp),INTENT(IN)     :: o3_vmr(NCX,nlev)
REAL(wp),INTENT(IN)     :: vmr2molm2(NCX,nlev)
LOGICAL,INTENT(IN)      :: ldown
REAL(wp),INTENT(INOUT)  :: o3_column(NCX,nlev)

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
