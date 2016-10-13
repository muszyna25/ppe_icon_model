SUBROUTINE cariolle_do3dt(jcb,jce,NCX,nlev,lat_weight_linear,pres_weight_linear,avi,do3dt)
USE mo_cariolle_kind,  ONLY: wp,wi
USE mo_cariolle_types, ONLY: t_avi,pvi,nlatx,nlevx
IMPLICIT NONE
INTEGER(wi),INTENT(IN) :: jcb,jce,NCX,nlev
EXTERNAL                  lat_weight_linear,pres_weight_linear
TYPE(t_avi),INTENT(IN) :: avi
REAL(wp),INTENT(INOUT) :: do3dt(NCX,nlev) ! tendency of VMR per second

INTEGER(wi)            :: ilev,ic
INTEGER(wi)            :: lev_temp,it_temp
REAL(wp)               :: o3_column(NCX,nlev)
REAL(wp)               :: wgt1_lat(NCX),wgt2_lat(NCX), &
                        & wgt1_p(NCX,nlev),wgt2_p(NCX,nlev)
INTEGER(wi)            :: inmw1_lat(NCX),inmw2_lat(NCX), &
                        & iw1_p(NCX,nlev),iw2_p(NCX,nlev)
REAL(wp)               :: wgt1,wgt2
INTEGER(wi)            :: iw1,iw2
REAL(wp)               :: al1(NCX,nlev), al2(NCX,nlev), &
                        & al3(NCX,nlev), al4(NCX,nlev), &
                        & al5(NCX,nlev), al6(NCX,nlev), &
                        & al7(NCX,nlev), al8(NCX,nlev)

CALL cariolle_o3_column(jcb,jce,NCX,nlev,avi%o3_vmr,avi%vmr2molm2,avi%ldown,o3_column)
write(*,*) 'o3_column(1,:)=',o3_column(1,:)
CALL lat_weight_linear( jcb,                 jce,                   NCX,               &
                      & avi%cell_center_lat, nlatx,                 pvi%rlat,          &
                      & pvi%delta_lat,       pvi%l_lat_sn,          wgt1_lat,          &
                      & wgt2_lat,            inmw1_lat,             inmw2_lat          )
CALL pres_weight_linear(jcb,jce,NCX,nlev,avi%pres,nlevx,pvi%plev,wgt1_p,wgt2_p,iw1_p,iw2_p)
lev_temp=10
it_temp=1
DO ilev=1,nlev
  DO ic=jcb,jce
    wgt1=wgt1_lat(ic)
    wgt2=wgt2_lat(ic)
    iw1 =inmw1_lat(ic)
    iw2 =inmw2_lat(ic)
    al1(ic,ilev) = wgt1*pvi%a1(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a1(iw2,lev_temp,it_temp)
    al2(ic,ilev) = wgt1*pvi%a2(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a2(iw2,lev_temp,it_temp)
    al3(ic,ilev) = wgt1*pvi%a3(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a3(iw2,lev_temp,it_temp)
    al4(ic,ilev) = wgt1*pvi%a4(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a4(iw2,lev_temp,it_temp)
    al5(ic,ilev) = wgt1*pvi%a5(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a5(iw2,lev_temp,it_temp)
    al6(ic,ilev) = wgt1*pvi%a6(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a6(iw2,lev_temp,it_temp)
    al7(ic,ilev) = wgt1*pvi%a7(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a7(iw2,lev_temp,it_temp)
    al8(ic,ilev) = wgt1*pvi%a8(iw1,lev_temp,it_temp) + &
                 & wgt2*pvi%a8(iw2,lev_temp,it_temp)
  END DO
END DO

DO ic=jcb,jce
  DO ilev=1,nlev
    do3dt(ic,ilev) = al1(ic,ilev) + &
                   & al2(ic,ilev)*(avi%o3_vmr(ic,ilev)-al3(ic,ilev)) + &
                   & al4(ic,ilev)*(avi%tmprt(ic,ilev)-al5(ic,ilev)) + &
                   & al6(ic,ilev)*(o3_column(ic,ilev)-al7(ic,ilev)) + &
                   & al8(ic,ilev)*avi%o3_vmr(ic,ilev)
  END DO
END DO
END SUBROUTINE cariolle_do3dt
