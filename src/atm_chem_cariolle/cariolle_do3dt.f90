SUBROUTINE cariolle_do3dt(jcb,                 jce,             NCX,               &
                         &nlev,                time_ip,         lat_weight_linear, &
                         &pres_weight_linear,  avi,             do3dt              )
USE mo_cariolle_kind,  ONLY: wp,wi
USE mo_cariolle_types, ONLY: t_time_interpolation,t_avi,pvi,nlatx,nlevx
IMPLICIT NONE
INTEGER(wi),INTENT(IN) :: jcb,jce,NCX,nlev
TYPE(t_time_interpolation), INTENT(IN) :: time_ip
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
REAL(wp)               :: wgt1,wgt2,wp1,wp2
INTEGER(wi)            :: iw1,iw2,ip1,ip2
REAL(wp)               :: a1_p1,a1_p2,a2_p1,a2_p2,a3_p1,a3_p2,a4_p1,a4_p2, &
                        & a5_p1,a5_p2,a6_p1,a6_p2,a7_p1,a7_p2,a8_p1,a8_p2
REAL(wp)               :: al1(NCX,nlev), al2(NCX,nlev), &
                        & al3(NCX,nlev), al4(NCX,nlev), &
                        & al5(NCX,nlev), al6(NCX,nlev), &
                        & al7(NCX,nlev), al8(NCX,nlev)
REAL(wp)               :: at1(0:nlatx+1,nlevx), at2(0:nlatx+1,nlevx), &
                        & at3(0:nlatx+1,nlevx), at4(0:nlatx+1,nlevx), &
                        & at5(0:nlatx+1,nlevx), at6(0:nlatx+1,nlevx), &
                        & at7(0:nlatx+1,nlevx), at8(0:nlatx+1,nlevx)

CALL cariolle_o3_column(jcb,jce,NCX,nlev,avi%o3_vmr,avi%vmr2molm2,avi%ldown,o3_column)
CALL lat_weight_linear( jcb,                 jce,                   NCX,               &
                      & avi%cell_center_lat, nlatx,                 pvi%rlat,          &
                      & pvi%delta_lat,       pvi%l_lat_sn,          wgt1_lat,          &
                      & wgt2_lat,            inmw1_lat,             inmw2_lat          )
CALL pres_weight_linear(jcb,                 jce,                   NCX,               &
                      & nlev,                avi%pres,              nlevx,             &
                      & pvi%plev,            wgt1_p,                wgt2_p,            &
                      & iw1_p,               iw2_p                                     )
! interpolate coefficients with respect to time
wp1=time_ip%weight1
wp2=time_ip%weight2
ip1=time_ip%imonth1
ip2=time_ip%imonth2
at1(:,:)=wp1*pvi%a1(:,:,ip1)+wp2*pvi%a1(:,:,ip2)
at2(:,:)=wp1*pvi%a2(:,:,ip1)+wp2*pvi%a2(:,:,ip2)
at3(:,:)=wp1*pvi%a3(:,:,ip1)+wp2*pvi%a3(:,:,ip2)
at4(:,:)=wp1*pvi%a4(:,:,ip1)+wp2*pvi%a4(:,:,ip2)
at5(:,:)=wp1*pvi%a5(:,:,ip1)+wp2*pvi%a5(:,:,ip2)
at6(:,:)=wp1*pvi%a6(:,:,ip1)+wp2*pvi%a6(:,:,ip2)
at7(:,:)=wp1*pvi%a7(:,:,ip1)+wp2*pvi%a7(:,:,ip2)
at8(:,:)=wp1*pvi%a8(:,:,ip1)+wp2*pvi%a8(:,:,ip2)
it_temp=1
DO ilev=1,nlev
  DO ic=jcb,jce
    wgt1=wgt1_lat(ic)
    wgt2=wgt2_lat(ic)
    iw1 =inmw1_lat(ic)
    iw2 =inmw2_lat(ic)
    wp1 =wgt1_p(ic,ilev)
    wp2 =wgt2_p(ic,ilev)
    ip1 =iw1_p(ic,ilev)
    ip2 =iw2_p(ic,ilev)
! latitude interpolation for pressure level 1 (iw1_p)
    a1_p1 = wgt1*at1(iw1,ip1) + &
          & wgt2*at1(iw2,ip1)
    a2_p1 = wgt1*at2(iw1,ip1) + &
          & wgt2*at2(iw2,ip1)
    a3_p1 = wgt1*at3(iw1,ip1) + &
          & wgt2*at3(iw2,ip1)
    a4_p1 = wgt1*at4(iw1,ip1) + &
          & wgt2*at4(iw2,ip1)
    a5_p1 = wgt1*at5(iw1,ip1) + &
          & wgt2*at5(iw2,ip1)
    a6_p1 = wgt1*at6(iw1,ip1) + &
          & wgt2*at6(iw2,ip1)
    a7_p1 = wgt1*at7(iw1,ip1) + &
          & wgt2*at7(iw2,ip1)
    a8_p1 = wgt1*at8(iw1,ip1) + &
          & wgt2*at8(iw2,ip1)
! latitude interpolation for pressure level 2 (iw2_p)
    a1_p2 = wgt1*at1(iw1,ip2) + &
          & wgt2*at1(iw2,ip2)
    a2_p2 = wgt1*at2(iw1,ip2) + &
          & wgt2*at2(iw2,ip2)
    a3_p2 = wgt1*at3(iw1,ip2) + &
          & wgt2*at3(iw2,ip2)
    a4_p2 = wgt1*at4(iw1,ip2) + &
          & wgt2*at4(iw2,ip2)
    a5_p2 = wgt1*at5(iw1,ip2) + &
          & wgt2*at5(iw2,ip2)
    a6_p2 = wgt1*at6(iw1,ip2) + &
          & wgt2*at6(iw2,ip2)
    a7_p2 = wgt1*at7(iw1,ip2) + &
          & wgt2*at7(iw2,ip2)
    a8_p2 = wgt1*at8(iw1,ip2) + &
          & wgt2*at8(iw2,ip2)
! pressure level interpolation
    al1(ic,ilev)=wp1*a1_p1+wp2*a1_p2
    al2(ic,ilev)=wp1*a2_p1+wp2*a2_p2
    al3(ic,ilev)=wp1*a3_p1+wp2*a3_p2
    al4(ic,ilev)=wp1*a4_p1+wp2*a4_p2
    al5(ic,ilev)=wp1*a5_p1+wp2*a5_p2
    al6(ic,ilev)=wp1*a6_p1+wp2*a6_p2
    al7(ic,ilev)=wp1*a7_p1+wp2*a7_p2
    al8(ic,ilev)=wp1*a8_p1+wp2*a8_p2
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
!!$write(*,*) 'do3dt(1,5)=',do3dt(1,5)
END SUBROUTINE cariolle_do3dt
