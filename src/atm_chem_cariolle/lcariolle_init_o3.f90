SUBROUTINE lcariolle_init_o3(                                   &
         & jcb,                 jce,         NCX,               &
         & nlev,                time_ip,     lat_intp_li,       &
         & pres_intp_li,        avi,         vmr_o3             )
USE mo_lcariolle_kind,  ONLY: wp,wi
USE mo_lcariolle_types, ONLY: t_time_interpolation,t_avi,pvi,nlatx,nlevx
IMPLICIT NONE
INTEGER(wi),INTENT(IN) :: jcb,jce,NCX,nlev
TYPE(t_time_interpolation), INTENT(IN) :: time_ip
EXTERNAL                  lat_intp_li,pres_intp_li
TYPE(t_avi),INTENT(IN) :: avi
REAL(wp),INTENT(INOUT) :: vmr_o3(NCX,nlev) ! tendency of VMR per second

INTEGER(wi)            :: ilev,ic
REAL(wp)               :: wgt1_lat(NCX),wgt2_lat(NCX), &
                        & wgt1_p(NCX,nlev),wgt2_p(NCX,nlev)
INTEGER(wi)            :: inmw1_lat(NCX),inmw2_lat(NCX), &
                        & iw1_p(NCX,nlev),iw2_p(NCX,nlev)
REAL(wp)               :: wgt1,wgt2,wp1,wp2
INTEGER(wi)            :: iw1,iw2,ip1,ip2
REAL(wp)               :: a3_p1,a3_p2
REAL(wp)               :: al3(NCX,nlev)
REAL(wp)               :: at3(0:nlatx+1,nlevx)

CALL lat_intp_li(                                           &
   & jcb,                 jce,                NCX,          &
   & avi%cell_center_lat, nlatx,              pvi%rlat,     &
   & pvi%delta_lat,       pvi%l_lat_sn,       wgt1_lat,     &
   & wgt2_lat,            inmw1_lat,          inmw2_lat     )
CALL pres_intp_li(                                          &
   & jcb,                 jce,                NCX,          &
   & nlev,                avi%pres,           nlevx,        &
   & pvi%plev,            wgt1_p,             wgt2_p,       &
   & iw1_p,               iw2_p                             )
! interpolate coefficients with respect to time
wp1=time_ip%weight1
wp2=time_ip%weight2
ip1=time_ip%imonth1
ip2=time_ip%imonth2
at3(:,:)=wp1*pvi%a3(:,:,ip1)+wp2*pvi%a3(:,:,ip2)
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
    a3_p1 = wgt1*at3(iw1,ip1) + &
          & wgt2*at3(iw2,ip1)
! latitude interpolation for pressure level 2 (iw2_p)
    a3_p2 = wgt1*at3(iw1,ip2) + &
          & wgt2*at3(iw2,ip2)
! pressure level interpolation
    vmr_o3(ic,ilev)=wp1*a3_p1+wp2*a3_p2
  END DO
END DO
END SUBROUTINE lcariolle_init_o3
