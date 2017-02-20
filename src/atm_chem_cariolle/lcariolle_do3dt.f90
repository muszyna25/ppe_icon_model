!>
!! @brief This subroutine calculates the ozone tendency according to equation (1)
!! of Cariolle et al.: Atmos. Chem. Phys. 7, 2183 (2007).
!! The subroutine can take any number of columns jcb,...,jce, but the columns must
!! comprise the full overhead ozone column.
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
SUBROUTINE lcariolle_do3dt(                                         &
         & jcb,                 jce,             NCX,               &
         & nlev,                time_ip,         lat_intp_li,       &
         & pres_intp_li,        avi,             do3dt              )
USE mo_lcariolle_kind,  ONLY: wp,wi
USE mo_lcariolle_types, ONLY: t_time_interpolation,t_avi,pvi,nlatx,nlevx
IMPLICIT NONE
INTEGER(wi),INTENT(IN) :: &
     & jcb,jce,   & !< begin, end index of column
     & NCX,       & !< first dim of fields as in calling subprogram
     & nlev         !< number of levels in column
TYPE(t_time_interpolation), INTENT(IN) :: time_ip   !< contains linear interpolation weights
EXTERNAL                  lat_intp_li, pres_intp_li !< subprograms for linear interpolation
                                                    !< wrt latitudes and pressure
TYPE(t_avi),INTENT(IN) :: avi                       !< derived type containing all variables
                                                    !< passed to submodel Cariolle
REAL(wp),INTENT(INOUT) :: do3dt(NCX,nlev)           !< tendency of ozone VMR per second

INTEGER(wi)            :: ilev,ic
INTEGER(wi)            :: lev_temp
REAL(wp)               :: o3_column(NCX,nlev)       !< overhead ozone column 
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

! calculate overhead ozone column for each model layer
CALL lcariolle_o3_column(                         &
   & jcb,           jce,           NCX,           &
   & nlev,          avi%o3_vmr,    avi%vmr2molm2, &
   & avi%ldown,     o3_column                     )
o3_column(jcb:jce,1:nlev)=o3_column(jcb:jce,1:nlev)*pvi%avogadro*1.e-4_wp
! calculate linear interpolation weights for latitude interpolation
CALL lat_intp_li(                                                   &
   & jcb,                 jce,                   NCX,               &
   & avi%cell_center_lat, nlatx,                 pvi%rlat,          &
   & pvi%delta_lat,       pvi%l_lat_sn,          wgt1_lat,          &
   & wgt2_lat,            inmw1_lat,             inmw2_lat          )
! calculate linear interpolation weights for pressure interpolation 
CALL pres_intp_li(                                                  &
   & jcb,                 jce,                   NCX,               &
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
! latitude and pressure interpolation of at1,...,at8
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
! calculate equation (1) of Cariolle et al., Atmos. Chem. Phys. 7, 2183 (2007).
! first step: all terms except the A_8 term for polar stratospheric clouds
DO ilev=1,nlev
  DO ic=jcb,jce
    do3dt(ic,ilev) = al1(ic,ilev) + &
                   & al2(ic,ilev)*(avi%o3_vmr(ic,ilev)-al3(ic,ilev)) + &
                   & al4(ic,ilev)*(avi%tmprt(ic,ilev)-al5(ic,ilev)) + &
                   & al6(ic,ilev)*(o3_column(ic,ilev)-al7(ic,ilev))
  END DO
END DO
! Add A_8 in case of polar stratospheric clouds and daylight for
! additional ozone destruction
DO ilev=1,nlev
  WHERE (avi%lday(jcb:jce).AND.avi%tmprt(jcb:jce,ilev)<=195._wp)
    do3dt(jcb:jce,ilev) = do3dt(jcb:jce,ilev) + &
                        & al8(jcb:jce,ilev)*avi%o3_vmr(jcb:jce,ilev)
  END WHERE
END DO
END SUBROUTINE lcariolle_do3dt
