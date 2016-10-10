SUBROUTINE cariolle_do3dt(jcb,jce,NCX,nlev,avi,do3dt)
USE mo_cariolle_kind,  ONLY: wp,wi
USE mo_cariolle_types, ONLY: t_avi,pvi
IMPLICIT NONE
INTEGER(wi),INTENT(IN) :: jcb,jce,NCX,nlev
TYPE(t_avi),INTENT(IN) :: avi
REAL(wp),INTENT(INOUT) :: do3dt(NCX,nlev) ! tendency of VMR per second

INTEGER(wi)            :: ilev,ic
INTEGER(wi)            :: lev_temp,it_temp
REAL(wp)               :: o3_column(NCX,nlev)
REAL(wp)               :: wgt1_lat(NCX),wgt2_lat(NCX)
INTEGER(wi)            :: inmw1_lat(NCX),inmw2_lat(NCX)
REAL(wp)               :: al1(NCX,nlev), al2(NCX,nlev), &
                        & al3(NCX,nlev), al4(NCX,nlev), &
                        & al5(NCX,nlev), al6(NCX,nlev), &
                        & al7(NCX,nlev), al8(NCX,nlev)

!!$write(*,*) 'avi%cell_center_lat=',avi%cell_center_lat(jcb:jce)
CALL cariolle_o3_column(jcb,jce,NCX,nlev,avi%o3_vmr,avi%vmr2molm2,avi%ldown,o3_column)
write(*,*) 'o3_column(1,:)=',o3_column(1,:)
!!$CALL lat_weight_li( jcb,                 jce,                   NCX,               &
!!$                  & avi%cell_center_lat, avi%l_cell_sn,         pvi%lat,           &
!!$                  & pvi%lat_shift,       pvi%delta_lat,         pvi%nlat,          &
!!$                  & pvi%l_lat_sn,        wgt1_lat,              wgt2_lat,          &
!!$                  & inmw1_lat,           inmw2_lat                                 )
!!$
!!$lev_temp=10
!!$it_temp=1
!!$DO ilev=1,nlev
!!$  DO ic=jcb,jce
!!$    al1(ic,ilev) = wgt1_lat(ic)*pvi%a1(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a1(inmw2_lat,lev_temp,it_temp)
!!$    al2(ic,ilev) = wgt1_lat(ic)*pvi%a2(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a2(inmw2_lat,lev_temp,it_temp)
!!$    al3(ic,ilev) = wgt1_lat(ic)*pvi%a3(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a3(inmw2_lat,lev_temp,it_temp)
!!$    al4(ic,ilev) = wgt1_lat(ic)*pvi%a4(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a4(inmw2_lat,lev_temp,it_temp)
!!$    al5(ic,ilev) = wgt1_lat(ic)*pvi%a5(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a5(inmw2_lat,lev_temp,it_temp)
!!$    al6(ic,ilev) = wgt1_lat(ic)*pvi%a6(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a6(inmw2_lat,lev_temp,it_temp)
!!$    al7(ic,ilev) = wgt1_lat(ic)*pvi%a7(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a7(inmw2_lat,lev_temp,it_temp)
!!$    al8(ic,ilev) = wgt1_lat(ic)*pvi%a8(inmw1_lat,lev_temp,it_temp) + &
!!$                 & wgt2_lat(ic)*pvi%a8(inmw2_lat,lev_temp,it_temp)
!!$  END DO
!!$END DO
!!$
!!$DO ic=jcb:jce
!!$  DO ilev=1,nlev
!!$    do3dt(ic,ilev) = al1(ic,ilev) + &
!!$                   & al2(ic,ilev)*(avi%o3_vmr(ic,ilev)-al3(ic,ilev) + &
!!$                   & al4(ic,ilev)*(avi%tmprt(ic,ilev)-al5(ic,ilev) + &
!!$                   & al6(ic,ilev)*(o3_column(ic,ilev)-al7(ic,ilev) + &
!!$                   & al8(ic,ilev)*avi%o3_vmr(ic,ilev)
!!$  END DO
!!$END DO
END SUBROUTINE cariolle_do3dt
