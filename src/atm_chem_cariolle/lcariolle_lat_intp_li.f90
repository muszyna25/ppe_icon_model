!>
!! @brief Interpolation weights for linear interpolation from a
!! regular grid between N pole and S pole. The N pole and S pole have to
!! be included, but the distance to the first (last) point of the regular
!! grid can be different. 
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
SUBROUTINE lcariolle_lat_intp_li(                                       &
         & jcb,                jce,                   NCX,              &
         & cell_center_lat,    n_lat_clim,            r_lat_clim,       &
         & r_delta_lat_clim,   l_lat_clim_sn,         wgt1_lat,         &
         & wgt2_lat,           inmw1_lat,             inmw2_lat         )
USE mo_lcariolle_kind,        ONLY: wi,wp
INTEGER(wi),INTENT(IN)    :: jcb, jce, & !< first, last point to interpolate on
                           & NCX         !< length of array in calling subprog.
REAL(wp),INTENT(IN)       :: cell_center_lat(NCX) !< latitudes in radiant
INTEGER(wi),INTENT(IN)    :: n_lat_clim           !< number of lats in climatol.
REAL(wp),INTENT(IN)       :: r_lat_clim(0:n_lat_clim+1), & !< lats in climatology
                                                           !< in radiant
                           & r_delta_lat_clim     !< spacing in climatology
LOGICAL,INTENT(IN)        :: l_lat_clim_sn        !< .true. if order is S->N,
                                                  !< .false. otherwise
REAL(wp),INTENT(OUT)      :: wgt1_lat(NCX),wgt2_lat(NCX)   !< interpolation weights
INTEGER(wi),INTENT(OUT)   :: inmw1_lat(NCX),inmw2_lat(NCX) !< corresponding indices 
INTEGER(wi)               :: n_order
REAL(wp)                  :: r_delta_lat_clim_i

IF (l_lat_clim_sn) THEN
  n_order=1
ELSE
  n_order=-1
END IF
r_delta_lat_clim_i=1._wp/r_delta_lat_clim
inmw1_lat(jcb:jce)=MAX(INT(n_order*(cell_center_lat(jcb:jce)-r_lat_clim(1))*r_delta_lat_clim_i+1._wp),0)
inmw2_lat(jcb:jce)=inmw1_lat(jcb:jce)+1
wgt2_lat(jcb:jce)=n_order*(cell_center_lat(jcb:jce)-r_lat_clim(inmw1_lat(jcb:jce)))*r_delta_lat_clim_i
wgt1_lat(jcb:jce)=1.0_wp-wgt2_lat(jcb:jce)
END SUBROUTINE lcariolle_lat_intp_li
