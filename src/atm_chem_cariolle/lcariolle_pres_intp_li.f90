!>
!! @brief This subroutine calculates interpolation weights for a linear
!! interpolation with respect to pressure layers.
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
SUBROUTINE lcariolle_pres_intp_li (                               &
         & jcb,              jce,                 NCX,            &
         & nlev,             plev,                nlev_clim,      &
         & plev_clim,        wgt1_p,              wgt2_p,         &
         & inmw1_p,          inmw2_p                              )
USE mo_lcariolle_kind,       ONLY: wi,wp
INTEGER(wi),INTENT(IN)     :: &
     & jcb,jce,  & !< begin, end index of column
     & NCX,      & !< first dim of fields as in calling subprogram
     & nlev        !< number of levels in column
REAL(wp),INTENT(IN)        :: plev(NCX,nlev) !< pressure to interpolate on
INTEGER(wi),INTENT(IN)     :: nlev_clim      !< number of pressures in climat.
REAL(wp),INTENT(IN)        :: plev_clim(nlev_clim) !< climatological pressures
REAL(wp),INTENT(OUT)       :: wgt1_p(NCX,nlev),wgt2_p(NCX,nlev) !< intp. weights
INTEGER(wi),INTENT(OUT)    :: inmw1_p(NCX,nlev),inmw2_p(NCX,nlev) !< indices
INTEGER(wi)                :: ic,ilev,ii,ilev_low
DO ic=jcb,jce
  ilev_low=1           
  DO ilev=1,nlev
    ii=ilev_low
    DO
! Pressures are assumed to be in irregular distances in the climatology
      IF (ii>nlev_clim) THEN
        inmw1_p(ic,ilev)=nlev_clim
        inmw2_p(ic,ilev)=nlev_clim
        wgt1_p(ic,ilev)=0.5_wp
        wgt2_p(ic,ilev)=0.5_wp
        EXIT
      END IF
      IF(plev(ic,ilev).GE.plev_clim(ii)) THEN
        ii=ii+1
      ELSE
        IF(ii==1) THEN
          inmw1_p(ic,ilev)=1
          inmw2_p(ic,ilev)=1
          wgt1_p(ic,ilev)=0.5_wp
          wgt2_p(ic,ilev)=0.5_wp
        ELSE
          inmw1_p(ic,ilev)=ii-1
          inmw2_p(ic,ilev)=ii
          wgt1_p(ic,ilev)=(plev_clim(ii)-plev(ic,ilev))/ &
                         &(plev_clim(ii)-plev_clim(ii-1))
          wgt2_p(ic,ilev)=1._wp-wgt1_p(ic,ilev)
        END IF
        EXIT
      END IF
    END DO
    ilev_low=ii
  END DO
END DO
END SUBROUTINE lcariolle_pres_intp_li
