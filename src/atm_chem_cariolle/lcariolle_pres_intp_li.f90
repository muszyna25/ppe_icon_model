SUBROUTINE lcariolle_pres_intp_li (                               &
         & jcb,              jce,                 NCX,            &
         & nlev,             plev,                nlev_clim,      &
         & plev_clim,        wgt1_p,              wgt2_p,         &
         & inmw1_p,          inmw2_p                              )
USE mo_lcariolle_kind,       ONLY: wi,wp
INTEGER(wi),INTENT(IN)         :: jcb,jce,NCX,nlev
REAL(wp),INTENT(IN)            :: plev(NCX,nlev)
INTEGER(wi),INTENT(IN)         :: nlev_clim
REAL(wp),INTENT(IN)            :: plev_clim(nlev_clim)
REAL(wp),INTENT(OUT)           :: wgt1_p(NCX,nlev),wgt2_p(NCX,nlev)
INTEGER(wi),INTENT(OUT)        :: inmw1_p(NCX,nlev),inmw2_p(NCX,nlev)
INTEGER(wi)                    :: ic,ilev,ii,ilev_low
DO ic=jcb,jce
  ilev_low=1           
  DO ilev=1,nlev
    ii=ilev_low
    DO
!      write(*,*) 'ii=',ii
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
