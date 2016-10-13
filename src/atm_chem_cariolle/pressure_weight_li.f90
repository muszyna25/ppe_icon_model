SUBROUTINE pressure_weight_li (jcb,jce,NCX,nlev,plev,nlev_clim,plev_clim,wgt1_p,wgt2_p,inmw1_p,inmw2_p)
USE mo_cariolle_kind,       ONLY: wi,wp
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
write(*,*) '---------------------------------------------------'
write(*,*) 'plev_clim=',plev_clim
write(*,*) '---------------------------------------------------'
write(*,*) 'plev(1,1)=',plev(1,1)
write(*,*) 'inmw1_p(1,1)=',inmw1_p(1,1),'inmw2_p(1,1)=',inmw2_p(1,1)
write(*,*) 'plev_clim(inmw1_p(1,1))=',plev_clim(inmw1_p(1,1)),&
          &'plev_clim(inmw2_p(1,1))',plev_clim(inmw2_p(1,1))
write(*,*) 'wgt1_p(1,1)=',wgt1_p(1,1),'wgt2(1,1)=',wgt2_p(1,1)
write(*,*) 'plev(1,5)=',plev(1,5)
write(*,*) 'inmw1_p(1,5)=',inmw1_p(1,5),'inmw2_p(1,5)=',inmw2_p(1,5)
write(*,*) 'plev_clim(inmw1_p(1,5))=',plev_clim(inmw1_p(1,5)),&
          &'plev_clim(inmw2_p(1,5))',plev_clim(inmw2_p(1,5))
write(*,*) 'wgt1_p(1,5)=',wgt1_p(1,5),'wgt2(1,5)=',wgt2_p(1,5)
write(*,*) 'plev(1,nlev)=',plev(1,nlev)
write(*,*) 'inmw1_p(1,nlev)=',inmw1_p(1,nlev),'inmw2_p(1,nlev)=',inmw2_p(1,nlev)
write(*,*) 'plev_clim(inmw1_p(1,nlev))=',plev_clim(inmw1_p(1,nlev)),&
          &'plev_clim(inmw2_p(1,nlev))',plev_clim(inmw2_p(1,nlev))
write(*,*) 'wgt1_p(1,nlev)=',wgt1_p(1,nlev),'wgt2(1,nlev)=',wgt2_p(1,nlev)
END SUBROUTINE pressure_weight_li
