#include "hamocc_omp_definitions.inc"

MODULE mo_settling


      USE mo_kind,    ONLY        : wp

      IMPLICIT NONE

      PRIVATE

      PUBLIC:: settling, settling_pdm


CONTAINS

!>
!! @file settling.f90
!! @brief compute settling of debris
      SUBROUTINE settling (klev,start_idx, end_idx,pddpo,za)

      USE mo_param1_bgc, ONLY     : icalc, iopal, kopex90,   &
       &                            idet, kcalex90, &
       &                            kprorca,kprcaca,   &
       &                            ksilpro,kcoex90, idust, &
       &                            kprodus, kopex1000, &
       &                            kopex2000,kcoex1000,&
       &                            kcoex2000, kcalex1000, &
       &                            kcalex2000

      USE mo_carbch,  ONLY        : bgctra, wpoc, bgcflux

      USE mo_biomod,  ONLY        : wopal, wcal, wdust, n90depth,&
      &                             n1000depth,n2000depth

      USE mo_sedmnt,  ONLY        : prorca, prcaca, silpro, produs, kbo

      USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs, dtbgc

      IMPLICIT NONE


      ! Arguments

      INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
      INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
      INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

      REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
      REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height


      ! Local variables
      INTEGER :: k, kpke,j


     ! implicit method:
     ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
     ! -->
     ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
     ! sedimentation=w*dt*C(ks,T+dt)
     !
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,kpke,k) HAMOCC_OMP_DEFAULT_SCHEDULE

       DO j=start_idx,end_idx
        kpke=klev(j)        

        IF(kpke > 0)THEN
        IF(pddpo(j,1) > 0.5_wp)THEN

         if(kpke>=n90depth)then
           bgcflux(j,kcoex90) = bgctra(j,n90depth,idet)*wpoc(n90depth)/dtbgc  
           bgcflux(j,kcalex90) = bgctra(j,n90depth,icalc)*wcal/dtbgc  
           bgcflux(j,kopex90) = bgctra(j,n90depth,iopal)*wopal/dtbgc  
         endif

         if(kpke>=n1000depth)then
           bgcflux(j,kcoex1000) = bgctra(j,n1000depth,idet)*wpoc(n1000depth)/dtbgc  
           bgcflux(j,kcalex1000) = bgctra(j,n1000depth,icalc)*wcal/dtbgc  
           bgcflux(j,kopex1000) = bgctra(j,n1000depth,iopal)*wopal/dtbgc  
         endif

         if(kpke>=n2000depth)then
           bgcflux(j,kcoex2000) = bgctra(j,n2000depth,idet)*wpoc(n2000depth)/dtbgc  
           bgcflux(j,kcalex2000) = bgctra(j,n2000depth,icalc)*wcal/dtbgc  
           bgcflux(j,kopex2000) = bgctra(j,n2000depth,iopal)*wopal/dtbgc  
         endif

                  ! -----------surface layer
         k=1 
           bgctra(j,k,idet)  = (bgctra(j,k,idet)*(pddpo(j,k)+za(j)))      &
                &              / (pddpo(j,k)+ za(j) +wpoc(k))

           bgctra(j,k,icalc)   = (bgctra(j,k,icalc)*(pddpo(j,k)+za(j)))   &
                &                / (pddpo(j,k)+za(j)+wcal)


           bgctra(j,k,iopal)   = (bgctra(j,k,iopal)*(pddpo(j,k)+za(j)))    &
                &                / (pddpo(j,k)+za(j)+wopal)

           bgctra(j,k,idust)   = (bgctra(j,k,idust)*(pddpo(j,k)+za(j)))    &
                &                / (pddpo(j,k)+za(j)+wdust)
          ENDIF

         DO k=2,kpke
          IF(pddpo(j,k) > 0.5_wp)THEN
          ! water column
              bgctra(j,k,idet)  = (bgctra(j,k  ,idet)*pddpo(j,k)    &
                   &        +  bgctra(j,k-1,idet)*wpoc(k-1))/  &
                   &                     (pddpo(j,k)+wpoc(k))
              bgctra(j,k,icalc)   = (bgctra(j,k  ,icalc)*pddpo(j,k)   &
                   &                +  bgctra(j,k-1,icalc)*wcal)/         &
                   &                           (pddpo(j,k)+wcal)

              bgctra(j,k,iopal)   = (bgctra(j,k  ,iopal)*pddpo(j,k)   &
                   &                +  bgctra(j,k-1,iopal)*wopal)/        &
                   &                           (pddpo(j,k)+wopal)

              bgctra(j,k,idust)   = (bgctra(j,k  ,idust)*pddpo(j,k)   &
                   &                +  bgctra(j,k-1,idust)*wdust)/        &
                   &                           (pddpo(j,k)+wdust)

          ENDIF
         ENDDO

        IF(pddpo(j,kbo(j)) > 0.5_wp)THEN
 
           ! sediment fluxes at the bottom
            prorca(j) = bgctra(j,kbo(j),idet )*wpoc(kbo(j))
            prcaca(j) = bgctra(j,kbo(j),icalc)*wcal
            silpro(j) = bgctra(j,kbo(j),iopal)*wopal
            produs(j) = bgctra(j,kbo(j),idust)*wdust
            bgcflux(j,kprorca) = prorca(j)/dtbgc
            bgcflux(j,kprcaca) = prcaca(j)/dtbgc
            bgcflux(j,ksilpro) = silpro(j)/dtbgc
            bgcflux(j,kprodus) = produs(j)/dtbgc
            !
        ENDIF
         ENDIF
      END DO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
      END SUBROUTINE settling 


!>
!! @file settling.f90
!! @brief compute settling of debris
      SUBROUTINE settling_pdm (klev,start_idx, end_idx,pddpo,za)

      USE mo_param1_bgc, ONLY     : icalc, iopal, kopex90,   &
       &                            idet, kcalex90, &
       &                            kprorca,kprcaca,   &
       &                            ksilpro,kcoex90, idust, &
       &                            kprodus, kopex1000, &
       &                            kopex2000,kcoex1000,&
       &                            kcoex2000, kcalex1000, &
       &                            kcalex2000

      USE mo_carbch,  ONLY        : bgctra, wpoc, bgcflux

      USE mo_biomod,  ONLY        : wopal, wcal, wdust, n90depth,&
      &                             n1000depth,n2000depth

      USE mo_sedmnt,  ONLY        : prorca, prcaca, silpro, produs, kbo

      USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs, inv_dtbgc

      IMPLICIT NONE


      ! Arguments

      INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
      INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
      INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

      REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
      REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height


      ! Local variables
      INTEGER :: k, kpke,j

      REAL(wp) :: det_flux(bgc_zlevs+1), calc_flux(bgc_zlevs+1)
      REAL(wp) :: opal_flux(bgc_zlevs+1), dust_flux(bgc_zlevs+1)
      REAL(wp) :: cn_det, cn_calc, cn_opal, cn_dust 
      REAL(wp) :: det_plim, calc_plim, opal_plim, dust_plim


!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,cn_det,cn_calc,cn_opal,cn_dust,&
!HAMOCC_OMP            det_plim, calc_plim, opal_plim, dust_plim,&
!HAMOCC_OMP            det_flux, calc_flux, opal_flux, dust_flux ) HAMOCC_OMP_DEFAULT_SCHEDULE

       DO j=start_idx,end_idx
        kpke=klev(j)        

        IF(kpke > 0)THEN
        IF(pddpo(j,1) > 0.5_wp)THEN
        
         DO k = 1,kpke + 1

           det_flux(k)  = 0._wp
           calc_flux(k) = 0._wp
           opal_flux(k) = 0._wp
           dust_flux(k) = 0._wp

         ENDDO

        ! flux from atmosphere to surface layer euqals zero (here)
        ! hence, starting calculation of flux with k=2

         IF( kbo(j) > 1 ) THEN ! only, when water column is larger than one box 
                             ! otherwise go immediately to the calculation 
                             ! for flux to sediment

            k = 2 ! flux at the bottom of first level
            cn_det  = wpoc(k-1)/pddpo(j,k-1) ! assuming both positive      
            cn_calc = wcal/pddpo(j,k-1)      ! assuming both positive
            cn_opal = wopal/pddpo(j,k-1)     ! assuming both positive
            cn_dust = wdust/pddpo(j,k-1)     ! assuming both positive
         
            det_plim    = 0._wp
            calc_plim   = 0._wp
            opal_plim   = 0._wp
            dust_plim   = 0._wp

            CALL plimiter(det_plim,  cn_det,  bgctra(j,k,idet), bgctra(j,k-1,idet),   0._wp)
            CALL plimiter(calc_plim, cn_calc, bgctra(j,k,icalc), bgctra(j,k-1,icalc),  0._wp)
            CALL plimiter(opal_plim, cn_opal, bgctra(j,k,iopal), bgctra(j,k-1,iopal),  0._wp)
            CALL plimiter(dust_plim, cn_dust, bgctra(j,k,idust),bgctra(j,k-1,idust), 0._wp)

            det_flux(k)  = wpoc(k-1) * (bgctra(j,k-1,idet)                    &
                              & + 0.5_wp * det_plim * (1._wp-cn_det)            &
                              & * (bgctra(j,k,idet) - bgctra(j,k-1,idet))) 

            calc_flux(k) = wcal      * (bgctra(j,k-1,icalc)                   &
                              & + 0.5_wp * calc_plim * (1._wp-cn_calc)          &
                              & * (bgctra(j,k,icalc) - bgctra(j,k-1,icalc)))

            opal_flux(k) = wopal     * (bgctra(j,k-1,iopal)                   &
                              & + 0.5_wp * opal_plim * (1._wp-cn_opal)          &
                              & * (bgctra(j,k,iopal) - bgctra(j,k-1,iopal)))

            dust_flux(k) = wdust * (bgctra(j,k-1,idust)                      &
                              & + 0.5_wp * dust_plim * (1._wp-cn_dust)          &
                            & * (bgctra(j,k,idust) - bgctra(j,k-1,idust)))
        
            DO k = 3, kbo(j) ! water column

              cn_det  = wpoc(k-1)/pddpo(j,k-1) ! assuming both positive      
              cn_calc = wcal/pddpo(j,k-1)      ! assuming both positive
              cn_opal = wopal/pddpo(j,k-1)     ! assuming both positive
              cn_dust = wdust/pddpo(j,k-1)     ! assuming both positive
         
              det_plim    = 0._wp
              calc_plim   = 0._wp
              opal_plim   = 0._wp
              dust_plim   = 0._wp

              CALL plimiter(det_plim, cn_det, bgctra(j,k,idet), bgctra(j,k-1,idet), bgctra(j,k-2,idet))
              CALL plimiter(calc_plim, cn_calc, bgctra(j,k,icalc), bgctra(j,k-1,icalc), bgctra(j,k-2,icalc))
              CALL plimiter(opal_plim, cn_opal, bgctra(j,k,iopal), bgctra(j,k-1,iopal), bgctra(j,k-2,iopal))
              CALL plimiter(dust_plim, cn_dust, bgctra(j,k,idust), bgctra(j,k-1,idust), bgctra(j,k-2,idust))

              det_flux(k)  = wpoc(k-1) * (bgctra(j,k-1,idet)                  &
                               & + 0.5_wp * det_plim  * (1._wp-cn_det)          &
                               & * (bgctra(j,k,idet) - bgctra(j,k-1,idet))) 

              calc_flux(k) = wcal      * (bgctra(j,k-1,icalc)                 &
                               & + 0.5_wp * calc_plim * (1._wp-cn_calc)         &
                              & * (bgctra(j,k,icalc) - bgctra(j,k-1,icalc)))

              opal_flux(k) = wopal     * (bgctra(j,k-1,iopal)                 &
                               & + 0.5_wp * opal_plim * (1._wp-cn_opal)         &
                              & * (bgctra(j,k,iopal) - bgctra(j,k-1,iopal)))

              dust_flux(k) = wdust * (bgctra(j,k-1,idust)                    &
                              & + 0.5_wp * dust_plim * (1._wp-cn_dust)          &
                            & * (bgctra(j,k,idust) - bgctra(j,k-1,idust)))

            ENDDO
         ENDIF !kbo > 1?

         k = kbo(j) + 1 ! flux to sediment
         det_flux(k)  = wpoc(k-1) * bgctra(j,k-1,idet) 
         calc_flux(k) = wcal      * bgctra(j,k-1,icalc)
         opal_flux(k) = wopal     * bgctra(j,k-1,iopal)
         dust_flux(k) = wdust     * bgctra(j,k-1,idust)
          
        
         ! calculating change of concentrations water column:
         DO k = 1,kbo(j)
             bgctra(j,k,idet) = bgctra(j,k,idet) + (det_flux(k) - det_flux(k+1))/pddpo(j,k)
             bgctra(j,k,icalc) = bgctra(j,k,icalc) + (calc_flux(k) - calc_flux(k+1))/pddpo(j,k)
             bgctra(j,k,iopal) = bgctra(j,k,iopal) + (opal_flux(k) - opal_flux(k+1))/pddpo(j,k)
             bgctra(j,k,idust) = bgctra(j,k,idust) + (dust_flux(k) - dust_flux(k+1))/pddpo(j,k)
         ENDDO

         ! here n90depth+1, since the interface from box k to k+1 is w-level k+1
         ! in the previous implicit scheme, it's defined as k-flux
         IF (kbo(j) > n90depth) THEN
           bgcflux(j,kcoex90) = det_flux(n90depth+1)*inv_dtbgc   
           bgcflux(j,kopex90) = opal_flux(n90depth+1)*inv_dtbgc 
           bgcflux(j,kcalex90) = calc_flux(n90depth+1)*inv_dtbgc 
         ENDIF

         ! write out fluxes at about 1000 m
         IF (kbo(j) > n1000depth) THEN
           bgcflux(j,kcoex1000) = det_flux(n1000depth+1)*inv_dtbgc     
           bgcflux(j,kopex1000) = opal_flux(n1000depth+1)*inv_dtbgc    
           bgcflux(j,kcalex1000) = calc_flux(n1000depth+1)*inv_dtbgc    
         ENDIF

         ! write out fluxes at about 1950 m
         IF (kbo(j) > n2000depth) THEN
           bgcflux(j,kcoex2000) = det_flux(n2000depth+1)*inv_dtbgc     
           bgcflux(j,kopex2000) = opal_flux(n2000depth+1)*inv_dtbgc    
           bgcflux(j,kcalex2000) = calc_flux(n2000depth+1)*inv_dtbgc    
         ENDIF

         IF (pddpo(j,1) > 0.5_wp) THEN

           prorca(j) = det_flux(kbo(j)+1)
           prcaca(j) = calc_flux(kbo(j)+1)
           silpro(j) = opal_flux(kbo(j)+1)
           produs(j) = dust_flux(kbo(j)+1)
           
           !
           !  write output
           !
           bgcflux(j,kprorca) = prorca(j)*inv_dtbgc     
           bgcflux(j,kprcaca) = prcaca(j)*inv_dtbgc     
           bgcflux(j,ksilpro) = silpro(j)*inv_dtbgc     
           bgcflux(j,kprodus) = produs(j)*inv_dtbgc     

         ENDIF 

        ENDIF ! ddpo > 0.5   
        ENDIF !    
       ENDDO ! j

!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

      END SUBROUTINE settling_pdm 


      SUBROUTINE plimiter(p_lim, cnnr, conc_p1, conc_m, conc_m1)

      IMPLICIT NONE

      REAL(wp), INTENT(inout) :: p_lim
      REAL(wp), INTENT(in)    :: cnnr
      REAL(wp), INTENT(in)    :: conc_p1
      REAL(wp), INTENT(in)    :: conc_m
      REAL(wp), INTENT(in)    :: conc_m1

      REAL(wp)                :: ratio_p, limiter, delta_conc

      delta_conc = conc_p1 - conc_m
      IF(ABS(delta_conc) > 0._wp) THEN ! unfortunately required, as it can rarely happen that it's equal 0
         ratio_p    = (conc_m - conc_m1) / delta_conc
      ELSE 
         ! range, where limiter 2._wp/(1._wp-cnnr)) should jump in
         ! don't have to care about (1-cn), since model is close to unstable anyway, when cn=1 
         ratio_p    = 1.e37_wp ! fix me - is there a value for the largest real value defined?
      ENDIF
      

      limiter   = 0.5_wp+1._wp/6._wp*(1._wp-2._wp*cnnr)  &
                   & +  (0.5_wp-1._wp/6._wp*(1._wp-2._wp*cnnr)) * ratio_p
      p_lim  = MAX(0._wp,MIN(MIN(limiter,2._wp/(1._wp-cnnr)),2._wp * ratio_p / cnnr))
     
      END SUBROUTINE plimiter
END MODULE
