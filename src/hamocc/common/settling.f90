!>
!! @file settling.f90
!! @brief compute settling of debris
#include "hamocc_omp_definitions.inc"

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

      USE mo_kind,    ONLY        : wp
      USE mo_sedmnt,  ONLY        : prorca, prcaca, silpro, produs

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

        IF(pddpo(j,kpke) > 0.5_wp)THEN
 
           ! sediment fluxes at the bottom
            prorca(j) = bgctra(j,kpke,idet )*wpoc(kpke)
            prcaca(j) = bgctra(j,kpke,icalc)*wcal
            silpro(j) = bgctra(j,kpke,iopal)*wopal
            produs(j) = bgctra(j,kpke,idust)*wdust

            bgcflux(j,kprorca) = prorca(j)
            bgcflux(j,kprcaca) = prcaca(j)
            bgcflux(j,ksilpro) = silpro(j)
            bgcflux(j,kprodus) = produs(j)
            !
        ENDIF
         ENDIF
      END DO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
      END SUBROUTINE settling 

