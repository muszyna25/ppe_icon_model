     !>
     !! @file settling.f90
     !! @brief compute settling of debris
  !!
  !! Notes:
  !!
  !!
  !!

  !!
      SUBROUTINE settling (klev,start_idx, end_idx,pddpo,za)

      USE mo_param1_bgc, ONLY     : icalc, iopal,   &
       &                            idet, &
       &                            kprorca,kprcaca,   &
       &                            ksilpro,kcoex90, idust, kprodus

      USE mo_carbch,  ONLY        : bgctra, wpoc, bgcflux

      USE mo_biomod,  ONLY        : wopal, wcal, wdust, kbo, n90depth

      USE mo_kind,    ONLY        : wp
      USE mo_sedmnt,  ONLY        : prorca, prcaca, silpro, produs

      USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs, dtbgc

      IMPLICIT NONE


      INTEGER, INTENT(in)    :: start_idx             !<  number of levels
      INTEGER, INTENT(in)    :: end_idx               !<  number of columns to be treated at once 
      INTEGER, INTENT(in), target    :: klev(bgc_nproma) 
      REAL(wp), INTENT(in), target   :: pddpo(bgc_nproma,bgc_zlevs)
      REAL(wp), INTENT(in), target   :: za(bgc_nproma)

      ! Local variables
      INTEGER :: k, kpke,j


     ! implicit method:
     ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
     ! -->
     ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
     ! sedimentation=w*dt*C(ks,T+dt)
     !


       DO j=start_idx,end_idx
        kpke=klev(j)        
        IF(kpke > 0)THEN
        IF(pddpo(j,1) > 0.5_wp)THEN

         if(kpke>=n90depth)bgcflux(j,kcoex90) = bgctra(j,n90depth,idet)*wpoc(n90depth)/dtbgc  

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

      END SUBROUTINE settling 

