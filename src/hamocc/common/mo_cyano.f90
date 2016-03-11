!>
!! @file mo_cyano.f90
!! @brief N2 fixation and cyanobateria dynamics.
!!
!! Contains computation of diagnostic and prognostic N2 fixation and.
!! cyanobacteria dyanmics (growth, decay, buouyancy) 
!!
!!

#include "hamocc_omp_definitions.inc"

MODULE mo_cyano

  USE mo_kind, ONLY           : wp
  USE mo_control_bgc, ONLY    : dtb, dtbgc, bgc_nproma, bgc_zlevs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cyano, &
            cyadyn

CONTAINS

SUBROUTINE cyano ( start_idx,end_idx,pddpo, za )

!! @brief diagostic N2 fixation

  USE mo_biomod, ONLY         : rnit, n2_fixation, rn2
  USE mo_carbch, ONLY         : bgctra, bgcflux, bgctend 
  USE mo_param1_bgc, ONLY     : iano3, iphosph, igasnit, &
       &                        ioxygen, ialkali, knfixd, &
       &                        kn2b

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in) :: start_idx            !< start index for j loop (ICON cells, MPIOM lat dir)   
  INTEGER, INTENT(in) :: end_idx              !< end index  for j loop  (ICON cells, MPIOM lat dir) 

  REAL(wp), INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd dimension) [m].
  REAL(wp), INTENT(in) :: za(bgc_nproma)              !< surface height [m].

  !! Local variables

  REAL(wp) :: oldnitrate
  INTEGER  :: j
  !
  ! --------------------------------------------------------------------
  !
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,oldnitrate) HAMOCC_OMP_DEFAULT_SCHEDULE
     DO j = start_idx, end_idx

        IF (pddpo(j,1) > 0.5_wp) THEN
           IF (bgctra(j,1,iano3) < (rnit*bgctra(j,1,iphosph))) THEN

              oldnitrate = bgctra(j,1,iano3)

              bgctra(j,1,iano3) = bgctra(j,1,iano3)                 &
                                     * (1._wp - n2_fixation*dtb)              &
                                     + n2_fixation * dtb * rnit * bgctra(j,1,iphosph)


              bgctra(j,1,igasnit)  = bgctra(j,1,igasnit)                  &
                                     - (bgctra(j,1,iano3)-oldnitrate) /rn2  


         
              bgctra(j,1,ioxygen)  = bgctra(j,1,ioxygen)                  &
                                     - (bgctra(j,1,iano3)-oldnitrate) * 1.5_wp  ! factor 1.5 results from: NO3 = N + 3O = N + 1.5*O2 
                                                                                 ! (O2 is used for the balancing)

              !alkalinity is decreased during n-fixation due to release of H+ 
              ! change in ALK is identical to change in NO3

              bgctra(j,1,ialkali)  = bgctra(j,1,ialkali)                  &
                   &                 - (bgctra(j,1,iano3)-oldnitrate)
             
              bgcflux(j,knfixd)    = (- oldnitrate + bgctra(j,1,iano3))/rn2 /dtbgc
              bgctend(j,1,kn2b)    = bgctend(j,1,kn2b) - (bgctra(j,1,iano3) - oldnitrate) * (pddpo(j,1) + za(j))
           ENDIF
        ENDIF

     ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
END SUBROUTINE cyano



SUBROUTINE cyadyn(klevs,start_idx,end_idx,pddpo,za,ptho)
!! @brief prognostic N2 fixation, cyanobacteria

      USE mo_biomod, ONLY         : cycdec, pi_alpha_cya,cya_growth_max,          &
       &                            Topt_cya,T1_cya,T2_cya,bkcya_N, bkcya_P,      &
       &                            fPAR, strahl, ro2ut,       &
       &                            doccya_fac, rnit, riron, rcar, rn2, &
       &                            strahl,bkcya_fe,   wcya, rnoi

      USE mo_carbch, ONLY         : bgctra, bgctend, swr_frac
      USE mo_param1_bgc, ONLY     : iano3, iphosph, igasnit, &
           &                        ioxygen, ialkali, icya,  &
           &                        isco212, idoccya, &
           &                        idet, iiron, knfix, &
           &                        kpho_cya, kcyaloss, kn2b


      IMPLICIT NONE

      INTEGER, INTENT(in) ::  start_idx                  !< 1st REAL of model grid
      INTEGER, INTENT(in) ::  end_idx                  !< 1st REAL of model grid
      INTEGER, INTENT(in) ::  klevs(bgc_nproma)                   !< 3rd (vertical) REAL of model grid.
  
      REAL(wp), INTENT(in) :: pddpo(bgc_nproma, bgc_zlevs) !< size of scalar grid cell (3rd dimension) [m].
      REAL(wp), INTENT(in) :: ptho(bgc_nproma, bgc_zlevs)  !< potential temperature [deg C]
      REAL(wp), INTENT(in) :: za(bgc_nproma)  !< potential temperature [deg C]

      !! Local variables
   
      INTEGER  :: j,k, kpke
      REAL(wp) :: oldigasnit                                                                          
      REAL(wp) :: cyapro,cyaloss
      REAL(wp) :: avanut,avcyabac                                             
      REAL(wp) :: pho,xn,avanfe,pho_fe,pho_p 
      REAL(wp) :: l_I, l_T
      REAL(wp) :: T_min_Topt,sgnT
      REAL(wp) :: xa_P, xa_fe, avnit,l_P,l_fe
      REAL(wp) :: xn_p,xn_fe
   
!HAMOCC_OMP_PARALLEL 
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,avcyabac,avanut,avanfe,avnit,l_I,T_min_Topt,&
!HAMOCC_OMP            sgnT,l_T,xa_p,l_P,xa_fe,pho_fe,pho_p,xn_p,xn_fe,pho,&
!HAMOCC_OMP            cyapro,oldigasnit,xn,cyaloss) HAMOCC_OMP_DEFAULT_SCHEDULE

  DO j = start_idx, end_idx
  
    kpke=klevs(j)

      DO k=1, kpke                     

            IF( pddpo(j,k) .GT. 0.5_wp ) THEN
                   

              avcyabac = MAX(1.e-11_wp,bgctra(j,k,icya))                !available cyanobacteria
              avanut = MAX(0._wp,bgctra(j,k,iphosph))                   !available phosphate    
              avanfe = MAX(0._wp,(bgctra(j,k,iiron)/riron))             !available iron
              avnit = MAX(0._wp,bgctra(j,k,iano3)/rnit)                 !available nitrate
                        
 
              l_I = (pi_alpha_cya*fPAR*strahl(j))*swr_frac(j,k) &     ! light limitation
                /SQRT(cya_growth_max**2 + (pi_alpha_cya**2)*(fPAR*strahl(j)*swr_frac(j,k))**2) 
          
              T_min_Topt = ptho(j,k)-Topt_cya                          
              sgnT = sign(1._wp,T_min_Topt)
              IF(T_min_Topt .eq. 0._wp) sgnT = 0._wp 

              l_T = exp(-((T_min_Topt**4)/(T1_cya-T2_cya*sgnT)**4))          !temperature limitation 

              xa_p = avanut                                                 
              l_P = xa_P / (bkcya_P + xa_P)                      !phosphate limitation 

              xa_fe = avanfe         
              l_fe = xa_fe / (bkcya_fe + xa_fe)                  !iron limitation

              pho_fe = dtb*cya_growth_max*l_I*l_T*l_fe  
              pho_p  = dtb*cya_growth_max*l_I*l_T*l_P

        
              pho=dtb*cya_growth_max*l_I*l_T*l_P*l_fe            !growth 

              xn_p = xa_p/ (1._wp + pho_fe*avcyabac/(bkcya_p + xa_p))
              xn_fe = xa_fe/ (1._wp + pho_p*avcyabac/(bkcya_fe + xa_fe))

              pho = MIN(xa_p-xn_p, xa_fe-xn_fe)

             
              cyapro = MIN(avnit-1.e-11_wp*rnoi, &
   &                       avcyabac*pho*avnit**2/(bkcya_N**2 + avnit**2))



              ! ---------- nutrient uptake

              bgctra(j,k,iphosph) = bgctra(j,k,iphosph) - pho    
              bgctra(j,k,iiron) = bgctra(j,k,iiron) - pho * riron  

              bgctra(j,k,igasnit) = bgctra(j,k,igasnit) - (pho - cyapro)*rnit/rn2  ! gasnit [N2]
  
           

              bgctend(j,k,knfix) =  (pho - cyapro)*rnit/rn2/dtbgc  ! N fixation
              bgctend(j,k,kpho_cya) =  pho/dtbgc   
              bgctend(j,k,kn2b) = bgctend(j,k,kn2b) - (pho -cyapro) * rnit
 
              bgctra(j,k,iano3) = bgctra(j,k,iano3) - cyapro*rnit                                                           

              ! ---------- change of alkalinity 
              !(only for production on nitrate, no change of alkalinity 
              !for N2 fixation Wolf-Gladrow et al.(2007))
              bgctra(j,k,ialkali) = bgctra(j,k,ialkali) + cyapro*rnit

              ! ---------- oxygen production
              ! O2 from cyano growth using NO3: cyapro [NO3] --> ro2ut 
              ! O2 from cyano growth fixing N2: (pho - cyapro) [NO3]
              ! -> ro2ut - 3 * rnit/rno2 = 172 - 3*16/2 = 148
              bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen) + (pho - cyapro)*148._wp &
            &                       + cyapro*ro2ut                                            

              ! --------- change of total CO2
              bgctra(j,k,isco212) = bgctra(j,k,isco212) - pho*rcar

              ! --------- decay of cyanobacteria
              xn = bgctra(j,k,icya)/(1._wp + cycdec *dtb)
              cyaloss = max(0._wp,bgctra(j,k,icya) - xn)                                

              ! --------- change of cyanobacteria
              bgctra(j,k,icya) = bgctra(j,k,icya) + pho - cyaloss
              bgctend(j,k,kcyaloss) =  cyaloss/dtbgc   

              ! --------- decaying cyanobacteria are distributed to DOCCYA and detritus  
              bgctra(j,k,idoccya) = bgctra(j,k,idoccya) + doccya_fac*cyaloss     
              bgctra(j,k,idet) = bgctra(j,k,idet) + (1.0_wp - doccya_fac)*cyaloss        




            ENDIF ! wet cells
      ENDDO ! 
  ENDDO ! 
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL


!HAMOCC_OMP_PARALLEL 
!HAMOCC_OMP_DO PRIVATE(j,kpke,k) HAMOCC_OMP_DEFAULT_SCHEDULE

! -------------- buoyancy of cyanobacteria----------------------------------------

  ! implicit method:
  ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
  ! -->
  ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
  ! sedimentation=w*dt*C(ks,T+dt)
  !

 DO j=start_idx,end_idx
    kpke=klevs(j)
    IF (kpke > 0)THEN
    IF(pddpo(j,kpke).GT.0.5_wp) THEN
           bgctra(j,kpke,icya)  = (bgctra(j,kpke,icya)*pddpo(j,kpke))      &
                   &              / (pddpo(j,kpke)+wcya)
             
    ENDIF

    do k=2,kpke-1
         ! water column
        if(pddpo(j,k+1).LE.0.5_wp)then ! last wet cell
              bgctra(j,k,icya)  = (bgctra(j,k,icya)*pddpo(j,k))      &
                   &              / (pddpo(j,k)+wcya)

         else
               bgctra(j,k,icya)    = (bgctra(j,k  ,icya)*pddpo(j,k)    &
                   &                +  bgctra(j,k+1,icya)*wcya)/          &
                   &                          (pddpo(j,k)+wcya)
         endif 
    
   ENDDO
    k=1
    IF((pddpo(j,k).GT.0.5_wp) .and. (pddpo(j,k+1).GT.0.5_wp) )then ! only if next cell also wet
         bgctra(j,k,icya)  =  bgctra(j,k,icya) + (wcya*bgctra(j,k+1,icya))/(pddpo(j,k)+za(j))
    endif   
   ENDIF
 ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

                  
 
END SUBROUTINE  cyadyn

END MODULE


