  !! @file ocprod.f90
  !! @brief Computes plankton dynamics, OM degradation
  !!
#include "hamocc_omp_definitions.inc"

  SUBROUTINE ocprod (klev,start_idx,end_idx, ptho,pddpo, za)
    
   USE mo_kind, ONLY           : wp
   
   USE mo_biomod, ONLY         :phytomi, grami, rnoi, riron, pi_alpha, &
       &                        fpar, bkphy,grazra, bkzoo, epsher,         &
       &                        zinges, drempoc, ro2ut, remido, dyphy, spemor,     &
       &                        gammaz, gammap, ecan, rnit, ropal, bkopal,         &
       &                        rcalc, rcar, dremopal, relaxfe, fesoly,            &
       &                        denitrification, nitdem, dremn2o,         &
       &                        n2prod, sulfate_reduction, strahl,                 &
       &                        remido_cya, thresh_aerob,   & 
       &                        thresh_sred, dmsp, calmax
 
   USE mo_carbch, ONLY         : satoxy, &
       &                         bgctra, swr_frac, bgctend


   USE mo_control_bgc, ONLY    : dtb, bgc_nproma, bgc_zlevs, dtbgc 
   
   USE mo_param1_bgc, ONLY     : icalc, iopal, ian2o, igasnit, idms, &
       &                         iphy, izoo, isilica, iphosph, &
       &                         iano3, ioxygen, idet, idoc, isco212, &
       &                         ialkali, idoccya, kphosy, kremin,   &
       &                         iiron, ksred, kdenit, kgraz, kbacfra, &
       &                         kn2b, kh2ob, kdelcar, kdelsil, kbacfrac, &
       &                         kdmsprod, kdmsuv, kdmsbac,keuexp, knlim, kflim,&
       &                         kplim, kgraton, kexudp, kexudz, kpdy,kzdy, kaou


   USE mo_hamocc_nml, ONLY     : l_cyadyn

    IMPLICIT NONE

    ! Arguments

    INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
    INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
    INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 
    REAL(wp), INTENT(in), TARGET   :: ptho(bgc_nproma,bgc_zlevs)       !<  potential temperature (degC)
    REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]

    REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height

   !  Local variables

    INTEGER                :: j,k, kpke
    REAL(wp) :: avphy,avanut,avanfe,xa,xn,ya,yn,phosy,                   &
       &      pho, phofa, temfa,                                         &
       &      avgra,grazing,avsil,graton,                                &
       &      gratpoc,grawa,bacfra,phymor,zoothresh,zoomor,excdoc,exud,  &
       &      export, delsil, delcar, remin,             &
       &      opalrem, remin2o, aou, refra,                      &
       &      bacfra_cya, tpremin, r_bacfra, r_bacfra_cya,   &
       &      r_remin, avoxy, rcyano

   REAL(wp) :: surface_height

   REAL(wp) :: dms_prod, dms_uv, dms_bac 
 
! check if cyanobacteria are switched on
 rcyano=merge(1._wp,0._wp,l_cyadyn)

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,surface_height,avphy,avgra,avsil,avanut,&
!HAMOCC_OMP            avanfe,phofa,temfa,pho,phosy,xa,xn,ya,yn,grazing,&
!HAMOCC_OMP            graton,gratpoc,grawa,phymor,zoomor,zoothresh,excdoc,&
!HAMOCC_OMP            exud,export,delsil,delcar,opalrem,dms_prod,dms_bac,&
!HAMOCC_OMP            dms_uv,avoxy,tpremin,r_bacfra,r_bacfra_cya,r_remin,&
!HAMOCC_OMP            bacfra,bacfra_cya,remin,aou,refra,remin2o) HAMOCC_OMP_DEFAULT_SCHEDULE

 DO j = start_idx, end_idx
  
     kpke=klev(j)
 
     DO k = 1, kpke

      IF(pddpo(j,k) > 0.5_wp) then

       surface_height = MERGE(za(j), 0._wp, k==1) 

       !=====PLANKTON DYNAMICS: ========================
       !===GROWTH
       !===GRAZING
       !===DYING
       !===(EXUDATION)
       !===(EXCRETION)

       avphy  = MAX(phytomi,bgctra(j,k,iphy))                    ! 'available' phytoplankton
       avgra  = MAX(grami, bgctra(j,k,izoo))                     ! 'available' zooplankton
       avsil  = MAX(0._wp, bgctra(j,k,isilica))                  ! available silicate
       avanut = MAX(0._wp, MIN(bgctra(j,k,iphosph),        &     ! available nutrients (phosphate   [kmol P /m3]
 &          rnoi*bgctra(j,k,iano3)))                 !                     + nitrate)
       avanfe = MAX(0._wp, MIN(avanut,bgctra(j,k,iiron)/riron))  ! available iron

       ! Nutrient limitation diagnostic
       !   Iron: flim=merge(1,0,(P<=N).and.(Fe<P))
        bgctend(j,k,kflim)=merge(1._wp,0._wp,avanfe.eq.bgctra(j,k,iiron)/riron)
       !Phosphate: plim=merge(1.,0.,(flim==0.).and.(P<=N))
        bgctend(j,k,kplim)=merge(1._wp,0._wp,(bgctend(j,k,kflim).eq.0._wp).and.(bgctra(j,k,iphosph).le.rnoi*bgctra(j,k,iano3)))
        ! Nitrate: nlim=merge(1.,0.,(flim==0.).and.(N<P))
        bgctend(j,k,knlim)=merge(1._wp,0._wp,(bgctend(j,k,kflim).eq.0._wp).and.(rnoi*bgctra(j,k,iano3).lt.bgctra(j,k,iphosph)))


       ! phytoplankton growth
       phofa = pi_alpha*fPAR*strahl(j)*swr_frac(j,k)
       temfa = 0.6_wp * 1.066_wp**ptho(j,k)
       pho   = dtb*phofa*temfa/SQRT(phofa**2 + temfa**2)

       xa    = avanfe
       xn    = xa / (1._wp + pho*avphy / (xa+bkphy) )                ! bkphy = half saturation constant
       phosy = MAX(0._wp, xa-xn)                                     ! photosynthesis

       ! zooplankton growth, phy grazing
       ya    = avphy + phosy                                         ! new phytoplankton concentration before grazing
       yn    = (ya+grazra*avgra*phytomi/(avphy+bkzoo))            &  ! grazing
&         / (1._wp + grazra * avgra / (avphy + bkzoo))
       grazing = MAX(0._wp, ya-yn)                                   
       graton  = epsher * (1._wp - zinges) * grazing                 ! "grazing to (re-dissolved) nutrients"
       gratpoc = (1._wp - epsher) * grazing                          ! epsher=0.8 "grazing to POC"
       grawa   = epsher*zinges*grazing                               ! grazer 'wachstum(?)'

       ! plankton dying
       phymor    = dyphy * MAX(0._wp, (bgctra(j,k,iphy) - 2._wp * phytomi))   ! phytoplankton mortality dyphy=.008*dt
       zoothresh = MAX(0._wp, (bgctra(j,k,izoo) - 2._wp * grami))
       zoomor    = spemor*zoothresh*zoothresh                             ! zooplankton mortality
       excdoc    = gammaz*zoothresh                                       ! excretion to DOC (zooplankton)
       exud      = gammap * MAX(0._wp, (bgctra(j,k,iphy) - 2._wp*phytomi))    ! exudation to DOC (phytoplankton)
   
       bgctra(j,k,iphosph) = bgctra(j,k,iphosph) &
                   &    - phosy + graton  &
                   &    + ecan*zoomor 

       bgctra(j,k,iano3) = bgctra(j,k,iano3) &
        &        + (-phosy+graton+ecan*zoomor)*rnit
        

       export = zoomor * (1._wp - ecan) + phymor + gratpoc               ! ecan=.95, gratpoc= .2*grazing [P-units]

       bgctra(j,k,idet) = bgctra(j,k,idet) + export                   ! 

       !=====SHELL PRODUCTIOIN ========================

       delsil = MIN(ropal * export * avsil / (avsil + bkopal), 0.5_wp * avsil)
       delcar = MIN(calmax * rcar *export, rcalc * export * bkopal/(avsil+bkopal))  ! 'detritus linked calcium carbonate '

       bgctra(j,k,isco212) = bgctra(j,k,isco212) - delcar         &                ! - CACO3 production
                   &  + rcar*(  - phosy + graton + ecan*zoomor)   ! + remineralization C-units


       bgctra(j,k,ialkali) = bgctra(j,k,ialkali) - 2._wp * delcar &
                      - rnit * ( - phosy + graton + ecan*zoomor)

       bgctra(j,k,iphy) = bgctra(j,k,iphy) + phosy - grazing       &
                   &             - phymor - exud
    
       bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)                   &
                   &   + ro2ut*phosy                    &
                   &   - (graton + ecan*zoomor)*ro2ut

      

       bgctra(j,k,izoo)  = bgctra(j,k,izoo) + grawa - excdoc - zoomor

       bgctra(j,k,idoc)  = bgctra(j,k,idoc) + excdoc + exud

       bgctra(j,k,icalc) = bgctra(j,k,icalc) + delcar


       opalrem = dremopal * 0.1_wp * (ptho(j,k) + 3.0_wp) * bgctra(j,k,iopal)

       bgctra(j,k,isilica) = bgctra(j,k,isilica) - delsil + opalrem

       bgctra(j,k,iopal) = bgctra(j,k,iopal) + delsil - opalrem

       bgctra(j,k,iiron) = bgctra(j,k,iiron) + (- phosy + graton + ecan*zoomor)*riron  &
                    &   - relaxfe * MAX(bgctra(j,k,iiron) - fesoly, 0._wp)


       bgctend(j,k,kphosy) = phosy / dtbgc 
       bgctend(j,k,kgraz) = grazing / dtbgc 
       bgctend(j,k,kgraton) = graton / dtbgc 
       bgctend(j,k,kexudp) = exud / dtbgc 
       bgctend(j,k,kexudz) = excdoc / dtbgc 
       bgctend(j,k,kzdy) = zoomor / dtbgc 
       bgctend(j,k,kpdy) =  phymor / dtbgc 
       bgctend(j,k,kdelsil) = delsil / dtbgc 
       bgctend(j,k,kdelcar) = delcar / dtbgc 
       bgctend(j,k,keuexp) = export / dtbgc 


      !===== DMS ===

       dms_prod = (dmsp(5)*delsil+dmsp(4)*delcar)                        &             ! production
         &  * (1._wp + 1._wp / (ptho(j,k) + dmsp(1))**2)

       dms_bac = dmsp(3) * dtb * ABS(ptho(j,k) + 3._wp) * bgctra(j,k,idms) & ! bacterial consumption
                   &  * (bgctra(j,k,idms)/(dmsp(6)+bgctra(j,k,idms)))

       dms_uv  = dmsp(2) * 4._wp * dtb * phofa * bgctra(j, k, idms)         ! decay due to UV-radiation


       bgctra(j,k,idms) = bgctra(j,k,idms)                      &
                   &             + dms_prod - dms_bac - dms_uv

  
       bgctend(j,k,kdmsprod) = dms_prod / dtbgc 
       bgctend(j,k,kdmsbac)  = dms_bac / dtbgc 
       bgctend(j,k,kdmsuv)   = dms_uv / dtbgc 
      

       avoxy = MAX(0._wp,bgctra(j,k,ioxygen))       ! available O2                       

       IF (avoxy >= thresh_aerob) THEN                          
           !=====AEROB REMINERALIZATION ========================
           tpremin = remido*bgctra(j,k,idoc)  &
              &         + rcyano *remido_cya*bgctra(j,k,idoccya) &
              &         + drempoc * bgctra(j,k,idet)
 
           ! Indiv. fractions of decomp.
           r_bacfra = remido*bgctra(j,k,idoc) / tpremin 
           r_bacfra_cya = rcyano * remido_cya*bgctra(j,k,idoccya) / tpremin
           r_remin  =  drempoc * bgctra(j,k,idet) / tpremin
        
           ! DOC decomposition
           xn=bgctra(j,k,idoc)/(1._wp+remido)
           bacfra=MAX(0._wp,bgctra(j,k,idoc) - xn)
           bacfra    = MIN(bacfra,   &                         !
               &     r_bacfra*(avoxy-thresh_aerob)/ro2ut)

           ! DOCcya decomposition
           xn=bgctra(j,k,idoccya)/(1._wp+remido_cya)
           bacfra_cya=MAX(0._wp,bgctra(j,k,idoccya) - xn )
           bacfra_cya = rcyano * MIN(bacfra_cya, &
              &     r_bacfra_cya*(avoxy-thresh_aerob)/ro2ut)
   
           ! POC decomposition
           xn=bgctra(j,k,idet)/(1._wp+drempoc)
           remin=MAX(0._wp,bgctra(j,k,idet)-xn)
           remin = MIN(remin,               &
                    &      r_remin*(avoxy-thresh_aerob)/ro2ut)      
          
           bgctra(j,k,idoc)  = bgctra(j,k,idoc) - bacfra 
       
           bgctra(j,k,idoccya) = bgctra(j,k,idoccya) - bacfra_cya

           bgctra(j,k,idet)  = bgctra(j,k,idet) - remin

           bgctra(j,k,iiron)  = bgctra(j,k,iiron) + (bacfra+bacfra_cya+remin)*riron  

           bgctra(j,k,iphosph) =  bgctra(j,k,iphosph) + bacfra +remin                    &
        &                     + bacfra_cya

           

           bgctra(j,k,isco212) = bgctra(j,k,isco212)             &  
        &                + rcar*( bacfra + bacfra_cya + remin)                       ! + remineralization C-units


           bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)               &
        &            -(bacfra + bacfra_cya +remin)*ro2ut     

          
           bgctra(j,k,ialkali) = bgctra(j,k,ialkali)       &
        &             - (bacfra + bacfra_cya + remin)*rnit

           bgctra(j,k,iano3)= bgctra(j,k,iano3)                  &
        &             +(bacfra + bacfra_cya +remin)*rnit


           !***********************************************************************
           !   There is about 1.e4 O2 on 1 N2O molecule (Broecker&Peng)
           !    refra : Tim Rixen, pers. communication
           !***********************************************************************

           aou   = satoxy(j,k) - bgctra(j,k,ioxygen)
           refra = 1._wp + 3._wp * (0.5_wp + SIGN(0.5_wp, aou - 1.97e-4_wp))

           bgctra(j,k,ian2o)   = bgctra(j,k,ian2o)                  &
                   &    + (remin + bacfra + bacfra_cya) * 1.e-4_wp * ro2ut * refra

           bgctra(j,k,igasnit) = bgctra(j,k,igasnit)                &
                   &     - (remin + bacfra + bacfra_cya) * 1.e-4_wp * ro2ut * refra

          

           bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)               &
                   &     - (remin+bacfra+bacfra_cya)*1.e-4_wp*ro2ut*refra*0.5_wp
         
           bgctend(j,k,kremin) = remin / dtbgc 
           bgctend(j,k,kbacfra) = bacfra / dtbgc 
           bgctend(j,k,kbacfrac) = bacfra_cya / dtbgc 
           bgctend(j,k,kdenit) = 0._wp 
       ENDIF   ! O2 >= thresh_aerob

       avoxy = MAX(0._wp,bgctra(j,k,ioxygen))       ! available O2                       

       IF (avoxy < thresh_aerob) THEN                          
            !=====DENITRIFICATION ========================

           ! NO3 reduction
           xa = MIN(bgctra(j,k,idet),  &
                    0.5_wp * bgctra(j,k,iano3) / nitdem)

           xn = xa / (1._wp + denitrification*drempoc)
         
           remin = MAX(0._wp,xa-xn)
          ! remin = denitrification*drempoc * MIN(bgctra(j,k,idet),  &
          !          0.5_wp * bgctra(j,k,iano3) / nitdem)

           bgctra(j,k,idet)      = bgctra(j,k,idet)      &
                   &                -  remin
           ! N2O reduction
           xa =  MIN(bgctra(j,k,idet),                    &  ! remineralization using N2O
                   &                  0.003_wp * bgctra(j,k,ian2o) / (2._wp*ro2ut))

           xn = xa/(1._wp + dremn2o)

           remin2o = MAX(0._wp, xa-xn)
          ! remin2o = dremn2o * MIN(bgctra(j,k,idet),                    &  ! remineralization using N2O
          !         &                  0.003_wp * bgctra(j,k,ian2o) / (2._wp*ro2ut))

           bgctra(j,k,idet)      = bgctra(j,k,idet)      &
                   &                -  remin2o
     
           bgctra(j,k,iphosph) = bgctra(j,k,iphosph) + remin + remin2o 

           bgctra(j,k,isco212)  = bgctra(j,k,isco212)  &
                   &                 + rcar*(remin + remin2o)



           bgctra(j,k,iano3)   = bgctra(j,k,iano3)     &
                   &          - nitdem*remin + rnit*remin2o

           bgctra(j,k,igasnit) = bgctra(j,k,igasnit)   &
                   &                + n2prod*remin + 2._wp*ro2ut*remin2o
          

           bgctra(j,k,iiron)   = bgctra(j,k,iiron)     &
                   &          + riron*(remin+remin2o)

           bgctra(j,k,ian2o)   = bgctra(j,k,ian2o)     &
                   &         - 2._wp * ro2ut * remin2o

           !alkalinity is increased during denitrification due to consumption of H+ (see Wolf-Gladrow etal,2007)
           bgctra(j,k,ialkali) = bgctra(j,k,ialkali)   &
                   &          + nitdem*remin - rnit*remin2o

          bgctend(j,k,kn2b) = bgctend(j,k,kn2b) + 2._wp * n2prod * remin * (pddpo(j,k) + surface_height) 
          !denitrification produces water (H2O), the corresponding O2 uptake is budgeted in h2obudget
          bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob) + 0.5_wp * n2prod * remin * (pddpo(j,k) + surface_height) 


           bgctend(j,k,kdenit) = remin / dtbgc 
           bgctend(j,k,kremin) = 0._wp 
           bgctend(j,k,kbacfra) = 0._wp 
           bgctend(j,k,kbacfrac) = 0._wp 

       ENDIF ! oxygen < thresh_aerob


       avoxy = MAX(0._wp,bgctra(j,k,ioxygen))       ! available O2                       
       bgctend(j,k,ksred) = 0._wp 

       !=====SULFATE REDUCTION ========================
       IF (bgctra(j,k, ioxygen) < thresh_sred) THEN

             xn=bgctra(j,k,idet)/(1._wp+sulfate_reduction)
             remin=MAX(0._wp,bgctra(j,k,idet)-xn)
 
             bgctra(j,k,idet)    = bgctra(j,k,idet)    -        remin
             bgctra(j,k,ialkali) = bgctra(j,k,ialkali) + rnit * remin 
             bgctra(j,k,isco212) = bgctra(j,k,isco212) + rcar * remin
             bgctra(j,k,iphosph) = bgctra(j,k,iphosph) +        remin

             bgctra(j,k,iano3)   = bgctra(j,k,iano3)   + rnit  * remin
             bgctra(j,k,iiron)   = bgctra(j,k,iiron)   + riron * remin
             
             bgctend(j,k,ksred) = remin / dtbgc 
     
             bgctend(j,k,kn2b) = bgctend(j,k,kn2b) + 2._wp * rnit * remin * (pddpo(j,k) + surface_height) 
             bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob) - ro2ut * remin * (pddpo(j,k) + surface_height) 

       ENDIF ! O2 < thresh_sred
       bgctend(j,k,kaou)   = satoxy(j,k) - bgctra(j,k,ioxygen)
  
      ENDIF ! wet cells
     ENDDO ! k=1,kpke
  ENDDO ! j=start_idx,end_idx 

!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  END SUBROUTINE ocprod
