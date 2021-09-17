  !! @file ocprod.f90
  !! @brief Computes plankton dynamics, OM degradation
  !!
#include "hamocc_omp_definitions.inc"

MODULE mo_ocprod

  USE mo_kind, ONLY           : wp
   
  USE mo_memory_bgc, ONLY     :phytomi, grami, rnoi, riron, pi_alpha, &
       &                        fpar, bkphy, bkzoo, epsher,         &
       &                        zinges, ro2ut, remido, dyphy, spemor,     &
       &                        gammaz, gammap, ecan, rnit, ropal, bkopal,         &
       &                        rcar, relaxfe, fesoly,            &
       &                        denitrification, nitdem, dremn2o,         &
       &                        n2prod, sulfate_reduction, strahl,                 &
       &                        thresh_aerob, thresh_o2, prodn2o, & 
       &                        thresh_sred, dmsp, &
       &                        satoxy, meanswr, ralk, bkh2sox, rh2sox,&
       &                        bgctra, bgctend,docmin, &
       &                        bkpo4, bkfe, bkno3, bknh4, ro2ammo, no2denit, &
       &                        anamoxra, bkno2, bkrad, nitriox, nitrira, bkfe, o2thresh, rno3no2, &
       &                        rno3nh4, rnh4no2, rno2no3, alk_nrn2, rno2n2, o2den_lim, swr_frac


  USE mo_control_bgc, ONLY    : dtb, bgc_nproma, bgc_zlevs, dtbgc 
  USE mo_param1_bgc, ONLY     : icalc, iopal, ian2o, igasnit, idms, &
       &                        iphy, izoo, isilica, iphosph, &
       &                        iano3, ioxygen, idet, idoc, isco212, &
       &                        ialkali, kphosy, kremin, ih2s,    &
       &                        iiron, ksred, kdenit, kgraz, kbacfra, &
       &                        kn2b, kh2ob, kdelcar, kdelsil, kbacfrac, &
       &                        kdmsprod, kdmsuv, kdmsbac,keuexp, knlim, kflim,&
       &                        kplim, kgraton, kexudp, kexudz, kpdy,kzdy,kaou,&
       &                        kh2sprod, kh2sloss, &
       &                        iammo, iano2, &
       &                        kgppnh, kdnrn, kdnra, kanam, kammox, knitox

    USE mo_hamocc_nml, ONLY    : grazra, calmax, dremopal, drempoc, &
       &                         l_N_cycle, no3nh4red, no3no2red 
  PUBLIC :: ocprod
           

CONTAINS


SUBROUTINE ocprod (klev,start_idx,end_idx, ptho,pddpo, za,ptiestu, l_dynamic_pi)
    

  IMPLICIT NONE

  ! Arguments

  INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 
  REAL(wp), INTENT(in), TARGET   :: ptho(bgc_nproma,bgc_zlevs)       !<  potential temperature (degC)
  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]

  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height
  REAL(wp), INTENT(in) :: ptiestu(bgc_nproma,bgc_zlevs) !< depth of scalar grid cell [m]

  LOGICAL, INTENT(in) :: l_dynamic_pi

 !  Local variables

  INTEGER                :: j,k, kpke
  REAL(wp) :: avphy,avanut,avanfe,xa,xn,ya,yn,phosy,                   &
       &      pho, phofa, temfa,                                         &
       &      avgra,grazing,avsil,graton,                                &
       &      gratpoc,grawa,bacfra,phymor,zoothresh,zoomor,excdoc,exud,  &
       &      export, delsil, delcar, remin,             &
       &      opalrem, remin2o, aou, refra,                      &
       &      o2lim, actn2o, avdet,   &
       &      maxn2o, avoxy, oxid

  REAL(wp) :: surface_height

  REAL(wp) :: dms_prod, dms_uv, dms_bac 

 ! N-cycle variables
  REAL(wp) :: limp, limf,limn, po4lim, felim, no3lim, nh4lim, hib
  REAL(wp) :: nfrac, detn, remin_nit, rdnrn, rdnra, fdnrn, fdnra, no3rmax
  REAL(wp) :: no3a, no3c_max, detc_max, detc_act, n2ormax
  REAL(wp) :: anam, nh4n, nh4a, annpot, fammox, fnitox, newammo
  REAL(wp) :: newnitr, oxpot, nitox, ammox, remsulf, detnew, ntotlim
  REAL(wp) :: dnrn, dnra, nlim, no2a, no2c_max, n2oa, nrn2, n2oc_max 
  REAL(wp) :: anamox, avo2, n2oprod, n2on2, no2rmax
 

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,surface_height,avphy,avgra,avsil,avanut,&
!HAMOCC_OMP            avanfe,phofa,temfa,pho,phosy,xa,xn,ya,yn,grazing,&
!HAMOCC_OMP            graton,gratpoc,grawa,phymor,zoomor,zoothresh,excdoc,&
!HAMOCC_OMP            exud,export,delsil,delcar,opalrem,dms_prod,dms_bac,&
!HAMOCC_OMP            dms_uv,avoxy,maxn2o,actn2o,o2lim, avdet, oxid,&
!HAMOCC_OMP            bacfra,remin,aou,refra,remin2o) HAMOCC_OMP_DEFAULT_SCHEDULE

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
      if(l_dynamic_pi)then
         phofa=(pi_alpha + 0.05_wp*ptiestu(j,k)/(ptiestu(j,k)+90._wp)) &
     &            *fPAR*strahl(j)*meanswr(j,k)
   
      else
       phofa = pi_alpha*fPAR*strahl(j)*meanswr(j,k)
      endif
       temfa = 0.6_wp * 1.066_wp**ptho(j,k)
       pho   = dtb*phofa*temfa/SQRT(phofa**2 + temfa**2)

       !!!! N cycle !!!!!!!!
       IF (l_N_cycle) THEN
          xa = MAX(0._wp, (bgctra(j,k,iphosph)))
          xn = xa/(1._wp + pho*avphy/(xa + bkpo4))
          po4lim = MAX(0._wp, xa - xn)
          limp = xn/(xa + bkpo4)             

          xa = MAX(0._wp, (bgctra(j,k,iiron)))
          xn = xa/(1._wp + pho*avphy*riron/(xa + bkfe))
          felim = MAX(0._wp,xa - xn)
          limf = xn/(xa + bkfe)              

          hib = 1._wp/(1._wp + bgctra(j,k,iammo)/bknh4) ! sigma_inhib

          xa = MAX(0._wp, (bgctra(j,k,iano3)))
          xn = xa/(1._wp + pho*avphy*rnit*hib/(xa + bkno3))
          no3lim = MAX(0._wp, xa - xn)
          limn = xn/(xa + bkno3)             

          xa = MAX(0._wp, (bgctra(j,k,iammo)))
          xn = xa/(1._wp + pho*avphy*rnit/(xa + bknh4))
          nh4lim = MAX(0._wp, xa - xn )       
          ntotlim = no3lim + nh4lim          
       
          phosy = MIN(po4lim, ntotlim/rnit, felim/riron)

          nfrac = 1._wp
          if(phosy .gt. 1.E-18_wp) nfrac= nh4lim/ntotlim     ! fraction of photosynthesis on NH4
          limn = limn*(1._wp - nfrac) + nfrac*xn/( xa + bknh4) 


          IF ( limf .le. limp .and. limf .le. limn) THEN
             bgctend(j,k,kflim) = 1._wp
             bgctend(j,k,knlim) = 0._wp
             bgctend(j,k,kplim) = 0._wp
          ELSE IF (limn .le. limp) THEN
             bgctend(j,k,kplim) = 0._wp
             bgctend(j,k,knlim) = 1._wp
             bgctend(j,k,kflim) = 0._wp
          ELSE
             bgctend(j,k,kplim) = 1._wp
             bgctend(j,k,knlim) = 0._wp
             bgctend(j,k,kflim) = 0._wp
          END IF

       ELSE

          xa    = avanfe
          xn    = xa / (1._wp + pho*avphy / (xa+bkphy) )                ! bkphy = half saturation constant
          phosy = MAX(0._wp, xa-xn)                                     ! photo synthesis 

       ENDIF
       !!!! N cycle !!!!!!!!

       phosy=MERGE(bgctra(j,k,isco212)/rcar,phosy,bgctra(j,k,isco212).le.rcar*phosy) ! limit phosy by available DIC
       if(pho < 0.000001_wp)phosy=0._wp 

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

       !!!! N cycle !!!!!!!!
       IF (l_N_cycle) THEN
          bgctra(j,k,iammo) =  bgctra(j,k,iammo)                            & 
                   &                + (graton + ecan*zoomor)*rnit                   & !  all remineralization products added to NH4
                   &                - nfrac*phosy*rnit
          
          bgctra(j,k,iano3) = bgctra(j,k,iano3) - (1._wp - nfrac)*phosy*rnit
       ELSE
          bgctra(j,k,iano3) = bgctra(j,k,iano3) &
        &        + (-phosy+graton+ecan*zoomor)*rnit
       ENDIF
       !!!! N cycle !!!!!!!!
        

       export = zoomor * (1._wp - ecan) + phymor + gratpoc               ! ecan=.95, gratpoc= .2*grazing [P-units]

       bgctra(j,k,idet) = bgctra(j,k,idet) + export                   ! 

       !=====SHELL PRODUCTIOIN ========================

       delsil = MIN(ropal * export * avsil / (avsil + bkopal), 0.5_wp * avsil)
       delcar = calmax * rcar *export * bkopal/(avsil+bkopal)  ! 'detritus linked calcium carbonate '

       bgctra(j,k,isco212) = bgctra(j,k,isco212) - delcar         &                ! - CACO3 production
                   &  + rcar*(  - phosy + graton + ecan*zoomor)   ! + remineralization C-units

       bgctra(j,k,iphy) = bgctra(j,k,iphy) + phosy - grazing       &
                   &             - phymor - exud

       !!!! N cycle !!!!!!!!
       IF (l_N_cycle) THEN

          ! LR: Why is rnit used here and not ralk ?!
          bgctra(j,k,ialkali) = bgctra(j,k,ialkali)                   &    ! ocean with NH4 - alkalinity change
                    &           - nfrac*phosy*rnit                      &    ! alk decrease if OM from NH4 
                    &           + (1.-nfrac)*phosy*rnit                     &    ! alk increase if OM from NO3
                    &           + rnit*(graton + ecan*zoomor)               &    ! remin all to NH4
                    &           - (graton - phosy + ecan*zoomor)              &    ! PO4 changes                      
                    &           - 2._wp*delcar

          bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)               &
                   &            + phosy*(ro2ut*(1._wp - nfrac) + ro2ammo*nfrac)       & ! phosy from NO3 produces ro2ut, from NH4 only ro2ammo
                   &            - (graton + ecan*zoomor)*ro2ammo                      ! since all re
 
       ELSE

          bgctra(j,k,ialkali) = bgctra(j,k,ialkali) - 2._wp * delcar &
                               - ralk * ( - phosy + graton + ecan*zoomor)

          bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)                   &
                            &   + ro2ut*phosy                    &
                            &   - (graton + ecan*zoomor)*ro2ut

       ENDIF
       !!!! N cycle !!!!!!!!     

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
       if (l_N_cycle) bgctend(j,k,kgppnh) = nfrac*phosy / dtbgc


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
      


       IF (bgctra(j,k,ioxygen) > thresh_aerob) THEN                          

           !=====AEROB REMINERALIZATION ========================
           avoxy = bgctra(j,k,ioxygen) - thresh_aerob      ! available O2                       
       
           o2lim = bgctra(j,k,ioxygen)/(thresh_o2 + bgctra(j,k,ioxygen))

           ! POC decomposition
           xn=bgctra(j,k,idet)/(1._wp+o2lim*drempoc)
           remin=MAX(0._wp,bgctra(j,k,idet)-xn)

           !!!! N cycle !!!!!!!!
           IF (l_N_cycle) THEN
              remin = MIN(remin, 0.5_wp*avoxy/ro2ammo)
           ELSE
              remin = MIN(remin, avoxy/ro2ut)      
           ENDIF
           !!!! N cycle !!!!!!!!

           ! DOC decomposition
           avoxy = bgctra(j,k,ioxygen) -remin*ro2ut - thresh_aerob      ! available O2                       
           xn=(bgctra(j,k,idoc)+remido*docmin)/(1._wp+remido)
           bacfra=MAX(0._wp,bgctra(j,k,idoc) - xn)

           !!!! N cycle !!!!!!!!
           IF (l_N_cycle) THEN
              bacfra = MERGE(-0._wp,bacfra, avoxy-bacfra*ro2ammo.lt.thresh_aerob) 
           ELSE
              bacfra = MERGE(-0._wp,bacfra, avoxy-bacfra*ro2ut.lt.thresh_aerob)
           ENDIF
           !!!! N cycle !!!!!!!!

           bgctra(j,k,idoc)  = bgctra(j,k,idoc) - bacfra 
       

           bgctra(j,k,idet)  = bgctra(j,k,idet) - remin

           bgctra(j,k,iiron)  = bgctra(j,k,iiron) + (bacfra+remin)*riron  

           bgctra(j,k,iphosph) =  bgctra(j,k,iphosph) + bacfra +remin                    

           

           bgctra(j,k,isco212) = bgctra(j,k,isco212)             &  
        &                + rcar*( bacfra +  remin)                       ! + remineralization C-units


           !!!! N cycle !!!!!!!!
           IF (l_N_cycle) THEN

              bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)             &
                         &            - bacfra*ro2ammo            & ! since all remineralization products go to NH4
                         &            - remin*ro2ammo               ! there is less oxygen demand ro2ammo=(ro2ut - 2*rnit) 

              ! LR: Why is rnit used here and not ralk ?!
              ! This is strange: why no prefactor for "PO4 update" ?
              bgctra(j,k,ialkali) = bgctra(j,k,ialkali)            &
                         &            -    (remin + bacfra)       & ! PO4 update
                         &            + rnit*(bacfra + remin)       ! NH

              bgctra(j,k,iammo)= bgctra(j,k,iammo)                 & !  all remineralization products in EU added to NH4
                         &            + (bacfra + remin)*rnit

           ELSE
       
              bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)       &
                        &            -(bacfra + remin)*ro2ut     

          
              bgctra(j,k,ialkali) = bgctra(j,k,ialkali)       &
                        &            -(bacfra +  remin)*ralk

              bgctra(j,k,iano3)= bgctra(j,k,iano3)            &
                        &            +(bacfra + remin)*rnit
           ENDIF
           !!!! N cycle !!!!!!!!


           !***********************************************************************
           !   There is about 1.e4 O2 on 1 N2O molecule (Broecker&Peng)
           !    refra : Tim Rixen, pers. communication
           !***********************************************************************

           aou   = satoxy(j,k) - bgctra(j,k,ioxygen)
           refra = 1._wp + 3._wp * (0.5_wp + SIGN(0.5_wp, aou - 1.97e-4_wp))
       

           maxn2o = (remin+bacfra)*prodn2o*ro2ut*refra*0.5_wp
           avoxy = max(0._wp, bgctra(j,k,ioxygen)-thresh_aerob)
           actn2o = min(avoxy,maxn2o)

           bgctra(j,k,ian2o)   = bgctra(j,k,ian2o)                  &
                   &    + 2._wp *actn2o

           bgctra(j,k,igasnit) = bgctra(j,k,igasnit)                &
                   &     - 2._wp *actn2o

           bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)               &
                   &     - actn2o
         
           bgctend(j,k,kremin) = remin / dtbgc 
           bgctend(j,k,kbacfra) = bacfra / dtbgc 
           bgctend(j,k,kdenit) = 0._wp 
       ENDIF   ! O2 >= thresh_aerob


       !!!! N-cycle !!!!!!!!
       IF (l_N_cycle) THEN
       IF (bgctra(j,k,idet) > 1.e-15_wp) THEN
       IF (bgctra(j,k,ioxygen) <= o2thresh) THEN  ! o2thresh = 1e-5 currently

          ! o2 limitation identical for all suboxic processes
          ! LR: this seems to be there for avoiding cases where oxygen is
          !     negative?
          o2lim = min(1._wp,1._wp - bgctra(j,k,ioxygen)/o2thresh)

          ! convert detritus in P-units to N-units for nitrogen cycle changes  
          detn = max(0._wp, bgctra(j,k,idet)*rnit)
          remin_nit = 0._wp

          ! denitrification on NO3
          rdnrn = o2lim*no3no2red *detn/(bgctra(j,k,iano3) + 3.E-5_wp)
          rdnra = o2lim*no3nh4red *detn/(bgctra(j,k,iano3) + 3.E-5_wp)

          no3rmax = rdnrn + rdnra                       ! max pot loss of NO3

          ! LR: no3rmax < 0 only occurs, if no3 is strongly negative?!
          IF (no3rmax > 0._wp) THEN
          
          fdnrn = rdnrn/no3rmax                       ! fraction each process
          fdnra = rdnra/no3rmax

          !< implicit formulation to avoid neg. nitrate concentration
          no3a = bgctra(j,k,iano3)/(1._wp + no3rmax)   ! max change in NO3  
          no3c_max = bgctra(j,k,iano3) - no3a         ! corresponding max NO3 loss 
          detc_max =  no3c_max*(fdnrn/rno3no2 + fdnra/rno3nh4)  ! corresponding max change in det 
          detc_act = min (0.9_wp*bgctra(j,k,idet), detc_max)

          dnrn = fdnrn*detc_act             ! in P units
          dnra = fdnra*detc_act             ! in P untis 

          remin_nit = dnrn + dnra      ! change for DIC and PO4

          bgctend(j,k,kdnrn) =  dnrn*rno3no2 / dtbgc    ! in N units
          bgctend(j,k,kdnra) =  dnra*rno3nh4 / dtbgc    ! in N units
          ! LR: rate diagnostics not yet implemented
          ! bgc_o_pro(i,j,k,kradnrn) = rdnrn*86400._wp*inv_dtbgc    ! in per day
          ! bgc_o_pro(i,j,k,kradnra) = rdnra*86400._wp*inv_dtbgc

          ! change from dnrn and dnra in other tracers 
          bgctra(j,k,idet) = bgctra(j,k,idet) - dnrn - dnra           ! change in detritus                        

          bgctra(j,k,iano3) = bgctra(j,k,iano3)    &    ! change in nitrate 
                        &     - rno3no2*dnrn       &    ! from DNRN
                        &     - rno3nh4*dnra            ! from DNRA

          bgctra(j,k,iammo) = bgctra(j,k,iammo)    &    ! change in ammonium 
                      &       + rnit*dnrn          &    ! from DNRN
                      &       + 86._wp*dnra             ! from DNRA

          bgctra(j,k,iano2) = bgctra(j,k,iano2)    &    ! change in nitrite 
                      &       + rno3no2*dnrn            ! from DNRN

          ! LR: Why rnit?!
          !     Where does the 86 come from and why is it used for changes in 
          !     ammo and alkali, even though units should be different?
          bgctra(j,k,ialkali) = bgctra(j,k,ialkali)  &   ! change from DNRN and DNRA
                                &          + rnit*dnrn             &   ! from DNRN
                                &          + 86._wp*dnra           &   ! from DNRA
                                &          - remin_nit                 ! PO4 update

          bgctend(j,k,kn2b) = bgctend(j,k,kn2b) - 70._wp*dnra*(pddpo(j,k) + surface_height)


          bgctra(j,k,iphosph) = bgctra(j,k,iphosph) + remin_nit
          bgctra(j,k,isco212) = bgctra(j,k,isco212) + rcar*remin_nit
          bgctra(j,k,iiron)   = bgctra(j,k,iiron)  + riron*remin_nit 

          ENDIF ! no3rmax > 0.0

          ! allow complete denitrification only at very low O2
          ! denitrification  on NO2
          if (bgctra(j,k,ioxygen) < o2den_lim) then
             detn = max(0._wp,bgctra(j,k,idet)*rnit)
             o2lim = 1._wp - bgctra(j,k,ioxygen)/o2thresh
             no2rmax = o2lim*no2denit*detn/(bgctra(j,k,iano2) + 0.1E-6_wp)

             ! implicit formulation to avoid neg. nitrite concentration
             no2a = bgctra(j,k,iano2)/(1._wp + no2rmax) 
             no2c_max = max(0._wp,bgctra(j,k,iano2) - no2a)         ! corresponding max NO2 loss 
             detc_max=  no2c_max/rno2n2                    ! corresponding max change in det, rno2n2 conversion to P 
             nrn2 = min (bgctra(j,k,idet),detc_max)   ! in P units
             nrn2 = max(0._wp,nrn2)

             ! change from nrn2 in other tracers 
             bgctra(j,k,idet) = bgctra(j,k,idet) - nrn2          ! change in detritus
             bgctra(j,k,iano2) = bgctra(j,k,iano2) - rno2n2*nrn2 ! change in nitrite      
             bgctra(j,k,iammo) = bgctra(j,k,iammo) + rnit*nrn2   ! change in ammonium 

             bgctra(j,k,igasnit) = bgctra(j,k,igasnit) & 
                           &       + rno2n2*nrn2/2._wp           ! from nitrite reduction      

             bgctra(j,k,ialkali) = bgctra(j,k,ialkali)  &  
                           &       + alk_nrn2*nrn2      &
                           &       - nrn2                   

             bgctend(j,k,kn2b) = bgctend(j,k,kn2b)       &
                           &       + (alk_nrn2 - rnit)*nrn2*(pddpo(j,k) + surface_height)


             bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob) + rno2n2*nrn2*0.25_wp         &
                           &      * (pddpo(j,k) + surface_height)  ! produces losses 560/3 o2

             bgctra(j,k,iphosph) = bgctra(j,k,iphosph) + nrn2
             bgctra(j,k,isco212) = bgctra(j,k,isco212) + rcar*nrn2
             bgctra(j,k,iiron)  = bgctra(j,k,iiron)  + riron*nrn2

             ! LR: not yet implemented
             ! bgcprod(i,j,k,kremnit) = bgcprod(i,j,k,kremnit) + nrn2*inv_dtbgc ! in P -units
             !bgc_o_pro(i,j,k,kradeni) = no2rmax*86400._wp*inv_dtbgc

             bgctend(j,k,kdenit)  = rno2n2*nrn2 / dtbgc  ! in N units

             ! denitrification using laughing gas
             ! with correct stoichometric
             detn = max(0._wp,bgctra(j,k,idet)*rnit)

             n2ormax = o2lim*dremn2o*detn/(bgctra(j,k,ian2o)+5.E-9_wp)

             n2oa = bgctra(j,k,ian2o)/(1._wp + n2ormax) 
             n2oc_max = max(0._wp, bgctra(j,k,ian2o) - n2oa)         ! corresponding max N2O loss 
             detc_max =  n2oc_max/280._wp                    ! corresponding max change in det, rno2n2 conversion to P 
             n2on2 = min (bgctra(j,k,idet), detc_max)   ! in P units
             n2on2 = max(0._wp, n2on2)   

             ! change from nrn2 in other tracers excl. DIC and PO4, done later
             bgctra(j,k,idet) = bgctra(j,k,idet) - n2on2              ! change in detritus
             bgctra(j,k,ian2o) = bgctra(j,k,ian2o) - 280._wp*n2on2    ! change in nitrous oxide
             bgctra(j,k,iammo) = bgctra(j,k,iammo) + rnit*n2on2       ! change in ammonium 

             bgctra(j,k,igasnit) = bgctra(j,k,igasnit) & 
                           &       + 280._wp*n2on2                    ! dinitrogen production

             ! LR: again, why rnit?
             bgctra(j,k,ialkali) = bgctra(j,k,ialkali)  &  
                           &       + rnit*n2on2                &
                           &       - n2on2                     

             bgctra(j,k,iphosph) = bgctra(j,k,iphosph) + n2on2
             bgctra(j,k,isco212) = bgctra(j,k,isco212) + rcar*n2on2
             bgctra(j,k,iiron)   = bgctra(j,k,iiron)  + riron*n2on2

             ! LR: not yet implemented
             !bgcprod(i,j,k,kremn2o) = n2on2*inv_dtbgc ! in P -units
          else
             bgctend(j,k,kdenit) = 0._wp
          end if ! O2 < o2denlim   
       end if ! O2 < o2thresh   
       ENDIF ! det > 1.e-15

       ELSE ! no extendend N-cycle

       IF (bgctra(j,k,ioxygen) <= o2den_lim) THEN                          
            !=====DENITRIFICATION ========================

           avdet = MAX(1.e-15_wp,bgctra(j,k,idet))

           ! NO3 reduction
           o2lim = 1._wp - bgctra(j,k,ioxygen)/thresh_o2

           remin  = denitrification*drempoc*o2lim  &
                   & * MIN(avdet,0.5_wp*bgctra(j,k,iano3)/nitdem)

           bgctra(j,k,idet)      = bgctra(j,k,idet)      &
                   &                -  remin
           ! N2O reduction

           remin2o = dremn2o * MIN(avdet,                    &  ! remineralization using N2O
                   &   0.003_wp * bgctra(j,k,ian2o) / (2._wp*ro2ut))


           bgctra(j,k,idet)      = bgctra(j,k,idet)      &
                   &                -  remin2o
     
           bgctra(j,k,iphosph) = bgctra(j,k,iphosph) + remin + remin2o 

           bgctra(j,k,isco212)  = bgctra(j,k,isco212)  &
                   &                 + rcar*(remin + remin2o)



           bgctra(j,k,iano3)   = bgctra(j,k,iano3)     &
                   &          - nitdem*remin + rnit*(remin + remin2o)

           bgctra(j,k,igasnit) = bgctra(j,k,igasnit)   &
                   &                + n2prod*remin + 2._wp*ro2ut*remin2o
          

           bgctra(j,k,iiron)   = bgctra(j,k,iiron)     &
                   &          + riron*(remin+remin2o)

           bgctra(j,k,ian2o)   = bgctra(j,k,ian2o)     &
                   &         - 2._wp * ro2ut * remin2o

           !alkalinity is increased during denitrification due to consumption of H+ (see Wolf-Gladrow etal,2007)
           bgctra(j,k,ialkali) = bgctra(j,k,ialkali)   &
                   &          + nitdem*remin - ralk*(remin+remin2o)

          bgctend(j,k,kn2b) = bgctend(j,k,kn2b) + nitdem * remin * (pddpo(j,k) + surface_height) 
          !denitrification produces water (H2O), the corresponding O2 uptake is budgeted in h2obudget
          bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob) + 0.5_wp * n2prod * remin * (pddpo(j,k) + surface_height) 


           bgctend(j,k,kdenit) = remin / dtbgc 
           bgctend(j,k,kremin) = 0._wp 
           bgctend(j,k,kbacfra) = 0._wp 
           bgctend(j,k,kbacfrac) = 0._wp 

       ENDIF ! oxygen < thresh_aerob
       ENDIF !!!! N-cycle !!!!!!!!


       !!!! N-cycle !!!!!!!!
       IF (l_N_cycle) THEN

          !!!! Anammox
          IF(bgctra(j,k,ioxygen) < o2thresh) THEN       ! in suboxic water

             ! anamox on NO2 and NH4
             nh4a = max(0._wp,bgctra(j,k,iammo) - 1.E-15_wp)
             no2a = max(0._wp,bgctra(j,k,iano2) - 1.E-15_wp)

             o2lim = 1._wp -max(0._wp,bgctra(j,k,ioxygen)/o2thresh)
             anam = o2lim*anamoxra*no2a/(no2a+bkno2)   

             nh4n = nh4a/(1._wp + anam)
             annpot = nh4a - nh4n            ! potential change in NH4 due to anamox
             anamox = min(annpot, no2a/1.3_wp)

             bgctra(j,k,iammo) = bgctra(j,k,iammo) - anamox           ! from anamox               
             bgctra(j,k,iano2) = bgctra(j,k,iano2) - 1.3_wp*anamox    ! from anamox               
             bgctra(j,k,igasnit) = bgctra(j,k,igasnit) + anamox       ! from anamox, gasnit in N2 
             bgctra(j,k,iano3) = bgctra(j,k,iano3) + 0.3_wp*anamox    ! change in nitrite 
             bgctra(j,k,ialkali) = bgctra(j,k,ialkali) - 0.3_wp*anamox

             ! loss of Os from NO2 - gain NO3 - 0.5 from NH4 =1.3 - 0.3*1.5- 0.5 = +0.35
             bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob)          &
                           &       + anamox*0.35_wp * (pddpo(j,k) + surface_height)          

             bgctend(j,k,kn2b) = bgctend(j,k,kn2b)  + 1.7_wp*anamox  &
                          &      * (pddpo(j,k) + surface_height)  ! alk is changed for NO3, but for change 
                                                               ! in NH4 we need n2bugdet change 

             bgctend(j,k,kanam)  = 2._wp*anamox / dtbgc !loss to N2 from NH4 and NO2, in kmolN /m3 s 

             ! LR: not yet implemented  
             ! bgc_o_pro(i,j,k,kraanam) = anam*86400._wp*inv_dtbgc
          ELSE
             bgctend(j,k,kanam)  = 0._wp
             ! bgc_o_pro(i,j,k,kraanam)= 0._wp
          ENDIF ! O2 < 0.5


          !!!! Oxidation of NH4 and NO2
          fammox = nitrira*bkrad/(bkrad + swr_frac(j,k)*strahl(j)) ! local rate of ammox per time step
          fnitox = nitriox*bkrad/(bkrad + swr_frac(j,k)*strahl(j)) ! local rate of nitrite ox

          ! LR: rates not yet implemented
          !bgc_o_pro(i,j,k,kraammox) = fammox*86400._wp*inv_dtbgc
          !bgc_o_pro(i,j,k,kranitox) = fnitox*86400._wp*inv_dtbgc

          nh4a  =  bgctra(j,k,iammo)                         
          newammo = nh4a/(1._wp + fammox)

          no2a =  max(0._wp, bgctra(j,k,iano2))  
          newnitr = no2a/ (1._wp+fnitox)                        ! change of nitrite

          oxpot = max(0._wp, (bgctra(j,k,ioxygen) - 0.5E-6_wp)/rno2no3)     ! max change in o2
          nitox = min(no2a - newnitr, oxpot)

          oxpot =  max(0._wp, (bgctra(j,k,ioxygen) - 0.5E-6_wp - rno2no3*nitox)/rnh4no2)
          ammox = min(nh4a - newammo, oxpot)

          bgctra(j,k,iammo) = nh4a - ammox
          bgctra(j,k,iano2) = bgctra(j,k,iano2) + ammox - nitox ! change of nitrite
          bgctra(j,k,iano3) = bgctra(j,k,iano3) + nitox
          bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen)               &
                   &            - rno2no3*nitox - rnh4no2*ammox        ! O2 will be used during nitrification 

          bgctra(j,k,ialkali) = bgctra(j,k,ialkali)               &    ! ocean with NH4 - alkalinity change
                    &           - 2._wp*ammox              ! according to Wolf-Gladrow etal Mar. Chem. (2007)

          bgctend(j,k,kammox) =  ammox/dtbgc                    
          bgctend(j,k,knitox) =  nitox/dtbgc                    

          ! N2O production relies on nitrification on ammonimum
          aou   = satoxy(j,k)-bgctra(j,k,ioxygen)
          refra = 1._wp + 3._wp * (0.5_wp + SIGN(0.5_wp, aou - 1.97e-4_wp))
          maxn2o = rnh4no2*ammox * prodn2o * refra
          avo2 = max(0._wp,bgctra(j,k,ioxygen)- 0.5e-6_wp)
          n2oprod =  min(avo2,maxn2o)
       
          bgctra(j,k,ian2o)   = bgctra(j,k,ian2o) + n2oprod
          bgctra(j,k,igasnit) = bgctra(j,k,igasnit) - n2oprod
          ! divide by 2 since it consumes 1 mol oxygen which is given in mol O2 
          bgctra(j,k,ioxygen) = bgctra(j,k,ioxygen) - n2oprod *0.5_wp

       ENDIF
       !!!! N-cycle !!!!!!!!


       bgctend(j,k,ksred) = 0._wp 

       !=====SULFATE REDUCTION ========================
       IF (bgctra(j,k, ioxygen) < thresh_sred) THEN

             IF (l_N_cycle) THEN

                o2lim = 1._wp - bgctra(j,k,ioxygen)/o2thresh
                nlim = max(0._wp,1._wp - (bgctra(j,k,iano3) + bgctra(j,k,iano2))/10.E-6_wp)
        
                avdet = max(0._wp,bgctra(j,k,idet))
                remsulf = o2lim*nlim*sulfate_reduction*dtb 
                detnew = avdet/(1._wp+remsulf)
                remin = avdet - detnew

                ! if we assume full oxdidation of ammo with sulfur according to Thamdrup (anammox tpye)
                ! we get N2 instead of ammonium and a smaller alkalinity change
                ! LR: Check use of rnit/ralk !!
                bgctra(j,k,ialkali) = bgctra(j,k,ialkali)  + 2._wp*rnit * remin + remin
                bgctra(j,k,igasnit)   = bgctra(j,k,igasnit)  + 0.5_wp *rnit  *remin

                bgctra(j,k,iphosph) = bgctra(j,k,iphosph) + remin
                bgctra(j,k,iiron)   = bgctra(j,k,iiron) + riron * remin
                bgctra(j,k,ih2s)    = bgctra(j,k,ih2s) + 2._wp*rnit * remin + remin!

                bgctend(j,k,kn2b) = bgctend(j,k,kn2b)  + 3._wp*rnit *remin*(pddpo(j,k) + surface_height)
                !sulphate reduction indirectly effects O2 bugdet, which is budgeted in h2obudget
                bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob) - (ro2ammo + 0.5_wp*rnit)*remin*(pddpo(j,k) + surface_height)
       
             ELSE
                o2lim = 1._wp - bgctra(j,k,ioxygen)/thresh_o2

                xa = max(0._wp,bgctra(j,k,idet))
                xn = xa / (1._wp + sulfate_reduction*o2lim)
                remin = MAX(0._wp, xa-xn)

 
                bgctra(j,k,idet)    = bgctra(j,k,idet)    -        remin
                bgctra(j,k,ialkali) = bgctra(j,k,ialkali) + ralk * remin 
                bgctra(j,k,ih2s)    = bgctra(j,k,ih2s) + ralk * remin 
                bgctra(j,k,isco212) = bgctra(j,k,isco212) + rcar * remin
                bgctra(j,k,iphosph) = bgctra(j,k,iphosph) +        remin

                bgctra(j,k,iano3)   = bgctra(j,k,iano3)   + rnit  * remin
                bgctra(j,k,iiron)   = bgctra(j,k,iiron)   + riron * remin
     
                bgctend(j,k,kn2b) = bgctend(j,k,kn2b) + 2._wp * ralk * remin * (pddpo(j,k) + surface_height) 
                bgctend(j,k,kh2ob) = bgctend(j,k,kh2ob) - ro2ut * remin * (pddpo(j,k) + surface_height)
             ENDIF
             !!!! N-cycle !!!!!!!!

             bgctend(j,k,ksred) = remin / dtbgc 
             bgctend(j,k,kh2sprod) =  ralk * remin /dtbgc
             bgctend(j,k,kh2sloss) =  0._wp

       else
             ! HS oxidation 
               o2lim = bgctra(j,k,ioxygen)/(bgctra(j,k,ioxygen) + bkh2sox)
               xa = max(0._wp,bgctra(j,k,ih2s))
               xn = xa / ( 1._wp + rh2sox*o2lim)
               oxid = max(0._wp, xa-xn)

               bgctra(j,k,ih2s) = bgctra(j,k,ih2s) - oxid

               bgctra(j,k,ialkali) = bgctra(j,k,ialkali) - 2._wp  * oxid
               bgctend(j,k,kn2b) = bgctend(j,k,kn2b) - 2._wp * oxid * (pddpo(j,k) + surface_height) 
   
               bgctend(j,k,kh2sprod) =  0._wp
               bgctend(j,k,kh2sloss) =  oxid /dtbgc


       ENDIF ! O2 < thresh_sred



       bgctend(j,k,kaou)   = satoxy(j,k) - bgctra(j,k,ioxygen)
  
      ENDIF ! wet cells
     ENDDO ! k=1,kpke
  ENDDO ! j=start_idx,end_idx 

!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

END SUBROUTINE ocprod
END MODULE
