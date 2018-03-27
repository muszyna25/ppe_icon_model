#ifndef __NO_ICON_OCEAN__
#include "hamocc_omp_definitions.inc"      
   MODULE mo_bgc_icon_comm

! icon specific routines for output, update etc

      USE mo_ocean_diagnostics_types,   ONLY:  t_ocean_regions

      USE mo_ocean_types,          ONLY:  t_hydro_ocean_prog

      USE mo_control_bgc,          ONLY: dtb,bgc_gin, bgc_arctic, bgc_lab, & 
       &                                 bgc_natl, bgc_atl, bgc_tatl, &
       &                                 bgc_tropac, &
       &                                 bgc_land, bgc_ind, &
       &                                 bgc_soce, bgc_npac, bgc_carb

      USE mo_sedmnt,               ONLY: powtra, sedlay, burial, sedhpl

      USE mo_exception, ONLY      : message, finish

      USE mo_model_domain,   ONLY: t_patch_3D, t_patch

      USE mo_ocean_nml,  ONLY: no_tracer, n_zlev

      USE mo_kind,     ONLY: wp

      USE mo_impl_constants,      ONLY: max_char_length


      USE mo_grid_config,         ONLY: n_dom


      USE mo_hamocc_types,       ONLY: t_hamocc_diag, t_hamocc_state, &
    &                                  t_hamocc_sed, t_hamocc_tend,   &
    &                                  t_hamocc_monitor                    
   

       
      USE mo_parallel_config,     ONLY: nproma

      USE mo_hamocc_nml,         ONLY: io_stdo_bgc,ks


      IMPLICIT NONE

      PUBLIC

      TYPE(t_hamocc_state), TARGET                  :: hamocc_state


      INTERFACE to_bgcout
        module procedure to_bgcout_real
        module procedure to_bgcout_int
        module procedure to_bgcout_logical
      END INTERFACE



      CONTAINS

!================================================================================== 
       SUBROUTINE ini_bgc_regions

        IMPLICIT NONE
        TYPE(t_ocean_regions)         :: ocean_regions


         bgc_land=ocean_regions%land                                !0
         bgc_gin= ocean_regions%greenland_iceland_norwegian_sea     !1
         bgc_arctic=ocean_regions%arctic_ocean                      !2
         bgc_lab= ocean_regions%labrador_sea                        !3
         bgc_natl= ocean_regions%north_atlantic                     !4
         bgc_tatl= ocean_regions%tropical_atlantic                  !5 
         bgc_soce=ocean_regions%southern_ocean                      !6
         bgc_ind=ocean_regions%indian_ocean                         !7
         bgc_tropac= ocean_regions%tropical_pacific                 !8
         bgc_npac=ocean_regions%north_pacific                       !9
         bgc_carb=ocean_regions%caribbean                           !-33
 
       END SUBROUTINE
!================================================================================== 
    
      SUBROUTINE update_icon(start_idx, end_idx, &
&             klevs, pddpo, ptracer)

      USE mo_memory_bgc, ONLY: bgctra
      USE mo_param1_bgc, ONLY: n_bgctra


      REAL(wp)     :: ptracer(nproma,n_zlev,no_tracer+n_bgctra)    
      INTEGER, INTENT(in)::klevs(nproma)
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]

      INTEGER :: jc, jk, kpke
      INTEGER :: start_idx, end_idx
      INTEGER :: itrac
      CHARACTER(LEN=max_char_length), PARAMETER :: &
                routine = 'update_icon'

     ! CALL message(TRIM(routine), 'start' )
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,kpke,jk,itrac) HAMOCC_OMP_DEFAULT_SCHEDULE
      DO jc=start_idx,end_idx 
        kpke=klevs(jc)
        IF (pddpo(jc, 1) .GT. 0.5_wp) THEN
        DO jk =1,kpke
          DO itrac=no_tracer+1,no_tracer+n_bgctra
             ptracer(jc,jk,itrac) = bgctra(jc,jk,itrac-no_tracer)
          ENDDO 
        ENDDO
        ENDIF
      ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
      END SUBROUTINE

!================================================================================== 
      SUBROUTINE update_bgc(start_index, end_index, &
&             klevs,pddpo,jb,ptracer,p_diag,p_sed,p_tend)

      USE mo_memory_bgc, ONLY: bgctra, co3, hi, bgctend, bgcflux, &
 &                         akw3,ak13,ak23,akb3,aksp,satoxy, &
 &                         satn2, satn2o, solco2,kbo,bolay
      USE MO_PARAM1_BGC, ONLY: n_bgctra, issso12,         &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipowasi, ipowafe, kn2b,    &
 &                             kh2ob, korginp, ksilinp,   &
&                              kcalinp,keuexp, ipowh2s


      USE mo_sedmnt,  ONLY: pown2bud, powh2obud

      REAL(wp)     :: ptracer(nproma,n_zlev,no_tracer+n_bgctra)    
      INTEGER, INTENT(in)::klevs(nproma), jb
      TYPE(t_hamocc_diag) :: p_diag
      TYPE(t_hamocc_sed) :: p_sed
      TYPE(t_hamocc_tend) :: p_tend
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]

      INTEGER :: jc, jk, kpke
      INTEGER :: start_index, end_index
      INTEGER :: itrac
      CHARACTER(LEN=max_char_length), PARAMETER :: &
                routine = 'update_icon'


!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,kpke,jk,itrac) HAMOCC_OMP_DEFAULT_SCHEDULE
      DO jc=start_index,end_index 
        kpke=klevs(jc)
        IF (pddpo(jc, 1) .GT. 0.5_wp) THEN

          satn2(jc)  = p_tend%satn2(jc,jb)    
          satn2o(jc) = p_tend%satn2o(jc,jb)    
          solco2(jc) = p_tend%solco2(jc,jb)    
        DO jk =1,kpke
          DO itrac=no_tracer+1,no_tracer+n_bgctra
             bgctra(jc,jk,itrac-no_tracer)=ptracer(jc,jk,itrac) 
          ENDDO 
          hi(jc, jk) = p_diag%hi(jc,jk,jb)    
          co3(jc, jk) = p_diag%co3(jc,jk,jb)    
          aksp(jc, jk) = p_tend%aksp(jc,jk,jb)    
          ak13(jc, jk) = p_tend%ak1(jc,jk,jb)    
          ak23(jc, jk) = p_tend%ak2(jc,jk,jb)    
          akb3(jc, jk) = p_tend%akb(jc,jk,jb)    
          akw3(jc, jk) = p_tend%akw(jc,jk,jb)    
          satoxy(jc, jk) = p_tend%satoxy(jc,jk,jb)    
          bgctend(jc,jk,kh2ob) =  p_tend%h2obudget(jc,jk,jb) 
          bgctend(jc,jk,kn2b) =  p_tend%n2budget(jc,jk,jb) 
        ENDDO
     ! Sediment
         kbo(jc) = p_sed%kbo(jc,jb)
         bolay(jc) = p_sed%bolay(jc,jb)
        ! Burial layers
         burial(jc,issso12) = p_sed%bo12(jc,jb) 
         burial(jc,isssc12) = p_sed%bc12(jc,jb)
         burial(jc,issssil) = p_sed%bsil(jc,jb)
         burial(jc,issster) = p_sed%bter(jc,jb)
        DO jk =1,ks
             ! Solid sediment
             sedlay(jc,jk,issso12) = p_sed%so12(jc,jk,jb) 
             sedlay(jc,jk,isssc12) = p_sed%sc12(jc,jk,jb) 
             sedlay(jc,jk,issssil) = p_sed%ssil(jc,jk,jb) 
             sedlay(jc,jk,issster) = p_sed%ster(jc,jk,jb)  
             ! Pore water
             powtra(jc,jk,ipowaic) = p_sed%pwic(jc,jk,jb)  
             powtra(jc,jk,ipowaal) = p_sed%pwal(jc,jk,jb) 
             powtra(jc,jk,ipowaph) = p_sed%pwph(jc,jk,jb) 
             powtra(jc,jk,ipowaox) = p_sed%pwox(jc,jk,jb)  
             powtra(jc,jk,ipowasi) = p_sed%pwsi(jc,jk,jb) 
             powtra(jc,jk,ipowafe) = p_sed%pwfe(jc,jk,jb) 
             powtra(jc,jk,ipown2)  = p_sed%pwn2(jc,jk,jb) 
             powtra(jc,jk,ipowno3) = p_sed%pwno3(jc,jk,jb) 
             powtra(jc,jk,ipowh2s) = p_sed%pwh2s(jc,jk,jb) 
             sedhpl(jc,jk)         = p_sed%sedhi(jc,jk,jb) 
             powh2obud(jc,jk)    = p_sed%pwh2ob(jc,jk,jb) 
             pown2bud(jc,jk)     = p_sed%pwn2b(jc,jk,jb) 
        ENDDO
       ENDIF

      ENDDO
!HAMICC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
      END SUBROUTINE

!================================================================================== 
  SUBROUTINE set_bgc_tendencies_output(start_idx, end_idx, &
&             klevs,pddpo,jb,p_tend, p_diag, p_sed)

      
      USE mo_memory_bgc, ONLY: bgctend, bgcflux, hi, co3, bgctra, sedfluxo, &
 &                         akw3, akb3, aksp, ak13, ak23, satoxy, satn2, &
 &                         satn2o, solco2, bolay

      USE mo_param1_bgc, ONLY: kphosy, ksred, kremin, kdenit, &
 &                             kcflux, koflux, knflux, knfixd, &
 &                             knfix, kgraz, ksilpro, kprorca, &
 &                             kn2oflux,korginp, ksilinp, kcalinp,  &
 &                             kprcaca, kcoex90, issso12, &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipowasi, ipowafe, kpho_cya, &
 &                             kcyaloss, kn2b, kh2ob, kprodus, &
&                              kbacfra, kdelsil, kdelcar,  &
&                              kdmsflux, kdmsprod, kdmsbac, kdmsuv, &
&                              keuexp, kplim, kflim, knlim,kcalex90,&
&                              kopex90, kgraton, kexudp, kexudz, &
&                              kzdy, kpdy,kcoex1000,kcoex2000, &
&                              kopex1000,kopex2000,kcalex1000,&
&                              kcalex2000, kaou, kcTlim, kcLlim, &
&                              kcPlim, kcFlim, ipowh2s,kh2sprod, &   
&                              kh2sloss
  
      USE mo_sedmnt, ONLY : pown2bud, powh2obud, sedtend, &
&                           isremino, isreminn, isremins

      TYPE(t_hamocc_tend) :: p_tend
      TYPE(t_hamocc_diag) :: p_diag
      TYPE(t_hamocc_sed) :: p_sed


      INTEGER, INTENT(in) :: klevs(nproma)
      INTEGER, INTENT(in) :: start_idx, end_idx,jb
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]


      INTEGER :: jc, jk, kpke
      INTEGER :: itrac

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,jk,kpke) HAMOCC_OMP_DEFAULT_SCHEDULE
      DO jc=start_idx,end_idx 
        IF (pddpo(jc, 1) .GT. 0.5_wp) THEN
        p_tend%cflux(jc,jb) = bgcflux(jc,kcflux)
        p_tend%oflux(jc,jb) = bgcflux(jc,koflux)
        p_tend%nflux(jc,jb) = bgcflux(jc,knflux)
        p_tend%dmsflux(jc,jb) = bgcflux(jc,kdmsflux)
        p_tend%orginp(jc,jb) = bgcflux(jc,korginp)
        p_tend%silinp(jc,jb) = bgcflux(jc,ksilinp)
        p_tend%calinp(jc,jb) = bgcflux(jc,kcalinp)
        p_tend%n2oflux(jc,jb) = bgcflux(jc,kn2oflux)
        p_tend%nfixd(jc,jb) = bgcflux(jc,knfixd)
        p_tend%prcaca(jc,jb) = bgcflux(jc,kprcaca)
        p_tend%prorca(jc,jb) = bgcflux(jc,kprorca)
        p_tend%silpro(jc,jb) = bgcflux(jc,ksilpro)
        p_tend%produs(jc,jb) = bgcflux(jc,kprodus)
        p_tend%coex90(jc,jb) = bgcflux(jc,kcoex90)
        p_tend%calex90(jc,jb) = bgcflux(jc,kcalex90)
        p_tend%opex90(jc,jb) = bgcflux(jc,kopex90)
        p_tend%coex1000(jc,jb) = bgcflux(jc,kcoex1000)
        p_tend%calex1000(jc,jb) = bgcflux(jc,kcalex1000)
        p_tend%opex1000(jc,jb) = bgcflux(jc,kopex1000)
        p_tend%coex2000(jc,jb) = bgcflux(jc,kcoex2000)
        p_tend%calex2000(jc,jb) = bgcflux(jc,kcalex2000)
        p_tend%opex2000(jc,jb) = bgcflux(jc,kopex2000)
        kpke=klevs(jc)
        DO jk =1,kpke
             p_tend%npp(jc,jk,jb) = bgctend(jc,jk,kphosy)
             p_tend%graz(jc,jk,jb) = bgctend(jc,jk,kgraz)
             p_tend%zoomor(jc,jk,jb) = bgctend(jc,jk,kzdy)
             p_tend%phymor(jc,jk,jb) = bgctend(jc,jk,kpdy)
             p_tend%exudz(jc,jk,jb) = bgctend(jc,jk,kexudz)
             p_tend%graton(jc,jk,jb) = bgctend(jc,jk,kgraton)
             p_tend%exud(jc,jk,jb) = bgctend(jc,jk,kexudp)
             p_tend%nfix(jc,jk,jb) = bgctend(jc,jk,knfix)
             p_tend%phoc(jc,jk,jb) = bgctend(jc,jk,kpho_cya)
             p_tend%cyloss(jc,jk,jb) = bgctend(jc,jk,kcyaloss)
             p_tend%bacfra(jc,jk,jb) = bgctend(jc,jk,kbacfra)
             p_tend%remina(jc,jk,jb) = bgctend(jc,jk,kremin)
             p_tend%remins(jc,jk,jb) = bgctend(jc,jk,ksred)
             p_tend%reminn(jc,jk,jb) = bgctend(jc,jk,kdenit)
             p_tend%h2obudget(jc,jk,jb) = bgctend(jc,jk,kh2ob)
             p_tend%n2budget(jc,jk,jb) = bgctend(jc,jk,kn2b)
             p_tend%delsil(jc,jk,jb) = bgctend(jc,jk,kdelsil)
             p_tend%delcar(jc,jk,jb) = bgctend(jc,jk,kdelcar)
             p_tend%h2sprod(jc,jk,jb) = bgctend(jc,jk,kh2sprod)
             p_tend%h2sloss(jc,jk,jb) = bgctend(jc,jk,kh2sloss)
             p_tend%dmsprod(jc,jk,jb) = bgctend(jc,jk,kdmsprod)
             p_tend%dmsbac(jc,jk,jb) = bgctend(jc,jk,kdmsbac)
             p_tend%dmsuv(jc,jk,jb) = bgctend(jc,jk,kdmsuv)
             p_tend%euexp(jc,jk,jb) = bgctend(jc,jk,keuexp)
             p_diag%hi(jc,jk,jb)     = hi(jc, jk)
             p_diag%co3(jc,jk,jb)    = co3(jc,jk) 
             p_tend%akb(jc,jk,jb)    = akb3(jc,jk) 
             p_tend%akw(jc,jk,jb)    = akw3(jc,jk) 
             p_tend%ak1(jc,jk,jb)    = ak13(jc,jk) 
             p_tend%ak2(jc,jk,jb)    = ak23(jc,jk) 
             p_tend%aksp(jc,jk,jb)   = aksp(jc,jk) 
             p_tend%flim(jc,jk,jb) = bgctend(jc,jk,kflim)
             p_tend%nlim(jc,jk,jb) = bgctend(jc,jk,knlim)
             p_tend%plim(jc,jk,jb) = bgctend(jc,jk,kplim)
             p_tend%cTlim(jc,jk,jb) = bgctend(jc,jk,kcTlim)
             p_tend%cLlim(jc,jk,jb) = bgctend(jc,jk,kcLlim)
             p_tend%cPlim(jc,jk,jb) = bgctend(jc,jk,kcPlim)
             p_tend%cFlim(jc,jk,jb) = bgctend(jc,jk,kcFlim)
             p_tend%satoxy(jc,jk,jb)    = satoxy(jc,jk) 
             p_tend%aou(jc,jk,jb) = bgctend(jc,jk,kaou)
        ENDDO
        p_tend%satn2(jc,jb)    = satn2(jc) 
        p_tend%satn2o(jc,jb)    = satn2o(jc) 
        p_tend%solco2(jc,jb)    = solco2(jc) 
        ! Sediment
        ! Burial layers
        p_sed%bo12(jc,jb) = burial(jc,issso12)
        p_sed%bc12(jc,jb) = burial(jc,isssc12)
        p_sed%bsil(jc,jb) = burial(jc,issssil)
        p_sed%bter(jc,jb) = burial(jc,issster) 
        ! Sediment-ocean fluxes
        p_tend%sedflic(jc,jb) = sedfluxo(jc,ipowaic) 
        p_tend%sedflal(jc,jb) = sedfluxo(jc,ipowaal) 
        p_tend%sedflph(jc,jb) = sedfluxo(jc,ipowaph) 
        p_tend%sedflox(jc,jb) = sedfluxo(jc,ipowaox) 
        p_tend%sedflsi(jc,jb) = sedfluxo(jc,ipowasi) 
        p_tend%sedflfe(jc,jb) = sedfluxo(jc,ipowafe) 
        p_tend%sedfln2(jc,jb) = sedfluxo(jc,ipown2) 
        p_tend%sedflno3(jc,jb) = sedfluxo(jc,ipowno3) 
        p_tend%sedflh2s(jc,jb) = sedfluxo(jc,ipowh2s) 
        DO jk =1,ks
             ! Solid sediment
             p_sed%so12(jc,jk,jb) = sedlay(jc,jk,issso12)
             p_sed%sc12(jc,jk,jb) = sedlay(jc,jk,isssc12)
             p_sed%ssil(jc,jk,jb) = sedlay(jc,jk,issssil)
             p_sed%ster(jc,jk,jb) = sedlay(jc,jk,issster)
             ! Pore water
             p_sed%pwic(jc,jk,jb) = powtra(jc,jk,ipowaic)
             p_sed%pwal(jc,jk,jb) = powtra(jc,jk,ipowaal)
             p_sed%pwph(jc,jk,jb) = powtra(jc,jk,ipowaph)
             p_sed%pwox(jc,jk,jb) = powtra(jc,jk,ipowaox)
             p_sed%pwsi(jc,jk,jb) = powtra(jc,jk,ipowasi)
             p_sed%pwfe(jc,jk,jb) = powtra(jc,jk,ipowafe)
             p_sed%pwn2(jc,jk,jb) = powtra(jc,jk,ipown2)
             p_sed%pwno3(jc,jk,jb) = powtra(jc,jk,ipowno3)
             p_sed%pwh2s(jc,jk,jb) = powtra(jc,jk,ipowh2s)
             p_sed%sedhi(jc,jk,jb) = sedhpl(jc,jk)
             p_sed%pwh2ob(jc,jk,jb) = powh2obud(jc,jk)
             p_sed%pwn2b(jc,jk,jb) = pown2bud(jc,jk)

             ! tendencies
             p_tend%sedro2(jc,ks,jb) = sedtend(jc,jk,isremino) 
             p_tend%sedrn(jc,ks,jb)  = sedtend(jc,jk,isreminn) 
             p_tend%sedrs(jc,ks,jb)  = sedtend(jc,jk,isremins) 
             
        ENDDO
      ENDIF
      ENDDO
!HAMOC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  END SUBROUTINE 

!================================================================================== 
    
      SUBROUTINE initial_update_icon(start_index, end_index, &
&             klevs,pddpo,jb,ptracer, p_sed,p_diag)

      USE mo_memory_bgc, ONLY: bgctra, hi, co3, kbo,bolay
      USE mo_param1_bgc, ONLY: n_bgctra, issso12, &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipowasi, ipowafe, ipowh2s
  


      REAL(wp)     :: ptracer(nproma,n_zlev,no_tracer+n_bgctra)    
      INTEGER, INTENT(in)::klevs(nproma)
      TYPE(t_hamocc_sed) :: p_sed
      TYPE(t_hamocc_diag) :: p_diag
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]

      INTEGER :: jc, jk, kpke,jb
      INTEGER :: start_index, end_index
      INTEGER :: itrac
      CHARACTER(LEN=max_char_length), PARAMETER :: &
                routine = 'update_icon'

     ! CALL message(TRIM(routine), 'start' )
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,kpke,jk) HAMOCC_OMP_DEFAULT_SCHEDULE
      DO jc=start_index,end_index 
        kpke=klevs(jc)
        IF (pddpo(jc, 1) .GT. 0.5_wp) THEN
        DO jk =1,kpke
          DO itrac=no_tracer+1,no_tracer+n_bgctra
             ptracer(jc,jk,itrac) = bgctra(jc,jk,itrac-no_tracer)
          ENDDO
             p_diag%hi(jc,jk,jb)  = hi(jc, jk)
             p_diag%co3(jc,jk,jb) = co3(jc,jk) 
        ENDDO
        ! Sediment
        ! Burial layers
        p_sed%bo12(jc,jb) = burial(jc,issso12)
        p_sed%bc12(jc,jb) = burial(jc,isssc12)
        p_sed%bsil(jc,jb) = burial(jc,issssil)
        p_sed%bter(jc,jb) = burial(jc,issster) 
        p_sed%bolay(jc,jb) = bolay(jc)
        p_sed%kbo(jc,jb) = kbo(jc)
        DO jk =1,ks
             ! Solid sediment
             p_sed%so12(jc,jk,jb) = sedlay(jc,jk,issso12)
             p_sed%sc12(jc,jk,jb) = sedlay(jc,jk,isssc12)
             p_sed%ssil(jc,jk,jb) = sedlay(jc,jk,issssil)
             p_sed%ster(jc,jk,jb) = sedlay(jc,jk,issster)
             ! Pore water
             p_sed%pwic(jc,jk,jb) = powtra(jc,jk,ipowaic)
             p_sed%pwal(jc,jk,jb) = powtra(jc,jk,ipowaal)
             p_sed%pwph(jc,jk,jb) = powtra(jc,jk,ipowaph)
             p_sed%pwox(jc,jk,jb) = powtra(jc,jk,ipowaox)
             p_sed%pwsi(jc,jk,jb) = powtra(jc,jk,ipowasi)
             p_sed%pwfe(jc,jk,jb) = powtra(jc,jk,ipowafe)
             p_sed%pwn2(jc,jk,jb) = powtra(jc,jk,ipown2)
             p_sed%pwno3(jc,jk,jb) = powtra(jc,jk,ipowno3)
             p_sed%pwh2s(jc,jk,jb) = powtra(jc,jk,ipowh2s)
             p_sed%sedhi(jc,jk,jb) = sedhpl(jc,jk)
        ENDDO
      ENDIF
 
      ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
      END SUBROUTINE

!================================================================================== 


  SUBROUTINE print_bgc_parameters
  USE mo_memory_bgc, ONLY      : phytomi, grami, remido, dyphy, zinges,        &
       &                     epsher, grazra, spemor, gammap, gammaz, ecan, &
       &                     pi_alpha, fpar, bkphy, bkzoo, bkopal,         &
       &                     drempoc,                                      &
       &                     dremopal, dremn2o, sulfate_reduction,         &
       &                     dremcalc, n2_fixation, ro2ut, rcar, rnit,     &
       &                     rnoi, nitdem, n2prod, rcalc, ropal, calmax,   &
       &                     perc_diron, riron, fesoly, relaxfe,     &
       &                     denitrification,             rn2,             &
       &                     cycdec, pi_alpha_cya,cya_growth_max,          &
       &                     Topt_cya,T1_cya,T2_cya,bkcya_N, bkcya_P, bkcya_fe, &
       &                     buoyancyspeed_cya, &
       &                     doccya_fac, thresh_aerob, thresh_sred 

   USE mo_hamocc_nml, ONLY: i_settling, l_cyadyn, denit_sed, disso_po, &
      &                 sinkspeed_opal, sinkspeed_calc


   USE mo_sedmnt, ONLY: disso_op, disso_cal,sred_sed


  CHARACTER(LEN=max_char_length) :: &
                cpara_name,cpara_val

   cpara_name='========PARAMETER SETUP'
   cpara_val="============="
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   ! Phytoplankton
   cpara_name='PHYTOPLANKTON'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("phytomi",phytomi)
   CALL to_bgcout("dyphy",dyphy)
   CALL to_bgcout("bkphy",bkphy)

   ! Zooplankton
   cpara_name='ZOOPLANKTON'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("grami",grami)
   CALL to_bgcout("zinges",zinges)
   CALL to_bgcout("grazra",grazra)
   CALL to_bgcout("bkzoo",bkzoo)

   ! Cyanobacteria
   cpara_name='CYANOBACTERIA'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("l_cyadyn",l_cyadyn)
   CALL to_bgcout("n2_fixation",n2_fixation)
   CALL to_bgcout("n2_fixation",n2_fixation)
   CALL to_bgcout("buoyancyspeed_cya",buoyancyspeed_cya)
   CALL to_bgcout("cycdec",cycdec)
   CALL to_bgcout("pi_alpha_cya",pi_alpha_cya)
   CALL to_bgcout("Topt_cya",Topt_cya)
   CALL to_bgcout("T1",T1_cya)
   CALL to_bgcout("T2",T2_cya)
   CALL to_bgcout("bkcya_p",bkcya_p)
   CALL to_bgcout("bkcya_n",bkcya_n)
   CALL to_bgcout("bkcya_fe",bkcya_fe)
   CALL to_bgcout("doccya_fac",doccya_fac)
  
   ! Detritus
   cpara_name='DETRITUS'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("i_settling",i_settling)
   CALL to_bgcout("drempoc",drempoc)
   CALL to_bgcout("denitrification",denitrification)
   CALL to_bgcout("sulfate_reduction",sulfate_reduction)
   CALL to_bgcout("thresh_aerob",thresh_aerob)
   CALL to_bgcout("thresh_sred",thresh_sred)


   ! DOC
   cpara_name='DOC'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("remido",remido)
  

   ! Opal
   cpara_name='Si/Opal'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("bkopal",bkopal)
   CALL to_bgcout("dremopal",dremopal)
   CALL to_bgcout("ropal",ropal)
   CALL to_bgcout("sinkspeed_opal",sinkspeed_opal)

   ! Iron
   cpara_name='Iron'
   cpara_val="========"
   CALL to_bgcout("perc_diron",perc_diron)
   CALL to_bgcout("fesoly",fesoly)
   CALL to_bgcout("riron",riron)
   CALL to_bgcout("relaxfe",relaxfe)

   ! Sediment
   cpara_name='Sediment'
   cpara_val="========"
   CALL to_bgcout("denit_sed",denit_sed)
   CALL to_bgcout("disso_po",disso_po)
   CALL to_bgcout("disso_op",disso_op)
   CALL to_bgcout("disso_cal",disso_cal)
   CALL to_bgcout("sred_sed",sred_sed)
   

   ! Calc
   cpara_name='Calc'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("dremcalc",dremcalc)
   CALL to_bgcout("sinkspeed_calc",sinkspeed_calc)

   cpara_name='======================='
   cpara_val="==========="
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  END SUBROUTINE print_bgc_parameters
!================================================================================== 

  SUBROUTINE print_wpoc
   USE mo_memory_bgc, ONLY: wpoc, wdust
    CHARACTER(LEN=max_char_length) :: &
                cpara_name,cpara_val

    INTEGER:: k

   CALL to_bgcout("wdust",wdust)
   
   cpara_name='========WPOC [m/d]'
   cpara_val="============="
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   
   DO k=1,n_zlev  
    write(cpara_name,'(i2)')k
    write(cpara_val,'(f6.2)')wpoc(k)/dtb
    CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   enddo

   cpara_name='======================='
   cpara_val="==========="
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  

  END SUBROUTINE print_wpoc

!================================================================================== 

SUBROUTINE to_bgcout_real(cname,val)
  REAL(wp),INTENT(in) ::val
  CHARACTER( LEN = * ):: cname
  CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

  cpara_name=cname
  IF(abs(val)<100._wp.and.abs(val)>0.1)then
    write(cpara_val, '(f9.2)') val
  ELSE
    write(cpara_val, '(ES22.15)') val
  ENDIF
  CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc)
END SUBROUTINE

SUBROUTINE to_bgcout_logical(cname,val)
  LOGICAL,INTENT(in) ::val
  CHARACTER( LEN = * ):: cname
  CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

  cpara_name=cname
  cpara_val=merge(".true. ",".false.",val)
  CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc)
END SUBROUTINE

SUBROUTINE to_bgcout_int(cname,val)
  INTEGER,INTENT(in) ::val
  CHARACTER( LEN = * ):: cname
  CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

  cpara_name=cname
  write(cpara_val, '(i0)') val
  CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc)
END SUBROUTINE


 END MODULE
#endif
