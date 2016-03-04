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

      USE mo_sedmnt,               ONLY: ks, powtra, sedlay, burial, sedhpl

      USE mo_exception, ONLY      : message, finish

      USE mo_model_domain,   ONLY: t_patch_3D, t_patch

      USE mo_ocean_nml,  ONLY: no_tracer, n_zlev

      USE mo_kind,     ONLY: wp

      USE mo_impl_constants,      ONLY: max_char_length


      USE mo_grid_config,         ONLY: n_dom


      USE mo_hamocc_types,       ONLY: t_hamocc_diag, t_hamocc_state, &
    &                                  t_hamocc_sed, t_hamocc_tend,   &
    &                                  t_hamocc_acc, t_hamocc_monitor                    
   





       
      USE mo_parallel_config,     ONLY: nproma

      USE mo_hamocc_nml,         ONLY: io_stdo_bgc


      IMPLICIT NONE

      PUBLIC

      !TYPE(t_hamocc_state), ALLOCATABLE, TARGET     :: hamocc_state(:)
      TYPE(t_hamocc_state), TARGET                  :: hamocc_state

      CONTAINS

!================================================================================== 

  function get_level_index_by_depth(depth,tiestw) result(level_index)

     real(wp), intent(in) :: depth,tiestw(n_zlev+1)
     
     integer :: level_index

     level_index = 1
     do while(tiestw(level_index+1) < depth .and. level_index < n_zlev)
        level_index = level_index + 1
     end do

  end function get_level_index_by_depth

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

      USE mo_carbch, ONLY: bgctra
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

      USE MO_CARBCH, ONLY: bgctra, co3, hi, bgctend, bgcflux, &
 &                         akw3,ak13,ak23,akb3,aksp 
      USE MO_PARAM1_BGC, ONLY: n_bgctra, issso12,         &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipowasi, ipowafe, kn2b,    &
 &                             kh2ob, korginp, ksilinp,   &
&                              kcalinp,keuexp

      USE MO_BIOMOD,  ONLY:kbo,bolay

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

!      CALL message(TRIM(routine), 'start' )

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,kpke,jk,itrac) HAMOCC_OMP_DEFAULT_SCHEDULE
      DO jc=start_index,end_index 
        kpke=klevs(jc)
        IF (pddpo(jc, 1) .GT. 0.5_wp) THEN

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
   SUBROUTINE set_bgc_output_pointers(timelevel,p_diag,p_prog)
    INTEGER, INTENT(in) :: timelevel
    TYPE(t_hamocc_diag) :: p_diag
    TYPE(t_hydro_ocean_prog) :: p_prog
    
    CHARACTER(LEN=max_char_length) :: timelevel_str
    !-------------------------------------------------------------------------
    WRITE(timelevel_str,'(a,i2.2)') '_TL',timelevel
    !write(0,*)'>>>>>>>>>>>>>>>> T timelevel_str:',TRIM(timelevel_str)
    p_diag%dic(:,:,:)        =  p_prog%tracer(:,:,:,3)
    p_diag%alk(:,:,:)        =  p_prog%tracer(:,:,:,4)
    p_diag%po4(:,:,:)        =  p_prog%tracer(:,:,:,5)
    p_diag%no3(:,:,:)        =  p_prog%tracer(:,:,:,6)
    p_diag%n2(:,:,:)         =  p_prog%tracer(:,:,:,7)
    p_diag%phy(:,:,:)        =  p_prog%tracer(:,:,:,8)
    p_diag%zoo(:,:,:)        =  p_prog%tracer(:,:,:,9)
    p_diag%cya(:,:,:)        =  p_prog%tracer(:,:,:,10)
    p_diag%o2(:,:,:)         =  p_prog%tracer(:,:,:,11)
    p_diag%si(:,:,:)         =  p_prog%tracer(:,:,:,12)
    p_diag%doc(:,:,:)        =  p_prog%tracer(:,:,:,13)
    p_diag%n2o(:,:,:)        =  p_prog%tracer(:,:,:,14)
    p_diag%det(:,:,:)        =  p_prog%tracer(:,:,:,15)
    p_diag%doccya(:,:,:)     =  p_prog%tracer(:,:,:,16)
    p_diag%iron(:,:,:)       =  p_prog%tracer(:,:,:,17)
    p_diag%dms(:,:,:)        =  p_prog%tracer(:,:,:,18)
    p_diag%calc(:,:,:)       =  p_prog%tracer(:,:,:,19)
    p_diag%opal(:,:,:)       =  p_prog%tracer(:,:,:,20)
    p_diag%dust(:,:,:)       =  p_prog%tracer(:,:,:,21)
   
  END SUBROUTINE 
!================================================================================== 
  SUBROUTINE set_bgc_tendencies_output(start_idx, end_idx, &
&             klevs,pddpo,jb,p_tend, p_diag, p_sed)

      
      USE mo_biomod, ONLY: bolay
      USE mo_carbch, ONLY: bgctend, bgcflux, hi, co3, bgctra, sedfluxo, &
 &                         akw3, akb3, aksp, ak13, ak23

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
&                              kbacfra, kdelsil, kdelcar, kbacfrac, &
&                              kdmsflux, kdmsprod, kdmsbac, kdmsuv, &
&                              keuexp, kplim, kflim, knlim,kcalex90,&
&                              kopex90, kgraton, kexudp, kexudz, &
&                              kzdy, kpdy,kcoex1000,kcoex2000, &
&                              kopex1000,kopex2000,kcalex1000,&
&                              kcalex2000
  
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
             p_tend%bacfrac(jc,jk,jb) = bgctend(jc,jk,kbacfrac)
             p_tend%remina(jc,jk,jb) = bgctend(jc,jk,kremin)
             p_tend%remins(jc,jk,jb) = bgctend(jc,jk,ksred)
             p_tend%reminn(jc,jk,jb) = bgctend(jc,jk,kdenit)
             p_tend%h2obudget(jc,jk,jb) = bgctend(jc,jk,kh2ob)
             p_tend%n2budget(jc,jk,jb) = bgctend(jc,jk,kn2b)
             p_tend%delsil(jc,jk,jb) = bgctend(jc,jk,kdelsil)
             p_tend%delcar(jc,jk,jb) = bgctend(jc,jk,kdelcar)
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
        ENDDO
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

      USE mo_carbch, ONLY: bgctra, hi, co3
      USE mo_biomod, ONLY: kbo,bolay
      USE mo_param1_bgc, ONLY: n_bgctra, issso12, &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipowasi, ipowafe
  


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
        p_sed%bc12(jc,jb) = burial(jc,issso12)
        p_sed%bsil(jc,jb) = burial(jc,issssil)
        p_sed%bter(jc,jb) = burial(jc,issster) 
        p_sed%bolay(jc,jb) = bolay(jc)
        p_sed%kbo(jc,jb) = kbo(jc)
        DO jk =1,ks
             ! Solid sediment
             p_sed%so12(jc,jk,jb) = sedlay(jc,jk,issso12)
             p_sed%sc12(jc,jk,jb) = sedlay(jc,jk,issso12)
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
             p_sed%sedhi(jc,jk,jb) = sedhpl(jc,jk)
        ENDDO
      ENDIF
 
      ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
      END SUBROUTINE

!================================================================================== 


  SUBROUTINE print_bgc_parameters
  USE mo_biomod, ONLY      : phytomi, grami, remido, dyphy, zinges,        &
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
       &                     remido_cya, dremdoc_cya, buoyancyspeed_cya, &
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
   cpara_name='phytomi'
   write(cpara_val,'(ES9.2)')phytomi
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='dyphy'
   write(cpara_val,'(f6.4)')dyphy
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='bkphy'
   write(cpara_val,'(ES9.2)')bkphy
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   ! Zooplankton
   cpara_name='ZOOPLANKTON'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='grami'
   write(cpara_val,'(ES9.2)')grami
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='zinges'
   write(cpara_val,'(f5.2)')zinges
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='grazra'
   write(cpara_val,'(f6.4)')grazra
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='bkzoo'
   write(cpara_val,'(ES9.2)')bkzoo
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   ! Cyanobacteria
   cpara_name='CYANOBACTERIA'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='l_cyadyn'
   write(cpara_val,'(L1)')l_cyadyn
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='n2_fixation'
   write(cpara_val,'(f6.4)')n2_fixation
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='buoyancyspeed_cya'
   write(cpara_val,'(f6.4)')buoyancyspeed_cya
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='cycdec'
   write(cpara_val,'(f6.4)')cycdec
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='pi_alpha_cya'
   write(cpara_val,'(f6.4)')pi_alpha_cya
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='Topt_cya'
   write(cpara_val,'(f6.2)')Topt_cya
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='T1_cya'
   write(cpara_val,'(f6.2)')T1_cya
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='T2_cya'
   write(cpara_val,'(f6.2)')T2_cya
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='bkcya_p'
   write(cpara_val,'(ES9.2)')bkcya_p
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='bkcya_n'
   write(cpara_val,'(ES9.2)')bkcya_n
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='bkcya_fe'
   write(cpara_val,'(ES9.2)')bkcya_fe
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='remido_cya'
   write(cpara_val,'(ES9.2)')remido_cya
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='doccya_cya'
   write(cpara_val,'(f6.4)')doccya_fac
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  
   ! Detritus
   cpara_name='DETRITUS'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='i_settling'
   write(cpara_val,'(i1)')i_settling
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='drempoc'
   write(cpara_val,'(f6.4)')drempoc
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='denitrification'
   write(cpara_val,'(f6.4)')denitrification
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='sulfate reduction'
   write(cpara_val,'(f6.4)')sulfate_reduction
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='thresh_aerob'
   write(cpara_val,'(ES9.2)')thresh_aerob
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='thresh_sred'
   write(cpara_val,'(ES9.2)')thresh_sred
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )


   ! DOC
   cpara_name='DOC'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='remido'
   write(cpara_val,'(f6.4)')remido
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  

   ! Opal
   cpara_name='Si/Opal'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='bkopal'
   write(cpara_val,'(ES9.2)')bkopal
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='dremopal'
   write(cpara_val,'(f6.4)')dremopal
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='ropal'
   write(cpara_val,'(f6.2)')ropal
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='sinkspeed_opal'
   write(cpara_val,'(f6.2)')sinkspeed_opal
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   ! Iron
   cpara_name='Iron'
   cpara_val="========"
   cpara_name='perc_diron'
   write(cpara_val,'(f6.4)')perc_diron
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='fesoly'
   write(cpara_val,'(ES9.2)')fesoly
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='riron'
   write(cpara_val,'(f6.4)')riron
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='relaxfe'
   write(cpara_val,'(ES9.2)')relaxfe
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   ! Sediment
   cpara_name='Sediment'
   cpara_val="========"
   cpara_name='denit_sed'
   write(cpara_val,'(f6.4)')denit_sed
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='disso_po'
   write(cpara_val,'(f6.4)')disso_po
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='disso_op'
   write(cpara_val,'(ES9.2)')disso_op
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='disso_cal'
   write(cpara_val,'(ES9.2)')disso_cal
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='sred_sed'
   write(cpara_val,'(ES9.2)')sred_sed
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   

   ! Calc
   cpara_name='Calc'
   cpara_val="========"
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='dremcalc'
   write(cpara_val,'(f6.4)')dremcalc
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   cpara_name='sinkspeed_calc'
   write(cpara_val,'(f6.2)')sinkspeed_calc
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   cpara_name='======================='
   cpara_val="==========="
   CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  END SUBROUTINE print_bgc_parameters

  SUBROUTINE print_wpoc
   USE mo_carbch, ONLY: wpoc    
    CHARACTER(LEN=max_char_length) :: &
                cpara_name,cpara_val

    INTEGER:: k

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

 END MODULE
