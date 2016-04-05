
!--------------------------------------
#include "hamocc_omp_definitions.inc"
!--------------------------------------
MODULE mo_hamocc_diagnostics

   USE mo_kind,     ONLY: wp
   USE mo_sync,     ONLY: global_sum_array
   USE mo_sedmnt,   ONLY: ks, porsol, porwat, seddw
   USE mo_exception, ONLY: message
   USE mo_impl_constants, ONLY: max_char_length
   USE mo_hamocc_types, ONLY: t_hamocc_state
   USE mo_ocean_types, ONLY: t_hydro_ocean_state
   USE mo_model_domain, ONLY: t_patch_3d, t_patch 
   USE mo_grid_subset, ONLY: t_subset_range, get_index_range
   USE mo_hamocc_nml, ONLY: io_stdo_bgc, l_cyadyn
   USE mo_dynamics_config,     ONLY: nold, nnew
   USE mo_ocean_nml,   ONLY: n_zlev,no_tracer
   USE mo_bgc_constants, ONLY:  s2year, n2tgn, c2gtc, kilo
   USE mo_biomod, ONLY: p2gtc
   USE mo_control_bgc, ONLY: dtbgc
   USE mo_carbch, ONLY: totalarea
   USE mo_bgc_icon_comm, ONLY: to_bgcout
   USE mo_param1_bgc, ONLY: isco212, ialkali, iphosph,iano3, igasnit, &
&                           iphy, izoo, icya, ioxygen, isilica, idoc, &
&                           ian2o, idet, idoccya, iiron, icalc, iopal,&
&                           idust, idms


IMPLICIT NONE

PRIVATE

PUBLIC:: get_inventories, get_monitoring


CONTAINS

SUBROUTINE get_monitoring(hamocc_state,ocean_state,p_patch_3d)

USE mo_biomod, ONLY: rcar, rn2, nitdem,doccya_fac
TYPE(t_hamocc_state) :: hamocc_state
TYPE(t_hydro_ocean_state) :: ocean_state
TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d

CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val
INTEGER:: i_time_stat
REAL(wp) :: glob_n2b, glob_pwn2b

i_time_stat=nold(1)


! First try:  do not accumulate monitoring itself, but put accumulated fields into it
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%npp(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%phosy(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%phoc(:,:,:), &
& i_time_stat, hamocc_state%p_tend%monitor%phosy_cya(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%graz(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%grazing(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%graton(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%graton(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%exud(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%exud(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%exudz(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%exudz(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%zoomor(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%zoomor(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%phymor(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%phymor(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%delsil(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%delsil(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%delcar(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%delcar(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%bacfra(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%bacfra(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%remina(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%remina(1))
if(l_cyadyn)then
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%nfix(:,:,:), i_time_stat, hamocc_state%p_tend%monitor%n2fix(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%cyloss(:,:,:), i_time_stat, &
& hamocc_state%p_tend%monitor%cyaldet(1))
else
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%nfixd(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%n2fix(1), 1, ocean_state)
endif
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%bacfrac(:,:,:), &
&i_time_stat, hamocc_state%p_tend%monitor%bacfrac(1))
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_acc%reminn(:,:,:), &
& i_time_stat, hamocc_state%p_tend%monitor%wcdenit(1))
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%cflux(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%net_co2_flux(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%coex90(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%omex90(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%calex90(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%calex90(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%opex90(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%opex90(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%coex1000(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%omex1000(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%calex1000(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%calex1000(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%opex1000(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%opex1000(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%coex2000(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%omex2000(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%calex2000(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%calex2000(1), -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%opex2000(:,:), i_time_stat,&
& hamocc_state%p_tend%monitor%opex2000(1), -2)
CALL calc_inventory2d(p_patch_3d, ocean_state%p_prog(i_time_stat)%tracer(:,1,:,ialkali+no_tracer), &
  i_time_stat,hamocc_state%p_tend%monitor%sfalk(1),-2)
CALL calc_inventory2d(p_patch_3d, ocean_state%p_prog(i_time_stat)%tracer(:,1,:,isco212+no_tracer), &
  i_time_stat,hamocc_state%p_tend%monitor%sfdic(1),-2)
CALL calc_inventory2d(p_patch_3d, ocean_state%p_prog(i_time_stat)%tracer(:,1,:,iphosph+no_tracer), &
  i_time_stat,hamocc_state%p_tend%monitor%sfphos(1),-2)
CALL calc_inventory2d(p_patch_3d, ocean_state%p_prog(i_time_stat)%tracer(:,1,:,isilica+no_tracer), &
  i_time_stat,hamocc_state%p_tend%monitor%sfsil(1),-2)

CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwn2b(:,:,:), porwat, glob_pwn2b)
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_tend%n2budget(:,:,:), i_time_stat, glob_n2b,.TRUE.)

CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_tend%sedrn(:,:,:), porwat, &
& hamocc_state%p_tend%monitor%seddenit(1))

! Unit conversion 
hamocc_state%p_tend%monitor%phosy(1) = hamocc_state%p_tend%monitor%phosy(1) * p2gtc
hamocc_state%p_tend%monitor%phosy_cya(1) = hamocc_state%p_tend%monitor%phosy_cya(1) * p2gtc
hamocc_state%p_tend%monitor%grazing(1) = hamocc_state%p_tend%monitor%grazing(1) * p2gtc
hamocc_state%p_tend%monitor%remina(1) = hamocc_state%p_tend%monitor%remina(1) * p2gtc
hamocc_state%p_tend%monitor%exud(1) = hamocc_state%p_tend%monitor%exud(1) * p2gtc
hamocc_state%p_tend%monitor%exudz(1) = hamocc_state%p_tend%monitor%exudz(1) * p2gtc
hamocc_state%p_tend%monitor%zoomor(1) = hamocc_state%p_tend%monitor%zoomor(1) * p2gtc
hamocc_state%p_tend%monitor%phymor(1) = hamocc_state%p_tend%monitor%phymor(1) * p2gtc
hamocc_state%p_tend%monitor%graton(1) = hamocc_state%p_tend%monitor%graton(1) * p2gtc
hamocc_state%p_tend%monitor%bacfra(1) = hamocc_state%p_tend%monitor%bacfra(1) * p2gtc
hamocc_state%p_tend%monitor%bacfrac(1) = hamocc_state%p_tend%monitor%bacfrac(1) * p2gtc
hamocc_state%p_tend%monitor%net_co2_flux(1) = hamocc_state%p_tend%monitor%net_co2_flux(1) * c2gtc
hamocc_state%p_tend%monitor%delcar(1) = hamocc_state%p_tend%monitor%delcar(1) * c2gtc
hamocc_state%p_tend%monitor%wcdenit(1) = hamocc_state%p_tend%monitor%wcdenit(1) * nitdem* n2tgn
hamocc_state%p_tend%monitor%n2fix(1) = hamocc_state%p_tend%monitor%n2fix(1) * n2tgn * rn2
hamocc_state%p_tend%monitor%omex90(1) = hamocc_state%p_tend%monitor%omex90(1) * p2gtc
hamocc_state%p_tend%monitor%calex90(1) = hamocc_state%p_tend%monitor%calex90(1) * c2gtc
hamocc_state%p_tend%monitor%omex1000(1) = hamocc_state%p_tend%monitor%omex1000(1) * p2gtc
hamocc_state%p_tend%monitor%calex1000(1) = hamocc_state%p_tend%monitor%calex1000(1) * c2gtc
hamocc_state%p_tend%monitor%omex2000(1) = hamocc_state%p_tend%monitor%omex2000(1) * p2gtc
hamocc_state%p_tend%monitor%calex2000(1) = hamocc_state%p_tend%monitor%calex2000(1) * c2gtc
hamocc_state%p_tend%monitor%seddenit(1) = hamocc_state%p_tend%monitor%seddenit(1) * nitdem*n2tgn
hamocc_state%p_tend%monitor%cyaldoc(1) = hamocc_state%p_tend%monitor%cyaldet(1) * p2gtc * doccya_fac
hamocc_state%p_tend%monitor%cyaldet(1) = hamocc_state%p_tend%monitor%cyaldet(1) * p2gtc *(1._wp - doccya_fac)

! mean values of surface concentrations
hamocc_state%p_tend%monitor%sfalk(1) = hamocc_state%p_tend%monitor%sfalk(1)/totalarea
hamocc_state%p_tend%monitor%sfdic(1) = hamocc_state%p_tend%monitor%sfdic(1)/totalarea
hamocc_state%p_tend%monitor%sfsil(1) = hamocc_state%p_tend%monitor%sfsil(1)/totalarea
hamocc_state%p_tend%monitor%sfphos(1) = hamocc_state%p_tend%monitor%sfphos(1)/totalarea
hamocc_state%p_tend%monitor%zalkn2(1) = glob_n2b + glob_pwn2b

END SUBROUTINE get_monitoring

SUBROUTINE get_inventories(hamocc_state,ocean_state,p_patch_3d,i_time_stat)

USE mo_biomod,      ONLY: rnit,rn2, ro2bal,rcar

INTEGER :: i_time_stat
TYPE(t_hydro_ocean_state) :: ocean_state
TYPE(t_hamocc_state) :: hamocc_state
TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d

! Local variables
REAL(wp) :: glob_det, glob_doc, glob_phy, glob_zoo
REAL(wp) :: glob_phos, glob_sedo12, glob_bo12
REAL(wp) :: glob_sedc12, glob_nit, glob_gnit, glob_pwn2b, glob_pwh2ob
REAL(wp) :: glob_n2o,glob_n2fl,glob_n2ofl, glob_orginp
REAL(wp) :: glob_calinp, glob_silinp, glob_alk, glob_calc
REAL(wp) :: glob_sil, glob_opal, glob_sedsi, glob_pwsi
REAL(wp) :: glob_bsil, glob_silpro, glob_n2b, glob_h2ob
REAL(wp) :: glob_prorca,  glob_cya, glob_doccya,glob_produs
REAL(wp) :: glob_pwn2, glob_pwno3, glob_prcaca, glob_ofl
REAL(wp) :: glob_pwic, glob_pwal, glob_pwph, glob_cfl
REAL(wp) :: glob_dic, glob_o2, glob_fe, glob_co3, glob_hi
REAL(wp) :: glob_pwox, glob_pwfe, glob_bclay, glob_sedclay
REAL(wp) :: total_ocean, watersum, sedsum
REAL(wp) :: rcyano, glob_bc12
CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

cpara_name='======================='
cpara_val="==========="
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

!Check if cyanobacteria are calculcated
rcyano = MERGE(1._wp,0._wp, l_cyadyn)

! Calculate global inventories of individual tracers
! Water column

CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,idoc+no_tracer), &
&                     i_time_stat, glob_doc)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,idet+no_tracer), &
&                     i_time_stat, glob_det)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,iphy+no_tracer), &
&                     i_time_stat, glob_phy)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,izoo+no_tracer), &
&                     i_time_stat, glob_zoo)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,iphosph+no_tracer),&
&                     i_time_stat, glob_phos)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,iano3+no_tracer),&
&                     i_time_stat, glob_nit)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,igasnit+no_tracer),&
&                     i_time_stat, glob_gnit)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,ian2o+no_tracer), &
&                     i_time_stat, glob_n2o)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,isilica+no_tracer), &
&                     i_time_stat, glob_sil)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,icalc+no_tracer), &
&                     i_time_stat, glob_calc)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,iopal+no_tracer), &
&                     i_time_stat, glob_opal)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,ialkali+no_tracer),&
&                     i_time_stat, glob_alk)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,isco212+no_tracer),&
&                     i_time_stat, glob_dic)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,ioxygen+no_tracer),&
&                     i_time_stat, glob_o2)
CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,iiron+no_tracer),&
&                     i_time_stat, glob_fe)
IF(l_cyadyn)THEN
 CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,icya+no_tracer),&
&                      i_time_stat, glob_cya)
 CALL calc_inventory3d(p_patch_3d, ocean_state, ocean_state%p_prog(i_time_stat)%tracer(:,:,:,idoccya+no_tracer),&
&                      i_time_stat, glob_doccya)
else
 glob_cya=0._wp
 glob_doccya=0._wp
ENDIF
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_diag%co3(:,:,:), i_time_stat, glob_co3)
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_diag%hi(:,:,:), i_time_stat, glob_hi)

! Sediment
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%so12(:,:,:), porsol,  glob_sedo12)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%sc12(:,:,:), porsol,  glob_sedc12)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%ster(:,:,:), porsol,  glob_sedclay)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%ssil(:,:,:), porsol,  glob_sedsi)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwsi(:,:,:), porwat,  glob_pwsi)
!CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwsi(:,:,:), porsol, glob_pwsi)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwal(:,:,:), porwat,  glob_pwal)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwic(:,:,:), porwat,  glob_pwic)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwph(:,:,:), porwat,  glob_pwph)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwox(:,:,:), porwat,  glob_pwox)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwfe(:,:,:), porwat,  glob_pwfe)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwn2b(:,:,:), porwat, glob_pwn2b)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwn2(:,:,:), porwat,  glob_pwn2)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwno3(:,:,:), porwat,  glob_pwno3)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwh2ob(:,:,:), porwat, glob_pwh2ob)

CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bo12(:,:), i_time_stat, glob_bo12,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bsil(:,:), i_time_stat, glob_bsil,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bc12(:,:), i_time_stat, glob_bc12,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bter(:,:), i_time_stat, glob_bclay,-2)

! Tendencies
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%prorca(:,:), i_time_stat, glob_prorca,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%prcaca(:,:), i_time_stat, glob_prcaca,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%silpro(:,:), i_time_stat, glob_silpro,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%produs(:,:), i_time_stat, glob_produs,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%cflux(:,:), i_time_stat, glob_cfl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%oflux(:,:), i_time_stat, glob_ofl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%nflux(:,:), i_time_stat, glob_n2fl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%n2oflux(:,:), i_time_stat, glob_n2ofl, 1)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%orginp(:,:), i_time_stat, glob_orginp, 1, ocean_state)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%silinp(:,:), i_time_stat, glob_silinp, 1, ocean_state)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_acc%calinp(:,:), i_time_stat, glob_calinp, 1, ocean_state)
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_tend%h2obudget(:,:,:), i_time_stat, glob_h2ob,.TRUE.)
CALL calc_inventory3d(p_patch_3d, ocean_state, hamocc_state%p_tend%n2budget(:,:,:), i_time_stat, glob_n2b,.TRUE.)

! Convert unit for fluxes
![kmol/s] to [kmol]

glob_cfl=glob_cfl*dtbgc
glob_ofl=glob_ofl*dtbgc
glob_n2fl=glob_n2fl*dtbgc
glob_n2ofl=glob_n2ofl*dtbgc


! Print tracer output
CALL message(' ', ' ', io_stdo_bgc)
CALL message('Global inventory of', 'ocean tracers', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('DIC',glob_dic)
CALL to_bgcout('Alkalinity',glob_alk)
CALL to_bgcout('Phosphate',glob_phos)
CALL to_bgcout('Oxygen',glob_o2)
CALL to_bgcout('N2',glob_gnit)
CALL to_bgcout('Nitrate',glob_nit)
CALL to_bgcout('Silicate',glob_sil)
CALL to_bgcout('DOC',glob_doc)
CALL to_bgcout('Phytoplankton',glob_phy)
CALL to_bgcout('Zooplankton',glob_zoo)
CALL to_bgcout('N2O',glob_n2o)
CALL to_bgcout('Iron',glob_fe)
CALL to_bgcout('Detritus',glob_det)
IF(l_cyadyn)THEN
 CALL to_bgcout('Cyanobacteria',glob_cya)
 CALL to_bgcout('DOC cya',glob_doccya)
ENDIF
CALL to_bgcout('Calc',glob_calc)
CALL to_bgcout('Opal',glob_opal)

CALL message('Global inventory of', 'additional tracers', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('hi',glob_hi)
CALL to_bgcout('co3',glob_co3)

CALL message(' ', ' ', io_stdo_bgc)
CALL message('Global inventory of', 'aqueous sediment tracers', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('DIC',glob_pwic)
CALL to_bgcout('Alkalinity',glob_pwal)
CALL to_bgcout('Phosphate',glob_pwph)
CALL to_bgcout('Oxygen',glob_pwox)
CALL to_bgcout('N2',glob_pwn2)
CALL to_bgcout('Nitrate',glob_pwno3)
CALL to_bgcout('Silicate',glob_pwsi)
CALL to_bgcout('Iron',glob_pwfe)

CALL message(' ', ' ', io_stdo_bgc)
CALL message('Global inventory of', 'solid sediment constituents', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('Solid C org',glob_sedo12)
CALL to_bgcout('Burial C org',glob_bo12)
CALL to_bgcout('Solid CaCO3',glob_sedc12)
CALL to_bgcout('Burial CaCO3',glob_bc12)
CALL to_bgcout('Solid opal',glob_sedsi)
CALL to_bgcout('Burial opal',glob_bsil)
CALL to_bgcout('Solid clay',glob_sedclay)
CALL to_bgcout('Burial opal',glob_bclay)

CALL message(' ', ' ', io_stdo_bgc)
cpara_name='======================='
cpara_val="==========="
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('Redfield (global)',glob_nit/glob_phos)

CALL message(' ', ' ', io_stdo_bgc)
CALL message('Global fluxes into', 'atmosphere [kmol]', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('CO2 flux',glob_cfl)
CALL to_bgcout('O2 flux',glob_ofl)
CALL to_bgcout('N2 flux',glob_n2fl)
CALL to_bgcout('N2O flux',glob_n2ofl)

CALL message(' ', ' ', io_stdo_bgc)
CALL message('Global fluxes into', 'sediment [kmol]', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('prcaca',glob_prcaca)
CALL to_bgcout('prorca',glob_prorca)
CALL to_bgcout('silpro',glob_silpro)
CALL to_bgcout('produs',glob_produs)

CALL to_bgcout('zalkn2',glob_n2b+glob_pwn2b)

CALL message(' ', ' ', io_stdo_bgc)
CALL message('Global waethering fluxes', ' [kmol]', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('orginp',glob_orginp)
CALL to_bgcout('silinp',glob_silinp)
CALL to_bgcout('calcinp',glob_calinp)

cpara_name='======================='
cpara_val="==========="
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL message(' ', ' ', io_stdo_bgc)


! Calculate total inventory diagnostics
! DO not consider prcaca, prorca, silpro and produs, as the call is after the loop
! and the fluxes are zeros after powach, the output in bgcflux however is not equal zero
!-------- Phosphate
watersum = glob_det + glob_doc + glob_phy + glob_zoo + glob_phos  &
     &     + rcyano*(glob_cya + glob_doccya)
     
sedsum =  glob_sedo12 + glob_bo12 +  glob_orginp + glob_pwph 

total_ocean = watersum + sedsum


CALL to_bgcout('Global water phosphate [kmol]',watersum)
CALL to_bgcout('Global sediment phosphate [kmol]',sedsum)
CALL to_bgcout('Global total phosphate [kmol]',total_ocean)
CALL message(' ', ' ', io_stdo_bgc)


!-------- Nitrate
watersum = rnit * (glob_det + glob_doc + glob_phy + glob_zoo  &
     &     + rcyano*(glob_cya + glob_doccya) ) + glob_nit     &
     &     + rn2 * (glob_gnit + glob_n2o + glob_n2fl + glob_n2ofl)&
     &     + glob_pwn2 + glob_pwno3  
     
sedsum =  rnit* (glob_sedo12 + glob_bo12)   & 
     &     - glob_orginp

total_ocean = watersum + sedsum

CALL to_bgcout('Global water nitrate [kmol]',watersum)
CALL to_bgcout('Global sediment nitrate [kmol]',sedsum)
CALL to_bgcout('Global total nitrate [kmol]',total_ocean)
CALL message(' ', ' ', io_stdo_bgc)

!-------- Silicate 
watersum =  glob_sil + glob_opal + glob_pwsi
     
sedsum = glob_sedsi + glob_bsil  - glob_silinp

total_ocean = watersum + sedsum

CALL to_bgcout('Global water silicate [kmol]',watersum)
CALL to_bgcout('Global sediment silicate [kmol]',sedsum)
CALL to_bgcout('Global total silicate [kmol]',total_ocean)
CALL message(' ', ' ', io_stdo_bgc)

! Alkalinity

watersum = glob_alk - rnit* (glob_det + glob_doc + glob_phy + glob_zoo &
  &        + rcyano*(glob_doccya + glob_cya)) - glob_n2b          &
  &        + 2._wp * glob_calc + rnit * glob_orginp - 2._wp * glob_calinp

sedsum =  glob_sedc12 + glob_bc12  

total_ocean = watersum + sedsum

CALL to_bgcout('Global water alkalinity [kmol]',watersum)
CALL to_bgcout('Global sediment alkalinity [kmol]',sedsum)
CALL to_bgcout('Global total alkalinity [kmol]',total_ocean)
CALL message(' ', ' ', io_stdo_bgc)

! Oxygen

watersum = (glob_det + glob_doc + glob_phy + glob_zoo +         &
  &         rcyano*glob_cya + rcyano * glob_doccya)*(-ro2bal) + &
  &         glob_o2 + glob_phos*2._wp + glob_dic + glob_calc +  &
  &         glob_nit * 1.5_wp + glob_n2o* 0.5_wp + glob_pwno3* 1.5 + &
  &         glob_pwic + glob_pwox + glob_pwph*2._wp + glob_ofl + &
  &         glob_n2ofl * 0.5_wp + glob_cfl + glob_h2ob +         &
  &         glob_orginp*ro2bal - glob_calinp  

sedsum =   (glob_sedo12 + glob_bo12)*(-ro2bal) + glob_sedc12 + glob_bc12 

total_ocean = watersum + sedsum

CALL to_bgcout('Global water oxygen [kmol]',watersum)
CALL to_bgcout('Global sediment oxygen [kmol]',sedsum)
CALL to_bgcout('Global total oxygen [kmol]',total_ocean)
CALL message(' ', ' ', io_stdo_bgc)


! Carbon

watersum = (glob_det + glob_doc + glob_phy + glob_zoo   &
     &     + rcyano*(glob_cya + glob_doccya)) *rcar + &
     &     glob_dic + glob_calc - glob_calinp - rcar * glob_orginp + &
     &     glob_cfl
     
sedsum =  (glob_sedo12 + glob_bo12 ) * rcar + glob_sedc12 + glob_bc12 


total_ocean = watersum + sedsum


CALL to_bgcout('Global water carbon [kmol]',watersum)
CALL to_bgcout('Global sediment carbon [kmol]',sedsum)
CALL to_bgcout('Global total carbon [kmol]',total_ocean)
CALL message(' ', ' ', io_stdo_bgc)


cpara_name='======================='
cpara_name='======================='
cpara_val="==========="
CALL message(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

END SUBROUTINE
!--------------------------------------------------------------------------------------------
SUBROUTINE calc_inventory3d(patch3D, ocean_state,pfield3d, i_time_stat,field_globsum,no_thick)
! Calculate inventory of the whole water column
! inv = sum (tracer * volume)
! if no_thick is given:
! inv = sum(tracer * area) --> thickness needs to be cons. elsewhere
REAL(wp), TARGET:: pfield3d(:,:,:)
TYPE(t_patch_3d), TARGET, INTENT(in) :: patch3D
TYPE(t_hydro_ocean_state), TARGET:: ocean_state
REAL(wp), INTENT(OUT):: field_globsum
LOGICAL, OPTIONAL:: no_thick
INTEGER:: i_time_stat

! Local
INTEGER:: jk
TYPE(t_patch), POINTER :: patch_2d
REAL(wp) :: ptmp


patch_2d => patch3D%p_patch_2d(1)

field_globsum=0._wp

IF(PRESENT(no_thick))then 
! if multiplied elsewhere, do not consider cell thickness and surface area here
 DO jk=1,n_zlev
  ptmp = global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,jk,:) * &
&   pfield3d(:,jk,:) &
)
 field_globsum = field_globsum + ptmp 
 ENDDO
ELSE
 DO jk=1,n_zlev
  ptmp = global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,jk,:) * &
&   MERGE(patch3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,jk,:)+ocean_state%p_prog(i_time_stat)%h(:,:),&
&         patch3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,jk,:),jk==1) * &
&   pfield3d(:,jk,:) &
)
 field_globsum = field_globsum + ptmp 
 ENDDO

ENDIF

END SUBROUTINE 


SUBROUTINE calc_inventory_sed(patch3D, pfield3d, sed_state, field_globsum)
! calculated sediment inventory
! sed_state: porsol or porwat (for solid or pore water tracer)
! inv = sum (tracer * {porsol or porwat} * area * seddw)
REAL(wp), TARGET:: pfield3d(:,:,:)
REAL(wp) :: sed_state(ks)
TYPE(t_patch_3d), TARGET, INTENT(in) :: patch3D
REAL(wp), INTENT(OUT) :: field_globsum

! Local
INTEGER:: jk
TYPE(t_patch), POINTER :: patch_2d
REAL(wp) :: ptmp

patch_2d => patch3D%p_patch_2d(1)

field_globsum=0._wp
DO jk =1,ks
  ptmp = global_sum_array( &
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,1,:) * &
&   seddw(jk) * sed_state(jk) * &
&   pfield3d(:,jk,:) &
 )
 field_globsum = field_globsum + ptmp 
ENDDO

END SUBROUTINE 


SUBROUTINE calc_inventory2d(patch3D, pfield2d, i_time_stat,field_globsum, jk, ocean_state)
! calculate inventory of 2d fields (optionally at given level)
! inv = (tracer * volume)
! volume = area* (cell thickness  + surface_height) : if ocean_state given
! volume = area* (cell thickness(level) ) : if ocean_state not given
REAL(wp), TARGET:: pfield2d(:,:)
TYPE(t_patch_3d), TARGET, INTENT(in) :: patch3D
TYPE(t_hydro_ocean_state), TARGET, OPTIONAL:: ocean_state
REAL(wp), INTENT(OUT) :: field_globsum
INTEGER,  OPTIONAL :: jk
INTEGER:: i_time_stat
! Local
TYPE(t_patch), POINTER :: patch_2d
INTEGER :: ik

ik = 1
IF(PRESENT(jk)) ik=jk

patch_2d => patch3D%p_patch_2d(1)
IF(PRESENT(ocean_state))THEN ! if needed add surface_height to cell thickness
field_globsum=global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,1,:) * &
&   (patch3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:)+ocean_state%p_prog(i_time_stat)%h(:,:))*&
&   pfield2d(:,:))
ELSEIF(ik < 0) THEN ! do not consider thickness
field_globsum=global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,1,:) * &
&   pfield2d(:,:))
ELSE
field_globsum=global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,ik,:) * &
&   patch3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,ik,:)*&
&   pfield2d(:,:))
ENDIF


END SUBROUTINE 

END MODULE mo_hamocc_diagnostics
