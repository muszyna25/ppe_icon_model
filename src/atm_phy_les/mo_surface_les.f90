!>
!! mo_surface_les
!!
!! Surface calculations for les physics using Businger Dyer relationship
!! 
!! @author Anurag Dipankar, MPI-M
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-03-07)
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_surface_les

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish,message_text
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_intp,                ONLY: verts2edges_scalar, edges2verts_scalar, &
                                    cells2verts_scalar, cells2edges_scalar, &
                                    edges2cells_scalar, verts2cells_scalar, &
                                    edges2cells_vector
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: ltimer, msg_level
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array_mult, sync_patch_array, global_max, &
                                    global_sum_array
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, alv, rd, rgrav, rd_o_cpd, vtmpc1, rcpd
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag 
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_les_config,          ONLY: les_config
  USE mo_math_constants,      ONLY: pi_2, ln2, dbl_eps
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_data_turbdiff,       ONLY: akt, alpha0, z0_ice  
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water
  USE mo_util_phys,           ONLY: nwp_dyn_gust
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_nh_testcases_nml,    ONLY: th_cbl  
  USE mo_util_dbg_prnt,       ONLY: dbg_print

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: surface_conditions

  REAL(wp), PARAMETER :: km_min      = 0.01_wp !min viscosity

  !Parameters for surface layer parameterizations: From RB Stull's book
  REAL(wp), PARAMETER :: bsm = 4.7_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 15._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 4.7_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 9._wp   !Businger Untable Heat
  REAL(wp), PARAMETER :: Pr  = 0.74_wp !Km/Kh factor  
  REAL(wp), PARAMETER :: min_wind = 0.01_wp !min wind

  !Parameters for surface parameterizations from COSMO docs
  REAL(wp), PARAMETER :: beta_10 = 0.042_wp
  REAL(wp), PARAMETER :: h_10    = 10._wp
  REAL(wp), PARAMETER :: B       = 5._wp 
  REAL(wp), PARAMETER :: C       = 5._wp 
  REAL(wp), PARAMETER :: D       = 5._wp 
  REAL(wp), PARAMETER :: zh_max  = 0.1_wp

  !Parameters for RICO case
  REAL(wp), PARAMETER :: c_m = 0.001229_wp  
  REAL(wp), PARAMETER :: c_h = 0.001094_wp
  REAL(wp), PARAMETER :: c_q = 0.001133_wp
  REAL(wp), PARAMETER :: th0_rico = 298.5_wp
  REAL(wp), PARAMETER :: psfc = 101540._wp

  CHARACTER(len=12)  :: str_module = 'surface_les'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 4           ! Determines level of detail for 1 line debug

  CONTAINS


  !>
  !! surface_conditions
  !!------------------------------------------------------------------------
  !! Calculate surface temperature and moisture given the fluxes using Businger
  !! Dyer relationships .OR. vice versa. All calculations are done at cell center
  !!  
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  SUBROUTINE  surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_int, &
                                 p_prog_lnd_now, p_prog_lnd_new, p_diag_lnd, &
                                 prm_diag, theta, qv, visc_sfc_c)

    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(in)        :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_new!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    REAL(wp),          INTENT(in)        :: theta(:,:,:)  !pot temp  
    REAL(wp),          INTENT(in)        :: qv(:,:,:)     !spec humidity
    REAL(wp),          INTENT(inout)     :: visc_sfc_c(:,:)!apparent visc at surface

    REAL(wp) :: rhos, th0_srf, obukhov_length, z_mc, ustar, inv_mwind, mwind, bflux
    REAL(wp) :: zrough, exner, var(nproma,p_patch%nblks_c), theta_nlev, qv_nlev
    REAL(wp), POINTER :: pres_sfc(:,:)
    REAL(wp) :: theta_sfc, shfl, lhfl, umfl, vmfl, bflx1, bflx2, theta_sfc1, diff
    REAL(wp) :: RIB, zh, tcn_mom, tcn_heat
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, isidx, isblk, rl
    INTEGER :: nlev, jg, itr, jkp1
    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_surface_les:surface_conditions'

    IF (msg_level >= 15) &
         CALL message(TRIM(routine), '')

    jg = p_patch%id

    pres_sfc => p_nh_diag%pres_sfc

    ! number of vertical levels
    nlev = p_patch%nlev
    jk   = nlev
    jkp1 =  jk+1

    i_nchdom   = MAX(1,p_patch%n_childdom)
 
    !AD(18-02-2014):Set some minimal value for visc_sfc_c in 
    !the boundary zone. 
    visc_sfc_c  = km_min
 
    !loop indices: exclude halo 
    rl_start   = grf_bdywidth_c+1 
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)         


    SELECT CASE(les_config(jg)%isrfc_type)

    !Default case with TERRA
    !Most of what follows has been taken from mo_nwp_turbtrans_interface 
    !for COSMO turbulence

    CASE(1)

     !For now the LES scheme uses the surface fluxes from land directly. Exchange coefficients
     !are calculated by calling turbtrans using GME approach. See comment on mo_interface_les
     !where it is called. Status as on 11.09.2013 (AD)

     !Calculate visc_sfc_c for this case (03-12-2013)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,mwind,rhos,z_mc),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

           !Get apparent viscosity at the surface from the surface fluxes: to be used
           !in mo_sgs_turbulence
           
           !Mean wind at nlev
           mwind  = MAX( min_wind,SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )

           !Get surface fluxes
           !rho at surface: no qc at suface
           rhos   =  pres_sfc(jc,jb)/( rd * &
                     p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

           z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)
           visc_sfc_c(jc,jb) = rhos * ABS(prm_diag%tch(jc,jb)*mwind*z_mc* &
                                          les_config(jg)%turb_prandtl)
         END DO  
      END DO
!$OMP END DO 
!$OMP END PARALLEL

      !change sign of surface fluxes
      prm_diag%shfl_s  = - prm_diag%shfl_s
      prm_diag%lhfl_s  = - prm_diag%lhfl_s
      prm_diag%umfl_s  = - prm_diag%umfl_s
      prm_diag%vmfl_s  = - prm_diag%vmfl_s
 
    !Prescribed latent/sensible heat fluxes: get ustar and surface temperature / moisture
    CASE(2)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,exner,zrough,th0_srf,bflux,mwind,z_mc, &
!$OMP            RIB,ustar,obukhov_length,theta_sfc,rhos),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Surface exner
            exner = EXP( rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref) )

            !Roughness length
            zrough = prm_diag%gz0(jc,jb) * rgrav

            !Get reference surface pot. temperature
            !First time step t_g takes value assigned in nwp_phy_init          
            th0_srf = p_prog_lnd_now%t_g(jc,jb) / exner

            !Buoyancy flux
            bflux = grav*(les_config(jg)%shflx +  &
                   vtmpc1*th0_srf*les_config(jg)%lhflx)/th0_srf

            !Mean wind at nlev
            mwind  = MAX( min_wind,SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            !Z height
            z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

            !Now diagnose friction velocity (ustar)
            IF(les_config(jg)%ufric<0._wp)THEN
              !Bulk Richardson no at first model level
              RIB = grav * (theta(jc,jk,jb)-th0_srf) * (z_mc - zrough) / (th0_srf * mwind**2)
              ustar = SQRT( diag_ustar_sq(z_mc,zrough,RIB,mwind) )
            ELSE
              ustar = les_config(jg)%ufric
            END IF

            !"-" sign in the begining because ustar*thstar = -shflx
            obukhov_length = -ustar**3 / (akt * bflux)
             
            theta_sfc   = theta(jc,jk,jb) + les_config(jg)%shflx / ustar * &
                          businger_heat(zrough,z_mc,obukhov_length) 

            p_prog_lnd_new%t_g(jc,jb) = theta_sfc * exner

            !Get surface qv
            p_diag_lnd%qv_s(jc,jb) = qv(jc,jk,jb) + les_config(jg)%lhflx / ustar * &
                                     businger_heat(zrough,z_mc,obukhov_length) 

            !Get surface fluxes
            !rho at surface: no qc at suface
            rhos   =  pres_sfc(jc,jb)/( rd * &
                      p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

            prm_diag%shfl_s(jc,jb)  = les_config(jg)%shflx * rhos * cpd
            prm_diag%lhfl_s(jc,jb)  = les_config(jg)%lhflx * rhos * alv
            prm_diag%umfl_s(jc,jb)  = ustar**2  * rhos * p_nh_diag%u(jc,jk,jb) / mwind
            prm_diag%vmfl_s(jc,jb)  = ustar**2  * rhos * p_nh_diag%v(jc,jk,jb) / mwind

           !Get apparent viscosity at the surface from the surface fluxes: to be used
           !in mo_sgs_turbulence
           visc_sfc_c(jc,jb) = rhos*ABS(les_config(jg)%shflx)*z_mc* &
                               les_config(jg)%turb_prandtl/ABS(theta(jc,jk,jb)-theta_sfc) 

         END DO  
      END DO
!$OMP END DO 
!$OMP END PARALLEL

    !Prescribed buoyancy flux and transfer coefficient at surface to get a uniform SST (Stevens 2007 JAS)
    !It uses fixed transfer coefficient and assumes that q_s is saturated
    !It assumes that surface temperature is same as theta at surface
    CASE(3)

      !Get mean theta and qv at first model level and mean surface pressure

      var(:,:) = theta(:,jk,:)
      WHERE(.NOT.p_patch%cells%decomp_info%owner_mask(:,:)) var(:,:) = 0._wp
      theta_nlev =  global_sum_array(var)/REAL(p_patch%n_patch_cells_g,wp)

      var(:,:) = qv(:,jk,:)
      WHERE(.NOT.p_patch%cells%decomp_info%owner_mask(:,:)) var(:,:) = 0._wp
      qv_nlev =  global_sum_array(var)/REAL(p_patch%n_patch_cells_g,wp)


      !Iterate to get surface temperature given buoyancy flux:following UCLA-LES

      !Note that t_g(:,:) is uniform so the first index is used below
      rl    = grf_bdywidth_c+1
      isblk = p_patch%cells%start_blk(rl,1)
      isidx = p_patch%cells%start_idx(rl,1)
      theta_sfc = p_prog_lnd_now%t_g(isidx,isblk)
      
      diff = 1._wp
      itr = 0 
      DO WHILE (diff > 1._wp-6 .AND. itr < 30)
         bflx1 = les_config(jg)%tran_coeff*( (theta_sfc-theta_nlev)+vtmpc1* &
                 theta_nlev*(spec_humi(sat_pres_water(theta_sfc),p0ref)- &
                 qv_nlev) )*grav/th_cbl(1)

         theta_sfc1 = theta_sfc + 0.01_wp 

         bflx2 = les_config(jg)%tran_coeff*( (theta_sfc1-theta_nlev)+vtmpc1* &
                 theta_nlev*(spec_humi(sat_pres_water(theta_sfc1),p0ref)- &
                 qv_nlev) )*grav/th_cbl(1)

         theta_sfc = theta_sfc1 + 0.01_wp*(les_config(jg)%bflux-bflx1)/(bflx2-bflx1) 

         diff = ABS(1._wp - theta_sfc/theta_sfc1)
         itr = itr + 1         
      END DO               
       
      !WRITE(message_text,'(i4,f14.6,f14.7)')itr,diff,bflx2
      !CALL message('FINAL ITR, RESID, AND BFLX:',message_text )

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,zrough,z_mc,mwind,RIB,ustar,rhos),ICON_OMP_RUNTIME_SCHEDULE 
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Mean wind at nlev
            mwind  = MAX( min_wind,SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

            !Now diagnose friction velocity (ustar)
            IF(les_config(jg)%ufric<0._wp)THEN
              !Roughness length
              zrough = prm_diag%gz0(jc,jb) * rgrav

              !Bulk Richardson no at first model level
              RIB = grav * (theta(jc,jk,jb)-theta_sfc) * (z_mc-zrough)/(theta_sfc*mwind**2)
              ustar = SQRT( diag_ustar_sq(z_mc,zrough,RIB,mwind) )
            ELSE
              ustar = les_config(jg)%ufric
            END IF

            !Surface temperature
            p_prog_lnd_new%t_g(jc,jb) = theta_sfc 

            !Get surface qv 
            p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(p_prog_lnd_new%t_g(jc,jb)),p0ref)

            !Get surface fluxes
            rhos   =  p0ref/( rd * &
                      p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

            prm_diag%shfl_s(jc,jb) = rhos*cpd*les_config(jg)%tran_coeff*(theta_sfc-theta(jc,jk,jb))
            prm_diag%lhfl_s(jc,jb) = rhos*alv*les_config(jg)%tran_coeff*(p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
            prm_diag%umfl_s(jc,jb)  = ustar**2  * rhos * p_nh_diag%u(jc,jk,jb) / mwind
            prm_diag%vmfl_s(jc,jb)  = ustar**2  * rhos * p_nh_diag%v(jc,jk,jb) / mwind

           !Get apparent viscosity at the surface from the surface fluxes: to be used
           !in mo_sgs_turbulence
           visc_sfc_c(jc,jb) = rhos*ABS(les_config(jg)%tran_coeff)*z_mc*&
                               les_config(jg)%turb_prandtl

         END DO  
      END DO
!$OMP END DO 
!$OMP END PARALLEL

   !Rico case
    CASE(4)

      !RICO case - bulk aerodynamic formulation with fixed exchange coef.
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
                              
            !Get surface qv and temperature
            p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(les_config(jg)%sst),psfc)
            p_prog_lnd_new%t_g(jc,jb) = les_config(jg)%sst
            
            !Mean wind at nlev
            mwind  = SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) 
 
            shfl  =   c_h * mwind * (th0_rico-theta(jc,jk,jb))
            lhfl  =   c_q * mwind * (p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
            umfl  =   c_m * mwind * p_nh_diag%u(jc,jk,jb)
            vmfl  =   c_m * mwind * p_nh_diag%v(jc,jk,jb)
                                 
            !Surface density 
            rhos   =  psfc/( rd * &
                      p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  
                      
            !Get surface fluxes                       
            prm_diag%shfl_s(jc,jb)  = - shfl * rhos * cpd
            prm_diag%lhfl_s(jc,jb)  = - lhfl * rhos * alv
            prm_diag%umfl_s(jc,jb)  = umfl * rhos
            prm_diag%vmfl_s(jc,jb)  = vmfl * rhos

            z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

           !Get apparent viscosity at the surface from the surface fluxes: to be used
           !in mo_sgs_turbulence
           visc_sfc_c(jc,jb) = rhos*ABS(shfl)*z_mc* &
                               les_config(jg)%turb_prandtl/ABS(theta(jc,jk,jb)-th0_rico) 

        END DO  
      END DO
 
    !Fix SST case
    CASE(5)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,zrough,theta_sfc,mwind,z_mc, &
!$OMP            RIB,tcn_mom,tcn_heat,zh,rhos),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

           !Roughness length
           zrough = prm_diag%gz0(jc,jb) * rgrav

           !Get surface pot. temperature and humidity
           p_prog_lnd_new%t_g(jc,jb) = les_config(jg)%sst
           theta_sfc = p_prog_lnd_new%t_g(jc,jb) / EXP( rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref) )

           p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(les_config(jg)%sst),pres_sfc(jc,jb))

           !Mean wind at nlev
           mwind  = MAX( min_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
          
           !Z height to be used as a reference height in surface layer
           z_mc   = p_nh_metrics%z_mc(jc,jk,jb) - p_nh_metrics%z_ifc(jc,jkp1,jb)

           !Bulk Richardson no at first model level
           RIB = grav * (theta(jc,jk,jb)-theta_sfc) * (z_mc-zrough) / (theta_sfc * mwind**2)

           !Momentum transfer coefficient
           tcn_mom             = (akt/LOG(z_mc/zrough))**2
           prm_diag%tcm(jc,jb) = tcn_mom * stability_function_mom(RIB,z_mc/zrough,tcn_mom)

           !Heat transfer coefficient
           zh = MIN(zrough,zh_max)
           tcn_heat            = akt**2/(LOG(z_mc/zrough)*LOG(z_mc/zh))
           prm_diag%tch(jc,jb) = tcn_heat * stability_function_heat(RIB,z_mc/zh,tcn_heat)

           !Get surface fluxes
           !rho at surface: no qc at suface
           rhos   =  pres_sfc(jc,jb)/( rd * &
                     p_prog_lnd_new%t_g(jc,jb)*(1._wp+vtmpc1*p_diag_lnd%qv_s(jc,jb)) )  

           prm_diag%shfl_s(jc,jb)  = rhos*cpd*prm_diag%tch(jc,jb)*mwind*(theta_sfc-theta(jc,jk,jb))
           prm_diag%lhfl_s(jc,jb)  = rhos*alv*prm_diag%tch(jc,jb)*mwind*(p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
           prm_diag%umfl_s(jc,jb)  = rhos*prm_diag%tcm(jc,jb)*mwind*p_nh_diag%u(jc,jk,jb) 
           prm_diag%vmfl_s(jc,jb)  = rhos*prm_diag%tcm(jc,jb)*mwind*p_nh_diag%v(jc,jk,jb) 

           !Get apparent viscosity at the surface from the surface fluxes: to be used
           !in mo_sgs_turbulence
           visc_sfc_c(jc,jb) = rhos*ABS(prm_diag%tch(jc,jb))*mwind*z_mc*les_config(jg)%turb_prandtl
           
         END DO
      END DO   
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !No fluxes
    CASE(0)
      prm_diag%shfl_s  = 0._wp
      prm_diag%lhfl_s  = 0._wp 
      prm_diag%umfl_s  = 0._wp 
      prm_diag%vmfl_s  = 0._wp 
      visc_sfc_c       = km_min

  END SELECT 


  !Sync is required for mom fluxes and visc_sfc_c
  CALL sync_patch_array(SYNC_C, p_patch, visc_sfc_c)
  CALL sync_patch_array(SYNC_C, p_patch, prm_diag%umfl_s)
  CALL sync_patch_array(SYNC_C, p_patch, prm_diag%vmfl_s)


  END SUBROUTINE surface_conditions

  !>
  !! factor_heat
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile for neutral and unstable case for scalars
  !! Stable case is still to be done
  !! Beniot, On the integral of the surface layer profile-gradient functions (JAM), 1977
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  FUNCTION businger_heat(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, lamda, psi
     REAL(wp) :: zeta0, lamda0, psi0

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L 
       zeta0  = z0/L 
       psi    = -bsh*zeta
       psi    = -bsh*zeta0
       factor = Pr * (LOG(z1/z0) - psi + psi0) / akt
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(1._wp - buh*zeta)  
       lamda0 = SQRT(1._wp - buh*zeta0)  
       psi    = 2._wp * ( LOG(1._wp+lamda) - ln2 )
       psi0   = 2._wp * ( LOG(1._wp+lamda0) - ln2 )
       factor = Pr * (LOG(z1/z0) - psi + psi0) / akt
     ELSE !Neutral
       factor = Pr * LOG(z1/z0) / akt
     END IF 

  END FUNCTION businger_heat 
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  FUNCTION phi_heat(z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z1, L
     REAL(wp) :: factor, zeta, lamda

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L 
       lamda  = bsh*zeta
       factor = Pr + lamda
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       lamda  = SQRT(1._wp - buh*zeta)  
       factor = Pr / lamda
     ELSE !neutral
       factor = Pr 
     END IF 

  END FUNCTION phi_heat 

  !>
  !! factor_mom
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile for neutral and unstable case for velocities
  !! Beniot, On the integral of the surface layer profile-gradient functions (JAM), 1977
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  FUNCTION businger_mom(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, psi, lamda
     REAL(wp) :: zeta0, psi0, lamda0

     IF(L > 0._wp)THEN !Stable
       zeta  = z1/L 
       zeta0 = z0/L 
       psi  = -bsm*zeta
       psi0 = -bsm*zeta0

       factor = ( LOG(z1/z0) - psi + psi0 ) / akt
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(SQRT(1._wp - bum*zeta))  
       lamda0 = SQRT(SQRT(1._wp - bum*zeta0))  

       psi    = 2._wp * LOG(1._wp+lamda) + LOG(1._wp+lamda*lamda) - &
                2._wp * ATAN(lamda) + pi_2 - 3._wp*ln2

       psi0   = 2._wp * LOG(1._wp+lamda0) + LOG(1._wp+lamda0*lamda0) - &
                2._wp * ATAN(lamda0) + pi_2 - 3._wp*ln2

       factor = ( LOG(z1/z0) - psi + psi0 ) / akt
     ELSE !neutral
       factor = LOG(z1/z0) / akt
     END IF

  END FUNCTION businger_mom 
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  FUNCTION phi_mom(z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z1, L
     REAL(wp) :: factor, zeta, lamda

     IF(L > 0._wp)THEN !Stable
       zeta   = z1/L 
       factor = 1._wp + bsm * zeta
     ELSEIF(L < 0._wp)THEN !unstable 
       zeta   = z1/L 
       lamda  = SQRT(SQRT(1._wp - bum*zeta))  
       factor = 1._wp / lamda
     ELSE !neutral
       factor  = 1._wp
     END IF

  END FUNCTION phi_mom 

  !>
  !! diagnose ustar
  !!------------------------------------------------------------------------
  FUNCTION diag_ustar_sq(h, z0, RIB, wind) RESULT(ustar_sq)
     REAL(wp), INTENT(IN) :: h, z0, RIB, wind

     REAL(wp) :: ustar_sq, tcn_mom

     tcn_mom  = (akt/LOG(h/z0))**2
     ustar_sq = tcn_mom*wind**2*stability_function_mom(RIB,h/z0,tcn_mom)
        
  END FUNCTION diag_ustar_sq
 
  !>
  !! diagnose wstar
  !!------------------------------------------------------------------------
  FUNCTION diag_wstar_sq(h, z0, RIB, wind) RESULT(wstar_sq)
     REAL(wp), INTENT(IN) :: h, z0, RIB, wind

     REAL(wp) :: wstar_sq, tcn_heat, zh
  
     zh = MIN(z0, zh_max)

     IF(RIB < 0._wp)THEN
       tcn_heat = akt**2/(LOG(h/z0)*LOG(h/zh))
       wstar_sq = ( tcn_heat * wind**3 * ABS(RIB) * &
                    stability_function_heat(RIB,h/zh,tcn_heat) )**(2._wp/3._wp)
     ELSE
       wstar_sq = 0._wp
     END IF      
        
  END FUNCTION diag_wstar_sq
 
  !>
  !! init_zrough
  !!------------------------------------------------------------------------
  FUNCTION init_gz(h, wind) RESULT(gz)
     REAL(wp), INTENT(IN) :: h, wind

     REAL(wp) :: z01, z02, gz
 

     z01 = alpha0 * wind**2 / ( (1._wp/beta_10) + LOG(h/h_10)/akt )**2
     !z02 = ( alpha0 * wind**2 * ABS(RIB) )**1.5_wp / (C * SQRT(grav*h))

     !gz = MAX(z01, z02) 
     gz = z01
        
  END FUNCTION init_gz

  !>
  !! stability_function_mom
  !! Taken from COSMO docs
  !!------------------------------------------------------------------------
  FUNCTION stability_function_mom(RIB, hz0, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hz0, tc

     REAL(wp) :: stab_fun, hz0_fac
 
     IF(RIB.GE.0._wp)THEN
       stab_fun = 1._wp / ( 1._wp + 2._wp*B*RIB/SQRT(1._wp+D*RIB) ) 
     ELSE
       hz0_fac = ( hz0**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (z0/h)**(1/3)<<1 giving hz0_fac=SQRT(h/z0)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 2._wp*B*ABS(RIB)/(1._wp + 3._wp*B*C*tc*hz0_fac*SQRT(ABS(RIB)))
     END IF 
        
  END FUNCTION stability_function_mom
  !>
  !! stability_function_heat
  !!------------------------------------------------------------------------
  FUNCTION stability_function_heat(RIB, hzh, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hzh, tc

     REAL(wp) :: stab_fun, hzh_fac
 
     IF(RIB.GE.0._wp)THEN
       stab_fun = 1._wp / ( 1._wp + 3._wp*B*RIB*SQRT(1._wp+D*RIB) ) 
     ELSE
       hzh_fac = ( hzh**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (zh/h)**(1/3)<<1 giving hzh_fac=SQRT(h/zh)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 2._wp*B*ABS(RIB)/(1._wp + 3._wp*B*C*tc*hzh_fac*SQRT(ABS(RIB)))
     END IF 
  END FUNCTION stability_function_heat

!-------------------------------------------------------------------------------

    
END MODULE mo_surface_les




