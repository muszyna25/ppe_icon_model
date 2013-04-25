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
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, global_max
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, alv, rd, rgrav, rd_o_cpd
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag 
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_les_config,          ONLY: les_config
  USE mo_math_constants,      ONLY: pi_2, ln2

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: surface_conditions, min_wind

  REAL(wp), PARAMETER :: min_wind = 0.1_wp
 
  !Parameters for surface layer parameterizations: From RB Stull's book
  REAL(wp), PARAMETER :: bsm = 4.7_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 15._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 4.7_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 9._wp   !Businger Untable Heat
  REAL(wp), PARAMETER :: Pr  = 0.74_wp !Km/Kh factor

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
                                 p_prog_lnd_now, p_diag_lnd, prm_diag,    &
                                 theta, theta_sfc, qv)

    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    REAL(wp),          INTENT(in)        :: theta(:,:,:)  !pot temp  
    REAL(wp),          INTENT(out)       :: theta_sfc(:,:)!sfc pot temp  
    REAL(wp),          INTENT(in)        :: qv(:,:,:)     !spec humidity

    REAL(wp) :: rhos, th0_srf, obukhov_length, z_mc, ustar, mwind, bflux
    REAL(wp) :: zrough, pres_sfc(nproma,p_patch%nblks_c), exner
    REAL(wp) :: theta_sfc_new, th_grad, qv_grad, theta_ic, qv_ic, phi
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc
    INTEGER :: nlev, jg, itr
    
    jg = p_patch%id

    ! number of vertical levels
    nlev = p_patch%nlev
    i_nchdom   = MAX(1,p_patch%n_childdom)
     
    !Pres_sfc from previous calls are not synced. Therefore, I sync it here locally
    IF (p_test_run) THEN
       pres_sfc(:,:) = 0._wp
    ENDIF

!$OMP PARALLEL WORKSHARE
    pres_sfc(:,:) = p_nh_diag%pres_sfc(:,:)
!$OMP END PARALLEL WORKSHARE

    CALL sync_patch_array(SYNC_C, p_patch, pres_sfc)

    SELECT CASE(les_config(jg)%isrfc_type)

    !Fix SST type
    CASE(1)

    !To be implemented

    !Prescribed latent/sensible heat fluxes: get ustar and surface temperature / moisture
    CASE(2)
      
      rl_start = 2
      rl_end   = min_rlcell_int-1
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

      jk = nlev
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,th0_srf,bflux,mwind,ustar,z_mc, &
!$OMP phi,th_grad,theta_ic,exner,zrough,obukhov_length,rhos),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Surface exner
            exner = EXP( rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref) )

            !Roughness length
            zrough = prm_diag%gz0(jc,jb) * rgrav

            !Get reference surface pot. temperature
            th0_srf = p_prog_lnd_now%t_g(jc,jb) / exner

            !Buoyancy flux
            bflux = grav*(les_config(jg)%shflx +  &
                   0.61_wp*th0_srf*les_config(jg)%lhflx)/th0_srf

            !Mean wind at nlev
            mwind  = MAX( min_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            !Z height
            z_mc = p_nh_metrics%z_mc(jc,jk,jb)             

            !Now diagnose friction velocity (ustar)
            IF(les_config(jg)%ufric<0._wp)THEN
              ustar = diag_ustar(z_mc,zrough,bflux,mwind)
            ELSE
              ustar = les_config(jg)%ufric
            END IF

            !"-" sign in the begining because ustar*thstar = -shflx
            obukhov_length = -ustar**3 * les_config(jg)%rkarman_constant / bflux
 
            !Following step is only a fix to get correct gradient for calculating
            !reasonable kh:

            !Get surface theta to satisfy the flux profile: set gradient
            phi     = phi_heat(z_mc,obukhov_length)

            th_grad = -les_config(jg)%shflx*phi*les_config(jg)%rkarman_constant/(ustar*z_mc)          
            
            theta_ic  = p_nh_metrics%wgtfac_c(jc,jk,jb) * theta(jc,jk,jb) + &
                        (1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb)) * theta(jc,jk-1,jb)

            theta_sfc(jc,jb) = theta_ic - th_grad * p_nh_metrics%ddqz_z_full(jc,jk,jb)

            !theta_sfc(jc,jb) = theta(jc,jk,jb) + les_config(jg)%shflx / ustar * &
            !                   businger_heat(zrough,z_mc,obukhov_length) 

            p_prog_lnd_now%t_g(jc,jb) = theta_sfc(jc,jb) * exner

            !Get surface qv normally
            p_diag_lnd%qv_s(jc,jb) = qv(jc,jk,jb) + les_config(jg)%lhflx / ustar * &
                                     businger_heat(zrough,z_mc,obukhov_length) 

            !Get surface fluxes
            rhos   =  pres_sfc(jc,jb)/( rd * &
                      p_prog_lnd_now%t_g(jc,jb)*(1._wp+0.61_wp*p_diag_lnd%qv_s(jc,jb)) )  

            prm_diag%shfl_s(jc,jb)  = les_config(jg)%shflx * rhos * cpd
            prm_diag%lhfl_s(jc,jb)  = les_config(jg)%lhflx * rhos * alv
            prm_diag%umfl_s(jc,jb)  = ustar**2  * rhos
          
         END DO  
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !Prescribed buoyancy flux and transfer coefficient at surface - Stevens 2007 JAS
    !It uses fixed transfer coefficient and assumes that q_s is saturated 
    CASE(3)

      rl_start = 2
      rl_end   = min_rlcell_int-1
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

      jk = nlev
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,exner,zrough,th0_srf,theta_sfc_new,itr, &
!$OMP  mwind,z_mc,ustar,obukhov_length,th_grad,theta_ic,rhos),ICON_OMP_RUNTIME_SCHEDULE 
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Surface exner
            exner = EXP( rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref) )

            !Roughness length
            zrough = prm_diag%gz0(jc,jb) * rgrav

            !Get reference surface pot. temperature
            th0_srf = p_prog_lnd_now%t_g(jc,jb) / exner

            !Iterate to get surface temperature given buoyancy flux
           
            !First guess
            theta_sfc(jc,jb) = th0_srf
            DO itr = 1 , 4
              theta_sfc_new =    ( th0_srf * rgrav * les_config(jg)%bflux - &
                  0.61_wp * th0_srf * les_config(jg)%tran_coeff *           &
                 (spec_humi(sat_pres_water(theta_sfc(jc,jb)),pres_sfc(jc,jb)) - &
                  qv(jc,jk,jb)) ) / les_config(jg)%tran_coeff + theta(jc,jk,jb)

              theta_sfc(jc,jb) = theta_sfc_new
            END DO               

            !Mean wind at nlev
            mwind  = MAX( min_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            !Z height
            z_mc = p_nh_metrics%z_mc(jc,jk,jb)             

            !Now diagnose friction velocity (ustar)
            IF(les_config(jg)%ufric<0._wp)THEN
              ustar = diag_ustar(z_mc,zrough,les_config(jg)%bflux,mwind)
            ELSE
              ustar = les_config(jg)%ufric
            END IF

            !------------------------------------------------------------------
            !Following step is only a fix to get correct gradient for calculating
            !reasonable kh:
            !Overwrite surface theta to satisfy the flux profile using the flux
            !calculated using theta_sfc calculated above
            !------------------------------------------------------------------

            !"-" sign in the begining because ustar*thstar = -shflx
            obukhov_length = -ustar**3 * les_config(jg)%rkarman_constant / les_config(jg)%bflux
 
            th_grad = -les_config(jg)%tran_coeff*(theta_sfc(jc,jb)-theta(jc,jk,jb))* &
                       phi_heat(z_mc,obukhov_length)*les_config(jg)%rkarman_constant/(ustar*z_mc)          
            
            theta_ic  = p_nh_metrics%wgtfac_c(jc,jk,jb) * theta(jc,jk,jb) + &
                        (1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb)) * theta(jc,jk-1,jb)

            theta_sfc(jc,jb) = theta_ic - th_grad * p_nh_metrics%ddqz_z_full(jc,jk,jb)
              
            !Surface temperature
            p_prog_lnd_now%t_g(jc,jb) = theta_sfc(jc,jb) * exner

            !Get surface qv 
            p_diag_lnd%qv_s(jc,jb) = spec_humi(sat_pres_water(theta_sfc(jc,jb)),pres_sfc(jc,jb))

            !Get surface fluxes
            rhos   =  pres_sfc(jc,jb)/( rd * &
                      p_prog_lnd_now%t_g(jc,jb)*(1._wp+0.61_wp*p_diag_lnd%qv_s(jc,jb)) )  
            prm_diag%shfl_s(jc,jb) = rhos*cpd*les_config(jg)%tran_coeff*(theta_sfc(jc,jb)-theta(jc,jk,jb))
            prm_diag%lhfl_s(jc,jb) = rhos*alv*les_config(jg)%tran_coeff*(p_diag_lnd%qv_s(jc,jb)-qv(jc,jk,jb))
            prm_diag%umfl_s(jc,jb) = ustar**2 *rhos

         END DO  
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   END SELECT 


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
       factor = Pr * (LOG(z1/z0) - psi + psi0) * les_config(1)%rkarman_constant
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(1._wp - buh*zeta)  
       lamda0 = SQRT(1._wp - buh*zeta0)  
       psi    = 2._wp * ( LOG(1._wp+lamda) - ln2 )
       psi0   = 2._wp * ( LOG(1._wp+lamda0) - ln2 )
       factor = Pr * (LOG(z1/z0) - psi + psi0) * les_config(1)%rkarman_constant
     ELSE !Neutral
       factor = Pr * LOG(z1/z0) * les_config(1)%rkarman_constant
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
       factor = Pr / SQRT(lamda)
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

       factor = les_config(1)%rkarman_constant * ( LOG(z1/z0) - psi + psi0 )
     ELSEIF(L < 0._wp)THEN !unstable
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(SQRT(1._wp - bum*zeta))  
       lamda0 = SQRT(SQRT(1._wp - bum*zeta0))  

       psi    = 2._wp * LOG(1._wp+lamda) + LOG(1._wp+lamda*lamda) - &
                2._wp * ATAN(lamda) + pi_2 - 3._wp*ln2

       psi0   = 2._wp * LOG(1._wp+lamda0) + LOG(1._wp+lamda0*lamda0) - &
                2._wp * ATAN(lamda0) + pi_2 - 3._wp*ln2

       factor = les_config(1)%rkarman_constant * ( LOG(z1/z0) - psi + psi0 )
     ELSE !neutral
       factor = les_config(1)%rkarman_constant * LOG(z1/z0) 
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
  FUNCTION  diag_ustar(z1, z0, bflux, wind) RESULT(ustar)
     REAL(wp), INTENT(IN) :: z1, z0, bflux, wind

     REAL(wp) :: ustar, L
     INTEGER  :: ITERATE      

     !First guess       
     ustar = wind * les_config(1)%karman_constant / LOG(z1/z0)

     IF(bflux > 0._wp)THEN   !Unstable case
       DO ITERATE = 1 , 4
          L = -ustar**3 * les_config(1)%rkarman_constant / bflux
          ustar = wind / businger_mom(z0,z1,L)
       END DO
     ELSEIF(bflux < 0._wp)THEN!Stable case
      
      !To be done

     END IF
        
  END FUNCTION diag_ustar
 
!-------------------------------------------------------------------------------
     
END MODULE mo_surface_les




