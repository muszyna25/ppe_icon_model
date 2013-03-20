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
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, alv, rd, rgrav
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag 
  USE mo_nh_torus_exp,        ONLY: shflx_cbl, lhflx_cbl, set_sst_cbl
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                               spec_humi !,qsat_rho !! Specific humidity
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: surface_conditions, min_wind, rkarman

  REAL(wp), PARAMETER :: min_wind = 0.1_wp
  REAL(wp), PARAMETER :: rkarman = 2.5_wp

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
                                 p_prog_lnd_now, p_diag_lnd, prm_diag, theta)

    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    REAL(wp),          INTENT(in)        :: theta(:,:,:)  !pot temp  

    REAL(wp) :: rhos, th0_srf, obukhov_length, z_mc, ustar, mwind, bflux
    REAL(wp) :: zrough, pres_sfc(nproma,p_patch%nblks_c)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc
    INTEGER  :: nlev              

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

    !SELECT CASE(isrfc_type)

    !Prescribed latent/sensible heat fluxes: get ustar and surface temperature / moisture
    !CASE(1)
      
      rl_start = 2
      rl_end   = min_rlcell_int-1
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

      jk = nlev
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,th0_srf,bflux,mwind,ustar,z_mc, &
!$OMP            zrough,obukhov_length,rhos),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx

            !Roughness length
            zrough = prm_diag%gz0(jc,jb) * rgrav

            !Get reference surface temperature
            th0_srf = p_prog_lnd_now%t_g(jc,jb)

            !Buoyancy flux
            bflux = grav*(shflx_cbl+0.61_wp* th0_srf*lhflx_cbl)/th0_srf

            !Mean wind at nlev
            mwind  = MAX( min_wind, SQRT(p_nh_diag%u(jc,jk,jb)**2+p_nh_diag%v(jc,jk,jb)**2) )
           
            !Z height
            z_mc = p_nh_metrics%z_mc(jc,jk,jb)             

            !Now diagnose friction velocity (ustar)
            ustar = diag_ustar(z_mc,zrough,bflux,mwind)

            !"-" sign in the begining because ustar*thstar = -shflx
            obukhov_length = -ustar**3 * rkarman / bflux
 
            p_prog_lnd_now%t_g(jc,jb) = theta(jc,jk,jb) + shflx_cbl / ustar * &
                                        businger_heat(zrough,z_mc,obukhov_length) 

            !Get surface qv using t_g: saturation value
            p_diag_lnd%qv_s(jc,jb) =    &
                spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),pres_sfc(jc,jb)) 

            !Get surface fluxes
            rhos   =  pres_sfc(jc,jb)/(rd * p_prog_lnd_now%t_g(jc,jb))  
            prm_diag%shfl_s(jc,jb)  = shflx_cbl * rhos * cpd
            prm_diag%lhfl_s(jc,jb)  = lhflx_cbl * rhos * alv
            prm_diag%umfl_s(jc,jb)  = ustar**2  * rhos

         END DO  
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF(lhflx_cbl==0._wp) p_diag_lnd%qv_s(:,:) = 0._wp
  
    !CASE(2)
     !Fix SST type


   !END SELECT 


  END SUBROUTINE surface_conditions

  !>
  !! factor_heat
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile for neutral and unstable case for scalars
  !! Stable case is still to be done
  !! Beniot, On the integral of the surface layer profile-gradient functions (JAM), 1977
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  FUNCTION businger_heat(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, zeta0, lamda, lamda0

     IF(L > 0._wp)THEN !Stable

       !To be done

     ELSE !unstable or neutral
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(1._wp - 9._wp*zeta)  
       lamda0 = SQRT(1._wp - 9._wp*zeta0)  
        
       factor = 0.74_wp * ( LOG(z1/z0) + 2._wp*(LOG(1._wp+lamda0)-LOG(1._wp+lamda)) ) &
                * rkarman
     END IF 

  END FUNCTION businger_heat 

  !>
  !! factor_mom
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile for neutral and unstable case for velocities
  !! Beniot, On the integral of the surface layer profile-gradient functions (JAM), 1977
  !!------------------------------------------------------------------------
  FUNCTION businger_mom(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, zeta0, lamda, lamda0

     IF(L > 0._wp)THEN !Stable

       !To be done

     ELSE !unstable or neutral
       zeta   = z1/L 
       zeta0  = z0/L 
       lamda  = SQRT(SQRT(1._wp - 15._wp*zeta))  
       lamda0 = SQRT(SQRT(1._wp - 15._wp*zeta0))  
        
       factor = rkarman * ( LOG(z1/z0) + LOG((1._wp+lamda0**2)*(1._wp+lamda0)**2) - &
                            LOG((1._wp+lamda**2)*(1._wp+lamda)**2) + 2._wp * (      &
                            ATAN(lamda) - ATAN(lamda0) )  )
     END IF

  END FUNCTION businger_mom 

  !>
  !! diagnose ustar
  !!------------------------------------------------------------------------
  FUNCTION  diag_ustar(z1, z0, bflux, wind) RESULT(ustar)
     REAL(wp), INTENT(IN) :: z1, z0, bflux, wind

     REAL(wp) :: ustar, L
     INTEGER  :: ITERATE      

     !First guess       
     ustar = wind / (rkarman * LOG(z1/z0))

     IF(bflux > 0._wp)THEN   !Unstable case
       DO ITERATE = 1 , 4
          L = -ustar**3 * rkarman / bflux
          ustar = wind / businger_mom(z0,z1,L)
       END DO
     ELSEIF(bflux < 0._wp)THEN!Stable case
      
      !To be done

     END IF
        
  END FUNCTION diag_ustar
 
!-------------------------------------------------------------------------------
     
END MODULE mo_surface_les




