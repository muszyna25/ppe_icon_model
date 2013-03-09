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
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag 
  USE mo_nh_torus_exp,        ONLY: shflx_cbl, lhflx_cbl, set_sst_cbl, ufric_cbl

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: surface_conditions, min_wind, karman

  REAL(wp), PARAMETER :: min_wind = 0.1_wp
  REAL(wp), PARAMETER :: zrough  = 0.0003_wp, karman = 0.4_wp

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
  SUBROUTINE  surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_nh_prog, p_int, &
                                 p_prog_lnd_now, p_diag_lnd)

    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_nh_diag),   INTENT(inout)     :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog  !< single nh prognostic state
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 

    REAL(wp) :: flux_up, flux_dn, tend_hor, tend_ver, sflux, cpd_o_cvd
    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e), Rkarman
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc
    INTEGER  :: nlev              

    ! number of vertical levels
    nlev = p_patch%nlev
    i_nchdom   = MAX(1,p_patch%n_childdom)

    Rkarman = 1._wp / karman

    !Get mean wind at cell center at first level near surface
    !CALL rbf_vec_interpol_cell(p_nh_diag%vn, p_patch, p_int, u_nlev, v_nlev, opt_slev=nlev, &
    !                           opt_elev=nlev, opt_rlend=min_rlcell_int-1)

    SELECT CASE(isrfc_type)

    !Prescribed friction velocity and latend/sensible heat fluxes: get ustar, theta_start, etc
    !and surface temperature / moisture
    CASE(1)
      rl_start = 2
      rl_end   = min_rlcell_int-1
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

      jk = nlev
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,theta,obukhov_length,zeta,zeta0),ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx
            p_nh_diag%u_star(jc,jb)  = ufric_cbl            
            p_nh_diag%th_star(jc,jb) = -shflx_cbl / p_nh_diag%u_star(jc,jb)
            p_nh_diag%qv_star(jc,jb) = -lhflx_cbl / p_nh_diag%u_star(jc,jb)

            !Get surface temperature using theta at nlev            
            th0_srf  = p_nh_metrics%theta_ref_ic(jc,nlev+1,jb)

            theta = p_nh_prog%theta_v(jc,jk,jb) / (1._wp + 0.61_wp * &
                    p_nh_prog%tracer(jc,jb,jk,iqv) - p_nh_prog%tracer(jc,jk,jb,iqc))

            obukhov_length = p_nh_diag%u_star(jc,jb)**2 * th0_srf * Rkarman / &
                     (grav*(p_nh_diag%th_star(jc,jb)+0.61_wp*th0_srf*p_nh_diag%qv_star(jc,jb)))
 
            p_prog_lnd_now%t_g(jc,jb) = theta - p_nh_diag%th_star(jc,jb) *  & 
                 factor_unstable_heat(zrough,p_nh_metrics%z_mc(jc,nlev,jb),obukhov_length) 

            !Get surface qv using t_g: saturation value
            p_diag_lnd%qv_s(jc,jb) =    &
                spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_nh_diag%pres_sfc(jc,jb)) 
         END DO  
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF(lhflx_cbl==0._wp) p_diag_lnd%qv_s(:,:) = 0._wp
  
    CASE(2)
     !Fix SST type


   END SELECT 


  END SUBROUTINE surface_conditions

  !>
  !! factor_unstable_heat
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile for neutral and unstable case:
  !! Beniot, On the integral of the surface layer profile-gradient functions (JAM), 1977
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-06)
  ELEMENTAL FUNCTION  factor_unstable_heat(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, zeta0, lamda, lambda0

     zeta   = z1/L 
     zeta0  = z0/L 
     lamda  = SQRT(1._wp - 9._wp*zeta)  
     lamda0 = SQRT(1._wp - 9._wp*zeta0)  
        
     factor = 0.74_wp * ( LOG(z1/z0) + 2._wp*(LOG(1._wp+lamda0)-LOG(1._wp-lambda)) ) &
              * Rkarman

  END FUNCTION factor_unstable_heat 
  !>
  !! factor_unstable_mom
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile for neutral and unstable case:
  !! Beniot, On the integral of the surface layer profile-gradient functions (JAM), 1977
  !!------------------------------------------------------------------------
  ELEMENTAL FUNCTION  factor_unstable_mom(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
     REAL(wp) :: factor, zeta, zeta0, lamda, lambda0

     zeta   = z1/L 
     zeta0  = z0/L 
     lamda  = SQRT(SQRT(1._wp - 15._wp*zeta))  
     lamda0 = SQRT(SQRT(1._wp - 15._wp*zeta0))  
        
     factor = Rkarman * ( LOG(z1/z0) + LOG((1._wp+lamda0**2)*(1._wp+lambda0)**2) - &
                          LOG((1._wp+lamda**2)*(1._wp+lambda)**2) + 2._wp * (      &
                          ATAN(lamda) - ATAN(lambda0) )  )

  END FUNCTION factor_unstable_mom
   
!-------------------------------------------------------------------------------
     
END MODULE mo_surface_les




