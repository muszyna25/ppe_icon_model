!>
!! mo_sgs_turbulence
!!
!! Calculates 3D subgrid-scale viscosity and diffusivity in the nonhydrostatic model
!! 
!! @author Anurag Dipankar, MPI-M
!!
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
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

! USE of #ifdef __LOOP_EXCHANGE is MISSING

MODULE mo_sgs_turbulence

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish,message_text, debug_messages_on
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_intp,                ONLY: verts2edges_scalar, edges2verts_scalar, &
                                    cells2verts_scalar, cells2edges_scalar, &
                                    edges2cells_scalar, verts2cells_scalar, &
                                    edges2cells_vector
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: iqv, iqc, msg_level
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_math_constants,      ONLY: dbl_eps, pi
  USE mo_math_utilities,      ONLY: tdma_solver
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
                                    sync_patch_array_mult
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, rcpd, alv
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag 
  USE mo_surface_les,         ONLY: min_wind, surface_conditions 
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_les_config,          ONLY: les_config
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_mpi,                 ONLY: p_pe

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  REAL(wp), PARAMETER :: km_min = 0.01_wp

  PUBLIC :: drive_subgrid_diffusion

  !Variables for the module
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: visc_smag_v, visc_smag_ie, diff_smag_e, &
                                             visc_smag_c, rho_e, DIV_c, u_vert, v_vert, &
                                             w_ie, w_vert
  
  CONTAINS

  !-------------------------------------------------------------------------------------
  !>
  !! drive_subgrid_diffusion
  !!------------------------------------------------------------------------
  !! Driver for computing the sgs viscosity and diffusivity using Smagorisnky model
  !! and evaluating diffusion term on triangles
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-03-05)
  SUBROUTINE drive_subgrid_diffusion(p_nh_prog, p_nh_prog_rcf, p_nh_diag, p_nh_metrics, p_patch, &
                                     p_int, p_prog_lnd_now, p_prog_lnd_new, p_diag_lnd, prm_diag,&
                                     prm_nwp_tend, dt)

    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog     !< single nh prognostic state
    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog_rcf !< rcf nh prognostic state 
    TYPE(t_nh_diag),   INTENT(inout)     :: p_nh_diag     !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch       !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int         !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(in)        :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_new!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
    REAL(wp),          INTENT(in)        :: dt

    REAL(wp), ALLOCATABLE :: theta(:,:,:), thetav(:,:,:)
    REAL(wp) :: visc_sfc_c(nproma,1,p_patch%nblks_c)

    INTEGER :: nlev, nlevp1
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, jg

    !CALL debug_messages_on

    jg = p_patch%id

    nlev   = p_patch%nlev
    nlevp1 = nlev+1
    i_nchdom   = MAX(1,p_patch%n_childdom)
 
   
    ALLOCATE( u_vert(nproma,nlev,p_patch%nblks_v),           &
              v_vert(nproma,nlev,p_patch%nblks_v),           &
              w_vert(nproma,nlevp1,p_patch%nblks_v),         &
              w_ie(nproma,nlevp1,p_patch%nblks_e),           &
              visc_smag_v(nproma,nlev,p_patch%nblks_v),      &
              visc_smag_c(nproma,nlev,p_patch%nblks_c),      &
              visc_smag_ie(nproma,nlevp1,p_patch%nblks_e),   &                              
              diff_smag_e(nproma,nlev,p_patch%nblks_e),      &
              theta(nproma,nlev,p_patch%nblks_c),            &
              thetav(nproma,nlev,p_patch%nblks_c),           &
              DIV_c(nproma,nlev,p_patch%nblks_c),            &
              rho_e(nproma,nlev,p_patch%nblks_e)             &
             )

    !Initialize
    IF(p_test_run)THEN
!ICON_OMP_WORKSHARE
      u_vert(:,:,:) = 0._wp; v_vert(:,:,:) = 0._wp; w_vert(:,:,:) = 0._wp 
!ICON_OMP_END_WORKSHARE
    END IF

    !Convert temperature/tempv to potential/thetav temperature: all routines within 
    !use theta and thetav
    !Note that tracers are NOT synced till this stage, therefore satad and hence
    !temperature, pressure are only calculated for interior nodes.
    !Syncing things that are needed locally is better than syncing tracers (ntracers!)
    !and it also avoids a lot of if-endif in nwp_phy_interface
   
    !Sync temp and relevant tracers here: tracer syncing is required for tracers
    !calculation in horizontal diffusion
    !-for now only qv and qc are needed and hence are synced here. Add more if needed
    !-also sync exncer because it is changed after calling satad_v_3D
    CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_nh_diag%temp, p_nh_diag%tempv, &
                               p_nh_prog%exner, p_nh_prog%tracer(:,:,:,iqv),        &
                               p_nh_prog%tracer(:,:,:,iqc))

  
    !Calculate theta (for diffusion) and theta_v for Ri calculation
    rl_start   = 2
    rl_end     = min_rlcell_int-2
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jc = i_startidx, i_endidx
         theta(jc,1:nlev,jb) = p_nh_diag%temp(jc,1:nlev,jb) / &
                               p_nh_prog%exner(jc,1:nlev,jb)

         thetav(jc,1:nlev,jb) = p_nh_diag%tempv(jc,1:nlev,jb) / &
                                p_nh_prog%exner(jc,1:nlev,jb)
       END DO
    END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    !Think about moving this call to mo_nh_interface_nwp where nwp_surface is called    
    CALL surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_int,     &
                            p_prog_lnd_now, p_prog_lnd_new, p_diag_lnd,  &
                            prm_diag, theta, p_nh_prog%tracer(:,:,:,iqv),&
                            visc_sfc_c(:,1,:))

    CALL smagorinsky_model(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag, &
                           thetav, visc_sfc_c)

    CALL diffuse_hori_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag, &
                               prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb, dt)

    !Vertical velocity is updated here
    CALL diffuse_vert_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, dt)

    CALL diffuse_scalar(theta, p_nh_metrics, p_patch, p_int, p_nh_diag,  &
                        prm_nwp_tend%ddt_temp_turb, p_nh_prog%exner,     &
                        prm_diag, p_nh_prog%rho, dt, 'theta')

    !For qv and qc
    IF(.NOT.les_config(jg)%is_dry_cbl)THEN
      CALL diffuse_scalar(p_nh_prog%tracer(:,:,:,iqv), p_nh_metrics, p_patch, p_int, &
                          p_nh_diag, prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv),        &
                          p_nh_prog%exner, prm_diag, p_nh_prog%rho, dt, 'qv')

      CALL diffuse_scalar(p_nh_prog%tracer(:,:,:,iqc), p_nh_metrics, p_patch, p_int, &
                          p_nh_diag, prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc),        &
                          p_nh_prog%exner, prm_diag, p_nh_prog%rho, dt, 'qc')
    ELSE 
!ICON_OMP_WORKSHARE
      prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc) = 0._wp
!ICON_OMP_END_WORKSHARE
    END IF

   DEALLOCATE(u_vert, v_vert, w_vert, w_ie, visc_smag_v, visc_smag_ie, diff_smag_e, &
              theta, visc_smag_c, rho_e, DIV_c)


  END SUBROUTINE drive_subgrid_diffusion
  !-------------------------------------------------------------------------------------

 
  !-------------------------------------------------------------------------------------
  !>
  !! smagorisnky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorisnky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j+du_j/dx_i
  !!
  !! and  K = cs*\Delta**2*D/sqrt(2) where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2du_1/dx_1
  !!        D_22 = 2d_u2/dx_2
  !!        D_33 = 2d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
  SUBROUTINE smagorinsky_model(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                                prm_diag, thetav, visc_sfc_c)

    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    REAL(wp),          INTENT(in)        :: thetav(:,:,:)! virtual potential temperature
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    REAL(wp),          INTENT(in)        :: visc_sfc_c(:,:,:)!surface sgs visc

    ! local variables
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: DD, div_of_stress, visc_smag_e, &
                                               vn_ie, vt_ie, thetav_e
    REAL(wp) :: visc_sfc_e(nproma,1,p_patch%nblks_e), z1, z2, d_thetav_ie
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, vt_vert1, vt_vert2, vt_vert3, &
                vt_vert4, w_full_c1
    REAL(wp) :: w_full_c2, w_full_v1, w_full_v2, brunt_vaisala_frq
    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33
    REAL(wp), POINTER :: diff_smag_ic(:,:,:)

    INTEGER  :: nlev, nlevp1             !< number of full levels
    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, je, ic, jkp1, jv, jkm1, jcn, jbn, jvn

    !--------------------------------------------------------------------------
    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = nlev+1
    i_nchdom   = MAX(1,p_patch%n_childdom)

    !Allocation
    ALLOCATE( vn_ie(nproma,nlevp1,p_patch%nblks_e),        &
              vt_ie(nproma,nlevp1,p_patch%nblks_e),        &
              DD(nproma,nlev,p_patch%nblks_e),             &
              visc_smag_e(nproma,nlev,p_patch%nblks_e),    &
              thetav_e(nproma,nlev,p_patch%nblks_e),       &
              div_of_stress(nproma,nlev,p_patch%nblks_e)   &
            )

    !Initialize
    IF(p_test_run)THEN
!ICON_OMP_WORKSHARE
      visc_smag_e(:,:,:) = 0._wp
!ICON_OMP_END_WORKSHARE
    END IF 

    diff_smag_ic => prm_diag%tkvh

    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------

    CALL cells2verts_scalar(p_nh_prog%w, p_patch, p_int%cells_aw_verts, w_vert, opt_rlend=min_rlvert_int)
    CALL sync_patch_array(SYNC_V, p_patch, w_vert)

    CALL cells2edges_scalar(p_nh_prog%w, p_patch, p_int%c_lin_e, w_ie, opt_rlend=min_rledge_int-2)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int, &
                                  u_vert, v_vert, opt_rlend=min_rlvert_int )
    CALL sync_patch_array(SYNC_V, p_patch, u_vert)
    CALL sync_patch_array(SYNC_V, p_patch, v_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation
    rl_start = 2
    rl_end   = min_rledge_int-3

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)
    
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2, nlev
         DO je = i_startidx, i_endidx
            vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb)*p_nh_prog%vn(je,jk,jb) &
                    + (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb))*p_nh_prog%vn(je,jk-1,jb)         
         END DO
       END DO
       DO je = i_startidx, i_endidx
          vn_ie(je,1,jb) =                              &
            p_nh_metrics%wgtfacq1_e(je,1,jb)*p_nh_prog%vn(je,1,jb) + &
            p_nh_metrics%wgtfacq1_e(je,2,jb)*p_nh_prog%vn(je,2,jb) + &
            p_nh_metrics%wgtfacq1_e(je,3,jb)*p_nh_prog%vn(je,3,jb)

          vn_ie(je,nlevp1,jb) =                           &
            p_nh_metrics%wgtfacq_e(je,1,jb)*p_nh_prog%vn(je,nlev,jb) +   &
            p_nh_metrics%wgtfacq_e(je,2,jb)*p_nh_prog%vn(je,nlev-1,jb) + &
            p_nh_metrics%wgtfacq_e(je,3,jb)*p_nh_prog%vn(je,nlev-2,jb)
       END DO
    END DO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
 
    CALL rbf_vec_interpol_edge(vn_ie, p_patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    rl_start = 4
    rl_end   = min_rledge_int-2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,  &
!ICON_OMP            vt_vert1,vt_vert2,vt_vert3,vt_vert4,w_full_c1,w_full_c2,w_full_v1,  &
!ICON_OMP            w_full_v2,D_11,D_12,D_13,D_22,D_23,D_33,jkp1)
    DO jb = i_startblk,i_endblk

       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1, nlev
         DO je = i_startidx, i_endidx
          
            jkp1 = jk + 1

            vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v2

            vt_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v2 

            vt_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v2

            vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v2

            vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v2

            vt_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v2 

            vt_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v2

            !W at full levels
            w_full_c1  = 0.5_wp * (                                        &
                         p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) + &
                         p_nh_prog%w(iecidx(je,jb,1),jkp1,iecblk(je,jb,1)) ) 
 
            w_full_c2  = 0.5_wp * (                                        &
                         p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) + &
                         p_nh_prog%w(iecidx(je,jb,2),jkp1,iecblk(je,jb,2)) ) 

            !W at full levels vertices from w at vertices at interface levels
            w_full_v1  = 0.5_wp * (                                    &
                         w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) +    &
                         w_vert(ividx(je,jb,1),jkp1,ivblk(je,jb,1)) ) 

            w_full_v2  = 0.5_wp * (                                    &
                         w_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) +    &
                         w_vert(ividx(je,jb,2),jkp1,ivblk(je,jb,2)) ) 

            !Strain rates at edge center
            D_11       =      2._wp * (vn_vert4-vn_vert3) * &       
                              p_patch%edges%inv_vert_vert_length(je,jb)

            D_12       =     p_patch%edges%system_orientation(je,jb) *  &
                            (vn_vert2-vn_vert1) *                       &
                             p_patch%edges%inv_primal_edge_length(je,jb)&
                             + (vt_vert4-vt_vert3)*                     &
                             p_patch%edges%inv_vert_vert_length(je,jb)   
                    
            D_13       =     (vn_ie(je,jk,jb)-vn_ie(je,jkp1,jb))* &
                              p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +   &
                             (w_full_c2 - w_full_c1) *                      &
                              p_patch%edges%inv_dual_edge_length(je,jb)   

            D_22       =     2._wp*(vt_vert2-vt_vert1)*p_patch%edges%system_orientation(je,jb) * &
                             p_patch%edges%inv_primal_edge_length(je,jb)

            D_23       =    (vt_ie(je,jk,jb) - vt_ie(je,jkp1,jb)) *          &
                             p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +     &
                             p_patch%edges%system_orientation(je,jb) *       &
                            (w_full_v2 - w_full_v1) *                        &
                             p_patch%edges%inv_primal_edge_length(je,jb)   

            D_33       =     2._wp * (w_ie(je,jk,jb) - w_ie(je,jkp1,jb)) *   &
                             p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)               
           
            DD(je,jk,jb)   = D_11**2 + D_22**2 + D_33**2  + 2._wp * ( D_12**2 + D_13**2 + D_23**2 )             

            !calculate divergence to get the deviatoric part of stress tensor in 
            !diffusion: D_11-1/3*(D_11+D_22+D_33)
            div_of_stress(je,jk,jb)  = 0.5_wp * ( D_11 + D_22 + D_33 )
 
         ENDDO
       ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL


    !div_of_stress from edge to cell-scalar interpolation
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int-1 !-1 for its use in hor diffusion

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           DIV_c(jc,jk,jb) =                                                               &
              (div_of_stress(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               div_of_stress(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               div_of_stress(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))
         END DO
       END DO
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    !--------------------------------------------------------------------------
    !3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !--------------------------------------------------------------------------

    ! 3(a)Calculate thetav at the edges for gradient Richardson number and 
    ! some additional interpolations. 
    
    CALL cells2edges_scalar(thetav, p_patch, p_int%c_lin_e, thetav_e,  &
                            opt_rlstart=4, opt_rlend= min_rledge_int)

    CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, rho_e, &
                            opt_rlstart=4, opt_rlend=min_rledge_int-1)

    !3(b) Calculate stability corrected turbulent viscosity
    ! visc = mixing_length_sq * SQRT(DD/2) * SQRT(1-Ri/Pr) where Ri = (g/thetav)*d_thetav_dz/(DD/2), 
    ! where Brunt_vaisala_freq = (g/thetav)*d_thetav_dz. After simplification: 
    ! visc = mixing_length_sq * SQRT[ DD/2 - (Brunt_vaisala_frq/Pr) ]

    CALL cells2edges_scalar(visc_sfc_c, p_patch, p_int%c_lin_e, visc_sfc_e, opt_slev=1, &
                            opt_elev=1, opt_rlstart=4, opt_rlend=min_rledge_int)


    rl_start = 4
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx,d_thetav_ie,brunt_vaisala_frq,z1,z2)
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev-1
         DO je = i_startidx, i_endidx

           d_thetav_ie = p_nh_metrics%wgtfac_e(je,jk,jb) * thetav_e(je,jk,jb) +       &
                     (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb)) * thetav_e(je,jk-1,jb) &
                     - (p_nh_metrics%wgtfac_e(je,jk+1,jb) * thetav_e(je,jk+1,jb) +    &
                     (1._wp - p_nh_metrics%wgtfac_e(je,jk+1,jb)) * thetav_e(je,jk,jb))

           brunt_vaisala_frq = grav * d_thetav_ie * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) &
                               / thetav_e(je,jk,jb)

           visc_smag_e(je,jk,jb) = rho_e(je,jk,jb) *                  &
               MAX( km_min, p_nh_metrics%mixing_length_sq(je,jk,jb) * &
               SQRT(MAX(0._wp, DD(je,jk,jb)*0.5_wp-les_config(1)%rturb_prandtl*brunt_vaisala_frq)) ) 

         END DO
       END DO
      DO je = i_startidx, i_endidx
        z1  = 1._wp / p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb)
        z2  = p_nh_metrics%ddqz_z_full_e(je,nlev,jb) * 0.5_wp
        visc_smag_e(je,nlev,jb) =  (visc_sfc_e(je,1,jb)*z1+visc_smag_e(je,nlev-1,jb)*z2)/(z1+z2)
        visc_smag_e(je,1,jb)    = rho_e(je,1,jb) * km_min
      END DO
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    CALL sync_patch_array(SYNC_E, p_patch, visc_smag_e)

    !--------------------------------------------------------------------------
    !4) Interpolate viscosity and strain rates to different locations: calculate them for 
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !4a) visc from edge to cell-scalar interpolation
    rl_start = grf_bdywidth_c
    rl_end   = min_rlcell_int-1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           visc_smag_c(jc,jk,jb) =                                                       &
              (visc_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               visc_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               visc_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))           
         END DO
       END DO
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    !4b) visc at vertices
    CALL edges2verts_scalar(visc_smag_e, p_patch, p_int%e_aw_v, visc_smag_v, &
                            opt_rlstart=5, opt_rlend=min_rlvert_int-1) 

   ! CALL cells2verts_scalar(visc_smag_c, p_patch, p_int%cells_aw_verts, visc_smag_v, &
   !                         opt_rlstart=5, opt_rlend=min_rlvert_int-1) 

    !4c) Now calculate visc_smag at half levels at edge. Boundary values not required

!ICON_OMP_PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int 

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev
         DO je = i_startidx, i_endidx
           visc_smag_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * visc_smag_e(je,jk,jb) + &
                       (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb)) * visc_smag_e(je,jk-1,jb)          
         END DO
       END DO     
    END DO      
!ICON_OMP_END_DO


    !--------------------------------------------------------------------------
    !5) Calculate turbulent diffusivity
    !--------------------------------------------------------------------------

    !Turbulent diffusivity at edge at full levels
    rl_start = grf_bdywidth_e
    rl_end   = min_rledge_int-1

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO je = i_startidx, i_endidx
           diff_smag_e(je,jk,jb) = visc_smag_e(je,jk,jb) * les_config(1)%rturb_prandtl
         END DO
       END DO
    END DO
!ICON_OMP_END_DO

    !Turbulent diffusivity at cell center at half levels and like visc_smag_ie 
    !it is also calculated for interior points only and boundary values not required
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!ICON_OMP_DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
           diff_smag_ic(jc,jk,jb) = ( p_nh_metrics%wgtfac_c(jc,jk,jb)*visc_smag_c(jc,jk,jb) + &
                                  (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb))*visc_smag_c(jc,jk-1,jb) &
                                  ) * les_config(1)%rturb_prandtl   
         END DO
       END DO     
    END DO      
!ICON_OMP_END_DO

    !Copy to prm_diag%tkvm for output: find another way to do it
!ICON_OMP_DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)

       !Boundary values or tkvh
       DO jc = i_startidx, i_endidx
         prm_diag%tkvh(jc,1,jb)     = prm_diag%tkvh(jc,2,jb) 
         prm_diag%tkvh(jc,nlevp1,jb)= 0._wp
       END DO
    
       !now calculate tkvm
       DO jk = 1 , nlevp1
         DO jc = i_startidx, i_endidx
           prm_diag%tkvm(jc,jk,jb) = prm_diag%tkvh(jc,jk,jb) * les_config(1)%turb_prandtl 
         END DO
       END DO     
    END DO      
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    !DEALLOCATE variables
    DEALLOCATE( DD, div_of_stress, visc_smag_e, vn_ie, vt_ie, thetav_e )
  
  END SUBROUTINE smagorinsky_model
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_hori_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for normal velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - The Euler scheme is generally known to be unconditionaly unstable for diffusion term
  !!   but it is not certain when implemented in time split manner
  !! - This needs to be investigated and if found unstable we will switch to MacKormack 
  !!   scheme which is known to be stable
  !!  
  !! d_vn/d_t =  d_tau_11/d_x1 + d_tau_12/d_x2 + d_tau_13/d_x3
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_hori_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                                   prm_diag, ddt_u, ddt_v, dt)

    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nwp_phy_diag),  INTENT(in)    :: prm_diag     !< atm phys vars
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_u(:,:,:) !< u tendency
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_v(:,:,:) !< v tendency
    REAL(wp),           INTENT(in)       :: dt           !< dt turb

    REAL(wp) :: flux_up_e, flux_dn_e, flux_up_v, flux_dn_v, flux_up_c, flux_dn_c
    REAL(wp) :: stress_uc, stress_vc, stress_c1n, stress_c2n, inv_mwind
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt
    REAL(wp) :: inv_rhoe(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: vn_new(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: unew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: vnew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jkp1, jb, je, jcn, jbn, jvn, jc
    INTEGER  :: nlev              

    ! number of vertical levels
    nlev     = p_patch%nlev
    i_nchdom = MAX(1,p_patch%n_childdom)
   
    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk
 
    !Some initializations

!ICON_OMP_WORKSHARE
    !total tendency
    tot_tend(:,:,:) = 0._wp

    !new vn
    vn_new(:,:,:) = p_nh_prog%vn(:,:,:)
!ICON_OMP_END_WORKSHARE


    !Inverse of density (global rho_e in this module)
    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
         inv_rhoe(je,jk,jb) = 1._wp / rho_e(je,jk,jb)
       END DO
      END DO
    END DO
!ICON_OMP_END_DO

    ! 1) First get the horizontal tendencies

!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!ICON_OMP   jcn,jbn,flux_up_c,flux_dn_c,jvn,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
 
         vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                    p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                    v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                    p_patch%edges%primal_normal_vert(je,jb,1)%v2

         vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                    p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                    v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                    p_patch%edges%primal_normal_vert(je,jb,2)%v2

         vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                    v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%primal_normal_vert(je,jb,3)%v2

         vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                    v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%primal_normal_vert(je,jb,4)%v2

         dvt      = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                    v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%dual_normal_vert(je,jb,4)%v2 - &
                    (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                    v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%dual_normal_vert(je,jb,3)%v2)
 

         !tendency in normal direction:         
         !flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

         jcn     = iecidx(je,jb,2)  
         jbn     = iecblk(je,jb,2)  
         flux_up_c = visc_smag_c(jcn,jk,jbn) * &
                   ( 4._wp*(vn_vert4-p_nh_prog%vn(je,jk,jb))   * &
                     p_patch%edges%inv_vert_vert_length(je,jb) - &
                     0.6666666_wp*DIV_c(jcn,jk,jbn) )

         jcn     = iecidx(je,jb,1)  
         jbn     = iecblk(je,jb,1)  
         flux_dn_c = visc_smag_c(jcn,jk,jbn) * &
                   ( 4._wp*(p_nh_prog%vn(je,jk,jb)-vn_vert3)   * &
                     p_patch%edges%inv_vert_vert_length(je,jb) - &
                     0.6666666_wp*DIV_c(jcn,jk,jbn) )

         !tendency in tangential direction

         !D_12 between edge center and the vertex: delta_v/(primal_edge_len/2) + 
         ! ((vt4+vt2)/2-(vt3+vt2)/2)/(distance_opp_edges)
         !flux = D_12*visc
          
         !Note that the tangential velocities at vertices are used in D_12 is an 
         !approximation for speed. Better way is to use vt reconsructed from vn at 
         !each edges. Also, visc at somewhere between edge mid point and the vertex 
         !should be used but this is a good approximation

         jvn     = ividx(je,jb,2)
         jbn     = ivblk(je,jb,2)
         flux_up_v = visc_smag_v(jvn,jk,jbn) * ( p_patch%edges%system_orientation(je,jb) *       & 
           (vn_vert2-p_nh_prog%vn(je,jk,jb))*p_patch%edges%inv_primal_edge_length(je,jb)*2._wp + &
            dvt*p_patch%edges%inv_vert_vert_length(je,jb) )  

         jvn     = ividx(je,jb,1)
         jbn     = ivblk(je,jb,1)
         flux_dn_v = visc_smag_v(jvn,jk,jbn) * ( p_patch%edges%system_orientation(je,jb)      *  & 
           (p_nh_prog%vn(je,jk,jb)-vn_vert1)*p_patch%edges%inv_primal_edge_length(je,jb)*2._wp + &
            dvt*p_patch%edges%inv_vert_vert_length(je,jb) )  


         tot_tend(je,jk,jb) = ( (flux_up_c-flux_dn_c)*p_patch%edges%inv_dual_edge_length(je,jb) + &
                        p_patch%edges%system_orientation(je,jb) * (flux_up_v-flux_dn_v)  * &
                        p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp ) * inv_rhoe(je,jk,jb)

       END DO
      END DO
    END DO                         
!ICON_OMP_END_DO

   ! 2) Vertical tendency

!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx,flux_up_e,flux_dn_e)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 2 , nlev-1
       DO je = i_startidx, i_endidx

         !tendency in vertical direction
         flux_up_e = visc_smag_ie(je,jk,jb) *                             &
                ( (p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb)) * &
                   p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) +           &
                  (p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                   p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) *   &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                   
         flux_dn_e = visc_smag_ie(je,jk+1,jb) *                           &
                ( (p_nh_prog%vn(je,jk,jb) - p_nh_prog%vn(je,jk+1,jb)) * &
                   p_nh_metrics%inv_ddqz_z_half_e(je,jk+1,jb) +         &
                  (p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) -  &
                   p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1))) * &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                  
        tot_tend(je,jk,jb) = tot_tend(je,jk,jb) +  (flux_up_e - flux_dn_e) * &
                          p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

       END DO
      END DO
    END DO  
!ICON_OMP_END_DO

   !-----------------------------------------------------------------
   !3) Boundary treatment in vertical: use surface fluxes from 
   !   surface_conditions
   !-----------------------------------------------------------------
   
   !-----------------------------------------------------------------
   !jk = 1
   !-----------------------------------------------------------------

!ICON_OMP_DO PRIVATE(jb,je,i_startidx,i_endidx,flux_up_e,flux_dn_e)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
       DO je = i_startidx, i_endidx

         flux_up_e = 0._wp
                   
         flux_dn_e = visc_smag_ie(je,2,jb) *                           &
                ( (p_nh_prog%vn(je,1,jb) - p_nh_prog%vn(je,2,jb)) * &
                   p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) +         &
                  (p_nh_prog%w(iecidx(je,jb,2),2,iecblk(je,jb,2)) -  &
                   p_nh_prog%w(iecidx(je,jb,1),2,iecblk(je,jb,1))) * &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                  
         tot_tend(je,1,jb) = tot_tend(je,1,jb) + (flux_up_e - flux_dn_e) *  &
                    p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)

       END DO
    END DO                         
!ICON_OMP_END_DO

   !-----------------------------------------------------------------
   ! jk = nlev
   !-----------------------------------------------------------------

!ICON_OMP_DO PRIVATE(jb,je,i_startidx,i_endidx,flux_up_e,flux_dn_e,stress_c1n,stress_c2n)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
       DO je = i_startidx, i_endidx

         flux_up_e = visc_smag_ie(je,nlev,jb) *                             &
                ( (p_nh_prog%vn(je,nlev-1,jb) - p_nh_prog%vn(je,nlev,jb)) * &
                   p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) +             &
                  (p_nh_prog%w(iecidx(je,jb,2),nlev,iecblk(je,jb,2)) -      &
                   p_nh_prog%w(iecidx(je,jb,1),nlev,iecblk(je,jb,1))) *     &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                  
         !Get net shear stress in the direction of vn at surface

         !shear stress in normal direction from cell 1 
         stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                      p_patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                      prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                      p_patch%edges%primal_normal_cell(je,jb,1)%v2

         !shear stress in normal direction from cell 2 
         stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                      p_patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                      prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                      p_patch%edges%primal_normal_cell(je,jb,2)%v2

         !Net stress at the edge
         flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + &
                        stress_c2n * p_int%c_lin_e(je,2,jb) 

         tot_tend(je,nlev,jb) = tot_tend(je,nlev,jb) + (flux_up_e - flux_dn_e) * &
              p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

       END DO
    END DO                         
!ICON_OMP_END_DO

   ! 4) Update vn: it makes more sense to first apply diffusion on vn
   ! and then get ddt_u/v than to interpolating tot_tend directly to 
   ! ddt_u/v. Proof: during the test phase it was found that the latter slowed
   ! down the computation by 15%, although the results look same.

!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
         vn_new(je,jk,jb) = p_nh_prog%vn(je,jk,jb) + dt * tot_tend(je,jk,jb) 
       END DO
      END DO
    END DO  
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    !5) Get turbulent tendency at cell center
    CALL sync_patch_array(SYNC_E, p_patch, vn_new)
    CALL rbf_vec_interpol_cell(vn_new, p_patch, p_int, unew, vnew, opt_rlend=min_rlcell_int)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
   
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1 , nlev
       DO jc = i_startidx, i_endidx
         ddt_u(jc,jk,jb) = ( unew(jc,jk,jb) - p_nh_diag%u(jc,jk,jb) ) / dt   
         ddt_v(jc,jk,jb) = ( vnew(jc,jk,jb) - p_nh_diag%v(jc,jk,jb) ) / dt   
       END DO
      END DO
    END DO  
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

!    CALL sync_patch_array(SYNC_E, p_patch, tot_tend)
!    CALL rbf_vec_interpol_cell(tot_tend, p_patch, p_int, ddt_u, ddt_v, opt_rlend=min_rlcell_int-1)

  END SUBROUTINE diffuse_hori_velocity
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_vert_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for vertical velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - The Euler scheme is generally known to be unconditionaly unstable for diffusion term
  !!   but it is not certain when implemented in time split manner
  !! - This needs to be investigated and if found unstable we will switch to MacKormack 
  !!   scheme which is known to be stable
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_vert_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, dt)

    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    REAL(wp),          INTENT(in)        :: dt

    REAL(wp) :: flux_up_c, flux_dn_c, dvn1, dvn2, dvt1, dvt2, flux_up_v, flux_dn_v
    REAL(wp) :: hor_tend(nproma,p_patch%nlev+1,p_patch%nblks_e)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: vt_e(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jkp1, jb, je, jkm1, jc, jcn, jbn, jvn
    INTEGER  :: nlev              

    ! number of vertical levels
    nlev     = p_patch%nlev
    i_nchdom = MAX(1,p_patch%n_childdom)
   
    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk
   
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    !Some initializations
 
    IF(p_test_run)THEN
!ICON_OMP_WORKSHARE
      tot_tend(:,:,:) = 0._wp
      hor_tend(:,:,:) = 0._wp
!ICON_OMP_END_WORKSHARE
    END IF

    CALL rbf_vec_interpol_edge(p_nh_prog%vn, p_patch, p_int, vt_e, opt_rlend=min_rledge_int-1)

    ! 1) First get the horizontal tendencies at the half level edges

    rl_start = grf_bdywidth_e
    rl_end   = min_rledge_int-1

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jkm1,jcn,jbn,dvn1,dvn2,flux_up_c,flux_dn_c,&
!ICON_OMP            jvn,dvt1,dvt2,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 2 , nlev
       DO je = i_startidx, i_endidx
        
         jkm1    = jk-1

         !tendency in normal direction
         !flux = visc_c * D_31_c where D_31(=D_13) is calculated at half level 
         !cell center
 
         jcn     = iecidx(je,jb,2)  
         jbn     = iecblk(je,jb,2)  

         dvn2  = p_nh_diag%u(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                 p_nh_diag%v(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v2 - &
                 p_nh_diag%u(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v1   - &
                 p_nh_diag%v(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v2

         flux_up_c = (visc_smag_c(jcn,jk,jbn)*p_nh_metrics%wgtfac_c(jcn,jk,jbn) +  &
                     visc_smag_c(jcn,jkm1,jbn)*(1._wp-p_nh_metrics%wgtfac_c(jcn,jk,jbn))) * &
                     dvn2*p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) + &
                     (w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))-w_ie(je,jk,jb)) * &
                     p_patch%edges%inv_vert_vert_length(je,jb)*2.0_wp

         jcn     = iecidx(je,jb,1)  
         jbn     = iecblk(je,jb,1)  

         dvn1  = p_nh_diag%u(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                 p_nh_diag%v(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v2 - &
                 p_nh_diag%u(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v1   - &
                 p_nh_diag%v(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v2
           

         flux_dn_c = (visc_smag_c(jcn,jk,jbn)*p_nh_metrics%wgtfac_c(jcn,jk,jbn) +  &
                     visc_smag_c(jcn,jkm1,jbn)*(1._wp-p_nh_metrics%wgtfac_c(jcn,jk,jbn))) * &
                     dvn1*p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) + &
                     (w_ie(je,jk,jb)-w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))) *  &
                     p_patch%edges%inv_vert_vert_length(je,jb)*2.0_wp


         !tendency in tangential direction
         !flux = visc_v * D_32_v where D_32(=D_23) is calculated at half level 
         !between vertex and edge center

         jvn     = ividx(je,jb,2)
         jbn     = ivblk(je,jb,2)

         dvt2    = ( u_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                     v_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v2 + &
                     vt_e(je,jkm1,jb) ) * 0.5_wp           - &
                   ( u_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                     v_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v2 + &
                     vt_e(je,jk,jb) ) * 0.5_wp

         flux_up_v = (visc_smag_v(jvn,jk,jbn)*p_nh_metrics%wgtfac_v(jvn,jk,jbn) + &
                     visc_smag_v(jvn,jkm1,jbn)*(1._wp-p_nh_metrics%wgtfac_v(jvn,jk,jbn))) * &
                     dvt2*p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) + &
                     p_patch%edges%system_orientation(je,jb)*(w_vert(jvn,jk,jbn)-w_ie(je,jk,jb)) / &
                     p_patch%edges%edge_cell_length(je,jb,2)


         jvn     = ividx(je,jb,1)
         jbn     = ivblk(je,jb,1)

         dvt1    = ( u_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                     v_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v2 + &
                     vt_e(je,jkm1,jb) ) * 0.5_wp           - &
                   ( u_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                     v_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v2 + &
                     vt_e(je,jk,jb) ) * 0.5_wp

         flux_dn_v = (visc_smag_v(jvn,jk,jbn)*p_nh_metrics%wgtfac_v(jvn,jk,jbn) + &
                     visc_smag_v(jvn,jkm1,jbn)*(1._wp-p_nh_metrics%wgtfac_v(jvn,jk,jbn))) * &
                     dvt1*p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) + &
                     p_patch%edges%system_orientation(je,jb)*(w_ie(je,jk,jb)-w_vert(jvn,jk,jbn)) / &
                     p_patch%edges%edge_cell_length(je,jb,1)

         hor_tend(je,jk,jb) = ( (flux_up_c - flux_dn_c) * p_patch%edges%inv_dual_edge_length(je,jb) + &
                                 p_patch%edges%system_orientation(je,jb) * (flux_up_v - flux_dn_v)  * &
                                 p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp ) / rho_e(je,jk,jb)

       END DO
      END DO
    END DO                         
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    !CALL sync_patch_array(SYNC_E, p_patch, hor_tend)

    !Interpolate horizontal tendencies to w point: except top and bottom boundaries
    !w==0 at these boundaries

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
      
           tot_tend(jc,jk,jb) =                                                       &
              (hor_tend(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               hor_tend(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               hor_tend(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))           
      
         END DO
       END DO
    END DO
!ICON_OMP_END_DO

   ! 
   ! 2) Vertical tendency: evaluated at w point
   !

!ICON_OMP_DO PRIVATE(jb,jk,je,i_startidx,i_endidx,flux_up_c,flux_dn_c, &
!ICON_OMP            jkm1)
    DO jb = i_startblk,i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 2 , nlev
       DO jc = i_startidx, i_endidx

         jkm1 = jk - 1

         flux_up_c = visc_smag_c(jc,jkm1,jb) * ( (p_nh_prog%w(jc,jkm1,jb)-p_nh_prog%w(jc,jk,jb))* &
                     p_nh_metrics%inv_ddqz_z_full(jc,jkm1,jb) - 0.333333_wp*DIV_c(jc,jkm1,jb) )

         flux_dn_c = visc_smag_c(jc,jk,jb) * ( (p_nh_prog%w(jc,jk,jb)-p_nh_prog%w(jc,jk+1,jb)) * &
                     p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) - 0.333333_wp*DIV_c(jc,jk,jb) )

         tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + 2._wp * (flux_up_c - flux_dn_c) *  &
                              p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) / p_nh_prog%rho(jc,jk,jb)

         p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) + dt * tot_tend(jc,jk,jb)

       END DO
      END DO
    END DO  
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

   CALL sync_patch_array(SYNC_C, p_patch, p_nh_prog%w)

  END SUBROUTINE diffuse_vert_velocity
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_scalar
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based scalars
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - The Euler scheme is generally known to be unconditionaly unstable for diffusion term
  !!   but it is not certain when implemented in time split manner
  !! - This needs to be investigated and if found unstable we will switch to MacKormack 
  !!   scheme which is known to be stable
  !!  
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_scalar(var, p_nh_metrics, p_patch, p_int, p_nh_diag, tot_tend,  &
                            exner, prm_diag, rho, dt, scalar_name)

    REAL(wp),          INTENT(in)        :: var(:,:,:)   !input scalar
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    REAL(wp),        INTENT(inout),TARGET:: tot_tend(:,:,:)!<total tendency
    REAL(wp),          INTENT(in)        :: exner(:,:,:)   !
    REAL(wp),          INTENT(in)        :: rho(:,:,:)     !density at cell center
    TYPE(t_nwp_phy_diag),  INTENT(in)    :: prm_diag     !< atm phys vars
    REAL(wp),          INTENT(in)        :: dt
    CHARACTER(*), INTENT(in)             :: scalar_name    

    !Local variables
    INTEGER, PARAMETER :: iexplicit = 1
    INTEGER, PARAMETER :: iimplicit = 2
    INTEGER, PARAMETER :: vert_scheme_type = iexplicit

    REAL(wp) :: flux_up, flux_dn, inv_dt
    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: fac(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: sflux(nproma,p_patch%nblks_c)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: a, b, c, rhs
    REAL(wp), ALLOCATABLE, DIMENSION(:)     :: var_new
    REAL(wp), POINTER :: diff_smag_ic(:,:,:)
     
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc
    INTEGER  :: nlev              

    ! number of vertical levels
    nlev = p_patch%nlev
    i_nchdom   = MAX(1,p_patch%n_childdom)
    inv_dt  = 1._wp / dt

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk
   
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    diff_smag_ic => prm_diag%tkvh

    !Special treatment for different scalars: includes
    !boundary treatment also

!ICON_OMP_PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
     rl_start = grf_bdywidth_c+1
     rl_end   = min_rlcell_int

     i_startblk = p_patch%cells%start_blk(rl_start,1)
     i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF(TRIM(scalar_name)=='theta')THEN
      !multiply by exner to convert from theta tend to temp tend
      !assuming that exner perturbation are small compared to temp
!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              fac(jc,jk,jb) = cpd * rcvd * exner(jc,jk,jb) / rho(jc,jk,jb)
            END DO
          END DO
          DO jc = i_startidx, i_endidx
             sflux(jc,jb) = prm_diag%shfl_s(jc,jb) * rcpd
          END DO
        END DO
!ICON_OMP_END_DO
    ELSEIF(TRIM(scalar_name)=='qv')THEN
!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              fac(jc,jk,jb) = 1._wp / rho(jc,jk,jb)
            END DO
          END DO
          DO jc = i_startidx, i_endidx
             sflux(jc,jb) = prm_diag%lhfl_s(jc,jb) / alv
          END DO
        END DO
!ICON_OMP_END_DO
    ELSEIF(TRIM(scalar_name)=='qc')THEN
!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              fac(jc,jk,jb) = 1._wp / rho(jc,jk,jb)
            END DO
          END DO
          DO jc = i_startidx, i_endidx
             sflux(jc,jb) = 0._wp
          END DO
        END DO
!ICON_OMP_END_DO
    END IF
      
    ! use conservative discretization div(k*grad(var))-horizontal part  
    ! following mo_nh_diffusion

    !include halo points and boundary points because these values will be 
    !used in next loop
    rl_start = grf_bdywidth_e
    rl_end   = min_rledge_int-1 

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!ICON_OMP_DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! compute kh_smag_e * grad(var) 
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
                nabla2_e(je,jk,jb) = diff_smag_e(je,jk,jb) *          &
                         p_patch%edges%inv_dual_edge_length(je,jb)*   &
                        (var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                         var(iecidx(je,jb,1),jk,iecblk(je,jb,1)))
            ENDDO
          ENDDO
        ENDDO
!ICON_OMP_END_DO

        ! now compute the divergence of the quantity above
        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              !horizontal tendency
              tot_tend(jc,jk,jb)  =  fac(jc,jk,jb) * (                                 &
                nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb) )
            ENDDO
          ENDDO
        ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

       !---------------------------------------------------------------
       !Vertical diffusion
       !---------------------------------------------------------------

       rl_start = grf_bdywidth_c+1
       rl_end   = min_rlcell_int

       i_startblk = p_patch%cells%start_blk(rl_start,1)
       i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


       SELECT CASE(vert_scheme_type)

       CASE(iexplicit)
       !Vertical tendency - explicit solver

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,flux_up,flux_dn)
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 2, nlev-1
            DO jc = i_startidx, i_endidx
              flux_up = diff_smag_ic(jc,jk,jb) * (var(jc,jk-1,jb) - var(jc,jk,jb)) * &
                        p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) 

              flux_dn = diff_smag_ic(jc,jk+1,jb) * (var(jc,jk,jb) - var(jc,jk+1,jb)) * &
                        p_nh_metrics%inv_ddqz_z_half(jc,jk+1,jb) 

              tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + (flux_up - flux_dn) *  &
                              p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * fac(jc,jk,jb)

            ENDDO
          ENDDO
        ENDDO
!ICON_OMP_END_DO

        !Boundary treatment
       
        !--------------------------------------------------------    
        !jk = 1        
        !--------------------------------------------------------    
!ICON_OMP_DO PRIVATE(jc,jb,i_startidx,i_endidx,flux_dn,flux_up)
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
            DO jc = i_startidx, i_endidx
              flux_up = 0._wp

              flux_dn = diff_smag_ic(jc,2,jb) * (var(jc,1,jb) - var(jc,2,jb)) * &
                        p_nh_metrics%inv_ddqz_z_half(jc,2,jb) 

              tot_tend(jc,1,jb) = tot_tend(jc,1,jb) + (flux_up- flux_dn) * &
                          p_nh_metrics%inv_ddqz_z_full(jc,1,jb) * fac(jc,1,jb)
            ENDDO
        ENDDO
!ICON_OMP_END_DO

        !--------------------------------------------------------    
        !jk = nlev
        !--------------------------------------------------------    

!ICON_OMP_DO PRIVATE(jc,jb,i_startidx,i_endidx,flux_up,flux_dn)
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
            DO jc = i_startidx, i_endidx
              flux_up = diff_smag_ic(jc,nlev,jb) * (var(jc,nlev-1,jb) - var(jc,nlev,jb)) * &
                        p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) 

              flux_dn  = -sflux(jc,jb)

              tot_tend(jc,nlev,jb) = tot_tend(jc,nlev,jb) + (flux_up - flux_dn) * & 
                      p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * fac(jc,nlev,jb)
            ENDDO
        ENDDO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
       
       CASE(iimplicit)

        !vertical tendency- implicit solver
        !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

       ! LL: The local arrays should be private for each block, then the allocation and the solver will be called inside the omp section
       ALLOCATE( a(nproma,nlev,p_patch%nblks_c), b(nproma,nlev,p_patch%nblks_c), var_new(nlev),  &
                 c(nproma,nlev,p_patch%nblks_c), rhs(nproma,nlev,p_patch%nblks_c) )

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
       DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 2,nlev-1
            DO jc = i_startidx, i_endidx
              a(jc,jk,jb)   = - diff_smag_ic(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * &
                              p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) 

              c(jc,jk,jb)   = - diff_smag_ic(jc,jk+1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * &
                              p_nh_metrics%inv_ddqz_z_half(jc,jk+1,jb)
           
              b(jc,jk,jb)   =  inv_dt - a(jc,jk,jb) - c(jc,jk,jb)

              rhs(jc,jk,jb) =   var(jc,jk,jb) * inv_dt
            END DO
          END DO
       END DO
!ICON_OMP_END_DO

      !TOP Boundary treatment
      ! jk = 1
      !
!ICON_OMP_DO PRIVATE(jc,jb,i_startidx,i_endidx)
       DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx
            a(jc,1,jb)   = 0._wp
            c(jc,1,jb)   = -diff_smag_ic(jc,2,jb) * p_nh_metrics%inv_ddqz_z_full(jc,1,jb) * &
                             p_nh_metrics%inv_ddqz_z_half(jc,2,jb) 
            b(jc,1,jb)   = inv_dt - a(jc,1,jb) - c(jc,1,jb)
            rhs(jc,1,jb) = var(jc,1,jb) * inv_dt 
          END DO
       END DO     
!ICON_OMP_END_DO

       !SFC Boundary treatment
       !jk = nlev
       !
!ICON_OMP_DO PRIVATE(jc,jb,i_startidx,i_endidx)
       DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx
            a(jc,nlev,jb)  = - diff_smag_ic(jc,nlev,jb) * p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) 
            c(jc,nlev,jb)  = 0._wp
            b(jc,nlev,jb)  = inv_dt - a(jc,nlev,jb) - c(jc,nlev,jb)
            rhs(jc,nlev,jb)= var(jc,nlev,jb) * inv_dt + sflux(jc,jb) * p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) 
          END DO
       END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

       !CALL TDMA

       DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jc = i_startidx, i_endidx
            CALL tdma_solver(a(jc,:,jb),b(jc,:,jb),c(jc,:,jb),rhs(jc,:,jb),nlev,var_new)
            tot_tend(jc,:,jb) = tot_tend(jc,:,jb) + (var_new(:)-var(jc,:,jb))*inv_dt*fac(jc,:,jb)
         END DO
       END DO

       DEALLOCATE( a, b, c, rhs, var_new )

    END SELECT !vert_scheme

  END SUBROUTINE diffuse_scalar

!---------------------------------------------------------------------------------

END MODULE mo_sgs_turbulence
