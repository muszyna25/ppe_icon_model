!>
!! This module contains subroutines for evaluating the right-hand side
!! of the primitive equations
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (MPI-M, 2009-11)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ha_dynamics_adv

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message,finish
  USE mo_model_domain,       ONLY: t_patch
  USE mo_math_gradients,     ONLY: grad_fd_norm
  USE mo_math_divrot,        ONLY: div, div_avg
  USE mo_dynamics_config,    ONLY: idiv_method, lshallow_water
  USE mo_io_config,          ONLY: l_diagtime
  USE mo_parallel_config,    ONLY: nproma
  USE mo_run_config,         ONLY: nlev, nlevm1, nlevp1
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_interpol_config,    ONLY: sick_a, sick_o
  USE mo_intp_rbf,           ONLY: rbf_vec_interpol_edge
  USE mo_intp,               ONLY: cells2edges_scalar, edges2cells_scalar, &
                                   verts2edges_scalar, cells2verts_scalar, &
                                   edges2verts_scalar, verts2cells_scalar, &
                                   edges2edges_scalar
  USE mo_interpol_config,    ONLY: i_cori_method, l_corner_vort
  USE mo_vertical_coord_table, ONLY: rdelpr, nplev, nplvp1
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,     ONLY: min_rledge
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_sync,               ONLY: SYNC_C, SYNC_V, sync_patch_array

  IMPLICIT NONE

  PUBLIC

CONTAINS

  !>
  !! Purpose: calculate the tendency of normal wind induced by vertical
  !! advection
  !!
  !! General comment:
  !! For nested domains, tendencies are interpolated from the parent domain
  !! on a boundary zone with a width of grf_bdywidth_c for cells and
  !! grf_bdywidth_e for edges, respectively. These tendencies must not be
  !! overwritten.
  !!
  !! Reference:
  !! Simmons and Burridge (1981, Mon. Wea. Rev.)
  !! Roeckner et al. (2005, MPI Report 349)
  !!
  !! @par Revision History
  !! Separated from subroutine dyn and re-written by Hui Wan (MPI-M, 2009-11-17)
  !!
  SUBROUTINE vn_adv_vertical( p_vn, p_weta_c, p_delp_e, &
                              pt_patch, pt_int_state,   &
                              p_ddt_vn )
  !! Arguments

  REAL(wp),INTENT(IN)    :: p_vn    (:,:,:)   !<  normal velocity
  REAL(wp),INTENT(INOUT) :: p_weta_c(:,:,:)   !<  cell-based vertical velocity
  REAL(wp),INTENT(IN)    :: p_delp_e(:,:,:)   !<  layer thickness at edges

  TYPE(t_patch),TARGET,INTENT(IN) :: pt_patch   !<  grid information
  TYPE(t_int_state),INTENT(IN) :: pt_int_state  !<  interpolation coefficients

  REAL(wp), INTENT(INOUT) :: p_ddt_vn(:,:,:)  !<  tendency of normal velocity

  !! Local variables

  INTEGER  :: nblks_e
  INTEGER  :: jb,jbs,is,ie,jk,jkp

  REAL(wp) :: z_weta_e( SIZE(p_vn,1),SIZE(p_weta_c,2),SIZE(p_vn,3) )
  REAL(wp) :: z_tmp_e ( SIZE(p_vn,1),SIZE(p_vn,    2),SIZE(p_vn,3) )

! Dimension parameter

   nblks_e  = pt_patch%nblks_e

! Interpolate vertical velocity (rho*eta-dot) from cell centers to edges

  ! LL: The vertical velocity is computed on the halos from the
  !     continuity
  !  CALL sync_patch_array( SYNC_C, pt_patch, p_weta_c )

   CALL cells2edges_scalar( p_weta_c, pt_patch, pt_int_state%c_lin_e, &! in
                            z_weta_e, opt_rlstart=4 )    ! out, optional in

! Tendency of velocity due to vertical advection.

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)

      z_tmp_e(is:ie,1,jb) = 0._wp

      DO jk = 1, nlevm1
         jkp = jk + 1
         z_tmp_e(is:ie,jkp,jb) = z_weta_e(is:ie,jkp,jb)  &
                                 *( p_vn(is:ie,jk,jb)- p_vn(is:ie,jkp,jb) )
         z_tmp_e(is:ie,jk ,jb) = z_tmp_e(is:ie,jk,jb) + z_tmp_e(is:ie,jkp,jb)
      ENDDO

     !Pressure levels
      DO jk = 1,nplev
        p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb) &
                              + z_tmp_e(is:ie,jk,jb)*rdelpr(jk)*0.5_wp
      ENDDO

     !Sigma and transition levels
      DO jk = nplvp1,nlev
        p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb) &
                              + z_tmp_e(is:ie,jk,jb)/p_delp_e(is:ie,jk,jb) &
                               *0.5_wp
      ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE vn_adv_vertical

  !>
  !! Vertical advection of temperature
  !!
  !! For nested domains, tendencies are interpolated from the parent domain
  !! on a boundary zone with a width of grf_bdywidth_c for cells and
  !! grf_bdywidth_e for edges, respectively. These tendencies must not be
  !! overwritten.
  !!
  !! Reference:
  !! Simmons and Burridge (1981, Mon. Wea. Rev.)
  !! Roeckner et al. (2005, MPI Report 349)
  !!
  !! @par Revision History
  !! Separation from subroutine dyn and rewriting by Hui Wan (MPI-M, 2009-11-18)
  !!
  SUBROUTINE temp_adv_vertical( p_temp, p_weta, p_rdelp, pt_patch, p_ddt_temp )

  !! Arguments

  REAL(wp),INTENT(IN)    :: p_temp    (:,:,:) !< temperature
  REAL(wp),INTENT(IN)    :: p_weta    (:,:,:) !< vertical velocity rho*eta-dot
  REAL(wp),INTENT(IN)    :: p_rdelp   (:,:,:) !< 1./layer_thickness (delp)
  REAL(wp),INTENT(INOUT) :: p_ddt_temp(:,:,:) !< temperature tendency
  TYPE(t_patch),INTENT(IN) :: pt_patch          !< grid information

  !! Local variables

  INTEGER  :: jb,jbs,is,ie,jk,jkp
  INTEGER  :: nblks_c
  REAL(wp) :: z_tmp_c (SIZE(p_temp,1),SIZE(p_temp,2),SIZE(p_temp,3))

! Dimension parameter

  nblks_c = pt_patch%nblks_c

! Vertical advection

!$OMP PARALLEL  PRIVATE(jbs)
  jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = jbs,nblks_c
     CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)

     z_tmp_c(is:ie,1,jb) = 0._wp

     DO jk = 1, nlevm1
        jkp = jk + 1
        z_tmp_c(is:ie,jkp,jb) = p_weta(is:ie,jkp,jb)  &
                               *( p_temp(is:ie,jk, jb) - p_temp(is:ie,jkp,jb) )
        z_tmp_c(is:ie,jk ,jb) = z_tmp_c(is:ie,jk,jb) + z_tmp_c(is:ie,jkp,jb)
     ENDDO

     p_ddt_temp(is:ie,1:nlev,jb) = p_ddt_temp(is:ie,1:nlev,jb)  &
                                 + z_tmp_c(is:ie,1:nlev,jb)     &
                                  *p_rdelp(is:ie,1:nlev,jb)*0.5_wp
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE temp_adv_vertical

  !>
  !! Horizontal advection of temperature
  !!
  !! Energy conserving formula [ div(T*delp*V) - T*div(delp*V) ]/delp.
  !! For nested domains, tendencies are interpolated from the parent domain
  !! on a boundary zone with a width of grf_bdywidth_c for cells and
  !! grf_bdywidth_e for edges, respectively. These tendencies must not be
  !! overwritten.
  !!
  !! @par Revision History
  !! Separation from subroutine dyn and rewriting by Hui Wan (MPI-M, 2009-11-18)
  !!
  SUBROUTINE temp_adv_horizontal( p_temp, p_mflux_e, p_rdelp_c, p_mdiv, &
                                  pt_patch, pt_int_state, p_ddt_temp )
  !! Arguments

  REAL(wp),INTENT(IN)    :: p_temp     (:,:,:)
  REAL(wp),INTENT(IN)    :: p_mflux_e  (:,:,:)
  REAL(wp),INTENT(IN)    :: p_rdelp_c  (:,:,:)
  REAL(wp),INTENT(IN)    :: p_mdiv     (:,:,:)
  REAL(wp),INTENT(INOUT) :: p_ddt_temp (:,:,:)

  TYPE(t_patch),TARGET,INTENT(IN) :: pt_patch
  TYPE(t_int_state),INTENT(IN) :: pt_int_state

  !! Local variables

  INTEGER :: jb,jbs,is,ie,nblks_e,nblks_c

  REAL(wp),DIMENSION(SIZE(p_mflux_e,1),SIZE(p_mflux_e,2),SIZE(p_mflux_e,3)) &
          :: z_temp_e, z_flux_e, z_tmdiv

! Dimension parameters
  nblks_e = pt_patch%nblks_e
  nblks_c = pt_patch%nblks_c

! Interpolate temperature from cells to edges
  CALL cells2edges_scalar( p_temp, pt_patch, pt_int_state%c_lin_e, z_temp_e )

! Calculate heat divergence at cell centres

   jbs = pt_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 2)
      z_flux_e(is:ie,:,jb) = z_temp_e(is:ie,:,jb)*p_mflux_e(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  SELECT CASE(idiv_method)
  CASE(1)
    CALL div( z_flux_e, pt_patch, pt_int_state, z_tmdiv, opt_rlstart=2 )
  CASE(2)
    CALL div_avg( z_flux_e, pt_patch, pt_int_state, pt_int_state%c_bln_avg, &
                  z_tmdiv, opt_rlstart=2 )
  END SELECT

! Temperature tendency due to horizontal advection

   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)

      p_ddt_temp(is:ie,:,jb) = p_ddt_temp(is:ie,:,jb)     &
                             - p_rdelp_c(is:ie,:,jb)      &
                              *( z_tmdiv(is:ie,:,jb)      &
                                -p_mdiv(is:ie,:,jb)*p_temp(is:ie,:,jb) )

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE temp_adv_horizontal


  !>
  !! Horizontal momentum advection and Coriolis force
  !! The algorithm on triangular grids is different from that on the
  !! hexagonal/pentagonal grids.
  !!
  !! @par Revision History
  !! Separated from mo_nonlinear_advection and re-written
  !! by Hui Wan (MPI-M, 2009-11-20)
  !!
  SUBROUTINE vn_adv_horizontal( p_vn, p_vort, p_delp_c,  &
                                pt_patch, pt_int,        &
                                p_ddt_vn, p_kin, p_vt,   &
                                p_delp_v,                &
                                opt_rlstart, opt_rlend  )
  !! Arguments
  IMPLICIT NONE

  REAL(wp),INTENT(IN) :: p_vn     (:,:,:)  !< normal velocity at edges
  REAL(wp),INTENT(IN) :: p_vort   (:,:,:)  !< relative vorticity at dual centers
  REAL(wp),INTENT(IN) :: p_delp_c (:,:,:)  !< pseudo-dencity at cell centers

  TYPE(t_patch),    TARGET,INTENT(IN) :: pt_patch
  TYPE(t_int_state),TARGET,INTENT(IN) :: pt_int

  REAL(wp),INTENT(INOUT) :: p_ddt_vn (:,:,:) !< tendency of normal wind
  REAL(wp),INTENT(INOUT) :: p_kin    (:,:,:) !< kinetic energy at cell centers
  REAL(wp),INTENT(INOUT) :: p_vt     (:,:,:) !< tangential wind at edges
  REAL(wp),INTENT(INOUT) :: p_delp_v (:,:,:) !< pseudo-dencity at dual grid

  INTEGER,INTENT(IN),OPTIONAL :: opt_rlstart, opt_rlend
                                 !< start and end values of refin_ctrl flag
  !! Local variables

  INTEGER  :: jb,jbs,jbe,is,ie,jk
  INTEGER  :: nblks_e, npromz_e, i_nchdom, nblks_c, npromz_c
  INTEGER  :: rl_start, rl_end

  REAL(wp) :: z_tmp_e    ( nproma, nlev, pt_patch%nblks_e )
  REAL(wp) :: z_vort_e   ( nproma, nlev, pt_patch%nblks_e )

! Dimension parameters

  nblks_e  = pt_patch%nblks_e
  npromz_e = pt_patch%npromz_e
  nblks_c  = pt_patch%nblks_c
  npromz_c = pt_patch%npromz_c

  i_nchdom = MAX(1,pt_patch%n_childdom)

! For diagnosing potential enstrophy in the
! shallow water model, interpolate the pseudo-dencity from triangle centers
! to vertices.
  IF (l_diagtime .AND. lshallow_water) THEN
    CALL cells2verts_scalar(p_delp_c, pt_patch, &
                            pt_int%cells_aw_verts, p_delp_v)
  ENDIF

!-----------------------------------------
! Kinetic energy and tangential velocity
!-----------------------------------------

! Reconstruct the tangential wind at edges from the normal velocity;
! Define the specific kinetic energy at edges as ke = 0.5*( vn*vn + vt*vt )
! The cell-based value needed for computing ke gradient is then obtained
! by interpolation.

  CALL rbf_vec_interpol_edge( p_vn, pt_patch, pt_int, p_vt)

  jbs = pt_patch%edges%start_blk(2,1) !for rbf_vec_dim_edge==4
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1,nblks_e
     CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 2)
     z_tmp_e(is:ie,:,jb) = 0.5_wp*( p_vt(is:ie,:,jb)*p_vt(is:ie,:,jb) &
                                   +p_vn(is:ie,:,jb)*p_vn(is:ie,:,jb) )
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CALL edges2cells_scalar( z_tmp_e, pt_patch, pt_int%e_bln_c_s,  &
                           p_kin, opt_rlstart=2 )

!--------------------------
! Absolute vorticity flux
!--------------------------
! Dimension parameters

  IF ( PRESENT(opt_rlstart) ) THEN
    IF ((opt_rlstart >= 1) .AND. (opt_rlstart <= 2)) THEN
      CALL finish ('mo_ha_dynamics_adv:vn_adv_horizontal',  &
                   'opt_rlstart must not be 1 or 2')
    ENDIF
    rl_start = opt_rlstart
  ELSE
    rl_start = 3
  END IF

  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rledge
  END IF

! Compute the absolute vorticity flux

  CALL verts2edges_scalar( p_vort, pt_patch, pt_int%v_1o2_e, &
                           z_vort_e, opt_rlstart=3)

  jbs = pt_patch%edges%start_blk(rl_start,1)
  jbe = pt_patch%edges%end_blk(rl_end,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = jbs,jbe
     CALL get_indices_e(pt_patch, jb,jbs,jbe, is,ie, rl_start,rl_end)

     DO jk = 1,nlev
        p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb)             &
                              - p_vt(is:ie,jk,jb)                 &
                               *( z_vort_e(is:ie,jk,jb)           &
                                 +pt_patch%edges%f_e(is:ie,jb) )
     ENDDO
     !!HW: minus sign here because the tangential direction in the code is
     !!    in fact the opposite from what is written in BR05.
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!=============================================================================
! Calculate the gradient of kinetic energy, and accumulate velocity tendency
!=============================================================================

  CALL grad_fd_norm( p_kin, pt_patch, z_tmp_e, opt_rlstart=4 )

  jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      p_ddt_vn(is:ie,:,jb) = p_ddt_vn(is:ie,:,jb) - z_tmp_e(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   END SUBROUTINE vn_adv_horizontal

END MODULE mo_ha_dynamics_adv
