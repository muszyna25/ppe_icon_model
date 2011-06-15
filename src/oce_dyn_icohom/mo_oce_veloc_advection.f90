!>
!! Contains the implementation of velocity advection in vector invariant form
!! that is used in the ocean model.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2010)
!!  Modified by Stephan Lorenz,     MPI-M (2010-11)
!!   - implementation of new PtP reconstruction
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_oce_veloc_advection
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_run_nml,             ONLY: nproma
USE mo_math_constants
USE mo_physical_constants
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, &
  &                               sea_boundary, sea, boundary, &
  &                               full_coriolis, beta_plane_coriolis,&
  &                               f_plane_coriolis, zero_coriolis
USE mo_model_domain,        ONLY: t_patch
USE mo_ocean_nml,           ONLY: n_zlev,iswm_oce, ab_beta, ab_gam, L_INVERSE_FLIP_FLOP, &
  &                               CORIOLIS_TYPE, basin_center_lat!, basin_center_lon,  &
!  &                               basin_width_deg, basin_height_deg  
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_oce_state,           ONLY: t_hydro_ocean_diag, t_hydro_ocean_aux
USE mo_oce_math_operators,  ONLY: rot_vertex_ocean,rot_vertex_ocean_origin,rot_vertex_ocean_total,&
 &                                grad_fd_norm_oce
USE mo_math_utilities,      ONLY: gvec2cvec, t_cartesian_coordinates!gc2cc,cc2gc

USE mo_scalar_product,      ONLY: map_cell2edges, dual_flip_flop, &
  &                               primal_map_c2e, map_edges2edges
USE mo_interpolation,       ONLY: t_int_state, rbf_vec_interpol_edge,&
  &                               edges2cells_scalar, verts2edges_scalar,&
  &                               rbf_vec_interpol_cell
!USE mo_oce_index,              ONLY: c_k, ne_b, ne_i, form4ar! , ldbg
IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: veloc_adv_horz_mimetic
PUBLIC :: veloc_adv_vert_mimetic
PUBLIC :: veloc_adv_horz_RBF
PUBLIC :: veloc_adv_vert_RBF
PUBLIC :: veloc_adv_diff_vert_RBF
!PUBLIC :: veloc_adv_diff_vert_RBF2

CONTAINS
!-------------------------------------------------------------------------
!
!
!>
!! Computes horizontal advection of a (edge based) vector field.
!!
!! Computes rotational term of a vector field given by its components in
!! the directions normal to triangle edges and the gradient of the kinetic energy
!! which is calculated using the reconstructed velocity at cell centers. Both
!! terms are combined and constitute the horizontal velocity advection.
!!
!!IMPORTANT: It is assumed that the reconstruction of the tangential velocity
!!           has been done before.
!1
!! input:  lives on edges (velocity points)
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE veloc_adv_horz_mimetic( p_patch, vn_old,vn_new, p_diag, veloc_adv_horz_e, p_int)
!
!
!  patch on which computation is performed
!
TYPE(t_patch), INTENT(in) :: p_patch

!
! normal and tangential velocity  of which advection is computed
!
REAL(wp), INTENT(inout) :: vn_old(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
REAL(wp), INTENT(inout) :: vn_new(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
!
!diagnostic ocean state stores horizontally advected velocity
!
TYPE(t_hydro_ocean_diag) :: p_diag
!
! variable in which horizontally advected velocity is stored
!
REAL(wp), INTENT(out) :: veloc_adv_horz_e(:,:,:)
!
!Interpolation necessary just for testing
TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int


INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jk, jb, jc, je!, jv, ile, ibe, ie, jev
INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: i_startblk_v, i_endblk_v!, i_startidx_v, i_endidx_v
INTEGER :: rl_start_e, rl_end_e, rl_start_v, rl_end_v, rl_start_c, rl_end_c
!INTEGER ::  i_v1_idx, i_v1_blk, i_v2_idx, i_v2_blk
REAL(wp) :: z_e  (nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_u_c(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_v_c(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_vort_glb(nproma,n_zlev,p_patch%nblks_v)
REAL(wp) :: z_vort_flx(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_grad_ekin_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_RBF_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_RBF_c(nproma,n_zlev,p_patch%nblks_c)
!REAL(wp) :: z_tmp
!REAL(wp) :: z_beta_plane_vort
!REAL(wp) :: z_f_plane_vort
!REAL(wp) :: z_lat_basin_center
!REAL(wp) :: z_y
!INTEGER :: i_v1_ctr, i_v2_ctr
!INTEGER :: i_ctr
!INTEGER :: il_c1, ib_c1, il_c2, ib_c2
!INTEGER :: il_e1, ib_e1,il_e2, ib_e2, il_e3, ib_e3
INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3
REAL(wp) :: z_weight_e1, z_weight_e2, z_weight_e3!, z_weight
!ARRAYS FOR TESTING
REAL(wp) :: z_vt(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_flx_RBF(nproma,n_zlev,p_patch%nblks_e)
!REAL(wp) :: z_veloc_adv_horz_e(nproma,n_zlev,p_patch%nblks_e)
LOGICAL, PARAMETER :: L_DEBUG = .TRUE. 
!-----------------------------------------------------------------------
! #slo# set local variable to zero due to nag -nan compiler-option
z_e             (:,:,:) = 0.0_wp
z_u_c           (:,:,:) = 0.0_wp
z_v_c           (:,:,:) = 0.0_wp
z_vort_glb      (:,:,:) = 0.0_wp
z_vort_flx      (:,:,:) = 0.0_wp
veloc_adv_horz_e(:,:,:) = 0.0_wp
z_vt            (:,:,:) = 0.0_wp
z_vort_e        (:,:,:) = 0.0_wp
z_vort_flx_RBF  (:,:,:) = 0.0_wp
z_kin_RBF_c     (:,:,:) = 0.0_wp
z_grad_ekin_RBF (:,:,:) = 0.0_wp

rl_start_v = 1
rl_end_v   = min_rlvert
rl_start_e = 1
rl_end_e   = min_rledge
rl_start_c = 1
rl_end_c   = min_rlcell

i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)

slev = 1
elev = n_zlev
!z_e = ab_gam*vn_new + (1.0_wp-ab_gam)*vn_old

! calculate global vorticity for all vertical layers
!CALL rot_vertex_ocean_total(p_diag%ptp_vn, p_patch, p_diag%vort)
CALL rot_vertex_ocean_total(vn_old, p_patch, p_diag%vort)

!calculate vorticity flux across dual edge
IF ( iswm_oce == 1 ) THEN
  z_vort_flx(:,:,:) = dual_flip_flop(p_patch, vn_old, vn_new, p_diag%vort,p_diag%thick_e)!p_diag%h_e
ELSEIF ( iswm_oce /= 1 ) THEN
  z_vort_flx(:,:,:) = dual_flip_flop(p_patch, vn_old, vn_new, p_diag%vort, p_diag%h_e)
ENDIF
 DO jk = slev, elev
 write(*,*)'max/min vorticity:              ', jk,MAXVAL(p_diag%vort(:,jk,:)),&
                                                 &MINVAL(p_diag%vort(:,jk,:)) 
 write(*,*)'max/min vort flux:              ', jk,MAXVAL(z_vort_flx(:,jk,:)),&
                                                 &MINVAL(z_vort_flx(:,jk,:)) 
 END DO

!calculate gradient of kinetic energy
!The kinetic energy is already calculated at beginning
!of actual timestep in sbr "calc_scalar_product_for_veloc"
CALL grad_fd_norm_oce( p_diag%kin, &
                 & p_patch,    &
                 & p_diag%grad,&
                 & opt_slev=slev,opt_elev=elev )
  write(*,*)'max/min kin energy:            ', MAXVAL(p_diag%kin(:,1,:)),&
                                              &MINVAL(p_diag%kin(:,1,:))
  write(*,*)'max/min grad kin energy:       ', MAXVAL(p_diag%grad(:,1,:)),&
                                              &MINVAL(p_diag%grad(:,1,:)) 

  write(987,*)'max/min kin energy:          ', MAXVAL(p_diag%kin(:,1,:)),&
                                              &MINVAL(p_diag%kin(:,1,:))
  write(987,*)'max/min grad kin energy:     ', MAXVAL(p_diag%grad(:,1,:)),&
                                              &MINVAL(p_diag%grad(:,1,:)) 

IF(L_INVERSE_FLIP_FLOP)THEN
  CALL map_edges2edges( p_patch,    &
                      & p_diag%grad,&
                      & z_e,        &
                      & p_diag%thick_e)
  p_diag%grad = z_e
ENDIF

!---------------------------for testing and comparison with RBF--------------------------
!----------nonlinear coriolis and grad of kinetic energy computed with RBFs--------------
IF (L_DEBUG) THEN

CALL rbf_vec_interpol_edge( vn_old,       &
                          & p_patch,  &
                          & p_int,    &
                          & p_diag%vt,&
                          & opt_slev=slev,opt_elev=elev)
DO jb = i_startblk_e,i_endblk_e
  CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
    &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
  DO jk = slev, elev
    DO je=i_startidx_e, i_endidx_e
      IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) == boundary ) THEN
        p_diag%vt(je,jk,jb) = 0.0_wp
        vn_old(je,jk,jb) = 0.0_wp
      ENDIF
    END DO
  END DO
END DO

CALL rot_vertex_ocean(vn_old, p_diag%vt, p_patch, p_diag%vort)
CALL verts2edges_scalar( p_diag%vort, p_patch, p_int%v_1o2_e, &
                         z_vort_e, opt_slev=slev,opt_elev=elev, opt_rlstart=3)

  DO jb = i_startblk_e,i_endblk_e
    CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
      &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
    DO jk = slev, elev
      DO je=i_startidx_e, i_endidx_e
        IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) == sea ) THEN
         z_vort_flx_RBF(je,jk,jb) =&
          & p_diag%vt(je,jk,jb)*(z_vort_e(je,jk,jb)+ p_patch%edges%f_e(je,jb))         !& p_diag%vt(je,jk,jb)*p_patch%edges%f_e(je,jb)
        ELSE
          z_vort_flx_RBF(je,jk,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
  ENDDO

CALL rbf_vec_interpol_cell( vn_old, p_patch, p_int, p_diag%u,  &
  &                         p_diag%v, opt_slev=slev, opt_elev=elev)
  !write(*,*)'max/min vort flux:', MAXVAL(z_vort_flx_RBF(:,1,:)),MINVAL(z_vort_flx_RBF(:,1,:)) 
DO jb = i_startblk_e,i_endblk_e
  CALL get_indices_e(p_patch, jb, i_startblk_e,i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e
      ! calculate kinetic energy at edges from normal and tangential comp.
      z_kin_RBF_e(je,jk,jb) =0.5_wp*(p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb)&
                            &          +    vn_old(je,jk,jb)*vn_old(je,jk,jb) )
    ENDDO
  ENDDO
ENDDO

DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
                     i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
  DO jk = slev, elev
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) > sea_boundary ) THEN
          z_kin_RBF_c(jc,jk,jb) = 0.0_wp
      ELSE

        ile1 = p_patch%cells%edge_idx(jc,jb,1)
        ibe1 = p_patch%cells%edge_blk(jc,jb,1)
        ile2 = p_patch%cells%edge_idx(jc,jb,2)
        ibe2 = p_patch%cells%edge_blk(jc,jb,2)
        ile3 = p_patch%cells%edge_idx(jc,jb,3)
        ibe3 = p_patch%cells%edge_blk(jc,jb,3)
        z_weight_e1 = 0.0_wp
        z_weight_e2 = 0.0_wp
        z_weight_e3 = 0.0_wp
        IF(p_patch%patch_oce%lsm_oce_e(ile1,jk,ibe1)<= boundary)THEN
          z_weight_e1 = p_patch%edges%area_edge(ile1,ibe1)
        ENDIF
        IF(p_patch%patch_oce%lsm_oce_e(ile2,jk,ibe2)<= boundary)THEN
          z_weight_e2 = p_patch%edges%area_edge(ile2,ibe2)
        ENDIF
        IF(p_patch%patch_oce%lsm_oce_e(ile3,jk,ibe3)<= boundary)THEN
          z_weight_e3 = p_patch%edges%area_edge(ile3,ibe3)
        ENDIF

        !write(*,*)'weights',jc,jk,jb,z_weight_e1,z_weight_e2,z_weight_e3
         z_kin_RBF_c(jc,jk,jb) = (z_kin_RBF_e(ile1,jk,ibe1)*z_weight_e1&
                             &+ z_kin_RBF_e(ile2,jk,ibe2)*z_weight_e2&
                             &+ z_kin_RBF_e(ile3,jk,ibe3)*z_weight_e3)&
                             &/(z_weight_e1+z_weight_e2+z_weight_e3) 
      ENDIF
    END DO
  END DO
END DO

CALL grad_fd_norm_oce( p_diag%kin, &
                     & p_patch,    &
                     & z_grad_ekin_RBF, opt_slev=slev,opt_elev=elev)

END IF
!--------------END OF TESTING----------------------------------------------------------

!Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
DO jb = i_startblk_e, i_endblk_e
  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e

      !z_veloc_adv_horz_e(je,jk,jb)  = z_vort_flx_RBF(je,jk,jb) + z_grad_ekin_RBF(je,jk,jb) 
      !veloc_adv_horz_e(je,jk,jb)= z_vort_flx(je,jk,jb) + z_grad_ekin_RBF(je,jk,jb)
      !veloc_adv_horz_e(je,jk,jb)    = z_vort_flx_RBF(je,jk,jb)+ p_diag%grad(je,jk,jb)
      veloc_adv_horz_e(je,jk,jb)    = z_vort_flx(je,jk,jb) + p_diag%grad(je,jk,jb)
!        write(*,*)'horz adv:vort-flx:vort-flx-RBF',je,jk,jb,z_vort_flx(je,1,jb), &
!          &        z_vort_flx_RBF(je,jk,jb),p_diag%grad(je,jk,jb),z_grad_ekin_RBF(je,jk,jb),&
!          &        veloc_adv_horz_e(je,jk,jb), z_veloc_adv_horz_e(je,jk,jb)
!       IF( (z_vort_flx(je,jk,jb)>0.0.AND.z_vort_flx_RBF(je,jk,jb)<0)&
!       &   .OR. (z_vort_flx(je,jk,jb)<0.0.AND.z_vort_flx_RBF(je,jk,jb)>0))THEN
!             Write(*,*)'wraning:differing signs',je,jk,jb,z_vort_flx(je,1,jb), &
!          &        z_vort_flx_RBF(je,jk,jb) 
!        ENDIF
!      write(*,*)'horz adv:grad ekin:vort-flx:',je,jk,jb, vn(je,1,jb),z_vort_flx(je,1,jb), &
!        &        veloc_adv_horz_e(je,1,jb) !, p_diag%grad(je,1,jb)
       !Store mapped value
!        z_tmp = veloc_adv_horz_e(je,jk,jb)
! 
!       !get adjacent triangles
!       il_c1 = p_patch%edges%cell_idx(je,jb,1)
!       ib_c1 = p_patch%edges%cell_blk(je,jb,1)
! 
!       il_c2 = p_patch%edges%cell_idx(je,jb,2)
!       ib_c2 = p_patch%edges%cell_blk(je,jb,2)
! 
!       !check if triangle 1 has a land edge
!       !if no do nothing, if yes replace ptp_vn by original values 
!       IF(    p_patch%patch_oce%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary&
!        &.AND.p_patch%patch_oce%lsm_oce_c(il_c2,jk,ib_c2) == sea)THEN
! 
!         !get index of cell1  edges
!         il_e1 = p_patch%cells%edge_idx(il_c1,ib_c1,1)
!         ib_e1 = p_patch%cells%edge_blk(il_c1,ib_c1,1)
! 
!         il_e2 = p_patch%cells%edge_idx(il_c1,ib_c1,2)
!         ib_e2 = p_patch%cells%edge_blk(il_c1,ib_c1,2)
! 
!         il_e3 = p_patch%cells%edge_idx(il_c1,ib_c1,3)
!         ib_e3 = p_patch%cells%edge_blk(il_c1,ib_c1,3)
!         !reset all mapped values to their original, except the 
!         !edge that connects both 
!         veloc_adv_horz_e(il_e1,jk,ib_e1) = z_vort_flx_RBF(il_e1,jk,ib_e1)
!         veloc_adv_horz_e(il_e2,jk,ib_e2) = z_vort_flx_RBF(il_e2,jk,ib_e2)
!         veloc_adv_horz_e(il_e3,jk,ib_e3) = z_vort_flx_RBF(il_e3,jk,ib_e3)
! 
!        veloc_adv_horz_e(je,jk,jb)=z_tmp
! 
!       ELSEIF(    p_patch%patch_oce%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary&
!        &.AND.p_patch%patch_oce%lsm_oce_c(il_c2,jk,ib_c2) == sea)THEN
! 
!         !get index of cell2  edges
!         il_e1 = p_patch%cells%edge_idx(il_c2,ib_c2,1)
!         ib_e1 = p_patch%cells%edge_blk(il_c2,ib_c2,1)
! 
!         il_e2 = p_patch%cells%edge_idx(il_c2,ib_c2,2)
!         ib_e2 = p_patch%cells%edge_blk(il_c2,ib_c2,2)
! 
!         il_e3 = p_patch%cells%edge_idx(il_c2,ib_c2,3)
!         ib_e3 = p_patch%cells%edge_blk(il_c2,ib_c2,3)
!         !reset all mapped values to their original, except the 
!         !edge that connects both 
!         veloc_adv_horz_e(il_e1,jk,ib_e1) =&
!         & z_vort_flx_RBF(il_e1,jk,ib_e1)+ z_grad_ekin_RBF(il_e1,jk,ib_e1)
!         veloc_adv_horz_e(il_e2,jk,ib_e2) =&
!         & z_vort_flx_RBF(il_e2,jk,ib_e2)+ z_grad_ekin_RBF(il_e2,jk,ib_e2)
!         veloc_adv_horz_e(il_e3,jk,ib_e3) =&
!         & z_vort_flx_RBF(il_e3,jk,ib_e3)+ z_grad_ekin_RBF(il_e3,jk,ib_e3)
! 
!        veloc_adv_horz_e(je,jk,jb)=z_tmp
! 
!       ELSEIF(    p_patch%patch_oce%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary&
!        &.AND.p_patch%patch_oce%lsm_oce_c(il_c2,jk,ib_c2) == sea_boundary)THEN
! 
!         il_e1 = p_patch%cells%edge_idx(il_c1,ib_c1,1)
!         ib_e1 = p_patch%cells%edge_blk(il_c1,ib_c1,1)
! 
!         il_e2 = p_patch%cells%edge_idx(il_c1,ib_c1,2)
!         ib_e2 = p_patch%cells%edge_blk(il_c1,ib_c1,2)
! 
!         il_e3 = p_patch%cells%edge_idx(il_c1,ib_c1,3)
!         ib_e3 = p_patch%cells%edge_blk(il_c1,ib_c1,3)
! 
!         veloc_adv_horz_e(il_e1,jk,ib_e1) =&
!         & z_vort_flx_RBF(il_e1,jk,ib_e1)+ z_grad_ekin_RBF(il_e1,jk,ib_e1)
!         veloc_adv_horz_e(il_e2,jk,ib_e2) =&
!         & z_vort_flx_RBF(il_e2,jk,ib_e2)+ z_grad_ekin_RBF(il_e2,jk,ib_e2)
!         veloc_adv_horz_e(il_e3,jk,ib_e3) =&
!         & z_vort_flx_RBF(il_e3,jk,ib_e3)+ z_grad_ekin_RBF(il_e3,jk,ib_e3)
! 
!         il_e1 = p_patch%cells%edge_idx(il_c2,ib_c2,1)
!         ib_e1 = p_patch%cells%edge_blk(il_c2,ib_c2,1)
! 
!         il_e2 = p_patch%cells%edge_idx(il_c2,ib_c2,2)
!         ib_e2 = p_patch%cells%edge_blk(il_c2,ib_c2,2)
! 
!         il_e3 = p_patch%cells%edge_idx(il_c2,ib_c2,3)
!         ib_e3 = p_patch%cells%edge_blk(il_c2,ib_c2,3)
! 
!         veloc_adv_horz_e(il_e1,jk,ib_e1) =&
!         & z_vort_flx_RBF(il_e1,jk,ib_e1)+ z_grad_ekin_RBF(il_e1,jk,ib_e1)
!         veloc_adv_horz_e(il_e2,jk,ib_e2) =&
!         & z_vort_flx_RBF(il_e2,jk,ib_e2)+ z_grad_ekin_RBF(il_e2,jk,ib_e2)
!         veloc_adv_horz_e(il_e3,jk,ib_e3) =&
!         & z_vort_flx_RBF(il_e3,jk,ib_e3)+ z_grad_ekin_RBF(il_e3,jk,ib_e3)
!      ENDIF

   END DO
  END DO
END DO
!    do jb=1,10!i_startblk_e, i_endblk_e
!     CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
!                        i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
!    do je= i_startidx_e, i_endidx_e
!   IF(je==5)THEN
! ! !     IF(vn(je,1,jb)/=0.0)THEN
!       write(*,*)'RBF: vn_in, grad_ekin, vort_flx:',je,jb, vn(je,1,jb), &
!        &             z_grad_ekin_RBF(je,1,jb), &
!        &            -z_vort_flx_RBF(je,1,jb), z_veloc_adv_horz_e(je,1,jb)
!       write(*,*)'PTP: vn_in, grad_ekin, vort_flx:',je,jb, vn(je,1,jb), &
!         &             p_diag%grad(je,1,jb), &
!         &             z_vort_flx(je,1,jb), veloc_adv_horz_e(je,1,jb)
!       ENDIF
! ! !  ENDIF
!    enddo
!    enddo
!  stop
END subroutine veloc_adv_horz_mimetic
!-------------------------------------------------------------------------
!
!
!>
!! Computes horizontal advection of a (edge based) vector field.
!!
!! Computes rotational term of a vector field given by its components in
!! the directions normal to triangle edges and the gradient of the kinetic energy
!! which is calculated using the reconstructed velocity at cell centers. Both
!! terms are combined and constitute the horizontal velocity advection.
!!
!!IMPORTANT: It is assumed that the reconstruction of the tangential velocity
!!           has been done before.
!1
!! input:  lives on edges (velocity points)
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE veloc_adv_horz_RBF( p_patch, vn, p_diag, veloc_adv_horz_e, p_int)
!
!
!  patch on which computation is performed
!
TYPE(t_patch), INTENT(in) :: p_patch

!
! normal and tangential velocity  of which advection is computed
!
REAL(wp), INTENT(inout) :: vn(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
!
!diagnostic ocean state stores horizontally advected velocity
!
TYPE(t_hydro_ocean_diag) :: p_diag
!
! variable in which horizontally advected velocity is stored
!
REAL(wp), INTENT(out) :: veloc_adv_horz_e(:,:,:)
!
!Interpolation necessary just for testing
TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int


INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jk, jb, je, jc
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: rl_start_e, rl_end_e, rl_start_c, rl_end_c!, rl_start_c, rl_end_c
INTEGER :: i_v1_idx, i_v1_blk, i_v2_idx, i_v2_blk
INTEGER :: jev, ile,ibe, i_v1_ctr, i_v2_ctr 
REAL(wp) :: z_e  (nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_glb(nproma,n_zlev,p_patch%nblks_v)
!REAL(wp) :: z_y
REAL(wp) :: z_grad_ekin_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_RBF_e(nproma,n_zlev,p_patch%nblks_e)
INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3
!ARRAYS FOR TESTING
REAL(wp) :: z_vort_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_flx_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_e_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_weight_e1,z_weight_e2, z_weight_e3!, z_weight
!REAL(wp) :: z_beta_plane_vort
!REAL(wp) :: z_f_plane_vort
!REAL(wp) :: z_lat_basin_center
!-----------------------------------------------------------------------
! #slo# set local variable to zero due to nag -nan compiler-option
z_e             (:,:,:) = 0.0_wp
z_vort_glb      (:,:,:) = 0.0_wp
!z_vort_flx      (:,:,:) = 0.0_wp
veloc_adv_horz_e(:,:,:) = 0.0_wp
z_vort_e        (:,:,:) = 0.0_wp
z_vort_flx_RBF  (:,:,:) = 0.0_wp
z_kin_e_RBF     (:,:,:) = 0.0_wp
z_grad_ekin_RBF (:,:,:) = 0.0_wp

rl_start_e = 1
rl_end_e   = min_rledge
rl_start_c = 1
rl_end_c   = min_rlcell

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

slev = 1
elev = n_zlev

CALL rbf_vec_interpol_edge( vn,       &
                          & p_patch,  &
                          & p_int,    &
                          & p_diag%vt,&
                          & opt_slev=slev,opt_elev=elev)
DO jb = i_startblk_e,i_endblk_e
  CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
    &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
  DO jk = slev, elev
    DO je=i_startidx_e, i_endidx_e
      IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) == boundary ) THEN
        p_diag%vt(je,jk,jb) = 0.0_wp
               vn(je,jk,jb) = 0.0_wp
      ENDIF
    END DO
  END DO
END DO

!CALL rot_vertex_ocean(vn, p_diag%vt, p_patch, p_diag%vort)
CALL rot_vertex_ocean_origin(vn, p_diag%vt, p_patch, p_diag%vort)
!CALL rot_vertex_ocean_total(vn, p_patch, p_diag%vort)
! CALL verts2edges_scalar( p_diag%vort, p_patch, p_int%v_1o2_e, &
!                          z_vort_e, opt_slev=slev,opt_elev=elev, opt_rlstart=3)

  DO jb = i_startblk_e,i_endblk_e
    CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
      &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
    DO jk = slev, slev
      DO je=i_startidx_e, i_endidx_e
        i_v1_idx = p_patch%edges%vertex_idx(je,jb,1)
        i_v1_blk = p_patch%edges%vertex_blk(je,jb,1)
        i_v2_idx = p_patch%edges%vertex_idx(je,jb,2)
        i_v2_blk = p_patch%edges%vertex_blk(je,jb,2)
        !count wet edges in vertex 1
        i_v1_ctr = 0
        DO jev = 1, p_patch%verts%num_edges(i_v1_idx,i_v1_blk)
          ile = p_patch%verts%edge_idx(i_v1_idx,i_v1_blk,jev)
          ibe = p_patch%verts%edge_blk(i_v1_idx,i_v1_blk,jev)
          IF ( p_patch%patch_oce%lsm_oce_e(ile,jk,ibe) == sea ) THEN
            i_v1_ctr = i_v1_ctr +1
          ENDIF
        END DO
        !count wet edges in vertex 2
        i_v2_ctr = 0
        DO jev = 1, p_patch%verts%num_edges(i_v2_idx,i_v2_blk)
          ile = p_patch%verts%edge_idx(i_v2_idx,i_v2_blk,jev)
          ibe = p_patch%verts%edge_blk(i_v2_idx,i_v2_blk,jev)
          IF ( p_patch%patch_oce%lsm_oce_e(ile,jk,ibe) == sea ) THEN
            i_v2_ctr = i_v2_ctr +1
          ENDIF
        END DO
         IF(   i_v1_ctr==p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
         &.AND.i_v2_ctr==p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN

           z_vort_e(je,jk,jb) =&
           & 0.5_wp*(p_diag%vort(i_v1_idx,jk,i_v1_blk)&
           &+        p_diag%vort(i_v2_idx,jk,i_v2_blk))

         ELSEIF(   i_v1_ctr==p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
         &.AND.i_v2_ctr <p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN

           z_vort_e(je,jk,jb) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
           &+ REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)

         ELSEIF(   i_v1_ctr<p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
         &.AND.i_v2_ctr==p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN

           z_vort_e(je,jk,jb) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
           &+ REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)

         ELSEIF(   i_v1_ctr<p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
         &.AND.i_v2_ctr<p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN

           z_vort_e(je,jk,jb) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
           &+ REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
         ELSE
           z_vort_e(je,jk,jb) = 0.0_wp
        ENDIF
!  IF(jk==1)THEN
!  IF(p_patch%patch_oce%lsm_oce_e(je,jk,jb)/=2)THEN
! ! IF(p_patch%patch_oce%lsm_oce_e(je,jk,jb)==0)THEN
! IF(i_v1_ctr==6.and.i_v2_ctr==6.and.p_patch%patch_oce%lsm_oce_e(je,jk,jb)==sea)THEN
! ELSE
!  write(101,*)'vert ctr', jk,je,jb, i_v1_ctr, i_v2_ctr, p_patch%patch_oce%lsm_oce_e(je,jk,jb)!,&
! ENDIF
! ! ENDIF
!  ENDIF
!  ENDIF
      END DO
    END DO
  ENDDO

  DO jb = i_startblk_e,i_endblk_e
    CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
      &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
    DO jk = slev, elev
      DO je=i_startidx_e, i_endidx_e
        IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) == sea ) THEN
         z_vort_flx_RBF(je,jk,jb) =&
          & p_diag%vt(je,jk,jb)*(z_vort_e(je,jk,jb)+ p_patch%edges%f_e(je,jb))
        ELSE
          z_vort_flx_RBF(je,jk,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
  ENDDO

CALL rbf_vec_interpol_cell( vn, p_patch, p_int, p_diag%u,  &
  &                         p_diag%v, opt_slev=slev, opt_elev=elev)
  !write(*,*)'max/min vort flux:', MAXVAL(z_vort_flx_RBF(:,1,:)),MINVAL(z_vort_flx_RBF(:,1,:)) 
DO jb = i_startblk_e,i_endblk_e
  CALL get_indices_e(p_patch, jb, i_startblk_e,i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e
      ! calculate kinetic energy at edges from normal and tangential comp.
      z_kin_RBF_e(je,jk,jb) =0.5_wp*(p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb)&
                            &          +    vn(je,jk,jb)*vn(je,jk,jb) )
    ENDDO
  ENDDO
ENDDO
!  write(*,*)'max/min kin energy edges:', MAXVAL(z_kin_RBF_e(:,1,:)),&
!                                          &MINVAL(z_kin_RBF_e(:,1,:)) 
!!$OMP END DO
!!$OMP END PARALLEL
  ! Bilinear interpolation of kinetic energy from the edges to the cells
!    CALL edges2cells_scalar( z_kin_RBF_e,     &
!                           & p_patch,         &
!                           & p_int%e_bln_c_s, &
!                           & p_diag%kin,      &
!                           & opt_slev=slev,opt_elev=elev)
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
                     i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
  DO jk = slev, elev
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) > sea_boundary ) THEN
          p_diag%kin(jc,jk,jb) = 0.0_wp
      ELSE

        ile1 = p_patch%cells%edge_idx(jc,jb,1)
        ibe1 = p_patch%cells%edge_blk(jc,jb,1)
        ile2 = p_patch%cells%edge_idx(jc,jb,2)
        ibe2 = p_patch%cells%edge_blk(jc,jb,2)
        ile3 = p_patch%cells%edge_idx(jc,jb,3)
        ibe3 = p_patch%cells%edge_blk(jc,jb,3)
        z_weight_e1 = 0.0_wp
        z_weight_e2 = 0.0_wp
        z_weight_e3 = 0.0_wp
        IF(p_patch%patch_oce%lsm_oce_e(ile1,jk,ibe1)<= boundary)THEN
          z_weight_e1 = p_patch%edges%area_edge(ile1,ibe1)
        ENDIF
        IF(p_patch%patch_oce%lsm_oce_e(ile2,jk,ibe2)<= boundary)THEN
          z_weight_e2 = p_patch%edges%area_edge(ile2,ibe2)
        ENDIF
        IF(p_patch%patch_oce%lsm_oce_e(ile3,jk,ibe3)<= boundary)THEN
          z_weight_e3 = p_patch%edges%area_edge(ile3,ibe3)
        ENDIF

        !write(*,*)'weights',jc,jk,jb,z_weight_e1,z_weight_e2,z_weight_e3
         p_diag%kin(jc,jk,jb) = (z_kin_RBF_e(ile1,jk,ibe1)*z_weight_e1&
                             &+ z_kin_RBF_e(ile2,jk,ibe2)*z_weight_e2&
                             &+ z_kin_RBF_e(ile3,jk,ibe3)*z_weight_e3)&
                             &/(z_weight_e1+z_weight_e2+z_weight_e3) 
!       p_diag%kin(jc,jk,jb) = 0.5_wp*(p_diag%u(jc,jk,jb)*p_diag%u(jc,jk,jb)&
!                                     &+p_diag%v(jc,jk,jb)*p_diag%v(jc,jk,jb))
      ENDIF
    END DO
  END DO
END DO

CALL grad_fd_norm_oce( p_diag%kin, &
                     & p_patch,    &
                     & z_grad_ekin_RBF, opt_slev=slev,opt_elev=elev)

!Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
DO jb = i_startblk_e, i_endblk_e
  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e
      IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
        veloc_adv_horz_e(je,jk,jb) =&
        & z_vort_flx_RBF(je,jk,jb) + z_grad_ekin_RBF(je,jk,jb) 
      ELSE
        veloc_adv_horz_e(je,jk,jb) = 0.0_wp
      ENDIF 
   END DO
  END DO
END DO

DO jk = slev, elev
 write(*,*)'max/min kin energy:             ', jk,MAXVAL(p_diag%kin(:,jk,:)),&
                                                 &MINVAL(p_diag%kin(:,jk,:)) 
 write(*,*)'max/min grad kin energy:        ', jk,MAXVAL(z_grad_ekin_RBF(:,jk,:)),&
                                                 &MINVAL(z_grad_ekin_RBF(:,jk,:)) 
 write(876,*)'max/min kin energy:           ', jk,MAXVAL(p_diag%kin(:,jk,:)),&
                                                 &MINVAL(p_diag%kin(:,jk,:)) 
 write(876,*)'max/min grad kin energy:      ', jk,MAXVAL(z_grad_ekin_RBF(:,jk,:)),&
                                                 &MINVAL(z_grad_ekin_RBF(:,jk,:)) 
 write(*,*)'max/min vorticity:              ', jk,MAXVAL(p_diag%vort(:,jk,:)),&
                                                 &MINVAL(p_diag%vort(:,jk,:)) 
 write(876,*)'max/min vorticity:            ', jk,MAXVAL(p_diag%vort(:,jk,:)),&
                                                 &MINVAL(p_diag%vort(:,jk,:)) 
!  write(*,*)'max/min vorticity edges:', jk,MAXVAL(z_vort_e(:,jk,:)),&
!                                          &MINVAL(z_vort_e(:,jk,:)) 
!write(876,*)'max/min vorticity edges:', jk,MAXVAL(z_vort_e(:,jk,:)),MINVAL(z_vort_e(:,jk,:)) 
END DO


DO jk = slev, elev
 write(*,*)'max/min vort flux:  advection:  ', jk,&
& MAXVAL(z_vort_flx_RBF(:,jk,:)),MINVAL(z_vort_flx_RBF(:,jk,:)),&
& MAXVAL(veloc_adv_horz_e(:,jk,:)),MINVAL(veloc_adv_horz_e(:,jk,:))

 write(876,*)'max/min vort flux:  advection:', jk,&
& MAXVAL(z_vort_flx_RBF(:,jk,:)),MINVAL(z_vort_flx_RBF(:,jk,:)),&
& MAXVAL(veloc_adv_horz_e(:,jk,:)),MINVAL(veloc_adv_horz_e(:,jk,:))
 END DO
END subroutine veloc_adv_horz_RBF
!-------------------------------------------------------------------------
!
!
!>
!! Computes vertical advection of a (edge based) horizontal vector field.
!! The vertical derivative of the velocity vector at circumcenters that
!! is reconstructed from edge data is calculated and then multiplied by
!! the vertical velocity. The product is mapped from top of the computational
!! prism to the middle (still at centers) via the transposed of vertical differentiation
!! and then transformed to edges.
!!
!! IMPORTANT: It is assumed that the velocity vector reconstruction from
!! edges to cells has been done before.
!!
!! input:  lives on cells (velocity points)
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE veloc_adv_vert_mimetic( p_patch, p_aux, p_diag,&
&                          top_bc_w_c,  bot_bc_w_c,&
&                          veloc_adv_vert_e)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch
TYPE(t_hydro_ocean_aux)           :: p_aux
TYPE(t_hydro_ocean_diag)          :: p_diag
!
REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! variable in which horizontally advected velocity is stored
REAL(wp), INTENT(inout) :: veloc_adv_vert_e(:,:,:)

!local variables
INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: z_dolic
!REAL(wp) :: dzi, dzi_m1, dzi_p1
!REAL(wp) :: dzm, dzm_m1, dzm_p1
TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &                              z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
!-----------------------------------------------------------------------
! #slo# set local variable to zero due to nag -nan compiler-option
! z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c)%x = 0.0_wp
! z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)%x   = 0.0_wp

! blocking
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev
DO jk = slev, elev
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      z_adv_u_i(jc,jk,jb)%x = 0.0_wp
      z_adv_u_m(jc,jk,jb)%x = 0.0_wp
    END DO
  END DO
END DO
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_adv_u_i(jc,elev+1,jb)%x = 0.0_wp
  END DO
END DO

!Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
!This requires appropriate boundary conditions
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    ! #slo# 2011-05-11 - where is  1a) bc_top_veloc?
    ! 1b) ocean interior 
    DO jk = slev+1, z_dolic-1
        z_adv_u_i(jc,jk,jb)%x&
          & = p_diag%w(jc,jk-1,jb)*(p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x)&
          & / p_patch%patch_oce%del_zlev_i(jk-1)
! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb),&
!&( p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x )
    END DO
    ! 1c) ocean bottom
    ! #slo# 2011-05-12 - is w(0) defined? Secure for dolic>1 at wet points
    IF ( z_dolic>0 ) &  ! wet points only 
      & z_adv_u_i(jc,z_dolic,jb)%x = p_diag%w(jc,z_dolic-1,jb)*p_aux%bc_bot_veloc_cc(jc,jb)%x

  END DO
! write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
! & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
END DO

! ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! ! This mapping is the transposed of the vertical differencing.
! 
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
    &                i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    ! 2b) ocean interior
    DO jk = slev,z_dolic-1
        z_adv_u_m(jc,jk,jb)%x &
        & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)%x&
        & +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_u_i(jc,jk,jb)%x) &
        & / (p_patch%patch_oce%del_zlev_i(jk+1)+p_patch%patch_oce%del_zlev_i(jk))
    END DO
    ! 2c) ocean bottom
    IF ( z_dolic>0 ) &  ! wet points only
      &  z_adv_u_m(jc,z_dolic,jb)%x = 0.0_wp
! #slo# 2011-05-11 - Attention: dolic (not elev) must be used here for non-zero values
!     & = (0.5_wp*p_patch%patch_oce%del_zlev_m(z_dolic)*z_adv_u_i(jc,z_dolic+1,jb)%x&
!     & +         p_patch%patch_oce%del_zlev_i(z_dolic)*z_adv_u_i(jc,z_dolic,  jb)%x)&
!     & / (2.0_wp*p_patch%patch_oce%del_zlev_m(z_dolic))
  END DO
! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
END DO

! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
CALL map_cell2edges( p_patch, z_adv_u_m, veloc_adv_vert_e, &
  &                  opt_slev=slev, opt_elev=elev )

!  DO jk=1,n_zlev
!    WRITE(*,*) 'max/min vert adv FINAL',jk, &
!      &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
!  END DO

END subroutine veloc_adv_vert_mimetic
! ! !-------------------------------------------------------------------------
!
!
!>
!! Computes vertical advection of a (edge based) horizontal vector field.
!! The vertical derivative of the velocity vector at circumcenters that
!! is reconstructed from edge data is calculated and then multiplied by
!! the vertical velocity. The product is mapped from top of the computational
!! prism to the middle (still at centers) via the transposed of vertical differentiation
!! and then transformed to edges.
!!
!! IMPORTANT: It is assumed that the velocity vector reconstruction from
!! edges to cells has been done before.
!!
!! input:  lives on cells (velocity points)
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE veloc_adv_diff_vert_RBF( p_patch, u_c, v_c, w_c, &
&                          top_bc_u_c, top_bc_v_c,          &
&                          bot_bc_u_c, bot_bc_v_c,          &
&                          A_v,                             & 
&                          veloc_adv_vert_e)
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch

!
! Components of cell based variable which is vertically advected
REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
!
! Top boundary condition for cell based variables
REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!
! Bottom boundary condition for cell based variables
REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: A_v(:,:,:)      ! dim: (nproma,n_zlev+1,nblks_c)
!
! variable in which horizontally advected velocity is stored
REAL(wp), INTENT(out) :: veloc_adv_vert_e(:,:,:)

!INTEGER, PARAMETER :: top=1
INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb, z_dolic
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
REAL(wp) :: dzi, dzi_m1!, dzi_p1
REAL(wp) :: dzm, dzm_m1!, dzm_p1

REAL(wp) :: z_adv_diff_u(nproma,n_zlev,p_patch%nblks_c),  &
  &         z_adv_diff_v(nproma,n_zlev,p_patch%nblks_c),  &
  &         z_diff_u(nproma,n_zlev+1,p_patch%nblks_c),    &
  &         z_diff_v(nproma,n_zlev+1,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: zu_cc(nproma,n_zlev,p_patch%nblks_c)
!-----------------------------------------------------------------------

! #slo# set local variable to zero due to nag -nan compiler-option
z_adv_diff_u(:,:,:) = 0.0_wp
z_adv_diff_v(:,:,:) = 0.0_wp

! blocking
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev

!Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
!This requires appropriate boundary conditions

DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    DO jk = slev+1, z_dolic-1
      !check if we have at least two layers of water
      !  #slo# - 2011-04-01 - Is this really intended here
      !  maybe this condition should be fulfilled everywhere
      !  then it must be calculated in fill_vertical_domain
      !  this condition could then be omitted here
      IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN 
      !IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        !ocean interior 
        dzi_m1 = p_patch%patch_oce%del_zlev_i(jk-1)

        z_diff_u(jc,jk,jb)&
        & = A_v(jc,jk,jb)*( u_c(jc,jk-1,jb) - u_c(jc,jk,jb) )&
        & / dzi_m1

        z_diff_v(jc,jk,jb)&
        & = A_v(jc,jk,jb)*( v_c(jc,jk-1,jb) - v_c(jc,jk,jb) )&
        & / dzi_m1
     ENDIF    ! at least 2 vertical layers
    END DO
  END DO
END DO

!Including top&bottom boundary condition
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)

    z_diff_u(jc,1,jb)       = top_bc_u_c(jc,jb)
    z_diff_v(jc,1,jb)       = top_bc_v_c(jc,jb)

    z_diff_u(jc,z_dolic,jb) = bot_bc_u_c(jc,jb)
    z_diff_v(jc,z_dolic,jb) = bot_bc_v_c(jc,jb)
  END DO
END DO

! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! This mapping is the transposed of the vertical differencing.

DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    DO jk = slev, z_dolic
      !distance between prism centers at levels jk and (jk-1)
      dzi   = p_patch%patch_oce%del_zlev_i(jk)
      !thickness of prism at level jk
      dzm   = p_patch%patch_oce%del_zlev_m(jk)

      !check if we are on land: To be replaced by 3D lsm  
      ! #slo# 2011-05-11 - replace by consistent formulation: vertical loop down to dolic
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
       IF ( jk == slev ) THEN
         z_adv_diff_u(jc,jk,jb) =&
         & (z_diff_u(jc,jk,jb)-z_diff_u(jc,jk+1,jb))/dzm

         z_adv_diff_v(jc,jk,jb) =&
         & (z_diff_v(jc,jk,jb)-z_diff_v(jc,jk+1,jb))/dzm

       ELSEIF(jk>slev .AND. jk < z_dolic ) THEN
         !thickness of prism at level (jk-1)
         dzm_m1= p_patch%patch_oce%del_zlev_m(jk-1)
         !distance between prism centers at levels (jk-1) and (jk-2)
         dzi_m1= p_patch%patch_oce%del_zlev_i(jk-1) 

         z_adv_diff_u(jc,jk,jb) &
         & = -(( u_c(jc,jk-1,jb)- u_c(jc,jk,jb))*w_c(jc,jk-1,jb)&
         &/dzi_m1                                               &
         & - ( u_c(jc,jk,jb) - u_c(jc,jk+1,jb))*w_c(jc,jk,jb)   &
         & /dzi)                                                &
         & / (2.0_wp*dzm)                                       &
         & + (z_diff_u(jc,jk,jb)-z_diff_u(jc,jk+1,jb))/dzm

         z_adv_diff_v(jc,jk,jb) &
         & = -(( v_c(jc,jk-1,jb)- v_c(jc,jk,jb))*w_c(jc,jk-1,jb)&
         &/dzi_m1                                               &
         & - ( v_c(jc,jk,jb) - v_c(jc,jk+1,jb))*w_c(jc,jk,jb)   &
         &/dzi)                                                 &
         & / (2.0_wp*dzm)                                       &
         & + (z_diff_v(jc,jk,jb)-z_diff_v(jc,jk+1,jb))/dzm

       ELSEIF ( jk == z_dolic ) THEN
         dzm_m1= p_patch%patch_oce%del_zlev_m(jk-1)
         dzi_m1= p_patch%patch_oce%del_zlev_i(jk-1)

         z_adv_diff_u(jc,jk,jb) =&
         &  (z_diff_u(jc,jk-1,jb)-z_diff_u(jc,jk,jb))/dzm_m1

         z_adv_diff_v(jc,jk,jb) =&
         &  (z_diff_v(jc,jk-1,jb)-z_diff_v(jc,jk,jb))/dzm_m1
       ENDIF
         CALL gvec2cvec( -z_adv_diff_u(jc,jk,jb), -z_adv_diff_v(jc,jk,jb), &
         &               p_patch%cells%center(jc,jb)%lon,                &
         &               p_patch%cells%center(jc,jb)%lat,                &
         &               zu_cc(jc,jk,jb)%x(1),zu_cc(jc,jk,jb)%x(2),zu_cc(jc,jk,jb)%x(3))
      ENDIF
    END DO
  END DO
!  write(*,*)'B max/min vert adv:',jk, maxval(z_adv_diff_u(:,jk,:)),&
!  & minval(z_adv_diff_u(:,jk,:)),&
!  & maxval(z_adv_diff_v(:,jk,:)),&
!  & minval(z_adv_diff_v(:,jk,:))
END DO

! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
CALL map_cell2edges( p_patch, zu_cc, veloc_adv_vert_e)

END subroutine veloc_adv_diff_vert_RBF
! ! !-------------------------------------------------------------------------
!
!
!>
!! Computes vertical advection of a (edge based) horizontal vector field.
!! The vertical derivative of the velocity vector at circumcenters that
!! is reconstructed from edge data is calculated and then multiplied by
!! the vertical velocity. The product is mapped from top of the computational
!! prism to the middle (still at centers) via the transposed of vertical differentiation
!! and then transformed to edges.
!!
!! IMPORTANT: It is assumed that the velocity vector reconstruction from
!! edges to cells has been done before.
!!
!! input:  lives on cells (velocity points)
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE veloc_adv_vert_RBF( p_patch, u_c, v_c, w_c, &
&                          top_bc_u_c, top_bc_v_c, &
&                          bot_bc_u_c,  bot_bc_v_c,&
&                          top_bc_w_c,  bot_bc_w_c,&
&                          veloc_adv_vert_e)
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch

!
! Components of cell based variable which is vertically advected
REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
!
! Top boundary condition for cell based variables
REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!
! Bottom boundary condition for cell based variables
REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!
REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)

! variable in which horizontally advected velocity is stored
REAL(wp), INTENT(out) :: veloc_adv_vert_e(:,:,:)

!INTEGER, PARAMETER :: top=1
INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb, i_dolic
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_v_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_u_m(nproma,n_zlev,p_patch%nblks_c),  &
  &         z_adv_v_m(nproma,n_zlev,p_patch%nblks_c)
!-----------------------------------------------------------------------

! #slo# set local variable to zero due to nag -nan compiler-option
z_adv_u_i(:,:,:) = 0.0_wp
z_adv_v_i(:,:,:) = 0.0_wp
z_adv_u_m(:,:,:) = 0.0_wp
z_adv_v_m(:,:,:) = 0.0_wp

! blocking
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev

!Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
!This requires appropriate boundary conditions
DO jk = slev, elev
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      !check if we have at least two layers of water
      !  #slo# - 2011-04-01 - Is this really intended here
      !  maybe this condition should be fulfilled everywhere
      !  then it must be calculated in fill_vertical_domain
      !  this condition could then be omitted here
      IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
        !1b) ocean bottom 
        IF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
          ! u,v-component
          z_adv_u_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_u_c(jc,jb)
          z_adv_v_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_v_c(jc,jb)

        !1c) ocean interior 
        ELSEIF( jk>slev .AND.  jk < p_patch%patch_oce%dolic_c(jc,jb))THEN
          ! u,v-component
          z_adv_u_i(jc,jk,jb)&
          & = w_c(jc,jk-1,jb) *( u_c(jc,jk-1,jb) - u_c(jc,jk,jb) )&
            & / p_patch%patch_oce%del_zlev_i(jk-1)

          z_adv_v_i(jc,jk,jb)&
          & = w_c(jc,jk-1,jb) *( v_c(jc,jk-1,jb) - v_c(jc,jk,jb) )&
          & / p_patch%patch_oce%del_zlev_i(jk-1) !&
! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb) *( v_c(jc,jk,jb) - v_c(jc,jk-1,jb) ),&
! &  v_c(jc,jk,jb), v_c(jc,jk-1,jb),p_patch%patch_oce%del_zlev_i(jk-1) 
! write(*,*)'vert adv:u: ',jk, jc,jb,w_c(jc,jk,jb) *( u_c(jc,jk,jb) - u_c(jc,jk-1,jb) ),&
! &  u_c(jc,jk,jb), u_c(jc,jk-1,jb)

        ENDIF  ! jk-condition
      ENDIF    ! at least 2 vertical layers
    END DO
  END DO
!  write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
!  & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
END DO

! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! This mapping is the transposed of the vertical differencing.

!1) From surface down to one layer before bottom
DO jk = slev, elev-1
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      !check if we are on land: To be replaced by 3D lsm  
      ! #slo# 2011-05-11 - replace by consistent formulation: vertical loop down to dolic
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

        z_adv_u_m(jc,jk,jb) &
        & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)&
        & +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_u_i(jc,jk,jb)) &
        & / (p_patch%patch_oce%del_zlev_i(jk+1)+p_patch%patch_oce%del_zlev_i(jk))

        z_adv_v_m(jc,jk,jb)&
        & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_v_i(jc,jk+1,jb)&
        &   +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_v_i(jc,jk,jb))&
        & / (p_patch%patch_oce%del_zlev_i(jk+1)+p_patch%patch_oce%del_zlev_i(jk))
      ENDIF
    END DO
  END DO
! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
END DO
!Bottom layer
!The value of p_patch%patch_oce%del_zlev_i at the botom is 0.5*p_patch%patch_oce%del_zlev_m
!The dimensioning of the firs arrays requires to seperate the vertical loop.
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    IF ( p_patch%patch_oce%dolic_c(jc,jb)>0 ) THEN  ! wet points only
       i_dolic = p_patch%patch_oce%dolic_c(jc,jb)
       z_adv_u_m(jc,i_dolic,jb) =0.0_wp!&
      !& = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,jb)&
      !& +        p_patch%patch_oce%del_zlev_i(elev)*z_adv_u_i(jc,elev,jb)) &
      !& / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))

      z_adv_v_m(jc,i_dolic,jb)=0.0_wp!&
      !& = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_v_i(jc,elev+1,jb)&
      !&   +        p_patch%patch_oce%del_zlev_i(elev)*z_adv_v_i(jc,elev,jb))&
      !& / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))

    END IF
  END DO
END DO

! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
CALL primal_map_c2e( p_patch,&
                  & z_adv_u_m, z_adv_v_m,&
                  & veloc_adv_vert_e )
! DO jk=1,n_zlev
!   WRITE(*,*) 'max/min vert adv FINAL',jk, &
!     &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
! END DO

END subroutine veloc_adv_vert_RBF
! !-----------------------------code below no longer used---------------------------------------------
! ! ! !-------------------------------------------------------------------------
! !
! !
! !>
! !! Computes vertical advection of a (edge based) horizontal vector field.
! !! The vertical derivative of the velocity vector at circumcenters that
! !! is reconstructed from edge data is calculated and then multiplied by
! !! the vertical velocity. The product is mapped from top of the computational
! !! prism to the middle (still at centers) via the transposed of vertical differentiation
! !! and then transformed to edges.
! !!
! !! IMPORTANT: It is assumed that the velocity vector reconstruction from
! !! edges to cells has been done before.
! !!
! !! input:  lives on cells (velocity points)
! !! output: lives on edges (velocity points)
! !!
! !! @par Revision History
! !! Developed  by  Peter Korn, MPI-M (2010).
! !!
! SUBROUTINE veloc_adv_vert_RBF2( p_patch, u_c, v_c, w_c, &
! &                          top_bc_u_c, top_bc_v_c, &
! &                          bot_bc_u_c,  bot_bc_v_c,&
! &                          top_bc_w_c,  bot_bc_w_c,&
! &                          veloc_adv_vert_e)
! !
! !  patch on which computation is performed
! !
! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! 
! !
! ! Components of cell based variable which is vertically advected
! REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
! !
! ! Top boundary condition for cell based variables
! REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! !
! ! Bottom boundary condition for cell based variables
! REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! !
! REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! 
! ! variable in which horizontally advected velocity is stored
! REAL(wp), INTENT(out) :: veloc_adv_vert_e(:,:,:)
! 
! !INTEGER, PARAMETER :: top=1
! INTEGER :: slev, elev     ! vertical start and end level
! INTEGER :: jc, jk, jb
! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! 
! REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
!   &         z_adv_v_i(nproma,n_zlev+1,p_patch%nblks_c),  &
!   &         z_adv_u_m(nproma,n_zlev,p_patch%nblks_c),  &
!   &         z_adv_v_m(nproma,n_zlev,p_patch%nblks_c)
! !-----------------------------------------------------------------------
! 
! ! #slo# set local variable to zero due to nag -nan compiler-option
! z_adv_u_i(:,:,:) = 0.0_wp
! z_adv_v_i(:,:,:) = 0.0_wp
! z_adv_u_m(:,:,:) = 0.0_wp
! z_adv_v_m(:,:,:) = 0.0_wp
! 
! ! blocking
! i_startblk = p_patch%cells%start_blk(1,1)
! i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
! slev = 1
! elev = n_zlev
! 
! !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
! !This requires appropriate boundary conditions
! DO jk = slev, elev
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
!       !check if we have at least two layers of water
!       !  #slo# - 2011-04-01 - Is this really intended here
!       !  maybe this condition should be fulfilled everywhere
!       !  then it must be calculated in fill_vertical_domain
!       !  this condition could then be omitted here
!       IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
! 
!         !1a) 0cean surface
!         IF(jk==slev)THEN
!           ! u,v-component
!           z_adv_u_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_u_c(jc,jb)
!           z_adv_v_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_v_c(jc,jb)
! !write(*,*)'vert adv: top:',jc,jb,jk,top_bc_w_c(jc,jb),top_bc_u_c(jc,jb),top_bc_v_c(jc,jb) 
!         !1b) ocean bottom 
!         ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
!           ! u,v-component
!           z_adv_u_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_u_c(jc,jb)
!           z_adv_v_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_v_c(jc,jb)
! 
!         !1c) ocean interior 
!         ELSEIF( jk>slev .AND.  jk < p_patch%patch_oce%dolic_c(jc,jb))THEN
!           ! u,v-component
!           z_adv_u_i(jc,jk,jb)&
!           & = w_c(jc,jk,jb) *( u_c(jc,jk-1,jb) + u_c(jc,jk,jb) )
! 
!           z_adv_v_i(jc,jk,jb)&
!           & = w_c(jc,jk,jb) *( v_c(jc,jk-1,jb) + v_c(jc,jk,jb) )
! ! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb) *( v_c(jc,jk,jb) - v_c(jc,jk-1,jb) ),&
! ! &  v_c(jc,jk,jb), v_c(jc,jk-1,jb),p_patch%patch_oce%del_zlev_i(jk-1) 
! ! write(*,*)'vert adv:u: ',jk, jc,jb,w_c(jc,jk,jb) *( u_c(jc,jk,jb) - u_c(jc,jk-1,jb) ),&
! ! &  u_c(jc,jk,jb), u_c(jc,jk-1,jb)
! 
!         ENDIF  ! jk-condition
!       ENDIF    ! at least 2 vertical layers
!     END DO
!   END DO
! !   write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
! !   & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
! END DO
! 
! ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! ! This mapping is the transposed of the vertical differencing.
! 
! !1) From surface down to one layer before bottom
! DO jk = slev, elev-1
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
!       !check if we are on land: To be replaced by 3D lsm  
!       IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!       !IF (p_patch%patch_oce%dolic_c(jc,jb) <= sea_boundary) THEN
! 
!         z_adv_u_m(jc,jk,jb) &
!         & = (z_adv_u_i(jc,jk,jb) -  z_adv_u_i(jc,jk+1,jb)) &
!         & / p_patch%patch_oce%del_zlev_m(jk)
! 
!         z_adv_v_m(jc,jk,jb)&
!         & = (z_adv_v_i(jc,jk,jb) - z_adv_v_i(jc,jk+1,jb))&
!         & / p_patch%patch_oce%del_zlev_m(jk)
!       ENDIF
!     END DO
!   END DO
! !  write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
! !  & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
! END DO
! !Bottom layer
! !The value of p_patch%patch_oce%del_zlev_i at the botom is 0.5*p_patch%patch_oce%del_zlev_m
! !The dimensioning of the firs arrays requires to seperate the vertical loop.
! DO jb = i_startblk, i_endblk
!   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                      i_startidx, i_endidx, 1,min_rlcell)
!   DO jc = i_startidx, i_endidx
!     !check if we are on land: To be replaced by 3D lsm      
!     IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!     !IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
! 
!       z_adv_u_m(jc,elev,jb) &
!       & = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,jb)&
!       & +  p_patch%patch_oce%del_zlev_i(elev)*z_adv_u_i(jc,elev,jb)) &
!       & / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))
! 
!       z_adv_v_m(jc,elev,jb)&
!       & = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_v_i(jc,elev+1,jb)&
!       &   +  p_patch%patch_oce%del_zlev_i(elev)*z_adv_v_i(jc,elev,jb))&
!       & / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))
! 
!       ENDIF
!     END DO
!   END DO
! 
! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
!  CALL primal_map_c2e( p_patch,&
!                    & z_adv_u_m, z_adv_v_m,&
!                    & veloc_adv_vert_e )
! ! DO jk=1,n_zlev
! !   WRITE(*,*) 'max/min vert adv FINAL',jk, &
! !     &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
! ! END DO
! 
! END subroutine veloc_adv_vert_RBF2
! ! ! !-------------------------------------------------------------------------
! !
! !
! !>
! !! Computes vertical advection of a (edge based) horizontal vector field.
! !! The vertical derivative of the velocity vector at circumcenters that
! !! is reconstructed from edge data is calculated and then multiplied by
! !! the vertical velocity. The product is mapped from top of the computational
! !! prism to the middle (still at centers) via the transposed of vertical differentiation
! !! and then transformed to edges.
! !!
! !! IMPORTANT: It is assumed that the velocity vector reconstruction from
! !! edges to cells has been done before.
! !!
! !! input:  lives on cells (velocity points)
! !! output: lives on edges (velocity points)
! !!
! !! @par Revision History
! !! Developed  by  Peter Korn, MPI-M (2010).
! !!
! SUBROUTINE veloc_adv_vert_mimetic_old( p_patch, p_aux, p_diag,&
! &                          top_bc_w_c,  bot_bc_w_c,&
! &                          veloc_adv_vert_e)
! 
! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! TYPE(t_hydro_ocean_aux)           :: p_aux
! TYPE(t_hydro_ocean_diag)          :: p_diag
! !
! REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! ! variable in which horizontally advected velocity is stored
! REAL(wp), INTENT(inout) :: veloc_adv_vert_e(:,:,:)
! 
! !local variables
! INTEGER :: slev, elev     ! vertical start and end level
! INTEGER :: jc, jk, jb
! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! 
! TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
!   &                              z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
! !-----------------------------------------------------------------------
! ! #slo# set local variable to zero due to nag -nan compiler-option
! ! z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c)%x = 0.0_wp
! ! z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)%x   = 0.0_wp
! 
! ! blocking
! i_startblk = p_patch%cells%start_blk(1,1)
! i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
! slev = 1
! elev = n_zlev
! DO jk = slev, elev
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
!       z_adv_u_i(jc,jk,jb)%x = 0.0_wp
!       z_adv_u_m(jc,jk,jb)%x = 0.0_wp
!     END DO
!   END DO
! END DO
! DO jb = i_startblk, i_endblk
!   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                      i_startidx, i_endidx, 1,min_rlcell)
!   DO jc = i_startidx, i_endidx
!     z_adv_u_i(jc,elev+1,jb)%x = 0.0_wp
!   END DO
! END DO
! 
! 
! !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
! !This requires appropriate boundary conditions
! DO jk = slev, elev
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
!       !check if we are on land: To be replaced by 3D lsm
!       IF (p_patch%patch_oce%dolic_c(jc,jb) /= 0) THEN
!         !1a) 0cean surface
!         IF(jk==slev)THEN
! 
! !          z_adv_u_i(jc,jk,jb)%x = top_bc_w_c(jc,jb)*p_aux%bc_top_veloc_cc(jc,jb)%x
!           z_adv_u_i(jc,jk,jb)%x = p_diag%w(jc,1,jb)*p_aux%bc_top_veloc_cc(jc,jb)%x
!         !1b) ocean bottom 
!         ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
! 
!           z_adv_u_i(jc,jk+1,jb)%x = bot_bc_w_c(jc,jb)*p_aux%bc_bot_veloc_cc(jc,jb)%x
! 
!         !1c) ocean interior 
!         ELSEIF( jk >slev .AND. jk < p_patch%patch_oce%dolic_c(jc,jb))THEN
! 
!           z_adv_u_i(jc,jk,jb)%x&
!             & = p_diag%w(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x)&
!             & / p_patch%patch_oce%del_zlev_i(jk-1)
! ! ! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb),&
! !&( p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x )
! 
!         ENDIF
!       ENDIF
!     END DO
!   END DO
! ! write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
! ! & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
! END DO
! 
! ! ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! ! ! This mapping is the transposed of the vertical differencing.
! ! 
! ! !1) From surface down to one layer before bottom
! DO jk = slev, elev-1
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
!       ! #slo# 2011-04-01 - dolic_c is now correct
!       ! compare also formulation in veloc_adv_vert_RBF l. 765-790
!       IF ( jk <= p_patch%patch_oce%dolic_c(jc,jb) ) THEN
!       !IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!         z_adv_u_m(jc,jk,jb)%x &
!         & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)%x&
!         & +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_u_i(jc,jk,jb)%x) &
!         & / (2.0_wp*p_patch%patch_oce%del_zlev_m(jk))
!       ! & / (p_patch%patch_oce%del_zlev_m(jk+1)+p_patch%patch_oce%del_zlev_m(jk))
!       ENDIF
!     END DO
!   END DO
! ! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
! ! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
! END DO
! 
! ! !Bottom layer
! ! !The value of p_patch%patch_oce%del_zlev_i at the botom is 0.5*p_patch%patch_oce%del_zlev_m
! ! !The dimensioning of the first arrays requires to seperate the vertical loop.
! DO jb = i_startblk, i_endblk
!   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                      i_startidx, i_endidx, 1,min_rlcell)
!   DO jc = i_startidx, i_endidx
!     ! #slo# 2011-04-01 - dolic_c is now correct
!     IF ( elev <= p_patch%patch_oce%dolic_c(jc,jb) ) THEN
!     !IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! 
!       z_adv_u_m(jc,elev,jb)%x &
!       ! #slo# 2011-04-01 - Attention to multiplying by dz(elev) in both terms
!       & = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,jb)%x&
!       & +  p_patch%patch_oce%del_zlev_i(elev)*z_adv_u_i(jc,elev,jb)%x) &
!       & / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))
!       ENDIF
!     END DO
!   END DO
! 
! ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
!   CALL map_cell2edges( p_patch, z_adv_u_m, veloc_adv_vert_e )
! 
! !  DO jk=1,n_zlev
! !    WRITE(*,*) 'max/min vert adv FINAL',jk, &
! !      &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
! !  END DO
! 
! END subroutine veloc_adv_vert_mimetic_old
! ! !-------------------------------------------------------------------------
!
!
!>
!! Computes vertical advection of a (edge based) horizontal vector field.
!! The vertical derivative of the velocity vector at circumcenters that
!! is reconstructed from edge data is calculated and then multiplied by
!! the vertical velocity. The product is mapped from top of the computational
!! prism to the middle (still at centers) via the transposed of vertical differentiation
!! and then transformed to edges.
!!
!! IMPORTANT: It is assumed that the velocity vector reconstruction from
!! edges to cells has been done before.
!!
!! input:  lives on cells (velocity points)
!! output: lives on edges (velocity points)
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE veloc_adv_vert_RBF_old( p_patch, u_c, v_c, w_c, &
&                          top_bc_u_c, top_bc_v_c, &
&                          bot_bc_u_c,  bot_bc_v_c,&
&                          top_bc_w_c,  bot_bc_w_c,&
&                          veloc_adv_vert_e)
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch

!
! Components of cell based variable which is vertically advected
REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
!
! Top boundary condition for cell based variables
REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!
! Bottom boundary condition for cell based variables
REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!
REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)

! variable in which horizontally advected velocity is stored
REAL(wp), INTENT(out) :: veloc_adv_vert_e(:,:,:)

!INTEGER, PARAMETER :: top=1
INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_v_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_u_m(nproma,n_zlev,p_patch%nblks_c),  &
  &         z_adv_v_m(nproma,n_zlev,p_patch%nblks_c)
!-----------------------------------------------------------------------

! #slo# set local variable to zero due to nag -nan compiler-option
z_adv_u_i(:,:,:) = 0.0_wp
z_adv_v_i(:,:,:) = 0.0_wp
z_adv_u_m(:,:,:) = 0.0_wp
z_adv_v_m(:,:,:) = 0.0_wp

! blocking
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev

!Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
!This requires appropriate boundary conditions
DO jk = slev, elev
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      !check if we have at least two layers of water
      !  #slo# - 2011-04-01 - Is this really intended here
      !  maybe this condition should be fulfilled everywhere
      !  then it must be calculated in fill_vertical_domain
      !  this condition could then be omitted here
      IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN

        !1a) 0cean surface
        IF(jk==slev)THEN
          ! u,v-component
!           z_adv_u_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_u_c(jc,jb)
!           z_adv_v_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_v_c(jc,jb)
          z_adv_u_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_u_c(jc,jb)
          z_adv_v_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_v_c(jc,jb)
!write(*,*)'vert adv: top:',jc,jb,jk,top_bc_w_c(jc,jb),top_bc_u_c(jc,jb),top_bc_v_c(jc,jb) 
        !1b) ocean bottom 
        ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
          ! u,v-component
          z_adv_u_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_u_c(jc,jb)
          z_adv_v_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_v_c(jc,jb)

        !1c) ocean interior 
        ELSEIF( jk>slev .AND.  jk < p_patch%patch_oce%dolic_c(jc,jb))THEN
          ! u,v-component
          z_adv_u_i(jc,jk,jb)&
          & = w_c(jc,jk,jb) *( u_c(jc,jk-1,jb) - u_c(jc,jk,jb) )&
            & / p_patch%patch_oce%del_zlev_i(jk-1)

          z_adv_v_i(jc,jk,jb)&
          & = w_c(jc,jk,jb) *( v_c(jc,jk-1,jb) - v_c(jc,jk,jb) )&
          & / p_patch%patch_oce%del_zlev_i(jk-1) !&
! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb) *( v_c(jc,jk,jb) - v_c(jc,jk-1,jb) ),&
! &  v_c(jc,jk,jb), v_c(jc,jk-1,jb),p_patch%patch_oce%del_zlev_i(jk-1) 
! write(*,*)'vert adv:u: ',jk, jc,jb,w_c(jc,jk,jb) *( u_c(jc,jk,jb) - u_c(jc,jk-1,jb) ),&
! &  u_c(jc,jk,jb), u_c(jc,jk-1,jb)

        ENDIF  ! jk-condition
      ENDIF    ! at least 2 vertical layers
    END DO
  END DO
!  write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
!  & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
END DO

! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! This mapping is the transposed of the vertical differencing.

!1) From surface down to one layer before bottom
DO jk = slev, elev-1
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      !check if we are on land: To be replaced by 3D lsm  
      ! #slo# 2011-05-11 - replace by consistent formulation: vertical loop down to dolic
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
      !IF (p_patch%patch_oce%dolic_c(jc,jb) <= sea_boundary) THEN

        z_adv_u_m(jc,jk,jb) &
        & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)&
        & +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_u_i(jc,jk,jb)) &
        & / (2.0_wp*p_patch%patch_oce%del_zlev_m(jk))

        z_adv_v_m(jc,jk,jb)&
        & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_v_i(jc,jk+1,jb)&
        &   +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_v_i(jc,jk,jb))&
        & / (2.0_wp*p_patch%patch_oce%del_zlev_m(jk))
      ENDIF
    END DO
  END DO
! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
END DO
!Bottom layer
!The value of p_patch%patch_oce%del_zlev_i at the botom is 0.5*p_patch%patch_oce%del_zlev_m
!The dimensioning of the firs arrays requires to seperate the vertical loop.
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    !check if we are on land: To be replaced by 3D lsm      
    ! #slo# 2011-05-11 - replace by consistent formulation: vertical loop down to dolic
    IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
    !IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN

      z_adv_u_m(jc,elev,jb) &
      & = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,jb)&
      & +  p_patch%patch_oce%del_zlev_i(elev)*z_adv_u_i(jc,elev,jb)) &
      & / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))

      z_adv_v_m(jc,elev,jb)&
      & = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_v_i(jc,elev+1,jb)&
      &   +  p_patch%patch_oce%del_zlev_i(elev)*z_adv_v_i(jc,elev,jb))&
      & / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))

      ENDIF
    END DO
  END DO

! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
 CALL primal_map_c2e( p_patch,&
                   & z_adv_u_m, z_adv_v_m,&
                   & veloc_adv_vert_e )
! DO jk=1,n_zlev
!   WRITE(*,*) 'max/min vert adv FINAL',jk, &
!     &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
! END DO

END subroutine veloc_adv_vert_RBF_old
! ! !-------------------------------------------------------------------------
! !
! !
! !>
! !! Computes vertical advection of a (edge based) horizontal vector field.
! !! The vertical derivative of the velocity vector at circumcenters that
! !! is reconstructed from edge data is calculated and then multiplied by
! !! the vertical velocity. The product is mapped from top of the computational
! !! prism to the middle (still at centers) via the transposed of vertical differentiation
! !! and then transformed to edges.
! !!
! !! IMPORTANT: It is assumed that the velocity vector reconstruction from
! !! edges to cells has been done before.
! !!
! !! input:  lives on cells (velocity points)
! !! output: lives on edges (velocity points)
! !!
! !! @par Revision History
! !! Developed  by  Peter Korn, MPI-M (2010).
! !!
! SUBROUTINE veloc_adv_diff_vert_RBF2( p_patch, u_c, v_c, w_c, &
! &                          top_bc_u_c, top_bc_v_c,          &
! &                          bot_bc_u_c, bot_bc_v_c,          &
! &                          A_v,                             & 
! &                          veloc_adv_vert_e)
! !
! !  patch on which computation is performed
! !
! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! 
! !
! ! Components of cell based variable which is vertically advected
! REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
! !
! ! Top boundary condition for cell based variables
! REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! !
! ! Bottom boundary condition for cell based variables
! REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! REAL(wp), INTENT(in) :: A_v(:,:,:)      ! dim: (nproma,n_zlev+1,nblks_c)
! !
! ! variable in which horizontally advected velocity is stored
! REAL(wp), INTENT(out) :: veloc_adv_vert_e(:,:,:)
! 
! !INTEGER, PARAMETER :: top=1
! INTEGER :: slev, elev     ! vertical start and end level
! INTEGER :: jc, jk, jb, z_dolic
! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! 
! REAL(wp) :: &
!   &         z_adv_diff_u_m(nproma,n_zlev,p_patch%nblks_c),  &
!   &         z_adv_diff_v_m(nproma,n_zlev,p_patch%nblks_c),  &
!   &         z_diff_u_i(nproma,n_zlev+1,p_patch%nblks_c),    &
!   &         z_diff_v_i(nproma,n_zlev+1,p_patch%nblks_c) 
! TYPE(t_cartesian_coordinates) :: zu_cc(nproma,n_zlev,p_patch%nblks_c)
! !-----------------------------------------------------------------------
! ! #slo# set local variable to zero due to nag -nan compiler-option
! z_adv_diff_u_m(:,:,:) = 0.0_wp
! z_adv_diff_v_m(:,:,:) = 0.0_wp
! 
! ! blocking
! i_startblk = p_patch%cells%start_blk(1,1)
! i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
! slev = 1
! elev = n_zlev
! 
! !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity 
! !This requires appropriate boundary conditions
! DO jk =slev+1, elev-1
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
!       !check if we have at least two layers of water
!       !  #slo# - 2011-04-01 - Is this really intended here
!       !  maybe this condition should be fulfilled everywhere
!       !  then it must be calculated in fill_vertical_domain
!       !  this condition could then be omitted here
!       IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
!         !ocean interior 
!         IF( jk <= p_patch%patch_oce%dolic_c(jc,jb))THEN
! 
!           z_diff_u_i(jc,jk,jb)&
!           & =- 0.5_wp*( u_c(jc,jk-1,jb) + u_c(jc,jk,jb))*w_c(jc,jk-1,jb)&
!           & + A_v(jc,jk,jb)*( u_c(jc,jk-1,jb) - u_c(jc,jk,jb) )/p_patch%patch_oce%del_zlev_i(jk-1)
! 
!           z_diff_v_i(jc,jk,jb)&
!           & =- 0.5_wp*( v_c(jc,jk-1,jb) + v_c(jc,jk,jb))*w_c(jc,jk-1,jb)&
!           & + A_v(jc,jk,jb)*( v_c(jc,jk-1,jb) - v_c(jc,jk,jb) )/p_patch%patch_oce%del_zlev_i(jk-1)
! 
! ! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb) *( v_c(jc,jk,jb) - v_c(jc,jk-1,jb) ),&
! ! &  v_c(jc,jk,jb), v_c(jc,jk-1,jb),p_patch%patch_oce%del_zlev_i(jk-1) 
! ! write(*,*)'vert adv:u: ',jk, jc,jb,w_c(jc,jk,jb) *( u_c(jc,jk,jb) - u_c(jc,jk-1,jb) ),&
! ! &  u_c(jc,jk,jb), u_c(jc,jk-1,jb)
!         ENDIF  ! jk-condition
!       ENDIF    ! at least 2 vertical layers
!     END DO
!   END DO
! !  write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
! !  & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
! END DO
! 
! !Including top&bottom boundary condition
! DO jb = i_startblk, i_endblk
!   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                      i_startidx, i_endidx, 1,min_rlcell)
!   DO jc = i_startidx, i_endidx
! 
!     z_diff_u_i(jc,1,jb) = top_bc_u_c(jc,jb)
!     z_diff_v_i(jc,1,jb) = top_bc_v_c(jc,jb)
! 
!     z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
! 
!     z_diff_u_i(jc,z_dolic,jb) = bot_bc_u_c(jc,jb)
!     z_diff_v_i(jc,z_dolic,jb) = bot_bc_v_c(jc,jb)
!   END DO
! END DO
! 
! ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! ! This mapping is the transposed of the vertical differencing.
! 
! DO jb = i_startblk, i_endblk
!   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!   DO jc = i_startidx, i_endidx
!     z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
!     IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!     DO jk = slev, z_dolic-1!elev-1
!         z_adv_diff_u_m(jc,jk,jb)                              &
!           & = (z_diff_u_i(jc,jk,jb) - z_diff_u_i(jc,jk+1,jb)) &
!           & / p_patch%patch_oce%del_zlev_m(jk) 
! 
!         z_adv_diff_v_m(jc,jk,jb)                              &
!           & = (z_diff_v_i(jc,jk,jb) - z_diff_v_i(jc,jk+1,jb)) &
!           & / p_patch%patch_oce%del_zlev_m(jk)
! 
!         CALL gvec2cvec( -z_adv_diff_u_m(jc,jk,jb), -z_adv_diff_v_m(jc,jk,jb), &
!           &             p_patch%cells%center(jc,jb)%lon,                    &
!           &             p_patch%cells%center(jc,jb)%lat,                    &
!           &             zu_cc(jc,jk,jb)%x(1),zu_cc(jc,jk,jb)%x(2),zu_cc(jc,jk,jb)%x(3))
! 
!     END DO
!       ENDIF
!   END DO
! !  write(*,*)'B max/min vert adv:',jk, maxval(z_adv_diff_u_m(:,jk,:)),&
! !  & minval(z_adv_diff_u_m(:,jk,:)),&
! !  & maxval(z_adv_diff_v_m(:,jk,:)),&
! !  & minval(z_adv_diff_v_m(:,jk,:))
! END DO
! 
! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers) 
!  CALL map_cell2edges( p_patch, zu_cc, veloc_adv_vert_e)
! 
! END subroutine veloc_adv_diff_vert_RBF2
! ! ! !-------------------------------------------------------------------------

END MODULE mo_oce_veloc_advection
