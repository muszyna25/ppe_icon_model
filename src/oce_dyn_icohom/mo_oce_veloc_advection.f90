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
                                & sea_boundary!,max_char_length
USE mo_model_domain,        ONLY: t_patch
USE mo_ocean_nml,           ONLY: n_zlev,iswm_oce, L_INVERSE_FLIP_FLOP!, idisc_scheme
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_oce_state,   ONLY: t_hydro_ocean_diag, t_hydro_ocean_aux
USE mo_oce_math_operators,  ONLY: rot_vertex_ocean, grad_fd_norm_oce
USE mo_scalar_product,      ONLY: map_cell2edges,  dual_flip_flop, &
  &                               primal_map_c2e, map_edges2edges
USE mo_interpolation,      ONLY: t_int_state, rbf_vec_interpol_edge,&
                                 &edges2cells_scalar, verts2edges_scalar
!USE mo_oce_index,              ONLY: c_k, ne_b, ne_i, form4ar! , ldbg

USE mo_math_utilities,     ONLY: t_cartesian_coordinates
IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


! PUBLIC :: veloc_adv_horz
! PUBLIC :: veloc_adv_vert

PUBLIC :: veloc_adv_horz_mimetic
PUBLIC :: veloc_adv_vert_mimetic
PUBLIC :: veloc_adv_horz_RBF
PUBLIC :: veloc_adv_vert_RBF


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
SUBROUTINE veloc_adv_horz_mimetic( p_patch, vn, p_diag, veloc_adv_horz_e, p_int)
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
INTEGER :: jk, jb, je, jv!,jc,ile,ibe, ie
!INTEGER :: i_startblk_c, i_endblk_c!, i_startidx_c, i_endidx_c
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
INTEGER :: rl_start_e, rl_end_e, rl_start_v, rl_end_v, rl_start_c, rl_end_c
!INTEGER ::  i_v1_idx, i_v1_blk, i_v2_idx, i_v2_blk
REAL(wp) :: z_e  (nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_u_c(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_v_c(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_vort_glb(nproma,n_zlev,p_patch%nblks_v)
REAL(wp) :: z_vort_flx(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_flx2(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_grad_ekin_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_RBF_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_RBF_c(nproma,n_zlev,p_patch%nblks_c)

!ARRAYS FOR TESTING
REAL(wp) :: z_vt(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_flx_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_veloc_adv_horz_e(nproma,n_zlev,p_patch%nblks_e)
LOGICAL, PARAMETER :: L_DEBUG = .TRUE. 
!-----------------------------------------------------------------------

! #slo# set local variable to zero due to nag -nan compiler-option
z_e             (:,:,:) = 0.0_wp
z_u_c           (:,:,:) = 0.0_wp
z_v_c           (:,:,:) = 0.0_wp
z_vort_glb      (:,:,:) = 0.0_wp
z_vort_flx      (:,:,:) = 0.0_wp
z_vort_flx2     (:,:,:) = 0.0_wp
veloc_adv_horz_e(:,:,:) = 0.0_wp
z_vt            (:,:,:) = 0.0_wp
z_vort_e        (:,:,:) = 0.0_wp
z_vort_flx_RBF  (:,:,:) = 0.0_wp
z_kin_RBF_c     (:,:,:) = 0.0_wp
z_grad_ekin_RBF (:,:,:) = 0.0_wp
rl_start_v = 1
rl_end_v = min_rlvert

rl_start_e = 1
rl_end_e  = min_rledge

rl_start_c = 1
rl_end_c  = min_rlcell

!i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
!i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)

slev = 1
elev = n_zlev


! calculate relative vorticity for all vertical layers
!  - tangential velocity (0) needed for no-slip boundary only
!CALL rot_vertex_ocean( vn, p_diag%vt, p_patch, p_diag%vort)
CALL rot_vertex_ocean(p_diag%ptp_vn, p_diag%vt, p_patch, p_diag%vort)

! add global to relative vorticity to obtain total vorticity
DO jb = i_startblk_v, i_endblk_v
  CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, &
                     i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
  DO jk = slev, elev
    DO jv = i_startidx_v, i_endidx_v
      z_vort_glb(jv,jk,jb) = p_diag%vort(jv,jk,jb) + p_patch%verts%f_v(jv,jb)
    ENDDO
  END DO
END DO

!calculate vorticity flux across dual edge
IF ( iswm_oce == 1 ) THEN
  z_vort_flx(:,:,:) = dual_flip_flop(p_patch, vn, z_vort_glb,p_diag%h_e)!p_diag%h_e
ELSEIF ( iswm_oce /= 1 ) THEN
  z_vort_flx(:,:,:) = dual_flip_flop(p_patch, vn, z_vort_glb, p_diag%h_e)
ENDIF
 DO jk = slev, elev
 write(*,*)'max/min vorticity:', jk,MAXVAL(p_diag%vort(:,jk,:)),MINVAL(p_diag%vort(:,jk,:)) 
 write(*,*)'max/min vort flux:', jk,MAXVAL(z_vort_flx(:,jk,:)),MINVAL(z_vort_flx(:,jk,:)) 
 END DO
! IF (L_DEBUG) THEN
!  DO jb = i_startblk_e, i_endblk_e
!  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
!                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e) 
!  DO je = i_startidx_e, i_endidx_e
!    write(*,*) 'vort flux:',je,jb,&!(jb-1)*nproma+je,jb,je, &
!      &               vn(je,1,jb),z_vort_flx(je,1,jb)
!  END DO
!  END DO
! ! ENDIF

!calculate gradient of kinetic energy
!The kinetic energy is already calculated at beginning
!of actual timestep in sbr "calc_scalar_product_for_veloc"
CALL grad_fd_norm_oce( p_diag%kin, &
                 & p_patch,    &
                 & p_diag%grad,&
                 & opt_slev=slev,opt_elev=elev )

IF(L_INVERSE_FLIP_FLOP)THEN
  z_e = p_diag%grad
  CALL map_edges2edges( p_patch,       &
                      & z_e,           &
                      & p_diag%grad,   &
                      & p_diag%thick_e )
ENDIF
!  write(*,*)'max/min kin energy:     ', MAXVAL(p_diag%kin(:,1,:)),MINVAL(p_diag%kin(:,1,:))
!  write(*,*)'max/min grad kin energy:', MAXVAL(p_diag%grad(:,1,:)),MINVAL(p_diag%grad(:,1,:)) 
! 

!---------------------------for testing and comparison with RBF--------------------------
!----------nonlinear coriolis and grad of kinetic energy computed with RBFs--------------
IF (L_DEBUG) THEN
 CALL rbf_vec_interpol_edge( vn, p_patch, p_int, p_diag%vt,&
                          & opt_slev=slev,opt_elev=elev)
 CALL rot_vertex_ocean(vn, p_diag%vt, p_patch, p_diag%vort,&
                          & opt_slev=slev,opt_elev=elev)
 CALL verts2edges_scalar( p_diag%vort, p_patch, p_int%v_1o2_e, &
                           z_vort_e, opt_slev=slev,opt_elev=elev, opt_rlstart=3)
!   DO jb = i_startblk_e,i_endblk_e
!     CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
!       &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
!     DO jk = slev, elev
!       DO je=i_startidx_e, i_endidx_e
!         i_v1_idx = p_patch%edges%vertex_idx(je,jb,1)
!         i_v1_blk = p_patch%edges%vertex_blk(je,jb,1)
!         i_v2_idx = p_patch%edges%vertex_idx(je,jb,2)
!         i_v2_blk = p_patch%edges%vertex_blk(je,jb,2)
!         z_vort_e(je,jk,jb) = 0.5_wp*(p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!                            &+        p_diag%vort(i_v2_idx,jk,i_v2_blk))
!       END DO
!     END DO
!   ENDDO

  DO jb = i_startblk_e,i_endblk_e
    CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
      &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
    DO jk = slev, elev
      DO je=i_startidx_e, i_endidx_e
        z_vort_flx_RBF(je,jk,jb) =&
      & -p_diag%vt(je,jk,jb)*(z_vort_e(je,jk,jb)+p_patch%edges%f_e(je,jb))

!write(*,*)'vort-flx:old:new:rbf',je,jb,&
!z_vort_flx(je,jk,jb), z_vort_flx_RBF(je,jk,jb)
      END DO
    END DO
  ENDDO
 !write(*,*)'max/min vort flux RBF:', MAXVAL(z_vort_flx_RBF(:,1,:)),MINVAL(z_vort_flx_RBF(:,1,:)) 
  DO jb = i_startblk_e,i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e,i_endblk_e, &
                       i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e
        ! calculate kinetic energy at edges from normal and tangential comp.
        z_kin_RBF_e(je,jk,jb) =0.5_wp*(p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb)+&
                          &               vn(je,jk,jb)*vn(je,jk,jb) )
      ENDDO
    ENDDO
  ENDDO

!!$OMP END DO
!!$OMP END PARALLEL
  ! Bilinear interpolation of kinetic energy from the edges to the cells
   CALL edges2cells_scalar( z_kin_RBF_e, p_patch, p_int%e_bln_c_s,  &
     &                      z_kin_RBF_c, opt_slev=slev, opt_elev=elev)

   CALL grad_fd_norm_oce( z_kin_RBF_c, p_patch, z_grad_ekin_RBF,opt_slev=slev,opt_elev=elev )
END IF
!--------------END OF TESTING----------------------------------------------------------

!Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e
      z_veloc_adv_horz_e(je,jk,jb)  = -z_vort_flx_RBF(je,jk,jb) - z_grad_ekin_RBF(je,jk,jb) 
!
!      veloc_adv_horz_e(je,jk,jb)    = -z_vort_flx_RBF(je,jk,jb) - p_diag%grad(je,jk,jb)
      veloc_adv_horz_e(je,jk,jb)    = z_vort_flx(je,jk,jb)- p_diag%grad(je,jk,jb)!
!       write(*,*)'horz adv:vort-flx:vort-flx-RBF',je,jk,jb,z_vort_flx(je,1,jb), &
!         &        -z_vort_flx_RBF(je,jk,jb),p_diag%grad(je,jk,jb),z_grad_ekin_RBF(je,jk,jb)
!      IF( (z_vort_flx(je,jk,jb)>0.0.AND.-z_vort_flx_RBF(je,jk,jb)<0)&
! &   .OR. (z_vort_flx(je,jk,jb)<0.0.AND.-z_vort_flx_RBF(je,jk,jb)>0))THEN
!           Write(*,*)'wraning:differing signs',je,jk,jb,z_vort_flx(je,1,jb), &
!        &        -z_vort_flx_RBF(je,jk,jb) 
! ENDIF
!      write(*,*)'horz adv:grad ekin:vort-flx:',je,jk,jb, vn(je,1,jb),z_vort_flx(je,1,jb), &
!        &        veloc_adv_horz_e(je,1,jb) !, p_diag%grad(je,1,jb)
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
INTEGER :: jk, jb, je
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e!, rl_start_v, rl_end_v, rl_start_c, rl_end_c
!INTEGER :: i_v1_idx, i_v1_blk, i_v2_idx, i_v2_blk 
REAL(wp) :: z_e  (nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_glb(nproma,n_zlev,p_patch%nblks_v)
!REAL(wp) :: z_vort_flx(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_grad_ekin_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_RBF_e(nproma,n_zlev,p_patch%nblks_e)

!ARRAYS FOR TESTING
REAL(wp) :: z_vort_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_vort_flx_RBF(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_kin_e_RBF(nproma,n_zlev,p_patch%nblks_e)
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
rl_end_e  = min_rledge

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

slev = 1
elev = n_zlev

  CALL rbf_vec_interpol_edge( vn,       &
                            & p_patch,  &
                            & p_int,    &
                            & p_diag%vt,&
                            & opt_slev=slev,opt_elev=elev)

  CALL rot_vertex_ocean(vn, p_diag%vt, p_patch, p_diag%vort)
  CALL verts2edges_scalar( p_diag%vort, p_patch, p_int%v_1o2_e, &
                           z_vort_e, opt_slev=slev,opt_elev=elev, opt_rlstart=3)

!   DO jb = i_startblk_e,i_endblk_e
!     CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
!       &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
!     DO jk = slev, elev
!       DO je=i_startidx_e, i_endidx_e
!         i_v1_idx = p_patch%edges%vertex_idx(je,jb,1)
!         i_v1_blk = p_patch%edges%vertex_blk(je,jb,1)
!         i_v2_idx = p_patch%edges%vertex_idx(je,jb,2)
!         i_v2_blk = p_patch%edges%vertex_blk(je,jb,2)
!         z_vort_e(je,jk,jb) =&
!         & 0.5_wp*(p_diag%vort( i_v1_idx,jk,i_v1_blk)&
!         &+        p_diag%vort(i_v1_idx,jk,i_v1_blk))
!       END DO
!     END DO
!   ENDDO

  DO jb = i_startblk_e,i_endblk_e
    CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, &
      &                              i_startidx_e, i_endidx_e, rl_start_e,rl_end_e)
    DO jk = slev, elev
      DO je=i_startidx_e, i_endidx_e
        z_vort_flx_RBF(je,jk,jb) =&
         & -p_diag%vt(je,jk,jb)*(z_vort_e(je,jk,jb)+p_patch%edges%f_e(je,jb))
!        & -p_diag%vt(je,jk,jb)*p_patch%edges%f_e(je,jb)
      END DO
    END DO
  ENDDO

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
!!$OMP END DO
!!$OMP END PARALLEL
  ! Bilinear interpolation of kinetic energy from the edges to the cells
   CALL edges2cells_scalar( z_kin_RBF_e,     &
                          & p_patch,         &
                          & p_int%e_bln_c_s, &
                          & p_diag%kin,      &
                          & opt_slev=slev,opt_elev=elev)

   CALL grad_fd_norm_oce( p_diag%kin, &
                        & p_patch,    &
                        & z_grad_ekin_RBF, opt_slev=slev,opt_elev=elev)


!Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e
      veloc_adv_horz_e(je,jk,jb) =&
      & -z_vort_flx_RBF(je,jk,jb) - z_grad_ekin_RBF(je,jk,jb) 
   END DO
  END DO
END DO

 DO jk = slev, elev
 write(*,*)'max/min vort flux:kin energy: advection', jk,&
& MAXVAL(z_vort_flx_RBF(:,jk,:)),MINVAL(z_vort_flx_RBF(:,jk,:)),&
& MAXVAL(z_grad_ekin_RBF(:,jk,:)),MINVAL(z_grad_ekin_RBF(:,jk,:)),&
& MAXVAL(veloc_adv_horz_e(:,jk,:)),MINVAL(veloc_adv_horz_e(:,jk,:))
 END DO


!   do jb=1,10!i_startblk_e, i_endblk_e
!    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
!                      i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
!   do je= i_startidx_e, i_endidx_e
!    write(*,*)'RBF: vn_in, grad_ekin, vort_flx:',je,jb, vn(je,1,jb), &
!       &       z_grad_ekin_RBF(je,1,jb), &
!       &       -z_vort_flx_RBF(je,1,jb), veloc_adv_horz_e(je,1,jb)
!   enddo
!   enddo
!stop
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
DO jk = slev, elev
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      !check if we are on land: To be replaced by 3D lsm
      IF (p_patch%patch_oce%dolic_c(jc,jb) /= 0) THEN
        !1a) 0cean surface
        IF(jk==slev)THEN

          z_adv_u_i(jc,jk,jb)%x = top_bc_w_c(jc,jb)*p_aux%bc_top_veloc_cc(jc,jb)%x

        !1b) ocean bottom 
        ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN

          z_adv_u_i(jc,jk+1,jb)%x = bot_bc_w_c(jc,jb)*p_aux%bc_bot_veloc_cc(jc,jb)%x

        !1c) ocean interior 
        ELSEIF( jk >slev .AND. jk < p_patch%patch_oce%dolic_c(jc,jb))THEN

          z_adv_u_i(jc,jk,jb)%x&
            & = p_diag%w(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x)&
            & / p_patch%patch_oce%del_zlev_i(jk-1)
! ! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb),&
!&( p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x )

        ENDIF
      ENDIF
    END DO
  END DO
! write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
! & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
END DO

! ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
! ! This mapping is the transposed of the vertical differencing.
! 
! !1) From surface down to one layer before bottom
DO jk = slev, elev-1
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx
      ! #slo# 2011-04-01 - dolic_c is now correct
      ! compare also formulation in veloc_adv_vert_RBF l. 765-790
      IF ( jk <= p_patch%patch_oce%dolic_c(jc,jb) ) THEN
      !IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        z_adv_u_m(jc,jk,jb)%x &
        & = (p_patch%patch_oce%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)%x&
        & +  p_patch%patch_oce%del_zlev_i(jk)*z_adv_u_i(jc,jk,jb)%x) &
        & / (2.0_wp*p_patch%patch_oce%del_zlev_m(jk))
      ! & / (p_patch%patch_oce%del_zlev_m(jk+1)+p_patch%patch_oce%del_zlev_m(jk))
      ENDIF
    END DO
  END DO
! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
END DO

! !Bottom layer
! !The value of p_patch%patch_oce%del_zlev_i at the botom is 0.5*p_patch%patch_oce%del_zlev_m
! !The dimensioning of the first arrays requires to seperate the vertical loop.
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    ! #slo# 2011-04-01 - dolic_c is now correct
    IF ( elev <= p_patch%patch_oce%dolic_c(jc,jb) ) THEN
    !IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

      z_adv_u_m(jc,elev,jb)%x &
      ! #slo# 2011-04-01 - Attention to multiplying by dz(elev) in both terms
      & = (0.5_wp*p_patch%patch_oce%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,jb)%x&
      & +  p_patch%patch_oce%del_zlev_i(elev)*z_adv_u_i(jc,elev,jb)%x) &
      & / (2.0_wp*p_patch%patch_oce%del_zlev_m(elev))
      ENDIF
    END DO
  END DO

! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
  CALL map_cell2edges( p_patch, z_adv_u_m, veloc_adv_vert_e )

 DO jk=1,n_zlev
   WRITE(*,*) 'max/min vert adv FINAL',jk, &
     &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
 END DO

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
          ! u-component
          z_adv_u_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_u_c(jc,jb)
          ! v-component
          z_adv_v_i(jc,jk,jb) =  top_bc_w_c(jc,jb)*top_bc_v_c(jc,jb)
!write(*,*)'vert adv: top:',jc,jb,jk,top_bc_w_c(jc,jb),top_bc_u_c(jc,jb),top_bc_v_c(jc,jb) 
        !1b) ocean bottom 
        ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
          ! u-component
          z_adv_u_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_u_c(jc,jb)
          ! v-component
          z_adv_v_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_v_c(jc,jb)

        !1c) ocean interior 
        ELSEIF( jk>slev .AND.  jk < p_patch%patch_oce%dolic_c(jc,jb))THEN
          ! u-component
          z_adv_u_i(jc,jk,jb)&
          & = w_c(jc,jk,jb) *( u_c(jc,jk-1,jb) - u_c(jc,jk,jb) )&
            & / p_patch%patch_oce%del_zlev_i(jk-1)
          ! v-component
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
 write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
 & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
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
DO jk=1,n_zlev
  WRITE(*,*) 'max/min vert adv FINAL',jk, &
    &        MAXVAL(veloc_adv_vert_e(:,jk,:)), MINVAL(veloc_adv_vert_e(:,jk,:))
END DO

END subroutine veloc_adv_vert_RBF
! ! !-------------------------------------------------------------------------


END MODULE mo_oce_veloc_advection
