!>
!! Contains the implementation of velocity advection in vector invariant form
!! that is used in the ocean model.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2010)
!!  Modified by Stephan Lorenz,     MPI-M (2010-11)
!!   - implementation of new PtP reconstruction
!!   mpi parallelized LL
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
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: sync_e, sync_c, sync_v, sync_patch_array
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_impl_constants,      ONLY: sea_boundary, sea, boundary, min_dolic!, &
  !  &                               min_rlcell, min_rledge, min_rlvert
  USE mo_ocean_nml,           ONLY: n_zlev!, iswm_oce, l_inverse_flip_flop, ab_beta, ab_gam
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_oce_state,           ONLY: t_hydro_ocean_diag!, v_base
  USE mo_oce_math_operators,  ONLY: rot_vertex_ocean_rbf, &
    &                               grad_fd_norm_oce_3d, &!grad_fd_norm_oce, div_oce_3d, &
    &                               rot_vertex_ocean_3d!, rot_vertex_ocean
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, vector_product!,cc2gc, gvec2cvec, gc2cc
  USE mo_scalar_product,      ONLY: map_cell2edges_3D, nonlinear_coriolis_3D,nonlinear_coriolis_3d_2!, map_edges2edges, map_edges2cell
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=12)  :: str_module = 'oceVelocAdv '  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
! CHARACTER(len=12)  :: str_module = '__FILE__'

  PUBLIC  :: veloc_adv_horz_mimetic
  PRIVATE :: veloc_adv_horz_mimetic_div
  PRIVATE :: veloc_adv_horz_mimetic_rot
  PUBLIC  :: veloc_adv_vert_mimetic
  PRIVATE :: veloc_adv_vert_mimetic_div
  PRIVATE :: veloc_adv_vert_mimetic_rot
  PRIVATE :: veloc_adv_vert_mimetic_rot_flux
  !PUBLIC  :: veloc_adv_horz_rbf
  !PUBLIC  :: veloc_adv_vert_rbf

  INTEGER, PARAMETER, PRIVATE :: rotational_form = 0
  INTEGER, PARAMETER, PRIVATE :: divergence_form = 1
  INTEGER, PARAMETER, PRIVATE :: velocity_advection_form=0

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Computes horizontal advection of a (edge based) vector field.
  !! either by using rotational/vector-invariant form of velocity advection
  !! or the divergence form
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!  
  !!   mpi parallelized LL
  SUBROUTINE veloc_adv_horz_mimetic( p_patch_3D,        &
    & vn_old,          &
    & vn_new,          &
    & p_diag,          &
    & veloc_adv_horz_e,&
    & p_op_coeff)
    !
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(inout)          :: vn_old(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(inout)          :: vn_new(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_hydro_ocean_diag)         :: p_diag
    REAL(wp), INTENT(out)            :: veloc_adv_horz_e(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff


    !Interpolation necessary just for testing
    !TYPE(t_int_state),TARGET,INTENT(in), OPTIONAL :: p_int
    !-----------------------------------------------------------------------

    IF (velocity_advection_form == rotational_form) THEN
      CALL veloc_adv_horz_mimetic_rot( p_patch_3D,        &
        & vn_old,          &
        & vn_new,          &
        & p_diag,          &
        & veloc_adv_horz_e,&
        & p_op_coeff)!,p_int)
    ELSEIF (velocity_advection_form == divergence_form) THEN
      CALL veloc_adv_horz_mimetic_div( p_patch_3D,      &
        & vn_old,        &
        & p_diag,        &
        & p_op_coeff,    &
        & veloc_adv_horz_e)
    ENDIF

  END SUBROUTINE veloc_adv_horz_mimetic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) vector field.
  !! either by using rotational/vector-invariant form of velocity advection
  !! or the divergence form
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!
  !!   mpi parallelized LL
  SUBROUTINE veloc_adv_vert_mimetic( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    !-----------------------------------------------------------------------

    IF (velocity_advection_form == rotational_form) THEN
      CALL veloc_adv_vert_mimetic_rot( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)!veloc_adv_vert_mim_rot_flux2
      !CALL veloc_adv_vert_mimetic_rot_flux2(p_diag,p_op_coeff, veloc_adv_vert_e)
      !CALL veloc_adv_vert_mimetic_rot_flux( p_diag,p_op_coeff, veloc_adv_vert_e)
    ELSEIF (velocity_advection_form == divergence_form) THEN
      CALL veloc_adv_vert_mimetic_div( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)
    ENDIF

  END SUBROUTINE veloc_adv_vert_mimetic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
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
  !!  mpi parallelized LL
  !!
  SUBROUTINE veloc_adv_horz_mimetic_rot( p_patch_3D,     &
    & vn_old,          &
    & vn_new,          &
    & p_diag,          &
    & veloc_adv_horz_e,&
    & p_op_coeff)
    !
    !
    !  patch on which computation is performed
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    !
    ! normal velocity  of which advection is computed
    REAL(wp), INTENT(inout) :: vn_old(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(inout) :: vn_new(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    !
    !diagnostic ocean state stores horizontally advected velocity
    TYPE(t_hydro_ocean_diag) :: p_diag

    ! variable in which horizontally advected velocity is stored
    REAL(wp), INTENT(out) :: veloc_adv_horz_e(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    !
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    !Interpolation necessary just for testing
    !TYPE(t_int_state),TARGET,INTENT(in), OPTIONAL :: p_int


    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jk, jb, jc, je!, jv, ile, ibe, ie, jev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_e, i_endidx_e
    REAL(wp) :: z_e            (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_vort_flx     (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_grad_ekin_rbf(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_kin_rbf_e    (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_kin_rbf_c    (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_vort_e       (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_vort_flx_rbf (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3
    REAL(wp) :: z_weight_e1, z_weight_e2, z_weight_e3
    LOGICAL, PARAMETER :: l_debug = .FALSE.
    !LOGICAL, PARAMETER :: L_ENSTROPHY_DISSIPATION=.FALSE.
    TYPE(t_subset_range), POINTER :: all_edges, all_cells
    TYPE(t_patch), POINTER         :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_edges => p_patch%edges%all
    all_cells => p_patch%cells%all
    !-----------------------------------------------------------------------

    z_e             (1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
    z_vort_flx      (1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
    veloc_adv_horz_e(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp

    slev = 1
    elev = n_zlev

    !calculate vorticity flux across dual edge
    !   LL: nonlinear_coriolis_3d is mpi parallized. p_diag%vort must have been synced
    CALL nonlinear_coriolis_3d( p_patch_3D, &
      & vn_old,          &
      & p_diag%p_vn_dual,&
      & p_diag%thick_e,  &
      & p_diag%vort,     &
      & p_op_coeff,      &
      & z_vort_flx)

!     CALL nonlinear_coriolis_3d_2( p_patch_3D, &
!       & vn_old,          &
!       & p_diag%p_vn_dual,&
!       & p_diag%thick_e,  &
!       & p_diag%vort,     &
!       & p_op_coeff,      &
!       & z_vort_flx)
    !-------------------------------------------------------------------------------
    ! IF(L_ENSTROPHY_DISSIPATION)THEN
    !  DO jk = slev, elev
    !   write(*,*)'max/min vort flux: ',MAXVAL(z_vort_flx(:,jk,:)),&
    !                                         &MINVAL(z_vort_flx(:,jk,:))
    !   END DO
    ! z_vort_flx=laplacian4vortex_flux(p_patch,z_vort_flx)
    ! ENDIF
    !-------------------------------------------------------------------------------

    !calculate gradient of kinetic energy
    CALL grad_fd_norm_oce_3d( p_diag%kin, &
      & p_patch_3D,           &
      & p_op_coeff%grad_coeff,&
      & p_diag%grad)

    CALL sync_patch_array(sync_e, p_patch, p_diag%grad(:,:,:))    

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimRot: kin energy'        ,p_diag%kin              ,str_module,idt_src)
    CALL dbg_print('HorzMimRot: vorticity'         ,p_diag%vort             ,str_module,idt_src)
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimRot: vorticity flux'    ,z_vort_flx              ,str_module,idt_src)
    !CALL dbg_print('HorzMimRot: p_vn%x(1)'         ,p_diag%p_vn%x(1)        ,str_module,idt_src)
    !CALL dbg_print('HorzMimRot: p_vn_dual%x(1)'    ,p_diag%p_vn%x(1)        ,str_module,idt_src)
    CALL dbg_print('HorzMimRot: grad kin en'       ,p_diag%grad             ,str_module,idt_src)
    !---------------------------------------------------------------------

    ! IF(L_INVERSE_FLIP_FLOP)THEN
    !   CALL map_edges2edges( p_patch,    &
    !                       & p_diag%grad,&
    !                       & z_e,        &
    !                       & p_diag%thick_e)
    !   p_diag%grad = z_e
    !   CALL sync_patch_array(sync_e, p_patch, p_diag%grad(:,jk,:))
    ! ENDIF

    !---------------------------for testing and comparison with RBF--------------------------
    !----------nonlinear coriolis and grad of kinetic energy computed with RBFs--------------
    !----------needs interpolation state

! !     IF (l_debug) THEN
! !       z_vort_flx_rbf  (:,:,:) = 0.0_wp
! !       z_kin_rbf_c     (:,:,:) = 0.0_wp
! !       z_grad_ekin_rbf (:,:,:) = 0.0_wp
! !       z_vort_e        (:,:,:) = 0.0_wp
! !       CALL rbf_vec_interpol_edge( vn_old,       &
! !         & p_patch,  &
! !         & p_int,    &
! !         & p_diag%vt,&
! !         & opt_slev=slev,opt_elev=elev)
! !       CALL sync_patch_array(SYNC_E, p_patch, p_diag%vt)
! ! 
! ! 
! !       CALL rot_vertex_ocean_rbf(p_patch, vn_old, p_diag%vt, p_diag%vort)
! !       CALL sync_patch_array(SYNC_V, p_patch, p_diag%vort)
! ! 
! !       CALL verts2edges_scalar( p_diag%vort, p_patch, p_int%v_1o2_e, &
! !         & z_vort_e, opt_slev=slev,opt_elev=elev, opt_rlstart=3)
! !       CALL sync_patch_array(SYNC_E, p_patch, z_vort_e)
! ! 
! !       DO jb = all_edges%start_block, all_edges%end_block
! !         CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
! !         DO jk = slev, elev
! !           DO je=i_startidx_e, i_endidx_e
! ! 
! !             IF ( v_base%lsm_e(je,jk,jb) == sea ) THEN  !  #slo# <= sea_boundary?
! !               z_vort_flx_rbf(je,jk,jb) = &
! !                 & p_diag%vt(je,jk,jb)*(z_vort_e(je,jk,jb)+ p_patch%edges%f_e(je,jb))
! !             ENDIF
! !           END DO
! !         END DO
! !       ENDDO
! ! 
! !       CALL rbf_vec_interpol_cell( vn_old, p_patch, p_int, p_diag%u,  &
! !         & p_diag%v, opt_slev=slev, opt_elev=elev)
! !       CALL sync_patch_array(SYNC_C, p_patch, p_diag%v)
! ! 
! !       DO jb = all_edges%start_block, all_edges%end_block
! !         CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
! !         DO jk = slev, elev
! !           DO je = i_startidx_e, i_endidx_e
! !             ! calculate kinetic energy at edges from normal and tangential comp.
! !             z_kin_rbf_e(je,jk,jb) =0.5_wp*(p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb)&
! !               & +    vn_old(je,jk,jb)*vn_old(je,jk,jb) )
! !           ENDDO
! !         ENDDO
! !       ENDDO
! ! 
! !       DO jb = all_cells%start_block, all_cells%end_block
! !         CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !         DO jk = slev, elev
! !           DO jc = i_startidx_c, i_endidx_c
! !             IF ( v_base%lsm_c(jc,jk,jb) > sea_boundary ) THEN
! !               z_kin_rbf_c(jc,jk,jb) = 0.0_wp
! !             ELSE
! ! 
! !               ile1 = p_patch%cells%edge_idx(jc,jb,1)
! !               ibe1 = p_patch%cells%edge_blk(jc,jb,1)
! !               ile2 = p_patch%cells%edge_idx(jc,jb,2)
! !               ibe2 = p_patch%cells%edge_blk(jc,jb,2)
! !               ile3 = p_patch%cells%edge_idx(jc,jb,3)
! !               ibe3 = p_patch%cells%edge_blk(jc,jb,3)
! !               z_weight_e1 = 0.0_wp
! !               z_weight_e2 = 0.0_wp
! !               z_weight_e3 = 0.0_wp
! !               IF(v_base%lsm_e(ile1,jk,ibe1)<= boundary)THEN
! !                 z_weight_e1 = p_patch%edges%area_edge(ile1,ibe1)
! !               ENDIF
! !               IF(v_base%lsm_e(ile2,jk,ibe2)<= boundary)THEN
! !                 z_weight_e2 = p_patch%edges%area_edge(ile2,ibe2)
! !               ENDIF
! !               IF(v_base%lsm_e(ile3,jk,ibe3)<= boundary)THEN
! !                 z_weight_e3 = p_patch%edges%area_edge(ile3,ibe3)
! !               ENDIF
! ! 
! !               !write(*,*)'weights',jc,jk,jb,z_weight_e1,z_weight_e2,z_weight_e3
! !               z_kin_rbf_c(jc,jk,jb) = (z_kin_rbf_e(ile1,jk,ibe1)*z_weight_e1&
! !                 & + z_kin_rbf_e(ile2,jk,ibe2)*z_weight_e2&
! !                 & + z_kin_rbf_e(ile3,jk,ibe3)*z_weight_e3)&
! !                 & /(z_weight_e1+z_weight_e2+z_weight_e3)
! !             ENDIF
! !           END DO
! !         END DO
! !       END DO
! ! 
! !       CALL grad_fd_norm_oce_3d( p_diag%kin, &
! !         & p_patch,    &
! !         & p_op_coeff%grad_coeff,&
! !         & z_grad_ekin_rbf)
! ! 
! !       CALL sync_patch_array(SYNC_E, p_patch, z_grad_ekin_rbf)      
! ! 
! !       !---------Debug Diagnostics-------------------------------------------
! !       idt_src=4  ! output print level (1-5, fix)
! !       CALL dbg_print('L_DEBUG: vort flux RBF'     ,z_vort_flx_rbf           ,str_module,idt_src)
! !       !CALL dbg_print('L_DEBUG: p_vn%x(1)'         ,p_diag%p_vn%x(1)         ,str_module,idt_src)
! !       !CALL dbg_print('L_DEBUG: p_vn_dual%x(1)'    ,p_diag%p_vn%x(1)         ,str_module,idt_src)
! !       !---------------------------------------------------------------------
! ! 
! !     END IF ! L_DEBUG
    !--------------END OF TESTING----------------------------------------------------------

    !Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = slev, elev
        DO je = i_startidx_e, i_endidx_e
         !IF(v_base%lsm_e(je,jk,jb)<= boundary)THEN
         IF(p_patch_3D%lsm_e(je,jk,jb)<= boundary)THEN
          !veloc_adv_horz_e(je,jk,jb)  = z_vort_flx_RBF(je,jk,jb) + z_grad_ekin_RBF(je,jk,jb)
          !veloc_adv_horz_e(je,jk,jb)= z_vort_flx(je,jk,jb) + z_grad_ekin_RBF(je,jk,jb)
          !veloc_adv_horz_e(je,jk,jb)    = z_vort_flx_RBF(je,jk,jb)+ p_diag%grad(je,jk,jb)
          veloc_adv_horz_e(je,jk,jb) = (z_vort_flx(je,jk,jb) + p_diag%grad(je,jk,jb))
          !        write(*,*)'horz adv:vort-flx:vort-flx-RBF',je,jk,jb,z_vort_flx(je,1,jb), &
          !          &        z_vort_flx_RBF(je,jk,jb)!,&!p_diag%grad(je,jk,jb),z_grad_ekin_RBF(je,jk,jb),&
          !&        veloc_adv_horz_e(je,jk,jb)!, z_veloc_adv_horz_e(je,jk,jb)
          !       IF( (z_vort_flx(je,jk,jb)>0.0.AND.z_vort_flx_RBF(je,jk,jb)<0)&
          !       &   .OR. (z_vort_flx(je,jk,jb)<0.0.AND.z_vort_flx_RBF(je,jk,jb)>0))THEN
          !             Write(*,*)'wraning:differing signs',je,jk,jb,z_vort_flx(je,1,jb), &
          !          &        z_vort_flx_RBF(je,jk,jb)
          !        ENDIF
          !      write(*,*)'horz adv:grad ekin:vort-flx:',je,jk,jb, vn(je,1,jb),z_vort_flx(je,1,jb), &
          !        &        veloc_adv_horz_e(je,1,jb) !, p_diag%grad(je,1,jb)
          ! IF ( v_base%lsm_e(je,jk,jb) == boundary ) THEN
          ! write(*,*)'vt ',je,jk,jb, p_diag%vt(je,jk,jb)
          ENDIF
        END DO
      END DO
    END DO

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimRot: final Vel.Adv.'    ,veloc_adv_horz_e        ,str_module,idt_src)
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_horz_mimetic_rot
  !-------------------------------------------------------------------------  

  !-------------------------------------------------------------------------
  !!  mpi parallelized LL
  SUBROUTINE veloc_adv_horz_mimetic_div( p_patch_3D,   &
    & vn,             &
    & p_diag,         &
    & p_op_coeff,     &
    & veloc_adv_horz_e)
    !

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(inout)          :: vn(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e) ! dim: (nproma,n_zlev,nblks_e)
    TYPE(t_hydro_ocean_diag)         :: p_diag
    TYPE(t_operator_coeff),INTENT(in):: p_op_coeff
    REAL(wp), INTENT(out)            :: veloc_adv_horz_e(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
    !
    !Local variables
    !
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jk, jb, je, jc
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e1, ib_e1, il_e2, ib_e2, il_e3, ib_e3

    REAL(wp)                      :: z_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: u_v_cc_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: u_v_cc_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_div_vec_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER         :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)

    ! #slo# set local variable to zero due to nag -nan compiler-option
    all_edges => p_patch%edges%all
    all_cells => p_patch%cells%all

    z_e             (1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
    veloc_adv_horz_e(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp

    slev = 1
    elev = n_zlev

    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = slev, elev
        DO je = i_startidx_e, i_endidx_e
          !IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
          IF(p_patch_3D%lsm_e(je,jk,jb)<= boundary)THEN
            !Neighbouring cells
            il_c1 = p_patch%edges%vertex_idx(je,jb,1)
            ib_c1 = p_patch%edges%vertex_blk(je,jb,1)
            il_c2 = p_patch%edges%vertex_idx(je,jb,2)
            ib_c2 = p_patch%edges%vertex_blk(je,jb,2)

            !velocity vector at edges
            u_v_cc_e(je,jk,jb)%x = 0.5_wp * (p_diag%p_vn_dual(il_c1,jk,ib_c1)%x + &
              & p_diag%p_vn_dual(il_c2,jk,ib_c2)%x)

            u_v_cc_e(je,jk,jb)%x = vn(je,jk,jb) * u_v_cc_e(je,jk,jb)%x
          ENDIF
        END DO
      END DO
    END DO


   DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jk = slev, elev
        DO jc = i_startidx_c, i_endidx_c

!           cc_c = gc2cc(p_patch%cells%center(jc,jb))p_op_coeff%cell_position_cc
          !cc_c = p_op_coeff%cell_position_cc(jc,jk,jb)
          u_v_cc_c(jc,jk,jb)= vector_product(p_op_coeff%cell_position_cc(jc,jk,jb),&
                                &p_diag%p_vn(jc,jk,jb))
          u_v_cc_c(jc,jk,jb)%x = p_patch%cells%f_c(jc,jb)*u_v_cc_c(jc,jk,jb)%x

         il_e1 = p_patch%cells%edge_idx(jc,jb,1)
         ib_e1 = p_patch%cells%edge_blk(jc,jb,1)

         il_e2 = p_patch%cells%edge_idx(jc,jb,2)
         ib_e2 = p_patch%cells%edge_blk(jc,jb,2)

         il_e3 = p_patch%cells%edge_idx(jc,jb,3)
         ib_e3 = p_patch%cells%edge_blk(jc,jb,3)

         z_div_vec_c(jc,jk,jb)%x =  &
              & u_v_cc_e(il_e1,jk,ib_e1)%x * p_op_coeff%div_coeff(jc,jk,jb,1) + &
              & u_v_cc_e(il_e2,jk,ib_e2)%x * p_op_coeff%div_coeff(jc,jk,jb,2) + &
              & u_v_cc_e(il_e3,jk,ib_e3)%x * p_op_coeff%div_coeff(jc,jk,jb,3) +&
             & u_v_cc_c(jc,jk,jb)%x

        END DO
      END DO
   END DO

!     CALL map_cell2edges( p_patch, z_div_vec_c, veloc_adv_horz_e, &
!       & opt_slev=slev, opt_elev=elev )
    CALL map_cell2edges_3D( p_patch_3D, z_div_vec_c, veloc_adv_horz_e,p_op_coeff)


    !calculates the curl. This is needed in Laplace-beltrami operator (velocity diffusion).
    !It is not needed for velocity advection.
      CALL rot_vertex_ocean_3D( p_patch_3D,          &
                              & vn,                  &
                              & p_diag%p_vn_dual,    &
                              & p_op_coeff,          &
                              & p_diag%vort)
    CALL sync_patch_array(SYNC_V, p_patch, p_diag%vort)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimDiv: final Vel.Adv.'    ,veloc_adv_horz_e        ,str_module,idt_src)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimDiv: vorticity'         ,p_diag%vort             ,str_module,idt_src)
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_horz_mimetic_div
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! ! FUNCTION laplacian4vortex_flux(p_patch, vort_flux_in) RESULT(vort_flux_out)
  ! ! TYPE(t_patch),TARGET, INTENT(in) :: p_patch
  ! ! REAL(wp), INTENT(inout) :: vort_flux_in(:,:,:)
  ! ! REAL(wp) ::vort_flux_out(SIZE(vort_flux_in,1), SIZE(vort_flux_in,2), SIZE(vort_flux_in,3))
  ! !
  ! ! !
  ! ! !Local Variables
  ! ! !
  ! ! INTEGER :: slev, elev
  ! ! INTEGER :: jk, jb, jc, je
  ! ! INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  ! ! INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  ! ! INTEGER :: rl_start_e, rl_end_e,rl_start_c, rl_end_c
  ! ! REAL(wp) :: z_tmp(nproma,n_zlev,p_patch%nblks_e)
  ! ! TYPE(t_cartesian_coordinates)    :: z_pv_cc(nproma,n_zlev,p_patch%nblks_c)
  ! ! INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  ! ! INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  ! ! TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,n_zlev,p_patch%nblks_e)
  ! ! TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,p_patch%nblks_c)
  ! ! !-----------------------------------------------------------------------
  ! ! rl_start_e = 1
  ! ! rl_end_e   = min_rledge
  ! ! rl_start_c = 1
  ! ! rl_end_c   = min_rlcell
  ! !
  ! ! i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
  ! ! i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
  ! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  ! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  ! !
  ! ! slev = 1
  ! ! elev = n_zlev
  ! !
  ! ! z_tmp = vort_flux_in
  ! !
  ! ! CALL map_edges2cell( p_patch, z_tmp, z_pv_cc)
  ! ! DO jb = i_startblk_c, i_endblk_c
  ! !
  ! !   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
  ! !                      i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
  ! !     DO jk = slev, elev
  ! !       DO jc = i_startidx_c, i_endidx_c
  ! !           z_div_grad_u(jc,jk,jb)%x =  0.0_wp
  ! !       END DO
  ! !     END DO
  ! !   END DO
  ! ! DO jb = i_startblk_e, i_endblk_e
  ! !
  ! !   CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e,&
  ! !                    &  i_startidx_e, i_endidx_e,&
  ! !                    &  rl_start_e, rl_end_e)
  ! !   DO jk = slev, elev
  ! !     DO je = i_startidx_e, i_endidx_e
  ! !       z_grad_u(je,jk,jb)%x = 0.0_wp
  ! !     ENDDO
  ! !   END DO
  ! ! END DO
  ! !
  ! ! !Step 1: Calculate gradient of cell velocity vector.
  ! ! !Result is a gradient vector, located at edges
  ! ! !Step 2: Multiply each component of gradient vector with mixing coefficients
  ! ! DO jb = i_startblk_e, i_endblk_e
  ! !
  ! !   CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e,&
  ! !                    &  i_startidx_e, i_endidx_e,&
  ! !                    &  rl_start_e, rl_end_e)
  ! !   DO jk = slev, elev
  ! !     DO je = i_startidx_e, i_endidx_e
  ! !
  ! !       !Get indices of two adjacent triangles
  ! !       il_c1 = p_patch%edges%cell_idx(je,jb,1)
  ! !       ib_c1 = p_patch%edges%cell_blk(je,jb,1)
  ! !       il_c2 = p_patch%edges%cell_idx(je,jb,2)
  ! !       ib_c2 = p_patch%edges%cell_blk(je,jb,2)
  ! !
  ! !     IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
  ! !       z_grad_u(je,jk,jb)%x = &
  ! !         &                  (z_pv_cc(il_c2,jk,ib_c2)%x &
  ! !         &                  - z_pv_cc(il_c1,jk,ib_c1)%x)              &
  ! !         &                  / p_patch%edges%dual_edge_length(je,jb)
  ! !     ELSE
  ! !       z_grad_u(je,jk,jb)%x = 0.0_wp
  ! !     ENDIF
  ! !     ENDDO
  ! !   END DO
  ! ! END DO
  ! !
  ! ! !Step 2: Apply divergence to each component of mixing times gradient vector
  ! ! iidx => p_patch%cells%edge_idx
  ! ! iblk => p_patch%cells%edge_blk
  ! !
  ! ! DO jb = i_startblk_c, i_endblk_c
  ! !
  ! !   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
  ! !                      i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
  ! !     DO jk = slev, elev
  ! !       DO jc = i_startidx_c, i_endidx_c
  ! !
  ! !          IF ( v_base%lsm_c(jc,jk,jb) >= boundary ) THEN
  ! !            z_div_grad_u(jc,jk,jb)%x = 0.0_wp
  ! !          ELSE
  ! !           z_div_grad_u(jc,jk,jb)%x =  &
  ! !             z_grad_u(iidx(jc,jb,1),jk,iblk(jc,jb,1))%x * p_int_state(1)%geofac_div(jc,1,jb) + &
  ! !             z_grad_u(iidx(jc,jb,2),jk,iblk(jc,jb,2))%x * p_int_state(1)%geofac_div(jc,2,jb) + &
  ! !             z_grad_u(iidx(jc,jb,3),jk,iblk(jc,jb,3))%x * p_int_state(1)%geofac_div(jc,3,jb)
  ! !         ENDIF
  ! !       END DO
  ! !     END DO
  ! !   END DO
  ! !
  ! ! !Step 3: Map divergence back to edges
  ! ! CALL map_cell2edges( p_patch, z_div_grad_u, z_tmp)
  ! ! vort_flux_out = vort_flux_in + 1000000.0_wp*z_tmp
  ! !
  ! !
  ! ! END FUNCTION laplacian4vortex_flux
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation.
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
  !!  mpi parallelized LL
  !!
  SUBROUTINE veloc_adv_vert_mimetic_rot( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag    
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER :: z_dolic
    REAL(wp), POINTER :: del_zlev_i(:)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER         :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
   !-----------------------------------------------------------------------

    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(3) = 0.0_wp

    z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3) = 0.0_wp

    slev = 1
    elev = n_zlev

    !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
    !This requires appropriate boundary conditions
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)

        IF(z_dolic>=min_dolic)THEN
          del_zlev_i => p_patch_3D%p_patch_1D(1)%prism_center_dist_c(jc,:,jb)
          !1a) ocean surface
          z_adv_u_i(jc,slev,jb)%x =               &
            & p_diag%w(jc,slev,jb)*&  !/v_base%del_zlev_i(slev)
            & (p_diag%p_vn(jc,slev,jb)%x - p_diag%p_vn(jc,slev+1,jb)%x)/del_zlev_i(slev)

          ! 1b) ocean interior
          DO jk = slev+1, z_dolic-1
            z_adv_u_i(jc,jk,jb)%x                  &
             & = p_diag%w(jc,jk,jb)* & !/ v_base%del_zlev_i(jk)
             &(p_diag%p_vn(jc,jk-1,jb)%x - p_diag%p_vn(jc,jk,jb)%x)/del_zlev_i(jk)
          END DO
          ! z_adv_u_i(jc,slev,jb)%x=0.0_wp!z_adv_u_i(jc,slev+1,jb)%x
        ENDIF
      END DO
    END DO

    ! ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
    ! ! This mapping is the transposed of the vertical differencing.!
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!  v_base%dolic_c(jc,jb)
        IF(z_dolic>=min_dolic)THEN
          del_zlev_i => p_patch_3D%p_patch_1D(1)%prism_center_dist_c(jc,:,jb)
          DO jk = slev+1,z_dolic-1
            !This seems to work well
            !z_adv_u_m(jc,jk,jb)%x &
            !& = 0.5_wp*(z_adv_u_i(jc,jk,jb)%x+z_adv_u_i(jc,jk+1,jb)%x)
            !Map back from cell top to cell center and weight according to thickness
            z_adv_u_m(jc,jk,jb)%x &
            & = (del_zlev_i(jk)  *z_adv_u_i(jc,jk,jb)%x    &
            & +  del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)%x) &
            & / (del_zlev_i(jk+1)+del_zlev_i(jk))
          END DO
          ! 2c) ocean bottom
          !z_adv_u_m(jc,z_dolic,jb)%x =0.0_wp!=  z_adv_u_i(jc,z_dolic,jb)%x
        ENDIF
      END DO
    END DO

    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( p_patch_3D, z_adv_u_m, veloc_adv_vert_e,p_op_coeff)

    CALL sync_patch_array(SYNC_E, p_patch, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('VertMimRot: V.Adv. Final'    ,veloc_adv_vert_e         ,str_module,idt_src)
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_vert_mimetic_rot
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation. The rotational form excludes
  !! the flux form of the vertical advection, instead it implies the use of
  !! $w \partial_z Pv$. In this subroutine the vertical advection is discretized
  !! by transforming into flux-form minus a correction term, 
  !! All calculations are carried out at cell centers and are 
  !! mapped to edges at the end.
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
  !!  mpi parallelized LL
  !!
  SUBROUTINE veloc_adv_vert_mim_rot_flux2( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag    
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER :: z_dolic
    REAL(wp), POINTER :: del_zlev_i(:)
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp)                      :: z_w_ave  (nproma,n_zlev,  p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp)                      :: z_w_diff (nproma,n_zlev-1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,  p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    !-----------------------------------------------------------------------

    slev = 1
    elev = n_zlev

    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(3) = 0.0_wp

    z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3) = 0.0_wp

    z_w_ave (1:nproma,1:n_zlev,  1:p_patch%nblks_c) = 0.0_wp
    z_w_diff(1:nproma,1:n_zlev-1,1:p_patch%nblks_c) = 0.0_wp


    !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
    !This requires appropriate boundary conditions
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)

        IF(z_dolic>=min_dolic)THEN
          !del_zlev_m=>p_diag%inv_prism_thick_c(jc,:,jb)
          del_zlev_m => p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,jb)
          DO jk = slev, z_dolic-1
            z_w_ave(jc,jk,jb) = 0.5_wp*       (p_diag%w(jc,jk,jb)+p_diag%w(jc,jk+1,jb))
            z_w_diff(jc,jk,jb)=del_zlev_m(jk)*(p_diag%w(jc,jk,jb)-p_diag%w(jc,jk+1,jb))
             !&p_diag%inv_prism_thick_c(jc,jk,jb)!/v_base%del_zlev_m(jk)
          END DO
          z_w_ave(jc,z_dolic,jb)= p_diag%w(jc,z_dolic,jb)
        !ELSE
        !  z_w_ave(jc,:,jb)=0.0_wp
        ENDIF
      END DO
    END DO

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
        IF(z_dolic>=min_dolic)THEN
          !del_zlev_i=>p_diag%inv_prism_center_dist_c(jc,:,jb)
          del_zlev_i => p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,:,jb)
          DO jk = slev,z_dolic-1
            !The last term is the correction term to the constructed flux form.
            !The result of this calculation lives on prism top or bottom.
            z_adv_u_i(jc,jk,jb)%x =                            &
              & (z_w_ave(jc,jk,jb)  *p_diag%p_vn(jc,jk,jb)%x   &
              & -z_w_ave(jc,jk+1,jb)*p_diag%p_vn(jc,jk+1,jb)%x)&
              &*del_zlev_i(jk)                                 &!!/v_base%del_zlev_i(jk)&
              & -z_w_diff(jc,jk,jb) * p_diag%p_vn(jc,jk,jb)%x
        END DO
        ENDIF
      END DO
    END DO
! write(*,*)'vert 2: A max/min vert adv:',jk,&
! & maxval(z_adv_u_i(:,jk,:)%x(1)), minval(z_adv_u_m(:,jk,:)%x(1)),&
! & maxval(z_adv_u_i(:,jk,:)%x(2)), minval(z_adv_u_m(:,jk,:)%x(2)),&
! & maxval(z_adv_u_i(:,jk,:)%x(3)), minval(z_adv_u_m(:,jk,:)%x(3))

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
        IF(z_dolic>=min_dolic)THEN
          !del_zlev_i=>p_diag%prism_center_dist_c(jc,:,jb)
          del_zlev_i=>p_patch_3D%p_patch_1D(1)%prism_center_dist_c(jc,:,jb)

          DO jk = slev,z_dolic-1!DO jk = slev+1,z_dolic-1
            !This seems to work well
            !z_adv_u_m(jc,jk,jb)%x &
            !& = 0.5_wp*(z_adv_u_i(jc,jk,jb)%x+z_adv_u_i(jc,jk+1,jb)%x)
            !Map back from cell top to cell center and weight according to thickness
            !  z_adv_u_m(jc,jk,jb)%x &
            !    & = (v_base%del_zlev_i(jk)  *z_adv_u_i(jc,jk,jb)%x&
            !    & +  v_base%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)%x) &
            !    & / (v_base%del_zlev_i(jk+1)+v_base%del_zlev_i(jk))
            z_adv_u_m(jc,jk,jb)%x &
              & = (del_zlev_i(jk)  *z_adv_u_i(jc,jk,jb)%x&
              & +  del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)%x) &
              & / (del_zlev_i(jk+1)+del_zlev_i(jk))
          END DO
          ! 2c) ocean bottom
          !z_adv_u_m(jc,z_dolic,jb)%x =0.0_wp!=  z_adv_u_i(jc,z_dolic,jb)%x
        ENDIF
      END DO
    END DO

    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( p_patch_3D, z_adv_u_m, veloc_adv_vert_e,p_op_coeff)


    CALL sync_patch_array(SYNC_E, p_patch, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    !CALL dbg_print('VertMimRot2: z_adv_u_m%x(1)' ,z_adv_u_m%x(1)           ,str_module,idt_src)
    CALL dbg_print('VertMimRot2: VelAdv Final'   ,veloc_adv_vert_e         ,str_module,idt_src)
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_vert_mim_rot_flux2
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation. The rotational form excludes
  !! the flux form of the vertical advection, instead it implies the use of
  !! $w \partial_z Pv$. In this subroutine the vertical advection is discretized
  !! by transforming into flux-form minus a correction term, 
  !! All calculations are carried out at cell centers and are 
  !! mapped to edges at the end.
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
  !!  mpi parallelized LL
  !!
  SUBROUTINE veloc_adv_vert_mimetic_rot_flux( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag    
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER :: z_dolic
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp)                      :: z_w_diff (nproma,n_zlev-1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    !-----------------------------------------------------------------------
    slev = 1
    elev = n_zlev

    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(3) = 0.0_wp

    z_w_diff(1:nproma,1:n_zlev-1,1:p_patch%nblks_c) = 0.0_wp

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
        IF(z_dolic>=min_dolic)THEN
          !del_zlev_m=>p_diag%inv_prism_thick_c(jc,:,jb)
          del_zlev_m=>p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,jb)

          DO jk = slev, z_dolic-1
            z_w_diff(jc,jk,jb)=del_zlev_m(jk)*(p_diag%w(jc,jk,jb)-p_diag%w(jc,jk+1,jb))
             !/v_base%del_zlev_m(jk)
          END DO

          !jk=slev          
          !z_adv_u_i(jc,jk,jb)%x = del_zlev_m(jk)*0.5_wp*&   
          !& (p_diag%w(jc,jk,jb)  *(p_diag%p_vn(jc,jk,jb)%x + p_diag%p_vn(jc,jk,jb)%x) &
          !& -p_diag%w(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x + p_diag%p_vn(jc,jk+1,jb)%x)) &
          !& -z_w_diff(jc,jk,jb) * p_diag%p_vn(jc,jk,jb)%x

          DO jk = slev+1,z_dolic-1
            !The last term is the correction term to the constructed flux form.
            !The result of this calculation lives on prism top or bottom.
            z_adv_u_i(jc,jk,jb)%x = del_zlev_m(jk)*0.5_wp*&   
              & (p_diag%w(jc,jk,jb)  *(p_diag%p_vn(jc,jk-1,jb)%x + p_diag%p_vn(jc,jk,jb)%x)    &
              & -p_diag%w(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x   + p_diag%p_vn(jc,jk+1,jb)%x)) &
              & -z_w_diff(jc,jk,jb) * p_diag%p_vn(jc,jk,jb)%x
          END DO
        ENDIF
      END DO
    END DO
    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( p_patch_3D, z_adv_u_i, veloc_adv_vert_e, p_op_coeff)

    CALL sync_patch_array(SYNC_E, p_patch, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    !CALL dbg_print('VertMimRot2: z_adv_u_m%x(1)' ,z_adv_u_m%x(1)           ,str_module,idt_src)
    CALL dbg_print('VertMimRot: VelAdv Final'   ,veloc_adv_vert_e(:,2,:)         ,str_module,idt_src)
    !---------------------------------------------------------------------
  END SUBROUTINE veloc_adv_vert_mimetic_rot_flux
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation. The rotational form excludes
  !! the flux form of the vertical advection, instead it implies the use of
  !! $w \partial_z Pv$. In this subroutine the vertical advection is discretized
  !! by transforming into flux-form minus a correction term, 
  !! All calculations are carried out at cell centers and are 
  !! mapped to edges at the end.
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
  !!  mpi parallelized LL
  !!
  SUBROUTINE veloc_adv_vert_mimetic_div( p_patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag    
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jc, jk, jb
    INTEGER :: i_startidx, i_endidx
    INTEGER :: z_dolic
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp)                      :: z_w_diff (nproma,n_zlev-1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    !-----------------------------------------------------------------------
    slev = 1
    elev = n_zlev

    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(3) = 0.0_wp

    z_w_diff(1:nproma,1:n_zlev-1,1:p_patch%nblks_c) = 0.0_wp

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
        IF(z_dolic>=min_dolic)THEN
          !del_zlev_m=>p_diag%inv_prism_thick_c(jc,:,jb)
          del_zlev_m=>p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,jb)

          DO jk = slev, z_dolic-1
            z_w_diff(jc,jk,jb)=del_zlev_m(jk)*(p_diag%w(jc,jk,jb)-p_diag%w(jc,jk+1,jb))
             !/v_base%del_zlev_m(jk)
          END DO

          jk=slev          
          z_adv_u_i(jc,jk,jb)%x = del_zlev_m(jk)*0.5_wp*&   
          & (p_diag%w(jc,jk,jb)*(p_diag%p_vn(jc,jk,jb)%x + p_diag%p_vn(jc,jk,jb)%x) &
          & -p_diag%w(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x + p_diag%p_vn(jc,jk+1,jb)%x))

          DO jk = slev+1,z_dolic-1
            !The last term is the correction term to the constructed flux form.
            !The result of this calculation lives on prism top or bottom.
            z_adv_u_i(jc,jk,jb)%x = del_zlev_m(jk)*0.5_wp*&   
              & (p_diag%w(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x + p_diag%p_vn(jc,jk,jb)%x)    &
              & -p_diag%w(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x + p_diag%p_vn(jc,jk+1,jb)%x)) 
          END DO
        ENDIF
      END DO
    END DO
    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( p_patch_3D, z_adv_u_i, veloc_adv_vert_e, p_op_coeff)

    CALL sync_patch_array(SYNC_E, p_patch, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    !CALL dbg_print('VertMimRot2: z_adv_u_m%x(1)' ,z_adv_u_m%x(1)           ,str_module,idt_src)
    CALL dbg_print('VertMimDiv: VelAdv Final'   ,veloc_adv_vert_e         ,str_module,idt_src)
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_vert_mimetic_div
  !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! Computes vertical advection of a (edge based) horizontal vector field that
!   !! suits to divergence form of velocity equation.
!   !! The vertical derivative of the velocity vector at circumcenters that
!   !! is reconstructed from edge data is calculated and then multiplied by
!   !! the vertical velocity. The product is mapped from top of the computational
!   !! prism to the middle (still at centers) via the transposed of vertical differentiation
!   !! and then transformed to edges.
!   !!
!   !! IMPORTANT: It is assumed that the velocity vector reconstruction from
!   !! edges to cells has been done before.
!   !!
!   !! input:  lives on cells (velocity points)
!   !! output: lives on edges (velocity points)
!   !!
!   !! @par Revision History
!   !! Developed  by  Peter Korn, MPI-M (2010).
!   !!  mpi parallelized LL
!   !!
!   SUBROUTINE veloc_adv_vert_mimetic_div( p_patch, p_diag,p_op_coeff, veloc_adv_vert_e)
! 
!     TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!     TYPE(t_hydro_ocean_diag)          :: p_diag
!     TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
!     REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,p_patch%nblks_e)
! 
!     !local variables
!     INTEGER :: slev, elev     ! vertical start and end level
!     INTEGER :: jc, jk, jb
!     INTEGER :: i_startidx, i_endidx
!     INTEGER :: z_dolic
!     TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c)
!     TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
!     TYPE(t_subset_range), POINTER :: all_cells
!     !-----------------------------------------------------------------------
!     all_cells => p_patch%cells%all
!     !-----------------------------------------------------------------------
!     slev = 1
!     elev = n_zlev
! 
!     z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(1) = 0.0_wp
!     z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(2) = 0.0_wp
!     z_adv_u_i(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(3) = 0.0_wp
! 
!     z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1) = 0.0_wp
!     z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2) = 0.0_wp
!     z_adv_u_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3) = 0.0_wp
! 
! 
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
!       DO jc = i_startidx, i_endidx
!         z_dolic = v_base%dolic_c(jc,jb)
! 
!         IF(z_dolic>=min_dolic)THEN
!           ! 1a) ocean surface  !Code below explodes: Use upper boundary condition for d_z u ?
!           z_adv_u_i(jc,slev,jb)%x = p_diag%w(jc,slev,jb)*p_diag%p_vn(jc,slev,jb)%x 
! 
!           ! 1b) ocean interior
!           DO jk = slev+1, z_dolic-1
!             z_adv_u_i(jc,jk,jb)%x&
!               & = 0.5_wp*p_diag%w(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x + p_diag%p_vn(jc,jk,jb)%x)
!           END DO
!         ENDIF
!       END DO
!     END DO
! 
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
!       DO jc = i_startidx, i_endidx
!         z_dolic = v_base%dolic_c(jc,jb)
!         IF(z_dolic>=min_dolic)THEN
!           DO jk = slev,z_dolic-1
!             z_adv_u_m(jc,jk,jb)%x &
!               & = (z_adv_u_i(jc,jk,jb)%x+z_adv_u_i(jc,jk+1,jb)%x)/v_base%del_zlev_m(jk)
!           END DO
!         ENDIF
!       END DO
!       ! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
!       ! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
!     END DO
! 
!     ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
!     CALL map_cell2edges_3D( p_patch, z_adv_u_m, veloc_adv_vert_e,p_op_coeff)
! 
!     CALL sync_patch_array(SYNC_E, p_patch, veloc_adv_vert_e)
! 
!     !---------Debug Diagnostics-------------------------------------------
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('VertMimDiv: z_adv_u_m%x(1)'  ,z_adv_u_m%x(1)           ,str_module,idt_src)
!     CALL dbg_print('VertMimDiv: VelAdv Final'    ,veloc_adv_vert_e         ,str_module,idt_src)
!     !---------------------------------------------------------------------
! 
!   END SUBROUTINE veloc_adv_vert_mimetic_div
!   !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! Computes horizontal advection of a (edge based) vector field.
!   !!
!   !! Computes rotational term of a vector field given by its components in
!   !! the directions normal to triangle edges and the gradient of the kinetic energy
!   !! which is calculated using the reconstructed velocity at cell centers. Both
!   !! terms are combined and constitute the horizontal velocity advection.
!   !!
!   !!IMPORTANT: It is assumed that the reconstruction of the tangential velocity
!   !!           has been done before.
!   !1
!   !! input:  lives on edges (velocity points)
!   !! output: lives on edges (velocity points)
!   !!
!   !! @par Revision History
!   !! Developed  by  Peter Korn, MPI-M (2010).
!   !!  mpi parallelized LL
!   !!
!   SUBROUTINE veloc_adv_horz_rbf( p_patch, vn, p_diag, grad_coeff, veloc_adv_horz_e, p_int)
!     !
!     !
!     !  patch on which computation is performed
!     !
!     TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! 
!     !
!     ! normal and tangential velocity  of which advection is computed
!     !
!     REAL(wp), INTENT(inout) :: vn(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
!     !
!     !diagnostic ocean state stores horizontally advected velocity
!     !
!     TYPE(t_hydro_ocean_diag) :: p_diag
! 
!     REAL(wp), INTENT(in)    :: grad_coeff(:,:,:)
!     !
!     ! variable in which horizontally advected velocity is stored
!     !
!     REAL(wp), INTENT(out) :: veloc_adv_horz_e(:,:,:)
!     !
!     !Interpolation necessary just for testing
!     TYPE(t_int_state),TARGET,INTENT(in)  :: p_int
! 
! 
!     INTEGER :: slev, elev     ! vertical start and end level
!     INTEGER :: jk, jb, je, jc
!     INTEGER :: i_startidx_e, i_endidx_e
!     INTEGER :: i_startidx_c, i_endidx_c
!     INTEGER :: i_v1_idx, i_v1_blk, i_v2_idx, i_v2_blk
!     INTEGER :: jev, ile,ibe, i_v1_ctr, i_v2_ctr
!     REAL(wp) :: z_e  (nproma,n_zlev,p_patch%nblks_e)
!     REAL(wp) :: z_vort_glb(nproma,n_zlev,p_patch%nblks_v)
!     REAL(wp) :: z_grad_ekin_rbf(nproma,n_zlev,p_patch%nblks_e)
!     REAL(wp) :: z_kin_rbf_e(nproma,n_zlev,p_patch%nblks_e)
!     INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3
!     REAL(wp) :: z_vort_e(nproma,n_zlev,p_patch%nblks_e)
!     REAL(wp) :: z_vort_flx_rbf(nproma,n_zlev,p_patch%nblks_e)
!     REAL(wp) :: z_kin_e_rbf(nproma,n_zlev,p_patch%nblks_e)
!     REAL(wp) :: z_weight_e1,z_weight_e2, z_weight_e3!, z_weight
!     TYPE(t_subset_range), POINTER :: all_edges, owned_edges, all_cells
!     !-----------------------------------------------------------------------
!     all_edges   => p_patch%edges%all
!     owned_edges => p_patch%edges%owned
!     all_cells   => p_patch%cells%all
! 
!     ! #slo# set local variable to zero due to nag -nan compiler-option
!     z_e             (:,:,:) = 0.0_wp
!     z_vort_glb      (:,:,:) = 0.0_wp
!     !z_vort_flx      (:,:,:) = 0.0_wp
!     veloc_adv_horz_e(:,:,:) = 0.0_wp
!     z_vort_e        (:,:,:) = 0.0_wp
!     z_vort_flx_rbf  (:,:,:) = 0.0_wp
!     z_kin_e_rbf     (:,:,:) = 0.0_wp
!     z_grad_ekin_rbf (:,:,:) = 0.0_wp
! 
!     slev = 1
!     elev = n_zlev
! 
!     CALL rbf_vec_interpol_edge( vn,       &
!       & p_patch,  &
!       & p_int,    &
!       & p_diag%vt,&
!       & opt_slev=slev,opt_elev=elev)
! 
!     CALL sync_patch_array(SYNC_E, p_patch, p_diag%v)
! 
! 
!     DO jb = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!       DO jk = slev, elev
!         DO je=i_startidx_e, i_endidx_e
!           IF ( v_base%lsm_e(je,jk,jb) == boundary ) THEN
!             p_diag%vt(je,jk,jb) = 0.0_wp
!             vn(je,jk,jb) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
!     END DO
! 
!     CALL rot_vertex_ocean_rbf(p_patch,vn, p_diag%vt, p_diag%vort)
!     ! CALL verts2edges_scalar( p_diag%vort, p_patch, p_int%v_1o2_e, &
!     !                          z_vort_e, opt_slev=slev,opt_elev=elev, opt_rlstart=3)
!     CALL sync_patch_array(SYNC_V, p_patch, p_diag%vort)
! 
! 
!     DO jb = owned_edges%start_block, owned_edges%end_block
!       CALL get_index_range(owned_edges, jb, i_startidx_e, i_endidx_e)
!       DO jk = slev, slev
!         DO je=i_startidx_e, i_endidx_e
!           i_v1_idx = p_patch%edges%vertex_idx(je,jb,1)
!           i_v1_blk = p_patch%edges%vertex_blk(je,jb,1)
!           i_v2_idx = p_patch%edges%vertex_idx(je,jb,2)
!           i_v2_blk = p_patch%edges%vertex_blk(je,jb,2)
!           !count wet edges in vertex 1
!           i_v1_ctr = 0
!           DO jev = 1, p_patch%verts%num_edges(i_v1_idx,i_v1_blk)
!             ile = p_patch%verts%edge_idx(i_v1_idx,i_v1_blk,jev)
!             ibe = p_patch%verts%edge_blk(i_v1_idx,i_v1_blk,jev)
!             IF ( v_base%lsm_e(ile,jk,ibe) == sea ) THEN
!               i_v1_ctr = i_v1_ctr +1
!             ENDIF
!           END DO
!           !count wet edges in vertex 2
!           i_v2_ctr = 0
!           DO jev = 1, p_patch%verts%num_edges(i_v2_idx,i_v2_blk)
!             ile = p_patch%verts%edge_idx(i_v2_idx,i_v2_blk,jev)
!             ibe = p_patch%verts%edge_blk(i_v2_idx,i_v2_blk,jev)
!             IF ( v_base%lsm_e(ile,jk,ibe) == sea ) THEN
!               i_v2_ctr = i_v2_ctr +1
!             ENDIF
!           END DO
!           IF(   i_v1_ctr==p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr==p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN
! 
!             z_vort_e(je,jk,jb) =&
!               & 0.5_wp*(p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & +        p_diag%vort(i_v2_idx,jk,i_v2_blk))
! 
!           ELSEIF(   i_v1_ctr==p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr <p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN
! 
!             z_vort_e(je,jk,jb) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & + REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
! 
!           ELSEIF(   i_v1_ctr<p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr==p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN
! 
!             z_vort_e(je,jk,jb) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & + REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
! 
!           ELSEIF(   i_v1_ctr<p_patch%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr<p_patch%verts%num_edges(i_v2_idx,i_v2_blk))THEN
! 
!             z_vort_e(je,jk,jb) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & + REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
!           ELSE
!             z_vort_e(je,jk,jb) = 0.0_wp
!           ENDIF
!           !  IF(jk==1)THEN
!           !  IF(v_base%lsm_e(je,jk,jb)/=2)THEN
!           ! ! IF(v_base%lsm_e(je,jk,jb)==0)THEN
!           ! IF(i_v1_ctr==6.and.i_v2_ctr==6.and.v_base%lsm_e(je,jk,jb)==sea)THEN
!           ! ELSE
!           !  write(101,*)'vert ctr', jk,je,jb, i_v1_ctr, i_v2_ctr, v_base%lsm_e(je,jk,jb)!,&
!           ! ENDIF
!           ! ! ENDIF
!           !  ENDIF
!           !  ENDIF
!         END DO
!       END DO
!     ENDDO
!     CALL sync_patch_array(SYNC_E, p_patch, z_vort_e)
! 
! 
!     DO jb = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!       DO jk = slev, elev
!         DO je=i_startidx_e, i_endidx_e
!           IF ( v_base%lsm_e(je,jk,jb) == sea ) THEN
!             z_vort_flx_rbf(je,jk,jb) =&
!               & p_diag%vt(je,jk,jb)*(p_patch%edges%f_e(je,jb)+z_vort_e(je,jk,jb))
!             !          & p_diag%vt(je,jk,jb)*p_patch%edges%f_e(je,jb)
!           ELSE
!             z_vort_flx_rbf(je,jk,jb) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
!     ENDDO
! 
!     CALL rbf_vec_interpol_cell( vn, p_patch, p_int, p_diag%u,  &
!       & p_diag%v, opt_slev=slev, opt_elev=elev)
!     CALL sync_patch_array(SYNC_C, p_patch, p_diag%v)
! 
!     !write(*,*)'max/min vort flux:', MAXVAL(z_vort_flx_RBF(:,1,:)),MINVAL(z_vort_flx_RBF(:,1,:))
!     DO jb = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!       DO jk = slev, elev
!         DO je = i_startidx_e, i_endidx_e
!           ! calculate kinetic energy at edges from normal and tangential comp.
!           z_kin_rbf_e(je,jk,jb) =0.5_wp*(p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb)&
!             & +    vn(je,jk,jb)*vn(je,jk,jb) )
!         ENDDO
!       ENDDO
!     ENDDO
! 
!     !!$OMP END DO
!     !!$OMP END PARALLEL
!     ! Bilinear interpolation of kinetic energy from the edges to the cells
!     !    CALL edges2cells_scalar( z_kin_RBF_e,     &
!     !                           & p_patch,         &
!     !                           & p_int%e_bln_c_s, &
!     !                           & p_diag%kin,      &
!     !                           & opt_slev=slev,opt_elev=elev)
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!       DO jk = slev, elev
!         DO jc = i_startidx_c, i_endidx_c
!           IF ( v_base%lsm_c(jc,jk,jb) > sea_boundary ) THEN
!             p_diag%kin(jc,jk,jb) = 0.0_wp
!           ELSE
! 
!             ile1 = p_patch%cells%edge_idx(jc,jb,1)
!             ibe1 = p_patch%cells%edge_blk(jc,jb,1)
!             ile2 = p_patch%cells%edge_idx(jc,jb,2)
!             ibe2 = p_patch%cells%edge_blk(jc,jb,2)
!             ile3 = p_patch%cells%edge_idx(jc,jb,3)
!             ibe3 = p_patch%cells%edge_blk(jc,jb,3)
!             z_weight_e1 = 0.0_wp
!             z_weight_e2 = 0.0_wp
!             z_weight_e3 = 0.0_wp
!             IF(v_base%lsm_e(ile1,jk,ibe1)<= boundary)THEN
!               z_weight_e1 = p_patch%edges%area_edge(ile1,ibe1)
!             ENDIF
!             IF(v_base%lsm_e(ile2,jk,ibe2)<= boundary)THEN
!               z_weight_e2 = p_patch%edges%area_edge(ile2,ibe2)
!             ENDIF
!             IF(v_base%lsm_e(ile3,jk,ibe3)<= boundary)THEN
!               z_weight_e3 = p_patch%edges%area_edge(ile3,ibe3)
!             ENDIF
! 
!             !write(*,*)'weights',jc,jk,jb,z_weight_e1,z_weight_e2,z_weight_e3
!             p_diag%kin(jc,jk,jb) = (z_kin_rbf_e(ile1,jk,ibe1)*z_weight_e1&
!               & + z_kin_rbf_e(ile2,jk,ibe2)*z_weight_e2&
!               & + z_kin_rbf_e(ile3,jk,ibe3)*z_weight_e3)&
!               & /(z_weight_e1+z_weight_e2+z_weight_e3)
!             !       p_diag%kin(jc,jk,jb) = 0.5_wp*(p_diag%u(jc,jk,jb)*p_diag%u(jc,jk,jb)&
!             !                                     &+p_diag%v(jc,jk,jb)*p_diag%v(jc,jk,jb))
!           ENDIF
!         END DO
!       END DO
!     END DO
! 
!    CALL grad_fd_norm_oce_3d( p_diag%kin, &
!       & p_patch,    &
!       & grad_coeff, &
!       & z_grad_ekin_rbf)
! 
! !     CALL grad_fd_norm_oce( p_diag%kin, &
! !       & p_patch,    &
! !       & z_grad_ekin_rbf, opt_slev=slev,opt_elev=elev)
! !     CALL sync_patch_array(SYNC_C, p_patch, z_grad_ekin_rbf)
! 
! 
!     !Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
!     DO jb = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!       DO jk = slev, elev
!         DO je = i_startidx_e, i_endidx_e
!           IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
!             veloc_adv_horz_e(je,jk,jb) =&
!               & z_vort_flx_rbf(je,jk,jb) + z_grad_ekin_rbf(je,jk,jb)
!           ELSE
!             veloc_adv_horz_e(je,jk,jb) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
!     END DO
! 
!     !---------Debug Diagnostics-------------------------------------------
!     idt_src=2  ! output print level (1-5, fix)
!     CALL dbg_print('HorzRBF: vorticity adv.'     ,veloc_adv_horz_e         ,str_module,idt_src)
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('HorzRBF: kin. energy'        ,p_diag%kin               ,str_module,idt_src)
!     CALL dbg_print('HorzRBF: vorticity'          ,p_diag%vort              ,str_module,idt_src)
!     idt_src=4  ! output print level (1-5, fix)
!     CALL dbg_print('HorzRBF: grad. kin. en.'     ,z_grad_ekin_rbf          ,str_module,idt_src)
!     CALL dbg_print('HorzRBF: vorticity_e'        ,z_vort_e                 ,str_module,idt_src)
!     CALL dbg_print('HorzRBF: vorticity flux'     ,z_vort_flx_rbf           ,str_module,idt_src)
!     !---------------------------------------------------------------------
! 
!   END SUBROUTINE veloc_adv_horz_rbf
!   !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! Computes vertical advection of a (edge based) horizontal vector field.
!   !! The vertical derivative of the velocity vector at circumcenters that
!   !! is reconstructed from edge data is calculated and then multiplied by
!   !! the vertical velocity. The product is mapped from top of the computational
!   !! prism to the middle (still at centers) via the transposed of vertical differentiation
!   !! and then transformed to edges.
!   !!
!   !! IMPORTANT: It is assumed that the velocity vector reconstruction from
!   !! edges to cells has been done before.
!   !!
!   !! input:  lives on cells (velocity points)
!   !! output: lives on edges (velocity points)
!   !!
!   !! @par Revision History
!   !! Developed  by  Peter Korn, MPI-M (2010).
!   !!  mpi parallelized LL
!    !!
!   SUBROUTINE veloc_adv_vert_rbf( p_patch, u_c, v_c, w_c, &
!     & top_bc_u_c, top_bc_v_c, &
!     & bot_bc_u_c,  bot_bc_v_c,&
!     & top_bc_w_c,  bot_bc_w_c,&
!     & veloc_adv_vert_e)
!     !
!     !  patch on which computation is performed
!     !
!     TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! 
!     !
!     ! Components of cell based variable which is vertically advected
!     REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
!     REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
!     REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
!     !
!     ! Top boundary condition for cell based variables
!     REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!     REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!     !
!     ! Bottom boundary condition for cell based variables
!     REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!     REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!     !
!     REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
!     REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
! 
!     ! variable in which horizontally advected velocity is stored
!     REAL(wp), INTENT(inout) :: veloc_adv_vert_e(:,:,:)
! 
!     INTEGER :: slev, elev     ! vertical start and end level
!     INTEGER :: jc, jk, jb, i_dolic
!     INTEGER :: i_startidx, i_endidx
! 
!     TYPE(t_subset_range), POINTER :: all_cells
! 
!     REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
!       & z_adv_v_i(nproma,n_zlev+1,p_patch%nblks_c),  &
!       & z_adv_u_m(nproma,n_zlev,p_patch%nblks_c),  &
!       & z_adv_v_m(nproma,n_zlev,p_patch%nblks_c)
!     !-----------------------------------------------------------------------
!     all_cells => p_patch%cells%all
! 
!     ! #slo# set local variable to zero due to nag -nan compiler-option
!     z_adv_u_i(:,:,:) = 0.0_wp
!     z_adv_v_i(:,:,:) = 0.0_wp
!     z_adv_u_m(:,:,:) = 0.0_wp
!     z_adv_v_m(:,:,:) = 0.0_wp
! 
!     slev = 1
!     elev = n_zlev
! 
!     !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
!     !This requires appropriate boundary conditions
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
!       DO jk = slev, elev
!         DO jc = i_startidx, i_endidx
!           !check if we have at least two layers of water
!           !  #slo# - 2011-04-01 - Is this really intended here
!           !  maybe this condition should be fulfilled everywhere
!           !  then it must be calculated in fill_vertical_domain
!           !  this condition could then be omitted here
!           IF (v_base%dolic_c(jc,jb) >= 2) THEN
!             !1b) ocean bottom
!             IF ( jk == v_base%dolic_c(jc,jb) ) THEN
!               ! u,v-component
!               z_adv_u_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_u_c(jc,jb)
!               z_adv_v_i(jc,jk+1,jb) = bot_bc_w_c(jc,jb)*bot_bc_v_c(jc,jb)
! 
!               !1c) ocean interior
!             ELSEIF( jk>slev .AND.  jk < v_base%dolic_c(jc,jb))THEN
!               ! u,v-component
!               z_adv_u_i(jc,jk,jb)&
!                 & = w_c(jc,jk,jb) *( u_c(jc,jk-1,jb) - u_c(jc,jk,jb) )&
!                 & / v_base%del_zlev_i(jk)
! 
!               z_adv_v_i(jc,jk,jb)&
!                 & = w_c(jc,jk,jb) *( v_c(jc,jk-1,jb) - v_c(jc,jk,jb) )&
!                 & / v_base%del_zlev_i(jk) !&
!               ! write(*,*)'vert adv:v: ',jk, jc,jb,w_c(jc,jk,jb) *( v_c(jc,jk,jb) - v_c(jc,jk-1,jb) ),&
!               ! &  v_c(jc,jk,jb), v_c(jc,jk-1,jb),v_base%del_zlev_i(jk-1)
!               ! write(*,*)'vert adv:u: ',jk, jc,jb,w_c(jc,jk,jb) *( u_c(jc,jk,jb) - u_c(jc,jk-1,jb) ),&
!               ! &  u_c(jc,jk,jb), u_c(jc,jk-1,jb)
! 
!             ENDIF  ! jk-condition
!           ENDIF    ! at least 2 vertical layers
!         END DO
!       END DO
!       !  write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
!       !  & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
!     END DO
! 
!     ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
!     ! This mapping is the transposed of the vertical differencing.
! 
!     !1) From surface down to one layer before bottom
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
!       DO jk = slev, elev-1
!         DO jc = i_startidx, i_endidx
!           !check if we are on land: To be replaced by 3D lsm
!           ! #slo# 2011-05-11 - replace by consistent formulation: vertical loop down to dolic
!           IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! 
!             z_adv_u_m(jc,jk,jb) &
!               & = (v_base%del_zlev_i(jk+1)*z_adv_u_i(jc,jk+1,jb)&
!               & +  v_base%del_zlev_i(jk)*z_adv_u_i(jc,jk,jb)) &
!               & / (v_base%del_zlev_i(jk+1)+v_base%del_zlev_i(jk))
! 
!             z_adv_v_m(jc,jk,jb)&
!               & = (v_base%del_zlev_i(jk+1)*z_adv_v_i(jc,jk+1,jb)&
!               & +  v_base%del_zlev_i(jk)*z_adv_v_i(jc,jk,jb))&
!               & / (v_base%del_zlev_i(jk+1)+v_base%del_zlev_i(jk))
!           ENDIF
!         END DO
!       END DO
!       ! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
!       ! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
!     END DO
!     !Bottom layer
!     !The value of v_base%del_zlev_i at the botom is 0.5*v_base%del_zlev_m
!     !The dimensioning of the firs arrays requires to seperate the vertical loop.
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
!       DO jc = i_startidx, i_endidx
!         IF ( v_base%dolic_c(jc,jb)>0 ) THEN  ! wet points only
!           i_dolic = v_base%dolic_c(jc,jb)
!           z_adv_u_m(jc,i_dolic,jb) =0.0_wp!&
!           !& = (0.5_wp*v_base%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,jb)&
!           !& +        v_base%del_zlev_i(elev)*z_adv_u_i(jc,elev,jb)) &
!           !& / (2.0_wp*v_base%del_zlev_m(elev))
! 
!           z_adv_v_m(jc,i_dolic,jb)=0.0_wp!&
!           !& = (0.5_wp*v_base%del_zlev_m(elev)*z_adv_v_i(jc,elev+1,jb)&
!           !&   +        v_base%del_zlev_i(elev)*z_adv_v_i(jc,elev,jb))&
!           !& / (2.0_wp*v_base%del_zlev_m(elev))
! 
!         END IF
!       END DO
!     END DO
! 
!     ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
! !     CALL primal_map_c2e( p_patch,&
! !       & z_adv_u_m, z_adv_v_m,&
! !       & veloc_adv_vert_e )
!     ! result is synced in the called funtion
! 
!   ! !---------Debug Diagnostics-------------------------------------------
!   !  NOT YET
!   ! idt_src=3  ! output print level (1-5, fix)
!   ! CALL dbg_print('VertRBF: kin. energy'        ,p_diag%kin               ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity'          ,p_diag%vort              ,str_module,idt_src)
!   ! idt_src=4  ! output print level (1-5, fix)
!   ! CALL dbg_print('VertRBF: grad. kin. en.'     ,z_grad_ekin_rbf          ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity_e'        ,z_vort_e                 ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity flux'     ,z_vort_flx_rbf           ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity adv.'     ,veloc_adv_horz_e         ,str_module,idt_src)
!   ! !---------------------------------------------------------------------
! 
!   END SUBROUTINE veloc_adv_vert_rbf
!   !-------------------------------------------------------------------------


END MODULE mo_oce_veloc_advection
