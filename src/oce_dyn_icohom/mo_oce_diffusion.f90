!>
!! Contains the implementation of velocity and tracer diffusion for the ICON ocean model.
!! 
!! 
!! @par Revision History
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
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
MODULE mo_oce_diffusion
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
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec !, gc2cc
USE mo_impl_constants,      ONLY: boundary, sea_boundary, MIN_DOLIC ! ,max_char_length
USE mo_parallel_config,     ONLY: nproma
USE mo_ocean_nml,           ONLY: n_zlev, iswm_oce, veloc_diffusion_order, veloc_diffusion_form
USE mo_run_config,          ONLY: dtime
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_oce_state,           ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
  &                               t_hydro_ocean_aux!, v_base
USE mo_model_domain,        ONLY: t_patch, t_patch_3D_oce
!USE mo_exception,           ONLY: message, finish!, message_text
USE mo_oce_physics,         ONLY: t_ho_params
USE mo_scalar_product,      ONLY: map_cell2edges_3D
USE mo_oce_math_operators,  ONLY: div_oce_3D, rot_vertex_ocean_3d,&
 &                                map_edges2vert_3D
USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array

IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module    = 'oceDiffusion'  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

!
! PUBLIC INTERFACE
!
INTEGER, PARAMETER  :: top=1
PUBLIC :: velocity_diffusion
PUBLIC :: veloc_diff_harmonic_div_grad, veloc_diff_biharmonic_div_grad
PUBLIC :: veloc_diff_harmonic_curl_curl, veloc_diff_biharmonic_curl_curl
PUBLIC :: velocity_diffusion_vert_mimetic
!PUBLIC :: velocity_diffusion_horz_rbf
!PUBLIC :: velocity_diffusion_vert_rbf
PUBLIC :: veloc_diffusion_vert_impl_hom
PUBLIC :: tracer_diffusion_horz
PUBLIC :: tracer_diffusion_vert_expl
PUBLIC :: tracer_diffusion_vert_impl_hom

CONTAINS

!-------------------------------------------------------------------------  
!
!  
!>
!! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
!!    implemented as P^T div( K_H grad P v).
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
!! mpi parallelized, sync required
SUBROUTINE velocity_diffusion(p_patch,p_patch_3D, vn_in, p_param, p_diag,p_op_coeff, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in)   :: p_patch
 TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  REAL(wp), INTENT(in)                :: vn_in(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_ho_params), INTENT(in)       :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag),INTENT(in) :: p_diag
  TYPE(t_operator_coeff),INTENT(in)   :: p_op_coeff
  REAL(wp), INTENT(INOUT)             :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)

  !Local variables
  !REAL(wp) :: z_lapl(nproma,n_zlev,p_patch%nblks_e)
  !INTEGER  :: jk
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        

  IF(veloc_diffusion_order==1)THEN

    !divgrad laplacian is chosen 
    IF(veloc_diffusion_form==2)THEN

      CALL veloc_diff_harmonic_div_grad( p_patch, p_patch_3D,  &
                                       & p_param,   &
                                       & p_diag,    &
                                       & p_op_coeff,&
                                       & laplacian_vn_out)

    ELSEIF(veloc_diffusion_form==1)THEN

       CALL veloc_diff_harmonic_curl_curl( vn_in, p_diag%vort, &
                                         & p_patch, p_patch_3D,&
                                         & p_param%K_veloc_h,  &
                                         & p_op_coeff,         &
                                         & p_diag%p_vn_dual,   &
                                         & laplacian_vn_out)
   ENDIF

  ELSEIF(veloc_diffusion_order==2)THEN

   IF(veloc_diffusion_form==2)THEN

      CALL veloc_diff_biharmonic_div_grad( p_patch,p_patch_3D,      &
                                         & p_param,      &
                                         & p_diag,       &
                                         & p_op_coeff,   &
                                         & laplacian_vn_out)
   ELSEIF(veloc_diffusion_form==1)THEN

       CALL veloc_diff_biharmonic_curl_curl( p_patch, p_patch_3D, &
                                           & vn_in,             &
                                           & p_diag%vort,       &
                                           & p_param%K_veloc_h, &
                                           & p_op_coeff,        &
                                           & p_diag%p_vn_dual,  &
                                           & laplacian_vn_out)
   ENDIF
 ENDIF

!    DO jk=1, n_zlev
!     write(*,*)'LAPLACIAN',jk,&
!     &maxval(laplacian_vn_out(:,jk,:)),minval(laplacian_vn_out(:,jk,:))!,&
!     &maxval(z_lapl(:,jk,:)),minval(z_lapl(:,jk,:))
!    END DO

 CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)


END SUBROUTINE velocity_diffusion
!-------------------------------------------------------------------------  
!
!  
!>
!! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
!!    implemented as P^T div( K_H grad P v).
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
!! mpi parallelized, sync required
SUBROUTINE veloc_diff_harmonic_div_grad(p_patch,p_patch_3D, p_param, p_diag,&
                                       &p_op_coeff, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  !REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
  REAL(wp), INTENT(INOUT)           :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)

  !Local variables
  INTEGER :: slev, elev
  INTEGER :: jk, jb, je,jc
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: i_startidx_c, i_endidx_c
  INTEGER :: i_startidx_e, i_endidx_e
  INTEGER :: idx_cartesian
  INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,p_patch%nblks_c)

  TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
!-------------------------------------------------------------------------------
  all_cells       => p_patch%cells%all
  all_edges       => p_patch%edges%all
  edges_in_domain => p_patch%edges%in_domain
!-------------------------------------------------------------------------------
  slev = 1
  elev = n_zlev

  laplacian_vn_out(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp

  ! loop over cells in local domain + halo
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = i_startidx_c, i_endidx_c
          z_div_grad_u(jc,jk,jb)%x =  0.0_wp
      END DO
    END DO
  END DO

  ! loop over edges in local domain + halo
  DO jb = all_edges%start_block, all_edges%end_block
    CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e
        z_grad_u(je,jk,jb)%x = 0.0_wp
      ENDDO
    END DO
  END DO

  !-------------------------------------------------------------------------------------------------------
  !Step 1: Calculate gradient of cell velocity vector.
  !Result is a gradient vector, located at edges
  !Step 2: Multiply each component of gradient vector with mixing coefficients
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)

    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e
        !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
        IF (p_patch_3D%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)

          z_grad_u(je,jk,jb)%x = p_param%K_veloc_h(je,jk,jb)*  &
            &                  ( p_diag%p_vn(il_c2,jk,ib_c2)%x &
            &                  - p_diag%p_vn(il_c1,jk,ib_c1)%x)&
            &                  / p_patch%edges%dual_edge_length(je,jb)
      ENDIF 
      ENDDO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_E, p_patch,z_grad_u(:,:,:)%x(idx_cartesian) )
  END DO

  !Step 2: Apply divergence to each component of mixing times gradient vector
  iidx => p_patch%cells%edge_idx
  iblk => p_patch%cells%edge_blk

  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = i_startidx_c, i_endidx_c

         !IF ( v_base%lsm_oce_c(jc,jk,jb) >= boundary ) THEN
         IF (p_patch_3D%lsm_oce_c(jc,jk,jb) >= boundary) THEN
           z_div_grad_u(jc,jk,jb)%x = 0.0_wp
         ELSE
          z_div_grad_u(jc,jk,jb)%x =  &
            z_grad_u(iidx(jc,jb,1),jk,iblk(jc,jb,1))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,1)+&
            z_grad_u(iidx(jc,jb,2),jk,iblk(jc,jb,2))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,2)+&
            z_grad_u(iidx(jc,jb,3),jk,iblk(jc,jb,3))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,3)
        ENDIF
      END DO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_C, p_patch,z_div_grad_u(:,:,:)%x(idx_cartesian) )
  END DO

  !Step 3: Map divergence back to edges
  CALL map_cell2edges_3D( p_patch,p_patch_3D, z_div_grad_u, laplacian_vn_out, p_op_coeff)

  ! write(*,*)'lapla',maxval(p_param%K_veloc_h(:,1,:)*p_patch%edges%primal_edge_length),&
  ! &minval(p_param%K_veloc_h(:,1,:)*p_patch%edges%primal_edge_length)
  ! !  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)
END SUBROUTINE veloc_diff_harmonic_div_grad
!-------------------------------------------------------------------------  
!
!  
!>
!! !  SUBROUTINE calculates horizontal diffusion of edge velocity via bilaplacian diffusion
!!    implemented as P^T (div K grad(divgrad P v)). Due to the position of the mixing matrix
!!    which is following MPI-OM, the operator can not be written as simple iteration of
!!    the laplacian in divgrad form.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
!! mpi parallelized, sync required
SUBROUTINE veloc_diff_biharmonic_div_grad(p_patch,p_patch_3D, p_param, p_diag,&
                                         &p_op_coeff, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
 TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  !REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
  REAL(wp), INTENT(INOUT)           :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)

  !Local variables
  REAL(wp) :: z_laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)
  INTEGER :: slev, elev
  INTEGER :: jk, jb, je,jc
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: i_startidx_c, i_endidx_c
  INTEGER :: i_startidx_e, i_endidx_e
  INTEGER :: idx_cartesian
  INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,p_patch%nblks_c)

  TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
!-------------------------------------------------------------------------------
  all_cells       => p_patch%cells%all
  all_edges       => p_patch%edges%all
  edges_in_domain => p_patch%edges%in_domain
!-------------------------------------------------------------------------------

  slev = 1
  elev = n_zlev


  laplacian_vn_out  (1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
  z_laplacian_vn_out(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp


  ! loop over cells in local domain + halo
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = i_startidx_c, i_endidx_c
          z_div_grad_u(jc,jk,jb)%x =  0.0_wp
      END DO
    END DO
  END DO

  ! loop over edges in local domain + halo
  DO jb = all_edges%start_block, all_edges%end_block
    CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e
        z_grad_u(je,jk,jb)%x = 0.0_wp
      ENDDO
    END DO
  END DO

  !-------------------------------------------------------------------------------------------------------
  !Step 1: Calculate gradient of cell velocity vector.
  !Result is a gradient vector, located at edges
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)

    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e

      !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
      IF (p_patch_3D%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        z_grad_u(je,jk,jb)%x = (p_diag%p_vn(il_c2,jk,ib_c2)%x &
          &                   - p_diag%p_vn(il_c1,jk,ib_c1)%x)&
          &                  / p_patch%edges%dual_edge_length(je,jb)
      ENDIF 
      ENDDO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_E, p_patch,z_grad_u(:,:,:)%x(idx_cartesian) )
  END DO

  !Step 2: Apply divergence to each component of gradient vector
  iidx => p_patch%cells%edge_idx
  iblk => p_patch%cells%edge_blk

  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = i_startidx_c, i_endidx_c

         !IF ( v_base%lsm_oce_c(jc,jk,jb) >= boundary ) THEN
         IF (p_patch_3D%lsm_oce_c(jc,jk,jb) >= boundary) THEN
           z_div_grad_u(jc,jk,jb)%x = 0.0_wp
         ELSE
          z_div_grad_u(jc,jk,jb)%x =  &
            z_grad_u(iidx(jc,jb,1),jk,iblk(jc,jb,1))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,1)+&
            z_grad_u(iidx(jc,jb,2),jk,iblk(jc,jb,2))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,2)+&
            z_grad_u(iidx(jc,jb,3),jk,iblk(jc,jb,3))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,3)
        ENDIF
      END DO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_C, p_patch,z_div_grad_u(:,:,:)%x(idx_cartesian) )
  END DO


  !Step 4: Repeat the application of div and grad and take the mixing coefficients into account
  !First the grad of previous result
  !now imes the mixiing/friction coefficient
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e

      !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
      IF (p_patch_3D%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        z_grad_u(je,jk,jb)%x = p_param%K_veloc_h(je,jk,jb)*&
          &                   (z_div_grad_u(il_c2,jk,ib_c2)%x &
          &                  - z_div_grad_u(il_c1,jk,ib_c1)%x)&
          &                  / p_patch%edges%dual_edge_length(je,jb)
      ENDIF 
      ENDDO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_E, p_patch,z_grad_u(:,:,:)%x(idx_cartesian) )
  END DO

  !Step 5: Apply divergence to each component of gradient vector
  iidx => p_patch%cells%edge_idx
  iblk => p_patch%cells%edge_blk

  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = i_startidx_c, i_endidx_c

         !IF ( v_base%lsm_oce_c(jc,jk,jb) >= boundary ) THEN
         IF (p_patch_3D%lsm_oce_c(jc,jk,jb) >= boundary) THEN
           z_div_grad_u(jc,jk,jb)%x = 0.0_wp
         ELSE
          z_div_grad_u(jc,jk,jb)%x =  &
            z_grad_u(iidx(jc,jb,1),jk,iblk(jc,jb,1))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,1)+&
            z_grad_u(iidx(jc,jb,2),jk,iblk(jc,jb,2))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,2)+&
            z_grad_u(iidx(jc,jb,3),jk,iblk(jc,jb,3))%x&
            & * p_op_coeff%div_coeff(jc,jk,jb,3)
        ENDIF
      END DO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_C, p_patch,z_div_grad_u(:,:,:)%x(idx_cartesian) )
  END DO


  !Step 6: Map divergence back to edges
  CALL map_cell2edges_3D( p_patch,p_patch_3D, z_div_grad_u,laplacian_vn_out,p_op_coeff)!,subset_range=all_cells)
  ! write(*,*)'lapla',maxval(p_param%K_veloc_h(:,1,:)*p_patch%edges%primal_edge_length),&
  ! &minval(p_param%K_veloc_h(:,1,:)*p_patch%edges%primal_edge_length)
  ! !  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)
END SUBROUTINE veloc_diff_biharmonic_div_grad
  !-------------------------------------------------------------------------
  !>
  !!  Computes  laplacian of a vector field in curl curl form.
  !!
  !! input:  lives on edges (velocity points)
  !! output: lives on edges
  !!
  !! @par Revision History
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE veloc_diff_harmonic_curl_curl( u_vec_e, vort, p_patch,p_patch_3D, k_h, p_op_coeff,&
    & p_vn_dual, nabla2_vec_e)
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
    REAL(wp), INTENT(in)                      :: u_vec_e(nproma,n_zlev,p_patch%nblks_e) 
    REAL(wp), INTENT(in)                      :: vort(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp), INTENT(in)                      :: k_h(:,:,:)
    TYPE(t_operator_coeff), INTENT(in)        :: p_op_coeff
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp), INTENT(inout)                   ::  nabla2_vec_e(nproma,n_zlev,p_patch%nblks_e)
    !
    !Local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    REAL(wp) ::  z_div_c(nproma,n_zlev,p_patch%nblks_c)!, &
    !REAL(wp) ::  z_vn_e(nproma,n_zlev,p_patch%nblks_e)
    !REAL(wp) ::  z_rot_v(nproma,n_zlev,p_patch%nblks_v)
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    edges_in_domain => p_patch%edges%in_domain
    !-----------------------------------------------------------------------

    slev = 1
    elev = n_zlev

    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    ! compute divergence of vector field
    CALL div_oce_3d( u_vec_e, p_patch, p_op_coeff%div_coeff, z_div_c)
    CALL sync_patch_array(SYNC_C,p_patch,z_div_c)

    ! compute rotation of vector field for the ocean
    !CALL rot_vertex_ocean_3D( p_patch, u_vec_e, p_vn_dual, p_op_coeff, z_rot_v)!
    !CALL sync_patch_array(SYNC_V,p_patch,z_rot_v)
    !z_rot_v=vort
    !
    !  loop through all patch edges (and blocks)
    !
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO je = i_startidx, i_endidx
        DO jk = slev, elev

          !IF(v_base%lsm_oce_e(je,jk,jb) < land_boundary)THEN
          nabla2_vec_e(je,jk,jb) = p_patch_3D%wet_e(je,jk,jb)*&    !v_base%wet_e(je,jk,jb)*&
            &k_h(je,jk,jb)*(   &
            & p_patch%edges%system_orientation(je,jb) *     &
            & ( vort(ividx(je,jb,2),jk,ivblk(je,jb,2))     &
            & - vort(ividx(je,jb,1),jk,ivblk(je,jb,1)) )   &
            & * p_patch%edges%inv_primal_edge_length(je,jb) &
            & + &
            & ( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))     &
            & - z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )   &
            & * p_patch%edges%inv_dual_edge_length(je,jb))
        END DO
      END DO
    END DO

  END SUBROUTINE veloc_diff_harmonic_curl_curl
  !-------------------------------------------------------------------------
  !>
  !!  Computes  Bilaplacian of a vector field.The placement of the mixing tensor 
  !! follows the practice in NEMO. There is no corresponding operator in MPI-OM
  !!
  !! input:  lives on edges (velocity points)
  !! output: lives on edges
  !!
  !! @par Revision History
  !! Developed by P. Korn, MPI-M(2012)
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE veloc_diff_biharmonic_curl_curl( p_patch,p_patch_3D,u_vec_e,vort, k_h, p_op_coeff,&
                                            & p_vn_dual, nabla4_vec_e)
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
    REAL(wp), INTENT(in)                      :: u_vec_e(nproma,n_zlev,p_patch%nblks_e) 
    REAL(wp), INTENT(in)                      :: vort(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp), INTENT(in)                      :: k_h(:,:,:)
    TYPE(t_operator_coeff), INTENT(in)        :: p_op_coeff
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp), INTENT(inout)                   :: nabla4_vec_e(nproma,n_zlev,p_patch%nblks_e)

    !
    !Local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    REAL(wp) ::  z_div_c(nproma,n_zlev,p_patch%nblks_c)
    REAL(wp) ::  z_rot_v(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp) ::  z_nabla2_e(nproma,n_zlev,p_patch%nblks_e)
    REAL(wp) :: h_e(nproma,p_patch%nblks_e)
    TYPE(t_cartesian_coordinates)  :: p_nabla2_dual(nproma,n_zlev,p_patch%nblks_v)
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    edges_in_domain => p_patch%edges%in_domain
    !-----------------------------------------------------------------------

    slev = 1
    elev = n_zlev

    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    CALL veloc_diff_harmonic_curl_curl( u_vec_e, vort, p_patch,p_patch_3D, k_h, p_op_coeff,&
                                      & p_vn_dual, z_nabla2_e)
    ! compute divergence of vector field
!     CALL div_oce_3d( u_vec_e, p_patch, p_op_coeff%div_coeff, z_div_c)
!     ! DO jk = slev, elev
!     ! write(*,*)'vort1:',jk,maxval(vort(:,jk,:)),minval(vort(:,jk,:))
!     ! END DO
!     ! compute rotation of vector field for the ocean
!     !CALL rot_vertex_ocean_3D( p_patch, u_vec_e, p_vn_dual, p_op_coeff, z_rot_v)!
!     !z_rot_v=vort
!     !
!     !  loop through all patch edges (and blocks)
!     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
!       DO je = i_startidx, i_endidx
!         DO jk = slev, elev
!           !DO je = i_startidx, i_endidx
!           !IF(v_base%lsm_oce_e(je,jk,jb) < land_boundary)THEN
!           z_nabla2_e(je,jk,jb) =  &
!             & v_base%wet_e(je,jk,jb)*     &
!             & (p_patch%edges%system_orientation(je,jb) *  &
!             & ( vort(ividx(je,jb,2),jk,ivblk(je,jb,2))  &
!             & - vort(ividx(je,jb,1),jk,ivblk(je,jb,1)) )  &
!             & * p_patch%edges%inv_primal_edge_length(je,jb))  &
!             & +v_base%wet_e(je,jk,jb)*&
!             & (( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))    &
!             & - z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )  &
!             & * p_patch%edges%inv_dual_edge_length(je,jb))
!            z_nabla2_e(je,jk,jb)=sqrt(k_h(je,jk,jb))*z_nabla2_e(je,jk,jb)
!         END DO
!       END DO
!     END DO

    ! compute divergence of vector field
    CALL div_oce_3d( z_nabla2_e, p_patch, p_op_coeff%div_coeff, z_div_c)

    ! compute rotation of vector field for the ocean
    CALL map_edges2vert_3D( p_patch, &
                          & z_nabla2_e,&
                          & p_op_coeff%edge2vert_coeff_cc,&
                          & p_nabla2_dual)!h_e dummy, not used. Delete in sbr map_edges2vert
    CALL rot_vertex_ocean_3D( p_patch, p_patch_3D, z_nabla2_e, p_nabla2_dual, p_op_coeff, z_rot_v)!


    !combine divergence and vorticity
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO je = i_startidx, i_endidx
        DO jk = slev, elev

          nabla4_vec_e(je,jk,jb) =  &
            & p_patch_3D%wet_e(je,jk,jb)*                  &!v_base%wet_e(je,jk,jb)*     &
            & (p_patch%edges%system_orientation(je,jb) *  &
            & ( z_rot_v(ividx(je,jb,2),jk,ivblk(je,jb,2))  &
            & - z_rot_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) )  &
            & * p_patch%edges%inv_primal_edge_length(je,jb))  &
            & +p_patch_3D%wet_e(je,jk,jb)*                     &!v_base%wet_e(je,jk,jb)*&
            & (( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))    &
            & - z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )  &
            & * p_patch%edges%inv_dual_edge_length(je,jb))

            nabla4_vec_e(je,jk,jb)=-nabla4_vec_e(je,jk,jb)
        END DO
      END DO
    END DO

  END SUBROUTINE veloc_diff_biharmonic_curl_curl
! ! !-------------------------------------------------------------------------  
! ! !
! ! !  
! ! !>
! ! !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
! ! !!    implemented as Average(div( K_H grad RBF v)).
! ! !! 
! ! !! @par Revision History
! ! !! Developed  by  Peter Korn, MPI-M (2010).
! ! !! 
! ! !! mpi parallelized, sync required
! ! SUBROUTINE velocity_diffusion_horz_rbf(p_patch, p_param, p_diag, laplacian_vn_out)
! !   !
! !   TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! !   !REAL(wp), INTENT(inout)           :: vn_in(nproma,n_zlev,p_patch%nblks_e)
! !   TYPE(t_ho_params), INTENT(in)     :: p_param
! !   TYPE(t_hydro_ocean_diag)          :: p_diag
! !   REAL(wp), INTENT(out)             :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)
! !   !
! !   !Local variables
! !   INTEGER :: slev, elev
! !   INTEGER :: jk, jb, je!,jc
! !   INTEGER :: i_startidx_e, i_endidx_e
! !   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
! !   TYPE(t_subset_range), POINTER :: edges_in_domain
! !   ! CHARACTER(len=max_char_length), PARAMETER :: &
! !   !        & routine = ('mo_ocean_semi_implicit_ab:velocity_diffusion_horz')
! !   !-------------------------------------------------------------------------------
! !   edges_in_domain => p_patch%edges%in_domain
! !   !-------------------------------------------------------------------------------
! !   ! #slo# set intent out variable to zero due to nag -nan compiler-option
! !   laplacian_vn_out(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
! ! 
! !   slev = 1
! !   elev = n_zlev
! ! 
! !   !Step 2: Multiply laplacian with mixing coefficients
! !   DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !     CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! ! 
! !     DO jk = slev, elev
! !       DO je = i_startidx_e, i_endidx_e
! ! 
! !         il_c1 = p_patch%edges%cell_idx(je,jb,1)
! !         ib_c1 = p_patch%edges%cell_blk(je,jb,1)
! !         il_c2 = p_patch%edges%cell_idx(je,jb,2)
! !         ib_c2 = p_patch%edges%cell_blk(je,jb,2)
! ! 
! !         IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
! ! 
! !           laplacian_vn_out(je,jk,jb) = p_param%K_veloc_h(je,jk,jb)*laplacian_vn_out(je,jk,jb)
! !         ENDIF
! !       ENDDO
! !     END DO
! !   END DO
! !   CALL sync_patch_array(SYNC_E, p_patch,laplacian_vn_out)
! ! 
! ! END SUBROUTINE velocity_diffusion_horz_rbf
!-------------------------------------------------------------------------
!>
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
!!  mpi parallelized LL (no sync required)
SUBROUTINE velocity_diffusion_vert_mimetic( p_patch, p_patch_3D, p_diag, p_aux,p_op_coeff,&
                                           &p_param, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
 TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_hydro_ocean_aux)           :: p_aux
  TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
  !REAL(wp), INTENT(in)              :: h_c(nproma,p_patch%nblks_c) 
  TYPE(t_ho_params), INTENT(in)     :: p_param
  REAL(wp)                          :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)

  INTEGER                       :: slev, elev     ! vertical start and end level
  INTEGER                       :: jc, jk, jb, z_dolic
  INTEGER                       :: i_startidx, i_endidx
  TYPE(t_cartesian_coordinates) :: z_u(nproma,n_zlev+1,p_patch%nblks_c)!,  &
  TYPE(t_subset_range), POINTER :: cells_in_domain! , all_cells
  !  &                              z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:veloc diffusion vert mimetic')
  !-----------------------------------------------------------------------
!  all_cells => p_patch%cells%all
  cells_in_domain => p_patch%cells%in_domain
  !-----------------------------------------------------------------------

  z_u(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(1) = 0.0_wp
  z_u(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(2) = 0.0_wp
  z_u(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)%x(3) = 0.0_wp

  slev       = 1
  elev       = n_zlev

  !1 Vertical derivative of cell velocity vector times horizontal velocity
  ! loop runs now from slev to z_dolic:
  !DO jk = slev, elev
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
      IF (z_dolic >= MIN_DOLIC) THEN
        !1a) ocean surface 
        jk = slev
   ! #slo 2011-09-08 - include effect of surface boundary condition on G_n
   !                   eliminate in vertical diffusion: laplacian_vert
   !    &( &
        z_u(jc,jk,jb)%x =                                                &
        & (p_aux%bc_top_veloc_cc(jc,jb)%x - p_param%A_veloc_v(jc,jk+1,jb)&
        &*(p_diag%p_vn(jc,jk,jb)%x        - p_diag%p_vn(jc,jk+1,jb)%x)   &
        & /p_patch_3D%p_patch_1D(1)%del_zlev_i(jk+1))/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)
!        &/v_base%del_zlev_i(jk+1))/v_base%del_zlev_m(jk)
        !1b) ocean bottom 
        !ELSEIF ( jk == v_base%dolic_c(jc,jb) ) THEN
        !ELSEIF ( jk == z_dolic ) THEN
        jk = z_dolic
  !        IF(v_base%n_zlev>=2)THEN
        z_u(jc,jk,jb)%x =                                      &
        &  (p_param%A_veloc_v(jc,jk,jb)                        &
        & *(p_diag%p_vn(jc,jk-1,jb)%x-p_diag%p_vn(jc,jk,jb)%x) &
!        & / v_base%del_zlev_i(jk)                              &
        & / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)              &
        & - p_aux%bc_bot_veloc_cc(jc,jb)%x)/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)

  !        ENDIF
        !1c) ocean interior 
        !ELSEIF( jk>slev .AND. jk < v_base%dolic_c(jc,jb) ) THEN
        DO jk = slev+1, z_dolic-1
           z_u(jc,jk,jb)%x = &
           &   ( p_param%A_veloc_v(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x - &
           &                                  p_diag%p_vn(jc,jk,jb)%x)    &
           !&   /v_base%del_zlev_i(jk)                                     &
           & / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)              &
           &  - p_param%A_veloc_v(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x -  &
           &                                   p_diag%p_vn(jc,jk+1,jb)%x) &
           & / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk+1))/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)            
           !&   /v_base%del_zlev_i(jk+1))/v_base%del_zlev_m(jk)
        ! ENDIF  ! jk-condition
       !ENDIF    ! at least 2 vertical layers
        END DO ! jk ocean interior
        ENDIF  ! dolic>0
      END DO
    END DO
  !END DO

  CALL sync_patch_array(SYNC_C, p_patch, z_u(:,:,:)%x(1))
  CALL sync_patch_array(SYNC_C, p_patch, z_u(:,:,:)%x(2))
  CALL sync_patch_array(SYNC_C, p_patch, z_u(:,:,:)%x(3))  
  ! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)  
  CALL map_cell2edges_3D( p_patch,p_patch_3D, z_u,laplacian_vn_out,p_op_coeff)
  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=4  ! output print level (1-5, fix)
  CALL dbg_print('VelDiffMim: Laplacian'     ,laplacian_vn_out         ,str_module,idt_src)
  !---------------------------------------------------------------------

END SUBROUTINE velocity_diffusion_vert_mimetic
!-------------------------------------------------------------------------  


! ! !-------------------------------------------------------------------------  
! ! !>
! ! !!
! ! !! IMPORTANT: It is assumed that the velocity vector reconstruction from
! ! !! edges to cells has been done before.
! ! !!
! ! !! input:  lives on cells (velocity points)
! ! !! output: lives on edges (velocity points)
! ! !!
! ! !! @par Revision History
! ! !! Developed  by  Peter Korn, MPI-M (2010).
! ! !!
! ! !!  mpi parallelized, sync required
! ! SUBROUTINE velocity_diffusion_vert_rbf( p_patch,p_patch_3D, u_c, v_c, top_bc_u_c, top_bc_v_c,&
! !                           &  bot_bc_u_c,  bot_bc_v_c,p_param,p_op_coeff, laplacian_vn_out)
! !   TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! !   TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
! !   ! Components of cell based variable which is vertically diffused
! !   REAL(wp), INTENT(inout)           :: u_c(nproma,n_zlev,p_patch%nblks_c)
! !   REAL(wp), INTENT(inout)           :: v_c(nproma,n_zlev,p_patch%nblks_c)
! !   !REAL(wp), INTENT(in)              :: h_c(nproma,p_patch%nblks_c)
! !   ! Top boundary condition for cell based variables
! !   REAL(wp), INTENT(in)              :: top_bc_u_c(nproma,p_patch%nblks_c)
! !   REAL(wp), INTENT(in)              :: top_bc_v_c(nproma,p_patch%nblks_c)
! !   ! Bottom boundary condition for cell based variables
! !   REAL(wp), INTENT(in)              :: bot_bc_u_c(nproma,p_patch%nblks_c)
! !   REAL(wp), INTENT(in)              :: bot_bc_v_c(nproma,p_patch%nblks_c)
! !   TYPE(t_ho_params), INTENT(in)     :: p_param
! !   TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
! !   ! variable in which horizontally advected velocity is stored
! !   REAL(wp)            :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)
! ! 
! !   INTEGER :: slev, elev     ! vertical start and end level
! !   INTEGER :: jc, jk, jb
! !   INTEGER :: i_startidx_c, i_endidx_c
! !   TYPE(t_cartesian_coordinates) :: zu_cc(nproma,n_zlev,p_patch%nblks_c)
! !   TYPE(t_subset_range), POINTER :: all_cells
! !   REAL(wp) :: z_u(nproma,n_zlev+1,p_patch%nblks_c),  &
! !     &         z_v(nproma,n_zlev+1,p_patch%nblks_c)
! !   ! CHARACTER(len=max_char_length), PARAMETER :: &
! !   !        & routine = ('mo_oce_diffusion:veloc diffusion vert rbf')
! !   !-----------------------------------------------------------------------
! !   all_cells => p_patch%cells%all
! !   !-----------------------------------------------------------------------
! ! 
! !   ! #slo# set local variable to zero due to nag -nan compiler-option
! !   z_u(1:nproma,1:n_zlev+1,1:p_patch%nblks_c) = 0.0_wp
! !   z_v(1:nproma,1:n_zlev+1,1:p_patch%nblks_c) = 0.0_wp
! ! 
! !   slev       = 1
! !   elev       = n_zlev
! ! 
! ! 
! !   !1 Vertical derivative of cell velocity vector times horizontal velocity
! !   DO jk = slev, elev
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !       DO jc = i_startidx_c, i_endidx_c
! !       !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !       IF (p_patch_3D%lsm_oce_c(jc,jk,jb) <= sea_boundary) THEN
! !         !1a) 0cean surface
! !         !check if we have at least two layers of water
! !        !IF (v_base%dolic_c(jc,jb) >= 2) THEN
! !         IF(jk==slev)THEN
! !             ! u-component
! !             z_u(jc,jk,jb) =      &
! !             &( top_bc_u_c(jc,jb)       &
! !             & -p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk,jb)-u_c(jc,jk+1,jb)) &
! !             & /v_base%del_zlev_i(jk+1))/v_base%del_zlev_m(jk)
! !             ! v-component
! !             z_v(jc,jk,jb) =      &
! !             &( top_bc_v_c(jc,jb)       &
! !             & -p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
! !             &/v_base%del_zlev_i(jk+1)) &
! !             & /v_base%del_zlev_m(jk)
! !             !write(*,*)'RBD:vert-diff:top',jc,jk,jb, z_u_i(jc,jk,jb),top_bc_u_c(jc,jb)
! ! 
! !         !1b) ocean bottom 
! !         !ELSEIF ( jk == v_base%dolic_c(jc,jb) ) THEN
! !         ELSEIF ( jk==p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)) THEN
! !           IF(v_base%n_zlev>=2)THEN
! !             ! u-component
! !             z_u(jc,jk,jb)              &
! !             & = (p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk-1,jb)-u_c(jc,jk,jb))&
! !             & /v_base%del_zlev_i(jk)&
! !             & - bot_bc_u_c(jc,jb))/v_base%del_zlev_m(jk)
! !             ! v-component
! !             z_v(jc,jk,jb)&
! !             & = (p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk-1,jb)-v_c(jc,jk,jb))&
! !             & /v_base%del_zlev_i(jk)&
! !             & -bot_bc_v_c(jc,jb))/v_base%del_zlev_m(jk)
! !           !write(*,*)'u-diff botom:',jk,u_c(jc,elev-1,jb),u_c(jc,elev,jb),u_c(jc,elev,jb),bot_bc_u_c(jc,jb), z_u_i(jc,jk,jb)
! !           ENDIF
! ! 
! !         !1c) ocean interior 
! !         ELSEIF(jk>slev .AND. jk < p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) ) THEN
! !           ! u-component
! !           z_u(jc,jk,jb)&
! !           & = &
! !           & ( p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk-1,jb)-u_c(jc,jk,jb))&
! !           &/v_base%del_zlev_i(jk)&
! !           & - p_param%A_veloc_v(jc,jk+1,jb)  *(u_c(jc,jk,jb)-u_c(jc,jk+1,jb))&
! !           &/v_base%del_zlev_i(jk+1)) &
! !           & / v_base%del_zlev_m(jk)
! !           ! v-component
! !           z_v(jc,jk,jb)&
! !           & = &
! !           & ( p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk-1,jb)-v_c(jc,jk,jb))&
! !           &/v_base%del_zlev_i(jk)&
! !           & - p_param%A_veloc_v(jc,jk,jb)  *(v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
! !           &/v_base%del_zlev_i(jk+1)) &
! !           & / v_base%del_zlev_m(jk)
! ! 
! !         !1c) ocean bottom 
! !         ENDIF  ! jk-condition
! !         ENDIF
! !         CALL gvec2cvec( z_u(jc,jk,jb), z_v(jc,jk,jb),       &
! !         &               p_patch%cells%center(jc,jb)%lon,                &
! !         &               p_patch%cells%center(jc,jb)%lat,                &
! !         &               zu_cc(jc,jk,jb)%x(1),zu_cc(jc,jk,jb)%x(2),zu_cc(jc,jk,jb)%x(3))
! ! 
! !       !ENDIF   ! at least 2 vertical layers
! !       END DO
! !     END DO
! !   END DO
! ! 
! !   ! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)
! ! !  CALL map_cell2edges_3D( p_patch, zu_cc,laplacian_vn_out,p_op_coeff)
! !   CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)
! ! 
! !   !---------DEBUG DIAGNOSTICS-------------------------------------------
! !   idt_src=4  ! output print level (1-5, fix)
! !   CALL dbg_print('VelDiffRbf: Laplacian'     ,laplacian_vn_out         ,str_module,idt_src)
! !   !---------------------------------------------------------------------
! ! 
! ! END SUBROUTINE velocity_diffusion_vert_rbf
!------------------------------------------------------------------------
! mpi parallelized, sync required
SUBROUTINE tracer_diffusion_horz(p_patch,p_patch_3D, trac_in, p_os, K_T, diff_flx, subset_range)

  !
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
  !
  ! Patch on which computation is performed
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,p_patch%nblks_c)
  TYPE(t_hydro_ocean_state), TARGET :: p_os
  REAL(wp), INTENT(in)              :: K_T(:,:,:) !mixing coefficient for tracer
  REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
  !
  !
  !Local variables
  INTEGER                       :: slev, elev     ! vertical start and end level
  INTEGER                       :: jk, jb, je
  INTEGER                       :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER                       :: i_startidx_e, i_endidx_e
  REAL(wp)                      :: delta_z
  TYPE(t_subset_range), POINTER :: edges_in_domain
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
  !-------------------------------------------------------------------------------
  edges_in_domain => p_patch%edges%in_domain
  !-------------------------------------------------------------------------------
  slev = 1
  elev = n_zlev

  IF ( iswm_oce /= 1) THEN
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      jk=1
      DO je = i_startidx_e, i_endidx_e

        !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
        IF (p_patch_3D%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)

          delta_z = p_patch_3D%p_patch_1D(1)%prism_thick_e(je,jk,jb) !v_base%del_zlev_m(slev) + p_os%p_diag%h_e(je,jb)

          diff_flx(je,jk,jb) = K_T(je,jk,jb)*delta_z&
                     &*(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
                     &/p_patch%edges%dual_edge_length(je,jb)
          ENDIF
        ENDDO

      DO jk = slev+1, elev

        !delta_z  = v_base%del_zlev_m(jk)

        DO je = i_startidx_e, i_endidx_e

          !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          IF (p_patch_3D%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
            !Get indices of two adjacent triangles
            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)

            diff_flx(je,jk,jb) = K_T(je,jk,jb)                           &
                       & *p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(je,jk,jb)&!  p_os%p_diag%prism_thick_flat_sfc_e(je,jk,jb)&
                       &*(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
                       &/p_patch%edges%dual_edge_length(je,jb)
          ENDIF
        ENDDO
      END DO
    END DO
  ELSEIF ( iswm_oce == 1) THEN
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO jk = slev, elev
        DO je = i_startidx_e, i_endidx_e
          !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          IF (p_patch_3D%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN

          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)

          delta_z  = p_os%p_diag%thick_e(je,jb)

          diff_flx(je,jk,jb) = K_T(je,jk,jb)*delta_z*p_patch_3D%wet_e(je,jk,jb)&!v_base%wet_e(je,jk,jb) &
            &          *(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
            &          /p_patch%edges%dual_edge_length(je,jb)
          ENDIF
        ENDDO
      END DO
    END DO
  ENDIF

  ! Apply divergence to mixing times gradient to get laplacian
  !CALL div_oce( diff_flx, p_patch, laplacian_trac_out)
  IF (PRESENT(subset_range)) THEN
    IF (.NOT. subset_range%is_in_domain) &
     & CALL sync_patch_array(SYNC_E, p_patch, diff_flx)
  ENDIF
      

END SUBROUTINE tracer_diffusion_horz
!-------------------------------------------------------------------------  
!
!!Subroutine computes the vertical diffusive flux of an arbitrary tracer.
!>
!!
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
!! mpi parallelized, no sync required
SUBROUTINE tracer_diffusion_vert_expl(p_patch,p_patch_3D,        &
                                    & trac_c,          &
                                    & top_bc_tracer,   & 
                                    & A_v,             &
                                    & div_diff_flx)

  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  REAL(wp), INTENT(in)              :: trac_c(nproma, n_zlev,p_patch%nblks_c)
  !REAL(wp), INTENT(in)              :: dz(nproma, n_zlev+1,p_patch%nblks_c)
  REAL(wp), INTENT(in)              :: top_bc_tracer(nproma, p_patch%nblks_c)
  !REAL(wp), INTENT(in)              :: bot_bc_tracer(nproma, p_patch%nblks_c)
  REAL(wp), INTENT(inout)           :: A_v(:,:,:) 
  REAL(wp), INTENT(out)             :: div_diff_flx(nproma, n_zlev,p_patch%nblks_c)
  !
  !Local variables
  INTEGER                       :: slev, elev
  INTEGER                       :: jc, jk, jb
  INTEGER                       :: i_startidx_c, i_endidx_c
  INTEGER                       :: z_dolic
  ! vertical diffusive tracer flux
  REAL(wp)                      :: z_diff_flx(nproma, n_zlev+1,p_patch%nblks_c)
  TYPE(t_subset_range), POINTER :: all_cells
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:tracer_diffusion_vert')
  !-----------------------------------------------------------------------
  all_cells => p_patch%cells%all
  !-----------------------------------------------------------------------
  slev              = 1
  elev              = n_zlev
  z_diff_flx(1:nproma, 1:n_zlev+1,1:p_patch%nblks_c) = 0.0_wp

  !1 Vertical derivative of tracer
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      z_dolic  = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
      IF ( z_dolic >=MIN_DOLIC ) THEN
        !1a) 0cean surface
        !diff_flx(jc,slev,jb) = A_v(slev)&
        !& *(top_bc_tracer(jc,jb)-trac_c(jc,slev,jb))/v_base%del_zlev_i(slev)
        z_diff_flx(jc,slev,jb) = top_bc_tracer(jc,jb)

        !1b) ocean interior 
        DO jk = slev+1, z_dolic
          !IF(dz(jc,jk,jb)/=0.0_wp)THEN
            ! #slo# 2011-09-15 - Corrections:
            !   - must be checked at all vertical processes in the model
            !   - correction only active for variable level thicknesses (distances)
            z_diff_flx(jc,jk,jb)&
            & = A_v(jc,jk,jb) &
            & * (trac_c(jc,jk-1,jb)-trac_c(jc,jk,jb))/ p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)!v_base%del_zlev_i(jk)
        END DO


        DO jk = 1, z_dolic
          ! positive vertical divergence in direction of w (upward positive)
           div_diff_flx(jc,jk,jb) = (z_diff_flx(jc,jk,jb) - z_diff_flx(jc,jk+1,jb))&
                                &/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)!v_base%del_zlev_m(jk) 
          END DO
        !1c) ocean bottom zero bottom boundary condition
        !diff_flx(jc,z_dolic+1,jb) = bot_bc_tracer(jc,jk)
      ENDIF
    END DO
  END DO

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=5  ! output print level (1-5, fix)
  CALL dbg_print('TrcDiffExpl: z_diff_flx'   ,z_diff_flx               ,str_module,idt_src)
  CALL dbg_print('TrcDiffExpl: div_diff_flx' ,div_diff_flx             ,str_module,idt_src)
  !---------------------------------------------------------------------

END SUBROUTINE tracer_diffusion_vert_expl
! !-------------------------------------------------------------------------  
! !
! !!Subroutine implements implicit vertical diffusion for scalar fields.
! !>
! !! sbr identical to sbr above but now with homogeneous boundary conditions
! !!
! !! @par Revision History
! !! Developed  by  Peter Korn, MPI-M (2011).
! !!
! !! mpi parallelized, no sync required
! SUBROUTINE tracer_diffusion_vert_impl_hom( p_patch,   &
!                                 & field_column,   &
!                                 & h_c,            &
!                                 & A_v,            &
!                                 & diff_column)
! 
!   TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!   REAL(wp), INTENT(inout)           :: field_column(:,:,:)
!   REAL(wp), INTENT(IN)              :: h_c(:,:)           !surface height, relevant for thickness of first cell 
!   REAL(wp), INTENT(in)              :: A_v(:,:,:)
!   REAL(wp), INTENT(out)             :: diff_column(:,:,:)
!   !
!   !Local variables
!   INTEGER :: slev
!   INTEGER :: jc, jk, jb
!   INTEGER :: i_startidx_c, i_endidx_c
!   REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
!   REAL(wp) :: z_tmp
!   REAL(wp) :: inv_zinv_i(1:n_zlev)
!   REAL(wp) :: inv_zinv_m(1:n_zlev)
!   !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
!   !REAL(wp) :: z_c1(nproma,1,p_patch%nblks_c)
!   INTEGER  :: z_dolic
!   TYPE(t_subset_range), POINTER :: all_cells
!   ! CHARACTER(len=max_char_length), PARAMETER :: &
!   !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!   !-----------------------------------------------------------------------
!   slev    = 1
!   !A_v    = 0.0001_wp
! 
!   diff_column(:,:,:)      = field_column(:,:,:)
!   a(slev:n_zlev)          = 0.0_wp
!   b(slev:n_zlev)          = 0.0_wp
!   c(slev:n_zlev)          = 0.0_wp
!   !bet(slev:n_zlev)       = 1.0_wp
!   !gam(slev:n_zlev)       = 0.0_wp
!   inv_zinv_i(slev:n_zlev) = 0.0_wp
!   inv_zinv_m(slev:n_zlev) = 0.0_wp
! 
!   all_cells => p_patch%cells%all
! 
!   DO jb = all_cells%start_block, all_cells%end_block
!     CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!     DO jc = i_startidx_c, i_endidx_c
!       z_dolic = v_base%dolic_c(jc,jb)
! 
!       IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
!         IF ( z_dolic >=MIN_DOLIC ) THEN
! 
!           inv_zinv_i(:) = 1.0_wp/v_base%del_zlev_i(:)
!           inv_zinv_m(:) = 1.0_wp/v_base%del_zlev_m(:)
! 
!            !inv_zinv_i(1) = 1.0_wp/(v_base%del_zlev_i(1)+h_c(jc,jb))
!            !inv_zinv_m(1) = 1.0_wp/(v_base%del_zlev_m(1)+h_c(jc,jb))
!           !inv_zinv_i(:) = 1.0_wp/h_c(jc,:,jb)
!           !inv_zinv_m(:) = 1.0_wp/h_c(jc,:,jb)
! 
!           !Fill triangular matrix
!           !b is diagonal a and c are upper and lower band
!           DO jk = slev+1, z_dolic-1
!             a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)*dtime
!             c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)*dtime
!             b(jk) = 1.0_wp-a(jk)-c(jk)
!           END DO
! 
!           ! The first row
!           c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)*dtime
!           a(slev) = 0.0_wp           
!           b(slev) = 1.0_wp- c(slev) - a(slev) 
! 
!           ! The last row
!           a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)*dtime
!           c(z_dolic) = 0.0_wp
!           b(z_dolic) = 1.0_wp - a(z_dolic) - c(z_dolic)
! 
!           DO jk=slev, z_dolic-1
!             IF(b(jk)/=0.0_wp)THEN
!               a(jk) = a(jk)/b(jk)
!               c(jk) = c(jk)/b(jk)
!               field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
!               b(jk)=1.0_wp
!             ENDIF
!           END DO
! 
!           DO jk=slev+1, z_dolic-1
!             b(jk)                  = b(jk)-a(jk)*c(jk-1)
!             field_column(jc,jk,jb) = field_column(jc,jk,jb)&
!                           &-a(jk)*field_column(jc,jk-1,jb)
!             c(jk)                  = c(jk)/b(jk)
!             field_column(jc,jk,jb) = field_column(jc,jk,jb)/b(jk)
!             b(jk)                  = 1.0_wp
!           END DO
! 
!           z_tmp = b(z_dolic)-a(z_dolic)*c(z_dolic-1)
!           z_tmp = (field_column(jc,z_dolic,jb)-a(z_dolic)*field_column(jc,z_dolic-1,jb))/z_tmp
! 
!           field_column(jc,z_dolic,jb) = z_tmp
!           DO jk = z_dolic-1,1,-1
!             field_column(jc,jk,jb) = field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
!           END DO
!           DO jk = 1,z_dolic!-1
!             diff_column(jc,jk,jb) = field_column(jc,jk,jb)
!           END DO
!         ELSEIF ( z_dolic < MIN_DOLIC ) THEN
!           diff_column(jc,:,jb) = 0.0_wp
!           field_column(jc,:,jb)= 0.0_wp
!         ENDIF
!       ELSEIF( v_base%lsm_oce_c(jc,1,jb) > sea_boundary ) THEN
!         diff_column(jc,:,jb) = field_column(jc,:,jb)
!       ENDIF
! 
! 
!     END DO
!   END DO
! 
! END SUBROUTINE tracer_diffusion_vert_impl_hom
!-------------------------------------------------------------------------  
!
!!Subroutine implements implicit vertical diffusion for scalar fields.
!>
!! sbr identical to sbr above but now with homogeneous boundary conditions
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
!! mpi parallelized, no sync required
SUBROUTINE tracer_diffusion_vert_impl_hom( p_patch,p_patch_3D,     &
                                         & field_column,&
                                         & h_c,         &
                                         & A_v,         &
                                         & p_op_coeff,  &
                                         & diff_column)

  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  REAL(wp), INTENT(inout)           :: field_column(1:nproma,1:n_zlev,1:p_patch%nblks_c)
  REAL(wp), INTENT(IN)              :: h_c(1:nproma,1:p_patch%nblks_c)  !surface height, relevant for thickness of first cell 
  REAL(wp), INTENT(in)              :: A_v(:,:,:)  
  TYPE(t_operator_coeff),TARGET     :: p_op_coeff
  REAL(wp), INTENT(inout)           :: diff_column(1:nproma,1:n_zlev,1:p_patch%nblks_c)
  !
  !Local variables
  INTEGER :: slev
  INTEGER :: jc, jk, jb
  INTEGER :: i_startidx_c, i_endidx_c
  !REAL(wp) :: a0(1:n_zlev), b0(1:n_zlev), c0(1:n_zlev)
  REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
  REAL(wp) :: z_tmp
  !REAL(wp) :: inv_zinv_i(1:n_zlev)
  !REAL(wp) :: inv_zinv_m(1:n_zlev)
  REAL(wp) :: dt_inv
  !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
  !REAL(wp) :: z_c1(nproma,1,p_patch%nblks_c)
  INTEGER  :: z_dolic
  TYPE(t_subset_range), POINTER :: all_cells
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
  !-----------------------------------------------------------------------
  all_cells => p_patch%cells%all
  !-----------------------------------------------------------------------
  slev    = 1
  !A_v    = 0.0001_wp
  dt_inv = 1.0_wp/dtime

  field_column(1:nproma,1:n_zlev,1:p_patch%nblks_c)&
  &=field_column(1:nproma,1:n_zlev,1:p_patch%nblks_c)*dt_inv

  CALL sync_patch_array(SYNC_C, p_patch, field_column)

  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)  !v_base%dolic_c(jc,jb)

      !IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
      IF (p_patch_3D%lsm_oce_c(jc,1,jb) <= sea_boundary) THEN
        IF ( z_dolic >=MIN_DOLIC ) THEN

          !inv_zinv_i(:) = 1.0_wp/v_base%del_zlev_i(:)
          !inv_zinv_m(:) = 1.0_wp/v_base%del_zlev_m(:)
           !inv_zinv_i(1) = 1.0_wp/(v_base%del_zlev_i(1)+h_c(jc,jb))
           !inv_zinv_m(1) = 1.0_wp/(v_base%del_zlev_m(1)+h_c(jc,jb))
          !inv_zinv_i(:) = 1.0_wp/h_c(jc,:,jb)
          !inv_zinv_m(:) = 1.0_wp/h_c(jc,:,jb)

          !Fill triangular matrix
          !b is diagonal a and c are upper and lower band
          !DO jk = slev+1, z_dolic-1
          !  a0(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)!*dtime
          !  c0(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)!*dtime
          !  b0(jk) = dt_inv-a0(jk)-c0(jk)
          !END DO
         !! The first row
          !c0(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)!*dtime
          !a0(slev) = 0.0_wp           
          !b0(slev) = dt_inv- c0(slev)! - a(slev) 
          !! The last row
          !a0(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)!*dtime
          !c0(z_dolic) = 0.0_wp
          !b0(z_dolic) = dt_inv - a0(z_dolic)! - c(z_dolic)

! ! ! !------------coefficient test--------------------------
! DO jk = slev, z_dolic
! IF(a0(jk)/=p_op_coeff%matrix_vert_diff_c(jc,jk,jb,1))THEN
! write(1239,*)'a wrong',jc,jk,jb,a0(jk), p_op_coeff%matrix_vert_diff_c(jc,jk,jb,1)
! !STOP
! ENDIF
! IF(b0(jk)/=p_op_coeff%matrix_vert_diff_c(jc,jk,jb,2))THEN
! write(1239,*)'b wrong',jc,jk,jb,b0(jk), p_op_coeff%matrix_vert_diff_c(jc,jk,jb,2)
! !STOP
! ENDIF
! IF(c0(jk)/=p_op_coeff%matrix_vert_diff_c(jc,jk,jb,3))THEN
! write(1239,*)'c wrong',jc,jk,jb,c0(jk),p_op_coeff%matrix_vert_diff_c(jc,jk,jb,3)
! !STOP
! ENDIF
! END DO
! ! ! !------------coefficient test--------------------------
           !a => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,1)
           !b => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,2)
           !c => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,3)
          !field_column(jc,1:z_dolic,jb)=field_column(jc,1:z_dolic,jb)*dt_inv
           a(1:z_dolic) = p_op_coeff%matrix_vert_diff_c(jc,1:z_dolic,jb,1)
           b(1:z_dolic) = p_op_coeff%matrix_vert_diff_c(jc,1:z_dolic,jb,2)
           c(1:z_dolic) = p_op_coeff%matrix_vert_diff_c(jc,1:z_dolic,jb,3)


          DO jk=slev, z_dolic-1
            IF(b(jk)/=0.0_wp)THEN
              a(jk) = a(jk)/b(jk)
              c(jk) = c(jk)/b(jk)
              field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
              b(jk)=1.0_wp
            ENDIF
          END DO

          !Apply the matrix
          DO jk=slev+1, z_dolic-1
            b(jk)                  = b(jk)-a(jk)*c(jk-1)
            field_column(jc,jk,jb) = field_column(jc,jk,jb)&
                          &-a(jk)*field_column(jc,jk-1,jb)
            c(jk)                  = c(jk)/b(jk)
            field_column(jc,jk,jb) = field_column(jc,jk,jb)/b(jk)
            b(jk)                  = 1.0_wp
          END DO

          z_tmp = b(z_dolic)-a(z_dolic)*c(z_dolic-1)
          z_tmp = (field_column(jc,z_dolic,jb)-a(z_dolic)*field_column(jc,z_dolic-1,jb))/z_tmp

          field_column(jc,z_dolic,jb) = z_tmp
          DO jk = z_dolic-1,1,-1
            field_column(jc,jk,jb) = field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
          END DO
           DO jk = 1,z_dolic!-1
             diff_column(jc,jk,jb) = field_column(jc,jk,jb)!*dtime
           END DO
        !ELSEIF ( z_dolic < MIN_DOLIC ) THEN
        !  diff_column(jc,:,jb) = 0.0_wp
        !  field_column(jc,:,jb)= 0.0_wp
        ENDIF
      !ELSEIF( v_base%lsm_oce_c(jc,1,jb) > sea_boundary ) THEN
      !  diff_column(jc,1:z_dolic,jb) = field_column(jc,1:z_dolic,jb)
      ENDIF
    END DO
  END DO

  CALL sync_patch_array(SYNC_C, p_patch, diff_column)

END SUBROUTINE tracer_diffusion_vert_impl_hom
!-------------------------------------------------------------------------  
!
!!Subroutine implements implicit vertical diffusion for horizontal velocity fields
!!by inverting a scalar field..
!>
!! sbr identical to previous one, except for homogeneous boundary conditions
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!  mpi parallelized LL (no sync required)
SUBROUTINE veloc_diffusion_vert_impl_hom( p_patch, p_patch_3D,      &
                                        & field_column,  &
                                        & h_e,           &
                                        & A_v,           &
                                        & p_op_coeff,    &
                                        & diff_column)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch 
  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  REAL(wp), INTENT(inout)           :: field_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)
  !surface height at edges, relevant for thickness of first cell 
  REAL(wp), INTENT(IN)              :: h_e(1:nproma,1:p_patch%nblks_e)
  REAL(wp), INTENT(inout)           :: A_v(:,:,:)   
  TYPE(t_operator_coeff), TARGET    :: p_op_coeff
  REAL(wp), INTENT(out)             :: diff_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)
  !
  !Local variables
  INTEGER :: slev
  INTEGER :: je, jk, jb
  INTEGER :: i_startidx, i_endidx
!   REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
  REAL(wp),POINTER :: a(:), b(:), c(:)
  !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
  REAL(wp) :: z_tmp
  REAL(wp) :: dt_inv
  INTEGER  :: z_dolic
  !REAL(wp) :: inv_zinv_i(1:n_zlev)
  !REAL(wp) :: inv_zinv_m(1:n_zlev)
  TYPE(t_subset_range), POINTER :: all_edges
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
  !-----------------------------------------------------------------------
  all_edges => p_patch%edges%all
  !-----------------------------------------------------------------------
  slev = 1
  dt_inv=1.0_wp/dtime

  !diff_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)= 0.0_wp

  field_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)&
  &=field_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)*dt_inv


  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=5  ! output print level (1-5, fix)
  CALL dbg_print('VelDifImplHomIn:field_col' ,field_column             ,str_module,idt_src)
  CALL dbg_print('VelDifImplHomIn: diff_col' ,diff_column              ,str_module,idt_src)
  !---------------------------------------------------------------------

  DO jb = all_edges%start_block, all_edges%end_block
    CALL get_index_range(all_edges, jb, i_startidx, i_endidx)
    DO je = i_startidx, i_endidx
      z_dolic = p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)!!v_base%dolic_e(je,jb)

      !IF ( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
      IF (p_patch_3D%lsm_oce_e(je,1,jb) <= sea_boundary) THEN
        IF ( z_dolic >= MIN_DOLIC ) THEN

          !inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
          !inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)
          !!Fill triangular matrix
          !!b is diagonal a and c are upper and lower band
          !DO jk = slev+1, z_dolic-1
          !  a(jk) = -A_v(je,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)!*dtime
          !  c(jk) = -A_v(je,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)!*dtime
          !  b(jk) = dt_inv-a(jk)-c(jk)
          !END DO
          !! The first row
          ! c(slev) = -A_v(je,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)!*dtime
          ! a(slev) = 0.0_wp           
          ! b(slev) = dt_inv- c(slev) !- a(slev) 
          !! The last row
          !a(z_dolic) = -A_v(je,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)!*dtime
          !c(z_dolic) = 0.0_wp
          !b(z_dolic) = dt_inv - a(z_dolic)! - c(z_dolic)

           a => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,1)
           b => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,2)
           c => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,3)
           !a(1:z_dolic) = p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,1)
           !b(1:z_dolic) = p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,2)
           !c(1:z_dolic) = p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,3)


          !Apply the matrix
          DO jk=slev, z_dolic-1
            IF(b(jk)/=0.0_wp)THEN
              a(jk) = a(jk)/b(jk)
              c(jk) = c(jk)/b(jk)
              field_column(je,jk,jb)=field_column(je,jk,jb)/b(jk)
              b(jk)=1.0_wp
            ENDIF
          END DO

          DO jk=slev+1, z_dolic-1
            b(jk)                  = b(jk)-a(jk)*c(jk-1)
            field_column(je,jk,jb) = field_column(je,jk,jb)&
                          &-a(jk)*field_column(je,jk-1,jb)
            c(jk)                  = c(jk)/b(jk)
            field_column(je,jk,jb) = field_column(je,jk,jb)/b(jk)
            b(jk)                  = 1.0_wp
          END DO

          z_tmp = b(z_dolic)-a(z_dolic)*c(z_dolic-1)
          z_tmp = (field_column(je,z_dolic,jb)-a(z_dolic)*field_column(je,z_dolic-1,jb))/z_tmp

          field_column(je,z_dolic,jb) = z_tmp
          DO jk = z_dolic-1,1,-1
            field_column(je,jk,jb) = field_column(je,jk,jb)-c(jk)*field_column(je,jk+1,jb)
          END DO
          DO jk = 1,z_dolic!-1
            diff_column(je,jk,jb) = field_column(je,jk,jb)
          END DO
        !ELSEIF ( z_dolic < MIN_DOLIC ) THEN
        !  diff_column(je,1:z_dolic,jb) = 0.0_wp
        !  field_column(je,1:z_dolic,jb)= 0.0_wp
        ENDIF
      !ELSEIF( v_base%lsm_oce_e(je,1,jb) > sea_boundary ) THEN
      !  diff_column(je,1:z_dolic,jb) = field_column(je,1:z_dolic,jb)
       !diff_column(je,:,jb) = 0.0_wp
       !field_column(je,:,jb)= 0.0_wp
      ENDIF
    END DO
  END DO

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=5  ! output print level (1-5, fix)
  CALL dbg_print('VelDifImplHom: field_col'  ,field_column             ,str_module,idt_src)
  CALL dbg_print('VelDifImplHom: diff_col'   ,diff_column              ,str_module,idt_src)
  !---------------------------------------------------------------------

END SUBROUTINE veloc_diffusion_vert_impl_hom
!------------------------------------------------------------------------
END MODULE mo_oce_diffusion
