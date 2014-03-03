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
USE mo_ocean_nml,           ONLY: n_zlev, iswm_oce, veloc_diffusion_order, veloc_diffusion_form, &
  & use_tracer_x_height
USE mo_run_config,          ONLY: dtime
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_oce_types,           ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, t_ocean_tracer, t_hydro_ocean_aux
USE mo_model_domain,        ONLY: t_patch, t_patch_3D
!USE mo_exception,           ONLY: message, finish!, message_text
USE mo_oce_physics,         ONLY: t_ho_params
USE mo_scalar_product,      ONLY: map_cell2edges_3D
USE mo_oce_math_operators,  ONLY: div_oce_3D, rot_vertex_ocean_3d,&
 &                                map_edges2vert_3D
USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
USE mo_exception,           ONLY: finish !, message_text, message

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
PUBLIC :: velocity_diffusion_vert_explicit
PUBLIC :: velocity_diffusion_vertical_implicit
PUBLIC :: tracer_diffusion_horz
PUBLIC :: tracer_diffusion_vert_explicit
PUBLIC :: tracer_diffusion_vertical_implicit

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
SUBROUTINE velocity_diffusion( p_patch_3D, vn_in, p_param, p_diag,p_op_coeff, laplacian_vn_out)

 TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  REAL(wp), INTENT(in)                :: vn_in(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_ho_params), INTENT(in)       :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag),INTENT(in) :: p_diag
  TYPE(t_operator_coeff),INTENT(in)   :: p_op_coeff
  REAL(wp), INTENT(INOUT)             :: laplacian_vn_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

  !Local variables
  !REAL(wp) :: z_lapl(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  !INTEGER  :: jk
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        

  IF(veloc_diffusion_order==1)THEN
     
 !divgrad laplacian is chosen 
    IF(veloc_diffusion_form==2)THEN

      CALL veloc_diff_harmonic_div_grad( p_patch_3D,&
                                       & p_param,   &
                                       & p_diag,    &
                                       & p_op_coeff,&
                                       & laplacian_vn_out)
      
    ELSEIF(veloc_diffusion_form==1)THEN

       CALL veloc_diff_harmonic_curl_curl( vn_in, p_diag%vort,&
                                         & p_patch_3D,           &
                                         & p_op_coeff,        &
                                         & p_diag%p_vn_dual,  &
                                         & laplacian_vn_out,  & 
                                         & p_param%K_veloc_h)
      
   ENDIF

  ELSEIF(veloc_diffusion_order==2)THEN

   IF(veloc_diffusion_form==2)THEN

      CALL veloc_diff_biharmonic_div_grad( p_patch_3D,   &
                                         & p_param,      &
                                         & p_diag,       &
                                         & p_op_coeff,   &
                                         & laplacian_vn_out)
   ELSEIF(veloc_diffusion_form==1)THEN

       CALL veloc_diff_biharmonic_curl_curl( p_patch_3D,        &
                                           & vn_in,             &
                                           & p_diag%vort,       &
                                           & p_param%K_veloc_h, &
                                           & p_op_coeff,        &
                                           & p_diag%p_vn_dual,  &
                                           & laplacian_vn_out)
   ENDIF
 ENDIF

 CALL sync_patch_array(SYNC_E, p_patch_3D%p_patch_2D(1), laplacian_vn_out)


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
SUBROUTINE veloc_diff_harmonic_div_grad( p_patch_3D, p_param, p_diag,&
                                       & p_op_coeff, laplacian_vn_out)
 
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
  REAL(wp), INTENT(INOUT)           :: laplacian_vn_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

  !Local variables
  INTEGER :: slev, elev
  INTEGER :: jk, jb, je,jc
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: start_index, end_index
  INTEGER :: i_startidx_e, i_endidx_e
  INTEGER :: idx_cartesian
  INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  TYPE(t_cartesian_coordinates) :: z_grad_u    (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
  TYPE(t_patch), POINTER        :: p_patch 
!-------------------------------------------------------------------------------
  p_patch         => p_patch_3D%p_patch_2D(1)
  all_cells       => p_patch%cells%all
  all_edges       => p_patch%edges%all
  edges_in_domain => p_patch%edges%in_domain
!-------------------------------------------------------------------------------
  slev = 1
  elev = n_zlev

  laplacian_vn_out(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e) = 0.0_wp

  ! loop over cells in local domain + halo
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, start_index, end_index)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = start_index, end_index
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
        !IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
        IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
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
    CALL get_index_range(all_cells, jb, start_index, end_index)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = start_index, end_index

         !IF ( v_base%lsm_c(jc,jk,jb) >= boundary ) THEN
         IF (p_patch_3D%lsm_c(jc,jk,jb) >= boundary) THEN
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
  CALL map_cell2edges_3D( p_patch_3D, z_div_grad_u, laplacian_vn_out, p_op_coeff)
  
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
SUBROUTINE veloc_diff_biharmonic_div_grad( p_patch_3D, p_param, p_diag,&
                                         & p_op_coeff, laplacian_vn_out)

 TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  !REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
  REAL(wp), INTENT(INOUT)           :: laplacian_vn_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

  !Local variables
  INTEGER :: slev, elev
  INTEGER :: jk, jb, je,jc
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: start_index, end_index
  INTEGER :: i_startidx_e, i_endidx_e
  INTEGER :: idx_cartesian
  INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  REAL(wp)                      :: z_laplacian_vn_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_cartesian_coordinates) :: z_grad_u          (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_cartesian_coordinates) :: z_div_grad_u      (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
  TYPE(t_patch), POINTER        :: p_patch 
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
!------------------------------------------------------------------------------- 
  p_patch         => p_patch_3D%p_patch_2D(1)
  all_cells       => p_patch%cells%all
  all_edges       => p_patch%edges%all
  edges_in_domain => p_patch%edges%in_domain
!-------------------------------------------------------------------------------
  slev = 1
  elev = n_zlev

  laplacian_vn_out  (1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e) = 0.0_wp
  z_laplacian_vn_out(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e) = 0.0_wp

  ! loop over cells in local domain + halo
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, start_index, end_index)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = start_index, end_index
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

      IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        z_grad_u(je,jk,jb)%x = (p_diag%p_vn(il_c2,jk,ib_c2)%x &
          &                   - p_diag%p_vn(il_c1,jk,ib_c1)%x)&
          &                   / p_patch%edges%dual_edge_length(je,jb)
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
    CALL get_index_range(all_cells, jb, start_index, end_index)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = start_index, end_index

         IF (p_patch_3D%lsm_c(jc,jk,jb) >= boundary) THEN
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

      IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
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
    CALL get_index_range(all_cells, jb, start_index, end_index)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = start_index, end_index

         IF (p_patch_3D%lsm_c(jc,jk,jb) >= boundary) THEN
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
  CALL map_cell2edges_3D( p_patch_3D, z_div_grad_u,laplacian_vn_out,p_op_coeff)!,subset_range=all_cells)

!   !laplacian_vn_out=-laplacian_vn_out
!   DO jk=1,n_zlev
!    write(*,*)'Bi-lapla',jk,maxval(laplacian_vn_out(:,jk,:)),&
!    &minval(laplacian_vn_out(:,jk,:))
!   END DO
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
  SUBROUTINE veloc_diff_harmonic_curl_curl( u_vec_e, vort, p_patch_3D, p_op_coeff,&
    & p_vn_dual, nabla2_vec_e,k_h )
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)  :: p_patch_3D
    REAL(wp), INTENT(in)                      :: u_vec_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e) 
    REAL(wp), INTENT(in)                      :: vort   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff), INTENT(in)        :: p_op_coeff
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual    (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(inout)                   :: nabla2_vec_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp),OPTIONAL, INTENT(in)             :: k_h(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)    
!
    !Local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    REAL(wp) ::  z_div_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)!, &
    !REAL(wp) ::  z_vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    !REAL(wp) ::  z_rot_v(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    p_patch         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => p_patch%edges%in_domain
    !-----------------------------------------------------------------------

    slev = 1
    elev = n_zlev

    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    ! compute divergence of vector field
    z_div_c(:,:,p_patch%alloc_cell_blocks) = 0.0_wp
    CALL div_oce_3d( u_vec_e, p_patch, p_op_coeff%div_coeff, z_div_c)
    CALL sync_patch_array(SYNC_C,p_patch,z_div_c)

    ! compute rotation of vector field for the ocean
    !CALL rot_vertex_ocean_3D( p_patch, u_vec_e, p_vn_dual, p_op_coeff, z_rot_v)!
    !CALL sync_patch_array(SYNC_V,p_patch,z_rot_v)
    !z_rot_v=vort
    nabla2_vec_e(:,:,:) = 0.0_wp

    IF(PRESENT(k_h))THEN
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO je = i_startidx, i_endidx
        ! DO jk = slev, elev
        DO jk = slev, p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)

!              write(0, *) "==============================="
!              write(0, *) "0",  je,jk,jb
!              write(0, *) "1",  p_patch_3D%wet_e(je,jk,jb)
!              write(0,*)  "2",  k_h(je,jk,jb)
!              write(0,*)  "3",  p_patch%edges%system_orientation(je,jb)
!              write(0,*)  "4",  vort(ividx(je,jb,2),jk,ivblk(je,jb,2))
!              write(0,*)  "5",  vort(ividx(je,jb,1),jk,ivblk(je,jb,1))
!              write(0,*)  "6",  p_patch%edges%inv_primal_edge_length(je,jb)
!              write(0,*)  "7",  z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))
!              write(0,*)  "8",  z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1))
!              write(0,*)  "9",  p_patch%edges%inv_dual_edge_length(je,jb)
          !IF(v_base%lsm_e(je,jk,jb) < land_boundary)THEN

          nabla2_vec_e(je,jk,jb) = p_patch_3D%wet_e(je,jk,jb)* &    !v_base%wet_e(je,jk,jb)*&
            & k_h(je,jk,jb) * (                                &
            & p_patch%edges%system_orientation(je,jb) *        &
            & ( vort(ividx(je,jb,2),jk,ivblk(je,jb,2))         &
            & - vort(ividx(je,jb,1),jk,ivblk(je,jb,1)) )       &
            & * p_patch%edges%inv_primal_edge_length(je,jb)    &
            & +                                                &
            & ( z_div_c(icidx(je,jb,2),jk,icblk(je,jb,2))      &
            & - z_div_c(icidx(je,jb,1),jk,icblk(je,jb,1)) )    &
            & * p_patch%edges%inv_dual_edge_length(je,jb))
        END DO
      END DO
    END DO
    ELSEIF(.NOT.PRESENT(k_h))THEN
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
          nabla2_vec_e(je,jk,jb) = p_patch_3D%wet_e(je,jk,jb)*&
            &(   &
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

    ENDIF
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
  SUBROUTINE veloc_diff_biharmonic_curl_curl( p_patch_3D,u_vec_e,vort, k_h, p_op_coeff,&
                                            & p_vn_dual, nabla4_vec_e)
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)  :: p_patch_3D
    REAL(wp), INTENT(in)                      :: u_vec_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e) 
    REAL(wp), INTENT(in)                      :: vort   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(in)                      :: k_h(:,:,:)
    TYPE(t_operator_coeff), INTENT(in)        :: p_op_coeff
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(inout)                   :: nabla4_vec_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !
    !Local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: i_startidx, i_endidx
    REAL(wp) ::  z_div_c   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) ::  z_rot_v   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp) ::  z_nabla2_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: h_e        (nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates)  :: p_nabla2_dual(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER,  DIMENSION(:,:,:), POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    p_patch         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => p_patch%edges%in_domain
    !-----------------------------------------------------------------------

    slev = 1
    elev = n_zlev

    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    z_nabla2_e(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
!TODO review: k_h missing in current version but used in oce-struct-dev branch
    CALL veloc_diff_harmonic_curl_curl( u_vec_e, vort, p_patch_3D, p_op_coeff,&
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
!           !IF(v_base%lsm_e(je,jk,jb) < land_boundary)THEN
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
    CALL map_edges2vert_3d( p_patch, &
                          & z_nabla2_e,&
                          & p_op_coeff%edge2vert_coeff_cc,&
                          & p_nabla2_dual)!h_e dummy, not used. Delete in sbr map_edges2vert
    CALL rot_vertex_ocean_3d( p_patch_3D, z_nabla2_e, p_nabla2_dual, p_op_coeff, z_rot_v)!


    !combine divergence and vorticity
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO je = i_startidx, i_endidx
        DO jk = slev, elev

          nabla4_vec_e(je,jk,jb) =  &
            & p_patch_3D%wet_e(je,jk,jb)*     &
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
SUBROUTINE velocity_diffusion_vert_explicit( p_patch_3D, p_diag, p_aux,p_op_coeff,&
                                           &p_param, laplacian_vn_out)

 TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_hydro_ocean_aux)           :: p_aux
  TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
  TYPE(t_ho_params), INTENT(in)     :: p_param
  REAL(wp)                          :: laplacian_vn_out(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

  INTEGER                       :: slev, elev     ! vertical start and end level
  INTEGER                       :: jc, jk, jb, z_dolic
  INTEGER                       :: i_startidx, i_endidx
  TYPE(t_cartesian_coordinates) :: z_u(nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)!,  &
  TYPE(t_subset_range), POINTER :: cells_in_domain! , all_cells
  TYPE(t_patch), POINTER         :: p_patch 
  !-----------------------------------------------------------------------
  p_patch         => p_patch_3D%p_patch_2D(1)
  cells_in_domain => p_patch%cells%in_domain
  !-----------------------------------------------------------------------

  z_u(1:nproma,1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)%x(1) = 0.0_wp
  z_u(1:nproma,1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)%x(2) = 0.0_wp
  z_u(1:nproma,1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)%x(3) = 0.0_wp

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
        z_u(jc,jk,jb)%x =                                                &
        & (p_aux%bc_top_veloc_cc(jc,jb)%x - p_param%A_veloc_v(jc,jk+1,jb)&
        &*(p_diag%p_vn(jc,jk,jb)%x        - p_diag%p_vn(jc,jk+1,jb)%x)   &
        & /p_patch_3D%p_patch_1D(1)%del_zlev_i(jk+1))/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)

        !1b) ocean bottom 
        jk = z_dolic
        z_u(jc,jk,jb)%x =                                      &
        &  (p_param%A_veloc_v(jc,jk,jb)                        &
        & *(p_diag%p_vn(jc,jk-1,jb)%x-p_diag%p_vn(jc,jk,jb)%x) &
        & / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)              &
        & - p_aux%bc_bot_veloc_cc(jc,jb)%x)/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)

        !1c) ocean interior 
        DO jk = slev+1, z_dolic-1
           z_u(jc,jk,jb)%x = &
           &   ( p_param%A_veloc_v(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x - &
           &                                  p_diag%p_vn(jc,jk,jb)%x)    &
           & / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)              &
           &  - p_param%A_veloc_v(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x -  &
           &                                   p_diag%p_vn(jc,jk+1,jb)%x) &
           & / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk+1))/p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)            
        END DO ! jk ocean interior
        ENDIF  ! dolic>0
      END DO
    END DO
  !END DO

  CALL sync_patch_array(SYNC_C, p_patch, z_u(:,:,:)%x(1))
  CALL sync_patch_array(SYNC_C, p_patch, z_u(:,:,:)%x(2))
  CALL sync_patch_array(SYNC_C, p_patch, z_u(:,:,:)%x(3))  
  ! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)  
  CALL map_cell2edges_3D( p_patch_3D, z_u,laplacian_vn_out,p_op_coeff)
  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=4  ! output print level (1-5, fix)
  CALL dbg_print('VelDiffMim: Laplacian'     ,laplacian_vn_out         ,str_module,idt_src)
  !---------------------------------------------------------------------

END SUBROUTINE velocity_diffusion_vert_explicit
!-------------------------------------------------------------------------  
! mpi parallelized, sync required
SUBROUTINE tracer_diffusion_horz(p_patch_3D, trac_in, p_os, K_T, diff_flx, subset_range)

  !
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
  !
  ! Patch on which computation is performed
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  TYPE(t_hydro_ocean_state), TARGET :: p_os
  REAL(wp), INTENT(in)              :: K_T(:,:,:) !mixing coefficient for tracer
  REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
  !
  !Local variables
  INTEGER                       :: slev, elev     ! vertical start and end level
  INTEGER                       :: jk, jb, je
  INTEGER                       :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER                       :: i_startidx_e, i_endidx_e
  REAL(wp)                      :: delta_z
  TYPE(t_subset_range), POINTER :: edges_in_domain
  TYPE(t_patch), POINTER        :: p_patch 
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
  !-------------------------------------------------------------------------------
   p_patch        => p_patch_3D%p_patch_2D(1)
  edges_in_domain => p_patch%edges%in_domain
  !-------------------------------------------------------------------------------
  slev = 1
  elev = n_zlev

!  IF ( iswm_oce /= 1) THEN
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      jk=1
      DO je = i_startidx_e, i_endidx_e

        !IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
        IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
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

          !IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
          IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
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
!   ELSEIF ( iswm_oce == 1) THEN
!     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
!       DO jk = slev, elev
!         DO je = i_startidx_e, i_endidx_e
!           !IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
!           IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
! 
!           !Get indices of two adjacent triangles
!           il_c1 = p_patch%edges%cell_idx(je,jb,1)
!           ib_c1 = p_patch%edges%cell_blk(je,jb,1)
!           il_c2 = p_patch%edges%cell_idx(je,jb,2)
!           ib_c2 = p_patch%edges%cell_blk(je,jb,2)
! 
!           delta_z  = p_os%p_diag%thick_e(je,jb)
! 
!           diff_flx(je,jk,jb) = K_T(je,jk,jb)*delta_z*p_patch_3D%wet_e(je,jk,jb)&!v_base%wet_e(je,jk,jb) &
!             &          *(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
!             &          /p_patch%edges%dual_edge_length(je,jb)
!           ENDIF
!         ENDDO
!       END DO
!     END DO
!   ENDIF

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
SUBROUTINE tracer_diffusion_vert_explicit(p_patch_3D,        &
                                    & trac_c,          &
                                    & top_bc_tracer,   & 
                                    & A_v,             &
                                    & div_diff_flx)

  TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
  REAL(wp), INTENT(in)              :: trac_c       (nproma, n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  REAL(wp), INTENT(in)              :: top_bc_tracer(nproma, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  REAL(wp), INTENT(in)              :: A_v(:,:,:) 
  REAL(wp), INTENT(inout)             :: div_diff_flx(nproma, n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  !
  !Local variables
  INTEGER                       :: slev, elev
  INTEGER                       :: jc, jk, jb
  INTEGER                       :: start_index, end_index
  INTEGER                       :: z_dolic
  ! vertical diffusive tracer flux
  REAL(wp)                      :: z_diff_flx(nproma, n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
  TYPE(t_subset_range), POINTER :: all_cells
  TYPE(t_patch), POINTER        :: p_patch 
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:tracer_diffusion_vert')
  !-----------------------------------------------------------------------
  p_patch   => p_patch_3D%p_patch_2D(1)
  all_cells => p_patch%cells%all
  !-----------------------------------------------------------------------
  slev              = 1
  elev              = n_zlev
  z_diff_flx(1:nproma, 1:n_zlev+1,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) = 0.0_wp

  !1 Vertical derivative of tracer
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, start_index, end_index)
    DO jc = start_index, end_index
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

END SUBROUTINE tracer_diffusion_vert_explicit

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !! Preconditioning by Leonidas Linardakis, MPI-M (2014)
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit( patch_3D,               &
                                           & ocean_tracer, A_v,   &
                                           & operators_coefficients ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff),TARGET        :: operators_coefficients
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: dt_inv, diagonal_product
    REAL(wp), POINTER   :: field_column(:,:,:)
    INTEGER  :: bottom_level
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (bottom_level < 1 ) CYCLE

        DO jk=1,bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, bottom_level-1
          a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(bottom_level) = -A_v(jc,bottom_level,jb) * inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
        b(bottom_level) = dt_inv - a(bottom_level)

        ! precondition: set diagonal equal to diagonal_product
        diagonal_product = PRODUCT(b(1:bottom_level))

        DO jk = 1, bottom_level
           fact(jk) = diagonal_product / b(jk)
           a(jk)  = a(jk)  * fact(jk)
           b(jk)  = diagonal_product
           c(jk)  = dt_inv * fact(jk) - a(jk) - b(jk)
           column_tracer(jk) = field_column(jc,jk,jb) * dt_inv * fact(jk)
        ENDDO
        c(bottom_level) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk=bottom_level-1, 1, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_tracer( jk ) = column_tracer( jk ) - fact(jk+1) * column_tracer( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
        DO jk =  2, bottom_level
          column_tracer( jk ) = ( column_tracer( jk ) - a( jk ) * column_tracer( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
        ENDDO

      END DO ! jc = start_index, end_index
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE tracer_diffusion_vertical_implicit
  !------------------------------------------------------------------------


!-------------------------------------------------------------------------  
! old version, will be removed
!
!!Subroutine implements implicit vertical diffusion for scalar fields.
!>
!! sbr identical to sbr above but now with homogeneous boundary conditions
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
!! The result diff_column is calculated on in_domain_cells
!  SUBROUTINE tracer_diffusion_vertical_implicit( p_patch_3D,               &
!                                           & ocean_tracer, h_c, A_v,   &
!                                           & p_op_coeff) !,  &
!                                          ! & diff_column)
!
!    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
!    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
!    REAL(wp), INTENT(inout)              :: h_c(:,:)  !surface height, relevant for thickness of first cell, in
!    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
!    TYPE(t_operator_coeff),TARGET     :: p_op_coeff
!    ! REAL(wp), INTENT(inout)           :: diff_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
!    !
!    !Local variables
!    INTEGER :: slev
!    INTEGER :: jc, jk, jb
!    INTEGER :: start_index, end_index
!    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
!    REAL(wp) :: fact
!    ! REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
!    ! REAL(wp) :: inv_prism_thickness(1:n_zlev)
!    REAL(wp) :: dt_inv
!    REAL(wp), POINTER   :: field_column(:,:,:)
!    REAL(wp) :: residual(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
!    REAL(wp) :: column_tracer(1:n_zlev), old_tracer_column(1:n_zlev)
!    REAL(wp) :: prism_thickness(1:n_zlev), inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
!    INTEGER  :: z_dolic
!    TYPE(t_subset_range), POINTER :: cells_in_domain
!    TYPE(t_patch), POINTER         :: p_patch
!    ! CHARACTER(len=max_char_length), PARAMETER :: &
!    !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!    !-----------------------------------------------------------------------
!    p_patch         => p_patch_3D%p_patch_2D(1)
!    cells_in_domain => p_patch%cells%in_domain
!    IF (use_tracer_x_height) THEN
!      field_column => ocean_tracer%concentration_x_height
!    ELSE
!      field_column => ocean_tracer%concentration
!    ENDIF
!    residual          = 0.0_wp
!    old_tracer_column = 0.0_wp
!    !-----------------------------------------------------------------------
!    slev   = 1
!    dt_inv = 1.0_wp/dtime
!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
!      DO jc = start_index, end_index
!        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
!
!        IF (z_dolic > 0 ) THEN
!
!          ! recalculate coefficients
!          prism_thickness(1)     = p_patch_3D%p_patch_1D(1)%del_zlev_m(1) + h_c(jc,jb)
!          inv_prism_thickness(1) = 1.0_wp / prism_thickness(1)
!          DO jk=2,z_dolic
!            prism_thickness(jk)            = p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)
!            ! inv_prism_thickness(jk)        = p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
!            inv_prism_thickness(jk)        = 1.0_wp / p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)
!            inv_prisms_center_distance(jk) = 2.0_wp / (prism_thickness(jk-1) + prism_thickness(jk))
!            !inv_prisms_center_distance(jk) = 1.0_wp / p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)
!          ENDDO
!
!          IF (use_tracer_x_height) THEN
!            ! top level
!            a(1) = 0.0_wp
!            c(1) = -A_v(jc,2,jb) * inv_prism_thickness(2) * inv_prisms_center_distance(2)         * dtime
!            b(1) = 1.0_wp + A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2) * dtime
!            !Fill triangular matrix
!            !b is diagonal a is the lower diagonal, c is the upper
!            DO jk = 2, z_dolic-1
!              a(jk) = -A_v(jc,jk,jb)   * inv_prism_thickness(jk-1) * inv_prisms_center_distance(jk)   * dtime
!              c(jk) = -A_v(jc,jk+1,jb) * inv_prism_thickness(jk+1) * inv_prisms_center_distance(jk+1) * dtime
!              b(jk) = 1.0_wp +                                         &
!                & (A_v(jc,jk,jb)  * inv_prisms_center_distance(jk) +   &
!                &  A_v(jc,jk+1,jb) * inv_prisms_center_distance(jk+1)) * inv_prism_thickness(jk) * dtime
!            END DO
!            !bottom
!            a(z_dolic) = -A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic-1) * inv_prisms_center_distance(z_dolic) * dtime
!            c(z_dolic) = 0.0_wp
!            b(z_dolic) = 1.0_wp + A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic) * inv_prisms_center_distance(z_dolic) * dtime
!
!            ! get locally the column tracer_x_height
!            column_tracer(1:z_dolic) = field_column(jc,1:z_dolic,jb)
!
!          !------------------------------------
!          ELSE ! not use_tracer_x_height
!            ! top level
!            a(1) = 0.0_wp
!            c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
!            b(1) = dt_inv - c(1)
!            !Fill triangular matrix
!            !b is diagonal a is the lower diagonal, c is the upper
!            DO jk = 2, z_dolic-1
!              a(jk) = -A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
!              c(jk) = -A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
!              b(jk) = dt_inv - a(jk) - c(jk)
!            END DO
!            !bottom
!            a(z_dolic) = -A_v(jc,z_dolic,jb) * inv_prism_thickness(z_dolic) * inv_prisms_center_distance(z_dolic)
!            c(z_dolic) = 0.0_wp
!            b(z_dolic) = dt_inv - a(z_dolic)
!
!            ! get locally the column tracer / dt
!            column_tracer(1:z_dolic) = field_column(jc,1:z_dolic,jb) * dt_inv
!
!          ENDIF ! use_tracer_x_height
!          !------------------------------------
!
!          ! solver from lapack
!          !
!          ! eliminate lower diagonal
!          DO jk=slev, z_dolic-1
!            fact = a( jk+1 ) / b( jk )
!            b( jk+1 ) = b( jk+1 ) - fact * c( jk )
!            column_tracer( jk+1 ) = column_tracer( jk+1 ) - fact * column_tracer( jk )
!          ENDDO
!  !        DO jk=slev+1, z_dolic
!  !          a(jk) = 0.0_wp
!  !        ENDDO
!
!          !     Back solve with the matrix U from the factorization.
!          column_tracer( z_dolic ) = column_tracer( z_dolic ) / b( z_dolic )
!          DO jk =  z_dolic-1, 1, -1
!            column_tracer( jk ) = ( column_tracer( jk ) - c( jk ) * column_tracer( jk+1 ) ) / b( jk )
!          ENDDO
!
!          IF (use_tracer_x_height) THEN
!            DO jk = 1, z_dolic
!              ocean_tracer%concentration_x_height(jc,jk,jb) = column_tracer(jk)
!              ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk) * inv_prism_thickness(jk)
!            ENDDO
!          ELSE
!            DO jk = 1, z_dolic
!              ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
!            ENDDO
!          ENDIF
!
!          ! check residual
!          CALL apply_triangular_matrix(z_dolic,c,b,a,column_tracer,old_tracer_column)
!          DO jk=1,z_dolic
!            residual(jc,jk,jb) = field_column(jc,jk,jb)*dt_inv - old_tracer_column(jk)
!          ENDDO
!
!        ENDIF ! z_dolic > 0
!
!      END DO ! jc = start_index, end_index
!    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block
!
!    CALL dbg_print('VertTracDiffImpl:residual' ,residual ,str_module,5, in_subset=p_patch%cells%owned)
!
!    ! CALL sync_patch_array(SYNC_C, p_patch, diff_column)
!
!  END SUBROUTINE tracer_diffusion_vertical_implicit


!-------------------------------------------------------------------------  
!
!!Subroutine implements implicit vertical diffusion for horizontal velocity fields
!!by inverting a scalar field..
!>
!! sbr identical to previous one, except for homogeneous boundary conditions
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!! Optimized ofr round-off erros, Leonidas Linardakis, MPI-M (2014)
  !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit( patch_3D,               &
                                           & velocity, A_v,   &
                                           & operators_coefficients ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    REAL(wp), INTENT(inout)              :: velocity(:,:,:)   ! on edges
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff),TARGET        :: operators_coefficients
    !
    REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev), diagonal_product
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)

    INTEGER :: bottom_level
    INTEGER :: edge_index, jk, edge_block
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER         :: patch_2D

    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_edges       => patch_2D%edges%all
    dt_inv = 1.0_wp/dtime
    !-----------------------------------------------------------------------

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_e(edge_index,edge_block)

        IF (bottom_level < 1 ) CYCLE

        ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
        DO jk=1, bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_e(edge_index,jk,edge_block)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1D(1)%inv_prism_center_dist_e(edge_index,jk,edge_block)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(edge_index,2,edge_block) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, bottom_level-1
          a(jk) = - A_v(edge_index,jk,edge_block)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(edge_index,jk+1,edge_block) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(bottom_level) = -A_v(edge_index,bottom_level,edge_block) * inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
        b(bottom_level) = dt_inv - a(bottom_level)

        ! precondition: set diagonal equal to diagonal_product
        diagonal_product = PRODUCT(b(1:bottom_level))

        DO jk = 1, bottom_level
           fact(jk) = diagonal_product / b(jk)
           a(jk)  = a(jk)  * fact(jk)
           b(jk)  = diagonal_product
           c(jk)  = dt_inv * fact(jk) - a(jk) - b(jk)
           column_velocity(jk) = velocity(edge_index,jk,edge_block) * dt_inv * fact(jk)
        ENDDO
        c(bottom_level) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk = bottom_level-1, 1, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_velocity( jk ) = column_velocity( jk ) - fact(jk+1) * column_velocity( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_velocity( 1 ) = column_velocity( 1 ) / b( 1 )
        DO jk = 2, bottom_level
          column_velocity( jk ) = ( column_velocity( jk ) - a( jk ) * column_velocity( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          velocity(edge_index,jk,edge_block) = column_velocity(jk)
        ENDDO

      END DO ! edge_index = start_index, end_index
    END DO ! edge_block = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE velocity_diffusion_vertical_implicit
  !------------------------------------------------------------------------


! OLD routine
!SUBROUTINE velocity_diffusion_vertical_implicit( p_patch_3D,    &
!                                        & field_column,  &
!                                        & A_v,           &
!                                        & p_op_coeff) !,    &
!                                        & diff_column)
!  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
!  REAL(wp), INTENT(inout)           :: field_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
!  !surface height at edges, relevant for thickness of first cell
!  REAL(wp), INTENT(inout)           :: A_v(:,:,:)
!  TYPE(t_operator_coeff), TARGET    :: p_op_coeff
!  REAL(wp), INTENT(inout)             :: diff_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
!  !
!  !Local variables
!  INTEGER :: slev
!  INTEGER :: je, jk, jb
!  INTEGER :: i_startidx, i_endidx
!   REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
!  !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
!  REAL(wp),POINTER :: a(:), b(:), c(:)
!  REAL(wp) :: z_tmp
!  REAL(wp) :: dt_inv
!  INTEGER  :: z_dolic
!  REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
!  REAL(wp) :: inv_prism_thickness(1:n_zlev)
!  TYPE(t_subset_range), POINTER :: all_edges
!  TYPE(t_patch), POINTER        :: p_patch
!  ! CHARACTER(len=max_char_length), PARAMETER :: &
!  !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!  !-----------------------------------------------------------------------
!  p_patch   => p_patch_3D%p_patch_2D(1)
!  all_edges => p_patch%edges%all
!  !-----------------------------------------------------------------------
!
!  slev = 1
!  dt_inv=1.0_wp/dtime
!
!  ! diff_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)= 0.0_wp
!
!  field_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)&
!  &=field_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)*dt_inv
!
!
!  !---------DEBUG DIAGNOSTICS-------------------------------------------
!  idt_src=5  ! output print level (1-5, fix)
!  CALL dbg_print('VelDifImplHomIn:field_col' ,field_column, str_module,idt_src, in_subset=p_patch%edges%owned)
!  !---------------------------------------------------------------------
!
!  DO jb = all_edges%start_block, all_edges%end_block
!    CALL get_index_range(all_edges, jb, i_startidx, i_endidx)
!    DO je = i_startidx, i_endidx
!      z_dolic = p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)!!v_base%dolic_e(je,jb)
!
!      IF (p_patch_3D%lsm_e(je,1,jb) <= sea_boundary) THEN
!        IF ( z_dolic >= MIN_DOLIC ) THEN
!
!
!           a => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,1)
!           b => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,2)
!           c => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,3)
!
!
!          !Apply the matrix
!          DO jk=slev, z_dolic-1
!            IF(b(jk)/=0.0_wp)THEN
!              a(jk) = a(jk)/b(jk)
!              c(jk) = c(jk)/b(jk)
!              field_column(je,jk,jb)=field_column(je,jk,jb)/b(jk)
!              b(jk)=1.0_wp
!            ENDIF
!          END DO
!
!          DO jk=slev+1, z_dolic-1
!            b(jk)                  = b(jk)-a(jk)*c(jk-1)
!            field_column(je,jk,jb) = field_column(je,jk,jb)&
!                          &-a(jk)*field_column(je,jk-1,jb)
!            c(jk)                  = c(jk)/b(jk)
!            field_column(je,jk,jb) = field_column(je,jk,jb)/b(jk)
!            b(jk)                  = 1.0_wp
!          END DO
!
!          z_tmp = b(z_dolic)-a(z_dolic)*c(z_dolic-1)
!          z_tmp = (field_column(je,z_dolic,jb)-a(z_dolic)*field_column(je,z_dolic-1,jb))/z_tmp
!
!          field_column(je,z_dolic,jb) = z_tmp
!          DO jk = z_dolic-1,1,-1
!            field_column(je,jk,jb) = field_column(je,jk,jb)-c(jk)*field_column(je,jk+1,jb)
!          END DO
!          DO jk = 1,z_dolic!-1
!            diff_column(je,jk,jb) = field_column(je,jk,jb)
!          END DO
!        !ELSEIF ( z_dolic < MIN_DOLIC ) THEN
!        !  diff_column(je,1:z_dolic,jb) = 0.0_wp
!        !  field_column(je,1:z_dolic,jb)= 0.0_wp
!        ENDIF
!      !ELSEIF( v_base%lsm_e(je,1,jb) > sea_boundary ) THEN
!      !  diff_column(je,1:z_dolic,jb) = field_column(je,1:z_dolic,jb)
!       !diff_column(je,:,jb) = 0.0_wp
!       !field_column(je,:,jb)= 0.0_wp
!      ENDIF
!    END DO
!  END DO
!
!  !---------DEBUG DIAGNOSTICS-------------------------------------------
!  idt_src=5  ! output print level (1-5, fix)
!  CALL dbg_print('VelDifImplHom: field_col'  ,field_column             ,str_module,idt_src, in_subset=p_patch%edges%owned)
!  CALL dbg_print('VelDifImplHom: diff_col'   ,diff_column              ,str_module,idt_src, in_subset=p_patch%edges%owned)
!  !---------------------------------------------------------------------
!
!END SUBROUTINE velocity_diffusion_vertical_implicit

!------------------------------------------------------------------------
END MODULE mo_oce_diffusion
