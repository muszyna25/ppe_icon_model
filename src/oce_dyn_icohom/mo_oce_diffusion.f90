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
USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec!, gc2cc
USE mo_impl_constants,      ONLY: boundary, sea,sea_boundary ,max_char_length, &
  &                               min_rlcell, min_rledge, min_rlcell, MIN_DOLIC
USE mo_parallel_config,     ONLY: nproma
USE mo_ocean_nml,           ONLY: n_zlev, iswm_oce
USE mo_run_config,          ONLY: dtime
USE mo_oce_index,           ONLY: print_mxmn, jkc, jkdim, ipl_src
USE mo_oce_state,           ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
  &                               t_hydro_ocean_aux, v_base
USE mo_model_domain,        ONLY: t_patch
!USE mo_exception,           ONLY: message, finish!, message_text
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_oce_physics,         ONLY: t_ho_params
USE mo_scalar_product,      ONLY: map_cell2edges, primal_map_c2e
USE mo_oce_math_operators,  ONLY: nabla2_vec_ocean, div_oce
USE mo_intp_data_strc,      ONLY: p_int_state
USE mo_util_subset,         ONLY: t_subset_range, get_index_range
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_idx, global_max

IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!
! PUBLIC INTERFACE
!
INTEGER, PARAMETER  :: top=1

PUBLIC :: velocity_diffusion_horz_mimetic
PUBLIC :: velocity_diffusion_vert_mimetic
PUBLIC :: velocity_diffusion_horz_rbf
PUBLIC :: velocity_diffusion_vert_rbf
PUBLIC :: tracer_diffusion_horz
PUBLIC :: tracer_diffusion_vert_expl
!PUBLIC :: tracer_diffusion_vert_impl
PUBLIC :: tracer_diffusion_vert_impl_hom
!PUBLIC :: veloc_diffusion_vert_impl
PUBLIC :: veloc_diffusion_vert_impl_hom
!INTEGER, PARAMETER :: MIN_DOLIC = 2

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
SUBROUTINE velocity_diffusion_horz_mimetic(p_patch, vn_in, p_param, p_diag, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag)          :: p_diag
  REAL(wp), INTENT(INOUT)           :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)

  !Local variables
  INTEGER :: slev, elev
  INTEGER :: jk, jb, je,jc
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e
  INTEGER :: idx_cartesian
  INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,p_patch%nblks_c)

  TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
!REAL(wp) :: z_grad_vn(nproma,n_zlev,p_patch%nblks_e)
!REAL(wp) :: z_div_grad_vn(nproma,n_zlev,p_patch%nblks_c)
!REAL(wp) :: laplacian_vn_out2(nproma,n_zlev,p_patch%nblks_e)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        

  laplacian_vn_out(:,:,:) = 0.0_wp

  slev = 1
  elev = n_zlev

  ! loop over cells in local domain + halo
  all_cells => p_patch%cells%all
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
  all_edges => p_patch%edges%all
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
  edges_in_domain => p_patch%edges%in_domain
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)

    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e

      IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        z_grad_u(je,jk,jb)%x = p_param%K_veloc_h(je,jk,jb)   &
          &                  *(p_diag%p_vn(il_c2,jk,ib_c2)%x &
          &                  - p_diag%p_vn(il_c1,jk,ib_c1)%x)&
          &                  / p_patch%edges%dual_edge_length(je,jb)
      ELSE
        z_grad_u(je,jk,jb)%x = 0.0_wp
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

         IF ( v_base%lsm_oce_c(jc,jk,jb) >= boundary ) THEN
           z_div_grad_u(jc,jk,jb)%x = 0.0_wp
         ELSE
          z_div_grad_u(jc,jk,jb)%x =  &
            z_grad_u(iidx(jc,jb,1),jk,iblk(jc,jb,1))%x * p_int_state(1)%geofac_div(jc,1,jb) + &
            z_grad_u(iidx(jc,jb,2),jk,iblk(jc,jb,2))%x * p_int_state(1)%geofac_div(jc,2,jb) + &
            z_grad_u(iidx(jc,jb,3),jk,iblk(jc,jb,3))%x * p_int_state(1)%geofac_div(jc,3,jb)
        ENDIF

      END DO
    END DO
  END DO
  DO idx_cartesian = 1,3
    CALL sync_patch_array(SYNC_C, p_patch,z_div_grad_u(:,:,:)%x(idx_cartesian) )
  END DO

  !Step 3: Map divergence back to edges
  CALL map_cell2edges( p_patch, z_div_grad_u, laplacian_vn_out)
  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)


END SUBROUTINE velocity_diffusion_horz_mimetic
!-------------------------------------------------------------------------  
!
!  
!>
!! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
!!    implemented as Average(div( K_H grad RBF v)).
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
!! mpi parallelized, sync required
SUBROUTINE velocity_diffusion_horz_rbf(p_patch, vn_in, p_param, p_diag, laplacian_vn_out)
  !
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  REAL(wp), INTENT(inout)            :: vn_in(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
  TYPE(t_hydro_ocean_diag)          :: p_diag
  REAL(wp), INTENT(out)             :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)
  !
  !Local variables
  INTEGER :: slev, elev
  INTEGER :: jk, jb, je!,jc
  INTEGER :: i_startidx_e, i_endidx_e
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  TYPE(t_subset_range), POINTER :: edges_in_domain        

  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_ocean_semi_implicit_ab:velocity_diffusion_horz')
  !-------------------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')        

  ! #slo# set intent out variable to zero due to nag -nan compiler-option
  laplacian_vn_out(:,:,:) = 0.0_wp

  slev = 1
  elev = n_zlev

  edges_in_domain => p_patch%edges%in_domain

  !Step 2: Multiply laplacian with mixing coefficients
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)

    DO jk = slev, elev
      DO je = i_startidx_e, i_endidx_e

        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          laplacian_vn_out(je,jk,jb) = p_param%K_veloc_h(je,jk,jb)*laplacian_vn_out(je,jk,jb)
        ELSE
          laplacian_vn_out(je,jk,jb) = 0.0_wp
        ENDIF
      ENDDO
    END DO
  END DO
  CALL sync_patch_array(SYNC_E, p_patch,laplacian_vn_out)

END SUBROUTINE velocity_diffusion_horz_rbf
!-------------------------------------------------------------------------
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
SUBROUTINE velocity_diffusion_vert_mimetic( p_patch, p_diag, p_aux,h_c,p_param, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  TYPE(t_hydro_ocean_diag)          :: p_diag
  TYPE(t_hydro_ocean_aux)           :: p_aux
  REAL(wp), INTENT(in)              :: h_c(:,:) 
  TYPE(t_ho_params), INTENT(in)     :: p_param
  REAL(wp)                          :: laplacian_vn_out(:,:,:)

  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: jc, jk, jb, z_dolic
  INTEGER :: i_startidx, i_endidx
  TYPE(t_cartesian_coordinates) :: z_u(nproma,n_zlev+1,p_patch%nblks_c)!,  &
  TYPE(t_subset_range), POINTER :: all_cells
  !  &                              z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:veloc diffusion vert mimetic')
  !-----------------------------------------------------------------------
  all_cells => p_patch%cells%all

  z_u(nproma,n_zlev+1,p_patch%nblks_c)%x = 0.0_wp
  slev       = 1
  elev       = n_zlev

  !1 Vertical derivative of cell velocity vector times horizontal velocity
  ! loop runs now from slev to z_dolic:
  !DO jk = slev, elev
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      z_dolic = v_base%dolic_c(jc,jb)
      IF (z_dolic >= MIN_DOLIC) THEN
        !1a) ocean surface 
        jk = slev
   ! #slo 2011-09-08 - include effect of surface boundary condition on G_n
   !                   eliminate in vertical diffusion: laplacian_vert
   !    &( &
        z_u(jc,jk,jb)%x =                                                &
        & (p_aux%bc_top_veloc_cc(jc,jb)%x - p_param%A_veloc_v(jc,jk+1,jb)&
        &*(p_diag%p_vn(jc,jk,jb)%x        - p_diag%p_vn(jc,jk+1,jb)%x)   &
        &/v_base%del_zlev_i(jk+1))/v_base%del_zlev_m(jk)
        !1b) ocean bottom 
        !ELSEIF ( jk == v_base%dolic_c(jc,jb) ) THEN
        !ELSEIF ( jk == z_dolic ) THEN
        jk = z_dolic
  !        IF(v_base%n_zlev>=2)THEN
        z_u(jc,jk,jb)%x =                                      &
        &  (p_param%A_veloc_v(jc,jk,jb)                        &
        & *(p_diag%p_vn(jc,jk-1,jb)%x-p_diag%p_vn(jc,jk,jb)%x) &
        & / v_base%del_zlev_i(jk)                              &
        & - p_aux%bc_bot_veloc_cc(jc,jb)%x)/v_base%del_zlev_m(jk)

  !        ENDIF
        !1c) ocean interior 
        !ELSEIF( jk>slev .AND. jk < v_base%dolic_c(jc,jb) ) THEN
        DO jk = slev+1, z_dolic-1
           z_u(jc,jk,jb)%x = &
           &   ( p_param%A_veloc_v(jc,jk,jb)*(p_diag%p_vn(jc,jk-1,jb)%x - &
           &                                  p_diag%p_vn(jc,jk,jb)%x)    &
           &   /v_base%del_zlev_i(jk)                                     &
           &  - p_param%A_veloc_v(jc,jk+1,jb)*(p_diag%p_vn(jc,jk,jb)%x -  &
           &                                   p_diag%p_vn(jc,jk+1,jb)%x) &
           &   /v_base%del_zlev_i(jk+1))/v_base%del_zlev_m(jk)
        ! ENDIF  ! jk-condition
       !ENDIF    ! at least 2 vertical layers
        END DO ! jk ocean interior
        ENDIF  ! dolic>0
      END DO
    END DO
  !END DO

  ! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)  
  CALL map_cell2edges( p_patch, z_u, laplacian_vn_out)
  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)

   DO jk=slev, elev
     ipl_src=4  ! output print level (1-5, fix)
     CALL print_mxmn('vert diffusion',jk,laplacian_vn_out(:,:,:),n_zlev, &
       &              p_patch%nblks_e,'dif',ipl_src)
   END DO

END SUBROUTINE velocity_diffusion_vert_mimetic
!-------------------------------------------------------------------------  


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
SUBROUTINE velocity_diffusion_vert_rbf( p_patch, u_c, v_c, h_c, top_bc_u_c, top_bc_v_c,&
                          &  bot_bc_u_c,  bot_bc_v_c,p_param, laplacian_vn_out)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  ! Components of cell based variable which is vertically advected
  REAL(wp), INTENT(inout) :: u_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(inout) :: v_c(:,:,:) ! dim: (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(in) :: h_c(:,:) ! dim: (nproma,nblks_c)
  !
  ! Top boundary condition for cell based variables
  REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
  !
  ! Bottom boundary condition for cell based variables
  REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,nblks_c)
  TYPE(t_ho_params), INTENT(in)     :: p_param
  ! variable in which horizontally advected velocity is stored
  REAL(wp)            :: laplacian_vn_out(:,:,:)

  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: jc, jk, jb
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
  TYPE(t_cartesian_coordinates) :: zu_cc(nproma,n_zlev,p_patch%nblks_c)
  REAL(wp) :: z_u(nproma,n_zlev+1,p_patch%nblks_c),  &
    &         z_v(nproma,n_zlev+1,p_patch%nblks_c)
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_oce_diffusion:veloc diffusion vert rbf')
  !-----------------------------------------------------------------------

  ! #slo# set local variable to zero due to nag -nan compiler-option
  z_u(:,:,:) = 0.0_wp
  z_v(:,:,:) = 0.0_wp

  slev       = 1
  elev       = n_zlev
  i_startblk = p_patch%cells%start_blk(1,1)
  i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
  !1 Vertical derivative of cell velocity vector times horizontal velocity
  DO jk = slev, elev

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1,min_rlcell)
      DO jc = i_startidx, i_endidx
      IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        !1a) 0cean surface
        !check if we have at least two layers of water
       !IF (v_base%dolic_c(jc,jb) >= 2) THEN
        IF(jk==slev)THEN
            ! u-component
            z_u(jc,jk,jb) =      &
            &( top_bc_u_c(jc,jb)       &
            & -p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk,jb)-u_c(jc,jk+1,jb)) &
            & /v_base%del_zlev_i(jk+1))/v_base%del_zlev_m(jk)
            ! v-component
            z_v(jc,jk,jb) =      &
            &( top_bc_v_c(jc,jb)       &
            & -p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
            &/v_base%del_zlev_i(jk+1)) &
            & /v_base%del_zlev_m(jk)
            !write(*,*)'RBD:vert-diff:top',jc,jk,jb, z_u_i(jc,jk,jb),top_bc_u_c(jc,jb)

        !1b) ocean bottom 
        ELSEIF ( jk == v_base%dolic_c(jc,jb) ) THEN

          IF(v_base%n_zlev>=2)THEN
            ! u-component
            z_u(jc,jk,jb)              &
            & = (p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk-1,jb)-u_c(jc,jk,jb))&
            & /v_base%del_zlev_i(jk)&
            & - bot_bc_u_c(jc,jb))/v_base%del_zlev_m(jk)
            ! v-component
            z_v(jc,jk,jb)&
            & = (p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk-1,jb)-v_c(jc,jk,jb))&
            & /v_base%del_zlev_i(jk)&
            & -bot_bc_v_c(jc,jb))/v_base%del_zlev_m(jk)
          !write(*,*)'u-diff botom:',jk,u_c(jc,elev-1,jb),u_c(jc,elev,jb),u_c(jc,elev,jb),bot_bc_u_c(jc,jb), z_u_i(jc,jk,jb)
          ENDIF

        !1c) ocean interior 
        ELSEIF(jk>slev .AND. jk < v_base%dolic_c(jc,jb) ) THEN
          ! u-component
          z_u(jc,jk,jb)&
          & = &
          & ( p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk-1,jb)-u_c(jc,jk,jb))&
          &/v_base%del_zlev_i(jk)&
          & - p_param%A_veloc_v(jc,jk+1,jb)  *(u_c(jc,jk,jb)-u_c(jc,jk+1,jb))&
          &/v_base%del_zlev_i(jk+1)) &
          & / v_base%del_zlev_m(jk)
          ! v-component
          z_v(jc,jk,jb)&
          & = &
          & ( p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk-1,jb)-v_c(jc,jk,jb))&
          &/v_base%del_zlev_i(jk)&
          & - p_param%A_veloc_v(jc,jk,jb)  *(v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
          &/v_base%del_zlev_i(jk+1)) &
          & / v_base%del_zlev_m(jk)

        !1c) ocean bottom 
        ENDIF  ! jk-condition
        ENDIF
        CALL gvec2cvec( z_u(jc,jk,jb), z_v(jc,jk,jb),       &
        &               p_patch%cells%center(jc,jb)%lon,                &
        &               p_patch%cells%center(jc,jb)%lat,                &
        &               zu_cc(jc,jk,jb)%x(1),zu_cc(jc,jk,jb)%x(2),zu_cc(jc,jk,jb)%x(3))

      !ENDIF   ! at least 2 vertical layers
      END DO
    END DO
  END DO

  ! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)
  CALL map_cell2edges( p_patch, zu_cc, laplacian_vn_out)
  CALL sync_patch_array(SYNC_E, p_patch, laplacian_vn_out)
   DO jk=slev, elev
     ipl_src=4  ! output print level (1-5, fix)
     CALL print_mxmn('vert diffusion',jk,laplacian_vn_out(:,:,:),n_zlev, &
       &              p_patch%nblks_e,'dif',ipl_src)
   END DO

END SUBROUTINE velocity_diffusion_vert_rbf

!------------------------------------------------------------------------
SUBROUTINE tracer_diffusion_horz(p_patch, trac_in, p_os, K_T, diff_flx)
  !
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
  !
  ! Patch on which computation is performed
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,p_patch%nblks_c)
  TYPE(t_hydro_ocean_state), TARGET :: p_os
  REAL(wp), INTENT(in)              :: K_T(:,:,:) !mixing coefficient for tracer
  REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,p_patch%nblks_e)
  !
  !
  !Local variables
  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: jk, jb, je
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: rl_start, rl_end
  REAL(wp) :: delta_z
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
  !-------------------------------------------------------------------------------
  rl_start = 1
  rl_end   = min_rledge

  slev = 1
  elev = n_zlev

  i_startblk_e = p_patch%edges%start_blk(1,1)
  i_endblk_e   = p_patch%edges%end_blk(min_rledge,1)

  IF ( iswm_oce /= 1) THEN
    DO jb = i_startblk_e, i_endblk_e
      CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e,&
                       &  i_startidx_e, i_endidx_e,&
                       &  rl_start, rl_end)
      DO jk = slev, elev

        delta_z  = v_base%del_zlev_m(jk)

        DO je = i_startidx_e, i_endidx_e

          IF(jk==slev) delta_z = v_base%del_zlev_m(slev) + p_os%p_diag%h_e(je,jb)

          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)

          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            diff_flx(je,jk,jb) = K_T(je,jk,jb)*delta_z&
                       &*(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
                       &/p_patch%edges%dual_edge_length(je,jb)
          ELSE
            diff_flx(je,jk,jb) = 0.0_wp
          ENDIF
        ENDDO
      END DO
    END DO
  ELSEIF ( iswm_oce == 1) THEN
    DO jb = i_startblk_e, i_endblk_e
      CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e,&
                       &  i_startidx_e, i_endidx_e,&
                       &  rl_start, rl_end)
      DO jk = slev, elev
        DO je = i_startidx_e, i_endidx_e


          !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)

          delta_z  = p_os%p_diag%thick_e(je,jb)

           diff_flx(je,jk,jb) = K_T(je,jk,jb)*delta_z*v_base%wet_e(je,jk,jb)&
                       &*(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
                       &/p_patch%edges%dual_edge_length(je,jb)
  !write(123,*)'trac diffusive flux',je,jk,jb,diff_flx(je,jk,jb),K_T(je,jk,jb),&
  !&trac_in(il_c2,jk,ib_c2),trac_in(il_c1,jk,ib_c1)
          !ELSE
          !  diff_flx(je,jk,jb) = 0.0_wp
          !ENDIF
        ENDDO
      END DO
    END DO
  ENDIF
  ! Apply divergence to mixing times gradient to get laplacian
  !CALL div_oce( diff_flx, p_patch, laplacian_trac_out)

  ! DO jk=slev, elev
  !   WRITE(*,*)'MAX/MIN horz tracer diffusion ',jk, &
  !     &        MAXVAL(diff_flx(:,jk,:)), MINVAL(diff_flx(:,jk,:))
  ! END DO


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
SUBROUTINE tracer_diffusion_vert_expl( p_patch,        &
                                    & trac_c,dz,       &
                                    & top_bc_tracer,   & 
                                    & bot_bc_tracer,   &
                                    & A_v,             &
                                    & div_diff_flx)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(inout)           :: trac_c(:,:,:)
REAL(wp), INTENT(in)              :: dz(:,:,:)
REAL(wp), INTENT(in)              :: top_bc_tracer(:,:)
REAL(wp), INTENT(in)              :: bot_bc_tracer(:,:)
REAL(wp), INTENT(inout)           :: A_v(:,:,:) 
REAL(wp), INTENT(out)             :: div_diff_flx(:,:,:)
!
!Local variables
INTEGER :: slev, elev
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: z_dolic
REAL(wp) :: z_diff_flx(nproma, n_zlev+1,p_patch%nblks_c)   ! vertical diffusive tracer flux
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:tracer_diffusion_vert')
!-----------------------------------------------------------------------
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev
z_diff_flx(:,:,:) = 0.0_wp
!A_v=0.0_wp
!1 Vertical derivative of tracer
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic  = v_base%dolic_c(jc,jb)
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
          & * (trac_c(jc,jk-1,jb)-trac_c(jc,jk,jb))/v_base%del_zlev_i(jk)
        !ELSE
        !  diff_flx(jc,jk,jb)= 0.0_wp
        !ENDIF
      END DO


      DO jk = 1, z_dolic
        ! positive vertical divergence in direction of w (upward positive)
         div_diff_flx(jc,jk,jb) = (z_diff_flx(jc,jk,jb) - z_diff_flx(jc,jk+1,jb))&
                              &/v_base%del_zlev_m(jk) 
        END DO


      !1c) ocean bottom zero bottom boundary condition
      !diff_flx(jc,z_dolic+1,jb) = bot_bc_tracer(jc,jk)
    ELSE
      div_diff_flx(jc,:,jb)= 0.0_wp
    ENDIF
  END DO
END DO
 DO jk=slev, elev
   ipl_src=4  ! output print level (1-5, fix)
   CALL print_mxmn('vert diffusion expl',jk,div_diff_flx(:,:,:),n_zlev, &
     &              p_patch%nblks_c,'dif',ipl_src)
   CALL print_mxmn('vrt.dif.expl.diff-flx',jk,z_diff_flx(:,:,:),n_zlev+1, &
     &              p_patch%nblks_c,'dif',ipl_src)
 END DO

END subroutine tracer_diffusion_vert_expl
! ! !-------------------------------------------------------------------------  
! ! !
! ! !!Subroutine implements implicit vertical diffusion for scalar fields.
! ! !>
! ! !!
! ! !!
! ! !! @par Revision History
! ! !! Developed  by  Peter Korn, MPI-M (2011).
! ! !!
! ! SUBROUTINE tracer_diffusion_vert_impl( p_patch,   &
! !                                 & field_column,   &
! !                                 & top_bc,         & 
! !                                 & bot_bc,         &
! !                                 & h_c,            &
! !                                 & A_v,            &
! !                                 & diff_column)
! ! 
! ! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! ! REAL(wp), INTENT(inout)           :: field_column(:,:,:)
! ! REAL(wp), INTENT(in)              :: top_bc(:,:)
! ! REAL(wp), INTENT(in)              :: bot_bc(:,:)  
! ! REAL(wp), INTENT(IN)              :: h_c(:,:)           !surface height, relevant for thickness of first cell 
! ! REAL(wp), INTENT(inout)           :: A_v(:,:,:) 
! ! REAL(wp), INTENT(out)             :: diff_column(:,:,:)
! ! !
! ! !Local variables
! ! INTEGER :: slev
! ! INTEGER :: jc, jk, jb
! ! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! ! REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
! ! REAL(wp) :: z_tmp
! ! REAL(wp) :: inv_zinv_i(1:n_zlev)
! ! REAL(wp) :: inv_zinv_m(1:n_zlev)
! ! !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
! ! REAL(wp) :: z_c1(nproma,n_zlev,p_patch%nblks_c)
! ! INTEGER  :: z_dolic
! ! ! CHARACTER(len=max_char_length), PARAMETER :: &
! ! !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
! ! !-----------------------------------------------------------------------
! ! i_startblk = p_patch%cells%start_blk(1,1)
! ! i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
! ! slev = 1
! ! !A_v=0.0001_wp
! ! !dt_inv = 1.0_wp/dtime
! ! !write(*,*)'impl vert trac diff: max/min top bc trac:',maxval(top_bc), minval(top_bc)
! ! !write(*,*)'impl vert trac diff: max/min bot bc trac:',maxval(bot_bc), minval(bot_bc)
! ! ipl_src=5  ! output print level (1-5, fix)
! ! z_c1(:,1,:)=top_bc(:,:)
! ! CALL print_mxmn('IMPL TRC: top bc',1,z_c1(:,:,:),1,p_patch%nblks_c,'dif',ipl_src)
! ! 
! ! diff_column(:,:,:)= field_column(:,:,:)
! ! z_c1= field_column
! ! 
! ! a(slev:n_zlev)    = 0.0_wp 
! ! b(slev:n_zlev)    = 0.0_wp 
! ! c(slev:n_zlev)    = 0.0_wp
! ! !bet(slev:n_zlev)  = 1.0_wp
! ! !gam(slev:n_zlev)  = 0.0_wp
! ! inv_zinv_i(slev:n_zlev)  = 0.0_wp
! ! inv_zinv_m(slev:n_zlev)  = 0.0_wp
! ! 
! ! DO jb = i_startblk, i_endblk
! !   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
! !                        i_startidx, i_endidx, 1,min_rlcell)
! !   DO jc = i_startidx, i_endidx
! !     z_dolic = v_base%dolic_c(jc,jb)
! ! 
! !     IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
! !       ! IF ( z_dolic >=MIN_DOLIC +1) THEN
! !       IF ( z_dolic >=MIN_DOLIC) THEN
! !         inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
! !         inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)
! ! 
! ! !         IF(h_c(jc,jb)/=0.0_wp)THEN
! ! !           inv_zinv_i(slev)=1.0_wp/(v_base%del_zlev_i(slev)+h_c(jc,jb))
! ! !           inv_zinv_m(slev)=1.0_wp/(v_base%del_zlev_m(slev)+h_c(jc,jb))
! ! ! !         ENDIF
! !         !Fill triangular matrix
! !         !b is diagonal a and c are upper and lower band
! !         DO jk = slev+1, z_dolic-1
! !           a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)*dtime
! !           c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1) *dtime
! !           b(jk) = 1.0_wp-a(jk)-c(jk)
! !         END DO
! ! 
! !         ! The first row
! !         c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1) *dtime
! !         a(slev) = 0.0_wp
! !         b(slev) = 1.0_wp- c(slev) !- a(slev) 
! ! 
! !         ! The last row
! !         a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic) *dtime
! !         c(z_dolic) = 0.0_wp
! !         b(z_dolic) = 1.0_wp - a(z_dolic)! - c(z_dolic)
! ! 
! !         field_column(jc,slev,jb)   = field_column(jc,slev,jb)   ! &
! !                                    !&+ top_bc(jc,jb)*inv_zinv_m(slev)*dtime
! !         !field_column(jc,z_dolic,jb)= field_column(jc,z_dolic,jb) !&
! !         !                           ! &- bot_bc(jc,jb)*inv_zinv_m(z_dolic)*dtime
! ! !----------------------------------------------------------- 
! ! ! !         bet(slev) = 1.0_wp/b(slev)
! ! ! !         DO jk=slev+1, z_dolic
! ! ! !           gam(jk) = a(jk-1)*bet(jk-1)
! ! ! !           IF((b(jk) - c(jk)*gam(jk))/=0.0_wp)THEN 
! ! ! !             bet(jk) = 1.0_wp/(b(jk) - c(jk)*gam(jk))
! ! ! !           ENDIF
! ! ! !         END DO
! ! ! !         diff_column(jc,slev,jb)=field_column(jc,slev,jb)*bet(slev)
! ! ! !         DO jk=slev+1, z_dolic
! ! ! !           diff_column(jc,jk,jb)=bet(jk)*&
! ! ! !           &(field_column(jc,jk,jb)-c(jk)*diff_column(jc,jk-1,jb))
! ! ! !         END DO
! ! ! !         ! Backward sweep
! ! ! !         DO jk=z_dolic-1,slev,-1
! ! ! !           diff_column(jc,jk,jb) = diff_column(jc,jk,jb)  &
! ! ! !           &            -gam(jk+1)*diff_column(jc,jk+1,jb)
! ! ! !         END DO
! ! ! !----------------------------------------------------------- 
! ! !         DO jk=slev, z_dolic-2
! ! !           IF(b(jk)/=0.0_wp)THEN
! ! !             a(jk) = a(jk)/b(jk)
! ! !             c(jk) = c(jk)/b(jk)
! ! !             field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! ! !             b(jk)=1.0_wp
! ! !           ENDIF
! ! !         END DO
! ! !         DO jk=slev+1, z_dolic-2
! ! !           b(jk)=b(jk)-a(jk)*c(jk-1)
! ! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)&
! ! !                         &-a(jk)*field_column(jc,jk-1,jb)
! ! !           c(jk)=c(jk)/b(jk)
! ! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! ! !           b(jk)=1.0_wp
! ! !         END DO
! ! !         z_tmp=b(z_dolic-1)-a(z_dolic-1)*c(z_dolic-2)
! ! !         z_tmp=(field_column(jc,z_dolic-1,jb)-a(z_dolic-1)*field_column(jc,z_dolic-2,jb))/z_tmp
! ! ! 
! ! !         field_column(jc,z_dolic-1,jb)=z_tmp
! ! !         DO jk=z_dolic-2,1,-1
! ! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
! ! !         END DO
! ! !         DO jk=1,z_dolic-1
! ! !           diff_column(jc,jk,jb)=field_column(jc,jk,jb)
! ! !         END DO
! ! !----------------------------------------------------------- 
! !         DO jk=slev, z_dolic-1
! !           IF(b(jk)/=0.0_wp)THEN
! !             a(jk) = a(jk)/b(jk)
! !             c(jk) = c(jk)/b(jk)
! !             field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! !             b(jk)=1.0_wp
! !           ENDIF
! !         END DO
! ! 
! !         DO jk=slev+1, z_dolic-1
! !           b(jk)=b(jk)-a(jk)*c(jk-1)
! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)&
! !                         &-a(jk)*field_column(jc,jk-1,jb)
! !           c(jk)=c(jk)/b(jk)
! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! !           b(jk)=1.0_wp
! !         END DO
! ! 
! !         z_tmp=b(z_dolic)-a(z_dolic)*c(z_dolic-1)
! !         z_tmp=(field_column(jc,z_dolic,jb)-a(z_dolic)*field_column(jc,z_dolic-1,jb))/z_tmp
! ! 
! !         field_column(jc,z_dolic,jb)=z_tmp
! !         DO jk=z_dolic-1,1,-1
! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
! !         END DO
! !         DO jk=1,z_dolic-1
! !           diff_column(jc,jk,jb)=field_column(jc,jk,jb)
! !         END DO
! ! !------------------------------------------
! ! !       IF(field_column(jc,1,jb)/=0.0_wp)THEN
! ! ! !      IF(field_column(jc,2,jb)/=diff_column(jc,2,jb))THEN 
! ! !      !write(234,*)'top bc',top_bc(jc,jb)/zinv
! ! !      !write(234,*)'coffs i',inv_zinv_i(:)
! ! !      !write(234,*)'coffs m',inv_zinv_m(:)
! ! !      write(234,*)'coffs A_V',jc,jb,A_v(jc,:,jb)
! ! !      write(234,*)'mat up  a:', a
! ! !      write(234,*)'mat dia b:', b
! ! !      write(234,*)'mat dow c:', c
! ! !       !   write(234,*)'mat gam:', gam
! ! !       write(234,*)'in :', z_c1(jc,:,jb)!field_column(jc,:,jb) 
! ! !       write(234,*)'out:',diff_column(jc,:,jb)     
! ! !      !write(234,*)'sum:',sum(field_column(jc,:,jb))/z_dolic,&
! ! !      ! & sum(diff_column(jc,:,jb))/z_dolic, &
! ! !      ! &(sum(diff_column(jc,:,jb))/z_dolic)/(sum(field_column(jc,:,jb))/z_dolic)
! ! !      write(234,*)
! ! !  !   ENDIF
! ! !    ENDIF
! !       ELSEIF ( z_dolic <MIN_DOLIC ) THEN
! !         diff_column(jc,:,jb) = 0.0_wp!field_column(jc,:,jb)
! !       ENDIF
! !     ELSEIF( v_base%lsm_oce_c(jc,1,jb) > sea_boundary ) THEN
! !       diff_column(jc,:,jb) = field_column(jc,:,jb)
! !     ENDIF
! ! 
! ! 
! !   END DO
! ! END DO
! ! 
! ! DO jk=slev, n_zlev
! !   ipl_src=5  ! output print level (1-5, fix)
! !   CALL print_mxmn('IMPL TRC: bef.vtrc.dif',jk,z_c1(:,:,:),n_zlev, &
! !     &              p_patch%nblks_c,'dif',ipl_src)
! ! ! WRITE(*,*)'IMPLICIT TRACER: MAX/MIN before vert tracer diffusion ',jk, &
! ! !   &        MAXVAL(field_column(:,jk,:)), MINVAL(field_column(:,jk,:))
! ! END DO
! ! DO jk=slev, n_zlev
! !   ipl_src=5  ! output print level (1-5, fix)
! !   CALL print_mxmn('IMPL TRC: aft.vtrc.dif',jk,diff_column(:,:,:),n_zlev, &
! !     &              p_patch%nblks_c,'dif',ipl_src)
! ! END DO
! ! 
! ! END subroutine tracer_diffusion_vert_impl
!-------------------------------------------------------------------------  
!
!!Subroutine implements implicit vertical diffusion for scalar fields.
!>
!! sbr identical to sbr above but now with homogeneous boundary conditions
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
SUBROUTINE tracer_diffusion_vert_impl_hom( p_patch,   &
                                & field_column,   &
                                & h_c,            &
                                & A_v,            &
                                & diff_column)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(inout)           :: field_column(:,:,:)
REAL(wp), INTENT(IN)              :: h_c(:,:)           !surface height, relevant for thickness of first cell 
REAL(wp), INTENT(inout)           :: A_v(:,:,:) 
REAL(wp), INTENT(out)             :: diff_column(:,:,:)
!
!Local variables
INTEGER :: slev
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
REAL(wp) :: z_tmp
REAL(wp) :: inv_zinv_i(1:n_zlev)
REAL(wp) :: inv_zinv_m(1:n_zlev)
!REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
!REAL(wp) :: z_c1(nproma,1,p_patch%nblks_c)
INTEGER  :: z_dolic
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!-----------------------------------------------------------------------
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
!A_v=0.0001_wp
!dt_inv = 1.0_wp/dtime
!write(*,*)'impl vert trac diff: max/min top bc trac:',maxval(top_bc), minval(top_bc)
!write(*,*)'impl vert trac diff: max/min bot bc trac:',maxval(bot_bc), minval(bot_bc)
ipl_src=5  ! output print level (1-5, fix)


diff_column(:,:,:)= field_column(:,:,:)
a(slev:n_zlev)    = 0.0_wp 
b(slev:n_zlev)    = 0.0_wp 
c(slev:n_zlev)    = 0.0_wp
!bet(slev:n_zlev)  = 1.0_wp
!gam(slev:n_zlev)  = 0.0_wp
inv_zinv_i(slev:n_zlev)  = 0.0_wp
inv_zinv_m(slev:n_zlev)  = 0.0_wp

DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    z_dolic             = v_base%dolic_c(jc,jb)

    IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
      !IF ( z_dolic >=MIN_DOLIC+1 ) THEN
      IF ( z_dolic >=MIN_DOLIC ) THEN

        inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
        inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)

        !Fill triangular matrix
        !b is diagonal a and c are upper and lower band
        DO jk = slev+1, z_dolic-1
          a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)*dtime
          c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)*dtime
          b(jk) = 1.0_wp-a(jk)-c(jk)
        END DO

        ! The first row
         c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)*dtime
         a(slev) = 0.0_wp           
         b(slev) = 1.0_wp- c(slev) !- a(slev) 

        ! The last row
        a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)*dtime
        c(z_dolic) = 0.0_wp
        b(z_dolic) = 1.0_wp - a(z_dolic)! - c(z_dolic)

        ! The matrix is now complete.
!------------------------------------------
! !         bet(slev) = 1.0_wp/b(slev)
! !         DO jk=slev+1, z_dolic
! !           gam(jk) = a(jk-1)*bet(jk-1)
! !           IF((b(jk) - c(jk)*gam(jk))/=0.0_wp)THEN 
! !             bet(jk) = 1.0_wp/(b(jk) - c(jk)*gam(jk))
! !           ENDIF
! !         END DO
! ! 
! !         diff_column(jc,slev,jb)=field_column(jc,slev,jb)!*bet(slev)
! !         diff_column(jc,z_dolic,jb) = field_column(jc,z_dolic,jb)
! !         DO jk=slev+1, z_dolic
! !           diff_column(jc,jk,jb)=bet(jk)*&
! !           &(field_column(jc,jk,jb)-c(jk)*diff_column(jc,jk-1,jb))
! !         END DO
! ! 
! !         ! Backward sweep
! !         DO jk=z_dolic-1,slev,-1
! !           diff_column(jc,jk,jb) = diff_column(jc,jk,jb)  &
! !           &            -gam(jk+1)*diff_column(jc,jk+1,jb)
! !         END DO
! ! 
!------------------------------------------
!         DO jk=slev, z_dolic-2
!           IF(b(jk)/=0.0_wp)THEN
!             a(jk) = a(jk)/b(jk)
!             c(jk) = c(jk)/b(jk)
!             field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
!             b(jk)=1.0_wp
!           ENDIF
!         END DO
!         DO jk=slev+1, z_dolic-2
!           b(jk)=b(jk)-a(jk)*c(jk-1)
!           field_column(jc,jk,jb)=field_column(jc,jk,jb)&
!                         &-a(jk)*field_column(jc,jk-1,jb)
!           c(jk)=c(jk)/b(jk)
!           field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
!           b(jk)=1.0_wp
!         END DO
!         z_tmp=b(z_dolic-1)-a(z_dolic-1)*c(z_dolic-2)
!         z_tmp=(field_column(jc,z_dolic-1,jb)-a(z_dolic-1)*field_column(jc,z_dolic-2,jb))/z_tmp
! 
!         field_column(jc,z_dolic-1,jb)=z_tmp
!         DO jk=z_dolic-2,1,-1
!           field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
!         END DO
!         DO jk=1,z_dolic-1
!           diff_column(jc,jk,jb)=field_column(jc,jk,jb)
!         END DO
!----------------------------------------------------------- 
        DO jk=slev, z_dolic-1
          IF(b(jk)/=0.0_wp)THEN
            a(jk) = a(jk)/b(jk)
            c(jk) = c(jk)/b(jk)
            field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
            b(jk)=1.0_wp
          ENDIF
        END DO

        DO jk=slev+1, z_dolic-1
          b(jk)=b(jk)-a(jk)*c(jk-1)
          field_column(jc,jk,jb)=field_column(jc,jk,jb)&
                        &-a(jk)*field_column(jc,jk-1,jb)
          c(jk)=c(jk)/b(jk)
          field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
          b(jk)=1.0_wp
        END DO

        z_tmp=b(z_dolic)-a(z_dolic)*c(z_dolic-1)
        z_tmp=(field_column(jc,z_dolic,jb)-a(z_dolic)*field_column(jc,z_dolic-1,jb))/z_tmp

        field_column(jc,z_dolic,jb)=z_tmp
        DO jk=z_dolic-1,1,-1
          field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
        END DO
        DO jk=1,z_dolic-1
          diff_column(jc,jk,jb)=field_column(jc,jk,jb)
        END DO
!      IF(field_column(jc,1,jb)/=0.0_wp)THEN
!     !write(234,*)'top bc',top_bc(jc,jb)/zinv
!     !write(234,*)'coffs i',inv_zinv_i(:)
!     !write(234,*)'coffs m',inv_zinv_m(:)
!     write(234,*)'coffs A_V',A_v(jc,:,jb)
!     write(234,*)'mat up  a:', a
!     write(234,*)'mat dia b:', b
!     write(234,*)'mat dow c:', c
!      write(234,*)'in :', field_column(jc,:,jb) 
!      write(234,*)'out:',diff_column(jc,:,jb)
!      !write(234,*)'sum:',sum(field_column(jc,:,jb))/z_dolic,&
!      !& sum(diff_column(jc,:,jb))/z_dolic, &
!      !&(sum(diff_column(jc,:,jb))/z_dolic)/(sum(field_column(jc,:,jb))/z_dolic)
!     write(234,*)
!    ENDIF
      ELSEIF ( z_dolic <MIN_DOLIC ) THEN
        diff_column(jc,:,jb) = 0.0_wp!field_column(jc,:,jb)
      ENDIF
    ELSEIF( v_base%lsm_oce_c(jc,1,jb) > sea_boundary ) THEN
      diff_column(jc,:,jb) = field_column(jc,:,jb)
    ENDIF


  END DO
END DO
! write(234,*)'-------'
DO jk=slev, n_zlev
  ipl_src=5  ! output print level (1-5, fix)
  CALL print_mxmn('IMPL TRC: aft.vtrc.dif',jk,diff_column(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'dif',ipl_src)
END DO

END subroutine tracer_diffusion_vert_impl_hom
! ! !-------------------------------------------------------------------------  
! ! !
! ! !!Subroutine implements implicit vertical diffusion for hrozontal velocity fields
! ! !!by inverting a scalar field..
! ! !>
! ! !!
! ! !!
! ! !! @par Revision History
! ! !! Developed  by  Peter Korn, MPI-M (2011).
! ! !!
! ! SUBROUTINE veloc_diffusion_vert_impl( p_patch,       &
! !                                     & field_column,  &
! !                                     & top_bc_vn,     & 
! !                                     & bot_bc_vn,     &
! !                                     & h_e,           &
! !                                     & A_v,           &
! !                                     & diff_column)
! ! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! ! REAL(wp), INTENT(inout)           :: field_column(:,:,:)
! ! REAL(wp), INTENT(in)              :: top_bc_vn(:,:)
! ! REAL(wp), INTENT(in)              :: bot_bc_vn(:,:)  
! ! REAL(wp), INTENT(IN)              :: h_e(:,:)           !surface height at edges, relevant for thickness of first cell 
! ! REAL(wp), INTENT(inout)           :: A_v(:,:,:) 
! ! REAL(wp), INTENT(out)             :: diff_column(:,:,:)
! ! !
! ! !Local variables
! ! INTEGER :: slev
! ! INTEGER :: jc, jk, jb
! ! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! ! REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
! ! !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
! ! !REAL(wp) :: z_rhs(1:n_zlev)
! ! REAL(wp) :: z_tmp!, zinv
! ! INTEGER  :: z_dolic
! ! REAL(wp) :: inv_zinv_i(1:n_zlev)
! ! REAL(wp) :: inv_zinv_m(1:n_zlev)
! ! REAL(wp) :: z_e1(nproma,n_zlev,p_patch%nblks_e)
! ! ! CHARACTER(len=max_char_length), PARAMETER :: &
! ! !        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
! ! !-----------------------------------------------------------------------
! ! i_startblk = p_patch%edges%start_blk(1,1)
! ! i_endblk   = p_patch%edges%end_blk(min_rledge,1)
! ! slev       = 1
! ! 
! ! ipl_src=5  ! output print level (1-5, fix)
! ! z_e1(:,1,:)=top_bc_vn(:,:)
! ! CALL print_mxmn('IMPL VEL: top bc',1,z_e1(:,:,:),1,p_patch%nblks_e,'dif',ipl_src)
! ! z_e1(:,1,:)=bot_bc_vn(:,:)
! ! CALL print_mxmn('IMPL VEL: bot bc',1,z_e1(:,:,:),1,p_patch%nblks_e,'dif',ipl_src)
! ! !write(*,*)'impl vert v-diff: max/min top bc vel:',maxval(top_bc_vn), minval(top_bc_vn)
! ! !write(*,*)'impl vert v-diff: max/min bot bc vel:',maxval(bot_bc_vn), minval(bot_bc_vn)
! ! 
! ! !gam(1:n_zlev)            = 0.0_wp
! ! a(slev:n_zlev)           = 0.0_wp 
! ! b(slev:n_zlev)           = 0.0_wp 
! ! c(slev:n_zlev)           = 0.0_wp
! ! !bet(slev:n_zlev)         = 1.0_wp
! ! inv_zinv_i(slev:n_zlev)  = 0.0_wp
! ! inv_zinv_m(slev:n_zlev)  = 0.0_wp
! ! diff_column=field_column
! ! z_e1=field_column
! ! DO jb = i_startblk, i_endblk
! !   CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
! !                        i_startidx, i_endidx, 1,min_rledge)
! !   DO jc = i_startidx, i_endidx
! ! 
! !     z_dolic = v_base%dolic_e(jc,jb)
! !     IF ( v_base%lsm_oce_e(jc,1,jb) < sea_boundary ) THEN
! !        IF ( z_dolic >=MIN_DOLIC ) THEN
! ! !      IF ( z_dolic >=MIN_DOLIC +1) THEN
! !         inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
! !         inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)
! ! 
! !         !Fill triangular matrix
! !         !b is diagonal a and c are upper and lower band
! !         DO jk = slev+1, z_dolic-1
! !           a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)*dtime
! !           c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1) *dtime
! !           b(jk) = 1.0_wp-a(jk)-c(jk)
! !         END DO
! ! 
! !         ! The first row
! !          c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)*dtime
! !          a(slev) = 0.0_wp           
! !          b(slev) = 1.0_wp- c(slev) !- a(slev) 
! ! 
! !         ! The last row
! !         a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic) *dtime
! !         c(z_dolic) = 0.0_wp
! !         b(z_dolic) = 1.0_wp - a(z_dolic)! - c(z_dolic)
! ! 
! !         ! The matrix is now complete, fill the rhs 
! !         ! The first row contains surface forcing, the last bottom boundary condition
! !         !These codelines are different to the homogeneous version of this sbr (see below) 
! !         field_column(jc,slev,jb)   = field_column(jc,slev,jb)    &
! !                                    &+ top_bc_vn(jc,jb)*inv_zinv_m(slev)*dtime
! !         field_column(jc,z_dolic,jb)= field_column(jc,z_dolic,jb) &
! !                                    &- bot_bc_vn(jc,jb)*inv_zinv_m(z_dolic)*dtime
! ! !----------------------------------------------------------- 
! ! ! !         bet(slev) = 1.0_wp/b(slev)
! ! ! !         DO jk=slev+1, z_dolic
! ! ! !           gam(jk) = a(jk-1)*bet(jk-1)
! ! ! !           IF((b(jk) - c(jk)*gam(jk))/=0.0_wp)THEN 
! ! ! !             bet(jk) = 1.0_wp/(b(jk) - c(jk)*gam(jk))
! ! ! !           ENDIF
! ! ! !         END DO
! ! ! !         diff_column(jc,slev,jb)=field_column(jc,slev,jb)*bet(slev)
! ! ! !         DO jk=slev+1, z_dolic
! ! ! !           diff_column(jc,jk,jb)=bet(jk)*&
! ! ! !           &(field_column(jc,jk,jb)-c(jk)*diff_column(jc,jk-1,jb))
! ! ! !         END DO
! ! ! !         ! Backward sweep
! ! ! !         DO jk=z_dolic-1,slev,-1
! ! ! !           diff_column(jc,jk,jb) = diff_column(jc,jk,jb)  &
! ! ! !           &            -gam(jk+1)*diff_column(jc,jk+1,jb)
! ! ! !         END DO
! ! ! !------------------------------------------
! ! !         DO jk=slev, z_dolic-2
! ! !           IF(b(jk)/=0.0_wp)THEN
! ! !             a(jk) = a(jk)/b(jk)
! ! !             c(jk) = c(jk)/b(jk)
! ! !             field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! ! !             b(jk)=1.0_wp
! ! !           ENDIF
! ! !         END DO
! ! !         DO jk=slev+1, z_dolic-2
! ! !           b(jk)=b(jk)-a(jk)*c(jk-1)
! ! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)&
! ! !                         &-a(jk)*field_column(jc,jk-1,jb)
! ! !           c(jk)=c(jk)/b(jk)
! ! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! ! !           b(jk)=1.0_wp
! ! !         END DO
! ! !         z_tmp=b(z_dolic-1)-a(z_dolic-1)*c(z_dolic-2)
! ! !         z_tmp=(field_column(jc,z_dolic-1,jb)-a(z_dolic-1)*field_column(jc,z_dolic-2,jb))/z_tmp
! ! ! 
! ! !         field_column(jc,z_dolic-1,jb)=z_tmp
! ! !         DO jk=z_dolic-2,1,-1
! ! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
! ! !         END DO
! ! !         DO jk=1,z_dolic-1
! ! !           diff_column(jc,jk,jb)=field_column(jc,jk,jb)
! ! !         END DO
! ! !----------------------------------------------------------- 
! !         DO jk=slev, z_dolic-1
! !           IF(b(jk)/=0.0_wp)THEN
! !             a(jk) = a(jk)/b(jk)
! !             c(jk) = c(jk)/b(jk)
! !             field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! !             b(jk)=1.0_wp
! !           ENDIF
! !         END DO
! ! 
! !         DO jk=slev+1, z_dolic-1
! !           b(jk)=b(jk)-a(jk)*c(jk-1)
! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)&
! !                         &-a(jk)*field_column(jc,jk-1,jb)
! !           c(jk)=c(jk)/b(jk)
! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
! !           b(jk)=1.0_wp
! !         END DO
! ! 
! !         z_tmp=b(z_dolic)-a(z_dolic)*c(z_dolic-1)
! !         z_tmp=(field_column(jc,z_dolic,jb)-a(z_dolic)*field_column(jc,z_dolic-1,jb))/z_tmp
! ! 
! !         field_column(jc,z_dolic,jb)=z_tmp
! !         DO jk=z_dolic-1,1,-1
! !           field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
! !         END DO
! !         DO jk=1,z_dolic
! ! !       DO jk=1,z_dolic-1
! !           diff_column(jc,jk,jb)=field_column(jc,jk,jb)
! !         END DO
! ! 
! ! !      IF(field_column(jc,1,jb)/=0.0_wp)THEN
! ! !      !write(234,*)'top bc',top_bc(jc,jb)/zinv
! ! !      write(234,*)'coffs',jc,jb,A_v(jc,:,jb)
! ! !      write(234,*)'mat up  a:', a
! ! !      write(234,*)'mat dia b:', b
! ! !      write(234,*)'mat dow c:', c
! ! !       write(234,*)'in :', z_e1(jc,:,jb) 
! ! !       write(234,*)'out:',diff_column(jc,:,jb)
! ! ! !      write(234,*)'sum:',sum(field_column(jc,:,jb))/z_dolic,&
! ! ! !      & sum(diff_column(jc,:,jb))/z_dolic, &
! ! ! !      &(sum(diff_column(jc,:,jb))/z_dolic)/(sum(field_column(jc,:,jb))/z_dolic)
! ! ! !     write(234,*)
! ! !      ENDIF
! ! 
! !       ELSEIF ( z_dolic <MIN_DOLIC ) THEN
! !         diff_column(jc,:,jb) = 0.0_wp
! !         field_column(jc,:,jb)= 0.0_wp
! !       ENDIF
! !     ELSEIF( v_base%lsm_oce_e(jc,1,jb) >= sea_boundary ) THEN
! !       diff_column(jc,:,jb) = 0.0_wp
! !       field_column(jc,:,jb)= 0.0_wp
! !     ENDIF
! !   END DO
! ! END DO
! ! 
! ! DO jk=slev, n_zlev
! !   ipl_src=5  ! output print level (1-5, fix)
! !   CALL print_mxmn('IMPL VEL: aft.vvel.dif',jk,diff_column(:,:,:),n_zlev, &
! !     &              p_patch%nblks_c,'dif',ipl_src)
! ! END DO
! ! 
! ! END subroutine veloc_diffusion_vert_impl
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
SUBROUTINE veloc_diffusion_vert_impl_hom( p_patch,       &
                                    & field_column,  &
                                    & h_e,           &
                                    & A_v,           &
                                    & diff_column)
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(inout)           :: field_column(:,:,:)
REAL(wp), INTENT(IN)              :: h_e(:,:)           !surface height at edges, relevant for thickness of first cell 
REAL(wp), INTENT(inout)           :: A_v(:,:,:) 
REAL(wp), INTENT(out)             :: diff_column(:,:,:)
!
!Local variables
INTEGER :: slev
INTEGER :: jc, jk, jb
INTEGER :: i_startidx, i_endidx
REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
!REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
REAL(wp) :: z_tmp
INTEGER  :: z_dolic
REAL(wp) :: inv_zinv_i(1:n_zlev)
REAL(wp) :: inv_zinv_m(1:n_zlev)
    
TYPE(t_subset_range), POINTER :: all_edges
!REAL(wp) :: z_e1(nproma,1,p_patch%nblks_e)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!-----------------------------------------------------------------------
all_edges => p_patch%edges%all

slev       = 1

!gam(1:n_zlev)     = 0.0_wp
a(slev:n_zlev)    = 0.0_wp 
b(slev:n_zlev)    = 0.0_wp 
c(slev:n_zlev)    = 0.0_wp
!bet(slev:n_zlev)  = 1.0_wp
inv_zinv_i(slev:n_zlev)  = 0.0_wp
inv_zinv_m(slev:n_zlev)  = 0.0_wp
diff_column=field_column

DO jb = all_edges%start_block, all_edges%end_block
  CALL get_index_range(all_edges, jb, i_startidx, i_endidx)
  DO jc = i_startidx, i_endidx

    z_dolic             = v_base%dolic_e(jc,jb)
    IF ( v_base%lsm_oce_e(jc,1,jb) < sea_boundary ) THEN
      IF ( z_dolic >MIN_DOLIC ) THEN

        inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
        inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)


        !Fill triangular matrix
        !b is diagonal a and c are upper and lower band
        DO jk = slev+1, z_dolic-1
          a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)*dtime
          c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)*dtime
          b(jk) = 1.0_wp-a(jk)-c(jk)
        END DO

        ! The first row
         c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)*dtime
         a(slev) = 0.0_wp           
         b(slev) = 1.0_wp- c(slev) !- a(slev) 

        ! The last row
        a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)
        c(z_dolic) = 0.0_wp
        b(z_dolic) = 1.0_wp - a(z_dolic)! - c(z_dolic)

!------------------------------------------
! !         bet(slev) = 1.0_wp/b(slev)
! !         DO jk=slev+1, z_dolic
! !           gam(jk) = a(jk-1)*bet(jk-1)
! !           IF((b(jk) - c(jk)*gam(jk))/=0.0_wp)THEN 
! !             bet(jk) = 1.0_wp/(b(jk) - c(jk)*gam(jk))
! !           ENDIF
! !         END DO
! !         diff_column(jc,slev,jb)=field_column(jc,slev,jb)!*bet(slev)
! !         DO jk=slev+1, z_dolic
! !           diff_column(jc,jk,jb)=bet(jk)*&
! !           &(field_column(jc,jk,jb)-c(jk)*diff_column(jc,jk-1,jb))
! !         END DO
! !         ! Backward sweep
! !         DO jk=z_dolic-1,slev,-1
! !           diff_column(jc,jk,jb) = diff_column(jc,jk,jb)  &
! !           &            -gam(jk+1)*diff_column(jc,jk+1,jb)
! !         END DO
!------------------------------------------
        DO jk=slev, z_dolic-2
          IF(b(jk)/=0.0_wp)THEN
            a(jk) = a(jk)/b(jk)
            c(jk) = c(jk)/b(jk)
            field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
            b(jk)=1.0_wp
          ENDIF
        END DO
        DO jk=slev+1, z_dolic-2
          b(jk)=b(jk)-a(jk)*c(jk-1)
          field_column(jc,jk,jb)=field_column(jc,jk,jb)&
                        &-a(jk)*field_column(jc,jk-1,jb)
          c(jk)=c(jk)/b(jk)
          field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
          b(jk)=1.0_wp
        END DO
        z_tmp=b(z_dolic-1)-a(z_dolic-1)*c(z_dolic-2)
        z_tmp=(field_column(jc,z_dolic-1,jb)-a(z_dolic-1)*field_column(jc,z_dolic-2,jb))/z_tmp

        field_column(jc,z_dolic-1,jb)=z_tmp
        DO jk=z_dolic-2,1,-1
          field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
        END DO
        DO jk=1,z_dolic
    !   DO jk=1,z_dolic-1
          diff_column(jc,jk,jb)=field_column(jc,jk,jb)
        END DO
!----------------------------------------------------------- 
!         DO jk=slev, z_dolic-1
!           IF(b(jk)/=0.0_wp)THEN
!             a(jk) = a(jk)/b(jk)
!             c(jk) = c(jk)/b(jk)
!             field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
!             b(jk)=1.0_wp
!           ENDIF
!         END DO
! 
!         DO jk=slev+1, z_dolic-1
!           b(jk)=b(jk)-a(jk)*c(jk-1)
!           field_column(jc,jk,jb)=field_column(jc,jk,jb)&
!                         &-a(jk)*field_column(jc,jk-1,jb)
!           c(jk)=c(jk)/b(jk)
!           field_column(jc,jk,jb)=field_column(jc,jk,jb)/b(jk)
!           b(jk)=1.0_wp
!         END DO
! 
!         z_tmp=b(z_dolic)-a(z_dolic)*c(z_dolic-1)
!         z_tmp=(field_column(jc,z_dolic,jb)-a(z_dolic)*field_column(jc,z_dolic-1,jb))/z_tmp
! 
!         field_column(jc,z_dolic,jb)=z_tmp
!         DO jk=z_dolic-1,1,-1
!           field_column(jc,jk,jb)=field_column(jc,jk,jb)-c(jk)*field_column(jc,jk+1,jb)
!         END DO
!         DO jk=1,z_dolic-1
!           diff_column(jc,jk,jb)=field_column(jc,jk,jb)
!         END DO




      ELSEIF ( z_dolic <MIN_DOLIC ) THEN
        diff_column(jc,:,jb) = 0.0_wp
        field_column(jc,:,jb)= 0.0_wp
      ENDIF
    ELSEIF( v_base%lsm_oce_e(jc,1,jb) >= sea_boundary ) THEN
      diff_column(jc,:,jb) = 0.0_wp
      field_column(jc,:,jb)= 0.0_wp
    ENDIF
  END DO
END DO

DO jk=slev, n_zlev
  ipl_src=5  ! output print level (1-5, fix)
  CALL print_mxmn('IMPL VEL: aft.vvel.dif',jk,diff_column(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'dif',ipl_src)
END DO



END subroutine veloc_diffusion_vert_impl_hom
!------------------------------------------------------------------------
END MODULE mo_oce_diffusion
