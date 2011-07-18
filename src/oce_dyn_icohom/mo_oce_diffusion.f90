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
USE mo_kind,                      ONLY: wp
!USE mo_mpi,                       ONLY: p_pe, p_io
USE mo_physical_constants
USE mo_math_utilities,            ONLY: t_cartesian_coordinates, gvec2cvec!, gc2cc
USE mo_impl_constants,            ONLY: boundary, sea_boundary ,max_char_length, &
  &                                     min_rlcell, min_rledge, min_rlcell,      &
  &                                     max_char_length
USE mo_parallel_config,  ONLY: nproma
USE mo_ocean_nml,                 ONLY: n_zlev, iswm_oce
USE mo_run_config,                ONLY: dtime
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, t_hydro_ocean_aux
USE mo_model_domain,              ONLY: t_patch
!USE mo_exception,                 ONLY: message, finish!, message_text
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
!USE mo_oce_boundcond,             ONLY: bot_bound_cond_horz_veloc, top_bound_cond_horz_veloc!,&
 ! &                                     bot_bound_cond_vert_veloc, top_bound_cond_vert_veloc
USE mo_oce_physics
USE mo_scalar_product,            ONLY: map_cell2edges, primal_map_c2e
USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d, &
  &                                     nabla2_vec_ocean
USE mo_interpolation,          ONLY: t_int_state
!USE mo_oce_index,              ONLY: c_i, c_b, c_k, ne_b, ne_i, nc_b, nc_i, form4ar, ldbg

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
PUBLIC :: tracer_diffusion_vert_impl
PUBLIC :: veloc_diffusion_vert_impl
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
SUBROUTINE velocity_diffusion_horz_mimetic(p_patch, vn_in, p_param, p_diag, laplacian_vn_out)
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch%nblks_e)
TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
TYPE(t_hydro_ocean_diag)          :: p_diag
REAL(wp), INTENT(out)             :: laplacian_vn_out(nproma,n_zlev,p_patch%nblks_e)
!
!Local variables
INTEGER :: slev, elev
INTEGER :: jk, jb, je,jc
INTEGER :: il_c1, ib_c1, il_c2, ib_c2
INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,n_zlev,p_patch%nblks_e)
TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        
rl_start_e = 1
rl_end_e  = min_rledge

rl_start_c = 1
rl_end_c  = min_rlcell

i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

slev = 1
elev = n_zlev

!Step 1: Calculate gradient of cell velocity vector.
!Result is a gradient vector, located at edges
!Step 2: Multiply each component of gradient vector with mixing coefficients
DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e,&
                   &  i_startidx_e, i_endidx_e,&
                   &  rl_start_e, rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e

      !Get indices of two adjacent triangles
      il_c1 = p_patch%edges%cell_idx(je,jb,1)
      ib_c1 = p_patch%edges%cell_blk(je,jb,1)
      il_c2 = p_patch%edges%cell_idx(je,jb,2)
      ib_c2 = p_patch%edges%cell_blk(je,jb,2)

    IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
      z_grad_u(je,jk,jb)%x = p_param%K_veloc_h(je,jk,jb)&
        &                  *(p_diag%p_vn(il_c2,jk,ib_c2)%x &
        &                  - p_diag%p_vn(il_c1,jk,ib_c1)%x)              &
        &                  / p_patch%edges%dual_edge_length(je,jb)
    ELSE
      z_grad_u(je,jk,jb)%x = 0.0_wp
    ENDIF 

    ENDDO
  END DO
END DO

!Step 2: Apply divergence to each component of mixing times gradient vector
iidx => p_patch%cells%edge_idx
iblk => p_patch%cells%edge_blk

DO jb = i_startblk_c, i_endblk_c

  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
                     i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = slev, elev
      DO jc = i_startidx_c, i_endidx_c

         IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) >= boundary ) THEN
           z_div_grad_u(jc,jk,jb)%x = 0.0_wp
         ELSE
          z_div_grad_u(jc,jk,jb)%x =  &
            z_grad_u(iidx(jc,jb,1),jk,iblk(jc,jb,1))%x * p_patch%patch_oce%geofac_div(jc,1,jb) + &
            z_grad_u(iidx(jc,jb,2),jk,iblk(jc,jb,2))%x * p_patch%patch_oce%geofac_div(jc,2,jb) + &
            z_grad_u(iidx(jc,jb,3),jk,iblk(jc,jb,3))%x * p_patch%patch_oce%geofac_div(jc,3,jb)
        ENDIF

      END DO
    END DO
  END DO

!Step 3: Map divergence back to edges
CALL map_cell2edges( p_patch, z_div_grad_u, laplacian_vn_out)

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
SUBROUTINE velocity_diffusion_horz_rbf(p_patch, vn_in, p_param, p_diag, laplacian_vn_out)
!
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
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e
INTEGER :: il_c1, ib_c1, il_c2, ib_c2
!REAL(wp) :: z_mixing_coeff
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_ocean_semi_implicit_ab:velocity_diffusion_horz')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        
rl_start_e = 1
rl_end_e  = min_rledge

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

slev = 1
elev = n_zlev

CALL nabla2_vec_ocean( vn_in,&
                    & p_diag%vt,&
                    & p_patch, &
                    &laplacian_vn_out,&
                    &opt_slev=slev,opt_elev=elev )

!Step 2: Multiply laplacian with mixing coefficients
DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e,&
                   &  i_startidx_e, i_endidx_e,&
                   &  rl_start_e, rl_end_e)
  DO jk = slev, elev
    DO je = i_startidx_e, i_endidx_e

      il_c1 = p_patch%edges%cell_idx(je,jb,1)
      ib_c1 = p_patch%edges%cell_blk(je,jb,1)
      il_c2 = p_patch%edges%cell_idx(je,jb,2)
      ib_c2 = p_patch%edges%cell_blk(je,jb,2)

      IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
        laplacian_vn_out(je,jk,jb) = p_param%K_veloc_h(je,jk,jb)*laplacian_vn_out(je,jk,jb)
      ELSE
        laplacian_vn_out(je,jk,jb) = 0.0_wp
      ENDIF
    ENDDO
  END DO
END DO

END SUBROUTINE velocity_diffusion_horz_rbf
! !-------------------------------------------------------------------------
!
!
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
SUBROUTINE velocity_diffusion_vert_mimetic( p_patch, p_diag, p_aux,h_c,p_param, laplacian_vn_out)
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
TYPE(t_hydro_ocean_diag)          :: p_diag
TYPE(t_hydro_ocean_aux)           :: p_aux
REAL(wp), INTENT(in)              :: h_c(:,:) 
TYPE(t_ho_params), INTENT(in)     :: p_param
REAL(wp)                          :: laplacian_vn_out(:,:,:)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb, z_dolic
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &                              z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:veloc diffusion vert mimetic')
!-----------------------------------------------------------------------
z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c)%x = 0.0_wp
z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)%x   = 0.0_wp

slev       = 1
elev       = n_zlev
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)

!1 Vertical derivative of cell velocity vector times horizontal velocity
! loop runs now from slev to z_dolic:
!DO jk = slev, elev
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx
    !check if we have at least two layers of water
    !  #slo# - 2011-04-01 - Is this really intended here
    !  #slo# - 2011-05-25 - with z_dolic and update in fill_vertical_domain now ok. 
    !IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    IF (z_dolic > 0) THEN
      !1a) ocean surface 
      !IF(jk==slev)THEN
      jk = slev
      z_adv_u_i(jc,jk,jb)%x =                              &
      &( p_aux%bc_top_veloc_cc(jc,jb)%x                    &
      &- p_param%A_veloc_v(jc,jk,jb)                       &
      &*(p_diag%p_vn(jc,jk,jb)%x-p_diag%p_vn(jc,jk+1,jb)%x)&
      &/p_patch%patch_oce%del_zlev_i(jk))/p_patch%patch_oce%del_zlev_m(jk)
      !1b) ocean bottom 
      !ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
      !ELSEIF ( jk == z_dolic ) THEN
      jk = z_dolic
!        IF(p_patch%patch_oce%n_zlev>=2)THEN
      z_adv_u_i(jc,jk,jb)%x&
      & = (p_param%A_veloc_v(jc,jk,jb)&
      &*(p_diag%p_vn(jc,jk-1,jb)%x-p_diag%p_vn(jc,jk,jb)%x) &
      & / p_patch%patch_oce%del_zlev_i(jk-1)&
      & - p_aux%bc_bot_veloc_cc(jc,jb)%x)/p_patch%patch_oce%del_zlev_m(jk)
      !write(*,*)'u-diff botom:',jk,u_c(jc,elev-1,jb),u_c(jc,elev,jb),u_c(jc,elev,jb), &
      !  &        bot_bc_u_c(jc,jb), z_adv_u_i(jc,jk,jb)
!        ENDIF
      !1c) ocean interior 
      !ELSEIF( jk>slev .AND. jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN
      DO jk = slev+1, z_dolic-1
         z_adv_u_i(jc,jk,jb)%x&
         & = ( p_param%A_veloc_v(jc,jk-1,jb)*(p_diag%p_vn(jc,jk-1,jb)%x-p_diag%p_vn(jc,jk,jb)%x)&
         &   /p_patch%patch_oce%del_zlev_i(jk-1)&
         &  - p_param%A_veloc_v(jc,jk,jb)*(p_diag%p_vn(jc,jk,jb)%x  -p_diag%p_vn(jc,jk+1,jb)%x)&
         &   /p_patch%patch_oce%del_zlev_i(jk))/&
         &    p_patch%patch_oce%del_zlev_m(jk)
! write(*,*)'u-diff interior:',jk,u_c(jc,jk-1,jb),u_c(jc,jk,jb),u_c(jc,jk,jb),u_c(jc,jk+1,jb), &
!   &        z_adv_u_i(jc,jk,jb), p_patch%patch_oce%dolic_c(jc,jb) 
      ! ENDIF  ! jk-condition
     !ENDIF    ! at least 2 vertical layers
      END DO ! jk ocean interior
      ENDIF  ! dolic>0
    END DO
  END DO
!END DO

! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)  
CALL map_cell2edges( p_patch, z_adv_u_i, laplacian_vn_out)

! DO jk=slev, elev
!   WRITE(*,*)'MAX/MIN vert diffusion ',jk, &
!     &        MAXVAL(laplacian_vn_out(:,jk,:)), MINVAL(laplacian_vn_out(:,jk,:))
! END DO

END subroutine velocity_diffusion_vert_mimetic
!-------------------------------------------------------------------------  
!
!
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
REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_v_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_u_m(nproma,n_zlev,p_patch%nblks_c),  &
  &         z_adv_v_m(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:veloc diffusion vert rbf')
!-----------------------------------------------------------------------

! #slo# set local variable to zero due to nag -nan compiler-option
z_adv_u_i(:,:,:) = 0.0_wp
z_adv_v_i(:,:,:) = 0.0_wp
z_adv_u_m(:,:,:) = 0.0_wp
z_adv_v_m(:,:,:) = 0.0_wp

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
    IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
      !1a) 0cean surface
      !check if we have at least two layers of water
     !IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
      IF(jk==slev)THEN
          ! u-component
          z_adv_u_i(jc,jk,jb) =      &
          &( top_bc_u_c(jc,jb)       &
          & -p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk,jb)-u_c(jc,jk+1,jb)) &
          & /p_patch%patch_oce%del_zlev_i(jk))/p_patch%patch_oce%del_zlev_m(jk)
          ! v-component
          z_adv_v_i(jc,jk,jb) =      &
          &( top_bc_v_c(jc,jb)       &
          & -p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
          &/p_patch%patch_oce%del_zlev_i(jk)) &
          & /p_patch%patch_oce%del_zlev_m(jk)
          !write(*,*)'RBD:vert-diff:top',jc,jk,jb, z_adv_u_i(jc,jk,jb),top_bc_u_c(jc,jb)

      !1b) ocean bottom 
      ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN

        IF(p_patch%patch_oce%n_zlev>=2)THEN
          ! u-component
          z_adv_u_i(jc,jk,jb)              &
          & = (p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk-1,jb)-u_c(jc,jk,jb))&
          & /p_patch%patch_oce%del_zlev_i(jk-1)&
          & - bot_bc_u_c(jc,jb))/p_patch%patch_oce%del_zlev_m(jk)
          ! v-component
          z_adv_v_i(jc,jk,jb)&
          & = (p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,jk-1,jb)-v_c(jc,jk,jb))&
          & /p_patch%patch_oce%del_zlev_i(jk-1)&
          & -bot_bc_v_c(jc,jb))/p_patch%patch_oce%del_zlev_m(jk)
        !write(*,*)'u-diff botom:',jk,u_c(jc,elev-1,jb),u_c(jc,elev,jb),u_c(jc,elev,jb),bot_bc_u_c(jc,jb), z_adv_u_i(jc,jk,jb)
        ENDIF

      !1c) ocean interior 
      ELSEIF(jk>slev .AND. jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN
        ! u-component
        z_adv_u_i(jc,jk,jb)&
        & = &
        & ( p_param%A_veloc_v(jc,jk-1,jb)*(u_c(jc,jk-1,jb)-u_c(jc,jk,jb))&
        &/p_patch%patch_oce%del_zlev_i(jk-1)&
        & - p_param%A_veloc_v(jc,jk,jb)  *(u_c(jc,jk,jb)-u_c(jc,jk+1,jb))&
        &/p_patch%patch_oce%del_zlev_i(jk)) &
        & / p_patch%patch_oce%del_zlev_m(jk)
        ! v-component
        z_adv_v_i(jc,jk,jb)&
        & = &
        & ( p_param%A_veloc_v(jc,jk-1,jb)*(v_c(jc,jk-1,jb)-v_c(jc,jk,jb))&
        &/p_patch%patch_oce%del_zlev_i(jk-1)&
        & - p_param%A_veloc_v(jc,jk,jb)  *(v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
        &/p_patch%patch_oce%del_zlev_i(jk)) &
        & / p_patch%patch_oce%del_zlev_m(jk)

! write(*,*)'u-diff interior:',jk,u_c(jc,jk-1,jb),u_c(jc,jk,jb),u_c(jc,jk,jb),u_c(jc,jk+1,jb), z_adv_u_i(jc,jk,jb),&
! &p_patch%patch_oce%dolic_c(jc,jb) 
      !1c) ocean bottom 

      ENDIF  ! jk-condition
      ENDIF
      CALL gvec2cvec( z_adv_u_i(jc,jk,jb), z_adv_v_i(jc,jk,jb),       &
      &               p_patch%cells%center(jc,jb)%lon,                &
      &               p_patch%cells%center(jc,jb)%lat,                &
      &               zu_cc(jc,jk,jb)%x(1),zu_cc(jc,jk,jb)%x(2),zu_cc(jc,jk,jb)%x(3))

    !ENDIF   ! at least 2 vertical layers
    END DO
  END DO
END DO

!  DO jk=slev, elev
!    WRITE(*,*)'MAX/MIN vert diff A:',jk, &
!      &        MAXVAL(z_adv_u_i(:,jk,:)), MINVAL(z_adv_u_i(:,jk,:)),&
!      &        MAXVAL(z_adv_v_i(:,jk,:)), MINVAL(z_adv_v_i(:,jk,:))
!  END DO

! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)
CALL map_cell2edges( p_patch, zu_cc, laplacian_vn_out)
 DO jk=slev, elev
   WRITE(*,*)'MAX/MIN vert diffusion ',jk, &
     &        MAXVAL(laplacian_vn_out(:,jk,:)), MINVAL(laplacian_vn_out(:,jk,:))
 END DO

END subroutine velocity_diffusion_vert_rbf
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

      delta_z  = p_patch%patch_oce%del_zlev_m(jk)

      DO je = i_startidx_e, i_endidx_e

        IF(jk==slev) delta_z = p_patch%patch_oce%del_zlev_m(slev) + p_os%p_diag%h_e(je,jb)

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
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

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

         delta_z  = p_os%p_diag%thick_e(je,jb)

          diff_flx(je,jk,jb) = K_T(je,jk,jb)*delta_z&
                     &*(trac_in(il_c2,jk,ib_c2)-trac_in(il_c1,jk,ib_c1))&
                     &/p_patch%edges%dual_edge_length(je,jb)
write(123,*)'trac diffusive flux',je,jk,jb,diff_flx(je,jk,jb),K_T(je,jk,jb),&
&trac_in(il_c2,jk,ib_c2),trac_in(il_c1,jk,ib_c1)
        ELSE
          diff_flx(je,jk,jb) = 0.0_wp
        ENDIF
      ENDDO
    END DO
  END DO
ENDIF
! Apply divergence to mixing times gradient to get laplacian
!CALL div_oce( z_grad, p_patch, laplacian_trac_out)

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
                                    & diff_flx)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(inout)           :: trac_c(:,:,:)
REAL(wp), INTENT(in)              :: dz(:,:,:)
REAL(wp), INTENT(in)              :: top_bc_tracer(:,:)
REAL(wp), INTENT(in)              :: bot_bc_tracer(:,:)
REAL(wp), INTENT(in)              :: A_v(:,:,:) 
REAL(wp), INTENT(out)             :: diff_flx(:,:,:)
!
!Local variables
INTEGER :: slev, elev
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:tracer_diffusion_vert')
!-----------------------------------------------------------------------
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev
diff_flx(:,:,:) = 0.0_wp
!1 Vertical derivative of tracer
DO jk = slev, elev+1
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
    DO jc = i_startidx, i_endidx

      !1a) 0cean surface
      IF(jk==slev)THEN

        IF (p_patch%patch_oce%dolic_c(jc,jb) /= 0) THEN
!           diff_flx(jc,jk,jb) = A_v(jk)&
!           & *(top_bc_tracer(jc,jb)-trac_c(jc,jk,jb))/p_patch%patch_oce%del_zlev_i(jk)
            diff_flx(jc,jk,jb) = top_bc_tracer(jc,jb)
        ENDIF

      !1b) ocean bottom 
      ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb)+1 ) THEN

        IF(p_patch%patch_oce%n_zlev>=2)THEN
          diff_flx(jc,jk,jb)= bot_bc_tracer(jc,jb)
        ENDIF
      !1c) ocean interior 
      ELSEIF ( jk>slev .AND.  jk <= p_patch%patch_oce%dolic_c(jc,jb) ) THEN
        IF(dz(jc,jk,jb)/=0.0_wp)THEN
          diff_flx(jc,jk,jb)&
          & = A_v(jc,jk,jb) &
          & * (trac_c(jc,jk-1,jb)-trac_c(jc,jk,jb))/dz(jc,jk,jb)! p_patch%patch_oce%del_zlev_i(jk)
        ELSE
          diff_flx(jc,jk,jb)= 0.0_wp
        ENDIF
      ENDIF
    END DO
  END DO
END DO
! DO jk=slev, elev
!   WRITE(*,*)'MAX/MIN vert tracer diffusion ',jk, &
!     &        MAXVAL(diff_flx(:,jk,:)), MINVAL(diff_flx(:,jk,:))
! END DO

END subroutine tracer_diffusion_vert_expl
!-------------------------------------------------------------------------  
!
!!Subroutine implements implicit vertical diffusion for scalar fields.
!>
!!
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
SUBROUTINE tracer_diffusion_vert_impl( p_patch,      &
                                & field_column,   &
                                & top_bc,         & 
                                & bot_bc,         &
                                & A_v,            &
                                & diff_column)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(inout)           :: field_column(:,:,:)
REAL(wp), INTENT(in)              :: top_bc(:,:)
REAL(wp), INTENT(in)              :: bot_bc(:,:)  
REAL(wp), INTENT(in)              :: A_v(:,:,:) 
REAL(wp), INTENT(out)             :: diff_column(:,:,:)
!
!Local variables
INTEGER :: slev, elev
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
REAL(wp) :: z_trac_rhs(1:n_zlev)
REAL(wp) :: Z(1:n_zlev), zbar(1:n_zlev+1)
REAL(wp) :: dt_inv, t1, zinv
INTEGER  :: z_dolic
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!-----------------------------------------------------------------------
i_startblk = p_patch%cells%start_blk(1,1)
i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
slev = 1
elev = n_zlev

dt_inv = 1.0_wp/dtime

diff_column(:,:,:)= 0.0_wp
a(slev:n_zlev)    = 0.0_wp 
b(slev:n_zlev)    = 0.0_wp 
c(slev:n_zlev)    = 0.0_wp
Z(slev:n_zlev)    = 0.0_wp 
zbar(slev:n_zlev) = 0.0_wp
z_trac_rhs(slev:n_zlev)= 0.0_wp

DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rlcell)
  DO jc = i_startidx, i_endidx

    z_dolic             = p_patch%patch_oce%dolic_c(jc,jb)
    Z(slev:n_zlev)      = p_patch%patch_oce%zlev_m(slev:n_zlev)
    zbar(slev:n_zlev+1) = p_patch%patch_oce%zlev_i(slev:n_zlev+1)

  IF ( z_dolic >=3 ) THEN
    !Fill triangular matrix
    !b is diagonal a nd c are lower and upper band
    DO jk = slev+1, z_dolic-1  
      zinv  = 1.0_wp/(zbar(jk)-zbar(jk+1))     
      a(jk) = -A_v(jc,jk,jb)/(Z(jk-1)-Z(jk))*zinv
      c(jk) = -A_v(jc,jk+1,jb)/(Z(jk)-Z(jk+1))*zinv
      b(jk) = -a(jk)-c(jk)+dt_inv
    END DO

    ! The first row
    zinv = 1.0_wp/(zbar(1)-zbar(2))
    c(1) = -A_v(jc,2,jb)/(Z(1)-Z(2))*zinv
    a(1) = top_bc(jc,jb)
    b(1) = -c(1)+dt_inv

    ! The last row
    zinv       = 1.0_wp/(zbar(z_dolic-1)-zbar(z_dolic))
    a(z_dolic) = -A_v(jc,z_dolic,jb)/(Z(z_dolic-1)-Z(z_dolic-2))*zinv
    b(z_dolic) = -a(z_dolic)+dt_inv
    c(z_dolic) = bot_bc(jc,jb)

    ! The matrix is now complete, fill the rhs 
    z_trac_rhs(2:z_dolic)=field_column(jc,2:z_dolic,jb)*dt_inv

    !  The first row contains also surface forcing
    z_trac_rhs(1)=dt_inv*field_column(jc,1,jb) !-1.0e-6*heat_flux(n)/(4.2_wp)+10.0_wp*alpha*(Tsurf(n)-TF(1,n)) 

    ! =============================================
    ! Prepare coefficients and rhs (scale with diagonal element)
    DO jk=slev, z_dolic-1
      a(jk)         = a(jk)/b(jk)
      c(jk)         = c(jk)/b(jk)
      z_trac_rhs(jk)= z_trac_rhs(jk)/b(jk)
      b(jk)         = 1.0_wp
    END DO

    DO jk=slev+1, z_dolic-1
      b(jk)          = b(jk)-a(jk)*c(jk-1)
      z_trac_rhs(jk) = z_trac_rhs(jk)-a(jk)*z_trac_rhs(jk-1)
      c(jk)          = c(jk)/b(jk)
      z_trac_rhs(jk) = z_trac_rhs(jk)/b(jk)
      b(jk)          = 1.0_wp
    END DO

    t1=b(z_dolic)-a(z_dolic)*c(z_dolic-1)
    t1=(z_trac_rhs(z_dolic)-a(z_dolic)*z_trac_rhs(z_dolic-1))/t1

    z_trac_rhs(z_dolic)=t1

    ! Backward sweep
    DO jk=z_dolic-1,1,-1
      z_trac_rhs(jk) = z_trac_rhs(jk)-c(jk)*z_trac_rhs(jk+1)
    END DO

    ! Now update of tracer fields
    DO jk=1,z_dolic
      diff_column(jc,jk,jb) = z_trac_rhs(jk)
    END DO
  ELSE
    diff_column(jc,:,jb) = field_column(jc,:,jb)
  ENDIF
  END DO
END DO
 DO jk=slev, elev
   WRITE(*,*)'MAX/MIN before vert tracer diffusion ',jk, &
     &        MAXVAL(field_column(:,jk,:)), MINVAL(field_column(:,jk,:))
 END DO


 DO jk=slev, elev
   WRITE(*,*)'MAX/MIN after vert tracer diffusion ',jk, &
     &        MAXVAL(diff_column(:,jk,:)), MINVAL(diff_column(:,jk,:))
 END DO
END subroutine tracer_diffusion_vert_impl
!-------------------------------------------------------------------------  
!
!!Subroutine implements implicit vertical diffusion for scalar fields.
!>
!!
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
SUBROUTINE veloc_diffusion_vert_impl( p_patch,       &
                                    & field_column,  &
                                    & p_os_aux,      &
                                    & A_v,           &
                                    & diff_column)
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(inout)           :: field_column(:,:,:)
TYPE(t_hydro_ocean_aux)           :: p_os_aux     
REAL(wp), INTENT(in)              :: A_v(:,:,:) 
REAL(wp), INTENT(out)             :: diff_column(:,:,:)
!
!Local variables
INTEGER :: slev, elev
INTEGER :: je, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
REAL(wp) :: z_trac_rhs(1:n_zlev)
REAL(wp) :: Z(1:n_zlev), zbar(1:n_zlev+1)
REAL(wp) :: dt_inv, zinv
REAL(wp) :: t1
INTEGER  :: z_dolic
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:tracer_diffusion_impl')
!-----------------------------------------------------------------------
i_startblk = p_patch%edges%start_blk(1,1)
i_endblk   = p_patch%edges%end_blk(min_rledge,1)
slev       = 1
elev       = n_zlev

dt_inv = 1.0_wp/dtime

diff_column(:,:,:)= 0.0_wp
a(slev:n_zlev)    = 0.0_wp 
b(slev:n_zlev)    = 0.0_wp 
c(slev:n_zlev)    = 0.0_wp
Z(slev:n_zlev)    = 0.0_wp 
zbar(slev:n_zlev) = 0.0_wp
z_trac_rhs(slev:n_zlev)= 0.0_wp

DO jb = i_startblk, i_endblk
  CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1,min_rledge)
  DO je = i_startidx, i_endidx

    z_dolic             = p_patch%patch_oce%dolic_e(je,jb)
    Z(slev:n_zlev)      = p_patch%patch_oce%zlev_m(slev:n_zlev)
    zbar(slev:n_zlev+1) = p_patch%patch_oce%zlev_i(slev:n_zlev+1)
   
  IF ( z_dolic >=3 ) THEN
  IF ( p_patch%patch_oce%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
    !Fill triangular matrix
    !b is diagonal a nd c are lower and upper band
    DO jk = slev+1, z_dolic-1  
      zinv  = 1.0_wp/(zbar(jk)-zbar(jk+1))     
      a(jk) = -A_v(je,jk,jb)/(Z(jk-1)-Z(jk))*zinv
      c(jk) = -A_v(je,jk+1,jb)/(Z(jk)-Z(jk+1))*zinv
      b(jk) = -a(jk)-c(jk)+dt_inv
    END DO

    ! The first row
    zinv = 1.0_wp/(zbar(1)-zbar(2))
    c(1) = -A_v(je,2,jb)/(Z(1)-Z(2))*zinv
    a(1) = p_os_aux%bc_top_vn(je,jb)
    b(1) = -c(1)+dt_inv

    ! The last row
    zinv       = 1.0_wp/(zbar(z_dolic-1)-zbar(z_dolic))
    a(z_dolic) = -A_v(je,z_dolic,jb)/(Z(z_dolic-1)-Z(z_dolic-2))*zinv
    b(z_dolic) = -a(z_dolic)+dt_inv
    c(z_dolic) = p_os_aux%bc_bot_vn(je,jb)
!write(*,*)'a(1) ',je,jb,a(1) 
    ! The matrix is now complete, fill the rhs 
    z_trac_rhs(2:z_dolic)=field_column(je,2:z_dolic,jb)*dt_inv

    !  The first row contains also surface forcing
    z_trac_rhs(1)=dt_inv*field_column(je,1,jb)

    ! =============================================
    ! Prepare coefficients and rhs (scale with diagonal element)
    DO jk=slev, z_dolic-1
      a(jk)         = a(jk)/b(jk)
      c(jk)         = c(jk)/b(jk)
      z_trac_rhs(jk)= z_trac_rhs(jk)/b(jk)
      b(jk)         = 1.0_wp
    END DO

    DO jk=slev+1, z_dolic-1
      b(jk)          = b(jk)-a(jk)*c(jk-1)
      z_trac_rhs(jk) = z_trac_rhs(jk)-a(jk)*z_trac_rhs(jk-1)
      c(jk)          = c(jk)/b(jk)
      z_trac_rhs(jk) = z_trac_rhs(jk)/b(jk)
      b(jk)          = 1.0_wp
    END DO

    t1=b(z_dolic)-a(z_dolic)*c(z_dolic-1)
    t1=(z_trac_rhs(z_dolic)-a(z_dolic)*z_trac_rhs(z_dolic-1))/t1

    z_trac_rhs(z_dolic)=t1

    ! Backward sweep
    DO jk=z_dolic-1,1,-1
      z_trac_rhs(jk) = z_trac_rhs(jk)-c(jk)*z_trac_rhs(jk+1)
    END DO

    ! Now update of fields
    DO jk=1,z_dolic
      diff_column(je,jk,jb) = z_trac_rhs(jk)
    END DO
  ELSE
    diff_column(je,:,jb) = 0.0_wp
  ENDIF
  ELSE
    diff_column(je,:,jb) = field_column(je,:,jb)
  ENDIF
  END DO
END DO
 DO jk=slev, elev
   WRITE(*,*)'MAX/MIN before vert veloc diffusion ',jk, &
     &        MAXVAL(field_column(:,jk,:)), MINVAL(field_column(:,jk,:))
 END DO


 DO jk=slev, elev
   WRITE(*,*)'MAX/MIN after vert veloc diffusion ',jk, &
     &        MAXVAL(diff_column(:,jk,:)), MINVAL(diff_column(:,jk,:))
 END DO
END subroutine veloc_diffusion_vert_impl
  !-------------------------------------------------------------------------
! !
! !
! !>
! !!
! !!
! !! @par Revision History
! !! Developed  by  Peter Korn, MPI-M (2010).
! !!
! SUBROUTINE tracer_diffusion_vert_old( p_patch,&
!                                 & h_c, trac_c,&
!                                 & top_bc_tracer, bot_bc_tracer,&
!                                 & laplacian_trac_vert_out)
! !
! !
! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! REAL(wp), INTENT(inout)           :: trac_c(:,:,:)
! REAL(wp), INTENT(in)              :: h_c(:,:)
! REAL(wp), INTENT(in)              :: top_bc_tracer(:,:) 
! REAL(wp), INTENT(in)              :: bot_bc_tracer(:,:) 
! REAL(wp), INTENT(out)             :: laplacian_trac_vert_out(:,:,:)
! !
! !Local variables
! INTEGER :: slev, elev
! INTEGER :: jc, jk, jb
! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! REAL(wp), PARAMETER :: K_T_v = 1.0e-2_wp !Preliminary solution should be taken from params
! ! REAL(wp) :: z_T_i(nproma,n_zlev+1,p_patch%nblks_c),  &
! !   &         z_T_m(nproma,n_zlev,p_patch%nblks_c)
! !-----------------------------------------------------------------------
! ! trac_c(:,1,:)=4.0_wp
! ! trac_c(:,2,:)=3.0_wp
! ! trac_c(:,3,:)=1.5_wp
! ! trac_c(:,4,:)=1.0_wp
! laplacian_trac_vert_out = 0.0_wp
! 
! i_startblk = p_patch%cells%start_blk(1,1)
! i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
! slev = 1
! elev = n_zlev
! 
! !1 Vertical derivative of tracer
! DO jk = slev, elev
!   DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
!                        i_startidx, i_endidx, 1,min_rlcell)
!     DO jc = i_startidx, i_endidx
! 
!       !1a) 0cean surface
!       IF(jk==slev)THEN
! 
!         IF (p_patch%patch_oce%dolic_c(jc,jb) /= 0) THEN
!           laplacian_trac_vert_out(jc,jk,jb) = &
!           &( top_bc_tracer(jc,jb)& !=0 zero flux
!           & -K_T_v*(trac_c(jc,jk,jb)-trac_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk+1))/&
!           & p_patch%patch_oce%del_zlev_m(jk)
! !           laplacian_trac_vert_out(jc,jk,jb) = &
! !           &(       (top_bc_tracer(jc,jb)-trac_c(jc,jk,jb))/(p_patch%patch_oce%del_zlev_i(jk)+h_c(jc,jb))&
! !           & -K_T_v*(trac_c(jc,jk,jb)-trac_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk+1))/p_patch%patch_oce%del_zlev_m(jk)
!         ENDIF
! 
!       !1b) ocean interior 
!       ELSEIF ( jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN
!         laplacian_trac_vert_out(jc,jk,jb)&
!         & = K_T_v *&
!         & ( (trac_c(jc,jk-1,jb)-trac_c(jc,jk,jb))/p_patch%patch_oce%del_zlev_i(jk-1)&
!         & - (trac_c(jc,jk,jb)-trac_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk))/&
!         & p_patch%patch_oce%del_zlev_m(jk)
! ! write(*,*)'T-diff interior:',jk,trac_c(jc,jk-1,jb),trac_c(jc,jk,jb),u_c(jc,jk,jb),trac_c(jc,jk+1,jb), z_T_i(jc,jk,jb),&
! ! &p_patch%patch_oce%dolic_c(jc,jb) 
!       !1c) ocean bottom 
!       ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
! 
!         IF(p_patch%patch_oce%n_zlev>=2)THEN
!           laplacian_trac_vert_out(jc,jk,jb)&
!           & = (K_T_v*(trac_c(jc,elev-1,jb)-trac_c(jc,elev,jb))/&
!           & p_patch%patch_oce%del_zlev_i(elev-1)&
!           & -bot_bc_tracer(jc,jb))/p_patch%patch_oce%del_zlev_m(elev)
! 
!  !write(*,*)'u-diff botom:',jk,trac_c(jc,elev-1,jb),trac_c(jc,elev,jb),trac_c(jc,elev,jb),bot_bc_u_c(jc,jb), z_adv_u_i(jc,jk,jb)
!         ENDIF
!       ENDIF
!     END DO
!   END DO
! END DO
! 
! 
! DO jk=slev, elev
!   WRITE(*,*)'MAX/MIN vert tracer diffusion ',jk, &
!     &        MAXVAL(laplacian_trac_vert_out(:,jk,:)), MINVAL(laplacian_trac_vert_out(:,jk,:))
! END DO
! 
! END subroutine tracer_diffusion_vert_old
!   !-------------------------------------------------------------------------
!
!
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
SUBROUTINE velocity_diffusion_vert_mim_old( p_patch, p_diag, p_aux,h_c,p_param, laplacian_vn_out)
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
TYPE(t_hydro_ocean_diag)          :: p_diag
TYPE(t_hydro_ocean_aux)           :: p_aux
REAL(wp), INTENT(in)              :: h_c(:,:) 
TYPE(t_ho_params), INTENT(in)     :: p_param
REAL(wp)                          :: laplacian_vn_out(:,:,:)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &                              z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:veloc diffusion vert mimetic')
!-----------------------------------------------------------------------
z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c)%x = 0.0_wp
z_adv_u_m(nproma,n_zlev,p_patch%nblks_c)%x   = 0.0_wp

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
      !check if we have at least two layers of water
      !  #slo# - 2011-04-01 - Is this really intended here
      !  maybe this condition should be fulfilled everywhere
      !  then it must be calculated in fill_vertical_domain
      !  this condition could then be omitted here
     IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN

      !1a) 0cean surface
      IF ( jk==slev ) THEN

!           z_adv_u_i(jc,jk,jb) = &
!           &(       (top_bc_u_c(jc,jb)-u_c(jc,jk,jb))/(p_patch%patch_oce%del_zlev_i(jk)+h_c(jc,jb))&
!           & -K_veloc_v*(u_c(jc,jk,jb)-u_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk+1))/p_patch%patch_oce%del_zlev_m(jk)
          z_adv_u_i(jc,jk,jb)%x = &
          &( p_aux%bc_top_veloc_cc(jc,jb)%x &
          & -p_param%A_veloc_v(jc,jk,jb)*(p_diag%p_vn(jc,jk,jb)%x-p_diag%p_vn(jc,jk+1,jb)%x) &
          & /p_patch%patch_oce%del_zlev_i(jk+1))&
          & /p_patch%patch_oce%del_zlev_m(jk)
!write(*,*)'MIM:vert-diff:top',jc,jk,jb, z_adv_u_i(jc,jk,jb)%x, p_aux%bc_top_veloc_cc(jc,jb)%x
!        ENDIF
      !1b) ocean bottom 
      ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN

!        IF(p_patch%patch_oce%n_zlev>=2)THEN
          z_adv_u_i(jc,jk,jb)%x&
          & = (p_param%A_veloc_v(jc,jk,jb)&
          &*(p_diag%p_vn(jc,elev-1,jb)%x-p_diag%p_vn(jc,elev,jb)%x) &
          & / p_patch%patch_oce%del_zlev_i(elev-1)&
          & - p_aux%bc_bot_veloc_cc(jc,jb)%x)/p_patch%patch_oce%del_zlev_m(elev)
          !write(*,*)'u-diff botom:',jk,u_c(jc,elev-1,jb),u_c(jc,elev,jb),u_c(jc,elev,jb),bot_bc_u_c(jc,jb), z_adv_u_i(jc,jk,jb)
!        ENDIF
      !1c) ocean interior 
      ELSEIF( jk>slev .AND. jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN

        z_adv_u_i(jc,jk,jb)%x&
        & = p_param%A_veloc_v(jc,jk,jb) *&
     & ( (p_diag%p_vn(jc,jk-1,jb)%x-p_diag%p_vn(jc,jk,jb)%x)  /p_patch%patch_oce%del_zlev_i(jk-1)&
     & - (p_diag%p_vn(jc,jk,jb)%x  -p_diag%p_vn(jc,jk+1,jb)%x)/p_patch%patch_oce%del_zlev_i(jk))/&
     & p_patch%patch_oce%del_zlev_m(jk)
! write(*,*)'u-diff interior:',jk,u_c(jc,jk-1,jb),u_c(jc,jk,jb),u_c(jc,jk,jb),u_c(jc,jk+1,jb), z_adv_u_i(jc,jk,jb),&
! &p_patch%patch_oce%dolic_c(jc,jb) 
      !1c) ocean bottom 

        ENDIF  ! jk-condition
      ENDIF    ! at least 2 vertical layers
    END DO
  END DO
END DO

! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)  
CALL map_cell2edges( p_patch, z_adv_u_i, laplacian_vn_out)

! DO jk=slev, elev
!   WRITE(*,*)'MAX/MIN vert diffusion ',jk, &
!     &        MAXVAL(laplacian_vn_out(:,jk,:)), MINVAL(laplacian_vn_out(:,jk,:))
! END DO

END subroutine velocity_diffusion_vert_mim_old
!-------------------------------------------------------------------------  
!
!
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
SUBROUTINE velocity_diffusion_vert_rbf_old( p_patch, u_c, v_c, h_c, top_bc_u_c, top_bc_v_c,&
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

REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_v_i(nproma,n_zlev+1,p_patch%nblks_c),  &
  &         z_adv_u_m(nproma,n_zlev,p_patch%nblks_c),  &
  &         z_adv_v_m(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diffusion:veloc diffusion vert rbf')
!-----------------------------------------------------------------------

! #slo# set local variable to zero due to nag -nan compiler-option
z_adv_u_i(:,:,:) = 0.0_wp
z_adv_v_i(:,:,:) = 0.0_wp
z_adv_u_m(:,:,:) = 0.0_wp
z_adv_v_m(:,:,:) = 0.0_wp

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
      !1a) 0cean surface
      !check if we have at least two layers of water
     !IF (p_patch%patch_oce%dolic_c(jc,jb) >= 2) THEN
      IF(jk==slev)THEN
          ! u-component
!           z_adv_u_i(jc,jk,jb) = &
!           &(       (top_bc_u_c(jc,jb)-u_c(jc,jk,jb))/(p_patch%patch_oce%del_zlev_i(jk)+h_c(jc,jb))&
!           & -K_veloc_v*(u_c(jc,jk,jb)-u_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk+1))/p_patch%patch_oce%del_zlev_m(jk)
          z_adv_u_i(jc,jk,jb) = &
          &(       top_bc_u_c(jc,jb) &
          & -p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,jk,jb)-u_c(jc,jk+1,jb)) &
          & /p_patch%patch_oce%del_zlev_i(jk))/p_patch%patch_oce%del_zlev_m(jk)
          ! v-component
!           z_adv_v_i(jc,jk,jb) = K_veloc_v*&
!           &( (top_bc_v_c(jc,jb)-v_c(jc,jk,jb))/(p_patch%patch_oce%del_zlev_i(jk)+h_c(jc,jb))&
!           & -(    v_c(jc,jk,jb)-v_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk+1))/p_patch%patch_oce%del_zlev_m(jk)
          z_adv_v_i(jc,jk,jb) = &
          &( top_bc_v_c(jc,jb)&
          & -p_param%A_veloc_v(jc,jk,jb)*(    v_c(jc,jk,jb)-v_c(jc,jk+1,jb))&
          &/p_patch%patch_oce%del_zlev_i(jk)) &
          & /p_patch%patch_oce%del_zlev_m(jk)
!write(*,*)'RBD:vert-diff:top',jc,jk,jb, z_adv_u_i(jc,jk,jb),top_bc_u_c(jc,jb)
      !1b) ocean bottom 
      ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN

        IF(p_patch%patch_oce%n_zlev>=2)THEN
          ! u-component
!           z_adv_u_i(jc,jk,jb)&
!           & = (K_veloc_v*(u_c(jc,elev-1,jb)-u_c(jc,elev,jb))/ p_patch%patch_oce%del_zlev_i(elev-1)&
!           & -(u_c(jc,elev,jb)-bot_bc_u_c(jc,jb))/ p_patch%patch_oce%del_zlev_i(elev))/p_patch%patch_oce%del_zlev_m(elev)
          z_adv_u_i(jc,jk,jb)&
        & = (p_param%A_veloc_v(jc,jk,jb)*(u_c(jc,elev-1,jb)-u_c(jc,elev,jb))/&
           & p_patch%patch_oce%del_zlev_i(elev-1)&
        & - bot_bc_u_c(jc,jb))/p_patch%patch_oce%del_zlev_m(elev)

!           z_adv_v_i(jc,jk,jb)&
!           & = (K_veloc_v*(v_c(jc,elev-1,jb)-v_c(jc,elev,jb))/ p_patch%patch_oce%del_zlev_i(elev-1)&
!           & -(v_c(jc,elev,jb)-bot_bc_v_c(jc,jb))/ p_patch%patch_oce%del_zlev_i(elev))/p_patch%patch_oce%del_zlev_m(elev)
          z_adv_v_i(jc,jk,jb)&
        & = (p_param%A_veloc_v(jc,jk,jb)*(v_c(jc,elev-1,jb)-v_c(jc,elev,jb))/&
        & p_patch%patch_oce%del_zlev_i(elev-1)&
        & -bot_bc_v_c(jc,jb))/p_patch%patch_oce%del_zlev_m(elev)

        !write(*,*)'u-diff botom:',jk,u_c(jc,elev-1,jb),u_c(jc,elev,jb),u_c(jc,elev,jb),bot_bc_u_c(jc,jb), z_adv_u_i(jc,jk,jb)
        ENDIF

      !1c) ocean interior 
      ELSEIF(jk>slev .AND. jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN
        ! u-component
        z_adv_u_i(jc,jk,jb)&
        & = p_param%A_veloc_v(jc,jk,jb) *&
        & ( (u_c(jc,jk-1,jb)-u_c(jc,jk,jb))/p_patch%patch_oce%del_zlev_i(jk-1)&
        & - (u_c(jc,jk,jb)-u_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk)) &
        & / p_patch%patch_oce%del_zlev_m(jk)
        ! v-component
        z_adv_v_i(jc,jk,jb)&
        & = p_param%A_veloc_v(jc,jk,jb)*&
        & ( (v_c(jc,jk-1,jb)-v_c(jc,jk,jb))/p_patch%patch_oce%del_zlev_i(jk-1)&
        & - (v_c(jc,jk,jb)-v_c(jc,jk+1,jb))/p_patch%patch_oce%del_zlev_i(jk)) &
        & / p_patch%patch_oce%del_zlev_m(jk)
! write(*,*)'u-diff interior:',jk,u_c(jc,jk-1,jb),u_c(jc,jk,jb),u_c(jc,jk,jb),u_c(jc,jk+1,jb), z_adv_u_i(jc,jk,jb),&
! &p_patch%patch_oce%dolic_c(jc,jb) 
      !1c) ocean bottom 

      ENDIF  ! jk-condition
    !ENDIF   ! at least 2 vertical layers
    END DO
  END DO
END DO

!  DO jk=slev, elev
!    WRITE(*,*)'MAX/MIN vert diff A:',jk, &
!      &        MAXVAL(z_adv_u_i(:,jk,:)), MINVAL(z_adv_u_i(:,jk,:)),&
!      &        MAXVAL(z_adv_v_i(:,jk,:)), MINVAL(z_adv_v_i(:,jk,:))
!  END DO

! Step 2: Map result of previous calculations from cell centers to edges (for all vertical layers)
CALL primal_map_c2e( p_patch,&
                   & z_adv_u_i,&
                   & z_adv_v_i,&
                   & laplacian_vn_out )

! DO jk=slev, elev
!   WRITE(*,*)'MAX/MIN vert diffusion ',jk, &
!     &        MAXVAL(laplacian_vn_out(:,jk,:)), MINVAL(laplacian_vn_out(:,jk,:))
! END DO

END subroutine velocity_diffusion_vert_rbf_old
!------------------------------------------------------------------------
END MODULE mo_oce_diffusion
