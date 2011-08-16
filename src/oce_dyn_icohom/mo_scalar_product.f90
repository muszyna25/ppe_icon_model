!>
!!               The module <i>mo_scalar_product</i>.
!!
!!               The module <i>mo_scalar_product</i>
!! implements discrete scalar products which depend on the grid geometry only
!! are used to formulate the primitive equations in weak form.
!! The coefficients used for projections from triangular edges to cell centers
!! and vice versa are contained in the ocean part of the model domain and
!! are calculated in <i>mo_ocean_topo</i>.
!!
!! @par Revision History
!! Initial version  by Peter Korn and Stephan Lorenz,  MPI-M, Hamburg, October 2010
!! Modification by Stephan Lorenz, MPI-M, (2010-11-02)
!! - initial primal_flip_flop divided into basic parts of primal_map_e2c and c2e
!! Modification by Stephan Lorenz, MPI-M, (2010-11-16)
!! - implementation as primal_map_e2c, primal_map_e2c_no_edge_height, primal_map_c2e
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_scalar_product
!-------------------------------------------------------------------------
!
USE mo_kind,               ONLY: wp
USE mo_parallel_config,    ONLY: nproma
USE mo_run_config,         ONLY: dtime
USE mo_impl_constants,     ONLY: sea_boundary,     &!max_char_length, &
!  &                             land, land_boundary, boundary, sea,  &
  &                              min_rlcell, min_rledge !, min_rlvert
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_model_domain,       ONLY: t_patch
USE mo_oce_state,          ONLY: t_hydro_ocean_diag, t_hydro_ocean_state, v_base
USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce !, ab_gam
USE mo_math_utilities,     ONLY: gvec2cvec, cvec2gvec, t_cartesian_coordinates,gc2cc
!USE mo_oce_index,          ONLY: ne_b, ne_i, nv_b, nv_i, form4ar, ldbg, c_k!, c_b, c_i
!USE mo_exception,          ONLY: message, message_text

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: calc_scalar_product_for_veloc
PUBLIC :: map_cell2edges
PUBLIC :: map_cell2edges_2D
PUBLIC :: map_edges2cell
PUBLIC :: primal_map_c2e
PUBLIC :: dual_flip_flop
PUBLIC :: map_edges2edges
PRIVATE :: map_cell2edges_upwind
PRIVATE :: map_cell2edges_upwind2

interface map_edges2cell

  module procedure map_edges2cell_with_height
  module procedure map_edges2cell_without_height

end interface


interface map_edges2edges

  module procedure map_edges2edges_with_height
  module procedure map_edges2edges_without_height

end interface

CONTAINS
 !-------------------------------------------------------------------------
  !
  !>
  !! Function implements discrete scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The functions uses coefficients that are calculated in sbr
  !! "init_scalar_product".
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-9)
  !!
  FUNCTION dual_flip_flop(p_patch, vn_old_e, vn_new_e, vort_v,p_vn, &
                        & h_e, opt_slev, opt_elev) RESULT(vn_out_e)

  TYPE(t_patch), INTENT(IN)    :: p_patch
  REAL(wp), INTENT(inout)      :: vn_old_e(:,:,:)
  REAL(wp), INTENT(inout)      :: vn_new_e(:,:,:)
  REAL(wp), INTENT(in)         :: vort_v(:,:,:)
  TYPE(t_cartesian_coordinates):: p_vn(nproma,n_zlev,p_patch%nblks_c)
  REAL(wp), INTENT(IN)      :: h_e(:,:)
  INTEGER, INTENT(IN), OPTIONAL ::  opt_slev  ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL ::  opt_elev  ! optional vertical end level
  REAL(wp)                  :: vn_out_e(SIZE(vn_old_e,1), SIZE(vn_old_e,2), SIZE(vn_old_e,3))

  !Local variables
  !INTEGER, PARAMETER :: no_vert_edges = 6
  INTEGER :: slev, elev
  !INTEGER :: rl_start_v, rl_end_v
  INTEGER :: rl_start_e, rl_end_e
  !INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: il_v1, il_v2, ib_v1, ib_v2, il_e, ib_e
  INTEGER :: je, jb, ie, jk, iie_v1, iie_v2!, jv
  INTEGER :: ic1, ib1, ic2,ib2

  REAL(wp) :: z_thick, z_weight!, z_tmp1, z_tmp2
  REAL(wp) :: vn_tmp!z_vn
  REAL(wp) :: z_mean
  REAL(wp) :: z_h_approx
  REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)
  TYPE(t_cartesian_coordinates) :: u_v1_cc!(SIZE(vort_v,1),SIZE(vort_v,3))
  TYPE(t_cartesian_coordinates) :: u_v2_cc!(SIZE(vort_v,1),SIZE(vort_v,3))
  INTEGER  :: i_v1_ctr, i_v2_ctr

  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:dual_flip_flop')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
   vn_out_e(:,:,:) = 0.0_wp

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

!  rl_start_v = 1
!  rl_end_v = min_rlvert
!  i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
!  i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)

 rl_start_e = 1
 rl_end_e = min_rledge
 i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
 i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

 z_vn_e = vn_old_e!ab_gam*vn_new_e + (1.0_wp-ab_gam)*vn_old_e

  LEVEL_LOOP: DO jk = slev, elev
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2)
    EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

      CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
                       & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)

      EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e

        IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

        !Get indices of two adjacent vertices
        il_v1 = p_patch%edges%vertex_idx(je,jb,1)
        ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
        il_v2 = p_patch%edges%vertex_idx(je,jb,2)
        ib_v2 = p_patch%edges%vertex_blk(je,jb,2)

        u_v1_cc%x  = 0.0_wp
        z_weight   = 0.0_wp
        i_v1_ctr = 0
        iie_v1   = 0
        z_h_approx = 0.0_wp
        DO ie=1, p_patch%verts%num_edges(il_v1,ib_v1) !no_vert_edges

          il_e = p_patch%verts%edge_idx(il_v1,ib_v1,ie)
          ib_e = p_patch%verts%edge_blk(il_v1,ib_v1,ie)

          IF(il_e==je.AND.ib_e==jb)THEN
            iie_v1 = ie
          ENDIF
          z_thick = v_base%del_zlev_m(jk)
           IF ( iswm_oce == 1 ) THEN
             z_thick = h_e(il_e,ib_e)
           ELSEIF( iswm_oce /= 1 ) THEN 
             IF (jk == 1 )THEN
               z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e) 
             ENDIF
           ENDIF 

           IF ( v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary ) THEN
             z_weight = z_weight + v_base%variable_dual_vol_norm(il_v1,ib_v1,ie)*z_thick

             z_h_approx = z_h_approx + z_thick
 
             u_v1_cc%x = u_v1_cc%x +                      &
             &           v_base%edge2vert_coeff_cc(il_v1,ib_v1,ie)%x * &
             &           z_vn_e(il_e,jk,ib_e)*z_thick
             i_v1_ctr=i_v1_ctr+1
           ENDIF
        END DO
        IF(z_weight/=0.0_wp)THEN
            z_h_approx = z_h_approx/i_v1_ctr
            u_v1_cc%x = u_v1_cc%x/(z_weight)!*z_h_approx)
        ENDIF
        !Multiply vector with global vorticity
        !u_v1_cc%x = u_v1_cc%x * vort_v(il_v1,jk,ib_v1)


        u_v2_cc%x   = 0.0_wp
        z_weight    = 0.0_wp
        i_v2_ctr    = 0
        iie_v2      = 0
        z_h_approx  = 0.0_wp
        DO ie=1, p_patch%verts%num_edges(il_v2,ib_v2) !no_vert_edges

          il_e = p_patch%verts%edge_idx(il_v2,ib_v2,ie)
          ib_e = p_patch%verts%edge_blk(il_v2,ib_v2,ie)

          IF(il_e==je.AND.ib_e==jb)THEN
            iie_v2 = ie
          ENDIF

          z_thick = v_base%del_zlev_m(jk)
           IF ( iswm_oce == 1 ) THEN
             z_thick = h_e(il_e,ib_e)
           ELSEIF( iswm_oce /= 1 ) THEN 
             IF (jk == 1 )THEN
               z_thick = v_base%del_zlev_m(jk) + h_e(il_e,ib_e) 
             ENDIF
           ENDIF 

           IF ( v_base%lsm_oce_e(il_e,jk,ib_e)  <= sea_boundary ) THEN
             z_weight = z_weight + v_base%variable_dual_vol_norm(il_v2,ib_v2,ie)*z_thick

             z_h_approx = z_h_approx + z_thick

             u_v2_cc%x = u_v2_cc%x +                      &
             &           v_base%edge2vert_coeff_cc(il_v2,ib_v2,ie)%x * &
             &           z_vn_e(il_e,jk,ib_e)*z_thick
             i_v2_ctr=i_v2_ctr+1
           ENDIF
        END DO
        IF(z_weight/=0.0_wp)THEN
            z_h_approx = z_h_approx/i_v2_ctr
            u_v2_cc%x = u_v2_cc%x/(z_weight)!*z_h_approx)
        ENDIF
        !Multiply vector with global vorticity
        !u_v2_cc%x = u_v2_cc%x * vort_v(il_v2,jk,ib_v2)

!  IF(i_v1_ctr/=6.or.i_v2_ctr/=6)then
!  write(*,*)'vertex 1:2 counter',i_v1_ctr, i_v2_ctr,&
! &p_patch%verts%num_edges(il_v1,ib_v1),&
! &p_patch%verts%num_edges(il_v2,ib_v2)
!  ENDIF
            vn_out_e(je,jk,jb) = &
           &- DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x)&
           &+ DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x)
 
  vn_out_e(je,jk,jb)  = vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
  &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))/(i_v1_ctr+i_v2_ctr))



!          IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
!          &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
! 
! 
! 
!     u_v1_cc%x=u_v1_cc%x!*(vort_v(il_v1,jk,ib_v1))!+p_patch%verts%f_v(il_v1,ib_v1))
!     u_v2_cc%x=u_v2_cc%x!*(vort_v(il_v2,jk,ib_v2))!+p_patch%verts%f_v(il_v2,ib_v2))
! 
!             vn_out_e(je,jk,jb) = &
!            &- DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x)&
!            &+ DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x)
!  
!   vn_out_e(je,jk,jb)  = vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
!   &+0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
! write(*,*)'comparison',je,jb,&
! & DOT_PRODUCT(u_v1_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x),&
! &  DOT_PRODUCT(u_v2_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x),&
! &- DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x),&
! &DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x),vn_out_e(je,jk,jb)
! ! IF(jk==1)THEN
! ! write(*,*)'dual-flip-flop: ',je,jb,vn_out_e(je,jk,jb),&
! ! &- 2.0_wp*DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x),&
! ! &  2.0_wp*DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x)
! ! write(*,*)'dual-flip-flop2:',je,jb,&
! ! &2.0_wp*DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_vector_cc(il_v1,ib_v1,iie_v1)%x),&
! ! &-2.0_wp*DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_vector_cc(il_v2,ib_v2,iie_v2)%x )
! ! ENDIF
! 
! 
!          ELSEIF( i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
!          &.AND.  i_v2_ctr< p_patch%verts%num_edges(il_v2,ib_v2))THEN
!    z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
!    &/REAL(i_v1_ctr+i_v2_ctr,wp)
!    u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
!    u_v2_cc%x=u_v2_cc%x*(z_mean+p_patch%verts%f_v(il_v2,ib_v2))
! 
!            vn_out_e(je,jk,jb) = &
!           & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x)
! 
!             vn_out_e(je,jk,jb) = DOT_PRODUCT(u_v1_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x)
! 
!           ELSEIF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
!           &.AND.  i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
! !z_mean = 0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v1,jk,ib_v1))
! !   z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
! !   &/REAL(i_v1_ctr+i_v2_ctr,wp)
! !   u_v1_cc%x=u_v1_cc%x*(z_mean+p_patch%verts%f_v(il_v1,ib_v1))
!     u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
! 
!  !           vn_out_e(je,jk,jb) = &
!  !          & -2.0_wp*DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x)
! 
!             vn_out_e(je,jk,jb) = DOT_PRODUCT(u_v2_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x)
! 
!           ELSEIF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
!           &.AND.  i_v2_ctr<p_patch%verts%num_edges(il_v2,ib_v2))THEN
! 
! ELSE
!    z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
!    &/REAL(i_v1_ctr+i_v2_ctr,wp)
!    u_v1_cc%x=u_v1_cc%x*(z_mean+p_patch%verts%f_v(il_v1,ib_v1))
!    u_v2_cc%x=u_v2_cc%x*(z_mean+p_patch%verts%f_v(il_v2,ib_v2))
!            vn_out_e(je,jk,jb) = &
!           &- DOT_PRODUCT(u_v2_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x)&
!           &+ DOT_PRODUCT(u_v1_cc%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x)
! 
!           ic1 = p_patch%edges%cell_idx(je,jb,1)
!           ib1 = p_patch%edges%cell_blk(je,jb,1)
!           ic2 = p_patch%edges%cell_idx(je,jb,2)
!           ib2 = p_patch%edges%cell_blk(je,jb,2)
! 
!           vn_out_e(je,jk,jb) = 0.5_wp*&
!           &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
!           &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
! 
!           vn_out_e(je,jk,jb) =&
!           & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
!           &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
!           &/(i_v1_ctr+i_v2_ctr)
! 
!           ENDIF
       ELSE
         vn_out_e(je,jk,jb)= 0.0_wp
       ENDIF
! write(*,*)'vertex vect 1:',il_v1,ib_v1,(jb-1)*nproma+je,jb,je, u_v_cc(il_v1,ib_v1)%x
! write(*,*)'vertex vect 2:',il_v2,ib_v2,(jb-1)*nproma+je,jb,je, u_v_cc(il_v2,ib_v2)%x

!write(*,*)'P^T coeff 1:',je,jb,1, v_base%edge2vert_coeff_t(je,jb,1)%x
!write(*,*)'P^T coeff 2:',je,jb,2, v_base%edge2vert_coeff_t(je,jb,2)%x
      END DO EDGE_IDX_LOOP
    END DO EDGE_BLK_LOOP
!$OMP END DO
!$OMP END PARALLEL
  END DO LEVEL_LOOP
!stop
  END FUNCTION dual_flip_flop
!  !-------------------------------------------------------------------------
!   !
!   !>
!   !! Function implements discrete scalar product on the primal grid. This
!   !! scalar product depends on the grid geometry only and  is used to formulate the primitive
!   !! equations in weak form. The functions uses coefficients that are calculated in sbr
!   !! "init_scalar_product".
!   !!
!   !! @par Revision History
!   !!  developed by Peter Korn, MPI-M (2010-9)
!   !!
!   FUNCTION dual_flip_flop(p_patch, vn_old_e, vn_new_e, vort_v,&
!                         & h_e, opt_slev, opt_elev) RESULT(vn_out_e)
! 
!   TYPE(t_patch), INTENT(IN) :: p_patch
!   REAL(wp), INTENT(inout)      :: vn_old_e(:,:,:)
!   REAL(wp), INTENT(inout)      :: vn_new_e(:,:,:)
!   REAL(wp), INTENT(in)      :: vort_v(:,:,:)
!   REAL(wp), INTENT(IN)      :: h_e(:,:)
!   INTEGER, INTENT(IN), OPTIONAL ::  opt_slev  ! optional vertical start level
!   INTEGER, INTENT(IN), OPTIONAL ::  opt_elev  ! optional vertical end level
!   REAL(wp)                  :: vn_out_e(SIZE(vn_old_e,1), SIZE(vn_old_e,2), SIZE(vn_old_e,3))
! 
!   !Local variables
!   INTEGER, PARAMETER :: no_vert_edges = 6
!   INTEGER :: slev, elev
!   INTEGER :: rl_start_v, rl_end_v
!   INTEGER :: rl_start_e, rl_end_e
!   INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
!   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
!   INTEGER :: il_v1, il_v2, ib_v1, ib_v2, il_e, ib_e
!   INTEGER :: je, jb, jv, ie, jk
! 
!   REAL(wp) :: z_thick, z_weight
!   REAL(wp) :: z_vn
!   REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)
!   TYPE(t_cartesian_coordinates) :: u_v_cc(SIZE(vort_v,1),SIZE(vort_v,3))
!   !TYPE(t_cartesian_coordinates) ::  u_tmp(SIZE(vort_v,1),SIZE(vort_v,3))
!   REAL(wp) :: u_v(SIZE(vort_v,1),SIZE(vort_v,3))
!   REAL(wp) :: v_v(SIZE(vort_v,1),SIZE(vort_v,3))
!   INTEGER  :: i_ctr
! 
!   !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!   !  & routine = ('mo_scalar_product:dual_flip_flop')
!   !-----------------------------------------------------------------------
!   !CALL message (TRIM(routine), 'start')
!    vn_out_e(:,:,:) = 0.0_wp
! 
!   ! check optional arguments
!   IF ( PRESENT(opt_slev) ) THEN
!     slev = opt_slev
!   ELSE
!     slev = 1
!   END IF
!   IF ( PRESENT(opt_elev) ) THEN
!     elev = opt_elev
!   ELSE
!     elev = n_zlev
!   END IF
! 
!  rl_start_v = 1
!  rl_end_v = min_rlvert
!  i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
!  i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)
! 
!  rl_start_e = 1
!  rl_end_e = min_rledge
!  i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
!  i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
! 
!  z_vn_e = ab_gam*vn_new_e + (1.0_wp-ab_gam)*vn_old_e
! 
!   LEVEL_LOOP: DO jk = slev, elev
! 
!   u_v(:,:) = 0.0_wp
!   v_v(:,:) = 0.0_wp
! 
!  !Step 1: edge to vertex mapping
! 
! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jv,ie,i_startidx_v,i_endidx_v,il_e,ib_e)
!     VERT_BLK_LOOP: DO jb = i_startblk_v, i_endblk_v
! 
!       CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v, &
!         &                rl_start_v, rl_end_v)
! 
!       VERT_IDX_LOOP: DO jv = i_startidx_v, i_endidx_v
!         u_v(jv,jb) = 0.0_wp
!         v_v(jv,jb) = 0.0_wp
!         u_v_cc(jv,jb)%x= 0.0_wp
!         z_weight   = 0.0_wp
!         i_ctr = 0
!         DO ie=1, p_patch%verts%num_edges(jv,jb) !no_vert_edges
! 
!           il_e = p_patch%verts%edge_idx(jv,jb,ie)
!           ib_e = p_patch%verts%edge_blk(jv,jb,ie)
! 
!           z_thick = v_base%del_zlev_m(jk)
!            IF ( iswm_oce == 1 ) THEN
!              z_thick = h_e(il_e,ib_e)
!            ELSEIF( iswm_oce /= 1 ) THEN 
!              IF (jk == 1 )THEN
!                z_thick = v_base%del_zlev_m(jk) + h_e(il_e,ib_e) 
!              ENDIF
!            ENDIF 
! 
!            IF ( v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary ) THEN
!              z_weight = z_weight + v_base%variable_dual_vol_norm(jv,jb,ie)*z_thick
!              u_v_cc(jv,jb)%x = u_v_cc(jv,jb)%x +                                  &
!              &                 v_base%edge2vert_coeff_cc(jv,jb,ie)%x * &
!              &                       z_vn_e(il_e,jk,ib_e)*z_thick
!              i_ctr=i_ctr+1
!            ENDIF
! ! write(*,*)'dual coeff',jv,jb,ie,&
! ! & v_base%edge2vert_coeff(jv,jb,ie,1),&
! ! & v_base%edge2vert_coeff(jv,jb,ie,2), z_thick
! 
!         END DO
!         IF(z_weight/=0.0_wp)THEN
!           IF(i_ctr<=1)THEN
!             u_v_cc(jv,jb)%x=0.0_wp
!           ELSE
!             u_v_cc(jv,jb)%x = u_v_cc(jv,jb)%x/z_weight
!           ENDIF
!         ELSE
!           u_v_cc(jv,jb)%x=0.0_wp
!         ENDIF
! 
!         !Multiply vector with global vorticity
!         u_v_cc(jv,jb)%x = u_v_cc(jv,jb)%x * vort_v(jv,jk,jb)
! ! IF(z_weight/=0.0_wp)THEN
! ! write(*,*)'compare',jv,jb,u_v_cc(jv,jb)%x
! ! write(*,*)'compare',jv,jb,u_tmp(jv,jb)%x 
! !ENDIF
!       END DO VERT_IDX_LOOP
!     END DO VERT_BLK_LOOP
! !$OMP END DO
! 
! 
! !$OMP DO PRIVATE(jb,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2)
!     EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
! 
!       CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
!                        & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
! 
!       EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
!         !Get indices of two adjacent vertices
!         il_v1 = p_patch%edges%vertex_idx(je,jb,1)
!         ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
!         il_v2 = p_patch%edges%vertex_idx(je,jb,2)
!         ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
! 
! !        IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
!           vn_out_e(je,jk,jb) = &
!    &- DOT_PRODUCT(u_v_cc(il_v2,ib_v2)%x,v_base%edge2vert_coeff_cc_t(je,jb,2)%x)&
!    &+ DOT_PRODUCT(u_v_cc(il_v1,ib_v1)%x,v_base%edge2vert_coeff_cc_t(je,jb,1)%x)
! !        ELSE
! !          vn_out_e(je,jk,jb)= 0.0_wp
! !        ENDIF
! ! write(*,*)'vertex vect 1:',il_v1,ib_v1,(jb-1)*nproma+je,jb,je, u_v_cc(il_v1,ib_v1)%x
! ! write(*,*)'vertex vect 2:',il_v2,ib_v2,(jb-1)*nproma+je,jb,je, u_v_cc(il_v2,ib_v2)%x
! 
! !write(*,*)'P^T coeff 1:',je,jb,1, v_base%edge2vert_coeff_t(je,jb,1)%x
! !write(*,*)'P^T coeff 2:',je,jb,2, v_base%edge2vert_coeff_t(je,jb,2)%x
!       END DO EDGE_IDX_LOOP
!     END DO EDGE_BLK_LOOP
! !$OMP END DO
! !$OMP END PARALLEL
!   END DO LEVEL_LOOP
!   END FUNCTION dual_flip_flop
!-------------------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE calc_scalar_product_for_veloc( p_patch, vn_e_old, vn_e_new,&
                                          & h_e, p_diag)

  TYPE(t_patch), INTENT(IN) :: p_patch            ! patch on which computation is performed
  REAL(wp), INTENT(IN)      :: vn_e_old(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
  REAL(wp), INTENT(IN)      :: vn_e_new(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
  REAL(wp), INTENT(IN)      :: h_e(:,:)           ! SW-case: h_e is thicknerss at edges ! 3D case: h_e is surface elevation at edges 
  TYPE(t_hydro_ocean_diag)  :: p_diag

  !Local variables
  !INTEGER, PARAMETER :: no_cell_edges = 3
  INTEGER :: slev, elev
  INTEGER :: rl_start_c, rl_end_c!, rl_start_e, rl_end_e 
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  !INTEGER :: il_e, ib_e
  !INTEGER :: il_e1, ib_e1,il_e2, ib_e2, il_e3, ib_e3
  INTEGER :: jc, jb, jk!, ie,je
  !REAL(wp) :: z_kin(nproma,n_zlev,p_patch%nblks_c) 
  !REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)    ! input vector (nproma,n_zlev,nblks_e)
  !REAL(wp) :: z_weight, z_tmp
   TYPE(t_cartesian_coordinates)    :: z_pv_cc(nproma,n_zlev,p_patch%nblks_c)
!   REAL(wp) :: vort_flux_e(nproma,1,p_patch%nblks_e)
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

slev = 1
elev = n_zlev

!rl_start_e = 1
!rl_end_e  = min_rledge

rl_start_c = 1
rl_end_c  = min_rlcell

i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

!i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
!i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

!z_vn_e = ab_gam*vn_e_new + (1.0_wp-ab_gam)*vn_e_old

!Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
CALL map_edges2cell( p_patch, vn_e_old, p_diag%p_vn, h_e )
!CALL map_edges2cell( p_patch, z_vn_e, p_diag%p_vn, h_e )

! !   !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
!   CELL_BLK_LOOP_SWM: DO jb = i_startblk_c, i_endblk_c
!     CALL get_indices_c( p_patch, jb,&
!                       & i_startblk_c, i_endblk_c,&
!                       & i_startidx_c, i_endidx_c,&
!                       & rl_start_c, rl_end_c)
!     LEVEL_LOOP_SWM: DO jk = slev, elev
!       CELL_IDX_LOOP_SWM: DO jc =  i_startidx_c, i_endidx_c
!         !calculate velocity reconstruction at cell center
!         z_weight           = 0.0_wp
!         z_pv_cc(jc,jk,jb)%x = 0.0_wp
!         DO ie=1, no_cell_edges
! 
!           il_e = p_patch%cells%edge_idx(jc,jb,ie)
!           ib_e = p_patch%cells%edge_blk(jc,jb,ie)
! 
!           z_weight = z_weight&
!           & + v_base%variable_vol_norm(jc,jb,ie)*p_diag%thick_e(il_e,ib_e)
!           z_pv_cc(jc,jk,jb)%x = z_pv_cc(jc,jk,jb)%x&
!                            & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
!                            & * z_vn_e(il_e,jk,ib_e)*p_diag%thick_e(il_e,ib_e)
!         END DO
! !        IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
!           z_pv_cc(jc,jk,jb)%x = z_pv_cc(jc,jk,jb)%x / z_weight
! !        ELSE
! !          z_pv_cc(jc,jk,jb)%x=0.0_wp
! !        ENDIF
!       END DO CELL_IDX_LOOP_SWM
!     END DO LEVEL_LOOP_SWM
!   END DO CELL_BLK_LOOP_SWM

CALL map_edges2cell( p_patch, vn_e_old, z_pv_cc)
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb,&
                    & i_startblk_c, i_endblk_c,&
                    & i_startidx_c, i_endidx_c,&
                    & rl_start_c, rl_end_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
  !calculate kinetic energy

  DO jk = slev, elev
    DO jc =  i_startidx_c, i_endidx_c

!      IF ( v_base%lsm_oce_c(jc,jk,jb) > sea_boundary ) THEN
!         p_diag%kin(jc,jk,jb) = 0.0_wp
!       ELSE
        ! p_diag%kin(jc,jk,jb) = 0.5_wp*DOT_PRODUCT(p_diag%p_vn(jc,jk,jb)%x,p_diag%p_vn(jc,jk,jb)%x)
      p_diag%kin(jc,jk,jb) = 0.5_wp*DOT_PRODUCT(z_pv_cc(jc,jk,jb)%x,z_pv_cc(jc,jk,jb)%x)
!       write(*,*)'kin energy',jc,jk,jb,p_diag%kin(jc,jk,jb), p_diag%p_vn(jc,jk,jb)%x!z_kin(jc,jk,jb)
!       ENDIF
    END DO
  END DO
END DO

! DO jk = slev, elev
!   write(*,*)'max/min kin energy:',maxval(p_diag%kin(:,jk,:)), minval(p_diag%kin(:,jk,:))!,&
! !&maxval(z_kin(:,1,:)), minval(z_kin(:,1,:))
! END DO
!convert cartesian velocity vector p_diag%p_vn(jc,jk,jb)%x to geographical coordinate system
!for output
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb,&
                      & i_startblk_c, i_endblk_c,&
                      & i_startidx_c, i_endidx_c,&
                      & rl_start_c, rl_end_c)
  DO jk = slev, elev
    DO jc =  i_startidx_c, i_endidx_c
      CALL cvec2gvec ( p_diag%p_vn(jc,jk,jb)%x(1),     &
                     & p_diag%p_vn(jc,jk,jb)%x(2),     &
                     & p_diag%p_vn(jc,jk,jb)%x(3),     &
                     & p_patch%cells%center(jc,jb)%lon,&
                     & p_patch%cells%center(jc,jb)%lat,&
                     & p_diag%u(jc,jk,jb), p_diag%v(jc,jk,jb) )
! DO ie=1,3
! il_e = p_patch%cells%edge_idx(jc,jb,ie)
! ib_e = p_patch%cells%edge_blk(jc,jb,ie)
! write(*,*)'vn',jc,jb,ie,&
! &p_diag%u(jc,jk,jb)*p_patch%edges%primal_normal(il_e,ib_e)%v1&
! &+p_diag%v(jc,jk,jb)*p_patch%edges%primal_normal(il_e,ib_e)%v2,&
! &DOT_PRODUCT(p_diag%p_vn(jc,jk,jb)%x,p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
! END DO
    END DO
  END DO
END DO

CALL map_cell2edges( p_patch, p_diag%p_vn, p_diag%ptp_vn)

!   EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
!     CALL get_indices_e(p_patch, jb,&
!                      & i_startblk_e, i_endblk_e,&
!                      & i_startidx_e, i_endidx_e,&
!                      & rl_start_e, rl_end_e)
!     LEVEL_LOOP_E: DO jk = slev, elev
!       EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
! 
!        !Store mapped value
!        z_tmp = p_diag%ptp_vn(je,jk,jb)
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
!       IF(    v_base%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary&
!        &.AND.v_base%lsm_oce_c(il_c2,jk,ib_c2) == sea)THEN
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
!         p_diag%ptp_vn(il_e1,jk,ib_e1) = z_vn_e(il_e1,jk,ib_e1)
!         p_diag%ptp_vn(il_e2,jk,ib_e2) = z_vn_e(il_e2,jk,ib_e2)
!         p_diag%ptp_vn(il_e3,jk,ib_e3) = z_vn_e(il_e3,jk,ib_e3)
! 
!        p_diag%ptp_vn(je,jk,jb)=z_tmp
! 
!       ELSEIF(    v_base%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary&
!        &.AND.v_base%lsm_oce_c(il_c2,jk,ib_c2) == sea)THEN
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
!         p_diag%ptp_vn(il_e1,jk,ib_e1) = z_vn_e(il_e1,jk,ib_e1)
!         p_diag%ptp_vn(il_e2,jk,ib_e2) = z_vn_e(il_e2,jk,ib_e2)
!         p_diag%ptp_vn(il_e3,jk,ib_e3) = z_vn_e(il_e3,jk,ib_e3)
! 
!        p_diag%ptp_vn(je,jk,jb)=z_tmp
! 
!       ELSEIF(    v_base%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary&
!        &.AND.v_base%lsm_oce_c(il_c2,jk,ib_c2) == sea_boundary)THEN
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
!         p_diag%ptp_vn(il_e1,jk,ib_e1) = z_vn_e(il_e1,jk,ib_e1)
!         p_diag%ptp_vn(il_e2,jk,ib_e2) = z_vn_e(il_e2,jk,ib_e2)
!         p_diag%ptp_vn(il_e3,jk,ib_e3) = z_vn_e(il_e3,jk,ib_e3)
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
!         p_diag%ptp_vn(il_e1,jk,ib_e1) = z_vn_e(il_e1,jk,ib_e1)
!         p_diag%ptp_vn(il_e2,jk,ib_e2) = z_vn_e(il_e2,jk,ib_e2)
!         p_diag%ptp_vn(il_e3,jk,ib_e3) = z_vn_e(il_e3,jk,ib_e3)
!      ENDIF
!      END DO EDGE_IDX_LOOP
!    END DO LEVEL_LOOP_E
!  END DO EDGE_BLK_LOOP
! DO jk = slev, elev
!   write(*,*)'max/min u,v:',maxval(p_diag%u(:,jk,:)), minval(p_diag%u(:,jk,:)),&
!   &maxval(p_diag%v(:,jk,:)), minval(p_diag%v(:,jk,:))
! END DO



!  EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
!    CALL get_indices_e(p_patch, jb,&
!                     & i_startblk_e, i_endblk_e,&
!                     & i_startidx_e, i_endidx_e,&
!                     & rl_start_e, rl_end_e)
!    LEVEL_LOOP_E: DO jk = slev, elev
!      EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
 !   write(*,*)'vn:ptp_vn:',je,jk,jb,z_vn_e(je,jk,jb), p_diag%ptp_vn(je,jk,jb),&
!    !&(z_vn_e(je,jk,jb)- p_diag%ptp_vn(je,jk,jb))
! ! & DOT_PRODUCT(p_diag%p_vn(il_c1,jk,ib_c1)%x, v_base%edge2cell_coeff_t(je,jb,1)%x),&
! ! &DOT_PRODUCT(p_diag%p_vn(il_c2,jk,ib_c2)%x, v_base%edge2cell_coeff_t(je,jb,2)%x)
!     END DO EDGE_IDX_LOOP
!   END DO LEVEL_LOOP_E
! END DO EDGE_BLK_LOOP
!stop
!write(90,*)'-----------------------'
! DO jk=1,n_zlev
! write(*,*)'NEW SCALAR: L infty difference:',jk,&
! &maxval(vn_e(:,jk,:)-p_diag%ptp_vn(:,jk,:))&
! &,minval(vn_e(:,jk,:)-p_diag%ptp_vn(:,jk,:)) 
! write(*,*)'NEW SCALAR: L2 difference    :',jk,&
! &SUM( (vn_e(:,jk,:)-p_diag%ptp_vn(:,jk,:))&
! &*(vn_e(:,jk,:)-p_diag%ptp_vn(:,jk,:)))
! write(*,*)'NEW SCALAR: norm',jk,&
! &SUM( vn_e(:,jk,:)*vn_e(:,jk,:)),&
! &SUM( p_diag%ptp_vn(:,jk,:)*p_diag%ptp_vn(:,jk,:))
! END DO
END SUBROUTINE calc_scalar_product_for_veloc
!-------------------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_edges2edges_with_height( p_patch, vn_in, vn_out,&
                                        & h_e,opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN) :: p_patch            ! patch on which computation is performed
  REAL(wp), INTENT(INOUT)   :: vn_in(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
  REAL(wp), INTENT(INOUT)   :: vn_out(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
  REAL(wp), INTENT(IN)      :: h_e(:,:)           ! SW-case: h_e is thicknerss at edges ! 3D case: h_e is surface elevation at edges 

  INTEGER, INTENT(IN), OPTIONAL ::  opt_slev  ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL ::  opt_elev  ! optional vertical end level

  !Local variables
!   INTEGER, PARAMETER :: no_cell_edges = 3
!   INTEGER :: slev, elev
!   INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e 
!   INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
!   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
!   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
!   INTEGER :: il_e, ib_e
!   INTEGER :: jc, jb, jk, ie,je
   TYPE(t_cartesian_coordinates)    :: z_p_vn_cc(nproma,n_zlev,p_patch%nblks_c)
!   REAL(wp) :: vort_flux_e(nproma,1,p_patch%nblks_e)
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')


  CALL map_edges2cell( p_patch, vn_in, z_p_vn_cc, h_e, opt_slev, opt_elev )
  CALL map_cell2edges( p_patch, z_p_vn_cc, vn_out, opt_slev, opt_elev )
! rl_start_e = 1
! rl_end_e  = min_rledge
! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
! 
! EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
!   CALL get_indices_e(p_patch, jb,&
!                    & i_startblk_e, i_endblk_e,&
!                    & i_startidx_e, i_endidx_e,&
!                    & rl_start_e, rl_end_e)
!   LEVEL_LOOP_E: DO jk = slev, elev
!     EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
!   write(*,*)'vn:ptp_vn:',je,jk,jb,vn_in(je,jk,jb), vn_out(je,jk,jb)!,&
!   &(vn_in(je,jk,jb)- vn_out(je,jk,jb))
!    END DO EDGE_IDX_LOOP
!  END DO LEVEL_LOOP_E
!END DO EDGE_BLK_LOOP
!stop
END SUBROUTINE map_edges2edges_with_height
 !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_edges2edges_without_height( p_patch, vn_in, vn_out,&
                                           & opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN) :: p_patch            ! patch on which computation is performed
  REAL(wp), INTENT(INOUT)   :: vn_in(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
  REAL(wp), INTENT(INOUT)   :: vn_out(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)

  INTEGER, INTENT(IN), OPTIONAL ::  opt_slev  ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL ::  opt_elev  ! optional vertical end level

  !Local variables
!   INTEGER, PARAMETER :: no_cell_edges = 3
!   INTEGER :: slev, elev
!   INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e 
!   INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
!   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
!   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
!   INTEGER :: il_e, ib_e
!   INTEGER :: jc, jb, jk, ie,je
   TYPE(t_cartesian_coordinates)    :: z_p_vn_cc(nproma,n_zlev,p_patch%nblks_c)
!   REAL(wp) :: vort_flux_e(nproma,1,p_patch%nblks_e)
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')


  CALL map_edges2cell( p_patch, vn_in, z_p_vn_cc, opt_slev, opt_elev)
  CALL map_cell2edges( p_patch, z_p_vn_cc, vn_out, opt_slev, opt_elev )
! rl_start_e = 1
! rl_end_e  = min_rledge
! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
! 
! EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
!   CALL get_indices_e(p_patch, jb,&
!                    & i_startblk_e, i_endblk_e,&
!                    & i_startidx_e, i_endidx_e,&
!                    & rl_start_e, rl_end_e)
!   LEVEL_LOOP_E: DO jk = slev, elev
!     EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
!   write(*,*)'vn:ptp_vn:',je,jk,jb,vn_in(je,jk,jb), vn_out(je,jk,jb)!,&
!   &(vn_in(je,jk,jb)- vn_out(je,jk,jb))
!    END DO EDGE_IDX_LOOP
!  END DO LEVEL_LOOP_E
!END DO EDGE_BLK_LOOP
!stop
END SUBROUTINE map_edges2edges_without_height
 !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_cell2edges( p_patch, p_vn_c, ptp_vn, opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN)                 :: p_patch          ! patch on which computation is performed
  TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vn_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(INOUT)                   :: ptp_vn(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)  
  INTEGER, INTENT(IN), OPTIONAL             :: opt_slev        ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL             :: opt_elev        ! optional vertical end level


  !Local variables
  INTEGER :: slev, elev
  INTEGER :: rl_start_e, rl_end_e 
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: je, jb, jk,i! ie,je
  INTEGER :: i_ctr_c1, i_ctr_c2, il_c1,ib_c1, il_c2,ib_c2, ile_c1 , ibe_c1, ile_c2 , ibe_c2
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

rl_start_e = 1
rl_end_e  = min_rledge

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)


! calculation of transposed P^TPv from Pv (incart coord)
EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb,&
                   & i_startblk_e, i_endblk_e,&
                   & i_startidx_e, i_endidx_e,&
                   & rl_start_e, rl_end_e)

  LEVEL_LOOP_E: DO jk = slev, elev
    EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e

      IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        i_ctr_c1 = 0
        i_ctr_c2 = 0
        DO i=1,3
          ile_c1 = p_patch%cells%edge_idx(il_c1,ib_c1,i)
          ibe_c1 = p_patch%cells%edge_blk(il_c1,ib_c1,i)

          ile_c2 = p_patch%cells%edge_idx(il_c2,ib_c2,i)
          ibe_c2 = p_patch%cells%edge_blk(il_c2,ib_c2,i)

          IF ( v_base%lsm_oce_e(ile_c1,jk,ibe_c1) <= sea_boundary ) THEN
            i_ctr_c1 = i_ctr_c1 +1
          ENDIF
          IF ( v_base%lsm_oce_e(ile_c2,jk,ibe_c2) <= sea_boundary ) THEN
            i_ctr_c2 = i_ctr_c2 +1
          ENDIF
        END DO
        IF( i_ctr_c1 == 3 .AND. i_ctr_c2==3)THEN
          ptp_vn(je,jk,jb) =&
        & DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
        &+DOT_PRODUCT(p_vn_c(il_c2,jk,ib_c2)%x,v_base%edge2cell_coeff_cc_t(je,jb,2)%x)
         ELSE!IF(i_ctr_c1 > i_ctr_c2)THEN
           ptp_vn(je,jk,jb) =&
           &i_ctr_c1*(&
 &DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
 &+i_ctr_c2&
 &*DOT_PRODUCT(p_vn_c(il_c2,jk,ib_c2)%x,v_base%edge2cell_coeff_cc_t(je,jb,2)%x))&
           &/(i_ctr_c1+i_ctr_c2)
 
! !        ELSEIF(i_ctr_c1 < i_ctr_c2)THEN
! 
         ENDIF

!write(*,*)'counter', i_ctr_c1,i_ctr_c2
       ELSE
         ptp_vn(je,jk,jb) = 0.0_wp
       ENDIF
! IF(jb==3)then
! write(*,*)'vn, ptp_vn:',je,jk,jb,&
! & il_c1,ib_c1, p_vn_c(il_c1,jk,ib_c1)%x,&
! & il_c1,ib_c1, p_vn_c(il_c2,jk,ib_c2)%x
! ENDIF
    END DO EDGE_IDX_LOOP
  END DO LEVEL_LOOP_E
END DO EDGE_BLK_LOOP
!stop
  END SUBROUTINE map_cell2edges
!-----------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_cell2edges_upwind( p_patch, p_vec_c, p_vn_e, p_os, opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN)                 :: p_patch          ! patch on which computation is performed
  TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vec_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(OUT)                     :: p_vn_e(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
  TYPE(t_hydro_ocean_state)                 :: p_os
  INTEGER, INTENT(IN), OPTIONAL             :: opt_slev        ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL             :: opt_elev        ! optional vertical end level


  !Local variables
  INTEGER :: slev, elev
  INTEGER :: rl_start_e, rl_end_e 
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: je, jb, jk!, ie,je
  REAL(wp) :: z_dot_c1, z_dot_c2, norm_c1_c2, norm
  TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2
  TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

rl_start_e = 1
rl_end_e  = min_rledge

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)


! calculation of transposed P^TPv from Pv (incart coord)
EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb,&
                   & i_startblk_e, i_endblk_e,&
                   & i_startidx_e, i_endidx_e,&
                   & rl_start_e, rl_end_e)

  LEVEL_LOOP_E: DO jk = slev, elev
    EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e

      IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        cc_e0 = gc2cc(p_patch%edges%center(je,jb))
        cc_c1 = gc2cc(p_patch%cells%center(il_c1,ib_c1))
        cc_c2 = gc2cc(p_patch%cells%center(il_c2,ib_c2))

        cv_c1_c2%x = cc_c1%x - cc_c2%x
        cv_c1_e0%x = cc_e0%x - cc_c1%x
        cv_c2_e0%x = cc_e0%x - cc_c2%x


        norm_c1_c2 = SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))+SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x)) 
 
        IF(p_os%p_diag%ptp_vn(je,jk,jb)>=0.0_wp)THEN !velocity is pointing
          norm=SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))
          p_vn_e(je,jk,jb) =&
          & DOT_PRODUCT(p_vec_c(il_c1,jk,ib_c1)%x,&
          &(v_base%edge2cell_coeff_cc_t(je,jb,1)%x)*(norm_c1_c2/norm))!&
          !&+&
          !& DOT_PRODUCT(p_vec_c(il_c2,jk,ib_c2)%x,&
          !&(v_base%edge2cell_coeff_cc_t(je,jb,2)%x)*2.0_wp)!(norm_c1_c2/norm))
        !&+DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)

        ELSEIF(p_os%p_diag%ptp_vn(je,jk,jb)<0.0_wp)THEN
          norm=SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x))
          p_vn_e(je,jk,jb) =&
          & DOT_PRODUCT(p_vec_c(il_c2,jk,ib_c2)%x,&
          & (v_base%edge2cell_coeff_cc_t(je,jb,2)%x)*(norm_c1_c2/norm))!&
          !&+&
          !& DOT_PRODUCT(p_vec_c(il_c1,jk,ib_c1)%x,&
          !&(v_base%edge2cell_coeff_cc_t(je,jb,1)%x)*2.0_wp)!(norm_c1_c2/norm))
          !write(*,*)'frac',(norm_c1_c2/norm)
        ENDIF

       ELSE
         p_vn_e(je,jk,jb) = 0.0_wp
       ENDIF
! IF(jb==3)then
! write(*,*)'vn, ptp_vn:',je,jk,jb,&
! & il_c1,ib_c1, p_vn_c(il_c1,jk,ib_c1)%x,&
! & il_c1,ib_c1, p_vn_c(il_c2,jk,ib_c2)%x
! ENDIF
    END DO EDGE_IDX_LOOP
  END DO LEVEL_LOOP_E
END DO EDGE_BLK_LOOP

  END SUBROUTINE map_cell2edges_upwind
!-----------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_cell2edges_upwind2( p_patch, p_vec_c, p_vn_e, p_os, opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN)                 :: p_patch          ! patch on which computation is performed
  TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vec_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(OUT)                     :: p_vn_e(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
  TYPE(t_hydro_ocean_state)                 :: p_os
  INTEGER, INTENT(IN), OPTIONAL             :: opt_slev        ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL             :: opt_elev        ! optional vertical end level


  !Local variables
  INTEGER :: slev, elev
  INTEGER :: rl_start_e, rl_end_e 
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: je, jb, jk!, ie,je
  REAL(wp) :: z_dot_c1, z_dot_c2, norm_c1_c2, norm
  TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2
  TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

rl_start_e = 1
rl_end_e  = min_rledge

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)


! calculation of transposed P^TPv from Pv (incart coord)
EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb,&
                   & i_startblk_e, i_endblk_e,&
                   & i_startidx_e, i_endidx_e,&
                   & rl_start_e, rl_end_e)

  LEVEL_LOOP_E: DO jk = slev, elev
    EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e

      IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

 
          p_vn_e(je,jk,jb) =&
          & DOT_PRODUCT(p_vec_c(il_c1,jk,ib_c1)%x,&
          &(v_base%edge2cell_coeff_cc_t(je,jb,1)%x))&
          &+&
          & DOT_PRODUCT(p_vec_c(il_c2,jk,ib_c2)%x,&
          &(v_base%edge2cell_coeff_cc_t(je,jb,2)%x))&
          &-&
          p_os%p_diag%ptp_vn(je,jk,jb)*dtime*&
          &SQRT(SUM( (p_vec_c(il_c2,jk,ib_c2)%x-p_vec_c(il_c1,jk,ib_c1)%x)&
          &*p_patch%edges%inv_dual_edge_length(je,jb)&
          &*(p_vec_c(il_c2,jk,ib_c2)%x- p_vec_c(il_c1,jk,ib_c1)%x)&
          &*p_patch%edges%inv_dual_edge_length(je,jb)))

       ELSE
         p_vn_e(je,jk,jb) = 0.0_wp
       ENDIF
! IF(jb==3)then
! write(*,*)'vn, ptp_vn:',je,jk,jb,&
! & il_c1,ib_c1, p_vn_c(il_c1,jk,ib_c1)%x,&
! & il_c1,ib_c1, p_vn_c(il_c2,jk,ib_c2)%x
! ENDIF
    END DO EDGE_IDX_LOOP
  END DO LEVEL_LOOP_E
END DO EDGE_BLK_LOOP

  END SUBROUTINE map_cell2edges_upwind2
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_cell2edges_2D( p_patch, p_vn_c, ptp_vn)

  TYPE(t_patch), INTENT(IN)                 :: p_patch          ! patch on which computation is performed
  TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(INOUT)                   :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)  


  !Local variables
  INTEGER :: rl_start_e, rl_end_e 
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: je, jb
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
  rl_start_e   = 1
  rl_end_e     = min_rledge
  i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)


  ! calculation of transposed P^TPv from Pv (incart coord)
  EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb,&
                   & i_startblk_e, i_endblk_e,&
                   & i_startidx_e, i_endidx_e,&
                   & rl_start_e, rl_end_e)

    EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e

      IF(v_base%lsm_oce_e(je,1,jb) <= sea_boundary)THEN

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        ptp_vn(je,jb) =&
        & DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
        &+DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x,v_base%edge2cell_coeff_cc_t(je,jb,2)%x)
       ELSE
         ptp_vn(je,jb) = 0.0_wp
       ENDIF
    END DO EDGE_IDX_LOOP
  END DO EDGE_BLK_LOOP

  END SUBROUTINE map_cell2edges_2D
  !-----------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  SUBROUTINE map_edges2cell_with_height( p_patch, vn_e, p_vn_c, h_e, opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN)                  :: p_patch        ! patch on which computation is performed
  REAL(wp), INTENT(IN)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
                                                               ! 3D case: h_e is surface elevation at edges
  TYPE(t_cartesian_coordinates),INTENT(INOUT):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,nblks_c)
  REAL(wp), INTENT(IN)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges 
  INTEGER, INTENT(IN), OPTIONAL              :: opt_slev       ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL              :: opt_elev       ! optional vertical end level

  !Local variables
  INTEGER, PARAMETER :: no_cell_edges = 3
  INTEGER :: slev, elev
  INTEGER :: rl_start_c, rl_end_c!, rl_start_e, rl_end_e, rl_start_v, rl_end_v 
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: il_e, ib_e
  INTEGER :: jc, jb, jk, ie!,je
  REAL(wp) :: z_weight
  REAL(wp) :: z_thick_e
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)


IF ( iswm_oce == 1 ) THEN

  !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
  CELL_BLK_LOOP_SWM: DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb,&
                      & i_startblk_c, i_endblk_c,&
                      & i_startidx_c, i_endidx_c,&
                      & rl_start_c, rl_end_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    LEVEL_LOOP_SWM: DO jk = slev, elev
      CELL_IDX_LOOP_SWM: DO jc =  i_startidx_c, i_endidx_c
        !calculate velocity reconstruction at cell center
        z_weight           = 0.0_wp
        p_vn_c(jc,jk,jb)%x = 0.0_wp
        DO ie=1, no_cell_edges

          il_e = p_patch%cells%edge_idx(jc,jb,ie)
          ib_e = p_patch%cells%edge_blk(jc,jb,ie)

          z_weight = z_weight + v_base%variable_vol_norm(jc,jb,ie)*h_e(il_e,ib_e)
          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                           & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
                           & * vn_e(il_e,jk,ib_e) * h_e(il_e,ib_e)
        END DO
        IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x / z_weight
        ELSE
          p_vn_c(jc,jk,jb)%x=0.0_wp
        ENDIF
      END DO CELL_IDX_LOOP_SWM
    END DO LEVEL_LOOP_SWM
  END DO CELL_BLK_LOOP_SWM

ELSEIF( iswm_oce /= 1 ) THEN 

  !Step 1: Calculation of Pv in cartesian coordinates 
  CELL_BLK_LOOP: DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb,&
                      & i_startblk_c, i_endblk_c,&
                      & i_startidx_c, i_endidx_c,&
                      & rl_start_c, rl_end_c)

    !We are dealing with the surface layer first
    CELL_IDX_LOOP_TOP: DO jc =  i_startidx_c, i_endidx_c
      z_weight             = 0.0_wp
      p_vn_c(jc,slev,jb)%x = 0.0_wp
      DO ie=1, no_cell_edges

        il_e = p_patch%cells%edge_idx(jc,jb,ie)
        ib_e = p_patch%cells%edge_blk(jc,jb,ie)

        z_thick_e = v_base%del_zlev_m(slev) + h_e(il_e,ib_e) 

        z_weight = z_weight + v_base%variable_vol_norm(jc,jb,ie)*z_thick_e
        p_vn_c(jc,slev,jb)%x = p_vn_c(jc,slev,jb)%x&
                           & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
                           & * vn_e(il_e,slev,ib_e) * z_thick_e
! IF(jb==3)THEN
! write(*,*)'map 2 Step:terms',jc,jb,ie,vn_e(il_e,slev,ib_e), z_thick_e,&
! &v_base%edge2cell_coeff_cc(jc,jb,ie)%x
! ENDIF
      END DO
! IF(jb==3)THEN
! write(*,*)'map 2 Step:Vec',jc,jb, p_vn_c(jc,slev,jb)%x 
! ENDIF

      IF(v_base%lsm_oce_c(jc,slev,jb) <= sea_boundary)THEN
        p_vn_c(jc,slev,jb)%x = p_vn_c(jc,slev,jb)%x / z_weight
      ELSE
       p_vn_c(jc,slev,jb)%x=0.0_wp 
      ENDIF
    END DO CELL_IDX_LOOP_TOP

    !Now we calculate at the levels below the surface
    LEVEL_LOOP: DO jk = slev+1, elev
      CELL_IDX_LOOP: DO jc =  i_startidx_c, i_endidx_c
        p_vn_c(jc,jk,jb)%x = 0.0_wp
        DO ie=1, no_cell_edges
          il_e = p_patch%cells%edge_idx(jc,jb,ie)
          ib_e = p_patch%cells%edge_blk(jc,jb,ie)

          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                             & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
                             & * vn_e(il_e,jk,ib_e)
        END DO
        IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x/v_base%fixed_vol_norm(jc,jb)
        ELSE
          p_vn_c(jc,jk,jb)%x = 0.0_wp
       ENDIF
      END DO CELL_IDX_LOOP
    END DO LEVEL_LOOP
  END DO CELL_BLK_LOOP
ENDIF
! IF(    v_base%lsm_oce_c(il_c1,jk,ib_c1) == sea_boundary
  END SUBROUTINE map_edges2cell_with_height
!----------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!map_edges2cell_without_height
  SUBROUTINE map_edges2cell_without_height( p_patch, vn_e, p_vn_c, opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN)                  :: p_patch        ! patch on which computation is performed
  REAL(wp)                                   :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
                                                               ! 3D case: h_e is surface elevation at edges
  TYPE(t_cartesian_coordinates),INTENT(INOUT):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,nblks_c)
  INTEGER, INTENT(IN), OPTIONAL              :: opt_slev       ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL              :: opt_elev       ! optional vertical end level

  !Local variables
  INTEGER, PARAMETER :: no_cell_edges = 3
  INTEGER :: slev, elev
  INTEGER :: rl_start_c, rl_end_c!, rl_start_e, rl_end_e, rl_start_v, rl_end_v 
  INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
  !INTEGER :: i_startblk_e, i_endblk_e!, i_startidx_e, i_endidx_e
  !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: il_e, ib_e
  INTEGER :: jc, jb, jk, ie
  !REAL(wp) :: z_weight
  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_e2c')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

  !Calculation of Pv in cartesian coordinates
  CELL_BLK_LOOP: DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb,&
                      & i_startblk_c, i_endblk_c,&
                      & i_startidx_c, i_endidx_c,&
                      & rl_start_c, rl_end_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    LEVEL_LOOP: DO jk = slev, elev
      CELL_IDX_LOOP: DO jc =  i_startidx_c, i_endidx_c
        !calculate velocity reconstruction at cell center
        p_vn_c(jc,jk,jb)%x = 0.0_wp
        DO ie=1, no_cell_edges
          il_e = p_patch%cells%edge_idx(jc,jb,ie)
          ib_e = p_patch%cells%edge_blk(jc,jb,ie)

          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                           & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
                           & * vn_e(il_e,jk,ib_e) 
        END DO
       IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x / v_base%fixed_vol_norm(jc,jb)
        ELSE
          p_vn_c(jc,jk,jb)%x=0.0_wp
        ENDIF
      END DO CELL_IDX_LOOP
    END DO LEVEL_LOOP
  END DO CELL_BLK_LOOP
  END SUBROUTINE map_edges2cell_without_height
  !-----------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based into edge-based vectors on the primal grid.
  !!
  !! Discrete mapping of cell-based into edge-based vectors on the primal grid.
  !! This mapping depends on the grid geometry only and is used to formulate
  !! the primitive equations in weak form.
  !! The coefficients are calculated in sbrt init_scalar_product and stored in
  !! the ocean part of the patch.
  !! Input lives on cells (2-dim, 3-dim), output lives on edges (3-dim)
  !!
  !! @par Revision History
  !!  developed by Stephan Lorenz, MPI-M (2010-11)
  !!
  SUBROUTINE primal_map_c2e( ppatch, u_c, v_c, vn_e, h_c, opt_slev, opt_elev )

  TYPE(t_patch), INTENT(IN) :: ppatch         ! patch on which computation is performed
  REAL(wp),    INTENT(IN)   :: u_c(:,:,:)     ! zonal input component (nproma,n_zlev,nblks_c)
  REAL(wp),    INTENT(IN)   :: v_c(:,:,:)     ! meridional input component (nproma,n_zlev,nblks_c)
  REAL(wp),    INTENT(OUT)  :: vn_e(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
  REAL(wp),INTENT(IN), OPTIONAL :: h_c(:,:)    ! SWE-case: thickness at cell centers
  INTEGER, INTENT(IN), OPTIONAL ::  opt_slev  ! optional vertical start level
  INTEGER, INTENT(IN), OPTIONAL ::  opt_elev  ! optional vertical end level

  !Local variables
  TYPE(t_cartesian_coordinates)   :: zu_cc(nproma,ppatch%nblks_c)
  TYPE(t_cartesian_coordinates)   :: zu_cc2(nproma,ppatch%nblks_c)
  INTEGER :: slev, elev
  INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: jb, jc, je, jk, ie

  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !  & routine = ('mo_scalar_product:primal_map_c2e')
  !-----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')

  DO ie = 1,3
    zu_cc(:,:)%x(ie) = 0.0_wp
  END DO

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = n_zlev
  END IF

  ! #slo# due to nag -nan compiler-option - remove for production
  vn_e(:,slev:elev,:) = 0.0_wp

  LEVEL_LOOP: DO jk = slev, elev

    rl_start = 1
    rl_end = min_rlcell

    i_startblk = ppatch%cells%start_blk(rl_start,1)
    i_endblk   = ppatch%cells%end_blk(rl_end,1)

    CELL_BLK_LOOP: DO jb = i_startblk, i_endblk

      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      CELL_IDX_LOOP: DO jc =  i_startidx, i_endidx

        CALL gvec2cvec( u_c(jc,jk,jb), v_c(jc,jk,jb),   &
          &             ppatch%cells%center(jc,jb)%lon, &
          &             ppatch%cells%center(jc,jb)%lat, &
          &             zu_cc(jc,jb)%x(1),zu_cc(jc,jb)%x(2), zu_cc(jc,jb)%x(3))

! IF(jk==1)THEN
! IF(jc==1.and.jb<=10)THEN
!  write(*,*)'P-cart-old:',jc,jb, zu_cc(jc,jb)%x
! ENDIF
! ENDIF
      END DO CELL_IDX_LOOP
    END DO CELL_BLK_LOOP


    rl_start = 1
    rl_end = min_rledge

    i_startblk = ppatch%edges%start_blk(rl_start,1)
    i_endblk   = ppatch%edges%end_blk(rl_end,1)

    EDGE_BLK_LOOP: DO jb = i_startblk, i_endblk

      CALL get_indices_e(ppatch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      EDGE_IDX_LOOP: DO je =  i_startidx, i_endidx

        !Get indices of two adjacent triangles
        il_c1 = ppatch%edges%cell_idx(je,jb,1)
        ib_c1 = ppatch%edges%cell_blk(je,jb,1)
        il_c2 = ppatch%edges%cell_idx(je,jb,2)
        ib_c2 = ppatch%edges%cell_blk(je,jb,2)

        IF(jk==1)THEN
          IF ( PRESENT(h_c) ) THEN
!IF(jb==17)then
!write(*,*)'c2e:vec',jb,je,zu_cc(il_c1,ib_c1)%x, zu_cc(il_c2,ib_c2)%x,&
!&h_c(il_c1,ib_c1),h_c(il_c2,ib_c2) 
!ENDIF
            !Multiply cartesian vector with height vaiable 
            zu_cc2(il_c1,ib_c1)%x = zu_cc(il_c1,ib_c1)%x*h_c(il_c1,ib_c1)
            zu_cc2(il_c2,ib_c2)%x = zu_cc(il_c2,ib_c2)%x*h_c(il_c2,ib_c2)

            IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
              vn_e(je,jk,jb) = &
       &  DOT_PRODUCT(zu_cc2(il_c1,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
       &+ DOT_PRODUCT(zu_cc2(il_c2,ib_c2)%x,v_base%edge2cell_coeff_cc_t(je,jb,2)%x)
            ELSE
              vn_e(je,jk,jb) = 0.0_wp
            ENDIF 

! !IF(je=12.and.jb==17)then
! IF(jb==17)then
! write(*,*)'c2e:vec*h',jb,je,zu_cc2(il_c1,ib_c1)%x, zu_cc2(il_c2,ib_c2)%x
! write(*,*)'c2e:',je,jb,vn_e(je,jk,jb),&
! &DOT_PRODUCT(zu_cc2(il_c1,ib_c1)%x, v_base%edge2cell_coeff_t(je,jb,1)%x),&
! &DOT_PRODUCT(zu_cc2(il_c2,ib_c2)%x, v_base%edge2cell_coeff_t(je,jb,2)%x)
! ENDIF
          ELSEIF(.NOT.PRESENT(h_c))THEN

            IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
              vn_e(je,jk,jb) = &
       &  DOT_PRODUCT(zu_cc(il_c1,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
       &+ DOT_PRODUCT(zu_cc(il_c2,ib_c2)%x,v_base%edge2cell_coeff_cc_t(je,jb,2)%x)
            ELSE
              vn_e(je,jk,jb) = 0.0_wp
            ENDIF
! IF(jb==17)then
! write(*,*)'c2e:vec',je,jb,zu_cc(il_c1,ib_c1)%x, zu_cc(il_c2,ib_c2)%x
!  write(*,*)'c2e:',je,jb,il_c1,ib_c1,il_c2,ib_c2, vn_e(je,jk,jb),&
!  &DOT_PRODUCT(zu_cc(il_c1,ib_c1)%x, v_base%edge2cell_coeff_t(je,jb,1)%x),&
!  &DOT_PRODUCT(zu_cc(il_c2,ib_c2)%x, v_base%edge2cell_coeff_t(je,jb,2)%x)
! ENDIF
          ENDIF
        ELSEIF(jk>1)THEN
          IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
            vn_e(je,jk,jb) =&
           &  DOT_PRODUCT(zu_cc(il_c1,ib_c1)%x, v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
           &+ DOT_PRODUCT(zu_cc(il_c2,ib_c2)%x, v_base%edge2cell_coeff_cc_t(je,jb,2)%x)
          ELSE
            vn_e(je,jk,jb) = 0.0_wp
          ENDIF 
        ENDIF!jk-condition
      END DO EDGE_IDX_LOOP
    END DO EDGE_BLK_LOOP
  END DO LEVEL_LOOP

  END SUBROUTINE primal_map_c2e


END MODULE mo_scalar_product

