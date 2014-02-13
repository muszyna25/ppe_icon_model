!>
!!   Contains the implementation of the mathematical operators for the ocean.
!!
!!   Contains the implementation of the mathematical operators for the ocean.
!!
!! @par Revision History
!!  Developed  by Peter Korn and Stephan Lorenz 2010-04
!!  Modified by Stephan Lorenz                  2011-02
!!    correct implementation of ocean boundaries
!!
!! @par To Do
!! Boundary exchange, nblks in presence of halos and dummy edge
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
MODULE mo_oce_operators_testbed
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
!  USE mo_exception,          ONLY: finish
  USE mo_run_config,         ONLY: ltimer
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_impl_constants,     ONLY: boundary, sea_boundary !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  !USE mo_exception,          ONLY: finish, message
#ifndef __SX__
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
#endif
  USE mo_oce_types,          ONLY: t_hydro_ocean_state
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, vector_product !, gc2cc
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  USE mo_oce_math_operators

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=12)           :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: operator_test

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE operator_test( patch_3d, p_os, p_op_coeff, vn_e, trac_c)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: p_op_coeff
    REAL(wp)                             :: vn_e(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp)                             :: trac_c(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    
    !
    !
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, je, jv, jk, jb,edge_index,edge_block
    REAL(wp) :: curl_integral(1:n_zlev), div_integral(1:n_zlev), lhs(1:n_zlev),rhs(1:n_zlev)
    REAL(wp) :: grad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp) :: curl(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp) :: curlgrad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_v)    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain, verts_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_cartesian_coordinates) :: vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    verts_in_domain => patch_2D%verts%in_domain
    !-------------------------------------------------------------------------------
    grad    (1:nproma,1:n_zlev,1:patch_2D%nblks_e)=0.0_wp
    div     (1:nproma,1:n_zlev,1:patch_2D%nblks_c)=0.0_wp
    curl    (1:nproma,1:n_zlev,1:patch_2D%nblks_v)=0.0_wp
    curlgrad(1:nproma,1:n_zlev,1:patch_2D%nblks_v)=0.0_wp
    
    
    !vn_e=0.0_wp
    !vn_e(1:5,1,1)=10.0_wp
    !trac_c=35.0_wp
    CALL map_edges2vert_3d(patch_2D, vn_e, p_op_coeff%edge2vert_coeff_cc, vn_dual)
    
    !calculate gradient of curl and curl
    CALL grad_fd_norm_oce_3d( trac_c, patch_3D, p_op_coeff%grad_coeff, grad)
    CALL rot_vertex_ocean_3d( patch_3D, grad, vn_dual, p_op_coeff, curlgrad)
    CALL rot_vertex_ocean_3d( patch_3D, vn_e, vn_dual, p_op_coeff, curl)
    
    !domain integrated curl
    curl_integral=0.0_wp
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, start_index, end_index)
      DO jk = 1, n_zlev
        DO jv = start_index, end_index
        curl_integral(jk) = curl_integral(jk) + curl(jv,jk,jb)*patch_2D%verts%dual_area(jv,jb)

        END DO
      END DO
    END DO
    write(0,*)'OPERATOR-TEST:curl-of-gradient:',maxval(curlgrad),minval(curlgrad)
    Do jk=1,n_zlev
    write(0,*)'OPERATOR-TEST:curl-integral:',jk,curl_integral(jk)
    END DO
    !calculate domain integrated div
    CALL div_oce_3d( vn_e, patch_2D,p_op_coeff%div_coeff, div)
    div_integral=0.0_wp
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          DO je=1,3
          edge_index = patch_2D%cells%edge_idx(jc,jb, je)
          edge_block = patch_2D%cells%edge_blk(jc,jb,je)
          
             div_integral(jk)= div_integral(jk)+(vn_e(edge_index,jk, edge_block)* &
             &patch_2D%edges%primal_edge_length(edge_index, edge_block)/patch_2D%cells%area(jc,jb)       * &
            & patch_2D%cells%edge_orientation(jc,jb, je))* patch_2D%cells%area(jc,jb)
IF(jk==1)&
&write(123,*)'details:',jk,jc,jb,je,edge_index,edge_block,vn_e(edge_index,jk, edge_block)* &
!             &patch_2D%edges%primal_edge_length(edge_index, edge_block)        * &
            & patch_2D%cells%edge_orientation(jc,jb, je)             
          END DO
          div_integral(jk)= div_integral(jk)
IF(jk==1)&          
&write(123,*)'details2',jk,jc,jb,div_integral(jk)         
!          div_integral(jk) = div_integral(jk) + div(jc,jk,jb)*patch_2D%cells%area(jc,jb)                            
        END DO
      END DO
    END DO
    write(123,*)'-----------------------'
    DO jk=1,n_zlev
    write(0,*)'OPERATOR-TEST: div-integral:',jk,div_integral(jk)
    END DO

    !product rule for divergence
    CALL div_oce_3d( vn_e, patch_2D,p_op_coeff%div_coeff, div)
    CALL grad_fd_norm_oce_3d( trac_c, patch_3D, p_op_coeff%grad_coeff, grad)

    lhs=0.0_wp
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jk = 1, n_zlev
        DO jc = start_index, end_index    
          lhs(jk)=lhs(jk)+trac_c(jc,jk,jb)*div(jc,jk,jb)*patch_2D%cells%area(jc,jb)        
        END DO
      END DO
    END DO
    rhs=0.0_wp
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      DO jk = 1, n_zlev
        DO je = start_index, end_index
        rhs(jk)=rhs(jk)+grad(je,jk,jb)*vn_e(je,jk,jb)*patch_2D%edges%area_edge(je,jb)
        !rhs=rhs+grad(je,jk,jb)*vn_e(je,jk,jb)&
        !&*patch_2D%edges%primal_edge_length(je,jb)*patch_2D%edges%dual_edge_length(je,jb)
        END DO
      END DO
    END DO
    Do jk=1,n_zlev
      write(0,*)'OPERATOR-TEST: integration-by-parts:',jk,lhs(jk),rhs(jk),lhs(jk)+rhs(jk)          
      write(0,*)'OPERATOR-TEST: integration-by-parts2:',jk,&
      &maxval(grad(:,jk,:)),minval(grad(:,jk,:)),maxval(vn_e(:,jk,:)),minval(vn_e(:,jk,:))          
      
    END DO  
  END SUBROUTINE operator_test
  !-------------------------------------------------------------------------------

END MODULE mo_oce_operators_testbed

