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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_testbed_operators
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_run_config,         ONLY: test_mode
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_impl_constants,     ONLY: boundary, sea, sea_boundary !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  !USE mo_exception,          ONLY: finish, message
#ifndef __SX__
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
#endif
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, vector_product 
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  USE mo_ocean_math_operators

  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_datetime,            ONLY: t_datetime
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_sea_ice_types,       ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, t_sea_ice
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=12)           :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: ocean_test_operators

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_operators( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data,                     &
    & surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),TARGET,   INTENT(inout)   :: operators_coefficients

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_test_operators"

    SELECT CASE (test_mode) ! 100 - 999

      CASE (101)
        CALL test_nonlinear_coriolis_3d(patch_3D, ocean_state(1), operators_coefficients)
      CASE (102)
          CALL operator_test(patch_3D, ocean_state(1), operators_coefficients)
      
      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT

  END SUBROUTINE ocean_test_operators
  !-------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE operator_test( patch_3d, ocean_state, operators_coefficients)!, vn_e, trac_c)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),TARGET, INTENT(in)   :: operators_coefficients
    
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
	REAL(wp), POINTER :: vn_e(:,:,:)!(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
	REAL(wp), POINTER :: trac_c(:,:,:)!(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
		
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    verts_in_domain => patch_2D%verts%in_domain
    !-------------------------------------------------------------------------------
	vn_e  => ocean_state%p_prog(1)%vn
	trac_c=>ocean_state%p_prog(1)%tracer(:,:,:,1) 	
	
    grad    (1:nproma,1:n_zlev,1:patch_2D%nblks_e)=0.0_wp
    div     (1:nproma,1:n_zlev,1:patch_2D%nblks_c)=0.0_wp
    curl    (1:nproma,1:n_zlev,1:patch_2D%nblks_v)=0.0_wp
    curlgrad(1:nproma,1:n_zlev,1:patch_2D%nblks_v)=0.0_wp
    
    
    !vn_e=0.0_wp
    vn_e(1:5,1,1)=10.0_wp
    trac_c=35.0_wp
    CALL map_edges2vert_3d(patch_2D, vn_e, operators_coefficients%edge2vert_coeff_cc, vn_dual)
    
    !calculate gradient of curl and curl
    CALL grad_fd_norm_oce_3d( trac_c, patch_3D, operators_coefficients%grad_coeff, grad)
    CALL rot_vertex_ocean_3d( patch_3D, grad, vn_dual, operators_coefficients, curlgrad)
    CALL rot_vertex_ocean_3d( patch_3D, vn_e, vn_dual, operators_coefficients, curl)
    
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
    CALL div_oce_3d( vn_e, patch_3D,operators_coefficients%div_coeff, div)
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
    CALL div_oce_3d( vn_e, patch_3D,operators_coefficients%div_coeff, div)
    CALL grad_fd_norm_oce_3d( trac_c, patch_3D, operators_coefficients%grad_coeff, grad)

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

  !-------------------------------------------------------------------------
  !>
  !! Note:at the momode it only works in sequential mode !
  SUBROUTINE test_nonlinear_coriolis_3d(patch_3D, ocean_state, operators_coefficients)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_operator_coeff),TARGET, INTENT(in)   :: operators_coefficients

    REAL(wp) :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), POINTER :: vn(:,:,:)
    TYPE(t_cartesian_coordinates), POINTER  :: vn_dual(:,:,:)
    TYPE(t_patch), POINTER :: patch_2d

    patch_2d   => patch_3d%p_patch_2d(1)
    vn      => ocean_state%p_prog(nold(1))%vn
    vn_dual => ocean_state%p_diag%p_vn_dual
    
    vn = 1.0_wp
    vort_v  = 0.0_wp
    patch_2d%verts%f_v = 1.0_wp    
    
    ! vn_dual = 0.0_wp
    ! CALL rot_vertex_ocean_3d( patch_3d, vn, p_vn_dual, operators_coefficients, vort_v)
    
    write(0,*) "-------------- original coriolis budget -----------------------"
    CALL testbed_nonlinear_coriolis_3d(patch_3d, vn, vort_v, &
       & operators_coefficients, vort_flux)
    CALL diagnose_coriolis_energy(patch_3d, vn, vort_flux)

    write(0,*) "-------------- fast coriolis budget -----------------------"
    CALL testbed_nonlinear_coriolis_3d_fast(patch_3d, vn, vort_v, &
       & operators_coefficients, vort_flux)
    CALL diagnose_coriolis_energy(patch_3d, vn, vort_flux)

  END SUBROUTINE test_nonlinear_coriolis_3d
  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o3 will produce very different results when using this version,
  !!  the model exhibits great sensitivity using the AtlanticBoxACC setup
  !<Optimize:inUse>
  SUBROUTINE testbed_nonlinear_coriolis_3d_fast(patch_3d, vn, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    ! TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel! , endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, vertex_edge
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    REAL(wp) :: this_vort_flux(n_zlev, 2) ! for each of the two vertices

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    ! endLevel    = n_zlev


!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk,  &
!ICON_OMP vertex_edge, ictr) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      ! vort_flux(:,:,blockNo) = 0.0_wp

      edge_idx_loop: DO je =  start_edge_index, end_edge_index

        this_vort_flux(:,:) = 0.0_wp

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

        ! vertex 1
        ictr = 0
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex1_idx,vertex1_blk)!no_dual_cell_edges

          ictr =ictr+1

          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

            this_vort_flux(level, 1) =  this_vort_flux(level, 1) + &
              & vn( patch_2d%verts%edge_idx(vertex1_idx,vertex1_blk,vertex_edge), level, &
              &      patch_2d%verts%edge_blk(vertex1_idx,vertex1_blk,vertex_edge))  * &
              &  operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)

          ENDDO

        END DO ! edges of this vertex

        ! vertex 2
        ictr = no_dual_edges
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex2_idx,vertex2_blk)!no_dual_cell_edges

          ictr =ictr+1

          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

            this_vort_flux(level, 2) =  this_vort_flux(level, 2) + &
              & vn( patch_2d%verts%edge_idx(vertex2_idx,vertex2_blk,vertex_edge), level, &
              &      patch_2d%verts%edge_blk(vertex2_idx,vertex2_blk,vertex_edge))  * &
              &  operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)

          ENDDO

        END DO ! edges of this vertex

        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          vort_flux(je,level,blockNo) =  &
            & this_vort_flux(level,1) * &
            & (vort_v(vertex1_idx, level, vertex1_blk) + patch_2d%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(level,2) * &
            & (vort_v(vertex2_idx, level, vertex2_blk) + patch_2d%verts%f_v(vertex2_idx, vertex2_blk))

        ENDDO

      END DO edge_idx_loop

    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE testbed_nonlinear_coriolis_3d_fast
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Note: it only works in sequential mode !
  SUBROUTINE testbed_nonlinear_coriolis_3d(patch_3d, vn, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    ! TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    REAL(wp) :: vort_global, thick_edge, thick_vert

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    endLevel    = n_zlev


    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      level_loop: DO level = startLevel, endLevel

        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN

            vort_flux(je,level,blockNo) = 0.0_wp

            DO neighbor=1,2
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = patch_2d%edges%vertex_idx(je,blockNo,neighbor)
              ib_v = patch_2d%edges%vertex_blk(je,blockNo,neighbor)

              vort_global = (vort_v(il_v,level,ib_v) + patch_2d%verts%f_v(il_v,ib_v))

              thick_vert=0.0_wp
              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)
                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)

                thick_edge = patch_3D%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
                thick_vert = thick_vert+thick_edge/patch_2d%verts%num_edges(il_v,ib_v)
              END DO

              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1

                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)

                thick_edge = patch_3D%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
                vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)+vn(il_e,level,ib_e)*vort_global&
                  & *operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)*(thick_edge/thick_vert)
              END DO
            END DO
            
          ELSE
            vort_flux(je,level,blockNo)= 0.0_wp
          ENDIF
          
        END DO edge_idx_loop
      END DO level_loop
    END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block
  END SUBROUTINE testbed_nonlinear_coriolis_3d
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Note: it only works in sequential mode !
  SUBROUTINE diagnose_coriolis_energy(patch_3d, vn, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in)        :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    REAL(wp) :: vort_flux_budget_perlevel(n_zlev), volume(n_zlev)
    REAL(wp) :: vort_flux_budget_perlevel_pos(n_zlev),vort_flux_budget_perlevel_neg(n_zlev)
    REAL(wp) :: in_product
    
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel    = 1
    endLevel    = n_zlev
    !-----------------------------------------------------------------------
    !Diagnostic of energetic neutrality
    vort_flux_budget_perlevel(1:n_zlev)=0.0_wp
    vort_flux_budget_perlevel_pos(1:n_zlev)=0.0_wp
    vort_flux_budget_perlevel_neg(1:n_zlev)=0.0_wp
    volume(1:n_zlev)=0.0_wp

    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      DO je =  start_edge_index, end_edge_index
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          in_product = vort_flux(je,level,blockNo) * vn(je,level,blockNo)
          
          vort_flux_budget_perlevel(level)=vort_flux_budget_perlevel(level)&
            & +in_product &
            & *patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
            & *patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)

          IF(in_product >0.0_wp) THEN
            vort_flux_budget_perlevel_pos(1:n_zlev)= vort_flux_budget_perlevel_pos(1:n_zlev)&
              & + in_product &
              &   * patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
              &   * patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
          ELSEIF(in_product <0.0_wp) THEN
            vort_flux_budget_perlevel_neg(1:n_zlev)= vort_flux_budget_perlevel_neg(1:n_zlev)&
              & + in_product &
              & * patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
              & * patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
          ENDIF
          
          volume(level)=volume(level)+patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
            & *patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
          
       END DO
      END DO
    END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block

    DO level = startLevel, endLevel
      write(0,*)'global vorticity budget', level, &
        & vort_flux_budget_perlevel(level)/volume(level)
      write(0,*)'global vorticity budget pos:neg', level,&
        & vort_flux_budget_perlevel_pos(level), vort_flux_budget_perlevel_neg(level), &
        & (vort_flux_budget_perlevel_pos(level) + vort_flux_budget_perlevel_neg(level))/volume(level)
    END DO
  END SUBROUTINE diagnose_coriolis_energy
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_operators

