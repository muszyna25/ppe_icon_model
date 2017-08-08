!>

!! Type definition for the GPU implementation of the dynamical core of ICONAM.
!!
!! @author William Sawyer (CSCS)
!!
!! @par Revision History
!! Initial release by William Sawyer (2015)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_nonhydro_gpu_types

  USE mo_kind,                 ONLY: wp, vp
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d
  USE mo_run_config,           ONLY: iforcing, ntracer, iqm_max,                &
    &                                iqv, iqc, iqi, iqr, iqs, iqt, iqtvar,      &
    &                                iqni, iqni_nuc, iqg, iqh, iqnr, iqns,      & 
    &                                iqng, iqnh, iqnc, inccn, ininpot, ininact, &
    &                                iqtke, nqtendphy, ltestcase, lart   
  USE mo_model_domain,         ONLY: t_patch, p_patch
  USE mo_nonhydro_state,       ONLY: p_nh_state
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_diag, t_nh_prog
  USE mo_nh_prepadv_types,     ONLY: t_prepare_adv

  IMPLICIT NONE

  PRIVATE 

#if defined( _OPENACC )

  PUBLIC :: my_patches
  PUBLIC :: save_patch_pointers, refresh_patch_pointers
  PUBLIC :: save_convenience_pointers, refresh_convenience_pointers
  PUBLIC :: h2d_solve_nonhydro, d2h_solve_nonhydro
  PUBLIC :: h2d_velocity_tendencies, d2h_velocity_tendencies

! POINTER information which needs to be stored and refreshed on the CPU

  TYPE my_patches
    TYPE ( t_patch ), POINTER :: all_cells_patch, all_edges_patch, all_verts_patch
    TYPE ( t_patch ), POINTER :: owned_cells_patch, owned_edges_patch, owned_verts_patch
    TYPE ( t_patch ), POINTER :: in_domain_cells_patch, in_domain_edges_patch, in_domain_verts_patch
    TYPE ( t_patch ), POINTER :: not_owned_cells_patch, not_owned_edges_patch, not_owned_verts_patch
    TYPE ( t_patch ), POINTER :: not_in_domain_cells_patch, not_in_domain_edges_patch, not_in_domain_verts_patch
    TYPE ( t_patch ), POINTER :: gradiscalculable_edges_patch
    TYPE ( t_patch ), POINTER :: one_edge_in_domain_cells_patch
  ENDTYPE my_patches

  TYPE my_nh_state
    TYPE ( t_nh_prog ), ALLOCATABLE  :: prog(:)
    TYPE ( t_nh_diag )               :: diag
  ENDTYPE my_nh_state

  TYPE ( my_patches ), ALLOCATABLE :: save_patches(:)

  TYPE ( my_nh_state ), ALLOCATABLE :: save_nh_state(:)

  CHARACTER(len=*), PARAMETER :: version = &
    & '$Id: mo_nonhydro_gpu_types.f90 16760 2014-04-01 09:04:22Z wsawyer $'

CONTAINS

     SUBROUTINE save_patch_pointers( p_patch, save_patch )

       TYPE ( t_patch ),    INTENT(INOUT) :: p_patch
       TYPE ( my_patches ), INTENT(INOUT) :: save_patch

     
!
! Save the critical aliases in one patch
!

       save_patch%all_cells_patch =>  p_patch%cells%all%patch
       save_patch%all_edges_patch =>  p_patch%edges%all%patch
       save_patch%all_verts_patch =>  p_patch%verts%all%patch

       NULLIFY( p_patch%cells%decomp_info%halo_level )
       NULLIFY( p_patch%edges%decomp_info%halo_level )
       NULLIFY( p_patch%verts%decomp_info%halo_level )

       save_patch%all_cells_patch =>  p_patch%cells%all%patch
       save_patch%all_edges_patch =>  p_patch%edges%all%patch
       save_patch%all_verts_patch =>  p_patch%verts%all%patch

       NULLIFY( p_patch%cells%all%patch )
       NULLIFY( p_patch%edges%all%patch )
       NULLIFY( p_patch%verts%all%patch )

       save_patch%owned_cells_patch =>  p_patch%cells%owned%patch
       save_patch%owned_edges_patch =>  p_patch%edges%owned%patch
       save_patch%owned_verts_patch =>  p_patch%verts%owned%patch

       NULLIFY( p_patch%cells%owned%patch )
       NULLIFY( p_patch%edges%owned%patch )
       NULLIFY( p_patch%verts%owned%patch )

       save_patch%in_domain_cells_patch =>  p_patch%cells%in_domain%patch
       save_patch%in_domain_edges_patch =>  p_patch%edges%in_domain%patch
       save_patch%in_domain_verts_patch =>  p_patch%verts%in_domain%patch

       NULLIFY( p_patch%cells%in_domain%patch )
       NULLIFY( p_patch%edges%in_domain%patch )
       NULLIFY( p_patch%verts%in_domain%patch )

       save_patch%gradiscalculable_edges_patch =>  p_patch%edges%gradiscalculable%patch

       NULLIFY( p_patch%edges%gradiscalculable%patch )

       save_patch%not_owned_cells_patch =>  p_patch%cells%not_owned%patch
       save_patch%not_owned_edges_patch =>  p_patch%edges%not_owned%patch
       save_patch%not_owned_verts_patch =>  p_patch%verts%not_owned%patch

       NULLIFY( p_patch%cells%not_owned%patch )
       NULLIFY( p_patch%edges%not_owned%patch )
       NULLIFY( p_patch%verts%not_owned%patch )

       save_patch%not_in_domain_cells_patch =>  p_patch%cells%not_in_domain%patch
       save_patch%not_in_domain_edges_patch =>  p_patch%edges%not_in_domain%patch
       save_patch%not_in_domain_verts_patch =>  p_patch%verts%not_in_domain%patch

       NULLIFY( p_patch%cells%not_in_domain%patch )
       NULLIFY( p_patch%edges%not_in_domain%patch )
       NULLIFY( p_patch%verts%not_in_domain%patch )

       save_patch%one_edge_in_domain_cells_patch =>  p_patch%cells%one_edge_in_domain%patch

       NULLIFY( p_patch%cells%one_edge_in_domain%patch )

     END SUBROUTINE save_patch_pointers

     SUBROUTINE refresh_patch_pointers( save_patch, p_patch )

       TYPE ( my_patches ), INTENT(INOUT) :: save_patch
       TYPE ( t_patch ),    INTENT(INOUT) :: p_patch

!
! Refresh the patch pointers only
!
       p_patch%cells%all%patch =>  save_patch%all_cells_patch
       p_patch%edges%all%patch =>  save_patch%all_edges_patch 
       p_patch%verts%all%patch =>  save_patch%all_verts_patch

       NULLIFY( save_patch%all_cells_patch )
       NULLIFY( save_patch%all_edges_patch )
       NULLIFY( save_patch%all_verts_patch )

       p_patch%cells%all%patch =>  save_patch%all_cells_patch
       p_patch%edges%all%patch =>  save_patch%all_edges_patch
       p_patch%verts%all%patch =>  save_patch%all_verts_patch

       NULLIFY( save_patch%all_cells_patch )
       NULLIFY( save_patch%all_edges_patch )
       NULLIFY( save_patch%all_verts_patch )

       p_patch%cells%owned%patch =>  save_patch%owned_cells_patch
       p_patch%edges%owned%patch =>  save_patch%owned_edges_patch
       p_patch%verts%owned%patch =>  save_patch%owned_verts_patch

       NULLIFY( save_patch%owned_cells_patch )
       NULLIFY( save_patch%owned_edges_patch )
       NULLIFY( save_patch%owned_verts_patch )

       p_patch%cells%in_domain%patch =>  save_patch%in_domain_cells_patch 
       p_patch%edges%in_domain%patch =>  save_patch%in_domain_edges_patch 
       p_patch%verts%in_domain%patch =>  save_patch%in_domain_verts_patch 

       NULLIFY( save_patch%in_domain_cells_patch )
       NULLIFY( save_patch%in_domain_edges_patch )
       NULLIFY( save_patch%in_domain_verts_patch )

       p_patch%edges%gradiscalculable%patch =>  save_patch%gradiscalculable_edges_patch

       NULLIFY( save_patch%gradiscalculable_edges_patch )

       p_patch%cells%not_owned%patch =>  save_patch%not_owned_cells_patch
       p_patch%edges%not_owned%patch =>  save_patch%not_owned_edges_patch
       p_patch%verts%not_owned%patch =>  save_patch%not_owned_verts_patch

       NULLIFY( save_patch%not_owned_cells_patch )
       NULLIFY( save_patch%not_owned_edges_patch )
       NULLIFY( save_patch%not_owned_verts_patch )

       p_patch%cells%not_in_domain%patch =>  save_patch%not_in_domain_cells_patch
       p_patch%edges%not_in_domain%patch =>  save_patch%not_in_domain_edges_patch
       p_patch%verts%not_in_domain%patch =>  save_patch%not_in_domain_verts_patch

       NULLIFY( save_patch%not_in_domain_cells_patch )
       NULLIFY( save_patch%not_in_domain_edges_patch )
       NULLIFY( save_patch%not_in_domain_verts_patch )

       p_patch%cells%one_edge_in_domain%patch =>  save_patch%one_edge_in_domain_cells_patch

       NULLIFY( save_patch%one_edge_in_domain_cells_patch )

     END SUBROUTINE refresh_patch_pointers

     SUBROUTINE save_convenience_pointers( )

       INTEGER :: jg, jt, jv, num_vars

!
! Save then nullify convenience pointers on the pHost, in order to perform a deep copy to the device
!

       ALLOCATE( save_patches( SIZE(p_patch) ) )
       DO jg = 1, SIZE(p_patch)
         CALL save_patch_pointers( p_patch(jg), save_patches(jg) )
       ENDDO

       ALLOCATE( save_nh_state( SIZE(p_nh_state) ) )

       DO jg = 1, SIZE(p_nh_state)

! ddt_grf_trc_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%ddt_grf_trc_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%ddt_grf_trc_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_2d => p_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_3d => p_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%ddt_grf_trc_ptr )  ! p_3d not used?
         ENDIF

! hfl_trc_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%hfl_trc_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%hfl_trc_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_2d => p_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_3d => p_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%hfl_trc_ptr )  ! p_3d not used?
         ENDIF

! vfl_trc_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%vfl_trc_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%vfl_trc_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_2d => p_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_3d => p_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%vfl_trc_ptr )  ! p_3d not used?
         ENDIF

! ddt_trc_adv_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%ddt_trc_adv_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%ddt_trc_adv_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%ddt_trc_adv_ptr( jv )%p_2d => p_nh_state(jg)%diag%ddt_trc_adv_ptr( jv )%p_2d 
             save_nh_state(jg)%diag%ddt_trc_adv_ptr( jv )%p_3d => p_nh_state(jg)%diag%ddt_trc_adv_ptr( jv )%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%ddt_trc_adv_ptr )  ! p_2d not used?
         ENDIF

! ddt_vn_adv_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%ddt_vn_adv_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%ddt_vn_adv_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%ddt_vn_adv_ptr( jv )%p_2d => p_nh_state(jg)%diag%ddt_vn_adv_ptr( jv )%p_2d 
             save_nh_state(jg)%diag%ddt_vn_adv_ptr( jv )%p_3d => p_nh_state(jg)%diag%ddt_vn_adv_ptr( jv )%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%ddt_vn_adv_ptr )  ! p_2d not used?
         ENDIF

! ddt_w_adv_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%ddt_w_adv_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%ddt_w_adv_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_2d => p_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_3d => p_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%ddt_w_adv_ptr )  ! p_2d not used?
         ENDIF

! q_int_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%q_int_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%q_int_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%q_int_ptr(jv)%p_2d => p_nh_state(jg)%diag%q_int_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%q_int_ptr(jv)%p_3d => p_nh_state(jg)%diag%q_int_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%q_int_ptr )  ! p_3d not used?
         ENDIF

! q_ubc_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%q_ubc_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%q_ubc_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%q_ubc_ptr(jv)%p_2d => p_nh_state(jg)%diag%q_ubc_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%q_ubc_ptr(jv)%p_3d => p_nh_state(jg)%diag%q_ubc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%q_ubc_ptr )  ! p_3d not used?
         ENDIF

! tracer_vi_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%tracer_vi_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%tracer_vi_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_2d => p_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_3d => p_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%tracer_vi_ptr )  ! p_3d not used?
         ENDIF

! tracer_vi_avg_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%tracer_vi_avg_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( save_nh_state(jg)%diag%tracer_vi_avg_ptr( num_vars ) )
           DO jv = 1, num_vars 
             save_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_2d => p_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_2d 
             save_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_3d => p_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( p_nh_state(jg)%diag%tracer_vi_avg_ptr )  ! p_3d not used?
         ENDIF

         ALLOCATE( save_nh_state(jg)%prog( SIZE( p_nh_state(jg)%prog ) ) )
         DO jt=1, SIZE( p_nh_state(jg)%prog )
           num_vars = SIZE(p_nh_state(jg)%prog(jt)%tracer_ptr)
           IF ( num_vars > 0 ) THEN
             ALLOCATE( save_nh_state(jg)%prog(jt)%tracer_ptr( num_vars ) )
             DO jv = 1, num_vars
               save_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_2d => p_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_2d 
               save_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_3d => p_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_3d 
             ENDDO
             DEALLOCATE( p_nh_state(jg)%prog(jt)%tracer_ptr )  ! p_2d not used?
           ENDIF
         ENDDO
       ENDDO

     END SUBROUTINE save_convenience_pointers

     SUBROUTINE refresh_convenience_pointers( )

       INTEGER :: jg, jt, jv, num_vars

!
! Save then nullify convenience pointers on the Host, in order to perform a deep copy to the device
!

       DO jg = 1, SIZE(p_patch)

         CALL refresh_patch_pointers( save_patches(jg), p_patch(jg) )

       ENDDO

       DEALLOCATE( save_patches )

       DO jg = 1, SIZE(save_nh_state)

! ddt_grf_trc_ptr
         num_vars = SIZE( save_nh_state(jg)%diag%ddt_grf_trc_ptr )
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%ddt_grf_trc_ptr( num_vars ) )  ! p_2d not used?
           DO jv = 1, num_vars
             p_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_2d  => save_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_3d  => save_nh_state(jg)%diag%ddt_grf_trc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%ddt_grf_trc_ptr )
         ENDIF

! hfl_trc_ptr
         num_vars = SIZE( save_nh_state(jg)%diag%hfl_trc_ptr )
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%hfl_trc_ptr( num_vars ) )  ! p_2d not used?
           DO jv = 1, num_vars
             p_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_2d  => save_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_3d  => save_nh_state(jg)%diag%hfl_trc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%hfl_trc_ptr )
         ENDIF

! vfl_trc_ptr
         num_vars = SIZE( save_nh_state(jg)%diag%vfl_trc_ptr )
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%vfl_trc_ptr( num_vars ) )  ! p_2d not used?
           DO jv = 1, num_vars
             p_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_2d  => save_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_3d  => save_nh_state(jg)%diag%vfl_trc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%vfl_trc_ptr )
         ENDIF

! ddt_trc_adv_ptr
         num_vars = SIZE( save_nh_state(jg)%diag%ddt_trc_adv_ptr )
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%ddt_trc_adv_ptr( num_vars ) )  ! p_2d not used?
           DO jv = 1, num_vars
             p_nh_state(jg)%diag%ddt_trc_adv_ptr(jv)%p_2d  => save_nh_state(jg)%diag%ddt_trc_adv_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%ddt_trc_adv_ptr(jv)%p_3d  => save_nh_state(jg)%diag%ddt_trc_adv_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%ddt_trc_adv_ptr )
         ENDIF

! ddt_vn_adv_ptr
         num_vars = SIZE( save_nh_state(jg)%diag%ddt_vn_adv_ptr )
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%ddt_vn_adv_ptr( num_vars ) )  ! p_2d not used?
           DO jv = 1, num_vars
             p_nh_state(jg)%diag%ddt_vn_adv_ptr(jv)%p_2d  => save_nh_state(jg)%diag%ddt_vn_adv_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%ddt_vn_adv_ptr(jv)%p_3d  => save_nh_state(jg)%diag%ddt_vn_adv_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%ddt_vn_adv_ptr )
         ENDIF

! ddt_w_adv_ptr
         num_vars = SIZE( save_nh_state(jg)%diag%ddt_w_adv_ptr )
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%ddt_w_adv_ptr( num_vars ) )  ! p_2d not used?
           DO jv = 1, num_vars
             p_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_2d  => save_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_3d  => save_nh_state(jg)%diag%ddt_w_adv_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%ddt_w_adv_ptr )
         ENDIF

! q_int_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%q_int_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%q_int_ptr(num_vars) )  ! p_3d not used?
           DO jv = 1, num_vars 
             p_nh_state(jg)%diag%q_int_ptr(jv)%p_2d => save_nh_state(jg)%diag%q_int_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%q_int_ptr(jv)%p_3d => save_nh_state(jg)%diag%q_int_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%q_int_ptr )
         ENDIF

! q_ubc_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%q_ubc_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%q_ubc_ptr(num_vars) )  ! p_3d not used?
           DO jv = 1, num_vars 
             p_nh_state(jg)%diag%q_ubc_ptr(jv)%p_2d => save_nh_state(jg)%diag%q_ubc_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%q_ubc_ptr(jv)%p_3d => save_nh_state(jg)%diag%q_ubc_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%q_ubc_ptr )
         ENDIF

! tracer_vi_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%tracer_vi_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%tracer_vi_ptr(num_vars) )  ! p_3d not used?
           DO jv = 1, num_vars 
             p_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_2d => save_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_3d => save_nh_state(jg)%diag%tracer_vi_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%tracer_vi_ptr )
         ENDIF

! tracer_vi_avg_ptr
         num_vars = SIZE(p_nh_state(jg)%diag%tracer_vi_avg_ptr)
         IF ( num_vars > 0 ) THEN
           ALLOCATE( p_nh_state(jg)%diag%tracer_vi_avg_ptr(num_vars) )  ! p_3d not used?
           DO jv = 1, num_vars 
             p_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_2d => save_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_2d 
             p_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_3d => save_nh_state(jg)%diag%tracer_vi_avg_ptr(jv)%p_3d 
           ENDDO
           DEALLOCATE( save_nh_state(jg)%diag%tracer_vi_avg_ptr )
         ENDIF

         DO jt=1, SIZE( p_nh_state(jg)%prog )
           num_vars = SIZE(p_nh_state(jg)%prog(jt)%tracer_ptr)
           IF ( num_vars > 0 ) THEN
             ALLOCATE( p_nh_state(jg)%prog(jt)%tracer_ptr( num_vars ) )
             DO jv = 1, num_vars
               p_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_2d => save_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_2d 
               p_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_3d => save_nh_state(jg)%prog(jt)%tracer_ptr(jv)%p_3d 
             ENDDO
             DEALLOCATE( save_nh_state(jg)%prog(jt)%tracer_ptr )
           ENDIF
         ENDDO
         DEALLOCATE( save_nh_state(jg)%prog )
       ENDDO

       DEALLOCATE( save_nh_state )


     END SUBROUTINE refresh_convenience_pointers

     SUBROUTINE h2d_solve_nonhydro( nnow, jstep, jg, idiv_method, grf_intmethod_e, lprep_adv, l_vert_nested, is_iau_active, p_nh, prep_adv )

       INTEGER, INTENT(IN)       :: nnow, jstep, jg, idiv_method, grf_intmethod_e
       LOGICAL, INTENT(IN)       :: l_vert_nested, lprep_adv, is_iau_active

       TYPE(t_nh_state),          INTENT(INOUT) :: p_nh
       TYPE(t_prepare_adv),       INTENT(INOUT) :: prep_adv

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp                 ! p_prog  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: theta_v_ic_tmp, rho_ic_tmp                                     ! p_diag  WP
       REAL(wp), DIMENSION(:,:),     POINTER  :: dvn_ie_ubc_tmp,  dtheta_v_ic_ubc_tmp, dw_ubc_tmp               ! p_diag  WP 2D
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: mass_fl_e_tmp,  mflx_ic_ubc_tmp, exner_pr_tmp                  ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: grf_bdy_mflx_tmp                                               ! p_diag  WP

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp, w_concorr_c_tmp, ddt_exner_phy_tmp          ! p_diag  VP
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: exner_dyn_incr_tmp, ddt_vn_phy_tmp                             ! p_diag  VP
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: rho_incr_tmp, exner_incr_tmp                                   ! p_diag  VP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp                  ! prep_adv WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_ref_tmp, w_ref_tmp                                          ! p_ref   WP

!
! OpenACC Implementation:  For testing in ACC_VALIDATE=.TRUE. mode, we would ultimately like to be able to run 
!                          this routine entirely on the accelerator with input on the host, and moving
!                          output back to the host.  I order to do this, an additional, far larger, number of fields
!                          must be updated here on the device:
!
! p_nh%prog(nnow)          All present (above)
! p_nh%diag:               ddt_exner_phy, ddt_vn_adv, ddt_vn_phy, ddt_w_adv
!                          vn_ref, dtheta_v_ic_ubc, dw_ubc, dvn_ie_ubc, mflx_ic_ubc
!                          rho_incr, exner_incr, vn_incr, exner_pr
!                          grf_tend_vn, grf_tend_mflx, grf_tend_rho, grf_tend_thv, grf_tend_w
!
! p_nh%metrics:            Entire structure (read-only)
!
! p_patch:                 Entire structure (read-only)
!  

       exner_tmp           => p_nh%prog(nnow)%exner 
       rho_tmp             => p_nh%prog(nnow)%rho
       theta_v_tmp         => p_nh%prog(nnow)%theta_v 
       vn_tmp              => p_nh%prog(nnow)%vn
       w_tmp               => p_nh%prog(nnow)%w
!$ACC UPDATE DEVICE ( exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp )

       rho_ic_tmp          => p_nh%diag%rho_ic
       theta_v_ic_tmp      => p_nh%diag%theta_v_ic
!$ACC UPDATE DEVICE ( rho_ic_tmp, theta_v_ic_tmp )

       vt_tmp              => p_nh%diag%vt
       vn_ie_tmp           => p_nh%diag%vn_ie
       w_concorr_c_tmp     => p_nh%diag%w_concorr_c
!$ACC UPDATE DEVICE ( vt_tmp, vn_ie_tmp, w_concorr_c_tmp )

       mass_fl_e_tmp       => p_nh%diag%mass_fl_e
       exner_pr_tmp        => p_nh%diag%exner_pr
       exner_dyn_incr_tmp  => p_nh%diag%exner_dyn_incr
!$ACC UPDATE DEVICE ( mass_fl_e_tmp, exner_pr_tmp, exner_dyn_incr_tmp )

       mflx_ic_ubc_tmp     => p_nh%diag%mflx_ic_ubc
       dvn_ie_ubc_tmp      => p_nh%diag%dvn_ie_ubc
       dtheta_v_ic_ubc_tmp => p_nh%diag%dtheta_v_ic_ubc
       dw_ubc_tmp          => p_nh%diag%dw_ubc
!$ACC UPDATE DEVICE ( mflx_ic_ubc_tmp, dvn_ie_ubc_tmp, dtheta_v_ic_ubc_tmp, dw_ubc_tmp ) IF( l_vert_nested )

       ddt_exner_phy_tmp   => p_nh%diag%ddt_exner_phy
       ddt_vn_phy_tmp      => p_nh%diag%ddt_vn_phy
!$ACC UPDATE DEVICE ( ddt_exner_phy_tmp,ddt_vn_phy_tmp )

       rho_incr_tmp        => p_nh%diag%rho_incr
       exner_incr_tmp      => p_nh%diag%exner_incr
!$ACC UPDATE DEVICE ( rho_incr_tmp, exner_incr_tmp )

       vn_traj_tmp       => prep_adv%vn_traj
       mass_flx_me_tmp   => prep_adv%mass_flx_me
       mass_flx_ic_tmp   => prep_adv%mass_flx_ic
!$ACC UPDATE DEVICE ( vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp ) IF( lprep_adv )

       vn_ref_tmp          => p_nh%ref%vn_ref
       w_ref_tmp           => p_nh%ref%w_ref
!$ACC UPDATE DEVICE ( vn_ref_tmp, w_ref_tmp )

       grf_bdy_mflx_tmp   => p_nh%diag%grf_bdy_mflx
!$ACC UPDATE DEVICE( grf_bdy_mflx_tmp ) IF( (jg > 1) .AND. (grf_intmethod_e >= 5) .AND. (idiv_method == 1) .AND. (jstep == 0) )

     END SUBROUTINE h2d_solve_nonhydro

     SUBROUTINE d2h_solve_nonhydro( nnew, jstep, jg, idyn_timestep, grf_intmethod_e, idiv_method, lsave_mflx, l_child_vertnest, lprep_adv, p_nh, prep_adv )

       INTEGER, INTENT(IN)       :: nnew, jstep, jg, idyn_timestep, grf_intmethod_e, idiv_method
       LOGICAL, INTENT(IN)       :: lsave_mflx, l_child_vertnest, lprep_adv

       TYPE(t_nh_state),          INTENT(INOUT) :: p_nh
       TYPE(t_prepare_adv),       INTENT(INOUT) :: prep_adv

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp                 ! p_prog  WP
       REAL(wp), DIMENSION(:,:),     POINTER  :: dvn_ie_int_tmp                                                 ! p_diag  WP 2D
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: theta_v_ic_tmp, rho_ic_tmp, dw_int_tmp                         ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: dtheta_v_ic_int_tmp,  grf_bdy_mflx_tmp                         ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: mass_fl_e_tmp,  mflx_ic_int_tmp, exner_pr_tmp                  ! p_diag  WP

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp, w_concorr_c_tmp                             ! p_diag  VP
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: mass_fl_e_sv_tmp, exner_dyn_incr_tmp                           ! p_diag  VP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp                  ! prep_adv WP


! The following code is necessary if the Dycore is to be run in isolation on the GPU
! Update all device output on host: the prognostic variables have shifted from nnow to nnew; diagnostics pointers set above

       exner_tmp           => p_nh%prog(nnew)%exner
       rho_tmp             => p_nh%prog(nnew)%rho
       theta_v_tmp         => p_nh%prog(nnew)%theta_v
       vn_tmp              => p_nh%prog(nnew)%vn
       w_tmp               => p_nh%prog(nnew)%w
!$ACC UPDATE HOST ( exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp )

       vt_tmp              => p_nh%diag%vt
       vn_ie_tmp           => p_nh%diag%vn_ie
       rho_ic_tmp          => p_nh%diag%rho_ic
       theta_v_ic_tmp      => p_nh%diag%theta_v_ic
       exner_pr_tmp        => p_nh%diag%exner_pr
!$ACC UPDATE HOST ( vt_tmp, vn_ie_tmp, rho_ic_tmp, theta_v_ic_tmp, exner_pr_tmp )

       w_concorr_c_tmp     => p_nh%diag%w_concorr_c
       mass_fl_e_tmp       => p_nh%diag%mass_fl_e
       exner_dyn_incr_tmp  => p_nh%diag%exner_dyn_incr
!$ACC UPDATE HOST ( w_concorr_c_tmp, mass_fl_e_tmp, exner_dyn_incr_tmp )

       mass_fl_e_sv_tmp    => p_nh%diag%mass_fl_e_sv
!$ACC UPDATE HOST ( mass_fl_e_sv_tmp ) IF( lsave_mflx )

       dw_int_tmp          => p_nh%diag%dw_int
       mflx_ic_int_tmp     => p_nh%diag%mflx_ic_int
       dtheta_v_ic_int_tmp => p_nh%diag%dtheta_v_ic_int
!$ACC UPDATE HOST ( dw_int_tmp, mflx_ic_int_tmp, dtheta_v_ic_int_tmp ) IF( l_child_vertnest )

      dvn_ie_int_tmp      => p_nh%diag%dvn_ie_int
!$ACC UPDATE HOST ( dvn_ie_int_tmp ) IF( idyn_timestep == 1 .AND. l_child_vertnest)

      grf_bdy_mflx_tmp    => p_nh%diag%grf_bdy_mflx
!$ACC UPDATE HOST ( grf_bdy_mflx_tmp ) IF( (jg > 1) .AND. (grf_intmethod_e >= 5) .AND. (idiv_method == 1) .AND. (jstep == 0) )

      vn_traj_tmp         => prep_adv%vn_traj
      mass_flx_me_tmp     => prep_adv%mass_flx_me
      mass_flx_ic_tmp     => prep_adv%mass_flx_ic
!$ACC UPDATE HOST ( vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp ) IF( lprep_adv )

     END SUBROUTINE d2h_solve_nonhydro


     SUBROUTINE h2d_velocity_tendencies( p_prog, p_diag, z_w_concorr_me, z_kin_hor_e, z_vt_ie )
       TYPE(t_nh_prog), INTENT(INOUT)            :: p_prog
       TYPE(t_nh_diag), INTENT(INOUT)            :: p_diag
       REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie


       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_tmp, w_tmp
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: w_concorr_c_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_adv_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_w_adv_tmp

       vn_tmp              => p_prog%vn
       w_tmp               => p_prog%w
       vt_tmp              => p_diag%vt
       vn_ie_tmp           => p_diag%vn_ie
       w_concorr_c_tmp     => p_diag%w_concorr_c
       ddt_vn_adv_tmp      => p_diag%ddt_vn_adv
       ddt_w_adv_tmp       => p_diag%ddt_w_adv

!$ACC UPDATE DEVICE ( vn_tmp, w_tmp, vt_tmp, vn_ie_tmp, w_concorr_c_tmp, ddt_vn_adv_tmp, ddt_w_adv_tmp )
!$ACC UPDATE DEVICE ( z_w_concorr_me, z_kin_hor_e, z_vt_ie )

     END SUBROUTINE h2d_velocity_tendencies

     SUBROUTINE d2h_velocity_tendencies( istep, ntnd, p_diag, z_w_concorr_me, z_kin_hor_e, z_vt_ie )

       INTEGER, INTENT(IN)                       :: istep, ntnd
       TYPE(t_nh_diag), INTENT(INOUT)            :: p_diag
       REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: w_concorr_c_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_adv_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_w_adv_tmp

       vt_tmp              => p_diag%vt
       vn_ie_tmp           => p_diag%vn_ie
       w_concorr_c_tmp     => p_diag%w_concorr_c
       ddt_vn_adv_tmp      => p_diag%ddt_vn_adv
       ddt_w_adv_tmp       => p_diag%ddt_w_adv

!$ACC UPDATE HOST( z_kin_hor_e, z_vt_ie, z_w_concorr_me, vt_tmp, vn_ie_tmp, w_concorr_c_tmp ), IF( istep==1 )
!$ACC UPDATE HOST( ddt_vn_adv_tmp, ddt_w_adv_tmp(:,:,:,ntnd) )

     END SUBROUTINE d2h_velocity_tendencies

#endif

END MODULE mo_nonhydro_gpu_types
