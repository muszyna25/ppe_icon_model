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

#endif

END MODULE mo_nonhydro_gpu_types
