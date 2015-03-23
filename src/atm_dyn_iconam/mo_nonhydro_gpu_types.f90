#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>

!! Type definition for the dynamical core of ICONAM.
!!
!! @author Almut Gassmann (MPI-M)
!! @author Daniel Reinert (DWD)
!! @author Guenther Zaengl (DWD)
!!
!! @par Revision History
!! Initial release by Daniel Reinert, DWD (2012-02-07)
!! - Moved here from mo_nonhydro_state to avoid circular dependencies
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
  USE mo_mpi,                  ONLY: my_process_is_work, i_am_accel_node
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d
  USE mo_model_domain,         ONLY: t_patch, p_patch
  USE mo_nonhydro_state,       ONLY: p_nh_state
  USE mo_nonhydro_types,       ONLY: t_nh_diag

  IMPLICIT NONE

  PRIVATE 

#if defined( _OPENACC )

  PUBLIC :: init_gpu_variables, finalize_gpu_variables
  PUBLIC :: save_convenience_pointers
  PUBLIC :: refresh_convenience_pointers

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
    TYPE ( t_ptr_2d3d ), ALLOCATABLE :: ddt_vn_adv_ptr(:)
    TYPE ( t_ptr_2d3d ), ALLOCATABLE :: ddt_w_adv_ptr(:)
  ENDTYPE my_nh_state

  TYPE ( my_patches ), ALLOCATABLE :: save_patches(:)

  TYPE ( my_nh_state ), ALLOCATABLE :: save_nh_state(:)

  CHARACTER(len=*), PARAMETER :: version = &
    & '$Id: mo_nonhydro_gpu_types.f90 16760 2014-04-01 09:04:22Z wsawyer $'

CONTAINS

     SUBROUTINE init_gpu_variables( )

       i_am_accel_node   = my_process_is_work()

     END SUBROUTINE init_gpu_variables


     SUBROUTINE finalize_gpu_variables( )

       i_am_accel_node = .FALSE.

     END SUBROUTINE finalize_gpu_variables


     SUBROUTINE save_convenience_pointers( )

       INTEGER :: jg, nl, num_levs

!
! Save then nullify convenience pointers on the Host, in order to perform a deep copy to the device
!

       ALLOCATE( save_patches( SIZE(p_patch) ) )
       DO jg = 1, SIZE(p_patch)
         NULLIFY( p_patch(jg)%cells%decomp_info%halo_level )
         NULLIFY( p_patch(jg)%edges%decomp_info%halo_level )
         NULLIFY( p_patch(jg)%verts%decomp_info%halo_level )

         save_patches(jg)%all_cells_patch =>  p_patch(jg)%cells%all%patch
         save_patches(jg)%all_edges_patch =>  p_patch(jg)%edges%all%patch
         save_patches(jg)%all_verts_patch =>  p_patch(jg)%verts%all%patch

         NULLIFY( p_patch(jg)%cells%all%patch )
         NULLIFY( p_patch(jg)%edges%all%patch )
         NULLIFY( p_patch(jg)%verts%all%patch )

         save_patches(jg)%owned_cells_patch =>  p_patch(jg)%cells%owned%patch
         save_patches(jg)%owned_edges_patch =>  p_patch(jg)%edges%owned%patch
         save_patches(jg)%owned_verts_patch =>  p_patch(jg)%verts%owned%patch

         NULLIFY( p_patch(jg)%cells%owned%patch )
         NULLIFY( p_patch(jg)%edges%owned%patch )
         NULLIFY( p_patch(jg)%verts%owned%patch )

         save_patches(jg)%in_domain_cells_patch =>  p_patch(jg)%cells%in_domain%patch
         save_patches(jg)%in_domain_edges_patch =>  p_patch(jg)%edges%in_domain%patch
         save_patches(jg)%in_domain_verts_patch =>  p_patch(jg)%verts%in_domain%patch

         NULLIFY( p_patch(jg)%cells%in_domain%patch )
         NULLIFY( p_patch(jg)%edges%in_domain%patch )
         NULLIFY( p_patch(jg)%verts%in_domain%patch )

         save_patches(jg)%gradiscalculable_edges_patch =>  p_patch(jg)%edges%gradiscalculable%patch

         NULLIFY( p_patch(jg)%edges%gradiscalculable%patch )

         save_patches(jg)%not_owned_cells_patch =>  p_patch(jg)%cells%not_owned%patch
         save_patches(jg)%not_owned_edges_patch =>  p_patch(jg)%edges%not_owned%patch
         save_patches(jg)%not_owned_verts_patch =>  p_patch(jg)%verts%not_owned%patch

         NULLIFY( p_patch(jg)%cells%not_owned%patch )
         NULLIFY( p_patch(jg)%edges%not_owned%patch )
         NULLIFY( p_patch(jg)%verts%not_owned%patch )

         save_patches(jg)%not_in_domain_cells_patch =>  p_patch(jg)%cells%not_in_domain%patch
         save_patches(jg)%not_in_domain_edges_patch =>  p_patch(jg)%edges%not_in_domain%patch
         save_patches(jg)%not_in_domain_verts_patch =>  p_patch(jg)%verts%not_in_domain%patch

         NULLIFY( p_patch(jg)%cells%not_in_domain%patch )
         NULLIFY( p_patch(jg)%edges%not_in_domain%patch )
         NULLIFY( p_patch(jg)%verts%not_in_domain%patch )

         save_patches(jg)%one_edge_in_domain_cells_patch =>  p_patch(jg)%cells%one_edge_in_domain%patch

         NULLIFY( p_patch(jg)%cells%one_edge_in_domain%patch )

       ENDDO

       ALLOCATE( save_nh_state( SIZE(p_nh_state) ) )

       DO jg = 1, SIZE(p_nh_state)
! ddt_vn_adv_ptr
         num_levs = SIZE(p_nh_state(jg)%diag%ddt_vn_adv_ptr)
         ALLOCATE( save_nh_state(jg)%ddt_vn_adv_ptr( num_levs ) )
         DO nl = 1, num_levs 
           save_nh_state(jg)%ddt_vn_adv_ptr( nl )%p_2d => p_nh_state(jg)%diag%ddt_vn_adv_ptr( nl )%p_2d 
           save_nh_state(jg)%ddt_vn_adv_ptr( nl )%p_3d => p_nh_state(jg)%diag%ddt_vn_adv_ptr( nl )%p_3d 
         ENDDO
         DEALLOCATE( p_nh_state(jg)%diag%ddt_vn_adv_ptr )  ! p_2d not used?

! ddt_w_adv_ptr
         num_levs = SIZE(p_nh_state(jg)%diag%ddt_w_adv_ptr)
         ALLOCATE( save_nh_state(jg)%ddt_w_adv_ptr( num_levs ) )
         DO nl = 1, num_levs 
           save_nh_state(jg)%ddt_w_adv_ptr( nl )%p_2d => p_nh_state(jg)%diag%ddt_w_adv_ptr( nl )%p_2d 
           save_nh_state(jg)%ddt_w_adv_ptr( nl )%p_3d => p_nh_state(jg)%diag%ddt_w_adv_ptr( nl )%p_3d 
         ENDDO
         DEALLOCATE( p_nh_state(jg)%diag%ddt_w_adv_ptr )  ! p_2d not used?

       ENDDO

     END SUBROUTINE save_convenience_pointers

     SUBROUTINE refresh_convenience_pointers( )

       INTEGER :: jg, nl, num_levs

!
! Save then nullify convenience pointers on the Host, in order to perform a deep copy to the device
!

       DO jg = 1, SIZE(p_patch)

         p_patch(jg)%cells%all%patch =>                   save_patches(jg)%all_cells_patch
         p_patch(jg)%edges%all%patch =>                   save_patches(jg)%all_edges_patch
         p_patch(jg)%verts%all%patch =>                   save_patches(jg)%all_verts_patch

         p_patch(jg)%cells%owned%patch =>                 save_patches(jg)%owned_cells_patch
         p_patch(jg)%edges%owned%patch =>                 save_patches(jg)%owned_edges_patch
         p_patch(jg)%verts%owned%patch =>                 save_patches(jg)%owned_verts_patch

         p_patch(jg)%cells%in_domain%patch =>             save_patches(jg)%in_domain_cells_patch
         p_patch(jg)%edges%in_domain%patch =>             save_patches(jg)%in_domain_edges_patch
         p_patch(jg)%verts%in_domain%patch =>             save_patches(jg)%in_domain_verts_patch

         p_patch(jg)%edges%gradiscalculable%patch =>       save_patches(jg)%gradiscalculable_edges_patch

         p_patch(jg)%cells%not_owned%patch =>             save_patches(jg)%not_owned_cells_patch
         p_patch(jg)%edges%not_owned%patch =>             save_patches(jg)%not_owned_edges_patch
         p_patch(jg)%verts%not_owned%patch =>             save_patches(jg)%not_owned_verts_patch

         p_patch(jg)%cells%not_in_domain%patch =>         save_patches(jg)%not_in_domain_cells_patch
         p_patch(jg)%edges%not_in_domain%patch =>         save_patches(jg)%not_in_domain_edges_patch
         p_patch(jg)%verts%not_in_domain%patch =>         save_patches(jg)%not_in_domain_verts_patch

         p_patch(jg)%cells%one_edge_in_domain%patch =>    save_patches(jg)%one_edge_in_domain_cells_patch

       ENDDO


       DEALLOCATE( save_patches )

       DO jg = 1, SIZE(save_nh_state)
! ddt_vn_adv_ptr
         num_levs = SIZE( save_nh_state(jg)%ddt_vn_adv_ptr )
         ALLOCATE( p_nh_state(jg)%diag%ddt_vn_adv_ptr( num_levs ) )  ! p_2d not used?
         DO nl = 1, num_levs
           p_nh_state(jg)%diag%ddt_vn_adv_ptr( nl )%p_2d  => save_nh_state(jg)%ddt_vn_adv_ptr( nl )%p_2d 
           p_nh_state(jg)%diag%ddt_vn_adv_ptr( nl )%p_3d  => save_nh_state(jg)%ddt_vn_adv_ptr( nl )%p_3d 
         ENDDO
         DEALLOCATE( save_nh_state(jg)%ddt_vn_adv_ptr )

! ddt_w_adv_ptr
         num_levs = SIZE( save_nh_state(jg)%ddt_w_adv_ptr )
         ALLOCATE( p_nh_state(jg)%diag%ddt_w_adv_ptr( num_levs ) )  ! p_2d not used?
         DO nl = 1, num_levs
           p_nh_state(jg)%diag%ddt_w_adv_ptr( nl )%p_2d  => save_nh_state(jg)%ddt_w_adv_ptr( nl )%p_2d 
           p_nh_state(jg)%diag%ddt_w_adv_ptr( nl )%p_3d  => save_nh_state(jg)%ddt_w_adv_ptr( nl )%p_3d 
         ENDDO
         DEALLOCATE( save_nh_state(jg)%ddt_w_adv_ptr )

       ENDDO

       DEALLOCATE( save_nh_state )


     END SUBROUTINE refresh_convenience_pointers

#endif

END MODULE mo_nonhydro_gpu_types
