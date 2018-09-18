!>
!! Provide an implementation of the ocean tracer nuding functionality.
!!
!!
!! @author Helmuth Haak, MPI
!!
!! @par Revision History
!!  Original version by Helmuth Haak, MPI-M (2018)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_nudging
!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,              &
    & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
    & type_3dimrelax_temp, para_3dimrelax_temp,                           &
    & type_3dimrelax_salt, para_3dimrelax_salt!,                           &
 !   & iswm_oce,                 use_none,                                 &
 !   & flux_calculation_horz, flux_calculation_vert, miura_order1,         &
 !   & l_with_vert_tracer_diffusion, l_with_vert_tracer_advection,         &
 !   & GMRedi_configuration,GMRedi_combined,  GM_only,Redi_only ,          &
 !   & Cartesian_Mixing, tracer_threshold_min, tracer_threshold_max,       &
 !   & namelist_tracer_name, tracer_update_mode, use_none, nbgcadv,        &
 !   & GMREDI_COMBINED_DIAGNOSTIC,GM_INDIVIDUAL_DIAGNOSTIC,REDI_INDIVIDUAL_DIAGNOSTIC
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer 
  USE mo_ocean_nudging_types,       ONLY: t_ocean_nudge
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
!  USE mo_ocean_boundcond,           ONLY: top_bound_cond_tracer
!  USE mo_ocean_physics_types,       ONLY: t_ho_params
!  USE mo_ocean_surface_types,       ONLY: t_ocean_surface
!  USE mo_ocean_diffusion,           ONLY: tracer_diffusion_vertical_implicit, tracer_diffusion_vert_explicit,tracer_diffusion_horz
!  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
!  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
!  USE mo_sync,                      ONLY: sync_c, sync_e, sync_patch_array
!  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_dif_vert, timer_extra30
!  USE mo_statistics,                ONLY: global_minmaxmean, print_value_location
!  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
!  USE mo_ocean_GM_Redi,             ONLY: calc_ocean_physics, prepare_ocean_physics
!  USE mo_ocean_math_operators,      ONLY: div_oce_3d, verticalDiv_scalar_onFullLevels! !verticalDiv_scalar_midlevel
!  USE mo_scalar_product,            ONLY: map_edges2edges_viacell_3d_const_z
!  USE mo_physical_constants,        ONLY: clw, rho_ref,sitodbar
!  USE mo_ocean_thermodyn,           ONLY: calculate_density, calc_potential_density
!  USE mo_ocean_pp_scheme,           ONLY: calculate_rho4GMRedi


  IMPLICIT NONE

  ! required for reading netcdf files
  INCLUDE 'netcdf.inc'


  TYPE(t_ocean_nudge)  :: ocean_nudge

  PRIVATE

  CHARACTER(LEN=12)           :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  ! Public interface
  PUBLIC :: nudge_ocean_tracers, ocean_nudge


CONTAINS

      !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE nudge saLINITY AND TEMERATURE
  !!
  !! @par Revision History
  !! Developed  by  Helmuth Haak, MPI-M (2018).
  !!
!<Optimize:inUse>
  SUBROUTINE nudge_ocean_tracers(patch_3d, p_os)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)      :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: p_os

    !Local variables

    REAL(wp) :: z_relax
    REAL(wp) :: z_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    patch_2d  => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%ALL

    ! Final step: 3-dim temperature relaxation
    !  - strict time constant, i.e. independent of layer thickness
    !  - additional forcing Term F_T = -1/tau(T-T*) [ K/s ]
    !    when using the sign convention
    !      dT/dt = Operators + F_T
    !    i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature if it is warmer than relaxation data)
    !  - discretized:
    !    tracer = tracer - 1/(para_3dimRelax_Temp[months]) * (tracer(1)-data_3dimRelax_Temp)
    IF (no_tracer>=1 .AND. type_3dimrelax_temp >0) THEN

      ! calculate relaxation term
      z_relax = 1.0_wp/(para_3dimrelax_temp*2.592e6_wp)
!      ocean_nudge%forc_3dimrelax_temp(:,:,:) = -z_relax * ocean_nudge%relax_3dim_coefficient(:,:,:) &
!        & * ( p_os%p_prog(nnew(1))%tracer(:,:,:,1) - ocean_nudge%data_3dimrelax_temp(:,:,:))


      ocean_nudge%forc_3dimrelax_temp(:,:,:) = z_relax * dtime &
         * (ocean_nudge%data_3dimrelax_temp(:,:,:) - p_os%p_prog(nnew(1))%tracer(:,:,:,1) )

      ! add relaxation term to new temperature
      p_os%p_prog(nnew(1))%tracer(:,:,:,1) = p_os%p_prog(nnew(1))%tracer(:,:,:,1) + &
                                             ocean_nudge%forc_3dimRelax_temp(:,:,:)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('3d_relax: tracer T forc', ocean_nudge%forc_3dimRelax_Temp, str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('3d_relax: tracer T data', ocean_nudge%data_3dimRelax_Temp, str_module,idt_src, in_subset=cells_in_domain)      
      z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,1)
      CALL dbg_print('3d_relax: tracer T trac', z_c, str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF

    ! Final step: 3-dim salinity relaxation
    !  - additional forcing Term F_S = -1/tau(S-S*) [ psu/s ]
    !    when using the sign convention
    !      dS/dt = Operators + F_S
    !    i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity if it is larger than relaxation data)
    !  - discretized:
    !    tracer = tracer - 1/(para_3dimRelax_salt[months]) * (tracer(1)-data_3dimRelax_salt)
    IF (no_tracer>=2 .AND. type_3dimrelax_salt >0) THEN

      ! calculate relaxation term
      z_relax = 1.0_wp/(para_3dimrelax_salt*2.592e6_wp)
!      ocean_nudge%forc_3dimrelax_salt(:,:,:) = -z_relax * ocean_nudge%relax_3dim_coefficient(:,:,:) &
!        & ( p_os%p_prog(nnew(1))%tracer(:,:,:,2) -       &
!        & ocean_nudge%forc_3dimrelax_salt(:,:,:))

      ocean_nudge%forc_3dimrelax_salt(:,:,:) = z_relax * dtime &
         * (ocean_nudge%data_3dimrelax_salt(:,:,:) - p_os%p_prog(nnew(1))%tracer(:,:,:,2) )

      ! add relaxation term to new temperature
      p_os%p_prog(nnew(1))%tracer(:,:,:,2) = p_os%p_prog(nnew(1))%tracer(:,:,:,2) + &
                                             ocean_nudge%forc_3dimRelax_salt(:,:,:)

      ! add relaxation term to new salinity

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('3d_relax: tracer S forc'  ,ocean_nudge%forc_3dimrelax_salt,str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('3d_relax: tracer S data'  ,ocean_nudge%data_3dimrelax_salt,str_module,idt_src, in_subset=cells_in_domain)
      z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,2)
      CALL dbg_print('3d_relax: tracer S trac'  ,z_c                           ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF    

  END SUBROUTINE nudge_ocean_tracers
  !-------------------------------------------------------------------------


END MODULE mo_ocean_nudging


