!>
!! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!!
!! @author Leonidas Linardakis, MPIM
!!
!
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!----------------------------
#include "icon_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_to_hamocc_interface
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length, success
  USE mo_parallel_config,        ONLY: nproma
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_coupling_config,        ONLY: is_coupled_run
  USE mo_ocean_nml,              ONLY: n_zlev, no_tracer, lhamocc, &
    &  Cartesian_Mixing, GMRedi_configuration
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_timer,                  ONLY: timer_start, timer_stop, &
    & timer_bgc_inv, timer_bgc_tot, ltimer
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, &
    & t_operator_coeff
  USE mo_ocean_tracer_transport_types, ONLY: t_tracer_collection, t_ocean_transport_state
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_physics_types,    ONLY: t_ho_params  
  USE mo_master_config,          ONLY: isRestart
  USE mo_hamocc_ocean_physics,   ONLY: tracer_biochemistry_transport
  USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_ocean_transport_state, t_hamocc_ocean_state
  USE mtime,                     ONLY: datetime    
  USE mo_construct_icon_hamocc,  ONLY: construct_icon_hamocc, destruct_icon_hamocc, init_icon_hamocc
   USE mo_ext_data_types,        ONLY: t_external_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: ocean_to_hamocc_construct, ocean_to_hamocc_init, ocean_to_hamocc_end, ocean_to_hamocc_interface
  !-------------------------------------------------------------------------

  TYPE(t_hamocc_ocean_state), TARGET   :: hamocc_ocean_state
  
CONTAINS
!-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE ocean_to_hamocc_construct(patch_3d, ext_data)

    TYPE(t_patch_3d ), INTENT(in)     :: patch_3d
    TYPE(t_external_data), TARGET, INTENT(inout) :: ext_data
   
    IF(.not. lhamocc) return
   
    CALL construct_icon_hamocc(patch_3d, ext_data)

  END SUBROUTINE ocean_to_hamocc_construct
  !--------------------------------------------datetime-----------------------------

!-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE ocean_to_hamocc_init(patch_3d, ocean_state, &
    & p_as, sea_ice, p_oce_sfc, p_phys_param)

    TYPE(t_patch_3d ), INTENT(in) , TARGET :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_atmos_for_ocean)               :: p_as
    TYPE (t_sea_ice)                      :: sea_ice
    TYPE(t_ocean_surface)                 :: p_oce_sfc
     TYPE (t_ho_params)                   :: p_phys_param
 
    TYPE(t_ocean_transport_state), TARGET   :: transport_state
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: alloc_cell_blocks, nblks_e

  
    IF(.not. lhamocc) return
   
    patch_2d => patch_3d%p_patch_2d(1)
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
   !----------------------------------------------
    ! dummy allocation for using the same communication patterns
    ALLOCATE(transport_state%h_new(nproma,alloc_cell_blocks), &
             transport_state%h_old(nproma,alloc_cell_blocks), &
             transport_state%vn(nproma,n_zlev,nblks_e),      &
             transport_state%w(nproma,n_zlev+1,alloc_cell_blocks))
    transport_state%h_new = 0.0_wp
    transport_state%h_old = 0.0_wp
    transport_state%vn    = 0.0_wp
    transport_state%w     = 0.0_wp
    transport_state%patch_3d => patch_3d
    !----------------------------------------------
    

    CALL fill_ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param)
      
    CALL init_icon_hamocc (hamocc_ocean_state)  
    
    DEALLOCATE(transport_state%h_new, &
             transport_state%h_old,   &
             transport_state%vn,      &
             transport_state%w)

  END SUBROUTINE ocean_to_hamocc_init
  !--------------------------------------------datetime-----------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE ocean_to_hamocc_end()

    IF(.not. lhamocc) return

    CALL destruct_icon_hamocc()
   
  END SUBROUTINE ocean_to_hamocc_end
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param, operators_coefficients, current_time)
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in)    :: ocean_state
    TYPE(t_ocean_transport_state), TARGET            :: transport_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
   
    IF(.not. lhamocc) return
    
    CALL fill_ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param)
    
    !------------------------------------------------------------------------
    CALL tracer_biochemistry_transport(hamocc_ocean_state, operators_coefficients, current_time)

  END SUBROUTINE ocean_to_hamocc_interface
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param)
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in)    :: ocean_state
    TYPE(t_ocean_transport_state), TARGET            :: transport_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
   
    TYPE(t_ocean_to_hamocc_state), POINTER           :: ocean_to_hamocc_state
    TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state
 
    ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
    hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state
    hamocc_ocean_state%ocean_transport_state => transport_state 
    hamocc_ocean_state%patch_3d => transport_state%patch_3d
    
    ocean_to_hamocc_state%top_dilution_coeff => p_oce_sfc%top_dilution_coeff
    ocean_to_hamocc_state%h_old              => transport_state%h_old
    ocean_to_hamocc_state%h_new              => transport_state%h_new
    ocean_to_hamocc_state%ice_concentration_sum => sea_ice%concSum    
    ocean_to_hamocc_state%temperature        => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    ocean_to_hamocc_state%salinity           => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
    ocean_to_hamocc_state%hor_diffusion_coeff => p_phys_param%TracerDiffusion_coeff(:,:,:,2)
    ocean_to_hamocc_state%ver_diffusion_coeff => p_phys_param%a_tracer_v(:,:,:,2)
    ocean_to_hamocc_state%short_wave_flux    => p_as%fswr  ! p_oce_sfc%HeatFlux_ShortWave
    ocean_to_hamocc_state%wind10m            => p_as%fu10
    ocean_to_hamocc_state%co2_mixing_ratio   => p_as%co2 
   
    hamocc_to_ocean_state%co2_flux           => p_as%co2flx
    
  END SUBROUTINE fill_ocean_to_hamocc_interface
     !------------------------------------------------------------------------

END MODULE mo_ocean_to_hamocc_interface
