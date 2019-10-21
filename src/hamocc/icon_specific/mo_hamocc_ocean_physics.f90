#ifndef __NO_ICON_OCEAN__
#include "hamocc_omp_definitions.inc"     
#include "icon_definitions.inc" 
   MODULE mo_hamocc_ocean_physics

    USE mo_kind,                         ONLY: wp
    USE mo_exception,                    ONLY: finish
    USE mo_param1_bgc,                   ONLY: n_bgctra
    USE mo_parallel_config,              ONLY: nproma
    USE mo_model_domain,                 ONLY: t_patch_3D, t_patch
    USE mo_grid_subset,                  ONLY: get_index_range, t_subset_range
    USE mo_timer,                        ONLY: timer_start, timer_stop, &
    &                                          timer_tracer_ab, timer_bgc_tot,&
    &                                          timers_level
    USE mo_ocean_tracer,                 ONLY: advect_ocean_tracers
    USE mo_ocean_tracer_GMRedi,          ONLY: advect_ocean_tracers_GMRedi
    USE mo_ocean_types,                  ONLY: t_hydro_ocean_state, &
    &                                          t_operator_coeff
    USE mo_ocean_tracer_transport_types, ONLY: t_tracer_collection, t_ocean_transport_state
    USE mo_ocean_surface_types,          ONLY: t_ocean_surface, t_atmos_for_ocean
    USE mo_sea_ice_types,                ONLY: t_sea_ice
    USE mo_run_config,                   ONLY: ltimer
    USE mo_ocean_nml,                    ONLY: Cartesian_Mixing, GMRedi_configuration
    USE mo_ocean_physics,                ONLY: update_ho_params
    USE mo_hamocc_types,                 ONLY: t_hamocc_prog, t_hamocc_state
    USE mo_ocean_physics_types,          ONLY: t_ho_params 
    USE mo_bgc_icon_comm,                ONLY: hamocc_state
    USE mo_dynamics_config,              ONLY: nold, nnew 
    USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
      & t_hamocc_ocean_state
    USE mo_hamocc_diagnostics, ONLY: get_monitoring 
    USE mo_bgc_bcond,              ONLY: ext_data_bgc, update_bgc_bcond
    USE mtime,                     ONLY: datetime
    
    IMPLICIT NONE
    PRIVATE

    PUBLIC:: tracer_biochemistry_transport

    CONTAINS

  SUBROUTINE tracer_biochemistry_transport(hamocc_ocean_state, operators_coefficients, current_time)

    TYPE(t_hamocc_ocean_state), TARGET               :: hamocc_ocean_state
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_to_hamocc_state), POINTER           :: ocean_to_hamocc_state
    TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state
    TYPE(t_patch_3d ),POINTER                        :: patch_3d
    TYPE(t_ocean_transport_state), POINTER           :: transport_state

    INTEGER :: i
    
    transport_state => hamocc_ocean_state%ocean_transport_state
    patch_3d => transport_state%patch_3d
    ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
    hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state

    CALL dilute_hamocc_tracers(patch_3d, ocean_to_hamocc_state%top_dilution_coeff, hamocc_state%p_prog(nold(1)))
    
    CALL update_bgc_bcond( patch_3d, ext_data_bgc,  current_time)
       
    !------------------------------------------------------------------------
    ! call HAMOCC
    if(ltimer) call timer_start(timer_bgc_tot)
    CALL bgc_icon(patch_3d, hamocc_ocean_state)
    if(ltimer) call timer_stop(timer_bgc_tot)

    !------------------------------------------------------------------------
    ! transport tracers and diffuse them
    ! fill diffusion coefficients
    old_tracer_collection => hamocc_state%p_prog(nold(1))%tracer_collection
    new_tracer_collection => hamocc_state%p_prog(nnew(1))%tracer_collection


    DO i = 1, old_tracer_collection%no_of_tracers
        old_tracer_collection%tracer(i)%hor_diffusion_coeff => ocean_to_hamocc_state%hor_diffusion_coeff
        old_tracer_collection%tracer(i)%ver_diffusion_coeff => ocean_to_hamocc_state%ver_diffusion_coeff
    ENDDO

    start_timer(timer_tracer_ab,1)

    IF (GMRedi_configuration==Cartesian_Mixing) THEN
        CALL advect_ocean_tracers(old_tracer_collection, new_tracer_collection, transport_state, operators_coefficients)
    ELSE
      CALL finish("tracer_biochemistry_transport", "GMRedi_configuration is not activated")
    ENDIF
        ! the GMRedi will be treated in the sequential case; for the moment we will skip it
!       ELSE
!         CALL  advect_ocean_tracers_GMRedi(old_tracer_collection, new_tracer_collection, &
!           &  ocean_state, transport_state, p_phys_param, operators_coefficients)
!       ENDIF

     stop_timer(timer_tracer_ab,1)
     !------------------------------------------------------------------------
    
     CALL get_monitoring( hamocc_state, hamocc_state%p_prog(nnew(1))%tracer, ocean_to_hamocc_state%h_new, patch_3d)
    !------------------------------------------------------------------------


    END SUBROUTINE tracer_biochemistry_transport



    SUBROUTINE DILUTE_HAMOCC_TRACERS(p_patch_3D, top_dilution_coeff, hamocc_state_prog)
    ! HAMOCC tracer dilution using old and new h

    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_hamocc_prog), INTENT(inout)       :: hamocc_state_prog
    REAL(wp), INTENT(in)                     :: top_dilution_coeff(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    ! Local variables
    INTEGER:: jb, jc, i_bgc_tra, i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch

    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all

    DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
            IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN
                DO i_bgc_tra = 1, n_bgctra
                     hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra)  = hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra) &
    &                                                               * top_dilution_coeff(jc,jb)
                ENDDO
            endif
        ENDDO
    ENDDO

   END SUBROUTINE DILUTE_HAMOCC_TRACERS
   END MODULE
#endif   
