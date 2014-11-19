!>
!! Contains basic diagnostics for ICON ocean model.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/02)
!!  Extended   by Stephan Lorenz,   MPI-M (2012)
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
!----------------------------
MODULE mo_ocean_diagnostics
  USE mo_kind,               ONLY: wp, dp, i8
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_grid_tools,         ONLY: get_oriented_edges_from_global_vertices, check_global_indexes
  USE mo_mpi,                ONLY: my_process_is_stdio, p_field_sum, get_my_mpi_work_id, &
    & p_comm_work_test, p_comm_work, p_io, p_bcast
  USE mo_sync,               ONLY: global_sum_array, disable_sync_checks, enable_sync_checks, &
    &                              sync_c, sync_e, sync_patch_array
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, cvec2gvec
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_math_constants,     ONLY: rad2deg, dbl_eps
  USE mo_impl_constants,     ONLY: sea_boundary,sea, &
    & min_rlcell, min_rledge, min_rlcell, &
    & max_char_length, min_dolic
  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer, &
    & gibraltar, &
    & denmark_strait, &
    & drake_passage, &
    & indonesian_throughflow,&
    & scotland_iceland, &
    & mozambique, &
    & framStrait, &
    & beringStrait, &
    & barentsOpening, &
    & agulhas, &
    & agulhas_long, &
    & agulhas_longer, &
    & ab_const, ab_beta, ab_gam, iswm_oce, discretization_scheme
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_physical_constants, ONLY: grav, rho_ref, rhos, rhoi,sice
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d,t_patch_vert, t_grid_edges
  USE mo_ocean_types,          ONLY: t_hydro_ocean_state, t_hydro_ocean_diag,&
    &                              t_ocean_regions, t_ocean_region_volumes, t_ocean_region_areas
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_sea_ice_types,      ONLY: t_sfc_flx, t_sea_ice
  USE mo_datetime,           ONLY: t_datetime, datetime_to_string, date_len
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_scalar_product,     ONLY: map_edges2cell_3d
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_util_file,          ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_statistics,         ONLY: subset_sum
  USE mo_fortran_tools,      ONLY: assign_if_present
  
  IMPLICIT NONE
  
  !PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceDiag     '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug
  
  INTEGER :: diag_unit = -1 ! file handle for the global timeseries output
  INTEGER :: moc_unit  = -1 ! file handle for the global timeseries output
  CHARACTER(LEN=max_char_length) :: diag_fname, moc_fname
  INTEGER, PARAMETER :: linecharacters  = 2048
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics
  PUBLIC :: construct_oce_diagnostics
  PUBLIC :: destruct_oce_diagnostics
  PUBLIC :: t_oce_monitor
  PUBLIC :: t_oce_timeseries
  PUBLIC :: calc_moc
  PUBLIC :: calc_psi
  PUBLIC :: calc_psi_vn
  
  TYPE t_oce_monitor
    REAL(wp) :: volume
    REAL(wp) :: kin_energy
    REAL(wp) :: pot_energy
    REAL(wp) :: total_energy
    REAL(wp) :: total_salt
    REAL(wp) :: vorticity
    REAL(wp) :: enstrophy
    REAL(wp) :: potential_enstrophy
    REAL(wp) :: absolute_vertical_velocity
    REAL(wp) :: HeatFlux_ShortWave
    REAL(wp) :: HeatFlux_LongWave
    REAL(wp) :: HeatFlux_Sensible
    REAL(wp) :: HeatFlux_Latent
    REAL(wp) :: HeatFlux_Total
    REAL(wp) :: FrshFlux_Precipitation
    REAL(wp) :: FrshFlux_SnowFall
    REAL(wp) :: FrshFlux_Evaporation
    REAL(wp) :: FrshFlux_Runoff
    REAL(wp) :: FrshFlux_TotalSalt
    REAL(wp) :: FrshFlux_TotalOcean
    REAL(wp) :: FrshFlux_TotalIce
    REAL(wp) :: FrshFlux_VolumeIce
    REAL(wp) :: FrshFlux_VolumeTotal
    REAL(wp) :: HeatFlux_Relax
    REAL(wp) :: FrshFlux_Relax
    REAL(wp) :: TempFlux_Relax
    REAL(wp) :: SaltFlux_Relax
    
    REAL(wp) :: ice_volume_nh !                                                           [km3]
    REAL(wp) :: ice_volume_sh !                                                           [km3]
    REAL(wp) :: ice_extent_nh !                                                           [km2]
    REAL(wp) :: ice_extent_sh !                                                           [km2]
    ! through flows {{{
    REAL(wp) :: gibraltar     ! though flow                                               [Sv]
    REAL(wp) :: denmark_strait! though flow                                               [Sv]
    REAL(wp) :: drake_passage ! though flow                                               [Sv]
    REAL(wp) :: indonesian_throughflow !                                                  [Sv]
    REAL(wp) :: scotland_iceland !                                                        [Sv]
    REAL(wp) :: mozambique
    REAL(wp) :: framStrait
    REAL(wp) :: beringStrait
    REAL(wp) :: barentsOpening
    REAL(wp) :: agulhas
    REAL(wp) :: agulhas_long
    REAL(wp) :: agulhas_longer
    ! }}}
    REAL(wp) :: t_mean_na_200m !                                                        [degC]
    REAL(wp) :: t_mean_na_800m !                                                        [degC]
    REAL(wp) :: ice_ocean_heat_budget
    REAL(wp) :: ice_ocean_salinity_budget
    REAL(wp) :: ice_ocean_volume_budget
    REAL(wp), ALLOCATABLE :: tracer_content(:)
    
  END TYPE t_oce_monitor
  
  TYPE t_oce_timeseries
    
    TYPE(t_oce_monitor), ALLOCATABLE :: oce_diagnostics(:)    ! time array of diagnostic values
    CHARACTER(LEN=40), DIMENSION(50)  :: names = (/ &
      & "total_volume                            ", &
      & "kin_energy                              ", &
      & "pot_energy                              ", &
      & "total_energy                            ", &
      & "total_salt                              ", &
      & "vorticity                               ", &
      & "enstrophy                               ", &
      & "potential_enstrophy                     ", &
      & "absolute_vertical_velocity              ", &
      & "HeatFlux_ShortWave                      ", &
      & "HeatFlux_LongWave                       ", &
      & "HeatFlux_Sensible                       ", &
      & "HeatFlux_Latent                         ", &
      & "HeatFlux_Total                          ", &
      & "FrshFlux_Precipitation                  ", &
      & "FrshFlux_SnowFall                       ", &
      & "FrshFlux_Evaporation                    ", &
      & "FrshFlux_Runoff                         ", &
      & "FrshFlux_TotalSalt                      ", &
      & "FrshFlux_TotalOcean                     ", &
      & "FrshFlux_TotalIce                       ", &
      & "FrshFlux_VolumeIce                      ", &
      & "FrshFlux_VolumeTotal                    ", &
      & "HeatFlux_Relax                          ", &
      & "FrshFlux_Relax                          ", &
      & "TempFlux_Relax                          ", &
      & "SaltFlux_Relax                          ", &
      & "ice_volume_nh                           ", &
      & "ice_volume_sh                           ", &
      & "ice_extent_nh                           ", &
      & "ice_extent_sh                           ", &
      & "gibraltar                               ", &
      & "denmark_strait                          ", &
      & "drake_passage                           ", &
      & "indonesian_throughflow                  ", &
      & "scotland_iceland                        ", &
      & "mozambique                              ", &
      & "framStrait                              ", &
      & "beringStrait                            ", &
      & "barentsOpening                          ", &
      & "agulhas                                 ", &
      & "agulhas_long                            ", &
      & "agulhas_longer                          ", &
      & "t_mean_NA_200m                          ", &
      & "t_mean_NA_800m                          ", &
      & "ice_ocean_heat_budget                   ", &
      & "ice_ocean_salinity_budget               ", &
      & "ice_ocean_volume_budget                 ", &
      & "total_temperature                       ", &
      & "total_salinity                          "/)
    
  END TYPE t_oce_timeseries
  
  TYPE t_oce_section
    TYPE(t_subset_indexed) :: subset
    REAL(wp), POINTER :: orientation(:)
  END TYPE t_oce_section
  
  INTEGER, PARAMETER  :: oce_section_count = 12
  PRIVATE             :: oce_section_count
  TYPE(t_oce_section) :: oce_sections(oce_section_count)
  PRIVATE             :: oce_sections
  
  TYPE(t_ocean_region_volumes),SAVE :: ocean_region_volumes
  PRIVATE                           :: ocean_region_volumes
  TYPE(t_ocean_region_areas),SAVE   :: ocean_region_areas
  PRIVATE                           :: ocean_region_areas

  TYPE(t_oce_timeseries),POINTER :: oce_ts
  
CONTAINS

  !-------------------------------------------------------------------------
  !  The constructor of the types related to ocean diagnostics
  !>
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!
!<Optimize:inUse>
  SUBROUTINE construct_oce_diagnostics( patch_3D, p_os, datestring )
    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    CHARACTER(LEN=32)                      :: datestring
    
    !local variable
    INTEGER :: i,ist
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_diagnostics:construct_oce_diagnostics')
    !-----------------------------------------------------------------------
    CHARACTER(LEN=linecharacters) :: headerline
    INTEGER  :: nblks_e,nblks_v,jb,jc,jk, region_index,start_index,end_index
    REAL(wp) :: surface_area, surface_height, prism_vol, prism_area, column_volume
    
    TYPE(t_patch), POINTER        :: p_patch
    TYPE(t_subset_range), POINTER :: owned_cells
    INTEGER, POINTER              :: regions(:,:)
    TYPE(t_ocean_regions)         :: ocean_regions
    !-----------------------------------------------------------------------
    p_patch => patch_3D%p_patch_2d(1)
    regions => patch_3D%regio_c
    !-----------------------------------------------------------------------
    owned_cells => p_patch%cells%owned
    !-----------------------------------------------------------------------
    
    CALL message (TRIM(routine), 'start')
    ALLOCATE(oce_ts)
    
    ALLOCATE(oce_ts%oce_diagnostics(0:nsteps))
    
    oce_ts%oce_diagnostics(0:nsteps)%volume                     = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%kin_energy                 = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%pot_energy                 = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%total_energy               = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%total_salt                 = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%vorticity                  = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%enstrophy                  = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%potential_enstrophy        = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%absolute_vertical_velocity = 0.0_wp

    oce_ts%oce_diagnostics(0:nsteps)%HeatFlux_ShortWave         = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%HeatFlux_LongWave          = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%HeatFlux_Sensible          = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%HeatFlux_Latent            = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%HeatFlux_Total             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_Precipitation     = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_SnowFall          = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_Evaporation       = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_Runoff            = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_TotalSalt         = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_TotalOcean        = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_TotalIce          = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_VolumeIce         = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_VolumeTotal       = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%HeatFlux_Relax             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%FrshFlux_Relax             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%TempFlux_Relax             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%SaltFlux_Relax             = 0.0_wp

    oce_ts%oce_diagnostics(0:nsteps)%ice_volume_nh              = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%ice_volume_sh              = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%ice_extent_nh              = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%ice_extent_sh              = 0.0_wp

    ! through flows
    oce_ts%oce_diagnostics(0:nsteps)%gibraltar                  = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%denmark_strait             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%drake_passage              = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%indonesian_throughflow     = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%scotland_iceland           = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%mozambique                 = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%framStrait                 = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%beringStrait               = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%barentsOpening             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%agulhas                    = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%agulhas_long               = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%agulhas_longer             = 0.0_wp

    oce_ts%oce_diagnostics(0:nsteps)%t_mean_na_200m             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%t_mean_na_800m             = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%ice_ocean_heat_budget      = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%ice_ocean_salinity_budget  = 0.0_wp
    oce_ts%oce_diagnostics(0:nsteps)%ice_ocean_volume_budget    = 0.0_wp
    
    DO i=0,nsteps
      ALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer))
      oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer) = 0.0_wp
    END DO
    
    ! open textfile for global timeseries
    diag_fname = 'oce_diagnostics-'//TRIM(datestring)//'.txt'
    diag_unit = find_next_free_unit(10,999)
    OPEN (UNIT=diag_unit,FILE=diag_fname,IOSTAT=ist,Recl=linecharacters)
    ! header of the text file
    headerline = ''
    ! * add timestep columns
    WRITE(headerline,'(a)') 'step date time'
    ! * add columne for each monitored variable
    DO i=1,SIZE(oce_ts%names)
      WRITE(headerline,'(a,a,a)')TRIM(headerline),' ',TRIM(oce_ts%names(i))
    END DO
    WRITE(diag_unit,'(a)')TRIM(headerline)

    ! CALL check_global_indexes(p_patch)
   
    ! open file for MOC - extraordinary at this time
    moc_fname='MOC.'//TRIM(datestring)
    moc_unit = find_next_free_unit(10,99)
    OPEN (moc_unit,FILE=moc_fname,FORM='unformatted')
    WRITE(message_text,'(2a)') ' MOC-file opened successfully, filename=',TRIM(moc_fname)
    CALL message (TRIM(routine), message_text)
    
    ! compute subsets for given sections path allong edges
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(1)%subset,      &
      & orientation = oce_sections(1)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array = gibraltar,            &
      & subset_name = 'gibraltar')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(2)%subset,      &
      & orientation = oce_sections(2)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =denmark_strait,        &
      & subset_name = 'denmark_strait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(3)%subset,      &
      & orientation = oce_sections(3)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =drake_passage,         &
      & subset_name = 'drake_passage')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(4)%subset,      &
      & orientation = oce_sections(4)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =indonesian_throughflow,&
      & subset_name = 'indonesian_throughflow')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(5)%subset,      &
      & orientation = oce_sections(5)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =scotland_iceland,      &
      & subset_name = 'scotland_iceland')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(6)%subset,      &
      & orientation = oce_sections(6)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =mozambique,      &
      & subset_name = 'mozambique')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(7)%subset,      &
      & orientation = oce_sections(7)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =framStrait,      &
      & subset_name = 'framStrait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(8)%subset,      &
      & orientation = oce_sections(8)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =beringStrait,      &
      & subset_name = 'beringStrait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(9)%subset,      &
      & orientation = oce_sections(9)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =barentsOpening,      &
      & subset_name = 'barentsOpening')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(10)%subset,      &
      & orientation = oce_sections(10)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas,      &
      & subset_name = 'agulhas')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(11)%subset,      &
      & orientation = oce_sections(11)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas_long,      &
      & subset_name = 'agulhas_long')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(12)%subset,      &
      & orientation = oce_sections(12)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas_longer,      &
      & subset_name = 'agulhas_longer')

!     CALL finish("","")
    
    surface_area   = 0.0_wp
    surface_height = 0.0_wp
    prism_vol      = 0.0_wp
    prism_area     = 0.0_wp
    ! compute regional ocean volumes
    DO jb = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
        
        ! area
        prism_area               = p_patch%cells%area(jc,jb)
        ocean_region_areas%total = ocean_region_areas%total + prism_area
        
        ! volume
        CALL compute_vertical_volume(jb,jc, &
          & prism_area, &
          & p_os%p_prog(nnew(1))%h(jc,jb), &
          & patch_3D%p_patch_1d(1)%prism_thick_c(jc,:,jb), &
          & patch_3D%p_patch_1d(1)%dolic_c(jc,jb), &
          & column_volume)
        ocean_region_volumes%total = ocean_region_volumes%total + column_volume
        
        region_index = regions(jc,jb)
        IF (ocean_regions%greenland_iceland_norwegian_sea == region_index) THEN
          ocean_region_volumes%greenland_iceland_norwegian_sea = &
            & ocean_region_volumes%greenland_iceland_norwegian_sea + column_volume
          ocean_region_areas%greenland_iceland_norwegian_sea   = &
            & ocean_region_areas%greenland_iceland_norwegian_sea + prism_area
          
        ELSEIF (ocean_regions%arctic_ocean == region_index) THEN
          ocean_region_volumes%arctic_ocean                    = ocean_region_volumes%arctic_ocean + column_volume
          ocean_region_areas%arctic_ocean                      = ocean_region_areas%arctic_ocean + prism_area
          
        ELSEIF (ocean_regions%labrador_sea == region_index) THEN
          ocean_region_volumes%labrador_sea                    = ocean_region_volumes%labrador_sea + column_volume
          ocean_region_areas%labrador_sea                      = ocean_region_areas%labrador_sea + prism_area
          
        ELSEIF (ocean_regions%north_atlantic == region_index) THEN
          ocean_region_volumes%north_atlantic                  = ocean_region_volumes%north_atlantic + column_volume
          ocean_region_areas%north_atlantic                    = ocean_region_areas%north_atlantic + prism_area
          
        ELSEIF (ocean_regions%tropical_atlantic == region_index) THEN
          ocean_region_volumes%tropical_atlantic               = ocean_region_volumes%tropical_atlantic + column_volume
          ocean_region_areas%tropical_atlantic                 = ocean_region_areas%tropical_atlantic + prism_area
          
        ELSEIF (ocean_regions%southern_ocean == region_index) THEN
          ocean_region_volumes%southern_ocean                  = ocean_region_volumes%southern_ocean + column_volume
          ocean_region_areas%southern_ocean                    = ocean_region_areas%southern_ocean + prism_area
          
        ELSEIF (ocean_regions%indian_ocean == region_index) THEN
          ocean_region_volumes%indian_ocean                    = ocean_region_volumes%indian_ocean + column_volume
          ocean_region_areas%indian_ocean                      = ocean_region_areas%indian_ocean + prism_area
          
        ELSEIF (ocean_regions%tropical_pacific == region_index) THEN
          ocean_region_volumes%tropical_pacific                = ocean_region_volumes%tropical_pacific + column_volume
          ocean_region_areas%tropical_pacific                  = ocean_region_areas%tropical_pacific + prism_area
          
        ELSEIF (ocean_regions%north_pacific == region_index) THEN
          ocean_region_volumes%north_pacific                   = ocean_region_volumes%north_pacific + column_volume
          ocean_region_areas%north_pacific                     = ocean_region_areas%north_pacific + prism_area
          
        ELSEIF (ocean_regions%caribbean == region_index) THEN
          ocean_region_volumes%caribbean                       = ocean_region_volumes%caribbean + column_volume
          ocean_region_areas%caribbean                         = ocean_region_areas%caribbean + prism_area
        END IF
        
      END DO
    END DO
    ! compute global values
    
    CALL disable_sync_checks()
    ocean_region_volumes%total = global_sum_array(ocean_region_volumes%total)
    ocean_region_areas%total   = global_sum_array(ocean_region_areas%total)
    CALL enable_sync_checks()
    
    CALL message (TRIM(routine), 'end')
  END SUBROUTINE construct_oce_diagnostics
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE compute_vertical_volume(jb,jc,prism_area,surface_height,thicknesses,max_vertical_level,volume)
    INTEGER,  INTENT(in)  :: jb,jc,max_vertical_level
    REAL(wp), INTENT(in)  :: prism_area, surface_height, thicknesses(:)
    REAL(wp), INTENT(inout) :: volume
    
    INTEGER :: jk
    REAL(wp) :: surface_height_,prism_vol_
    
    volume  = 0.0_wp
    DO jk = 1,max_vertical_level
      !local volume
      surface_height_ = MERGE(surface_height,0.0_wp, 1 == jk)
      prism_vol_      = prism_area * (thicknesses(jk) + surface_height_)
      !Fluid volume wrt lsm
      volume          = volume + prism_vol_
    END DO
  END SUBROUTINE compute_vertical_volume
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !  !The destructor of the types related to ocean diagnostics
  !>
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_oce_diagnostics()
    !
    !
    !local variables
    INTEGER :: i,iret
    CHARACTER(LEN=max_char_length)  :: linkname
    CHARACTER(LEN=max_char_length)  :: message_text
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_diagnostics:destruct_oce_diagnostics')
    !-----------------------------------------------------------------------
    DO i=0,nsteps
      DEALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content)
    END DO
    DEALLOCATE(oce_ts%oce_diagnostics)
    DEALLOCATE(oce_ts)
    ! close the global diagnostics text file and the SRV MOC file
    CLOSE(UNIT=diag_unit)
    CLOSE(UNIT=moc_unit)
    ! create a link to the last diagnostics file
    linkname = 'oce_diagnostics.txt'
    IF (util_islink(TRIM(linkname))) THEN
      iret = util_unlink(TRIM(linkname))
    ENDIF
    iret = util_symlink(TRIM(diag_fname),TRIM(linkname))
    WRITE(message_text,'(t1,a,t50,a)') TRIM(diag_fname), TRIM(linkname)
    CALL message('',message_text)
    
    CALL message (TRIM(routine), 'end')
  END SUBROUTINE destruct_oce_diagnostics
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! !  calculate_oce_diagnostics
  !
  ! @par Revision History
  ! Developed  by  Peter Korn, MPI-M (2010).
  !
  SUBROUTINE calc_slow_oce_diagnostics(patch_3D, p_os, p_sfc_flx, p_ice, &
    & timestep, datetime)
    TYPE(t_patch_3d ),TARGET, INTENT(in)    :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_sfc_flx),    INTENT(in)          :: p_sfc_flx
    TYPE (t_sea_ice),   INTENT(in)          :: p_ice
    INTEGER, INTENT(in) :: timestep
    TYPE(t_datetime), INTENT(in)            :: datetime
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index!,i_startblk_c, i_endblk_c,
    INTEGER :: jk,jc,jb!,je
    INTEGER :: i_no_t, i
    REAL(wp):: prism_vol, surface_height, prism_area, surface_area, z_w
    INTEGER :: reference_timestep
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp) :: sflux
    
    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_oce_monitor),  POINTER :: monitor
    CHARACTER(LEN=linecharacters) :: line, nvars
    CHARACTER(LEN=linecharacters) :: fmt_string, real_fmt
    CHARACTER(LEN=date_len)       :: datestring
    REAL(wp), PARAMETER :: equator = 0.00001_wp
    TYPE(t_ocean_regions)         :: ocean_regions
    
    !-----------------------------------------------------------------------
    p_patch        => patch_3D%p_patch_2d(1)
    owned_cells    => p_patch%cells%owned
    !-----------------------------------------------------------------------
    monitor        => oce_ts%oce_diagnostics(timestep)
    surface_area   = 0.0_wp
    surface_height = 0.0_wp
    prism_vol      = 0.0_wp
    prism_area     = 0.0_wp
    z_w            = 0.0_wp
    CALL datetime_to_string(datestring, datetime, plain=.TRUE.)
    
    !cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
    SELECT CASE (iswm_oce)
    CASE (1) ! shallow water mode
      !Potential energy in SW-casep_patch%patch_oce%del_zlev_m(1)
      DO jb = owned_cells%start_block,owned_cells%end_block
        CALL get_index_range(owned_cells, jb, start_cell_index, end_cell_index)
        
        DO jc =  start_cell_index, end_cell_index
          IF ( patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            prism_vol      = p_patch%cells%area(jc,jb)*patch_3D%p_patch_1d(1)%prism_thick_c(jc,1,jb)
            monitor%volume = monitor%volume + prism_vol
            
            
            monitor%pot_energy = monitor%pot_energy &
              & + 0.5_wp*grav* prism_vol*patch_3D%p_patch_1d(1)%prism_thick_c(jc,1,jb)
            
            monitor%kin_energy = monitor%kin_energy + p_os%p_diag%kin(jc,1,jb)*prism_vol
            
            monitor%total_energy=monitor%kin_energy+monitor%pot_energy
            DO i_no_t=1, no_tracer
              monitor%tracer_content(i_no_t) = monitor%tracer_content(i_no_t)&
                & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,1,jb,i_no_t)
            END DO
          ENDIF
        END DO
      END DO
      
    CASE default !3D model
      DO jb = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, jb, start_cell_index, end_cell_index)
        !We are dealing with the surface layer first
        DO jc =  start_cell_index, end_cell_index
          
          ! area
          prism_area   = p_patch%cells%area(jc,jb)
          surface_area = surface_area + prism_area
          ! sum of top layer vertical velocities abolsute values
          monitor%absolute_vertical_velocity = &
            & monitor%absolute_vertical_velocity + ABS(p_os%p_diag%w(jc,1,jb))*prism_area
          
          monitor%HeatFlux_ShortWave   = monitor%HeatFlux_ShortWave   + p_sfc_flx%HeatFlux_ShortWave(jc,jb)*prism_area
          monitor%HeatFlux_LongWave    = monitor%HeatFlux_LongWave    + p_sfc_flx%HeatFlux_LongWave (jc,jb)*prism_area
          monitor%HeatFlux_Sensible    = monitor%HeatFlux_Sensible    + p_sfc_flx%HeatFlux_Sensible (jc,jb)*prism_area
          monitor%HeatFlux_Latent      = monitor%HeatFlux_Latent      + p_sfc_flx%HeatFlux_Latent   (jc,jb)*prism_area
          monitor%HeatFlux_Total       = monitor%HeatFlux_Total       + p_sfc_flx%HeatFlux_Total    (jc,jb)*prism_area
          monitor%FrshFlux_Precipitation  = monitor%FrshFlux_Precipitation  + p_sfc_flx%FrshFlux_Precipitation(jc,jb)*prism_area
          monitor%FrshFlux_SnowFall    = monitor%FrshFlux_SnowFall    + p_sfc_flx%FrshFlux_SnowFall(jc,jb)*prism_area
          monitor%FrshFlux_Evaporation    = monitor%FrshFlux_Evaporation    + p_sfc_flx%FrshFlux_Evaporation(jc,jb)*prism_area
          monitor%FrshFlux_Runoff  = monitor%FrshFlux_Runoff  + p_sfc_flx%FrshFlux_Runoff(jc,jb)*prism_area
          monitor%FrshFlux_TotalSalt   = monitor%FrshFlux_TotalSalt   + p_sfc_flx%FrshFlux_TotalSalt(jc,jb)*prism_area
          monitor%FrshFlux_TotalOcean    = monitor%FrshFlux_TotalOcean    + p_sfc_flx%FrshFlux_TotalOcean(jc,jb)*prism_area
          monitor%FrshFlux_TotalIce    = monitor%FrshFlux_TotalIce    + p_sfc_flx%FrshFlux_TotalIce(jc,jb)*prism_area
          monitor%FrshFlux_VolumeIce   = monitor%FrshFlux_VolumeIce   + p_sfc_flx%FrshFlux_VolumeIce (jc,jb)*prism_area
          monitor%FrshFlux_VolumeTotal  = monitor%FrshFlux_VolumeTotal  + p_sfc_flx%FrshFlux_VolumeTotal(jc,jb)*prism_area
          monitor%HeatFlux_Relax = monitor%HeatFlux_Relax + p_sfc_flx%HeatFlux_Relax(jc,jb)*prism_area
          monitor%FrshFlux_Relax = monitor%FrshFlux_Relax + p_sfc_flx%FrshFlux_Relax(jc,jb)*prism_area
          monitor%TempFlux_Relax = monitor%TempFlux_Relax + p_sfc_flx%TempFlux_Relax(jc,jb)*prism_area
          monitor%SaltFlux_Relax = monitor%SaltFlux_Relax + p_sfc_flx%SaltFlux_Relax(jc,jb)*prism_area
          
          ! northern hemisphere
          IF (p_patch%cells%center(jc,jb)%lat > equator) THEN
            monitor%ice_volume_nh  = monitor%ice_volume_nh + prism_area*SUM(p_ice%vol(jc,:,jb)*p_ice%conc(jc,:,jb))
            monitor%ice_extent_nh  = monitor%ice_extent_nh + p_ice%concsum(jc,jb)*prism_area
          ELSE
            ! southern hemisphere
            monitor%ice_volume_sh  = monitor%ice_volume_sh + prism_area*SUM(p_ice%vol(jc,:,jb)*p_ice%conc(jc,:,jb))
            monitor%ice_extent_sh  = monitor%ice_extent_sh + p_ice%concsum(jc,jb)*prism_area
          END IF

          ! ice budgets
          ! heat
          ! 
          ! salinity
          ! volume
          
          DO jk = 1,patch_3D%p_patch_1d(1)%dolic_c(jc,jb)
            
            !local volume
            surface_height = MERGE(p_os%p_prog(nnew(1))%h(jc,jb),0.0_wp, 1 == jk)
            prism_vol      = prism_area * (patch_3D%p_patch_1d(1)%prism_thick_c(jc,jk,jb) + surface_height)
            
            !Fluid volume
            monitor%volume = monitor%volume + prism_vol
            
            !kinetic energy
            monitor%kin_energy = monitor%kin_energy + p_os%p_diag%kin(jc,jk,jb)*prism_vol
            
            !Potential energy
            IF(jk==1)THEN
              z_w = (p_os%p_diag%w(jc,jk,jb)*p_os%p_prog(nold(1))%h(jc,jb)&
                & +p_os%p_diag%w(jc,jk+1,jb)*0.5_wp*patch_3D%p_patch_1d(1)%del_zlev_i(jk))&
                & /(0.5_wp*patch_3D%p_patch_1d(1)%del_zlev_i(jk)+p_os%p_prog(nold(1))%h(jc,jb))
            ELSEIF(jk>1.AND.jk<n_zlev)THEN
              z_w = (p_os%p_diag%w(jc,jk,jb)*patch_3D%p_patch_1d(1)%del_zlev_i(jk)&
                & +p_os%p_diag%w(jc,jk+1,jb)*patch_3D%p_patch_1d(1)%del_zlev_i(jk+1))&
                & /(patch_3D%p_patch_1d(1)%del_zlev_i(jk)+patch_3D%p_patch_1d(1)%del_zlev_i(jk+1))
            ENDIF
            monitor%pot_energy = monitor%pot_energy + grav*z_w* p_os%p_diag%rho(jc,jk,jb)* prism_vol
            
            !Tracer content
            DO i_no_t=1, no_tracer
              monitor%tracer_content(i_no_t) = &
                & monitor%tracer_content(i_no_t) + prism_vol*p_os%p_prog(nold(1))%tracer(jc,jk,jb,i_no_t)
            END DO
          END DO
        END DO
      END DO
      
    END SELECT
    
    ! compute global sums {
    CALL disable_sync_checks()
    monitor%volume                     = global_sum_array(monitor%volume)
    surface_area                       = global_sum_array(surface_area)
    monitor%kin_energy                 = global_sum_array(monitor%kin_energy)/monitor%volume
    monitor%pot_energy                 = global_sum_array(monitor%pot_energy)/monitor%volume
    monitor%total_energy               = global_sum_array(monitor%total_energy)/monitor%volume
!   monitor%total_salt                 = calc_total_salt_content(p_patch, patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:),&
!       &                                                        p_ice, p_os,p_sfc_flx,p_ice%zUnderIce)
    monitor%vorticity                  = global_sum_array(monitor%vorticity)
    monitor%enstrophy                  = global_sum_array(monitor%enstrophy)
    monitor%potential_enstrophy        = global_sum_array(monitor%potential_enstrophy)
    monitor%absolute_vertical_velocity = global_sum_array(monitor%absolute_vertical_velocity)/surface_area
    monitor%HeatFlux_ShortWave         = global_sum_array(monitor%HeatFlux_ShortWave)/surface_area
    monitor%HeatFlux_LongWave          = global_sum_array(monitor%HeatFlux_LongWave )/surface_area
    monitor%HeatFlux_Sensible          = global_sum_array(monitor%HeatFlux_Sensible )/surface_area
    monitor%HeatFlux_Latent            = global_sum_array(monitor%HeatFlux_Latent   )/surface_area
    monitor%HeatFlux_Total             = global_sum_array(monitor%HeatFlux_Total    )/surface_area
    monitor%FrshFlux_Precipitation     = global_sum_array(monitor%FrshFlux_Precipitation)/surface_area
    monitor%FrshFlux_SnowFall          = global_sum_array(monitor%FrshFlux_SnowFall)/surface_area
    monitor%FrshFlux_Evaporation       = global_sum_array(monitor%FrshFlux_Evaporation)/surface_area
    monitor%FrshFlux_Runoff            = global_sum_array(monitor%FrshFlux_Runoff)/surface_area
    monitor%FrshFlux_TotalSalt         = global_sum_array(monitor%FrshFlux_TotalSalt)/surface_area
    monitor%FrshFlux_TotalOcean        = global_sum_array(monitor%FrshFlux_TotalOcean)/surface_area
    monitor%FrshFlux_TotalIce          = global_sum_array(monitor%FrshFlux_TotalIce)/surface_area
    monitor%FrshFlux_VolumeIce         = global_sum_array(monitor%FrshFlux_VolumeIce)/surface_area
    monitor%FrshFlux_VolumeTotal       = global_sum_array(monitor%FrshFlux_VolumeTotal)/surface_area
    monitor%HeatFlux_Relax             = global_sum_array(monitor%HeatFlux_Relax)/surface_area
    monitor%FrshFlux_Relax             = global_sum_array(monitor%FrshFlux_Relax)/surface_area
    monitor%SaltFlux_Relax             = global_sum_array(monitor%SaltFlux_Relax)/surface_area
    monitor%TempFlux_Relax             = global_sum_array(monitor%TempFlux_Relax)/surface_area
    monitor%ice_volume_nh              = global_sum_array(monitor%ice_volume_nh)/1.0e9_wp
    monitor%ice_volume_sh              = global_sum_array(monitor%ice_volume_sh)/1.0e9_wp
    monitor%ice_extent_nh              = global_sum_array(monitor%ice_extent_nh)/1.0e6_wp
    monitor%ice_extent_sh              = global_sum_array(monitor%ice_extent_sh)/1.0e6_wp
    DO i_no_t=1,no_tracer
      monitor%tracer_content(i_no_t) = global_sum_array(monitor%tracer_content(i_no_t))
    END DO
    CALL enable_sync_checks()
    ! fluxes through given paths
    IF (my_process_is_stdio()) &
      & WRITE(0,*) "---------------  fluxes --------------------------------"
    DO i=1,oce_section_count
      sflux = section_flux(oce_sections(i), p_os%p_prog(nnew(1))%vn)
      !
      ! #slo# disabled since subset%block is not allocated (#3759, HPC_sun_debug)
      ! #ifdef NOMPI
      !     IF (my_process_is_stdio()) &
      !       & write(0,*) oce_sections(i)%subset%name, ":", sflux, 'at edges:',oce_sections(i)%subset%block
      ! #else
      IF (my_process_is_stdio()) &
        & WRITE(0,*) oce_sections(i)%subset%name, ":", sflux
      ! #endif
      
      SELECT CASE (i)
      CASE (1)
        monitor%gibraltar              = sflux*rho_ref
      CASE (2)
        monitor%denmark_strait         = sflux*rho_ref
      CASE (3)
        monitor%drake_passage          = sflux*rho_ref
      CASE (4)
        monitor%indonesian_throughflow = sflux*rho_ref
      CASE (5)
        monitor%scotland_iceland       = sflux*rho_ref
      CASE (6)
        monitor%mozambique             = sflux*rho_ref
      CASE (7)
        monitor%framStrait             = sflux*rho_ref
      CASE (8)
        monitor%beringStrait           = sflux*rho_ref
      CASE (9)
        monitor%barentsOpening         = sflux*rho_ref
      CASE (10)
        monitor%agulhas                = sflux*rho_ref
      CASE (11)
        monitor%agulhas_long           = sflux*rho_ref
      CASE (12)
        monitor%agulhas_longer         = sflux*rho_ref
      END SELECT
    ENDDO
    IF (my_process_is_stdio()) &
      & WRITE(0,*) "---------------  end fluxes ----------------------------"
    
    
    IF (my_process_is_stdio()) THEN
      ! write things to diagnostics output file
      real_fmt   = 'es26.18'
      ! * number of non-tracer diag. variables
      WRITE(nvars,'(i3)') SIZE(oce_ts%names)-no_tracer
      WRITE(fmt_string,'(a)') '(i15.15,1x,a,1x,'//TRIM(ADJUSTL(nvars))//TRIM(real_fmt)//')'
      ! create date and time string
      ! * non-tracer diags
      WRITE(line,fmt_string) &
        & timestep, &
        & TRIM(datestring), &
        & monitor%volume, &
        & monitor%kin_energy, &
        & monitor%pot_energy, &
        & monitor%total_energy, &
        & monitor%total_salt, &
        & monitor%vorticity, &
        & monitor%enstrophy, &
        & monitor%potential_enstrophy, &
        & monitor%absolute_vertical_velocity, &
        & monitor%HeatFlux_ShortWave, &
        & monitor%HeatFlux_LongWave , &
        & monitor%HeatFlux_Sensible , &
        & monitor%HeatFlux_Latent   , &
        & monitor%HeatFlux_Total,     &
        & monitor%FrshFlux_Precipitation, &
        & monitor%FrshFlux_SnowFall, &
        & monitor%FrshFlux_Evaporation, &
        & monitor%FrshFlux_Runoff, &
        & monitor%FrshFlux_TotalSalt, &
        & monitor%FrshFlux_TotalOcean, &
        & monitor%FrshFlux_TotalIce, &
        & monitor%FrshFlux_VolumeIce, &
        & monitor%FrshFlux_VolumeTotal, &
        & monitor%HeatFlux_Relax, &
        & monitor%FrshFlux_Relax, &
        & monitor%TempFlux_Relax, &
        & monitor%SaltFlux_Relax, &
        & monitor%ice_volume_nh, &
        & monitor%ice_volume_sh, &
        & monitor%ice_extent_nh, &
        & monitor%ice_extent_sh, &
        & monitor%gibraltar, &
        & monitor%denmark_strait, &
        & monitor%drake_passage,  &
        & monitor%indonesian_throughflow, &
        & monitor%scotland_iceland, &
        & monitor%mozambique    , &
        & monitor%framStrait    , &
        & monitor%beringStrait  , &
        & monitor%barentsOpening, &
        & monitor%agulhas       , &
        & monitor%agulhas_long  , &
        & monitor%agulhas_longer, &
        & monitor%t_mean_na_200m, &
        & monitor%t_mean_na_800m, &
        & monitor%ice_ocean_heat_budget, &
        & monitor%ice_ocean_salinity_budget, &
        & monitor%ice_ocean_volume_budget

      ! * tracers
      DO i_no_t=1,no_tracer
        WRITE(line,'(a,'//TRIM(real_fmt)//')') TRIM(line),monitor%tracer_content(i_no_t)
      END DO

      WRITE(diag_unit,'(a)') TRIM(line)
    END IF
    
  END SUBROUTINE calc_slow_oce_diagnostics
  
!<Optimize:inUse>
  SUBROUTINE calc_fast_oce_diagnostics(p_patch, dolic, prism_thickness,depths, p_diag)
    TYPE(t_patch ),TARGET :: p_patch
    INTEGER,  POINTER :: dolic(:,:)
    REAL(wp), POINTER :: prism_thickness(:,:,:)
    REAL(wp), INTENT(in)              :: depths(:)
    TYPE(t_hydro_ocean_diag), TARGET :: p_diag
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index!,i_startblk_c, i_endblk_c,
    INTEGER :: jk,jc,jb!,je
    
    TYPE(t_subset_range), POINTER :: owned_cells
    !-----------------------------------------------------------------------
    owned_cells    => p_patch%cells%owned
    
    !cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
    SELECT CASE (iswm_oce)
    CASE (1) ! shallow water mode
      
    CASE default !3D model
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index) SCHEDULE(dynamic)
      DO jb = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, jb, start_cell_index, end_cell_index)
        !We are dealing with the surface layer first
        DO jc =  start_cell_index, end_cell_index
          
         !p_diag%mld(jc,jb) = calc_mixed_layer_depth(p_diag%zgrad_rho(jc,:,jb),&
         !  & 0.125_wp, &
         !  & dolic(jc,jb), &
         !  & prism_thickness(jc,:,jb), &
         !  & depths(1))
          ! save the maximal mixed layer depth  - RESET TO 0 ECHT OUTPUT TIMESTEP
          p_diag%mld(jc,jb) = MAX(calc_mixed_layer_depth(p_diag%zgrad_rho(jc,:,jb),&
            & 0.125_wp, &
            & dolic(jc,jb), &
            & prism_thickness(jc,:,jb), &
            & depths(1)),p_diag%mld(jc,jb))
          p_diag%condep(jc,jb) = calc_condep(p_diag%zgrad_rho(jc,:,jb), dolic(jc,jb))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO 
      CALL dbg_print('Diag: mld',p_diag%mld,str_module,4,in_subset=owned_cells)
    END SELECT
  END SUBROUTINE calc_fast_oce_diagnostics
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION section_flux(in_oce_section, velocity_values)
    TYPE(t_oce_section) :: in_oce_section
    REAL(wp), POINTER :: velocity_values(:,:,:)
    
    INTEGER :: i, k, edge_idx, edge_block
    REAL(wp) :: oriented_length
    REAL(wp), ALLOCATABLE :: flux_weights(:,:)
    TYPE(t_grid_edges), POINTER ::  edges
    TYPE(t_patch_vert),POINTER :: patch_vertical
    
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_ocean_diagnostics:section_flux'
    
    edges          => in_oce_section%subset%patch%edges
    patch_vertical => in_oce_section%subset%patch_3d%p_patch_1d(1)
    
    ! calculate weights
    ! flux_weights can also be preallocated
    ALLOCATE(flux_weights(n_zlev, MAX(in_oce_section%subset%SIZE, 1)))
    flux_weights(:,:) = 0.0_wp
    DO i=1, in_oce_section%subset%SIZE
      
      edge_idx   = in_oce_section%subset%idx(i)
      edge_block = in_oce_section%subset%BLOCK(i)
      oriented_length = edges%primal_edge_length(edge_idx, edge_block) * &
        & in_oce_section%orientation(i) ! this can also be pre-calculated and stored in in_oce_section%orientation
      
      !write(0,*) "oriented_length:",  oriented_length
      
      DO k=1, n_zlev
        flux_weights(k, i) = patch_vertical%prism_thick_e(edge_idx, k, edge_block) * oriented_length ! maybe also use slm
        !write(0,*) i, k, in_oce_section%subset%name, " flux_weights:",  flux_weights(k, i), &
        !  & patch_vertical%prism_thick_e(edge_idx, k, edge_block)
        !write(0,*) i, k, in_oce_section%subset%name, " velocity_value:", velocity_values(edge_idx, k, edge_block)
      ENDDO
      
    ENDDO
    
    
    section_flux = subset_sum(                           &
      & values                 = velocity_values,        &
      & indexed_subset         = in_oce_section%subset,  &
      & subset_indexed_weights = flux_weights)
    
    DEALLOCATE(flux_weights)
    
    !write(0,*) get_my_mpi_work_id(), ": section_flux on subset ", in_oce_section%subset%name, ":", &
    !  & section_flux, in_oce_section%subset%size
    
  END FUNCTION section_flux
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of meridional overturning circulation (MOC)
  !
  !   Calculation of meridional overturning circulation for different basins
  !   (Atlantic, Pacific, Indian, global)
  !>
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2012).
  !!  based on code from MPIOM
  !
  ! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
  ! TODO: calculate the 1 deg resolution meridional distance
  !!
  SUBROUTINE calc_moc (p_patch, patch_3D, w, datetime)
    
    TYPE(t_patch), TARGET, INTENT(in)  :: p_patch
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: w(:,:,:)   ! vertical velocity at cell centers
    ! dims: (nproma,nlev+1,alloc_cell_blocks)
    TYPE(t_datetime), INTENT(in)       :: datetime
    !
    ! local variables
    ! INTEGER :: i
    INTEGER, PARAMETER ::  jbrei=3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
    INTEGER :: jb, jc, jk, start_index, end_index !, il_e, ib_e
    INTEGER :: lbrei, lbr, idate, itime
    INTEGER :: mpi_comm
    INTEGER(i8) :: i1,i2,i3,i4
    
    REAL(wp) :: z_lat, z_lat_deg, z_lat_dim
    REAL(wp) :: global_moc(180,n_zlev), atlant_moc(180,n_zlev), pacind_moc(180,n_zlev)
    REAL(dp) :: local_moc(180), res_moc(180)
    
    TYPE(t_subset_range), POINTER :: dom_cells
    
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_moc')
    
    !-----------------------------------------------------------------------
    
    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF
    
    global_moc(:,:) = 0.0_wp
    pacind_moc(:,:) = 0.0_wp
    atlant_moc(:,:) = 0.0_wp
    
    ! set barrier:
    ! CALL MPI_BARRIER(0)
    
    ! with all cells no sync is necessary
    !owned_cells => p_patch%cells%owned
    dom_cells   => p_patch%cells%in_domain
    
    !write(81,*) 'MOC: datetime:',datetime
    
    DO jk = 1, n_zlev   !  not yet on intermediate levels
      DO jb = dom_cells%start_block, dom_cells%end_block
        CALL get_index_range(dom_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          
          !  could be replaced by vertical loop to bottom
          IF ( patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            
            ! lbrei: corresponding latitude row of 1 deg extension
            !       1 south pole
            !     180 north pole
            z_lat = p_patch%cells%center(jc,jb)%lat
            z_lat_deg = z_lat*rad2deg
            lbrei = NINT(90.0_wp + z_lat_deg)
            lbrei = MAX(lbrei,1)
            lbrei = MIN(lbrei,180)
            
            ! get neighbor edge for scaling
            !   il_e = p_patch%cells%edge_idx(jc,jb,1)
            !   ib_e = p_patch%cells%edge_blk(jc,jb,1)
            
            ! z_lat_dim: scale to 1 deg resolution
            ! z_lat_dim: latitudinal extent of triangle divided by latitudinal smoothing extent
            !   z_lat_dim = p_patch%edges%primal_edge_length(il_e,ib_e) / &
            !     & (REAL(2*jbrei, wp) * 111111._wp*1.3_wp)
            z_lat_dim = 1.0_wp
            
            ! distribute MOC over (2*jbrei)+1 latitude rows
            !  - no weighting with latitudes done
            !  - lbrei: index of 180 X 1 deg meridional resolution
            DO lbr = -jbrei, jbrei
              lbrei = NINT(90.0_wp + z_lat_deg + REAL(lbr, wp) * z_lat_dim)
              lbrei = MAX(lbrei,1)
              lbrei = MIN(lbrei,180)
              
              global_moc(lbrei,jk) = global_moc(lbrei,jk) - &
              !  multiply with wet (or loop to bottom)
                & p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) * &
                & patch_3D%wet_c(jc,jk,jb) / &
                & REAL(2*jbrei + 1, wp)
              
              IF (patch_3D%basin_c(jc,jb) == 1) THEN         !  1: Atlantic; 0: Land
                
                atlant_moc(lbrei,jk) = atlant_moc(lbrei,jk) - &
                  & p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) * &
                  & patch_3D%wet_c(jc,jk,jb) / &
                  & REAL(2*jbrei + 1, wp)
              ELSE IF (patch_3D%basin_c(jc,jb) >= 2) THEN   !  2: Indian; 4: Pacific
                pacind_moc(lbrei,jk) = pacind_moc(lbrei,jk) - &
                  & p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) * &
                  & patch_3D%wet_c(jc,jk,jb) / &
                  & REAL(2*jbrei + 1, wp)
              END IF
              
            END DO
            
          END IF
        END DO
      END DO
      
      ! test parallelization:
      ! function field_sum_all using mpi_allreduce and working precisions wp does not exist
      ! res_moc(:) = p_field_sum_all_wp(global_moc(:,jk))
      ! res_moc(:) = p_field_sum_all_wp(atlant_moc(:,jk))
      ! res_moc(:) = p_field_sum_all_wp(pacind_moc(:,jk))
      
      ! function field_sum using mpi_reduce, then broadcast
      local_moc(:)     = REAL(global_moc(:,jk),dp)
      res_moc(:)       = p_field_sum(local_moc, mpi_comm)
      CALL p_bcast(res_moc(:), p_io, mpi_comm)
      global_moc(:,jk) = REAL(res_moc(:),wp)
      
      local_moc(:)     = REAL(atlant_moc(:,jk),dp)
      res_moc(:)       = p_field_sum(local_moc, mpi_comm)
      CALL p_bcast(res_moc(:), p_io, mpi_comm)
      atlant_moc(:,jk) = REAL(res_moc(:),wp)
      
      local_moc(:)     = REAL(pacind_moc(:,jk),dp)
      res_moc(:)       = p_field_sum(local_moc, mpi_comm)
      CALL p_bcast(res_moc(:), p_io, mpi_comm)
      pacind_moc(:,jk) = REAL(res_moc(:),wp)
      
    END DO  ! n_zlev-loop
    
    IF (my_process_is_stdio()) THEN
      DO lbr=179,1,-1   ! fixed to 1 deg meridional resolution
        
        global_moc(lbr,:)=global_moc(lbr+1,:)+global_moc(lbr,:)
        atlant_moc(lbr,:)=atlant_moc(lbr+1,:)+atlant_moc(lbr,:)
        pacind_moc(lbr,:)=pacind_moc(lbr+1,:)+pacind_moc(lbr,:)
        
      END DO
      
      ! write out MOC in extra format, file opened in mo_hydro_ocean_run  - integer*8
      !  - correct date in extra format - i.e YYYYMMDD - no time info
      !idate=datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute
      idate = datetime%year*10000+datetime%month*100+datetime%day
      itime = datetime%hour*100+datetime%minute
      WRITE(message_text,*) 'Write MOC at year =',datetime%year,', date =',idate,' time =', itime
      CALL message (TRIM(routine), message_text)
      
      DO jk = 1,n_zlev
        i1=INT(idate,i8)
        i2 = INT(777,i8)
        i3 = INT(patch_3D%p_patch_1d(1)%zlev_i(jk),i8)
        i4 = INT(180,i8)
        WRITE(moc_unit) i1,i2,i3,i4
        WRITE(moc_unit) (global_moc(lbr,jk),lbr=1,180)
        i2 = INT(778,i8)
        WRITE(moc_unit) i1,i2,i3,i4
        WRITE(moc_unit) (atlant_moc(lbr,jk),lbr=1,180)
        i2 = INT(779,i8)
        WRITE(moc_unit) i1,i2,i3,i4
        WRITE(moc_unit) (pacind_moc(lbr,jk),lbr=1,180)
        
      END DO
    END IF
    
  END SUBROUTINE calc_moc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of horizontal stream function
  !
  !>
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2012).
  !!  based on code from MPIOM
  !
  ! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
  !!
!<Optimize:inUse>
  SUBROUTINE calc_psi (patch_3D, u, prism_thickness, u_vint, datetime)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: u(:,:,:)     ! zonal velocity at cell centers
    REAL(wp), INTENT(in)               :: prism_thickness(:,:,:)       ! elevation on cell centers
    ! dims: (nproma,nlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)            :: u_vint(:,:)  ! barotropic zonal velocity on icon grid
    TYPE(t_datetime), INTENT(in)       :: datetime
    !
    ! local variables
    ! INTEGER :: i
    
    ! switch for writing stream function (not yet in namelist); 1: icon-grid; 2: regular grid output
    INTEGER, PARAMETER ::  idiag_psi = 1
    
    INTEGER, PARAMETER ::  nlat = 180                    ! meridional dimension of regular grid
    INTEGER, PARAMETER ::  nlon = 360                    ! zonal dimension of regular grid
    
    ! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
    INTEGER, PARAMETER ::  jsmth = 3
    INTEGER :: jb, jc, jk, start_index, end_index
    INTEGER :: jlat, jlon, jlt, jln, jltx, jlnx, jsmth2
    INTEGER(i8)        :: idate, iextra(4)
    
    
    REAL(wp) :: z_lat_deg, z_lon_deg, z_lat_dist, delta_z, rsmth
    REAL(wp) :: z_uint_reg(nlon,nlat)                     ! vertical integral on regular grid
    REAL(wp) :: psi_reg(nlon,nlat)                        ! horizontal stream function
    
    TYPE(t_patch), POINTER  :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, dom_cells
    
    !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_psi')
    
    !-----------------------------------------------------------------------
    
    psi_reg(:,:)    = 0.0_wp
    z_uint_reg(:,:) = 0.0_wp
    
    jsmth2          = 2*jsmth + 1
    rsmth           = REAL(jsmth2*jsmth2, wp)
    
    ! with all cells no sync is necessary
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    dom_cells => patch_2d%cells%in_domain
    
    ! (1) barotropic system:
    !     vertical integration of zonal velocity times vertical layer thickness [m/s*m]
!ICON_OMP_PARALLEL_DO PRIVATE(jc, jk, start_index, end_index) SCHEDULE(dynamic)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      u_vint(:,jb)     = 0.0_wp
      DO jc = start_index, end_index
      
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        
          u_vint(jc,jb) = u_vint(jc,jb) - u(jc,jk,jb) * prism_thickness(jc,jk,jb)
          
        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
                
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi:    u_vint    ', u_vint        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------
    
    IF (idiag_psi == 1) RETURN
    
    ! (2) distribute integrated zonal velocity (u*dz) on 1x1 deg grid
    !     this code is not mature yet
    
    ! in domain: count all cells only once
    DO jb = dom_cells%start_block, dom_cells%end_block
      CALL get_index_range(dom_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
        z_lat_deg = patch_2d%cells%center(jc,jb)%lat * rad2deg
        z_lon_deg = patch_2d%cells%center(jc,jb)%lon * rad2deg
        
        !  ! 0 <= lon <= 360 deg
        !  z_lon_deg = z_lon_deg + 180.0_wp
        
        ! jlat/jlon: corresponding latitude/longitude coordinates of 1 deg extension
        ! jlat: 1 = south of 89.0S; 89 = 1S-Eq.; 90 = Eq-1N;  180 = north of 89N
        ! jlon: 1 = 180W-179W; 180 = 1-0 deg west; 360 = 179E-180E
        
        jlat = NINT(91.0_wp + z_lat_deg)
        jlon = NINT(z_lon_deg + 180.5_wp)
        
        ! distribute stream function over rsmth=(2*jsmth+1)**2 lat/lon regular grid points
        !  - no weighting with latitudes done
        !  - no correction with regular lsm done
        DO jltx = jlat-jsmth, jlat+jsmth
          
          jlt = jltx
          IF (jlt <    1) jlt =      1-jlt  ! apply equatorwards
          IF (jlt > nlat) jlt = 2*nlat-jlt  ! apply equatorwards
          DO jlnx = jlon-jsmth, jlon+jsmth
            
            jln = jlnx
            IF (jln <    1) jln = jln+nlon  ! circular boundary
            IF (jln > nlon) jln = jln-nlon  ! circular boundary
            
            z_uint_reg(jln,jlt) = z_uint_reg(jln,jlt) + u_vint(jc,jb) / rsmth
            
            ! 99 format('J lat=',f8.2,' lon=',f8.2,' jlat=',i4,' jlon=',i4,' lsm=',i3, &
            !      &    ' jlt=',i4,  ' jln=',i4,' uint=',1p10e12.3)
            ! 98 format(' lat=',f8.2,' lon=',f8.2,' jlat=',i4,' jlon=',i4,' lsm=',i3, &
            !      &    ' uint=',1p10e12.3)
            !    if ((jlat==101 .and. jlon==270) &
            !      & write(82,99) z_lat_deg,z_lon_deg,jlat,jlon,v_base%lsm_c(jc,1,jb), &
            !      &              jlt,jln,z_uint_reg(jln,jlt)
            
          END DO
        END DO
        !    write(82,98) z_lat_deg,z_lon_deg,jlat,jlon,v_base%lsm_c(jc,1,jb),z_uint_reg(jlon,jlat)
        
      END DO
    END DO
    
    ! (3) calculate meridional integral on regular grid starting from south pole:
    
    DO jlt = nlat-1, 1, -1
      z_uint_reg(:,jlt) = z_uint_reg(:,jlt) + z_uint_reg(:,jlt+1)
    END DO
    
    ! (4) calculate stream function: scale with length of 1 deg*rho [m/s*m*m*kg/m3=kg/s]
    
    ! meridional distance of 1 deg
    ! ATTENTION - fixed 1 deg resolution should be related to icon-resolution
    z_lat_dist = 111111.0_wp  ! * 1.3_wp ??
    
    psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * rho_ref
    
    ! stream function on icon grid without calculation of meridional integral
    !  - tbd after interpolation to regular grid externally
    !  psi    (:,:) = u_vint    (:,:)              * rho_ref
    
    
    ! write out in extra format - integer*8
    idate = INT(datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute,i8)
    WRITE(0,*) 'write global PSI at iyear, idate:',datetime%year, idate
    
    iextra(1) = INT(idate,i8)
    iextra(2) = INT(780,i8)
    iextra(3) = INT(0,i8)
    iextra(4) = INT(nlon*nlat,i8)
    
    WRITE(80) (iextra(jb),jb=1,4)
    WRITE(80) ((psi_reg(jln,jlt),jln=1,nlon),jlt=1,nlat)
    
    DO jlat=1,nlat
      WRITE(82,*) 'jlat=',jlat
      WRITE(82,'(1p10e12.3)') (psi_reg(jlon,jlat),jlon=1,nlon)
    ENDDO
                
    !---------DEBUG DIAGNOSTICS-------------------------------------------
 !  idt_src=3  ! output print level (1-5, fix)
 !  CALL dbg_print('calc_psi:    psi_reg   ', psi_reg       , str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE calc_psi
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of horizontal stream function using normal velocity on edges
  !
  !>
  !>  Calculation of horizontal stream function using normal velocity on edges
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_psi_vn (patch_3D, vn, prism_thickness_e, op_coeff, u_vint, v_vint, datetime)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: vn(:,:,:)                 ! normal velocity at cell edges
    REAL(wp), INTENT(in)               :: prism_thickness_e(:,:,:)  ! elevation on cell edges
    ! dims: (nproma,nlev,alloc_edge_blocks)
    REAL(wp), INTENT(inout)            :: u_vint(:,:)               ! barotropic zonal velocity on cell centers
    REAL(wp), INTENT(inout)            :: v_vint(:,:)               ! barotropic meridional velocity on cell centers
    TYPE(t_operator_coeff),INTENT(in)  :: op_coeff
    TYPE(t_datetime), INTENT(in)       :: datetime
    !
    INTEGER  :: jb, jc, je, jk, start_index, end_index
    ! vertical integral vn on edges and in cartesian coordinates - no 2-dim mapping is available
    REAL(wp) :: vn_vint(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: u_2d(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !  scratch arrays for test
    REAL(wp) :: v_2d(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !  scratch arrays for test
    TYPE(t_cartesian_coordinates) :: vint_cc(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    
    TYPE(t_patch), POINTER  :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges, all_cells
    
    !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_psi_vn')
    
    !-----------------------------------------------------------------------
    
    patch_2d  => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL
    all_cells => patch_2d%cells%ALL
    
    vn_vint  (:,:,:)    = 0.0_wp
    vint_cc(:,:,:)%x(1) = 0.0_wp
    vint_cc(:,:,:)%x(2) = 0.0_wp
    vint_cc(:,:,:)%x(3) = 0.0_wp
    u_2d       (:,:)    = 0.0_wp
    v_2d       (:,:)    = 0.0_wp
    
    ! (1) barotropic system:
    !     vertical integration of normal velocity times vertical layer thickness [m/s*m]
!ICON_OMP_PARALLEL_DO PRIVATE(je, jk, start_index, end_index) SCHEDULE(dynamic)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_index, end_index)
      vn_vint(:,1,jb) = 0.0_wp
      DO je = start_index, end_index
      
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
        
          vn_vint(je,1,jb) = vn_vint(je,1,jb) - vn(je,jk,jb) * prism_thickness_e(je,jk,jb)
          
        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
    
    ! (2) remapping normal velocity to zonal and meridional velocity at cell centers
    CALL sync_patch_array(sync_e, patch_2D, vn_vint)       

    CALL map_edges2cell_3d(patch_3D, vn_vint, op_coeff, vint_cc)

    CALL sync_patch_array(sync_c, patch_2D, vint_cc(:,:,:)%x(1))
    CALL sync_patch_array(sync_c, patch_2D, vint_cc(:,:,:)%x(2))
    CALL sync_patch_array(sync_c, patch_2D, vint_cc(:,:,:)%x(3))

    ! calculate zonal and meridional velocity:
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
        CALL cvec2gvec(vint_cc(jc,1,jb)%x(1), vint_cc(jc,1,jb)%x(2), vint_cc(jc,1,jb)%x(3), &
          &            patch_2d%cells%center(jc,jb)%lon, patch_2d%cells%center(jc,jb)%lat,  &
          &            u_2d(jc,jb), v_2d(jc,jb))
!         &            u_vint(jc,jb), v_vint(jc,jb))
      END DO
    END DO
    !CALL sync_patch_array(sync_c, patch_2D, u_vint)       
    !CALL sync_patch_array(sync_c, patch_2D, v_vint)       
    CALL sync_patch_array(sync_c, patch_2D, u_2d)       
    CALL sync_patch_array(sync_c, patch_2D, v_2d)       

    ! hack for test: calc_psy for u_vint, calc_psi_vn for v_vint - accumulated and written out
    v_vint(:,:) = u_2d(:,:)
                
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi_vn: u_2d        ', u_2d        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi_vn: v_2d        ', v_2d        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    CALL dbg_print('calc_psi_vn: vint_cc%x(1)' ,vint_cc%x(1), str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------

    
  END SUBROUTINE calc_psi_vn
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!taken from MPIOM
!<Optimize:inUse>
  FUNCTION calc_condep(vertical_density_gradient,max_lev) result(condep)
    REAL(wp),INTENT(in)  :: vertical_density_gradient(n_zlev)
    INTEGER, INTENT(in)  :: max_lev
    INTEGER :: condep
    
    INTEGER :: jk
    INTEGER :: maxcondep     !< maximum convective penetration level
    REAL(wp) :: masked_vertical_density_gradient(n_zlev)
    
    condep = 1

    ! remove dbl_eps, which  is added in the vertical gradient computation
    masked_vertical_density_gradient = MAX(vertical_density_gradient - dbl_eps,0.0_wp)
    
    !! diagnose maximum convection level
    !! condep = maximum model level penetrated by vertically continous
    !! convection from the surface downward
    !! calculated over integration period ; it should be written out
    !! as snapshot at the end of the run
   maxcondep=1
    DO jk=2,max_lev
      IF (masked_vertical_density_gradient(jk) .ne. 0.0_wp) THEN
        maxcondep = jk
        EXIT
      ENDIF
    ENDDO
    
    condep = max0(maxcondep,condep)
  END FUNCTION calc_condep
!<Optimize:inUse>
  FUNCTION calc_mixed_layer_depth(vertical_density_gradient,critical_value,max_lev,thickness, depth_of_first_layer) &
    & result(mixed_layer_depth)
    REAL(wp), TARGET :: vertical_density_gradient(n_zlev)
    REAL(wp), INTENT(in)  :: critical_value
    INTEGER,  INTENT(in)  :: max_lev
    REAL(wp), INTENT(in)  :: thickness(n_zlev)
    REAL(wp), INTENT(in)  :: depth_of_first_layer
    
    REAL(wp) :: sigh        ,zzz
    REAL(wp) :: mixed_layer_depth
    REAL(wp) :: masked_vertical_density_gradient(n_zlev)
    INTEGER :: jk
    
    sigh              = critical_value
    mixed_layer_depth = depth_of_first_layer
    masked_vertical_density_gradient = MAX(vertical_density_gradient,0.0_wp)
    
    ! This diagnostic calculates the mixed layer depth.
    ! It uses the incremental density increase between two
    ! levels and substracts it from the initial density criterion (sigcrit)
    ! and adds the level thickness (zzz) to zmld. This is done till
    ! the accumulated density increase between the surface and
    ! layer k is sigcrit or sigh = O, respectively.
    
    ! stabio(k) = insitu density gradient
    ! sigh = remaining density difference
    
    DO jk = 2, max_lev
      IF (sigh .GT. 1.e-6_wp) THEN
        zzz               = MIN(sigh/(ABS(masked_vertical_density_gradient(jk))+1.0E-19_wp),thickness(jk))
        sigh              = MAX(0._wp, sigh-zzz*masked_vertical_density_gradient(jk))
        mixed_layer_depth = mixed_layer_depth + zzz
      ELSE
        sigh = 0._wp
      ENDIF
    ENDDO
    
  END FUNCTION calc_mixed_layer_depth

  FUNCTION calc_total_salt_content(p_patch, thickness, p_ice, p_os, surface_fluxes, zUnderIce, &
      & computation_type) RESULT(total_salt_content)
    TYPE(t_patch),POINTER                                 :: p_patch
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks),INTENT(IN) :: thickness,zUnderIce
    TYPE (t_sea_ice),       INTENT(IN)                    :: p_ice
    TYPE(t_hydro_ocean_state)                             :: p_os
    TYPE(t_sfc_flx)                                       :: surface_fluxes
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type
    REAL(wp)                                              :: total_salt_content

    REAL(wp), DIMENSION(nproma,n_zlev,p_patch%alloc_cell_blocks) :: salt

    salt               = calc_salt_content(p_patch, thickness, p_ice, p_os, surface_fluxes, zUnderIce, computation_type)
    total_salt_content = global_sum_array(salt)
  END FUNCTION calc_total_salt_content

  FUNCTION calc_salt_content(p_patch, thickness, p_ice, p_os, surface_fluxes, zUnderIce, &
      & computation_type) &
      & RESULT(salt)
    TYPE(t_patch),POINTER                                 :: p_patch
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks),INTENT(IN) :: zUnderIce
    REAL(wp),DIMENSION(nproma,n_zlev,p_patch%alloc_cell_blocks),INTENT(IN) :: thickness
    TYPE (t_sea_ice),       INTENT(IN)                    :: p_ice
    TYPE(t_hydro_ocean_state)                             :: p_os
    TYPE(t_sfc_flx)                                       :: surface_fluxes
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type

    ! locals
    REAL(wp), DIMENSION(nproma,n_zlev,p_patch%alloc_cell_blocks) :: salt

    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: saltInSeaice, saltInLiquidWater
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd, level
    INTEGER                                               :: my_computation_type

    my_computation_type = 0
    salt         = 0.0_wp

    CALL assign_if_present(my_computation_type, computation_type)
    subset => p_patch%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
        ! surface:
          saltInSeaice(cell,block)    = sice*rhoi &
            &                         * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                         * p_patch%cells%area(cell,block)
        !!DN This is no longer needed since we now update surface salinity
        !directly
        !!DN   p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIce(cell,block) &
        !!DN   &                                            -   dtime &
        !!DN   &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
        !!DN   &                                              * p_os%p_prog(nold(1))%tracer(cell,1,block,2)) &
        !!DN   &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        CASE (1) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
        ! surface:
          saltInSeaice(cell,block)    = sice*rhoi &
            &                         * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                         * p_patch%cells%area(cell,block)
        !!DN This is no longer needed since we now update surface salinity directly
          p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIce(cell,block) &
          &                                            -   dtime &
          &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
          &                                              * p_os%p_prog(nold(1))%tracer(cell,1,block,2)) &
          &                                           /zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        END SELECT

        salt(cell,1,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)
        DO level=2,subset%vertical_levels(cell,block)
          salt(cell,level,block) = p_os%p_prog(nold(1))%tracer(cell,level,block,2) &
            &                    * thickness(cell,level,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        END DO

        ! rest of the underwater world
      END DO ! cell
    END DO !block
  END FUNCTION calc_salt_content

END MODULE mo_ocean_diagnostics
